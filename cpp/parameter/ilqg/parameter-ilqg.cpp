#include "../parameter.h"

#include "util/matrix.h"

#include "ilqg.h"

#include <time.h>

#include <Python.h>
#include <boost/python.hpp>
#include <boost/filesystem.hpp>

namespace py = boost::python;


SymmetricMatrix<U_DIM> Rint;
SymmetricMatrix<X_DIM> Qint;
SymmetricMatrix<X_DIM> QGoal, QGoalVariance, QintVariance;
SymmetricMatrix<X_DIM> Sigma0;
Matrix<X_DIM,X_DIM> SqrtSigma0;
Matrix<X_DIM,X_DIM> SqrtTemp;

Matrix<X_DIM> x0, xGoal;

inline Matrix<X_DIM> f(const Matrix<X_DIM>& x, const Matrix<U_DIM>& u)
{
	return dynfunc(x, u,zeros<Q_DIM,1>());
}

// Observation model
inline Matrix<Z_DIM> h(const Matrix<X_DIM>& x)
{
	return obsfunc(x,zeros<R_DIM,1>());
}

// Jacobian df/dx(x,u)
template <size_t _xDim, size_t _uDim>
inline Matrix<_xDim,_xDim> dfdx(Matrix<_xDim> (*f)(const Matrix<_xDim>&, const Matrix<_uDim>&), const Matrix<_xDim>& x, const Matrix<_uDim>& u) 
{
	Matrix<_xDim,_xDim> A;
	Matrix<_xDim> xr(x), xl(x);
	for (size_t i = 0; i < _xDim; ++i) {
		xr[i] += diffEps; xl[i] -= diffEps;
		A.insert(0,i, (f(xr, u) - f(xl, u)) / (2.0*diffEps));
		xr[i] = xl[i] = x[i];
	}
	return A;
}

// Jacobian df/du(x,u)
template <size_t _xDim, size_t _uDim>
inline Matrix<_xDim,_uDim> dfdu(Matrix<_xDim> (*f)(const Matrix<_xDim>&, const Matrix<_uDim>&),	const Matrix<_xDim>& x, const Matrix<_uDim>& u) 
{
	Matrix<_xDim,_uDim> B;
	Matrix<_uDim> ur(u), ul(u);
	for (size_t i = 0; i < _uDim; ++i) {
		ur[i] += diffEps; ul[i] -= diffEps;
		B.insert(0,i, (f(x, ur) - f(x, ul)) / (2.0*diffEps));
		ur[i] = ul[i] = u[i];
	}
	return B;
}

// Jacobian dh/dx(x)
template <size_t _xDim, size_t _zDim>
inline Matrix<_zDim,_xDim> dhdx(Matrix<_zDim> (*h)(const Matrix<_xDim>&), const Matrix<_xDim>& x) 
{
	Matrix<_zDim,_xDim> H;
	Matrix<_xDim> xr(x), xl(x);
	for (size_t i = 0; i < _xDim; ++i) {
		xr[i] += diffEps; xl[i] -= diffEps;
		H.insert(0,i, (h(xr) - h(xl)) / (2.0*diffEps));
		xr[i] = xl[i] = x[i];
	}
	return H;
}

inline SymmetricMatrix<X_DIM> varM(const Matrix<X_DIM>& x, const Matrix<U_DIM>& u)
{
	SymmetricMatrix<X_DIM> S = identity<X_DIM>();
	S(0,0) = 0.01;
	S(1,1) = 0.01;
	S(2,2) = 0.01;
	S(3,3) = 0.01;
	S(4,4) = 1e-6;
	S(5,5) = 1e-6;
	S(6,6) = 1e-6;
	S(7,7) = 1e-6;
	return S;
}

inline SymmetricMatrix<Z_DIM> varN(const Matrix<X_DIM>& x)
{
	SymmetricMatrix<Z_DIM> S = identity<Z_DIM>();
	S(0,0) = 1e-2;
	S(1,1) = 1e-2;
	S(2,2) = 1e-3;
	S(3,3) = 1e-3;
	return S;
}

// Compute closed form quadratic finalCost function around around b
inline void quadratizeFinalCost(const Matrix<X_DIM>& xBar, const SymmetricMatrix<X_DIM>& SigmaBar, double& s, SymmetricMatrix<X_DIM>& S, Matrix<1,X_DIM>& sT, Matrix<1,S_DIM>& tT, unsigned int flag) 
{
	if (flag & COMPUTE_S) S = QGoal;
	if (flag & COMPUTE_sT) sT = ~(xBar - xGoal)*QGoal;
	if (flag & COMPUTE_s) s = 0.5*scalar(~(xBar - xGoal)*QGoal*(xBar - xGoal)) + scalar(vecTh(QGoalVariance)*vectorize(SigmaBar)); // TODO: this cost is huge, why?
	if (flag & COMPUTE_tT) tT = vecTh(QGoalVariance);
}

// Compute closed form quadratic cost function around around b and u
inline bool quadratizeCost(const Matrix<X_DIM>& xBar, const SymmetricMatrix<X_DIM>& SigmaBar, const Matrix<U_DIM>& uBar, double& q, SymmetricMatrix<X_DIM>& Q, SymmetricMatrix<U_DIM>& R, 
						   Matrix<U_DIM, X_DIM>& P, Matrix<1,X_DIM>& qT, Matrix<1,U_DIM>& rT, Matrix<1,S_DIM>& pT, unsigned int flag) 
{
	if (flag & COMPUTE_Q) Q = Qint; // penalize uncertainty
	if (flag & COMPUTE_R) R = Rint;
	if (flag & COMPUTE_P) P = zeros<U_DIM,X_DIM>();
	if (flag & COMPUTE_qT) qT = ~(xBar - xGoal)*Qint;
	if (flag & COMPUTE_rT) rT = ~uBar*Rint;
	if (flag & COMPUTE_pT) pT = vecTh(QintVariance);
	if (flag & COMPUTE_q) q = 0.5*scalar(~(xBar - xGoal)*Qint*(xBar - xGoal)) + 0.5*scalar(~uBar*Rint*uBar) + scalar(vecTh(QintVariance)*vectorize(SigmaBar));
	return true;
}

inline void linearizeDynamics(const Matrix<X_DIM>& xBar, const Matrix<U_DIM>& uBar, Matrix<X_DIM>& c, Matrix<X_DIM, X_DIM>& A, Matrix<X_DIM, U_DIM>& B, SymmetricMatrix<X_DIM>& M, unsigned int flag)
{
	if (flag & COMPUTE_c) c = f(xBar, uBar);
	if (flag & COMPUTE_A) A = dfdx(f, xBar, uBar); //identity<X_DIM>();
	if (flag & COMPUTE_B) B = dfdu(f, xBar, uBar); //identity<U_DIM>();
	if (flag & COMPUTE_M) M = varM(xBar, uBar);
}

inline void linearizeObservation(const Matrix<X_DIM>& xBar, Matrix<Z_DIM, X_DIM>& H, SymmetricMatrix<Z_DIM>& N)
{
	H = dhdx(h, xBar); // identity<X_DIM>();
	N = varN(xBar);
}


double costfunc(const std::vector< Matrix<B_DIM> >& B, const std::vector< Matrix<U_DIM> >& U)
{
	double cost = 0;
	Matrix<X_DIM> x;
	Matrix<X_DIM, X_DIM> SqrtSigma;

	for(int t = 0; t < T-1; ++t) {
		unVec(B[t], x, SqrtSigma);
		cost += alpha_belief*tr(SqrtSigma*SqrtSigma) + alpha_control*tr(~U[t]*U[t]);
	}
	unVec(B[T-1], x, SqrtSigma);
	cost += alpha_final_belief*tr(SqrtSigma*SqrtSigma);
	return cost;
}

Matrix<X_DIM> SigmaDiag(const Matrix<X_DIM,X_DIM>& Sigma) {
	Matrix<X_DIM> diag;
	for(int i = 0; i < X_DIM; ++i) {
		diag[i] = Sigma(i,i);
	}
	return diag;
}

int main(int argc, char* argv[])
{
	double length1_est = .3,
			length2_est = .7,
			mass1_est = .3,
			mass2_est = .35;

	// position, then velocity
	x0[0] = M_PI*0.5; x0[1] = M_PI*0.5; x0[2] = 0; x0[3] = 0;
	// parameter start estimates (alphabetical, then numerical order)
	x0[4] = 1/length1_est; x0[5] = 1/length2_est; x0[6] = 1/mass1_est; x0[7] = 1/mass2_est;


	Matrix<X_DIM> x_real;
	x_real[0] = M_PI*0.45; x_real[1] = M_PI*0.55; x_real[2] = -0.01; x_real[3] = 0.01;
	x_real[4] = 1/dynamics::length1; x_real[5] = 1/dynamics::length2; x_real[6] = 1/dynamics::mass1; x_real[7] = 1/dynamics::mass2;



	// init controls from start to goal -- straight line trajectory
	// TODO: possibly switch to random controls
	std::vector< Matrix<U_DIM> > uBar((T-1), (x0.subMatrix<U_DIM,1>(0,0) - x0.subMatrix<U_DIM,1>(0,0))/(double)(T-1));

	Rint = 1e-3*identity<U_DIM>();

	Qint(0,0) = 10.0;
	Qint(1,1) = 10.0;
	Qint(2,2) = 10.0;
	Qint(3,3) = 10.0;
	Qint(4,4) = 0.0;
	Qint(5,5) = 0.0;
	Qint(6,6) = 0.0;
	Qint(7,7) = 0.0;

	QintVariance(0,0) = 0.0;
	QintVariance(1,1) = 0.0;
	QintVariance(2,2) = 0.0;
	QintVariance(3,3) = 0.0;
	QintVariance(4,4) = 10.0;
	QintVariance(5,5) = 10.0;
	QintVariance(6,6) = 10.0;
	QintVariance(7,7) = 10.0;

	QGoal(0,0) = 0.0;
	QGoal(1,1) = 0.0;
	QGoal(2,2) = 0.0;
	QGoal(3,3) = 0.0;
	QGoal(4,4) = 0.0;
	QGoal(5,5) = 0.0;
	QGoal(6,6) = 0.0;
	QGoal(7,7) = 0.0;

	QGoalVariance(0,0) = 0.0;
	QGoalVariance(1,1) = 0.0;
	QGoalVariance(2,2) = 0.0;
	QGoalVariance(3,3) = 0.0;
	QGoalVariance(4,4) = 0.0;
	QGoalVariance(5,5) = 0.0;
	QGoalVariance(6,6) = 0.0;
	QGoalVariance(7,7) = 0.0;

	std::vector< Matrix<U_DIM, X_DIM> > L;
	std::vector< Matrix<X_DIM> > xBar;
	std::vector< SymmetricMatrix<X_DIM> > SigmaBar;
	
	SqrtSigma0(0,0) = 0.1;
	SqrtSigma0(1,1) = 0.1;
	SqrtSigma0(2,2) = 0.05;
	SqrtSigma0(3,3) = 0.05;
	SqrtSigma0(4,4) = 0.5;
	SqrtSigma0(5,5) = 0.5;
	SqrtSigma0(6,6) = 0.5;
	SqrtSigma0(7,7) = 0.5;
	for(int i=0; i < X_DIM; ++i) { Sigma0(i,i) = SqrtSigma0(i,i)*SqrtSigma0(i,i); } 
	xBar.push_back(x0);
	SigmaBar.push_back(Sigma0);

	std::vector<Matrix<U_DIM> > HistoryU(HORIZON);
	std::vector<Matrix<B_DIM> > HistoryB(HORIZON);

	std::vector< Matrix<B_DIM> > Bpomdp(T);
	std::vector< Matrix<B_DIM> > Bekf(T);
	vec(xBar[0], sqrt(SigmaBar[0]), Bekf[0]);
	std::vector< Matrix<B_DIM> > Binitial(T);
	vec(x0, sqrt(Sigma0), Binitial[0]);
	for (int t = 0; t < T - 1; ++t)
	{
		Binitial[t+1] = beliefDynamics(Binitial[t], uBar[t]);
	}
	for (int h=0; h<HORIZON; h++){
		
		double cost_initial = costfunc(Binitial, uBar);

		std::cout << "solvePOMDP " <<h<< std::endl;

		solvePOMDP(linearizeDynamics, linearizeObservation, quadratizeFinalCost, quadratizeCost, xBar, SigmaBar, uBar, L);

		std::cout << "POMDP solved" << std::endl;

	//for (size_t i = 0; i < xBar.size(); ++i) {
	//	std::cout << ~(xBar[i]);
	//}
		for(int t = 0; t < T; ++t) {
			vec(xBar[t], sqrt(SigmaBar[t]), Bpomdp[t]);
		}

		
		

		HistoryU[h] = uBar[0];
		HistoryB[h] = Bekf[0];
		Bekf[0] = executeControlStep(x_real, Bekf[0], uBar[0]);
		std::cout<<"U "<<~uBar[0]<<"\n";
		std::cout<<"X "<<~xBar[0]<<"\n";
		for (int t = 0; t < T - 1; ++t)
		{
			Bekf[t+1] = beliefDynamics(Bekf[t], uBar[t]);
		}



		double cost_pomdp = costfunc(Bpomdp, uBar);
		double cost_ekf = costfunc(Bekf, uBar);

		std::cout << "cost initial: " << cost_initial << std::endl;
		std::cout << "cost pomdp: " << cost_pomdp << std::endl;
		std::cout << "cost ekf: " << cost_ekf << std::endl;

		//Update XBar, SigmaBar, uBar
		xBar.clear(); 
		SigmaBar.clear(); 

		unVec(Bekf[0], x0, SqrtSigma0);

		xBar.push_back(x0); 

		for(int i=0; i < X_DIM; ++i) { Sigma0(i,i) = SqrtSigma0(i,i)*SqrtSigma0(i,i); } 

		SigmaBar.push_back(Sigma0); 

		for(int t = 0; t < T-2; ++t) {
		
			uBar[t] = zeros<U_DIM,1>();
		}


	}
	pythonDisplayHistory(HistoryU,HistoryB, SqrtSigma0, x0, HORIZON);
	Matrix<X_DIM> xpomdp, xekf;
	Matrix<X_DIM,X_DIM> SqrtSigmapomdp, SqrtSigmaekf;
	for(int t = 0; t < T; ++t) {
		unVec(Bpomdp[t], xpomdp, SqrtSigmapomdp);
		unVec(Bekf[t], xekf, SqrtSigmaekf);
		//std::cout << "t: " << t << " xpomdp" << std::endl;
		//std::cout << ~xpomdp;
		//std::cout << "t: " << t << " xekf" << std::endl;
		//std::cout << ~xekf << std::endl;
		std::cout << "t: " << t << " Sigmapomdp" << std::endl;
		std::cout << ~SigmaDiag(SqrtSigmapomdp);
		std::cout << "t: " << t << " Sigmaekf" << std::endl;
		std::cout << ~SigmaDiag(SqrtSigmaekf) << std::endl << std::endl;
	}


#ifdef CPP_PLOT

	Py_Initialize();

	py::list Bvec;
	for(int j=0; j < B_DIM; j++) {
		for(int i=0; i < T; i++) {
			Bvec.append(Bekf[i][j]);
		}
	}

	py::list Uvec;
	for(int j=0; j < U_DIM; j++) {
		for(int i=0; i < T-1; i++) {
			Uvec.append(uBar[i][j]);
		}
	}

	std::string workingDir = boost::filesystem::current_path().normalize().string();

	py::object main_module = py::import("__main__");
	py::object main_namespace = main_module.attr("__dict__");
	py::exec("import sys, os", main_namespace);
	py::exec(py::str("sys.path.append('"+workingDir+"/parameter')"), main_namespace);
	py::object plot_mod = py::import("plot_parameter");
	py::object plot_traj = plot_mod.attr("plot_parameter_trajectory");

	plot_traj(Bvec, Uvec, B_DIM, X_DIM, U_DIM, T);

	//int k;
	//std::cin >> k;
#endif


	return 0;
}

