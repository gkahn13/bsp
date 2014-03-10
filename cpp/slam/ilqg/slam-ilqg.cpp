#include <vector>
#include <iomanip>

#include "util/Timer.h"

#include "../slam.h"
#include "ilqg.h"

const double alpha_belief = 10, alpha_final_belief = 10, alpha_control = 1, alpha_goal_state = 1;

// iLQG params
SymmetricMatrix<U_DIM> Rint;
SymmetricMatrix<X_DIM> Qint;
SymmetricMatrix<X_DIM> QGoal;
SymmetricMatrix<X_DIM> Sigma0;

Matrix<X_DIM> f(const Matrix<X_DIM>& x, const Matrix<U_DIM>& u)
{
	return dynfunc(x, u, zeros<Q_DIM,1>());
}

// Observation model
Matrix<Z_DIM> h(const Matrix<X_DIM>& x)
{
	return obsfunc(x, zeros<R_DIM,1>());
}

/*
double computeLQGMPcost(const std::vector<Matrix<X_DIM> >& X, const std::vector<Matrix<U_DIM> >& U)
{
	std::vector<Matrix<G_DIM, X_DIM> > J;
	std::vector<Matrix<X_DIM, X_DIM> > Sigma;
	preprocess(X, J, Sigma); // TODO: need to do same preprocess?

	double cost = 0;
	for(int t = 0; t < T-1; ++t) {
		cost += alpha_belief*tr(J[t]*Sigma[t]*~J[t]) + alpha_control*tr(~U[t]*U[t]);
	}
	cost += alpha_final_belief*tr(J[T-1]*Sigma[T-1]*~J[T-1]);
	return cost;
}
*/

double computeCost(const std::vector< Matrix<B_DIM> >& B, const std::vector< Matrix<U_DIM> >& U)
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

// Jacobian df/dx(x,u)
template <size_t _xDim, size_t _uDim>
inline Matrix<_xDim,_xDim> dfdx(Matrix<_xDim> (*f)(const Matrix<_xDim>&, const Matrix<_uDim>&), const Matrix<_xDim>& x, const Matrix<_uDim>& u) 
{
	Matrix<_xDim,_xDim> A;
	Matrix<_xDim> xr(x), xl(x);
	for (size_t i = 0; i < _xDim; ++i) {
		xr[i] += step; xl[i] -= step;
		A.insert(0,i, (f(xr, u) - f(xl, u)) / (2.0*step));
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
		ur[i] += step; ul[i] -= step;
		B.insert(0,i, (f(x, ur) - f(x, ul)) / (2.0*step));
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
		xr[i] += step; xl[i] -= step;
		H.insert(0,i, (h(xr) - h(xl)) / (2.0*step));
		xr[i] = xl[i] = x[i];
	}
	return H;
}

inline SymmetricMatrix<X_DIM> varM(const Matrix<X_DIM>& x, const Matrix<U_DIM>& u)
{
	Matrix<C_DIM,C_DIM> Acar;
	Matrix<C_DIM,Q_DIM> Mcar;
	linearizeDynamics(x, u, zeros<Q_DIM,1>(), Acar, Mcar);

	Matrix<X_DIM,Q_DIM> M = zeros<X_DIM,Q_DIM>();
	M.insert<C_DIM, 2>(0, 0, Mcar);

	SymmetricMatrix<X_DIM> S = SymProd(M*Q,~M);

	return S;
}

inline SymmetricMatrix<Z_DIM> varN(const Matrix<X_DIM>& x)
{
	return R;
}

// TODO: needs to be changed?
// Compute closed form quadratic finalCost function around around b
inline void quadratizeFinalCost(const Matrix<X_DIM>& xBar, const SymmetricMatrix<X_DIM>& SigmaBar, double& s, SymmetricMatrix<X_DIM>& S, Matrix<1,X_DIM>& sT, Matrix<1,S_DIM>& tT, unsigned int flag) 
{
	if (flag & COMPUTE_S) S = QGoal;
	if (flag & COMPUTE_sT) sT = ~(xBar - xGoal)*QGoal;
	if (flag & COMPUTE_s) s = 0.5*scalar(~(xBar - xGoal)*QGoal*(xBar - xGoal)) + scalar(vecTh(QGoal)*vec(SigmaBar));
	//if (flag & COMPUTE_s) s = 2.0*(0.5*scalar(~(xBar - xGoal)*QGoal*(xBar - xGoal)) + scalar(vecTh(QGoal)*vec(SigmaBar)));
	if (flag & COMPUTE_tT) tT = vecTh(QGoal);
}

// TODO: needs to be changed?
// Compute closed form quadratic cost function around around b and u
inline bool quadratizeCost(const Matrix<X_DIM>& xBar, const SymmetricMatrix<X_DIM>& SigmaBar, const Matrix<U_DIM>& uBar, double& q, SymmetricMatrix<X_DIM>& Q, SymmetricMatrix<U_DIM>& R, 
						   Matrix<U_DIM, X_DIM>& P, Matrix<1,X_DIM>& qT, Matrix<1,U_DIM>& rT, Matrix<1,S_DIM>& pT, unsigned int flag) 
{
	if (flag & COMPUTE_Q) Q = zeros<X_DIM>(); // penalize uncertainty
	if (flag & COMPUTE_R) R = Rint;
	if (flag & COMPUTE_P) P = zeros<U_DIM,X_DIM>();
	if (flag & COMPUTE_qT) qT = zeros<1, X_DIM>();
	if (flag & COMPUTE_rT) rT = ~uBar*Rint;
	if (flag & COMPUTE_pT) pT = vecTh(Qint);
	if (flag & COMPUTE_q) q = 0.5*scalar(~uBar*Rint*uBar) + scalar(vecTh(Qint)*vec(SigmaBar));
	//if (flag & COMPUTE_q) q = 2*(0.5*scalar(~uBar*Rint*uBar) + scalar(vecTh(Qint)*vec(SigmaBar)));

	return true;
}

inline void linearizeDynamics(const Matrix<X_DIM>& xBar, const Matrix<U_DIM>& uBar, Matrix<X_DIM>& c, Matrix<X_DIM, X_DIM>& A, Matrix<X_DIM, U_DIM>& B, SymmetricMatrix<X_DIM>& M, unsigned int flag)
{
	if (flag & COMPUTE_c) c = f(xBar, uBar);
	if (flag & COMPUTE_A) A = dfdx(f, xBar, uBar);
	if (flag & COMPUTE_B) B = dfdu(f, xBar, uBar);
	if (flag & COMPUTE_M) M = varM(xBar, uBar);
}

inline void linearizeObservation(const Matrix<X_DIM>& xBar, Matrix<Z_DIM, X_DIM>& H, SymmetricMatrix<Z_DIM>& N)
{
	H = dhdx(h, xBar);
	N = varN(xBar);
}

/*
bool testInitializationFeasibility(const std::vector<Matrix<B_DIM> >& B, const std::vector<Matrix<U_DIM> >& U)
{
	//std::cout << "X initial" << std::endl;
	for (int t = 0; t < T; ++t) { 
		const Matrix<U_DIM>& xt = B[t].subMatrix<X_DIM,1>(0,0);
		for (int i = 0; i < X_DIM; ++i) {
			if (xt[i] > xMax[i] || xt[i] < xMin[i]) {
				LOG_ERROR("Joint angle limit violated at joint %d and time %d", i,t);
				return false;
			}
		}

		//std::cout << std::setprecision(8) << ~xt; 
	}

	//std::cout << "U initial" << std::endl;
	for (int t = 0; t < T-1; ++t) { 
		const Matrix<U_DIM>& ut = U[t];
		for (int i = 0; i < U_DIM; ++i) {
			if (ut[i] > uMax[i] || ut[i] < uMin[i]) {
				LOG_ERROR("Control limit violated at joint %d and time %d", i,t);
				return false;
			}
		}
		//std::cout << std::setprecision(8) << ~ut; 
	}

	return true;
}
*/

int main(int argc, char* argv[])
{
	std::vector<Matrix<P_DIM> > l(NUM_LANDMARKS);
	l[0][0] = 30; l[0][1] = -10;
	l[1][0] = 70; l[1][1] = 12.5;
	l[2][0] = 20; l[2][1] = 10;
	initProblemParams(l);

	Q *= 1;
	R *= 1;

	Rint = alpha_control*identity<U_DIM>();
	Qint = alpha_belief*identity<X_DIM>();
	QGoal = alpha_final_belief*identity<X_DIM>();


	std::vector< Matrix<U_DIM, X_DIM> > L;
	std::vector< Matrix<X_DIM> > xBar;
	std::vector< SymmetricMatrix<X_DIM> > SigmaBar;
	
	xGoal = x0;
	xGoal[0] = 60;
	xGoal[1] = 0;
	xGoal[2] = 0;

	// TODO: should I initialize using slam-traj?
	Matrix<U_DIM> uinit;
	uinit[0] = (xGoal[0] - x0[0])/((double)(T-1)*DT);
	uinit[1] = 0;
	std::vector<Matrix<U_DIM> > uBar(T-1, uinit);

	uBar[0][1] = 3.14/6;
	uBar[6][1] = -2*3.14/6;

	Sigma0 = SymProd(SqrtSigma0,SqrtSigma0); //identity<X_DIM>();
	xBar.push_back(x0);

	SigmaBar.push_back(Sigma0);

//	std::cout << "Rint\n" << Rint;
//	std::cout << "Qint\n" << Qint;
//	std::cout << "QGoal\n" << QGoal;
//	std::cout << "uinit: " << ~uinit;
//	std::cout << "x0: " << ~x0;
//	std::cout << "xGoal: " << ~xGoal;
//	std::cout << "varM\n" << varM(x0, uinit);
//	std::cout << "varN\n" << varN(x0);

	for(int t=0; t < T-1; ++t) {
		std::cout << ~uBar[t];
	}

	//pythonDisplayTrajectory(uBar, T, true);

	util::Timer solveTimer;
	Timer_tic(&solveTimer);

	// solve ilqg here
	solvePOMDP(linearizeDynamics, linearizeObservation, quadratizeFinalCost, quadratizeCost, xBar, SigmaBar, uBar, L);

	double solvetime = util::Timer_toc(&solveTimer);

	LOG_INFO("Solve time: %5.3f ms", solvetime*1000);
	
	pythonDisplayTrajectory(uBar, T, true);

	for(int t=0; t < T-1; ++t) {
		std::cout << ~uBar[t];
	}


	return 0;
}
