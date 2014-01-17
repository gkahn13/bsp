#include <vector>

#include "util/matrix.h"
#include "util/utils.h"
#include "util/Timer.h"
#include "util/logging.h"

#include <Python.h>
#include <boost/python.hpp>
#include <numpy/ndarrayobject.h>
#include <boost/filesystem.hpp>

namespace py = boost::python;

extern "C" {
#include "beliefMPC.h"
}

#define DT 1.0
#define X_DIM 2
#define U_DIM 2
#define Z_DIM 2
#define Q_DIM 2
#define R_DIM 2

#define S_DIM (((X_DIM+1)*X_DIM)/2)
#define B_DIM (X_DIM+S_DIM)

const double step = 0.0078125*0.0078125;

Matrix<X_DIM> x0;
Matrix<X_DIM,X_DIM> SqrtSigma0;
Matrix<X_DIM> xGoal;
Matrix<X_DIM> xMin, xMax;
Matrix<U_DIM> uMin, uMax;

const int T = 15;
const double INFTY = 1e10;
const double alpha_belief = 10, alpha_final_belief = 10, alpha_control = 1;

// beliefMPC vars
beliefMPC_FLOAT **lb, **ub, **C, **e, **z;

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

inline Matrix<X_DIM> f(const Matrix<X_DIM>& x, const Matrix<U_DIM>& u, const Matrix<U_DIM>& q)
{  
	Matrix<X_DIM> xNew = x + u*DT + 0.01*q;
	return xNew;
}

// Observation model
inline Matrix<Z_DIM> h(const Matrix<X_DIM>& x, const Matrix<R_DIM>& r)
{
	double intensity = sqrt(sqr(0.5*x[0]) + 1e-6);
	Matrix<Z_DIM> z = x + intensity*r;
	return z;
}

// Jacobians: df(x,u,q)/dx, df(x,u,q)/dq
inline void linearizeDynamics(const Matrix<X_DIM>& x, const Matrix<U_DIM>& u, const Matrix<Q_DIM>& q, Matrix<X_DIM,X_DIM>& A, Matrix<X_DIM,Q_DIM>& M)
{
	A.reset();
	Matrix<X_DIM> xr(x), xl(x);
	for (size_t i = 0; i < X_DIM; ++i) {
		xr[i] += step; xl[i] -= step;
		A.insert(0,i, (f(xr, u, q) - f(xl, u, q)) / (xr[i] - xl[i]));
		xr[i] = x[i]; xl[i] = x[i];
	}

	M.reset();
	Matrix<Q_DIM> qr(q), ql(q);
	for (size_t i = 0; i < Q_DIM; ++i) {
		qr[i] += step; ql[i] -= step;
		M.insert(0,i, (f(x, u, qr) - f(x, u, ql)) / (qr[i] - ql[i]));
		qr[i] = q[i]; ql[i] = q[i];
	}
}

// Jacobians: dh(x,r)/dx, dh(x,r)/dr
inline void linearizeObservation(const Matrix<X_DIM>& x, const Matrix<R_DIM>& r, Matrix<Z_DIM,X_DIM>& H, Matrix<Z_DIM,R_DIM>& N)
{
	H.reset();
	Matrix<X_DIM> xr(x), xl(x);
	for (size_t i = 0; i < X_DIM; ++i) {
		xr[i] += step; xl[i] -= step;
		H.insert(0,i, (h(xr, r) - h(xl, r)) / (xr[i] - xl[i]));
		xr[i] = x[i]; xl[i] = x[i];
	}

	N.reset();
	Matrix<R_DIM> rr(r), rl(r);
	for (size_t i = 0; i < R_DIM; ++i) {
		rr[i] += step; rl[i] -= step;
		N.insert(0,i, (h(x, rr) - h(x, rl)) / (rr[i] - rl[i]));
		rr[i] = r[i]; rl[i] = r[i];
	}
}

// Switch between belief vector and matrices
inline void unVec(const Matrix<B_DIM>& b, Matrix<X_DIM>& x, Matrix<X_DIM,X_DIM>& SqrtSigma) {
	x = b.subMatrix<X_DIM,1>(0,0);
	size_t idx = X_DIM;
	for (size_t j = 0; j < X_DIM; ++j) {
		for (size_t i = j; i < X_DIM; ++i) {
			SqrtSigma(i,j) = b[idx];
			SqrtSigma(j,i) = b[idx];
			++idx;
		}
	}
}

inline void vec(const Matrix<X_DIM>& x, const Matrix<X_DIM,X_DIM>& SqrtSigma, Matrix<B_DIM>& b) {
	b.insert(0,0,x);
	size_t idx = X_DIM;
	for (size_t j = 0; j < X_DIM; ++j) {
		for (size_t i = j; i < X_DIM; ++i) {
			b[idx] = 0.5 * (SqrtSigma(i,j) + SqrtSigma(j,i));
			++idx;
		}
	}
}

// Belief dynamics
inline Matrix<B_DIM> beliefDynamics(const Matrix<B_DIM>& b, const Matrix<U_DIM>& u) {
	Matrix<X_DIM> x;
	Matrix<X_DIM,X_DIM> SqrtSigma;
	unVec(b, x, SqrtSigma);

	Matrix<X_DIM,X_DIM> Sigma = SqrtSigma*SqrtSigma;

	Matrix<X_DIM,X_DIM> A;
	Matrix<X_DIM,Q_DIM> M;
	linearizeDynamics(x, u, zeros<Q_DIM,1>(), A, M);

	x = f(x, u, zeros<Q_DIM,1>());
	Sigma = A*Sigma*~A + M*~M;

	Matrix<Z_DIM,X_DIM> H;
	Matrix<Z_DIM,R_DIM> N;
	linearizeObservation(x, zeros<R_DIM,1>(), H, N);

	Matrix<X_DIM,Z_DIM> K = Sigma*~H/(H*Sigma*~H + N*~N);

	Sigma = (identity<X_DIM>() - K*H)*Sigma;

	Matrix<B_DIM> g;
	vec(x, sqrt(Sigma), g);

	return g;
}

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

// Jacobians: dg(b,u)/db, dg(b,u)/du
void linearizeBeliefDynamics(const Matrix<B_DIM>& b, const Matrix<U_DIM>& u, Matrix<B_DIM,B_DIM>& F, Matrix<B_DIM,U_DIM>& G, Matrix<B_DIM>& h)
{
	F.reset();
	Matrix<B_DIM> br(b), bl(b);
	for (size_t i = 0; i < B_DIM; ++i) {
		br[i] += step; bl[i] -= step;
		//std::cout << "bplus: " << ~(beliefDynamics(br, u)) << std::endl;
		//std::cout << "bminus: " << ~(beliefDynamics(bl, u)) << std::endl;
		F.insert(0,i, (beliefDynamics(br, u) - beliefDynamics(bl, u)) / (br[i] - bl[i]));
		br[i] = b[i]; bl[i] = b[i];
	}

	G.reset();
	Matrix<U_DIM> ur(u), ul(u);
	for (size_t i = 0; i < U_DIM; ++i) {
		ur[i] += step; ul[i] -= step;
		G.insert(0,i, (beliefDynamics(b, ur) - beliefDynamics(b, ul)) / (ur[i] - ul[i]));
		ur[i] = u[i]; ul[i] = u[i];
	}

	h = beliefDynamics(b, u);
}

// TODO: Find better way to do this using macro expansions?
void setupBeliefMPCVars(beliefMPC_params& problem, beliefMPC_output& output)
{
	lb = new beliefMPC_FLOAT*[T];
	ub = new beliefMPC_FLOAT*[T];
	C = new beliefMPC_FLOAT*[T-1];
	e = new beliefMPC_FLOAT*[T-1];
	z = new beliefMPC_FLOAT*[T];

	lb[0] = problem.lb01; ub[0] = problem.ub01; C[0] = problem.C01; e[0] = problem.e01;
	lb[1] = problem.lb02; ub[1] = problem.ub02; C[1] = problem.C02; e[1] = problem.e02;
	lb[2] = problem.lb03; ub[2] = problem.ub03; C[2] = problem.C03; e[2] = problem.e03;
	lb[3] = problem.lb04; ub[3] = problem.ub04; C[3] = problem.C04; e[3] = problem.e04;
	lb[4] = problem.lb05; ub[4] = problem.ub05; C[4] = problem.C05; e[4] = problem.e05;
	lb[5] = problem.lb06; ub[5] = problem.ub06; C[5] = problem.C06; e[5] = problem.e06;
	lb[6] = problem.lb07; ub[6] = problem.ub07; C[6] = problem.C07; e[6] = problem.e07;
	lb[7] = problem.lb08; ub[7] = problem.ub08; C[7] = problem.C08; e[7] = problem.e08;
	lb[8] = problem.lb09; ub[8] = problem.ub09; C[8] = problem.C09; e[8] = problem.e09;
	lb[9] = problem.lb10; ub[9] = problem.ub10; C[9] = problem.C10; e[9] = problem.e10;
	lb[10] = problem.lb11; ub[10] = problem.ub11; C[10] = problem.C11; e[10] = problem.e11;
	lb[11] = problem.lb12; ub[11] = problem.ub12; C[11] = problem.C12; e[11] = problem.e12;
	lb[12] = problem.lb13; ub[12] = problem.ub13; C[12] = problem.C13; e[12] = problem.e13;
	lb[13] = problem.lb14; ub[13] = problem.ub14; C[13] = problem.C14; e[13] = problem.e14;
	lb[14] = problem.lb15; ub[14] = problem.ub15;

	z[0] = output.z1; z[1] = output.z2; z[2] = output.z3; z[3] = output.z4; z[4] = output.z5;
	z[5] = output.z6; z[6] = output.z7; z[7] = output.z8; z[8] = output.z9; z[9] = output.z10; 
	z[10] = output.z11; z[11] = output.z12; z[12] = output.z13; z[13] = output.z14; z[14] = output.z15; 
}

void cleanupBeliefMPCVars()
{
	delete[] lb; 
	delete[] ub; 
	delete[] C;
	delete[] e;
	delete[] z;
}

double beliefCollocation(std::vector< Matrix<B_DIM> >& B, std::vector< Matrix<U_DIM> >& U, beliefMPC_params& problem, beliefMPC_output& output, beliefMPC_info& info)
{
	int maxIter = 10;
	double Beps = 1;
	double Ueps = 1;

	// box constraint around goal
	double delta = 0.01;

	Matrix<B_DIM,1> b0 = B[0];

	std::vector< Matrix<B_DIM,B_DIM> > F(T-1);
	std::vector< Matrix<B_DIM,U_DIM> > G(T-1);
	std::vector< Matrix<B_DIM> > h(T-1);

	double prevcost = computeCost(B, U);
	double optcost;

	//std::cout << "Initialization trajectory cost: " << std::setprecision(10) << prevcost << std::endl;

	for(int it = 0; it < maxIter; ++it) 
	{
		//std::cout << "Iter: " << it << std::endl;

		// linearize belief dynamics constraint here
		for (int t = 0; t < T-1; ++t) 
		{
			Matrix<B_DIM>& bt = B[t];
			Matrix<U_DIM>& ut = U[t];

			linearizeBeliefDynamics(bt, ut, F[t], G[t], h[t]);

			// Fill in lb, ub, C, e
			lb[t][0] = MAX(xMin[0], bt[0] - Beps);
			lb[t][1] = MAX(xMin[1], bt[1] - Beps);
			lb[t][2] = bt[2] - Beps;
			lb[t][3] = bt[3] - Beps;
			lb[t][4] = bt[4] - Beps;
			lb[t][5] = MAX(uMin[0], ut[0] - Ueps);
			lb[t][6] = MAX(uMin[1], ut[1] - Ueps);

			ub[t][0] = MIN(xMax[0], bt[0] + Beps);
			ub[t][1] = MIN(xMax[1], bt[1] + Beps);
			ub[t][2] = bt[2] + Beps;
			ub[t][3] = bt[3] + Beps;
			ub[t][4] = bt[4] + Beps;
			ub[t][5] = MIN(uMax[0], ut[0] + Ueps);
			ub[t][6] = MIN(uMax[1], ut[1] + Ueps);

			if (t > 0) {
				Matrix<B_DIM,B_DIM+U_DIM> CMat;
				Matrix<B_DIM> eVec;

				CMat.insert<B_DIM,B_DIM>(0,0,F[t]);
				CMat.insert<B_DIM,U_DIM>(0,B_DIM,G[t]);
				Matrix<B_DIM+U_DIM, B_DIM> CMatT = ~CMat;
				int idx = 0;
				for(int c = 0; c < (B_DIM+U_DIM); ++c) {
					for(int r = 0; r < B_DIM; ++r) {
						C[t][idx++] = CMat[c + r*(B_DIM+U_DIM)];
					}
				}
				eVec = -h[t] + F[t]*bt + G[t]*ut;
				for(int i = 0; i < B_DIM; ++i) {
					e[t][i] = eVec[i];
				}
			}
			else {
				Matrix<2*B_DIM,B_DIM+U_DIM> CMat;
				Matrix<2*B_DIM> eVec;

				CMat.insert<B_DIM,B_DIM>(0,0,identity<B_DIM>());
				CMat.insert<B_DIM,U_DIM>(0,B_DIM,zeros<B_DIM,U_DIM>());
				CMat.insert<B_DIM,B_DIM>(B_DIM,0,F[t]);
				CMat.insert<B_DIM,U_DIM>(B_DIM,B_DIM,G[t]);
				int idx = 0;
				for(int c = 0; c < (B_DIM+U_DIM); ++c) {
					for(int r = 0; r < 2*B_DIM; ++r) {
						C[t][idx++] = CMat[c + r*(B_DIM+U_DIM)];
					}
				}
				eVec.insert<B_DIM,1>(0,0,bt);
				eVec.insert<B_DIM,1>(B_DIM,0,zeros<B_DIM,1>());
				for(int i = 0; i < 2*B_DIM; ++i) {
					e[t][i] = eVec[i];
				}
			}
		} //setting up problem

		Matrix<B_DIM>& bT = B[T-1];

		// Fill in lb, ub, C, e
		lb[T-1][0] = MAX(xGoal[0] - delta, bT[0] - Beps);
		lb[T-1][1] = MAX(xGoal[1] - delta, bT[1] - Beps);
		lb[T-1][2] = bT[2] - Beps;
		lb[T-1][3] = bT[3] - Beps;
		lb[T-1][4] = bT[4] - Beps;

		ub[T-1][0] = MIN(xGoal[0] + delta, bT[0] + Beps);
		ub[T-1][1] = MIN(xGoal[1] + delta, bT[1] + Beps);
		ub[T-1][2] = bT[2] + Beps;
		ub[T-1][3] = bT[3] + Beps;
		ub[T-1][4] = bT[4] + Beps;

		// Verify problem inputs

		//int num;
		//std::cin >> num;
		
		int exitflag = beliefMPC_solve(&problem, &output, &info);
		if (exitflag == 1) {
			for(int t = 0; t < T-1; ++t) {
				Matrix<B_DIM>& bt = B[t];
				Matrix<U_DIM>& ut = U[t];

				for(int i = 0; i < B_DIM; ++i) {
					bt[i] = z[t][i];
				}
				for(int i = 0; i < U_DIM; ++i) {
					ut[i] = z[t][B_DIM+i];
				}
				optcost = info.pobj;
			}
		}
		else {
			LOG_ERROR("Some problem in solver");
			std::exit(-1);
		}
		LOG_DEBUG("Optimized cost: %4.10f", optcost);

		if ((optcost > prevcost) || (fabs(optcost - prevcost)/prevcost < 0.01))
			break; 
		else {
			prevcost = optcost;
			// TODO: integrate trajectory?
			// TODO: plot trajectory
		}
		
		//int num;
		//std::cin >> num;
	}
	return computeCost(B, U);

}

// default for unix
// requires path to Python bsp already be on PYTHONPATH
void pythonDisplayTrajectory(std::vector< Matrix<B_DIM> >& B, std::vector< Matrix<U_DIM> >& U)
{
	for (int t = 0; t < T-1; ++t) {
		B[t+1] = beliefDynamics(B[t], U[t]);
	}

	py::list Bvec;
	for(int j=0; j < B_DIM; j++) {
		for(int i=0; i < T; i++) {
			Bvec.append(B[i][j]);
		}
	}

	py::list Uvec;
		for(int j=0; j < U_DIM; j++) {
			for(int i=0; i < T-1; i++) {
			Uvec.append(U[i][j]);
		}
	}

	std::string workingDir = boost::filesystem::current_path().normalize().string();
	std::string bspDir = workingDir.substr(0,workingDir.find("bsp"));

	try
	{
		Py_Initialize();
		py::object main_module = py::import("__main__");
		py::object main_namespace = main_module.attr("__dict__");
		py::exec("import sys, os", main_namespace);
		py::exec(py::str("sys.path.append('"+bspDir+"bsp/python')"), main_namespace);
		py::exec("from bsp_light_dark import LightDarkModel", main_namespace);
		py::object model = py::eval("LightDarkModel()", main_namespace);
		py::object plot_mod = py::import("plot");
		py::object plot_traj = plot_mod.attr("plot_belief_trajectory");

		plot_traj(Bvec, Uvec, model);
	}
	catch(py::error_already_set const &)
	{
		PyErr_Print();
	}

}

int main(int argc, char* argv[])
{
	x0[0] = -3.5; x0[1] = 2;
	SqrtSigma0 = identity<X_DIM>();
	xGoal[0] = -3.5; xGoal[1] = -2;

	xMin[0] = -5; xMin[1] = -3; 
	xMax[0] = 5; xMax[1] = 3;
	uMin[0] = -1; uMin[1] = -1;
	uMax[0] = 1; uMax[1] = 1;

	Matrix<U_DIM> uinit;
	uinit[0] = (xGoal[0] - x0[0]) / (T-1);
	uinit[1] = (xGoal[1] - x0[1]) / (T-1);
	
	std::vector<Matrix<U_DIM> > U(T-1, uinit); 

	std::vector<Matrix<B_DIM> > B(T);

	vec(x0, SqrtSigma0, B[0]);
	for (size_t t = 0; t < T-1; ++t) {
		B[t+1] = beliefDynamics(B[t], U[t]);
		//std::cout << ~B[t] << std::endl;
	}

	//for (size_t t = 0; t < T; ++t) {
	//	std::cout << ~B[t];
	//}

	beliefMPC_params problem;
	beliefMPC_output output;
	beliefMPC_info info;

	setupBeliefMPCVars(problem, output);

	util::Timer solveTimer;
	util::Timer_tic(&solveTimer);
	
	// B&U optimized in-place
	double cost = beliefCollocation(B, U, problem, output, info);

	double solvetime = util::Timer_toc(&solveTimer);
	LOG_INFO("Optimized cost: %4.10f", cost);
	LOG_INFO("Solve time: %5.3f ms", solvetime*1000);
	
	cleanupBeliefMPCVars();

	//pythonDisplayTrajectory(B, U);

	/*
	for (size_t t = 0; t < T; ++t) {
		std::cout << ~B[t] << std::endl;
	}
	*/

	//int k;
	//std::cin >> k;

	//CAL_End();
	return 0;
}
