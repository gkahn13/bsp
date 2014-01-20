#ifndef __POINT_H__
#define __POINT_H__

#include <fstream>

#include "util/matrix.h"
//extern "C" {
#include "util/utils.h"
//}
//#include "util/Timer.h"
#include "util/logging.h"

#include <Python.h>
//#include <pythonrun.h>
#include <boost/python.hpp>
#include <boost/filesystem.hpp>

namespace py = boost::python;

#define TIMESTEPS 15
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

const int T = TIMESTEPS;
const double INFTY = 1e10;

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

const double alpha_belief = 10, alpha_final_belief = 10, alpha_control = 1, alpha_goal_state = 1;

double *inputVars, *vars;
std::vector<int> maskIndices;


Matrix<X_DIM> dynfunc(const Matrix<X_DIM>& x, const Matrix<U_DIM>& u, const Matrix<U_DIM>& q)
{
	Matrix<X_DIM> xNew = x + u*DT + 0.01*q;
	return xNew;
}

// Observation model
Matrix<Z_DIM> obsfunc(const Matrix<X_DIM>& x, const Matrix<R_DIM>& r)
{
	double intensity = sqrt(sqr(0.5*x[0]) + 1e-6);
	Matrix<Z_DIM> z = x + intensity*r;
	return z;
}

// Jacobians: df(x,u,q)/dx, df(x,u,q)/dq
void linearizeDynamics(const Matrix<X_DIM>& x, const Matrix<U_DIM>& u, const Matrix<Q_DIM>& q, Matrix<X_DIM,X_DIM>& A, Matrix<X_DIM,Q_DIM>& M)
{
	A.reset();
	Matrix<X_DIM> xr(x), xl(x);
	for (size_t i = 0; i < X_DIM; ++i) {
		xr[i] += step; xl[i] -= step;
		A.insert(0,i, (dynfunc(xr, u, q) - dynfunc(xl, u, q)) / (xr[i] - xl[i]));
		xr[i] = x[i]; xl[i] = x[i];
	}

	M.reset();
	Matrix<Q_DIM> qr(q), ql(q);
	for (size_t i = 0; i < Q_DIM; ++i) {
		qr[i] += step; ql[i] -= step;
		M.insert(0,i, (dynfunc(x, u, qr) - dynfunc(x, u, ql)) / (qr[i] - ql[i]));
		qr[i] = q[i]; ql[i] = q[i];
	}
}

// Jacobians: dh(x,r)/dx, dh(x,r)/dr
void linearizeObservation(const Matrix<X_DIM>& x, const Matrix<R_DIM>& r, Matrix<Z_DIM,X_DIM>& H, Matrix<Z_DIM,R_DIM>& N)
{
	H.reset();
	Matrix<X_DIM> xr(x), xl(x);
	for (size_t i = 0; i < X_DIM; ++i) {
		xr[i] += step; xl[i] -= step;
		H.insert(0,i, (obsfunc(xr, r) - obsfunc(xl, r)) / (xr[i] - xl[i]));
		xr[i] = x[i]; xl[i] = x[i];
	}

	N.reset();
	Matrix<R_DIM> rr(r), rl(r);
	for (size_t i = 0; i < R_DIM; ++i) {
		rr[i] += step; rl[i] -= step;
		N.insert(0,i, (obsfunc(x, rr) - obsfunc(x, rl)) / (rr[i] - rl[i]));
		rr[i] = r[i]; rl[i] = r[i];
	}
}

// Switch between belief vector and matrices
void unVec(const Matrix<B_DIM>& b, Matrix<X_DIM>& x, Matrix<X_DIM,X_DIM>& SqrtSigma) {
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

void vec(const Matrix<X_DIM>& x, const Matrix<X_DIM,X_DIM>& SqrtSigma, Matrix<B_DIM>& b) {
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
Matrix<B_DIM> beliefDynamics(const Matrix<B_DIM>& b, const Matrix<U_DIM>& u) {
	Matrix<X_DIM> x;
	Matrix<X_DIM,X_DIM> SqrtSigma;
	unVec(b, x, SqrtSigma);

	Matrix<X_DIM,X_DIM> Sigma = SqrtSigma*SqrtSigma;

	Matrix<X_DIM,X_DIM> A;
	Matrix<X_DIM,Q_DIM> M;
	linearizeDynamics(x, u, zeros<Q_DIM,1>(), A, M);

	x = dynfunc(x, u, zeros<Q_DIM,1>());
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


void setupDstarInterface(std::string mask) {
	std::stringstream ss(mask);
	int val, i=0;
	while (ss >> val) {
		if (val == 1) {
			maskIndices.push_back(i);
		}
		i++;
	}

	inputVars = new double[i];
}

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

	py::list x0_list, xGoal_list;
	for(int i=0; i < X_DIM; i++) {
		x0_list.append(x0[i]);
		xGoal_list.append(xGoal[i]);
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
		py::object plot_traj = plot_mod.attr("plot_belief_trajectory_cpp");

		plot_traj(Bvec, Uvec, model, x0_list, xGoal_list, T);
	}
	catch(py::error_already_set const &)
	{
		PyErr_Print();
	}

}

#endif
