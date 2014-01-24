#ifndef __POINT_DOOUBLE_H__
#define __POINT_DOUBLE_H__

#include <fstream>
#include <tgmath.h>

/*
#ifndef USE_ADOLC
#include "util/dMatrix.h"
#else
#include <adolc/adolc.h>
#include "util/adMatrix.h"
#endif
*/
//#include <adolc/adolc.h>
#include "util/dmatrix.h"


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

/*
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

dMatrix<X_DIM> x0;
dMatrix<X_DIM,X_DIM> SqrtSigma0;
dMatrix<X_DIM> xGoal;
dMatrix<X_DIM> xMin, xMax;
dMatrix<U_DIM> uMin, uMax;

const int T = TIMESTEPS;
const double INFTY = 1e10;

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

const double alpha_belief = 10, alpha_final_belief = 10, alpha_control = 1, alpha_goal_state = 1;

double *inputVars, *vars;
std::vector<int> maskIndices;
*/


dMatrix<X_DIM> dynfunc_d(const dMatrix<X_DIM>& x, const dMatrix<U_DIM>& u, const dMatrix<U_DIM>& q)
{
	dMatrix<X_DIM> xNew = x + u*DT + 0.01*q;
	return xNew;
}

// Observation model
dMatrix<Z_DIM> obsfunc_d(const dMatrix<X_DIM>& x, const dMatrix<R_DIM>& r)
{
	//double intensity = sqrt(0.5*0.5*x[0]*x[0] + 1e-6);
	dMatrix<Z_DIM> z = x + sqrt(0.5*0.5*x[0]*x[0] + 1e-6)*r;
	return z;
}

// Jacobians: df(x,u,q)/dx, df(x,u,q)/dq
void linearizeDynamics_d(const dMatrix<X_DIM>& x, const dMatrix<U_DIM>& u, const dMatrix<Q_DIM>& q, dMatrix<X_DIM,X_DIM>& A, dMatrix<X_DIM,Q_DIM>& M)
{
	A.reset();
	dMatrix<X_DIM> xr(x), xl(x);
	for (size_t i = 0; i < X_DIM; ++i) {
		xr[i] += step; xl[i] -= step;
		A.insert(0,i, (dynfunc_d(xr, u, q) - dynfunc_d(xl, u, q)) / (xr[i] - xl[i]));
		xr[i] = x[i]; xl[i] = x[i];
	}

	M.reset();
	dMatrix<Q_DIM> qr(q), ql(q);
	for (size_t i = 0; i < Q_DIM; ++i) {
		qr[i] += step; ql[i] -= step;
		M.insert(0,i, (dynfunc_d(x, u, qr) - dynfunc_d(x, u, ql)) / (qr[i] - ql[i]));
		qr[i] = q[i]; ql[i] = q[i];
	}
}

// Jacobians: dh(x,r)/dx, dh(x,r)/dr
void linearizeObservation_d(const dMatrix<X_DIM>& x, const dMatrix<R_DIM>& r, dMatrix<Z_DIM,X_DIM>& H, dMatrix<Z_DIM,R_DIM>& N)
{
	H.reset();
	dMatrix<X_DIM> xr(x), xl(x);
	for (size_t i = 0; i < X_DIM; ++i) {
		xr[i] += step; xl[i] -= step;
		H.insert(0,i, (obsfunc_d(xr, r) - obsfunc_d(xl, r)) / (xr[i] - xl[i]));
		xr[i] = x[i]; xl[i] = x[i];
	}

	N.reset();
	dMatrix<R_DIM> rr(r), rl(r);
	for (size_t i = 0; i < R_DIM; ++i) {
		rr[i] += step; rl[i] -= step;
		N.insert(0,i, (obsfunc_d(x, rr) - obsfunc_d(x, rl)) / (rr[i] - rl[i]));
		rr[i] = r[i]; rl[i] = r[i];
	}
}

// Switch between belief vec_dtor and matrices
void unVec_d(const dMatrix<B_DIM>& b, dMatrix<X_DIM>& x, dMatrix<X_DIM,X_DIM>& SqrtSigma) {
	x = b.subdMatrix<X_DIM,1>(0,0);
	size_t idx = X_DIM;
	for (size_t j = 0; j < X_DIM; ++j) {
		for (size_t i = j; i < X_DIM; ++i) {
			SqrtSigma(i,j) = b[idx];
			SqrtSigma(j,i) = b[idx];
			++idx;
		}
	}
}

void vec_d(const dMatrix<X_DIM>& x, const dMatrix<X_DIM,X_DIM>& SqrtSigma, dMatrix<B_DIM>& b) {
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
dMatrix<B_DIM> beliefDynamics_d(const dMatrix<B_DIM>& b, const dMatrix<U_DIM>& u) {
	dMatrix<X_DIM> x;
	dMatrix<X_DIM,X_DIM> SqrtSigma;
	unVec_d(b, x, SqrtSigma);

	dMatrix<X_DIM,X_DIM> Sigma = SqrtSigma*SqrtSigma;

	dMatrix<X_DIM,X_DIM> A = dIdentity<X_DIM>();
	dMatrix<X_DIM,Q_DIM> M = .01*dIdentity<U_DIM>();
	//linearizeDynamics_d(x, u, zeros<Q_DIM,1>(), A, M);

	x = dynfunc_d(x, u, dZeros<Q_DIM,1>());
	Sigma = A*Sigma*~A + M*~M;

	dMatrix<Z_DIM,X_DIM> H = dZeros<Z_DIM,X_DIM>();
	dMatrix<Z_DIM,R_DIM> N = dZeros<Z_DIM,R_DIM>();
	H(0,0) = 1; H(1,1) = 1;
	N(0,0) = sqrt(x(0,0) * x(0,0) * 0.5 * 0.5 + 1e-6);
	N(1,1) = sqrt(x(0,0) * x(0,0) * 0.5 * 0.5 + 1e-6);
	//linearizeObservation_d(x, zeros<R_DIM,1>(), H, N);

	dMatrix<X_DIM,Z_DIM> K = Sigma*~H/(H*Sigma*~H + N*~N);

	Sigma = (dIdentity<X_DIM>() - K*H)*Sigma;

	dMatrix<B_DIM> g;
	vec_d(x, sqrt(Sigma), g);

	return g;
}


void linearizeObservation_d(std::string mask) {
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

void pythonDisplayTrajectory_d(std::vector< dMatrix<B_DIM> >& B, std::vector< dMatrix<U_DIM> >& U)
{
	for (int t = 0; t < T-1; ++t) {
		B[t+1] = beliefDynamics_d(B[t], U[t]);
	}

	py::list Bvec_d;
	for(int j=0; j < B_DIM; j++) {
		for(int i=0; i < T; i++) {
			Bvec_d.append(B[i][j]);
		}
	}

	py::list Uvec_d;
		for(int j=0; j < U_DIM; j++) {
			for(int i=0; i < T-1; i++) {
			Uvec_d.append(U[i][j]);
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

		plot_traj(Bvec_d, Uvec_d, model, x0_list, xGoal_list, T);
	}
	catch(py::error_already_set const &)
	{
		PyErr_Print();
	}

}

#endif
