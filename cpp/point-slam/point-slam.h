#ifndef __POINT_H__
#define __POINT_H__

#include <fstream>

#include "util/matrix.h"


#include "util/utils.h"
#include "util/logging.h"

#include <Python.h>
//#include <pythonrun.h>
#include <boost/python.hpp>
#include <boost/filesystem.hpp>

namespace py = boost::python;

#define TIMESTEPS 15
#define DT 1.0

#define NUM_LANDMARKS 4
#define NUM_WAYPOINTS 5

#define P_DIM 2 // robot state size (x, y)
#define L_DIM (2*NUM_LANDMARKS) // landmark state size (x, y)

#define X_DIM (P_DIM + L_DIM)
#define U_DIM P_DIM
#define Z_DIM L_DIM
#define Q_DIM P_DIM
#define R_DIM L_DIM

#define S_DIM (((X_DIM+1)*X_DIM)/2)
#define B_DIM (X_DIM+S_DIM)
#define XU_DIM (X_DIM*T+U_DIM*(T-1))

const double step = 0.0078125*0.0078125;
const double range = 5;

const int T = TIMESTEPS;
const double INFTY = 1e10;

#define MIN(a,b) (*((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

const double alpha_belief = 10, alpha_final_belief = 10, alpha_control = 1, alpha_goal_state = 1;

const double initial_sigma_factor = 1;

double *inputVars, *vars;
std::vector<int> maskIndices;


Matrix<X_DIM> dynfunc(const Matrix<X_DIM>& x, const Matrix<U_DIM>& u, const Matrix<U_DIM>& q)
{
	Matrix<P_DIM> pNew;
	Matrix<X_DIM> xNew;

	pNew = x.subMatrix<P_DIM>(0,0) + u*DT + q;
	xNew.insert(0, 0, pNew);
	xNew.insert(P_DIM, 0, x.subMatrix<L_DIM>(P_DIM,0));

	/*
	Matrix<L_DIM> l = x.subMatrix<L_DIM>(P_DIM,0);
	double x0 = x[0], x1 = x[1], l0, l1, dist;

	for(int i = 0; i < L_DIM; i += 2) {
		l0 = l[i];
		l1 = l[i+1];

		dist = sqrt((x0 - l0)*(x0 - l0) + (x1 - l1)*(x1 - l1));

		if (dist > range) {
			xNew[P_DIM+i] += q[0];
			xNew[P_DIM+i+1] += q[1];
		}
	}
	*/

	return xNew;
}

// Observation model
Matrix<Z_DIM> obsfunc(const Matrix<X_DIM>& x, const Matrix<R_DIM>& r)
{
	Matrix<Z_DIM> z;
	Matrix<L_DIM> l = x.subMatrix<L_DIM>(P_DIM,0);
	double x0 = x[0], x1 = x[1], l0, l1, dist;

	double alpha = 10;
	for(int i = 0; i < L_DIM; i += 2) {
		l0 = l[i];
		l1 = l[i+1];

		dist = sqrt((x0 - l0)*(x0 - l0) + (x1 - l1)*(x1 - l1));

		z[i] = (x0 - l0) + (8*abs(range-dist)*r[i])/ (1+exp(alpha*(range-dist)));
		z[i+1] = (x1 - l1) + (8*abs(range-dist)*r[i+1])/ (1+exp(alpha*(range-dist)));

		std::cout << 1/(1+exp(alpha*(range-dist))) << " " << std::endl;
	}

	//exit(0);

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
void unVec(const Matrix<B_DIM>& b, Matrix<X_DIM>& x, Matrix<X_DIM,X_DIM>& S) {
	x = b.subMatrix<X_DIM,1>(0,0);
	size_t idx = X_DIM;
	for (size_t j = 0; j < X_DIM; ++j) {
		for (size_t i = j; i < X_DIM; ++i) {
			S(i,j) = b[idx];
			S(j,i) = b[idx];
			++idx;
		}
	}
}

void vec(const Matrix<X_DIM>& x, const Matrix<X_DIM,X_DIM>& S, Matrix<B_DIM>& b) {
	b.insert(0,0,x);
	size_t idx = X_DIM;
	for (size_t j = 0; j < X_DIM; ++j) {
		for (size_t i = j; i < X_DIM; ++i) {
			b[idx] = 0.5 * (S(i,j) + S(j,i));
			++idx;
		}
	}
}


// Belief dynamics
Matrix<B_DIM> beliefDynamics(const Matrix<B_DIM>& b, const Matrix<U_DIM>& u) {
	Matrix<X_DIM> x;
	Matrix<X_DIM,X_DIM> Sigma;
	unVec(b, x, Sigma);

	Sigma = Sigma*Sigma;

	//Matrix<X_DIM,X_DIM> A = identity<X_DIM>();
	//Matrix<X_DIM,Q_DIM> M = .01*identity<U_DIM>();
	Matrix<X_DIM,X_DIM> A;
	Matrix<X_DIM,Q_DIM> M;

	Matrix<Q_DIM,Q_DIM> QC = identity<Q_DIM>();

	linearizeDynamics(x, u, zeros<Q_DIM,1>(), A, M);
	x = dynfunc(x, u, zeros<Q_DIM,1>());

	Sigma = A*Sigma*~A + M*QC*~M;

	Matrix<Z_DIM,X_DIM> H = zeros<Z_DIM,X_DIM>();
	Matrix<Z_DIM,R_DIM> N = zeros<Z_DIM,R_DIM>();
	linearizeObservation(x, zeros<R_DIM,1>(), H, N);

	Matrix<R_DIM,R_DIM> RC = 10*identity<R_DIM>();

	Matrix<X_DIM,Z_DIM> K = Sigma*~H/(H*Sigma*~H + N*RC*~N);

	std::cout << "I - KH" << std::endl << identity<X_DIM>() - K*H << std::endl;

	Sigma = (identity<X_DIM>() - K*H)*Sigma;
	
	Matrix<B_DIM> g;
	vec(x, sqrt(Sigma), g);

	return g;
}


void pythonDisplayTrajectory(std::vector< Matrix<B_DIM> >& B, std::vector< Matrix<P_DIM> >& waypoints, int time_steps)
{

	// B_vec is only for the robot, not the landmarks
	py::list B_vec;
	for(int j=0; j < P_DIM; j++) {
		for(int i=0; i < time_steps; i++) {
			B_vec.append(B[i][j]);
		}
	}

	Matrix<X_DIM> x;
	Matrix<X_DIM,X_DIM> SqrtSigma;
	for(int i=0; i < time_steps; i++) {
		unVec(B[i], x, SqrtSigma);
		B_vec.append(SqrtSigma[0,0]);
	}
	for(int i=0; i < time_steps; i++) {
		unVec(B[i], x, SqrtSigma);
		B_vec.append(SqrtSigma[1,0]);
	}
	for(int i=0; i < time_steps; i++) {
		unVec(B[i], x, SqrtSigma);
		B_vec.append(SqrtSigma[1,1]);
	}

	py::list waypoints_vec;
	for(int j=0; j < 2; j++) {
		for(int i=0; i < NUM_WAYPOINTS; i++) {
			waypoints_vec.append(waypoints[i][j]);
		}
	}

	std::string workingDir = boost::filesystem::current_path().normalize().string();

	try
	{
		Py_Initialize();
		py::object main_module = py::import("__main__");
		py::object main_namespace = main_module.attr("__dict__");
		py::exec("import sys, os", main_namespace);
		py::exec(py::str("sys.path.append('"+workingDir+"/point-slam')"), main_namespace);
		py::object plot_mod = py::import("plot_point_slam");
		py::object plot_traj = plot_mod.attr("plot_point_trajectory");

		plot_traj(B_vec, waypoints_vec, time_steps);
	}
	catch(py::error_already_set const &)
	{
		PyErr_Print();
	}

}

#endif
