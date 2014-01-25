#ifndef __POINT_H__
#define __POINT_H__

#include <fstream>
#include <tgmath.h>

#include "util/matrix.h"


#include "util/utils.h"
#include "util/logging.h"

#include <Python.h>
//#include <pythonrun.h>
#include <boost/python.hpp>
#include <boost/filesystem.hpp>

namespace py = boost::python;

#define TIMESTEPS 10
#define DT 1
#define X_DIM 6
#define U_DIM 6
#define Z_DIM 4
#define Q_DIM 6
#define R_DIM 4
#define G_DIM 3 // dimension of g(x) function, which is 3d position

#define S_DIM (((X_DIM+1)*X_DIM)/2)
#define B_DIM (X_DIM+S_DIM)
#define XU_DIM (X_DIM*T+U_DIM*(T-1))

const double l4 = 2.375;
const double l3 = 10.375;
const double l2 = 8;
const double l1 = 7.25;

const double step = 0.0078125*0.0078125;

Matrix<X_DIM> x0;
Matrix<X_DIM,X_DIM> SqrtSigma0;
Matrix<G_DIM> posGoal;
Matrix<X_DIM> xGoal; // TODO: temporary, since goal is a vector of joints
Matrix<X_DIM> xMin, xMax;
Matrix<U_DIM> uMin, uMax;

Matrix<G_DIM> cam0, cam1;


const int T = TIMESTEPS;
const double INFTY = 1e10;

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

const double alpha_belief = 10, alpha_final_belief = 10, alpha_control = 1, alpha_goal_state = 1;

//double *inputVars, *vars;
//std::vector<int> maskIndices;


Matrix<X_DIM> dynfunc(const Matrix<X_DIM>& x, const Matrix<U_DIM>& u, const Matrix<U_DIM>& q)
{
	Matrix<X_DIM> xNew = x + (u + q)*DT;
	return xNew;
}

// joint angles -> end effector position
Matrix<G_DIM> g(const Matrix<X_DIM>& x)
{
    double a0 = x[0], a1 = x[1], a2 = x[2], a3 = x[3], a4 = x[4], a5 = x[5];

    Matrix<G_DIM> p;
    p[0] = sin(a0) * (cos(a1) * (sin(a2) * (cos(a4) * l4 + l3) + cos(a2) * cos(a3) * sin(a4) * l4) + sin(a1) * (cos(a2) * (cos(a4) * l4 + l3) - sin(a2) * cos(a3) * sin(a4) * l4 + l2)) + cos(a0) * sin(a3) * sin(a4) * l4;
    p[1] = -sin(a1) * (sin(a2) * (cos(a4) * l4 + l3) + cos(a2) * cos(a3) * sin(a4) * l4) + cos(a1) * (cos(a2) * (cos(a4) * l4 + l3) - sin(a2) * cos(a3) * sin(a4) * l4 + l2) + l1;
    p[2] = cos(a0) * (cos(a1) * (sin(a2) * (cos(a4) * l4 + l3) + cos(a2) * cos(a3) * sin(a4) * l4) + sin(a1) * (cos(a2) * (cos(a4) * l4 + l3) - sin(a2) * cos(a3) * sin(a4) * l4 + l2)) - sin(a0) * sin(a3) * sin(a4) * l4;

    return p;
}


// Observation model
Matrix<Z_DIM> obsfunc(const Matrix<X_DIM>& x, const Matrix<R_DIM>& r)
{
	Matrix<G_DIM> ee_pos = g(x);

	Matrix<Z_DIM> obs;
	obs[0] = (ee_pos[0] - cam0[0]) / (ee_pos[1] - cam0[1]) + r[0];
	obs[1] = (ee_pos[2] - cam0[2]) / (ee_pos[1] - cam0[1]) + r[1];
	obs[2] = (ee_pos[0] - cam1[0]) / (ee_pos[1] - cam1[1]) + r[2];
    obs[3] = (ee_pos[2] - cam1[2]) / (ee_pos[1] - cam1[1]) + r[3];

    return obs;
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

	Matrix<X_DIM,X_DIM> A = identity<X_DIM>();
	Matrix<X_DIM,Q_DIM> M = DT*identity<U_DIM>();
	//linearizeDynamics(x, u, zeros<Q_DIM,1>(), A, M);

	x = dynfunc(x, u, zeros<Q_DIM,1>());
	Sigma = A*Sigma*~A + M*~M;

	// TODO: fill in with symbolically computed linearizations
	Matrix<Z_DIM,X_DIM> H = zeros<Z_DIM,X_DIM>();
	Matrix<Z_DIM,R_DIM> N = zeros<Z_DIM,R_DIM>();
	linearizeObservation(x, zeros<R_DIM>(), H, N);

	// TODO: currently fails on second 7th sqp iteration b/c H is huge
	Matrix<X_DIM,Z_DIM> K = Sigma*~H/(H*Sigma*~H + N*~N);

	Sigma = (identity<X_DIM>() - K*H)*Sigma;

	Matrix<B_DIM> g;
	vec(x, sqrt(Sigma), g);

	return g;
}



#endif
