#ifndef __POINT_ADOUBLE_H__
#define __POINT_ADOUBLE_H__

#include <fstream>
#include <tgmath.h>

#include <adolc/adolc.h>

#include "util/amatrix.h"
#include "point.h"

#include "util/utils.h"
#include "util/logging.h"

#include <Python.h>
//#include <pythonrun.h>
#include <boost/python.hpp>
//#include <boost/filesystem.hpp>

namespace py = boost::python;


aMatrix<X_DIM> dynfunc_a(const aMatrix<X_DIM>& x, const aMatrix<U_DIM>& u, const aMatrix<U_DIM>& q)
{
	aMatrix<X_DIM> xNew = x + u*DT + 0.01*q;
	return xNew;
}

// Observation model
aMatrix<Z_DIM> obsfunc_a(const aMatrix<X_DIM>& x, const aMatrix<R_DIM>& r)
{
	//double intensity = sqrt(0.5*0.5*x[0]*x[0] + 1e-6);
	aMatrix<Z_DIM> z = x + sqrt(0.5*0.5*x[0]*x[0] + 1e-6)*r;
	return z;
}

// Jacobians: df(x,u,q)/dx, df(x,u,q)/dq
void linearizeDynamics_a(const aMatrix<X_DIM>& x, const aMatrix<U_DIM>& u, const aMatrix<Q_DIM>& q, aMatrix<X_DIM,X_DIM>& A, aMatrix<X_DIM,Q_DIM>& M)
{
	A.reset();
	aMatrix<X_DIM> xr(x), xl(x);
	for (size_t i = 0; i < X_DIM; ++i) {
		xr[i] += step; xl[i] -= step;
		A.insert(0,i, (dynfunc_a(xr, u, q) - dynfunc_a(xl, u, q)) / (xr[i] - xl[i]));
		xr[i] = x[i]; xl[i] = x[i];
	}

	M.reset();
	aMatrix<Q_DIM> qr(q), ql(q);
	for (size_t i = 0; i < Q_DIM; ++i) {
		qr[i] += step; ql[i] -= step;
		M.insert(0,i, (dynfunc_a(x, u, qr) - dynfunc_a(x, u, ql)) / (qr[i] - ql[i]));
		qr[i] = q[i]; ql[i] = q[i];
	}
}

// Jacobians: dh(x,r)/dx, dh(x,r)/dr
void linearizeObservation_a(const aMatrix<X_DIM>& x, const aMatrix<R_DIM>& r, aMatrix<Z_DIM,X_DIM>& H, aMatrix<Z_DIM,R_DIM>& N)
{
	H.reset();
	aMatrix<X_DIM> xr(x), xl(x);
	for (size_t i = 0; i < X_DIM; ++i) {
		xr[i] += step; xl[i] -= step;
		H.insert(0,i, (obsfunc_a(xr, r) - obsfunc_a(xl, r)) / (xr[i] - xl[i]));
		xr[i] = x[i]; xl[i] = x[i];
	}

	N.reset();
	aMatrix<R_DIM> rr(r), rl(r);
	for (size_t i = 0; i < R_DIM; ++i) {
		rr[i] += step; rl[i] -= step;
		N.insert(0,i, (obsfunc_a(x, rr) - obsfunc_a(x, rl)) / (rr[i] - rl[i]));
		rr[i] = r[i]; rl[i] = r[i];
	}
}

// Switch between belief vec_ator and matrices
void unVec_a(const aMatrix<B_DIM>& b, aMatrix<X_DIM>& x, aMatrix<X_DIM,X_DIM>& SqrtSigma) {
	x = b.subaMatrix<X_DIM,1>(0,0);
	size_t idx = X_DIM;
	for (size_t j = 0; j < X_DIM; ++j) {
		for (size_t i = j; i < X_DIM; ++i) {
			SqrtSigma(i,j) = b[idx];
			SqrtSigma(j,i) = b[idx];
			++idx;
		}
	}
}

void vec_a(const aMatrix<X_DIM>& x, const aMatrix<X_DIM,X_DIM>& SqrtSigma, aMatrix<B_DIM>& b) {
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
aMatrix<B_DIM> beliefDynamics_a(const aMatrix<B_DIM>& b, const aMatrix<U_DIM>& u) {
	aMatrix<X_DIM> x;
	aMatrix<X_DIM,X_DIM> SqrtSigma;
	unVec_a(b, x, SqrtSigma);

	aMatrix<X_DIM,X_DIM> Sigma = SqrtSigma*SqrtSigma;

	aMatrix<X_DIM,X_DIM> A = aIdentity<X_DIM>();
	aMatrix<X_DIM,Q_DIM> M = .01*aIdentity<U_DIM>();
	//linearizeDynamics_a(x, u, zeros<Q_DIM,1>(), A, M);

	x = dynfunc_a(x, u, aZeros<Q_DIM,1>());
	Sigma = A*Sigma*~A + M*~M;

	aMatrix<Z_DIM,X_DIM> H = aZeros<Z_DIM,X_DIM>();
	aMatrix<Z_DIM,R_DIM> N = aZeros<Z_DIM,R_DIM>();
	H(0,0) = 1; H(1,1) = 1;
	N(0,0) = sqrt(x(0,0) * x(0,0) * 0.5 * 0.5 + 1e-6);
	N(1,1) = sqrt(x(0,0) * x(0,0) * 0.5 * 0.5 + 1e-6);
	//linearizeObservation_a(x, zeros<R_DIM,1>(), H, N);

	aMatrix<X_DIM,Z_DIM> K = Sigma*~H/(H*Sigma*~H + N*~N);

	Sigma = (aIdentity<X_DIM>() - K*H)*Sigma;

	aMatrix<B_DIM> g;
	vec_a(x, sqrt(Sigma), g);

	return g;
}



#endif
