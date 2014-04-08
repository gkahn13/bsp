#ifndef __POINT_PLATT_H__
#define __POINT_PLATT_H__

#include <vector>
#include <cmath>

#include "../util/matrix.h"
#include "../util/utils.h"
#include "../util/logging.h"

#include <Python.h>
#include <boost/python.hpp>
#include <boost/filesystem.hpp>

namespace py = boost::python;

#define TIMESTEPS 3
#define DT 1.0
#define X_DIM 2
#define U_DIM 2
#define Z_DIM 2
#define Q_DIM 2
#define R_DIM 2

const int T = TIMESTEPS;

const double alpha_control = 1;

SymmetricMatrix<Q_DIM> Q;
SymmetricMatrix<R_DIM> R;

namespace point_platt {

Matrix<X_DIM> dynfunc(const Matrix<X_DIM>& x, const Matrix<U_DIM>& u, const Matrix<Q_DIM>& q)
{
	Matrix<X_DIM> xNew = x + u*DT + .01*q;
	return xNew;
}

Matrix<Z_DIM> obsfunc(const Matrix<X_DIM>& x, const Matrix<R_DIM>& r)
{
	Matrix<Z_DIM> z = x + sqrt(0.5*0.5*x[0]*x[0] + 1e-6)*r;
	return z;
}

// assume P[t][0] is most likely particle
float costfunc(const std::vector<std::vector<Matrix<X_DIM> > >& P, const std::vector<Matrix<U_DIM> >& U) {
	float cost = 0;
	int M = P[0].size();

	Matrix<X_DIM> P_t_mle, P_t_other;
	Matrix<Z_DIM> h_t_mle, h_t_other, diff;
	for(int t=1; t < T; ++t) {
		P_t_mle = P[t][0];
		h_t_mle = obsfunc(P_t_mle, zeros<R_DIM,1>());
		for(int m=1; m < M; ++m) {
			P_t_other = P[t][m];
			h_t_other = obsfunc(P_t_other, zeros<R_DIM,1>());
			diff = h_t_mle - h_t_other;
			cost += tr(~diff*diff);
		}
	}

	for(int t=0; t < T-1; ++t) {
		cost += alpha_control*tr(~U[t]*U[t]);
	}

	return cost;
}

}

#endif
