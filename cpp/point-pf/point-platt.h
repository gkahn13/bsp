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

#define TIMESTEPS 15
#define PARTICLES 5
#define DT 1.0
#define X_DIM 2
#define U_DIM 2
#define Z_DIM 2
#define Q_DIM 2
#define R_DIM 2

const int T = TIMESTEPS;
const int M = PARTICLES;
const int TOTAL_VARS = T*M*X_DIM + (T-1)*U_DIM;

const double alpha_control = 1;

SymmetricMatrix<Q_DIM> Q;
SymmetricMatrix<R_DIM> R;

std::vector<Matrix<X_DIM>> P0(M);
Matrix<X_DIM> xGoal;
Matrix<X_DIM> xMin, xMax;
Matrix<U_DIM> uMin, uMax;

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

const double step = 0.0078125*0.0078125;

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

Matrix<TOTAL_VARS> deriv_costfunc(const std::vector<std::vector<Matrix<X_DIM> > >& P, const std::vector<Matrix<U_DIM> >& U) {
	Matrix<TOTAL_VARS> d;

	float orig, cost_p, cost_l;
	int index = 0;
	for(int t=0; t < T; ++t) {
		for(int m=0; m < M; ++m) {
			for(int i=0; i < X_DIM; ++i) {
				orig = P[t][m][i];

				P[t][m][i] = orig + step;
				cost_p = costfunc(P,U);

				P[t][m][i] = orig - step;
				cost_l = costfunc(P,U);

				P[t][m][i] = orig;
				d[index++] = (cost_p - cost_l)/(2*step);
			}
		}

		if (t < T-1) {
			for(int i=0; i < U_DIM; ++i) {
				orig = U[t][i];

				U[t][i] = orig + step;
				cost_p = costfunc(P,U);

				U[t][i] = orig - step;
				cost_l = costfunc(P,U);

				U[t][i] = orig;
				d[index++] = (cost_p - cost_l)/(2*step);
			}
		}
	}

	return d;
}

}

#endif
