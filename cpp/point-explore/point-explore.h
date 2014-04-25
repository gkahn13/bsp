#ifndef __POINT_EXPLORE_H__
#define __POINT_EXPLORE_H__

#include <Python.h>

#include <vector>
#include <cmath>

#include "../util/matrix.h"
#include "../util/utils.h"
#include "../util/logging.h"

#include <boost/python.hpp>
#include <boost/filesystem.hpp>

#include <symbolic/casadi.hpp>
#include <symbolic/stl_vector_tools.hpp>
#include <cstdlib>

namespace py = boost::python;

#define TIMESTEPS 10
#define PARTICLES 100
#define AGENTS 2
#define DT 1.0 // Note: if you change this, must change the FORCES matlab file
#define X_DIM 2
#define U_DIM 2
#define Z_DIM 1
#define Q_DIM 2
#define R_DIM 1

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

const int T = TIMESTEPS;
const int M = PARTICLES;
const int N = AGENTS;
const int TOTAL_VARS = T*N*X_DIM + (T-1)*N*U_DIM;

const double step = 0.0078125*0.0078125;
const double INFTY = 1e10;

const double alpha = 2.5;
const double max_range = 0.25;

const double alpha_control_norm = 0; // 1e-2
const double alpha_control_smooth = 0; // 1e-2
const double alpha_separation = 0; // 1e-3



extern SymmetricMatrix<N*R_DIM> R;

extern Matrix<N*X_DIM> x0;
extern Matrix<X_DIM> xMin, xMax;
extern Matrix<U_DIM> uMin, uMax;

extern Matrix<X_DIM> target;

#include "casadi/casadi-point-explore.h"
namespace AD = CasADi;
extern AD::SXFunction casadi_differential_entropy_func;
extern AD::SXFunction casadi_grad_differential_entropy_func;
extern AD::SXFunction casadi_diaghess_differential_entropy_func;

inline double uniform(double low, double high) {
	return (high - low)*(rand() / double(RAND_MAX)) + low;
}

template<int _dim>
inline double dist(const Matrix<_dim>& a, const Matrix<_dim>& b) {
	return sqrt(tr(~(a-b)*(a-b)));
}


template <class C>
inline void subsample(const std::vector<C>& full, std::vector<C>& sub) {
	std::shuffle(full.begin(), full.end(), std::default_random_engine(time(0)));
	sub = std::vector<C>(full.begin(), full.begin() + sub.size());
}


namespace point_explore {


void initialize();

Matrix<N*X_DIM> dynfunc(const Matrix<N*X_DIM>& x, const Matrix<N*U_DIM>& u);

Matrix<N*Z_DIM> obsfunc(const Matrix<N*X_DIM>& x, const Matrix<X_DIM>& t, const Matrix<N*R_DIM>& r);

template <int _vDim, int _sDim>
double gaussLikelihood(const Matrix<_vDim>& v, const SymmetricMatrix<_sDim>& S);

std::vector<Matrix<X_DIM> > lowVarianceSampler(const std::vector<Matrix<X_DIM> >& P, const std::vector<double>& W, double r);


void updateStateAndParticles(const Matrix<N*X_DIM>& x_t, const std::vector<Matrix<X_DIM> >& P_t, const Matrix<N*U_DIM>& u_t,
								Matrix<N*X_DIM>& x_tp1, std::vector<Matrix<X_DIM> >& P_tp1);


double differential_entropy(const std::vector<Matrix<N*X_DIM> >& X, const std::vector<Matrix<N*U_DIM> >& U,
					 const std::vector<Matrix<X_DIM> >& P);

double casadi_differential_entropy(const std::vector<Matrix<N*X_DIM> >& X, const std::vector<Matrix<N*U_DIM> >& U,
					 const std::vector<Matrix<X_DIM> >& P);

Matrix<TOTAL_VARS> grad_differential_entropy(std::vector<Matrix<N*X_DIM> >& X, std::vector<Matrix<N*U_DIM> >& U,
												const std::vector<Matrix<X_DIM> >& P);

Matrix<TOTAL_VARS> casadi_grad_differential_entropy(const std::vector<Matrix<N*X_DIM> >& X, const std::vector<Matrix<N*U_DIM> >& U,
					 const std::vector<Matrix<X_DIM> >& P);

//Matrix<TOTAL_VARS> diaghess_differential_entropy(std::vector<Matrix<X_DIM> >& X, std::vector<Matrix<U_DIM> >& U,
//												const std::vector<Matrix<X_DIM> >& P);
//
//Matrix<TOTAL_VARS> casadi_diaghess_differential_entropy(const std::vector<Matrix<X_DIM> >& X, const std::vector<Matrix<U_DIM> >& U,
//					 const std::vector<Matrix<X_DIM> >& P);



void pythonDisplayStatesAndParticles(const std::vector<Matrix<N*X_DIM>>& X, const std::vector<Matrix<X_DIM> >& P, const Matrix<X_DIM>& targ);

}


#endif
