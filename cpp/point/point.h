#ifndef __POINT_H__
#define __POINT_H__

#include "util/matrix.h"
extern "C" {
#include "util/utils.h"
}
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

extern const double step;

extern Matrix<X_DIM> x0;
extern Matrix<X_DIM,X_DIM> SqrtSigma0;
extern Matrix<X_DIM> xGoal;
extern Matrix<X_DIM> xMin, xMax;
extern Matrix<U_DIM> uMin, uMax;

extern const int T;
extern const double INFTY;

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

extern double *inputVars, *vars;
extern std::vector<int> maskIndices;


// dynamics function
extern inline Matrix<X_DIM> dynfunc(const Matrix<X_DIM>& x, const Matrix<U_DIM>& u, const Matrix<U_DIM>& q);

// observation function
extern inline Matrix<Z_DIM> obsfunc(const Matrix<X_DIM>& x, const Matrix<R_DIM>& r);

// Jacobians: df(x,u,q)/dx, df(x,u,q)/dq
extern inline void linearizeDynamics(const Matrix<X_DIM>& x, const Matrix<U_DIM>& u, const Matrix<Q_DIM>& q, Matrix<X_DIM,X_DIM>& A, Matrix<X_DIM,Q_DIM>& M);

// Jacobians: dh(x,r)/dx, dh(x,r)/dr
extern inline void linearizeObservation(const Matrix<X_DIM>& x, const Matrix<R_DIM>& r, Matrix<Z_DIM,X_DIM>& H, Matrix<Z_DIM,R_DIM>& N);

// Switch between belief vector and matrices
extern inline void unVec(const Matrix<B_DIM>& b, Matrix<X_DIM>& x, Matrix<X_DIM,X_DIM>& SqrtSigma);

extern inline void vec(const Matrix<X_DIM>& x, const Matrix<X_DIM,X_DIM>& SqrtSigma, Matrix<B_DIM>& b);

// Belief dynamics
extern inline Matrix<B_DIM> beliefDynamics(const Matrix<B_DIM>& b, const Matrix<U_DIM>& u);

extern void pythonDisplayTrajectory(std::vector< Matrix<B_DIM> >& B, std::vector< Matrix<U_DIM> >& U);

#endif
