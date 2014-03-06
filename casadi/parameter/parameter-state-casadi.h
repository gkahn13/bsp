#ifndef __PARAMETER_STATE_CASADI_H__
#define __PARAMETER_STATE_CASADI_H__

#include <math.h>

inline double sq(double x){ return x*x;}
#define sq(x) ((x)*(x))

inline double sign(double x){ return x<0 ? -1 : x>0 ? 1 : x;}

int initDyn(int *n_in, int *n_out);
int getDynSparsity(int i, int *nrow, int *ncol, int **rowindouble, int **col);
void evaluateDyn(const double* x0,const double* x1,const double* x2,double* r0);
int evaluateDynWrap(const double** x, double** r);


int initObs(int *n_in, int *n_out);
int getObsSparsity(int i, int *nrow, int *ncol, int **rowindouble, int **col);
void evaluateObs(const double* x0,const double* x1,const double* x2,double* r0,double* r1);
int evaluateObsWrap(const double** x, double** r);


//#include "parameter-state-grad.c"
//#include "parameter-state-cost.c"

#endif
