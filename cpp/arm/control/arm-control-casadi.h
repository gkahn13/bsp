#ifndef __ARM_CONTROL_CASADI_H__
#define __ARM_CONTROL_CASADI_H__

#include <math.h>

inline double sq(double x){ return x*x;}
#define sq(x) ((x)*(x))

inline double sign(double x){ return x<0 ? -1 : x>0 ? 1 : x;}

int initCost(int *n_in, int *n_out);
int getCostSparsity(int i, int *nrow, int *ncol, int **rowindouble, int **col);
void evaluateCost(const double* x0,const double* x1,const double* x2,const double* x3,const double* x4,double* r0);
int evaluateCostWrap(const double** x, double** r);

int initCostGrad(int *n_in, int *n_out);
int getCostGradSparsity(int i, int *nrow, int *ncol, int **rowindouble, int **col);
void evaluateCostGrad(const double* x0,const double* x1,const double* x2,const double* x3,const double* x4,double* r0, double* r1);
int evaluateCostGradWrap(const double** x, double** r);

#endif
