#include <math.h>

double sq(double x){ return x*x;}

double sign(double x){ return x<0 ? -1 : x>0 ? 1 : x;}

int initCost(int *n_in, int *n_out);
int getCostSparsity(int i, int *nrow, int *ncol, int **rowindouble, int **col);
double evaluateCost(const double* x0,const double* x1,const double* x2,const double* x3,const double* x4,double* r0);
int evaluateCostWrap(const double** x, double** r);

int initGrad(int *n_in, int *n_out);
int getGradSparsity(int i, int *nrow, int *ncol, int **rowindouble, int **col);
void evaluateGrad(const double* x0,const double* x1,const double* x2,const double* x3,const double* x4,double* r0);
int evaluateGradWrap(const double** x, double** r);