// Code generated from D* symbolic/automatic differentiation

#ifndef __symeval_H__
#define __symeval_H__

// For 15 timesteps for the point example

// Function definitions
void evalCost(double *result, double *vars);

void evalCostGrad(double *result, double *vars);

void evalCostGradDiagHess(double *result, double *vars);

#endif
