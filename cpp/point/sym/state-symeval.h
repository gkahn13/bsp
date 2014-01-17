// Code generated from D* symbolic/automatic differentiation
// common across all time steps
#ifndef __state_symeval_H__
#define __state_symeval_H__

// Function definitions
void evalCost(double *result, double *vars);

void evalCostGrad(double *result, double *vars);

void evalCostGradDiagHess(double *result, double *vars);

#endif