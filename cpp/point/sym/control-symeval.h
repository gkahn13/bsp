// Code generated from D* symbolic/automatic differentiation
// common across all time steps
#ifndef __control_symeval_H__
#define __control_symeval_H__

// Function definitions
void evalCost(double *result, double *vars);

void evalCostGradDiagHess(double *result, double *vars);

char* getMask();

#endif
