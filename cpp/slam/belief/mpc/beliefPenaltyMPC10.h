/*
FORCES - Fast interior point code generation for multistage problems.
Copyright (C) 2011-14 Alexander Domahidi [domahidi@control.ee.ethz.ch],
Automatic Control Laboratory, ETH Zurich.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef __beliefPenaltyMPC_H__
#define __beliefPenaltyMPC_H__


/* DATA TYPE ------------------------------------------------------------*/
typedef double beliefPenaltyMPC_FLOAT;


/* SOLVER SETTINGS ------------------------------------------------------*/
/* print level */
#ifndef beliefPenaltyMPC_SET_PRINTLEVEL
#define beliefPenaltyMPC_SET_PRINTLEVEL    (1)
#endif

/* timing */
#ifndef beliefPenaltyMPC_SET_TIMING
#define beliefPenaltyMPC_SET_TIMING    (1)
#endif

/* Numeric Warnings */
/* #define PRINTNUMERICALWARNINGS */

/* maximum number of iterations  */
#define beliefPenaltyMPC_SET_MAXIT         (30)	

/* scaling factor of line search (affine direction) */
#define beliefPenaltyMPC_SET_LS_SCALE_AFF  (0.9)      

/* scaling factor of line search (combined direction) */
#define beliefPenaltyMPC_SET_LS_SCALE      (0.95)  

/* minimum required step size in each iteration */
#define beliefPenaltyMPC_SET_LS_MINSTEP    (1E-08)

/* maximum step size (combined direction) */
#define beliefPenaltyMPC_SET_LS_MAXSTEP    (0.995)

/* desired relative duality gap */
#define beliefPenaltyMPC_SET_ACC_RDGAP     (0.0001)

/* desired maximum residual on equality constraints */
#define beliefPenaltyMPC_SET_ACC_RESEQ     (1E-06)

/* desired maximum residual on inequality constraints */
#define beliefPenaltyMPC_SET_ACC_RESINEQ   (1E-06)

/* desired maximum violation of complementarity */
#define beliefPenaltyMPC_SET_ACC_KKTCOMPL  (1E-06)


/* RETURN CODES----------------------------------------------------------*/
/* solver has converged within desired accuracy */
#define beliefPenaltyMPC_OPTIMAL      (1)

/* maximum number of iterations has been reached */
#define beliefPenaltyMPC_MAXITREACHED (0)

/* no progress in line search possible */
#define beliefPenaltyMPC_NOPROGRESS   (-7)




/* PARAMETERS -----------------------------------------------------------*/
/* fill this with data before calling the solver! */
typedef struct beliefPenaltyMPC_params
{
    /* diagonal matrix of size [325 x 325] (only the diagonal is stored) */
    beliefPenaltyMPC_FLOAT H1[325];

    /* vector of size 325 */
    beliefPenaltyMPC_FLOAT f1[325];

    /* vector of size 325 */
    beliefPenaltyMPC_FLOAT lb1[325];

    /* vector of size 117 */
    beliefPenaltyMPC_FLOAT ub1[117];

    /* matrix of size [208 x 325] (column major format) */
    beliefPenaltyMPC_FLOAT C1[67600];

    /* vector of size 208 */
    beliefPenaltyMPC_FLOAT e1[208];

    /* diagonal matrix of size [325 x 325] (only the diagonal is stored) */
    beliefPenaltyMPC_FLOAT H2[325];

    /* vector of size 325 */
    beliefPenaltyMPC_FLOAT f2[325];

    /* vector of size 325 */
    beliefPenaltyMPC_FLOAT lb2[325];

    /* vector of size 117 */
    beliefPenaltyMPC_FLOAT ub2[117];

    /* matrix of size [104 x 325] (column major format) */
    beliefPenaltyMPC_FLOAT C2[33800];

    /* vector of size 104 */
    beliefPenaltyMPC_FLOAT e2[104];

    /* matrix of size [208 x 325] (column major format) */
    beliefPenaltyMPC_FLOAT D2[67600];

    /* diagonal matrix of size [325 x 325] (only the diagonal is stored) */
    beliefPenaltyMPC_FLOAT H3[325];

    /* vector of size 325 */
    beliefPenaltyMPC_FLOAT f3[325];

    /* vector of size 325 */
    beliefPenaltyMPC_FLOAT lb3[325];

    /* vector of size 117 */
    beliefPenaltyMPC_FLOAT ub3[117];

    /* matrix of size [104 x 325] (column major format) */
    beliefPenaltyMPC_FLOAT C3[33800];

    /* vector of size 104 */
    beliefPenaltyMPC_FLOAT e3[104];

    /* matrix of size [104 x 325] (column major format) */
    beliefPenaltyMPC_FLOAT D3[33800];

    /* diagonal matrix of size [325 x 325] (only the diagonal is stored) */
    beliefPenaltyMPC_FLOAT H4[325];

    /* vector of size 325 */
    beliefPenaltyMPC_FLOAT f4[325];

    /* vector of size 325 */
    beliefPenaltyMPC_FLOAT lb4[325];

    /* vector of size 117 */
    beliefPenaltyMPC_FLOAT ub4[117];

    /* matrix of size [104 x 325] (column major format) */
    beliefPenaltyMPC_FLOAT C4[33800];

    /* vector of size 104 */
    beliefPenaltyMPC_FLOAT e4[104];

    /* matrix of size [104 x 325] (column major format) */
    beliefPenaltyMPC_FLOAT D4[33800];

    /* diagonal matrix of size [325 x 325] (only the diagonal is stored) */
    beliefPenaltyMPC_FLOAT H5[325];

    /* vector of size 325 */
    beliefPenaltyMPC_FLOAT f5[325];

    /* vector of size 325 */
    beliefPenaltyMPC_FLOAT lb5[325];

    /* vector of size 117 */
    beliefPenaltyMPC_FLOAT ub5[117];

    /* matrix of size [104 x 325] (column major format) */
    beliefPenaltyMPC_FLOAT C5[33800];

    /* vector of size 104 */
    beliefPenaltyMPC_FLOAT e5[104];

    /* matrix of size [104 x 325] (column major format) */
    beliefPenaltyMPC_FLOAT D5[33800];

    /* diagonal matrix of size [325 x 325] (only the diagonal is stored) */
    beliefPenaltyMPC_FLOAT H6[325];

    /* vector of size 325 */
    beliefPenaltyMPC_FLOAT f6[325];

    /* vector of size 325 */
    beliefPenaltyMPC_FLOAT lb6[325];

    /* vector of size 117 */
    beliefPenaltyMPC_FLOAT ub6[117];

    /* matrix of size [104 x 325] (column major format) */
    beliefPenaltyMPC_FLOAT C6[33800];

    /* vector of size 104 */
    beliefPenaltyMPC_FLOAT e6[104];

    /* matrix of size [104 x 325] (column major format) */
    beliefPenaltyMPC_FLOAT D6[33800];

    /* diagonal matrix of size [325 x 325] (only the diagonal is stored) */
    beliefPenaltyMPC_FLOAT H7[325];

    /* vector of size 325 */
    beliefPenaltyMPC_FLOAT f7[325];

    /* vector of size 325 */
    beliefPenaltyMPC_FLOAT lb7[325];

    /* vector of size 117 */
    beliefPenaltyMPC_FLOAT ub7[117];

    /* matrix of size [104 x 325] (column major format) */
    beliefPenaltyMPC_FLOAT C7[33800];

    /* vector of size 104 */
    beliefPenaltyMPC_FLOAT e7[104];

    /* matrix of size [104 x 325] (column major format) */
    beliefPenaltyMPC_FLOAT D7[33800];

    /* diagonal matrix of size [325 x 325] (only the diagonal is stored) */
    beliefPenaltyMPC_FLOAT H8[325];

    /* vector of size 325 */
    beliefPenaltyMPC_FLOAT f8[325];

    /* vector of size 325 */
    beliefPenaltyMPC_FLOAT lb8[325];

    /* vector of size 117 */
    beliefPenaltyMPC_FLOAT ub8[117];

    /* matrix of size [104 x 325] (column major format) */
    beliefPenaltyMPC_FLOAT C8[33800];

    /* vector of size 104 */
    beliefPenaltyMPC_FLOAT e8[104];

    /* matrix of size [104 x 325] (column major format) */
    beliefPenaltyMPC_FLOAT D8[33800];

    /* diagonal matrix of size [325 x 325] (only the diagonal is stored) */
    beliefPenaltyMPC_FLOAT H9[325];

    /* vector of size 325 */
    beliefPenaltyMPC_FLOAT f9[325];

    /* vector of size 325 */
    beliefPenaltyMPC_FLOAT lb9[325];

    /* vector of size 117 */
    beliefPenaltyMPC_FLOAT ub9[117];

    /* matrix of size [104 x 325] (column major format) */
    beliefPenaltyMPC_FLOAT C9[33800];

    /* vector of size 104 */
    beliefPenaltyMPC_FLOAT e9[104];

    /* matrix of size [104 x 325] (column major format) */
    beliefPenaltyMPC_FLOAT D9[33800];

    /* diagonal matrix of size [104 x 104] (only the diagonal is stored) */
    beliefPenaltyMPC_FLOAT H10[104];

    /* vector of size 104 */
    beliefPenaltyMPC_FLOAT lb10[104];

    /* vector of size 104 */
    beliefPenaltyMPC_FLOAT ub10[104];

    /* matrix of size [104 x 104] (column major format) */
    beliefPenaltyMPC_FLOAT D10[10816];

} beliefPenaltyMPC_params;


/* OUTPUTS --------------------------------------------------------------*/
/* the desired variables are put here by the solver */
typedef struct beliefPenaltyMPC_output
{
    /* vector of size 117 */
    beliefPenaltyMPC_FLOAT z1[117];

    /* vector of size 117 */
    beliefPenaltyMPC_FLOAT z2[117];

    /* vector of size 117 */
    beliefPenaltyMPC_FLOAT z3[117];

    /* vector of size 117 */
    beliefPenaltyMPC_FLOAT z4[117];

    /* vector of size 117 */
    beliefPenaltyMPC_FLOAT z5[117];

    /* vector of size 117 */
    beliefPenaltyMPC_FLOAT z6[117];

    /* vector of size 117 */
    beliefPenaltyMPC_FLOAT z7[117];

    /* vector of size 117 */
    beliefPenaltyMPC_FLOAT z8[117];

    /* vector of size 117 */
    beliefPenaltyMPC_FLOAT z9[117];

    /* vector of size 104 */
    beliefPenaltyMPC_FLOAT z10[104];

} beliefPenaltyMPC_output;


/* SOLVER INFO ----------------------------------------------------------*/
/* diagnostic data from last interior point step */
typedef struct beliefPenaltyMPC_info
{
    /* iteration number */
    int it;
	
    /* inf-norm of equality constraint residuals */
    beliefPenaltyMPC_FLOAT res_eq;
	
    /* inf-norm of inequality constraint residuals */
    beliefPenaltyMPC_FLOAT res_ineq;

    /* primal objective */
    beliefPenaltyMPC_FLOAT pobj;	
	
    /* dual objective */
    beliefPenaltyMPC_FLOAT dobj;	

    /* duality gap := pobj - dobj */
    beliefPenaltyMPC_FLOAT dgap;		
	
    /* relative duality gap := |dgap / pobj | */
    beliefPenaltyMPC_FLOAT rdgap;		

    /* duality measure */
    beliefPenaltyMPC_FLOAT mu;

	/* duality measure (after affine step) */
    beliefPenaltyMPC_FLOAT mu_aff;
	
    /* centering parameter */
    beliefPenaltyMPC_FLOAT sigma;
	
    /* number of backtracking line search steps (affine direction) */
    int lsit_aff;
    
    /* number of backtracking line search steps (combined direction) */
    int lsit_cc;
    
    /* step size (affine direction) */
    beliefPenaltyMPC_FLOAT step_aff;
    
    /* step size (combined direction) */
    beliefPenaltyMPC_FLOAT step_cc;    

	/* solvertime */
	beliefPenaltyMPC_FLOAT solvetime;   

} beliefPenaltyMPC_info;


/* SOLVER FUNCTION DEFINITION -------------------------------------------*/
/* examine exitflag before using the result! */
int beliefPenaltyMPC_solve(beliefPenaltyMPC_params* params, beliefPenaltyMPC_output* output, beliefPenaltyMPC_info* info);


#endif