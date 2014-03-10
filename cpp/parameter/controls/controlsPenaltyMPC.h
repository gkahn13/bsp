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

#ifndef __controlsPenaltyMPC_H__
#define __controlsPenaltyMPC_H__


/* DATA TYPE ------------------------------------------------------------*/
typedef double controlsPenaltyMPC_FLOAT;


/* SOLVER SETTINGS ------------------------------------------------------*/
/* print level */
#ifndef controlsPenaltyMPC_SET_PRINTLEVEL
#define controlsPenaltyMPC_SET_PRINTLEVEL    (0)
#endif

/* timing */
#ifndef controlsPenaltyMPC_SET_TIMING
#define controlsPenaltyMPC_SET_TIMING    (0)
#endif

/* Numeric Warnings */
/* #define PRINTNUMERICALWARNINGS */

/* maximum number of iterations  */
#define controlsPenaltyMPC_SET_MAXIT         (60)	

/* scaling factor of line search (affine direction) */
#define controlsPenaltyMPC_SET_LS_SCALE_AFF  (0.9)      

/* scaling factor of line search (combined direction) */
#define controlsPenaltyMPC_SET_LS_SCALE      (0.95)  

/* minimum required step size in each iteration */
#define controlsPenaltyMPC_SET_LS_MINSTEP    (1E-08)

/* maximum step size (combined direction) */
#define controlsPenaltyMPC_SET_LS_MAXSTEP    (0.995)

/* desired relative duality gap */
#define controlsPenaltyMPC_SET_ACC_RDGAP     (0.0001)

/* desired maximum residual on equality constraints */
#define controlsPenaltyMPC_SET_ACC_RESEQ     (1E-06)

/* desired maximum residual on inequality constraints */
#define controlsPenaltyMPC_SET_ACC_RESINEQ   (1E-06)

/* desired maximum violation of complementarity */
#define controlsPenaltyMPC_SET_ACC_KKTCOMPL  (1E-06)


/* RETURN CODES----------------------------------------------------------*/
/* solver has converged within desired accuracy */
#define controlsPenaltyMPC_OPTIMAL      (1)

/* maximum number of iterations has been reached */
#define controlsPenaltyMPC_MAXITREACHED (0)

/* no progress in line search possible */
#define controlsPenaltyMPC_NOPROGRESS   (-7)




/* PARAMETERS -----------------------------------------------------------*/
/* fill this with data before calling the solver! */
typedef struct controlsPenaltyMPC_params
{
    /* diagonal matrix of size [2 x 2] (only the diagonal is stored) */
    controlsPenaltyMPC_FLOAT H1[2];

    /* vector of size 2 */
    controlsPenaltyMPC_FLOAT f1[2];

    /* vector of size 2 */
    controlsPenaltyMPC_FLOAT lb1[2];

    /* vector of size 2 */
    controlsPenaltyMPC_FLOAT ub1[2];

    /* diagonal matrix of size [2 x 2] (only the diagonal is stored) */
    controlsPenaltyMPC_FLOAT H2[2];

    /* vector of size 2 */
    controlsPenaltyMPC_FLOAT f2[2];

    /* vector of size 2 */
    controlsPenaltyMPC_FLOAT lb2[2];

    /* vector of size 2 */
    controlsPenaltyMPC_FLOAT ub2[2];

    /* diagonal matrix of size [2 x 2] (only the diagonal is stored) */
    controlsPenaltyMPC_FLOAT H3[2];

    /* vector of size 2 */
    controlsPenaltyMPC_FLOAT f3[2];

    /* vector of size 2 */
    controlsPenaltyMPC_FLOAT lb3[2];

    /* vector of size 2 */
    controlsPenaltyMPC_FLOAT ub3[2];

    /* diagonal matrix of size [2 x 2] (only the diagonal is stored) */
    controlsPenaltyMPC_FLOAT H4[2];

    /* vector of size 2 */
    controlsPenaltyMPC_FLOAT f4[2];

    /* vector of size 2 */
    controlsPenaltyMPC_FLOAT lb4[2];

    /* vector of size 2 */
    controlsPenaltyMPC_FLOAT ub4[2];

    /* diagonal matrix of size [2 x 2] (only the diagonal is stored) */
    controlsPenaltyMPC_FLOAT H5[2];

    /* vector of size 2 */
    controlsPenaltyMPC_FLOAT f5[2];

    /* vector of size 2 */
    controlsPenaltyMPC_FLOAT lb5[2];

    /* vector of size 2 */
    controlsPenaltyMPC_FLOAT ub5[2];

    /* diagonal matrix of size [2 x 2] (only the diagonal is stored) */
    controlsPenaltyMPC_FLOAT H6[2];

    /* vector of size 2 */
    controlsPenaltyMPC_FLOAT f6[2];

    /* vector of size 2 */
    controlsPenaltyMPC_FLOAT lb6[2];

    /* vector of size 2 */
    controlsPenaltyMPC_FLOAT ub6[2];

    /* diagonal matrix of size [2 x 2] (only the diagonal is stored) */
    controlsPenaltyMPC_FLOAT H7[2];

    /* vector of size 2 */
    controlsPenaltyMPC_FLOAT f7[2];

    /* vector of size 2 */
    controlsPenaltyMPC_FLOAT lb7[2];

    /* vector of size 2 */
    controlsPenaltyMPC_FLOAT ub7[2];

    /* diagonal matrix of size [2 x 2] (only the diagonal is stored) */
    controlsPenaltyMPC_FLOAT H8[2];

    /* vector of size 2 */
    controlsPenaltyMPC_FLOAT f8[2];

    /* vector of size 2 */
    controlsPenaltyMPC_FLOAT lb8[2];

    /* vector of size 2 */
    controlsPenaltyMPC_FLOAT ub8[2];

    /* diagonal matrix of size [2 x 2] (only the diagonal is stored) */
    controlsPenaltyMPC_FLOAT H9[2];

    /* vector of size 2 */
    controlsPenaltyMPC_FLOAT f9[2];

    /* vector of size 2 */
    controlsPenaltyMPC_FLOAT lb9[2];

    /* vector of size 2 */
    controlsPenaltyMPC_FLOAT ub9[2];

    /* diagonal matrix of size [2 x 2] (only the diagonal is stored) */
    controlsPenaltyMPC_FLOAT H10[2];

    /* vector of size 2 */
    controlsPenaltyMPC_FLOAT f10[2];

    /* vector of size 2 */
    controlsPenaltyMPC_FLOAT lb10[2];

    /* vector of size 2 */
    controlsPenaltyMPC_FLOAT ub10[2];

    /* diagonal matrix of size [2 x 2] (only the diagonal is stored) */
    controlsPenaltyMPC_FLOAT H11[2];

    /* vector of size 2 */
    controlsPenaltyMPC_FLOAT f11[2];

    /* vector of size 2 */
    controlsPenaltyMPC_FLOAT lb11[2];

    /* vector of size 2 */
    controlsPenaltyMPC_FLOAT ub11[2];

    /* diagonal matrix of size [2 x 2] (only the diagonal is stored) */
    controlsPenaltyMPC_FLOAT H12[2];

    /* vector of size 2 */
    controlsPenaltyMPC_FLOAT f12[2];

    /* vector of size 2 */
    controlsPenaltyMPC_FLOAT lb12[2];

    /* vector of size 2 */
    controlsPenaltyMPC_FLOAT ub12[2];

    /* diagonal matrix of size [2 x 2] (only the diagonal is stored) */
    controlsPenaltyMPC_FLOAT H13[2];

    /* vector of size 2 */
    controlsPenaltyMPC_FLOAT f13[2];

    /* vector of size 2 */
    controlsPenaltyMPC_FLOAT lb13[2];

    /* vector of size 2 */
    controlsPenaltyMPC_FLOAT ub13[2];

    /* diagonal matrix of size [2 x 2] (only the diagonal is stored) */
    controlsPenaltyMPC_FLOAT H14[2];

    /* vector of size 2 */
    controlsPenaltyMPC_FLOAT f14[2];

    /* vector of size 2 */
    controlsPenaltyMPC_FLOAT lb14[2];

    /* vector of size 2 */
    controlsPenaltyMPC_FLOAT ub14[2];

} controlsPenaltyMPC_params;


/* OUTPUTS --------------------------------------------------------------*/
/* the desired variables are put here by the solver */
typedef struct controlsPenaltyMPC_output
{
    /* vector of size 2 */
    controlsPenaltyMPC_FLOAT z1[2];

    /* vector of size 2 */
    controlsPenaltyMPC_FLOAT z2[2];

    /* vector of size 2 */
    controlsPenaltyMPC_FLOAT z3[2];

    /* vector of size 2 */
    controlsPenaltyMPC_FLOAT z4[2];

    /* vector of size 2 */
    controlsPenaltyMPC_FLOAT z5[2];

    /* vector of size 2 */
    controlsPenaltyMPC_FLOAT z6[2];

    /* vector of size 2 */
    controlsPenaltyMPC_FLOAT z7[2];

    /* vector of size 2 */
    controlsPenaltyMPC_FLOAT z8[2];

    /* vector of size 2 */
    controlsPenaltyMPC_FLOAT z9[2];

    /* vector of size 2 */
    controlsPenaltyMPC_FLOAT z10[2];

    /* vector of size 2 */
    controlsPenaltyMPC_FLOAT z11[2];

    /* vector of size 2 */
    controlsPenaltyMPC_FLOAT z12[2];

    /* vector of size 2 */
    controlsPenaltyMPC_FLOAT z13[2];

    /* vector of size 2 */
    controlsPenaltyMPC_FLOAT z14[2];

} controlsPenaltyMPC_output;


/* SOLVER INFO ----------------------------------------------------------*/
/* diagnostic data from last interior point step */
typedef struct controlsPenaltyMPC_info
{
    /* iteration number */
    int it;
	
    /* inf-norm of equality constraint residuals */
    controlsPenaltyMPC_FLOAT res_eq;
	
    /* inf-norm of inequality constraint residuals */
    controlsPenaltyMPC_FLOAT res_ineq;

    /* primal objective */
    controlsPenaltyMPC_FLOAT pobj;	
	
    /* dual objective */
    controlsPenaltyMPC_FLOAT dobj;	

    /* duality gap := pobj - dobj */
    controlsPenaltyMPC_FLOAT dgap;		
	
    /* relative duality gap := |dgap / pobj | */
    controlsPenaltyMPC_FLOAT rdgap;		

    /* duality measure */
    controlsPenaltyMPC_FLOAT mu;

	/* duality measure (after affine step) */
    controlsPenaltyMPC_FLOAT mu_aff;
	
    /* centering parameter */
    controlsPenaltyMPC_FLOAT sigma;
	
    /* number of backtracking line search steps (affine direction) */
    int lsit_aff;
    
    /* number of backtracking line search steps (combined direction) */
    int lsit_cc;
    
    /* step size (affine direction) */
    controlsPenaltyMPC_FLOAT step_aff;
    
    /* step size (combined direction) */
    controlsPenaltyMPC_FLOAT step_cc;    

	/* solvertime */
	controlsPenaltyMPC_FLOAT solvetime;   

} controlsPenaltyMPC_info;


/* SOLVER FUNCTION DEFINITION -------------------------------------------*/
/* examine exitflag before using the result! */
int controlsPenaltyMPC_solve(controlsPenaltyMPC_params* params, controlsPenaltyMPC_output* output, controlsPenaltyMPC_info* info);


#endif