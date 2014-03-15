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

#ifndef __controlPenaltyMPC_H__
#define __controlPenaltyMPC_H__


/* DATA TYPE ------------------------------------------------------------*/
typedef double controlPenaltyMPC_FLOAT;


/* SOLVER SETTINGS ------------------------------------------------------*/
/* print level */
#ifndef controlPenaltyMPC_SET_PRINTLEVEL
#define controlPenaltyMPC_SET_PRINTLEVEL    (0)
#endif

/* timing */
#ifndef controlPenaltyMPC_SET_TIMING
#define controlPenaltyMPC_SET_TIMING    (0)
#endif

/* Numeric Warnings */
/* #define PRINTNUMERICALWARNINGS */

/* maximum number of iterations  */
#define controlPenaltyMPC_SET_MAXIT         (40)	

/* scaling factor of line search (affine direction) */
#define controlPenaltyMPC_SET_LS_SCALE_AFF  (0.9)      

/* scaling factor of line search (combined direction) */
#define controlPenaltyMPC_SET_LS_SCALE      (0.95)  

/* minimum required step size in each iteration */
#define controlPenaltyMPC_SET_LS_MINSTEP    (1E-08)

/* maximum step size (combined direction) */
#define controlPenaltyMPC_SET_LS_MAXSTEP    (0.995)

/* desired relative duality gap */
#define controlPenaltyMPC_SET_ACC_RDGAP     (0.0001)

/* desired maximum residual on equality constraints */
#define controlPenaltyMPC_SET_ACC_RESEQ     (1E-06)

/* desired maximum residual on inequality constraints */
#define controlPenaltyMPC_SET_ACC_RESINEQ   (1E-06)

/* desired maximum violation of complementarity */
#define controlPenaltyMPC_SET_ACC_KKTCOMPL  (1E-06)


/* RETURN CODES----------------------------------------------------------*/
/* solver has converged within desired accuracy */
#define controlPenaltyMPC_OPTIMAL      (1)

/* maximum number of iterations has been reached */
#define controlPenaltyMPC_MAXITREACHED (0)

/* no progress in line search possible */
#define controlPenaltyMPC_NOPROGRESS   (-7)




/* PARAMETERS -----------------------------------------------------------*/
/* fill this with data before calling the solver! */
typedef struct controlPenaltyMPC_params
{
    /* diagonal matrix of size [12 x 12] (only the diagonal is stored) */
    controlPenaltyMPC_FLOAT Q1[12];

    /* vector of size 12 */
    controlPenaltyMPC_FLOAT f1[12];

    /* vector of size 12 */
    controlPenaltyMPC_FLOAT lb1[12];

    /* vector of size 12 */
    controlPenaltyMPC_FLOAT ub1[12];

    /* vector of size 6 */
    controlPenaltyMPC_FLOAT e1[6];

    /* diagonal matrix of size [12 x 12] (only the diagonal is stored) */
    controlPenaltyMPC_FLOAT Q2[12];

    /* vector of size 12 */
    controlPenaltyMPC_FLOAT f2[12];

    /* vector of size 12 */
    controlPenaltyMPC_FLOAT lb2[12];

    /* vector of size 12 */
    controlPenaltyMPC_FLOAT ub2[12];

    /* diagonal matrix of size [12 x 12] (only the diagonal is stored) */
    controlPenaltyMPC_FLOAT Q3[12];

    /* vector of size 12 */
    controlPenaltyMPC_FLOAT f3[12];

    /* vector of size 12 */
    controlPenaltyMPC_FLOAT lb3[12];

    /* vector of size 12 */
    controlPenaltyMPC_FLOAT ub3[12];

    /* diagonal matrix of size [12 x 12] (only the diagonal is stored) */
    controlPenaltyMPC_FLOAT Q4[12];

    /* vector of size 12 */
    controlPenaltyMPC_FLOAT f4[12];

    /* vector of size 12 */
    controlPenaltyMPC_FLOAT lb4[12];

    /* vector of size 12 */
    controlPenaltyMPC_FLOAT ub4[12];

    /* diagonal matrix of size [12 x 12] (only the diagonal is stored) */
    controlPenaltyMPC_FLOAT Q5[12];

    /* vector of size 12 */
    controlPenaltyMPC_FLOAT f5[12];

    /* vector of size 12 */
    controlPenaltyMPC_FLOAT lb5[12];

    /* vector of size 12 */
    controlPenaltyMPC_FLOAT ub5[12];

    /* diagonal matrix of size [12 x 12] (only the diagonal is stored) */
    controlPenaltyMPC_FLOAT Q6[12];

    /* vector of size 12 */
    controlPenaltyMPC_FLOAT f6[12];

    /* vector of size 12 */
    controlPenaltyMPC_FLOAT lb6[12];

    /* vector of size 12 */
    controlPenaltyMPC_FLOAT ub6[12];

    /* diagonal matrix of size [12 x 12] (only the diagonal is stored) */
    controlPenaltyMPC_FLOAT Q7[12];

    /* vector of size 12 */
    controlPenaltyMPC_FLOAT f7[12];

    /* vector of size 12 */
    controlPenaltyMPC_FLOAT lb7[12];

    /* vector of size 12 */
    controlPenaltyMPC_FLOAT ub7[12];

    /* diagonal matrix of size [12 x 12] (only the diagonal is stored) */
    controlPenaltyMPC_FLOAT Q8[12];

    /* vector of size 12 */
    controlPenaltyMPC_FLOAT f8[12];

    /* vector of size 12 */
    controlPenaltyMPC_FLOAT lb8[12];

    /* vector of size 12 */
    controlPenaltyMPC_FLOAT ub8[12];

    /* diagonal matrix of size [12 x 12] (only the diagonal is stored) */
    controlPenaltyMPC_FLOAT Q9[12];

    /* vector of size 12 */
    controlPenaltyMPC_FLOAT f9[12];

    /* vector of size 12 */
    controlPenaltyMPC_FLOAT lb9[12];

    /* vector of size 12 */
    controlPenaltyMPC_FLOAT ub9[12];

    /* diagonal matrix of size [12 x 12] (only the diagonal is stored) */
    controlPenaltyMPC_FLOAT Q10[12];

    /* vector of size 12 */
    controlPenaltyMPC_FLOAT f10[12];

    /* vector of size 12 */
    controlPenaltyMPC_FLOAT lb10[12];

    /* vector of size 12 */
    controlPenaltyMPC_FLOAT ub10[12];

    /* diagonal matrix of size [12 x 12] (only the diagonal is stored) */
    controlPenaltyMPC_FLOAT Q11[12];

    /* vector of size 12 */
    controlPenaltyMPC_FLOAT f11[12];

    /* vector of size 12 */
    controlPenaltyMPC_FLOAT lb11[12];

    /* vector of size 12 */
    controlPenaltyMPC_FLOAT ub11[12];

    /* diagonal matrix of size [12 x 12] (only the diagonal is stored) */
    controlPenaltyMPC_FLOAT Q12[12];

    /* vector of size 12 */
    controlPenaltyMPC_FLOAT f12[12];

    /* vector of size 12 */
    controlPenaltyMPC_FLOAT lb12[12];

    /* vector of size 12 */
    controlPenaltyMPC_FLOAT ub12[12];

    /* diagonal matrix of size [12 x 12] (only the diagonal is stored) */
    controlPenaltyMPC_FLOAT Q13[12];

    /* vector of size 12 */
    controlPenaltyMPC_FLOAT f13[12];

    /* vector of size 12 */
    controlPenaltyMPC_FLOAT lb13[12];

    /* vector of size 12 */
    controlPenaltyMPC_FLOAT ub13[12];

    /* diagonal matrix of size [12 x 12] (only the diagonal is stored) */
    controlPenaltyMPC_FLOAT Q14[12];

    /* vector of size 12 */
    controlPenaltyMPC_FLOAT f14[12];

    /* vector of size 12 */
    controlPenaltyMPC_FLOAT lb14[12];

    /* vector of size 12 */
    controlPenaltyMPC_FLOAT ub14[12];

    /* diagonal matrix of size [6 x 6] (only the diagonal is stored) */
    controlPenaltyMPC_FLOAT Q15[6];

    /* vector of size 6 */
    controlPenaltyMPC_FLOAT f15[6];

    /* vector of size 6 */
    controlPenaltyMPC_FLOAT lb15[6];

    /* vector of size 6 */
    controlPenaltyMPC_FLOAT ub15[6];

} controlPenaltyMPC_params;


/* OUTPUTS --------------------------------------------------------------*/
/* the desired variables are put here by the solver */
typedef struct controlPenaltyMPC_output
{
    /* vector of size 12 */
    controlPenaltyMPC_FLOAT z1[12];

    /* vector of size 12 */
    controlPenaltyMPC_FLOAT z2[12];

    /* vector of size 12 */
    controlPenaltyMPC_FLOAT z3[12];

    /* vector of size 12 */
    controlPenaltyMPC_FLOAT z4[12];

    /* vector of size 12 */
    controlPenaltyMPC_FLOAT z5[12];

    /* vector of size 12 */
    controlPenaltyMPC_FLOAT z6[12];

    /* vector of size 12 */
    controlPenaltyMPC_FLOAT z7[12];

    /* vector of size 12 */
    controlPenaltyMPC_FLOAT z8[12];

    /* vector of size 12 */
    controlPenaltyMPC_FLOAT z9[12];

    /* vector of size 12 */
    controlPenaltyMPC_FLOAT z10[12];

    /* vector of size 12 */
    controlPenaltyMPC_FLOAT z11[12];

    /* vector of size 12 */
    controlPenaltyMPC_FLOAT z12[12];

    /* vector of size 12 */
    controlPenaltyMPC_FLOAT z13[12];

    /* vector of size 12 */
    controlPenaltyMPC_FLOAT z14[12];

    /* vector of size 6 */
    controlPenaltyMPC_FLOAT z15[6];

} controlPenaltyMPC_output;


/* SOLVER INFO ----------------------------------------------------------*/
/* diagnostic data from last interior point step */
typedef struct controlPenaltyMPC_info
{
    /* iteration number */
    int it;
	
    /* inf-norm of equality constraint residuals */
    controlPenaltyMPC_FLOAT res_eq;
	
    /* inf-norm of inequality constraint residuals */
    controlPenaltyMPC_FLOAT res_ineq;

    /* primal objective */
    controlPenaltyMPC_FLOAT pobj;	
	
    /* dual objective */
    controlPenaltyMPC_FLOAT dobj;	

    /* duality gap := pobj - dobj */
    controlPenaltyMPC_FLOAT dgap;		
	
    /* relative duality gap := |dgap / pobj | */
    controlPenaltyMPC_FLOAT rdgap;		

    /* duality measure */
    controlPenaltyMPC_FLOAT mu;

	/* duality measure (after affine step) */
    controlPenaltyMPC_FLOAT mu_aff;
	
    /* centering parameter */
    controlPenaltyMPC_FLOAT sigma;
	
    /* number of backtracking line search steps (affine direction) */
    int lsit_aff;
    
    /* number of backtracking line search steps (combined direction) */
    int lsit_cc;
    
    /* step size (affine direction) */
    controlPenaltyMPC_FLOAT step_aff;
    
    /* step size (combined direction) */
    controlPenaltyMPC_FLOAT step_cc;    

	/* solvertime */
	controlPenaltyMPC_FLOAT solvetime;   

} controlPenaltyMPC_info;


/* SOLVER FUNCTION DEFINITION -------------------------------------------*/
/* examine exitflag before using the result! */
int controlPenaltyMPC_solve(controlPenaltyMPC_params* params, controlPenaltyMPC_output* output, controlPenaltyMPC_info* info);


#endif