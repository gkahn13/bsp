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

#ifndef __planarMPC_H__
#define __planarMPC_H__


/* DATA TYPE ------------------------------------------------------------*/
typedef double planarMPC_FLOAT;


/* SOLVER SETTINGS ------------------------------------------------------*/
/* print level */
#ifndef planarMPC_SET_PRINTLEVEL
#define planarMPC_SET_PRINTLEVEL    (1)
#endif

/* timing */
#ifndef planarMPC_SET_TIMING
#define planarMPC_SET_TIMING    (0)
#endif

/* Numeric Warnings */
/* #define PRINTNUMERICALWARNINGS */

/* maximum number of iterations  */
#define planarMPC_SET_MAXIT         (100)	

/* scaling factor of line search (affine direction) */
#define planarMPC_SET_LS_SCALE_AFF  (0.9)      

/* scaling factor of line search (combined direction) */
#define planarMPC_SET_LS_SCALE      (0.95)  

/* minimum required step size in each iteration */
#define planarMPC_SET_LS_MINSTEP    (1E-08)

/* maximum step size (combined direction) */
#define planarMPC_SET_LS_MAXSTEP    (0.995)

/* desired relative duality gap */
#define planarMPC_SET_ACC_RDGAP     (0.0001)

/* desired maximum residual on equality constraints */
#define planarMPC_SET_ACC_RESEQ     (1E-06)

/* desired maximum residual on inequality constraints */
#define planarMPC_SET_ACC_RESINEQ   (1E-06)

/* desired maximum violation of complementarity */
#define planarMPC_SET_ACC_KKTCOMPL  (1E-06)


/* RETURN CODES----------------------------------------------------------*/
/* solver has converged within desired accuracy */
#define planarMPC_OPTIMAL      (1)

/* maximum number of iterations has been reached */
#define planarMPC_MAXITREACHED (0)

/* no progress in line search possible */
#define planarMPC_NOPROGRESS   (-7)




/* PARAMETERS -----------------------------------------------------------*/
/* fill this with data before calling the solver! */
typedef struct planarMPC_params
{
    /* diagonal matrix of size [8 x 8] (only the diagonal is stored) */
    planarMPC_FLOAT H1[8];

    /* vector of size 8 */
    planarMPC_FLOAT f1[8];

    /* vector of size 8 */
    planarMPC_FLOAT lb1[8];

    /* vector of size 8 */
    planarMPC_FLOAT ub1[8];

    /* vector of size 4 */
    planarMPC_FLOAT c1[4];

    /* diagonal matrix of size [8 x 8] (only the diagonal is stored) */
    planarMPC_FLOAT H2[8];

    /* vector of size 8 */
    planarMPC_FLOAT f2[8];

    /* vector of size 8 */
    planarMPC_FLOAT lb2[8];

    /* vector of size 8 */
    planarMPC_FLOAT ub2[8];

    /* diagonal matrix of size [8 x 8] (only the diagonal is stored) */
    planarMPC_FLOAT H3[8];

    /* vector of size 8 */
    planarMPC_FLOAT f3[8];

    /* vector of size 8 */
    planarMPC_FLOAT lb3[8];

    /* vector of size 8 */
    planarMPC_FLOAT ub3[8];

    /* diagonal matrix of size [8 x 8] (only the diagonal is stored) */
    planarMPC_FLOAT H4[8];

    /* vector of size 8 */
    planarMPC_FLOAT f4[8];

    /* vector of size 8 */
    planarMPC_FLOAT lb4[8];

    /* vector of size 8 */
    planarMPC_FLOAT ub4[8];

    /* diagonal matrix of size [8 x 8] (only the diagonal is stored) */
    planarMPC_FLOAT H5[8];

    /* vector of size 8 */
    planarMPC_FLOAT f5[8];

    /* vector of size 8 */
    planarMPC_FLOAT lb5[8];

    /* vector of size 8 */
    planarMPC_FLOAT ub5[8];

    /* diagonal matrix of size [8 x 8] (only the diagonal is stored) */
    planarMPC_FLOAT H6[8];

    /* vector of size 8 */
    planarMPC_FLOAT f6[8];

    /* vector of size 8 */
    planarMPC_FLOAT lb6[8];

    /* vector of size 8 */
    planarMPC_FLOAT ub6[8];

    /* diagonal matrix of size [8 x 8] (only the diagonal is stored) */
    planarMPC_FLOAT H7[8];

    /* vector of size 8 */
    planarMPC_FLOAT f7[8];

    /* vector of size 8 */
    planarMPC_FLOAT lb7[8];

    /* vector of size 8 */
    planarMPC_FLOAT ub7[8];

    /* diagonal matrix of size [8 x 8] (only the diagonal is stored) */
    planarMPC_FLOAT H8[8];

    /* vector of size 8 */
    planarMPC_FLOAT f8[8];

    /* vector of size 8 */
    planarMPC_FLOAT lb8[8];

    /* vector of size 8 */
    planarMPC_FLOAT ub8[8];

    /* diagonal matrix of size [8 x 8] (only the diagonal is stored) */
    planarMPC_FLOAT H9[8];

    /* vector of size 8 */
    planarMPC_FLOAT f9[8];

    /* vector of size 8 */
    planarMPC_FLOAT lb9[8];

    /* vector of size 8 */
    planarMPC_FLOAT ub9[8];

    /* diagonal matrix of size [4 x 4] (only the diagonal is stored) */
    planarMPC_FLOAT H10[4];

    /* vector of size 4 */
    planarMPC_FLOAT f10[4];

    /* vector of size 4 */
    planarMPC_FLOAT lb10[4];

    /* vector of size 4 */
    planarMPC_FLOAT ub10[4];

    /* matrix of size [4 x 4] (column major format) */
    planarMPC_FLOAT A10[16];

    /* vector of size 4 */
    planarMPC_FLOAT b10[4];

} planarMPC_params;


/* OUTPUTS --------------------------------------------------------------*/
/* the desired variables are put here by the solver */
typedef struct planarMPC_output
{
    /* vector of size 8 */
    planarMPC_FLOAT z1[8];

    /* vector of size 8 */
    planarMPC_FLOAT z2[8];

    /* vector of size 8 */
    planarMPC_FLOAT z3[8];

    /* vector of size 8 */
    planarMPC_FLOAT z4[8];

    /* vector of size 8 */
    planarMPC_FLOAT z5[8];

    /* vector of size 8 */
    planarMPC_FLOAT z6[8];

    /* vector of size 8 */
    planarMPC_FLOAT z7[8];

    /* vector of size 8 */
    planarMPC_FLOAT z8[8];

    /* vector of size 8 */
    planarMPC_FLOAT z9[8];

    /* vector of size 4 */
    planarMPC_FLOAT z10[4];

} planarMPC_output;


/* SOLVER INFO ----------------------------------------------------------*/
/* diagnostic data from last interior point step */
typedef struct planarMPC_info
{
    /* iteration number */
    int it;
	
    /* inf-norm of equality constraint residuals */
    planarMPC_FLOAT res_eq;
	
    /* inf-norm of inequality constraint residuals */
    planarMPC_FLOAT res_ineq;

    /* primal objective */
    planarMPC_FLOAT pobj;	
	
    /* dual objective */
    planarMPC_FLOAT dobj;	

    /* duality gap := pobj - dobj */
    planarMPC_FLOAT dgap;		
	
    /* relative duality gap := |dgap / pobj | */
    planarMPC_FLOAT rdgap;		

    /* duality measure */
    planarMPC_FLOAT mu;

	/* duality measure (after affine step) */
    planarMPC_FLOAT mu_aff;
	
    /* centering parameter */
    planarMPC_FLOAT sigma;
	
    /* number of backtracking line search steps (affine direction) */
    int lsit_aff;
    
    /* number of backtracking line search steps (combined direction) */
    int lsit_cc;
    
    /* step size (affine direction) */
    planarMPC_FLOAT step_aff;
    
    /* step size (combined direction) */
    planarMPC_FLOAT step_cc;    

	/* solvertime */
	planarMPC_FLOAT solvetime;   

} planarMPC_info;


/* SOLVER FUNCTION DEFINITION -------------------------------------------*/
/* examine exitflag before using the result! */
int planarMPC_solve(planarMPC_params* params, planarMPC_output* output, planarMPC_info* info);


#endif