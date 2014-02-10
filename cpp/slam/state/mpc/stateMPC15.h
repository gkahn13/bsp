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

#ifndef __stateMPC_H__
#define __stateMPC_H__


/* DATA TYPE ------------------------------------------------------------*/
typedef double stateMPC_FLOAT;


/* SOLVER SETTINGS ------------------------------------------------------*/
/* print level */
#ifndef stateMPC_SET_PRINTLEVEL
#define stateMPC_SET_PRINTLEVEL    (2)
#endif

/* timing */
#ifndef stateMPC_SET_TIMING
#define stateMPC_SET_TIMING    (0)
#endif

/* Numeric Warnings */
/* #define PRINTNUMERICALWARNINGS */

/* maximum number of iterations  */
#define stateMPC_SET_MAXIT         (50)	

/* scaling factor of line search (affine direction) */
#define stateMPC_SET_LS_SCALE_AFF  (0.9)      

/* scaling factor of line search (combined direction) */
#define stateMPC_SET_LS_SCALE      (0.95)  

/* minimum required step size in each iteration */
#define stateMPC_SET_LS_MINSTEP    (1E-08)

/* maximum step size (combined direction) */
#define stateMPC_SET_LS_MAXSTEP    (0.995)

/* desired relative duality gap */
#define stateMPC_SET_ACC_RDGAP     (0.0001)

/* desired maximum residual on equality constraints */
#define stateMPC_SET_ACC_RESEQ     (1E-06)

/* desired maximum residual on inequality constraints */
#define stateMPC_SET_ACC_RESINEQ   (1E-06)

/* desired maximum violation of complementarity */
#define stateMPC_SET_ACC_KKTCOMPL  (1E-06)


/* RETURN CODES----------------------------------------------------------*/
/* solver has converged within desired accuracy */
#define stateMPC_OPTIMAL      (1)

/* maximum number of iterations has been reached */
#define stateMPC_MAXITREACHED (0)

/* no progress in line search possible */
#define stateMPC_NOPROGRESS   (-7)




/* PARAMETERS -----------------------------------------------------------*/
/* fill this with data before calling the solver! */
typedef struct stateMPC_params
{
    /* diagonal matrix of size [23 x 23] (only the diagonal is stored) */
    stateMPC_FLOAT H1[23];

    /* vector of size 23 */
    stateMPC_FLOAT f1[23];

    /* vector of size 23 */
    stateMPC_FLOAT lb1[23];

    /* vector of size 9 */
    stateMPC_FLOAT ub1[9];

    /* matrix of size [7 x 23] (column major format) */
    stateMPC_FLOAT C1[161];

    /* vector of size 7 */
    stateMPC_FLOAT e1[7];

    /* diagonal matrix of size [23 x 23] (only the diagonal is stored) */
    stateMPC_FLOAT H2[23];

    /* vector of size 23 */
    stateMPC_FLOAT f2[23];

    /* vector of size 23 */
    stateMPC_FLOAT lb2[23];

    /* vector of size 9 */
    stateMPC_FLOAT ub2[9];

    /* matrix of size [7 x 23] (column major format) */
    stateMPC_FLOAT C2[161];

    /* vector of size 7 */
    stateMPC_FLOAT e2[7];

    /* diagonal matrix of size [23 x 23] (only the diagonal is stored) */
    stateMPC_FLOAT H3[23];

    /* vector of size 23 */
    stateMPC_FLOAT f3[23];

    /* vector of size 23 */
    stateMPC_FLOAT lb3[23];

    /* vector of size 9 */
    stateMPC_FLOAT ub3[9];

    /* matrix of size [7 x 23] (column major format) */
    stateMPC_FLOAT C3[161];

    /* vector of size 7 */
    stateMPC_FLOAT e3[7];

    /* diagonal matrix of size [23 x 23] (only the diagonal is stored) */
    stateMPC_FLOAT H4[23];

    /* vector of size 23 */
    stateMPC_FLOAT f4[23];

    /* vector of size 23 */
    stateMPC_FLOAT lb4[23];

    /* vector of size 9 */
    stateMPC_FLOAT ub4[9];

    /* matrix of size [7 x 23] (column major format) */
    stateMPC_FLOAT C4[161];

    /* vector of size 7 */
    stateMPC_FLOAT e4[7];

    /* diagonal matrix of size [23 x 23] (only the diagonal is stored) */
    stateMPC_FLOAT H5[23];

    /* vector of size 23 */
    stateMPC_FLOAT f5[23];

    /* vector of size 23 */
    stateMPC_FLOAT lb5[23];

    /* vector of size 9 */
    stateMPC_FLOAT ub5[9];

    /* matrix of size [7 x 23] (column major format) */
    stateMPC_FLOAT C5[161];

    /* vector of size 7 */
    stateMPC_FLOAT e5[7];

    /* diagonal matrix of size [23 x 23] (only the diagonal is stored) */
    stateMPC_FLOAT H6[23];

    /* vector of size 23 */
    stateMPC_FLOAT f6[23];

    /* vector of size 23 */
    stateMPC_FLOAT lb6[23];

    /* vector of size 9 */
    stateMPC_FLOAT ub6[9];

    /* matrix of size [7 x 23] (column major format) */
    stateMPC_FLOAT C6[161];

    /* vector of size 7 */
    stateMPC_FLOAT e6[7];

    /* diagonal matrix of size [23 x 23] (only the diagonal is stored) */
    stateMPC_FLOAT H7[23];

    /* vector of size 23 */
    stateMPC_FLOAT f7[23];

    /* vector of size 23 */
    stateMPC_FLOAT lb7[23];

    /* vector of size 9 */
    stateMPC_FLOAT ub7[9];

    /* matrix of size [7 x 23] (column major format) */
    stateMPC_FLOAT C7[161];

    /* vector of size 7 */
    stateMPC_FLOAT e7[7];

    /* diagonal matrix of size [23 x 23] (only the diagonal is stored) */
    stateMPC_FLOAT H8[23];

    /* vector of size 23 */
    stateMPC_FLOAT f8[23];

    /* vector of size 23 */
    stateMPC_FLOAT lb8[23];

    /* vector of size 9 */
    stateMPC_FLOAT ub8[9];

    /* matrix of size [7 x 23] (column major format) */
    stateMPC_FLOAT C8[161];

    /* vector of size 7 */
    stateMPC_FLOAT e8[7];

    /* diagonal matrix of size [23 x 23] (only the diagonal is stored) */
    stateMPC_FLOAT H9[23];

    /* vector of size 23 */
    stateMPC_FLOAT f9[23];

    /* vector of size 23 */
    stateMPC_FLOAT lb9[23];

    /* vector of size 9 */
    stateMPC_FLOAT ub9[9];

    /* matrix of size [7 x 23] (column major format) */
    stateMPC_FLOAT C9[161];

    /* vector of size 7 */
    stateMPC_FLOAT e9[7];

    /* diagonal matrix of size [23 x 23] (only the diagonal is stored) */
    stateMPC_FLOAT H10[23];

    /* vector of size 23 */
    stateMPC_FLOAT f10[23];

    /* vector of size 23 */
    stateMPC_FLOAT lb10[23];

    /* vector of size 9 */
    stateMPC_FLOAT ub10[9];

    /* matrix of size [7 x 23] (column major format) */
    stateMPC_FLOAT C10[161];

    /* vector of size 7 */
    stateMPC_FLOAT e10[7];

    /* diagonal matrix of size [23 x 23] (only the diagonal is stored) */
    stateMPC_FLOAT H11[23];

    /* vector of size 23 */
    stateMPC_FLOAT f11[23];

    /* vector of size 23 */
    stateMPC_FLOAT lb11[23];

    /* vector of size 9 */
    stateMPC_FLOAT ub11[9];

    /* matrix of size [7 x 23] (column major format) */
    stateMPC_FLOAT C11[161];

    /* vector of size 7 */
    stateMPC_FLOAT e11[7];

    /* diagonal matrix of size [23 x 23] (only the diagonal is stored) */
    stateMPC_FLOAT H12[23];

    /* vector of size 23 */
    stateMPC_FLOAT f12[23];

    /* vector of size 23 */
    stateMPC_FLOAT lb12[23];

    /* vector of size 9 */
    stateMPC_FLOAT ub12[9];

    /* matrix of size [7 x 23] (column major format) */
    stateMPC_FLOAT C12[161];

    /* vector of size 7 */
    stateMPC_FLOAT e12[7];

    /* diagonal matrix of size [23 x 23] (only the diagonal is stored) */
    stateMPC_FLOAT H13[23];

    /* vector of size 23 */
    stateMPC_FLOAT f13[23];

    /* vector of size 23 */
    stateMPC_FLOAT lb13[23];

    /* vector of size 9 */
    stateMPC_FLOAT ub13[9];

    /* matrix of size [7 x 23] (column major format) */
    stateMPC_FLOAT C13[161];

    /* vector of size 7 */
    stateMPC_FLOAT e13[7];

    /* diagonal matrix of size [23 x 23] (only the diagonal is stored) */
    stateMPC_FLOAT H14[23];

    /* vector of size 23 */
    stateMPC_FLOAT f14[23];

    /* vector of size 23 */
    stateMPC_FLOAT lb14[23];

    /* vector of size 9 */
    stateMPC_FLOAT ub14[9];

    /* matrix of size [7 x 23] (column major format) */
    stateMPC_FLOAT C14[161];

    /* vector of size 7 */
    stateMPC_FLOAT e14[7];

    /* diagonal matrix of size [7 x 7] (only the diagonal is stored) */
    stateMPC_FLOAT H15[7];

    /* vector of size 7 */
    stateMPC_FLOAT f15[7];

    /* vector of size 7 */
    stateMPC_FLOAT lb15[7];

    /* vector of size 7 */
    stateMPC_FLOAT ub15[7];

    /* vector of size 7 */
    stateMPC_FLOAT e15[7];

} stateMPC_params;


/* OUTPUTS --------------------------------------------------------------*/
/* the desired variables are put here by the solver */
typedef struct stateMPC_output
{
    /* vector of size 9 */
    stateMPC_FLOAT z1[9];

    /* vector of size 9 */
    stateMPC_FLOAT z2[9];

    /* vector of size 9 */
    stateMPC_FLOAT z3[9];

    /* vector of size 9 */
    stateMPC_FLOAT z4[9];

    /* vector of size 9 */
    stateMPC_FLOAT z5[9];

    /* vector of size 9 */
    stateMPC_FLOAT z6[9];

    /* vector of size 9 */
    stateMPC_FLOAT z7[9];

    /* vector of size 9 */
    stateMPC_FLOAT z8[9];

    /* vector of size 9 */
    stateMPC_FLOAT z9[9];

    /* vector of size 9 */
    stateMPC_FLOAT z10[9];

    /* vector of size 9 */
    stateMPC_FLOAT z11[9];

    /* vector of size 9 */
    stateMPC_FLOAT z12[9];

    /* vector of size 9 */
    stateMPC_FLOAT z13[9];

    /* vector of size 9 */
    stateMPC_FLOAT z14[9];

    /* vector of size 7 */
    stateMPC_FLOAT z15[7];

} stateMPC_output;


/* SOLVER INFO ----------------------------------------------------------*/
/* diagnostic data from last interior point step */
typedef struct stateMPC_info
{
    /* iteration number */
    int it;
	
    /* inf-norm of equality constraint residuals */
    stateMPC_FLOAT res_eq;
	
    /* inf-norm of inequality constraint residuals */
    stateMPC_FLOAT res_ineq;

    /* primal objective */
    stateMPC_FLOAT pobj;	
	
    /* dual objective */
    stateMPC_FLOAT dobj;	

    /* duality gap := pobj - dobj */
    stateMPC_FLOAT dgap;		
	
    /* relative duality gap := |dgap / pobj | */
    stateMPC_FLOAT rdgap;		

    /* duality measure */
    stateMPC_FLOAT mu;

	/* duality measure (after affine step) */
    stateMPC_FLOAT mu_aff;
	
    /* centering parameter */
    stateMPC_FLOAT sigma;
	
    /* number of backtracking line search steps (affine direction) */
    int lsit_aff;
    
    /* number of backtracking line search steps (combined direction) */
    int lsit_cc;
    
    /* step size (affine direction) */
    stateMPC_FLOAT step_aff;
    
    /* step size (combined direction) */
    stateMPC_FLOAT step_cc;    

	/* solvertime */
	stateMPC_FLOAT solvetime;   

} stateMPC_info;


/* SOLVER FUNCTION DEFINITION -------------------------------------------*/
/* examine exitflag before using the result! */
int stateMPC_solve(stateMPC_params* params, stateMPC_output* output, stateMPC_info* info);


#endif