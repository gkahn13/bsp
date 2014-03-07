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
#define stateMPC_SET_PRINTLEVEL    (0)
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
    /* diagonal matrix of size [11 x 11] (only the diagonal is stored) */
    stateMPC_FLOAT H1[11];

    /* vector of size 11 */
    stateMPC_FLOAT f1[11];

    /* vector of size 11 */
    stateMPC_FLOAT lb1[11];

    /* vector of size 5 */
    stateMPC_FLOAT ub1[5];

    /* matrix of size [3 x 11] (column major format) */
    stateMPC_FLOAT C1[33];

    /* vector of size 3 */
    stateMPC_FLOAT e1[3];

    /* diagonal matrix of size [11 x 11] (only the diagonal is stored) */
    stateMPC_FLOAT H2[11];

    /* vector of size 11 */
    stateMPC_FLOAT f2[11];

    /* vector of size 11 */
    stateMPC_FLOAT lb2[11];

    /* vector of size 5 */
    stateMPC_FLOAT ub2[5];

    /* matrix of size [3 x 11] (column major format) */
    stateMPC_FLOAT C2[33];

    /* vector of size 3 */
    stateMPC_FLOAT e2[3];

    /* diagonal matrix of size [11 x 11] (only the diagonal is stored) */
    stateMPC_FLOAT H3[11];

    /* vector of size 11 */
    stateMPC_FLOAT f3[11];

    /* vector of size 11 */
    stateMPC_FLOAT lb3[11];

    /* vector of size 5 */
    stateMPC_FLOAT ub3[5];

    /* matrix of size [3 x 11] (column major format) */
    stateMPC_FLOAT C3[33];

    /* vector of size 3 */
    stateMPC_FLOAT e3[3];

    /* diagonal matrix of size [11 x 11] (only the diagonal is stored) */
    stateMPC_FLOAT H4[11];

    /* vector of size 11 */
    stateMPC_FLOAT f4[11];

    /* vector of size 11 */
    stateMPC_FLOAT lb4[11];

    /* vector of size 5 */
    stateMPC_FLOAT ub4[5];

    /* matrix of size [3 x 11] (column major format) */
    stateMPC_FLOAT C4[33];

    /* vector of size 3 */
    stateMPC_FLOAT e4[3];

    /* diagonal matrix of size [11 x 11] (only the diagonal is stored) */
    stateMPC_FLOAT H5[11];

    /* vector of size 11 */
    stateMPC_FLOAT f5[11];

    /* vector of size 11 */
    stateMPC_FLOAT lb5[11];

    /* vector of size 5 */
    stateMPC_FLOAT ub5[5];

    /* matrix of size [3 x 11] (column major format) */
    stateMPC_FLOAT C5[33];

    /* vector of size 3 */
    stateMPC_FLOAT e5[3];

    /* diagonal matrix of size [11 x 11] (only the diagonal is stored) */
    stateMPC_FLOAT H6[11];

    /* vector of size 11 */
    stateMPC_FLOAT f6[11];

    /* vector of size 11 */
    stateMPC_FLOAT lb6[11];

    /* vector of size 5 */
    stateMPC_FLOAT ub6[5];

    /* matrix of size [3 x 11] (column major format) */
    stateMPC_FLOAT C6[33];

    /* vector of size 3 */
    stateMPC_FLOAT e6[3];

    /* diagonal matrix of size [11 x 11] (only the diagonal is stored) */
    stateMPC_FLOAT H7[11];

    /* vector of size 11 */
    stateMPC_FLOAT f7[11];

    /* vector of size 11 */
    stateMPC_FLOAT lb7[11];

    /* vector of size 5 */
    stateMPC_FLOAT ub7[5];

    /* matrix of size [3 x 11] (column major format) */
    stateMPC_FLOAT C7[33];

    /* vector of size 3 */
    stateMPC_FLOAT e7[3];

    /* diagonal matrix of size [11 x 11] (only the diagonal is stored) */
    stateMPC_FLOAT H8[11];

    /* vector of size 11 */
    stateMPC_FLOAT f8[11];

    /* vector of size 11 */
    stateMPC_FLOAT lb8[11];

    /* vector of size 5 */
    stateMPC_FLOAT ub8[5];

    /* matrix of size [3 x 11] (column major format) */
    stateMPC_FLOAT C8[33];

    /* vector of size 3 */
    stateMPC_FLOAT e8[3];

    /* diagonal matrix of size [11 x 11] (only the diagonal is stored) */
    stateMPC_FLOAT H9[11];

    /* vector of size 11 */
    stateMPC_FLOAT f9[11];

    /* vector of size 11 */
    stateMPC_FLOAT lb9[11];

    /* vector of size 5 */
    stateMPC_FLOAT ub9[5];

    /* matrix of size [3 x 11] (column major format) */
    stateMPC_FLOAT C9[33];

    /* vector of size 3 */
    stateMPC_FLOAT e9[3];

    /* diagonal matrix of size [11 x 11] (only the diagonal is stored) */
    stateMPC_FLOAT H10[11];

    /* vector of size 11 */
    stateMPC_FLOAT f10[11];

    /* vector of size 11 */
    stateMPC_FLOAT lb10[11];

    /* vector of size 5 */
    stateMPC_FLOAT ub10[5];

    /* matrix of size [3 x 11] (column major format) */
    stateMPC_FLOAT C10[33];

    /* vector of size 3 */
    stateMPC_FLOAT e10[3];

    /* diagonal matrix of size [11 x 11] (only the diagonal is stored) */
    stateMPC_FLOAT H11[11];

    /* vector of size 11 */
    stateMPC_FLOAT f11[11];

    /* vector of size 11 */
    stateMPC_FLOAT lb11[11];

    /* vector of size 5 */
    stateMPC_FLOAT ub11[5];

    /* matrix of size [3 x 11] (column major format) */
    stateMPC_FLOAT C11[33];

    /* vector of size 3 */
    stateMPC_FLOAT e11[3];

    /* diagonal matrix of size [11 x 11] (only the diagonal is stored) */
    stateMPC_FLOAT H12[11];

    /* vector of size 11 */
    stateMPC_FLOAT f12[11];

    /* vector of size 11 */
    stateMPC_FLOAT lb12[11];

    /* vector of size 5 */
    stateMPC_FLOAT ub12[5];

    /* matrix of size [3 x 11] (column major format) */
    stateMPC_FLOAT C12[33];

    /* vector of size 3 */
    stateMPC_FLOAT e12[3];

    /* diagonal matrix of size [11 x 11] (only the diagonal is stored) */
    stateMPC_FLOAT H13[11];

    /* vector of size 11 */
    stateMPC_FLOAT f13[11];

    /* vector of size 11 */
    stateMPC_FLOAT lb13[11];

    /* vector of size 5 */
    stateMPC_FLOAT ub13[5];

    /* matrix of size [3 x 11] (column major format) */
    stateMPC_FLOAT C13[33];

    /* vector of size 3 */
    stateMPC_FLOAT e13[3];

    /* diagonal matrix of size [11 x 11] (only the diagonal is stored) */
    stateMPC_FLOAT H14[11];

    /* vector of size 11 */
    stateMPC_FLOAT f14[11];

    /* vector of size 11 */
    stateMPC_FLOAT lb14[11];

    /* vector of size 5 */
    stateMPC_FLOAT ub14[5];

    /* matrix of size [3 x 11] (column major format) */
    stateMPC_FLOAT C14[33];

    /* vector of size 3 */
    stateMPC_FLOAT e14[3];

    /* diagonal matrix of size [3 x 3] (only the diagonal is stored) */
    stateMPC_FLOAT H15[3];

    /* vector of size 3 */
    stateMPC_FLOAT f15[3];

    /* vector of size 3 */
    stateMPC_FLOAT lb15[3];

    /* vector of size 3 */
    stateMPC_FLOAT ub15[3];

    /* vector of size 3 */
    stateMPC_FLOAT e15[3];

} stateMPC_params;


/* OUTPUTS --------------------------------------------------------------*/
/* the desired variables are put here by the solver */
typedef struct stateMPC_output
{
    /* vector of size 5 */
    stateMPC_FLOAT z1[5];

    /* vector of size 5 */
    stateMPC_FLOAT z2[5];

    /* vector of size 5 */
    stateMPC_FLOAT z3[5];

    /* vector of size 5 */
    stateMPC_FLOAT z4[5];

    /* vector of size 5 */
    stateMPC_FLOAT z5[5];

    /* vector of size 5 */
    stateMPC_FLOAT z6[5];

    /* vector of size 5 */
    stateMPC_FLOAT z7[5];

    /* vector of size 5 */
    stateMPC_FLOAT z8[5];

    /* vector of size 5 */
    stateMPC_FLOAT z9[5];

    /* vector of size 5 */
    stateMPC_FLOAT z10[5];

    /* vector of size 5 */
    stateMPC_FLOAT z11[5];

    /* vector of size 5 */
    stateMPC_FLOAT z12[5];

    /* vector of size 5 */
    stateMPC_FLOAT z13[5];

    /* vector of size 5 */
    stateMPC_FLOAT z14[5];

    /* vector of size 3 */
    stateMPC_FLOAT z15[3];

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