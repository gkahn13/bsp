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
#define stateMPC_SET_PRINTLEVEL    (1)
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
    /* diagonal matrix of size [29 x 29] (only the diagonal is stored) */
    stateMPC_FLOAT H1[29];

    /* vector of size 29 */
    stateMPC_FLOAT f1[29];

    /* vector of size 29 */
    stateMPC_FLOAT lb1[29];

    /* vector of size 11 */
    stateMPC_FLOAT ub1[11];

    /* matrix of size [9 x 29] (column major format) */
    stateMPC_FLOAT C1[261];

    /* vector of size 9 */
    stateMPC_FLOAT e1[9];

    /* diagonal matrix of size [29 x 29] (only the diagonal is stored) */
    stateMPC_FLOAT H2[29];

    /* vector of size 29 */
    stateMPC_FLOAT f2[29];

    /* vector of size 29 */
    stateMPC_FLOAT lb2[29];

    /* vector of size 11 */
    stateMPC_FLOAT ub2[11];

    /* matrix of size [9 x 29] (column major format) */
    stateMPC_FLOAT C2[261];

    /* vector of size 9 */
    stateMPC_FLOAT e2[9];

    /* diagonal matrix of size [29 x 29] (only the diagonal is stored) */
    stateMPC_FLOAT H3[29];

    /* vector of size 29 */
    stateMPC_FLOAT f3[29];

    /* vector of size 29 */
    stateMPC_FLOAT lb3[29];

    /* vector of size 11 */
    stateMPC_FLOAT ub3[11];

    /* matrix of size [9 x 29] (column major format) */
    stateMPC_FLOAT C3[261];

    /* vector of size 9 */
    stateMPC_FLOAT e3[9];

    /* diagonal matrix of size [29 x 29] (only the diagonal is stored) */
    stateMPC_FLOAT H4[29];

    /* vector of size 29 */
    stateMPC_FLOAT f4[29];

    /* vector of size 29 */
    stateMPC_FLOAT lb4[29];

    /* vector of size 11 */
    stateMPC_FLOAT ub4[11];

    /* matrix of size [9 x 29] (column major format) */
    stateMPC_FLOAT C4[261];

    /* vector of size 9 */
    stateMPC_FLOAT e4[9];

    /* diagonal matrix of size [29 x 29] (only the diagonal is stored) */
    stateMPC_FLOAT H5[29];

    /* vector of size 29 */
    stateMPC_FLOAT f5[29];

    /* vector of size 29 */
    stateMPC_FLOAT lb5[29];

    /* vector of size 11 */
    stateMPC_FLOAT ub5[11];

    /* matrix of size [9 x 29] (column major format) */
    stateMPC_FLOAT C5[261];

    /* vector of size 9 */
    stateMPC_FLOAT e5[9];

    /* diagonal matrix of size [29 x 29] (only the diagonal is stored) */
    stateMPC_FLOAT H6[29];

    /* vector of size 29 */
    stateMPC_FLOAT f6[29];

    /* vector of size 29 */
    stateMPC_FLOAT lb6[29];

    /* vector of size 11 */
    stateMPC_FLOAT ub6[11];

    /* matrix of size [9 x 29] (column major format) */
    stateMPC_FLOAT C6[261];

    /* vector of size 9 */
    stateMPC_FLOAT e6[9];

    /* diagonal matrix of size [29 x 29] (only the diagonal is stored) */
    stateMPC_FLOAT H7[29];

    /* vector of size 29 */
    stateMPC_FLOAT f7[29];

    /* vector of size 29 */
    stateMPC_FLOAT lb7[29];

    /* vector of size 11 */
    stateMPC_FLOAT ub7[11];

    /* matrix of size [9 x 29] (column major format) */
    stateMPC_FLOAT C7[261];

    /* vector of size 9 */
    stateMPC_FLOAT e7[9];

    /* diagonal matrix of size [29 x 29] (only the diagonal is stored) */
    stateMPC_FLOAT H8[29];

    /* vector of size 29 */
    stateMPC_FLOAT f8[29];

    /* vector of size 29 */
    stateMPC_FLOAT lb8[29];

    /* vector of size 11 */
    stateMPC_FLOAT ub8[11];

    /* matrix of size [9 x 29] (column major format) */
    stateMPC_FLOAT C8[261];

    /* vector of size 9 */
    stateMPC_FLOAT e8[9];

    /* diagonal matrix of size [29 x 29] (only the diagonal is stored) */
    stateMPC_FLOAT H9[29];

    /* vector of size 29 */
    stateMPC_FLOAT f9[29];

    /* vector of size 29 */
    stateMPC_FLOAT lb9[29];

    /* vector of size 11 */
    stateMPC_FLOAT ub9[11];

    /* matrix of size [9 x 29] (column major format) */
    stateMPC_FLOAT C9[261];

    /* vector of size 9 */
    stateMPC_FLOAT e9[9];

    /* diagonal matrix of size [29 x 29] (only the diagonal is stored) */
    stateMPC_FLOAT H10[29];

    /* vector of size 29 */
    stateMPC_FLOAT f10[29];

    /* vector of size 29 */
    stateMPC_FLOAT lb10[29];

    /* vector of size 11 */
    stateMPC_FLOAT ub10[11];

    /* matrix of size [9 x 29] (column major format) */
    stateMPC_FLOAT C10[261];

    /* vector of size 9 */
    stateMPC_FLOAT e10[9];

    /* diagonal matrix of size [29 x 29] (only the diagonal is stored) */
    stateMPC_FLOAT H11[29];

    /* vector of size 29 */
    stateMPC_FLOAT f11[29];

    /* vector of size 29 */
    stateMPC_FLOAT lb11[29];

    /* vector of size 11 */
    stateMPC_FLOAT ub11[11];

    /* matrix of size [9 x 29] (column major format) */
    stateMPC_FLOAT C11[261];

    /* vector of size 9 */
    stateMPC_FLOAT e11[9];

    /* diagonal matrix of size [29 x 29] (only the diagonal is stored) */
    stateMPC_FLOAT H12[29];

    /* vector of size 29 */
    stateMPC_FLOAT f12[29];

    /* vector of size 29 */
    stateMPC_FLOAT lb12[29];

    /* vector of size 11 */
    stateMPC_FLOAT ub12[11];

    /* matrix of size [9 x 29] (column major format) */
    stateMPC_FLOAT C12[261];

    /* vector of size 9 */
    stateMPC_FLOAT e12[9];

    /* diagonal matrix of size [29 x 29] (only the diagonal is stored) */
    stateMPC_FLOAT H13[29];

    /* vector of size 29 */
    stateMPC_FLOAT f13[29];

    /* vector of size 29 */
    stateMPC_FLOAT lb13[29];

    /* vector of size 11 */
    stateMPC_FLOAT ub13[11];

    /* matrix of size [9 x 29] (column major format) */
    stateMPC_FLOAT C13[261];

    /* vector of size 9 */
    stateMPC_FLOAT e13[9];

    /* diagonal matrix of size [29 x 29] (only the diagonal is stored) */
    stateMPC_FLOAT H14[29];

    /* vector of size 29 */
    stateMPC_FLOAT f14[29];

    /* vector of size 29 */
    stateMPC_FLOAT lb14[29];

    /* vector of size 11 */
    stateMPC_FLOAT ub14[11];

    /* matrix of size [9 x 29] (column major format) */
    stateMPC_FLOAT C14[261];

    /* vector of size 9 */
    stateMPC_FLOAT e14[9];

    /* diagonal matrix of size [9 x 9] (only the diagonal is stored) */
    stateMPC_FLOAT H15[9];

    /* vector of size 9 */
    stateMPC_FLOAT f15[9];

    /* vector of size 9 */
    stateMPC_FLOAT lb15[9];

    /* vector of size 9 */
    stateMPC_FLOAT ub15[9];

    /* vector of size 9 */
    stateMPC_FLOAT e15[9];

} stateMPC_params;


/* OUTPUTS --------------------------------------------------------------*/
/* the desired variables are put here by the solver */
typedef struct stateMPC_output
{
    /* vector of size 11 */
    stateMPC_FLOAT z1[11];

    /* vector of size 11 */
    stateMPC_FLOAT z2[11];

    /* vector of size 11 */
    stateMPC_FLOAT z3[11];

    /* vector of size 11 */
    stateMPC_FLOAT z4[11];

    /* vector of size 11 */
    stateMPC_FLOAT z5[11];

    /* vector of size 11 */
    stateMPC_FLOAT z6[11];

    /* vector of size 11 */
    stateMPC_FLOAT z7[11];

    /* vector of size 11 */
    stateMPC_FLOAT z8[11];

    /* vector of size 11 */
    stateMPC_FLOAT z9[11];

    /* vector of size 11 */
    stateMPC_FLOAT z10[11];

    /* vector of size 11 */
    stateMPC_FLOAT z11[11];

    /* vector of size 11 */
    stateMPC_FLOAT z12[11];

    /* vector of size 11 */
    stateMPC_FLOAT z13[11];

    /* vector of size 11 */
    stateMPC_FLOAT z14[11];

    /* vector of size 9 */
    stateMPC_FLOAT z15[9];

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