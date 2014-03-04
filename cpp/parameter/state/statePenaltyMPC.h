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

#ifndef __statePenaltyMPC_H__
#define __statePenaltyMPC_H__


/* DATA TYPE ------------------------------------------------------------*/
typedef double statePenaltyMPC_FLOAT;


/* SOLVER SETTINGS ------------------------------------------------------*/
/* print level */
#ifndef statePenaltyMPC_SET_PRINTLEVEL
#define statePenaltyMPC_SET_PRINTLEVEL    (0)
#endif

/* timing */
#ifndef statePenaltyMPC_SET_TIMING
#define statePenaltyMPC_SET_TIMING    (0)
#endif

/* Numeric Warnings */
/* #define PRINTNUMERICALWARNINGS */

/* maximum number of iterations  */
#define statePenaltyMPC_SET_MAXIT         (50)	

/* scaling factor of line search (affine direction) */
#define statePenaltyMPC_SET_LS_SCALE_AFF  (0.9)      

/* scaling factor of line search (combined direction) */
#define statePenaltyMPC_SET_LS_SCALE      (0.95)  

/* minimum required step size in each iteration */
#define statePenaltyMPC_SET_LS_MINSTEP    (1E-08)

/* maximum step size (combined direction) */
#define statePenaltyMPC_SET_LS_MAXSTEP    (0.995)

/* desired relative duality gap */
#define statePenaltyMPC_SET_ACC_RDGAP     (0.0001)

/* desired maximum residual on equality constraints */
#define statePenaltyMPC_SET_ACC_RESEQ     (1E-06)

/* desired maximum residual on inequality constraints */
#define statePenaltyMPC_SET_ACC_RESINEQ   (1E-06)

/* desired maximum violation of complementarity */
#define statePenaltyMPC_SET_ACC_KKTCOMPL  (1E-06)


/* RETURN CODES----------------------------------------------------------*/
/* solver has converged within desired accuracy */
#define statePenaltyMPC_OPTIMAL      (1)

/* maximum number of iterations has been reached */
#define statePenaltyMPC_MAXITREACHED (0)

/* no progress in line search possible */
#define statePenaltyMPC_NOPROGRESS   (-7)




/* PARAMETERS -----------------------------------------------------------*/
/* fill this with data before calling the solver! */
typedef struct statePenaltyMPC_params
{
    /* diagonal matrix of size [26 x 26] (only the diagonal is stored) */
    statePenaltyMPC_FLOAT H1[26];

    /* vector of size 26 */
    statePenaltyMPC_FLOAT f1[26];

    /* vector of size 26 */
    statePenaltyMPC_FLOAT lb1[26];

    /* vector of size 10 */
    statePenaltyMPC_FLOAT ub1[10];

    /* matrix of size [8 x 26] (column major format) */
    statePenaltyMPC_FLOAT C1[208];

    /* vector of size 8 */
    statePenaltyMPC_FLOAT e1[8];

    /* diagonal matrix of size [26 x 26] (only the diagonal is stored) */
    statePenaltyMPC_FLOAT H2[26];

    /* vector of size 26 */
    statePenaltyMPC_FLOAT f2[26];

    /* vector of size 26 */
    statePenaltyMPC_FLOAT lb2[26];

    /* vector of size 10 */
    statePenaltyMPC_FLOAT ub2[10];

    /* matrix of size [8 x 26] (column major format) */
    statePenaltyMPC_FLOAT C2[208];

    /* vector of size 8 */
    statePenaltyMPC_FLOAT e2[8];

    /* diagonal matrix of size [26 x 26] (only the diagonal is stored) */
    statePenaltyMPC_FLOAT H3[26];

    /* vector of size 26 */
    statePenaltyMPC_FLOAT f3[26];

    /* vector of size 26 */
    statePenaltyMPC_FLOAT lb3[26];

    /* vector of size 10 */
    statePenaltyMPC_FLOAT ub3[10];

    /* matrix of size [8 x 26] (column major format) */
    statePenaltyMPC_FLOAT C3[208];

    /* vector of size 8 */
    statePenaltyMPC_FLOAT e3[8];

    /* diagonal matrix of size [26 x 26] (only the diagonal is stored) */
    statePenaltyMPC_FLOAT H4[26];

    /* vector of size 26 */
    statePenaltyMPC_FLOAT f4[26];

    /* vector of size 26 */
    statePenaltyMPC_FLOAT lb4[26];

    /* vector of size 10 */
    statePenaltyMPC_FLOAT ub4[10];

    /* matrix of size [8 x 26] (column major format) */
    statePenaltyMPC_FLOAT C4[208];

    /* vector of size 8 */
    statePenaltyMPC_FLOAT e4[8];

    /* diagonal matrix of size [26 x 26] (only the diagonal is stored) */
    statePenaltyMPC_FLOAT H5[26];

    /* vector of size 26 */
    statePenaltyMPC_FLOAT f5[26];

    /* vector of size 26 */
    statePenaltyMPC_FLOAT lb5[26];

    /* vector of size 10 */
    statePenaltyMPC_FLOAT ub5[10];

    /* matrix of size [8 x 26] (column major format) */
    statePenaltyMPC_FLOAT C5[208];

    /* vector of size 8 */
    statePenaltyMPC_FLOAT e5[8];

    /* diagonal matrix of size [26 x 26] (only the diagonal is stored) */
    statePenaltyMPC_FLOAT H6[26];

    /* vector of size 26 */
    statePenaltyMPC_FLOAT f6[26];

    /* vector of size 26 */
    statePenaltyMPC_FLOAT lb6[26];

    /* vector of size 10 */
    statePenaltyMPC_FLOAT ub6[10];

    /* matrix of size [8 x 26] (column major format) */
    statePenaltyMPC_FLOAT C6[208];

    /* vector of size 8 */
    statePenaltyMPC_FLOAT e6[8];

    /* diagonal matrix of size [26 x 26] (only the diagonal is stored) */
    statePenaltyMPC_FLOAT H7[26];

    /* vector of size 26 */
    statePenaltyMPC_FLOAT f7[26];

    /* vector of size 26 */
    statePenaltyMPC_FLOAT lb7[26];

    /* vector of size 10 */
    statePenaltyMPC_FLOAT ub7[10];

    /* matrix of size [8 x 26] (column major format) */
    statePenaltyMPC_FLOAT C7[208];

    /* vector of size 8 */
    statePenaltyMPC_FLOAT e7[8];

    /* diagonal matrix of size [26 x 26] (only the diagonal is stored) */
    statePenaltyMPC_FLOAT H8[26];

    /* vector of size 26 */
    statePenaltyMPC_FLOAT f8[26];

    /* vector of size 26 */
    statePenaltyMPC_FLOAT lb8[26];

    /* vector of size 10 */
    statePenaltyMPC_FLOAT ub8[10];

    /* matrix of size [8 x 26] (column major format) */
    statePenaltyMPC_FLOAT C8[208];

    /* vector of size 8 */
    statePenaltyMPC_FLOAT e8[8];

    /* diagonal matrix of size [26 x 26] (only the diagonal is stored) */
    statePenaltyMPC_FLOAT H9[26];

    /* vector of size 26 */
    statePenaltyMPC_FLOAT f9[26];

    /* vector of size 26 */
    statePenaltyMPC_FLOAT lb9[26];

    /* vector of size 10 */
    statePenaltyMPC_FLOAT ub9[10];

    /* matrix of size [8 x 26] (column major format) */
    statePenaltyMPC_FLOAT C9[208];

    /* vector of size 8 */
    statePenaltyMPC_FLOAT e9[8];

    /* diagonal matrix of size [26 x 26] (only the diagonal is stored) */
    statePenaltyMPC_FLOAT H10[26];

    /* vector of size 26 */
    statePenaltyMPC_FLOAT f10[26];

    /* vector of size 26 */
    statePenaltyMPC_FLOAT lb10[26];

    /* vector of size 10 */
    statePenaltyMPC_FLOAT ub10[10];

    /* matrix of size [8 x 26] (column major format) */
    statePenaltyMPC_FLOAT C10[208];

    /* vector of size 8 */
    statePenaltyMPC_FLOAT e10[8];

    /* diagonal matrix of size [26 x 26] (only the diagonal is stored) */
    statePenaltyMPC_FLOAT H11[26];

    /* vector of size 26 */
    statePenaltyMPC_FLOAT f11[26];

    /* vector of size 26 */
    statePenaltyMPC_FLOAT lb11[26];

    /* vector of size 10 */
    statePenaltyMPC_FLOAT ub11[10];

    /* matrix of size [8 x 26] (column major format) */
    statePenaltyMPC_FLOAT C11[208];

    /* vector of size 8 */
    statePenaltyMPC_FLOAT e11[8];

    /* diagonal matrix of size [26 x 26] (only the diagonal is stored) */
    statePenaltyMPC_FLOAT H12[26];

    /* vector of size 26 */
    statePenaltyMPC_FLOAT f12[26];

    /* vector of size 26 */
    statePenaltyMPC_FLOAT lb12[26];

    /* vector of size 10 */
    statePenaltyMPC_FLOAT ub12[10];

    /* matrix of size [8 x 26] (column major format) */
    statePenaltyMPC_FLOAT C12[208];

    /* vector of size 8 */
    statePenaltyMPC_FLOAT e12[8];

    /* diagonal matrix of size [26 x 26] (only the diagonal is stored) */
    statePenaltyMPC_FLOAT H13[26];

    /* vector of size 26 */
    statePenaltyMPC_FLOAT f13[26];

    /* vector of size 26 */
    statePenaltyMPC_FLOAT lb13[26];

    /* vector of size 10 */
    statePenaltyMPC_FLOAT ub13[10];

    /* matrix of size [8 x 26] (column major format) */
    statePenaltyMPC_FLOAT C13[208];

    /* vector of size 8 */
    statePenaltyMPC_FLOAT e13[8];

    /* diagonal matrix of size [26 x 26] (only the diagonal is stored) */
    statePenaltyMPC_FLOAT H14[26];

    /* vector of size 26 */
    statePenaltyMPC_FLOAT f14[26];

    /* vector of size 26 */
    statePenaltyMPC_FLOAT lb14[26];

    /* vector of size 10 */
    statePenaltyMPC_FLOAT ub14[10];

    /* matrix of size [8 x 26] (column major format) */
    statePenaltyMPC_FLOAT C14[208];

    /* vector of size 8 */
    statePenaltyMPC_FLOAT e14[8];

    /* diagonal matrix of size [8 x 8] (only the diagonal is stored) */
    statePenaltyMPC_FLOAT H15[8];

    /* vector of size 8 */
    statePenaltyMPC_FLOAT f15[8];

    /* vector of size 8 */
    statePenaltyMPC_FLOAT lb15[8];

    /* vector of size 8 */
    statePenaltyMPC_FLOAT ub15[8];

    /* vector of size 8 */
    statePenaltyMPC_FLOAT e15[8];

} statePenaltyMPC_params;


/* OUTPUTS --------------------------------------------------------------*/
/* the desired variables are put here by the solver */
typedef struct statePenaltyMPC_output
{
    /* vector of size 10 */
    statePenaltyMPC_FLOAT z1[10];

    /* vector of size 10 */
    statePenaltyMPC_FLOAT z2[10];

    /* vector of size 10 */
    statePenaltyMPC_FLOAT z3[10];

    /* vector of size 10 */
    statePenaltyMPC_FLOAT z4[10];

    /* vector of size 10 */
    statePenaltyMPC_FLOAT z5[10];

    /* vector of size 10 */
    statePenaltyMPC_FLOAT z6[10];

    /* vector of size 10 */
    statePenaltyMPC_FLOAT z7[10];

    /* vector of size 10 */
    statePenaltyMPC_FLOAT z8[10];

    /* vector of size 10 */
    statePenaltyMPC_FLOAT z9[10];

    /* vector of size 10 */
    statePenaltyMPC_FLOAT z10[10];

    /* vector of size 10 */
    statePenaltyMPC_FLOAT z11[10];

    /* vector of size 10 */
    statePenaltyMPC_FLOAT z12[10];

    /* vector of size 10 */
    statePenaltyMPC_FLOAT z13[10];

    /* vector of size 10 */
    statePenaltyMPC_FLOAT z14[10];

    /* vector of size 8 */
    statePenaltyMPC_FLOAT z15[8];

} statePenaltyMPC_output;


/* SOLVER INFO ----------------------------------------------------------*/
/* diagnostic data from last interior point step */
typedef struct statePenaltyMPC_info
{
    /* iteration number */
    int it;
	
    /* inf-norm of equality constraint residuals */
    statePenaltyMPC_FLOAT res_eq;
	
    /* inf-norm of inequality constraint residuals */
    statePenaltyMPC_FLOAT res_ineq;

    /* primal objective */
    statePenaltyMPC_FLOAT pobj;	
	
    /* dual objective */
    statePenaltyMPC_FLOAT dobj;	

    /* duality gap := pobj - dobj */
    statePenaltyMPC_FLOAT dgap;		
	
    /* relative duality gap := |dgap / pobj | */
    statePenaltyMPC_FLOAT rdgap;		

    /* duality measure */
    statePenaltyMPC_FLOAT mu;

	/* duality measure (after affine step) */
    statePenaltyMPC_FLOAT mu_aff;
	
    /* centering parameter */
    statePenaltyMPC_FLOAT sigma;
	
    /* number of backtracking line search steps (affine direction) */
    int lsit_aff;
    
    /* number of backtracking line search steps (combined direction) */
    int lsit_cc;
    
    /* step size (affine direction) */
    statePenaltyMPC_FLOAT step_aff;
    
    /* step size (combined direction) */
    statePenaltyMPC_FLOAT step_cc;    

	/* solvertime */
	statePenaltyMPC_FLOAT solvetime;   

} statePenaltyMPC_info;


/* SOLVER FUNCTION DEFINITION -------------------------------------------*/
/* examine exitflag before using the result! */
int statePenaltyMPC_solve(statePenaltyMPC_params* params, statePenaltyMPC_output* output, statePenaltyMPC_info* info);


#endif