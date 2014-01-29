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
#define beliefPenaltyMPC_SET_PRINTLEVEL    (0)
#endif

/* timing */
#ifndef beliefPenaltyMPC_SET_TIMING
#define beliefPenaltyMPC_SET_TIMING    (0)
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
    /* diagonal matrix of size [197 x 197] (only the diagonal is stored) */
    beliefPenaltyMPC_FLOAT H1[197];

    /* vector of size 197 */
    beliefPenaltyMPC_FLOAT f1[197];

    /* vector of size 197 */
    beliefPenaltyMPC_FLOAT lb1[197];

    /* vector of size 67 */
    beliefPenaltyMPC_FLOAT ub1[67];

    /* matrix of size [130 x 197] (column major format) */
    beliefPenaltyMPC_FLOAT C1[25610];

    /* vector of size 130 */
    beliefPenaltyMPC_FLOAT e1[130];

    /* diagonal matrix of size [197 x 197] (only the diagonal is stored) */
    beliefPenaltyMPC_FLOAT H2[197];

    /* vector of size 197 */
    beliefPenaltyMPC_FLOAT f2[197];

    /* vector of size 197 */
    beliefPenaltyMPC_FLOAT lb2[197];

    /* vector of size 67 */
    beliefPenaltyMPC_FLOAT ub2[67];

    /* matrix of size [65 x 197] (column major format) */
    beliefPenaltyMPC_FLOAT C2[12805];

    /* vector of size 65 */
    beliefPenaltyMPC_FLOAT e2[65];

    /* diagonal matrix of size [197 x 197] (only the diagonal is stored) */
    beliefPenaltyMPC_FLOAT H3[197];

    /* vector of size 197 */
    beliefPenaltyMPC_FLOAT f3[197];

    /* vector of size 197 */
    beliefPenaltyMPC_FLOAT lb3[197];

    /* vector of size 67 */
    beliefPenaltyMPC_FLOAT ub3[67];

    /* matrix of size [65 x 197] (column major format) */
    beliefPenaltyMPC_FLOAT C3[12805];

    /* vector of size 65 */
    beliefPenaltyMPC_FLOAT e3[65];

    /* diagonal matrix of size [197 x 197] (only the diagonal is stored) */
    beliefPenaltyMPC_FLOAT H4[197];

    /* vector of size 197 */
    beliefPenaltyMPC_FLOAT f4[197];

    /* vector of size 197 */
    beliefPenaltyMPC_FLOAT lb4[197];

    /* vector of size 67 */
    beliefPenaltyMPC_FLOAT ub4[67];

    /* matrix of size [65 x 197] (column major format) */
    beliefPenaltyMPC_FLOAT C4[12805];

    /* vector of size 65 */
    beliefPenaltyMPC_FLOAT e4[65];

    /* diagonal matrix of size [197 x 197] (only the diagonal is stored) */
    beliefPenaltyMPC_FLOAT H5[197];

    /* vector of size 197 */
    beliefPenaltyMPC_FLOAT f5[197];

    /* vector of size 197 */
    beliefPenaltyMPC_FLOAT lb5[197];

    /* vector of size 67 */
    beliefPenaltyMPC_FLOAT ub5[67];

    /* matrix of size [65 x 197] (column major format) */
    beliefPenaltyMPC_FLOAT C5[12805];

    /* vector of size 65 */
    beliefPenaltyMPC_FLOAT e5[65];

    /* diagonal matrix of size [197 x 197] (only the diagonal is stored) */
    beliefPenaltyMPC_FLOAT H6[197];

    /* vector of size 197 */
    beliefPenaltyMPC_FLOAT f6[197];

    /* vector of size 197 */
    beliefPenaltyMPC_FLOAT lb6[197];

    /* vector of size 67 */
    beliefPenaltyMPC_FLOAT ub6[67];

    /* matrix of size [65 x 197] (column major format) */
    beliefPenaltyMPC_FLOAT C6[12805];

    /* vector of size 65 */
    beliefPenaltyMPC_FLOAT e6[65];

    /* diagonal matrix of size [197 x 197] (only the diagonal is stored) */
    beliefPenaltyMPC_FLOAT H7[197];

    /* vector of size 197 */
    beliefPenaltyMPC_FLOAT f7[197];

    /* vector of size 197 */
    beliefPenaltyMPC_FLOAT lb7[197];

    /* vector of size 67 */
    beliefPenaltyMPC_FLOAT ub7[67];

    /* matrix of size [65 x 197] (column major format) */
    beliefPenaltyMPC_FLOAT C7[12805];

    /* vector of size 65 */
    beliefPenaltyMPC_FLOAT e7[65];

    /* diagonal matrix of size [197 x 197] (only the diagonal is stored) */
    beliefPenaltyMPC_FLOAT H8[197];

    /* vector of size 197 */
    beliefPenaltyMPC_FLOAT f8[197];

    /* vector of size 197 */
    beliefPenaltyMPC_FLOAT lb8[197];

    /* vector of size 67 */
    beliefPenaltyMPC_FLOAT ub8[67];

    /* matrix of size [65 x 197] (column major format) */
    beliefPenaltyMPC_FLOAT C8[12805];

    /* vector of size 65 */
    beliefPenaltyMPC_FLOAT e8[65];

    /* diagonal matrix of size [197 x 197] (only the diagonal is stored) */
    beliefPenaltyMPC_FLOAT H9[197];

    /* vector of size 197 */
    beliefPenaltyMPC_FLOAT f9[197];

    /* vector of size 197 */
    beliefPenaltyMPC_FLOAT lb9[197];

    /* vector of size 67 */
    beliefPenaltyMPC_FLOAT ub9[67];

    /* matrix of size [65 x 197] (column major format) */
    beliefPenaltyMPC_FLOAT C9[12805];

    /* vector of size 65 */
    beliefPenaltyMPC_FLOAT e9[65];

    /* diagonal matrix of size [197 x 197] (only the diagonal is stored) */
    beliefPenaltyMPC_FLOAT H10[197];

    /* vector of size 197 */
    beliefPenaltyMPC_FLOAT f10[197];

    /* vector of size 197 */
    beliefPenaltyMPC_FLOAT lb10[197];

    /* vector of size 67 */
    beliefPenaltyMPC_FLOAT ub10[67];

    /* matrix of size [65 x 197] (column major format) */
    beliefPenaltyMPC_FLOAT C10[12805];

    /* vector of size 65 */
    beliefPenaltyMPC_FLOAT e10[65];

    /* diagonal matrix of size [197 x 197] (only the diagonal is stored) */
    beliefPenaltyMPC_FLOAT H11[197];

    /* vector of size 197 */
    beliefPenaltyMPC_FLOAT f11[197];

    /* vector of size 197 */
    beliefPenaltyMPC_FLOAT lb11[197];

    /* vector of size 67 */
    beliefPenaltyMPC_FLOAT ub11[67];

    /* matrix of size [65 x 197] (column major format) */
    beliefPenaltyMPC_FLOAT C11[12805];

    /* vector of size 65 */
    beliefPenaltyMPC_FLOAT e11[65];

    /* diagonal matrix of size [197 x 197] (only the diagonal is stored) */
    beliefPenaltyMPC_FLOAT H12[197];

    /* vector of size 197 */
    beliefPenaltyMPC_FLOAT f12[197];

    /* vector of size 197 */
    beliefPenaltyMPC_FLOAT lb12[197];

    /* vector of size 67 */
    beliefPenaltyMPC_FLOAT ub12[67];

    /* matrix of size [65 x 197] (column major format) */
    beliefPenaltyMPC_FLOAT C12[12805];

    /* vector of size 65 */
    beliefPenaltyMPC_FLOAT e12[65];

    /* diagonal matrix of size [197 x 197] (only the diagonal is stored) */
    beliefPenaltyMPC_FLOAT H13[197];

    /* vector of size 197 */
    beliefPenaltyMPC_FLOAT f13[197];

    /* vector of size 197 */
    beliefPenaltyMPC_FLOAT lb13[197];

    /* vector of size 67 */
    beliefPenaltyMPC_FLOAT ub13[67];

    /* matrix of size [65 x 197] (column major format) */
    beliefPenaltyMPC_FLOAT C13[12805];

    /* vector of size 65 */
    beliefPenaltyMPC_FLOAT e13[65];

    /* diagonal matrix of size [197 x 197] (only the diagonal is stored) */
    beliefPenaltyMPC_FLOAT H14[197];

    /* vector of size 197 */
    beliefPenaltyMPC_FLOAT f14[197];

    /* vector of size 197 */
    beliefPenaltyMPC_FLOAT lb14[197];

    /* vector of size 67 */
    beliefPenaltyMPC_FLOAT ub14[67];

    /* matrix of size [65 x 197] (column major format) */
    beliefPenaltyMPC_FLOAT C14[12805];

    /* vector of size 65 */
    beliefPenaltyMPC_FLOAT e14[65];

    /* diagonal matrix of size [65 x 65] (only the diagonal is stored) */
    beliefPenaltyMPC_FLOAT H15[65];

    /* vector of size 65 */
    beliefPenaltyMPC_FLOAT lb15[65];

    /* vector of size 65 */
    beliefPenaltyMPC_FLOAT ub15[65];

} beliefPenaltyMPC_params;


/* OUTPUTS --------------------------------------------------------------*/
/* the desired variables are put here by the solver */
typedef struct beliefPenaltyMPC_output
{
    /* vector of size 67 */
    beliefPenaltyMPC_FLOAT z1[67];

    /* vector of size 67 */
    beliefPenaltyMPC_FLOAT z2[67];

    /* vector of size 67 */
    beliefPenaltyMPC_FLOAT z3[67];

    /* vector of size 67 */
    beliefPenaltyMPC_FLOAT z4[67];

    /* vector of size 67 */
    beliefPenaltyMPC_FLOAT z5[67];

    /* vector of size 67 */
    beliefPenaltyMPC_FLOAT z6[67];

    /* vector of size 67 */
    beliefPenaltyMPC_FLOAT z7[67];

    /* vector of size 67 */
    beliefPenaltyMPC_FLOAT z8[67];

    /* vector of size 67 */
    beliefPenaltyMPC_FLOAT z9[67];

    /* vector of size 67 */
    beliefPenaltyMPC_FLOAT z10[67];

    /* vector of size 67 */
    beliefPenaltyMPC_FLOAT z11[67];

    /* vector of size 67 */
    beliefPenaltyMPC_FLOAT z12[67];

    /* vector of size 67 */
    beliefPenaltyMPC_FLOAT z13[67];

    /* vector of size 67 */
    beliefPenaltyMPC_FLOAT z14[67];

    /* vector of size 65 */
    beliefPenaltyMPC_FLOAT z15[65];

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