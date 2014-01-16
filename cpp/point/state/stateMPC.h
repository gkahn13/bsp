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
#define stateMPC_SET_MAXIT         (30)	

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
    /* diagonal matrix of size [4 x 4] (only the diagonal is stored) */
    stateMPC_FLOAT H01[4];

    /* vector of size 4 */
    stateMPC_FLOAT f01[4];

    /* vector of size 4 */
    stateMPC_FLOAT lb01[4];

    /* vector of size 4 */
    stateMPC_FLOAT ub01[4];

    /* matrix of size [2 x 4] (column major format) */
    stateMPC_FLOAT C01[8];

    /* vector of size 2 */
    stateMPC_FLOAT e01[2];

    /* diagonal matrix of size [4 x 4] (only the diagonal is stored) */
    stateMPC_FLOAT H02[4];

    /* vector of size 4 */
    stateMPC_FLOAT f02[4];

    /* vector of size 4 */
    stateMPC_FLOAT lb02[4];

    /* vector of size 4 */
    stateMPC_FLOAT ub02[4];

    /* matrix of size [2 x 4] (column major format) */
    stateMPC_FLOAT C02[8];

    /* vector of size 2 */
    stateMPC_FLOAT e02[2];

    /* diagonal matrix of size [4 x 4] (only the diagonal is stored) */
    stateMPC_FLOAT H03[4];

    /* vector of size 4 */
    stateMPC_FLOAT f03[4];

    /* vector of size 4 */
    stateMPC_FLOAT lb03[4];

    /* vector of size 4 */
    stateMPC_FLOAT ub03[4];

    /* matrix of size [2 x 4] (column major format) */
    stateMPC_FLOAT C03[8];

    /* vector of size 2 */
    stateMPC_FLOAT e03[2];

    /* diagonal matrix of size [4 x 4] (only the diagonal is stored) */
    stateMPC_FLOAT H04[4];

    /* vector of size 4 */
    stateMPC_FLOAT f04[4];

    /* vector of size 4 */
    stateMPC_FLOAT lb04[4];

    /* vector of size 4 */
    stateMPC_FLOAT ub04[4];

    /* matrix of size [2 x 4] (column major format) */
    stateMPC_FLOAT C04[8];

    /* vector of size 2 */
    stateMPC_FLOAT e04[2];

    /* diagonal matrix of size [4 x 4] (only the diagonal is stored) */
    stateMPC_FLOAT H05[4];

    /* vector of size 4 */
    stateMPC_FLOAT f05[4];

    /* vector of size 4 */
    stateMPC_FLOAT lb05[4];

    /* vector of size 4 */
    stateMPC_FLOAT ub05[4];

    /* matrix of size [2 x 4] (column major format) */
    stateMPC_FLOAT C05[8];

    /* vector of size 2 */
    stateMPC_FLOAT e05[2];

    /* diagonal matrix of size [4 x 4] (only the diagonal is stored) */
    stateMPC_FLOAT H06[4];

    /* vector of size 4 */
    stateMPC_FLOAT f06[4];

    /* vector of size 4 */
    stateMPC_FLOAT lb06[4];

    /* vector of size 4 */
    stateMPC_FLOAT ub06[4];

    /* matrix of size [2 x 4] (column major format) */
    stateMPC_FLOAT C06[8];

    /* vector of size 2 */
    stateMPC_FLOAT e06[2];

    /* diagonal matrix of size [4 x 4] (only the diagonal is stored) */
    stateMPC_FLOAT H07[4];

    /* vector of size 4 */
    stateMPC_FLOAT f07[4];

    /* vector of size 4 */
    stateMPC_FLOAT lb07[4];

    /* vector of size 4 */
    stateMPC_FLOAT ub07[4];

    /* matrix of size [2 x 4] (column major format) */
    stateMPC_FLOAT C07[8];

    /* vector of size 2 */
    stateMPC_FLOAT e07[2];

    /* diagonal matrix of size [4 x 4] (only the diagonal is stored) */
    stateMPC_FLOAT H08[4];

    /* vector of size 4 */
    stateMPC_FLOAT f08[4];

    /* vector of size 4 */
    stateMPC_FLOAT lb08[4];

    /* vector of size 4 */
    stateMPC_FLOAT ub08[4];

    /* matrix of size [2 x 4] (column major format) */
    stateMPC_FLOAT C08[8];

    /* vector of size 2 */
    stateMPC_FLOAT e08[2];

    /* diagonal matrix of size [4 x 4] (only the diagonal is stored) */
    stateMPC_FLOAT H09[4];

    /* vector of size 4 */
    stateMPC_FLOAT f09[4];

    /* vector of size 4 */
    stateMPC_FLOAT lb09[4];

    /* vector of size 4 */
    stateMPC_FLOAT ub09[4];

    /* matrix of size [2 x 4] (column major format) */
    stateMPC_FLOAT C09[8];

    /* vector of size 2 */
    stateMPC_FLOAT e09[2];

    /* diagonal matrix of size [4 x 4] (only the diagonal is stored) */
    stateMPC_FLOAT H10[4];

    /* vector of size 4 */
    stateMPC_FLOAT f10[4];

    /* vector of size 4 */
    stateMPC_FLOAT lb10[4];

    /* vector of size 4 */
    stateMPC_FLOAT ub10[4];

    /* matrix of size [2 x 4] (column major format) */
    stateMPC_FLOAT C10[8];

    /* vector of size 2 */
    stateMPC_FLOAT e10[2];

    /* diagonal matrix of size [4 x 4] (only the diagonal is stored) */
    stateMPC_FLOAT H11[4];

    /* vector of size 4 */
    stateMPC_FLOAT f11[4];

    /* vector of size 4 */
    stateMPC_FLOAT lb11[4];

    /* vector of size 4 */
    stateMPC_FLOAT ub11[4];

    /* matrix of size [2 x 4] (column major format) */
    stateMPC_FLOAT C11[8];

    /* vector of size 2 */
    stateMPC_FLOAT e11[2];

    /* diagonal matrix of size [4 x 4] (only the diagonal is stored) */
    stateMPC_FLOAT H12[4];

    /* vector of size 4 */
    stateMPC_FLOAT f12[4];

    /* vector of size 4 */
    stateMPC_FLOAT lb12[4];

    /* vector of size 4 */
    stateMPC_FLOAT ub12[4];

    /* matrix of size [2 x 4] (column major format) */
    stateMPC_FLOAT C12[8];

    /* vector of size 2 */
    stateMPC_FLOAT e12[2];

    /* diagonal matrix of size [4 x 4] (only the diagonal is stored) */
    stateMPC_FLOAT H13[4];

    /* vector of size 4 */
    stateMPC_FLOAT f13[4];

    /* vector of size 4 */
    stateMPC_FLOAT lb13[4];

    /* vector of size 4 */
    stateMPC_FLOAT ub13[4];

    /* matrix of size [2 x 4] (column major format) */
    stateMPC_FLOAT C13[8];

    /* vector of size 2 */
    stateMPC_FLOAT e13[2];

    /* diagonal matrix of size [4 x 4] (only the diagonal is stored) */
    stateMPC_FLOAT H14[4];

    /* vector of size 4 */
    stateMPC_FLOAT f14[4];

    /* vector of size 4 */
    stateMPC_FLOAT lb14[4];

    /* vector of size 4 */
    stateMPC_FLOAT ub14[4];

    /* matrix of size [2 x 4] (column major format) */
    stateMPC_FLOAT C14[8];

    /* vector of size 2 */
    stateMPC_FLOAT e14[2];

    /* diagonal matrix of size [2 x 2] (only the diagonal is stored) */
    stateMPC_FLOAT H15[2];

    /* vector of size 2 */
    stateMPC_FLOAT f15[2];

    /* vector of size 2 */
    stateMPC_FLOAT lb15[2];

    /* vector of size 2 */
    stateMPC_FLOAT ub15[2];

    /* vector of size 2 */
    stateMPC_FLOAT e15[2];

} stateMPC_params;


/* OUTPUTS --------------------------------------------------------------*/
/* the desired variables are put here by the solver */
typedef struct stateMPC_output
{
    /* vector of size 4 */
    stateMPC_FLOAT z1[4];

    /* vector of size 4 */
    stateMPC_FLOAT z2[4];

    /* vector of size 4 */
    stateMPC_FLOAT z3[4];

    /* vector of size 4 */
    stateMPC_FLOAT z4[4];

    /* vector of size 4 */
    stateMPC_FLOAT z5[4];

    /* vector of size 4 */
    stateMPC_FLOAT z6[4];

    /* vector of size 4 */
    stateMPC_FLOAT z7[4];

    /* vector of size 4 */
    stateMPC_FLOAT z8[4];

    /* vector of size 4 */
    stateMPC_FLOAT z9[4];

    /* vector of size 4 */
    stateMPC_FLOAT z10[4];

    /* vector of size 4 */
    stateMPC_FLOAT z11[4];

    /* vector of size 4 */
    stateMPC_FLOAT z12[4];

    /* vector of size 4 */
    stateMPC_FLOAT z13[4];

    /* vector of size 4 */
    stateMPC_FLOAT z14[4];

    /* vector of size 2 */
    stateMPC_FLOAT z15[2];

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
