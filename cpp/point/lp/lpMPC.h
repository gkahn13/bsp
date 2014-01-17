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

#ifndef __lpMPC_H__
#define __lpMPC_H__


/* DATA TYPE ------------------------------------------------------------*/
typedef double lpMPC_FLOAT;


/* SOLVER SETTINGS ------------------------------------------------------*/
/* print level */
#ifndef lpMPC_SET_PRINTLEVEL
#define lpMPC_SET_PRINTLEVEL    (0)
#endif

/* timing */
#ifndef lpMPC_SET_TIMING
#define lpMPC_SET_TIMING    (0)
#endif

/* Numeric Warnings */
/* #define PRINTNUMERICALWARNINGS */

/* maximum number of iterations  */
#define lpMPC_SET_MAXIT         (30)	

/* scaling factor of line search (affine direction) */
#define lpMPC_SET_LS_SCALE_AFF  (0.9)      

/* scaling factor of line search (combined direction) */
#define lpMPC_SET_LS_SCALE      (0.95)  

/* minimum required step size in each iteration */
#define lpMPC_SET_LS_MINSTEP    (1E-08)

/* maximum step size (combined direction) */
#define lpMPC_SET_LS_MAXSTEP    (0.995)

/* desired relative duality gap */
#define lpMPC_SET_ACC_RDGAP     (0.0001)

/* desired maximum residual on equality constraints */
#define lpMPC_SET_ACC_RESEQ     (1E-06)

/* desired maximum residual on inequality constraints */
#define lpMPC_SET_ACC_RESINEQ   (1E-06)

/* desired maximum violation of complementarity */
#define lpMPC_SET_ACC_KKTCOMPL  (1E-06)


/* RETURN CODES----------------------------------------------------------*/
/* solver has converged within desired accuracy */
#define lpMPC_OPTIMAL      (1)

/* maximum number of iterations has been reached */
#define lpMPC_MAXITREACHED (0)

/* no progress in line search possible */
#define lpMPC_NOPROGRESS   (-7)




/* PARAMETERS -----------------------------------------------------------*/
/* fill this with data before calling the solver! */
typedef struct lpMPC_params
{
    /* vector of size 4 */
    lpMPC_FLOAT f01[4];

    /* vector of size 4 */
    lpMPC_FLOAT lb01[4];

    /* vector of size 4 */
    lpMPC_FLOAT ub01[4];

    /* matrix of size [2 x 4] (column major format) */
    lpMPC_FLOAT C01[8];

    /* vector of size 2 */
    lpMPC_FLOAT e01[2];

    /* vector of size 4 */
    lpMPC_FLOAT f02[4];

    /* vector of size 4 */
    lpMPC_FLOAT lb02[4];

    /* vector of size 4 */
    lpMPC_FLOAT ub02[4];

    /* matrix of size [2 x 4] (column major format) */
    lpMPC_FLOAT C02[8];

    /* vector of size 2 */
    lpMPC_FLOAT e02[2];

    /* vector of size 4 */
    lpMPC_FLOAT f03[4];

    /* vector of size 4 */
    lpMPC_FLOAT lb03[4];

    /* vector of size 4 */
    lpMPC_FLOAT ub03[4];

    /* matrix of size [2 x 4] (column major format) */
    lpMPC_FLOAT C03[8];

    /* vector of size 2 */
    lpMPC_FLOAT e03[2];

    /* vector of size 4 */
    lpMPC_FLOAT f04[4];

    /* vector of size 4 */
    lpMPC_FLOAT lb04[4];

    /* vector of size 4 */
    lpMPC_FLOAT ub04[4];

    /* matrix of size [2 x 4] (column major format) */
    lpMPC_FLOAT C04[8];

    /* vector of size 2 */
    lpMPC_FLOAT e04[2];

    /* vector of size 4 */
    lpMPC_FLOAT f05[4];

    /* vector of size 4 */
    lpMPC_FLOAT lb05[4];

    /* vector of size 4 */
    lpMPC_FLOAT ub05[4];

    /* matrix of size [2 x 4] (column major format) */
    lpMPC_FLOAT C05[8];

    /* vector of size 2 */
    lpMPC_FLOAT e05[2];

    /* vector of size 4 */
    lpMPC_FLOAT f06[4];

    /* vector of size 4 */
    lpMPC_FLOAT lb06[4];

    /* vector of size 4 */
    lpMPC_FLOAT ub06[4];

    /* matrix of size [2 x 4] (column major format) */
    lpMPC_FLOAT C06[8];

    /* vector of size 2 */
    lpMPC_FLOAT e06[2];

    /* vector of size 4 */
    lpMPC_FLOAT f07[4];

    /* vector of size 4 */
    lpMPC_FLOAT lb07[4];

    /* vector of size 4 */
    lpMPC_FLOAT ub07[4];

    /* matrix of size [2 x 4] (column major format) */
    lpMPC_FLOAT C07[8];

    /* vector of size 2 */
    lpMPC_FLOAT e07[2];

    /* vector of size 4 */
    lpMPC_FLOAT f08[4];

    /* vector of size 4 */
    lpMPC_FLOAT lb08[4];

    /* vector of size 4 */
    lpMPC_FLOAT ub08[4];

    /* matrix of size [2 x 4] (column major format) */
    lpMPC_FLOAT C08[8];

    /* vector of size 2 */
    lpMPC_FLOAT e08[2];

    /* vector of size 4 */
    lpMPC_FLOAT f09[4];

    /* vector of size 4 */
    lpMPC_FLOAT lb09[4];

    /* vector of size 4 */
    lpMPC_FLOAT ub09[4];

    /* matrix of size [2 x 4] (column major format) */
    lpMPC_FLOAT C09[8];

    /* vector of size 2 */
    lpMPC_FLOAT e09[2];

    /* vector of size 4 */
    lpMPC_FLOAT f10[4];

    /* vector of size 4 */
    lpMPC_FLOAT lb10[4];

    /* vector of size 4 */
    lpMPC_FLOAT ub10[4];

    /* matrix of size [2 x 4] (column major format) */
    lpMPC_FLOAT C10[8];

    /* vector of size 2 */
    lpMPC_FLOAT e10[2];

    /* vector of size 4 */
    lpMPC_FLOAT f11[4];

    /* vector of size 4 */
    lpMPC_FLOAT lb11[4];

    /* vector of size 4 */
    lpMPC_FLOAT ub11[4];

    /* matrix of size [2 x 4] (column major format) */
    lpMPC_FLOAT C11[8];

    /* vector of size 2 */
    lpMPC_FLOAT e11[2];

    /* vector of size 4 */
    lpMPC_FLOAT f12[4];

    /* vector of size 4 */
    lpMPC_FLOAT lb12[4];

    /* vector of size 4 */
    lpMPC_FLOAT ub12[4];

    /* matrix of size [2 x 4] (column major format) */
    lpMPC_FLOAT C12[8];

    /* vector of size 2 */
    lpMPC_FLOAT e12[2];

    /* vector of size 4 */
    lpMPC_FLOAT f13[4];

    /* vector of size 4 */
    lpMPC_FLOAT lb13[4];

    /* vector of size 4 */
    lpMPC_FLOAT ub13[4];

    /* matrix of size [2 x 4] (column major format) */
    lpMPC_FLOAT C13[8];

    /* vector of size 2 */
    lpMPC_FLOAT e13[2];

    /* vector of size 4 */
    lpMPC_FLOAT f14[4];

    /* vector of size 4 */
    lpMPC_FLOAT lb14[4];

    /* vector of size 4 */
    lpMPC_FLOAT ub14[4];

    /* matrix of size [2 x 4] (column major format) */
    lpMPC_FLOAT C14[8];

    /* vector of size 2 */
    lpMPC_FLOAT e14[2];

    /* vector of size 2 */
    lpMPC_FLOAT f15[2];

    /* vector of size 2 */
    lpMPC_FLOAT lb15[2];

    /* vector of size 2 */
    lpMPC_FLOAT ub15[2];

    /* vector of size 2 */
    lpMPC_FLOAT e15[2];

} lpMPC_params;


/* OUTPUTS --------------------------------------------------------------*/
/* the desired variables are put here by the solver */
typedef struct lpMPC_output
{
    /* vector of size 4 */
    lpMPC_FLOAT z1[4];

    /* vector of size 4 */
    lpMPC_FLOAT z2[4];

    /* vector of size 4 */
    lpMPC_FLOAT z3[4];

    /* vector of size 4 */
    lpMPC_FLOAT z4[4];

    /* vector of size 4 */
    lpMPC_FLOAT z5[4];

    /* vector of size 4 */
    lpMPC_FLOAT z6[4];

    /* vector of size 4 */
    lpMPC_FLOAT z7[4];

    /* vector of size 4 */
    lpMPC_FLOAT z8[4];

    /* vector of size 4 */
    lpMPC_FLOAT z9[4];

    /* vector of size 4 */
    lpMPC_FLOAT z10[4];

    /* vector of size 4 */
    lpMPC_FLOAT z11[4];

    /* vector of size 4 */
    lpMPC_FLOAT z12[4];

    /* vector of size 4 */
    lpMPC_FLOAT z13[4];

    /* vector of size 4 */
    lpMPC_FLOAT z14[4];

    /* vector of size 2 */
    lpMPC_FLOAT z15[2];

} lpMPC_output;


/* SOLVER INFO ----------------------------------------------------------*/
/* diagnostic data from last interior point step */
typedef struct lpMPC_info
{
    /* iteration number */
    int it;
	
    /* inf-norm of equality constraint residuals */
    lpMPC_FLOAT res_eq;
	
    /* inf-norm of inequality constraint residuals */
    lpMPC_FLOAT res_ineq;

    /* primal objective */
    lpMPC_FLOAT pobj;	
	
    /* dual objective */
    lpMPC_FLOAT dobj;	

    /* duality gap := pobj - dobj */
    lpMPC_FLOAT dgap;		
	
    /* relative duality gap := |dgap / pobj | */
    lpMPC_FLOAT rdgap;		

    /* duality measure */
    lpMPC_FLOAT mu;

	/* duality measure (after affine step) */
    lpMPC_FLOAT mu_aff;
	
    /* centering parameter */
    lpMPC_FLOAT sigma;
	
    /* number of backtracking line search steps (affine direction) */
    int lsit_aff;
    
    /* number of backtracking line search steps (combined direction) */
    int lsit_cc;
    
    /* step size (affine direction) */
    lpMPC_FLOAT step_aff;
    
    /* step size (combined direction) */
    lpMPC_FLOAT step_cc;    

	/* solvertime */
	lpMPC_FLOAT solvetime;   

} lpMPC_info;


/* SOLVER FUNCTION DEFINITION -------------------------------------------*/
/* examine exitflag before using the result! */
int lpMPC_solve(lpMPC_params* params, lpMPC_output* output, lpMPC_info* info);


#endif