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

#ifndef __beliefMPC_H__
#define __beliefMPC_H__


/* DATA TYPE ------------------------------------------------------------*/
typedef double beliefMPC_FLOAT;


/* SOLVER SETTINGS ------------------------------------------------------*/
/* print level */
#ifndef beliefMPC_SET_PRINTLEVEL
#define beliefMPC_SET_PRINTLEVEL    (0)
#endif

/* timing */
#ifndef beliefMPC_SET_TIMING
#define beliefMPC_SET_TIMING    (0)
#endif

/* Numeric Warnings */
/* #define PRINTNUMERICALWARNINGS */

/* maximum number of iterations  */
#define beliefMPC_SET_MAXIT         (30)	

/* scaling factor of line search (affine direction) */
#define beliefMPC_SET_LS_SCALE_AFF  (0.9)      

/* scaling factor of line search (combined direction) */
#define beliefMPC_SET_LS_SCALE      (0.95)  

/* minimum required step size in each iteration */
#define beliefMPC_SET_LS_MINSTEP    (1E-08)

/* maximum step size (combined direction) */
#define beliefMPC_SET_LS_MAXSTEP    (0.995)

/* desired relative duality gap */
#define beliefMPC_SET_ACC_RDGAP     (0.0001)

/* desired maximum residual on equality constraints */
#define beliefMPC_SET_ACC_RESEQ     (1E-06)

/* desired maximum residual on inequality constraints */
#define beliefMPC_SET_ACC_RESINEQ   (1E-06)

/* desired maximum violation of complementarity */
#define beliefMPC_SET_ACC_KKTCOMPL  (1E-06)


/* RETURN CODES----------------------------------------------------------*/
/* solver has converged within desired accuracy */
#define beliefMPC_OPTIMAL      (1)

/* maximum number of iterations has been reached */
#define beliefMPC_MAXITREACHED (0)

/* no progress in line search possible */
#define beliefMPC_NOPROGRESS   (-7)




/* PARAMETERS -----------------------------------------------------------*/
/* fill this with data before calling the solver! */
typedef struct beliefMPC_params
{
    /* vector of size 7 */
    beliefMPC_FLOAT lb01[7];

    /* vector of size 7 */
    beliefMPC_FLOAT ub01[7];

    /* matrix of size [10 x 7] (column major format) */
    beliefMPC_FLOAT C01[70];

    /* vector of size 10 */
    beliefMPC_FLOAT e01[10];

    /* vector of size 7 */
    beliefMPC_FLOAT lb02[7];

    /* vector of size 7 */
    beliefMPC_FLOAT ub02[7];

    /* matrix of size [5 x 7] (column major format) */
    beliefMPC_FLOAT C02[35];

    /* vector of size 5 */
    beliefMPC_FLOAT e02[5];

    /* vector of size 7 */
    beliefMPC_FLOAT lb03[7];

    /* vector of size 7 */
    beliefMPC_FLOAT ub03[7];

    /* matrix of size [5 x 7] (column major format) */
    beliefMPC_FLOAT C03[35];

    /* vector of size 5 */
    beliefMPC_FLOAT e03[5];

    /* vector of size 7 */
    beliefMPC_FLOAT lb04[7];

    /* vector of size 7 */
    beliefMPC_FLOAT ub04[7];

    /* matrix of size [5 x 7] (column major format) */
    beliefMPC_FLOAT C04[35];

    /* vector of size 5 */
    beliefMPC_FLOAT e04[5];

    /* vector of size 7 */
    beliefMPC_FLOAT lb05[7];

    /* vector of size 7 */
    beliefMPC_FLOAT ub05[7];

    /* matrix of size [5 x 7] (column major format) */
    beliefMPC_FLOAT C05[35];

    /* vector of size 5 */
    beliefMPC_FLOAT e05[5];

    /* vector of size 7 */
    beliefMPC_FLOAT lb06[7];

    /* vector of size 7 */
    beliefMPC_FLOAT ub06[7];

    /* matrix of size [5 x 7] (column major format) */
    beliefMPC_FLOAT C06[35];

    /* vector of size 5 */
    beliefMPC_FLOAT e06[5];

    /* vector of size 7 */
    beliefMPC_FLOAT lb07[7];

    /* vector of size 7 */
    beliefMPC_FLOAT ub07[7];

    /* matrix of size [5 x 7] (column major format) */
    beliefMPC_FLOAT C07[35];

    /* vector of size 5 */
    beliefMPC_FLOAT e07[5];

    /* vector of size 7 */
    beliefMPC_FLOAT lb08[7];

    /* vector of size 7 */
    beliefMPC_FLOAT ub08[7];

    /* matrix of size [5 x 7] (column major format) */
    beliefMPC_FLOAT C08[35];

    /* vector of size 5 */
    beliefMPC_FLOAT e08[5];

    /* vector of size 7 */
    beliefMPC_FLOAT lb09[7];

    /* vector of size 7 */
    beliefMPC_FLOAT ub09[7];

    /* matrix of size [5 x 7] (column major format) */
    beliefMPC_FLOAT C09[35];

    /* vector of size 5 */
    beliefMPC_FLOAT e09[5];

    /* vector of size 7 */
    beliefMPC_FLOAT lb10[7];

    /* vector of size 7 */
    beliefMPC_FLOAT ub10[7];

    /* matrix of size [5 x 7] (column major format) */
    beliefMPC_FLOAT C10[35];

    /* vector of size 5 */
    beliefMPC_FLOAT e10[5];

    /* vector of size 7 */
    beliefMPC_FLOAT lb11[7];

    /* vector of size 7 */
    beliefMPC_FLOAT ub11[7];

    /* matrix of size [5 x 7] (column major format) */
    beliefMPC_FLOAT C11[35];

    /* vector of size 5 */
    beliefMPC_FLOAT e11[5];

    /* vector of size 7 */
    beliefMPC_FLOAT lb12[7];

    /* vector of size 7 */
    beliefMPC_FLOAT ub12[7];

    /* matrix of size [5 x 7] (column major format) */
    beliefMPC_FLOAT C12[35];

    /* vector of size 5 */
    beliefMPC_FLOAT e12[5];

    /* vector of size 7 */
    beliefMPC_FLOAT lb13[7];

    /* vector of size 7 */
    beliefMPC_FLOAT ub13[7];

    /* matrix of size [5 x 7] (column major format) */
    beliefMPC_FLOAT C13[35];

    /* vector of size 5 */
    beliefMPC_FLOAT e13[5];

    /* vector of size 7 */
    beliefMPC_FLOAT lb14[7];

    /* vector of size 7 */
    beliefMPC_FLOAT ub14[7];

    /* matrix of size [5 x 7] (column major format) */
    beliefMPC_FLOAT C14[35];

    /* vector of size 5 */
    beliefMPC_FLOAT e14[5];

    /* vector of size 5 */
    beliefMPC_FLOAT lb15[5];

    /* vector of size 5 */
    beliefMPC_FLOAT ub15[5];

} beliefMPC_params;


/* OUTPUTS --------------------------------------------------------------*/
/* the desired variables are put here by the solver */
typedef struct beliefMPC_output
{
    /* vector of size 7 */
    beliefMPC_FLOAT z1[7];

    /* vector of size 7 */
    beliefMPC_FLOAT z2[7];

    /* vector of size 7 */
    beliefMPC_FLOAT z3[7];

    /* vector of size 7 */
    beliefMPC_FLOAT z4[7];

    /* vector of size 7 */
    beliefMPC_FLOAT z5[7];

    /* vector of size 7 */
    beliefMPC_FLOAT z6[7];

    /* vector of size 7 */
    beliefMPC_FLOAT z7[7];

    /* vector of size 7 */
    beliefMPC_FLOAT z8[7];

    /* vector of size 7 */
    beliefMPC_FLOAT z9[7];

    /* vector of size 7 */
    beliefMPC_FLOAT z10[7];

    /* vector of size 7 */
    beliefMPC_FLOAT z11[7];

    /* vector of size 7 */
    beliefMPC_FLOAT z12[7];

    /* vector of size 7 */
    beliefMPC_FLOAT z13[7];

    /* vector of size 7 */
    beliefMPC_FLOAT z14[7];

    /* vector of size 5 */
    beliefMPC_FLOAT z15[5];

} beliefMPC_output;


/* SOLVER INFO ----------------------------------------------------------*/
/* diagnostic data from last interior point step */
typedef struct beliefMPC_info
{
    /* iteration number */
    int it;
	
    /* inf-norm of equality constraint residuals */
    beliefMPC_FLOAT res_eq;
	
    /* inf-norm of inequality constraint residuals */
    beliefMPC_FLOAT res_ineq;

    /* primal objective */
    beliefMPC_FLOAT pobj;	
	
    /* dual objective */
    beliefMPC_FLOAT dobj;	

    /* duality gap := pobj - dobj */
    beliefMPC_FLOAT dgap;		
	
    /* relative duality gap := |dgap / pobj | */
    beliefMPC_FLOAT rdgap;		

    /* duality measure */
    beliefMPC_FLOAT mu;

	/* duality measure (after affine step) */
    beliefMPC_FLOAT mu_aff;
	
    /* centering parameter */
    beliefMPC_FLOAT sigma;
	
    /* number of backtracking line search steps (affine direction) */
    int lsit_aff;
    
    /* number of backtracking line search steps (combined direction) */
    int lsit_cc;
    
    /* step size (affine direction) */
    beliefMPC_FLOAT step_aff;
    
    /* step size (combined direction) */
    beliefMPC_FLOAT step_cc;    

	/* solvertime */
	beliefMPC_FLOAT solvetime;   

} beliefMPC_info;


/* SOLVER FUNCTION DEFINITION -------------------------------------------*/
/* examine exitflag before using the result! */
int beliefMPC_solve(beliefMPC_params* params, beliefMPC_output* output, beliefMPC_info* info);


#endif