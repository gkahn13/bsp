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
    beliefMPC_FLOAT lb1[7];

    /* vector of size 7 */
    beliefMPC_FLOAT ub1[7];

    /* matrix of size [10 x 7] (column major format) */
    beliefMPC_FLOAT C1[70];

    /* vector of size 10 */
    beliefMPC_FLOAT e1[10];

    /* vector of size 7 */
    beliefMPC_FLOAT lb2[7];

    /* vector of size 7 */
    beliefMPC_FLOAT ub2[7];

    /* matrix of size [5 x 7] (column major format) */
    beliefMPC_FLOAT C2[35];

    /* vector of size 5 */
    beliefMPC_FLOAT e2[5];

    /* vector of size 7 */
    beliefMPC_FLOAT lb3[7];

    /* vector of size 7 */
    beliefMPC_FLOAT ub3[7];

    /* matrix of size [5 x 7] (column major format) */
    beliefMPC_FLOAT C3[35];

    /* vector of size 5 */
    beliefMPC_FLOAT e3[5];

    /* vector of size 7 */
    beliefMPC_FLOAT lb4[7];

    /* vector of size 7 */
    beliefMPC_FLOAT ub4[7];

    /* matrix of size [5 x 7] (column major format) */
    beliefMPC_FLOAT C4[35];

    /* vector of size 5 */
    beliefMPC_FLOAT e4[5];

    /* vector of size 7 */
    beliefMPC_FLOAT lb5[7];

    /* vector of size 7 */
    beliefMPC_FLOAT ub5[7];

    /* matrix of size [5 x 7] (column major format) */
    beliefMPC_FLOAT C5[35];

    /* vector of size 5 */
    beliefMPC_FLOAT e5[5];

    /* vector of size 7 */
    beliefMPC_FLOAT lb6[7];

    /* vector of size 7 */
    beliefMPC_FLOAT ub6[7];

    /* matrix of size [5 x 7] (column major format) */
    beliefMPC_FLOAT C6[35];

    /* vector of size 5 */
    beliefMPC_FLOAT e6[5];

    /* vector of size 7 */
    beliefMPC_FLOAT lb7[7];

    /* vector of size 7 */
    beliefMPC_FLOAT ub7[7];

    /* matrix of size [5 x 7] (column major format) */
    beliefMPC_FLOAT C7[35];

    /* vector of size 5 */
    beliefMPC_FLOAT e7[5];

    /* vector of size 7 */
    beliefMPC_FLOAT lb8[7];

    /* vector of size 7 */
    beliefMPC_FLOAT ub8[7];

    /* matrix of size [5 x 7] (column major format) */
    beliefMPC_FLOAT C8[35];

    /* vector of size 5 */
    beliefMPC_FLOAT e8[5];

    /* vector of size 7 */
    beliefMPC_FLOAT lb9[7];

    /* vector of size 7 */
    beliefMPC_FLOAT ub9[7];

    /* matrix of size [5 x 7] (column major format) */
    beliefMPC_FLOAT C9[35];

    /* vector of size 5 */
    beliefMPC_FLOAT e9[5];

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

    /* vector of size 7 */
    beliefMPC_FLOAT lb15[7];

    /* vector of size 7 */
    beliefMPC_FLOAT ub15[7];

    /* matrix of size [5 x 7] (column major format) */
    beliefMPC_FLOAT C15[35];

    /* vector of size 5 */
    beliefMPC_FLOAT e15[5];

    /* vector of size 7 */
    beliefMPC_FLOAT lb16[7];

    /* vector of size 7 */
    beliefMPC_FLOAT ub16[7];

    /* matrix of size [5 x 7] (column major format) */
    beliefMPC_FLOAT C16[35];

    /* vector of size 5 */
    beliefMPC_FLOAT e16[5];

    /* vector of size 7 */
    beliefMPC_FLOAT lb17[7];

    /* vector of size 7 */
    beliefMPC_FLOAT ub17[7];

    /* matrix of size [5 x 7] (column major format) */
    beliefMPC_FLOAT C17[35];

    /* vector of size 5 */
    beliefMPC_FLOAT e17[5];

    /* vector of size 7 */
    beliefMPC_FLOAT lb18[7];

    /* vector of size 7 */
    beliefMPC_FLOAT ub18[7];

    /* matrix of size [5 x 7] (column major format) */
    beliefMPC_FLOAT C18[35];

    /* vector of size 5 */
    beliefMPC_FLOAT e18[5];

    /* vector of size 7 */
    beliefMPC_FLOAT lb19[7];

    /* vector of size 7 */
    beliefMPC_FLOAT ub19[7];

    /* matrix of size [5 x 7] (column major format) */
    beliefMPC_FLOAT C19[35];

    /* vector of size 5 */
    beliefMPC_FLOAT e19[5];

    /* vector of size 5 */
    beliefMPC_FLOAT lb20[5];

    /* vector of size 5 */
    beliefMPC_FLOAT ub20[5];

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

    /* vector of size 7 */
    beliefMPC_FLOAT z15[7];

    /* vector of size 7 */
    beliefMPC_FLOAT z16[7];

    /* vector of size 7 */
    beliefMPC_FLOAT z17[7];

    /* vector of size 7 */
    beliefMPC_FLOAT z18[7];

    /* vector of size 7 */
    beliefMPC_FLOAT z19[7];

    /* vector of size 5 */
    beliefMPC_FLOAT z20[5];

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