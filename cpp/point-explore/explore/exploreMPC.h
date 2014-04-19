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

#ifndef __exploreMPC_H__
#define __exploreMPC_H__


/* DATA TYPE ------------------------------------------------------------*/
typedef double exploreMPC_FLOAT;


/* SOLVER SETTINGS ------------------------------------------------------*/
/* print level */
#ifndef exploreMPC_SET_PRINTLEVEL
#define exploreMPC_SET_PRINTLEVEL    (0)
#endif

/* timing */
#ifndef exploreMPC_SET_TIMING
#define exploreMPC_SET_TIMING    (0)
#endif

/* Numeric Warnings */
/* #define PRINTNUMERICALWARNINGS */

/* maximum number of iterations  */
#define exploreMPC_SET_MAXIT         (50)	

/* scaling factor of line search (affine direction) */
#define exploreMPC_SET_LS_SCALE_AFF  (0.9)      

/* scaling factor of line search (combined direction) */
#define exploreMPC_SET_LS_SCALE      (0.95)  

/* minimum required step size in each iteration */
#define exploreMPC_SET_LS_MINSTEP    (1E-08)

/* maximum step size (combined direction) */
#define exploreMPC_SET_LS_MAXSTEP    (0.995)

/* desired relative duality gap */
#define exploreMPC_SET_ACC_RDGAP     (0.0001)

/* desired maximum residual on equality constraints */
#define exploreMPC_SET_ACC_RESEQ     (1E-06)

/* desired maximum residual on inequality constraints */
#define exploreMPC_SET_ACC_RESINEQ   (1E-06)

/* desired maximum violation of complementarity */
#define exploreMPC_SET_ACC_KKTCOMPL  (1E-06)


/* RETURN CODES----------------------------------------------------------*/
/* solver has converged within desired accuracy */
#define exploreMPC_OPTIMAL      (1)

/* maximum number of iterations has been reached */
#define exploreMPC_MAXITREACHED (0)

/* no progress in line search possible */
#define exploreMPC_NOPROGRESS   (-7)




/* PARAMETERS -----------------------------------------------------------*/
/* fill this with data before calling the solver! */
typedef struct exploreMPC_params
{
    /* diagonal matrix of size [4 x 4] (only the diagonal is stored) */
    exploreMPC_FLOAT H1[4];

    /* vector of size 4 */
    exploreMPC_FLOAT f1[4];

    /* vector of size 4 */
    exploreMPC_FLOAT lb1[4];

    /* vector of size 4 */
    exploreMPC_FLOAT ub1[4];

    /* vector of size 2 */
    exploreMPC_FLOAT c1[2];

    /* diagonal matrix of size [4 x 4] (only the diagonal is stored) */
    exploreMPC_FLOAT H2[4];

    /* vector of size 4 */
    exploreMPC_FLOAT f2[4];

    /* vector of size 4 */
    exploreMPC_FLOAT lb2[4];

    /* vector of size 4 */
    exploreMPC_FLOAT ub2[4];

    /* diagonal matrix of size [4 x 4] (only the diagonal is stored) */
    exploreMPC_FLOAT H3[4];

    /* vector of size 4 */
    exploreMPC_FLOAT f3[4];

    /* vector of size 4 */
    exploreMPC_FLOAT lb3[4];

    /* vector of size 4 */
    exploreMPC_FLOAT ub3[4];

    /* diagonal matrix of size [4 x 4] (only the diagonal is stored) */
    exploreMPC_FLOAT H4[4];

    /* vector of size 4 */
    exploreMPC_FLOAT f4[4];

    /* vector of size 4 */
    exploreMPC_FLOAT lb4[4];

    /* vector of size 4 */
    exploreMPC_FLOAT ub4[4];

    /* diagonal matrix of size [4 x 4] (only the diagonal is stored) */
    exploreMPC_FLOAT H5[4];

    /* vector of size 4 */
    exploreMPC_FLOAT f5[4];

    /* vector of size 4 */
    exploreMPC_FLOAT lb5[4];

    /* vector of size 4 */
    exploreMPC_FLOAT ub5[4];

    /* diagonal matrix of size [4 x 4] (only the diagonal is stored) */
    exploreMPC_FLOAT H6[4];

    /* vector of size 4 */
    exploreMPC_FLOAT f6[4];

    /* vector of size 4 */
    exploreMPC_FLOAT lb6[4];

    /* vector of size 4 */
    exploreMPC_FLOAT ub6[4];

    /* diagonal matrix of size [4 x 4] (only the diagonal is stored) */
    exploreMPC_FLOAT H7[4];

    /* vector of size 4 */
    exploreMPC_FLOAT f7[4];

    /* vector of size 4 */
    exploreMPC_FLOAT lb7[4];

    /* vector of size 4 */
    exploreMPC_FLOAT ub7[4];

    /* diagonal matrix of size [4 x 4] (only the diagonal is stored) */
    exploreMPC_FLOAT H8[4];

    /* vector of size 4 */
    exploreMPC_FLOAT f8[4];

    /* vector of size 4 */
    exploreMPC_FLOAT lb8[4];

    /* vector of size 4 */
    exploreMPC_FLOAT ub8[4];

    /* diagonal matrix of size [4 x 4] (only the diagonal is stored) */
    exploreMPC_FLOAT H9[4];

    /* vector of size 4 */
    exploreMPC_FLOAT f9[4];

    /* vector of size 4 */
    exploreMPC_FLOAT lb9[4];

    /* vector of size 4 */
    exploreMPC_FLOAT ub9[4];

    /* diagonal matrix of size [2 x 2] (only the diagonal is stored) */
    exploreMPC_FLOAT H10[2];

    /* vector of size 2 */
    exploreMPC_FLOAT f10[2];

    /* vector of size 2 */
    exploreMPC_FLOAT lb10[2];

    /* vector of size 2 */
    exploreMPC_FLOAT ub10[2];

} exploreMPC_params;


/* OUTPUTS --------------------------------------------------------------*/
/* the desired variables are put here by the solver */
typedef struct exploreMPC_output
{
    /* vector of size 4 */
    exploreMPC_FLOAT z1[4];

    /* vector of size 4 */
    exploreMPC_FLOAT z2[4];

    /* vector of size 4 */
    exploreMPC_FLOAT z3[4];

    /* vector of size 4 */
    exploreMPC_FLOAT z4[4];

    /* vector of size 4 */
    exploreMPC_FLOAT z5[4];

    /* vector of size 4 */
    exploreMPC_FLOAT z6[4];

    /* vector of size 4 */
    exploreMPC_FLOAT z7[4];

    /* vector of size 4 */
    exploreMPC_FLOAT z8[4];

    /* vector of size 4 */
    exploreMPC_FLOAT z9[4];

    /* vector of size 2 */
    exploreMPC_FLOAT z10[2];

} exploreMPC_output;


/* SOLVER INFO ----------------------------------------------------------*/
/* diagnostic data from last interior point step */
typedef struct exploreMPC_info
{
    /* iteration number */
    int it;
	
    /* inf-norm of equality constraint residuals */
    exploreMPC_FLOAT res_eq;
	
    /* inf-norm of inequality constraint residuals */
    exploreMPC_FLOAT res_ineq;

    /* primal objective */
    exploreMPC_FLOAT pobj;	
	
    /* dual objective */
    exploreMPC_FLOAT dobj;	

    /* duality gap := pobj - dobj */
    exploreMPC_FLOAT dgap;		
	
    /* relative duality gap := |dgap / pobj | */
    exploreMPC_FLOAT rdgap;		

    /* duality measure */
    exploreMPC_FLOAT mu;

	/* duality measure (after affine step) */
    exploreMPC_FLOAT mu_aff;
	
    /* centering parameter */
    exploreMPC_FLOAT sigma;
	
    /* number of backtracking line search steps (affine direction) */
    int lsit_aff;
    
    /* number of backtracking line search steps (combined direction) */
    int lsit_cc;
    
    /* step size (affine direction) */
    exploreMPC_FLOAT step_aff;
    
    /* step size (combined direction) */
    exploreMPC_FLOAT step_cc;    

	/* solvertime */
	exploreMPC_FLOAT solvetime;   

} exploreMPC_info;


/* SOLVER FUNCTION DEFINITION -------------------------------------------*/
/* examine exitflag before using the result! */
int exploreMPC_solve(exploreMPC_params* params, exploreMPC_output* output, exploreMPC_info* info);


#endif