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

#ifndef __controlMPC_H__
#define __controlMPC_H__


/* DATA TYPE ------------------------------------------------------------*/
typedef double controlMPC_FLOAT;


/* SOLVER SETTINGS ------------------------------------------------------*/
/* print level */
#ifndef controlMPC_SET_PRINTLEVEL
#define controlMPC_SET_PRINTLEVEL    (0)
#endif

/* timing */
#ifndef controlMPC_SET_TIMING
#define controlMPC_SET_TIMING    (0)
#endif

/* Numeric Warnings */
/* #define PRINTNUMERICALWARNINGS */

/* maximum number of iterations  */
#define controlMPC_SET_MAXIT         (30)	

/* scaling factor of line search (affine direction) */
#define controlMPC_SET_LS_SCALE_AFF  (0.9)      

/* scaling factor of line search (combined direction) */
#define controlMPC_SET_LS_SCALE      (0.95)  

/* minimum required step size in each iteration */
#define controlMPC_SET_LS_MINSTEP    (1E-08)

/* maximum step size (combined direction) */
#define controlMPC_SET_LS_MAXSTEP    (0.995)

/* desired relative duality gap */
#define controlMPC_SET_ACC_RDGAP     (0.0001)

/* desired maximum residual on equality constraints */
#define controlMPC_SET_ACC_RESEQ     (1E-06)

/* desired maximum residual on inequality constraints */
#define controlMPC_SET_ACC_RESINEQ   (1E-06)

/* desired maximum violation of complementarity */
#define controlMPC_SET_ACC_KKTCOMPL  (1E-06)


/* RETURN CODES----------------------------------------------------------*/
/* solver has converged within desired accuracy */
#define controlMPC_OPTIMAL      (1)

/* maximum number of iterations has been reached */
#define controlMPC_MAXITREACHED (0)

/* no progress in line search possible */
#define controlMPC_NOPROGRESS   (-7)




/* PARAMETERS -----------------------------------------------------------*/
/* fill this with data before calling the solver! */
typedef struct controlMPC_params
{
    /* diagonal matrix of size [2 x 2] (only the diagonal is stored) */
    controlMPC_FLOAT H1[2];

    /* vector of size 2 */
    controlMPC_FLOAT f1[2];

    /* vector of size 2 */
    controlMPC_FLOAT lb1[2];

    /* vector of size 2 */
    controlMPC_FLOAT ub1[2];

    /* diagonal matrix of size [2 x 2] (only the diagonal is stored) */
    controlMPC_FLOAT H2[2];

    /* vector of size 2 */
    controlMPC_FLOAT f2[2];

    /* vector of size 2 */
    controlMPC_FLOAT lb2[2];

    /* vector of size 2 */
    controlMPC_FLOAT ub2[2];

    /* diagonal matrix of size [2 x 2] (only the diagonal is stored) */
    controlMPC_FLOAT H3[2];

    /* vector of size 2 */
    controlMPC_FLOAT f3[2];

    /* vector of size 2 */
    controlMPC_FLOAT lb3[2];

    /* vector of size 2 */
    controlMPC_FLOAT ub3[2];

    /* diagonal matrix of size [2 x 2] (only the diagonal is stored) */
    controlMPC_FLOAT H4[2];

    /* vector of size 2 */
    controlMPC_FLOAT f4[2];

    /* vector of size 2 */
    controlMPC_FLOAT lb4[2];

    /* vector of size 2 */
    controlMPC_FLOAT ub4[2];

    /* diagonal matrix of size [2 x 2] (only the diagonal is stored) */
    controlMPC_FLOAT H5[2];

    /* vector of size 2 */
    controlMPC_FLOAT f5[2];

    /* vector of size 2 */
    controlMPC_FLOAT lb5[2];

    /* vector of size 2 */
    controlMPC_FLOAT ub5[2];

    /* diagonal matrix of size [2 x 2] (only the diagonal is stored) */
    controlMPC_FLOAT H6[2];

    /* vector of size 2 */
    controlMPC_FLOAT f6[2];

    /* vector of size 2 */
    controlMPC_FLOAT lb6[2];

    /* vector of size 2 */
    controlMPC_FLOAT ub6[2];

    /* diagonal matrix of size [2 x 2] (only the diagonal is stored) */
    controlMPC_FLOAT H7[2];

    /* vector of size 2 */
    controlMPC_FLOAT f7[2];

    /* vector of size 2 */
    controlMPC_FLOAT lb7[2];

    /* vector of size 2 */
    controlMPC_FLOAT ub7[2];

    /* diagonal matrix of size [2 x 2] (only the diagonal is stored) */
    controlMPC_FLOAT H8[2];

    /* vector of size 2 */
    controlMPC_FLOAT f8[2];

    /* vector of size 2 */
    controlMPC_FLOAT lb8[2];

    /* vector of size 2 */
    controlMPC_FLOAT ub8[2];

    /* diagonal matrix of size [2 x 2] (only the diagonal is stored) */
    controlMPC_FLOAT H9[2];

    /* vector of size 2 */
    controlMPC_FLOAT f9[2];

    /* vector of size 2 */
    controlMPC_FLOAT lb9[2];

    /* vector of size 2 */
    controlMPC_FLOAT ub9[2];

    /* diagonal matrix of size [2 x 2] (only the diagonal is stored) */
    controlMPC_FLOAT H10[2];

    /* vector of size 2 */
    controlMPC_FLOAT f10[2];

    /* vector of size 2 */
    controlMPC_FLOAT lb10[2];

    /* vector of size 2 */
    controlMPC_FLOAT ub10[2];

    /* diagonal matrix of size [2 x 2] (only the diagonal is stored) */
    controlMPC_FLOAT H11[2];

    /* vector of size 2 */
    controlMPC_FLOAT f11[2];

    /* vector of size 2 */
    controlMPC_FLOAT lb11[2];

    /* vector of size 2 */
    controlMPC_FLOAT ub11[2];

    /* diagonal matrix of size [2 x 2] (only the diagonal is stored) */
    controlMPC_FLOAT H12[2];

    /* vector of size 2 */
    controlMPC_FLOAT f12[2];

    /* vector of size 2 */
    controlMPC_FLOAT lb12[2];

    /* vector of size 2 */
    controlMPC_FLOAT ub12[2];

    /* diagonal matrix of size [2 x 2] (only the diagonal is stored) */
    controlMPC_FLOAT H13[2];

    /* vector of size 2 */
    controlMPC_FLOAT f13[2];

    /* vector of size 2 */
    controlMPC_FLOAT lb13[2];

    /* vector of size 2 */
    controlMPC_FLOAT ub13[2];

    /* diagonal matrix of size [2 x 2] (only the diagonal is stored) */
    controlMPC_FLOAT H14[2];

    /* vector of size 2 */
    controlMPC_FLOAT f14[2];

    /* vector of size 2 */
    controlMPC_FLOAT lb14[2];

    /* vector of size 2 */
    controlMPC_FLOAT ub14[2];

    /* diagonal matrix of size [2 x 2] (only the diagonal is stored) */
    controlMPC_FLOAT H15[2];

    /* vector of size 2 */
    controlMPC_FLOAT f15[2];

    /* vector of size 2 */
    controlMPC_FLOAT lb15[2];

    /* vector of size 2 */
    controlMPC_FLOAT ub15[2];

    /* diagonal matrix of size [2 x 2] (only the diagonal is stored) */
    controlMPC_FLOAT H16[2];

    /* vector of size 2 */
    controlMPC_FLOAT f16[2];

    /* vector of size 2 */
    controlMPC_FLOAT lb16[2];

    /* vector of size 2 */
    controlMPC_FLOAT ub16[2];

    /* diagonal matrix of size [2 x 2] (only the diagonal is stored) */
    controlMPC_FLOAT H17[2];

    /* vector of size 2 */
    controlMPC_FLOAT f17[2];

    /* vector of size 2 */
    controlMPC_FLOAT lb17[2];

    /* vector of size 2 */
    controlMPC_FLOAT ub17[2];

    /* diagonal matrix of size [2 x 2] (only the diagonal is stored) */
    controlMPC_FLOAT H18[2];

    /* vector of size 2 */
    controlMPC_FLOAT f18[2];

    /* vector of size 2 */
    controlMPC_FLOAT lb18[2];

    /* vector of size 2 */
    controlMPC_FLOAT ub18[2];

    /* diagonal matrix of size [2 x 2] (only the diagonal is stored) */
    controlMPC_FLOAT H19[2];

    /* vector of size 2 */
    controlMPC_FLOAT f19[2];

    /* vector of size 2 */
    controlMPC_FLOAT lb19[2];

    /* vector of size 2 */
    controlMPC_FLOAT ub19[2];

    /* diagonal matrix of size [2 x 2] (only the diagonal is stored) */
    controlMPC_FLOAT H20[2];

    /* vector of size 2 */
    controlMPC_FLOAT f20[2];

    /* vector of size 2 */
    controlMPC_FLOAT lb20[2];

    /* vector of size 2 */
    controlMPC_FLOAT ub20[2];

    /* diagonal matrix of size [2 x 2] (only the diagonal is stored) */
    controlMPC_FLOAT H21[2];

    /* vector of size 2 */
    controlMPC_FLOAT f21[2];

    /* vector of size 2 */
    controlMPC_FLOAT lb21[2];

    /* vector of size 2 */
    controlMPC_FLOAT ub21[2];

    /* diagonal matrix of size [2 x 2] (only the diagonal is stored) */
    controlMPC_FLOAT H22[2];

    /* vector of size 2 */
    controlMPC_FLOAT f22[2];

    /* vector of size 2 */
    controlMPC_FLOAT lb22[2];

    /* vector of size 2 */
    controlMPC_FLOAT ub22[2];

    /* diagonal matrix of size [2 x 2] (only the diagonal is stored) */
    controlMPC_FLOAT H23[2];

    /* vector of size 2 */
    controlMPC_FLOAT f23[2];

    /* vector of size 2 */
    controlMPC_FLOAT lb23[2];

    /* vector of size 2 */
    controlMPC_FLOAT ub23[2];

    /* diagonal matrix of size [2 x 2] (only the diagonal is stored) */
    controlMPC_FLOAT H24[2];

    /* vector of size 2 */
    controlMPC_FLOAT f24[2];

    /* vector of size 2 */
    controlMPC_FLOAT lb24[2];

    /* vector of size 2 */
    controlMPC_FLOAT ub24[2];

} controlMPC_params;


/* OUTPUTS --------------------------------------------------------------*/
/* the desired variables are put here by the solver */
typedef struct controlMPC_output
{
    /* vector of size 2 */
    controlMPC_FLOAT z1[2];

    /* vector of size 2 */
    controlMPC_FLOAT z2[2];

    /* vector of size 2 */
    controlMPC_FLOAT z3[2];

    /* vector of size 2 */
    controlMPC_FLOAT z4[2];

    /* vector of size 2 */
    controlMPC_FLOAT z5[2];

    /* vector of size 2 */
    controlMPC_FLOAT z6[2];

    /* vector of size 2 */
    controlMPC_FLOAT z7[2];

    /* vector of size 2 */
    controlMPC_FLOAT z8[2];

    /* vector of size 2 */
    controlMPC_FLOAT z9[2];

    /* vector of size 2 */
    controlMPC_FLOAT z10[2];

    /* vector of size 2 */
    controlMPC_FLOAT z11[2];

    /* vector of size 2 */
    controlMPC_FLOAT z12[2];

    /* vector of size 2 */
    controlMPC_FLOAT z13[2];

    /* vector of size 2 */
    controlMPC_FLOAT z14[2];

    /* vector of size 2 */
    controlMPC_FLOAT z15[2];

    /* vector of size 2 */
    controlMPC_FLOAT z16[2];

    /* vector of size 2 */
    controlMPC_FLOAT z17[2];

    /* vector of size 2 */
    controlMPC_FLOAT z18[2];

    /* vector of size 2 */
    controlMPC_FLOAT z19[2];

    /* vector of size 2 */
    controlMPC_FLOAT z20[2];

    /* vector of size 2 */
    controlMPC_FLOAT z21[2];

    /* vector of size 2 */
    controlMPC_FLOAT z22[2];

    /* vector of size 2 */
    controlMPC_FLOAT z23[2];

    /* vector of size 2 */
    controlMPC_FLOAT z24[2];

} controlMPC_output;


/* SOLVER INFO ----------------------------------------------------------*/
/* diagnostic data from last interior point step */
typedef struct controlMPC_info
{
    /* iteration number */
    int it;
	
    /* inf-norm of equality constraint residuals */
    controlMPC_FLOAT res_eq;
	
    /* inf-norm of inequality constraint residuals */
    controlMPC_FLOAT res_ineq;

    /* primal objective */
    controlMPC_FLOAT pobj;	
	
    /* dual objective */
    controlMPC_FLOAT dobj;	

    /* duality gap := pobj - dobj */
    controlMPC_FLOAT dgap;		
	
    /* relative duality gap := |dgap / pobj | */
    controlMPC_FLOAT rdgap;		

    /* duality measure */
    controlMPC_FLOAT mu;

	/* duality measure (after affine step) */
    controlMPC_FLOAT mu_aff;
	
    /* centering parameter */
    controlMPC_FLOAT sigma;
	
    /* number of backtracking line search steps (affine direction) */
    int lsit_aff;
    
    /* number of backtracking line search steps (combined direction) */
    int lsit_cc;
    
    /* step size (affine direction) */
    controlMPC_FLOAT step_aff;
    
    /* step size (combined direction) */
    controlMPC_FLOAT step_cc;    

	/* solvertime */
	controlMPC_FLOAT solvetime;   

} controlMPC_info;


/* SOLVER FUNCTION DEFINITION -------------------------------------------*/
/* examine exitflag before using the result! */
int controlMPC_solve(controlMPC_params* params, controlMPC_output* output, controlMPC_info* info);


#endif