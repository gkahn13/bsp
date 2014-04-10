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

#ifndef __plattMPC_H__
#define __plattMPC_H__


/* DATA TYPE ------------------------------------------------------------*/
typedef double plattMPC_FLOAT;


/* SOLVER SETTINGS ------------------------------------------------------*/
/* print level */
#ifndef plattMPC_SET_PRINTLEVEL
#define plattMPC_SET_PRINTLEVEL    (2)
#endif

/* timing */
#ifndef plattMPC_SET_TIMING
#define plattMPC_SET_TIMING    (0)
#endif

/* Numeric Warnings */
/* #define PRINTNUMERICALWARNINGS */

/* maximum number of iterations  */
#define plattMPC_SET_MAXIT         (100)	

/* scaling factor of line search (affine direction) */
#define plattMPC_SET_LS_SCALE_AFF  (0.9)      

/* scaling factor of line search (combined direction) */
#define plattMPC_SET_LS_SCALE      (0.95)  

/* minimum required step size in each iteration */
#define plattMPC_SET_LS_MINSTEP    (1E-08)

/* maximum step size (combined direction) */
#define plattMPC_SET_LS_MAXSTEP    (0.995)

/* desired relative duality gap */
#define plattMPC_SET_ACC_RDGAP     (0.0001)

/* desired maximum residual on equality constraints */
#define plattMPC_SET_ACC_RESEQ     (1E-06)

/* desired maximum residual on inequality constraints */
#define plattMPC_SET_ACC_RESINEQ   (1E-06)

/* desired maximum violation of complementarity */
#define plattMPC_SET_ACC_KKTCOMPL  (1E-06)


/* RETURN CODES----------------------------------------------------------*/
/* solver has converged within desired accuracy */
#define plattMPC_OPTIMAL      (1)

/* maximum number of iterations has been reached */
#define plattMPC_MAXITREACHED (0)

/* no progress in line search possible */
#define plattMPC_NOPROGRESS   (-7)




/* PARAMETERS -----------------------------------------------------------*/
/* fill this with data before calling the solver! */
typedef struct plattMPC_params
{
    /* diagonal matrix of size [12 x 12] (only the diagonal is stored) */
    plattMPC_FLOAT H1[12];

    /* vector of size 12 */
    plattMPC_FLOAT f1[12];

    /* vector of size 12 */
    plattMPC_FLOAT lb1[12];

    /* vector of size 12 */
    plattMPC_FLOAT ub1[12];

    /* vector of size 10 */
    plattMPC_FLOAT c1[10];

    /* diagonal matrix of size [12 x 12] (only the diagonal is stored) */
    plattMPC_FLOAT H2[12];

    /* vector of size 12 */
    plattMPC_FLOAT f2[12];

    /* vector of size 12 */
    plattMPC_FLOAT lb2[12];

    /* vector of size 12 */
    plattMPC_FLOAT ub2[12];

    /* diagonal matrix of size [12 x 12] (only the diagonal is stored) */
    plattMPC_FLOAT H3[12];

    /* vector of size 12 */
    plattMPC_FLOAT f3[12];

    /* vector of size 12 */
    plattMPC_FLOAT lb3[12];

    /* vector of size 12 */
    plattMPC_FLOAT ub3[12];

    /* diagonal matrix of size [12 x 12] (only the diagonal is stored) */
    plattMPC_FLOAT H4[12];

    /* vector of size 12 */
    plattMPC_FLOAT f4[12];

    /* vector of size 12 */
    plattMPC_FLOAT lb4[12];

    /* vector of size 12 */
    plattMPC_FLOAT ub4[12];

    /* diagonal matrix of size [12 x 12] (only the diagonal is stored) */
    plattMPC_FLOAT H5[12];

    /* vector of size 12 */
    plattMPC_FLOAT f5[12];

    /* vector of size 12 */
    plattMPC_FLOAT lb5[12];

    /* vector of size 12 */
    plattMPC_FLOAT ub5[12];

    /* diagonal matrix of size [12 x 12] (only the diagonal is stored) */
    plattMPC_FLOAT H6[12];

    /* vector of size 12 */
    plattMPC_FLOAT f6[12];

    /* vector of size 12 */
    plattMPC_FLOAT lb6[12];

    /* vector of size 12 */
    plattMPC_FLOAT ub6[12];

    /* diagonal matrix of size [12 x 12] (only the diagonal is stored) */
    plattMPC_FLOAT H7[12];

    /* vector of size 12 */
    plattMPC_FLOAT f7[12];

    /* vector of size 12 */
    plattMPC_FLOAT lb7[12];

    /* vector of size 12 */
    plattMPC_FLOAT ub7[12];

    /* diagonal matrix of size [12 x 12] (only the diagonal is stored) */
    plattMPC_FLOAT H8[12];

    /* vector of size 12 */
    plattMPC_FLOAT f8[12];

    /* vector of size 12 */
    plattMPC_FLOAT lb8[12];

    /* vector of size 12 */
    plattMPC_FLOAT ub8[12];

    /* diagonal matrix of size [12 x 12] (only the diagonal is stored) */
    plattMPC_FLOAT H9[12];

    /* vector of size 12 */
    plattMPC_FLOAT f9[12];

    /* vector of size 12 */
    plattMPC_FLOAT lb9[12];

    /* vector of size 12 */
    plattMPC_FLOAT ub9[12];

    /* diagonal matrix of size [12 x 12] (only the diagonal is stored) */
    plattMPC_FLOAT H10[12];

    /* vector of size 12 */
    plattMPC_FLOAT f10[12];

    /* vector of size 12 */
    plattMPC_FLOAT lb10[12];

    /* vector of size 12 */
    plattMPC_FLOAT ub10[12];

    /* diagonal matrix of size [12 x 12] (only the diagonal is stored) */
    plattMPC_FLOAT H11[12];

    /* vector of size 12 */
    plattMPC_FLOAT f11[12];

    /* vector of size 12 */
    plattMPC_FLOAT lb11[12];

    /* vector of size 12 */
    plattMPC_FLOAT ub11[12];

    /* diagonal matrix of size [12 x 12] (only the diagonal is stored) */
    plattMPC_FLOAT H12[12];

    /* vector of size 12 */
    plattMPC_FLOAT f12[12];

    /* vector of size 12 */
    plattMPC_FLOAT lb12[12];

    /* vector of size 12 */
    plattMPC_FLOAT ub12[12];

    /* diagonal matrix of size [12 x 12] (only the diagonal is stored) */
    plattMPC_FLOAT H13[12];

    /* vector of size 12 */
    plattMPC_FLOAT f13[12];

    /* vector of size 12 */
    plattMPC_FLOAT lb13[12];

    /* vector of size 12 */
    plattMPC_FLOAT ub13[12];

    /* diagonal matrix of size [12 x 12] (only the diagonal is stored) */
    plattMPC_FLOAT H14[12];

    /* vector of size 12 */
    plattMPC_FLOAT f14[12];

    /* vector of size 12 */
    plattMPC_FLOAT lb14[12];

    /* vector of size 12 */
    plattMPC_FLOAT ub14[12];

    /* diagonal matrix of size [10 x 10] (only the diagonal is stored) */
    plattMPC_FLOAT H15[10];

    /* vector of size 10 */
    plattMPC_FLOAT f15[10];

    /* vector of size 10 */
    plattMPC_FLOAT lb15[10];

    /* vector of size 10 */
    plattMPC_FLOAT ub15[10];

} plattMPC_params;


/* OUTPUTS --------------------------------------------------------------*/
/* the desired variables are put here by the solver */
typedef struct plattMPC_output
{
    /* vector of size 12 */
    plattMPC_FLOAT z1[12];

    /* vector of size 12 */
    plattMPC_FLOAT z2[12];

    /* vector of size 12 */
    plattMPC_FLOAT z3[12];

    /* vector of size 12 */
    plattMPC_FLOAT z4[12];

    /* vector of size 12 */
    plattMPC_FLOAT z5[12];

    /* vector of size 12 */
    plattMPC_FLOAT z6[12];

    /* vector of size 12 */
    plattMPC_FLOAT z7[12];

    /* vector of size 12 */
    plattMPC_FLOAT z8[12];

    /* vector of size 12 */
    plattMPC_FLOAT z9[12];

    /* vector of size 12 */
    plattMPC_FLOAT z10[12];

    /* vector of size 12 */
    plattMPC_FLOAT z11[12];

    /* vector of size 12 */
    plattMPC_FLOAT z12[12];

    /* vector of size 12 */
    plattMPC_FLOAT z13[12];

    /* vector of size 12 */
    plattMPC_FLOAT z14[12];

    /* vector of size 10 */
    plattMPC_FLOAT z15[10];

} plattMPC_output;


/* SOLVER INFO ----------------------------------------------------------*/
/* diagnostic data from last interior point step */
typedef struct plattMPC_info
{
    /* iteration number */
    int it;
	
    /* inf-norm of equality constraint residuals */
    plattMPC_FLOAT res_eq;
	
    /* inf-norm of inequality constraint residuals */
    plattMPC_FLOAT res_ineq;

    /* primal objective */
    plattMPC_FLOAT pobj;	
	
    /* dual objective */
    plattMPC_FLOAT dobj;	

    /* duality gap := pobj - dobj */
    plattMPC_FLOAT dgap;		
	
    /* relative duality gap := |dgap / pobj | */
    plattMPC_FLOAT rdgap;		

    /* duality measure */
    plattMPC_FLOAT mu;

	/* duality measure (after affine step) */
    plattMPC_FLOAT mu_aff;
	
    /* centering parameter */
    plattMPC_FLOAT sigma;
	
    /* number of backtracking line search steps (affine direction) */
    int lsit_aff;
    
    /* number of backtracking line search steps (combined direction) */
    int lsit_cc;
    
    /* step size (affine direction) */
    plattMPC_FLOAT step_aff;
    
    /* step size (combined direction) */
    plattMPC_FLOAT step_cc;    

	/* solvertime */
	plattMPC_FLOAT solvetime;   

} plattMPC_info;


/* SOLVER FUNCTION DEFINITION -------------------------------------------*/
/* examine exitflag before using the result! */
int plattMPC_solve(plattMPC_params* params, plattMPC_output* output, plattMPC_info* info);


#endif