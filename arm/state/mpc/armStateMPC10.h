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

#ifndef __armStateMPC_H__
#define __armStateMPC_H__


/* DATA TYPE ------------------------------------------------------------*/
typedef double armStateMPC_FLOAT;


/* SOLVER SETTINGS ------------------------------------------------------*/
/* print level */
#ifndef armStateMPC_SET_PRINTLEVEL
#define armStateMPC_SET_PRINTLEVEL    (0)
#endif

/* timing */
#ifndef armStateMPC_SET_TIMING
#define armStateMPC_SET_TIMING    (0)
#endif

/* Numeric Warnings */
/* #define PRINTNUMERICALWARNINGS */

/* maximum number of iterations  */
#define armStateMPC_SET_MAXIT         (30)	

/* scaling factor of line search (affine direction) */
#define armStateMPC_SET_LS_SCALE_AFF  (0.9)      

/* scaling factor of line search (combined direction) */
#define armStateMPC_SET_LS_SCALE      (0.95)  

/* minimum required step size in each iteration */
#define armStateMPC_SET_LS_MINSTEP    (1E-08)

/* maximum step size (combined direction) */
#define armStateMPC_SET_LS_MAXSTEP    (0.995)

/* desired relative duality gap */
#define armStateMPC_SET_ACC_RDGAP     (0.0001)

/* desired maximum residual on equality constraints */
#define armStateMPC_SET_ACC_RESEQ     (1E-06)

/* desired maximum residual on inequality constraints */
#define armStateMPC_SET_ACC_RESINEQ   (1E-06)

/* desired maximum violation of complementarity */
#define armStateMPC_SET_ACC_KKTCOMPL  (1E-06)


/* RETURN CODES----------------------------------------------------------*/
/* solver has converged within desired accuracy */
#define armStateMPC_OPTIMAL      (1)

/* maximum number of iterations has been reached */
#define armStateMPC_MAXITREACHED (0)

/* no progress in line search possible */
#define armStateMPC_NOPROGRESS   (-7)




/* PARAMETERS -----------------------------------------------------------*/
/* fill this with data before calling the solver! */
typedef struct armStateMPC_params
{
    /* diagonal matrix of size [12 x 12] (only the diagonal is stored) */
    armStateMPC_FLOAT H1[12];

    /* vector of size 12 */
    armStateMPC_FLOAT f1[12];

    /* vector of size 12 */
    armStateMPC_FLOAT lb1[12];

    /* vector of size 12 */
    armStateMPC_FLOAT ub1[12];

    /* matrix of size [6 x 12] (column major format) */
    armStateMPC_FLOAT C1[72];

    /* vector of size 6 */
    armStateMPC_FLOAT e1[6];

    /* diagonal matrix of size [12 x 12] (only the diagonal is stored) */
    armStateMPC_FLOAT H2[12];

    /* vector of size 12 */
    armStateMPC_FLOAT f2[12];

    /* vector of size 12 */
    armStateMPC_FLOAT lb2[12];

    /* vector of size 12 */
    armStateMPC_FLOAT ub2[12];

    /* matrix of size [6 x 12] (column major format) */
    armStateMPC_FLOAT C2[72];

    /* vector of size 6 */
    armStateMPC_FLOAT e2[6];

    /* diagonal matrix of size [12 x 12] (only the diagonal is stored) */
    armStateMPC_FLOAT H3[12];

    /* vector of size 12 */
    armStateMPC_FLOAT f3[12];

    /* vector of size 12 */
    armStateMPC_FLOAT lb3[12];

    /* vector of size 12 */
    armStateMPC_FLOAT ub3[12];

    /* matrix of size [6 x 12] (column major format) */
    armStateMPC_FLOAT C3[72];

    /* vector of size 6 */
    armStateMPC_FLOAT e3[6];

    /* diagonal matrix of size [12 x 12] (only the diagonal is stored) */
    armStateMPC_FLOAT H4[12];

    /* vector of size 12 */
    armStateMPC_FLOAT f4[12];

    /* vector of size 12 */
    armStateMPC_FLOAT lb4[12];

    /* vector of size 12 */
    armStateMPC_FLOAT ub4[12];

    /* matrix of size [6 x 12] (column major format) */
    armStateMPC_FLOAT C4[72];

    /* vector of size 6 */
    armStateMPC_FLOAT e4[6];

    /* diagonal matrix of size [12 x 12] (only the diagonal is stored) */
    armStateMPC_FLOAT H5[12];

    /* vector of size 12 */
    armStateMPC_FLOAT f5[12];

    /* vector of size 12 */
    armStateMPC_FLOAT lb5[12];

    /* vector of size 12 */
    armStateMPC_FLOAT ub5[12];

    /* matrix of size [6 x 12] (column major format) */
    armStateMPC_FLOAT C5[72];

    /* vector of size 6 */
    armStateMPC_FLOAT e5[6];

    /* diagonal matrix of size [12 x 12] (only the diagonal is stored) */
    armStateMPC_FLOAT H6[12];

    /* vector of size 12 */
    armStateMPC_FLOAT f6[12];

    /* vector of size 12 */
    armStateMPC_FLOAT lb6[12];

    /* vector of size 12 */
    armStateMPC_FLOAT ub6[12];

    /* matrix of size [6 x 12] (column major format) */
    armStateMPC_FLOAT C6[72];

    /* vector of size 6 */
    armStateMPC_FLOAT e6[6];

    /* diagonal matrix of size [12 x 12] (only the diagonal is stored) */
    armStateMPC_FLOAT H7[12];

    /* vector of size 12 */
    armStateMPC_FLOAT f7[12];

    /* vector of size 12 */
    armStateMPC_FLOAT lb7[12];

    /* vector of size 12 */
    armStateMPC_FLOAT ub7[12];

    /* matrix of size [6 x 12] (column major format) */
    armStateMPC_FLOAT C7[72];

    /* vector of size 6 */
    armStateMPC_FLOAT e7[6];

    /* diagonal matrix of size [12 x 12] (only the diagonal is stored) */
    armStateMPC_FLOAT H8[12];

    /* vector of size 12 */
    armStateMPC_FLOAT f8[12];

    /* vector of size 12 */
    armStateMPC_FLOAT lb8[12];

    /* vector of size 12 */
    armStateMPC_FLOAT ub8[12];

    /* matrix of size [6 x 12] (column major format) */
    armStateMPC_FLOAT C8[72];

    /* vector of size 6 */
    armStateMPC_FLOAT e8[6];

    /* diagonal matrix of size [12 x 12] (only the diagonal is stored) */
    armStateMPC_FLOAT H9[12];

    /* vector of size 12 */
    armStateMPC_FLOAT f9[12];

    /* vector of size 12 */
    armStateMPC_FLOAT lb9[12];

    /* vector of size 12 */
    armStateMPC_FLOAT ub9[12];

    /* matrix of size [6 x 12] (column major format) */
    armStateMPC_FLOAT C9[72];

    /* vector of size 6 */
    armStateMPC_FLOAT e9[6];

    /* diagonal matrix of size [6 x 6] (only the diagonal is stored) */
    armStateMPC_FLOAT H10[6];

    /* vector of size 6 */
    armStateMPC_FLOAT f10[6];

    /* vector of size 6 */
    armStateMPC_FLOAT lb10[6];

    /* vector of size 6 */
    armStateMPC_FLOAT ub10[6];

    /* vector of size 6 */
    armStateMPC_FLOAT e10[6];

} armStateMPC_params;


/* OUTPUTS --------------------------------------------------------------*/
/* the desired variables are put here by the solver */
typedef struct armStateMPC_output
{
    /* vector of size 12 */
    armStateMPC_FLOAT z1[12];

    /* vector of size 12 */
    armStateMPC_FLOAT z2[12];

    /* vector of size 12 */
    armStateMPC_FLOAT z3[12];

    /* vector of size 12 */
    armStateMPC_FLOAT z4[12];

    /* vector of size 12 */
    armStateMPC_FLOAT z5[12];

    /* vector of size 12 */
    armStateMPC_FLOAT z6[12];

    /* vector of size 12 */
    armStateMPC_FLOAT z7[12];

    /* vector of size 12 */
    armStateMPC_FLOAT z8[12];

    /* vector of size 12 */
    armStateMPC_FLOAT z9[12];

    /* vector of size 6 */
    armStateMPC_FLOAT z10[6];

} armStateMPC_output;


/* SOLVER INFO ----------------------------------------------------------*/
/* diagnostic data from last interior point step */
typedef struct armStateMPC_info
{
    /* iteration number */
    int it;
	
    /* inf-norm of equality constraint residuals */
    armStateMPC_FLOAT res_eq;
	
    /* inf-norm of inequality constraint residuals */
    armStateMPC_FLOAT res_ineq;

    /* primal objective */
    armStateMPC_FLOAT pobj;	
	
    /* dual objective */
    armStateMPC_FLOAT dobj;	

    /* duality gap := pobj - dobj */
    armStateMPC_FLOAT dgap;		
	
    /* relative duality gap := |dgap / pobj | */
    armStateMPC_FLOAT rdgap;		

    /* duality measure */
    armStateMPC_FLOAT mu;

	/* duality measure (after affine step) */
    armStateMPC_FLOAT mu_aff;
	
    /* centering parameter */
    armStateMPC_FLOAT sigma;
	
    /* number of backtracking line search steps (affine direction) */
    int lsit_aff;
    
    /* number of backtracking line search steps (combined direction) */
    int lsit_cc;
    
    /* step size (affine direction) */
    armStateMPC_FLOAT step_aff;
    
    /* step size (combined direction) */
    armStateMPC_FLOAT step_cc;    

	/* solvertime */
	armStateMPC_FLOAT solvetime;   

} armStateMPC_info;


/* SOLVER FUNCTION DEFINITION -------------------------------------------*/
/* examine exitflag before using the result! */
int armStateMPC_solve(armStateMPC_params* params, armStateMPC_output* output, armStateMPC_info* info);


#endif