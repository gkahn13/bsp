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
#define beliefPenaltyMPC_SET_PRINTLEVEL    (1)
#endif

/* timing */
#ifndef beliefPenaltyMPC_SET_TIMING
#define beliefPenaltyMPC_SET_TIMING    (1)
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
    /* matrix of size [87 x 87] (column major format) */
    beliefPenaltyMPC_FLOAT H1[7569];

    /* vector of size 87 */
    beliefPenaltyMPC_FLOAT f1[87];

    /* vector of size 87 */
    beliefPenaltyMPC_FLOAT lb1[87];

    /* vector of size 33 */
    beliefPenaltyMPC_FLOAT ub1[33];

    /* matrix of size [27 x 87] (column major format) */
    beliefPenaltyMPC_FLOAT C1[2349];

    /* vector of size 27 */
    beliefPenaltyMPC_FLOAT e1[27];

    /* matrix of size [87 x 87] (column major format) */
    beliefPenaltyMPC_FLOAT H2[7569];

    /* vector of size 87 */
    beliefPenaltyMPC_FLOAT f2[87];

    /* vector of size 87 */
    beliefPenaltyMPC_FLOAT lb2[87];

    /* vector of size 33 */
    beliefPenaltyMPC_FLOAT ub2[33];

    /* matrix of size [27 x 87] (column major format) */
    beliefPenaltyMPC_FLOAT C2[2349];

    /* vector of size 27 */
    beliefPenaltyMPC_FLOAT e2[27];

    /* matrix of size [87 x 87] (column major format) */
    beliefPenaltyMPC_FLOAT H3[7569];

    /* vector of size 87 */
    beliefPenaltyMPC_FLOAT f3[87];

    /* vector of size 87 */
    beliefPenaltyMPC_FLOAT lb3[87];

    /* vector of size 33 */
    beliefPenaltyMPC_FLOAT ub3[33];

    /* matrix of size [27 x 87] (column major format) */
    beliefPenaltyMPC_FLOAT C3[2349];

    /* vector of size 27 */
    beliefPenaltyMPC_FLOAT e3[27];

    /* matrix of size [87 x 87] (column major format) */
    beliefPenaltyMPC_FLOAT H4[7569];

    /* vector of size 87 */
    beliefPenaltyMPC_FLOAT f4[87];

    /* vector of size 87 */
    beliefPenaltyMPC_FLOAT lb4[87];

    /* vector of size 33 */
    beliefPenaltyMPC_FLOAT ub4[33];

    /* matrix of size [27 x 87] (column major format) */
    beliefPenaltyMPC_FLOAT C4[2349];

    /* vector of size 27 */
    beliefPenaltyMPC_FLOAT e4[27];

    /* matrix of size [87 x 87] (column major format) */
    beliefPenaltyMPC_FLOAT H5[7569];

    /* vector of size 87 */
    beliefPenaltyMPC_FLOAT f5[87];

    /* vector of size 87 */
    beliefPenaltyMPC_FLOAT lb5[87];

    /* vector of size 33 */
    beliefPenaltyMPC_FLOAT ub5[33];

    /* matrix of size [27 x 87] (column major format) */
    beliefPenaltyMPC_FLOAT C5[2349];

    /* vector of size 27 */
    beliefPenaltyMPC_FLOAT e5[27];

    /* matrix of size [87 x 87] (column major format) */
    beliefPenaltyMPC_FLOAT H6[7569];

    /* vector of size 87 */
    beliefPenaltyMPC_FLOAT f6[87];

    /* vector of size 87 */
    beliefPenaltyMPC_FLOAT lb6[87];

    /* vector of size 33 */
    beliefPenaltyMPC_FLOAT ub6[33];

    /* matrix of size [27 x 87] (column major format) */
    beliefPenaltyMPC_FLOAT C6[2349];

    /* vector of size 27 */
    beliefPenaltyMPC_FLOAT e6[27];

    /* matrix of size [87 x 87] (column major format) */
    beliefPenaltyMPC_FLOAT H7[7569];

    /* vector of size 87 */
    beliefPenaltyMPC_FLOAT f7[87];

    /* vector of size 87 */
    beliefPenaltyMPC_FLOAT lb7[87];

    /* vector of size 33 */
    beliefPenaltyMPC_FLOAT ub7[33];

    /* matrix of size [27 x 87] (column major format) */
    beliefPenaltyMPC_FLOAT C7[2349];

    /* vector of size 27 */
    beliefPenaltyMPC_FLOAT e7[27];

    /* matrix of size [87 x 87] (column major format) */
    beliefPenaltyMPC_FLOAT H8[7569];

    /* vector of size 87 */
    beliefPenaltyMPC_FLOAT f8[87];

    /* vector of size 87 */
    beliefPenaltyMPC_FLOAT lb8[87];

    /* vector of size 33 */
    beliefPenaltyMPC_FLOAT ub8[33];

    /* matrix of size [27 x 87] (column major format) */
    beliefPenaltyMPC_FLOAT C8[2349];

    /* vector of size 27 */
    beliefPenaltyMPC_FLOAT e8[27];

    /* matrix of size [87 x 87] (column major format) */
    beliefPenaltyMPC_FLOAT H9[7569];

    /* vector of size 87 */
    beliefPenaltyMPC_FLOAT f9[87];

    /* vector of size 87 */
    beliefPenaltyMPC_FLOAT lb9[87];

    /* vector of size 33 */
    beliefPenaltyMPC_FLOAT ub9[33];

    /* matrix of size [27 x 87] (column major format) */
    beliefPenaltyMPC_FLOAT C9[2349];

    /* vector of size 27 */
    beliefPenaltyMPC_FLOAT e9[27];

    /* matrix of size [87 x 87] (column major format) */
    beliefPenaltyMPC_FLOAT H10[7569];

    /* vector of size 87 */
    beliefPenaltyMPC_FLOAT f10[87];

    /* vector of size 87 */
    beliefPenaltyMPC_FLOAT lb10[87];

    /* vector of size 33 */
    beliefPenaltyMPC_FLOAT ub10[33];

    /* matrix of size [27 x 87] (column major format) */
    beliefPenaltyMPC_FLOAT C10[2349];

    /* vector of size 27 */
    beliefPenaltyMPC_FLOAT e10[27];

    /* matrix of size [87 x 87] (column major format) */
    beliefPenaltyMPC_FLOAT H11[7569];

    /* vector of size 87 */
    beliefPenaltyMPC_FLOAT f11[87];

    /* vector of size 87 */
    beliefPenaltyMPC_FLOAT lb11[87];

    /* vector of size 33 */
    beliefPenaltyMPC_FLOAT ub11[33];

    /* matrix of size [27 x 87] (column major format) */
    beliefPenaltyMPC_FLOAT C11[2349];

    /* vector of size 27 */
    beliefPenaltyMPC_FLOAT e11[27];

    /* matrix of size [87 x 87] (column major format) */
    beliefPenaltyMPC_FLOAT H12[7569];

    /* vector of size 87 */
    beliefPenaltyMPC_FLOAT f12[87];

    /* vector of size 87 */
    beliefPenaltyMPC_FLOAT lb12[87];

    /* vector of size 33 */
    beliefPenaltyMPC_FLOAT ub12[33];

    /* matrix of size [27 x 87] (column major format) */
    beliefPenaltyMPC_FLOAT C12[2349];

    /* vector of size 27 */
    beliefPenaltyMPC_FLOAT e12[27];

    /* matrix of size [87 x 87] (column major format) */
    beliefPenaltyMPC_FLOAT H13[7569];

    /* vector of size 87 */
    beliefPenaltyMPC_FLOAT f13[87];

    /* vector of size 87 */
    beliefPenaltyMPC_FLOAT lb13[87];

    /* vector of size 33 */
    beliefPenaltyMPC_FLOAT ub13[33];

    /* matrix of size [27 x 87] (column major format) */
    beliefPenaltyMPC_FLOAT C13[2349];

    /* vector of size 27 */
    beliefPenaltyMPC_FLOAT e13[27];

    /* matrix of size [87 x 87] (column major format) */
    beliefPenaltyMPC_FLOAT H14[7569];

    /* vector of size 87 */
    beliefPenaltyMPC_FLOAT f14[87];

    /* vector of size 87 */
    beliefPenaltyMPC_FLOAT lb14[87];

    /* vector of size 33 */
    beliefPenaltyMPC_FLOAT ub14[33];

    /* matrix of size [30 x 87] (column major format) */
    beliefPenaltyMPC_FLOAT C14[2610];

    /* vector of size 30 */
    beliefPenaltyMPC_FLOAT e14[30];

    /* matrix of size [33 x 33] (column major format) */
    beliefPenaltyMPC_FLOAT H15[1089];

    /* vector of size 33 */
    beliefPenaltyMPC_FLOAT f15[33];

    /* vector of size 33 */
    beliefPenaltyMPC_FLOAT lb15[33];

    /* vector of size 27 */
    beliefPenaltyMPC_FLOAT ub15[27];

    /* vector of size 30 */
    beliefPenaltyMPC_FLOAT e15[30];

    /* matrix of size [30 x 33] (column major format) */
    beliefPenaltyMPC_FLOAT D15[990];

} beliefPenaltyMPC_params;


/* OUTPUTS --------------------------------------------------------------*/
/* the desired variables are put here by the solver */
typedef struct beliefPenaltyMPC_output
{
    /* vector of size 33 */
    beliefPenaltyMPC_FLOAT z1[33];

    /* vector of size 33 */
    beliefPenaltyMPC_FLOAT z2[33];

    /* vector of size 33 */
    beliefPenaltyMPC_FLOAT z3[33];

    /* vector of size 33 */
    beliefPenaltyMPC_FLOAT z4[33];

    /* vector of size 33 */
    beliefPenaltyMPC_FLOAT z5[33];

    /* vector of size 33 */
    beliefPenaltyMPC_FLOAT z6[33];

    /* vector of size 33 */
    beliefPenaltyMPC_FLOAT z7[33];

    /* vector of size 33 */
    beliefPenaltyMPC_FLOAT z8[33];

    /* vector of size 33 */
    beliefPenaltyMPC_FLOAT z9[33];

    /* vector of size 33 */
    beliefPenaltyMPC_FLOAT z10[33];

    /* vector of size 33 */
    beliefPenaltyMPC_FLOAT z11[33];

    /* vector of size 33 */
    beliefPenaltyMPC_FLOAT z12[33];

    /* vector of size 33 */
    beliefPenaltyMPC_FLOAT z13[33];

    /* vector of size 33 */
    beliefPenaltyMPC_FLOAT z14[33];

    /* vector of size 27 */
    beliefPenaltyMPC_FLOAT z15[27];

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