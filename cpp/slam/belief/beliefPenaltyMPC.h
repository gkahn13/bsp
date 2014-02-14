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
#define beliefPenaltyMPC_SET_MAXIT         (50)	

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
    /* diagonal matrix of size [164 x 164] (only the diagonal is stored) */
    beliefPenaltyMPC_FLOAT H1[164];

    /* vector of size 164 */
    beliefPenaltyMPC_FLOAT f1[164];

    /* vector of size 164 */
    beliefPenaltyMPC_FLOAT lb1[164];

    /* vector of size 56 */
    beliefPenaltyMPC_FLOAT ub1[56];

    /* matrix of size [54 x 164] (column major format) */
    beliefPenaltyMPC_FLOAT C1[8856];

    /* vector of size 54 */
    beliefPenaltyMPC_FLOAT e1[54];

    /* diagonal matrix of size [164 x 164] (only the diagonal is stored) */
    beliefPenaltyMPC_FLOAT H2[164];

    /* vector of size 164 */
    beliefPenaltyMPC_FLOAT f2[164];

    /* vector of size 164 */
    beliefPenaltyMPC_FLOAT lb2[164];

    /* vector of size 56 */
    beliefPenaltyMPC_FLOAT ub2[56];

    /* matrix of size [54 x 164] (column major format) */
    beliefPenaltyMPC_FLOAT C2[8856];

    /* vector of size 54 */
    beliefPenaltyMPC_FLOAT e2[54];

    /* diagonal matrix of size [164 x 164] (only the diagonal is stored) */
    beliefPenaltyMPC_FLOAT H3[164];

    /* vector of size 164 */
    beliefPenaltyMPC_FLOAT f3[164];

    /* vector of size 164 */
    beliefPenaltyMPC_FLOAT lb3[164];

    /* vector of size 56 */
    beliefPenaltyMPC_FLOAT ub3[56];

    /* matrix of size [54 x 164] (column major format) */
    beliefPenaltyMPC_FLOAT C3[8856];

    /* vector of size 54 */
    beliefPenaltyMPC_FLOAT e3[54];

    /* diagonal matrix of size [164 x 164] (only the diagonal is stored) */
    beliefPenaltyMPC_FLOAT H4[164];

    /* vector of size 164 */
    beliefPenaltyMPC_FLOAT f4[164];

    /* vector of size 164 */
    beliefPenaltyMPC_FLOAT lb4[164];

    /* vector of size 56 */
    beliefPenaltyMPC_FLOAT ub4[56];

    /* matrix of size [54 x 164] (column major format) */
    beliefPenaltyMPC_FLOAT C4[8856];

    /* vector of size 54 */
    beliefPenaltyMPC_FLOAT e4[54];

    /* diagonal matrix of size [164 x 164] (only the diagonal is stored) */
    beliefPenaltyMPC_FLOAT H5[164];

    /* vector of size 164 */
    beliefPenaltyMPC_FLOAT f5[164];

    /* vector of size 164 */
    beliefPenaltyMPC_FLOAT lb5[164];

    /* vector of size 56 */
    beliefPenaltyMPC_FLOAT ub5[56];

    /* matrix of size [54 x 164] (column major format) */
    beliefPenaltyMPC_FLOAT C5[8856];

    /* vector of size 54 */
    beliefPenaltyMPC_FLOAT e5[54];

    /* diagonal matrix of size [164 x 164] (only the diagonal is stored) */
    beliefPenaltyMPC_FLOAT H6[164];

    /* vector of size 164 */
    beliefPenaltyMPC_FLOAT f6[164];

    /* vector of size 164 */
    beliefPenaltyMPC_FLOAT lb6[164];

    /* vector of size 56 */
    beliefPenaltyMPC_FLOAT ub6[56];

    /* matrix of size [54 x 164] (column major format) */
    beliefPenaltyMPC_FLOAT C6[8856];

    /* vector of size 54 */
    beliefPenaltyMPC_FLOAT e6[54];

    /* diagonal matrix of size [164 x 164] (only the diagonal is stored) */
    beliefPenaltyMPC_FLOAT H7[164];

    /* vector of size 164 */
    beliefPenaltyMPC_FLOAT f7[164];

    /* vector of size 164 */
    beliefPenaltyMPC_FLOAT lb7[164];

    /* vector of size 56 */
    beliefPenaltyMPC_FLOAT ub7[56];

    /* matrix of size [54 x 164] (column major format) */
    beliefPenaltyMPC_FLOAT C7[8856];

    /* vector of size 54 */
    beliefPenaltyMPC_FLOAT e7[54];

    /* diagonal matrix of size [164 x 164] (only the diagonal is stored) */
    beliefPenaltyMPC_FLOAT H8[164];

    /* vector of size 164 */
    beliefPenaltyMPC_FLOAT f8[164];

    /* vector of size 164 */
    beliefPenaltyMPC_FLOAT lb8[164];

    /* vector of size 56 */
    beliefPenaltyMPC_FLOAT ub8[56];

    /* matrix of size [54 x 164] (column major format) */
    beliefPenaltyMPC_FLOAT C8[8856];

    /* vector of size 54 */
    beliefPenaltyMPC_FLOAT e8[54];

    /* diagonal matrix of size [164 x 164] (only the diagonal is stored) */
    beliefPenaltyMPC_FLOAT H9[164];

    /* vector of size 164 */
    beliefPenaltyMPC_FLOAT f9[164];

    /* vector of size 164 */
    beliefPenaltyMPC_FLOAT lb9[164];

    /* vector of size 56 */
    beliefPenaltyMPC_FLOAT ub9[56];

    /* matrix of size [54 x 164] (column major format) */
    beliefPenaltyMPC_FLOAT C9[8856];

    /* vector of size 54 */
    beliefPenaltyMPC_FLOAT e9[54];

    /* diagonal matrix of size [164 x 164] (only the diagonal is stored) */
    beliefPenaltyMPC_FLOAT H10[164];

    /* vector of size 164 */
    beliefPenaltyMPC_FLOAT f10[164];

    /* vector of size 164 */
    beliefPenaltyMPC_FLOAT lb10[164];

    /* vector of size 56 */
    beliefPenaltyMPC_FLOAT ub10[56];

    /* matrix of size [54 x 164] (column major format) */
    beliefPenaltyMPC_FLOAT C10[8856];

    /* vector of size 54 */
    beliefPenaltyMPC_FLOAT e10[54];

    /* diagonal matrix of size [164 x 164] (only the diagonal is stored) */
    beliefPenaltyMPC_FLOAT H11[164];

    /* vector of size 164 */
    beliefPenaltyMPC_FLOAT f11[164];

    /* vector of size 164 */
    beliefPenaltyMPC_FLOAT lb11[164];

    /* vector of size 56 */
    beliefPenaltyMPC_FLOAT ub11[56];

    /* matrix of size [54 x 164] (column major format) */
    beliefPenaltyMPC_FLOAT C11[8856];

    /* vector of size 54 */
    beliefPenaltyMPC_FLOAT e11[54];

    /* diagonal matrix of size [164 x 164] (only the diagonal is stored) */
    beliefPenaltyMPC_FLOAT H12[164];

    /* vector of size 164 */
    beliefPenaltyMPC_FLOAT f12[164];

    /* vector of size 164 */
    beliefPenaltyMPC_FLOAT lb12[164];

    /* vector of size 56 */
    beliefPenaltyMPC_FLOAT ub12[56];

    /* matrix of size [54 x 164] (column major format) */
    beliefPenaltyMPC_FLOAT C12[8856];

    /* vector of size 54 */
    beliefPenaltyMPC_FLOAT e12[54];

    /* diagonal matrix of size [164 x 164] (only the diagonal is stored) */
    beliefPenaltyMPC_FLOAT H13[164];

    /* vector of size 164 */
    beliefPenaltyMPC_FLOAT f13[164];

    /* vector of size 164 */
    beliefPenaltyMPC_FLOAT lb13[164];

    /* vector of size 56 */
    beliefPenaltyMPC_FLOAT ub13[56];

    /* matrix of size [54 x 164] (column major format) */
    beliefPenaltyMPC_FLOAT C13[8856];

    /* vector of size 54 */
    beliefPenaltyMPC_FLOAT e13[54];

    /* diagonal matrix of size [164 x 164] (only the diagonal is stored) */
    beliefPenaltyMPC_FLOAT H14[164];

    /* vector of size 164 */
    beliefPenaltyMPC_FLOAT f14[164];

    /* vector of size 164 */
    beliefPenaltyMPC_FLOAT lb14[164];

    /* vector of size 56 */
    beliefPenaltyMPC_FLOAT ub14[56];

    /* matrix of size [54 x 164] (column major format) */
    beliefPenaltyMPC_FLOAT C14[8856];

    /* vector of size 54 */
    beliefPenaltyMPC_FLOAT e14[54];

    /* diagonal matrix of size [54 x 54] (only the diagonal is stored) */
    beliefPenaltyMPC_FLOAT H15[54];

    /* vector of size 54 */
    beliefPenaltyMPC_FLOAT lb15[54];

    /* vector of size 54 */
    beliefPenaltyMPC_FLOAT ub15[54];

    /* vector of size 54 */
    beliefPenaltyMPC_FLOAT e15[54];

} beliefPenaltyMPC_params;


/* OUTPUTS --------------------------------------------------------------*/
/* the desired variables are put here by the solver */
typedef struct beliefPenaltyMPC_output
{
    /* vector of size 56 */
    beliefPenaltyMPC_FLOAT z1[56];

    /* vector of size 56 */
    beliefPenaltyMPC_FLOAT z2[56];

    /* vector of size 56 */
    beliefPenaltyMPC_FLOAT z3[56];

    /* vector of size 56 */
    beliefPenaltyMPC_FLOAT z4[56];

    /* vector of size 56 */
    beliefPenaltyMPC_FLOAT z5[56];

    /* vector of size 56 */
    beliefPenaltyMPC_FLOAT z6[56];

    /* vector of size 56 */
    beliefPenaltyMPC_FLOAT z7[56];

    /* vector of size 56 */
    beliefPenaltyMPC_FLOAT z8[56];

    /* vector of size 56 */
    beliefPenaltyMPC_FLOAT z9[56];

    /* vector of size 56 */
    beliefPenaltyMPC_FLOAT z10[56];

    /* vector of size 56 */
    beliefPenaltyMPC_FLOAT z11[56];

    /* vector of size 56 */
    beliefPenaltyMPC_FLOAT z12[56];

    /* vector of size 56 */
    beliefPenaltyMPC_FLOAT z13[56];

    /* vector of size 56 */
    beliefPenaltyMPC_FLOAT z14[56];

    /* vector of size 54 */
    beliefPenaltyMPC_FLOAT z15[54];

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