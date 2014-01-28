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

#ifndef __statePenaltyMPC_H__
#define __statePenaltyMPC_H__


/* DATA TYPE ------------------------------------------------------------*/
typedef double statePenaltyMPC_FLOAT;


/* SOLVER SETTINGS ------------------------------------------------------*/
/* print level */
#ifndef statePenaltyMPC_SET_PRINTLEVEL
#define statePenaltyMPC_SET_PRINTLEVEL    (2)
#endif

/* timing */
#ifndef statePenaltyMPC_SET_TIMING
#define statePenaltyMPC_SET_TIMING    (1)
#endif

/* Numeric Warnings */
/* #define PRINTNUMERICALWARNINGS */

/* maximum number of iterations  */
#define statePenaltyMPC_SET_MAXIT         (40)	

/* scaling factor of line search (affine direction) */
#define statePenaltyMPC_SET_LS_SCALE_AFF  (0.9)      

/* scaling factor of line search (combined direction) */
#define statePenaltyMPC_SET_LS_SCALE      (0.95)  

/* minimum required step size in each iteration */
#define statePenaltyMPC_SET_LS_MINSTEP    (1E-08)

/* maximum step size (combined direction) */
#define statePenaltyMPC_SET_LS_MAXSTEP    (0.995)

/* desired relative duality gap */
#define statePenaltyMPC_SET_ACC_RDGAP     (0.0001)

/* desired maximum residual on equality constraints */
#define statePenaltyMPC_SET_ACC_RESEQ     (1E-06)

/* desired maximum residual on inequality constraints */
#define statePenaltyMPC_SET_ACC_RESINEQ   (1E-06)

/* desired maximum violation of complementarity */
#define statePenaltyMPC_SET_ACC_KKTCOMPL  (1E-06)


/* RETURN CODES----------------------------------------------------------*/
/* solver has converged within desired accuracy */
#define statePenaltyMPC_OPTIMAL      (1)

/* maximum number of iterations has been reached */
#define statePenaltyMPC_MAXITREACHED (0)

/* no progress in line search possible */
#define statePenaltyMPC_NOPROGRESS   (-7)




/* PARAMETERS -----------------------------------------------------------*/
/* fill this with data before calling the solver! */
typedef struct statePenaltyMPC_params
{
    /* matrix of size [12 x 12] (column major format) */
    statePenaltyMPC_FLOAT Q1[144];

    /* vector of size 12 */
    statePenaltyMPC_FLOAT f1[12];

    /* vector of size 12 */
    statePenaltyMPC_FLOAT lb1[12];

    /* vector of size 12 */
    statePenaltyMPC_FLOAT ub1[12];

    /* vector of size 6 */
    statePenaltyMPC_FLOAT e1[6];

    /* matrix of size [12 x 12] (column major format) */
    statePenaltyMPC_FLOAT Q2[144];

    /* vector of size 12 */
    statePenaltyMPC_FLOAT f2[12];

    /* vector of size 12 */
    statePenaltyMPC_FLOAT lb2[12];

    /* vector of size 12 */
    statePenaltyMPC_FLOAT ub2[12];

    /* matrix of size [12 x 12] (column major format) */
    statePenaltyMPC_FLOAT Q3[144];

    /* vector of size 12 */
    statePenaltyMPC_FLOAT f3[12];

    /* vector of size 12 */
    statePenaltyMPC_FLOAT lb3[12];

    /* vector of size 12 */
    statePenaltyMPC_FLOAT ub3[12];

    /* matrix of size [12 x 12] (column major format) */
    statePenaltyMPC_FLOAT Q4[144];

    /* vector of size 12 */
    statePenaltyMPC_FLOAT f4[12];

    /* vector of size 12 */
    statePenaltyMPC_FLOAT lb4[12];

    /* vector of size 12 */
    statePenaltyMPC_FLOAT ub4[12];

    /* matrix of size [12 x 12] (column major format) */
    statePenaltyMPC_FLOAT Q5[144];

    /* vector of size 12 */
    statePenaltyMPC_FLOAT f5[12];

    /* vector of size 12 */
    statePenaltyMPC_FLOAT lb5[12];

    /* vector of size 12 */
    statePenaltyMPC_FLOAT ub5[12];

    /* matrix of size [12 x 12] (column major format) */
    statePenaltyMPC_FLOAT Q6[144];

    /* vector of size 12 */
    statePenaltyMPC_FLOAT f6[12];

    /* vector of size 12 */
    statePenaltyMPC_FLOAT lb6[12];

    /* vector of size 12 */
    statePenaltyMPC_FLOAT ub6[12];

    /* matrix of size [12 x 12] (column major format) */
    statePenaltyMPC_FLOAT Q7[144];

    /* vector of size 12 */
    statePenaltyMPC_FLOAT f7[12];

    /* vector of size 12 */
    statePenaltyMPC_FLOAT lb7[12];

    /* vector of size 12 */
    statePenaltyMPC_FLOAT ub7[12];

    /* matrix of size [12 x 12] (column major format) */
    statePenaltyMPC_FLOAT Q8[144];

    /* vector of size 12 */
    statePenaltyMPC_FLOAT f8[12];

    /* vector of size 12 */
    statePenaltyMPC_FLOAT lb8[12];

    /* vector of size 12 */
    statePenaltyMPC_FLOAT ub8[12];

    /* matrix of size [12 x 12] (column major format) */
    statePenaltyMPC_FLOAT Q9[144];

    /* vector of size 12 */
    statePenaltyMPC_FLOAT f9[12];

    /* vector of size 12 */
    statePenaltyMPC_FLOAT lb9[12];

    /* vector of size 12 */
    statePenaltyMPC_FLOAT ub9[12];

    /* matrix of size [12 x 12] (column major format) */
    statePenaltyMPC_FLOAT Q10[144];

    /* vector of size 12 */
    statePenaltyMPC_FLOAT f10[12];

    /* vector of size 12 */
    statePenaltyMPC_FLOAT lb10[12];

    /* vector of size 12 */
    statePenaltyMPC_FLOAT ub10[12];

    /* matrix of size [12 x 12] (column major format) */
    statePenaltyMPC_FLOAT Q11[144];

    /* vector of size 12 */
    statePenaltyMPC_FLOAT f11[12];

    /* vector of size 12 */
    statePenaltyMPC_FLOAT lb11[12];

    /* vector of size 12 */
    statePenaltyMPC_FLOAT ub11[12];

    /* matrix of size [12 x 12] (column major format) */
    statePenaltyMPC_FLOAT Q12[144];

    /* vector of size 12 */
    statePenaltyMPC_FLOAT f12[12];

    /* vector of size 12 */
    statePenaltyMPC_FLOAT lb12[12];

    /* vector of size 12 */
    statePenaltyMPC_FLOAT ub12[12];

    /* matrix of size [12 x 12] (column major format) */
    statePenaltyMPC_FLOAT Q13[144];

    /* vector of size 12 */
    statePenaltyMPC_FLOAT f13[12];

    /* vector of size 12 */
    statePenaltyMPC_FLOAT lb13[12];

    /* vector of size 12 */
    statePenaltyMPC_FLOAT ub13[12];

    /* matrix of size [12 x 12] (column major format) */
    statePenaltyMPC_FLOAT Q14[144];

    /* vector of size 12 */
    statePenaltyMPC_FLOAT f14[12];

    /* vector of size 12 */
    statePenaltyMPC_FLOAT lb14[12];

    /* vector of size 12 */
    statePenaltyMPC_FLOAT ub14[12];

    /* matrix of size [12 x 12] (column major format) */
    statePenaltyMPC_FLOAT Q15[144];

    /* vector of size 12 */
    statePenaltyMPC_FLOAT f15[12];

    /* vector of size 12 */
    statePenaltyMPC_FLOAT lb15[12];

    /* vector of size 6 */
    statePenaltyMPC_FLOAT ub15[6];

    /* matrix of size [12 x 12] (column major format) */
    statePenaltyMPC_FLOAT A15[144];

    /* vector of size 12 */
    statePenaltyMPC_FLOAT b15[12];

} statePenaltyMPC_params;


/* OUTPUTS --------------------------------------------------------------*/
/* the desired variables are put here by the solver */
typedef struct statePenaltyMPC_output
{
    /* vector of size 12 */
    statePenaltyMPC_FLOAT z1[12];

    /* vector of size 12 */
    statePenaltyMPC_FLOAT z2[12];

    /* vector of size 12 */
    statePenaltyMPC_FLOAT z3[12];

    /* vector of size 12 */
    statePenaltyMPC_FLOAT z4[12];

    /* vector of size 12 */
    statePenaltyMPC_FLOAT z5[12];

    /* vector of size 12 */
    statePenaltyMPC_FLOAT z6[12];

    /* vector of size 12 */
    statePenaltyMPC_FLOAT z7[12];

    /* vector of size 12 */
    statePenaltyMPC_FLOAT z8[12];

    /* vector of size 12 */
    statePenaltyMPC_FLOAT z9[12];

    /* vector of size 12 */
    statePenaltyMPC_FLOAT z10[12];

    /* vector of size 12 */
    statePenaltyMPC_FLOAT z11[12];

    /* vector of size 12 */
    statePenaltyMPC_FLOAT z12[12];

    /* vector of size 12 */
    statePenaltyMPC_FLOAT z13[12];

    /* vector of size 12 */
    statePenaltyMPC_FLOAT z14[12];

    /* vector of size 6 */
    statePenaltyMPC_FLOAT z15[6];

} statePenaltyMPC_output;


/* SOLVER INFO ----------------------------------------------------------*/
/* diagnostic data from last interior point step */
typedef struct statePenaltyMPC_info
{
    /* iteration number */
    int it;
	
    /* inf-norm of equality constraint residuals */
    statePenaltyMPC_FLOAT res_eq;
	
    /* inf-norm of inequality constraint residuals */
    statePenaltyMPC_FLOAT res_ineq;

    /* primal objective */
    statePenaltyMPC_FLOAT pobj;	
	
    /* dual objective */
    statePenaltyMPC_FLOAT dobj;	

    /* duality gap := pobj - dobj */
    statePenaltyMPC_FLOAT dgap;		
	
    /* relative duality gap := |dgap / pobj | */
    statePenaltyMPC_FLOAT rdgap;		

    /* duality measure */
    statePenaltyMPC_FLOAT mu;

	/* duality measure (after affine step) */
    statePenaltyMPC_FLOAT mu_aff;
	
    /* centering parameter */
    statePenaltyMPC_FLOAT sigma;
	
    /* number of backtracking line search steps (affine direction) */
    int lsit_aff;
    
    /* number of backtracking line search steps (combined direction) */
    int lsit_cc;
    
    /* step size (affine direction) */
    statePenaltyMPC_FLOAT step_aff;
    
    /* step size (combined direction) */
    statePenaltyMPC_FLOAT step_cc;    

	/* solvertime */
	statePenaltyMPC_FLOAT solvetime;   

} statePenaltyMPC_info;


/* SOLVER FUNCTION DEFINITION -------------------------------------------*/
/* examine exitflag before using the result! */
int statePenaltyMPC_solve(statePenaltyMPC_params* params, statePenaltyMPC_output* output, statePenaltyMPC_info* info);


#endif