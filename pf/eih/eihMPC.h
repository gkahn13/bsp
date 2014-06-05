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

#ifndef __eihMPC_H__
#define __eihMPC_H__


/* DATA TYPE ------------------------------------------------------------*/
typedef double eihMPC_FLOAT;


/* SOLVER SETTINGS ------------------------------------------------------*/
/* print level */
#ifndef eihMPC_SET_PRINTLEVEL
#define eihMPC_SET_PRINTLEVEL    (0)
#endif

/* timing */
#ifndef eihMPC_SET_TIMING
#define eihMPC_SET_TIMING    (0)
#endif

/* Numeric Warnings */
/* #define PRINTNUMERICALWARNINGS */

/* maximum number of iterations  */
#define eihMPC_SET_MAXIT         (50)	

/* scaling factor of line search (affine direction) */
#define eihMPC_SET_LS_SCALE_AFF  (0.9)      

/* scaling factor of line search (combined direction) */
#define eihMPC_SET_LS_SCALE      (0.95)  

/* minimum required step size in each iteration */
#define eihMPC_SET_LS_MINSTEP    (1E-08)

/* maximum step size (combined direction) */
#define eihMPC_SET_LS_MAXSTEP    (0.995)

/* desired relative duality gap */
#define eihMPC_SET_ACC_RDGAP     (0.0001)

/* desired maximum residual on equality constraints */
#define eihMPC_SET_ACC_RESEQ     (1E-06)

/* desired maximum residual on inequality constraints */
#define eihMPC_SET_ACC_RESINEQ   (1E-06)

/* desired maximum violation of complementarity */
#define eihMPC_SET_ACC_KKTCOMPL  (1E-06)


/* RETURN CODES----------------------------------------------------------*/
/* solver has converged within desired accuracy */
#define eihMPC_OPTIMAL      (1)

/* maximum number of iterations has been reached */
#define eihMPC_MAXITREACHED (0)

/* no progress in line search possible */
#define eihMPC_NOPROGRESS   (-7)




/* PARAMETERS -----------------------------------------------------------*/
/* fill this with data before calling the solver! */
typedef struct eihMPC_params
{
    /* diagonal matrix of size [14 x 14] (only the diagonal is stored) */
    eihMPC_FLOAT H1[14];

    /* vector of size 14 */
    eihMPC_FLOAT f1[14];

    /* vector of size 14 */
    eihMPC_FLOAT lb1[14];

    /* vector of size 14 */
    eihMPC_FLOAT ub1[14];

    /* vector of size 7 */
    eihMPC_FLOAT c1[7];

    /* diagonal matrix of size [14 x 14] (only the diagonal is stored) */
    eihMPC_FLOAT H2[14];

    /* vector of size 14 */
    eihMPC_FLOAT f2[14];

    /* vector of size 14 */
    eihMPC_FLOAT lb2[14];

    /* vector of size 14 */
    eihMPC_FLOAT ub2[14];

    /* diagonal matrix of size [14 x 14] (only the diagonal is stored) */
    eihMPC_FLOAT H3[14];

    /* vector of size 14 */
    eihMPC_FLOAT f3[14];

    /* vector of size 14 */
    eihMPC_FLOAT lb3[14];

    /* vector of size 14 */
    eihMPC_FLOAT ub3[14];

    /* diagonal matrix of size [14 x 14] (only the diagonal is stored) */
    eihMPC_FLOAT H4[14];

    /* vector of size 14 */
    eihMPC_FLOAT f4[14];

    /* vector of size 14 */
    eihMPC_FLOAT lb4[14];

    /* vector of size 14 */
    eihMPC_FLOAT ub4[14];

    /* diagonal matrix of size [7 x 7] (only the diagonal is stored) */
    eihMPC_FLOAT H5[7];

    /* vector of size 7 */
    eihMPC_FLOAT f5[7];

    /* vector of size 7 */
    eihMPC_FLOAT lb5[7];

    /* vector of size 7 */
    eihMPC_FLOAT ub5[7];

} eihMPC_params;


/* OUTPUTS --------------------------------------------------------------*/
/* the desired variables are put here by the solver */
typedef struct eihMPC_output
{
    /* vector of size 14 */
    eihMPC_FLOAT z1[14];

    /* vector of size 14 */
    eihMPC_FLOAT z2[14];

    /* vector of size 14 */
    eihMPC_FLOAT z3[14];

    /* vector of size 14 */
    eihMPC_FLOAT z4[14];

    /* vector of size 7 */
    eihMPC_FLOAT z5[7];

} eihMPC_output;


/* SOLVER INFO ----------------------------------------------------------*/
/* diagnostic data from last interior point step */
typedef struct eihMPC_info
{
    /* iteration number */
    int it;
	
    /* inf-norm of equality constraint residuals */
    eihMPC_FLOAT res_eq;
	
    /* inf-norm of inequality constraint residuals */
    eihMPC_FLOAT res_ineq;

    /* primal objective */
    eihMPC_FLOAT pobj;	
	
    /* dual objective */
    eihMPC_FLOAT dobj;	

    /* duality gap := pobj - dobj */
    eihMPC_FLOAT dgap;		
	
    /* relative duality gap := |dgap / pobj | */
    eihMPC_FLOAT rdgap;		

    /* duality measure */
    eihMPC_FLOAT mu;

	/* duality measure (after affine step) */
    eihMPC_FLOAT mu_aff;
	
    /* centering parameter */
    eihMPC_FLOAT sigma;
	
    /* number of backtracking line search steps (affine direction) */
    int lsit_aff;
    
    /* number of backtracking line search steps (combined direction) */
    int lsit_cc;
    
    /* step size (affine direction) */
    eihMPC_FLOAT step_aff;
    
    /* step size (combined direction) */
    eihMPC_FLOAT step_cc;    

	/* solvertime */
	eihMPC_FLOAT solvetime;   

} eihMPC_info;


/* SOLVER FUNCTION DEFINITION -------------------------------------------*/
/* examine exitflag before using the result! */
int eihMPC_solve(eihMPC_params* params, eihMPC_output* output, eihMPC_info* info);


#endif