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

#ifndef __smoothMPC_H__
#define __smoothMPC_H__


/* DATA TYPE ------------------------------------------------------------*/
typedef double smoothMPC_FLOAT;


/* SOLVER SETTINGS ------------------------------------------------------*/
/* print level */
#ifndef smoothMPC_SET_PRINTLEVEL
#define smoothMPC_SET_PRINTLEVEL    (0)
#endif

/* timing */
#ifndef smoothMPC_SET_TIMING
#define smoothMPC_SET_TIMING    (0)
#endif

/* Numeric Warnings */
/* #define PRINTNUMERICALWARNINGS */

/* maximum number of iterations  */
#define smoothMPC_SET_MAXIT         (50)	

/* scaling factor of line search (affine direction) */
#define smoothMPC_SET_LS_SCALE_AFF  (0.9)      

/* scaling factor of line search (combined direction) */
#define smoothMPC_SET_LS_SCALE      (0.95)  

/* minimum required step size in each iteration */
#define smoothMPC_SET_LS_MINSTEP    (1E-08)

/* maximum step size (combined direction) */
#define smoothMPC_SET_LS_MAXSTEP    (0.995)

/* desired relative duality gap */
#define smoothMPC_SET_ACC_RDGAP     (0.0001)

/* desired maximum residual on equality constraints */
#define smoothMPC_SET_ACC_RESEQ     (1E-06)

/* desired maximum residual on inequality constraints */
#define smoothMPC_SET_ACC_RESINEQ   (1E-06)

/* desired maximum violation of complementarity */
#define smoothMPC_SET_ACC_KKTCOMPL  (1E-06)


/* RETURN CODES----------------------------------------------------------*/
/* solver has converged within desired accuracy */
#define smoothMPC_OPTIMAL      (1)

/* maximum number of iterations has been reached */
#define smoothMPC_MAXITREACHED (0)

/* no progress in line search possible */
#define smoothMPC_NOPROGRESS   (-7)




/* PARAMETERS -----------------------------------------------------------*/
/* fill this with data before calling the solver! */
typedef struct smoothMPC_params
{
    /* diagonal matrix of size [11 x 11] (only the diagonal is stored) */
    smoothMPC_FLOAT H1[11];

    /* vector of size 11 */
    smoothMPC_FLOAT f1[11];

    /* vector of size 11 */
    smoothMPC_FLOAT lb1[11];

    /* vector of size 5 */
    smoothMPC_FLOAT ub1[5];

    /* matrix of size [3 x 11] (column major format) */
    smoothMPC_FLOAT C1[33];

    /* vector of size 3 */
    smoothMPC_FLOAT e1[3];

    /* diagonal matrix of size [11 x 11] (only the diagonal is stored) */
    smoothMPC_FLOAT H2[11];

    /* vector of size 11 */
    smoothMPC_FLOAT f2[11];

    /* vector of size 11 */
    smoothMPC_FLOAT lb2[11];

    /* vector of size 5 */
    smoothMPC_FLOAT ub2[5];

    /* matrix of size [3 x 11] (column major format) */
    smoothMPC_FLOAT C2[33];

    /* vector of size 3 */
    smoothMPC_FLOAT e2[3];

    /* diagonal matrix of size [11 x 11] (only the diagonal is stored) */
    smoothMPC_FLOAT H3[11];

    /* vector of size 11 */
    smoothMPC_FLOAT f3[11];

    /* vector of size 11 */
    smoothMPC_FLOAT lb3[11];

    /* vector of size 5 */
    smoothMPC_FLOAT ub3[5];

    /* matrix of size [3 x 11] (column major format) */
    smoothMPC_FLOAT C3[33];

    /* vector of size 3 */
    smoothMPC_FLOAT e3[3];

    /* diagonal matrix of size [11 x 11] (only the diagonal is stored) */
    smoothMPC_FLOAT H4[11];

    /* vector of size 11 */
    smoothMPC_FLOAT f4[11];

    /* vector of size 11 */
    smoothMPC_FLOAT lb4[11];

    /* vector of size 5 */
    smoothMPC_FLOAT ub4[5];

    /* matrix of size [3 x 11] (column major format) */
    smoothMPC_FLOAT C4[33];

    /* vector of size 3 */
    smoothMPC_FLOAT e4[3];

    /* diagonal matrix of size [11 x 11] (only the diagonal is stored) */
    smoothMPC_FLOAT H5[11];

    /* vector of size 11 */
    smoothMPC_FLOAT f5[11];

    /* vector of size 11 */
    smoothMPC_FLOAT lb5[11];

    /* vector of size 5 */
    smoothMPC_FLOAT ub5[5];

    /* matrix of size [3 x 11] (column major format) */
    smoothMPC_FLOAT C5[33];

    /* vector of size 3 */
    smoothMPC_FLOAT e5[3];

    /* diagonal matrix of size [11 x 11] (only the diagonal is stored) */
    smoothMPC_FLOAT H6[11];

    /* vector of size 11 */
    smoothMPC_FLOAT f6[11];

    /* vector of size 11 */
    smoothMPC_FLOAT lb6[11];

    /* vector of size 5 */
    smoothMPC_FLOAT ub6[5];

    /* matrix of size [3 x 11] (column major format) */
    smoothMPC_FLOAT C6[33];

    /* vector of size 3 */
    smoothMPC_FLOAT e6[3];

    /* diagonal matrix of size [11 x 11] (only the diagonal is stored) */
    smoothMPC_FLOAT H7[11];

    /* vector of size 11 */
    smoothMPC_FLOAT f7[11];

    /* vector of size 11 */
    smoothMPC_FLOAT lb7[11];

    /* vector of size 5 */
    smoothMPC_FLOAT ub7[5];

    /* matrix of size [3 x 11] (column major format) */
    smoothMPC_FLOAT C7[33];

    /* vector of size 3 */
    smoothMPC_FLOAT e7[3];

    /* diagonal matrix of size [11 x 11] (only the diagonal is stored) */
    smoothMPC_FLOAT H8[11];

    /* vector of size 11 */
    smoothMPC_FLOAT f8[11];

    /* vector of size 11 */
    smoothMPC_FLOAT lb8[11];

    /* vector of size 5 */
    smoothMPC_FLOAT ub8[5];

    /* matrix of size [3 x 11] (column major format) */
    smoothMPC_FLOAT C8[33];

    /* vector of size 3 */
    smoothMPC_FLOAT e8[3];

    /* diagonal matrix of size [11 x 11] (only the diagonal is stored) */
    smoothMPC_FLOAT H9[11];

    /* vector of size 11 */
    smoothMPC_FLOAT f9[11];

    /* vector of size 11 */
    smoothMPC_FLOAT lb9[11];

    /* vector of size 5 */
    smoothMPC_FLOAT ub9[5];

    /* matrix of size [3 x 11] (column major format) */
    smoothMPC_FLOAT C9[33];

    /* vector of size 3 */
    smoothMPC_FLOAT e9[3];

    /* diagonal matrix of size [11 x 11] (only the diagonal is stored) */
    smoothMPC_FLOAT H10[11];

    /* vector of size 11 */
    smoothMPC_FLOAT f10[11];

    /* vector of size 11 */
    smoothMPC_FLOAT lb10[11];

    /* vector of size 5 */
    smoothMPC_FLOAT ub10[5];

    /* matrix of size [3 x 11] (column major format) */
    smoothMPC_FLOAT C10[33];

    /* vector of size 3 */
    smoothMPC_FLOAT e10[3];

    /* diagonal matrix of size [11 x 11] (only the diagonal is stored) */
    smoothMPC_FLOAT H11[11];

    /* vector of size 11 */
    smoothMPC_FLOAT f11[11];

    /* vector of size 11 */
    smoothMPC_FLOAT lb11[11];

    /* vector of size 5 */
    smoothMPC_FLOAT ub11[5];

    /* matrix of size [3 x 11] (column major format) */
    smoothMPC_FLOAT C11[33];

    /* vector of size 3 */
    smoothMPC_FLOAT e11[3];

    /* diagonal matrix of size [11 x 11] (only the diagonal is stored) */
    smoothMPC_FLOAT H12[11];

    /* vector of size 11 */
    smoothMPC_FLOAT f12[11];

    /* vector of size 11 */
    smoothMPC_FLOAT lb12[11];

    /* vector of size 5 */
    smoothMPC_FLOAT ub12[5];

    /* matrix of size [3 x 11] (column major format) */
    smoothMPC_FLOAT C12[33];

    /* vector of size 3 */
    smoothMPC_FLOAT e12[3];

    /* diagonal matrix of size [11 x 11] (only the diagonal is stored) */
    smoothMPC_FLOAT H13[11];

    /* vector of size 11 */
    smoothMPC_FLOAT f13[11];

    /* vector of size 11 */
    smoothMPC_FLOAT lb13[11];

    /* vector of size 5 */
    smoothMPC_FLOAT ub13[5];

    /* matrix of size [3 x 11] (column major format) */
    smoothMPC_FLOAT C13[33];

    /* vector of size 3 */
    smoothMPC_FLOAT e13[3];

    /* diagonal matrix of size [11 x 11] (only the diagonal is stored) */
    smoothMPC_FLOAT H14[11];

    /* vector of size 11 */
    smoothMPC_FLOAT f14[11];

    /* vector of size 11 */
    smoothMPC_FLOAT lb14[11];

    /* vector of size 5 */
    smoothMPC_FLOAT ub14[5];

    /* matrix of size [3 x 11] (column major format) */
    smoothMPC_FLOAT C14[33];

    /* vector of size 3 */
    smoothMPC_FLOAT e14[3];

    /* diagonal matrix of size [11 x 11] (only the diagonal is stored) */
    smoothMPC_FLOAT H15[11];

    /* vector of size 11 */
    smoothMPC_FLOAT f15[11];

    /* vector of size 11 */
    smoothMPC_FLOAT lb15[11];

    /* vector of size 5 */
    smoothMPC_FLOAT ub15[5];

    /* matrix of size [3 x 11] (column major format) */
    smoothMPC_FLOAT C15[33];

    /* vector of size 3 */
    smoothMPC_FLOAT e15[3];

    /* diagonal matrix of size [11 x 11] (only the diagonal is stored) */
    smoothMPC_FLOAT H16[11];

    /* vector of size 11 */
    smoothMPC_FLOAT f16[11];

    /* vector of size 11 */
    smoothMPC_FLOAT lb16[11];

    /* vector of size 5 */
    smoothMPC_FLOAT ub16[5];

    /* matrix of size [3 x 11] (column major format) */
    smoothMPC_FLOAT C16[33];

    /* vector of size 3 */
    smoothMPC_FLOAT e16[3];

    /* diagonal matrix of size [11 x 11] (only the diagonal is stored) */
    smoothMPC_FLOAT H17[11];

    /* vector of size 11 */
    smoothMPC_FLOAT f17[11];

    /* vector of size 11 */
    smoothMPC_FLOAT lb17[11];

    /* vector of size 5 */
    smoothMPC_FLOAT ub17[5];

    /* matrix of size [3 x 11] (column major format) */
    smoothMPC_FLOAT C17[33];

    /* vector of size 3 */
    smoothMPC_FLOAT e17[3];

    /* diagonal matrix of size [11 x 11] (only the diagonal is stored) */
    smoothMPC_FLOAT H18[11];

    /* vector of size 11 */
    smoothMPC_FLOAT f18[11];

    /* vector of size 11 */
    smoothMPC_FLOAT lb18[11];

    /* vector of size 5 */
    smoothMPC_FLOAT ub18[5];

    /* matrix of size [3 x 11] (column major format) */
    smoothMPC_FLOAT C18[33];

    /* vector of size 3 */
    smoothMPC_FLOAT e18[3];

    /* diagonal matrix of size [11 x 11] (only the diagonal is stored) */
    smoothMPC_FLOAT H19[11];

    /* vector of size 11 */
    smoothMPC_FLOAT f19[11];

    /* vector of size 11 */
    smoothMPC_FLOAT lb19[11];

    /* vector of size 5 */
    smoothMPC_FLOAT ub19[5];

    /* matrix of size [3 x 11] (column major format) */
    smoothMPC_FLOAT C19[33];

    /* vector of size 3 */
    smoothMPC_FLOAT e19[3];

    /* diagonal matrix of size [11 x 11] (only the diagonal is stored) */
    smoothMPC_FLOAT H20[11];

    /* vector of size 11 */
    smoothMPC_FLOAT f20[11];

    /* vector of size 11 */
    smoothMPC_FLOAT lb20[11];

    /* vector of size 5 */
    smoothMPC_FLOAT ub20[5];

    /* matrix of size [3 x 11] (column major format) */
    smoothMPC_FLOAT C20[33];

    /* vector of size 3 */
    smoothMPC_FLOAT e20[3];

    /* diagonal matrix of size [11 x 11] (only the diagonal is stored) */
    smoothMPC_FLOAT H21[11];

    /* vector of size 11 */
    smoothMPC_FLOAT f21[11];

    /* vector of size 11 */
    smoothMPC_FLOAT lb21[11];

    /* vector of size 5 */
    smoothMPC_FLOAT ub21[5];

    /* matrix of size [3 x 11] (column major format) */
    smoothMPC_FLOAT C21[33];

    /* vector of size 3 */
    smoothMPC_FLOAT e21[3];

    /* diagonal matrix of size [11 x 11] (only the diagonal is stored) */
    smoothMPC_FLOAT H22[11];

    /* vector of size 11 */
    smoothMPC_FLOAT f22[11];

    /* vector of size 11 */
    smoothMPC_FLOAT lb22[11];

    /* vector of size 5 */
    smoothMPC_FLOAT ub22[5];

    /* matrix of size [3 x 11] (column major format) */
    smoothMPC_FLOAT C22[33];

    /* vector of size 3 */
    smoothMPC_FLOAT e22[3];

    /* diagonal matrix of size [11 x 11] (only the diagonal is stored) */
    smoothMPC_FLOAT H23[11];

    /* vector of size 11 */
    smoothMPC_FLOAT f23[11];

    /* vector of size 11 */
    smoothMPC_FLOAT lb23[11];

    /* vector of size 5 */
    smoothMPC_FLOAT ub23[5];

    /* matrix of size [3 x 11] (column major format) */
    smoothMPC_FLOAT C23[33];

    /* vector of size 3 */
    smoothMPC_FLOAT e23[3];

    /* diagonal matrix of size [11 x 11] (only the diagonal is stored) */
    smoothMPC_FLOAT H24[11];

    /* vector of size 11 */
    smoothMPC_FLOAT f24[11];

    /* vector of size 11 */
    smoothMPC_FLOAT lb24[11];

    /* vector of size 5 */
    smoothMPC_FLOAT ub24[5];

    /* matrix of size [3 x 11] (column major format) */
    smoothMPC_FLOAT C24[33];

    /* vector of size 3 */
    smoothMPC_FLOAT e24[3];

    /* diagonal matrix of size [11 x 11] (only the diagonal is stored) */
    smoothMPC_FLOAT H25[11];

    /* vector of size 11 */
    smoothMPC_FLOAT f25[11];

    /* vector of size 11 */
    smoothMPC_FLOAT lb25[11];

    /* vector of size 5 */
    smoothMPC_FLOAT ub25[5];

    /* matrix of size [3 x 11] (column major format) */
    smoothMPC_FLOAT C25[33];

    /* vector of size 3 */
    smoothMPC_FLOAT e25[3];

    /* diagonal matrix of size [11 x 11] (only the diagonal is stored) */
    smoothMPC_FLOAT H26[11];

    /* vector of size 11 */
    smoothMPC_FLOAT f26[11];

    /* vector of size 11 */
    smoothMPC_FLOAT lb26[11];

    /* vector of size 5 */
    smoothMPC_FLOAT ub26[5];

    /* matrix of size [3 x 11] (column major format) */
    smoothMPC_FLOAT C26[33];

    /* vector of size 3 */
    smoothMPC_FLOAT e26[3];

    /* diagonal matrix of size [11 x 11] (only the diagonal is stored) */
    smoothMPC_FLOAT H27[11];

    /* vector of size 11 */
    smoothMPC_FLOAT f27[11];

    /* vector of size 11 */
    smoothMPC_FLOAT lb27[11];

    /* vector of size 5 */
    smoothMPC_FLOAT ub27[5];

    /* matrix of size [3 x 11] (column major format) */
    smoothMPC_FLOAT C27[33];

    /* vector of size 3 */
    smoothMPC_FLOAT e27[3];

    /* diagonal matrix of size [11 x 11] (only the diagonal is stored) */
    smoothMPC_FLOAT H28[11];

    /* vector of size 11 */
    smoothMPC_FLOAT f28[11];

    /* vector of size 11 */
    smoothMPC_FLOAT lb28[11];

    /* vector of size 5 */
    smoothMPC_FLOAT ub28[5];

    /* matrix of size [3 x 11] (column major format) */
    smoothMPC_FLOAT C28[33];

    /* vector of size 3 */
    smoothMPC_FLOAT e28[3];

    /* diagonal matrix of size [11 x 11] (only the diagonal is stored) */
    smoothMPC_FLOAT H29[11];

    /* vector of size 11 */
    smoothMPC_FLOAT f29[11];

    /* vector of size 11 */
    smoothMPC_FLOAT lb29[11];

    /* vector of size 5 */
    smoothMPC_FLOAT ub29[5];

    /* matrix of size [3 x 11] (column major format) */
    smoothMPC_FLOAT C29[33];

    /* vector of size 3 */
    smoothMPC_FLOAT e29[3];

    /* diagonal matrix of size [11 x 11] (only the diagonal is stored) */
    smoothMPC_FLOAT H30[11];

    /* vector of size 11 */
    smoothMPC_FLOAT f30[11];

    /* vector of size 11 */
    smoothMPC_FLOAT lb30[11];

    /* vector of size 5 */
    smoothMPC_FLOAT ub30[5];

    /* matrix of size [3 x 11] (column major format) */
    smoothMPC_FLOAT C30[33];

    /* vector of size 3 */
    smoothMPC_FLOAT e30[3];

    /* diagonal matrix of size [11 x 11] (only the diagonal is stored) */
    smoothMPC_FLOAT H31[11];

    /* vector of size 11 */
    smoothMPC_FLOAT f31[11];

    /* vector of size 11 */
    smoothMPC_FLOAT lb31[11];

    /* vector of size 5 */
    smoothMPC_FLOAT ub31[5];

    /* matrix of size [3 x 11] (column major format) */
    smoothMPC_FLOAT C31[33];

    /* vector of size 3 */
    smoothMPC_FLOAT e31[3];

    /* diagonal matrix of size [11 x 11] (only the diagonal is stored) */
    smoothMPC_FLOAT H32[11];

    /* vector of size 11 */
    smoothMPC_FLOAT f32[11];

    /* vector of size 11 */
    smoothMPC_FLOAT lb32[11];

    /* vector of size 5 */
    smoothMPC_FLOAT ub32[5];

    /* matrix of size [3 x 11] (column major format) */
    smoothMPC_FLOAT C32[33];

    /* vector of size 3 */
    smoothMPC_FLOAT e32[3];

    /* diagonal matrix of size [11 x 11] (only the diagonal is stored) */
    smoothMPC_FLOAT H33[11];

    /* vector of size 11 */
    smoothMPC_FLOAT f33[11];

    /* vector of size 11 */
    smoothMPC_FLOAT lb33[11];

    /* vector of size 5 */
    smoothMPC_FLOAT ub33[5];

    /* matrix of size [3 x 11] (column major format) */
    smoothMPC_FLOAT C33[33];

    /* vector of size 3 */
    smoothMPC_FLOAT e33[3];

    /* diagonal matrix of size [11 x 11] (only the diagonal is stored) */
    smoothMPC_FLOAT H34[11];

    /* vector of size 11 */
    smoothMPC_FLOAT f34[11];

    /* vector of size 11 */
    smoothMPC_FLOAT lb34[11];

    /* vector of size 5 */
    smoothMPC_FLOAT ub34[5];

    /* matrix of size [3 x 11] (column major format) */
    smoothMPC_FLOAT C34[33];

    /* vector of size 3 */
    smoothMPC_FLOAT e34[3];

    /* diagonal matrix of size [11 x 11] (only the diagonal is stored) */
    smoothMPC_FLOAT H35[11];

    /* vector of size 11 */
    smoothMPC_FLOAT f35[11];

    /* vector of size 11 */
    smoothMPC_FLOAT lb35[11];

    /* vector of size 5 */
    smoothMPC_FLOAT ub35[5];

    /* matrix of size [3 x 11] (column major format) */
    smoothMPC_FLOAT C35[33];

    /* vector of size 3 */
    smoothMPC_FLOAT e35[3];

    /* diagonal matrix of size [11 x 11] (only the diagonal is stored) */
    smoothMPC_FLOAT H36[11];

    /* vector of size 11 */
    smoothMPC_FLOAT f36[11];

    /* vector of size 11 */
    smoothMPC_FLOAT lb36[11];

    /* vector of size 5 */
    smoothMPC_FLOAT ub36[5];

    /* matrix of size [3 x 11] (column major format) */
    smoothMPC_FLOAT C36[33];

    /* vector of size 3 */
    smoothMPC_FLOAT e36[3];

    /* diagonal matrix of size [11 x 11] (only the diagonal is stored) */
    smoothMPC_FLOAT H37[11];

    /* vector of size 11 */
    smoothMPC_FLOAT f37[11];

    /* vector of size 11 */
    smoothMPC_FLOAT lb37[11];

    /* vector of size 5 */
    smoothMPC_FLOAT ub37[5];

    /* matrix of size [3 x 11] (column major format) */
    smoothMPC_FLOAT C37[33];

    /* vector of size 3 */
    smoothMPC_FLOAT e37[3];

    /* diagonal matrix of size [11 x 11] (only the diagonal is stored) */
    smoothMPC_FLOAT H38[11];

    /* vector of size 11 */
    smoothMPC_FLOAT f38[11];

    /* vector of size 11 */
    smoothMPC_FLOAT lb38[11];

    /* vector of size 5 */
    smoothMPC_FLOAT ub38[5];

    /* matrix of size [3 x 11] (column major format) */
    smoothMPC_FLOAT C38[33];

    /* vector of size 3 */
    smoothMPC_FLOAT e38[3];

    /* diagonal matrix of size [11 x 11] (only the diagonal is stored) */
    smoothMPC_FLOAT H39[11];

    /* vector of size 11 */
    smoothMPC_FLOAT f39[11];

    /* vector of size 11 */
    smoothMPC_FLOAT lb39[11];

    /* vector of size 5 */
    smoothMPC_FLOAT ub39[5];

    /* matrix of size [3 x 11] (column major format) */
    smoothMPC_FLOAT C39[33];

    /* vector of size 3 */
    smoothMPC_FLOAT e39[3];

    /* diagonal matrix of size [11 x 11] (only the diagonal is stored) */
    smoothMPC_FLOAT H40[11];

    /* vector of size 11 */
    smoothMPC_FLOAT f40[11];

    /* vector of size 11 */
    smoothMPC_FLOAT lb40[11];

    /* vector of size 5 */
    smoothMPC_FLOAT ub40[5];

    /* matrix of size [3 x 11] (column major format) */
    smoothMPC_FLOAT C40[33];

    /* vector of size 3 */
    smoothMPC_FLOAT e40[3];

    /* diagonal matrix of size [11 x 11] (only the diagonal is stored) */
    smoothMPC_FLOAT H41[11];

    /* vector of size 11 */
    smoothMPC_FLOAT f41[11];

    /* vector of size 11 */
    smoothMPC_FLOAT lb41[11];

    /* vector of size 5 */
    smoothMPC_FLOAT ub41[5];

    /* matrix of size [3 x 11] (column major format) */
    smoothMPC_FLOAT C41[33];

    /* vector of size 3 */
    smoothMPC_FLOAT e41[3];

    /* diagonal matrix of size [11 x 11] (only the diagonal is stored) */
    smoothMPC_FLOAT H42[11];

    /* vector of size 11 */
    smoothMPC_FLOAT f42[11];

    /* vector of size 11 */
    smoothMPC_FLOAT lb42[11];

    /* vector of size 5 */
    smoothMPC_FLOAT ub42[5];

    /* matrix of size [3 x 11] (column major format) */
    smoothMPC_FLOAT C42[33];

    /* vector of size 3 */
    smoothMPC_FLOAT e42[3];

    /* diagonal matrix of size [11 x 11] (only the diagonal is stored) */
    smoothMPC_FLOAT H43[11];

    /* vector of size 11 */
    smoothMPC_FLOAT f43[11];

    /* vector of size 11 */
    smoothMPC_FLOAT lb43[11];

    /* vector of size 5 */
    smoothMPC_FLOAT ub43[5];

    /* matrix of size [3 x 11] (column major format) */
    smoothMPC_FLOAT C43[33];

    /* vector of size 3 */
    smoothMPC_FLOAT e43[3];

    /* diagonal matrix of size [11 x 11] (only the diagonal is stored) */
    smoothMPC_FLOAT H44[11];

    /* vector of size 11 */
    smoothMPC_FLOAT f44[11];

    /* vector of size 11 */
    smoothMPC_FLOAT lb44[11];

    /* vector of size 5 */
    smoothMPC_FLOAT ub44[5];

    /* matrix of size [3 x 11] (column major format) */
    smoothMPC_FLOAT C44[33];

    /* vector of size 3 */
    smoothMPC_FLOAT e44[3];

    /* diagonal matrix of size [11 x 11] (only the diagonal is stored) */
    smoothMPC_FLOAT H45[11];

    /* vector of size 11 */
    smoothMPC_FLOAT f45[11];

    /* vector of size 11 */
    smoothMPC_FLOAT lb45[11];

    /* vector of size 5 */
    smoothMPC_FLOAT ub45[5];

    /* matrix of size [3 x 11] (column major format) */
    smoothMPC_FLOAT C45[33];

    /* vector of size 3 */
    smoothMPC_FLOAT e45[3];

    /* diagonal matrix of size [11 x 11] (only the diagonal is stored) */
    smoothMPC_FLOAT H46[11];

    /* vector of size 11 */
    smoothMPC_FLOAT f46[11];

    /* vector of size 11 */
    smoothMPC_FLOAT lb46[11];

    /* vector of size 5 */
    smoothMPC_FLOAT ub46[5];

    /* matrix of size [3 x 11] (column major format) */
    smoothMPC_FLOAT C46[33];

    /* vector of size 3 */
    smoothMPC_FLOAT e46[3];

    /* diagonal matrix of size [11 x 11] (only the diagonal is stored) */
    smoothMPC_FLOAT H47[11];

    /* vector of size 11 */
    smoothMPC_FLOAT f47[11];

    /* vector of size 11 */
    smoothMPC_FLOAT lb47[11];

    /* vector of size 5 */
    smoothMPC_FLOAT ub47[5];

    /* matrix of size [3 x 11] (column major format) */
    smoothMPC_FLOAT C47[33];

    /* vector of size 3 */
    smoothMPC_FLOAT e47[3];

    /* diagonal matrix of size [11 x 11] (only the diagonal is stored) */
    smoothMPC_FLOAT H48[11];

    /* vector of size 11 */
    smoothMPC_FLOAT f48[11];

    /* vector of size 11 */
    smoothMPC_FLOAT lb48[11];

    /* vector of size 5 */
    smoothMPC_FLOAT ub48[5];

    /* matrix of size [3 x 11] (column major format) */
    smoothMPC_FLOAT C48[33];

    /* vector of size 3 */
    smoothMPC_FLOAT e48[3];

    /* diagonal matrix of size [11 x 11] (only the diagonal is stored) */
    smoothMPC_FLOAT H49[11];

    /* vector of size 11 */
    smoothMPC_FLOAT f49[11];

    /* vector of size 11 */
    smoothMPC_FLOAT lb49[11];

    /* vector of size 5 */
    smoothMPC_FLOAT ub49[5];

    /* matrix of size [3 x 11] (column major format) */
    smoothMPC_FLOAT C49[33];

    /* vector of size 3 */
    smoothMPC_FLOAT e49[3];

    /* diagonal matrix of size [11 x 11] (only the diagonal is stored) */
    smoothMPC_FLOAT H50[11];

    /* vector of size 11 */
    smoothMPC_FLOAT f50[11];

    /* vector of size 11 */
    smoothMPC_FLOAT lb50[11];

    /* vector of size 5 */
    smoothMPC_FLOAT ub50[5];

    /* matrix of size [3 x 11] (column major format) */
    smoothMPC_FLOAT C50[33];

    /* vector of size 3 */
    smoothMPC_FLOAT e50[3];

    /* diagonal matrix of size [11 x 11] (only the diagonal is stored) */
    smoothMPC_FLOAT H51[11];

    /* vector of size 11 */
    smoothMPC_FLOAT f51[11];

    /* vector of size 11 */
    smoothMPC_FLOAT lb51[11];

    /* vector of size 5 */
    smoothMPC_FLOAT ub51[5];

    /* matrix of size [3 x 11] (column major format) */
    smoothMPC_FLOAT C51[33];

    /* vector of size 3 */
    smoothMPC_FLOAT e51[3];

    /* diagonal matrix of size [11 x 11] (only the diagonal is stored) */
    smoothMPC_FLOAT H52[11];

    /* vector of size 11 */
    smoothMPC_FLOAT f52[11];

    /* vector of size 11 */
    smoothMPC_FLOAT lb52[11];

    /* vector of size 5 */
    smoothMPC_FLOAT ub52[5];

    /* matrix of size [3 x 11] (column major format) */
    smoothMPC_FLOAT C52[33];

    /* vector of size 3 */
    smoothMPC_FLOAT e52[3];

    /* diagonal matrix of size [11 x 11] (only the diagonal is stored) */
    smoothMPC_FLOAT H53[11];

    /* vector of size 11 */
    smoothMPC_FLOAT f53[11];

    /* vector of size 11 */
    smoothMPC_FLOAT lb53[11];

    /* vector of size 5 */
    smoothMPC_FLOAT ub53[5];

    /* matrix of size [3 x 11] (column major format) */
    smoothMPC_FLOAT C53[33];

    /* vector of size 3 */
    smoothMPC_FLOAT e53[3];

    /* diagonal matrix of size [11 x 11] (only the diagonal is stored) */
    smoothMPC_FLOAT H54[11];

    /* vector of size 11 */
    smoothMPC_FLOAT f54[11];

    /* vector of size 11 */
    smoothMPC_FLOAT lb54[11];

    /* vector of size 5 */
    smoothMPC_FLOAT ub54[5];

    /* matrix of size [3 x 11] (column major format) */
    smoothMPC_FLOAT C54[33];

    /* vector of size 3 */
    smoothMPC_FLOAT e54[3];

    /* diagonal matrix of size [11 x 11] (only the diagonal is stored) */
    smoothMPC_FLOAT H55[11];

    /* vector of size 11 */
    smoothMPC_FLOAT f55[11];

    /* vector of size 11 */
    smoothMPC_FLOAT lb55[11];

    /* vector of size 5 */
    smoothMPC_FLOAT ub55[5];

    /* matrix of size [3 x 11] (column major format) */
    smoothMPC_FLOAT C55[33];

    /* vector of size 3 */
    smoothMPC_FLOAT e55[3];

    /* diagonal matrix of size [11 x 11] (only the diagonal is stored) */
    smoothMPC_FLOAT H56[11];

    /* vector of size 11 */
    smoothMPC_FLOAT f56[11];

    /* vector of size 11 */
    smoothMPC_FLOAT lb56[11];

    /* vector of size 5 */
    smoothMPC_FLOAT ub56[5];

    /* matrix of size [3 x 11] (column major format) */
    smoothMPC_FLOAT C56[33];

    /* vector of size 3 */
    smoothMPC_FLOAT e56[3];

    /* diagonal matrix of size [11 x 11] (only the diagonal is stored) */
    smoothMPC_FLOAT H57[11];

    /* vector of size 11 */
    smoothMPC_FLOAT f57[11];

    /* vector of size 11 */
    smoothMPC_FLOAT lb57[11];

    /* vector of size 5 */
    smoothMPC_FLOAT ub57[5];

    /* matrix of size [3 x 11] (column major format) */
    smoothMPC_FLOAT C57[33];

    /* vector of size 3 */
    smoothMPC_FLOAT e57[3];

    /* diagonal matrix of size [11 x 11] (only the diagonal is stored) */
    smoothMPC_FLOAT H58[11];

    /* vector of size 11 */
    smoothMPC_FLOAT f58[11];

    /* vector of size 11 */
    smoothMPC_FLOAT lb58[11];

    /* vector of size 5 */
    smoothMPC_FLOAT ub58[5];

    /* matrix of size [3 x 11] (column major format) */
    smoothMPC_FLOAT C58[33];

    /* vector of size 3 */
    smoothMPC_FLOAT e58[3];

    /* diagonal matrix of size [11 x 11] (only the diagonal is stored) */
    smoothMPC_FLOAT H59[11];

    /* vector of size 11 */
    smoothMPC_FLOAT f59[11];

    /* vector of size 11 */
    smoothMPC_FLOAT lb59[11];

    /* vector of size 5 */
    smoothMPC_FLOAT ub59[5];

    /* matrix of size [3 x 11] (column major format) */
    smoothMPC_FLOAT C59[33];

    /* vector of size 3 */
    smoothMPC_FLOAT e59[3];

    /* diagonal matrix of size [3 x 3] (only the diagonal is stored) */
    smoothMPC_FLOAT H60[3];

    /* vector of size 3 */
    smoothMPC_FLOAT f60[3];

    /* vector of size 3 */
    smoothMPC_FLOAT lb60[3];

    /* vector of size 3 */
    smoothMPC_FLOAT ub60[3];

    /* vector of size 3 */
    smoothMPC_FLOAT e60[3];

} smoothMPC_params;


/* OUTPUTS --------------------------------------------------------------*/
/* the desired variables are put here by the solver */
typedef struct smoothMPC_output
{
    /* vector of size 5 */
    smoothMPC_FLOAT z1[5];

    /* vector of size 5 */
    smoothMPC_FLOAT z2[5];

    /* vector of size 5 */
    smoothMPC_FLOAT z3[5];

    /* vector of size 5 */
    smoothMPC_FLOAT z4[5];

    /* vector of size 5 */
    smoothMPC_FLOAT z5[5];

    /* vector of size 5 */
    smoothMPC_FLOAT z6[5];

    /* vector of size 5 */
    smoothMPC_FLOAT z7[5];

    /* vector of size 5 */
    smoothMPC_FLOAT z8[5];

    /* vector of size 5 */
    smoothMPC_FLOAT z9[5];

    /* vector of size 5 */
    smoothMPC_FLOAT z10[5];

    /* vector of size 5 */
    smoothMPC_FLOAT z11[5];

    /* vector of size 5 */
    smoothMPC_FLOAT z12[5];

    /* vector of size 5 */
    smoothMPC_FLOAT z13[5];

    /* vector of size 5 */
    smoothMPC_FLOAT z14[5];

    /* vector of size 5 */
    smoothMPC_FLOAT z15[5];

    /* vector of size 5 */
    smoothMPC_FLOAT z16[5];

    /* vector of size 5 */
    smoothMPC_FLOAT z17[5];

    /* vector of size 5 */
    smoothMPC_FLOAT z18[5];

    /* vector of size 5 */
    smoothMPC_FLOAT z19[5];

    /* vector of size 5 */
    smoothMPC_FLOAT z20[5];

    /* vector of size 5 */
    smoothMPC_FLOAT z21[5];

    /* vector of size 5 */
    smoothMPC_FLOAT z22[5];

    /* vector of size 5 */
    smoothMPC_FLOAT z23[5];

    /* vector of size 5 */
    smoothMPC_FLOAT z24[5];

    /* vector of size 5 */
    smoothMPC_FLOAT z25[5];

    /* vector of size 5 */
    smoothMPC_FLOAT z26[5];

    /* vector of size 5 */
    smoothMPC_FLOAT z27[5];

    /* vector of size 5 */
    smoothMPC_FLOAT z28[5];

    /* vector of size 5 */
    smoothMPC_FLOAT z29[5];

    /* vector of size 5 */
    smoothMPC_FLOAT z30[5];

    /* vector of size 5 */
    smoothMPC_FLOAT z31[5];

    /* vector of size 5 */
    smoothMPC_FLOAT z32[5];

    /* vector of size 5 */
    smoothMPC_FLOAT z33[5];

    /* vector of size 5 */
    smoothMPC_FLOAT z34[5];

    /* vector of size 5 */
    smoothMPC_FLOAT z35[5];

    /* vector of size 5 */
    smoothMPC_FLOAT z36[5];

    /* vector of size 5 */
    smoothMPC_FLOAT z37[5];

    /* vector of size 5 */
    smoothMPC_FLOAT z38[5];

    /* vector of size 5 */
    smoothMPC_FLOAT z39[5];

    /* vector of size 5 */
    smoothMPC_FLOAT z40[5];

    /* vector of size 5 */
    smoothMPC_FLOAT z41[5];

    /* vector of size 5 */
    smoothMPC_FLOAT z42[5];

    /* vector of size 5 */
    smoothMPC_FLOAT z43[5];

    /* vector of size 5 */
    smoothMPC_FLOAT z44[5];

    /* vector of size 5 */
    smoothMPC_FLOAT z45[5];

    /* vector of size 5 */
    smoothMPC_FLOAT z46[5];

    /* vector of size 5 */
    smoothMPC_FLOAT z47[5];

    /* vector of size 5 */
    smoothMPC_FLOAT z48[5];

    /* vector of size 5 */
    smoothMPC_FLOAT z49[5];

    /* vector of size 5 */
    smoothMPC_FLOAT z50[5];

    /* vector of size 5 */
    smoothMPC_FLOAT z51[5];

    /* vector of size 5 */
    smoothMPC_FLOAT z52[5];

    /* vector of size 5 */
    smoothMPC_FLOAT z53[5];

    /* vector of size 5 */
    smoothMPC_FLOAT z54[5];

    /* vector of size 5 */
    smoothMPC_FLOAT z55[5];

    /* vector of size 5 */
    smoothMPC_FLOAT z56[5];

    /* vector of size 5 */
    smoothMPC_FLOAT z57[5];

    /* vector of size 5 */
    smoothMPC_FLOAT z58[5];

    /* vector of size 5 */
    smoothMPC_FLOAT z59[5];

    /* vector of size 3 */
    smoothMPC_FLOAT z60[3];

} smoothMPC_output;


/* SOLVER INFO ----------------------------------------------------------*/
/* diagnostic data from last interior point step */
typedef struct smoothMPC_info
{
    /* iteration number */
    int it;
	
    /* inf-norm of equality constraint residuals */
    smoothMPC_FLOAT res_eq;
	
    /* inf-norm of inequality constraint residuals */
    smoothMPC_FLOAT res_ineq;

    /* primal objective */
    smoothMPC_FLOAT pobj;	
	
    /* dual objective */
    smoothMPC_FLOAT dobj;	

    /* duality gap := pobj - dobj */
    smoothMPC_FLOAT dgap;		
	
    /* relative duality gap := |dgap / pobj | */
    smoothMPC_FLOAT rdgap;		

    /* duality measure */
    smoothMPC_FLOAT mu;

	/* duality measure (after affine step) */
    smoothMPC_FLOAT mu_aff;
	
    /* centering parameter */
    smoothMPC_FLOAT sigma;
	
    /* number of backtracking line search steps (affine direction) */
    int lsit_aff;
    
    /* number of backtracking line search steps (combined direction) */
    int lsit_cc;
    
    /* step size (affine direction) */
    smoothMPC_FLOAT step_aff;
    
    /* step size (combined direction) */
    smoothMPC_FLOAT step_cc;    

	/* solvertime */
	smoothMPC_FLOAT solvetime;   

} smoothMPC_info;


/* SOLVER FUNCTION DEFINITION -------------------------------------------*/
/* examine exitflag before using the result! */
int smoothMPC_solve(smoothMPC_params* params, smoothMPC_output* output, smoothMPC_info* info);


#endif