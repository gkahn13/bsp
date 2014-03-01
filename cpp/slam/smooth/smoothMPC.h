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

    /* vector of size 11 */
    smoothMPC_FLOAT f60[11];

    /* vector of size 11 */
    smoothMPC_FLOAT lb60[11];

    /* vector of size 5 */
    smoothMPC_FLOAT ub60[5];

    /* matrix of size [3 x 11] (column major format) */
    smoothMPC_FLOAT C60[33];

    /* vector of size 3 */
    smoothMPC_FLOAT e60[3];

    /* vector of size 11 */
    smoothMPC_FLOAT f61[11];

    /* vector of size 11 */
    smoothMPC_FLOAT lb61[11];

    /* vector of size 5 */
    smoothMPC_FLOAT ub61[5];

    /* matrix of size [3 x 11] (column major format) */
    smoothMPC_FLOAT C61[33];

    /* vector of size 3 */
    smoothMPC_FLOAT e61[3];

    /* vector of size 11 */
    smoothMPC_FLOAT f62[11];

    /* vector of size 11 */
    smoothMPC_FLOAT lb62[11];

    /* vector of size 5 */
    smoothMPC_FLOAT ub62[5];

    /* matrix of size [3 x 11] (column major format) */
    smoothMPC_FLOAT C62[33];

    /* vector of size 3 */
    smoothMPC_FLOAT e62[3];

    /* vector of size 11 */
    smoothMPC_FLOAT f63[11];

    /* vector of size 11 */
    smoothMPC_FLOAT lb63[11];

    /* vector of size 5 */
    smoothMPC_FLOAT ub63[5];

    /* matrix of size [3 x 11] (column major format) */
    smoothMPC_FLOAT C63[33];

    /* vector of size 3 */
    smoothMPC_FLOAT e63[3];

    /* vector of size 11 */
    smoothMPC_FLOAT f64[11];

    /* vector of size 11 */
    smoothMPC_FLOAT lb64[11];

    /* vector of size 5 */
    smoothMPC_FLOAT ub64[5];

    /* matrix of size [3 x 11] (column major format) */
    smoothMPC_FLOAT C64[33];

    /* vector of size 3 */
    smoothMPC_FLOAT e64[3];

    /* vector of size 11 */
    smoothMPC_FLOAT f65[11];

    /* vector of size 11 */
    smoothMPC_FLOAT lb65[11];

    /* vector of size 5 */
    smoothMPC_FLOAT ub65[5];

    /* matrix of size [3 x 11] (column major format) */
    smoothMPC_FLOAT C65[33];

    /* vector of size 3 */
    smoothMPC_FLOAT e65[3];

    /* vector of size 11 */
    smoothMPC_FLOAT f66[11];

    /* vector of size 11 */
    smoothMPC_FLOAT lb66[11];

    /* vector of size 5 */
    smoothMPC_FLOAT ub66[5];

    /* matrix of size [3 x 11] (column major format) */
    smoothMPC_FLOAT C66[33];

    /* vector of size 3 */
    smoothMPC_FLOAT e66[3];

    /* vector of size 11 */
    smoothMPC_FLOAT f67[11];

    /* vector of size 11 */
    smoothMPC_FLOAT lb67[11];

    /* vector of size 5 */
    smoothMPC_FLOAT ub67[5];

    /* matrix of size [3 x 11] (column major format) */
    smoothMPC_FLOAT C67[33];

    /* vector of size 3 */
    smoothMPC_FLOAT e67[3];

    /* vector of size 11 */
    smoothMPC_FLOAT f68[11];

    /* vector of size 11 */
    smoothMPC_FLOAT lb68[11];

    /* vector of size 5 */
    smoothMPC_FLOAT ub68[5];

    /* matrix of size [3 x 11] (column major format) */
    smoothMPC_FLOAT C68[33];

    /* vector of size 3 */
    smoothMPC_FLOAT e68[3];

    /* vector of size 11 */
    smoothMPC_FLOAT f69[11];

    /* vector of size 11 */
    smoothMPC_FLOAT lb69[11];

    /* vector of size 5 */
    smoothMPC_FLOAT ub69[5];

    /* matrix of size [3 x 11] (column major format) */
    smoothMPC_FLOAT C69[33];

    /* vector of size 3 */
    smoothMPC_FLOAT e69[3];

    /* vector of size 11 */
    smoothMPC_FLOAT f70[11];

    /* vector of size 11 */
    smoothMPC_FLOAT lb70[11];

    /* vector of size 5 */
    smoothMPC_FLOAT ub70[5];

    /* matrix of size [3 x 11] (column major format) */
    smoothMPC_FLOAT C70[33];

    /* vector of size 3 */
    smoothMPC_FLOAT e70[3];

    /* vector of size 11 */
    smoothMPC_FLOAT f71[11];

    /* vector of size 11 */
    smoothMPC_FLOAT lb71[11];

    /* vector of size 5 */
    smoothMPC_FLOAT ub71[5];

    /* matrix of size [3 x 11] (column major format) */
    smoothMPC_FLOAT C71[33];

    /* vector of size 3 */
    smoothMPC_FLOAT e71[3];

    /* vector of size 11 */
    smoothMPC_FLOAT f72[11];

    /* vector of size 11 */
    smoothMPC_FLOAT lb72[11];

    /* vector of size 5 */
    smoothMPC_FLOAT ub72[5];

    /* matrix of size [3 x 11] (column major format) */
    smoothMPC_FLOAT C72[33];

    /* vector of size 3 */
    smoothMPC_FLOAT e72[3];

    /* vector of size 11 */
    smoothMPC_FLOAT f73[11];

    /* vector of size 11 */
    smoothMPC_FLOAT lb73[11];

    /* vector of size 5 */
    smoothMPC_FLOAT ub73[5];

    /* matrix of size [3 x 11] (column major format) */
    smoothMPC_FLOAT C73[33];

    /* vector of size 3 */
    smoothMPC_FLOAT e73[3];

    /* vector of size 11 */
    smoothMPC_FLOAT f74[11];

    /* vector of size 11 */
    smoothMPC_FLOAT lb74[11];

    /* vector of size 5 */
    smoothMPC_FLOAT ub74[5];

    /* matrix of size [3 x 11] (column major format) */
    smoothMPC_FLOAT C74[33];

    /* vector of size 3 */
    smoothMPC_FLOAT e74[3];

    /* vector of size 11 */
    smoothMPC_FLOAT f75[11];

    /* vector of size 11 */
    smoothMPC_FLOAT lb75[11];

    /* vector of size 5 */
    smoothMPC_FLOAT ub75[5];

    /* matrix of size [3 x 11] (column major format) */
    smoothMPC_FLOAT C75[33];

    /* vector of size 3 */
    smoothMPC_FLOAT e75[3];

    /* vector of size 11 */
    smoothMPC_FLOAT f76[11];

    /* vector of size 11 */
    smoothMPC_FLOAT lb76[11];

    /* vector of size 5 */
    smoothMPC_FLOAT ub76[5];

    /* matrix of size [3 x 11] (column major format) */
    smoothMPC_FLOAT C76[33];

    /* vector of size 3 */
    smoothMPC_FLOAT e76[3];

    /* vector of size 11 */
    smoothMPC_FLOAT f77[11];

    /* vector of size 11 */
    smoothMPC_FLOAT lb77[11];

    /* vector of size 5 */
    smoothMPC_FLOAT ub77[5];

    /* matrix of size [3 x 11] (column major format) */
    smoothMPC_FLOAT C77[33];

    /* vector of size 3 */
    smoothMPC_FLOAT e77[3];

    /* vector of size 11 */
    smoothMPC_FLOAT f78[11];

    /* vector of size 11 */
    smoothMPC_FLOAT lb78[11];

    /* vector of size 5 */
    smoothMPC_FLOAT ub78[5];

    /* matrix of size [3 x 11] (column major format) */
    smoothMPC_FLOAT C78[33];

    /* vector of size 3 */
    smoothMPC_FLOAT e78[3];

    /* vector of size 11 */
    smoothMPC_FLOAT f79[11];

    /* vector of size 11 */
    smoothMPC_FLOAT lb79[11];

    /* vector of size 5 */
    smoothMPC_FLOAT ub79[5];

    /* matrix of size [3 x 11] (column major format) */
    smoothMPC_FLOAT C79[33];

    /* vector of size 3 */
    smoothMPC_FLOAT e79[3];

    /* vector of size 11 */
    smoothMPC_FLOAT f80[11];

    /* vector of size 11 */
    smoothMPC_FLOAT lb80[11];

    /* vector of size 5 */
    smoothMPC_FLOAT ub80[5];

    /* matrix of size [3 x 11] (column major format) */
    smoothMPC_FLOAT C80[33];

    /* vector of size 3 */
    smoothMPC_FLOAT e80[3];

    /* vector of size 11 */
    smoothMPC_FLOAT f81[11];

    /* vector of size 11 */
    smoothMPC_FLOAT lb81[11];

    /* vector of size 5 */
    smoothMPC_FLOAT ub81[5];

    /* matrix of size [3 x 11] (column major format) */
    smoothMPC_FLOAT C81[33];

    /* vector of size 3 */
    smoothMPC_FLOAT e81[3];

    /* vector of size 11 */
    smoothMPC_FLOAT f82[11];

    /* vector of size 11 */
    smoothMPC_FLOAT lb82[11];

    /* vector of size 5 */
    smoothMPC_FLOAT ub82[5];

    /* matrix of size [3 x 11] (column major format) */
    smoothMPC_FLOAT C82[33];

    /* vector of size 3 */
    smoothMPC_FLOAT e82[3];

    /* vector of size 11 */
    smoothMPC_FLOAT f83[11];

    /* vector of size 11 */
    smoothMPC_FLOAT lb83[11];

    /* vector of size 5 */
    smoothMPC_FLOAT ub83[5];

    /* matrix of size [3 x 11] (column major format) */
    smoothMPC_FLOAT C83[33];

    /* vector of size 3 */
    smoothMPC_FLOAT e83[3];

    /* vector of size 11 */
    smoothMPC_FLOAT f84[11];

    /* vector of size 11 */
    smoothMPC_FLOAT lb84[11];

    /* vector of size 5 */
    smoothMPC_FLOAT ub84[5];

    /* matrix of size [3 x 11] (column major format) */
    smoothMPC_FLOAT C84[33];

    /* vector of size 3 */
    smoothMPC_FLOAT e84[3];

    /* vector of size 11 */
    smoothMPC_FLOAT f85[11];

    /* vector of size 11 */
    smoothMPC_FLOAT lb85[11];

    /* vector of size 5 */
    smoothMPC_FLOAT ub85[5];

    /* matrix of size [3 x 11] (column major format) */
    smoothMPC_FLOAT C85[33];

    /* vector of size 3 */
    smoothMPC_FLOAT e85[3];

    /* vector of size 11 */
    smoothMPC_FLOAT f86[11];

    /* vector of size 11 */
    smoothMPC_FLOAT lb86[11];

    /* vector of size 5 */
    smoothMPC_FLOAT ub86[5];

    /* matrix of size [3 x 11] (column major format) */
    smoothMPC_FLOAT C86[33];

    /* vector of size 3 */
    smoothMPC_FLOAT e86[3];

    /* vector of size 11 */
    smoothMPC_FLOAT f87[11];

    /* vector of size 11 */
    smoothMPC_FLOAT lb87[11];

    /* vector of size 5 */
    smoothMPC_FLOAT ub87[5];

    /* matrix of size [3 x 11] (column major format) */
    smoothMPC_FLOAT C87[33];

    /* vector of size 3 */
    smoothMPC_FLOAT e87[3];

    /* vector of size 11 */
    smoothMPC_FLOAT f88[11];

    /* vector of size 11 */
    smoothMPC_FLOAT lb88[11];

    /* vector of size 5 */
    smoothMPC_FLOAT ub88[5];

    /* matrix of size [3 x 11] (column major format) */
    smoothMPC_FLOAT C88[33];

    /* vector of size 3 */
    smoothMPC_FLOAT e88[3];

    /* vector of size 11 */
    smoothMPC_FLOAT f89[11];

    /* vector of size 11 */
    smoothMPC_FLOAT lb89[11];

    /* vector of size 5 */
    smoothMPC_FLOAT ub89[5];

    /* matrix of size [3 x 11] (column major format) */
    smoothMPC_FLOAT C89[33];

    /* vector of size 3 */
    smoothMPC_FLOAT e89[3];

    /* vector of size 11 */
    smoothMPC_FLOAT f90[11];

    /* vector of size 11 */
    smoothMPC_FLOAT lb90[11];

    /* vector of size 5 */
    smoothMPC_FLOAT ub90[5];

    /* matrix of size [3 x 11] (column major format) */
    smoothMPC_FLOAT C90[33];

    /* vector of size 3 */
    smoothMPC_FLOAT e90[3];

    /* vector of size 11 */
    smoothMPC_FLOAT f91[11];

    /* vector of size 11 */
    smoothMPC_FLOAT lb91[11];

    /* vector of size 5 */
    smoothMPC_FLOAT ub91[5];

    /* matrix of size [3 x 11] (column major format) */
    smoothMPC_FLOAT C91[33];

    /* vector of size 3 */
    smoothMPC_FLOAT e91[3];

    /* vector of size 11 */
    smoothMPC_FLOAT f92[11];

    /* vector of size 11 */
    smoothMPC_FLOAT lb92[11];

    /* vector of size 5 */
    smoothMPC_FLOAT ub92[5];

    /* matrix of size [3 x 11] (column major format) */
    smoothMPC_FLOAT C92[33];

    /* vector of size 3 */
    smoothMPC_FLOAT e92[3];

    /* vector of size 11 */
    smoothMPC_FLOAT f93[11];

    /* vector of size 11 */
    smoothMPC_FLOAT lb93[11];

    /* vector of size 5 */
    smoothMPC_FLOAT ub93[5];

    /* matrix of size [3 x 11] (column major format) */
    smoothMPC_FLOAT C93[33];

    /* vector of size 3 */
    smoothMPC_FLOAT e93[3];

    /* vector of size 11 */
    smoothMPC_FLOAT f94[11];

    /* vector of size 11 */
    smoothMPC_FLOAT lb94[11];

    /* vector of size 5 */
    smoothMPC_FLOAT ub94[5];

    /* matrix of size [3 x 11] (column major format) */
    smoothMPC_FLOAT C94[33];

    /* vector of size 3 */
    smoothMPC_FLOAT e94[3];

    /* vector of size 11 */
    smoothMPC_FLOAT f95[11];

    /* vector of size 11 */
    smoothMPC_FLOAT lb95[11];

    /* vector of size 5 */
    smoothMPC_FLOAT ub95[5];

    /* matrix of size [3 x 11] (column major format) */
    smoothMPC_FLOAT C95[33];

    /* vector of size 3 */
    smoothMPC_FLOAT e95[3];

    /* vector of size 11 */
    smoothMPC_FLOAT f96[11];

    /* vector of size 11 */
    smoothMPC_FLOAT lb96[11];

    /* vector of size 5 */
    smoothMPC_FLOAT ub96[5];

    /* matrix of size [3 x 11] (column major format) */
    smoothMPC_FLOAT C96[33];

    /* vector of size 3 */
    smoothMPC_FLOAT e96[3];

    /* vector of size 11 */
    smoothMPC_FLOAT f97[11];

    /* vector of size 11 */
    smoothMPC_FLOAT lb97[11];

    /* vector of size 5 */
    smoothMPC_FLOAT ub97[5];

    /* matrix of size [3 x 11] (column major format) */
    smoothMPC_FLOAT C97[33];

    /* vector of size 3 */
    smoothMPC_FLOAT e97[3];

    /* vector of size 11 */
    smoothMPC_FLOAT f98[11];

    /* vector of size 11 */
    smoothMPC_FLOAT lb98[11];

    /* vector of size 5 */
    smoothMPC_FLOAT ub98[5];

    /* matrix of size [3 x 11] (column major format) */
    smoothMPC_FLOAT C98[33];

    /* vector of size 3 */
    smoothMPC_FLOAT e98[3];

    /* vector of size 11 */
    smoothMPC_FLOAT f99[11];

    /* vector of size 11 */
    smoothMPC_FLOAT lb99[11];

    /* vector of size 5 */
    smoothMPC_FLOAT ub99[5];

    /* matrix of size [3 x 11] (column major format) */
    smoothMPC_FLOAT C99[33];

    /* vector of size 3 */
    smoothMPC_FLOAT e99[3];

    /* vector of size 3 */
    smoothMPC_FLOAT f100[3];

    /* vector of size 3 */
    smoothMPC_FLOAT lb100[3];

    /* vector of size 3 */
    smoothMPC_FLOAT ub100[3];

    /* vector of size 3 */
    smoothMPC_FLOAT e100[3];

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

    /* vector of size 5 */
    smoothMPC_FLOAT z60[5];

    /* vector of size 5 */
    smoothMPC_FLOAT z61[5];

    /* vector of size 5 */
    smoothMPC_FLOAT z62[5];

    /* vector of size 5 */
    smoothMPC_FLOAT z63[5];

    /* vector of size 5 */
    smoothMPC_FLOAT z64[5];

    /* vector of size 5 */
    smoothMPC_FLOAT z65[5];

    /* vector of size 5 */
    smoothMPC_FLOAT z66[5];

    /* vector of size 5 */
    smoothMPC_FLOAT z67[5];

    /* vector of size 5 */
    smoothMPC_FLOAT z68[5];

    /* vector of size 5 */
    smoothMPC_FLOAT z69[5];

    /* vector of size 5 */
    smoothMPC_FLOAT z70[5];

    /* vector of size 5 */
    smoothMPC_FLOAT z71[5];

    /* vector of size 5 */
    smoothMPC_FLOAT z72[5];

    /* vector of size 5 */
    smoothMPC_FLOAT z73[5];

    /* vector of size 5 */
    smoothMPC_FLOAT z74[5];

    /* vector of size 5 */
    smoothMPC_FLOAT z75[5];

    /* vector of size 5 */
    smoothMPC_FLOAT z76[5];

    /* vector of size 5 */
    smoothMPC_FLOAT z77[5];

    /* vector of size 5 */
    smoothMPC_FLOAT z78[5];

    /* vector of size 5 */
    smoothMPC_FLOAT z79[5];

    /* vector of size 5 */
    smoothMPC_FLOAT z80[5];

    /* vector of size 5 */
    smoothMPC_FLOAT z81[5];

    /* vector of size 5 */
    smoothMPC_FLOAT z82[5];

    /* vector of size 5 */
    smoothMPC_FLOAT z83[5];

    /* vector of size 5 */
    smoothMPC_FLOAT z84[5];

    /* vector of size 5 */
    smoothMPC_FLOAT z85[5];

    /* vector of size 5 */
    smoothMPC_FLOAT z86[5];

    /* vector of size 5 */
    smoothMPC_FLOAT z87[5];

    /* vector of size 5 */
    smoothMPC_FLOAT z88[5];

    /* vector of size 5 */
    smoothMPC_FLOAT z89[5];

    /* vector of size 5 */
    smoothMPC_FLOAT z90[5];

    /* vector of size 5 */
    smoothMPC_FLOAT z91[5];

    /* vector of size 5 */
    smoothMPC_FLOAT z92[5];

    /* vector of size 5 */
    smoothMPC_FLOAT z93[5];

    /* vector of size 5 */
    smoothMPC_FLOAT z94[5];

    /* vector of size 5 */
    smoothMPC_FLOAT z95[5];

    /* vector of size 5 */
    smoothMPC_FLOAT z96[5];

    /* vector of size 5 */
    smoothMPC_FLOAT z97[5];

    /* vector of size 5 */
    smoothMPC_FLOAT z98[5];

    /* vector of size 5 */
    smoothMPC_FLOAT z99[5];

    /* vector of size 3 */
    smoothMPC_FLOAT z100[3];

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