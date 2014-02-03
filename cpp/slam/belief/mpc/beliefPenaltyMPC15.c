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

#include "beliefPenaltyMPC.h"

/* for square root */
#include <math.h> 

/* SAFE DIVISION ------------------------------------------------------- */
#define MAX(X,Y)  ((X) < (Y) ? (Y) : (X))
#define MIN(X,Y)  ((X) < (Y) ? (X) : (Y))
/*#define SAFEDIV_POS(X,Y)  ( (Y) < EPS ? ((X)/EPS) : (X)/(Y) ) 
#define EPS (1.0000E-013) */
#define BIGM (1E30)
#define BIGMM (1E60)

/* includes for parallel computation if necessary */


/* SYSTEM INCLUDES FOR PRINTING ---------------------------------------- */
#ifndef USEMEXPRINTS
#include <stdio.h>
#define PRINTTEXT printf
#else
#include "mex.h"
#define PRINTTEXT mexPrintf
#endif



/* LINEAR ALGEBRA LIBRARY ---------------------------------------------- */
/*
 * Initializes a vector of length 888 with a value.
 */
void beliefPenaltyMPC_LA_INITIALIZEVECTOR_888(beliefPenaltyMPC_FLOAT* vec, beliefPenaltyMPC_FLOAT value)
{
	int i;
	for( i=0; i<888; i++ )
	{
		vec[i] = value;
	}
}


/*
 * Initializes a vector of length 300 with a value.
 */
void beliefPenaltyMPC_LA_INITIALIZEVECTOR_300(beliefPenaltyMPC_FLOAT* vec, beliefPenaltyMPC_FLOAT value)
{
	int i;
	for( i=0; i<300; i++ )
	{
		vec[i] = value;
	}
}


/*
 * Initializes a vector of length 1216 with a value.
 */
void beliefPenaltyMPC_LA_INITIALIZEVECTOR_1216(beliefPenaltyMPC_FLOAT* vec, beliefPenaltyMPC_FLOAT value)
{
	int i;
	for( i=0; i<1216; i++ )
	{
		vec[i] = value;
	}
}


/* 
 * Calculates a dot product and adds it to a variable: z += x'*y; 
 * This function is for vectors of length 1216.
 */
void beliefPenaltyMPC_LA_DOTACC_1216(beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *y, beliefPenaltyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<1216; i++ ){
		*z += x[i]*y[i];
	}
}


/*
 * Calculates the gradient and the value for a quadratic function 0.5*z'*H*z + f'*z
 *
 * INPUTS:     H  - Symmetric Hessian, diag matrix of size [62 x 62]
 *             f  - column vector of size 62
 *             z  - column vector of size 62
 *
 * OUTPUTS: grad  - gradient at z (= H*z + f), column vector of size 62
 *          value <-- value + 0.5*z'*H*z + f'*z (value will be modified)
 */
void beliefPenaltyMPC_LA_DIAG_QUADFCN_62(beliefPenaltyMPC_FLOAT* H, beliefPenaltyMPC_FLOAT* f, beliefPenaltyMPC_FLOAT* z, beliefPenaltyMPC_FLOAT* grad, beliefPenaltyMPC_FLOAT* value)
{
	int i;
	beliefPenaltyMPC_FLOAT hz;	
	for( i=0; i<62; i++){
		hz = H[i]*z[i];
		grad[i] = hz + f[i];
		*value += 0.5*hz*z[i] + f[i]*z[i];
	}
}


/*
 * Calculates the gradient and the value for a quadratic function 0.5*z'*H*z + f'*z
 *
 * INPUTS:     H  - Symmetric Hessian, diag matrix of size [20 x 20]
 *             f  - column vector of size 20
 *             z  - column vector of size 20
 *
 * OUTPUTS: grad  - gradient at z (= H*z + f), column vector of size 20
 *          value <-- value + 0.5*z'*H*z + f'*z (value will be modified)
 */
void beliefPenaltyMPC_LA_DIAG_QUADFCN_20(beliefPenaltyMPC_FLOAT* H, beliefPenaltyMPC_FLOAT* f, beliefPenaltyMPC_FLOAT* z, beliefPenaltyMPC_FLOAT* grad, beliefPenaltyMPC_FLOAT* value)
{
	int i;
	beliefPenaltyMPC_FLOAT hz;	
	for( i=0; i<20; i++){
		hz = H[i]*z[i];
		grad[i] = hz + f[i];
		*value += 0.5*hz*z[i] + f[i]*z[i];
	}
}


/* 
 * Computes r = B*u - b
 * and      y = max([norm(r,inf), y])
 * and      z -= l'*r
 * where B is stored in diabzero format
 */
void beliefPenaltyMPC_LA_DIAGZERO_MVMSUB6_20(beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *b, beliefPenaltyMPC_FLOAT *l, beliefPenaltyMPC_FLOAT *r, beliefPenaltyMPC_FLOAT *z, beliefPenaltyMPC_FLOAT *y)
{
	int i;
	beliefPenaltyMPC_FLOAT Bu[20];
	beliefPenaltyMPC_FLOAT norm = *y;
	beliefPenaltyMPC_FLOAT lr = 0;

	/* do A*x + B*u first */
	for( i=0; i<20; i++ ){
		Bu[i] = B[i]*u[i];
	}	

	for( i=0; i<20; i++ ){
		r[i] = Bu[i] - b[i];
		lr += l[i]*r[i];
		if( r[i] > norm ){
			norm = r[i];
		}
		if( -r[i] > norm ){
			norm = -r[i];
		}
	}
	*y = norm;
	*z -= lr;
}


/* 
 * Computes r = A*x + B*u - b
 * and      y = max([norm(r,inf), y])
 * and      z -= l'*r
 * where A is stored in column major format
 */
void beliefPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_20_62_62(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *b, beliefPenaltyMPC_FLOAT *l, beliefPenaltyMPC_FLOAT *r, beliefPenaltyMPC_FLOAT *z, beliefPenaltyMPC_FLOAT *y)
{
	int i;
	int j;
	int k = 0;
	beliefPenaltyMPC_FLOAT AxBu[20];
	beliefPenaltyMPC_FLOAT norm = *y;
	beliefPenaltyMPC_FLOAT lr = 0;

	/* do A*x + B*u first */
	for( i=0; i<20; i++ ){
		AxBu[i] = A[k++]*x[0] + B[i]*u[i];
	}	

	for( j=1; j<62; j++ ){		
		for( i=0; i<20; i++ ){
			AxBu[i] += A[k++]*x[j];
		}
	}

	for( i=0; i<20; i++ ){
		r[i] = AxBu[i] - b[i];
		lr += l[i]*r[i];
		if( r[i] > norm ){
			norm = r[i];
		}
		if( -r[i] > norm ){
			norm = -r[i];
		}
	}
	*y = norm;
	*z -= lr;
}


/* 
 * Computes r = A*x + B*u - b
 * and      y = max([norm(r,inf), y])
 * and      z -= l'*r
 * where A is stored in column major format
 */
void beliefPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_20_62_20(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *b, beliefPenaltyMPC_FLOAT *l, beliefPenaltyMPC_FLOAT *r, beliefPenaltyMPC_FLOAT *z, beliefPenaltyMPC_FLOAT *y)
{
	int i;
	int j;
	int k = 0;
	beliefPenaltyMPC_FLOAT AxBu[20];
	beliefPenaltyMPC_FLOAT norm = *y;
	beliefPenaltyMPC_FLOAT lr = 0;

	/* do A*x + B*u first */
	for( i=0; i<20; i++ ){
		AxBu[i] = A[k++]*x[0] + B[i]*u[i];
	}	

	for( j=1; j<62; j++ ){		
		for( i=0; i<20; i++ ){
			AxBu[i] += A[k++]*x[j];
		}
	}

	for( i=0; i<20; i++ ){
		r[i] = AxBu[i] - b[i];
		lr += l[i]*r[i];
		if( r[i] > norm ){
			norm = r[i];
		}
		if( -r[i] > norm ){
			norm = -r[i];
		}
	}
	*y = norm;
	*z -= lr;
}


/*
 * Matrix vector multiplication z = A'*x + B'*y 
 * where A is of size [20 x 62] and stored in column major format.
 * and B is of size [20 x 62] and stored in diagzero format
 * Note the transposes of A and B!
 */
void beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_20_62_20(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *y, beliefPenaltyMPC_FLOAT *z)
{
	int i;
	int j;
	int k = 0;
	for( i=0; i<20; i++ ){
		z[i] = 0;
		for( j=0; j<20; j++ ){
			z[i] += A[k++]*x[j];
		}
		z[i] += B[i]*y[i];
	}
	for( i=20 ;i<62; i++ ){
		z[i] = 0;
		for( j=0; j<20; j++ ){
			z[i] += A[k++]*x[j];
		}
	}
}


/*
 * Matrix vector multiplication y = M'*x where M is of size [20 x 20]
 * and stored in diagzero format. Note the transpose of M!
 */
void beliefPenaltyMPC_LA_DIAGZERO_MTVM_20_20(beliefPenaltyMPC_FLOAT *M, beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *y)
{
	int i;
	for( i=0; i<20; i++ ){
		y[i] = M[i]*x[i];
	}
}


/*
 * Vector subtraction and addition.
 *	 Input: five vectors t, tidx, u, v, w and two scalars z and r
 *	 Output: y = t(tidx) - u + w
 *           z = z - v'*x;
 *           r = max([norm(y,inf), z]);
 * for vectors of length 62. Output z is of course scalar.
 */
void beliefPenaltyMPC_LA_VSUBADD3_62(beliefPenaltyMPC_FLOAT* t, beliefPenaltyMPC_FLOAT* u, int* uidx, beliefPenaltyMPC_FLOAT* v, beliefPenaltyMPC_FLOAT* w, beliefPenaltyMPC_FLOAT* y, beliefPenaltyMPC_FLOAT* z, beliefPenaltyMPC_FLOAT* r)
{
	int i;
	beliefPenaltyMPC_FLOAT norm = *r;
	beliefPenaltyMPC_FLOAT vx = 0;
	beliefPenaltyMPC_FLOAT x;
	for( i=0; i<62; i++){
		x = t[i] - u[uidx[i]];
		y[i] = x + w[i];
		vx += v[i]*x;
		if( y[i] > norm ){
			norm = y[i];
		}
		if( -y[i] > norm ){
			norm = -y[i];
		}
	}
	*z -= vx;
	*r = norm;
}


/*
 * Vector subtraction and addition.
 *	 Input: five vectors t, tidx, u, v, w and two scalars z and r
 *	 Output: y = t(tidx) - u + w
 *           z = z - v'*x;
 *           r = max([norm(y,inf), z]);
 * for vectors of length 22. Output z is of course scalar.
 */
void beliefPenaltyMPC_LA_VSUBADD2_22(beliefPenaltyMPC_FLOAT* t, int* tidx, beliefPenaltyMPC_FLOAT* u, beliefPenaltyMPC_FLOAT* v, beliefPenaltyMPC_FLOAT* w, beliefPenaltyMPC_FLOAT* y, beliefPenaltyMPC_FLOAT* z, beliefPenaltyMPC_FLOAT* r)
{
	int i;
	beliefPenaltyMPC_FLOAT norm = *r;
	beliefPenaltyMPC_FLOAT vx = 0;
	beliefPenaltyMPC_FLOAT x;
	for( i=0; i<22; i++){
		x = t[tidx[i]] - u[i];
		y[i] = x + w[i];
		vx += v[i]*x;
		if( y[i] > norm ){
			norm = y[i];
		}
		if( -y[i] > norm ){
			norm = -y[i];
		}
	}
	*z -= vx;
	*r = norm;
}


/*
 * Vector subtraction and addition.
 *	 Input: five vectors t, tidx, u, v, w and two scalars z and r
 *	 Output: y = t(tidx) - u + w
 *           z = z - v'*x;
 *           r = max([norm(y,inf), z]);
 * for vectors of length 20. Output z is of course scalar.
 */
void beliefPenaltyMPC_LA_VSUBADD3_20(beliefPenaltyMPC_FLOAT* t, beliefPenaltyMPC_FLOAT* u, int* uidx, beliefPenaltyMPC_FLOAT* v, beliefPenaltyMPC_FLOAT* w, beliefPenaltyMPC_FLOAT* y, beliefPenaltyMPC_FLOAT* z, beliefPenaltyMPC_FLOAT* r)
{
	int i;
	beliefPenaltyMPC_FLOAT norm = *r;
	beliefPenaltyMPC_FLOAT vx = 0;
	beliefPenaltyMPC_FLOAT x;
	for( i=0; i<20; i++){
		x = t[i] - u[uidx[i]];
		y[i] = x + w[i];
		vx += v[i]*x;
		if( y[i] > norm ){
			norm = y[i];
		}
		if( -y[i] > norm ){
			norm = -y[i];
		}
	}
	*z -= vx;
	*r = norm;
}


/*
 * Vector subtraction and addition.
 *	 Input: five vectors t, tidx, u, v, w and two scalars z and r
 *	 Output: y = t(tidx) - u + w
 *           z = z - v'*x;
 *           r = max([norm(y,inf), z]);
 * for vectors of length 20. Output z is of course scalar.
 */
void beliefPenaltyMPC_LA_VSUBADD2_20(beliefPenaltyMPC_FLOAT* t, int* tidx, beliefPenaltyMPC_FLOAT* u, beliefPenaltyMPC_FLOAT* v, beliefPenaltyMPC_FLOAT* w, beliefPenaltyMPC_FLOAT* y, beliefPenaltyMPC_FLOAT* z, beliefPenaltyMPC_FLOAT* r)
{
	int i;
	beliefPenaltyMPC_FLOAT norm = *r;
	beliefPenaltyMPC_FLOAT vx = 0;
	beliefPenaltyMPC_FLOAT x;
	for( i=0; i<20; i++){
		x = t[tidx[i]] - u[i];
		y[i] = x + w[i];
		vx += v[i]*x;
		if( y[i] > norm ){
			norm = y[i];
		}
		if( -y[i] > norm ){
			norm = -y[i];
		}
	}
	*z -= vx;
	*r = norm;
}


/*
 * Computes inequality constraints gradient-
 * Special function for box constraints of length 62
 * Returns also L/S, a value that is often used elsewhere.
 */
void beliefPenaltyMPC_LA_INEQ_B_GRAD_62_62_22(beliefPenaltyMPC_FLOAT *lu, beliefPenaltyMPC_FLOAT *su, beliefPenaltyMPC_FLOAT *ru, beliefPenaltyMPC_FLOAT *ll, beliefPenaltyMPC_FLOAT *sl, beliefPenaltyMPC_FLOAT *rl, int* lbIdx, int* ubIdx, beliefPenaltyMPC_FLOAT *grad, beliefPenaltyMPC_FLOAT *lubysu, beliefPenaltyMPC_FLOAT *llbysl)
{
	int i;
	for( i=0; i<62; i++ ){
		grad[i] = 0;
	}
	for( i=0; i<62; i++ ){		
		llbysl[i] = ll[i] / sl[i];
		grad[lbIdx[i]] -= llbysl[i]*rl[i];
	}
	for( i=0; i<22; i++ ){
		lubysu[i] = lu[i] / su[i];
		grad[ubIdx[i]] += lubysu[i]*ru[i];
	}
}


/*
 * Computes inequality constraints gradient-
 * Special function for box constraints of length 20
 * Returns also L/S, a value that is often used elsewhere.
 */
void beliefPenaltyMPC_LA_INEQ_B_GRAD_20_20_20(beliefPenaltyMPC_FLOAT *lu, beliefPenaltyMPC_FLOAT *su, beliefPenaltyMPC_FLOAT *ru, beliefPenaltyMPC_FLOAT *ll, beliefPenaltyMPC_FLOAT *sl, beliefPenaltyMPC_FLOAT *rl, int* lbIdx, int* ubIdx, beliefPenaltyMPC_FLOAT *grad, beliefPenaltyMPC_FLOAT *lubysu, beliefPenaltyMPC_FLOAT *llbysl)
{
	int i;
	for( i=0; i<20; i++ ){
		grad[i] = 0;
	}
	for( i=0; i<20; i++ ){		
		llbysl[i] = ll[i] / sl[i];
		grad[lbIdx[i]] -= llbysl[i]*rl[i];
	}
	for( i=0; i<20; i++ ){
		lubysu[i] = lu[i] / su[i];
		grad[ubIdx[i]] += lubysu[i]*ru[i];
	}
}


/*
 * Addition of three vectors  z = u + w + v
 * of length 888.
 */
void beliefPenaltyMPC_LA_VVADD3_888(beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *v, beliefPenaltyMPC_FLOAT *w, beliefPenaltyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<888; i++ ){
		z[i] = u[i] + v[i] + w[i];
	}
}


/*
 * Special function to compute the diagonal cholesky factorization of the 
 * positive definite augmented Hessian for block size 62.
 *
 * Inputs: - H = diagonal cost Hessian in diagonal storage format
 *         - llbysl = L / S of lower bounds
 *         - lubysu = L / S of upper bounds
 *
 * Output: Phi = sqrt(H + diag(llbysl) + diag(lubysu))
 * where Phi is stored in diagonal storage format
 */
void beliefPenaltyMPC_LA_DIAG_CHOL_LBUB_62_62_22(beliefPenaltyMPC_FLOAT *H, beliefPenaltyMPC_FLOAT *llbysl, int* lbIdx, beliefPenaltyMPC_FLOAT *lubysu, int* ubIdx, beliefPenaltyMPC_FLOAT *Phi)


{
	int i;
	
	/* copy  H into PHI */
	for( i=0; i<62; i++ ){
		Phi[i] = H[i];
	}

	/* add llbysl onto Phi where necessary */
	for( i=0; i<62; i++ ){
		Phi[lbIdx[i]] += llbysl[i];
	}

	/* add lubysu onto Phi where necessary */
	for( i=0; i<22; i++){
		Phi[ubIdx[i]] +=  lubysu[i];
	}
	
	/* compute cholesky */
	for(i=0; i<62; i++)
	{
#if beliefPenaltyMPC_SET_PRINTLEVEL > 0 && defined PRINTNUMERICALWARNINGS
		if( Phi[i] < 1.0000000000000000E-013 )
		{
            PRINTTEXT("WARNING: small pivot in Cholesky fact. (=%3.1e < eps=%3.1e), regularizing to %3.1e\n",Phi[i],1.0000000000000000E-013,4.0000000000000002E-004);
			Phi[i] = 2.0000000000000000E-002;
		}
		else
		{
			Phi[i] = sqrt(Phi[i]);
		}
#else
		Phi[i] = Phi[i] < 1.0000000000000000E-013 ? 2.0000000000000000E-002 : sqrt(Phi[i]);
#endif
	}

}


/**
 * Forward substitution for the matrix equation A*L' = B
 * where A is to be computed and is of size [20 x 62],
 * B is given and of size [20 x 62], L is a diagonal
 * matrix of size 20 stored in diagonal matrix 
 * storage format. Note the transpose of L has no impact!
 *
 * Result: A in column major storage format.
 *
 */
void beliefPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_20_62(beliefPenaltyMPC_FLOAT *L, beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *A)
{
    int i,j;
	 int k = 0;

	for( j=0; j<62; j++){
		for( i=0; i<20; i++){
			A[k] = B[k]/L[j];
			k++;
		}
	}

}


/**
 * Forward substitution for the matrix equation A*L' = B
 * where A is to be computed and is of size [20 x 62],
 * B is given and of size [20 x 62], L is a diagonal
 *  matrix of size 62 stored in diagonal 
 * storage format. Note the transpose of L!
 *
 * Result: A in diagonalzero storage format.
 *
 */
void beliefPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_20_62(beliefPenaltyMPC_FLOAT *L, beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *A)
{
	int j;
    for( j=0; j<62; j++ ){   
		A[j] = B[j]/L[j];
     }
}


/**
 * Compute C = A*B' where 
 *
 *	size(A) = [20 x 62]
 *  size(B) = [20 x 62] in diagzero format
 * 
 * A and C matrices are stored in column major format.
 * 
 * 
 */
void beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMTM_20_62_20(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *C)
{
    int i, j;
	
	for( i=0; i<20; i++ ){
		for( j=0; j<20; j++){
			C[j*20+i] = B[i*20+j]*A[i];
		}
	}

}


/**
 * Forward substitution to solve L*y = b where L is a
 * diagonal matrix in vector storage format.
 * 
 * The dimensions involved are 62.
 */
void beliefPenaltyMPC_LA_DIAG_FORWARDSUB_62(beliefPenaltyMPC_FLOAT *L, beliefPenaltyMPC_FLOAT *b, beliefPenaltyMPC_FLOAT *y)
{
    int i;

    for( i=0; i<62; i++ ){
		y[i] = b[i]/L[i];
    }
}


/*
 * Special function to compute the diagonal cholesky factorization of the 
 * positive definite augmented Hessian for block size 20.
 *
 * Inputs: - H = diagonal cost Hessian in diagonal storage format
 *         - llbysl = L / S of lower bounds
 *         - lubysu = L / S of upper bounds
 *
 * Output: Phi = sqrt(H + diag(llbysl) + diag(lubysu))
 * where Phi is stored in diagonal storage format
 */
void beliefPenaltyMPC_LA_DIAG_CHOL_ONELOOP_LBUB_20_20_20(beliefPenaltyMPC_FLOAT *H, beliefPenaltyMPC_FLOAT *llbysl, int* lbIdx, beliefPenaltyMPC_FLOAT *lubysu, int* ubIdx, beliefPenaltyMPC_FLOAT *Phi)


{
	int i;
	
	/* compute cholesky */
	for( i=0; i<20; i++ ){
		Phi[i] = H[i] + llbysl[i] + lubysu[i];

#if beliefPenaltyMPC_SET_PRINTLEVEL > 0 && defined PRINTNUMERICALWARNINGS
		if( Phi[i] < 1.0000000000000000E-013 )
		{
            PRINTTEXT("WARNING: small pivot in Cholesky fact. (=%3.1e < eps=%3.1e), regularizing to %3.1e\n",Phi[i],1.0000000000000000E-013,4.0000000000000002E-004);
			Phi[i] = 2.0000000000000000E-002;
		}
		else
		{
			Phi[i] = sqrt(Phi[i]);
		}
#else
		Phi[i] = Phi[i] < 1.0000000000000000E-013 ? 2.0000000000000000E-002 : sqrt(Phi[i]);
#endif
	}
	
}


/**
 * Forward substitution for the matrix equation A*L' = B
 * where A is to be computed and is of size [20 x 20],
 * B is given and of size [20 x 20], L is a diagonal
 *  matrix of size 20 stored in diagonal 
 * storage format. Note the transpose of L!
 *
 * Result: A in diagonalzero storage format.
 *
 */
void beliefPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_20_20(beliefPenaltyMPC_FLOAT *L, beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *A)
{
	int j;
    for( j=0; j<20; j++ ){   
		A[j] = B[j]/L[j];
     }
}


/**
 * Forward substitution to solve L*y = b where L is a
 * diagonal matrix in vector storage format.
 * 
 * The dimensions involved are 20.
 */
void beliefPenaltyMPC_LA_DIAG_FORWARDSUB_20(beliefPenaltyMPC_FLOAT *L, beliefPenaltyMPC_FLOAT *b, beliefPenaltyMPC_FLOAT *y)
{
    int i;

    for( i=0; i<20; i++ ){
		y[i] = b[i]/L[i];
    }
}


/**
 * Compute L = B*B', where L is lower triangular of size NXp1
 * and B is a diagzero matrix of size [20 x 62] in column
 * storage format.
 * 
 */
void beliefPenaltyMPC_LA_DIAGZERO_MMT_20(beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *L)
{
    int i, ii, di;
    
    ii = 0; di = 0;
    for( i=0; i<20; i++ ){        
		L[ii+i] = B[i]*B[i];
        ii += ++di;
    }
}


/* 
 * Computes r = b - B*u
 * B is stored in diagzero format
 */
void beliefPenaltyMPC_LA_DIAGZERO_MVMSUB7_20(beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *b, beliefPenaltyMPC_FLOAT *r)
{
	int i;

	for( i=0; i<20; i++ ){
		r[i] = b[i] - B[i]*u[i];
	}	
	
}


/**
 * Compute L = A*A' + B*B', where L is lower triangular of size NXp1
 * and A is a dense matrix of size [20 x 62] in column
 * storage format, and B is of size [20 x 62] diagonalzero
 * storage format.
 * 
 * THIS ONE HAS THE WORST ACCES PATTERN POSSIBLE. 
 * POSSIBKE FIX: PUT A AND B INTO ROW MAJOR FORMAT FIRST.
 * 
 */
void beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMT2_20_62_62(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *L)
{
    int i, j, k, ii, di;
    beliefPenaltyMPC_FLOAT ltemp;
    
    ii = 0; di = 0;
    for( i=0; i<20; i++ ){        
        for( j=0; j<=i; j++ ){
            ltemp = 0; 
            for( k=0; k<62; k++ ){
                ltemp += A[k*20+i]*A[k*20+j];
            }		
            L[ii+j] = ltemp;
        }
		/* work on the diagonal
		 * there might be i == j, but j has already been incremented so it is i == j-1 */
		L[ii+i] += B[i]*B[i];
        ii += ++di;
    }
}


/* 
 * Computes r = b - A*x - B*u
 * where A is stored in column major format
 * and B is stored in diagzero format
 */
void beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_20_62_62(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *b, beliefPenaltyMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;

	for( i=0; i<20; i++ ){
		r[i] = b[i] - A[k++]*x[0] - B[i]*u[i];
	}	

	for( j=1; j<62; j++ ){		
		for( i=0; i<20; i++ ){
			r[i] -= A[k++]*x[j];
		}
	}
	
}


/**
 * Compute L = A*A' + B*B', where L is lower triangular of size NXp1
 * and A is a dense matrix of size [20 x 62] in column
 * storage format, and B is of size [20 x 20] diagonalzero
 * storage format.
 * 
 * THIS ONE HAS THE WORST ACCES PATTERN POSSIBLE. 
 * POSSIBKE FIX: PUT A AND B INTO ROW MAJOR FORMAT FIRST.
 * 
 */
void beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMT2_20_62_20(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *L)
{
    int i, j, k, ii, di;
    beliefPenaltyMPC_FLOAT ltemp;
    
    ii = 0; di = 0;
    for( i=0; i<20; i++ ){        
        for( j=0; j<=i; j++ ){
            ltemp = 0; 
            for( k=0; k<62; k++ ){
                ltemp += A[k*20+i]*A[k*20+j];
            }		
            L[ii+j] = ltemp;
        }
		/* work on the diagonal
		 * there might be i == j, but j has already been incremented so it is i == j-1 */
		L[ii+i] += B[i]*B[i];
        ii += ++di;
    }
}


/* 
 * Computes r = b - A*x - B*u
 * where A is stored in column major format
 * and B is stored in diagzero format
 */
void beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_20_62_20(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *b, beliefPenaltyMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;

	for( i=0; i<20; i++ ){
		r[i] = b[i] - A[k++]*x[0] - B[i]*u[i];
	}	

	for( j=1; j<62; j++ ){		
		for( i=0; i<20; i++ ){
			r[i] -= A[k++]*x[j];
		}
	}
	
}


/**
 * Cholesky factorization as above, but working on a matrix in 
 * lower triangular storage format of size 20 and outputting
 * the Cholesky factor to matrix L in lower triangular format.
 */
void beliefPenaltyMPC_LA_DENSE_CHOL_20(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *L)
{
    int i, j, k, di, dj;
	 int ii, jj;

    beliefPenaltyMPC_FLOAT l;
    beliefPenaltyMPC_FLOAT Mii;

	/* copy A to L first and then operate on L */
	/* COULD BE OPTIMIZED */
	ii=0; di=0;
	for( i=0; i<20; i++ ){
		for( j=0; j<=i; j++ ){
			L[ii+j] = A[ii+j];
		}
		ii += ++di;
	}    
	
	/* factor L */
	ii=0; di=0;
    for( i=0; i<20; i++ ){
        l = 0;
        for( k=0; k<i; k++ ){
            l += L[ii+k]*L[ii+k];
        }        
        
        Mii = L[ii+i] - l;
        
#if beliefPenaltyMPC_SET_PRINTLEVEL > 0 && defined PRINTNUMERICALWARNINGS
        if( Mii < 1.0000000000000000E-013 ){
             PRINTTEXT("WARNING (CHOL): small %d-th pivot in Cholesky fact. (=%3.1e < eps=%3.1e), regularizing to %3.1e\n",i,Mii,1.0000000000000000E-013,4.0000000000000002E-004);
			 L[ii+i] = 2.0000000000000000E-002;
		} else
		{
			L[ii+i] = sqrt(Mii);
		}
#else
		L[ii+i] = Mii < 1.0000000000000000E-013 ? 2.0000000000000000E-002 : sqrt(Mii);
#endif

		jj = ((i+1)*(i+2))/2; dj = i+1;
        for( j=i+1; j<20; j++ ){
            l = 0;            
            for( k=0; k<i; k++ ){
                l += L[jj+k]*L[ii+k];
            }

			/* saturate values for numerical stability */
			l = MIN(l,  BIGMM);
			l = MAX(l, -BIGMM);

            L[jj+i] = (L[jj+i] - l)/L[ii+i];            
			jj += ++dj;
        }
		ii += ++di;
    }	
}


/**
 * Forward substitution to solve L*y = b where L is a
 * lower triangular matrix in triangular storage format.
 * 
 * The dimensions involved are 20.
 */
void beliefPenaltyMPC_LA_DENSE_FORWARDSUB_20(beliefPenaltyMPC_FLOAT *L, beliefPenaltyMPC_FLOAT *b, beliefPenaltyMPC_FLOAT *y)
{
    int i,j,ii,di;
    beliefPenaltyMPC_FLOAT yel;
            
    ii = 0; di = 0;
    for( i=0; i<20; i++ ){
        yel = b[i];        
        for( j=0; j<i; j++ ){
            yel -= y[j]*L[ii+j];
        }

		/* saturate for numerical stability  */
		yel = MIN(yel, BIGM);
		yel = MAX(yel, -BIGM);

        y[i] = yel / L[ii+i];
        ii += ++di;
    }
}


/** 
 * Forward substitution for the matrix equation A*L' = B'
 * where A is to be computed and is of size [20 x 20],
 * B is given and of size [20 x 20], L is a lower tri-
 * angular matrix of size 20 stored in lower triangular 
 * storage format. Note the transpose of L AND B!
 *
 * Result: A in column major storage format.
 *
 */
void beliefPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_20_20(beliefPenaltyMPC_FLOAT *L, beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *A)
{
    int i,j,k,ii,di;
    beliefPenaltyMPC_FLOAT a;
    
    ii=0; di=0;
    for( j=0; j<20; j++ ){        
        for( i=0; i<20; i++ ){
            a = B[i*20+j];
            for( k=0; k<j; k++ ){
                a -= A[k*20+i]*L[ii+k];
            }    

			/* saturate for numerical stability */
			a = MIN(a, BIGM);
			a = MAX(a, -BIGM); 

			A[j*20+i] = a/L[ii+j];			
        }
        ii += ++di;
    }
}


/**
 * Compute L = L - A*A', where L is lower triangular of size 20
 * and A is a dense matrix of size [20 x 20] in column
 * storage format.
 * 
 * THIS ONE HAS THE WORST ACCES PATTERN POSSIBLE. 
 * POSSIBKE FIX: PUT A INTO ROW MAJOR FORMAT FIRST.
 * 
 */
void beliefPenaltyMPC_LA_DENSE_MMTSUB_20_20(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *L)
{
    int i, j, k, ii, di;
    beliefPenaltyMPC_FLOAT ltemp;
    
    ii = 0; di = 0;
    for( i=0; i<20; i++ ){        
        for( j=0; j<=i; j++ ){
            ltemp = 0; 
            for( k=0; k<20; k++ ){
                ltemp += A[k*20+i]*A[k*20+j];
            }						
            L[ii+j] -= ltemp;
        }
        ii += ++di;
    }
}


/* 
 * Computes r = b - A*x
 * where A is stored in column major format
 */
void beliefPenaltyMPC_LA_DENSE_MVMSUB1_20_20(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *b, beliefPenaltyMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;

	for( i=0; i<20; i++ ){
		r[i] = b[i] - A[k++]*x[0];
	}	
	for( j=1; j<20; j++ ){		
		for( i=0; i<20; i++ ){
			r[i] -= A[k++]*x[j];
		}
	}
}


/**
 * Backward Substitution to solve L^T*x = y where L is a
 * lower triangular matrix in triangular storage format.
 * 
 * All involved dimensions are 20.
 */
void beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_20(beliefPenaltyMPC_FLOAT *L, beliefPenaltyMPC_FLOAT *y, beliefPenaltyMPC_FLOAT *x)
{
    int i, ii, di, j, jj, dj;
    beliefPenaltyMPC_FLOAT xel;    
	int start = 190;
    
    /* now solve L^T*x = y by backward substitution */
    ii = start; di = 19;
    for( i=19; i>=0; i-- ){        
        xel = y[i];        
        jj = start; dj = 19;
        for( j=19; j>i; j-- ){
            xel -= x[j]*L[jj+i];
            jj -= dj--;
        }

		/* saturate for numerical stability */
		xel = MIN(xel, BIGM);
		xel = MAX(xel, -BIGM); 

        x[i] = xel / L[ii+i];
        ii -= di--;
    }
}


/*
 * Matrix vector multiplication y = b - M'*x where M is of size [20 x 20]
 * and stored in column major format. Note the transpose of M!
 */
void beliefPenaltyMPC_LA_DENSE_MTVMSUB_20_20(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *b, beliefPenaltyMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0; 
	for( i=0; i<20; i++ ){
		r[i] = b[i];
		for( j=0; j<20; j++ ){
			r[i] -= A[k++]*x[j];
		}
	}
}


/*
 * Vector subtraction z = -x - y for vectors of length 888.
 */
void beliefPenaltyMPC_LA_VSUB2_888(beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *y, beliefPenaltyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<888; i++){
		z[i] = -x[i] - y[i];
	}
}


/**
 * Forward-Backward-Substitution to solve L*L^T*x = b where L is a
 * diagonal matrix of size 62 in vector
 * storage format.
 */
void beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_62(beliefPenaltyMPC_FLOAT *L, beliefPenaltyMPC_FLOAT *b, beliefPenaltyMPC_FLOAT *x)
{
    int i;
            
    /* solve Ly = b by forward and backward substitution */
    for( i=0; i<62; i++ ){
		x[i] = b[i]/(L[i]*L[i]);
    }
    
}


/**
 * Forward-Backward-Substitution to solve L*L^T*x = b where L is a
 * diagonal matrix of size 20 in vector
 * storage format.
 */
void beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_20(beliefPenaltyMPC_FLOAT *L, beliefPenaltyMPC_FLOAT *b, beliefPenaltyMPC_FLOAT *x)
{
    int i;
            
    /* solve Ly = b by forward and backward substitution */
    for( i=0; i<20; i++ ){
		x[i] = b[i]/(L[i]*L[i]);
    }
    
}


/*
 * Vector subtraction z = x(xidx) - y where y, z and xidx are of length 62,
 * and x has length 62 and is indexed through yidx.
 */
void beliefPenaltyMPC_LA_VSUB_INDEXED_62(beliefPenaltyMPC_FLOAT *x, int* xidx, beliefPenaltyMPC_FLOAT *y, beliefPenaltyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<62; i++){
		z[i] = x[xidx[i]] - y[i];
	}
}


/*
 * Vector subtraction x = -u.*v - w for vectors of length 62.
 */
void beliefPenaltyMPC_LA_VSUB3_62(beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *v, beliefPenaltyMPC_FLOAT *w, beliefPenaltyMPC_FLOAT *x)
{
	int i;
	for( i=0; i<62; i++){
		x[i] = -u[i]*v[i] - w[i];
	}
}


/*
 * Vector subtraction z = -x - y(yidx) where y is of length 62
 * and z, x and yidx are of length 22.
 */
void beliefPenaltyMPC_LA_VSUB2_INDEXED_22(beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *y, int* yidx, beliefPenaltyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<22; i++){
		z[i] = -x[i] - y[yidx[i]];
	}
}


/*
 * Vector subtraction x = -u.*v - w for vectors of length 22.
 */
void beliefPenaltyMPC_LA_VSUB3_22(beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *v, beliefPenaltyMPC_FLOAT *w, beliefPenaltyMPC_FLOAT *x)
{
	int i;
	for( i=0; i<22; i++){
		x[i] = -u[i]*v[i] - w[i];
	}
}


/*
 * Vector subtraction z = x(xidx) - y where y, z and xidx are of length 20,
 * and x has length 20 and is indexed through yidx.
 */
void beliefPenaltyMPC_LA_VSUB_INDEXED_20(beliefPenaltyMPC_FLOAT *x, int* xidx, beliefPenaltyMPC_FLOAT *y, beliefPenaltyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<20; i++){
		z[i] = x[xidx[i]] - y[i];
	}
}


/*
 * Vector subtraction x = -u.*v - w for vectors of length 20.
 */
void beliefPenaltyMPC_LA_VSUB3_20(beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *v, beliefPenaltyMPC_FLOAT *w, beliefPenaltyMPC_FLOAT *x)
{
	int i;
	for( i=0; i<20; i++){
		x[i] = -u[i]*v[i] - w[i];
	}
}


/*
 * Vector subtraction z = -x - y(yidx) where y is of length 20
 * and z, x and yidx are of length 20.
 */
void beliefPenaltyMPC_LA_VSUB2_INDEXED_20(beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *y, int* yidx, beliefPenaltyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<20; i++){
		z[i] = -x[i] - y[yidx[i]];
	}
}


/**
 * Backtracking line search.
 * 
 * First determine the maximum line length by a feasibility line
 * search, i.e. a ~= argmax{ a \in [0...1] s.t. l+a*dl >= 0 and s+a*ds >= 0}.
 *
 * The function returns either the number of iterations or exits the error code
 * beliefPenaltyMPC_NOPROGRESS (should be negative).
 */
int beliefPenaltyMPC_LINESEARCH_BACKTRACKING_AFFINE(beliefPenaltyMPC_FLOAT *l, beliefPenaltyMPC_FLOAT *s, beliefPenaltyMPC_FLOAT *dl, beliefPenaltyMPC_FLOAT *ds, beliefPenaltyMPC_FLOAT *a, beliefPenaltyMPC_FLOAT *mu_aff)
{
    int i;
	int lsIt=1;    
    beliefPenaltyMPC_FLOAT dltemp;
    beliefPenaltyMPC_FLOAT dstemp;
    beliefPenaltyMPC_FLOAT mya = 1.0;
    beliefPenaltyMPC_FLOAT mymu;
        
    while( 1 ){                        

        /* 
         * Compute both snew and wnew together.
         * We compute also mu_affine along the way here, as the
         * values might be in registers, so it should be cheaper.
         */
        mymu = 0;
        for( i=0; i<1216; i++ ){
            dltemp = l[i] + mya*dl[i];
            dstemp = s[i] + mya*ds[i];
            if( dltemp < 0 || dstemp < 0 ){
                lsIt++;
                break;
            } else {                
                mymu += dstemp*dltemp;
            }
        }
        
        /* 
         * If no early termination of the for-loop above occurred, we
         * found the required value of a and we can quit the while loop.
         */
        if( i == 1216 ){
            break;
        } else {
            mya *= beliefPenaltyMPC_SET_LS_SCALE_AFF;
            if( mya < beliefPenaltyMPC_SET_LS_MINSTEP ){
                return beliefPenaltyMPC_NOPROGRESS;
            }
        }
    }
    
    /* return new values and iteration counter */
    *a = mya;
    *mu_aff = mymu / (beliefPenaltyMPC_FLOAT)1216;
    return lsIt;
}


/*
 * Vector subtraction x = u.*v - a where a is a scalar
*  and x,u,v are vectors of length 1216.
 */
void beliefPenaltyMPC_LA_VSUB5_1216(beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *v, beliefPenaltyMPC_FLOAT a, beliefPenaltyMPC_FLOAT *x)
{
	int i;
	for( i=0; i<1216; i++){
		x[i] = u[i]*v[i] - a;
	}
}


/*
 * Computes x=0; x(uidx) += u/su; x(vidx) -= v/sv where x is of length 62,
 * u, su, uidx are of length 22 and v, sv, vidx are of length 62.
 */
void beliefPenaltyMPC_LA_VSUB6_INDEXED_62_22_62(beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *su, int* uidx, beliefPenaltyMPC_FLOAT *v, beliefPenaltyMPC_FLOAT *sv, int* vidx, beliefPenaltyMPC_FLOAT *x)
{
	int i;
	for( i=0; i<62; i++ ){
		x[i] = 0;
	}
	for( i=0; i<22; i++){
		x[uidx[i]] += u[i]/su[i];
	}
	for( i=0; i<62; i++){
		x[vidx[i]] -= v[i]/sv[i];
	}
}


/* 
 * Computes r =  B*u
 * where B is stored in diagzero format
 */
void beliefPenaltyMPC_LA_DIAGZERO_MVM_20(beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *r)
{
	int i;

	for( i=0; i<20; i++ ){
		r[i] = B[i]*u[i];
	}	
	
}


/* 
 * Computes r = A*x + B*u
 * where A is stored in column major format
 * and B is stored in diagzero format
 */
void beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_20_62_62(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;

	for( i=0; i<20; i++ ){
		r[i] = A[k++]*x[0] + B[i]*u[i];
	}	

	for( j=1; j<62; j++ ){		
		for( i=0; i<20; i++ ){
			r[i] += A[k++]*x[j];
		}
	}
	
}


/*
 * Computes x=0; x(uidx) += u/su; x(vidx) -= v/sv where x is of length 20,
 * u, su, uidx are of length 20 and v, sv, vidx are of length 20.
 */
void beliefPenaltyMPC_LA_VSUB6_INDEXED_20_20_20(beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *su, int* uidx, beliefPenaltyMPC_FLOAT *v, beliefPenaltyMPC_FLOAT *sv, int* vidx, beliefPenaltyMPC_FLOAT *x)
{
	int i;
	for( i=0; i<20; i++ ){
		x[i] = 0;
	}
	for( i=0; i<20; i++){
		x[uidx[i]] += u[i]/su[i];
	}
	for( i=0; i<20; i++){
		x[vidx[i]] -= v[i]/sv[i];
	}
}


/* 
 * Computes r = A*x + B*u
 * where A is stored in column major format
 * and B is stored in diagzero format
 */
void beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_20_62_20(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;

	for( i=0; i<20; i++ ){
		r[i] = A[k++]*x[0] + B[i]*u[i];
	}	

	for( j=1; j<62; j++ ){		
		for( i=0; i<20; i++ ){
			r[i] += A[k++]*x[j];
		}
	}
	
}


/*
 * Vector subtraction z = x - y for vectors of length 888.
 */
void beliefPenaltyMPC_LA_VSUB_888(beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *y, beliefPenaltyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<888; i++){
		z[i] = x[i] - y[i];
	}
}


/** 
 * Computes z = -r./s - u.*y(y)
 * where all vectors except of y are of length 62 (length of y >= 62).
 */
void beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_62(beliefPenaltyMPC_FLOAT *r, beliefPenaltyMPC_FLOAT *s, beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *y, int* yidx, beliefPenaltyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<62; i++ ){
		z[i] = -r[i]/s[i] - u[i]*y[yidx[i]];
	}
}


/** 
 * Computes z = -r./s + u.*y(y)
 * where all vectors except of y are of length 22 (length of y >= 22).
 */
void beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_22(beliefPenaltyMPC_FLOAT *r, beliefPenaltyMPC_FLOAT *s, beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *y, int* yidx, beliefPenaltyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<22; i++ ){
		z[i] = -r[i]/s[i] + u[i]*y[yidx[i]];
	}
}


/** 
 * Computes z = -r./s - u.*y(y)
 * where all vectors except of y are of length 20 (length of y >= 20).
 */
void beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_20(beliefPenaltyMPC_FLOAT *r, beliefPenaltyMPC_FLOAT *s, beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *y, int* yidx, beliefPenaltyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<20; i++ ){
		z[i] = -r[i]/s[i] - u[i]*y[yidx[i]];
	}
}


/** 
 * Computes z = -r./s + u.*y(y)
 * where all vectors except of y are of length 20 (length of y >= 20).
 */
void beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_20(beliefPenaltyMPC_FLOAT *r, beliefPenaltyMPC_FLOAT *s, beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *y, int* yidx, beliefPenaltyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<20; i++ ){
		z[i] = -r[i]/s[i] + u[i]*y[yidx[i]];
	}
}


/*
 * Computes ds = -l.\(r + s.*dl) for vectors of length 1216.
 */
void beliefPenaltyMPC_LA_VSUB7_1216(beliefPenaltyMPC_FLOAT *l, beliefPenaltyMPC_FLOAT *r, beliefPenaltyMPC_FLOAT *s, beliefPenaltyMPC_FLOAT *dl, beliefPenaltyMPC_FLOAT *ds)
{
	int i;
	for( i=0; i<1216; i++){
		ds[i] = -(r[i] + s[i]*dl[i])/l[i];
	}
}


/*
 * Vector addition x = x + y for vectors of length 888.
 */
void beliefPenaltyMPC_LA_VADD_888(beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *y)
{
	int i;
	for( i=0; i<888; i++){
		x[i] += y[i];
	}
}


/*
 * Vector addition x = x + y for vectors of length 300.
 */
void beliefPenaltyMPC_LA_VADD_300(beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *y)
{
	int i;
	for( i=0; i<300; i++){
		x[i] += y[i];
	}
}


/*
 * Vector addition x = x + y for vectors of length 1216.
 */
void beliefPenaltyMPC_LA_VADD_1216(beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *y)
{
	int i;
	for( i=0; i<1216; i++){
		x[i] += y[i];
	}
}


/**
 * Backtracking line search for combined predictor/corrector step.
 * Update on variables with safety factor gamma (to keep us away from
 * boundary).
 */
int beliefPenaltyMPC_LINESEARCH_BACKTRACKING_COMBINED(beliefPenaltyMPC_FLOAT *z, beliefPenaltyMPC_FLOAT *v, beliefPenaltyMPC_FLOAT *l, beliefPenaltyMPC_FLOAT *s, beliefPenaltyMPC_FLOAT *dz, beliefPenaltyMPC_FLOAT *dv, beliefPenaltyMPC_FLOAT *dl, beliefPenaltyMPC_FLOAT *ds, beliefPenaltyMPC_FLOAT *a, beliefPenaltyMPC_FLOAT *mu)
{
    int i, lsIt=1;       
    beliefPenaltyMPC_FLOAT dltemp;
    beliefPenaltyMPC_FLOAT dstemp;    
    beliefPenaltyMPC_FLOAT a_gamma;
            
    *a = 1.0;
    while( 1 ){                        

        /* check whether search criterion is fulfilled */
        for( i=0; i<1216; i++ ){
            dltemp = l[i] + (*a)*dl[i];
            dstemp = s[i] + (*a)*ds[i];
            if( dltemp < 0 || dstemp < 0 ){
                lsIt++;
                break;
            }
        }
        
        /* 
         * If no early termination of the for-loop above occurred, we
         * found the required value of a and we can quit the while loop.
         */
        if( i == 1216 ){
            break;
        } else {
            *a *= beliefPenaltyMPC_SET_LS_SCALE;
            if( *a < beliefPenaltyMPC_SET_LS_MINSTEP ){
                return beliefPenaltyMPC_NOPROGRESS;
            }
        }
    }
    
    /* update variables with safety margin */
    a_gamma = (*a)*beliefPenaltyMPC_SET_LS_MAXSTEP;
    
    /* primal variables */
    for( i=0; i<888; i++ ){
        z[i] += a_gamma*dz[i];
    }
    
    /* equality constraint multipliers */
    for( i=0; i<300; i++ ){
        v[i] += a_gamma*dv[i];
    }
    
    /* inequality constraint multipliers & slacks, also update mu */
    *mu = 0;
    for( i=0; i<1216; i++ ){
        dltemp = l[i] + a_gamma*dl[i]; l[i] = dltemp;
        dstemp = s[i] + a_gamma*ds[i]; s[i] = dstemp;
        *mu += dltemp*dstemp;
    }
    
    *a = a_gamma;
    *mu /= (beliefPenaltyMPC_FLOAT)1216;
    return lsIt;
}




/* VARIABLE DEFINITIONS ------------------------------------------------ */
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_z[888];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_v[300];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_dz_aff[888];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_dv_aff[300];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_grad_cost[888];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_grad_eq[888];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_rd[888];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_l[1216];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_s[1216];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_lbys[1216];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_dl_aff[1216];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ds_aff[1216];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_dz_cc[888];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_dv_cc[300];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_dl_cc[1216];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ds_cc[1216];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ccrhs[1216];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_grad_ineq[888];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_H00[62] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 2.0000000000000000E+001, 2.0000000000000000E+001, 2.0000000000000000E+001, 2.0000000000000000E+001, 2.0000000000000000E+001, 2.0000000000000000E+001, 2.0000000000000000E+001, 2.0000000000000000E+001, 2.0000000000000000E+001, 2.0000000000000000E+001, 2.0000000000000000E+001, 2.0000000000000000E+001, 2.0000000000000000E+001, 2.0000000000000000E+001, 2.0000000000000000E+001, 2.0000000000000000E-002, 2.0000000000000000E-002, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_z00 = beliefPenaltyMPC_z + 0;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzaff00 = beliefPenaltyMPC_dz_aff + 0;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzcc00 = beliefPenaltyMPC_dz_cc + 0;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_rd00 = beliefPenaltyMPC_rd + 0;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lbyrd00[62];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_cost00 = beliefPenaltyMPC_grad_cost + 0;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_eq00 = beliefPenaltyMPC_grad_eq + 0;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_ineq00 = beliefPenaltyMPC_grad_ineq + 0;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ctv00[62];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_v00 = beliefPenaltyMPC_v + 0;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_re00[20];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_beta00[20];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_betacc00[20];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvaff00 = beliefPenaltyMPC_dv_aff + 0;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvcc00 = beliefPenaltyMPC_dv_cc + 0;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_V00[1240];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Yd00[210];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ld00[210];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_yy00[20];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_bmy00[20];
int beliefPenaltyMPC_lbIdx00[62] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llb00 = beliefPenaltyMPC_l + 0;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_slb00 = beliefPenaltyMPC_s + 0;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llbbyslb00 = beliefPenaltyMPC_lbys + 0;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_rilb00[62];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbaff00 = beliefPenaltyMPC_dl_aff + 0;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbaff00 = beliefPenaltyMPC_ds_aff + 0;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbcc00 = beliefPenaltyMPC_dl_cc + 0;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbcc00 = beliefPenaltyMPC_ds_cc + 0;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsl00 = beliefPenaltyMPC_ccrhs + 0;
int beliefPenaltyMPC_ubIdx00[22] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lub00 = beliefPenaltyMPC_l + 62;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_sub00 = beliefPenaltyMPC_s + 62;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lubbysub00 = beliefPenaltyMPC_lbys + 62;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_riub00[22];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubaff00 = beliefPenaltyMPC_dl_aff + 62;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubaff00 = beliefPenaltyMPC_ds_aff + 62;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubcc00 = beliefPenaltyMPC_dl_cc + 62;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubcc00 = beliefPenaltyMPC_ds_cc + 62;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsub00 = beliefPenaltyMPC_ccrhs + 62;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Phi00[62];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_D00[62] = {1.0000000000000000E+000, 
1.0000000000000000E+000, 
1.0000000000000000E+000, 
1.0000000000000000E+000, 
1.0000000000000000E+000, 
1.0000000000000000E+000, 
1.0000000000000000E+000, 
1.0000000000000000E+000, 
1.0000000000000000E+000, 
1.0000000000000000E+000, 
1.0000000000000000E+000, 
1.0000000000000000E+000, 
1.0000000000000000E+000, 
1.0000000000000000E+000, 
1.0000000000000000E+000, 
1.0000000000000000E+000, 
1.0000000000000000E+000, 
1.0000000000000000E+000, 
1.0000000000000000E+000, 
1.0000000000000000E+000};
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_W00[62];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_z01 = beliefPenaltyMPC_z + 62;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzaff01 = beliefPenaltyMPC_dz_aff + 62;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzcc01 = beliefPenaltyMPC_dz_cc + 62;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_rd01 = beliefPenaltyMPC_rd + 62;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lbyrd01[62];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_cost01 = beliefPenaltyMPC_grad_cost + 62;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_eq01 = beliefPenaltyMPC_grad_eq + 62;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_ineq01 = beliefPenaltyMPC_grad_ineq + 62;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ctv01[62];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_v01 = beliefPenaltyMPC_v + 20;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_re01[20];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_beta01[20];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_betacc01[20];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvaff01 = beliefPenaltyMPC_dv_aff + 20;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvcc01 = beliefPenaltyMPC_dv_cc + 20;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_V01[1240];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Yd01[210];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ld01[210];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_yy01[20];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_bmy01[20];
int beliefPenaltyMPC_lbIdx01[62] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llb01 = beliefPenaltyMPC_l + 84;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_slb01 = beliefPenaltyMPC_s + 84;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llbbyslb01 = beliefPenaltyMPC_lbys + 84;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_rilb01[62];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbaff01 = beliefPenaltyMPC_dl_aff + 84;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbaff01 = beliefPenaltyMPC_ds_aff + 84;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbcc01 = beliefPenaltyMPC_dl_cc + 84;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbcc01 = beliefPenaltyMPC_ds_cc + 84;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsl01 = beliefPenaltyMPC_ccrhs + 84;
int beliefPenaltyMPC_ubIdx01[22] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lub01 = beliefPenaltyMPC_l + 146;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_sub01 = beliefPenaltyMPC_s + 146;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lubbysub01 = beliefPenaltyMPC_lbys + 146;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_riub01[22];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubaff01 = beliefPenaltyMPC_dl_aff + 146;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubaff01 = beliefPenaltyMPC_ds_aff + 146;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubcc01 = beliefPenaltyMPC_dl_cc + 146;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubcc01 = beliefPenaltyMPC_ds_cc + 146;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsub01 = beliefPenaltyMPC_ccrhs + 146;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Phi01[62];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_D01[62] = {-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000};
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_W01[62];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ysd01[400];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lsd01[400];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_z02 = beliefPenaltyMPC_z + 124;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzaff02 = beliefPenaltyMPC_dz_aff + 124;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzcc02 = beliefPenaltyMPC_dz_cc + 124;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_rd02 = beliefPenaltyMPC_rd + 124;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lbyrd02[62];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_cost02 = beliefPenaltyMPC_grad_cost + 124;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_eq02 = beliefPenaltyMPC_grad_eq + 124;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_ineq02 = beliefPenaltyMPC_grad_ineq + 124;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ctv02[62];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_v02 = beliefPenaltyMPC_v + 40;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_re02[20];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_beta02[20];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_betacc02[20];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvaff02 = beliefPenaltyMPC_dv_aff + 40;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvcc02 = beliefPenaltyMPC_dv_cc + 40;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_V02[1240];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Yd02[210];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ld02[210];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_yy02[20];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_bmy02[20];
int beliefPenaltyMPC_lbIdx02[62] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llb02 = beliefPenaltyMPC_l + 168;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_slb02 = beliefPenaltyMPC_s + 168;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llbbyslb02 = beliefPenaltyMPC_lbys + 168;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_rilb02[62];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbaff02 = beliefPenaltyMPC_dl_aff + 168;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbaff02 = beliefPenaltyMPC_ds_aff + 168;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbcc02 = beliefPenaltyMPC_dl_cc + 168;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbcc02 = beliefPenaltyMPC_ds_cc + 168;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsl02 = beliefPenaltyMPC_ccrhs + 168;
int beliefPenaltyMPC_ubIdx02[22] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lub02 = beliefPenaltyMPC_l + 230;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_sub02 = beliefPenaltyMPC_s + 230;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lubbysub02 = beliefPenaltyMPC_lbys + 230;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_riub02[22];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubaff02 = beliefPenaltyMPC_dl_aff + 230;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubaff02 = beliefPenaltyMPC_ds_aff + 230;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubcc02 = beliefPenaltyMPC_dl_cc + 230;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubcc02 = beliefPenaltyMPC_ds_cc + 230;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsub02 = beliefPenaltyMPC_ccrhs + 230;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Phi02[62];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_W02[62];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ysd02[400];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lsd02[400];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_z03 = beliefPenaltyMPC_z + 186;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzaff03 = beliefPenaltyMPC_dz_aff + 186;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzcc03 = beliefPenaltyMPC_dz_cc + 186;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_rd03 = beliefPenaltyMPC_rd + 186;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lbyrd03[62];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_cost03 = beliefPenaltyMPC_grad_cost + 186;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_eq03 = beliefPenaltyMPC_grad_eq + 186;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_ineq03 = beliefPenaltyMPC_grad_ineq + 186;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ctv03[62];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_v03 = beliefPenaltyMPC_v + 60;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_re03[20];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_beta03[20];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_betacc03[20];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvaff03 = beliefPenaltyMPC_dv_aff + 60;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvcc03 = beliefPenaltyMPC_dv_cc + 60;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_V03[1240];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Yd03[210];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ld03[210];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_yy03[20];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_bmy03[20];
int beliefPenaltyMPC_lbIdx03[62] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llb03 = beliefPenaltyMPC_l + 252;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_slb03 = beliefPenaltyMPC_s + 252;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llbbyslb03 = beliefPenaltyMPC_lbys + 252;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_rilb03[62];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbaff03 = beliefPenaltyMPC_dl_aff + 252;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbaff03 = beliefPenaltyMPC_ds_aff + 252;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbcc03 = beliefPenaltyMPC_dl_cc + 252;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbcc03 = beliefPenaltyMPC_ds_cc + 252;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsl03 = beliefPenaltyMPC_ccrhs + 252;
int beliefPenaltyMPC_ubIdx03[22] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lub03 = beliefPenaltyMPC_l + 314;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_sub03 = beliefPenaltyMPC_s + 314;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lubbysub03 = beliefPenaltyMPC_lbys + 314;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_riub03[22];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubaff03 = beliefPenaltyMPC_dl_aff + 314;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubaff03 = beliefPenaltyMPC_ds_aff + 314;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubcc03 = beliefPenaltyMPC_dl_cc + 314;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubcc03 = beliefPenaltyMPC_ds_cc + 314;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsub03 = beliefPenaltyMPC_ccrhs + 314;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Phi03[62];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_W03[62];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ysd03[400];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lsd03[400];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_z04 = beliefPenaltyMPC_z + 248;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzaff04 = beliefPenaltyMPC_dz_aff + 248;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzcc04 = beliefPenaltyMPC_dz_cc + 248;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_rd04 = beliefPenaltyMPC_rd + 248;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lbyrd04[62];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_cost04 = beliefPenaltyMPC_grad_cost + 248;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_eq04 = beliefPenaltyMPC_grad_eq + 248;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_ineq04 = beliefPenaltyMPC_grad_ineq + 248;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ctv04[62];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_v04 = beliefPenaltyMPC_v + 80;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_re04[20];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_beta04[20];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_betacc04[20];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvaff04 = beliefPenaltyMPC_dv_aff + 80;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvcc04 = beliefPenaltyMPC_dv_cc + 80;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_V04[1240];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Yd04[210];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ld04[210];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_yy04[20];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_bmy04[20];
int beliefPenaltyMPC_lbIdx04[62] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llb04 = beliefPenaltyMPC_l + 336;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_slb04 = beliefPenaltyMPC_s + 336;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llbbyslb04 = beliefPenaltyMPC_lbys + 336;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_rilb04[62];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbaff04 = beliefPenaltyMPC_dl_aff + 336;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbaff04 = beliefPenaltyMPC_ds_aff + 336;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbcc04 = beliefPenaltyMPC_dl_cc + 336;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbcc04 = beliefPenaltyMPC_ds_cc + 336;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsl04 = beliefPenaltyMPC_ccrhs + 336;
int beliefPenaltyMPC_ubIdx04[22] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lub04 = beliefPenaltyMPC_l + 398;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_sub04 = beliefPenaltyMPC_s + 398;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lubbysub04 = beliefPenaltyMPC_lbys + 398;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_riub04[22];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubaff04 = beliefPenaltyMPC_dl_aff + 398;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubaff04 = beliefPenaltyMPC_ds_aff + 398;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubcc04 = beliefPenaltyMPC_dl_cc + 398;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubcc04 = beliefPenaltyMPC_ds_cc + 398;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsub04 = beliefPenaltyMPC_ccrhs + 398;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Phi04[62];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_W04[62];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ysd04[400];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lsd04[400];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_z05 = beliefPenaltyMPC_z + 310;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzaff05 = beliefPenaltyMPC_dz_aff + 310;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzcc05 = beliefPenaltyMPC_dz_cc + 310;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_rd05 = beliefPenaltyMPC_rd + 310;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lbyrd05[62];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_cost05 = beliefPenaltyMPC_grad_cost + 310;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_eq05 = beliefPenaltyMPC_grad_eq + 310;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_ineq05 = beliefPenaltyMPC_grad_ineq + 310;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ctv05[62];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_v05 = beliefPenaltyMPC_v + 100;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_re05[20];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_beta05[20];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_betacc05[20];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvaff05 = beliefPenaltyMPC_dv_aff + 100;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvcc05 = beliefPenaltyMPC_dv_cc + 100;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_V05[1240];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Yd05[210];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ld05[210];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_yy05[20];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_bmy05[20];
int beliefPenaltyMPC_lbIdx05[62] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llb05 = beliefPenaltyMPC_l + 420;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_slb05 = beliefPenaltyMPC_s + 420;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llbbyslb05 = beliefPenaltyMPC_lbys + 420;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_rilb05[62];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbaff05 = beliefPenaltyMPC_dl_aff + 420;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbaff05 = beliefPenaltyMPC_ds_aff + 420;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbcc05 = beliefPenaltyMPC_dl_cc + 420;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbcc05 = beliefPenaltyMPC_ds_cc + 420;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsl05 = beliefPenaltyMPC_ccrhs + 420;
int beliefPenaltyMPC_ubIdx05[22] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lub05 = beliefPenaltyMPC_l + 482;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_sub05 = beliefPenaltyMPC_s + 482;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lubbysub05 = beliefPenaltyMPC_lbys + 482;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_riub05[22];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubaff05 = beliefPenaltyMPC_dl_aff + 482;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubaff05 = beliefPenaltyMPC_ds_aff + 482;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubcc05 = beliefPenaltyMPC_dl_cc + 482;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubcc05 = beliefPenaltyMPC_ds_cc + 482;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsub05 = beliefPenaltyMPC_ccrhs + 482;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Phi05[62];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_W05[62];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ysd05[400];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lsd05[400];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_z06 = beliefPenaltyMPC_z + 372;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzaff06 = beliefPenaltyMPC_dz_aff + 372;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzcc06 = beliefPenaltyMPC_dz_cc + 372;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_rd06 = beliefPenaltyMPC_rd + 372;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lbyrd06[62];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_cost06 = beliefPenaltyMPC_grad_cost + 372;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_eq06 = beliefPenaltyMPC_grad_eq + 372;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_ineq06 = beliefPenaltyMPC_grad_ineq + 372;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ctv06[62];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_v06 = beliefPenaltyMPC_v + 120;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_re06[20];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_beta06[20];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_betacc06[20];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvaff06 = beliefPenaltyMPC_dv_aff + 120;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvcc06 = beliefPenaltyMPC_dv_cc + 120;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_V06[1240];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Yd06[210];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ld06[210];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_yy06[20];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_bmy06[20];
int beliefPenaltyMPC_lbIdx06[62] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llb06 = beliefPenaltyMPC_l + 504;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_slb06 = beliefPenaltyMPC_s + 504;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llbbyslb06 = beliefPenaltyMPC_lbys + 504;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_rilb06[62];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbaff06 = beliefPenaltyMPC_dl_aff + 504;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbaff06 = beliefPenaltyMPC_ds_aff + 504;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbcc06 = beliefPenaltyMPC_dl_cc + 504;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbcc06 = beliefPenaltyMPC_ds_cc + 504;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsl06 = beliefPenaltyMPC_ccrhs + 504;
int beliefPenaltyMPC_ubIdx06[22] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lub06 = beliefPenaltyMPC_l + 566;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_sub06 = beliefPenaltyMPC_s + 566;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lubbysub06 = beliefPenaltyMPC_lbys + 566;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_riub06[22];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubaff06 = beliefPenaltyMPC_dl_aff + 566;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubaff06 = beliefPenaltyMPC_ds_aff + 566;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubcc06 = beliefPenaltyMPC_dl_cc + 566;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubcc06 = beliefPenaltyMPC_ds_cc + 566;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsub06 = beliefPenaltyMPC_ccrhs + 566;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Phi06[62];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_W06[62];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ysd06[400];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lsd06[400];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_z07 = beliefPenaltyMPC_z + 434;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzaff07 = beliefPenaltyMPC_dz_aff + 434;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzcc07 = beliefPenaltyMPC_dz_cc + 434;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_rd07 = beliefPenaltyMPC_rd + 434;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lbyrd07[62];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_cost07 = beliefPenaltyMPC_grad_cost + 434;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_eq07 = beliefPenaltyMPC_grad_eq + 434;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_ineq07 = beliefPenaltyMPC_grad_ineq + 434;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ctv07[62];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_v07 = beliefPenaltyMPC_v + 140;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_re07[20];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_beta07[20];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_betacc07[20];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvaff07 = beliefPenaltyMPC_dv_aff + 140;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvcc07 = beliefPenaltyMPC_dv_cc + 140;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_V07[1240];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Yd07[210];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ld07[210];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_yy07[20];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_bmy07[20];
int beliefPenaltyMPC_lbIdx07[62] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llb07 = beliefPenaltyMPC_l + 588;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_slb07 = beliefPenaltyMPC_s + 588;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llbbyslb07 = beliefPenaltyMPC_lbys + 588;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_rilb07[62];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbaff07 = beliefPenaltyMPC_dl_aff + 588;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbaff07 = beliefPenaltyMPC_ds_aff + 588;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbcc07 = beliefPenaltyMPC_dl_cc + 588;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbcc07 = beliefPenaltyMPC_ds_cc + 588;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsl07 = beliefPenaltyMPC_ccrhs + 588;
int beliefPenaltyMPC_ubIdx07[22] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lub07 = beliefPenaltyMPC_l + 650;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_sub07 = beliefPenaltyMPC_s + 650;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lubbysub07 = beliefPenaltyMPC_lbys + 650;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_riub07[22];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubaff07 = beliefPenaltyMPC_dl_aff + 650;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubaff07 = beliefPenaltyMPC_ds_aff + 650;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubcc07 = beliefPenaltyMPC_dl_cc + 650;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubcc07 = beliefPenaltyMPC_ds_cc + 650;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsub07 = beliefPenaltyMPC_ccrhs + 650;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Phi07[62];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_W07[62];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ysd07[400];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lsd07[400];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_z08 = beliefPenaltyMPC_z + 496;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzaff08 = beliefPenaltyMPC_dz_aff + 496;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzcc08 = beliefPenaltyMPC_dz_cc + 496;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_rd08 = beliefPenaltyMPC_rd + 496;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lbyrd08[62];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_cost08 = beliefPenaltyMPC_grad_cost + 496;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_eq08 = beliefPenaltyMPC_grad_eq + 496;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_ineq08 = beliefPenaltyMPC_grad_ineq + 496;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ctv08[62];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_v08 = beliefPenaltyMPC_v + 160;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_re08[20];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_beta08[20];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_betacc08[20];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvaff08 = beliefPenaltyMPC_dv_aff + 160;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvcc08 = beliefPenaltyMPC_dv_cc + 160;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_V08[1240];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Yd08[210];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ld08[210];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_yy08[20];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_bmy08[20];
int beliefPenaltyMPC_lbIdx08[62] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llb08 = beliefPenaltyMPC_l + 672;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_slb08 = beliefPenaltyMPC_s + 672;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llbbyslb08 = beliefPenaltyMPC_lbys + 672;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_rilb08[62];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbaff08 = beliefPenaltyMPC_dl_aff + 672;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbaff08 = beliefPenaltyMPC_ds_aff + 672;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbcc08 = beliefPenaltyMPC_dl_cc + 672;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbcc08 = beliefPenaltyMPC_ds_cc + 672;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsl08 = beliefPenaltyMPC_ccrhs + 672;
int beliefPenaltyMPC_ubIdx08[22] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lub08 = beliefPenaltyMPC_l + 734;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_sub08 = beliefPenaltyMPC_s + 734;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lubbysub08 = beliefPenaltyMPC_lbys + 734;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_riub08[22];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubaff08 = beliefPenaltyMPC_dl_aff + 734;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubaff08 = beliefPenaltyMPC_ds_aff + 734;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubcc08 = beliefPenaltyMPC_dl_cc + 734;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubcc08 = beliefPenaltyMPC_ds_cc + 734;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsub08 = beliefPenaltyMPC_ccrhs + 734;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Phi08[62];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_W08[62];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ysd08[400];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lsd08[400];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_z09 = beliefPenaltyMPC_z + 558;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzaff09 = beliefPenaltyMPC_dz_aff + 558;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzcc09 = beliefPenaltyMPC_dz_cc + 558;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_rd09 = beliefPenaltyMPC_rd + 558;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lbyrd09[62];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_cost09 = beliefPenaltyMPC_grad_cost + 558;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_eq09 = beliefPenaltyMPC_grad_eq + 558;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_ineq09 = beliefPenaltyMPC_grad_ineq + 558;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ctv09[62];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_v09 = beliefPenaltyMPC_v + 180;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_re09[20];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_beta09[20];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_betacc09[20];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvaff09 = beliefPenaltyMPC_dv_aff + 180;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvcc09 = beliefPenaltyMPC_dv_cc + 180;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_V09[1240];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Yd09[210];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ld09[210];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_yy09[20];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_bmy09[20];
int beliefPenaltyMPC_lbIdx09[62] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llb09 = beliefPenaltyMPC_l + 756;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_slb09 = beliefPenaltyMPC_s + 756;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llbbyslb09 = beliefPenaltyMPC_lbys + 756;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_rilb09[62];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbaff09 = beliefPenaltyMPC_dl_aff + 756;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbaff09 = beliefPenaltyMPC_ds_aff + 756;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbcc09 = beliefPenaltyMPC_dl_cc + 756;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbcc09 = beliefPenaltyMPC_ds_cc + 756;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsl09 = beliefPenaltyMPC_ccrhs + 756;
int beliefPenaltyMPC_ubIdx09[22] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lub09 = beliefPenaltyMPC_l + 818;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_sub09 = beliefPenaltyMPC_s + 818;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lubbysub09 = beliefPenaltyMPC_lbys + 818;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_riub09[22];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubaff09 = beliefPenaltyMPC_dl_aff + 818;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubaff09 = beliefPenaltyMPC_ds_aff + 818;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubcc09 = beliefPenaltyMPC_dl_cc + 818;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubcc09 = beliefPenaltyMPC_ds_cc + 818;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsub09 = beliefPenaltyMPC_ccrhs + 818;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Phi09[62];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_W09[62];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ysd09[400];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lsd09[400];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_z10 = beliefPenaltyMPC_z + 620;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzaff10 = beliefPenaltyMPC_dz_aff + 620;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzcc10 = beliefPenaltyMPC_dz_cc + 620;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_rd10 = beliefPenaltyMPC_rd + 620;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lbyrd10[62];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_cost10 = beliefPenaltyMPC_grad_cost + 620;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_eq10 = beliefPenaltyMPC_grad_eq + 620;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_ineq10 = beliefPenaltyMPC_grad_ineq + 620;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ctv10[62];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_v10 = beliefPenaltyMPC_v + 200;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_re10[20];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_beta10[20];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_betacc10[20];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvaff10 = beliefPenaltyMPC_dv_aff + 200;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvcc10 = beliefPenaltyMPC_dv_cc + 200;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_V10[1240];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Yd10[210];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ld10[210];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_yy10[20];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_bmy10[20];
int beliefPenaltyMPC_lbIdx10[62] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llb10 = beliefPenaltyMPC_l + 840;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_slb10 = beliefPenaltyMPC_s + 840;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llbbyslb10 = beliefPenaltyMPC_lbys + 840;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_rilb10[62];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbaff10 = beliefPenaltyMPC_dl_aff + 840;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbaff10 = beliefPenaltyMPC_ds_aff + 840;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbcc10 = beliefPenaltyMPC_dl_cc + 840;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbcc10 = beliefPenaltyMPC_ds_cc + 840;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsl10 = beliefPenaltyMPC_ccrhs + 840;
int beliefPenaltyMPC_ubIdx10[22] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lub10 = beliefPenaltyMPC_l + 902;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_sub10 = beliefPenaltyMPC_s + 902;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lubbysub10 = beliefPenaltyMPC_lbys + 902;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_riub10[22];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubaff10 = beliefPenaltyMPC_dl_aff + 902;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubaff10 = beliefPenaltyMPC_ds_aff + 902;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubcc10 = beliefPenaltyMPC_dl_cc + 902;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubcc10 = beliefPenaltyMPC_ds_cc + 902;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsub10 = beliefPenaltyMPC_ccrhs + 902;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Phi10[62];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_W10[62];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ysd10[400];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lsd10[400];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_z11 = beliefPenaltyMPC_z + 682;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzaff11 = beliefPenaltyMPC_dz_aff + 682;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzcc11 = beliefPenaltyMPC_dz_cc + 682;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_rd11 = beliefPenaltyMPC_rd + 682;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lbyrd11[62];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_cost11 = beliefPenaltyMPC_grad_cost + 682;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_eq11 = beliefPenaltyMPC_grad_eq + 682;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_ineq11 = beliefPenaltyMPC_grad_ineq + 682;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ctv11[62];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_v11 = beliefPenaltyMPC_v + 220;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_re11[20];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_beta11[20];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_betacc11[20];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvaff11 = beliefPenaltyMPC_dv_aff + 220;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvcc11 = beliefPenaltyMPC_dv_cc + 220;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_V11[1240];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Yd11[210];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ld11[210];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_yy11[20];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_bmy11[20];
int beliefPenaltyMPC_lbIdx11[62] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llb11 = beliefPenaltyMPC_l + 924;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_slb11 = beliefPenaltyMPC_s + 924;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llbbyslb11 = beliefPenaltyMPC_lbys + 924;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_rilb11[62];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbaff11 = beliefPenaltyMPC_dl_aff + 924;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbaff11 = beliefPenaltyMPC_ds_aff + 924;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbcc11 = beliefPenaltyMPC_dl_cc + 924;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbcc11 = beliefPenaltyMPC_ds_cc + 924;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsl11 = beliefPenaltyMPC_ccrhs + 924;
int beliefPenaltyMPC_ubIdx11[22] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lub11 = beliefPenaltyMPC_l + 986;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_sub11 = beliefPenaltyMPC_s + 986;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lubbysub11 = beliefPenaltyMPC_lbys + 986;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_riub11[22];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubaff11 = beliefPenaltyMPC_dl_aff + 986;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubaff11 = beliefPenaltyMPC_ds_aff + 986;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubcc11 = beliefPenaltyMPC_dl_cc + 986;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubcc11 = beliefPenaltyMPC_ds_cc + 986;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsub11 = beliefPenaltyMPC_ccrhs + 986;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Phi11[62];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_W11[62];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ysd11[400];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lsd11[400];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_z12 = beliefPenaltyMPC_z + 744;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzaff12 = beliefPenaltyMPC_dz_aff + 744;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzcc12 = beliefPenaltyMPC_dz_cc + 744;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_rd12 = beliefPenaltyMPC_rd + 744;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lbyrd12[62];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_cost12 = beliefPenaltyMPC_grad_cost + 744;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_eq12 = beliefPenaltyMPC_grad_eq + 744;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_ineq12 = beliefPenaltyMPC_grad_ineq + 744;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ctv12[62];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_v12 = beliefPenaltyMPC_v + 240;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_re12[20];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_beta12[20];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_betacc12[20];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvaff12 = beliefPenaltyMPC_dv_aff + 240;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvcc12 = beliefPenaltyMPC_dv_cc + 240;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_V12[1240];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Yd12[210];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ld12[210];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_yy12[20];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_bmy12[20];
int beliefPenaltyMPC_lbIdx12[62] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llb12 = beliefPenaltyMPC_l + 1008;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_slb12 = beliefPenaltyMPC_s + 1008;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llbbyslb12 = beliefPenaltyMPC_lbys + 1008;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_rilb12[62];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbaff12 = beliefPenaltyMPC_dl_aff + 1008;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbaff12 = beliefPenaltyMPC_ds_aff + 1008;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbcc12 = beliefPenaltyMPC_dl_cc + 1008;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbcc12 = beliefPenaltyMPC_ds_cc + 1008;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsl12 = beliefPenaltyMPC_ccrhs + 1008;
int beliefPenaltyMPC_ubIdx12[22] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lub12 = beliefPenaltyMPC_l + 1070;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_sub12 = beliefPenaltyMPC_s + 1070;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lubbysub12 = beliefPenaltyMPC_lbys + 1070;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_riub12[22];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubaff12 = beliefPenaltyMPC_dl_aff + 1070;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubaff12 = beliefPenaltyMPC_ds_aff + 1070;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubcc12 = beliefPenaltyMPC_dl_cc + 1070;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubcc12 = beliefPenaltyMPC_ds_cc + 1070;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsub12 = beliefPenaltyMPC_ccrhs + 1070;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Phi12[62];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_W12[62];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ysd12[400];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lsd12[400];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_z13 = beliefPenaltyMPC_z + 806;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzaff13 = beliefPenaltyMPC_dz_aff + 806;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzcc13 = beliefPenaltyMPC_dz_cc + 806;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_rd13 = beliefPenaltyMPC_rd + 806;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lbyrd13[62];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_cost13 = beliefPenaltyMPC_grad_cost + 806;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_eq13 = beliefPenaltyMPC_grad_eq + 806;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_ineq13 = beliefPenaltyMPC_grad_ineq + 806;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ctv13[62];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_v13 = beliefPenaltyMPC_v + 260;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_re13[20];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_beta13[20];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_betacc13[20];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvaff13 = beliefPenaltyMPC_dv_aff + 260;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvcc13 = beliefPenaltyMPC_dv_cc + 260;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_V13[1240];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Yd13[210];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ld13[210];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_yy13[20];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_bmy13[20];
int beliefPenaltyMPC_lbIdx13[62] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llb13 = beliefPenaltyMPC_l + 1092;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_slb13 = beliefPenaltyMPC_s + 1092;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llbbyslb13 = beliefPenaltyMPC_lbys + 1092;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_rilb13[62];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbaff13 = beliefPenaltyMPC_dl_aff + 1092;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbaff13 = beliefPenaltyMPC_ds_aff + 1092;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbcc13 = beliefPenaltyMPC_dl_cc + 1092;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbcc13 = beliefPenaltyMPC_ds_cc + 1092;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsl13 = beliefPenaltyMPC_ccrhs + 1092;
int beliefPenaltyMPC_ubIdx13[22] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lub13 = beliefPenaltyMPC_l + 1154;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_sub13 = beliefPenaltyMPC_s + 1154;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lubbysub13 = beliefPenaltyMPC_lbys + 1154;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_riub13[22];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubaff13 = beliefPenaltyMPC_dl_aff + 1154;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubaff13 = beliefPenaltyMPC_ds_aff + 1154;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubcc13 = beliefPenaltyMPC_dl_cc + 1154;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubcc13 = beliefPenaltyMPC_ds_cc + 1154;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsub13 = beliefPenaltyMPC_ccrhs + 1154;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Phi13[62];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_W13[62];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ysd13[400];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lsd13[400];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_H14[20] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 2.0000000000000000E+001, 2.0000000000000000E+001, 2.0000000000000000E+001, 2.0000000000000000E+001, 2.0000000000000000E+001, 2.0000000000000000E+001, 2.0000000000000000E+001, 2.0000000000000000E+001, 2.0000000000000000E+001, 2.0000000000000000E+001, 2.0000000000000000E+001, 2.0000000000000000E+001, 2.0000000000000000E+001, 2.0000000000000000E+001, 2.0000000000000000E+001};
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_f14[20] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_z14 = beliefPenaltyMPC_z + 868;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzaff14 = beliefPenaltyMPC_dz_aff + 868;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzcc14 = beliefPenaltyMPC_dz_cc + 868;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_rd14 = beliefPenaltyMPC_rd + 868;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lbyrd14[20];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_cost14 = beliefPenaltyMPC_grad_cost + 868;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_eq14 = beliefPenaltyMPC_grad_eq + 868;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_ineq14 = beliefPenaltyMPC_grad_ineq + 868;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ctv14[20];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_v14 = beliefPenaltyMPC_v + 280;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_re14[20];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_beta14[20];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_betacc14[20];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvaff14 = beliefPenaltyMPC_dv_aff + 280;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvcc14 = beliefPenaltyMPC_dv_cc + 280;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_V14[400];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Yd14[210];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ld14[210];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_yy14[20];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_bmy14[20];
int beliefPenaltyMPC_lbIdx14[20] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llb14 = beliefPenaltyMPC_l + 1176;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_slb14 = beliefPenaltyMPC_s + 1176;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llbbyslb14 = beliefPenaltyMPC_lbys + 1176;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_rilb14[20];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbaff14 = beliefPenaltyMPC_dl_aff + 1176;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbaff14 = beliefPenaltyMPC_ds_aff + 1176;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbcc14 = beliefPenaltyMPC_dl_cc + 1176;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbcc14 = beliefPenaltyMPC_ds_cc + 1176;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsl14 = beliefPenaltyMPC_ccrhs + 1176;
int beliefPenaltyMPC_ubIdx14[20] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lub14 = beliefPenaltyMPC_l + 1196;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_sub14 = beliefPenaltyMPC_s + 1196;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lubbysub14 = beliefPenaltyMPC_lbys + 1196;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_riub14[20];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubaff14 = beliefPenaltyMPC_dl_aff + 1196;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubaff14 = beliefPenaltyMPC_ds_aff + 1196;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubcc14 = beliefPenaltyMPC_dl_cc + 1196;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubcc14 = beliefPenaltyMPC_ds_cc + 1196;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsub14 = beliefPenaltyMPC_ccrhs + 1196;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Phi14[20];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_D14[20] = {-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000};
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_W14[20];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ysd14[400];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lsd14[400];
beliefPenaltyMPC_FLOAT musigma;
beliefPenaltyMPC_FLOAT sigma_3rdroot;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Diag1_0[62];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Diag2_0[62];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_L_0[1891];




/* SOLVER CODE --------------------------------------------------------- */
int beliefPenaltyMPC_solve(beliefPenaltyMPC_params* params, beliefPenaltyMPC_output* output, beliefPenaltyMPC_info* info)
{	
int exitcode;

#if beliefPenaltyMPC_SET_TIMING == 1
	beliefPenaltyMPC_timer solvertimer;
	beliefPenaltyMPC_tic(&solvertimer);
#endif
/* FUNCTION CALLS INTO LA LIBRARY -------------------------------------- */
info->it = 0;
beliefPenaltyMPC_LA_INITIALIZEVECTOR_888(beliefPenaltyMPC_z, 0);
beliefPenaltyMPC_LA_INITIALIZEVECTOR_300(beliefPenaltyMPC_v, 1);
beliefPenaltyMPC_LA_INITIALIZEVECTOR_1216(beliefPenaltyMPC_l, 1);
beliefPenaltyMPC_LA_INITIALIZEVECTOR_1216(beliefPenaltyMPC_s, 1);
info->mu = 0;
beliefPenaltyMPC_LA_DOTACC_1216(beliefPenaltyMPC_l, beliefPenaltyMPC_s, &info->mu);
info->mu /= 1216;
PRINTTEXT("This is beliefPenaltyMPC, a solver generated by FORCES (forces.ethz.ch).\n");
PRINTTEXT("(c) Alexander Domahidi, Automatic Control Laboratory, ETH Zurich, 2011-2014.\n");
PRINTTEXT("\n  #it  res_eq   res_ineq     pobj         dobj       dgap     rdgap     mu\n");
PRINTTEXT("  ---------------------------------------------------------------------------\n");
while( 1 ){
info->pobj = 0;
beliefPenaltyMPC_LA_DIAG_QUADFCN_62(beliefPenaltyMPC_H00, params->f1, beliefPenaltyMPC_z00, beliefPenaltyMPC_grad_cost00, &info->pobj);
beliefPenaltyMPC_LA_DIAG_QUADFCN_62(beliefPenaltyMPC_H00, params->f2, beliefPenaltyMPC_z01, beliefPenaltyMPC_grad_cost01, &info->pobj);
beliefPenaltyMPC_LA_DIAG_QUADFCN_62(beliefPenaltyMPC_H00, params->f3, beliefPenaltyMPC_z02, beliefPenaltyMPC_grad_cost02, &info->pobj);
beliefPenaltyMPC_LA_DIAG_QUADFCN_62(beliefPenaltyMPC_H00, params->f4, beliefPenaltyMPC_z03, beliefPenaltyMPC_grad_cost03, &info->pobj);
beliefPenaltyMPC_LA_DIAG_QUADFCN_62(beliefPenaltyMPC_H00, params->f5, beliefPenaltyMPC_z04, beliefPenaltyMPC_grad_cost04, &info->pobj);
beliefPenaltyMPC_LA_DIAG_QUADFCN_62(beliefPenaltyMPC_H00, params->f6, beliefPenaltyMPC_z05, beliefPenaltyMPC_grad_cost05, &info->pobj);
beliefPenaltyMPC_LA_DIAG_QUADFCN_62(beliefPenaltyMPC_H00, params->f7, beliefPenaltyMPC_z06, beliefPenaltyMPC_grad_cost06, &info->pobj);
beliefPenaltyMPC_LA_DIAG_QUADFCN_62(beliefPenaltyMPC_H00, params->f8, beliefPenaltyMPC_z07, beliefPenaltyMPC_grad_cost07, &info->pobj);
beliefPenaltyMPC_LA_DIAG_QUADFCN_62(beliefPenaltyMPC_H00, params->f9, beliefPenaltyMPC_z08, beliefPenaltyMPC_grad_cost08, &info->pobj);
beliefPenaltyMPC_LA_DIAG_QUADFCN_62(beliefPenaltyMPC_H00, params->f10, beliefPenaltyMPC_z09, beliefPenaltyMPC_grad_cost09, &info->pobj);
beliefPenaltyMPC_LA_DIAG_QUADFCN_62(beliefPenaltyMPC_H00, params->f11, beliefPenaltyMPC_z10, beliefPenaltyMPC_grad_cost10, &info->pobj);
beliefPenaltyMPC_LA_DIAG_QUADFCN_62(beliefPenaltyMPC_H00, params->f12, beliefPenaltyMPC_z11, beliefPenaltyMPC_grad_cost11, &info->pobj);
beliefPenaltyMPC_LA_DIAG_QUADFCN_62(beliefPenaltyMPC_H00, params->f13, beliefPenaltyMPC_z12, beliefPenaltyMPC_grad_cost12, &info->pobj);
beliefPenaltyMPC_LA_DIAG_QUADFCN_62(beliefPenaltyMPC_H00, params->f14, beliefPenaltyMPC_z13, beliefPenaltyMPC_grad_cost13, &info->pobj);
beliefPenaltyMPC_LA_DIAG_QUADFCN_20(beliefPenaltyMPC_H14, beliefPenaltyMPC_f14, beliefPenaltyMPC_z14, beliefPenaltyMPC_grad_cost14, &info->pobj);
info->res_eq = 0;
info->dgap = 0;
beliefPenaltyMPC_LA_DIAGZERO_MVMSUB6_20(beliefPenaltyMPC_D00, beliefPenaltyMPC_z00, params->e1, beliefPenaltyMPC_v00, beliefPenaltyMPC_re00, &info->dgap, &info->res_eq);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_20_62_62(params->C1, beliefPenaltyMPC_z00, beliefPenaltyMPC_D01, beliefPenaltyMPC_z01, params->e2, beliefPenaltyMPC_v01, beliefPenaltyMPC_re01, &info->dgap, &info->res_eq);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_20_62_62(params->C2, beliefPenaltyMPC_z01, beliefPenaltyMPC_D01, beliefPenaltyMPC_z02, params->e3, beliefPenaltyMPC_v02, beliefPenaltyMPC_re02, &info->dgap, &info->res_eq);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_20_62_62(params->C3, beliefPenaltyMPC_z02, beliefPenaltyMPC_D01, beliefPenaltyMPC_z03, params->e4, beliefPenaltyMPC_v03, beliefPenaltyMPC_re03, &info->dgap, &info->res_eq);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_20_62_62(params->C4, beliefPenaltyMPC_z03, beliefPenaltyMPC_D01, beliefPenaltyMPC_z04, params->e5, beliefPenaltyMPC_v04, beliefPenaltyMPC_re04, &info->dgap, &info->res_eq);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_20_62_62(params->C5, beliefPenaltyMPC_z04, beliefPenaltyMPC_D01, beliefPenaltyMPC_z05, params->e6, beliefPenaltyMPC_v05, beliefPenaltyMPC_re05, &info->dgap, &info->res_eq);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_20_62_62(params->C6, beliefPenaltyMPC_z05, beliefPenaltyMPC_D01, beliefPenaltyMPC_z06, params->e7, beliefPenaltyMPC_v06, beliefPenaltyMPC_re06, &info->dgap, &info->res_eq);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_20_62_62(params->C7, beliefPenaltyMPC_z06, beliefPenaltyMPC_D01, beliefPenaltyMPC_z07, params->e8, beliefPenaltyMPC_v07, beliefPenaltyMPC_re07, &info->dgap, &info->res_eq);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_20_62_62(params->C8, beliefPenaltyMPC_z07, beliefPenaltyMPC_D01, beliefPenaltyMPC_z08, params->e9, beliefPenaltyMPC_v08, beliefPenaltyMPC_re08, &info->dgap, &info->res_eq);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_20_62_62(params->C9, beliefPenaltyMPC_z08, beliefPenaltyMPC_D01, beliefPenaltyMPC_z09, params->e10, beliefPenaltyMPC_v09, beliefPenaltyMPC_re09, &info->dgap, &info->res_eq);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_20_62_62(params->C10, beliefPenaltyMPC_z09, beliefPenaltyMPC_D01, beliefPenaltyMPC_z10, params->e11, beliefPenaltyMPC_v10, beliefPenaltyMPC_re10, &info->dgap, &info->res_eq);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_20_62_62(params->C11, beliefPenaltyMPC_z10, beliefPenaltyMPC_D01, beliefPenaltyMPC_z11, params->e12, beliefPenaltyMPC_v11, beliefPenaltyMPC_re11, &info->dgap, &info->res_eq);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_20_62_62(params->C12, beliefPenaltyMPC_z11, beliefPenaltyMPC_D01, beliefPenaltyMPC_z12, params->e13, beliefPenaltyMPC_v12, beliefPenaltyMPC_re12, &info->dgap, &info->res_eq);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_20_62_62(params->C13, beliefPenaltyMPC_z12, beliefPenaltyMPC_D01, beliefPenaltyMPC_z13, params->e14, beliefPenaltyMPC_v13, beliefPenaltyMPC_re13, &info->dgap, &info->res_eq);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_20_62_20(params->C14, beliefPenaltyMPC_z13, beliefPenaltyMPC_D14, beliefPenaltyMPC_z14, params->e15, beliefPenaltyMPC_v14, beliefPenaltyMPC_re14, &info->dgap, &info->res_eq);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_20_62_20(params->C1, beliefPenaltyMPC_v01, beliefPenaltyMPC_D00, beliefPenaltyMPC_v00, beliefPenaltyMPC_grad_eq00);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_20_62_20(params->C2, beliefPenaltyMPC_v02, beliefPenaltyMPC_D01, beliefPenaltyMPC_v01, beliefPenaltyMPC_grad_eq01);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_20_62_20(params->C3, beliefPenaltyMPC_v03, beliefPenaltyMPC_D01, beliefPenaltyMPC_v02, beliefPenaltyMPC_grad_eq02);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_20_62_20(params->C4, beliefPenaltyMPC_v04, beliefPenaltyMPC_D01, beliefPenaltyMPC_v03, beliefPenaltyMPC_grad_eq03);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_20_62_20(params->C5, beliefPenaltyMPC_v05, beliefPenaltyMPC_D01, beliefPenaltyMPC_v04, beliefPenaltyMPC_grad_eq04);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_20_62_20(params->C6, beliefPenaltyMPC_v06, beliefPenaltyMPC_D01, beliefPenaltyMPC_v05, beliefPenaltyMPC_grad_eq05);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_20_62_20(params->C7, beliefPenaltyMPC_v07, beliefPenaltyMPC_D01, beliefPenaltyMPC_v06, beliefPenaltyMPC_grad_eq06);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_20_62_20(params->C8, beliefPenaltyMPC_v08, beliefPenaltyMPC_D01, beliefPenaltyMPC_v07, beliefPenaltyMPC_grad_eq07);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_20_62_20(params->C9, beliefPenaltyMPC_v09, beliefPenaltyMPC_D01, beliefPenaltyMPC_v08, beliefPenaltyMPC_grad_eq08);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_20_62_20(params->C10, beliefPenaltyMPC_v10, beliefPenaltyMPC_D01, beliefPenaltyMPC_v09, beliefPenaltyMPC_grad_eq09);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_20_62_20(params->C11, beliefPenaltyMPC_v11, beliefPenaltyMPC_D01, beliefPenaltyMPC_v10, beliefPenaltyMPC_grad_eq10);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_20_62_20(params->C12, beliefPenaltyMPC_v12, beliefPenaltyMPC_D01, beliefPenaltyMPC_v11, beliefPenaltyMPC_grad_eq11);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_20_62_20(params->C13, beliefPenaltyMPC_v13, beliefPenaltyMPC_D01, beliefPenaltyMPC_v12, beliefPenaltyMPC_grad_eq12);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_20_62_20(params->C14, beliefPenaltyMPC_v14, beliefPenaltyMPC_D01, beliefPenaltyMPC_v13, beliefPenaltyMPC_grad_eq13);
beliefPenaltyMPC_LA_DIAGZERO_MTVM_20_20(beliefPenaltyMPC_D14, beliefPenaltyMPC_v14, beliefPenaltyMPC_grad_eq14);
info->res_ineq = 0;
beliefPenaltyMPC_LA_VSUBADD3_62(params->lb1, beliefPenaltyMPC_z00, beliefPenaltyMPC_lbIdx00, beliefPenaltyMPC_llb00, beliefPenaltyMPC_slb00, beliefPenaltyMPC_rilb00, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD2_22(beliefPenaltyMPC_z00, beliefPenaltyMPC_ubIdx00, params->ub1, beliefPenaltyMPC_lub00, beliefPenaltyMPC_sub00, beliefPenaltyMPC_riub00, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD3_62(params->lb2, beliefPenaltyMPC_z01, beliefPenaltyMPC_lbIdx01, beliefPenaltyMPC_llb01, beliefPenaltyMPC_slb01, beliefPenaltyMPC_rilb01, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD2_22(beliefPenaltyMPC_z01, beliefPenaltyMPC_ubIdx01, params->ub2, beliefPenaltyMPC_lub01, beliefPenaltyMPC_sub01, beliefPenaltyMPC_riub01, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD3_62(params->lb3, beliefPenaltyMPC_z02, beliefPenaltyMPC_lbIdx02, beliefPenaltyMPC_llb02, beliefPenaltyMPC_slb02, beliefPenaltyMPC_rilb02, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD2_22(beliefPenaltyMPC_z02, beliefPenaltyMPC_ubIdx02, params->ub3, beliefPenaltyMPC_lub02, beliefPenaltyMPC_sub02, beliefPenaltyMPC_riub02, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD3_62(params->lb4, beliefPenaltyMPC_z03, beliefPenaltyMPC_lbIdx03, beliefPenaltyMPC_llb03, beliefPenaltyMPC_slb03, beliefPenaltyMPC_rilb03, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD2_22(beliefPenaltyMPC_z03, beliefPenaltyMPC_ubIdx03, params->ub4, beliefPenaltyMPC_lub03, beliefPenaltyMPC_sub03, beliefPenaltyMPC_riub03, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD3_62(params->lb5, beliefPenaltyMPC_z04, beliefPenaltyMPC_lbIdx04, beliefPenaltyMPC_llb04, beliefPenaltyMPC_slb04, beliefPenaltyMPC_rilb04, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD2_22(beliefPenaltyMPC_z04, beliefPenaltyMPC_ubIdx04, params->ub5, beliefPenaltyMPC_lub04, beliefPenaltyMPC_sub04, beliefPenaltyMPC_riub04, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD3_62(params->lb6, beliefPenaltyMPC_z05, beliefPenaltyMPC_lbIdx05, beliefPenaltyMPC_llb05, beliefPenaltyMPC_slb05, beliefPenaltyMPC_rilb05, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD2_22(beliefPenaltyMPC_z05, beliefPenaltyMPC_ubIdx05, params->ub6, beliefPenaltyMPC_lub05, beliefPenaltyMPC_sub05, beliefPenaltyMPC_riub05, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD3_62(params->lb7, beliefPenaltyMPC_z06, beliefPenaltyMPC_lbIdx06, beliefPenaltyMPC_llb06, beliefPenaltyMPC_slb06, beliefPenaltyMPC_rilb06, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD2_22(beliefPenaltyMPC_z06, beliefPenaltyMPC_ubIdx06, params->ub7, beliefPenaltyMPC_lub06, beliefPenaltyMPC_sub06, beliefPenaltyMPC_riub06, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD3_62(params->lb8, beliefPenaltyMPC_z07, beliefPenaltyMPC_lbIdx07, beliefPenaltyMPC_llb07, beliefPenaltyMPC_slb07, beliefPenaltyMPC_rilb07, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD2_22(beliefPenaltyMPC_z07, beliefPenaltyMPC_ubIdx07, params->ub8, beliefPenaltyMPC_lub07, beliefPenaltyMPC_sub07, beliefPenaltyMPC_riub07, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD3_62(params->lb9, beliefPenaltyMPC_z08, beliefPenaltyMPC_lbIdx08, beliefPenaltyMPC_llb08, beliefPenaltyMPC_slb08, beliefPenaltyMPC_rilb08, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD2_22(beliefPenaltyMPC_z08, beliefPenaltyMPC_ubIdx08, params->ub9, beliefPenaltyMPC_lub08, beliefPenaltyMPC_sub08, beliefPenaltyMPC_riub08, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD3_62(params->lb10, beliefPenaltyMPC_z09, beliefPenaltyMPC_lbIdx09, beliefPenaltyMPC_llb09, beliefPenaltyMPC_slb09, beliefPenaltyMPC_rilb09, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD2_22(beliefPenaltyMPC_z09, beliefPenaltyMPC_ubIdx09, params->ub10, beliefPenaltyMPC_lub09, beliefPenaltyMPC_sub09, beliefPenaltyMPC_riub09, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD3_62(params->lb11, beliefPenaltyMPC_z10, beliefPenaltyMPC_lbIdx10, beliefPenaltyMPC_llb10, beliefPenaltyMPC_slb10, beliefPenaltyMPC_rilb10, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD2_22(beliefPenaltyMPC_z10, beliefPenaltyMPC_ubIdx10, params->ub11, beliefPenaltyMPC_lub10, beliefPenaltyMPC_sub10, beliefPenaltyMPC_riub10, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD3_62(params->lb12, beliefPenaltyMPC_z11, beliefPenaltyMPC_lbIdx11, beliefPenaltyMPC_llb11, beliefPenaltyMPC_slb11, beliefPenaltyMPC_rilb11, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD2_22(beliefPenaltyMPC_z11, beliefPenaltyMPC_ubIdx11, params->ub12, beliefPenaltyMPC_lub11, beliefPenaltyMPC_sub11, beliefPenaltyMPC_riub11, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD3_62(params->lb13, beliefPenaltyMPC_z12, beliefPenaltyMPC_lbIdx12, beliefPenaltyMPC_llb12, beliefPenaltyMPC_slb12, beliefPenaltyMPC_rilb12, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD2_22(beliefPenaltyMPC_z12, beliefPenaltyMPC_ubIdx12, params->ub13, beliefPenaltyMPC_lub12, beliefPenaltyMPC_sub12, beliefPenaltyMPC_riub12, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD3_62(params->lb14, beliefPenaltyMPC_z13, beliefPenaltyMPC_lbIdx13, beliefPenaltyMPC_llb13, beliefPenaltyMPC_slb13, beliefPenaltyMPC_rilb13, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD2_22(beliefPenaltyMPC_z13, beliefPenaltyMPC_ubIdx13, params->ub14, beliefPenaltyMPC_lub13, beliefPenaltyMPC_sub13, beliefPenaltyMPC_riub13, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD3_20(params->lb15, beliefPenaltyMPC_z14, beliefPenaltyMPC_lbIdx14, beliefPenaltyMPC_llb14, beliefPenaltyMPC_slb14, beliefPenaltyMPC_rilb14, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD2_20(beliefPenaltyMPC_z14, beliefPenaltyMPC_ubIdx14, params->ub15, beliefPenaltyMPC_lub14, beliefPenaltyMPC_sub14, beliefPenaltyMPC_riub14, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_INEQ_B_GRAD_62_62_22(beliefPenaltyMPC_lub00, beliefPenaltyMPC_sub00, beliefPenaltyMPC_riub00, beliefPenaltyMPC_llb00, beliefPenaltyMPC_slb00, beliefPenaltyMPC_rilb00, beliefPenaltyMPC_lbIdx00, beliefPenaltyMPC_ubIdx00, beliefPenaltyMPC_grad_ineq00, beliefPenaltyMPC_lubbysub00, beliefPenaltyMPC_llbbyslb00);
beliefPenaltyMPC_LA_INEQ_B_GRAD_62_62_22(beliefPenaltyMPC_lub01, beliefPenaltyMPC_sub01, beliefPenaltyMPC_riub01, beliefPenaltyMPC_llb01, beliefPenaltyMPC_slb01, beliefPenaltyMPC_rilb01, beliefPenaltyMPC_lbIdx01, beliefPenaltyMPC_ubIdx01, beliefPenaltyMPC_grad_ineq01, beliefPenaltyMPC_lubbysub01, beliefPenaltyMPC_llbbyslb01);
beliefPenaltyMPC_LA_INEQ_B_GRAD_62_62_22(beliefPenaltyMPC_lub02, beliefPenaltyMPC_sub02, beliefPenaltyMPC_riub02, beliefPenaltyMPC_llb02, beliefPenaltyMPC_slb02, beliefPenaltyMPC_rilb02, beliefPenaltyMPC_lbIdx02, beliefPenaltyMPC_ubIdx02, beliefPenaltyMPC_grad_ineq02, beliefPenaltyMPC_lubbysub02, beliefPenaltyMPC_llbbyslb02);
beliefPenaltyMPC_LA_INEQ_B_GRAD_62_62_22(beliefPenaltyMPC_lub03, beliefPenaltyMPC_sub03, beliefPenaltyMPC_riub03, beliefPenaltyMPC_llb03, beliefPenaltyMPC_slb03, beliefPenaltyMPC_rilb03, beliefPenaltyMPC_lbIdx03, beliefPenaltyMPC_ubIdx03, beliefPenaltyMPC_grad_ineq03, beliefPenaltyMPC_lubbysub03, beliefPenaltyMPC_llbbyslb03);
beliefPenaltyMPC_LA_INEQ_B_GRAD_62_62_22(beliefPenaltyMPC_lub04, beliefPenaltyMPC_sub04, beliefPenaltyMPC_riub04, beliefPenaltyMPC_llb04, beliefPenaltyMPC_slb04, beliefPenaltyMPC_rilb04, beliefPenaltyMPC_lbIdx04, beliefPenaltyMPC_ubIdx04, beliefPenaltyMPC_grad_ineq04, beliefPenaltyMPC_lubbysub04, beliefPenaltyMPC_llbbyslb04);
beliefPenaltyMPC_LA_INEQ_B_GRAD_62_62_22(beliefPenaltyMPC_lub05, beliefPenaltyMPC_sub05, beliefPenaltyMPC_riub05, beliefPenaltyMPC_llb05, beliefPenaltyMPC_slb05, beliefPenaltyMPC_rilb05, beliefPenaltyMPC_lbIdx05, beliefPenaltyMPC_ubIdx05, beliefPenaltyMPC_grad_ineq05, beliefPenaltyMPC_lubbysub05, beliefPenaltyMPC_llbbyslb05);
beliefPenaltyMPC_LA_INEQ_B_GRAD_62_62_22(beliefPenaltyMPC_lub06, beliefPenaltyMPC_sub06, beliefPenaltyMPC_riub06, beliefPenaltyMPC_llb06, beliefPenaltyMPC_slb06, beliefPenaltyMPC_rilb06, beliefPenaltyMPC_lbIdx06, beliefPenaltyMPC_ubIdx06, beliefPenaltyMPC_grad_ineq06, beliefPenaltyMPC_lubbysub06, beliefPenaltyMPC_llbbyslb06);
beliefPenaltyMPC_LA_INEQ_B_GRAD_62_62_22(beliefPenaltyMPC_lub07, beliefPenaltyMPC_sub07, beliefPenaltyMPC_riub07, beliefPenaltyMPC_llb07, beliefPenaltyMPC_slb07, beliefPenaltyMPC_rilb07, beliefPenaltyMPC_lbIdx07, beliefPenaltyMPC_ubIdx07, beliefPenaltyMPC_grad_ineq07, beliefPenaltyMPC_lubbysub07, beliefPenaltyMPC_llbbyslb07);
beliefPenaltyMPC_LA_INEQ_B_GRAD_62_62_22(beliefPenaltyMPC_lub08, beliefPenaltyMPC_sub08, beliefPenaltyMPC_riub08, beliefPenaltyMPC_llb08, beliefPenaltyMPC_slb08, beliefPenaltyMPC_rilb08, beliefPenaltyMPC_lbIdx08, beliefPenaltyMPC_ubIdx08, beliefPenaltyMPC_grad_ineq08, beliefPenaltyMPC_lubbysub08, beliefPenaltyMPC_llbbyslb08);
beliefPenaltyMPC_LA_INEQ_B_GRAD_62_62_22(beliefPenaltyMPC_lub09, beliefPenaltyMPC_sub09, beliefPenaltyMPC_riub09, beliefPenaltyMPC_llb09, beliefPenaltyMPC_slb09, beliefPenaltyMPC_rilb09, beliefPenaltyMPC_lbIdx09, beliefPenaltyMPC_ubIdx09, beliefPenaltyMPC_grad_ineq09, beliefPenaltyMPC_lubbysub09, beliefPenaltyMPC_llbbyslb09);
beliefPenaltyMPC_LA_INEQ_B_GRAD_62_62_22(beliefPenaltyMPC_lub10, beliefPenaltyMPC_sub10, beliefPenaltyMPC_riub10, beliefPenaltyMPC_llb10, beliefPenaltyMPC_slb10, beliefPenaltyMPC_rilb10, beliefPenaltyMPC_lbIdx10, beliefPenaltyMPC_ubIdx10, beliefPenaltyMPC_grad_ineq10, beliefPenaltyMPC_lubbysub10, beliefPenaltyMPC_llbbyslb10);
beliefPenaltyMPC_LA_INEQ_B_GRAD_62_62_22(beliefPenaltyMPC_lub11, beliefPenaltyMPC_sub11, beliefPenaltyMPC_riub11, beliefPenaltyMPC_llb11, beliefPenaltyMPC_slb11, beliefPenaltyMPC_rilb11, beliefPenaltyMPC_lbIdx11, beliefPenaltyMPC_ubIdx11, beliefPenaltyMPC_grad_ineq11, beliefPenaltyMPC_lubbysub11, beliefPenaltyMPC_llbbyslb11);
beliefPenaltyMPC_LA_INEQ_B_GRAD_62_62_22(beliefPenaltyMPC_lub12, beliefPenaltyMPC_sub12, beliefPenaltyMPC_riub12, beliefPenaltyMPC_llb12, beliefPenaltyMPC_slb12, beliefPenaltyMPC_rilb12, beliefPenaltyMPC_lbIdx12, beliefPenaltyMPC_ubIdx12, beliefPenaltyMPC_grad_ineq12, beliefPenaltyMPC_lubbysub12, beliefPenaltyMPC_llbbyslb12);
beliefPenaltyMPC_LA_INEQ_B_GRAD_62_62_22(beliefPenaltyMPC_lub13, beliefPenaltyMPC_sub13, beliefPenaltyMPC_riub13, beliefPenaltyMPC_llb13, beliefPenaltyMPC_slb13, beliefPenaltyMPC_rilb13, beliefPenaltyMPC_lbIdx13, beliefPenaltyMPC_ubIdx13, beliefPenaltyMPC_grad_ineq13, beliefPenaltyMPC_lubbysub13, beliefPenaltyMPC_llbbyslb13);
beliefPenaltyMPC_LA_INEQ_B_GRAD_20_20_20(beliefPenaltyMPC_lub14, beliefPenaltyMPC_sub14, beliefPenaltyMPC_riub14, beliefPenaltyMPC_llb14, beliefPenaltyMPC_slb14, beliefPenaltyMPC_rilb14, beliefPenaltyMPC_lbIdx14, beliefPenaltyMPC_ubIdx14, beliefPenaltyMPC_grad_ineq14, beliefPenaltyMPC_lubbysub14, beliefPenaltyMPC_llbbyslb14);
info->dobj = info->pobj - info->dgap;
info->rdgap = info->pobj ? info->dgap / info->pobj : 1e6;
if( info->rdgap < 0 ) info->rdgap = -info->rdgap;
PRINTTEXT("  %3d  %3.1e  %3.1e  %+6.4e  %+6.4e  %+3.1e  %3.1e  %3.1e\n",info->it, info->res_eq, info->res_ineq, info->pobj, info->dobj, info->dgap, info->rdgap, info->mu);
if( info->mu < beliefPenaltyMPC_SET_ACC_KKTCOMPL
    && (info->rdgap < beliefPenaltyMPC_SET_ACC_RDGAP || info->dgap < beliefPenaltyMPC_SET_ACC_KKTCOMPL)
    && info->res_eq < beliefPenaltyMPC_SET_ACC_RESEQ
    && info->res_ineq < beliefPenaltyMPC_SET_ACC_RESINEQ ){
PRINTTEXT("OPTIMAL (within RESEQ=%2.1e, RESINEQ=%2.1e, (R)DGAP=(%2.1e)%2.1e, MU=%2.1e).\n",beliefPenaltyMPC_SET_ACC_RESEQ, beliefPenaltyMPC_SET_ACC_RESINEQ,beliefPenaltyMPC_SET_ACC_KKTCOMPL,beliefPenaltyMPC_SET_ACC_RDGAP,beliefPenaltyMPC_SET_ACC_KKTCOMPL);
exitcode = beliefPenaltyMPC_OPTIMAL; break; }
if( info->it == beliefPenaltyMPC_SET_MAXIT ){
PRINTTEXT("Maximum number of iterations reached, exiting.\n");
exitcode = beliefPenaltyMPC_MAXITREACHED; break; }
beliefPenaltyMPC_LA_VVADD3_888(beliefPenaltyMPC_grad_cost, beliefPenaltyMPC_grad_eq, beliefPenaltyMPC_grad_ineq, beliefPenaltyMPC_rd);
beliefPenaltyMPC_LA_DIAG_CHOL_LBUB_62_62_22(beliefPenaltyMPC_H00, beliefPenaltyMPC_llbbyslb00, beliefPenaltyMPC_lbIdx00, beliefPenaltyMPC_lubbysub00, beliefPenaltyMPC_ubIdx00, beliefPenaltyMPC_Phi00);
beliefPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_20_62(beliefPenaltyMPC_Phi00, params->C1, beliefPenaltyMPC_V00);
beliefPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_20_62(beliefPenaltyMPC_Phi00, beliefPenaltyMPC_D00, beliefPenaltyMPC_W00);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMTM_20_62_20(beliefPenaltyMPC_W00, beliefPenaltyMPC_V00, beliefPenaltyMPC_Ysd01);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_62(beliefPenaltyMPC_Phi00, beliefPenaltyMPC_rd00, beliefPenaltyMPC_Lbyrd00);
beliefPenaltyMPC_LA_DIAG_CHOL_LBUB_62_62_22(beliefPenaltyMPC_H00, beliefPenaltyMPC_llbbyslb01, beliefPenaltyMPC_lbIdx01, beliefPenaltyMPC_lubbysub01, beliefPenaltyMPC_ubIdx01, beliefPenaltyMPC_Phi01);
beliefPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_20_62(beliefPenaltyMPC_Phi01, params->C2, beliefPenaltyMPC_V01);
beliefPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_20_62(beliefPenaltyMPC_Phi01, beliefPenaltyMPC_D01, beliefPenaltyMPC_W01);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMTM_20_62_20(beliefPenaltyMPC_W01, beliefPenaltyMPC_V01, beliefPenaltyMPC_Ysd02);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_62(beliefPenaltyMPC_Phi01, beliefPenaltyMPC_rd01, beliefPenaltyMPC_Lbyrd01);
beliefPenaltyMPC_LA_DIAG_CHOL_LBUB_62_62_22(beliefPenaltyMPC_H00, beliefPenaltyMPC_llbbyslb02, beliefPenaltyMPC_lbIdx02, beliefPenaltyMPC_lubbysub02, beliefPenaltyMPC_ubIdx02, beliefPenaltyMPC_Phi02);
beliefPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_20_62(beliefPenaltyMPC_Phi02, params->C3, beliefPenaltyMPC_V02);
beliefPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_20_62(beliefPenaltyMPC_Phi02, beliefPenaltyMPC_D01, beliefPenaltyMPC_W02);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMTM_20_62_20(beliefPenaltyMPC_W02, beliefPenaltyMPC_V02, beliefPenaltyMPC_Ysd03);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_62(beliefPenaltyMPC_Phi02, beliefPenaltyMPC_rd02, beliefPenaltyMPC_Lbyrd02);
beliefPenaltyMPC_LA_DIAG_CHOL_LBUB_62_62_22(beliefPenaltyMPC_H00, beliefPenaltyMPC_llbbyslb03, beliefPenaltyMPC_lbIdx03, beliefPenaltyMPC_lubbysub03, beliefPenaltyMPC_ubIdx03, beliefPenaltyMPC_Phi03);
beliefPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_20_62(beliefPenaltyMPC_Phi03, params->C4, beliefPenaltyMPC_V03);
beliefPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_20_62(beliefPenaltyMPC_Phi03, beliefPenaltyMPC_D01, beliefPenaltyMPC_W03);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMTM_20_62_20(beliefPenaltyMPC_W03, beliefPenaltyMPC_V03, beliefPenaltyMPC_Ysd04);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_62(beliefPenaltyMPC_Phi03, beliefPenaltyMPC_rd03, beliefPenaltyMPC_Lbyrd03);
beliefPenaltyMPC_LA_DIAG_CHOL_LBUB_62_62_22(beliefPenaltyMPC_H00, beliefPenaltyMPC_llbbyslb04, beliefPenaltyMPC_lbIdx04, beliefPenaltyMPC_lubbysub04, beliefPenaltyMPC_ubIdx04, beliefPenaltyMPC_Phi04);
beliefPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_20_62(beliefPenaltyMPC_Phi04, params->C5, beliefPenaltyMPC_V04);
beliefPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_20_62(beliefPenaltyMPC_Phi04, beliefPenaltyMPC_D01, beliefPenaltyMPC_W04);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMTM_20_62_20(beliefPenaltyMPC_W04, beliefPenaltyMPC_V04, beliefPenaltyMPC_Ysd05);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_62(beliefPenaltyMPC_Phi04, beliefPenaltyMPC_rd04, beliefPenaltyMPC_Lbyrd04);
beliefPenaltyMPC_LA_DIAG_CHOL_LBUB_62_62_22(beliefPenaltyMPC_H00, beliefPenaltyMPC_llbbyslb05, beliefPenaltyMPC_lbIdx05, beliefPenaltyMPC_lubbysub05, beliefPenaltyMPC_ubIdx05, beliefPenaltyMPC_Phi05);
beliefPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_20_62(beliefPenaltyMPC_Phi05, params->C6, beliefPenaltyMPC_V05);
beliefPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_20_62(beliefPenaltyMPC_Phi05, beliefPenaltyMPC_D01, beliefPenaltyMPC_W05);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMTM_20_62_20(beliefPenaltyMPC_W05, beliefPenaltyMPC_V05, beliefPenaltyMPC_Ysd06);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_62(beliefPenaltyMPC_Phi05, beliefPenaltyMPC_rd05, beliefPenaltyMPC_Lbyrd05);
beliefPenaltyMPC_LA_DIAG_CHOL_LBUB_62_62_22(beliefPenaltyMPC_H00, beliefPenaltyMPC_llbbyslb06, beliefPenaltyMPC_lbIdx06, beliefPenaltyMPC_lubbysub06, beliefPenaltyMPC_ubIdx06, beliefPenaltyMPC_Phi06);
beliefPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_20_62(beliefPenaltyMPC_Phi06, params->C7, beliefPenaltyMPC_V06);
beliefPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_20_62(beliefPenaltyMPC_Phi06, beliefPenaltyMPC_D01, beliefPenaltyMPC_W06);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMTM_20_62_20(beliefPenaltyMPC_W06, beliefPenaltyMPC_V06, beliefPenaltyMPC_Ysd07);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_62(beliefPenaltyMPC_Phi06, beliefPenaltyMPC_rd06, beliefPenaltyMPC_Lbyrd06);
beliefPenaltyMPC_LA_DIAG_CHOL_LBUB_62_62_22(beliefPenaltyMPC_H00, beliefPenaltyMPC_llbbyslb07, beliefPenaltyMPC_lbIdx07, beliefPenaltyMPC_lubbysub07, beliefPenaltyMPC_ubIdx07, beliefPenaltyMPC_Phi07);
beliefPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_20_62(beliefPenaltyMPC_Phi07, params->C8, beliefPenaltyMPC_V07);
beliefPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_20_62(beliefPenaltyMPC_Phi07, beliefPenaltyMPC_D01, beliefPenaltyMPC_W07);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMTM_20_62_20(beliefPenaltyMPC_W07, beliefPenaltyMPC_V07, beliefPenaltyMPC_Ysd08);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_62(beliefPenaltyMPC_Phi07, beliefPenaltyMPC_rd07, beliefPenaltyMPC_Lbyrd07);
beliefPenaltyMPC_LA_DIAG_CHOL_LBUB_62_62_22(beliefPenaltyMPC_H00, beliefPenaltyMPC_llbbyslb08, beliefPenaltyMPC_lbIdx08, beliefPenaltyMPC_lubbysub08, beliefPenaltyMPC_ubIdx08, beliefPenaltyMPC_Phi08);
beliefPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_20_62(beliefPenaltyMPC_Phi08, params->C9, beliefPenaltyMPC_V08);
beliefPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_20_62(beliefPenaltyMPC_Phi08, beliefPenaltyMPC_D01, beliefPenaltyMPC_W08);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMTM_20_62_20(beliefPenaltyMPC_W08, beliefPenaltyMPC_V08, beliefPenaltyMPC_Ysd09);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_62(beliefPenaltyMPC_Phi08, beliefPenaltyMPC_rd08, beliefPenaltyMPC_Lbyrd08);
beliefPenaltyMPC_LA_DIAG_CHOL_LBUB_62_62_22(beliefPenaltyMPC_H00, beliefPenaltyMPC_llbbyslb09, beliefPenaltyMPC_lbIdx09, beliefPenaltyMPC_lubbysub09, beliefPenaltyMPC_ubIdx09, beliefPenaltyMPC_Phi09);
beliefPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_20_62(beliefPenaltyMPC_Phi09, params->C10, beliefPenaltyMPC_V09);
beliefPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_20_62(beliefPenaltyMPC_Phi09, beliefPenaltyMPC_D01, beliefPenaltyMPC_W09);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMTM_20_62_20(beliefPenaltyMPC_W09, beliefPenaltyMPC_V09, beliefPenaltyMPC_Ysd10);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_62(beliefPenaltyMPC_Phi09, beliefPenaltyMPC_rd09, beliefPenaltyMPC_Lbyrd09);
beliefPenaltyMPC_LA_DIAG_CHOL_LBUB_62_62_22(beliefPenaltyMPC_H00, beliefPenaltyMPC_llbbyslb10, beliefPenaltyMPC_lbIdx10, beliefPenaltyMPC_lubbysub10, beliefPenaltyMPC_ubIdx10, beliefPenaltyMPC_Phi10);
beliefPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_20_62(beliefPenaltyMPC_Phi10, params->C11, beliefPenaltyMPC_V10);
beliefPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_20_62(beliefPenaltyMPC_Phi10, beliefPenaltyMPC_D01, beliefPenaltyMPC_W10);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMTM_20_62_20(beliefPenaltyMPC_W10, beliefPenaltyMPC_V10, beliefPenaltyMPC_Ysd11);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_62(beliefPenaltyMPC_Phi10, beliefPenaltyMPC_rd10, beliefPenaltyMPC_Lbyrd10);
beliefPenaltyMPC_LA_DIAG_CHOL_LBUB_62_62_22(beliefPenaltyMPC_H00, beliefPenaltyMPC_llbbyslb11, beliefPenaltyMPC_lbIdx11, beliefPenaltyMPC_lubbysub11, beliefPenaltyMPC_ubIdx11, beliefPenaltyMPC_Phi11);
beliefPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_20_62(beliefPenaltyMPC_Phi11, params->C12, beliefPenaltyMPC_V11);
beliefPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_20_62(beliefPenaltyMPC_Phi11, beliefPenaltyMPC_D01, beliefPenaltyMPC_W11);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMTM_20_62_20(beliefPenaltyMPC_W11, beliefPenaltyMPC_V11, beliefPenaltyMPC_Ysd12);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_62(beliefPenaltyMPC_Phi11, beliefPenaltyMPC_rd11, beliefPenaltyMPC_Lbyrd11);
beliefPenaltyMPC_LA_DIAG_CHOL_LBUB_62_62_22(beliefPenaltyMPC_H00, beliefPenaltyMPC_llbbyslb12, beliefPenaltyMPC_lbIdx12, beliefPenaltyMPC_lubbysub12, beliefPenaltyMPC_ubIdx12, beliefPenaltyMPC_Phi12);
beliefPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_20_62(beliefPenaltyMPC_Phi12, params->C13, beliefPenaltyMPC_V12);
beliefPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_20_62(beliefPenaltyMPC_Phi12, beliefPenaltyMPC_D01, beliefPenaltyMPC_W12);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMTM_20_62_20(beliefPenaltyMPC_W12, beliefPenaltyMPC_V12, beliefPenaltyMPC_Ysd13);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_62(beliefPenaltyMPC_Phi12, beliefPenaltyMPC_rd12, beliefPenaltyMPC_Lbyrd12);
beliefPenaltyMPC_LA_DIAG_CHOL_LBUB_62_62_22(beliefPenaltyMPC_H00, beliefPenaltyMPC_llbbyslb13, beliefPenaltyMPC_lbIdx13, beliefPenaltyMPC_lubbysub13, beliefPenaltyMPC_ubIdx13, beliefPenaltyMPC_Phi13);
beliefPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_20_62(beliefPenaltyMPC_Phi13, params->C14, beliefPenaltyMPC_V13);
beliefPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_20_62(beliefPenaltyMPC_Phi13, beliefPenaltyMPC_D01, beliefPenaltyMPC_W13);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMTM_20_62_20(beliefPenaltyMPC_W13, beliefPenaltyMPC_V13, beliefPenaltyMPC_Ysd14);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_62(beliefPenaltyMPC_Phi13, beliefPenaltyMPC_rd13, beliefPenaltyMPC_Lbyrd13);
beliefPenaltyMPC_LA_DIAG_CHOL_ONELOOP_LBUB_20_20_20(beliefPenaltyMPC_H14, beliefPenaltyMPC_llbbyslb14, beliefPenaltyMPC_lbIdx14, beliefPenaltyMPC_lubbysub14, beliefPenaltyMPC_ubIdx14, beliefPenaltyMPC_Phi14);
beliefPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_20_20(beliefPenaltyMPC_Phi14, beliefPenaltyMPC_D14, beliefPenaltyMPC_W14);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_20(beliefPenaltyMPC_Phi14, beliefPenaltyMPC_rd14, beliefPenaltyMPC_Lbyrd14);
beliefPenaltyMPC_LA_DIAGZERO_MMT_20(beliefPenaltyMPC_W00, beliefPenaltyMPC_Yd00);
beliefPenaltyMPC_LA_DIAGZERO_MVMSUB7_20(beliefPenaltyMPC_W00, beliefPenaltyMPC_Lbyrd00, beliefPenaltyMPC_re00, beliefPenaltyMPC_beta00);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMT2_20_62_62(beliefPenaltyMPC_V00, beliefPenaltyMPC_W01, beliefPenaltyMPC_Yd01);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_20_62_62(beliefPenaltyMPC_V00, beliefPenaltyMPC_Lbyrd00, beliefPenaltyMPC_W01, beliefPenaltyMPC_Lbyrd01, beliefPenaltyMPC_re01, beliefPenaltyMPC_beta01);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMT2_20_62_62(beliefPenaltyMPC_V01, beliefPenaltyMPC_W02, beliefPenaltyMPC_Yd02);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_20_62_62(beliefPenaltyMPC_V01, beliefPenaltyMPC_Lbyrd01, beliefPenaltyMPC_W02, beliefPenaltyMPC_Lbyrd02, beliefPenaltyMPC_re02, beliefPenaltyMPC_beta02);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMT2_20_62_62(beliefPenaltyMPC_V02, beliefPenaltyMPC_W03, beliefPenaltyMPC_Yd03);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_20_62_62(beliefPenaltyMPC_V02, beliefPenaltyMPC_Lbyrd02, beliefPenaltyMPC_W03, beliefPenaltyMPC_Lbyrd03, beliefPenaltyMPC_re03, beliefPenaltyMPC_beta03);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMT2_20_62_62(beliefPenaltyMPC_V03, beliefPenaltyMPC_W04, beliefPenaltyMPC_Yd04);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_20_62_62(beliefPenaltyMPC_V03, beliefPenaltyMPC_Lbyrd03, beliefPenaltyMPC_W04, beliefPenaltyMPC_Lbyrd04, beliefPenaltyMPC_re04, beliefPenaltyMPC_beta04);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMT2_20_62_62(beliefPenaltyMPC_V04, beliefPenaltyMPC_W05, beliefPenaltyMPC_Yd05);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_20_62_62(beliefPenaltyMPC_V04, beliefPenaltyMPC_Lbyrd04, beliefPenaltyMPC_W05, beliefPenaltyMPC_Lbyrd05, beliefPenaltyMPC_re05, beliefPenaltyMPC_beta05);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMT2_20_62_62(beliefPenaltyMPC_V05, beliefPenaltyMPC_W06, beliefPenaltyMPC_Yd06);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_20_62_62(beliefPenaltyMPC_V05, beliefPenaltyMPC_Lbyrd05, beliefPenaltyMPC_W06, beliefPenaltyMPC_Lbyrd06, beliefPenaltyMPC_re06, beliefPenaltyMPC_beta06);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMT2_20_62_62(beliefPenaltyMPC_V06, beliefPenaltyMPC_W07, beliefPenaltyMPC_Yd07);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_20_62_62(beliefPenaltyMPC_V06, beliefPenaltyMPC_Lbyrd06, beliefPenaltyMPC_W07, beliefPenaltyMPC_Lbyrd07, beliefPenaltyMPC_re07, beliefPenaltyMPC_beta07);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMT2_20_62_62(beliefPenaltyMPC_V07, beliefPenaltyMPC_W08, beliefPenaltyMPC_Yd08);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_20_62_62(beliefPenaltyMPC_V07, beliefPenaltyMPC_Lbyrd07, beliefPenaltyMPC_W08, beliefPenaltyMPC_Lbyrd08, beliefPenaltyMPC_re08, beliefPenaltyMPC_beta08);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMT2_20_62_62(beliefPenaltyMPC_V08, beliefPenaltyMPC_W09, beliefPenaltyMPC_Yd09);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_20_62_62(beliefPenaltyMPC_V08, beliefPenaltyMPC_Lbyrd08, beliefPenaltyMPC_W09, beliefPenaltyMPC_Lbyrd09, beliefPenaltyMPC_re09, beliefPenaltyMPC_beta09);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMT2_20_62_62(beliefPenaltyMPC_V09, beliefPenaltyMPC_W10, beliefPenaltyMPC_Yd10);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_20_62_62(beliefPenaltyMPC_V09, beliefPenaltyMPC_Lbyrd09, beliefPenaltyMPC_W10, beliefPenaltyMPC_Lbyrd10, beliefPenaltyMPC_re10, beliefPenaltyMPC_beta10);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMT2_20_62_62(beliefPenaltyMPC_V10, beliefPenaltyMPC_W11, beliefPenaltyMPC_Yd11);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_20_62_62(beliefPenaltyMPC_V10, beliefPenaltyMPC_Lbyrd10, beliefPenaltyMPC_W11, beliefPenaltyMPC_Lbyrd11, beliefPenaltyMPC_re11, beliefPenaltyMPC_beta11);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMT2_20_62_62(beliefPenaltyMPC_V11, beliefPenaltyMPC_W12, beliefPenaltyMPC_Yd12);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_20_62_62(beliefPenaltyMPC_V11, beliefPenaltyMPC_Lbyrd11, beliefPenaltyMPC_W12, beliefPenaltyMPC_Lbyrd12, beliefPenaltyMPC_re12, beliefPenaltyMPC_beta12);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMT2_20_62_62(beliefPenaltyMPC_V12, beliefPenaltyMPC_W13, beliefPenaltyMPC_Yd13);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_20_62_62(beliefPenaltyMPC_V12, beliefPenaltyMPC_Lbyrd12, beliefPenaltyMPC_W13, beliefPenaltyMPC_Lbyrd13, beliefPenaltyMPC_re13, beliefPenaltyMPC_beta13);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMT2_20_62_20(beliefPenaltyMPC_V13, beliefPenaltyMPC_W14, beliefPenaltyMPC_Yd14);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_20_62_20(beliefPenaltyMPC_V13, beliefPenaltyMPC_Lbyrd13, beliefPenaltyMPC_W14, beliefPenaltyMPC_Lbyrd14, beliefPenaltyMPC_re14, beliefPenaltyMPC_beta14);
beliefPenaltyMPC_LA_DENSE_CHOL_20(beliefPenaltyMPC_Yd00, beliefPenaltyMPC_Ld00);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_20(beliefPenaltyMPC_Ld00, beliefPenaltyMPC_beta00, beliefPenaltyMPC_yy00);
beliefPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_20_20(beliefPenaltyMPC_Ld00, beliefPenaltyMPC_Ysd01, beliefPenaltyMPC_Lsd01);
beliefPenaltyMPC_LA_DENSE_MMTSUB_20_20(beliefPenaltyMPC_Lsd01, beliefPenaltyMPC_Yd01);
beliefPenaltyMPC_LA_DENSE_CHOL_20(beliefPenaltyMPC_Yd01, beliefPenaltyMPC_Ld01);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_20_20(beliefPenaltyMPC_Lsd01, beliefPenaltyMPC_yy00, beliefPenaltyMPC_beta01, beliefPenaltyMPC_bmy01);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_20(beliefPenaltyMPC_Ld01, beliefPenaltyMPC_bmy01, beliefPenaltyMPC_yy01);
beliefPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_20_20(beliefPenaltyMPC_Ld01, beliefPenaltyMPC_Ysd02, beliefPenaltyMPC_Lsd02);
beliefPenaltyMPC_LA_DENSE_MMTSUB_20_20(beliefPenaltyMPC_Lsd02, beliefPenaltyMPC_Yd02);
beliefPenaltyMPC_LA_DENSE_CHOL_20(beliefPenaltyMPC_Yd02, beliefPenaltyMPC_Ld02);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_20_20(beliefPenaltyMPC_Lsd02, beliefPenaltyMPC_yy01, beliefPenaltyMPC_beta02, beliefPenaltyMPC_bmy02);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_20(beliefPenaltyMPC_Ld02, beliefPenaltyMPC_bmy02, beliefPenaltyMPC_yy02);
beliefPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_20_20(beliefPenaltyMPC_Ld02, beliefPenaltyMPC_Ysd03, beliefPenaltyMPC_Lsd03);
beliefPenaltyMPC_LA_DENSE_MMTSUB_20_20(beliefPenaltyMPC_Lsd03, beliefPenaltyMPC_Yd03);
beliefPenaltyMPC_LA_DENSE_CHOL_20(beliefPenaltyMPC_Yd03, beliefPenaltyMPC_Ld03);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_20_20(beliefPenaltyMPC_Lsd03, beliefPenaltyMPC_yy02, beliefPenaltyMPC_beta03, beliefPenaltyMPC_bmy03);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_20(beliefPenaltyMPC_Ld03, beliefPenaltyMPC_bmy03, beliefPenaltyMPC_yy03);
beliefPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_20_20(beliefPenaltyMPC_Ld03, beliefPenaltyMPC_Ysd04, beliefPenaltyMPC_Lsd04);
beliefPenaltyMPC_LA_DENSE_MMTSUB_20_20(beliefPenaltyMPC_Lsd04, beliefPenaltyMPC_Yd04);
beliefPenaltyMPC_LA_DENSE_CHOL_20(beliefPenaltyMPC_Yd04, beliefPenaltyMPC_Ld04);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_20_20(beliefPenaltyMPC_Lsd04, beliefPenaltyMPC_yy03, beliefPenaltyMPC_beta04, beliefPenaltyMPC_bmy04);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_20(beliefPenaltyMPC_Ld04, beliefPenaltyMPC_bmy04, beliefPenaltyMPC_yy04);
beliefPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_20_20(beliefPenaltyMPC_Ld04, beliefPenaltyMPC_Ysd05, beliefPenaltyMPC_Lsd05);
beliefPenaltyMPC_LA_DENSE_MMTSUB_20_20(beliefPenaltyMPC_Lsd05, beliefPenaltyMPC_Yd05);
beliefPenaltyMPC_LA_DENSE_CHOL_20(beliefPenaltyMPC_Yd05, beliefPenaltyMPC_Ld05);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_20_20(beliefPenaltyMPC_Lsd05, beliefPenaltyMPC_yy04, beliefPenaltyMPC_beta05, beliefPenaltyMPC_bmy05);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_20(beliefPenaltyMPC_Ld05, beliefPenaltyMPC_bmy05, beliefPenaltyMPC_yy05);
beliefPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_20_20(beliefPenaltyMPC_Ld05, beliefPenaltyMPC_Ysd06, beliefPenaltyMPC_Lsd06);
beliefPenaltyMPC_LA_DENSE_MMTSUB_20_20(beliefPenaltyMPC_Lsd06, beliefPenaltyMPC_Yd06);
beliefPenaltyMPC_LA_DENSE_CHOL_20(beliefPenaltyMPC_Yd06, beliefPenaltyMPC_Ld06);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_20_20(beliefPenaltyMPC_Lsd06, beliefPenaltyMPC_yy05, beliefPenaltyMPC_beta06, beliefPenaltyMPC_bmy06);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_20(beliefPenaltyMPC_Ld06, beliefPenaltyMPC_bmy06, beliefPenaltyMPC_yy06);
beliefPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_20_20(beliefPenaltyMPC_Ld06, beliefPenaltyMPC_Ysd07, beliefPenaltyMPC_Lsd07);
beliefPenaltyMPC_LA_DENSE_MMTSUB_20_20(beliefPenaltyMPC_Lsd07, beliefPenaltyMPC_Yd07);
beliefPenaltyMPC_LA_DENSE_CHOL_20(beliefPenaltyMPC_Yd07, beliefPenaltyMPC_Ld07);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_20_20(beliefPenaltyMPC_Lsd07, beliefPenaltyMPC_yy06, beliefPenaltyMPC_beta07, beliefPenaltyMPC_bmy07);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_20(beliefPenaltyMPC_Ld07, beliefPenaltyMPC_bmy07, beliefPenaltyMPC_yy07);
beliefPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_20_20(beliefPenaltyMPC_Ld07, beliefPenaltyMPC_Ysd08, beliefPenaltyMPC_Lsd08);
beliefPenaltyMPC_LA_DENSE_MMTSUB_20_20(beliefPenaltyMPC_Lsd08, beliefPenaltyMPC_Yd08);
beliefPenaltyMPC_LA_DENSE_CHOL_20(beliefPenaltyMPC_Yd08, beliefPenaltyMPC_Ld08);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_20_20(beliefPenaltyMPC_Lsd08, beliefPenaltyMPC_yy07, beliefPenaltyMPC_beta08, beliefPenaltyMPC_bmy08);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_20(beliefPenaltyMPC_Ld08, beliefPenaltyMPC_bmy08, beliefPenaltyMPC_yy08);
beliefPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_20_20(beliefPenaltyMPC_Ld08, beliefPenaltyMPC_Ysd09, beliefPenaltyMPC_Lsd09);
beliefPenaltyMPC_LA_DENSE_MMTSUB_20_20(beliefPenaltyMPC_Lsd09, beliefPenaltyMPC_Yd09);
beliefPenaltyMPC_LA_DENSE_CHOL_20(beliefPenaltyMPC_Yd09, beliefPenaltyMPC_Ld09);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_20_20(beliefPenaltyMPC_Lsd09, beliefPenaltyMPC_yy08, beliefPenaltyMPC_beta09, beliefPenaltyMPC_bmy09);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_20(beliefPenaltyMPC_Ld09, beliefPenaltyMPC_bmy09, beliefPenaltyMPC_yy09);
beliefPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_20_20(beliefPenaltyMPC_Ld09, beliefPenaltyMPC_Ysd10, beliefPenaltyMPC_Lsd10);
beliefPenaltyMPC_LA_DENSE_MMTSUB_20_20(beliefPenaltyMPC_Lsd10, beliefPenaltyMPC_Yd10);
beliefPenaltyMPC_LA_DENSE_CHOL_20(beliefPenaltyMPC_Yd10, beliefPenaltyMPC_Ld10);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_20_20(beliefPenaltyMPC_Lsd10, beliefPenaltyMPC_yy09, beliefPenaltyMPC_beta10, beliefPenaltyMPC_bmy10);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_20(beliefPenaltyMPC_Ld10, beliefPenaltyMPC_bmy10, beliefPenaltyMPC_yy10);
beliefPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_20_20(beliefPenaltyMPC_Ld10, beliefPenaltyMPC_Ysd11, beliefPenaltyMPC_Lsd11);
beliefPenaltyMPC_LA_DENSE_MMTSUB_20_20(beliefPenaltyMPC_Lsd11, beliefPenaltyMPC_Yd11);
beliefPenaltyMPC_LA_DENSE_CHOL_20(beliefPenaltyMPC_Yd11, beliefPenaltyMPC_Ld11);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_20_20(beliefPenaltyMPC_Lsd11, beliefPenaltyMPC_yy10, beliefPenaltyMPC_beta11, beliefPenaltyMPC_bmy11);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_20(beliefPenaltyMPC_Ld11, beliefPenaltyMPC_bmy11, beliefPenaltyMPC_yy11);
beliefPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_20_20(beliefPenaltyMPC_Ld11, beliefPenaltyMPC_Ysd12, beliefPenaltyMPC_Lsd12);
beliefPenaltyMPC_LA_DENSE_MMTSUB_20_20(beliefPenaltyMPC_Lsd12, beliefPenaltyMPC_Yd12);
beliefPenaltyMPC_LA_DENSE_CHOL_20(beliefPenaltyMPC_Yd12, beliefPenaltyMPC_Ld12);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_20_20(beliefPenaltyMPC_Lsd12, beliefPenaltyMPC_yy11, beliefPenaltyMPC_beta12, beliefPenaltyMPC_bmy12);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_20(beliefPenaltyMPC_Ld12, beliefPenaltyMPC_bmy12, beliefPenaltyMPC_yy12);
beliefPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_20_20(beliefPenaltyMPC_Ld12, beliefPenaltyMPC_Ysd13, beliefPenaltyMPC_Lsd13);
beliefPenaltyMPC_LA_DENSE_MMTSUB_20_20(beliefPenaltyMPC_Lsd13, beliefPenaltyMPC_Yd13);
beliefPenaltyMPC_LA_DENSE_CHOL_20(beliefPenaltyMPC_Yd13, beliefPenaltyMPC_Ld13);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_20_20(beliefPenaltyMPC_Lsd13, beliefPenaltyMPC_yy12, beliefPenaltyMPC_beta13, beliefPenaltyMPC_bmy13);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_20(beliefPenaltyMPC_Ld13, beliefPenaltyMPC_bmy13, beliefPenaltyMPC_yy13);
beliefPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_20_20(beliefPenaltyMPC_Ld13, beliefPenaltyMPC_Ysd14, beliefPenaltyMPC_Lsd14);
beliefPenaltyMPC_LA_DENSE_MMTSUB_20_20(beliefPenaltyMPC_Lsd14, beliefPenaltyMPC_Yd14);
beliefPenaltyMPC_LA_DENSE_CHOL_20(beliefPenaltyMPC_Yd14, beliefPenaltyMPC_Ld14);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_20_20(beliefPenaltyMPC_Lsd14, beliefPenaltyMPC_yy13, beliefPenaltyMPC_beta14, beliefPenaltyMPC_bmy14);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_20(beliefPenaltyMPC_Ld14, beliefPenaltyMPC_bmy14, beliefPenaltyMPC_yy14);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_20(beliefPenaltyMPC_Ld14, beliefPenaltyMPC_yy14, beliefPenaltyMPC_dvaff14);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_20_20(beliefPenaltyMPC_Lsd14, beliefPenaltyMPC_dvaff14, beliefPenaltyMPC_yy13, beliefPenaltyMPC_bmy13);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_20(beliefPenaltyMPC_Ld13, beliefPenaltyMPC_bmy13, beliefPenaltyMPC_dvaff13);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_20_20(beliefPenaltyMPC_Lsd13, beliefPenaltyMPC_dvaff13, beliefPenaltyMPC_yy12, beliefPenaltyMPC_bmy12);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_20(beliefPenaltyMPC_Ld12, beliefPenaltyMPC_bmy12, beliefPenaltyMPC_dvaff12);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_20_20(beliefPenaltyMPC_Lsd12, beliefPenaltyMPC_dvaff12, beliefPenaltyMPC_yy11, beliefPenaltyMPC_bmy11);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_20(beliefPenaltyMPC_Ld11, beliefPenaltyMPC_bmy11, beliefPenaltyMPC_dvaff11);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_20_20(beliefPenaltyMPC_Lsd11, beliefPenaltyMPC_dvaff11, beliefPenaltyMPC_yy10, beliefPenaltyMPC_bmy10);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_20(beliefPenaltyMPC_Ld10, beliefPenaltyMPC_bmy10, beliefPenaltyMPC_dvaff10);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_20_20(beliefPenaltyMPC_Lsd10, beliefPenaltyMPC_dvaff10, beliefPenaltyMPC_yy09, beliefPenaltyMPC_bmy09);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_20(beliefPenaltyMPC_Ld09, beliefPenaltyMPC_bmy09, beliefPenaltyMPC_dvaff09);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_20_20(beliefPenaltyMPC_Lsd09, beliefPenaltyMPC_dvaff09, beliefPenaltyMPC_yy08, beliefPenaltyMPC_bmy08);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_20(beliefPenaltyMPC_Ld08, beliefPenaltyMPC_bmy08, beliefPenaltyMPC_dvaff08);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_20_20(beliefPenaltyMPC_Lsd08, beliefPenaltyMPC_dvaff08, beliefPenaltyMPC_yy07, beliefPenaltyMPC_bmy07);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_20(beliefPenaltyMPC_Ld07, beliefPenaltyMPC_bmy07, beliefPenaltyMPC_dvaff07);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_20_20(beliefPenaltyMPC_Lsd07, beliefPenaltyMPC_dvaff07, beliefPenaltyMPC_yy06, beliefPenaltyMPC_bmy06);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_20(beliefPenaltyMPC_Ld06, beliefPenaltyMPC_bmy06, beliefPenaltyMPC_dvaff06);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_20_20(beliefPenaltyMPC_Lsd06, beliefPenaltyMPC_dvaff06, beliefPenaltyMPC_yy05, beliefPenaltyMPC_bmy05);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_20(beliefPenaltyMPC_Ld05, beliefPenaltyMPC_bmy05, beliefPenaltyMPC_dvaff05);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_20_20(beliefPenaltyMPC_Lsd05, beliefPenaltyMPC_dvaff05, beliefPenaltyMPC_yy04, beliefPenaltyMPC_bmy04);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_20(beliefPenaltyMPC_Ld04, beliefPenaltyMPC_bmy04, beliefPenaltyMPC_dvaff04);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_20_20(beliefPenaltyMPC_Lsd04, beliefPenaltyMPC_dvaff04, beliefPenaltyMPC_yy03, beliefPenaltyMPC_bmy03);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_20(beliefPenaltyMPC_Ld03, beliefPenaltyMPC_bmy03, beliefPenaltyMPC_dvaff03);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_20_20(beliefPenaltyMPC_Lsd03, beliefPenaltyMPC_dvaff03, beliefPenaltyMPC_yy02, beliefPenaltyMPC_bmy02);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_20(beliefPenaltyMPC_Ld02, beliefPenaltyMPC_bmy02, beliefPenaltyMPC_dvaff02);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_20_20(beliefPenaltyMPC_Lsd02, beliefPenaltyMPC_dvaff02, beliefPenaltyMPC_yy01, beliefPenaltyMPC_bmy01);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_20(beliefPenaltyMPC_Ld01, beliefPenaltyMPC_bmy01, beliefPenaltyMPC_dvaff01);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_20_20(beliefPenaltyMPC_Lsd01, beliefPenaltyMPC_dvaff01, beliefPenaltyMPC_yy00, beliefPenaltyMPC_bmy00);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_20(beliefPenaltyMPC_Ld00, beliefPenaltyMPC_bmy00, beliefPenaltyMPC_dvaff00);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_20_62_20(params->C1, beliefPenaltyMPC_dvaff01, beliefPenaltyMPC_D00, beliefPenaltyMPC_dvaff00, beliefPenaltyMPC_grad_eq00);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_20_62_20(params->C2, beliefPenaltyMPC_dvaff02, beliefPenaltyMPC_D01, beliefPenaltyMPC_dvaff01, beliefPenaltyMPC_grad_eq01);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_20_62_20(params->C3, beliefPenaltyMPC_dvaff03, beliefPenaltyMPC_D01, beliefPenaltyMPC_dvaff02, beliefPenaltyMPC_grad_eq02);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_20_62_20(params->C4, beliefPenaltyMPC_dvaff04, beliefPenaltyMPC_D01, beliefPenaltyMPC_dvaff03, beliefPenaltyMPC_grad_eq03);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_20_62_20(params->C5, beliefPenaltyMPC_dvaff05, beliefPenaltyMPC_D01, beliefPenaltyMPC_dvaff04, beliefPenaltyMPC_grad_eq04);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_20_62_20(params->C6, beliefPenaltyMPC_dvaff06, beliefPenaltyMPC_D01, beliefPenaltyMPC_dvaff05, beliefPenaltyMPC_grad_eq05);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_20_62_20(params->C7, beliefPenaltyMPC_dvaff07, beliefPenaltyMPC_D01, beliefPenaltyMPC_dvaff06, beliefPenaltyMPC_grad_eq06);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_20_62_20(params->C8, beliefPenaltyMPC_dvaff08, beliefPenaltyMPC_D01, beliefPenaltyMPC_dvaff07, beliefPenaltyMPC_grad_eq07);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_20_62_20(params->C9, beliefPenaltyMPC_dvaff09, beliefPenaltyMPC_D01, beliefPenaltyMPC_dvaff08, beliefPenaltyMPC_grad_eq08);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_20_62_20(params->C10, beliefPenaltyMPC_dvaff10, beliefPenaltyMPC_D01, beliefPenaltyMPC_dvaff09, beliefPenaltyMPC_grad_eq09);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_20_62_20(params->C11, beliefPenaltyMPC_dvaff11, beliefPenaltyMPC_D01, beliefPenaltyMPC_dvaff10, beliefPenaltyMPC_grad_eq10);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_20_62_20(params->C12, beliefPenaltyMPC_dvaff12, beliefPenaltyMPC_D01, beliefPenaltyMPC_dvaff11, beliefPenaltyMPC_grad_eq11);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_20_62_20(params->C13, beliefPenaltyMPC_dvaff13, beliefPenaltyMPC_D01, beliefPenaltyMPC_dvaff12, beliefPenaltyMPC_grad_eq12);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_20_62_20(params->C14, beliefPenaltyMPC_dvaff14, beliefPenaltyMPC_D01, beliefPenaltyMPC_dvaff13, beliefPenaltyMPC_grad_eq13);
beliefPenaltyMPC_LA_DIAGZERO_MTVM_20_20(beliefPenaltyMPC_D14, beliefPenaltyMPC_dvaff14, beliefPenaltyMPC_grad_eq14);
beliefPenaltyMPC_LA_VSUB2_888(beliefPenaltyMPC_rd, beliefPenaltyMPC_grad_eq, beliefPenaltyMPC_rd);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_62(beliefPenaltyMPC_Phi00, beliefPenaltyMPC_rd00, beliefPenaltyMPC_dzaff00);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_62(beliefPenaltyMPC_Phi01, beliefPenaltyMPC_rd01, beliefPenaltyMPC_dzaff01);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_62(beliefPenaltyMPC_Phi02, beliefPenaltyMPC_rd02, beliefPenaltyMPC_dzaff02);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_62(beliefPenaltyMPC_Phi03, beliefPenaltyMPC_rd03, beliefPenaltyMPC_dzaff03);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_62(beliefPenaltyMPC_Phi04, beliefPenaltyMPC_rd04, beliefPenaltyMPC_dzaff04);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_62(beliefPenaltyMPC_Phi05, beliefPenaltyMPC_rd05, beliefPenaltyMPC_dzaff05);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_62(beliefPenaltyMPC_Phi06, beliefPenaltyMPC_rd06, beliefPenaltyMPC_dzaff06);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_62(beliefPenaltyMPC_Phi07, beliefPenaltyMPC_rd07, beliefPenaltyMPC_dzaff07);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_62(beliefPenaltyMPC_Phi08, beliefPenaltyMPC_rd08, beliefPenaltyMPC_dzaff08);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_62(beliefPenaltyMPC_Phi09, beliefPenaltyMPC_rd09, beliefPenaltyMPC_dzaff09);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_62(beliefPenaltyMPC_Phi10, beliefPenaltyMPC_rd10, beliefPenaltyMPC_dzaff10);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_62(beliefPenaltyMPC_Phi11, beliefPenaltyMPC_rd11, beliefPenaltyMPC_dzaff11);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_62(beliefPenaltyMPC_Phi12, beliefPenaltyMPC_rd12, beliefPenaltyMPC_dzaff12);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_62(beliefPenaltyMPC_Phi13, beliefPenaltyMPC_rd13, beliefPenaltyMPC_dzaff13);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_20(beliefPenaltyMPC_Phi14, beliefPenaltyMPC_rd14, beliefPenaltyMPC_dzaff14);
beliefPenaltyMPC_LA_VSUB_INDEXED_62(beliefPenaltyMPC_dzaff00, beliefPenaltyMPC_lbIdx00, beliefPenaltyMPC_rilb00, beliefPenaltyMPC_dslbaff00);
beliefPenaltyMPC_LA_VSUB3_62(beliefPenaltyMPC_llbbyslb00, beliefPenaltyMPC_dslbaff00, beliefPenaltyMPC_llb00, beliefPenaltyMPC_dllbaff00);
beliefPenaltyMPC_LA_VSUB2_INDEXED_22(beliefPenaltyMPC_riub00, beliefPenaltyMPC_dzaff00, beliefPenaltyMPC_ubIdx00, beliefPenaltyMPC_dsubaff00);
beliefPenaltyMPC_LA_VSUB3_22(beliefPenaltyMPC_lubbysub00, beliefPenaltyMPC_dsubaff00, beliefPenaltyMPC_lub00, beliefPenaltyMPC_dlubaff00);
beliefPenaltyMPC_LA_VSUB_INDEXED_62(beliefPenaltyMPC_dzaff01, beliefPenaltyMPC_lbIdx01, beliefPenaltyMPC_rilb01, beliefPenaltyMPC_dslbaff01);
beliefPenaltyMPC_LA_VSUB3_62(beliefPenaltyMPC_llbbyslb01, beliefPenaltyMPC_dslbaff01, beliefPenaltyMPC_llb01, beliefPenaltyMPC_dllbaff01);
beliefPenaltyMPC_LA_VSUB2_INDEXED_22(beliefPenaltyMPC_riub01, beliefPenaltyMPC_dzaff01, beliefPenaltyMPC_ubIdx01, beliefPenaltyMPC_dsubaff01);
beliefPenaltyMPC_LA_VSUB3_22(beliefPenaltyMPC_lubbysub01, beliefPenaltyMPC_dsubaff01, beliefPenaltyMPC_lub01, beliefPenaltyMPC_dlubaff01);
beliefPenaltyMPC_LA_VSUB_INDEXED_62(beliefPenaltyMPC_dzaff02, beliefPenaltyMPC_lbIdx02, beliefPenaltyMPC_rilb02, beliefPenaltyMPC_dslbaff02);
beliefPenaltyMPC_LA_VSUB3_62(beliefPenaltyMPC_llbbyslb02, beliefPenaltyMPC_dslbaff02, beliefPenaltyMPC_llb02, beliefPenaltyMPC_dllbaff02);
beliefPenaltyMPC_LA_VSUB2_INDEXED_22(beliefPenaltyMPC_riub02, beliefPenaltyMPC_dzaff02, beliefPenaltyMPC_ubIdx02, beliefPenaltyMPC_dsubaff02);
beliefPenaltyMPC_LA_VSUB3_22(beliefPenaltyMPC_lubbysub02, beliefPenaltyMPC_dsubaff02, beliefPenaltyMPC_lub02, beliefPenaltyMPC_dlubaff02);
beliefPenaltyMPC_LA_VSUB_INDEXED_62(beliefPenaltyMPC_dzaff03, beliefPenaltyMPC_lbIdx03, beliefPenaltyMPC_rilb03, beliefPenaltyMPC_dslbaff03);
beliefPenaltyMPC_LA_VSUB3_62(beliefPenaltyMPC_llbbyslb03, beliefPenaltyMPC_dslbaff03, beliefPenaltyMPC_llb03, beliefPenaltyMPC_dllbaff03);
beliefPenaltyMPC_LA_VSUB2_INDEXED_22(beliefPenaltyMPC_riub03, beliefPenaltyMPC_dzaff03, beliefPenaltyMPC_ubIdx03, beliefPenaltyMPC_dsubaff03);
beliefPenaltyMPC_LA_VSUB3_22(beliefPenaltyMPC_lubbysub03, beliefPenaltyMPC_dsubaff03, beliefPenaltyMPC_lub03, beliefPenaltyMPC_dlubaff03);
beliefPenaltyMPC_LA_VSUB_INDEXED_62(beliefPenaltyMPC_dzaff04, beliefPenaltyMPC_lbIdx04, beliefPenaltyMPC_rilb04, beliefPenaltyMPC_dslbaff04);
beliefPenaltyMPC_LA_VSUB3_62(beliefPenaltyMPC_llbbyslb04, beliefPenaltyMPC_dslbaff04, beliefPenaltyMPC_llb04, beliefPenaltyMPC_dllbaff04);
beliefPenaltyMPC_LA_VSUB2_INDEXED_22(beliefPenaltyMPC_riub04, beliefPenaltyMPC_dzaff04, beliefPenaltyMPC_ubIdx04, beliefPenaltyMPC_dsubaff04);
beliefPenaltyMPC_LA_VSUB3_22(beliefPenaltyMPC_lubbysub04, beliefPenaltyMPC_dsubaff04, beliefPenaltyMPC_lub04, beliefPenaltyMPC_dlubaff04);
beliefPenaltyMPC_LA_VSUB_INDEXED_62(beliefPenaltyMPC_dzaff05, beliefPenaltyMPC_lbIdx05, beliefPenaltyMPC_rilb05, beliefPenaltyMPC_dslbaff05);
beliefPenaltyMPC_LA_VSUB3_62(beliefPenaltyMPC_llbbyslb05, beliefPenaltyMPC_dslbaff05, beliefPenaltyMPC_llb05, beliefPenaltyMPC_dllbaff05);
beliefPenaltyMPC_LA_VSUB2_INDEXED_22(beliefPenaltyMPC_riub05, beliefPenaltyMPC_dzaff05, beliefPenaltyMPC_ubIdx05, beliefPenaltyMPC_dsubaff05);
beliefPenaltyMPC_LA_VSUB3_22(beliefPenaltyMPC_lubbysub05, beliefPenaltyMPC_dsubaff05, beliefPenaltyMPC_lub05, beliefPenaltyMPC_dlubaff05);
beliefPenaltyMPC_LA_VSUB_INDEXED_62(beliefPenaltyMPC_dzaff06, beliefPenaltyMPC_lbIdx06, beliefPenaltyMPC_rilb06, beliefPenaltyMPC_dslbaff06);
beliefPenaltyMPC_LA_VSUB3_62(beliefPenaltyMPC_llbbyslb06, beliefPenaltyMPC_dslbaff06, beliefPenaltyMPC_llb06, beliefPenaltyMPC_dllbaff06);
beliefPenaltyMPC_LA_VSUB2_INDEXED_22(beliefPenaltyMPC_riub06, beliefPenaltyMPC_dzaff06, beliefPenaltyMPC_ubIdx06, beliefPenaltyMPC_dsubaff06);
beliefPenaltyMPC_LA_VSUB3_22(beliefPenaltyMPC_lubbysub06, beliefPenaltyMPC_dsubaff06, beliefPenaltyMPC_lub06, beliefPenaltyMPC_dlubaff06);
beliefPenaltyMPC_LA_VSUB_INDEXED_62(beliefPenaltyMPC_dzaff07, beliefPenaltyMPC_lbIdx07, beliefPenaltyMPC_rilb07, beliefPenaltyMPC_dslbaff07);
beliefPenaltyMPC_LA_VSUB3_62(beliefPenaltyMPC_llbbyslb07, beliefPenaltyMPC_dslbaff07, beliefPenaltyMPC_llb07, beliefPenaltyMPC_dllbaff07);
beliefPenaltyMPC_LA_VSUB2_INDEXED_22(beliefPenaltyMPC_riub07, beliefPenaltyMPC_dzaff07, beliefPenaltyMPC_ubIdx07, beliefPenaltyMPC_dsubaff07);
beliefPenaltyMPC_LA_VSUB3_22(beliefPenaltyMPC_lubbysub07, beliefPenaltyMPC_dsubaff07, beliefPenaltyMPC_lub07, beliefPenaltyMPC_dlubaff07);
beliefPenaltyMPC_LA_VSUB_INDEXED_62(beliefPenaltyMPC_dzaff08, beliefPenaltyMPC_lbIdx08, beliefPenaltyMPC_rilb08, beliefPenaltyMPC_dslbaff08);
beliefPenaltyMPC_LA_VSUB3_62(beliefPenaltyMPC_llbbyslb08, beliefPenaltyMPC_dslbaff08, beliefPenaltyMPC_llb08, beliefPenaltyMPC_dllbaff08);
beliefPenaltyMPC_LA_VSUB2_INDEXED_22(beliefPenaltyMPC_riub08, beliefPenaltyMPC_dzaff08, beliefPenaltyMPC_ubIdx08, beliefPenaltyMPC_dsubaff08);
beliefPenaltyMPC_LA_VSUB3_22(beliefPenaltyMPC_lubbysub08, beliefPenaltyMPC_dsubaff08, beliefPenaltyMPC_lub08, beliefPenaltyMPC_dlubaff08);
beliefPenaltyMPC_LA_VSUB_INDEXED_62(beliefPenaltyMPC_dzaff09, beliefPenaltyMPC_lbIdx09, beliefPenaltyMPC_rilb09, beliefPenaltyMPC_dslbaff09);
beliefPenaltyMPC_LA_VSUB3_62(beliefPenaltyMPC_llbbyslb09, beliefPenaltyMPC_dslbaff09, beliefPenaltyMPC_llb09, beliefPenaltyMPC_dllbaff09);
beliefPenaltyMPC_LA_VSUB2_INDEXED_22(beliefPenaltyMPC_riub09, beliefPenaltyMPC_dzaff09, beliefPenaltyMPC_ubIdx09, beliefPenaltyMPC_dsubaff09);
beliefPenaltyMPC_LA_VSUB3_22(beliefPenaltyMPC_lubbysub09, beliefPenaltyMPC_dsubaff09, beliefPenaltyMPC_lub09, beliefPenaltyMPC_dlubaff09);
beliefPenaltyMPC_LA_VSUB_INDEXED_62(beliefPenaltyMPC_dzaff10, beliefPenaltyMPC_lbIdx10, beliefPenaltyMPC_rilb10, beliefPenaltyMPC_dslbaff10);
beliefPenaltyMPC_LA_VSUB3_62(beliefPenaltyMPC_llbbyslb10, beliefPenaltyMPC_dslbaff10, beliefPenaltyMPC_llb10, beliefPenaltyMPC_dllbaff10);
beliefPenaltyMPC_LA_VSUB2_INDEXED_22(beliefPenaltyMPC_riub10, beliefPenaltyMPC_dzaff10, beliefPenaltyMPC_ubIdx10, beliefPenaltyMPC_dsubaff10);
beliefPenaltyMPC_LA_VSUB3_22(beliefPenaltyMPC_lubbysub10, beliefPenaltyMPC_dsubaff10, beliefPenaltyMPC_lub10, beliefPenaltyMPC_dlubaff10);
beliefPenaltyMPC_LA_VSUB_INDEXED_62(beliefPenaltyMPC_dzaff11, beliefPenaltyMPC_lbIdx11, beliefPenaltyMPC_rilb11, beliefPenaltyMPC_dslbaff11);
beliefPenaltyMPC_LA_VSUB3_62(beliefPenaltyMPC_llbbyslb11, beliefPenaltyMPC_dslbaff11, beliefPenaltyMPC_llb11, beliefPenaltyMPC_dllbaff11);
beliefPenaltyMPC_LA_VSUB2_INDEXED_22(beliefPenaltyMPC_riub11, beliefPenaltyMPC_dzaff11, beliefPenaltyMPC_ubIdx11, beliefPenaltyMPC_dsubaff11);
beliefPenaltyMPC_LA_VSUB3_22(beliefPenaltyMPC_lubbysub11, beliefPenaltyMPC_dsubaff11, beliefPenaltyMPC_lub11, beliefPenaltyMPC_dlubaff11);
beliefPenaltyMPC_LA_VSUB_INDEXED_62(beliefPenaltyMPC_dzaff12, beliefPenaltyMPC_lbIdx12, beliefPenaltyMPC_rilb12, beliefPenaltyMPC_dslbaff12);
beliefPenaltyMPC_LA_VSUB3_62(beliefPenaltyMPC_llbbyslb12, beliefPenaltyMPC_dslbaff12, beliefPenaltyMPC_llb12, beliefPenaltyMPC_dllbaff12);
beliefPenaltyMPC_LA_VSUB2_INDEXED_22(beliefPenaltyMPC_riub12, beliefPenaltyMPC_dzaff12, beliefPenaltyMPC_ubIdx12, beliefPenaltyMPC_dsubaff12);
beliefPenaltyMPC_LA_VSUB3_22(beliefPenaltyMPC_lubbysub12, beliefPenaltyMPC_dsubaff12, beliefPenaltyMPC_lub12, beliefPenaltyMPC_dlubaff12);
beliefPenaltyMPC_LA_VSUB_INDEXED_62(beliefPenaltyMPC_dzaff13, beliefPenaltyMPC_lbIdx13, beliefPenaltyMPC_rilb13, beliefPenaltyMPC_dslbaff13);
beliefPenaltyMPC_LA_VSUB3_62(beliefPenaltyMPC_llbbyslb13, beliefPenaltyMPC_dslbaff13, beliefPenaltyMPC_llb13, beliefPenaltyMPC_dllbaff13);
beliefPenaltyMPC_LA_VSUB2_INDEXED_22(beliefPenaltyMPC_riub13, beliefPenaltyMPC_dzaff13, beliefPenaltyMPC_ubIdx13, beliefPenaltyMPC_dsubaff13);
beliefPenaltyMPC_LA_VSUB3_22(beliefPenaltyMPC_lubbysub13, beliefPenaltyMPC_dsubaff13, beliefPenaltyMPC_lub13, beliefPenaltyMPC_dlubaff13);
beliefPenaltyMPC_LA_VSUB_INDEXED_20(beliefPenaltyMPC_dzaff14, beliefPenaltyMPC_lbIdx14, beliefPenaltyMPC_rilb14, beliefPenaltyMPC_dslbaff14);
beliefPenaltyMPC_LA_VSUB3_20(beliefPenaltyMPC_llbbyslb14, beliefPenaltyMPC_dslbaff14, beliefPenaltyMPC_llb14, beliefPenaltyMPC_dllbaff14);
beliefPenaltyMPC_LA_VSUB2_INDEXED_20(beliefPenaltyMPC_riub14, beliefPenaltyMPC_dzaff14, beliefPenaltyMPC_ubIdx14, beliefPenaltyMPC_dsubaff14);
beliefPenaltyMPC_LA_VSUB3_20(beliefPenaltyMPC_lubbysub14, beliefPenaltyMPC_dsubaff14, beliefPenaltyMPC_lub14, beliefPenaltyMPC_dlubaff14);
info->lsit_aff = beliefPenaltyMPC_LINESEARCH_BACKTRACKING_AFFINE(beliefPenaltyMPC_l, beliefPenaltyMPC_s, beliefPenaltyMPC_dl_aff, beliefPenaltyMPC_ds_aff, &info->step_aff, &info->mu_aff);
if( info->lsit_aff == beliefPenaltyMPC_NOPROGRESS ){
PRINTTEXT("Affine line search could not proceed at iteration %d.\nThe problem might be infeasible -- exiting.\n",info->it+1);
exitcode = beliefPenaltyMPC_NOPROGRESS; break;
}
sigma_3rdroot = info->mu_aff / info->mu;
info->sigma = sigma_3rdroot*sigma_3rdroot*sigma_3rdroot;
musigma = info->mu * info->sigma;
beliefPenaltyMPC_LA_VSUB5_1216(beliefPenaltyMPC_ds_aff, beliefPenaltyMPC_dl_aff, musigma, beliefPenaltyMPC_ccrhs);
beliefPenaltyMPC_LA_VSUB6_INDEXED_62_22_62(beliefPenaltyMPC_ccrhsub00, beliefPenaltyMPC_sub00, beliefPenaltyMPC_ubIdx00, beliefPenaltyMPC_ccrhsl00, beliefPenaltyMPC_slb00, beliefPenaltyMPC_lbIdx00, beliefPenaltyMPC_rd00);
beliefPenaltyMPC_LA_VSUB6_INDEXED_62_22_62(beliefPenaltyMPC_ccrhsub01, beliefPenaltyMPC_sub01, beliefPenaltyMPC_ubIdx01, beliefPenaltyMPC_ccrhsl01, beliefPenaltyMPC_slb01, beliefPenaltyMPC_lbIdx01, beliefPenaltyMPC_rd01);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_62(beliefPenaltyMPC_Phi00, beliefPenaltyMPC_rd00, beliefPenaltyMPC_Lbyrd00);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_62(beliefPenaltyMPC_Phi01, beliefPenaltyMPC_rd01, beliefPenaltyMPC_Lbyrd01);
beliefPenaltyMPC_LA_DIAGZERO_MVM_20(beliefPenaltyMPC_W00, beliefPenaltyMPC_Lbyrd00, beliefPenaltyMPC_beta00);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_20(beliefPenaltyMPC_Ld00, beliefPenaltyMPC_beta00, beliefPenaltyMPC_yy00);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_20_62_62(beliefPenaltyMPC_V00, beliefPenaltyMPC_Lbyrd00, beliefPenaltyMPC_W01, beliefPenaltyMPC_Lbyrd01, beliefPenaltyMPC_beta01);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_20_20(beliefPenaltyMPC_Lsd01, beliefPenaltyMPC_yy00, beliefPenaltyMPC_beta01, beliefPenaltyMPC_bmy01);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_20(beliefPenaltyMPC_Ld01, beliefPenaltyMPC_bmy01, beliefPenaltyMPC_yy01);
beliefPenaltyMPC_LA_VSUB6_INDEXED_62_22_62(beliefPenaltyMPC_ccrhsub02, beliefPenaltyMPC_sub02, beliefPenaltyMPC_ubIdx02, beliefPenaltyMPC_ccrhsl02, beliefPenaltyMPC_slb02, beliefPenaltyMPC_lbIdx02, beliefPenaltyMPC_rd02);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_62(beliefPenaltyMPC_Phi02, beliefPenaltyMPC_rd02, beliefPenaltyMPC_Lbyrd02);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_20_62_62(beliefPenaltyMPC_V01, beliefPenaltyMPC_Lbyrd01, beliefPenaltyMPC_W02, beliefPenaltyMPC_Lbyrd02, beliefPenaltyMPC_beta02);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_20_20(beliefPenaltyMPC_Lsd02, beliefPenaltyMPC_yy01, beliefPenaltyMPC_beta02, beliefPenaltyMPC_bmy02);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_20(beliefPenaltyMPC_Ld02, beliefPenaltyMPC_bmy02, beliefPenaltyMPC_yy02);
beliefPenaltyMPC_LA_VSUB6_INDEXED_62_22_62(beliefPenaltyMPC_ccrhsub03, beliefPenaltyMPC_sub03, beliefPenaltyMPC_ubIdx03, beliefPenaltyMPC_ccrhsl03, beliefPenaltyMPC_slb03, beliefPenaltyMPC_lbIdx03, beliefPenaltyMPC_rd03);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_62(beliefPenaltyMPC_Phi03, beliefPenaltyMPC_rd03, beliefPenaltyMPC_Lbyrd03);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_20_62_62(beliefPenaltyMPC_V02, beliefPenaltyMPC_Lbyrd02, beliefPenaltyMPC_W03, beliefPenaltyMPC_Lbyrd03, beliefPenaltyMPC_beta03);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_20_20(beliefPenaltyMPC_Lsd03, beliefPenaltyMPC_yy02, beliefPenaltyMPC_beta03, beliefPenaltyMPC_bmy03);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_20(beliefPenaltyMPC_Ld03, beliefPenaltyMPC_bmy03, beliefPenaltyMPC_yy03);
beliefPenaltyMPC_LA_VSUB6_INDEXED_62_22_62(beliefPenaltyMPC_ccrhsub04, beliefPenaltyMPC_sub04, beliefPenaltyMPC_ubIdx04, beliefPenaltyMPC_ccrhsl04, beliefPenaltyMPC_slb04, beliefPenaltyMPC_lbIdx04, beliefPenaltyMPC_rd04);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_62(beliefPenaltyMPC_Phi04, beliefPenaltyMPC_rd04, beliefPenaltyMPC_Lbyrd04);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_20_62_62(beliefPenaltyMPC_V03, beliefPenaltyMPC_Lbyrd03, beliefPenaltyMPC_W04, beliefPenaltyMPC_Lbyrd04, beliefPenaltyMPC_beta04);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_20_20(beliefPenaltyMPC_Lsd04, beliefPenaltyMPC_yy03, beliefPenaltyMPC_beta04, beliefPenaltyMPC_bmy04);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_20(beliefPenaltyMPC_Ld04, beliefPenaltyMPC_bmy04, beliefPenaltyMPC_yy04);
beliefPenaltyMPC_LA_VSUB6_INDEXED_62_22_62(beliefPenaltyMPC_ccrhsub05, beliefPenaltyMPC_sub05, beliefPenaltyMPC_ubIdx05, beliefPenaltyMPC_ccrhsl05, beliefPenaltyMPC_slb05, beliefPenaltyMPC_lbIdx05, beliefPenaltyMPC_rd05);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_62(beliefPenaltyMPC_Phi05, beliefPenaltyMPC_rd05, beliefPenaltyMPC_Lbyrd05);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_20_62_62(beliefPenaltyMPC_V04, beliefPenaltyMPC_Lbyrd04, beliefPenaltyMPC_W05, beliefPenaltyMPC_Lbyrd05, beliefPenaltyMPC_beta05);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_20_20(beliefPenaltyMPC_Lsd05, beliefPenaltyMPC_yy04, beliefPenaltyMPC_beta05, beliefPenaltyMPC_bmy05);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_20(beliefPenaltyMPC_Ld05, beliefPenaltyMPC_bmy05, beliefPenaltyMPC_yy05);
beliefPenaltyMPC_LA_VSUB6_INDEXED_62_22_62(beliefPenaltyMPC_ccrhsub06, beliefPenaltyMPC_sub06, beliefPenaltyMPC_ubIdx06, beliefPenaltyMPC_ccrhsl06, beliefPenaltyMPC_slb06, beliefPenaltyMPC_lbIdx06, beliefPenaltyMPC_rd06);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_62(beliefPenaltyMPC_Phi06, beliefPenaltyMPC_rd06, beliefPenaltyMPC_Lbyrd06);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_20_62_62(beliefPenaltyMPC_V05, beliefPenaltyMPC_Lbyrd05, beliefPenaltyMPC_W06, beliefPenaltyMPC_Lbyrd06, beliefPenaltyMPC_beta06);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_20_20(beliefPenaltyMPC_Lsd06, beliefPenaltyMPC_yy05, beliefPenaltyMPC_beta06, beliefPenaltyMPC_bmy06);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_20(beliefPenaltyMPC_Ld06, beliefPenaltyMPC_bmy06, beliefPenaltyMPC_yy06);
beliefPenaltyMPC_LA_VSUB6_INDEXED_62_22_62(beliefPenaltyMPC_ccrhsub07, beliefPenaltyMPC_sub07, beliefPenaltyMPC_ubIdx07, beliefPenaltyMPC_ccrhsl07, beliefPenaltyMPC_slb07, beliefPenaltyMPC_lbIdx07, beliefPenaltyMPC_rd07);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_62(beliefPenaltyMPC_Phi07, beliefPenaltyMPC_rd07, beliefPenaltyMPC_Lbyrd07);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_20_62_62(beliefPenaltyMPC_V06, beliefPenaltyMPC_Lbyrd06, beliefPenaltyMPC_W07, beliefPenaltyMPC_Lbyrd07, beliefPenaltyMPC_beta07);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_20_20(beliefPenaltyMPC_Lsd07, beliefPenaltyMPC_yy06, beliefPenaltyMPC_beta07, beliefPenaltyMPC_bmy07);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_20(beliefPenaltyMPC_Ld07, beliefPenaltyMPC_bmy07, beliefPenaltyMPC_yy07);
beliefPenaltyMPC_LA_VSUB6_INDEXED_62_22_62(beliefPenaltyMPC_ccrhsub08, beliefPenaltyMPC_sub08, beliefPenaltyMPC_ubIdx08, beliefPenaltyMPC_ccrhsl08, beliefPenaltyMPC_slb08, beliefPenaltyMPC_lbIdx08, beliefPenaltyMPC_rd08);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_62(beliefPenaltyMPC_Phi08, beliefPenaltyMPC_rd08, beliefPenaltyMPC_Lbyrd08);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_20_62_62(beliefPenaltyMPC_V07, beliefPenaltyMPC_Lbyrd07, beliefPenaltyMPC_W08, beliefPenaltyMPC_Lbyrd08, beliefPenaltyMPC_beta08);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_20_20(beliefPenaltyMPC_Lsd08, beliefPenaltyMPC_yy07, beliefPenaltyMPC_beta08, beliefPenaltyMPC_bmy08);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_20(beliefPenaltyMPC_Ld08, beliefPenaltyMPC_bmy08, beliefPenaltyMPC_yy08);
beliefPenaltyMPC_LA_VSUB6_INDEXED_62_22_62(beliefPenaltyMPC_ccrhsub09, beliefPenaltyMPC_sub09, beliefPenaltyMPC_ubIdx09, beliefPenaltyMPC_ccrhsl09, beliefPenaltyMPC_slb09, beliefPenaltyMPC_lbIdx09, beliefPenaltyMPC_rd09);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_62(beliefPenaltyMPC_Phi09, beliefPenaltyMPC_rd09, beliefPenaltyMPC_Lbyrd09);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_20_62_62(beliefPenaltyMPC_V08, beliefPenaltyMPC_Lbyrd08, beliefPenaltyMPC_W09, beliefPenaltyMPC_Lbyrd09, beliefPenaltyMPC_beta09);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_20_20(beliefPenaltyMPC_Lsd09, beliefPenaltyMPC_yy08, beliefPenaltyMPC_beta09, beliefPenaltyMPC_bmy09);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_20(beliefPenaltyMPC_Ld09, beliefPenaltyMPC_bmy09, beliefPenaltyMPC_yy09);
beliefPenaltyMPC_LA_VSUB6_INDEXED_62_22_62(beliefPenaltyMPC_ccrhsub10, beliefPenaltyMPC_sub10, beliefPenaltyMPC_ubIdx10, beliefPenaltyMPC_ccrhsl10, beliefPenaltyMPC_slb10, beliefPenaltyMPC_lbIdx10, beliefPenaltyMPC_rd10);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_62(beliefPenaltyMPC_Phi10, beliefPenaltyMPC_rd10, beliefPenaltyMPC_Lbyrd10);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_20_62_62(beliefPenaltyMPC_V09, beliefPenaltyMPC_Lbyrd09, beliefPenaltyMPC_W10, beliefPenaltyMPC_Lbyrd10, beliefPenaltyMPC_beta10);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_20_20(beliefPenaltyMPC_Lsd10, beliefPenaltyMPC_yy09, beliefPenaltyMPC_beta10, beliefPenaltyMPC_bmy10);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_20(beliefPenaltyMPC_Ld10, beliefPenaltyMPC_bmy10, beliefPenaltyMPC_yy10);
beliefPenaltyMPC_LA_VSUB6_INDEXED_62_22_62(beliefPenaltyMPC_ccrhsub11, beliefPenaltyMPC_sub11, beliefPenaltyMPC_ubIdx11, beliefPenaltyMPC_ccrhsl11, beliefPenaltyMPC_slb11, beliefPenaltyMPC_lbIdx11, beliefPenaltyMPC_rd11);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_62(beliefPenaltyMPC_Phi11, beliefPenaltyMPC_rd11, beliefPenaltyMPC_Lbyrd11);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_20_62_62(beliefPenaltyMPC_V10, beliefPenaltyMPC_Lbyrd10, beliefPenaltyMPC_W11, beliefPenaltyMPC_Lbyrd11, beliefPenaltyMPC_beta11);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_20_20(beliefPenaltyMPC_Lsd11, beliefPenaltyMPC_yy10, beliefPenaltyMPC_beta11, beliefPenaltyMPC_bmy11);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_20(beliefPenaltyMPC_Ld11, beliefPenaltyMPC_bmy11, beliefPenaltyMPC_yy11);
beliefPenaltyMPC_LA_VSUB6_INDEXED_62_22_62(beliefPenaltyMPC_ccrhsub12, beliefPenaltyMPC_sub12, beliefPenaltyMPC_ubIdx12, beliefPenaltyMPC_ccrhsl12, beliefPenaltyMPC_slb12, beliefPenaltyMPC_lbIdx12, beliefPenaltyMPC_rd12);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_62(beliefPenaltyMPC_Phi12, beliefPenaltyMPC_rd12, beliefPenaltyMPC_Lbyrd12);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_20_62_62(beliefPenaltyMPC_V11, beliefPenaltyMPC_Lbyrd11, beliefPenaltyMPC_W12, beliefPenaltyMPC_Lbyrd12, beliefPenaltyMPC_beta12);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_20_20(beliefPenaltyMPC_Lsd12, beliefPenaltyMPC_yy11, beliefPenaltyMPC_beta12, beliefPenaltyMPC_bmy12);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_20(beliefPenaltyMPC_Ld12, beliefPenaltyMPC_bmy12, beliefPenaltyMPC_yy12);
beliefPenaltyMPC_LA_VSUB6_INDEXED_62_22_62(beliefPenaltyMPC_ccrhsub13, beliefPenaltyMPC_sub13, beliefPenaltyMPC_ubIdx13, beliefPenaltyMPC_ccrhsl13, beliefPenaltyMPC_slb13, beliefPenaltyMPC_lbIdx13, beliefPenaltyMPC_rd13);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_62(beliefPenaltyMPC_Phi13, beliefPenaltyMPC_rd13, beliefPenaltyMPC_Lbyrd13);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_20_62_62(beliefPenaltyMPC_V12, beliefPenaltyMPC_Lbyrd12, beliefPenaltyMPC_W13, beliefPenaltyMPC_Lbyrd13, beliefPenaltyMPC_beta13);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_20_20(beliefPenaltyMPC_Lsd13, beliefPenaltyMPC_yy12, beliefPenaltyMPC_beta13, beliefPenaltyMPC_bmy13);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_20(beliefPenaltyMPC_Ld13, beliefPenaltyMPC_bmy13, beliefPenaltyMPC_yy13);
beliefPenaltyMPC_LA_VSUB6_INDEXED_20_20_20(beliefPenaltyMPC_ccrhsub14, beliefPenaltyMPC_sub14, beliefPenaltyMPC_ubIdx14, beliefPenaltyMPC_ccrhsl14, beliefPenaltyMPC_slb14, beliefPenaltyMPC_lbIdx14, beliefPenaltyMPC_rd14);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_20(beliefPenaltyMPC_Phi14, beliefPenaltyMPC_rd14, beliefPenaltyMPC_Lbyrd14);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_20_62_20(beliefPenaltyMPC_V13, beliefPenaltyMPC_Lbyrd13, beliefPenaltyMPC_W14, beliefPenaltyMPC_Lbyrd14, beliefPenaltyMPC_beta14);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_20_20(beliefPenaltyMPC_Lsd14, beliefPenaltyMPC_yy13, beliefPenaltyMPC_beta14, beliefPenaltyMPC_bmy14);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_20(beliefPenaltyMPC_Ld14, beliefPenaltyMPC_bmy14, beliefPenaltyMPC_yy14);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_20(beliefPenaltyMPC_Ld14, beliefPenaltyMPC_yy14, beliefPenaltyMPC_dvcc14);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_20_20(beliefPenaltyMPC_Lsd14, beliefPenaltyMPC_dvcc14, beliefPenaltyMPC_yy13, beliefPenaltyMPC_bmy13);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_20(beliefPenaltyMPC_Ld13, beliefPenaltyMPC_bmy13, beliefPenaltyMPC_dvcc13);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_20_20(beliefPenaltyMPC_Lsd13, beliefPenaltyMPC_dvcc13, beliefPenaltyMPC_yy12, beliefPenaltyMPC_bmy12);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_20(beliefPenaltyMPC_Ld12, beliefPenaltyMPC_bmy12, beliefPenaltyMPC_dvcc12);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_20_20(beliefPenaltyMPC_Lsd12, beliefPenaltyMPC_dvcc12, beliefPenaltyMPC_yy11, beliefPenaltyMPC_bmy11);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_20(beliefPenaltyMPC_Ld11, beliefPenaltyMPC_bmy11, beliefPenaltyMPC_dvcc11);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_20_20(beliefPenaltyMPC_Lsd11, beliefPenaltyMPC_dvcc11, beliefPenaltyMPC_yy10, beliefPenaltyMPC_bmy10);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_20(beliefPenaltyMPC_Ld10, beliefPenaltyMPC_bmy10, beliefPenaltyMPC_dvcc10);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_20_20(beliefPenaltyMPC_Lsd10, beliefPenaltyMPC_dvcc10, beliefPenaltyMPC_yy09, beliefPenaltyMPC_bmy09);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_20(beliefPenaltyMPC_Ld09, beliefPenaltyMPC_bmy09, beliefPenaltyMPC_dvcc09);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_20_20(beliefPenaltyMPC_Lsd09, beliefPenaltyMPC_dvcc09, beliefPenaltyMPC_yy08, beliefPenaltyMPC_bmy08);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_20(beliefPenaltyMPC_Ld08, beliefPenaltyMPC_bmy08, beliefPenaltyMPC_dvcc08);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_20_20(beliefPenaltyMPC_Lsd08, beliefPenaltyMPC_dvcc08, beliefPenaltyMPC_yy07, beliefPenaltyMPC_bmy07);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_20(beliefPenaltyMPC_Ld07, beliefPenaltyMPC_bmy07, beliefPenaltyMPC_dvcc07);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_20_20(beliefPenaltyMPC_Lsd07, beliefPenaltyMPC_dvcc07, beliefPenaltyMPC_yy06, beliefPenaltyMPC_bmy06);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_20(beliefPenaltyMPC_Ld06, beliefPenaltyMPC_bmy06, beliefPenaltyMPC_dvcc06);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_20_20(beliefPenaltyMPC_Lsd06, beliefPenaltyMPC_dvcc06, beliefPenaltyMPC_yy05, beliefPenaltyMPC_bmy05);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_20(beliefPenaltyMPC_Ld05, beliefPenaltyMPC_bmy05, beliefPenaltyMPC_dvcc05);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_20_20(beliefPenaltyMPC_Lsd05, beliefPenaltyMPC_dvcc05, beliefPenaltyMPC_yy04, beliefPenaltyMPC_bmy04);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_20(beliefPenaltyMPC_Ld04, beliefPenaltyMPC_bmy04, beliefPenaltyMPC_dvcc04);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_20_20(beliefPenaltyMPC_Lsd04, beliefPenaltyMPC_dvcc04, beliefPenaltyMPC_yy03, beliefPenaltyMPC_bmy03);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_20(beliefPenaltyMPC_Ld03, beliefPenaltyMPC_bmy03, beliefPenaltyMPC_dvcc03);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_20_20(beliefPenaltyMPC_Lsd03, beliefPenaltyMPC_dvcc03, beliefPenaltyMPC_yy02, beliefPenaltyMPC_bmy02);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_20(beliefPenaltyMPC_Ld02, beliefPenaltyMPC_bmy02, beliefPenaltyMPC_dvcc02);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_20_20(beliefPenaltyMPC_Lsd02, beliefPenaltyMPC_dvcc02, beliefPenaltyMPC_yy01, beliefPenaltyMPC_bmy01);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_20(beliefPenaltyMPC_Ld01, beliefPenaltyMPC_bmy01, beliefPenaltyMPC_dvcc01);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_20_20(beliefPenaltyMPC_Lsd01, beliefPenaltyMPC_dvcc01, beliefPenaltyMPC_yy00, beliefPenaltyMPC_bmy00);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_20(beliefPenaltyMPC_Ld00, beliefPenaltyMPC_bmy00, beliefPenaltyMPC_dvcc00);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_20_62_20(params->C1, beliefPenaltyMPC_dvcc01, beliefPenaltyMPC_D00, beliefPenaltyMPC_dvcc00, beliefPenaltyMPC_grad_eq00);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_20_62_20(params->C2, beliefPenaltyMPC_dvcc02, beliefPenaltyMPC_D01, beliefPenaltyMPC_dvcc01, beliefPenaltyMPC_grad_eq01);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_20_62_20(params->C3, beliefPenaltyMPC_dvcc03, beliefPenaltyMPC_D01, beliefPenaltyMPC_dvcc02, beliefPenaltyMPC_grad_eq02);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_20_62_20(params->C4, beliefPenaltyMPC_dvcc04, beliefPenaltyMPC_D01, beliefPenaltyMPC_dvcc03, beliefPenaltyMPC_grad_eq03);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_20_62_20(params->C5, beliefPenaltyMPC_dvcc05, beliefPenaltyMPC_D01, beliefPenaltyMPC_dvcc04, beliefPenaltyMPC_grad_eq04);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_20_62_20(params->C6, beliefPenaltyMPC_dvcc06, beliefPenaltyMPC_D01, beliefPenaltyMPC_dvcc05, beliefPenaltyMPC_grad_eq05);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_20_62_20(params->C7, beliefPenaltyMPC_dvcc07, beliefPenaltyMPC_D01, beliefPenaltyMPC_dvcc06, beliefPenaltyMPC_grad_eq06);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_20_62_20(params->C8, beliefPenaltyMPC_dvcc08, beliefPenaltyMPC_D01, beliefPenaltyMPC_dvcc07, beliefPenaltyMPC_grad_eq07);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_20_62_20(params->C9, beliefPenaltyMPC_dvcc09, beliefPenaltyMPC_D01, beliefPenaltyMPC_dvcc08, beliefPenaltyMPC_grad_eq08);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_20_62_20(params->C10, beliefPenaltyMPC_dvcc10, beliefPenaltyMPC_D01, beliefPenaltyMPC_dvcc09, beliefPenaltyMPC_grad_eq09);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_20_62_20(params->C11, beliefPenaltyMPC_dvcc11, beliefPenaltyMPC_D01, beliefPenaltyMPC_dvcc10, beliefPenaltyMPC_grad_eq10);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_20_62_20(params->C12, beliefPenaltyMPC_dvcc12, beliefPenaltyMPC_D01, beliefPenaltyMPC_dvcc11, beliefPenaltyMPC_grad_eq11);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_20_62_20(params->C13, beliefPenaltyMPC_dvcc13, beliefPenaltyMPC_D01, beliefPenaltyMPC_dvcc12, beliefPenaltyMPC_grad_eq12);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_20_62_20(params->C14, beliefPenaltyMPC_dvcc14, beliefPenaltyMPC_D01, beliefPenaltyMPC_dvcc13, beliefPenaltyMPC_grad_eq13);
beliefPenaltyMPC_LA_DIAGZERO_MTVM_20_20(beliefPenaltyMPC_D14, beliefPenaltyMPC_dvcc14, beliefPenaltyMPC_grad_eq14);
beliefPenaltyMPC_LA_VSUB_888(beliefPenaltyMPC_rd, beliefPenaltyMPC_grad_eq, beliefPenaltyMPC_rd);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_62(beliefPenaltyMPC_Phi00, beliefPenaltyMPC_rd00, beliefPenaltyMPC_dzcc00);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_62(beliefPenaltyMPC_Phi01, beliefPenaltyMPC_rd01, beliefPenaltyMPC_dzcc01);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_62(beliefPenaltyMPC_Phi02, beliefPenaltyMPC_rd02, beliefPenaltyMPC_dzcc02);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_62(beliefPenaltyMPC_Phi03, beliefPenaltyMPC_rd03, beliefPenaltyMPC_dzcc03);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_62(beliefPenaltyMPC_Phi04, beliefPenaltyMPC_rd04, beliefPenaltyMPC_dzcc04);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_62(beliefPenaltyMPC_Phi05, beliefPenaltyMPC_rd05, beliefPenaltyMPC_dzcc05);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_62(beliefPenaltyMPC_Phi06, beliefPenaltyMPC_rd06, beliefPenaltyMPC_dzcc06);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_62(beliefPenaltyMPC_Phi07, beliefPenaltyMPC_rd07, beliefPenaltyMPC_dzcc07);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_62(beliefPenaltyMPC_Phi08, beliefPenaltyMPC_rd08, beliefPenaltyMPC_dzcc08);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_62(beliefPenaltyMPC_Phi09, beliefPenaltyMPC_rd09, beliefPenaltyMPC_dzcc09);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_62(beliefPenaltyMPC_Phi10, beliefPenaltyMPC_rd10, beliefPenaltyMPC_dzcc10);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_62(beliefPenaltyMPC_Phi11, beliefPenaltyMPC_rd11, beliefPenaltyMPC_dzcc11);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_62(beliefPenaltyMPC_Phi12, beliefPenaltyMPC_rd12, beliefPenaltyMPC_dzcc12);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_62(beliefPenaltyMPC_Phi13, beliefPenaltyMPC_rd13, beliefPenaltyMPC_dzcc13);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_20(beliefPenaltyMPC_Phi14, beliefPenaltyMPC_rd14, beliefPenaltyMPC_dzcc14);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_62(beliefPenaltyMPC_ccrhsl00, beliefPenaltyMPC_slb00, beliefPenaltyMPC_llbbyslb00, beliefPenaltyMPC_dzcc00, beliefPenaltyMPC_lbIdx00, beliefPenaltyMPC_dllbcc00);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_22(beliefPenaltyMPC_ccrhsub00, beliefPenaltyMPC_sub00, beliefPenaltyMPC_lubbysub00, beliefPenaltyMPC_dzcc00, beliefPenaltyMPC_ubIdx00, beliefPenaltyMPC_dlubcc00);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_62(beliefPenaltyMPC_ccrhsl01, beliefPenaltyMPC_slb01, beliefPenaltyMPC_llbbyslb01, beliefPenaltyMPC_dzcc01, beliefPenaltyMPC_lbIdx01, beliefPenaltyMPC_dllbcc01);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_22(beliefPenaltyMPC_ccrhsub01, beliefPenaltyMPC_sub01, beliefPenaltyMPC_lubbysub01, beliefPenaltyMPC_dzcc01, beliefPenaltyMPC_ubIdx01, beliefPenaltyMPC_dlubcc01);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_62(beliefPenaltyMPC_ccrhsl02, beliefPenaltyMPC_slb02, beliefPenaltyMPC_llbbyslb02, beliefPenaltyMPC_dzcc02, beliefPenaltyMPC_lbIdx02, beliefPenaltyMPC_dllbcc02);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_22(beliefPenaltyMPC_ccrhsub02, beliefPenaltyMPC_sub02, beliefPenaltyMPC_lubbysub02, beliefPenaltyMPC_dzcc02, beliefPenaltyMPC_ubIdx02, beliefPenaltyMPC_dlubcc02);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_62(beliefPenaltyMPC_ccrhsl03, beliefPenaltyMPC_slb03, beliefPenaltyMPC_llbbyslb03, beliefPenaltyMPC_dzcc03, beliefPenaltyMPC_lbIdx03, beliefPenaltyMPC_dllbcc03);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_22(beliefPenaltyMPC_ccrhsub03, beliefPenaltyMPC_sub03, beliefPenaltyMPC_lubbysub03, beliefPenaltyMPC_dzcc03, beliefPenaltyMPC_ubIdx03, beliefPenaltyMPC_dlubcc03);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_62(beliefPenaltyMPC_ccrhsl04, beliefPenaltyMPC_slb04, beliefPenaltyMPC_llbbyslb04, beliefPenaltyMPC_dzcc04, beliefPenaltyMPC_lbIdx04, beliefPenaltyMPC_dllbcc04);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_22(beliefPenaltyMPC_ccrhsub04, beliefPenaltyMPC_sub04, beliefPenaltyMPC_lubbysub04, beliefPenaltyMPC_dzcc04, beliefPenaltyMPC_ubIdx04, beliefPenaltyMPC_dlubcc04);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_62(beliefPenaltyMPC_ccrhsl05, beliefPenaltyMPC_slb05, beliefPenaltyMPC_llbbyslb05, beliefPenaltyMPC_dzcc05, beliefPenaltyMPC_lbIdx05, beliefPenaltyMPC_dllbcc05);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_22(beliefPenaltyMPC_ccrhsub05, beliefPenaltyMPC_sub05, beliefPenaltyMPC_lubbysub05, beliefPenaltyMPC_dzcc05, beliefPenaltyMPC_ubIdx05, beliefPenaltyMPC_dlubcc05);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_62(beliefPenaltyMPC_ccrhsl06, beliefPenaltyMPC_slb06, beliefPenaltyMPC_llbbyslb06, beliefPenaltyMPC_dzcc06, beliefPenaltyMPC_lbIdx06, beliefPenaltyMPC_dllbcc06);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_22(beliefPenaltyMPC_ccrhsub06, beliefPenaltyMPC_sub06, beliefPenaltyMPC_lubbysub06, beliefPenaltyMPC_dzcc06, beliefPenaltyMPC_ubIdx06, beliefPenaltyMPC_dlubcc06);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_62(beliefPenaltyMPC_ccrhsl07, beliefPenaltyMPC_slb07, beliefPenaltyMPC_llbbyslb07, beliefPenaltyMPC_dzcc07, beliefPenaltyMPC_lbIdx07, beliefPenaltyMPC_dllbcc07);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_22(beliefPenaltyMPC_ccrhsub07, beliefPenaltyMPC_sub07, beliefPenaltyMPC_lubbysub07, beliefPenaltyMPC_dzcc07, beliefPenaltyMPC_ubIdx07, beliefPenaltyMPC_dlubcc07);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_62(beliefPenaltyMPC_ccrhsl08, beliefPenaltyMPC_slb08, beliefPenaltyMPC_llbbyslb08, beliefPenaltyMPC_dzcc08, beliefPenaltyMPC_lbIdx08, beliefPenaltyMPC_dllbcc08);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_22(beliefPenaltyMPC_ccrhsub08, beliefPenaltyMPC_sub08, beliefPenaltyMPC_lubbysub08, beliefPenaltyMPC_dzcc08, beliefPenaltyMPC_ubIdx08, beliefPenaltyMPC_dlubcc08);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_62(beliefPenaltyMPC_ccrhsl09, beliefPenaltyMPC_slb09, beliefPenaltyMPC_llbbyslb09, beliefPenaltyMPC_dzcc09, beliefPenaltyMPC_lbIdx09, beliefPenaltyMPC_dllbcc09);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_22(beliefPenaltyMPC_ccrhsub09, beliefPenaltyMPC_sub09, beliefPenaltyMPC_lubbysub09, beliefPenaltyMPC_dzcc09, beliefPenaltyMPC_ubIdx09, beliefPenaltyMPC_dlubcc09);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_62(beliefPenaltyMPC_ccrhsl10, beliefPenaltyMPC_slb10, beliefPenaltyMPC_llbbyslb10, beliefPenaltyMPC_dzcc10, beliefPenaltyMPC_lbIdx10, beliefPenaltyMPC_dllbcc10);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_22(beliefPenaltyMPC_ccrhsub10, beliefPenaltyMPC_sub10, beliefPenaltyMPC_lubbysub10, beliefPenaltyMPC_dzcc10, beliefPenaltyMPC_ubIdx10, beliefPenaltyMPC_dlubcc10);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_62(beliefPenaltyMPC_ccrhsl11, beliefPenaltyMPC_slb11, beliefPenaltyMPC_llbbyslb11, beliefPenaltyMPC_dzcc11, beliefPenaltyMPC_lbIdx11, beliefPenaltyMPC_dllbcc11);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_22(beliefPenaltyMPC_ccrhsub11, beliefPenaltyMPC_sub11, beliefPenaltyMPC_lubbysub11, beliefPenaltyMPC_dzcc11, beliefPenaltyMPC_ubIdx11, beliefPenaltyMPC_dlubcc11);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_62(beliefPenaltyMPC_ccrhsl12, beliefPenaltyMPC_slb12, beliefPenaltyMPC_llbbyslb12, beliefPenaltyMPC_dzcc12, beliefPenaltyMPC_lbIdx12, beliefPenaltyMPC_dllbcc12);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_22(beliefPenaltyMPC_ccrhsub12, beliefPenaltyMPC_sub12, beliefPenaltyMPC_lubbysub12, beliefPenaltyMPC_dzcc12, beliefPenaltyMPC_ubIdx12, beliefPenaltyMPC_dlubcc12);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_62(beliefPenaltyMPC_ccrhsl13, beliefPenaltyMPC_slb13, beliefPenaltyMPC_llbbyslb13, beliefPenaltyMPC_dzcc13, beliefPenaltyMPC_lbIdx13, beliefPenaltyMPC_dllbcc13);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_22(beliefPenaltyMPC_ccrhsub13, beliefPenaltyMPC_sub13, beliefPenaltyMPC_lubbysub13, beliefPenaltyMPC_dzcc13, beliefPenaltyMPC_ubIdx13, beliefPenaltyMPC_dlubcc13);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_20(beliefPenaltyMPC_ccrhsl14, beliefPenaltyMPC_slb14, beliefPenaltyMPC_llbbyslb14, beliefPenaltyMPC_dzcc14, beliefPenaltyMPC_lbIdx14, beliefPenaltyMPC_dllbcc14);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_20(beliefPenaltyMPC_ccrhsub14, beliefPenaltyMPC_sub14, beliefPenaltyMPC_lubbysub14, beliefPenaltyMPC_dzcc14, beliefPenaltyMPC_ubIdx14, beliefPenaltyMPC_dlubcc14);
beliefPenaltyMPC_LA_VSUB7_1216(beliefPenaltyMPC_l, beliefPenaltyMPC_ccrhs, beliefPenaltyMPC_s, beliefPenaltyMPC_dl_cc, beliefPenaltyMPC_ds_cc);
beliefPenaltyMPC_LA_VADD_888(beliefPenaltyMPC_dz_cc, beliefPenaltyMPC_dz_aff);
beliefPenaltyMPC_LA_VADD_300(beliefPenaltyMPC_dv_cc, beliefPenaltyMPC_dv_aff);
beliefPenaltyMPC_LA_VADD_1216(beliefPenaltyMPC_dl_cc, beliefPenaltyMPC_dl_aff);
beliefPenaltyMPC_LA_VADD_1216(beliefPenaltyMPC_ds_cc, beliefPenaltyMPC_ds_aff);
info->lsit_cc = beliefPenaltyMPC_LINESEARCH_BACKTRACKING_COMBINED(beliefPenaltyMPC_z, beliefPenaltyMPC_v, beliefPenaltyMPC_l, beliefPenaltyMPC_s, beliefPenaltyMPC_dz_cc, beliefPenaltyMPC_dv_cc, beliefPenaltyMPC_dl_cc, beliefPenaltyMPC_ds_cc, &info->step_cc, &info->mu);
if( info->lsit_cc == beliefPenaltyMPC_NOPROGRESS ){
PRINTTEXT("Line search could not proceed at iteration %d, exiting.\n",info->it+1);
exitcode = beliefPenaltyMPC_NOPROGRESS; break;
}
info->it++;
}
output->z1[0] = beliefPenaltyMPC_z00[0];
output->z1[1] = beliefPenaltyMPC_z00[1];
output->z1[2] = beliefPenaltyMPC_z00[2];
output->z1[3] = beliefPenaltyMPC_z00[3];
output->z1[4] = beliefPenaltyMPC_z00[4];
output->z1[5] = beliefPenaltyMPC_z00[5];
output->z1[6] = beliefPenaltyMPC_z00[6];
output->z1[7] = beliefPenaltyMPC_z00[7];
output->z1[8] = beliefPenaltyMPC_z00[8];
output->z1[9] = beliefPenaltyMPC_z00[9];
output->z1[10] = beliefPenaltyMPC_z00[10];
output->z1[11] = beliefPenaltyMPC_z00[11];
output->z1[12] = beliefPenaltyMPC_z00[12];
output->z1[13] = beliefPenaltyMPC_z00[13];
output->z1[14] = beliefPenaltyMPC_z00[14];
output->z1[15] = beliefPenaltyMPC_z00[15];
output->z1[16] = beliefPenaltyMPC_z00[16];
output->z1[17] = beliefPenaltyMPC_z00[17];
output->z1[18] = beliefPenaltyMPC_z00[18];
output->z1[19] = beliefPenaltyMPC_z00[19];
output->z1[20] = beliefPenaltyMPC_z00[20];
output->z1[21] = beliefPenaltyMPC_z00[21];
output->z2[0] = beliefPenaltyMPC_z01[0];
output->z2[1] = beliefPenaltyMPC_z01[1];
output->z2[2] = beliefPenaltyMPC_z01[2];
output->z2[3] = beliefPenaltyMPC_z01[3];
output->z2[4] = beliefPenaltyMPC_z01[4];
output->z2[5] = beliefPenaltyMPC_z01[5];
output->z2[6] = beliefPenaltyMPC_z01[6];
output->z2[7] = beliefPenaltyMPC_z01[7];
output->z2[8] = beliefPenaltyMPC_z01[8];
output->z2[9] = beliefPenaltyMPC_z01[9];
output->z2[10] = beliefPenaltyMPC_z01[10];
output->z2[11] = beliefPenaltyMPC_z01[11];
output->z2[12] = beliefPenaltyMPC_z01[12];
output->z2[13] = beliefPenaltyMPC_z01[13];
output->z2[14] = beliefPenaltyMPC_z01[14];
output->z2[15] = beliefPenaltyMPC_z01[15];
output->z2[16] = beliefPenaltyMPC_z01[16];
output->z2[17] = beliefPenaltyMPC_z01[17];
output->z2[18] = beliefPenaltyMPC_z01[18];
output->z2[19] = beliefPenaltyMPC_z01[19];
output->z2[20] = beliefPenaltyMPC_z01[20];
output->z2[21] = beliefPenaltyMPC_z01[21];
output->z3[0] = beliefPenaltyMPC_z02[0];
output->z3[1] = beliefPenaltyMPC_z02[1];
output->z3[2] = beliefPenaltyMPC_z02[2];
output->z3[3] = beliefPenaltyMPC_z02[3];
output->z3[4] = beliefPenaltyMPC_z02[4];
output->z3[5] = beliefPenaltyMPC_z02[5];
output->z3[6] = beliefPenaltyMPC_z02[6];
output->z3[7] = beliefPenaltyMPC_z02[7];
output->z3[8] = beliefPenaltyMPC_z02[8];
output->z3[9] = beliefPenaltyMPC_z02[9];
output->z3[10] = beliefPenaltyMPC_z02[10];
output->z3[11] = beliefPenaltyMPC_z02[11];
output->z3[12] = beliefPenaltyMPC_z02[12];
output->z3[13] = beliefPenaltyMPC_z02[13];
output->z3[14] = beliefPenaltyMPC_z02[14];
output->z3[15] = beliefPenaltyMPC_z02[15];
output->z3[16] = beliefPenaltyMPC_z02[16];
output->z3[17] = beliefPenaltyMPC_z02[17];
output->z3[18] = beliefPenaltyMPC_z02[18];
output->z3[19] = beliefPenaltyMPC_z02[19];
output->z3[20] = beliefPenaltyMPC_z02[20];
output->z3[21] = beliefPenaltyMPC_z02[21];
output->z4[0] = beliefPenaltyMPC_z03[0];
output->z4[1] = beliefPenaltyMPC_z03[1];
output->z4[2] = beliefPenaltyMPC_z03[2];
output->z4[3] = beliefPenaltyMPC_z03[3];
output->z4[4] = beliefPenaltyMPC_z03[4];
output->z4[5] = beliefPenaltyMPC_z03[5];
output->z4[6] = beliefPenaltyMPC_z03[6];
output->z4[7] = beliefPenaltyMPC_z03[7];
output->z4[8] = beliefPenaltyMPC_z03[8];
output->z4[9] = beliefPenaltyMPC_z03[9];
output->z4[10] = beliefPenaltyMPC_z03[10];
output->z4[11] = beliefPenaltyMPC_z03[11];
output->z4[12] = beliefPenaltyMPC_z03[12];
output->z4[13] = beliefPenaltyMPC_z03[13];
output->z4[14] = beliefPenaltyMPC_z03[14];
output->z4[15] = beliefPenaltyMPC_z03[15];
output->z4[16] = beliefPenaltyMPC_z03[16];
output->z4[17] = beliefPenaltyMPC_z03[17];
output->z4[18] = beliefPenaltyMPC_z03[18];
output->z4[19] = beliefPenaltyMPC_z03[19];
output->z4[20] = beliefPenaltyMPC_z03[20];
output->z4[21] = beliefPenaltyMPC_z03[21];
output->z5[0] = beliefPenaltyMPC_z04[0];
output->z5[1] = beliefPenaltyMPC_z04[1];
output->z5[2] = beliefPenaltyMPC_z04[2];
output->z5[3] = beliefPenaltyMPC_z04[3];
output->z5[4] = beliefPenaltyMPC_z04[4];
output->z5[5] = beliefPenaltyMPC_z04[5];
output->z5[6] = beliefPenaltyMPC_z04[6];
output->z5[7] = beliefPenaltyMPC_z04[7];
output->z5[8] = beliefPenaltyMPC_z04[8];
output->z5[9] = beliefPenaltyMPC_z04[9];
output->z5[10] = beliefPenaltyMPC_z04[10];
output->z5[11] = beliefPenaltyMPC_z04[11];
output->z5[12] = beliefPenaltyMPC_z04[12];
output->z5[13] = beliefPenaltyMPC_z04[13];
output->z5[14] = beliefPenaltyMPC_z04[14];
output->z5[15] = beliefPenaltyMPC_z04[15];
output->z5[16] = beliefPenaltyMPC_z04[16];
output->z5[17] = beliefPenaltyMPC_z04[17];
output->z5[18] = beliefPenaltyMPC_z04[18];
output->z5[19] = beliefPenaltyMPC_z04[19];
output->z5[20] = beliefPenaltyMPC_z04[20];
output->z5[21] = beliefPenaltyMPC_z04[21];
output->z6[0] = beliefPenaltyMPC_z05[0];
output->z6[1] = beliefPenaltyMPC_z05[1];
output->z6[2] = beliefPenaltyMPC_z05[2];
output->z6[3] = beliefPenaltyMPC_z05[3];
output->z6[4] = beliefPenaltyMPC_z05[4];
output->z6[5] = beliefPenaltyMPC_z05[5];
output->z6[6] = beliefPenaltyMPC_z05[6];
output->z6[7] = beliefPenaltyMPC_z05[7];
output->z6[8] = beliefPenaltyMPC_z05[8];
output->z6[9] = beliefPenaltyMPC_z05[9];
output->z6[10] = beliefPenaltyMPC_z05[10];
output->z6[11] = beliefPenaltyMPC_z05[11];
output->z6[12] = beliefPenaltyMPC_z05[12];
output->z6[13] = beliefPenaltyMPC_z05[13];
output->z6[14] = beliefPenaltyMPC_z05[14];
output->z6[15] = beliefPenaltyMPC_z05[15];
output->z6[16] = beliefPenaltyMPC_z05[16];
output->z6[17] = beliefPenaltyMPC_z05[17];
output->z6[18] = beliefPenaltyMPC_z05[18];
output->z6[19] = beliefPenaltyMPC_z05[19];
output->z6[20] = beliefPenaltyMPC_z05[20];
output->z6[21] = beliefPenaltyMPC_z05[21];
output->z7[0] = beliefPenaltyMPC_z06[0];
output->z7[1] = beliefPenaltyMPC_z06[1];
output->z7[2] = beliefPenaltyMPC_z06[2];
output->z7[3] = beliefPenaltyMPC_z06[3];
output->z7[4] = beliefPenaltyMPC_z06[4];
output->z7[5] = beliefPenaltyMPC_z06[5];
output->z7[6] = beliefPenaltyMPC_z06[6];
output->z7[7] = beliefPenaltyMPC_z06[7];
output->z7[8] = beliefPenaltyMPC_z06[8];
output->z7[9] = beliefPenaltyMPC_z06[9];
output->z7[10] = beliefPenaltyMPC_z06[10];
output->z7[11] = beliefPenaltyMPC_z06[11];
output->z7[12] = beliefPenaltyMPC_z06[12];
output->z7[13] = beliefPenaltyMPC_z06[13];
output->z7[14] = beliefPenaltyMPC_z06[14];
output->z7[15] = beliefPenaltyMPC_z06[15];
output->z7[16] = beliefPenaltyMPC_z06[16];
output->z7[17] = beliefPenaltyMPC_z06[17];
output->z7[18] = beliefPenaltyMPC_z06[18];
output->z7[19] = beliefPenaltyMPC_z06[19];
output->z7[20] = beliefPenaltyMPC_z06[20];
output->z7[21] = beliefPenaltyMPC_z06[21];
output->z8[0] = beliefPenaltyMPC_z07[0];
output->z8[1] = beliefPenaltyMPC_z07[1];
output->z8[2] = beliefPenaltyMPC_z07[2];
output->z8[3] = beliefPenaltyMPC_z07[3];
output->z8[4] = beliefPenaltyMPC_z07[4];
output->z8[5] = beliefPenaltyMPC_z07[5];
output->z8[6] = beliefPenaltyMPC_z07[6];
output->z8[7] = beliefPenaltyMPC_z07[7];
output->z8[8] = beliefPenaltyMPC_z07[8];
output->z8[9] = beliefPenaltyMPC_z07[9];
output->z8[10] = beliefPenaltyMPC_z07[10];
output->z8[11] = beliefPenaltyMPC_z07[11];
output->z8[12] = beliefPenaltyMPC_z07[12];
output->z8[13] = beliefPenaltyMPC_z07[13];
output->z8[14] = beliefPenaltyMPC_z07[14];
output->z8[15] = beliefPenaltyMPC_z07[15];
output->z8[16] = beliefPenaltyMPC_z07[16];
output->z8[17] = beliefPenaltyMPC_z07[17];
output->z8[18] = beliefPenaltyMPC_z07[18];
output->z8[19] = beliefPenaltyMPC_z07[19];
output->z8[20] = beliefPenaltyMPC_z07[20];
output->z8[21] = beliefPenaltyMPC_z07[21];
output->z9[0] = beliefPenaltyMPC_z08[0];
output->z9[1] = beliefPenaltyMPC_z08[1];
output->z9[2] = beliefPenaltyMPC_z08[2];
output->z9[3] = beliefPenaltyMPC_z08[3];
output->z9[4] = beliefPenaltyMPC_z08[4];
output->z9[5] = beliefPenaltyMPC_z08[5];
output->z9[6] = beliefPenaltyMPC_z08[6];
output->z9[7] = beliefPenaltyMPC_z08[7];
output->z9[8] = beliefPenaltyMPC_z08[8];
output->z9[9] = beliefPenaltyMPC_z08[9];
output->z9[10] = beliefPenaltyMPC_z08[10];
output->z9[11] = beliefPenaltyMPC_z08[11];
output->z9[12] = beliefPenaltyMPC_z08[12];
output->z9[13] = beliefPenaltyMPC_z08[13];
output->z9[14] = beliefPenaltyMPC_z08[14];
output->z9[15] = beliefPenaltyMPC_z08[15];
output->z9[16] = beliefPenaltyMPC_z08[16];
output->z9[17] = beliefPenaltyMPC_z08[17];
output->z9[18] = beliefPenaltyMPC_z08[18];
output->z9[19] = beliefPenaltyMPC_z08[19];
output->z9[20] = beliefPenaltyMPC_z08[20];
output->z9[21] = beliefPenaltyMPC_z08[21];
output->z10[0] = beliefPenaltyMPC_z09[0];
output->z10[1] = beliefPenaltyMPC_z09[1];
output->z10[2] = beliefPenaltyMPC_z09[2];
output->z10[3] = beliefPenaltyMPC_z09[3];
output->z10[4] = beliefPenaltyMPC_z09[4];
output->z10[5] = beliefPenaltyMPC_z09[5];
output->z10[6] = beliefPenaltyMPC_z09[6];
output->z10[7] = beliefPenaltyMPC_z09[7];
output->z10[8] = beliefPenaltyMPC_z09[8];
output->z10[9] = beliefPenaltyMPC_z09[9];
output->z10[10] = beliefPenaltyMPC_z09[10];
output->z10[11] = beliefPenaltyMPC_z09[11];
output->z10[12] = beliefPenaltyMPC_z09[12];
output->z10[13] = beliefPenaltyMPC_z09[13];
output->z10[14] = beliefPenaltyMPC_z09[14];
output->z10[15] = beliefPenaltyMPC_z09[15];
output->z10[16] = beliefPenaltyMPC_z09[16];
output->z10[17] = beliefPenaltyMPC_z09[17];
output->z10[18] = beliefPenaltyMPC_z09[18];
output->z10[19] = beliefPenaltyMPC_z09[19];
output->z10[20] = beliefPenaltyMPC_z09[20];
output->z10[21] = beliefPenaltyMPC_z09[21];
output->z11[0] = beliefPenaltyMPC_z10[0];
output->z11[1] = beliefPenaltyMPC_z10[1];
output->z11[2] = beliefPenaltyMPC_z10[2];
output->z11[3] = beliefPenaltyMPC_z10[3];
output->z11[4] = beliefPenaltyMPC_z10[4];
output->z11[5] = beliefPenaltyMPC_z10[5];
output->z11[6] = beliefPenaltyMPC_z10[6];
output->z11[7] = beliefPenaltyMPC_z10[7];
output->z11[8] = beliefPenaltyMPC_z10[8];
output->z11[9] = beliefPenaltyMPC_z10[9];
output->z11[10] = beliefPenaltyMPC_z10[10];
output->z11[11] = beliefPenaltyMPC_z10[11];
output->z11[12] = beliefPenaltyMPC_z10[12];
output->z11[13] = beliefPenaltyMPC_z10[13];
output->z11[14] = beliefPenaltyMPC_z10[14];
output->z11[15] = beliefPenaltyMPC_z10[15];
output->z11[16] = beliefPenaltyMPC_z10[16];
output->z11[17] = beliefPenaltyMPC_z10[17];
output->z11[18] = beliefPenaltyMPC_z10[18];
output->z11[19] = beliefPenaltyMPC_z10[19];
output->z11[20] = beliefPenaltyMPC_z10[20];
output->z11[21] = beliefPenaltyMPC_z10[21];
output->z12[0] = beliefPenaltyMPC_z11[0];
output->z12[1] = beliefPenaltyMPC_z11[1];
output->z12[2] = beliefPenaltyMPC_z11[2];
output->z12[3] = beliefPenaltyMPC_z11[3];
output->z12[4] = beliefPenaltyMPC_z11[4];
output->z12[5] = beliefPenaltyMPC_z11[5];
output->z12[6] = beliefPenaltyMPC_z11[6];
output->z12[7] = beliefPenaltyMPC_z11[7];
output->z12[8] = beliefPenaltyMPC_z11[8];
output->z12[9] = beliefPenaltyMPC_z11[9];
output->z12[10] = beliefPenaltyMPC_z11[10];
output->z12[11] = beliefPenaltyMPC_z11[11];
output->z12[12] = beliefPenaltyMPC_z11[12];
output->z12[13] = beliefPenaltyMPC_z11[13];
output->z12[14] = beliefPenaltyMPC_z11[14];
output->z12[15] = beliefPenaltyMPC_z11[15];
output->z12[16] = beliefPenaltyMPC_z11[16];
output->z12[17] = beliefPenaltyMPC_z11[17];
output->z12[18] = beliefPenaltyMPC_z11[18];
output->z12[19] = beliefPenaltyMPC_z11[19];
output->z12[20] = beliefPenaltyMPC_z11[20];
output->z12[21] = beliefPenaltyMPC_z11[21];
output->z13[0] = beliefPenaltyMPC_z12[0];
output->z13[1] = beliefPenaltyMPC_z12[1];
output->z13[2] = beliefPenaltyMPC_z12[2];
output->z13[3] = beliefPenaltyMPC_z12[3];
output->z13[4] = beliefPenaltyMPC_z12[4];
output->z13[5] = beliefPenaltyMPC_z12[5];
output->z13[6] = beliefPenaltyMPC_z12[6];
output->z13[7] = beliefPenaltyMPC_z12[7];
output->z13[8] = beliefPenaltyMPC_z12[8];
output->z13[9] = beliefPenaltyMPC_z12[9];
output->z13[10] = beliefPenaltyMPC_z12[10];
output->z13[11] = beliefPenaltyMPC_z12[11];
output->z13[12] = beliefPenaltyMPC_z12[12];
output->z13[13] = beliefPenaltyMPC_z12[13];
output->z13[14] = beliefPenaltyMPC_z12[14];
output->z13[15] = beliefPenaltyMPC_z12[15];
output->z13[16] = beliefPenaltyMPC_z12[16];
output->z13[17] = beliefPenaltyMPC_z12[17];
output->z13[18] = beliefPenaltyMPC_z12[18];
output->z13[19] = beliefPenaltyMPC_z12[19];
output->z13[20] = beliefPenaltyMPC_z12[20];
output->z13[21] = beliefPenaltyMPC_z12[21];
output->z14[0] = beliefPenaltyMPC_z13[0];
output->z14[1] = beliefPenaltyMPC_z13[1];
output->z14[2] = beliefPenaltyMPC_z13[2];
output->z14[3] = beliefPenaltyMPC_z13[3];
output->z14[4] = beliefPenaltyMPC_z13[4];
output->z14[5] = beliefPenaltyMPC_z13[5];
output->z14[6] = beliefPenaltyMPC_z13[6];
output->z14[7] = beliefPenaltyMPC_z13[7];
output->z14[8] = beliefPenaltyMPC_z13[8];
output->z14[9] = beliefPenaltyMPC_z13[9];
output->z14[10] = beliefPenaltyMPC_z13[10];
output->z14[11] = beliefPenaltyMPC_z13[11];
output->z14[12] = beliefPenaltyMPC_z13[12];
output->z14[13] = beliefPenaltyMPC_z13[13];
output->z14[14] = beliefPenaltyMPC_z13[14];
output->z14[15] = beliefPenaltyMPC_z13[15];
output->z14[16] = beliefPenaltyMPC_z13[16];
output->z14[17] = beliefPenaltyMPC_z13[17];
output->z14[18] = beliefPenaltyMPC_z13[18];
output->z14[19] = beliefPenaltyMPC_z13[19];
output->z14[20] = beliefPenaltyMPC_z13[20];
output->z14[21] = beliefPenaltyMPC_z13[21];
output->z15[0] = beliefPenaltyMPC_z14[0];
output->z15[1] = beliefPenaltyMPC_z14[1];
output->z15[2] = beliefPenaltyMPC_z14[2];
output->z15[3] = beliefPenaltyMPC_z14[3];
output->z15[4] = beliefPenaltyMPC_z14[4];
output->z15[5] = beliefPenaltyMPC_z14[5];
output->z15[6] = beliefPenaltyMPC_z14[6];
output->z15[7] = beliefPenaltyMPC_z14[7];
output->z15[8] = beliefPenaltyMPC_z14[8];
output->z15[9] = beliefPenaltyMPC_z14[9];
output->z15[10] = beliefPenaltyMPC_z14[10];
output->z15[11] = beliefPenaltyMPC_z14[11];
output->z15[12] = beliefPenaltyMPC_z14[12];
output->z15[13] = beliefPenaltyMPC_z14[13];
output->z15[14] = beliefPenaltyMPC_z14[14];
output->z15[15] = beliefPenaltyMPC_z14[15];
output->z15[16] = beliefPenaltyMPC_z14[16];
output->z15[17] = beliefPenaltyMPC_z14[17];
output->z15[18] = beliefPenaltyMPC_z14[18];
output->z15[19] = beliefPenaltyMPC_z14[19];

#if beliefPenaltyMPC_SET_TIMING == 1
info->solvetime = beliefPenaltyMPC_toc(&solvertimer);
#if beliefPenaltyMPC_SET_PRINTLEVEL > 0 && beliefPenaltyMPC_SET_TIMING == 1
if( info->it > 1 ){
	PRINTTEXT("Solve time: %5.3f ms (%d iterations)\n\n", info->solvetime*1000, info->it);
} else {
	PRINTTEXT("Solve time: %5.3f ms (%d iteration)\n\n", info->solvetime*1000, info->it);
}
#endif
#else
info->solvetime = -1;
#endif
return exitcode;
}
