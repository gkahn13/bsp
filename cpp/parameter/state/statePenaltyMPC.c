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

#include "statePenaltyMPC.h"

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
 * Initializes a vector of length 372 with a value.
 */
void statePenaltyMPC_LA_INITIALIZEVECTOR_372(statePenaltyMPC_FLOAT* vec, statePenaltyMPC_FLOAT value)
{
	int i;
	for( i=0; i<372; i++ )
	{
		vec[i] = value;
	}
}


/*
 * Initializes a vector of length 120 with a value.
 */
void statePenaltyMPC_LA_INITIALIZEVECTOR_120(statePenaltyMPC_FLOAT* vec, statePenaltyMPC_FLOAT value)
{
	int i;
	for( i=0; i<120; i++ )
	{
		vec[i] = value;
	}
}


/*
 * Initializes a vector of length 520 with a value.
 */
void statePenaltyMPC_LA_INITIALIZEVECTOR_520(statePenaltyMPC_FLOAT* vec, statePenaltyMPC_FLOAT value)
{
	int i;
	for( i=0; i<520; i++ )
	{
		vec[i] = value;
	}
}


/* 
 * Calculates a dot product and adds it to a variable: z += x'*y; 
 * This function is for vectors of length 520.
 */
void statePenaltyMPC_LA_DOTACC_520(statePenaltyMPC_FLOAT *x, statePenaltyMPC_FLOAT *y, statePenaltyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<520; i++ ){
		*z += x[i]*y[i];
	}
}


/*
 * Calculates the gradient and the value for a quadratic function 0.5*z'*H*z + f'*z
 *
 * INPUTS:     H  - Symmetric Hessian, diag matrix of size [26 x 26]
 *             f  - column vector of size 26
 *             z  - column vector of size 26
 *
 * OUTPUTS: grad  - gradient at z (= H*z + f), column vector of size 26
 *          value <-- value + 0.5*z'*H*z + f'*z (value will be modified)
 */
void statePenaltyMPC_LA_DIAG_QUADFCN_26(statePenaltyMPC_FLOAT* H, statePenaltyMPC_FLOAT* f, statePenaltyMPC_FLOAT* z, statePenaltyMPC_FLOAT* grad, statePenaltyMPC_FLOAT* value)
{
	int i;
	statePenaltyMPC_FLOAT hz;	
	for( i=0; i<26; i++){
		hz = H[i]*z[i];
		grad[i] = hz + f[i];
		*value += 0.5*hz*z[i] + f[i]*z[i];
	}
}


/*
 * Calculates the gradient and the value for a quadratic function 0.5*z'*H*z + f'*z
 *
 * INPUTS:     H  - Symmetric Hessian, diag matrix of size [8 x 8]
 *             f  - column vector of size 8
 *             z  - column vector of size 8
 *
 * OUTPUTS: grad  - gradient at z (= H*z + f), column vector of size 8
 *          value <-- value + 0.5*z'*H*z + f'*z (value will be modified)
 */
void statePenaltyMPC_LA_DIAG_QUADFCN_8(statePenaltyMPC_FLOAT* H, statePenaltyMPC_FLOAT* f, statePenaltyMPC_FLOAT* z, statePenaltyMPC_FLOAT* grad, statePenaltyMPC_FLOAT* value)
{
	int i;
	statePenaltyMPC_FLOAT hz;	
	for( i=0; i<8; i++){
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
void statePenaltyMPC_LA_DIAGZERO_MVMSUB6_8(statePenaltyMPC_FLOAT *B, statePenaltyMPC_FLOAT *u, statePenaltyMPC_FLOAT *b, statePenaltyMPC_FLOAT *l, statePenaltyMPC_FLOAT *r, statePenaltyMPC_FLOAT *z, statePenaltyMPC_FLOAT *y)
{
	int i;
	statePenaltyMPC_FLOAT Bu[8];
	statePenaltyMPC_FLOAT norm = *y;
	statePenaltyMPC_FLOAT lr = 0;

	/* do A*x + B*u first */
	for( i=0; i<8; i++ ){
		Bu[i] = B[i]*u[i];
	}	

	for( i=0; i<8; i++ ){
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
void statePenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_8_26_26(statePenaltyMPC_FLOAT *A, statePenaltyMPC_FLOAT *x, statePenaltyMPC_FLOAT *B, statePenaltyMPC_FLOAT *u, statePenaltyMPC_FLOAT *b, statePenaltyMPC_FLOAT *l, statePenaltyMPC_FLOAT *r, statePenaltyMPC_FLOAT *z, statePenaltyMPC_FLOAT *y)
{
	int i;
	int j;
	int k = 0;
	statePenaltyMPC_FLOAT AxBu[8];
	statePenaltyMPC_FLOAT norm = *y;
	statePenaltyMPC_FLOAT lr = 0;

	/* do A*x + B*u first */
	for( i=0; i<8; i++ ){
		AxBu[i] = A[k++]*x[0] + B[i]*u[i];
	}	

	for( j=1; j<26; j++ ){		
		for( i=0; i<8; i++ ){
			AxBu[i] += A[k++]*x[j];
		}
	}

	for( i=0; i<8; i++ ){
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
void statePenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_8_26_8(statePenaltyMPC_FLOAT *A, statePenaltyMPC_FLOAT *x, statePenaltyMPC_FLOAT *B, statePenaltyMPC_FLOAT *u, statePenaltyMPC_FLOAT *b, statePenaltyMPC_FLOAT *l, statePenaltyMPC_FLOAT *r, statePenaltyMPC_FLOAT *z, statePenaltyMPC_FLOAT *y)
{
	int i;
	int j;
	int k = 0;
	statePenaltyMPC_FLOAT AxBu[8];
	statePenaltyMPC_FLOAT norm = *y;
	statePenaltyMPC_FLOAT lr = 0;

	/* do A*x + B*u first */
	for( i=0; i<8; i++ ){
		AxBu[i] = A[k++]*x[0] + B[i]*u[i];
	}	

	for( j=1; j<26; j++ ){		
		for( i=0; i<8; i++ ){
			AxBu[i] += A[k++]*x[j];
		}
	}

	for( i=0; i<8; i++ ){
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
 * where A is of size [8 x 26] and stored in column major format.
 * and B is of size [8 x 26] and stored in diagzero format
 * Note the transposes of A and B!
 */
void statePenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_8_26_8(statePenaltyMPC_FLOAT *A, statePenaltyMPC_FLOAT *x, statePenaltyMPC_FLOAT *B, statePenaltyMPC_FLOAT *y, statePenaltyMPC_FLOAT *z)
{
	int i;
	int j;
	int k = 0;
	for( i=0; i<8; i++ ){
		z[i] = 0;
		for( j=0; j<8; j++ ){
			z[i] += A[k++]*x[j];
		}
		z[i] += B[i]*y[i];
	}
	for( i=8 ;i<26; i++ ){
		z[i] = 0;
		for( j=0; j<8; j++ ){
			z[i] += A[k++]*x[j];
		}
	}
}


/*
 * Matrix vector multiplication y = M'*x where M is of size [8 x 8]
 * and stored in diagzero format. Note the transpose of M!
 */
void statePenaltyMPC_LA_DIAGZERO_MTVM_8_8(statePenaltyMPC_FLOAT *M, statePenaltyMPC_FLOAT *x, statePenaltyMPC_FLOAT *y)
{
	int i;
	for( i=0; i<8; i++ ){
		y[i] = M[i]*x[i];
	}
}


/*
 * Vector subtraction and addition.
 *	 Input: five vectors t, tidx, u, v, w and two scalars z and r
 *	 Output: y = t(tidx) - u + w
 *           z = z - v'*x;
 *           r = max([norm(y,inf), z]);
 * for vectors of length 26. Output z is of course scalar.
 */
void statePenaltyMPC_LA_VSUBADD3_26(statePenaltyMPC_FLOAT* t, statePenaltyMPC_FLOAT* u, int* uidx, statePenaltyMPC_FLOAT* v, statePenaltyMPC_FLOAT* w, statePenaltyMPC_FLOAT* y, statePenaltyMPC_FLOAT* z, statePenaltyMPC_FLOAT* r)
{
	int i;
	statePenaltyMPC_FLOAT norm = *r;
	statePenaltyMPC_FLOAT vx = 0;
	statePenaltyMPC_FLOAT x;
	for( i=0; i<26; i++){
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
 * for vectors of length 10. Output z is of course scalar.
 */
void statePenaltyMPC_LA_VSUBADD2_10(statePenaltyMPC_FLOAT* t, int* tidx, statePenaltyMPC_FLOAT* u, statePenaltyMPC_FLOAT* v, statePenaltyMPC_FLOAT* w, statePenaltyMPC_FLOAT* y, statePenaltyMPC_FLOAT* z, statePenaltyMPC_FLOAT* r)
{
	int i;
	statePenaltyMPC_FLOAT norm = *r;
	statePenaltyMPC_FLOAT vx = 0;
	statePenaltyMPC_FLOAT x;
	for( i=0; i<10; i++){
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
 * for vectors of length 8. Output z is of course scalar.
 */
void statePenaltyMPC_LA_VSUBADD3_8(statePenaltyMPC_FLOAT* t, statePenaltyMPC_FLOAT* u, int* uidx, statePenaltyMPC_FLOAT* v, statePenaltyMPC_FLOAT* w, statePenaltyMPC_FLOAT* y, statePenaltyMPC_FLOAT* z, statePenaltyMPC_FLOAT* r)
{
	int i;
	statePenaltyMPC_FLOAT norm = *r;
	statePenaltyMPC_FLOAT vx = 0;
	statePenaltyMPC_FLOAT x;
	for( i=0; i<8; i++){
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
 * for vectors of length 8. Output z is of course scalar.
 */
void statePenaltyMPC_LA_VSUBADD2_8(statePenaltyMPC_FLOAT* t, int* tidx, statePenaltyMPC_FLOAT* u, statePenaltyMPC_FLOAT* v, statePenaltyMPC_FLOAT* w, statePenaltyMPC_FLOAT* y, statePenaltyMPC_FLOAT* z, statePenaltyMPC_FLOAT* r)
{
	int i;
	statePenaltyMPC_FLOAT norm = *r;
	statePenaltyMPC_FLOAT vx = 0;
	statePenaltyMPC_FLOAT x;
	for( i=0; i<8; i++){
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
 * Special function for box constraints of length 26
 * Returns also L/S, a value that is often used elsewhere.
 */
void statePenaltyMPC_LA_INEQ_B_GRAD_26_26_10(statePenaltyMPC_FLOAT *lu, statePenaltyMPC_FLOAT *su, statePenaltyMPC_FLOAT *ru, statePenaltyMPC_FLOAT *ll, statePenaltyMPC_FLOAT *sl, statePenaltyMPC_FLOAT *rl, int* lbIdx, int* ubIdx, statePenaltyMPC_FLOAT *grad, statePenaltyMPC_FLOAT *lubysu, statePenaltyMPC_FLOAT *llbysl)
{
	int i;
	for( i=0; i<26; i++ ){
		grad[i] = 0;
	}
	for( i=0; i<26; i++ ){		
		llbysl[i] = ll[i] / sl[i];
		grad[lbIdx[i]] -= llbysl[i]*rl[i];
	}
	for( i=0; i<10; i++ ){
		lubysu[i] = lu[i] / su[i];
		grad[ubIdx[i]] += lubysu[i]*ru[i];
	}
}


/*
 * Computes inequality constraints gradient-
 * Special function for box constraints of length 8
 * Returns also L/S, a value that is often used elsewhere.
 */
void statePenaltyMPC_LA_INEQ_B_GRAD_8_8_8(statePenaltyMPC_FLOAT *lu, statePenaltyMPC_FLOAT *su, statePenaltyMPC_FLOAT *ru, statePenaltyMPC_FLOAT *ll, statePenaltyMPC_FLOAT *sl, statePenaltyMPC_FLOAT *rl, int* lbIdx, int* ubIdx, statePenaltyMPC_FLOAT *grad, statePenaltyMPC_FLOAT *lubysu, statePenaltyMPC_FLOAT *llbysl)
{
	int i;
	for( i=0; i<8; i++ ){
		grad[i] = 0;
	}
	for( i=0; i<8; i++ ){		
		llbysl[i] = ll[i] / sl[i];
		grad[lbIdx[i]] -= llbysl[i]*rl[i];
	}
	for( i=0; i<8; i++ ){
		lubysu[i] = lu[i] / su[i];
		grad[ubIdx[i]] += lubysu[i]*ru[i];
	}
}


/*
 * Addition of three vectors  z = u + w + v
 * of length 372.
 */
void statePenaltyMPC_LA_VVADD3_372(statePenaltyMPC_FLOAT *u, statePenaltyMPC_FLOAT *v, statePenaltyMPC_FLOAT *w, statePenaltyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<372; i++ ){
		z[i] = u[i] + v[i] + w[i];
	}
}


/*
 * Special function to compute the diagonal cholesky factorization of the 
 * positive definite augmented Hessian for block size 26.
 *
 * Inputs: - H = diagonal cost Hessian in diagonal storage format
 *         - llbysl = L / S of lower bounds
 *         - lubysu = L / S of upper bounds
 *
 * Output: Phi = sqrt(H + diag(llbysl) + diag(lubysu))
 * where Phi is stored in diagonal storage format
 */
void statePenaltyMPC_LA_DIAG_CHOL_LBUB_26_26_10(statePenaltyMPC_FLOAT *H, statePenaltyMPC_FLOAT *llbysl, int* lbIdx, statePenaltyMPC_FLOAT *lubysu, int* ubIdx, statePenaltyMPC_FLOAT *Phi)


{
	int i;
	
	/* copy  H into PHI */
	for( i=0; i<26; i++ ){
		Phi[i] = H[i];
	}

	/* add llbysl onto Phi where necessary */
	for( i=0; i<26; i++ ){
		Phi[lbIdx[i]] += llbysl[i];
	}

	/* add lubysu onto Phi where necessary */
	for( i=0; i<10; i++){
		Phi[ubIdx[i]] +=  lubysu[i];
	}
	
	/* compute cholesky */
	for(i=0; i<26; i++)
	{
#if statePenaltyMPC_SET_PRINTLEVEL > 0 && defined PRINTNUMERICALWARNINGS
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
 * where A is to be computed and is of size [8 x 26],
 * B is given and of size [8 x 26], L is a diagonal
 * matrix of size 8 stored in diagonal matrix 
 * storage format. Note the transpose of L has no impact!
 *
 * Result: A in column major storage format.
 *
 */
void statePenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_8_26(statePenaltyMPC_FLOAT *L, statePenaltyMPC_FLOAT *B, statePenaltyMPC_FLOAT *A)
{
    int i,j;
	 int k = 0;

	for( j=0; j<26; j++){
		for( i=0; i<8; i++){
			A[k] = B[k]/L[j];
			k++;
		}
	}

}


/**
 * Forward substitution for the matrix equation A*L' = B
 * where A is to be computed and is of size [8 x 26],
 * B is given and of size [8 x 26], L is a diagonal
 *  matrix of size 26 stored in diagonal 
 * storage format. Note the transpose of L!
 *
 * Result: A in diagonalzero storage format.
 *
 */
void statePenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_8_26(statePenaltyMPC_FLOAT *L, statePenaltyMPC_FLOAT *B, statePenaltyMPC_FLOAT *A)
{
	int j;
    for( j=0; j<26; j++ ){   
		A[j] = B[j]/L[j];
     }
}


/**
 * Compute C = A*B' where 
 *
 *	size(A) = [8 x 26]
 *  size(B) = [8 x 26] in diagzero format
 * 
 * A and C matrices are stored in column major format.
 * 
 * 
 */
void statePenaltyMPC_LA_DENSE_DIAGZERO_MMTM_8_26_8(statePenaltyMPC_FLOAT *A, statePenaltyMPC_FLOAT *B, statePenaltyMPC_FLOAT *C)
{
    int i, j;
	
	for( i=0; i<8; i++ ){
		for( j=0; j<8; j++){
			C[j*8+i] = B[i*8+j]*A[i];
		}
	}

}


/**
 * Forward substitution to solve L*y = b where L is a
 * diagonal matrix in vector storage format.
 * 
 * The dimensions involved are 26.
 */
void statePenaltyMPC_LA_DIAG_FORWARDSUB_26(statePenaltyMPC_FLOAT *L, statePenaltyMPC_FLOAT *b, statePenaltyMPC_FLOAT *y)
{
    int i;

    for( i=0; i<26; i++ ){
		y[i] = b[i]/L[i];
    }
}


/*
 * Special function to compute the diagonal cholesky factorization of the 
 * positive definite augmented Hessian for block size 8.
 *
 * Inputs: - H = diagonal cost Hessian in diagonal storage format
 *         - llbysl = L / S of lower bounds
 *         - lubysu = L / S of upper bounds
 *
 * Output: Phi = sqrt(H + diag(llbysl) + diag(lubysu))
 * where Phi is stored in diagonal storage format
 */
void statePenaltyMPC_LA_DIAG_CHOL_ONELOOP_LBUB_8_8_8(statePenaltyMPC_FLOAT *H, statePenaltyMPC_FLOAT *llbysl, int* lbIdx, statePenaltyMPC_FLOAT *lubysu, int* ubIdx, statePenaltyMPC_FLOAT *Phi)


{
	int i;
	
	/* compute cholesky */
	for( i=0; i<8; i++ ){
		Phi[i] = H[i] + llbysl[i] + lubysu[i];

#if statePenaltyMPC_SET_PRINTLEVEL > 0 && defined PRINTNUMERICALWARNINGS
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
 * where A is to be computed and is of size [8 x 8],
 * B is given and of size [8 x 8], L is a diagonal
 *  matrix of size 8 stored in diagonal 
 * storage format. Note the transpose of L!
 *
 * Result: A in diagonalzero storage format.
 *
 */
void statePenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_8_8(statePenaltyMPC_FLOAT *L, statePenaltyMPC_FLOAT *B, statePenaltyMPC_FLOAT *A)
{
	int j;
    for( j=0; j<8; j++ ){   
		A[j] = B[j]/L[j];
     }
}


/**
 * Forward substitution to solve L*y = b where L is a
 * diagonal matrix in vector storage format.
 * 
 * The dimensions involved are 8.
 */
void statePenaltyMPC_LA_DIAG_FORWARDSUB_8(statePenaltyMPC_FLOAT *L, statePenaltyMPC_FLOAT *b, statePenaltyMPC_FLOAT *y)
{
    int i;

    for( i=0; i<8; i++ ){
		y[i] = b[i]/L[i];
    }
}


/**
 * Compute L = B*B', where L is lower triangular of size NXp1
 * and B is a diagzero matrix of size [8 x 26] in column
 * storage format.
 * 
 */
void statePenaltyMPC_LA_DIAGZERO_MMT_8(statePenaltyMPC_FLOAT *B, statePenaltyMPC_FLOAT *L)
{
    int i, ii, di;
    
    ii = 0; di = 0;
    for( i=0; i<8; i++ ){        
		L[ii+i] = B[i]*B[i];
        ii += ++di;
    }
}


/* 
 * Computes r = b - B*u
 * B is stored in diagzero format
 */
void statePenaltyMPC_LA_DIAGZERO_MVMSUB7_8(statePenaltyMPC_FLOAT *B, statePenaltyMPC_FLOAT *u, statePenaltyMPC_FLOAT *b, statePenaltyMPC_FLOAT *r)
{
	int i;

	for( i=0; i<8; i++ ){
		r[i] = b[i] - B[i]*u[i];
	}	
	
}


/**
 * Compute L = A*A' + B*B', where L is lower triangular of size NXp1
 * and A is a dense matrix of size [8 x 26] in column
 * storage format, and B is of size [8 x 26] diagonalzero
 * storage format.
 * 
 * THIS ONE HAS THE WORST ACCES PATTERN POSSIBLE. 
 * POSSIBKE FIX: PUT A AND B INTO ROW MAJOR FORMAT FIRST.
 * 
 */
void statePenaltyMPC_LA_DENSE_DIAGZERO_MMT2_8_26_26(statePenaltyMPC_FLOAT *A, statePenaltyMPC_FLOAT *B, statePenaltyMPC_FLOAT *L)
{
    int i, j, k, ii, di;
    statePenaltyMPC_FLOAT ltemp;
    
    ii = 0; di = 0;
    for( i=0; i<8; i++ ){        
        for( j=0; j<=i; j++ ){
            ltemp = 0; 
            for( k=0; k<26; k++ ){
                ltemp += A[k*8+i]*A[k*8+j];
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
void statePenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_8_26_26(statePenaltyMPC_FLOAT *A, statePenaltyMPC_FLOAT *x, statePenaltyMPC_FLOAT *B, statePenaltyMPC_FLOAT *u, statePenaltyMPC_FLOAT *b, statePenaltyMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;

	for( i=0; i<8; i++ ){
		r[i] = b[i] - A[k++]*x[0] - B[i]*u[i];
	}	

	for( j=1; j<26; j++ ){		
		for( i=0; i<8; i++ ){
			r[i] -= A[k++]*x[j];
		}
	}
	
}


/**
 * Compute L = A*A' + B*B', where L is lower triangular of size NXp1
 * and A is a dense matrix of size [8 x 26] in column
 * storage format, and B is of size [8 x 8] diagonalzero
 * storage format.
 * 
 * THIS ONE HAS THE WORST ACCES PATTERN POSSIBLE. 
 * POSSIBKE FIX: PUT A AND B INTO ROW MAJOR FORMAT FIRST.
 * 
 */
void statePenaltyMPC_LA_DENSE_DIAGZERO_MMT2_8_26_8(statePenaltyMPC_FLOAT *A, statePenaltyMPC_FLOAT *B, statePenaltyMPC_FLOAT *L)
{
    int i, j, k, ii, di;
    statePenaltyMPC_FLOAT ltemp;
    
    ii = 0; di = 0;
    for( i=0; i<8; i++ ){        
        for( j=0; j<=i; j++ ){
            ltemp = 0; 
            for( k=0; k<26; k++ ){
                ltemp += A[k*8+i]*A[k*8+j];
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
void statePenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_8_26_8(statePenaltyMPC_FLOAT *A, statePenaltyMPC_FLOAT *x, statePenaltyMPC_FLOAT *B, statePenaltyMPC_FLOAT *u, statePenaltyMPC_FLOAT *b, statePenaltyMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;

	for( i=0; i<8; i++ ){
		r[i] = b[i] - A[k++]*x[0] - B[i]*u[i];
	}	

	for( j=1; j<26; j++ ){		
		for( i=0; i<8; i++ ){
			r[i] -= A[k++]*x[j];
		}
	}
	
}


/**
 * Cholesky factorization as above, but working on a matrix in 
 * lower triangular storage format of size 8 and outputting
 * the Cholesky factor to matrix L in lower triangular format.
 */
void statePenaltyMPC_LA_DENSE_CHOL_8(statePenaltyMPC_FLOAT *A, statePenaltyMPC_FLOAT *L)
{
    int i, j, k, di, dj;
	 int ii, jj;

    statePenaltyMPC_FLOAT l;
    statePenaltyMPC_FLOAT Mii;

	/* copy A to L first and then operate on L */
	/* COULD BE OPTIMIZED */
	ii=0; di=0;
	for( i=0; i<8; i++ ){
		for( j=0; j<=i; j++ ){
			L[ii+j] = A[ii+j];
		}
		ii += ++di;
	}    
	
	/* factor L */
	ii=0; di=0;
    for( i=0; i<8; i++ ){
        l = 0;
        for( k=0; k<i; k++ ){
            l += L[ii+k]*L[ii+k];
        }        
        
        Mii = L[ii+i] - l;
        
#if statePenaltyMPC_SET_PRINTLEVEL > 0 && defined PRINTNUMERICALWARNINGS
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
        for( j=i+1; j<8; j++ ){
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
 * The dimensions involved are 8.
 */
void statePenaltyMPC_LA_DENSE_FORWARDSUB_8(statePenaltyMPC_FLOAT *L, statePenaltyMPC_FLOAT *b, statePenaltyMPC_FLOAT *y)
{
    int i,j,ii,di;
    statePenaltyMPC_FLOAT yel;
            
    ii = 0; di = 0;
    for( i=0; i<8; i++ ){
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
 * where A is to be computed and is of size [8 x 8],
 * B is given and of size [8 x 8], L is a lower tri-
 * angular matrix of size 8 stored in lower triangular 
 * storage format. Note the transpose of L AND B!
 *
 * Result: A in column major storage format.
 *
 */
void statePenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_8_8(statePenaltyMPC_FLOAT *L, statePenaltyMPC_FLOAT *B, statePenaltyMPC_FLOAT *A)
{
    int i,j,k,ii,di;
    statePenaltyMPC_FLOAT a;
    
    ii=0; di=0;
    for( j=0; j<8; j++ ){        
        for( i=0; i<8; i++ ){
            a = B[i*8+j];
            for( k=0; k<j; k++ ){
                a -= A[k*8+i]*L[ii+k];
            }    

			/* saturate for numerical stability */
			a = MIN(a, BIGM);
			a = MAX(a, -BIGM); 

			A[j*8+i] = a/L[ii+j];			
        }
        ii += ++di;
    }
}


/**
 * Compute L = L - A*A', where L is lower triangular of size 8
 * and A is a dense matrix of size [8 x 8] in column
 * storage format.
 * 
 * THIS ONE HAS THE WORST ACCES PATTERN POSSIBLE. 
 * POSSIBKE FIX: PUT A INTO ROW MAJOR FORMAT FIRST.
 * 
 */
void statePenaltyMPC_LA_DENSE_MMTSUB_8_8(statePenaltyMPC_FLOAT *A, statePenaltyMPC_FLOAT *L)
{
    int i, j, k, ii, di;
    statePenaltyMPC_FLOAT ltemp;
    
    ii = 0; di = 0;
    for( i=0; i<8; i++ ){        
        for( j=0; j<=i; j++ ){
            ltemp = 0; 
            for( k=0; k<8; k++ ){
                ltemp += A[k*8+i]*A[k*8+j];
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
void statePenaltyMPC_LA_DENSE_MVMSUB1_8_8(statePenaltyMPC_FLOAT *A, statePenaltyMPC_FLOAT *x, statePenaltyMPC_FLOAT *b, statePenaltyMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;

	for( i=0; i<8; i++ ){
		r[i] = b[i] - A[k++]*x[0];
	}	
	for( j=1; j<8; j++ ){		
		for( i=0; i<8; i++ ){
			r[i] -= A[k++]*x[j];
		}
	}
}


/**
 * Backward Substitution to solve L^T*x = y where L is a
 * lower triangular matrix in triangular storage format.
 * 
 * All involved dimensions are 8.
 */
void statePenaltyMPC_LA_DENSE_BACKWARDSUB_8(statePenaltyMPC_FLOAT *L, statePenaltyMPC_FLOAT *y, statePenaltyMPC_FLOAT *x)
{
    int i, ii, di, j, jj, dj;
    statePenaltyMPC_FLOAT xel;    
	int start = 28;
    
    /* now solve L^T*x = y by backward substitution */
    ii = start; di = 7;
    for( i=7; i>=0; i-- ){        
        xel = y[i];        
        jj = start; dj = 7;
        for( j=7; j>i; j-- ){
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
 * Matrix vector multiplication y = b - M'*x where M is of size [8 x 8]
 * and stored in column major format. Note the transpose of M!
 */
void statePenaltyMPC_LA_DENSE_MTVMSUB_8_8(statePenaltyMPC_FLOAT *A, statePenaltyMPC_FLOAT *x, statePenaltyMPC_FLOAT *b, statePenaltyMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0; 
	for( i=0; i<8; i++ ){
		r[i] = b[i];
		for( j=0; j<8; j++ ){
			r[i] -= A[k++]*x[j];
		}
	}
}


/*
 * Vector subtraction z = -x - y for vectors of length 372.
 */
void statePenaltyMPC_LA_VSUB2_372(statePenaltyMPC_FLOAT *x, statePenaltyMPC_FLOAT *y, statePenaltyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<372; i++){
		z[i] = -x[i] - y[i];
	}
}


/**
 * Forward-Backward-Substitution to solve L*L^T*x = b where L is a
 * diagonal matrix of size 26 in vector
 * storage format.
 */
void statePenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_26(statePenaltyMPC_FLOAT *L, statePenaltyMPC_FLOAT *b, statePenaltyMPC_FLOAT *x)
{
    int i;
            
    /* solve Ly = b by forward and backward substitution */
    for( i=0; i<26; i++ ){
		x[i] = b[i]/(L[i]*L[i]);
    }
    
}


/**
 * Forward-Backward-Substitution to solve L*L^T*x = b where L is a
 * diagonal matrix of size 8 in vector
 * storage format.
 */
void statePenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_8(statePenaltyMPC_FLOAT *L, statePenaltyMPC_FLOAT *b, statePenaltyMPC_FLOAT *x)
{
    int i;
            
    /* solve Ly = b by forward and backward substitution */
    for( i=0; i<8; i++ ){
		x[i] = b[i]/(L[i]*L[i]);
    }
    
}


/*
 * Vector subtraction z = x(xidx) - y where y, z and xidx are of length 26,
 * and x has length 26 and is indexed through yidx.
 */
void statePenaltyMPC_LA_VSUB_INDEXED_26(statePenaltyMPC_FLOAT *x, int* xidx, statePenaltyMPC_FLOAT *y, statePenaltyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<26; i++){
		z[i] = x[xidx[i]] - y[i];
	}
}


/*
 * Vector subtraction x = -u.*v - w for vectors of length 26.
 */
void statePenaltyMPC_LA_VSUB3_26(statePenaltyMPC_FLOAT *u, statePenaltyMPC_FLOAT *v, statePenaltyMPC_FLOAT *w, statePenaltyMPC_FLOAT *x)
{
	int i;
	for( i=0; i<26; i++){
		x[i] = -u[i]*v[i] - w[i];
	}
}


/*
 * Vector subtraction z = -x - y(yidx) where y is of length 26
 * and z, x and yidx are of length 10.
 */
void statePenaltyMPC_LA_VSUB2_INDEXED_10(statePenaltyMPC_FLOAT *x, statePenaltyMPC_FLOAT *y, int* yidx, statePenaltyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<10; i++){
		z[i] = -x[i] - y[yidx[i]];
	}
}


/*
 * Vector subtraction x = -u.*v - w for vectors of length 10.
 */
void statePenaltyMPC_LA_VSUB3_10(statePenaltyMPC_FLOAT *u, statePenaltyMPC_FLOAT *v, statePenaltyMPC_FLOAT *w, statePenaltyMPC_FLOAT *x)
{
	int i;
	for( i=0; i<10; i++){
		x[i] = -u[i]*v[i] - w[i];
	}
}


/*
 * Vector subtraction z = x(xidx) - y where y, z and xidx are of length 8,
 * and x has length 8 and is indexed through yidx.
 */
void statePenaltyMPC_LA_VSUB_INDEXED_8(statePenaltyMPC_FLOAT *x, int* xidx, statePenaltyMPC_FLOAT *y, statePenaltyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<8; i++){
		z[i] = x[xidx[i]] - y[i];
	}
}


/*
 * Vector subtraction x = -u.*v - w for vectors of length 8.
 */
void statePenaltyMPC_LA_VSUB3_8(statePenaltyMPC_FLOAT *u, statePenaltyMPC_FLOAT *v, statePenaltyMPC_FLOAT *w, statePenaltyMPC_FLOAT *x)
{
	int i;
	for( i=0; i<8; i++){
		x[i] = -u[i]*v[i] - w[i];
	}
}


/*
 * Vector subtraction z = -x - y(yidx) where y is of length 8
 * and z, x and yidx are of length 8.
 */
void statePenaltyMPC_LA_VSUB2_INDEXED_8(statePenaltyMPC_FLOAT *x, statePenaltyMPC_FLOAT *y, int* yidx, statePenaltyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<8; i++){
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
 * statePenaltyMPC_NOPROGRESS (should be negative).
 */
int statePenaltyMPC_LINESEARCH_BACKTRACKING_AFFINE(statePenaltyMPC_FLOAT *l, statePenaltyMPC_FLOAT *s, statePenaltyMPC_FLOAT *dl, statePenaltyMPC_FLOAT *ds, statePenaltyMPC_FLOAT *a, statePenaltyMPC_FLOAT *mu_aff)
{
    int i;
	int lsIt=1;    
    statePenaltyMPC_FLOAT dltemp;
    statePenaltyMPC_FLOAT dstemp;
    statePenaltyMPC_FLOAT mya = 1.0;
    statePenaltyMPC_FLOAT mymu;
        
    while( 1 ){                        

        /* 
         * Compute both snew and wnew together.
         * We compute also mu_affine along the way here, as the
         * values might be in registers, so it should be cheaper.
         */
        mymu = 0;
        for( i=0; i<520; i++ ){
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
        if( i == 520 ){
            break;
        } else {
            mya *= statePenaltyMPC_SET_LS_SCALE_AFF;
            if( mya < statePenaltyMPC_SET_LS_MINSTEP ){
                return statePenaltyMPC_NOPROGRESS;
            }
        }
    }
    
    /* return new values and iteration counter */
    *a = mya;
    *mu_aff = mymu / (statePenaltyMPC_FLOAT)520;
    return lsIt;
}


/*
 * Vector subtraction x = u.*v - a where a is a scalar
*  and x,u,v are vectors of length 520.
 */
void statePenaltyMPC_LA_VSUB5_520(statePenaltyMPC_FLOAT *u, statePenaltyMPC_FLOAT *v, statePenaltyMPC_FLOAT a, statePenaltyMPC_FLOAT *x)
{
	int i;
	for( i=0; i<520; i++){
		x[i] = u[i]*v[i] - a;
	}
}


/*
 * Computes x=0; x(uidx) += u/su; x(vidx) -= v/sv where x is of length 26,
 * u, su, uidx are of length 10 and v, sv, vidx are of length 26.
 */
void statePenaltyMPC_LA_VSUB6_INDEXED_26_10_26(statePenaltyMPC_FLOAT *u, statePenaltyMPC_FLOAT *su, int* uidx, statePenaltyMPC_FLOAT *v, statePenaltyMPC_FLOAT *sv, int* vidx, statePenaltyMPC_FLOAT *x)
{
	int i;
	for( i=0; i<26; i++ ){
		x[i] = 0;
	}
	for( i=0; i<10; i++){
		x[uidx[i]] += u[i]/su[i];
	}
	for( i=0; i<26; i++){
		x[vidx[i]] -= v[i]/sv[i];
	}
}


/* 
 * Computes r =  B*u
 * where B is stored in diagzero format
 */
void statePenaltyMPC_LA_DIAGZERO_MVM_8(statePenaltyMPC_FLOAT *B, statePenaltyMPC_FLOAT *u, statePenaltyMPC_FLOAT *r)
{
	int i;

	for( i=0; i<8; i++ ){
		r[i] = B[i]*u[i];
	}	
	
}


/* 
 * Computes r = A*x + B*u
 * where A is stored in column major format
 * and B is stored in diagzero format
 */
void statePenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_8_26_26(statePenaltyMPC_FLOAT *A, statePenaltyMPC_FLOAT *x, statePenaltyMPC_FLOAT *B, statePenaltyMPC_FLOAT *u, statePenaltyMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;

	for( i=0; i<8; i++ ){
		r[i] = A[k++]*x[0] + B[i]*u[i];
	}	

	for( j=1; j<26; j++ ){		
		for( i=0; i<8; i++ ){
			r[i] += A[k++]*x[j];
		}
	}
	
}


/*
 * Computes x=0; x(uidx) += u/su; x(vidx) -= v/sv where x is of length 8,
 * u, su, uidx are of length 8 and v, sv, vidx are of length 8.
 */
void statePenaltyMPC_LA_VSUB6_INDEXED_8_8_8(statePenaltyMPC_FLOAT *u, statePenaltyMPC_FLOAT *su, int* uidx, statePenaltyMPC_FLOAT *v, statePenaltyMPC_FLOAT *sv, int* vidx, statePenaltyMPC_FLOAT *x)
{
	int i;
	for( i=0; i<8; i++ ){
		x[i] = 0;
	}
	for( i=0; i<8; i++){
		x[uidx[i]] += u[i]/su[i];
	}
	for( i=0; i<8; i++){
		x[vidx[i]] -= v[i]/sv[i];
	}
}


/* 
 * Computes r = A*x + B*u
 * where A is stored in column major format
 * and B is stored in diagzero format
 */
void statePenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_8_26_8(statePenaltyMPC_FLOAT *A, statePenaltyMPC_FLOAT *x, statePenaltyMPC_FLOAT *B, statePenaltyMPC_FLOAT *u, statePenaltyMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;

	for( i=0; i<8; i++ ){
		r[i] = A[k++]*x[0] + B[i]*u[i];
	}	

	for( j=1; j<26; j++ ){		
		for( i=0; i<8; i++ ){
			r[i] += A[k++]*x[j];
		}
	}
	
}


/*
 * Vector subtraction z = x - y for vectors of length 372.
 */
void statePenaltyMPC_LA_VSUB_372(statePenaltyMPC_FLOAT *x, statePenaltyMPC_FLOAT *y, statePenaltyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<372; i++){
		z[i] = x[i] - y[i];
	}
}


/** 
 * Computes z = -r./s - u.*y(y)
 * where all vectors except of y are of length 26 (length of y >= 26).
 */
void statePenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_26(statePenaltyMPC_FLOAT *r, statePenaltyMPC_FLOAT *s, statePenaltyMPC_FLOAT *u, statePenaltyMPC_FLOAT *y, int* yidx, statePenaltyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<26; i++ ){
		z[i] = -r[i]/s[i] - u[i]*y[yidx[i]];
	}
}


/** 
 * Computes z = -r./s + u.*y(y)
 * where all vectors except of y are of length 10 (length of y >= 10).
 */
void statePenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_10(statePenaltyMPC_FLOAT *r, statePenaltyMPC_FLOAT *s, statePenaltyMPC_FLOAT *u, statePenaltyMPC_FLOAT *y, int* yidx, statePenaltyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<10; i++ ){
		z[i] = -r[i]/s[i] + u[i]*y[yidx[i]];
	}
}


/** 
 * Computes z = -r./s - u.*y(y)
 * where all vectors except of y are of length 8 (length of y >= 8).
 */
void statePenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_8(statePenaltyMPC_FLOAT *r, statePenaltyMPC_FLOAT *s, statePenaltyMPC_FLOAT *u, statePenaltyMPC_FLOAT *y, int* yidx, statePenaltyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<8; i++ ){
		z[i] = -r[i]/s[i] - u[i]*y[yidx[i]];
	}
}


/** 
 * Computes z = -r./s + u.*y(y)
 * where all vectors except of y are of length 8 (length of y >= 8).
 */
void statePenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_8(statePenaltyMPC_FLOAT *r, statePenaltyMPC_FLOAT *s, statePenaltyMPC_FLOAT *u, statePenaltyMPC_FLOAT *y, int* yidx, statePenaltyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<8; i++ ){
		z[i] = -r[i]/s[i] + u[i]*y[yidx[i]];
	}
}


/*
 * Computes ds = -l.\(r + s.*dl) for vectors of length 520.
 */
void statePenaltyMPC_LA_VSUB7_520(statePenaltyMPC_FLOAT *l, statePenaltyMPC_FLOAT *r, statePenaltyMPC_FLOAT *s, statePenaltyMPC_FLOAT *dl, statePenaltyMPC_FLOAT *ds)
{
	int i;
	for( i=0; i<520; i++){
		ds[i] = -(r[i] + s[i]*dl[i])/l[i];
	}
}


/*
 * Vector addition x = x + y for vectors of length 372.
 */
void statePenaltyMPC_LA_VADD_372(statePenaltyMPC_FLOAT *x, statePenaltyMPC_FLOAT *y)
{
	int i;
	for( i=0; i<372; i++){
		x[i] += y[i];
	}
}


/*
 * Vector addition x = x + y for vectors of length 120.
 */
void statePenaltyMPC_LA_VADD_120(statePenaltyMPC_FLOAT *x, statePenaltyMPC_FLOAT *y)
{
	int i;
	for( i=0; i<120; i++){
		x[i] += y[i];
	}
}


/*
 * Vector addition x = x + y for vectors of length 520.
 */
void statePenaltyMPC_LA_VADD_520(statePenaltyMPC_FLOAT *x, statePenaltyMPC_FLOAT *y)
{
	int i;
	for( i=0; i<520; i++){
		x[i] += y[i];
	}
}


/**
 * Backtracking line search for combined predictor/corrector step.
 * Update on variables with safety factor gamma (to keep us away from
 * boundary).
 */
int statePenaltyMPC_LINESEARCH_BACKTRACKING_COMBINED(statePenaltyMPC_FLOAT *z, statePenaltyMPC_FLOAT *v, statePenaltyMPC_FLOAT *l, statePenaltyMPC_FLOAT *s, statePenaltyMPC_FLOAT *dz, statePenaltyMPC_FLOAT *dv, statePenaltyMPC_FLOAT *dl, statePenaltyMPC_FLOAT *ds, statePenaltyMPC_FLOAT *a, statePenaltyMPC_FLOAT *mu)
{
    int i, lsIt=1;       
    statePenaltyMPC_FLOAT dltemp;
    statePenaltyMPC_FLOAT dstemp;    
    statePenaltyMPC_FLOAT a_gamma;
            
    *a = 1.0;
    while( 1 ){                        

        /* check whether search criterion is fulfilled */
        for( i=0; i<520; i++ ){
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
        if( i == 520 ){
            break;
        } else {
            *a *= statePenaltyMPC_SET_LS_SCALE;
            if( *a < statePenaltyMPC_SET_LS_MINSTEP ){
                return statePenaltyMPC_NOPROGRESS;
            }
        }
    }
    
    /* update variables with safety margin */
    a_gamma = (*a)*statePenaltyMPC_SET_LS_MAXSTEP;
    
    /* primal variables */
    for( i=0; i<372; i++ ){
        z[i] += a_gamma*dz[i];
    }
    
    /* equality constraint multipliers */
    for( i=0; i<120; i++ ){
        v[i] += a_gamma*dv[i];
    }
    
    /* inequality constraint multipliers & slacks, also update mu */
    *mu = 0;
    for( i=0; i<520; i++ ){
        dltemp = l[i] + a_gamma*dl[i]; l[i] = dltemp;
        dstemp = s[i] + a_gamma*ds[i]; s[i] = dstemp;
        *mu += dltemp*dstemp;
    }
    
    *a = a_gamma;
    *mu /= (statePenaltyMPC_FLOAT)520;
    return lsIt;
}




/* VARIABLE DEFINITIONS ------------------------------------------------ */
statePenaltyMPC_FLOAT statePenaltyMPC_z[372];
statePenaltyMPC_FLOAT statePenaltyMPC_v[120];
statePenaltyMPC_FLOAT statePenaltyMPC_dz_aff[372];
statePenaltyMPC_FLOAT statePenaltyMPC_dv_aff[120];
statePenaltyMPC_FLOAT statePenaltyMPC_grad_cost[372];
statePenaltyMPC_FLOAT statePenaltyMPC_grad_eq[372];
statePenaltyMPC_FLOAT statePenaltyMPC_rd[372];
statePenaltyMPC_FLOAT statePenaltyMPC_l[520];
statePenaltyMPC_FLOAT statePenaltyMPC_s[520];
statePenaltyMPC_FLOAT statePenaltyMPC_lbys[520];
statePenaltyMPC_FLOAT statePenaltyMPC_dl_aff[520];
statePenaltyMPC_FLOAT statePenaltyMPC_ds_aff[520];
statePenaltyMPC_FLOAT statePenaltyMPC_dz_cc[372];
statePenaltyMPC_FLOAT statePenaltyMPC_dv_cc[120];
statePenaltyMPC_FLOAT statePenaltyMPC_dl_cc[520];
statePenaltyMPC_FLOAT statePenaltyMPC_ds_cc[520];
statePenaltyMPC_FLOAT statePenaltyMPC_ccrhs[520];
statePenaltyMPC_FLOAT statePenaltyMPC_grad_ineq[372];
statePenaltyMPC_FLOAT* statePenaltyMPC_z00 = statePenaltyMPC_z + 0;
statePenaltyMPC_FLOAT* statePenaltyMPC_dzaff00 = statePenaltyMPC_dz_aff + 0;
statePenaltyMPC_FLOAT* statePenaltyMPC_dzcc00 = statePenaltyMPC_dz_cc + 0;
statePenaltyMPC_FLOAT* statePenaltyMPC_rd00 = statePenaltyMPC_rd + 0;
statePenaltyMPC_FLOAT statePenaltyMPC_Lbyrd00[26];
statePenaltyMPC_FLOAT* statePenaltyMPC_grad_cost00 = statePenaltyMPC_grad_cost + 0;
statePenaltyMPC_FLOAT* statePenaltyMPC_grad_eq00 = statePenaltyMPC_grad_eq + 0;
statePenaltyMPC_FLOAT* statePenaltyMPC_grad_ineq00 = statePenaltyMPC_grad_ineq + 0;
statePenaltyMPC_FLOAT statePenaltyMPC_ctv00[26];
statePenaltyMPC_FLOAT* statePenaltyMPC_v00 = statePenaltyMPC_v + 0;
statePenaltyMPC_FLOAT statePenaltyMPC_re00[8];
statePenaltyMPC_FLOAT statePenaltyMPC_beta00[8];
statePenaltyMPC_FLOAT statePenaltyMPC_betacc00[8];
statePenaltyMPC_FLOAT* statePenaltyMPC_dvaff00 = statePenaltyMPC_dv_aff + 0;
statePenaltyMPC_FLOAT* statePenaltyMPC_dvcc00 = statePenaltyMPC_dv_cc + 0;
statePenaltyMPC_FLOAT statePenaltyMPC_V00[208];
statePenaltyMPC_FLOAT statePenaltyMPC_Yd00[36];
statePenaltyMPC_FLOAT statePenaltyMPC_Ld00[36];
statePenaltyMPC_FLOAT statePenaltyMPC_yy00[8];
statePenaltyMPC_FLOAT statePenaltyMPC_bmy00[8];
int statePenaltyMPC_lbIdx00[26] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25};
statePenaltyMPC_FLOAT* statePenaltyMPC_llb00 = statePenaltyMPC_l + 0;
statePenaltyMPC_FLOAT* statePenaltyMPC_slb00 = statePenaltyMPC_s + 0;
statePenaltyMPC_FLOAT* statePenaltyMPC_llbbyslb00 = statePenaltyMPC_lbys + 0;
statePenaltyMPC_FLOAT statePenaltyMPC_rilb00[26];
statePenaltyMPC_FLOAT* statePenaltyMPC_dllbaff00 = statePenaltyMPC_dl_aff + 0;
statePenaltyMPC_FLOAT* statePenaltyMPC_dslbaff00 = statePenaltyMPC_ds_aff + 0;
statePenaltyMPC_FLOAT* statePenaltyMPC_dllbcc00 = statePenaltyMPC_dl_cc + 0;
statePenaltyMPC_FLOAT* statePenaltyMPC_dslbcc00 = statePenaltyMPC_ds_cc + 0;
statePenaltyMPC_FLOAT* statePenaltyMPC_ccrhsl00 = statePenaltyMPC_ccrhs + 0;
int statePenaltyMPC_ubIdx00[10] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
statePenaltyMPC_FLOAT* statePenaltyMPC_lub00 = statePenaltyMPC_l + 26;
statePenaltyMPC_FLOAT* statePenaltyMPC_sub00 = statePenaltyMPC_s + 26;
statePenaltyMPC_FLOAT* statePenaltyMPC_lubbysub00 = statePenaltyMPC_lbys + 26;
statePenaltyMPC_FLOAT statePenaltyMPC_riub00[10];
statePenaltyMPC_FLOAT* statePenaltyMPC_dlubaff00 = statePenaltyMPC_dl_aff + 26;
statePenaltyMPC_FLOAT* statePenaltyMPC_dsubaff00 = statePenaltyMPC_ds_aff + 26;
statePenaltyMPC_FLOAT* statePenaltyMPC_dlubcc00 = statePenaltyMPC_dl_cc + 26;
statePenaltyMPC_FLOAT* statePenaltyMPC_dsubcc00 = statePenaltyMPC_ds_cc + 26;
statePenaltyMPC_FLOAT* statePenaltyMPC_ccrhsub00 = statePenaltyMPC_ccrhs + 26;
statePenaltyMPC_FLOAT statePenaltyMPC_Phi00[26];
statePenaltyMPC_FLOAT statePenaltyMPC_D00[26] = {1.0000000000000000E+000, 
1.0000000000000000E+000, 
1.0000000000000000E+000, 
1.0000000000000000E+000, 
1.0000000000000000E+000, 
1.0000000000000000E+000, 
1.0000000000000000E+000, 
1.0000000000000000E+000};
statePenaltyMPC_FLOAT statePenaltyMPC_W00[26];
statePenaltyMPC_FLOAT* statePenaltyMPC_z01 = statePenaltyMPC_z + 26;
statePenaltyMPC_FLOAT* statePenaltyMPC_dzaff01 = statePenaltyMPC_dz_aff + 26;
statePenaltyMPC_FLOAT* statePenaltyMPC_dzcc01 = statePenaltyMPC_dz_cc + 26;
statePenaltyMPC_FLOAT* statePenaltyMPC_rd01 = statePenaltyMPC_rd + 26;
statePenaltyMPC_FLOAT statePenaltyMPC_Lbyrd01[26];
statePenaltyMPC_FLOAT* statePenaltyMPC_grad_cost01 = statePenaltyMPC_grad_cost + 26;
statePenaltyMPC_FLOAT* statePenaltyMPC_grad_eq01 = statePenaltyMPC_grad_eq + 26;
statePenaltyMPC_FLOAT* statePenaltyMPC_grad_ineq01 = statePenaltyMPC_grad_ineq + 26;
statePenaltyMPC_FLOAT statePenaltyMPC_ctv01[26];
statePenaltyMPC_FLOAT* statePenaltyMPC_v01 = statePenaltyMPC_v + 8;
statePenaltyMPC_FLOAT statePenaltyMPC_re01[8];
statePenaltyMPC_FLOAT statePenaltyMPC_beta01[8];
statePenaltyMPC_FLOAT statePenaltyMPC_betacc01[8];
statePenaltyMPC_FLOAT* statePenaltyMPC_dvaff01 = statePenaltyMPC_dv_aff + 8;
statePenaltyMPC_FLOAT* statePenaltyMPC_dvcc01 = statePenaltyMPC_dv_cc + 8;
statePenaltyMPC_FLOAT statePenaltyMPC_V01[208];
statePenaltyMPC_FLOAT statePenaltyMPC_Yd01[36];
statePenaltyMPC_FLOAT statePenaltyMPC_Ld01[36];
statePenaltyMPC_FLOAT statePenaltyMPC_yy01[8];
statePenaltyMPC_FLOAT statePenaltyMPC_bmy01[8];
int statePenaltyMPC_lbIdx01[26] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25};
statePenaltyMPC_FLOAT* statePenaltyMPC_llb01 = statePenaltyMPC_l + 36;
statePenaltyMPC_FLOAT* statePenaltyMPC_slb01 = statePenaltyMPC_s + 36;
statePenaltyMPC_FLOAT* statePenaltyMPC_llbbyslb01 = statePenaltyMPC_lbys + 36;
statePenaltyMPC_FLOAT statePenaltyMPC_rilb01[26];
statePenaltyMPC_FLOAT* statePenaltyMPC_dllbaff01 = statePenaltyMPC_dl_aff + 36;
statePenaltyMPC_FLOAT* statePenaltyMPC_dslbaff01 = statePenaltyMPC_ds_aff + 36;
statePenaltyMPC_FLOAT* statePenaltyMPC_dllbcc01 = statePenaltyMPC_dl_cc + 36;
statePenaltyMPC_FLOAT* statePenaltyMPC_dslbcc01 = statePenaltyMPC_ds_cc + 36;
statePenaltyMPC_FLOAT* statePenaltyMPC_ccrhsl01 = statePenaltyMPC_ccrhs + 36;
int statePenaltyMPC_ubIdx01[10] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
statePenaltyMPC_FLOAT* statePenaltyMPC_lub01 = statePenaltyMPC_l + 62;
statePenaltyMPC_FLOAT* statePenaltyMPC_sub01 = statePenaltyMPC_s + 62;
statePenaltyMPC_FLOAT* statePenaltyMPC_lubbysub01 = statePenaltyMPC_lbys + 62;
statePenaltyMPC_FLOAT statePenaltyMPC_riub01[10];
statePenaltyMPC_FLOAT* statePenaltyMPC_dlubaff01 = statePenaltyMPC_dl_aff + 62;
statePenaltyMPC_FLOAT* statePenaltyMPC_dsubaff01 = statePenaltyMPC_ds_aff + 62;
statePenaltyMPC_FLOAT* statePenaltyMPC_dlubcc01 = statePenaltyMPC_dl_cc + 62;
statePenaltyMPC_FLOAT* statePenaltyMPC_dsubcc01 = statePenaltyMPC_ds_cc + 62;
statePenaltyMPC_FLOAT* statePenaltyMPC_ccrhsub01 = statePenaltyMPC_ccrhs + 62;
statePenaltyMPC_FLOAT statePenaltyMPC_Phi01[26];
statePenaltyMPC_FLOAT statePenaltyMPC_D01[26] = {-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000};
statePenaltyMPC_FLOAT statePenaltyMPC_W01[26];
statePenaltyMPC_FLOAT statePenaltyMPC_Ysd01[64];
statePenaltyMPC_FLOAT statePenaltyMPC_Lsd01[64];
statePenaltyMPC_FLOAT* statePenaltyMPC_z02 = statePenaltyMPC_z + 52;
statePenaltyMPC_FLOAT* statePenaltyMPC_dzaff02 = statePenaltyMPC_dz_aff + 52;
statePenaltyMPC_FLOAT* statePenaltyMPC_dzcc02 = statePenaltyMPC_dz_cc + 52;
statePenaltyMPC_FLOAT* statePenaltyMPC_rd02 = statePenaltyMPC_rd + 52;
statePenaltyMPC_FLOAT statePenaltyMPC_Lbyrd02[26];
statePenaltyMPC_FLOAT* statePenaltyMPC_grad_cost02 = statePenaltyMPC_grad_cost + 52;
statePenaltyMPC_FLOAT* statePenaltyMPC_grad_eq02 = statePenaltyMPC_grad_eq + 52;
statePenaltyMPC_FLOAT* statePenaltyMPC_grad_ineq02 = statePenaltyMPC_grad_ineq + 52;
statePenaltyMPC_FLOAT statePenaltyMPC_ctv02[26];
statePenaltyMPC_FLOAT* statePenaltyMPC_v02 = statePenaltyMPC_v + 16;
statePenaltyMPC_FLOAT statePenaltyMPC_re02[8];
statePenaltyMPC_FLOAT statePenaltyMPC_beta02[8];
statePenaltyMPC_FLOAT statePenaltyMPC_betacc02[8];
statePenaltyMPC_FLOAT* statePenaltyMPC_dvaff02 = statePenaltyMPC_dv_aff + 16;
statePenaltyMPC_FLOAT* statePenaltyMPC_dvcc02 = statePenaltyMPC_dv_cc + 16;
statePenaltyMPC_FLOAT statePenaltyMPC_V02[208];
statePenaltyMPC_FLOAT statePenaltyMPC_Yd02[36];
statePenaltyMPC_FLOAT statePenaltyMPC_Ld02[36];
statePenaltyMPC_FLOAT statePenaltyMPC_yy02[8];
statePenaltyMPC_FLOAT statePenaltyMPC_bmy02[8];
int statePenaltyMPC_lbIdx02[26] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25};
statePenaltyMPC_FLOAT* statePenaltyMPC_llb02 = statePenaltyMPC_l + 72;
statePenaltyMPC_FLOAT* statePenaltyMPC_slb02 = statePenaltyMPC_s + 72;
statePenaltyMPC_FLOAT* statePenaltyMPC_llbbyslb02 = statePenaltyMPC_lbys + 72;
statePenaltyMPC_FLOAT statePenaltyMPC_rilb02[26];
statePenaltyMPC_FLOAT* statePenaltyMPC_dllbaff02 = statePenaltyMPC_dl_aff + 72;
statePenaltyMPC_FLOAT* statePenaltyMPC_dslbaff02 = statePenaltyMPC_ds_aff + 72;
statePenaltyMPC_FLOAT* statePenaltyMPC_dllbcc02 = statePenaltyMPC_dl_cc + 72;
statePenaltyMPC_FLOAT* statePenaltyMPC_dslbcc02 = statePenaltyMPC_ds_cc + 72;
statePenaltyMPC_FLOAT* statePenaltyMPC_ccrhsl02 = statePenaltyMPC_ccrhs + 72;
int statePenaltyMPC_ubIdx02[10] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
statePenaltyMPC_FLOAT* statePenaltyMPC_lub02 = statePenaltyMPC_l + 98;
statePenaltyMPC_FLOAT* statePenaltyMPC_sub02 = statePenaltyMPC_s + 98;
statePenaltyMPC_FLOAT* statePenaltyMPC_lubbysub02 = statePenaltyMPC_lbys + 98;
statePenaltyMPC_FLOAT statePenaltyMPC_riub02[10];
statePenaltyMPC_FLOAT* statePenaltyMPC_dlubaff02 = statePenaltyMPC_dl_aff + 98;
statePenaltyMPC_FLOAT* statePenaltyMPC_dsubaff02 = statePenaltyMPC_ds_aff + 98;
statePenaltyMPC_FLOAT* statePenaltyMPC_dlubcc02 = statePenaltyMPC_dl_cc + 98;
statePenaltyMPC_FLOAT* statePenaltyMPC_dsubcc02 = statePenaltyMPC_ds_cc + 98;
statePenaltyMPC_FLOAT* statePenaltyMPC_ccrhsub02 = statePenaltyMPC_ccrhs + 98;
statePenaltyMPC_FLOAT statePenaltyMPC_Phi02[26];
statePenaltyMPC_FLOAT statePenaltyMPC_W02[26];
statePenaltyMPC_FLOAT statePenaltyMPC_Ysd02[64];
statePenaltyMPC_FLOAT statePenaltyMPC_Lsd02[64];
statePenaltyMPC_FLOAT* statePenaltyMPC_z03 = statePenaltyMPC_z + 78;
statePenaltyMPC_FLOAT* statePenaltyMPC_dzaff03 = statePenaltyMPC_dz_aff + 78;
statePenaltyMPC_FLOAT* statePenaltyMPC_dzcc03 = statePenaltyMPC_dz_cc + 78;
statePenaltyMPC_FLOAT* statePenaltyMPC_rd03 = statePenaltyMPC_rd + 78;
statePenaltyMPC_FLOAT statePenaltyMPC_Lbyrd03[26];
statePenaltyMPC_FLOAT* statePenaltyMPC_grad_cost03 = statePenaltyMPC_grad_cost + 78;
statePenaltyMPC_FLOAT* statePenaltyMPC_grad_eq03 = statePenaltyMPC_grad_eq + 78;
statePenaltyMPC_FLOAT* statePenaltyMPC_grad_ineq03 = statePenaltyMPC_grad_ineq + 78;
statePenaltyMPC_FLOAT statePenaltyMPC_ctv03[26];
statePenaltyMPC_FLOAT* statePenaltyMPC_v03 = statePenaltyMPC_v + 24;
statePenaltyMPC_FLOAT statePenaltyMPC_re03[8];
statePenaltyMPC_FLOAT statePenaltyMPC_beta03[8];
statePenaltyMPC_FLOAT statePenaltyMPC_betacc03[8];
statePenaltyMPC_FLOAT* statePenaltyMPC_dvaff03 = statePenaltyMPC_dv_aff + 24;
statePenaltyMPC_FLOAT* statePenaltyMPC_dvcc03 = statePenaltyMPC_dv_cc + 24;
statePenaltyMPC_FLOAT statePenaltyMPC_V03[208];
statePenaltyMPC_FLOAT statePenaltyMPC_Yd03[36];
statePenaltyMPC_FLOAT statePenaltyMPC_Ld03[36];
statePenaltyMPC_FLOAT statePenaltyMPC_yy03[8];
statePenaltyMPC_FLOAT statePenaltyMPC_bmy03[8];
int statePenaltyMPC_lbIdx03[26] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25};
statePenaltyMPC_FLOAT* statePenaltyMPC_llb03 = statePenaltyMPC_l + 108;
statePenaltyMPC_FLOAT* statePenaltyMPC_slb03 = statePenaltyMPC_s + 108;
statePenaltyMPC_FLOAT* statePenaltyMPC_llbbyslb03 = statePenaltyMPC_lbys + 108;
statePenaltyMPC_FLOAT statePenaltyMPC_rilb03[26];
statePenaltyMPC_FLOAT* statePenaltyMPC_dllbaff03 = statePenaltyMPC_dl_aff + 108;
statePenaltyMPC_FLOAT* statePenaltyMPC_dslbaff03 = statePenaltyMPC_ds_aff + 108;
statePenaltyMPC_FLOAT* statePenaltyMPC_dllbcc03 = statePenaltyMPC_dl_cc + 108;
statePenaltyMPC_FLOAT* statePenaltyMPC_dslbcc03 = statePenaltyMPC_ds_cc + 108;
statePenaltyMPC_FLOAT* statePenaltyMPC_ccrhsl03 = statePenaltyMPC_ccrhs + 108;
int statePenaltyMPC_ubIdx03[10] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
statePenaltyMPC_FLOAT* statePenaltyMPC_lub03 = statePenaltyMPC_l + 134;
statePenaltyMPC_FLOAT* statePenaltyMPC_sub03 = statePenaltyMPC_s + 134;
statePenaltyMPC_FLOAT* statePenaltyMPC_lubbysub03 = statePenaltyMPC_lbys + 134;
statePenaltyMPC_FLOAT statePenaltyMPC_riub03[10];
statePenaltyMPC_FLOAT* statePenaltyMPC_dlubaff03 = statePenaltyMPC_dl_aff + 134;
statePenaltyMPC_FLOAT* statePenaltyMPC_dsubaff03 = statePenaltyMPC_ds_aff + 134;
statePenaltyMPC_FLOAT* statePenaltyMPC_dlubcc03 = statePenaltyMPC_dl_cc + 134;
statePenaltyMPC_FLOAT* statePenaltyMPC_dsubcc03 = statePenaltyMPC_ds_cc + 134;
statePenaltyMPC_FLOAT* statePenaltyMPC_ccrhsub03 = statePenaltyMPC_ccrhs + 134;
statePenaltyMPC_FLOAT statePenaltyMPC_Phi03[26];
statePenaltyMPC_FLOAT statePenaltyMPC_W03[26];
statePenaltyMPC_FLOAT statePenaltyMPC_Ysd03[64];
statePenaltyMPC_FLOAT statePenaltyMPC_Lsd03[64];
statePenaltyMPC_FLOAT* statePenaltyMPC_z04 = statePenaltyMPC_z + 104;
statePenaltyMPC_FLOAT* statePenaltyMPC_dzaff04 = statePenaltyMPC_dz_aff + 104;
statePenaltyMPC_FLOAT* statePenaltyMPC_dzcc04 = statePenaltyMPC_dz_cc + 104;
statePenaltyMPC_FLOAT* statePenaltyMPC_rd04 = statePenaltyMPC_rd + 104;
statePenaltyMPC_FLOAT statePenaltyMPC_Lbyrd04[26];
statePenaltyMPC_FLOAT* statePenaltyMPC_grad_cost04 = statePenaltyMPC_grad_cost + 104;
statePenaltyMPC_FLOAT* statePenaltyMPC_grad_eq04 = statePenaltyMPC_grad_eq + 104;
statePenaltyMPC_FLOAT* statePenaltyMPC_grad_ineq04 = statePenaltyMPC_grad_ineq + 104;
statePenaltyMPC_FLOAT statePenaltyMPC_ctv04[26];
statePenaltyMPC_FLOAT* statePenaltyMPC_v04 = statePenaltyMPC_v + 32;
statePenaltyMPC_FLOAT statePenaltyMPC_re04[8];
statePenaltyMPC_FLOAT statePenaltyMPC_beta04[8];
statePenaltyMPC_FLOAT statePenaltyMPC_betacc04[8];
statePenaltyMPC_FLOAT* statePenaltyMPC_dvaff04 = statePenaltyMPC_dv_aff + 32;
statePenaltyMPC_FLOAT* statePenaltyMPC_dvcc04 = statePenaltyMPC_dv_cc + 32;
statePenaltyMPC_FLOAT statePenaltyMPC_V04[208];
statePenaltyMPC_FLOAT statePenaltyMPC_Yd04[36];
statePenaltyMPC_FLOAT statePenaltyMPC_Ld04[36];
statePenaltyMPC_FLOAT statePenaltyMPC_yy04[8];
statePenaltyMPC_FLOAT statePenaltyMPC_bmy04[8];
int statePenaltyMPC_lbIdx04[26] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25};
statePenaltyMPC_FLOAT* statePenaltyMPC_llb04 = statePenaltyMPC_l + 144;
statePenaltyMPC_FLOAT* statePenaltyMPC_slb04 = statePenaltyMPC_s + 144;
statePenaltyMPC_FLOAT* statePenaltyMPC_llbbyslb04 = statePenaltyMPC_lbys + 144;
statePenaltyMPC_FLOAT statePenaltyMPC_rilb04[26];
statePenaltyMPC_FLOAT* statePenaltyMPC_dllbaff04 = statePenaltyMPC_dl_aff + 144;
statePenaltyMPC_FLOAT* statePenaltyMPC_dslbaff04 = statePenaltyMPC_ds_aff + 144;
statePenaltyMPC_FLOAT* statePenaltyMPC_dllbcc04 = statePenaltyMPC_dl_cc + 144;
statePenaltyMPC_FLOAT* statePenaltyMPC_dslbcc04 = statePenaltyMPC_ds_cc + 144;
statePenaltyMPC_FLOAT* statePenaltyMPC_ccrhsl04 = statePenaltyMPC_ccrhs + 144;
int statePenaltyMPC_ubIdx04[10] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
statePenaltyMPC_FLOAT* statePenaltyMPC_lub04 = statePenaltyMPC_l + 170;
statePenaltyMPC_FLOAT* statePenaltyMPC_sub04 = statePenaltyMPC_s + 170;
statePenaltyMPC_FLOAT* statePenaltyMPC_lubbysub04 = statePenaltyMPC_lbys + 170;
statePenaltyMPC_FLOAT statePenaltyMPC_riub04[10];
statePenaltyMPC_FLOAT* statePenaltyMPC_dlubaff04 = statePenaltyMPC_dl_aff + 170;
statePenaltyMPC_FLOAT* statePenaltyMPC_dsubaff04 = statePenaltyMPC_ds_aff + 170;
statePenaltyMPC_FLOAT* statePenaltyMPC_dlubcc04 = statePenaltyMPC_dl_cc + 170;
statePenaltyMPC_FLOAT* statePenaltyMPC_dsubcc04 = statePenaltyMPC_ds_cc + 170;
statePenaltyMPC_FLOAT* statePenaltyMPC_ccrhsub04 = statePenaltyMPC_ccrhs + 170;
statePenaltyMPC_FLOAT statePenaltyMPC_Phi04[26];
statePenaltyMPC_FLOAT statePenaltyMPC_W04[26];
statePenaltyMPC_FLOAT statePenaltyMPC_Ysd04[64];
statePenaltyMPC_FLOAT statePenaltyMPC_Lsd04[64];
statePenaltyMPC_FLOAT* statePenaltyMPC_z05 = statePenaltyMPC_z + 130;
statePenaltyMPC_FLOAT* statePenaltyMPC_dzaff05 = statePenaltyMPC_dz_aff + 130;
statePenaltyMPC_FLOAT* statePenaltyMPC_dzcc05 = statePenaltyMPC_dz_cc + 130;
statePenaltyMPC_FLOAT* statePenaltyMPC_rd05 = statePenaltyMPC_rd + 130;
statePenaltyMPC_FLOAT statePenaltyMPC_Lbyrd05[26];
statePenaltyMPC_FLOAT* statePenaltyMPC_grad_cost05 = statePenaltyMPC_grad_cost + 130;
statePenaltyMPC_FLOAT* statePenaltyMPC_grad_eq05 = statePenaltyMPC_grad_eq + 130;
statePenaltyMPC_FLOAT* statePenaltyMPC_grad_ineq05 = statePenaltyMPC_grad_ineq + 130;
statePenaltyMPC_FLOAT statePenaltyMPC_ctv05[26];
statePenaltyMPC_FLOAT* statePenaltyMPC_v05 = statePenaltyMPC_v + 40;
statePenaltyMPC_FLOAT statePenaltyMPC_re05[8];
statePenaltyMPC_FLOAT statePenaltyMPC_beta05[8];
statePenaltyMPC_FLOAT statePenaltyMPC_betacc05[8];
statePenaltyMPC_FLOAT* statePenaltyMPC_dvaff05 = statePenaltyMPC_dv_aff + 40;
statePenaltyMPC_FLOAT* statePenaltyMPC_dvcc05 = statePenaltyMPC_dv_cc + 40;
statePenaltyMPC_FLOAT statePenaltyMPC_V05[208];
statePenaltyMPC_FLOAT statePenaltyMPC_Yd05[36];
statePenaltyMPC_FLOAT statePenaltyMPC_Ld05[36];
statePenaltyMPC_FLOAT statePenaltyMPC_yy05[8];
statePenaltyMPC_FLOAT statePenaltyMPC_bmy05[8];
int statePenaltyMPC_lbIdx05[26] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25};
statePenaltyMPC_FLOAT* statePenaltyMPC_llb05 = statePenaltyMPC_l + 180;
statePenaltyMPC_FLOAT* statePenaltyMPC_slb05 = statePenaltyMPC_s + 180;
statePenaltyMPC_FLOAT* statePenaltyMPC_llbbyslb05 = statePenaltyMPC_lbys + 180;
statePenaltyMPC_FLOAT statePenaltyMPC_rilb05[26];
statePenaltyMPC_FLOAT* statePenaltyMPC_dllbaff05 = statePenaltyMPC_dl_aff + 180;
statePenaltyMPC_FLOAT* statePenaltyMPC_dslbaff05 = statePenaltyMPC_ds_aff + 180;
statePenaltyMPC_FLOAT* statePenaltyMPC_dllbcc05 = statePenaltyMPC_dl_cc + 180;
statePenaltyMPC_FLOAT* statePenaltyMPC_dslbcc05 = statePenaltyMPC_ds_cc + 180;
statePenaltyMPC_FLOAT* statePenaltyMPC_ccrhsl05 = statePenaltyMPC_ccrhs + 180;
int statePenaltyMPC_ubIdx05[10] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
statePenaltyMPC_FLOAT* statePenaltyMPC_lub05 = statePenaltyMPC_l + 206;
statePenaltyMPC_FLOAT* statePenaltyMPC_sub05 = statePenaltyMPC_s + 206;
statePenaltyMPC_FLOAT* statePenaltyMPC_lubbysub05 = statePenaltyMPC_lbys + 206;
statePenaltyMPC_FLOAT statePenaltyMPC_riub05[10];
statePenaltyMPC_FLOAT* statePenaltyMPC_dlubaff05 = statePenaltyMPC_dl_aff + 206;
statePenaltyMPC_FLOAT* statePenaltyMPC_dsubaff05 = statePenaltyMPC_ds_aff + 206;
statePenaltyMPC_FLOAT* statePenaltyMPC_dlubcc05 = statePenaltyMPC_dl_cc + 206;
statePenaltyMPC_FLOAT* statePenaltyMPC_dsubcc05 = statePenaltyMPC_ds_cc + 206;
statePenaltyMPC_FLOAT* statePenaltyMPC_ccrhsub05 = statePenaltyMPC_ccrhs + 206;
statePenaltyMPC_FLOAT statePenaltyMPC_Phi05[26];
statePenaltyMPC_FLOAT statePenaltyMPC_W05[26];
statePenaltyMPC_FLOAT statePenaltyMPC_Ysd05[64];
statePenaltyMPC_FLOAT statePenaltyMPC_Lsd05[64];
statePenaltyMPC_FLOAT* statePenaltyMPC_z06 = statePenaltyMPC_z + 156;
statePenaltyMPC_FLOAT* statePenaltyMPC_dzaff06 = statePenaltyMPC_dz_aff + 156;
statePenaltyMPC_FLOAT* statePenaltyMPC_dzcc06 = statePenaltyMPC_dz_cc + 156;
statePenaltyMPC_FLOAT* statePenaltyMPC_rd06 = statePenaltyMPC_rd + 156;
statePenaltyMPC_FLOAT statePenaltyMPC_Lbyrd06[26];
statePenaltyMPC_FLOAT* statePenaltyMPC_grad_cost06 = statePenaltyMPC_grad_cost + 156;
statePenaltyMPC_FLOAT* statePenaltyMPC_grad_eq06 = statePenaltyMPC_grad_eq + 156;
statePenaltyMPC_FLOAT* statePenaltyMPC_grad_ineq06 = statePenaltyMPC_grad_ineq + 156;
statePenaltyMPC_FLOAT statePenaltyMPC_ctv06[26];
statePenaltyMPC_FLOAT* statePenaltyMPC_v06 = statePenaltyMPC_v + 48;
statePenaltyMPC_FLOAT statePenaltyMPC_re06[8];
statePenaltyMPC_FLOAT statePenaltyMPC_beta06[8];
statePenaltyMPC_FLOAT statePenaltyMPC_betacc06[8];
statePenaltyMPC_FLOAT* statePenaltyMPC_dvaff06 = statePenaltyMPC_dv_aff + 48;
statePenaltyMPC_FLOAT* statePenaltyMPC_dvcc06 = statePenaltyMPC_dv_cc + 48;
statePenaltyMPC_FLOAT statePenaltyMPC_V06[208];
statePenaltyMPC_FLOAT statePenaltyMPC_Yd06[36];
statePenaltyMPC_FLOAT statePenaltyMPC_Ld06[36];
statePenaltyMPC_FLOAT statePenaltyMPC_yy06[8];
statePenaltyMPC_FLOAT statePenaltyMPC_bmy06[8];
int statePenaltyMPC_lbIdx06[26] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25};
statePenaltyMPC_FLOAT* statePenaltyMPC_llb06 = statePenaltyMPC_l + 216;
statePenaltyMPC_FLOAT* statePenaltyMPC_slb06 = statePenaltyMPC_s + 216;
statePenaltyMPC_FLOAT* statePenaltyMPC_llbbyslb06 = statePenaltyMPC_lbys + 216;
statePenaltyMPC_FLOAT statePenaltyMPC_rilb06[26];
statePenaltyMPC_FLOAT* statePenaltyMPC_dllbaff06 = statePenaltyMPC_dl_aff + 216;
statePenaltyMPC_FLOAT* statePenaltyMPC_dslbaff06 = statePenaltyMPC_ds_aff + 216;
statePenaltyMPC_FLOAT* statePenaltyMPC_dllbcc06 = statePenaltyMPC_dl_cc + 216;
statePenaltyMPC_FLOAT* statePenaltyMPC_dslbcc06 = statePenaltyMPC_ds_cc + 216;
statePenaltyMPC_FLOAT* statePenaltyMPC_ccrhsl06 = statePenaltyMPC_ccrhs + 216;
int statePenaltyMPC_ubIdx06[10] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
statePenaltyMPC_FLOAT* statePenaltyMPC_lub06 = statePenaltyMPC_l + 242;
statePenaltyMPC_FLOAT* statePenaltyMPC_sub06 = statePenaltyMPC_s + 242;
statePenaltyMPC_FLOAT* statePenaltyMPC_lubbysub06 = statePenaltyMPC_lbys + 242;
statePenaltyMPC_FLOAT statePenaltyMPC_riub06[10];
statePenaltyMPC_FLOAT* statePenaltyMPC_dlubaff06 = statePenaltyMPC_dl_aff + 242;
statePenaltyMPC_FLOAT* statePenaltyMPC_dsubaff06 = statePenaltyMPC_ds_aff + 242;
statePenaltyMPC_FLOAT* statePenaltyMPC_dlubcc06 = statePenaltyMPC_dl_cc + 242;
statePenaltyMPC_FLOAT* statePenaltyMPC_dsubcc06 = statePenaltyMPC_ds_cc + 242;
statePenaltyMPC_FLOAT* statePenaltyMPC_ccrhsub06 = statePenaltyMPC_ccrhs + 242;
statePenaltyMPC_FLOAT statePenaltyMPC_Phi06[26];
statePenaltyMPC_FLOAT statePenaltyMPC_W06[26];
statePenaltyMPC_FLOAT statePenaltyMPC_Ysd06[64];
statePenaltyMPC_FLOAT statePenaltyMPC_Lsd06[64];
statePenaltyMPC_FLOAT* statePenaltyMPC_z07 = statePenaltyMPC_z + 182;
statePenaltyMPC_FLOAT* statePenaltyMPC_dzaff07 = statePenaltyMPC_dz_aff + 182;
statePenaltyMPC_FLOAT* statePenaltyMPC_dzcc07 = statePenaltyMPC_dz_cc + 182;
statePenaltyMPC_FLOAT* statePenaltyMPC_rd07 = statePenaltyMPC_rd + 182;
statePenaltyMPC_FLOAT statePenaltyMPC_Lbyrd07[26];
statePenaltyMPC_FLOAT* statePenaltyMPC_grad_cost07 = statePenaltyMPC_grad_cost + 182;
statePenaltyMPC_FLOAT* statePenaltyMPC_grad_eq07 = statePenaltyMPC_grad_eq + 182;
statePenaltyMPC_FLOAT* statePenaltyMPC_grad_ineq07 = statePenaltyMPC_grad_ineq + 182;
statePenaltyMPC_FLOAT statePenaltyMPC_ctv07[26];
statePenaltyMPC_FLOAT* statePenaltyMPC_v07 = statePenaltyMPC_v + 56;
statePenaltyMPC_FLOAT statePenaltyMPC_re07[8];
statePenaltyMPC_FLOAT statePenaltyMPC_beta07[8];
statePenaltyMPC_FLOAT statePenaltyMPC_betacc07[8];
statePenaltyMPC_FLOAT* statePenaltyMPC_dvaff07 = statePenaltyMPC_dv_aff + 56;
statePenaltyMPC_FLOAT* statePenaltyMPC_dvcc07 = statePenaltyMPC_dv_cc + 56;
statePenaltyMPC_FLOAT statePenaltyMPC_V07[208];
statePenaltyMPC_FLOAT statePenaltyMPC_Yd07[36];
statePenaltyMPC_FLOAT statePenaltyMPC_Ld07[36];
statePenaltyMPC_FLOAT statePenaltyMPC_yy07[8];
statePenaltyMPC_FLOAT statePenaltyMPC_bmy07[8];
int statePenaltyMPC_lbIdx07[26] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25};
statePenaltyMPC_FLOAT* statePenaltyMPC_llb07 = statePenaltyMPC_l + 252;
statePenaltyMPC_FLOAT* statePenaltyMPC_slb07 = statePenaltyMPC_s + 252;
statePenaltyMPC_FLOAT* statePenaltyMPC_llbbyslb07 = statePenaltyMPC_lbys + 252;
statePenaltyMPC_FLOAT statePenaltyMPC_rilb07[26];
statePenaltyMPC_FLOAT* statePenaltyMPC_dllbaff07 = statePenaltyMPC_dl_aff + 252;
statePenaltyMPC_FLOAT* statePenaltyMPC_dslbaff07 = statePenaltyMPC_ds_aff + 252;
statePenaltyMPC_FLOAT* statePenaltyMPC_dllbcc07 = statePenaltyMPC_dl_cc + 252;
statePenaltyMPC_FLOAT* statePenaltyMPC_dslbcc07 = statePenaltyMPC_ds_cc + 252;
statePenaltyMPC_FLOAT* statePenaltyMPC_ccrhsl07 = statePenaltyMPC_ccrhs + 252;
int statePenaltyMPC_ubIdx07[10] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
statePenaltyMPC_FLOAT* statePenaltyMPC_lub07 = statePenaltyMPC_l + 278;
statePenaltyMPC_FLOAT* statePenaltyMPC_sub07 = statePenaltyMPC_s + 278;
statePenaltyMPC_FLOAT* statePenaltyMPC_lubbysub07 = statePenaltyMPC_lbys + 278;
statePenaltyMPC_FLOAT statePenaltyMPC_riub07[10];
statePenaltyMPC_FLOAT* statePenaltyMPC_dlubaff07 = statePenaltyMPC_dl_aff + 278;
statePenaltyMPC_FLOAT* statePenaltyMPC_dsubaff07 = statePenaltyMPC_ds_aff + 278;
statePenaltyMPC_FLOAT* statePenaltyMPC_dlubcc07 = statePenaltyMPC_dl_cc + 278;
statePenaltyMPC_FLOAT* statePenaltyMPC_dsubcc07 = statePenaltyMPC_ds_cc + 278;
statePenaltyMPC_FLOAT* statePenaltyMPC_ccrhsub07 = statePenaltyMPC_ccrhs + 278;
statePenaltyMPC_FLOAT statePenaltyMPC_Phi07[26];
statePenaltyMPC_FLOAT statePenaltyMPC_W07[26];
statePenaltyMPC_FLOAT statePenaltyMPC_Ysd07[64];
statePenaltyMPC_FLOAT statePenaltyMPC_Lsd07[64];
statePenaltyMPC_FLOAT* statePenaltyMPC_z08 = statePenaltyMPC_z + 208;
statePenaltyMPC_FLOAT* statePenaltyMPC_dzaff08 = statePenaltyMPC_dz_aff + 208;
statePenaltyMPC_FLOAT* statePenaltyMPC_dzcc08 = statePenaltyMPC_dz_cc + 208;
statePenaltyMPC_FLOAT* statePenaltyMPC_rd08 = statePenaltyMPC_rd + 208;
statePenaltyMPC_FLOAT statePenaltyMPC_Lbyrd08[26];
statePenaltyMPC_FLOAT* statePenaltyMPC_grad_cost08 = statePenaltyMPC_grad_cost + 208;
statePenaltyMPC_FLOAT* statePenaltyMPC_grad_eq08 = statePenaltyMPC_grad_eq + 208;
statePenaltyMPC_FLOAT* statePenaltyMPC_grad_ineq08 = statePenaltyMPC_grad_ineq + 208;
statePenaltyMPC_FLOAT statePenaltyMPC_ctv08[26];
statePenaltyMPC_FLOAT* statePenaltyMPC_v08 = statePenaltyMPC_v + 64;
statePenaltyMPC_FLOAT statePenaltyMPC_re08[8];
statePenaltyMPC_FLOAT statePenaltyMPC_beta08[8];
statePenaltyMPC_FLOAT statePenaltyMPC_betacc08[8];
statePenaltyMPC_FLOAT* statePenaltyMPC_dvaff08 = statePenaltyMPC_dv_aff + 64;
statePenaltyMPC_FLOAT* statePenaltyMPC_dvcc08 = statePenaltyMPC_dv_cc + 64;
statePenaltyMPC_FLOAT statePenaltyMPC_V08[208];
statePenaltyMPC_FLOAT statePenaltyMPC_Yd08[36];
statePenaltyMPC_FLOAT statePenaltyMPC_Ld08[36];
statePenaltyMPC_FLOAT statePenaltyMPC_yy08[8];
statePenaltyMPC_FLOAT statePenaltyMPC_bmy08[8];
int statePenaltyMPC_lbIdx08[26] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25};
statePenaltyMPC_FLOAT* statePenaltyMPC_llb08 = statePenaltyMPC_l + 288;
statePenaltyMPC_FLOAT* statePenaltyMPC_slb08 = statePenaltyMPC_s + 288;
statePenaltyMPC_FLOAT* statePenaltyMPC_llbbyslb08 = statePenaltyMPC_lbys + 288;
statePenaltyMPC_FLOAT statePenaltyMPC_rilb08[26];
statePenaltyMPC_FLOAT* statePenaltyMPC_dllbaff08 = statePenaltyMPC_dl_aff + 288;
statePenaltyMPC_FLOAT* statePenaltyMPC_dslbaff08 = statePenaltyMPC_ds_aff + 288;
statePenaltyMPC_FLOAT* statePenaltyMPC_dllbcc08 = statePenaltyMPC_dl_cc + 288;
statePenaltyMPC_FLOAT* statePenaltyMPC_dslbcc08 = statePenaltyMPC_ds_cc + 288;
statePenaltyMPC_FLOAT* statePenaltyMPC_ccrhsl08 = statePenaltyMPC_ccrhs + 288;
int statePenaltyMPC_ubIdx08[10] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
statePenaltyMPC_FLOAT* statePenaltyMPC_lub08 = statePenaltyMPC_l + 314;
statePenaltyMPC_FLOAT* statePenaltyMPC_sub08 = statePenaltyMPC_s + 314;
statePenaltyMPC_FLOAT* statePenaltyMPC_lubbysub08 = statePenaltyMPC_lbys + 314;
statePenaltyMPC_FLOAT statePenaltyMPC_riub08[10];
statePenaltyMPC_FLOAT* statePenaltyMPC_dlubaff08 = statePenaltyMPC_dl_aff + 314;
statePenaltyMPC_FLOAT* statePenaltyMPC_dsubaff08 = statePenaltyMPC_ds_aff + 314;
statePenaltyMPC_FLOAT* statePenaltyMPC_dlubcc08 = statePenaltyMPC_dl_cc + 314;
statePenaltyMPC_FLOAT* statePenaltyMPC_dsubcc08 = statePenaltyMPC_ds_cc + 314;
statePenaltyMPC_FLOAT* statePenaltyMPC_ccrhsub08 = statePenaltyMPC_ccrhs + 314;
statePenaltyMPC_FLOAT statePenaltyMPC_Phi08[26];
statePenaltyMPC_FLOAT statePenaltyMPC_W08[26];
statePenaltyMPC_FLOAT statePenaltyMPC_Ysd08[64];
statePenaltyMPC_FLOAT statePenaltyMPC_Lsd08[64];
statePenaltyMPC_FLOAT* statePenaltyMPC_z09 = statePenaltyMPC_z + 234;
statePenaltyMPC_FLOAT* statePenaltyMPC_dzaff09 = statePenaltyMPC_dz_aff + 234;
statePenaltyMPC_FLOAT* statePenaltyMPC_dzcc09 = statePenaltyMPC_dz_cc + 234;
statePenaltyMPC_FLOAT* statePenaltyMPC_rd09 = statePenaltyMPC_rd + 234;
statePenaltyMPC_FLOAT statePenaltyMPC_Lbyrd09[26];
statePenaltyMPC_FLOAT* statePenaltyMPC_grad_cost09 = statePenaltyMPC_grad_cost + 234;
statePenaltyMPC_FLOAT* statePenaltyMPC_grad_eq09 = statePenaltyMPC_grad_eq + 234;
statePenaltyMPC_FLOAT* statePenaltyMPC_grad_ineq09 = statePenaltyMPC_grad_ineq + 234;
statePenaltyMPC_FLOAT statePenaltyMPC_ctv09[26];
statePenaltyMPC_FLOAT* statePenaltyMPC_v09 = statePenaltyMPC_v + 72;
statePenaltyMPC_FLOAT statePenaltyMPC_re09[8];
statePenaltyMPC_FLOAT statePenaltyMPC_beta09[8];
statePenaltyMPC_FLOAT statePenaltyMPC_betacc09[8];
statePenaltyMPC_FLOAT* statePenaltyMPC_dvaff09 = statePenaltyMPC_dv_aff + 72;
statePenaltyMPC_FLOAT* statePenaltyMPC_dvcc09 = statePenaltyMPC_dv_cc + 72;
statePenaltyMPC_FLOAT statePenaltyMPC_V09[208];
statePenaltyMPC_FLOAT statePenaltyMPC_Yd09[36];
statePenaltyMPC_FLOAT statePenaltyMPC_Ld09[36];
statePenaltyMPC_FLOAT statePenaltyMPC_yy09[8];
statePenaltyMPC_FLOAT statePenaltyMPC_bmy09[8];
int statePenaltyMPC_lbIdx09[26] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25};
statePenaltyMPC_FLOAT* statePenaltyMPC_llb09 = statePenaltyMPC_l + 324;
statePenaltyMPC_FLOAT* statePenaltyMPC_slb09 = statePenaltyMPC_s + 324;
statePenaltyMPC_FLOAT* statePenaltyMPC_llbbyslb09 = statePenaltyMPC_lbys + 324;
statePenaltyMPC_FLOAT statePenaltyMPC_rilb09[26];
statePenaltyMPC_FLOAT* statePenaltyMPC_dllbaff09 = statePenaltyMPC_dl_aff + 324;
statePenaltyMPC_FLOAT* statePenaltyMPC_dslbaff09 = statePenaltyMPC_ds_aff + 324;
statePenaltyMPC_FLOAT* statePenaltyMPC_dllbcc09 = statePenaltyMPC_dl_cc + 324;
statePenaltyMPC_FLOAT* statePenaltyMPC_dslbcc09 = statePenaltyMPC_ds_cc + 324;
statePenaltyMPC_FLOAT* statePenaltyMPC_ccrhsl09 = statePenaltyMPC_ccrhs + 324;
int statePenaltyMPC_ubIdx09[10] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
statePenaltyMPC_FLOAT* statePenaltyMPC_lub09 = statePenaltyMPC_l + 350;
statePenaltyMPC_FLOAT* statePenaltyMPC_sub09 = statePenaltyMPC_s + 350;
statePenaltyMPC_FLOAT* statePenaltyMPC_lubbysub09 = statePenaltyMPC_lbys + 350;
statePenaltyMPC_FLOAT statePenaltyMPC_riub09[10];
statePenaltyMPC_FLOAT* statePenaltyMPC_dlubaff09 = statePenaltyMPC_dl_aff + 350;
statePenaltyMPC_FLOAT* statePenaltyMPC_dsubaff09 = statePenaltyMPC_ds_aff + 350;
statePenaltyMPC_FLOAT* statePenaltyMPC_dlubcc09 = statePenaltyMPC_dl_cc + 350;
statePenaltyMPC_FLOAT* statePenaltyMPC_dsubcc09 = statePenaltyMPC_ds_cc + 350;
statePenaltyMPC_FLOAT* statePenaltyMPC_ccrhsub09 = statePenaltyMPC_ccrhs + 350;
statePenaltyMPC_FLOAT statePenaltyMPC_Phi09[26];
statePenaltyMPC_FLOAT statePenaltyMPC_W09[26];
statePenaltyMPC_FLOAT statePenaltyMPC_Ysd09[64];
statePenaltyMPC_FLOAT statePenaltyMPC_Lsd09[64];
statePenaltyMPC_FLOAT* statePenaltyMPC_z10 = statePenaltyMPC_z + 260;
statePenaltyMPC_FLOAT* statePenaltyMPC_dzaff10 = statePenaltyMPC_dz_aff + 260;
statePenaltyMPC_FLOAT* statePenaltyMPC_dzcc10 = statePenaltyMPC_dz_cc + 260;
statePenaltyMPC_FLOAT* statePenaltyMPC_rd10 = statePenaltyMPC_rd + 260;
statePenaltyMPC_FLOAT statePenaltyMPC_Lbyrd10[26];
statePenaltyMPC_FLOAT* statePenaltyMPC_grad_cost10 = statePenaltyMPC_grad_cost + 260;
statePenaltyMPC_FLOAT* statePenaltyMPC_grad_eq10 = statePenaltyMPC_grad_eq + 260;
statePenaltyMPC_FLOAT* statePenaltyMPC_grad_ineq10 = statePenaltyMPC_grad_ineq + 260;
statePenaltyMPC_FLOAT statePenaltyMPC_ctv10[26];
statePenaltyMPC_FLOAT* statePenaltyMPC_v10 = statePenaltyMPC_v + 80;
statePenaltyMPC_FLOAT statePenaltyMPC_re10[8];
statePenaltyMPC_FLOAT statePenaltyMPC_beta10[8];
statePenaltyMPC_FLOAT statePenaltyMPC_betacc10[8];
statePenaltyMPC_FLOAT* statePenaltyMPC_dvaff10 = statePenaltyMPC_dv_aff + 80;
statePenaltyMPC_FLOAT* statePenaltyMPC_dvcc10 = statePenaltyMPC_dv_cc + 80;
statePenaltyMPC_FLOAT statePenaltyMPC_V10[208];
statePenaltyMPC_FLOAT statePenaltyMPC_Yd10[36];
statePenaltyMPC_FLOAT statePenaltyMPC_Ld10[36];
statePenaltyMPC_FLOAT statePenaltyMPC_yy10[8];
statePenaltyMPC_FLOAT statePenaltyMPC_bmy10[8];
int statePenaltyMPC_lbIdx10[26] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25};
statePenaltyMPC_FLOAT* statePenaltyMPC_llb10 = statePenaltyMPC_l + 360;
statePenaltyMPC_FLOAT* statePenaltyMPC_slb10 = statePenaltyMPC_s + 360;
statePenaltyMPC_FLOAT* statePenaltyMPC_llbbyslb10 = statePenaltyMPC_lbys + 360;
statePenaltyMPC_FLOAT statePenaltyMPC_rilb10[26];
statePenaltyMPC_FLOAT* statePenaltyMPC_dllbaff10 = statePenaltyMPC_dl_aff + 360;
statePenaltyMPC_FLOAT* statePenaltyMPC_dslbaff10 = statePenaltyMPC_ds_aff + 360;
statePenaltyMPC_FLOAT* statePenaltyMPC_dllbcc10 = statePenaltyMPC_dl_cc + 360;
statePenaltyMPC_FLOAT* statePenaltyMPC_dslbcc10 = statePenaltyMPC_ds_cc + 360;
statePenaltyMPC_FLOAT* statePenaltyMPC_ccrhsl10 = statePenaltyMPC_ccrhs + 360;
int statePenaltyMPC_ubIdx10[10] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
statePenaltyMPC_FLOAT* statePenaltyMPC_lub10 = statePenaltyMPC_l + 386;
statePenaltyMPC_FLOAT* statePenaltyMPC_sub10 = statePenaltyMPC_s + 386;
statePenaltyMPC_FLOAT* statePenaltyMPC_lubbysub10 = statePenaltyMPC_lbys + 386;
statePenaltyMPC_FLOAT statePenaltyMPC_riub10[10];
statePenaltyMPC_FLOAT* statePenaltyMPC_dlubaff10 = statePenaltyMPC_dl_aff + 386;
statePenaltyMPC_FLOAT* statePenaltyMPC_dsubaff10 = statePenaltyMPC_ds_aff + 386;
statePenaltyMPC_FLOAT* statePenaltyMPC_dlubcc10 = statePenaltyMPC_dl_cc + 386;
statePenaltyMPC_FLOAT* statePenaltyMPC_dsubcc10 = statePenaltyMPC_ds_cc + 386;
statePenaltyMPC_FLOAT* statePenaltyMPC_ccrhsub10 = statePenaltyMPC_ccrhs + 386;
statePenaltyMPC_FLOAT statePenaltyMPC_Phi10[26];
statePenaltyMPC_FLOAT statePenaltyMPC_W10[26];
statePenaltyMPC_FLOAT statePenaltyMPC_Ysd10[64];
statePenaltyMPC_FLOAT statePenaltyMPC_Lsd10[64];
statePenaltyMPC_FLOAT* statePenaltyMPC_z11 = statePenaltyMPC_z + 286;
statePenaltyMPC_FLOAT* statePenaltyMPC_dzaff11 = statePenaltyMPC_dz_aff + 286;
statePenaltyMPC_FLOAT* statePenaltyMPC_dzcc11 = statePenaltyMPC_dz_cc + 286;
statePenaltyMPC_FLOAT* statePenaltyMPC_rd11 = statePenaltyMPC_rd + 286;
statePenaltyMPC_FLOAT statePenaltyMPC_Lbyrd11[26];
statePenaltyMPC_FLOAT* statePenaltyMPC_grad_cost11 = statePenaltyMPC_grad_cost + 286;
statePenaltyMPC_FLOAT* statePenaltyMPC_grad_eq11 = statePenaltyMPC_grad_eq + 286;
statePenaltyMPC_FLOAT* statePenaltyMPC_grad_ineq11 = statePenaltyMPC_grad_ineq + 286;
statePenaltyMPC_FLOAT statePenaltyMPC_ctv11[26];
statePenaltyMPC_FLOAT* statePenaltyMPC_v11 = statePenaltyMPC_v + 88;
statePenaltyMPC_FLOAT statePenaltyMPC_re11[8];
statePenaltyMPC_FLOAT statePenaltyMPC_beta11[8];
statePenaltyMPC_FLOAT statePenaltyMPC_betacc11[8];
statePenaltyMPC_FLOAT* statePenaltyMPC_dvaff11 = statePenaltyMPC_dv_aff + 88;
statePenaltyMPC_FLOAT* statePenaltyMPC_dvcc11 = statePenaltyMPC_dv_cc + 88;
statePenaltyMPC_FLOAT statePenaltyMPC_V11[208];
statePenaltyMPC_FLOAT statePenaltyMPC_Yd11[36];
statePenaltyMPC_FLOAT statePenaltyMPC_Ld11[36];
statePenaltyMPC_FLOAT statePenaltyMPC_yy11[8];
statePenaltyMPC_FLOAT statePenaltyMPC_bmy11[8];
int statePenaltyMPC_lbIdx11[26] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25};
statePenaltyMPC_FLOAT* statePenaltyMPC_llb11 = statePenaltyMPC_l + 396;
statePenaltyMPC_FLOAT* statePenaltyMPC_slb11 = statePenaltyMPC_s + 396;
statePenaltyMPC_FLOAT* statePenaltyMPC_llbbyslb11 = statePenaltyMPC_lbys + 396;
statePenaltyMPC_FLOAT statePenaltyMPC_rilb11[26];
statePenaltyMPC_FLOAT* statePenaltyMPC_dllbaff11 = statePenaltyMPC_dl_aff + 396;
statePenaltyMPC_FLOAT* statePenaltyMPC_dslbaff11 = statePenaltyMPC_ds_aff + 396;
statePenaltyMPC_FLOAT* statePenaltyMPC_dllbcc11 = statePenaltyMPC_dl_cc + 396;
statePenaltyMPC_FLOAT* statePenaltyMPC_dslbcc11 = statePenaltyMPC_ds_cc + 396;
statePenaltyMPC_FLOAT* statePenaltyMPC_ccrhsl11 = statePenaltyMPC_ccrhs + 396;
int statePenaltyMPC_ubIdx11[10] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
statePenaltyMPC_FLOAT* statePenaltyMPC_lub11 = statePenaltyMPC_l + 422;
statePenaltyMPC_FLOAT* statePenaltyMPC_sub11 = statePenaltyMPC_s + 422;
statePenaltyMPC_FLOAT* statePenaltyMPC_lubbysub11 = statePenaltyMPC_lbys + 422;
statePenaltyMPC_FLOAT statePenaltyMPC_riub11[10];
statePenaltyMPC_FLOAT* statePenaltyMPC_dlubaff11 = statePenaltyMPC_dl_aff + 422;
statePenaltyMPC_FLOAT* statePenaltyMPC_dsubaff11 = statePenaltyMPC_ds_aff + 422;
statePenaltyMPC_FLOAT* statePenaltyMPC_dlubcc11 = statePenaltyMPC_dl_cc + 422;
statePenaltyMPC_FLOAT* statePenaltyMPC_dsubcc11 = statePenaltyMPC_ds_cc + 422;
statePenaltyMPC_FLOAT* statePenaltyMPC_ccrhsub11 = statePenaltyMPC_ccrhs + 422;
statePenaltyMPC_FLOAT statePenaltyMPC_Phi11[26];
statePenaltyMPC_FLOAT statePenaltyMPC_W11[26];
statePenaltyMPC_FLOAT statePenaltyMPC_Ysd11[64];
statePenaltyMPC_FLOAT statePenaltyMPC_Lsd11[64];
statePenaltyMPC_FLOAT* statePenaltyMPC_z12 = statePenaltyMPC_z + 312;
statePenaltyMPC_FLOAT* statePenaltyMPC_dzaff12 = statePenaltyMPC_dz_aff + 312;
statePenaltyMPC_FLOAT* statePenaltyMPC_dzcc12 = statePenaltyMPC_dz_cc + 312;
statePenaltyMPC_FLOAT* statePenaltyMPC_rd12 = statePenaltyMPC_rd + 312;
statePenaltyMPC_FLOAT statePenaltyMPC_Lbyrd12[26];
statePenaltyMPC_FLOAT* statePenaltyMPC_grad_cost12 = statePenaltyMPC_grad_cost + 312;
statePenaltyMPC_FLOAT* statePenaltyMPC_grad_eq12 = statePenaltyMPC_grad_eq + 312;
statePenaltyMPC_FLOAT* statePenaltyMPC_grad_ineq12 = statePenaltyMPC_grad_ineq + 312;
statePenaltyMPC_FLOAT statePenaltyMPC_ctv12[26];
statePenaltyMPC_FLOAT* statePenaltyMPC_v12 = statePenaltyMPC_v + 96;
statePenaltyMPC_FLOAT statePenaltyMPC_re12[8];
statePenaltyMPC_FLOAT statePenaltyMPC_beta12[8];
statePenaltyMPC_FLOAT statePenaltyMPC_betacc12[8];
statePenaltyMPC_FLOAT* statePenaltyMPC_dvaff12 = statePenaltyMPC_dv_aff + 96;
statePenaltyMPC_FLOAT* statePenaltyMPC_dvcc12 = statePenaltyMPC_dv_cc + 96;
statePenaltyMPC_FLOAT statePenaltyMPC_V12[208];
statePenaltyMPC_FLOAT statePenaltyMPC_Yd12[36];
statePenaltyMPC_FLOAT statePenaltyMPC_Ld12[36];
statePenaltyMPC_FLOAT statePenaltyMPC_yy12[8];
statePenaltyMPC_FLOAT statePenaltyMPC_bmy12[8];
int statePenaltyMPC_lbIdx12[26] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25};
statePenaltyMPC_FLOAT* statePenaltyMPC_llb12 = statePenaltyMPC_l + 432;
statePenaltyMPC_FLOAT* statePenaltyMPC_slb12 = statePenaltyMPC_s + 432;
statePenaltyMPC_FLOAT* statePenaltyMPC_llbbyslb12 = statePenaltyMPC_lbys + 432;
statePenaltyMPC_FLOAT statePenaltyMPC_rilb12[26];
statePenaltyMPC_FLOAT* statePenaltyMPC_dllbaff12 = statePenaltyMPC_dl_aff + 432;
statePenaltyMPC_FLOAT* statePenaltyMPC_dslbaff12 = statePenaltyMPC_ds_aff + 432;
statePenaltyMPC_FLOAT* statePenaltyMPC_dllbcc12 = statePenaltyMPC_dl_cc + 432;
statePenaltyMPC_FLOAT* statePenaltyMPC_dslbcc12 = statePenaltyMPC_ds_cc + 432;
statePenaltyMPC_FLOAT* statePenaltyMPC_ccrhsl12 = statePenaltyMPC_ccrhs + 432;
int statePenaltyMPC_ubIdx12[10] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
statePenaltyMPC_FLOAT* statePenaltyMPC_lub12 = statePenaltyMPC_l + 458;
statePenaltyMPC_FLOAT* statePenaltyMPC_sub12 = statePenaltyMPC_s + 458;
statePenaltyMPC_FLOAT* statePenaltyMPC_lubbysub12 = statePenaltyMPC_lbys + 458;
statePenaltyMPC_FLOAT statePenaltyMPC_riub12[10];
statePenaltyMPC_FLOAT* statePenaltyMPC_dlubaff12 = statePenaltyMPC_dl_aff + 458;
statePenaltyMPC_FLOAT* statePenaltyMPC_dsubaff12 = statePenaltyMPC_ds_aff + 458;
statePenaltyMPC_FLOAT* statePenaltyMPC_dlubcc12 = statePenaltyMPC_dl_cc + 458;
statePenaltyMPC_FLOAT* statePenaltyMPC_dsubcc12 = statePenaltyMPC_ds_cc + 458;
statePenaltyMPC_FLOAT* statePenaltyMPC_ccrhsub12 = statePenaltyMPC_ccrhs + 458;
statePenaltyMPC_FLOAT statePenaltyMPC_Phi12[26];
statePenaltyMPC_FLOAT statePenaltyMPC_W12[26];
statePenaltyMPC_FLOAT statePenaltyMPC_Ysd12[64];
statePenaltyMPC_FLOAT statePenaltyMPC_Lsd12[64];
statePenaltyMPC_FLOAT* statePenaltyMPC_z13 = statePenaltyMPC_z + 338;
statePenaltyMPC_FLOAT* statePenaltyMPC_dzaff13 = statePenaltyMPC_dz_aff + 338;
statePenaltyMPC_FLOAT* statePenaltyMPC_dzcc13 = statePenaltyMPC_dz_cc + 338;
statePenaltyMPC_FLOAT* statePenaltyMPC_rd13 = statePenaltyMPC_rd + 338;
statePenaltyMPC_FLOAT statePenaltyMPC_Lbyrd13[26];
statePenaltyMPC_FLOAT* statePenaltyMPC_grad_cost13 = statePenaltyMPC_grad_cost + 338;
statePenaltyMPC_FLOAT* statePenaltyMPC_grad_eq13 = statePenaltyMPC_grad_eq + 338;
statePenaltyMPC_FLOAT* statePenaltyMPC_grad_ineq13 = statePenaltyMPC_grad_ineq + 338;
statePenaltyMPC_FLOAT statePenaltyMPC_ctv13[26];
statePenaltyMPC_FLOAT* statePenaltyMPC_v13 = statePenaltyMPC_v + 104;
statePenaltyMPC_FLOAT statePenaltyMPC_re13[8];
statePenaltyMPC_FLOAT statePenaltyMPC_beta13[8];
statePenaltyMPC_FLOAT statePenaltyMPC_betacc13[8];
statePenaltyMPC_FLOAT* statePenaltyMPC_dvaff13 = statePenaltyMPC_dv_aff + 104;
statePenaltyMPC_FLOAT* statePenaltyMPC_dvcc13 = statePenaltyMPC_dv_cc + 104;
statePenaltyMPC_FLOAT statePenaltyMPC_V13[208];
statePenaltyMPC_FLOAT statePenaltyMPC_Yd13[36];
statePenaltyMPC_FLOAT statePenaltyMPC_Ld13[36];
statePenaltyMPC_FLOAT statePenaltyMPC_yy13[8];
statePenaltyMPC_FLOAT statePenaltyMPC_bmy13[8];
int statePenaltyMPC_lbIdx13[26] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25};
statePenaltyMPC_FLOAT* statePenaltyMPC_llb13 = statePenaltyMPC_l + 468;
statePenaltyMPC_FLOAT* statePenaltyMPC_slb13 = statePenaltyMPC_s + 468;
statePenaltyMPC_FLOAT* statePenaltyMPC_llbbyslb13 = statePenaltyMPC_lbys + 468;
statePenaltyMPC_FLOAT statePenaltyMPC_rilb13[26];
statePenaltyMPC_FLOAT* statePenaltyMPC_dllbaff13 = statePenaltyMPC_dl_aff + 468;
statePenaltyMPC_FLOAT* statePenaltyMPC_dslbaff13 = statePenaltyMPC_ds_aff + 468;
statePenaltyMPC_FLOAT* statePenaltyMPC_dllbcc13 = statePenaltyMPC_dl_cc + 468;
statePenaltyMPC_FLOAT* statePenaltyMPC_dslbcc13 = statePenaltyMPC_ds_cc + 468;
statePenaltyMPC_FLOAT* statePenaltyMPC_ccrhsl13 = statePenaltyMPC_ccrhs + 468;
int statePenaltyMPC_ubIdx13[10] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
statePenaltyMPC_FLOAT* statePenaltyMPC_lub13 = statePenaltyMPC_l + 494;
statePenaltyMPC_FLOAT* statePenaltyMPC_sub13 = statePenaltyMPC_s + 494;
statePenaltyMPC_FLOAT* statePenaltyMPC_lubbysub13 = statePenaltyMPC_lbys + 494;
statePenaltyMPC_FLOAT statePenaltyMPC_riub13[10];
statePenaltyMPC_FLOAT* statePenaltyMPC_dlubaff13 = statePenaltyMPC_dl_aff + 494;
statePenaltyMPC_FLOAT* statePenaltyMPC_dsubaff13 = statePenaltyMPC_ds_aff + 494;
statePenaltyMPC_FLOAT* statePenaltyMPC_dlubcc13 = statePenaltyMPC_dl_cc + 494;
statePenaltyMPC_FLOAT* statePenaltyMPC_dsubcc13 = statePenaltyMPC_ds_cc + 494;
statePenaltyMPC_FLOAT* statePenaltyMPC_ccrhsub13 = statePenaltyMPC_ccrhs + 494;
statePenaltyMPC_FLOAT statePenaltyMPC_Phi13[26];
statePenaltyMPC_FLOAT statePenaltyMPC_W13[26];
statePenaltyMPC_FLOAT statePenaltyMPC_Ysd13[64];
statePenaltyMPC_FLOAT statePenaltyMPC_Lsd13[64];
statePenaltyMPC_FLOAT statePenaltyMPC_f14[8] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
statePenaltyMPC_FLOAT* statePenaltyMPC_z14 = statePenaltyMPC_z + 364;
statePenaltyMPC_FLOAT* statePenaltyMPC_dzaff14 = statePenaltyMPC_dz_aff + 364;
statePenaltyMPC_FLOAT* statePenaltyMPC_dzcc14 = statePenaltyMPC_dz_cc + 364;
statePenaltyMPC_FLOAT* statePenaltyMPC_rd14 = statePenaltyMPC_rd + 364;
statePenaltyMPC_FLOAT statePenaltyMPC_Lbyrd14[8];
statePenaltyMPC_FLOAT* statePenaltyMPC_grad_cost14 = statePenaltyMPC_grad_cost + 364;
statePenaltyMPC_FLOAT* statePenaltyMPC_grad_eq14 = statePenaltyMPC_grad_eq + 364;
statePenaltyMPC_FLOAT* statePenaltyMPC_grad_ineq14 = statePenaltyMPC_grad_ineq + 364;
statePenaltyMPC_FLOAT statePenaltyMPC_ctv14[8];
statePenaltyMPC_FLOAT* statePenaltyMPC_v14 = statePenaltyMPC_v + 112;
statePenaltyMPC_FLOAT statePenaltyMPC_re14[8];
statePenaltyMPC_FLOAT statePenaltyMPC_beta14[8];
statePenaltyMPC_FLOAT statePenaltyMPC_betacc14[8];
statePenaltyMPC_FLOAT* statePenaltyMPC_dvaff14 = statePenaltyMPC_dv_aff + 112;
statePenaltyMPC_FLOAT* statePenaltyMPC_dvcc14 = statePenaltyMPC_dv_cc + 112;
statePenaltyMPC_FLOAT statePenaltyMPC_V14[64];
statePenaltyMPC_FLOAT statePenaltyMPC_Yd14[36];
statePenaltyMPC_FLOAT statePenaltyMPC_Ld14[36];
statePenaltyMPC_FLOAT statePenaltyMPC_yy14[8];
statePenaltyMPC_FLOAT statePenaltyMPC_bmy14[8];
int statePenaltyMPC_lbIdx14[8] = {0, 1, 2, 3, 4, 5, 6, 7};
statePenaltyMPC_FLOAT* statePenaltyMPC_llb14 = statePenaltyMPC_l + 504;
statePenaltyMPC_FLOAT* statePenaltyMPC_slb14 = statePenaltyMPC_s + 504;
statePenaltyMPC_FLOAT* statePenaltyMPC_llbbyslb14 = statePenaltyMPC_lbys + 504;
statePenaltyMPC_FLOAT statePenaltyMPC_rilb14[8];
statePenaltyMPC_FLOAT* statePenaltyMPC_dllbaff14 = statePenaltyMPC_dl_aff + 504;
statePenaltyMPC_FLOAT* statePenaltyMPC_dslbaff14 = statePenaltyMPC_ds_aff + 504;
statePenaltyMPC_FLOAT* statePenaltyMPC_dllbcc14 = statePenaltyMPC_dl_cc + 504;
statePenaltyMPC_FLOAT* statePenaltyMPC_dslbcc14 = statePenaltyMPC_ds_cc + 504;
statePenaltyMPC_FLOAT* statePenaltyMPC_ccrhsl14 = statePenaltyMPC_ccrhs + 504;
int statePenaltyMPC_ubIdx14[8] = {0, 1, 2, 3, 4, 5, 6, 7};
statePenaltyMPC_FLOAT* statePenaltyMPC_lub14 = statePenaltyMPC_l + 512;
statePenaltyMPC_FLOAT* statePenaltyMPC_sub14 = statePenaltyMPC_s + 512;
statePenaltyMPC_FLOAT* statePenaltyMPC_lubbysub14 = statePenaltyMPC_lbys + 512;
statePenaltyMPC_FLOAT statePenaltyMPC_riub14[8];
statePenaltyMPC_FLOAT* statePenaltyMPC_dlubaff14 = statePenaltyMPC_dl_aff + 512;
statePenaltyMPC_FLOAT* statePenaltyMPC_dsubaff14 = statePenaltyMPC_ds_aff + 512;
statePenaltyMPC_FLOAT* statePenaltyMPC_dlubcc14 = statePenaltyMPC_dl_cc + 512;
statePenaltyMPC_FLOAT* statePenaltyMPC_dsubcc14 = statePenaltyMPC_ds_cc + 512;
statePenaltyMPC_FLOAT* statePenaltyMPC_ccrhsub14 = statePenaltyMPC_ccrhs + 512;
statePenaltyMPC_FLOAT statePenaltyMPC_Phi14[8];
statePenaltyMPC_FLOAT statePenaltyMPC_D14[8] = {-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000};
statePenaltyMPC_FLOAT statePenaltyMPC_W14[8];
statePenaltyMPC_FLOAT statePenaltyMPC_Ysd14[64];
statePenaltyMPC_FLOAT statePenaltyMPC_Lsd14[64];
statePenaltyMPC_FLOAT musigma;
statePenaltyMPC_FLOAT sigma_3rdroot;
statePenaltyMPC_FLOAT statePenaltyMPC_Diag1_0[26];
statePenaltyMPC_FLOAT statePenaltyMPC_Diag2_0[26];
statePenaltyMPC_FLOAT statePenaltyMPC_L_0[325];




/* SOLVER CODE --------------------------------------------------------- */
int statePenaltyMPC_solve(statePenaltyMPC_params* params, statePenaltyMPC_output* output, statePenaltyMPC_info* info)
{	
int exitcode;

#if statePenaltyMPC_SET_TIMING == 1
	statePenaltyMPC_timer solvertimer;
	statePenaltyMPC_tic(&solvertimer);
#endif
/* FUNCTION CALLS INTO LA LIBRARY -------------------------------------- */
info->it = 0;
statePenaltyMPC_LA_INITIALIZEVECTOR_372(statePenaltyMPC_z, 0);
statePenaltyMPC_LA_INITIALIZEVECTOR_120(statePenaltyMPC_v, 1);
statePenaltyMPC_LA_INITIALIZEVECTOR_520(statePenaltyMPC_l, 1);
statePenaltyMPC_LA_INITIALIZEVECTOR_520(statePenaltyMPC_s, 1);
info->mu = 0;
statePenaltyMPC_LA_DOTACC_520(statePenaltyMPC_l, statePenaltyMPC_s, &info->mu);
info->mu /= 520;
PRINTTEXT("This is statePenaltyMPC, a solver generated by FORCES (forces.ethz.ch).\n");
PRINTTEXT("(c) Alexander Domahidi, Automatic Control Laboratory, ETH Zurich, 2011-2014.\n");
PRINTTEXT("\n  #it  res_eq   res_ineq     pobj         dobj       dgap     rdgap     mu\n");
PRINTTEXT("  ---------------------------------------------------------------------------\n");
while( 1 ){
info->pobj = 0;
statePenaltyMPC_LA_DIAG_QUADFCN_26(params->H1, params->f1, statePenaltyMPC_z00, statePenaltyMPC_grad_cost00, &info->pobj);
statePenaltyMPC_LA_DIAG_QUADFCN_26(params->H2, params->f2, statePenaltyMPC_z01, statePenaltyMPC_grad_cost01, &info->pobj);
statePenaltyMPC_LA_DIAG_QUADFCN_26(params->H3, params->f3, statePenaltyMPC_z02, statePenaltyMPC_grad_cost02, &info->pobj);
statePenaltyMPC_LA_DIAG_QUADFCN_26(params->H4, params->f4, statePenaltyMPC_z03, statePenaltyMPC_grad_cost03, &info->pobj);
statePenaltyMPC_LA_DIAG_QUADFCN_26(params->H5, params->f5, statePenaltyMPC_z04, statePenaltyMPC_grad_cost04, &info->pobj);
statePenaltyMPC_LA_DIAG_QUADFCN_26(params->H6, params->f6, statePenaltyMPC_z05, statePenaltyMPC_grad_cost05, &info->pobj);
statePenaltyMPC_LA_DIAG_QUADFCN_26(params->H7, params->f7, statePenaltyMPC_z06, statePenaltyMPC_grad_cost06, &info->pobj);
statePenaltyMPC_LA_DIAG_QUADFCN_26(params->H8, params->f8, statePenaltyMPC_z07, statePenaltyMPC_grad_cost07, &info->pobj);
statePenaltyMPC_LA_DIAG_QUADFCN_26(params->H9, params->f9, statePenaltyMPC_z08, statePenaltyMPC_grad_cost08, &info->pobj);
statePenaltyMPC_LA_DIAG_QUADFCN_26(params->H10, params->f10, statePenaltyMPC_z09, statePenaltyMPC_grad_cost09, &info->pobj);
statePenaltyMPC_LA_DIAG_QUADFCN_26(params->H11, params->f11, statePenaltyMPC_z10, statePenaltyMPC_grad_cost10, &info->pobj);
statePenaltyMPC_LA_DIAG_QUADFCN_26(params->H12, params->f12, statePenaltyMPC_z11, statePenaltyMPC_grad_cost11, &info->pobj);
statePenaltyMPC_LA_DIAG_QUADFCN_26(params->H13, params->f13, statePenaltyMPC_z12, statePenaltyMPC_grad_cost12, &info->pobj);
statePenaltyMPC_LA_DIAG_QUADFCN_26(params->H14, params->f14, statePenaltyMPC_z13, statePenaltyMPC_grad_cost13, &info->pobj);
statePenaltyMPC_LA_DIAG_QUADFCN_8(params->H15, statePenaltyMPC_f14, statePenaltyMPC_z14, statePenaltyMPC_grad_cost14, &info->pobj);
info->res_eq = 0;
info->dgap = 0;
statePenaltyMPC_LA_DIAGZERO_MVMSUB6_8(statePenaltyMPC_D00, statePenaltyMPC_z00, params->e1, statePenaltyMPC_v00, statePenaltyMPC_re00, &info->dgap, &info->res_eq);
statePenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_8_26_26(params->C1, statePenaltyMPC_z00, statePenaltyMPC_D01, statePenaltyMPC_z01, params->e2, statePenaltyMPC_v01, statePenaltyMPC_re01, &info->dgap, &info->res_eq);
statePenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_8_26_26(params->C2, statePenaltyMPC_z01, statePenaltyMPC_D01, statePenaltyMPC_z02, params->e3, statePenaltyMPC_v02, statePenaltyMPC_re02, &info->dgap, &info->res_eq);
statePenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_8_26_26(params->C3, statePenaltyMPC_z02, statePenaltyMPC_D01, statePenaltyMPC_z03, params->e4, statePenaltyMPC_v03, statePenaltyMPC_re03, &info->dgap, &info->res_eq);
statePenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_8_26_26(params->C4, statePenaltyMPC_z03, statePenaltyMPC_D01, statePenaltyMPC_z04, params->e5, statePenaltyMPC_v04, statePenaltyMPC_re04, &info->dgap, &info->res_eq);
statePenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_8_26_26(params->C5, statePenaltyMPC_z04, statePenaltyMPC_D01, statePenaltyMPC_z05, params->e6, statePenaltyMPC_v05, statePenaltyMPC_re05, &info->dgap, &info->res_eq);
statePenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_8_26_26(params->C6, statePenaltyMPC_z05, statePenaltyMPC_D01, statePenaltyMPC_z06, params->e7, statePenaltyMPC_v06, statePenaltyMPC_re06, &info->dgap, &info->res_eq);
statePenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_8_26_26(params->C7, statePenaltyMPC_z06, statePenaltyMPC_D01, statePenaltyMPC_z07, params->e8, statePenaltyMPC_v07, statePenaltyMPC_re07, &info->dgap, &info->res_eq);
statePenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_8_26_26(params->C8, statePenaltyMPC_z07, statePenaltyMPC_D01, statePenaltyMPC_z08, params->e9, statePenaltyMPC_v08, statePenaltyMPC_re08, &info->dgap, &info->res_eq);
statePenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_8_26_26(params->C9, statePenaltyMPC_z08, statePenaltyMPC_D01, statePenaltyMPC_z09, params->e10, statePenaltyMPC_v09, statePenaltyMPC_re09, &info->dgap, &info->res_eq);
statePenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_8_26_26(params->C10, statePenaltyMPC_z09, statePenaltyMPC_D01, statePenaltyMPC_z10, params->e11, statePenaltyMPC_v10, statePenaltyMPC_re10, &info->dgap, &info->res_eq);
statePenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_8_26_26(params->C11, statePenaltyMPC_z10, statePenaltyMPC_D01, statePenaltyMPC_z11, params->e12, statePenaltyMPC_v11, statePenaltyMPC_re11, &info->dgap, &info->res_eq);
statePenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_8_26_26(params->C12, statePenaltyMPC_z11, statePenaltyMPC_D01, statePenaltyMPC_z12, params->e13, statePenaltyMPC_v12, statePenaltyMPC_re12, &info->dgap, &info->res_eq);
statePenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_8_26_26(params->C13, statePenaltyMPC_z12, statePenaltyMPC_D01, statePenaltyMPC_z13, params->e14, statePenaltyMPC_v13, statePenaltyMPC_re13, &info->dgap, &info->res_eq);
statePenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_8_26_8(params->C14, statePenaltyMPC_z13, statePenaltyMPC_D14, statePenaltyMPC_z14, params->e15, statePenaltyMPC_v14, statePenaltyMPC_re14, &info->dgap, &info->res_eq);
statePenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_8_26_8(params->C1, statePenaltyMPC_v01, statePenaltyMPC_D00, statePenaltyMPC_v00, statePenaltyMPC_grad_eq00);
statePenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_8_26_8(params->C2, statePenaltyMPC_v02, statePenaltyMPC_D01, statePenaltyMPC_v01, statePenaltyMPC_grad_eq01);
statePenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_8_26_8(params->C3, statePenaltyMPC_v03, statePenaltyMPC_D01, statePenaltyMPC_v02, statePenaltyMPC_grad_eq02);
statePenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_8_26_8(params->C4, statePenaltyMPC_v04, statePenaltyMPC_D01, statePenaltyMPC_v03, statePenaltyMPC_grad_eq03);
statePenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_8_26_8(params->C5, statePenaltyMPC_v05, statePenaltyMPC_D01, statePenaltyMPC_v04, statePenaltyMPC_grad_eq04);
statePenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_8_26_8(params->C6, statePenaltyMPC_v06, statePenaltyMPC_D01, statePenaltyMPC_v05, statePenaltyMPC_grad_eq05);
statePenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_8_26_8(params->C7, statePenaltyMPC_v07, statePenaltyMPC_D01, statePenaltyMPC_v06, statePenaltyMPC_grad_eq06);
statePenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_8_26_8(params->C8, statePenaltyMPC_v08, statePenaltyMPC_D01, statePenaltyMPC_v07, statePenaltyMPC_grad_eq07);
statePenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_8_26_8(params->C9, statePenaltyMPC_v09, statePenaltyMPC_D01, statePenaltyMPC_v08, statePenaltyMPC_grad_eq08);
statePenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_8_26_8(params->C10, statePenaltyMPC_v10, statePenaltyMPC_D01, statePenaltyMPC_v09, statePenaltyMPC_grad_eq09);
statePenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_8_26_8(params->C11, statePenaltyMPC_v11, statePenaltyMPC_D01, statePenaltyMPC_v10, statePenaltyMPC_grad_eq10);
statePenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_8_26_8(params->C12, statePenaltyMPC_v12, statePenaltyMPC_D01, statePenaltyMPC_v11, statePenaltyMPC_grad_eq11);
statePenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_8_26_8(params->C13, statePenaltyMPC_v13, statePenaltyMPC_D01, statePenaltyMPC_v12, statePenaltyMPC_grad_eq12);
statePenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_8_26_8(params->C14, statePenaltyMPC_v14, statePenaltyMPC_D01, statePenaltyMPC_v13, statePenaltyMPC_grad_eq13);
statePenaltyMPC_LA_DIAGZERO_MTVM_8_8(statePenaltyMPC_D14, statePenaltyMPC_v14, statePenaltyMPC_grad_eq14);
info->res_ineq = 0;
statePenaltyMPC_LA_VSUBADD3_26(params->lb1, statePenaltyMPC_z00, statePenaltyMPC_lbIdx00, statePenaltyMPC_llb00, statePenaltyMPC_slb00, statePenaltyMPC_rilb00, &info->dgap, &info->res_ineq);
statePenaltyMPC_LA_VSUBADD2_10(statePenaltyMPC_z00, statePenaltyMPC_ubIdx00, params->ub1, statePenaltyMPC_lub00, statePenaltyMPC_sub00, statePenaltyMPC_riub00, &info->dgap, &info->res_ineq);
statePenaltyMPC_LA_VSUBADD3_26(params->lb2, statePenaltyMPC_z01, statePenaltyMPC_lbIdx01, statePenaltyMPC_llb01, statePenaltyMPC_slb01, statePenaltyMPC_rilb01, &info->dgap, &info->res_ineq);
statePenaltyMPC_LA_VSUBADD2_10(statePenaltyMPC_z01, statePenaltyMPC_ubIdx01, params->ub2, statePenaltyMPC_lub01, statePenaltyMPC_sub01, statePenaltyMPC_riub01, &info->dgap, &info->res_ineq);
statePenaltyMPC_LA_VSUBADD3_26(params->lb3, statePenaltyMPC_z02, statePenaltyMPC_lbIdx02, statePenaltyMPC_llb02, statePenaltyMPC_slb02, statePenaltyMPC_rilb02, &info->dgap, &info->res_ineq);
statePenaltyMPC_LA_VSUBADD2_10(statePenaltyMPC_z02, statePenaltyMPC_ubIdx02, params->ub3, statePenaltyMPC_lub02, statePenaltyMPC_sub02, statePenaltyMPC_riub02, &info->dgap, &info->res_ineq);
statePenaltyMPC_LA_VSUBADD3_26(params->lb4, statePenaltyMPC_z03, statePenaltyMPC_lbIdx03, statePenaltyMPC_llb03, statePenaltyMPC_slb03, statePenaltyMPC_rilb03, &info->dgap, &info->res_ineq);
statePenaltyMPC_LA_VSUBADD2_10(statePenaltyMPC_z03, statePenaltyMPC_ubIdx03, params->ub4, statePenaltyMPC_lub03, statePenaltyMPC_sub03, statePenaltyMPC_riub03, &info->dgap, &info->res_ineq);
statePenaltyMPC_LA_VSUBADD3_26(params->lb5, statePenaltyMPC_z04, statePenaltyMPC_lbIdx04, statePenaltyMPC_llb04, statePenaltyMPC_slb04, statePenaltyMPC_rilb04, &info->dgap, &info->res_ineq);
statePenaltyMPC_LA_VSUBADD2_10(statePenaltyMPC_z04, statePenaltyMPC_ubIdx04, params->ub5, statePenaltyMPC_lub04, statePenaltyMPC_sub04, statePenaltyMPC_riub04, &info->dgap, &info->res_ineq);
statePenaltyMPC_LA_VSUBADD3_26(params->lb6, statePenaltyMPC_z05, statePenaltyMPC_lbIdx05, statePenaltyMPC_llb05, statePenaltyMPC_slb05, statePenaltyMPC_rilb05, &info->dgap, &info->res_ineq);
statePenaltyMPC_LA_VSUBADD2_10(statePenaltyMPC_z05, statePenaltyMPC_ubIdx05, params->ub6, statePenaltyMPC_lub05, statePenaltyMPC_sub05, statePenaltyMPC_riub05, &info->dgap, &info->res_ineq);
statePenaltyMPC_LA_VSUBADD3_26(params->lb7, statePenaltyMPC_z06, statePenaltyMPC_lbIdx06, statePenaltyMPC_llb06, statePenaltyMPC_slb06, statePenaltyMPC_rilb06, &info->dgap, &info->res_ineq);
statePenaltyMPC_LA_VSUBADD2_10(statePenaltyMPC_z06, statePenaltyMPC_ubIdx06, params->ub7, statePenaltyMPC_lub06, statePenaltyMPC_sub06, statePenaltyMPC_riub06, &info->dgap, &info->res_ineq);
statePenaltyMPC_LA_VSUBADD3_26(params->lb8, statePenaltyMPC_z07, statePenaltyMPC_lbIdx07, statePenaltyMPC_llb07, statePenaltyMPC_slb07, statePenaltyMPC_rilb07, &info->dgap, &info->res_ineq);
statePenaltyMPC_LA_VSUBADD2_10(statePenaltyMPC_z07, statePenaltyMPC_ubIdx07, params->ub8, statePenaltyMPC_lub07, statePenaltyMPC_sub07, statePenaltyMPC_riub07, &info->dgap, &info->res_ineq);
statePenaltyMPC_LA_VSUBADD3_26(params->lb9, statePenaltyMPC_z08, statePenaltyMPC_lbIdx08, statePenaltyMPC_llb08, statePenaltyMPC_slb08, statePenaltyMPC_rilb08, &info->dgap, &info->res_ineq);
statePenaltyMPC_LA_VSUBADD2_10(statePenaltyMPC_z08, statePenaltyMPC_ubIdx08, params->ub9, statePenaltyMPC_lub08, statePenaltyMPC_sub08, statePenaltyMPC_riub08, &info->dgap, &info->res_ineq);
statePenaltyMPC_LA_VSUBADD3_26(params->lb10, statePenaltyMPC_z09, statePenaltyMPC_lbIdx09, statePenaltyMPC_llb09, statePenaltyMPC_slb09, statePenaltyMPC_rilb09, &info->dgap, &info->res_ineq);
statePenaltyMPC_LA_VSUBADD2_10(statePenaltyMPC_z09, statePenaltyMPC_ubIdx09, params->ub10, statePenaltyMPC_lub09, statePenaltyMPC_sub09, statePenaltyMPC_riub09, &info->dgap, &info->res_ineq);
statePenaltyMPC_LA_VSUBADD3_26(params->lb11, statePenaltyMPC_z10, statePenaltyMPC_lbIdx10, statePenaltyMPC_llb10, statePenaltyMPC_slb10, statePenaltyMPC_rilb10, &info->dgap, &info->res_ineq);
statePenaltyMPC_LA_VSUBADD2_10(statePenaltyMPC_z10, statePenaltyMPC_ubIdx10, params->ub11, statePenaltyMPC_lub10, statePenaltyMPC_sub10, statePenaltyMPC_riub10, &info->dgap, &info->res_ineq);
statePenaltyMPC_LA_VSUBADD3_26(params->lb12, statePenaltyMPC_z11, statePenaltyMPC_lbIdx11, statePenaltyMPC_llb11, statePenaltyMPC_slb11, statePenaltyMPC_rilb11, &info->dgap, &info->res_ineq);
statePenaltyMPC_LA_VSUBADD2_10(statePenaltyMPC_z11, statePenaltyMPC_ubIdx11, params->ub12, statePenaltyMPC_lub11, statePenaltyMPC_sub11, statePenaltyMPC_riub11, &info->dgap, &info->res_ineq);
statePenaltyMPC_LA_VSUBADD3_26(params->lb13, statePenaltyMPC_z12, statePenaltyMPC_lbIdx12, statePenaltyMPC_llb12, statePenaltyMPC_slb12, statePenaltyMPC_rilb12, &info->dgap, &info->res_ineq);
statePenaltyMPC_LA_VSUBADD2_10(statePenaltyMPC_z12, statePenaltyMPC_ubIdx12, params->ub13, statePenaltyMPC_lub12, statePenaltyMPC_sub12, statePenaltyMPC_riub12, &info->dgap, &info->res_ineq);
statePenaltyMPC_LA_VSUBADD3_26(params->lb14, statePenaltyMPC_z13, statePenaltyMPC_lbIdx13, statePenaltyMPC_llb13, statePenaltyMPC_slb13, statePenaltyMPC_rilb13, &info->dgap, &info->res_ineq);
statePenaltyMPC_LA_VSUBADD2_10(statePenaltyMPC_z13, statePenaltyMPC_ubIdx13, params->ub14, statePenaltyMPC_lub13, statePenaltyMPC_sub13, statePenaltyMPC_riub13, &info->dgap, &info->res_ineq);
statePenaltyMPC_LA_VSUBADD3_8(params->lb15, statePenaltyMPC_z14, statePenaltyMPC_lbIdx14, statePenaltyMPC_llb14, statePenaltyMPC_slb14, statePenaltyMPC_rilb14, &info->dgap, &info->res_ineq);
statePenaltyMPC_LA_VSUBADD2_8(statePenaltyMPC_z14, statePenaltyMPC_ubIdx14, params->ub15, statePenaltyMPC_lub14, statePenaltyMPC_sub14, statePenaltyMPC_riub14, &info->dgap, &info->res_ineq);
statePenaltyMPC_LA_INEQ_B_GRAD_26_26_10(statePenaltyMPC_lub00, statePenaltyMPC_sub00, statePenaltyMPC_riub00, statePenaltyMPC_llb00, statePenaltyMPC_slb00, statePenaltyMPC_rilb00, statePenaltyMPC_lbIdx00, statePenaltyMPC_ubIdx00, statePenaltyMPC_grad_ineq00, statePenaltyMPC_lubbysub00, statePenaltyMPC_llbbyslb00);
statePenaltyMPC_LA_INEQ_B_GRAD_26_26_10(statePenaltyMPC_lub01, statePenaltyMPC_sub01, statePenaltyMPC_riub01, statePenaltyMPC_llb01, statePenaltyMPC_slb01, statePenaltyMPC_rilb01, statePenaltyMPC_lbIdx01, statePenaltyMPC_ubIdx01, statePenaltyMPC_grad_ineq01, statePenaltyMPC_lubbysub01, statePenaltyMPC_llbbyslb01);
statePenaltyMPC_LA_INEQ_B_GRAD_26_26_10(statePenaltyMPC_lub02, statePenaltyMPC_sub02, statePenaltyMPC_riub02, statePenaltyMPC_llb02, statePenaltyMPC_slb02, statePenaltyMPC_rilb02, statePenaltyMPC_lbIdx02, statePenaltyMPC_ubIdx02, statePenaltyMPC_grad_ineq02, statePenaltyMPC_lubbysub02, statePenaltyMPC_llbbyslb02);
statePenaltyMPC_LA_INEQ_B_GRAD_26_26_10(statePenaltyMPC_lub03, statePenaltyMPC_sub03, statePenaltyMPC_riub03, statePenaltyMPC_llb03, statePenaltyMPC_slb03, statePenaltyMPC_rilb03, statePenaltyMPC_lbIdx03, statePenaltyMPC_ubIdx03, statePenaltyMPC_grad_ineq03, statePenaltyMPC_lubbysub03, statePenaltyMPC_llbbyslb03);
statePenaltyMPC_LA_INEQ_B_GRAD_26_26_10(statePenaltyMPC_lub04, statePenaltyMPC_sub04, statePenaltyMPC_riub04, statePenaltyMPC_llb04, statePenaltyMPC_slb04, statePenaltyMPC_rilb04, statePenaltyMPC_lbIdx04, statePenaltyMPC_ubIdx04, statePenaltyMPC_grad_ineq04, statePenaltyMPC_lubbysub04, statePenaltyMPC_llbbyslb04);
statePenaltyMPC_LA_INEQ_B_GRAD_26_26_10(statePenaltyMPC_lub05, statePenaltyMPC_sub05, statePenaltyMPC_riub05, statePenaltyMPC_llb05, statePenaltyMPC_slb05, statePenaltyMPC_rilb05, statePenaltyMPC_lbIdx05, statePenaltyMPC_ubIdx05, statePenaltyMPC_grad_ineq05, statePenaltyMPC_lubbysub05, statePenaltyMPC_llbbyslb05);
statePenaltyMPC_LA_INEQ_B_GRAD_26_26_10(statePenaltyMPC_lub06, statePenaltyMPC_sub06, statePenaltyMPC_riub06, statePenaltyMPC_llb06, statePenaltyMPC_slb06, statePenaltyMPC_rilb06, statePenaltyMPC_lbIdx06, statePenaltyMPC_ubIdx06, statePenaltyMPC_grad_ineq06, statePenaltyMPC_lubbysub06, statePenaltyMPC_llbbyslb06);
statePenaltyMPC_LA_INEQ_B_GRAD_26_26_10(statePenaltyMPC_lub07, statePenaltyMPC_sub07, statePenaltyMPC_riub07, statePenaltyMPC_llb07, statePenaltyMPC_slb07, statePenaltyMPC_rilb07, statePenaltyMPC_lbIdx07, statePenaltyMPC_ubIdx07, statePenaltyMPC_grad_ineq07, statePenaltyMPC_lubbysub07, statePenaltyMPC_llbbyslb07);
statePenaltyMPC_LA_INEQ_B_GRAD_26_26_10(statePenaltyMPC_lub08, statePenaltyMPC_sub08, statePenaltyMPC_riub08, statePenaltyMPC_llb08, statePenaltyMPC_slb08, statePenaltyMPC_rilb08, statePenaltyMPC_lbIdx08, statePenaltyMPC_ubIdx08, statePenaltyMPC_grad_ineq08, statePenaltyMPC_lubbysub08, statePenaltyMPC_llbbyslb08);
statePenaltyMPC_LA_INEQ_B_GRAD_26_26_10(statePenaltyMPC_lub09, statePenaltyMPC_sub09, statePenaltyMPC_riub09, statePenaltyMPC_llb09, statePenaltyMPC_slb09, statePenaltyMPC_rilb09, statePenaltyMPC_lbIdx09, statePenaltyMPC_ubIdx09, statePenaltyMPC_grad_ineq09, statePenaltyMPC_lubbysub09, statePenaltyMPC_llbbyslb09);
statePenaltyMPC_LA_INEQ_B_GRAD_26_26_10(statePenaltyMPC_lub10, statePenaltyMPC_sub10, statePenaltyMPC_riub10, statePenaltyMPC_llb10, statePenaltyMPC_slb10, statePenaltyMPC_rilb10, statePenaltyMPC_lbIdx10, statePenaltyMPC_ubIdx10, statePenaltyMPC_grad_ineq10, statePenaltyMPC_lubbysub10, statePenaltyMPC_llbbyslb10);
statePenaltyMPC_LA_INEQ_B_GRAD_26_26_10(statePenaltyMPC_lub11, statePenaltyMPC_sub11, statePenaltyMPC_riub11, statePenaltyMPC_llb11, statePenaltyMPC_slb11, statePenaltyMPC_rilb11, statePenaltyMPC_lbIdx11, statePenaltyMPC_ubIdx11, statePenaltyMPC_grad_ineq11, statePenaltyMPC_lubbysub11, statePenaltyMPC_llbbyslb11);
statePenaltyMPC_LA_INEQ_B_GRAD_26_26_10(statePenaltyMPC_lub12, statePenaltyMPC_sub12, statePenaltyMPC_riub12, statePenaltyMPC_llb12, statePenaltyMPC_slb12, statePenaltyMPC_rilb12, statePenaltyMPC_lbIdx12, statePenaltyMPC_ubIdx12, statePenaltyMPC_grad_ineq12, statePenaltyMPC_lubbysub12, statePenaltyMPC_llbbyslb12);
statePenaltyMPC_LA_INEQ_B_GRAD_26_26_10(statePenaltyMPC_lub13, statePenaltyMPC_sub13, statePenaltyMPC_riub13, statePenaltyMPC_llb13, statePenaltyMPC_slb13, statePenaltyMPC_rilb13, statePenaltyMPC_lbIdx13, statePenaltyMPC_ubIdx13, statePenaltyMPC_grad_ineq13, statePenaltyMPC_lubbysub13, statePenaltyMPC_llbbyslb13);
statePenaltyMPC_LA_INEQ_B_GRAD_8_8_8(statePenaltyMPC_lub14, statePenaltyMPC_sub14, statePenaltyMPC_riub14, statePenaltyMPC_llb14, statePenaltyMPC_slb14, statePenaltyMPC_rilb14, statePenaltyMPC_lbIdx14, statePenaltyMPC_ubIdx14, statePenaltyMPC_grad_ineq14, statePenaltyMPC_lubbysub14, statePenaltyMPC_llbbyslb14);
info->dobj = info->pobj - info->dgap;
info->rdgap = info->pobj ? info->dgap / info->pobj : 1e6;
if( info->rdgap < 0 ) info->rdgap = -info->rdgap;
PRINTTEXT("  %3d  %3.1e  %3.1e  %+6.4e  %+6.4e  %+3.1e  %3.1e  %3.1e\n",info->it, info->res_eq, info->res_ineq, info->pobj, info->dobj, info->dgap, info->rdgap, info->mu);
if( info->mu < statePenaltyMPC_SET_ACC_KKTCOMPL
    && (info->rdgap < statePenaltyMPC_SET_ACC_RDGAP || info->dgap < statePenaltyMPC_SET_ACC_KKTCOMPL)
    && info->res_eq < statePenaltyMPC_SET_ACC_RESEQ
    && info->res_ineq < statePenaltyMPC_SET_ACC_RESINEQ ){
PRINTTEXT("OPTIMAL (within RESEQ=%2.1e, RESINEQ=%2.1e, (R)DGAP=(%2.1e)%2.1e, MU=%2.1e).\n",statePenaltyMPC_SET_ACC_RESEQ, statePenaltyMPC_SET_ACC_RESINEQ,statePenaltyMPC_SET_ACC_KKTCOMPL,statePenaltyMPC_SET_ACC_RDGAP,statePenaltyMPC_SET_ACC_KKTCOMPL);
exitcode = statePenaltyMPC_OPTIMAL; break; }
if( info->it == statePenaltyMPC_SET_MAXIT ){
PRINTTEXT("Maximum number of iterations reached, exiting.\n");
exitcode = statePenaltyMPC_MAXITREACHED; break; }
statePenaltyMPC_LA_VVADD3_372(statePenaltyMPC_grad_cost, statePenaltyMPC_grad_eq, statePenaltyMPC_grad_ineq, statePenaltyMPC_rd);
statePenaltyMPC_LA_DIAG_CHOL_LBUB_26_26_10(params->H1, statePenaltyMPC_llbbyslb00, statePenaltyMPC_lbIdx00, statePenaltyMPC_lubbysub00, statePenaltyMPC_ubIdx00, statePenaltyMPC_Phi00);
statePenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_8_26(statePenaltyMPC_Phi00, params->C1, statePenaltyMPC_V00);
statePenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_8_26(statePenaltyMPC_Phi00, statePenaltyMPC_D00, statePenaltyMPC_W00);
statePenaltyMPC_LA_DENSE_DIAGZERO_MMTM_8_26_8(statePenaltyMPC_W00, statePenaltyMPC_V00, statePenaltyMPC_Ysd01);
statePenaltyMPC_LA_DIAG_FORWARDSUB_26(statePenaltyMPC_Phi00, statePenaltyMPC_rd00, statePenaltyMPC_Lbyrd00);
statePenaltyMPC_LA_DIAG_CHOL_LBUB_26_26_10(params->H2, statePenaltyMPC_llbbyslb01, statePenaltyMPC_lbIdx01, statePenaltyMPC_lubbysub01, statePenaltyMPC_ubIdx01, statePenaltyMPC_Phi01);
statePenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_8_26(statePenaltyMPC_Phi01, params->C2, statePenaltyMPC_V01);
statePenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_8_26(statePenaltyMPC_Phi01, statePenaltyMPC_D01, statePenaltyMPC_W01);
statePenaltyMPC_LA_DENSE_DIAGZERO_MMTM_8_26_8(statePenaltyMPC_W01, statePenaltyMPC_V01, statePenaltyMPC_Ysd02);
statePenaltyMPC_LA_DIAG_FORWARDSUB_26(statePenaltyMPC_Phi01, statePenaltyMPC_rd01, statePenaltyMPC_Lbyrd01);
statePenaltyMPC_LA_DIAG_CHOL_LBUB_26_26_10(params->H3, statePenaltyMPC_llbbyslb02, statePenaltyMPC_lbIdx02, statePenaltyMPC_lubbysub02, statePenaltyMPC_ubIdx02, statePenaltyMPC_Phi02);
statePenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_8_26(statePenaltyMPC_Phi02, params->C3, statePenaltyMPC_V02);
statePenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_8_26(statePenaltyMPC_Phi02, statePenaltyMPC_D01, statePenaltyMPC_W02);
statePenaltyMPC_LA_DENSE_DIAGZERO_MMTM_8_26_8(statePenaltyMPC_W02, statePenaltyMPC_V02, statePenaltyMPC_Ysd03);
statePenaltyMPC_LA_DIAG_FORWARDSUB_26(statePenaltyMPC_Phi02, statePenaltyMPC_rd02, statePenaltyMPC_Lbyrd02);
statePenaltyMPC_LA_DIAG_CHOL_LBUB_26_26_10(params->H4, statePenaltyMPC_llbbyslb03, statePenaltyMPC_lbIdx03, statePenaltyMPC_lubbysub03, statePenaltyMPC_ubIdx03, statePenaltyMPC_Phi03);
statePenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_8_26(statePenaltyMPC_Phi03, params->C4, statePenaltyMPC_V03);
statePenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_8_26(statePenaltyMPC_Phi03, statePenaltyMPC_D01, statePenaltyMPC_W03);
statePenaltyMPC_LA_DENSE_DIAGZERO_MMTM_8_26_8(statePenaltyMPC_W03, statePenaltyMPC_V03, statePenaltyMPC_Ysd04);
statePenaltyMPC_LA_DIAG_FORWARDSUB_26(statePenaltyMPC_Phi03, statePenaltyMPC_rd03, statePenaltyMPC_Lbyrd03);
statePenaltyMPC_LA_DIAG_CHOL_LBUB_26_26_10(params->H5, statePenaltyMPC_llbbyslb04, statePenaltyMPC_lbIdx04, statePenaltyMPC_lubbysub04, statePenaltyMPC_ubIdx04, statePenaltyMPC_Phi04);
statePenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_8_26(statePenaltyMPC_Phi04, params->C5, statePenaltyMPC_V04);
statePenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_8_26(statePenaltyMPC_Phi04, statePenaltyMPC_D01, statePenaltyMPC_W04);
statePenaltyMPC_LA_DENSE_DIAGZERO_MMTM_8_26_8(statePenaltyMPC_W04, statePenaltyMPC_V04, statePenaltyMPC_Ysd05);
statePenaltyMPC_LA_DIAG_FORWARDSUB_26(statePenaltyMPC_Phi04, statePenaltyMPC_rd04, statePenaltyMPC_Lbyrd04);
statePenaltyMPC_LA_DIAG_CHOL_LBUB_26_26_10(params->H6, statePenaltyMPC_llbbyslb05, statePenaltyMPC_lbIdx05, statePenaltyMPC_lubbysub05, statePenaltyMPC_ubIdx05, statePenaltyMPC_Phi05);
statePenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_8_26(statePenaltyMPC_Phi05, params->C6, statePenaltyMPC_V05);
statePenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_8_26(statePenaltyMPC_Phi05, statePenaltyMPC_D01, statePenaltyMPC_W05);
statePenaltyMPC_LA_DENSE_DIAGZERO_MMTM_8_26_8(statePenaltyMPC_W05, statePenaltyMPC_V05, statePenaltyMPC_Ysd06);
statePenaltyMPC_LA_DIAG_FORWARDSUB_26(statePenaltyMPC_Phi05, statePenaltyMPC_rd05, statePenaltyMPC_Lbyrd05);
statePenaltyMPC_LA_DIAG_CHOL_LBUB_26_26_10(params->H7, statePenaltyMPC_llbbyslb06, statePenaltyMPC_lbIdx06, statePenaltyMPC_lubbysub06, statePenaltyMPC_ubIdx06, statePenaltyMPC_Phi06);
statePenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_8_26(statePenaltyMPC_Phi06, params->C7, statePenaltyMPC_V06);
statePenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_8_26(statePenaltyMPC_Phi06, statePenaltyMPC_D01, statePenaltyMPC_W06);
statePenaltyMPC_LA_DENSE_DIAGZERO_MMTM_8_26_8(statePenaltyMPC_W06, statePenaltyMPC_V06, statePenaltyMPC_Ysd07);
statePenaltyMPC_LA_DIAG_FORWARDSUB_26(statePenaltyMPC_Phi06, statePenaltyMPC_rd06, statePenaltyMPC_Lbyrd06);
statePenaltyMPC_LA_DIAG_CHOL_LBUB_26_26_10(params->H8, statePenaltyMPC_llbbyslb07, statePenaltyMPC_lbIdx07, statePenaltyMPC_lubbysub07, statePenaltyMPC_ubIdx07, statePenaltyMPC_Phi07);
statePenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_8_26(statePenaltyMPC_Phi07, params->C8, statePenaltyMPC_V07);
statePenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_8_26(statePenaltyMPC_Phi07, statePenaltyMPC_D01, statePenaltyMPC_W07);
statePenaltyMPC_LA_DENSE_DIAGZERO_MMTM_8_26_8(statePenaltyMPC_W07, statePenaltyMPC_V07, statePenaltyMPC_Ysd08);
statePenaltyMPC_LA_DIAG_FORWARDSUB_26(statePenaltyMPC_Phi07, statePenaltyMPC_rd07, statePenaltyMPC_Lbyrd07);
statePenaltyMPC_LA_DIAG_CHOL_LBUB_26_26_10(params->H9, statePenaltyMPC_llbbyslb08, statePenaltyMPC_lbIdx08, statePenaltyMPC_lubbysub08, statePenaltyMPC_ubIdx08, statePenaltyMPC_Phi08);
statePenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_8_26(statePenaltyMPC_Phi08, params->C9, statePenaltyMPC_V08);
statePenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_8_26(statePenaltyMPC_Phi08, statePenaltyMPC_D01, statePenaltyMPC_W08);
statePenaltyMPC_LA_DENSE_DIAGZERO_MMTM_8_26_8(statePenaltyMPC_W08, statePenaltyMPC_V08, statePenaltyMPC_Ysd09);
statePenaltyMPC_LA_DIAG_FORWARDSUB_26(statePenaltyMPC_Phi08, statePenaltyMPC_rd08, statePenaltyMPC_Lbyrd08);
statePenaltyMPC_LA_DIAG_CHOL_LBUB_26_26_10(params->H10, statePenaltyMPC_llbbyslb09, statePenaltyMPC_lbIdx09, statePenaltyMPC_lubbysub09, statePenaltyMPC_ubIdx09, statePenaltyMPC_Phi09);
statePenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_8_26(statePenaltyMPC_Phi09, params->C10, statePenaltyMPC_V09);
statePenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_8_26(statePenaltyMPC_Phi09, statePenaltyMPC_D01, statePenaltyMPC_W09);
statePenaltyMPC_LA_DENSE_DIAGZERO_MMTM_8_26_8(statePenaltyMPC_W09, statePenaltyMPC_V09, statePenaltyMPC_Ysd10);
statePenaltyMPC_LA_DIAG_FORWARDSUB_26(statePenaltyMPC_Phi09, statePenaltyMPC_rd09, statePenaltyMPC_Lbyrd09);
statePenaltyMPC_LA_DIAG_CHOL_LBUB_26_26_10(params->H11, statePenaltyMPC_llbbyslb10, statePenaltyMPC_lbIdx10, statePenaltyMPC_lubbysub10, statePenaltyMPC_ubIdx10, statePenaltyMPC_Phi10);
statePenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_8_26(statePenaltyMPC_Phi10, params->C11, statePenaltyMPC_V10);
statePenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_8_26(statePenaltyMPC_Phi10, statePenaltyMPC_D01, statePenaltyMPC_W10);
statePenaltyMPC_LA_DENSE_DIAGZERO_MMTM_8_26_8(statePenaltyMPC_W10, statePenaltyMPC_V10, statePenaltyMPC_Ysd11);
statePenaltyMPC_LA_DIAG_FORWARDSUB_26(statePenaltyMPC_Phi10, statePenaltyMPC_rd10, statePenaltyMPC_Lbyrd10);
statePenaltyMPC_LA_DIAG_CHOL_LBUB_26_26_10(params->H12, statePenaltyMPC_llbbyslb11, statePenaltyMPC_lbIdx11, statePenaltyMPC_lubbysub11, statePenaltyMPC_ubIdx11, statePenaltyMPC_Phi11);
statePenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_8_26(statePenaltyMPC_Phi11, params->C12, statePenaltyMPC_V11);
statePenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_8_26(statePenaltyMPC_Phi11, statePenaltyMPC_D01, statePenaltyMPC_W11);
statePenaltyMPC_LA_DENSE_DIAGZERO_MMTM_8_26_8(statePenaltyMPC_W11, statePenaltyMPC_V11, statePenaltyMPC_Ysd12);
statePenaltyMPC_LA_DIAG_FORWARDSUB_26(statePenaltyMPC_Phi11, statePenaltyMPC_rd11, statePenaltyMPC_Lbyrd11);
statePenaltyMPC_LA_DIAG_CHOL_LBUB_26_26_10(params->H13, statePenaltyMPC_llbbyslb12, statePenaltyMPC_lbIdx12, statePenaltyMPC_lubbysub12, statePenaltyMPC_ubIdx12, statePenaltyMPC_Phi12);
statePenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_8_26(statePenaltyMPC_Phi12, params->C13, statePenaltyMPC_V12);
statePenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_8_26(statePenaltyMPC_Phi12, statePenaltyMPC_D01, statePenaltyMPC_W12);
statePenaltyMPC_LA_DENSE_DIAGZERO_MMTM_8_26_8(statePenaltyMPC_W12, statePenaltyMPC_V12, statePenaltyMPC_Ysd13);
statePenaltyMPC_LA_DIAG_FORWARDSUB_26(statePenaltyMPC_Phi12, statePenaltyMPC_rd12, statePenaltyMPC_Lbyrd12);
statePenaltyMPC_LA_DIAG_CHOL_LBUB_26_26_10(params->H14, statePenaltyMPC_llbbyslb13, statePenaltyMPC_lbIdx13, statePenaltyMPC_lubbysub13, statePenaltyMPC_ubIdx13, statePenaltyMPC_Phi13);
statePenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_8_26(statePenaltyMPC_Phi13, params->C14, statePenaltyMPC_V13);
statePenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_8_26(statePenaltyMPC_Phi13, statePenaltyMPC_D01, statePenaltyMPC_W13);
statePenaltyMPC_LA_DENSE_DIAGZERO_MMTM_8_26_8(statePenaltyMPC_W13, statePenaltyMPC_V13, statePenaltyMPC_Ysd14);
statePenaltyMPC_LA_DIAG_FORWARDSUB_26(statePenaltyMPC_Phi13, statePenaltyMPC_rd13, statePenaltyMPC_Lbyrd13);
statePenaltyMPC_LA_DIAG_CHOL_ONELOOP_LBUB_8_8_8(params->H15, statePenaltyMPC_llbbyslb14, statePenaltyMPC_lbIdx14, statePenaltyMPC_lubbysub14, statePenaltyMPC_ubIdx14, statePenaltyMPC_Phi14);
statePenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_8_8(statePenaltyMPC_Phi14, statePenaltyMPC_D14, statePenaltyMPC_W14);
statePenaltyMPC_LA_DIAG_FORWARDSUB_8(statePenaltyMPC_Phi14, statePenaltyMPC_rd14, statePenaltyMPC_Lbyrd14);
statePenaltyMPC_LA_DIAGZERO_MMT_8(statePenaltyMPC_W00, statePenaltyMPC_Yd00);
statePenaltyMPC_LA_DIAGZERO_MVMSUB7_8(statePenaltyMPC_W00, statePenaltyMPC_Lbyrd00, statePenaltyMPC_re00, statePenaltyMPC_beta00);
statePenaltyMPC_LA_DENSE_DIAGZERO_MMT2_8_26_26(statePenaltyMPC_V00, statePenaltyMPC_W01, statePenaltyMPC_Yd01);
statePenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_8_26_26(statePenaltyMPC_V00, statePenaltyMPC_Lbyrd00, statePenaltyMPC_W01, statePenaltyMPC_Lbyrd01, statePenaltyMPC_re01, statePenaltyMPC_beta01);
statePenaltyMPC_LA_DENSE_DIAGZERO_MMT2_8_26_26(statePenaltyMPC_V01, statePenaltyMPC_W02, statePenaltyMPC_Yd02);
statePenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_8_26_26(statePenaltyMPC_V01, statePenaltyMPC_Lbyrd01, statePenaltyMPC_W02, statePenaltyMPC_Lbyrd02, statePenaltyMPC_re02, statePenaltyMPC_beta02);
statePenaltyMPC_LA_DENSE_DIAGZERO_MMT2_8_26_26(statePenaltyMPC_V02, statePenaltyMPC_W03, statePenaltyMPC_Yd03);
statePenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_8_26_26(statePenaltyMPC_V02, statePenaltyMPC_Lbyrd02, statePenaltyMPC_W03, statePenaltyMPC_Lbyrd03, statePenaltyMPC_re03, statePenaltyMPC_beta03);
statePenaltyMPC_LA_DENSE_DIAGZERO_MMT2_8_26_26(statePenaltyMPC_V03, statePenaltyMPC_W04, statePenaltyMPC_Yd04);
statePenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_8_26_26(statePenaltyMPC_V03, statePenaltyMPC_Lbyrd03, statePenaltyMPC_W04, statePenaltyMPC_Lbyrd04, statePenaltyMPC_re04, statePenaltyMPC_beta04);
statePenaltyMPC_LA_DENSE_DIAGZERO_MMT2_8_26_26(statePenaltyMPC_V04, statePenaltyMPC_W05, statePenaltyMPC_Yd05);
statePenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_8_26_26(statePenaltyMPC_V04, statePenaltyMPC_Lbyrd04, statePenaltyMPC_W05, statePenaltyMPC_Lbyrd05, statePenaltyMPC_re05, statePenaltyMPC_beta05);
statePenaltyMPC_LA_DENSE_DIAGZERO_MMT2_8_26_26(statePenaltyMPC_V05, statePenaltyMPC_W06, statePenaltyMPC_Yd06);
statePenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_8_26_26(statePenaltyMPC_V05, statePenaltyMPC_Lbyrd05, statePenaltyMPC_W06, statePenaltyMPC_Lbyrd06, statePenaltyMPC_re06, statePenaltyMPC_beta06);
statePenaltyMPC_LA_DENSE_DIAGZERO_MMT2_8_26_26(statePenaltyMPC_V06, statePenaltyMPC_W07, statePenaltyMPC_Yd07);
statePenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_8_26_26(statePenaltyMPC_V06, statePenaltyMPC_Lbyrd06, statePenaltyMPC_W07, statePenaltyMPC_Lbyrd07, statePenaltyMPC_re07, statePenaltyMPC_beta07);
statePenaltyMPC_LA_DENSE_DIAGZERO_MMT2_8_26_26(statePenaltyMPC_V07, statePenaltyMPC_W08, statePenaltyMPC_Yd08);
statePenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_8_26_26(statePenaltyMPC_V07, statePenaltyMPC_Lbyrd07, statePenaltyMPC_W08, statePenaltyMPC_Lbyrd08, statePenaltyMPC_re08, statePenaltyMPC_beta08);
statePenaltyMPC_LA_DENSE_DIAGZERO_MMT2_8_26_26(statePenaltyMPC_V08, statePenaltyMPC_W09, statePenaltyMPC_Yd09);
statePenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_8_26_26(statePenaltyMPC_V08, statePenaltyMPC_Lbyrd08, statePenaltyMPC_W09, statePenaltyMPC_Lbyrd09, statePenaltyMPC_re09, statePenaltyMPC_beta09);
statePenaltyMPC_LA_DENSE_DIAGZERO_MMT2_8_26_26(statePenaltyMPC_V09, statePenaltyMPC_W10, statePenaltyMPC_Yd10);
statePenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_8_26_26(statePenaltyMPC_V09, statePenaltyMPC_Lbyrd09, statePenaltyMPC_W10, statePenaltyMPC_Lbyrd10, statePenaltyMPC_re10, statePenaltyMPC_beta10);
statePenaltyMPC_LA_DENSE_DIAGZERO_MMT2_8_26_26(statePenaltyMPC_V10, statePenaltyMPC_W11, statePenaltyMPC_Yd11);
statePenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_8_26_26(statePenaltyMPC_V10, statePenaltyMPC_Lbyrd10, statePenaltyMPC_W11, statePenaltyMPC_Lbyrd11, statePenaltyMPC_re11, statePenaltyMPC_beta11);
statePenaltyMPC_LA_DENSE_DIAGZERO_MMT2_8_26_26(statePenaltyMPC_V11, statePenaltyMPC_W12, statePenaltyMPC_Yd12);
statePenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_8_26_26(statePenaltyMPC_V11, statePenaltyMPC_Lbyrd11, statePenaltyMPC_W12, statePenaltyMPC_Lbyrd12, statePenaltyMPC_re12, statePenaltyMPC_beta12);
statePenaltyMPC_LA_DENSE_DIAGZERO_MMT2_8_26_26(statePenaltyMPC_V12, statePenaltyMPC_W13, statePenaltyMPC_Yd13);
statePenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_8_26_26(statePenaltyMPC_V12, statePenaltyMPC_Lbyrd12, statePenaltyMPC_W13, statePenaltyMPC_Lbyrd13, statePenaltyMPC_re13, statePenaltyMPC_beta13);
statePenaltyMPC_LA_DENSE_DIAGZERO_MMT2_8_26_8(statePenaltyMPC_V13, statePenaltyMPC_W14, statePenaltyMPC_Yd14);
statePenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_8_26_8(statePenaltyMPC_V13, statePenaltyMPC_Lbyrd13, statePenaltyMPC_W14, statePenaltyMPC_Lbyrd14, statePenaltyMPC_re14, statePenaltyMPC_beta14);
statePenaltyMPC_LA_DENSE_CHOL_8(statePenaltyMPC_Yd00, statePenaltyMPC_Ld00);
statePenaltyMPC_LA_DENSE_FORWARDSUB_8(statePenaltyMPC_Ld00, statePenaltyMPC_beta00, statePenaltyMPC_yy00);
statePenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_8_8(statePenaltyMPC_Ld00, statePenaltyMPC_Ysd01, statePenaltyMPC_Lsd01);
statePenaltyMPC_LA_DENSE_MMTSUB_8_8(statePenaltyMPC_Lsd01, statePenaltyMPC_Yd01);
statePenaltyMPC_LA_DENSE_CHOL_8(statePenaltyMPC_Yd01, statePenaltyMPC_Ld01);
statePenaltyMPC_LA_DENSE_MVMSUB1_8_8(statePenaltyMPC_Lsd01, statePenaltyMPC_yy00, statePenaltyMPC_beta01, statePenaltyMPC_bmy01);
statePenaltyMPC_LA_DENSE_FORWARDSUB_8(statePenaltyMPC_Ld01, statePenaltyMPC_bmy01, statePenaltyMPC_yy01);
statePenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_8_8(statePenaltyMPC_Ld01, statePenaltyMPC_Ysd02, statePenaltyMPC_Lsd02);
statePenaltyMPC_LA_DENSE_MMTSUB_8_8(statePenaltyMPC_Lsd02, statePenaltyMPC_Yd02);
statePenaltyMPC_LA_DENSE_CHOL_8(statePenaltyMPC_Yd02, statePenaltyMPC_Ld02);
statePenaltyMPC_LA_DENSE_MVMSUB1_8_8(statePenaltyMPC_Lsd02, statePenaltyMPC_yy01, statePenaltyMPC_beta02, statePenaltyMPC_bmy02);
statePenaltyMPC_LA_DENSE_FORWARDSUB_8(statePenaltyMPC_Ld02, statePenaltyMPC_bmy02, statePenaltyMPC_yy02);
statePenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_8_8(statePenaltyMPC_Ld02, statePenaltyMPC_Ysd03, statePenaltyMPC_Lsd03);
statePenaltyMPC_LA_DENSE_MMTSUB_8_8(statePenaltyMPC_Lsd03, statePenaltyMPC_Yd03);
statePenaltyMPC_LA_DENSE_CHOL_8(statePenaltyMPC_Yd03, statePenaltyMPC_Ld03);
statePenaltyMPC_LA_DENSE_MVMSUB1_8_8(statePenaltyMPC_Lsd03, statePenaltyMPC_yy02, statePenaltyMPC_beta03, statePenaltyMPC_bmy03);
statePenaltyMPC_LA_DENSE_FORWARDSUB_8(statePenaltyMPC_Ld03, statePenaltyMPC_bmy03, statePenaltyMPC_yy03);
statePenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_8_8(statePenaltyMPC_Ld03, statePenaltyMPC_Ysd04, statePenaltyMPC_Lsd04);
statePenaltyMPC_LA_DENSE_MMTSUB_8_8(statePenaltyMPC_Lsd04, statePenaltyMPC_Yd04);
statePenaltyMPC_LA_DENSE_CHOL_8(statePenaltyMPC_Yd04, statePenaltyMPC_Ld04);
statePenaltyMPC_LA_DENSE_MVMSUB1_8_8(statePenaltyMPC_Lsd04, statePenaltyMPC_yy03, statePenaltyMPC_beta04, statePenaltyMPC_bmy04);
statePenaltyMPC_LA_DENSE_FORWARDSUB_8(statePenaltyMPC_Ld04, statePenaltyMPC_bmy04, statePenaltyMPC_yy04);
statePenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_8_8(statePenaltyMPC_Ld04, statePenaltyMPC_Ysd05, statePenaltyMPC_Lsd05);
statePenaltyMPC_LA_DENSE_MMTSUB_8_8(statePenaltyMPC_Lsd05, statePenaltyMPC_Yd05);
statePenaltyMPC_LA_DENSE_CHOL_8(statePenaltyMPC_Yd05, statePenaltyMPC_Ld05);
statePenaltyMPC_LA_DENSE_MVMSUB1_8_8(statePenaltyMPC_Lsd05, statePenaltyMPC_yy04, statePenaltyMPC_beta05, statePenaltyMPC_bmy05);
statePenaltyMPC_LA_DENSE_FORWARDSUB_8(statePenaltyMPC_Ld05, statePenaltyMPC_bmy05, statePenaltyMPC_yy05);
statePenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_8_8(statePenaltyMPC_Ld05, statePenaltyMPC_Ysd06, statePenaltyMPC_Lsd06);
statePenaltyMPC_LA_DENSE_MMTSUB_8_8(statePenaltyMPC_Lsd06, statePenaltyMPC_Yd06);
statePenaltyMPC_LA_DENSE_CHOL_8(statePenaltyMPC_Yd06, statePenaltyMPC_Ld06);
statePenaltyMPC_LA_DENSE_MVMSUB1_8_8(statePenaltyMPC_Lsd06, statePenaltyMPC_yy05, statePenaltyMPC_beta06, statePenaltyMPC_bmy06);
statePenaltyMPC_LA_DENSE_FORWARDSUB_8(statePenaltyMPC_Ld06, statePenaltyMPC_bmy06, statePenaltyMPC_yy06);
statePenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_8_8(statePenaltyMPC_Ld06, statePenaltyMPC_Ysd07, statePenaltyMPC_Lsd07);
statePenaltyMPC_LA_DENSE_MMTSUB_8_8(statePenaltyMPC_Lsd07, statePenaltyMPC_Yd07);
statePenaltyMPC_LA_DENSE_CHOL_8(statePenaltyMPC_Yd07, statePenaltyMPC_Ld07);
statePenaltyMPC_LA_DENSE_MVMSUB1_8_8(statePenaltyMPC_Lsd07, statePenaltyMPC_yy06, statePenaltyMPC_beta07, statePenaltyMPC_bmy07);
statePenaltyMPC_LA_DENSE_FORWARDSUB_8(statePenaltyMPC_Ld07, statePenaltyMPC_bmy07, statePenaltyMPC_yy07);
statePenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_8_8(statePenaltyMPC_Ld07, statePenaltyMPC_Ysd08, statePenaltyMPC_Lsd08);
statePenaltyMPC_LA_DENSE_MMTSUB_8_8(statePenaltyMPC_Lsd08, statePenaltyMPC_Yd08);
statePenaltyMPC_LA_DENSE_CHOL_8(statePenaltyMPC_Yd08, statePenaltyMPC_Ld08);
statePenaltyMPC_LA_DENSE_MVMSUB1_8_8(statePenaltyMPC_Lsd08, statePenaltyMPC_yy07, statePenaltyMPC_beta08, statePenaltyMPC_bmy08);
statePenaltyMPC_LA_DENSE_FORWARDSUB_8(statePenaltyMPC_Ld08, statePenaltyMPC_bmy08, statePenaltyMPC_yy08);
statePenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_8_8(statePenaltyMPC_Ld08, statePenaltyMPC_Ysd09, statePenaltyMPC_Lsd09);
statePenaltyMPC_LA_DENSE_MMTSUB_8_8(statePenaltyMPC_Lsd09, statePenaltyMPC_Yd09);
statePenaltyMPC_LA_DENSE_CHOL_8(statePenaltyMPC_Yd09, statePenaltyMPC_Ld09);
statePenaltyMPC_LA_DENSE_MVMSUB1_8_8(statePenaltyMPC_Lsd09, statePenaltyMPC_yy08, statePenaltyMPC_beta09, statePenaltyMPC_bmy09);
statePenaltyMPC_LA_DENSE_FORWARDSUB_8(statePenaltyMPC_Ld09, statePenaltyMPC_bmy09, statePenaltyMPC_yy09);
statePenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_8_8(statePenaltyMPC_Ld09, statePenaltyMPC_Ysd10, statePenaltyMPC_Lsd10);
statePenaltyMPC_LA_DENSE_MMTSUB_8_8(statePenaltyMPC_Lsd10, statePenaltyMPC_Yd10);
statePenaltyMPC_LA_DENSE_CHOL_8(statePenaltyMPC_Yd10, statePenaltyMPC_Ld10);
statePenaltyMPC_LA_DENSE_MVMSUB1_8_8(statePenaltyMPC_Lsd10, statePenaltyMPC_yy09, statePenaltyMPC_beta10, statePenaltyMPC_bmy10);
statePenaltyMPC_LA_DENSE_FORWARDSUB_8(statePenaltyMPC_Ld10, statePenaltyMPC_bmy10, statePenaltyMPC_yy10);
statePenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_8_8(statePenaltyMPC_Ld10, statePenaltyMPC_Ysd11, statePenaltyMPC_Lsd11);
statePenaltyMPC_LA_DENSE_MMTSUB_8_8(statePenaltyMPC_Lsd11, statePenaltyMPC_Yd11);
statePenaltyMPC_LA_DENSE_CHOL_8(statePenaltyMPC_Yd11, statePenaltyMPC_Ld11);
statePenaltyMPC_LA_DENSE_MVMSUB1_8_8(statePenaltyMPC_Lsd11, statePenaltyMPC_yy10, statePenaltyMPC_beta11, statePenaltyMPC_bmy11);
statePenaltyMPC_LA_DENSE_FORWARDSUB_8(statePenaltyMPC_Ld11, statePenaltyMPC_bmy11, statePenaltyMPC_yy11);
statePenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_8_8(statePenaltyMPC_Ld11, statePenaltyMPC_Ysd12, statePenaltyMPC_Lsd12);
statePenaltyMPC_LA_DENSE_MMTSUB_8_8(statePenaltyMPC_Lsd12, statePenaltyMPC_Yd12);
statePenaltyMPC_LA_DENSE_CHOL_8(statePenaltyMPC_Yd12, statePenaltyMPC_Ld12);
statePenaltyMPC_LA_DENSE_MVMSUB1_8_8(statePenaltyMPC_Lsd12, statePenaltyMPC_yy11, statePenaltyMPC_beta12, statePenaltyMPC_bmy12);
statePenaltyMPC_LA_DENSE_FORWARDSUB_8(statePenaltyMPC_Ld12, statePenaltyMPC_bmy12, statePenaltyMPC_yy12);
statePenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_8_8(statePenaltyMPC_Ld12, statePenaltyMPC_Ysd13, statePenaltyMPC_Lsd13);
statePenaltyMPC_LA_DENSE_MMTSUB_8_8(statePenaltyMPC_Lsd13, statePenaltyMPC_Yd13);
statePenaltyMPC_LA_DENSE_CHOL_8(statePenaltyMPC_Yd13, statePenaltyMPC_Ld13);
statePenaltyMPC_LA_DENSE_MVMSUB1_8_8(statePenaltyMPC_Lsd13, statePenaltyMPC_yy12, statePenaltyMPC_beta13, statePenaltyMPC_bmy13);
statePenaltyMPC_LA_DENSE_FORWARDSUB_8(statePenaltyMPC_Ld13, statePenaltyMPC_bmy13, statePenaltyMPC_yy13);
statePenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_8_8(statePenaltyMPC_Ld13, statePenaltyMPC_Ysd14, statePenaltyMPC_Lsd14);
statePenaltyMPC_LA_DENSE_MMTSUB_8_8(statePenaltyMPC_Lsd14, statePenaltyMPC_Yd14);
statePenaltyMPC_LA_DENSE_CHOL_8(statePenaltyMPC_Yd14, statePenaltyMPC_Ld14);
statePenaltyMPC_LA_DENSE_MVMSUB1_8_8(statePenaltyMPC_Lsd14, statePenaltyMPC_yy13, statePenaltyMPC_beta14, statePenaltyMPC_bmy14);
statePenaltyMPC_LA_DENSE_FORWARDSUB_8(statePenaltyMPC_Ld14, statePenaltyMPC_bmy14, statePenaltyMPC_yy14);
statePenaltyMPC_LA_DENSE_BACKWARDSUB_8(statePenaltyMPC_Ld14, statePenaltyMPC_yy14, statePenaltyMPC_dvaff14);
statePenaltyMPC_LA_DENSE_MTVMSUB_8_8(statePenaltyMPC_Lsd14, statePenaltyMPC_dvaff14, statePenaltyMPC_yy13, statePenaltyMPC_bmy13);
statePenaltyMPC_LA_DENSE_BACKWARDSUB_8(statePenaltyMPC_Ld13, statePenaltyMPC_bmy13, statePenaltyMPC_dvaff13);
statePenaltyMPC_LA_DENSE_MTVMSUB_8_8(statePenaltyMPC_Lsd13, statePenaltyMPC_dvaff13, statePenaltyMPC_yy12, statePenaltyMPC_bmy12);
statePenaltyMPC_LA_DENSE_BACKWARDSUB_8(statePenaltyMPC_Ld12, statePenaltyMPC_bmy12, statePenaltyMPC_dvaff12);
statePenaltyMPC_LA_DENSE_MTVMSUB_8_8(statePenaltyMPC_Lsd12, statePenaltyMPC_dvaff12, statePenaltyMPC_yy11, statePenaltyMPC_bmy11);
statePenaltyMPC_LA_DENSE_BACKWARDSUB_8(statePenaltyMPC_Ld11, statePenaltyMPC_bmy11, statePenaltyMPC_dvaff11);
statePenaltyMPC_LA_DENSE_MTVMSUB_8_8(statePenaltyMPC_Lsd11, statePenaltyMPC_dvaff11, statePenaltyMPC_yy10, statePenaltyMPC_bmy10);
statePenaltyMPC_LA_DENSE_BACKWARDSUB_8(statePenaltyMPC_Ld10, statePenaltyMPC_bmy10, statePenaltyMPC_dvaff10);
statePenaltyMPC_LA_DENSE_MTVMSUB_8_8(statePenaltyMPC_Lsd10, statePenaltyMPC_dvaff10, statePenaltyMPC_yy09, statePenaltyMPC_bmy09);
statePenaltyMPC_LA_DENSE_BACKWARDSUB_8(statePenaltyMPC_Ld09, statePenaltyMPC_bmy09, statePenaltyMPC_dvaff09);
statePenaltyMPC_LA_DENSE_MTVMSUB_8_8(statePenaltyMPC_Lsd09, statePenaltyMPC_dvaff09, statePenaltyMPC_yy08, statePenaltyMPC_bmy08);
statePenaltyMPC_LA_DENSE_BACKWARDSUB_8(statePenaltyMPC_Ld08, statePenaltyMPC_bmy08, statePenaltyMPC_dvaff08);
statePenaltyMPC_LA_DENSE_MTVMSUB_8_8(statePenaltyMPC_Lsd08, statePenaltyMPC_dvaff08, statePenaltyMPC_yy07, statePenaltyMPC_bmy07);
statePenaltyMPC_LA_DENSE_BACKWARDSUB_8(statePenaltyMPC_Ld07, statePenaltyMPC_bmy07, statePenaltyMPC_dvaff07);
statePenaltyMPC_LA_DENSE_MTVMSUB_8_8(statePenaltyMPC_Lsd07, statePenaltyMPC_dvaff07, statePenaltyMPC_yy06, statePenaltyMPC_bmy06);
statePenaltyMPC_LA_DENSE_BACKWARDSUB_8(statePenaltyMPC_Ld06, statePenaltyMPC_bmy06, statePenaltyMPC_dvaff06);
statePenaltyMPC_LA_DENSE_MTVMSUB_8_8(statePenaltyMPC_Lsd06, statePenaltyMPC_dvaff06, statePenaltyMPC_yy05, statePenaltyMPC_bmy05);
statePenaltyMPC_LA_DENSE_BACKWARDSUB_8(statePenaltyMPC_Ld05, statePenaltyMPC_bmy05, statePenaltyMPC_dvaff05);
statePenaltyMPC_LA_DENSE_MTVMSUB_8_8(statePenaltyMPC_Lsd05, statePenaltyMPC_dvaff05, statePenaltyMPC_yy04, statePenaltyMPC_bmy04);
statePenaltyMPC_LA_DENSE_BACKWARDSUB_8(statePenaltyMPC_Ld04, statePenaltyMPC_bmy04, statePenaltyMPC_dvaff04);
statePenaltyMPC_LA_DENSE_MTVMSUB_8_8(statePenaltyMPC_Lsd04, statePenaltyMPC_dvaff04, statePenaltyMPC_yy03, statePenaltyMPC_bmy03);
statePenaltyMPC_LA_DENSE_BACKWARDSUB_8(statePenaltyMPC_Ld03, statePenaltyMPC_bmy03, statePenaltyMPC_dvaff03);
statePenaltyMPC_LA_DENSE_MTVMSUB_8_8(statePenaltyMPC_Lsd03, statePenaltyMPC_dvaff03, statePenaltyMPC_yy02, statePenaltyMPC_bmy02);
statePenaltyMPC_LA_DENSE_BACKWARDSUB_8(statePenaltyMPC_Ld02, statePenaltyMPC_bmy02, statePenaltyMPC_dvaff02);
statePenaltyMPC_LA_DENSE_MTVMSUB_8_8(statePenaltyMPC_Lsd02, statePenaltyMPC_dvaff02, statePenaltyMPC_yy01, statePenaltyMPC_bmy01);
statePenaltyMPC_LA_DENSE_BACKWARDSUB_8(statePenaltyMPC_Ld01, statePenaltyMPC_bmy01, statePenaltyMPC_dvaff01);
statePenaltyMPC_LA_DENSE_MTVMSUB_8_8(statePenaltyMPC_Lsd01, statePenaltyMPC_dvaff01, statePenaltyMPC_yy00, statePenaltyMPC_bmy00);
statePenaltyMPC_LA_DENSE_BACKWARDSUB_8(statePenaltyMPC_Ld00, statePenaltyMPC_bmy00, statePenaltyMPC_dvaff00);
statePenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_8_26_8(params->C1, statePenaltyMPC_dvaff01, statePenaltyMPC_D00, statePenaltyMPC_dvaff00, statePenaltyMPC_grad_eq00);
statePenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_8_26_8(params->C2, statePenaltyMPC_dvaff02, statePenaltyMPC_D01, statePenaltyMPC_dvaff01, statePenaltyMPC_grad_eq01);
statePenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_8_26_8(params->C3, statePenaltyMPC_dvaff03, statePenaltyMPC_D01, statePenaltyMPC_dvaff02, statePenaltyMPC_grad_eq02);
statePenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_8_26_8(params->C4, statePenaltyMPC_dvaff04, statePenaltyMPC_D01, statePenaltyMPC_dvaff03, statePenaltyMPC_grad_eq03);
statePenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_8_26_8(params->C5, statePenaltyMPC_dvaff05, statePenaltyMPC_D01, statePenaltyMPC_dvaff04, statePenaltyMPC_grad_eq04);
statePenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_8_26_8(params->C6, statePenaltyMPC_dvaff06, statePenaltyMPC_D01, statePenaltyMPC_dvaff05, statePenaltyMPC_grad_eq05);
statePenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_8_26_8(params->C7, statePenaltyMPC_dvaff07, statePenaltyMPC_D01, statePenaltyMPC_dvaff06, statePenaltyMPC_grad_eq06);
statePenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_8_26_8(params->C8, statePenaltyMPC_dvaff08, statePenaltyMPC_D01, statePenaltyMPC_dvaff07, statePenaltyMPC_grad_eq07);
statePenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_8_26_8(params->C9, statePenaltyMPC_dvaff09, statePenaltyMPC_D01, statePenaltyMPC_dvaff08, statePenaltyMPC_grad_eq08);
statePenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_8_26_8(params->C10, statePenaltyMPC_dvaff10, statePenaltyMPC_D01, statePenaltyMPC_dvaff09, statePenaltyMPC_grad_eq09);
statePenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_8_26_8(params->C11, statePenaltyMPC_dvaff11, statePenaltyMPC_D01, statePenaltyMPC_dvaff10, statePenaltyMPC_grad_eq10);
statePenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_8_26_8(params->C12, statePenaltyMPC_dvaff12, statePenaltyMPC_D01, statePenaltyMPC_dvaff11, statePenaltyMPC_grad_eq11);
statePenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_8_26_8(params->C13, statePenaltyMPC_dvaff13, statePenaltyMPC_D01, statePenaltyMPC_dvaff12, statePenaltyMPC_grad_eq12);
statePenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_8_26_8(params->C14, statePenaltyMPC_dvaff14, statePenaltyMPC_D01, statePenaltyMPC_dvaff13, statePenaltyMPC_grad_eq13);
statePenaltyMPC_LA_DIAGZERO_MTVM_8_8(statePenaltyMPC_D14, statePenaltyMPC_dvaff14, statePenaltyMPC_grad_eq14);
statePenaltyMPC_LA_VSUB2_372(statePenaltyMPC_rd, statePenaltyMPC_grad_eq, statePenaltyMPC_rd);
statePenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_26(statePenaltyMPC_Phi00, statePenaltyMPC_rd00, statePenaltyMPC_dzaff00);
statePenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_26(statePenaltyMPC_Phi01, statePenaltyMPC_rd01, statePenaltyMPC_dzaff01);
statePenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_26(statePenaltyMPC_Phi02, statePenaltyMPC_rd02, statePenaltyMPC_dzaff02);
statePenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_26(statePenaltyMPC_Phi03, statePenaltyMPC_rd03, statePenaltyMPC_dzaff03);
statePenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_26(statePenaltyMPC_Phi04, statePenaltyMPC_rd04, statePenaltyMPC_dzaff04);
statePenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_26(statePenaltyMPC_Phi05, statePenaltyMPC_rd05, statePenaltyMPC_dzaff05);
statePenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_26(statePenaltyMPC_Phi06, statePenaltyMPC_rd06, statePenaltyMPC_dzaff06);
statePenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_26(statePenaltyMPC_Phi07, statePenaltyMPC_rd07, statePenaltyMPC_dzaff07);
statePenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_26(statePenaltyMPC_Phi08, statePenaltyMPC_rd08, statePenaltyMPC_dzaff08);
statePenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_26(statePenaltyMPC_Phi09, statePenaltyMPC_rd09, statePenaltyMPC_dzaff09);
statePenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_26(statePenaltyMPC_Phi10, statePenaltyMPC_rd10, statePenaltyMPC_dzaff10);
statePenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_26(statePenaltyMPC_Phi11, statePenaltyMPC_rd11, statePenaltyMPC_dzaff11);
statePenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_26(statePenaltyMPC_Phi12, statePenaltyMPC_rd12, statePenaltyMPC_dzaff12);
statePenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_26(statePenaltyMPC_Phi13, statePenaltyMPC_rd13, statePenaltyMPC_dzaff13);
statePenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_8(statePenaltyMPC_Phi14, statePenaltyMPC_rd14, statePenaltyMPC_dzaff14);
statePenaltyMPC_LA_VSUB_INDEXED_26(statePenaltyMPC_dzaff00, statePenaltyMPC_lbIdx00, statePenaltyMPC_rilb00, statePenaltyMPC_dslbaff00);
statePenaltyMPC_LA_VSUB3_26(statePenaltyMPC_llbbyslb00, statePenaltyMPC_dslbaff00, statePenaltyMPC_llb00, statePenaltyMPC_dllbaff00);
statePenaltyMPC_LA_VSUB2_INDEXED_10(statePenaltyMPC_riub00, statePenaltyMPC_dzaff00, statePenaltyMPC_ubIdx00, statePenaltyMPC_dsubaff00);
statePenaltyMPC_LA_VSUB3_10(statePenaltyMPC_lubbysub00, statePenaltyMPC_dsubaff00, statePenaltyMPC_lub00, statePenaltyMPC_dlubaff00);
statePenaltyMPC_LA_VSUB_INDEXED_26(statePenaltyMPC_dzaff01, statePenaltyMPC_lbIdx01, statePenaltyMPC_rilb01, statePenaltyMPC_dslbaff01);
statePenaltyMPC_LA_VSUB3_26(statePenaltyMPC_llbbyslb01, statePenaltyMPC_dslbaff01, statePenaltyMPC_llb01, statePenaltyMPC_dllbaff01);
statePenaltyMPC_LA_VSUB2_INDEXED_10(statePenaltyMPC_riub01, statePenaltyMPC_dzaff01, statePenaltyMPC_ubIdx01, statePenaltyMPC_dsubaff01);
statePenaltyMPC_LA_VSUB3_10(statePenaltyMPC_lubbysub01, statePenaltyMPC_dsubaff01, statePenaltyMPC_lub01, statePenaltyMPC_dlubaff01);
statePenaltyMPC_LA_VSUB_INDEXED_26(statePenaltyMPC_dzaff02, statePenaltyMPC_lbIdx02, statePenaltyMPC_rilb02, statePenaltyMPC_dslbaff02);
statePenaltyMPC_LA_VSUB3_26(statePenaltyMPC_llbbyslb02, statePenaltyMPC_dslbaff02, statePenaltyMPC_llb02, statePenaltyMPC_dllbaff02);
statePenaltyMPC_LA_VSUB2_INDEXED_10(statePenaltyMPC_riub02, statePenaltyMPC_dzaff02, statePenaltyMPC_ubIdx02, statePenaltyMPC_dsubaff02);
statePenaltyMPC_LA_VSUB3_10(statePenaltyMPC_lubbysub02, statePenaltyMPC_dsubaff02, statePenaltyMPC_lub02, statePenaltyMPC_dlubaff02);
statePenaltyMPC_LA_VSUB_INDEXED_26(statePenaltyMPC_dzaff03, statePenaltyMPC_lbIdx03, statePenaltyMPC_rilb03, statePenaltyMPC_dslbaff03);
statePenaltyMPC_LA_VSUB3_26(statePenaltyMPC_llbbyslb03, statePenaltyMPC_dslbaff03, statePenaltyMPC_llb03, statePenaltyMPC_dllbaff03);
statePenaltyMPC_LA_VSUB2_INDEXED_10(statePenaltyMPC_riub03, statePenaltyMPC_dzaff03, statePenaltyMPC_ubIdx03, statePenaltyMPC_dsubaff03);
statePenaltyMPC_LA_VSUB3_10(statePenaltyMPC_lubbysub03, statePenaltyMPC_dsubaff03, statePenaltyMPC_lub03, statePenaltyMPC_dlubaff03);
statePenaltyMPC_LA_VSUB_INDEXED_26(statePenaltyMPC_dzaff04, statePenaltyMPC_lbIdx04, statePenaltyMPC_rilb04, statePenaltyMPC_dslbaff04);
statePenaltyMPC_LA_VSUB3_26(statePenaltyMPC_llbbyslb04, statePenaltyMPC_dslbaff04, statePenaltyMPC_llb04, statePenaltyMPC_dllbaff04);
statePenaltyMPC_LA_VSUB2_INDEXED_10(statePenaltyMPC_riub04, statePenaltyMPC_dzaff04, statePenaltyMPC_ubIdx04, statePenaltyMPC_dsubaff04);
statePenaltyMPC_LA_VSUB3_10(statePenaltyMPC_lubbysub04, statePenaltyMPC_dsubaff04, statePenaltyMPC_lub04, statePenaltyMPC_dlubaff04);
statePenaltyMPC_LA_VSUB_INDEXED_26(statePenaltyMPC_dzaff05, statePenaltyMPC_lbIdx05, statePenaltyMPC_rilb05, statePenaltyMPC_dslbaff05);
statePenaltyMPC_LA_VSUB3_26(statePenaltyMPC_llbbyslb05, statePenaltyMPC_dslbaff05, statePenaltyMPC_llb05, statePenaltyMPC_dllbaff05);
statePenaltyMPC_LA_VSUB2_INDEXED_10(statePenaltyMPC_riub05, statePenaltyMPC_dzaff05, statePenaltyMPC_ubIdx05, statePenaltyMPC_dsubaff05);
statePenaltyMPC_LA_VSUB3_10(statePenaltyMPC_lubbysub05, statePenaltyMPC_dsubaff05, statePenaltyMPC_lub05, statePenaltyMPC_dlubaff05);
statePenaltyMPC_LA_VSUB_INDEXED_26(statePenaltyMPC_dzaff06, statePenaltyMPC_lbIdx06, statePenaltyMPC_rilb06, statePenaltyMPC_dslbaff06);
statePenaltyMPC_LA_VSUB3_26(statePenaltyMPC_llbbyslb06, statePenaltyMPC_dslbaff06, statePenaltyMPC_llb06, statePenaltyMPC_dllbaff06);
statePenaltyMPC_LA_VSUB2_INDEXED_10(statePenaltyMPC_riub06, statePenaltyMPC_dzaff06, statePenaltyMPC_ubIdx06, statePenaltyMPC_dsubaff06);
statePenaltyMPC_LA_VSUB3_10(statePenaltyMPC_lubbysub06, statePenaltyMPC_dsubaff06, statePenaltyMPC_lub06, statePenaltyMPC_dlubaff06);
statePenaltyMPC_LA_VSUB_INDEXED_26(statePenaltyMPC_dzaff07, statePenaltyMPC_lbIdx07, statePenaltyMPC_rilb07, statePenaltyMPC_dslbaff07);
statePenaltyMPC_LA_VSUB3_26(statePenaltyMPC_llbbyslb07, statePenaltyMPC_dslbaff07, statePenaltyMPC_llb07, statePenaltyMPC_dllbaff07);
statePenaltyMPC_LA_VSUB2_INDEXED_10(statePenaltyMPC_riub07, statePenaltyMPC_dzaff07, statePenaltyMPC_ubIdx07, statePenaltyMPC_dsubaff07);
statePenaltyMPC_LA_VSUB3_10(statePenaltyMPC_lubbysub07, statePenaltyMPC_dsubaff07, statePenaltyMPC_lub07, statePenaltyMPC_dlubaff07);
statePenaltyMPC_LA_VSUB_INDEXED_26(statePenaltyMPC_dzaff08, statePenaltyMPC_lbIdx08, statePenaltyMPC_rilb08, statePenaltyMPC_dslbaff08);
statePenaltyMPC_LA_VSUB3_26(statePenaltyMPC_llbbyslb08, statePenaltyMPC_dslbaff08, statePenaltyMPC_llb08, statePenaltyMPC_dllbaff08);
statePenaltyMPC_LA_VSUB2_INDEXED_10(statePenaltyMPC_riub08, statePenaltyMPC_dzaff08, statePenaltyMPC_ubIdx08, statePenaltyMPC_dsubaff08);
statePenaltyMPC_LA_VSUB3_10(statePenaltyMPC_lubbysub08, statePenaltyMPC_dsubaff08, statePenaltyMPC_lub08, statePenaltyMPC_dlubaff08);
statePenaltyMPC_LA_VSUB_INDEXED_26(statePenaltyMPC_dzaff09, statePenaltyMPC_lbIdx09, statePenaltyMPC_rilb09, statePenaltyMPC_dslbaff09);
statePenaltyMPC_LA_VSUB3_26(statePenaltyMPC_llbbyslb09, statePenaltyMPC_dslbaff09, statePenaltyMPC_llb09, statePenaltyMPC_dllbaff09);
statePenaltyMPC_LA_VSUB2_INDEXED_10(statePenaltyMPC_riub09, statePenaltyMPC_dzaff09, statePenaltyMPC_ubIdx09, statePenaltyMPC_dsubaff09);
statePenaltyMPC_LA_VSUB3_10(statePenaltyMPC_lubbysub09, statePenaltyMPC_dsubaff09, statePenaltyMPC_lub09, statePenaltyMPC_dlubaff09);
statePenaltyMPC_LA_VSUB_INDEXED_26(statePenaltyMPC_dzaff10, statePenaltyMPC_lbIdx10, statePenaltyMPC_rilb10, statePenaltyMPC_dslbaff10);
statePenaltyMPC_LA_VSUB3_26(statePenaltyMPC_llbbyslb10, statePenaltyMPC_dslbaff10, statePenaltyMPC_llb10, statePenaltyMPC_dllbaff10);
statePenaltyMPC_LA_VSUB2_INDEXED_10(statePenaltyMPC_riub10, statePenaltyMPC_dzaff10, statePenaltyMPC_ubIdx10, statePenaltyMPC_dsubaff10);
statePenaltyMPC_LA_VSUB3_10(statePenaltyMPC_lubbysub10, statePenaltyMPC_dsubaff10, statePenaltyMPC_lub10, statePenaltyMPC_dlubaff10);
statePenaltyMPC_LA_VSUB_INDEXED_26(statePenaltyMPC_dzaff11, statePenaltyMPC_lbIdx11, statePenaltyMPC_rilb11, statePenaltyMPC_dslbaff11);
statePenaltyMPC_LA_VSUB3_26(statePenaltyMPC_llbbyslb11, statePenaltyMPC_dslbaff11, statePenaltyMPC_llb11, statePenaltyMPC_dllbaff11);
statePenaltyMPC_LA_VSUB2_INDEXED_10(statePenaltyMPC_riub11, statePenaltyMPC_dzaff11, statePenaltyMPC_ubIdx11, statePenaltyMPC_dsubaff11);
statePenaltyMPC_LA_VSUB3_10(statePenaltyMPC_lubbysub11, statePenaltyMPC_dsubaff11, statePenaltyMPC_lub11, statePenaltyMPC_dlubaff11);
statePenaltyMPC_LA_VSUB_INDEXED_26(statePenaltyMPC_dzaff12, statePenaltyMPC_lbIdx12, statePenaltyMPC_rilb12, statePenaltyMPC_dslbaff12);
statePenaltyMPC_LA_VSUB3_26(statePenaltyMPC_llbbyslb12, statePenaltyMPC_dslbaff12, statePenaltyMPC_llb12, statePenaltyMPC_dllbaff12);
statePenaltyMPC_LA_VSUB2_INDEXED_10(statePenaltyMPC_riub12, statePenaltyMPC_dzaff12, statePenaltyMPC_ubIdx12, statePenaltyMPC_dsubaff12);
statePenaltyMPC_LA_VSUB3_10(statePenaltyMPC_lubbysub12, statePenaltyMPC_dsubaff12, statePenaltyMPC_lub12, statePenaltyMPC_dlubaff12);
statePenaltyMPC_LA_VSUB_INDEXED_26(statePenaltyMPC_dzaff13, statePenaltyMPC_lbIdx13, statePenaltyMPC_rilb13, statePenaltyMPC_dslbaff13);
statePenaltyMPC_LA_VSUB3_26(statePenaltyMPC_llbbyslb13, statePenaltyMPC_dslbaff13, statePenaltyMPC_llb13, statePenaltyMPC_dllbaff13);
statePenaltyMPC_LA_VSUB2_INDEXED_10(statePenaltyMPC_riub13, statePenaltyMPC_dzaff13, statePenaltyMPC_ubIdx13, statePenaltyMPC_dsubaff13);
statePenaltyMPC_LA_VSUB3_10(statePenaltyMPC_lubbysub13, statePenaltyMPC_dsubaff13, statePenaltyMPC_lub13, statePenaltyMPC_dlubaff13);
statePenaltyMPC_LA_VSUB_INDEXED_8(statePenaltyMPC_dzaff14, statePenaltyMPC_lbIdx14, statePenaltyMPC_rilb14, statePenaltyMPC_dslbaff14);
statePenaltyMPC_LA_VSUB3_8(statePenaltyMPC_llbbyslb14, statePenaltyMPC_dslbaff14, statePenaltyMPC_llb14, statePenaltyMPC_dllbaff14);
statePenaltyMPC_LA_VSUB2_INDEXED_8(statePenaltyMPC_riub14, statePenaltyMPC_dzaff14, statePenaltyMPC_ubIdx14, statePenaltyMPC_dsubaff14);
statePenaltyMPC_LA_VSUB3_8(statePenaltyMPC_lubbysub14, statePenaltyMPC_dsubaff14, statePenaltyMPC_lub14, statePenaltyMPC_dlubaff14);
info->lsit_aff = statePenaltyMPC_LINESEARCH_BACKTRACKING_AFFINE(statePenaltyMPC_l, statePenaltyMPC_s, statePenaltyMPC_dl_aff, statePenaltyMPC_ds_aff, &info->step_aff, &info->mu_aff);
if( info->lsit_aff == statePenaltyMPC_NOPROGRESS ){
PRINTTEXT("Affine line search could not proceed at iteration %d.\nThe problem might be infeasible -- exiting.\n",info->it+1);
exitcode = statePenaltyMPC_NOPROGRESS; break;
}
sigma_3rdroot = info->mu_aff / info->mu;
info->sigma = sigma_3rdroot*sigma_3rdroot*sigma_3rdroot;
musigma = info->mu * info->sigma;
statePenaltyMPC_LA_VSUB5_520(statePenaltyMPC_ds_aff, statePenaltyMPC_dl_aff, musigma, statePenaltyMPC_ccrhs);
statePenaltyMPC_LA_VSUB6_INDEXED_26_10_26(statePenaltyMPC_ccrhsub00, statePenaltyMPC_sub00, statePenaltyMPC_ubIdx00, statePenaltyMPC_ccrhsl00, statePenaltyMPC_slb00, statePenaltyMPC_lbIdx00, statePenaltyMPC_rd00);
statePenaltyMPC_LA_VSUB6_INDEXED_26_10_26(statePenaltyMPC_ccrhsub01, statePenaltyMPC_sub01, statePenaltyMPC_ubIdx01, statePenaltyMPC_ccrhsl01, statePenaltyMPC_slb01, statePenaltyMPC_lbIdx01, statePenaltyMPC_rd01);
statePenaltyMPC_LA_DIAG_FORWARDSUB_26(statePenaltyMPC_Phi00, statePenaltyMPC_rd00, statePenaltyMPC_Lbyrd00);
statePenaltyMPC_LA_DIAG_FORWARDSUB_26(statePenaltyMPC_Phi01, statePenaltyMPC_rd01, statePenaltyMPC_Lbyrd01);
statePenaltyMPC_LA_DIAGZERO_MVM_8(statePenaltyMPC_W00, statePenaltyMPC_Lbyrd00, statePenaltyMPC_beta00);
statePenaltyMPC_LA_DENSE_FORWARDSUB_8(statePenaltyMPC_Ld00, statePenaltyMPC_beta00, statePenaltyMPC_yy00);
statePenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_8_26_26(statePenaltyMPC_V00, statePenaltyMPC_Lbyrd00, statePenaltyMPC_W01, statePenaltyMPC_Lbyrd01, statePenaltyMPC_beta01);
statePenaltyMPC_LA_DENSE_MVMSUB1_8_8(statePenaltyMPC_Lsd01, statePenaltyMPC_yy00, statePenaltyMPC_beta01, statePenaltyMPC_bmy01);
statePenaltyMPC_LA_DENSE_FORWARDSUB_8(statePenaltyMPC_Ld01, statePenaltyMPC_bmy01, statePenaltyMPC_yy01);
statePenaltyMPC_LA_VSUB6_INDEXED_26_10_26(statePenaltyMPC_ccrhsub02, statePenaltyMPC_sub02, statePenaltyMPC_ubIdx02, statePenaltyMPC_ccrhsl02, statePenaltyMPC_slb02, statePenaltyMPC_lbIdx02, statePenaltyMPC_rd02);
statePenaltyMPC_LA_DIAG_FORWARDSUB_26(statePenaltyMPC_Phi02, statePenaltyMPC_rd02, statePenaltyMPC_Lbyrd02);
statePenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_8_26_26(statePenaltyMPC_V01, statePenaltyMPC_Lbyrd01, statePenaltyMPC_W02, statePenaltyMPC_Lbyrd02, statePenaltyMPC_beta02);
statePenaltyMPC_LA_DENSE_MVMSUB1_8_8(statePenaltyMPC_Lsd02, statePenaltyMPC_yy01, statePenaltyMPC_beta02, statePenaltyMPC_bmy02);
statePenaltyMPC_LA_DENSE_FORWARDSUB_8(statePenaltyMPC_Ld02, statePenaltyMPC_bmy02, statePenaltyMPC_yy02);
statePenaltyMPC_LA_VSUB6_INDEXED_26_10_26(statePenaltyMPC_ccrhsub03, statePenaltyMPC_sub03, statePenaltyMPC_ubIdx03, statePenaltyMPC_ccrhsl03, statePenaltyMPC_slb03, statePenaltyMPC_lbIdx03, statePenaltyMPC_rd03);
statePenaltyMPC_LA_DIAG_FORWARDSUB_26(statePenaltyMPC_Phi03, statePenaltyMPC_rd03, statePenaltyMPC_Lbyrd03);
statePenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_8_26_26(statePenaltyMPC_V02, statePenaltyMPC_Lbyrd02, statePenaltyMPC_W03, statePenaltyMPC_Lbyrd03, statePenaltyMPC_beta03);
statePenaltyMPC_LA_DENSE_MVMSUB1_8_8(statePenaltyMPC_Lsd03, statePenaltyMPC_yy02, statePenaltyMPC_beta03, statePenaltyMPC_bmy03);
statePenaltyMPC_LA_DENSE_FORWARDSUB_8(statePenaltyMPC_Ld03, statePenaltyMPC_bmy03, statePenaltyMPC_yy03);
statePenaltyMPC_LA_VSUB6_INDEXED_26_10_26(statePenaltyMPC_ccrhsub04, statePenaltyMPC_sub04, statePenaltyMPC_ubIdx04, statePenaltyMPC_ccrhsl04, statePenaltyMPC_slb04, statePenaltyMPC_lbIdx04, statePenaltyMPC_rd04);
statePenaltyMPC_LA_DIAG_FORWARDSUB_26(statePenaltyMPC_Phi04, statePenaltyMPC_rd04, statePenaltyMPC_Lbyrd04);
statePenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_8_26_26(statePenaltyMPC_V03, statePenaltyMPC_Lbyrd03, statePenaltyMPC_W04, statePenaltyMPC_Lbyrd04, statePenaltyMPC_beta04);
statePenaltyMPC_LA_DENSE_MVMSUB1_8_8(statePenaltyMPC_Lsd04, statePenaltyMPC_yy03, statePenaltyMPC_beta04, statePenaltyMPC_bmy04);
statePenaltyMPC_LA_DENSE_FORWARDSUB_8(statePenaltyMPC_Ld04, statePenaltyMPC_bmy04, statePenaltyMPC_yy04);
statePenaltyMPC_LA_VSUB6_INDEXED_26_10_26(statePenaltyMPC_ccrhsub05, statePenaltyMPC_sub05, statePenaltyMPC_ubIdx05, statePenaltyMPC_ccrhsl05, statePenaltyMPC_slb05, statePenaltyMPC_lbIdx05, statePenaltyMPC_rd05);
statePenaltyMPC_LA_DIAG_FORWARDSUB_26(statePenaltyMPC_Phi05, statePenaltyMPC_rd05, statePenaltyMPC_Lbyrd05);
statePenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_8_26_26(statePenaltyMPC_V04, statePenaltyMPC_Lbyrd04, statePenaltyMPC_W05, statePenaltyMPC_Lbyrd05, statePenaltyMPC_beta05);
statePenaltyMPC_LA_DENSE_MVMSUB1_8_8(statePenaltyMPC_Lsd05, statePenaltyMPC_yy04, statePenaltyMPC_beta05, statePenaltyMPC_bmy05);
statePenaltyMPC_LA_DENSE_FORWARDSUB_8(statePenaltyMPC_Ld05, statePenaltyMPC_bmy05, statePenaltyMPC_yy05);
statePenaltyMPC_LA_VSUB6_INDEXED_26_10_26(statePenaltyMPC_ccrhsub06, statePenaltyMPC_sub06, statePenaltyMPC_ubIdx06, statePenaltyMPC_ccrhsl06, statePenaltyMPC_slb06, statePenaltyMPC_lbIdx06, statePenaltyMPC_rd06);
statePenaltyMPC_LA_DIAG_FORWARDSUB_26(statePenaltyMPC_Phi06, statePenaltyMPC_rd06, statePenaltyMPC_Lbyrd06);
statePenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_8_26_26(statePenaltyMPC_V05, statePenaltyMPC_Lbyrd05, statePenaltyMPC_W06, statePenaltyMPC_Lbyrd06, statePenaltyMPC_beta06);
statePenaltyMPC_LA_DENSE_MVMSUB1_8_8(statePenaltyMPC_Lsd06, statePenaltyMPC_yy05, statePenaltyMPC_beta06, statePenaltyMPC_bmy06);
statePenaltyMPC_LA_DENSE_FORWARDSUB_8(statePenaltyMPC_Ld06, statePenaltyMPC_bmy06, statePenaltyMPC_yy06);
statePenaltyMPC_LA_VSUB6_INDEXED_26_10_26(statePenaltyMPC_ccrhsub07, statePenaltyMPC_sub07, statePenaltyMPC_ubIdx07, statePenaltyMPC_ccrhsl07, statePenaltyMPC_slb07, statePenaltyMPC_lbIdx07, statePenaltyMPC_rd07);
statePenaltyMPC_LA_DIAG_FORWARDSUB_26(statePenaltyMPC_Phi07, statePenaltyMPC_rd07, statePenaltyMPC_Lbyrd07);
statePenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_8_26_26(statePenaltyMPC_V06, statePenaltyMPC_Lbyrd06, statePenaltyMPC_W07, statePenaltyMPC_Lbyrd07, statePenaltyMPC_beta07);
statePenaltyMPC_LA_DENSE_MVMSUB1_8_8(statePenaltyMPC_Lsd07, statePenaltyMPC_yy06, statePenaltyMPC_beta07, statePenaltyMPC_bmy07);
statePenaltyMPC_LA_DENSE_FORWARDSUB_8(statePenaltyMPC_Ld07, statePenaltyMPC_bmy07, statePenaltyMPC_yy07);
statePenaltyMPC_LA_VSUB6_INDEXED_26_10_26(statePenaltyMPC_ccrhsub08, statePenaltyMPC_sub08, statePenaltyMPC_ubIdx08, statePenaltyMPC_ccrhsl08, statePenaltyMPC_slb08, statePenaltyMPC_lbIdx08, statePenaltyMPC_rd08);
statePenaltyMPC_LA_DIAG_FORWARDSUB_26(statePenaltyMPC_Phi08, statePenaltyMPC_rd08, statePenaltyMPC_Lbyrd08);
statePenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_8_26_26(statePenaltyMPC_V07, statePenaltyMPC_Lbyrd07, statePenaltyMPC_W08, statePenaltyMPC_Lbyrd08, statePenaltyMPC_beta08);
statePenaltyMPC_LA_DENSE_MVMSUB1_8_8(statePenaltyMPC_Lsd08, statePenaltyMPC_yy07, statePenaltyMPC_beta08, statePenaltyMPC_bmy08);
statePenaltyMPC_LA_DENSE_FORWARDSUB_8(statePenaltyMPC_Ld08, statePenaltyMPC_bmy08, statePenaltyMPC_yy08);
statePenaltyMPC_LA_VSUB6_INDEXED_26_10_26(statePenaltyMPC_ccrhsub09, statePenaltyMPC_sub09, statePenaltyMPC_ubIdx09, statePenaltyMPC_ccrhsl09, statePenaltyMPC_slb09, statePenaltyMPC_lbIdx09, statePenaltyMPC_rd09);
statePenaltyMPC_LA_DIAG_FORWARDSUB_26(statePenaltyMPC_Phi09, statePenaltyMPC_rd09, statePenaltyMPC_Lbyrd09);
statePenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_8_26_26(statePenaltyMPC_V08, statePenaltyMPC_Lbyrd08, statePenaltyMPC_W09, statePenaltyMPC_Lbyrd09, statePenaltyMPC_beta09);
statePenaltyMPC_LA_DENSE_MVMSUB1_8_8(statePenaltyMPC_Lsd09, statePenaltyMPC_yy08, statePenaltyMPC_beta09, statePenaltyMPC_bmy09);
statePenaltyMPC_LA_DENSE_FORWARDSUB_8(statePenaltyMPC_Ld09, statePenaltyMPC_bmy09, statePenaltyMPC_yy09);
statePenaltyMPC_LA_VSUB6_INDEXED_26_10_26(statePenaltyMPC_ccrhsub10, statePenaltyMPC_sub10, statePenaltyMPC_ubIdx10, statePenaltyMPC_ccrhsl10, statePenaltyMPC_slb10, statePenaltyMPC_lbIdx10, statePenaltyMPC_rd10);
statePenaltyMPC_LA_DIAG_FORWARDSUB_26(statePenaltyMPC_Phi10, statePenaltyMPC_rd10, statePenaltyMPC_Lbyrd10);
statePenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_8_26_26(statePenaltyMPC_V09, statePenaltyMPC_Lbyrd09, statePenaltyMPC_W10, statePenaltyMPC_Lbyrd10, statePenaltyMPC_beta10);
statePenaltyMPC_LA_DENSE_MVMSUB1_8_8(statePenaltyMPC_Lsd10, statePenaltyMPC_yy09, statePenaltyMPC_beta10, statePenaltyMPC_bmy10);
statePenaltyMPC_LA_DENSE_FORWARDSUB_8(statePenaltyMPC_Ld10, statePenaltyMPC_bmy10, statePenaltyMPC_yy10);
statePenaltyMPC_LA_VSUB6_INDEXED_26_10_26(statePenaltyMPC_ccrhsub11, statePenaltyMPC_sub11, statePenaltyMPC_ubIdx11, statePenaltyMPC_ccrhsl11, statePenaltyMPC_slb11, statePenaltyMPC_lbIdx11, statePenaltyMPC_rd11);
statePenaltyMPC_LA_DIAG_FORWARDSUB_26(statePenaltyMPC_Phi11, statePenaltyMPC_rd11, statePenaltyMPC_Lbyrd11);
statePenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_8_26_26(statePenaltyMPC_V10, statePenaltyMPC_Lbyrd10, statePenaltyMPC_W11, statePenaltyMPC_Lbyrd11, statePenaltyMPC_beta11);
statePenaltyMPC_LA_DENSE_MVMSUB1_8_8(statePenaltyMPC_Lsd11, statePenaltyMPC_yy10, statePenaltyMPC_beta11, statePenaltyMPC_bmy11);
statePenaltyMPC_LA_DENSE_FORWARDSUB_8(statePenaltyMPC_Ld11, statePenaltyMPC_bmy11, statePenaltyMPC_yy11);
statePenaltyMPC_LA_VSUB6_INDEXED_26_10_26(statePenaltyMPC_ccrhsub12, statePenaltyMPC_sub12, statePenaltyMPC_ubIdx12, statePenaltyMPC_ccrhsl12, statePenaltyMPC_slb12, statePenaltyMPC_lbIdx12, statePenaltyMPC_rd12);
statePenaltyMPC_LA_DIAG_FORWARDSUB_26(statePenaltyMPC_Phi12, statePenaltyMPC_rd12, statePenaltyMPC_Lbyrd12);
statePenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_8_26_26(statePenaltyMPC_V11, statePenaltyMPC_Lbyrd11, statePenaltyMPC_W12, statePenaltyMPC_Lbyrd12, statePenaltyMPC_beta12);
statePenaltyMPC_LA_DENSE_MVMSUB1_8_8(statePenaltyMPC_Lsd12, statePenaltyMPC_yy11, statePenaltyMPC_beta12, statePenaltyMPC_bmy12);
statePenaltyMPC_LA_DENSE_FORWARDSUB_8(statePenaltyMPC_Ld12, statePenaltyMPC_bmy12, statePenaltyMPC_yy12);
statePenaltyMPC_LA_VSUB6_INDEXED_26_10_26(statePenaltyMPC_ccrhsub13, statePenaltyMPC_sub13, statePenaltyMPC_ubIdx13, statePenaltyMPC_ccrhsl13, statePenaltyMPC_slb13, statePenaltyMPC_lbIdx13, statePenaltyMPC_rd13);
statePenaltyMPC_LA_DIAG_FORWARDSUB_26(statePenaltyMPC_Phi13, statePenaltyMPC_rd13, statePenaltyMPC_Lbyrd13);
statePenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_8_26_26(statePenaltyMPC_V12, statePenaltyMPC_Lbyrd12, statePenaltyMPC_W13, statePenaltyMPC_Lbyrd13, statePenaltyMPC_beta13);
statePenaltyMPC_LA_DENSE_MVMSUB1_8_8(statePenaltyMPC_Lsd13, statePenaltyMPC_yy12, statePenaltyMPC_beta13, statePenaltyMPC_bmy13);
statePenaltyMPC_LA_DENSE_FORWARDSUB_8(statePenaltyMPC_Ld13, statePenaltyMPC_bmy13, statePenaltyMPC_yy13);
statePenaltyMPC_LA_VSUB6_INDEXED_8_8_8(statePenaltyMPC_ccrhsub14, statePenaltyMPC_sub14, statePenaltyMPC_ubIdx14, statePenaltyMPC_ccrhsl14, statePenaltyMPC_slb14, statePenaltyMPC_lbIdx14, statePenaltyMPC_rd14);
statePenaltyMPC_LA_DIAG_FORWARDSUB_8(statePenaltyMPC_Phi14, statePenaltyMPC_rd14, statePenaltyMPC_Lbyrd14);
statePenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_8_26_8(statePenaltyMPC_V13, statePenaltyMPC_Lbyrd13, statePenaltyMPC_W14, statePenaltyMPC_Lbyrd14, statePenaltyMPC_beta14);
statePenaltyMPC_LA_DENSE_MVMSUB1_8_8(statePenaltyMPC_Lsd14, statePenaltyMPC_yy13, statePenaltyMPC_beta14, statePenaltyMPC_bmy14);
statePenaltyMPC_LA_DENSE_FORWARDSUB_8(statePenaltyMPC_Ld14, statePenaltyMPC_bmy14, statePenaltyMPC_yy14);
statePenaltyMPC_LA_DENSE_BACKWARDSUB_8(statePenaltyMPC_Ld14, statePenaltyMPC_yy14, statePenaltyMPC_dvcc14);
statePenaltyMPC_LA_DENSE_MTVMSUB_8_8(statePenaltyMPC_Lsd14, statePenaltyMPC_dvcc14, statePenaltyMPC_yy13, statePenaltyMPC_bmy13);
statePenaltyMPC_LA_DENSE_BACKWARDSUB_8(statePenaltyMPC_Ld13, statePenaltyMPC_bmy13, statePenaltyMPC_dvcc13);
statePenaltyMPC_LA_DENSE_MTVMSUB_8_8(statePenaltyMPC_Lsd13, statePenaltyMPC_dvcc13, statePenaltyMPC_yy12, statePenaltyMPC_bmy12);
statePenaltyMPC_LA_DENSE_BACKWARDSUB_8(statePenaltyMPC_Ld12, statePenaltyMPC_bmy12, statePenaltyMPC_dvcc12);
statePenaltyMPC_LA_DENSE_MTVMSUB_8_8(statePenaltyMPC_Lsd12, statePenaltyMPC_dvcc12, statePenaltyMPC_yy11, statePenaltyMPC_bmy11);
statePenaltyMPC_LA_DENSE_BACKWARDSUB_8(statePenaltyMPC_Ld11, statePenaltyMPC_bmy11, statePenaltyMPC_dvcc11);
statePenaltyMPC_LA_DENSE_MTVMSUB_8_8(statePenaltyMPC_Lsd11, statePenaltyMPC_dvcc11, statePenaltyMPC_yy10, statePenaltyMPC_bmy10);
statePenaltyMPC_LA_DENSE_BACKWARDSUB_8(statePenaltyMPC_Ld10, statePenaltyMPC_bmy10, statePenaltyMPC_dvcc10);
statePenaltyMPC_LA_DENSE_MTVMSUB_8_8(statePenaltyMPC_Lsd10, statePenaltyMPC_dvcc10, statePenaltyMPC_yy09, statePenaltyMPC_bmy09);
statePenaltyMPC_LA_DENSE_BACKWARDSUB_8(statePenaltyMPC_Ld09, statePenaltyMPC_bmy09, statePenaltyMPC_dvcc09);
statePenaltyMPC_LA_DENSE_MTVMSUB_8_8(statePenaltyMPC_Lsd09, statePenaltyMPC_dvcc09, statePenaltyMPC_yy08, statePenaltyMPC_bmy08);
statePenaltyMPC_LA_DENSE_BACKWARDSUB_8(statePenaltyMPC_Ld08, statePenaltyMPC_bmy08, statePenaltyMPC_dvcc08);
statePenaltyMPC_LA_DENSE_MTVMSUB_8_8(statePenaltyMPC_Lsd08, statePenaltyMPC_dvcc08, statePenaltyMPC_yy07, statePenaltyMPC_bmy07);
statePenaltyMPC_LA_DENSE_BACKWARDSUB_8(statePenaltyMPC_Ld07, statePenaltyMPC_bmy07, statePenaltyMPC_dvcc07);
statePenaltyMPC_LA_DENSE_MTVMSUB_8_8(statePenaltyMPC_Lsd07, statePenaltyMPC_dvcc07, statePenaltyMPC_yy06, statePenaltyMPC_bmy06);
statePenaltyMPC_LA_DENSE_BACKWARDSUB_8(statePenaltyMPC_Ld06, statePenaltyMPC_bmy06, statePenaltyMPC_dvcc06);
statePenaltyMPC_LA_DENSE_MTVMSUB_8_8(statePenaltyMPC_Lsd06, statePenaltyMPC_dvcc06, statePenaltyMPC_yy05, statePenaltyMPC_bmy05);
statePenaltyMPC_LA_DENSE_BACKWARDSUB_8(statePenaltyMPC_Ld05, statePenaltyMPC_bmy05, statePenaltyMPC_dvcc05);
statePenaltyMPC_LA_DENSE_MTVMSUB_8_8(statePenaltyMPC_Lsd05, statePenaltyMPC_dvcc05, statePenaltyMPC_yy04, statePenaltyMPC_bmy04);
statePenaltyMPC_LA_DENSE_BACKWARDSUB_8(statePenaltyMPC_Ld04, statePenaltyMPC_bmy04, statePenaltyMPC_dvcc04);
statePenaltyMPC_LA_DENSE_MTVMSUB_8_8(statePenaltyMPC_Lsd04, statePenaltyMPC_dvcc04, statePenaltyMPC_yy03, statePenaltyMPC_bmy03);
statePenaltyMPC_LA_DENSE_BACKWARDSUB_8(statePenaltyMPC_Ld03, statePenaltyMPC_bmy03, statePenaltyMPC_dvcc03);
statePenaltyMPC_LA_DENSE_MTVMSUB_8_8(statePenaltyMPC_Lsd03, statePenaltyMPC_dvcc03, statePenaltyMPC_yy02, statePenaltyMPC_bmy02);
statePenaltyMPC_LA_DENSE_BACKWARDSUB_8(statePenaltyMPC_Ld02, statePenaltyMPC_bmy02, statePenaltyMPC_dvcc02);
statePenaltyMPC_LA_DENSE_MTVMSUB_8_8(statePenaltyMPC_Lsd02, statePenaltyMPC_dvcc02, statePenaltyMPC_yy01, statePenaltyMPC_bmy01);
statePenaltyMPC_LA_DENSE_BACKWARDSUB_8(statePenaltyMPC_Ld01, statePenaltyMPC_bmy01, statePenaltyMPC_dvcc01);
statePenaltyMPC_LA_DENSE_MTVMSUB_8_8(statePenaltyMPC_Lsd01, statePenaltyMPC_dvcc01, statePenaltyMPC_yy00, statePenaltyMPC_bmy00);
statePenaltyMPC_LA_DENSE_BACKWARDSUB_8(statePenaltyMPC_Ld00, statePenaltyMPC_bmy00, statePenaltyMPC_dvcc00);
statePenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_8_26_8(params->C1, statePenaltyMPC_dvcc01, statePenaltyMPC_D00, statePenaltyMPC_dvcc00, statePenaltyMPC_grad_eq00);
statePenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_8_26_8(params->C2, statePenaltyMPC_dvcc02, statePenaltyMPC_D01, statePenaltyMPC_dvcc01, statePenaltyMPC_grad_eq01);
statePenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_8_26_8(params->C3, statePenaltyMPC_dvcc03, statePenaltyMPC_D01, statePenaltyMPC_dvcc02, statePenaltyMPC_grad_eq02);
statePenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_8_26_8(params->C4, statePenaltyMPC_dvcc04, statePenaltyMPC_D01, statePenaltyMPC_dvcc03, statePenaltyMPC_grad_eq03);
statePenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_8_26_8(params->C5, statePenaltyMPC_dvcc05, statePenaltyMPC_D01, statePenaltyMPC_dvcc04, statePenaltyMPC_grad_eq04);
statePenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_8_26_8(params->C6, statePenaltyMPC_dvcc06, statePenaltyMPC_D01, statePenaltyMPC_dvcc05, statePenaltyMPC_grad_eq05);
statePenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_8_26_8(params->C7, statePenaltyMPC_dvcc07, statePenaltyMPC_D01, statePenaltyMPC_dvcc06, statePenaltyMPC_grad_eq06);
statePenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_8_26_8(params->C8, statePenaltyMPC_dvcc08, statePenaltyMPC_D01, statePenaltyMPC_dvcc07, statePenaltyMPC_grad_eq07);
statePenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_8_26_8(params->C9, statePenaltyMPC_dvcc09, statePenaltyMPC_D01, statePenaltyMPC_dvcc08, statePenaltyMPC_grad_eq08);
statePenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_8_26_8(params->C10, statePenaltyMPC_dvcc10, statePenaltyMPC_D01, statePenaltyMPC_dvcc09, statePenaltyMPC_grad_eq09);
statePenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_8_26_8(params->C11, statePenaltyMPC_dvcc11, statePenaltyMPC_D01, statePenaltyMPC_dvcc10, statePenaltyMPC_grad_eq10);
statePenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_8_26_8(params->C12, statePenaltyMPC_dvcc12, statePenaltyMPC_D01, statePenaltyMPC_dvcc11, statePenaltyMPC_grad_eq11);
statePenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_8_26_8(params->C13, statePenaltyMPC_dvcc13, statePenaltyMPC_D01, statePenaltyMPC_dvcc12, statePenaltyMPC_grad_eq12);
statePenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_8_26_8(params->C14, statePenaltyMPC_dvcc14, statePenaltyMPC_D01, statePenaltyMPC_dvcc13, statePenaltyMPC_grad_eq13);
statePenaltyMPC_LA_DIAGZERO_MTVM_8_8(statePenaltyMPC_D14, statePenaltyMPC_dvcc14, statePenaltyMPC_grad_eq14);
statePenaltyMPC_LA_VSUB_372(statePenaltyMPC_rd, statePenaltyMPC_grad_eq, statePenaltyMPC_rd);
statePenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_26(statePenaltyMPC_Phi00, statePenaltyMPC_rd00, statePenaltyMPC_dzcc00);
statePenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_26(statePenaltyMPC_Phi01, statePenaltyMPC_rd01, statePenaltyMPC_dzcc01);
statePenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_26(statePenaltyMPC_Phi02, statePenaltyMPC_rd02, statePenaltyMPC_dzcc02);
statePenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_26(statePenaltyMPC_Phi03, statePenaltyMPC_rd03, statePenaltyMPC_dzcc03);
statePenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_26(statePenaltyMPC_Phi04, statePenaltyMPC_rd04, statePenaltyMPC_dzcc04);
statePenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_26(statePenaltyMPC_Phi05, statePenaltyMPC_rd05, statePenaltyMPC_dzcc05);
statePenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_26(statePenaltyMPC_Phi06, statePenaltyMPC_rd06, statePenaltyMPC_dzcc06);
statePenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_26(statePenaltyMPC_Phi07, statePenaltyMPC_rd07, statePenaltyMPC_dzcc07);
statePenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_26(statePenaltyMPC_Phi08, statePenaltyMPC_rd08, statePenaltyMPC_dzcc08);
statePenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_26(statePenaltyMPC_Phi09, statePenaltyMPC_rd09, statePenaltyMPC_dzcc09);
statePenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_26(statePenaltyMPC_Phi10, statePenaltyMPC_rd10, statePenaltyMPC_dzcc10);
statePenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_26(statePenaltyMPC_Phi11, statePenaltyMPC_rd11, statePenaltyMPC_dzcc11);
statePenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_26(statePenaltyMPC_Phi12, statePenaltyMPC_rd12, statePenaltyMPC_dzcc12);
statePenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_26(statePenaltyMPC_Phi13, statePenaltyMPC_rd13, statePenaltyMPC_dzcc13);
statePenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_8(statePenaltyMPC_Phi14, statePenaltyMPC_rd14, statePenaltyMPC_dzcc14);
statePenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_26(statePenaltyMPC_ccrhsl00, statePenaltyMPC_slb00, statePenaltyMPC_llbbyslb00, statePenaltyMPC_dzcc00, statePenaltyMPC_lbIdx00, statePenaltyMPC_dllbcc00);
statePenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_10(statePenaltyMPC_ccrhsub00, statePenaltyMPC_sub00, statePenaltyMPC_lubbysub00, statePenaltyMPC_dzcc00, statePenaltyMPC_ubIdx00, statePenaltyMPC_dlubcc00);
statePenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_26(statePenaltyMPC_ccrhsl01, statePenaltyMPC_slb01, statePenaltyMPC_llbbyslb01, statePenaltyMPC_dzcc01, statePenaltyMPC_lbIdx01, statePenaltyMPC_dllbcc01);
statePenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_10(statePenaltyMPC_ccrhsub01, statePenaltyMPC_sub01, statePenaltyMPC_lubbysub01, statePenaltyMPC_dzcc01, statePenaltyMPC_ubIdx01, statePenaltyMPC_dlubcc01);
statePenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_26(statePenaltyMPC_ccrhsl02, statePenaltyMPC_slb02, statePenaltyMPC_llbbyslb02, statePenaltyMPC_dzcc02, statePenaltyMPC_lbIdx02, statePenaltyMPC_dllbcc02);
statePenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_10(statePenaltyMPC_ccrhsub02, statePenaltyMPC_sub02, statePenaltyMPC_lubbysub02, statePenaltyMPC_dzcc02, statePenaltyMPC_ubIdx02, statePenaltyMPC_dlubcc02);
statePenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_26(statePenaltyMPC_ccrhsl03, statePenaltyMPC_slb03, statePenaltyMPC_llbbyslb03, statePenaltyMPC_dzcc03, statePenaltyMPC_lbIdx03, statePenaltyMPC_dllbcc03);
statePenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_10(statePenaltyMPC_ccrhsub03, statePenaltyMPC_sub03, statePenaltyMPC_lubbysub03, statePenaltyMPC_dzcc03, statePenaltyMPC_ubIdx03, statePenaltyMPC_dlubcc03);
statePenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_26(statePenaltyMPC_ccrhsl04, statePenaltyMPC_slb04, statePenaltyMPC_llbbyslb04, statePenaltyMPC_dzcc04, statePenaltyMPC_lbIdx04, statePenaltyMPC_dllbcc04);
statePenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_10(statePenaltyMPC_ccrhsub04, statePenaltyMPC_sub04, statePenaltyMPC_lubbysub04, statePenaltyMPC_dzcc04, statePenaltyMPC_ubIdx04, statePenaltyMPC_dlubcc04);
statePenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_26(statePenaltyMPC_ccrhsl05, statePenaltyMPC_slb05, statePenaltyMPC_llbbyslb05, statePenaltyMPC_dzcc05, statePenaltyMPC_lbIdx05, statePenaltyMPC_dllbcc05);
statePenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_10(statePenaltyMPC_ccrhsub05, statePenaltyMPC_sub05, statePenaltyMPC_lubbysub05, statePenaltyMPC_dzcc05, statePenaltyMPC_ubIdx05, statePenaltyMPC_dlubcc05);
statePenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_26(statePenaltyMPC_ccrhsl06, statePenaltyMPC_slb06, statePenaltyMPC_llbbyslb06, statePenaltyMPC_dzcc06, statePenaltyMPC_lbIdx06, statePenaltyMPC_dllbcc06);
statePenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_10(statePenaltyMPC_ccrhsub06, statePenaltyMPC_sub06, statePenaltyMPC_lubbysub06, statePenaltyMPC_dzcc06, statePenaltyMPC_ubIdx06, statePenaltyMPC_dlubcc06);
statePenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_26(statePenaltyMPC_ccrhsl07, statePenaltyMPC_slb07, statePenaltyMPC_llbbyslb07, statePenaltyMPC_dzcc07, statePenaltyMPC_lbIdx07, statePenaltyMPC_dllbcc07);
statePenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_10(statePenaltyMPC_ccrhsub07, statePenaltyMPC_sub07, statePenaltyMPC_lubbysub07, statePenaltyMPC_dzcc07, statePenaltyMPC_ubIdx07, statePenaltyMPC_dlubcc07);
statePenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_26(statePenaltyMPC_ccrhsl08, statePenaltyMPC_slb08, statePenaltyMPC_llbbyslb08, statePenaltyMPC_dzcc08, statePenaltyMPC_lbIdx08, statePenaltyMPC_dllbcc08);
statePenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_10(statePenaltyMPC_ccrhsub08, statePenaltyMPC_sub08, statePenaltyMPC_lubbysub08, statePenaltyMPC_dzcc08, statePenaltyMPC_ubIdx08, statePenaltyMPC_dlubcc08);
statePenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_26(statePenaltyMPC_ccrhsl09, statePenaltyMPC_slb09, statePenaltyMPC_llbbyslb09, statePenaltyMPC_dzcc09, statePenaltyMPC_lbIdx09, statePenaltyMPC_dllbcc09);
statePenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_10(statePenaltyMPC_ccrhsub09, statePenaltyMPC_sub09, statePenaltyMPC_lubbysub09, statePenaltyMPC_dzcc09, statePenaltyMPC_ubIdx09, statePenaltyMPC_dlubcc09);
statePenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_26(statePenaltyMPC_ccrhsl10, statePenaltyMPC_slb10, statePenaltyMPC_llbbyslb10, statePenaltyMPC_dzcc10, statePenaltyMPC_lbIdx10, statePenaltyMPC_dllbcc10);
statePenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_10(statePenaltyMPC_ccrhsub10, statePenaltyMPC_sub10, statePenaltyMPC_lubbysub10, statePenaltyMPC_dzcc10, statePenaltyMPC_ubIdx10, statePenaltyMPC_dlubcc10);
statePenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_26(statePenaltyMPC_ccrhsl11, statePenaltyMPC_slb11, statePenaltyMPC_llbbyslb11, statePenaltyMPC_dzcc11, statePenaltyMPC_lbIdx11, statePenaltyMPC_dllbcc11);
statePenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_10(statePenaltyMPC_ccrhsub11, statePenaltyMPC_sub11, statePenaltyMPC_lubbysub11, statePenaltyMPC_dzcc11, statePenaltyMPC_ubIdx11, statePenaltyMPC_dlubcc11);
statePenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_26(statePenaltyMPC_ccrhsl12, statePenaltyMPC_slb12, statePenaltyMPC_llbbyslb12, statePenaltyMPC_dzcc12, statePenaltyMPC_lbIdx12, statePenaltyMPC_dllbcc12);
statePenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_10(statePenaltyMPC_ccrhsub12, statePenaltyMPC_sub12, statePenaltyMPC_lubbysub12, statePenaltyMPC_dzcc12, statePenaltyMPC_ubIdx12, statePenaltyMPC_dlubcc12);
statePenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_26(statePenaltyMPC_ccrhsl13, statePenaltyMPC_slb13, statePenaltyMPC_llbbyslb13, statePenaltyMPC_dzcc13, statePenaltyMPC_lbIdx13, statePenaltyMPC_dllbcc13);
statePenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_10(statePenaltyMPC_ccrhsub13, statePenaltyMPC_sub13, statePenaltyMPC_lubbysub13, statePenaltyMPC_dzcc13, statePenaltyMPC_ubIdx13, statePenaltyMPC_dlubcc13);
statePenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_8(statePenaltyMPC_ccrhsl14, statePenaltyMPC_slb14, statePenaltyMPC_llbbyslb14, statePenaltyMPC_dzcc14, statePenaltyMPC_lbIdx14, statePenaltyMPC_dllbcc14);
statePenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_8(statePenaltyMPC_ccrhsub14, statePenaltyMPC_sub14, statePenaltyMPC_lubbysub14, statePenaltyMPC_dzcc14, statePenaltyMPC_ubIdx14, statePenaltyMPC_dlubcc14);
statePenaltyMPC_LA_VSUB7_520(statePenaltyMPC_l, statePenaltyMPC_ccrhs, statePenaltyMPC_s, statePenaltyMPC_dl_cc, statePenaltyMPC_ds_cc);
statePenaltyMPC_LA_VADD_372(statePenaltyMPC_dz_cc, statePenaltyMPC_dz_aff);
statePenaltyMPC_LA_VADD_120(statePenaltyMPC_dv_cc, statePenaltyMPC_dv_aff);
statePenaltyMPC_LA_VADD_520(statePenaltyMPC_dl_cc, statePenaltyMPC_dl_aff);
statePenaltyMPC_LA_VADD_520(statePenaltyMPC_ds_cc, statePenaltyMPC_ds_aff);
info->lsit_cc = statePenaltyMPC_LINESEARCH_BACKTRACKING_COMBINED(statePenaltyMPC_z, statePenaltyMPC_v, statePenaltyMPC_l, statePenaltyMPC_s, statePenaltyMPC_dz_cc, statePenaltyMPC_dv_cc, statePenaltyMPC_dl_cc, statePenaltyMPC_ds_cc, &info->step_cc, &info->mu);
if( info->lsit_cc == statePenaltyMPC_NOPROGRESS ){
PRINTTEXT("Line search could not proceed at iteration %d, exiting.\n",info->it+1);
exitcode = statePenaltyMPC_NOPROGRESS; break;
}
info->it++;
}
output->z1[0] = statePenaltyMPC_z00[0];
output->z1[1] = statePenaltyMPC_z00[1];
output->z1[2] = statePenaltyMPC_z00[2];
output->z1[3] = statePenaltyMPC_z00[3];
output->z1[4] = statePenaltyMPC_z00[4];
output->z1[5] = statePenaltyMPC_z00[5];
output->z1[6] = statePenaltyMPC_z00[6];
output->z1[7] = statePenaltyMPC_z00[7];
output->z1[8] = statePenaltyMPC_z00[8];
output->z1[9] = statePenaltyMPC_z00[9];
output->z2[0] = statePenaltyMPC_z01[0];
output->z2[1] = statePenaltyMPC_z01[1];
output->z2[2] = statePenaltyMPC_z01[2];
output->z2[3] = statePenaltyMPC_z01[3];
output->z2[4] = statePenaltyMPC_z01[4];
output->z2[5] = statePenaltyMPC_z01[5];
output->z2[6] = statePenaltyMPC_z01[6];
output->z2[7] = statePenaltyMPC_z01[7];
output->z2[8] = statePenaltyMPC_z01[8];
output->z2[9] = statePenaltyMPC_z01[9];
output->z3[0] = statePenaltyMPC_z02[0];
output->z3[1] = statePenaltyMPC_z02[1];
output->z3[2] = statePenaltyMPC_z02[2];
output->z3[3] = statePenaltyMPC_z02[3];
output->z3[4] = statePenaltyMPC_z02[4];
output->z3[5] = statePenaltyMPC_z02[5];
output->z3[6] = statePenaltyMPC_z02[6];
output->z3[7] = statePenaltyMPC_z02[7];
output->z3[8] = statePenaltyMPC_z02[8];
output->z3[9] = statePenaltyMPC_z02[9];
output->z4[0] = statePenaltyMPC_z03[0];
output->z4[1] = statePenaltyMPC_z03[1];
output->z4[2] = statePenaltyMPC_z03[2];
output->z4[3] = statePenaltyMPC_z03[3];
output->z4[4] = statePenaltyMPC_z03[4];
output->z4[5] = statePenaltyMPC_z03[5];
output->z4[6] = statePenaltyMPC_z03[6];
output->z4[7] = statePenaltyMPC_z03[7];
output->z4[8] = statePenaltyMPC_z03[8];
output->z4[9] = statePenaltyMPC_z03[9];
output->z5[0] = statePenaltyMPC_z04[0];
output->z5[1] = statePenaltyMPC_z04[1];
output->z5[2] = statePenaltyMPC_z04[2];
output->z5[3] = statePenaltyMPC_z04[3];
output->z5[4] = statePenaltyMPC_z04[4];
output->z5[5] = statePenaltyMPC_z04[5];
output->z5[6] = statePenaltyMPC_z04[6];
output->z5[7] = statePenaltyMPC_z04[7];
output->z5[8] = statePenaltyMPC_z04[8];
output->z5[9] = statePenaltyMPC_z04[9];
output->z6[0] = statePenaltyMPC_z05[0];
output->z6[1] = statePenaltyMPC_z05[1];
output->z6[2] = statePenaltyMPC_z05[2];
output->z6[3] = statePenaltyMPC_z05[3];
output->z6[4] = statePenaltyMPC_z05[4];
output->z6[5] = statePenaltyMPC_z05[5];
output->z6[6] = statePenaltyMPC_z05[6];
output->z6[7] = statePenaltyMPC_z05[7];
output->z6[8] = statePenaltyMPC_z05[8];
output->z6[9] = statePenaltyMPC_z05[9];
output->z7[0] = statePenaltyMPC_z06[0];
output->z7[1] = statePenaltyMPC_z06[1];
output->z7[2] = statePenaltyMPC_z06[2];
output->z7[3] = statePenaltyMPC_z06[3];
output->z7[4] = statePenaltyMPC_z06[4];
output->z7[5] = statePenaltyMPC_z06[5];
output->z7[6] = statePenaltyMPC_z06[6];
output->z7[7] = statePenaltyMPC_z06[7];
output->z7[8] = statePenaltyMPC_z06[8];
output->z7[9] = statePenaltyMPC_z06[9];
output->z8[0] = statePenaltyMPC_z07[0];
output->z8[1] = statePenaltyMPC_z07[1];
output->z8[2] = statePenaltyMPC_z07[2];
output->z8[3] = statePenaltyMPC_z07[3];
output->z8[4] = statePenaltyMPC_z07[4];
output->z8[5] = statePenaltyMPC_z07[5];
output->z8[6] = statePenaltyMPC_z07[6];
output->z8[7] = statePenaltyMPC_z07[7];
output->z8[8] = statePenaltyMPC_z07[8];
output->z8[9] = statePenaltyMPC_z07[9];
output->z9[0] = statePenaltyMPC_z08[0];
output->z9[1] = statePenaltyMPC_z08[1];
output->z9[2] = statePenaltyMPC_z08[2];
output->z9[3] = statePenaltyMPC_z08[3];
output->z9[4] = statePenaltyMPC_z08[4];
output->z9[5] = statePenaltyMPC_z08[5];
output->z9[6] = statePenaltyMPC_z08[6];
output->z9[7] = statePenaltyMPC_z08[7];
output->z9[8] = statePenaltyMPC_z08[8];
output->z9[9] = statePenaltyMPC_z08[9];
output->z10[0] = statePenaltyMPC_z09[0];
output->z10[1] = statePenaltyMPC_z09[1];
output->z10[2] = statePenaltyMPC_z09[2];
output->z10[3] = statePenaltyMPC_z09[3];
output->z10[4] = statePenaltyMPC_z09[4];
output->z10[5] = statePenaltyMPC_z09[5];
output->z10[6] = statePenaltyMPC_z09[6];
output->z10[7] = statePenaltyMPC_z09[7];
output->z10[8] = statePenaltyMPC_z09[8];
output->z10[9] = statePenaltyMPC_z09[9];
output->z11[0] = statePenaltyMPC_z10[0];
output->z11[1] = statePenaltyMPC_z10[1];
output->z11[2] = statePenaltyMPC_z10[2];
output->z11[3] = statePenaltyMPC_z10[3];
output->z11[4] = statePenaltyMPC_z10[4];
output->z11[5] = statePenaltyMPC_z10[5];
output->z11[6] = statePenaltyMPC_z10[6];
output->z11[7] = statePenaltyMPC_z10[7];
output->z11[8] = statePenaltyMPC_z10[8];
output->z11[9] = statePenaltyMPC_z10[9];
output->z12[0] = statePenaltyMPC_z11[0];
output->z12[1] = statePenaltyMPC_z11[1];
output->z12[2] = statePenaltyMPC_z11[2];
output->z12[3] = statePenaltyMPC_z11[3];
output->z12[4] = statePenaltyMPC_z11[4];
output->z12[5] = statePenaltyMPC_z11[5];
output->z12[6] = statePenaltyMPC_z11[6];
output->z12[7] = statePenaltyMPC_z11[7];
output->z12[8] = statePenaltyMPC_z11[8];
output->z12[9] = statePenaltyMPC_z11[9];
output->z13[0] = statePenaltyMPC_z12[0];
output->z13[1] = statePenaltyMPC_z12[1];
output->z13[2] = statePenaltyMPC_z12[2];
output->z13[3] = statePenaltyMPC_z12[3];
output->z13[4] = statePenaltyMPC_z12[4];
output->z13[5] = statePenaltyMPC_z12[5];
output->z13[6] = statePenaltyMPC_z12[6];
output->z13[7] = statePenaltyMPC_z12[7];
output->z13[8] = statePenaltyMPC_z12[8];
output->z13[9] = statePenaltyMPC_z12[9];
output->z14[0] = statePenaltyMPC_z13[0];
output->z14[1] = statePenaltyMPC_z13[1];
output->z14[2] = statePenaltyMPC_z13[2];
output->z14[3] = statePenaltyMPC_z13[3];
output->z14[4] = statePenaltyMPC_z13[4];
output->z14[5] = statePenaltyMPC_z13[5];
output->z14[6] = statePenaltyMPC_z13[6];
output->z14[7] = statePenaltyMPC_z13[7];
output->z14[8] = statePenaltyMPC_z13[8];
output->z14[9] = statePenaltyMPC_z13[9];
output->z15[0] = statePenaltyMPC_z14[0];
output->z15[1] = statePenaltyMPC_z14[1];
output->z15[2] = statePenaltyMPC_z14[2];
output->z15[3] = statePenaltyMPC_z14[3];
output->z15[4] = statePenaltyMPC_z14[4];
output->z15[5] = statePenaltyMPC_z14[5];
output->z15[6] = statePenaltyMPC_z14[6];
output->z15[7] = statePenaltyMPC_z14[7];

#if statePenaltyMPC_SET_TIMING == 1
info->solvetime = statePenaltyMPC_toc(&solvertimer);
#if statePenaltyMPC_SET_PRINTLEVEL > 0 && statePenaltyMPC_SET_TIMING == 1
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
