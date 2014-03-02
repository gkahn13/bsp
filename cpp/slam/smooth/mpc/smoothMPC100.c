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

#include "smoothMPC.h"

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
 * Initializes a vector of length 1092 with a value.
 */
void smoothMPC_LA_INITIALIZEVECTOR_1092(smoothMPC_FLOAT* vec, smoothMPC_FLOAT value)
{
	int i;
	for( i=0; i<1092; i++ )
	{
		vec[i] = value;
	}
}


/*
 * Initializes a vector of length 300 with a value.
 */
void smoothMPC_LA_INITIALIZEVECTOR_300(smoothMPC_FLOAT* vec, smoothMPC_FLOAT value)
{
	int i;
	for( i=0; i<300; i++ )
	{
		vec[i] = value;
	}
}


/*
 * Initializes a vector of length 1590 with a value.
 */
void smoothMPC_LA_INITIALIZEVECTOR_1590(smoothMPC_FLOAT* vec, smoothMPC_FLOAT value)
{
	int i;
	for( i=0; i<1590; i++ )
	{
		vec[i] = value;
	}
}


/* 
 * Calculates a dot product and adds it to a variable: z += x'*y; 
 * This function is for vectors of length 1590.
 */
void smoothMPC_LA_DOTACC_1590(smoothMPC_FLOAT *x, smoothMPC_FLOAT *y, smoothMPC_FLOAT *z)
{
	int i;
	for( i=0; i<1590; i++ ){
		*z += x[i]*y[i];
	}
}


/*
 * Calculates the gradient and the value for a quadratic function 0.5*z'*H*z + f'*z
 *
 * INPUTS:     H  - Symmetric Hessian, diag matrix of size [11 x 11]
 *             f  - column vector of size 11
 *             z  - column vector of size 11
 *
 * OUTPUTS: grad  - gradient at z (= H*z + f), column vector of size 11
 *          value <-- value + 0.5*z'*H*z + f'*z (value will be modified)
 */
void smoothMPC_LA_DIAG_QUADFCN_11(smoothMPC_FLOAT* H, smoothMPC_FLOAT* f, smoothMPC_FLOAT* z, smoothMPC_FLOAT* grad, smoothMPC_FLOAT* value)
{
	int i;
	smoothMPC_FLOAT hz;	
	for( i=0; i<11; i++){
		hz = H[i]*z[i];
		grad[i] = hz + f[i];
		*value += 0.5*hz*z[i] + f[i]*z[i];
	}
}


/*
 * Calculates the gradient and the value for a quadratic function 0.5*z'*H*z + f'*z
 *
 * INPUTS:     H  - Symmetric Hessian, diag matrix of size [3 x 3]
 *             f  - column vector of size 3
 *             z  - column vector of size 3
 *
 * OUTPUTS: grad  - gradient at z (= H*z + f), column vector of size 3
 *          value <-- value + 0.5*z'*H*z + f'*z (value will be modified)
 */
void smoothMPC_LA_DIAG_QUADFCN_3(smoothMPC_FLOAT* H, smoothMPC_FLOAT* f, smoothMPC_FLOAT* z, smoothMPC_FLOAT* grad, smoothMPC_FLOAT* value)
{
	int i;
	smoothMPC_FLOAT hz;	
	for( i=0; i<3; i++){
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
void smoothMPC_LA_DIAGZERO_MVMSUB6_3(smoothMPC_FLOAT *B, smoothMPC_FLOAT *u, smoothMPC_FLOAT *b, smoothMPC_FLOAT *l, smoothMPC_FLOAT *r, smoothMPC_FLOAT *z, smoothMPC_FLOAT *y)
{
	int i;
	smoothMPC_FLOAT Bu[3];
	smoothMPC_FLOAT norm = *y;
	smoothMPC_FLOAT lr = 0;

	/* do A*x + B*u first */
	for( i=0; i<3; i++ ){
		Bu[i] = B[i]*u[i];
	}	

	for( i=0; i<3; i++ ){
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
void smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(smoothMPC_FLOAT *A, smoothMPC_FLOAT *x, smoothMPC_FLOAT *B, smoothMPC_FLOAT *u, smoothMPC_FLOAT *b, smoothMPC_FLOAT *l, smoothMPC_FLOAT *r, smoothMPC_FLOAT *z, smoothMPC_FLOAT *y)
{
	int i;
	int j;
	int k = 0;
	smoothMPC_FLOAT AxBu[3];
	smoothMPC_FLOAT norm = *y;
	smoothMPC_FLOAT lr = 0;

	/* do A*x + B*u first */
	for( i=0; i<3; i++ ){
		AxBu[i] = A[k++]*x[0] + B[i]*u[i];
	}	

	for( j=1; j<11; j++ ){		
		for( i=0; i<3; i++ ){
			AxBu[i] += A[k++]*x[j];
		}
	}

	for( i=0; i<3; i++ ){
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
void smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_3(smoothMPC_FLOAT *A, smoothMPC_FLOAT *x, smoothMPC_FLOAT *B, smoothMPC_FLOAT *u, smoothMPC_FLOAT *b, smoothMPC_FLOAT *l, smoothMPC_FLOAT *r, smoothMPC_FLOAT *z, smoothMPC_FLOAT *y)
{
	int i;
	int j;
	int k = 0;
	smoothMPC_FLOAT AxBu[3];
	smoothMPC_FLOAT norm = *y;
	smoothMPC_FLOAT lr = 0;

	/* do A*x + B*u first */
	for( i=0; i<3; i++ ){
		AxBu[i] = A[k++]*x[0] + B[i]*u[i];
	}	

	for( j=1; j<11; j++ ){		
		for( i=0; i<3; i++ ){
			AxBu[i] += A[k++]*x[j];
		}
	}

	for( i=0; i<3; i++ ){
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
 * where A is of size [3 x 11] and stored in column major format.
 * and B is of size [3 x 11] and stored in diagzero format
 * Note the transposes of A and B!
 */
void smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(smoothMPC_FLOAT *A, smoothMPC_FLOAT *x, smoothMPC_FLOAT *B, smoothMPC_FLOAT *y, smoothMPC_FLOAT *z)
{
	int i;
	int j;
	int k = 0;
	for( i=0; i<3; i++ ){
		z[i] = 0;
		for( j=0; j<3; j++ ){
			z[i] += A[k++]*x[j];
		}
		z[i] += B[i]*y[i];
	}
	for( i=3 ;i<11; i++ ){
		z[i] = 0;
		for( j=0; j<3; j++ ){
			z[i] += A[k++]*x[j];
		}
	}
}


/*
 * Matrix vector multiplication y = M'*x where M is of size [3 x 3]
 * and stored in diagzero format. Note the transpose of M!
 */
void smoothMPC_LA_DIAGZERO_MTVM_3_3(smoothMPC_FLOAT *M, smoothMPC_FLOAT *x, smoothMPC_FLOAT *y)
{
	int i;
	for( i=0; i<3; i++ ){
		y[i] = M[i]*x[i];
	}
}


/*
 * Vector subtraction and addition.
 *	 Input: five vectors t, tidx, u, v, w and two scalars z and r
 *	 Output: y = t(tidx) - u + w
 *           z = z - v'*x;
 *           r = max([norm(y,inf), z]);
 * for vectors of length 11. Output z is of course scalar.
 */
void smoothMPC_LA_VSUBADD3_11(smoothMPC_FLOAT* t, smoothMPC_FLOAT* u, int* uidx, smoothMPC_FLOAT* v, smoothMPC_FLOAT* w, smoothMPC_FLOAT* y, smoothMPC_FLOAT* z, smoothMPC_FLOAT* r)
{
	int i;
	smoothMPC_FLOAT norm = *r;
	smoothMPC_FLOAT vx = 0;
	smoothMPC_FLOAT x;
	for( i=0; i<11; i++){
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
 * for vectors of length 5. Output z is of course scalar.
 */
void smoothMPC_LA_VSUBADD2_5(smoothMPC_FLOAT* t, int* tidx, smoothMPC_FLOAT* u, smoothMPC_FLOAT* v, smoothMPC_FLOAT* w, smoothMPC_FLOAT* y, smoothMPC_FLOAT* z, smoothMPC_FLOAT* r)
{
	int i;
	smoothMPC_FLOAT norm = *r;
	smoothMPC_FLOAT vx = 0;
	smoothMPC_FLOAT x;
	for( i=0; i<5; i++){
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
 * for vectors of length 3. Output z is of course scalar.
 */
void smoothMPC_LA_VSUBADD3_3(smoothMPC_FLOAT* t, smoothMPC_FLOAT* u, int* uidx, smoothMPC_FLOAT* v, smoothMPC_FLOAT* w, smoothMPC_FLOAT* y, smoothMPC_FLOAT* z, smoothMPC_FLOAT* r)
{
	int i;
	smoothMPC_FLOAT norm = *r;
	smoothMPC_FLOAT vx = 0;
	smoothMPC_FLOAT x;
	for( i=0; i<3; i++){
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
 * for vectors of length 3. Output z is of course scalar.
 */
void smoothMPC_LA_VSUBADD2_3(smoothMPC_FLOAT* t, int* tidx, smoothMPC_FLOAT* u, smoothMPC_FLOAT* v, smoothMPC_FLOAT* w, smoothMPC_FLOAT* y, smoothMPC_FLOAT* z, smoothMPC_FLOAT* r)
{
	int i;
	smoothMPC_FLOAT norm = *r;
	smoothMPC_FLOAT vx = 0;
	smoothMPC_FLOAT x;
	for( i=0; i<3; i++){
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
 * Special function for box constraints of length 11
 * Returns also L/S, a value that is often used elsewhere.
 */
void smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_FLOAT *lu, smoothMPC_FLOAT *su, smoothMPC_FLOAT *ru, smoothMPC_FLOAT *ll, smoothMPC_FLOAT *sl, smoothMPC_FLOAT *rl, int* lbIdx, int* ubIdx, smoothMPC_FLOAT *grad, smoothMPC_FLOAT *lubysu, smoothMPC_FLOAT *llbysl)
{
	int i;
	for( i=0; i<11; i++ ){
		grad[i] = 0;
	}
	for( i=0; i<11; i++ ){		
		llbysl[i] = ll[i] / sl[i];
		grad[lbIdx[i]] -= llbysl[i]*rl[i];
	}
	for( i=0; i<5; i++ ){
		lubysu[i] = lu[i] / su[i];
		grad[ubIdx[i]] += lubysu[i]*ru[i];
	}
}


/*
 * Computes inequality constraints gradient-
 * Special function for box constraints of length 3
 * Returns also L/S, a value that is often used elsewhere.
 */
void smoothMPC_LA_INEQ_B_GRAD_3_3_3(smoothMPC_FLOAT *lu, smoothMPC_FLOAT *su, smoothMPC_FLOAT *ru, smoothMPC_FLOAT *ll, smoothMPC_FLOAT *sl, smoothMPC_FLOAT *rl, int* lbIdx, int* ubIdx, smoothMPC_FLOAT *grad, smoothMPC_FLOAT *lubysu, smoothMPC_FLOAT *llbysl)
{
	int i;
	for( i=0; i<3; i++ ){
		grad[i] = 0;
	}
	for( i=0; i<3; i++ ){		
		llbysl[i] = ll[i] / sl[i];
		grad[lbIdx[i]] -= llbysl[i]*rl[i];
	}
	for( i=0; i<3; i++ ){
		lubysu[i] = lu[i] / su[i];
		grad[ubIdx[i]] += lubysu[i]*ru[i];
	}
}


/*
 * Addition of three vectors  z = u + w + v
 * of length 1092.
 */
void smoothMPC_LA_VVADD3_1092(smoothMPC_FLOAT *u, smoothMPC_FLOAT *v, smoothMPC_FLOAT *w, smoothMPC_FLOAT *z)
{
	int i;
	for( i=0; i<1092; i++ ){
		z[i] = u[i] + v[i] + w[i];
	}
}


/*
 * Special function to compute the diagonal cholesky factorization of the 
 * positive definite augmented Hessian for block size 11.
 *
 * Inputs: - H = diagonal cost Hessian in diagonal storage format
 *         - llbysl = L / S of lower bounds
 *         - lubysu = L / S of upper bounds
 *
 * Output: Phi = sqrt(H + diag(llbysl) + diag(lubysu))
 * where Phi is stored in diagonal storage format
 */
void smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(smoothMPC_FLOAT *H, smoothMPC_FLOAT *llbysl, int* lbIdx, smoothMPC_FLOAT *lubysu, int* ubIdx, smoothMPC_FLOAT *Phi)


{
	int i;
	
	/* copy  H into PHI */
	for( i=0; i<11; i++ ){
		Phi[i] = H[i];
	}

	/* add llbysl onto Phi where necessary */
	for( i=0; i<11; i++ ){
		Phi[lbIdx[i]] += llbysl[i];
	}

	/* add lubysu onto Phi where necessary */
	for( i=0; i<5; i++){
		Phi[ubIdx[i]] +=  lubysu[i];
	}
	
	/* compute cholesky */
	for(i=0; i<11; i++)
	{
#if smoothMPC_SET_PRINTLEVEL > 0 && defined PRINTNUMERICALWARNINGS
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
 * where A is to be computed and is of size [3 x 11],
 * B is given and of size [3 x 11], L is a diagonal
 * matrix of size 3 stored in diagonal matrix 
 * storage format. Note the transpose of L has no impact!
 *
 * Result: A in column major storage format.
 *
 */
void smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_FLOAT *L, smoothMPC_FLOAT *B, smoothMPC_FLOAT *A)
{
    int i,j;
	 int k = 0;

	for( j=0; j<11; j++){
		for( i=0; i<3; i++){
			A[k] = B[k]/L[j];
			k++;
		}
	}

}


/**
 * Forward substitution for the matrix equation A*L' = B
 * where A is to be computed and is of size [3 x 11],
 * B is given and of size [3 x 11], L is a diagonal
 *  matrix of size 11 stored in diagonal 
 * storage format. Note the transpose of L!
 *
 * Result: A in diagonalzero storage format.
 *
 */
void smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_FLOAT *L, smoothMPC_FLOAT *B, smoothMPC_FLOAT *A)
{
	int j;
    for( j=0; j<11; j++ ){   
		A[j] = B[j]/L[j];
     }
}


/**
 * Compute C = A*B' where 
 *
 *	size(A) = [3 x 11]
 *  size(B) = [3 x 11] in diagzero format
 * 
 * A and C matrices are stored in column major format.
 * 
 * 
 */
void smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_FLOAT *A, smoothMPC_FLOAT *B, smoothMPC_FLOAT *C)
{
    int i, j;
	
	for( i=0; i<3; i++ ){
		for( j=0; j<3; j++){
			C[j*3+i] = B[i*3+j]*A[i];
		}
	}

}


/**
 * Forward substitution to solve L*y = b where L is a
 * diagonal matrix in vector storage format.
 * 
 * The dimensions involved are 11.
 */
void smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_FLOAT *L, smoothMPC_FLOAT *b, smoothMPC_FLOAT *y)
{
    int i;

    for( i=0; i<11; i++ ){
		y[i] = b[i]/L[i];
    }
}


/*
 * Special function to compute the diagonal cholesky factorization of the 
 * positive definite augmented Hessian for block size 3.
 *
 * Inputs: - H = diagonal cost Hessian in diagonal storage format
 *         - llbysl = L / S of lower bounds
 *         - lubysu = L / S of upper bounds
 *
 * Output: Phi = sqrt(H + diag(llbysl) + diag(lubysu))
 * where Phi is stored in diagonal storage format
 */
void smoothMPC_LA_DIAG_CHOL_ONELOOP_LBUB_3_3_3(smoothMPC_FLOAT *H, smoothMPC_FLOAT *llbysl, int* lbIdx, smoothMPC_FLOAT *lubysu, int* ubIdx, smoothMPC_FLOAT *Phi)


{
	int i;
	
	/* compute cholesky */
	for( i=0; i<3; i++ ){
		Phi[i] = H[i] + llbysl[i] + lubysu[i];

#if smoothMPC_SET_PRINTLEVEL > 0 && defined PRINTNUMERICALWARNINGS
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
 * where A is to be computed and is of size [3 x 3],
 * B is given and of size [3 x 3], L is a diagonal
 *  matrix of size 3 stored in diagonal 
 * storage format. Note the transpose of L!
 *
 * Result: A in diagonalzero storage format.
 *
 */
void smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_3(smoothMPC_FLOAT *L, smoothMPC_FLOAT *B, smoothMPC_FLOAT *A)
{
	int j;
    for( j=0; j<3; j++ ){   
		A[j] = B[j]/L[j];
     }
}


/**
 * Forward substitution to solve L*y = b where L is a
 * diagonal matrix in vector storage format.
 * 
 * The dimensions involved are 3.
 */
void smoothMPC_LA_DIAG_FORWARDSUB_3(smoothMPC_FLOAT *L, smoothMPC_FLOAT *b, smoothMPC_FLOAT *y)
{
    int i;

    for( i=0; i<3; i++ ){
		y[i] = b[i]/L[i];
    }
}


/**
 * Compute L = B*B', where L is lower triangular of size NXp1
 * and B is a diagzero matrix of size [3 x 11] in column
 * storage format.
 * 
 */
void smoothMPC_LA_DIAGZERO_MMT_3(smoothMPC_FLOAT *B, smoothMPC_FLOAT *L)
{
    int i, ii, di;
    
    ii = 0; di = 0;
    for( i=0; i<3; i++ ){        
		L[ii+i] = B[i]*B[i];
        ii += ++di;
    }
}


/* 
 * Computes r = b - B*u
 * B is stored in diagzero format
 */
void smoothMPC_LA_DIAGZERO_MVMSUB7_3(smoothMPC_FLOAT *B, smoothMPC_FLOAT *u, smoothMPC_FLOAT *b, smoothMPC_FLOAT *r)
{
	int i;

	for( i=0; i<3; i++ ){
		r[i] = b[i] - B[i]*u[i];
	}	
	
}


/**
 * Compute L = A*A' + B*B', where L is lower triangular of size NXp1
 * and A is a dense matrix of size [3 x 11] in column
 * storage format, and B is of size [3 x 11] diagonalzero
 * storage format.
 * 
 * THIS ONE HAS THE WORST ACCES PATTERN POSSIBLE. 
 * POSSIBKE FIX: PUT A AND B INTO ROW MAJOR FORMAT FIRST.
 * 
 */
void smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_FLOAT *A, smoothMPC_FLOAT *B, smoothMPC_FLOAT *L)
{
    int i, j, k, ii, di;
    smoothMPC_FLOAT ltemp;
    
    ii = 0; di = 0;
    for( i=0; i<3; i++ ){        
        for( j=0; j<=i; j++ ){
            ltemp = 0; 
            for( k=0; k<11; k++ ){
                ltemp += A[k*3+i]*A[k*3+j];
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
void smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_FLOAT *A, smoothMPC_FLOAT *x, smoothMPC_FLOAT *B, smoothMPC_FLOAT *u, smoothMPC_FLOAT *b, smoothMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;

	for( i=0; i<3; i++ ){
		r[i] = b[i] - A[k++]*x[0] - B[i]*u[i];
	}	

	for( j=1; j<11; j++ ){		
		for( i=0; i<3; i++ ){
			r[i] -= A[k++]*x[j];
		}
	}
	
}


/**
 * Compute L = A*A' + B*B', where L is lower triangular of size NXp1
 * and A is a dense matrix of size [3 x 11] in column
 * storage format, and B is of size [3 x 3] diagonalzero
 * storage format.
 * 
 * THIS ONE HAS THE WORST ACCES PATTERN POSSIBLE. 
 * POSSIBKE FIX: PUT A AND B INTO ROW MAJOR FORMAT FIRST.
 * 
 */
void smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_3(smoothMPC_FLOAT *A, smoothMPC_FLOAT *B, smoothMPC_FLOAT *L)
{
    int i, j, k, ii, di;
    smoothMPC_FLOAT ltemp;
    
    ii = 0; di = 0;
    for( i=0; i<3; i++ ){        
        for( j=0; j<=i; j++ ){
            ltemp = 0; 
            for( k=0; k<11; k++ ){
                ltemp += A[k*3+i]*A[k*3+j];
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
void smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_3(smoothMPC_FLOAT *A, smoothMPC_FLOAT *x, smoothMPC_FLOAT *B, smoothMPC_FLOAT *u, smoothMPC_FLOAT *b, smoothMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;

	for( i=0; i<3; i++ ){
		r[i] = b[i] - A[k++]*x[0] - B[i]*u[i];
	}	

	for( j=1; j<11; j++ ){		
		for( i=0; i<3; i++ ){
			r[i] -= A[k++]*x[j];
		}
	}
	
}


/**
 * Cholesky factorization as above, but working on a matrix in 
 * lower triangular storage format of size 3 and outputting
 * the Cholesky factor to matrix L in lower triangular format.
 */
void smoothMPC_LA_DENSE_CHOL_3(smoothMPC_FLOAT *A, smoothMPC_FLOAT *L)
{
    int i, j, k, di, dj;
	 int ii, jj;

    smoothMPC_FLOAT l;
    smoothMPC_FLOAT Mii;

	/* copy A to L first and then operate on L */
	/* COULD BE OPTIMIZED */
	ii=0; di=0;
	for( i=0; i<3; i++ ){
		for( j=0; j<=i; j++ ){
			L[ii+j] = A[ii+j];
		}
		ii += ++di;
	}    
	
	/* factor L */
	ii=0; di=0;
    for( i=0; i<3; i++ ){
        l = 0;
        for( k=0; k<i; k++ ){
            l += L[ii+k]*L[ii+k];
        }        
        
        Mii = L[ii+i] - l;
        
#if smoothMPC_SET_PRINTLEVEL > 0 && defined PRINTNUMERICALWARNINGS
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
        for( j=i+1; j<3; j++ ){
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
 * The dimensions involved are 3.
 */
void smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_FLOAT *L, smoothMPC_FLOAT *b, smoothMPC_FLOAT *y)
{
    int i,j,ii,di;
    smoothMPC_FLOAT yel;
            
    ii = 0; di = 0;
    for( i=0; i<3; i++ ){
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
 * where A is to be computed and is of size [3 x 3],
 * B is given and of size [3 x 3], L is a lower tri-
 * angular matrix of size 3 stored in lower triangular 
 * storage format. Note the transpose of L AND B!
 *
 * Result: A in column major storage format.
 *
 */
void smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_FLOAT *L, smoothMPC_FLOAT *B, smoothMPC_FLOAT *A)
{
    int i,j,k,ii,di;
    smoothMPC_FLOAT a;
    
    ii=0; di=0;
    for( j=0; j<3; j++ ){        
        for( i=0; i<3; i++ ){
            a = B[i*3+j];
            for( k=0; k<j; k++ ){
                a -= A[k*3+i]*L[ii+k];
            }    

			/* saturate for numerical stability */
			a = MIN(a, BIGM);
			a = MAX(a, -BIGM); 

			A[j*3+i] = a/L[ii+j];			
        }
        ii += ++di;
    }
}


/**
 * Compute L = L - A*A', where L is lower triangular of size 3
 * and A is a dense matrix of size [3 x 3] in column
 * storage format.
 * 
 * THIS ONE HAS THE WORST ACCES PATTERN POSSIBLE. 
 * POSSIBKE FIX: PUT A INTO ROW MAJOR FORMAT FIRST.
 * 
 */
void smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_FLOAT *A, smoothMPC_FLOAT *L)
{
    int i, j, k, ii, di;
    smoothMPC_FLOAT ltemp;
    
    ii = 0; di = 0;
    for( i=0; i<3; i++ ){        
        for( j=0; j<=i; j++ ){
            ltemp = 0; 
            for( k=0; k<3; k++ ){
                ltemp += A[k*3+i]*A[k*3+j];
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
void smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_FLOAT *A, smoothMPC_FLOAT *x, smoothMPC_FLOAT *b, smoothMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;

	for( i=0; i<3; i++ ){
		r[i] = b[i] - A[k++]*x[0];
	}	
	for( j=1; j<3; j++ ){		
		for( i=0; i<3; i++ ){
			r[i] -= A[k++]*x[j];
		}
	}
}


/**
 * Backward Substitution to solve L^T*x = y where L is a
 * lower triangular matrix in triangular storage format.
 * 
 * All involved dimensions are 3.
 */
void smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_FLOAT *L, smoothMPC_FLOAT *y, smoothMPC_FLOAT *x)
{
    int i, ii, di, j, jj, dj;
    smoothMPC_FLOAT xel;    
	int start = 3;
    
    /* now solve L^T*x = y by backward substitution */
    ii = start; di = 2;
    for( i=2; i>=0; i-- ){        
        xel = y[i];        
        jj = start; dj = 2;
        for( j=2; j>i; j-- ){
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
 * Matrix vector multiplication y = b - M'*x where M is of size [3 x 3]
 * and stored in column major format. Note the transpose of M!
 */
void smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_FLOAT *A, smoothMPC_FLOAT *x, smoothMPC_FLOAT *b, smoothMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0; 
	for( i=0; i<3; i++ ){
		r[i] = b[i];
		for( j=0; j<3; j++ ){
			r[i] -= A[k++]*x[j];
		}
	}
}


/*
 * Vector subtraction z = -x - y for vectors of length 1092.
 */
void smoothMPC_LA_VSUB2_1092(smoothMPC_FLOAT *x, smoothMPC_FLOAT *y, smoothMPC_FLOAT *z)
{
	int i;
	for( i=0; i<1092; i++){
		z[i] = -x[i] - y[i];
	}
}


/**
 * Forward-Backward-Substitution to solve L*L^T*x = b where L is a
 * diagonal matrix of size 11 in vector
 * storage format.
 */
void smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_FLOAT *L, smoothMPC_FLOAT *b, smoothMPC_FLOAT *x)
{
    int i;
            
    /* solve Ly = b by forward and backward substitution */
    for( i=0; i<11; i++ ){
		x[i] = b[i]/(L[i]*L[i]);
    }
    
}


/**
 * Forward-Backward-Substitution to solve L*L^T*x = b where L is a
 * diagonal matrix of size 3 in vector
 * storage format.
 */
void smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_3(smoothMPC_FLOAT *L, smoothMPC_FLOAT *b, smoothMPC_FLOAT *x)
{
    int i;
            
    /* solve Ly = b by forward and backward substitution */
    for( i=0; i<3; i++ ){
		x[i] = b[i]/(L[i]*L[i]);
    }
    
}


/*
 * Vector subtraction z = x(xidx) - y where y, z and xidx are of length 11,
 * and x has length 11 and is indexed through yidx.
 */
void smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_FLOAT *x, int* xidx, smoothMPC_FLOAT *y, smoothMPC_FLOAT *z)
{
	int i;
	for( i=0; i<11; i++){
		z[i] = x[xidx[i]] - y[i];
	}
}


/*
 * Vector subtraction x = -u.*v - w for vectors of length 11.
 */
void smoothMPC_LA_VSUB3_11(smoothMPC_FLOAT *u, smoothMPC_FLOAT *v, smoothMPC_FLOAT *w, smoothMPC_FLOAT *x)
{
	int i;
	for( i=0; i<11; i++){
		x[i] = -u[i]*v[i] - w[i];
	}
}


/*
 * Vector subtraction z = -x - y(yidx) where y is of length 11
 * and z, x and yidx are of length 5.
 */
void smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_FLOAT *x, smoothMPC_FLOAT *y, int* yidx, smoothMPC_FLOAT *z)
{
	int i;
	for( i=0; i<5; i++){
		z[i] = -x[i] - y[yidx[i]];
	}
}


/*
 * Vector subtraction x = -u.*v - w for vectors of length 5.
 */
void smoothMPC_LA_VSUB3_5(smoothMPC_FLOAT *u, smoothMPC_FLOAT *v, smoothMPC_FLOAT *w, smoothMPC_FLOAT *x)
{
	int i;
	for( i=0; i<5; i++){
		x[i] = -u[i]*v[i] - w[i];
	}
}


/*
 * Vector subtraction z = x(xidx) - y where y, z and xidx are of length 3,
 * and x has length 3 and is indexed through yidx.
 */
void smoothMPC_LA_VSUB_INDEXED_3(smoothMPC_FLOAT *x, int* xidx, smoothMPC_FLOAT *y, smoothMPC_FLOAT *z)
{
	int i;
	for( i=0; i<3; i++){
		z[i] = x[xidx[i]] - y[i];
	}
}


/*
 * Vector subtraction x = -u.*v - w for vectors of length 3.
 */
void smoothMPC_LA_VSUB3_3(smoothMPC_FLOAT *u, smoothMPC_FLOAT *v, smoothMPC_FLOAT *w, smoothMPC_FLOAT *x)
{
	int i;
	for( i=0; i<3; i++){
		x[i] = -u[i]*v[i] - w[i];
	}
}


/*
 * Vector subtraction z = -x - y(yidx) where y is of length 3
 * and z, x and yidx are of length 3.
 */
void smoothMPC_LA_VSUB2_INDEXED_3(smoothMPC_FLOAT *x, smoothMPC_FLOAT *y, int* yidx, smoothMPC_FLOAT *z)
{
	int i;
	for( i=0; i<3; i++){
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
 * smoothMPC_NOPROGRESS (should be negative).
 */
int smoothMPC_LINESEARCH_BACKTRACKING_AFFINE(smoothMPC_FLOAT *l, smoothMPC_FLOAT *s, smoothMPC_FLOAT *dl, smoothMPC_FLOAT *ds, smoothMPC_FLOAT *a, smoothMPC_FLOAT *mu_aff)
{
    int i;
	int lsIt=1;    
    smoothMPC_FLOAT dltemp;
    smoothMPC_FLOAT dstemp;
    smoothMPC_FLOAT mya = 1.0;
    smoothMPC_FLOAT mymu;
        
    while( 1 ){                        

        /* 
         * Compute both snew and wnew together.
         * We compute also mu_affine along the way here, as the
         * values might be in registers, so it should be cheaper.
         */
        mymu = 0;
        for( i=0; i<1590; i++ ){
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
        if( i == 1590 ){
            break;
        } else {
            mya *= smoothMPC_SET_LS_SCALE_AFF;
            if( mya < smoothMPC_SET_LS_MINSTEP ){
                return smoothMPC_NOPROGRESS;
            }
        }
    }
    
    /* return new values and iteration counter */
    *a = mya;
    *mu_aff = mymu / (smoothMPC_FLOAT)1590;
    return lsIt;
}


/*
 * Vector subtraction x = u.*v - a where a is a scalar
*  and x,u,v are vectors of length 1590.
 */
void smoothMPC_LA_VSUB5_1590(smoothMPC_FLOAT *u, smoothMPC_FLOAT *v, smoothMPC_FLOAT a, smoothMPC_FLOAT *x)
{
	int i;
	for( i=0; i<1590; i++){
		x[i] = u[i]*v[i] - a;
	}
}


/*
 * Computes x=0; x(uidx) += u/su; x(vidx) -= v/sv where x is of length 11,
 * u, su, uidx are of length 5 and v, sv, vidx are of length 11.
 */
void smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_FLOAT *u, smoothMPC_FLOAT *su, int* uidx, smoothMPC_FLOAT *v, smoothMPC_FLOAT *sv, int* vidx, smoothMPC_FLOAT *x)
{
	int i;
	for( i=0; i<11; i++ ){
		x[i] = 0;
	}
	for( i=0; i<5; i++){
		x[uidx[i]] += u[i]/su[i];
	}
	for( i=0; i<11; i++){
		x[vidx[i]] -= v[i]/sv[i];
	}
}


/* 
 * Computes r =  B*u
 * where B is stored in diagzero format
 */
void smoothMPC_LA_DIAGZERO_MVM_3(smoothMPC_FLOAT *B, smoothMPC_FLOAT *u, smoothMPC_FLOAT *r)
{
	int i;

	for( i=0; i<3; i++ ){
		r[i] = B[i]*u[i];
	}	
	
}


/* 
 * Computes r = A*x + B*u
 * where A is stored in column major format
 * and B is stored in diagzero format
 */
void smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_FLOAT *A, smoothMPC_FLOAT *x, smoothMPC_FLOAT *B, smoothMPC_FLOAT *u, smoothMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;

	for( i=0; i<3; i++ ){
		r[i] = A[k++]*x[0] + B[i]*u[i];
	}	

	for( j=1; j<11; j++ ){		
		for( i=0; i<3; i++ ){
			r[i] += A[k++]*x[j];
		}
	}
	
}


/*
 * Computes x=0; x(uidx) += u/su; x(vidx) -= v/sv where x is of length 3,
 * u, su, uidx are of length 3 and v, sv, vidx are of length 3.
 */
void smoothMPC_LA_VSUB6_INDEXED_3_3_3(smoothMPC_FLOAT *u, smoothMPC_FLOAT *su, int* uidx, smoothMPC_FLOAT *v, smoothMPC_FLOAT *sv, int* vidx, smoothMPC_FLOAT *x)
{
	int i;
	for( i=0; i<3; i++ ){
		x[i] = 0;
	}
	for( i=0; i<3; i++){
		x[uidx[i]] += u[i]/su[i];
	}
	for( i=0; i<3; i++){
		x[vidx[i]] -= v[i]/sv[i];
	}
}


/* 
 * Computes r = A*x + B*u
 * where A is stored in column major format
 * and B is stored in diagzero format
 */
void smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_3(smoothMPC_FLOAT *A, smoothMPC_FLOAT *x, smoothMPC_FLOAT *B, smoothMPC_FLOAT *u, smoothMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;

	for( i=0; i<3; i++ ){
		r[i] = A[k++]*x[0] + B[i]*u[i];
	}	

	for( j=1; j<11; j++ ){		
		for( i=0; i<3; i++ ){
			r[i] += A[k++]*x[j];
		}
	}
	
}


/*
 * Vector subtraction z = x - y for vectors of length 1092.
 */
void smoothMPC_LA_VSUB_1092(smoothMPC_FLOAT *x, smoothMPC_FLOAT *y, smoothMPC_FLOAT *z)
{
	int i;
	for( i=0; i<1092; i++){
		z[i] = x[i] - y[i];
	}
}


/** 
 * Computes z = -r./s - u.*y(y)
 * where all vectors except of y are of length 11 (length of y >= 11).
 */
void smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_FLOAT *r, smoothMPC_FLOAT *s, smoothMPC_FLOAT *u, smoothMPC_FLOAT *y, int* yidx, smoothMPC_FLOAT *z)
{
	int i;
	for( i=0; i<11; i++ ){
		z[i] = -r[i]/s[i] - u[i]*y[yidx[i]];
	}
}


/** 
 * Computes z = -r./s + u.*y(y)
 * where all vectors except of y are of length 5 (length of y >= 5).
 */
void smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_FLOAT *r, smoothMPC_FLOAT *s, smoothMPC_FLOAT *u, smoothMPC_FLOAT *y, int* yidx, smoothMPC_FLOAT *z)
{
	int i;
	for( i=0; i<5; i++ ){
		z[i] = -r[i]/s[i] + u[i]*y[yidx[i]];
	}
}


/** 
 * Computes z = -r./s - u.*y(y)
 * where all vectors except of y are of length 3 (length of y >= 3).
 */
void smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_3(smoothMPC_FLOAT *r, smoothMPC_FLOAT *s, smoothMPC_FLOAT *u, smoothMPC_FLOAT *y, int* yidx, smoothMPC_FLOAT *z)
{
	int i;
	for( i=0; i<3; i++ ){
		z[i] = -r[i]/s[i] - u[i]*y[yidx[i]];
	}
}


/** 
 * Computes z = -r./s + u.*y(y)
 * where all vectors except of y are of length 3 (length of y >= 3).
 */
void smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_3(smoothMPC_FLOAT *r, smoothMPC_FLOAT *s, smoothMPC_FLOAT *u, smoothMPC_FLOAT *y, int* yidx, smoothMPC_FLOAT *z)
{
	int i;
	for( i=0; i<3; i++ ){
		z[i] = -r[i]/s[i] + u[i]*y[yidx[i]];
	}
}


/*
 * Computes ds = -l.\(r + s.*dl) for vectors of length 1590.
 */
void smoothMPC_LA_VSUB7_1590(smoothMPC_FLOAT *l, smoothMPC_FLOAT *r, smoothMPC_FLOAT *s, smoothMPC_FLOAT *dl, smoothMPC_FLOAT *ds)
{
	int i;
	for( i=0; i<1590; i++){
		ds[i] = -(r[i] + s[i]*dl[i])/l[i];
	}
}


/*
 * Vector addition x = x + y for vectors of length 1092.
 */
void smoothMPC_LA_VADD_1092(smoothMPC_FLOAT *x, smoothMPC_FLOAT *y)
{
	int i;
	for( i=0; i<1092; i++){
		x[i] += y[i];
	}
}


/*
 * Vector addition x = x + y for vectors of length 300.
 */
void smoothMPC_LA_VADD_300(smoothMPC_FLOAT *x, smoothMPC_FLOAT *y)
{
	int i;
	for( i=0; i<300; i++){
		x[i] += y[i];
	}
}


/*
 * Vector addition x = x + y for vectors of length 1590.
 */
void smoothMPC_LA_VADD_1590(smoothMPC_FLOAT *x, smoothMPC_FLOAT *y)
{
	int i;
	for( i=0; i<1590; i++){
		x[i] += y[i];
	}
}


/**
 * Backtracking line search for combined predictor/corrector step.
 * Update on variables with safety factor gamma (to keep us away from
 * boundary).
 */
int smoothMPC_LINESEARCH_BACKTRACKING_COMBINED(smoothMPC_FLOAT *z, smoothMPC_FLOAT *v, smoothMPC_FLOAT *l, smoothMPC_FLOAT *s, smoothMPC_FLOAT *dz, smoothMPC_FLOAT *dv, smoothMPC_FLOAT *dl, smoothMPC_FLOAT *ds, smoothMPC_FLOAT *a, smoothMPC_FLOAT *mu)
{
    int i, lsIt=1;       
    smoothMPC_FLOAT dltemp;
    smoothMPC_FLOAT dstemp;    
    smoothMPC_FLOAT a_gamma;
            
    *a = 1.0;
    while( 1 ){                        

        /* check whether search criterion is fulfilled */
        for( i=0; i<1590; i++ ){
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
        if( i == 1590 ){
            break;
        } else {
            *a *= smoothMPC_SET_LS_SCALE;
            if( *a < smoothMPC_SET_LS_MINSTEP ){
                return smoothMPC_NOPROGRESS;
            }
        }
    }
    
    /* update variables with safety margin */
    a_gamma = (*a)*smoothMPC_SET_LS_MAXSTEP;
    
    /* primal variables */
    for( i=0; i<1092; i++ ){
        z[i] += a_gamma*dz[i];
    }
    
    /* equality constraint multipliers */
    for( i=0; i<300; i++ ){
        v[i] += a_gamma*dv[i];
    }
    
    /* inequality constraint multipliers & slacks, also update mu */
    *mu = 0;
    for( i=0; i<1590; i++ ){
        dltemp = l[i] + a_gamma*dl[i]; l[i] = dltemp;
        dstemp = s[i] + a_gamma*ds[i]; s[i] = dstemp;
        *mu += dltemp*dstemp;
    }
    
    *a = a_gamma;
    *mu /= (smoothMPC_FLOAT)1590;
    return lsIt;
}




/* VARIABLE DEFINITIONS ------------------------------------------------ */
smoothMPC_FLOAT smoothMPC_z[1092];
smoothMPC_FLOAT smoothMPC_v[300];
smoothMPC_FLOAT smoothMPC_dz_aff[1092];
smoothMPC_FLOAT smoothMPC_dv_aff[300];
smoothMPC_FLOAT smoothMPC_grad_cost[1092];
smoothMPC_FLOAT smoothMPC_grad_eq[1092];
smoothMPC_FLOAT smoothMPC_rd[1092];
smoothMPC_FLOAT smoothMPC_l[1590];
smoothMPC_FLOAT smoothMPC_s[1590];
smoothMPC_FLOAT smoothMPC_lbys[1590];
smoothMPC_FLOAT smoothMPC_dl_aff[1590];
smoothMPC_FLOAT smoothMPC_ds_aff[1590];
smoothMPC_FLOAT smoothMPC_dz_cc[1092];
smoothMPC_FLOAT smoothMPC_dv_cc[300];
smoothMPC_FLOAT smoothMPC_dl_cc[1590];
smoothMPC_FLOAT smoothMPC_ds_cc[1590];
smoothMPC_FLOAT smoothMPC_ccrhs[1590];
smoothMPC_FLOAT smoothMPC_grad_ineq[1092];
smoothMPC_FLOAT* smoothMPC_z00 = smoothMPC_z + 0;
smoothMPC_FLOAT* smoothMPC_dzaff00 = smoothMPC_dz_aff + 0;
smoothMPC_FLOAT* smoothMPC_dzcc00 = smoothMPC_dz_cc + 0;
smoothMPC_FLOAT* smoothMPC_rd00 = smoothMPC_rd + 0;
smoothMPC_FLOAT smoothMPC_Lbyrd00[11];
smoothMPC_FLOAT* smoothMPC_grad_cost00 = smoothMPC_grad_cost + 0;
smoothMPC_FLOAT* smoothMPC_grad_eq00 = smoothMPC_grad_eq + 0;
smoothMPC_FLOAT* smoothMPC_grad_ineq00 = smoothMPC_grad_ineq + 0;
smoothMPC_FLOAT smoothMPC_ctv00[11];
smoothMPC_FLOAT* smoothMPC_v00 = smoothMPC_v + 0;
smoothMPC_FLOAT smoothMPC_re00[3];
smoothMPC_FLOAT smoothMPC_beta00[3];
smoothMPC_FLOAT smoothMPC_betacc00[3];
smoothMPC_FLOAT* smoothMPC_dvaff00 = smoothMPC_dv_aff + 0;
smoothMPC_FLOAT* smoothMPC_dvcc00 = smoothMPC_dv_cc + 0;
smoothMPC_FLOAT smoothMPC_V00[33];
smoothMPC_FLOAT smoothMPC_Yd00[6];
smoothMPC_FLOAT smoothMPC_Ld00[6];
smoothMPC_FLOAT smoothMPC_yy00[3];
smoothMPC_FLOAT smoothMPC_bmy00[3];
int smoothMPC_lbIdx00[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb00 = smoothMPC_l + 0;
smoothMPC_FLOAT* smoothMPC_slb00 = smoothMPC_s + 0;
smoothMPC_FLOAT* smoothMPC_llbbyslb00 = smoothMPC_lbys + 0;
smoothMPC_FLOAT smoothMPC_rilb00[11];
smoothMPC_FLOAT* smoothMPC_dllbaff00 = smoothMPC_dl_aff + 0;
smoothMPC_FLOAT* smoothMPC_dslbaff00 = smoothMPC_ds_aff + 0;
smoothMPC_FLOAT* smoothMPC_dllbcc00 = smoothMPC_dl_cc + 0;
smoothMPC_FLOAT* smoothMPC_dslbcc00 = smoothMPC_ds_cc + 0;
smoothMPC_FLOAT* smoothMPC_ccrhsl00 = smoothMPC_ccrhs + 0;
int smoothMPC_ubIdx00[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub00 = smoothMPC_l + 11;
smoothMPC_FLOAT* smoothMPC_sub00 = smoothMPC_s + 11;
smoothMPC_FLOAT* smoothMPC_lubbysub00 = smoothMPC_lbys + 11;
smoothMPC_FLOAT smoothMPC_riub00[5];
smoothMPC_FLOAT* smoothMPC_dlubaff00 = smoothMPC_dl_aff + 11;
smoothMPC_FLOAT* smoothMPC_dsubaff00 = smoothMPC_ds_aff + 11;
smoothMPC_FLOAT* smoothMPC_dlubcc00 = smoothMPC_dl_cc + 11;
smoothMPC_FLOAT* smoothMPC_dsubcc00 = smoothMPC_ds_cc + 11;
smoothMPC_FLOAT* smoothMPC_ccrhsub00 = smoothMPC_ccrhs + 11;
smoothMPC_FLOAT smoothMPC_Phi00[11];
smoothMPC_FLOAT smoothMPC_D00[11] = {1.0000000000000000E+000, 
1.0000000000000000E+000, 
1.0000000000000000E+000};
smoothMPC_FLOAT smoothMPC_W00[11];
smoothMPC_FLOAT* smoothMPC_z01 = smoothMPC_z + 11;
smoothMPC_FLOAT* smoothMPC_dzaff01 = smoothMPC_dz_aff + 11;
smoothMPC_FLOAT* smoothMPC_dzcc01 = smoothMPC_dz_cc + 11;
smoothMPC_FLOAT* smoothMPC_rd01 = smoothMPC_rd + 11;
smoothMPC_FLOAT smoothMPC_Lbyrd01[11];
smoothMPC_FLOAT* smoothMPC_grad_cost01 = smoothMPC_grad_cost + 11;
smoothMPC_FLOAT* smoothMPC_grad_eq01 = smoothMPC_grad_eq + 11;
smoothMPC_FLOAT* smoothMPC_grad_ineq01 = smoothMPC_grad_ineq + 11;
smoothMPC_FLOAT smoothMPC_ctv01[11];
smoothMPC_FLOAT* smoothMPC_v01 = smoothMPC_v + 3;
smoothMPC_FLOAT smoothMPC_re01[3];
smoothMPC_FLOAT smoothMPC_beta01[3];
smoothMPC_FLOAT smoothMPC_betacc01[3];
smoothMPC_FLOAT* smoothMPC_dvaff01 = smoothMPC_dv_aff + 3;
smoothMPC_FLOAT* smoothMPC_dvcc01 = smoothMPC_dv_cc + 3;
smoothMPC_FLOAT smoothMPC_V01[33];
smoothMPC_FLOAT smoothMPC_Yd01[6];
smoothMPC_FLOAT smoothMPC_Ld01[6];
smoothMPC_FLOAT smoothMPC_yy01[3];
smoothMPC_FLOAT smoothMPC_bmy01[3];
int smoothMPC_lbIdx01[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb01 = smoothMPC_l + 16;
smoothMPC_FLOAT* smoothMPC_slb01 = smoothMPC_s + 16;
smoothMPC_FLOAT* smoothMPC_llbbyslb01 = smoothMPC_lbys + 16;
smoothMPC_FLOAT smoothMPC_rilb01[11];
smoothMPC_FLOAT* smoothMPC_dllbaff01 = smoothMPC_dl_aff + 16;
smoothMPC_FLOAT* smoothMPC_dslbaff01 = smoothMPC_ds_aff + 16;
smoothMPC_FLOAT* smoothMPC_dllbcc01 = smoothMPC_dl_cc + 16;
smoothMPC_FLOAT* smoothMPC_dslbcc01 = smoothMPC_ds_cc + 16;
smoothMPC_FLOAT* smoothMPC_ccrhsl01 = smoothMPC_ccrhs + 16;
int smoothMPC_ubIdx01[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub01 = smoothMPC_l + 27;
smoothMPC_FLOAT* smoothMPC_sub01 = smoothMPC_s + 27;
smoothMPC_FLOAT* smoothMPC_lubbysub01 = smoothMPC_lbys + 27;
smoothMPC_FLOAT smoothMPC_riub01[5];
smoothMPC_FLOAT* smoothMPC_dlubaff01 = smoothMPC_dl_aff + 27;
smoothMPC_FLOAT* smoothMPC_dsubaff01 = smoothMPC_ds_aff + 27;
smoothMPC_FLOAT* smoothMPC_dlubcc01 = smoothMPC_dl_cc + 27;
smoothMPC_FLOAT* smoothMPC_dsubcc01 = smoothMPC_ds_cc + 27;
smoothMPC_FLOAT* smoothMPC_ccrhsub01 = smoothMPC_ccrhs + 27;
smoothMPC_FLOAT smoothMPC_Phi01[11];
smoothMPC_FLOAT smoothMPC_D01[11] = {-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000};
smoothMPC_FLOAT smoothMPC_W01[11];
smoothMPC_FLOAT smoothMPC_Ysd01[9];
smoothMPC_FLOAT smoothMPC_Lsd01[9];
smoothMPC_FLOAT* smoothMPC_z02 = smoothMPC_z + 22;
smoothMPC_FLOAT* smoothMPC_dzaff02 = smoothMPC_dz_aff + 22;
smoothMPC_FLOAT* smoothMPC_dzcc02 = smoothMPC_dz_cc + 22;
smoothMPC_FLOAT* smoothMPC_rd02 = smoothMPC_rd + 22;
smoothMPC_FLOAT smoothMPC_Lbyrd02[11];
smoothMPC_FLOAT* smoothMPC_grad_cost02 = smoothMPC_grad_cost + 22;
smoothMPC_FLOAT* smoothMPC_grad_eq02 = smoothMPC_grad_eq + 22;
smoothMPC_FLOAT* smoothMPC_grad_ineq02 = smoothMPC_grad_ineq + 22;
smoothMPC_FLOAT smoothMPC_ctv02[11];
smoothMPC_FLOAT* smoothMPC_v02 = smoothMPC_v + 6;
smoothMPC_FLOAT smoothMPC_re02[3];
smoothMPC_FLOAT smoothMPC_beta02[3];
smoothMPC_FLOAT smoothMPC_betacc02[3];
smoothMPC_FLOAT* smoothMPC_dvaff02 = smoothMPC_dv_aff + 6;
smoothMPC_FLOAT* smoothMPC_dvcc02 = smoothMPC_dv_cc + 6;
smoothMPC_FLOAT smoothMPC_V02[33];
smoothMPC_FLOAT smoothMPC_Yd02[6];
smoothMPC_FLOAT smoothMPC_Ld02[6];
smoothMPC_FLOAT smoothMPC_yy02[3];
smoothMPC_FLOAT smoothMPC_bmy02[3];
int smoothMPC_lbIdx02[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb02 = smoothMPC_l + 32;
smoothMPC_FLOAT* smoothMPC_slb02 = smoothMPC_s + 32;
smoothMPC_FLOAT* smoothMPC_llbbyslb02 = smoothMPC_lbys + 32;
smoothMPC_FLOAT smoothMPC_rilb02[11];
smoothMPC_FLOAT* smoothMPC_dllbaff02 = smoothMPC_dl_aff + 32;
smoothMPC_FLOAT* smoothMPC_dslbaff02 = smoothMPC_ds_aff + 32;
smoothMPC_FLOAT* smoothMPC_dllbcc02 = smoothMPC_dl_cc + 32;
smoothMPC_FLOAT* smoothMPC_dslbcc02 = smoothMPC_ds_cc + 32;
smoothMPC_FLOAT* smoothMPC_ccrhsl02 = smoothMPC_ccrhs + 32;
int smoothMPC_ubIdx02[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub02 = smoothMPC_l + 43;
smoothMPC_FLOAT* smoothMPC_sub02 = smoothMPC_s + 43;
smoothMPC_FLOAT* smoothMPC_lubbysub02 = smoothMPC_lbys + 43;
smoothMPC_FLOAT smoothMPC_riub02[5];
smoothMPC_FLOAT* smoothMPC_dlubaff02 = smoothMPC_dl_aff + 43;
smoothMPC_FLOAT* smoothMPC_dsubaff02 = smoothMPC_ds_aff + 43;
smoothMPC_FLOAT* smoothMPC_dlubcc02 = smoothMPC_dl_cc + 43;
smoothMPC_FLOAT* smoothMPC_dsubcc02 = smoothMPC_ds_cc + 43;
smoothMPC_FLOAT* smoothMPC_ccrhsub02 = smoothMPC_ccrhs + 43;
smoothMPC_FLOAT smoothMPC_Phi02[11];
smoothMPC_FLOAT smoothMPC_W02[11];
smoothMPC_FLOAT smoothMPC_Ysd02[9];
smoothMPC_FLOAT smoothMPC_Lsd02[9];
smoothMPC_FLOAT* smoothMPC_z03 = smoothMPC_z + 33;
smoothMPC_FLOAT* smoothMPC_dzaff03 = smoothMPC_dz_aff + 33;
smoothMPC_FLOAT* smoothMPC_dzcc03 = smoothMPC_dz_cc + 33;
smoothMPC_FLOAT* smoothMPC_rd03 = smoothMPC_rd + 33;
smoothMPC_FLOAT smoothMPC_Lbyrd03[11];
smoothMPC_FLOAT* smoothMPC_grad_cost03 = smoothMPC_grad_cost + 33;
smoothMPC_FLOAT* smoothMPC_grad_eq03 = smoothMPC_grad_eq + 33;
smoothMPC_FLOAT* smoothMPC_grad_ineq03 = smoothMPC_grad_ineq + 33;
smoothMPC_FLOAT smoothMPC_ctv03[11];
smoothMPC_FLOAT* smoothMPC_v03 = smoothMPC_v + 9;
smoothMPC_FLOAT smoothMPC_re03[3];
smoothMPC_FLOAT smoothMPC_beta03[3];
smoothMPC_FLOAT smoothMPC_betacc03[3];
smoothMPC_FLOAT* smoothMPC_dvaff03 = smoothMPC_dv_aff + 9;
smoothMPC_FLOAT* smoothMPC_dvcc03 = smoothMPC_dv_cc + 9;
smoothMPC_FLOAT smoothMPC_V03[33];
smoothMPC_FLOAT smoothMPC_Yd03[6];
smoothMPC_FLOAT smoothMPC_Ld03[6];
smoothMPC_FLOAT smoothMPC_yy03[3];
smoothMPC_FLOAT smoothMPC_bmy03[3];
int smoothMPC_lbIdx03[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb03 = smoothMPC_l + 48;
smoothMPC_FLOAT* smoothMPC_slb03 = smoothMPC_s + 48;
smoothMPC_FLOAT* smoothMPC_llbbyslb03 = smoothMPC_lbys + 48;
smoothMPC_FLOAT smoothMPC_rilb03[11];
smoothMPC_FLOAT* smoothMPC_dllbaff03 = smoothMPC_dl_aff + 48;
smoothMPC_FLOAT* smoothMPC_dslbaff03 = smoothMPC_ds_aff + 48;
smoothMPC_FLOAT* smoothMPC_dllbcc03 = smoothMPC_dl_cc + 48;
smoothMPC_FLOAT* smoothMPC_dslbcc03 = smoothMPC_ds_cc + 48;
smoothMPC_FLOAT* smoothMPC_ccrhsl03 = smoothMPC_ccrhs + 48;
int smoothMPC_ubIdx03[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub03 = smoothMPC_l + 59;
smoothMPC_FLOAT* smoothMPC_sub03 = smoothMPC_s + 59;
smoothMPC_FLOAT* smoothMPC_lubbysub03 = smoothMPC_lbys + 59;
smoothMPC_FLOAT smoothMPC_riub03[5];
smoothMPC_FLOAT* smoothMPC_dlubaff03 = smoothMPC_dl_aff + 59;
smoothMPC_FLOAT* smoothMPC_dsubaff03 = smoothMPC_ds_aff + 59;
smoothMPC_FLOAT* smoothMPC_dlubcc03 = smoothMPC_dl_cc + 59;
smoothMPC_FLOAT* smoothMPC_dsubcc03 = smoothMPC_ds_cc + 59;
smoothMPC_FLOAT* smoothMPC_ccrhsub03 = smoothMPC_ccrhs + 59;
smoothMPC_FLOAT smoothMPC_Phi03[11];
smoothMPC_FLOAT smoothMPC_W03[11];
smoothMPC_FLOAT smoothMPC_Ysd03[9];
smoothMPC_FLOAT smoothMPC_Lsd03[9];
smoothMPC_FLOAT* smoothMPC_z04 = smoothMPC_z + 44;
smoothMPC_FLOAT* smoothMPC_dzaff04 = smoothMPC_dz_aff + 44;
smoothMPC_FLOAT* smoothMPC_dzcc04 = smoothMPC_dz_cc + 44;
smoothMPC_FLOAT* smoothMPC_rd04 = smoothMPC_rd + 44;
smoothMPC_FLOAT smoothMPC_Lbyrd04[11];
smoothMPC_FLOAT* smoothMPC_grad_cost04 = smoothMPC_grad_cost + 44;
smoothMPC_FLOAT* smoothMPC_grad_eq04 = smoothMPC_grad_eq + 44;
smoothMPC_FLOAT* smoothMPC_grad_ineq04 = smoothMPC_grad_ineq + 44;
smoothMPC_FLOAT smoothMPC_ctv04[11];
smoothMPC_FLOAT* smoothMPC_v04 = smoothMPC_v + 12;
smoothMPC_FLOAT smoothMPC_re04[3];
smoothMPC_FLOAT smoothMPC_beta04[3];
smoothMPC_FLOAT smoothMPC_betacc04[3];
smoothMPC_FLOAT* smoothMPC_dvaff04 = smoothMPC_dv_aff + 12;
smoothMPC_FLOAT* smoothMPC_dvcc04 = smoothMPC_dv_cc + 12;
smoothMPC_FLOAT smoothMPC_V04[33];
smoothMPC_FLOAT smoothMPC_Yd04[6];
smoothMPC_FLOAT smoothMPC_Ld04[6];
smoothMPC_FLOAT smoothMPC_yy04[3];
smoothMPC_FLOAT smoothMPC_bmy04[3];
int smoothMPC_lbIdx04[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb04 = smoothMPC_l + 64;
smoothMPC_FLOAT* smoothMPC_slb04 = smoothMPC_s + 64;
smoothMPC_FLOAT* smoothMPC_llbbyslb04 = smoothMPC_lbys + 64;
smoothMPC_FLOAT smoothMPC_rilb04[11];
smoothMPC_FLOAT* smoothMPC_dllbaff04 = smoothMPC_dl_aff + 64;
smoothMPC_FLOAT* smoothMPC_dslbaff04 = smoothMPC_ds_aff + 64;
smoothMPC_FLOAT* smoothMPC_dllbcc04 = smoothMPC_dl_cc + 64;
smoothMPC_FLOAT* smoothMPC_dslbcc04 = smoothMPC_ds_cc + 64;
smoothMPC_FLOAT* smoothMPC_ccrhsl04 = smoothMPC_ccrhs + 64;
int smoothMPC_ubIdx04[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub04 = smoothMPC_l + 75;
smoothMPC_FLOAT* smoothMPC_sub04 = smoothMPC_s + 75;
smoothMPC_FLOAT* smoothMPC_lubbysub04 = smoothMPC_lbys + 75;
smoothMPC_FLOAT smoothMPC_riub04[5];
smoothMPC_FLOAT* smoothMPC_dlubaff04 = smoothMPC_dl_aff + 75;
smoothMPC_FLOAT* smoothMPC_dsubaff04 = smoothMPC_ds_aff + 75;
smoothMPC_FLOAT* smoothMPC_dlubcc04 = smoothMPC_dl_cc + 75;
smoothMPC_FLOAT* smoothMPC_dsubcc04 = smoothMPC_ds_cc + 75;
smoothMPC_FLOAT* smoothMPC_ccrhsub04 = smoothMPC_ccrhs + 75;
smoothMPC_FLOAT smoothMPC_Phi04[11];
smoothMPC_FLOAT smoothMPC_W04[11];
smoothMPC_FLOAT smoothMPC_Ysd04[9];
smoothMPC_FLOAT smoothMPC_Lsd04[9];
smoothMPC_FLOAT* smoothMPC_z05 = smoothMPC_z + 55;
smoothMPC_FLOAT* smoothMPC_dzaff05 = smoothMPC_dz_aff + 55;
smoothMPC_FLOAT* smoothMPC_dzcc05 = smoothMPC_dz_cc + 55;
smoothMPC_FLOAT* smoothMPC_rd05 = smoothMPC_rd + 55;
smoothMPC_FLOAT smoothMPC_Lbyrd05[11];
smoothMPC_FLOAT* smoothMPC_grad_cost05 = smoothMPC_grad_cost + 55;
smoothMPC_FLOAT* smoothMPC_grad_eq05 = smoothMPC_grad_eq + 55;
smoothMPC_FLOAT* smoothMPC_grad_ineq05 = smoothMPC_grad_ineq + 55;
smoothMPC_FLOAT smoothMPC_ctv05[11];
smoothMPC_FLOAT* smoothMPC_v05 = smoothMPC_v + 15;
smoothMPC_FLOAT smoothMPC_re05[3];
smoothMPC_FLOAT smoothMPC_beta05[3];
smoothMPC_FLOAT smoothMPC_betacc05[3];
smoothMPC_FLOAT* smoothMPC_dvaff05 = smoothMPC_dv_aff + 15;
smoothMPC_FLOAT* smoothMPC_dvcc05 = smoothMPC_dv_cc + 15;
smoothMPC_FLOAT smoothMPC_V05[33];
smoothMPC_FLOAT smoothMPC_Yd05[6];
smoothMPC_FLOAT smoothMPC_Ld05[6];
smoothMPC_FLOAT smoothMPC_yy05[3];
smoothMPC_FLOAT smoothMPC_bmy05[3];
int smoothMPC_lbIdx05[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb05 = smoothMPC_l + 80;
smoothMPC_FLOAT* smoothMPC_slb05 = smoothMPC_s + 80;
smoothMPC_FLOAT* smoothMPC_llbbyslb05 = smoothMPC_lbys + 80;
smoothMPC_FLOAT smoothMPC_rilb05[11];
smoothMPC_FLOAT* smoothMPC_dllbaff05 = smoothMPC_dl_aff + 80;
smoothMPC_FLOAT* smoothMPC_dslbaff05 = smoothMPC_ds_aff + 80;
smoothMPC_FLOAT* smoothMPC_dllbcc05 = smoothMPC_dl_cc + 80;
smoothMPC_FLOAT* smoothMPC_dslbcc05 = smoothMPC_ds_cc + 80;
smoothMPC_FLOAT* smoothMPC_ccrhsl05 = smoothMPC_ccrhs + 80;
int smoothMPC_ubIdx05[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub05 = smoothMPC_l + 91;
smoothMPC_FLOAT* smoothMPC_sub05 = smoothMPC_s + 91;
smoothMPC_FLOAT* smoothMPC_lubbysub05 = smoothMPC_lbys + 91;
smoothMPC_FLOAT smoothMPC_riub05[5];
smoothMPC_FLOAT* smoothMPC_dlubaff05 = smoothMPC_dl_aff + 91;
smoothMPC_FLOAT* smoothMPC_dsubaff05 = smoothMPC_ds_aff + 91;
smoothMPC_FLOAT* smoothMPC_dlubcc05 = smoothMPC_dl_cc + 91;
smoothMPC_FLOAT* smoothMPC_dsubcc05 = smoothMPC_ds_cc + 91;
smoothMPC_FLOAT* smoothMPC_ccrhsub05 = smoothMPC_ccrhs + 91;
smoothMPC_FLOAT smoothMPC_Phi05[11];
smoothMPC_FLOAT smoothMPC_W05[11];
smoothMPC_FLOAT smoothMPC_Ysd05[9];
smoothMPC_FLOAT smoothMPC_Lsd05[9];
smoothMPC_FLOAT* smoothMPC_z06 = smoothMPC_z + 66;
smoothMPC_FLOAT* smoothMPC_dzaff06 = smoothMPC_dz_aff + 66;
smoothMPC_FLOAT* smoothMPC_dzcc06 = smoothMPC_dz_cc + 66;
smoothMPC_FLOAT* smoothMPC_rd06 = smoothMPC_rd + 66;
smoothMPC_FLOAT smoothMPC_Lbyrd06[11];
smoothMPC_FLOAT* smoothMPC_grad_cost06 = smoothMPC_grad_cost + 66;
smoothMPC_FLOAT* smoothMPC_grad_eq06 = smoothMPC_grad_eq + 66;
smoothMPC_FLOAT* smoothMPC_grad_ineq06 = smoothMPC_grad_ineq + 66;
smoothMPC_FLOAT smoothMPC_ctv06[11];
smoothMPC_FLOAT* smoothMPC_v06 = smoothMPC_v + 18;
smoothMPC_FLOAT smoothMPC_re06[3];
smoothMPC_FLOAT smoothMPC_beta06[3];
smoothMPC_FLOAT smoothMPC_betacc06[3];
smoothMPC_FLOAT* smoothMPC_dvaff06 = smoothMPC_dv_aff + 18;
smoothMPC_FLOAT* smoothMPC_dvcc06 = smoothMPC_dv_cc + 18;
smoothMPC_FLOAT smoothMPC_V06[33];
smoothMPC_FLOAT smoothMPC_Yd06[6];
smoothMPC_FLOAT smoothMPC_Ld06[6];
smoothMPC_FLOAT smoothMPC_yy06[3];
smoothMPC_FLOAT smoothMPC_bmy06[3];
int smoothMPC_lbIdx06[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb06 = smoothMPC_l + 96;
smoothMPC_FLOAT* smoothMPC_slb06 = smoothMPC_s + 96;
smoothMPC_FLOAT* smoothMPC_llbbyslb06 = smoothMPC_lbys + 96;
smoothMPC_FLOAT smoothMPC_rilb06[11];
smoothMPC_FLOAT* smoothMPC_dllbaff06 = smoothMPC_dl_aff + 96;
smoothMPC_FLOAT* smoothMPC_dslbaff06 = smoothMPC_ds_aff + 96;
smoothMPC_FLOAT* smoothMPC_dllbcc06 = smoothMPC_dl_cc + 96;
smoothMPC_FLOAT* smoothMPC_dslbcc06 = smoothMPC_ds_cc + 96;
smoothMPC_FLOAT* smoothMPC_ccrhsl06 = smoothMPC_ccrhs + 96;
int smoothMPC_ubIdx06[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub06 = smoothMPC_l + 107;
smoothMPC_FLOAT* smoothMPC_sub06 = smoothMPC_s + 107;
smoothMPC_FLOAT* smoothMPC_lubbysub06 = smoothMPC_lbys + 107;
smoothMPC_FLOAT smoothMPC_riub06[5];
smoothMPC_FLOAT* smoothMPC_dlubaff06 = smoothMPC_dl_aff + 107;
smoothMPC_FLOAT* smoothMPC_dsubaff06 = smoothMPC_ds_aff + 107;
smoothMPC_FLOAT* smoothMPC_dlubcc06 = smoothMPC_dl_cc + 107;
smoothMPC_FLOAT* smoothMPC_dsubcc06 = smoothMPC_ds_cc + 107;
smoothMPC_FLOAT* smoothMPC_ccrhsub06 = smoothMPC_ccrhs + 107;
smoothMPC_FLOAT smoothMPC_Phi06[11];
smoothMPC_FLOAT smoothMPC_W06[11];
smoothMPC_FLOAT smoothMPC_Ysd06[9];
smoothMPC_FLOAT smoothMPC_Lsd06[9];
smoothMPC_FLOAT* smoothMPC_z07 = smoothMPC_z + 77;
smoothMPC_FLOAT* smoothMPC_dzaff07 = smoothMPC_dz_aff + 77;
smoothMPC_FLOAT* smoothMPC_dzcc07 = smoothMPC_dz_cc + 77;
smoothMPC_FLOAT* smoothMPC_rd07 = smoothMPC_rd + 77;
smoothMPC_FLOAT smoothMPC_Lbyrd07[11];
smoothMPC_FLOAT* smoothMPC_grad_cost07 = smoothMPC_grad_cost + 77;
smoothMPC_FLOAT* smoothMPC_grad_eq07 = smoothMPC_grad_eq + 77;
smoothMPC_FLOAT* smoothMPC_grad_ineq07 = smoothMPC_grad_ineq + 77;
smoothMPC_FLOAT smoothMPC_ctv07[11];
smoothMPC_FLOAT* smoothMPC_v07 = smoothMPC_v + 21;
smoothMPC_FLOAT smoothMPC_re07[3];
smoothMPC_FLOAT smoothMPC_beta07[3];
smoothMPC_FLOAT smoothMPC_betacc07[3];
smoothMPC_FLOAT* smoothMPC_dvaff07 = smoothMPC_dv_aff + 21;
smoothMPC_FLOAT* smoothMPC_dvcc07 = smoothMPC_dv_cc + 21;
smoothMPC_FLOAT smoothMPC_V07[33];
smoothMPC_FLOAT smoothMPC_Yd07[6];
smoothMPC_FLOAT smoothMPC_Ld07[6];
smoothMPC_FLOAT smoothMPC_yy07[3];
smoothMPC_FLOAT smoothMPC_bmy07[3];
int smoothMPC_lbIdx07[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb07 = smoothMPC_l + 112;
smoothMPC_FLOAT* smoothMPC_slb07 = smoothMPC_s + 112;
smoothMPC_FLOAT* smoothMPC_llbbyslb07 = smoothMPC_lbys + 112;
smoothMPC_FLOAT smoothMPC_rilb07[11];
smoothMPC_FLOAT* smoothMPC_dllbaff07 = smoothMPC_dl_aff + 112;
smoothMPC_FLOAT* smoothMPC_dslbaff07 = smoothMPC_ds_aff + 112;
smoothMPC_FLOAT* smoothMPC_dllbcc07 = smoothMPC_dl_cc + 112;
smoothMPC_FLOAT* smoothMPC_dslbcc07 = smoothMPC_ds_cc + 112;
smoothMPC_FLOAT* smoothMPC_ccrhsl07 = smoothMPC_ccrhs + 112;
int smoothMPC_ubIdx07[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub07 = smoothMPC_l + 123;
smoothMPC_FLOAT* smoothMPC_sub07 = smoothMPC_s + 123;
smoothMPC_FLOAT* smoothMPC_lubbysub07 = smoothMPC_lbys + 123;
smoothMPC_FLOAT smoothMPC_riub07[5];
smoothMPC_FLOAT* smoothMPC_dlubaff07 = smoothMPC_dl_aff + 123;
smoothMPC_FLOAT* smoothMPC_dsubaff07 = smoothMPC_ds_aff + 123;
smoothMPC_FLOAT* smoothMPC_dlubcc07 = smoothMPC_dl_cc + 123;
smoothMPC_FLOAT* smoothMPC_dsubcc07 = smoothMPC_ds_cc + 123;
smoothMPC_FLOAT* smoothMPC_ccrhsub07 = smoothMPC_ccrhs + 123;
smoothMPC_FLOAT smoothMPC_Phi07[11];
smoothMPC_FLOAT smoothMPC_W07[11];
smoothMPC_FLOAT smoothMPC_Ysd07[9];
smoothMPC_FLOAT smoothMPC_Lsd07[9];
smoothMPC_FLOAT* smoothMPC_z08 = smoothMPC_z + 88;
smoothMPC_FLOAT* smoothMPC_dzaff08 = smoothMPC_dz_aff + 88;
smoothMPC_FLOAT* smoothMPC_dzcc08 = smoothMPC_dz_cc + 88;
smoothMPC_FLOAT* smoothMPC_rd08 = smoothMPC_rd + 88;
smoothMPC_FLOAT smoothMPC_Lbyrd08[11];
smoothMPC_FLOAT* smoothMPC_grad_cost08 = smoothMPC_grad_cost + 88;
smoothMPC_FLOAT* smoothMPC_grad_eq08 = smoothMPC_grad_eq + 88;
smoothMPC_FLOAT* smoothMPC_grad_ineq08 = smoothMPC_grad_ineq + 88;
smoothMPC_FLOAT smoothMPC_ctv08[11];
smoothMPC_FLOAT* smoothMPC_v08 = smoothMPC_v + 24;
smoothMPC_FLOAT smoothMPC_re08[3];
smoothMPC_FLOAT smoothMPC_beta08[3];
smoothMPC_FLOAT smoothMPC_betacc08[3];
smoothMPC_FLOAT* smoothMPC_dvaff08 = smoothMPC_dv_aff + 24;
smoothMPC_FLOAT* smoothMPC_dvcc08 = smoothMPC_dv_cc + 24;
smoothMPC_FLOAT smoothMPC_V08[33];
smoothMPC_FLOAT smoothMPC_Yd08[6];
smoothMPC_FLOAT smoothMPC_Ld08[6];
smoothMPC_FLOAT smoothMPC_yy08[3];
smoothMPC_FLOAT smoothMPC_bmy08[3];
int smoothMPC_lbIdx08[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb08 = smoothMPC_l + 128;
smoothMPC_FLOAT* smoothMPC_slb08 = smoothMPC_s + 128;
smoothMPC_FLOAT* smoothMPC_llbbyslb08 = smoothMPC_lbys + 128;
smoothMPC_FLOAT smoothMPC_rilb08[11];
smoothMPC_FLOAT* smoothMPC_dllbaff08 = smoothMPC_dl_aff + 128;
smoothMPC_FLOAT* smoothMPC_dslbaff08 = smoothMPC_ds_aff + 128;
smoothMPC_FLOAT* smoothMPC_dllbcc08 = smoothMPC_dl_cc + 128;
smoothMPC_FLOAT* smoothMPC_dslbcc08 = smoothMPC_ds_cc + 128;
smoothMPC_FLOAT* smoothMPC_ccrhsl08 = smoothMPC_ccrhs + 128;
int smoothMPC_ubIdx08[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub08 = smoothMPC_l + 139;
smoothMPC_FLOAT* smoothMPC_sub08 = smoothMPC_s + 139;
smoothMPC_FLOAT* smoothMPC_lubbysub08 = smoothMPC_lbys + 139;
smoothMPC_FLOAT smoothMPC_riub08[5];
smoothMPC_FLOAT* smoothMPC_dlubaff08 = smoothMPC_dl_aff + 139;
smoothMPC_FLOAT* smoothMPC_dsubaff08 = smoothMPC_ds_aff + 139;
smoothMPC_FLOAT* smoothMPC_dlubcc08 = smoothMPC_dl_cc + 139;
smoothMPC_FLOAT* smoothMPC_dsubcc08 = smoothMPC_ds_cc + 139;
smoothMPC_FLOAT* smoothMPC_ccrhsub08 = smoothMPC_ccrhs + 139;
smoothMPC_FLOAT smoothMPC_Phi08[11];
smoothMPC_FLOAT smoothMPC_W08[11];
smoothMPC_FLOAT smoothMPC_Ysd08[9];
smoothMPC_FLOAT smoothMPC_Lsd08[9];
smoothMPC_FLOAT* smoothMPC_z09 = smoothMPC_z + 99;
smoothMPC_FLOAT* smoothMPC_dzaff09 = smoothMPC_dz_aff + 99;
smoothMPC_FLOAT* smoothMPC_dzcc09 = smoothMPC_dz_cc + 99;
smoothMPC_FLOAT* smoothMPC_rd09 = smoothMPC_rd + 99;
smoothMPC_FLOAT smoothMPC_Lbyrd09[11];
smoothMPC_FLOAT* smoothMPC_grad_cost09 = smoothMPC_grad_cost + 99;
smoothMPC_FLOAT* smoothMPC_grad_eq09 = smoothMPC_grad_eq + 99;
smoothMPC_FLOAT* smoothMPC_grad_ineq09 = smoothMPC_grad_ineq + 99;
smoothMPC_FLOAT smoothMPC_ctv09[11];
smoothMPC_FLOAT* smoothMPC_v09 = smoothMPC_v + 27;
smoothMPC_FLOAT smoothMPC_re09[3];
smoothMPC_FLOAT smoothMPC_beta09[3];
smoothMPC_FLOAT smoothMPC_betacc09[3];
smoothMPC_FLOAT* smoothMPC_dvaff09 = smoothMPC_dv_aff + 27;
smoothMPC_FLOAT* smoothMPC_dvcc09 = smoothMPC_dv_cc + 27;
smoothMPC_FLOAT smoothMPC_V09[33];
smoothMPC_FLOAT smoothMPC_Yd09[6];
smoothMPC_FLOAT smoothMPC_Ld09[6];
smoothMPC_FLOAT smoothMPC_yy09[3];
smoothMPC_FLOAT smoothMPC_bmy09[3];
int smoothMPC_lbIdx09[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb09 = smoothMPC_l + 144;
smoothMPC_FLOAT* smoothMPC_slb09 = smoothMPC_s + 144;
smoothMPC_FLOAT* smoothMPC_llbbyslb09 = smoothMPC_lbys + 144;
smoothMPC_FLOAT smoothMPC_rilb09[11];
smoothMPC_FLOAT* smoothMPC_dllbaff09 = smoothMPC_dl_aff + 144;
smoothMPC_FLOAT* smoothMPC_dslbaff09 = smoothMPC_ds_aff + 144;
smoothMPC_FLOAT* smoothMPC_dllbcc09 = smoothMPC_dl_cc + 144;
smoothMPC_FLOAT* smoothMPC_dslbcc09 = smoothMPC_ds_cc + 144;
smoothMPC_FLOAT* smoothMPC_ccrhsl09 = smoothMPC_ccrhs + 144;
int smoothMPC_ubIdx09[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub09 = smoothMPC_l + 155;
smoothMPC_FLOAT* smoothMPC_sub09 = smoothMPC_s + 155;
smoothMPC_FLOAT* smoothMPC_lubbysub09 = smoothMPC_lbys + 155;
smoothMPC_FLOAT smoothMPC_riub09[5];
smoothMPC_FLOAT* smoothMPC_dlubaff09 = smoothMPC_dl_aff + 155;
smoothMPC_FLOAT* smoothMPC_dsubaff09 = smoothMPC_ds_aff + 155;
smoothMPC_FLOAT* smoothMPC_dlubcc09 = smoothMPC_dl_cc + 155;
smoothMPC_FLOAT* smoothMPC_dsubcc09 = smoothMPC_ds_cc + 155;
smoothMPC_FLOAT* smoothMPC_ccrhsub09 = smoothMPC_ccrhs + 155;
smoothMPC_FLOAT smoothMPC_Phi09[11];
smoothMPC_FLOAT smoothMPC_W09[11];
smoothMPC_FLOAT smoothMPC_Ysd09[9];
smoothMPC_FLOAT smoothMPC_Lsd09[9];
smoothMPC_FLOAT* smoothMPC_z10 = smoothMPC_z + 110;
smoothMPC_FLOAT* smoothMPC_dzaff10 = smoothMPC_dz_aff + 110;
smoothMPC_FLOAT* smoothMPC_dzcc10 = smoothMPC_dz_cc + 110;
smoothMPC_FLOAT* smoothMPC_rd10 = smoothMPC_rd + 110;
smoothMPC_FLOAT smoothMPC_Lbyrd10[11];
smoothMPC_FLOAT* smoothMPC_grad_cost10 = smoothMPC_grad_cost + 110;
smoothMPC_FLOAT* smoothMPC_grad_eq10 = smoothMPC_grad_eq + 110;
smoothMPC_FLOAT* smoothMPC_grad_ineq10 = smoothMPC_grad_ineq + 110;
smoothMPC_FLOAT smoothMPC_ctv10[11];
smoothMPC_FLOAT* smoothMPC_v10 = smoothMPC_v + 30;
smoothMPC_FLOAT smoothMPC_re10[3];
smoothMPC_FLOAT smoothMPC_beta10[3];
smoothMPC_FLOAT smoothMPC_betacc10[3];
smoothMPC_FLOAT* smoothMPC_dvaff10 = smoothMPC_dv_aff + 30;
smoothMPC_FLOAT* smoothMPC_dvcc10 = smoothMPC_dv_cc + 30;
smoothMPC_FLOAT smoothMPC_V10[33];
smoothMPC_FLOAT smoothMPC_Yd10[6];
smoothMPC_FLOAT smoothMPC_Ld10[6];
smoothMPC_FLOAT smoothMPC_yy10[3];
smoothMPC_FLOAT smoothMPC_bmy10[3];
int smoothMPC_lbIdx10[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb10 = smoothMPC_l + 160;
smoothMPC_FLOAT* smoothMPC_slb10 = smoothMPC_s + 160;
smoothMPC_FLOAT* smoothMPC_llbbyslb10 = smoothMPC_lbys + 160;
smoothMPC_FLOAT smoothMPC_rilb10[11];
smoothMPC_FLOAT* smoothMPC_dllbaff10 = smoothMPC_dl_aff + 160;
smoothMPC_FLOAT* smoothMPC_dslbaff10 = smoothMPC_ds_aff + 160;
smoothMPC_FLOAT* smoothMPC_dllbcc10 = smoothMPC_dl_cc + 160;
smoothMPC_FLOAT* smoothMPC_dslbcc10 = smoothMPC_ds_cc + 160;
smoothMPC_FLOAT* smoothMPC_ccrhsl10 = smoothMPC_ccrhs + 160;
int smoothMPC_ubIdx10[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub10 = smoothMPC_l + 171;
smoothMPC_FLOAT* smoothMPC_sub10 = smoothMPC_s + 171;
smoothMPC_FLOAT* smoothMPC_lubbysub10 = smoothMPC_lbys + 171;
smoothMPC_FLOAT smoothMPC_riub10[5];
smoothMPC_FLOAT* smoothMPC_dlubaff10 = smoothMPC_dl_aff + 171;
smoothMPC_FLOAT* smoothMPC_dsubaff10 = smoothMPC_ds_aff + 171;
smoothMPC_FLOAT* smoothMPC_dlubcc10 = smoothMPC_dl_cc + 171;
smoothMPC_FLOAT* smoothMPC_dsubcc10 = smoothMPC_ds_cc + 171;
smoothMPC_FLOAT* smoothMPC_ccrhsub10 = smoothMPC_ccrhs + 171;
smoothMPC_FLOAT smoothMPC_Phi10[11];
smoothMPC_FLOAT smoothMPC_W10[11];
smoothMPC_FLOAT smoothMPC_Ysd10[9];
smoothMPC_FLOAT smoothMPC_Lsd10[9];
smoothMPC_FLOAT* smoothMPC_z11 = smoothMPC_z + 121;
smoothMPC_FLOAT* smoothMPC_dzaff11 = smoothMPC_dz_aff + 121;
smoothMPC_FLOAT* smoothMPC_dzcc11 = smoothMPC_dz_cc + 121;
smoothMPC_FLOAT* smoothMPC_rd11 = smoothMPC_rd + 121;
smoothMPC_FLOAT smoothMPC_Lbyrd11[11];
smoothMPC_FLOAT* smoothMPC_grad_cost11 = smoothMPC_grad_cost + 121;
smoothMPC_FLOAT* smoothMPC_grad_eq11 = smoothMPC_grad_eq + 121;
smoothMPC_FLOAT* smoothMPC_grad_ineq11 = smoothMPC_grad_ineq + 121;
smoothMPC_FLOAT smoothMPC_ctv11[11];
smoothMPC_FLOAT* smoothMPC_v11 = smoothMPC_v + 33;
smoothMPC_FLOAT smoothMPC_re11[3];
smoothMPC_FLOAT smoothMPC_beta11[3];
smoothMPC_FLOAT smoothMPC_betacc11[3];
smoothMPC_FLOAT* smoothMPC_dvaff11 = smoothMPC_dv_aff + 33;
smoothMPC_FLOAT* smoothMPC_dvcc11 = smoothMPC_dv_cc + 33;
smoothMPC_FLOAT smoothMPC_V11[33];
smoothMPC_FLOAT smoothMPC_Yd11[6];
smoothMPC_FLOAT smoothMPC_Ld11[6];
smoothMPC_FLOAT smoothMPC_yy11[3];
smoothMPC_FLOAT smoothMPC_bmy11[3];
int smoothMPC_lbIdx11[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb11 = smoothMPC_l + 176;
smoothMPC_FLOAT* smoothMPC_slb11 = smoothMPC_s + 176;
smoothMPC_FLOAT* smoothMPC_llbbyslb11 = smoothMPC_lbys + 176;
smoothMPC_FLOAT smoothMPC_rilb11[11];
smoothMPC_FLOAT* smoothMPC_dllbaff11 = smoothMPC_dl_aff + 176;
smoothMPC_FLOAT* smoothMPC_dslbaff11 = smoothMPC_ds_aff + 176;
smoothMPC_FLOAT* smoothMPC_dllbcc11 = smoothMPC_dl_cc + 176;
smoothMPC_FLOAT* smoothMPC_dslbcc11 = smoothMPC_ds_cc + 176;
smoothMPC_FLOAT* smoothMPC_ccrhsl11 = smoothMPC_ccrhs + 176;
int smoothMPC_ubIdx11[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub11 = smoothMPC_l + 187;
smoothMPC_FLOAT* smoothMPC_sub11 = smoothMPC_s + 187;
smoothMPC_FLOAT* smoothMPC_lubbysub11 = smoothMPC_lbys + 187;
smoothMPC_FLOAT smoothMPC_riub11[5];
smoothMPC_FLOAT* smoothMPC_dlubaff11 = smoothMPC_dl_aff + 187;
smoothMPC_FLOAT* smoothMPC_dsubaff11 = smoothMPC_ds_aff + 187;
smoothMPC_FLOAT* smoothMPC_dlubcc11 = smoothMPC_dl_cc + 187;
smoothMPC_FLOAT* smoothMPC_dsubcc11 = smoothMPC_ds_cc + 187;
smoothMPC_FLOAT* smoothMPC_ccrhsub11 = smoothMPC_ccrhs + 187;
smoothMPC_FLOAT smoothMPC_Phi11[11];
smoothMPC_FLOAT smoothMPC_W11[11];
smoothMPC_FLOAT smoothMPC_Ysd11[9];
smoothMPC_FLOAT smoothMPC_Lsd11[9];
smoothMPC_FLOAT* smoothMPC_z12 = smoothMPC_z + 132;
smoothMPC_FLOAT* smoothMPC_dzaff12 = smoothMPC_dz_aff + 132;
smoothMPC_FLOAT* smoothMPC_dzcc12 = smoothMPC_dz_cc + 132;
smoothMPC_FLOAT* smoothMPC_rd12 = smoothMPC_rd + 132;
smoothMPC_FLOAT smoothMPC_Lbyrd12[11];
smoothMPC_FLOAT* smoothMPC_grad_cost12 = smoothMPC_grad_cost + 132;
smoothMPC_FLOAT* smoothMPC_grad_eq12 = smoothMPC_grad_eq + 132;
smoothMPC_FLOAT* smoothMPC_grad_ineq12 = smoothMPC_grad_ineq + 132;
smoothMPC_FLOAT smoothMPC_ctv12[11];
smoothMPC_FLOAT* smoothMPC_v12 = smoothMPC_v + 36;
smoothMPC_FLOAT smoothMPC_re12[3];
smoothMPC_FLOAT smoothMPC_beta12[3];
smoothMPC_FLOAT smoothMPC_betacc12[3];
smoothMPC_FLOAT* smoothMPC_dvaff12 = smoothMPC_dv_aff + 36;
smoothMPC_FLOAT* smoothMPC_dvcc12 = smoothMPC_dv_cc + 36;
smoothMPC_FLOAT smoothMPC_V12[33];
smoothMPC_FLOAT smoothMPC_Yd12[6];
smoothMPC_FLOAT smoothMPC_Ld12[6];
smoothMPC_FLOAT smoothMPC_yy12[3];
smoothMPC_FLOAT smoothMPC_bmy12[3];
int smoothMPC_lbIdx12[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb12 = smoothMPC_l + 192;
smoothMPC_FLOAT* smoothMPC_slb12 = smoothMPC_s + 192;
smoothMPC_FLOAT* smoothMPC_llbbyslb12 = smoothMPC_lbys + 192;
smoothMPC_FLOAT smoothMPC_rilb12[11];
smoothMPC_FLOAT* smoothMPC_dllbaff12 = smoothMPC_dl_aff + 192;
smoothMPC_FLOAT* smoothMPC_dslbaff12 = smoothMPC_ds_aff + 192;
smoothMPC_FLOAT* smoothMPC_dllbcc12 = smoothMPC_dl_cc + 192;
smoothMPC_FLOAT* smoothMPC_dslbcc12 = smoothMPC_ds_cc + 192;
smoothMPC_FLOAT* smoothMPC_ccrhsl12 = smoothMPC_ccrhs + 192;
int smoothMPC_ubIdx12[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub12 = smoothMPC_l + 203;
smoothMPC_FLOAT* smoothMPC_sub12 = smoothMPC_s + 203;
smoothMPC_FLOAT* smoothMPC_lubbysub12 = smoothMPC_lbys + 203;
smoothMPC_FLOAT smoothMPC_riub12[5];
smoothMPC_FLOAT* smoothMPC_dlubaff12 = smoothMPC_dl_aff + 203;
smoothMPC_FLOAT* smoothMPC_dsubaff12 = smoothMPC_ds_aff + 203;
smoothMPC_FLOAT* smoothMPC_dlubcc12 = smoothMPC_dl_cc + 203;
smoothMPC_FLOAT* smoothMPC_dsubcc12 = smoothMPC_ds_cc + 203;
smoothMPC_FLOAT* smoothMPC_ccrhsub12 = smoothMPC_ccrhs + 203;
smoothMPC_FLOAT smoothMPC_Phi12[11];
smoothMPC_FLOAT smoothMPC_W12[11];
smoothMPC_FLOAT smoothMPC_Ysd12[9];
smoothMPC_FLOAT smoothMPC_Lsd12[9];
smoothMPC_FLOAT* smoothMPC_z13 = smoothMPC_z + 143;
smoothMPC_FLOAT* smoothMPC_dzaff13 = smoothMPC_dz_aff + 143;
smoothMPC_FLOAT* smoothMPC_dzcc13 = smoothMPC_dz_cc + 143;
smoothMPC_FLOAT* smoothMPC_rd13 = smoothMPC_rd + 143;
smoothMPC_FLOAT smoothMPC_Lbyrd13[11];
smoothMPC_FLOAT* smoothMPC_grad_cost13 = smoothMPC_grad_cost + 143;
smoothMPC_FLOAT* smoothMPC_grad_eq13 = smoothMPC_grad_eq + 143;
smoothMPC_FLOAT* smoothMPC_grad_ineq13 = smoothMPC_grad_ineq + 143;
smoothMPC_FLOAT smoothMPC_ctv13[11];
smoothMPC_FLOAT* smoothMPC_v13 = smoothMPC_v + 39;
smoothMPC_FLOAT smoothMPC_re13[3];
smoothMPC_FLOAT smoothMPC_beta13[3];
smoothMPC_FLOAT smoothMPC_betacc13[3];
smoothMPC_FLOAT* smoothMPC_dvaff13 = smoothMPC_dv_aff + 39;
smoothMPC_FLOAT* smoothMPC_dvcc13 = smoothMPC_dv_cc + 39;
smoothMPC_FLOAT smoothMPC_V13[33];
smoothMPC_FLOAT smoothMPC_Yd13[6];
smoothMPC_FLOAT smoothMPC_Ld13[6];
smoothMPC_FLOAT smoothMPC_yy13[3];
smoothMPC_FLOAT smoothMPC_bmy13[3];
int smoothMPC_lbIdx13[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb13 = smoothMPC_l + 208;
smoothMPC_FLOAT* smoothMPC_slb13 = smoothMPC_s + 208;
smoothMPC_FLOAT* smoothMPC_llbbyslb13 = smoothMPC_lbys + 208;
smoothMPC_FLOAT smoothMPC_rilb13[11];
smoothMPC_FLOAT* smoothMPC_dllbaff13 = smoothMPC_dl_aff + 208;
smoothMPC_FLOAT* smoothMPC_dslbaff13 = smoothMPC_ds_aff + 208;
smoothMPC_FLOAT* smoothMPC_dllbcc13 = smoothMPC_dl_cc + 208;
smoothMPC_FLOAT* smoothMPC_dslbcc13 = smoothMPC_ds_cc + 208;
smoothMPC_FLOAT* smoothMPC_ccrhsl13 = smoothMPC_ccrhs + 208;
int smoothMPC_ubIdx13[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub13 = smoothMPC_l + 219;
smoothMPC_FLOAT* smoothMPC_sub13 = smoothMPC_s + 219;
smoothMPC_FLOAT* smoothMPC_lubbysub13 = smoothMPC_lbys + 219;
smoothMPC_FLOAT smoothMPC_riub13[5];
smoothMPC_FLOAT* smoothMPC_dlubaff13 = smoothMPC_dl_aff + 219;
smoothMPC_FLOAT* smoothMPC_dsubaff13 = smoothMPC_ds_aff + 219;
smoothMPC_FLOAT* smoothMPC_dlubcc13 = smoothMPC_dl_cc + 219;
smoothMPC_FLOAT* smoothMPC_dsubcc13 = smoothMPC_ds_cc + 219;
smoothMPC_FLOAT* smoothMPC_ccrhsub13 = smoothMPC_ccrhs + 219;
smoothMPC_FLOAT smoothMPC_Phi13[11];
smoothMPC_FLOAT smoothMPC_W13[11];
smoothMPC_FLOAT smoothMPC_Ysd13[9];
smoothMPC_FLOAT smoothMPC_Lsd13[9];
smoothMPC_FLOAT* smoothMPC_z14 = smoothMPC_z + 154;
smoothMPC_FLOAT* smoothMPC_dzaff14 = smoothMPC_dz_aff + 154;
smoothMPC_FLOAT* smoothMPC_dzcc14 = smoothMPC_dz_cc + 154;
smoothMPC_FLOAT* smoothMPC_rd14 = smoothMPC_rd + 154;
smoothMPC_FLOAT smoothMPC_Lbyrd14[11];
smoothMPC_FLOAT* smoothMPC_grad_cost14 = smoothMPC_grad_cost + 154;
smoothMPC_FLOAT* smoothMPC_grad_eq14 = smoothMPC_grad_eq + 154;
smoothMPC_FLOAT* smoothMPC_grad_ineq14 = smoothMPC_grad_ineq + 154;
smoothMPC_FLOAT smoothMPC_ctv14[11];
smoothMPC_FLOAT* smoothMPC_v14 = smoothMPC_v + 42;
smoothMPC_FLOAT smoothMPC_re14[3];
smoothMPC_FLOAT smoothMPC_beta14[3];
smoothMPC_FLOAT smoothMPC_betacc14[3];
smoothMPC_FLOAT* smoothMPC_dvaff14 = smoothMPC_dv_aff + 42;
smoothMPC_FLOAT* smoothMPC_dvcc14 = smoothMPC_dv_cc + 42;
smoothMPC_FLOAT smoothMPC_V14[33];
smoothMPC_FLOAT smoothMPC_Yd14[6];
smoothMPC_FLOAT smoothMPC_Ld14[6];
smoothMPC_FLOAT smoothMPC_yy14[3];
smoothMPC_FLOAT smoothMPC_bmy14[3];
int smoothMPC_lbIdx14[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb14 = smoothMPC_l + 224;
smoothMPC_FLOAT* smoothMPC_slb14 = smoothMPC_s + 224;
smoothMPC_FLOAT* smoothMPC_llbbyslb14 = smoothMPC_lbys + 224;
smoothMPC_FLOAT smoothMPC_rilb14[11];
smoothMPC_FLOAT* smoothMPC_dllbaff14 = smoothMPC_dl_aff + 224;
smoothMPC_FLOAT* smoothMPC_dslbaff14 = smoothMPC_ds_aff + 224;
smoothMPC_FLOAT* smoothMPC_dllbcc14 = smoothMPC_dl_cc + 224;
smoothMPC_FLOAT* smoothMPC_dslbcc14 = smoothMPC_ds_cc + 224;
smoothMPC_FLOAT* smoothMPC_ccrhsl14 = smoothMPC_ccrhs + 224;
int smoothMPC_ubIdx14[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub14 = smoothMPC_l + 235;
smoothMPC_FLOAT* smoothMPC_sub14 = smoothMPC_s + 235;
smoothMPC_FLOAT* smoothMPC_lubbysub14 = smoothMPC_lbys + 235;
smoothMPC_FLOAT smoothMPC_riub14[5];
smoothMPC_FLOAT* smoothMPC_dlubaff14 = smoothMPC_dl_aff + 235;
smoothMPC_FLOAT* smoothMPC_dsubaff14 = smoothMPC_ds_aff + 235;
smoothMPC_FLOAT* smoothMPC_dlubcc14 = smoothMPC_dl_cc + 235;
smoothMPC_FLOAT* smoothMPC_dsubcc14 = smoothMPC_ds_cc + 235;
smoothMPC_FLOAT* smoothMPC_ccrhsub14 = smoothMPC_ccrhs + 235;
smoothMPC_FLOAT smoothMPC_Phi14[11];
smoothMPC_FLOAT smoothMPC_W14[11];
smoothMPC_FLOAT smoothMPC_Ysd14[9];
smoothMPC_FLOAT smoothMPC_Lsd14[9];
smoothMPC_FLOAT* smoothMPC_z15 = smoothMPC_z + 165;
smoothMPC_FLOAT* smoothMPC_dzaff15 = smoothMPC_dz_aff + 165;
smoothMPC_FLOAT* smoothMPC_dzcc15 = smoothMPC_dz_cc + 165;
smoothMPC_FLOAT* smoothMPC_rd15 = smoothMPC_rd + 165;
smoothMPC_FLOAT smoothMPC_Lbyrd15[11];
smoothMPC_FLOAT* smoothMPC_grad_cost15 = smoothMPC_grad_cost + 165;
smoothMPC_FLOAT* smoothMPC_grad_eq15 = smoothMPC_grad_eq + 165;
smoothMPC_FLOAT* smoothMPC_grad_ineq15 = smoothMPC_grad_ineq + 165;
smoothMPC_FLOAT smoothMPC_ctv15[11];
smoothMPC_FLOAT* smoothMPC_v15 = smoothMPC_v + 45;
smoothMPC_FLOAT smoothMPC_re15[3];
smoothMPC_FLOAT smoothMPC_beta15[3];
smoothMPC_FLOAT smoothMPC_betacc15[3];
smoothMPC_FLOAT* smoothMPC_dvaff15 = smoothMPC_dv_aff + 45;
smoothMPC_FLOAT* smoothMPC_dvcc15 = smoothMPC_dv_cc + 45;
smoothMPC_FLOAT smoothMPC_V15[33];
smoothMPC_FLOAT smoothMPC_Yd15[6];
smoothMPC_FLOAT smoothMPC_Ld15[6];
smoothMPC_FLOAT smoothMPC_yy15[3];
smoothMPC_FLOAT smoothMPC_bmy15[3];
int smoothMPC_lbIdx15[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb15 = smoothMPC_l + 240;
smoothMPC_FLOAT* smoothMPC_slb15 = smoothMPC_s + 240;
smoothMPC_FLOAT* smoothMPC_llbbyslb15 = smoothMPC_lbys + 240;
smoothMPC_FLOAT smoothMPC_rilb15[11];
smoothMPC_FLOAT* smoothMPC_dllbaff15 = smoothMPC_dl_aff + 240;
smoothMPC_FLOAT* smoothMPC_dslbaff15 = smoothMPC_ds_aff + 240;
smoothMPC_FLOAT* smoothMPC_dllbcc15 = smoothMPC_dl_cc + 240;
smoothMPC_FLOAT* smoothMPC_dslbcc15 = smoothMPC_ds_cc + 240;
smoothMPC_FLOAT* smoothMPC_ccrhsl15 = smoothMPC_ccrhs + 240;
int smoothMPC_ubIdx15[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub15 = smoothMPC_l + 251;
smoothMPC_FLOAT* smoothMPC_sub15 = smoothMPC_s + 251;
smoothMPC_FLOAT* smoothMPC_lubbysub15 = smoothMPC_lbys + 251;
smoothMPC_FLOAT smoothMPC_riub15[5];
smoothMPC_FLOAT* smoothMPC_dlubaff15 = smoothMPC_dl_aff + 251;
smoothMPC_FLOAT* smoothMPC_dsubaff15 = smoothMPC_ds_aff + 251;
smoothMPC_FLOAT* smoothMPC_dlubcc15 = smoothMPC_dl_cc + 251;
smoothMPC_FLOAT* smoothMPC_dsubcc15 = smoothMPC_ds_cc + 251;
smoothMPC_FLOAT* smoothMPC_ccrhsub15 = smoothMPC_ccrhs + 251;
smoothMPC_FLOAT smoothMPC_Phi15[11];
smoothMPC_FLOAT smoothMPC_W15[11];
smoothMPC_FLOAT smoothMPC_Ysd15[9];
smoothMPC_FLOAT smoothMPC_Lsd15[9];
smoothMPC_FLOAT* smoothMPC_z16 = smoothMPC_z + 176;
smoothMPC_FLOAT* smoothMPC_dzaff16 = smoothMPC_dz_aff + 176;
smoothMPC_FLOAT* smoothMPC_dzcc16 = smoothMPC_dz_cc + 176;
smoothMPC_FLOAT* smoothMPC_rd16 = smoothMPC_rd + 176;
smoothMPC_FLOAT smoothMPC_Lbyrd16[11];
smoothMPC_FLOAT* smoothMPC_grad_cost16 = smoothMPC_grad_cost + 176;
smoothMPC_FLOAT* smoothMPC_grad_eq16 = smoothMPC_grad_eq + 176;
smoothMPC_FLOAT* smoothMPC_grad_ineq16 = smoothMPC_grad_ineq + 176;
smoothMPC_FLOAT smoothMPC_ctv16[11];
smoothMPC_FLOAT* smoothMPC_v16 = smoothMPC_v + 48;
smoothMPC_FLOAT smoothMPC_re16[3];
smoothMPC_FLOAT smoothMPC_beta16[3];
smoothMPC_FLOAT smoothMPC_betacc16[3];
smoothMPC_FLOAT* smoothMPC_dvaff16 = smoothMPC_dv_aff + 48;
smoothMPC_FLOAT* smoothMPC_dvcc16 = smoothMPC_dv_cc + 48;
smoothMPC_FLOAT smoothMPC_V16[33];
smoothMPC_FLOAT smoothMPC_Yd16[6];
smoothMPC_FLOAT smoothMPC_Ld16[6];
smoothMPC_FLOAT smoothMPC_yy16[3];
smoothMPC_FLOAT smoothMPC_bmy16[3];
int smoothMPC_lbIdx16[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb16 = smoothMPC_l + 256;
smoothMPC_FLOAT* smoothMPC_slb16 = smoothMPC_s + 256;
smoothMPC_FLOAT* smoothMPC_llbbyslb16 = smoothMPC_lbys + 256;
smoothMPC_FLOAT smoothMPC_rilb16[11];
smoothMPC_FLOAT* smoothMPC_dllbaff16 = smoothMPC_dl_aff + 256;
smoothMPC_FLOAT* smoothMPC_dslbaff16 = smoothMPC_ds_aff + 256;
smoothMPC_FLOAT* smoothMPC_dllbcc16 = smoothMPC_dl_cc + 256;
smoothMPC_FLOAT* smoothMPC_dslbcc16 = smoothMPC_ds_cc + 256;
smoothMPC_FLOAT* smoothMPC_ccrhsl16 = smoothMPC_ccrhs + 256;
int smoothMPC_ubIdx16[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub16 = smoothMPC_l + 267;
smoothMPC_FLOAT* smoothMPC_sub16 = smoothMPC_s + 267;
smoothMPC_FLOAT* smoothMPC_lubbysub16 = smoothMPC_lbys + 267;
smoothMPC_FLOAT smoothMPC_riub16[5];
smoothMPC_FLOAT* smoothMPC_dlubaff16 = smoothMPC_dl_aff + 267;
smoothMPC_FLOAT* smoothMPC_dsubaff16 = smoothMPC_ds_aff + 267;
smoothMPC_FLOAT* smoothMPC_dlubcc16 = smoothMPC_dl_cc + 267;
smoothMPC_FLOAT* smoothMPC_dsubcc16 = smoothMPC_ds_cc + 267;
smoothMPC_FLOAT* smoothMPC_ccrhsub16 = smoothMPC_ccrhs + 267;
smoothMPC_FLOAT smoothMPC_Phi16[11];
smoothMPC_FLOAT smoothMPC_W16[11];
smoothMPC_FLOAT smoothMPC_Ysd16[9];
smoothMPC_FLOAT smoothMPC_Lsd16[9];
smoothMPC_FLOAT* smoothMPC_z17 = smoothMPC_z + 187;
smoothMPC_FLOAT* smoothMPC_dzaff17 = smoothMPC_dz_aff + 187;
smoothMPC_FLOAT* smoothMPC_dzcc17 = smoothMPC_dz_cc + 187;
smoothMPC_FLOAT* smoothMPC_rd17 = smoothMPC_rd + 187;
smoothMPC_FLOAT smoothMPC_Lbyrd17[11];
smoothMPC_FLOAT* smoothMPC_grad_cost17 = smoothMPC_grad_cost + 187;
smoothMPC_FLOAT* smoothMPC_grad_eq17 = smoothMPC_grad_eq + 187;
smoothMPC_FLOAT* smoothMPC_grad_ineq17 = smoothMPC_grad_ineq + 187;
smoothMPC_FLOAT smoothMPC_ctv17[11];
smoothMPC_FLOAT* smoothMPC_v17 = smoothMPC_v + 51;
smoothMPC_FLOAT smoothMPC_re17[3];
smoothMPC_FLOAT smoothMPC_beta17[3];
smoothMPC_FLOAT smoothMPC_betacc17[3];
smoothMPC_FLOAT* smoothMPC_dvaff17 = smoothMPC_dv_aff + 51;
smoothMPC_FLOAT* smoothMPC_dvcc17 = smoothMPC_dv_cc + 51;
smoothMPC_FLOAT smoothMPC_V17[33];
smoothMPC_FLOAT smoothMPC_Yd17[6];
smoothMPC_FLOAT smoothMPC_Ld17[6];
smoothMPC_FLOAT smoothMPC_yy17[3];
smoothMPC_FLOAT smoothMPC_bmy17[3];
int smoothMPC_lbIdx17[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb17 = smoothMPC_l + 272;
smoothMPC_FLOAT* smoothMPC_slb17 = smoothMPC_s + 272;
smoothMPC_FLOAT* smoothMPC_llbbyslb17 = smoothMPC_lbys + 272;
smoothMPC_FLOAT smoothMPC_rilb17[11];
smoothMPC_FLOAT* smoothMPC_dllbaff17 = smoothMPC_dl_aff + 272;
smoothMPC_FLOAT* smoothMPC_dslbaff17 = smoothMPC_ds_aff + 272;
smoothMPC_FLOAT* smoothMPC_dllbcc17 = smoothMPC_dl_cc + 272;
smoothMPC_FLOAT* smoothMPC_dslbcc17 = smoothMPC_ds_cc + 272;
smoothMPC_FLOAT* smoothMPC_ccrhsl17 = smoothMPC_ccrhs + 272;
int smoothMPC_ubIdx17[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub17 = smoothMPC_l + 283;
smoothMPC_FLOAT* smoothMPC_sub17 = smoothMPC_s + 283;
smoothMPC_FLOAT* smoothMPC_lubbysub17 = smoothMPC_lbys + 283;
smoothMPC_FLOAT smoothMPC_riub17[5];
smoothMPC_FLOAT* smoothMPC_dlubaff17 = smoothMPC_dl_aff + 283;
smoothMPC_FLOAT* smoothMPC_dsubaff17 = smoothMPC_ds_aff + 283;
smoothMPC_FLOAT* smoothMPC_dlubcc17 = smoothMPC_dl_cc + 283;
smoothMPC_FLOAT* smoothMPC_dsubcc17 = smoothMPC_ds_cc + 283;
smoothMPC_FLOAT* smoothMPC_ccrhsub17 = smoothMPC_ccrhs + 283;
smoothMPC_FLOAT smoothMPC_Phi17[11];
smoothMPC_FLOAT smoothMPC_W17[11];
smoothMPC_FLOAT smoothMPC_Ysd17[9];
smoothMPC_FLOAT smoothMPC_Lsd17[9];
smoothMPC_FLOAT* smoothMPC_z18 = smoothMPC_z + 198;
smoothMPC_FLOAT* smoothMPC_dzaff18 = smoothMPC_dz_aff + 198;
smoothMPC_FLOAT* smoothMPC_dzcc18 = smoothMPC_dz_cc + 198;
smoothMPC_FLOAT* smoothMPC_rd18 = smoothMPC_rd + 198;
smoothMPC_FLOAT smoothMPC_Lbyrd18[11];
smoothMPC_FLOAT* smoothMPC_grad_cost18 = smoothMPC_grad_cost + 198;
smoothMPC_FLOAT* smoothMPC_grad_eq18 = smoothMPC_grad_eq + 198;
smoothMPC_FLOAT* smoothMPC_grad_ineq18 = smoothMPC_grad_ineq + 198;
smoothMPC_FLOAT smoothMPC_ctv18[11];
smoothMPC_FLOAT* smoothMPC_v18 = smoothMPC_v + 54;
smoothMPC_FLOAT smoothMPC_re18[3];
smoothMPC_FLOAT smoothMPC_beta18[3];
smoothMPC_FLOAT smoothMPC_betacc18[3];
smoothMPC_FLOAT* smoothMPC_dvaff18 = smoothMPC_dv_aff + 54;
smoothMPC_FLOAT* smoothMPC_dvcc18 = smoothMPC_dv_cc + 54;
smoothMPC_FLOAT smoothMPC_V18[33];
smoothMPC_FLOAT smoothMPC_Yd18[6];
smoothMPC_FLOAT smoothMPC_Ld18[6];
smoothMPC_FLOAT smoothMPC_yy18[3];
smoothMPC_FLOAT smoothMPC_bmy18[3];
int smoothMPC_lbIdx18[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb18 = smoothMPC_l + 288;
smoothMPC_FLOAT* smoothMPC_slb18 = smoothMPC_s + 288;
smoothMPC_FLOAT* smoothMPC_llbbyslb18 = smoothMPC_lbys + 288;
smoothMPC_FLOAT smoothMPC_rilb18[11];
smoothMPC_FLOAT* smoothMPC_dllbaff18 = smoothMPC_dl_aff + 288;
smoothMPC_FLOAT* smoothMPC_dslbaff18 = smoothMPC_ds_aff + 288;
smoothMPC_FLOAT* smoothMPC_dllbcc18 = smoothMPC_dl_cc + 288;
smoothMPC_FLOAT* smoothMPC_dslbcc18 = smoothMPC_ds_cc + 288;
smoothMPC_FLOAT* smoothMPC_ccrhsl18 = smoothMPC_ccrhs + 288;
int smoothMPC_ubIdx18[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub18 = smoothMPC_l + 299;
smoothMPC_FLOAT* smoothMPC_sub18 = smoothMPC_s + 299;
smoothMPC_FLOAT* smoothMPC_lubbysub18 = smoothMPC_lbys + 299;
smoothMPC_FLOAT smoothMPC_riub18[5];
smoothMPC_FLOAT* smoothMPC_dlubaff18 = smoothMPC_dl_aff + 299;
smoothMPC_FLOAT* smoothMPC_dsubaff18 = smoothMPC_ds_aff + 299;
smoothMPC_FLOAT* smoothMPC_dlubcc18 = smoothMPC_dl_cc + 299;
smoothMPC_FLOAT* smoothMPC_dsubcc18 = smoothMPC_ds_cc + 299;
smoothMPC_FLOAT* smoothMPC_ccrhsub18 = smoothMPC_ccrhs + 299;
smoothMPC_FLOAT smoothMPC_Phi18[11];
smoothMPC_FLOAT smoothMPC_W18[11];
smoothMPC_FLOAT smoothMPC_Ysd18[9];
smoothMPC_FLOAT smoothMPC_Lsd18[9];
smoothMPC_FLOAT* smoothMPC_z19 = smoothMPC_z + 209;
smoothMPC_FLOAT* smoothMPC_dzaff19 = smoothMPC_dz_aff + 209;
smoothMPC_FLOAT* smoothMPC_dzcc19 = smoothMPC_dz_cc + 209;
smoothMPC_FLOAT* smoothMPC_rd19 = smoothMPC_rd + 209;
smoothMPC_FLOAT smoothMPC_Lbyrd19[11];
smoothMPC_FLOAT* smoothMPC_grad_cost19 = smoothMPC_grad_cost + 209;
smoothMPC_FLOAT* smoothMPC_grad_eq19 = smoothMPC_grad_eq + 209;
smoothMPC_FLOAT* smoothMPC_grad_ineq19 = smoothMPC_grad_ineq + 209;
smoothMPC_FLOAT smoothMPC_ctv19[11];
smoothMPC_FLOAT* smoothMPC_v19 = smoothMPC_v + 57;
smoothMPC_FLOAT smoothMPC_re19[3];
smoothMPC_FLOAT smoothMPC_beta19[3];
smoothMPC_FLOAT smoothMPC_betacc19[3];
smoothMPC_FLOAT* smoothMPC_dvaff19 = smoothMPC_dv_aff + 57;
smoothMPC_FLOAT* smoothMPC_dvcc19 = smoothMPC_dv_cc + 57;
smoothMPC_FLOAT smoothMPC_V19[33];
smoothMPC_FLOAT smoothMPC_Yd19[6];
smoothMPC_FLOAT smoothMPC_Ld19[6];
smoothMPC_FLOAT smoothMPC_yy19[3];
smoothMPC_FLOAT smoothMPC_bmy19[3];
int smoothMPC_lbIdx19[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb19 = smoothMPC_l + 304;
smoothMPC_FLOAT* smoothMPC_slb19 = smoothMPC_s + 304;
smoothMPC_FLOAT* smoothMPC_llbbyslb19 = smoothMPC_lbys + 304;
smoothMPC_FLOAT smoothMPC_rilb19[11];
smoothMPC_FLOAT* smoothMPC_dllbaff19 = smoothMPC_dl_aff + 304;
smoothMPC_FLOAT* smoothMPC_dslbaff19 = smoothMPC_ds_aff + 304;
smoothMPC_FLOAT* smoothMPC_dllbcc19 = smoothMPC_dl_cc + 304;
smoothMPC_FLOAT* smoothMPC_dslbcc19 = smoothMPC_ds_cc + 304;
smoothMPC_FLOAT* smoothMPC_ccrhsl19 = smoothMPC_ccrhs + 304;
int smoothMPC_ubIdx19[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub19 = smoothMPC_l + 315;
smoothMPC_FLOAT* smoothMPC_sub19 = smoothMPC_s + 315;
smoothMPC_FLOAT* smoothMPC_lubbysub19 = smoothMPC_lbys + 315;
smoothMPC_FLOAT smoothMPC_riub19[5];
smoothMPC_FLOAT* smoothMPC_dlubaff19 = smoothMPC_dl_aff + 315;
smoothMPC_FLOAT* smoothMPC_dsubaff19 = smoothMPC_ds_aff + 315;
smoothMPC_FLOAT* smoothMPC_dlubcc19 = smoothMPC_dl_cc + 315;
smoothMPC_FLOAT* smoothMPC_dsubcc19 = smoothMPC_ds_cc + 315;
smoothMPC_FLOAT* smoothMPC_ccrhsub19 = smoothMPC_ccrhs + 315;
smoothMPC_FLOAT smoothMPC_Phi19[11];
smoothMPC_FLOAT smoothMPC_W19[11];
smoothMPC_FLOAT smoothMPC_Ysd19[9];
smoothMPC_FLOAT smoothMPC_Lsd19[9];
smoothMPC_FLOAT* smoothMPC_z20 = smoothMPC_z + 220;
smoothMPC_FLOAT* smoothMPC_dzaff20 = smoothMPC_dz_aff + 220;
smoothMPC_FLOAT* smoothMPC_dzcc20 = smoothMPC_dz_cc + 220;
smoothMPC_FLOAT* smoothMPC_rd20 = smoothMPC_rd + 220;
smoothMPC_FLOAT smoothMPC_Lbyrd20[11];
smoothMPC_FLOAT* smoothMPC_grad_cost20 = smoothMPC_grad_cost + 220;
smoothMPC_FLOAT* smoothMPC_grad_eq20 = smoothMPC_grad_eq + 220;
smoothMPC_FLOAT* smoothMPC_grad_ineq20 = smoothMPC_grad_ineq + 220;
smoothMPC_FLOAT smoothMPC_ctv20[11];
smoothMPC_FLOAT* smoothMPC_v20 = smoothMPC_v + 60;
smoothMPC_FLOAT smoothMPC_re20[3];
smoothMPC_FLOAT smoothMPC_beta20[3];
smoothMPC_FLOAT smoothMPC_betacc20[3];
smoothMPC_FLOAT* smoothMPC_dvaff20 = smoothMPC_dv_aff + 60;
smoothMPC_FLOAT* smoothMPC_dvcc20 = smoothMPC_dv_cc + 60;
smoothMPC_FLOAT smoothMPC_V20[33];
smoothMPC_FLOAT smoothMPC_Yd20[6];
smoothMPC_FLOAT smoothMPC_Ld20[6];
smoothMPC_FLOAT smoothMPC_yy20[3];
smoothMPC_FLOAT smoothMPC_bmy20[3];
int smoothMPC_lbIdx20[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb20 = smoothMPC_l + 320;
smoothMPC_FLOAT* smoothMPC_slb20 = smoothMPC_s + 320;
smoothMPC_FLOAT* smoothMPC_llbbyslb20 = smoothMPC_lbys + 320;
smoothMPC_FLOAT smoothMPC_rilb20[11];
smoothMPC_FLOAT* smoothMPC_dllbaff20 = smoothMPC_dl_aff + 320;
smoothMPC_FLOAT* smoothMPC_dslbaff20 = smoothMPC_ds_aff + 320;
smoothMPC_FLOAT* smoothMPC_dllbcc20 = smoothMPC_dl_cc + 320;
smoothMPC_FLOAT* smoothMPC_dslbcc20 = smoothMPC_ds_cc + 320;
smoothMPC_FLOAT* smoothMPC_ccrhsl20 = smoothMPC_ccrhs + 320;
int smoothMPC_ubIdx20[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub20 = smoothMPC_l + 331;
smoothMPC_FLOAT* smoothMPC_sub20 = smoothMPC_s + 331;
smoothMPC_FLOAT* smoothMPC_lubbysub20 = smoothMPC_lbys + 331;
smoothMPC_FLOAT smoothMPC_riub20[5];
smoothMPC_FLOAT* smoothMPC_dlubaff20 = smoothMPC_dl_aff + 331;
smoothMPC_FLOAT* smoothMPC_dsubaff20 = smoothMPC_ds_aff + 331;
smoothMPC_FLOAT* smoothMPC_dlubcc20 = smoothMPC_dl_cc + 331;
smoothMPC_FLOAT* smoothMPC_dsubcc20 = smoothMPC_ds_cc + 331;
smoothMPC_FLOAT* smoothMPC_ccrhsub20 = smoothMPC_ccrhs + 331;
smoothMPC_FLOAT smoothMPC_Phi20[11];
smoothMPC_FLOAT smoothMPC_W20[11];
smoothMPC_FLOAT smoothMPC_Ysd20[9];
smoothMPC_FLOAT smoothMPC_Lsd20[9];
smoothMPC_FLOAT* smoothMPC_z21 = smoothMPC_z + 231;
smoothMPC_FLOAT* smoothMPC_dzaff21 = smoothMPC_dz_aff + 231;
smoothMPC_FLOAT* smoothMPC_dzcc21 = smoothMPC_dz_cc + 231;
smoothMPC_FLOAT* smoothMPC_rd21 = smoothMPC_rd + 231;
smoothMPC_FLOAT smoothMPC_Lbyrd21[11];
smoothMPC_FLOAT* smoothMPC_grad_cost21 = smoothMPC_grad_cost + 231;
smoothMPC_FLOAT* smoothMPC_grad_eq21 = smoothMPC_grad_eq + 231;
smoothMPC_FLOAT* smoothMPC_grad_ineq21 = smoothMPC_grad_ineq + 231;
smoothMPC_FLOAT smoothMPC_ctv21[11];
smoothMPC_FLOAT* smoothMPC_v21 = smoothMPC_v + 63;
smoothMPC_FLOAT smoothMPC_re21[3];
smoothMPC_FLOAT smoothMPC_beta21[3];
smoothMPC_FLOAT smoothMPC_betacc21[3];
smoothMPC_FLOAT* smoothMPC_dvaff21 = smoothMPC_dv_aff + 63;
smoothMPC_FLOAT* smoothMPC_dvcc21 = smoothMPC_dv_cc + 63;
smoothMPC_FLOAT smoothMPC_V21[33];
smoothMPC_FLOAT smoothMPC_Yd21[6];
smoothMPC_FLOAT smoothMPC_Ld21[6];
smoothMPC_FLOAT smoothMPC_yy21[3];
smoothMPC_FLOAT smoothMPC_bmy21[3];
int smoothMPC_lbIdx21[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb21 = smoothMPC_l + 336;
smoothMPC_FLOAT* smoothMPC_slb21 = smoothMPC_s + 336;
smoothMPC_FLOAT* smoothMPC_llbbyslb21 = smoothMPC_lbys + 336;
smoothMPC_FLOAT smoothMPC_rilb21[11];
smoothMPC_FLOAT* smoothMPC_dllbaff21 = smoothMPC_dl_aff + 336;
smoothMPC_FLOAT* smoothMPC_dslbaff21 = smoothMPC_ds_aff + 336;
smoothMPC_FLOAT* smoothMPC_dllbcc21 = smoothMPC_dl_cc + 336;
smoothMPC_FLOAT* smoothMPC_dslbcc21 = smoothMPC_ds_cc + 336;
smoothMPC_FLOAT* smoothMPC_ccrhsl21 = smoothMPC_ccrhs + 336;
int smoothMPC_ubIdx21[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub21 = smoothMPC_l + 347;
smoothMPC_FLOAT* smoothMPC_sub21 = smoothMPC_s + 347;
smoothMPC_FLOAT* smoothMPC_lubbysub21 = smoothMPC_lbys + 347;
smoothMPC_FLOAT smoothMPC_riub21[5];
smoothMPC_FLOAT* smoothMPC_dlubaff21 = smoothMPC_dl_aff + 347;
smoothMPC_FLOAT* smoothMPC_dsubaff21 = smoothMPC_ds_aff + 347;
smoothMPC_FLOAT* smoothMPC_dlubcc21 = smoothMPC_dl_cc + 347;
smoothMPC_FLOAT* smoothMPC_dsubcc21 = smoothMPC_ds_cc + 347;
smoothMPC_FLOAT* smoothMPC_ccrhsub21 = smoothMPC_ccrhs + 347;
smoothMPC_FLOAT smoothMPC_Phi21[11];
smoothMPC_FLOAT smoothMPC_W21[11];
smoothMPC_FLOAT smoothMPC_Ysd21[9];
smoothMPC_FLOAT smoothMPC_Lsd21[9];
smoothMPC_FLOAT* smoothMPC_z22 = smoothMPC_z + 242;
smoothMPC_FLOAT* smoothMPC_dzaff22 = smoothMPC_dz_aff + 242;
smoothMPC_FLOAT* smoothMPC_dzcc22 = smoothMPC_dz_cc + 242;
smoothMPC_FLOAT* smoothMPC_rd22 = smoothMPC_rd + 242;
smoothMPC_FLOAT smoothMPC_Lbyrd22[11];
smoothMPC_FLOAT* smoothMPC_grad_cost22 = smoothMPC_grad_cost + 242;
smoothMPC_FLOAT* smoothMPC_grad_eq22 = smoothMPC_grad_eq + 242;
smoothMPC_FLOAT* smoothMPC_grad_ineq22 = smoothMPC_grad_ineq + 242;
smoothMPC_FLOAT smoothMPC_ctv22[11];
smoothMPC_FLOAT* smoothMPC_v22 = smoothMPC_v + 66;
smoothMPC_FLOAT smoothMPC_re22[3];
smoothMPC_FLOAT smoothMPC_beta22[3];
smoothMPC_FLOAT smoothMPC_betacc22[3];
smoothMPC_FLOAT* smoothMPC_dvaff22 = smoothMPC_dv_aff + 66;
smoothMPC_FLOAT* smoothMPC_dvcc22 = smoothMPC_dv_cc + 66;
smoothMPC_FLOAT smoothMPC_V22[33];
smoothMPC_FLOAT smoothMPC_Yd22[6];
smoothMPC_FLOAT smoothMPC_Ld22[6];
smoothMPC_FLOAT smoothMPC_yy22[3];
smoothMPC_FLOAT smoothMPC_bmy22[3];
int smoothMPC_lbIdx22[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb22 = smoothMPC_l + 352;
smoothMPC_FLOAT* smoothMPC_slb22 = smoothMPC_s + 352;
smoothMPC_FLOAT* smoothMPC_llbbyslb22 = smoothMPC_lbys + 352;
smoothMPC_FLOAT smoothMPC_rilb22[11];
smoothMPC_FLOAT* smoothMPC_dllbaff22 = smoothMPC_dl_aff + 352;
smoothMPC_FLOAT* smoothMPC_dslbaff22 = smoothMPC_ds_aff + 352;
smoothMPC_FLOAT* smoothMPC_dllbcc22 = smoothMPC_dl_cc + 352;
smoothMPC_FLOAT* smoothMPC_dslbcc22 = smoothMPC_ds_cc + 352;
smoothMPC_FLOAT* smoothMPC_ccrhsl22 = smoothMPC_ccrhs + 352;
int smoothMPC_ubIdx22[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub22 = smoothMPC_l + 363;
smoothMPC_FLOAT* smoothMPC_sub22 = smoothMPC_s + 363;
smoothMPC_FLOAT* smoothMPC_lubbysub22 = smoothMPC_lbys + 363;
smoothMPC_FLOAT smoothMPC_riub22[5];
smoothMPC_FLOAT* smoothMPC_dlubaff22 = smoothMPC_dl_aff + 363;
smoothMPC_FLOAT* smoothMPC_dsubaff22 = smoothMPC_ds_aff + 363;
smoothMPC_FLOAT* smoothMPC_dlubcc22 = smoothMPC_dl_cc + 363;
smoothMPC_FLOAT* smoothMPC_dsubcc22 = smoothMPC_ds_cc + 363;
smoothMPC_FLOAT* smoothMPC_ccrhsub22 = smoothMPC_ccrhs + 363;
smoothMPC_FLOAT smoothMPC_Phi22[11];
smoothMPC_FLOAT smoothMPC_W22[11];
smoothMPC_FLOAT smoothMPC_Ysd22[9];
smoothMPC_FLOAT smoothMPC_Lsd22[9];
smoothMPC_FLOAT* smoothMPC_z23 = smoothMPC_z + 253;
smoothMPC_FLOAT* smoothMPC_dzaff23 = smoothMPC_dz_aff + 253;
smoothMPC_FLOAT* smoothMPC_dzcc23 = smoothMPC_dz_cc + 253;
smoothMPC_FLOAT* smoothMPC_rd23 = smoothMPC_rd + 253;
smoothMPC_FLOAT smoothMPC_Lbyrd23[11];
smoothMPC_FLOAT* smoothMPC_grad_cost23 = smoothMPC_grad_cost + 253;
smoothMPC_FLOAT* smoothMPC_grad_eq23 = smoothMPC_grad_eq + 253;
smoothMPC_FLOAT* smoothMPC_grad_ineq23 = smoothMPC_grad_ineq + 253;
smoothMPC_FLOAT smoothMPC_ctv23[11];
smoothMPC_FLOAT* smoothMPC_v23 = smoothMPC_v + 69;
smoothMPC_FLOAT smoothMPC_re23[3];
smoothMPC_FLOAT smoothMPC_beta23[3];
smoothMPC_FLOAT smoothMPC_betacc23[3];
smoothMPC_FLOAT* smoothMPC_dvaff23 = smoothMPC_dv_aff + 69;
smoothMPC_FLOAT* smoothMPC_dvcc23 = smoothMPC_dv_cc + 69;
smoothMPC_FLOAT smoothMPC_V23[33];
smoothMPC_FLOAT smoothMPC_Yd23[6];
smoothMPC_FLOAT smoothMPC_Ld23[6];
smoothMPC_FLOAT smoothMPC_yy23[3];
smoothMPC_FLOAT smoothMPC_bmy23[3];
int smoothMPC_lbIdx23[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb23 = smoothMPC_l + 368;
smoothMPC_FLOAT* smoothMPC_slb23 = smoothMPC_s + 368;
smoothMPC_FLOAT* smoothMPC_llbbyslb23 = smoothMPC_lbys + 368;
smoothMPC_FLOAT smoothMPC_rilb23[11];
smoothMPC_FLOAT* smoothMPC_dllbaff23 = smoothMPC_dl_aff + 368;
smoothMPC_FLOAT* smoothMPC_dslbaff23 = smoothMPC_ds_aff + 368;
smoothMPC_FLOAT* smoothMPC_dllbcc23 = smoothMPC_dl_cc + 368;
smoothMPC_FLOAT* smoothMPC_dslbcc23 = smoothMPC_ds_cc + 368;
smoothMPC_FLOAT* smoothMPC_ccrhsl23 = smoothMPC_ccrhs + 368;
int smoothMPC_ubIdx23[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub23 = smoothMPC_l + 379;
smoothMPC_FLOAT* smoothMPC_sub23 = smoothMPC_s + 379;
smoothMPC_FLOAT* smoothMPC_lubbysub23 = smoothMPC_lbys + 379;
smoothMPC_FLOAT smoothMPC_riub23[5];
smoothMPC_FLOAT* smoothMPC_dlubaff23 = smoothMPC_dl_aff + 379;
smoothMPC_FLOAT* smoothMPC_dsubaff23 = smoothMPC_ds_aff + 379;
smoothMPC_FLOAT* smoothMPC_dlubcc23 = smoothMPC_dl_cc + 379;
smoothMPC_FLOAT* smoothMPC_dsubcc23 = smoothMPC_ds_cc + 379;
smoothMPC_FLOAT* smoothMPC_ccrhsub23 = smoothMPC_ccrhs + 379;
smoothMPC_FLOAT smoothMPC_Phi23[11];
smoothMPC_FLOAT smoothMPC_W23[11];
smoothMPC_FLOAT smoothMPC_Ysd23[9];
smoothMPC_FLOAT smoothMPC_Lsd23[9];
smoothMPC_FLOAT* smoothMPC_z24 = smoothMPC_z + 264;
smoothMPC_FLOAT* smoothMPC_dzaff24 = smoothMPC_dz_aff + 264;
smoothMPC_FLOAT* smoothMPC_dzcc24 = smoothMPC_dz_cc + 264;
smoothMPC_FLOAT* smoothMPC_rd24 = smoothMPC_rd + 264;
smoothMPC_FLOAT smoothMPC_Lbyrd24[11];
smoothMPC_FLOAT* smoothMPC_grad_cost24 = smoothMPC_grad_cost + 264;
smoothMPC_FLOAT* smoothMPC_grad_eq24 = smoothMPC_grad_eq + 264;
smoothMPC_FLOAT* smoothMPC_grad_ineq24 = smoothMPC_grad_ineq + 264;
smoothMPC_FLOAT smoothMPC_ctv24[11];
smoothMPC_FLOAT* smoothMPC_v24 = smoothMPC_v + 72;
smoothMPC_FLOAT smoothMPC_re24[3];
smoothMPC_FLOAT smoothMPC_beta24[3];
smoothMPC_FLOAT smoothMPC_betacc24[3];
smoothMPC_FLOAT* smoothMPC_dvaff24 = smoothMPC_dv_aff + 72;
smoothMPC_FLOAT* smoothMPC_dvcc24 = smoothMPC_dv_cc + 72;
smoothMPC_FLOAT smoothMPC_V24[33];
smoothMPC_FLOAT smoothMPC_Yd24[6];
smoothMPC_FLOAT smoothMPC_Ld24[6];
smoothMPC_FLOAT smoothMPC_yy24[3];
smoothMPC_FLOAT smoothMPC_bmy24[3];
int smoothMPC_lbIdx24[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb24 = smoothMPC_l + 384;
smoothMPC_FLOAT* smoothMPC_slb24 = smoothMPC_s + 384;
smoothMPC_FLOAT* smoothMPC_llbbyslb24 = smoothMPC_lbys + 384;
smoothMPC_FLOAT smoothMPC_rilb24[11];
smoothMPC_FLOAT* smoothMPC_dllbaff24 = smoothMPC_dl_aff + 384;
smoothMPC_FLOAT* smoothMPC_dslbaff24 = smoothMPC_ds_aff + 384;
smoothMPC_FLOAT* smoothMPC_dllbcc24 = smoothMPC_dl_cc + 384;
smoothMPC_FLOAT* smoothMPC_dslbcc24 = smoothMPC_ds_cc + 384;
smoothMPC_FLOAT* smoothMPC_ccrhsl24 = smoothMPC_ccrhs + 384;
int smoothMPC_ubIdx24[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub24 = smoothMPC_l + 395;
smoothMPC_FLOAT* smoothMPC_sub24 = smoothMPC_s + 395;
smoothMPC_FLOAT* smoothMPC_lubbysub24 = smoothMPC_lbys + 395;
smoothMPC_FLOAT smoothMPC_riub24[5];
smoothMPC_FLOAT* smoothMPC_dlubaff24 = smoothMPC_dl_aff + 395;
smoothMPC_FLOAT* smoothMPC_dsubaff24 = smoothMPC_ds_aff + 395;
smoothMPC_FLOAT* smoothMPC_dlubcc24 = smoothMPC_dl_cc + 395;
smoothMPC_FLOAT* smoothMPC_dsubcc24 = smoothMPC_ds_cc + 395;
smoothMPC_FLOAT* smoothMPC_ccrhsub24 = smoothMPC_ccrhs + 395;
smoothMPC_FLOAT smoothMPC_Phi24[11];
smoothMPC_FLOAT smoothMPC_W24[11];
smoothMPC_FLOAT smoothMPC_Ysd24[9];
smoothMPC_FLOAT smoothMPC_Lsd24[9];
smoothMPC_FLOAT* smoothMPC_z25 = smoothMPC_z + 275;
smoothMPC_FLOAT* smoothMPC_dzaff25 = smoothMPC_dz_aff + 275;
smoothMPC_FLOAT* smoothMPC_dzcc25 = smoothMPC_dz_cc + 275;
smoothMPC_FLOAT* smoothMPC_rd25 = smoothMPC_rd + 275;
smoothMPC_FLOAT smoothMPC_Lbyrd25[11];
smoothMPC_FLOAT* smoothMPC_grad_cost25 = smoothMPC_grad_cost + 275;
smoothMPC_FLOAT* smoothMPC_grad_eq25 = smoothMPC_grad_eq + 275;
smoothMPC_FLOAT* smoothMPC_grad_ineq25 = smoothMPC_grad_ineq + 275;
smoothMPC_FLOAT smoothMPC_ctv25[11];
smoothMPC_FLOAT* smoothMPC_v25 = smoothMPC_v + 75;
smoothMPC_FLOAT smoothMPC_re25[3];
smoothMPC_FLOAT smoothMPC_beta25[3];
smoothMPC_FLOAT smoothMPC_betacc25[3];
smoothMPC_FLOAT* smoothMPC_dvaff25 = smoothMPC_dv_aff + 75;
smoothMPC_FLOAT* smoothMPC_dvcc25 = smoothMPC_dv_cc + 75;
smoothMPC_FLOAT smoothMPC_V25[33];
smoothMPC_FLOAT smoothMPC_Yd25[6];
smoothMPC_FLOAT smoothMPC_Ld25[6];
smoothMPC_FLOAT smoothMPC_yy25[3];
smoothMPC_FLOAT smoothMPC_bmy25[3];
int smoothMPC_lbIdx25[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb25 = smoothMPC_l + 400;
smoothMPC_FLOAT* smoothMPC_slb25 = smoothMPC_s + 400;
smoothMPC_FLOAT* smoothMPC_llbbyslb25 = smoothMPC_lbys + 400;
smoothMPC_FLOAT smoothMPC_rilb25[11];
smoothMPC_FLOAT* smoothMPC_dllbaff25 = smoothMPC_dl_aff + 400;
smoothMPC_FLOAT* smoothMPC_dslbaff25 = smoothMPC_ds_aff + 400;
smoothMPC_FLOAT* smoothMPC_dllbcc25 = smoothMPC_dl_cc + 400;
smoothMPC_FLOAT* smoothMPC_dslbcc25 = smoothMPC_ds_cc + 400;
smoothMPC_FLOAT* smoothMPC_ccrhsl25 = smoothMPC_ccrhs + 400;
int smoothMPC_ubIdx25[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub25 = smoothMPC_l + 411;
smoothMPC_FLOAT* smoothMPC_sub25 = smoothMPC_s + 411;
smoothMPC_FLOAT* smoothMPC_lubbysub25 = smoothMPC_lbys + 411;
smoothMPC_FLOAT smoothMPC_riub25[5];
smoothMPC_FLOAT* smoothMPC_dlubaff25 = smoothMPC_dl_aff + 411;
smoothMPC_FLOAT* smoothMPC_dsubaff25 = smoothMPC_ds_aff + 411;
smoothMPC_FLOAT* smoothMPC_dlubcc25 = smoothMPC_dl_cc + 411;
smoothMPC_FLOAT* smoothMPC_dsubcc25 = smoothMPC_ds_cc + 411;
smoothMPC_FLOAT* smoothMPC_ccrhsub25 = smoothMPC_ccrhs + 411;
smoothMPC_FLOAT smoothMPC_Phi25[11];
smoothMPC_FLOAT smoothMPC_W25[11];
smoothMPC_FLOAT smoothMPC_Ysd25[9];
smoothMPC_FLOAT smoothMPC_Lsd25[9];
smoothMPC_FLOAT* smoothMPC_z26 = smoothMPC_z + 286;
smoothMPC_FLOAT* smoothMPC_dzaff26 = smoothMPC_dz_aff + 286;
smoothMPC_FLOAT* smoothMPC_dzcc26 = smoothMPC_dz_cc + 286;
smoothMPC_FLOAT* smoothMPC_rd26 = smoothMPC_rd + 286;
smoothMPC_FLOAT smoothMPC_Lbyrd26[11];
smoothMPC_FLOAT* smoothMPC_grad_cost26 = smoothMPC_grad_cost + 286;
smoothMPC_FLOAT* smoothMPC_grad_eq26 = smoothMPC_grad_eq + 286;
smoothMPC_FLOAT* smoothMPC_grad_ineq26 = smoothMPC_grad_ineq + 286;
smoothMPC_FLOAT smoothMPC_ctv26[11];
smoothMPC_FLOAT* smoothMPC_v26 = smoothMPC_v + 78;
smoothMPC_FLOAT smoothMPC_re26[3];
smoothMPC_FLOAT smoothMPC_beta26[3];
smoothMPC_FLOAT smoothMPC_betacc26[3];
smoothMPC_FLOAT* smoothMPC_dvaff26 = smoothMPC_dv_aff + 78;
smoothMPC_FLOAT* smoothMPC_dvcc26 = smoothMPC_dv_cc + 78;
smoothMPC_FLOAT smoothMPC_V26[33];
smoothMPC_FLOAT smoothMPC_Yd26[6];
smoothMPC_FLOAT smoothMPC_Ld26[6];
smoothMPC_FLOAT smoothMPC_yy26[3];
smoothMPC_FLOAT smoothMPC_bmy26[3];
int smoothMPC_lbIdx26[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb26 = smoothMPC_l + 416;
smoothMPC_FLOAT* smoothMPC_slb26 = smoothMPC_s + 416;
smoothMPC_FLOAT* smoothMPC_llbbyslb26 = smoothMPC_lbys + 416;
smoothMPC_FLOAT smoothMPC_rilb26[11];
smoothMPC_FLOAT* smoothMPC_dllbaff26 = smoothMPC_dl_aff + 416;
smoothMPC_FLOAT* smoothMPC_dslbaff26 = smoothMPC_ds_aff + 416;
smoothMPC_FLOAT* smoothMPC_dllbcc26 = smoothMPC_dl_cc + 416;
smoothMPC_FLOAT* smoothMPC_dslbcc26 = smoothMPC_ds_cc + 416;
smoothMPC_FLOAT* smoothMPC_ccrhsl26 = smoothMPC_ccrhs + 416;
int smoothMPC_ubIdx26[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub26 = smoothMPC_l + 427;
smoothMPC_FLOAT* smoothMPC_sub26 = smoothMPC_s + 427;
smoothMPC_FLOAT* smoothMPC_lubbysub26 = smoothMPC_lbys + 427;
smoothMPC_FLOAT smoothMPC_riub26[5];
smoothMPC_FLOAT* smoothMPC_dlubaff26 = smoothMPC_dl_aff + 427;
smoothMPC_FLOAT* smoothMPC_dsubaff26 = smoothMPC_ds_aff + 427;
smoothMPC_FLOAT* smoothMPC_dlubcc26 = smoothMPC_dl_cc + 427;
smoothMPC_FLOAT* smoothMPC_dsubcc26 = smoothMPC_ds_cc + 427;
smoothMPC_FLOAT* smoothMPC_ccrhsub26 = smoothMPC_ccrhs + 427;
smoothMPC_FLOAT smoothMPC_Phi26[11];
smoothMPC_FLOAT smoothMPC_W26[11];
smoothMPC_FLOAT smoothMPC_Ysd26[9];
smoothMPC_FLOAT smoothMPC_Lsd26[9];
smoothMPC_FLOAT* smoothMPC_z27 = smoothMPC_z + 297;
smoothMPC_FLOAT* smoothMPC_dzaff27 = smoothMPC_dz_aff + 297;
smoothMPC_FLOAT* smoothMPC_dzcc27 = smoothMPC_dz_cc + 297;
smoothMPC_FLOAT* smoothMPC_rd27 = smoothMPC_rd + 297;
smoothMPC_FLOAT smoothMPC_Lbyrd27[11];
smoothMPC_FLOAT* smoothMPC_grad_cost27 = smoothMPC_grad_cost + 297;
smoothMPC_FLOAT* smoothMPC_grad_eq27 = smoothMPC_grad_eq + 297;
smoothMPC_FLOAT* smoothMPC_grad_ineq27 = smoothMPC_grad_ineq + 297;
smoothMPC_FLOAT smoothMPC_ctv27[11];
smoothMPC_FLOAT* smoothMPC_v27 = smoothMPC_v + 81;
smoothMPC_FLOAT smoothMPC_re27[3];
smoothMPC_FLOAT smoothMPC_beta27[3];
smoothMPC_FLOAT smoothMPC_betacc27[3];
smoothMPC_FLOAT* smoothMPC_dvaff27 = smoothMPC_dv_aff + 81;
smoothMPC_FLOAT* smoothMPC_dvcc27 = smoothMPC_dv_cc + 81;
smoothMPC_FLOAT smoothMPC_V27[33];
smoothMPC_FLOAT smoothMPC_Yd27[6];
smoothMPC_FLOAT smoothMPC_Ld27[6];
smoothMPC_FLOAT smoothMPC_yy27[3];
smoothMPC_FLOAT smoothMPC_bmy27[3];
int smoothMPC_lbIdx27[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb27 = smoothMPC_l + 432;
smoothMPC_FLOAT* smoothMPC_slb27 = smoothMPC_s + 432;
smoothMPC_FLOAT* smoothMPC_llbbyslb27 = smoothMPC_lbys + 432;
smoothMPC_FLOAT smoothMPC_rilb27[11];
smoothMPC_FLOAT* smoothMPC_dllbaff27 = smoothMPC_dl_aff + 432;
smoothMPC_FLOAT* smoothMPC_dslbaff27 = smoothMPC_ds_aff + 432;
smoothMPC_FLOAT* smoothMPC_dllbcc27 = smoothMPC_dl_cc + 432;
smoothMPC_FLOAT* smoothMPC_dslbcc27 = smoothMPC_ds_cc + 432;
smoothMPC_FLOAT* smoothMPC_ccrhsl27 = smoothMPC_ccrhs + 432;
int smoothMPC_ubIdx27[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub27 = smoothMPC_l + 443;
smoothMPC_FLOAT* smoothMPC_sub27 = smoothMPC_s + 443;
smoothMPC_FLOAT* smoothMPC_lubbysub27 = smoothMPC_lbys + 443;
smoothMPC_FLOAT smoothMPC_riub27[5];
smoothMPC_FLOAT* smoothMPC_dlubaff27 = smoothMPC_dl_aff + 443;
smoothMPC_FLOAT* smoothMPC_dsubaff27 = smoothMPC_ds_aff + 443;
smoothMPC_FLOAT* smoothMPC_dlubcc27 = smoothMPC_dl_cc + 443;
smoothMPC_FLOAT* smoothMPC_dsubcc27 = smoothMPC_ds_cc + 443;
smoothMPC_FLOAT* smoothMPC_ccrhsub27 = smoothMPC_ccrhs + 443;
smoothMPC_FLOAT smoothMPC_Phi27[11];
smoothMPC_FLOAT smoothMPC_W27[11];
smoothMPC_FLOAT smoothMPC_Ysd27[9];
smoothMPC_FLOAT smoothMPC_Lsd27[9];
smoothMPC_FLOAT* smoothMPC_z28 = smoothMPC_z + 308;
smoothMPC_FLOAT* smoothMPC_dzaff28 = smoothMPC_dz_aff + 308;
smoothMPC_FLOAT* smoothMPC_dzcc28 = smoothMPC_dz_cc + 308;
smoothMPC_FLOAT* smoothMPC_rd28 = smoothMPC_rd + 308;
smoothMPC_FLOAT smoothMPC_Lbyrd28[11];
smoothMPC_FLOAT* smoothMPC_grad_cost28 = smoothMPC_grad_cost + 308;
smoothMPC_FLOAT* smoothMPC_grad_eq28 = smoothMPC_grad_eq + 308;
smoothMPC_FLOAT* smoothMPC_grad_ineq28 = smoothMPC_grad_ineq + 308;
smoothMPC_FLOAT smoothMPC_ctv28[11];
smoothMPC_FLOAT* smoothMPC_v28 = smoothMPC_v + 84;
smoothMPC_FLOAT smoothMPC_re28[3];
smoothMPC_FLOAT smoothMPC_beta28[3];
smoothMPC_FLOAT smoothMPC_betacc28[3];
smoothMPC_FLOAT* smoothMPC_dvaff28 = smoothMPC_dv_aff + 84;
smoothMPC_FLOAT* smoothMPC_dvcc28 = smoothMPC_dv_cc + 84;
smoothMPC_FLOAT smoothMPC_V28[33];
smoothMPC_FLOAT smoothMPC_Yd28[6];
smoothMPC_FLOAT smoothMPC_Ld28[6];
smoothMPC_FLOAT smoothMPC_yy28[3];
smoothMPC_FLOAT smoothMPC_bmy28[3];
int smoothMPC_lbIdx28[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb28 = smoothMPC_l + 448;
smoothMPC_FLOAT* smoothMPC_slb28 = smoothMPC_s + 448;
smoothMPC_FLOAT* smoothMPC_llbbyslb28 = smoothMPC_lbys + 448;
smoothMPC_FLOAT smoothMPC_rilb28[11];
smoothMPC_FLOAT* smoothMPC_dllbaff28 = smoothMPC_dl_aff + 448;
smoothMPC_FLOAT* smoothMPC_dslbaff28 = smoothMPC_ds_aff + 448;
smoothMPC_FLOAT* smoothMPC_dllbcc28 = smoothMPC_dl_cc + 448;
smoothMPC_FLOAT* smoothMPC_dslbcc28 = smoothMPC_ds_cc + 448;
smoothMPC_FLOAT* smoothMPC_ccrhsl28 = smoothMPC_ccrhs + 448;
int smoothMPC_ubIdx28[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub28 = smoothMPC_l + 459;
smoothMPC_FLOAT* smoothMPC_sub28 = smoothMPC_s + 459;
smoothMPC_FLOAT* smoothMPC_lubbysub28 = smoothMPC_lbys + 459;
smoothMPC_FLOAT smoothMPC_riub28[5];
smoothMPC_FLOAT* smoothMPC_dlubaff28 = smoothMPC_dl_aff + 459;
smoothMPC_FLOAT* smoothMPC_dsubaff28 = smoothMPC_ds_aff + 459;
smoothMPC_FLOAT* smoothMPC_dlubcc28 = smoothMPC_dl_cc + 459;
smoothMPC_FLOAT* smoothMPC_dsubcc28 = smoothMPC_ds_cc + 459;
smoothMPC_FLOAT* smoothMPC_ccrhsub28 = smoothMPC_ccrhs + 459;
smoothMPC_FLOAT smoothMPC_Phi28[11];
smoothMPC_FLOAT smoothMPC_W28[11];
smoothMPC_FLOAT smoothMPC_Ysd28[9];
smoothMPC_FLOAT smoothMPC_Lsd28[9];
smoothMPC_FLOAT* smoothMPC_z29 = smoothMPC_z + 319;
smoothMPC_FLOAT* smoothMPC_dzaff29 = smoothMPC_dz_aff + 319;
smoothMPC_FLOAT* smoothMPC_dzcc29 = smoothMPC_dz_cc + 319;
smoothMPC_FLOAT* smoothMPC_rd29 = smoothMPC_rd + 319;
smoothMPC_FLOAT smoothMPC_Lbyrd29[11];
smoothMPC_FLOAT* smoothMPC_grad_cost29 = smoothMPC_grad_cost + 319;
smoothMPC_FLOAT* smoothMPC_grad_eq29 = smoothMPC_grad_eq + 319;
smoothMPC_FLOAT* smoothMPC_grad_ineq29 = smoothMPC_grad_ineq + 319;
smoothMPC_FLOAT smoothMPC_ctv29[11];
smoothMPC_FLOAT* smoothMPC_v29 = smoothMPC_v + 87;
smoothMPC_FLOAT smoothMPC_re29[3];
smoothMPC_FLOAT smoothMPC_beta29[3];
smoothMPC_FLOAT smoothMPC_betacc29[3];
smoothMPC_FLOAT* smoothMPC_dvaff29 = smoothMPC_dv_aff + 87;
smoothMPC_FLOAT* smoothMPC_dvcc29 = smoothMPC_dv_cc + 87;
smoothMPC_FLOAT smoothMPC_V29[33];
smoothMPC_FLOAT smoothMPC_Yd29[6];
smoothMPC_FLOAT smoothMPC_Ld29[6];
smoothMPC_FLOAT smoothMPC_yy29[3];
smoothMPC_FLOAT smoothMPC_bmy29[3];
int smoothMPC_lbIdx29[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb29 = smoothMPC_l + 464;
smoothMPC_FLOAT* smoothMPC_slb29 = smoothMPC_s + 464;
smoothMPC_FLOAT* smoothMPC_llbbyslb29 = smoothMPC_lbys + 464;
smoothMPC_FLOAT smoothMPC_rilb29[11];
smoothMPC_FLOAT* smoothMPC_dllbaff29 = smoothMPC_dl_aff + 464;
smoothMPC_FLOAT* smoothMPC_dslbaff29 = smoothMPC_ds_aff + 464;
smoothMPC_FLOAT* smoothMPC_dllbcc29 = smoothMPC_dl_cc + 464;
smoothMPC_FLOAT* smoothMPC_dslbcc29 = smoothMPC_ds_cc + 464;
smoothMPC_FLOAT* smoothMPC_ccrhsl29 = smoothMPC_ccrhs + 464;
int smoothMPC_ubIdx29[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub29 = smoothMPC_l + 475;
smoothMPC_FLOAT* smoothMPC_sub29 = smoothMPC_s + 475;
smoothMPC_FLOAT* smoothMPC_lubbysub29 = smoothMPC_lbys + 475;
smoothMPC_FLOAT smoothMPC_riub29[5];
smoothMPC_FLOAT* smoothMPC_dlubaff29 = smoothMPC_dl_aff + 475;
smoothMPC_FLOAT* smoothMPC_dsubaff29 = smoothMPC_ds_aff + 475;
smoothMPC_FLOAT* smoothMPC_dlubcc29 = smoothMPC_dl_cc + 475;
smoothMPC_FLOAT* smoothMPC_dsubcc29 = smoothMPC_ds_cc + 475;
smoothMPC_FLOAT* smoothMPC_ccrhsub29 = smoothMPC_ccrhs + 475;
smoothMPC_FLOAT smoothMPC_Phi29[11];
smoothMPC_FLOAT smoothMPC_W29[11];
smoothMPC_FLOAT smoothMPC_Ysd29[9];
smoothMPC_FLOAT smoothMPC_Lsd29[9];
smoothMPC_FLOAT* smoothMPC_z30 = smoothMPC_z + 330;
smoothMPC_FLOAT* smoothMPC_dzaff30 = smoothMPC_dz_aff + 330;
smoothMPC_FLOAT* smoothMPC_dzcc30 = smoothMPC_dz_cc + 330;
smoothMPC_FLOAT* smoothMPC_rd30 = smoothMPC_rd + 330;
smoothMPC_FLOAT smoothMPC_Lbyrd30[11];
smoothMPC_FLOAT* smoothMPC_grad_cost30 = smoothMPC_grad_cost + 330;
smoothMPC_FLOAT* smoothMPC_grad_eq30 = smoothMPC_grad_eq + 330;
smoothMPC_FLOAT* smoothMPC_grad_ineq30 = smoothMPC_grad_ineq + 330;
smoothMPC_FLOAT smoothMPC_ctv30[11];
smoothMPC_FLOAT* smoothMPC_v30 = smoothMPC_v + 90;
smoothMPC_FLOAT smoothMPC_re30[3];
smoothMPC_FLOAT smoothMPC_beta30[3];
smoothMPC_FLOAT smoothMPC_betacc30[3];
smoothMPC_FLOAT* smoothMPC_dvaff30 = smoothMPC_dv_aff + 90;
smoothMPC_FLOAT* smoothMPC_dvcc30 = smoothMPC_dv_cc + 90;
smoothMPC_FLOAT smoothMPC_V30[33];
smoothMPC_FLOAT smoothMPC_Yd30[6];
smoothMPC_FLOAT smoothMPC_Ld30[6];
smoothMPC_FLOAT smoothMPC_yy30[3];
smoothMPC_FLOAT smoothMPC_bmy30[3];
int smoothMPC_lbIdx30[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb30 = smoothMPC_l + 480;
smoothMPC_FLOAT* smoothMPC_slb30 = smoothMPC_s + 480;
smoothMPC_FLOAT* smoothMPC_llbbyslb30 = smoothMPC_lbys + 480;
smoothMPC_FLOAT smoothMPC_rilb30[11];
smoothMPC_FLOAT* smoothMPC_dllbaff30 = smoothMPC_dl_aff + 480;
smoothMPC_FLOAT* smoothMPC_dslbaff30 = smoothMPC_ds_aff + 480;
smoothMPC_FLOAT* smoothMPC_dllbcc30 = smoothMPC_dl_cc + 480;
smoothMPC_FLOAT* smoothMPC_dslbcc30 = smoothMPC_ds_cc + 480;
smoothMPC_FLOAT* smoothMPC_ccrhsl30 = smoothMPC_ccrhs + 480;
int smoothMPC_ubIdx30[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub30 = smoothMPC_l + 491;
smoothMPC_FLOAT* smoothMPC_sub30 = smoothMPC_s + 491;
smoothMPC_FLOAT* smoothMPC_lubbysub30 = smoothMPC_lbys + 491;
smoothMPC_FLOAT smoothMPC_riub30[5];
smoothMPC_FLOAT* smoothMPC_dlubaff30 = smoothMPC_dl_aff + 491;
smoothMPC_FLOAT* smoothMPC_dsubaff30 = smoothMPC_ds_aff + 491;
smoothMPC_FLOAT* smoothMPC_dlubcc30 = smoothMPC_dl_cc + 491;
smoothMPC_FLOAT* smoothMPC_dsubcc30 = smoothMPC_ds_cc + 491;
smoothMPC_FLOAT* smoothMPC_ccrhsub30 = smoothMPC_ccrhs + 491;
smoothMPC_FLOAT smoothMPC_Phi30[11];
smoothMPC_FLOAT smoothMPC_W30[11];
smoothMPC_FLOAT smoothMPC_Ysd30[9];
smoothMPC_FLOAT smoothMPC_Lsd30[9];
smoothMPC_FLOAT* smoothMPC_z31 = smoothMPC_z + 341;
smoothMPC_FLOAT* smoothMPC_dzaff31 = smoothMPC_dz_aff + 341;
smoothMPC_FLOAT* smoothMPC_dzcc31 = smoothMPC_dz_cc + 341;
smoothMPC_FLOAT* smoothMPC_rd31 = smoothMPC_rd + 341;
smoothMPC_FLOAT smoothMPC_Lbyrd31[11];
smoothMPC_FLOAT* smoothMPC_grad_cost31 = smoothMPC_grad_cost + 341;
smoothMPC_FLOAT* smoothMPC_grad_eq31 = smoothMPC_grad_eq + 341;
smoothMPC_FLOAT* smoothMPC_grad_ineq31 = smoothMPC_grad_ineq + 341;
smoothMPC_FLOAT smoothMPC_ctv31[11];
smoothMPC_FLOAT* smoothMPC_v31 = smoothMPC_v + 93;
smoothMPC_FLOAT smoothMPC_re31[3];
smoothMPC_FLOAT smoothMPC_beta31[3];
smoothMPC_FLOAT smoothMPC_betacc31[3];
smoothMPC_FLOAT* smoothMPC_dvaff31 = smoothMPC_dv_aff + 93;
smoothMPC_FLOAT* smoothMPC_dvcc31 = smoothMPC_dv_cc + 93;
smoothMPC_FLOAT smoothMPC_V31[33];
smoothMPC_FLOAT smoothMPC_Yd31[6];
smoothMPC_FLOAT smoothMPC_Ld31[6];
smoothMPC_FLOAT smoothMPC_yy31[3];
smoothMPC_FLOAT smoothMPC_bmy31[3];
int smoothMPC_lbIdx31[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb31 = smoothMPC_l + 496;
smoothMPC_FLOAT* smoothMPC_slb31 = smoothMPC_s + 496;
smoothMPC_FLOAT* smoothMPC_llbbyslb31 = smoothMPC_lbys + 496;
smoothMPC_FLOAT smoothMPC_rilb31[11];
smoothMPC_FLOAT* smoothMPC_dllbaff31 = smoothMPC_dl_aff + 496;
smoothMPC_FLOAT* smoothMPC_dslbaff31 = smoothMPC_ds_aff + 496;
smoothMPC_FLOAT* smoothMPC_dllbcc31 = smoothMPC_dl_cc + 496;
smoothMPC_FLOAT* smoothMPC_dslbcc31 = smoothMPC_ds_cc + 496;
smoothMPC_FLOAT* smoothMPC_ccrhsl31 = smoothMPC_ccrhs + 496;
int smoothMPC_ubIdx31[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub31 = smoothMPC_l + 507;
smoothMPC_FLOAT* smoothMPC_sub31 = smoothMPC_s + 507;
smoothMPC_FLOAT* smoothMPC_lubbysub31 = smoothMPC_lbys + 507;
smoothMPC_FLOAT smoothMPC_riub31[5];
smoothMPC_FLOAT* smoothMPC_dlubaff31 = smoothMPC_dl_aff + 507;
smoothMPC_FLOAT* smoothMPC_dsubaff31 = smoothMPC_ds_aff + 507;
smoothMPC_FLOAT* smoothMPC_dlubcc31 = smoothMPC_dl_cc + 507;
smoothMPC_FLOAT* smoothMPC_dsubcc31 = smoothMPC_ds_cc + 507;
smoothMPC_FLOAT* smoothMPC_ccrhsub31 = smoothMPC_ccrhs + 507;
smoothMPC_FLOAT smoothMPC_Phi31[11];
smoothMPC_FLOAT smoothMPC_W31[11];
smoothMPC_FLOAT smoothMPC_Ysd31[9];
smoothMPC_FLOAT smoothMPC_Lsd31[9];
smoothMPC_FLOAT* smoothMPC_z32 = smoothMPC_z + 352;
smoothMPC_FLOAT* smoothMPC_dzaff32 = smoothMPC_dz_aff + 352;
smoothMPC_FLOAT* smoothMPC_dzcc32 = smoothMPC_dz_cc + 352;
smoothMPC_FLOAT* smoothMPC_rd32 = smoothMPC_rd + 352;
smoothMPC_FLOAT smoothMPC_Lbyrd32[11];
smoothMPC_FLOAT* smoothMPC_grad_cost32 = smoothMPC_grad_cost + 352;
smoothMPC_FLOAT* smoothMPC_grad_eq32 = smoothMPC_grad_eq + 352;
smoothMPC_FLOAT* smoothMPC_grad_ineq32 = smoothMPC_grad_ineq + 352;
smoothMPC_FLOAT smoothMPC_ctv32[11];
smoothMPC_FLOAT* smoothMPC_v32 = smoothMPC_v + 96;
smoothMPC_FLOAT smoothMPC_re32[3];
smoothMPC_FLOAT smoothMPC_beta32[3];
smoothMPC_FLOAT smoothMPC_betacc32[3];
smoothMPC_FLOAT* smoothMPC_dvaff32 = smoothMPC_dv_aff + 96;
smoothMPC_FLOAT* smoothMPC_dvcc32 = smoothMPC_dv_cc + 96;
smoothMPC_FLOAT smoothMPC_V32[33];
smoothMPC_FLOAT smoothMPC_Yd32[6];
smoothMPC_FLOAT smoothMPC_Ld32[6];
smoothMPC_FLOAT smoothMPC_yy32[3];
smoothMPC_FLOAT smoothMPC_bmy32[3];
int smoothMPC_lbIdx32[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb32 = smoothMPC_l + 512;
smoothMPC_FLOAT* smoothMPC_slb32 = smoothMPC_s + 512;
smoothMPC_FLOAT* smoothMPC_llbbyslb32 = smoothMPC_lbys + 512;
smoothMPC_FLOAT smoothMPC_rilb32[11];
smoothMPC_FLOAT* smoothMPC_dllbaff32 = smoothMPC_dl_aff + 512;
smoothMPC_FLOAT* smoothMPC_dslbaff32 = smoothMPC_ds_aff + 512;
smoothMPC_FLOAT* smoothMPC_dllbcc32 = smoothMPC_dl_cc + 512;
smoothMPC_FLOAT* smoothMPC_dslbcc32 = smoothMPC_ds_cc + 512;
smoothMPC_FLOAT* smoothMPC_ccrhsl32 = smoothMPC_ccrhs + 512;
int smoothMPC_ubIdx32[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub32 = smoothMPC_l + 523;
smoothMPC_FLOAT* smoothMPC_sub32 = smoothMPC_s + 523;
smoothMPC_FLOAT* smoothMPC_lubbysub32 = smoothMPC_lbys + 523;
smoothMPC_FLOAT smoothMPC_riub32[5];
smoothMPC_FLOAT* smoothMPC_dlubaff32 = smoothMPC_dl_aff + 523;
smoothMPC_FLOAT* smoothMPC_dsubaff32 = smoothMPC_ds_aff + 523;
smoothMPC_FLOAT* smoothMPC_dlubcc32 = smoothMPC_dl_cc + 523;
smoothMPC_FLOAT* smoothMPC_dsubcc32 = smoothMPC_ds_cc + 523;
smoothMPC_FLOAT* smoothMPC_ccrhsub32 = smoothMPC_ccrhs + 523;
smoothMPC_FLOAT smoothMPC_Phi32[11];
smoothMPC_FLOAT smoothMPC_W32[11];
smoothMPC_FLOAT smoothMPC_Ysd32[9];
smoothMPC_FLOAT smoothMPC_Lsd32[9];
smoothMPC_FLOAT* smoothMPC_z33 = smoothMPC_z + 363;
smoothMPC_FLOAT* smoothMPC_dzaff33 = smoothMPC_dz_aff + 363;
smoothMPC_FLOAT* smoothMPC_dzcc33 = smoothMPC_dz_cc + 363;
smoothMPC_FLOAT* smoothMPC_rd33 = smoothMPC_rd + 363;
smoothMPC_FLOAT smoothMPC_Lbyrd33[11];
smoothMPC_FLOAT* smoothMPC_grad_cost33 = smoothMPC_grad_cost + 363;
smoothMPC_FLOAT* smoothMPC_grad_eq33 = smoothMPC_grad_eq + 363;
smoothMPC_FLOAT* smoothMPC_grad_ineq33 = smoothMPC_grad_ineq + 363;
smoothMPC_FLOAT smoothMPC_ctv33[11];
smoothMPC_FLOAT* smoothMPC_v33 = smoothMPC_v + 99;
smoothMPC_FLOAT smoothMPC_re33[3];
smoothMPC_FLOAT smoothMPC_beta33[3];
smoothMPC_FLOAT smoothMPC_betacc33[3];
smoothMPC_FLOAT* smoothMPC_dvaff33 = smoothMPC_dv_aff + 99;
smoothMPC_FLOAT* smoothMPC_dvcc33 = smoothMPC_dv_cc + 99;
smoothMPC_FLOAT smoothMPC_V33[33];
smoothMPC_FLOAT smoothMPC_Yd33[6];
smoothMPC_FLOAT smoothMPC_Ld33[6];
smoothMPC_FLOAT smoothMPC_yy33[3];
smoothMPC_FLOAT smoothMPC_bmy33[3];
int smoothMPC_lbIdx33[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb33 = smoothMPC_l + 528;
smoothMPC_FLOAT* smoothMPC_slb33 = smoothMPC_s + 528;
smoothMPC_FLOAT* smoothMPC_llbbyslb33 = smoothMPC_lbys + 528;
smoothMPC_FLOAT smoothMPC_rilb33[11];
smoothMPC_FLOAT* smoothMPC_dllbaff33 = smoothMPC_dl_aff + 528;
smoothMPC_FLOAT* smoothMPC_dslbaff33 = smoothMPC_ds_aff + 528;
smoothMPC_FLOAT* smoothMPC_dllbcc33 = smoothMPC_dl_cc + 528;
smoothMPC_FLOAT* smoothMPC_dslbcc33 = smoothMPC_ds_cc + 528;
smoothMPC_FLOAT* smoothMPC_ccrhsl33 = smoothMPC_ccrhs + 528;
int smoothMPC_ubIdx33[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub33 = smoothMPC_l + 539;
smoothMPC_FLOAT* smoothMPC_sub33 = smoothMPC_s + 539;
smoothMPC_FLOAT* smoothMPC_lubbysub33 = smoothMPC_lbys + 539;
smoothMPC_FLOAT smoothMPC_riub33[5];
smoothMPC_FLOAT* smoothMPC_dlubaff33 = smoothMPC_dl_aff + 539;
smoothMPC_FLOAT* smoothMPC_dsubaff33 = smoothMPC_ds_aff + 539;
smoothMPC_FLOAT* smoothMPC_dlubcc33 = smoothMPC_dl_cc + 539;
smoothMPC_FLOAT* smoothMPC_dsubcc33 = smoothMPC_ds_cc + 539;
smoothMPC_FLOAT* smoothMPC_ccrhsub33 = smoothMPC_ccrhs + 539;
smoothMPC_FLOAT smoothMPC_Phi33[11];
smoothMPC_FLOAT smoothMPC_W33[11];
smoothMPC_FLOAT smoothMPC_Ysd33[9];
smoothMPC_FLOAT smoothMPC_Lsd33[9];
smoothMPC_FLOAT* smoothMPC_z34 = smoothMPC_z + 374;
smoothMPC_FLOAT* smoothMPC_dzaff34 = smoothMPC_dz_aff + 374;
smoothMPC_FLOAT* smoothMPC_dzcc34 = smoothMPC_dz_cc + 374;
smoothMPC_FLOAT* smoothMPC_rd34 = smoothMPC_rd + 374;
smoothMPC_FLOAT smoothMPC_Lbyrd34[11];
smoothMPC_FLOAT* smoothMPC_grad_cost34 = smoothMPC_grad_cost + 374;
smoothMPC_FLOAT* smoothMPC_grad_eq34 = smoothMPC_grad_eq + 374;
smoothMPC_FLOAT* smoothMPC_grad_ineq34 = smoothMPC_grad_ineq + 374;
smoothMPC_FLOAT smoothMPC_ctv34[11];
smoothMPC_FLOAT* smoothMPC_v34 = smoothMPC_v + 102;
smoothMPC_FLOAT smoothMPC_re34[3];
smoothMPC_FLOAT smoothMPC_beta34[3];
smoothMPC_FLOAT smoothMPC_betacc34[3];
smoothMPC_FLOAT* smoothMPC_dvaff34 = smoothMPC_dv_aff + 102;
smoothMPC_FLOAT* smoothMPC_dvcc34 = smoothMPC_dv_cc + 102;
smoothMPC_FLOAT smoothMPC_V34[33];
smoothMPC_FLOAT smoothMPC_Yd34[6];
smoothMPC_FLOAT smoothMPC_Ld34[6];
smoothMPC_FLOAT smoothMPC_yy34[3];
smoothMPC_FLOAT smoothMPC_bmy34[3];
int smoothMPC_lbIdx34[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb34 = smoothMPC_l + 544;
smoothMPC_FLOAT* smoothMPC_slb34 = smoothMPC_s + 544;
smoothMPC_FLOAT* smoothMPC_llbbyslb34 = smoothMPC_lbys + 544;
smoothMPC_FLOAT smoothMPC_rilb34[11];
smoothMPC_FLOAT* smoothMPC_dllbaff34 = smoothMPC_dl_aff + 544;
smoothMPC_FLOAT* smoothMPC_dslbaff34 = smoothMPC_ds_aff + 544;
smoothMPC_FLOAT* smoothMPC_dllbcc34 = smoothMPC_dl_cc + 544;
smoothMPC_FLOAT* smoothMPC_dslbcc34 = smoothMPC_ds_cc + 544;
smoothMPC_FLOAT* smoothMPC_ccrhsl34 = smoothMPC_ccrhs + 544;
int smoothMPC_ubIdx34[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub34 = smoothMPC_l + 555;
smoothMPC_FLOAT* smoothMPC_sub34 = smoothMPC_s + 555;
smoothMPC_FLOAT* smoothMPC_lubbysub34 = smoothMPC_lbys + 555;
smoothMPC_FLOAT smoothMPC_riub34[5];
smoothMPC_FLOAT* smoothMPC_dlubaff34 = smoothMPC_dl_aff + 555;
smoothMPC_FLOAT* smoothMPC_dsubaff34 = smoothMPC_ds_aff + 555;
smoothMPC_FLOAT* smoothMPC_dlubcc34 = smoothMPC_dl_cc + 555;
smoothMPC_FLOAT* smoothMPC_dsubcc34 = smoothMPC_ds_cc + 555;
smoothMPC_FLOAT* smoothMPC_ccrhsub34 = smoothMPC_ccrhs + 555;
smoothMPC_FLOAT smoothMPC_Phi34[11];
smoothMPC_FLOAT smoothMPC_W34[11];
smoothMPC_FLOAT smoothMPC_Ysd34[9];
smoothMPC_FLOAT smoothMPC_Lsd34[9];
smoothMPC_FLOAT* smoothMPC_z35 = smoothMPC_z + 385;
smoothMPC_FLOAT* smoothMPC_dzaff35 = smoothMPC_dz_aff + 385;
smoothMPC_FLOAT* smoothMPC_dzcc35 = smoothMPC_dz_cc + 385;
smoothMPC_FLOAT* smoothMPC_rd35 = smoothMPC_rd + 385;
smoothMPC_FLOAT smoothMPC_Lbyrd35[11];
smoothMPC_FLOAT* smoothMPC_grad_cost35 = smoothMPC_grad_cost + 385;
smoothMPC_FLOAT* smoothMPC_grad_eq35 = smoothMPC_grad_eq + 385;
smoothMPC_FLOAT* smoothMPC_grad_ineq35 = smoothMPC_grad_ineq + 385;
smoothMPC_FLOAT smoothMPC_ctv35[11];
smoothMPC_FLOAT* smoothMPC_v35 = smoothMPC_v + 105;
smoothMPC_FLOAT smoothMPC_re35[3];
smoothMPC_FLOAT smoothMPC_beta35[3];
smoothMPC_FLOAT smoothMPC_betacc35[3];
smoothMPC_FLOAT* smoothMPC_dvaff35 = smoothMPC_dv_aff + 105;
smoothMPC_FLOAT* smoothMPC_dvcc35 = smoothMPC_dv_cc + 105;
smoothMPC_FLOAT smoothMPC_V35[33];
smoothMPC_FLOAT smoothMPC_Yd35[6];
smoothMPC_FLOAT smoothMPC_Ld35[6];
smoothMPC_FLOAT smoothMPC_yy35[3];
smoothMPC_FLOAT smoothMPC_bmy35[3];
int smoothMPC_lbIdx35[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb35 = smoothMPC_l + 560;
smoothMPC_FLOAT* smoothMPC_slb35 = smoothMPC_s + 560;
smoothMPC_FLOAT* smoothMPC_llbbyslb35 = smoothMPC_lbys + 560;
smoothMPC_FLOAT smoothMPC_rilb35[11];
smoothMPC_FLOAT* smoothMPC_dllbaff35 = smoothMPC_dl_aff + 560;
smoothMPC_FLOAT* smoothMPC_dslbaff35 = smoothMPC_ds_aff + 560;
smoothMPC_FLOAT* smoothMPC_dllbcc35 = smoothMPC_dl_cc + 560;
smoothMPC_FLOAT* smoothMPC_dslbcc35 = smoothMPC_ds_cc + 560;
smoothMPC_FLOAT* smoothMPC_ccrhsl35 = smoothMPC_ccrhs + 560;
int smoothMPC_ubIdx35[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub35 = smoothMPC_l + 571;
smoothMPC_FLOAT* smoothMPC_sub35 = smoothMPC_s + 571;
smoothMPC_FLOAT* smoothMPC_lubbysub35 = smoothMPC_lbys + 571;
smoothMPC_FLOAT smoothMPC_riub35[5];
smoothMPC_FLOAT* smoothMPC_dlubaff35 = smoothMPC_dl_aff + 571;
smoothMPC_FLOAT* smoothMPC_dsubaff35 = smoothMPC_ds_aff + 571;
smoothMPC_FLOAT* smoothMPC_dlubcc35 = smoothMPC_dl_cc + 571;
smoothMPC_FLOAT* smoothMPC_dsubcc35 = smoothMPC_ds_cc + 571;
smoothMPC_FLOAT* smoothMPC_ccrhsub35 = smoothMPC_ccrhs + 571;
smoothMPC_FLOAT smoothMPC_Phi35[11];
smoothMPC_FLOAT smoothMPC_W35[11];
smoothMPC_FLOAT smoothMPC_Ysd35[9];
smoothMPC_FLOAT smoothMPC_Lsd35[9];
smoothMPC_FLOAT* smoothMPC_z36 = smoothMPC_z + 396;
smoothMPC_FLOAT* smoothMPC_dzaff36 = smoothMPC_dz_aff + 396;
smoothMPC_FLOAT* smoothMPC_dzcc36 = smoothMPC_dz_cc + 396;
smoothMPC_FLOAT* smoothMPC_rd36 = smoothMPC_rd + 396;
smoothMPC_FLOAT smoothMPC_Lbyrd36[11];
smoothMPC_FLOAT* smoothMPC_grad_cost36 = smoothMPC_grad_cost + 396;
smoothMPC_FLOAT* smoothMPC_grad_eq36 = smoothMPC_grad_eq + 396;
smoothMPC_FLOAT* smoothMPC_grad_ineq36 = smoothMPC_grad_ineq + 396;
smoothMPC_FLOAT smoothMPC_ctv36[11];
smoothMPC_FLOAT* smoothMPC_v36 = smoothMPC_v + 108;
smoothMPC_FLOAT smoothMPC_re36[3];
smoothMPC_FLOAT smoothMPC_beta36[3];
smoothMPC_FLOAT smoothMPC_betacc36[3];
smoothMPC_FLOAT* smoothMPC_dvaff36 = smoothMPC_dv_aff + 108;
smoothMPC_FLOAT* smoothMPC_dvcc36 = smoothMPC_dv_cc + 108;
smoothMPC_FLOAT smoothMPC_V36[33];
smoothMPC_FLOAT smoothMPC_Yd36[6];
smoothMPC_FLOAT smoothMPC_Ld36[6];
smoothMPC_FLOAT smoothMPC_yy36[3];
smoothMPC_FLOAT smoothMPC_bmy36[3];
int smoothMPC_lbIdx36[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb36 = smoothMPC_l + 576;
smoothMPC_FLOAT* smoothMPC_slb36 = smoothMPC_s + 576;
smoothMPC_FLOAT* smoothMPC_llbbyslb36 = smoothMPC_lbys + 576;
smoothMPC_FLOAT smoothMPC_rilb36[11];
smoothMPC_FLOAT* smoothMPC_dllbaff36 = smoothMPC_dl_aff + 576;
smoothMPC_FLOAT* smoothMPC_dslbaff36 = smoothMPC_ds_aff + 576;
smoothMPC_FLOAT* smoothMPC_dllbcc36 = smoothMPC_dl_cc + 576;
smoothMPC_FLOAT* smoothMPC_dslbcc36 = smoothMPC_ds_cc + 576;
smoothMPC_FLOAT* smoothMPC_ccrhsl36 = smoothMPC_ccrhs + 576;
int smoothMPC_ubIdx36[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub36 = smoothMPC_l + 587;
smoothMPC_FLOAT* smoothMPC_sub36 = smoothMPC_s + 587;
smoothMPC_FLOAT* smoothMPC_lubbysub36 = smoothMPC_lbys + 587;
smoothMPC_FLOAT smoothMPC_riub36[5];
smoothMPC_FLOAT* smoothMPC_dlubaff36 = smoothMPC_dl_aff + 587;
smoothMPC_FLOAT* smoothMPC_dsubaff36 = smoothMPC_ds_aff + 587;
smoothMPC_FLOAT* smoothMPC_dlubcc36 = smoothMPC_dl_cc + 587;
smoothMPC_FLOAT* smoothMPC_dsubcc36 = smoothMPC_ds_cc + 587;
smoothMPC_FLOAT* smoothMPC_ccrhsub36 = smoothMPC_ccrhs + 587;
smoothMPC_FLOAT smoothMPC_Phi36[11];
smoothMPC_FLOAT smoothMPC_W36[11];
smoothMPC_FLOAT smoothMPC_Ysd36[9];
smoothMPC_FLOAT smoothMPC_Lsd36[9];
smoothMPC_FLOAT* smoothMPC_z37 = smoothMPC_z + 407;
smoothMPC_FLOAT* smoothMPC_dzaff37 = smoothMPC_dz_aff + 407;
smoothMPC_FLOAT* smoothMPC_dzcc37 = smoothMPC_dz_cc + 407;
smoothMPC_FLOAT* smoothMPC_rd37 = smoothMPC_rd + 407;
smoothMPC_FLOAT smoothMPC_Lbyrd37[11];
smoothMPC_FLOAT* smoothMPC_grad_cost37 = smoothMPC_grad_cost + 407;
smoothMPC_FLOAT* smoothMPC_grad_eq37 = smoothMPC_grad_eq + 407;
smoothMPC_FLOAT* smoothMPC_grad_ineq37 = smoothMPC_grad_ineq + 407;
smoothMPC_FLOAT smoothMPC_ctv37[11];
smoothMPC_FLOAT* smoothMPC_v37 = smoothMPC_v + 111;
smoothMPC_FLOAT smoothMPC_re37[3];
smoothMPC_FLOAT smoothMPC_beta37[3];
smoothMPC_FLOAT smoothMPC_betacc37[3];
smoothMPC_FLOAT* smoothMPC_dvaff37 = smoothMPC_dv_aff + 111;
smoothMPC_FLOAT* smoothMPC_dvcc37 = smoothMPC_dv_cc + 111;
smoothMPC_FLOAT smoothMPC_V37[33];
smoothMPC_FLOAT smoothMPC_Yd37[6];
smoothMPC_FLOAT smoothMPC_Ld37[6];
smoothMPC_FLOAT smoothMPC_yy37[3];
smoothMPC_FLOAT smoothMPC_bmy37[3];
int smoothMPC_lbIdx37[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb37 = smoothMPC_l + 592;
smoothMPC_FLOAT* smoothMPC_slb37 = smoothMPC_s + 592;
smoothMPC_FLOAT* smoothMPC_llbbyslb37 = smoothMPC_lbys + 592;
smoothMPC_FLOAT smoothMPC_rilb37[11];
smoothMPC_FLOAT* smoothMPC_dllbaff37 = smoothMPC_dl_aff + 592;
smoothMPC_FLOAT* smoothMPC_dslbaff37 = smoothMPC_ds_aff + 592;
smoothMPC_FLOAT* smoothMPC_dllbcc37 = smoothMPC_dl_cc + 592;
smoothMPC_FLOAT* smoothMPC_dslbcc37 = smoothMPC_ds_cc + 592;
smoothMPC_FLOAT* smoothMPC_ccrhsl37 = smoothMPC_ccrhs + 592;
int smoothMPC_ubIdx37[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub37 = smoothMPC_l + 603;
smoothMPC_FLOAT* smoothMPC_sub37 = smoothMPC_s + 603;
smoothMPC_FLOAT* smoothMPC_lubbysub37 = smoothMPC_lbys + 603;
smoothMPC_FLOAT smoothMPC_riub37[5];
smoothMPC_FLOAT* smoothMPC_dlubaff37 = smoothMPC_dl_aff + 603;
smoothMPC_FLOAT* smoothMPC_dsubaff37 = smoothMPC_ds_aff + 603;
smoothMPC_FLOAT* smoothMPC_dlubcc37 = smoothMPC_dl_cc + 603;
smoothMPC_FLOAT* smoothMPC_dsubcc37 = smoothMPC_ds_cc + 603;
smoothMPC_FLOAT* smoothMPC_ccrhsub37 = smoothMPC_ccrhs + 603;
smoothMPC_FLOAT smoothMPC_Phi37[11];
smoothMPC_FLOAT smoothMPC_W37[11];
smoothMPC_FLOAT smoothMPC_Ysd37[9];
smoothMPC_FLOAT smoothMPC_Lsd37[9];
smoothMPC_FLOAT* smoothMPC_z38 = smoothMPC_z + 418;
smoothMPC_FLOAT* smoothMPC_dzaff38 = smoothMPC_dz_aff + 418;
smoothMPC_FLOAT* smoothMPC_dzcc38 = smoothMPC_dz_cc + 418;
smoothMPC_FLOAT* smoothMPC_rd38 = smoothMPC_rd + 418;
smoothMPC_FLOAT smoothMPC_Lbyrd38[11];
smoothMPC_FLOAT* smoothMPC_grad_cost38 = smoothMPC_grad_cost + 418;
smoothMPC_FLOAT* smoothMPC_grad_eq38 = smoothMPC_grad_eq + 418;
smoothMPC_FLOAT* smoothMPC_grad_ineq38 = smoothMPC_grad_ineq + 418;
smoothMPC_FLOAT smoothMPC_ctv38[11];
smoothMPC_FLOAT* smoothMPC_v38 = smoothMPC_v + 114;
smoothMPC_FLOAT smoothMPC_re38[3];
smoothMPC_FLOAT smoothMPC_beta38[3];
smoothMPC_FLOAT smoothMPC_betacc38[3];
smoothMPC_FLOAT* smoothMPC_dvaff38 = smoothMPC_dv_aff + 114;
smoothMPC_FLOAT* smoothMPC_dvcc38 = smoothMPC_dv_cc + 114;
smoothMPC_FLOAT smoothMPC_V38[33];
smoothMPC_FLOAT smoothMPC_Yd38[6];
smoothMPC_FLOAT smoothMPC_Ld38[6];
smoothMPC_FLOAT smoothMPC_yy38[3];
smoothMPC_FLOAT smoothMPC_bmy38[3];
int smoothMPC_lbIdx38[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb38 = smoothMPC_l + 608;
smoothMPC_FLOAT* smoothMPC_slb38 = smoothMPC_s + 608;
smoothMPC_FLOAT* smoothMPC_llbbyslb38 = smoothMPC_lbys + 608;
smoothMPC_FLOAT smoothMPC_rilb38[11];
smoothMPC_FLOAT* smoothMPC_dllbaff38 = smoothMPC_dl_aff + 608;
smoothMPC_FLOAT* smoothMPC_dslbaff38 = smoothMPC_ds_aff + 608;
smoothMPC_FLOAT* smoothMPC_dllbcc38 = smoothMPC_dl_cc + 608;
smoothMPC_FLOAT* smoothMPC_dslbcc38 = smoothMPC_ds_cc + 608;
smoothMPC_FLOAT* smoothMPC_ccrhsl38 = smoothMPC_ccrhs + 608;
int smoothMPC_ubIdx38[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub38 = smoothMPC_l + 619;
smoothMPC_FLOAT* smoothMPC_sub38 = smoothMPC_s + 619;
smoothMPC_FLOAT* smoothMPC_lubbysub38 = smoothMPC_lbys + 619;
smoothMPC_FLOAT smoothMPC_riub38[5];
smoothMPC_FLOAT* smoothMPC_dlubaff38 = smoothMPC_dl_aff + 619;
smoothMPC_FLOAT* smoothMPC_dsubaff38 = smoothMPC_ds_aff + 619;
smoothMPC_FLOAT* smoothMPC_dlubcc38 = smoothMPC_dl_cc + 619;
smoothMPC_FLOAT* smoothMPC_dsubcc38 = smoothMPC_ds_cc + 619;
smoothMPC_FLOAT* smoothMPC_ccrhsub38 = smoothMPC_ccrhs + 619;
smoothMPC_FLOAT smoothMPC_Phi38[11];
smoothMPC_FLOAT smoothMPC_W38[11];
smoothMPC_FLOAT smoothMPC_Ysd38[9];
smoothMPC_FLOAT smoothMPC_Lsd38[9];
smoothMPC_FLOAT* smoothMPC_z39 = smoothMPC_z + 429;
smoothMPC_FLOAT* smoothMPC_dzaff39 = smoothMPC_dz_aff + 429;
smoothMPC_FLOAT* smoothMPC_dzcc39 = smoothMPC_dz_cc + 429;
smoothMPC_FLOAT* smoothMPC_rd39 = smoothMPC_rd + 429;
smoothMPC_FLOAT smoothMPC_Lbyrd39[11];
smoothMPC_FLOAT* smoothMPC_grad_cost39 = smoothMPC_grad_cost + 429;
smoothMPC_FLOAT* smoothMPC_grad_eq39 = smoothMPC_grad_eq + 429;
smoothMPC_FLOAT* smoothMPC_grad_ineq39 = smoothMPC_grad_ineq + 429;
smoothMPC_FLOAT smoothMPC_ctv39[11];
smoothMPC_FLOAT* smoothMPC_v39 = smoothMPC_v + 117;
smoothMPC_FLOAT smoothMPC_re39[3];
smoothMPC_FLOAT smoothMPC_beta39[3];
smoothMPC_FLOAT smoothMPC_betacc39[3];
smoothMPC_FLOAT* smoothMPC_dvaff39 = smoothMPC_dv_aff + 117;
smoothMPC_FLOAT* smoothMPC_dvcc39 = smoothMPC_dv_cc + 117;
smoothMPC_FLOAT smoothMPC_V39[33];
smoothMPC_FLOAT smoothMPC_Yd39[6];
smoothMPC_FLOAT smoothMPC_Ld39[6];
smoothMPC_FLOAT smoothMPC_yy39[3];
smoothMPC_FLOAT smoothMPC_bmy39[3];
int smoothMPC_lbIdx39[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb39 = smoothMPC_l + 624;
smoothMPC_FLOAT* smoothMPC_slb39 = smoothMPC_s + 624;
smoothMPC_FLOAT* smoothMPC_llbbyslb39 = smoothMPC_lbys + 624;
smoothMPC_FLOAT smoothMPC_rilb39[11];
smoothMPC_FLOAT* smoothMPC_dllbaff39 = smoothMPC_dl_aff + 624;
smoothMPC_FLOAT* smoothMPC_dslbaff39 = smoothMPC_ds_aff + 624;
smoothMPC_FLOAT* smoothMPC_dllbcc39 = smoothMPC_dl_cc + 624;
smoothMPC_FLOAT* smoothMPC_dslbcc39 = smoothMPC_ds_cc + 624;
smoothMPC_FLOAT* smoothMPC_ccrhsl39 = smoothMPC_ccrhs + 624;
int smoothMPC_ubIdx39[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub39 = smoothMPC_l + 635;
smoothMPC_FLOAT* smoothMPC_sub39 = smoothMPC_s + 635;
smoothMPC_FLOAT* smoothMPC_lubbysub39 = smoothMPC_lbys + 635;
smoothMPC_FLOAT smoothMPC_riub39[5];
smoothMPC_FLOAT* smoothMPC_dlubaff39 = smoothMPC_dl_aff + 635;
smoothMPC_FLOAT* smoothMPC_dsubaff39 = smoothMPC_ds_aff + 635;
smoothMPC_FLOAT* smoothMPC_dlubcc39 = smoothMPC_dl_cc + 635;
smoothMPC_FLOAT* smoothMPC_dsubcc39 = smoothMPC_ds_cc + 635;
smoothMPC_FLOAT* smoothMPC_ccrhsub39 = smoothMPC_ccrhs + 635;
smoothMPC_FLOAT smoothMPC_Phi39[11];
smoothMPC_FLOAT smoothMPC_W39[11];
smoothMPC_FLOAT smoothMPC_Ysd39[9];
smoothMPC_FLOAT smoothMPC_Lsd39[9];
smoothMPC_FLOAT* smoothMPC_z40 = smoothMPC_z + 440;
smoothMPC_FLOAT* smoothMPC_dzaff40 = smoothMPC_dz_aff + 440;
smoothMPC_FLOAT* smoothMPC_dzcc40 = smoothMPC_dz_cc + 440;
smoothMPC_FLOAT* smoothMPC_rd40 = smoothMPC_rd + 440;
smoothMPC_FLOAT smoothMPC_Lbyrd40[11];
smoothMPC_FLOAT* smoothMPC_grad_cost40 = smoothMPC_grad_cost + 440;
smoothMPC_FLOAT* smoothMPC_grad_eq40 = smoothMPC_grad_eq + 440;
smoothMPC_FLOAT* smoothMPC_grad_ineq40 = smoothMPC_grad_ineq + 440;
smoothMPC_FLOAT smoothMPC_ctv40[11];
smoothMPC_FLOAT* smoothMPC_v40 = smoothMPC_v + 120;
smoothMPC_FLOAT smoothMPC_re40[3];
smoothMPC_FLOAT smoothMPC_beta40[3];
smoothMPC_FLOAT smoothMPC_betacc40[3];
smoothMPC_FLOAT* smoothMPC_dvaff40 = smoothMPC_dv_aff + 120;
smoothMPC_FLOAT* smoothMPC_dvcc40 = smoothMPC_dv_cc + 120;
smoothMPC_FLOAT smoothMPC_V40[33];
smoothMPC_FLOAT smoothMPC_Yd40[6];
smoothMPC_FLOAT smoothMPC_Ld40[6];
smoothMPC_FLOAT smoothMPC_yy40[3];
smoothMPC_FLOAT smoothMPC_bmy40[3];
int smoothMPC_lbIdx40[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb40 = smoothMPC_l + 640;
smoothMPC_FLOAT* smoothMPC_slb40 = smoothMPC_s + 640;
smoothMPC_FLOAT* smoothMPC_llbbyslb40 = smoothMPC_lbys + 640;
smoothMPC_FLOAT smoothMPC_rilb40[11];
smoothMPC_FLOAT* smoothMPC_dllbaff40 = smoothMPC_dl_aff + 640;
smoothMPC_FLOAT* smoothMPC_dslbaff40 = smoothMPC_ds_aff + 640;
smoothMPC_FLOAT* smoothMPC_dllbcc40 = smoothMPC_dl_cc + 640;
smoothMPC_FLOAT* smoothMPC_dslbcc40 = smoothMPC_ds_cc + 640;
smoothMPC_FLOAT* smoothMPC_ccrhsl40 = smoothMPC_ccrhs + 640;
int smoothMPC_ubIdx40[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub40 = smoothMPC_l + 651;
smoothMPC_FLOAT* smoothMPC_sub40 = smoothMPC_s + 651;
smoothMPC_FLOAT* smoothMPC_lubbysub40 = smoothMPC_lbys + 651;
smoothMPC_FLOAT smoothMPC_riub40[5];
smoothMPC_FLOAT* smoothMPC_dlubaff40 = smoothMPC_dl_aff + 651;
smoothMPC_FLOAT* smoothMPC_dsubaff40 = smoothMPC_ds_aff + 651;
smoothMPC_FLOAT* smoothMPC_dlubcc40 = smoothMPC_dl_cc + 651;
smoothMPC_FLOAT* smoothMPC_dsubcc40 = smoothMPC_ds_cc + 651;
smoothMPC_FLOAT* smoothMPC_ccrhsub40 = smoothMPC_ccrhs + 651;
smoothMPC_FLOAT smoothMPC_Phi40[11];
smoothMPC_FLOAT smoothMPC_W40[11];
smoothMPC_FLOAT smoothMPC_Ysd40[9];
smoothMPC_FLOAT smoothMPC_Lsd40[9];
smoothMPC_FLOAT* smoothMPC_z41 = smoothMPC_z + 451;
smoothMPC_FLOAT* smoothMPC_dzaff41 = smoothMPC_dz_aff + 451;
smoothMPC_FLOAT* smoothMPC_dzcc41 = smoothMPC_dz_cc + 451;
smoothMPC_FLOAT* smoothMPC_rd41 = smoothMPC_rd + 451;
smoothMPC_FLOAT smoothMPC_Lbyrd41[11];
smoothMPC_FLOAT* smoothMPC_grad_cost41 = smoothMPC_grad_cost + 451;
smoothMPC_FLOAT* smoothMPC_grad_eq41 = smoothMPC_grad_eq + 451;
smoothMPC_FLOAT* smoothMPC_grad_ineq41 = smoothMPC_grad_ineq + 451;
smoothMPC_FLOAT smoothMPC_ctv41[11];
smoothMPC_FLOAT* smoothMPC_v41 = smoothMPC_v + 123;
smoothMPC_FLOAT smoothMPC_re41[3];
smoothMPC_FLOAT smoothMPC_beta41[3];
smoothMPC_FLOAT smoothMPC_betacc41[3];
smoothMPC_FLOAT* smoothMPC_dvaff41 = smoothMPC_dv_aff + 123;
smoothMPC_FLOAT* smoothMPC_dvcc41 = smoothMPC_dv_cc + 123;
smoothMPC_FLOAT smoothMPC_V41[33];
smoothMPC_FLOAT smoothMPC_Yd41[6];
smoothMPC_FLOAT smoothMPC_Ld41[6];
smoothMPC_FLOAT smoothMPC_yy41[3];
smoothMPC_FLOAT smoothMPC_bmy41[3];
int smoothMPC_lbIdx41[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb41 = smoothMPC_l + 656;
smoothMPC_FLOAT* smoothMPC_slb41 = smoothMPC_s + 656;
smoothMPC_FLOAT* smoothMPC_llbbyslb41 = smoothMPC_lbys + 656;
smoothMPC_FLOAT smoothMPC_rilb41[11];
smoothMPC_FLOAT* smoothMPC_dllbaff41 = smoothMPC_dl_aff + 656;
smoothMPC_FLOAT* smoothMPC_dslbaff41 = smoothMPC_ds_aff + 656;
smoothMPC_FLOAT* smoothMPC_dllbcc41 = smoothMPC_dl_cc + 656;
smoothMPC_FLOAT* smoothMPC_dslbcc41 = smoothMPC_ds_cc + 656;
smoothMPC_FLOAT* smoothMPC_ccrhsl41 = smoothMPC_ccrhs + 656;
int smoothMPC_ubIdx41[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub41 = smoothMPC_l + 667;
smoothMPC_FLOAT* smoothMPC_sub41 = smoothMPC_s + 667;
smoothMPC_FLOAT* smoothMPC_lubbysub41 = smoothMPC_lbys + 667;
smoothMPC_FLOAT smoothMPC_riub41[5];
smoothMPC_FLOAT* smoothMPC_dlubaff41 = smoothMPC_dl_aff + 667;
smoothMPC_FLOAT* smoothMPC_dsubaff41 = smoothMPC_ds_aff + 667;
smoothMPC_FLOAT* smoothMPC_dlubcc41 = smoothMPC_dl_cc + 667;
smoothMPC_FLOAT* smoothMPC_dsubcc41 = smoothMPC_ds_cc + 667;
smoothMPC_FLOAT* smoothMPC_ccrhsub41 = smoothMPC_ccrhs + 667;
smoothMPC_FLOAT smoothMPC_Phi41[11];
smoothMPC_FLOAT smoothMPC_W41[11];
smoothMPC_FLOAT smoothMPC_Ysd41[9];
smoothMPC_FLOAT smoothMPC_Lsd41[9];
smoothMPC_FLOAT* smoothMPC_z42 = smoothMPC_z + 462;
smoothMPC_FLOAT* smoothMPC_dzaff42 = smoothMPC_dz_aff + 462;
smoothMPC_FLOAT* smoothMPC_dzcc42 = smoothMPC_dz_cc + 462;
smoothMPC_FLOAT* smoothMPC_rd42 = smoothMPC_rd + 462;
smoothMPC_FLOAT smoothMPC_Lbyrd42[11];
smoothMPC_FLOAT* smoothMPC_grad_cost42 = smoothMPC_grad_cost + 462;
smoothMPC_FLOAT* smoothMPC_grad_eq42 = smoothMPC_grad_eq + 462;
smoothMPC_FLOAT* smoothMPC_grad_ineq42 = smoothMPC_grad_ineq + 462;
smoothMPC_FLOAT smoothMPC_ctv42[11];
smoothMPC_FLOAT* smoothMPC_v42 = smoothMPC_v + 126;
smoothMPC_FLOAT smoothMPC_re42[3];
smoothMPC_FLOAT smoothMPC_beta42[3];
smoothMPC_FLOAT smoothMPC_betacc42[3];
smoothMPC_FLOAT* smoothMPC_dvaff42 = smoothMPC_dv_aff + 126;
smoothMPC_FLOAT* smoothMPC_dvcc42 = smoothMPC_dv_cc + 126;
smoothMPC_FLOAT smoothMPC_V42[33];
smoothMPC_FLOAT smoothMPC_Yd42[6];
smoothMPC_FLOAT smoothMPC_Ld42[6];
smoothMPC_FLOAT smoothMPC_yy42[3];
smoothMPC_FLOAT smoothMPC_bmy42[3];
int smoothMPC_lbIdx42[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb42 = smoothMPC_l + 672;
smoothMPC_FLOAT* smoothMPC_slb42 = smoothMPC_s + 672;
smoothMPC_FLOAT* smoothMPC_llbbyslb42 = smoothMPC_lbys + 672;
smoothMPC_FLOAT smoothMPC_rilb42[11];
smoothMPC_FLOAT* smoothMPC_dllbaff42 = smoothMPC_dl_aff + 672;
smoothMPC_FLOAT* smoothMPC_dslbaff42 = smoothMPC_ds_aff + 672;
smoothMPC_FLOAT* smoothMPC_dllbcc42 = smoothMPC_dl_cc + 672;
smoothMPC_FLOAT* smoothMPC_dslbcc42 = smoothMPC_ds_cc + 672;
smoothMPC_FLOAT* smoothMPC_ccrhsl42 = smoothMPC_ccrhs + 672;
int smoothMPC_ubIdx42[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub42 = smoothMPC_l + 683;
smoothMPC_FLOAT* smoothMPC_sub42 = smoothMPC_s + 683;
smoothMPC_FLOAT* smoothMPC_lubbysub42 = smoothMPC_lbys + 683;
smoothMPC_FLOAT smoothMPC_riub42[5];
smoothMPC_FLOAT* smoothMPC_dlubaff42 = smoothMPC_dl_aff + 683;
smoothMPC_FLOAT* smoothMPC_dsubaff42 = smoothMPC_ds_aff + 683;
smoothMPC_FLOAT* smoothMPC_dlubcc42 = smoothMPC_dl_cc + 683;
smoothMPC_FLOAT* smoothMPC_dsubcc42 = smoothMPC_ds_cc + 683;
smoothMPC_FLOAT* smoothMPC_ccrhsub42 = smoothMPC_ccrhs + 683;
smoothMPC_FLOAT smoothMPC_Phi42[11];
smoothMPC_FLOAT smoothMPC_W42[11];
smoothMPC_FLOAT smoothMPC_Ysd42[9];
smoothMPC_FLOAT smoothMPC_Lsd42[9];
smoothMPC_FLOAT* smoothMPC_z43 = smoothMPC_z + 473;
smoothMPC_FLOAT* smoothMPC_dzaff43 = smoothMPC_dz_aff + 473;
smoothMPC_FLOAT* smoothMPC_dzcc43 = smoothMPC_dz_cc + 473;
smoothMPC_FLOAT* smoothMPC_rd43 = smoothMPC_rd + 473;
smoothMPC_FLOAT smoothMPC_Lbyrd43[11];
smoothMPC_FLOAT* smoothMPC_grad_cost43 = smoothMPC_grad_cost + 473;
smoothMPC_FLOAT* smoothMPC_grad_eq43 = smoothMPC_grad_eq + 473;
smoothMPC_FLOAT* smoothMPC_grad_ineq43 = smoothMPC_grad_ineq + 473;
smoothMPC_FLOAT smoothMPC_ctv43[11];
smoothMPC_FLOAT* smoothMPC_v43 = smoothMPC_v + 129;
smoothMPC_FLOAT smoothMPC_re43[3];
smoothMPC_FLOAT smoothMPC_beta43[3];
smoothMPC_FLOAT smoothMPC_betacc43[3];
smoothMPC_FLOAT* smoothMPC_dvaff43 = smoothMPC_dv_aff + 129;
smoothMPC_FLOAT* smoothMPC_dvcc43 = smoothMPC_dv_cc + 129;
smoothMPC_FLOAT smoothMPC_V43[33];
smoothMPC_FLOAT smoothMPC_Yd43[6];
smoothMPC_FLOAT smoothMPC_Ld43[6];
smoothMPC_FLOAT smoothMPC_yy43[3];
smoothMPC_FLOAT smoothMPC_bmy43[3];
int smoothMPC_lbIdx43[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb43 = smoothMPC_l + 688;
smoothMPC_FLOAT* smoothMPC_slb43 = smoothMPC_s + 688;
smoothMPC_FLOAT* smoothMPC_llbbyslb43 = smoothMPC_lbys + 688;
smoothMPC_FLOAT smoothMPC_rilb43[11];
smoothMPC_FLOAT* smoothMPC_dllbaff43 = smoothMPC_dl_aff + 688;
smoothMPC_FLOAT* smoothMPC_dslbaff43 = smoothMPC_ds_aff + 688;
smoothMPC_FLOAT* smoothMPC_dllbcc43 = smoothMPC_dl_cc + 688;
smoothMPC_FLOAT* smoothMPC_dslbcc43 = smoothMPC_ds_cc + 688;
smoothMPC_FLOAT* smoothMPC_ccrhsl43 = smoothMPC_ccrhs + 688;
int smoothMPC_ubIdx43[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub43 = smoothMPC_l + 699;
smoothMPC_FLOAT* smoothMPC_sub43 = smoothMPC_s + 699;
smoothMPC_FLOAT* smoothMPC_lubbysub43 = smoothMPC_lbys + 699;
smoothMPC_FLOAT smoothMPC_riub43[5];
smoothMPC_FLOAT* smoothMPC_dlubaff43 = smoothMPC_dl_aff + 699;
smoothMPC_FLOAT* smoothMPC_dsubaff43 = smoothMPC_ds_aff + 699;
smoothMPC_FLOAT* smoothMPC_dlubcc43 = smoothMPC_dl_cc + 699;
smoothMPC_FLOAT* smoothMPC_dsubcc43 = smoothMPC_ds_cc + 699;
smoothMPC_FLOAT* smoothMPC_ccrhsub43 = smoothMPC_ccrhs + 699;
smoothMPC_FLOAT smoothMPC_Phi43[11];
smoothMPC_FLOAT smoothMPC_W43[11];
smoothMPC_FLOAT smoothMPC_Ysd43[9];
smoothMPC_FLOAT smoothMPC_Lsd43[9];
smoothMPC_FLOAT* smoothMPC_z44 = smoothMPC_z + 484;
smoothMPC_FLOAT* smoothMPC_dzaff44 = smoothMPC_dz_aff + 484;
smoothMPC_FLOAT* smoothMPC_dzcc44 = smoothMPC_dz_cc + 484;
smoothMPC_FLOAT* smoothMPC_rd44 = smoothMPC_rd + 484;
smoothMPC_FLOAT smoothMPC_Lbyrd44[11];
smoothMPC_FLOAT* smoothMPC_grad_cost44 = smoothMPC_grad_cost + 484;
smoothMPC_FLOAT* smoothMPC_grad_eq44 = smoothMPC_grad_eq + 484;
smoothMPC_FLOAT* smoothMPC_grad_ineq44 = smoothMPC_grad_ineq + 484;
smoothMPC_FLOAT smoothMPC_ctv44[11];
smoothMPC_FLOAT* smoothMPC_v44 = smoothMPC_v + 132;
smoothMPC_FLOAT smoothMPC_re44[3];
smoothMPC_FLOAT smoothMPC_beta44[3];
smoothMPC_FLOAT smoothMPC_betacc44[3];
smoothMPC_FLOAT* smoothMPC_dvaff44 = smoothMPC_dv_aff + 132;
smoothMPC_FLOAT* smoothMPC_dvcc44 = smoothMPC_dv_cc + 132;
smoothMPC_FLOAT smoothMPC_V44[33];
smoothMPC_FLOAT smoothMPC_Yd44[6];
smoothMPC_FLOAT smoothMPC_Ld44[6];
smoothMPC_FLOAT smoothMPC_yy44[3];
smoothMPC_FLOAT smoothMPC_bmy44[3];
int smoothMPC_lbIdx44[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb44 = smoothMPC_l + 704;
smoothMPC_FLOAT* smoothMPC_slb44 = smoothMPC_s + 704;
smoothMPC_FLOAT* smoothMPC_llbbyslb44 = smoothMPC_lbys + 704;
smoothMPC_FLOAT smoothMPC_rilb44[11];
smoothMPC_FLOAT* smoothMPC_dllbaff44 = smoothMPC_dl_aff + 704;
smoothMPC_FLOAT* smoothMPC_dslbaff44 = smoothMPC_ds_aff + 704;
smoothMPC_FLOAT* smoothMPC_dllbcc44 = smoothMPC_dl_cc + 704;
smoothMPC_FLOAT* smoothMPC_dslbcc44 = smoothMPC_ds_cc + 704;
smoothMPC_FLOAT* smoothMPC_ccrhsl44 = smoothMPC_ccrhs + 704;
int smoothMPC_ubIdx44[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub44 = smoothMPC_l + 715;
smoothMPC_FLOAT* smoothMPC_sub44 = smoothMPC_s + 715;
smoothMPC_FLOAT* smoothMPC_lubbysub44 = smoothMPC_lbys + 715;
smoothMPC_FLOAT smoothMPC_riub44[5];
smoothMPC_FLOAT* smoothMPC_dlubaff44 = smoothMPC_dl_aff + 715;
smoothMPC_FLOAT* smoothMPC_dsubaff44 = smoothMPC_ds_aff + 715;
smoothMPC_FLOAT* smoothMPC_dlubcc44 = smoothMPC_dl_cc + 715;
smoothMPC_FLOAT* smoothMPC_dsubcc44 = smoothMPC_ds_cc + 715;
smoothMPC_FLOAT* smoothMPC_ccrhsub44 = smoothMPC_ccrhs + 715;
smoothMPC_FLOAT smoothMPC_Phi44[11];
smoothMPC_FLOAT smoothMPC_W44[11];
smoothMPC_FLOAT smoothMPC_Ysd44[9];
smoothMPC_FLOAT smoothMPC_Lsd44[9];
smoothMPC_FLOAT* smoothMPC_z45 = smoothMPC_z + 495;
smoothMPC_FLOAT* smoothMPC_dzaff45 = smoothMPC_dz_aff + 495;
smoothMPC_FLOAT* smoothMPC_dzcc45 = smoothMPC_dz_cc + 495;
smoothMPC_FLOAT* smoothMPC_rd45 = smoothMPC_rd + 495;
smoothMPC_FLOAT smoothMPC_Lbyrd45[11];
smoothMPC_FLOAT* smoothMPC_grad_cost45 = smoothMPC_grad_cost + 495;
smoothMPC_FLOAT* smoothMPC_grad_eq45 = smoothMPC_grad_eq + 495;
smoothMPC_FLOAT* smoothMPC_grad_ineq45 = smoothMPC_grad_ineq + 495;
smoothMPC_FLOAT smoothMPC_ctv45[11];
smoothMPC_FLOAT* smoothMPC_v45 = smoothMPC_v + 135;
smoothMPC_FLOAT smoothMPC_re45[3];
smoothMPC_FLOAT smoothMPC_beta45[3];
smoothMPC_FLOAT smoothMPC_betacc45[3];
smoothMPC_FLOAT* smoothMPC_dvaff45 = smoothMPC_dv_aff + 135;
smoothMPC_FLOAT* smoothMPC_dvcc45 = smoothMPC_dv_cc + 135;
smoothMPC_FLOAT smoothMPC_V45[33];
smoothMPC_FLOAT smoothMPC_Yd45[6];
smoothMPC_FLOAT smoothMPC_Ld45[6];
smoothMPC_FLOAT smoothMPC_yy45[3];
smoothMPC_FLOAT smoothMPC_bmy45[3];
int smoothMPC_lbIdx45[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb45 = smoothMPC_l + 720;
smoothMPC_FLOAT* smoothMPC_slb45 = smoothMPC_s + 720;
smoothMPC_FLOAT* smoothMPC_llbbyslb45 = smoothMPC_lbys + 720;
smoothMPC_FLOAT smoothMPC_rilb45[11];
smoothMPC_FLOAT* smoothMPC_dllbaff45 = smoothMPC_dl_aff + 720;
smoothMPC_FLOAT* smoothMPC_dslbaff45 = smoothMPC_ds_aff + 720;
smoothMPC_FLOAT* smoothMPC_dllbcc45 = smoothMPC_dl_cc + 720;
smoothMPC_FLOAT* smoothMPC_dslbcc45 = smoothMPC_ds_cc + 720;
smoothMPC_FLOAT* smoothMPC_ccrhsl45 = smoothMPC_ccrhs + 720;
int smoothMPC_ubIdx45[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub45 = smoothMPC_l + 731;
smoothMPC_FLOAT* smoothMPC_sub45 = smoothMPC_s + 731;
smoothMPC_FLOAT* smoothMPC_lubbysub45 = smoothMPC_lbys + 731;
smoothMPC_FLOAT smoothMPC_riub45[5];
smoothMPC_FLOAT* smoothMPC_dlubaff45 = smoothMPC_dl_aff + 731;
smoothMPC_FLOAT* smoothMPC_dsubaff45 = smoothMPC_ds_aff + 731;
smoothMPC_FLOAT* smoothMPC_dlubcc45 = smoothMPC_dl_cc + 731;
smoothMPC_FLOAT* smoothMPC_dsubcc45 = smoothMPC_ds_cc + 731;
smoothMPC_FLOAT* smoothMPC_ccrhsub45 = smoothMPC_ccrhs + 731;
smoothMPC_FLOAT smoothMPC_Phi45[11];
smoothMPC_FLOAT smoothMPC_W45[11];
smoothMPC_FLOAT smoothMPC_Ysd45[9];
smoothMPC_FLOAT smoothMPC_Lsd45[9];
smoothMPC_FLOAT* smoothMPC_z46 = smoothMPC_z + 506;
smoothMPC_FLOAT* smoothMPC_dzaff46 = smoothMPC_dz_aff + 506;
smoothMPC_FLOAT* smoothMPC_dzcc46 = smoothMPC_dz_cc + 506;
smoothMPC_FLOAT* smoothMPC_rd46 = smoothMPC_rd + 506;
smoothMPC_FLOAT smoothMPC_Lbyrd46[11];
smoothMPC_FLOAT* smoothMPC_grad_cost46 = smoothMPC_grad_cost + 506;
smoothMPC_FLOAT* smoothMPC_grad_eq46 = smoothMPC_grad_eq + 506;
smoothMPC_FLOAT* smoothMPC_grad_ineq46 = smoothMPC_grad_ineq + 506;
smoothMPC_FLOAT smoothMPC_ctv46[11];
smoothMPC_FLOAT* smoothMPC_v46 = smoothMPC_v + 138;
smoothMPC_FLOAT smoothMPC_re46[3];
smoothMPC_FLOAT smoothMPC_beta46[3];
smoothMPC_FLOAT smoothMPC_betacc46[3];
smoothMPC_FLOAT* smoothMPC_dvaff46 = smoothMPC_dv_aff + 138;
smoothMPC_FLOAT* smoothMPC_dvcc46 = smoothMPC_dv_cc + 138;
smoothMPC_FLOAT smoothMPC_V46[33];
smoothMPC_FLOAT smoothMPC_Yd46[6];
smoothMPC_FLOAT smoothMPC_Ld46[6];
smoothMPC_FLOAT smoothMPC_yy46[3];
smoothMPC_FLOAT smoothMPC_bmy46[3];
int smoothMPC_lbIdx46[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb46 = smoothMPC_l + 736;
smoothMPC_FLOAT* smoothMPC_slb46 = smoothMPC_s + 736;
smoothMPC_FLOAT* smoothMPC_llbbyslb46 = smoothMPC_lbys + 736;
smoothMPC_FLOAT smoothMPC_rilb46[11];
smoothMPC_FLOAT* smoothMPC_dllbaff46 = smoothMPC_dl_aff + 736;
smoothMPC_FLOAT* smoothMPC_dslbaff46 = smoothMPC_ds_aff + 736;
smoothMPC_FLOAT* smoothMPC_dllbcc46 = smoothMPC_dl_cc + 736;
smoothMPC_FLOAT* smoothMPC_dslbcc46 = smoothMPC_ds_cc + 736;
smoothMPC_FLOAT* smoothMPC_ccrhsl46 = smoothMPC_ccrhs + 736;
int smoothMPC_ubIdx46[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub46 = smoothMPC_l + 747;
smoothMPC_FLOAT* smoothMPC_sub46 = smoothMPC_s + 747;
smoothMPC_FLOAT* smoothMPC_lubbysub46 = smoothMPC_lbys + 747;
smoothMPC_FLOAT smoothMPC_riub46[5];
smoothMPC_FLOAT* smoothMPC_dlubaff46 = smoothMPC_dl_aff + 747;
smoothMPC_FLOAT* smoothMPC_dsubaff46 = smoothMPC_ds_aff + 747;
smoothMPC_FLOAT* smoothMPC_dlubcc46 = smoothMPC_dl_cc + 747;
smoothMPC_FLOAT* smoothMPC_dsubcc46 = smoothMPC_ds_cc + 747;
smoothMPC_FLOAT* smoothMPC_ccrhsub46 = smoothMPC_ccrhs + 747;
smoothMPC_FLOAT smoothMPC_Phi46[11];
smoothMPC_FLOAT smoothMPC_W46[11];
smoothMPC_FLOAT smoothMPC_Ysd46[9];
smoothMPC_FLOAT smoothMPC_Lsd46[9];
smoothMPC_FLOAT* smoothMPC_z47 = smoothMPC_z + 517;
smoothMPC_FLOAT* smoothMPC_dzaff47 = smoothMPC_dz_aff + 517;
smoothMPC_FLOAT* smoothMPC_dzcc47 = smoothMPC_dz_cc + 517;
smoothMPC_FLOAT* smoothMPC_rd47 = smoothMPC_rd + 517;
smoothMPC_FLOAT smoothMPC_Lbyrd47[11];
smoothMPC_FLOAT* smoothMPC_grad_cost47 = smoothMPC_grad_cost + 517;
smoothMPC_FLOAT* smoothMPC_grad_eq47 = smoothMPC_grad_eq + 517;
smoothMPC_FLOAT* smoothMPC_grad_ineq47 = smoothMPC_grad_ineq + 517;
smoothMPC_FLOAT smoothMPC_ctv47[11];
smoothMPC_FLOAT* smoothMPC_v47 = smoothMPC_v + 141;
smoothMPC_FLOAT smoothMPC_re47[3];
smoothMPC_FLOAT smoothMPC_beta47[3];
smoothMPC_FLOAT smoothMPC_betacc47[3];
smoothMPC_FLOAT* smoothMPC_dvaff47 = smoothMPC_dv_aff + 141;
smoothMPC_FLOAT* smoothMPC_dvcc47 = smoothMPC_dv_cc + 141;
smoothMPC_FLOAT smoothMPC_V47[33];
smoothMPC_FLOAT smoothMPC_Yd47[6];
smoothMPC_FLOAT smoothMPC_Ld47[6];
smoothMPC_FLOAT smoothMPC_yy47[3];
smoothMPC_FLOAT smoothMPC_bmy47[3];
int smoothMPC_lbIdx47[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb47 = smoothMPC_l + 752;
smoothMPC_FLOAT* smoothMPC_slb47 = smoothMPC_s + 752;
smoothMPC_FLOAT* smoothMPC_llbbyslb47 = smoothMPC_lbys + 752;
smoothMPC_FLOAT smoothMPC_rilb47[11];
smoothMPC_FLOAT* smoothMPC_dllbaff47 = smoothMPC_dl_aff + 752;
smoothMPC_FLOAT* smoothMPC_dslbaff47 = smoothMPC_ds_aff + 752;
smoothMPC_FLOAT* smoothMPC_dllbcc47 = smoothMPC_dl_cc + 752;
smoothMPC_FLOAT* smoothMPC_dslbcc47 = smoothMPC_ds_cc + 752;
smoothMPC_FLOAT* smoothMPC_ccrhsl47 = smoothMPC_ccrhs + 752;
int smoothMPC_ubIdx47[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub47 = smoothMPC_l + 763;
smoothMPC_FLOAT* smoothMPC_sub47 = smoothMPC_s + 763;
smoothMPC_FLOAT* smoothMPC_lubbysub47 = smoothMPC_lbys + 763;
smoothMPC_FLOAT smoothMPC_riub47[5];
smoothMPC_FLOAT* smoothMPC_dlubaff47 = smoothMPC_dl_aff + 763;
smoothMPC_FLOAT* smoothMPC_dsubaff47 = smoothMPC_ds_aff + 763;
smoothMPC_FLOAT* smoothMPC_dlubcc47 = smoothMPC_dl_cc + 763;
smoothMPC_FLOAT* smoothMPC_dsubcc47 = smoothMPC_ds_cc + 763;
smoothMPC_FLOAT* smoothMPC_ccrhsub47 = smoothMPC_ccrhs + 763;
smoothMPC_FLOAT smoothMPC_Phi47[11];
smoothMPC_FLOAT smoothMPC_W47[11];
smoothMPC_FLOAT smoothMPC_Ysd47[9];
smoothMPC_FLOAT smoothMPC_Lsd47[9];
smoothMPC_FLOAT* smoothMPC_z48 = smoothMPC_z + 528;
smoothMPC_FLOAT* smoothMPC_dzaff48 = smoothMPC_dz_aff + 528;
smoothMPC_FLOAT* smoothMPC_dzcc48 = smoothMPC_dz_cc + 528;
smoothMPC_FLOAT* smoothMPC_rd48 = smoothMPC_rd + 528;
smoothMPC_FLOAT smoothMPC_Lbyrd48[11];
smoothMPC_FLOAT* smoothMPC_grad_cost48 = smoothMPC_grad_cost + 528;
smoothMPC_FLOAT* smoothMPC_grad_eq48 = smoothMPC_grad_eq + 528;
smoothMPC_FLOAT* smoothMPC_grad_ineq48 = smoothMPC_grad_ineq + 528;
smoothMPC_FLOAT smoothMPC_ctv48[11];
smoothMPC_FLOAT* smoothMPC_v48 = smoothMPC_v + 144;
smoothMPC_FLOAT smoothMPC_re48[3];
smoothMPC_FLOAT smoothMPC_beta48[3];
smoothMPC_FLOAT smoothMPC_betacc48[3];
smoothMPC_FLOAT* smoothMPC_dvaff48 = smoothMPC_dv_aff + 144;
smoothMPC_FLOAT* smoothMPC_dvcc48 = smoothMPC_dv_cc + 144;
smoothMPC_FLOAT smoothMPC_V48[33];
smoothMPC_FLOAT smoothMPC_Yd48[6];
smoothMPC_FLOAT smoothMPC_Ld48[6];
smoothMPC_FLOAT smoothMPC_yy48[3];
smoothMPC_FLOAT smoothMPC_bmy48[3];
int smoothMPC_lbIdx48[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb48 = smoothMPC_l + 768;
smoothMPC_FLOAT* smoothMPC_slb48 = smoothMPC_s + 768;
smoothMPC_FLOAT* smoothMPC_llbbyslb48 = smoothMPC_lbys + 768;
smoothMPC_FLOAT smoothMPC_rilb48[11];
smoothMPC_FLOAT* smoothMPC_dllbaff48 = smoothMPC_dl_aff + 768;
smoothMPC_FLOAT* smoothMPC_dslbaff48 = smoothMPC_ds_aff + 768;
smoothMPC_FLOAT* smoothMPC_dllbcc48 = smoothMPC_dl_cc + 768;
smoothMPC_FLOAT* smoothMPC_dslbcc48 = smoothMPC_ds_cc + 768;
smoothMPC_FLOAT* smoothMPC_ccrhsl48 = smoothMPC_ccrhs + 768;
int smoothMPC_ubIdx48[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub48 = smoothMPC_l + 779;
smoothMPC_FLOAT* smoothMPC_sub48 = smoothMPC_s + 779;
smoothMPC_FLOAT* smoothMPC_lubbysub48 = smoothMPC_lbys + 779;
smoothMPC_FLOAT smoothMPC_riub48[5];
smoothMPC_FLOAT* smoothMPC_dlubaff48 = smoothMPC_dl_aff + 779;
smoothMPC_FLOAT* smoothMPC_dsubaff48 = smoothMPC_ds_aff + 779;
smoothMPC_FLOAT* smoothMPC_dlubcc48 = smoothMPC_dl_cc + 779;
smoothMPC_FLOAT* smoothMPC_dsubcc48 = smoothMPC_ds_cc + 779;
smoothMPC_FLOAT* smoothMPC_ccrhsub48 = smoothMPC_ccrhs + 779;
smoothMPC_FLOAT smoothMPC_Phi48[11];
smoothMPC_FLOAT smoothMPC_W48[11];
smoothMPC_FLOAT smoothMPC_Ysd48[9];
smoothMPC_FLOAT smoothMPC_Lsd48[9];
smoothMPC_FLOAT* smoothMPC_z49 = smoothMPC_z + 539;
smoothMPC_FLOAT* smoothMPC_dzaff49 = smoothMPC_dz_aff + 539;
smoothMPC_FLOAT* smoothMPC_dzcc49 = smoothMPC_dz_cc + 539;
smoothMPC_FLOAT* smoothMPC_rd49 = smoothMPC_rd + 539;
smoothMPC_FLOAT smoothMPC_Lbyrd49[11];
smoothMPC_FLOAT* smoothMPC_grad_cost49 = smoothMPC_grad_cost + 539;
smoothMPC_FLOAT* smoothMPC_grad_eq49 = smoothMPC_grad_eq + 539;
smoothMPC_FLOAT* smoothMPC_grad_ineq49 = smoothMPC_grad_ineq + 539;
smoothMPC_FLOAT smoothMPC_ctv49[11];
smoothMPC_FLOAT* smoothMPC_v49 = smoothMPC_v + 147;
smoothMPC_FLOAT smoothMPC_re49[3];
smoothMPC_FLOAT smoothMPC_beta49[3];
smoothMPC_FLOAT smoothMPC_betacc49[3];
smoothMPC_FLOAT* smoothMPC_dvaff49 = smoothMPC_dv_aff + 147;
smoothMPC_FLOAT* smoothMPC_dvcc49 = smoothMPC_dv_cc + 147;
smoothMPC_FLOAT smoothMPC_V49[33];
smoothMPC_FLOAT smoothMPC_Yd49[6];
smoothMPC_FLOAT smoothMPC_Ld49[6];
smoothMPC_FLOAT smoothMPC_yy49[3];
smoothMPC_FLOAT smoothMPC_bmy49[3];
int smoothMPC_lbIdx49[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb49 = smoothMPC_l + 784;
smoothMPC_FLOAT* smoothMPC_slb49 = smoothMPC_s + 784;
smoothMPC_FLOAT* smoothMPC_llbbyslb49 = smoothMPC_lbys + 784;
smoothMPC_FLOAT smoothMPC_rilb49[11];
smoothMPC_FLOAT* smoothMPC_dllbaff49 = smoothMPC_dl_aff + 784;
smoothMPC_FLOAT* smoothMPC_dslbaff49 = smoothMPC_ds_aff + 784;
smoothMPC_FLOAT* smoothMPC_dllbcc49 = smoothMPC_dl_cc + 784;
smoothMPC_FLOAT* smoothMPC_dslbcc49 = smoothMPC_ds_cc + 784;
smoothMPC_FLOAT* smoothMPC_ccrhsl49 = smoothMPC_ccrhs + 784;
int smoothMPC_ubIdx49[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub49 = smoothMPC_l + 795;
smoothMPC_FLOAT* smoothMPC_sub49 = smoothMPC_s + 795;
smoothMPC_FLOAT* smoothMPC_lubbysub49 = smoothMPC_lbys + 795;
smoothMPC_FLOAT smoothMPC_riub49[5];
smoothMPC_FLOAT* smoothMPC_dlubaff49 = smoothMPC_dl_aff + 795;
smoothMPC_FLOAT* smoothMPC_dsubaff49 = smoothMPC_ds_aff + 795;
smoothMPC_FLOAT* smoothMPC_dlubcc49 = smoothMPC_dl_cc + 795;
smoothMPC_FLOAT* smoothMPC_dsubcc49 = smoothMPC_ds_cc + 795;
smoothMPC_FLOAT* smoothMPC_ccrhsub49 = smoothMPC_ccrhs + 795;
smoothMPC_FLOAT smoothMPC_Phi49[11];
smoothMPC_FLOAT smoothMPC_W49[11];
smoothMPC_FLOAT smoothMPC_Ysd49[9];
smoothMPC_FLOAT smoothMPC_Lsd49[9];
smoothMPC_FLOAT* smoothMPC_z50 = smoothMPC_z + 550;
smoothMPC_FLOAT* smoothMPC_dzaff50 = smoothMPC_dz_aff + 550;
smoothMPC_FLOAT* smoothMPC_dzcc50 = smoothMPC_dz_cc + 550;
smoothMPC_FLOAT* smoothMPC_rd50 = smoothMPC_rd + 550;
smoothMPC_FLOAT smoothMPC_Lbyrd50[11];
smoothMPC_FLOAT* smoothMPC_grad_cost50 = smoothMPC_grad_cost + 550;
smoothMPC_FLOAT* smoothMPC_grad_eq50 = smoothMPC_grad_eq + 550;
smoothMPC_FLOAT* smoothMPC_grad_ineq50 = smoothMPC_grad_ineq + 550;
smoothMPC_FLOAT smoothMPC_ctv50[11];
smoothMPC_FLOAT* smoothMPC_v50 = smoothMPC_v + 150;
smoothMPC_FLOAT smoothMPC_re50[3];
smoothMPC_FLOAT smoothMPC_beta50[3];
smoothMPC_FLOAT smoothMPC_betacc50[3];
smoothMPC_FLOAT* smoothMPC_dvaff50 = smoothMPC_dv_aff + 150;
smoothMPC_FLOAT* smoothMPC_dvcc50 = smoothMPC_dv_cc + 150;
smoothMPC_FLOAT smoothMPC_V50[33];
smoothMPC_FLOAT smoothMPC_Yd50[6];
smoothMPC_FLOAT smoothMPC_Ld50[6];
smoothMPC_FLOAT smoothMPC_yy50[3];
smoothMPC_FLOAT smoothMPC_bmy50[3];
int smoothMPC_lbIdx50[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb50 = smoothMPC_l + 800;
smoothMPC_FLOAT* smoothMPC_slb50 = smoothMPC_s + 800;
smoothMPC_FLOAT* smoothMPC_llbbyslb50 = smoothMPC_lbys + 800;
smoothMPC_FLOAT smoothMPC_rilb50[11];
smoothMPC_FLOAT* smoothMPC_dllbaff50 = smoothMPC_dl_aff + 800;
smoothMPC_FLOAT* smoothMPC_dslbaff50 = smoothMPC_ds_aff + 800;
smoothMPC_FLOAT* smoothMPC_dllbcc50 = smoothMPC_dl_cc + 800;
smoothMPC_FLOAT* smoothMPC_dslbcc50 = smoothMPC_ds_cc + 800;
smoothMPC_FLOAT* smoothMPC_ccrhsl50 = smoothMPC_ccrhs + 800;
int smoothMPC_ubIdx50[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub50 = smoothMPC_l + 811;
smoothMPC_FLOAT* smoothMPC_sub50 = smoothMPC_s + 811;
smoothMPC_FLOAT* smoothMPC_lubbysub50 = smoothMPC_lbys + 811;
smoothMPC_FLOAT smoothMPC_riub50[5];
smoothMPC_FLOAT* smoothMPC_dlubaff50 = smoothMPC_dl_aff + 811;
smoothMPC_FLOAT* smoothMPC_dsubaff50 = smoothMPC_ds_aff + 811;
smoothMPC_FLOAT* smoothMPC_dlubcc50 = smoothMPC_dl_cc + 811;
smoothMPC_FLOAT* smoothMPC_dsubcc50 = smoothMPC_ds_cc + 811;
smoothMPC_FLOAT* smoothMPC_ccrhsub50 = smoothMPC_ccrhs + 811;
smoothMPC_FLOAT smoothMPC_Phi50[11];
smoothMPC_FLOAT smoothMPC_W50[11];
smoothMPC_FLOAT smoothMPC_Ysd50[9];
smoothMPC_FLOAT smoothMPC_Lsd50[9];
smoothMPC_FLOAT* smoothMPC_z51 = smoothMPC_z + 561;
smoothMPC_FLOAT* smoothMPC_dzaff51 = smoothMPC_dz_aff + 561;
smoothMPC_FLOAT* smoothMPC_dzcc51 = smoothMPC_dz_cc + 561;
smoothMPC_FLOAT* smoothMPC_rd51 = smoothMPC_rd + 561;
smoothMPC_FLOAT smoothMPC_Lbyrd51[11];
smoothMPC_FLOAT* smoothMPC_grad_cost51 = smoothMPC_grad_cost + 561;
smoothMPC_FLOAT* smoothMPC_grad_eq51 = smoothMPC_grad_eq + 561;
smoothMPC_FLOAT* smoothMPC_grad_ineq51 = smoothMPC_grad_ineq + 561;
smoothMPC_FLOAT smoothMPC_ctv51[11];
smoothMPC_FLOAT* smoothMPC_v51 = smoothMPC_v + 153;
smoothMPC_FLOAT smoothMPC_re51[3];
smoothMPC_FLOAT smoothMPC_beta51[3];
smoothMPC_FLOAT smoothMPC_betacc51[3];
smoothMPC_FLOAT* smoothMPC_dvaff51 = smoothMPC_dv_aff + 153;
smoothMPC_FLOAT* smoothMPC_dvcc51 = smoothMPC_dv_cc + 153;
smoothMPC_FLOAT smoothMPC_V51[33];
smoothMPC_FLOAT smoothMPC_Yd51[6];
smoothMPC_FLOAT smoothMPC_Ld51[6];
smoothMPC_FLOAT smoothMPC_yy51[3];
smoothMPC_FLOAT smoothMPC_bmy51[3];
int smoothMPC_lbIdx51[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb51 = smoothMPC_l + 816;
smoothMPC_FLOAT* smoothMPC_slb51 = smoothMPC_s + 816;
smoothMPC_FLOAT* smoothMPC_llbbyslb51 = smoothMPC_lbys + 816;
smoothMPC_FLOAT smoothMPC_rilb51[11];
smoothMPC_FLOAT* smoothMPC_dllbaff51 = smoothMPC_dl_aff + 816;
smoothMPC_FLOAT* smoothMPC_dslbaff51 = smoothMPC_ds_aff + 816;
smoothMPC_FLOAT* smoothMPC_dllbcc51 = smoothMPC_dl_cc + 816;
smoothMPC_FLOAT* smoothMPC_dslbcc51 = smoothMPC_ds_cc + 816;
smoothMPC_FLOAT* smoothMPC_ccrhsl51 = smoothMPC_ccrhs + 816;
int smoothMPC_ubIdx51[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub51 = smoothMPC_l + 827;
smoothMPC_FLOAT* smoothMPC_sub51 = smoothMPC_s + 827;
smoothMPC_FLOAT* smoothMPC_lubbysub51 = smoothMPC_lbys + 827;
smoothMPC_FLOAT smoothMPC_riub51[5];
smoothMPC_FLOAT* smoothMPC_dlubaff51 = smoothMPC_dl_aff + 827;
smoothMPC_FLOAT* smoothMPC_dsubaff51 = smoothMPC_ds_aff + 827;
smoothMPC_FLOAT* smoothMPC_dlubcc51 = smoothMPC_dl_cc + 827;
smoothMPC_FLOAT* smoothMPC_dsubcc51 = smoothMPC_ds_cc + 827;
smoothMPC_FLOAT* smoothMPC_ccrhsub51 = smoothMPC_ccrhs + 827;
smoothMPC_FLOAT smoothMPC_Phi51[11];
smoothMPC_FLOAT smoothMPC_W51[11];
smoothMPC_FLOAT smoothMPC_Ysd51[9];
smoothMPC_FLOAT smoothMPC_Lsd51[9];
smoothMPC_FLOAT* smoothMPC_z52 = smoothMPC_z + 572;
smoothMPC_FLOAT* smoothMPC_dzaff52 = smoothMPC_dz_aff + 572;
smoothMPC_FLOAT* smoothMPC_dzcc52 = smoothMPC_dz_cc + 572;
smoothMPC_FLOAT* smoothMPC_rd52 = smoothMPC_rd + 572;
smoothMPC_FLOAT smoothMPC_Lbyrd52[11];
smoothMPC_FLOAT* smoothMPC_grad_cost52 = smoothMPC_grad_cost + 572;
smoothMPC_FLOAT* smoothMPC_grad_eq52 = smoothMPC_grad_eq + 572;
smoothMPC_FLOAT* smoothMPC_grad_ineq52 = smoothMPC_grad_ineq + 572;
smoothMPC_FLOAT smoothMPC_ctv52[11];
smoothMPC_FLOAT* smoothMPC_v52 = smoothMPC_v + 156;
smoothMPC_FLOAT smoothMPC_re52[3];
smoothMPC_FLOAT smoothMPC_beta52[3];
smoothMPC_FLOAT smoothMPC_betacc52[3];
smoothMPC_FLOAT* smoothMPC_dvaff52 = smoothMPC_dv_aff + 156;
smoothMPC_FLOAT* smoothMPC_dvcc52 = smoothMPC_dv_cc + 156;
smoothMPC_FLOAT smoothMPC_V52[33];
smoothMPC_FLOAT smoothMPC_Yd52[6];
smoothMPC_FLOAT smoothMPC_Ld52[6];
smoothMPC_FLOAT smoothMPC_yy52[3];
smoothMPC_FLOAT smoothMPC_bmy52[3];
int smoothMPC_lbIdx52[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb52 = smoothMPC_l + 832;
smoothMPC_FLOAT* smoothMPC_slb52 = smoothMPC_s + 832;
smoothMPC_FLOAT* smoothMPC_llbbyslb52 = smoothMPC_lbys + 832;
smoothMPC_FLOAT smoothMPC_rilb52[11];
smoothMPC_FLOAT* smoothMPC_dllbaff52 = smoothMPC_dl_aff + 832;
smoothMPC_FLOAT* smoothMPC_dslbaff52 = smoothMPC_ds_aff + 832;
smoothMPC_FLOAT* smoothMPC_dllbcc52 = smoothMPC_dl_cc + 832;
smoothMPC_FLOAT* smoothMPC_dslbcc52 = smoothMPC_ds_cc + 832;
smoothMPC_FLOAT* smoothMPC_ccrhsl52 = smoothMPC_ccrhs + 832;
int smoothMPC_ubIdx52[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub52 = smoothMPC_l + 843;
smoothMPC_FLOAT* smoothMPC_sub52 = smoothMPC_s + 843;
smoothMPC_FLOAT* smoothMPC_lubbysub52 = smoothMPC_lbys + 843;
smoothMPC_FLOAT smoothMPC_riub52[5];
smoothMPC_FLOAT* smoothMPC_dlubaff52 = smoothMPC_dl_aff + 843;
smoothMPC_FLOAT* smoothMPC_dsubaff52 = smoothMPC_ds_aff + 843;
smoothMPC_FLOAT* smoothMPC_dlubcc52 = smoothMPC_dl_cc + 843;
smoothMPC_FLOAT* smoothMPC_dsubcc52 = smoothMPC_ds_cc + 843;
smoothMPC_FLOAT* smoothMPC_ccrhsub52 = smoothMPC_ccrhs + 843;
smoothMPC_FLOAT smoothMPC_Phi52[11];
smoothMPC_FLOAT smoothMPC_W52[11];
smoothMPC_FLOAT smoothMPC_Ysd52[9];
smoothMPC_FLOAT smoothMPC_Lsd52[9];
smoothMPC_FLOAT* smoothMPC_z53 = smoothMPC_z + 583;
smoothMPC_FLOAT* smoothMPC_dzaff53 = smoothMPC_dz_aff + 583;
smoothMPC_FLOAT* smoothMPC_dzcc53 = smoothMPC_dz_cc + 583;
smoothMPC_FLOAT* smoothMPC_rd53 = smoothMPC_rd + 583;
smoothMPC_FLOAT smoothMPC_Lbyrd53[11];
smoothMPC_FLOAT* smoothMPC_grad_cost53 = smoothMPC_grad_cost + 583;
smoothMPC_FLOAT* smoothMPC_grad_eq53 = smoothMPC_grad_eq + 583;
smoothMPC_FLOAT* smoothMPC_grad_ineq53 = smoothMPC_grad_ineq + 583;
smoothMPC_FLOAT smoothMPC_ctv53[11];
smoothMPC_FLOAT* smoothMPC_v53 = smoothMPC_v + 159;
smoothMPC_FLOAT smoothMPC_re53[3];
smoothMPC_FLOAT smoothMPC_beta53[3];
smoothMPC_FLOAT smoothMPC_betacc53[3];
smoothMPC_FLOAT* smoothMPC_dvaff53 = smoothMPC_dv_aff + 159;
smoothMPC_FLOAT* smoothMPC_dvcc53 = smoothMPC_dv_cc + 159;
smoothMPC_FLOAT smoothMPC_V53[33];
smoothMPC_FLOAT smoothMPC_Yd53[6];
smoothMPC_FLOAT smoothMPC_Ld53[6];
smoothMPC_FLOAT smoothMPC_yy53[3];
smoothMPC_FLOAT smoothMPC_bmy53[3];
int smoothMPC_lbIdx53[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb53 = smoothMPC_l + 848;
smoothMPC_FLOAT* smoothMPC_slb53 = smoothMPC_s + 848;
smoothMPC_FLOAT* smoothMPC_llbbyslb53 = smoothMPC_lbys + 848;
smoothMPC_FLOAT smoothMPC_rilb53[11];
smoothMPC_FLOAT* smoothMPC_dllbaff53 = smoothMPC_dl_aff + 848;
smoothMPC_FLOAT* smoothMPC_dslbaff53 = smoothMPC_ds_aff + 848;
smoothMPC_FLOAT* smoothMPC_dllbcc53 = smoothMPC_dl_cc + 848;
smoothMPC_FLOAT* smoothMPC_dslbcc53 = smoothMPC_ds_cc + 848;
smoothMPC_FLOAT* smoothMPC_ccrhsl53 = smoothMPC_ccrhs + 848;
int smoothMPC_ubIdx53[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub53 = smoothMPC_l + 859;
smoothMPC_FLOAT* smoothMPC_sub53 = smoothMPC_s + 859;
smoothMPC_FLOAT* smoothMPC_lubbysub53 = smoothMPC_lbys + 859;
smoothMPC_FLOAT smoothMPC_riub53[5];
smoothMPC_FLOAT* smoothMPC_dlubaff53 = smoothMPC_dl_aff + 859;
smoothMPC_FLOAT* smoothMPC_dsubaff53 = smoothMPC_ds_aff + 859;
smoothMPC_FLOAT* smoothMPC_dlubcc53 = smoothMPC_dl_cc + 859;
smoothMPC_FLOAT* smoothMPC_dsubcc53 = smoothMPC_ds_cc + 859;
smoothMPC_FLOAT* smoothMPC_ccrhsub53 = smoothMPC_ccrhs + 859;
smoothMPC_FLOAT smoothMPC_Phi53[11];
smoothMPC_FLOAT smoothMPC_W53[11];
smoothMPC_FLOAT smoothMPC_Ysd53[9];
smoothMPC_FLOAT smoothMPC_Lsd53[9];
smoothMPC_FLOAT* smoothMPC_z54 = smoothMPC_z + 594;
smoothMPC_FLOAT* smoothMPC_dzaff54 = smoothMPC_dz_aff + 594;
smoothMPC_FLOAT* smoothMPC_dzcc54 = smoothMPC_dz_cc + 594;
smoothMPC_FLOAT* smoothMPC_rd54 = smoothMPC_rd + 594;
smoothMPC_FLOAT smoothMPC_Lbyrd54[11];
smoothMPC_FLOAT* smoothMPC_grad_cost54 = smoothMPC_grad_cost + 594;
smoothMPC_FLOAT* smoothMPC_grad_eq54 = smoothMPC_grad_eq + 594;
smoothMPC_FLOAT* smoothMPC_grad_ineq54 = smoothMPC_grad_ineq + 594;
smoothMPC_FLOAT smoothMPC_ctv54[11];
smoothMPC_FLOAT* smoothMPC_v54 = smoothMPC_v + 162;
smoothMPC_FLOAT smoothMPC_re54[3];
smoothMPC_FLOAT smoothMPC_beta54[3];
smoothMPC_FLOAT smoothMPC_betacc54[3];
smoothMPC_FLOAT* smoothMPC_dvaff54 = smoothMPC_dv_aff + 162;
smoothMPC_FLOAT* smoothMPC_dvcc54 = smoothMPC_dv_cc + 162;
smoothMPC_FLOAT smoothMPC_V54[33];
smoothMPC_FLOAT smoothMPC_Yd54[6];
smoothMPC_FLOAT smoothMPC_Ld54[6];
smoothMPC_FLOAT smoothMPC_yy54[3];
smoothMPC_FLOAT smoothMPC_bmy54[3];
int smoothMPC_lbIdx54[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb54 = smoothMPC_l + 864;
smoothMPC_FLOAT* smoothMPC_slb54 = smoothMPC_s + 864;
smoothMPC_FLOAT* smoothMPC_llbbyslb54 = smoothMPC_lbys + 864;
smoothMPC_FLOAT smoothMPC_rilb54[11];
smoothMPC_FLOAT* smoothMPC_dllbaff54 = smoothMPC_dl_aff + 864;
smoothMPC_FLOAT* smoothMPC_dslbaff54 = smoothMPC_ds_aff + 864;
smoothMPC_FLOAT* smoothMPC_dllbcc54 = smoothMPC_dl_cc + 864;
smoothMPC_FLOAT* smoothMPC_dslbcc54 = smoothMPC_ds_cc + 864;
smoothMPC_FLOAT* smoothMPC_ccrhsl54 = smoothMPC_ccrhs + 864;
int smoothMPC_ubIdx54[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub54 = smoothMPC_l + 875;
smoothMPC_FLOAT* smoothMPC_sub54 = smoothMPC_s + 875;
smoothMPC_FLOAT* smoothMPC_lubbysub54 = smoothMPC_lbys + 875;
smoothMPC_FLOAT smoothMPC_riub54[5];
smoothMPC_FLOAT* smoothMPC_dlubaff54 = smoothMPC_dl_aff + 875;
smoothMPC_FLOAT* smoothMPC_dsubaff54 = smoothMPC_ds_aff + 875;
smoothMPC_FLOAT* smoothMPC_dlubcc54 = smoothMPC_dl_cc + 875;
smoothMPC_FLOAT* smoothMPC_dsubcc54 = smoothMPC_ds_cc + 875;
smoothMPC_FLOAT* smoothMPC_ccrhsub54 = smoothMPC_ccrhs + 875;
smoothMPC_FLOAT smoothMPC_Phi54[11];
smoothMPC_FLOAT smoothMPC_W54[11];
smoothMPC_FLOAT smoothMPC_Ysd54[9];
smoothMPC_FLOAT smoothMPC_Lsd54[9];
smoothMPC_FLOAT* smoothMPC_z55 = smoothMPC_z + 605;
smoothMPC_FLOAT* smoothMPC_dzaff55 = smoothMPC_dz_aff + 605;
smoothMPC_FLOAT* smoothMPC_dzcc55 = smoothMPC_dz_cc + 605;
smoothMPC_FLOAT* smoothMPC_rd55 = smoothMPC_rd + 605;
smoothMPC_FLOAT smoothMPC_Lbyrd55[11];
smoothMPC_FLOAT* smoothMPC_grad_cost55 = smoothMPC_grad_cost + 605;
smoothMPC_FLOAT* smoothMPC_grad_eq55 = smoothMPC_grad_eq + 605;
smoothMPC_FLOAT* smoothMPC_grad_ineq55 = smoothMPC_grad_ineq + 605;
smoothMPC_FLOAT smoothMPC_ctv55[11];
smoothMPC_FLOAT* smoothMPC_v55 = smoothMPC_v + 165;
smoothMPC_FLOAT smoothMPC_re55[3];
smoothMPC_FLOAT smoothMPC_beta55[3];
smoothMPC_FLOAT smoothMPC_betacc55[3];
smoothMPC_FLOAT* smoothMPC_dvaff55 = smoothMPC_dv_aff + 165;
smoothMPC_FLOAT* smoothMPC_dvcc55 = smoothMPC_dv_cc + 165;
smoothMPC_FLOAT smoothMPC_V55[33];
smoothMPC_FLOAT smoothMPC_Yd55[6];
smoothMPC_FLOAT smoothMPC_Ld55[6];
smoothMPC_FLOAT smoothMPC_yy55[3];
smoothMPC_FLOAT smoothMPC_bmy55[3];
int smoothMPC_lbIdx55[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb55 = smoothMPC_l + 880;
smoothMPC_FLOAT* smoothMPC_slb55 = smoothMPC_s + 880;
smoothMPC_FLOAT* smoothMPC_llbbyslb55 = smoothMPC_lbys + 880;
smoothMPC_FLOAT smoothMPC_rilb55[11];
smoothMPC_FLOAT* smoothMPC_dllbaff55 = smoothMPC_dl_aff + 880;
smoothMPC_FLOAT* smoothMPC_dslbaff55 = smoothMPC_ds_aff + 880;
smoothMPC_FLOAT* smoothMPC_dllbcc55 = smoothMPC_dl_cc + 880;
smoothMPC_FLOAT* smoothMPC_dslbcc55 = smoothMPC_ds_cc + 880;
smoothMPC_FLOAT* smoothMPC_ccrhsl55 = smoothMPC_ccrhs + 880;
int smoothMPC_ubIdx55[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub55 = smoothMPC_l + 891;
smoothMPC_FLOAT* smoothMPC_sub55 = smoothMPC_s + 891;
smoothMPC_FLOAT* smoothMPC_lubbysub55 = smoothMPC_lbys + 891;
smoothMPC_FLOAT smoothMPC_riub55[5];
smoothMPC_FLOAT* smoothMPC_dlubaff55 = smoothMPC_dl_aff + 891;
smoothMPC_FLOAT* smoothMPC_dsubaff55 = smoothMPC_ds_aff + 891;
smoothMPC_FLOAT* smoothMPC_dlubcc55 = smoothMPC_dl_cc + 891;
smoothMPC_FLOAT* smoothMPC_dsubcc55 = smoothMPC_ds_cc + 891;
smoothMPC_FLOAT* smoothMPC_ccrhsub55 = smoothMPC_ccrhs + 891;
smoothMPC_FLOAT smoothMPC_Phi55[11];
smoothMPC_FLOAT smoothMPC_W55[11];
smoothMPC_FLOAT smoothMPC_Ysd55[9];
smoothMPC_FLOAT smoothMPC_Lsd55[9];
smoothMPC_FLOAT* smoothMPC_z56 = smoothMPC_z + 616;
smoothMPC_FLOAT* smoothMPC_dzaff56 = smoothMPC_dz_aff + 616;
smoothMPC_FLOAT* smoothMPC_dzcc56 = smoothMPC_dz_cc + 616;
smoothMPC_FLOAT* smoothMPC_rd56 = smoothMPC_rd + 616;
smoothMPC_FLOAT smoothMPC_Lbyrd56[11];
smoothMPC_FLOAT* smoothMPC_grad_cost56 = smoothMPC_grad_cost + 616;
smoothMPC_FLOAT* smoothMPC_grad_eq56 = smoothMPC_grad_eq + 616;
smoothMPC_FLOAT* smoothMPC_grad_ineq56 = smoothMPC_grad_ineq + 616;
smoothMPC_FLOAT smoothMPC_ctv56[11];
smoothMPC_FLOAT* smoothMPC_v56 = smoothMPC_v + 168;
smoothMPC_FLOAT smoothMPC_re56[3];
smoothMPC_FLOAT smoothMPC_beta56[3];
smoothMPC_FLOAT smoothMPC_betacc56[3];
smoothMPC_FLOAT* smoothMPC_dvaff56 = smoothMPC_dv_aff + 168;
smoothMPC_FLOAT* smoothMPC_dvcc56 = smoothMPC_dv_cc + 168;
smoothMPC_FLOAT smoothMPC_V56[33];
smoothMPC_FLOAT smoothMPC_Yd56[6];
smoothMPC_FLOAT smoothMPC_Ld56[6];
smoothMPC_FLOAT smoothMPC_yy56[3];
smoothMPC_FLOAT smoothMPC_bmy56[3];
int smoothMPC_lbIdx56[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb56 = smoothMPC_l + 896;
smoothMPC_FLOAT* smoothMPC_slb56 = smoothMPC_s + 896;
smoothMPC_FLOAT* smoothMPC_llbbyslb56 = smoothMPC_lbys + 896;
smoothMPC_FLOAT smoothMPC_rilb56[11];
smoothMPC_FLOAT* smoothMPC_dllbaff56 = smoothMPC_dl_aff + 896;
smoothMPC_FLOAT* smoothMPC_dslbaff56 = smoothMPC_ds_aff + 896;
smoothMPC_FLOAT* smoothMPC_dllbcc56 = smoothMPC_dl_cc + 896;
smoothMPC_FLOAT* smoothMPC_dslbcc56 = smoothMPC_ds_cc + 896;
smoothMPC_FLOAT* smoothMPC_ccrhsl56 = smoothMPC_ccrhs + 896;
int smoothMPC_ubIdx56[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub56 = smoothMPC_l + 907;
smoothMPC_FLOAT* smoothMPC_sub56 = smoothMPC_s + 907;
smoothMPC_FLOAT* smoothMPC_lubbysub56 = smoothMPC_lbys + 907;
smoothMPC_FLOAT smoothMPC_riub56[5];
smoothMPC_FLOAT* smoothMPC_dlubaff56 = smoothMPC_dl_aff + 907;
smoothMPC_FLOAT* smoothMPC_dsubaff56 = smoothMPC_ds_aff + 907;
smoothMPC_FLOAT* smoothMPC_dlubcc56 = smoothMPC_dl_cc + 907;
smoothMPC_FLOAT* smoothMPC_dsubcc56 = smoothMPC_ds_cc + 907;
smoothMPC_FLOAT* smoothMPC_ccrhsub56 = smoothMPC_ccrhs + 907;
smoothMPC_FLOAT smoothMPC_Phi56[11];
smoothMPC_FLOAT smoothMPC_W56[11];
smoothMPC_FLOAT smoothMPC_Ysd56[9];
smoothMPC_FLOAT smoothMPC_Lsd56[9];
smoothMPC_FLOAT* smoothMPC_z57 = smoothMPC_z + 627;
smoothMPC_FLOAT* smoothMPC_dzaff57 = smoothMPC_dz_aff + 627;
smoothMPC_FLOAT* smoothMPC_dzcc57 = smoothMPC_dz_cc + 627;
smoothMPC_FLOAT* smoothMPC_rd57 = smoothMPC_rd + 627;
smoothMPC_FLOAT smoothMPC_Lbyrd57[11];
smoothMPC_FLOAT* smoothMPC_grad_cost57 = smoothMPC_grad_cost + 627;
smoothMPC_FLOAT* smoothMPC_grad_eq57 = smoothMPC_grad_eq + 627;
smoothMPC_FLOAT* smoothMPC_grad_ineq57 = smoothMPC_grad_ineq + 627;
smoothMPC_FLOAT smoothMPC_ctv57[11];
smoothMPC_FLOAT* smoothMPC_v57 = smoothMPC_v + 171;
smoothMPC_FLOAT smoothMPC_re57[3];
smoothMPC_FLOAT smoothMPC_beta57[3];
smoothMPC_FLOAT smoothMPC_betacc57[3];
smoothMPC_FLOAT* smoothMPC_dvaff57 = smoothMPC_dv_aff + 171;
smoothMPC_FLOAT* smoothMPC_dvcc57 = smoothMPC_dv_cc + 171;
smoothMPC_FLOAT smoothMPC_V57[33];
smoothMPC_FLOAT smoothMPC_Yd57[6];
smoothMPC_FLOAT smoothMPC_Ld57[6];
smoothMPC_FLOAT smoothMPC_yy57[3];
smoothMPC_FLOAT smoothMPC_bmy57[3];
int smoothMPC_lbIdx57[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb57 = smoothMPC_l + 912;
smoothMPC_FLOAT* smoothMPC_slb57 = smoothMPC_s + 912;
smoothMPC_FLOAT* smoothMPC_llbbyslb57 = smoothMPC_lbys + 912;
smoothMPC_FLOAT smoothMPC_rilb57[11];
smoothMPC_FLOAT* smoothMPC_dllbaff57 = smoothMPC_dl_aff + 912;
smoothMPC_FLOAT* smoothMPC_dslbaff57 = smoothMPC_ds_aff + 912;
smoothMPC_FLOAT* smoothMPC_dllbcc57 = smoothMPC_dl_cc + 912;
smoothMPC_FLOAT* smoothMPC_dslbcc57 = smoothMPC_ds_cc + 912;
smoothMPC_FLOAT* smoothMPC_ccrhsl57 = smoothMPC_ccrhs + 912;
int smoothMPC_ubIdx57[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub57 = smoothMPC_l + 923;
smoothMPC_FLOAT* smoothMPC_sub57 = smoothMPC_s + 923;
smoothMPC_FLOAT* smoothMPC_lubbysub57 = smoothMPC_lbys + 923;
smoothMPC_FLOAT smoothMPC_riub57[5];
smoothMPC_FLOAT* smoothMPC_dlubaff57 = smoothMPC_dl_aff + 923;
smoothMPC_FLOAT* smoothMPC_dsubaff57 = smoothMPC_ds_aff + 923;
smoothMPC_FLOAT* smoothMPC_dlubcc57 = smoothMPC_dl_cc + 923;
smoothMPC_FLOAT* smoothMPC_dsubcc57 = smoothMPC_ds_cc + 923;
smoothMPC_FLOAT* smoothMPC_ccrhsub57 = smoothMPC_ccrhs + 923;
smoothMPC_FLOAT smoothMPC_Phi57[11];
smoothMPC_FLOAT smoothMPC_W57[11];
smoothMPC_FLOAT smoothMPC_Ysd57[9];
smoothMPC_FLOAT smoothMPC_Lsd57[9];
smoothMPC_FLOAT* smoothMPC_z58 = smoothMPC_z + 638;
smoothMPC_FLOAT* smoothMPC_dzaff58 = smoothMPC_dz_aff + 638;
smoothMPC_FLOAT* smoothMPC_dzcc58 = smoothMPC_dz_cc + 638;
smoothMPC_FLOAT* smoothMPC_rd58 = smoothMPC_rd + 638;
smoothMPC_FLOAT smoothMPC_Lbyrd58[11];
smoothMPC_FLOAT* smoothMPC_grad_cost58 = smoothMPC_grad_cost + 638;
smoothMPC_FLOAT* smoothMPC_grad_eq58 = smoothMPC_grad_eq + 638;
smoothMPC_FLOAT* smoothMPC_grad_ineq58 = smoothMPC_grad_ineq + 638;
smoothMPC_FLOAT smoothMPC_ctv58[11];
smoothMPC_FLOAT* smoothMPC_v58 = smoothMPC_v + 174;
smoothMPC_FLOAT smoothMPC_re58[3];
smoothMPC_FLOAT smoothMPC_beta58[3];
smoothMPC_FLOAT smoothMPC_betacc58[3];
smoothMPC_FLOAT* smoothMPC_dvaff58 = smoothMPC_dv_aff + 174;
smoothMPC_FLOAT* smoothMPC_dvcc58 = smoothMPC_dv_cc + 174;
smoothMPC_FLOAT smoothMPC_V58[33];
smoothMPC_FLOAT smoothMPC_Yd58[6];
smoothMPC_FLOAT smoothMPC_Ld58[6];
smoothMPC_FLOAT smoothMPC_yy58[3];
smoothMPC_FLOAT smoothMPC_bmy58[3];
int smoothMPC_lbIdx58[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb58 = smoothMPC_l + 928;
smoothMPC_FLOAT* smoothMPC_slb58 = smoothMPC_s + 928;
smoothMPC_FLOAT* smoothMPC_llbbyslb58 = smoothMPC_lbys + 928;
smoothMPC_FLOAT smoothMPC_rilb58[11];
smoothMPC_FLOAT* smoothMPC_dllbaff58 = smoothMPC_dl_aff + 928;
smoothMPC_FLOAT* smoothMPC_dslbaff58 = smoothMPC_ds_aff + 928;
smoothMPC_FLOAT* smoothMPC_dllbcc58 = smoothMPC_dl_cc + 928;
smoothMPC_FLOAT* smoothMPC_dslbcc58 = smoothMPC_ds_cc + 928;
smoothMPC_FLOAT* smoothMPC_ccrhsl58 = smoothMPC_ccrhs + 928;
int smoothMPC_ubIdx58[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub58 = smoothMPC_l + 939;
smoothMPC_FLOAT* smoothMPC_sub58 = smoothMPC_s + 939;
smoothMPC_FLOAT* smoothMPC_lubbysub58 = smoothMPC_lbys + 939;
smoothMPC_FLOAT smoothMPC_riub58[5];
smoothMPC_FLOAT* smoothMPC_dlubaff58 = smoothMPC_dl_aff + 939;
smoothMPC_FLOAT* smoothMPC_dsubaff58 = smoothMPC_ds_aff + 939;
smoothMPC_FLOAT* smoothMPC_dlubcc58 = smoothMPC_dl_cc + 939;
smoothMPC_FLOAT* smoothMPC_dsubcc58 = smoothMPC_ds_cc + 939;
smoothMPC_FLOAT* smoothMPC_ccrhsub58 = smoothMPC_ccrhs + 939;
smoothMPC_FLOAT smoothMPC_Phi58[11];
smoothMPC_FLOAT smoothMPC_W58[11];
smoothMPC_FLOAT smoothMPC_Ysd58[9];
smoothMPC_FLOAT smoothMPC_Lsd58[9];
smoothMPC_FLOAT* smoothMPC_z59 = smoothMPC_z + 649;
smoothMPC_FLOAT* smoothMPC_dzaff59 = smoothMPC_dz_aff + 649;
smoothMPC_FLOAT* smoothMPC_dzcc59 = smoothMPC_dz_cc + 649;
smoothMPC_FLOAT* smoothMPC_rd59 = smoothMPC_rd + 649;
smoothMPC_FLOAT smoothMPC_Lbyrd59[11];
smoothMPC_FLOAT* smoothMPC_grad_cost59 = smoothMPC_grad_cost + 649;
smoothMPC_FLOAT* smoothMPC_grad_eq59 = smoothMPC_grad_eq + 649;
smoothMPC_FLOAT* smoothMPC_grad_ineq59 = smoothMPC_grad_ineq + 649;
smoothMPC_FLOAT smoothMPC_ctv59[11];
smoothMPC_FLOAT* smoothMPC_v59 = smoothMPC_v + 177;
smoothMPC_FLOAT smoothMPC_re59[3];
smoothMPC_FLOAT smoothMPC_beta59[3];
smoothMPC_FLOAT smoothMPC_betacc59[3];
smoothMPC_FLOAT* smoothMPC_dvaff59 = smoothMPC_dv_aff + 177;
smoothMPC_FLOAT* smoothMPC_dvcc59 = smoothMPC_dv_cc + 177;
smoothMPC_FLOAT smoothMPC_V59[33];
smoothMPC_FLOAT smoothMPC_Yd59[6];
smoothMPC_FLOAT smoothMPC_Ld59[6];
smoothMPC_FLOAT smoothMPC_yy59[3];
smoothMPC_FLOAT smoothMPC_bmy59[3];
int smoothMPC_lbIdx59[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb59 = smoothMPC_l + 944;
smoothMPC_FLOAT* smoothMPC_slb59 = smoothMPC_s + 944;
smoothMPC_FLOAT* smoothMPC_llbbyslb59 = smoothMPC_lbys + 944;
smoothMPC_FLOAT smoothMPC_rilb59[11];
smoothMPC_FLOAT* smoothMPC_dllbaff59 = smoothMPC_dl_aff + 944;
smoothMPC_FLOAT* smoothMPC_dslbaff59 = smoothMPC_ds_aff + 944;
smoothMPC_FLOAT* smoothMPC_dllbcc59 = smoothMPC_dl_cc + 944;
smoothMPC_FLOAT* smoothMPC_dslbcc59 = smoothMPC_ds_cc + 944;
smoothMPC_FLOAT* smoothMPC_ccrhsl59 = smoothMPC_ccrhs + 944;
int smoothMPC_ubIdx59[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub59 = smoothMPC_l + 955;
smoothMPC_FLOAT* smoothMPC_sub59 = smoothMPC_s + 955;
smoothMPC_FLOAT* smoothMPC_lubbysub59 = smoothMPC_lbys + 955;
smoothMPC_FLOAT smoothMPC_riub59[5];
smoothMPC_FLOAT* smoothMPC_dlubaff59 = smoothMPC_dl_aff + 955;
smoothMPC_FLOAT* smoothMPC_dsubaff59 = smoothMPC_ds_aff + 955;
smoothMPC_FLOAT* smoothMPC_dlubcc59 = smoothMPC_dl_cc + 955;
smoothMPC_FLOAT* smoothMPC_dsubcc59 = smoothMPC_ds_cc + 955;
smoothMPC_FLOAT* smoothMPC_ccrhsub59 = smoothMPC_ccrhs + 955;
smoothMPC_FLOAT smoothMPC_Phi59[11];
smoothMPC_FLOAT smoothMPC_W59[11];
smoothMPC_FLOAT smoothMPC_Ysd59[9];
smoothMPC_FLOAT smoothMPC_Lsd59[9];
smoothMPC_FLOAT* smoothMPC_z60 = smoothMPC_z + 660;
smoothMPC_FLOAT* smoothMPC_dzaff60 = smoothMPC_dz_aff + 660;
smoothMPC_FLOAT* smoothMPC_dzcc60 = smoothMPC_dz_cc + 660;
smoothMPC_FLOAT* smoothMPC_rd60 = smoothMPC_rd + 660;
smoothMPC_FLOAT smoothMPC_Lbyrd60[11];
smoothMPC_FLOAT* smoothMPC_grad_cost60 = smoothMPC_grad_cost + 660;
smoothMPC_FLOAT* smoothMPC_grad_eq60 = smoothMPC_grad_eq + 660;
smoothMPC_FLOAT* smoothMPC_grad_ineq60 = smoothMPC_grad_ineq + 660;
smoothMPC_FLOAT smoothMPC_ctv60[11];
smoothMPC_FLOAT* smoothMPC_v60 = smoothMPC_v + 180;
smoothMPC_FLOAT smoothMPC_re60[3];
smoothMPC_FLOAT smoothMPC_beta60[3];
smoothMPC_FLOAT smoothMPC_betacc60[3];
smoothMPC_FLOAT* smoothMPC_dvaff60 = smoothMPC_dv_aff + 180;
smoothMPC_FLOAT* smoothMPC_dvcc60 = smoothMPC_dv_cc + 180;
smoothMPC_FLOAT smoothMPC_V60[33];
smoothMPC_FLOAT smoothMPC_Yd60[6];
smoothMPC_FLOAT smoothMPC_Ld60[6];
smoothMPC_FLOAT smoothMPC_yy60[3];
smoothMPC_FLOAT smoothMPC_bmy60[3];
int smoothMPC_lbIdx60[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb60 = smoothMPC_l + 960;
smoothMPC_FLOAT* smoothMPC_slb60 = smoothMPC_s + 960;
smoothMPC_FLOAT* smoothMPC_llbbyslb60 = smoothMPC_lbys + 960;
smoothMPC_FLOAT smoothMPC_rilb60[11];
smoothMPC_FLOAT* smoothMPC_dllbaff60 = smoothMPC_dl_aff + 960;
smoothMPC_FLOAT* smoothMPC_dslbaff60 = smoothMPC_ds_aff + 960;
smoothMPC_FLOAT* smoothMPC_dllbcc60 = smoothMPC_dl_cc + 960;
smoothMPC_FLOAT* smoothMPC_dslbcc60 = smoothMPC_ds_cc + 960;
smoothMPC_FLOAT* smoothMPC_ccrhsl60 = smoothMPC_ccrhs + 960;
int smoothMPC_ubIdx60[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub60 = smoothMPC_l + 971;
smoothMPC_FLOAT* smoothMPC_sub60 = smoothMPC_s + 971;
smoothMPC_FLOAT* smoothMPC_lubbysub60 = smoothMPC_lbys + 971;
smoothMPC_FLOAT smoothMPC_riub60[5];
smoothMPC_FLOAT* smoothMPC_dlubaff60 = smoothMPC_dl_aff + 971;
smoothMPC_FLOAT* smoothMPC_dsubaff60 = smoothMPC_ds_aff + 971;
smoothMPC_FLOAT* smoothMPC_dlubcc60 = smoothMPC_dl_cc + 971;
smoothMPC_FLOAT* smoothMPC_dsubcc60 = smoothMPC_ds_cc + 971;
smoothMPC_FLOAT* smoothMPC_ccrhsub60 = smoothMPC_ccrhs + 971;
smoothMPC_FLOAT smoothMPC_Phi60[11];
smoothMPC_FLOAT smoothMPC_W60[11];
smoothMPC_FLOAT smoothMPC_Ysd60[9];
smoothMPC_FLOAT smoothMPC_Lsd60[9];
smoothMPC_FLOAT* smoothMPC_z61 = smoothMPC_z + 671;
smoothMPC_FLOAT* smoothMPC_dzaff61 = smoothMPC_dz_aff + 671;
smoothMPC_FLOAT* smoothMPC_dzcc61 = smoothMPC_dz_cc + 671;
smoothMPC_FLOAT* smoothMPC_rd61 = smoothMPC_rd + 671;
smoothMPC_FLOAT smoothMPC_Lbyrd61[11];
smoothMPC_FLOAT* smoothMPC_grad_cost61 = smoothMPC_grad_cost + 671;
smoothMPC_FLOAT* smoothMPC_grad_eq61 = smoothMPC_grad_eq + 671;
smoothMPC_FLOAT* smoothMPC_grad_ineq61 = smoothMPC_grad_ineq + 671;
smoothMPC_FLOAT smoothMPC_ctv61[11];
smoothMPC_FLOAT* smoothMPC_v61 = smoothMPC_v + 183;
smoothMPC_FLOAT smoothMPC_re61[3];
smoothMPC_FLOAT smoothMPC_beta61[3];
smoothMPC_FLOAT smoothMPC_betacc61[3];
smoothMPC_FLOAT* smoothMPC_dvaff61 = smoothMPC_dv_aff + 183;
smoothMPC_FLOAT* smoothMPC_dvcc61 = smoothMPC_dv_cc + 183;
smoothMPC_FLOAT smoothMPC_V61[33];
smoothMPC_FLOAT smoothMPC_Yd61[6];
smoothMPC_FLOAT smoothMPC_Ld61[6];
smoothMPC_FLOAT smoothMPC_yy61[3];
smoothMPC_FLOAT smoothMPC_bmy61[3];
int smoothMPC_lbIdx61[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb61 = smoothMPC_l + 976;
smoothMPC_FLOAT* smoothMPC_slb61 = smoothMPC_s + 976;
smoothMPC_FLOAT* smoothMPC_llbbyslb61 = smoothMPC_lbys + 976;
smoothMPC_FLOAT smoothMPC_rilb61[11];
smoothMPC_FLOAT* smoothMPC_dllbaff61 = smoothMPC_dl_aff + 976;
smoothMPC_FLOAT* smoothMPC_dslbaff61 = smoothMPC_ds_aff + 976;
smoothMPC_FLOAT* smoothMPC_dllbcc61 = smoothMPC_dl_cc + 976;
smoothMPC_FLOAT* smoothMPC_dslbcc61 = smoothMPC_ds_cc + 976;
smoothMPC_FLOAT* smoothMPC_ccrhsl61 = smoothMPC_ccrhs + 976;
int smoothMPC_ubIdx61[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub61 = smoothMPC_l + 987;
smoothMPC_FLOAT* smoothMPC_sub61 = smoothMPC_s + 987;
smoothMPC_FLOAT* smoothMPC_lubbysub61 = smoothMPC_lbys + 987;
smoothMPC_FLOAT smoothMPC_riub61[5];
smoothMPC_FLOAT* smoothMPC_dlubaff61 = smoothMPC_dl_aff + 987;
smoothMPC_FLOAT* smoothMPC_dsubaff61 = smoothMPC_ds_aff + 987;
smoothMPC_FLOAT* smoothMPC_dlubcc61 = smoothMPC_dl_cc + 987;
smoothMPC_FLOAT* smoothMPC_dsubcc61 = smoothMPC_ds_cc + 987;
smoothMPC_FLOAT* smoothMPC_ccrhsub61 = smoothMPC_ccrhs + 987;
smoothMPC_FLOAT smoothMPC_Phi61[11];
smoothMPC_FLOAT smoothMPC_W61[11];
smoothMPC_FLOAT smoothMPC_Ysd61[9];
smoothMPC_FLOAT smoothMPC_Lsd61[9];
smoothMPC_FLOAT* smoothMPC_z62 = smoothMPC_z + 682;
smoothMPC_FLOAT* smoothMPC_dzaff62 = smoothMPC_dz_aff + 682;
smoothMPC_FLOAT* smoothMPC_dzcc62 = smoothMPC_dz_cc + 682;
smoothMPC_FLOAT* smoothMPC_rd62 = smoothMPC_rd + 682;
smoothMPC_FLOAT smoothMPC_Lbyrd62[11];
smoothMPC_FLOAT* smoothMPC_grad_cost62 = smoothMPC_grad_cost + 682;
smoothMPC_FLOAT* smoothMPC_grad_eq62 = smoothMPC_grad_eq + 682;
smoothMPC_FLOAT* smoothMPC_grad_ineq62 = smoothMPC_grad_ineq + 682;
smoothMPC_FLOAT smoothMPC_ctv62[11];
smoothMPC_FLOAT* smoothMPC_v62 = smoothMPC_v + 186;
smoothMPC_FLOAT smoothMPC_re62[3];
smoothMPC_FLOAT smoothMPC_beta62[3];
smoothMPC_FLOAT smoothMPC_betacc62[3];
smoothMPC_FLOAT* smoothMPC_dvaff62 = smoothMPC_dv_aff + 186;
smoothMPC_FLOAT* smoothMPC_dvcc62 = smoothMPC_dv_cc + 186;
smoothMPC_FLOAT smoothMPC_V62[33];
smoothMPC_FLOAT smoothMPC_Yd62[6];
smoothMPC_FLOAT smoothMPC_Ld62[6];
smoothMPC_FLOAT smoothMPC_yy62[3];
smoothMPC_FLOAT smoothMPC_bmy62[3];
int smoothMPC_lbIdx62[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb62 = smoothMPC_l + 992;
smoothMPC_FLOAT* smoothMPC_slb62 = smoothMPC_s + 992;
smoothMPC_FLOAT* smoothMPC_llbbyslb62 = smoothMPC_lbys + 992;
smoothMPC_FLOAT smoothMPC_rilb62[11];
smoothMPC_FLOAT* smoothMPC_dllbaff62 = smoothMPC_dl_aff + 992;
smoothMPC_FLOAT* smoothMPC_dslbaff62 = smoothMPC_ds_aff + 992;
smoothMPC_FLOAT* smoothMPC_dllbcc62 = smoothMPC_dl_cc + 992;
smoothMPC_FLOAT* smoothMPC_dslbcc62 = smoothMPC_ds_cc + 992;
smoothMPC_FLOAT* smoothMPC_ccrhsl62 = smoothMPC_ccrhs + 992;
int smoothMPC_ubIdx62[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub62 = smoothMPC_l + 1003;
smoothMPC_FLOAT* smoothMPC_sub62 = smoothMPC_s + 1003;
smoothMPC_FLOAT* smoothMPC_lubbysub62 = smoothMPC_lbys + 1003;
smoothMPC_FLOAT smoothMPC_riub62[5];
smoothMPC_FLOAT* smoothMPC_dlubaff62 = smoothMPC_dl_aff + 1003;
smoothMPC_FLOAT* smoothMPC_dsubaff62 = smoothMPC_ds_aff + 1003;
smoothMPC_FLOAT* smoothMPC_dlubcc62 = smoothMPC_dl_cc + 1003;
smoothMPC_FLOAT* smoothMPC_dsubcc62 = smoothMPC_ds_cc + 1003;
smoothMPC_FLOAT* smoothMPC_ccrhsub62 = smoothMPC_ccrhs + 1003;
smoothMPC_FLOAT smoothMPC_Phi62[11];
smoothMPC_FLOAT smoothMPC_W62[11];
smoothMPC_FLOAT smoothMPC_Ysd62[9];
smoothMPC_FLOAT smoothMPC_Lsd62[9];
smoothMPC_FLOAT* smoothMPC_z63 = smoothMPC_z + 693;
smoothMPC_FLOAT* smoothMPC_dzaff63 = smoothMPC_dz_aff + 693;
smoothMPC_FLOAT* smoothMPC_dzcc63 = smoothMPC_dz_cc + 693;
smoothMPC_FLOAT* smoothMPC_rd63 = smoothMPC_rd + 693;
smoothMPC_FLOAT smoothMPC_Lbyrd63[11];
smoothMPC_FLOAT* smoothMPC_grad_cost63 = smoothMPC_grad_cost + 693;
smoothMPC_FLOAT* smoothMPC_grad_eq63 = smoothMPC_grad_eq + 693;
smoothMPC_FLOAT* smoothMPC_grad_ineq63 = smoothMPC_grad_ineq + 693;
smoothMPC_FLOAT smoothMPC_ctv63[11];
smoothMPC_FLOAT* smoothMPC_v63 = smoothMPC_v + 189;
smoothMPC_FLOAT smoothMPC_re63[3];
smoothMPC_FLOAT smoothMPC_beta63[3];
smoothMPC_FLOAT smoothMPC_betacc63[3];
smoothMPC_FLOAT* smoothMPC_dvaff63 = smoothMPC_dv_aff + 189;
smoothMPC_FLOAT* smoothMPC_dvcc63 = smoothMPC_dv_cc + 189;
smoothMPC_FLOAT smoothMPC_V63[33];
smoothMPC_FLOAT smoothMPC_Yd63[6];
smoothMPC_FLOAT smoothMPC_Ld63[6];
smoothMPC_FLOAT smoothMPC_yy63[3];
smoothMPC_FLOAT smoothMPC_bmy63[3];
int smoothMPC_lbIdx63[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb63 = smoothMPC_l + 1008;
smoothMPC_FLOAT* smoothMPC_slb63 = smoothMPC_s + 1008;
smoothMPC_FLOAT* smoothMPC_llbbyslb63 = smoothMPC_lbys + 1008;
smoothMPC_FLOAT smoothMPC_rilb63[11];
smoothMPC_FLOAT* smoothMPC_dllbaff63 = smoothMPC_dl_aff + 1008;
smoothMPC_FLOAT* smoothMPC_dslbaff63 = smoothMPC_ds_aff + 1008;
smoothMPC_FLOAT* smoothMPC_dllbcc63 = smoothMPC_dl_cc + 1008;
smoothMPC_FLOAT* smoothMPC_dslbcc63 = smoothMPC_ds_cc + 1008;
smoothMPC_FLOAT* smoothMPC_ccrhsl63 = smoothMPC_ccrhs + 1008;
int smoothMPC_ubIdx63[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub63 = smoothMPC_l + 1019;
smoothMPC_FLOAT* smoothMPC_sub63 = smoothMPC_s + 1019;
smoothMPC_FLOAT* smoothMPC_lubbysub63 = smoothMPC_lbys + 1019;
smoothMPC_FLOAT smoothMPC_riub63[5];
smoothMPC_FLOAT* smoothMPC_dlubaff63 = smoothMPC_dl_aff + 1019;
smoothMPC_FLOAT* smoothMPC_dsubaff63 = smoothMPC_ds_aff + 1019;
smoothMPC_FLOAT* smoothMPC_dlubcc63 = smoothMPC_dl_cc + 1019;
smoothMPC_FLOAT* smoothMPC_dsubcc63 = smoothMPC_ds_cc + 1019;
smoothMPC_FLOAT* smoothMPC_ccrhsub63 = smoothMPC_ccrhs + 1019;
smoothMPC_FLOAT smoothMPC_Phi63[11];
smoothMPC_FLOAT smoothMPC_W63[11];
smoothMPC_FLOAT smoothMPC_Ysd63[9];
smoothMPC_FLOAT smoothMPC_Lsd63[9];
smoothMPC_FLOAT* smoothMPC_z64 = smoothMPC_z + 704;
smoothMPC_FLOAT* smoothMPC_dzaff64 = smoothMPC_dz_aff + 704;
smoothMPC_FLOAT* smoothMPC_dzcc64 = smoothMPC_dz_cc + 704;
smoothMPC_FLOAT* smoothMPC_rd64 = smoothMPC_rd + 704;
smoothMPC_FLOAT smoothMPC_Lbyrd64[11];
smoothMPC_FLOAT* smoothMPC_grad_cost64 = smoothMPC_grad_cost + 704;
smoothMPC_FLOAT* smoothMPC_grad_eq64 = smoothMPC_grad_eq + 704;
smoothMPC_FLOAT* smoothMPC_grad_ineq64 = smoothMPC_grad_ineq + 704;
smoothMPC_FLOAT smoothMPC_ctv64[11];
smoothMPC_FLOAT* smoothMPC_v64 = smoothMPC_v + 192;
smoothMPC_FLOAT smoothMPC_re64[3];
smoothMPC_FLOAT smoothMPC_beta64[3];
smoothMPC_FLOAT smoothMPC_betacc64[3];
smoothMPC_FLOAT* smoothMPC_dvaff64 = smoothMPC_dv_aff + 192;
smoothMPC_FLOAT* smoothMPC_dvcc64 = smoothMPC_dv_cc + 192;
smoothMPC_FLOAT smoothMPC_V64[33];
smoothMPC_FLOAT smoothMPC_Yd64[6];
smoothMPC_FLOAT smoothMPC_Ld64[6];
smoothMPC_FLOAT smoothMPC_yy64[3];
smoothMPC_FLOAT smoothMPC_bmy64[3];
int smoothMPC_lbIdx64[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb64 = smoothMPC_l + 1024;
smoothMPC_FLOAT* smoothMPC_slb64 = smoothMPC_s + 1024;
smoothMPC_FLOAT* smoothMPC_llbbyslb64 = smoothMPC_lbys + 1024;
smoothMPC_FLOAT smoothMPC_rilb64[11];
smoothMPC_FLOAT* smoothMPC_dllbaff64 = smoothMPC_dl_aff + 1024;
smoothMPC_FLOAT* smoothMPC_dslbaff64 = smoothMPC_ds_aff + 1024;
smoothMPC_FLOAT* smoothMPC_dllbcc64 = smoothMPC_dl_cc + 1024;
smoothMPC_FLOAT* smoothMPC_dslbcc64 = smoothMPC_ds_cc + 1024;
smoothMPC_FLOAT* smoothMPC_ccrhsl64 = smoothMPC_ccrhs + 1024;
int smoothMPC_ubIdx64[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub64 = smoothMPC_l + 1035;
smoothMPC_FLOAT* smoothMPC_sub64 = smoothMPC_s + 1035;
smoothMPC_FLOAT* smoothMPC_lubbysub64 = smoothMPC_lbys + 1035;
smoothMPC_FLOAT smoothMPC_riub64[5];
smoothMPC_FLOAT* smoothMPC_dlubaff64 = smoothMPC_dl_aff + 1035;
smoothMPC_FLOAT* smoothMPC_dsubaff64 = smoothMPC_ds_aff + 1035;
smoothMPC_FLOAT* smoothMPC_dlubcc64 = smoothMPC_dl_cc + 1035;
smoothMPC_FLOAT* smoothMPC_dsubcc64 = smoothMPC_ds_cc + 1035;
smoothMPC_FLOAT* smoothMPC_ccrhsub64 = smoothMPC_ccrhs + 1035;
smoothMPC_FLOAT smoothMPC_Phi64[11];
smoothMPC_FLOAT smoothMPC_W64[11];
smoothMPC_FLOAT smoothMPC_Ysd64[9];
smoothMPC_FLOAT smoothMPC_Lsd64[9];
smoothMPC_FLOAT* smoothMPC_z65 = smoothMPC_z + 715;
smoothMPC_FLOAT* smoothMPC_dzaff65 = smoothMPC_dz_aff + 715;
smoothMPC_FLOAT* smoothMPC_dzcc65 = smoothMPC_dz_cc + 715;
smoothMPC_FLOAT* smoothMPC_rd65 = smoothMPC_rd + 715;
smoothMPC_FLOAT smoothMPC_Lbyrd65[11];
smoothMPC_FLOAT* smoothMPC_grad_cost65 = smoothMPC_grad_cost + 715;
smoothMPC_FLOAT* smoothMPC_grad_eq65 = smoothMPC_grad_eq + 715;
smoothMPC_FLOAT* smoothMPC_grad_ineq65 = smoothMPC_grad_ineq + 715;
smoothMPC_FLOAT smoothMPC_ctv65[11];
smoothMPC_FLOAT* smoothMPC_v65 = smoothMPC_v + 195;
smoothMPC_FLOAT smoothMPC_re65[3];
smoothMPC_FLOAT smoothMPC_beta65[3];
smoothMPC_FLOAT smoothMPC_betacc65[3];
smoothMPC_FLOAT* smoothMPC_dvaff65 = smoothMPC_dv_aff + 195;
smoothMPC_FLOAT* smoothMPC_dvcc65 = smoothMPC_dv_cc + 195;
smoothMPC_FLOAT smoothMPC_V65[33];
smoothMPC_FLOAT smoothMPC_Yd65[6];
smoothMPC_FLOAT smoothMPC_Ld65[6];
smoothMPC_FLOAT smoothMPC_yy65[3];
smoothMPC_FLOAT smoothMPC_bmy65[3];
int smoothMPC_lbIdx65[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb65 = smoothMPC_l + 1040;
smoothMPC_FLOAT* smoothMPC_slb65 = smoothMPC_s + 1040;
smoothMPC_FLOAT* smoothMPC_llbbyslb65 = smoothMPC_lbys + 1040;
smoothMPC_FLOAT smoothMPC_rilb65[11];
smoothMPC_FLOAT* smoothMPC_dllbaff65 = smoothMPC_dl_aff + 1040;
smoothMPC_FLOAT* smoothMPC_dslbaff65 = smoothMPC_ds_aff + 1040;
smoothMPC_FLOAT* smoothMPC_dllbcc65 = smoothMPC_dl_cc + 1040;
smoothMPC_FLOAT* smoothMPC_dslbcc65 = smoothMPC_ds_cc + 1040;
smoothMPC_FLOAT* smoothMPC_ccrhsl65 = smoothMPC_ccrhs + 1040;
int smoothMPC_ubIdx65[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub65 = smoothMPC_l + 1051;
smoothMPC_FLOAT* smoothMPC_sub65 = smoothMPC_s + 1051;
smoothMPC_FLOAT* smoothMPC_lubbysub65 = smoothMPC_lbys + 1051;
smoothMPC_FLOAT smoothMPC_riub65[5];
smoothMPC_FLOAT* smoothMPC_dlubaff65 = smoothMPC_dl_aff + 1051;
smoothMPC_FLOAT* smoothMPC_dsubaff65 = smoothMPC_ds_aff + 1051;
smoothMPC_FLOAT* smoothMPC_dlubcc65 = smoothMPC_dl_cc + 1051;
smoothMPC_FLOAT* smoothMPC_dsubcc65 = smoothMPC_ds_cc + 1051;
smoothMPC_FLOAT* smoothMPC_ccrhsub65 = smoothMPC_ccrhs + 1051;
smoothMPC_FLOAT smoothMPC_Phi65[11];
smoothMPC_FLOAT smoothMPC_W65[11];
smoothMPC_FLOAT smoothMPC_Ysd65[9];
smoothMPC_FLOAT smoothMPC_Lsd65[9];
smoothMPC_FLOAT* smoothMPC_z66 = smoothMPC_z + 726;
smoothMPC_FLOAT* smoothMPC_dzaff66 = smoothMPC_dz_aff + 726;
smoothMPC_FLOAT* smoothMPC_dzcc66 = smoothMPC_dz_cc + 726;
smoothMPC_FLOAT* smoothMPC_rd66 = smoothMPC_rd + 726;
smoothMPC_FLOAT smoothMPC_Lbyrd66[11];
smoothMPC_FLOAT* smoothMPC_grad_cost66 = smoothMPC_grad_cost + 726;
smoothMPC_FLOAT* smoothMPC_grad_eq66 = smoothMPC_grad_eq + 726;
smoothMPC_FLOAT* smoothMPC_grad_ineq66 = smoothMPC_grad_ineq + 726;
smoothMPC_FLOAT smoothMPC_ctv66[11];
smoothMPC_FLOAT* smoothMPC_v66 = smoothMPC_v + 198;
smoothMPC_FLOAT smoothMPC_re66[3];
smoothMPC_FLOAT smoothMPC_beta66[3];
smoothMPC_FLOAT smoothMPC_betacc66[3];
smoothMPC_FLOAT* smoothMPC_dvaff66 = smoothMPC_dv_aff + 198;
smoothMPC_FLOAT* smoothMPC_dvcc66 = smoothMPC_dv_cc + 198;
smoothMPC_FLOAT smoothMPC_V66[33];
smoothMPC_FLOAT smoothMPC_Yd66[6];
smoothMPC_FLOAT smoothMPC_Ld66[6];
smoothMPC_FLOAT smoothMPC_yy66[3];
smoothMPC_FLOAT smoothMPC_bmy66[3];
int smoothMPC_lbIdx66[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb66 = smoothMPC_l + 1056;
smoothMPC_FLOAT* smoothMPC_slb66 = smoothMPC_s + 1056;
smoothMPC_FLOAT* smoothMPC_llbbyslb66 = smoothMPC_lbys + 1056;
smoothMPC_FLOAT smoothMPC_rilb66[11];
smoothMPC_FLOAT* smoothMPC_dllbaff66 = smoothMPC_dl_aff + 1056;
smoothMPC_FLOAT* smoothMPC_dslbaff66 = smoothMPC_ds_aff + 1056;
smoothMPC_FLOAT* smoothMPC_dllbcc66 = smoothMPC_dl_cc + 1056;
smoothMPC_FLOAT* smoothMPC_dslbcc66 = smoothMPC_ds_cc + 1056;
smoothMPC_FLOAT* smoothMPC_ccrhsl66 = smoothMPC_ccrhs + 1056;
int smoothMPC_ubIdx66[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub66 = smoothMPC_l + 1067;
smoothMPC_FLOAT* smoothMPC_sub66 = smoothMPC_s + 1067;
smoothMPC_FLOAT* smoothMPC_lubbysub66 = smoothMPC_lbys + 1067;
smoothMPC_FLOAT smoothMPC_riub66[5];
smoothMPC_FLOAT* smoothMPC_dlubaff66 = smoothMPC_dl_aff + 1067;
smoothMPC_FLOAT* smoothMPC_dsubaff66 = smoothMPC_ds_aff + 1067;
smoothMPC_FLOAT* smoothMPC_dlubcc66 = smoothMPC_dl_cc + 1067;
smoothMPC_FLOAT* smoothMPC_dsubcc66 = smoothMPC_ds_cc + 1067;
smoothMPC_FLOAT* smoothMPC_ccrhsub66 = smoothMPC_ccrhs + 1067;
smoothMPC_FLOAT smoothMPC_Phi66[11];
smoothMPC_FLOAT smoothMPC_W66[11];
smoothMPC_FLOAT smoothMPC_Ysd66[9];
smoothMPC_FLOAT smoothMPC_Lsd66[9];
smoothMPC_FLOAT* smoothMPC_z67 = smoothMPC_z + 737;
smoothMPC_FLOAT* smoothMPC_dzaff67 = smoothMPC_dz_aff + 737;
smoothMPC_FLOAT* smoothMPC_dzcc67 = smoothMPC_dz_cc + 737;
smoothMPC_FLOAT* smoothMPC_rd67 = smoothMPC_rd + 737;
smoothMPC_FLOAT smoothMPC_Lbyrd67[11];
smoothMPC_FLOAT* smoothMPC_grad_cost67 = smoothMPC_grad_cost + 737;
smoothMPC_FLOAT* smoothMPC_grad_eq67 = smoothMPC_grad_eq + 737;
smoothMPC_FLOAT* smoothMPC_grad_ineq67 = smoothMPC_grad_ineq + 737;
smoothMPC_FLOAT smoothMPC_ctv67[11];
smoothMPC_FLOAT* smoothMPC_v67 = smoothMPC_v + 201;
smoothMPC_FLOAT smoothMPC_re67[3];
smoothMPC_FLOAT smoothMPC_beta67[3];
smoothMPC_FLOAT smoothMPC_betacc67[3];
smoothMPC_FLOAT* smoothMPC_dvaff67 = smoothMPC_dv_aff + 201;
smoothMPC_FLOAT* smoothMPC_dvcc67 = smoothMPC_dv_cc + 201;
smoothMPC_FLOAT smoothMPC_V67[33];
smoothMPC_FLOAT smoothMPC_Yd67[6];
smoothMPC_FLOAT smoothMPC_Ld67[6];
smoothMPC_FLOAT smoothMPC_yy67[3];
smoothMPC_FLOAT smoothMPC_bmy67[3];
int smoothMPC_lbIdx67[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb67 = smoothMPC_l + 1072;
smoothMPC_FLOAT* smoothMPC_slb67 = smoothMPC_s + 1072;
smoothMPC_FLOAT* smoothMPC_llbbyslb67 = smoothMPC_lbys + 1072;
smoothMPC_FLOAT smoothMPC_rilb67[11];
smoothMPC_FLOAT* smoothMPC_dllbaff67 = smoothMPC_dl_aff + 1072;
smoothMPC_FLOAT* smoothMPC_dslbaff67 = smoothMPC_ds_aff + 1072;
smoothMPC_FLOAT* smoothMPC_dllbcc67 = smoothMPC_dl_cc + 1072;
smoothMPC_FLOAT* smoothMPC_dslbcc67 = smoothMPC_ds_cc + 1072;
smoothMPC_FLOAT* smoothMPC_ccrhsl67 = smoothMPC_ccrhs + 1072;
int smoothMPC_ubIdx67[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub67 = smoothMPC_l + 1083;
smoothMPC_FLOAT* smoothMPC_sub67 = smoothMPC_s + 1083;
smoothMPC_FLOAT* smoothMPC_lubbysub67 = smoothMPC_lbys + 1083;
smoothMPC_FLOAT smoothMPC_riub67[5];
smoothMPC_FLOAT* smoothMPC_dlubaff67 = smoothMPC_dl_aff + 1083;
smoothMPC_FLOAT* smoothMPC_dsubaff67 = smoothMPC_ds_aff + 1083;
smoothMPC_FLOAT* smoothMPC_dlubcc67 = smoothMPC_dl_cc + 1083;
smoothMPC_FLOAT* smoothMPC_dsubcc67 = smoothMPC_ds_cc + 1083;
smoothMPC_FLOAT* smoothMPC_ccrhsub67 = smoothMPC_ccrhs + 1083;
smoothMPC_FLOAT smoothMPC_Phi67[11];
smoothMPC_FLOAT smoothMPC_W67[11];
smoothMPC_FLOAT smoothMPC_Ysd67[9];
smoothMPC_FLOAT smoothMPC_Lsd67[9];
smoothMPC_FLOAT* smoothMPC_z68 = smoothMPC_z + 748;
smoothMPC_FLOAT* smoothMPC_dzaff68 = smoothMPC_dz_aff + 748;
smoothMPC_FLOAT* smoothMPC_dzcc68 = smoothMPC_dz_cc + 748;
smoothMPC_FLOAT* smoothMPC_rd68 = smoothMPC_rd + 748;
smoothMPC_FLOAT smoothMPC_Lbyrd68[11];
smoothMPC_FLOAT* smoothMPC_grad_cost68 = smoothMPC_grad_cost + 748;
smoothMPC_FLOAT* smoothMPC_grad_eq68 = smoothMPC_grad_eq + 748;
smoothMPC_FLOAT* smoothMPC_grad_ineq68 = smoothMPC_grad_ineq + 748;
smoothMPC_FLOAT smoothMPC_ctv68[11];
smoothMPC_FLOAT* smoothMPC_v68 = smoothMPC_v + 204;
smoothMPC_FLOAT smoothMPC_re68[3];
smoothMPC_FLOAT smoothMPC_beta68[3];
smoothMPC_FLOAT smoothMPC_betacc68[3];
smoothMPC_FLOAT* smoothMPC_dvaff68 = smoothMPC_dv_aff + 204;
smoothMPC_FLOAT* smoothMPC_dvcc68 = smoothMPC_dv_cc + 204;
smoothMPC_FLOAT smoothMPC_V68[33];
smoothMPC_FLOAT smoothMPC_Yd68[6];
smoothMPC_FLOAT smoothMPC_Ld68[6];
smoothMPC_FLOAT smoothMPC_yy68[3];
smoothMPC_FLOAT smoothMPC_bmy68[3];
int smoothMPC_lbIdx68[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb68 = smoothMPC_l + 1088;
smoothMPC_FLOAT* smoothMPC_slb68 = smoothMPC_s + 1088;
smoothMPC_FLOAT* smoothMPC_llbbyslb68 = smoothMPC_lbys + 1088;
smoothMPC_FLOAT smoothMPC_rilb68[11];
smoothMPC_FLOAT* smoothMPC_dllbaff68 = smoothMPC_dl_aff + 1088;
smoothMPC_FLOAT* smoothMPC_dslbaff68 = smoothMPC_ds_aff + 1088;
smoothMPC_FLOAT* smoothMPC_dllbcc68 = smoothMPC_dl_cc + 1088;
smoothMPC_FLOAT* smoothMPC_dslbcc68 = smoothMPC_ds_cc + 1088;
smoothMPC_FLOAT* smoothMPC_ccrhsl68 = smoothMPC_ccrhs + 1088;
int smoothMPC_ubIdx68[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub68 = smoothMPC_l + 1099;
smoothMPC_FLOAT* smoothMPC_sub68 = smoothMPC_s + 1099;
smoothMPC_FLOAT* smoothMPC_lubbysub68 = smoothMPC_lbys + 1099;
smoothMPC_FLOAT smoothMPC_riub68[5];
smoothMPC_FLOAT* smoothMPC_dlubaff68 = smoothMPC_dl_aff + 1099;
smoothMPC_FLOAT* smoothMPC_dsubaff68 = smoothMPC_ds_aff + 1099;
smoothMPC_FLOAT* smoothMPC_dlubcc68 = smoothMPC_dl_cc + 1099;
smoothMPC_FLOAT* smoothMPC_dsubcc68 = smoothMPC_ds_cc + 1099;
smoothMPC_FLOAT* smoothMPC_ccrhsub68 = smoothMPC_ccrhs + 1099;
smoothMPC_FLOAT smoothMPC_Phi68[11];
smoothMPC_FLOAT smoothMPC_W68[11];
smoothMPC_FLOAT smoothMPC_Ysd68[9];
smoothMPC_FLOAT smoothMPC_Lsd68[9];
smoothMPC_FLOAT* smoothMPC_z69 = smoothMPC_z + 759;
smoothMPC_FLOAT* smoothMPC_dzaff69 = smoothMPC_dz_aff + 759;
smoothMPC_FLOAT* smoothMPC_dzcc69 = smoothMPC_dz_cc + 759;
smoothMPC_FLOAT* smoothMPC_rd69 = smoothMPC_rd + 759;
smoothMPC_FLOAT smoothMPC_Lbyrd69[11];
smoothMPC_FLOAT* smoothMPC_grad_cost69 = smoothMPC_grad_cost + 759;
smoothMPC_FLOAT* smoothMPC_grad_eq69 = smoothMPC_grad_eq + 759;
smoothMPC_FLOAT* smoothMPC_grad_ineq69 = smoothMPC_grad_ineq + 759;
smoothMPC_FLOAT smoothMPC_ctv69[11];
smoothMPC_FLOAT* smoothMPC_v69 = smoothMPC_v + 207;
smoothMPC_FLOAT smoothMPC_re69[3];
smoothMPC_FLOAT smoothMPC_beta69[3];
smoothMPC_FLOAT smoothMPC_betacc69[3];
smoothMPC_FLOAT* smoothMPC_dvaff69 = smoothMPC_dv_aff + 207;
smoothMPC_FLOAT* smoothMPC_dvcc69 = smoothMPC_dv_cc + 207;
smoothMPC_FLOAT smoothMPC_V69[33];
smoothMPC_FLOAT smoothMPC_Yd69[6];
smoothMPC_FLOAT smoothMPC_Ld69[6];
smoothMPC_FLOAT smoothMPC_yy69[3];
smoothMPC_FLOAT smoothMPC_bmy69[3];
int smoothMPC_lbIdx69[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb69 = smoothMPC_l + 1104;
smoothMPC_FLOAT* smoothMPC_slb69 = smoothMPC_s + 1104;
smoothMPC_FLOAT* smoothMPC_llbbyslb69 = smoothMPC_lbys + 1104;
smoothMPC_FLOAT smoothMPC_rilb69[11];
smoothMPC_FLOAT* smoothMPC_dllbaff69 = smoothMPC_dl_aff + 1104;
smoothMPC_FLOAT* smoothMPC_dslbaff69 = smoothMPC_ds_aff + 1104;
smoothMPC_FLOAT* smoothMPC_dllbcc69 = smoothMPC_dl_cc + 1104;
smoothMPC_FLOAT* smoothMPC_dslbcc69 = smoothMPC_ds_cc + 1104;
smoothMPC_FLOAT* smoothMPC_ccrhsl69 = smoothMPC_ccrhs + 1104;
int smoothMPC_ubIdx69[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub69 = smoothMPC_l + 1115;
smoothMPC_FLOAT* smoothMPC_sub69 = smoothMPC_s + 1115;
smoothMPC_FLOAT* smoothMPC_lubbysub69 = smoothMPC_lbys + 1115;
smoothMPC_FLOAT smoothMPC_riub69[5];
smoothMPC_FLOAT* smoothMPC_dlubaff69 = smoothMPC_dl_aff + 1115;
smoothMPC_FLOAT* smoothMPC_dsubaff69 = smoothMPC_ds_aff + 1115;
smoothMPC_FLOAT* smoothMPC_dlubcc69 = smoothMPC_dl_cc + 1115;
smoothMPC_FLOAT* smoothMPC_dsubcc69 = smoothMPC_ds_cc + 1115;
smoothMPC_FLOAT* smoothMPC_ccrhsub69 = smoothMPC_ccrhs + 1115;
smoothMPC_FLOAT smoothMPC_Phi69[11];
smoothMPC_FLOAT smoothMPC_W69[11];
smoothMPC_FLOAT smoothMPC_Ysd69[9];
smoothMPC_FLOAT smoothMPC_Lsd69[9];
smoothMPC_FLOAT* smoothMPC_z70 = smoothMPC_z + 770;
smoothMPC_FLOAT* smoothMPC_dzaff70 = smoothMPC_dz_aff + 770;
smoothMPC_FLOAT* smoothMPC_dzcc70 = smoothMPC_dz_cc + 770;
smoothMPC_FLOAT* smoothMPC_rd70 = smoothMPC_rd + 770;
smoothMPC_FLOAT smoothMPC_Lbyrd70[11];
smoothMPC_FLOAT* smoothMPC_grad_cost70 = smoothMPC_grad_cost + 770;
smoothMPC_FLOAT* smoothMPC_grad_eq70 = smoothMPC_grad_eq + 770;
smoothMPC_FLOAT* smoothMPC_grad_ineq70 = smoothMPC_grad_ineq + 770;
smoothMPC_FLOAT smoothMPC_ctv70[11];
smoothMPC_FLOAT* smoothMPC_v70 = smoothMPC_v + 210;
smoothMPC_FLOAT smoothMPC_re70[3];
smoothMPC_FLOAT smoothMPC_beta70[3];
smoothMPC_FLOAT smoothMPC_betacc70[3];
smoothMPC_FLOAT* smoothMPC_dvaff70 = smoothMPC_dv_aff + 210;
smoothMPC_FLOAT* smoothMPC_dvcc70 = smoothMPC_dv_cc + 210;
smoothMPC_FLOAT smoothMPC_V70[33];
smoothMPC_FLOAT smoothMPC_Yd70[6];
smoothMPC_FLOAT smoothMPC_Ld70[6];
smoothMPC_FLOAT smoothMPC_yy70[3];
smoothMPC_FLOAT smoothMPC_bmy70[3];
int smoothMPC_lbIdx70[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb70 = smoothMPC_l + 1120;
smoothMPC_FLOAT* smoothMPC_slb70 = smoothMPC_s + 1120;
smoothMPC_FLOAT* smoothMPC_llbbyslb70 = smoothMPC_lbys + 1120;
smoothMPC_FLOAT smoothMPC_rilb70[11];
smoothMPC_FLOAT* smoothMPC_dllbaff70 = smoothMPC_dl_aff + 1120;
smoothMPC_FLOAT* smoothMPC_dslbaff70 = smoothMPC_ds_aff + 1120;
smoothMPC_FLOAT* smoothMPC_dllbcc70 = smoothMPC_dl_cc + 1120;
smoothMPC_FLOAT* smoothMPC_dslbcc70 = smoothMPC_ds_cc + 1120;
smoothMPC_FLOAT* smoothMPC_ccrhsl70 = smoothMPC_ccrhs + 1120;
int smoothMPC_ubIdx70[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub70 = smoothMPC_l + 1131;
smoothMPC_FLOAT* smoothMPC_sub70 = smoothMPC_s + 1131;
smoothMPC_FLOAT* smoothMPC_lubbysub70 = smoothMPC_lbys + 1131;
smoothMPC_FLOAT smoothMPC_riub70[5];
smoothMPC_FLOAT* smoothMPC_dlubaff70 = smoothMPC_dl_aff + 1131;
smoothMPC_FLOAT* smoothMPC_dsubaff70 = smoothMPC_ds_aff + 1131;
smoothMPC_FLOAT* smoothMPC_dlubcc70 = smoothMPC_dl_cc + 1131;
smoothMPC_FLOAT* smoothMPC_dsubcc70 = smoothMPC_ds_cc + 1131;
smoothMPC_FLOAT* smoothMPC_ccrhsub70 = smoothMPC_ccrhs + 1131;
smoothMPC_FLOAT smoothMPC_Phi70[11];
smoothMPC_FLOAT smoothMPC_W70[11];
smoothMPC_FLOAT smoothMPC_Ysd70[9];
smoothMPC_FLOAT smoothMPC_Lsd70[9];
smoothMPC_FLOAT* smoothMPC_z71 = smoothMPC_z + 781;
smoothMPC_FLOAT* smoothMPC_dzaff71 = smoothMPC_dz_aff + 781;
smoothMPC_FLOAT* smoothMPC_dzcc71 = smoothMPC_dz_cc + 781;
smoothMPC_FLOAT* smoothMPC_rd71 = smoothMPC_rd + 781;
smoothMPC_FLOAT smoothMPC_Lbyrd71[11];
smoothMPC_FLOAT* smoothMPC_grad_cost71 = smoothMPC_grad_cost + 781;
smoothMPC_FLOAT* smoothMPC_grad_eq71 = smoothMPC_grad_eq + 781;
smoothMPC_FLOAT* smoothMPC_grad_ineq71 = smoothMPC_grad_ineq + 781;
smoothMPC_FLOAT smoothMPC_ctv71[11];
smoothMPC_FLOAT* smoothMPC_v71 = smoothMPC_v + 213;
smoothMPC_FLOAT smoothMPC_re71[3];
smoothMPC_FLOAT smoothMPC_beta71[3];
smoothMPC_FLOAT smoothMPC_betacc71[3];
smoothMPC_FLOAT* smoothMPC_dvaff71 = smoothMPC_dv_aff + 213;
smoothMPC_FLOAT* smoothMPC_dvcc71 = smoothMPC_dv_cc + 213;
smoothMPC_FLOAT smoothMPC_V71[33];
smoothMPC_FLOAT smoothMPC_Yd71[6];
smoothMPC_FLOAT smoothMPC_Ld71[6];
smoothMPC_FLOAT smoothMPC_yy71[3];
smoothMPC_FLOAT smoothMPC_bmy71[3];
int smoothMPC_lbIdx71[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb71 = smoothMPC_l + 1136;
smoothMPC_FLOAT* smoothMPC_slb71 = smoothMPC_s + 1136;
smoothMPC_FLOAT* smoothMPC_llbbyslb71 = smoothMPC_lbys + 1136;
smoothMPC_FLOAT smoothMPC_rilb71[11];
smoothMPC_FLOAT* smoothMPC_dllbaff71 = smoothMPC_dl_aff + 1136;
smoothMPC_FLOAT* smoothMPC_dslbaff71 = smoothMPC_ds_aff + 1136;
smoothMPC_FLOAT* smoothMPC_dllbcc71 = smoothMPC_dl_cc + 1136;
smoothMPC_FLOAT* smoothMPC_dslbcc71 = smoothMPC_ds_cc + 1136;
smoothMPC_FLOAT* smoothMPC_ccrhsl71 = smoothMPC_ccrhs + 1136;
int smoothMPC_ubIdx71[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub71 = smoothMPC_l + 1147;
smoothMPC_FLOAT* smoothMPC_sub71 = smoothMPC_s + 1147;
smoothMPC_FLOAT* smoothMPC_lubbysub71 = smoothMPC_lbys + 1147;
smoothMPC_FLOAT smoothMPC_riub71[5];
smoothMPC_FLOAT* smoothMPC_dlubaff71 = smoothMPC_dl_aff + 1147;
smoothMPC_FLOAT* smoothMPC_dsubaff71 = smoothMPC_ds_aff + 1147;
smoothMPC_FLOAT* smoothMPC_dlubcc71 = smoothMPC_dl_cc + 1147;
smoothMPC_FLOAT* smoothMPC_dsubcc71 = smoothMPC_ds_cc + 1147;
smoothMPC_FLOAT* smoothMPC_ccrhsub71 = smoothMPC_ccrhs + 1147;
smoothMPC_FLOAT smoothMPC_Phi71[11];
smoothMPC_FLOAT smoothMPC_W71[11];
smoothMPC_FLOAT smoothMPC_Ysd71[9];
smoothMPC_FLOAT smoothMPC_Lsd71[9];
smoothMPC_FLOAT* smoothMPC_z72 = smoothMPC_z + 792;
smoothMPC_FLOAT* smoothMPC_dzaff72 = smoothMPC_dz_aff + 792;
smoothMPC_FLOAT* smoothMPC_dzcc72 = smoothMPC_dz_cc + 792;
smoothMPC_FLOAT* smoothMPC_rd72 = smoothMPC_rd + 792;
smoothMPC_FLOAT smoothMPC_Lbyrd72[11];
smoothMPC_FLOAT* smoothMPC_grad_cost72 = smoothMPC_grad_cost + 792;
smoothMPC_FLOAT* smoothMPC_grad_eq72 = smoothMPC_grad_eq + 792;
smoothMPC_FLOAT* smoothMPC_grad_ineq72 = smoothMPC_grad_ineq + 792;
smoothMPC_FLOAT smoothMPC_ctv72[11];
smoothMPC_FLOAT* smoothMPC_v72 = smoothMPC_v + 216;
smoothMPC_FLOAT smoothMPC_re72[3];
smoothMPC_FLOAT smoothMPC_beta72[3];
smoothMPC_FLOAT smoothMPC_betacc72[3];
smoothMPC_FLOAT* smoothMPC_dvaff72 = smoothMPC_dv_aff + 216;
smoothMPC_FLOAT* smoothMPC_dvcc72 = smoothMPC_dv_cc + 216;
smoothMPC_FLOAT smoothMPC_V72[33];
smoothMPC_FLOAT smoothMPC_Yd72[6];
smoothMPC_FLOAT smoothMPC_Ld72[6];
smoothMPC_FLOAT smoothMPC_yy72[3];
smoothMPC_FLOAT smoothMPC_bmy72[3];
int smoothMPC_lbIdx72[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb72 = smoothMPC_l + 1152;
smoothMPC_FLOAT* smoothMPC_slb72 = smoothMPC_s + 1152;
smoothMPC_FLOAT* smoothMPC_llbbyslb72 = smoothMPC_lbys + 1152;
smoothMPC_FLOAT smoothMPC_rilb72[11];
smoothMPC_FLOAT* smoothMPC_dllbaff72 = smoothMPC_dl_aff + 1152;
smoothMPC_FLOAT* smoothMPC_dslbaff72 = smoothMPC_ds_aff + 1152;
smoothMPC_FLOAT* smoothMPC_dllbcc72 = smoothMPC_dl_cc + 1152;
smoothMPC_FLOAT* smoothMPC_dslbcc72 = smoothMPC_ds_cc + 1152;
smoothMPC_FLOAT* smoothMPC_ccrhsl72 = smoothMPC_ccrhs + 1152;
int smoothMPC_ubIdx72[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub72 = smoothMPC_l + 1163;
smoothMPC_FLOAT* smoothMPC_sub72 = smoothMPC_s + 1163;
smoothMPC_FLOAT* smoothMPC_lubbysub72 = smoothMPC_lbys + 1163;
smoothMPC_FLOAT smoothMPC_riub72[5];
smoothMPC_FLOAT* smoothMPC_dlubaff72 = smoothMPC_dl_aff + 1163;
smoothMPC_FLOAT* smoothMPC_dsubaff72 = smoothMPC_ds_aff + 1163;
smoothMPC_FLOAT* smoothMPC_dlubcc72 = smoothMPC_dl_cc + 1163;
smoothMPC_FLOAT* smoothMPC_dsubcc72 = smoothMPC_ds_cc + 1163;
smoothMPC_FLOAT* smoothMPC_ccrhsub72 = smoothMPC_ccrhs + 1163;
smoothMPC_FLOAT smoothMPC_Phi72[11];
smoothMPC_FLOAT smoothMPC_W72[11];
smoothMPC_FLOAT smoothMPC_Ysd72[9];
smoothMPC_FLOAT smoothMPC_Lsd72[9];
smoothMPC_FLOAT* smoothMPC_z73 = smoothMPC_z + 803;
smoothMPC_FLOAT* smoothMPC_dzaff73 = smoothMPC_dz_aff + 803;
smoothMPC_FLOAT* smoothMPC_dzcc73 = smoothMPC_dz_cc + 803;
smoothMPC_FLOAT* smoothMPC_rd73 = smoothMPC_rd + 803;
smoothMPC_FLOAT smoothMPC_Lbyrd73[11];
smoothMPC_FLOAT* smoothMPC_grad_cost73 = smoothMPC_grad_cost + 803;
smoothMPC_FLOAT* smoothMPC_grad_eq73 = smoothMPC_grad_eq + 803;
smoothMPC_FLOAT* smoothMPC_grad_ineq73 = smoothMPC_grad_ineq + 803;
smoothMPC_FLOAT smoothMPC_ctv73[11];
smoothMPC_FLOAT* smoothMPC_v73 = smoothMPC_v + 219;
smoothMPC_FLOAT smoothMPC_re73[3];
smoothMPC_FLOAT smoothMPC_beta73[3];
smoothMPC_FLOAT smoothMPC_betacc73[3];
smoothMPC_FLOAT* smoothMPC_dvaff73 = smoothMPC_dv_aff + 219;
smoothMPC_FLOAT* smoothMPC_dvcc73 = smoothMPC_dv_cc + 219;
smoothMPC_FLOAT smoothMPC_V73[33];
smoothMPC_FLOAT smoothMPC_Yd73[6];
smoothMPC_FLOAT smoothMPC_Ld73[6];
smoothMPC_FLOAT smoothMPC_yy73[3];
smoothMPC_FLOAT smoothMPC_bmy73[3];
int smoothMPC_lbIdx73[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb73 = smoothMPC_l + 1168;
smoothMPC_FLOAT* smoothMPC_slb73 = smoothMPC_s + 1168;
smoothMPC_FLOAT* smoothMPC_llbbyslb73 = smoothMPC_lbys + 1168;
smoothMPC_FLOAT smoothMPC_rilb73[11];
smoothMPC_FLOAT* smoothMPC_dllbaff73 = smoothMPC_dl_aff + 1168;
smoothMPC_FLOAT* smoothMPC_dslbaff73 = smoothMPC_ds_aff + 1168;
smoothMPC_FLOAT* smoothMPC_dllbcc73 = smoothMPC_dl_cc + 1168;
smoothMPC_FLOAT* smoothMPC_dslbcc73 = smoothMPC_ds_cc + 1168;
smoothMPC_FLOAT* smoothMPC_ccrhsl73 = smoothMPC_ccrhs + 1168;
int smoothMPC_ubIdx73[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub73 = smoothMPC_l + 1179;
smoothMPC_FLOAT* smoothMPC_sub73 = smoothMPC_s + 1179;
smoothMPC_FLOAT* smoothMPC_lubbysub73 = smoothMPC_lbys + 1179;
smoothMPC_FLOAT smoothMPC_riub73[5];
smoothMPC_FLOAT* smoothMPC_dlubaff73 = smoothMPC_dl_aff + 1179;
smoothMPC_FLOAT* smoothMPC_dsubaff73 = smoothMPC_ds_aff + 1179;
smoothMPC_FLOAT* smoothMPC_dlubcc73 = smoothMPC_dl_cc + 1179;
smoothMPC_FLOAT* smoothMPC_dsubcc73 = smoothMPC_ds_cc + 1179;
smoothMPC_FLOAT* smoothMPC_ccrhsub73 = smoothMPC_ccrhs + 1179;
smoothMPC_FLOAT smoothMPC_Phi73[11];
smoothMPC_FLOAT smoothMPC_W73[11];
smoothMPC_FLOAT smoothMPC_Ysd73[9];
smoothMPC_FLOAT smoothMPC_Lsd73[9];
smoothMPC_FLOAT* smoothMPC_z74 = smoothMPC_z + 814;
smoothMPC_FLOAT* smoothMPC_dzaff74 = smoothMPC_dz_aff + 814;
smoothMPC_FLOAT* smoothMPC_dzcc74 = smoothMPC_dz_cc + 814;
smoothMPC_FLOAT* smoothMPC_rd74 = smoothMPC_rd + 814;
smoothMPC_FLOAT smoothMPC_Lbyrd74[11];
smoothMPC_FLOAT* smoothMPC_grad_cost74 = smoothMPC_grad_cost + 814;
smoothMPC_FLOAT* smoothMPC_grad_eq74 = smoothMPC_grad_eq + 814;
smoothMPC_FLOAT* smoothMPC_grad_ineq74 = smoothMPC_grad_ineq + 814;
smoothMPC_FLOAT smoothMPC_ctv74[11];
smoothMPC_FLOAT* smoothMPC_v74 = smoothMPC_v + 222;
smoothMPC_FLOAT smoothMPC_re74[3];
smoothMPC_FLOAT smoothMPC_beta74[3];
smoothMPC_FLOAT smoothMPC_betacc74[3];
smoothMPC_FLOAT* smoothMPC_dvaff74 = smoothMPC_dv_aff + 222;
smoothMPC_FLOAT* smoothMPC_dvcc74 = smoothMPC_dv_cc + 222;
smoothMPC_FLOAT smoothMPC_V74[33];
smoothMPC_FLOAT smoothMPC_Yd74[6];
smoothMPC_FLOAT smoothMPC_Ld74[6];
smoothMPC_FLOAT smoothMPC_yy74[3];
smoothMPC_FLOAT smoothMPC_bmy74[3];
int smoothMPC_lbIdx74[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb74 = smoothMPC_l + 1184;
smoothMPC_FLOAT* smoothMPC_slb74 = smoothMPC_s + 1184;
smoothMPC_FLOAT* smoothMPC_llbbyslb74 = smoothMPC_lbys + 1184;
smoothMPC_FLOAT smoothMPC_rilb74[11];
smoothMPC_FLOAT* smoothMPC_dllbaff74 = smoothMPC_dl_aff + 1184;
smoothMPC_FLOAT* smoothMPC_dslbaff74 = smoothMPC_ds_aff + 1184;
smoothMPC_FLOAT* smoothMPC_dllbcc74 = smoothMPC_dl_cc + 1184;
smoothMPC_FLOAT* smoothMPC_dslbcc74 = smoothMPC_ds_cc + 1184;
smoothMPC_FLOAT* smoothMPC_ccrhsl74 = smoothMPC_ccrhs + 1184;
int smoothMPC_ubIdx74[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub74 = smoothMPC_l + 1195;
smoothMPC_FLOAT* smoothMPC_sub74 = smoothMPC_s + 1195;
smoothMPC_FLOAT* smoothMPC_lubbysub74 = smoothMPC_lbys + 1195;
smoothMPC_FLOAT smoothMPC_riub74[5];
smoothMPC_FLOAT* smoothMPC_dlubaff74 = smoothMPC_dl_aff + 1195;
smoothMPC_FLOAT* smoothMPC_dsubaff74 = smoothMPC_ds_aff + 1195;
smoothMPC_FLOAT* smoothMPC_dlubcc74 = smoothMPC_dl_cc + 1195;
smoothMPC_FLOAT* smoothMPC_dsubcc74 = smoothMPC_ds_cc + 1195;
smoothMPC_FLOAT* smoothMPC_ccrhsub74 = smoothMPC_ccrhs + 1195;
smoothMPC_FLOAT smoothMPC_Phi74[11];
smoothMPC_FLOAT smoothMPC_W74[11];
smoothMPC_FLOAT smoothMPC_Ysd74[9];
smoothMPC_FLOAT smoothMPC_Lsd74[9];
smoothMPC_FLOAT* smoothMPC_z75 = smoothMPC_z + 825;
smoothMPC_FLOAT* smoothMPC_dzaff75 = smoothMPC_dz_aff + 825;
smoothMPC_FLOAT* smoothMPC_dzcc75 = smoothMPC_dz_cc + 825;
smoothMPC_FLOAT* smoothMPC_rd75 = smoothMPC_rd + 825;
smoothMPC_FLOAT smoothMPC_Lbyrd75[11];
smoothMPC_FLOAT* smoothMPC_grad_cost75 = smoothMPC_grad_cost + 825;
smoothMPC_FLOAT* smoothMPC_grad_eq75 = smoothMPC_grad_eq + 825;
smoothMPC_FLOAT* smoothMPC_grad_ineq75 = smoothMPC_grad_ineq + 825;
smoothMPC_FLOAT smoothMPC_ctv75[11];
smoothMPC_FLOAT* smoothMPC_v75 = smoothMPC_v + 225;
smoothMPC_FLOAT smoothMPC_re75[3];
smoothMPC_FLOAT smoothMPC_beta75[3];
smoothMPC_FLOAT smoothMPC_betacc75[3];
smoothMPC_FLOAT* smoothMPC_dvaff75 = smoothMPC_dv_aff + 225;
smoothMPC_FLOAT* smoothMPC_dvcc75 = smoothMPC_dv_cc + 225;
smoothMPC_FLOAT smoothMPC_V75[33];
smoothMPC_FLOAT smoothMPC_Yd75[6];
smoothMPC_FLOAT smoothMPC_Ld75[6];
smoothMPC_FLOAT smoothMPC_yy75[3];
smoothMPC_FLOAT smoothMPC_bmy75[3];
int smoothMPC_lbIdx75[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb75 = smoothMPC_l + 1200;
smoothMPC_FLOAT* smoothMPC_slb75 = smoothMPC_s + 1200;
smoothMPC_FLOAT* smoothMPC_llbbyslb75 = smoothMPC_lbys + 1200;
smoothMPC_FLOAT smoothMPC_rilb75[11];
smoothMPC_FLOAT* smoothMPC_dllbaff75 = smoothMPC_dl_aff + 1200;
smoothMPC_FLOAT* smoothMPC_dslbaff75 = smoothMPC_ds_aff + 1200;
smoothMPC_FLOAT* smoothMPC_dllbcc75 = smoothMPC_dl_cc + 1200;
smoothMPC_FLOAT* smoothMPC_dslbcc75 = smoothMPC_ds_cc + 1200;
smoothMPC_FLOAT* smoothMPC_ccrhsl75 = smoothMPC_ccrhs + 1200;
int smoothMPC_ubIdx75[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub75 = smoothMPC_l + 1211;
smoothMPC_FLOAT* smoothMPC_sub75 = smoothMPC_s + 1211;
smoothMPC_FLOAT* smoothMPC_lubbysub75 = smoothMPC_lbys + 1211;
smoothMPC_FLOAT smoothMPC_riub75[5];
smoothMPC_FLOAT* smoothMPC_dlubaff75 = smoothMPC_dl_aff + 1211;
smoothMPC_FLOAT* smoothMPC_dsubaff75 = smoothMPC_ds_aff + 1211;
smoothMPC_FLOAT* smoothMPC_dlubcc75 = smoothMPC_dl_cc + 1211;
smoothMPC_FLOAT* smoothMPC_dsubcc75 = smoothMPC_ds_cc + 1211;
smoothMPC_FLOAT* smoothMPC_ccrhsub75 = smoothMPC_ccrhs + 1211;
smoothMPC_FLOAT smoothMPC_Phi75[11];
smoothMPC_FLOAT smoothMPC_W75[11];
smoothMPC_FLOAT smoothMPC_Ysd75[9];
smoothMPC_FLOAT smoothMPC_Lsd75[9];
smoothMPC_FLOAT* smoothMPC_z76 = smoothMPC_z + 836;
smoothMPC_FLOAT* smoothMPC_dzaff76 = smoothMPC_dz_aff + 836;
smoothMPC_FLOAT* smoothMPC_dzcc76 = smoothMPC_dz_cc + 836;
smoothMPC_FLOAT* smoothMPC_rd76 = smoothMPC_rd + 836;
smoothMPC_FLOAT smoothMPC_Lbyrd76[11];
smoothMPC_FLOAT* smoothMPC_grad_cost76 = smoothMPC_grad_cost + 836;
smoothMPC_FLOAT* smoothMPC_grad_eq76 = smoothMPC_grad_eq + 836;
smoothMPC_FLOAT* smoothMPC_grad_ineq76 = smoothMPC_grad_ineq + 836;
smoothMPC_FLOAT smoothMPC_ctv76[11];
smoothMPC_FLOAT* smoothMPC_v76 = smoothMPC_v + 228;
smoothMPC_FLOAT smoothMPC_re76[3];
smoothMPC_FLOAT smoothMPC_beta76[3];
smoothMPC_FLOAT smoothMPC_betacc76[3];
smoothMPC_FLOAT* smoothMPC_dvaff76 = smoothMPC_dv_aff + 228;
smoothMPC_FLOAT* smoothMPC_dvcc76 = smoothMPC_dv_cc + 228;
smoothMPC_FLOAT smoothMPC_V76[33];
smoothMPC_FLOAT smoothMPC_Yd76[6];
smoothMPC_FLOAT smoothMPC_Ld76[6];
smoothMPC_FLOAT smoothMPC_yy76[3];
smoothMPC_FLOAT smoothMPC_bmy76[3];
int smoothMPC_lbIdx76[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb76 = smoothMPC_l + 1216;
smoothMPC_FLOAT* smoothMPC_slb76 = smoothMPC_s + 1216;
smoothMPC_FLOAT* smoothMPC_llbbyslb76 = smoothMPC_lbys + 1216;
smoothMPC_FLOAT smoothMPC_rilb76[11];
smoothMPC_FLOAT* smoothMPC_dllbaff76 = smoothMPC_dl_aff + 1216;
smoothMPC_FLOAT* smoothMPC_dslbaff76 = smoothMPC_ds_aff + 1216;
smoothMPC_FLOAT* smoothMPC_dllbcc76 = smoothMPC_dl_cc + 1216;
smoothMPC_FLOAT* smoothMPC_dslbcc76 = smoothMPC_ds_cc + 1216;
smoothMPC_FLOAT* smoothMPC_ccrhsl76 = smoothMPC_ccrhs + 1216;
int smoothMPC_ubIdx76[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub76 = smoothMPC_l + 1227;
smoothMPC_FLOAT* smoothMPC_sub76 = smoothMPC_s + 1227;
smoothMPC_FLOAT* smoothMPC_lubbysub76 = smoothMPC_lbys + 1227;
smoothMPC_FLOAT smoothMPC_riub76[5];
smoothMPC_FLOAT* smoothMPC_dlubaff76 = smoothMPC_dl_aff + 1227;
smoothMPC_FLOAT* smoothMPC_dsubaff76 = smoothMPC_ds_aff + 1227;
smoothMPC_FLOAT* smoothMPC_dlubcc76 = smoothMPC_dl_cc + 1227;
smoothMPC_FLOAT* smoothMPC_dsubcc76 = smoothMPC_ds_cc + 1227;
smoothMPC_FLOAT* smoothMPC_ccrhsub76 = smoothMPC_ccrhs + 1227;
smoothMPC_FLOAT smoothMPC_Phi76[11];
smoothMPC_FLOAT smoothMPC_W76[11];
smoothMPC_FLOAT smoothMPC_Ysd76[9];
smoothMPC_FLOAT smoothMPC_Lsd76[9];
smoothMPC_FLOAT* smoothMPC_z77 = smoothMPC_z + 847;
smoothMPC_FLOAT* smoothMPC_dzaff77 = smoothMPC_dz_aff + 847;
smoothMPC_FLOAT* smoothMPC_dzcc77 = smoothMPC_dz_cc + 847;
smoothMPC_FLOAT* smoothMPC_rd77 = smoothMPC_rd + 847;
smoothMPC_FLOAT smoothMPC_Lbyrd77[11];
smoothMPC_FLOAT* smoothMPC_grad_cost77 = smoothMPC_grad_cost + 847;
smoothMPC_FLOAT* smoothMPC_grad_eq77 = smoothMPC_grad_eq + 847;
smoothMPC_FLOAT* smoothMPC_grad_ineq77 = smoothMPC_grad_ineq + 847;
smoothMPC_FLOAT smoothMPC_ctv77[11];
smoothMPC_FLOAT* smoothMPC_v77 = smoothMPC_v + 231;
smoothMPC_FLOAT smoothMPC_re77[3];
smoothMPC_FLOAT smoothMPC_beta77[3];
smoothMPC_FLOAT smoothMPC_betacc77[3];
smoothMPC_FLOAT* smoothMPC_dvaff77 = smoothMPC_dv_aff + 231;
smoothMPC_FLOAT* smoothMPC_dvcc77 = smoothMPC_dv_cc + 231;
smoothMPC_FLOAT smoothMPC_V77[33];
smoothMPC_FLOAT smoothMPC_Yd77[6];
smoothMPC_FLOAT smoothMPC_Ld77[6];
smoothMPC_FLOAT smoothMPC_yy77[3];
smoothMPC_FLOAT smoothMPC_bmy77[3];
int smoothMPC_lbIdx77[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb77 = smoothMPC_l + 1232;
smoothMPC_FLOAT* smoothMPC_slb77 = smoothMPC_s + 1232;
smoothMPC_FLOAT* smoothMPC_llbbyslb77 = smoothMPC_lbys + 1232;
smoothMPC_FLOAT smoothMPC_rilb77[11];
smoothMPC_FLOAT* smoothMPC_dllbaff77 = smoothMPC_dl_aff + 1232;
smoothMPC_FLOAT* smoothMPC_dslbaff77 = smoothMPC_ds_aff + 1232;
smoothMPC_FLOAT* smoothMPC_dllbcc77 = smoothMPC_dl_cc + 1232;
smoothMPC_FLOAT* smoothMPC_dslbcc77 = smoothMPC_ds_cc + 1232;
smoothMPC_FLOAT* smoothMPC_ccrhsl77 = smoothMPC_ccrhs + 1232;
int smoothMPC_ubIdx77[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub77 = smoothMPC_l + 1243;
smoothMPC_FLOAT* smoothMPC_sub77 = smoothMPC_s + 1243;
smoothMPC_FLOAT* smoothMPC_lubbysub77 = smoothMPC_lbys + 1243;
smoothMPC_FLOAT smoothMPC_riub77[5];
smoothMPC_FLOAT* smoothMPC_dlubaff77 = smoothMPC_dl_aff + 1243;
smoothMPC_FLOAT* smoothMPC_dsubaff77 = smoothMPC_ds_aff + 1243;
smoothMPC_FLOAT* smoothMPC_dlubcc77 = smoothMPC_dl_cc + 1243;
smoothMPC_FLOAT* smoothMPC_dsubcc77 = smoothMPC_ds_cc + 1243;
smoothMPC_FLOAT* smoothMPC_ccrhsub77 = smoothMPC_ccrhs + 1243;
smoothMPC_FLOAT smoothMPC_Phi77[11];
smoothMPC_FLOAT smoothMPC_W77[11];
smoothMPC_FLOAT smoothMPC_Ysd77[9];
smoothMPC_FLOAT smoothMPC_Lsd77[9];
smoothMPC_FLOAT* smoothMPC_z78 = smoothMPC_z + 858;
smoothMPC_FLOAT* smoothMPC_dzaff78 = smoothMPC_dz_aff + 858;
smoothMPC_FLOAT* smoothMPC_dzcc78 = smoothMPC_dz_cc + 858;
smoothMPC_FLOAT* smoothMPC_rd78 = smoothMPC_rd + 858;
smoothMPC_FLOAT smoothMPC_Lbyrd78[11];
smoothMPC_FLOAT* smoothMPC_grad_cost78 = smoothMPC_grad_cost + 858;
smoothMPC_FLOAT* smoothMPC_grad_eq78 = smoothMPC_grad_eq + 858;
smoothMPC_FLOAT* smoothMPC_grad_ineq78 = smoothMPC_grad_ineq + 858;
smoothMPC_FLOAT smoothMPC_ctv78[11];
smoothMPC_FLOAT* smoothMPC_v78 = smoothMPC_v + 234;
smoothMPC_FLOAT smoothMPC_re78[3];
smoothMPC_FLOAT smoothMPC_beta78[3];
smoothMPC_FLOAT smoothMPC_betacc78[3];
smoothMPC_FLOAT* smoothMPC_dvaff78 = smoothMPC_dv_aff + 234;
smoothMPC_FLOAT* smoothMPC_dvcc78 = smoothMPC_dv_cc + 234;
smoothMPC_FLOAT smoothMPC_V78[33];
smoothMPC_FLOAT smoothMPC_Yd78[6];
smoothMPC_FLOAT smoothMPC_Ld78[6];
smoothMPC_FLOAT smoothMPC_yy78[3];
smoothMPC_FLOAT smoothMPC_bmy78[3];
int smoothMPC_lbIdx78[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb78 = smoothMPC_l + 1248;
smoothMPC_FLOAT* smoothMPC_slb78 = smoothMPC_s + 1248;
smoothMPC_FLOAT* smoothMPC_llbbyslb78 = smoothMPC_lbys + 1248;
smoothMPC_FLOAT smoothMPC_rilb78[11];
smoothMPC_FLOAT* smoothMPC_dllbaff78 = smoothMPC_dl_aff + 1248;
smoothMPC_FLOAT* smoothMPC_dslbaff78 = smoothMPC_ds_aff + 1248;
smoothMPC_FLOAT* smoothMPC_dllbcc78 = smoothMPC_dl_cc + 1248;
smoothMPC_FLOAT* smoothMPC_dslbcc78 = smoothMPC_ds_cc + 1248;
smoothMPC_FLOAT* smoothMPC_ccrhsl78 = smoothMPC_ccrhs + 1248;
int smoothMPC_ubIdx78[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub78 = smoothMPC_l + 1259;
smoothMPC_FLOAT* smoothMPC_sub78 = smoothMPC_s + 1259;
smoothMPC_FLOAT* smoothMPC_lubbysub78 = smoothMPC_lbys + 1259;
smoothMPC_FLOAT smoothMPC_riub78[5];
smoothMPC_FLOAT* smoothMPC_dlubaff78 = smoothMPC_dl_aff + 1259;
smoothMPC_FLOAT* smoothMPC_dsubaff78 = smoothMPC_ds_aff + 1259;
smoothMPC_FLOAT* smoothMPC_dlubcc78 = smoothMPC_dl_cc + 1259;
smoothMPC_FLOAT* smoothMPC_dsubcc78 = smoothMPC_ds_cc + 1259;
smoothMPC_FLOAT* smoothMPC_ccrhsub78 = smoothMPC_ccrhs + 1259;
smoothMPC_FLOAT smoothMPC_Phi78[11];
smoothMPC_FLOAT smoothMPC_W78[11];
smoothMPC_FLOAT smoothMPC_Ysd78[9];
smoothMPC_FLOAT smoothMPC_Lsd78[9];
smoothMPC_FLOAT* smoothMPC_z79 = smoothMPC_z + 869;
smoothMPC_FLOAT* smoothMPC_dzaff79 = smoothMPC_dz_aff + 869;
smoothMPC_FLOAT* smoothMPC_dzcc79 = smoothMPC_dz_cc + 869;
smoothMPC_FLOAT* smoothMPC_rd79 = smoothMPC_rd + 869;
smoothMPC_FLOAT smoothMPC_Lbyrd79[11];
smoothMPC_FLOAT* smoothMPC_grad_cost79 = smoothMPC_grad_cost + 869;
smoothMPC_FLOAT* smoothMPC_grad_eq79 = smoothMPC_grad_eq + 869;
smoothMPC_FLOAT* smoothMPC_grad_ineq79 = smoothMPC_grad_ineq + 869;
smoothMPC_FLOAT smoothMPC_ctv79[11];
smoothMPC_FLOAT* smoothMPC_v79 = smoothMPC_v + 237;
smoothMPC_FLOAT smoothMPC_re79[3];
smoothMPC_FLOAT smoothMPC_beta79[3];
smoothMPC_FLOAT smoothMPC_betacc79[3];
smoothMPC_FLOAT* smoothMPC_dvaff79 = smoothMPC_dv_aff + 237;
smoothMPC_FLOAT* smoothMPC_dvcc79 = smoothMPC_dv_cc + 237;
smoothMPC_FLOAT smoothMPC_V79[33];
smoothMPC_FLOAT smoothMPC_Yd79[6];
smoothMPC_FLOAT smoothMPC_Ld79[6];
smoothMPC_FLOAT smoothMPC_yy79[3];
smoothMPC_FLOAT smoothMPC_bmy79[3];
int smoothMPC_lbIdx79[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb79 = smoothMPC_l + 1264;
smoothMPC_FLOAT* smoothMPC_slb79 = smoothMPC_s + 1264;
smoothMPC_FLOAT* smoothMPC_llbbyslb79 = smoothMPC_lbys + 1264;
smoothMPC_FLOAT smoothMPC_rilb79[11];
smoothMPC_FLOAT* smoothMPC_dllbaff79 = smoothMPC_dl_aff + 1264;
smoothMPC_FLOAT* smoothMPC_dslbaff79 = smoothMPC_ds_aff + 1264;
smoothMPC_FLOAT* smoothMPC_dllbcc79 = smoothMPC_dl_cc + 1264;
smoothMPC_FLOAT* smoothMPC_dslbcc79 = smoothMPC_ds_cc + 1264;
smoothMPC_FLOAT* smoothMPC_ccrhsl79 = smoothMPC_ccrhs + 1264;
int smoothMPC_ubIdx79[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub79 = smoothMPC_l + 1275;
smoothMPC_FLOAT* smoothMPC_sub79 = smoothMPC_s + 1275;
smoothMPC_FLOAT* smoothMPC_lubbysub79 = smoothMPC_lbys + 1275;
smoothMPC_FLOAT smoothMPC_riub79[5];
smoothMPC_FLOAT* smoothMPC_dlubaff79 = smoothMPC_dl_aff + 1275;
smoothMPC_FLOAT* smoothMPC_dsubaff79 = smoothMPC_ds_aff + 1275;
smoothMPC_FLOAT* smoothMPC_dlubcc79 = smoothMPC_dl_cc + 1275;
smoothMPC_FLOAT* smoothMPC_dsubcc79 = smoothMPC_ds_cc + 1275;
smoothMPC_FLOAT* smoothMPC_ccrhsub79 = smoothMPC_ccrhs + 1275;
smoothMPC_FLOAT smoothMPC_Phi79[11];
smoothMPC_FLOAT smoothMPC_W79[11];
smoothMPC_FLOAT smoothMPC_Ysd79[9];
smoothMPC_FLOAT smoothMPC_Lsd79[9];
smoothMPC_FLOAT* smoothMPC_z80 = smoothMPC_z + 880;
smoothMPC_FLOAT* smoothMPC_dzaff80 = smoothMPC_dz_aff + 880;
smoothMPC_FLOAT* smoothMPC_dzcc80 = smoothMPC_dz_cc + 880;
smoothMPC_FLOAT* smoothMPC_rd80 = smoothMPC_rd + 880;
smoothMPC_FLOAT smoothMPC_Lbyrd80[11];
smoothMPC_FLOAT* smoothMPC_grad_cost80 = smoothMPC_grad_cost + 880;
smoothMPC_FLOAT* smoothMPC_grad_eq80 = smoothMPC_grad_eq + 880;
smoothMPC_FLOAT* smoothMPC_grad_ineq80 = smoothMPC_grad_ineq + 880;
smoothMPC_FLOAT smoothMPC_ctv80[11];
smoothMPC_FLOAT* smoothMPC_v80 = smoothMPC_v + 240;
smoothMPC_FLOAT smoothMPC_re80[3];
smoothMPC_FLOAT smoothMPC_beta80[3];
smoothMPC_FLOAT smoothMPC_betacc80[3];
smoothMPC_FLOAT* smoothMPC_dvaff80 = smoothMPC_dv_aff + 240;
smoothMPC_FLOAT* smoothMPC_dvcc80 = smoothMPC_dv_cc + 240;
smoothMPC_FLOAT smoothMPC_V80[33];
smoothMPC_FLOAT smoothMPC_Yd80[6];
smoothMPC_FLOAT smoothMPC_Ld80[6];
smoothMPC_FLOAT smoothMPC_yy80[3];
smoothMPC_FLOAT smoothMPC_bmy80[3];
int smoothMPC_lbIdx80[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb80 = smoothMPC_l + 1280;
smoothMPC_FLOAT* smoothMPC_slb80 = smoothMPC_s + 1280;
smoothMPC_FLOAT* smoothMPC_llbbyslb80 = smoothMPC_lbys + 1280;
smoothMPC_FLOAT smoothMPC_rilb80[11];
smoothMPC_FLOAT* smoothMPC_dllbaff80 = smoothMPC_dl_aff + 1280;
smoothMPC_FLOAT* smoothMPC_dslbaff80 = smoothMPC_ds_aff + 1280;
smoothMPC_FLOAT* smoothMPC_dllbcc80 = smoothMPC_dl_cc + 1280;
smoothMPC_FLOAT* smoothMPC_dslbcc80 = smoothMPC_ds_cc + 1280;
smoothMPC_FLOAT* smoothMPC_ccrhsl80 = smoothMPC_ccrhs + 1280;
int smoothMPC_ubIdx80[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub80 = smoothMPC_l + 1291;
smoothMPC_FLOAT* smoothMPC_sub80 = smoothMPC_s + 1291;
smoothMPC_FLOAT* smoothMPC_lubbysub80 = smoothMPC_lbys + 1291;
smoothMPC_FLOAT smoothMPC_riub80[5];
smoothMPC_FLOAT* smoothMPC_dlubaff80 = smoothMPC_dl_aff + 1291;
smoothMPC_FLOAT* smoothMPC_dsubaff80 = smoothMPC_ds_aff + 1291;
smoothMPC_FLOAT* smoothMPC_dlubcc80 = smoothMPC_dl_cc + 1291;
smoothMPC_FLOAT* smoothMPC_dsubcc80 = smoothMPC_ds_cc + 1291;
smoothMPC_FLOAT* smoothMPC_ccrhsub80 = smoothMPC_ccrhs + 1291;
smoothMPC_FLOAT smoothMPC_Phi80[11];
smoothMPC_FLOAT smoothMPC_W80[11];
smoothMPC_FLOAT smoothMPC_Ysd80[9];
smoothMPC_FLOAT smoothMPC_Lsd80[9];
smoothMPC_FLOAT* smoothMPC_z81 = smoothMPC_z + 891;
smoothMPC_FLOAT* smoothMPC_dzaff81 = smoothMPC_dz_aff + 891;
smoothMPC_FLOAT* smoothMPC_dzcc81 = smoothMPC_dz_cc + 891;
smoothMPC_FLOAT* smoothMPC_rd81 = smoothMPC_rd + 891;
smoothMPC_FLOAT smoothMPC_Lbyrd81[11];
smoothMPC_FLOAT* smoothMPC_grad_cost81 = smoothMPC_grad_cost + 891;
smoothMPC_FLOAT* smoothMPC_grad_eq81 = smoothMPC_grad_eq + 891;
smoothMPC_FLOAT* smoothMPC_grad_ineq81 = smoothMPC_grad_ineq + 891;
smoothMPC_FLOAT smoothMPC_ctv81[11];
smoothMPC_FLOAT* smoothMPC_v81 = smoothMPC_v + 243;
smoothMPC_FLOAT smoothMPC_re81[3];
smoothMPC_FLOAT smoothMPC_beta81[3];
smoothMPC_FLOAT smoothMPC_betacc81[3];
smoothMPC_FLOAT* smoothMPC_dvaff81 = smoothMPC_dv_aff + 243;
smoothMPC_FLOAT* smoothMPC_dvcc81 = smoothMPC_dv_cc + 243;
smoothMPC_FLOAT smoothMPC_V81[33];
smoothMPC_FLOAT smoothMPC_Yd81[6];
smoothMPC_FLOAT smoothMPC_Ld81[6];
smoothMPC_FLOAT smoothMPC_yy81[3];
smoothMPC_FLOAT smoothMPC_bmy81[3];
int smoothMPC_lbIdx81[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb81 = smoothMPC_l + 1296;
smoothMPC_FLOAT* smoothMPC_slb81 = smoothMPC_s + 1296;
smoothMPC_FLOAT* smoothMPC_llbbyslb81 = smoothMPC_lbys + 1296;
smoothMPC_FLOAT smoothMPC_rilb81[11];
smoothMPC_FLOAT* smoothMPC_dllbaff81 = smoothMPC_dl_aff + 1296;
smoothMPC_FLOAT* smoothMPC_dslbaff81 = smoothMPC_ds_aff + 1296;
smoothMPC_FLOAT* smoothMPC_dllbcc81 = smoothMPC_dl_cc + 1296;
smoothMPC_FLOAT* smoothMPC_dslbcc81 = smoothMPC_ds_cc + 1296;
smoothMPC_FLOAT* smoothMPC_ccrhsl81 = smoothMPC_ccrhs + 1296;
int smoothMPC_ubIdx81[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub81 = smoothMPC_l + 1307;
smoothMPC_FLOAT* smoothMPC_sub81 = smoothMPC_s + 1307;
smoothMPC_FLOAT* smoothMPC_lubbysub81 = smoothMPC_lbys + 1307;
smoothMPC_FLOAT smoothMPC_riub81[5];
smoothMPC_FLOAT* smoothMPC_dlubaff81 = smoothMPC_dl_aff + 1307;
smoothMPC_FLOAT* smoothMPC_dsubaff81 = smoothMPC_ds_aff + 1307;
smoothMPC_FLOAT* smoothMPC_dlubcc81 = smoothMPC_dl_cc + 1307;
smoothMPC_FLOAT* smoothMPC_dsubcc81 = smoothMPC_ds_cc + 1307;
smoothMPC_FLOAT* smoothMPC_ccrhsub81 = smoothMPC_ccrhs + 1307;
smoothMPC_FLOAT smoothMPC_Phi81[11];
smoothMPC_FLOAT smoothMPC_W81[11];
smoothMPC_FLOAT smoothMPC_Ysd81[9];
smoothMPC_FLOAT smoothMPC_Lsd81[9];
smoothMPC_FLOAT* smoothMPC_z82 = smoothMPC_z + 902;
smoothMPC_FLOAT* smoothMPC_dzaff82 = smoothMPC_dz_aff + 902;
smoothMPC_FLOAT* smoothMPC_dzcc82 = smoothMPC_dz_cc + 902;
smoothMPC_FLOAT* smoothMPC_rd82 = smoothMPC_rd + 902;
smoothMPC_FLOAT smoothMPC_Lbyrd82[11];
smoothMPC_FLOAT* smoothMPC_grad_cost82 = smoothMPC_grad_cost + 902;
smoothMPC_FLOAT* smoothMPC_grad_eq82 = smoothMPC_grad_eq + 902;
smoothMPC_FLOAT* smoothMPC_grad_ineq82 = smoothMPC_grad_ineq + 902;
smoothMPC_FLOAT smoothMPC_ctv82[11];
smoothMPC_FLOAT* smoothMPC_v82 = smoothMPC_v + 246;
smoothMPC_FLOAT smoothMPC_re82[3];
smoothMPC_FLOAT smoothMPC_beta82[3];
smoothMPC_FLOAT smoothMPC_betacc82[3];
smoothMPC_FLOAT* smoothMPC_dvaff82 = smoothMPC_dv_aff + 246;
smoothMPC_FLOAT* smoothMPC_dvcc82 = smoothMPC_dv_cc + 246;
smoothMPC_FLOAT smoothMPC_V82[33];
smoothMPC_FLOAT smoothMPC_Yd82[6];
smoothMPC_FLOAT smoothMPC_Ld82[6];
smoothMPC_FLOAT smoothMPC_yy82[3];
smoothMPC_FLOAT smoothMPC_bmy82[3];
int smoothMPC_lbIdx82[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb82 = smoothMPC_l + 1312;
smoothMPC_FLOAT* smoothMPC_slb82 = smoothMPC_s + 1312;
smoothMPC_FLOAT* smoothMPC_llbbyslb82 = smoothMPC_lbys + 1312;
smoothMPC_FLOAT smoothMPC_rilb82[11];
smoothMPC_FLOAT* smoothMPC_dllbaff82 = smoothMPC_dl_aff + 1312;
smoothMPC_FLOAT* smoothMPC_dslbaff82 = smoothMPC_ds_aff + 1312;
smoothMPC_FLOAT* smoothMPC_dllbcc82 = smoothMPC_dl_cc + 1312;
smoothMPC_FLOAT* smoothMPC_dslbcc82 = smoothMPC_ds_cc + 1312;
smoothMPC_FLOAT* smoothMPC_ccrhsl82 = smoothMPC_ccrhs + 1312;
int smoothMPC_ubIdx82[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub82 = smoothMPC_l + 1323;
smoothMPC_FLOAT* smoothMPC_sub82 = smoothMPC_s + 1323;
smoothMPC_FLOAT* smoothMPC_lubbysub82 = smoothMPC_lbys + 1323;
smoothMPC_FLOAT smoothMPC_riub82[5];
smoothMPC_FLOAT* smoothMPC_dlubaff82 = smoothMPC_dl_aff + 1323;
smoothMPC_FLOAT* smoothMPC_dsubaff82 = smoothMPC_ds_aff + 1323;
smoothMPC_FLOAT* smoothMPC_dlubcc82 = smoothMPC_dl_cc + 1323;
smoothMPC_FLOAT* smoothMPC_dsubcc82 = smoothMPC_ds_cc + 1323;
smoothMPC_FLOAT* smoothMPC_ccrhsub82 = smoothMPC_ccrhs + 1323;
smoothMPC_FLOAT smoothMPC_Phi82[11];
smoothMPC_FLOAT smoothMPC_W82[11];
smoothMPC_FLOAT smoothMPC_Ysd82[9];
smoothMPC_FLOAT smoothMPC_Lsd82[9];
smoothMPC_FLOAT* smoothMPC_z83 = smoothMPC_z + 913;
smoothMPC_FLOAT* smoothMPC_dzaff83 = smoothMPC_dz_aff + 913;
smoothMPC_FLOAT* smoothMPC_dzcc83 = smoothMPC_dz_cc + 913;
smoothMPC_FLOAT* smoothMPC_rd83 = smoothMPC_rd + 913;
smoothMPC_FLOAT smoothMPC_Lbyrd83[11];
smoothMPC_FLOAT* smoothMPC_grad_cost83 = smoothMPC_grad_cost + 913;
smoothMPC_FLOAT* smoothMPC_grad_eq83 = smoothMPC_grad_eq + 913;
smoothMPC_FLOAT* smoothMPC_grad_ineq83 = smoothMPC_grad_ineq + 913;
smoothMPC_FLOAT smoothMPC_ctv83[11];
smoothMPC_FLOAT* smoothMPC_v83 = smoothMPC_v + 249;
smoothMPC_FLOAT smoothMPC_re83[3];
smoothMPC_FLOAT smoothMPC_beta83[3];
smoothMPC_FLOAT smoothMPC_betacc83[3];
smoothMPC_FLOAT* smoothMPC_dvaff83 = smoothMPC_dv_aff + 249;
smoothMPC_FLOAT* smoothMPC_dvcc83 = smoothMPC_dv_cc + 249;
smoothMPC_FLOAT smoothMPC_V83[33];
smoothMPC_FLOAT smoothMPC_Yd83[6];
smoothMPC_FLOAT smoothMPC_Ld83[6];
smoothMPC_FLOAT smoothMPC_yy83[3];
smoothMPC_FLOAT smoothMPC_bmy83[3];
int smoothMPC_lbIdx83[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb83 = smoothMPC_l + 1328;
smoothMPC_FLOAT* smoothMPC_slb83 = smoothMPC_s + 1328;
smoothMPC_FLOAT* smoothMPC_llbbyslb83 = smoothMPC_lbys + 1328;
smoothMPC_FLOAT smoothMPC_rilb83[11];
smoothMPC_FLOAT* smoothMPC_dllbaff83 = smoothMPC_dl_aff + 1328;
smoothMPC_FLOAT* smoothMPC_dslbaff83 = smoothMPC_ds_aff + 1328;
smoothMPC_FLOAT* smoothMPC_dllbcc83 = smoothMPC_dl_cc + 1328;
smoothMPC_FLOAT* smoothMPC_dslbcc83 = smoothMPC_ds_cc + 1328;
smoothMPC_FLOAT* smoothMPC_ccrhsl83 = smoothMPC_ccrhs + 1328;
int smoothMPC_ubIdx83[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub83 = smoothMPC_l + 1339;
smoothMPC_FLOAT* smoothMPC_sub83 = smoothMPC_s + 1339;
smoothMPC_FLOAT* smoothMPC_lubbysub83 = smoothMPC_lbys + 1339;
smoothMPC_FLOAT smoothMPC_riub83[5];
smoothMPC_FLOAT* smoothMPC_dlubaff83 = smoothMPC_dl_aff + 1339;
smoothMPC_FLOAT* smoothMPC_dsubaff83 = smoothMPC_ds_aff + 1339;
smoothMPC_FLOAT* smoothMPC_dlubcc83 = smoothMPC_dl_cc + 1339;
smoothMPC_FLOAT* smoothMPC_dsubcc83 = smoothMPC_ds_cc + 1339;
smoothMPC_FLOAT* smoothMPC_ccrhsub83 = smoothMPC_ccrhs + 1339;
smoothMPC_FLOAT smoothMPC_Phi83[11];
smoothMPC_FLOAT smoothMPC_W83[11];
smoothMPC_FLOAT smoothMPC_Ysd83[9];
smoothMPC_FLOAT smoothMPC_Lsd83[9];
smoothMPC_FLOAT* smoothMPC_z84 = smoothMPC_z + 924;
smoothMPC_FLOAT* smoothMPC_dzaff84 = smoothMPC_dz_aff + 924;
smoothMPC_FLOAT* smoothMPC_dzcc84 = smoothMPC_dz_cc + 924;
smoothMPC_FLOAT* smoothMPC_rd84 = smoothMPC_rd + 924;
smoothMPC_FLOAT smoothMPC_Lbyrd84[11];
smoothMPC_FLOAT* smoothMPC_grad_cost84 = smoothMPC_grad_cost + 924;
smoothMPC_FLOAT* smoothMPC_grad_eq84 = smoothMPC_grad_eq + 924;
smoothMPC_FLOAT* smoothMPC_grad_ineq84 = smoothMPC_grad_ineq + 924;
smoothMPC_FLOAT smoothMPC_ctv84[11];
smoothMPC_FLOAT* smoothMPC_v84 = smoothMPC_v + 252;
smoothMPC_FLOAT smoothMPC_re84[3];
smoothMPC_FLOAT smoothMPC_beta84[3];
smoothMPC_FLOAT smoothMPC_betacc84[3];
smoothMPC_FLOAT* smoothMPC_dvaff84 = smoothMPC_dv_aff + 252;
smoothMPC_FLOAT* smoothMPC_dvcc84 = smoothMPC_dv_cc + 252;
smoothMPC_FLOAT smoothMPC_V84[33];
smoothMPC_FLOAT smoothMPC_Yd84[6];
smoothMPC_FLOAT smoothMPC_Ld84[6];
smoothMPC_FLOAT smoothMPC_yy84[3];
smoothMPC_FLOAT smoothMPC_bmy84[3];
int smoothMPC_lbIdx84[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb84 = smoothMPC_l + 1344;
smoothMPC_FLOAT* smoothMPC_slb84 = smoothMPC_s + 1344;
smoothMPC_FLOAT* smoothMPC_llbbyslb84 = smoothMPC_lbys + 1344;
smoothMPC_FLOAT smoothMPC_rilb84[11];
smoothMPC_FLOAT* smoothMPC_dllbaff84 = smoothMPC_dl_aff + 1344;
smoothMPC_FLOAT* smoothMPC_dslbaff84 = smoothMPC_ds_aff + 1344;
smoothMPC_FLOAT* smoothMPC_dllbcc84 = smoothMPC_dl_cc + 1344;
smoothMPC_FLOAT* smoothMPC_dslbcc84 = smoothMPC_ds_cc + 1344;
smoothMPC_FLOAT* smoothMPC_ccrhsl84 = smoothMPC_ccrhs + 1344;
int smoothMPC_ubIdx84[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub84 = smoothMPC_l + 1355;
smoothMPC_FLOAT* smoothMPC_sub84 = smoothMPC_s + 1355;
smoothMPC_FLOAT* smoothMPC_lubbysub84 = smoothMPC_lbys + 1355;
smoothMPC_FLOAT smoothMPC_riub84[5];
smoothMPC_FLOAT* smoothMPC_dlubaff84 = smoothMPC_dl_aff + 1355;
smoothMPC_FLOAT* smoothMPC_dsubaff84 = smoothMPC_ds_aff + 1355;
smoothMPC_FLOAT* smoothMPC_dlubcc84 = smoothMPC_dl_cc + 1355;
smoothMPC_FLOAT* smoothMPC_dsubcc84 = smoothMPC_ds_cc + 1355;
smoothMPC_FLOAT* smoothMPC_ccrhsub84 = smoothMPC_ccrhs + 1355;
smoothMPC_FLOAT smoothMPC_Phi84[11];
smoothMPC_FLOAT smoothMPC_W84[11];
smoothMPC_FLOAT smoothMPC_Ysd84[9];
smoothMPC_FLOAT smoothMPC_Lsd84[9];
smoothMPC_FLOAT* smoothMPC_z85 = smoothMPC_z + 935;
smoothMPC_FLOAT* smoothMPC_dzaff85 = smoothMPC_dz_aff + 935;
smoothMPC_FLOAT* smoothMPC_dzcc85 = smoothMPC_dz_cc + 935;
smoothMPC_FLOAT* smoothMPC_rd85 = smoothMPC_rd + 935;
smoothMPC_FLOAT smoothMPC_Lbyrd85[11];
smoothMPC_FLOAT* smoothMPC_grad_cost85 = smoothMPC_grad_cost + 935;
smoothMPC_FLOAT* smoothMPC_grad_eq85 = smoothMPC_grad_eq + 935;
smoothMPC_FLOAT* smoothMPC_grad_ineq85 = smoothMPC_grad_ineq + 935;
smoothMPC_FLOAT smoothMPC_ctv85[11];
smoothMPC_FLOAT* smoothMPC_v85 = smoothMPC_v + 255;
smoothMPC_FLOAT smoothMPC_re85[3];
smoothMPC_FLOAT smoothMPC_beta85[3];
smoothMPC_FLOAT smoothMPC_betacc85[3];
smoothMPC_FLOAT* smoothMPC_dvaff85 = smoothMPC_dv_aff + 255;
smoothMPC_FLOAT* smoothMPC_dvcc85 = smoothMPC_dv_cc + 255;
smoothMPC_FLOAT smoothMPC_V85[33];
smoothMPC_FLOAT smoothMPC_Yd85[6];
smoothMPC_FLOAT smoothMPC_Ld85[6];
smoothMPC_FLOAT smoothMPC_yy85[3];
smoothMPC_FLOAT smoothMPC_bmy85[3];
int smoothMPC_lbIdx85[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb85 = smoothMPC_l + 1360;
smoothMPC_FLOAT* smoothMPC_slb85 = smoothMPC_s + 1360;
smoothMPC_FLOAT* smoothMPC_llbbyslb85 = smoothMPC_lbys + 1360;
smoothMPC_FLOAT smoothMPC_rilb85[11];
smoothMPC_FLOAT* smoothMPC_dllbaff85 = smoothMPC_dl_aff + 1360;
smoothMPC_FLOAT* smoothMPC_dslbaff85 = smoothMPC_ds_aff + 1360;
smoothMPC_FLOAT* smoothMPC_dllbcc85 = smoothMPC_dl_cc + 1360;
smoothMPC_FLOAT* smoothMPC_dslbcc85 = smoothMPC_ds_cc + 1360;
smoothMPC_FLOAT* smoothMPC_ccrhsl85 = smoothMPC_ccrhs + 1360;
int smoothMPC_ubIdx85[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub85 = smoothMPC_l + 1371;
smoothMPC_FLOAT* smoothMPC_sub85 = smoothMPC_s + 1371;
smoothMPC_FLOAT* smoothMPC_lubbysub85 = smoothMPC_lbys + 1371;
smoothMPC_FLOAT smoothMPC_riub85[5];
smoothMPC_FLOAT* smoothMPC_dlubaff85 = smoothMPC_dl_aff + 1371;
smoothMPC_FLOAT* smoothMPC_dsubaff85 = smoothMPC_ds_aff + 1371;
smoothMPC_FLOAT* smoothMPC_dlubcc85 = smoothMPC_dl_cc + 1371;
smoothMPC_FLOAT* smoothMPC_dsubcc85 = smoothMPC_ds_cc + 1371;
smoothMPC_FLOAT* smoothMPC_ccrhsub85 = smoothMPC_ccrhs + 1371;
smoothMPC_FLOAT smoothMPC_Phi85[11];
smoothMPC_FLOAT smoothMPC_W85[11];
smoothMPC_FLOAT smoothMPC_Ysd85[9];
smoothMPC_FLOAT smoothMPC_Lsd85[9];
smoothMPC_FLOAT* smoothMPC_z86 = smoothMPC_z + 946;
smoothMPC_FLOAT* smoothMPC_dzaff86 = smoothMPC_dz_aff + 946;
smoothMPC_FLOAT* smoothMPC_dzcc86 = smoothMPC_dz_cc + 946;
smoothMPC_FLOAT* smoothMPC_rd86 = smoothMPC_rd + 946;
smoothMPC_FLOAT smoothMPC_Lbyrd86[11];
smoothMPC_FLOAT* smoothMPC_grad_cost86 = smoothMPC_grad_cost + 946;
smoothMPC_FLOAT* smoothMPC_grad_eq86 = smoothMPC_grad_eq + 946;
smoothMPC_FLOAT* smoothMPC_grad_ineq86 = smoothMPC_grad_ineq + 946;
smoothMPC_FLOAT smoothMPC_ctv86[11];
smoothMPC_FLOAT* smoothMPC_v86 = smoothMPC_v + 258;
smoothMPC_FLOAT smoothMPC_re86[3];
smoothMPC_FLOAT smoothMPC_beta86[3];
smoothMPC_FLOAT smoothMPC_betacc86[3];
smoothMPC_FLOAT* smoothMPC_dvaff86 = smoothMPC_dv_aff + 258;
smoothMPC_FLOAT* smoothMPC_dvcc86 = smoothMPC_dv_cc + 258;
smoothMPC_FLOAT smoothMPC_V86[33];
smoothMPC_FLOAT smoothMPC_Yd86[6];
smoothMPC_FLOAT smoothMPC_Ld86[6];
smoothMPC_FLOAT smoothMPC_yy86[3];
smoothMPC_FLOAT smoothMPC_bmy86[3];
int smoothMPC_lbIdx86[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb86 = smoothMPC_l + 1376;
smoothMPC_FLOAT* smoothMPC_slb86 = smoothMPC_s + 1376;
smoothMPC_FLOAT* smoothMPC_llbbyslb86 = smoothMPC_lbys + 1376;
smoothMPC_FLOAT smoothMPC_rilb86[11];
smoothMPC_FLOAT* smoothMPC_dllbaff86 = smoothMPC_dl_aff + 1376;
smoothMPC_FLOAT* smoothMPC_dslbaff86 = smoothMPC_ds_aff + 1376;
smoothMPC_FLOAT* smoothMPC_dllbcc86 = smoothMPC_dl_cc + 1376;
smoothMPC_FLOAT* smoothMPC_dslbcc86 = smoothMPC_ds_cc + 1376;
smoothMPC_FLOAT* smoothMPC_ccrhsl86 = smoothMPC_ccrhs + 1376;
int smoothMPC_ubIdx86[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub86 = smoothMPC_l + 1387;
smoothMPC_FLOAT* smoothMPC_sub86 = smoothMPC_s + 1387;
smoothMPC_FLOAT* smoothMPC_lubbysub86 = smoothMPC_lbys + 1387;
smoothMPC_FLOAT smoothMPC_riub86[5];
smoothMPC_FLOAT* smoothMPC_dlubaff86 = smoothMPC_dl_aff + 1387;
smoothMPC_FLOAT* smoothMPC_dsubaff86 = smoothMPC_ds_aff + 1387;
smoothMPC_FLOAT* smoothMPC_dlubcc86 = smoothMPC_dl_cc + 1387;
smoothMPC_FLOAT* smoothMPC_dsubcc86 = smoothMPC_ds_cc + 1387;
smoothMPC_FLOAT* smoothMPC_ccrhsub86 = smoothMPC_ccrhs + 1387;
smoothMPC_FLOAT smoothMPC_Phi86[11];
smoothMPC_FLOAT smoothMPC_W86[11];
smoothMPC_FLOAT smoothMPC_Ysd86[9];
smoothMPC_FLOAT smoothMPC_Lsd86[9];
smoothMPC_FLOAT* smoothMPC_z87 = smoothMPC_z + 957;
smoothMPC_FLOAT* smoothMPC_dzaff87 = smoothMPC_dz_aff + 957;
smoothMPC_FLOAT* smoothMPC_dzcc87 = smoothMPC_dz_cc + 957;
smoothMPC_FLOAT* smoothMPC_rd87 = smoothMPC_rd + 957;
smoothMPC_FLOAT smoothMPC_Lbyrd87[11];
smoothMPC_FLOAT* smoothMPC_grad_cost87 = smoothMPC_grad_cost + 957;
smoothMPC_FLOAT* smoothMPC_grad_eq87 = smoothMPC_grad_eq + 957;
smoothMPC_FLOAT* smoothMPC_grad_ineq87 = smoothMPC_grad_ineq + 957;
smoothMPC_FLOAT smoothMPC_ctv87[11];
smoothMPC_FLOAT* smoothMPC_v87 = smoothMPC_v + 261;
smoothMPC_FLOAT smoothMPC_re87[3];
smoothMPC_FLOAT smoothMPC_beta87[3];
smoothMPC_FLOAT smoothMPC_betacc87[3];
smoothMPC_FLOAT* smoothMPC_dvaff87 = smoothMPC_dv_aff + 261;
smoothMPC_FLOAT* smoothMPC_dvcc87 = smoothMPC_dv_cc + 261;
smoothMPC_FLOAT smoothMPC_V87[33];
smoothMPC_FLOAT smoothMPC_Yd87[6];
smoothMPC_FLOAT smoothMPC_Ld87[6];
smoothMPC_FLOAT smoothMPC_yy87[3];
smoothMPC_FLOAT smoothMPC_bmy87[3];
int smoothMPC_lbIdx87[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb87 = smoothMPC_l + 1392;
smoothMPC_FLOAT* smoothMPC_slb87 = smoothMPC_s + 1392;
smoothMPC_FLOAT* smoothMPC_llbbyslb87 = smoothMPC_lbys + 1392;
smoothMPC_FLOAT smoothMPC_rilb87[11];
smoothMPC_FLOAT* smoothMPC_dllbaff87 = smoothMPC_dl_aff + 1392;
smoothMPC_FLOAT* smoothMPC_dslbaff87 = smoothMPC_ds_aff + 1392;
smoothMPC_FLOAT* smoothMPC_dllbcc87 = smoothMPC_dl_cc + 1392;
smoothMPC_FLOAT* smoothMPC_dslbcc87 = smoothMPC_ds_cc + 1392;
smoothMPC_FLOAT* smoothMPC_ccrhsl87 = smoothMPC_ccrhs + 1392;
int smoothMPC_ubIdx87[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub87 = smoothMPC_l + 1403;
smoothMPC_FLOAT* smoothMPC_sub87 = smoothMPC_s + 1403;
smoothMPC_FLOAT* smoothMPC_lubbysub87 = smoothMPC_lbys + 1403;
smoothMPC_FLOAT smoothMPC_riub87[5];
smoothMPC_FLOAT* smoothMPC_dlubaff87 = smoothMPC_dl_aff + 1403;
smoothMPC_FLOAT* smoothMPC_dsubaff87 = smoothMPC_ds_aff + 1403;
smoothMPC_FLOAT* smoothMPC_dlubcc87 = smoothMPC_dl_cc + 1403;
smoothMPC_FLOAT* smoothMPC_dsubcc87 = smoothMPC_ds_cc + 1403;
smoothMPC_FLOAT* smoothMPC_ccrhsub87 = smoothMPC_ccrhs + 1403;
smoothMPC_FLOAT smoothMPC_Phi87[11];
smoothMPC_FLOAT smoothMPC_W87[11];
smoothMPC_FLOAT smoothMPC_Ysd87[9];
smoothMPC_FLOAT smoothMPC_Lsd87[9];
smoothMPC_FLOAT* smoothMPC_z88 = smoothMPC_z + 968;
smoothMPC_FLOAT* smoothMPC_dzaff88 = smoothMPC_dz_aff + 968;
smoothMPC_FLOAT* smoothMPC_dzcc88 = smoothMPC_dz_cc + 968;
smoothMPC_FLOAT* smoothMPC_rd88 = smoothMPC_rd + 968;
smoothMPC_FLOAT smoothMPC_Lbyrd88[11];
smoothMPC_FLOAT* smoothMPC_grad_cost88 = smoothMPC_grad_cost + 968;
smoothMPC_FLOAT* smoothMPC_grad_eq88 = smoothMPC_grad_eq + 968;
smoothMPC_FLOAT* smoothMPC_grad_ineq88 = smoothMPC_grad_ineq + 968;
smoothMPC_FLOAT smoothMPC_ctv88[11];
smoothMPC_FLOAT* smoothMPC_v88 = smoothMPC_v + 264;
smoothMPC_FLOAT smoothMPC_re88[3];
smoothMPC_FLOAT smoothMPC_beta88[3];
smoothMPC_FLOAT smoothMPC_betacc88[3];
smoothMPC_FLOAT* smoothMPC_dvaff88 = smoothMPC_dv_aff + 264;
smoothMPC_FLOAT* smoothMPC_dvcc88 = smoothMPC_dv_cc + 264;
smoothMPC_FLOAT smoothMPC_V88[33];
smoothMPC_FLOAT smoothMPC_Yd88[6];
smoothMPC_FLOAT smoothMPC_Ld88[6];
smoothMPC_FLOAT smoothMPC_yy88[3];
smoothMPC_FLOAT smoothMPC_bmy88[3];
int smoothMPC_lbIdx88[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb88 = smoothMPC_l + 1408;
smoothMPC_FLOAT* smoothMPC_slb88 = smoothMPC_s + 1408;
smoothMPC_FLOAT* smoothMPC_llbbyslb88 = smoothMPC_lbys + 1408;
smoothMPC_FLOAT smoothMPC_rilb88[11];
smoothMPC_FLOAT* smoothMPC_dllbaff88 = smoothMPC_dl_aff + 1408;
smoothMPC_FLOAT* smoothMPC_dslbaff88 = smoothMPC_ds_aff + 1408;
smoothMPC_FLOAT* smoothMPC_dllbcc88 = smoothMPC_dl_cc + 1408;
smoothMPC_FLOAT* smoothMPC_dslbcc88 = smoothMPC_ds_cc + 1408;
smoothMPC_FLOAT* smoothMPC_ccrhsl88 = smoothMPC_ccrhs + 1408;
int smoothMPC_ubIdx88[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub88 = smoothMPC_l + 1419;
smoothMPC_FLOAT* smoothMPC_sub88 = smoothMPC_s + 1419;
smoothMPC_FLOAT* smoothMPC_lubbysub88 = smoothMPC_lbys + 1419;
smoothMPC_FLOAT smoothMPC_riub88[5];
smoothMPC_FLOAT* smoothMPC_dlubaff88 = smoothMPC_dl_aff + 1419;
smoothMPC_FLOAT* smoothMPC_dsubaff88 = smoothMPC_ds_aff + 1419;
smoothMPC_FLOAT* smoothMPC_dlubcc88 = smoothMPC_dl_cc + 1419;
smoothMPC_FLOAT* smoothMPC_dsubcc88 = smoothMPC_ds_cc + 1419;
smoothMPC_FLOAT* smoothMPC_ccrhsub88 = smoothMPC_ccrhs + 1419;
smoothMPC_FLOAT smoothMPC_Phi88[11];
smoothMPC_FLOAT smoothMPC_W88[11];
smoothMPC_FLOAT smoothMPC_Ysd88[9];
smoothMPC_FLOAT smoothMPC_Lsd88[9];
smoothMPC_FLOAT* smoothMPC_z89 = smoothMPC_z + 979;
smoothMPC_FLOAT* smoothMPC_dzaff89 = smoothMPC_dz_aff + 979;
smoothMPC_FLOAT* smoothMPC_dzcc89 = smoothMPC_dz_cc + 979;
smoothMPC_FLOAT* smoothMPC_rd89 = smoothMPC_rd + 979;
smoothMPC_FLOAT smoothMPC_Lbyrd89[11];
smoothMPC_FLOAT* smoothMPC_grad_cost89 = smoothMPC_grad_cost + 979;
smoothMPC_FLOAT* smoothMPC_grad_eq89 = smoothMPC_grad_eq + 979;
smoothMPC_FLOAT* smoothMPC_grad_ineq89 = smoothMPC_grad_ineq + 979;
smoothMPC_FLOAT smoothMPC_ctv89[11];
smoothMPC_FLOAT* smoothMPC_v89 = smoothMPC_v + 267;
smoothMPC_FLOAT smoothMPC_re89[3];
smoothMPC_FLOAT smoothMPC_beta89[3];
smoothMPC_FLOAT smoothMPC_betacc89[3];
smoothMPC_FLOAT* smoothMPC_dvaff89 = smoothMPC_dv_aff + 267;
smoothMPC_FLOAT* smoothMPC_dvcc89 = smoothMPC_dv_cc + 267;
smoothMPC_FLOAT smoothMPC_V89[33];
smoothMPC_FLOAT smoothMPC_Yd89[6];
smoothMPC_FLOAT smoothMPC_Ld89[6];
smoothMPC_FLOAT smoothMPC_yy89[3];
smoothMPC_FLOAT smoothMPC_bmy89[3];
int smoothMPC_lbIdx89[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb89 = smoothMPC_l + 1424;
smoothMPC_FLOAT* smoothMPC_slb89 = smoothMPC_s + 1424;
smoothMPC_FLOAT* smoothMPC_llbbyslb89 = smoothMPC_lbys + 1424;
smoothMPC_FLOAT smoothMPC_rilb89[11];
smoothMPC_FLOAT* smoothMPC_dllbaff89 = smoothMPC_dl_aff + 1424;
smoothMPC_FLOAT* smoothMPC_dslbaff89 = smoothMPC_ds_aff + 1424;
smoothMPC_FLOAT* smoothMPC_dllbcc89 = smoothMPC_dl_cc + 1424;
smoothMPC_FLOAT* smoothMPC_dslbcc89 = smoothMPC_ds_cc + 1424;
smoothMPC_FLOAT* smoothMPC_ccrhsl89 = smoothMPC_ccrhs + 1424;
int smoothMPC_ubIdx89[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub89 = smoothMPC_l + 1435;
smoothMPC_FLOAT* smoothMPC_sub89 = smoothMPC_s + 1435;
smoothMPC_FLOAT* smoothMPC_lubbysub89 = smoothMPC_lbys + 1435;
smoothMPC_FLOAT smoothMPC_riub89[5];
smoothMPC_FLOAT* smoothMPC_dlubaff89 = smoothMPC_dl_aff + 1435;
smoothMPC_FLOAT* smoothMPC_dsubaff89 = smoothMPC_ds_aff + 1435;
smoothMPC_FLOAT* smoothMPC_dlubcc89 = smoothMPC_dl_cc + 1435;
smoothMPC_FLOAT* smoothMPC_dsubcc89 = smoothMPC_ds_cc + 1435;
smoothMPC_FLOAT* smoothMPC_ccrhsub89 = smoothMPC_ccrhs + 1435;
smoothMPC_FLOAT smoothMPC_Phi89[11];
smoothMPC_FLOAT smoothMPC_W89[11];
smoothMPC_FLOAT smoothMPC_Ysd89[9];
smoothMPC_FLOAT smoothMPC_Lsd89[9];
smoothMPC_FLOAT* smoothMPC_z90 = smoothMPC_z + 990;
smoothMPC_FLOAT* smoothMPC_dzaff90 = smoothMPC_dz_aff + 990;
smoothMPC_FLOAT* smoothMPC_dzcc90 = smoothMPC_dz_cc + 990;
smoothMPC_FLOAT* smoothMPC_rd90 = smoothMPC_rd + 990;
smoothMPC_FLOAT smoothMPC_Lbyrd90[11];
smoothMPC_FLOAT* smoothMPC_grad_cost90 = smoothMPC_grad_cost + 990;
smoothMPC_FLOAT* smoothMPC_grad_eq90 = smoothMPC_grad_eq + 990;
smoothMPC_FLOAT* smoothMPC_grad_ineq90 = smoothMPC_grad_ineq + 990;
smoothMPC_FLOAT smoothMPC_ctv90[11];
smoothMPC_FLOAT* smoothMPC_v90 = smoothMPC_v + 270;
smoothMPC_FLOAT smoothMPC_re90[3];
smoothMPC_FLOAT smoothMPC_beta90[3];
smoothMPC_FLOAT smoothMPC_betacc90[3];
smoothMPC_FLOAT* smoothMPC_dvaff90 = smoothMPC_dv_aff + 270;
smoothMPC_FLOAT* smoothMPC_dvcc90 = smoothMPC_dv_cc + 270;
smoothMPC_FLOAT smoothMPC_V90[33];
smoothMPC_FLOAT smoothMPC_Yd90[6];
smoothMPC_FLOAT smoothMPC_Ld90[6];
smoothMPC_FLOAT smoothMPC_yy90[3];
smoothMPC_FLOAT smoothMPC_bmy90[3];
int smoothMPC_lbIdx90[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb90 = smoothMPC_l + 1440;
smoothMPC_FLOAT* smoothMPC_slb90 = smoothMPC_s + 1440;
smoothMPC_FLOAT* smoothMPC_llbbyslb90 = smoothMPC_lbys + 1440;
smoothMPC_FLOAT smoothMPC_rilb90[11];
smoothMPC_FLOAT* smoothMPC_dllbaff90 = smoothMPC_dl_aff + 1440;
smoothMPC_FLOAT* smoothMPC_dslbaff90 = smoothMPC_ds_aff + 1440;
smoothMPC_FLOAT* smoothMPC_dllbcc90 = smoothMPC_dl_cc + 1440;
smoothMPC_FLOAT* smoothMPC_dslbcc90 = smoothMPC_ds_cc + 1440;
smoothMPC_FLOAT* smoothMPC_ccrhsl90 = smoothMPC_ccrhs + 1440;
int smoothMPC_ubIdx90[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub90 = smoothMPC_l + 1451;
smoothMPC_FLOAT* smoothMPC_sub90 = smoothMPC_s + 1451;
smoothMPC_FLOAT* smoothMPC_lubbysub90 = smoothMPC_lbys + 1451;
smoothMPC_FLOAT smoothMPC_riub90[5];
smoothMPC_FLOAT* smoothMPC_dlubaff90 = smoothMPC_dl_aff + 1451;
smoothMPC_FLOAT* smoothMPC_dsubaff90 = smoothMPC_ds_aff + 1451;
smoothMPC_FLOAT* smoothMPC_dlubcc90 = smoothMPC_dl_cc + 1451;
smoothMPC_FLOAT* smoothMPC_dsubcc90 = smoothMPC_ds_cc + 1451;
smoothMPC_FLOAT* smoothMPC_ccrhsub90 = smoothMPC_ccrhs + 1451;
smoothMPC_FLOAT smoothMPC_Phi90[11];
smoothMPC_FLOAT smoothMPC_W90[11];
smoothMPC_FLOAT smoothMPC_Ysd90[9];
smoothMPC_FLOAT smoothMPC_Lsd90[9];
smoothMPC_FLOAT* smoothMPC_z91 = smoothMPC_z + 1001;
smoothMPC_FLOAT* smoothMPC_dzaff91 = smoothMPC_dz_aff + 1001;
smoothMPC_FLOAT* smoothMPC_dzcc91 = smoothMPC_dz_cc + 1001;
smoothMPC_FLOAT* smoothMPC_rd91 = smoothMPC_rd + 1001;
smoothMPC_FLOAT smoothMPC_Lbyrd91[11];
smoothMPC_FLOAT* smoothMPC_grad_cost91 = smoothMPC_grad_cost + 1001;
smoothMPC_FLOAT* smoothMPC_grad_eq91 = smoothMPC_grad_eq + 1001;
smoothMPC_FLOAT* smoothMPC_grad_ineq91 = smoothMPC_grad_ineq + 1001;
smoothMPC_FLOAT smoothMPC_ctv91[11];
smoothMPC_FLOAT* smoothMPC_v91 = smoothMPC_v + 273;
smoothMPC_FLOAT smoothMPC_re91[3];
smoothMPC_FLOAT smoothMPC_beta91[3];
smoothMPC_FLOAT smoothMPC_betacc91[3];
smoothMPC_FLOAT* smoothMPC_dvaff91 = smoothMPC_dv_aff + 273;
smoothMPC_FLOAT* smoothMPC_dvcc91 = smoothMPC_dv_cc + 273;
smoothMPC_FLOAT smoothMPC_V91[33];
smoothMPC_FLOAT smoothMPC_Yd91[6];
smoothMPC_FLOAT smoothMPC_Ld91[6];
smoothMPC_FLOAT smoothMPC_yy91[3];
smoothMPC_FLOAT smoothMPC_bmy91[3];
int smoothMPC_lbIdx91[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb91 = smoothMPC_l + 1456;
smoothMPC_FLOAT* smoothMPC_slb91 = smoothMPC_s + 1456;
smoothMPC_FLOAT* smoothMPC_llbbyslb91 = smoothMPC_lbys + 1456;
smoothMPC_FLOAT smoothMPC_rilb91[11];
smoothMPC_FLOAT* smoothMPC_dllbaff91 = smoothMPC_dl_aff + 1456;
smoothMPC_FLOAT* smoothMPC_dslbaff91 = smoothMPC_ds_aff + 1456;
smoothMPC_FLOAT* smoothMPC_dllbcc91 = smoothMPC_dl_cc + 1456;
smoothMPC_FLOAT* smoothMPC_dslbcc91 = smoothMPC_ds_cc + 1456;
smoothMPC_FLOAT* smoothMPC_ccrhsl91 = smoothMPC_ccrhs + 1456;
int smoothMPC_ubIdx91[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub91 = smoothMPC_l + 1467;
smoothMPC_FLOAT* smoothMPC_sub91 = smoothMPC_s + 1467;
smoothMPC_FLOAT* smoothMPC_lubbysub91 = smoothMPC_lbys + 1467;
smoothMPC_FLOAT smoothMPC_riub91[5];
smoothMPC_FLOAT* smoothMPC_dlubaff91 = smoothMPC_dl_aff + 1467;
smoothMPC_FLOAT* smoothMPC_dsubaff91 = smoothMPC_ds_aff + 1467;
smoothMPC_FLOAT* smoothMPC_dlubcc91 = smoothMPC_dl_cc + 1467;
smoothMPC_FLOAT* smoothMPC_dsubcc91 = smoothMPC_ds_cc + 1467;
smoothMPC_FLOAT* smoothMPC_ccrhsub91 = smoothMPC_ccrhs + 1467;
smoothMPC_FLOAT smoothMPC_Phi91[11];
smoothMPC_FLOAT smoothMPC_W91[11];
smoothMPC_FLOAT smoothMPC_Ysd91[9];
smoothMPC_FLOAT smoothMPC_Lsd91[9];
smoothMPC_FLOAT* smoothMPC_z92 = smoothMPC_z + 1012;
smoothMPC_FLOAT* smoothMPC_dzaff92 = smoothMPC_dz_aff + 1012;
smoothMPC_FLOAT* smoothMPC_dzcc92 = smoothMPC_dz_cc + 1012;
smoothMPC_FLOAT* smoothMPC_rd92 = smoothMPC_rd + 1012;
smoothMPC_FLOAT smoothMPC_Lbyrd92[11];
smoothMPC_FLOAT* smoothMPC_grad_cost92 = smoothMPC_grad_cost + 1012;
smoothMPC_FLOAT* smoothMPC_grad_eq92 = smoothMPC_grad_eq + 1012;
smoothMPC_FLOAT* smoothMPC_grad_ineq92 = smoothMPC_grad_ineq + 1012;
smoothMPC_FLOAT smoothMPC_ctv92[11];
smoothMPC_FLOAT* smoothMPC_v92 = smoothMPC_v + 276;
smoothMPC_FLOAT smoothMPC_re92[3];
smoothMPC_FLOAT smoothMPC_beta92[3];
smoothMPC_FLOAT smoothMPC_betacc92[3];
smoothMPC_FLOAT* smoothMPC_dvaff92 = smoothMPC_dv_aff + 276;
smoothMPC_FLOAT* smoothMPC_dvcc92 = smoothMPC_dv_cc + 276;
smoothMPC_FLOAT smoothMPC_V92[33];
smoothMPC_FLOAT smoothMPC_Yd92[6];
smoothMPC_FLOAT smoothMPC_Ld92[6];
smoothMPC_FLOAT smoothMPC_yy92[3];
smoothMPC_FLOAT smoothMPC_bmy92[3];
int smoothMPC_lbIdx92[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb92 = smoothMPC_l + 1472;
smoothMPC_FLOAT* smoothMPC_slb92 = smoothMPC_s + 1472;
smoothMPC_FLOAT* smoothMPC_llbbyslb92 = smoothMPC_lbys + 1472;
smoothMPC_FLOAT smoothMPC_rilb92[11];
smoothMPC_FLOAT* smoothMPC_dllbaff92 = smoothMPC_dl_aff + 1472;
smoothMPC_FLOAT* smoothMPC_dslbaff92 = smoothMPC_ds_aff + 1472;
smoothMPC_FLOAT* smoothMPC_dllbcc92 = smoothMPC_dl_cc + 1472;
smoothMPC_FLOAT* smoothMPC_dslbcc92 = smoothMPC_ds_cc + 1472;
smoothMPC_FLOAT* smoothMPC_ccrhsl92 = smoothMPC_ccrhs + 1472;
int smoothMPC_ubIdx92[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub92 = smoothMPC_l + 1483;
smoothMPC_FLOAT* smoothMPC_sub92 = smoothMPC_s + 1483;
smoothMPC_FLOAT* smoothMPC_lubbysub92 = smoothMPC_lbys + 1483;
smoothMPC_FLOAT smoothMPC_riub92[5];
smoothMPC_FLOAT* smoothMPC_dlubaff92 = smoothMPC_dl_aff + 1483;
smoothMPC_FLOAT* smoothMPC_dsubaff92 = smoothMPC_ds_aff + 1483;
smoothMPC_FLOAT* smoothMPC_dlubcc92 = smoothMPC_dl_cc + 1483;
smoothMPC_FLOAT* smoothMPC_dsubcc92 = smoothMPC_ds_cc + 1483;
smoothMPC_FLOAT* smoothMPC_ccrhsub92 = smoothMPC_ccrhs + 1483;
smoothMPC_FLOAT smoothMPC_Phi92[11];
smoothMPC_FLOAT smoothMPC_W92[11];
smoothMPC_FLOAT smoothMPC_Ysd92[9];
smoothMPC_FLOAT smoothMPC_Lsd92[9];
smoothMPC_FLOAT* smoothMPC_z93 = smoothMPC_z + 1023;
smoothMPC_FLOAT* smoothMPC_dzaff93 = smoothMPC_dz_aff + 1023;
smoothMPC_FLOAT* smoothMPC_dzcc93 = smoothMPC_dz_cc + 1023;
smoothMPC_FLOAT* smoothMPC_rd93 = smoothMPC_rd + 1023;
smoothMPC_FLOAT smoothMPC_Lbyrd93[11];
smoothMPC_FLOAT* smoothMPC_grad_cost93 = smoothMPC_grad_cost + 1023;
smoothMPC_FLOAT* smoothMPC_grad_eq93 = smoothMPC_grad_eq + 1023;
smoothMPC_FLOAT* smoothMPC_grad_ineq93 = smoothMPC_grad_ineq + 1023;
smoothMPC_FLOAT smoothMPC_ctv93[11];
smoothMPC_FLOAT* smoothMPC_v93 = smoothMPC_v + 279;
smoothMPC_FLOAT smoothMPC_re93[3];
smoothMPC_FLOAT smoothMPC_beta93[3];
smoothMPC_FLOAT smoothMPC_betacc93[3];
smoothMPC_FLOAT* smoothMPC_dvaff93 = smoothMPC_dv_aff + 279;
smoothMPC_FLOAT* smoothMPC_dvcc93 = smoothMPC_dv_cc + 279;
smoothMPC_FLOAT smoothMPC_V93[33];
smoothMPC_FLOAT smoothMPC_Yd93[6];
smoothMPC_FLOAT smoothMPC_Ld93[6];
smoothMPC_FLOAT smoothMPC_yy93[3];
smoothMPC_FLOAT smoothMPC_bmy93[3];
int smoothMPC_lbIdx93[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb93 = smoothMPC_l + 1488;
smoothMPC_FLOAT* smoothMPC_slb93 = smoothMPC_s + 1488;
smoothMPC_FLOAT* smoothMPC_llbbyslb93 = smoothMPC_lbys + 1488;
smoothMPC_FLOAT smoothMPC_rilb93[11];
smoothMPC_FLOAT* smoothMPC_dllbaff93 = smoothMPC_dl_aff + 1488;
smoothMPC_FLOAT* smoothMPC_dslbaff93 = smoothMPC_ds_aff + 1488;
smoothMPC_FLOAT* smoothMPC_dllbcc93 = smoothMPC_dl_cc + 1488;
smoothMPC_FLOAT* smoothMPC_dslbcc93 = smoothMPC_ds_cc + 1488;
smoothMPC_FLOAT* smoothMPC_ccrhsl93 = smoothMPC_ccrhs + 1488;
int smoothMPC_ubIdx93[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub93 = smoothMPC_l + 1499;
smoothMPC_FLOAT* smoothMPC_sub93 = smoothMPC_s + 1499;
smoothMPC_FLOAT* smoothMPC_lubbysub93 = smoothMPC_lbys + 1499;
smoothMPC_FLOAT smoothMPC_riub93[5];
smoothMPC_FLOAT* smoothMPC_dlubaff93 = smoothMPC_dl_aff + 1499;
smoothMPC_FLOAT* smoothMPC_dsubaff93 = smoothMPC_ds_aff + 1499;
smoothMPC_FLOAT* smoothMPC_dlubcc93 = smoothMPC_dl_cc + 1499;
smoothMPC_FLOAT* smoothMPC_dsubcc93 = smoothMPC_ds_cc + 1499;
smoothMPC_FLOAT* smoothMPC_ccrhsub93 = smoothMPC_ccrhs + 1499;
smoothMPC_FLOAT smoothMPC_Phi93[11];
smoothMPC_FLOAT smoothMPC_W93[11];
smoothMPC_FLOAT smoothMPC_Ysd93[9];
smoothMPC_FLOAT smoothMPC_Lsd93[9];
smoothMPC_FLOAT* smoothMPC_z94 = smoothMPC_z + 1034;
smoothMPC_FLOAT* smoothMPC_dzaff94 = smoothMPC_dz_aff + 1034;
smoothMPC_FLOAT* smoothMPC_dzcc94 = smoothMPC_dz_cc + 1034;
smoothMPC_FLOAT* smoothMPC_rd94 = smoothMPC_rd + 1034;
smoothMPC_FLOAT smoothMPC_Lbyrd94[11];
smoothMPC_FLOAT* smoothMPC_grad_cost94 = smoothMPC_grad_cost + 1034;
smoothMPC_FLOAT* smoothMPC_grad_eq94 = smoothMPC_grad_eq + 1034;
smoothMPC_FLOAT* smoothMPC_grad_ineq94 = smoothMPC_grad_ineq + 1034;
smoothMPC_FLOAT smoothMPC_ctv94[11];
smoothMPC_FLOAT* smoothMPC_v94 = smoothMPC_v + 282;
smoothMPC_FLOAT smoothMPC_re94[3];
smoothMPC_FLOAT smoothMPC_beta94[3];
smoothMPC_FLOAT smoothMPC_betacc94[3];
smoothMPC_FLOAT* smoothMPC_dvaff94 = smoothMPC_dv_aff + 282;
smoothMPC_FLOAT* smoothMPC_dvcc94 = smoothMPC_dv_cc + 282;
smoothMPC_FLOAT smoothMPC_V94[33];
smoothMPC_FLOAT smoothMPC_Yd94[6];
smoothMPC_FLOAT smoothMPC_Ld94[6];
smoothMPC_FLOAT smoothMPC_yy94[3];
smoothMPC_FLOAT smoothMPC_bmy94[3];
int smoothMPC_lbIdx94[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb94 = smoothMPC_l + 1504;
smoothMPC_FLOAT* smoothMPC_slb94 = smoothMPC_s + 1504;
smoothMPC_FLOAT* smoothMPC_llbbyslb94 = smoothMPC_lbys + 1504;
smoothMPC_FLOAT smoothMPC_rilb94[11];
smoothMPC_FLOAT* smoothMPC_dllbaff94 = smoothMPC_dl_aff + 1504;
smoothMPC_FLOAT* smoothMPC_dslbaff94 = smoothMPC_ds_aff + 1504;
smoothMPC_FLOAT* smoothMPC_dllbcc94 = smoothMPC_dl_cc + 1504;
smoothMPC_FLOAT* smoothMPC_dslbcc94 = smoothMPC_ds_cc + 1504;
smoothMPC_FLOAT* smoothMPC_ccrhsl94 = smoothMPC_ccrhs + 1504;
int smoothMPC_ubIdx94[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub94 = smoothMPC_l + 1515;
smoothMPC_FLOAT* smoothMPC_sub94 = smoothMPC_s + 1515;
smoothMPC_FLOAT* smoothMPC_lubbysub94 = smoothMPC_lbys + 1515;
smoothMPC_FLOAT smoothMPC_riub94[5];
smoothMPC_FLOAT* smoothMPC_dlubaff94 = smoothMPC_dl_aff + 1515;
smoothMPC_FLOAT* smoothMPC_dsubaff94 = smoothMPC_ds_aff + 1515;
smoothMPC_FLOAT* smoothMPC_dlubcc94 = smoothMPC_dl_cc + 1515;
smoothMPC_FLOAT* smoothMPC_dsubcc94 = smoothMPC_ds_cc + 1515;
smoothMPC_FLOAT* smoothMPC_ccrhsub94 = smoothMPC_ccrhs + 1515;
smoothMPC_FLOAT smoothMPC_Phi94[11];
smoothMPC_FLOAT smoothMPC_W94[11];
smoothMPC_FLOAT smoothMPC_Ysd94[9];
smoothMPC_FLOAT smoothMPC_Lsd94[9];
smoothMPC_FLOAT* smoothMPC_z95 = smoothMPC_z + 1045;
smoothMPC_FLOAT* smoothMPC_dzaff95 = smoothMPC_dz_aff + 1045;
smoothMPC_FLOAT* smoothMPC_dzcc95 = smoothMPC_dz_cc + 1045;
smoothMPC_FLOAT* smoothMPC_rd95 = smoothMPC_rd + 1045;
smoothMPC_FLOAT smoothMPC_Lbyrd95[11];
smoothMPC_FLOAT* smoothMPC_grad_cost95 = smoothMPC_grad_cost + 1045;
smoothMPC_FLOAT* smoothMPC_grad_eq95 = smoothMPC_grad_eq + 1045;
smoothMPC_FLOAT* smoothMPC_grad_ineq95 = smoothMPC_grad_ineq + 1045;
smoothMPC_FLOAT smoothMPC_ctv95[11];
smoothMPC_FLOAT* smoothMPC_v95 = smoothMPC_v + 285;
smoothMPC_FLOAT smoothMPC_re95[3];
smoothMPC_FLOAT smoothMPC_beta95[3];
smoothMPC_FLOAT smoothMPC_betacc95[3];
smoothMPC_FLOAT* smoothMPC_dvaff95 = smoothMPC_dv_aff + 285;
smoothMPC_FLOAT* smoothMPC_dvcc95 = smoothMPC_dv_cc + 285;
smoothMPC_FLOAT smoothMPC_V95[33];
smoothMPC_FLOAT smoothMPC_Yd95[6];
smoothMPC_FLOAT smoothMPC_Ld95[6];
smoothMPC_FLOAT smoothMPC_yy95[3];
smoothMPC_FLOAT smoothMPC_bmy95[3];
int smoothMPC_lbIdx95[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb95 = smoothMPC_l + 1520;
smoothMPC_FLOAT* smoothMPC_slb95 = smoothMPC_s + 1520;
smoothMPC_FLOAT* smoothMPC_llbbyslb95 = smoothMPC_lbys + 1520;
smoothMPC_FLOAT smoothMPC_rilb95[11];
smoothMPC_FLOAT* smoothMPC_dllbaff95 = smoothMPC_dl_aff + 1520;
smoothMPC_FLOAT* smoothMPC_dslbaff95 = smoothMPC_ds_aff + 1520;
smoothMPC_FLOAT* smoothMPC_dllbcc95 = smoothMPC_dl_cc + 1520;
smoothMPC_FLOAT* smoothMPC_dslbcc95 = smoothMPC_ds_cc + 1520;
smoothMPC_FLOAT* smoothMPC_ccrhsl95 = smoothMPC_ccrhs + 1520;
int smoothMPC_ubIdx95[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub95 = smoothMPC_l + 1531;
smoothMPC_FLOAT* smoothMPC_sub95 = smoothMPC_s + 1531;
smoothMPC_FLOAT* smoothMPC_lubbysub95 = smoothMPC_lbys + 1531;
smoothMPC_FLOAT smoothMPC_riub95[5];
smoothMPC_FLOAT* smoothMPC_dlubaff95 = smoothMPC_dl_aff + 1531;
smoothMPC_FLOAT* smoothMPC_dsubaff95 = smoothMPC_ds_aff + 1531;
smoothMPC_FLOAT* smoothMPC_dlubcc95 = smoothMPC_dl_cc + 1531;
smoothMPC_FLOAT* smoothMPC_dsubcc95 = smoothMPC_ds_cc + 1531;
smoothMPC_FLOAT* smoothMPC_ccrhsub95 = smoothMPC_ccrhs + 1531;
smoothMPC_FLOAT smoothMPC_Phi95[11];
smoothMPC_FLOAT smoothMPC_W95[11];
smoothMPC_FLOAT smoothMPC_Ysd95[9];
smoothMPC_FLOAT smoothMPC_Lsd95[9];
smoothMPC_FLOAT* smoothMPC_z96 = smoothMPC_z + 1056;
smoothMPC_FLOAT* smoothMPC_dzaff96 = smoothMPC_dz_aff + 1056;
smoothMPC_FLOAT* smoothMPC_dzcc96 = smoothMPC_dz_cc + 1056;
smoothMPC_FLOAT* smoothMPC_rd96 = smoothMPC_rd + 1056;
smoothMPC_FLOAT smoothMPC_Lbyrd96[11];
smoothMPC_FLOAT* smoothMPC_grad_cost96 = smoothMPC_grad_cost + 1056;
smoothMPC_FLOAT* smoothMPC_grad_eq96 = smoothMPC_grad_eq + 1056;
smoothMPC_FLOAT* smoothMPC_grad_ineq96 = smoothMPC_grad_ineq + 1056;
smoothMPC_FLOAT smoothMPC_ctv96[11];
smoothMPC_FLOAT* smoothMPC_v96 = smoothMPC_v + 288;
smoothMPC_FLOAT smoothMPC_re96[3];
smoothMPC_FLOAT smoothMPC_beta96[3];
smoothMPC_FLOAT smoothMPC_betacc96[3];
smoothMPC_FLOAT* smoothMPC_dvaff96 = smoothMPC_dv_aff + 288;
smoothMPC_FLOAT* smoothMPC_dvcc96 = smoothMPC_dv_cc + 288;
smoothMPC_FLOAT smoothMPC_V96[33];
smoothMPC_FLOAT smoothMPC_Yd96[6];
smoothMPC_FLOAT smoothMPC_Ld96[6];
smoothMPC_FLOAT smoothMPC_yy96[3];
smoothMPC_FLOAT smoothMPC_bmy96[3];
int smoothMPC_lbIdx96[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb96 = smoothMPC_l + 1536;
smoothMPC_FLOAT* smoothMPC_slb96 = smoothMPC_s + 1536;
smoothMPC_FLOAT* smoothMPC_llbbyslb96 = smoothMPC_lbys + 1536;
smoothMPC_FLOAT smoothMPC_rilb96[11];
smoothMPC_FLOAT* smoothMPC_dllbaff96 = smoothMPC_dl_aff + 1536;
smoothMPC_FLOAT* smoothMPC_dslbaff96 = smoothMPC_ds_aff + 1536;
smoothMPC_FLOAT* smoothMPC_dllbcc96 = smoothMPC_dl_cc + 1536;
smoothMPC_FLOAT* smoothMPC_dslbcc96 = smoothMPC_ds_cc + 1536;
smoothMPC_FLOAT* smoothMPC_ccrhsl96 = smoothMPC_ccrhs + 1536;
int smoothMPC_ubIdx96[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub96 = smoothMPC_l + 1547;
smoothMPC_FLOAT* smoothMPC_sub96 = smoothMPC_s + 1547;
smoothMPC_FLOAT* smoothMPC_lubbysub96 = smoothMPC_lbys + 1547;
smoothMPC_FLOAT smoothMPC_riub96[5];
smoothMPC_FLOAT* smoothMPC_dlubaff96 = smoothMPC_dl_aff + 1547;
smoothMPC_FLOAT* smoothMPC_dsubaff96 = smoothMPC_ds_aff + 1547;
smoothMPC_FLOAT* smoothMPC_dlubcc96 = smoothMPC_dl_cc + 1547;
smoothMPC_FLOAT* smoothMPC_dsubcc96 = smoothMPC_ds_cc + 1547;
smoothMPC_FLOAT* smoothMPC_ccrhsub96 = smoothMPC_ccrhs + 1547;
smoothMPC_FLOAT smoothMPC_Phi96[11];
smoothMPC_FLOAT smoothMPC_W96[11];
smoothMPC_FLOAT smoothMPC_Ysd96[9];
smoothMPC_FLOAT smoothMPC_Lsd96[9];
smoothMPC_FLOAT* smoothMPC_z97 = smoothMPC_z + 1067;
smoothMPC_FLOAT* smoothMPC_dzaff97 = smoothMPC_dz_aff + 1067;
smoothMPC_FLOAT* smoothMPC_dzcc97 = smoothMPC_dz_cc + 1067;
smoothMPC_FLOAT* smoothMPC_rd97 = smoothMPC_rd + 1067;
smoothMPC_FLOAT smoothMPC_Lbyrd97[11];
smoothMPC_FLOAT* smoothMPC_grad_cost97 = smoothMPC_grad_cost + 1067;
smoothMPC_FLOAT* smoothMPC_grad_eq97 = smoothMPC_grad_eq + 1067;
smoothMPC_FLOAT* smoothMPC_grad_ineq97 = smoothMPC_grad_ineq + 1067;
smoothMPC_FLOAT smoothMPC_ctv97[11];
smoothMPC_FLOAT* smoothMPC_v97 = smoothMPC_v + 291;
smoothMPC_FLOAT smoothMPC_re97[3];
smoothMPC_FLOAT smoothMPC_beta97[3];
smoothMPC_FLOAT smoothMPC_betacc97[3];
smoothMPC_FLOAT* smoothMPC_dvaff97 = smoothMPC_dv_aff + 291;
smoothMPC_FLOAT* smoothMPC_dvcc97 = smoothMPC_dv_cc + 291;
smoothMPC_FLOAT smoothMPC_V97[33];
smoothMPC_FLOAT smoothMPC_Yd97[6];
smoothMPC_FLOAT smoothMPC_Ld97[6];
smoothMPC_FLOAT smoothMPC_yy97[3];
smoothMPC_FLOAT smoothMPC_bmy97[3];
int smoothMPC_lbIdx97[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb97 = smoothMPC_l + 1552;
smoothMPC_FLOAT* smoothMPC_slb97 = smoothMPC_s + 1552;
smoothMPC_FLOAT* smoothMPC_llbbyslb97 = smoothMPC_lbys + 1552;
smoothMPC_FLOAT smoothMPC_rilb97[11];
smoothMPC_FLOAT* smoothMPC_dllbaff97 = smoothMPC_dl_aff + 1552;
smoothMPC_FLOAT* smoothMPC_dslbaff97 = smoothMPC_ds_aff + 1552;
smoothMPC_FLOAT* smoothMPC_dllbcc97 = smoothMPC_dl_cc + 1552;
smoothMPC_FLOAT* smoothMPC_dslbcc97 = smoothMPC_ds_cc + 1552;
smoothMPC_FLOAT* smoothMPC_ccrhsl97 = smoothMPC_ccrhs + 1552;
int smoothMPC_ubIdx97[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub97 = smoothMPC_l + 1563;
smoothMPC_FLOAT* smoothMPC_sub97 = smoothMPC_s + 1563;
smoothMPC_FLOAT* smoothMPC_lubbysub97 = smoothMPC_lbys + 1563;
smoothMPC_FLOAT smoothMPC_riub97[5];
smoothMPC_FLOAT* smoothMPC_dlubaff97 = smoothMPC_dl_aff + 1563;
smoothMPC_FLOAT* smoothMPC_dsubaff97 = smoothMPC_ds_aff + 1563;
smoothMPC_FLOAT* smoothMPC_dlubcc97 = smoothMPC_dl_cc + 1563;
smoothMPC_FLOAT* smoothMPC_dsubcc97 = smoothMPC_ds_cc + 1563;
smoothMPC_FLOAT* smoothMPC_ccrhsub97 = smoothMPC_ccrhs + 1563;
smoothMPC_FLOAT smoothMPC_Phi97[11];
smoothMPC_FLOAT smoothMPC_W97[11];
smoothMPC_FLOAT smoothMPC_Ysd97[9];
smoothMPC_FLOAT smoothMPC_Lsd97[9];
smoothMPC_FLOAT* smoothMPC_z98 = smoothMPC_z + 1078;
smoothMPC_FLOAT* smoothMPC_dzaff98 = smoothMPC_dz_aff + 1078;
smoothMPC_FLOAT* smoothMPC_dzcc98 = smoothMPC_dz_cc + 1078;
smoothMPC_FLOAT* smoothMPC_rd98 = smoothMPC_rd + 1078;
smoothMPC_FLOAT smoothMPC_Lbyrd98[11];
smoothMPC_FLOAT* smoothMPC_grad_cost98 = smoothMPC_grad_cost + 1078;
smoothMPC_FLOAT* smoothMPC_grad_eq98 = smoothMPC_grad_eq + 1078;
smoothMPC_FLOAT* smoothMPC_grad_ineq98 = smoothMPC_grad_ineq + 1078;
smoothMPC_FLOAT smoothMPC_ctv98[11];
smoothMPC_FLOAT* smoothMPC_v98 = smoothMPC_v + 294;
smoothMPC_FLOAT smoothMPC_re98[3];
smoothMPC_FLOAT smoothMPC_beta98[3];
smoothMPC_FLOAT smoothMPC_betacc98[3];
smoothMPC_FLOAT* smoothMPC_dvaff98 = smoothMPC_dv_aff + 294;
smoothMPC_FLOAT* smoothMPC_dvcc98 = smoothMPC_dv_cc + 294;
smoothMPC_FLOAT smoothMPC_V98[33];
smoothMPC_FLOAT smoothMPC_Yd98[6];
smoothMPC_FLOAT smoothMPC_Ld98[6];
smoothMPC_FLOAT smoothMPC_yy98[3];
smoothMPC_FLOAT smoothMPC_bmy98[3];
int smoothMPC_lbIdx98[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb98 = smoothMPC_l + 1568;
smoothMPC_FLOAT* smoothMPC_slb98 = smoothMPC_s + 1568;
smoothMPC_FLOAT* smoothMPC_llbbyslb98 = smoothMPC_lbys + 1568;
smoothMPC_FLOAT smoothMPC_rilb98[11];
smoothMPC_FLOAT* smoothMPC_dllbaff98 = smoothMPC_dl_aff + 1568;
smoothMPC_FLOAT* smoothMPC_dslbaff98 = smoothMPC_ds_aff + 1568;
smoothMPC_FLOAT* smoothMPC_dllbcc98 = smoothMPC_dl_cc + 1568;
smoothMPC_FLOAT* smoothMPC_dslbcc98 = smoothMPC_ds_cc + 1568;
smoothMPC_FLOAT* smoothMPC_ccrhsl98 = smoothMPC_ccrhs + 1568;
int smoothMPC_ubIdx98[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub98 = smoothMPC_l + 1579;
smoothMPC_FLOAT* smoothMPC_sub98 = smoothMPC_s + 1579;
smoothMPC_FLOAT* smoothMPC_lubbysub98 = smoothMPC_lbys + 1579;
smoothMPC_FLOAT smoothMPC_riub98[5];
smoothMPC_FLOAT* smoothMPC_dlubaff98 = smoothMPC_dl_aff + 1579;
smoothMPC_FLOAT* smoothMPC_dsubaff98 = smoothMPC_ds_aff + 1579;
smoothMPC_FLOAT* smoothMPC_dlubcc98 = smoothMPC_dl_cc + 1579;
smoothMPC_FLOAT* smoothMPC_dsubcc98 = smoothMPC_ds_cc + 1579;
smoothMPC_FLOAT* smoothMPC_ccrhsub98 = smoothMPC_ccrhs + 1579;
smoothMPC_FLOAT smoothMPC_Phi98[11];
smoothMPC_FLOAT smoothMPC_W98[11];
smoothMPC_FLOAT smoothMPC_Ysd98[9];
smoothMPC_FLOAT smoothMPC_Lsd98[9];
smoothMPC_FLOAT* smoothMPC_z99 = smoothMPC_z + 1089;
smoothMPC_FLOAT* smoothMPC_dzaff99 = smoothMPC_dz_aff + 1089;
smoothMPC_FLOAT* smoothMPC_dzcc99 = smoothMPC_dz_cc + 1089;
smoothMPC_FLOAT* smoothMPC_rd99 = smoothMPC_rd + 1089;
smoothMPC_FLOAT smoothMPC_Lbyrd99[3];
smoothMPC_FLOAT* smoothMPC_grad_cost99 = smoothMPC_grad_cost + 1089;
smoothMPC_FLOAT* smoothMPC_grad_eq99 = smoothMPC_grad_eq + 1089;
smoothMPC_FLOAT* smoothMPC_grad_ineq99 = smoothMPC_grad_ineq + 1089;
smoothMPC_FLOAT smoothMPC_ctv99[3];
smoothMPC_FLOAT* smoothMPC_v99 = smoothMPC_v + 297;
smoothMPC_FLOAT smoothMPC_re99[3];
smoothMPC_FLOAT smoothMPC_beta99[3];
smoothMPC_FLOAT smoothMPC_betacc99[3];
smoothMPC_FLOAT* smoothMPC_dvaff99 = smoothMPC_dv_aff + 297;
smoothMPC_FLOAT* smoothMPC_dvcc99 = smoothMPC_dv_cc + 297;
smoothMPC_FLOAT smoothMPC_V99[9];
smoothMPC_FLOAT smoothMPC_Yd99[6];
smoothMPC_FLOAT smoothMPC_Ld99[6];
smoothMPC_FLOAT smoothMPC_yy99[3];
smoothMPC_FLOAT smoothMPC_bmy99[3];
int smoothMPC_lbIdx99[3] = {0, 1, 2};
smoothMPC_FLOAT* smoothMPC_llb99 = smoothMPC_l + 1584;
smoothMPC_FLOAT* smoothMPC_slb99 = smoothMPC_s + 1584;
smoothMPC_FLOAT* smoothMPC_llbbyslb99 = smoothMPC_lbys + 1584;
smoothMPC_FLOAT smoothMPC_rilb99[3];
smoothMPC_FLOAT* smoothMPC_dllbaff99 = smoothMPC_dl_aff + 1584;
smoothMPC_FLOAT* smoothMPC_dslbaff99 = smoothMPC_ds_aff + 1584;
smoothMPC_FLOAT* smoothMPC_dllbcc99 = smoothMPC_dl_cc + 1584;
smoothMPC_FLOAT* smoothMPC_dslbcc99 = smoothMPC_ds_cc + 1584;
smoothMPC_FLOAT* smoothMPC_ccrhsl99 = smoothMPC_ccrhs + 1584;
int smoothMPC_ubIdx99[3] = {0, 1, 2};
smoothMPC_FLOAT* smoothMPC_lub99 = smoothMPC_l + 1587;
smoothMPC_FLOAT* smoothMPC_sub99 = smoothMPC_s + 1587;
smoothMPC_FLOAT* smoothMPC_lubbysub99 = smoothMPC_lbys + 1587;
smoothMPC_FLOAT smoothMPC_riub99[3];
smoothMPC_FLOAT* smoothMPC_dlubaff99 = smoothMPC_dl_aff + 1587;
smoothMPC_FLOAT* smoothMPC_dsubaff99 = smoothMPC_ds_aff + 1587;
smoothMPC_FLOAT* smoothMPC_dlubcc99 = smoothMPC_dl_cc + 1587;
smoothMPC_FLOAT* smoothMPC_dsubcc99 = smoothMPC_ds_cc + 1587;
smoothMPC_FLOAT* smoothMPC_ccrhsub99 = smoothMPC_ccrhs + 1587;
smoothMPC_FLOAT smoothMPC_Phi99[3];
smoothMPC_FLOAT smoothMPC_D99[3] = {-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000};
smoothMPC_FLOAT smoothMPC_W99[3];
smoothMPC_FLOAT smoothMPC_Ysd99[9];
smoothMPC_FLOAT smoothMPC_Lsd99[9];
smoothMPC_FLOAT musigma;
smoothMPC_FLOAT sigma_3rdroot;
smoothMPC_FLOAT smoothMPC_Diag1_0[11];
smoothMPC_FLOAT smoothMPC_Diag2_0[11];
smoothMPC_FLOAT smoothMPC_L_0[55];




/* SOLVER CODE --------------------------------------------------------- */
int smoothMPC_solve(smoothMPC_params* params, smoothMPC_output* output, smoothMPC_info* info)
{	
int exitcode;

#if smoothMPC_SET_TIMING == 1
	smoothMPC_timer solvertimer;
	smoothMPC_tic(&solvertimer);
#endif
/* FUNCTION CALLS INTO LA LIBRARY -------------------------------------- */
info->it = 0;
smoothMPC_LA_INITIALIZEVECTOR_1092(smoothMPC_z, 0);
smoothMPC_LA_INITIALIZEVECTOR_300(smoothMPC_v, 1);
smoothMPC_LA_INITIALIZEVECTOR_1590(smoothMPC_l, 1);
smoothMPC_LA_INITIALIZEVECTOR_1590(smoothMPC_s, 1);
info->mu = 0;
smoothMPC_LA_DOTACC_1590(smoothMPC_l, smoothMPC_s, &info->mu);
info->mu /= 1590;
PRINTTEXT("This is smoothMPC, a solver generated by FORCES (forces.ethz.ch).\n");
PRINTTEXT("(c) Alexander Domahidi, Automatic Control Laboratory, ETH Zurich, 2011-2014.\n");
PRINTTEXT("\n  #it  res_eq   res_ineq     pobj         dobj       dgap     rdgap     mu\n");
PRINTTEXT("  ---------------------------------------------------------------------------\n");
while( 1 ){
info->pobj = 0;
smoothMPC_LA_DIAG_QUADFCN_11(params->H1, params->f1, smoothMPC_z00, smoothMPC_grad_cost00, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H2, params->f2, smoothMPC_z01, smoothMPC_grad_cost01, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H3, params->f3, smoothMPC_z02, smoothMPC_grad_cost02, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H4, params->f4, smoothMPC_z03, smoothMPC_grad_cost03, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H5, params->f5, smoothMPC_z04, smoothMPC_grad_cost04, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H6, params->f6, smoothMPC_z05, smoothMPC_grad_cost05, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H7, params->f7, smoothMPC_z06, smoothMPC_grad_cost06, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H8, params->f8, smoothMPC_z07, smoothMPC_grad_cost07, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H9, params->f9, smoothMPC_z08, smoothMPC_grad_cost08, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H10, params->f10, smoothMPC_z09, smoothMPC_grad_cost09, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H11, params->f11, smoothMPC_z10, smoothMPC_grad_cost10, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H12, params->f12, smoothMPC_z11, smoothMPC_grad_cost11, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H13, params->f13, smoothMPC_z12, smoothMPC_grad_cost12, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H14, params->f14, smoothMPC_z13, smoothMPC_grad_cost13, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H15, params->f15, smoothMPC_z14, smoothMPC_grad_cost14, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H16, params->f16, smoothMPC_z15, smoothMPC_grad_cost15, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H17, params->f17, smoothMPC_z16, smoothMPC_grad_cost16, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H18, params->f18, smoothMPC_z17, smoothMPC_grad_cost17, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H19, params->f19, smoothMPC_z18, smoothMPC_grad_cost18, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H20, params->f20, smoothMPC_z19, smoothMPC_grad_cost19, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H21, params->f21, smoothMPC_z20, smoothMPC_grad_cost20, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H22, params->f22, smoothMPC_z21, smoothMPC_grad_cost21, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H23, params->f23, smoothMPC_z22, smoothMPC_grad_cost22, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H24, params->f24, smoothMPC_z23, smoothMPC_grad_cost23, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H25, params->f25, smoothMPC_z24, smoothMPC_grad_cost24, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H26, params->f26, smoothMPC_z25, smoothMPC_grad_cost25, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H27, params->f27, smoothMPC_z26, smoothMPC_grad_cost26, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H28, params->f28, smoothMPC_z27, smoothMPC_grad_cost27, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H29, params->f29, smoothMPC_z28, smoothMPC_grad_cost28, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H30, params->f30, smoothMPC_z29, smoothMPC_grad_cost29, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H31, params->f31, smoothMPC_z30, smoothMPC_grad_cost30, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H32, params->f32, smoothMPC_z31, smoothMPC_grad_cost31, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H33, params->f33, smoothMPC_z32, smoothMPC_grad_cost32, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H34, params->f34, smoothMPC_z33, smoothMPC_grad_cost33, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H35, params->f35, smoothMPC_z34, smoothMPC_grad_cost34, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H36, params->f36, smoothMPC_z35, smoothMPC_grad_cost35, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H37, params->f37, smoothMPC_z36, smoothMPC_grad_cost36, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H38, params->f38, smoothMPC_z37, smoothMPC_grad_cost37, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H39, params->f39, smoothMPC_z38, smoothMPC_grad_cost38, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H40, params->f40, smoothMPC_z39, smoothMPC_grad_cost39, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H41, params->f41, smoothMPC_z40, smoothMPC_grad_cost40, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H42, params->f42, smoothMPC_z41, smoothMPC_grad_cost41, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H43, params->f43, smoothMPC_z42, smoothMPC_grad_cost42, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H44, params->f44, smoothMPC_z43, smoothMPC_grad_cost43, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H45, params->f45, smoothMPC_z44, smoothMPC_grad_cost44, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H46, params->f46, smoothMPC_z45, smoothMPC_grad_cost45, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H47, params->f47, smoothMPC_z46, smoothMPC_grad_cost46, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H48, params->f48, smoothMPC_z47, smoothMPC_grad_cost47, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H49, params->f49, smoothMPC_z48, smoothMPC_grad_cost48, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H50, params->f50, smoothMPC_z49, smoothMPC_grad_cost49, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H51, params->f51, smoothMPC_z50, smoothMPC_grad_cost50, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H52, params->f52, smoothMPC_z51, smoothMPC_grad_cost51, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H53, params->f53, smoothMPC_z52, smoothMPC_grad_cost52, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H54, params->f54, smoothMPC_z53, smoothMPC_grad_cost53, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H55, params->f55, smoothMPC_z54, smoothMPC_grad_cost54, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H56, params->f56, smoothMPC_z55, smoothMPC_grad_cost55, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H57, params->f57, smoothMPC_z56, smoothMPC_grad_cost56, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H58, params->f58, smoothMPC_z57, smoothMPC_grad_cost57, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H59, params->f59, smoothMPC_z58, smoothMPC_grad_cost58, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H60, params->f60, smoothMPC_z59, smoothMPC_grad_cost59, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H61, params->f61, smoothMPC_z60, smoothMPC_grad_cost60, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H62, params->f62, smoothMPC_z61, smoothMPC_grad_cost61, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H63, params->f63, smoothMPC_z62, smoothMPC_grad_cost62, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H64, params->f64, smoothMPC_z63, smoothMPC_grad_cost63, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H65, params->f65, smoothMPC_z64, smoothMPC_grad_cost64, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H66, params->f66, smoothMPC_z65, smoothMPC_grad_cost65, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H67, params->f67, smoothMPC_z66, smoothMPC_grad_cost66, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H68, params->f68, smoothMPC_z67, smoothMPC_grad_cost67, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H69, params->f69, smoothMPC_z68, smoothMPC_grad_cost68, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H70, params->f70, smoothMPC_z69, smoothMPC_grad_cost69, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H71, params->f71, smoothMPC_z70, smoothMPC_grad_cost70, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H72, params->f72, smoothMPC_z71, smoothMPC_grad_cost71, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H73, params->f73, smoothMPC_z72, smoothMPC_grad_cost72, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H74, params->f74, smoothMPC_z73, smoothMPC_grad_cost73, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H75, params->f75, smoothMPC_z74, smoothMPC_grad_cost74, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H76, params->f76, smoothMPC_z75, smoothMPC_grad_cost75, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H77, params->f77, smoothMPC_z76, smoothMPC_grad_cost76, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H78, params->f78, smoothMPC_z77, smoothMPC_grad_cost77, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H79, params->f79, smoothMPC_z78, smoothMPC_grad_cost78, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H80, params->f80, smoothMPC_z79, smoothMPC_grad_cost79, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H81, params->f81, smoothMPC_z80, smoothMPC_grad_cost80, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H82, params->f82, smoothMPC_z81, smoothMPC_grad_cost81, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H83, params->f83, smoothMPC_z82, smoothMPC_grad_cost82, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H84, params->f84, smoothMPC_z83, smoothMPC_grad_cost83, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H85, params->f85, smoothMPC_z84, smoothMPC_grad_cost84, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H86, params->f86, smoothMPC_z85, smoothMPC_grad_cost85, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H87, params->f87, smoothMPC_z86, smoothMPC_grad_cost86, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H88, params->f88, smoothMPC_z87, smoothMPC_grad_cost87, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H89, params->f89, smoothMPC_z88, smoothMPC_grad_cost88, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H90, params->f90, smoothMPC_z89, smoothMPC_grad_cost89, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H91, params->f91, smoothMPC_z90, smoothMPC_grad_cost90, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H92, params->f92, smoothMPC_z91, smoothMPC_grad_cost91, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H93, params->f93, smoothMPC_z92, smoothMPC_grad_cost92, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H94, params->f94, smoothMPC_z93, smoothMPC_grad_cost93, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H95, params->f95, smoothMPC_z94, smoothMPC_grad_cost94, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H96, params->f96, smoothMPC_z95, smoothMPC_grad_cost95, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H97, params->f97, smoothMPC_z96, smoothMPC_grad_cost96, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H98, params->f98, smoothMPC_z97, smoothMPC_grad_cost97, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H99, params->f99, smoothMPC_z98, smoothMPC_grad_cost98, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_3(params->H100, params->f100, smoothMPC_z99, smoothMPC_grad_cost99, &info->pobj);
info->res_eq = 0;
info->dgap = 0;
smoothMPC_LA_DIAGZERO_MVMSUB6_3(smoothMPC_D00, smoothMPC_z00, params->e1, smoothMPC_v00, smoothMPC_re00, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C1, smoothMPC_z00, smoothMPC_D01, smoothMPC_z01, params->e2, smoothMPC_v01, smoothMPC_re01, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C2, smoothMPC_z01, smoothMPC_D01, smoothMPC_z02, params->e3, smoothMPC_v02, smoothMPC_re02, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C3, smoothMPC_z02, smoothMPC_D01, smoothMPC_z03, params->e4, smoothMPC_v03, smoothMPC_re03, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C4, smoothMPC_z03, smoothMPC_D01, smoothMPC_z04, params->e5, smoothMPC_v04, smoothMPC_re04, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C5, smoothMPC_z04, smoothMPC_D01, smoothMPC_z05, params->e6, smoothMPC_v05, smoothMPC_re05, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C6, smoothMPC_z05, smoothMPC_D01, smoothMPC_z06, params->e7, smoothMPC_v06, smoothMPC_re06, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C7, smoothMPC_z06, smoothMPC_D01, smoothMPC_z07, params->e8, smoothMPC_v07, smoothMPC_re07, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C8, smoothMPC_z07, smoothMPC_D01, smoothMPC_z08, params->e9, smoothMPC_v08, smoothMPC_re08, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C9, smoothMPC_z08, smoothMPC_D01, smoothMPC_z09, params->e10, smoothMPC_v09, smoothMPC_re09, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C10, smoothMPC_z09, smoothMPC_D01, smoothMPC_z10, params->e11, smoothMPC_v10, smoothMPC_re10, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C11, smoothMPC_z10, smoothMPC_D01, smoothMPC_z11, params->e12, smoothMPC_v11, smoothMPC_re11, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C12, smoothMPC_z11, smoothMPC_D01, smoothMPC_z12, params->e13, smoothMPC_v12, smoothMPC_re12, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C13, smoothMPC_z12, smoothMPC_D01, smoothMPC_z13, params->e14, smoothMPC_v13, smoothMPC_re13, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C14, smoothMPC_z13, smoothMPC_D01, smoothMPC_z14, params->e15, smoothMPC_v14, smoothMPC_re14, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C15, smoothMPC_z14, smoothMPC_D01, smoothMPC_z15, params->e16, smoothMPC_v15, smoothMPC_re15, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C16, smoothMPC_z15, smoothMPC_D01, smoothMPC_z16, params->e17, smoothMPC_v16, smoothMPC_re16, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C17, smoothMPC_z16, smoothMPC_D01, smoothMPC_z17, params->e18, smoothMPC_v17, smoothMPC_re17, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C18, smoothMPC_z17, smoothMPC_D01, smoothMPC_z18, params->e19, smoothMPC_v18, smoothMPC_re18, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C19, smoothMPC_z18, smoothMPC_D01, smoothMPC_z19, params->e20, smoothMPC_v19, smoothMPC_re19, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C20, smoothMPC_z19, smoothMPC_D01, smoothMPC_z20, params->e21, smoothMPC_v20, smoothMPC_re20, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C21, smoothMPC_z20, smoothMPC_D01, smoothMPC_z21, params->e22, smoothMPC_v21, smoothMPC_re21, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C22, smoothMPC_z21, smoothMPC_D01, smoothMPC_z22, params->e23, smoothMPC_v22, smoothMPC_re22, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C23, smoothMPC_z22, smoothMPC_D01, smoothMPC_z23, params->e24, smoothMPC_v23, smoothMPC_re23, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C24, smoothMPC_z23, smoothMPC_D01, smoothMPC_z24, params->e25, smoothMPC_v24, smoothMPC_re24, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C25, smoothMPC_z24, smoothMPC_D01, smoothMPC_z25, params->e26, smoothMPC_v25, smoothMPC_re25, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C26, smoothMPC_z25, smoothMPC_D01, smoothMPC_z26, params->e27, smoothMPC_v26, smoothMPC_re26, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C27, smoothMPC_z26, smoothMPC_D01, smoothMPC_z27, params->e28, smoothMPC_v27, smoothMPC_re27, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C28, smoothMPC_z27, smoothMPC_D01, smoothMPC_z28, params->e29, smoothMPC_v28, smoothMPC_re28, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C29, smoothMPC_z28, smoothMPC_D01, smoothMPC_z29, params->e30, smoothMPC_v29, smoothMPC_re29, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C30, smoothMPC_z29, smoothMPC_D01, smoothMPC_z30, params->e31, smoothMPC_v30, smoothMPC_re30, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C31, smoothMPC_z30, smoothMPC_D01, smoothMPC_z31, params->e32, smoothMPC_v31, smoothMPC_re31, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C32, smoothMPC_z31, smoothMPC_D01, smoothMPC_z32, params->e33, smoothMPC_v32, smoothMPC_re32, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C33, smoothMPC_z32, smoothMPC_D01, smoothMPC_z33, params->e34, smoothMPC_v33, smoothMPC_re33, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C34, smoothMPC_z33, smoothMPC_D01, smoothMPC_z34, params->e35, smoothMPC_v34, smoothMPC_re34, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C35, smoothMPC_z34, smoothMPC_D01, smoothMPC_z35, params->e36, smoothMPC_v35, smoothMPC_re35, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C36, smoothMPC_z35, smoothMPC_D01, smoothMPC_z36, params->e37, smoothMPC_v36, smoothMPC_re36, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C37, smoothMPC_z36, smoothMPC_D01, smoothMPC_z37, params->e38, smoothMPC_v37, smoothMPC_re37, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C38, smoothMPC_z37, smoothMPC_D01, smoothMPC_z38, params->e39, smoothMPC_v38, smoothMPC_re38, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C39, smoothMPC_z38, smoothMPC_D01, smoothMPC_z39, params->e40, smoothMPC_v39, smoothMPC_re39, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C40, smoothMPC_z39, smoothMPC_D01, smoothMPC_z40, params->e41, smoothMPC_v40, smoothMPC_re40, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C41, smoothMPC_z40, smoothMPC_D01, smoothMPC_z41, params->e42, smoothMPC_v41, smoothMPC_re41, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C42, smoothMPC_z41, smoothMPC_D01, smoothMPC_z42, params->e43, smoothMPC_v42, smoothMPC_re42, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C43, smoothMPC_z42, smoothMPC_D01, smoothMPC_z43, params->e44, smoothMPC_v43, smoothMPC_re43, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C44, smoothMPC_z43, smoothMPC_D01, smoothMPC_z44, params->e45, smoothMPC_v44, smoothMPC_re44, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C45, smoothMPC_z44, smoothMPC_D01, smoothMPC_z45, params->e46, smoothMPC_v45, smoothMPC_re45, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C46, smoothMPC_z45, smoothMPC_D01, smoothMPC_z46, params->e47, smoothMPC_v46, smoothMPC_re46, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C47, smoothMPC_z46, smoothMPC_D01, smoothMPC_z47, params->e48, smoothMPC_v47, smoothMPC_re47, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C48, smoothMPC_z47, smoothMPC_D01, smoothMPC_z48, params->e49, smoothMPC_v48, smoothMPC_re48, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C49, smoothMPC_z48, smoothMPC_D01, smoothMPC_z49, params->e50, smoothMPC_v49, smoothMPC_re49, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C50, smoothMPC_z49, smoothMPC_D01, smoothMPC_z50, params->e51, smoothMPC_v50, smoothMPC_re50, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C51, smoothMPC_z50, smoothMPC_D01, smoothMPC_z51, params->e52, smoothMPC_v51, smoothMPC_re51, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C52, smoothMPC_z51, smoothMPC_D01, smoothMPC_z52, params->e53, smoothMPC_v52, smoothMPC_re52, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C53, smoothMPC_z52, smoothMPC_D01, smoothMPC_z53, params->e54, smoothMPC_v53, smoothMPC_re53, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C54, smoothMPC_z53, smoothMPC_D01, smoothMPC_z54, params->e55, smoothMPC_v54, smoothMPC_re54, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C55, smoothMPC_z54, smoothMPC_D01, smoothMPC_z55, params->e56, smoothMPC_v55, smoothMPC_re55, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C56, smoothMPC_z55, smoothMPC_D01, smoothMPC_z56, params->e57, smoothMPC_v56, smoothMPC_re56, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C57, smoothMPC_z56, smoothMPC_D01, smoothMPC_z57, params->e58, smoothMPC_v57, smoothMPC_re57, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C58, smoothMPC_z57, smoothMPC_D01, smoothMPC_z58, params->e59, smoothMPC_v58, smoothMPC_re58, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C59, smoothMPC_z58, smoothMPC_D01, smoothMPC_z59, params->e60, smoothMPC_v59, smoothMPC_re59, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C60, smoothMPC_z59, smoothMPC_D01, smoothMPC_z60, params->e61, smoothMPC_v60, smoothMPC_re60, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C61, smoothMPC_z60, smoothMPC_D01, smoothMPC_z61, params->e62, smoothMPC_v61, smoothMPC_re61, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C62, smoothMPC_z61, smoothMPC_D01, smoothMPC_z62, params->e63, smoothMPC_v62, smoothMPC_re62, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C63, smoothMPC_z62, smoothMPC_D01, smoothMPC_z63, params->e64, smoothMPC_v63, smoothMPC_re63, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C64, smoothMPC_z63, smoothMPC_D01, smoothMPC_z64, params->e65, smoothMPC_v64, smoothMPC_re64, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C65, smoothMPC_z64, smoothMPC_D01, smoothMPC_z65, params->e66, smoothMPC_v65, smoothMPC_re65, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C66, smoothMPC_z65, smoothMPC_D01, smoothMPC_z66, params->e67, smoothMPC_v66, smoothMPC_re66, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C67, smoothMPC_z66, smoothMPC_D01, smoothMPC_z67, params->e68, smoothMPC_v67, smoothMPC_re67, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C68, smoothMPC_z67, smoothMPC_D01, smoothMPC_z68, params->e69, smoothMPC_v68, smoothMPC_re68, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C69, smoothMPC_z68, smoothMPC_D01, smoothMPC_z69, params->e70, smoothMPC_v69, smoothMPC_re69, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C70, smoothMPC_z69, smoothMPC_D01, smoothMPC_z70, params->e71, smoothMPC_v70, smoothMPC_re70, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C71, smoothMPC_z70, smoothMPC_D01, smoothMPC_z71, params->e72, smoothMPC_v71, smoothMPC_re71, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C72, smoothMPC_z71, smoothMPC_D01, smoothMPC_z72, params->e73, smoothMPC_v72, smoothMPC_re72, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C73, smoothMPC_z72, smoothMPC_D01, smoothMPC_z73, params->e74, smoothMPC_v73, smoothMPC_re73, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C74, smoothMPC_z73, smoothMPC_D01, smoothMPC_z74, params->e75, smoothMPC_v74, smoothMPC_re74, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C75, smoothMPC_z74, smoothMPC_D01, smoothMPC_z75, params->e76, smoothMPC_v75, smoothMPC_re75, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C76, smoothMPC_z75, smoothMPC_D01, smoothMPC_z76, params->e77, smoothMPC_v76, smoothMPC_re76, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C77, smoothMPC_z76, smoothMPC_D01, smoothMPC_z77, params->e78, smoothMPC_v77, smoothMPC_re77, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C78, smoothMPC_z77, smoothMPC_D01, smoothMPC_z78, params->e79, smoothMPC_v78, smoothMPC_re78, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C79, smoothMPC_z78, smoothMPC_D01, smoothMPC_z79, params->e80, smoothMPC_v79, smoothMPC_re79, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C80, smoothMPC_z79, smoothMPC_D01, smoothMPC_z80, params->e81, smoothMPC_v80, smoothMPC_re80, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C81, smoothMPC_z80, smoothMPC_D01, smoothMPC_z81, params->e82, smoothMPC_v81, smoothMPC_re81, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C82, smoothMPC_z81, smoothMPC_D01, smoothMPC_z82, params->e83, smoothMPC_v82, smoothMPC_re82, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C83, smoothMPC_z82, smoothMPC_D01, smoothMPC_z83, params->e84, smoothMPC_v83, smoothMPC_re83, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C84, smoothMPC_z83, smoothMPC_D01, smoothMPC_z84, params->e85, smoothMPC_v84, smoothMPC_re84, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C85, smoothMPC_z84, smoothMPC_D01, smoothMPC_z85, params->e86, smoothMPC_v85, smoothMPC_re85, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C86, smoothMPC_z85, smoothMPC_D01, smoothMPC_z86, params->e87, smoothMPC_v86, smoothMPC_re86, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C87, smoothMPC_z86, smoothMPC_D01, smoothMPC_z87, params->e88, smoothMPC_v87, smoothMPC_re87, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C88, smoothMPC_z87, smoothMPC_D01, smoothMPC_z88, params->e89, smoothMPC_v88, smoothMPC_re88, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C89, smoothMPC_z88, smoothMPC_D01, smoothMPC_z89, params->e90, smoothMPC_v89, smoothMPC_re89, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C90, smoothMPC_z89, smoothMPC_D01, smoothMPC_z90, params->e91, smoothMPC_v90, smoothMPC_re90, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C91, smoothMPC_z90, smoothMPC_D01, smoothMPC_z91, params->e92, smoothMPC_v91, smoothMPC_re91, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C92, smoothMPC_z91, smoothMPC_D01, smoothMPC_z92, params->e93, smoothMPC_v92, smoothMPC_re92, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C93, smoothMPC_z92, smoothMPC_D01, smoothMPC_z93, params->e94, smoothMPC_v93, smoothMPC_re93, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C94, smoothMPC_z93, smoothMPC_D01, smoothMPC_z94, params->e95, smoothMPC_v94, smoothMPC_re94, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C95, smoothMPC_z94, smoothMPC_D01, smoothMPC_z95, params->e96, smoothMPC_v95, smoothMPC_re95, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C96, smoothMPC_z95, smoothMPC_D01, smoothMPC_z96, params->e97, smoothMPC_v96, smoothMPC_re96, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C97, smoothMPC_z96, smoothMPC_D01, smoothMPC_z97, params->e98, smoothMPC_v97, smoothMPC_re97, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C98, smoothMPC_z97, smoothMPC_D01, smoothMPC_z98, params->e99, smoothMPC_v98, smoothMPC_re98, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_3(params->C99, smoothMPC_z98, smoothMPC_D99, smoothMPC_z99, params->e100, smoothMPC_v99, smoothMPC_re99, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C1, smoothMPC_v01, smoothMPC_D00, smoothMPC_v00, smoothMPC_grad_eq00);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C2, smoothMPC_v02, smoothMPC_D01, smoothMPC_v01, smoothMPC_grad_eq01);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C3, smoothMPC_v03, smoothMPC_D01, smoothMPC_v02, smoothMPC_grad_eq02);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C4, smoothMPC_v04, smoothMPC_D01, smoothMPC_v03, smoothMPC_grad_eq03);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C5, smoothMPC_v05, smoothMPC_D01, smoothMPC_v04, smoothMPC_grad_eq04);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C6, smoothMPC_v06, smoothMPC_D01, smoothMPC_v05, smoothMPC_grad_eq05);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C7, smoothMPC_v07, smoothMPC_D01, smoothMPC_v06, smoothMPC_grad_eq06);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C8, smoothMPC_v08, smoothMPC_D01, smoothMPC_v07, smoothMPC_grad_eq07);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C9, smoothMPC_v09, smoothMPC_D01, smoothMPC_v08, smoothMPC_grad_eq08);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C10, smoothMPC_v10, smoothMPC_D01, smoothMPC_v09, smoothMPC_grad_eq09);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C11, smoothMPC_v11, smoothMPC_D01, smoothMPC_v10, smoothMPC_grad_eq10);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C12, smoothMPC_v12, smoothMPC_D01, smoothMPC_v11, smoothMPC_grad_eq11);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C13, smoothMPC_v13, smoothMPC_D01, smoothMPC_v12, smoothMPC_grad_eq12);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C14, smoothMPC_v14, smoothMPC_D01, smoothMPC_v13, smoothMPC_grad_eq13);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C15, smoothMPC_v15, smoothMPC_D01, smoothMPC_v14, smoothMPC_grad_eq14);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C16, smoothMPC_v16, smoothMPC_D01, smoothMPC_v15, smoothMPC_grad_eq15);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C17, smoothMPC_v17, smoothMPC_D01, smoothMPC_v16, smoothMPC_grad_eq16);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C18, smoothMPC_v18, smoothMPC_D01, smoothMPC_v17, smoothMPC_grad_eq17);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C19, smoothMPC_v19, smoothMPC_D01, smoothMPC_v18, smoothMPC_grad_eq18);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C20, smoothMPC_v20, smoothMPC_D01, smoothMPC_v19, smoothMPC_grad_eq19);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C21, smoothMPC_v21, smoothMPC_D01, smoothMPC_v20, smoothMPC_grad_eq20);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C22, smoothMPC_v22, smoothMPC_D01, smoothMPC_v21, smoothMPC_grad_eq21);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C23, smoothMPC_v23, smoothMPC_D01, smoothMPC_v22, smoothMPC_grad_eq22);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C24, smoothMPC_v24, smoothMPC_D01, smoothMPC_v23, smoothMPC_grad_eq23);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C25, smoothMPC_v25, smoothMPC_D01, smoothMPC_v24, smoothMPC_grad_eq24);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C26, smoothMPC_v26, smoothMPC_D01, smoothMPC_v25, smoothMPC_grad_eq25);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C27, smoothMPC_v27, smoothMPC_D01, smoothMPC_v26, smoothMPC_grad_eq26);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C28, smoothMPC_v28, smoothMPC_D01, smoothMPC_v27, smoothMPC_grad_eq27);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C29, smoothMPC_v29, smoothMPC_D01, smoothMPC_v28, smoothMPC_grad_eq28);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C30, smoothMPC_v30, smoothMPC_D01, smoothMPC_v29, smoothMPC_grad_eq29);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C31, smoothMPC_v31, smoothMPC_D01, smoothMPC_v30, smoothMPC_grad_eq30);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C32, smoothMPC_v32, smoothMPC_D01, smoothMPC_v31, smoothMPC_grad_eq31);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C33, smoothMPC_v33, smoothMPC_D01, smoothMPC_v32, smoothMPC_grad_eq32);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C34, smoothMPC_v34, smoothMPC_D01, smoothMPC_v33, smoothMPC_grad_eq33);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C35, smoothMPC_v35, smoothMPC_D01, smoothMPC_v34, smoothMPC_grad_eq34);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C36, smoothMPC_v36, smoothMPC_D01, smoothMPC_v35, smoothMPC_grad_eq35);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C37, smoothMPC_v37, smoothMPC_D01, smoothMPC_v36, smoothMPC_grad_eq36);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C38, smoothMPC_v38, smoothMPC_D01, smoothMPC_v37, smoothMPC_grad_eq37);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C39, smoothMPC_v39, smoothMPC_D01, smoothMPC_v38, smoothMPC_grad_eq38);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C40, smoothMPC_v40, smoothMPC_D01, smoothMPC_v39, smoothMPC_grad_eq39);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C41, smoothMPC_v41, smoothMPC_D01, smoothMPC_v40, smoothMPC_grad_eq40);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C42, smoothMPC_v42, smoothMPC_D01, smoothMPC_v41, smoothMPC_grad_eq41);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C43, smoothMPC_v43, smoothMPC_D01, smoothMPC_v42, smoothMPC_grad_eq42);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C44, smoothMPC_v44, smoothMPC_D01, smoothMPC_v43, smoothMPC_grad_eq43);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C45, smoothMPC_v45, smoothMPC_D01, smoothMPC_v44, smoothMPC_grad_eq44);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C46, smoothMPC_v46, smoothMPC_D01, smoothMPC_v45, smoothMPC_grad_eq45);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C47, smoothMPC_v47, smoothMPC_D01, smoothMPC_v46, smoothMPC_grad_eq46);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C48, smoothMPC_v48, smoothMPC_D01, smoothMPC_v47, smoothMPC_grad_eq47);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C49, smoothMPC_v49, smoothMPC_D01, smoothMPC_v48, smoothMPC_grad_eq48);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C50, smoothMPC_v50, smoothMPC_D01, smoothMPC_v49, smoothMPC_grad_eq49);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C51, smoothMPC_v51, smoothMPC_D01, smoothMPC_v50, smoothMPC_grad_eq50);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C52, smoothMPC_v52, smoothMPC_D01, smoothMPC_v51, smoothMPC_grad_eq51);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C53, smoothMPC_v53, smoothMPC_D01, smoothMPC_v52, smoothMPC_grad_eq52);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C54, smoothMPC_v54, smoothMPC_D01, smoothMPC_v53, smoothMPC_grad_eq53);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C55, smoothMPC_v55, smoothMPC_D01, smoothMPC_v54, smoothMPC_grad_eq54);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C56, smoothMPC_v56, smoothMPC_D01, smoothMPC_v55, smoothMPC_grad_eq55);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C57, smoothMPC_v57, smoothMPC_D01, smoothMPC_v56, smoothMPC_grad_eq56);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C58, smoothMPC_v58, smoothMPC_D01, smoothMPC_v57, smoothMPC_grad_eq57);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C59, smoothMPC_v59, smoothMPC_D01, smoothMPC_v58, smoothMPC_grad_eq58);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C60, smoothMPC_v60, smoothMPC_D01, smoothMPC_v59, smoothMPC_grad_eq59);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C61, smoothMPC_v61, smoothMPC_D01, smoothMPC_v60, smoothMPC_grad_eq60);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C62, smoothMPC_v62, smoothMPC_D01, smoothMPC_v61, smoothMPC_grad_eq61);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C63, smoothMPC_v63, smoothMPC_D01, smoothMPC_v62, smoothMPC_grad_eq62);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C64, smoothMPC_v64, smoothMPC_D01, smoothMPC_v63, smoothMPC_grad_eq63);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C65, smoothMPC_v65, smoothMPC_D01, smoothMPC_v64, smoothMPC_grad_eq64);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C66, smoothMPC_v66, smoothMPC_D01, smoothMPC_v65, smoothMPC_grad_eq65);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C67, smoothMPC_v67, smoothMPC_D01, smoothMPC_v66, smoothMPC_grad_eq66);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C68, smoothMPC_v68, smoothMPC_D01, smoothMPC_v67, smoothMPC_grad_eq67);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C69, smoothMPC_v69, smoothMPC_D01, smoothMPC_v68, smoothMPC_grad_eq68);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C70, smoothMPC_v70, smoothMPC_D01, smoothMPC_v69, smoothMPC_grad_eq69);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C71, smoothMPC_v71, smoothMPC_D01, smoothMPC_v70, smoothMPC_grad_eq70);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C72, smoothMPC_v72, smoothMPC_D01, smoothMPC_v71, smoothMPC_grad_eq71);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C73, smoothMPC_v73, smoothMPC_D01, smoothMPC_v72, smoothMPC_grad_eq72);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C74, smoothMPC_v74, smoothMPC_D01, smoothMPC_v73, smoothMPC_grad_eq73);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C75, smoothMPC_v75, smoothMPC_D01, smoothMPC_v74, smoothMPC_grad_eq74);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C76, smoothMPC_v76, smoothMPC_D01, smoothMPC_v75, smoothMPC_grad_eq75);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C77, smoothMPC_v77, smoothMPC_D01, smoothMPC_v76, smoothMPC_grad_eq76);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C78, smoothMPC_v78, smoothMPC_D01, smoothMPC_v77, smoothMPC_grad_eq77);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C79, smoothMPC_v79, smoothMPC_D01, smoothMPC_v78, smoothMPC_grad_eq78);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C80, smoothMPC_v80, smoothMPC_D01, smoothMPC_v79, smoothMPC_grad_eq79);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C81, smoothMPC_v81, smoothMPC_D01, smoothMPC_v80, smoothMPC_grad_eq80);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C82, smoothMPC_v82, smoothMPC_D01, smoothMPC_v81, smoothMPC_grad_eq81);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C83, smoothMPC_v83, smoothMPC_D01, smoothMPC_v82, smoothMPC_grad_eq82);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C84, smoothMPC_v84, smoothMPC_D01, smoothMPC_v83, smoothMPC_grad_eq83);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C85, smoothMPC_v85, smoothMPC_D01, smoothMPC_v84, smoothMPC_grad_eq84);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C86, smoothMPC_v86, smoothMPC_D01, smoothMPC_v85, smoothMPC_grad_eq85);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C87, smoothMPC_v87, smoothMPC_D01, smoothMPC_v86, smoothMPC_grad_eq86);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C88, smoothMPC_v88, smoothMPC_D01, smoothMPC_v87, smoothMPC_grad_eq87);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C89, smoothMPC_v89, smoothMPC_D01, smoothMPC_v88, smoothMPC_grad_eq88);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C90, smoothMPC_v90, smoothMPC_D01, smoothMPC_v89, smoothMPC_grad_eq89);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C91, smoothMPC_v91, smoothMPC_D01, smoothMPC_v90, smoothMPC_grad_eq90);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C92, smoothMPC_v92, smoothMPC_D01, smoothMPC_v91, smoothMPC_grad_eq91);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C93, smoothMPC_v93, smoothMPC_D01, smoothMPC_v92, smoothMPC_grad_eq92);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C94, smoothMPC_v94, smoothMPC_D01, smoothMPC_v93, smoothMPC_grad_eq93);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C95, smoothMPC_v95, smoothMPC_D01, smoothMPC_v94, smoothMPC_grad_eq94);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C96, smoothMPC_v96, smoothMPC_D01, smoothMPC_v95, smoothMPC_grad_eq95);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C97, smoothMPC_v97, smoothMPC_D01, smoothMPC_v96, smoothMPC_grad_eq96);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C98, smoothMPC_v98, smoothMPC_D01, smoothMPC_v97, smoothMPC_grad_eq97);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C99, smoothMPC_v99, smoothMPC_D01, smoothMPC_v98, smoothMPC_grad_eq98);
smoothMPC_LA_DIAGZERO_MTVM_3_3(smoothMPC_D99, smoothMPC_v99, smoothMPC_grad_eq99);
info->res_ineq = 0;
smoothMPC_LA_VSUBADD3_11(params->lb1, smoothMPC_z00, smoothMPC_lbIdx00, smoothMPC_llb00, smoothMPC_slb00, smoothMPC_rilb00, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z00, smoothMPC_ubIdx00, params->ub1, smoothMPC_lub00, smoothMPC_sub00, smoothMPC_riub00, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb2, smoothMPC_z01, smoothMPC_lbIdx01, smoothMPC_llb01, smoothMPC_slb01, smoothMPC_rilb01, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z01, smoothMPC_ubIdx01, params->ub2, smoothMPC_lub01, smoothMPC_sub01, smoothMPC_riub01, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb3, smoothMPC_z02, smoothMPC_lbIdx02, smoothMPC_llb02, smoothMPC_slb02, smoothMPC_rilb02, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z02, smoothMPC_ubIdx02, params->ub3, smoothMPC_lub02, smoothMPC_sub02, smoothMPC_riub02, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb4, smoothMPC_z03, smoothMPC_lbIdx03, smoothMPC_llb03, smoothMPC_slb03, smoothMPC_rilb03, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z03, smoothMPC_ubIdx03, params->ub4, smoothMPC_lub03, smoothMPC_sub03, smoothMPC_riub03, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb5, smoothMPC_z04, smoothMPC_lbIdx04, smoothMPC_llb04, smoothMPC_slb04, smoothMPC_rilb04, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z04, smoothMPC_ubIdx04, params->ub5, smoothMPC_lub04, smoothMPC_sub04, smoothMPC_riub04, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb6, smoothMPC_z05, smoothMPC_lbIdx05, smoothMPC_llb05, smoothMPC_slb05, smoothMPC_rilb05, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z05, smoothMPC_ubIdx05, params->ub6, smoothMPC_lub05, smoothMPC_sub05, smoothMPC_riub05, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb7, smoothMPC_z06, smoothMPC_lbIdx06, smoothMPC_llb06, smoothMPC_slb06, smoothMPC_rilb06, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z06, smoothMPC_ubIdx06, params->ub7, smoothMPC_lub06, smoothMPC_sub06, smoothMPC_riub06, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb8, smoothMPC_z07, smoothMPC_lbIdx07, smoothMPC_llb07, smoothMPC_slb07, smoothMPC_rilb07, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z07, smoothMPC_ubIdx07, params->ub8, smoothMPC_lub07, smoothMPC_sub07, smoothMPC_riub07, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb9, smoothMPC_z08, smoothMPC_lbIdx08, smoothMPC_llb08, smoothMPC_slb08, smoothMPC_rilb08, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z08, smoothMPC_ubIdx08, params->ub9, smoothMPC_lub08, smoothMPC_sub08, smoothMPC_riub08, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb10, smoothMPC_z09, smoothMPC_lbIdx09, smoothMPC_llb09, smoothMPC_slb09, smoothMPC_rilb09, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z09, smoothMPC_ubIdx09, params->ub10, smoothMPC_lub09, smoothMPC_sub09, smoothMPC_riub09, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb11, smoothMPC_z10, smoothMPC_lbIdx10, smoothMPC_llb10, smoothMPC_slb10, smoothMPC_rilb10, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z10, smoothMPC_ubIdx10, params->ub11, smoothMPC_lub10, smoothMPC_sub10, smoothMPC_riub10, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb12, smoothMPC_z11, smoothMPC_lbIdx11, smoothMPC_llb11, smoothMPC_slb11, smoothMPC_rilb11, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z11, smoothMPC_ubIdx11, params->ub12, smoothMPC_lub11, smoothMPC_sub11, smoothMPC_riub11, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb13, smoothMPC_z12, smoothMPC_lbIdx12, smoothMPC_llb12, smoothMPC_slb12, smoothMPC_rilb12, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z12, smoothMPC_ubIdx12, params->ub13, smoothMPC_lub12, smoothMPC_sub12, smoothMPC_riub12, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb14, smoothMPC_z13, smoothMPC_lbIdx13, smoothMPC_llb13, smoothMPC_slb13, smoothMPC_rilb13, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z13, smoothMPC_ubIdx13, params->ub14, smoothMPC_lub13, smoothMPC_sub13, smoothMPC_riub13, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb15, smoothMPC_z14, smoothMPC_lbIdx14, smoothMPC_llb14, smoothMPC_slb14, smoothMPC_rilb14, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z14, smoothMPC_ubIdx14, params->ub15, smoothMPC_lub14, smoothMPC_sub14, smoothMPC_riub14, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb16, smoothMPC_z15, smoothMPC_lbIdx15, smoothMPC_llb15, smoothMPC_slb15, smoothMPC_rilb15, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z15, smoothMPC_ubIdx15, params->ub16, smoothMPC_lub15, smoothMPC_sub15, smoothMPC_riub15, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb17, smoothMPC_z16, smoothMPC_lbIdx16, smoothMPC_llb16, smoothMPC_slb16, smoothMPC_rilb16, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z16, smoothMPC_ubIdx16, params->ub17, smoothMPC_lub16, smoothMPC_sub16, smoothMPC_riub16, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb18, smoothMPC_z17, smoothMPC_lbIdx17, smoothMPC_llb17, smoothMPC_slb17, smoothMPC_rilb17, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z17, smoothMPC_ubIdx17, params->ub18, smoothMPC_lub17, smoothMPC_sub17, smoothMPC_riub17, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb19, smoothMPC_z18, smoothMPC_lbIdx18, smoothMPC_llb18, smoothMPC_slb18, smoothMPC_rilb18, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z18, smoothMPC_ubIdx18, params->ub19, smoothMPC_lub18, smoothMPC_sub18, smoothMPC_riub18, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb20, smoothMPC_z19, smoothMPC_lbIdx19, smoothMPC_llb19, smoothMPC_slb19, smoothMPC_rilb19, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z19, smoothMPC_ubIdx19, params->ub20, smoothMPC_lub19, smoothMPC_sub19, smoothMPC_riub19, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb21, smoothMPC_z20, smoothMPC_lbIdx20, smoothMPC_llb20, smoothMPC_slb20, smoothMPC_rilb20, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z20, smoothMPC_ubIdx20, params->ub21, smoothMPC_lub20, smoothMPC_sub20, smoothMPC_riub20, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb22, smoothMPC_z21, smoothMPC_lbIdx21, smoothMPC_llb21, smoothMPC_slb21, smoothMPC_rilb21, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z21, smoothMPC_ubIdx21, params->ub22, smoothMPC_lub21, smoothMPC_sub21, smoothMPC_riub21, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb23, smoothMPC_z22, smoothMPC_lbIdx22, smoothMPC_llb22, smoothMPC_slb22, smoothMPC_rilb22, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z22, smoothMPC_ubIdx22, params->ub23, smoothMPC_lub22, smoothMPC_sub22, smoothMPC_riub22, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb24, smoothMPC_z23, smoothMPC_lbIdx23, smoothMPC_llb23, smoothMPC_slb23, smoothMPC_rilb23, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z23, smoothMPC_ubIdx23, params->ub24, smoothMPC_lub23, smoothMPC_sub23, smoothMPC_riub23, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb25, smoothMPC_z24, smoothMPC_lbIdx24, smoothMPC_llb24, smoothMPC_slb24, smoothMPC_rilb24, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z24, smoothMPC_ubIdx24, params->ub25, smoothMPC_lub24, smoothMPC_sub24, smoothMPC_riub24, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb26, smoothMPC_z25, smoothMPC_lbIdx25, smoothMPC_llb25, smoothMPC_slb25, smoothMPC_rilb25, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z25, smoothMPC_ubIdx25, params->ub26, smoothMPC_lub25, smoothMPC_sub25, smoothMPC_riub25, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb27, smoothMPC_z26, smoothMPC_lbIdx26, smoothMPC_llb26, smoothMPC_slb26, smoothMPC_rilb26, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z26, smoothMPC_ubIdx26, params->ub27, smoothMPC_lub26, smoothMPC_sub26, smoothMPC_riub26, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb28, smoothMPC_z27, smoothMPC_lbIdx27, smoothMPC_llb27, smoothMPC_slb27, smoothMPC_rilb27, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z27, smoothMPC_ubIdx27, params->ub28, smoothMPC_lub27, smoothMPC_sub27, smoothMPC_riub27, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb29, smoothMPC_z28, smoothMPC_lbIdx28, smoothMPC_llb28, smoothMPC_slb28, smoothMPC_rilb28, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z28, smoothMPC_ubIdx28, params->ub29, smoothMPC_lub28, smoothMPC_sub28, smoothMPC_riub28, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb30, smoothMPC_z29, smoothMPC_lbIdx29, smoothMPC_llb29, smoothMPC_slb29, smoothMPC_rilb29, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z29, smoothMPC_ubIdx29, params->ub30, smoothMPC_lub29, smoothMPC_sub29, smoothMPC_riub29, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb31, smoothMPC_z30, smoothMPC_lbIdx30, smoothMPC_llb30, smoothMPC_slb30, smoothMPC_rilb30, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z30, smoothMPC_ubIdx30, params->ub31, smoothMPC_lub30, smoothMPC_sub30, smoothMPC_riub30, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb32, smoothMPC_z31, smoothMPC_lbIdx31, smoothMPC_llb31, smoothMPC_slb31, smoothMPC_rilb31, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z31, smoothMPC_ubIdx31, params->ub32, smoothMPC_lub31, smoothMPC_sub31, smoothMPC_riub31, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb33, smoothMPC_z32, smoothMPC_lbIdx32, smoothMPC_llb32, smoothMPC_slb32, smoothMPC_rilb32, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z32, smoothMPC_ubIdx32, params->ub33, smoothMPC_lub32, smoothMPC_sub32, smoothMPC_riub32, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb34, smoothMPC_z33, smoothMPC_lbIdx33, smoothMPC_llb33, smoothMPC_slb33, smoothMPC_rilb33, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z33, smoothMPC_ubIdx33, params->ub34, smoothMPC_lub33, smoothMPC_sub33, smoothMPC_riub33, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb35, smoothMPC_z34, smoothMPC_lbIdx34, smoothMPC_llb34, smoothMPC_slb34, smoothMPC_rilb34, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z34, smoothMPC_ubIdx34, params->ub35, smoothMPC_lub34, smoothMPC_sub34, smoothMPC_riub34, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb36, smoothMPC_z35, smoothMPC_lbIdx35, smoothMPC_llb35, smoothMPC_slb35, smoothMPC_rilb35, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z35, smoothMPC_ubIdx35, params->ub36, smoothMPC_lub35, smoothMPC_sub35, smoothMPC_riub35, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb37, smoothMPC_z36, smoothMPC_lbIdx36, smoothMPC_llb36, smoothMPC_slb36, smoothMPC_rilb36, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z36, smoothMPC_ubIdx36, params->ub37, smoothMPC_lub36, smoothMPC_sub36, smoothMPC_riub36, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb38, smoothMPC_z37, smoothMPC_lbIdx37, smoothMPC_llb37, smoothMPC_slb37, smoothMPC_rilb37, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z37, smoothMPC_ubIdx37, params->ub38, smoothMPC_lub37, smoothMPC_sub37, smoothMPC_riub37, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb39, smoothMPC_z38, smoothMPC_lbIdx38, smoothMPC_llb38, smoothMPC_slb38, smoothMPC_rilb38, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z38, smoothMPC_ubIdx38, params->ub39, smoothMPC_lub38, smoothMPC_sub38, smoothMPC_riub38, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb40, smoothMPC_z39, smoothMPC_lbIdx39, smoothMPC_llb39, smoothMPC_slb39, smoothMPC_rilb39, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z39, smoothMPC_ubIdx39, params->ub40, smoothMPC_lub39, smoothMPC_sub39, smoothMPC_riub39, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb41, smoothMPC_z40, smoothMPC_lbIdx40, smoothMPC_llb40, smoothMPC_slb40, smoothMPC_rilb40, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z40, smoothMPC_ubIdx40, params->ub41, smoothMPC_lub40, smoothMPC_sub40, smoothMPC_riub40, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb42, smoothMPC_z41, smoothMPC_lbIdx41, smoothMPC_llb41, smoothMPC_slb41, smoothMPC_rilb41, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z41, smoothMPC_ubIdx41, params->ub42, smoothMPC_lub41, smoothMPC_sub41, smoothMPC_riub41, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb43, smoothMPC_z42, smoothMPC_lbIdx42, smoothMPC_llb42, smoothMPC_slb42, smoothMPC_rilb42, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z42, smoothMPC_ubIdx42, params->ub43, smoothMPC_lub42, smoothMPC_sub42, smoothMPC_riub42, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb44, smoothMPC_z43, smoothMPC_lbIdx43, smoothMPC_llb43, smoothMPC_slb43, smoothMPC_rilb43, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z43, smoothMPC_ubIdx43, params->ub44, smoothMPC_lub43, smoothMPC_sub43, smoothMPC_riub43, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb45, smoothMPC_z44, smoothMPC_lbIdx44, smoothMPC_llb44, smoothMPC_slb44, smoothMPC_rilb44, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z44, smoothMPC_ubIdx44, params->ub45, smoothMPC_lub44, smoothMPC_sub44, smoothMPC_riub44, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb46, smoothMPC_z45, smoothMPC_lbIdx45, smoothMPC_llb45, smoothMPC_slb45, smoothMPC_rilb45, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z45, smoothMPC_ubIdx45, params->ub46, smoothMPC_lub45, smoothMPC_sub45, smoothMPC_riub45, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb47, smoothMPC_z46, smoothMPC_lbIdx46, smoothMPC_llb46, smoothMPC_slb46, smoothMPC_rilb46, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z46, smoothMPC_ubIdx46, params->ub47, smoothMPC_lub46, smoothMPC_sub46, smoothMPC_riub46, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb48, smoothMPC_z47, smoothMPC_lbIdx47, smoothMPC_llb47, smoothMPC_slb47, smoothMPC_rilb47, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z47, smoothMPC_ubIdx47, params->ub48, smoothMPC_lub47, smoothMPC_sub47, smoothMPC_riub47, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb49, smoothMPC_z48, smoothMPC_lbIdx48, smoothMPC_llb48, smoothMPC_slb48, smoothMPC_rilb48, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z48, smoothMPC_ubIdx48, params->ub49, smoothMPC_lub48, smoothMPC_sub48, smoothMPC_riub48, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb50, smoothMPC_z49, smoothMPC_lbIdx49, smoothMPC_llb49, smoothMPC_slb49, smoothMPC_rilb49, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z49, smoothMPC_ubIdx49, params->ub50, smoothMPC_lub49, smoothMPC_sub49, smoothMPC_riub49, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb51, smoothMPC_z50, smoothMPC_lbIdx50, smoothMPC_llb50, smoothMPC_slb50, smoothMPC_rilb50, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z50, smoothMPC_ubIdx50, params->ub51, smoothMPC_lub50, smoothMPC_sub50, smoothMPC_riub50, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb52, smoothMPC_z51, smoothMPC_lbIdx51, smoothMPC_llb51, smoothMPC_slb51, smoothMPC_rilb51, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z51, smoothMPC_ubIdx51, params->ub52, smoothMPC_lub51, smoothMPC_sub51, smoothMPC_riub51, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb53, smoothMPC_z52, smoothMPC_lbIdx52, smoothMPC_llb52, smoothMPC_slb52, smoothMPC_rilb52, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z52, smoothMPC_ubIdx52, params->ub53, smoothMPC_lub52, smoothMPC_sub52, smoothMPC_riub52, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb54, smoothMPC_z53, smoothMPC_lbIdx53, smoothMPC_llb53, smoothMPC_slb53, smoothMPC_rilb53, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z53, smoothMPC_ubIdx53, params->ub54, smoothMPC_lub53, smoothMPC_sub53, smoothMPC_riub53, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb55, smoothMPC_z54, smoothMPC_lbIdx54, smoothMPC_llb54, smoothMPC_slb54, smoothMPC_rilb54, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z54, smoothMPC_ubIdx54, params->ub55, smoothMPC_lub54, smoothMPC_sub54, smoothMPC_riub54, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb56, smoothMPC_z55, smoothMPC_lbIdx55, smoothMPC_llb55, smoothMPC_slb55, smoothMPC_rilb55, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z55, smoothMPC_ubIdx55, params->ub56, smoothMPC_lub55, smoothMPC_sub55, smoothMPC_riub55, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb57, smoothMPC_z56, smoothMPC_lbIdx56, smoothMPC_llb56, smoothMPC_slb56, smoothMPC_rilb56, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z56, smoothMPC_ubIdx56, params->ub57, smoothMPC_lub56, smoothMPC_sub56, smoothMPC_riub56, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb58, smoothMPC_z57, smoothMPC_lbIdx57, smoothMPC_llb57, smoothMPC_slb57, smoothMPC_rilb57, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z57, smoothMPC_ubIdx57, params->ub58, smoothMPC_lub57, smoothMPC_sub57, smoothMPC_riub57, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb59, smoothMPC_z58, smoothMPC_lbIdx58, smoothMPC_llb58, smoothMPC_slb58, smoothMPC_rilb58, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z58, smoothMPC_ubIdx58, params->ub59, smoothMPC_lub58, smoothMPC_sub58, smoothMPC_riub58, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb60, smoothMPC_z59, smoothMPC_lbIdx59, smoothMPC_llb59, smoothMPC_slb59, smoothMPC_rilb59, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z59, smoothMPC_ubIdx59, params->ub60, smoothMPC_lub59, smoothMPC_sub59, smoothMPC_riub59, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb61, smoothMPC_z60, smoothMPC_lbIdx60, smoothMPC_llb60, smoothMPC_slb60, smoothMPC_rilb60, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z60, smoothMPC_ubIdx60, params->ub61, smoothMPC_lub60, smoothMPC_sub60, smoothMPC_riub60, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb62, smoothMPC_z61, smoothMPC_lbIdx61, smoothMPC_llb61, smoothMPC_slb61, smoothMPC_rilb61, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z61, smoothMPC_ubIdx61, params->ub62, smoothMPC_lub61, smoothMPC_sub61, smoothMPC_riub61, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb63, smoothMPC_z62, smoothMPC_lbIdx62, smoothMPC_llb62, smoothMPC_slb62, smoothMPC_rilb62, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z62, smoothMPC_ubIdx62, params->ub63, smoothMPC_lub62, smoothMPC_sub62, smoothMPC_riub62, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb64, smoothMPC_z63, smoothMPC_lbIdx63, smoothMPC_llb63, smoothMPC_slb63, smoothMPC_rilb63, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z63, smoothMPC_ubIdx63, params->ub64, smoothMPC_lub63, smoothMPC_sub63, smoothMPC_riub63, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb65, smoothMPC_z64, smoothMPC_lbIdx64, smoothMPC_llb64, smoothMPC_slb64, smoothMPC_rilb64, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z64, smoothMPC_ubIdx64, params->ub65, smoothMPC_lub64, smoothMPC_sub64, smoothMPC_riub64, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb66, smoothMPC_z65, smoothMPC_lbIdx65, smoothMPC_llb65, smoothMPC_slb65, smoothMPC_rilb65, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z65, smoothMPC_ubIdx65, params->ub66, smoothMPC_lub65, smoothMPC_sub65, smoothMPC_riub65, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb67, smoothMPC_z66, smoothMPC_lbIdx66, smoothMPC_llb66, smoothMPC_slb66, smoothMPC_rilb66, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z66, smoothMPC_ubIdx66, params->ub67, smoothMPC_lub66, smoothMPC_sub66, smoothMPC_riub66, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb68, smoothMPC_z67, smoothMPC_lbIdx67, smoothMPC_llb67, smoothMPC_slb67, smoothMPC_rilb67, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z67, smoothMPC_ubIdx67, params->ub68, smoothMPC_lub67, smoothMPC_sub67, smoothMPC_riub67, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb69, smoothMPC_z68, smoothMPC_lbIdx68, smoothMPC_llb68, smoothMPC_slb68, smoothMPC_rilb68, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z68, smoothMPC_ubIdx68, params->ub69, smoothMPC_lub68, smoothMPC_sub68, smoothMPC_riub68, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb70, smoothMPC_z69, smoothMPC_lbIdx69, smoothMPC_llb69, smoothMPC_slb69, smoothMPC_rilb69, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z69, smoothMPC_ubIdx69, params->ub70, smoothMPC_lub69, smoothMPC_sub69, smoothMPC_riub69, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb71, smoothMPC_z70, smoothMPC_lbIdx70, smoothMPC_llb70, smoothMPC_slb70, smoothMPC_rilb70, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z70, smoothMPC_ubIdx70, params->ub71, smoothMPC_lub70, smoothMPC_sub70, smoothMPC_riub70, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb72, smoothMPC_z71, smoothMPC_lbIdx71, smoothMPC_llb71, smoothMPC_slb71, smoothMPC_rilb71, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z71, smoothMPC_ubIdx71, params->ub72, smoothMPC_lub71, smoothMPC_sub71, smoothMPC_riub71, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb73, smoothMPC_z72, smoothMPC_lbIdx72, smoothMPC_llb72, smoothMPC_slb72, smoothMPC_rilb72, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z72, smoothMPC_ubIdx72, params->ub73, smoothMPC_lub72, smoothMPC_sub72, smoothMPC_riub72, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb74, smoothMPC_z73, smoothMPC_lbIdx73, smoothMPC_llb73, smoothMPC_slb73, smoothMPC_rilb73, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z73, smoothMPC_ubIdx73, params->ub74, smoothMPC_lub73, smoothMPC_sub73, smoothMPC_riub73, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb75, smoothMPC_z74, smoothMPC_lbIdx74, smoothMPC_llb74, smoothMPC_slb74, smoothMPC_rilb74, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z74, smoothMPC_ubIdx74, params->ub75, smoothMPC_lub74, smoothMPC_sub74, smoothMPC_riub74, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb76, smoothMPC_z75, smoothMPC_lbIdx75, smoothMPC_llb75, smoothMPC_slb75, smoothMPC_rilb75, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z75, smoothMPC_ubIdx75, params->ub76, smoothMPC_lub75, smoothMPC_sub75, smoothMPC_riub75, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb77, smoothMPC_z76, smoothMPC_lbIdx76, smoothMPC_llb76, smoothMPC_slb76, smoothMPC_rilb76, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z76, smoothMPC_ubIdx76, params->ub77, smoothMPC_lub76, smoothMPC_sub76, smoothMPC_riub76, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb78, smoothMPC_z77, smoothMPC_lbIdx77, smoothMPC_llb77, smoothMPC_slb77, smoothMPC_rilb77, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z77, smoothMPC_ubIdx77, params->ub78, smoothMPC_lub77, smoothMPC_sub77, smoothMPC_riub77, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb79, smoothMPC_z78, smoothMPC_lbIdx78, smoothMPC_llb78, smoothMPC_slb78, smoothMPC_rilb78, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z78, smoothMPC_ubIdx78, params->ub79, smoothMPC_lub78, smoothMPC_sub78, smoothMPC_riub78, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb80, smoothMPC_z79, smoothMPC_lbIdx79, smoothMPC_llb79, smoothMPC_slb79, smoothMPC_rilb79, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z79, smoothMPC_ubIdx79, params->ub80, smoothMPC_lub79, smoothMPC_sub79, smoothMPC_riub79, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb81, smoothMPC_z80, smoothMPC_lbIdx80, smoothMPC_llb80, smoothMPC_slb80, smoothMPC_rilb80, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z80, smoothMPC_ubIdx80, params->ub81, smoothMPC_lub80, smoothMPC_sub80, smoothMPC_riub80, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb82, smoothMPC_z81, smoothMPC_lbIdx81, smoothMPC_llb81, smoothMPC_slb81, smoothMPC_rilb81, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z81, smoothMPC_ubIdx81, params->ub82, smoothMPC_lub81, smoothMPC_sub81, smoothMPC_riub81, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb83, smoothMPC_z82, smoothMPC_lbIdx82, smoothMPC_llb82, smoothMPC_slb82, smoothMPC_rilb82, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z82, smoothMPC_ubIdx82, params->ub83, smoothMPC_lub82, smoothMPC_sub82, smoothMPC_riub82, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb84, smoothMPC_z83, smoothMPC_lbIdx83, smoothMPC_llb83, smoothMPC_slb83, smoothMPC_rilb83, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z83, smoothMPC_ubIdx83, params->ub84, smoothMPC_lub83, smoothMPC_sub83, smoothMPC_riub83, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb85, smoothMPC_z84, smoothMPC_lbIdx84, smoothMPC_llb84, smoothMPC_slb84, smoothMPC_rilb84, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z84, smoothMPC_ubIdx84, params->ub85, smoothMPC_lub84, smoothMPC_sub84, smoothMPC_riub84, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb86, smoothMPC_z85, smoothMPC_lbIdx85, smoothMPC_llb85, smoothMPC_slb85, smoothMPC_rilb85, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z85, smoothMPC_ubIdx85, params->ub86, smoothMPC_lub85, smoothMPC_sub85, smoothMPC_riub85, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb87, smoothMPC_z86, smoothMPC_lbIdx86, smoothMPC_llb86, smoothMPC_slb86, smoothMPC_rilb86, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z86, smoothMPC_ubIdx86, params->ub87, smoothMPC_lub86, smoothMPC_sub86, smoothMPC_riub86, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb88, smoothMPC_z87, smoothMPC_lbIdx87, smoothMPC_llb87, smoothMPC_slb87, smoothMPC_rilb87, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z87, smoothMPC_ubIdx87, params->ub88, smoothMPC_lub87, smoothMPC_sub87, smoothMPC_riub87, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb89, smoothMPC_z88, smoothMPC_lbIdx88, smoothMPC_llb88, smoothMPC_slb88, smoothMPC_rilb88, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z88, smoothMPC_ubIdx88, params->ub89, smoothMPC_lub88, smoothMPC_sub88, smoothMPC_riub88, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb90, smoothMPC_z89, smoothMPC_lbIdx89, smoothMPC_llb89, smoothMPC_slb89, smoothMPC_rilb89, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z89, smoothMPC_ubIdx89, params->ub90, smoothMPC_lub89, smoothMPC_sub89, smoothMPC_riub89, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb91, smoothMPC_z90, smoothMPC_lbIdx90, smoothMPC_llb90, smoothMPC_slb90, smoothMPC_rilb90, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z90, smoothMPC_ubIdx90, params->ub91, smoothMPC_lub90, smoothMPC_sub90, smoothMPC_riub90, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb92, smoothMPC_z91, smoothMPC_lbIdx91, smoothMPC_llb91, smoothMPC_slb91, smoothMPC_rilb91, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z91, smoothMPC_ubIdx91, params->ub92, smoothMPC_lub91, smoothMPC_sub91, smoothMPC_riub91, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb93, smoothMPC_z92, smoothMPC_lbIdx92, smoothMPC_llb92, smoothMPC_slb92, smoothMPC_rilb92, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z92, smoothMPC_ubIdx92, params->ub93, smoothMPC_lub92, smoothMPC_sub92, smoothMPC_riub92, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb94, smoothMPC_z93, smoothMPC_lbIdx93, smoothMPC_llb93, smoothMPC_slb93, smoothMPC_rilb93, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z93, smoothMPC_ubIdx93, params->ub94, smoothMPC_lub93, smoothMPC_sub93, smoothMPC_riub93, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb95, smoothMPC_z94, smoothMPC_lbIdx94, smoothMPC_llb94, smoothMPC_slb94, smoothMPC_rilb94, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z94, smoothMPC_ubIdx94, params->ub95, smoothMPC_lub94, smoothMPC_sub94, smoothMPC_riub94, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb96, smoothMPC_z95, smoothMPC_lbIdx95, smoothMPC_llb95, smoothMPC_slb95, smoothMPC_rilb95, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z95, smoothMPC_ubIdx95, params->ub96, smoothMPC_lub95, smoothMPC_sub95, smoothMPC_riub95, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb97, smoothMPC_z96, smoothMPC_lbIdx96, smoothMPC_llb96, smoothMPC_slb96, smoothMPC_rilb96, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z96, smoothMPC_ubIdx96, params->ub97, smoothMPC_lub96, smoothMPC_sub96, smoothMPC_riub96, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb98, smoothMPC_z97, smoothMPC_lbIdx97, smoothMPC_llb97, smoothMPC_slb97, smoothMPC_rilb97, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z97, smoothMPC_ubIdx97, params->ub98, smoothMPC_lub97, smoothMPC_sub97, smoothMPC_riub97, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb99, smoothMPC_z98, smoothMPC_lbIdx98, smoothMPC_llb98, smoothMPC_slb98, smoothMPC_rilb98, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z98, smoothMPC_ubIdx98, params->ub99, smoothMPC_lub98, smoothMPC_sub98, smoothMPC_riub98, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_3(params->lb100, smoothMPC_z99, smoothMPC_lbIdx99, smoothMPC_llb99, smoothMPC_slb99, smoothMPC_rilb99, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_3(smoothMPC_z99, smoothMPC_ubIdx99, params->ub100, smoothMPC_lub99, smoothMPC_sub99, smoothMPC_riub99, &info->dgap, &info->res_ineq);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub00, smoothMPC_sub00, smoothMPC_riub00, smoothMPC_llb00, smoothMPC_slb00, smoothMPC_rilb00, smoothMPC_lbIdx00, smoothMPC_ubIdx00, smoothMPC_grad_ineq00, smoothMPC_lubbysub00, smoothMPC_llbbyslb00);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub01, smoothMPC_sub01, smoothMPC_riub01, smoothMPC_llb01, smoothMPC_slb01, smoothMPC_rilb01, smoothMPC_lbIdx01, smoothMPC_ubIdx01, smoothMPC_grad_ineq01, smoothMPC_lubbysub01, smoothMPC_llbbyslb01);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub02, smoothMPC_sub02, smoothMPC_riub02, smoothMPC_llb02, smoothMPC_slb02, smoothMPC_rilb02, smoothMPC_lbIdx02, smoothMPC_ubIdx02, smoothMPC_grad_ineq02, smoothMPC_lubbysub02, smoothMPC_llbbyslb02);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub03, smoothMPC_sub03, smoothMPC_riub03, smoothMPC_llb03, smoothMPC_slb03, smoothMPC_rilb03, smoothMPC_lbIdx03, smoothMPC_ubIdx03, smoothMPC_grad_ineq03, smoothMPC_lubbysub03, smoothMPC_llbbyslb03);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub04, smoothMPC_sub04, smoothMPC_riub04, smoothMPC_llb04, smoothMPC_slb04, smoothMPC_rilb04, smoothMPC_lbIdx04, smoothMPC_ubIdx04, smoothMPC_grad_ineq04, smoothMPC_lubbysub04, smoothMPC_llbbyslb04);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub05, smoothMPC_sub05, smoothMPC_riub05, smoothMPC_llb05, smoothMPC_slb05, smoothMPC_rilb05, smoothMPC_lbIdx05, smoothMPC_ubIdx05, smoothMPC_grad_ineq05, smoothMPC_lubbysub05, smoothMPC_llbbyslb05);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub06, smoothMPC_sub06, smoothMPC_riub06, smoothMPC_llb06, smoothMPC_slb06, smoothMPC_rilb06, smoothMPC_lbIdx06, smoothMPC_ubIdx06, smoothMPC_grad_ineq06, smoothMPC_lubbysub06, smoothMPC_llbbyslb06);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub07, smoothMPC_sub07, smoothMPC_riub07, smoothMPC_llb07, smoothMPC_slb07, smoothMPC_rilb07, smoothMPC_lbIdx07, smoothMPC_ubIdx07, smoothMPC_grad_ineq07, smoothMPC_lubbysub07, smoothMPC_llbbyslb07);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub08, smoothMPC_sub08, smoothMPC_riub08, smoothMPC_llb08, smoothMPC_slb08, smoothMPC_rilb08, smoothMPC_lbIdx08, smoothMPC_ubIdx08, smoothMPC_grad_ineq08, smoothMPC_lubbysub08, smoothMPC_llbbyslb08);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub09, smoothMPC_sub09, smoothMPC_riub09, smoothMPC_llb09, smoothMPC_slb09, smoothMPC_rilb09, smoothMPC_lbIdx09, smoothMPC_ubIdx09, smoothMPC_grad_ineq09, smoothMPC_lubbysub09, smoothMPC_llbbyslb09);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub10, smoothMPC_sub10, smoothMPC_riub10, smoothMPC_llb10, smoothMPC_slb10, smoothMPC_rilb10, smoothMPC_lbIdx10, smoothMPC_ubIdx10, smoothMPC_grad_ineq10, smoothMPC_lubbysub10, smoothMPC_llbbyslb10);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub11, smoothMPC_sub11, smoothMPC_riub11, smoothMPC_llb11, smoothMPC_slb11, smoothMPC_rilb11, smoothMPC_lbIdx11, smoothMPC_ubIdx11, smoothMPC_grad_ineq11, smoothMPC_lubbysub11, smoothMPC_llbbyslb11);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub12, smoothMPC_sub12, smoothMPC_riub12, smoothMPC_llb12, smoothMPC_slb12, smoothMPC_rilb12, smoothMPC_lbIdx12, smoothMPC_ubIdx12, smoothMPC_grad_ineq12, smoothMPC_lubbysub12, smoothMPC_llbbyslb12);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub13, smoothMPC_sub13, smoothMPC_riub13, smoothMPC_llb13, smoothMPC_slb13, smoothMPC_rilb13, smoothMPC_lbIdx13, smoothMPC_ubIdx13, smoothMPC_grad_ineq13, smoothMPC_lubbysub13, smoothMPC_llbbyslb13);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub14, smoothMPC_sub14, smoothMPC_riub14, smoothMPC_llb14, smoothMPC_slb14, smoothMPC_rilb14, smoothMPC_lbIdx14, smoothMPC_ubIdx14, smoothMPC_grad_ineq14, smoothMPC_lubbysub14, smoothMPC_llbbyslb14);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub15, smoothMPC_sub15, smoothMPC_riub15, smoothMPC_llb15, smoothMPC_slb15, smoothMPC_rilb15, smoothMPC_lbIdx15, smoothMPC_ubIdx15, smoothMPC_grad_ineq15, smoothMPC_lubbysub15, smoothMPC_llbbyslb15);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub16, smoothMPC_sub16, smoothMPC_riub16, smoothMPC_llb16, smoothMPC_slb16, smoothMPC_rilb16, smoothMPC_lbIdx16, smoothMPC_ubIdx16, smoothMPC_grad_ineq16, smoothMPC_lubbysub16, smoothMPC_llbbyslb16);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub17, smoothMPC_sub17, smoothMPC_riub17, smoothMPC_llb17, smoothMPC_slb17, smoothMPC_rilb17, smoothMPC_lbIdx17, smoothMPC_ubIdx17, smoothMPC_grad_ineq17, smoothMPC_lubbysub17, smoothMPC_llbbyslb17);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub18, smoothMPC_sub18, smoothMPC_riub18, smoothMPC_llb18, smoothMPC_slb18, smoothMPC_rilb18, smoothMPC_lbIdx18, smoothMPC_ubIdx18, smoothMPC_grad_ineq18, smoothMPC_lubbysub18, smoothMPC_llbbyslb18);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub19, smoothMPC_sub19, smoothMPC_riub19, smoothMPC_llb19, smoothMPC_slb19, smoothMPC_rilb19, smoothMPC_lbIdx19, smoothMPC_ubIdx19, smoothMPC_grad_ineq19, smoothMPC_lubbysub19, smoothMPC_llbbyslb19);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub20, smoothMPC_sub20, smoothMPC_riub20, smoothMPC_llb20, smoothMPC_slb20, smoothMPC_rilb20, smoothMPC_lbIdx20, smoothMPC_ubIdx20, smoothMPC_grad_ineq20, smoothMPC_lubbysub20, smoothMPC_llbbyslb20);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub21, smoothMPC_sub21, smoothMPC_riub21, smoothMPC_llb21, smoothMPC_slb21, smoothMPC_rilb21, smoothMPC_lbIdx21, smoothMPC_ubIdx21, smoothMPC_grad_ineq21, smoothMPC_lubbysub21, smoothMPC_llbbyslb21);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub22, smoothMPC_sub22, smoothMPC_riub22, smoothMPC_llb22, smoothMPC_slb22, smoothMPC_rilb22, smoothMPC_lbIdx22, smoothMPC_ubIdx22, smoothMPC_grad_ineq22, smoothMPC_lubbysub22, smoothMPC_llbbyslb22);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub23, smoothMPC_sub23, smoothMPC_riub23, smoothMPC_llb23, smoothMPC_slb23, smoothMPC_rilb23, smoothMPC_lbIdx23, smoothMPC_ubIdx23, smoothMPC_grad_ineq23, smoothMPC_lubbysub23, smoothMPC_llbbyslb23);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub24, smoothMPC_sub24, smoothMPC_riub24, smoothMPC_llb24, smoothMPC_slb24, smoothMPC_rilb24, smoothMPC_lbIdx24, smoothMPC_ubIdx24, smoothMPC_grad_ineq24, smoothMPC_lubbysub24, smoothMPC_llbbyslb24);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub25, smoothMPC_sub25, smoothMPC_riub25, smoothMPC_llb25, smoothMPC_slb25, smoothMPC_rilb25, smoothMPC_lbIdx25, smoothMPC_ubIdx25, smoothMPC_grad_ineq25, smoothMPC_lubbysub25, smoothMPC_llbbyslb25);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub26, smoothMPC_sub26, smoothMPC_riub26, smoothMPC_llb26, smoothMPC_slb26, smoothMPC_rilb26, smoothMPC_lbIdx26, smoothMPC_ubIdx26, smoothMPC_grad_ineq26, smoothMPC_lubbysub26, smoothMPC_llbbyslb26);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub27, smoothMPC_sub27, smoothMPC_riub27, smoothMPC_llb27, smoothMPC_slb27, smoothMPC_rilb27, smoothMPC_lbIdx27, smoothMPC_ubIdx27, smoothMPC_grad_ineq27, smoothMPC_lubbysub27, smoothMPC_llbbyslb27);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub28, smoothMPC_sub28, smoothMPC_riub28, smoothMPC_llb28, smoothMPC_slb28, smoothMPC_rilb28, smoothMPC_lbIdx28, smoothMPC_ubIdx28, smoothMPC_grad_ineq28, smoothMPC_lubbysub28, smoothMPC_llbbyslb28);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub29, smoothMPC_sub29, smoothMPC_riub29, smoothMPC_llb29, smoothMPC_slb29, smoothMPC_rilb29, smoothMPC_lbIdx29, smoothMPC_ubIdx29, smoothMPC_grad_ineq29, smoothMPC_lubbysub29, smoothMPC_llbbyslb29);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub30, smoothMPC_sub30, smoothMPC_riub30, smoothMPC_llb30, smoothMPC_slb30, smoothMPC_rilb30, smoothMPC_lbIdx30, smoothMPC_ubIdx30, smoothMPC_grad_ineq30, smoothMPC_lubbysub30, smoothMPC_llbbyslb30);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub31, smoothMPC_sub31, smoothMPC_riub31, smoothMPC_llb31, smoothMPC_slb31, smoothMPC_rilb31, smoothMPC_lbIdx31, smoothMPC_ubIdx31, smoothMPC_grad_ineq31, smoothMPC_lubbysub31, smoothMPC_llbbyslb31);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub32, smoothMPC_sub32, smoothMPC_riub32, smoothMPC_llb32, smoothMPC_slb32, smoothMPC_rilb32, smoothMPC_lbIdx32, smoothMPC_ubIdx32, smoothMPC_grad_ineq32, smoothMPC_lubbysub32, smoothMPC_llbbyslb32);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub33, smoothMPC_sub33, smoothMPC_riub33, smoothMPC_llb33, smoothMPC_slb33, smoothMPC_rilb33, smoothMPC_lbIdx33, smoothMPC_ubIdx33, smoothMPC_grad_ineq33, smoothMPC_lubbysub33, smoothMPC_llbbyslb33);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub34, smoothMPC_sub34, smoothMPC_riub34, smoothMPC_llb34, smoothMPC_slb34, smoothMPC_rilb34, smoothMPC_lbIdx34, smoothMPC_ubIdx34, smoothMPC_grad_ineq34, smoothMPC_lubbysub34, smoothMPC_llbbyslb34);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub35, smoothMPC_sub35, smoothMPC_riub35, smoothMPC_llb35, smoothMPC_slb35, smoothMPC_rilb35, smoothMPC_lbIdx35, smoothMPC_ubIdx35, smoothMPC_grad_ineq35, smoothMPC_lubbysub35, smoothMPC_llbbyslb35);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub36, smoothMPC_sub36, smoothMPC_riub36, smoothMPC_llb36, smoothMPC_slb36, smoothMPC_rilb36, smoothMPC_lbIdx36, smoothMPC_ubIdx36, smoothMPC_grad_ineq36, smoothMPC_lubbysub36, smoothMPC_llbbyslb36);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub37, smoothMPC_sub37, smoothMPC_riub37, smoothMPC_llb37, smoothMPC_slb37, smoothMPC_rilb37, smoothMPC_lbIdx37, smoothMPC_ubIdx37, smoothMPC_grad_ineq37, smoothMPC_lubbysub37, smoothMPC_llbbyslb37);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub38, smoothMPC_sub38, smoothMPC_riub38, smoothMPC_llb38, smoothMPC_slb38, smoothMPC_rilb38, smoothMPC_lbIdx38, smoothMPC_ubIdx38, smoothMPC_grad_ineq38, smoothMPC_lubbysub38, smoothMPC_llbbyslb38);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub39, smoothMPC_sub39, smoothMPC_riub39, smoothMPC_llb39, smoothMPC_slb39, smoothMPC_rilb39, smoothMPC_lbIdx39, smoothMPC_ubIdx39, smoothMPC_grad_ineq39, smoothMPC_lubbysub39, smoothMPC_llbbyslb39);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub40, smoothMPC_sub40, smoothMPC_riub40, smoothMPC_llb40, smoothMPC_slb40, smoothMPC_rilb40, smoothMPC_lbIdx40, smoothMPC_ubIdx40, smoothMPC_grad_ineq40, smoothMPC_lubbysub40, smoothMPC_llbbyslb40);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub41, smoothMPC_sub41, smoothMPC_riub41, smoothMPC_llb41, smoothMPC_slb41, smoothMPC_rilb41, smoothMPC_lbIdx41, smoothMPC_ubIdx41, smoothMPC_grad_ineq41, smoothMPC_lubbysub41, smoothMPC_llbbyslb41);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub42, smoothMPC_sub42, smoothMPC_riub42, smoothMPC_llb42, smoothMPC_slb42, smoothMPC_rilb42, smoothMPC_lbIdx42, smoothMPC_ubIdx42, smoothMPC_grad_ineq42, smoothMPC_lubbysub42, smoothMPC_llbbyslb42);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub43, smoothMPC_sub43, smoothMPC_riub43, smoothMPC_llb43, smoothMPC_slb43, smoothMPC_rilb43, smoothMPC_lbIdx43, smoothMPC_ubIdx43, smoothMPC_grad_ineq43, smoothMPC_lubbysub43, smoothMPC_llbbyslb43);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub44, smoothMPC_sub44, smoothMPC_riub44, smoothMPC_llb44, smoothMPC_slb44, smoothMPC_rilb44, smoothMPC_lbIdx44, smoothMPC_ubIdx44, smoothMPC_grad_ineq44, smoothMPC_lubbysub44, smoothMPC_llbbyslb44);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub45, smoothMPC_sub45, smoothMPC_riub45, smoothMPC_llb45, smoothMPC_slb45, smoothMPC_rilb45, smoothMPC_lbIdx45, smoothMPC_ubIdx45, smoothMPC_grad_ineq45, smoothMPC_lubbysub45, smoothMPC_llbbyslb45);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub46, smoothMPC_sub46, smoothMPC_riub46, smoothMPC_llb46, smoothMPC_slb46, smoothMPC_rilb46, smoothMPC_lbIdx46, smoothMPC_ubIdx46, smoothMPC_grad_ineq46, smoothMPC_lubbysub46, smoothMPC_llbbyslb46);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub47, smoothMPC_sub47, smoothMPC_riub47, smoothMPC_llb47, smoothMPC_slb47, smoothMPC_rilb47, smoothMPC_lbIdx47, smoothMPC_ubIdx47, smoothMPC_grad_ineq47, smoothMPC_lubbysub47, smoothMPC_llbbyslb47);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub48, smoothMPC_sub48, smoothMPC_riub48, smoothMPC_llb48, smoothMPC_slb48, smoothMPC_rilb48, smoothMPC_lbIdx48, smoothMPC_ubIdx48, smoothMPC_grad_ineq48, smoothMPC_lubbysub48, smoothMPC_llbbyslb48);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub49, smoothMPC_sub49, smoothMPC_riub49, smoothMPC_llb49, smoothMPC_slb49, smoothMPC_rilb49, smoothMPC_lbIdx49, smoothMPC_ubIdx49, smoothMPC_grad_ineq49, smoothMPC_lubbysub49, smoothMPC_llbbyslb49);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub50, smoothMPC_sub50, smoothMPC_riub50, smoothMPC_llb50, smoothMPC_slb50, smoothMPC_rilb50, smoothMPC_lbIdx50, smoothMPC_ubIdx50, smoothMPC_grad_ineq50, smoothMPC_lubbysub50, smoothMPC_llbbyslb50);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub51, smoothMPC_sub51, smoothMPC_riub51, smoothMPC_llb51, smoothMPC_slb51, smoothMPC_rilb51, smoothMPC_lbIdx51, smoothMPC_ubIdx51, smoothMPC_grad_ineq51, smoothMPC_lubbysub51, smoothMPC_llbbyslb51);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub52, smoothMPC_sub52, smoothMPC_riub52, smoothMPC_llb52, smoothMPC_slb52, smoothMPC_rilb52, smoothMPC_lbIdx52, smoothMPC_ubIdx52, smoothMPC_grad_ineq52, smoothMPC_lubbysub52, smoothMPC_llbbyslb52);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub53, smoothMPC_sub53, smoothMPC_riub53, smoothMPC_llb53, smoothMPC_slb53, smoothMPC_rilb53, smoothMPC_lbIdx53, smoothMPC_ubIdx53, smoothMPC_grad_ineq53, smoothMPC_lubbysub53, smoothMPC_llbbyslb53);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub54, smoothMPC_sub54, smoothMPC_riub54, smoothMPC_llb54, smoothMPC_slb54, smoothMPC_rilb54, smoothMPC_lbIdx54, smoothMPC_ubIdx54, smoothMPC_grad_ineq54, smoothMPC_lubbysub54, smoothMPC_llbbyslb54);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub55, smoothMPC_sub55, smoothMPC_riub55, smoothMPC_llb55, smoothMPC_slb55, smoothMPC_rilb55, smoothMPC_lbIdx55, smoothMPC_ubIdx55, smoothMPC_grad_ineq55, smoothMPC_lubbysub55, smoothMPC_llbbyslb55);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub56, smoothMPC_sub56, smoothMPC_riub56, smoothMPC_llb56, smoothMPC_slb56, smoothMPC_rilb56, smoothMPC_lbIdx56, smoothMPC_ubIdx56, smoothMPC_grad_ineq56, smoothMPC_lubbysub56, smoothMPC_llbbyslb56);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub57, smoothMPC_sub57, smoothMPC_riub57, smoothMPC_llb57, smoothMPC_slb57, smoothMPC_rilb57, smoothMPC_lbIdx57, smoothMPC_ubIdx57, smoothMPC_grad_ineq57, smoothMPC_lubbysub57, smoothMPC_llbbyslb57);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub58, smoothMPC_sub58, smoothMPC_riub58, smoothMPC_llb58, smoothMPC_slb58, smoothMPC_rilb58, smoothMPC_lbIdx58, smoothMPC_ubIdx58, smoothMPC_grad_ineq58, smoothMPC_lubbysub58, smoothMPC_llbbyslb58);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub59, smoothMPC_sub59, smoothMPC_riub59, smoothMPC_llb59, smoothMPC_slb59, smoothMPC_rilb59, smoothMPC_lbIdx59, smoothMPC_ubIdx59, smoothMPC_grad_ineq59, smoothMPC_lubbysub59, smoothMPC_llbbyslb59);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub60, smoothMPC_sub60, smoothMPC_riub60, smoothMPC_llb60, smoothMPC_slb60, smoothMPC_rilb60, smoothMPC_lbIdx60, smoothMPC_ubIdx60, smoothMPC_grad_ineq60, smoothMPC_lubbysub60, smoothMPC_llbbyslb60);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub61, smoothMPC_sub61, smoothMPC_riub61, smoothMPC_llb61, smoothMPC_slb61, smoothMPC_rilb61, smoothMPC_lbIdx61, smoothMPC_ubIdx61, smoothMPC_grad_ineq61, smoothMPC_lubbysub61, smoothMPC_llbbyslb61);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub62, smoothMPC_sub62, smoothMPC_riub62, smoothMPC_llb62, smoothMPC_slb62, smoothMPC_rilb62, smoothMPC_lbIdx62, smoothMPC_ubIdx62, smoothMPC_grad_ineq62, smoothMPC_lubbysub62, smoothMPC_llbbyslb62);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub63, smoothMPC_sub63, smoothMPC_riub63, smoothMPC_llb63, smoothMPC_slb63, smoothMPC_rilb63, smoothMPC_lbIdx63, smoothMPC_ubIdx63, smoothMPC_grad_ineq63, smoothMPC_lubbysub63, smoothMPC_llbbyslb63);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub64, smoothMPC_sub64, smoothMPC_riub64, smoothMPC_llb64, smoothMPC_slb64, smoothMPC_rilb64, smoothMPC_lbIdx64, smoothMPC_ubIdx64, smoothMPC_grad_ineq64, smoothMPC_lubbysub64, smoothMPC_llbbyslb64);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub65, smoothMPC_sub65, smoothMPC_riub65, smoothMPC_llb65, smoothMPC_slb65, smoothMPC_rilb65, smoothMPC_lbIdx65, smoothMPC_ubIdx65, smoothMPC_grad_ineq65, smoothMPC_lubbysub65, smoothMPC_llbbyslb65);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub66, smoothMPC_sub66, smoothMPC_riub66, smoothMPC_llb66, smoothMPC_slb66, smoothMPC_rilb66, smoothMPC_lbIdx66, smoothMPC_ubIdx66, smoothMPC_grad_ineq66, smoothMPC_lubbysub66, smoothMPC_llbbyslb66);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub67, smoothMPC_sub67, smoothMPC_riub67, smoothMPC_llb67, smoothMPC_slb67, smoothMPC_rilb67, smoothMPC_lbIdx67, smoothMPC_ubIdx67, smoothMPC_grad_ineq67, smoothMPC_lubbysub67, smoothMPC_llbbyslb67);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub68, smoothMPC_sub68, smoothMPC_riub68, smoothMPC_llb68, smoothMPC_slb68, smoothMPC_rilb68, smoothMPC_lbIdx68, smoothMPC_ubIdx68, smoothMPC_grad_ineq68, smoothMPC_lubbysub68, smoothMPC_llbbyslb68);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub69, smoothMPC_sub69, smoothMPC_riub69, smoothMPC_llb69, smoothMPC_slb69, smoothMPC_rilb69, smoothMPC_lbIdx69, smoothMPC_ubIdx69, smoothMPC_grad_ineq69, smoothMPC_lubbysub69, smoothMPC_llbbyslb69);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub70, smoothMPC_sub70, smoothMPC_riub70, smoothMPC_llb70, smoothMPC_slb70, smoothMPC_rilb70, smoothMPC_lbIdx70, smoothMPC_ubIdx70, smoothMPC_grad_ineq70, smoothMPC_lubbysub70, smoothMPC_llbbyslb70);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub71, smoothMPC_sub71, smoothMPC_riub71, smoothMPC_llb71, smoothMPC_slb71, smoothMPC_rilb71, smoothMPC_lbIdx71, smoothMPC_ubIdx71, smoothMPC_grad_ineq71, smoothMPC_lubbysub71, smoothMPC_llbbyslb71);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub72, smoothMPC_sub72, smoothMPC_riub72, smoothMPC_llb72, smoothMPC_slb72, smoothMPC_rilb72, smoothMPC_lbIdx72, smoothMPC_ubIdx72, smoothMPC_grad_ineq72, smoothMPC_lubbysub72, smoothMPC_llbbyslb72);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub73, smoothMPC_sub73, smoothMPC_riub73, smoothMPC_llb73, smoothMPC_slb73, smoothMPC_rilb73, smoothMPC_lbIdx73, smoothMPC_ubIdx73, smoothMPC_grad_ineq73, smoothMPC_lubbysub73, smoothMPC_llbbyslb73);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub74, smoothMPC_sub74, smoothMPC_riub74, smoothMPC_llb74, smoothMPC_slb74, smoothMPC_rilb74, smoothMPC_lbIdx74, smoothMPC_ubIdx74, smoothMPC_grad_ineq74, smoothMPC_lubbysub74, smoothMPC_llbbyslb74);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub75, smoothMPC_sub75, smoothMPC_riub75, smoothMPC_llb75, smoothMPC_slb75, smoothMPC_rilb75, smoothMPC_lbIdx75, smoothMPC_ubIdx75, smoothMPC_grad_ineq75, smoothMPC_lubbysub75, smoothMPC_llbbyslb75);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub76, smoothMPC_sub76, smoothMPC_riub76, smoothMPC_llb76, smoothMPC_slb76, smoothMPC_rilb76, smoothMPC_lbIdx76, smoothMPC_ubIdx76, smoothMPC_grad_ineq76, smoothMPC_lubbysub76, smoothMPC_llbbyslb76);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub77, smoothMPC_sub77, smoothMPC_riub77, smoothMPC_llb77, smoothMPC_slb77, smoothMPC_rilb77, smoothMPC_lbIdx77, smoothMPC_ubIdx77, smoothMPC_grad_ineq77, smoothMPC_lubbysub77, smoothMPC_llbbyslb77);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub78, smoothMPC_sub78, smoothMPC_riub78, smoothMPC_llb78, smoothMPC_slb78, smoothMPC_rilb78, smoothMPC_lbIdx78, smoothMPC_ubIdx78, smoothMPC_grad_ineq78, smoothMPC_lubbysub78, smoothMPC_llbbyslb78);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub79, smoothMPC_sub79, smoothMPC_riub79, smoothMPC_llb79, smoothMPC_slb79, smoothMPC_rilb79, smoothMPC_lbIdx79, smoothMPC_ubIdx79, smoothMPC_grad_ineq79, smoothMPC_lubbysub79, smoothMPC_llbbyslb79);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub80, smoothMPC_sub80, smoothMPC_riub80, smoothMPC_llb80, smoothMPC_slb80, smoothMPC_rilb80, smoothMPC_lbIdx80, smoothMPC_ubIdx80, smoothMPC_grad_ineq80, smoothMPC_lubbysub80, smoothMPC_llbbyslb80);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub81, smoothMPC_sub81, smoothMPC_riub81, smoothMPC_llb81, smoothMPC_slb81, smoothMPC_rilb81, smoothMPC_lbIdx81, smoothMPC_ubIdx81, smoothMPC_grad_ineq81, smoothMPC_lubbysub81, smoothMPC_llbbyslb81);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub82, smoothMPC_sub82, smoothMPC_riub82, smoothMPC_llb82, smoothMPC_slb82, smoothMPC_rilb82, smoothMPC_lbIdx82, smoothMPC_ubIdx82, smoothMPC_grad_ineq82, smoothMPC_lubbysub82, smoothMPC_llbbyslb82);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub83, smoothMPC_sub83, smoothMPC_riub83, smoothMPC_llb83, smoothMPC_slb83, smoothMPC_rilb83, smoothMPC_lbIdx83, smoothMPC_ubIdx83, smoothMPC_grad_ineq83, smoothMPC_lubbysub83, smoothMPC_llbbyslb83);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub84, smoothMPC_sub84, smoothMPC_riub84, smoothMPC_llb84, smoothMPC_slb84, smoothMPC_rilb84, smoothMPC_lbIdx84, smoothMPC_ubIdx84, smoothMPC_grad_ineq84, smoothMPC_lubbysub84, smoothMPC_llbbyslb84);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub85, smoothMPC_sub85, smoothMPC_riub85, smoothMPC_llb85, smoothMPC_slb85, smoothMPC_rilb85, smoothMPC_lbIdx85, smoothMPC_ubIdx85, smoothMPC_grad_ineq85, smoothMPC_lubbysub85, smoothMPC_llbbyslb85);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub86, smoothMPC_sub86, smoothMPC_riub86, smoothMPC_llb86, smoothMPC_slb86, smoothMPC_rilb86, smoothMPC_lbIdx86, smoothMPC_ubIdx86, smoothMPC_grad_ineq86, smoothMPC_lubbysub86, smoothMPC_llbbyslb86);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub87, smoothMPC_sub87, smoothMPC_riub87, smoothMPC_llb87, smoothMPC_slb87, smoothMPC_rilb87, smoothMPC_lbIdx87, smoothMPC_ubIdx87, smoothMPC_grad_ineq87, smoothMPC_lubbysub87, smoothMPC_llbbyslb87);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub88, smoothMPC_sub88, smoothMPC_riub88, smoothMPC_llb88, smoothMPC_slb88, smoothMPC_rilb88, smoothMPC_lbIdx88, smoothMPC_ubIdx88, smoothMPC_grad_ineq88, smoothMPC_lubbysub88, smoothMPC_llbbyslb88);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub89, smoothMPC_sub89, smoothMPC_riub89, smoothMPC_llb89, smoothMPC_slb89, smoothMPC_rilb89, smoothMPC_lbIdx89, smoothMPC_ubIdx89, smoothMPC_grad_ineq89, smoothMPC_lubbysub89, smoothMPC_llbbyslb89);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub90, smoothMPC_sub90, smoothMPC_riub90, smoothMPC_llb90, smoothMPC_slb90, smoothMPC_rilb90, smoothMPC_lbIdx90, smoothMPC_ubIdx90, smoothMPC_grad_ineq90, smoothMPC_lubbysub90, smoothMPC_llbbyslb90);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub91, smoothMPC_sub91, smoothMPC_riub91, smoothMPC_llb91, smoothMPC_slb91, smoothMPC_rilb91, smoothMPC_lbIdx91, smoothMPC_ubIdx91, smoothMPC_grad_ineq91, smoothMPC_lubbysub91, smoothMPC_llbbyslb91);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub92, smoothMPC_sub92, smoothMPC_riub92, smoothMPC_llb92, smoothMPC_slb92, smoothMPC_rilb92, smoothMPC_lbIdx92, smoothMPC_ubIdx92, smoothMPC_grad_ineq92, smoothMPC_lubbysub92, smoothMPC_llbbyslb92);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub93, smoothMPC_sub93, smoothMPC_riub93, smoothMPC_llb93, smoothMPC_slb93, smoothMPC_rilb93, smoothMPC_lbIdx93, smoothMPC_ubIdx93, smoothMPC_grad_ineq93, smoothMPC_lubbysub93, smoothMPC_llbbyslb93);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub94, smoothMPC_sub94, smoothMPC_riub94, smoothMPC_llb94, smoothMPC_slb94, smoothMPC_rilb94, smoothMPC_lbIdx94, smoothMPC_ubIdx94, smoothMPC_grad_ineq94, smoothMPC_lubbysub94, smoothMPC_llbbyslb94);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub95, smoothMPC_sub95, smoothMPC_riub95, smoothMPC_llb95, smoothMPC_slb95, smoothMPC_rilb95, smoothMPC_lbIdx95, smoothMPC_ubIdx95, smoothMPC_grad_ineq95, smoothMPC_lubbysub95, smoothMPC_llbbyslb95);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub96, smoothMPC_sub96, smoothMPC_riub96, smoothMPC_llb96, smoothMPC_slb96, smoothMPC_rilb96, smoothMPC_lbIdx96, smoothMPC_ubIdx96, smoothMPC_grad_ineq96, smoothMPC_lubbysub96, smoothMPC_llbbyslb96);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub97, smoothMPC_sub97, smoothMPC_riub97, smoothMPC_llb97, smoothMPC_slb97, smoothMPC_rilb97, smoothMPC_lbIdx97, smoothMPC_ubIdx97, smoothMPC_grad_ineq97, smoothMPC_lubbysub97, smoothMPC_llbbyslb97);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub98, smoothMPC_sub98, smoothMPC_riub98, smoothMPC_llb98, smoothMPC_slb98, smoothMPC_rilb98, smoothMPC_lbIdx98, smoothMPC_ubIdx98, smoothMPC_grad_ineq98, smoothMPC_lubbysub98, smoothMPC_llbbyslb98);
smoothMPC_LA_INEQ_B_GRAD_3_3_3(smoothMPC_lub99, smoothMPC_sub99, smoothMPC_riub99, smoothMPC_llb99, smoothMPC_slb99, smoothMPC_rilb99, smoothMPC_lbIdx99, smoothMPC_ubIdx99, smoothMPC_grad_ineq99, smoothMPC_lubbysub99, smoothMPC_llbbyslb99);
info->dobj = info->pobj - info->dgap;
info->rdgap = info->pobj ? info->dgap / info->pobj : 1e6;
if( info->rdgap < 0 ) info->rdgap = -info->rdgap;
PRINTTEXT("  %3d  %3.1e  %3.1e  %+6.4e  %+6.4e  %+3.1e  %3.1e  %3.1e\n",info->it, info->res_eq, info->res_ineq, info->pobj, info->dobj, info->dgap, info->rdgap, info->mu);
if( info->mu < smoothMPC_SET_ACC_KKTCOMPL
    && (info->rdgap < smoothMPC_SET_ACC_RDGAP || info->dgap < smoothMPC_SET_ACC_KKTCOMPL)
    && info->res_eq < smoothMPC_SET_ACC_RESEQ
    && info->res_ineq < smoothMPC_SET_ACC_RESINEQ ){
PRINTTEXT("OPTIMAL (within RESEQ=%2.1e, RESINEQ=%2.1e, (R)DGAP=(%2.1e)%2.1e, MU=%2.1e).\n",smoothMPC_SET_ACC_RESEQ, smoothMPC_SET_ACC_RESINEQ,smoothMPC_SET_ACC_KKTCOMPL,smoothMPC_SET_ACC_RDGAP,smoothMPC_SET_ACC_KKTCOMPL);
exitcode = smoothMPC_OPTIMAL; break; }
if( info->it == smoothMPC_SET_MAXIT ){
PRINTTEXT("Maximum number of iterations reached, exiting.\n");
exitcode = smoothMPC_MAXITREACHED; break; }
smoothMPC_LA_VVADD3_1092(smoothMPC_grad_cost, smoothMPC_grad_eq, smoothMPC_grad_ineq, smoothMPC_rd);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H1, smoothMPC_llbbyslb00, smoothMPC_lbIdx00, smoothMPC_lubbysub00, smoothMPC_ubIdx00, smoothMPC_Phi00);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi00, params->C1, smoothMPC_V00);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi00, smoothMPC_D00, smoothMPC_W00);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W00, smoothMPC_V00, smoothMPC_Ysd01);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi00, smoothMPC_rd00, smoothMPC_Lbyrd00);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H2, smoothMPC_llbbyslb01, smoothMPC_lbIdx01, smoothMPC_lubbysub01, smoothMPC_ubIdx01, smoothMPC_Phi01);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi01, params->C2, smoothMPC_V01);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi01, smoothMPC_D01, smoothMPC_W01);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W01, smoothMPC_V01, smoothMPC_Ysd02);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi01, smoothMPC_rd01, smoothMPC_Lbyrd01);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H3, smoothMPC_llbbyslb02, smoothMPC_lbIdx02, smoothMPC_lubbysub02, smoothMPC_ubIdx02, smoothMPC_Phi02);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi02, params->C3, smoothMPC_V02);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi02, smoothMPC_D01, smoothMPC_W02);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W02, smoothMPC_V02, smoothMPC_Ysd03);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi02, smoothMPC_rd02, smoothMPC_Lbyrd02);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H4, smoothMPC_llbbyslb03, smoothMPC_lbIdx03, smoothMPC_lubbysub03, smoothMPC_ubIdx03, smoothMPC_Phi03);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi03, params->C4, smoothMPC_V03);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi03, smoothMPC_D01, smoothMPC_W03);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W03, smoothMPC_V03, smoothMPC_Ysd04);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi03, smoothMPC_rd03, smoothMPC_Lbyrd03);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H5, smoothMPC_llbbyslb04, smoothMPC_lbIdx04, smoothMPC_lubbysub04, smoothMPC_ubIdx04, smoothMPC_Phi04);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi04, params->C5, smoothMPC_V04);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi04, smoothMPC_D01, smoothMPC_W04);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W04, smoothMPC_V04, smoothMPC_Ysd05);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi04, smoothMPC_rd04, smoothMPC_Lbyrd04);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H6, smoothMPC_llbbyslb05, smoothMPC_lbIdx05, smoothMPC_lubbysub05, smoothMPC_ubIdx05, smoothMPC_Phi05);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi05, params->C6, smoothMPC_V05);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi05, smoothMPC_D01, smoothMPC_W05);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W05, smoothMPC_V05, smoothMPC_Ysd06);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi05, smoothMPC_rd05, smoothMPC_Lbyrd05);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H7, smoothMPC_llbbyslb06, smoothMPC_lbIdx06, smoothMPC_lubbysub06, smoothMPC_ubIdx06, smoothMPC_Phi06);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi06, params->C7, smoothMPC_V06);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi06, smoothMPC_D01, smoothMPC_W06);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W06, smoothMPC_V06, smoothMPC_Ysd07);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi06, smoothMPC_rd06, smoothMPC_Lbyrd06);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H8, smoothMPC_llbbyslb07, smoothMPC_lbIdx07, smoothMPC_lubbysub07, smoothMPC_ubIdx07, smoothMPC_Phi07);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi07, params->C8, smoothMPC_V07);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi07, smoothMPC_D01, smoothMPC_W07);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W07, smoothMPC_V07, smoothMPC_Ysd08);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi07, smoothMPC_rd07, smoothMPC_Lbyrd07);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H9, smoothMPC_llbbyslb08, smoothMPC_lbIdx08, smoothMPC_lubbysub08, smoothMPC_ubIdx08, smoothMPC_Phi08);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi08, params->C9, smoothMPC_V08);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi08, smoothMPC_D01, smoothMPC_W08);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W08, smoothMPC_V08, smoothMPC_Ysd09);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi08, smoothMPC_rd08, smoothMPC_Lbyrd08);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H10, smoothMPC_llbbyslb09, smoothMPC_lbIdx09, smoothMPC_lubbysub09, smoothMPC_ubIdx09, smoothMPC_Phi09);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi09, params->C10, smoothMPC_V09);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi09, smoothMPC_D01, smoothMPC_W09);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W09, smoothMPC_V09, smoothMPC_Ysd10);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi09, smoothMPC_rd09, smoothMPC_Lbyrd09);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H11, smoothMPC_llbbyslb10, smoothMPC_lbIdx10, smoothMPC_lubbysub10, smoothMPC_ubIdx10, smoothMPC_Phi10);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi10, params->C11, smoothMPC_V10);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi10, smoothMPC_D01, smoothMPC_W10);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W10, smoothMPC_V10, smoothMPC_Ysd11);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi10, smoothMPC_rd10, smoothMPC_Lbyrd10);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H12, smoothMPC_llbbyslb11, smoothMPC_lbIdx11, smoothMPC_lubbysub11, smoothMPC_ubIdx11, smoothMPC_Phi11);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi11, params->C12, smoothMPC_V11);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi11, smoothMPC_D01, smoothMPC_W11);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W11, smoothMPC_V11, smoothMPC_Ysd12);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi11, smoothMPC_rd11, smoothMPC_Lbyrd11);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H13, smoothMPC_llbbyslb12, smoothMPC_lbIdx12, smoothMPC_lubbysub12, smoothMPC_ubIdx12, smoothMPC_Phi12);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi12, params->C13, smoothMPC_V12);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi12, smoothMPC_D01, smoothMPC_W12);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W12, smoothMPC_V12, smoothMPC_Ysd13);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi12, smoothMPC_rd12, smoothMPC_Lbyrd12);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H14, smoothMPC_llbbyslb13, smoothMPC_lbIdx13, smoothMPC_lubbysub13, smoothMPC_ubIdx13, smoothMPC_Phi13);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi13, params->C14, smoothMPC_V13);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi13, smoothMPC_D01, smoothMPC_W13);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W13, smoothMPC_V13, smoothMPC_Ysd14);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi13, smoothMPC_rd13, smoothMPC_Lbyrd13);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H15, smoothMPC_llbbyslb14, smoothMPC_lbIdx14, smoothMPC_lubbysub14, smoothMPC_ubIdx14, smoothMPC_Phi14);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi14, params->C15, smoothMPC_V14);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi14, smoothMPC_D01, smoothMPC_W14);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W14, smoothMPC_V14, smoothMPC_Ysd15);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi14, smoothMPC_rd14, smoothMPC_Lbyrd14);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H16, smoothMPC_llbbyslb15, smoothMPC_lbIdx15, smoothMPC_lubbysub15, smoothMPC_ubIdx15, smoothMPC_Phi15);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi15, params->C16, smoothMPC_V15);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi15, smoothMPC_D01, smoothMPC_W15);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W15, smoothMPC_V15, smoothMPC_Ysd16);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi15, smoothMPC_rd15, smoothMPC_Lbyrd15);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H17, smoothMPC_llbbyslb16, smoothMPC_lbIdx16, smoothMPC_lubbysub16, smoothMPC_ubIdx16, smoothMPC_Phi16);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi16, params->C17, smoothMPC_V16);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi16, smoothMPC_D01, smoothMPC_W16);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W16, smoothMPC_V16, smoothMPC_Ysd17);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi16, smoothMPC_rd16, smoothMPC_Lbyrd16);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H18, smoothMPC_llbbyslb17, smoothMPC_lbIdx17, smoothMPC_lubbysub17, smoothMPC_ubIdx17, smoothMPC_Phi17);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi17, params->C18, smoothMPC_V17);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi17, smoothMPC_D01, smoothMPC_W17);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W17, smoothMPC_V17, smoothMPC_Ysd18);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi17, smoothMPC_rd17, smoothMPC_Lbyrd17);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H19, smoothMPC_llbbyslb18, smoothMPC_lbIdx18, smoothMPC_lubbysub18, smoothMPC_ubIdx18, smoothMPC_Phi18);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi18, params->C19, smoothMPC_V18);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi18, smoothMPC_D01, smoothMPC_W18);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W18, smoothMPC_V18, smoothMPC_Ysd19);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi18, smoothMPC_rd18, smoothMPC_Lbyrd18);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H20, smoothMPC_llbbyslb19, smoothMPC_lbIdx19, smoothMPC_lubbysub19, smoothMPC_ubIdx19, smoothMPC_Phi19);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi19, params->C20, smoothMPC_V19);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi19, smoothMPC_D01, smoothMPC_W19);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W19, smoothMPC_V19, smoothMPC_Ysd20);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi19, smoothMPC_rd19, smoothMPC_Lbyrd19);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H21, smoothMPC_llbbyslb20, smoothMPC_lbIdx20, smoothMPC_lubbysub20, smoothMPC_ubIdx20, smoothMPC_Phi20);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi20, params->C21, smoothMPC_V20);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi20, smoothMPC_D01, smoothMPC_W20);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W20, smoothMPC_V20, smoothMPC_Ysd21);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi20, smoothMPC_rd20, smoothMPC_Lbyrd20);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H22, smoothMPC_llbbyslb21, smoothMPC_lbIdx21, smoothMPC_lubbysub21, smoothMPC_ubIdx21, smoothMPC_Phi21);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi21, params->C22, smoothMPC_V21);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi21, smoothMPC_D01, smoothMPC_W21);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W21, smoothMPC_V21, smoothMPC_Ysd22);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi21, smoothMPC_rd21, smoothMPC_Lbyrd21);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H23, smoothMPC_llbbyslb22, smoothMPC_lbIdx22, smoothMPC_lubbysub22, smoothMPC_ubIdx22, smoothMPC_Phi22);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi22, params->C23, smoothMPC_V22);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi22, smoothMPC_D01, smoothMPC_W22);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W22, smoothMPC_V22, smoothMPC_Ysd23);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi22, smoothMPC_rd22, smoothMPC_Lbyrd22);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H24, smoothMPC_llbbyslb23, smoothMPC_lbIdx23, smoothMPC_lubbysub23, smoothMPC_ubIdx23, smoothMPC_Phi23);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi23, params->C24, smoothMPC_V23);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi23, smoothMPC_D01, smoothMPC_W23);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W23, smoothMPC_V23, smoothMPC_Ysd24);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi23, smoothMPC_rd23, smoothMPC_Lbyrd23);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H25, smoothMPC_llbbyslb24, smoothMPC_lbIdx24, smoothMPC_lubbysub24, smoothMPC_ubIdx24, smoothMPC_Phi24);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi24, params->C25, smoothMPC_V24);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi24, smoothMPC_D01, smoothMPC_W24);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W24, smoothMPC_V24, smoothMPC_Ysd25);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi24, smoothMPC_rd24, smoothMPC_Lbyrd24);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H26, smoothMPC_llbbyslb25, smoothMPC_lbIdx25, smoothMPC_lubbysub25, smoothMPC_ubIdx25, smoothMPC_Phi25);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi25, params->C26, smoothMPC_V25);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi25, smoothMPC_D01, smoothMPC_W25);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W25, smoothMPC_V25, smoothMPC_Ysd26);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi25, smoothMPC_rd25, smoothMPC_Lbyrd25);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H27, smoothMPC_llbbyslb26, smoothMPC_lbIdx26, smoothMPC_lubbysub26, smoothMPC_ubIdx26, smoothMPC_Phi26);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi26, params->C27, smoothMPC_V26);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi26, smoothMPC_D01, smoothMPC_W26);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W26, smoothMPC_V26, smoothMPC_Ysd27);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi26, smoothMPC_rd26, smoothMPC_Lbyrd26);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H28, smoothMPC_llbbyslb27, smoothMPC_lbIdx27, smoothMPC_lubbysub27, smoothMPC_ubIdx27, smoothMPC_Phi27);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi27, params->C28, smoothMPC_V27);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi27, smoothMPC_D01, smoothMPC_W27);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W27, smoothMPC_V27, smoothMPC_Ysd28);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi27, smoothMPC_rd27, smoothMPC_Lbyrd27);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H29, smoothMPC_llbbyslb28, smoothMPC_lbIdx28, smoothMPC_lubbysub28, smoothMPC_ubIdx28, smoothMPC_Phi28);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi28, params->C29, smoothMPC_V28);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi28, smoothMPC_D01, smoothMPC_W28);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W28, smoothMPC_V28, smoothMPC_Ysd29);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi28, smoothMPC_rd28, smoothMPC_Lbyrd28);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H30, smoothMPC_llbbyslb29, smoothMPC_lbIdx29, smoothMPC_lubbysub29, smoothMPC_ubIdx29, smoothMPC_Phi29);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi29, params->C30, smoothMPC_V29);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi29, smoothMPC_D01, smoothMPC_W29);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W29, smoothMPC_V29, smoothMPC_Ysd30);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi29, smoothMPC_rd29, smoothMPC_Lbyrd29);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H31, smoothMPC_llbbyslb30, smoothMPC_lbIdx30, smoothMPC_lubbysub30, smoothMPC_ubIdx30, smoothMPC_Phi30);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi30, params->C31, smoothMPC_V30);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi30, smoothMPC_D01, smoothMPC_W30);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W30, smoothMPC_V30, smoothMPC_Ysd31);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi30, smoothMPC_rd30, smoothMPC_Lbyrd30);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H32, smoothMPC_llbbyslb31, smoothMPC_lbIdx31, smoothMPC_lubbysub31, smoothMPC_ubIdx31, smoothMPC_Phi31);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi31, params->C32, smoothMPC_V31);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi31, smoothMPC_D01, smoothMPC_W31);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W31, smoothMPC_V31, smoothMPC_Ysd32);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi31, smoothMPC_rd31, smoothMPC_Lbyrd31);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H33, smoothMPC_llbbyslb32, smoothMPC_lbIdx32, smoothMPC_lubbysub32, smoothMPC_ubIdx32, smoothMPC_Phi32);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi32, params->C33, smoothMPC_V32);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi32, smoothMPC_D01, smoothMPC_W32);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W32, smoothMPC_V32, smoothMPC_Ysd33);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi32, smoothMPC_rd32, smoothMPC_Lbyrd32);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H34, smoothMPC_llbbyslb33, smoothMPC_lbIdx33, smoothMPC_lubbysub33, smoothMPC_ubIdx33, smoothMPC_Phi33);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi33, params->C34, smoothMPC_V33);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi33, smoothMPC_D01, smoothMPC_W33);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W33, smoothMPC_V33, smoothMPC_Ysd34);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi33, smoothMPC_rd33, smoothMPC_Lbyrd33);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H35, smoothMPC_llbbyslb34, smoothMPC_lbIdx34, smoothMPC_lubbysub34, smoothMPC_ubIdx34, smoothMPC_Phi34);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi34, params->C35, smoothMPC_V34);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi34, smoothMPC_D01, smoothMPC_W34);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W34, smoothMPC_V34, smoothMPC_Ysd35);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi34, smoothMPC_rd34, smoothMPC_Lbyrd34);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H36, smoothMPC_llbbyslb35, smoothMPC_lbIdx35, smoothMPC_lubbysub35, smoothMPC_ubIdx35, smoothMPC_Phi35);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi35, params->C36, smoothMPC_V35);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi35, smoothMPC_D01, smoothMPC_W35);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W35, smoothMPC_V35, smoothMPC_Ysd36);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi35, smoothMPC_rd35, smoothMPC_Lbyrd35);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H37, smoothMPC_llbbyslb36, smoothMPC_lbIdx36, smoothMPC_lubbysub36, smoothMPC_ubIdx36, smoothMPC_Phi36);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi36, params->C37, smoothMPC_V36);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi36, smoothMPC_D01, smoothMPC_W36);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W36, smoothMPC_V36, smoothMPC_Ysd37);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi36, smoothMPC_rd36, smoothMPC_Lbyrd36);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H38, smoothMPC_llbbyslb37, smoothMPC_lbIdx37, smoothMPC_lubbysub37, smoothMPC_ubIdx37, smoothMPC_Phi37);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi37, params->C38, smoothMPC_V37);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi37, smoothMPC_D01, smoothMPC_W37);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W37, smoothMPC_V37, smoothMPC_Ysd38);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi37, smoothMPC_rd37, smoothMPC_Lbyrd37);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H39, smoothMPC_llbbyslb38, smoothMPC_lbIdx38, smoothMPC_lubbysub38, smoothMPC_ubIdx38, smoothMPC_Phi38);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi38, params->C39, smoothMPC_V38);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi38, smoothMPC_D01, smoothMPC_W38);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W38, smoothMPC_V38, smoothMPC_Ysd39);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi38, smoothMPC_rd38, smoothMPC_Lbyrd38);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H40, smoothMPC_llbbyslb39, smoothMPC_lbIdx39, smoothMPC_lubbysub39, smoothMPC_ubIdx39, smoothMPC_Phi39);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi39, params->C40, smoothMPC_V39);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi39, smoothMPC_D01, smoothMPC_W39);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W39, smoothMPC_V39, smoothMPC_Ysd40);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi39, smoothMPC_rd39, smoothMPC_Lbyrd39);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H41, smoothMPC_llbbyslb40, smoothMPC_lbIdx40, smoothMPC_lubbysub40, smoothMPC_ubIdx40, smoothMPC_Phi40);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi40, params->C41, smoothMPC_V40);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi40, smoothMPC_D01, smoothMPC_W40);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W40, smoothMPC_V40, smoothMPC_Ysd41);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi40, smoothMPC_rd40, smoothMPC_Lbyrd40);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H42, smoothMPC_llbbyslb41, smoothMPC_lbIdx41, smoothMPC_lubbysub41, smoothMPC_ubIdx41, smoothMPC_Phi41);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi41, params->C42, smoothMPC_V41);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi41, smoothMPC_D01, smoothMPC_W41);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W41, smoothMPC_V41, smoothMPC_Ysd42);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi41, smoothMPC_rd41, smoothMPC_Lbyrd41);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H43, smoothMPC_llbbyslb42, smoothMPC_lbIdx42, smoothMPC_lubbysub42, smoothMPC_ubIdx42, smoothMPC_Phi42);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi42, params->C43, smoothMPC_V42);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi42, smoothMPC_D01, smoothMPC_W42);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W42, smoothMPC_V42, smoothMPC_Ysd43);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi42, smoothMPC_rd42, smoothMPC_Lbyrd42);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H44, smoothMPC_llbbyslb43, smoothMPC_lbIdx43, smoothMPC_lubbysub43, smoothMPC_ubIdx43, smoothMPC_Phi43);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi43, params->C44, smoothMPC_V43);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi43, smoothMPC_D01, smoothMPC_W43);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W43, smoothMPC_V43, smoothMPC_Ysd44);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi43, smoothMPC_rd43, smoothMPC_Lbyrd43);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H45, smoothMPC_llbbyslb44, smoothMPC_lbIdx44, smoothMPC_lubbysub44, smoothMPC_ubIdx44, smoothMPC_Phi44);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi44, params->C45, smoothMPC_V44);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi44, smoothMPC_D01, smoothMPC_W44);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W44, smoothMPC_V44, smoothMPC_Ysd45);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi44, smoothMPC_rd44, smoothMPC_Lbyrd44);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H46, smoothMPC_llbbyslb45, smoothMPC_lbIdx45, smoothMPC_lubbysub45, smoothMPC_ubIdx45, smoothMPC_Phi45);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi45, params->C46, smoothMPC_V45);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi45, smoothMPC_D01, smoothMPC_W45);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W45, smoothMPC_V45, smoothMPC_Ysd46);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi45, smoothMPC_rd45, smoothMPC_Lbyrd45);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H47, smoothMPC_llbbyslb46, smoothMPC_lbIdx46, smoothMPC_lubbysub46, smoothMPC_ubIdx46, smoothMPC_Phi46);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi46, params->C47, smoothMPC_V46);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi46, smoothMPC_D01, smoothMPC_W46);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W46, smoothMPC_V46, smoothMPC_Ysd47);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi46, smoothMPC_rd46, smoothMPC_Lbyrd46);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H48, smoothMPC_llbbyslb47, smoothMPC_lbIdx47, smoothMPC_lubbysub47, smoothMPC_ubIdx47, smoothMPC_Phi47);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi47, params->C48, smoothMPC_V47);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi47, smoothMPC_D01, smoothMPC_W47);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W47, smoothMPC_V47, smoothMPC_Ysd48);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi47, smoothMPC_rd47, smoothMPC_Lbyrd47);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H49, smoothMPC_llbbyslb48, smoothMPC_lbIdx48, smoothMPC_lubbysub48, smoothMPC_ubIdx48, smoothMPC_Phi48);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi48, params->C49, smoothMPC_V48);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi48, smoothMPC_D01, smoothMPC_W48);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W48, smoothMPC_V48, smoothMPC_Ysd49);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi48, smoothMPC_rd48, smoothMPC_Lbyrd48);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H50, smoothMPC_llbbyslb49, smoothMPC_lbIdx49, smoothMPC_lubbysub49, smoothMPC_ubIdx49, smoothMPC_Phi49);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi49, params->C50, smoothMPC_V49);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi49, smoothMPC_D01, smoothMPC_W49);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W49, smoothMPC_V49, smoothMPC_Ysd50);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi49, smoothMPC_rd49, smoothMPC_Lbyrd49);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H51, smoothMPC_llbbyslb50, smoothMPC_lbIdx50, smoothMPC_lubbysub50, smoothMPC_ubIdx50, smoothMPC_Phi50);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi50, params->C51, smoothMPC_V50);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi50, smoothMPC_D01, smoothMPC_W50);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W50, smoothMPC_V50, smoothMPC_Ysd51);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi50, smoothMPC_rd50, smoothMPC_Lbyrd50);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H52, smoothMPC_llbbyslb51, smoothMPC_lbIdx51, smoothMPC_lubbysub51, smoothMPC_ubIdx51, smoothMPC_Phi51);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi51, params->C52, smoothMPC_V51);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi51, smoothMPC_D01, smoothMPC_W51);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W51, smoothMPC_V51, smoothMPC_Ysd52);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi51, smoothMPC_rd51, smoothMPC_Lbyrd51);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H53, smoothMPC_llbbyslb52, smoothMPC_lbIdx52, smoothMPC_lubbysub52, smoothMPC_ubIdx52, smoothMPC_Phi52);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi52, params->C53, smoothMPC_V52);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi52, smoothMPC_D01, smoothMPC_W52);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W52, smoothMPC_V52, smoothMPC_Ysd53);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi52, smoothMPC_rd52, smoothMPC_Lbyrd52);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H54, smoothMPC_llbbyslb53, smoothMPC_lbIdx53, smoothMPC_lubbysub53, smoothMPC_ubIdx53, smoothMPC_Phi53);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi53, params->C54, smoothMPC_V53);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi53, smoothMPC_D01, smoothMPC_W53);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W53, smoothMPC_V53, smoothMPC_Ysd54);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi53, smoothMPC_rd53, smoothMPC_Lbyrd53);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H55, smoothMPC_llbbyslb54, smoothMPC_lbIdx54, smoothMPC_lubbysub54, smoothMPC_ubIdx54, smoothMPC_Phi54);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi54, params->C55, smoothMPC_V54);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi54, smoothMPC_D01, smoothMPC_W54);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W54, smoothMPC_V54, smoothMPC_Ysd55);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi54, smoothMPC_rd54, smoothMPC_Lbyrd54);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H56, smoothMPC_llbbyslb55, smoothMPC_lbIdx55, smoothMPC_lubbysub55, smoothMPC_ubIdx55, smoothMPC_Phi55);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi55, params->C56, smoothMPC_V55);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi55, smoothMPC_D01, smoothMPC_W55);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W55, smoothMPC_V55, smoothMPC_Ysd56);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi55, smoothMPC_rd55, smoothMPC_Lbyrd55);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H57, smoothMPC_llbbyslb56, smoothMPC_lbIdx56, smoothMPC_lubbysub56, smoothMPC_ubIdx56, smoothMPC_Phi56);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi56, params->C57, smoothMPC_V56);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi56, smoothMPC_D01, smoothMPC_W56);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W56, smoothMPC_V56, smoothMPC_Ysd57);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi56, smoothMPC_rd56, smoothMPC_Lbyrd56);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H58, smoothMPC_llbbyslb57, smoothMPC_lbIdx57, smoothMPC_lubbysub57, smoothMPC_ubIdx57, smoothMPC_Phi57);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi57, params->C58, smoothMPC_V57);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi57, smoothMPC_D01, smoothMPC_W57);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W57, smoothMPC_V57, smoothMPC_Ysd58);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi57, smoothMPC_rd57, smoothMPC_Lbyrd57);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H59, smoothMPC_llbbyslb58, smoothMPC_lbIdx58, smoothMPC_lubbysub58, smoothMPC_ubIdx58, smoothMPC_Phi58);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi58, params->C59, smoothMPC_V58);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi58, smoothMPC_D01, smoothMPC_W58);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W58, smoothMPC_V58, smoothMPC_Ysd59);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi58, smoothMPC_rd58, smoothMPC_Lbyrd58);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H60, smoothMPC_llbbyslb59, smoothMPC_lbIdx59, smoothMPC_lubbysub59, smoothMPC_ubIdx59, smoothMPC_Phi59);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi59, params->C60, smoothMPC_V59);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi59, smoothMPC_D01, smoothMPC_W59);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W59, smoothMPC_V59, smoothMPC_Ysd60);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi59, smoothMPC_rd59, smoothMPC_Lbyrd59);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H61, smoothMPC_llbbyslb60, smoothMPC_lbIdx60, smoothMPC_lubbysub60, smoothMPC_ubIdx60, smoothMPC_Phi60);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi60, params->C61, smoothMPC_V60);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi60, smoothMPC_D01, smoothMPC_W60);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W60, smoothMPC_V60, smoothMPC_Ysd61);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi60, smoothMPC_rd60, smoothMPC_Lbyrd60);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H62, smoothMPC_llbbyslb61, smoothMPC_lbIdx61, smoothMPC_lubbysub61, smoothMPC_ubIdx61, smoothMPC_Phi61);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi61, params->C62, smoothMPC_V61);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi61, smoothMPC_D01, smoothMPC_W61);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W61, smoothMPC_V61, smoothMPC_Ysd62);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi61, smoothMPC_rd61, smoothMPC_Lbyrd61);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H63, smoothMPC_llbbyslb62, smoothMPC_lbIdx62, smoothMPC_lubbysub62, smoothMPC_ubIdx62, smoothMPC_Phi62);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi62, params->C63, smoothMPC_V62);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi62, smoothMPC_D01, smoothMPC_W62);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W62, smoothMPC_V62, smoothMPC_Ysd63);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi62, smoothMPC_rd62, smoothMPC_Lbyrd62);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H64, smoothMPC_llbbyslb63, smoothMPC_lbIdx63, smoothMPC_lubbysub63, smoothMPC_ubIdx63, smoothMPC_Phi63);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi63, params->C64, smoothMPC_V63);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi63, smoothMPC_D01, smoothMPC_W63);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W63, smoothMPC_V63, smoothMPC_Ysd64);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi63, smoothMPC_rd63, smoothMPC_Lbyrd63);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H65, smoothMPC_llbbyslb64, smoothMPC_lbIdx64, smoothMPC_lubbysub64, smoothMPC_ubIdx64, smoothMPC_Phi64);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi64, params->C65, smoothMPC_V64);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi64, smoothMPC_D01, smoothMPC_W64);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W64, smoothMPC_V64, smoothMPC_Ysd65);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi64, smoothMPC_rd64, smoothMPC_Lbyrd64);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H66, smoothMPC_llbbyslb65, smoothMPC_lbIdx65, smoothMPC_lubbysub65, smoothMPC_ubIdx65, smoothMPC_Phi65);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi65, params->C66, smoothMPC_V65);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi65, smoothMPC_D01, smoothMPC_W65);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W65, smoothMPC_V65, smoothMPC_Ysd66);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi65, smoothMPC_rd65, smoothMPC_Lbyrd65);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H67, smoothMPC_llbbyslb66, smoothMPC_lbIdx66, smoothMPC_lubbysub66, smoothMPC_ubIdx66, smoothMPC_Phi66);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi66, params->C67, smoothMPC_V66);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi66, smoothMPC_D01, smoothMPC_W66);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W66, smoothMPC_V66, smoothMPC_Ysd67);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi66, smoothMPC_rd66, smoothMPC_Lbyrd66);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H68, smoothMPC_llbbyslb67, smoothMPC_lbIdx67, smoothMPC_lubbysub67, smoothMPC_ubIdx67, smoothMPC_Phi67);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi67, params->C68, smoothMPC_V67);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi67, smoothMPC_D01, smoothMPC_W67);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W67, smoothMPC_V67, smoothMPC_Ysd68);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi67, smoothMPC_rd67, smoothMPC_Lbyrd67);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H69, smoothMPC_llbbyslb68, smoothMPC_lbIdx68, smoothMPC_lubbysub68, smoothMPC_ubIdx68, smoothMPC_Phi68);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi68, params->C69, smoothMPC_V68);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi68, smoothMPC_D01, smoothMPC_W68);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W68, smoothMPC_V68, smoothMPC_Ysd69);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi68, smoothMPC_rd68, smoothMPC_Lbyrd68);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H70, smoothMPC_llbbyslb69, smoothMPC_lbIdx69, smoothMPC_lubbysub69, smoothMPC_ubIdx69, smoothMPC_Phi69);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi69, params->C70, smoothMPC_V69);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi69, smoothMPC_D01, smoothMPC_W69);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W69, smoothMPC_V69, smoothMPC_Ysd70);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi69, smoothMPC_rd69, smoothMPC_Lbyrd69);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H71, smoothMPC_llbbyslb70, smoothMPC_lbIdx70, smoothMPC_lubbysub70, smoothMPC_ubIdx70, smoothMPC_Phi70);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi70, params->C71, smoothMPC_V70);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi70, smoothMPC_D01, smoothMPC_W70);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W70, smoothMPC_V70, smoothMPC_Ysd71);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi70, smoothMPC_rd70, smoothMPC_Lbyrd70);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H72, smoothMPC_llbbyslb71, smoothMPC_lbIdx71, smoothMPC_lubbysub71, smoothMPC_ubIdx71, smoothMPC_Phi71);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi71, params->C72, smoothMPC_V71);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi71, smoothMPC_D01, smoothMPC_W71);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W71, smoothMPC_V71, smoothMPC_Ysd72);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi71, smoothMPC_rd71, smoothMPC_Lbyrd71);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H73, smoothMPC_llbbyslb72, smoothMPC_lbIdx72, smoothMPC_lubbysub72, smoothMPC_ubIdx72, smoothMPC_Phi72);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi72, params->C73, smoothMPC_V72);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi72, smoothMPC_D01, smoothMPC_W72);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W72, smoothMPC_V72, smoothMPC_Ysd73);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi72, smoothMPC_rd72, smoothMPC_Lbyrd72);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H74, smoothMPC_llbbyslb73, smoothMPC_lbIdx73, smoothMPC_lubbysub73, smoothMPC_ubIdx73, smoothMPC_Phi73);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi73, params->C74, smoothMPC_V73);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi73, smoothMPC_D01, smoothMPC_W73);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W73, smoothMPC_V73, smoothMPC_Ysd74);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi73, smoothMPC_rd73, smoothMPC_Lbyrd73);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H75, smoothMPC_llbbyslb74, smoothMPC_lbIdx74, smoothMPC_lubbysub74, smoothMPC_ubIdx74, smoothMPC_Phi74);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi74, params->C75, smoothMPC_V74);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi74, smoothMPC_D01, smoothMPC_W74);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W74, smoothMPC_V74, smoothMPC_Ysd75);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi74, smoothMPC_rd74, smoothMPC_Lbyrd74);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H76, smoothMPC_llbbyslb75, smoothMPC_lbIdx75, smoothMPC_lubbysub75, smoothMPC_ubIdx75, smoothMPC_Phi75);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi75, params->C76, smoothMPC_V75);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi75, smoothMPC_D01, smoothMPC_W75);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W75, smoothMPC_V75, smoothMPC_Ysd76);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi75, smoothMPC_rd75, smoothMPC_Lbyrd75);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H77, smoothMPC_llbbyslb76, smoothMPC_lbIdx76, smoothMPC_lubbysub76, smoothMPC_ubIdx76, smoothMPC_Phi76);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi76, params->C77, smoothMPC_V76);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi76, smoothMPC_D01, smoothMPC_W76);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W76, smoothMPC_V76, smoothMPC_Ysd77);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi76, smoothMPC_rd76, smoothMPC_Lbyrd76);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H78, smoothMPC_llbbyslb77, smoothMPC_lbIdx77, smoothMPC_lubbysub77, smoothMPC_ubIdx77, smoothMPC_Phi77);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi77, params->C78, smoothMPC_V77);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi77, smoothMPC_D01, smoothMPC_W77);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W77, smoothMPC_V77, smoothMPC_Ysd78);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi77, smoothMPC_rd77, smoothMPC_Lbyrd77);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H79, smoothMPC_llbbyslb78, smoothMPC_lbIdx78, smoothMPC_lubbysub78, smoothMPC_ubIdx78, smoothMPC_Phi78);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi78, params->C79, smoothMPC_V78);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi78, smoothMPC_D01, smoothMPC_W78);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W78, smoothMPC_V78, smoothMPC_Ysd79);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi78, smoothMPC_rd78, smoothMPC_Lbyrd78);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H80, smoothMPC_llbbyslb79, smoothMPC_lbIdx79, smoothMPC_lubbysub79, smoothMPC_ubIdx79, smoothMPC_Phi79);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi79, params->C80, smoothMPC_V79);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi79, smoothMPC_D01, smoothMPC_W79);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W79, smoothMPC_V79, smoothMPC_Ysd80);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi79, smoothMPC_rd79, smoothMPC_Lbyrd79);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H81, smoothMPC_llbbyslb80, smoothMPC_lbIdx80, smoothMPC_lubbysub80, smoothMPC_ubIdx80, smoothMPC_Phi80);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi80, params->C81, smoothMPC_V80);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi80, smoothMPC_D01, smoothMPC_W80);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W80, smoothMPC_V80, smoothMPC_Ysd81);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi80, smoothMPC_rd80, smoothMPC_Lbyrd80);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H82, smoothMPC_llbbyslb81, smoothMPC_lbIdx81, smoothMPC_lubbysub81, smoothMPC_ubIdx81, smoothMPC_Phi81);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi81, params->C82, smoothMPC_V81);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi81, smoothMPC_D01, smoothMPC_W81);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W81, smoothMPC_V81, smoothMPC_Ysd82);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi81, smoothMPC_rd81, smoothMPC_Lbyrd81);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H83, smoothMPC_llbbyslb82, smoothMPC_lbIdx82, smoothMPC_lubbysub82, smoothMPC_ubIdx82, smoothMPC_Phi82);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi82, params->C83, smoothMPC_V82);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi82, smoothMPC_D01, smoothMPC_W82);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W82, smoothMPC_V82, smoothMPC_Ysd83);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi82, smoothMPC_rd82, smoothMPC_Lbyrd82);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H84, smoothMPC_llbbyslb83, smoothMPC_lbIdx83, smoothMPC_lubbysub83, smoothMPC_ubIdx83, smoothMPC_Phi83);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi83, params->C84, smoothMPC_V83);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi83, smoothMPC_D01, smoothMPC_W83);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W83, smoothMPC_V83, smoothMPC_Ysd84);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi83, smoothMPC_rd83, smoothMPC_Lbyrd83);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H85, smoothMPC_llbbyslb84, smoothMPC_lbIdx84, smoothMPC_lubbysub84, smoothMPC_ubIdx84, smoothMPC_Phi84);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi84, params->C85, smoothMPC_V84);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi84, smoothMPC_D01, smoothMPC_W84);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W84, smoothMPC_V84, smoothMPC_Ysd85);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi84, smoothMPC_rd84, smoothMPC_Lbyrd84);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H86, smoothMPC_llbbyslb85, smoothMPC_lbIdx85, smoothMPC_lubbysub85, smoothMPC_ubIdx85, smoothMPC_Phi85);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi85, params->C86, smoothMPC_V85);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi85, smoothMPC_D01, smoothMPC_W85);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W85, smoothMPC_V85, smoothMPC_Ysd86);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi85, smoothMPC_rd85, smoothMPC_Lbyrd85);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H87, smoothMPC_llbbyslb86, smoothMPC_lbIdx86, smoothMPC_lubbysub86, smoothMPC_ubIdx86, smoothMPC_Phi86);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi86, params->C87, smoothMPC_V86);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi86, smoothMPC_D01, smoothMPC_W86);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W86, smoothMPC_V86, smoothMPC_Ysd87);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi86, smoothMPC_rd86, smoothMPC_Lbyrd86);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H88, smoothMPC_llbbyslb87, smoothMPC_lbIdx87, smoothMPC_lubbysub87, smoothMPC_ubIdx87, smoothMPC_Phi87);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi87, params->C88, smoothMPC_V87);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi87, smoothMPC_D01, smoothMPC_W87);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W87, smoothMPC_V87, smoothMPC_Ysd88);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi87, smoothMPC_rd87, smoothMPC_Lbyrd87);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H89, smoothMPC_llbbyslb88, smoothMPC_lbIdx88, smoothMPC_lubbysub88, smoothMPC_ubIdx88, smoothMPC_Phi88);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi88, params->C89, smoothMPC_V88);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi88, smoothMPC_D01, smoothMPC_W88);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W88, smoothMPC_V88, smoothMPC_Ysd89);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi88, smoothMPC_rd88, smoothMPC_Lbyrd88);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H90, smoothMPC_llbbyslb89, smoothMPC_lbIdx89, smoothMPC_lubbysub89, smoothMPC_ubIdx89, smoothMPC_Phi89);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi89, params->C90, smoothMPC_V89);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi89, smoothMPC_D01, smoothMPC_W89);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W89, smoothMPC_V89, smoothMPC_Ysd90);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi89, smoothMPC_rd89, smoothMPC_Lbyrd89);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H91, smoothMPC_llbbyslb90, smoothMPC_lbIdx90, smoothMPC_lubbysub90, smoothMPC_ubIdx90, smoothMPC_Phi90);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi90, params->C91, smoothMPC_V90);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi90, smoothMPC_D01, smoothMPC_W90);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W90, smoothMPC_V90, smoothMPC_Ysd91);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi90, smoothMPC_rd90, smoothMPC_Lbyrd90);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H92, smoothMPC_llbbyslb91, smoothMPC_lbIdx91, smoothMPC_lubbysub91, smoothMPC_ubIdx91, smoothMPC_Phi91);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi91, params->C92, smoothMPC_V91);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi91, smoothMPC_D01, smoothMPC_W91);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W91, smoothMPC_V91, smoothMPC_Ysd92);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi91, smoothMPC_rd91, smoothMPC_Lbyrd91);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H93, smoothMPC_llbbyslb92, smoothMPC_lbIdx92, smoothMPC_lubbysub92, smoothMPC_ubIdx92, smoothMPC_Phi92);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi92, params->C93, smoothMPC_V92);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi92, smoothMPC_D01, smoothMPC_W92);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W92, smoothMPC_V92, smoothMPC_Ysd93);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi92, smoothMPC_rd92, smoothMPC_Lbyrd92);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H94, smoothMPC_llbbyslb93, smoothMPC_lbIdx93, smoothMPC_lubbysub93, smoothMPC_ubIdx93, smoothMPC_Phi93);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi93, params->C94, smoothMPC_V93);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi93, smoothMPC_D01, smoothMPC_W93);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W93, smoothMPC_V93, smoothMPC_Ysd94);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi93, smoothMPC_rd93, smoothMPC_Lbyrd93);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H95, smoothMPC_llbbyslb94, smoothMPC_lbIdx94, smoothMPC_lubbysub94, smoothMPC_ubIdx94, smoothMPC_Phi94);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi94, params->C95, smoothMPC_V94);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi94, smoothMPC_D01, smoothMPC_W94);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W94, smoothMPC_V94, smoothMPC_Ysd95);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi94, smoothMPC_rd94, smoothMPC_Lbyrd94);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H96, smoothMPC_llbbyslb95, smoothMPC_lbIdx95, smoothMPC_lubbysub95, smoothMPC_ubIdx95, smoothMPC_Phi95);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi95, params->C96, smoothMPC_V95);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi95, smoothMPC_D01, smoothMPC_W95);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W95, smoothMPC_V95, smoothMPC_Ysd96);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi95, smoothMPC_rd95, smoothMPC_Lbyrd95);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H97, smoothMPC_llbbyslb96, smoothMPC_lbIdx96, smoothMPC_lubbysub96, smoothMPC_ubIdx96, smoothMPC_Phi96);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi96, params->C97, smoothMPC_V96);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi96, smoothMPC_D01, smoothMPC_W96);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W96, smoothMPC_V96, smoothMPC_Ysd97);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi96, smoothMPC_rd96, smoothMPC_Lbyrd96);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H98, smoothMPC_llbbyslb97, smoothMPC_lbIdx97, smoothMPC_lubbysub97, smoothMPC_ubIdx97, smoothMPC_Phi97);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi97, params->C98, smoothMPC_V97);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi97, smoothMPC_D01, smoothMPC_W97);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W97, smoothMPC_V97, smoothMPC_Ysd98);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi97, smoothMPC_rd97, smoothMPC_Lbyrd97);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H99, smoothMPC_llbbyslb98, smoothMPC_lbIdx98, smoothMPC_lubbysub98, smoothMPC_ubIdx98, smoothMPC_Phi98);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi98, params->C99, smoothMPC_V98);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi98, smoothMPC_D01, smoothMPC_W98);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W98, smoothMPC_V98, smoothMPC_Ysd99);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi98, smoothMPC_rd98, smoothMPC_Lbyrd98);
smoothMPC_LA_DIAG_CHOL_ONELOOP_LBUB_3_3_3(params->H100, smoothMPC_llbbyslb99, smoothMPC_lbIdx99, smoothMPC_lubbysub99, smoothMPC_ubIdx99, smoothMPC_Phi99);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_3(smoothMPC_Phi99, smoothMPC_D99, smoothMPC_W99);
smoothMPC_LA_DIAG_FORWARDSUB_3(smoothMPC_Phi99, smoothMPC_rd99, smoothMPC_Lbyrd99);
smoothMPC_LA_DIAGZERO_MMT_3(smoothMPC_W00, smoothMPC_Yd00);
smoothMPC_LA_DIAGZERO_MVMSUB7_3(smoothMPC_W00, smoothMPC_Lbyrd00, smoothMPC_re00, smoothMPC_beta00);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V00, smoothMPC_W01, smoothMPC_Yd01);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V00, smoothMPC_Lbyrd00, smoothMPC_W01, smoothMPC_Lbyrd01, smoothMPC_re01, smoothMPC_beta01);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V01, smoothMPC_W02, smoothMPC_Yd02);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V01, smoothMPC_Lbyrd01, smoothMPC_W02, smoothMPC_Lbyrd02, smoothMPC_re02, smoothMPC_beta02);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V02, smoothMPC_W03, smoothMPC_Yd03);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V02, smoothMPC_Lbyrd02, smoothMPC_W03, smoothMPC_Lbyrd03, smoothMPC_re03, smoothMPC_beta03);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V03, smoothMPC_W04, smoothMPC_Yd04);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V03, smoothMPC_Lbyrd03, smoothMPC_W04, smoothMPC_Lbyrd04, smoothMPC_re04, smoothMPC_beta04);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V04, smoothMPC_W05, smoothMPC_Yd05);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V04, smoothMPC_Lbyrd04, smoothMPC_W05, smoothMPC_Lbyrd05, smoothMPC_re05, smoothMPC_beta05);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V05, smoothMPC_W06, smoothMPC_Yd06);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V05, smoothMPC_Lbyrd05, smoothMPC_W06, smoothMPC_Lbyrd06, smoothMPC_re06, smoothMPC_beta06);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V06, smoothMPC_W07, smoothMPC_Yd07);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V06, smoothMPC_Lbyrd06, smoothMPC_W07, smoothMPC_Lbyrd07, smoothMPC_re07, smoothMPC_beta07);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V07, smoothMPC_W08, smoothMPC_Yd08);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V07, smoothMPC_Lbyrd07, smoothMPC_W08, smoothMPC_Lbyrd08, smoothMPC_re08, smoothMPC_beta08);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V08, smoothMPC_W09, smoothMPC_Yd09);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V08, smoothMPC_Lbyrd08, smoothMPC_W09, smoothMPC_Lbyrd09, smoothMPC_re09, smoothMPC_beta09);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V09, smoothMPC_W10, smoothMPC_Yd10);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V09, smoothMPC_Lbyrd09, smoothMPC_W10, smoothMPC_Lbyrd10, smoothMPC_re10, smoothMPC_beta10);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V10, smoothMPC_W11, smoothMPC_Yd11);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V10, smoothMPC_Lbyrd10, smoothMPC_W11, smoothMPC_Lbyrd11, smoothMPC_re11, smoothMPC_beta11);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V11, smoothMPC_W12, smoothMPC_Yd12);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V11, smoothMPC_Lbyrd11, smoothMPC_W12, smoothMPC_Lbyrd12, smoothMPC_re12, smoothMPC_beta12);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V12, smoothMPC_W13, smoothMPC_Yd13);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V12, smoothMPC_Lbyrd12, smoothMPC_W13, smoothMPC_Lbyrd13, smoothMPC_re13, smoothMPC_beta13);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V13, smoothMPC_W14, smoothMPC_Yd14);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V13, smoothMPC_Lbyrd13, smoothMPC_W14, smoothMPC_Lbyrd14, smoothMPC_re14, smoothMPC_beta14);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V14, smoothMPC_W15, smoothMPC_Yd15);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V14, smoothMPC_Lbyrd14, smoothMPC_W15, smoothMPC_Lbyrd15, smoothMPC_re15, smoothMPC_beta15);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V15, smoothMPC_W16, smoothMPC_Yd16);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V15, smoothMPC_Lbyrd15, smoothMPC_W16, smoothMPC_Lbyrd16, smoothMPC_re16, smoothMPC_beta16);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V16, smoothMPC_W17, smoothMPC_Yd17);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V16, smoothMPC_Lbyrd16, smoothMPC_W17, smoothMPC_Lbyrd17, smoothMPC_re17, smoothMPC_beta17);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V17, smoothMPC_W18, smoothMPC_Yd18);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V17, smoothMPC_Lbyrd17, smoothMPC_W18, smoothMPC_Lbyrd18, smoothMPC_re18, smoothMPC_beta18);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V18, smoothMPC_W19, smoothMPC_Yd19);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V18, smoothMPC_Lbyrd18, smoothMPC_W19, smoothMPC_Lbyrd19, smoothMPC_re19, smoothMPC_beta19);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V19, smoothMPC_W20, smoothMPC_Yd20);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V19, smoothMPC_Lbyrd19, smoothMPC_W20, smoothMPC_Lbyrd20, smoothMPC_re20, smoothMPC_beta20);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V20, smoothMPC_W21, smoothMPC_Yd21);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V20, smoothMPC_Lbyrd20, smoothMPC_W21, smoothMPC_Lbyrd21, smoothMPC_re21, smoothMPC_beta21);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V21, smoothMPC_W22, smoothMPC_Yd22);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V21, smoothMPC_Lbyrd21, smoothMPC_W22, smoothMPC_Lbyrd22, smoothMPC_re22, smoothMPC_beta22);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V22, smoothMPC_W23, smoothMPC_Yd23);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V22, smoothMPC_Lbyrd22, smoothMPC_W23, smoothMPC_Lbyrd23, smoothMPC_re23, smoothMPC_beta23);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V23, smoothMPC_W24, smoothMPC_Yd24);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V23, smoothMPC_Lbyrd23, smoothMPC_W24, smoothMPC_Lbyrd24, smoothMPC_re24, smoothMPC_beta24);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V24, smoothMPC_W25, smoothMPC_Yd25);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V24, smoothMPC_Lbyrd24, smoothMPC_W25, smoothMPC_Lbyrd25, smoothMPC_re25, smoothMPC_beta25);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V25, smoothMPC_W26, smoothMPC_Yd26);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V25, smoothMPC_Lbyrd25, smoothMPC_W26, smoothMPC_Lbyrd26, smoothMPC_re26, smoothMPC_beta26);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V26, smoothMPC_W27, smoothMPC_Yd27);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V26, smoothMPC_Lbyrd26, smoothMPC_W27, smoothMPC_Lbyrd27, smoothMPC_re27, smoothMPC_beta27);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V27, smoothMPC_W28, smoothMPC_Yd28);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V27, smoothMPC_Lbyrd27, smoothMPC_W28, smoothMPC_Lbyrd28, smoothMPC_re28, smoothMPC_beta28);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V28, smoothMPC_W29, smoothMPC_Yd29);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V28, smoothMPC_Lbyrd28, smoothMPC_W29, smoothMPC_Lbyrd29, smoothMPC_re29, smoothMPC_beta29);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V29, smoothMPC_W30, smoothMPC_Yd30);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V29, smoothMPC_Lbyrd29, smoothMPC_W30, smoothMPC_Lbyrd30, smoothMPC_re30, smoothMPC_beta30);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V30, smoothMPC_W31, smoothMPC_Yd31);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V30, smoothMPC_Lbyrd30, smoothMPC_W31, smoothMPC_Lbyrd31, smoothMPC_re31, smoothMPC_beta31);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V31, smoothMPC_W32, smoothMPC_Yd32);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V31, smoothMPC_Lbyrd31, smoothMPC_W32, smoothMPC_Lbyrd32, smoothMPC_re32, smoothMPC_beta32);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V32, smoothMPC_W33, smoothMPC_Yd33);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V32, smoothMPC_Lbyrd32, smoothMPC_W33, smoothMPC_Lbyrd33, smoothMPC_re33, smoothMPC_beta33);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V33, smoothMPC_W34, smoothMPC_Yd34);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V33, smoothMPC_Lbyrd33, smoothMPC_W34, smoothMPC_Lbyrd34, smoothMPC_re34, smoothMPC_beta34);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V34, smoothMPC_W35, smoothMPC_Yd35);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V34, smoothMPC_Lbyrd34, smoothMPC_W35, smoothMPC_Lbyrd35, smoothMPC_re35, smoothMPC_beta35);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V35, smoothMPC_W36, smoothMPC_Yd36);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V35, smoothMPC_Lbyrd35, smoothMPC_W36, smoothMPC_Lbyrd36, smoothMPC_re36, smoothMPC_beta36);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V36, smoothMPC_W37, smoothMPC_Yd37);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V36, smoothMPC_Lbyrd36, smoothMPC_W37, smoothMPC_Lbyrd37, smoothMPC_re37, smoothMPC_beta37);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V37, smoothMPC_W38, smoothMPC_Yd38);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V37, smoothMPC_Lbyrd37, smoothMPC_W38, smoothMPC_Lbyrd38, smoothMPC_re38, smoothMPC_beta38);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V38, smoothMPC_W39, smoothMPC_Yd39);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V38, smoothMPC_Lbyrd38, smoothMPC_W39, smoothMPC_Lbyrd39, smoothMPC_re39, smoothMPC_beta39);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V39, smoothMPC_W40, smoothMPC_Yd40);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V39, smoothMPC_Lbyrd39, smoothMPC_W40, smoothMPC_Lbyrd40, smoothMPC_re40, smoothMPC_beta40);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V40, smoothMPC_W41, smoothMPC_Yd41);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V40, smoothMPC_Lbyrd40, smoothMPC_W41, smoothMPC_Lbyrd41, smoothMPC_re41, smoothMPC_beta41);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V41, smoothMPC_W42, smoothMPC_Yd42);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V41, smoothMPC_Lbyrd41, smoothMPC_W42, smoothMPC_Lbyrd42, smoothMPC_re42, smoothMPC_beta42);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V42, smoothMPC_W43, smoothMPC_Yd43);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V42, smoothMPC_Lbyrd42, smoothMPC_W43, smoothMPC_Lbyrd43, smoothMPC_re43, smoothMPC_beta43);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V43, smoothMPC_W44, smoothMPC_Yd44);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V43, smoothMPC_Lbyrd43, smoothMPC_W44, smoothMPC_Lbyrd44, smoothMPC_re44, smoothMPC_beta44);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V44, smoothMPC_W45, smoothMPC_Yd45);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V44, smoothMPC_Lbyrd44, smoothMPC_W45, smoothMPC_Lbyrd45, smoothMPC_re45, smoothMPC_beta45);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V45, smoothMPC_W46, smoothMPC_Yd46);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V45, smoothMPC_Lbyrd45, smoothMPC_W46, smoothMPC_Lbyrd46, smoothMPC_re46, smoothMPC_beta46);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V46, smoothMPC_W47, smoothMPC_Yd47);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V46, smoothMPC_Lbyrd46, smoothMPC_W47, smoothMPC_Lbyrd47, smoothMPC_re47, smoothMPC_beta47);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V47, smoothMPC_W48, smoothMPC_Yd48);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V47, smoothMPC_Lbyrd47, smoothMPC_W48, smoothMPC_Lbyrd48, smoothMPC_re48, smoothMPC_beta48);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V48, smoothMPC_W49, smoothMPC_Yd49);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V48, smoothMPC_Lbyrd48, smoothMPC_W49, smoothMPC_Lbyrd49, smoothMPC_re49, smoothMPC_beta49);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V49, smoothMPC_W50, smoothMPC_Yd50);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V49, smoothMPC_Lbyrd49, smoothMPC_W50, smoothMPC_Lbyrd50, smoothMPC_re50, smoothMPC_beta50);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V50, smoothMPC_W51, smoothMPC_Yd51);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V50, smoothMPC_Lbyrd50, smoothMPC_W51, smoothMPC_Lbyrd51, smoothMPC_re51, smoothMPC_beta51);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V51, smoothMPC_W52, smoothMPC_Yd52);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V51, smoothMPC_Lbyrd51, smoothMPC_W52, smoothMPC_Lbyrd52, smoothMPC_re52, smoothMPC_beta52);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V52, smoothMPC_W53, smoothMPC_Yd53);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V52, smoothMPC_Lbyrd52, smoothMPC_W53, smoothMPC_Lbyrd53, smoothMPC_re53, smoothMPC_beta53);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V53, smoothMPC_W54, smoothMPC_Yd54);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V53, smoothMPC_Lbyrd53, smoothMPC_W54, smoothMPC_Lbyrd54, smoothMPC_re54, smoothMPC_beta54);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V54, smoothMPC_W55, smoothMPC_Yd55);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V54, smoothMPC_Lbyrd54, smoothMPC_W55, smoothMPC_Lbyrd55, smoothMPC_re55, smoothMPC_beta55);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V55, smoothMPC_W56, smoothMPC_Yd56);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V55, smoothMPC_Lbyrd55, smoothMPC_W56, smoothMPC_Lbyrd56, smoothMPC_re56, smoothMPC_beta56);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V56, smoothMPC_W57, smoothMPC_Yd57);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V56, smoothMPC_Lbyrd56, smoothMPC_W57, smoothMPC_Lbyrd57, smoothMPC_re57, smoothMPC_beta57);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V57, smoothMPC_W58, smoothMPC_Yd58);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V57, smoothMPC_Lbyrd57, smoothMPC_W58, smoothMPC_Lbyrd58, smoothMPC_re58, smoothMPC_beta58);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V58, smoothMPC_W59, smoothMPC_Yd59);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V58, smoothMPC_Lbyrd58, smoothMPC_W59, smoothMPC_Lbyrd59, smoothMPC_re59, smoothMPC_beta59);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V59, smoothMPC_W60, smoothMPC_Yd60);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V59, smoothMPC_Lbyrd59, smoothMPC_W60, smoothMPC_Lbyrd60, smoothMPC_re60, smoothMPC_beta60);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V60, smoothMPC_W61, smoothMPC_Yd61);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V60, smoothMPC_Lbyrd60, smoothMPC_W61, smoothMPC_Lbyrd61, smoothMPC_re61, smoothMPC_beta61);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V61, smoothMPC_W62, smoothMPC_Yd62);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V61, smoothMPC_Lbyrd61, smoothMPC_W62, smoothMPC_Lbyrd62, smoothMPC_re62, smoothMPC_beta62);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V62, smoothMPC_W63, smoothMPC_Yd63);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V62, smoothMPC_Lbyrd62, smoothMPC_W63, smoothMPC_Lbyrd63, smoothMPC_re63, smoothMPC_beta63);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V63, smoothMPC_W64, smoothMPC_Yd64);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V63, smoothMPC_Lbyrd63, smoothMPC_W64, smoothMPC_Lbyrd64, smoothMPC_re64, smoothMPC_beta64);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V64, smoothMPC_W65, smoothMPC_Yd65);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V64, smoothMPC_Lbyrd64, smoothMPC_W65, smoothMPC_Lbyrd65, smoothMPC_re65, smoothMPC_beta65);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V65, smoothMPC_W66, smoothMPC_Yd66);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V65, smoothMPC_Lbyrd65, smoothMPC_W66, smoothMPC_Lbyrd66, smoothMPC_re66, smoothMPC_beta66);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V66, smoothMPC_W67, smoothMPC_Yd67);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V66, smoothMPC_Lbyrd66, smoothMPC_W67, smoothMPC_Lbyrd67, smoothMPC_re67, smoothMPC_beta67);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V67, smoothMPC_W68, smoothMPC_Yd68);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V67, smoothMPC_Lbyrd67, smoothMPC_W68, smoothMPC_Lbyrd68, smoothMPC_re68, smoothMPC_beta68);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V68, smoothMPC_W69, smoothMPC_Yd69);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V68, smoothMPC_Lbyrd68, smoothMPC_W69, smoothMPC_Lbyrd69, smoothMPC_re69, smoothMPC_beta69);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V69, smoothMPC_W70, smoothMPC_Yd70);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V69, smoothMPC_Lbyrd69, smoothMPC_W70, smoothMPC_Lbyrd70, smoothMPC_re70, smoothMPC_beta70);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V70, smoothMPC_W71, smoothMPC_Yd71);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V70, smoothMPC_Lbyrd70, smoothMPC_W71, smoothMPC_Lbyrd71, smoothMPC_re71, smoothMPC_beta71);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V71, smoothMPC_W72, smoothMPC_Yd72);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V71, smoothMPC_Lbyrd71, smoothMPC_W72, smoothMPC_Lbyrd72, smoothMPC_re72, smoothMPC_beta72);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V72, smoothMPC_W73, smoothMPC_Yd73);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V72, smoothMPC_Lbyrd72, smoothMPC_W73, smoothMPC_Lbyrd73, smoothMPC_re73, smoothMPC_beta73);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V73, smoothMPC_W74, smoothMPC_Yd74);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V73, smoothMPC_Lbyrd73, smoothMPC_W74, smoothMPC_Lbyrd74, smoothMPC_re74, smoothMPC_beta74);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V74, smoothMPC_W75, smoothMPC_Yd75);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V74, smoothMPC_Lbyrd74, smoothMPC_W75, smoothMPC_Lbyrd75, smoothMPC_re75, smoothMPC_beta75);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V75, smoothMPC_W76, smoothMPC_Yd76);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V75, smoothMPC_Lbyrd75, smoothMPC_W76, smoothMPC_Lbyrd76, smoothMPC_re76, smoothMPC_beta76);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V76, smoothMPC_W77, smoothMPC_Yd77);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V76, smoothMPC_Lbyrd76, smoothMPC_W77, smoothMPC_Lbyrd77, smoothMPC_re77, smoothMPC_beta77);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V77, smoothMPC_W78, smoothMPC_Yd78);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V77, smoothMPC_Lbyrd77, smoothMPC_W78, smoothMPC_Lbyrd78, smoothMPC_re78, smoothMPC_beta78);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V78, smoothMPC_W79, smoothMPC_Yd79);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V78, smoothMPC_Lbyrd78, smoothMPC_W79, smoothMPC_Lbyrd79, smoothMPC_re79, smoothMPC_beta79);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V79, smoothMPC_W80, smoothMPC_Yd80);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V79, smoothMPC_Lbyrd79, smoothMPC_W80, smoothMPC_Lbyrd80, smoothMPC_re80, smoothMPC_beta80);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V80, smoothMPC_W81, smoothMPC_Yd81);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V80, smoothMPC_Lbyrd80, smoothMPC_W81, smoothMPC_Lbyrd81, smoothMPC_re81, smoothMPC_beta81);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V81, smoothMPC_W82, smoothMPC_Yd82);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V81, smoothMPC_Lbyrd81, smoothMPC_W82, smoothMPC_Lbyrd82, smoothMPC_re82, smoothMPC_beta82);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V82, smoothMPC_W83, smoothMPC_Yd83);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V82, smoothMPC_Lbyrd82, smoothMPC_W83, smoothMPC_Lbyrd83, smoothMPC_re83, smoothMPC_beta83);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V83, smoothMPC_W84, smoothMPC_Yd84);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V83, smoothMPC_Lbyrd83, smoothMPC_W84, smoothMPC_Lbyrd84, smoothMPC_re84, smoothMPC_beta84);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V84, smoothMPC_W85, smoothMPC_Yd85);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V84, smoothMPC_Lbyrd84, smoothMPC_W85, smoothMPC_Lbyrd85, smoothMPC_re85, smoothMPC_beta85);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V85, smoothMPC_W86, smoothMPC_Yd86);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V85, smoothMPC_Lbyrd85, smoothMPC_W86, smoothMPC_Lbyrd86, smoothMPC_re86, smoothMPC_beta86);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V86, smoothMPC_W87, smoothMPC_Yd87);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V86, smoothMPC_Lbyrd86, smoothMPC_W87, smoothMPC_Lbyrd87, smoothMPC_re87, smoothMPC_beta87);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V87, smoothMPC_W88, smoothMPC_Yd88);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V87, smoothMPC_Lbyrd87, smoothMPC_W88, smoothMPC_Lbyrd88, smoothMPC_re88, smoothMPC_beta88);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V88, smoothMPC_W89, smoothMPC_Yd89);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V88, smoothMPC_Lbyrd88, smoothMPC_W89, smoothMPC_Lbyrd89, smoothMPC_re89, smoothMPC_beta89);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V89, smoothMPC_W90, smoothMPC_Yd90);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V89, smoothMPC_Lbyrd89, smoothMPC_W90, smoothMPC_Lbyrd90, smoothMPC_re90, smoothMPC_beta90);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V90, smoothMPC_W91, smoothMPC_Yd91);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V90, smoothMPC_Lbyrd90, smoothMPC_W91, smoothMPC_Lbyrd91, smoothMPC_re91, smoothMPC_beta91);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V91, smoothMPC_W92, smoothMPC_Yd92);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V91, smoothMPC_Lbyrd91, smoothMPC_W92, smoothMPC_Lbyrd92, smoothMPC_re92, smoothMPC_beta92);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V92, smoothMPC_W93, smoothMPC_Yd93);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V92, smoothMPC_Lbyrd92, smoothMPC_W93, smoothMPC_Lbyrd93, smoothMPC_re93, smoothMPC_beta93);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V93, smoothMPC_W94, smoothMPC_Yd94);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V93, smoothMPC_Lbyrd93, smoothMPC_W94, smoothMPC_Lbyrd94, smoothMPC_re94, smoothMPC_beta94);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V94, smoothMPC_W95, smoothMPC_Yd95);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V94, smoothMPC_Lbyrd94, smoothMPC_W95, smoothMPC_Lbyrd95, smoothMPC_re95, smoothMPC_beta95);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V95, smoothMPC_W96, smoothMPC_Yd96);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V95, smoothMPC_Lbyrd95, smoothMPC_W96, smoothMPC_Lbyrd96, smoothMPC_re96, smoothMPC_beta96);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V96, smoothMPC_W97, smoothMPC_Yd97);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V96, smoothMPC_Lbyrd96, smoothMPC_W97, smoothMPC_Lbyrd97, smoothMPC_re97, smoothMPC_beta97);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V97, smoothMPC_W98, smoothMPC_Yd98);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V97, smoothMPC_Lbyrd97, smoothMPC_W98, smoothMPC_Lbyrd98, smoothMPC_re98, smoothMPC_beta98);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_3(smoothMPC_V98, smoothMPC_W99, smoothMPC_Yd99);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_3(smoothMPC_V98, smoothMPC_Lbyrd98, smoothMPC_W99, smoothMPC_Lbyrd99, smoothMPC_re99, smoothMPC_beta99);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd00, smoothMPC_Ld00);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld00, smoothMPC_beta00, smoothMPC_yy00);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld00, smoothMPC_Ysd01, smoothMPC_Lsd01);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd01, smoothMPC_Yd01);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd01, smoothMPC_Ld01);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd01, smoothMPC_yy00, smoothMPC_beta01, smoothMPC_bmy01);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld01, smoothMPC_bmy01, smoothMPC_yy01);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld01, smoothMPC_Ysd02, smoothMPC_Lsd02);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd02, smoothMPC_Yd02);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd02, smoothMPC_Ld02);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd02, smoothMPC_yy01, smoothMPC_beta02, smoothMPC_bmy02);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld02, smoothMPC_bmy02, smoothMPC_yy02);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld02, smoothMPC_Ysd03, smoothMPC_Lsd03);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd03, smoothMPC_Yd03);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd03, smoothMPC_Ld03);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd03, smoothMPC_yy02, smoothMPC_beta03, smoothMPC_bmy03);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld03, smoothMPC_bmy03, smoothMPC_yy03);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld03, smoothMPC_Ysd04, smoothMPC_Lsd04);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd04, smoothMPC_Yd04);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd04, smoothMPC_Ld04);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd04, smoothMPC_yy03, smoothMPC_beta04, smoothMPC_bmy04);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld04, smoothMPC_bmy04, smoothMPC_yy04);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld04, smoothMPC_Ysd05, smoothMPC_Lsd05);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd05, smoothMPC_Yd05);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd05, smoothMPC_Ld05);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd05, smoothMPC_yy04, smoothMPC_beta05, smoothMPC_bmy05);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld05, smoothMPC_bmy05, smoothMPC_yy05);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld05, smoothMPC_Ysd06, smoothMPC_Lsd06);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd06, smoothMPC_Yd06);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd06, smoothMPC_Ld06);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd06, smoothMPC_yy05, smoothMPC_beta06, smoothMPC_bmy06);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld06, smoothMPC_bmy06, smoothMPC_yy06);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld06, smoothMPC_Ysd07, smoothMPC_Lsd07);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd07, smoothMPC_Yd07);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd07, smoothMPC_Ld07);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd07, smoothMPC_yy06, smoothMPC_beta07, smoothMPC_bmy07);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld07, smoothMPC_bmy07, smoothMPC_yy07);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld07, smoothMPC_Ysd08, smoothMPC_Lsd08);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd08, smoothMPC_Yd08);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd08, smoothMPC_Ld08);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd08, smoothMPC_yy07, smoothMPC_beta08, smoothMPC_bmy08);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld08, smoothMPC_bmy08, smoothMPC_yy08);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld08, smoothMPC_Ysd09, smoothMPC_Lsd09);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd09, smoothMPC_Yd09);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd09, smoothMPC_Ld09);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd09, smoothMPC_yy08, smoothMPC_beta09, smoothMPC_bmy09);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld09, smoothMPC_bmy09, smoothMPC_yy09);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld09, smoothMPC_Ysd10, smoothMPC_Lsd10);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd10, smoothMPC_Yd10);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd10, smoothMPC_Ld10);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd10, smoothMPC_yy09, smoothMPC_beta10, smoothMPC_bmy10);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld10, smoothMPC_bmy10, smoothMPC_yy10);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld10, smoothMPC_Ysd11, smoothMPC_Lsd11);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd11, smoothMPC_Yd11);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd11, smoothMPC_Ld11);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd11, smoothMPC_yy10, smoothMPC_beta11, smoothMPC_bmy11);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld11, smoothMPC_bmy11, smoothMPC_yy11);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld11, smoothMPC_Ysd12, smoothMPC_Lsd12);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd12, smoothMPC_Yd12);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd12, smoothMPC_Ld12);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd12, smoothMPC_yy11, smoothMPC_beta12, smoothMPC_bmy12);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld12, smoothMPC_bmy12, smoothMPC_yy12);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld12, smoothMPC_Ysd13, smoothMPC_Lsd13);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd13, smoothMPC_Yd13);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd13, smoothMPC_Ld13);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd13, smoothMPC_yy12, smoothMPC_beta13, smoothMPC_bmy13);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld13, smoothMPC_bmy13, smoothMPC_yy13);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld13, smoothMPC_Ysd14, smoothMPC_Lsd14);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd14, smoothMPC_Yd14);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd14, smoothMPC_Ld14);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd14, smoothMPC_yy13, smoothMPC_beta14, smoothMPC_bmy14);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld14, smoothMPC_bmy14, smoothMPC_yy14);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld14, smoothMPC_Ysd15, smoothMPC_Lsd15);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd15, smoothMPC_Yd15);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd15, smoothMPC_Ld15);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd15, smoothMPC_yy14, smoothMPC_beta15, smoothMPC_bmy15);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld15, smoothMPC_bmy15, smoothMPC_yy15);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld15, smoothMPC_Ysd16, smoothMPC_Lsd16);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd16, smoothMPC_Yd16);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd16, smoothMPC_Ld16);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd16, smoothMPC_yy15, smoothMPC_beta16, smoothMPC_bmy16);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld16, smoothMPC_bmy16, smoothMPC_yy16);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld16, smoothMPC_Ysd17, smoothMPC_Lsd17);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd17, smoothMPC_Yd17);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd17, smoothMPC_Ld17);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd17, smoothMPC_yy16, smoothMPC_beta17, smoothMPC_bmy17);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld17, smoothMPC_bmy17, smoothMPC_yy17);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld17, smoothMPC_Ysd18, smoothMPC_Lsd18);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd18, smoothMPC_Yd18);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd18, smoothMPC_Ld18);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd18, smoothMPC_yy17, smoothMPC_beta18, smoothMPC_bmy18);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld18, smoothMPC_bmy18, smoothMPC_yy18);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld18, smoothMPC_Ysd19, smoothMPC_Lsd19);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd19, smoothMPC_Yd19);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd19, smoothMPC_Ld19);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd19, smoothMPC_yy18, smoothMPC_beta19, smoothMPC_bmy19);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld19, smoothMPC_bmy19, smoothMPC_yy19);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld19, smoothMPC_Ysd20, smoothMPC_Lsd20);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd20, smoothMPC_Yd20);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd20, smoothMPC_Ld20);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd20, smoothMPC_yy19, smoothMPC_beta20, smoothMPC_bmy20);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld20, smoothMPC_bmy20, smoothMPC_yy20);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld20, smoothMPC_Ysd21, smoothMPC_Lsd21);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd21, smoothMPC_Yd21);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd21, smoothMPC_Ld21);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd21, smoothMPC_yy20, smoothMPC_beta21, smoothMPC_bmy21);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld21, smoothMPC_bmy21, smoothMPC_yy21);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld21, smoothMPC_Ysd22, smoothMPC_Lsd22);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd22, smoothMPC_Yd22);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd22, smoothMPC_Ld22);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd22, smoothMPC_yy21, smoothMPC_beta22, smoothMPC_bmy22);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld22, smoothMPC_bmy22, smoothMPC_yy22);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld22, smoothMPC_Ysd23, smoothMPC_Lsd23);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd23, smoothMPC_Yd23);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd23, smoothMPC_Ld23);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd23, smoothMPC_yy22, smoothMPC_beta23, smoothMPC_bmy23);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld23, smoothMPC_bmy23, smoothMPC_yy23);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld23, smoothMPC_Ysd24, smoothMPC_Lsd24);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd24, smoothMPC_Yd24);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd24, smoothMPC_Ld24);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd24, smoothMPC_yy23, smoothMPC_beta24, smoothMPC_bmy24);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld24, smoothMPC_bmy24, smoothMPC_yy24);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld24, smoothMPC_Ysd25, smoothMPC_Lsd25);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd25, smoothMPC_Yd25);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd25, smoothMPC_Ld25);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd25, smoothMPC_yy24, smoothMPC_beta25, smoothMPC_bmy25);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld25, smoothMPC_bmy25, smoothMPC_yy25);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld25, smoothMPC_Ysd26, smoothMPC_Lsd26);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd26, smoothMPC_Yd26);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd26, smoothMPC_Ld26);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd26, smoothMPC_yy25, smoothMPC_beta26, smoothMPC_bmy26);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld26, smoothMPC_bmy26, smoothMPC_yy26);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld26, smoothMPC_Ysd27, smoothMPC_Lsd27);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd27, smoothMPC_Yd27);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd27, smoothMPC_Ld27);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd27, smoothMPC_yy26, smoothMPC_beta27, smoothMPC_bmy27);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld27, smoothMPC_bmy27, smoothMPC_yy27);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld27, smoothMPC_Ysd28, smoothMPC_Lsd28);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd28, smoothMPC_Yd28);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd28, smoothMPC_Ld28);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd28, smoothMPC_yy27, smoothMPC_beta28, smoothMPC_bmy28);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld28, smoothMPC_bmy28, smoothMPC_yy28);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld28, smoothMPC_Ysd29, smoothMPC_Lsd29);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd29, smoothMPC_Yd29);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd29, smoothMPC_Ld29);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd29, smoothMPC_yy28, smoothMPC_beta29, smoothMPC_bmy29);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld29, smoothMPC_bmy29, smoothMPC_yy29);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld29, smoothMPC_Ysd30, smoothMPC_Lsd30);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd30, smoothMPC_Yd30);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd30, smoothMPC_Ld30);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd30, smoothMPC_yy29, smoothMPC_beta30, smoothMPC_bmy30);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld30, smoothMPC_bmy30, smoothMPC_yy30);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld30, smoothMPC_Ysd31, smoothMPC_Lsd31);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd31, smoothMPC_Yd31);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd31, smoothMPC_Ld31);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd31, smoothMPC_yy30, smoothMPC_beta31, smoothMPC_bmy31);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld31, smoothMPC_bmy31, smoothMPC_yy31);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld31, smoothMPC_Ysd32, smoothMPC_Lsd32);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd32, smoothMPC_Yd32);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd32, smoothMPC_Ld32);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd32, smoothMPC_yy31, smoothMPC_beta32, smoothMPC_bmy32);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld32, smoothMPC_bmy32, smoothMPC_yy32);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld32, smoothMPC_Ysd33, smoothMPC_Lsd33);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd33, smoothMPC_Yd33);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd33, smoothMPC_Ld33);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd33, smoothMPC_yy32, smoothMPC_beta33, smoothMPC_bmy33);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld33, smoothMPC_bmy33, smoothMPC_yy33);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld33, smoothMPC_Ysd34, smoothMPC_Lsd34);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd34, smoothMPC_Yd34);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd34, smoothMPC_Ld34);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd34, smoothMPC_yy33, smoothMPC_beta34, smoothMPC_bmy34);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld34, smoothMPC_bmy34, smoothMPC_yy34);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld34, smoothMPC_Ysd35, smoothMPC_Lsd35);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd35, smoothMPC_Yd35);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd35, smoothMPC_Ld35);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd35, smoothMPC_yy34, smoothMPC_beta35, smoothMPC_bmy35);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld35, smoothMPC_bmy35, smoothMPC_yy35);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld35, smoothMPC_Ysd36, smoothMPC_Lsd36);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd36, smoothMPC_Yd36);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd36, smoothMPC_Ld36);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd36, smoothMPC_yy35, smoothMPC_beta36, smoothMPC_bmy36);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld36, smoothMPC_bmy36, smoothMPC_yy36);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld36, smoothMPC_Ysd37, smoothMPC_Lsd37);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd37, smoothMPC_Yd37);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd37, smoothMPC_Ld37);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd37, smoothMPC_yy36, smoothMPC_beta37, smoothMPC_bmy37);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld37, smoothMPC_bmy37, smoothMPC_yy37);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld37, smoothMPC_Ysd38, smoothMPC_Lsd38);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd38, smoothMPC_Yd38);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd38, smoothMPC_Ld38);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd38, smoothMPC_yy37, smoothMPC_beta38, smoothMPC_bmy38);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld38, smoothMPC_bmy38, smoothMPC_yy38);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld38, smoothMPC_Ysd39, smoothMPC_Lsd39);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd39, smoothMPC_Yd39);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd39, smoothMPC_Ld39);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd39, smoothMPC_yy38, smoothMPC_beta39, smoothMPC_bmy39);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld39, smoothMPC_bmy39, smoothMPC_yy39);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld39, smoothMPC_Ysd40, smoothMPC_Lsd40);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd40, smoothMPC_Yd40);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd40, smoothMPC_Ld40);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd40, smoothMPC_yy39, smoothMPC_beta40, smoothMPC_bmy40);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld40, smoothMPC_bmy40, smoothMPC_yy40);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld40, smoothMPC_Ysd41, smoothMPC_Lsd41);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd41, smoothMPC_Yd41);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd41, smoothMPC_Ld41);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd41, smoothMPC_yy40, smoothMPC_beta41, smoothMPC_bmy41);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld41, smoothMPC_bmy41, smoothMPC_yy41);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld41, smoothMPC_Ysd42, smoothMPC_Lsd42);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd42, smoothMPC_Yd42);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd42, smoothMPC_Ld42);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd42, smoothMPC_yy41, smoothMPC_beta42, smoothMPC_bmy42);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld42, smoothMPC_bmy42, smoothMPC_yy42);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld42, smoothMPC_Ysd43, smoothMPC_Lsd43);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd43, smoothMPC_Yd43);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd43, smoothMPC_Ld43);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd43, smoothMPC_yy42, smoothMPC_beta43, smoothMPC_bmy43);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld43, smoothMPC_bmy43, smoothMPC_yy43);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld43, smoothMPC_Ysd44, smoothMPC_Lsd44);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd44, smoothMPC_Yd44);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd44, smoothMPC_Ld44);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd44, smoothMPC_yy43, smoothMPC_beta44, smoothMPC_bmy44);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld44, smoothMPC_bmy44, smoothMPC_yy44);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld44, smoothMPC_Ysd45, smoothMPC_Lsd45);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd45, smoothMPC_Yd45);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd45, smoothMPC_Ld45);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd45, smoothMPC_yy44, smoothMPC_beta45, smoothMPC_bmy45);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld45, smoothMPC_bmy45, smoothMPC_yy45);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld45, smoothMPC_Ysd46, smoothMPC_Lsd46);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd46, smoothMPC_Yd46);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd46, smoothMPC_Ld46);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd46, smoothMPC_yy45, smoothMPC_beta46, smoothMPC_bmy46);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld46, smoothMPC_bmy46, smoothMPC_yy46);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld46, smoothMPC_Ysd47, smoothMPC_Lsd47);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd47, smoothMPC_Yd47);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd47, smoothMPC_Ld47);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd47, smoothMPC_yy46, smoothMPC_beta47, smoothMPC_bmy47);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld47, smoothMPC_bmy47, smoothMPC_yy47);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld47, smoothMPC_Ysd48, smoothMPC_Lsd48);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd48, smoothMPC_Yd48);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd48, smoothMPC_Ld48);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd48, smoothMPC_yy47, smoothMPC_beta48, smoothMPC_bmy48);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld48, smoothMPC_bmy48, smoothMPC_yy48);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld48, smoothMPC_Ysd49, smoothMPC_Lsd49);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd49, smoothMPC_Yd49);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd49, smoothMPC_Ld49);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd49, smoothMPC_yy48, smoothMPC_beta49, smoothMPC_bmy49);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld49, smoothMPC_bmy49, smoothMPC_yy49);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld49, smoothMPC_Ysd50, smoothMPC_Lsd50);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd50, smoothMPC_Yd50);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd50, smoothMPC_Ld50);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd50, smoothMPC_yy49, smoothMPC_beta50, smoothMPC_bmy50);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld50, smoothMPC_bmy50, smoothMPC_yy50);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld50, smoothMPC_Ysd51, smoothMPC_Lsd51);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd51, smoothMPC_Yd51);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd51, smoothMPC_Ld51);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd51, smoothMPC_yy50, smoothMPC_beta51, smoothMPC_bmy51);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld51, smoothMPC_bmy51, smoothMPC_yy51);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld51, smoothMPC_Ysd52, smoothMPC_Lsd52);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd52, smoothMPC_Yd52);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd52, smoothMPC_Ld52);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd52, smoothMPC_yy51, smoothMPC_beta52, smoothMPC_bmy52);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld52, smoothMPC_bmy52, smoothMPC_yy52);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld52, smoothMPC_Ysd53, smoothMPC_Lsd53);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd53, smoothMPC_Yd53);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd53, smoothMPC_Ld53);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd53, smoothMPC_yy52, smoothMPC_beta53, smoothMPC_bmy53);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld53, smoothMPC_bmy53, smoothMPC_yy53);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld53, smoothMPC_Ysd54, smoothMPC_Lsd54);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd54, smoothMPC_Yd54);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd54, smoothMPC_Ld54);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd54, smoothMPC_yy53, smoothMPC_beta54, smoothMPC_bmy54);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld54, smoothMPC_bmy54, smoothMPC_yy54);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld54, smoothMPC_Ysd55, smoothMPC_Lsd55);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd55, smoothMPC_Yd55);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd55, smoothMPC_Ld55);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd55, smoothMPC_yy54, smoothMPC_beta55, smoothMPC_bmy55);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld55, smoothMPC_bmy55, smoothMPC_yy55);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld55, smoothMPC_Ysd56, smoothMPC_Lsd56);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd56, smoothMPC_Yd56);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd56, smoothMPC_Ld56);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd56, smoothMPC_yy55, smoothMPC_beta56, smoothMPC_bmy56);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld56, smoothMPC_bmy56, smoothMPC_yy56);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld56, smoothMPC_Ysd57, smoothMPC_Lsd57);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd57, smoothMPC_Yd57);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd57, smoothMPC_Ld57);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd57, smoothMPC_yy56, smoothMPC_beta57, smoothMPC_bmy57);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld57, smoothMPC_bmy57, smoothMPC_yy57);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld57, smoothMPC_Ysd58, smoothMPC_Lsd58);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd58, smoothMPC_Yd58);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd58, smoothMPC_Ld58);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd58, smoothMPC_yy57, smoothMPC_beta58, smoothMPC_bmy58);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld58, smoothMPC_bmy58, smoothMPC_yy58);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld58, smoothMPC_Ysd59, smoothMPC_Lsd59);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd59, smoothMPC_Yd59);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd59, smoothMPC_Ld59);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd59, smoothMPC_yy58, smoothMPC_beta59, smoothMPC_bmy59);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld59, smoothMPC_bmy59, smoothMPC_yy59);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld59, smoothMPC_Ysd60, smoothMPC_Lsd60);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd60, smoothMPC_Yd60);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd60, smoothMPC_Ld60);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd60, smoothMPC_yy59, smoothMPC_beta60, smoothMPC_bmy60);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld60, smoothMPC_bmy60, smoothMPC_yy60);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld60, smoothMPC_Ysd61, smoothMPC_Lsd61);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd61, smoothMPC_Yd61);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd61, smoothMPC_Ld61);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd61, smoothMPC_yy60, smoothMPC_beta61, smoothMPC_bmy61);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld61, smoothMPC_bmy61, smoothMPC_yy61);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld61, smoothMPC_Ysd62, smoothMPC_Lsd62);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd62, smoothMPC_Yd62);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd62, smoothMPC_Ld62);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd62, smoothMPC_yy61, smoothMPC_beta62, smoothMPC_bmy62);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld62, smoothMPC_bmy62, smoothMPC_yy62);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld62, smoothMPC_Ysd63, smoothMPC_Lsd63);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd63, smoothMPC_Yd63);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd63, smoothMPC_Ld63);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd63, smoothMPC_yy62, smoothMPC_beta63, smoothMPC_bmy63);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld63, smoothMPC_bmy63, smoothMPC_yy63);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld63, smoothMPC_Ysd64, smoothMPC_Lsd64);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd64, smoothMPC_Yd64);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd64, smoothMPC_Ld64);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd64, smoothMPC_yy63, smoothMPC_beta64, smoothMPC_bmy64);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld64, smoothMPC_bmy64, smoothMPC_yy64);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld64, smoothMPC_Ysd65, smoothMPC_Lsd65);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd65, smoothMPC_Yd65);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd65, smoothMPC_Ld65);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd65, smoothMPC_yy64, smoothMPC_beta65, smoothMPC_bmy65);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld65, smoothMPC_bmy65, smoothMPC_yy65);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld65, smoothMPC_Ysd66, smoothMPC_Lsd66);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd66, smoothMPC_Yd66);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd66, smoothMPC_Ld66);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd66, smoothMPC_yy65, smoothMPC_beta66, smoothMPC_bmy66);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld66, smoothMPC_bmy66, smoothMPC_yy66);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld66, smoothMPC_Ysd67, smoothMPC_Lsd67);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd67, smoothMPC_Yd67);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd67, smoothMPC_Ld67);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd67, smoothMPC_yy66, smoothMPC_beta67, smoothMPC_bmy67);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld67, smoothMPC_bmy67, smoothMPC_yy67);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld67, smoothMPC_Ysd68, smoothMPC_Lsd68);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd68, smoothMPC_Yd68);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd68, smoothMPC_Ld68);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd68, smoothMPC_yy67, smoothMPC_beta68, smoothMPC_bmy68);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld68, smoothMPC_bmy68, smoothMPC_yy68);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld68, smoothMPC_Ysd69, smoothMPC_Lsd69);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd69, smoothMPC_Yd69);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd69, smoothMPC_Ld69);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd69, smoothMPC_yy68, smoothMPC_beta69, smoothMPC_bmy69);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld69, smoothMPC_bmy69, smoothMPC_yy69);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld69, smoothMPC_Ysd70, smoothMPC_Lsd70);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd70, smoothMPC_Yd70);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd70, smoothMPC_Ld70);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd70, smoothMPC_yy69, smoothMPC_beta70, smoothMPC_bmy70);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld70, smoothMPC_bmy70, smoothMPC_yy70);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld70, smoothMPC_Ysd71, smoothMPC_Lsd71);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd71, smoothMPC_Yd71);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd71, smoothMPC_Ld71);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd71, smoothMPC_yy70, smoothMPC_beta71, smoothMPC_bmy71);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld71, smoothMPC_bmy71, smoothMPC_yy71);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld71, smoothMPC_Ysd72, smoothMPC_Lsd72);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd72, smoothMPC_Yd72);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd72, smoothMPC_Ld72);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd72, smoothMPC_yy71, smoothMPC_beta72, smoothMPC_bmy72);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld72, smoothMPC_bmy72, smoothMPC_yy72);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld72, smoothMPC_Ysd73, smoothMPC_Lsd73);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd73, smoothMPC_Yd73);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd73, smoothMPC_Ld73);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd73, smoothMPC_yy72, smoothMPC_beta73, smoothMPC_bmy73);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld73, smoothMPC_bmy73, smoothMPC_yy73);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld73, smoothMPC_Ysd74, smoothMPC_Lsd74);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd74, smoothMPC_Yd74);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd74, smoothMPC_Ld74);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd74, smoothMPC_yy73, smoothMPC_beta74, smoothMPC_bmy74);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld74, smoothMPC_bmy74, smoothMPC_yy74);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld74, smoothMPC_Ysd75, smoothMPC_Lsd75);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd75, smoothMPC_Yd75);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd75, smoothMPC_Ld75);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd75, smoothMPC_yy74, smoothMPC_beta75, smoothMPC_bmy75);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld75, smoothMPC_bmy75, smoothMPC_yy75);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld75, smoothMPC_Ysd76, smoothMPC_Lsd76);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd76, smoothMPC_Yd76);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd76, smoothMPC_Ld76);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd76, smoothMPC_yy75, smoothMPC_beta76, smoothMPC_bmy76);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld76, smoothMPC_bmy76, smoothMPC_yy76);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld76, smoothMPC_Ysd77, smoothMPC_Lsd77);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd77, smoothMPC_Yd77);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd77, smoothMPC_Ld77);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd77, smoothMPC_yy76, smoothMPC_beta77, smoothMPC_bmy77);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld77, smoothMPC_bmy77, smoothMPC_yy77);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld77, smoothMPC_Ysd78, smoothMPC_Lsd78);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd78, smoothMPC_Yd78);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd78, smoothMPC_Ld78);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd78, smoothMPC_yy77, smoothMPC_beta78, smoothMPC_bmy78);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld78, smoothMPC_bmy78, smoothMPC_yy78);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld78, smoothMPC_Ysd79, smoothMPC_Lsd79);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd79, smoothMPC_Yd79);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd79, smoothMPC_Ld79);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd79, smoothMPC_yy78, smoothMPC_beta79, smoothMPC_bmy79);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld79, smoothMPC_bmy79, smoothMPC_yy79);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld79, smoothMPC_Ysd80, smoothMPC_Lsd80);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd80, smoothMPC_Yd80);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd80, smoothMPC_Ld80);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd80, smoothMPC_yy79, smoothMPC_beta80, smoothMPC_bmy80);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld80, smoothMPC_bmy80, smoothMPC_yy80);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld80, smoothMPC_Ysd81, smoothMPC_Lsd81);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd81, smoothMPC_Yd81);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd81, smoothMPC_Ld81);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd81, smoothMPC_yy80, smoothMPC_beta81, smoothMPC_bmy81);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld81, smoothMPC_bmy81, smoothMPC_yy81);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld81, smoothMPC_Ysd82, smoothMPC_Lsd82);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd82, smoothMPC_Yd82);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd82, smoothMPC_Ld82);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd82, smoothMPC_yy81, smoothMPC_beta82, smoothMPC_bmy82);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld82, smoothMPC_bmy82, smoothMPC_yy82);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld82, smoothMPC_Ysd83, smoothMPC_Lsd83);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd83, smoothMPC_Yd83);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd83, smoothMPC_Ld83);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd83, smoothMPC_yy82, smoothMPC_beta83, smoothMPC_bmy83);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld83, smoothMPC_bmy83, smoothMPC_yy83);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld83, smoothMPC_Ysd84, smoothMPC_Lsd84);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd84, smoothMPC_Yd84);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd84, smoothMPC_Ld84);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd84, smoothMPC_yy83, smoothMPC_beta84, smoothMPC_bmy84);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld84, smoothMPC_bmy84, smoothMPC_yy84);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld84, smoothMPC_Ysd85, smoothMPC_Lsd85);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd85, smoothMPC_Yd85);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd85, smoothMPC_Ld85);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd85, smoothMPC_yy84, smoothMPC_beta85, smoothMPC_bmy85);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld85, smoothMPC_bmy85, smoothMPC_yy85);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld85, smoothMPC_Ysd86, smoothMPC_Lsd86);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd86, smoothMPC_Yd86);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd86, smoothMPC_Ld86);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd86, smoothMPC_yy85, smoothMPC_beta86, smoothMPC_bmy86);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld86, smoothMPC_bmy86, smoothMPC_yy86);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld86, smoothMPC_Ysd87, smoothMPC_Lsd87);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd87, smoothMPC_Yd87);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd87, smoothMPC_Ld87);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd87, smoothMPC_yy86, smoothMPC_beta87, smoothMPC_bmy87);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld87, smoothMPC_bmy87, smoothMPC_yy87);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld87, smoothMPC_Ysd88, smoothMPC_Lsd88);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd88, smoothMPC_Yd88);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd88, smoothMPC_Ld88);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd88, smoothMPC_yy87, smoothMPC_beta88, smoothMPC_bmy88);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld88, smoothMPC_bmy88, smoothMPC_yy88);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld88, smoothMPC_Ysd89, smoothMPC_Lsd89);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd89, smoothMPC_Yd89);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd89, smoothMPC_Ld89);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd89, smoothMPC_yy88, smoothMPC_beta89, smoothMPC_bmy89);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld89, smoothMPC_bmy89, smoothMPC_yy89);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld89, smoothMPC_Ysd90, smoothMPC_Lsd90);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd90, smoothMPC_Yd90);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd90, smoothMPC_Ld90);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd90, smoothMPC_yy89, smoothMPC_beta90, smoothMPC_bmy90);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld90, smoothMPC_bmy90, smoothMPC_yy90);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld90, smoothMPC_Ysd91, smoothMPC_Lsd91);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd91, smoothMPC_Yd91);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd91, smoothMPC_Ld91);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd91, smoothMPC_yy90, smoothMPC_beta91, smoothMPC_bmy91);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld91, smoothMPC_bmy91, smoothMPC_yy91);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld91, smoothMPC_Ysd92, smoothMPC_Lsd92);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd92, smoothMPC_Yd92);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd92, smoothMPC_Ld92);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd92, smoothMPC_yy91, smoothMPC_beta92, smoothMPC_bmy92);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld92, smoothMPC_bmy92, smoothMPC_yy92);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld92, smoothMPC_Ysd93, smoothMPC_Lsd93);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd93, smoothMPC_Yd93);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd93, smoothMPC_Ld93);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd93, smoothMPC_yy92, smoothMPC_beta93, smoothMPC_bmy93);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld93, smoothMPC_bmy93, smoothMPC_yy93);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld93, smoothMPC_Ysd94, smoothMPC_Lsd94);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd94, smoothMPC_Yd94);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd94, smoothMPC_Ld94);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd94, smoothMPC_yy93, smoothMPC_beta94, smoothMPC_bmy94);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld94, smoothMPC_bmy94, smoothMPC_yy94);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld94, smoothMPC_Ysd95, smoothMPC_Lsd95);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd95, smoothMPC_Yd95);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd95, smoothMPC_Ld95);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd95, smoothMPC_yy94, smoothMPC_beta95, smoothMPC_bmy95);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld95, smoothMPC_bmy95, smoothMPC_yy95);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld95, smoothMPC_Ysd96, smoothMPC_Lsd96);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd96, smoothMPC_Yd96);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd96, smoothMPC_Ld96);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd96, smoothMPC_yy95, smoothMPC_beta96, smoothMPC_bmy96);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld96, smoothMPC_bmy96, smoothMPC_yy96);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld96, smoothMPC_Ysd97, smoothMPC_Lsd97);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd97, smoothMPC_Yd97);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd97, smoothMPC_Ld97);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd97, smoothMPC_yy96, smoothMPC_beta97, smoothMPC_bmy97);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld97, smoothMPC_bmy97, smoothMPC_yy97);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld97, smoothMPC_Ysd98, smoothMPC_Lsd98);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd98, smoothMPC_Yd98);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd98, smoothMPC_Ld98);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd98, smoothMPC_yy97, smoothMPC_beta98, smoothMPC_bmy98);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld98, smoothMPC_bmy98, smoothMPC_yy98);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld98, smoothMPC_Ysd99, smoothMPC_Lsd99);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd99, smoothMPC_Yd99);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd99, smoothMPC_Ld99);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd99, smoothMPC_yy98, smoothMPC_beta99, smoothMPC_bmy99);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld99, smoothMPC_bmy99, smoothMPC_yy99);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld99, smoothMPC_yy99, smoothMPC_dvaff99);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd99, smoothMPC_dvaff99, smoothMPC_yy98, smoothMPC_bmy98);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld98, smoothMPC_bmy98, smoothMPC_dvaff98);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd98, smoothMPC_dvaff98, smoothMPC_yy97, smoothMPC_bmy97);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld97, smoothMPC_bmy97, smoothMPC_dvaff97);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd97, smoothMPC_dvaff97, smoothMPC_yy96, smoothMPC_bmy96);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld96, smoothMPC_bmy96, smoothMPC_dvaff96);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd96, smoothMPC_dvaff96, smoothMPC_yy95, smoothMPC_bmy95);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld95, smoothMPC_bmy95, smoothMPC_dvaff95);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd95, smoothMPC_dvaff95, smoothMPC_yy94, smoothMPC_bmy94);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld94, smoothMPC_bmy94, smoothMPC_dvaff94);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd94, smoothMPC_dvaff94, smoothMPC_yy93, smoothMPC_bmy93);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld93, smoothMPC_bmy93, smoothMPC_dvaff93);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd93, smoothMPC_dvaff93, smoothMPC_yy92, smoothMPC_bmy92);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld92, smoothMPC_bmy92, smoothMPC_dvaff92);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd92, smoothMPC_dvaff92, smoothMPC_yy91, smoothMPC_bmy91);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld91, smoothMPC_bmy91, smoothMPC_dvaff91);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd91, smoothMPC_dvaff91, smoothMPC_yy90, smoothMPC_bmy90);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld90, smoothMPC_bmy90, smoothMPC_dvaff90);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd90, smoothMPC_dvaff90, smoothMPC_yy89, smoothMPC_bmy89);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld89, smoothMPC_bmy89, smoothMPC_dvaff89);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd89, smoothMPC_dvaff89, smoothMPC_yy88, smoothMPC_bmy88);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld88, smoothMPC_bmy88, smoothMPC_dvaff88);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd88, smoothMPC_dvaff88, smoothMPC_yy87, smoothMPC_bmy87);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld87, smoothMPC_bmy87, smoothMPC_dvaff87);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd87, smoothMPC_dvaff87, smoothMPC_yy86, smoothMPC_bmy86);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld86, smoothMPC_bmy86, smoothMPC_dvaff86);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd86, smoothMPC_dvaff86, smoothMPC_yy85, smoothMPC_bmy85);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld85, smoothMPC_bmy85, smoothMPC_dvaff85);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd85, smoothMPC_dvaff85, smoothMPC_yy84, smoothMPC_bmy84);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld84, smoothMPC_bmy84, smoothMPC_dvaff84);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd84, smoothMPC_dvaff84, smoothMPC_yy83, smoothMPC_bmy83);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld83, smoothMPC_bmy83, smoothMPC_dvaff83);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd83, smoothMPC_dvaff83, smoothMPC_yy82, smoothMPC_bmy82);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld82, smoothMPC_bmy82, smoothMPC_dvaff82);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd82, smoothMPC_dvaff82, smoothMPC_yy81, smoothMPC_bmy81);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld81, smoothMPC_bmy81, smoothMPC_dvaff81);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd81, smoothMPC_dvaff81, smoothMPC_yy80, smoothMPC_bmy80);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld80, smoothMPC_bmy80, smoothMPC_dvaff80);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd80, smoothMPC_dvaff80, smoothMPC_yy79, smoothMPC_bmy79);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld79, smoothMPC_bmy79, smoothMPC_dvaff79);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd79, smoothMPC_dvaff79, smoothMPC_yy78, smoothMPC_bmy78);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld78, smoothMPC_bmy78, smoothMPC_dvaff78);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd78, smoothMPC_dvaff78, smoothMPC_yy77, smoothMPC_bmy77);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld77, smoothMPC_bmy77, smoothMPC_dvaff77);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd77, smoothMPC_dvaff77, smoothMPC_yy76, smoothMPC_bmy76);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld76, smoothMPC_bmy76, smoothMPC_dvaff76);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd76, smoothMPC_dvaff76, smoothMPC_yy75, smoothMPC_bmy75);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld75, smoothMPC_bmy75, smoothMPC_dvaff75);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd75, smoothMPC_dvaff75, smoothMPC_yy74, smoothMPC_bmy74);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld74, smoothMPC_bmy74, smoothMPC_dvaff74);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd74, smoothMPC_dvaff74, smoothMPC_yy73, smoothMPC_bmy73);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld73, smoothMPC_bmy73, smoothMPC_dvaff73);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd73, smoothMPC_dvaff73, smoothMPC_yy72, smoothMPC_bmy72);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld72, smoothMPC_bmy72, smoothMPC_dvaff72);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd72, smoothMPC_dvaff72, smoothMPC_yy71, smoothMPC_bmy71);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld71, smoothMPC_bmy71, smoothMPC_dvaff71);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd71, smoothMPC_dvaff71, smoothMPC_yy70, smoothMPC_bmy70);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld70, smoothMPC_bmy70, smoothMPC_dvaff70);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd70, smoothMPC_dvaff70, smoothMPC_yy69, smoothMPC_bmy69);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld69, smoothMPC_bmy69, smoothMPC_dvaff69);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd69, smoothMPC_dvaff69, smoothMPC_yy68, smoothMPC_bmy68);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld68, smoothMPC_bmy68, smoothMPC_dvaff68);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd68, smoothMPC_dvaff68, smoothMPC_yy67, smoothMPC_bmy67);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld67, smoothMPC_bmy67, smoothMPC_dvaff67);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd67, smoothMPC_dvaff67, smoothMPC_yy66, smoothMPC_bmy66);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld66, smoothMPC_bmy66, smoothMPC_dvaff66);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd66, smoothMPC_dvaff66, smoothMPC_yy65, smoothMPC_bmy65);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld65, smoothMPC_bmy65, smoothMPC_dvaff65);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd65, smoothMPC_dvaff65, smoothMPC_yy64, smoothMPC_bmy64);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld64, smoothMPC_bmy64, smoothMPC_dvaff64);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd64, smoothMPC_dvaff64, smoothMPC_yy63, smoothMPC_bmy63);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld63, smoothMPC_bmy63, smoothMPC_dvaff63);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd63, smoothMPC_dvaff63, smoothMPC_yy62, smoothMPC_bmy62);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld62, smoothMPC_bmy62, smoothMPC_dvaff62);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd62, smoothMPC_dvaff62, smoothMPC_yy61, smoothMPC_bmy61);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld61, smoothMPC_bmy61, smoothMPC_dvaff61);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd61, smoothMPC_dvaff61, smoothMPC_yy60, smoothMPC_bmy60);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld60, smoothMPC_bmy60, smoothMPC_dvaff60);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd60, smoothMPC_dvaff60, smoothMPC_yy59, smoothMPC_bmy59);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld59, smoothMPC_bmy59, smoothMPC_dvaff59);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd59, smoothMPC_dvaff59, smoothMPC_yy58, smoothMPC_bmy58);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld58, smoothMPC_bmy58, smoothMPC_dvaff58);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd58, smoothMPC_dvaff58, smoothMPC_yy57, smoothMPC_bmy57);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld57, smoothMPC_bmy57, smoothMPC_dvaff57);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd57, smoothMPC_dvaff57, smoothMPC_yy56, smoothMPC_bmy56);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld56, smoothMPC_bmy56, smoothMPC_dvaff56);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd56, smoothMPC_dvaff56, smoothMPC_yy55, smoothMPC_bmy55);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld55, smoothMPC_bmy55, smoothMPC_dvaff55);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd55, smoothMPC_dvaff55, smoothMPC_yy54, smoothMPC_bmy54);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld54, smoothMPC_bmy54, smoothMPC_dvaff54);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd54, smoothMPC_dvaff54, smoothMPC_yy53, smoothMPC_bmy53);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld53, smoothMPC_bmy53, smoothMPC_dvaff53);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd53, smoothMPC_dvaff53, smoothMPC_yy52, smoothMPC_bmy52);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld52, smoothMPC_bmy52, smoothMPC_dvaff52);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd52, smoothMPC_dvaff52, smoothMPC_yy51, smoothMPC_bmy51);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld51, smoothMPC_bmy51, smoothMPC_dvaff51);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd51, smoothMPC_dvaff51, smoothMPC_yy50, smoothMPC_bmy50);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld50, smoothMPC_bmy50, smoothMPC_dvaff50);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd50, smoothMPC_dvaff50, smoothMPC_yy49, smoothMPC_bmy49);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld49, smoothMPC_bmy49, smoothMPC_dvaff49);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd49, smoothMPC_dvaff49, smoothMPC_yy48, smoothMPC_bmy48);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld48, smoothMPC_bmy48, smoothMPC_dvaff48);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd48, smoothMPC_dvaff48, smoothMPC_yy47, smoothMPC_bmy47);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld47, smoothMPC_bmy47, smoothMPC_dvaff47);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd47, smoothMPC_dvaff47, smoothMPC_yy46, smoothMPC_bmy46);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld46, smoothMPC_bmy46, smoothMPC_dvaff46);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd46, smoothMPC_dvaff46, smoothMPC_yy45, smoothMPC_bmy45);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld45, smoothMPC_bmy45, smoothMPC_dvaff45);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd45, smoothMPC_dvaff45, smoothMPC_yy44, smoothMPC_bmy44);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld44, smoothMPC_bmy44, smoothMPC_dvaff44);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd44, smoothMPC_dvaff44, smoothMPC_yy43, smoothMPC_bmy43);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld43, smoothMPC_bmy43, smoothMPC_dvaff43);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd43, smoothMPC_dvaff43, smoothMPC_yy42, smoothMPC_bmy42);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld42, smoothMPC_bmy42, smoothMPC_dvaff42);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd42, smoothMPC_dvaff42, smoothMPC_yy41, smoothMPC_bmy41);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld41, smoothMPC_bmy41, smoothMPC_dvaff41);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd41, smoothMPC_dvaff41, smoothMPC_yy40, smoothMPC_bmy40);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld40, smoothMPC_bmy40, smoothMPC_dvaff40);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd40, smoothMPC_dvaff40, smoothMPC_yy39, smoothMPC_bmy39);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld39, smoothMPC_bmy39, smoothMPC_dvaff39);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd39, smoothMPC_dvaff39, smoothMPC_yy38, smoothMPC_bmy38);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld38, smoothMPC_bmy38, smoothMPC_dvaff38);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd38, smoothMPC_dvaff38, smoothMPC_yy37, smoothMPC_bmy37);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld37, smoothMPC_bmy37, smoothMPC_dvaff37);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd37, smoothMPC_dvaff37, smoothMPC_yy36, smoothMPC_bmy36);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld36, smoothMPC_bmy36, smoothMPC_dvaff36);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd36, smoothMPC_dvaff36, smoothMPC_yy35, smoothMPC_bmy35);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld35, smoothMPC_bmy35, smoothMPC_dvaff35);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd35, smoothMPC_dvaff35, smoothMPC_yy34, smoothMPC_bmy34);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld34, smoothMPC_bmy34, smoothMPC_dvaff34);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd34, smoothMPC_dvaff34, smoothMPC_yy33, smoothMPC_bmy33);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld33, smoothMPC_bmy33, smoothMPC_dvaff33);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd33, smoothMPC_dvaff33, smoothMPC_yy32, smoothMPC_bmy32);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld32, smoothMPC_bmy32, smoothMPC_dvaff32);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd32, smoothMPC_dvaff32, smoothMPC_yy31, smoothMPC_bmy31);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld31, smoothMPC_bmy31, smoothMPC_dvaff31);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd31, smoothMPC_dvaff31, smoothMPC_yy30, smoothMPC_bmy30);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld30, smoothMPC_bmy30, smoothMPC_dvaff30);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd30, smoothMPC_dvaff30, smoothMPC_yy29, smoothMPC_bmy29);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld29, smoothMPC_bmy29, smoothMPC_dvaff29);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd29, smoothMPC_dvaff29, smoothMPC_yy28, smoothMPC_bmy28);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld28, smoothMPC_bmy28, smoothMPC_dvaff28);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd28, smoothMPC_dvaff28, smoothMPC_yy27, smoothMPC_bmy27);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld27, smoothMPC_bmy27, smoothMPC_dvaff27);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd27, smoothMPC_dvaff27, smoothMPC_yy26, smoothMPC_bmy26);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld26, smoothMPC_bmy26, smoothMPC_dvaff26);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd26, smoothMPC_dvaff26, smoothMPC_yy25, smoothMPC_bmy25);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld25, smoothMPC_bmy25, smoothMPC_dvaff25);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd25, smoothMPC_dvaff25, smoothMPC_yy24, smoothMPC_bmy24);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld24, smoothMPC_bmy24, smoothMPC_dvaff24);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd24, smoothMPC_dvaff24, smoothMPC_yy23, smoothMPC_bmy23);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld23, smoothMPC_bmy23, smoothMPC_dvaff23);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd23, smoothMPC_dvaff23, smoothMPC_yy22, smoothMPC_bmy22);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld22, smoothMPC_bmy22, smoothMPC_dvaff22);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd22, smoothMPC_dvaff22, smoothMPC_yy21, smoothMPC_bmy21);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld21, smoothMPC_bmy21, smoothMPC_dvaff21);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd21, smoothMPC_dvaff21, smoothMPC_yy20, smoothMPC_bmy20);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld20, smoothMPC_bmy20, smoothMPC_dvaff20);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd20, smoothMPC_dvaff20, smoothMPC_yy19, smoothMPC_bmy19);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld19, smoothMPC_bmy19, smoothMPC_dvaff19);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd19, smoothMPC_dvaff19, smoothMPC_yy18, smoothMPC_bmy18);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld18, smoothMPC_bmy18, smoothMPC_dvaff18);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd18, smoothMPC_dvaff18, smoothMPC_yy17, smoothMPC_bmy17);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld17, smoothMPC_bmy17, smoothMPC_dvaff17);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd17, smoothMPC_dvaff17, smoothMPC_yy16, smoothMPC_bmy16);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld16, smoothMPC_bmy16, smoothMPC_dvaff16);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd16, smoothMPC_dvaff16, smoothMPC_yy15, smoothMPC_bmy15);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld15, smoothMPC_bmy15, smoothMPC_dvaff15);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd15, smoothMPC_dvaff15, smoothMPC_yy14, smoothMPC_bmy14);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld14, smoothMPC_bmy14, smoothMPC_dvaff14);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd14, smoothMPC_dvaff14, smoothMPC_yy13, smoothMPC_bmy13);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld13, smoothMPC_bmy13, smoothMPC_dvaff13);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd13, smoothMPC_dvaff13, smoothMPC_yy12, smoothMPC_bmy12);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld12, smoothMPC_bmy12, smoothMPC_dvaff12);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd12, smoothMPC_dvaff12, smoothMPC_yy11, smoothMPC_bmy11);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld11, smoothMPC_bmy11, smoothMPC_dvaff11);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd11, smoothMPC_dvaff11, smoothMPC_yy10, smoothMPC_bmy10);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld10, smoothMPC_bmy10, smoothMPC_dvaff10);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd10, smoothMPC_dvaff10, smoothMPC_yy09, smoothMPC_bmy09);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld09, smoothMPC_bmy09, smoothMPC_dvaff09);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd09, smoothMPC_dvaff09, smoothMPC_yy08, smoothMPC_bmy08);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld08, smoothMPC_bmy08, smoothMPC_dvaff08);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd08, smoothMPC_dvaff08, smoothMPC_yy07, smoothMPC_bmy07);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld07, smoothMPC_bmy07, smoothMPC_dvaff07);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd07, smoothMPC_dvaff07, smoothMPC_yy06, smoothMPC_bmy06);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld06, smoothMPC_bmy06, smoothMPC_dvaff06);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd06, smoothMPC_dvaff06, smoothMPC_yy05, smoothMPC_bmy05);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld05, smoothMPC_bmy05, smoothMPC_dvaff05);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd05, smoothMPC_dvaff05, smoothMPC_yy04, smoothMPC_bmy04);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld04, smoothMPC_bmy04, smoothMPC_dvaff04);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd04, smoothMPC_dvaff04, smoothMPC_yy03, smoothMPC_bmy03);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld03, smoothMPC_bmy03, smoothMPC_dvaff03);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd03, smoothMPC_dvaff03, smoothMPC_yy02, smoothMPC_bmy02);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld02, smoothMPC_bmy02, smoothMPC_dvaff02);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd02, smoothMPC_dvaff02, smoothMPC_yy01, smoothMPC_bmy01);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld01, smoothMPC_bmy01, smoothMPC_dvaff01);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd01, smoothMPC_dvaff01, smoothMPC_yy00, smoothMPC_bmy00);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld00, smoothMPC_bmy00, smoothMPC_dvaff00);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C1, smoothMPC_dvaff01, smoothMPC_D00, smoothMPC_dvaff00, smoothMPC_grad_eq00);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C2, smoothMPC_dvaff02, smoothMPC_D01, smoothMPC_dvaff01, smoothMPC_grad_eq01);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C3, smoothMPC_dvaff03, smoothMPC_D01, smoothMPC_dvaff02, smoothMPC_grad_eq02);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C4, smoothMPC_dvaff04, smoothMPC_D01, smoothMPC_dvaff03, smoothMPC_grad_eq03);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C5, smoothMPC_dvaff05, smoothMPC_D01, smoothMPC_dvaff04, smoothMPC_grad_eq04);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C6, smoothMPC_dvaff06, smoothMPC_D01, smoothMPC_dvaff05, smoothMPC_grad_eq05);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C7, smoothMPC_dvaff07, smoothMPC_D01, smoothMPC_dvaff06, smoothMPC_grad_eq06);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C8, smoothMPC_dvaff08, smoothMPC_D01, smoothMPC_dvaff07, smoothMPC_grad_eq07);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C9, smoothMPC_dvaff09, smoothMPC_D01, smoothMPC_dvaff08, smoothMPC_grad_eq08);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C10, smoothMPC_dvaff10, smoothMPC_D01, smoothMPC_dvaff09, smoothMPC_grad_eq09);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C11, smoothMPC_dvaff11, smoothMPC_D01, smoothMPC_dvaff10, smoothMPC_grad_eq10);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C12, smoothMPC_dvaff12, smoothMPC_D01, smoothMPC_dvaff11, smoothMPC_grad_eq11);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C13, smoothMPC_dvaff13, smoothMPC_D01, smoothMPC_dvaff12, smoothMPC_grad_eq12);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C14, smoothMPC_dvaff14, smoothMPC_D01, smoothMPC_dvaff13, smoothMPC_grad_eq13);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C15, smoothMPC_dvaff15, smoothMPC_D01, smoothMPC_dvaff14, smoothMPC_grad_eq14);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C16, smoothMPC_dvaff16, smoothMPC_D01, smoothMPC_dvaff15, smoothMPC_grad_eq15);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C17, smoothMPC_dvaff17, smoothMPC_D01, smoothMPC_dvaff16, smoothMPC_grad_eq16);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C18, smoothMPC_dvaff18, smoothMPC_D01, smoothMPC_dvaff17, smoothMPC_grad_eq17);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C19, smoothMPC_dvaff19, smoothMPC_D01, smoothMPC_dvaff18, smoothMPC_grad_eq18);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C20, smoothMPC_dvaff20, smoothMPC_D01, smoothMPC_dvaff19, smoothMPC_grad_eq19);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C21, smoothMPC_dvaff21, smoothMPC_D01, smoothMPC_dvaff20, smoothMPC_grad_eq20);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C22, smoothMPC_dvaff22, smoothMPC_D01, smoothMPC_dvaff21, smoothMPC_grad_eq21);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C23, smoothMPC_dvaff23, smoothMPC_D01, smoothMPC_dvaff22, smoothMPC_grad_eq22);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C24, smoothMPC_dvaff24, smoothMPC_D01, smoothMPC_dvaff23, smoothMPC_grad_eq23);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C25, smoothMPC_dvaff25, smoothMPC_D01, smoothMPC_dvaff24, smoothMPC_grad_eq24);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C26, smoothMPC_dvaff26, smoothMPC_D01, smoothMPC_dvaff25, smoothMPC_grad_eq25);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C27, smoothMPC_dvaff27, smoothMPC_D01, smoothMPC_dvaff26, smoothMPC_grad_eq26);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C28, smoothMPC_dvaff28, smoothMPC_D01, smoothMPC_dvaff27, smoothMPC_grad_eq27);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C29, smoothMPC_dvaff29, smoothMPC_D01, smoothMPC_dvaff28, smoothMPC_grad_eq28);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C30, smoothMPC_dvaff30, smoothMPC_D01, smoothMPC_dvaff29, smoothMPC_grad_eq29);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C31, smoothMPC_dvaff31, smoothMPC_D01, smoothMPC_dvaff30, smoothMPC_grad_eq30);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C32, smoothMPC_dvaff32, smoothMPC_D01, smoothMPC_dvaff31, smoothMPC_grad_eq31);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C33, smoothMPC_dvaff33, smoothMPC_D01, smoothMPC_dvaff32, smoothMPC_grad_eq32);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C34, smoothMPC_dvaff34, smoothMPC_D01, smoothMPC_dvaff33, smoothMPC_grad_eq33);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C35, smoothMPC_dvaff35, smoothMPC_D01, smoothMPC_dvaff34, smoothMPC_grad_eq34);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C36, smoothMPC_dvaff36, smoothMPC_D01, smoothMPC_dvaff35, smoothMPC_grad_eq35);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C37, smoothMPC_dvaff37, smoothMPC_D01, smoothMPC_dvaff36, smoothMPC_grad_eq36);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C38, smoothMPC_dvaff38, smoothMPC_D01, smoothMPC_dvaff37, smoothMPC_grad_eq37);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C39, smoothMPC_dvaff39, smoothMPC_D01, smoothMPC_dvaff38, smoothMPC_grad_eq38);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C40, smoothMPC_dvaff40, smoothMPC_D01, smoothMPC_dvaff39, smoothMPC_grad_eq39);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C41, smoothMPC_dvaff41, smoothMPC_D01, smoothMPC_dvaff40, smoothMPC_grad_eq40);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C42, smoothMPC_dvaff42, smoothMPC_D01, smoothMPC_dvaff41, smoothMPC_grad_eq41);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C43, smoothMPC_dvaff43, smoothMPC_D01, smoothMPC_dvaff42, smoothMPC_grad_eq42);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C44, smoothMPC_dvaff44, smoothMPC_D01, smoothMPC_dvaff43, smoothMPC_grad_eq43);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C45, smoothMPC_dvaff45, smoothMPC_D01, smoothMPC_dvaff44, smoothMPC_grad_eq44);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C46, smoothMPC_dvaff46, smoothMPC_D01, smoothMPC_dvaff45, smoothMPC_grad_eq45);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C47, smoothMPC_dvaff47, smoothMPC_D01, smoothMPC_dvaff46, smoothMPC_grad_eq46);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C48, smoothMPC_dvaff48, smoothMPC_D01, smoothMPC_dvaff47, smoothMPC_grad_eq47);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C49, smoothMPC_dvaff49, smoothMPC_D01, smoothMPC_dvaff48, smoothMPC_grad_eq48);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C50, smoothMPC_dvaff50, smoothMPC_D01, smoothMPC_dvaff49, smoothMPC_grad_eq49);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C51, smoothMPC_dvaff51, smoothMPC_D01, smoothMPC_dvaff50, smoothMPC_grad_eq50);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C52, smoothMPC_dvaff52, smoothMPC_D01, smoothMPC_dvaff51, smoothMPC_grad_eq51);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C53, smoothMPC_dvaff53, smoothMPC_D01, smoothMPC_dvaff52, smoothMPC_grad_eq52);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C54, smoothMPC_dvaff54, smoothMPC_D01, smoothMPC_dvaff53, smoothMPC_grad_eq53);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C55, smoothMPC_dvaff55, smoothMPC_D01, smoothMPC_dvaff54, smoothMPC_grad_eq54);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C56, smoothMPC_dvaff56, smoothMPC_D01, smoothMPC_dvaff55, smoothMPC_grad_eq55);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C57, smoothMPC_dvaff57, smoothMPC_D01, smoothMPC_dvaff56, smoothMPC_grad_eq56);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C58, smoothMPC_dvaff58, smoothMPC_D01, smoothMPC_dvaff57, smoothMPC_grad_eq57);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C59, smoothMPC_dvaff59, smoothMPC_D01, smoothMPC_dvaff58, smoothMPC_grad_eq58);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C60, smoothMPC_dvaff60, smoothMPC_D01, smoothMPC_dvaff59, smoothMPC_grad_eq59);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C61, smoothMPC_dvaff61, smoothMPC_D01, smoothMPC_dvaff60, smoothMPC_grad_eq60);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C62, smoothMPC_dvaff62, smoothMPC_D01, smoothMPC_dvaff61, smoothMPC_grad_eq61);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C63, smoothMPC_dvaff63, smoothMPC_D01, smoothMPC_dvaff62, smoothMPC_grad_eq62);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C64, smoothMPC_dvaff64, smoothMPC_D01, smoothMPC_dvaff63, smoothMPC_grad_eq63);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C65, smoothMPC_dvaff65, smoothMPC_D01, smoothMPC_dvaff64, smoothMPC_grad_eq64);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C66, smoothMPC_dvaff66, smoothMPC_D01, smoothMPC_dvaff65, smoothMPC_grad_eq65);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C67, smoothMPC_dvaff67, smoothMPC_D01, smoothMPC_dvaff66, smoothMPC_grad_eq66);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C68, smoothMPC_dvaff68, smoothMPC_D01, smoothMPC_dvaff67, smoothMPC_grad_eq67);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C69, smoothMPC_dvaff69, smoothMPC_D01, smoothMPC_dvaff68, smoothMPC_grad_eq68);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C70, smoothMPC_dvaff70, smoothMPC_D01, smoothMPC_dvaff69, smoothMPC_grad_eq69);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C71, smoothMPC_dvaff71, smoothMPC_D01, smoothMPC_dvaff70, smoothMPC_grad_eq70);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C72, smoothMPC_dvaff72, smoothMPC_D01, smoothMPC_dvaff71, smoothMPC_grad_eq71);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C73, smoothMPC_dvaff73, smoothMPC_D01, smoothMPC_dvaff72, smoothMPC_grad_eq72);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C74, smoothMPC_dvaff74, smoothMPC_D01, smoothMPC_dvaff73, smoothMPC_grad_eq73);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C75, smoothMPC_dvaff75, smoothMPC_D01, smoothMPC_dvaff74, smoothMPC_grad_eq74);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C76, smoothMPC_dvaff76, smoothMPC_D01, smoothMPC_dvaff75, smoothMPC_grad_eq75);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C77, smoothMPC_dvaff77, smoothMPC_D01, smoothMPC_dvaff76, smoothMPC_grad_eq76);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C78, smoothMPC_dvaff78, smoothMPC_D01, smoothMPC_dvaff77, smoothMPC_grad_eq77);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C79, smoothMPC_dvaff79, smoothMPC_D01, smoothMPC_dvaff78, smoothMPC_grad_eq78);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C80, smoothMPC_dvaff80, smoothMPC_D01, smoothMPC_dvaff79, smoothMPC_grad_eq79);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C81, smoothMPC_dvaff81, smoothMPC_D01, smoothMPC_dvaff80, smoothMPC_grad_eq80);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C82, smoothMPC_dvaff82, smoothMPC_D01, smoothMPC_dvaff81, smoothMPC_grad_eq81);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C83, smoothMPC_dvaff83, smoothMPC_D01, smoothMPC_dvaff82, smoothMPC_grad_eq82);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C84, smoothMPC_dvaff84, smoothMPC_D01, smoothMPC_dvaff83, smoothMPC_grad_eq83);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C85, smoothMPC_dvaff85, smoothMPC_D01, smoothMPC_dvaff84, smoothMPC_grad_eq84);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C86, smoothMPC_dvaff86, smoothMPC_D01, smoothMPC_dvaff85, smoothMPC_grad_eq85);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C87, smoothMPC_dvaff87, smoothMPC_D01, smoothMPC_dvaff86, smoothMPC_grad_eq86);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C88, smoothMPC_dvaff88, smoothMPC_D01, smoothMPC_dvaff87, smoothMPC_grad_eq87);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C89, smoothMPC_dvaff89, smoothMPC_D01, smoothMPC_dvaff88, smoothMPC_grad_eq88);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C90, smoothMPC_dvaff90, smoothMPC_D01, smoothMPC_dvaff89, smoothMPC_grad_eq89);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C91, smoothMPC_dvaff91, smoothMPC_D01, smoothMPC_dvaff90, smoothMPC_grad_eq90);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C92, smoothMPC_dvaff92, smoothMPC_D01, smoothMPC_dvaff91, smoothMPC_grad_eq91);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C93, smoothMPC_dvaff93, smoothMPC_D01, smoothMPC_dvaff92, smoothMPC_grad_eq92);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C94, smoothMPC_dvaff94, smoothMPC_D01, smoothMPC_dvaff93, smoothMPC_grad_eq93);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C95, smoothMPC_dvaff95, smoothMPC_D01, smoothMPC_dvaff94, smoothMPC_grad_eq94);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C96, smoothMPC_dvaff96, smoothMPC_D01, smoothMPC_dvaff95, smoothMPC_grad_eq95);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C97, smoothMPC_dvaff97, smoothMPC_D01, smoothMPC_dvaff96, smoothMPC_grad_eq96);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C98, smoothMPC_dvaff98, smoothMPC_D01, smoothMPC_dvaff97, smoothMPC_grad_eq97);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C99, smoothMPC_dvaff99, smoothMPC_D01, smoothMPC_dvaff98, smoothMPC_grad_eq98);
smoothMPC_LA_DIAGZERO_MTVM_3_3(smoothMPC_D99, smoothMPC_dvaff99, smoothMPC_grad_eq99);
smoothMPC_LA_VSUB2_1092(smoothMPC_rd, smoothMPC_grad_eq, smoothMPC_rd);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi00, smoothMPC_rd00, smoothMPC_dzaff00);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi01, smoothMPC_rd01, smoothMPC_dzaff01);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi02, smoothMPC_rd02, smoothMPC_dzaff02);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi03, smoothMPC_rd03, smoothMPC_dzaff03);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi04, smoothMPC_rd04, smoothMPC_dzaff04);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi05, smoothMPC_rd05, smoothMPC_dzaff05);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi06, smoothMPC_rd06, smoothMPC_dzaff06);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi07, smoothMPC_rd07, smoothMPC_dzaff07);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi08, smoothMPC_rd08, smoothMPC_dzaff08);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi09, smoothMPC_rd09, smoothMPC_dzaff09);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi10, smoothMPC_rd10, smoothMPC_dzaff10);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi11, smoothMPC_rd11, smoothMPC_dzaff11);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi12, smoothMPC_rd12, smoothMPC_dzaff12);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi13, smoothMPC_rd13, smoothMPC_dzaff13);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi14, smoothMPC_rd14, smoothMPC_dzaff14);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi15, smoothMPC_rd15, smoothMPC_dzaff15);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi16, smoothMPC_rd16, smoothMPC_dzaff16);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi17, smoothMPC_rd17, smoothMPC_dzaff17);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi18, smoothMPC_rd18, smoothMPC_dzaff18);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi19, smoothMPC_rd19, smoothMPC_dzaff19);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi20, smoothMPC_rd20, smoothMPC_dzaff20);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi21, smoothMPC_rd21, smoothMPC_dzaff21);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi22, smoothMPC_rd22, smoothMPC_dzaff22);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi23, smoothMPC_rd23, smoothMPC_dzaff23);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi24, smoothMPC_rd24, smoothMPC_dzaff24);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi25, smoothMPC_rd25, smoothMPC_dzaff25);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi26, smoothMPC_rd26, smoothMPC_dzaff26);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi27, smoothMPC_rd27, smoothMPC_dzaff27);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi28, smoothMPC_rd28, smoothMPC_dzaff28);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi29, smoothMPC_rd29, smoothMPC_dzaff29);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi30, smoothMPC_rd30, smoothMPC_dzaff30);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi31, smoothMPC_rd31, smoothMPC_dzaff31);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi32, smoothMPC_rd32, smoothMPC_dzaff32);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi33, smoothMPC_rd33, smoothMPC_dzaff33);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi34, smoothMPC_rd34, smoothMPC_dzaff34);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi35, smoothMPC_rd35, smoothMPC_dzaff35);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi36, smoothMPC_rd36, smoothMPC_dzaff36);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi37, smoothMPC_rd37, smoothMPC_dzaff37);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi38, smoothMPC_rd38, smoothMPC_dzaff38);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi39, smoothMPC_rd39, smoothMPC_dzaff39);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi40, smoothMPC_rd40, smoothMPC_dzaff40);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi41, smoothMPC_rd41, smoothMPC_dzaff41);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi42, smoothMPC_rd42, smoothMPC_dzaff42);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi43, smoothMPC_rd43, smoothMPC_dzaff43);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi44, smoothMPC_rd44, smoothMPC_dzaff44);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi45, smoothMPC_rd45, smoothMPC_dzaff45);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi46, smoothMPC_rd46, smoothMPC_dzaff46);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi47, smoothMPC_rd47, smoothMPC_dzaff47);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi48, smoothMPC_rd48, smoothMPC_dzaff48);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi49, smoothMPC_rd49, smoothMPC_dzaff49);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi50, smoothMPC_rd50, smoothMPC_dzaff50);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi51, smoothMPC_rd51, smoothMPC_dzaff51);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi52, smoothMPC_rd52, smoothMPC_dzaff52);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi53, smoothMPC_rd53, smoothMPC_dzaff53);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi54, smoothMPC_rd54, smoothMPC_dzaff54);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi55, smoothMPC_rd55, smoothMPC_dzaff55);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi56, smoothMPC_rd56, smoothMPC_dzaff56);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi57, smoothMPC_rd57, smoothMPC_dzaff57);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi58, smoothMPC_rd58, smoothMPC_dzaff58);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi59, smoothMPC_rd59, smoothMPC_dzaff59);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi60, smoothMPC_rd60, smoothMPC_dzaff60);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi61, smoothMPC_rd61, smoothMPC_dzaff61);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi62, smoothMPC_rd62, smoothMPC_dzaff62);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi63, smoothMPC_rd63, smoothMPC_dzaff63);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi64, smoothMPC_rd64, smoothMPC_dzaff64);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi65, smoothMPC_rd65, smoothMPC_dzaff65);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi66, smoothMPC_rd66, smoothMPC_dzaff66);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi67, smoothMPC_rd67, smoothMPC_dzaff67);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi68, smoothMPC_rd68, smoothMPC_dzaff68);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi69, smoothMPC_rd69, smoothMPC_dzaff69);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi70, smoothMPC_rd70, smoothMPC_dzaff70);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi71, smoothMPC_rd71, smoothMPC_dzaff71);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi72, smoothMPC_rd72, smoothMPC_dzaff72);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi73, smoothMPC_rd73, smoothMPC_dzaff73);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi74, smoothMPC_rd74, smoothMPC_dzaff74);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi75, smoothMPC_rd75, smoothMPC_dzaff75);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi76, smoothMPC_rd76, smoothMPC_dzaff76);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi77, smoothMPC_rd77, smoothMPC_dzaff77);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi78, smoothMPC_rd78, smoothMPC_dzaff78);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi79, smoothMPC_rd79, smoothMPC_dzaff79);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi80, smoothMPC_rd80, smoothMPC_dzaff80);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi81, smoothMPC_rd81, smoothMPC_dzaff81);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi82, smoothMPC_rd82, smoothMPC_dzaff82);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi83, smoothMPC_rd83, smoothMPC_dzaff83);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi84, smoothMPC_rd84, smoothMPC_dzaff84);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi85, smoothMPC_rd85, smoothMPC_dzaff85);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi86, smoothMPC_rd86, smoothMPC_dzaff86);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi87, smoothMPC_rd87, smoothMPC_dzaff87);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi88, smoothMPC_rd88, smoothMPC_dzaff88);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi89, smoothMPC_rd89, smoothMPC_dzaff89);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi90, smoothMPC_rd90, smoothMPC_dzaff90);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi91, smoothMPC_rd91, smoothMPC_dzaff91);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi92, smoothMPC_rd92, smoothMPC_dzaff92);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi93, smoothMPC_rd93, smoothMPC_dzaff93);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi94, smoothMPC_rd94, smoothMPC_dzaff94);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi95, smoothMPC_rd95, smoothMPC_dzaff95);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi96, smoothMPC_rd96, smoothMPC_dzaff96);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi97, smoothMPC_rd97, smoothMPC_dzaff97);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi98, smoothMPC_rd98, smoothMPC_dzaff98);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_3(smoothMPC_Phi99, smoothMPC_rd99, smoothMPC_dzaff99);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff00, smoothMPC_lbIdx00, smoothMPC_rilb00, smoothMPC_dslbaff00);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb00, smoothMPC_dslbaff00, smoothMPC_llb00, smoothMPC_dllbaff00);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub00, smoothMPC_dzaff00, smoothMPC_ubIdx00, smoothMPC_dsubaff00);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub00, smoothMPC_dsubaff00, smoothMPC_lub00, smoothMPC_dlubaff00);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff01, smoothMPC_lbIdx01, smoothMPC_rilb01, smoothMPC_dslbaff01);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb01, smoothMPC_dslbaff01, smoothMPC_llb01, smoothMPC_dllbaff01);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub01, smoothMPC_dzaff01, smoothMPC_ubIdx01, smoothMPC_dsubaff01);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub01, smoothMPC_dsubaff01, smoothMPC_lub01, smoothMPC_dlubaff01);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff02, smoothMPC_lbIdx02, smoothMPC_rilb02, smoothMPC_dslbaff02);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb02, smoothMPC_dslbaff02, smoothMPC_llb02, smoothMPC_dllbaff02);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub02, smoothMPC_dzaff02, smoothMPC_ubIdx02, smoothMPC_dsubaff02);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub02, smoothMPC_dsubaff02, smoothMPC_lub02, smoothMPC_dlubaff02);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff03, smoothMPC_lbIdx03, smoothMPC_rilb03, smoothMPC_dslbaff03);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb03, smoothMPC_dslbaff03, smoothMPC_llb03, smoothMPC_dllbaff03);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub03, smoothMPC_dzaff03, smoothMPC_ubIdx03, smoothMPC_dsubaff03);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub03, smoothMPC_dsubaff03, smoothMPC_lub03, smoothMPC_dlubaff03);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff04, smoothMPC_lbIdx04, smoothMPC_rilb04, smoothMPC_dslbaff04);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb04, smoothMPC_dslbaff04, smoothMPC_llb04, smoothMPC_dllbaff04);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub04, smoothMPC_dzaff04, smoothMPC_ubIdx04, smoothMPC_dsubaff04);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub04, smoothMPC_dsubaff04, smoothMPC_lub04, smoothMPC_dlubaff04);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff05, smoothMPC_lbIdx05, smoothMPC_rilb05, smoothMPC_dslbaff05);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb05, smoothMPC_dslbaff05, smoothMPC_llb05, smoothMPC_dllbaff05);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub05, smoothMPC_dzaff05, smoothMPC_ubIdx05, smoothMPC_dsubaff05);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub05, smoothMPC_dsubaff05, smoothMPC_lub05, smoothMPC_dlubaff05);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff06, smoothMPC_lbIdx06, smoothMPC_rilb06, smoothMPC_dslbaff06);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb06, smoothMPC_dslbaff06, smoothMPC_llb06, smoothMPC_dllbaff06);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub06, smoothMPC_dzaff06, smoothMPC_ubIdx06, smoothMPC_dsubaff06);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub06, smoothMPC_dsubaff06, smoothMPC_lub06, smoothMPC_dlubaff06);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff07, smoothMPC_lbIdx07, smoothMPC_rilb07, smoothMPC_dslbaff07);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb07, smoothMPC_dslbaff07, smoothMPC_llb07, smoothMPC_dllbaff07);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub07, smoothMPC_dzaff07, smoothMPC_ubIdx07, smoothMPC_dsubaff07);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub07, smoothMPC_dsubaff07, smoothMPC_lub07, smoothMPC_dlubaff07);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff08, smoothMPC_lbIdx08, smoothMPC_rilb08, smoothMPC_dslbaff08);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb08, smoothMPC_dslbaff08, smoothMPC_llb08, smoothMPC_dllbaff08);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub08, smoothMPC_dzaff08, smoothMPC_ubIdx08, smoothMPC_dsubaff08);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub08, smoothMPC_dsubaff08, smoothMPC_lub08, smoothMPC_dlubaff08);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff09, smoothMPC_lbIdx09, smoothMPC_rilb09, smoothMPC_dslbaff09);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb09, smoothMPC_dslbaff09, smoothMPC_llb09, smoothMPC_dllbaff09);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub09, smoothMPC_dzaff09, smoothMPC_ubIdx09, smoothMPC_dsubaff09);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub09, smoothMPC_dsubaff09, smoothMPC_lub09, smoothMPC_dlubaff09);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff10, smoothMPC_lbIdx10, smoothMPC_rilb10, smoothMPC_dslbaff10);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb10, smoothMPC_dslbaff10, smoothMPC_llb10, smoothMPC_dllbaff10);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub10, smoothMPC_dzaff10, smoothMPC_ubIdx10, smoothMPC_dsubaff10);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub10, smoothMPC_dsubaff10, smoothMPC_lub10, smoothMPC_dlubaff10);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff11, smoothMPC_lbIdx11, smoothMPC_rilb11, smoothMPC_dslbaff11);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb11, smoothMPC_dslbaff11, smoothMPC_llb11, smoothMPC_dllbaff11);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub11, smoothMPC_dzaff11, smoothMPC_ubIdx11, smoothMPC_dsubaff11);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub11, smoothMPC_dsubaff11, smoothMPC_lub11, smoothMPC_dlubaff11);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff12, smoothMPC_lbIdx12, smoothMPC_rilb12, smoothMPC_dslbaff12);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb12, smoothMPC_dslbaff12, smoothMPC_llb12, smoothMPC_dllbaff12);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub12, smoothMPC_dzaff12, smoothMPC_ubIdx12, smoothMPC_dsubaff12);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub12, smoothMPC_dsubaff12, smoothMPC_lub12, smoothMPC_dlubaff12);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff13, smoothMPC_lbIdx13, smoothMPC_rilb13, smoothMPC_dslbaff13);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb13, smoothMPC_dslbaff13, smoothMPC_llb13, smoothMPC_dllbaff13);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub13, smoothMPC_dzaff13, smoothMPC_ubIdx13, smoothMPC_dsubaff13);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub13, smoothMPC_dsubaff13, smoothMPC_lub13, smoothMPC_dlubaff13);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff14, smoothMPC_lbIdx14, smoothMPC_rilb14, smoothMPC_dslbaff14);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb14, smoothMPC_dslbaff14, smoothMPC_llb14, smoothMPC_dllbaff14);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub14, smoothMPC_dzaff14, smoothMPC_ubIdx14, smoothMPC_dsubaff14);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub14, smoothMPC_dsubaff14, smoothMPC_lub14, smoothMPC_dlubaff14);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff15, smoothMPC_lbIdx15, smoothMPC_rilb15, smoothMPC_dslbaff15);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb15, smoothMPC_dslbaff15, smoothMPC_llb15, smoothMPC_dllbaff15);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub15, smoothMPC_dzaff15, smoothMPC_ubIdx15, smoothMPC_dsubaff15);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub15, smoothMPC_dsubaff15, smoothMPC_lub15, smoothMPC_dlubaff15);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff16, smoothMPC_lbIdx16, smoothMPC_rilb16, smoothMPC_dslbaff16);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb16, smoothMPC_dslbaff16, smoothMPC_llb16, smoothMPC_dllbaff16);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub16, smoothMPC_dzaff16, smoothMPC_ubIdx16, smoothMPC_dsubaff16);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub16, smoothMPC_dsubaff16, smoothMPC_lub16, smoothMPC_dlubaff16);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff17, smoothMPC_lbIdx17, smoothMPC_rilb17, smoothMPC_dslbaff17);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb17, smoothMPC_dslbaff17, smoothMPC_llb17, smoothMPC_dllbaff17);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub17, smoothMPC_dzaff17, smoothMPC_ubIdx17, smoothMPC_dsubaff17);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub17, smoothMPC_dsubaff17, smoothMPC_lub17, smoothMPC_dlubaff17);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff18, smoothMPC_lbIdx18, smoothMPC_rilb18, smoothMPC_dslbaff18);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb18, smoothMPC_dslbaff18, smoothMPC_llb18, smoothMPC_dllbaff18);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub18, smoothMPC_dzaff18, smoothMPC_ubIdx18, smoothMPC_dsubaff18);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub18, smoothMPC_dsubaff18, smoothMPC_lub18, smoothMPC_dlubaff18);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff19, smoothMPC_lbIdx19, smoothMPC_rilb19, smoothMPC_dslbaff19);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb19, smoothMPC_dslbaff19, smoothMPC_llb19, smoothMPC_dllbaff19);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub19, smoothMPC_dzaff19, smoothMPC_ubIdx19, smoothMPC_dsubaff19);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub19, smoothMPC_dsubaff19, smoothMPC_lub19, smoothMPC_dlubaff19);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff20, smoothMPC_lbIdx20, smoothMPC_rilb20, smoothMPC_dslbaff20);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb20, smoothMPC_dslbaff20, smoothMPC_llb20, smoothMPC_dllbaff20);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub20, smoothMPC_dzaff20, smoothMPC_ubIdx20, smoothMPC_dsubaff20);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub20, smoothMPC_dsubaff20, smoothMPC_lub20, smoothMPC_dlubaff20);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff21, smoothMPC_lbIdx21, smoothMPC_rilb21, smoothMPC_dslbaff21);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb21, smoothMPC_dslbaff21, smoothMPC_llb21, smoothMPC_dllbaff21);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub21, smoothMPC_dzaff21, smoothMPC_ubIdx21, smoothMPC_dsubaff21);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub21, smoothMPC_dsubaff21, smoothMPC_lub21, smoothMPC_dlubaff21);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff22, smoothMPC_lbIdx22, smoothMPC_rilb22, smoothMPC_dslbaff22);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb22, smoothMPC_dslbaff22, smoothMPC_llb22, smoothMPC_dllbaff22);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub22, smoothMPC_dzaff22, smoothMPC_ubIdx22, smoothMPC_dsubaff22);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub22, smoothMPC_dsubaff22, smoothMPC_lub22, smoothMPC_dlubaff22);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff23, smoothMPC_lbIdx23, smoothMPC_rilb23, smoothMPC_dslbaff23);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb23, smoothMPC_dslbaff23, smoothMPC_llb23, smoothMPC_dllbaff23);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub23, smoothMPC_dzaff23, smoothMPC_ubIdx23, smoothMPC_dsubaff23);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub23, smoothMPC_dsubaff23, smoothMPC_lub23, smoothMPC_dlubaff23);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff24, smoothMPC_lbIdx24, smoothMPC_rilb24, smoothMPC_dslbaff24);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb24, smoothMPC_dslbaff24, smoothMPC_llb24, smoothMPC_dllbaff24);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub24, smoothMPC_dzaff24, smoothMPC_ubIdx24, smoothMPC_dsubaff24);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub24, smoothMPC_dsubaff24, smoothMPC_lub24, smoothMPC_dlubaff24);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff25, smoothMPC_lbIdx25, smoothMPC_rilb25, smoothMPC_dslbaff25);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb25, smoothMPC_dslbaff25, smoothMPC_llb25, smoothMPC_dllbaff25);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub25, smoothMPC_dzaff25, smoothMPC_ubIdx25, smoothMPC_dsubaff25);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub25, smoothMPC_dsubaff25, smoothMPC_lub25, smoothMPC_dlubaff25);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff26, smoothMPC_lbIdx26, smoothMPC_rilb26, smoothMPC_dslbaff26);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb26, smoothMPC_dslbaff26, smoothMPC_llb26, smoothMPC_dllbaff26);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub26, smoothMPC_dzaff26, smoothMPC_ubIdx26, smoothMPC_dsubaff26);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub26, smoothMPC_dsubaff26, smoothMPC_lub26, smoothMPC_dlubaff26);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff27, smoothMPC_lbIdx27, smoothMPC_rilb27, smoothMPC_dslbaff27);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb27, smoothMPC_dslbaff27, smoothMPC_llb27, smoothMPC_dllbaff27);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub27, smoothMPC_dzaff27, smoothMPC_ubIdx27, smoothMPC_dsubaff27);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub27, smoothMPC_dsubaff27, smoothMPC_lub27, smoothMPC_dlubaff27);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff28, smoothMPC_lbIdx28, smoothMPC_rilb28, smoothMPC_dslbaff28);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb28, smoothMPC_dslbaff28, smoothMPC_llb28, smoothMPC_dllbaff28);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub28, smoothMPC_dzaff28, smoothMPC_ubIdx28, smoothMPC_dsubaff28);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub28, smoothMPC_dsubaff28, smoothMPC_lub28, smoothMPC_dlubaff28);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff29, smoothMPC_lbIdx29, smoothMPC_rilb29, smoothMPC_dslbaff29);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb29, smoothMPC_dslbaff29, smoothMPC_llb29, smoothMPC_dllbaff29);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub29, smoothMPC_dzaff29, smoothMPC_ubIdx29, smoothMPC_dsubaff29);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub29, smoothMPC_dsubaff29, smoothMPC_lub29, smoothMPC_dlubaff29);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff30, smoothMPC_lbIdx30, smoothMPC_rilb30, smoothMPC_dslbaff30);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb30, smoothMPC_dslbaff30, smoothMPC_llb30, smoothMPC_dllbaff30);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub30, smoothMPC_dzaff30, smoothMPC_ubIdx30, smoothMPC_dsubaff30);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub30, smoothMPC_dsubaff30, smoothMPC_lub30, smoothMPC_dlubaff30);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff31, smoothMPC_lbIdx31, smoothMPC_rilb31, smoothMPC_dslbaff31);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb31, smoothMPC_dslbaff31, smoothMPC_llb31, smoothMPC_dllbaff31);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub31, smoothMPC_dzaff31, smoothMPC_ubIdx31, smoothMPC_dsubaff31);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub31, smoothMPC_dsubaff31, smoothMPC_lub31, smoothMPC_dlubaff31);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff32, smoothMPC_lbIdx32, smoothMPC_rilb32, smoothMPC_dslbaff32);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb32, smoothMPC_dslbaff32, smoothMPC_llb32, smoothMPC_dllbaff32);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub32, smoothMPC_dzaff32, smoothMPC_ubIdx32, smoothMPC_dsubaff32);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub32, smoothMPC_dsubaff32, smoothMPC_lub32, smoothMPC_dlubaff32);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff33, smoothMPC_lbIdx33, smoothMPC_rilb33, smoothMPC_dslbaff33);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb33, smoothMPC_dslbaff33, smoothMPC_llb33, smoothMPC_dllbaff33);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub33, smoothMPC_dzaff33, smoothMPC_ubIdx33, smoothMPC_dsubaff33);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub33, smoothMPC_dsubaff33, smoothMPC_lub33, smoothMPC_dlubaff33);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff34, smoothMPC_lbIdx34, smoothMPC_rilb34, smoothMPC_dslbaff34);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb34, smoothMPC_dslbaff34, smoothMPC_llb34, smoothMPC_dllbaff34);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub34, smoothMPC_dzaff34, smoothMPC_ubIdx34, smoothMPC_dsubaff34);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub34, smoothMPC_dsubaff34, smoothMPC_lub34, smoothMPC_dlubaff34);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff35, smoothMPC_lbIdx35, smoothMPC_rilb35, smoothMPC_dslbaff35);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb35, smoothMPC_dslbaff35, smoothMPC_llb35, smoothMPC_dllbaff35);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub35, smoothMPC_dzaff35, smoothMPC_ubIdx35, smoothMPC_dsubaff35);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub35, smoothMPC_dsubaff35, smoothMPC_lub35, smoothMPC_dlubaff35);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff36, smoothMPC_lbIdx36, smoothMPC_rilb36, smoothMPC_dslbaff36);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb36, smoothMPC_dslbaff36, smoothMPC_llb36, smoothMPC_dllbaff36);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub36, smoothMPC_dzaff36, smoothMPC_ubIdx36, smoothMPC_dsubaff36);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub36, smoothMPC_dsubaff36, smoothMPC_lub36, smoothMPC_dlubaff36);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff37, smoothMPC_lbIdx37, smoothMPC_rilb37, smoothMPC_dslbaff37);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb37, smoothMPC_dslbaff37, smoothMPC_llb37, smoothMPC_dllbaff37);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub37, smoothMPC_dzaff37, smoothMPC_ubIdx37, smoothMPC_dsubaff37);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub37, smoothMPC_dsubaff37, smoothMPC_lub37, smoothMPC_dlubaff37);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff38, smoothMPC_lbIdx38, smoothMPC_rilb38, smoothMPC_dslbaff38);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb38, smoothMPC_dslbaff38, smoothMPC_llb38, smoothMPC_dllbaff38);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub38, smoothMPC_dzaff38, smoothMPC_ubIdx38, smoothMPC_dsubaff38);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub38, smoothMPC_dsubaff38, smoothMPC_lub38, smoothMPC_dlubaff38);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff39, smoothMPC_lbIdx39, smoothMPC_rilb39, smoothMPC_dslbaff39);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb39, smoothMPC_dslbaff39, smoothMPC_llb39, smoothMPC_dllbaff39);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub39, smoothMPC_dzaff39, smoothMPC_ubIdx39, smoothMPC_dsubaff39);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub39, smoothMPC_dsubaff39, smoothMPC_lub39, smoothMPC_dlubaff39);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff40, smoothMPC_lbIdx40, smoothMPC_rilb40, smoothMPC_dslbaff40);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb40, smoothMPC_dslbaff40, smoothMPC_llb40, smoothMPC_dllbaff40);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub40, smoothMPC_dzaff40, smoothMPC_ubIdx40, smoothMPC_dsubaff40);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub40, smoothMPC_dsubaff40, smoothMPC_lub40, smoothMPC_dlubaff40);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff41, smoothMPC_lbIdx41, smoothMPC_rilb41, smoothMPC_dslbaff41);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb41, smoothMPC_dslbaff41, smoothMPC_llb41, smoothMPC_dllbaff41);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub41, smoothMPC_dzaff41, smoothMPC_ubIdx41, smoothMPC_dsubaff41);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub41, smoothMPC_dsubaff41, smoothMPC_lub41, smoothMPC_dlubaff41);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff42, smoothMPC_lbIdx42, smoothMPC_rilb42, smoothMPC_dslbaff42);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb42, smoothMPC_dslbaff42, smoothMPC_llb42, smoothMPC_dllbaff42);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub42, smoothMPC_dzaff42, smoothMPC_ubIdx42, smoothMPC_dsubaff42);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub42, smoothMPC_dsubaff42, smoothMPC_lub42, smoothMPC_dlubaff42);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff43, smoothMPC_lbIdx43, smoothMPC_rilb43, smoothMPC_dslbaff43);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb43, smoothMPC_dslbaff43, smoothMPC_llb43, smoothMPC_dllbaff43);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub43, smoothMPC_dzaff43, smoothMPC_ubIdx43, smoothMPC_dsubaff43);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub43, smoothMPC_dsubaff43, smoothMPC_lub43, smoothMPC_dlubaff43);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff44, smoothMPC_lbIdx44, smoothMPC_rilb44, smoothMPC_dslbaff44);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb44, smoothMPC_dslbaff44, smoothMPC_llb44, smoothMPC_dllbaff44);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub44, smoothMPC_dzaff44, smoothMPC_ubIdx44, smoothMPC_dsubaff44);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub44, smoothMPC_dsubaff44, smoothMPC_lub44, smoothMPC_dlubaff44);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff45, smoothMPC_lbIdx45, smoothMPC_rilb45, smoothMPC_dslbaff45);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb45, smoothMPC_dslbaff45, smoothMPC_llb45, smoothMPC_dllbaff45);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub45, smoothMPC_dzaff45, smoothMPC_ubIdx45, smoothMPC_dsubaff45);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub45, smoothMPC_dsubaff45, smoothMPC_lub45, smoothMPC_dlubaff45);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff46, smoothMPC_lbIdx46, smoothMPC_rilb46, smoothMPC_dslbaff46);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb46, smoothMPC_dslbaff46, smoothMPC_llb46, smoothMPC_dllbaff46);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub46, smoothMPC_dzaff46, smoothMPC_ubIdx46, smoothMPC_dsubaff46);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub46, smoothMPC_dsubaff46, smoothMPC_lub46, smoothMPC_dlubaff46);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff47, smoothMPC_lbIdx47, smoothMPC_rilb47, smoothMPC_dslbaff47);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb47, smoothMPC_dslbaff47, smoothMPC_llb47, smoothMPC_dllbaff47);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub47, smoothMPC_dzaff47, smoothMPC_ubIdx47, smoothMPC_dsubaff47);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub47, smoothMPC_dsubaff47, smoothMPC_lub47, smoothMPC_dlubaff47);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff48, smoothMPC_lbIdx48, smoothMPC_rilb48, smoothMPC_dslbaff48);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb48, smoothMPC_dslbaff48, smoothMPC_llb48, smoothMPC_dllbaff48);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub48, smoothMPC_dzaff48, smoothMPC_ubIdx48, smoothMPC_dsubaff48);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub48, smoothMPC_dsubaff48, smoothMPC_lub48, smoothMPC_dlubaff48);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff49, smoothMPC_lbIdx49, smoothMPC_rilb49, smoothMPC_dslbaff49);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb49, smoothMPC_dslbaff49, smoothMPC_llb49, smoothMPC_dllbaff49);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub49, smoothMPC_dzaff49, smoothMPC_ubIdx49, smoothMPC_dsubaff49);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub49, smoothMPC_dsubaff49, smoothMPC_lub49, smoothMPC_dlubaff49);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff50, smoothMPC_lbIdx50, smoothMPC_rilb50, smoothMPC_dslbaff50);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb50, smoothMPC_dslbaff50, smoothMPC_llb50, smoothMPC_dllbaff50);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub50, smoothMPC_dzaff50, smoothMPC_ubIdx50, smoothMPC_dsubaff50);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub50, smoothMPC_dsubaff50, smoothMPC_lub50, smoothMPC_dlubaff50);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff51, smoothMPC_lbIdx51, smoothMPC_rilb51, smoothMPC_dslbaff51);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb51, smoothMPC_dslbaff51, smoothMPC_llb51, smoothMPC_dllbaff51);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub51, smoothMPC_dzaff51, smoothMPC_ubIdx51, smoothMPC_dsubaff51);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub51, smoothMPC_dsubaff51, smoothMPC_lub51, smoothMPC_dlubaff51);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff52, smoothMPC_lbIdx52, smoothMPC_rilb52, smoothMPC_dslbaff52);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb52, smoothMPC_dslbaff52, smoothMPC_llb52, smoothMPC_dllbaff52);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub52, smoothMPC_dzaff52, smoothMPC_ubIdx52, smoothMPC_dsubaff52);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub52, smoothMPC_dsubaff52, smoothMPC_lub52, smoothMPC_dlubaff52);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff53, smoothMPC_lbIdx53, smoothMPC_rilb53, smoothMPC_dslbaff53);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb53, smoothMPC_dslbaff53, smoothMPC_llb53, smoothMPC_dllbaff53);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub53, smoothMPC_dzaff53, smoothMPC_ubIdx53, smoothMPC_dsubaff53);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub53, smoothMPC_dsubaff53, smoothMPC_lub53, smoothMPC_dlubaff53);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff54, smoothMPC_lbIdx54, smoothMPC_rilb54, smoothMPC_dslbaff54);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb54, smoothMPC_dslbaff54, smoothMPC_llb54, smoothMPC_dllbaff54);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub54, smoothMPC_dzaff54, smoothMPC_ubIdx54, smoothMPC_dsubaff54);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub54, smoothMPC_dsubaff54, smoothMPC_lub54, smoothMPC_dlubaff54);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff55, smoothMPC_lbIdx55, smoothMPC_rilb55, smoothMPC_dslbaff55);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb55, smoothMPC_dslbaff55, smoothMPC_llb55, smoothMPC_dllbaff55);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub55, smoothMPC_dzaff55, smoothMPC_ubIdx55, smoothMPC_dsubaff55);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub55, smoothMPC_dsubaff55, smoothMPC_lub55, smoothMPC_dlubaff55);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff56, smoothMPC_lbIdx56, smoothMPC_rilb56, smoothMPC_dslbaff56);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb56, smoothMPC_dslbaff56, smoothMPC_llb56, smoothMPC_dllbaff56);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub56, smoothMPC_dzaff56, smoothMPC_ubIdx56, smoothMPC_dsubaff56);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub56, smoothMPC_dsubaff56, smoothMPC_lub56, smoothMPC_dlubaff56);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff57, smoothMPC_lbIdx57, smoothMPC_rilb57, smoothMPC_dslbaff57);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb57, smoothMPC_dslbaff57, smoothMPC_llb57, smoothMPC_dllbaff57);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub57, smoothMPC_dzaff57, smoothMPC_ubIdx57, smoothMPC_dsubaff57);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub57, smoothMPC_dsubaff57, smoothMPC_lub57, smoothMPC_dlubaff57);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff58, smoothMPC_lbIdx58, smoothMPC_rilb58, smoothMPC_dslbaff58);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb58, smoothMPC_dslbaff58, smoothMPC_llb58, smoothMPC_dllbaff58);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub58, smoothMPC_dzaff58, smoothMPC_ubIdx58, smoothMPC_dsubaff58);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub58, smoothMPC_dsubaff58, smoothMPC_lub58, smoothMPC_dlubaff58);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff59, smoothMPC_lbIdx59, smoothMPC_rilb59, smoothMPC_dslbaff59);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb59, smoothMPC_dslbaff59, smoothMPC_llb59, smoothMPC_dllbaff59);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub59, smoothMPC_dzaff59, smoothMPC_ubIdx59, smoothMPC_dsubaff59);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub59, smoothMPC_dsubaff59, smoothMPC_lub59, smoothMPC_dlubaff59);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff60, smoothMPC_lbIdx60, smoothMPC_rilb60, smoothMPC_dslbaff60);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb60, smoothMPC_dslbaff60, smoothMPC_llb60, smoothMPC_dllbaff60);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub60, smoothMPC_dzaff60, smoothMPC_ubIdx60, smoothMPC_dsubaff60);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub60, smoothMPC_dsubaff60, smoothMPC_lub60, smoothMPC_dlubaff60);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff61, smoothMPC_lbIdx61, smoothMPC_rilb61, smoothMPC_dslbaff61);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb61, smoothMPC_dslbaff61, smoothMPC_llb61, smoothMPC_dllbaff61);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub61, smoothMPC_dzaff61, smoothMPC_ubIdx61, smoothMPC_dsubaff61);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub61, smoothMPC_dsubaff61, smoothMPC_lub61, smoothMPC_dlubaff61);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff62, smoothMPC_lbIdx62, smoothMPC_rilb62, smoothMPC_dslbaff62);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb62, smoothMPC_dslbaff62, smoothMPC_llb62, smoothMPC_dllbaff62);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub62, smoothMPC_dzaff62, smoothMPC_ubIdx62, smoothMPC_dsubaff62);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub62, smoothMPC_dsubaff62, smoothMPC_lub62, smoothMPC_dlubaff62);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff63, smoothMPC_lbIdx63, smoothMPC_rilb63, smoothMPC_dslbaff63);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb63, smoothMPC_dslbaff63, smoothMPC_llb63, smoothMPC_dllbaff63);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub63, smoothMPC_dzaff63, smoothMPC_ubIdx63, smoothMPC_dsubaff63);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub63, smoothMPC_dsubaff63, smoothMPC_lub63, smoothMPC_dlubaff63);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff64, smoothMPC_lbIdx64, smoothMPC_rilb64, smoothMPC_dslbaff64);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb64, smoothMPC_dslbaff64, smoothMPC_llb64, smoothMPC_dllbaff64);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub64, smoothMPC_dzaff64, smoothMPC_ubIdx64, smoothMPC_dsubaff64);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub64, smoothMPC_dsubaff64, smoothMPC_lub64, smoothMPC_dlubaff64);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff65, smoothMPC_lbIdx65, smoothMPC_rilb65, smoothMPC_dslbaff65);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb65, smoothMPC_dslbaff65, smoothMPC_llb65, smoothMPC_dllbaff65);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub65, smoothMPC_dzaff65, smoothMPC_ubIdx65, smoothMPC_dsubaff65);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub65, smoothMPC_dsubaff65, smoothMPC_lub65, smoothMPC_dlubaff65);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff66, smoothMPC_lbIdx66, smoothMPC_rilb66, smoothMPC_dslbaff66);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb66, smoothMPC_dslbaff66, smoothMPC_llb66, smoothMPC_dllbaff66);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub66, smoothMPC_dzaff66, smoothMPC_ubIdx66, smoothMPC_dsubaff66);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub66, smoothMPC_dsubaff66, smoothMPC_lub66, smoothMPC_dlubaff66);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff67, smoothMPC_lbIdx67, smoothMPC_rilb67, smoothMPC_dslbaff67);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb67, smoothMPC_dslbaff67, smoothMPC_llb67, smoothMPC_dllbaff67);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub67, smoothMPC_dzaff67, smoothMPC_ubIdx67, smoothMPC_dsubaff67);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub67, smoothMPC_dsubaff67, smoothMPC_lub67, smoothMPC_dlubaff67);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff68, smoothMPC_lbIdx68, smoothMPC_rilb68, smoothMPC_dslbaff68);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb68, smoothMPC_dslbaff68, smoothMPC_llb68, smoothMPC_dllbaff68);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub68, smoothMPC_dzaff68, smoothMPC_ubIdx68, smoothMPC_dsubaff68);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub68, smoothMPC_dsubaff68, smoothMPC_lub68, smoothMPC_dlubaff68);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff69, smoothMPC_lbIdx69, smoothMPC_rilb69, smoothMPC_dslbaff69);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb69, smoothMPC_dslbaff69, smoothMPC_llb69, smoothMPC_dllbaff69);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub69, smoothMPC_dzaff69, smoothMPC_ubIdx69, smoothMPC_dsubaff69);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub69, smoothMPC_dsubaff69, smoothMPC_lub69, smoothMPC_dlubaff69);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff70, smoothMPC_lbIdx70, smoothMPC_rilb70, smoothMPC_dslbaff70);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb70, smoothMPC_dslbaff70, smoothMPC_llb70, smoothMPC_dllbaff70);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub70, smoothMPC_dzaff70, smoothMPC_ubIdx70, smoothMPC_dsubaff70);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub70, smoothMPC_dsubaff70, smoothMPC_lub70, smoothMPC_dlubaff70);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff71, smoothMPC_lbIdx71, smoothMPC_rilb71, smoothMPC_dslbaff71);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb71, smoothMPC_dslbaff71, smoothMPC_llb71, smoothMPC_dllbaff71);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub71, smoothMPC_dzaff71, smoothMPC_ubIdx71, smoothMPC_dsubaff71);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub71, smoothMPC_dsubaff71, smoothMPC_lub71, smoothMPC_dlubaff71);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff72, smoothMPC_lbIdx72, smoothMPC_rilb72, smoothMPC_dslbaff72);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb72, smoothMPC_dslbaff72, smoothMPC_llb72, smoothMPC_dllbaff72);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub72, smoothMPC_dzaff72, smoothMPC_ubIdx72, smoothMPC_dsubaff72);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub72, smoothMPC_dsubaff72, smoothMPC_lub72, smoothMPC_dlubaff72);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff73, smoothMPC_lbIdx73, smoothMPC_rilb73, smoothMPC_dslbaff73);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb73, smoothMPC_dslbaff73, smoothMPC_llb73, smoothMPC_dllbaff73);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub73, smoothMPC_dzaff73, smoothMPC_ubIdx73, smoothMPC_dsubaff73);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub73, smoothMPC_dsubaff73, smoothMPC_lub73, smoothMPC_dlubaff73);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff74, smoothMPC_lbIdx74, smoothMPC_rilb74, smoothMPC_dslbaff74);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb74, smoothMPC_dslbaff74, smoothMPC_llb74, smoothMPC_dllbaff74);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub74, smoothMPC_dzaff74, smoothMPC_ubIdx74, smoothMPC_dsubaff74);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub74, smoothMPC_dsubaff74, smoothMPC_lub74, smoothMPC_dlubaff74);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff75, smoothMPC_lbIdx75, smoothMPC_rilb75, smoothMPC_dslbaff75);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb75, smoothMPC_dslbaff75, smoothMPC_llb75, smoothMPC_dllbaff75);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub75, smoothMPC_dzaff75, smoothMPC_ubIdx75, smoothMPC_dsubaff75);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub75, smoothMPC_dsubaff75, smoothMPC_lub75, smoothMPC_dlubaff75);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff76, smoothMPC_lbIdx76, smoothMPC_rilb76, smoothMPC_dslbaff76);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb76, smoothMPC_dslbaff76, smoothMPC_llb76, smoothMPC_dllbaff76);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub76, smoothMPC_dzaff76, smoothMPC_ubIdx76, smoothMPC_dsubaff76);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub76, smoothMPC_dsubaff76, smoothMPC_lub76, smoothMPC_dlubaff76);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff77, smoothMPC_lbIdx77, smoothMPC_rilb77, smoothMPC_dslbaff77);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb77, smoothMPC_dslbaff77, smoothMPC_llb77, smoothMPC_dllbaff77);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub77, smoothMPC_dzaff77, smoothMPC_ubIdx77, smoothMPC_dsubaff77);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub77, smoothMPC_dsubaff77, smoothMPC_lub77, smoothMPC_dlubaff77);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff78, smoothMPC_lbIdx78, smoothMPC_rilb78, smoothMPC_dslbaff78);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb78, smoothMPC_dslbaff78, smoothMPC_llb78, smoothMPC_dllbaff78);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub78, smoothMPC_dzaff78, smoothMPC_ubIdx78, smoothMPC_dsubaff78);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub78, smoothMPC_dsubaff78, smoothMPC_lub78, smoothMPC_dlubaff78);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff79, smoothMPC_lbIdx79, smoothMPC_rilb79, smoothMPC_dslbaff79);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb79, smoothMPC_dslbaff79, smoothMPC_llb79, smoothMPC_dllbaff79);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub79, smoothMPC_dzaff79, smoothMPC_ubIdx79, smoothMPC_dsubaff79);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub79, smoothMPC_dsubaff79, smoothMPC_lub79, smoothMPC_dlubaff79);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff80, smoothMPC_lbIdx80, smoothMPC_rilb80, smoothMPC_dslbaff80);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb80, smoothMPC_dslbaff80, smoothMPC_llb80, smoothMPC_dllbaff80);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub80, smoothMPC_dzaff80, smoothMPC_ubIdx80, smoothMPC_dsubaff80);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub80, smoothMPC_dsubaff80, smoothMPC_lub80, smoothMPC_dlubaff80);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff81, smoothMPC_lbIdx81, smoothMPC_rilb81, smoothMPC_dslbaff81);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb81, smoothMPC_dslbaff81, smoothMPC_llb81, smoothMPC_dllbaff81);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub81, smoothMPC_dzaff81, smoothMPC_ubIdx81, smoothMPC_dsubaff81);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub81, smoothMPC_dsubaff81, smoothMPC_lub81, smoothMPC_dlubaff81);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff82, smoothMPC_lbIdx82, smoothMPC_rilb82, smoothMPC_dslbaff82);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb82, smoothMPC_dslbaff82, smoothMPC_llb82, smoothMPC_dllbaff82);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub82, smoothMPC_dzaff82, smoothMPC_ubIdx82, smoothMPC_dsubaff82);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub82, smoothMPC_dsubaff82, smoothMPC_lub82, smoothMPC_dlubaff82);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff83, smoothMPC_lbIdx83, smoothMPC_rilb83, smoothMPC_dslbaff83);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb83, smoothMPC_dslbaff83, smoothMPC_llb83, smoothMPC_dllbaff83);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub83, smoothMPC_dzaff83, smoothMPC_ubIdx83, smoothMPC_dsubaff83);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub83, smoothMPC_dsubaff83, smoothMPC_lub83, smoothMPC_dlubaff83);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff84, smoothMPC_lbIdx84, smoothMPC_rilb84, smoothMPC_dslbaff84);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb84, smoothMPC_dslbaff84, smoothMPC_llb84, smoothMPC_dllbaff84);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub84, smoothMPC_dzaff84, smoothMPC_ubIdx84, smoothMPC_dsubaff84);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub84, smoothMPC_dsubaff84, smoothMPC_lub84, smoothMPC_dlubaff84);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff85, smoothMPC_lbIdx85, smoothMPC_rilb85, smoothMPC_dslbaff85);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb85, smoothMPC_dslbaff85, smoothMPC_llb85, smoothMPC_dllbaff85);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub85, smoothMPC_dzaff85, smoothMPC_ubIdx85, smoothMPC_dsubaff85);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub85, smoothMPC_dsubaff85, smoothMPC_lub85, smoothMPC_dlubaff85);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff86, smoothMPC_lbIdx86, smoothMPC_rilb86, smoothMPC_dslbaff86);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb86, smoothMPC_dslbaff86, smoothMPC_llb86, smoothMPC_dllbaff86);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub86, smoothMPC_dzaff86, smoothMPC_ubIdx86, smoothMPC_dsubaff86);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub86, smoothMPC_dsubaff86, smoothMPC_lub86, smoothMPC_dlubaff86);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff87, smoothMPC_lbIdx87, smoothMPC_rilb87, smoothMPC_dslbaff87);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb87, smoothMPC_dslbaff87, smoothMPC_llb87, smoothMPC_dllbaff87);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub87, smoothMPC_dzaff87, smoothMPC_ubIdx87, smoothMPC_dsubaff87);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub87, smoothMPC_dsubaff87, smoothMPC_lub87, smoothMPC_dlubaff87);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff88, smoothMPC_lbIdx88, smoothMPC_rilb88, smoothMPC_dslbaff88);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb88, smoothMPC_dslbaff88, smoothMPC_llb88, smoothMPC_dllbaff88);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub88, smoothMPC_dzaff88, smoothMPC_ubIdx88, smoothMPC_dsubaff88);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub88, smoothMPC_dsubaff88, smoothMPC_lub88, smoothMPC_dlubaff88);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff89, smoothMPC_lbIdx89, smoothMPC_rilb89, smoothMPC_dslbaff89);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb89, smoothMPC_dslbaff89, smoothMPC_llb89, smoothMPC_dllbaff89);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub89, smoothMPC_dzaff89, smoothMPC_ubIdx89, smoothMPC_dsubaff89);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub89, smoothMPC_dsubaff89, smoothMPC_lub89, smoothMPC_dlubaff89);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff90, smoothMPC_lbIdx90, smoothMPC_rilb90, smoothMPC_dslbaff90);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb90, smoothMPC_dslbaff90, smoothMPC_llb90, smoothMPC_dllbaff90);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub90, smoothMPC_dzaff90, smoothMPC_ubIdx90, smoothMPC_dsubaff90);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub90, smoothMPC_dsubaff90, smoothMPC_lub90, smoothMPC_dlubaff90);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff91, smoothMPC_lbIdx91, smoothMPC_rilb91, smoothMPC_dslbaff91);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb91, smoothMPC_dslbaff91, smoothMPC_llb91, smoothMPC_dllbaff91);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub91, smoothMPC_dzaff91, smoothMPC_ubIdx91, smoothMPC_dsubaff91);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub91, smoothMPC_dsubaff91, smoothMPC_lub91, smoothMPC_dlubaff91);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff92, smoothMPC_lbIdx92, smoothMPC_rilb92, smoothMPC_dslbaff92);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb92, smoothMPC_dslbaff92, smoothMPC_llb92, smoothMPC_dllbaff92);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub92, smoothMPC_dzaff92, smoothMPC_ubIdx92, smoothMPC_dsubaff92);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub92, smoothMPC_dsubaff92, smoothMPC_lub92, smoothMPC_dlubaff92);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff93, smoothMPC_lbIdx93, smoothMPC_rilb93, smoothMPC_dslbaff93);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb93, smoothMPC_dslbaff93, smoothMPC_llb93, smoothMPC_dllbaff93);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub93, smoothMPC_dzaff93, smoothMPC_ubIdx93, smoothMPC_dsubaff93);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub93, smoothMPC_dsubaff93, smoothMPC_lub93, smoothMPC_dlubaff93);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff94, smoothMPC_lbIdx94, smoothMPC_rilb94, smoothMPC_dslbaff94);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb94, smoothMPC_dslbaff94, smoothMPC_llb94, smoothMPC_dllbaff94);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub94, smoothMPC_dzaff94, smoothMPC_ubIdx94, smoothMPC_dsubaff94);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub94, smoothMPC_dsubaff94, smoothMPC_lub94, smoothMPC_dlubaff94);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff95, smoothMPC_lbIdx95, smoothMPC_rilb95, smoothMPC_dslbaff95);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb95, smoothMPC_dslbaff95, smoothMPC_llb95, smoothMPC_dllbaff95);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub95, smoothMPC_dzaff95, smoothMPC_ubIdx95, smoothMPC_dsubaff95);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub95, smoothMPC_dsubaff95, smoothMPC_lub95, smoothMPC_dlubaff95);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff96, smoothMPC_lbIdx96, smoothMPC_rilb96, smoothMPC_dslbaff96);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb96, smoothMPC_dslbaff96, smoothMPC_llb96, smoothMPC_dllbaff96);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub96, smoothMPC_dzaff96, smoothMPC_ubIdx96, smoothMPC_dsubaff96);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub96, smoothMPC_dsubaff96, smoothMPC_lub96, smoothMPC_dlubaff96);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff97, smoothMPC_lbIdx97, smoothMPC_rilb97, smoothMPC_dslbaff97);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb97, smoothMPC_dslbaff97, smoothMPC_llb97, smoothMPC_dllbaff97);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub97, smoothMPC_dzaff97, smoothMPC_ubIdx97, smoothMPC_dsubaff97);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub97, smoothMPC_dsubaff97, smoothMPC_lub97, smoothMPC_dlubaff97);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff98, smoothMPC_lbIdx98, smoothMPC_rilb98, smoothMPC_dslbaff98);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb98, smoothMPC_dslbaff98, smoothMPC_llb98, smoothMPC_dllbaff98);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub98, smoothMPC_dzaff98, smoothMPC_ubIdx98, smoothMPC_dsubaff98);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub98, smoothMPC_dsubaff98, smoothMPC_lub98, smoothMPC_dlubaff98);
smoothMPC_LA_VSUB_INDEXED_3(smoothMPC_dzaff99, smoothMPC_lbIdx99, smoothMPC_rilb99, smoothMPC_dslbaff99);
smoothMPC_LA_VSUB3_3(smoothMPC_llbbyslb99, smoothMPC_dslbaff99, smoothMPC_llb99, smoothMPC_dllbaff99);
smoothMPC_LA_VSUB2_INDEXED_3(smoothMPC_riub99, smoothMPC_dzaff99, smoothMPC_ubIdx99, smoothMPC_dsubaff99);
smoothMPC_LA_VSUB3_3(smoothMPC_lubbysub99, smoothMPC_dsubaff99, smoothMPC_lub99, smoothMPC_dlubaff99);
info->lsit_aff = smoothMPC_LINESEARCH_BACKTRACKING_AFFINE(smoothMPC_l, smoothMPC_s, smoothMPC_dl_aff, smoothMPC_ds_aff, &info->step_aff, &info->mu_aff);
if( info->lsit_aff == smoothMPC_NOPROGRESS ){
PRINTTEXT("Affine line search could not proceed at iteration %d.\nThe problem might be infeasible -- exiting.\n",info->it+1);
exitcode = smoothMPC_NOPROGRESS; break;
}
sigma_3rdroot = info->mu_aff / info->mu;
info->sigma = sigma_3rdroot*sigma_3rdroot*sigma_3rdroot;
musigma = info->mu * info->sigma;
smoothMPC_LA_VSUB5_1590(smoothMPC_ds_aff, smoothMPC_dl_aff, musigma, smoothMPC_ccrhs);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub00, smoothMPC_sub00, smoothMPC_ubIdx00, smoothMPC_ccrhsl00, smoothMPC_slb00, smoothMPC_lbIdx00, smoothMPC_rd00);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub01, smoothMPC_sub01, smoothMPC_ubIdx01, smoothMPC_ccrhsl01, smoothMPC_slb01, smoothMPC_lbIdx01, smoothMPC_rd01);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi00, smoothMPC_rd00, smoothMPC_Lbyrd00);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi01, smoothMPC_rd01, smoothMPC_Lbyrd01);
smoothMPC_LA_DIAGZERO_MVM_3(smoothMPC_W00, smoothMPC_Lbyrd00, smoothMPC_beta00);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld00, smoothMPC_beta00, smoothMPC_yy00);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V00, smoothMPC_Lbyrd00, smoothMPC_W01, smoothMPC_Lbyrd01, smoothMPC_beta01);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd01, smoothMPC_yy00, smoothMPC_beta01, smoothMPC_bmy01);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld01, smoothMPC_bmy01, smoothMPC_yy01);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub02, smoothMPC_sub02, smoothMPC_ubIdx02, smoothMPC_ccrhsl02, smoothMPC_slb02, smoothMPC_lbIdx02, smoothMPC_rd02);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi02, smoothMPC_rd02, smoothMPC_Lbyrd02);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V01, smoothMPC_Lbyrd01, smoothMPC_W02, smoothMPC_Lbyrd02, smoothMPC_beta02);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd02, smoothMPC_yy01, smoothMPC_beta02, smoothMPC_bmy02);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld02, smoothMPC_bmy02, smoothMPC_yy02);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub03, smoothMPC_sub03, smoothMPC_ubIdx03, smoothMPC_ccrhsl03, smoothMPC_slb03, smoothMPC_lbIdx03, smoothMPC_rd03);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi03, smoothMPC_rd03, smoothMPC_Lbyrd03);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V02, smoothMPC_Lbyrd02, smoothMPC_W03, smoothMPC_Lbyrd03, smoothMPC_beta03);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd03, smoothMPC_yy02, smoothMPC_beta03, smoothMPC_bmy03);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld03, smoothMPC_bmy03, smoothMPC_yy03);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub04, smoothMPC_sub04, smoothMPC_ubIdx04, smoothMPC_ccrhsl04, smoothMPC_slb04, smoothMPC_lbIdx04, smoothMPC_rd04);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi04, smoothMPC_rd04, smoothMPC_Lbyrd04);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V03, smoothMPC_Lbyrd03, smoothMPC_W04, smoothMPC_Lbyrd04, smoothMPC_beta04);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd04, smoothMPC_yy03, smoothMPC_beta04, smoothMPC_bmy04);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld04, smoothMPC_bmy04, smoothMPC_yy04);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub05, smoothMPC_sub05, smoothMPC_ubIdx05, smoothMPC_ccrhsl05, smoothMPC_slb05, smoothMPC_lbIdx05, smoothMPC_rd05);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi05, smoothMPC_rd05, smoothMPC_Lbyrd05);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V04, smoothMPC_Lbyrd04, smoothMPC_W05, smoothMPC_Lbyrd05, smoothMPC_beta05);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd05, smoothMPC_yy04, smoothMPC_beta05, smoothMPC_bmy05);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld05, smoothMPC_bmy05, smoothMPC_yy05);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub06, smoothMPC_sub06, smoothMPC_ubIdx06, smoothMPC_ccrhsl06, smoothMPC_slb06, smoothMPC_lbIdx06, smoothMPC_rd06);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi06, smoothMPC_rd06, smoothMPC_Lbyrd06);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V05, smoothMPC_Lbyrd05, smoothMPC_W06, smoothMPC_Lbyrd06, smoothMPC_beta06);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd06, smoothMPC_yy05, smoothMPC_beta06, smoothMPC_bmy06);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld06, smoothMPC_bmy06, smoothMPC_yy06);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub07, smoothMPC_sub07, smoothMPC_ubIdx07, smoothMPC_ccrhsl07, smoothMPC_slb07, smoothMPC_lbIdx07, smoothMPC_rd07);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi07, smoothMPC_rd07, smoothMPC_Lbyrd07);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V06, smoothMPC_Lbyrd06, smoothMPC_W07, smoothMPC_Lbyrd07, smoothMPC_beta07);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd07, smoothMPC_yy06, smoothMPC_beta07, smoothMPC_bmy07);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld07, smoothMPC_bmy07, smoothMPC_yy07);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub08, smoothMPC_sub08, smoothMPC_ubIdx08, smoothMPC_ccrhsl08, smoothMPC_slb08, smoothMPC_lbIdx08, smoothMPC_rd08);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi08, smoothMPC_rd08, smoothMPC_Lbyrd08);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V07, smoothMPC_Lbyrd07, smoothMPC_W08, smoothMPC_Lbyrd08, smoothMPC_beta08);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd08, smoothMPC_yy07, smoothMPC_beta08, smoothMPC_bmy08);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld08, smoothMPC_bmy08, smoothMPC_yy08);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub09, smoothMPC_sub09, smoothMPC_ubIdx09, smoothMPC_ccrhsl09, smoothMPC_slb09, smoothMPC_lbIdx09, smoothMPC_rd09);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi09, smoothMPC_rd09, smoothMPC_Lbyrd09);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V08, smoothMPC_Lbyrd08, smoothMPC_W09, smoothMPC_Lbyrd09, smoothMPC_beta09);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd09, smoothMPC_yy08, smoothMPC_beta09, smoothMPC_bmy09);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld09, smoothMPC_bmy09, smoothMPC_yy09);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub10, smoothMPC_sub10, smoothMPC_ubIdx10, smoothMPC_ccrhsl10, smoothMPC_slb10, smoothMPC_lbIdx10, smoothMPC_rd10);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi10, smoothMPC_rd10, smoothMPC_Lbyrd10);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V09, smoothMPC_Lbyrd09, smoothMPC_W10, smoothMPC_Lbyrd10, smoothMPC_beta10);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd10, smoothMPC_yy09, smoothMPC_beta10, smoothMPC_bmy10);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld10, smoothMPC_bmy10, smoothMPC_yy10);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub11, smoothMPC_sub11, smoothMPC_ubIdx11, smoothMPC_ccrhsl11, smoothMPC_slb11, smoothMPC_lbIdx11, smoothMPC_rd11);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi11, smoothMPC_rd11, smoothMPC_Lbyrd11);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V10, smoothMPC_Lbyrd10, smoothMPC_W11, smoothMPC_Lbyrd11, smoothMPC_beta11);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd11, smoothMPC_yy10, smoothMPC_beta11, smoothMPC_bmy11);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld11, smoothMPC_bmy11, smoothMPC_yy11);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub12, smoothMPC_sub12, smoothMPC_ubIdx12, smoothMPC_ccrhsl12, smoothMPC_slb12, smoothMPC_lbIdx12, smoothMPC_rd12);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi12, smoothMPC_rd12, smoothMPC_Lbyrd12);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V11, smoothMPC_Lbyrd11, smoothMPC_W12, smoothMPC_Lbyrd12, smoothMPC_beta12);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd12, smoothMPC_yy11, smoothMPC_beta12, smoothMPC_bmy12);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld12, smoothMPC_bmy12, smoothMPC_yy12);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub13, smoothMPC_sub13, smoothMPC_ubIdx13, smoothMPC_ccrhsl13, smoothMPC_slb13, smoothMPC_lbIdx13, smoothMPC_rd13);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi13, smoothMPC_rd13, smoothMPC_Lbyrd13);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V12, smoothMPC_Lbyrd12, smoothMPC_W13, smoothMPC_Lbyrd13, smoothMPC_beta13);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd13, smoothMPC_yy12, smoothMPC_beta13, smoothMPC_bmy13);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld13, smoothMPC_bmy13, smoothMPC_yy13);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub14, smoothMPC_sub14, smoothMPC_ubIdx14, smoothMPC_ccrhsl14, smoothMPC_slb14, smoothMPC_lbIdx14, smoothMPC_rd14);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi14, smoothMPC_rd14, smoothMPC_Lbyrd14);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V13, smoothMPC_Lbyrd13, smoothMPC_W14, smoothMPC_Lbyrd14, smoothMPC_beta14);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd14, smoothMPC_yy13, smoothMPC_beta14, smoothMPC_bmy14);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld14, smoothMPC_bmy14, smoothMPC_yy14);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub15, smoothMPC_sub15, smoothMPC_ubIdx15, smoothMPC_ccrhsl15, smoothMPC_slb15, smoothMPC_lbIdx15, smoothMPC_rd15);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi15, smoothMPC_rd15, smoothMPC_Lbyrd15);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V14, smoothMPC_Lbyrd14, smoothMPC_W15, smoothMPC_Lbyrd15, smoothMPC_beta15);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd15, smoothMPC_yy14, smoothMPC_beta15, smoothMPC_bmy15);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld15, smoothMPC_bmy15, smoothMPC_yy15);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub16, smoothMPC_sub16, smoothMPC_ubIdx16, smoothMPC_ccrhsl16, smoothMPC_slb16, smoothMPC_lbIdx16, smoothMPC_rd16);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi16, smoothMPC_rd16, smoothMPC_Lbyrd16);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V15, smoothMPC_Lbyrd15, smoothMPC_W16, smoothMPC_Lbyrd16, smoothMPC_beta16);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd16, smoothMPC_yy15, smoothMPC_beta16, smoothMPC_bmy16);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld16, smoothMPC_bmy16, smoothMPC_yy16);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub17, smoothMPC_sub17, smoothMPC_ubIdx17, smoothMPC_ccrhsl17, smoothMPC_slb17, smoothMPC_lbIdx17, smoothMPC_rd17);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi17, smoothMPC_rd17, smoothMPC_Lbyrd17);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V16, smoothMPC_Lbyrd16, smoothMPC_W17, smoothMPC_Lbyrd17, smoothMPC_beta17);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd17, smoothMPC_yy16, smoothMPC_beta17, smoothMPC_bmy17);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld17, smoothMPC_bmy17, smoothMPC_yy17);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub18, smoothMPC_sub18, smoothMPC_ubIdx18, smoothMPC_ccrhsl18, smoothMPC_slb18, smoothMPC_lbIdx18, smoothMPC_rd18);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi18, smoothMPC_rd18, smoothMPC_Lbyrd18);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V17, smoothMPC_Lbyrd17, smoothMPC_W18, smoothMPC_Lbyrd18, smoothMPC_beta18);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd18, smoothMPC_yy17, smoothMPC_beta18, smoothMPC_bmy18);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld18, smoothMPC_bmy18, smoothMPC_yy18);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub19, smoothMPC_sub19, smoothMPC_ubIdx19, smoothMPC_ccrhsl19, smoothMPC_slb19, smoothMPC_lbIdx19, smoothMPC_rd19);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi19, smoothMPC_rd19, smoothMPC_Lbyrd19);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V18, smoothMPC_Lbyrd18, smoothMPC_W19, smoothMPC_Lbyrd19, smoothMPC_beta19);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd19, smoothMPC_yy18, smoothMPC_beta19, smoothMPC_bmy19);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld19, smoothMPC_bmy19, smoothMPC_yy19);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub20, smoothMPC_sub20, smoothMPC_ubIdx20, smoothMPC_ccrhsl20, smoothMPC_slb20, smoothMPC_lbIdx20, smoothMPC_rd20);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi20, smoothMPC_rd20, smoothMPC_Lbyrd20);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V19, smoothMPC_Lbyrd19, smoothMPC_W20, smoothMPC_Lbyrd20, smoothMPC_beta20);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd20, smoothMPC_yy19, smoothMPC_beta20, smoothMPC_bmy20);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld20, smoothMPC_bmy20, smoothMPC_yy20);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub21, smoothMPC_sub21, smoothMPC_ubIdx21, smoothMPC_ccrhsl21, smoothMPC_slb21, smoothMPC_lbIdx21, smoothMPC_rd21);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi21, smoothMPC_rd21, smoothMPC_Lbyrd21);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V20, smoothMPC_Lbyrd20, smoothMPC_W21, smoothMPC_Lbyrd21, smoothMPC_beta21);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd21, smoothMPC_yy20, smoothMPC_beta21, smoothMPC_bmy21);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld21, smoothMPC_bmy21, smoothMPC_yy21);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub22, smoothMPC_sub22, smoothMPC_ubIdx22, smoothMPC_ccrhsl22, smoothMPC_slb22, smoothMPC_lbIdx22, smoothMPC_rd22);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi22, smoothMPC_rd22, smoothMPC_Lbyrd22);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V21, smoothMPC_Lbyrd21, smoothMPC_W22, smoothMPC_Lbyrd22, smoothMPC_beta22);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd22, smoothMPC_yy21, smoothMPC_beta22, smoothMPC_bmy22);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld22, smoothMPC_bmy22, smoothMPC_yy22);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub23, smoothMPC_sub23, smoothMPC_ubIdx23, smoothMPC_ccrhsl23, smoothMPC_slb23, smoothMPC_lbIdx23, smoothMPC_rd23);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi23, smoothMPC_rd23, smoothMPC_Lbyrd23);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V22, smoothMPC_Lbyrd22, smoothMPC_W23, smoothMPC_Lbyrd23, smoothMPC_beta23);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd23, smoothMPC_yy22, smoothMPC_beta23, smoothMPC_bmy23);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld23, smoothMPC_bmy23, smoothMPC_yy23);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub24, smoothMPC_sub24, smoothMPC_ubIdx24, smoothMPC_ccrhsl24, smoothMPC_slb24, smoothMPC_lbIdx24, smoothMPC_rd24);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi24, smoothMPC_rd24, smoothMPC_Lbyrd24);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V23, smoothMPC_Lbyrd23, smoothMPC_W24, smoothMPC_Lbyrd24, smoothMPC_beta24);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd24, smoothMPC_yy23, smoothMPC_beta24, smoothMPC_bmy24);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld24, smoothMPC_bmy24, smoothMPC_yy24);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub25, smoothMPC_sub25, smoothMPC_ubIdx25, smoothMPC_ccrhsl25, smoothMPC_slb25, smoothMPC_lbIdx25, smoothMPC_rd25);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi25, smoothMPC_rd25, smoothMPC_Lbyrd25);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V24, smoothMPC_Lbyrd24, smoothMPC_W25, smoothMPC_Lbyrd25, smoothMPC_beta25);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd25, smoothMPC_yy24, smoothMPC_beta25, smoothMPC_bmy25);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld25, smoothMPC_bmy25, smoothMPC_yy25);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub26, smoothMPC_sub26, smoothMPC_ubIdx26, smoothMPC_ccrhsl26, smoothMPC_slb26, smoothMPC_lbIdx26, smoothMPC_rd26);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi26, smoothMPC_rd26, smoothMPC_Lbyrd26);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V25, smoothMPC_Lbyrd25, smoothMPC_W26, smoothMPC_Lbyrd26, smoothMPC_beta26);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd26, smoothMPC_yy25, smoothMPC_beta26, smoothMPC_bmy26);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld26, smoothMPC_bmy26, smoothMPC_yy26);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub27, smoothMPC_sub27, smoothMPC_ubIdx27, smoothMPC_ccrhsl27, smoothMPC_slb27, smoothMPC_lbIdx27, smoothMPC_rd27);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi27, smoothMPC_rd27, smoothMPC_Lbyrd27);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V26, smoothMPC_Lbyrd26, smoothMPC_W27, smoothMPC_Lbyrd27, smoothMPC_beta27);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd27, smoothMPC_yy26, smoothMPC_beta27, smoothMPC_bmy27);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld27, smoothMPC_bmy27, smoothMPC_yy27);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub28, smoothMPC_sub28, smoothMPC_ubIdx28, smoothMPC_ccrhsl28, smoothMPC_slb28, smoothMPC_lbIdx28, smoothMPC_rd28);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi28, smoothMPC_rd28, smoothMPC_Lbyrd28);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V27, smoothMPC_Lbyrd27, smoothMPC_W28, smoothMPC_Lbyrd28, smoothMPC_beta28);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd28, smoothMPC_yy27, smoothMPC_beta28, smoothMPC_bmy28);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld28, smoothMPC_bmy28, smoothMPC_yy28);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub29, smoothMPC_sub29, smoothMPC_ubIdx29, smoothMPC_ccrhsl29, smoothMPC_slb29, smoothMPC_lbIdx29, smoothMPC_rd29);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi29, smoothMPC_rd29, smoothMPC_Lbyrd29);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V28, smoothMPC_Lbyrd28, smoothMPC_W29, smoothMPC_Lbyrd29, smoothMPC_beta29);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd29, smoothMPC_yy28, smoothMPC_beta29, smoothMPC_bmy29);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld29, smoothMPC_bmy29, smoothMPC_yy29);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub30, smoothMPC_sub30, smoothMPC_ubIdx30, smoothMPC_ccrhsl30, smoothMPC_slb30, smoothMPC_lbIdx30, smoothMPC_rd30);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi30, smoothMPC_rd30, smoothMPC_Lbyrd30);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V29, smoothMPC_Lbyrd29, smoothMPC_W30, smoothMPC_Lbyrd30, smoothMPC_beta30);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd30, smoothMPC_yy29, smoothMPC_beta30, smoothMPC_bmy30);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld30, smoothMPC_bmy30, smoothMPC_yy30);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub31, smoothMPC_sub31, smoothMPC_ubIdx31, smoothMPC_ccrhsl31, smoothMPC_slb31, smoothMPC_lbIdx31, smoothMPC_rd31);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi31, smoothMPC_rd31, smoothMPC_Lbyrd31);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V30, smoothMPC_Lbyrd30, smoothMPC_W31, smoothMPC_Lbyrd31, smoothMPC_beta31);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd31, smoothMPC_yy30, smoothMPC_beta31, smoothMPC_bmy31);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld31, smoothMPC_bmy31, smoothMPC_yy31);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub32, smoothMPC_sub32, smoothMPC_ubIdx32, smoothMPC_ccrhsl32, smoothMPC_slb32, smoothMPC_lbIdx32, smoothMPC_rd32);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi32, smoothMPC_rd32, smoothMPC_Lbyrd32);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V31, smoothMPC_Lbyrd31, smoothMPC_W32, smoothMPC_Lbyrd32, smoothMPC_beta32);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd32, smoothMPC_yy31, smoothMPC_beta32, smoothMPC_bmy32);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld32, smoothMPC_bmy32, smoothMPC_yy32);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub33, smoothMPC_sub33, smoothMPC_ubIdx33, smoothMPC_ccrhsl33, smoothMPC_slb33, smoothMPC_lbIdx33, smoothMPC_rd33);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi33, smoothMPC_rd33, smoothMPC_Lbyrd33);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V32, smoothMPC_Lbyrd32, smoothMPC_W33, smoothMPC_Lbyrd33, smoothMPC_beta33);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd33, smoothMPC_yy32, smoothMPC_beta33, smoothMPC_bmy33);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld33, smoothMPC_bmy33, smoothMPC_yy33);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub34, smoothMPC_sub34, smoothMPC_ubIdx34, smoothMPC_ccrhsl34, smoothMPC_slb34, smoothMPC_lbIdx34, smoothMPC_rd34);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi34, smoothMPC_rd34, smoothMPC_Lbyrd34);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V33, smoothMPC_Lbyrd33, smoothMPC_W34, smoothMPC_Lbyrd34, smoothMPC_beta34);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd34, smoothMPC_yy33, smoothMPC_beta34, smoothMPC_bmy34);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld34, smoothMPC_bmy34, smoothMPC_yy34);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub35, smoothMPC_sub35, smoothMPC_ubIdx35, smoothMPC_ccrhsl35, smoothMPC_slb35, smoothMPC_lbIdx35, smoothMPC_rd35);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi35, smoothMPC_rd35, smoothMPC_Lbyrd35);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V34, smoothMPC_Lbyrd34, smoothMPC_W35, smoothMPC_Lbyrd35, smoothMPC_beta35);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd35, smoothMPC_yy34, smoothMPC_beta35, smoothMPC_bmy35);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld35, smoothMPC_bmy35, smoothMPC_yy35);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub36, smoothMPC_sub36, smoothMPC_ubIdx36, smoothMPC_ccrhsl36, smoothMPC_slb36, smoothMPC_lbIdx36, smoothMPC_rd36);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi36, smoothMPC_rd36, smoothMPC_Lbyrd36);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V35, smoothMPC_Lbyrd35, smoothMPC_W36, smoothMPC_Lbyrd36, smoothMPC_beta36);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd36, smoothMPC_yy35, smoothMPC_beta36, smoothMPC_bmy36);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld36, smoothMPC_bmy36, smoothMPC_yy36);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub37, smoothMPC_sub37, smoothMPC_ubIdx37, smoothMPC_ccrhsl37, smoothMPC_slb37, smoothMPC_lbIdx37, smoothMPC_rd37);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi37, smoothMPC_rd37, smoothMPC_Lbyrd37);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V36, smoothMPC_Lbyrd36, smoothMPC_W37, smoothMPC_Lbyrd37, smoothMPC_beta37);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd37, smoothMPC_yy36, smoothMPC_beta37, smoothMPC_bmy37);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld37, smoothMPC_bmy37, smoothMPC_yy37);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub38, smoothMPC_sub38, smoothMPC_ubIdx38, smoothMPC_ccrhsl38, smoothMPC_slb38, smoothMPC_lbIdx38, smoothMPC_rd38);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi38, smoothMPC_rd38, smoothMPC_Lbyrd38);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V37, smoothMPC_Lbyrd37, smoothMPC_W38, smoothMPC_Lbyrd38, smoothMPC_beta38);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd38, smoothMPC_yy37, smoothMPC_beta38, smoothMPC_bmy38);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld38, smoothMPC_bmy38, smoothMPC_yy38);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub39, smoothMPC_sub39, smoothMPC_ubIdx39, smoothMPC_ccrhsl39, smoothMPC_slb39, smoothMPC_lbIdx39, smoothMPC_rd39);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi39, smoothMPC_rd39, smoothMPC_Lbyrd39);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V38, smoothMPC_Lbyrd38, smoothMPC_W39, smoothMPC_Lbyrd39, smoothMPC_beta39);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd39, smoothMPC_yy38, smoothMPC_beta39, smoothMPC_bmy39);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld39, smoothMPC_bmy39, smoothMPC_yy39);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub40, smoothMPC_sub40, smoothMPC_ubIdx40, smoothMPC_ccrhsl40, smoothMPC_slb40, smoothMPC_lbIdx40, smoothMPC_rd40);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi40, smoothMPC_rd40, smoothMPC_Lbyrd40);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V39, smoothMPC_Lbyrd39, smoothMPC_W40, smoothMPC_Lbyrd40, smoothMPC_beta40);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd40, smoothMPC_yy39, smoothMPC_beta40, smoothMPC_bmy40);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld40, smoothMPC_bmy40, smoothMPC_yy40);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub41, smoothMPC_sub41, smoothMPC_ubIdx41, smoothMPC_ccrhsl41, smoothMPC_slb41, smoothMPC_lbIdx41, smoothMPC_rd41);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi41, smoothMPC_rd41, smoothMPC_Lbyrd41);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V40, smoothMPC_Lbyrd40, smoothMPC_W41, smoothMPC_Lbyrd41, smoothMPC_beta41);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd41, smoothMPC_yy40, smoothMPC_beta41, smoothMPC_bmy41);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld41, smoothMPC_bmy41, smoothMPC_yy41);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub42, smoothMPC_sub42, smoothMPC_ubIdx42, smoothMPC_ccrhsl42, smoothMPC_slb42, smoothMPC_lbIdx42, smoothMPC_rd42);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi42, smoothMPC_rd42, smoothMPC_Lbyrd42);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V41, smoothMPC_Lbyrd41, smoothMPC_W42, smoothMPC_Lbyrd42, smoothMPC_beta42);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd42, smoothMPC_yy41, smoothMPC_beta42, smoothMPC_bmy42);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld42, smoothMPC_bmy42, smoothMPC_yy42);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub43, smoothMPC_sub43, smoothMPC_ubIdx43, smoothMPC_ccrhsl43, smoothMPC_slb43, smoothMPC_lbIdx43, smoothMPC_rd43);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi43, smoothMPC_rd43, smoothMPC_Lbyrd43);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V42, smoothMPC_Lbyrd42, smoothMPC_W43, smoothMPC_Lbyrd43, smoothMPC_beta43);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd43, smoothMPC_yy42, smoothMPC_beta43, smoothMPC_bmy43);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld43, smoothMPC_bmy43, smoothMPC_yy43);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub44, smoothMPC_sub44, smoothMPC_ubIdx44, smoothMPC_ccrhsl44, smoothMPC_slb44, smoothMPC_lbIdx44, smoothMPC_rd44);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi44, smoothMPC_rd44, smoothMPC_Lbyrd44);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V43, smoothMPC_Lbyrd43, smoothMPC_W44, smoothMPC_Lbyrd44, smoothMPC_beta44);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd44, smoothMPC_yy43, smoothMPC_beta44, smoothMPC_bmy44);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld44, smoothMPC_bmy44, smoothMPC_yy44);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub45, smoothMPC_sub45, smoothMPC_ubIdx45, smoothMPC_ccrhsl45, smoothMPC_slb45, smoothMPC_lbIdx45, smoothMPC_rd45);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi45, smoothMPC_rd45, smoothMPC_Lbyrd45);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V44, smoothMPC_Lbyrd44, smoothMPC_W45, smoothMPC_Lbyrd45, smoothMPC_beta45);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd45, smoothMPC_yy44, smoothMPC_beta45, smoothMPC_bmy45);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld45, smoothMPC_bmy45, smoothMPC_yy45);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub46, smoothMPC_sub46, smoothMPC_ubIdx46, smoothMPC_ccrhsl46, smoothMPC_slb46, smoothMPC_lbIdx46, smoothMPC_rd46);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi46, smoothMPC_rd46, smoothMPC_Lbyrd46);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V45, smoothMPC_Lbyrd45, smoothMPC_W46, smoothMPC_Lbyrd46, smoothMPC_beta46);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd46, smoothMPC_yy45, smoothMPC_beta46, smoothMPC_bmy46);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld46, smoothMPC_bmy46, smoothMPC_yy46);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub47, smoothMPC_sub47, smoothMPC_ubIdx47, smoothMPC_ccrhsl47, smoothMPC_slb47, smoothMPC_lbIdx47, smoothMPC_rd47);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi47, smoothMPC_rd47, smoothMPC_Lbyrd47);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V46, smoothMPC_Lbyrd46, smoothMPC_W47, smoothMPC_Lbyrd47, smoothMPC_beta47);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd47, smoothMPC_yy46, smoothMPC_beta47, smoothMPC_bmy47);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld47, smoothMPC_bmy47, smoothMPC_yy47);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub48, smoothMPC_sub48, smoothMPC_ubIdx48, smoothMPC_ccrhsl48, smoothMPC_slb48, smoothMPC_lbIdx48, smoothMPC_rd48);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi48, smoothMPC_rd48, smoothMPC_Lbyrd48);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V47, smoothMPC_Lbyrd47, smoothMPC_W48, smoothMPC_Lbyrd48, smoothMPC_beta48);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd48, smoothMPC_yy47, smoothMPC_beta48, smoothMPC_bmy48);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld48, smoothMPC_bmy48, smoothMPC_yy48);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub49, smoothMPC_sub49, smoothMPC_ubIdx49, smoothMPC_ccrhsl49, smoothMPC_slb49, smoothMPC_lbIdx49, smoothMPC_rd49);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi49, smoothMPC_rd49, smoothMPC_Lbyrd49);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V48, smoothMPC_Lbyrd48, smoothMPC_W49, smoothMPC_Lbyrd49, smoothMPC_beta49);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd49, smoothMPC_yy48, smoothMPC_beta49, smoothMPC_bmy49);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld49, smoothMPC_bmy49, smoothMPC_yy49);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub50, smoothMPC_sub50, smoothMPC_ubIdx50, smoothMPC_ccrhsl50, smoothMPC_slb50, smoothMPC_lbIdx50, smoothMPC_rd50);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi50, smoothMPC_rd50, smoothMPC_Lbyrd50);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V49, smoothMPC_Lbyrd49, smoothMPC_W50, smoothMPC_Lbyrd50, smoothMPC_beta50);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd50, smoothMPC_yy49, smoothMPC_beta50, smoothMPC_bmy50);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld50, smoothMPC_bmy50, smoothMPC_yy50);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub51, smoothMPC_sub51, smoothMPC_ubIdx51, smoothMPC_ccrhsl51, smoothMPC_slb51, smoothMPC_lbIdx51, smoothMPC_rd51);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi51, smoothMPC_rd51, smoothMPC_Lbyrd51);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V50, smoothMPC_Lbyrd50, smoothMPC_W51, smoothMPC_Lbyrd51, smoothMPC_beta51);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd51, smoothMPC_yy50, smoothMPC_beta51, smoothMPC_bmy51);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld51, smoothMPC_bmy51, smoothMPC_yy51);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub52, smoothMPC_sub52, smoothMPC_ubIdx52, smoothMPC_ccrhsl52, smoothMPC_slb52, smoothMPC_lbIdx52, smoothMPC_rd52);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi52, smoothMPC_rd52, smoothMPC_Lbyrd52);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V51, smoothMPC_Lbyrd51, smoothMPC_W52, smoothMPC_Lbyrd52, smoothMPC_beta52);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd52, smoothMPC_yy51, smoothMPC_beta52, smoothMPC_bmy52);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld52, smoothMPC_bmy52, smoothMPC_yy52);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub53, smoothMPC_sub53, smoothMPC_ubIdx53, smoothMPC_ccrhsl53, smoothMPC_slb53, smoothMPC_lbIdx53, smoothMPC_rd53);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi53, smoothMPC_rd53, smoothMPC_Lbyrd53);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V52, smoothMPC_Lbyrd52, smoothMPC_W53, smoothMPC_Lbyrd53, smoothMPC_beta53);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd53, smoothMPC_yy52, smoothMPC_beta53, smoothMPC_bmy53);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld53, smoothMPC_bmy53, smoothMPC_yy53);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub54, smoothMPC_sub54, smoothMPC_ubIdx54, smoothMPC_ccrhsl54, smoothMPC_slb54, smoothMPC_lbIdx54, smoothMPC_rd54);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi54, smoothMPC_rd54, smoothMPC_Lbyrd54);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V53, smoothMPC_Lbyrd53, smoothMPC_W54, smoothMPC_Lbyrd54, smoothMPC_beta54);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd54, smoothMPC_yy53, smoothMPC_beta54, smoothMPC_bmy54);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld54, smoothMPC_bmy54, smoothMPC_yy54);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub55, smoothMPC_sub55, smoothMPC_ubIdx55, smoothMPC_ccrhsl55, smoothMPC_slb55, smoothMPC_lbIdx55, smoothMPC_rd55);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi55, smoothMPC_rd55, smoothMPC_Lbyrd55);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V54, smoothMPC_Lbyrd54, smoothMPC_W55, smoothMPC_Lbyrd55, smoothMPC_beta55);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd55, smoothMPC_yy54, smoothMPC_beta55, smoothMPC_bmy55);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld55, smoothMPC_bmy55, smoothMPC_yy55);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub56, smoothMPC_sub56, smoothMPC_ubIdx56, smoothMPC_ccrhsl56, smoothMPC_slb56, smoothMPC_lbIdx56, smoothMPC_rd56);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi56, smoothMPC_rd56, smoothMPC_Lbyrd56);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V55, smoothMPC_Lbyrd55, smoothMPC_W56, smoothMPC_Lbyrd56, smoothMPC_beta56);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd56, smoothMPC_yy55, smoothMPC_beta56, smoothMPC_bmy56);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld56, smoothMPC_bmy56, smoothMPC_yy56);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub57, smoothMPC_sub57, smoothMPC_ubIdx57, smoothMPC_ccrhsl57, smoothMPC_slb57, smoothMPC_lbIdx57, smoothMPC_rd57);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi57, smoothMPC_rd57, smoothMPC_Lbyrd57);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V56, smoothMPC_Lbyrd56, smoothMPC_W57, smoothMPC_Lbyrd57, smoothMPC_beta57);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd57, smoothMPC_yy56, smoothMPC_beta57, smoothMPC_bmy57);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld57, smoothMPC_bmy57, smoothMPC_yy57);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub58, smoothMPC_sub58, smoothMPC_ubIdx58, smoothMPC_ccrhsl58, smoothMPC_slb58, smoothMPC_lbIdx58, smoothMPC_rd58);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi58, smoothMPC_rd58, smoothMPC_Lbyrd58);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V57, smoothMPC_Lbyrd57, smoothMPC_W58, smoothMPC_Lbyrd58, smoothMPC_beta58);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd58, smoothMPC_yy57, smoothMPC_beta58, smoothMPC_bmy58);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld58, smoothMPC_bmy58, smoothMPC_yy58);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub59, smoothMPC_sub59, smoothMPC_ubIdx59, smoothMPC_ccrhsl59, smoothMPC_slb59, smoothMPC_lbIdx59, smoothMPC_rd59);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi59, smoothMPC_rd59, smoothMPC_Lbyrd59);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V58, smoothMPC_Lbyrd58, smoothMPC_W59, smoothMPC_Lbyrd59, smoothMPC_beta59);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd59, smoothMPC_yy58, smoothMPC_beta59, smoothMPC_bmy59);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld59, smoothMPC_bmy59, smoothMPC_yy59);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub60, smoothMPC_sub60, smoothMPC_ubIdx60, smoothMPC_ccrhsl60, smoothMPC_slb60, smoothMPC_lbIdx60, smoothMPC_rd60);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi60, smoothMPC_rd60, smoothMPC_Lbyrd60);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V59, smoothMPC_Lbyrd59, smoothMPC_W60, smoothMPC_Lbyrd60, smoothMPC_beta60);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd60, smoothMPC_yy59, smoothMPC_beta60, smoothMPC_bmy60);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld60, smoothMPC_bmy60, smoothMPC_yy60);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub61, smoothMPC_sub61, smoothMPC_ubIdx61, smoothMPC_ccrhsl61, smoothMPC_slb61, smoothMPC_lbIdx61, smoothMPC_rd61);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi61, smoothMPC_rd61, smoothMPC_Lbyrd61);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V60, smoothMPC_Lbyrd60, smoothMPC_W61, smoothMPC_Lbyrd61, smoothMPC_beta61);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd61, smoothMPC_yy60, smoothMPC_beta61, smoothMPC_bmy61);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld61, smoothMPC_bmy61, smoothMPC_yy61);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub62, smoothMPC_sub62, smoothMPC_ubIdx62, smoothMPC_ccrhsl62, smoothMPC_slb62, smoothMPC_lbIdx62, smoothMPC_rd62);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi62, smoothMPC_rd62, smoothMPC_Lbyrd62);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V61, smoothMPC_Lbyrd61, smoothMPC_W62, smoothMPC_Lbyrd62, smoothMPC_beta62);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd62, smoothMPC_yy61, smoothMPC_beta62, smoothMPC_bmy62);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld62, smoothMPC_bmy62, smoothMPC_yy62);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub63, smoothMPC_sub63, smoothMPC_ubIdx63, smoothMPC_ccrhsl63, smoothMPC_slb63, smoothMPC_lbIdx63, smoothMPC_rd63);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi63, smoothMPC_rd63, smoothMPC_Lbyrd63);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V62, smoothMPC_Lbyrd62, smoothMPC_W63, smoothMPC_Lbyrd63, smoothMPC_beta63);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd63, smoothMPC_yy62, smoothMPC_beta63, smoothMPC_bmy63);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld63, smoothMPC_bmy63, smoothMPC_yy63);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub64, smoothMPC_sub64, smoothMPC_ubIdx64, smoothMPC_ccrhsl64, smoothMPC_slb64, smoothMPC_lbIdx64, smoothMPC_rd64);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi64, smoothMPC_rd64, smoothMPC_Lbyrd64);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V63, smoothMPC_Lbyrd63, smoothMPC_W64, smoothMPC_Lbyrd64, smoothMPC_beta64);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd64, smoothMPC_yy63, smoothMPC_beta64, smoothMPC_bmy64);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld64, smoothMPC_bmy64, smoothMPC_yy64);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub65, smoothMPC_sub65, smoothMPC_ubIdx65, smoothMPC_ccrhsl65, smoothMPC_slb65, smoothMPC_lbIdx65, smoothMPC_rd65);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi65, smoothMPC_rd65, smoothMPC_Lbyrd65);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V64, smoothMPC_Lbyrd64, smoothMPC_W65, smoothMPC_Lbyrd65, smoothMPC_beta65);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd65, smoothMPC_yy64, smoothMPC_beta65, smoothMPC_bmy65);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld65, smoothMPC_bmy65, smoothMPC_yy65);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub66, smoothMPC_sub66, smoothMPC_ubIdx66, smoothMPC_ccrhsl66, smoothMPC_slb66, smoothMPC_lbIdx66, smoothMPC_rd66);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi66, smoothMPC_rd66, smoothMPC_Lbyrd66);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V65, smoothMPC_Lbyrd65, smoothMPC_W66, smoothMPC_Lbyrd66, smoothMPC_beta66);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd66, smoothMPC_yy65, smoothMPC_beta66, smoothMPC_bmy66);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld66, smoothMPC_bmy66, smoothMPC_yy66);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub67, smoothMPC_sub67, smoothMPC_ubIdx67, smoothMPC_ccrhsl67, smoothMPC_slb67, smoothMPC_lbIdx67, smoothMPC_rd67);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi67, smoothMPC_rd67, smoothMPC_Lbyrd67);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V66, smoothMPC_Lbyrd66, smoothMPC_W67, smoothMPC_Lbyrd67, smoothMPC_beta67);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd67, smoothMPC_yy66, smoothMPC_beta67, smoothMPC_bmy67);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld67, smoothMPC_bmy67, smoothMPC_yy67);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub68, smoothMPC_sub68, smoothMPC_ubIdx68, smoothMPC_ccrhsl68, smoothMPC_slb68, smoothMPC_lbIdx68, smoothMPC_rd68);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi68, smoothMPC_rd68, smoothMPC_Lbyrd68);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V67, smoothMPC_Lbyrd67, smoothMPC_W68, smoothMPC_Lbyrd68, smoothMPC_beta68);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd68, smoothMPC_yy67, smoothMPC_beta68, smoothMPC_bmy68);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld68, smoothMPC_bmy68, smoothMPC_yy68);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub69, smoothMPC_sub69, smoothMPC_ubIdx69, smoothMPC_ccrhsl69, smoothMPC_slb69, smoothMPC_lbIdx69, smoothMPC_rd69);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi69, smoothMPC_rd69, smoothMPC_Lbyrd69);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V68, smoothMPC_Lbyrd68, smoothMPC_W69, smoothMPC_Lbyrd69, smoothMPC_beta69);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd69, smoothMPC_yy68, smoothMPC_beta69, smoothMPC_bmy69);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld69, smoothMPC_bmy69, smoothMPC_yy69);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub70, smoothMPC_sub70, smoothMPC_ubIdx70, smoothMPC_ccrhsl70, smoothMPC_slb70, smoothMPC_lbIdx70, smoothMPC_rd70);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi70, smoothMPC_rd70, smoothMPC_Lbyrd70);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V69, smoothMPC_Lbyrd69, smoothMPC_W70, smoothMPC_Lbyrd70, smoothMPC_beta70);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd70, smoothMPC_yy69, smoothMPC_beta70, smoothMPC_bmy70);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld70, smoothMPC_bmy70, smoothMPC_yy70);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub71, smoothMPC_sub71, smoothMPC_ubIdx71, smoothMPC_ccrhsl71, smoothMPC_slb71, smoothMPC_lbIdx71, smoothMPC_rd71);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi71, smoothMPC_rd71, smoothMPC_Lbyrd71);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V70, smoothMPC_Lbyrd70, smoothMPC_W71, smoothMPC_Lbyrd71, smoothMPC_beta71);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd71, smoothMPC_yy70, smoothMPC_beta71, smoothMPC_bmy71);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld71, smoothMPC_bmy71, smoothMPC_yy71);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub72, smoothMPC_sub72, smoothMPC_ubIdx72, smoothMPC_ccrhsl72, smoothMPC_slb72, smoothMPC_lbIdx72, smoothMPC_rd72);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi72, smoothMPC_rd72, smoothMPC_Lbyrd72);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V71, smoothMPC_Lbyrd71, smoothMPC_W72, smoothMPC_Lbyrd72, smoothMPC_beta72);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd72, smoothMPC_yy71, smoothMPC_beta72, smoothMPC_bmy72);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld72, smoothMPC_bmy72, smoothMPC_yy72);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub73, smoothMPC_sub73, smoothMPC_ubIdx73, smoothMPC_ccrhsl73, smoothMPC_slb73, smoothMPC_lbIdx73, smoothMPC_rd73);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi73, smoothMPC_rd73, smoothMPC_Lbyrd73);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V72, smoothMPC_Lbyrd72, smoothMPC_W73, smoothMPC_Lbyrd73, smoothMPC_beta73);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd73, smoothMPC_yy72, smoothMPC_beta73, smoothMPC_bmy73);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld73, smoothMPC_bmy73, smoothMPC_yy73);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub74, smoothMPC_sub74, smoothMPC_ubIdx74, smoothMPC_ccrhsl74, smoothMPC_slb74, smoothMPC_lbIdx74, smoothMPC_rd74);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi74, smoothMPC_rd74, smoothMPC_Lbyrd74);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V73, smoothMPC_Lbyrd73, smoothMPC_W74, smoothMPC_Lbyrd74, smoothMPC_beta74);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd74, smoothMPC_yy73, smoothMPC_beta74, smoothMPC_bmy74);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld74, smoothMPC_bmy74, smoothMPC_yy74);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub75, smoothMPC_sub75, smoothMPC_ubIdx75, smoothMPC_ccrhsl75, smoothMPC_slb75, smoothMPC_lbIdx75, smoothMPC_rd75);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi75, smoothMPC_rd75, smoothMPC_Lbyrd75);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V74, smoothMPC_Lbyrd74, smoothMPC_W75, smoothMPC_Lbyrd75, smoothMPC_beta75);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd75, smoothMPC_yy74, smoothMPC_beta75, smoothMPC_bmy75);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld75, smoothMPC_bmy75, smoothMPC_yy75);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub76, smoothMPC_sub76, smoothMPC_ubIdx76, smoothMPC_ccrhsl76, smoothMPC_slb76, smoothMPC_lbIdx76, smoothMPC_rd76);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi76, smoothMPC_rd76, smoothMPC_Lbyrd76);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V75, smoothMPC_Lbyrd75, smoothMPC_W76, smoothMPC_Lbyrd76, smoothMPC_beta76);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd76, smoothMPC_yy75, smoothMPC_beta76, smoothMPC_bmy76);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld76, smoothMPC_bmy76, smoothMPC_yy76);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub77, smoothMPC_sub77, smoothMPC_ubIdx77, smoothMPC_ccrhsl77, smoothMPC_slb77, smoothMPC_lbIdx77, smoothMPC_rd77);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi77, smoothMPC_rd77, smoothMPC_Lbyrd77);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V76, smoothMPC_Lbyrd76, smoothMPC_W77, smoothMPC_Lbyrd77, smoothMPC_beta77);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd77, smoothMPC_yy76, smoothMPC_beta77, smoothMPC_bmy77);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld77, smoothMPC_bmy77, smoothMPC_yy77);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub78, smoothMPC_sub78, smoothMPC_ubIdx78, smoothMPC_ccrhsl78, smoothMPC_slb78, smoothMPC_lbIdx78, smoothMPC_rd78);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi78, smoothMPC_rd78, smoothMPC_Lbyrd78);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V77, smoothMPC_Lbyrd77, smoothMPC_W78, smoothMPC_Lbyrd78, smoothMPC_beta78);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd78, smoothMPC_yy77, smoothMPC_beta78, smoothMPC_bmy78);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld78, smoothMPC_bmy78, smoothMPC_yy78);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub79, smoothMPC_sub79, smoothMPC_ubIdx79, smoothMPC_ccrhsl79, smoothMPC_slb79, smoothMPC_lbIdx79, smoothMPC_rd79);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi79, smoothMPC_rd79, smoothMPC_Lbyrd79);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V78, smoothMPC_Lbyrd78, smoothMPC_W79, smoothMPC_Lbyrd79, smoothMPC_beta79);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd79, smoothMPC_yy78, smoothMPC_beta79, smoothMPC_bmy79);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld79, smoothMPC_bmy79, smoothMPC_yy79);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub80, smoothMPC_sub80, smoothMPC_ubIdx80, smoothMPC_ccrhsl80, smoothMPC_slb80, smoothMPC_lbIdx80, smoothMPC_rd80);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi80, smoothMPC_rd80, smoothMPC_Lbyrd80);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V79, smoothMPC_Lbyrd79, smoothMPC_W80, smoothMPC_Lbyrd80, smoothMPC_beta80);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd80, smoothMPC_yy79, smoothMPC_beta80, smoothMPC_bmy80);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld80, smoothMPC_bmy80, smoothMPC_yy80);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub81, smoothMPC_sub81, smoothMPC_ubIdx81, smoothMPC_ccrhsl81, smoothMPC_slb81, smoothMPC_lbIdx81, smoothMPC_rd81);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi81, smoothMPC_rd81, smoothMPC_Lbyrd81);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V80, smoothMPC_Lbyrd80, smoothMPC_W81, smoothMPC_Lbyrd81, smoothMPC_beta81);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd81, smoothMPC_yy80, smoothMPC_beta81, smoothMPC_bmy81);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld81, smoothMPC_bmy81, smoothMPC_yy81);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub82, smoothMPC_sub82, smoothMPC_ubIdx82, smoothMPC_ccrhsl82, smoothMPC_slb82, smoothMPC_lbIdx82, smoothMPC_rd82);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi82, smoothMPC_rd82, smoothMPC_Lbyrd82);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V81, smoothMPC_Lbyrd81, smoothMPC_W82, smoothMPC_Lbyrd82, smoothMPC_beta82);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd82, smoothMPC_yy81, smoothMPC_beta82, smoothMPC_bmy82);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld82, smoothMPC_bmy82, smoothMPC_yy82);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub83, smoothMPC_sub83, smoothMPC_ubIdx83, smoothMPC_ccrhsl83, smoothMPC_slb83, smoothMPC_lbIdx83, smoothMPC_rd83);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi83, smoothMPC_rd83, smoothMPC_Lbyrd83);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V82, smoothMPC_Lbyrd82, smoothMPC_W83, smoothMPC_Lbyrd83, smoothMPC_beta83);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd83, smoothMPC_yy82, smoothMPC_beta83, smoothMPC_bmy83);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld83, smoothMPC_bmy83, smoothMPC_yy83);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub84, smoothMPC_sub84, smoothMPC_ubIdx84, smoothMPC_ccrhsl84, smoothMPC_slb84, smoothMPC_lbIdx84, smoothMPC_rd84);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi84, smoothMPC_rd84, smoothMPC_Lbyrd84);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V83, smoothMPC_Lbyrd83, smoothMPC_W84, smoothMPC_Lbyrd84, smoothMPC_beta84);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd84, smoothMPC_yy83, smoothMPC_beta84, smoothMPC_bmy84);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld84, smoothMPC_bmy84, smoothMPC_yy84);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub85, smoothMPC_sub85, smoothMPC_ubIdx85, smoothMPC_ccrhsl85, smoothMPC_slb85, smoothMPC_lbIdx85, smoothMPC_rd85);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi85, smoothMPC_rd85, smoothMPC_Lbyrd85);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V84, smoothMPC_Lbyrd84, smoothMPC_W85, smoothMPC_Lbyrd85, smoothMPC_beta85);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd85, smoothMPC_yy84, smoothMPC_beta85, smoothMPC_bmy85);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld85, smoothMPC_bmy85, smoothMPC_yy85);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub86, smoothMPC_sub86, smoothMPC_ubIdx86, smoothMPC_ccrhsl86, smoothMPC_slb86, smoothMPC_lbIdx86, smoothMPC_rd86);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi86, smoothMPC_rd86, smoothMPC_Lbyrd86);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V85, smoothMPC_Lbyrd85, smoothMPC_W86, smoothMPC_Lbyrd86, smoothMPC_beta86);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd86, smoothMPC_yy85, smoothMPC_beta86, smoothMPC_bmy86);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld86, smoothMPC_bmy86, smoothMPC_yy86);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub87, smoothMPC_sub87, smoothMPC_ubIdx87, smoothMPC_ccrhsl87, smoothMPC_slb87, smoothMPC_lbIdx87, smoothMPC_rd87);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi87, smoothMPC_rd87, smoothMPC_Lbyrd87);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V86, smoothMPC_Lbyrd86, smoothMPC_W87, smoothMPC_Lbyrd87, smoothMPC_beta87);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd87, smoothMPC_yy86, smoothMPC_beta87, smoothMPC_bmy87);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld87, smoothMPC_bmy87, smoothMPC_yy87);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub88, smoothMPC_sub88, smoothMPC_ubIdx88, smoothMPC_ccrhsl88, smoothMPC_slb88, smoothMPC_lbIdx88, smoothMPC_rd88);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi88, smoothMPC_rd88, smoothMPC_Lbyrd88);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V87, smoothMPC_Lbyrd87, smoothMPC_W88, smoothMPC_Lbyrd88, smoothMPC_beta88);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd88, smoothMPC_yy87, smoothMPC_beta88, smoothMPC_bmy88);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld88, smoothMPC_bmy88, smoothMPC_yy88);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub89, smoothMPC_sub89, smoothMPC_ubIdx89, smoothMPC_ccrhsl89, smoothMPC_slb89, smoothMPC_lbIdx89, smoothMPC_rd89);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi89, smoothMPC_rd89, smoothMPC_Lbyrd89);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V88, smoothMPC_Lbyrd88, smoothMPC_W89, smoothMPC_Lbyrd89, smoothMPC_beta89);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd89, smoothMPC_yy88, smoothMPC_beta89, smoothMPC_bmy89);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld89, smoothMPC_bmy89, smoothMPC_yy89);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub90, smoothMPC_sub90, smoothMPC_ubIdx90, smoothMPC_ccrhsl90, smoothMPC_slb90, smoothMPC_lbIdx90, smoothMPC_rd90);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi90, smoothMPC_rd90, smoothMPC_Lbyrd90);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V89, smoothMPC_Lbyrd89, smoothMPC_W90, smoothMPC_Lbyrd90, smoothMPC_beta90);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd90, smoothMPC_yy89, smoothMPC_beta90, smoothMPC_bmy90);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld90, smoothMPC_bmy90, smoothMPC_yy90);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub91, smoothMPC_sub91, smoothMPC_ubIdx91, smoothMPC_ccrhsl91, smoothMPC_slb91, smoothMPC_lbIdx91, smoothMPC_rd91);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi91, smoothMPC_rd91, smoothMPC_Lbyrd91);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V90, smoothMPC_Lbyrd90, smoothMPC_W91, smoothMPC_Lbyrd91, smoothMPC_beta91);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd91, smoothMPC_yy90, smoothMPC_beta91, smoothMPC_bmy91);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld91, smoothMPC_bmy91, smoothMPC_yy91);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub92, smoothMPC_sub92, smoothMPC_ubIdx92, smoothMPC_ccrhsl92, smoothMPC_slb92, smoothMPC_lbIdx92, smoothMPC_rd92);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi92, smoothMPC_rd92, smoothMPC_Lbyrd92);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V91, smoothMPC_Lbyrd91, smoothMPC_W92, smoothMPC_Lbyrd92, smoothMPC_beta92);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd92, smoothMPC_yy91, smoothMPC_beta92, smoothMPC_bmy92);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld92, smoothMPC_bmy92, smoothMPC_yy92);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub93, smoothMPC_sub93, smoothMPC_ubIdx93, smoothMPC_ccrhsl93, smoothMPC_slb93, smoothMPC_lbIdx93, smoothMPC_rd93);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi93, smoothMPC_rd93, smoothMPC_Lbyrd93);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V92, smoothMPC_Lbyrd92, smoothMPC_W93, smoothMPC_Lbyrd93, smoothMPC_beta93);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd93, smoothMPC_yy92, smoothMPC_beta93, smoothMPC_bmy93);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld93, smoothMPC_bmy93, smoothMPC_yy93);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub94, smoothMPC_sub94, smoothMPC_ubIdx94, smoothMPC_ccrhsl94, smoothMPC_slb94, smoothMPC_lbIdx94, smoothMPC_rd94);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi94, smoothMPC_rd94, smoothMPC_Lbyrd94);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V93, smoothMPC_Lbyrd93, smoothMPC_W94, smoothMPC_Lbyrd94, smoothMPC_beta94);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd94, smoothMPC_yy93, smoothMPC_beta94, smoothMPC_bmy94);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld94, smoothMPC_bmy94, smoothMPC_yy94);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub95, smoothMPC_sub95, smoothMPC_ubIdx95, smoothMPC_ccrhsl95, smoothMPC_slb95, smoothMPC_lbIdx95, smoothMPC_rd95);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi95, smoothMPC_rd95, smoothMPC_Lbyrd95);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V94, smoothMPC_Lbyrd94, smoothMPC_W95, smoothMPC_Lbyrd95, smoothMPC_beta95);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd95, smoothMPC_yy94, smoothMPC_beta95, smoothMPC_bmy95);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld95, smoothMPC_bmy95, smoothMPC_yy95);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub96, smoothMPC_sub96, smoothMPC_ubIdx96, smoothMPC_ccrhsl96, smoothMPC_slb96, smoothMPC_lbIdx96, smoothMPC_rd96);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi96, smoothMPC_rd96, smoothMPC_Lbyrd96);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V95, smoothMPC_Lbyrd95, smoothMPC_W96, smoothMPC_Lbyrd96, smoothMPC_beta96);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd96, smoothMPC_yy95, smoothMPC_beta96, smoothMPC_bmy96);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld96, smoothMPC_bmy96, smoothMPC_yy96);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub97, smoothMPC_sub97, smoothMPC_ubIdx97, smoothMPC_ccrhsl97, smoothMPC_slb97, smoothMPC_lbIdx97, smoothMPC_rd97);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi97, smoothMPC_rd97, smoothMPC_Lbyrd97);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V96, smoothMPC_Lbyrd96, smoothMPC_W97, smoothMPC_Lbyrd97, smoothMPC_beta97);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd97, smoothMPC_yy96, smoothMPC_beta97, smoothMPC_bmy97);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld97, smoothMPC_bmy97, smoothMPC_yy97);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub98, smoothMPC_sub98, smoothMPC_ubIdx98, smoothMPC_ccrhsl98, smoothMPC_slb98, smoothMPC_lbIdx98, smoothMPC_rd98);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi98, smoothMPC_rd98, smoothMPC_Lbyrd98);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V97, smoothMPC_Lbyrd97, smoothMPC_W98, smoothMPC_Lbyrd98, smoothMPC_beta98);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd98, smoothMPC_yy97, smoothMPC_beta98, smoothMPC_bmy98);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld98, smoothMPC_bmy98, smoothMPC_yy98);
smoothMPC_LA_VSUB6_INDEXED_3_3_3(smoothMPC_ccrhsub99, smoothMPC_sub99, smoothMPC_ubIdx99, smoothMPC_ccrhsl99, smoothMPC_slb99, smoothMPC_lbIdx99, smoothMPC_rd99);
smoothMPC_LA_DIAG_FORWARDSUB_3(smoothMPC_Phi99, smoothMPC_rd99, smoothMPC_Lbyrd99);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_3(smoothMPC_V98, smoothMPC_Lbyrd98, smoothMPC_W99, smoothMPC_Lbyrd99, smoothMPC_beta99);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd99, smoothMPC_yy98, smoothMPC_beta99, smoothMPC_bmy99);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld99, smoothMPC_bmy99, smoothMPC_yy99);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld99, smoothMPC_yy99, smoothMPC_dvcc99);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd99, smoothMPC_dvcc99, smoothMPC_yy98, smoothMPC_bmy98);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld98, smoothMPC_bmy98, smoothMPC_dvcc98);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd98, smoothMPC_dvcc98, smoothMPC_yy97, smoothMPC_bmy97);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld97, smoothMPC_bmy97, smoothMPC_dvcc97);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd97, smoothMPC_dvcc97, smoothMPC_yy96, smoothMPC_bmy96);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld96, smoothMPC_bmy96, smoothMPC_dvcc96);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd96, smoothMPC_dvcc96, smoothMPC_yy95, smoothMPC_bmy95);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld95, smoothMPC_bmy95, smoothMPC_dvcc95);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd95, smoothMPC_dvcc95, smoothMPC_yy94, smoothMPC_bmy94);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld94, smoothMPC_bmy94, smoothMPC_dvcc94);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd94, smoothMPC_dvcc94, smoothMPC_yy93, smoothMPC_bmy93);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld93, smoothMPC_bmy93, smoothMPC_dvcc93);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd93, smoothMPC_dvcc93, smoothMPC_yy92, smoothMPC_bmy92);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld92, smoothMPC_bmy92, smoothMPC_dvcc92);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd92, smoothMPC_dvcc92, smoothMPC_yy91, smoothMPC_bmy91);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld91, smoothMPC_bmy91, smoothMPC_dvcc91);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd91, smoothMPC_dvcc91, smoothMPC_yy90, smoothMPC_bmy90);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld90, smoothMPC_bmy90, smoothMPC_dvcc90);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd90, smoothMPC_dvcc90, smoothMPC_yy89, smoothMPC_bmy89);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld89, smoothMPC_bmy89, smoothMPC_dvcc89);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd89, smoothMPC_dvcc89, smoothMPC_yy88, smoothMPC_bmy88);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld88, smoothMPC_bmy88, smoothMPC_dvcc88);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd88, smoothMPC_dvcc88, smoothMPC_yy87, smoothMPC_bmy87);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld87, smoothMPC_bmy87, smoothMPC_dvcc87);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd87, smoothMPC_dvcc87, smoothMPC_yy86, smoothMPC_bmy86);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld86, smoothMPC_bmy86, smoothMPC_dvcc86);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd86, smoothMPC_dvcc86, smoothMPC_yy85, smoothMPC_bmy85);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld85, smoothMPC_bmy85, smoothMPC_dvcc85);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd85, smoothMPC_dvcc85, smoothMPC_yy84, smoothMPC_bmy84);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld84, smoothMPC_bmy84, smoothMPC_dvcc84);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd84, smoothMPC_dvcc84, smoothMPC_yy83, smoothMPC_bmy83);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld83, smoothMPC_bmy83, smoothMPC_dvcc83);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd83, smoothMPC_dvcc83, smoothMPC_yy82, smoothMPC_bmy82);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld82, smoothMPC_bmy82, smoothMPC_dvcc82);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd82, smoothMPC_dvcc82, smoothMPC_yy81, smoothMPC_bmy81);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld81, smoothMPC_bmy81, smoothMPC_dvcc81);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd81, smoothMPC_dvcc81, smoothMPC_yy80, smoothMPC_bmy80);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld80, smoothMPC_bmy80, smoothMPC_dvcc80);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd80, smoothMPC_dvcc80, smoothMPC_yy79, smoothMPC_bmy79);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld79, smoothMPC_bmy79, smoothMPC_dvcc79);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd79, smoothMPC_dvcc79, smoothMPC_yy78, smoothMPC_bmy78);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld78, smoothMPC_bmy78, smoothMPC_dvcc78);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd78, smoothMPC_dvcc78, smoothMPC_yy77, smoothMPC_bmy77);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld77, smoothMPC_bmy77, smoothMPC_dvcc77);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd77, smoothMPC_dvcc77, smoothMPC_yy76, smoothMPC_bmy76);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld76, smoothMPC_bmy76, smoothMPC_dvcc76);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd76, smoothMPC_dvcc76, smoothMPC_yy75, smoothMPC_bmy75);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld75, smoothMPC_bmy75, smoothMPC_dvcc75);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd75, smoothMPC_dvcc75, smoothMPC_yy74, smoothMPC_bmy74);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld74, smoothMPC_bmy74, smoothMPC_dvcc74);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd74, smoothMPC_dvcc74, smoothMPC_yy73, smoothMPC_bmy73);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld73, smoothMPC_bmy73, smoothMPC_dvcc73);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd73, smoothMPC_dvcc73, smoothMPC_yy72, smoothMPC_bmy72);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld72, smoothMPC_bmy72, smoothMPC_dvcc72);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd72, smoothMPC_dvcc72, smoothMPC_yy71, smoothMPC_bmy71);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld71, smoothMPC_bmy71, smoothMPC_dvcc71);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd71, smoothMPC_dvcc71, smoothMPC_yy70, smoothMPC_bmy70);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld70, smoothMPC_bmy70, smoothMPC_dvcc70);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd70, smoothMPC_dvcc70, smoothMPC_yy69, smoothMPC_bmy69);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld69, smoothMPC_bmy69, smoothMPC_dvcc69);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd69, smoothMPC_dvcc69, smoothMPC_yy68, smoothMPC_bmy68);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld68, smoothMPC_bmy68, smoothMPC_dvcc68);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd68, smoothMPC_dvcc68, smoothMPC_yy67, smoothMPC_bmy67);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld67, smoothMPC_bmy67, smoothMPC_dvcc67);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd67, smoothMPC_dvcc67, smoothMPC_yy66, smoothMPC_bmy66);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld66, smoothMPC_bmy66, smoothMPC_dvcc66);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd66, smoothMPC_dvcc66, smoothMPC_yy65, smoothMPC_bmy65);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld65, smoothMPC_bmy65, smoothMPC_dvcc65);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd65, smoothMPC_dvcc65, smoothMPC_yy64, smoothMPC_bmy64);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld64, smoothMPC_bmy64, smoothMPC_dvcc64);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd64, smoothMPC_dvcc64, smoothMPC_yy63, smoothMPC_bmy63);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld63, smoothMPC_bmy63, smoothMPC_dvcc63);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd63, smoothMPC_dvcc63, smoothMPC_yy62, smoothMPC_bmy62);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld62, smoothMPC_bmy62, smoothMPC_dvcc62);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd62, smoothMPC_dvcc62, smoothMPC_yy61, smoothMPC_bmy61);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld61, smoothMPC_bmy61, smoothMPC_dvcc61);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd61, smoothMPC_dvcc61, smoothMPC_yy60, smoothMPC_bmy60);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld60, smoothMPC_bmy60, smoothMPC_dvcc60);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd60, smoothMPC_dvcc60, smoothMPC_yy59, smoothMPC_bmy59);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld59, smoothMPC_bmy59, smoothMPC_dvcc59);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd59, smoothMPC_dvcc59, smoothMPC_yy58, smoothMPC_bmy58);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld58, smoothMPC_bmy58, smoothMPC_dvcc58);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd58, smoothMPC_dvcc58, smoothMPC_yy57, smoothMPC_bmy57);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld57, smoothMPC_bmy57, smoothMPC_dvcc57);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd57, smoothMPC_dvcc57, smoothMPC_yy56, smoothMPC_bmy56);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld56, smoothMPC_bmy56, smoothMPC_dvcc56);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd56, smoothMPC_dvcc56, smoothMPC_yy55, smoothMPC_bmy55);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld55, smoothMPC_bmy55, smoothMPC_dvcc55);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd55, smoothMPC_dvcc55, smoothMPC_yy54, smoothMPC_bmy54);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld54, smoothMPC_bmy54, smoothMPC_dvcc54);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd54, smoothMPC_dvcc54, smoothMPC_yy53, smoothMPC_bmy53);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld53, smoothMPC_bmy53, smoothMPC_dvcc53);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd53, smoothMPC_dvcc53, smoothMPC_yy52, smoothMPC_bmy52);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld52, smoothMPC_bmy52, smoothMPC_dvcc52);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd52, smoothMPC_dvcc52, smoothMPC_yy51, smoothMPC_bmy51);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld51, smoothMPC_bmy51, smoothMPC_dvcc51);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd51, smoothMPC_dvcc51, smoothMPC_yy50, smoothMPC_bmy50);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld50, smoothMPC_bmy50, smoothMPC_dvcc50);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd50, smoothMPC_dvcc50, smoothMPC_yy49, smoothMPC_bmy49);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld49, smoothMPC_bmy49, smoothMPC_dvcc49);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd49, smoothMPC_dvcc49, smoothMPC_yy48, smoothMPC_bmy48);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld48, smoothMPC_bmy48, smoothMPC_dvcc48);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd48, smoothMPC_dvcc48, smoothMPC_yy47, smoothMPC_bmy47);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld47, smoothMPC_bmy47, smoothMPC_dvcc47);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd47, smoothMPC_dvcc47, smoothMPC_yy46, smoothMPC_bmy46);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld46, smoothMPC_bmy46, smoothMPC_dvcc46);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd46, smoothMPC_dvcc46, smoothMPC_yy45, smoothMPC_bmy45);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld45, smoothMPC_bmy45, smoothMPC_dvcc45);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd45, smoothMPC_dvcc45, smoothMPC_yy44, smoothMPC_bmy44);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld44, smoothMPC_bmy44, smoothMPC_dvcc44);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd44, smoothMPC_dvcc44, smoothMPC_yy43, smoothMPC_bmy43);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld43, smoothMPC_bmy43, smoothMPC_dvcc43);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd43, smoothMPC_dvcc43, smoothMPC_yy42, smoothMPC_bmy42);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld42, smoothMPC_bmy42, smoothMPC_dvcc42);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd42, smoothMPC_dvcc42, smoothMPC_yy41, smoothMPC_bmy41);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld41, smoothMPC_bmy41, smoothMPC_dvcc41);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd41, smoothMPC_dvcc41, smoothMPC_yy40, smoothMPC_bmy40);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld40, smoothMPC_bmy40, smoothMPC_dvcc40);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd40, smoothMPC_dvcc40, smoothMPC_yy39, smoothMPC_bmy39);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld39, smoothMPC_bmy39, smoothMPC_dvcc39);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd39, smoothMPC_dvcc39, smoothMPC_yy38, smoothMPC_bmy38);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld38, smoothMPC_bmy38, smoothMPC_dvcc38);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd38, smoothMPC_dvcc38, smoothMPC_yy37, smoothMPC_bmy37);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld37, smoothMPC_bmy37, smoothMPC_dvcc37);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd37, smoothMPC_dvcc37, smoothMPC_yy36, smoothMPC_bmy36);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld36, smoothMPC_bmy36, smoothMPC_dvcc36);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd36, smoothMPC_dvcc36, smoothMPC_yy35, smoothMPC_bmy35);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld35, smoothMPC_bmy35, smoothMPC_dvcc35);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd35, smoothMPC_dvcc35, smoothMPC_yy34, smoothMPC_bmy34);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld34, smoothMPC_bmy34, smoothMPC_dvcc34);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd34, smoothMPC_dvcc34, smoothMPC_yy33, smoothMPC_bmy33);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld33, smoothMPC_bmy33, smoothMPC_dvcc33);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd33, smoothMPC_dvcc33, smoothMPC_yy32, smoothMPC_bmy32);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld32, smoothMPC_bmy32, smoothMPC_dvcc32);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd32, smoothMPC_dvcc32, smoothMPC_yy31, smoothMPC_bmy31);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld31, smoothMPC_bmy31, smoothMPC_dvcc31);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd31, smoothMPC_dvcc31, smoothMPC_yy30, smoothMPC_bmy30);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld30, smoothMPC_bmy30, smoothMPC_dvcc30);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd30, smoothMPC_dvcc30, smoothMPC_yy29, smoothMPC_bmy29);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld29, smoothMPC_bmy29, smoothMPC_dvcc29);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd29, smoothMPC_dvcc29, smoothMPC_yy28, smoothMPC_bmy28);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld28, smoothMPC_bmy28, smoothMPC_dvcc28);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd28, smoothMPC_dvcc28, smoothMPC_yy27, smoothMPC_bmy27);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld27, smoothMPC_bmy27, smoothMPC_dvcc27);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd27, smoothMPC_dvcc27, smoothMPC_yy26, smoothMPC_bmy26);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld26, smoothMPC_bmy26, smoothMPC_dvcc26);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd26, smoothMPC_dvcc26, smoothMPC_yy25, smoothMPC_bmy25);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld25, smoothMPC_bmy25, smoothMPC_dvcc25);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd25, smoothMPC_dvcc25, smoothMPC_yy24, smoothMPC_bmy24);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld24, smoothMPC_bmy24, smoothMPC_dvcc24);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd24, smoothMPC_dvcc24, smoothMPC_yy23, smoothMPC_bmy23);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld23, smoothMPC_bmy23, smoothMPC_dvcc23);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd23, smoothMPC_dvcc23, smoothMPC_yy22, smoothMPC_bmy22);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld22, smoothMPC_bmy22, smoothMPC_dvcc22);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd22, smoothMPC_dvcc22, smoothMPC_yy21, smoothMPC_bmy21);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld21, smoothMPC_bmy21, smoothMPC_dvcc21);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd21, smoothMPC_dvcc21, smoothMPC_yy20, smoothMPC_bmy20);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld20, smoothMPC_bmy20, smoothMPC_dvcc20);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd20, smoothMPC_dvcc20, smoothMPC_yy19, smoothMPC_bmy19);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld19, smoothMPC_bmy19, smoothMPC_dvcc19);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd19, smoothMPC_dvcc19, smoothMPC_yy18, smoothMPC_bmy18);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld18, smoothMPC_bmy18, smoothMPC_dvcc18);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd18, smoothMPC_dvcc18, smoothMPC_yy17, smoothMPC_bmy17);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld17, smoothMPC_bmy17, smoothMPC_dvcc17);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd17, smoothMPC_dvcc17, smoothMPC_yy16, smoothMPC_bmy16);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld16, smoothMPC_bmy16, smoothMPC_dvcc16);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd16, smoothMPC_dvcc16, smoothMPC_yy15, smoothMPC_bmy15);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld15, smoothMPC_bmy15, smoothMPC_dvcc15);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd15, smoothMPC_dvcc15, smoothMPC_yy14, smoothMPC_bmy14);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld14, smoothMPC_bmy14, smoothMPC_dvcc14);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd14, smoothMPC_dvcc14, smoothMPC_yy13, smoothMPC_bmy13);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld13, smoothMPC_bmy13, smoothMPC_dvcc13);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd13, smoothMPC_dvcc13, smoothMPC_yy12, smoothMPC_bmy12);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld12, smoothMPC_bmy12, smoothMPC_dvcc12);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd12, smoothMPC_dvcc12, smoothMPC_yy11, smoothMPC_bmy11);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld11, smoothMPC_bmy11, smoothMPC_dvcc11);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd11, smoothMPC_dvcc11, smoothMPC_yy10, smoothMPC_bmy10);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld10, smoothMPC_bmy10, smoothMPC_dvcc10);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd10, smoothMPC_dvcc10, smoothMPC_yy09, smoothMPC_bmy09);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld09, smoothMPC_bmy09, smoothMPC_dvcc09);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd09, smoothMPC_dvcc09, smoothMPC_yy08, smoothMPC_bmy08);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld08, smoothMPC_bmy08, smoothMPC_dvcc08);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd08, smoothMPC_dvcc08, smoothMPC_yy07, smoothMPC_bmy07);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld07, smoothMPC_bmy07, smoothMPC_dvcc07);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd07, smoothMPC_dvcc07, smoothMPC_yy06, smoothMPC_bmy06);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld06, smoothMPC_bmy06, smoothMPC_dvcc06);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd06, smoothMPC_dvcc06, smoothMPC_yy05, smoothMPC_bmy05);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld05, smoothMPC_bmy05, smoothMPC_dvcc05);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd05, smoothMPC_dvcc05, smoothMPC_yy04, smoothMPC_bmy04);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld04, smoothMPC_bmy04, smoothMPC_dvcc04);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd04, smoothMPC_dvcc04, smoothMPC_yy03, smoothMPC_bmy03);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld03, smoothMPC_bmy03, smoothMPC_dvcc03);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd03, smoothMPC_dvcc03, smoothMPC_yy02, smoothMPC_bmy02);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld02, smoothMPC_bmy02, smoothMPC_dvcc02);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd02, smoothMPC_dvcc02, smoothMPC_yy01, smoothMPC_bmy01);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld01, smoothMPC_bmy01, smoothMPC_dvcc01);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd01, smoothMPC_dvcc01, smoothMPC_yy00, smoothMPC_bmy00);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld00, smoothMPC_bmy00, smoothMPC_dvcc00);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C1, smoothMPC_dvcc01, smoothMPC_D00, smoothMPC_dvcc00, smoothMPC_grad_eq00);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C2, smoothMPC_dvcc02, smoothMPC_D01, smoothMPC_dvcc01, smoothMPC_grad_eq01);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C3, smoothMPC_dvcc03, smoothMPC_D01, smoothMPC_dvcc02, smoothMPC_grad_eq02);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C4, smoothMPC_dvcc04, smoothMPC_D01, smoothMPC_dvcc03, smoothMPC_grad_eq03);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C5, smoothMPC_dvcc05, smoothMPC_D01, smoothMPC_dvcc04, smoothMPC_grad_eq04);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C6, smoothMPC_dvcc06, smoothMPC_D01, smoothMPC_dvcc05, smoothMPC_grad_eq05);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C7, smoothMPC_dvcc07, smoothMPC_D01, smoothMPC_dvcc06, smoothMPC_grad_eq06);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C8, smoothMPC_dvcc08, smoothMPC_D01, smoothMPC_dvcc07, smoothMPC_grad_eq07);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C9, smoothMPC_dvcc09, smoothMPC_D01, smoothMPC_dvcc08, smoothMPC_grad_eq08);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C10, smoothMPC_dvcc10, smoothMPC_D01, smoothMPC_dvcc09, smoothMPC_grad_eq09);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C11, smoothMPC_dvcc11, smoothMPC_D01, smoothMPC_dvcc10, smoothMPC_grad_eq10);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C12, smoothMPC_dvcc12, smoothMPC_D01, smoothMPC_dvcc11, smoothMPC_grad_eq11);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C13, smoothMPC_dvcc13, smoothMPC_D01, smoothMPC_dvcc12, smoothMPC_grad_eq12);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C14, smoothMPC_dvcc14, smoothMPC_D01, smoothMPC_dvcc13, smoothMPC_grad_eq13);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C15, smoothMPC_dvcc15, smoothMPC_D01, smoothMPC_dvcc14, smoothMPC_grad_eq14);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C16, smoothMPC_dvcc16, smoothMPC_D01, smoothMPC_dvcc15, smoothMPC_grad_eq15);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C17, smoothMPC_dvcc17, smoothMPC_D01, smoothMPC_dvcc16, smoothMPC_grad_eq16);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C18, smoothMPC_dvcc18, smoothMPC_D01, smoothMPC_dvcc17, smoothMPC_grad_eq17);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C19, smoothMPC_dvcc19, smoothMPC_D01, smoothMPC_dvcc18, smoothMPC_grad_eq18);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C20, smoothMPC_dvcc20, smoothMPC_D01, smoothMPC_dvcc19, smoothMPC_grad_eq19);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C21, smoothMPC_dvcc21, smoothMPC_D01, smoothMPC_dvcc20, smoothMPC_grad_eq20);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C22, smoothMPC_dvcc22, smoothMPC_D01, smoothMPC_dvcc21, smoothMPC_grad_eq21);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C23, smoothMPC_dvcc23, smoothMPC_D01, smoothMPC_dvcc22, smoothMPC_grad_eq22);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C24, smoothMPC_dvcc24, smoothMPC_D01, smoothMPC_dvcc23, smoothMPC_grad_eq23);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C25, smoothMPC_dvcc25, smoothMPC_D01, smoothMPC_dvcc24, smoothMPC_grad_eq24);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C26, smoothMPC_dvcc26, smoothMPC_D01, smoothMPC_dvcc25, smoothMPC_grad_eq25);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C27, smoothMPC_dvcc27, smoothMPC_D01, smoothMPC_dvcc26, smoothMPC_grad_eq26);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C28, smoothMPC_dvcc28, smoothMPC_D01, smoothMPC_dvcc27, smoothMPC_grad_eq27);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C29, smoothMPC_dvcc29, smoothMPC_D01, smoothMPC_dvcc28, smoothMPC_grad_eq28);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C30, smoothMPC_dvcc30, smoothMPC_D01, smoothMPC_dvcc29, smoothMPC_grad_eq29);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C31, smoothMPC_dvcc31, smoothMPC_D01, smoothMPC_dvcc30, smoothMPC_grad_eq30);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C32, smoothMPC_dvcc32, smoothMPC_D01, smoothMPC_dvcc31, smoothMPC_grad_eq31);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C33, smoothMPC_dvcc33, smoothMPC_D01, smoothMPC_dvcc32, smoothMPC_grad_eq32);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C34, smoothMPC_dvcc34, smoothMPC_D01, smoothMPC_dvcc33, smoothMPC_grad_eq33);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C35, smoothMPC_dvcc35, smoothMPC_D01, smoothMPC_dvcc34, smoothMPC_grad_eq34);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C36, smoothMPC_dvcc36, smoothMPC_D01, smoothMPC_dvcc35, smoothMPC_grad_eq35);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C37, smoothMPC_dvcc37, smoothMPC_D01, smoothMPC_dvcc36, smoothMPC_grad_eq36);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C38, smoothMPC_dvcc38, smoothMPC_D01, smoothMPC_dvcc37, smoothMPC_grad_eq37);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C39, smoothMPC_dvcc39, smoothMPC_D01, smoothMPC_dvcc38, smoothMPC_grad_eq38);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C40, smoothMPC_dvcc40, smoothMPC_D01, smoothMPC_dvcc39, smoothMPC_grad_eq39);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C41, smoothMPC_dvcc41, smoothMPC_D01, smoothMPC_dvcc40, smoothMPC_grad_eq40);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C42, smoothMPC_dvcc42, smoothMPC_D01, smoothMPC_dvcc41, smoothMPC_grad_eq41);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C43, smoothMPC_dvcc43, smoothMPC_D01, smoothMPC_dvcc42, smoothMPC_grad_eq42);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C44, smoothMPC_dvcc44, smoothMPC_D01, smoothMPC_dvcc43, smoothMPC_grad_eq43);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C45, smoothMPC_dvcc45, smoothMPC_D01, smoothMPC_dvcc44, smoothMPC_grad_eq44);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C46, smoothMPC_dvcc46, smoothMPC_D01, smoothMPC_dvcc45, smoothMPC_grad_eq45);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C47, smoothMPC_dvcc47, smoothMPC_D01, smoothMPC_dvcc46, smoothMPC_grad_eq46);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C48, smoothMPC_dvcc48, smoothMPC_D01, smoothMPC_dvcc47, smoothMPC_grad_eq47);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C49, smoothMPC_dvcc49, smoothMPC_D01, smoothMPC_dvcc48, smoothMPC_grad_eq48);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C50, smoothMPC_dvcc50, smoothMPC_D01, smoothMPC_dvcc49, smoothMPC_grad_eq49);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C51, smoothMPC_dvcc51, smoothMPC_D01, smoothMPC_dvcc50, smoothMPC_grad_eq50);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C52, smoothMPC_dvcc52, smoothMPC_D01, smoothMPC_dvcc51, smoothMPC_grad_eq51);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C53, smoothMPC_dvcc53, smoothMPC_D01, smoothMPC_dvcc52, smoothMPC_grad_eq52);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C54, smoothMPC_dvcc54, smoothMPC_D01, smoothMPC_dvcc53, smoothMPC_grad_eq53);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C55, smoothMPC_dvcc55, smoothMPC_D01, smoothMPC_dvcc54, smoothMPC_grad_eq54);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C56, smoothMPC_dvcc56, smoothMPC_D01, smoothMPC_dvcc55, smoothMPC_grad_eq55);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C57, smoothMPC_dvcc57, smoothMPC_D01, smoothMPC_dvcc56, smoothMPC_grad_eq56);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C58, smoothMPC_dvcc58, smoothMPC_D01, smoothMPC_dvcc57, smoothMPC_grad_eq57);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C59, smoothMPC_dvcc59, smoothMPC_D01, smoothMPC_dvcc58, smoothMPC_grad_eq58);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C60, smoothMPC_dvcc60, smoothMPC_D01, smoothMPC_dvcc59, smoothMPC_grad_eq59);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C61, smoothMPC_dvcc61, smoothMPC_D01, smoothMPC_dvcc60, smoothMPC_grad_eq60);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C62, smoothMPC_dvcc62, smoothMPC_D01, smoothMPC_dvcc61, smoothMPC_grad_eq61);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C63, smoothMPC_dvcc63, smoothMPC_D01, smoothMPC_dvcc62, smoothMPC_grad_eq62);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C64, smoothMPC_dvcc64, smoothMPC_D01, smoothMPC_dvcc63, smoothMPC_grad_eq63);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C65, smoothMPC_dvcc65, smoothMPC_D01, smoothMPC_dvcc64, smoothMPC_grad_eq64);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C66, smoothMPC_dvcc66, smoothMPC_D01, smoothMPC_dvcc65, smoothMPC_grad_eq65);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C67, smoothMPC_dvcc67, smoothMPC_D01, smoothMPC_dvcc66, smoothMPC_grad_eq66);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C68, smoothMPC_dvcc68, smoothMPC_D01, smoothMPC_dvcc67, smoothMPC_grad_eq67);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C69, smoothMPC_dvcc69, smoothMPC_D01, smoothMPC_dvcc68, smoothMPC_grad_eq68);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C70, smoothMPC_dvcc70, smoothMPC_D01, smoothMPC_dvcc69, smoothMPC_grad_eq69);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C71, smoothMPC_dvcc71, smoothMPC_D01, smoothMPC_dvcc70, smoothMPC_grad_eq70);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C72, smoothMPC_dvcc72, smoothMPC_D01, smoothMPC_dvcc71, smoothMPC_grad_eq71);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C73, smoothMPC_dvcc73, smoothMPC_D01, smoothMPC_dvcc72, smoothMPC_grad_eq72);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C74, smoothMPC_dvcc74, smoothMPC_D01, smoothMPC_dvcc73, smoothMPC_grad_eq73);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C75, smoothMPC_dvcc75, smoothMPC_D01, smoothMPC_dvcc74, smoothMPC_grad_eq74);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C76, smoothMPC_dvcc76, smoothMPC_D01, smoothMPC_dvcc75, smoothMPC_grad_eq75);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C77, smoothMPC_dvcc77, smoothMPC_D01, smoothMPC_dvcc76, smoothMPC_grad_eq76);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C78, smoothMPC_dvcc78, smoothMPC_D01, smoothMPC_dvcc77, smoothMPC_grad_eq77);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C79, smoothMPC_dvcc79, smoothMPC_D01, smoothMPC_dvcc78, smoothMPC_grad_eq78);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C80, smoothMPC_dvcc80, smoothMPC_D01, smoothMPC_dvcc79, smoothMPC_grad_eq79);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C81, smoothMPC_dvcc81, smoothMPC_D01, smoothMPC_dvcc80, smoothMPC_grad_eq80);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C82, smoothMPC_dvcc82, smoothMPC_D01, smoothMPC_dvcc81, smoothMPC_grad_eq81);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C83, smoothMPC_dvcc83, smoothMPC_D01, smoothMPC_dvcc82, smoothMPC_grad_eq82);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C84, smoothMPC_dvcc84, smoothMPC_D01, smoothMPC_dvcc83, smoothMPC_grad_eq83);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C85, smoothMPC_dvcc85, smoothMPC_D01, smoothMPC_dvcc84, smoothMPC_grad_eq84);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C86, smoothMPC_dvcc86, smoothMPC_D01, smoothMPC_dvcc85, smoothMPC_grad_eq85);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C87, smoothMPC_dvcc87, smoothMPC_D01, smoothMPC_dvcc86, smoothMPC_grad_eq86);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C88, smoothMPC_dvcc88, smoothMPC_D01, smoothMPC_dvcc87, smoothMPC_grad_eq87);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C89, smoothMPC_dvcc89, smoothMPC_D01, smoothMPC_dvcc88, smoothMPC_grad_eq88);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C90, smoothMPC_dvcc90, smoothMPC_D01, smoothMPC_dvcc89, smoothMPC_grad_eq89);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C91, smoothMPC_dvcc91, smoothMPC_D01, smoothMPC_dvcc90, smoothMPC_grad_eq90);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C92, smoothMPC_dvcc92, smoothMPC_D01, smoothMPC_dvcc91, smoothMPC_grad_eq91);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C93, smoothMPC_dvcc93, smoothMPC_D01, smoothMPC_dvcc92, smoothMPC_grad_eq92);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C94, smoothMPC_dvcc94, smoothMPC_D01, smoothMPC_dvcc93, smoothMPC_grad_eq93);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C95, smoothMPC_dvcc95, smoothMPC_D01, smoothMPC_dvcc94, smoothMPC_grad_eq94);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C96, smoothMPC_dvcc96, smoothMPC_D01, smoothMPC_dvcc95, smoothMPC_grad_eq95);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C97, smoothMPC_dvcc97, smoothMPC_D01, smoothMPC_dvcc96, smoothMPC_grad_eq96);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C98, smoothMPC_dvcc98, smoothMPC_D01, smoothMPC_dvcc97, smoothMPC_grad_eq97);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C99, smoothMPC_dvcc99, smoothMPC_D01, smoothMPC_dvcc98, smoothMPC_grad_eq98);
smoothMPC_LA_DIAGZERO_MTVM_3_3(smoothMPC_D99, smoothMPC_dvcc99, smoothMPC_grad_eq99);
smoothMPC_LA_VSUB_1092(smoothMPC_rd, smoothMPC_grad_eq, smoothMPC_rd);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi00, smoothMPC_rd00, smoothMPC_dzcc00);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi01, smoothMPC_rd01, smoothMPC_dzcc01);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi02, smoothMPC_rd02, smoothMPC_dzcc02);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi03, smoothMPC_rd03, smoothMPC_dzcc03);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi04, smoothMPC_rd04, smoothMPC_dzcc04);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi05, smoothMPC_rd05, smoothMPC_dzcc05);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi06, smoothMPC_rd06, smoothMPC_dzcc06);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi07, smoothMPC_rd07, smoothMPC_dzcc07);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi08, smoothMPC_rd08, smoothMPC_dzcc08);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi09, smoothMPC_rd09, smoothMPC_dzcc09);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi10, smoothMPC_rd10, smoothMPC_dzcc10);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi11, smoothMPC_rd11, smoothMPC_dzcc11);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi12, smoothMPC_rd12, smoothMPC_dzcc12);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi13, smoothMPC_rd13, smoothMPC_dzcc13);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi14, smoothMPC_rd14, smoothMPC_dzcc14);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi15, smoothMPC_rd15, smoothMPC_dzcc15);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi16, smoothMPC_rd16, smoothMPC_dzcc16);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi17, smoothMPC_rd17, smoothMPC_dzcc17);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi18, smoothMPC_rd18, smoothMPC_dzcc18);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi19, smoothMPC_rd19, smoothMPC_dzcc19);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi20, smoothMPC_rd20, smoothMPC_dzcc20);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi21, smoothMPC_rd21, smoothMPC_dzcc21);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi22, smoothMPC_rd22, smoothMPC_dzcc22);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi23, smoothMPC_rd23, smoothMPC_dzcc23);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi24, smoothMPC_rd24, smoothMPC_dzcc24);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi25, smoothMPC_rd25, smoothMPC_dzcc25);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi26, smoothMPC_rd26, smoothMPC_dzcc26);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi27, smoothMPC_rd27, smoothMPC_dzcc27);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi28, smoothMPC_rd28, smoothMPC_dzcc28);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi29, smoothMPC_rd29, smoothMPC_dzcc29);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi30, smoothMPC_rd30, smoothMPC_dzcc30);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi31, smoothMPC_rd31, smoothMPC_dzcc31);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi32, smoothMPC_rd32, smoothMPC_dzcc32);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi33, smoothMPC_rd33, smoothMPC_dzcc33);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi34, smoothMPC_rd34, smoothMPC_dzcc34);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi35, smoothMPC_rd35, smoothMPC_dzcc35);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi36, smoothMPC_rd36, smoothMPC_dzcc36);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi37, smoothMPC_rd37, smoothMPC_dzcc37);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi38, smoothMPC_rd38, smoothMPC_dzcc38);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi39, smoothMPC_rd39, smoothMPC_dzcc39);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi40, smoothMPC_rd40, smoothMPC_dzcc40);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi41, smoothMPC_rd41, smoothMPC_dzcc41);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi42, smoothMPC_rd42, smoothMPC_dzcc42);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi43, smoothMPC_rd43, smoothMPC_dzcc43);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi44, smoothMPC_rd44, smoothMPC_dzcc44);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi45, smoothMPC_rd45, smoothMPC_dzcc45);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi46, smoothMPC_rd46, smoothMPC_dzcc46);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi47, smoothMPC_rd47, smoothMPC_dzcc47);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi48, smoothMPC_rd48, smoothMPC_dzcc48);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi49, smoothMPC_rd49, smoothMPC_dzcc49);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi50, smoothMPC_rd50, smoothMPC_dzcc50);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi51, smoothMPC_rd51, smoothMPC_dzcc51);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi52, smoothMPC_rd52, smoothMPC_dzcc52);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi53, smoothMPC_rd53, smoothMPC_dzcc53);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi54, smoothMPC_rd54, smoothMPC_dzcc54);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi55, smoothMPC_rd55, smoothMPC_dzcc55);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi56, smoothMPC_rd56, smoothMPC_dzcc56);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi57, smoothMPC_rd57, smoothMPC_dzcc57);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi58, smoothMPC_rd58, smoothMPC_dzcc58);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi59, smoothMPC_rd59, smoothMPC_dzcc59);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi60, smoothMPC_rd60, smoothMPC_dzcc60);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi61, smoothMPC_rd61, smoothMPC_dzcc61);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi62, smoothMPC_rd62, smoothMPC_dzcc62);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi63, smoothMPC_rd63, smoothMPC_dzcc63);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi64, smoothMPC_rd64, smoothMPC_dzcc64);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi65, smoothMPC_rd65, smoothMPC_dzcc65);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi66, smoothMPC_rd66, smoothMPC_dzcc66);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi67, smoothMPC_rd67, smoothMPC_dzcc67);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi68, smoothMPC_rd68, smoothMPC_dzcc68);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi69, smoothMPC_rd69, smoothMPC_dzcc69);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi70, smoothMPC_rd70, smoothMPC_dzcc70);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi71, smoothMPC_rd71, smoothMPC_dzcc71);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi72, smoothMPC_rd72, smoothMPC_dzcc72);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi73, smoothMPC_rd73, smoothMPC_dzcc73);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi74, smoothMPC_rd74, smoothMPC_dzcc74);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi75, smoothMPC_rd75, smoothMPC_dzcc75);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi76, smoothMPC_rd76, smoothMPC_dzcc76);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi77, smoothMPC_rd77, smoothMPC_dzcc77);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi78, smoothMPC_rd78, smoothMPC_dzcc78);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi79, smoothMPC_rd79, smoothMPC_dzcc79);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi80, smoothMPC_rd80, smoothMPC_dzcc80);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi81, smoothMPC_rd81, smoothMPC_dzcc81);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi82, smoothMPC_rd82, smoothMPC_dzcc82);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi83, smoothMPC_rd83, smoothMPC_dzcc83);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi84, smoothMPC_rd84, smoothMPC_dzcc84);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi85, smoothMPC_rd85, smoothMPC_dzcc85);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi86, smoothMPC_rd86, smoothMPC_dzcc86);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi87, smoothMPC_rd87, smoothMPC_dzcc87);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi88, smoothMPC_rd88, smoothMPC_dzcc88);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi89, smoothMPC_rd89, smoothMPC_dzcc89);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi90, smoothMPC_rd90, smoothMPC_dzcc90);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi91, smoothMPC_rd91, smoothMPC_dzcc91);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi92, smoothMPC_rd92, smoothMPC_dzcc92);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi93, smoothMPC_rd93, smoothMPC_dzcc93);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi94, smoothMPC_rd94, smoothMPC_dzcc94);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi95, smoothMPC_rd95, smoothMPC_dzcc95);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi96, smoothMPC_rd96, smoothMPC_dzcc96);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi97, smoothMPC_rd97, smoothMPC_dzcc97);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi98, smoothMPC_rd98, smoothMPC_dzcc98);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_3(smoothMPC_Phi99, smoothMPC_rd99, smoothMPC_dzcc99);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl00, smoothMPC_slb00, smoothMPC_llbbyslb00, smoothMPC_dzcc00, smoothMPC_lbIdx00, smoothMPC_dllbcc00);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub00, smoothMPC_sub00, smoothMPC_lubbysub00, smoothMPC_dzcc00, smoothMPC_ubIdx00, smoothMPC_dlubcc00);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl01, smoothMPC_slb01, smoothMPC_llbbyslb01, smoothMPC_dzcc01, smoothMPC_lbIdx01, smoothMPC_dllbcc01);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub01, smoothMPC_sub01, smoothMPC_lubbysub01, smoothMPC_dzcc01, smoothMPC_ubIdx01, smoothMPC_dlubcc01);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl02, smoothMPC_slb02, smoothMPC_llbbyslb02, smoothMPC_dzcc02, smoothMPC_lbIdx02, smoothMPC_dllbcc02);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub02, smoothMPC_sub02, smoothMPC_lubbysub02, smoothMPC_dzcc02, smoothMPC_ubIdx02, smoothMPC_dlubcc02);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl03, smoothMPC_slb03, smoothMPC_llbbyslb03, smoothMPC_dzcc03, smoothMPC_lbIdx03, smoothMPC_dllbcc03);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub03, smoothMPC_sub03, smoothMPC_lubbysub03, smoothMPC_dzcc03, smoothMPC_ubIdx03, smoothMPC_dlubcc03);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl04, smoothMPC_slb04, smoothMPC_llbbyslb04, smoothMPC_dzcc04, smoothMPC_lbIdx04, smoothMPC_dllbcc04);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub04, smoothMPC_sub04, smoothMPC_lubbysub04, smoothMPC_dzcc04, smoothMPC_ubIdx04, smoothMPC_dlubcc04);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl05, smoothMPC_slb05, smoothMPC_llbbyslb05, smoothMPC_dzcc05, smoothMPC_lbIdx05, smoothMPC_dllbcc05);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub05, smoothMPC_sub05, smoothMPC_lubbysub05, smoothMPC_dzcc05, smoothMPC_ubIdx05, smoothMPC_dlubcc05);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl06, smoothMPC_slb06, smoothMPC_llbbyslb06, smoothMPC_dzcc06, smoothMPC_lbIdx06, smoothMPC_dllbcc06);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub06, smoothMPC_sub06, smoothMPC_lubbysub06, smoothMPC_dzcc06, smoothMPC_ubIdx06, smoothMPC_dlubcc06);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl07, smoothMPC_slb07, smoothMPC_llbbyslb07, smoothMPC_dzcc07, smoothMPC_lbIdx07, smoothMPC_dllbcc07);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub07, smoothMPC_sub07, smoothMPC_lubbysub07, smoothMPC_dzcc07, smoothMPC_ubIdx07, smoothMPC_dlubcc07);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl08, smoothMPC_slb08, smoothMPC_llbbyslb08, smoothMPC_dzcc08, smoothMPC_lbIdx08, smoothMPC_dllbcc08);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub08, smoothMPC_sub08, smoothMPC_lubbysub08, smoothMPC_dzcc08, smoothMPC_ubIdx08, smoothMPC_dlubcc08);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl09, smoothMPC_slb09, smoothMPC_llbbyslb09, smoothMPC_dzcc09, smoothMPC_lbIdx09, smoothMPC_dllbcc09);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub09, smoothMPC_sub09, smoothMPC_lubbysub09, smoothMPC_dzcc09, smoothMPC_ubIdx09, smoothMPC_dlubcc09);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl10, smoothMPC_slb10, smoothMPC_llbbyslb10, smoothMPC_dzcc10, smoothMPC_lbIdx10, smoothMPC_dllbcc10);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub10, smoothMPC_sub10, smoothMPC_lubbysub10, smoothMPC_dzcc10, smoothMPC_ubIdx10, smoothMPC_dlubcc10);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl11, smoothMPC_slb11, smoothMPC_llbbyslb11, smoothMPC_dzcc11, smoothMPC_lbIdx11, smoothMPC_dllbcc11);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub11, smoothMPC_sub11, smoothMPC_lubbysub11, smoothMPC_dzcc11, smoothMPC_ubIdx11, smoothMPC_dlubcc11);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl12, smoothMPC_slb12, smoothMPC_llbbyslb12, smoothMPC_dzcc12, smoothMPC_lbIdx12, smoothMPC_dllbcc12);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub12, smoothMPC_sub12, smoothMPC_lubbysub12, smoothMPC_dzcc12, smoothMPC_ubIdx12, smoothMPC_dlubcc12);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl13, smoothMPC_slb13, smoothMPC_llbbyslb13, smoothMPC_dzcc13, smoothMPC_lbIdx13, smoothMPC_dllbcc13);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub13, smoothMPC_sub13, smoothMPC_lubbysub13, smoothMPC_dzcc13, smoothMPC_ubIdx13, smoothMPC_dlubcc13);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl14, smoothMPC_slb14, smoothMPC_llbbyslb14, smoothMPC_dzcc14, smoothMPC_lbIdx14, smoothMPC_dllbcc14);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub14, smoothMPC_sub14, smoothMPC_lubbysub14, smoothMPC_dzcc14, smoothMPC_ubIdx14, smoothMPC_dlubcc14);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl15, smoothMPC_slb15, smoothMPC_llbbyslb15, smoothMPC_dzcc15, smoothMPC_lbIdx15, smoothMPC_dllbcc15);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub15, smoothMPC_sub15, smoothMPC_lubbysub15, smoothMPC_dzcc15, smoothMPC_ubIdx15, smoothMPC_dlubcc15);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl16, smoothMPC_slb16, smoothMPC_llbbyslb16, smoothMPC_dzcc16, smoothMPC_lbIdx16, smoothMPC_dllbcc16);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub16, smoothMPC_sub16, smoothMPC_lubbysub16, smoothMPC_dzcc16, smoothMPC_ubIdx16, smoothMPC_dlubcc16);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl17, smoothMPC_slb17, smoothMPC_llbbyslb17, smoothMPC_dzcc17, smoothMPC_lbIdx17, smoothMPC_dllbcc17);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub17, smoothMPC_sub17, smoothMPC_lubbysub17, smoothMPC_dzcc17, smoothMPC_ubIdx17, smoothMPC_dlubcc17);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl18, smoothMPC_slb18, smoothMPC_llbbyslb18, smoothMPC_dzcc18, smoothMPC_lbIdx18, smoothMPC_dllbcc18);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub18, smoothMPC_sub18, smoothMPC_lubbysub18, smoothMPC_dzcc18, smoothMPC_ubIdx18, smoothMPC_dlubcc18);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl19, smoothMPC_slb19, smoothMPC_llbbyslb19, smoothMPC_dzcc19, smoothMPC_lbIdx19, smoothMPC_dllbcc19);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub19, smoothMPC_sub19, smoothMPC_lubbysub19, smoothMPC_dzcc19, smoothMPC_ubIdx19, smoothMPC_dlubcc19);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl20, smoothMPC_slb20, smoothMPC_llbbyslb20, smoothMPC_dzcc20, smoothMPC_lbIdx20, smoothMPC_dllbcc20);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub20, smoothMPC_sub20, smoothMPC_lubbysub20, smoothMPC_dzcc20, smoothMPC_ubIdx20, smoothMPC_dlubcc20);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl21, smoothMPC_slb21, smoothMPC_llbbyslb21, smoothMPC_dzcc21, smoothMPC_lbIdx21, smoothMPC_dllbcc21);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub21, smoothMPC_sub21, smoothMPC_lubbysub21, smoothMPC_dzcc21, smoothMPC_ubIdx21, smoothMPC_dlubcc21);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl22, smoothMPC_slb22, smoothMPC_llbbyslb22, smoothMPC_dzcc22, smoothMPC_lbIdx22, smoothMPC_dllbcc22);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub22, smoothMPC_sub22, smoothMPC_lubbysub22, smoothMPC_dzcc22, smoothMPC_ubIdx22, smoothMPC_dlubcc22);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl23, smoothMPC_slb23, smoothMPC_llbbyslb23, smoothMPC_dzcc23, smoothMPC_lbIdx23, smoothMPC_dllbcc23);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub23, smoothMPC_sub23, smoothMPC_lubbysub23, smoothMPC_dzcc23, smoothMPC_ubIdx23, smoothMPC_dlubcc23);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl24, smoothMPC_slb24, smoothMPC_llbbyslb24, smoothMPC_dzcc24, smoothMPC_lbIdx24, smoothMPC_dllbcc24);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub24, smoothMPC_sub24, smoothMPC_lubbysub24, smoothMPC_dzcc24, smoothMPC_ubIdx24, smoothMPC_dlubcc24);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl25, smoothMPC_slb25, smoothMPC_llbbyslb25, smoothMPC_dzcc25, smoothMPC_lbIdx25, smoothMPC_dllbcc25);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub25, smoothMPC_sub25, smoothMPC_lubbysub25, smoothMPC_dzcc25, smoothMPC_ubIdx25, smoothMPC_dlubcc25);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl26, smoothMPC_slb26, smoothMPC_llbbyslb26, smoothMPC_dzcc26, smoothMPC_lbIdx26, smoothMPC_dllbcc26);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub26, smoothMPC_sub26, smoothMPC_lubbysub26, smoothMPC_dzcc26, smoothMPC_ubIdx26, smoothMPC_dlubcc26);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl27, smoothMPC_slb27, smoothMPC_llbbyslb27, smoothMPC_dzcc27, smoothMPC_lbIdx27, smoothMPC_dllbcc27);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub27, smoothMPC_sub27, smoothMPC_lubbysub27, smoothMPC_dzcc27, smoothMPC_ubIdx27, smoothMPC_dlubcc27);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl28, smoothMPC_slb28, smoothMPC_llbbyslb28, smoothMPC_dzcc28, smoothMPC_lbIdx28, smoothMPC_dllbcc28);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub28, smoothMPC_sub28, smoothMPC_lubbysub28, smoothMPC_dzcc28, smoothMPC_ubIdx28, smoothMPC_dlubcc28);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl29, smoothMPC_slb29, smoothMPC_llbbyslb29, smoothMPC_dzcc29, smoothMPC_lbIdx29, smoothMPC_dllbcc29);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub29, smoothMPC_sub29, smoothMPC_lubbysub29, smoothMPC_dzcc29, smoothMPC_ubIdx29, smoothMPC_dlubcc29);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl30, smoothMPC_slb30, smoothMPC_llbbyslb30, smoothMPC_dzcc30, smoothMPC_lbIdx30, smoothMPC_dllbcc30);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub30, smoothMPC_sub30, smoothMPC_lubbysub30, smoothMPC_dzcc30, smoothMPC_ubIdx30, smoothMPC_dlubcc30);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl31, smoothMPC_slb31, smoothMPC_llbbyslb31, smoothMPC_dzcc31, smoothMPC_lbIdx31, smoothMPC_dllbcc31);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub31, smoothMPC_sub31, smoothMPC_lubbysub31, smoothMPC_dzcc31, smoothMPC_ubIdx31, smoothMPC_dlubcc31);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl32, smoothMPC_slb32, smoothMPC_llbbyslb32, smoothMPC_dzcc32, smoothMPC_lbIdx32, smoothMPC_dllbcc32);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub32, smoothMPC_sub32, smoothMPC_lubbysub32, smoothMPC_dzcc32, smoothMPC_ubIdx32, smoothMPC_dlubcc32);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl33, smoothMPC_slb33, smoothMPC_llbbyslb33, smoothMPC_dzcc33, smoothMPC_lbIdx33, smoothMPC_dllbcc33);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub33, smoothMPC_sub33, smoothMPC_lubbysub33, smoothMPC_dzcc33, smoothMPC_ubIdx33, smoothMPC_dlubcc33);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl34, smoothMPC_slb34, smoothMPC_llbbyslb34, smoothMPC_dzcc34, smoothMPC_lbIdx34, smoothMPC_dllbcc34);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub34, smoothMPC_sub34, smoothMPC_lubbysub34, smoothMPC_dzcc34, smoothMPC_ubIdx34, smoothMPC_dlubcc34);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl35, smoothMPC_slb35, smoothMPC_llbbyslb35, smoothMPC_dzcc35, smoothMPC_lbIdx35, smoothMPC_dllbcc35);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub35, smoothMPC_sub35, smoothMPC_lubbysub35, smoothMPC_dzcc35, smoothMPC_ubIdx35, smoothMPC_dlubcc35);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl36, smoothMPC_slb36, smoothMPC_llbbyslb36, smoothMPC_dzcc36, smoothMPC_lbIdx36, smoothMPC_dllbcc36);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub36, smoothMPC_sub36, smoothMPC_lubbysub36, smoothMPC_dzcc36, smoothMPC_ubIdx36, smoothMPC_dlubcc36);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl37, smoothMPC_slb37, smoothMPC_llbbyslb37, smoothMPC_dzcc37, smoothMPC_lbIdx37, smoothMPC_dllbcc37);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub37, smoothMPC_sub37, smoothMPC_lubbysub37, smoothMPC_dzcc37, smoothMPC_ubIdx37, smoothMPC_dlubcc37);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl38, smoothMPC_slb38, smoothMPC_llbbyslb38, smoothMPC_dzcc38, smoothMPC_lbIdx38, smoothMPC_dllbcc38);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub38, smoothMPC_sub38, smoothMPC_lubbysub38, smoothMPC_dzcc38, smoothMPC_ubIdx38, smoothMPC_dlubcc38);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl39, smoothMPC_slb39, smoothMPC_llbbyslb39, smoothMPC_dzcc39, smoothMPC_lbIdx39, smoothMPC_dllbcc39);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub39, smoothMPC_sub39, smoothMPC_lubbysub39, smoothMPC_dzcc39, smoothMPC_ubIdx39, smoothMPC_dlubcc39);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl40, smoothMPC_slb40, smoothMPC_llbbyslb40, smoothMPC_dzcc40, smoothMPC_lbIdx40, smoothMPC_dllbcc40);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub40, smoothMPC_sub40, smoothMPC_lubbysub40, smoothMPC_dzcc40, smoothMPC_ubIdx40, smoothMPC_dlubcc40);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl41, smoothMPC_slb41, smoothMPC_llbbyslb41, smoothMPC_dzcc41, smoothMPC_lbIdx41, smoothMPC_dllbcc41);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub41, smoothMPC_sub41, smoothMPC_lubbysub41, smoothMPC_dzcc41, smoothMPC_ubIdx41, smoothMPC_dlubcc41);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl42, smoothMPC_slb42, smoothMPC_llbbyslb42, smoothMPC_dzcc42, smoothMPC_lbIdx42, smoothMPC_dllbcc42);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub42, smoothMPC_sub42, smoothMPC_lubbysub42, smoothMPC_dzcc42, smoothMPC_ubIdx42, smoothMPC_dlubcc42);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl43, smoothMPC_slb43, smoothMPC_llbbyslb43, smoothMPC_dzcc43, smoothMPC_lbIdx43, smoothMPC_dllbcc43);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub43, smoothMPC_sub43, smoothMPC_lubbysub43, smoothMPC_dzcc43, smoothMPC_ubIdx43, smoothMPC_dlubcc43);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl44, smoothMPC_slb44, smoothMPC_llbbyslb44, smoothMPC_dzcc44, smoothMPC_lbIdx44, smoothMPC_dllbcc44);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub44, smoothMPC_sub44, smoothMPC_lubbysub44, smoothMPC_dzcc44, smoothMPC_ubIdx44, smoothMPC_dlubcc44);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl45, smoothMPC_slb45, smoothMPC_llbbyslb45, smoothMPC_dzcc45, smoothMPC_lbIdx45, smoothMPC_dllbcc45);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub45, smoothMPC_sub45, smoothMPC_lubbysub45, smoothMPC_dzcc45, smoothMPC_ubIdx45, smoothMPC_dlubcc45);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl46, smoothMPC_slb46, smoothMPC_llbbyslb46, smoothMPC_dzcc46, smoothMPC_lbIdx46, smoothMPC_dllbcc46);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub46, smoothMPC_sub46, smoothMPC_lubbysub46, smoothMPC_dzcc46, smoothMPC_ubIdx46, smoothMPC_dlubcc46);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl47, smoothMPC_slb47, smoothMPC_llbbyslb47, smoothMPC_dzcc47, smoothMPC_lbIdx47, smoothMPC_dllbcc47);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub47, smoothMPC_sub47, smoothMPC_lubbysub47, smoothMPC_dzcc47, smoothMPC_ubIdx47, smoothMPC_dlubcc47);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl48, smoothMPC_slb48, smoothMPC_llbbyslb48, smoothMPC_dzcc48, smoothMPC_lbIdx48, smoothMPC_dllbcc48);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub48, smoothMPC_sub48, smoothMPC_lubbysub48, smoothMPC_dzcc48, smoothMPC_ubIdx48, smoothMPC_dlubcc48);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl49, smoothMPC_slb49, smoothMPC_llbbyslb49, smoothMPC_dzcc49, smoothMPC_lbIdx49, smoothMPC_dllbcc49);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub49, smoothMPC_sub49, smoothMPC_lubbysub49, smoothMPC_dzcc49, smoothMPC_ubIdx49, smoothMPC_dlubcc49);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl50, smoothMPC_slb50, smoothMPC_llbbyslb50, smoothMPC_dzcc50, smoothMPC_lbIdx50, smoothMPC_dllbcc50);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub50, smoothMPC_sub50, smoothMPC_lubbysub50, smoothMPC_dzcc50, smoothMPC_ubIdx50, smoothMPC_dlubcc50);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl51, smoothMPC_slb51, smoothMPC_llbbyslb51, smoothMPC_dzcc51, smoothMPC_lbIdx51, smoothMPC_dllbcc51);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub51, smoothMPC_sub51, smoothMPC_lubbysub51, smoothMPC_dzcc51, smoothMPC_ubIdx51, smoothMPC_dlubcc51);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl52, smoothMPC_slb52, smoothMPC_llbbyslb52, smoothMPC_dzcc52, smoothMPC_lbIdx52, smoothMPC_dllbcc52);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub52, smoothMPC_sub52, smoothMPC_lubbysub52, smoothMPC_dzcc52, smoothMPC_ubIdx52, smoothMPC_dlubcc52);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl53, smoothMPC_slb53, smoothMPC_llbbyslb53, smoothMPC_dzcc53, smoothMPC_lbIdx53, smoothMPC_dllbcc53);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub53, smoothMPC_sub53, smoothMPC_lubbysub53, smoothMPC_dzcc53, smoothMPC_ubIdx53, smoothMPC_dlubcc53);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl54, smoothMPC_slb54, smoothMPC_llbbyslb54, smoothMPC_dzcc54, smoothMPC_lbIdx54, smoothMPC_dllbcc54);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub54, smoothMPC_sub54, smoothMPC_lubbysub54, smoothMPC_dzcc54, smoothMPC_ubIdx54, smoothMPC_dlubcc54);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl55, smoothMPC_slb55, smoothMPC_llbbyslb55, smoothMPC_dzcc55, smoothMPC_lbIdx55, smoothMPC_dllbcc55);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub55, smoothMPC_sub55, smoothMPC_lubbysub55, smoothMPC_dzcc55, smoothMPC_ubIdx55, smoothMPC_dlubcc55);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl56, smoothMPC_slb56, smoothMPC_llbbyslb56, smoothMPC_dzcc56, smoothMPC_lbIdx56, smoothMPC_dllbcc56);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub56, smoothMPC_sub56, smoothMPC_lubbysub56, smoothMPC_dzcc56, smoothMPC_ubIdx56, smoothMPC_dlubcc56);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl57, smoothMPC_slb57, smoothMPC_llbbyslb57, smoothMPC_dzcc57, smoothMPC_lbIdx57, smoothMPC_dllbcc57);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub57, smoothMPC_sub57, smoothMPC_lubbysub57, smoothMPC_dzcc57, smoothMPC_ubIdx57, smoothMPC_dlubcc57);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl58, smoothMPC_slb58, smoothMPC_llbbyslb58, smoothMPC_dzcc58, smoothMPC_lbIdx58, smoothMPC_dllbcc58);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub58, smoothMPC_sub58, smoothMPC_lubbysub58, smoothMPC_dzcc58, smoothMPC_ubIdx58, smoothMPC_dlubcc58);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl59, smoothMPC_slb59, smoothMPC_llbbyslb59, smoothMPC_dzcc59, smoothMPC_lbIdx59, smoothMPC_dllbcc59);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub59, smoothMPC_sub59, smoothMPC_lubbysub59, smoothMPC_dzcc59, smoothMPC_ubIdx59, smoothMPC_dlubcc59);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl60, smoothMPC_slb60, smoothMPC_llbbyslb60, smoothMPC_dzcc60, smoothMPC_lbIdx60, smoothMPC_dllbcc60);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub60, smoothMPC_sub60, smoothMPC_lubbysub60, smoothMPC_dzcc60, smoothMPC_ubIdx60, smoothMPC_dlubcc60);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl61, smoothMPC_slb61, smoothMPC_llbbyslb61, smoothMPC_dzcc61, smoothMPC_lbIdx61, smoothMPC_dllbcc61);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub61, smoothMPC_sub61, smoothMPC_lubbysub61, smoothMPC_dzcc61, smoothMPC_ubIdx61, smoothMPC_dlubcc61);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl62, smoothMPC_slb62, smoothMPC_llbbyslb62, smoothMPC_dzcc62, smoothMPC_lbIdx62, smoothMPC_dllbcc62);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub62, smoothMPC_sub62, smoothMPC_lubbysub62, smoothMPC_dzcc62, smoothMPC_ubIdx62, smoothMPC_dlubcc62);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl63, smoothMPC_slb63, smoothMPC_llbbyslb63, smoothMPC_dzcc63, smoothMPC_lbIdx63, smoothMPC_dllbcc63);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub63, smoothMPC_sub63, smoothMPC_lubbysub63, smoothMPC_dzcc63, smoothMPC_ubIdx63, smoothMPC_dlubcc63);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl64, smoothMPC_slb64, smoothMPC_llbbyslb64, smoothMPC_dzcc64, smoothMPC_lbIdx64, smoothMPC_dllbcc64);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub64, smoothMPC_sub64, smoothMPC_lubbysub64, smoothMPC_dzcc64, smoothMPC_ubIdx64, smoothMPC_dlubcc64);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl65, smoothMPC_slb65, smoothMPC_llbbyslb65, smoothMPC_dzcc65, smoothMPC_lbIdx65, smoothMPC_dllbcc65);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub65, smoothMPC_sub65, smoothMPC_lubbysub65, smoothMPC_dzcc65, smoothMPC_ubIdx65, smoothMPC_dlubcc65);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl66, smoothMPC_slb66, smoothMPC_llbbyslb66, smoothMPC_dzcc66, smoothMPC_lbIdx66, smoothMPC_dllbcc66);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub66, smoothMPC_sub66, smoothMPC_lubbysub66, smoothMPC_dzcc66, smoothMPC_ubIdx66, smoothMPC_dlubcc66);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl67, smoothMPC_slb67, smoothMPC_llbbyslb67, smoothMPC_dzcc67, smoothMPC_lbIdx67, smoothMPC_dllbcc67);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub67, smoothMPC_sub67, smoothMPC_lubbysub67, smoothMPC_dzcc67, smoothMPC_ubIdx67, smoothMPC_dlubcc67);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl68, smoothMPC_slb68, smoothMPC_llbbyslb68, smoothMPC_dzcc68, smoothMPC_lbIdx68, smoothMPC_dllbcc68);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub68, smoothMPC_sub68, smoothMPC_lubbysub68, smoothMPC_dzcc68, smoothMPC_ubIdx68, smoothMPC_dlubcc68);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl69, smoothMPC_slb69, smoothMPC_llbbyslb69, smoothMPC_dzcc69, smoothMPC_lbIdx69, smoothMPC_dllbcc69);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub69, smoothMPC_sub69, smoothMPC_lubbysub69, smoothMPC_dzcc69, smoothMPC_ubIdx69, smoothMPC_dlubcc69);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl70, smoothMPC_slb70, smoothMPC_llbbyslb70, smoothMPC_dzcc70, smoothMPC_lbIdx70, smoothMPC_dllbcc70);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub70, smoothMPC_sub70, smoothMPC_lubbysub70, smoothMPC_dzcc70, smoothMPC_ubIdx70, smoothMPC_dlubcc70);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl71, smoothMPC_slb71, smoothMPC_llbbyslb71, smoothMPC_dzcc71, smoothMPC_lbIdx71, smoothMPC_dllbcc71);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub71, smoothMPC_sub71, smoothMPC_lubbysub71, smoothMPC_dzcc71, smoothMPC_ubIdx71, smoothMPC_dlubcc71);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl72, smoothMPC_slb72, smoothMPC_llbbyslb72, smoothMPC_dzcc72, smoothMPC_lbIdx72, smoothMPC_dllbcc72);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub72, smoothMPC_sub72, smoothMPC_lubbysub72, smoothMPC_dzcc72, smoothMPC_ubIdx72, smoothMPC_dlubcc72);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl73, smoothMPC_slb73, smoothMPC_llbbyslb73, smoothMPC_dzcc73, smoothMPC_lbIdx73, smoothMPC_dllbcc73);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub73, smoothMPC_sub73, smoothMPC_lubbysub73, smoothMPC_dzcc73, smoothMPC_ubIdx73, smoothMPC_dlubcc73);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl74, smoothMPC_slb74, smoothMPC_llbbyslb74, smoothMPC_dzcc74, smoothMPC_lbIdx74, smoothMPC_dllbcc74);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub74, smoothMPC_sub74, smoothMPC_lubbysub74, smoothMPC_dzcc74, smoothMPC_ubIdx74, smoothMPC_dlubcc74);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl75, smoothMPC_slb75, smoothMPC_llbbyslb75, smoothMPC_dzcc75, smoothMPC_lbIdx75, smoothMPC_dllbcc75);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub75, smoothMPC_sub75, smoothMPC_lubbysub75, smoothMPC_dzcc75, smoothMPC_ubIdx75, smoothMPC_dlubcc75);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl76, smoothMPC_slb76, smoothMPC_llbbyslb76, smoothMPC_dzcc76, smoothMPC_lbIdx76, smoothMPC_dllbcc76);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub76, smoothMPC_sub76, smoothMPC_lubbysub76, smoothMPC_dzcc76, smoothMPC_ubIdx76, smoothMPC_dlubcc76);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl77, smoothMPC_slb77, smoothMPC_llbbyslb77, smoothMPC_dzcc77, smoothMPC_lbIdx77, smoothMPC_dllbcc77);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub77, smoothMPC_sub77, smoothMPC_lubbysub77, smoothMPC_dzcc77, smoothMPC_ubIdx77, smoothMPC_dlubcc77);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl78, smoothMPC_slb78, smoothMPC_llbbyslb78, smoothMPC_dzcc78, smoothMPC_lbIdx78, smoothMPC_dllbcc78);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub78, smoothMPC_sub78, smoothMPC_lubbysub78, smoothMPC_dzcc78, smoothMPC_ubIdx78, smoothMPC_dlubcc78);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl79, smoothMPC_slb79, smoothMPC_llbbyslb79, smoothMPC_dzcc79, smoothMPC_lbIdx79, smoothMPC_dllbcc79);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub79, smoothMPC_sub79, smoothMPC_lubbysub79, smoothMPC_dzcc79, smoothMPC_ubIdx79, smoothMPC_dlubcc79);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl80, smoothMPC_slb80, smoothMPC_llbbyslb80, smoothMPC_dzcc80, smoothMPC_lbIdx80, smoothMPC_dllbcc80);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub80, smoothMPC_sub80, smoothMPC_lubbysub80, smoothMPC_dzcc80, smoothMPC_ubIdx80, smoothMPC_dlubcc80);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl81, smoothMPC_slb81, smoothMPC_llbbyslb81, smoothMPC_dzcc81, smoothMPC_lbIdx81, smoothMPC_dllbcc81);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub81, smoothMPC_sub81, smoothMPC_lubbysub81, smoothMPC_dzcc81, smoothMPC_ubIdx81, smoothMPC_dlubcc81);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl82, smoothMPC_slb82, smoothMPC_llbbyslb82, smoothMPC_dzcc82, smoothMPC_lbIdx82, smoothMPC_dllbcc82);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub82, smoothMPC_sub82, smoothMPC_lubbysub82, smoothMPC_dzcc82, smoothMPC_ubIdx82, smoothMPC_dlubcc82);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl83, smoothMPC_slb83, smoothMPC_llbbyslb83, smoothMPC_dzcc83, smoothMPC_lbIdx83, smoothMPC_dllbcc83);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub83, smoothMPC_sub83, smoothMPC_lubbysub83, smoothMPC_dzcc83, smoothMPC_ubIdx83, smoothMPC_dlubcc83);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl84, smoothMPC_slb84, smoothMPC_llbbyslb84, smoothMPC_dzcc84, smoothMPC_lbIdx84, smoothMPC_dllbcc84);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub84, smoothMPC_sub84, smoothMPC_lubbysub84, smoothMPC_dzcc84, smoothMPC_ubIdx84, smoothMPC_dlubcc84);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl85, smoothMPC_slb85, smoothMPC_llbbyslb85, smoothMPC_dzcc85, smoothMPC_lbIdx85, smoothMPC_dllbcc85);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub85, smoothMPC_sub85, smoothMPC_lubbysub85, smoothMPC_dzcc85, smoothMPC_ubIdx85, smoothMPC_dlubcc85);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl86, smoothMPC_slb86, smoothMPC_llbbyslb86, smoothMPC_dzcc86, smoothMPC_lbIdx86, smoothMPC_dllbcc86);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub86, smoothMPC_sub86, smoothMPC_lubbysub86, smoothMPC_dzcc86, smoothMPC_ubIdx86, smoothMPC_dlubcc86);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl87, smoothMPC_slb87, smoothMPC_llbbyslb87, smoothMPC_dzcc87, smoothMPC_lbIdx87, smoothMPC_dllbcc87);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub87, smoothMPC_sub87, smoothMPC_lubbysub87, smoothMPC_dzcc87, smoothMPC_ubIdx87, smoothMPC_dlubcc87);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl88, smoothMPC_slb88, smoothMPC_llbbyslb88, smoothMPC_dzcc88, smoothMPC_lbIdx88, smoothMPC_dllbcc88);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub88, smoothMPC_sub88, smoothMPC_lubbysub88, smoothMPC_dzcc88, smoothMPC_ubIdx88, smoothMPC_dlubcc88);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl89, smoothMPC_slb89, smoothMPC_llbbyslb89, smoothMPC_dzcc89, smoothMPC_lbIdx89, smoothMPC_dllbcc89);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub89, smoothMPC_sub89, smoothMPC_lubbysub89, smoothMPC_dzcc89, smoothMPC_ubIdx89, smoothMPC_dlubcc89);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl90, smoothMPC_slb90, smoothMPC_llbbyslb90, smoothMPC_dzcc90, smoothMPC_lbIdx90, smoothMPC_dllbcc90);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub90, smoothMPC_sub90, smoothMPC_lubbysub90, smoothMPC_dzcc90, smoothMPC_ubIdx90, smoothMPC_dlubcc90);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl91, smoothMPC_slb91, smoothMPC_llbbyslb91, smoothMPC_dzcc91, smoothMPC_lbIdx91, smoothMPC_dllbcc91);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub91, smoothMPC_sub91, smoothMPC_lubbysub91, smoothMPC_dzcc91, smoothMPC_ubIdx91, smoothMPC_dlubcc91);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl92, smoothMPC_slb92, smoothMPC_llbbyslb92, smoothMPC_dzcc92, smoothMPC_lbIdx92, smoothMPC_dllbcc92);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub92, smoothMPC_sub92, smoothMPC_lubbysub92, smoothMPC_dzcc92, smoothMPC_ubIdx92, smoothMPC_dlubcc92);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl93, smoothMPC_slb93, smoothMPC_llbbyslb93, smoothMPC_dzcc93, smoothMPC_lbIdx93, smoothMPC_dllbcc93);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub93, smoothMPC_sub93, smoothMPC_lubbysub93, smoothMPC_dzcc93, smoothMPC_ubIdx93, smoothMPC_dlubcc93);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl94, smoothMPC_slb94, smoothMPC_llbbyslb94, smoothMPC_dzcc94, smoothMPC_lbIdx94, smoothMPC_dllbcc94);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub94, smoothMPC_sub94, smoothMPC_lubbysub94, smoothMPC_dzcc94, smoothMPC_ubIdx94, smoothMPC_dlubcc94);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl95, smoothMPC_slb95, smoothMPC_llbbyslb95, smoothMPC_dzcc95, smoothMPC_lbIdx95, smoothMPC_dllbcc95);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub95, smoothMPC_sub95, smoothMPC_lubbysub95, smoothMPC_dzcc95, smoothMPC_ubIdx95, smoothMPC_dlubcc95);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl96, smoothMPC_slb96, smoothMPC_llbbyslb96, smoothMPC_dzcc96, smoothMPC_lbIdx96, smoothMPC_dllbcc96);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub96, smoothMPC_sub96, smoothMPC_lubbysub96, smoothMPC_dzcc96, smoothMPC_ubIdx96, smoothMPC_dlubcc96);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl97, smoothMPC_slb97, smoothMPC_llbbyslb97, smoothMPC_dzcc97, smoothMPC_lbIdx97, smoothMPC_dllbcc97);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub97, smoothMPC_sub97, smoothMPC_lubbysub97, smoothMPC_dzcc97, smoothMPC_ubIdx97, smoothMPC_dlubcc97);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl98, smoothMPC_slb98, smoothMPC_llbbyslb98, smoothMPC_dzcc98, smoothMPC_lbIdx98, smoothMPC_dllbcc98);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub98, smoothMPC_sub98, smoothMPC_lubbysub98, smoothMPC_dzcc98, smoothMPC_ubIdx98, smoothMPC_dlubcc98);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_3(smoothMPC_ccrhsl99, smoothMPC_slb99, smoothMPC_llbbyslb99, smoothMPC_dzcc99, smoothMPC_lbIdx99, smoothMPC_dllbcc99);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_3(smoothMPC_ccrhsub99, smoothMPC_sub99, smoothMPC_lubbysub99, smoothMPC_dzcc99, smoothMPC_ubIdx99, smoothMPC_dlubcc99);
smoothMPC_LA_VSUB7_1590(smoothMPC_l, smoothMPC_ccrhs, smoothMPC_s, smoothMPC_dl_cc, smoothMPC_ds_cc);
smoothMPC_LA_VADD_1092(smoothMPC_dz_cc, smoothMPC_dz_aff);
smoothMPC_LA_VADD_300(smoothMPC_dv_cc, smoothMPC_dv_aff);
smoothMPC_LA_VADD_1590(smoothMPC_dl_cc, smoothMPC_dl_aff);
smoothMPC_LA_VADD_1590(smoothMPC_ds_cc, smoothMPC_ds_aff);
info->lsit_cc = smoothMPC_LINESEARCH_BACKTRACKING_COMBINED(smoothMPC_z, smoothMPC_v, smoothMPC_l, smoothMPC_s, smoothMPC_dz_cc, smoothMPC_dv_cc, smoothMPC_dl_cc, smoothMPC_ds_cc, &info->step_cc, &info->mu);
if( info->lsit_cc == smoothMPC_NOPROGRESS ){
PRINTTEXT("Line search could not proceed at iteration %d, exiting.\n",info->it+1);
exitcode = smoothMPC_NOPROGRESS; break;
}
info->it++;
}
output->z1[0] = smoothMPC_z00[0];
output->z1[1] = smoothMPC_z00[1];
output->z1[2] = smoothMPC_z00[2];
output->z1[3] = smoothMPC_z00[3];
output->z1[4] = smoothMPC_z00[4];
output->z2[0] = smoothMPC_z01[0];
output->z2[1] = smoothMPC_z01[1];
output->z2[2] = smoothMPC_z01[2];
output->z2[3] = smoothMPC_z01[3];
output->z2[4] = smoothMPC_z01[4];
output->z3[0] = smoothMPC_z02[0];
output->z3[1] = smoothMPC_z02[1];
output->z3[2] = smoothMPC_z02[2];
output->z3[3] = smoothMPC_z02[3];
output->z3[4] = smoothMPC_z02[4];
output->z4[0] = smoothMPC_z03[0];
output->z4[1] = smoothMPC_z03[1];
output->z4[2] = smoothMPC_z03[2];
output->z4[3] = smoothMPC_z03[3];
output->z4[4] = smoothMPC_z03[4];
output->z5[0] = smoothMPC_z04[0];
output->z5[1] = smoothMPC_z04[1];
output->z5[2] = smoothMPC_z04[2];
output->z5[3] = smoothMPC_z04[3];
output->z5[4] = smoothMPC_z04[4];
output->z6[0] = smoothMPC_z05[0];
output->z6[1] = smoothMPC_z05[1];
output->z6[2] = smoothMPC_z05[2];
output->z6[3] = smoothMPC_z05[3];
output->z6[4] = smoothMPC_z05[4];
output->z7[0] = smoothMPC_z06[0];
output->z7[1] = smoothMPC_z06[1];
output->z7[2] = smoothMPC_z06[2];
output->z7[3] = smoothMPC_z06[3];
output->z7[4] = smoothMPC_z06[4];
output->z8[0] = smoothMPC_z07[0];
output->z8[1] = smoothMPC_z07[1];
output->z8[2] = smoothMPC_z07[2];
output->z8[3] = smoothMPC_z07[3];
output->z8[4] = smoothMPC_z07[4];
output->z9[0] = smoothMPC_z08[0];
output->z9[1] = smoothMPC_z08[1];
output->z9[2] = smoothMPC_z08[2];
output->z9[3] = smoothMPC_z08[3];
output->z9[4] = smoothMPC_z08[4];
output->z10[0] = smoothMPC_z09[0];
output->z10[1] = smoothMPC_z09[1];
output->z10[2] = smoothMPC_z09[2];
output->z10[3] = smoothMPC_z09[3];
output->z10[4] = smoothMPC_z09[4];
output->z11[0] = smoothMPC_z10[0];
output->z11[1] = smoothMPC_z10[1];
output->z11[2] = smoothMPC_z10[2];
output->z11[3] = smoothMPC_z10[3];
output->z11[4] = smoothMPC_z10[4];
output->z12[0] = smoothMPC_z11[0];
output->z12[1] = smoothMPC_z11[1];
output->z12[2] = smoothMPC_z11[2];
output->z12[3] = smoothMPC_z11[3];
output->z12[4] = smoothMPC_z11[4];
output->z13[0] = smoothMPC_z12[0];
output->z13[1] = smoothMPC_z12[1];
output->z13[2] = smoothMPC_z12[2];
output->z13[3] = smoothMPC_z12[3];
output->z13[4] = smoothMPC_z12[4];
output->z14[0] = smoothMPC_z13[0];
output->z14[1] = smoothMPC_z13[1];
output->z14[2] = smoothMPC_z13[2];
output->z14[3] = smoothMPC_z13[3];
output->z14[4] = smoothMPC_z13[4];
output->z15[0] = smoothMPC_z14[0];
output->z15[1] = smoothMPC_z14[1];
output->z15[2] = smoothMPC_z14[2];
output->z15[3] = smoothMPC_z14[3];
output->z15[4] = smoothMPC_z14[4];
output->z16[0] = smoothMPC_z15[0];
output->z16[1] = smoothMPC_z15[1];
output->z16[2] = smoothMPC_z15[2];
output->z16[3] = smoothMPC_z15[3];
output->z16[4] = smoothMPC_z15[4];
output->z17[0] = smoothMPC_z16[0];
output->z17[1] = smoothMPC_z16[1];
output->z17[2] = smoothMPC_z16[2];
output->z17[3] = smoothMPC_z16[3];
output->z17[4] = smoothMPC_z16[4];
output->z18[0] = smoothMPC_z17[0];
output->z18[1] = smoothMPC_z17[1];
output->z18[2] = smoothMPC_z17[2];
output->z18[3] = smoothMPC_z17[3];
output->z18[4] = smoothMPC_z17[4];
output->z19[0] = smoothMPC_z18[0];
output->z19[1] = smoothMPC_z18[1];
output->z19[2] = smoothMPC_z18[2];
output->z19[3] = smoothMPC_z18[3];
output->z19[4] = smoothMPC_z18[4];
output->z20[0] = smoothMPC_z19[0];
output->z20[1] = smoothMPC_z19[1];
output->z20[2] = smoothMPC_z19[2];
output->z20[3] = smoothMPC_z19[3];
output->z20[4] = smoothMPC_z19[4];
output->z21[0] = smoothMPC_z20[0];
output->z21[1] = smoothMPC_z20[1];
output->z21[2] = smoothMPC_z20[2];
output->z21[3] = smoothMPC_z20[3];
output->z21[4] = smoothMPC_z20[4];
output->z22[0] = smoothMPC_z21[0];
output->z22[1] = smoothMPC_z21[1];
output->z22[2] = smoothMPC_z21[2];
output->z22[3] = smoothMPC_z21[3];
output->z22[4] = smoothMPC_z21[4];
output->z23[0] = smoothMPC_z22[0];
output->z23[1] = smoothMPC_z22[1];
output->z23[2] = smoothMPC_z22[2];
output->z23[3] = smoothMPC_z22[3];
output->z23[4] = smoothMPC_z22[4];
output->z24[0] = smoothMPC_z23[0];
output->z24[1] = smoothMPC_z23[1];
output->z24[2] = smoothMPC_z23[2];
output->z24[3] = smoothMPC_z23[3];
output->z24[4] = smoothMPC_z23[4];
output->z25[0] = smoothMPC_z24[0];
output->z25[1] = smoothMPC_z24[1];
output->z25[2] = smoothMPC_z24[2];
output->z25[3] = smoothMPC_z24[3];
output->z25[4] = smoothMPC_z24[4];
output->z26[0] = smoothMPC_z25[0];
output->z26[1] = smoothMPC_z25[1];
output->z26[2] = smoothMPC_z25[2];
output->z26[3] = smoothMPC_z25[3];
output->z26[4] = smoothMPC_z25[4];
output->z27[0] = smoothMPC_z26[0];
output->z27[1] = smoothMPC_z26[1];
output->z27[2] = smoothMPC_z26[2];
output->z27[3] = smoothMPC_z26[3];
output->z27[4] = smoothMPC_z26[4];
output->z28[0] = smoothMPC_z27[0];
output->z28[1] = smoothMPC_z27[1];
output->z28[2] = smoothMPC_z27[2];
output->z28[3] = smoothMPC_z27[3];
output->z28[4] = smoothMPC_z27[4];
output->z29[0] = smoothMPC_z28[0];
output->z29[1] = smoothMPC_z28[1];
output->z29[2] = smoothMPC_z28[2];
output->z29[3] = smoothMPC_z28[3];
output->z29[4] = smoothMPC_z28[4];
output->z30[0] = smoothMPC_z29[0];
output->z30[1] = smoothMPC_z29[1];
output->z30[2] = smoothMPC_z29[2];
output->z30[3] = smoothMPC_z29[3];
output->z30[4] = smoothMPC_z29[4];
output->z31[0] = smoothMPC_z30[0];
output->z31[1] = smoothMPC_z30[1];
output->z31[2] = smoothMPC_z30[2];
output->z31[3] = smoothMPC_z30[3];
output->z31[4] = smoothMPC_z30[4];
output->z32[0] = smoothMPC_z31[0];
output->z32[1] = smoothMPC_z31[1];
output->z32[2] = smoothMPC_z31[2];
output->z32[3] = smoothMPC_z31[3];
output->z32[4] = smoothMPC_z31[4];
output->z33[0] = smoothMPC_z32[0];
output->z33[1] = smoothMPC_z32[1];
output->z33[2] = smoothMPC_z32[2];
output->z33[3] = smoothMPC_z32[3];
output->z33[4] = smoothMPC_z32[4];
output->z34[0] = smoothMPC_z33[0];
output->z34[1] = smoothMPC_z33[1];
output->z34[2] = smoothMPC_z33[2];
output->z34[3] = smoothMPC_z33[3];
output->z34[4] = smoothMPC_z33[4];
output->z35[0] = smoothMPC_z34[0];
output->z35[1] = smoothMPC_z34[1];
output->z35[2] = smoothMPC_z34[2];
output->z35[3] = smoothMPC_z34[3];
output->z35[4] = smoothMPC_z34[4];
output->z36[0] = smoothMPC_z35[0];
output->z36[1] = smoothMPC_z35[1];
output->z36[2] = smoothMPC_z35[2];
output->z36[3] = smoothMPC_z35[3];
output->z36[4] = smoothMPC_z35[4];
output->z37[0] = smoothMPC_z36[0];
output->z37[1] = smoothMPC_z36[1];
output->z37[2] = smoothMPC_z36[2];
output->z37[3] = smoothMPC_z36[3];
output->z37[4] = smoothMPC_z36[4];
output->z38[0] = smoothMPC_z37[0];
output->z38[1] = smoothMPC_z37[1];
output->z38[2] = smoothMPC_z37[2];
output->z38[3] = smoothMPC_z37[3];
output->z38[4] = smoothMPC_z37[4];
output->z39[0] = smoothMPC_z38[0];
output->z39[1] = smoothMPC_z38[1];
output->z39[2] = smoothMPC_z38[2];
output->z39[3] = smoothMPC_z38[3];
output->z39[4] = smoothMPC_z38[4];
output->z40[0] = smoothMPC_z39[0];
output->z40[1] = smoothMPC_z39[1];
output->z40[2] = smoothMPC_z39[2];
output->z40[3] = smoothMPC_z39[3];
output->z40[4] = smoothMPC_z39[4];
output->z41[0] = smoothMPC_z40[0];
output->z41[1] = smoothMPC_z40[1];
output->z41[2] = smoothMPC_z40[2];
output->z41[3] = smoothMPC_z40[3];
output->z41[4] = smoothMPC_z40[4];
output->z42[0] = smoothMPC_z41[0];
output->z42[1] = smoothMPC_z41[1];
output->z42[2] = smoothMPC_z41[2];
output->z42[3] = smoothMPC_z41[3];
output->z42[4] = smoothMPC_z41[4];
output->z43[0] = smoothMPC_z42[0];
output->z43[1] = smoothMPC_z42[1];
output->z43[2] = smoothMPC_z42[2];
output->z43[3] = smoothMPC_z42[3];
output->z43[4] = smoothMPC_z42[4];
output->z44[0] = smoothMPC_z43[0];
output->z44[1] = smoothMPC_z43[1];
output->z44[2] = smoothMPC_z43[2];
output->z44[3] = smoothMPC_z43[3];
output->z44[4] = smoothMPC_z43[4];
output->z45[0] = smoothMPC_z44[0];
output->z45[1] = smoothMPC_z44[1];
output->z45[2] = smoothMPC_z44[2];
output->z45[3] = smoothMPC_z44[3];
output->z45[4] = smoothMPC_z44[4];
output->z46[0] = smoothMPC_z45[0];
output->z46[1] = smoothMPC_z45[1];
output->z46[2] = smoothMPC_z45[2];
output->z46[3] = smoothMPC_z45[3];
output->z46[4] = smoothMPC_z45[4];
output->z47[0] = smoothMPC_z46[0];
output->z47[1] = smoothMPC_z46[1];
output->z47[2] = smoothMPC_z46[2];
output->z47[3] = smoothMPC_z46[3];
output->z47[4] = smoothMPC_z46[4];
output->z48[0] = smoothMPC_z47[0];
output->z48[1] = smoothMPC_z47[1];
output->z48[2] = smoothMPC_z47[2];
output->z48[3] = smoothMPC_z47[3];
output->z48[4] = smoothMPC_z47[4];
output->z49[0] = smoothMPC_z48[0];
output->z49[1] = smoothMPC_z48[1];
output->z49[2] = smoothMPC_z48[2];
output->z49[3] = smoothMPC_z48[3];
output->z49[4] = smoothMPC_z48[4];
output->z50[0] = smoothMPC_z49[0];
output->z50[1] = smoothMPC_z49[1];
output->z50[2] = smoothMPC_z49[2];
output->z50[3] = smoothMPC_z49[3];
output->z50[4] = smoothMPC_z49[4];
output->z51[0] = smoothMPC_z50[0];
output->z51[1] = smoothMPC_z50[1];
output->z51[2] = smoothMPC_z50[2];
output->z51[3] = smoothMPC_z50[3];
output->z51[4] = smoothMPC_z50[4];
output->z52[0] = smoothMPC_z51[0];
output->z52[1] = smoothMPC_z51[1];
output->z52[2] = smoothMPC_z51[2];
output->z52[3] = smoothMPC_z51[3];
output->z52[4] = smoothMPC_z51[4];
output->z53[0] = smoothMPC_z52[0];
output->z53[1] = smoothMPC_z52[1];
output->z53[2] = smoothMPC_z52[2];
output->z53[3] = smoothMPC_z52[3];
output->z53[4] = smoothMPC_z52[4];
output->z54[0] = smoothMPC_z53[0];
output->z54[1] = smoothMPC_z53[1];
output->z54[2] = smoothMPC_z53[2];
output->z54[3] = smoothMPC_z53[3];
output->z54[4] = smoothMPC_z53[4];
output->z55[0] = smoothMPC_z54[0];
output->z55[1] = smoothMPC_z54[1];
output->z55[2] = smoothMPC_z54[2];
output->z55[3] = smoothMPC_z54[3];
output->z55[4] = smoothMPC_z54[4];
output->z56[0] = smoothMPC_z55[0];
output->z56[1] = smoothMPC_z55[1];
output->z56[2] = smoothMPC_z55[2];
output->z56[3] = smoothMPC_z55[3];
output->z56[4] = smoothMPC_z55[4];
output->z57[0] = smoothMPC_z56[0];
output->z57[1] = smoothMPC_z56[1];
output->z57[2] = smoothMPC_z56[2];
output->z57[3] = smoothMPC_z56[3];
output->z57[4] = smoothMPC_z56[4];
output->z58[0] = smoothMPC_z57[0];
output->z58[1] = smoothMPC_z57[1];
output->z58[2] = smoothMPC_z57[2];
output->z58[3] = smoothMPC_z57[3];
output->z58[4] = smoothMPC_z57[4];
output->z59[0] = smoothMPC_z58[0];
output->z59[1] = smoothMPC_z58[1];
output->z59[2] = smoothMPC_z58[2];
output->z59[3] = smoothMPC_z58[3];
output->z59[4] = smoothMPC_z58[4];
output->z60[0] = smoothMPC_z59[0];
output->z60[1] = smoothMPC_z59[1];
output->z60[2] = smoothMPC_z59[2];
output->z60[3] = smoothMPC_z59[3];
output->z60[4] = smoothMPC_z59[4];
output->z61[0] = smoothMPC_z60[0];
output->z61[1] = smoothMPC_z60[1];
output->z61[2] = smoothMPC_z60[2];
output->z61[3] = smoothMPC_z60[3];
output->z61[4] = smoothMPC_z60[4];
output->z62[0] = smoothMPC_z61[0];
output->z62[1] = smoothMPC_z61[1];
output->z62[2] = smoothMPC_z61[2];
output->z62[3] = smoothMPC_z61[3];
output->z62[4] = smoothMPC_z61[4];
output->z63[0] = smoothMPC_z62[0];
output->z63[1] = smoothMPC_z62[1];
output->z63[2] = smoothMPC_z62[2];
output->z63[3] = smoothMPC_z62[3];
output->z63[4] = smoothMPC_z62[4];
output->z64[0] = smoothMPC_z63[0];
output->z64[1] = smoothMPC_z63[1];
output->z64[2] = smoothMPC_z63[2];
output->z64[3] = smoothMPC_z63[3];
output->z64[4] = smoothMPC_z63[4];
output->z65[0] = smoothMPC_z64[0];
output->z65[1] = smoothMPC_z64[1];
output->z65[2] = smoothMPC_z64[2];
output->z65[3] = smoothMPC_z64[3];
output->z65[4] = smoothMPC_z64[4];
output->z66[0] = smoothMPC_z65[0];
output->z66[1] = smoothMPC_z65[1];
output->z66[2] = smoothMPC_z65[2];
output->z66[3] = smoothMPC_z65[3];
output->z66[4] = smoothMPC_z65[4];
output->z67[0] = smoothMPC_z66[0];
output->z67[1] = smoothMPC_z66[1];
output->z67[2] = smoothMPC_z66[2];
output->z67[3] = smoothMPC_z66[3];
output->z67[4] = smoothMPC_z66[4];
output->z68[0] = smoothMPC_z67[0];
output->z68[1] = smoothMPC_z67[1];
output->z68[2] = smoothMPC_z67[2];
output->z68[3] = smoothMPC_z67[3];
output->z68[4] = smoothMPC_z67[4];
output->z69[0] = smoothMPC_z68[0];
output->z69[1] = smoothMPC_z68[1];
output->z69[2] = smoothMPC_z68[2];
output->z69[3] = smoothMPC_z68[3];
output->z69[4] = smoothMPC_z68[4];
output->z70[0] = smoothMPC_z69[0];
output->z70[1] = smoothMPC_z69[1];
output->z70[2] = smoothMPC_z69[2];
output->z70[3] = smoothMPC_z69[3];
output->z70[4] = smoothMPC_z69[4];
output->z71[0] = smoothMPC_z70[0];
output->z71[1] = smoothMPC_z70[1];
output->z71[2] = smoothMPC_z70[2];
output->z71[3] = smoothMPC_z70[3];
output->z71[4] = smoothMPC_z70[4];
output->z72[0] = smoothMPC_z71[0];
output->z72[1] = smoothMPC_z71[1];
output->z72[2] = smoothMPC_z71[2];
output->z72[3] = smoothMPC_z71[3];
output->z72[4] = smoothMPC_z71[4];
output->z73[0] = smoothMPC_z72[0];
output->z73[1] = smoothMPC_z72[1];
output->z73[2] = smoothMPC_z72[2];
output->z73[3] = smoothMPC_z72[3];
output->z73[4] = smoothMPC_z72[4];
output->z74[0] = smoothMPC_z73[0];
output->z74[1] = smoothMPC_z73[1];
output->z74[2] = smoothMPC_z73[2];
output->z74[3] = smoothMPC_z73[3];
output->z74[4] = smoothMPC_z73[4];
output->z75[0] = smoothMPC_z74[0];
output->z75[1] = smoothMPC_z74[1];
output->z75[2] = smoothMPC_z74[2];
output->z75[3] = smoothMPC_z74[3];
output->z75[4] = smoothMPC_z74[4];
output->z76[0] = smoothMPC_z75[0];
output->z76[1] = smoothMPC_z75[1];
output->z76[2] = smoothMPC_z75[2];
output->z76[3] = smoothMPC_z75[3];
output->z76[4] = smoothMPC_z75[4];
output->z77[0] = smoothMPC_z76[0];
output->z77[1] = smoothMPC_z76[1];
output->z77[2] = smoothMPC_z76[2];
output->z77[3] = smoothMPC_z76[3];
output->z77[4] = smoothMPC_z76[4];
output->z78[0] = smoothMPC_z77[0];
output->z78[1] = smoothMPC_z77[1];
output->z78[2] = smoothMPC_z77[2];
output->z78[3] = smoothMPC_z77[3];
output->z78[4] = smoothMPC_z77[4];
output->z79[0] = smoothMPC_z78[0];
output->z79[1] = smoothMPC_z78[1];
output->z79[2] = smoothMPC_z78[2];
output->z79[3] = smoothMPC_z78[3];
output->z79[4] = smoothMPC_z78[4];
output->z80[0] = smoothMPC_z79[0];
output->z80[1] = smoothMPC_z79[1];
output->z80[2] = smoothMPC_z79[2];
output->z80[3] = smoothMPC_z79[3];
output->z80[4] = smoothMPC_z79[4];
output->z81[0] = smoothMPC_z80[0];
output->z81[1] = smoothMPC_z80[1];
output->z81[2] = smoothMPC_z80[2];
output->z81[3] = smoothMPC_z80[3];
output->z81[4] = smoothMPC_z80[4];
output->z82[0] = smoothMPC_z81[0];
output->z82[1] = smoothMPC_z81[1];
output->z82[2] = smoothMPC_z81[2];
output->z82[3] = smoothMPC_z81[3];
output->z82[4] = smoothMPC_z81[4];
output->z83[0] = smoothMPC_z82[0];
output->z83[1] = smoothMPC_z82[1];
output->z83[2] = smoothMPC_z82[2];
output->z83[3] = smoothMPC_z82[3];
output->z83[4] = smoothMPC_z82[4];
output->z84[0] = smoothMPC_z83[0];
output->z84[1] = smoothMPC_z83[1];
output->z84[2] = smoothMPC_z83[2];
output->z84[3] = smoothMPC_z83[3];
output->z84[4] = smoothMPC_z83[4];
output->z85[0] = smoothMPC_z84[0];
output->z85[1] = smoothMPC_z84[1];
output->z85[2] = smoothMPC_z84[2];
output->z85[3] = smoothMPC_z84[3];
output->z85[4] = smoothMPC_z84[4];
output->z86[0] = smoothMPC_z85[0];
output->z86[1] = smoothMPC_z85[1];
output->z86[2] = smoothMPC_z85[2];
output->z86[3] = smoothMPC_z85[3];
output->z86[4] = smoothMPC_z85[4];
output->z87[0] = smoothMPC_z86[0];
output->z87[1] = smoothMPC_z86[1];
output->z87[2] = smoothMPC_z86[2];
output->z87[3] = smoothMPC_z86[3];
output->z87[4] = smoothMPC_z86[4];
output->z88[0] = smoothMPC_z87[0];
output->z88[1] = smoothMPC_z87[1];
output->z88[2] = smoothMPC_z87[2];
output->z88[3] = smoothMPC_z87[3];
output->z88[4] = smoothMPC_z87[4];
output->z89[0] = smoothMPC_z88[0];
output->z89[1] = smoothMPC_z88[1];
output->z89[2] = smoothMPC_z88[2];
output->z89[3] = smoothMPC_z88[3];
output->z89[4] = smoothMPC_z88[4];
output->z90[0] = smoothMPC_z89[0];
output->z90[1] = smoothMPC_z89[1];
output->z90[2] = smoothMPC_z89[2];
output->z90[3] = smoothMPC_z89[3];
output->z90[4] = smoothMPC_z89[4];
output->z91[0] = smoothMPC_z90[0];
output->z91[1] = smoothMPC_z90[1];
output->z91[2] = smoothMPC_z90[2];
output->z91[3] = smoothMPC_z90[3];
output->z91[4] = smoothMPC_z90[4];
output->z92[0] = smoothMPC_z91[0];
output->z92[1] = smoothMPC_z91[1];
output->z92[2] = smoothMPC_z91[2];
output->z92[3] = smoothMPC_z91[3];
output->z92[4] = smoothMPC_z91[4];
output->z93[0] = smoothMPC_z92[0];
output->z93[1] = smoothMPC_z92[1];
output->z93[2] = smoothMPC_z92[2];
output->z93[3] = smoothMPC_z92[3];
output->z93[4] = smoothMPC_z92[4];
output->z94[0] = smoothMPC_z93[0];
output->z94[1] = smoothMPC_z93[1];
output->z94[2] = smoothMPC_z93[2];
output->z94[3] = smoothMPC_z93[3];
output->z94[4] = smoothMPC_z93[4];
output->z95[0] = smoothMPC_z94[0];
output->z95[1] = smoothMPC_z94[1];
output->z95[2] = smoothMPC_z94[2];
output->z95[3] = smoothMPC_z94[3];
output->z95[4] = smoothMPC_z94[4];
output->z96[0] = smoothMPC_z95[0];
output->z96[1] = smoothMPC_z95[1];
output->z96[2] = smoothMPC_z95[2];
output->z96[3] = smoothMPC_z95[3];
output->z96[4] = smoothMPC_z95[4];
output->z97[0] = smoothMPC_z96[0];
output->z97[1] = smoothMPC_z96[1];
output->z97[2] = smoothMPC_z96[2];
output->z97[3] = smoothMPC_z96[3];
output->z97[4] = smoothMPC_z96[4];
output->z98[0] = smoothMPC_z97[0];
output->z98[1] = smoothMPC_z97[1];
output->z98[2] = smoothMPC_z97[2];
output->z98[3] = smoothMPC_z97[3];
output->z98[4] = smoothMPC_z97[4];
output->z99[0] = smoothMPC_z98[0];
output->z99[1] = smoothMPC_z98[1];
output->z99[2] = smoothMPC_z98[2];
output->z99[3] = smoothMPC_z98[3];
output->z99[4] = smoothMPC_z98[4];
output->z100[0] = smoothMPC_z99[0];
output->z100[1] = smoothMPC_z99[1];
output->z100[2] = smoothMPC_z99[2];

#if smoothMPC_SET_TIMING == 1
info->solvetime = smoothMPC_toc(&solvertimer);
#if smoothMPC_SET_PRINTLEVEL > 0 && smoothMPC_SET_TIMING == 1
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
