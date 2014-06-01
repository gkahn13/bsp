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

#include "entropyMPC.h"

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




/* LINEAR ALGEBRA LIBRARY ---------------------------------------------- */
/*
 * Initializes a vector of length 178 with a value.
 */
void entropyMPC_LA_INITIALIZEVECTOR_178(entropyMPC_FLOAT* vec, entropyMPC_FLOAT value)
{
	int i;
	for( i=0; i<178; i++ )
	{
		vec[i] = value;
	}
}


/*
 * Initializes a vector of length 150 with a value.
 */
void entropyMPC_LA_INITIALIZEVECTOR_150(entropyMPC_FLOAT* vec, entropyMPC_FLOAT value)
{
	int i;
	for( i=0; i<150; i++ )
	{
		vec[i] = value;
	}
}


/*
 * Initializes a vector of length 356 with a value.
 */
void entropyMPC_LA_INITIALIZEVECTOR_356(entropyMPC_FLOAT* vec, entropyMPC_FLOAT value)
{
	int i;
	for( i=0; i<356; i++ )
	{
		vec[i] = value;
	}
}


/* 
 * Calculates a dot product and adds it to a variable: z += x'*y; 
 * This function is for vectors of length 356.
 */
void entropyMPC_LA_DOTACC_356(entropyMPC_FLOAT *x, entropyMPC_FLOAT *y, entropyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<356; i++ ){
		*z += x[i]*y[i];
	}
}


/*
 * Calculates the gradient and the value for a quadratic function 0.5*z'*H*z + f'*z
 *
 * INPUTS:     H  - Symmetric Hessian, diag matrix of size [12 x 12]
 *             f  - column vector of size 12
 *             z  - column vector of size 12
 *
 * OUTPUTS: grad  - gradient at z (= H*z + f), column vector of size 12
 *          value <-- value + 0.5*z'*H*z + f'*z (value will be modified)
 */
void entropyMPC_LA_DIAG_QUADFCN_12(entropyMPC_FLOAT* H, entropyMPC_FLOAT* f, entropyMPC_FLOAT* z, entropyMPC_FLOAT* grad, entropyMPC_FLOAT* value)
{
	int i;
	entropyMPC_FLOAT hz;	
	for( i=0; i<12; i++){
		hz = H[i]*z[i];
		grad[i] = hz + f[i];
		*value += 0.5*hz*z[i] + f[i]*z[i];
	}
}


/*
 * Calculates the gradient and the value for a quadratic function 0.5*z'*H*z + f'*z
 *
 * INPUTS:     H  - Symmetric Hessian, diag matrix of size [10 x 10]
 *             f  - column vector of size 10
 *             z  - column vector of size 10
 *
 * OUTPUTS: grad  - gradient at z (= H*z + f), column vector of size 10
 *          value <-- value + 0.5*z'*H*z + f'*z (value will be modified)
 */
void entropyMPC_LA_DIAG_QUADFCN_10(entropyMPC_FLOAT* H, entropyMPC_FLOAT* f, entropyMPC_FLOAT* z, entropyMPC_FLOAT* grad, entropyMPC_FLOAT* value)
{
	int i;
	entropyMPC_FLOAT hz;	
	for( i=0; i<10; i++){
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
void entropyMPC_LA_DIAGZERO_MVMSUB6_10(entropyMPC_FLOAT *B, entropyMPC_FLOAT *u, entropyMPC_FLOAT *b, entropyMPC_FLOAT *l, entropyMPC_FLOAT *r, entropyMPC_FLOAT *z, entropyMPC_FLOAT *y)
{
	int i;
	entropyMPC_FLOAT Bu[10];
	entropyMPC_FLOAT norm = *y;
	entropyMPC_FLOAT lr = 0;

	/* do A*x + B*u first */
	for( i=0; i<10; i++ ){
		Bu[i] = B[i]*u[i];
	}	

	for( i=0; i<10; i++ ){
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
void entropyMPC_LA_DENSE_DIAGZERO_MVMSUB3_10_12_12(entropyMPC_FLOAT *A, entropyMPC_FLOAT *x, entropyMPC_FLOAT *B, entropyMPC_FLOAT *u, entropyMPC_FLOAT *b, entropyMPC_FLOAT *l, entropyMPC_FLOAT *r, entropyMPC_FLOAT *z, entropyMPC_FLOAT *y)
{
	int i;
	int j;
	int k = 0;
	entropyMPC_FLOAT AxBu[10];
	entropyMPC_FLOAT norm = *y;
	entropyMPC_FLOAT lr = 0;

	/* do A*x + B*u first */
	for( i=0; i<10; i++ ){
		AxBu[i] = A[k++]*x[0] + B[i]*u[i];
	}	

	for( j=1; j<12; j++ ){		
		for( i=0; i<10; i++ ){
			AxBu[i] += A[k++]*x[j];
		}
	}

	for( i=0; i<10; i++ ){
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
void entropyMPC_LA_DENSE_DIAGZERO_MVMSUB3_10_12_10(entropyMPC_FLOAT *A, entropyMPC_FLOAT *x, entropyMPC_FLOAT *B, entropyMPC_FLOAT *u, entropyMPC_FLOAT *b, entropyMPC_FLOAT *l, entropyMPC_FLOAT *r, entropyMPC_FLOAT *z, entropyMPC_FLOAT *y)
{
	int i;
	int j;
	int k = 0;
	entropyMPC_FLOAT AxBu[10];
	entropyMPC_FLOAT norm = *y;
	entropyMPC_FLOAT lr = 0;

	/* do A*x + B*u first */
	for( i=0; i<10; i++ ){
		AxBu[i] = A[k++]*x[0] + B[i]*u[i];
	}	

	for( j=1; j<12; j++ ){		
		for( i=0; i<10; i++ ){
			AxBu[i] += A[k++]*x[j];
		}
	}

	for( i=0; i<10; i++ ){
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
 * where A is of size [10 x 12] and stored in column major format.
 * and B is of size [10 x 12] and stored in diagzero format
 * Note the transposes of A and B!
 */
void entropyMPC_LA_DENSE_DIAGZERO_MTVM2_10_12_10(entropyMPC_FLOAT *A, entropyMPC_FLOAT *x, entropyMPC_FLOAT *B, entropyMPC_FLOAT *y, entropyMPC_FLOAT *z)
{
	int i;
	int j;
	int k = 0;
	for( i=0; i<10; i++ ){
		z[i] = 0;
		for( j=0; j<10; j++ ){
			z[i] += A[k++]*x[j];
		}
		z[i] += B[i]*y[i];
	}
	for( i=10 ;i<12; i++ ){
		z[i] = 0;
		for( j=0; j<10; j++ ){
			z[i] += A[k++]*x[j];
		}
	}
}


/*
 * Matrix vector multiplication y = M'*x where M is of size [10 x 10]
 * and stored in diagzero format. Note the transpose of M!
 */
void entropyMPC_LA_DIAGZERO_MTVM_10_10(entropyMPC_FLOAT *M, entropyMPC_FLOAT *x, entropyMPC_FLOAT *y)
{
	int i;
	for( i=0; i<10; i++ ){
		y[i] = M[i]*x[i];
	}
}


/*
 * Vector subtraction and addition.
 *	 Input: five vectors t, tidx, u, v, w and two scalars z and r
 *	 Output: y = t(tidx) - u + w
 *           z = z - v'*x;
 *           r = max([norm(y,inf), z]);
 * for vectors of length 12. Output z is of course scalar.
 */
void entropyMPC_LA_VSUBADD3_12(entropyMPC_FLOAT* t, entropyMPC_FLOAT* u, int* uidx, entropyMPC_FLOAT* v, entropyMPC_FLOAT* w, entropyMPC_FLOAT* y, entropyMPC_FLOAT* z, entropyMPC_FLOAT* r)
{
	int i;
	entropyMPC_FLOAT norm = *r;
	entropyMPC_FLOAT vx = 0;
	entropyMPC_FLOAT x;
	for( i=0; i<12; i++){
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
 * for vectors of length 12. Output z is of course scalar.
 */
void entropyMPC_LA_VSUBADD2_12(entropyMPC_FLOAT* t, int* tidx, entropyMPC_FLOAT* u, entropyMPC_FLOAT* v, entropyMPC_FLOAT* w, entropyMPC_FLOAT* y, entropyMPC_FLOAT* z, entropyMPC_FLOAT* r)
{
	int i;
	entropyMPC_FLOAT norm = *r;
	entropyMPC_FLOAT vx = 0;
	entropyMPC_FLOAT x;
	for( i=0; i<12; i++){
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
 * for vectors of length 10. Output z is of course scalar.
 */
void entropyMPC_LA_VSUBADD3_10(entropyMPC_FLOAT* t, entropyMPC_FLOAT* u, int* uidx, entropyMPC_FLOAT* v, entropyMPC_FLOAT* w, entropyMPC_FLOAT* y, entropyMPC_FLOAT* z, entropyMPC_FLOAT* r)
{
	int i;
	entropyMPC_FLOAT norm = *r;
	entropyMPC_FLOAT vx = 0;
	entropyMPC_FLOAT x;
	for( i=0; i<10; i++){
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
void entropyMPC_LA_VSUBADD2_10(entropyMPC_FLOAT* t, int* tidx, entropyMPC_FLOAT* u, entropyMPC_FLOAT* v, entropyMPC_FLOAT* w, entropyMPC_FLOAT* y, entropyMPC_FLOAT* z, entropyMPC_FLOAT* r)
{
	int i;
	entropyMPC_FLOAT norm = *r;
	entropyMPC_FLOAT vx = 0;
	entropyMPC_FLOAT x;
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
 * Computes inequality constraints gradient-
 * Special function for box constraints of length 12
 * Returns also L/S, a value that is often used elsewhere.
 */
void entropyMPC_LA_INEQ_B_GRAD_12_12_12(entropyMPC_FLOAT *lu, entropyMPC_FLOAT *su, entropyMPC_FLOAT *ru, entropyMPC_FLOAT *ll, entropyMPC_FLOAT *sl, entropyMPC_FLOAT *rl, int* lbIdx, int* ubIdx, entropyMPC_FLOAT *grad, entropyMPC_FLOAT *lubysu, entropyMPC_FLOAT *llbysl)
{
	int i;
	for( i=0; i<12; i++ ){
		grad[i] = 0;
	}
	for( i=0; i<12; i++ ){		
		llbysl[i] = ll[i] / sl[i];
		grad[lbIdx[i]] -= llbysl[i]*rl[i];
	}
	for( i=0; i<12; i++ ){
		lubysu[i] = lu[i] / su[i];
		grad[ubIdx[i]] += lubysu[i]*ru[i];
	}
}


/*
 * Computes inequality constraints gradient-
 * Special function for box constraints of length 10
 * Returns also L/S, a value that is often used elsewhere.
 */
void entropyMPC_LA_INEQ_B_GRAD_10_10_10(entropyMPC_FLOAT *lu, entropyMPC_FLOAT *su, entropyMPC_FLOAT *ru, entropyMPC_FLOAT *ll, entropyMPC_FLOAT *sl, entropyMPC_FLOAT *rl, int* lbIdx, int* ubIdx, entropyMPC_FLOAT *grad, entropyMPC_FLOAT *lubysu, entropyMPC_FLOAT *llbysl)
{
	int i;
	for( i=0; i<10; i++ ){
		grad[i] = 0;
	}
	for( i=0; i<10; i++ ){		
		llbysl[i] = ll[i] / sl[i];
		grad[lbIdx[i]] -= llbysl[i]*rl[i];
	}
	for( i=0; i<10; i++ ){
		lubysu[i] = lu[i] / su[i];
		grad[ubIdx[i]] += lubysu[i]*ru[i];
	}
}


/*
 * Addition of three vectors  z = u + w + v
 * of length 178.
 */
void entropyMPC_LA_VVADD3_178(entropyMPC_FLOAT *u, entropyMPC_FLOAT *v, entropyMPC_FLOAT *w, entropyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<178; i++ ){
		z[i] = u[i] + v[i] + w[i];
	}
}


/*
 * Special function to compute the diagonal cholesky factorization of the 
 * positive definite augmented Hessian for block size 12.
 *
 * Inputs: - H = diagonal cost Hessian in diagonal storage format
 *         - llbysl = L / S of lower bounds
 *         - lubysu = L / S of upper bounds
 *
 * Output: Phi = sqrt(H + diag(llbysl) + diag(lubysu))
 * where Phi is stored in diagonal storage format
 */
void entropyMPC_LA_DIAG_CHOL_ONELOOP_LBUB_12_12_12(entropyMPC_FLOAT *H, entropyMPC_FLOAT *llbysl, int* lbIdx, entropyMPC_FLOAT *lubysu, int* ubIdx, entropyMPC_FLOAT *Phi)


{
	int i;
	
	/* compute cholesky */
	for( i=0; i<12; i++ ){
		Phi[i] = H[i] + llbysl[i] + lubysu[i];

#if entropyMPC_SET_PRINTLEVEL > 0 && defined PRINTNUMERICALWARNINGS
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
 * where A is to be computed and is of size [10 x 12],
 * B is given and of size [10 x 12], L is a diagonal
 * matrix of size 10 stored in diagonal matrix 
 * storage format. Note the transpose of L has no impact!
 *
 * Result: A in column major storage format.
 *
 */
void entropyMPC_LA_DIAG_MATRIXFORWARDSUB_10_12(entropyMPC_FLOAT *L, entropyMPC_FLOAT *B, entropyMPC_FLOAT *A)
{
    int i,j;
	 int k = 0;

	for( j=0; j<12; j++){
		for( i=0; i<10; i++){
			A[k] = B[k]/L[j];
			k++;
		}
	}

}


/**
 * Forward substitution for the matrix equation A*L' = B
 * where A is to be computed and is of size [10 x 12],
 * B is given and of size [10 x 12], L is a diagonal
 *  matrix of size 12 stored in diagonal 
 * storage format. Note the transpose of L!
 *
 * Result: A in diagonalzero storage format.
 *
 */
void entropyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_10_12(entropyMPC_FLOAT *L, entropyMPC_FLOAT *B, entropyMPC_FLOAT *A)
{
	int j;
    for( j=0; j<12; j++ ){   
		A[j] = B[j]/L[j];
     }
}


/**
 * Compute C = A*B' where 
 *
 *	size(A) = [10 x 12]
 *  size(B) = [10 x 12] in diagzero format
 * 
 * A and C matrices are stored in column major format.
 * 
 * 
 */
void entropyMPC_LA_DENSE_DIAGZERO_MMTM_10_12_10(entropyMPC_FLOAT *A, entropyMPC_FLOAT *B, entropyMPC_FLOAT *C)
{
    int i, j;
	
	for( i=0; i<10; i++ ){
		for( j=0; j<10; j++){
			C[j*10+i] = B[i*10+j]*A[i];
		}
	}

}


/**
 * Forward substitution to solve L*y = b where L is a
 * diagonal matrix in vector storage format.
 * 
 * The dimensions involved are 12.
 */
void entropyMPC_LA_DIAG_FORWARDSUB_12(entropyMPC_FLOAT *L, entropyMPC_FLOAT *b, entropyMPC_FLOAT *y)
{
    int i;

    for( i=0; i<12; i++ ){
		y[i] = b[i]/L[i];
    }
}


/*
 * Special function to compute the diagonal cholesky factorization of the 
 * positive definite augmented Hessian for block size 10.
 *
 * Inputs: - H = diagonal cost Hessian in diagonal storage format
 *         - llbysl = L / S of lower bounds
 *         - lubysu = L / S of upper bounds
 *
 * Output: Phi = sqrt(H + diag(llbysl) + diag(lubysu))
 * where Phi is stored in diagonal storage format
 */
void entropyMPC_LA_DIAG_CHOL_ONELOOP_LBUB_10_10_10(entropyMPC_FLOAT *H, entropyMPC_FLOAT *llbysl, int* lbIdx, entropyMPC_FLOAT *lubysu, int* ubIdx, entropyMPC_FLOAT *Phi)


{
	int i;
	
	/* compute cholesky */
	for( i=0; i<10; i++ ){
		Phi[i] = H[i] + llbysl[i] + lubysu[i];

#if entropyMPC_SET_PRINTLEVEL > 0 && defined PRINTNUMERICALWARNINGS
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
 * where A is to be computed and is of size [10 x 10],
 * B is given and of size [10 x 10], L is a diagonal
 *  matrix of size 10 stored in diagonal 
 * storage format. Note the transpose of L!
 *
 * Result: A in diagonalzero storage format.
 *
 */
void entropyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_10_10(entropyMPC_FLOAT *L, entropyMPC_FLOAT *B, entropyMPC_FLOAT *A)
{
	int j;
    for( j=0; j<10; j++ ){   
		A[j] = B[j]/L[j];
     }
}


/**
 * Forward substitution to solve L*y = b where L is a
 * diagonal matrix in vector storage format.
 * 
 * The dimensions involved are 10.
 */
void entropyMPC_LA_DIAG_FORWARDSUB_10(entropyMPC_FLOAT *L, entropyMPC_FLOAT *b, entropyMPC_FLOAT *y)
{
    int i;

    for( i=0; i<10; i++ ){
		y[i] = b[i]/L[i];
    }
}


/**
 * Compute L = B*B', where L is lower triangular of size NXp1
 * and B is a diagzero matrix of size [10 x 12] in column
 * storage format.
 * 
 */
void entropyMPC_LA_DIAGZERO_MMT_10(entropyMPC_FLOAT *B, entropyMPC_FLOAT *L)
{
    int i, ii, di;
    
    ii = 0; di = 0;
    for( i=0; i<10; i++ ){        
		L[ii+i] = B[i]*B[i];
        ii += ++di;
    }
}


/* 
 * Computes r = b - B*u
 * B is stored in diagzero format
 */
void entropyMPC_LA_DIAGZERO_MVMSUB7_10(entropyMPC_FLOAT *B, entropyMPC_FLOAT *u, entropyMPC_FLOAT *b, entropyMPC_FLOAT *r)
{
	int i;

	for( i=0; i<10; i++ ){
		r[i] = b[i] - B[i]*u[i];
	}	
	
}


/**
 * Compute L = A*A' + B*B', where L is lower triangular of size NXp1
 * and A is a dense matrix of size [10 x 12] in column
 * storage format, and B is of size [10 x 12] diagonalzero
 * storage format.
 * 
 * THIS ONE HAS THE WORST ACCES PATTERN POSSIBLE. 
 * POSSIBKE FIX: PUT A AND B INTO ROW MAJOR FORMAT FIRST.
 * 
 */
void entropyMPC_LA_DENSE_DIAGZERO_MMT2_10_12_12(entropyMPC_FLOAT *A, entropyMPC_FLOAT *B, entropyMPC_FLOAT *L)
{
    int i, j, k, ii, di;
    entropyMPC_FLOAT ltemp;
    
    ii = 0; di = 0;
    for( i=0; i<10; i++ ){        
        for( j=0; j<=i; j++ ){
            ltemp = 0; 
            for( k=0; k<12; k++ ){
                ltemp += A[k*10+i]*A[k*10+j];
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
void entropyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_10_12_12(entropyMPC_FLOAT *A, entropyMPC_FLOAT *x, entropyMPC_FLOAT *B, entropyMPC_FLOAT *u, entropyMPC_FLOAT *b, entropyMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;

	for( i=0; i<10; i++ ){
		r[i] = b[i] - A[k++]*x[0] - B[i]*u[i];
	}	

	for( j=1; j<12; j++ ){		
		for( i=0; i<10; i++ ){
			r[i] -= A[k++]*x[j];
		}
	}
	
}


/**
 * Compute L = A*A' + B*B', where L is lower triangular of size NXp1
 * and A is a dense matrix of size [10 x 12] in column
 * storage format, and B is of size [10 x 10] diagonalzero
 * storage format.
 * 
 * THIS ONE HAS THE WORST ACCES PATTERN POSSIBLE. 
 * POSSIBKE FIX: PUT A AND B INTO ROW MAJOR FORMAT FIRST.
 * 
 */
void entropyMPC_LA_DENSE_DIAGZERO_MMT2_10_12_10(entropyMPC_FLOAT *A, entropyMPC_FLOAT *B, entropyMPC_FLOAT *L)
{
    int i, j, k, ii, di;
    entropyMPC_FLOAT ltemp;
    
    ii = 0; di = 0;
    for( i=0; i<10; i++ ){        
        for( j=0; j<=i; j++ ){
            ltemp = 0; 
            for( k=0; k<12; k++ ){
                ltemp += A[k*10+i]*A[k*10+j];
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
void entropyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_10_12_10(entropyMPC_FLOAT *A, entropyMPC_FLOAT *x, entropyMPC_FLOAT *B, entropyMPC_FLOAT *u, entropyMPC_FLOAT *b, entropyMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;

	for( i=0; i<10; i++ ){
		r[i] = b[i] - A[k++]*x[0] - B[i]*u[i];
	}	

	for( j=1; j<12; j++ ){		
		for( i=0; i<10; i++ ){
			r[i] -= A[k++]*x[j];
		}
	}
	
}


/**
 * Cholesky factorization as above, but working on a matrix in 
 * lower triangular storage format of size 10 and outputting
 * the Cholesky factor to matrix L in lower triangular format.
 */
void entropyMPC_LA_DENSE_CHOL_10(entropyMPC_FLOAT *A, entropyMPC_FLOAT *L)
{
    int i, j, k, di, dj;
	 int ii, jj;

    entropyMPC_FLOAT l;
    entropyMPC_FLOAT Mii;

	/* copy A to L first and then operate on L */
	/* COULD BE OPTIMIZED */
	ii=0; di=0;
	for( i=0; i<10; i++ ){
		for( j=0; j<=i; j++ ){
			L[ii+j] = A[ii+j];
		}
		ii += ++di;
	}    
	
	/* factor L */
	ii=0; di=0;
    for( i=0; i<10; i++ ){
        l = 0;
        for( k=0; k<i; k++ ){
            l += L[ii+k]*L[ii+k];
        }        
        
        Mii = L[ii+i] - l;
        
#if entropyMPC_SET_PRINTLEVEL > 0 && defined PRINTNUMERICALWARNINGS
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
        for( j=i+1; j<10; j++ ){
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
 * The dimensions involved are 10.
 */
void entropyMPC_LA_DENSE_FORWARDSUB_10(entropyMPC_FLOAT *L, entropyMPC_FLOAT *b, entropyMPC_FLOAT *y)
{
    int i,j,ii,di;
    entropyMPC_FLOAT yel;
            
    ii = 0; di = 0;
    for( i=0; i<10; i++ ){
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
 * where A is to be computed and is of size [10 x 10],
 * B is given and of size [10 x 10], L is a lower tri-
 * angular matrix of size 10 stored in lower triangular 
 * storage format. Note the transpose of L AND B!
 *
 * Result: A in column major storage format.
 *
 */
void entropyMPC_LA_DENSE_MATRIXTFORWARDSUB_10_10(entropyMPC_FLOAT *L, entropyMPC_FLOAT *B, entropyMPC_FLOAT *A)
{
    int i,j,k,ii,di;
    entropyMPC_FLOAT a;
    
    ii=0; di=0;
    for( j=0; j<10; j++ ){        
        for( i=0; i<10; i++ ){
            a = B[i*10+j];
            for( k=0; k<j; k++ ){
                a -= A[k*10+i]*L[ii+k];
            }    

			/* saturate for numerical stability */
			a = MIN(a, BIGM);
			a = MAX(a, -BIGM); 

			A[j*10+i] = a/L[ii+j];			
        }
        ii += ++di;
    }
}


/**
 * Compute L = L - A*A', where L is lower triangular of size 10
 * and A is a dense matrix of size [10 x 10] in column
 * storage format.
 * 
 * THIS ONE HAS THE WORST ACCES PATTERN POSSIBLE. 
 * POSSIBKE FIX: PUT A INTO ROW MAJOR FORMAT FIRST.
 * 
 */
void entropyMPC_LA_DENSE_MMTSUB_10_10(entropyMPC_FLOAT *A, entropyMPC_FLOAT *L)
{
    int i, j, k, ii, di;
    entropyMPC_FLOAT ltemp;
    
    ii = 0; di = 0;
    for( i=0; i<10; i++ ){        
        for( j=0; j<=i; j++ ){
            ltemp = 0; 
            for( k=0; k<10; k++ ){
                ltemp += A[k*10+i]*A[k*10+j];
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
void entropyMPC_LA_DENSE_MVMSUB1_10_10(entropyMPC_FLOAT *A, entropyMPC_FLOAT *x, entropyMPC_FLOAT *b, entropyMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;

	for( i=0; i<10; i++ ){
		r[i] = b[i] - A[k++]*x[0];
	}	
	for( j=1; j<10; j++ ){		
		for( i=0; i<10; i++ ){
			r[i] -= A[k++]*x[j];
		}
	}
}


/**
 * Backward Substitution to solve L^T*x = y where L is a
 * lower triangular matrix in triangular storage format.
 * 
 * All involved dimensions are 10.
 */
void entropyMPC_LA_DENSE_BACKWARDSUB_10(entropyMPC_FLOAT *L, entropyMPC_FLOAT *y, entropyMPC_FLOAT *x)
{
    int i, ii, di, j, jj, dj;
    entropyMPC_FLOAT xel;    
	int start = 45;
    
    /* now solve L^T*x = y by backward substitution */
    ii = start; di = 9;
    for( i=9; i>=0; i-- ){        
        xel = y[i];        
        jj = start; dj = 9;
        for( j=9; j>i; j-- ){
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
 * Matrix vector multiplication y = b - M'*x where M is of size [10 x 10]
 * and stored in column major format. Note the transpose of M!
 */
void entropyMPC_LA_DENSE_MTVMSUB_10_10(entropyMPC_FLOAT *A, entropyMPC_FLOAT *x, entropyMPC_FLOAT *b, entropyMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0; 
	for( i=0; i<10; i++ ){
		r[i] = b[i];
		for( j=0; j<10; j++ ){
			r[i] -= A[k++]*x[j];
		}
	}
}


/*
 * Vector subtraction z = -x - y for vectors of length 178.
 */
void entropyMPC_LA_VSUB2_178(entropyMPC_FLOAT *x, entropyMPC_FLOAT *y, entropyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<178; i++){
		z[i] = -x[i] - y[i];
	}
}


/**
 * Forward-Backward-Substitution to solve L*L^T*x = b where L is a
 * diagonal matrix of size 12 in vector
 * storage format.
 */
void entropyMPC_LA_DIAG_FORWARDBACKWARDSUB_12(entropyMPC_FLOAT *L, entropyMPC_FLOAT *b, entropyMPC_FLOAT *x)
{
    int i;
            
    /* solve Ly = b by forward and backward substitution */
    for( i=0; i<12; i++ ){
		x[i] = b[i]/(L[i]*L[i]);
    }
    
}


/**
 * Forward-Backward-Substitution to solve L*L^T*x = b where L is a
 * diagonal matrix of size 10 in vector
 * storage format.
 */
void entropyMPC_LA_DIAG_FORWARDBACKWARDSUB_10(entropyMPC_FLOAT *L, entropyMPC_FLOAT *b, entropyMPC_FLOAT *x)
{
    int i;
            
    /* solve Ly = b by forward and backward substitution */
    for( i=0; i<10; i++ ){
		x[i] = b[i]/(L[i]*L[i]);
    }
    
}


/*
 * Vector subtraction z = x(xidx) - y where y, z and xidx are of length 12,
 * and x has length 12 and is indexed through yidx.
 */
void entropyMPC_LA_VSUB_INDEXED_12(entropyMPC_FLOAT *x, int* xidx, entropyMPC_FLOAT *y, entropyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<12; i++){
		z[i] = x[xidx[i]] - y[i];
	}
}


/*
 * Vector subtraction x = -u.*v - w for vectors of length 12.
 */
void entropyMPC_LA_VSUB3_12(entropyMPC_FLOAT *u, entropyMPC_FLOAT *v, entropyMPC_FLOAT *w, entropyMPC_FLOAT *x)
{
	int i;
	for( i=0; i<12; i++){
		x[i] = -u[i]*v[i] - w[i];
	}
}


/*
 * Vector subtraction z = -x - y(yidx) where y is of length 12
 * and z, x and yidx are of length 12.
 */
void entropyMPC_LA_VSUB2_INDEXED_12(entropyMPC_FLOAT *x, entropyMPC_FLOAT *y, int* yidx, entropyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<12; i++){
		z[i] = -x[i] - y[yidx[i]];
	}
}


/*
 * Vector subtraction z = x(xidx) - y where y, z and xidx are of length 10,
 * and x has length 10 and is indexed through yidx.
 */
void entropyMPC_LA_VSUB_INDEXED_10(entropyMPC_FLOAT *x, int* xidx, entropyMPC_FLOAT *y, entropyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<10; i++){
		z[i] = x[xidx[i]] - y[i];
	}
}


/*
 * Vector subtraction x = -u.*v - w for vectors of length 10.
 */
void entropyMPC_LA_VSUB3_10(entropyMPC_FLOAT *u, entropyMPC_FLOAT *v, entropyMPC_FLOAT *w, entropyMPC_FLOAT *x)
{
	int i;
	for( i=0; i<10; i++){
		x[i] = -u[i]*v[i] - w[i];
	}
}


/*
 * Vector subtraction z = -x - y(yidx) where y is of length 10
 * and z, x and yidx are of length 10.
 */
void entropyMPC_LA_VSUB2_INDEXED_10(entropyMPC_FLOAT *x, entropyMPC_FLOAT *y, int* yidx, entropyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<10; i++){
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
 * entropyMPC_NOPROGRESS (should be negative).
 */
int entropyMPC_LINESEARCH_BACKTRACKING_AFFINE(entropyMPC_FLOAT *l, entropyMPC_FLOAT *s, entropyMPC_FLOAT *dl, entropyMPC_FLOAT *ds, entropyMPC_FLOAT *a, entropyMPC_FLOAT *mu_aff)
{
    int i;
	int lsIt=1;    
    entropyMPC_FLOAT dltemp;
    entropyMPC_FLOAT dstemp;
    entropyMPC_FLOAT mya = 1.0;
    entropyMPC_FLOAT mymu;
        
    while( 1 ){                        

        /* 
         * Compute both snew and wnew together.
         * We compute also mu_affine along the way here, as the
         * values might be in registers, so it should be cheaper.
         */
        mymu = 0;
        for( i=0; i<356; i++ ){
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
        if( i == 356 ){
            break;
        } else {
            mya *= entropyMPC_SET_LS_SCALE_AFF;
            if( mya < entropyMPC_SET_LS_MINSTEP ){
                return entropyMPC_NOPROGRESS;
            }
        }
    }
    
    /* return new values and iteration counter */
    *a = mya;
    *mu_aff = mymu / (entropyMPC_FLOAT)356;
    return lsIt;
}


/*
 * Vector subtraction x = u.*v - a where a is a scalar
*  and x,u,v are vectors of length 356.
 */
void entropyMPC_LA_VSUB5_356(entropyMPC_FLOAT *u, entropyMPC_FLOAT *v, entropyMPC_FLOAT a, entropyMPC_FLOAT *x)
{
	int i;
	for( i=0; i<356; i++){
		x[i] = u[i]*v[i] - a;
	}
}


/*
 * Computes x=0; x(uidx) += u/su; x(vidx) -= v/sv where x is of length 12,
 * u, su, uidx are of length 12 and v, sv, vidx are of length 12.
 */
void entropyMPC_LA_VSUB6_INDEXED_12_12_12(entropyMPC_FLOAT *u, entropyMPC_FLOAT *su, int* uidx, entropyMPC_FLOAT *v, entropyMPC_FLOAT *sv, int* vidx, entropyMPC_FLOAT *x)
{
	int i;
	for( i=0; i<12; i++ ){
		x[i] = 0;
	}
	for( i=0; i<12; i++){
		x[uidx[i]] += u[i]/su[i];
	}
	for( i=0; i<12; i++){
		x[vidx[i]] -= v[i]/sv[i];
	}
}


/* 
 * Computes r =  B*u
 * where B is stored in diagzero format
 */
void entropyMPC_LA_DIAGZERO_MVM_10(entropyMPC_FLOAT *B, entropyMPC_FLOAT *u, entropyMPC_FLOAT *r)
{
	int i;

	for( i=0; i<10; i++ ){
		r[i] = B[i]*u[i];
	}	
	
}


/* 
 * Computes r = A*x + B*u
 * where A is stored in column major format
 * and B is stored in diagzero format
 */
void entropyMPC_LA_DENSE_DIAGZERO_2MVMADD_10_12_12(entropyMPC_FLOAT *A, entropyMPC_FLOAT *x, entropyMPC_FLOAT *B, entropyMPC_FLOAT *u, entropyMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;

	for( i=0; i<10; i++ ){
		r[i] = A[k++]*x[0] + B[i]*u[i];
	}	

	for( j=1; j<12; j++ ){		
		for( i=0; i<10; i++ ){
			r[i] += A[k++]*x[j];
		}
	}
	
}


/*
 * Computes x=0; x(uidx) += u/su; x(vidx) -= v/sv where x is of length 10,
 * u, su, uidx are of length 10 and v, sv, vidx are of length 10.
 */
void entropyMPC_LA_VSUB6_INDEXED_10_10_10(entropyMPC_FLOAT *u, entropyMPC_FLOAT *su, int* uidx, entropyMPC_FLOAT *v, entropyMPC_FLOAT *sv, int* vidx, entropyMPC_FLOAT *x)
{
	int i;
	for( i=0; i<10; i++ ){
		x[i] = 0;
	}
	for( i=0; i<10; i++){
		x[uidx[i]] += u[i]/su[i];
	}
	for( i=0; i<10; i++){
		x[vidx[i]] -= v[i]/sv[i];
	}
}


/* 
 * Computes r = A*x + B*u
 * where A is stored in column major format
 * and B is stored in diagzero format
 */
void entropyMPC_LA_DENSE_DIAGZERO_2MVMADD_10_12_10(entropyMPC_FLOAT *A, entropyMPC_FLOAT *x, entropyMPC_FLOAT *B, entropyMPC_FLOAT *u, entropyMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;

	for( i=0; i<10; i++ ){
		r[i] = A[k++]*x[0] + B[i]*u[i];
	}	

	for( j=1; j<12; j++ ){		
		for( i=0; i<10; i++ ){
			r[i] += A[k++]*x[j];
		}
	}
	
}


/*
 * Vector subtraction z = x - y for vectors of length 178.
 */
void entropyMPC_LA_VSUB_178(entropyMPC_FLOAT *x, entropyMPC_FLOAT *y, entropyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<178; i++){
		z[i] = x[i] - y[i];
	}
}


/** 
 * Computes z = -r./s - u.*y(y)
 * where all vectors except of y are of length 12 (length of y >= 12).
 */
void entropyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_12(entropyMPC_FLOAT *r, entropyMPC_FLOAT *s, entropyMPC_FLOAT *u, entropyMPC_FLOAT *y, int* yidx, entropyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<12; i++ ){
		z[i] = -r[i]/s[i] - u[i]*y[yidx[i]];
	}
}


/** 
 * Computes z = -r./s + u.*y(y)
 * where all vectors except of y are of length 12 (length of y >= 12).
 */
void entropyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_12(entropyMPC_FLOAT *r, entropyMPC_FLOAT *s, entropyMPC_FLOAT *u, entropyMPC_FLOAT *y, int* yidx, entropyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<12; i++ ){
		z[i] = -r[i]/s[i] + u[i]*y[yidx[i]];
	}
}


/** 
 * Computes z = -r./s - u.*y(y)
 * where all vectors except of y are of length 10 (length of y >= 10).
 */
void entropyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_10(entropyMPC_FLOAT *r, entropyMPC_FLOAT *s, entropyMPC_FLOAT *u, entropyMPC_FLOAT *y, int* yidx, entropyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<10; i++ ){
		z[i] = -r[i]/s[i] - u[i]*y[yidx[i]];
	}
}


/** 
 * Computes z = -r./s + u.*y(y)
 * where all vectors except of y are of length 10 (length of y >= 10).
 */
void entropyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_10(entropyMPC_FLOAT *r, entropyMPC_FLOAT *s, entropyMPC_FLOAT *u, entropyMPC_FLOAT *y, int* yidx, entropyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<10; i++ ){
		z[i] = -r[i]/s[i] + u[i]*y[yidx[i]];
	}
}


/*
 * Computes ds = -l.\(r + s.*dl) for vectors of length 356.
 */
void entropyMPC_LA_VSUB7_356(entropyMPC_FLOAT *l, entropyMPC_FLOAT *r, entropyMPC_FLOAT *s, entropyMPC_FLOAT *dl, entropyMPC_FLOAT *ds)
{
	int i;
	for( i=0; i<356; i++){
		ds[i] = -(r[i] + s[i]*dl[i])/l[i];
	}
}


/*
 * Vector addition x = x + y for vectors of length 178.
 */
void entropyMPC_LA_VADD_178(entropyMPC_FLOAT *x, entropyMPC_FLOAT *y)
{
	int i;
	for( i=0; i<178; i++){
		x[i] += y[i];
	}
}


/*
 * Vector addition x = x + y for vectors of length 150.
 */
void entropyMPC_LA_VADD_150(entropyMPC_FLOAT *x, entropyMPC_FLOAT *y)
{
	int i;
	for( i=0; i<150; i++){
		x[i] += y[i];
	}
}


/*
 * Vector addition x = x + y for vectors of length 356.
 */
void entropyMPC_LA_VADD_356(entropyMPC_FLOAT *x, entropyMPC_FLOAT *y)
{
	int i;
	for( i=0; i<356; i++){
		x[i] += y[i];
	}
}


/**
 * Backtracking line search for combined predictor/corrector step.
 * Update on variables with safety factor gamma (to keep us away from
 * boundary).
 */
int entropyMPC_LINESEARCH_BACKTRACKING_COMBINED(entropyMPC_FLOAT *z, entropyMPC_FLOAT *v, entropyMPC_FLOAT *l, entropyMPC_FLOAT *s, entropyMPC_FLOAT *dz, entropyMPC_FLOAT *dv, entropyMPC_FLOAT *dl, entropyMPC_FLOAT *ds, entropyMPC_FLOAT *a, entropyMPC_FLOAT *mu)
{
    int i, lsIt=1;       
    entropyMPC_FLOAT dltemp;
    entropyMPC_FLOAT dstemp;    
    entropyMPC_FLOAT a_gamma;
            
    *a = 1.0;
    while( 1 ){                        

        /* check whether search criterion is fulfilled */
        for( i=0; i<356; i++ ){
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
        if( i == 356 ){
            break;
        } else {
            *a *= entropyMPC_SET_LS_SCALE;
            if( *a < entropyMPC_SET_LS_MINSTEP ){
                return entropyMPC_NOPROGRESS;
            }
        }
    }
    
    /* update variables with safety margin */
    a_gamma = (*a)*entropyMPC_SET_LS_MAXSTEP;
    
    /* primal variables */
    for( i=0; i<178; i++ ){
        z[i] += a_gamma*dz[i];
    }
    
    /* equality constraint multipliers */
    for( i=0; i<150; i++ ){
        v[i] += a_gamma*dv[i];
    }
    
    /* inequality constraint multipliers & slacks, also update mu */
    *mu = 0;
    for( i=0; i<356; i++ ){
        dltemp = l[i] + a_gamma*dl[i]; l[i] = dltemp;
        dstemp = s[i] + a_gamma*ds[i]; s[i] = dstemp;
        *mu += dltemp*dstemp;
    }
    
    *a = a_gamma;
    *mu /= (entropyMPC_FLOAT)356;
    return lsIt;
}




/* VARIABLE DEFINITIONS ------------------------------------------------ */
entropyMPC_FLOAT entropyMPC_z[178];
entropyMPC_FLOAT entropyMPC_v[150];
entropyMPC_FLOAT entropyMPC_dz_aff[178];
entropyMPC_FLOAT entropyMPC_dv_aff[150];
entropyMPC_FLOAT entropyMPC_grad_cost[178];
entropyMPC_FLOAT entropyMPC_grad_eq[178];
entropyMPC_FLOAT entropyMPC_rd[178];
entropyMPC_FLOAT entropyMPC_l[356];
entropyMPC_FLOAT entropyMPC_s[356];
entropyMPC_FLOAT entropyMPC_lbys[356];
entropyMPC_FLOAT entropyMPC_dl_aff[356];
entropyMPC_FLOAT entropyMPC_ds_aff[356];
entropyMPC_FLOAT entropyMPC_dz_cc[178];
entropyMPC_FLOAT entropyMPC_dv_cc[150];
entropyMPC_FLOAT entropyMPC_dl_cc[356];
entropyMPC_FLOAT entropyMPC_ds_cc[356];
entropyMPC_FLOAT entropyMPC_ccrhs[356];
entropyMPC_FLOAT entropyMPC_grad_ineq[178];
entropyMPC_FLOAT* entropyMPC_z00 = entropyMPC_z + 0;
entropyMPC_FLOAT* entropyMPC_dzaff00 = entropyMPC_dz_aff + 0;
entropyMPC_FLOAT* entropyMPC_dzcc00 = entropyMPC_dz_cc + 0;
entropyMPC_FLOAT* entropyMPC_rd00 = entropyMPC_rd + 0;
entropyMPC_FLOAT entropyMPC_Lbyrd00[12];
entropyMPC_FLOAT* entropyMPC_grad_cost00 = entropyMPC_grad_cost + 0;
entropyMPC_FLOAT* entropyMPC_grad_eq00 = entropyMPC_grad_eq + 0;
entropyMPC_FLOAT* entropyMPC_grad_ineq00 = entropyMPC_grad_ineq + 0;
entropyMPC_FLOAT entropyMPC_ctv00[12];
entropyMPC_FLOAT entropyMPC_C00[120] = {1.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 
0.0000000000000000E+000, 1.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 
0.0000000000000000E+000, 0.0000000000000000E+000, 1.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 
0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 1.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 
0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 1.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 
0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 1.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 
0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 1.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 
0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 1.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 
0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 1.0000000000000000E+000, 0.0000000000000000E+000, 
0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 1.0000000000000000E+000, 
1.0000000000000000E+000, 0.0000000000000000E+000, 1.0000000000000000E+000, 0.0000000000000000E+000, 1.0000000000000000E+000, 0.0000000000000000E+000, 1.0000000000000000E+000, 0.0000000000000000E+000, 1.0000000000000000E+000, 0.0000000000000000E+000, 
0.0000000000000000E+000, 1.0000000000000000E+000, 0.0000000000000000E+000, 1.0000000000000000E+000, 0.0000000000000000E+000, 1.0000000000000000E+000, 0.0000000000000000E+000, 1.0000000000000000E+000, 0.0000000000000000E+000, 1.0000000000000000E+000};
entropyMPC_FLOAT* entropyMPC_v00 = entropyMPC_v + 0;
entropyMPC_FLOAT entropyMPC_re00[10];
entropyMPC_FLOAT entropyMPC_beta00[10];
entropyMPC_FLOAT entropyMPC_betacc00[10];
entropyMPC_FLOAT* entropyMPC_dvaff00 = entropyMPC_dv_aff + 0;
entropyMPC_FLOAT* entropyMPC_dvcc00 = entropyMPC_dv_cc + 0;
entropyMPC_FLOAT entropyMPC_V00[120];
entropyMPC_FLOAT entropyMPC_Yd00[55];
entropyMPC_FLOAT entropyMPC_Ld00[55];
entropyMPC_FLOAT entropyMPC_yy00[10];
entropyMPC_FLOAT entropyMPC_bmy00[10];
int entropyMPC_lbIdx00[12] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
entropyMPC_FLOAT* entropyMPC_llb00 = entropyMPC_l + 0;
entropyMPC_FLOAT* entropyMPC_slb00 = entropyMPC_s + 0;
entropyMPC_FLOAT* entropyMPC_llbbyslb00 = entropyMPC_lbys + 0;
entropyMPC_FLOAT entropyMPC_rilb00[12];
entropyMPC_FLOAT* entropyMPC_dllbaff00 = entropyMPC_dl_aff + 0;
entropyMPC_FLOAT* entropyMPC_dslbaff00 = entropyMPC_ds_aff + 0;
entropyMPC_FLOAT* entropyMPC_dllbcc00 = entropyMPC_dl_cc + 0;
entropyMPC_FLOAT* entropyMPC_dslbcc00 = entropyMPC_ds_cc + 0;
entropyMPC_FLOAT* entropyMPC_ccrhsl00 = entropyMPC_ccrhs + 0;
int entropyMPC_ubIdx00[12] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
entropyMPC_FLOAT* entropyMPC_lub00 = entropyMPC_l + 12;
entropyMPC_FLOAT* entropyMPC_sub00 = entropyMPC_s + 12;
entropyMPC_FLOAT* entropyMPC_lubbysub00 = entropyMPC_lbys + 12;
entropyMPC_FLOAT entropyMPC_riub00[12];
entropyMPC_FLOAT* entropyMPC_dlubaff00 = entropyMPC_dl_aff + 12;
entropyMPC_FLOAT* entropyMPC_dsubaff00 = entropyMPC_ds_aff + 12;
entropyMPC_FLOAT* entropyMPC_dlubcc00 = entropyMPC_dl_cc + 12;
entropyMPC_FLOAT* entropyMPC_dsubcc00 = entropyMPC_ds_cc + 12;
entropyMPC_FLOAT* entropyMPC_ccrhsub00 = entropyMPC_ccrhs + 12;
entropyMPC_FLOAT entropyMPC_Phi00[12];
entropyMPC_FLOAT entropyMPC_D00[12] = {1.0000000000000000E+000, 
1.0000000000000000E+000, 
1.0000000000000000E+000, 
1.0000000000000000E+000, 
1.0000000000000000E+000, 
1.0000000000000000E+000, 
1.0000000000000000E+000, 
1.0000000000000000E+000, 
1.0000000000000000E+000, 
1.0000000000000000E+000};
entropyMPC_FLOAT entropyMPC_W00[12];
entropyMPC_FLOAT* entropyMPC_z01 = entropyMPC_z + 12;
entropyMPC_FLOAT* entropyMPC_dzaff01 = entropyMPC_dz_aff + 12;
entropyMPC_FLOAT* entropyMPC_dzcc01 = entropyMPC_dz_cc + 12;
entropyMPC_FLOAT* entropyMPC_rd01 = entropyMPC_rd + 12;
entropyMPC_FLOAT entropyMPC_Lbyrd01[12];
entropyMPC_FLOAT* entropyMPC_grad_cost01 = entropyMPC_grad_cost + 12;
entropyMPC_FLOAT* entropyMPC_grad_eq01 = entropyMPC_grad_eq + 12;
entropyMPC_FLOAT* entropyMPC_grad_ineq01 = entropyMPC_grad_ineq + 12;
entropyMPC_FLOAT entropyMPC_ctv01[12];
entropyMPC_FLOAT* entropyMPC_v01 = entropyMPC_v + 10;
entropyMPC_FLOAT entropyMPC_re01[10];
entropyMPC_FLOAT entropyMPC_beta01[10];
entropyMPC_FLOAT entropyMPC_betacc01[10];
entropyMPC_FLOAT* entropyMPC_dvaff01 = entropyMPC_dv_aff + 10;
entropyMPC_FLOAT* entropyMPC_dvcc01 = entropyMPC_dv_cc + 10;
entropyMPC_FLOAT entropyMPC_V01[120];
entropyMPC_FLOAT entropyMPC_Yd01[55];
entropyMPC_FLOAT entropyMPC_Ld01[55];
entropyMPC_FLOAT entropyMPC_yy01[10];
entropyMPC_FLOAT entropyMPC_bmy01[10];
entropyMPC_FLOAT entropyMPC_c01[10] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
int entropyMPC_lbIdx01[12] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
entropyMPC_FLOAT* entropyMPC_llb01 = entropyMPC_l + 24;
entropyMPC_FLOAT* entropyMPC_slb01 = entropyMPC_s + 24;
entropyMPC_FLOAT* entropyMPC_llbbyslb01 = entropyMPC_lbys + 24;
entropyMPC_FLOAT entropyMPC_rilb01[12];
entropyMPC_FLOAT* entropyMPC_dllbaff01 = entropyMPC_dl_aff + 24;
entropyMPC_FLOAT* entropyMPC_dslbaff01 = entropyMPC_ds_aff + 24;
entropyMPC_FLOAT* entropyMPC_dllbcc01 = entropyMPC_dl_cc + 24;
entropyMPC_FLOAT* entropyMPC_dslbcc01 = entropyMPC_ds_cc + 24;
entropyMPC_FLOAT* entropyMPC_ccrhsl01 = entropyMPC_ccrhs + 24;
int entropyMPC_ubIdx01[12] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
entropyMPC_FLOAT* entropyMPC_lub01 = entropyMPC_l + 36;
entropyMPC_FLOAT* entropyMPC_sub01 = entropyMPC_s + 36;
entropyMPC_FLOAT* entropyMPC_lubbysub01 = entropyMPC_lbys + 36;
entropyMPC_FLOAT entropyMPC_riub01[12];
entropyMPC_FLOAT* entropyMPC_dlubaff01 = entropyMPC_dl_aff + 36;
entropyMPC_FLOAT* entropyMPC_dsubaff01 = entropyMPC_ds_aff + 36;
entropyMPC_FLOAT* entropyMPC_dlubcc01 = entropyMPC_dl_cc + 36;
entropyMPC_FLOAT* entropyMPC_dsubcc01 = entropyMPC_ds_cc + 36;
entropyMPC_FLOAT* entropyMPC_ccrhsub01 = entropyMPC_ccrhs + 36;
entropyMPC_FLOAT entropyMPC_Phi01[12];
entropyMPC_FLOAT entropyMPC_D01[12] = {-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000};
entropyMPC_FLOAT entropyMPC_W01[12];
entropyMPC_FLOAT entropyMPC_Ysd01[100];
entropyMPC_FLOAT entropyMPC_Lsd01[100];
entropyMPC_FLOAT* entropyMPC_z02 = entropyMPC_z + 24;
entropyMPC_FLOAT* entropyMPC_dzaff02 = entropyMPC_dz_aff + 24;
entropyMPC_FLOAT* entropyMPC_dzcc02 = entropyMPC_dz_cc + 24;
entropyMPC_FLOAT* entropyMPC_rd02 = entropyMPC_rd + 24;
entropyMPC_FLOAT entropyMPC_Lbyrd02[12];
entropyMPC_FLOAT* entropyMPC_grad_cost02 = entropyMPC_grad_cost + 24;
entropyMPC_FLOAT* entropyMPC_grad_eq02 = entropyMPC_grad_eq + 24;
entropyMPC_FLOAT* entropyMPC_grad_ineq02 = entropyMPC_grad_ineq + 24;
entropyMPC_FLOAT entropyMPC_ctv02[12];
entropyMPC_FLOAT* entropyMPC_v02 = entropyMPC_v + 20;
entropyMPC_FLOAT entropyMPC_re02[10];
entropyMPC_FLOAT entropyMPC_beta02[10];
entropyMPC_FLOAT entropyMPC_betacc02[10];
entropyMPC_FLOAT* entropyMPC_dvaff02 = entropyMPC_dv_aff + 20;
entropyMPC_FLOAT* entropyMPC_dvcc02 = entropyMPC_dv_cc + 20;
entropyMPC_FLOAT entropyMPC_V02[120];
entropyMPC_FLOAT entropyMPC_Yd02[55];
entropyMPC_FLOAT entropyMPC_Ld02[55];
entropyMPC_FLOAT entropyMPC_yy02[10];
entropyMPC_FLOAT entropyMPC_bmy02[10];
entropyMPC_FLOAT entropyMPC_c02[10] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
int entropyMPC_lbIdx02[12] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
entropyMPC_FLOAT* entropyMPC_llb02 = entropyMPC_l + 48;
entropyMPC_FLOAT* entropyMPC_slb02 = entropyMPC_s + 48;
entropyMPC_FLOAT* entropyMPC_llbbyslb02 = entropyMPC_lbys + 48;
entropyMPC_FLOAT entropyMPC_rilb02[12];
entropyMPC_FLOAT* entropyMPC_dllbaff02 = entropyMPC_dl_aff + 48;
entropyMPC_FLOAT* entropyMPC_dslbaff02 = entropyMPC_ds_aff + 48;
entropyMPC_FLOAT* entropyMPC_dllbcc02 = entropyMPC_dl_cc + 48;
entropyMPC_FLOAT* entropyMPC_dslbcc02 = entropyMPC_ds_cc + 48;
entropyMPC_FLOAT* entropyMPC_ccrhsl02 = entropyMPC_ccrhs + 48;
int entropyMPC_ubIdx02[12] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
entropyMPC_FLOAT* entropyMPC_lub02 = entropyMPC_l + 60;
entropyMPC_FLOAT* entropyMPC_sub02 = entropyMPC_s + 60;
entropyMPC_FLOAT* entropyMPC_lubbysub02 = entropyMPC_lbys + 60;
entropyMPC_FLOAT entropyMPC_riub02[12];
entropyMPC_FLOAT* entropyMPC_dlubaff02 = entropyMPC_dl_aff + 60;
entropyMPC_FLOAT* entropyMPC_dsubaff02 = entropyMPC_ds_aff + 60;
entropyMPC_FLOAT* entropyMPC_dlubcc02 = entropyMPC_dl_cc + 60;
entropyMPC_FLOAT* entropyMPC_dsubcc02 = entropyMPC_ds_cc + 60;
entropyMPC_FLOAT* entropyMPC_ccrhsub02 = entropyMPC_ccrhs + 60;
entropyMPC_FLOAT entropyMPC_Phi02[12];
entropyMPC_FLOAT entropyMPC_W02[12];
entropyMPC_FLOAT entropyMPC_Ysd02[100];
entropyMPC_FLOAT entropyMPC_Lsd02[100];
entropyMPC_FLOAT* entropyMPC_z03 = entropyMPC_z + 36;
entropyMPC_FLOAT* entropyMPC_dzaff03 = entropyMPC_dz_aff + 36;
entropyMPC_FLOAT* entropyMPC_dzcc03 = entropyMPC_dz_cc + 36;
entropyMPC_FLOAT* entropyMPC_rd03 = entropyMPC_rd + 36;
entropyMPC_FLOAT entropyMPC_Lbyrd03[12];
entropyMPC_FLOAT* entropyMPC_grad_cost03 = entropyMPC_grad_cost + 36;
entropyMPC_FLOAT* entropyMPC_grad_eq03 = entropyMPC_grad_eq + 36;
entropyMPC_FLOAT* entropyMPC_grad_ineq03 = entropyMPC_grad_ineq + 36;
entropyMPC_FLOAT entropyMPC_ctv03[12];
entropyMPC_FLOAT* entropyMPC_v03 = entropyMPC_v + 30;
entropyMPC_FLOAT entropyMPC_re03[10];
entropyMPC_FLOAT entropyMPC_beta03[10];
entropyMPC_FLOAT entropyMPC_betacc03[10];
entropyMPC_FLOAT* entropyMPC_dvaff03 = entropyMPC_dv_aff + 30;
entropyMPC_FLOAT* entropyMPC_dvcc03 = entropyMPC_dv_cc + 30;
entropyMPC_FLOAT entropyMPC_V03[120];
entropyMPC_FLOAT entropyMPC_Yd03[55];
entropyMPC_FLOAT entropyMPC_Ld03[55];
entropyMPC_FLOAT entropyMPC_yy03[10];
entropyMPC_FLOAT entropyMPC_bmy03[10];
entropyMPC_FLOAT entropyMPC_c03[10] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
int entropyMPC_lbIdx03[12] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
entropyMPC_FLOAT* entropyMPC_llb03 = entropyMPC_l + 72;
entropyMPC_FLOAT* entropyMPC_slb03 = entropyMPC_s + 72;
entropyMPC_FLOAT* entropyMPC_llbbyslb03 = entropyMPC_lbys + 72;
entropyMPC_FLOAT entropyMPC_rilb03[12];
entropyMPC_FLOAT* entropyMPC_dllbaff03 = entropyMPC_dl_aff + 72;
entropyMPC_FLOAT* entropyMPC_dslbaff03 = entropyMPC_ds_aff + 72;
entropyMPC_FLOAT* entropyMPC_dllbcc03 = entropyMPC_dl_cc + 72;
entropyMPC_FLOAT* entropyMPC_dslbcc03 = entropyMPC_ds_cc + 72;
entropyMPC_FLOAT* entropyMPC_ccrhsl03 = entropyMPC_ccrhs + 72;
int entropyMPC_ubIdx03[12] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
entropyMPC_FLOAT* entropyMPC_lub03 = entropyMPC_l + 84;
entropyMPC_FLOAT* entropyMPC_sub03 = entropyMPC_s + 84;
entropyMPC_FLOAT* entropyMPC_lubbysub03 = entropyMPC_lbys + 84;
entropyMPC_FLOAT entropyMPC_riub03[12];
entropyMPC_FLOAT* entropyMPC_dlubaff03 = entropyMPC_dl_aff + 84;
entropyMPC_FLOAT* entropyMPC_dsubaff03 = entropyMPC_ds_aff + 84;
entropyMPC_FLOAT* entropyMPC_dlubcc03 = entropyMPC_dl_cc + 84;
entropyMPC_FLOAT* entropyMPC_dsubcc03 = entropyMPC_ds_cc + 84;
entropyMPC_FLOAT* entropyMPC_ccrhsub03 = entropyMPC_ccrhs + 84;
entropyMPC_FLOAT entropyMPC_Phi03[12];
entropyMPC_FLOAT entropyMPC_W03[12];
entropyMPC_FLOAT entropyMPC_Ysd03[100];
entropyMPC_FLOAT entropyMPC_Lsd03[100];
entropyMPC_FLOAT* entropyMPC_z04 = entropyMPC_z + 48;
entropyMPC_FLOAT* entropyMPC_dzaff04 = entropyMPC_dz_aff + 48;
entropyMPC_FLOAT* entropyMPC_dzcc04 = entropyMPC_dz_cc + 48;
entropyMPC_FLOAT* entropyMPC_rd04 = entropyMPC_rd + 48;
entropyMPC_FLOAT entropyMPC_Lbyrd04[12];
entropyMPC_FLOAT* entropyMPC_grad_cost04 = entropyMPC_grad_cost + 48;
entropyMPC_FLOAT* entropyMPC_grad_eq04 = entropyMPC_grad_eq + 48;
entropyMPC_FLOAT* entropyMPC_grad_ineq04 = entropyMPC_grad_ineq + 48;
entropyMPC_FLOAT entropyMPC_ctv04[12];
entropyMPC_FLOAT* entropyMPC_v04 = entropyMPC_v + 40;
entropyMPC_FLOAT entropyMPC_re04[10];
entropyMPC_FLOAT entropyMPC_beta04[10];
entropyMPC_FLOAT entropyMPC_betacc04[10];
entropyMPC_FLOAT* entropyMPC_dvaff04 = entropyMPC_dv_aff + 40;
entropyMPC_FLOAT* entropyMPC_dvcc04 = entropyMPC_dv_cc + 40;
entropyMPC_FLOAT entropyMPC_V04[120];
entropyMPC_FLOAT entropyMPC_Yd04[55];
entropyMPC_FLOAT entropyMPC_Ld04[55];
entropyMPC_FLOAT entropyMPC_yy04[10];
entropyMPC_FLOAT entropyMPC_bmy04[10];
entropyMPC_FLOAT entropyMPC_c04[10] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
int entropyMPC_lbIdx04[12] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
entropyMPC_FLOAT* entropyMPC_llb04 = entropyMPC_l + 96;
entropyMPC_FLOAT* entropyMPC_slb04 = entropyMPC_s + 96;
entropyMPC_FLOAT* entropyMPC_llbbyslb04 = entropyMPC_lbys + 96;
entropyMPC_FLOAT entropyMPC_rilb04[12];
entropyMPC_FLOAT* entropyMPC_dllbaff04 = entropyMPC_dl_aff + 96;
entropyMPC_FLOAT* entropyMPC_dslbaff04 = entropyMPC_ds_aff + 96;
entropyMPC_FLOAT* entropyMPC_dllbcc04 = entropyMPC_dl_cc + 96;
entropyMPC_FLOAT* entropyMPC_dslbcc04 = entropyMPC_ds_cc + 96;
entropyMPC_FLOAT* entropyMPC_ccrhsl04 = entropyMPC_ccrhs + 96;
int entropyMPC_ubIdx04[12] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
entropyMPC_FLOAT* entropyMPC_lub04 = entropyMPC_l + 108;
entropyMPC_FLOAT* entropyMPC_sub04 = entropyMPC_s + 108;
entropyMPC_FLOAT* entropyMPC_lubbysub04 = entropyMPC_lbys + 108;
entropyMPC_FLOAT entropyMPC_riub04[12];
entropyMPC_FLOAT* entropyMPC_dlubaff04 = entropyMPC_dl_aff + 108;
entropyMPC_FLOAT* entropyMPC_dsubaff04 = entropyMPC_ds_aff + 108;
entropyMPC_FLOAT* entropyMPC_dlubcc04 = entropyMPC_dl_cc + 108;
entropyMPC_FLOAT* entropyMPC_dsubcc04 = entropyMPC_ds_cc + 108;
entropyMPC_FLOAT* entropyMPC_ccrhsub04 = entropyMPC_ccrhs + 108;
entropyMPC_FLOAT entropyMPC_Phi04[12];
entropyMPC_FLOAT entropyMPC_W04[12];
entropyMPC_FLOAT entropyMPC_Ysd04[100];
entropyMPC_FLOAT entropyMPC_Lsd04[100];
entropyMPC_FLOAT* entropyMPC_z05 = entropyMPC_z + 60;
entropyMPC_FLOAT* entropyMPC_dzaff05 = entropyMPC_dz_aff + 60;
entropyMPC_FLOAT* entropyMPC_dzcc05 = entropyMPC_dz_cc + 60;
entropyMPC_FLOAT* entropyMPC_rd05 = entropyMPC_rd + 60;
entropyMPC_FLOAT entropyMPC_Lbyrd05[12];
entropyMPC_FLOAT* entropyMPC_grad_cost05 = entropyMPC_grad_cost + 60;
entropyMPC_FLOAT* entropyMPC_grad_eq05 = entropyMPC_grad_eq + 60;
entropyMPC_FLOAT* entropyMPC_grad_ineq05 = entropyMPC_grad_ineq + 60;
entropyMPC_FLOAT entropyMPC_ctv05[12];
entropyMPC_FLOAT* entropyMPC_v05 = entropyMPC_v + 50;
entropyMPC_FLOAT entropyMPC_re05[10];
entropyMPC_FLOAT entropyMPC_beta05[10];
entropyMPC_FLOAT entropyMPC_betacc05[10];
entropyMPC_FLOAT* entropyMPC_dvaff05 = entropyMPC_dv_aff + 50;
entropyMPC_FLOAT* entropyMPC_dvcc05 = entropyMPC_dv_cc + 50;
entropyMPC_FLOAT entropyMPC_V05[120];
entropyMPC_FLOAT entropyMPC_Yd05[55];
entropyMPC_FLOAT entropyMPC_Ld05[55];
entropyMPC_FLOAT entropyMPC_yy05[10];
entropyMPC_FLOAT entropyMPC_bmy05[10];
entropyMPC_FLOAT entropyMPC_c05[10] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
int entropyMPC_lbIdx05[12] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
entropyMPC_FLOAT* entropyMPC_llb05 = entropyMPC_l + 120;
entropyMPC_FLOAT* entropyMPC_slb05 = entropyMPC_s + 120;
entropyMPC_FLOAT* entropyMPC_llbbyslb05 = entropyMPC_lbys + 120;
entropyMPC_FLOAT entropyMPC_rilb05[12];
entropyMPC_FLOAT* entropyMPC_dllbaff05 = entropyMPC_dl_aff + 120;
entropyMPC_FLOAT* entropyMPC_dslbaff05 = entropyMPC_ds_aff + 120;
entropyMPC_FLOAT* entropyMPC_dllbcc05 = entropyMPC_dl_cc + 120;
entropyMPC_FLOAT* entropyMPC_dslbcc05 = entropyMPC_ds_cc + 120;
entropyMPC_FLOAT* entropyMPC_ccrhsl05 = entropyMPC_ccrhs + 120;
int entropyMPC_ubIdx05[12] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
entropyMPC_FLOAT* entropyMPC_lub05 = entropyMPC_l + 132;
entropyMPC_FLOAT* entropyMPC_sub05 = entropyMPC_s + 132;
entropyMPC_FLOAT* entropyMPC_lubbysub05 = entropyMPC_lbys + 132;
entropyMPC_FLOAT entropyMPC_riub05[12];
entropyMPC_FLOAT* entropyMPC_dlubaff05 = entropyMPC_dl_aff + 132;
entropyMPC_FLOAT* entropyMPC_dsubaff05 = entropyMPC_ds_aff + 132;
entropyMPC_FLOAT* entropyMPC_dlubcc05 = entropyMPC_dl_cc + 132;
entropyMPC_FLOAT* entropyMPC_dsubcc05 = entropyMPC_ds_cc + 132;
entropyMPC_FLOAT* entropyMPC_ccrhsub05 = entropyMPC_ccrhs + 132;
entropyMPC_FLOAT entropyMPC_Phi05[12];
entropyMPC_FLOAT entropyMPC_W05[12];
entropyMPC_FLOAT entropyMPC_Ysd05[100];
entropyMPC_FLOAT entropyMPC_Lsd05[100];
entropyMPC_FLOAT* entropyMPC_z06 = entropyMPC_z + 72;
entropyMPC_FLOAT* entropyMPC_dzaff06 = entropyMPC_dz_aff + 72;
entropyMPC_FLOAT* entropyMPC_dzcc06 = entropyMPC_dz_cc + 72;
entropyMPC_FLOAT* entropyMPC_rd06 = entropyMPC_rd + 72;
entropyMPC_FLOAT entropyMPC_Lbyrd06[12];
entropyMPC_FLOAT* entropyMPC_grad_cost06 = entropyMPC_grad_cost + 72;
entropyMPC_FLOAT* entropyMPC_grad_eq06 = entropyMPC_grad_eq + 72;
entropyMPC_FLOAT* entropyMPC_grad_ineq06 = entropyMPC_grad_ineq + 72;
entropyMPC_FLOAT entropyMPC_ctv06[12];
entropyMPC_FLOAT* entropyMPC_v06 = entropyMPC_v + 60;
entropyMPC_FLOAT entropyMPC_re06[10];
entropyMPC_FLOAT entropyMPC_beta06[10];
entropyMPC_FLOAT entropyMPC_betacc06[10];
entropyMPC_FLOAT* entropyMPC_dvaff06 = entropyMPC_dv_aff + 60;
entropyMPC_FLOAT* entropyMPC_dvcc06 = entropyMPC_dv_cc + 60;
entropyMPC_FLOAT entropyMPC_V06[120];
entropyMPC_FLOAT entropyMPC_Yd06[55];
entropyMPC_FLOAT entropyMPC_Ld06[55];
entropyMPC_FLOAT entropyMPC_yy06[10];
entropyMPC_FLOAT entropyMPC_bmy06[10];
entropyMPC_FLOAT entropyMPC_c06[10] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
int entropyMPC_lbIdx06[12] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
entropyMPC_FLOAT* entropyMPC_llb06 = entropyMPC_l + 144;
entropyMPC_FLOAT* entropyMPC_slb06 = entropyMPC_s + 144;
entropyMPC_FLOAT* entropyMPC_llbbyslb06 = entropyMPC_lbys + 144;
entropyMPC_FLOAT entropyMPC_rilb06[12];
entropyMPC_FLOAT* entropyMPC_dllbaff06 = entropyMPC_dl_aff + 144;
entropyMPC_FLOAT* entropyMPC_dslbaff06 = entropyMPC_ds_aff + 144;
entropyMPC_FLOAT* entropyMPC_dllbcc06 = entropyMPC_dl_cc + 144;
entropyMPC_FLOAT* entropyMPC_dslbcc06 = entropyMPC_ds_cc + 144;
entropyMPC_FLOAT* entropyMPC_ccrhsl06 = entropyMPC_ccrhs + 144;
int entropyMPC_ubIdx06[12] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
entropyMPC_FLOAT* entropyMPC_lub06 = entropyMPC_l + 156;
entropyMPC_FLOAT* entropyMPC_sub06 = entropyMPC_s + 156;
entropyMPC_FLOAT* entropyMPC_lubbysub06 = entropyMPC_lbys + 156;
entropyMPC_FLOAT entropyMPC_riub06[12];
entropyMPC_FLOAT* entropyMPC_dlubaff06 = entropyMPC_dl_aff + 156;
entropyMPC_FLOAT* entropyMPC_dsubaff06 = entropyMPC_ds_aff + 156;
entropyMPC_FLOAT* entropyMPC_dlubcc06 = entropyMPC_dl_cc + 156;
entropyMPC_FLOAT* entropyMPC_dsubcc06 = entropyMPC_ds_cc + 156;
entropyMPC_FLOAT* entropyMPC_ccrhsub06 = entropyMPC_ccrhs + 156;
entropyMPC_FLOAT entropyMPC_Phi06[12];
entropyMPC_FLOAT entropyMPC_W06[12];
entropyMPC_FLOAT entropyMPC_Ysd06[100];
entropyMPC_FLOAT entropyMPC_Lsd06[100];
entropyMPC_FLOAT* entropyMPC_z07 = entropyMPC_z + 84;
entropyMPC_FLOAT* entropyMPC_dzaff07 = entropyMPC_dz_aff + 84;
entropyMPC_FLOAT* entropyMPC_dzcc07 = entropyMPC_dz_cc + 84;
entropyMPC_FLOAT* entropyMPC_rd07 = entropyMPC_rd + 84;
entropyMPC_FLOAT entropyMPC_Lbyrd07[12];
entropyMPC_FLOAT* entropyMPC_grad_cost07 = entropyMPC_grad_cost + 84;
entropyMPC_FLOAT* entropyMPC_grad_eq07 = entropyMPC_grad_eq + 84;
entropyMPC_FLOAT* entropyMPC_grad_ineq07 = entropyMPC_grad_ineq + 84;
entropyMPC_FLOAT entropyMPC_ctv07[12];
entropyMPC_FLOAT* entropyMPC_v07 = entropyMPC_v + 70;
entropyMPC_FLOAT entropyMPC_re07[10];
entropyMPC_FLOAT entropyMPC_beta07[10];
entropyMPC_FLOAT entropyMPC_betacc07[10];
entropyMPC_FLOAT* entropyMPC_dvaff07 = entropyMPC_dv_aff + 70;
entropyMPC_FLOAT* entropyMPC_dvcc07 = entropyMPC_dv_cc + 70;
entropyMPC_FLOAT entropyMPC_V07[120];
entropyMPC_FLOAT entropyMPC_Yd07[55];
entropyMPC_FLOAT entropyMPC_Ld07[55];
entropyMPC_FLOAT entropyMPC_yy07[10];
entropyMPC_FLOAT entropyMPC_bmy07[10];
entropyMPC_FLOAT entropyMPC_c07[10] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
int entropyMPC_lbIdx07[12] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
entropyMPC_FLOAT* entropyMPC_llb07 = entropyMPC_l + 168;
entropyMPC_FLOAT* entropyMPC_slb07 = entropyMPC_s + 168;
entropyMPC_FLOAT* entropyMPC_llbbyslb07 = entropyMPC_lbys + 168;
entropyMPC_FLOAT entropyMPC_rilb07[12];
entropyMPC_FLOAT* entropyMPC_dllbaff07 = entropyMPC_dl_aff + 168;
entropyMPC_FLOAT* entropyMPC_dslbaff07 = entropyMPC_ds_aff + 168;
entropyMPC_FLOAT* entropyMPC_dllbcc07 = entropyMPC_dl_cc + 168;
entropyMPC_FLOAT* entropyMPC_dslbcc07 = entropyMPC_ds_cc + 168;
entropyMPC_FLOAT* entropyMPC_ccrhsl07 = entropyMPC_ccrhs + 168;
int entropyMPC_ubIdx07[12] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
entropyMPC_FLOAT* entropyMPC_lub07 = entropyMPC_l + 180;
entropyMPC_FLOAT* entropyMPC_sub07 = entropyMPC_s + 180;
entropyMPC_FLOAT* entropyMPC_lubbysub07 = entropyMPC_lbys + 180;
entropyMPC_FLOAT entropyMPC_riub07[12];
entropyMPC_FLOAT* entropyMPC_dlubaff07 = entropyMPC_dl_aff + 180;
entropyMPC_FLOAT* entropyMPC_dsubaff07 = entropyMPC_ds_aff + 180;
entropyMPC_FLOAT* entropyMPC_dlubcc07 = entropyMPC_dl_cc + 180;
entropyMPC_FLOAT* entropyMPC_dsubcc07 = entropyMPC_ds_cc + 180;
entropyMPC_FLOAT* entropyMPC_ccrhsub07 = entropyMPC_ccrhs + 180;
entropyMPC_FLOAT entropyMPC_Phi07[12];
entropyMPC_FLOAT entropyMPC_W07[12];
entropyMPC_FLOAT entropyMPC_Ysd07[100];
entropyMPC_FLOAT entropyMPC_Lsd07[100];
entropyMPC_FLOAT* entropyMPC_z08 = entropyMPC_z + 96;
entropyMPC_FLOAT* entropyMPC_dzaff08 = entropyMPC_dz_aff + 96;
entropyMPC_FLOAT* entropyMPC_dzcc08 = entropyMPC_dz_cc + 96;
entropyMPC_FLOAT* entropyMPC_rd08 = entropyMPC_rd + 96;
entropyMPC_FLOAT entropyMPC_Lbyrd08[12];
entropyMPC_FLOAT* entropyMPC_grad_cost08 = entropyMPC_grad_cost + 96;
entropyMPC_FLOAT* entropyMPC_grad_eq08 = entropyMPC_grad_eq + 96;
entropyMPC_FLOAT* entropyMPC_grad_ineq08 = entropyMPC_grad_ineq + 96;
entropyMPC_FLOAT entropyMPC_ctv08[12];
entropyMPC_FLOAT* entropyMPC_v08 = entropyMPC_v + 80;
entropyMPC_FLOAT entropyMPC_re08[10];
entropyMPC_FLOAT entropyMPC_beta08[10];
entropyMPC_FLOAT entropyMPC_betacc08[10];
entropyMPC_FLOAT* entropyMPC_dvaff08 = entropyMPC_dv_aff + 80;
entropyMPC_FLOAT* entropyMPC_dvcc08 = entropyMPC_dv_cc + 80;
entropyMPC_FLOAT entropyMPC_V08[120];
entropyMPC_FLOAT entropyMPC_Yd08[55];
entropyMPC_FLOAT entropyMPC_Ld08[55];
entropyMPC_FLOAT entropyMPC_yy08[10];
entropyMPC_FLOAT entropyMPC_bmy08[10];
entropyMPC_FLOAT entropyMPC_c08[10] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
int entropyMPC_lbIdx08[12] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
entropyMPC_FLOAT* entropyMPC_llb08 = entropyMPC_l + 192;
entropyMPC_FLOAT* entropyMPC_slb08 = entropyMPC_s + 192;
entropyMPC_FLOAT* entropyMPC_llbbyslb08 = entropyMPC_lbys + 192;
entropyMPC_FLOAT entropyMPC_rilb08[12];
entropyMPC_FLOAT* entropyMPC_dllbaff08 = entropyMPC_dl_aff + 192;
entropyMPC_FLOAT* entropyMPC_dslbaff08 = entropyMPC_ds_aff + 192;
entropyMPC_FLOAT* entropyMPC_dllbcc08 = entropyMPC_dl_cc + 192;
entropyMPC_FLOAT* entropyMPC_dslbcc08 = entropyMPC_ds_cc + 192;
entropyMPC_FLOAT* entropyMPC_ccrhsl08 = entropyMPC_ccrhs + 192;
int entropyMPC_ubIdx08[12] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
entropyMPC_FLOAT* entropyMPC_lub08 = entropyMPC_l + 204;
entropyMPC_FLOAT* entropyMPC_sub08 = entropyMPC_s + 204;
entropyMPC_FLOAT* entropyMPC_lubbysub08 = entropyMPC_lbys + 204;
entropyMPC_FLOAT entropyMPC_riub08[12];
entropyMPC_FLOAT* entropyMPC_dlubaff08 = entropyMPC_dl_aff + 204;
entropyMPC_FLOAT* entropyMPC_dsubaff08 = entropyMPC_ds_aff + 204;
entropyMPC_FLOAT* entropyMPC_dlubcc08 = entropyMPC_dl_cc + 204;
entropyMPC_FLOAT* entropyMPC_dsubcc08 = entropyMPC_ds_cc + 204;
entropyMPC_FLOAT* entropyMPC_ccrhsub08 = entropyMPC_ccrhs + 204;
entropyMPC_FLOAT entropyMPC_Phi08[12];
entropyMPC_FLOAT entropyMPC_W08[12];
entropyMPC_FLOAT entropyMPC_Ysd08[100];
entropyMPC_FLOAT entropyMPC_Lsd08[100];
entropyMPC_FLOAT* entropyMPC_z09 = entropyMPC_z + 108;
entropyMPC_FLOAT* entropyMPC_dzaff09 = entropyMPC_dz_aff + 108;
entropyMPC_FLOAT* entropyMPC_dzcc09 = entropyMPC_dz_cc + 108;
entropyMPC_FLOAT* entropyMPC_rd09 = entropyMPC_rd + 108;
entropyMPC_FLOAT entropyMPC_Lbyrd09[12];
entropyMPC_FLOAT* entropyMPC_grad_cost09 = entropyMPC_grad_cost + 108;
entropyMPC_FLOAT* entropyMPC_grad_eq09 = entropyMPC_grad_eq + 108;
entropyMPC_FLOAT* entropyMPC_grad_ineq09 = entropyMPC_grad_ineq + 108;
entropyMPC_FLOAT entropyMPC_ctv09[12];
entropyMPC_FLOAT* entropyMPC_v09 = entropyMPC_v + 90;
entropyMPC_FLOAT entropyMPC_re09[10];
entropyMPC_FLOAT entropyMPC_beta09[10];
entropyMPC_FLOAT entropyMPC_betacc09[10];
entropyMPC_FLOAT* entropyMPC_dvaff09 = entropyMPC_dv_aff + 90;
entropyMPC_FLOAT* entropyMPC_dvcc09 = entropyMPC_dv_cc + 90;
entropyMPC_FLOAT entropyMPC_V09[120];
entropyMPC_FLOAT entropyMPC_Yd09[55];
entropyMPC_FLOAT entropyMPC_Ld09[55];
entropyMPC_FLOAT entropyMPC_yy09[10];
entropyMPC_FLOAT entropyMPC_bmy09[10];
entropyMPC_FLOAT entropyMPC_c09[10] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
int entropyMPC_lbIdx09[12] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
entropyMPC_FLOAT* entropyMPC_llb09 = entropyMPC_l + 216;
entropyMPC_FLOAT* entropyMPC_slb09 = entropyMPC_s + 216;
entropyMPC_FLOAT* entropyMPC_llbbyslb09 = entropyMPC_lbys + 216;
entropyMPC_FLOAT entropyMPC_rilb09[12];
entropyMPC_FLOAT* entropyMPC_dllbaff09 = entropyMPC_dl_aff + 216;
entropyMPC_FLOAT* entropyMPC_dslbaff09 = entropyMPC_ds_aff + 216;
entropyMPC_FLOAT* entropyMPC_dllbcc09 = entropyMPC_dl_cc + 216;
entropyMPC_FLOAT* entropyMPC_dslbcc09 = entropyMPC_ds_cc + 216;
entropyMPC_FLOAT* entropyMPC_ccrhsl09 = entropyMPC_ccrhs + 216;
int entropyMPC_ubIdx09[12] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
entropyMPC_FLOAT* entropyMPC_lub09 = entropyMPC_l + 228;
entropyMPC_FLOAT* entropyMPC_sub09 = entropyMPC_s + 228;
entropyMPC_FLOAT* entropyMPC_lubbysub09 = entropyMPC_lbys + 228;
entropyMPC_FLOAT entropyMPC_riub09[12];
entropyMPC_FLOAT* entropyMPC_dlubaff09 = entropyMPC_dl_aff + 228;
entropyMPC_FLOAT* entropyMPC_dsubaff09 = entropyMPC_ds_aff + 228;
entropyMPC_FLOAT* entropyMPC_dlubcc09 = entropyMPC_dl_cc + 228;
entropyMPC_FLOAT* entropyMPC_dsubcc09 = entropyMPC_ds_cc + 228;
entropyMPC_FLOAT* entropyMPC_ccrhsub09 = entropyMPC_ccrhs + 228;
entropyMPC_FLOAT entropyMPC_Phi09[12];
entropyMPC_FLOAT entropyMPC_W09[12];
entropyMPC_FLOAT entropyMPC_Ysd09[100];
entropyMPC_FLOAT entropyMPC_Lsd09[100];
entropyMPC_FLOAT* entropyMPC_z10 = entropyMPC_z + 120;
entropyMPC_FLOAT* entropyMPC_dzaff10 = entropyMPC_dz_aff + 120;
entropyMPC_FLOAT* entropyMPC_dzcc10 = entropyMPC_dz_cc + 120;
entropyMPC_FLOAT* entropyMPC_rd10 = entropyMPC_rd + 120;
entropyMPC_FLOAT entropyMPC_Lbyrd10[12];
entropyMPC_FLOAT* entropyMPC_grad_cost10 = entropyMPC_grad_cost + 120;
entropyMPC_FLOAT* entropyMPC_grad_eq10 = entropyMPC_grad_eq + 120;
entropyMPC_FLOAT* entropyMPC_grad_ineq10 = entropyMPC_grad_ineq + 120;
entropyMPC_FLOAT entropyMPC_ctv10[12];
entropyMPC_FLOAT* entropyMPC_v10 = entropyMPC_v + 100;
entropyMPC_FLOAT entropyMPC_re10[10];
entropyMPC_FLOAT entropyMPC_beta10[10];
entropyMPC_FLOAT entropyMPC_betacc10[10];
entropyMPC_FLOAT* entropyMPC_dvaff10 = entropyMPC_dv_aff + 100;
entropyMPC_FLOAT* entropyMPC_dvcc10 = entropyMPC_dv_cc + 100;
entropyMPC_FLOAT entropyMPC_V10[120];
entropyMPC_FLOAT entropyMPC_Yd10[55];
entropyMPC_FLOAT entropyMPC_Ld10[55];
entropyMPC_FLOAT entropyMPC_yy10[10];
entropyMPC_FLOAT entropyMPC_bmy10[10];
entropyMPC_FLOAT entropyMPC_c10[10] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
int entropyMPC_lbIdx10[12] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
entropyMPC_FLOAT* entropyMPC_llb10 = entropyMPC_l + 240;
entropyMPC_FLOAT* entropyMPC_slb10 = entropyMPC_s + 240;
entropyMPC_FLOAT* entropyMPC_llbbyslb10 = entropyMPC_lbys + 240;
entropyMPC_FLOAT entropyMPC_rilb10[12];
entropyMPC_FLOAT* entropyMPC_dllbaff10 = entropyMPC_dl_aff + 240;
entropyMPC_FLOAT* entropyMPC_dslbaff10 = entropyMPC_ds_aff + 240;
entropyMPC_FLOAT* entropyMPC_dllbcc10 = entropyMPC_dl_cc + 240;
entropyMPC_FLOAT* entropyMPC_dslbcc10 = entropyMPC_ds_cc + 240;
entropyMPC_FLOAT* entropyMPC_ccrhsl10 = entropyMPC_ccrhs + 240;
int entropyMPC_ubIdx10[12] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
entropyMPC_FLOAT* entropyMPC_lub10 = entropyMPC_l + 252;
entropyMPC_FLOAT* entropyMPC_sub10 = entropyMPC_s + 252;
entropyMPC_FLOAT* entropyMPC_lubbysub10 = entropyMPC_lbys + 252;
entropyMPC_FLOAT entropyMPC_riub10[12];
entropyMPC_FLOAT* entropyMPC_dlubaff10 = entropyMPC_dl_aff + 252;
entropyMPC_FLOAT* entropyMPC_dsubaff10 = entropyMPC_ds_aff + 252;
entropyMPC_FLOAT* entropyMPC_dlubcc10 = entropyMPC_dl_cc + 252;
entropyMPC_FLOAT* entropyMPC_dsubcc10 = entropyMPC_ds_cc + 252;
entropyMPC_FLOAT* entropyMPC_ccrhsub10 = entropyMPC_ccrhs + 252;
entropyMPC_FLOAT entropyMPC_Phi10[12];
entropyMPC_FLOAT entropyMPC_W10[12];
entropyMPC_FLOAT entropyMPC_Ysd10[100];
entropyMPC_FLOAT entropyMPC_Lsd10[100];
entropyMPC_FLOAT* entropyMPC_z11 = entropyMPC_z + 132;
entropyMPC_FLOAT* entropyMPC_dzaff11 = entropyMPC_dz_aff + 132;
entropyMPC_FLOAT* entropyMPC_dzcc11 = entropyMPC_dz_cc + 132;
entropyMPC_FLOAT* entropyMPC_rd11 = entropyMPC_rd + 132;
entropyMPC_FLOAT entropyMPC_Lbyrd11[12];
entropyMPC_FLOAT* entropyMPC_grad_cost11 = entropyMPC_grad_cost + 132;
entropyMPC_FLOAT* entropyMPC_grad_eq11 = entropyMPC_grad_eq + 132;
entropyMPC_FLOAT* entropyMPC_grad_ineq11 = entropyMPC_grad_ineq + 132;
entropyMPC_FLOAT entropyMPC_ctv11[12];
entropyMPC_FLOAT* entropyMPC_v11 = entropyMPC_v + 110;
entropyMPC_FLOAT entropyMPC_re11[10];
entropyMPC_FLOAT entropyMPC_beta11[10];
entropyMPC_FLOAT entropyMPC_betacc11[10];
entropyMPC_FLOAT* entropyMPC_dvaff11 = entropyMPC_dv_aff + 110;
entropyMPC_FLOAT* entropyMPC_dvcc11 = entropyMPC_dv_cc + 110;
entropyMPC_FLOAT entropyMPC_V11[120];
entropyMPC_FLOAT entropyMPC_Yd11[55];
entropyMPC_FLOAT entropyMPC_Ld11[55];
entropyMPC_FLOAT entropyMPC_yy11[10];
entropyMPC_FLOAT entropyMPC_bmy11[10];
entropyMPC_FLOAT entropyMPC_c11[10] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
int entropyMPC_lbIdx11[12] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
entropyMPC_FLOAT* entropyMPC_llb11 = entropyMPC_l + 264;
entropyMPC_FLOAT* entropyMPC_slb11 = entropyMPC_s + 264;
entropyMPC_FLOAT* entropyMPC_llbbyslb11 = entropyMPC_lbys + 264;
entropyMPC_FLOAT entropyMPC_rilb11[12];
entropyMPC_FLOAT* entropyMPC_dllbaff11 = entropyMPC_dl_aff + 264;
entropyMPC_FLOAT* entropyMPC_dslbaff11 = entropyMPC_ds_aff + 264;
entropyMPC_FLOAT* entropyMPC_dllbcc11 = entropyMPC_dl_cc + 264;
entropyMPC_FLOAT* entropyMPC_dslbcc11 = entropyMPC_ds_cc + 264;
entropyMPC_FLOAT* entropyMPC_ccrhsl11 = entropyMPC_ccrhs + 264;
int entropyMPC_ubIdx11[12] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
entropyMPC_FLOAT* entropyMPC_lub11 = entropyMPC_l + 276;
entropyMPC_FLOAT* entropyMPC_sub11 = entropyMPC_s + 276;
entropyMPC_FLOAT* entropyMPC_lubbysub11 = entropyMPC_lbys + 276;
entropyMPC_FLOAT entropyMPC_riub11[12];
entropyMPC_FLOAT* entropyMPC_dlubaff11 = entropyMPC_dl_aff + 276;
entropyMPC_FLOAT* entropyMPC_dsubaff11 = entropyMPC_ds_aff + 276;
entropyMPC_FLOAT* entropyMPC_dlubcc11 = entropyMPC_dl_cc + 276;
entropyMPC_FLOAT* entropyMPC_dsubcc11 = entropyMPC_ds_cc + 276;
entropyMPC_FLOAT* entropyMPC_ccrhsub11 = entropyMPC_ccrhs + 276;
entropyMPC_FLOAT entropyMPC_Phi11[12];
entropyMPC_FLOAT entropyMPC_W11[12];
entropyMPC_FLOAT entropyMPC_Ysd11[100];
entropyMPC_FLOAT entropyMPC_Lsd11[100];
entropyMPC_FLOAT* entropyMPC_z12 = entropyMPC_z + 144;
entropyMPC_FLOAT* entropyMPC_dzaff12 = entropyMPC_dz_aff + 144;
entropyMPC_FLOAT* entropyMPC_dzcc12 = entropyMPC_dz_cc + 144;
entropyMPC_FLOAT* entropyMPC_rd12 = entropyMPC_rd + 144;
entropyMPC_FLOAT entropyMPC_Lbyrd12[12];
entropyMPC_FLOAT* entropyMPC_grad_cost12 = entropyMPC_grad_cost + 144;
entropyMPC_FLOAT* entropyMPC_grad_eq12 = entropyMPC_grad_eq + 144;
entropyMPC_FLOAT* entropyMPC_grad_ineq12 = entropyMPC_grad_ineq + 144;
entropyMPC_FLOAT entropyMPC_ctv12[12];
entropyMPC_FLOAT* entropyMPC_v12 = entropyMPC_v + 120;
entropyMPC_FLOAT entropyMPC_re12[10];
entropyMPC_FLOAT entropyMPC_beta12[10];
entropyMPC_FLOAT entropyMPC_betacc12[10];
entropyMPC_FLOAT* entropyMPC_dvaff12 = entropyMPC_dv_aff + 120;
entropyMPC_FLOAT* entropyMPC_dvcc12 = entropyMPC_dv_cc + 120;
entropyMPC_FLOAT entropyMPC_V12[120];
entropyMPC_FLOAT entropyMPC_Yd12[55];
entropyMPC_FLOAT entropyMPC_Ld12[55];
entropyMPC_FLOAT entropyMPC_yy12[10];
entropyMPC_FLOAT entropyMPC_bmy12[10];
entropyMPC_FLOAT entropyMPC_c12[10] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
int entropyMPC_lbIdx12[12] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
entropyMPC_FLOAT* entropyMPC_llb12 = entropyMPC_l + 288;
entropyMPC_FLOAT* entropyMPC_slb12 = entropyMPC_s + 288;
entropyMPC_FLOAT* entropyMPC_llbbyslb12 = entropyMPC_lbys + 288;
entropyMPC_FLOAT entropyMPC_rilb12[12];
entropyMPC_FLOAT* entropyMPC_dllbaff12 = entropyMPC_dl_aff + 288;
entropyMPC_FLOAT* entropyMPC_dslbaff12 = entropyMPC_ds_aff + 288;
entropyMPC_FLOAT* entropyMPC_dllbcc12 = entropyMPC_dl_cc + 288;
entropyMPC_FLOAT* entropyMPC_dslbcc12 = entropyMPC_ds_cc + 288;
entropyMPC_FLOAT* entropyMPC_ccrhsl12 = entropyMPC_ccrhs + 288;
int entropyMPC_ubIdx12[12] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
entropyMPC_FLOAT* entropyMPC_lub12 = entropyMPC_l + 300;
entropyMPC_FLOAT* entropyMPC_sub12 = entropyMPC_s + 300;
entropyMPC_FLOAT* entropyMPC_lubbysub12 = entropyMPC_lbys + 300;
entropyMPC_FLOAT entropyMPC_riub12[12];
entropyMPC_FLOAT* entropyMPC_dlubaff12 = entropyMPC_dl_aff + 300;
entropyMPC_FLOAT* entropyMPC_dsubaff12 = entropyMPC_ds_aff + 300;
entropyMPC_FLOAT* entropyMPC_dlubcc12 = entropyMPC_dl_cc + 300;
entropyMPC_FLOAT* entropyMPC_dsubcc12 = entropyMPC_ds_cc + 300;
entropyMPC_FLOAT* entropyMPC_ccrhsub12 = entropyMPC_ccrhs + 300;
entropyMPC_FLOAT entropyMPC_Phi12[12];
entropyMPC_FLOAT entropyMPC_W12[12];
entropyMPC_FLOAT entropyMPC_Ysd12[100];
entropyMPC_FLOAT entropyMPC_Lsd12[100];
entropyMPC_FLOAT* entropyMPC_z13 = entropyMPC_z + 156;
entropyMPC_FLOAT* entropyMPC_dzaff13 = entropyMPC_dz_aff + 156;
entropyMPC_FLOAT* entropyMPC_dzcc13 = entropyMPC_dz_cc + 156;
entropyMPC_FLOAT* entropyMPC_rd13 = entropyMPC_rd + 156;
entropyMPC_FLOAT entropyMPC_Lbyrd13[12];
entropyMPC_FLOAT* entropyMPC_grad_cost13 = entropyMPC_grad_cost + 156;
entropyMPC_FLOAT* entropyMPC_grad_eq13 = entropyMPC_grad_eq + 156;
entropyMPC_FLOAT* entropyMPC_grad_ineq13 = entropyMPC_grad_ineq + 156;
entropyMPC_FLOAT entropyMPC_ctv13[12];
entropyMPC_FLOAT* entropyMPC_v13 = entropyMPC_v + 130;
entropyMPC_FLOAT entropyMPC_re13[10];
entropyMPC_FLOAT entropyMPC_beta13[10];
entropyMPC_FLOAT entropyMPC_betacc13[10];
entropyMPC_FLOAT* entropyMPC_dvaff13 = entropyMPC_dv_aff + 130;
entropyMPC_FLOAT* entropyMPC_dvcc13 = entropyMPC_dv_cc + 130;
entropyMPC_FLOAT entropyMPC_V13[120];
entropyMPC_FLOAT entropyMPC_Yd13[55];
entropyMPC_FLOAT entropyMPC_Ld13[55];
entropyMPC_FLOAT entropyMPC_yy13[10];
entropyMPC_FLOAT entropyMPC_bmy13[10];
entropyMPC_FLOAT entropyMPC_c13[10] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
int entropyMPC_lbIdx13[12] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
entropyMPC_FLOAT* entropyMPC_llb13 = entropyMPC_l + 312;
entropyMPC_FLOAT* entropyMPC_slb13 = entropyMPC_s + 312;
entropyMPC_FLOAT* entropyMPC_llbbyslb13 = entropyMPC_lbys + 312;
entropyMPC_FLOAT entropyMPC_rilb13[12];
entropyMPC_FLOAT* entropyMPC_dllbaff13 = entropyMPC_dl_aff + 312;
entropyMPC_FLOAT* entropyMPC_dslbaff13 = entropyMPC_ds_aff + 312;
entropyMPC_FLOAT* entropyMPC_dllbcc13 = entropyMPC_dl_cc + 312;
entropyMPC_FLOAT* entropyMPC_dslbcc13 = entropyMPC_ds_cc + 312;
entropyMPC_FLOAT* entropyMPC_ccrhsl13 = entropyMPC_ccrhs + 312;
int entropyMPC_ubIdx13[12] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
entropyMPC_FLOAT* entropyMPC_lub13 = entropyMPC_l + 324;
entropyMPC_FLOAT* entropyMPC_sub13 = entropyMPC_s + 324;
entropyMPC_FLOAT* entropyMPC_lubbysub13 = entropyMPC_lbys + 324;
entropyMPC_FLOAT entropyMPC_riub13[12];
entropyMPC_FLOAT* entropyMPC_dlubaff13 = entropyMPC_dl_aff + 324;
entropyMPC_FLOAT* entropyMPC_dsubaff13 = entropyMPC_ds_aff + 324;
entropyMPC_FLOAT* entropyMPC_dlubcc13 = entropyMPC_dl_cc + 324;
entropyMPC_FLOAT* entropyMPC_dsubcc13 = entropyMPC_ds_cc + 324;
entropyMPC_FLOAT* entropyMPC_ccrhsub13 = entropyMPC_ccrhs + 324;
entropyMPC_FLOAT entropyMPC_Phi13[12];
entropyMPC_FLOAT entropyMPC_W13[12];
entropyMPC_FLOAT entropyMPC_Ysd13[100];
entropyMPC_FLOAT entropyMPC_Lsd13[100];
entropyMPC_FLOAT* entropyMPC_z14 = entropyMPC_z + 168;
entropyMPC_FLOAT* entropyMPC_dzaff14 = entropyMPC_dz_aff + 168;
entropyMPC_FLOAT* entropyMPC_dzcc14 = entropyMPC_dz_cc + 168;
entropyMPC_FLOAT* entropyMPC_rd14 = entropyMPC_rd + 168;
entropyMPC_FLOAT entropyMPC_Lbyrd14[10];
entropyMPC_FLOAT* entropyMPC_grad_cost14 = entropyMPC_grad_cost + 168;
entropyMPC_FLOAT* entropyMPC_grad_eq14 = entropyMPC_grad_eq + 168;
entropyMPC_FLOAT* entropyMPC_grad_ineq14 = entropyMPC_grad_ineq + 168;
entropyMPC_FLOAT entropyMPC_ctv14[10];
entropyMPC_FLOAT* entropyMPC_v14 = entropyMPC_v + 140;
entropyMPC_FLOAT entropyMPC_re14[10];
entropyMPC_FLOAT entropyMPC_beta14[10];
entropyMPC_FLOAT entropyMPC_betacc14[10];
entropyMPC_FLOAT* entropyMPC_dvaff14 = entropyMPC_dv_aff + 140;
entropyMPC_FLOAT* entropyMPC_dvcc14 = entropyMPC_dv_cc + 140;
entropyMPC_FLOAT entropyMPC_V14[100];
entropyMPC_FLOAT entropyMPC_Yd14[55];
entropyMPC_FLOAT entropyMPC_Ld14[55];
entropyMPC_FLOAT entropyMPC_yy14[10];
entropyMPC_FLOAT entropyMPC_bmy14[10];
entropyMPC_FLOAT entropyMPC_c14[10] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
int entropyMPC_lbIdx14[10] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
entropyMPC_FLOAT* entropyMPC_llb14 = entropyMPC_l + 336;
entropyMPC_FLOAT* entropyMPC_slb14 = entropyMPC_s + 336;
entropyMPC_FLOAT* entropyMPC_llbbyslb14 = entropyMPC_lbys + 336;
entropyMPC_FLOAT entropyMPC_rilb14[10];
entropyMPC_FLOAT* entropyMPC_dllbaff14 = entropyMPC_dl_aff + 336;
entropyMPC_FLOAT* entropyMPC_dslbaff14 = entropyMPC_ds_aff + 336;
entropyMPC_FLOAT* entropyMPC_dllbcc14 = entropyMPC_dl_cc + 336;
entropyMPC_FLOAT* entropyMPC_dslbcc14 = entropyMPC_ds_cc + 336;
entropyMPC_FLOAT* entropyMPC_ccrhsl14 = entropyMPC_ccrhs + 336;
int entropyMPC_ubIdx14[10] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
entropyMPC_FLOAT* entropyMPC_lub14 = entropyMPC_l + 346;
entropyMPC_FLOAT* entropyMPC_sub14 = entropyMPC_s + 346;
entropyMPC_FLOAT* entropyMPC_lubbysub14 = entropyMPC_lbys + 346;
entropyMPC_FLOAT entropyMPC_riub14[10];
entropyMPC_FLOAT* entropyMPC_dlubaff14 = entropyMPC_dl_aff + 346;
entropyMPC_FLOAT* entropyMPC_dsubaff14 = entropyMPC_ds_aff + 346;
entropyMPC_FLOAT* entropyMPC_dlubcc14 = entropyMPC_dl_cc + 346;
entropyMPC_FLOAT* entropyMPC_dsubcc14 = entropyMPC_ds_cc + 346;
entropyMPC_FLOAT* entropyMPC_ccrhsub14 = entropyMPC_ccrhs + 346;
entropyMPC_FLOAT entropyMPC_Phi14[10];
entropyMPC_FLOAT entropyMPC_D14[10] = {-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000};
entropyMPC_FLOAT entropyMPC_W14[10];
entropyMPC_FLOAT entropyMPC_Ysd14[100];
entropyMPC_FLOAT entropyMPC_Lsd14[100];
entropyMPC_FLOAT musigma;
entropyMPC_FLOAT sigma_3rdroot;
entropyMPC_FLOAT entropyMPC_Diag1_0[12];
entropyMPC_FLOAT entropyMPC_Diag2_0[12];
entropyMPC_FLOAT entropyMPC_L_0[66];




/* SOLVER CODE --------------------------------------------------------- */
int entropyMPC_solve(entropyMPC_params* params, entropyMPC_output* output, entropyMPC_info* info)
{	
int exitcode;

#if entropyMPC_SET_TIMING == 1
	entropyMPC_timer solvertimer;
	entropyMPC_tic(&solvertimer);
#endif
/* FUNCTION CALLS INTO LA LIBRARY -------------------------------------- */
info->it = 0;
entropyMPC_LA_INITIALIZEVECTOR_178(entropyMPC_z, 0);
entropyMPC_LA_INITIALIZEVECTOR_150(entropyMPC_v, 1);
entropyMPC_LA_INITIALIZEVECTOR_356(entropyMPC_l, 1);
entropyMPC_LA_INITIALIZEVECTOR_356(entropyMPC_s, 1);
info->mu = 0;
entropyMPC_LA_DOTACC_356(entropyMPC_l, entropyMPC_s, &info->mu);
info->mu /= 356;
while( 1 ){
info->pobj = 0;
entropyMPC_LA_DIAG_QUADFCN_12(params->H1, params->f1, entropyMPC_z00, entropyMPC_grad_cost00, &info->pobj);
entropyMPC_LA_DIAG_QUADFCN_12(params->H2, params->f2, entropyMPC_z01, entropyMPC_grad_cost01, &info->pobj);
entropyMPC_LA_DIAG_QUADFCN_12(params->H3, params->f3, entropyMPC_z02, entropyMPC_grad_cost02, &info->pobj);
entropyMPC_LA_DIAG_QUADFCN_12(params->H4, params->f4, entropyMPC_z03, entropyMPC_grad_cost03, &info->pobj);
entropyMPC_LA_DIAG_QUADFCN_12(params->H5, params->f5, entropyMPC_z04, entropyMPC_grad_cost04, &info->pobj);
entropyMPC_LA_DIAG_QUADFCN_12(params->H6, params->f6, entropyMPC_z05, entropyMPC_grad_cost05, &info->pobj);
entropyMPC_LA_DIAG_QUADFCN_12(params->H7, params->f7, entropyMPC_z06, entropyMPC_grad_cost06, &info->pobj);
entropyMPC_LA_DIAG_QUADFCN_12(params->H8, params->f8, entropyMPC_z07, entropyMPC_grad_cost07, &info->pobj);
entropyMPC_LA_DIAG_QUADFCN_12(params->H9, params->f9, entropyMPC_z08, entropyMPC_grad_cost08, &info->pobj);
entropyMPC_LA_DIAG_QUADFCN_12(params->H10, params->f10, entropyMPC_z09, entropyMPC_grad_cost09, &info->pobj);
entropyMPC_LA_DIAG_QUADFCN_12(params->H11, params->f11, entropyMPC_z10, entropyMPC_grad_cost10, &info->pobj);
entropyMPC_LA_DIAG_QUADFCN_12(params->H12, params->f12, entropyMPC_z11, entropyMPC_grad_cost11, &info->pobj);
entropyMPC_LA_DIAG_QUADFCN_12(params->H13, params->f13, entropyMPC_z12, entropyMPC_grad_cost12, &info->pobj);
entropyMPC_LA_DIAG_QUADFCN_12(params->H14, params->f14, entropyMPC_z13, entropyMPC_grad_cost13, &info->pobj);
entropyMPC_LA_DIAG_QUADFCN_10(params->H15, params->f15, entropyMPC_z14, entropyMPC_grad_cost14, &info->pobj);
info->res_eq = 0;
info->dgap = 0;
entropyMPC_LA_DIAGZERO_MVMSUB6_10(entropyMPC_D00, entropyMPC_z00, params->c1, entropyMPC_v00, entropyMPC_re00, &info->dgap, &info->res_eq);
entropyMPC_LA_DENSE_DIAGZERO_MVMSUB3_10_12_12(entropyMPC_C00, entropyMPC_z00, entropyMPC_D01, entropyMPC_z01, entropyMPC_c01, entropyMPC_v01, entropyMPC_re01, &info->dgap, &info->res_eq);
entropyMPC_LA_DENSE_DIAGZERO_MVMSUB3_10_12_12(entropyMPC_C00, entropyMPC_z01, entropyMPC_D01, entropyMPC_z02, entropyMPC_c02, entropyMPC_v02, entropyMPC_re02, &info->dgap, &info->res_eq);
entropyMPC_LA_DENSE_DIAGZERO_MVMSUB3_10_12_12(entropyMPC_C00, entropyMPC_z02, entropyMPC_D01, entropyMPC_z03, entropyMPC_c03, entropyMPC_v03, entropyMPC_re03, &info->dgap, &info->res_eq);
entropyMPC_LA_DENSE_DIAGZERO_MVMSUB3_10_12_12(entropyMPC_C00, entropyMPC_z03, entropyMPC_D01, entropyMPC_z04, entropyMPC_c04, entropyMPC_v04, entropyMPC_re04, &info->dgap, &info->res_eq);
entropyMPC_LA_DENSE_DIAGZERO_MVMSUB3_10_12_12(entropyMPC_C00, entropyMPC_z04, entropyMPC_D01, entropyMPC_z05, entropyMPC_c05, entropyMPC_v05, entropyMPC_re05, &info->dgap, &info->res_eq);
entropyMPC_LA_DENSE_DIAGZERO_MVMSUB3_10_12_12(entropyMPC_C00, entropyMPC_z05, entropyMPC_D01, entropyMPC_z06, entropyMPC_c06, entropyMPC_v06, entropyMPC_re06, &info->dgap, &info->res_eq);
entropyMPC_LA_DENSE_DIAGZERO_MVMSUB3_10_12_12(entropyMPC_C00, entropyMPC_z06, entropyMPC_D01, entropyMPC_z07, entropyMPC_c07, entropyMPC_v07, entropyMPC_re07, &info->dgap, &info->res_eq);
entropyMPC_LA_DENSE_DIAGZERO_MVMSUB3_10_12_12(entropyMPC_C00, entropyMPC_z07, entropyMPC_D01, entropyMPC_z08, entropyMPC_c08, entropyMPC_v08, entropyMPC_re08, &info->dgap, &info->res_eq);
entropyMPC_LA_DENSE_DIAGZERO_MVMSUB3_10_12_12(entropyMPC_C00, entropyMPC_z08, entropyMPC_D01, entropyMPC_z09, entropyMPC_c09, entropyMPC_v09, entropyMPC_re09, &info->dgap, &info->res_eq);
entropyMPC_LA_DENSE_DIAGZERO_MVMSUB3_10_12_12(entropyMPC_C00, entropyMPC_z09, entropyMPC_D01, entropyMPC_z10, entropyMPC_c10, entropyMPC_v10, entropyMPC_re10, &info->dgap, &info->res_eq);
entropyMPC_LA_DENSE_DIAGZERO_MVMSUB3_10_12_12(entropyMPC_C00, entropyMPC_z10, entropyMPC_D01, entropyMPC_z11, entropyMPC_c11, entropyMPC_v11, entropyMPC_re11, &info->dgap, &info->res_eq);
entropyMPC_LA_DENSE_DIAGZERO_MVMSUB3_10_12_12(entropyMPC_C00, entropyMPC_z11, entropyMPC_D01, entropyMPC_z12, entropyMPC_c12, entropyMPC_v12, entropyMPC_re12, &info->dgap, &info->res_eq);
entropyMPC_LA_DENSE_DIAGZERO_MVMSUB3_10_12_12(entropyMPC_C00, entropyMPC_z12, entropyMPC_D01, entropyMPC_z13, entropyMPC_c13, entropyMPC_v13, entropyMPC_re13, &info->dgap, &info->res_eq);
entropyMPC_LA_DENSE_DIAGZERO_MVMSUB3_10_12_10(entropyMPC_C00, entropyMPC_z13, entropyMPC_D14, entropyMPC_z14, entropyMPC_c14, entropyMPC_v14, entropyMPC_re14, &info->dgap, &info->res_eq);
entropyMPC_LA_DENSE_DIAGZERO_MTVM2_10_12_10(entropyMPC_C00, entropyMPC_v01, entropyMPC_D00, entropyMPC_v00, entropyMPC_grad_eq00);
entropyMPC_LA_DENSE_DIAGZERO_MTVM2_10_12_10(entropyMPC_C00, entropyMPC_v02, entropyMPC_D01, entropyMPC_v01, entropyMPC_grad_eq01);
entropyMPC_LA_DENSE_DIAGZERO_MTVM2_10_12_10(entropyMPC_C00, entropyMPC_v03, entropyMPC_D01, entropyMPC_v02, entropyMPC_grad_eq02);
entropyMPC_LA_DENSE_DIAGZERO_MTVM2_10_12_10(entropyMPC_C00, entropyMPC_v04, entropyMPC_D01, entropyMPC_v03, entropyMPC_grad_eq03);
entropyMPC_LA_DENSE_DIAGZERO_MTVM2_10_12_10(entropyMPC_C00, entropyMPC_v05, entropyMPC_D01, entropyMPC_v04, entropyMPC_grad_eq04);
entropyMPC_LA_DENSE_DIAGZERO_MTVM2_10_12_10(entropyMPC_C00, entropyMPC_v06, entropyMPC_D01, entropyMPC_v05, entropyMPC_grad_eq05);
entropyMPC_LA_DENSE_DIAGZERO_MTVM2_10_12_10(entropyMPC_C00, entropyMPC_v07, entropyMPC_D01, entropyMPC_v06, entropyMPC_grad_eq06);
entropyMPC_LA_DENSE_DIAGZERO_MTVM2_10_12_10(entropyMPC_C00, entropyMPC_v08, entropyMPC_D01, entropyMPC_v07, entropyMPC_grad_eq07);
entropyMPC_LA_DENSE_DIAGZERO_MTVM2_10_12_10(entropyMPC_C00, entropyMPC_v09, entropyMPC_D01, entropyMPC_v08, entropyMPC_grad_eq08);
entropyMPC_LA_DENSE_DIAGZERO_MTVM2_10_12_10(entropyMPC_C00, entropyMPC_v10, entropyMPC_D01, entropyMPC_v09, entropyMPC_grad_eq09);
entropyMPC_LA_DENSE_DIAGZERO_MTVM2_10_12_10(entropyMPC_C00, entropyMPC_v11, entropyMPC_D01, entropyMPC_v10, entropyMPC_grad_eq10);
entropyMPC_LA_DENSE_DIAGZERO_MTVM2_10_12_10(entropyMPC_C00, entropyMPC_v12, entropyMPC_D01, entropyMPC_v11, entropyMPC_grad_eq11);
entropyMPC_LA_DENSE_DIAGZERO_MTVM2_10_12_10(entropyMPC_C00, entropyMPC_v13, entropyMPC_D01, entropyMPC_v12, entropyMPC_grad_eq12);
entropyMPC_LA_DENSE_DIAGZERO_MTVM2_10_12_10(entropyMPC_C00, entropyMPC_v14, entropyMPC_D01, entropyMPC_v13, entropyMPC_grad_eq13);
entropyMPC_LA_DIAGZERO_MTVM_10_10(entropyMPC_D14, entropyMPC_v14, entropyMPC_grad_eq14);
info->res_ineq = 0;
entropyMPC_LA_VSUBADD3_12(params->lb1, entropyMPC_z00, entropyMPC_lbIdx00, entropyMPC_llb00, entropyMPC_slb00, entropyMPC_rilb00, &info->dgap, &info->res_ineq);
entropyMPC_LA_VSUBADD2_12(entropyMPC_z00, entropyMPC_ubIdx00, params->ub1, entropyMPC_lub00, entropyMPC_sub00, entropyMPC_riub00, &info->dgap, &info->res_ineq);
entropyMPC_LA_VSUBADD3_12(params->lb2, entropyMPC_z01, entropyMPC_lbIdx01, entropyMPC_llb01, entropyMPC_slb01, entropyMPC_rilb01, &info->dgap, &info->res_ineq);
entropyMPC_LA_VSUBADD2_12(entropyMPC_z01, entropyMPC_ubIdx01, params->ub2, entropyMPC_lub01, entropyMPC_sub01, entropyMPC_riub01, &info->dgap, &info->res_ineq);
entropyMPC_LA_VSUBADD3_12(params->lb3, entropyMPC_z02, entropyMPC_lbIdx02, entropyMPC_llb02, entropyMPC_slb02, entropyMPC_rilb02, &info->dgap, &info->res_ineq);
entropyMPC_LA_VSUBADD2_12(entropyMPC_z02, entropyMPC_ubIdx02, params->ub3, entropyMPC_lub02, entropyMPC_sub02, entropyMPC_riub02, &info->dgap, &info->res_ineq);
entropyMPC_LA_VSUBADD3_12(params->lb4, entropyMPC_z03, entropyMPC_lbIdx03, entropyMPC_llb03, entropyMPC_slb03, entropyMPC_rilb03, &info->dgap, &info->res_ineq);
entropyMPC_LA_VSUBADD2_12(entropyMPC_z03, entropyMPC_ubIdx03, params->ub4, entropyMPC_lub03, entropyMPC_sub03, entropyMPC_riub03, &info->dgap, &info->res_ineq);
entropyMPC_LA_VSUBADD3_12(params->lb5, entropyMPC_z04, entropyMPC_lbIdx04, entropyMPC_llb04, entropyMPC_slb04, entropyMPC_rilb04, &info->dgap, &info->res_ineq);
entropyMPC_LA_VSUBADD2_12(entropyMPC_z04, entropyMPC_ubIdx04, params->ub5, entropyMPC_lub04, entropyMPC_sub04, entropyMPC_riub04, &info->dgap, &info->res_ineq);
entropyMPC_LA_VSUBADD3_12(params->lb6, entropyMPC_z05, entropyMPC_lbIdx05, entropyMPC_llb05, entropyMPC_slb05, entropyMPC_rilb05, &info->dgap, &info->res_ineq);
entropyMPC_LA_VSUBADD2_12(entropyMPC_z05, entropyMPC_ubIdx05, params->ub6, entropyMPC_lub05, entropyMPC_sub05, entropyMPC_riub05, &info->dgap, &info->res_ineq);
entropyMPC_LA_VSUBADD3_12(params->lb7, entropyMPC_z06, entropyMPC_lbIdx06, entropyMPC_llb06, entropyMPC_slb06, entropyMPC_rilb06, &info->dgap, &info->res_ineq);
entropyMPC_LA_VSUBADD2_12(entropyMPC_z06, entropyMPC_ubIdx06, params->ub7, entropyMPC_lub06, entropyMPC_sub06, entropyMPC_riub06, &info->dgap, &info->res_ineq);
entropyMPC_LA_VSUBADD3_12(params->lb8, entropyMPC_z07, entropyMPC_lbIdx07, entropyMPC_llb07, entropyMPC_slb07, entropyMPC_rilb07, &info->dgap, &info->res_ineq);
entropyMPC_LA_VSUBADD2_12(entropyMPC_z07, entropyMPC_ubIdx07, params->ub8, entropyMPC_lub07, entropyMPC_sub07, entropyMPC_riub07, &info->dgap, &info->res_ineq);
entropyMPC_LA_VSUBADD3_12(params->lb9, entropyMPC_z08, entropyMPC_lbIdx08, entropyMPC_llb08, entropyMPC_slb08, entropyMPC_rilb08, &info->dgap, &info->res_ineq);
entropyMPC_LA_VSUBADD2_12(entropyMPC_z08, entropyMPC_ubIdx08, params->ub9, entropyMPC_lub08, entropyMPC_sub08, entropyMPC_riub08, &info->dgap, &info->res_ineq);
entropyMPC_LA_VSUBADD3_12(params->lb10, entropyMPC_z09, entropyMPC_lbIdx09, entropyMPC_llb09, entropyMPC_slb09, entropyMPC_rilb09, &info->dgap, &info->res_ineq);
entropyMPC_LA_VSUBADD2_12(entropyMPC_z09, entropyMPC_ubIdx09, params->ub10, entropyMPC_lub09, entropyMPC_sub09, entropyMPC_riub09, &info->dgap, &info->res_ineq);
entropyMPC_LA_VSUBADD3_12(params->lb11, entropyMPC_z10, entropyMPC_lbIdx10, entropyMPC_llb10, entropyMPC_slb10, entropyMPC_rilb10, &info->dgap, &info->res_ineq);
entropyMPC_LA_VSUBADD2_12(entropyMPC_z10, entropyMPC_ubIdx10, params->ub11, entropyMPC_lub10, entropyMPC_sub10, entropyMPC_riub10, &info->dgap, &info->res_ineq);
entropyMPC_LA_VSUBADD3_12(params->lb12, entropyMPC_z11, entropyMPC_lbIdx11, entropyMPC_llb11, entropyMPC_slb11, entropyMPC_rilb11, &info->dgap, &info->res_ineq);
entropyMPC_LA_VSUBADD2_12(entropyMPC_z11, entropyMPC_ubIdx11, params->ub12, entropyMPC_lub11, entropyMPC_sub11, entropyMPC_riub11, &info->dgap, &info->res_ineq);
entropyMPC_LA_VSUBADD3_12(params->lb13, entropyMPC_z12, entropyMPC_lbIdx12, entropyMPC_llb12, entropyMPC_slb12, entropyMPC_rilb12, &info->dgap, &info->res_ineq);
entropyMPC_LA_VSUBADD2_12(entropyMPC_z12, entropyMPC_ubIdx12, params->ub13, entropyMPC_lub12, entropyMPC_sub12, entropyMPC_riub12, &info->dgap, &info->res_ineq);
entropyMPC_LA_VSUBADD3_12(params->lb14, entropyMPC_z13, entropyMPC_lbIdx13, entropyMPC_llb13, entropyMPC_slb13, entropyMPC_rilb13, &info->dgap, &info->res_ineq);
entropyMPC_LA_VSUBADD2_12(entropyMPC_z13, entropyMPC_ubIdx13, params->ub14, entropyMPC_lub13, entropyMPC_sub13, entropyMPC_riub13, &info->dgap, &info->res_ineq);
entropyMPC_LA_VSUBADD3_10(params->lb15, entropyMPC_z14, entropyMPC_lbIdx14, entropyMPC_llb14, entropyMPC_slb14, entropyMPC_rilb14, &info->dgap, &info->res_ineq);
entropyMPC_LA_VSUBADD2_10(entropyMPC_z14, entropyMPC_ubIdx14, params->ub15, entropyMPC_lub14, entropyMPC_sub14, entropyMPC_riub14, &info->dgap, &info->res_ineq);
entropyMPC_LA_INEQ_B_GRAD_12_12_12(entropyMPC_lub00, entropyMPC_sub00, entropyMPC_riub00, entropyMPC_llb00, entropyMPC_slb00, entropyMPC_rilb00, entropyMPC_lbIdx00, entropyMPC_ubIdx00, entropyMPC_grad_ineq00, entropyMPC_lubbysub00, entropyMPC_llbbyslb00);
entropyMPC_LA_INEQ_B_GRAD_12_12_12(entropyMPC_lub01, entropyMPC_sub01, entropyMPC_riub01, entropyMPC_llb01, entropyMPC_slb01, entropyMPC_rilb01, entropyMPC_lbIdx01, entropyMPC_ubIdx01, entropyMPC_grad_ineq01, entropyMPC_lubbysub01, entropyMPC_llbbyslb01);
entropyMPC_LA_INEQ_B_GRAD_12_12_12(entropyMPC_lub02, entropyMPC_sub02, entropyMPC_riub02, entropyMPC_llb02, entropyMPC_slb02, entropyMPC_rilb02, entropyMPC_lbIdx02, entropyMPC_ubIdx02, entropyMPC_grad_ineq02, entropyMPC_lubbysub02, entropyMPC_llbbyslb02);
entropyMPC_LA_INEQ_B_GRAD_12_12_12(entropyMPC_lub03, entropyMPC_sub03, entropyMPC_riub03, entropyMPC_llb03, entropyMPC_slb03, entropyMPC_rilb03, entropyMPC_lbIdx03, entropyMPC_ubIdx03, entropyMPC_grad_ineq03, entropyMPC_lubbysub03, entropyMPC_llbbyslb03);
entropyMPC_LA_INEQ_B_GRAD_12_12_12(entropyMPC_lub04, entropyMPC_sub04, entropyMPC_riub04, entropyMPC_llb04, entropyMPC_slb04, entropyMPC_rilb04, entropyMPC_lbIdx04, entropyMPC_ubIdx04, entropyMPC_grad_ineq04, entropyMPC_lubbysub04, entropyMPC_llbbyslb04);
entropyMPC_LA_INEQ_B_GRAD_12_12_12(entropyMPC_lub05, entropyMPC_sub05, entropyMPC_riub05, entropyMPC_llb05, entropyMPC_slb05, entropyMPC_rilb05, entropyMPC_lbIdx05, entropyMPC_ubIdx05, entropyMPC_grad_ineq05, entropyMPC_lubbysub05, entropyMPC_llbbyslb05);
entropyMPC_LA_INEQ_B_GRAD_12_12_12(entropyMPC_lub06, entropyMPC_sub06, entropyMPC_riub06, entropyMPC_llb06, entropyMPC_slb06, entropyMPC_rilb06, entropyMPC_lbIdx06, entropyMPC_ubIdx06, entropyMPC_grad_ineq06, entropyMPC_lubbysub06, entropyMPC_llbbyslb06);
entropyMPC_LA_INEQ_B_GRAD_12_12_12(entropyMPC_lub07, entropyMPC_sub07, entropyMPC_riub07, entropyMPC_llb07, entropyMPC_slb07, entropyMPC_rilb07, entropyMPC_lbIdx07, entropyMPC_ubIdx07, entropyMPC_grad_ineq07, entropyMPC_lubbysub07, entropyMPC_llbbyslb07);
entropyMPC_LA_INEQ_B_GRAD_12_12_12(entropyMPC_lub08, entropyMPC_sub08, entropyMPC_riub08, entropyMPC_llb08, entropyMPC_slb08, entropyMPC_rilb08, entropyMPC_lbIdx08, entropyMPC_ubIdx08, entropyMPC_grad_ineq08, entropyMPC_lubbysub08, entropyMPC_llbbyslb08);
entropyMPC_LA_INEQ_B_GRAD_12_12_12(entropyMPC_lub09, entropyMPC_sub09, entropyMPC_riub09, entropyMPC_llb09, entropyMPC_slb09, entropyMPC_rilb09, entropyMPC_lbIdx09, entropyMPC_ubIdx09, entropyMPC_grad_ineq09, entropyMPC_lubbysub09, entropyMPC_llbbyslb09);
entropyMPC_LA_INEQ_B_GRAD_12_12_12(entropyMPC_lub10, entropyMPC_sub10, entropyMPC_riub10, entropyMPC_llb10, entropyMPC_slb10, entropyMPC_rilb10, entropyMPC_lbIdx10, entropyMPC_ubIdx10, entropyMPC_grad_ineq10, entropyMPC_lubbysub10, entropyMPC_llbbyslb10);
entropyMPC_LA_INEQ_B_GRAD_12_12_12(entropyMPC_lub11, entropyMPC_sub11, entropyMPC_riub11, entropyMPC_llb11, entropyMPC_slb11, entropyMPC_rilb11, entropyMPC_lbIdx11, entropyMPC_ubIdx11, entropyMPC_grad_ineq11, entropyMPC_lubbysub11, entropyMPC_llbbyslb11);
entropyMPC_LA_INEQ_B_GRAD_12_12_12(entropyMPC_lub12, entropyMPC_sub12, entropyMPC_riub12, entropyMPC_llb12, entropyMPC_slb12, entropyMPC_rilb12, entropyMPC_lbIdx12, entropyMPC_ubIdx12, entropyMPC_grad_ineq12, entropyMPC_lubbysub12, entropyMPC_llbbyslb12);
entropyMPC_LA_INEQ_B_GRAD_12_12_12(entropyMPC_lub13, entropyMPC_sub13, entropyMPC_riub13, entropyMPC_llb13, entropyMPC_slb13, entropyMPC_rilb13, entropyMPC_lbIdx13, entropyMPC_ubIdx13, entropyMPC_grad_ineq13, entropyMPC_lubbysub13, entropyMPC_llbbyslb13);
entropyMPC_LA_INEQ_B_GRAD_10_10_10(entropyMPC_lub14, entropyMPC_sub14, entropyMPC_riub14, entropyMPC_llb14, entropyMPC_slb14, entropyMPC_rilb14, entropyMPC_lbIdx14, entropyMPC_ubIdx14, entropyMPC_grad_ineq14, entropyMPC_lubbysub14, entropyMPC_llbbyslb14);
info->dobj = info->pobj - info->dgap;
info->rdgap = info->pobj ? info->dgap / info->pobj : 1e6;
if( info->rdgap < 0 ) info->rdgap = -info->rdgap;
if( info->mu < entropyMPC_SET_ACC_KKTCOMPL
    && (info->rdgap < entropyMPC_SET_ACC_RDGAP || info->dgap < entropyMPC_SET_ACC_KKTCOMPL)
    && info->res_eq < entropyMPC_SET_ACC_RESEQ
    && info->res_ineq < entropyMPC_SET_ACC_RESINEQ ){
exitcode = entropyMPC_OPTIMAL; break; }
if( info->it == entropyMPC_SET_MAXIT ){
exitcode = entropyMPC_MAXITREACHED; break; }
entropyMPC_LA_VVADD3_178(entropyMPC_grad_cost, entropyMPC_grad_eq, entropyMPC_grad_ineq, entropyMPC_rd);
entropyMPC_LA_DIAG_CHOL_ONELOOP_LBUB_12_12_12(params->H1, entropyMPC_llbbyslb00, entropyMPC_lbIdx00, entropyMPC_lubbysub00, entropyMPC_ubIdx00, entropyMPC_Phi00);
entropyMPC_LA_DIAG_MATRIXFORWARDSUB_10_12(entropyMPC_Phi00, entropyMPC_C00, entropyMPC_V00);
entropyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_10_12(entropyMPC_Phi00, entropyMPC_D00, entropyMPC_W00);
entropyMPC_LA_DENSE_DIAGZERO_MMTM_10_12_10(entropyMPC_W00, entropyMPC_V00, entropyMPC_Ysd01);
entropyMPC_LA_DIAG_FORWARDSUB_12(entropyMPC_Phi00, entropyMPC_rd00, entropyMPC_Lbyrd00);
entropyMPC_LA_DIAG_CHOL_ONELOOP_LBUB_12_12_12(params->H2, entropyMPC_llbbyslb01, entropyMPC_lbIdx01, entropyMPC_lubbysub01, entropyMPC_ubIdx01, entropyMPC_Phi01);
entropyMPC_LA_DIAG_MATRIXFORWARDSUB_10_12(entropyMPC_Phi01, entropyMPC_C00, entropyMPC_V01);
entropyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_10_12(entropyMPC_Phi01, entropyMPC_D01, entropyMPC_W01);
entropyMPC_LA_DENSE_DIAGZERO_MMTM_10_12_10(entropyMPC_W01, entropyMPC_V01, entropyMPC_Ysd02);
entropyMPC_LA_DIAG_FORWARDSUB_12(entropyMPC_Phi01, entropyMPC_rd01, entropyMPC_Lbyrd01);
entropyMPC_LA_DIAG_CHOL_ONELOOP_LBUB_12_12_12(params->H3, entropyMPC_llbbyslb02, entropyMPC_lbIdx02, entropyMPC_lubbysub02, entropyMPC_ubIdx02, entropyMPC_Phi02);
entropyMPC_LA_DIAG_MATRIXFORWARDSUB_10_12(entropyMPC_Phi02, entropyMPC_C00, entropyMPC_V02);
entropyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_10_12(entropyMPC_Phi02, entropyMPC_D01, entropyMPC_W02);
entropyMPC_LA_DENSE_DIAGZERO_MMTM_10_12_10(entropyMPC_W02, entropyMPC_V02, entropyMPC_Ysd03);
entropyMPC_LA_DIAG_FORWARDSUB_12(entropyMPC_Phi02, entropyMPC_rd02, entropyMPC_Lbyrd02);
entropyMPC_LA_DIAG_CHOL_ONELOOP_LBUB_12_12_12(params->H4, entropyMPC_llbbyslb03, entropyMPC_lbIdx03, entropyMPC_lubbysub03, entropyMPC_ubIdx03, entropyMPC_Phi03);
entropyMPC_LA_DIAG_MATRIXFORWARDSUB_10_12(entropyMPC_Phi03, entropyMPC_C00, entropyMPC_V03);
entropyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_10_12(entropyMPC_Phi03, entropyMPC_D01, entropyMPC_W03);
entropyMPC_LA_DENSE_DIAGZERO_MMTM_10_12_10(entropyMPC_W03, entropyMPC_V03, entropyMPC_Ysd04);
entropyMPC_LA_DIAG_FORWARDSUB_12(entropyMPC_Phi03, entropyMPC_rd03, entropyMPC_Lbyrd03);
entropyMPC_LA_DIAG_CHOL_ONELOOP_LBUB_12_12_12(params->H5, entropyMPC_llbbyslb04, entropyMPC_lbIdx04, entropyMPC_lubbysub04, entropyMPC_ubIdx04, entropyMPC_Phi04);
entropyMPC_LA_DIAG_MATRIXFORWARDSUB_10_12(entropyMPC_Phi04, entropyMPC_C00, entropyMPC_V04);
entropyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_10_12(entropyMPC_Phi04, entropyMPC_D01, entropyMPC_W04);
entropyMPC_LA_DENSE_DIAGZERO_MMTM_10_12_10(entropyMPC_W04, entropyMPC_V04, entropyMPC_Ysd05);
entropyMPC_LA_DIAG_FORWARDSUB_12(entropyMPC_Phi04, entropyMPC_rd04, entropyMPC_Lbyrd04);
entropyMPC_LA_DIAG_CHOL_ONELOOP_LBUB_12_12_12(params->H6, entropyMPC_llbbyslb05, entropyMPC_lbIdx05, entropyMPC_lubbysub05, entropyMPC_ubIdx05, entropyMPC_Phi05);
entropyMPC_LA_DIAG_MATRIXFORWARDSUB_10_12(entropyMPC_Phi05, entropyMPC_C00, entropyMPC_V05);
entropyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_10_12(entropyMPC_Phi05, entropyMPC_D01, entropyMPC_W05);
entropyMPC_LA_DENSE_DIAGZERO_MMTM_10_12_10(entropyMPC_W05, entropyMPC_V05, entropyMPC_Ysd06);
entropyMPC_LA_DIAG_FORWARDSUB_12(entropyMPC_Phi05, entropyMPC_rd05, entropyMPC_Lbyrd05);
entropyMPC_LA_DIAG_CHOL_ONELOOP_LBUB_12_12_12(params->H7, entropyMPC_llbbyslb06, entropyMPC_lbIdx06, entropyMPC_lubbysub06, entropyMPC_ubIdx06, entropyMPC_Phi06);
entropyMPC_LA_DIAG_MATRIXFORWARDSUB_10_12(entropyMPC_Phi06, entropyMPC_C00, entropyMPC_V06);
entropyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_10_12(entropyMPC_Phi06, entropyMPC_D01, entropyMPC_W06);
entropyMPC_LA_DENSE_DIAGZERO_MMTM_10_12_10(entropyMPC_W06, entropyMPC_V06, entropyMPC_Ysd07);
entropyMPC_LA_DIAG_FORWARDSUB_12(entropyMPC_Phi06, entropyMPC_rd06, entropyMPC_Lbyrd06);
entropyMPC_LA_DIAG_CHOL_ONELOOP_LBUB_12_12_12(params->H8, entropyMPC_llbbyslb07, entropyMPC_lbIdx07, entropyMPC_lubbysub07, entropyMPC_ubIdx07, entropyMPC_Phi07);
entropyMPC_LA_DIAG_MATRIXFORWARDSUB_10_12(entropyMPC_Phi07, entropyMPC_C00, entropyMPC_V07);
entropyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_10_12(entropyMPC_Phi07, entropyMPC_D01, entropyMPC_W07);
entropyMPC_LA_DENSE_DIAGZERO_MMTM_10_12_10(entropyMPC_W07, entropyMPC_V07, entropyMPC_Ysd08);
entropyMPC_LA_DIAG_FORWARDSUB_12(entropyMPC_Phi07, entropyMPC_rd07, entropyMPC_Lbyrd07);
entropyMPC_LA_DIAG_CHOL_ONELOOP_LBUB_12_12_12(params->H9, entropyMPC_llbbyslb08, entropyMPC_lbIdx08, entropyMPC_lubbysub08, entropyMPC_ubIdx08, entropyMPC_Phi08);
entropyMPC_LA_DIAG_MATRIXFORWARDSUB_10_12(entropyMPC_Phi08, entropyMPC_C00, entropyMPC_V08);
entropyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_10_12(entropyMPC_Phi08, entropyMPC_D01, entropyMPC_W08);
entropyMPC_LA_DENSE_DIAGZERO_MMTM_10_12_10(entropyMPC_W08, entropyMPC_V08, entropyMPC_Ysd09);
entropyMPC_LA_DIAG_FORWARDSUB_12(entropyMPC_Phi08, entropyMPC_rd08, entropyMPC_Lbyrd08);
entropyMPC_LA_DIAG_CHOL_ONELOOP_LBUB_12_12_12(params->H10, entropyMPC_llbbyslb09, entropyMPC_lbIdx09, entropyMPC_lubbysub09, entropyMPC_ubIdx09, entropyMPC_Phi09);
entropyMPC_LA_DIAG_MATRIXFORWARDSUB_10_12(entropyMPC_Phi09, entropyMPC_C00, entropyMPC_V09);
entropyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_10_12(entropyMPC_Phi09, entropyMPC_D01, entropyMPC_W09);
entropyMPC_LA_DENSE_DIAGZERO_MMTM_10_12_10(entropyMPC_W09, entropyMPC_V09, entropyMPC_Ysd10);
entropyMPC_LA_DIAG_FORWARDSUB_12(entropyMPC_Phi09, entropyMPC_rd09, entropyMPC_Lbyrd09);
entropyMPC_LA_DIAG_CHOL_ONELOOP_LBUB_12_12_12(params->H11, entropyMPC_llbbyslb10, entropyMPC_lbIdx10, entropyMPC_lubbysub10, entropyMPC_ubIdx10, entropyMPC_Phi10);
entropyMPC_LA_DIAG_MATRIXFORWARDSUB_10_12(entropyMPC_Phi10, entropyMPC_C00, entropyMPC_V10);
entropyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_10_12(entropyMPC_Phi10, entropyMPC_D01, entropyMPC_W10);
entropyMPC_LA_DENSE_DIAGZERO_MMTM_10_12_10(entropyMPC_W10, entropyMPC_V10, entropyMPC_Ysd11);
entropyMPC_LA_DIAG_FORWARDSUB_12(entropyMPC_Phi10, entropyMPC_rd10, entropyMPC_Lbyrd10);
entropyMPC_LA_DIAG_CHOL_ONELOOP_LBUB_12_12_12(params->H12, entropyMPC_llbbyslb11, entropyMPC_lbIdx11, entropyMPC_lubbysub11, entropyMPC_ubIdx11, entropyMPC_Phi11);
entropyMPC_LA_DIAG_MATRIXFORWARDSUB_10_12(entropyMPC_Phi11, entropyMPC_C00, entropyMPC_V11);
entropyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_10_12(entropyMPC_Phi11, entropyMPC_D01, entropyMPC_W11);
entropyMPC_LA_DENSE_DIAGZERO_MMTM_10_12_10(entropyMPC_W11, entropyMPC_V11, entropyMPC_Ysd12);
entropyMPC_LA_DIAG_FORWARDSUB_12(entropyMPC_Phi11, entropyMPC_rd11, entropyMPC_Lbyrd11);
entropyMPC_LA_DIAG_CHOL_ONELOOP_LBUB_12_12_12(params->H13, entropyMPC_llbbyslb12, entropyMPC_lbIdx12, entropyMPC_lubbysub12, entropyMPC_ubIdx12, entropyMPC_Phi12);
entropyMPC_LA_DIAG_MATRIXFORWARDSUB_10_12(entropyMPC_Phi12, entropyMPC_C00, entropyMPC_V12);
entropyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_10_12(entropyMPC_Phi12, entropyMPC_D01, entropyMPC_W12);
entropyMPC_LA_DENSE_DIAGZERO_MMTM_10_12_10(entropyMPC_W12, entropyMPC_V12, entropyMPC_Ysd13);
entropyMPC_LA_DIAG_FORWARDSUB_12(entropyMPC_Phi12, entropyMPC_rd12, entropyMPC_Lbyrd12);
entropyMPC_LA_DIAG_CHOL_ONELOOP_LBUB_12_12_12(params->H14, entropyMPC_llbbyslb13, entropyMPC_lbIdx13, entropyMPC_lubbysub13, entropyMPC_ubIdx13, entropyMPC_Phi13);
entropyMPC_LA_DIAG_MATRIXFORWARDSUB_10_12(entropyMPC_Phi13, entropyMPC_C00, entropyMPC_V13);
entropyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_10_12(entropyMPC_Phi13, entropyMPC_D01, entropyMPC_W13);
entropyMPC_LA_DENSE_DIAGZERO_MMTM_10_12_10(entropyMPC_W13, entropyMPC_V13, entropyMPC_Ysd14);
entropyMPC_LA_DIAG_FORWARDSUB_12(entropyMPC_Phi13, entropyMPC_rd13, entropyMPC_Lbyrd13);
entropyMPC_LA_DIAG_CHOL_ONELOOP_LBUB_10_10_10(params->H15, entropyMPC_llbbyslb14, entropyMPC_lbIdx14, entropyMPC_lubbysub14, entropyMPC_ubIdx14, entropyMPC_Phi14);
entropyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_10_10(entropyMPC_Phi14, entropyMPC_D14, entropyMPC_W14);
entropyMPC_LA_DIAG_FORWARDSUB_10(entropyMPC_Phi14, entropyMPC_rd14, entropyMPC_Lbyrd14);
entropyMPC_LA_DIAGZERO_MMT_10(entropyMPC_W00, entropyMPC_Yd00);
entropyMPC_LA_DIAGZERO_MVMSUB7_10(entropyMPC_W00, entropyMPC_Lbyrd00, entropyMPC_re00, entropyMPC_beta00);
entropyMPC_LA_DENSE_DIAGZERO_MMT2_10_12_12(entropyMPC_V00, entropyMPC_W01, entropyMPC_Yd01);
entropyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_10_12_12(entropyMPC_V00, entropyMPC_Lbyrd00, entropyMPC_W01, entropyMPC_Lbyrd01, entropyMPC_re01, entropyMPC_beta01);
entropyMPC_LA_DENSE_DIAGZERO_MMT2_10_12_12(entropyMPC_V01, entropyMPC_W02, entropyMPC_Yd02);
entropyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_10_12_12(entropyMPC_V01, entropyMPC_Lbyrd01, entropyMPC_W02, entropyMPC_Lbyrd02, entropyMPC_re02, entropyMPC_beta02);
entropyMPC_LA_DENSE_DIAGZERO_MMT2_10_12_12(entropyMPC_V02, entropyMPC_W03, entropyMPC_Yd03);
entropyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_10_12_12(entropyMPC_V02, entropyMPC_Lbyrd02, entropyMPC_W03, entropyMPC_Lbyrd03, entropyMPC_re03, entropyMPC_beta03);
entropyMPC_LA_DENSE_DIAGZERO_MMT2_10_12_12(entropyMPC_V03, entropyMPC_W04, entropyMPC_Yd04);
entropyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_10_12_12(entropyMPC_V03, entropyMPC_Lbyrd03, entropyMPC_W04, entropyMPC_Lbyrd04, entropyMPC_re04, entropyMPC_beta04);
entropyMPC_LA_DENSE_DIAGZERO_MMT2_10_12_12(entropyMPC_V04, entropyMPC_W05, entropyMPC_Yd05);
entropyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_10_12_12(entropyMPC_V04, entropyMPC_Lbyrd04, entropyMPC_W05, entropyMPC_Lbyrd05, entropyMPC_re05, entropyMPC_beta05);
entropyMPC_LA_DENSE_DIAGZERO_MMT2_10_12_12(entropyMPC_V05, entropyMPC_W06, entropyMPC_Yd06);
entropyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_10_12_12(entropyMPC_V05, entropyMPC_Lbyrd05, entropyMPC_W06, entropyMPC_Lbyrd06, entropyMPC_re06, entropyMPC_beta06);
entropyMPC_LA_DENSE_DIAGZERO_MMT2_10_12_12(entropyMPC_V06, entropyMPC_W07, entropyMPC_Yd07);
entropyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_10_12_12(entropyMPC_V06, entropyMPC_Lbyrd06, entropyMPC_W07, entropyMPC_Lbyrd07, entropyMPC_re07, entropyMPC_beta07);
entropyMPC_LA_DENSE_DIAGZERO_MMT2_10_12_12(entropyMPC_V07, entropyMPC_W08, entropyMPC_Yd08);
entropyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_10_12_12(entropyMPC_V07, entropyMPC_Lbyrd07, entropyMPC_W08, entropyMPC_Lbyrd08, entropyMPC_re08, entropyMPC_beta08);
entropyMPC_LA_DENSE_DIAGZERO_MMT2_10_12_12(entropyMPC_V08, entropyMPC_W09, entropyMPC_Yd09);
entropyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_10_12_12(entropyMPC_V08, entropyMPC_Lbyrd08, entropyMPC_W09, entropyMPC_Lbyrd09, entropyMPC_re09, entropyMPC_beta09);
entropyMPC_LA_DENSE_DIAGZERO_MMT2_10_12_12(entropyMPC_V09, entropyMPC_W10, entropyMPC_Yd10);
entropyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_10_12_12(entropyMPC_V09, entropyMPC_Lbyrd09, entropyMPC_W10, entropyMPC_Lbyrd10, entropyMPC_re10, entropyMPC_beta10);
entropyMPC_LA_DENSE_DIAGZERO_MMT2_10_12_12(entropyMPC_V10, entropyMPC_W11, entropyMPC_Yd11);
entropyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_10_12_12(entropyMPC_V10, entropyMPC_Lbyrd10, entropyMPC_W11, entropyMPC_Lbyrd11, entropyMPC_re11, entropyMPC_beta11);
entropyMPC_LA_DENSE_DIAGZERO_MMT2_10_12_12(entropyMPC_V11, entropyMPC_W12, entropyMPC_Yd12);
entropyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_10_12_12(entropyMPC_V11, entropyMPC_Lbyrd11, entropyMPC_W12, entropyMPC_Lbyrd12, entropyMPC_re12, entropyMPC_beta12);
entropyMPC_LA_DENSE_DIAGZERO_MMT2_10_12_12(entropyMPC_V12, entropyMPC_W13, entropyMPC_Yd13);
entropyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_10_12_12(entropyMPC_V12, entropyMPC_Lbyrd12, entropyMPC_W13, entropyMPC_Lbyrd13, entropyMPC_re13, entropyMPC_beta13);
entropyMPC_LA_DENSE_DIAGZERO_MMT2_10_12_10(entropyMPC_V13, entropyMPC_W14, entropyMPC_Yd14);
entropyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_10_12_10(entropyMPC_V13, entropyMPC_Lbyrd13, entropyMPC_W14, entropyMPC_Lbyrd14, entropyMPC_re14, entropyMPC_beta14);
entropyMPC_LA_DENSE_CHOL_10(entropyMPC_Yd00, entropyMPC_Ld00);
entropyMPC_LA_DENSE_FORWARDSUB_10(entropyMPC_Ld00, entropyMPC_beta00, entropyMPC_yy00);
entropyMPC_LA_DENSE_MATRIXTFORWARDSUB_10_10(entropyMPC_Ld00, entropyMPC_Ysd01, entropyMPC_Lsd01);
entropyMPC_LA_DENSE_MMTSUB_10_10(entropyMPC_Lsd01, entropyMPC_Yd01);
entropyMPC_LA_DENSE_CHOL_10(entropyMPC_Yd01, entropyMPC_Ld01);
entropyMPC_LA_DENSE_MVMSUB1_10_10(entropyMPC_Lsd01, entropyMPC_yy00, entropyMPC_beta01, entropyMPC_bmy01);
entropyMPC_LA_DENSE_FORWARDSUB_10(entropyMPC_Ld01, entropyMPC_bmy01, entropyMPC_yy01);
entropyMPC_LA_DENSE_MATRIXTFORWARDSUB_10_10(entropyMPC_Ld01, entropyMPC_Ysd02, entropyMPC_Lsd02);
entropyMPC_LA_DENSE_MMTSUB_10_10(entropyMPC_Lsd02, entropyMPC_Yd02);
entropyMPC_LA_DENSE_CHOL_10(entropyMPC_Yd02, entropyMPC_Ld02);
entropyMPC_LA_DENSE_MVMSUB1_10_10(entropyMPC_Lsd02, entropyMPC_yy01, entropyMPC_beta02, entropyMPC_bmy02);
entropyMPC_LA_DENSE_FORWARDSUB_10(entropyMPC_Ld02, entropyMPC_bmy02, entropyMPC_yy02);
entropyMPC_LA_DENSE_MATRIXTFORWARDSUB_10_10(entropyMPC_Ld02, entropyMPC_Ysd03, entropyMPC_Lsd03);
entropyMPC_LA_DENSE_MMTSUB_10_10(entropyMPC_Lsd03, entropyMPC_Yd03);
entropyMPC_LA_DENSE_CHOL_10(entropyMPC_Yd03, entropyMPC_Ld03);
entropyMPC_LA_DENSE_MVMSUB1_10_10(entropyMPC_Lsd03, entropyMPC_yy02, entropyMPC_beta03, entropyMPC_bmy03);
entropyMPC_LA_DENSE_FORWARDSUB_10(entropyMPC_Ld03, entropyMPC_bmy03, entropyMPC_yy03);
entropyMPC_LA_DENSE_MATRIXTFORWARDSUB_10_10(entropyMPC_Ld03, entropyMPC_Ysd04, entropyMPC_Lsd04);
entropyMPC_LA_DENSE_MMTSUB_10_10(entropyMPC_Lsd04, entropyMPC_Yd04);
entropyMPC_LA_DENSE_CHOL_10(entropyMPC_Yd04, entropyMPC_Ld04);
entropyMPC_LA_DENSE_MVMSUB1_10_10(entropyMPC_Lsd04, entropyMPC_yy03, entropyMPC_beta04, entropyMPC_bmy04);
entropyMPC_LA_DENSE_FORWARDSUB_10(entropyMPC_Ld04, entropyMPC_bmy04, entropyMPC_yy04);
entropyMPC_LA_DENSE_MATRIXTFORWARDSUB_10_10(entropyMPC_Ld04, entropyMPC_Ysd05, entropyMPC_Lsd05);
entropyMPC_LA_DENSE_MMTSUB_10_10(entropyMPC_Lsd05, entropyMPC_Yd05);
entropyMPC_LA_DENSE_CHOL_10(entropyMPC_Yd05, entropyMPC_Ld05);
entropyMPC_LA_DENSE_MVMSUB1_10_10(entropyMPC_Lsd05, entropyMPC_yy04, entropyMPC_beta05, entropyMPC_bmy05);
entropyMPC_LA_DENSE_FORWARDSUB_10(entropyMPC_Ld05, entropyMPC_bmy05, entropyMPC_yy05);
entropyMPC_LA_DENSE_MATRIXTFORWARDSUB_10_10(entropyMPC_Ld05, entropyMPC_Ysd06, entropyMPC_Lsd06);
entropyMPC_LA_DENSE_MMTSUB_10_10(entropyMPC_Lsd06, entropyMPC_Yd06);
entropyMPC_LA_DENSE_CHOL_10(entropyMPC_Yd06, entropyMPC_Ld06);
entropyMPC_LA_DENSE_MVMSUB1_10_10(entropyMPC_Lsd06, entropyMPC_yy05, entropyMPC_beta06, entropyMPC_bmy06);
entropyMPC_LA_DENSE_FORWARDSUB_10(entropyMPC_Ld06, entropyMPC_bmy06, entropyMPC_yy06);
entropyMPC_LA_DENSE_MATRIXTFORWARDSUB_10_10(entropyMPC_Ld06, entropyMPC_Ysd07, entropyMPC_Lsd07);
entropyMPC_LA_DENSE_MMTSUB_10_10(entropyMPC_Lsd07, entropyMPC_Yd07);
entropyMPC_LA_DENSE_CHOL_10(entropyMPC_Yd07, entropyMPC_Ld07);
entropyMPC_LA_DENSE_MVMSUB1_10_10(entropyMPC_Lsd07, entropyMPC_yy06, entropyMPC_beta07, entropyMPC_bmy07);
entropyMPC_LA_DENSE_FORWARDSUB_10(entropyMPC_Ld07, entropyMPC_bmy07, entropyMPC_yy07);
entropyMPC_LA_DENSE_MATRIXTFORWARDSUB_10_10(entropyMPC_Ld07, entropyMPC_Ysd08, entropyMPC_Lsd08);
entropyMPC_LA_DENSE_MMTSUB_10_10(entropyMPC_Lsd08, entropyMPC_Yd08);
entropyMPC_LA_DENSE_CHOL_10(entropyMPC_Yd08, entropyMPC_Ld08);
entropyMPC_LA_DENSE_MVMSUB1_10_10(entropyMPC_Lsd08, entropyMPC_yy07, entropyMPC_beta08, entropyMPC_bmy08);
entropyMPC_LA_DENSE_FORWARDSUB_10(entropyMPC_Ld08, entropyMPC_bmy08, entropyMPC_yy08);
entropyMPC_LA_DENSE_MATRIXTFORWARDSUB_10_10(entropyMPC_Ld08, entropyMPC_Ysd09, entropyMPC_Lsd09);
entropyMPC_LA_DENSE_MMTSUB_10_10(entropyMPC_Lsd09, entropyMPC_Yd09);
entropyMPC_LA_DENSE_CHOL_10(entropyMPC_Yd09, entropyMPC_Ld09);
entropyMPC_LA_DENSE_MVMSUB1_10_10(entropyMPC_Lsd09, entropyMPC_yy08, entropyMPC_beta09, entropyMPC_bmy09);
entropyMPC_LA_DENSE_FORWARDSUB_10(entropyMPC_Ld09, entropyMPC_bmy09, entropyMPC_yy09);
entropyMPC_LA_DENSE_MATRIXTFORWARDSUB_10_10(entropyMPC_Ld09, entropyMPC_Ysd10, entropyMPC_Lsd10);
entropyMPC_LA_DENSE_MMTSUB_10_10(entropyMPC_Lsd10, entropyMPC_Yd10);
entropyMPC_LA_DENSE_CHOL_10(entropyMPC_Yd10, entropyMPC_Ld10);
entropyMPC_LA_DENSE_MVMSUB1_10_10(entropyMPC_Lsd10, entropyMPC_yy09, entropyMPC_beta10, entropyMPC_bmy10);
entropyMPC_LA_DENSE_FORWARDSUB_10(entropyMPC_Ld10, entropyMPC_bmy10, entropyMPC_yy10);
entropyMPC_LA_DENSE_MATRIXTFORWARDSUB_10_10(entropyMPC_Ld10, entropyMPC_Ysd11, entropyMPC_Lsd11);
entropyMPC_LA_DENSE_MMTSUB_10_10(entropyMPC_Lsd11, entropyMPC_Yd11);
entropyMPC_LA_DENSE_CHOL_10(entropyMPC_Yd11, entropyMPC_Ld11);
entropyMPC_LA_DENSE_MVMSUB1_10_10(entropyMPC_Lsd11, entropyMPC_yy10, entropyMPC_beta11, entropyMPC_bmy11);
entropyMPC_LA_DENSE_FORWARDSUB_10(entropyMPC_Ld11, entropyMPC_bmy11, entropyMPC_yy11);
entropyMPC_LA_DENSE_MATRIXTFORWARDSUB_10_10(entropyMPC_Ld11, entropyMPC_Ysd12, entropyMPC_Lsd12);
entropyMPC_LA_DENSE_MMTSUB_10_10(entropyMPC_Lsd12, entropyMPC_Yd12);
entropyMPC_LA_DENSE_CHOL_10(entropyMPC_Yd12, entropyMPC_Ld12);
entropyMPC_LA_DENSE_MVMSUB1_10_10(entropyMPC_Lsd12, entropyMPC_yy11, entropyMPC_beta12, entropyMPC_bmy12);
entropyMPC_LA_DENSE_FORWARDSUB_10(entropyMPC_Ld12, entropyMPC_bmy12, entropyMPC_yy12);
entropyMPC_LA_DENSE_MATRIXTFORWARDSUB_10_10(entropyMPC_Ld12, entropyMPC_Ysd13, entropyMPC_Lsd13);
entropyMPC_LA_DENSE_MMTSUB_10_10(entropyMPC_Lsd13, entropyMPC_Yd13);
entropyMPC_LA_DENSE_CHOL_10(entropyMPC_Yd13, entropyMPC_Ld13);
entropyMPC_LA_DENSE_MVMSUB1_10_10(entropyMPC_Lsd13, entropyMPC_yy12, entropyMPC_beta13, entropyMPC_bmy13);
entropyMPC_LA_DENSE_FORWARDSUB_10(entropyMPC_Ld13, entropyMPC_bmy13, entropyMPC_yy13);
entropyMPC_LA_DENSE_MATRIXTFORWARDSUB_10_10(entropyMPC_Ld13, entropyMPC_Ysd14, entropyMPC_Lsd14);
entropyMPC_LA_DENSE_MMTSUB_10_10(entropyMPC_Lsd14, entropyMPC_Yd14);
entropyMPC_LA_DENSE_CHOL_10(entropyMPC_Yd14, entropyMPC_Ld14);
entropyMPC_LA_DENSE_MVMSUB1_10_10(entropyMPC_Lsd14, entropyMPC_yy13, entropyMPC_beta14, entropyMPC_bmy14);
entropyMPC_LA_DENSE_FORWARDSUB_10(entropyMPC_Ld14, entropyMPC_bmy14, entropyMPC_yy14);
entropyMPC_LA_DENSE_BACKWARDSUB_10(entropyMPC_Ld14, entropyMPC_yy14, entropyMPC_dvaff14);
entropyMPC_LA_DENSE_MTVMSUB_10_10(entropyMPC_Lsd14, entropyMPC_dvaff14, entropyMPC_yy13, entropyMPC_bmy13);
entropyMPC_LA_DENSE_BACKWARDSUB_10(entropyMPC_Ld13, entropyMPC_bmy13, entropyMPC_dvaff13);
entropyMPC_LA_DENSE_MTVMSUB_10_10(entropyMPC_Lsd13, entropyMPC_dvaff13, entropyMPC_yy12, entropyMPC_bmy12);
entropyMPC_LA_DENSE_BACKWARDSUB_10(entropyMPC_Ld12, entropyMPC_bmy12, entropyMPC_dvaff12);
entropyMPC_LA_DENSE_MTVMSUB_10_10(entropyMPC_Lsd12, entropyMPC_dvaff12, entropyMPC_yy11, entropyMPC_bmy11);
entropyMPC_LA_DENSE_BACKWARDSUB_10(entropyMPC_Ld11, entropyMPC_bmy11, entropyMPC_dvaff11);
entropyMPC_LA_DENSE_MTVMSUB_10_10(entropyMPC_Lsd11, entropyMPC_dvaff11, entropyMPC_yy10, entropyMPC_bmy10);
entropyMPC_LA_DENSE_BACKWARDSUB_10(entropyMPC_Ld10, entropyMPC_bmy10, entropyMPC_dvaff10);
entropyMPC_LA_DENSE_MTVMSUB_10_10(entropyMPC_Lsd10, entropyMPC_dvaff10, entropyMPC_yy09, entropyMPC_bmy09);
entropyMPC_LA_DENSE_BACKWARDSUB_10(entropyMPC_Ld09, entropyMPC_bmy09, entropyMPC_dvaff09);
entropyMPC_LA_DENSE_MTVMSUB_10_10(entropyMPC_Lsd09, entropyMPC_dvaff09, entropyMPC_yy08, entropyMPC_bmy08);
entropyMPC_LA_DENSE_BACKWARDSUB_10(entropyMPC_Ld08, entropyMPC_bmy08, entropyMPC_dvaff08);
entropyMPC_LA_DENSE_MTVMSUB_10_10(entropyMPC_Lsd08, entropyMPC_dvaff08, entropyMPC_yy07, entropyMPC_bmy07);
entropyMPC_LA_DENSE_BACKWARDSUB_10(entropyMPC_Ld07, entropyMPC_bmy07, entropyMPC_dvaff07);
entropyMPC_LA_DENSE_MTVMSUB_10_10(entropyMPC_Lsd07, entropyMPC_dvaff07, entropyMPC_yy06, entropyMPC_bmy06);
entropyMPC_LA_DENSE_BACKWARDSUB_10(entropyMPC_Ld06, entropyMPC_bmy06, entropyMPC_dvaff06);
entropyMPC_LA_DENSE_MTVMSUB_10_10(entropyMPC_Lsd06, entropyMPC_dvaff06, entropyMPC_yy05, entropyMPC_bmy05);
entropyMPC_LA_DENSE_BACKWARDSUB_10(entropyMPC_Ld05, entropyMPC_bmy05, entropyMPC_dvaff05);
entropyMPC_LA_DENSE_MTVMSUB_10_10(entropyMPC_Lsd05, entropyMPC_dvaff05, entropyMPC_yy04, entropyMPC_bmy04);
entropyMPC_LA_DENSE_BACKWARDSUB_10(entropyMPC_Ld04, entropyMPC_bmy04, entropyMPC_dvaff04);
entropyMPC_LA_DENSE_MTVMSUB_10_10(entropyMPC_Lsd04, entropyMPC_dvaff04, entropyMPC_yy03, entropyMPC_bmy03);
entropyMPC_LA_DENSE_BACKWARDSUB_10(entropyMPC_Ld03, entropyMPC_bmy03, entropyMPC_dvaff03);
entropyMPC_LA_DENSE_MTVMSUB_10_10(entropyMPC_Lsd03, entropyMPC_dvaff03, entropyMPC_yy02, entropyMPC_bmy02);
entropyMPC_LA_DENSE_BACKWARDSUB_10(entropyMPC_Ld02, entropyMPC_bmy02, entropyMPC_dvaff02);
entropyMPC_LA_DENSE_MTVMSUB_10_10(entropyMPC_Lsd02, entropyMPC_dvaff02, entropyMPC_yy01, entropyMPC_bmy01);
entropyMPC_LA_DENSE_BACKWARDSUB_10(entropyMPC_Ld01, entropyMPC_bmy01, entropyMPC_dvaff01);
entropyMPC_LA_DENSE_MTVMSUB_10_10(entropyMPC_Lsd01, entropyMPC_dvaff01, entropyMPC_yy00, entropyMPC_bmy00);
entropyMPC_LA_DENSE_BACKWARDSUB_10(entropyMPC_Ld00, entropyMPC_bmy00, entropyMPC_dvaff00);
entropyMPC_LA_DENSE_DIAGZERO_MTVM2_10_12_10(entropyMPC_C00, entropyMPC_dvaff01, entropyMPC_D00, entropyMPC_dvaff00, entropyMPC_grad_eq00);
entropyMPC_LA_DENSE_DIAGZERO_MTVM2_10_12_10(entropyMPC_C00, entropyMPC_dvaff02, entropyMPC_D01, entropyMPC_dvaff01, entropyMPC_grad_eq01);
entropyMPC_LA_DENSE_DIAGZERO_MTVM2_10_12_10(entropyMPC_C00, entropyMPC_dvaff03, entropyMPC_D01, entropyMPC_dvaff02, entropyMPC_grad_eq02);
entropyMPC_LA_DENSE_DIAGZERO_MTVM2_10_12_10(entropyMPC_C00, entropyMPC_dvaff04, entropyMPC_D01, entropyMPC_dvaff03, entropyMPC_grad_eq03);
entropyMPC_LA_DENSE_DIAGZERO_MTVM2_10_12_10(entropyMPC_C00, entropyMPC_dvaff05, entropyMPC_D01, entropyMPC_dvaff04, entropyMPC_grad_eq04);
entropyMPC_LA_DENSE_DIAGZERO_MTVM2_10_12_10(entropyMPC_C00, entropyMPC_dvaff06, entropyMPC_D01, entropyMPC_dvaff05, entropyMPC_grad_eq05);
entropyMPC_LA_DENSE_DIAGZERO_MTVM2_10_12_10(entropyMPC_C00, entropyMPC_dvaff07, entropyMPC_D01, entropyMPC_dvaff06, entropyMPC_grad_eq06);
entropyMPC_LA_DENSE_DIAGZERO_MTVM2_10_12_10(entropyMPC_C00, entropyMPC_dvaff08, entropyMPC_D01, entropyMPC_dvaff07, entropyMPC_grad_eq07);
entropyMPC_LA_DENSE_DIAGZERO_MTVM2_10_12_10(entropyMPC_C00, entropyMPC_dvaff09, entropyMPC_D01, entropyMPC_dvaff08, entropyMPC_grad_eq08);
entropyMPC_LA_DENSE_DIAGZERO_MTVM2_10_12_10(entropyMPC_C00, entropyMPC_dvaff10, entropyMPC_D01, entropyMPC_dvaff09, entropyMPC_grad_eq09);
entropyMPC_LA_DENSE_DIAGZERO_MTVM2_10_12_10(entropyMPC_C00, entropyMPC_dvaff11, entropyMPC_D01, entropyMPC_dvaff10, entropyMPC_grad_eq10);
entropyMPC_LA_DENSE_DIAGZERO_MTVM2_10_12_10(entropyMPC_C00, entropyMPC_dvaff12, entropyMPC_D01, entropyMPC_dvaff11, entropyMPC_grad_eq11);
entropyMPC_LA_DENSE_DIAGZERO_MTVM2_10_12_10(entropyMPC_C00, entropyMPC_dvaff13, entropyMPC_D01, entropyMPC_dvaff12, entropyMPC_grad_eq12);
entropyMPC_LA_DENSE_DIAGZERO_MTVM2_10_12_10(entropyMPC_C00, entropyMPC_dvaff14, entropyMPC_D01, entropyMPC_dvaff13, entropyMPC_grad_eq13);
entropyMPC_LA_DIAGZERO_MTVM_10_10(entropyMPC_D14, entropyMPC_dvaff14, entropyMPC_grad_eq14);
entropyMPC_LA_VSUB2_178(entropyMPC_rd, entropyMPC_grad_eq, entropyMPC_rd);
entropyMPC_LA_DIAG_FORWARDBACKWARDSUB_12(entropyMPC_Phi00, entropyMPC_rd00, entropyMPC_dzaff00);
entropyMPC_LA_DIAG_FORWARDBACKWARDSUB_12(entropyMPC_Phi01, entropyMPC_rd01, entropyMPC_dzaff01);
entropyMPC_LA_DIAG_FORWARDBACKWARDSUB_12(entropyMPC_Phi02, entropyMPC_rd02, entropyMPC_dzaff02);
entropyMPC_LA_DIAG_FORWARDBACKWARDSUB_12(entropyMPC_Phi03, entropyMPC_rd03, entropyMPC_dzaff03);
entropyMPC_LA_DIAG_FORWARDBACKWARDSUB_12(entropyMPC_Phi04, entropyMPC_rd04, entropyMPC_dzaff04);
entropyMPC_LA_DIAG_FORWARDBACKWARDSUB_12(entropyMPC_Phi05, entropyMPC_rd05, entropyMPC_dzaff05);
entropyMPC_LA_DIAG_FORWARDBACKWARDSUB_12(entropyMPC_Phi06, entropyMPC_rd06, entropyMPC_dzaff06);
entropyMPC_LA_DIAG_FORWARDBACKWARDSUB_12(entropyMPC_Phi07, entropyMPC_rd07, entropyMPC_dzaff07);
entropyMPC_LA_DIAG_FORWARDBACKWARDSUB_12(entropyMPC_Phi08, entropyMPC_rd08, entropyMPC_dzaff08);
entropyMPC_LA_DIAG_FORWARDBACKWARDSUB_12(entropyMPC_Phi09, entropyMPC_rd09, entropyMPC_dzaff09);
entropyMPC_LA_DIAG_FORWARDBACKWARDSUB_12(entropyMPC_Phi10, entropyMPC_rd10, entropyMPC_dzaff10);
entropyMPC_LA_DIAG_FORWARDBACKWARDSUB_12(entropyMPC_Phi11, entropyMPC_rd11, entropyMPC_dzaff11);
entropyMPC_LA_DIAG_FORWARDBACKWARDSUB_12(entropyMPC_Phi12, entropyMPC_rd12, entropyMPC_dzaff12);
entropyMPC_LA_DIAG_FORWARDBACKWARDSUB_12(entropyMPC_Phi13, entropyMPC_rd13, entropyMPC_dzaff13);
entropyMPC_LA_DIAG_FORWARDBACKWARDSUB_10(entropyMPC_Phi14, entropyMPC_rd14, entropyMPC_dzaff14);
entropyMPC_LA_VSUB_INDEXED_12(entropyMPC_dzaff00, entropyMPC_lbIdx00, entropyMPC_rilb00, entropyMPC_dslbaff00);
entropyMPC_LA_VSUB3_12(entropyMPC_llbbyslb00, entropyMPC_dslbaff00, entropyMPC_llb00, entropyMPC_dllbaff00);
entropyMPC_LA_VSUB2_INDEXED_12(entropyMPC_riub00, entropyMPC_dzaff00, entropyMPC_ubIdx00, entropyMPC_dsubaff00);
entropyMPC_LA_VSUB3_12(entropyMPC_lubbysub00, entropyMPC_dsubaff00, entropyMPC_lub00, entropyMPC_dlubaff00);
entropyMPC_LA_VSUB_INDEXED_12(entropyMPC_dzaff01, entropyMPC_lbIdx01, entropyMPC_rilb01, entropyMPC_dslbaff01);
entropyMPC_LA_VSUB3_12(entropyMPC_llbbyslb01, entropyMPC_dslbaff01, entropyMPC_llb01, entropyMPC_dllbaff01);
entropyMPC_LA_VSUB2_INDEXED_12(entropyMPC_riub01, entropyMPC_dzaff01, entropyMPC_ubIdx01, entropyMPC_dsubaff01);
entropyMPC_LA_VSUB3_12(entropyMPC_lubbysub01, entropyMPC_dsubaff01, entropyMPC_lub01, entropyMPC_dlubaff01);
entropyMPC_LA_VSUB_INDEXED_12(entropyMPC_dzaff02, entropyMPC_lbIdx02, entropyMPC_rilb02, entropyMPC_dslbaff02);
entropyMPC_LA_VSUB3_12(entropyMPC_llbbyslb02, entropyMPC_dslbaff02, entropyMPC_llb02, entropyMPC_dllbaff02);
entropyMPC_LA_VSUB2_INDEXED_12(entropyMPC_riub02, entropyMPC_dzaff02, entropyMPC_ubIdx02, entropyMPC_dsubaff02);
entropyMPC_LA_VSUB3_12(entropyMPC_lubbysub02, entropyMPC_dsubaff02, entropyMPC_lub02, entropyMPC_dlubaff02);
entropyMPC_LA_VSUB_INDEXED_12(entropyMPC_dzaff03, entropyMPC_lbIdx03, entropyMPC_rilb03, entropyMPC_dslbaff03);
entropyMPC_LA_VSUB3_12(entropyMPC_llbbyslb03, entropyMPC_dslbaff03, entropyMPC_llb03, entropyMPC_dllbaff03);
entropyMPC_LA_VSUB2_INDEXED_12(entropyMPC_riub03, entropyMPC_dzaff03, entropyMPC_ubIdx03, entropyMPC_dsubaff03);
entropyMPC_LA_VSUB3_12(entropyMPC_lubbysub03, entropyMPC_dsubaff03, entropyMPC_lub03, entropyMPC_dlubaff03);
entropyMPC_LA_VSUB_INDEXED_12(entropyMPC_dzaff04, entropyMPC_lbIdx04, entropyMPC_rilb04, entropyMPC_dslbaff04);
entropyMPC_LA_VSUB3_12(entropyMPC_llbbyslb04, entropyMPC_dslbaff04, entropyMPC_llb04, entropyMPC_dllbaff04);
entropyMPC_LA_VSUB2_INDEXED_12(entropyMPC_riub04, entropyMPC_dzaff04, entropyMPC_ubIdx04, entropyMPC_dsubaff04);
entropyMPC_LA_VSUB3_12(entropyMPC_lubbysub04, entropyMPC_dsubaff04, entropyMPC_lub04, entropyMPC_dlubaff04);
entropyMPC_LA_VSUB_INDEXED_12(entropyMPC_dzaff05, entropyMPC_lbIdx05, entropyMPC_rilb05, entropyMPC_dslbaff05);
entropyMPC_LA_VSUB3_12(entropyMPC_llbbyslb05, entropyMPC_dslbaff05, entropyMPC_llb05, entropyMPC_dllbaff05);
entropyMPC_LA_VSUB2_INDEXED_12(entropyMPC_riub05, entropyMPC_dzaff05, entropyMPC_ubIdx05, entropyMPC_dsubaff05);
entropyMPC_LA_VSUB3_12(entropyMPC_lubbysub05, entropyMPC_dsubaff05, entropyMPC_lub05, entropyMPC_dlubaff05);
entropyMPC_LA_VSUB_INDEXED_12(entropyMPC_dzaff06, entropyMPC_lbIdx06, entropyMPC_rilb06, entropyMPC_dslbaff06);
entropyMPC_LA_VSUB3_12(entropyMPC_llbbyslb06, entropyMPC_dslbaff06, entropyMPC_llb06, entropyMPC_dllbaff06);
entropyMPC_LA_VSUB2_INDEXED_12(entropyMPC_riub06, entropyMPC_dzaff06, entropyMPC_ubIdx06, entropyMPC_dsubaff06);
entropyMPC_LA_VSUB3_12(entropyMPC_lubbysub06, entropyMPC_dsubaff06, entropyMPC_lub06, entropyMPC_dlubaff06);
entropyMPC_LA_VSUB_INDEXED_12(entropyMPC_dzaff07, entropyMPC_lbIdx07, entropyMPC_rilb07, entropyMPC_dslbaff07);
entropyMPC_LA_VSUB3_12(entropyMPC_llbbyslb07, entropyMPC_dslbaff07, entropyMPC_llb07, entropyMPC_dllbaff07);
entropyMPC_LA_VSUB2_INDEXED_12(entropyMPC_riub07, entropyMPC_dzaff07, entropyMPC_ubIdx07, entropyMPC_dsubaff07);
entropyMPC_LA_VSUB3_12(entropyMPC_lubbysub07, entropyMPC_dsubaff07, entropyMPC_lub07, entropyMPC_dlubaff07);
entropyMPC_LA_VSUB_INDEXED_12(entropyMPC_dzaff08, entropyMPC_lbIdx08, entropyMPC_rilb08, entropyMPC_dslbaff08);
entropyMPC_LA_VSUB3_12(entropyMPC_llbbyslb08, entropyMPC_dslbaff08, entropyMPC_llb08, entropyMPC_dllbaff08);
entropyMPC_LA_VSUB2_INDEXED_12(entropyMPC_riub08, entropyMPC_dzaff08, entropyMPC_ubIdx08, entropyMPC_dsubaff08);
entropyMPC_LA_VSUB3_12(entropyMPC_lubbysub08, entropyMPC_dsubaff08, entropyMPC_lub08, entropyMPC_dlubaff08);
entropyMPC_LA_VSUB_INDEXED_12(entropyMPC_dzaff09, entropyMPC_lbIdx09, entropyMPC_rilb09, entropyMPC_dslbaff09);
entropyMPC_LA_VSUB3_12(entropyMPC_llbbyslb09, entropyMPC_dslbaff09, entropyMPC_llb09, entropyMPC_dllbaff09);
entropyMPC_LA_VSUB2_INDEXED_12(entropyMPC_riub09, entropyMPC_dzaff09, entropyMPC_ubIdx09, entropyMPC_dsubaff09);
entropyMPC_LA_VSUB3_12(entropyMPC_lubbysub09, entropyMPC_dsubaff09, entropyMPC_lub09, entropyMPC_dlubaff09);
entropyMPC_LA_VSUB_INDEXED_12(entropyMPC_dzaff10, entropyMPC_lbIdx10, entropyMPC_rilb10, entropyMPC_dslbaff10);
entropyMPC_LA_VSUB3_12(entropyMPC_llbbyslb10, entropyMPC_dslbaff10, entropyMPC_llb10, entropyMPC_dllbaff10);
entropyMPC_LA_VSUB2_INDEXED_12(entropyMPC_riub10, entropyMPC_dzaff10, entropyMPC_ubIdx10, entropyMPC_dsubaff10);
entropyMPC_LA_VSUB3_12(entropyMPC_lubbysub10, entropyMPC_dsubaff10, entropyMPC_lub10, entropyMPC_dlubaff10);
entropyMPC_LA_VSUB_INDEXED_12(entropyMPC_dzaff11, entropyMPC_lbIdx11, entropyMPC_rilb11, entropyMPC_dslbaff11);
entropyMPC_LA_VSUB3_12(entropyMPC_llbbyslb11, entropyMPC_dslbaff11, entropyMPC_llb11, entropyMPC_dllbaff11);
entropyMPC_LA_VSUB2_INDEXED_12(entropyMPC_riub11, entropyMPC_dzaff11, entropyMPC_ubIdx11, entropyMPC_dsubaff11);
entropyMPC_LA_VSUB3_12(entropyMPC_lubbysub11, entropyMPC_dsubaff11, entropyMPC_lub11, entropyMPC_dlubaff11);
entropyMPC_LA_VSUB_INDEXED_12(entropyMPC_dzaff12, entropyMPC_lbIdx12, entropyMPC_rilb12, entropyMPC_dslbaff12);
entropyMPC_LA_VSUB3_12(entropyMPC_llbbyslb12, entropyMPC_dslbaff12, entropyMPC_llb12, entropyMPC_dllbaff12);
entropyMPC_LA_VSUB2_INDEXED_12(entropyMPC_riub12, entropyMPC_dzaff12, entropyMPC_ubIdx12, entropyMPC_dsubaff12);
entropyMPC_LA_VSUB3_12(entropyMPC_lubbysub12, entropyMPC_dsubaff12, entropyMPC_lub12, entropyMPC_dlubaff12);
entropyMPC_LA_VSUB_INDEXED_12(entropyMPC_dzaff13, entropyMPC_lbIdx13, entropyMPC_rilb13, entropyMPC_dslbaff13);
entropyMPC_LA_VSUB3_12(entropyMPC_llbbyslb13, entropyMPC_dslbaff13, entropyMPC_llb13, entropyMPC_dllbaff13);
entropyMPC_LA_VSUB2_INDEXED_12(entropyMPC_riub13, entropyMPC_dzaff13, entropyMPC_ubIdx13, entropyMPC_dsubaff13);
entropyMPC_LA_VSUB3_12(entropyMPC_lubbysub13, entropyMPC_dsubaff13, entropyMPC_lub13, entropyMPC_dlubaff13);
entropyMPC_LA_VSUB_INDEXED_10(entropyMPC_dzaff14, entropyMPC_lbIdx14, entropyMPC_rilb14, entropyMPC_dslbaff14);
entropyMPC_LA_VSUB3_10(entropyMPC_llbbyslb14, entropyMPC_dslbaff14, entropyMPC_llb14, entropyMPC_dllbaff14);
entropyMPC_LA_VSUB2_INDEXED_10(entropyMPC_riub14, entropyMPC_dzaff14, entropyMPC_ubIdx14, entropyMPC_dsubaff14);
entropyMPC_LA_VSUB3_10(entropyMPC_lubbysub14, entropyMPC_dsubaff14, entropyMPC_lub14, entropyMPC_dlubaff14);
info->lsit_aff = entropyMPC_LINESEARCH_BACKTRACKING_AFFINE(entropyMPC_l, entropyMPC_s, entropyMPC_dl_aff, entropyMPC_ds_aff, &info->step_aff, &info->mu_aff);
if( info->lsit_aff == entropyMPC_NOPROGRESS ){
exitcode = entropyMPC_NOPROGRESS; break;
}
sigma_3rdroot = info->mu_aff / info->mu;
info->sigma = sigma_3rdroot*sigma_3rdroot*sigma_3rdroot;
musigma = info->mu * info->sigma;
entropyMPC_LA_VSUB5_356(entropyMPC_ds_aff, entropyMPC_dl_aff, musigma, entropyMPC_ccrhs);
entropyMPC_LA_VSUB6_INDEXED_12_12_12(entropyMPC_ccrhsub00, entropyMPC_sub00, entropyMPC_ubIdx00, entropyMPC_ccrhsl00, entropyMPC_slb00, entropyMPC_lbIdx00, entropyMPC_rd00);
entropyMPC_LA_VSUB6_INDEXED_12_12_12(entropyMPC_ccrhsub01, entropyMPC_sub01, entropyMPC_ubIdx01, entropyMPC_ccrhsl01, entropyMPC_slb01, entropyMPC_lbIdx01, entropyMPC_rd01);
entropyMPC_LA_DIAG_FORWARDSUB_12(entropyMPC_Phi00, entropyMPC_rd00, entropyMPC_Lbyrd00);
entropyMPC_LA_DIAG_FORWARDSUB_12(entropyMPC_Phi01, entropyMPC_rd01, entropyMPC_Lbyrd01);
entropyMPC_LA_DIAGZERO_MVM_10(entropyMPC_W00, entropyMPC_Lbyrd00, entropyMPC_beta00);
entropyMPC_LA_DENSE_FORWARDSUB_10(entropyMPC_Ld00, entropyMPC_beta00, entropyMPC_yy00);
entropyMPC_LA_DENSE_DIAGZERO_2MVMADD_10_12_12(entropyMPC_V00, entropyMPC_Lbyrd00, entropyMPC_W01, entropyMPC_Lbyrd01, entropyMPC_beta01);
entropyMPC_LA_DENSE_MVMSUB1_10_10(entropyMPC_Lsd01, entropyMPC_yy00, entropyMPC_beta01, entropyMPC_bmy01);
entropyMPC_LA_DENSE_FORWARDSUB_10(entropyMPC_Ld01, entropyMPC_bmy01, entropyMPC_yy01);
entropyMPC_LA_VSUB6_INDEXED_12_12_12(entropyMPC_ccrhsub02, entropyMPC_sub02, entropyMPC_ubIdx02, entropyMPC_ccrhsl02, entropyMPC_slb02, entropyMPC_lbIdx02, entropyMPC_rd02);
entropyMPC_LA_DIAG_FORWARDSUB_12(entropyMPC_Phi02, entropyMPC_rd02, entropyMPC_Lbyrd02);
entropyMPC_LA_DENSE_DIAGZERO_2MVMADD_10_12_12(entropyMPC_V01, entropyMPC_Lbyrd01, entropyMPC_W02, entropyMPC_Lbyrd02, entropyMPC_beta02);
entropyMPC_LA_DENSE_MVMSUB1_10_10(entropyMPC_Lsd02, entropyMPC_yy01, entropyMPC_beta02, entropyMPC_bmy02);
entropyMPC_LA_DENSE_FORWARDSUB_10(entropyMPC_Ld02, entropyMPC_bmy02, entropyMPC_yy02);
entropyMPC_LA_VSUB6_INDEXED_12_12_12(entropyMPC_ccrhsub03, entropyMPC_sub03, entropyMPC_ubIdx03, entropyMPC_ccrhsl03, entropyMPC_slb03, entropyMPC_lbIdx03, entropyMPC_rd03);
entropyMPC_LA_DIAG_FORWARDSUB_12(entropyMPC_Phi03, entropyMPC_rd03, entropyMPC_Lbyrd03);
entropyMPC_LA_DENSE_DIAGZERO_2MVMADD_10_12_12(entropyMPC_V02, entropyMPC_Lbyrd02, entropyMPC_W03, entropyMPC_Lbyrd03, entropyMPC_beta03);
entropyMPC_LA_DENSE_MVMSUB1_10_10(entropyMPC_Lsd03, entropyMPC_yy02, entropyMPC_beta03, entropyMPC_bmy03);
entropyMPC_LA_DENSE_FORWARDSUB_10(entropyMPC_Ld03, entropyMPC_bmy03, entropyMPC_yy03);
entropyMPC_LA_VSUB6_INDEXED_12_12_12(entropyMPC_ccrhsub04, entropyMPC_sub04, entropyMPC_ubIdx04, entropyMPC_ccrhsl04, entropyMPC_slb04, entropyMPC_lbIdx04, entropyMPC_rd04);
entropyMPC_LA_DIAG_FORWARDSUB_12(entropyMPC_Phi04, entropyMPC_rd04, entropyMPC_Lbyrd04);
entropyMPC_LA_DENSE_DIAGZERO_2MVMADD_10_12_12(entropyMPC_V03, entropyMPC_Lbyrd03, entropyMPC_W04, entropyMPC_Lbyrd04, entropyMPC_beta04);
entropyMPC_LA_DENSE_MVMSUB1_10_10(entropyMPC_Lsd04, entropyMPC_yy03, entropyMPC_beta04, entropyMPC_bmy04);
entropyMPC_LA_DENSE_FORWARDSUB_10(entropyMPC_Ld04, entropyMPC_bmy04, entropyMPC_yy04);
entropyMPC_LA_VSUB6_INDEXED_12_12_12(entropyMPC_ccrhsub05, entropyMPC_sub05, entropyMPC_ubIdx05, entropyMPC_ccrhsl05, entropyMPC_slb05, entropyMPC_lbIdx05, entropyMPC_rd05);
entropyMPC_LA_DIAG_FORWARDSUB_12(entropyMPC_Phi05, entropyMPC_rd05, entropyMPC_Lbyrd05);
entropyMPC_LA_DENSE_DIAGZERO_2MVMADD_10_12_12(entropyMPC_V04, entropyMPC_Lbyrd04, entropyMPC_W05, entropyMPC_Lbyrd05, entropyMPC_beta05);
entropyMPC_LA_DENSE_MVMSUB1_10_10(entropyMPC_Lsd05, entropyMPC_yy04, entropyMPC_beta05, entropyMPC_bmy05);
entropyMPC_LA_DENSE_FORWARDSUB_10(entropyMPC_Ld05, entropyMPC_bmy05, entropyMPC_yy05);
entropyMPC_LA_VSUB6_INDEXED_12_12_12(entropyMPC_ccrhsub06, entropyMPC_sub06, entropyMPC_ubIdx06, entropyMPC_ccrhsl06, entropyMPC_slb06, entropyMPC_lbIdx06, entropyMPC_rd06);
entropyMPC_LA_DIAG_FORWARDSUB_12(entropyMPC_Phi06, entropyMPC_rd06, entropyMPC_Lbyrd06);
entropyMPC_LA_DENSE_DIAGZERO_2MVMADD_10_12_12(entropyMPC_V05, entropyMPC_Lbyrd05, entropyMPC_W06, entropyMPC_Lbyrd06, entropyMPC_beta06);
entropyMPC_LA_DENSE_MVMSUB1_10_10(entropyMPC_Lsd06, entropyMPC_yy05, entropyMPC_beta06, entropyMPC_bmy06);
entropyMPC_LA_DENSE_FORWARDSUB_10(entropyMPC_Ld06, entropyMPC_bmy06, entropyMPC_yy06);
entropyMPC_LA_VSUB6_INDEXED_12_12_12(entropyMPC_ccrhsub07, entropyMPC_sub07, entropyMPC_ubIdx07, entropyMPC_ccrhsl07, entropyMPC_slb07, entropyMPC_lbIdx07, entropyMPC_rd07);
entropyMPC_LA_DIAG_FORWARDSUB_12(entropyMPC_Phi07, entropyMPC_rd07, entropyMPC_Lbyrd07);
entropyMPC_LA_DENSE_DIAGZERO_2MVMADD_10_12_12(entropyMPC_V06, entropyMPC_Lbyrd06, entropyMPC_W07, entropyMPC_Lbyrd07, entropyMPC_beta07);
entropyMPC_LA_DENSE_MVMSUB1_10_10(entropyMPC_Lsd07, entropyMPC_yy06, entropyMPC_beta07, entropyMPC_bmy07);
entropyMPC_LA_DENSE_FORWARDSUB_10(entropyMPC_Ld07, entropyMPC_bmy07, entropyMPC_yy07);
entropyMPC_LA_VSUB6_INDEXED_12_12_12(entropyMPC_ccrhsub08, entropyMPC_sub08, entropyMPC_ubIdx08, entropyMPC_ccrhsl08, entropyMPC_slb08, entropyMPC_lbIdx08, entropyMPC_rd08);
entropyMPC_LA_DIAG_FORWARDSUB_12(entropyMPC_Phi08, entropyMPC_rd08, entropyMPC_Lbyrd08);
entropyMPC_LA_DENSE_DIAGZERO_2MVMADD_10_12_12(entropyMPC_V07, entropyMPC_Lbyrd07, entropyMPC_W08, entropyMPC_Lbyrd08, entropyMPC_beta08);
entropyMPC_LA_DENSE_MVMSUB1_10_10(entropyMPC_Lsd08, entropyMPC_yy07, entropyMPC_beta08, entropyMPC_bmy08);
entropyMPC_LA_DENSE_FORWARDSUB_10(entropyMPC_Ld08, entropyMPC_bmy08, entropyMPC_yy08);
entropyMPC_LA_VSUB6_INDEXED_12_12_12(entropyMPC_ccrhsub09, entropyMPC_sub09, entropyMPC_ubIdx09, entropyMPC_ccrhsl09, entropyMPC_slb09, entropyMPC_lbIdx09, entropyMPC_rd09);
entropyMPC_LA_DIAG_FORWARDSUB_12(entropyMPC_Phi09, entropyMPC_rd09, entropyMPC_Lbyrd09);
entropyMPC_LA_DENSE_DIAGZERO_2MVMADD_10_12_12(entropyMPC_V08, entropyMPC_Lbyrd08, entropyMPC_W09, entropyMPC_Lbyrd09, entropyMPC_beta09);
entropyMPC_LA_DENSE_MVMSUB1_10_10(entropyMPC_Lsd09, entropyMPC_yy08, entropyMPC_beta09, entropyMPC_bmy09);
entropyMPC_LA_DENSE_FORWARDSUB_10(entropyMPC_Ld09, entropyMPC_bmy09, entropyMPC_yy09);
entropyMPC_LA_VSUB6_INDEXED_12_12_12(entropyMPC_ccrhsub10, entropyMPC_sub10, entropyMPC_ubIdx10, entropyMPC_ccrhsl10, entropyMPC_slb10, entropyMPC_lbIdx10, entropyMPC_rd10);
entropyMPC_LA_DIAG_FORWARDSUB_12(entropyMPC_Phi10, entropyMPC_rd10, entropyMPC_Lbyrd10);
entropyMPC_LA_DENSE_DIAGZERO_2MVMADD_10_12_12(entropyMPC_V09, entropyMPC_Lbyrd09, entropyMPC_W10, entropyMPC_Lbyrd10, entropyMPC_beta10);
entropyMPC_LA_DENSE_MVMSUB1_10_10(entropyMPC_Lsd10, entropyMPC_yy09, entropyMPC_beta10, entropyMPC_bmy10);
entropyMPC_LA_DENSE_FORWARDSUB_10(entropyMPC_Ld10, entropyMPC_bmy10, entropyMPC_yy10);
entropyMPC_LA_VSUB6_INDEXED_12_12_12(entropyMPC_ccrhsub11, entropyMPC_sub11, entropyMPC_ubIdx11, entropyMPC_ccrhsl11, entropyMPC_slb11, entropyMPC_lbIdx11, entropyMPC_rd11);
entropyMPC_LA_DIAG_FORWARDSUB_12(entropyMPC_Phi11, entropyMPC_rd11, entropyMPC_Lbyrd11);
entropyMPC_LA_DENSE_DIAGZERO_2MVMADD_10_12_12(entropyMPC_V10, entropyMPC_Lbyrd10, entropyMPC_W11, entropyMPC_Lbyrd11, entropyMPC_beta11);
entropyMPC_LA_DENSE_MVMSUB1_10_10(entropyMPC_Lsd11, entropyMPC_yy10, entropyMPC_beta11, entropyMPC_bmy11);
entropyMPC_LA_DENSE_FORWARDSUB_10(entropyMPC_Ld11, entropyMPC_bmy11, entropyMPC_yy11);
entropyMPC_LA_VSUB6_INDEXED_12_12_12(entropyMPC_ccrhsub12, entropyMPC_sub12, entropyMPC_ubIdx12, entropyMPC_ccrhsl12, entropyMPC_slb12, entropyMPC_lbIdx12, entropyMPC_rd12);
entropyMPC_LA_DIAG_FORWARDSUB_12(entropyMPC_Phi12, entropyMPC_rd12, entropyMPC_Lbyrd12);
entropyMPC_LA_DENSE_DIAGZERO_2MVMADD_10_12_12(entropyMPC_V11, entropyMPC_Lbyrd11, entropyMPC_W12, entropyMPC_Lbyrd12, entropyMPC_beta12);
entropyMPC_LA_DENSE_MVMSUB1_10_10(entropyMPC_Lsd12, entropyMPC_yy11, entropyMPC_beta12, entropyMPC_bmy12);
entropyMPC_LA_DENSE_FORWARDSUB_10(entropyMPC_Ld12, entropyMPC_bmy12, entropyMPC_yy12);
entropyMPC_LA_VSUB6_INDEXED_12_12_12(entropyMPC_ccrhsub13, entropyMPC_sub13, entropyMPC_ubIdx13, entropyMPC_ccrhsl13, entropyMPC_slb13, entropyMPC_lbIdx13, entropyMPC_rd13);
entropyMPC_LA_DIAG_FORWARDSUB_12(entropyMPC_Phi13, entropyMPC_rd13, entropyMPC_Lbyrd13);
entropyMPC_LA_DENSE_DIAGZERO_2MVMADD_10_12_12(entropyMPC_V12, entropyMPC_Lbyrd12, entropyMPC_W13, entropyMPC_Lbyrd13, entropyMPC_beta13);
entropyMPC_LA_DENSE_MVMSUB1_10_10(entropyMPC_Lsd13, entropyMPC_yy12, entropyMPC_beta13, entropyMPC_bmy13);
entropyMPC_LA_DENSE_FORWARDSUB_10(entropyMPC_Ld13, entropyMPC_bmy13, entropyMPC_yy13);
entropyMPC_LA_VSUB6_INDEXED_10_10_10(entropyMPC_ccrhsub14, entropyMPC_sub14, entropyMPC_ubIdx14, entropyMPC_ccrhsl14, entropyMPC_slb14, entropyMPC_lbIdx14, entropyMPC_rd14);
entropyMPC_LA_DIAG_FORWARDSUB_10(entropyMPC_Phi14, entropyMPC_rd14, entropyMPC_Lbyrd14);
entropyMPC_LA_DENSE_DIAGZERO_2MVMADD_10_12_10(entropyMPC_V13, entropyMPC_Lbyrd13, entropyMPC_W14, entropyMPC_Lbyrd14, entropyMPC_beta14);
entropyMPC_LA_DENSE_MVMSUB1_10_10(entropyMPC_Lsd14, entropyMPC_yy13, entropyMPC_beta14, entropyMPC_bmy14);
entropyMPC_LA_DENSE_FORWARDSUB_10(entropyMPC_Ld14, entropyMPC_bmy14, entropyMPC_yy14);
entropyMPC_LA_DENSE_BACKWARDSUB_10(entropyMPC_Ld14, entropyMPC_yy14, entropyMPC_dvcc14);
entropyMPC_LA_DENSE_MTVMSUB_10_10(entropyMPC_Lsd14, entropyMPC_dvcc14, entropyMPC_yy13, entropyMPC_bmy13);
entropyMPC_LA_DENSE_BACKWARDSUB_10(entropyMPC_Ld13, entropyMPC_bmy13, entropyMPC_dvcc13);
entropyMPC_LA_DENSE_MTVMSUB_10_10(entropyMPC_Lsd13, entropyMPC_dvcc13, entropyMPC_yy12, entropyMPC_bmy12);
entropyMPC_LA_DENSE_BACKWARDSUB_10(entropyMPC_Ld12, entropyMPC_bmy12, entropyMPC_dvcc12);
entropyMPC_LA_DENSE_MTVMSUB_10_10(entropyMPC_Lsd12, entropyMPC_dvcc12, entropyMPC_yy11, entropyMPC_bmy11);
entropyMPC_LA_DENSE_BACKWARDSUB_10(entropyMPC_Ld11, entropyMPC_bmy11, entropyMPC_dvcc11);
entropyMPC_LA_DENSE_MTVMSUB_10_10(entropyMPC_Lsd11, entropyMPC_dvcc11, entropyMPC_yy10, entropyMPC_bmy10);
entropyMPC_LA_DENSE_BACKWARDSUB_10(entropyMPC_Ld10, entropyMPC_bmy10, entropyMPC_dvcc10);
entropyMPC_LA_DENSE_MTVMSUB_10_10(entropyMPC_Lsd10, entropyMPC_dvcc10, entropyMPC_yy09, entropyMPC_bmy09);
entropyMPC_LA_DENSE_BACKWARDSUB_10(entropyMPC_Ld09, entropyMPC_bmy09, entropyMPC_dvcc09);
entropyMPC_LA_DENSE_MTVMSUB_10_10(entropyMPC_Lsd09, entropyMPC_dvcc09, entropyMPC_yy08, entropyMPC_bmy08);
entropyMPC_LA_DENSE_BACKWARDSUB_10(entropyMPC_Ld08, entropyMPC_bmy08, entropyMPC_dvcc08);
entropyMPC_LA_DENSE_MTVMSUB_10_10(entropyMPC_Lsd08, entropyMPC_dvcc08, entropyMPC_yy07, entropyMPC_bmy07);
entropyMPC_LA_DENSE_BACKWARDSUB_10(entropyMPC_Ld07, entropyMPC_bmy07, entropyMPC_dvcc07);
entropyMPC_LA_DENSE_MTVMSUB_10_10(entropyMPC_Lsd07, entropyMPC_dvcc07, entropyMPC_yy06, entropyMPC_bmy06);
entropyMPC_LA_DENSE_BACKWARDSUB_10(entropyMPC_Ld06, entropyMPC_bmy06, entropyMPC_dvcc06);
entropyMPC_LA_DENSE_MTVMSUB_10_10(entropyMPC_Lsd06, entropyMPC_dvcc06, entropyMPC_yy05, entropyMPC_bmy05);
entropyMPC_LA_DENSE_BACKWARDSUB_10(entropyMPC_Ld05, entropyMPC_bmy05, entropyMPC_dvcc05);
entropyMPC_LA_DENSE_MTVMSUB_10_10(entropyMPC_Lsd05, entropyMPC_dvcc05, entropyMPC_yy04, entropyMPC_bmy04);
entropyMPC_LA_DENSE_BACKWARDSUB_10(entropyMPC_Ld04, entropyMPC_bmy04, entropyMPC_dvcc04);
entropyMPC_LA_DENSE_MTVMSUB_10_10(entropyMPC_Lsd04, entropyMPC_dvcc04, entropyMPC_yy03, entropyMPC_bmy03);
entropyMPC_LA_DENSE_BACKWARDSUB_10(entropyMPC_Ld03, entropyMPC_bmy03, entropyMPC_dvcc03);
entropyMPC_LA_DENSE_MTVMSUB_10_10(entropyMPC_Lsd03, entropyMPC_dvcc03, entropyMPC_yy02, entropyMPC_bmy02);
entropyMPC_LA_DENSE_BACKWARDSUB_10(entropyMPC_Ld02, entropyMPC_bmy02, entropyMPC_dvcc02);
entropyMPC_LA_DENSE_MTVMSUB_10_10(entropyMPC_Lsd02, entropyMPC_dvcc02, entropyMPC_yy01, entropyMPC_bmy01);
entropyMPC_LA_DENSE_BACKWARDSUB_10(entropyMPC_Ld01, entropyMPC_bmy01, entropyMPC_dvcc01);
entropyMPC_LA_DENSE_MTVMSUB_10_10(entropyMPC_Lsd01, entropyMPC_dvcc01, entropyMPC_yy00, entropyMPC_bmy00);
entropyMPC_LA_DENSE_BACKWARDSUB_10(entropyMPC_Ld00, entropyMPC_bmy00, entropyMPC_dvcc00);
entropyMPC_LA_DENSE_DIAGZERO_MTVM2_10_12_10(entropyMPC_C00, entropyMPC_dvcc01, entropyMPC_D00, entropyMPC_dvcc00, entropyMPC_grad_eq00);
entropyMPC_LA_DENSE_DIAGZERO_MTVM2_10_12_10(entropyMPC_C00, entropyMPC_dvcc02, entropyMPC_D01, entropyMPC_dvcc01, entropyMPC_grad_eq01);
entropyMPC_LA_DENSE_DIAGZERO_MTVM2_10_12_10(entropyMPC_C00, entropyMPC_dvcc03, entropyMPC_D01, entropyMPC_dvcc02, entropyMPC_grad_eq02);
entropyMPC_LA_DENSE_DIAGZERO_MTVM2_10_12_10(entropyMPC_C00, entropyMPC_dvcc04, entropyMPC_D01, entropyMPC_dvcc03, entropyMPC_grad_eq03);
entropyMPC_LA_DENSE_DIAGZERO_MTVM2_10_12_10(entropyMPC_C00, entropyMPC_dvcc05, entropyMPC_D01, entropyMPC_dvcc04, entropyMPC_grad_eq04);
entropyMPC_LA_DENSE_DIAGZERO_MTVM2_10_12_10(entropyMPC_C00, entropyMPC_dvcc06, entropyMPC_D01, entropyMPC_dvcc05, entropyMPC_grad_eq05);
entropyMPC_LA_DENSE_DIAGZERO_MTVM2_10_12_10(entropyMPC_C00, entropyMPC_dvcc07, entropyMPC_D01, entropyMPC_dvcc06, entropyMPC_grad_eq06);
entropyMPC_LA_DENSE_DIAGZERO_MTVM2_10_12_10(entropyMPC_C00, entropyMPC_dvcc08, entropyMPC_D01, entropyMPC_dvcc07, entropyMPC_grad_eq07);
entropyMPC_LA_DENSE_DIAGZERO_MTVM2_10_12_10(entropyMPC_C00, entropyMPC_dvcc09, entropyMPC_D01, entropyMPC_dvcc08, entropyMPC_grad_eq08);
entropyMPC_LA_DENSE_DIAGZERO_MTVM2_10_12_10(entropyMPC_C00, entropyMPC_dvcc10, entropyMPC_D01, entropyMPC_dvcc09, entropyMPC_grad_eq09);
entropyMPC_LA_DENSE_DIAGZERO_MTVM2_10_12_10(entropyMPC_C00, entropyMPC_dvcc11, entropyMPC_D01, entropyMPC_dvcc10, entropyMPC_grad_eq10);
entropyMPC_LA_DENSE_DIAGZERO_MTVM2_10_12_10(entropyMPC_C00, entropyMPC_dvcc12, entropyMPC_D01, entropyMPC_dvcc11, entropyMPC_grad_eq11);
entropyMPC_LA_DENSE_DIAGZERO_MTVM2_10_12_10(entropyMPC_C00, entropyMPC_dvcc13, entropyMPC_D01, entropyMPC_dvcc12, entropyMPC_grad_eq12);
entropyMPC_LA_DENSE_DIAGZERO_MTVM2_10_12_10(entropyMPC_C00, entropyMPC_dvcc14, entropyMPC_D01, entropyMPC_dvcc13, entropyMPC_grad_eq13);
entropyMPC_LA_DIAGZERO_MTVM_10_10(entropyMPC_D14, entropyMPC_dvcc14, entropyMPC_grad_eq14);
entropyMPC_LA_VSUB_178(entropyMPC_rd, entropyMPC_grad_eq, entropyMPC_rd);
entropyMPC_LA_DIAG_FORWARDBACKWARDSUB_12(entropyMPC_Phi00, entropyMPC_rd00, entropyMPC_dzcc00);
entropyMPC_LA_DIAG_FORWARDBACKWARDSUB_12(entropyMPC_Phi01, entropyMPC_rd01, entropyMPC_dzcc01);
entropyMPC_LA_DIAG_FORWARDBACKWARDSUB_12(entropyMPC_Phi02, entropyMPC_rd02, entropyMPC_dzcc02);
entropyMPC_LA_DIAG_FORWARDBACKWARDSUB_12(entropyMPC_Phi03, entropyMPC_rd03, entropyMPC_dzcc03);
entropyMPC_LA_DIAG_FORWARDBACKWARDSUB_12(entropyMPC_Phi04, entropyMPC_rd04, entropyMPC_dzcc04);
entropyMPC_LA_DIAG_FORWARDBACKWARDSUB_12(entropyMPC_Phi05, entropyMPC_rd05, entropyMPC_dzcc05);
entropyMPC_LA_DIAG_FORWARDBACKWARDSUB_12(entropyMPC_Phi06, entropyMPC_rd06, entropyMPC_dzcc06);
entropyMPC_LA_DIAG_FORWARDBACKWARDSUB_12(entropyMPC_Phi07, entropyMPC_rd07, entropyMPC_dzcc07);
entropyMPC_LA_DIAG_FORWARDBACKWARDSUB_12(entropyMPC_Phi08, entropyMPC_rd08, entropyMPC_dzcc08);
entropyMPC_LA_DIAG_FORWARDBACKWARDSUB_12(entropyMPC_Phi09, entropyMPC_rd09, entropyMPC_dzcc09);
entropyMPC_LA_DIAG_FORWARDBACKWARDSUB_12(entropyMPC_Phi10, entropyMPC_rd10, entropyMPC_dzcc10);
entropyMPC_LA_DIAG_FORWARDBACKWARDSUB_12(entropyMPC_Phi11, entropyMPC_rd11, entropyMPC_dzcc11);
entropyMPC_LA_DIAG_FORWARDBACKWARDSUB_12(entropyMPC_Phi12, entropyMPC_rd12, entropyMPC_dzcc12);
entropyMPC_LA_DIAG_FORWARDBACKWARDSUB_12(entropyMPC_Phi13, entropyMPC_rd13, entropyMPC_dzcc13);
entropyMPC_LA_DIAG_FORWARDBACKWARDSUB_10(entropyMPC_Phi14, entropyMPC_rd14, entropyMPC_dzcc14);
entropyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_12(entropyMPC_ccrhsl00, entropyMPC_slb00, entropyMPC_llbbyslb00, entropyMPC_dzcc00, entropyMPC_lbIdx00, entropyMPC_dllbcc00);
entropyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_12(entropyMPC_ccrhsub00, entropyMPC_sub00, entropyMPC_lubbysub00, entropyMPC_dzcc00, entropyMPC_ubIdx00, entropyMPC_dlubcc00);
entropyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_12(entropyMPC_ccrhsl01, entropyMPC_slb01, entropyMPC_llbbyslb01, entropyMPC_dzcc01, entropyMPC_lbIdx01, entropyMPC_dllbcc01);
entropyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_12(entropyMPC_ccrhsub01, entropyMPC_sub01, entropyMPC_lubbysub01, entropyMPC_dzcc01, entropyMPC_ubIdx01, entropyMPC_dlubcc01);
entropyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_12(entropyMPC_ccrhsl02, entropyMPC_slb02, entropyMPC_llbbyslb02, entropyMPC_dzcc02, entropyMPC_lbIdx02, entropyMPC_dllbcc02);
entropyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_12(entropyMPC_ccrhsub02, entropyMPC_sub02, entropyMPC_lubbysub02, entropyMPC_dzcc02, entropyMPC_ubIdx02, entropyMPC_dlubcc02);
entropyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_12(entropyMPC_ccrhsl03, entropyMPC_slb03, entropyMPC_llbbyslb03, entropyMPC_dzcc03, entropyMPC_lbIdx03, entropyMPC_dllbcc03);
entropyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_12(entropyMPC_ccrhsub03, entropyMPC_sub03, entropyMPC_lubbysub03, entropyMPC_dzcc03, entropyMPC_ubIdx03, entropyMPC_dlubcc03);
entropyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_12(entropyMPC_ccrhsl04, entropyMPC_slb04, entropyMPC_llbbyslb04, entropyMPC_dzcc04, entropyMPC_lbIdx04, entropyMPC_dllbcc04);
entropyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_12(entropyMPC_ccrhsub04, entropyMPC_sub04, entropyMPC_lubbysub04, entropyMPC_dzcc04, entropyMPC_ubIdx04, entropyMPC_dlubcc04);
entropyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_12(entropyMPC_ccrhsl05, entropyMPC_slb05, entropyMPC_llbbyslb05, entropyMPC_dzcc05, entropyMPC_lbIdx05, entropyMPC_dllbcc05);
entropyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_12(entropyMPC_ccrhsub05, entropyMPC_sub05, entropyMPC_lubbysub05, entropyMPC_dzcc05, entropyMPC_ubIdx05, entropyMPC_dlubcc05);
entropyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_12(entropyMPC_ccrhsl06, entropyMPC_slb06, entropyMPC_llbbyslb06, entropyMPC_dzcc06, entropyMPC_lbIdx06, entropyMPC_dllbcc06);
entropyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_12(entropyMPC_ccrhsub06, entropyMPC_sub06, entropyMPC_lubbysub06, entropyMPC_dzcc06, entropyMPC_ubIdx06, entropyMPC_dlubcc06);
entropyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_12(entropyMPC_ccrhsl07, entropyMPC_slb07, entropyMPC_llbbyslb07, entropyMPC_dzcc07, entropyMPC_lbIdx07, entropyMPC_dllbcc07);
entropyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_12(entropyMPC_ccrhsub07, entropyMPC_sub07, entropyMPC_lubbysub07, entropyMPC_dzcc07, entropyMPC_ubIdx07, entropyMPC_dlubcc07);
entropyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_12(entropyMPC_ccrhsl08, entropyMPC_slb08, entropyMPC_llbbyslb08, entropyMPC_dzcc08, entropyMPC_lbIdx08, entropyMPC_dllbcc08);
entropyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_12(entropyMPC_ccrhsub08, entropyMPC_sub08, entropyMPC_lubbysub08, entropyMPC_dzcc08, entropyMPC_ubIdx08, entropyMPC_dlubcc08);
entropyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_12(entropyMPC_ccrhsl09, entropyMPC_slb09, entropyMPC_llbbyslb09, entropyMPC_dzcc09, entropyMPC_lbIdx09, entropyMPC_dllbcc09);
entropyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_12(entropyMPC_ccrhsub09, entropyMPC_sub09, entropyMPC_lubbysub09, entropyMPC_dzcc09, entropyMPC_ubIdx09, entropyMPC_dlubcc09);
entropyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_12(entropyMPC_ccrhsl10, entropyMPC_slb10, entropyMPC_llbbyslb10, entropyMPC_dzcc10, entropyMPC_lbIdx10, entropyMPC_dllbcc10);
entropyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_12(entropyMPC_ccrhsub10, entropyMPC_sub10, entropyMPC_lubbysub10, entropyMPC_dzcc10, entropyMPC_ubIdx10, entropyMPC_dlubcc10);
entropyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_12(entropyMPC_ccrhsl11, entropyMPC_slb11, entropyMPC_llbbyslb11, entropyMPC_dzcc11, entropyMPC_lbIdx11, entropyMPC_dllbcc11);
entropyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_12(entropyMPC_ccrhsub11, entropyMPC_sub11, entropyMPC_lubbysub11, entropyMPC_dzcc11, entropyMPC_ubIdx11, entropyMPC_dlubcc11);
entropyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_12(entropyMPC_ccrhsl12, entropyMPC_slb12, entropyMPC_llbbyslb12, entropyMPC_dzcc12, entropyMPC_lbIdx12, entropyMPC_dllbcc12);
entropyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_12(entropyMPC_ccrhsub12, entropyMPC_sub12, entropyMPC_lubbysub12, entropyMPC_dzcc12, entropyMPC_ubIdx12, entropyMPC_dlubcc12);
entropyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_12(entropyMPC_ccrhsl13, entropyMPC_slb13, entropyMPC_llbbyslb13, entropyMPC_dzcc13, entropyMPC_lbIdx13, entropyMPC_dllbcc13);
entropyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_12(entropyMPC_ccrhsub13, entropyMPC_sub13, entropyMPC_lubbysub13, entropyMPC_dzcc13, entropyMPC_ubIdx13, entropyMPC_dlubcc13);
entropyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_10(entropyMPC_ccrhsl14, entropyMPC_slb14, entropyMPC_llbbyslb14, entropyMPC_dzcc14, entropyMPC_lbIdx14, entropyMPC_dllbcc14);
entropyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_10(entropyMPC_ccrhsub14, entropyMPC_sub14, entropyMPC_lubbysub14, entropyMPC_dzcc14, entropyMPC_ubIdx14, entropyMPC_dlubcc14);
entropyMPC_LA_VSUB7_356(entropyMPC_l, entropyMPC_ccrhs, entropyMPC_s, entropyMPC_dl_cc, entropyMPC_ds_cc);
entropyMPC_LA_VADD_178(entropyMPC_dz_cc, entropyMPC_dz_aff);
entropyMPC_LA_VADD_150(entropyMPC_dv_cc, entropyMPC_dv_aff);
entropyMPC_LA_VADD_356(entropyMPC_dl_cc, entropyMPC_dl_aff);
entropyMPC_LA_VADD_356(entropyMPC_ds_cc, entropyMPC_ds_aff);
info->lsit_cc = entropyMPC_LINESEARCH_BACKTRACKING_COMBINED(entropyMPC_z, entropyMPC_v, entropyMPC_l, entropyMPC_s, entropyMPC_dz_cc, entropyMPC_dv_cc, entropyMPC_dl_cc, entropyMPC_ds_cc, &info->step_cc, &info->mu);
if( info->lsit_cc == entropyMPC_NOPROGRESS ){
exitcode = entropyMPC_NOPROGRESS; break;
}
info->it++;
}
output->z1[0] = entropyMPC_z00[0];
output->z1[1] = entropyMPC_z00[1];
output->z1[2] = entropyMPC_z00[2];
output->z1[3] = entropyMPC_z00[3];
output->z1[4] = entropyMPC_z00[4];
output->z1[5] = entropyMPC_z00[5];
output->z1[6] = entropyMPC_z00[6];
output->z1[7] = entropyMPC_z00[7];
output->z1[8] = entropyMPC_z00[8];
output->z1[9] = entropyMPC_z00[9];
output->z1[10] = entropyMPC_z00[10];
output->z1[11] = entropyMPC_z00[11];
output->z2[0] = entropyMPC_z01[0];
output->z2[1] = entropyMPC_z01[1];
output->z2[2] = entropyMPC_z01[2];
output->z2[3] = entropyMPC_z01[3];
output->z2[4] = entropyMPC_z01[4];
output->z2[5] = entropyMPC_z01[5];
output->z2[6] = entropyMPC_z01[6];
output->z2[7] = entropyMPC_z01[7];
output->z2[8] = entropyMPC_z01[8];
output->z2[9] = entropyMPC_z01[9];
output->z2[10] = entropyMPC_z01[10];
output->z2[11] = entropyMPC_z01[11];
output->z3[0] = entropyMPC_z02[0];
output->z3[1] = entropyMPC_z02[1];
output->z3[2] = entropyMPC_z02[2];
output->z3[3] = entropyMPC_z02[3];
output->z3[4] = entropyMPC_z02[4];
output->z3[5] = entropyMPC_z02[5];
output->z3[6] = entropyMPC_z02[6];
output->z3[7] = entropyMPC_z02[7];
output->z3[8] = entropyMPC_z02[8];
output->z3[9] = entropyMPC_z02[9];
output->z3[10] = entropyMPC_z02[10];
output->z3[11] = entropyMPC_z02[11];
output->z4[0] = entropyMPC_z03[0];
output->z4[1] = entropyMPC_z03[1];
output->z4[2] = entropyMPC_z03[2];
output->z4[3] = entropyMPC_z03[3];
output->z4[4] = entropyMPC_z03[4];
output->z4[5] = entropyMPC_z03[5];
output->z4[6] = entropyMPC_z03[6];
output->z4[7] = entropyMPC_z03[7];
output->z4[8] = entropyMPC_z03[8];
output->z4[9] = entropyMPC_z03[9];
output->z4[10] = entropyMPC_z03[10];
output->z4[11] = entropyMPC_z03[11];
output->z5[0] = entropyMPC_z04[0];
output->z5[1] = entropyMPC_z04[1];
output->z5[2] = entropyMPC_z04[2];
output->z5[3] = entropyMPC_z04[3];
output->z5[4] = entropyMPC_z04[4];
output->z5[5] = entropyMPC_z04[5];
output->z5[6] = entropyMPC_z04[6];
output->z5[7] = entropyMPC_z04[7];
output->z5[8] = entropyMPC_z04[8];
output->z5[9] = entropyMPC_z04[9];
output->z5[10] = entropyMPC_z04[10];
output->z5[11] = entropyMPC_z04[11];
output->z6[0] = entropyMPC_z05[0];
output->z6[1] = entropyMPC_z05[1];
output->z6[2] = entropyMPC_z05[2];
output->z6[3] = entropyMPC_z05[3];
output->z6[4] = entropyMPC_z05[4];
output->z6[5] = entropyMPC_z05[5];
output->z6[6] = entropyMPC_z05[6];
output->z6[7] = entropyMPC_z05[7];
output->z6[8] = entropyMPC_z05[8];
output->z6[9] = entropyMPC_z05[9];
output->z6[10] = entropyMPC_z05[10];
output->z6[11] = entropyMPC_z05[11];
output->z7[0] = entropyMPC_z06[0];
output->z7[1] = entropyMPC_z06[1];
output->z7[2] = entropyMPC_z06[2];
output->z7[3] = entropyMPC_z06[3];
output->z7[4] = entropyMPC_z06[4];
output->z7[5] = entropyMPC_z06[5];
output->z7[6] = entropyMPC_z06[6];
output->z7[7] = entropyMPC_z06[7];
output->z7[8] = entropyMPC_z06[8];
output->z7[9] = entropyMPC_z06[9];
output->z7[10] = entropyMPC_z06[10];
output->z7[11] = entropyMPC_z06[11];
output->z8[0] = entropyMPC_z07[0];
output->z8[1] = entropyMPC_z07[1];
output->z8[2] = entropyMPC_z07[2];
output->z8[3] = entropyMPC_z07[3];
output->z8[4] = entropyMPC_z07[4];
output->z8[5] = entropyMPC_z07[5];
output->z8[6] = entropyMPC_z07[6];
output->z8[7] = entropyMPC_z07[7];
output->z8[8] = entropyMPC_z07[8];
output->z8[9] = entropyMPC_z07[9];
output->z8[10] = entropyMPC_z07[10];
output->z8[11] = entropyMPC_z07[11];
output->z9[0] = entropyMPC_z08[0];
output->z9[1] = entropyMPC_z08[1];
output->z9[2] = entropyMPC_z08[2];
output->z9[3] = entropyMPC_z08[3];
output->z9[4] = entropyMPC_z08[4];
output->z9[5] = entropyMPC_z08[5];
output->z9[6] = entropyMPC_z08[6];
output->z9[7] = entropyMPC_z08[7];
output->z9[8] = entropyMPC_z08[8];
output->z9[9] = entropyMPC_z08[9];
output->z9[10] = entropyMPC_z08[10];
output->z9[11] = entropyMPC_z08[11];
output->z10[0] = entropyMPC_z09[0];
output->z10[1] = entropyMPC_z09[1];
output->z10[2] = entropyMPC_z09[2];
output->z10[3] = entropyMPC_z09[3];
output->z10[4] = entropyMPC_z09[4];
output->z10[5] = entropyMPC_z09[5];
output->z10[6] = entropyMPC_z09[6];
output->z10[7] = entropyMPC_z09[7];
output->z10[8] = entropyMPC_z09[8];
output->z10[9] = entropyMPC_z09[9];
output->z10[10] = entropyMPC_z09[10];
output->z10[11] = entropyMPC_z09[11];
output->z11[0] = entropyMPC_z10[0];
output->z11[1] = entropyMPC_z10[1];
output->z11[2] = entropyMPC_z10[2];
output->z11[3] = entropyMPC_z10[3];
output->z11[4] = entropyMPC_z10[4];
output->z11[5] = entropyMPC_z10[5];
output->z11[6] = entropyMPC_z10[6];
output->z11[7] = entropyMPC_z10[7];
output->z11[8] = entropyMPC_z10[8];
output->z11[9] = entropyMPC_z10[9];
output->z11[10] = entropyMPC_z10[10];
output->z11[11] = entropyMPC_z10[11];
output->z12[0] = entropyMPC_z11[0];
output->z12[1] = entropyMPC_z11[1];
output->z12[2] = entropyMPC_z11[2];
output->z12[3] = entropyMPC_z11[3];
output->z12[4] = entropyMPC_z11[4];
output->z12[5] = entropyMPC_z11[5];
output->z12[6] = entropyMPC_z11[6];
output->z12[7] = entropyMPC_z11[7];
output->z12[8] = entropyMPC_z11[8];
output->z12[9] = entropyMPC_z11[9];
output->z12[10] = entropyMPC_z11[10];
output->z12[11] = entropyMPC_z11[11];
output->z13[0] = entropyMPC_z12[0];
output->z13[1] = entropyMPC_z12[1];
output->z13[2] = entropyMPC_z12[2];
output->z13[3] = entropyMPC_z12[3];
output->z13[4] = entropyMPC_z12[4];
output->z13[5] = entropyMPC_z12[5];
output->z13[6] = entropyMPC_z12[6];
output->z13[7] = entropyMPC_z12[7];
output->z13[8] = entropyMPC_z12[8];
output->z13[9] = entropyMPC_z12[9];
output->z13[10] = entropyMPC_z12[10];
output->z13[11] = entropyMPC_z12[11];
output->z14[0] = entropyMPC_z13[0];
output->z14[1] = entropyMPC_z13[1];
output->z14[2] = entropyMPC_z13[2];
output->z14[3] = entropyMPC_z13[3];
output->z14[4] = entropyMPC_z13[4];
output->z14[5] = entropyMPC_z13[5];
output->z14[6] = entropyMPC_z13[6];
output->z14[7] = entropyMPC_z13[7];
output->z14[8] = entropyMPC_z13[8];
output->z14[9] = entropyMPC_z13[9];
output->z14[10] = entropyMPC_z13[10];
output->z14[11] = entropyMPC_z13[11];
output->z15[0] = entropyMPC_z14[0];
output->z15[1] = entropyMPC_z14[1];
output->z15[2] = entropyMPC_z14[2];
output->z15[3] = entropyMPC_z14[3];
output->z15[4] = entropyMPC_z14[4];
output->z15[5] = entropyMPC_z14[5];
output->z15[6] = entropyMPC_z14[6];
output->z15[7] = entropyMPC_z14[7];
output->z15[8] = entropyMPC_z14[8];
output->z15[9] = entropyMPC_z14[9];

#if entropyMPC_SET_TIMING == 1
info->solvetime = entropyMPC_toc(&solvertimer);
#if entropyMPC_SET_PRINTLEVEL > 0 && entropyMPC_SET_TIMING == 1
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
