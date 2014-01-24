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

#include "lpMPC.h"

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
 * Initializes a vector of length 58 with a value.
 */
void lpMPC_LA_INITIALIZEVECTOR_58(lpMPC_FLOAT* vec, lpMPC_FLOAT value)
{
	int i;
	for( i=0; i<58; i++ )
	{
		vec[i] = value;
	}
}


/*
 * Initializes a vector of length 30 with a value.
 */
void lpMPC_LA_INITIALIZEVECTOR_30(lpMPC_FLOAT* vec, lpMPC_FLOAT value)
{
	int i;
	for( i=0; i<30; i++ )
	{
		vec[i] = value;
	}
}


/*
 * Initializes a vector of length 116 with a value.
 */
void lpMPC_LA_INITIALIZEVECTOR_116(lpMPC_FLOAT* vec, lpMPC_FLOAT value)
{
	int i;
	for( i=0; i<116; i++ )
	{
		vec[i] = value;
	}
}


/* 
 * Calculates a dot product and adds it to a variable: z += x'*y; 
 * This function is for vectors of length 116.
 */
void lpMPC_LA_DOTACC_116(lpMPC_FLOAT *x, lpMPC_FLOAT *y, lpMPC_FLOAT *z)
{
	int i;
	for( i=0; i<116; i++ ){
		*z += x[i]*y[i];
	}
}


/*
 * Calculates the gradient and the value for a quadratic function 0.5*z'*H*z + f'*z
 *
 * INPUTS:     H  - Symmetric Hessian, diag matrix of size [4 x 4]
 *             f  - column vector of size 4
 *             z  - column vector of size 4
 *
 * OUTPUTS: grad  - gradient at z (= H*z + f), column vector of size 4
 *          value <-- value + 0.5*z'*H*z + f'*z (value will be modified)
 */
void lpMPC_LA_DIAG_QUADFCN_4(lpMPC_FLOAT* H, lpMPC_FLOAT* f, lpMPC_FLOAT* z, lpMPC_FLOAT* grad, lpMPC_FLOAT* value)
{
	int i;
	lpMPC_FLOAT hz;	
	for( i=0; i<4; i++){
		hz = H[i]*z[i];
		grad[i] = hz + f[i];
		*value += 0.5*hz*z[i] + f[i]*z[i];
	}
}


/*
 * Calculates the gradient and the value for a quadratic function 0.5*z'*H*z + f'*z
 *
 * INPUTS:     H  - Symmetric Hessian, diag matrix of size [2 x 2]
 *             f  - column vector of size 2
 *             z  - column vector of size 2
 *
 * OUTPUTS: grad  - gradient at z (= H*z + f), column vector of size 2
 *          value <-- value + 0.5*z'*H*z + f'*z (value will be modified)
 */
void lpMPC_LA_DIAG_QUADFCN_2(lpMPC_FLOAT* H, lpMPC_FLOAT* f, lpMPC_FLOAT* z, lpMPC_FLOAT* grad, lpMPC_FLOAT* value)
{
	int i;
	lpMPC_FLOAT hz;	
	for( i=0; i<2; i++){
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
void lpMPC_LA_DIAGZERO_MVMSUB6_2(lpMPC_FLOAT *B, lpMPC_FLOAT *u, lpMPC_FLOAT *b, lpMPC_FLOAT *l, lpMPC_FLOAT *r, lpMPC_FLOAT *z, lpMPC_FLOAT *y)
{
	int i;
	lpMPC_FLOAT Bu[2];
	lpMPC_FLOAT norm = *y;
	lpMPC_FLOAT lr = 0;

	/* do A*x + B*u first */
	for( i=0; i<2; i++ ){
		Bu[i] = B[i]*u[i];
	}	

	for( i=0; i<2; i++ ){
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
void lpMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_4_4(lpMPC_FLOAT *A, lpMPC_FLOAT *x, lpMPC_FLOAT *B, lpMPC_FLOAT *u, lpMPC_FLOAT *b, lpMPC_FLOAT *l, lpMPC_FLOAT *r, lpMPC_FLOAT *z, lpMPC_FLOAT *y)
{
	int i;
	int j;
	int k = 0;
	lpMPC_FLOAT AxBu[2];
	lpMPC_FLOAT norm = *y;
	lpMPC_FLOAT lr = 0;

	/* do A*x + B*u first */
	for( i=0; i<2; i++ ){
		AxBu[i] = A[k++]*x[0] + B[i]*u[i];
	}	

	for( j=1; j<4; j++ ){		
		for( i=0; i<2; i++ ){
			AxBu[i] += A[k++]*x[j];
		}
	}

	for( i=0; i<2; i++ ){
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
void lpMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_4_2(lpMPC_FLOAT *A, lpMPC_FLOAT *x, lpMPC_FLOAT *B, lpMPC_FLOAT *u, lpMPC_FLOAT *b, lpMPC_FLOAT *l, lpMPC_FLOAT *r, lpMPC_FLOAT *z, lpMPC_FLOAT *y)
{
	int i;
	int j;
	int k = 0;
	lpMPC_FLOAT AxBu[2];
	lpMPC_FLOAT norm = *y;
	lpMPC_FLOAT lr = 0;

	/* do A*x + B*u first */
	for( i=0; i<2; i++ ){
		AxBu[i] = A[k++]*x[0] + B[i]*u[i];
	}	

	for( j=1; j<4; j++ ){		
		for( i=0; i<2; i++ ){
			AxBu[i] += A[k++]*x[j];
		}
	}

	for( i=0; i<2; i++ ){
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
 * where A is of size [2 x 4] and stored in column major format.
 * and B is of size [2 x 4] and stored in diagzero format
 * Note the transposes of A and B!
 */
void lpMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(lpMPC_FLOAT *A, lpMPC_FLOAT *x, lpMPC_FLOAT *B, lpMPC_FLOAT *y, lpMPC_FLOAT *z)
{
	int i;
	int j;
	int k = 0;
	for( i=0; i<2; i++ ){
		z[i] = 0;
		for( j=0; j<2; j++ ){
			z[i] += A[k++]*x[j];
		}
		z[i] += B[i]*y[i];
	}
	for( i=2 ;i<4; i++ ){
		z[i] = 0;
		for( j=0; j<2; j++ ){
			z[i] += A[k++]*x[j];
		}
	}
}


/*
 * Matrix vector multiplication y = M'*x where M is of size [2 x 2]
 * and stored in diagzero format. Note the transpose of M!
 */
void lpMPC_LA_DIAGZERO_MTVM_2_2(lpMPC_FLOAT *M, lpMPC_FLOAT *x, lpMPC_FLOAT *y)
{
	int i;
	for( i=0; i<2; i++ ){
		y[i] = M[i]*x[i];
	}
}


/*
 * Vector subtraction and addition.
 *	 Input: five vectors t, tidx, u, v, w and two scalars z and r
 *	 Output: y = t(tidx) - u + w
 *           z = z - v'*x;
 *           r = max([norm(y,inf), z]);
 * for vectors of length 4. Output z is of course scalar.
 */
void lpMPC_LA_VSUBADD3_4(lpMPC_FLOAT* t, lpMPC_FLOAT* u, int* uidx, lpMPC_FLOAT* v, lpMPC_FLOAT* w, lpMPC_FLOAT* y, lpMPC_FLOAT* z, lpMPC_FLOAT* r)
{
	int i;
	lpMPC_FLOAT norm = *r;
	lpMPC_FLOAT vx = 0;
	lpMPC_FLOAT x;
	for( i=0; i<4; i++){
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
 * for vectors of length 4. Output z is of course scalar.
 */
void lpMPC_LA_VSUBADD2_4(lpMPC_FLOAT* t, int* tidx, lpMPC_FLOAT* u, lpMPC_FLOAT* v, lpMPC_FLOAT* w, lpMPC_FLOAT* y, lpMPC_FLOAT* z, lpMPC_FLOAT* r)
{
	int i;
	lpMPC_FLOAT norm = *r;
	lpMPC_FLOAT vx = 0;
	lpMPC_FLOAT x;
	for( i=0; i<4; i++){
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
 * for vectors of length 2. Output z is of course scalar.
 */
void lpMPC_LA_VSUBADD3_2(lpMPC_FLOAT* t, lpMPC_FLOAT* u, int* uidx, lpMPC_FLOAT* v, lpMPC_FLOAT* w, lpMPC_FLOAT* y, lpMPC_FLOAT* z, lpMPC_FLOAT* r)
{
	int i;
	lpMPC_FLOAT norm = *r;
	lpMPC_FLOAT vx = 0;
	lpMPC_FLOAT x;
	for( i=0; i<2; i++){
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
 * for vectors of length 2. Output z is of course scalar.
 */
void lpMPC_LA_VSUBADD2_2(lpMPC_FLOAT* t, int* tidx, lpMPC_FLOAT* u, lpMPC_FLOAT* v, lpMPC_FLOAT* w, lpMPC_FLOAT* y, lpMPC_FLOAT* z, lpMPC_FLOAT* r)
{
	int i;
	lpMPC_FLOAT norm = *r;
	lpMPC_FLOAT vx = 0;
	lpMPC_FLOAT x;
	for( i=0; i<2; i++){
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
 * Special function for box constraints of length 4
 * Returns also L/S, a value that is often used elsewhere.
 */
void lpMPC_LA_INEQ_B_GRAD_4_4_4(lpMPC_FLOAT *lu, lpMPC_FLOAT *su, lpMPC_FLOAT *ru, lpMPC_FLOAT *ll, lpMPC_FLOAT *sl, lpMPC_FLOAT *rl, int* lbIdx, int* ubIdx, lpMPC_FLOAT *grad, lpMPC_FLOAT *lubysu, lpMPC_FLOAT *llbysl)
{
	int i;
	for( i=0; i<4; i++ ){
		grad[i] = 0;
	}
	for( i=0; i<4; i++ ){		
		llbysl[i] = ll[i] / sl[i];
		grad[lbIdx[i]] -= llbysl[i]*rl[i];
	}
	for( i=0; i<4; i++ ){
		lubysu[i] = lu[i] / su[i];
		grad[ubIdx[i]] += lubysu[i]*ru[i];
	}
}


/*
 * Computes inequality constraints gradient-
 * Special function for box constraints of length 2
 * Returns also L/S, a value that is often used elsewhere.
 */
void lpMPC_LA_INEQ_B_GRAD_2_2_2(lpMPC_FLOAT *lu, lpMPC_FLOAT *su, lpMPC_FLOAT *ru, lpMPC_FLOAT *ll, lpMPC_FLOAT *sl, lpMPC_FLOAT *rl, int* lbIdx, int* ubIdx, lpMPC_FLOAT *grad, lpMPC_FLOAT *lubysu, lpMPC_FLOAT *llbysl)
{
	int i;
	for( i=0; i<2; i++ ){
		grad[i] = 0;
	}
	for( i=0; i<2; i++ ){		
		llbysl[i] = ll[i] / sl[i];
		grad[lbIdx[i]] -= llbysl[i]*rl[i];
	}
	for( i=0; i<2; i++ ){
		lubysu[i] = lu[i] / su[i];
		grad[ubIdx[i]] += lubysu[i]*ru[i];
	}
}


/*
 * Addition of three vectors  z = u + w + v
 * of length 58.
 */
void lpMPC_LA_VVADD3_58(lpMPC_FLOAT *u, lpMPC_FLOAT *v, lpMPC_FLOAT *w, lpMPC_FLOAT *z)
{
	int i;
	for( i=0; i<58; i++ ){
		z[i] = u[i] + v[i] + w[i];
	}
}


/*
 * Special function to compute the diagonal cholesky factorization of the 
 * positive definite augmented Hessian for block size 4.
 *
 * Inputs: - H = diagonal cost Hessian in diagonal storage format
 *         - llbysl = L / S of lower bounds
 *         - lubysu = L / S of upper bounds
 *
 * Output: Phi = sqrt(H + diag(llbysl) + diag(lubysu))
 * where Phi is stored in diagonal storage format
 */
void lpMPC_LA_DIAG_CHOL_ONELOOP_LBUB_4_4_4(lpMPC_FLOAT *H, lpMPC_FLOAT *llbysl, int* lbIdx, lpMPC_FLOAT *lubysu, int* ubIdx, lpMPC_FLOAT *Phi)


{
	int i;
	
	/* compute cholesky */
	for( i=0; i<4; i++ ){
		Phi[i] = H[i] + llbysl[i] + lubysu[i];

#if lpMPC_SET_PRINTLEVEL > 0 && defined PRINTNUMERICALWARNINGS
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
 * where A is to be computed and is of size [2 x 4],
 * B is given and of size [2 x 4], L is a diagonal
 * matrix of size 2 stored in diagonal matrix 
 * storage format. Note the transpose of L has no impact!
 *
 * Result: A in column major storage format.
 *
 */
void lpMPC_LA_DIAG_MATRIXFORWARDSUB_2_4(lpMPC_FLOAT *L, lpMPC_FLOAT *B, lpMPC_FLOAT *A)
{
    int i,j;
	 int k = 0;

	for( j=0; j<4; j++){
		for( i=0; i<2; i++){
			A[k] = B[k]/L[j];
			k++;
		}
	}

}


/**
 * Forward substitution for the matrix equation A*L' = B
 * where A is to be computed and is of size [2 x 4],
 * B is given and of size [2 x 4], L is a diagonal
 *  matrix of size 4 stored in diagonal 
 * storage format. Note the transpose of L!
 *
 * Result: A in diagonalzero storage format.
 *
 */
void lpMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_4(lpMPC_FLOAT *L, lpMPC_FLOAT *B, lpMPC_FLOAT *A)
{
	int j;
    for( j=0; j<4; j++ ){   
		A[j] = B[j]/L[j];
     }
}


/**
 * Compute C = A*B' where 
 *
 *	size(A) = [2 x 4]
 *  size(B) = [2 x 4] in diagzero format
 * 
 * A and C matrices are stored in column major format.
 * 
 * 
 */
void lpMPC_LA_DENSE_DIAGZERO_MMTM_2_4_2(lpMPC_FLOAT *A, lpMPC_FLOAT *B, lpMPC_FLOAT *C)
{
    int i, j;
	
	for( i=0; i<2; i++ ){
		for( j=0; j<2; j++){
			C[j*2+i] = B[i*2+j]*A[i];
		}
	}

}


/**
 * Forward substitution to solve L*y = b where L is a
 * diagonal matrix in vector storage format.
 * 
 * The dimensions involved are 4.
 */
void lpMPC_LA_DIAG_FORWARDSUB_4(lpMPC_FLOAT *L, lpMPC_FLOAT *b, lpMPC_FLOAT *y)
{
    int i;

    for( i=0; i<4; i++ ){
		y[i] = b[i]/L[i];
    }
}


/*
 * Special function to compute the diagonal cholesky factorization of the 
 * positive definite augmented Hessian for block size 2.
 *
 * Inputs: - H = diagonal cost Hessian in diagonal storage format
 *         - llbysl = L / S of lower bounds
 *         - lubysu = L / S of upper bounds
 *
 * Output: Phi = sqrt(H + diag(llbysl) + diag(lubysu))
 * where Phi is stored in diagonal storage format
 */
void lpMPC_LA_DIAG_CHOL_ONELOOP_LBUB_2_2_2(lpMPC_FLOAT *H, lpMPC_FLOAT *llbysl, int* lbIdx, lpMPC_FLOAT *lubysu, int* ubIdx, lpMPC_FLOAT *Phi)


{
	int i;
	
	/* compute cholesky */
	for( i=0; i<2; i++ ){
		Phi[i] = H[i] + llbysl[i] + lubysu[i];

#if lpMPC_SET_PRINTLEVEL > 0 && defined PRINTNUMERICALWARNINGS
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
 * where A is to be computed and is of size [2 x 2],
 * B is given and of size [2 x 2], L is a diagonal
 *  matrix of size 2 stored in diagonal 
 * storage format. Note the transpose of L!
 *
 * Result: A in diagonalzero storage format.
 *
 */
void lpMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_2(lpMPC_FLOAT *L, lpMPC_FLOAT *B, lpMPC_FLOAT *A)
{
	int j;
    for( j=0; j<2; j++ ){   
		A[j] = B[j]/L[j];
     }
}


/**
 * Forward substitution to solve L*y = b where L is a
 * diagonal matrix in vector storage format.
 * 
 * The dimensions involved are 2.
 */
void lpMPC_LA_DIAG_FORWARDSUB_2(lpMPC_FLOAT *L, lpMPC_FLOAT *b, lpMPC_FLOAT *y)
{
    int i;

    for( i=0; i<2; i++ ){
		y[i] = b[i]/L[i];
    }
}


/**
 * Compute L = B*B', where L is lower triangular of size NXp1
 * and B is a diagzero matrix of size [2 x 4] in column
 * storage format.
 * 
 */
void lpMPC_LA_DIAGZERO_MMT_2(lpMPC_FLOAT *B, lpMPC_FLOAT *L)
{
    int i, ii, di;
    
    ii = 0; di = 0;
    for( i=0; i<2; i++ ){        
		L[ii+i] = B[i]*B[i];
        ii += ++di;
    }
}


/* 
 * Computes r = b - B*u
 * B is stored in diagzero format
 */
void lpMPC_LA_DIAGZERO_MVMSUB7_2(lpMPC_FLOAT *B, lpMPC_FLOAT *u, lpMPC_FLOAT *b, lpMPC_FLOAT *r)
{
	int i;

	for( i=0; i<2; i++ ){
		r[i] = b[i] - B[i]*u[i];
	}	
	
}


/**
 * Compute L = A*A' + B*B', where L is lower triangular of size NXp1
 * and A is a dense matrix of size [2 x 4] in column
 * storage format, and B is of size [2 x 4] diagonalzero
 * storage format.
 * 
 * THIS ONE HAS THE WORST ACCES PATTERN POSSIBLE. 
 * POSSIBKE FIX: PUT A AND B INTO ROW MAJOR FORMAT FIRST.
 * 
 */
void lpMPC_LA_DENSE_DIAGZERO_MMT2_2_4_4(lpMPC_FLOAT *A, lpMPC_FLOAT *B, lpMPC_FLOAT *L)
{
    int i, j, k, ii, di;
    lpMPC_FLOAT ltemp;
    
    ii = 0; di = 0;
    for( i=0; i<2; i++ ){        
        for( j=0; j<=i; j++ ){
            ltemp = 0; 
            for( k=0; k<4; k++ ){
                ltemp += A[k*2+i]*A[k*2+j];
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
void lpMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_4_4(lpMPC_FLOAT *A, lpMPC_FLOAT *x, lpMPC_FLOAT *B, lpMPC_FLOAT *u, lpMPC_FLOAT *b, lpMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;

	for( i=0; i<2; i++ ){
		r[i] = b[i] - A[k++]*x[0] - B[i]*u[i];
	}	

	for( j=1; j<4; j++ ){		
		for( i=0; i<2; i++ ){
			r[i] -= A[k++]*x[j];
		}
	}
	
}


/**
 * Compute L = A*A' + B*B', where L is lower triangular of size NXp1
 * and A is a dense matrix of size [2 x 4] in column
 * storage format, and B is of size [2 x 2] diagonalzero
 * storage format.
 * 
 * THIS ONE HAS THE WORST ACCES PATTERN POSSIBLE. 
 * POSSIBKE FIX: PUT A AND B INTO ROW MAJOR FORMAT FIRST.
 * 
 */
void lpMPC_LA_DENSE_DIAGZERO_MMT2_2_4_2(lpMPC_FLOAT *A, lpMPC_FLOAT *B, lpMPC_FLOAT *L)
{
    int i, j, k, ii, di;
    lpMPC_FLOAT ltemp;
    
    ii = 0; di = 0;
    for( i=0; i<2; i++ ){        
        for( j=0; j<=i; j++ ){
            ltemp = 0; 
            for( k=0; k<4; k++ ){
                ltemp += A[k*2+i]*A[k*2+j];
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
void lpMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_4_2(lpMPC_FLOAT *A, lpMPC_FLOAT *x, lpMPC_FLOAT *B, lpMPC_FLOAT *u, lpMPC_FLOAT *b, lpMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;

	for( i=0; i<2; i++ ){
		r[i] = b[i] - A[k++]*x[0] - B[i]*u[i];
	}	

	for( j=1; j<4; j++ ){		
		for( i=0; i<2; i++ ){
			r[i] -= A[k++]*x[j];
		}
	}
	
}


/**
 * Cholesky factorization as above, but working on a matrix in 
 * lower triangular storage format of size 2 and outputting
 * the Cholesky factor to matrix L in lower triangular format.
 */
void lpMPC_LA_DENSE_CHOL_2(lpMPC_FLOAT *A, lpMPC_FLOAT *L)
{
    int i, j, k, di, dj;
	 int ii, jj;

    lpMPC_FLOAT l;
    lpMPC_FLOAT Mii;

	/* copy A to L first and then operate on L */
	/* COULD BE OPTIMIZED */
	ii=0; di=0;
	for( i=0; i<2; i++ ){
		for( j=0; j<=i; j++ ){
			L[ii+j] = A[ii+j];
		}
		ii += ++di;
	}    
	
	/* factor L */
	ii=0; di=0;
    for( i=0; i<2; i++ ){
        l = 0;
        for( k=0; k<i; k++ ){
            l += L[ii+k]*L[ii+k];
        }        
        
        Mii = L[ii+i] - l;
        
#if lpMPC_SET_PRINTLEVEL > 0 && defined PRINTNUMERICALWARNINGS
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
        for( j=i+1; j<2; j++ ){
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
 * The dimensions involved are 2.
 */
void lpMPC_LA_DENSE_FORWARDSUB_2(lpMPC_FLOAT *L, lpMPC_FLOAT *b, lpMPC_FLOAT *y)
{
    int i,j,ii,di;
    lpMPC_FLOAT yel;
            
    ii = 0; di = 0;
    for( i=0; i<2; i++ ){
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
 * where A is to be computed and is of size [2 x 2],
 * B is given and of size [2 x 2], L is a lower tri-
 * angular matrix of size 2 stored in lower triangular 
 * storage format. Note the transpose of L AND B!
 *
 * Result: A in column major storage format.
 *
 */
void lpMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(lpMPC_FLOAT *L, lpMPC_FLOAT *B, lpMPC_FLOAT *A)
{
    int i,j,k,ii,di;
    lpMPC_FLOAT a;
    
    ii=0; di=0;
    for( j=0; j<2; j++ ){        
        for( i=0; i<2; i++ ){
            a = B[i*2+j];
            for( k=0; k<j; k++ ){
                a -= A[k*2+i]*L[ii+k];
            }    

			/* saturate for numerical stability */
			a = MIN(a, BIGM);
			a = MAX(a, -BIGM); 

			A[j*2+i] = a/L[ii+j];			
        }
        ii += ++di;
    }
}


/**
 * Compute L = L - A*A', where L is lower triangular of size 2
 * and A is a dense matrix of size [2 x 2] in column
 * storage format.
 * 
 * THIS ONE HAS THE WORST ACCES PATTERN POSSIBLE. 
 * POSSIBKE FIX: PUT A INTO ROW MAJOR FORMAT FIRST.
 * 
 */
void lpMPC_LA_DENSE_MMTSUB_2_2(lpMPC_FLOAT *A, lpMPC_FLOAT *L)
{
    int i, j, k, ii, di;
    lpMPC_FLOAT ltemp;
    
    ii = 0; di = 0;
    for( i=0; i<2; i++ ){        
        for( j=0; j<=i; j++ ){
            ltemp = 0; 
            for( k=0; k<2; k++ ){
                ltemp += A[k*2+i]*A[k*2+j];
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
void lpMPC_LA_DENSE_MVMSUB1_2_2(lpMPC_FLOAT *A, lpMPC_FLOAT *x, lpMPC_FLOAT *b, lpMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;

	for( i=0; i<2; i++ ){
		r[i] = b[i] - A[k++]*x[0];
	}	
	for( j=1; j<2; j++ ){		
		for( i=0; i<2; i++ ){
			r[i] -= A[k++]*x[j];
		}
	}
}


/**
 * Backward Substitution to solve L^T*x = y where L is a
 * lower triangular matrix in triangular storage format.
 * 
 * All involved dimensions are 2.
 */
void lpMPC_LA_DENSE_BACKWARDSUB_2(lpMPC_FLOAT *L, lpMPC_FLOAT *y, lpMPC_FLOAT *x)
{
    int i, ii, di, j, jj, dj;
    lpMPC_FLOAT xel;    
	int start = 1;
    
    /* now solve L^T*x = y by backward substitution */
    ii = start; di = 1;
    for( i=1; i>=0; i-- ){        
        xel = y[i];        
        jj = start; dj = 1;
        for( j=1; j>i; j-- ){
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
 * Matrix vector multiplication y = b - M'*x where M is of size [2 x 2]
 * and stored in column major format. Note the transpose of M!
 */
void lpMPC_LA_DENSE_MTVMSUB_2_2(lpMPC_FLOAT *A, lpMPC_FLOAT *x, lpMPC_FLOAT *b, lpMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0; 
	for( i=0; i<2; i++ ){
		r[i] = b[i];
		for( j=0; j<2; j++ ){
			r[i] -= A[k++]*x[j];
		}
	}
}


/*
 * Vector subtraction z = -x - y for vectors of length 58.
 */
void lpMPC_LA_VSUB2_58(lpMPC_FLOAT *x, lpMPC_FLOAT *y, lpMPC_FLOAT *z)
{
	int i;
	for( i=0; i<58; i++){
		z[i] = -x[i] - y[i];
	}
}


/**
 * Forward-Backward-Substitution to solve L*L^T*x = b where L is a
 * diagonal matrix of size 4 in vector
 * storage format.
 */
void lpMPC_LA_DIAG_FORWARDBACKWARDSUB_4(lpMPC_FLOAT *L, lpMPC_FLOAT *b, lpMPC_FLOAT *x)
{
    int i;
            
    /* solve Ly = b by forward and backward substitution */
    for( i=0; i<4; i++ ){
		x[i] = b[i]/(L[i]*L[i]);
    }
    
}


/**
 * Forward-Backward-Substitution to solve L*L^T*x = b where L is a
 * diagonal matrix of size 2 in vector
 * storage format.
 */
void lpMPC_LA_DIAG_FORWARDBACKWARDSUB_2(lpMPC_FLOAT *L, lpMPC_FLOAT *b, lpMPC_FLOAT *x)
{
    int i;
            
    /* solve Ly = b by forward and backward substitution */
    for( i=0; i<2; i++ ){
		x[i] = b[i]/(L[i]*L[i]);
    }
    
}


/*
 * Vector subtraction z = x(xidx) - y where y, z and xidx are of length 4,
 * and x has length 4 and is indexed through yidx.
 */
void lpMPC_LA_VSUB_INDEXED_4(lpMPC_FLOAT *x, int* xidx, lpMPC_FLOAT *y, lpMPC_FLOAT *z)
{
	int i;
	for( i=0; i<4; i++){
		z[i] = x[xidx[i]] - y[i];
	}
}


/*
 * Vector subtraction x = -u.*v - w for vectors of length 4.
 */
void lpMPC_LA_VSUB3_4(lpMPC_FLOAT *u, lpMPC_FLOAT *v, lpMPC_FLOAT *w, lpMPC_FLOAT *x)
{
	int i;
	for( i=0; i<4; i++){
		x[i] = -u[i]*v[i] - w[i];
	}
}


/*
 * Vector subtraction z = -x - y(yidx) where y is of length 4
 * and z, x and yidx are of length 4.
 */
void lpMPC_LA_VSUB2_INDEXED_4(lpMPC_FLOAT *x, lpMPC_FLOAT *y, int* yidx, lpMPC_FLOAT *z)
{
	int i;
	for( i=0; i<4; i++){
		z[i] = -x[i] - y[yidx[i]];
	}
}


/*
 * Vector subtraction z = x(xidx) - y where y, z and xidx are of length 2,
 * and x has length 2 and is indexed through yidx.
 */
void lpMPC_LA_VSUB_INDEXED_2(lpMPC_FLOAT *x, int* xidx, lpMPC_FLOAT *y, lpMPC_FLOAT *z)
{
	int i;
	for( i=0; i<2; i++){
		z[i] = x[xidx[i]] - y[i];
	}
}


/*
 * Vector subtraction x = -u.*v - w for vectors of length 2.
 */
void lpMPC_LA_VSUB3_2(lpMPC_FLOAT *u, lpMPC_FLOAT *v, lpMPC_FLOAT *w, lpMPC_FLOAT *x)
{
	int i;
	for( i=0; i<2; i++){
		x[i] = -u[i]*v[i] - w[i];
	}
}


/*
 * Vector subtraction z = -x - y(yidx) where y is of length 2
 * and z, x and yidx are of length 2.
 */
void lpMPC_LA_VSUB2_INDEXED_2(lpMPC_FLOAT *x, lpMPC_FLOAT *y, int* yidx, lpMPC_FLOAT *z)
{
	int i;
	for( i=0; i<2; i++){
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
 * lpMPC_NOPROGRESS (should be negative).
 */
int lpMPC_LINESEARCH_BACKTRACKING_AFFINE(lpMPC_FLOAT *l, lpMPC_FLOAT *s, lpMPC_FLOAT *dl, lpMPC_FLOAT *ds, lpMPC_FLOAT *a, lpMPC_FLOAT *mu_aff)
{
    int i;
	int lsIt=1;    
    lpMPC_FLOAT dltemp;
    lpMPC_FLOAT dstemp;
    lpMPC_FLOAT mya = 1.0;
    lpMPC_FLOAT mymu;
        
    while( 1 ){                        

        /* 
         * Compute both snew and wnew together.
         * We compute also mu_affine along the way here, as the
         * values might be in registers, so it should be cheaper.
         */
        mymu = 0;
        for( i=0; i<116; i++ ){
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
        if( i == 116 ){
            break;
        } else {
            mya *= lpMPC_SET_LS_SCALE_AFF;
            if( mya < lpMPC_SET_LS_MINSTEP ){
                return lpMPC_NOPROGRESS;
            }
        }
    }
    
    /* return new values and iteration counter */
    *a = mya;
    *mu_aff = mymu / (lpMPC_FLOAT)116;
    return lsIt;
}


/*
 * Vector subtraction x = u.*v - a where a is a scalar
*  and x,u,v are vectors of length 116.
 */
void lpMPC_LA_VSUB5_116(lpMPC_FLOAT *u, lpMPC_FLOAT *v, lpMPC_FLOAT a, lpMPC_FLOAT *x)
{
	int i;
	for( i=0; i<116; i++){
		x[i] = u[i]*v[i] - a;
	}
}


/*
 * Computes x=0; x(uidx) += u/su; x(vidx) -= v/sv where x is of length 4,
 * u, su, uidx are of length 4 and v, sv, vidx are of length 4.
 */
void lpMPC_LA_VSUB6_INDEXED_4_4_4(lpMPC_FLOAT *u, lpMPC_FLOAT *su, int* uidx, lpMPC_FLOAT *v, lpMPC_FLOAT *sv, int* vidx, lpMPC_FLOAT *x)
{
	int i;
	for( i=0; i<4; i++ ){
		x[i] = 0;
	}
	for( i=0; i<4; i++){
		x[uidx[i]] += u[i]/su[i];
	}
	for( i=0; i<4; i++){
		x[vidx[i]] -= v[i]/sv[i];
	}
}


/* 
 * Computes r =  B*u
 * where B is stored in diagzero format
 */
void lpMPC_LA_DIAGZERO_MVM_2(lpMPC_FLOAT *B, lpMPC_FLOAT *u, lpMPC_FLOAT *r)
{
	int i;

	for( i=0; i<2; i++ ){
		r[i] = B[i]*u[i];
	}	
	
}


/* 
 * Computes r = A*x + B*u
 * where A is stored in column major format
 * and B is stored in diagzero format
 */
void lpMPC_LA_DENSE_DIAGZERO_2MVMADD_2_4_4(lpMPC_FLOAT *A, lpMPC_FLOAT *x, lpMPC_FLOAT *B, lpMPC_FLOAT *u, lpMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;

	for( i=0; i<2; i++ ){
		r[i] = A[k++]*x[0] + B[i]*u[i];
	}	

	for( j=1; j<4; j++ ){		
		for( i=0; i<2; i++ ){
			r[i] += A[k++]*x[j];
		}
	}
	
}


/*
 * Computes x=0; x(uidx) += u/su; x(vidx) -= v/sv where x is of length 2,
 * u, su, uidx are of length 2 and v, sv, vidx are of length 2.
 */
void lpMPC_LA_VSUB6_INDEXED_2_2_2(lpMPC_FLOAT *u, lpMPC_FLOAT *su, int* uidx, lpMPC_FLOAT *v, lpMPC_FLOAT *sv, int* vidx, lpMPC_FLOAT *x)
{
	int i;
	for( i=0; i<2; i++ ){
		x[i] = 0;
	}
	for( i=0; i<2; i++){
		x[uidx[i]] += u[i]/su[i];
	}
	for( i=0; i<2; i++){
		x[vidx[i]] -= v[i]/sv[i];
	}
}


/* 
 * Computes r = A*x + B*u
 * where A is stored in column major format
 * and B is stored in diagzero format
 */
void lpMPC_LA_DENSE_DIAGZERO_2MVMADD_2_4_2(lpMPC_FLOAT *A, lpMPC_FLOAT *x, lpMPC_FLOAT *B, lpMPC_FLOAT *u, lpMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;

	for( i=0; i<2; i++ ){
		r[i] = A[k++]*x[0] + B[i]*u[i];
	}	

	for( j=1; j<4; j++ ){		
		for( i=0; i<2; i++ ){
			r[i] += A[k++]*x[j];
		}
	}
	
}


/*
 * Vector subtraction z = x - y for vectors of length 58.
 */
void lpMPC_LA_VSUB_58(lpMPC_FLOAT *x, lpMPC_FLOAT *y, lpMPC_FLOAT *z)
{
	int i;
	for( i=0; i<58; i++){
		z[i] = x[i] - y[i];
	}
}


/** 
 * Computes z = -r./s - u.*y(y)
 * where all vectors except of y are of length 4 (length of y >= 4).
 */
void lpMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_4(lpMPC_FLOAT *r, lpMPC_FLOAT *s, lpMPC_FLOAT *u, lpMPC_FLOAT *y, int* yidx, lpMPC_FLOAT *z)
{
	int i;
	for( i=0; i<4; i++ ){
		z[i] = -r[i]/s[i] - u[i]*y[yidx[i]];
	}
}


/** 
 * Computes z = -r./s + u.*y(y)
 * where all vectors except of y are of length 4 (length of y >= 4).
 */
void lpMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_4(lpMPC_FLOAT *r, lpMPC_FLOAT *s, lpMPC_FLOAT *u, lpMPC_FLOAT *y, int* yidx, lpMPC_FLOAT *z)
{
	int i;
	for( i=0; i<4; i++ ){
		z[i] = -r[i]/s[i] + u[i]*y[yidx[i]];
	}
}


/** 
 * Computes z = -r./s - u.*y(y)
 * where all vectors except of y are of length 2 (length of y >= 2).
 */
void lpMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_2(lpMPC_FLOAT *r, lpMPC_FLOAT *s, lpMPC_FLOAT *u, lpMPC_FLOAT *y, int* yidx, lpMPC_FLOAT *z)
{
	int i;
	for( i=0; i<2; i++ ){
		z[i] = -r[i]/s[i] - u[i]*y[yidx[i]];
	}
}


/** 
 * Computes z = -r./s + u.*y(y)
 * where all vectors except of y are of length 2 (length of y >= 2).
 */
void lpMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_2(lpMPC_FLOAT *r, lpMPC_FLOAT *s, lpMPC_FLOAT *u, lpMPC_FLOAT *y, int* yidx, lpMPC_FLOAT *z)
{
	int i;
	for( i=0; i<2; i++ ){
		z[i] = -r[i]/s[i] + u[i]*y[yidx[i]];
	}
}


/*
 * Computes ds = -l.\(r + s.*dl) for vectors of length 116.
 */
void lpMPC_LA_VSUB7_116(lpMPC_FLOAT *l, lpMPC_FLOAT *r, lpMPC_FLOAT *s, lpMPC_FLOAT *dl, lpMPC_FLOAT *ds)
{
	int i;
	for( i=0; i<116; i++){
		ds[i] = -(r[i] + s[i]*dl[i])/l[i];
	}
}


/*
 * Vector addition x = x + y for vectors of length 58.
 */
void lpMPC_LA_VADD_58(lpMPC_FLOAT *x, lpMPC_FLOAT *y)
{
	int i;
	for( i=0; i<58; i++){
		x[i] += y[i];
	}
}


/*
 * Vector addition x = x + y for vectors of length 30.
 */
void lpMPC_LA_VADD_30(lpMPC_FLOAT *x, lpMPC_FLOAT *y)
{
	int i;
	for( i=0; i<30; i++){
		x[i] += y[i];
	}
}


/*
 * Vector addition x = x + y for vectors of length 116.
 */
void lpMPC_LA_VADD_116(lpMPC_FLOAT *x, lpMPC_FLOAT *y)
{
	int i;
	for( i=0; i<116; i++){
		x[i] += y[i];
	}
}


/**
 * Backtracking line search for combined predictor/corrector step.
 * Update on variables with safety factor gamma (to keep us away from
 * boundary).
 */
int lpMPC_LINESEARCH_BACKTRACKING_COMBINED(lpMPC_FLOAT *z, lpMPC_FLOAT *v, lpMPC_FLOAT *l, lpMPC_FLOAT *s, lpMPC_FLOAT *dz, lpMPC_FLOAT *dv, lpMPC_FLOAT *dl, lpMPC_FLOAT *ds, lpMPC_FLOAT *a, lpMPC_FLOAT *mu)
{
    int i, lsIt=1;       
    lpMPC_FLOAT dltemp;
    lpMPC_FLOAT dstemp;    
    lpMPC_FLOAT a_gamma;
            
    *a = 1.0;
    while( 1 ){                        

        /* check whether search criterion is fulfilled */
        for( i=0; i<116; i++ ){
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
        if( i == 116 ){
            break;
        } else {
            *a *= lpMPC_SET_LS_SCALE;
            if( *a < lpMPC_SET_LS_MINSTEP ){
                return lpMPC_NOPROGRESS;
            }
        }
    }
    
    /* update variables with safety margin */
    a_gamma = (*a)*lpMPC_SET_LS_MAXSTEP;
    
    /* primal variables */
    for( i=0; i<58; i++ ){
        z[i] += a_gamma*dz[i];
    }
    
    /* equality constraint multipliers */
    for( i=0; i<30; i++ ){
        v[i] += a_gamma*dv[i];
    }
    
    /* inequality constraint multipliers & slacks, also update mu */
    *mu = 0;
    for( i=0; i<116; i++ ){
        dltemp = l[i] + a_gamma*dl[i]; l[i] = dltemp;
        dstemp = s[i] + a_gamma*ds[i]; s[i] = dstemp;
        *mu += dltemp*dstemp;
    }
    
    *a = a_gamma;
    *mu /= (lpMPC_FLOAT)116;
    return lsIt;
}




/* VARIABLE DEFINITIONS ------------------------------------------------ */
lpMPC_FLOAT lpMPC_z[58];
lpMPC_FLOAT lpMPC_v[30];
lpMPC_FLOAT lpMPC_dz_aff[58];
lpMPC_FLOAT lpMPC_dv_aff[30];
lpMPC_FLOAT lpMPC_grad_cost[58];
lpMPC_FLOAT lpMPC_grad_eq[58];
lpMPC_FLOAT lpMPC_rd[58];
lpMPC_FLOAT lpMPC_l[116];
lpMPC_FLOAT lpMPC_s[116];
lpMPC_FLOAT lpMPC_lbys[116];
lpMPC_FLOAT lpMPC_dl_aff[116];
lpMPC_FLOAT lpMPC_ds_aff[116];
lpMPC_FLOAT lpMPC_dz_cc[58];
lpMPC_FLOAT lpMPC_dv_cc[30];
lpMPC_FLOAT lpMPC_dl_cc[116];
lpMPC_FLOAT lpMPC_ds_cc[116];
lpMPC_FLOAT lpMPC_ccrhs[116];
lpMPC_FLOAT lpMPC_grad_ineq[58];
lpMPC_FLOAT lpMPC_H00[4] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
lpMPC_FLOAT* lpMPC_z00 = lpMPC_z + 0;
lpMPC_FLOAT* lpMPC_dzaff00 = lpMPC_dz_aff + 0;
lpMPC_FLOAT* lpMPC_dzcc00 = lpMPC_dz_cc + 0;
lpMPC_FLOAT* lpMPC_rd00 = lpMPC_rd + 0;
lpMPC_FLOAT lpMPC_Lbyrd00[4];
lpMPC_FLOAT* lpMPC_grad_cost00 = lpMPC_grad_cost + 0;
lpMPC_FLOAT* lpMPC_grad_eq00 = lpMPC_grad_eq + 0;
lpMPC_FLOAT* lpMPC_grad_ineq00 = lpMPC_grad_ineq + 0;
lpMPC_FLOAT lpMPC_ctv00[4];
lpMPC_FLOAT* lpMPC_v00 = lpMPC_v + 0;
lpMPC_FLOAT lpMPC_re00[2];
lpMPC_FLOAT lpMPC_beta00[2];
lpMPC_FLOAT lpMPC_betacc00[2];
lpMPC_FLOAT* lpMPC_dvaff00 = lpMPC_dv_aff + 0;
lpMPC_FLOAT* lpMPC_dvcc00 = lpMPC_dv_cc + 0;
lpMPC_FLOAT lpMPC_V00[8];
lpMPC_FLOAT lpMPC_Yd00[3];
lpMPC_FLOAT lpMPC_Ld00[3];
lpMPC_FLOAT lpMPC_yy00[2];
lpMPC_FLOAT lpMPC_bmy00[2];
int lpMPC_lbIdx00[4] = {0, 1, 2, 3};
lpMPC_FLOAT* lpMPC_llb00 = lpMPC_l + 0;
lpMPC_FLOAT* lpMPC_slb00 = lpMPC_s + 0;
lpMPC_FLOAT* lpMPC_llbbyslb00 = lpMPC_lbys + 0;
lpMPC_FLOAT lpMPC_rilb00[4];
lpMPC_FLOAT* lpMPC_dllbaff00 = lpMPC_dl_aff + 0;
lpMPC_FLOAT* lpMPC_dslbaff00 = lpMPC_ds_aff + 0;
lpMPC_FLOAT* lpMPC_dllbcc00 = lpMPC_dl_cc + 0;
lpMPC_FLOAT* lpMPC_dslbcc00 = lpMPC_ds_cc + 0;
lpMPC_FLOAT* lpMPC_ccrhsl00 = lpMPC_ccrhs + 0;
int lpMPC_ubIdx00[4] = {0, 1, 2, 3};
lpMPC_FLOAT* lpMPC_lub00 = lpMPC_l + 4;
lpMPC_FLOAT* lpMPC_sub00 = lpMPC_s + 4;
lpMPC_FLOAT* lpMPC_lubbysub00 = lpMPC_lbys + 4;
lpMPC_FLOAT lpMPC_riub00[4];
lpMPC_FLOAT* lpMPC_dlubaff00 = lpMPC_dl_aff + 4;
lpMPC_FLOAT* lpMPC_dsubaff00 = lpMPC_ds_aff + 4;
lpMPC_FLOAT* lpMPC_dlubcc00 = lpMPC_dl_cc + 4;
lpMPC_FLOAT* lpMPC_dsubcc00 = lpMPC_ds_cc + 4;
lpMPC_FLOAT* lpMPC_ccrhsub00 = lpMPC_ccrhs + 4;
lpMPC_FLOAT lpMPC_Phi00[4];
lpMPC_FLOAT lpMPC_D00[4] = {1.0000000000000000E+000, 
1.0000000000000000E+000};
lpMPC_FLOAT lpMPC_W00[4];
lpMPC_FLOAT* lpMPC_z01 = lpMPC_z + 4;
lpMPC_FLOAT* lpMPC_dzaff01 = lpMPC_dz_aff + 4;
lpMPC_FLOAT* lpMPC_dzcc01 = lpMPC_dz_cc + 4;
lpMPC_FLOAT* lpMPC_rd01 = lpMPC_rd + 4;
lpMPC_FLOAT lpMPC_Lbyrd01[4];
lpMPC_FLOAT* lpMPC_grad_cost01 = lpMPC_grad_cost + 4;
lpMPC_FLOAT* lpMPC_grad_eq01 = lpMPC_grad_eq + 4;
lpMPC_FLOAT* lpMPC_grad_ineq01 = lpMPC_grad_ineq + 4;
lpMPC_FLOAT lpMPC_ctv01[4];
lpMPC_FLOAT* lpMPC_v01 = lpMPC_v + 2;
lpMPC_FLOAT lpMPC_re01[2];
lpMPC_FLOAT lpMPC_beta01[2];
lpMPC_FLOAT lpMPC_betacc01[2];
lpMPC_FLOAT* lpMPC_dvaff01 = lpMPC_dv_aff + 2;
lpMPC_FLOAT* lpMPC_dvcc01 = lpMPC_dv_cc + 2;
lpMPC_FLOAT lpMPC_V01[8];
lpMPC_FLOAT lpMPC_Yd01[3];
lpMPC_FLOAT lpMPC_Ld01[3];
lpMPC_FLOAT lpMPC_yy01[2];
lpMPC_FLOAT lpMPC_bmy01[2];
int lpMPC_lbIdx01[4] = {0, 1, 2, 3};
lpMPC_FLOAT* lpMPC_llb01 = lpMPC_l + 8;
lpMPC_FLOAT* lpMPC_slb01 = lpMPC_s + 8;
lpMPC_FLOAT* lpMPC_llbbyslb01 = lpMPC_lbys + 8;
lpMPC_FLOAT lpMPC_rilb01[4];
lpMPC_FLOAT* lpMPC_dllbaff01 = lpMPC_dl_aff + 8;
lpMPC_FLOAT* lpMPC_dslbaff01 = lpMPC_ds_aff + 8;
lpMPC_FLOAT* lpMPC_dllbcc01 = lpMPC_dl_cc + 8;
lpMPC_FLOAT* lpMPC_dslbcc01 = lpMPC_ds_cc + 8;
lpMPC_FLOAT* lpMPC_ccrhsl01 = lpMPC_ccrhs + 8;
int lpMPC_ubIdx01[4] = {0, 1, 2, 3};
lpMPC_FLOAT* lpMPC_lub01 = lpMPC_l + 12;
lpMPC_FLOAT* lpMPC_sub01 = lpMPC_s + 12;
lpMPC_FLOAT* lpMPC_lubbysub01 = lpMPC_lbys + 12;
lpMPC_FLOAT lpMPC_riub01[4];
lpMPC_FLOAT* lpMPC_dlubaff01 = lpMPC_dl_aff + 12;
lpMPC_FLOAT* lpMPC_dsubaff01 = lpMPC_ds_aff + 12;
lpMPC_FLOAT* lpMPC_dlubcc01 = lpMPC_dl_cc + 12;
lpMPC_FLOAT* lpMPC_dsubcc01 = lpMPC_ds_cc + 12;
lpMPC_FLOAT* lpMPC_ccrhsub01 = lpMPC_ccrhs + 12;
lpMPC_FLOAT lpMPC_Phi01[4];
lpMPC_FLOAT lpMPC_D01[4] = {-1.0000000000000000E+000, 
-1.0000000000000000E+000};
lpMPC_FLOAT lpMPC_W01[4];
lpMPC_FLOAT lpMPC_Ysd01[4];
lpMPC_FLOAT lpMPC_Lsd01[4];
lpMPC_FLOAT* lpMPC_z02 = lpMPC_z + 8;
lpMPC_FLOAT* lpMPC_dzaff02 = lpMPC_dz_aff + 8;
lpMPC_FLOAT* lpMPC_dzcc02 = lpMPC_dz_cc + 8;
lpMPC_FLOAT* lpMPC_rd02 = lpMPC_rd + 8;
lpMPC_FLOAT lpMPC_Lbyrd02[4];
lpMPC_FLOAT* lpMPC_grad_cost02 = lpMPC_grad_cost + 8;
lpMPC_FLOAT* lpMPC_grad_eq02 = lpMPC_grad_eq + 8;
lpMPC_FLOAT* lpMPC_grad_ineq02 = lpMPC_grad_ineq + 8;
lpMPC_FLOAT lpMPC_ctv02[4];
lpMPC_FLOAT* lpMPC_v02 = lpMPC_v + 4;
lpMPC_FLOAT lpMPC_re02[2];
lpMPC_FLOAT lpMPC_beta02[2];
lpMPC_FLOAT lpMPC_betacc02[2];
lpMPC_FLOAT* lpMPC_dvaff02 = lpMPC_dv_aff + 4;
lpMPC_FLOAT* lpMPC_dvcc02 = lpMPC_dv_cc + 4;
lpMPC_FLOAT lpMPC_V02[8];
lpMPC_FLOAT lpMPC_Yd02[3];
lpMPC_FLOAT lpMPC_Ld02[3];
lpMPC_FLOAT lpMPC_yy02[2];
lpMPC_FLOAT lpMPC_bmy02[2];
int lpMPC_lbIdx02[4] = {0, 1, 2, 3};
lpMPC_FLOAT* lpMPC_llb02 = lpMPC_l + 16;
lpMPC_FLOAT* lpMPC_slb02 = lpMPC_s + 16;
lpMPC_FLOAT* lpMPC_llbbyslb02 = lpMPC_lbys + 16;
lpMPC_FLOAT lpMPC_rilb02[4];
lpMPC_FLOAT* lpMPC_dllbaff02 = lpMPC_dl_aff + 16;
lpMPC_FLOAT* lpMPC_dslbaff02 = lpMPC_ds_aff + 16;
lpMPC_FLOAT* lpMPC_dllbcc02 = lpMPC_dl_cc + 16;
lpMPC_FLOAT* lpMPC_dslbcc02 = lpMPC_ds_cc + 16;
lpMPC_FLOAT* lpMPC_ccrhsl02 = lpMPC_ccrhs + 16;
int lpMPC_ubIdx02[4] = {0, 1, 2, 3};
lpMPC_FLOAT* lpMPC_lub02 = lpMPC_l + 20;
lpMPC_FLOAT* lpMPC_sub02 = lpMPC_s + 20;
lpMPC_FLOAT* lpMPC_lubbysub02 = lpMPC_lbys + 20;
lpMPC_FLOAT lpMPC_riub02[4];
lpMPC_FLOAT* lpMPC_dlubaff02 = lpMPC_dl_aff + 20;
lpMPC_FLOAT* lpMPC_dsubaff02 = lpMPC_ds_aff + 20;
lpMPC_FLOAT* lpMPC_dlubcc02 = lpMPC_dl_cc + 20;
lpMPC_FLOAT* lpMPC_dsubcc02 = lpMPC_ds_cc + 20;
lpMPC_FLOAT* lpMPC_ccrhsub02 = lpMPC_ccrhs + 20;
lpMPC_FLOAT lpMPC_Phi02[4];
lpMPC_FLOAT lpMPC_W02[4];
lpMPC_FLOAT lpMPC_Ysd02[4];
lpMPC_FLOAT lpMPC_Lsd02[4];
lpMPC_FLOAT* lpMPC_z03 = lpMPC_z + 12;
lpMPC_FLOAT* lpMPC_dzaff03 = lpMPC_dz_aff + 12;
lpMPC_FLOAT* lpMPC_dzcc03 = lpMPC_dz_cc + 12;
lpMPC_FLOAT* lpMPC_rd03 = lpMPC_rd + 12;
lpMPC_FLOAT lpMPC_Lbyrd03[4];
lpMPC_FLOAT* lpMPC_grad_cost03 = lpMPC_grad_cost + 12;
lpMPC_FLOAT* lpMPC_grad_eq03 = lpMPC_grad_eq + 12;
lpMPC_FLOAT* lpMPC_grad_ineq03 = lpMPC_grad_ineq + 12;
lpMPC_FLOAT lpMPC_ctv03[4];
lpMPC_FLOAT* lpMPC_v03 = lpMPC_v + 6;
lpMPC_FLOAT lpMPC_re03[2];
lpMPC_FLOAT lpMPC_beta03[2];
lpMPC_FLOAT lpMPC_betacc03[2];
lpMPC_FLOAT* lpMPC_dvaff03 = lpMPC_dv_aff + 6;
lpMPC_FLOAT* lpMPC_dvcc03 = lpMPC_dv_cc + 6;
lpMPC_FLOAT lpMPC_V03[8];
lpMPC_FLOAT lpMPC_Yd03[3];
lpMPC_FLOAT lpMPC_Ld03[3];
lpMPC_FLOAT lpMPC_yy03[2];
lpMPC_FLOAT lpMPC_bmy03[2];
int lpMPC_lbIdx03[4] = {0, 1, 2, 3};
lpMPC_FLOAT* lpMPC_llb03 = lpMPC_l + 24;
lpMPC_FLOAT* lpMPC_slb03 = lpMPC_s + 24;
lpMPC_FLOAT* lpMPC_llbbyslb03 = lpMPC_lbys + 24;
lpMPC_FLOAT lpMPC_rilb03[4];
lpMPC_FLOAT* lpMPC_dllbaff03 = lpMPC_dl_aff + 24;
lpMPC_FLOAT* lpMPC_dslbaff03 = lpMPC_ds_aff + 24;
lpMPC_FLOAT* lpMPC_dllbcc03 = lpMPC_dl_cc + 24;
lpMPC_FLOAT* lpMPC_dslbcc03 = lpMPC_ds_cc + 24;
lpMPC_FLOAT* lpMPC_ccrhsl03 = lpMPC_ccrhs + 24;
int lpMPC_ubIdx03[4] = {0, 1, 2, 3};
lpMPC_FLOAT* lpMPC_lub03 = lpMPC_l + 28;
lpMPC_FLOAT* lpMPC_sub03 = lpMPC_s + 28;
lpMPC_FLOAT* lpMPC_lubbysub03 = lpMPC_lbys + 28;
lpMPC_FLOAT lpMPC_riub03[4];
lpMPC_FLOAT* lpMPC_dlubaff03 = lpMPC_dl_aff + 28;
lpMPC_FLOAT* lpMPC_dsubaff03 = lpMPC_ds_aff + 28;
lpMPC_FLOAT* lpMPC_dlubcc03 = lpMPC_dl_cc + 28;
lpMPC_FLOAT* lpMPC_dsubcc03 = lpMPC_ds_cc + 28;
lpMPC_FLOAT* lpMPC_ccrhsub03 = lpMPC_ccrhs + 28;
lpMPC_FLOAT lpMPC_Phi03[4];
lpMPC_FLOAT lpMPC_W03[4];
lpMPC_FLOAT lpMPC_Ysd03[4];
lpMPC_FLOAT lpMPC_Lsd03[4];
lpMPC_FLOAT* lpMPC_z04 = lpMPC_z + 16;
lpMPC_FLOAT* lpMPC_dzaff04 = lpMPC_dz_aff + 16;
lpMPC_FLOAT* lpMPC_dzcc04 = lpMPC_dz_cc + 16;
lpMPC_FLOAT* lpMPC_rd04 = lpMPC_rd + 16;
lpMPC_FLOAT lpMPC_Lbyrd04[4];
lpMPC_FLOAT* lpMPC_grad_cost04 = lpMPC_grad_cost + 16;
lpMPC_FLOAT* lpMPC_grad_eq04 = lpMPC_grad_eq + 16;
lpMPC_FLOAT* lpMPC_grad_ineq04 = lpMPC_grad_ineq + 16;
lpMPC_FLOAT lpMPC_ctv04[4];
lpMPC_FLOAT* lpMPC_v04 = lpMPC_v + 8;
lpMPC_FLOAT lpMPC_re04[2];
lpMPC_FLOAT lpMPC_beta04[2];
lpMPC_FLOAT lpMPC_betacc04[2];
lpMPC_FLOAT* lpMPC_dvaff04 = lpMPC_dv_aff + 8;
lpMPC_FLOAT* lpMPC_dvcc04 = lpMPC_dv_cc + 8;
lpMPC_FLOAT lpMPC_V04[8];
lpMPC_FLOAT lpMPC_Yd04[3];
lpMPC_FLOAT lpMPC_Ld04[3];
lpMPC_FLOAT lpMPC_yy04[2];
lpMPC_FLOAT lpMPC_bmy04[2];
int lpMPC_lbIdx04[4] = {0, 1, 2, 3};
lpMPC_FLOAT* lpMPC_llb04 = lpMPC_l + 32;
lpMPC_FLOAT* lpMPC_slb04 = lpMPC_s + 32;
lpMPC_FLOAT* lpMPC_llbbyslb04 = lpMPC_lbys + 32;
lpMPC_FLOAT lpMPC_rilb04[4];
lpMPC_FLOAT* lpMPC_dllbaff04 = lpMPC_dl_aff + 32;
lpMPC_FLOAT* lpMPC_dslbaff04 = lpMPC_ds_aff + 32;
lpMPC_FLOAT* lpMPC_dllbcc04 = lpMPC_dl_cc + 32;
lpMPC_FLOAT* lpMPC_dslbcc04 = lpMPC_ds_cc + 32;
lpMPC_FLOAT* lpMPC_ccrhsl04 = lpMPC_ccrhs + 32;
int lpMPC_ubIdx04[4] = {0, 1, 2, 3};
lpMPC_FLOAT* lpMPC_lub04 = lpMPC_l + 36;
lpMPC_FLOAT* lpMPC_sub04 = lpMPC_s + 36;
lpMPC_FLOAT* lpMPC_lubbysub04 = lpMPC_lbys + 36;
lpMPC_FLOAT lpMPC_riub04[4];
lpMPC_FLOAT* lpMPC_dlubaff04 = lpMPC_dl_aff + 36;
lpMPC_FLOAT* lpMPC_dsubaff04 = lpMPC_ds_aff + 36;
lpMPC_FLOAT* lpMPC_dlubcc04 = lpMPC_dl_cc + 36;
lpMPC_FLOAT* lpMPC_dsubcc04 = lpMPC_ds_cc + 36;
lpMPC_FLOAT* lpMPC_ccrhsub04 = lpMPC_ccrhs + 36;
lpMPC_FLOAT lpMPC_Phi04[4];
lpMPC_FLOAT lpMPC_W04[4];
lpMPC_FLOAT lpMPC_Ysd04[4];
lpMPC_FLOAT lpMPC_Lsd04[4];
lpMPC_FLOAT* lpMPC_z05 = lpMPC_z + 20;
lpMPC_FLOAT* lpMPC_dzaff05 = lpMPC_dz_aff + 20;
lpMPC_FLOAT* lpMPC_dzcc05 = lpMPC_dz_cc + 20;
lpMPC_FLOAT* lpMPC_rd05 = lpMPC_rd + 20;
lpMPC_FLOAT lpMPC_Lbyrd05[4];
lpMPC_FLOAT* lpMPC_grad_cost05 = lpMPC_grad_cost + 20;
lpMPC_FLOAT* lpMPC_grad_eq05 = lpMPC_grad_eq + 20;
lpMPC_FLOAT* lpMPC_grad_ineq05 = lpMPC_grad_ineq + 20;
lpMPC_FLOAT lpMPC_ctv05[4];
lpMPC_FLOAT* lpMPC_v05 = lpMPC_v + 10;
lpMPC_FLOAT lpMPC_re05[2];
lpMPC_FLOAT lpMPC_beta05[2];
lpMPC_FLOAT lpMPC_betacc05[2];
lpMPC_FLOAT* lpMPC_dvaff05 = lpMPC_dv_aff + 10;
lpMPC_FLOAT* lpMPC_dvcc05 = lpMPC_dv_cc + 10;
lpMPC_FLOAT lpMPC_V05[8];
lpMPC_FLOAT lpMPC_Yd05[3];
lpMPC_FLOAT lpMPC_Ld05[3];
lpMPC_FLOAT lpMPC_yy05[2];
lpMPC_FLOAT lpMPC_bmy05[2];
int lpMPC_lbIdx05[4] = {0, 1, 2, 3};
lpMPC_FLOAT* lpMPC_llb05 = lpMPC_l + 40;
lpMPC_FLOAT* lpMPC_slb05 = lpMPC_s + 40;
lpMPC_FLOAT* lpMPC_llbbyslb05 = lpMPC_lbys + 40;
lpMPC_FLOAT lpMPC_rilb05[4];
lpMPC_FLOAT* lpMPC_dllbaff05 = lpMPC_dl_aff + 40;
lpMPC_FLOAT* lpMPC_dslbaff05 = lpMPC_ds_aff + 40;
lpMPC_FLOAT* lpMPC_dllbcc05 = lpMPC_dl_cc + 40;
lpMPC_FLOAT* lpMPC_dslbcc05 = lpMPC_ds_cc + 40;
lpMPC_FLOAT* lpMPC_ccrhsl05 = lpMPC_ccrhs + 40;
int lpMPC_ubIdx05[4] = {0, 1, 2, 3};
lpMPC_FLOAT* lpMPC_lub05 = lpMPC_l + 44;
lpMPC_FLOAT* lpMPC_sub05 = lpMPC_s + 44;
lpMPC_FLOAT* lpMPC_lubbysub05 = lpMPC_lbys + 44;
lpMPC_FLOAT lpMPC_riub05[4];
lpMPC_FLOAT* lpMPC_dlubaff05 = lpMPC_dl_aff + 44;
lpMPC_FLOAT* lpMPC_dsubaff05 = lpMPC_ds_aff + 44;
lpMPC_FLOAT* lpMPC_dlubcc05 = lpMPC_dl_cc + 44;
lpMPC_FLOAT* lpMPC_dsubcc05 = lpMPC_ds_cc + 44;
lpMPC_FLOAT* lpMPC_ccrhsub05 = lpMPC_ccrhs + 44;
lpMPC_FLOAT lpMPC_Phi05[4];
lpMPC_FLOAT lpMPC_W05[4];
lpMPC_FLOAT lpMPC_Ysd05[4];
lpMPC_FLOAT lpMPC_Lsd05[4];
lpMPC_FLOAT* lpMPC_z06 = lpMPC_z + 24;
lpMPC_FLOAT* lpMPC_dzaff06 = lpMPC_dz_aff + 24;
lpMPC_FLOAT* lpMPC_dzcc06 = lpMPC_dz_cc + 24;
lpMPC_FLOAT* lpMPC_rd06 = lpMPC_rd + 24;
lpMPC_FLOAT lpMPC_Lbyrd06[4];
lpMPC_FLOAT* lpMPC_grad_cost06 = lpMPC_grad_cost + 24;
lpMPC_FLOAT* lpMPC_grad_eq06 = lpMPC_grad_eq + 24;
lpMPC_FLOAT* lpMPC_grad_ineq06 = lpMPC_grad_ineq + 24;
lpMPC_FLOAT lpMPC_ctv06[4];
lpMPC_FLOAT* lpMPC_v06 = lpMPC_v + 12;
lpMPC_FLOAT lpMPC_re06[2];
lpMPC_FLOAT lpMPC_beta06[2];
lpMPC_FLOAT lpMPC_betacc06[2];
lpMPC_FLOAT* lpMPC_dvaff06 = lpMPC_dv_aff + 12;
lpMPC_FLOAT* lpMPC_dvcc06 = lpMPC_dv_cc + 12;
lpMPC_FLOAT lpMPC_V06[8];
lpMPC_FLOAT lpMPC_Yd06[3];
lpMPC_FLOAT lpMPC_Ld06[3];
lpMPC_FLOAT lpMPC_yy06[2];
lpMPC_FLOAT lpMPC_bmy06[2];
int lpMPC_lbIdx06[4] = {0, 1, 2, 3};
lpMPC_FLOAT* lpMPC_llb06 = lpMPC_l + 48;
lpMPC_FLOAT* lpMPC_slb06 = lpMPC_s + 48;
lpMPC_FLOAT* lpMPC_llbbyslb06 = lpMPC_lbys + 48;
lpMPC_FLOAT lpMPC_rilb06[4];
lpMPC_FLOAT* lpMPC_dllbaff06 = lpMPC_dl_aff + 48;
lpMPC_FLOAT* lpMPC_dslbaff06 = lpMPC_ds_aff + 48;
lpMPC_FLOAT* lpMPC_dllbcc06 = lpMPC_dl_cc + 48;
lpMPC_FLOAT* lpMPC_dslbcc06 = lpMPC_ds_cc + 48;
lpMPC_FLOAT* lpMPC_ccrhsl06 = lpMPC_ccrhs + 48;
int lpMPC_ubIdx06[4] = {0, 1, 2, 3};
lpMPC_FLOAT* lpMPC_lub06 = lpMPC_l + 52;
lpMPC_FLOAT* lpMPC_sub06 = lpMPC_s + 52;
lpMPC_FLOAT* lpMPC_lubbysub06 = lpMPC_lbys + 52;
lpMPC_FLOAT lpMPC_riub06[4];
lpMPC_FLOAT* lpMPC_dlubaff06 = lpMPC_dl_aff + 52;
lpMPC_FLOAT* lpMPC_dsubaff06 = lpMPC_ds_aff + 52;
lpMPC_FLOAT* lpMPC_dlubcc06 = lpMPC_dl_cc + 52;
lpMPC_FLOAT* lpMPC_dsubcc06 = lpMPC_ds_cc + 52;
lpMPC_FLOAT* lpMPC_ccrhsub06 = lpMPC_ccrhs + 52;
lpMPC_FLOAT lpMPC_Phi06[4];
lpMPC_FLOAT lpMPC_W06[4];
lpMPC_FLOAT lpMPC_Ysd06[4];
lpMPC_FLOAT lpMPC_Lsd06[4];
lpMPC_FLOAT* lpMPC_z07 = lpMPC_z + 28;
lpMPC_FLOAT* lpMPC_dzaff07 = lpMPC_dz_aff + 28;
lpMPC_FLOAT* lpMPC_dzcc07 = lpMPC_dz_cc + 28;
lpMPC_FLOAT* lpMPC_rd07 = lpMPC_rd + 28;
lpMPC_FLOAT lpMPC_Lbyrd07[4];
lpMPC_FLOAT* lpMPC_grad_cost07 = lpMPC_grad_cost + 28;
lpMPC_FLOAT* lpMPC_grad_eq07 = lpMPC_grad_eq + 28;
lpMPC_FLOAT* lpMPC_grad_ineq07 = lpMPC_grad_ineq + 28;
lpMPC_FLOAT lpMPC_ctv07[4];
lpMPC_FLOAT* lpMPC_v07 = lpMPC_v + 14;
lpMPC_FLOAT lpMPC_re07[2];
lpMPC_FLOAT lpMPC_beta07[2];
lpMPC_FLOAT lpMPC_betacc07[2];
lpMPC_FLOAT* lpMPC_dvaff07 = lpMPC_dv_aff + 14;
lpMPC_FLOAT* lpMPC_dvcc07 = lpMPC_dv_cc + 14;
lpMPC_FLOAT lpMPC_V07[8];
lpMPC_FLOAT lpMPC_Yd07[3];
lpMPC_FLOAT lpMPC_Ld07[3];
lpMPC_FLOAT lpMPC_yy07[2];
lpMPC_FLOAT lpMPC_bmy07[2];
int lpMPC_lbIdx07[4] = {0, 1, 2, 3};
lpMPC_FLOAT* lpMPC_llb07 = lpMPC_l + 56;
lpMPC_FLOAT* lpMPC_slb07 = lpMPC_s + 56;
lpMPC_FLOAT* lpMPC_llbbyslb07 = lpMPC_lbys + 56;
lpMPC_FLOAT lpMPC_rilb07[4];
lpMPC_FLOAT* lpMPC_dllbaff07 = lpMPC_dl_aff + 56;
lpMPC_FLOAT* lpMPC_dslbaff07 = lpMPC_ds_aff + 56;
lpMPC_FLOAT* lpMPC_dllbcc07 = lpMPC_dl_cc + 56;
lpMPC_FLOAT* lpMPC_dslbcc07 = lpMPC_ds_cc + 56;
lpMPC_FLOAT* lpMPC_ccrhsl07 = lpMPC_ccrhs + 56;
int lpMPC_ubIdx07[4] = {0, 1, 2, 3};
lpMPC_FLOAT* lpMPC_lub07 = lpMPC_l + 60;
lpMPC_FLOAT* lpMPC_sub07 = lpMPC_s + 60;
lpMPC_FLOAT* lpMPC_lubbysub07 = lpMPC_lbys + 60;
lpMPC_FLOAT lpMPC_riub07[4];
lpMPC_FLOAT* lpMPC_dlubaff07 = lpMPC_dl_aff + 60;
lpMPC_FLOAT* lpMPC_dsubaff07 = lpMPC_ds_aff + 60;
lpMPC_FLOAT* lpMPC_dlubcc07 = lpMPC_dl_cc + 60;
lpMPC_FLOAT* lpMPC_dsubcc07 = lpMPC_ds_cc + 60;
lpMPC_FLOAT* lpMPC_ccrhsub07 = lpMPC_ccrhs + 60;
lpMPC_FLOAT lpMPC_Phi07[4];
lpMPC_FLOAT lpMPC_W07[4];
lpMPC_FLOAT lpMPC_Ysd07[4];
lpMPC_FLOAT lpMPC_Lsd07[4];
lpMPC_FLOAT* lpMPC_z08 = lpMPC_z + 32;
lpMPC_FLOAT* lpMPC_dzaff08 = lpMPC_dz_aff + 32;
lpMPC_FLOAT* lpMPC_dzcc08 = lpMPC_dz_cc + 32;
lpMPC_FLOAT* lpMPC_rd08 = lpMPC_rd + 32;
lpMPC_FLOAT lpMPC_Lbyrd08[4];
lpMPC_FLOAT* lpMPC_grad_cost08 = lpMPC_grad_cost + 32;
lpMPC_FLOAT* lpMPC_grad_eq08 = lpMPC_grad_eq + 32;
lpMPC_FLOAT* lpMPC_grad_ineq08 = lpMPC_grad_ineq + 32;
lpMPC_FLOAT lpMPC_ctv08[4];
lpMPC_FLOAT* lpMPC_v08 = lpMPC_v + 16;
lpMPC_FLOAT lpMPC_re08[2];
lpMPC_FLOAT lpMPC_beta08[2];
lpMPC_FLOAT lpMPC_betacc08[2];
lpMPC_FLOAT* lpMPC_dvaff08 = lpMPC_dv_aff + 16;
lpMPC_FLOAT* lpMPC_dvcc08 = lpMPC_dv_cc + 16;
lpMPC_FLOAT lpMPC_V08[8];
lpMPC_FLOAT lpMPC_Yd08[3];
lpMPC_FLOAT lpMPC_Ld08[3];
lpMPC_FLOAT lpMPC_yy08[2];
lpMPC_FLOAT lpMPC_bmy08[2];
int lpMPC_lbIdx08[4] = {0, 1, 2, 3};
lpMPC_FLOAT* lpMPC_llb08 = lpMPC_l + 64;
lpMPC_FLOAT* lpMPC_slb08 = lpMPC_s + 64;
lpMPC_FLOAT* lpMPC_llbbyslb08 = lpMPC_lbys + 64;
lpMPC_FLOAT lpMPC_rilb08[4];
lpMPC_FLOAT* lpMPC_dllbaff08 = lpMPC_dl_aff + 64;
lpMPC_FLOAT* lpMPC_dslbaff08 = lpMPC_ds_aff + 64;
lpMPC_FLOAT* lpMPC_dllbcc08 = lpMPC_dl_cc + 64;
lpMPC_FLOAT* lpMPC_dslbcc08 = lpMPC_ds_cc + 64;
lpMPC_FLOAT* lpMPC_ccrhsl08 = lpMPC_ccrhs + 64;
int lpMPC_ubIdx08[4] = {0, 1, 2, 3};
lpMPC_FLOAT* lpMPC_lub08 = lpMPC_l + 68;
lpMPC_FLOAT* lpMPC_sub08 = lpMPC_s + 68;
lpMPC_FLOAT* lpMPC_lubbysub08 = lpMPC_lbys + 68;
lpMPC_FLOAT lpMPC_riub08[4];
lpMPC_FLOAT* lpMPC_dlubaff08 = lpMPC_dl_aff + 68;
lpMPC_FLOAT* lpMPC_dsubaff08 = lpMPC_ds_aff + 68;
lpMPC_FLOAT* lpMPC_dlubcc08 = lpMPC_dl_cc + 68;
lpMPC_FLOAT* lpMPC_dsubcc08 = lpMPC_ds_cc + 68;
lpMPC_FLOAT* lpMPC_ccrhsub08 = lpMPC_ccrhs + 68;
lpMPC_FLOAT lpMPC_Phi08[4];
lpMPC_FLOAT lpMPC_W08[4];
lpMPC_FLOAT lpMPC_Ysd08[4];
lpMPC_FLOAT lpMPC_Lsd08[4];
lpMPC_FLOAT* lpMPC_z09 = lpMPC_z + 36;
lpMPC_FLOAT* lpMPC_dzaff09 = lpMPC_dz_aff + 36;
lpMPC_FLOAT* lpMPC_dzcc09 = lpMPC_dz_cc + 36;
lpMPC_FLOAT* lpMPC_rd09 = lpMPC_rd + 36;
lpMPC_FLOAT lpMPC_Lbyrd09[4];
lpMPC_FLOAT* lpMPC_grad_cost09 = lpMPC_grad_cost + 36;
lpMPC_FLOAT* lpMPC_grad_eq09 = lpMPC_grad_eq + 36;
lpMPC_FLOAT* lpMPC_grad_ineq09 = lpMPC_grad_ineq + 36;
lpMPC_FLOAT lpMPC_ctv09[4];
lpMPC_FLOAT* lpMPC_v09 = lpMPC_v + 18;
lpMPC_FLOAT lpMPC_re09[2];
lpMPC_FLOAT lpMPC_beta09[2];
lpMPC_FLOAT lpMPC_betacc09[2];
lpMPC_FLOAT* lpMPC_dvaff09 = lpMPC_dv_aff + 18;
lpMPC_FLOAT* lpMPC_dvcc09 = lpMPC_dv_cc + 18;
lpMPC_FLOAT lpMPC_V09[8];
lpMPC_FLOAT lpMPC_Yd09[3];
lpMPC_FLOAT lpMPC_Ld09[3];
lpMPC_FLOAT lpMPC_yy09[2];
lpMPC_FLOAT lpMPC_bmy09[2];
int lpMPC_lbIdx09[4] = {0, 1, 2, 3};
lpMPC_FLOAT* lpMPC_llb09 = lpMPC_l + 72;
lpMPC_FLOAT* lpMPC_slb09 = lpMPC_s + 72;
lpMPC_FLOAT* lpMPC_llbbyslb09 = lpMPC_lbys + 72;
lpMPC_FLOAT lpMPC_rilb09[4];
lpMPC_FLOAT* lpMPC_dllbaff09 = lpMPC_dl_aff + 72;
lpMPC_FLOAT* lpMPC_dslbaff09 = lpMPC_ds_aff + 72;
lpMPC_FLOAT* lpMPC_dllbcc09 = lpMPC_dl_cc + 72;
lpMPC_FLOAT* lpMPC_dslbcc09 = lpMPC_ds_cc + 72;
lpMPC_FLOAT* lpMPC_ccrhsl09 = lpMPC_ccrhs + 72;
int lpMPC_ubIdx09[4] = {0, 1, 2, 3};
lpMPC_FLOAT* lpMPC_lub09 = lpMPC_l + 76;
lpMPC_FLOAT* lpMPC_sub09 = lpMPC_s + 76;
lpMPC_FLOAT* lpMPC_lubbysub09 = lpMPC_lbys + 76;
lpMPC_FLOAT lpMPC_riub09[4];
lpMPC_FLOAT* lpMPC_dlubaff09 = lpMPC_dl_aff + 76;
lpMPC_FLOAT* lpMPC_dsubaff09 = lpMPC_ds_aff + 76;
lpMPC_FLOAT* lpMPC_dlubcc09 = lpMPC_dl_cc + 76;
lpMPC_FLOAT* lpMPC_dsubcc09 = lpMPC_ds_cc + 76;
lpMPC_FLOAT* lpMPC_ccrhsub09 = lpMPC_ccrhs + 76;
lpMPC_FLOAT lpMPC_Phi09[4];
lpMPC_FLOAT lpMPC_W09[4];
lpMPC_FLOAT lpMPC_Ysd09[4];
lpMPC_FLOAT lpMPC_Lsd09[4];
lpMPC_FLOAT* lpMPC_z10 = lpMPC_z + 40;
lpMPC_FLOAT* lpMPC_dzaff10 = lpMPC_dz_aff + 40;
lpMPC_FLOAT* lpMPC_dzcc10 = lpMPC_dz_cc + 40;
lpMPC_FLOAT* lpMPC_rd10 = lpMPC_rd + 40;
lpMPC_FLOAT lpMPC_Lbyrd10[4];
lpMPC_FLOAT* lpMPC_grad_cost10 = lpMPC_grad_cost + 40;
lpMPC_FLOAT* lpMPC_grad_eq10 = lpMPC_grad_eq + 40;
lpMPC_FLOAT* lpMPC_grad_ineq10 = lpMPC_grad_ineq + 40;
lpMPC_FLOAT lpMPC_ctv10[4];
lpMPC_FLOAT* lpMPC_v10 = lpMPC_v + 20;
lpMPC_FLOAT lpMPC_re10[2];
lpMPC_FLOAT lpMPC_beta10[2];
lpMPC_FLOAT lpMPC_betacc10[2];
lpMPC_FLOAT* lpMPC_dvaff10 = lpMPC_dv_aff + 20;
lpMPC_FLOAT* lpMPC_dvcc10 = lpMPC_dv_cc + 20;
lpMPC_FLOAT lpMPC_V10[8];
lpMPC_FLOAT lpMPC_Yd10[3];
lpMPC_FLOAT lpMPC_Ld10[3];
lpMPC_FLOAT lpMPC_yy10[2];
lpMPC_FLOAT lpMPC_bmy10[2];
int lpMPC_lbIdx10[4] = {0, 1, 2, 3};
lpMPC_FLOAT* lpMPC_llb10 = lpMPC_l + 80;
lpMPC_FLOAT* lpMPC_slb10 = lpMPC_s + 80;
lpMPC_FLOAT* lpMPC_llbbyslb10 = lpMPC_lbys + 80;
lpMPC_FLOAT lpMPC_rilb10[4];
lpMPC_FLOAT* lpMPC_dllbaff10 = lpMPC_dl_aff + 80;
lpMPC_FLOAT* lpMPC_dslbaff10 = lpMPC_ds_aff + 80;
lpMPC_FLOAT* lpMPC_dllbcc10 = lpMPC_dl_cc + 80;
lpMPC_FLOAT* lpMPC_dslbcc10 = lpMPC_ds_cc + 80;
lpMPC_FLOAT* lpMPC_ccrhsl10 = lpMPC_ccrhs + 80;
int lpMPC_ubIdx10[4] = {0, 1, 2, 3};
lpMPC_FLOAT* lpMPC_lub10 = lpMPC_l + 84;
lpMPC_FLOAT* lpMPC_sub10 = lpMPC_s + 84;
lpMPC_FLOAT* lpMPC_lubbysub10 = lpMPC_lbys + 84;
lpMPC_FLOAT lpMPC_riub10[4];
lpMPC_FLOAT* lpMPC_dlubaff10 = lpMPC_dl_aff + 84;
lpMPC_FLOAT* lpMPC_dsubaff10 = lpMPC_ds_aff + 84;
lpMPC_FLOAT* lpMPC_dlubcc10 = lpMPC_dl_cc + 84;
lpMPC_FLOAT* lpMPC_dsubcc10 = lpMPC_ds_cc + 84;
lpMPC_FLOAT* lpMPC_ccrhsub10 = lpMPC_ccrhs + 84;
lpMPC_FLOAT lpMPC_Phi10[4];
lpMPC_FLOAT lpMPC_W10[4];
lpMPC_FLOAT lpMPC_Ysd10[4];
lpMPC_FLOAT lpMPC_Lsd10[4];
lpMPC_FLOAT* lpMPC_z11 = lpMPC_z + 44;
lpMPC_FLOAT* lpMPC_dzaff11 = lpMPC_dz_aff + 44;
lpMPC_FLOAT* lpMPC_dzcc11 = lpMPC_dz_cc + 44;
lpMPC_FLOAT* lpMPC_rd11 = lpMPC_rd + 44;
lpMPC_FLOAT lpMPC_Lbyrd11[4];
lpMPC_FLOAT* lpMPC_grad_cost11 = lpMPC_grad_cost + 44;
lpMPC_FLOAT* lpMPC_grad_eq11 = lpMPC_grad_eq + 44;
lpMPC_FLOAT* lpMPC_grad_ineq11 = lpMPC_grad_ineq + 44;
lpMPC_FLOAT lpMPC_ctv11[4];
lpMPC_FLOAT* lpMPC_v11 = lpMPC_v + 22;
lpMPC_FLOAT lpMPC_re11[2];
lpMPC_FLOAT lpMPC_beta11[2];
lpMPC_FLOAT lpMPC_betacc11[2];
lpMPC_FLOAT* lpMPC_dvaff11 = lpMPC_dv_aff + 22;
lpMPC_FLOAT* lpMPC_dvcc11 = lpMPC_dv_cc + 22;
lpMPC_FLOAT lpMPC_V11[8];
lpMPC_FLOAT lpMPC_Yd11[3];
lpMPC_FLOAT lpMPC_Ld11[3];
lpMPC_FLOAT lpMPC_yy11[2];
lpMPC_FLOAT lpMPC_bmy11[2];
int lpMPC_lbIdx11[4] = {0, 1, 2, 3};
lpMPC_FLOAT* lpMPC_llb11 = lpMPC_l + 88;
lpMPC_FLOAT* lpMPC_slb11 = lpMPC_s + 88;
lpMPC_FLOAT* lpMPC_llbbyslb11 = lpMPC_lbys + 88;
lpMPC_FLOAT lpMPC_rilb11[4];
lpMPC_FLOAT* lpMPC_dllbaff11 = lpMPC_dl_aff + 88;
lpMPC_FLOAT* lpMPC_dslbaff11 = lpMPC_ds_aff + 88;
lpMPC_FLOAT* lpMPC_dllbcc11 = lpMPC_dl_cc + 88;
lpMPC_FLOAT* lpMPC_dslbcc11 = lpMPC_ds_cc + 88;
lpMPC_FLOAT* lpMPC_ccrhsl11 = lpMPC_ccrhs + 88;
int lpMPC_ubIdx11[4] = {0, 1, 2, 3};
lpMPC_FLOAT* lpMPC_lub11 = lpMPC_l + 92;
lpMPC_FLOAT* lpMPC_sub11 = lpMPC_s + 92;
lpMPC_FLOAT* lpMPC_lubbysub11 = lpMPC_lbys + 92;
lpMPC_FLOAT lpMPC_riub11[4];
lpMPC_FLOAT* lpMPC_dlubaff11 = lpMPC_dl_aff + 92;
lpMPC_FLOAT* lpMPC_dsubaff11 = lpMPC_ds_aff + 92;
lpMPC_FLOAT* lpMPC_dlubcc11 = lpMPC_dl_cc + 92;
lpMPC_FLOAT* lpMPC_dsubcc11 = lpMPC_ds_cc + 92;
lpMPC_FLOAT* lpMPC_ccrhsub11 = lpMPC_ccrhs + 92;
lpMPC_FLOAT lpMPC_Phi11[4];
lpMPC_FLOAT lpMPC_W11[4];
lpMPC_FLOAT lpMPC_Ysd11[4];
lpMPC_FLOAT lpMPC_Lsd11[4];
lpMPC_FLOAT* lpMPC_z12 = lpMPC_z + 48;
lpMPC_FLOAT* lpMPC_dzaff12 = lpMPC_dz_aff + 48;
lpMPC_FLOAT* lpMPC_dzcc12 = lpMPC_dz_cc + 48;
lpMPC_FLOAT* lpMPC_rd12 = lpMPC_rd + 48;
lpMPC_FLOAT lpMPC_Lbyrd12[4];
lpMPC_FLOAT* lpMPC_grad_cost12 = lpMPC_grad_cost + 48;
lpMPC_FLOAT* lpMPC_grad_eq12 = lpMPC_grad_eq + 48;
lpMPC_FLOAT* lpMPC_grad_ineq12 = lpMPC_grad_ineq + 48;
lpMPC_FLOAT lpMPC_ctv12[4];
lpMPC_FLOAT* lpMPC_v12 = lpMPC_v + 24;
lpMPC_FLOAT lpMPC_re12[2];
lpMPC_FLOAT lpMPC_beta12[2];
lpMPC_FLOAT lpMPC_betacc12[2];
lpMPC_FLOAT* lpMPC_dvaff12 = lpMPC_dv_aff + 24;
lpMPC_FLOAT* lpMPC_dvcc12 = lpMPC_dv_cc + 24;
lpMPC_FLOAT lpMPC_V12[8];
lpMPC_FLOAT lpMPC_Yd12[3];
lpMPC_FLOAT lpMPC_Ld12[3];
lpMPC_FLOAT lpMPC_yy12[2];
lpMPC_FLOAT lpMPC_bmy12[2];
int lpMPC_lbIdx12[4] = {0, 1, 2, 3};
lpMPC_FLOAT* lpMPC_llb12 = lpMPC_l + 96;
lpMPC_FLOAT* lpMPC_slb12 = lpMPC_s + 96;
lpMPC_FLOAT* lpMPC_llbbyslb12 = lpMPC_lbys + 96;
lpMPC_FLOAT lpMPC_rilb12[4];
lpMPC_FLOAT* lpMPC_dllbaff12 = lpMPC_dl_aff + 96;
lpMPC_FLOAT* lpMPC_dslbaff12 = lpMPC_ds_aff + 96;
lpMPC_FLOAT* lpMPC_dllbcc12 = lpMPC_dl_cc + 96;
lpMPC_FLOAT* lpMPC_dslbcc12 = lpMPC_ds_cc + 96;
lpMPC_FLOAT* lpMPC_ccrhsl12 = lpMPC_ccrhs + 96;
int lpMPC_ubIdx12[4] = {0, 1, 2, 3};
lpMPC_FLOAT* lpMPC_lub12 = lpMPC_l + 100;
lpMPC_FLOAT* lpMPC_sub12 = lpMPC_s + 100;
lpMPC_FLOAT* lpMPC_lubbysub12 = lpMPC_lbys + 100;
lpMPC_FLOAT lpMPC_riub12[4];
lpMPC_FLOAT* lpMPC_dlubaff12 = lpMPC_dl_aff + 100;
lpMPC_FLOAT* lpMPC_dsubaff12 = lpMPC_ds_aff + 100;
lpMPC_FLOAT* lpMPC_dlubcc12 = lpMPC_dl_cc + 100;
lpMPC_FLOAT* lpMPC_dsubcc12 = lpMPC_ds_cc + 100;
lpMPC_FLOAT* lpMPC_ccrhsub12 = lpMPC_ccrhs + 100;
lpMPC_FLOAT lpMPC_Phi12[4];
lpMPC_FLOAT lpMPC_W12[4];
lpMPC_FLOAT lpMPC_Ysd12[4];
lpMPC_FLOAT lpMPC_Lsd12[4];
lpMPC_FLOAT* lpMPC_z13 = lpMPC_z + 52;
lpMPC_FLOAT* lpMPC_dzaff13 = lpMPC_dz_aff + 52;
lpMPC_FLOAT* lpMPC_dzcc13 = lpMPC_dz_cc + 52;
lpMPC_FLOAT* lpMPC_rd13 = lpMPC_rd + 52;
lpMPC_FLOAT lpMPC_Lbyrd13[4];
lpMPC_FLOAT* lpMPC_grad_cost13 = lpMPC_grad_cost + 52;
lpMPC_FLOAT* lpMPC_grad_eq13 = lpMPC_grad_eq + 52;
lpMPC_FLOAT* lpMPC_grad_ineq13 = lpMPC_grad_ineq + 52;
lpMPC_FLOAT lpMPC_ctv13[4];
lpMPC_FLOAT* lpMPC_v13 = lpMPC_v + 26;
lpMPC_FLOAT lpMPC_re13[2];
lpMPC_FLOAT lpMPC_beta13[2];
lpMPC_FLOAT lpMPC_betacc13[2];
lpMPC_FLOAT* lpMPC_dvaff13 = lpMPC_dv_aff + 26;
lpMPC_FLOAT* lpMPC_dvcc13 = lpMPC_dv_cc + 26;
lpMPC_FLOAT lpMPC_V13[8];
lpMPC_FLOAT lpMPC_Yd13[3];
lpMPC_FLOAT lpMPC_Ld13[3];
lpMPC_FLOAT lpMPC_yy13[2];
lpMPC_FLOAT lpMPC_bmy13[2];
int lpMPC_lbIdx13[4] = {0, 1, 2, 3};
lpMPC_FLOAT* lpMPC_llb13 = lpMPC_l + 104;
lpMPC_FLOAT* lpMPC_slb13 = lpMPC_s + 104;
lpMPC_FLOAT* lpMPC_llbbyslb13 = lpMPC_lbys + 104;
lpMPC_FLOAT lpMPC_rilb13[4];
lpMPC_FLOAT* lpMPC_dllbaff13 = lpMPC_dl_aff + 104;
lpMPC_FLOAT* lpMPC_dslbaff13 = lpMPC_ds_aff + 104;
lpMPC_FLOAT* lpMPC_dllbcc13 = lpMPC_dl_cc + 104;
lpMPC_FLOAT* lpMPC_dslbcc13 = lpMPC_ds_cc + 104;
lpMPC_FLOAT* lpMPC_ccrhsl13 = lpMPC_ccrhs + 104;
int lpMPC_ubIdx13[4] = {0, 1, 2, 3};
lpMPC_FLOAT* lpMPC_lub13 = lpMPC_l + 108;
lpMPC_FLOAT* lpMPC_sub13 = lpMPC_s + 108;
lpMPC_FLOAT* lpMPC_lubbysub13 = lpMPC_lbys + 108;
lpMPC_FLOAT lpMPC_riub13[4];
lpMPC_FLOAT* lpMPC_dlubaff13 = lpMPC_dl_aff + 108;
lpMPC_FLOAT* lpMPC_dsubaff13 = lpMPC_ds_aff + 108;
lpMPC_FLOAT* lpMPC_dlubcc13 = lpMPC_dl_cc + 108;
lpMPC_FLOAT* lpMPC_dsubcc13 = lpMPC_ds_cc + 108;
lpMPC_FLOAT* lpMPC_ccrhsub13 = lpMPC_ccrhs + 108;
lpMPC_FLOAT lpMPC_Phi13[4];
lpMPC_FLOAT lpMPC_W13[4];
lpMPC_FLOAT lpMPC_Ysd13[4];
lpMPC_FLOAT lpMPC_Lsd13[4];
lpMPC_FLOAT* lpMPC_z14 = lpMPC_z + 56;
lpMPC_FLOAT* lpMPC_dzaff14 = lpMPC_dz_aff + 56;
lpMPC_FLOAT* lpMPC_dzcc14 = lpMPC_dz_cc + 56;
lpMPC_FLOAT* lpMPC_rd14 = lpMPC_rd + 56;
lpMPC_FLOAT lpMPC_Lbyrd14[2];
lpMPC_FLOAT* lpMPC_grad_cost14 = lpMPC_grad_cost + 56;
lpMPC_FLOAT* lpMPC_grad_eq14 = lpMPC_grad_eq + 56;
lpMPC_FLOAT* lpMPC_grad_ineq14 = lpMPC_grad_ineq + 56;
lpMPC_FLOAT lpMPC_ctv14[2];
lpMPC_FLOAT* lpMPC_v14 = lpMPC_v + 28;
lpMPC_FLOAT lpMPC_re14[2];
lpMPC_FLOAT lpMPC_beta14[2];
lpMPC_FLOAT lpMPC_betacc14[2];
lpMPC_FLOAT* lpMPC_dvaff14 = lpMPC_dv_aff + 28;
lpMPC_FLOAT* lpMPC_dvcc14 = lpMPC_dv_cc + 28;
lpMPC_FLOAT lpMPC_V14[4];
lpMPC_FLOAT lpMPC_Yd14[3];
lpMPC_FLOAT lpMPC_Ld14[3];
lpMPC_FLOAT lpMPC_yy14[2];
lpMPC_FLOAT lpMPC_bmy14[2];
int lpMPC_lbIdx14[2] = {0, 1};
lpMPC_FLOAT* lpMPC_llb14 = lpMPC_l + 112;
lpMPC_FLOAT* lpMPC_slb14 = lpMPC_s + 112;
lpMPC_FLOAT* lpMPC_llbbyslb14 = lpMPC_lbys + 112;
lpMPC_FLOAT lpMPC_rilb14[2];
lpMPC_FLOAT* lpMPC_dllbaff14 = lpMPC_dl_aff + 112;
lpMPC_FLOAT* lpMPC_dslbaff14 = lpMPC_ds_aff + 112;
lpMPC_FLOAT* lpMPC_dllbcc14 = lpMPC_dl_cc + 112;
lpMPC_FLOAT* lpMPC_dslbcc14 = lpMPC_ds_cc + 112;
lpMPC_FLOAT* lpMPC_ccrhsl14 = lpMPC_ccrhs + 112;
int lpMPC_ubIdx14[2] = {0, 1};
lpMPC_FLOAT* lpMPC_lub14 = lpMPC_l + 114;
lpMPC_FLOAT* lpMPC_sub14 = lpMPC_s + 114;
lpMPC_FLOAT* lpMPC_lubbysub14 = lpMPC_lbys + 114;
lpMPC_FLOAT lpMPC_riub14[2];
lpMPC_FLOAT* lpMPC_dlubaff14 = lpMPC_dl_aff + 114;
lpMPC_FLOAT* lpMPC_dsubaff14 = lpMPC_ds_aff + 114;
lpMPC_FLOAT* lpMPC_dlubcc14 = lpMPC_dl_cc + 114;
lpMPC_FLOAT* lpMPC_dsubcc14 = lpMPC_ds_cc + 114;
lpMPC_FLOAT* lpMPC_ccrhsub14 = lpMPC_ccrhs + 114;
lpMPC_FLOAT lpMPC_Phi14[2];
lpMPC_FLOAT lpMPC_D14[2] = {-1.0000000000000000E+000, 
-1.0000000000000000E+000};
lpMPC_FLOAT lpMPC_W14[2];
lpMPC_FLOAT lpMPC_Ysd14[4];
lpMPC_FLOAT lpMPC_Lsd14[4];
lpMPC_FLOAT musigma;
lpMPC_FLOAT sigma_3rdroot;
lpMPC_FLOAT lpMPC_Diag1_0[4];
lpMPC_FLOAT lpMPC_Diag2_0[4];
lpMPC_FLOAT lpMPC_L_0[6];




/* SOLVER CODE --------------------------------------------------------- */
int lpMPC_solve(lpMPC_params* params, lpMPC_output* output, lpMPC_info* info)
{	
int exitcode;

#if lpMPC_SET_TIMING == 1
	lpMPC_timer solvertimer;
	lpMPC_tic(&solvertimer);
#endif
/* FUNCTION CALLS INTO LA LIBRARY -------------------------------------- */
info->it = 0;
lpMPC_LA_INITIALIZEVECTOR_58(lpMPC_z, 0);
lpMPC_LA_INITIALIZEVECTOR_30(lpMPC_v, 1);
lpMPC_LA_INITIALIZEVECTOR_116(lpMPC_l, 1);
lpMPC_LA_INITIALIZEVECTOR_116(lpMPC_s, 1);
info->mu = 0;
lpMPC_LA_DOTACC_116(lpMPC_l, lpMPC_s, &info->mu);
info->mu /= 116;
while( 1 ){
info->pobj = 0;
lpMPC_LA_DIAG_QUADFCN_4(lpMPC_H00, params->f1, lpMPC_z00, lpMPC_grad_cost00, &info->pobj);
lpMPC_LA_DIAG_QUADFCN_4(lpMPC_H00, params->f2, lpMPC_z01, lpMPC_grad_cost01, &info->pobj);
lpMPC_LA_DIAG_QUADFCN_4(lpMPC_H00, params->f3, lpMPC_z02, lpMPC_grad_cost02, &info->pobj);
lpMPC_LA_DIAG_QUADFCN_4(lpMPC_H00, params->f4, lpMPC_z03, lpMPC_grad_cost03, &info->pobj);
lpMPC_LA_DIAG_QUADFCN_4(lpMPC_H00, params->f5, lpMPC_z04, lpMPC_grad_cost04, &info->pobj);
lpMPC_LA_DIAG_QUADFCN_4(lpMPC_H00, params->f6, lpMPC_z05, lpMPC_grad_cost05, &info->pobj);
lpMPC_LA_DIAG_QUADFCN_4(lpMPC_H00, params->f7, lpMPC_z06, lpMPC_grad_cost06, &info->pobj);
lpMPC_LA_DIAG_QUADFCN_4(lpMPC_H00, params->f8, lpMPC_z07, lpMPC_grad_cost07, &info->pobj);
lpMPC_LA_DIAG_QUADFCN_4(lpMPC_H00, params->f9, lpMPC_z08, lpMPC_grad_cost08, &info->pobj);
lpMPC_LA_DIAG_QUADFCN_4(lpMPC_H00, params->f10, lpMPC_z09, lpMPC_grad_cost09, &info->pobj);
lpMPC_LA_DIAG_QUADFCN_4(lpMPC_H00, params->f11, lpMPC_z10, lpMPC_grad_cost10, &info->pobj);
lpMPC_LA_DIAG_QUADFCN_4(lpMPC_H00, params->f12, lpMPC_z11, lpMPC_grad_cost11, &info->pobj);
lpMPC_LA_DIAG_QUADFCN_4(lpMPC_H00, params->f13, lpMPC_z12, lpMPC_grad_cost12, &info->pobj);
lpMPC_LA_DIAG_QUADFCN_4(lpMPC_H00, params->f14, lpMPC_z13, lpMPC_grad_cost13, &info->pobj);
lpMPC_LA_DIAG_QUADFCN_2(lpMPC_H00, params->f15, lpMPC_z14, lpMPC_grad_cost14, &info->pobj);
info->res_eq = 0;
info->dgap = 0;
lpMPC_LA_DIAGZERO_MVMSUB6_2(lpMPC_D00, lpMPC_z00, params->e1, lpMPC_v00, lpMPC_re00, &info->dgap, &info->res_eq);
lpMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_4_4(params->C1, lpMPC_z00, lpMPC_D01, lpMPC_z01, params->e2, lpMPC_v01, lpMPC_re01, &info->dgap, &info->res_eq);
lpMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_4_4(params->C2, lpMPC_z01, lpMPC_D01, lpMPC_z02, params->e3, lpMPC_v02, lpMPC_re02, &info->dgap, &info->res_eq);
lpMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_4_4(params->C3, lpMPC_z02, lpMPC_D01, lpMPC_z03, params->e4, lpMPC_v03, lpMPC_re03, &info->dgap, &info->res_eq);
lpMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_4_4(params->C4, lpMPC_z03, lpMPC_D01, lpMPC_z04, params->e5, lpMPC_v04, lpMPC_re04, &info->dgap, &info->res_eq);
lpMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_4_4(params->C5, lpMPC_z04, lpMPC_D01, lpMPC_z05, params->e6, lpMPC_v05, lpMPC_re05, &info->dgap, &info->res_eq);
lpMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_4_4(params->C6, lpMPC_z05, lpMPC_D01, lpMPC_z06, params->e7, lpMPC_v06, lpMPC_re06, &info->dgap, &info->res_eq);
lpMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_4_4(params->C7, lpMPC_z06, lpMPC_D01, lpMPC_z07, params->e8, lpMPC_v07, lpMPC_re07, &info->dgap, &info->res_eq);
lpMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_4_4(params->C8, lpMPC_z07, lpMPC_D01, lpMPC_z08, params->e9, lpMPC_v08, lpMPC_re08, &info->dgap, &info->res_eq);
lpMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_4_4(params->C9, lpMPC_z08, lpMPC_D01, lpMPC_z09, params->e10, lpMPC_v09, lpMPC_re09, &info->dgap, &info->res_eq);
lpMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_4_4(params->C10, lpMPC_z09, lpMPC_D01, lpMPC_z10, params->e11, lpMPC_v10, lpMPC_re10, &info->dgap, &info->res_eq);
lpMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_4_4(params->C11, lpMPC_z10, lpMPC_D01, lpMPC_z11, params->e12, lpMPC_v11, lpMPC_re11, &info->dgap, &info->res_eq);
lpMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_4_4(params->C12, lpMPC_z11, lpMPC_D01, lpMPC_z12, params->e13, lpMPC_v12, lpMPC_re12, &info->dgap, &info->res_eq);
lpMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_4_4(params->C13, lpMPC_z12, lpMPC_D01, lpMPC_z13, params->e14, lpMPC_v13, lpMPC_re13, &info->dgap, &info->res_eq);
lpMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_4_2(params->C14, lpMPC_z13, lpMPC_D14, lpMPC_z14, params->e15, lpMPC_v14, lpMPC_re14, &info->dgap, &info->res_eq);
lpMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C1, lpMPC_v01, lpMPC_D00, lpMPC_v00, lpMPC_grad_eq00);
lpMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C2, lpMPC_v02, lpMPC_D01, lpMPC_v01, lpMPC_grad_eq01);
lpMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C3, lpMPC_v03, lpMPC_D01, lpMPC_v02, lpMPC_grad_eq02);
lpMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C4, lpMPC_v04, lpMPC_D01, lpMPC_v03, lpMPC_grad_eq03);
lpMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C5, lpMPC_v05, lpMPC_D01, lpMPC_v04, lpMPC_grad_eq04);
lpMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C6, lpMPC_v06, lpMPC_D01, lpMPC_v05, lpMPC_grad_eq05);
lpMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C7, lpMPC_v07, lpMPC_D01, lpMPC_v06, lpMPC_grad_eq06);
lpMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C8, lpMPC_v08, lpMPC_D01, lpMPC_v07, lpMPC_grad_eq07);
lpMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C9, lpMPC_v09, lpMPC_D01, lpMPC_v08, lpMPC_grad_eq08);
lpMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C10, lpMPC_v10, lpMPC_D01, lpMPC_v09, lpMPC_grad_eq09);
lpMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C11, lpMPC_v11, lpMPC_D01, lpMPC_v10, lpMPC_grad_eq10);
lpMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C12, lpMPC_v12, lpMPC_D01, lpMPC_v11, lpMPC_grad_eq11);
lpMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C13, lpMPC_v13, lpMPC_D01, lpMPC_v12, lpMPC_grad_eq12);
lpMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C14, lpMPC_v14, lpMPC_D01, lpMPC_v13, lpMPC_grad_eq13);
lpMPC_LA_DIAGZERO_MTVM_2_2(lpMPC_D14, lpMPC_v14, lpMPC_grad_eq14);
info->res_ineq = 0;
lpMPC_LA_VSUBADD3_4(params->lb1, lpMPC_z00, lpMPC_lbIdx00, lpMPC_llb00, lpMPC_slb00, lpMPC_rilb00, &info->dgap, &info->res_ineq);
lpMPC_LA_VSUBADD2_4(lpMPC_z00, lpMPC_ubIdx00, params->ub1, lpMPC_lub00, lpMPC_sub00, lpMPC_riub00, &info->dgap, &info->res_ineq);
lpMPC_LA_VSUBADD3_4(params->lb2, lpMPC_z01, lpMPC_lbIdx01, lpMPC_llb01, lpMPC_slb01, lpMPC_rilb01, &info->dgap, &info->res_ineq);
lpMPC_LA_VSUBADD2_4(lpMPC_z01, lpMPC_ubIdx01, params->ub2, lpMPC_lub01, lpMPC_sub01, lpMPC_riub01, &info->dgap, &info->res_ineq);
lpMPC_LA_VSUBADD3_4(params->lb3, lpMPC_z02, lpMPC_lbIdx02, lpMPC_llb02, lpMPC_slb02, lpMPC_rilb02, &info->dgap, &info->res_ineq);
lpMPC_LA_VSUBADD2_4(lpMPC_z02, lpMPC_ubIdx02, params->ub3, lpMPC_lub02, lpMPC_sub02, lpMPC_riub02, &info->dgap, &info->res_ineq);
lpMPC_LA_VSUBADD3_4(params->lb4, lpMPC_z03, lpMPC_lbIdx03, lpMPC_llb03, lpMPC_slb03, lpMPC_rilb03, &info->dgap, &info->res_ineq);
lpMPC_LA_VSUBADD2_4(lpMPC_z03, lpMPC_ubIdx03, params->ub4, lpMPC_lub03, lpMPC_sub03, lpMPC_riub03, &info->dgap, &info->res_ineq);
lpMPC_LA_VSUBADD3_4(params->lb5, lpMPC_z04, lpMPC_lbIdx04, lpMPC_llb04, lpMPC_slb04, lpMPC_rilb04, &info->dgap, &info->res_ineq);
lpMPC_LA_VSUBADD2_4(lpMPC_z04, lpMPC_ubIdx04, params->ub5, lpMPC_lub04, lpMPC_sub04, lpMPC_riub04, &info->dgap, &info->res_ineq);
lpMPC_LA_VSUBADD3_4(params->lb6, lpMPC_z05, lpMPC_lbIdx05, lpMPC_llb05, lpMPC_slb05, lpMPC_rilb05, &info->dgap, &info->res_ineq);
lpMPC_LA_VSUBADD2_4(lpMPC_z05, lpMPC_ubIdx05, params->ub6, lpMPC_lub05, lpMPC_sub05, lpMPC_riub05, &info->dgap, &info->res_ineq);
lpMPC_LA_VSUBADD3_4(params->lb7, lpMPC_z06, lpMPC_lbIdx06, lpMPC_llb06, lpMPC_slb06, lpMPC_rilb06, &info->dgap, &info->res_ineq);
lpMPC_LA_VSUBADD2_4(lpMPC_z06, lpMPC_ubIdx06, params->ub7, lpMPC_lub06, lpMPC_sub06, lpMPC_riub06, &info->dgap, &info->res_ineq);
lpMPC_LA_VSUBADD3_4(params->lb8, lpMPC_z07, lpMPC_lbIdx07, lpMPC_llb07, lpMPC_slb07, lpMPC_rilb07, &info->dgap, &info->res_ineq);
lpMPC_LA_VSUBADD2_4(lpMPC_z07, lpMPC_ubIdx07, params->ub8, lpMPC_lub07, lpMPC_sub07, lpMPC_riub07, &info->dgap, &info->res_ineq);
lpMPC_LA_VSUBADD3_4(params->lb9, lpMPC_z08, lpMPC_lbIdx08, lpMPC_llb08, lpMPC_slb08, lpMPC_rilb08, &info->dgap, &info->res_ineq);
lpMPC_LA_VSUBADD2_4(lpMPC_z08, lpMPC_ubIdx08, params->ub9, lpMPC_lub08, lpMPC_sub08, lpMPC_riub08, &info->dgap, &info->res_ineq);
lpMPC_LA_VSUBADD3_4(params->lb10, lpMPC_z09, lpMPC_lbIdx09, lpMPC_llb09, lpMPC_slb09, lpMPC_rilb09, &info->dgap, &info->res_ineq);
lpMPC_LA_VSUBADD2_4(lpMPC_z09, lpMPC_ubIdx09, params->ub10, lpMPC_lub09, lpMPC_sub09, lpMPC_riub09, &info->dgap, &info->res_ineq);
lpMPC_LA_VSUBADD3_4(params->lb11, lpMPC_z10, lpMPC_lbIdx10, lpMPC_llb10, lpMPC_slb10, lpMPC_rilb10, &info->dgap, &info->res_ineq);
lpMPC_LA_VSUBADD2_4(lpMPC_z10, lpMPC_ubIdx10, params->ub11, lpMPC_lub10, lpMPC_sub10, lpMPC_riub10, &info->dgap, &info->res_ineq);
lpMPC_LA_VSUBADD3_4(params->lb12, lpMPC_z11, lpMPC_lbIdx11, lpMPC_llb11, lpMPC_slb11, lpMPC_rilb11, &info->dgap, &info->res_ineq);
lpMPC_LA_VSUBADD2_4(lpMPC_z11, lpMPC_ubIdx11, params->ub12, lpMPC_lub11, lpMPC_sub11, lpMPC_riub11, &info->dgap, &info->res_ineq);
lpMPC_LA_VSUBADD3_4(params->lb13, lpMPC_z12, lpMPC_lbIdx12, lpMPC_llb12, lpMPC_slb12, lpMPC_rilb12, &info->dgap, &info->res_ineq);
lpMPC_LA_VSUBADD2_4(lpMPC_z12, lpMPC_ubIdx12, params->ub13, lpMPC_lub12, lpMPC_sub12, lpMPC_riub12, &info->dgap, &info->res_ineq);
lpMPC_LA_VSUBADD3_4(params->lb14, lpMPC_z13, lpMPC_lbIdx13, lpMPC_llb13, lpMPC_slb13, lpMPC_rilb13, &info->dgap, &info->res_ineq);
lpMPC_LA_VSUBADD2_4(lpMPC_z13, lpMPC_ubIdx13, params->ub14, lpMPC_lub13, lpMPC_sub13, lpMPC_riub13, &info->dgap, &info->res_ineq);
lpMPC_LA_VSUBADD3_2(params->lb15, lpMPC_z14, lpMPC_lbIdx14, lpMPC_llb14, lpMPC_slb14, lpMPC_rilb14, &info->dgap, &info->res_ineq);
lpMPC_LA_VSUBADD2_2(lpMPC_z14, lpMPC_ubIdx14, params->ub15, lpMPC_lub14, lpMPC_sub14, lpMPC_riub14, &info->dgap, &info->res_ineq);
lpMPC_LA_INEQ_B_GRAD_4_4_4(lpMPC_lub00, lpMPC_sub00, lpMPC_riub00, lpMPC_llb00, lpMPC_slb00, lpMPC_rilb00, lpMPC_lbIdx00, lpMPC_ubIdx00, lpMPC_grad_ineq00, lpMPC_lubbysub00, lpMPC_llbbyslb00);
lpMPC_LA_INEQ_B_GRAD_4_4_4(lpMPC_lub01, lpMPC_sub01, lpMPC_riub01, lpMPC_llb01, lpMPC_slb01, lpMPC_rilb01, lpMPC_lbIdx01, lpMPC_ubIdx01, lpMPC_grad_ineq01, lpMPC_lubbysub01, lpMPC_llbbyslb01);
lpMPC_LA_INEQ_B_GRAD_4_4_4(lpMPC_lub02, lpMPC_sub02, lpMPC_riub02, lpMPC_llb02, lpMPC_slb02, lpMPC_rilb02, lpMPC_lbIdx02, lpMPC_ubIdx02, lpMPC_grad_ineq02, lpMPC_lubbysub02, lpMPC_llbbyslb02);
lpMPC_LA_INEQ_B_GRAD_4_4_4(lpMPC_lub03, lpMPC_sub03, lpMPC_riub03, lpMPC_llb03, lpMPC_slb03, lpMPC_rilb03, lpMPC_lbIdx03, lpMPC_ubIdx03, lpMPC_grad_ineq03, lpMPC_lubbysub03, lpMPC_llbbyslb03);
lpMPC_LA_INEQ_B_GRAD_4_4_4(lpMPC_lub04, lpMPC_sub04, lpMPC_riub04, lpMPC_llb04, lpMPC_slb04, lpMPC_rilb04, lpMPC_lbIdx04, lpMPC_ubIdx04, lpMPC_grad_ineq04, lpMPC_lubbysub04, lpMPC_llbbyslb04);
lpMPC_LA_INEQ_B_GRAD_4_4_4(lpMPC_lub05, lpMPC_sub05, lpMPC_riub05, lpMPC_llb05, lpMPC_slb05, lpMPC_rilb05, lpMPC_lbIdx05, lpMPC_ubIdx05, lpMPC_grad_ineq05, lpMPC_lubbysub05, lpMPC_llbbyslb05);
lpMPC_LA_INEQ_B_GRAD_4_4_4(lpMPC_lub06, lpMPC_sub06, lpMPC_riub06, lpMPC_llb06, lpMPC_slb06, lpMPC_rilb06, lpMPC_lbIdx06, lpMPC_ubIdx06, lpMPC_grad_ineq06, lpMPC_lubbysub06, lpMPC_llbbyslb06);
lpMPC_LA_INEQ_B_GRAD_4_4_4(lpMPC_lub07, lpMPC_sub07, lpMPC_riub07, lpMPC_llb07, lpMPC_slb07, lpMPC_rilb07, lpMPC_lbIdx07, lpMPC_ubIdx07, lpMPC_grad_ineq07, lpMPC_lubbysub07, lpMPC_llbbyslb07);
lpMPC_LA_INEQ_B_GRAD_4_4_4(lpMPC_lub08, lpMPC_sub08, lpMPC_riub08, lpMPC_llb08, lpMPC_slb08, lpMPC_rilb08, lpMPC_lbIdx08, lpMPC_ubIdx08, lpMPC_grad_ineq08, lpMPC_lubbysub08, lpMPC_llbbyslb08);
lpMPC_LA_INEQ_B_GRAD_4_4_4(lpMPC_lub09, lpMPC_sub09, lpMPC_riub09, lpMPC_llb09, lpMPC_slb09, lpMPC_rilb09, lpMPC_lbIdx09, lpMPC_ubIdx09, lpMPC_grad_ineq09, lpMPC_lubbysub09, lpMPC_llbbyslb09);
lpMPC_LA_INEQ_B_GRAD_4_4_4(lpMPC_lub10, lpMPC_sub10, lpMPC_riub10, lpMPC_llb10, lpMPC_slb10, lpMPC_rilb10, lpMPC_lbIdx10, lpMPC_ubIdx10, lpMPC_grad_ineq10, lpMPC_lubbysub10, lpMPC_llbbyslb10);
lpMPC_LA_INEQ_B_GRAD_4_4_4(lpMPC_lub11, lpMPC_sub11, lpMPC_riub11, lpMPC_llb11, lpMPC_slb11, lpMPC_rilb11, lpMPC_lbIdx11, lpMPC_ubIdx11, lpMPC_grad_ineq11, lpMPC_lubbysub11, lpMPC_llbbyslb11);
lpMPC_LA_INEQ_B_GRAD_4_4_4(lpMPC_lub12, lpMPC_sub12, lpMPC_riub12, lpMPC_llb12, lpMPC_slb12, lpMPC_rilb12, lpMPC_lbIdx12, lpMPC_ubIdx12, lpMPC_grad_ineq12, lpMPC_lubbysub12, lpMPC_llbbyslb12);
lpMPC_LA_INEQ_B_GRAD_4_4_4(lpMPC_lub13, lpMPC_sub13, lpMPC_riub13, lpMPC_llb13, lpMPC_slb13, lpMPC_rilb13, lpMPC_lbIdx13, lpMPC_ubIdx13, lpMPC_grad_ineq13, lpMPC_lubbysub13, lpMPC_llbbyslb13);
lpMPC_LA_INEQ_B_GRAD_2_2_2(lpMPC_lub14, lpMPC_sub14, lpMPC_riub14, lpMPC_llb14, lpMPC_slb14, lpMPC_rilb14, lpMPC_lbIdx14, lpMPC_ubIdx14, lpMPC_grad_ineq14, lpMPC_lubbysub14, lpMPC_llbbyslb14);
info->dobj = info->pobj - info->dgap;
info->rdgap = info->pobj ? info->dgap / info->pobj : 1e6;
if( info->rdgap < 0 ) info->rdgap = -info->rdgap;
if( info->mu < lpMPC_SET_ACC_KKTCOMPL
    && (info->rdgap < lpMPC_SET_ACC_RDGAP || info->dgap < lpMPC_SET_ACC_KKTCOMPL)
    && info->res_eq < lpMPC_SET_ACC_RESEQ
    && info->res_ineq < lpMPC_SET_ACC_RESINEQ ){
exitcode = lpMPC_OPTIMAL; break; }
if( info->it == lpMPC_SET_MAXIT ){
exitcode = lpMPC_MAXITREACHED; break; }
lpMPC_LA_VVADD3_58(lpMPC_grad_cost, lpMPC_grad_eq, lpMPC_grad_ineq, lpMPC_rd);
lpMPC_LA_DIAG_CHOL_ONELOOP_LBUB_4_4_4(lpMPC_H00, lpMPC_llbbyslb00, lpMPC_lbIdx00, lpMPC_lubbysub00, lpMPC_ubIdx00, lpMPC_Phi00);
lpMPC_LA_DIAG_MATRIXFORWARDSUB_2_4(lpMPC_Phi00, params->C1, lpMPC_V00);
lpMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_4(lpMPC_Phi00, lpMPC_D00, lpMPC_W00);
lpMPC_LA_DENSE_DIAGZERO_MMTM_2_4_2(lpMPC_W00, lpMPC_V00, lpMPC_Ysd01);
lpMPC_LA_DIAG_FORWARDSUB_4(lpMPC_Phi00, lpMPC_rd00, lpMPC_Lbyrd00);
lpMPC_LA_DIAG_CHOL_ONELOOP_LBUB_4_4_4(lpMPC_H00, lpMPC_llbbyslb01, lpMPC_lbIdx01, lpMPC_lubbysub01, lpMPC_ubIdx01, lpMPC_Phi01);
lpMPC_LA_DIAG_MATRIXFORWARDSUB_2_4(lpMPC_Phi01, params->C2, lpMPC_V01);
lpMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_4(lpMPC_Phi01, lpMPC_D01, lpMPC_W01);
lpMPC_LA_DENSE_DIAGZERO_MMTM_2_4_2(lpMPC_W01, lpMPC_V01, lpMPC_Ysd02);
lpMPC_LA_DIAG_FORWARDSUB_4(lpMPC_Phi01, lpMPC_rd01, lpMPC_Lbyrd01);
lpMPC_LA_DIAG_CHOL_ONELOOP_LBUB_4_4_4(lpMPC_H00, lpMPC_llbbyslb02, lpMPC_lbIdx02, lpMPC_lubbysub02, lpMPC_ubIdx02, lpMPC_Phi02);
lpMPC_LA_DIAG_MATRIXFORWARDSUB_2_4(lpMPC_Phi02, params->C3, lpMPC_V02);
lpMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_4(lpMPC_Phi02, lpMPC_D01, lpMPC_W02);
lpMPC_LA_DENSE_DIAGZERO_MMTM_2_4_2(lpMPC_W02, lpMPC_V02, lpMPC_Ysd03);
lpMPC_LA_DIAG_FORWARDSUB_4(lpMPC_Phi02, lpMPC_rd02, lpMPC_Lbyrd02);
lpMPC_LA_DIAG_CHOL_ONELOOP_LBUB_4_4_4(lpMPC_H00, lpMPC_llbbyslb03, lpMPC_lbIdx03, lpMPC_lubbysub03, lpMPC_ubIdx03, lpMPC_Phi03);
lpMPC_LA_DIAG_MATRIXFORWARDSUB_2_4(lpMPC_Phi03, params->C4, lpMPC_V03);
lpMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_4(lpMPC_Phi03, lpMPC_D01, lpMPC_W03);
lpMPC_LA_DENSE_DIAGZERO_MMTM_2_4_2(lpMPC_W03, lpMPC_V03, lpMPC_Ysd04);
lpMPC_LA_DIAG_FORWARDSUB_4(lpMPC_Phi03, lpMPC_rd03, lpMPC_Lbyrd03);
lpMPC_LA_DIAG_CHOL_ONELOOP_LBUB_4_4_4(lpMPC_H00, lpMPC_llbbyslb04, lpMPC_lbIdx04, lpMPC_lubbysub04, lpMPC_ubIdx04, lpMPC_Phi04);
lpMPC_LA_DIAG_MATRIXFORWARDSUB_2_4(lpMPC_Phi04, params->C5, lpMPC_V04);
lpMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_4(lpMPC_Phi04, lpMPC_D01, lpMPC_W04);
lpMPC_LA_DENSE_DIAGZERO_MMTM_2_4_2(lpMPC_W04, lpMPC_V04, lpMPC_Ysd05);
lpMPC_LA_DIAG_FORWARDSUB_4(lpMPC_Phi04, lpMPC_rd04, lpMPC_Lbyrd04);
lpMPC_LA_DIAG_CHOL_ONELOOP_LBUB_4_4_4(lpMPC_H00, lpMPC_llbbyslb05, lpMPC_lbIdx05, lpMPC_lubbysub05, lpMPC_ubIdx05, lpMPC_Phi05);
lpMPC_LA_DIAG_MATRIXFORWARDSUB_2_4(lpMPC_Phi05, params->C6, lpMPC_V05);
lpMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_4(lpMPC_Phi05, lpMPC_D01, lpMPC_W05);
lpMPC_LA_DENSE_DIAGZERO_MMTM_2_4_2(lpMPC_W05, lpMPC_V05, lpMPC_Ysd06);
lpMPC_LA_DIAG_FORWARDSUB_4(lpMPC_Phi05, lpMPC_rd05, lpMPC_Lbyrd05);
lpMPC_LA_DIAG_CHOL_ONELOOP_LBUB_4_4_4(lpMPC_H00, lpMPC_llbbyslb06, lpMPC_lbIdx06, lpMPC_lubbysub06, lpMPC_ubIdx06, lpMPC_Phi06);
lpMPC_LA_DIAG_MATRIXFORWARDSUB_2_4(lpMPC_Phi06, params->C7, lpMPC_V06);
lpMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_4(lpMPC_Phi06, lpMPC_D01, lpMPC_W06);
lpMPC_LA_DENSE_DIAGZERO_MMTM_2_4_2(lpMPC_W06, lpMPC_V06, lpMPC_Ysd07);
lpMPC_LA_DIAG_FORWARDSUB_4(lpMPC_Phi06, lpMPC_rd06, lpMPC_Lbyrd06);
lpMPC_LA_DIAG_CHOL_ONELOOP_LBUB_4_4_4(lpMPC_H00, lpMPC_llbbyslb07, lpMPC_lbIdx07, lpMPC_lubbysub07, lpMPC_ubIdx07, lpMPC_Phi07);
lpMPC_LA_DIAG_MATRIXFORWARDSUB_2_4(lpMPC_Phi07, params->C8, lpMPC_V07);
lpMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_4(lpMPC_Phi07, lpMPC_D01, lpMPC_W07);
lpMPC_LA_DENSE_DIAGZERO_MMTM_2_4_2(lpMPC_W07, lpMPC_V07, lpMPC_Ysd08);
lpMPC_LA_DIAG_FORWARDSUB_4(lpMPC_Phi07, lpMPC_rd07, lpMPC_Lbyrd07);
lpMPC_LA_DIAG_CHOL_ONELOOP_LBUB_4_4_4(lpMPC_H00, lpMPC_llbbyslb08, lpMPC_lbIdx08, lpMPC_lubbysub08, lpMPC_ubIdx08, lpMPC_Phi08);
lpMPC_LA_DIAG_MATRIXFORWARDSUB_2_4(lpMPC_Phi08, params->C9, lpMPC_V08);
lpMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_4(lpMPC_Phi08, lpMPC_D01, lpMPC_W08);
lpMPC_LA_DENSE_DIAGZERO_MMTM_2_4_2(lpMPC_W08, lpMPC_V08, lpMPC_Ysd09);
lpMPC_LA_DIAG_FORWARDSUB_4(lpMPC_Phi08, lpMPC_rd08, lpMPC_Lbyrd08);
lpMPC_LA_DIAG_CHOL_ONELOOP_LBUB_4_4_4(lpMPC_H00, lpMPC_llbbyslb09, lpMPC_lbIdx09, lpMPC_lubbysub09, lpMPC_ubIdx09, lpMPC_Phi09);
lpMPC_LA_DIAG_MATRIXFORWARDSUB_2_4(lpMPC_Phi09, params->C10, lpMPC_V09);
lpMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_4(lpMPC_Phi09, lpMPC_D01, lpMPC_W09);
lpMPC_LA_DENSE_DIAGZERO_MMTM_2_4_2(lpMPC_W09, lpMPC_V09, lpMPC_Ysd10);
lpMPC_LA_DIAG_FORWARDSUB_4(lpMPC_Phi09, lpMPC_rd09, lpMPC_Lbyrd09);
lpMPC_LA_DIAG_CHOL_ONELOOP_LBUB_4_4_4(lpMPC_H00, lpMPC_llbbyslb10, lpMPC_lbIdx10, lpMPC_lubbysub10, lpMPC_ubIdx10, lpMPC_Phi10);
lpMPC_LA_DIAG_MATRIXFORWARDSUB_2_4(lpMPC_Phi10, params->C11, lpMPC_V10);
lpMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_4(lpMPC_Phi10, lpMPC_D01, lpMPC_W10);
lpMPC_LA_DENSE_DIAGZERO_MMTM_2_4_2(lpMPC_W10, lpMPC_V10, lpMPC_Ysd11);
lpMPC_LA_DIAG_FORWARDSUB_4(lpMPC_Phi10, lpMPC_rd10, lpMPC_Lbyrd10);
lpMPC_LA_DIAG_CHOL_ONELOOP_LBUB_4_4_4(lpMPC_H00, lpMPC_llbbyslb11, lpMPC_lbIdx11, lpMPC_lubbysub11, lpMPC_ubIdx11, lpMPC_Phi11);
lpMPC_LA_DIAG_MATRIXFORWARDSUB_2_4(lpMPC_Phi11, params->C12, lpMPC_V11);
lpMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_4(lpMPC_Phi11, lpMPC_D01, lpMPC_W11);
lpMPC_LA_DENSE_DIAGZERO_MMTM_2_4_2(lpMPC_W11, lpMPC_V11, lpMPC_Ysd12);
lpMPC_LA_DIAG_FORWARDSUB_4(lpMPC_Phi11, lpMPC_rd11, lpMPC_Lbyrd11);
lpMPC_LA_DIAG_CHOL_ONELOOP_LBUB_4_4_4(lpMPC_H00, lpMPC_llbbyslb12, lpMPC_lbIdx12, lpMPC_lubbysub12, lpMPC_ubIdx12, lpMPC_Phi12);
lpMPC_LA_DIAG_MATRIXFORWARDSUB_2_4(lpMPC_Phi12, params->C13, lpMPC_V12);
lpMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_4(lpMPC_Phi12, lpMPC_D01, lpMPC_W12);
lpMPC_LA_DENSE_DIAGZERO_MMTM_2_4_2(lpMPC_W12, lpMPC_V12, lpMPC_Ysd13);
lpMPC_LA_DIAG_FORWARDSUB_4(lpMPC_Phi12, lpMPC_rd12, lpMPC_Lbyrd12);
lpMPC_LA_DIAG_CHOL_ONELOOP_LBUB_4_4_4(lpMPC_H00, lpMPC_llbbyslb13, lpMPC_lbIdx13, lpMPC_lubbysub13, lpMPC_ubIdx13, lpMPC_Phi13);
lpMPC_LA_DIAG_MATRIXFORWARDSUB_2_4(lpMPC_Phi13, params->C14, lpMPC_V13);
lpMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_4(lpMPC_Phi13, lpMPC_D01, lpMPC_W13);
lpMPC_LA_DENSE_DIAGZERO_MMTM_2_4_2(lpMPC_W13, lpMPC_V13, lpMPC_Ysd14);
lpMPC_LA_DIAG_FORWARDSUB_4(lpMPC_Phi13, lpMPC_rd13, lpMPC_Lbyrd13);
lpMPC_LA_DIAG_CHOL_ONELOOP_LBUB_2_2_2(lpMPC_H00, lpMPC_llbbyslb14, lpMPC_lbIdx14, lpMPC_lubbysub14, lpMPC_ubIdx14, lpMPC_Phi14);
lpMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_2(lpMPC_Phi14, lpMPC_D14, lpMPC_W14);
lpMPC_LA_DIAG_FORWARDSUB_2(lpMPC_Phi14, lpMPC_rd14, lpMPC_Lbyrd14);
lpMPC_LA_DIAGZERO_MMT_2(lpMPC_W00, lpMPC_Yd00);
lpMPC_LA_DIAGZERO_MVMSUB7_2(lpMPC_W00, lpMPC_Lbyrd00, lpMPC_re00, lpMPC_beta00);
lpMPC_LA_DENSE_DIAGZERO_MMT2_2_4_4(lpMPC_V00, lpMPC_W01, lpMPC_Yd01);
lpMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_4_4(lpMPC_V00, lpMPC_Lbyrd00, lpMPC_W01, lpMPC_Lbyrd01, lpMPC_re01, lpMPC_beta01);
lpMPC_LA_DENSE_DIAGZERO_MMT2_2_4_4(lpMPC_V01, lpMPC_W02, lpMPC_Yd02);
lpMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_4_4(lpMPC_V01, lpMPC_Lbyrd01, lpMPC_W02, lpMPC_Lbyrd02, lpMPC_re02, lpMPC_beta02);
lpMPC_LA_DENSE_DIAGZERO_MMT2_2_4_4(lpMPC_V02, lpMPC_W03, lpMPC_Yd03);
lpMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_4_4(lpMPC_V02, lpMPC_Lbyrd02, lpMPC_W03, lpMPC_Lbyrd03, lpMPC_re03, lpMPC_beta03);
lpMPC_LA_DENSE_DIAGZERO_MMT2_2_4_4(lpMPC_V03, lpMPC_W04, lpMPC_Yd04);
lpMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_4_4(lpMPC_V03, lpMPC_Lbyrd03, lpMPC_W04, lpMPC_Lbyrd04, lpMPC_re04, lpMPC_beta04);
lpMPC_LA_DENSE_DIAGZERO_MMT2_2_4_4(lpMPC_V04, lpMPC_W05, lpMPC_Yd05);
lpMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_4_4(lpMPC_V04, lpMPC_Lbyrd04, lpMPC_W05, lpMPC_Lbyrd05, lpMPC_re05, lpMPC_beta05);
lpMPC_LA_DENSE_DIAGZERO_MMT2_2_4_4(lpMPC_V05, lpMPC_W06, lpMPC_Yd06);
lpMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_4_4(lpMPC_V05, lpMPC_Lbyrd05, lpMPC_W06, lpMPC_Lbyrd06, lpMPC_re06, lpMPC_beta06);
lpMPC_LA_DENSE_DIAGZERO_MMT2_2_4_4(lpMPC_V06, lpMPC_W07, lpMPC_Yd07);
lpMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_4_4(lpMPC_V06, lpMPC_Lbyrd06, lpMPC_W07, lpMPC_Lbyrd07, lpMPC_re07, lpMPC_beta07);
lpMPC_LA_DENSE_DIAGZERO_MMT2_2_4_4(lpMPC_V07, lpMPC_W08, lpMPC_Yd08);
lpMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_4_4(lpMPC_V07, lpMPC_Lbyrd07, lpMPC_W08, lpMPC_Lbyrd08, lpMPC_re08, lpMPC_beta08);
lpMPC_LA_DENSE_DIAGZERO_MMT2_2_4_4(lpMPC_V08, lpMPC_W09, lpMPC_Yd09);
lpMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_4_4(lpMPC_V08, lpMPC_Lbyrd08, lpMPC_W09, lpMPC_Lbyrd09, lpMPC_re09, lpMPC_beta09);
lpMPC_LA_DENSE_DIAGZERO_MMT2_2_4_4(lpMPC_V09, lpMPC_W10, lpMPC_Yd10);
lpMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_4_4(lpMPC_V09, lpMPC_Lbyrd09, lpMPC_W10, lpMPC_Lbyrd10, lpMPC_re10, lpMPC_beta10);
lpMPC_LA_DENSE_DIAGZERO_MMT2_2_4_4(lpMPC_V10, lpMPC_W11, lpMPC_Yd11);
lpMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_4_4(lpMPC_V10, lpMPC_Lbyrd10, lpMPC_W11, lpMPC_Lbyrd11, lpMPC_re11, lpMPC_beta11);
lpMPC_LA_DENSE_DIAGZERO_MMT2_2_4_4(lpMPC_V11, lpMPC_W12, lpMPC_Yd12);
lpMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_4_4(lpMPC_V11, lpMPC_Lbyrd11, lpMPC_W12, lpMPC_Lbyrd12, lpMPC_re12, lpMPC_beta12);
lpMPC_LA_DENSE_DIAGZERO_MMT2_2_4_4(lpMPC_V12, lpMPC_W13, lpMPC_Yd13);
lpMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_4_4(lpMPC_V12, lpMPC_Lbyrd12, lpMPC_W13, lpMPC_Lbyrd13, lpMPC_re13, lpMPC_beta13);
lpMPC_LA_DENSE_DIAGZERO_MMT2_2_4_2(lpMPC_V13, lpMPC_W14, lpMPC_Yd14);
lpMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_4_2(lpMPC_V13, lpMPC_Lbyrd13, lpMPC_W14, lpMPC_Lbyrd14, lpMPC_re14, lpMPC_beta14);
lpMPC_LA_DENSE_CHOL_2(lpMPC_Yd00, lpMPC_Ld00);
lpMPC_LA_DENSE_FORWARDSUB_2(lpMPC_Ld00, lpMPC_beta00, lpMPC_yy00);
lpMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(lpMPC_Ld00, lpMPC_Ysd01, lpMPC_Lsd01);
lpMPC_LA_DENSE_MMTSUB_2_2(lpMPC_Lsd01, lpMPC_Yd01);
lpMPC_LA_DENSE_CHOL_2(lpMPC_Yd01, lpMPC_Ld01);
lpMPC_LA_DENSE_MVMSUB1_2_2(lpMPC_Lsd01, lpMPC_yy00, lpMPC_beta01, lpMPC_bmy01);
lpMPC_LA_DENSE_FORWARDSUB_2(lpMPC_Ld01, lpMPC_bmy01, lpMPC_yy01);
lpMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(lpMPC_Ld01, lpMPC_Ysd02, lpMPC_Lsd02);
lpMPC_LA_DENSE_MMTSUB_2_2(lpMPC_Lsd02, lpMPC_Yd02);
lpMPC_LA_DENSE_CHOL_2(lpMPC_Yd02, lpMPC_Ld02);
lpMPC_LA_DENSE_MVMSUB1_2_2(lpMPC_Lsd02, lpMPC_yy01, lpMPC_beta02, lpMPC_bmy02);
lpMPC_LA_DENSE_FORWARDSUB_2(lpMPC_Ld02, lpMPC_bmy02, lpMPC_yy02);
lpMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(lpMPC_Ld02, lpMPC_Ysd03, lpMPC_Lsd03);
lpMPC_LA_DENSE_MMTSUB_2_2(lpMPC_Lsd03, lpMPC_Yd03);
lpMPC_LA_DENSE_CHOL_2(lpMPC_Yd03, lpMPC_Ld03);
lpMPC_LA_DENSE_MVMSUB1_2_2(lpMPC_Lsd03, lpMPC_yy02, lpMPC_beta03, lpMPC_bmy03);
lpMPC_LA_DENSE_FORWARDSUB_2(lpMPC_Ld03, lpMPC_bmy03, lpMPC_yy03);
lpMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(lpMPC_Ld03, lpMPC_Ysd04, lpMPC_Lsd04);
lpMPC_LA_DENSE_MMTSUB_2_2(lpMPC_Lsd04, lpMPC_Yd04);
lpMPC_LA_DENSE_CHOL_2(lpMPC_Yd04, lpMPC_Ld04);
lpMPC_LA_DENSE_MVMSUB1_2_2(lpMPC_Lsd04, lpMPC_yy03, lpMPC_beta04, lpMPC_bmy04);
lpMPC_LA_DENSE_FORWARDSUB_2(lpMPC_Ld04, lpMPC_bmy04, lpMPC_yy04);
lpMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(lpMPC_Ld04, lpMPC_Ysd05, lpMPC_Lsd05);
lpMPC_LA_DENSE_MMTSUB_2_2(lpMPC_Lsd05, lpMPC_Yd05);
lpMPC_LA_DENSE_CHOL_2(lpMPC_Yd05, lpMPC_Ld05);
lpMPC_LA_DENSE_MVMSUB1_2_2(lpMPC_Lsd05, lpMPC_yy04, lpMPC_beta05, lpMPC_bmy05);
lpMPC_LA_DENSE_FORWARDSUB_2(lpMPC_Ld05, lpMPC_bmy05, lpMPC_yy05);
lpMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(lpMPC_Ld05, lpMPC_Ysd06, lpMPC_Lsd06);
lpMPC_LA_DENSE_MMTSUB_2_2(lpMPC_Lsd06, lpMPC_Yd06);
lpMPC_LA_DENSE_CHOL_2(lpMPC_Yd06, lpMPC_Ld06);
lpMPC_LA_DENSE_MVMSUB1_2_2(lpMPC_Lsd06, lpMPC_yy05, lpMPC_beta06, lpMPC_bmy06);
lpMPC_LA_DENSE_FORWARDSUB_2(lpMPC_Ld06, lpMPC_bmy06, lpMPC_yy06);
lpMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(lpMPC_Ld06, lpMPC_Ysd07, lpMPC_Lsd07);
lpMPC_LA_DENSE_MMTSUB_2_2(lpMPC_Lsd07, lpMPC_Yd07);
lpMPC_LA_DENSE_CHOL_2(lpMPC_Yd07, lpMPC_Ld07);
lpMPC_LA_DENSE_MVMSUB1_2_2(lpMPC_Lsd07, lpMPC_yy06, lpMPC_beta07, lpMPC_bmy07);
lpMPC_LA_DENSE_FORWARDSUB_2(lpMPC_Ld07, lpMPC_bmy07, lpMPC_yy07);
lpMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(lpMPC_Ld07, lpMPC_Ysd08, lpMPC_Lsd08);
lpMPC_LA_DENSE_MMTSUB_2_2(lpMPC_Lsd08, lpMPC_Yd08);
lpMPC_LA_DENSE_CHOL_2(lpMPC_Yd08, lpMPC_Ld08);
lpMPC_LA_DENSE_MVMSUB1_2_2(lpMPC_Lsd08, lpMPC_yy07, lpMPC_beta08, lpMPC_bmy08);
lpMPC_LA_DENSE_FORWARDSUB_2(lpMPC_Ld08, lpMPC_bmy08, lpMPC_yy08);
lpMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(lpMPC_Ld08, lpMPC_Ysd09, lpMPC_Lsd09);
lpMPC_LA_DENSE_MMTSUB_2_2(lpMPC_Lsd09, lpMPC_Yd09);
lpMPC_LA_DENSE_CHOL_2(lpMPC_Yd09, lpMPC_Ld09);
lpMPC_LA_DENSE_MVMSUB1_2_2(lpMPC_Lsd09, lpMPC_yy08, lpMPC_beta09, lpMPC_bmy09);
lpMPC_LA_DENSE_FORWARDSUB_2(lpMPC_Ld09, lpMPC_bmy09, lpMPC_yy09);
lpMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(lpMPC_Ld09, lpMPC_Ysd10, lpMPC_Lsd10);
lpMPC_LA_DENSE_MMTSUB_2_2(lpMPC_Lsd10, lpMPC_Yd10);
lpMPC_LA_DENSE_CHOL_2(lpMPC_Yd10, lpMPC_Ld10);
lpMPC_LA_DENSE_MVMSUB1_2_2(lpMPC_Lsd10, lpMPC_yy09, lpMPC_beta10, lpMPC_bmy10);
lpMPC_LA_DENSE_FORWARDSUB_2(lpMPC_Ld10, lpMPC_bmy10, lpMPC_yy10);
lpMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(lpMPC_Ld10, lpMPC_Ysd11, lpMPC_Lsd11);
lpMPC_LA_DENSE_MMTSUB_2_2(lpMPC_Lsd11, lpMPC_Yd11);
lpMPC_LA_DENSE_CHOL_2(lpMPC_Yd11, lpMPC_Ld11);
lpMPC_LA_DENSE_MVMSUB1_2_2(lpMPC_Lsd11, lpMPC_yy10, lpMPC_beta11, lpMPC_bmy11);
lpMPC_LA_DENSE_FORWARDSUB_2(lpMPC_Ld11, lpMPC_bmy11, lpMPC_yy11);
lpMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(lpMPC_Ld11, lpMPC_Ysd12, lpMPC_Lsd12);
lpMPC_LA_DENSE_MMTSUB_2_2(lpMPC_Lsd12, lpMPC_Yd12);
lpMPC_LA_DENSE_CHOL_2(lpMPC_Yd12, lpMPC_Ld12);
lpMPC_LA_DENSE_MVMSUB1_2_2(lpMPC_Lsd12, lpMPC_yy11, lpMPC_beta12, lpMPC_bmy12);
lpMPC_LA_DENSE_FORWARDSUB_2(lpMPC_Ld12, lpMPC_bmy12, lpMPC_yy12);
lpMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(lpMPC_Ld12, lpMPC_Ysd13, lpMPC_Lsd13);
lpMPC_LA_DENSE_MMTSUB_2_2(lpMPC_Lsd13, lpMPC_Yd13);
lpMPC_LA_DENSE_CHOL_2(lpMPC_Yd13, lpMPC_Ld13);
lpMPC_LA_DENSE_MVMSUB1_2_2(lpMPC_Lsd13, lpMPC_yy12, lpMPC_beta13, lpMPC_bmy13);
lpMPC_LA_DENSE_FORWARDSUB_2(lpMPC_Ld13, lpMPC_bmy13, lpMPC_yy13);
lpMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(lpMPC_Ld13, lpMPC_Ysd14, lpMPC_Lsd14);
lpMPC_LA_DENSE_MMTSUB_2_2(lpMPC_Lsd14, lpMPC_Yd14);
lpMPC_LA_DENSE_CHOL_2(lpMPC_Yd14, lpMPC_Ld14);
lpMPC_LA_DENSE_MVMSUB1_2_2(lpMPC_Lsd14, lpMPC_yy13, lpMPC_beta14, lpMPC_bmy14);
lpMPC_LA_DENSE_FORWARDSUB_2(lpMPC_Ld14, lpMPC_bmy14, lpMPC_yy14);
lpMPC_LA_DENSE_BACKWARDSUB_2(lpMPC_Ld14, lpMPC_yy14, lpMPC_dvaff14);
lpMPC_LA_DENSE_MTVMSUB_2_2(lpMPC_Lsd14, lpMPC_dvaff14, lpMPC_yy13, lpMPC_bmy13);
lpMPC_LA_DENSE_BACKWARDSUB_2(lpMPC_Ld13, lpMPC_bmy13, lpMPC_dvaff13);
lpMPC_LA_DENSE_MTVMSUB_2_2(lpMPC_Lsd13, lpMPC_dvaff13, lpMPC_yy12, lpMPC_bmy12);
lpMPC_LA_DENSE_BACKWARDSUB_2(lpMPC_Ld12, lpMPC_bmy12, lpMPC_dvaff12);
lpMPC_LA_DENSE_MTVMSUB_2_2(lpMPC_Lsd12, lpMPC_dvaff12, lpMPC_yy11, lpMPC_bmy11);
lpMPC_LA_DENSE_BACKWARDSUB_2(lpMPC_Ld11, lpMPC_bmy11, lpMPC_dvaff11);
lpMPC_LA_DENSE_MTVMSUB_2_2(lpMPC_Lsd11, lpMPC_dvaff11, lpMPC_yy10, lpMPC_bmy10);
lpMPC_LA_DENSE_BACKWARDSUB_2(lpMPC_Ld10, lpMPC_bmy10, lpMPC_dvaff10);
lpMPC_LA_DENSE_MTVMSUB_2_2(lpMPC_Lsd10, lpMPC_dvaff10, lpMPC_yy09, lpMPC_bmy09);
lpMPC_LA_DENSE_BACKWARDSUB_2(lpMPC_Ld09, lpMPC_bmy09, lpMPC_dvaff09);
lpMPC_LA_DENSE_MTVMSUB_2_2(lpMPC_Lsd09, lpMPC_dvaff09, lpMPC_yy08, lpMPC_bmy08);
lpMPC_LA_DENSE_BACKWARDSUB_2(lpMPC_Ld08, lpMPC_bmy08, lpMPC_dvaff08);
lpMPC_LA_DENSE_MTVMSUB_2_2(lpMPC_Lsd08, lpMPC_dvaff08, lpMPC_yy07, lpMPC_bmy07);
lpMPC_LA_DENSE_BACKWARDSUB_2(lpMPC_Ld07, lpMPC_bmy07, lpMPC_dvaff07);
lpMPC_LA_DENSE_MTVMSUB_2_2(lpMPC_Lsd07, lpMPC_dvaff07, lpMPC_yy06, lpMPC_bmy06);
lpMPC_LA_DENSE_BACKWARDSUB_2(lpMPC_Ld06, lpMPC_bmy06, lpMPC_dvaff06);
lpMPC_LA_DENSE_MTVMSUB_2_2(lpMPC_Lsd06, lpMPC_dvaff06, lpMPC_yy05, lpMPC_bmy05);
lpMPC_LA_DENSE_BACKWARDSUB_2(lpMPC_Ld05, lpMPC_bmy05, lpMPC_dvaff05);
lpMPC_LA_DENSE_MTVMSUB_2_2(lpMPC_Lsd05, lpMPC_dvaff05, lpMPC_yy04, lpMPC_bmy04);
lpMPC_LA_DENSE_BACKWARDSUB_2(lpMPC_Ld04, lpMPC_bmy04, lpMPC_dvaff04);
lpMPC_LA_DENSE_MTVMSUB_2_2(lpMPC_Lsd04, lpMPC_dvaff04, lpMPC_yy03, lpMPC_bmy03);
lpMPC_LA_DENSE_BACKWARDSUB_2(lpMPC_Ld03, lpMPC_bmy03, lpMPC_dvaff03);
lpMPC_LA_DENSE_MTVMSUB_2_2(lpMPC_Lsd03, lpMPC_dvaff03, lpMPC_yy02, lpMPC_bmy02);
lpMPC_LA_DENSE_BACKWARDSUB_2(lpMPC_Ld02, lpMPC_bmy02, lpMPC_dvaff02);
lpMPC_LA_DENSE_MTVMSUB_2_2(lpMPC_Lsd02, lpMPC_dvaff02, lpMPC_yy01, lpMPC_bmy01);
lpMPC_LA_DENSE_BACKWARDSUB_2(lpMPC_Ld01, lpMPC_bmy01, lpMPC_dvaff01);
lpMPC_LA_DENSE_MTVMSUB_2_2(lpMPC_Lsd01, lpMPC_dvaff01, lpMPC_yy00, lpMPC_bmy00);
lpMPC_LA_DENSE_BACKWARDSUB_2(lpMPC_Ld00, lpMPC_bmy00, lpMPC_dvaff00);
lpMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C1, lpMPC_dvaff01, lpMPC_D00, lpMPC_dvaff00, lpMPC_grad_eq00);
lpMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C2, lpMPC_dvaff02, lpMPC_D01, lpMPC_dvaff01, lpMPC_grad_eq01);
lpMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C3, lpMPC_dvaff03, lpMPC_D01, lpMPC_dvaff02, lpMPC_grad_eq02);
lpMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C4, lpMPC_dvaff04, lpMPC_D01, lpMPC_dvaff03, lpMPC_grad_eq03);
lpMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C5, lpMPC_dvaff05, lpMPC_D01, lpMPC_dvaff04, lpMPC_grad_eq04);
lpMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C6, lpMPC_dvaff06, lpMPC_D01, lpMPC_dvaff05, lpMPC_grad_eq05);
lpMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C7, lpMPC_dvaff07, lpMPC_D01, lpMPC_dvaff06, lpMPC_grad_eq06);
lpMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C8, lpMPC_dvaff08, lpMPC_D01, lpMPC_dvaff07, lpMPC_grad_eq07);
lpMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C9, lpMPC_dvaff09, lpMPC_D01, lpMPC_dvaff08, lpMPC_grad_eq08);
lpMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C10, lpMPC_dvaff10, lpMPC_D01, lpMPC_dvaff09, lpMPC_grad_eq09);
lpMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C11, lpMPC_dvaff11, lpMPC_D01, lpMPC_dvaff10, lpMPC_grad_eq10);
lpMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C12, lpMPC_dvaff12, lpMPC_D01, lpMPC_dvaff11, lpMPC_grad_eq11);
lpMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C13, lpMPC_dvaff13, lpMPC_D01, lpMPC_dvaff12, lpMPC_grad_eq12);
lpMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C14, lpMPC_dvaff14, lpMPC_D01, lpMPC_dvaff13, lpMPC_grad_eq13);
lpMPC_LA_DIAGZERO_MTVM_2_2(lpMPC_D14, lpMPC_dvaff14, lpMPC_grad_eq14);
lpMPC_LA_VSUB2_58(lpMPC_rd, lpMPC_grad_eq, lpMPC_rd);
lpMPC_LA_DIAG_FORWARDBACKWARDSUB_4(lpMPC_Phi00, lpMPC_rd00, lpMPC_dzaff00);
lpMPC_LA_DIAG_FORWARDBACKWARDSUB_4(lpMPC_Phi01, lpMPC_rd01, lpMPC_dzaff01);
lpMPC_LA_DIAG_FORWARDBACKWARDSUB_4(lpMPC_Phi02, lpMPC_rd02, lpMPC_dzaff02);
lpMPC_LA_DIAG_FORWARDBACKWARDSUB_4(lpMPC_Phi03, lpMPC_rd03, lpMPC_dzaff03);
lpMPC_LA_DIAG_FORWARDBACKWARDSUB_4(lpMPC_Phi04, lpMPC_rd04, lpMPC_dzaff04);
lpMPC_LA_DIAG_FORWARDBACKWARDSUB_4(lpMPC_Phi05, lpMPC_rd05, lpMPC_dzaff05);
lpMPC_LA_DIAG_FORWARDBACKWARDSUB_4(lpMPC_Phi06, lpMPC_rd06, lpMPC_dzaff06);
lpMPC_LA_DIAG_FORWARDBACKWARDSUB_4(lpMPC_Phi07, lpMPC_rd07, lpMPC_dzaff07);
lpMPC_LA_DIAG_FORWARDBACKWARDSUB_4(lpMPC_Phi08, lpMPC_rd08, lpMPC_dzaff08);
lpMPC_LA_DIAG_FORWARDBACKWARDSUB_4(lpMPC_Phi09, lpMPC_rd09, lpMPC_dzaff09);
lpMPC_LA_DIAG_FORWARDBACKWARDSUB_4(lpMPC_Phi10, lpMPC_rd10, lpMPC_dzaff10);
lpMPC_LA_DIAG_FORWARDBACKWARDSUB_4(lpMPC_Phi11, lpMPC_rd11, lpMPC_dzaff11);
lpMPC_LA_DIAG_FORWARDBACKWARDSUB_4(lpMPC_Phi12, lpMPC_rd12, lpMPC_dzaff12);
lpMPC_LA_DIAG_FORWARDBACKWARDSUB_4(lpMPC_Phi13, lpMPC_rd13, lpMPC_dzaff13);
lpMPC_LA_DIAG_FORWARDBACKWARDSUB_2(lpMPC_Phi14, lpMPC_rd14, lpMPC_dzaff14);
lpMPC_LA_VSUB_INDEXED_4(lpMPC_dzaff00, lpMPC_lbIdx00, lpMPC_rilb00, lpMPC_dslbaff00);
lpMPC_LA_VSUB3_4(lpMPC_llbbyslb00, lpMPC_dslbaff00, lpMPC_llb00, lpMPC_dllbaff00);
lpMPC_LA_VSUB2_INDEXED_4(lpMPC_riub00, lpMPC_dzaff00, lpMPC_ubIdx00, lpMPC_dsubaff00);
lpMPC_LA_VSUB3_4(lpMPC_lubbysub00, lpMPC_dsubaff00, lpMPC_lub00, lpMPC_dlubaff00);
lpMPC_LA_VSUB_INDEXED_4(lpMPC_dzaff01, lpMPC_lbIdx01, lpMPC_rilb01, lpMPC_dslbaff01);
lpMPC_LA_VSUB3_4(lpMPC_llbbyslb01, lpMPC_dslbaff01, lpMPC_llb01, lpMPC_dllbaff01);
lpMPC_LA_VSUB2_INDEXED_4(lpMPC_riub01, lpMPC_dzaff01, lpMPC_ubIdx01, lpMPC_dsubaff01);
lpMPC_LA_VSUB3_4(lpMPC_lubbysub01, lpMPC_dsubaff01, lpMPC_lub01, lpMPC_dlubaff01);
lpMPC_LA_VSUB_INDEXED_4(lpMPC_dzaff02, lpMPC_lbIdx02, lpMPC_rilb02, lpMPC_dslbaff02);
lpMPC_LA_VSUB3_4(lpMPC_llbbyslb02, lpMPC_dslbaff02, lpMPC_llb02, lpMPC_dllbaff02);
lpMPC_LA_VSUB2_INDEXED_4(lpMPC_riub02, lpMPC_dzaff02, lpMPC_ubIdx02, lpMPC_dsubaff02);
lpMPC_LA_VSUB3_4(lpMPC_lubbysub02, lpMPC_dsubaff02, lpMPC_lub02, lpMPC_dlubaff02);
lpMPC_LA_VSUB_INDEXED_4(lpMPC_dzaff03, lpMPC_lbIdx03, lpMPC_rilb03, lpMPC_dslbaff03);
lpMPC_LA_VSUB3_4(lpMPC_llbbyslb03, lpMPC_dslbaff03, lpMPC_llb03, lpMPC_dllbaff03);
lpMPC_LA_VSUB2_INDEXED_4(lpMPC_riub03, lpMPC_dzaff03, lpMPC_ubIdx03, lpMPC_dsubaff03);
lpMPC_LA_VSUB3_4(lpMPC_lubbysub03, lpMPC_dsubaff03, lpMPC_lub03, lpMPC_dlubaff03);
lpMPC_LA_VSUB_INDEXED_4(lpMPC_dzaff04, lpMPC_lbIdx04, lpMPC_rilb04, lpMPC_dslbaff04);
lpMPC_LA_VSUB3_4(lpMPC_llbbyslb04, lpMPC_dslbaff04, lpMPC_llb04, lpMPC_dllbaff04);
lpMPC_LA_VSUB2_INDEXED_4(lpMPC_riub04, lpMPC_dzaff04, lpMPC_ubIdx04, lpMPC_dsubaff04);
lpMPC_LA_VSUB3_4(lpMPC_lubbysub04, lpMPC_dsubaff04, lpMPC_lub04, lpMPC_dlubaff04);
lpMPC_LA_VSUB_INDEXED_4(lpMPC_dzaff05, lpMPC_lbIdx05, lpMPC_rilb05, lpMPC_dslbaff05);
lpMPC_LA_VSUB3_4(lpMPC_llbbyslb05, lpMPC_dslbaff05, lpMPC_llb05, lpMPC_dllbaff05);
lpMPC_LA_VSUB2_INDEXED_4(lpMPC_riub05, lpMPC_dzaff05, lpMPC_ubIdx05, lpMPC_dsubaff05);
lpMPC_LA_VSUB3_4(lpMPC_lubbysub05, lpMPC_dsubaff05, lpMPC_lub05, lpMPC_dlubaff05);
lpMPC_LA_VSUB_INDEXED_4(lpMPC_dzaff06, lpMPC_lbIdx06, lpMPC_rilb06, lpMPC_dslbaff06);
lpMPC_LA_VSUB3_4(lpMPC_llbbyslb06, lpMPC_dslbaff06, lpMPC_llb06, lpMPC_dllbaff06);
lpMPC_LA_VSUB2_INDEXED_4(lpMPC_riub06, lpMPC_dzaff06, lpMPC_ubIdx06, lpMPC_dsubaff06);
lpMPC_LA_VSUB3_4(lpMPC_lubbysub06, lpMPC_dsubaff06, lpMPC_lub06, lpMPC_dlubaff06);
lpMPC_LA_VSUB_INDEXED_4(lpMPC_dzaff07, lpMPC_lbIdx07, lpMPC_rilb07, lpMPC_dslbaff07);
lpMPC_LA_VSUB3_4(lpMPC_llbbyslb07, lpMPC_dslbaff07, lpMPC_llb07, lpMPC_dllbaff07);
lpMPC_LA_VSUB2_INDEXED_4(lpMPC_riub07, lpMPC_dzaff07, lpMPC_ubIdx07, lpMPC_dsubaff07);
lpMPC_LA_VSUB3_4(lpMPC_lubbysub07, lpMPC_dsubaff07, lpMPC_lub07, lpMPC_dlubaff07);
lpMPC_LA_VSUB_INDEXED_4(lpMPC_dzaff08, lpMPC_lbIdx08, lpMPC_rilb08, lpMPC_dslbaff08);
lpMPC_LA_VSUB3_4(lpMPC_llbbyslb08, lpMPC_dslbaff08, lpMPC_llb08, lpMPC_dllbaff08);
lpMPC_LA_VSUB2_INDEXED_4(lpMPC_riub08, lpMPC_dzaff08, lpMPC_ubIdx08, lpMPC_dsubaff08);
lpMPC_LA_VSUB3_4(lpMPC_lubbysub08, lpMPC_dsubaff08, lpMPC_lub08, lpMPC_dlubaff08);
lpMPC_LA_VSUB_INDEXED_4(lpMPC_dzaff09, lpMPC_lbIdx09, lpMPC_rilb09, lpMPC_dslbaff09);
lpMPC_LA_VSUB3_4(lpMPC_llbbyslb09, lpMPC_dslbaff09, lpMPC_llb09, lpMPC_dllbaff09);
lpMPC_LA_VSUB2_INDEXED_4(lpMPC_riub09, lpMPC_dzaff09, lpMPC_ubIdx09, lpMPC_dsubaff09);
lpMPC_LA_VSUB3_4(lpMPC_lubbysub09, lpMPC_dsubaff09, lpMPC_lub09, lpMPC_dlubaff09);
lpMPC_LA_VSUB_INDEXED_4(lpMPC_dzaff10, lpMPC_lbIdx10, lpMPC_rilb10, lpMPC_dslbaff10);
lpMPC_LA_VSUB3_4(lpMPC_llbbyslb10, lpMPC_dslbaff10, lpMPC_llb10, lpMPC_dllbaff10);
lpMPC_LA_VSUB2_INDEXED_4(lpMPC_riub10, lpMPC_dzaff10, lpMPC_ubIdx10, lpMPC_dsubaff10);
lpMPC_LA_VSUB3_4(lpMPC_lubbysub10, lpMPC_dsubaff10, lpMPC_lub10, lpMPC_dlubaff10);
lpMPC_LA_VSUB_INDEXED_4(lpMPC_dzaff11, lpMPC_lbIdx11, lpMPC_rilb11, lpMPC_dslbaff11);
lpMPC_LA_VSUB3_4(lpMPC_llbbyslb11, lpMPC_dslbaff11, lpMPC_llb11, lpMPC_dllbaff11);
lpMPC_LA_VSUB2_INDEXED_4(lpMPC_riub11, lpMPC_dzaff11, lpMPC_ubIdx11, lpMPC_dsubaff11);
lpMPC_LA_VSUB3_4(lpMPC_lubbysub11, lpMPC_dsubaff11, lpMPC_lub11, lpMPC_dlubaff11);
lpMPC_LA_VSUB_INDEXED_4(lpMPC_dzaff12, lpMPC_lbIdx12, lpMPC_rilb12, lpMPC_dslbaff12);
lpMPC_LA_VSUB3_4(lpMPC_llbbyslb12, lpMPC_dslbaff12, lpMPC_llb12, lpMPC_dllbaff12);
lpMPC_LA_VSUB2_INDEXED_4(lpMPC_riub12, lpMPC_dzaff12, lpMPC_ubIdx12, lpMPC_dsubaff12);
lpMPC_LA_VSUB3_4(lpMPC_lubbysub12, lpMPC_dsubaff12, lpMPC_lub12, lpMPC_dlubaff12);
lpMPC_LA_VSUB_INDEXED_4(lpMPC_dzaff13, lpMPC_lbIdx13, lpMPC_rilb13, lpMPC_dslbaff13);
lpMPC_LA_VSUB3_4(lpMPC_llbbyslb13, lpMPC_dslbaff13, lpMPC_llb13, lpMPC_dllbaff13);
lpMPC_LA_VSUB2_INDEXED_4(lpMPC_riub13, lpMPC_dzaff13, lpMPC_ubIdx13, lpMPC_dsubaff13);
lpMPC_LA_VSUB3_4(lpMPC_lubbysub13, lpMPC_dsubaff13, lpMPC_lub13, lpMPC_dlubaff13);
lpMPC_LA_VSUB_INDEXED_2(lpMPC_dzaff14, lpMPC_lbIdx14, lpMPC_rilb14, lpMPC_dslbaff14);
lpMPC_LA_VSUB3_2(lpMPC_llbbyslb14, lpMPC_dslbaff14, lpMPC_llb14, lpMPC_dllbaff14);
lpMPC_LA_VSUB2_INDEXED_2(lpMPC_riub14, lpMPC_dzaff14, lpMPC_ubIdx14, lpMPC_dsubaff14);
lpMPC_LA_VSUB3_2(lpMPC_lubbysub14, lpMPC_dsubaff14, lpMPC_lub14, lpMPC_dlubaff14);
info->lsit_aff = lpMPC_LINESEARCH_BACKTRACKING_AFFINE(lpMPC_l, lpMPC_s, lpMPC_dl_aff, lpMPC_ds_aff, &info->step_aff, &info->mu_aff);
if( info->lsit_aff == lpMPC_NOPROGRESS ){
exitcode = lpMPC_NOPROGRESS; break;
}
sigma_3rdroot = info->mu_aff / info->mu;
info->sigma = sigma_3rdroot*sigma_3rdroot*sigma_3rdroot;
musigma = info->mu * info->sigma;
lpMPC_LA_VSUB5_116(lpMPC_ds_aff, lpMPC_dl_aff, musigma, lpMPC_ccrhs);
lpMPC_LA_VSUB6_INDEXED_4_4_4(lpMPC_ccrhsub00, lpMPC_sub00, lpMPC_ubIdx00, lpMPC_ccrhsl00, lpMPC_slb00, lpMPC_lbIdx00, lpMPC_rd00);
lpMPC_LA_VSUB6_INDEXED_4_4_4(lpMPC_ccrhsub01, lpMPC_sub01, lpMPC_ubIdx01, lpMPC_ccrhsl01, lpMPC_slb01, lpMPC_lbIdx01, lpMPC_rd01);
lpMPC_LA_DIAG_FORWARDSUB_4(lpMPC_Phi00, lpMPC_rd00, lpMPC_Lbyrd00);
lpMPC_LA_DIAG_FORWARDSUB_4(lpMPC_Phi01, lpMPC_rd01, lpMPC_Lbyrd01);
lpMPC_LA_DIAGZERO_MVM_2(lpMPC_W00, lpMPC_Lbyrd00, lpMPC_beta00);
lpMPC_LA_DENSE_FORWARDSUB_2(lpMPC_Ld00, lpMPC_beta00, lpMPC_yy00);
lpMPC_LA_DENSE_DIAGZERO_2MVMADD_2_4_4(lpMPC_V00, lpMPC_Lbyrd00, lpMPC_W01, lpMPC_Lbyrd01, lpMPC_beta01);
lpMPC_LA_DENSE_MVMSUB1_2_2(lpMPC_Lsd01, lpMPC_yy00, lpMPC_beta01, lpMPC_bmy01);
lpMPC_LA_DENSE_FORWARDSUB_2(lpMPC_Ld01, lpMPC_bmy01, lpMPC_yy01);
lpMPC_LA_VSUB6_INDEXED_4_4_4(lpMPC_ccrhsub02, lpMPC_sub02, lpMPC_ubIdx02, lpMPC_ccrhsl02, lpMPC_slb02, lpMPC_lbIdx02, lpMPC_rd02);
lpMPC_LA_DIAG_FORWARDSUB_4(lpMPC_Phi02, lpMPC_rd02, lpMPC_Lbyrd02);
lpMPC_LA_DENSE_DIAGZERO_2MVMADD_2_4_4(lpMPC_V01, lpMPC_Lbyrd01, lpMPC_W02, lpMPC_Lbyrd02, lpMPC_beta02);
lpMPC_LA_DENSE_MVMSUB1_2_2(lpMPC_Lsd02, lpMPC_yy01, lpMPC_beta02, lpMPC_bmy02);
lpMPC_LA_DENSE_FORWARDSUB_2(lpMPC_Ld02, lpMPC_bmy02, lpMPC_yy02);
lpMPC_LA_VSUB6_INDEXED_4_4_4(lpMPC_ccrhsub03, lpMPC_sub03, lpMPC_ubIdx03, lpMPC_ccrhsl03, lpMPC_slb03, lpMPC_lbIdx03, lpMPC_rd03);
lpMPC_LA_DIAG_FORWARDSUB_4(lpMPC_Phi03, lpMPC_rd03, lpMPC_Lbyrd03);
lpMPC_LA_DENSE_DIAGZERO_2MVMADD_2_4_4(lpMPC_V02, lpMPC_Lbyrd02, lpMPC_W03, lpMPC_Lbyrd03, lpMPC_beta03);
lpMPC_LA_DENSE_MVMSUB1_2_2(lpMPC_Lsd03, lpMPC_yy02, lpMPC_beta03, lpMPC_bmy03);
lpMPC_LA_DENSE_FORWARDSUB_2(lpMPC_Ld03, lpMPC_bmy03, lpMPC_yy03);
lpMPC_LA_VSUB6_INDEXED_4_4_4(lpMPC_ccrhsub04, lpMPC_sub04, lpMPC_ubIdx04, lpMPC_ccrhsl04, lpMPC_slb04, lpMPC_lbIdx04, lpMPC_rd04);
lpMPC_LA_DIAG_FORWARDSUB_4(lpMPC_Phi04, lpMPC_rd04, lpMPC_Lbyrd04);
lpMPC_LA_DENSE_DIAGZERO_2MVMADD_2_4_4(lpMPC_V03, lpMPC_Lbyrd03, lpMPC_W04, lpMPC_Lbyrd04, lpMPC_beta04);
lpMPC_LA_DENSE_MVMSUB1_2_2(lpMPC_Lsd04, lpMPC_yy03, lpMPC_beta04, lpMPC_bmy04);
lpMPC_LA_DENSE_FORWARDSUB_2(lpMPC_Ld04, lpMPC_bmy04, lpMPC_yy04);
lpMPC_LA_VSUB6_INDEXED_4_4_4(lpMPC_ccrhsub05, lpMPC_sub05, lpMPC_ubIdx05, lpMPC_ccrhsl05, lpMPC_slb05, lpMPC_lbIdx05, lpMPC_rd05);
lpMPC_LA_DIAG_FORWARDSUB_4(lpMPC_Phi05, lpMPC_rd05, lpMPC_Lbyrd05);
lpMPC_LA_DENSE_DIAGZERO_2MVMADD_2_4_4(lpMPC_V04, lpMPC_Lbyrd04, lpMPC_W05, lpMPC_Lbyrd05, lpMPC_beta05);
lpMPC_LA_DENSE_MVMSUB1_2_2(lpMPC_Lsd05, lpMPC_yy04, lpMPC_beta05, lpMPC_bmy05);
lpMPC_LA_DENSE_FORWARDSUB_2(lpMPC_Ld05, lpMPC_bmy05, lpMPC_yy05);
lpMPC_LA_VSUB6_INDEXED_4_4_4(lpMPC_ccrhsub06, lpMPC_sub06, lpMPC_ubIdx06, lpMPC_ccrhsl06, lpMPC_slb06, lpMPC_lbIdx06, lpMPC_rd06);
lpMPC_LA_DIAG_FORWARDSUB_4(lpMPC_Phi06, lpMPC_rd06, lpMPC_Lbyrd06);
lpMPC_LA_DENSE_DIAGZERO_2MVMADD_2_4_4(lpMPC_V05, lpMPC_Lbyrd05, lpMPC_W06, lpMPC_Lbyrd06, lpMPC_beta06);
lpMPC_LA_DENSE_MVMSUB1_2_2(lpMPC_Lsd06, lpMPC_yy05, lpMPC_beta06, lpMPC_bmy06);
lpMPC_LA_DENSE_FORWARDSUB_2(lpMPC_Ld06, lpMPC_bmy06, lpMPC_yy06);
lpMPC_LA_VSUB6_INDEXED_4_4_4(lpMPC_ccrhsub07, lpMPC_sub07, lpMPC_ubIdx07, lpMPC_ccrhsl07, lpMPC_slb07, lpMPC_lbIdx07, lpMPC_rd07);
lpMPC_LA_DIAG_FORWARDSUB_4(lpMPC_Phi07, lpMPC_rd07, lpMPC_Lbyrd07);
lpMPC_LA_DENSE_DIAGZERO_2MVMADD_2_4_4(lpMPC_V06, lpMPC_Lbyrd06, lpMPC_W07, lpMPC_Lbyrd07, lpMPC_beta07);
lpMPC_LA_DENSE_MVMSUB1_2_2(lpMPC_Lsd07, lpMPC_yy06, lpMPC_beta07, lpMPC_bmy07);
lpMPC_LA_DENSE_FORWARDSUB_2(lpMPC_Ld07, lpMPC_bmy07, lpMPC_yy07);
lpMPC_LA_VSUB6_INDEXED_4_4_4(lpMPC_ccrhsub08, lpMPC_sub08, lpMPC_ubIdx08, lpMPC_ccrhsl08, lpMPC_slb08, lpMPC_lbIdx08, lpMPC_rd08);
lpMPC_LA_DIAG_FORWARDSUB_4(lpMPC_Phi08, lpMPC_rd08, lpMPC_Lbyrd08);
lpMPC_LA_DENSE_DIAGZERO_2MVMADD_2_4_4(lpMPC_V07, lpMPC_Lbyrd07, lpMPC_W08, lpMPC_Lbyrd08, lpMPC_beta08);
lpMPC_LA_DENSE_MVMSUB1_2_2(lpMPC_Lsd08, lpMPC_yy07, lpMPC_beta08, lpMPC_bmy08);
lpMPC_LA_DENSE_FORWARDSUB_2(lpMPC_Ld08, lpMPC_bmy08, lpMPC_yy08);
lpMPC_LA_VSUB6_INDEXED_4_4_4(lpMPC_ccrhsub09, lpMPC_sub09, lpMPC_ubIdx09, lpMPC_ccrhsl09, lpMPC_slb09, lpMPC_lbIdx09, lpMPC_rd09);
lpMPC_LA_DIAG_FORWARDSUB_4(lpMPC_Phi09, lpMPC_rd09, lpMPC_Lbyrd09);
lpMPC_LA_DENSE_DIAGZERO_2MVMADD_2_4_4(lpMPC_V08, lpMPC_Lbyrd08, lpMPC_W09, lpMPC_Lbyrd09, lpMPC_beta09);
lpMPC_LA_DENSE_MVMSUB1_2_2(lpMPC_Lsd09, lpMPC_yy08, lpMPC_beta09, lpMPC_bmy09);
lpMPC_LA_DENSE_FORWARDSUB_2(lpMPC_Ld09, lpMPC_bmy09, lpMPC_yy09);
lpMPC_LA_VSUB6_INDEXED_4_4_4(lpMPC_ccrhsub10, lpMPC_sub10, lpMPC_ubIdx10, lpMPC_ccrhsl10, lpMPC_slb10, lpMPC_lbIdx10, lpMPC_rd10);
lpMPC_LA_DIAG_FORWARDSUB_4(lpMPC_Phi10, lpMPC_rd10, lpMPC_Lbyrd10);
lpMPC_LA_DENSE_DIAGZERO_2MVMADD_2_4_4(lpMPC_V09, lpMPC_Lbyrd09, lpMPC_W10, lpMPC_Lbyrd10, lpMPC_beta10);
lpMPC_LA_DENSE_MVMSUB1_2_2(lpMPC_Lsd10, lpMPC_yy09, lpMPC_beta10, lpMPC_bmy10);
lpMPC_LA_DENSE_FORWARDSUB_2(lpMPC_Ld10, lpMPC_bmy10, lpMPC_yy10);
lpMPC_LA_VSUB6_INDEXED_4_4_4(lpMPC_ccrhsub11, lpMPC_sub11, lpMPC_ubIdx11, lpMPC_ccrhsl11, lpMPC_slb11, lpMPC_lbIdx11, lpMPC_rd11);
lpMPC_LA_DIAG_FORWARDSUB_4(lpMPC_Phi11, lpMPC_rd11, lpMPC_Lbyrd11);
lpMPC_LA_DENSE_DIAGZERO_2MVMADD_2_4_4(lpMPC_V10, lpMPC_Lbyrd10, lpMPC_W11, lpMPC_Lbyrd11, lpMPC_beta11);
lpMPC_LA_DENSE_MVMSUB1_2_2(lpMPC_Lsd11, lpMPC_yy10, lpMPC_beta11, lpMPC_bmy11);
lpMPC_LA_DENSE_FORWARDSUB_2(lpMPC_Ld11, lpMPC_bmy11, lpMPC_yy11);
lpMPC_LA_VSUB6_INDEXED_4_4_4(lpMPC_ccrhsub12, lpMPC_sub12, lpMPC_ubIdx12, lpMPC_ccrhsl12, lpMPC_slb12, lpMPC_lbIdx12, lpMPC_rd12);
lpMPC_LA_DIAG_FORWARDSUB_4(lpMPC_Phi12, lpMPC_rd12, lpMPC_Lbyrd12);
lpMPC_LA_DENSE_DIAGZERO_2MVMADD_2_4_4(lpMPC_V11, lpMPC_Lbyrd11, lpMPC_W12, lpMPC_Lbyrd12, lpMPC_beta12);
lpMPC_LA_DENSE_MVMSUB1_2_2(lpMPC_Lsd12, lpMPC_yy11, lpMPC_beta12, lpMPC_bmy12);
lpMPC_LA_DENSE_FORWARDSUB_2(lpMPC_Ld12, lpMPC_bmy12, lpMPC_yy12);
lpMPC_LA_VSUB6_INDEXED_4_4_4(lpMPC_ccrhsub13, lpMPC_sub13, lpMPC_ubIdx13, lpMPC_ccrhsl13, lpMPC_slb13, lpMPC_lbIdx13, lpMPC_rd13);
lpMPC_LA_DIAG_FORWARDSUB_4(lpMPC_Phi13, lpMPC_rd13, lpMPC_Lbyrd13);
lpMPC_LA_DENSE_DIAGZERO_2MVMADD_2_4_4(lpMPC_V12, lpMPC_Lbyrd12, lpMPC_W13, lpMPC_Lbyrd13, lpMPC_beta13);
lpMPC_LA_DENSE_MVMSUB1_2_2(lpMPC_Lsd13, lpMPC_yy12, lpMPC_beta13, lpMPC_bmy13);
lpMPC_LA_DENSE_FORWARDSUB_2(lpMPC_Ld13, lpMPC_bmy13, lpMPC_yy13);
lpMPC_LA_VSUB6_INDEXED_2_2_2(lpMPC_ccrhsub14, lpMPC_sub14, lpMPC_ubIdx14, lpMPC_ccrhsl14, lpMPC_slb14, lpMPC_lbIdx14, lpMPC_rd14);
lpMPC_LA_DIAG_FORWARDSUB_2(lpMPC_Phi14, lpMPC_rd14, lpMPC_Lbyrd14);
lpMPC_LA_DENSE_DIAGZERO_2MVMADD_2_4_2(lpMPC_V13, lpMPC_Lbyrd13, lpMPC_W14, lpMPC_Lbyrd14, lpMPC_beta14);
lpMPC_LA_DENSE_MVMSUB1_2_2(lpMPC_Lsd14, lpMPC_yy13, lpMPC_beta14, lpMPC_bmy14);
lpMPC_LA_DENSE_FORWARDSUB_2(lpMPC_Ld14, lpMPC_bmy14, lpMPC_yy14);
lpMPC_LA_DENSE_BACKWARDSUB_2(lpMPC_Ld14, lpMPC_yy14, lpMPC_dvcc14);
lpMPC_LA_DENSE_MTVMSUB_2_2(lpMPC_Lsd14, lpMPC_dvcc14, lpMPC_yy13, lpMPC_bmy13);
lpMPC_LA_DENSE_BACKWARDSUB_2(lpMPC_Ld13, lpMPC_bmy13, lpMPC_dvcc13);
lpMPC_LA_DENSE_MTVMSUB_2_2(lpMPC_Lsd13, lpMPC_dvcc13, lpMPC_yy12, lpMPC_bmy12);
lpMPC_LA_DENSE_BACKWARDSUB_2(lpMPC_Ld12, lpMPC_bmy12, lpMPC_dvcc12);
lpMPC_LA_DENSE_MTVMSUB_2_2(lpMPC_Lsd12, lpMPC_dvcc12, lpMPC_yy11, lpMPC_bmy11);
lpMPC_LA_DENSE_BACKWARDSUB_2(lpMPC_Ld11, lpMPC_bmy11, lpMPC_dvcc11);
lpMPC_LA_DENSE_MTVMSUB_2_2(lpMPC_Lsd11, lpMPC_dvcc11, lpMPC_yy10, lpMPC_bmy10);
lpMPC_LA_DENSE_BACKWARDSUB_2(lpMPC_Ld10, lpMPC_bmy10, lpMPC_dvcc10);
lpMPC_LA_DENSE_MTVMSUB_2_2(lpMPC_Lsd10, lpMPC_dvcc10, lpMPC_yy09, lpMPC_bmy09);
lpMPC_LA_DENSE_BACKWARDSUB_2(lpMPC_Ld09, lpMPC_bmy09, lpMPC_dvcc09);
lpMPC_LA_DENSE_MTVMSUB_2_2(lpMPC_Lsd09, lpMPC_dvcc09, lpMPC_yy08, lpMPC_bmy08);
lpMPC_LA_DENSE_BACKWARDSUB_2(lpMPC_Ld08, lpMPC_bmy08, lpMPC_dvcc08);
lpMPC_LA_DENSE_MTVMSUB_2_2(lpMPC_Lsd08, lpMPC_dvcc08, lpMPC_yy07, lpMPC_bmy07);
lpMPC_LA_DENSE_BACKWARDSUB_2(lpMPC_Ld07, lpMPC_bmy07, lpMPC_dvcc07);
lpMPC_LA_DENSE_MTVMSUB_2_2(lpMPC_Lsd07, lpMPC_dvcc07, lpMPC_yy06, lpMPC_bmy06);
lpMPC_LA_DENSE_BACKWARDSUB_2(lpMPC_Ld06, lpMPC_bmy06, lpMPC_dvcc06);
lpMPC_LA_DENSE_MTVMSUB_2_2(lpMPC_Lsd06, lpMPC_dvcc06, lpMPC_yy05, lpMPC_bmy05);
lpMPC_LA_DENSE_BACKWARDSUB_2(lpMPC_Ld05, lpMPC_bmy05, lpMPC_dvcc05);
lpMPC_LA_DENSE_MTVMSUB_2_2(lpMPC_Lsd05, lpMPC_dvcc05, lpMPC_yy04, lpMPC_bmy04);
lpMPC_LA_DENSE_BACKWARDSUB_2(lpMPC_Ld04, lpMPC_bmy04, lpMPC_dvcc04);
lpMPC_LA_DENSE_MTVMSUB_2_2(lpMPC_Lsd04, lpMPC_dvcc04, lpMPC_yy03, lpMPC_bmy03);
lpMPC_LA_DENSE_BACKWARDSUB_2(lpMPC_Ld03, lpMPC_bmy03, lpMPC_dvcc03);
lpMPC_LA_DENSE_MTVMSUB_2_2(lpMPC_Lsd03, lpMPC_dvcc03, lpMPC_yy02, lpMPC_bmy02);
lpMPC_LA_DENSE_BACKWARDSUB_2(lpMPC_Ld02, lpMPC_bmy02, lpMPC_dvcc02);
lpMPC_LA_DENSE_MTVMSUB_2_2(lpMPC_Lsd02, lpMPC_dvcc02, lpMPC_yy01, lpMPC_bmy01);
lpMPC_LA_DENSE_BACKWARDSUB_2(lpMPC_Ld01, lpMPC_bmy01, lpMPC_dvcc01);
lpMPC_LA_DENSE_MTVMSUB_2_2(lpMPC_Lsd01, lpMPC_dvcc01, lpMPC_yy00, lpMPC_bmy00);
lpMPC_LA_DENSE_BACKWARDSUB_2(lpMPC_Ld00, lpMPC_bmy00, lpMPC_dvcc00);
lpMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C1, lpMPC_dvcc01, lpMPC_D00, lpMPC_dvcc00, lpMPC_grad_eq00);
lpMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C2, lpMPC_dvcc02, lpMPC_D01, lpMPC_dvcc01, lpMPC_grad_eq01);
lpMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C3, lpMPC_dvcc03, lpMPC_D01, lpMPC_dvcc02, lpMPC_grad_eq02);
lpMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C4, lpMPC_dvcc04, lpMPC_D01, lpMPC_dvcc03, lpMPC_grad_eq03);
lpMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C5, lpMPC_dvcc05, lpMPC_D01, lpMPC_dvcc04, lpMPC_grad_eq04);
lpMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C6, lpMPC_dvcc06, lpMPC_D01, lpMPC_dvcc05, lpMPC_grad_eq05);
lpMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C7, lpMPC_dvcc07, lpMPC_D01, lpMPC_dvcc06, lpMPC_grad_eq06);
lpMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C8, lpMPC_dvcc08, lpMPC_D01, lpMPC_dvcc07, lpMPC_grad_eq07);
lpMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C9, lpMPC_dvcc09, lpMPC_D01, lpMPC_dvcc08, lpMPC_grad_eq08);
lpMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C10, lpMPC_dvcc10, lpMPC_D01, lpMPC_dvcc09, lpMPC_grad_eq09);
lpMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C11, lpMPC_dvcc11, lpMPC_D01, lpMPC_dvcc10, lpMPC_grad_eq10);
lpMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C12, lpMPC_dvcc12, lpMPC_D01, lpMPC_dvcc11, lpMPC_grad_eq11);
lpMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C13, lpMPC_dvcc13, lpMPC_D01, lpMPC_dvcc12, lpMPC_grad_eq12);
lpMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C14, lpMPC_dvcc14, lpMPC_D01, lpMPC_dvcc13, lpMPC_grad_eq13);
lpMPC_LA_DIAGZERO_MTVM_2_2(lpMPC_D14, lpMPC_dvcc14, lpMPC_grad_eq14);
lpMPC_LA_VSUB_58(lpMPC_rd, lpMPC_grad_eq, lpMPC_rd);
lpMPC_LA_DIAG_FORWARDBACKWARDSUB_4(lpMPC_Phi00, lpMPC_rd00, lpMPC_dzcc00);
lpMPC_LA_DIAG_FORWARDBACKWARDSUB_4(lpMPC_Phi01, lpMPC_rd01, lpMPC_dzcc01);
lpMPC_LA_DIAG_FORWARDBACKWARDSUB_4(lpMPC_Phi02, lpMPC_rd02, lpMPC_dzcc02);
lpMPC_LA_DIAG_FORWARDBACKWARDSUB_4(lpMPC_Phi03, lpMPC_rd03, lpMPC_dzcc03);
lpMPC_LA_DIAG_FORWARDBACKWARDSUB_4(lpMPC_Phi04, lpMPC_rd04, lpMPC_dzcc04);
lpMPC_LA_DIAG_FORWARDBACKWARDSUB_4(lpMPC_Phi05, lpMPC_rd05, lpMPC_dzcc05);
lpMPC_LA_DIAG_FORWARDBACKWARDSUB_4(lpMPC_Phi06, lpMPC_rd06, lpMPC_dzcc06);
lpMPC_LA_DIAG_FORWARDBACKWARDSUB_4(lpMPC_Phi07, lpMPC_rd07, lpMPC_dzcc07);
lpMPC_LA_DIAG_FORWARDBACKWARDSUB_4(lpMPC_Phi08, lpMPC_rd08, lpMPC_dzcc08);
lpMPC_LA_DIAG_FORWARDBACKWARDSUB_4(lpMPC_Phi09, lpMPC_rd09, lpMPC_dzcc09);
lpMPC_LA_DIAG_FORWARDBACKWARDSUB_4(lpMPC_Phi10, lpMPC_rd10, lpMPC_dzcc10);
lpMPC_LA_DIAG_FORWARDBACKWARDSUB_4(lpMPC_Phi11, lpMPC_rd11, lpMPC_dzcc11);
lpMPC_LA_DIAG_FORWARDBACKWARDSUB_4(lpMPC_Phi12, lpMPC_rd12, lpMPC_dzcc12);
lpMPC_LA_DIAG_FORWARDBACKWARDSUB_4(lpMPC_Phi13, lpMPC_rd13, lpMPC_dzcc13);
lpMPC_LA_DIAG_FORWARDBACKWARDSUB_2(lpMPC_Phi14, lpMPC_rd14, lpMPC_dzcc14);
lpMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_4(lpMPC_ccrhsl00, lpMPC_slb00, lpMPC_llbbyslb00, lpMPC_dzcc00, lpMPC_lbIdx00, lpMPC_dllbcc00);
lpMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_4(lpMPC_ccrhsub00, lpMPC_sub00, lpMPC_lubbysub00, lpMPC_dzcc00, lpMPC_ubIdx00, lpMPC_dlubcc00);
lpMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_4(lpMPC_ccrhsl01, lpMPC_slb01, lpMPC_llbbyslb01, lpMPC_dzcc01, lpMPC_lbIdx01, lpMPC_dllbcc01);
lpMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_4(lpMPC_ccrhsub01, lpMPC_sub01, lpMPC_lubbysub01, lpMPC_dzcc01, lpMPC_ubIdx01, lpMPC_dlubcc01);
lpMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_4(lpMPC_ccrhsl02, lpMPC_slb02, lpMPC_llbbyslb02, lpMPC_dzcc02, lpMPC_lbIdx02, lpMPC_dllbcc02);
lpMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_4(lpMPC_ccrhsub02, lpMPC_sub02, lpMPC_lubbysub02, lpMPC_dzcc02, lpMPC_ubIdx02, lpMPC_dlubcc02);
lpMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_4(lpMPC_ccrhsl03, lpMPC_slb03, lpMPC_llbbyslb03, lpMPC_dzcc03, lpMPC_lbIdx03, lpMPC_dllbcc03);
lpMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_4(lpMPC_ccrhsub03, lpMPC_sub03, lpMPC_lubbysub03, lpMPC_dzcc03, lpMPC_ubIdx03, lpMPC_dlubcc03);
lpMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_4(lpMPC_ccrhsl04, lpMPC_slb04, lpMPC_llbbyslb04, lpMPC_dzcc04, lpMPC_lbIdx04, lpMPC_dllbcc04);
lpMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_4(lpMPC_ccrhsub04, lpMPC_sub04, lpMPC_lubbysub04, lpMPC_dzcc04, lpMPC_ubIdx04, lpMPC_dlubcc04);
lpMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_4(lpMPC_ccrhsl05, lpMPC_slb05, lpMPC_llbbyslb05, lpMPC_dzcc05, lpMPC_lbIdx05, lpMPC_dllbcc05);
lpMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_4(lpMPC_ccrhsub05, lpMPC_sub05, lpMPC_lubbysub05, lpMPC_dzcc05, lpMPC_ubIdx05, lpMPC_dlubcc05);
lpMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_4(lpMPC_ccrhsl06, lpMPC_slb06, lpMPC_llbbyslb06, lpMPC_dzcc06, lpMPC_lbIdx06, lpMPC_dllbcc06);
lpMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_4(lpMPC_ccrhsub06, lpMPC_sub06, lpMPC_lubbysub06, lpMPC_dzcc06, lpMPC_ubIdx06, lpMPC_dlubcc06);
lpMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_4(lpMPC_ccrhsl07, lpMPC_slb07, lpMPC_llbbyslb07, lpMPC_dzcc07, lpMPC_lbIdx07, lpMPC_dllbcc07);
lpMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_4(lpMPC_ccrhsub07, lpMPC_sub07, lpMPC_lubbysub07, lpMPC_dzcc07, lpMPC_ubIdx07, lpMPC_dlubcc07);
lpMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_4(lpMPC_ccrhsl08, lpMPC_slb08, lpMPC_llbbyslb08, lpMPC_dzcc08, lpMPC_lbIdx08, lpMPC_dllbcc08);
lpMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_4(lpMPC_ccrhsub08, lpMPC_sub08, lpMPC_lubbysub08, lpMPC_dzcc08, lpMPC_ubIdx08, lpMPC_dlubcc08);
lpMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_4(lpMPC_ccrhsl09, lpMPC_slb09, lpMPC_llbbyslb09, lpMPC_dzcc09, lpMPC_lbIdx09, lpMPC_dllbcc09);
lpMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_4(lpMPC_ccrhsub09, lpMPC_sub09, lpMPC_lubbysub09, lpMPC_dzcc09, lpMPC_ubIdx09, lpMPC_dlubcc09);
lpMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_4(lpMPC_ccrhsl10, lpMPC_slb10, lpMPC_llbbyslb10, lpMPC_dzcc10, lpMPC_lbIdx10, lpMPC_dllbcc10);
lpMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_4(lpMPC_ccrhsub10, lpMPC_sub10, lpMPC_lubbysub10, lpMPC_dzcc10, lpMPC_ubIdx10, lpMPC_dlubcc10);
lpMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_4(lpMPC_ccrhsl11, lpMPC_slb11, lpMPC_llbbyslb11, lpMPC_dzcc11, lpMPC_lbIdx11, lpMPC_dllbcc11);
lpMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_4(lpMPC_ccrhsub11, lpMPC_sub11, lpMPC_lubbysub11, lpMPC_dzcc11, lpMPC_ubIdx11, lpMPC_dlubcc11);
lpMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_4(lpMPC_ccrhsl12, lpMPC_slb12, lpMPC_llbbyslb12, lpMPC_dzcc12, lpMPC_lbIdx12, lpMPC_dllbcc12);
lpMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_4(lpMPC_ccrhsub12, lpMPC_sub12, lpMPC_lubbysub12, lpMPC_dzcc12, lpMPC_ubIdx12, lpMPC_dlubcc12);
lpMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_4(lpMPC_ccrhsl13, lpMPC_slb13, lpMPC_llbbyslb13, lpMPC_dzcc13, lpMPC_lbIdx13, lpMPC_dllbcc13);
lpMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_4(lpMPC_ccrhsub13, lpMPC_sub13, lpMPC_lubbysub13, lpMPC_dzcc13, lpMPC_ubIdx13, lpMPC_dlubcc13);
lpMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_2(lpMPC_ccrhsl14, lpMPC_slb14, lpMPC_llbbyslb14, lpMPC_dzcc14, lpMPC_lbIdx14, lpMPC_dllbcc14);
lpMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_2(lpMPC_ccrhsub14, lpMPC_sub14, lpMPC_lubbysub14, lpMPC_dzcc14, lpMPC_ubIdx14, lpMPC_dlubcc14);
lpMPC_LA_VSUB7_116(lpMPC_l, lpMPC_ccrhs, lpMPC_s, lpMPC_dl_cc, lpMPC_ds_cc);
lpMPC_LA_VADD_58(lpMPC_dz_cc, lpMPC_dz_aff);
lpMPC_LA_VADD_30(lpMPC_dv_cc, lpMPC_dv_aff);
lpMPC_LA_VADD_116(lpMPC_dl_cc, lpMPC_dl_aff);
lpMPC_LA_VADD_116(lpMPC_ds_cc, lpMPC_ds_aff);
info->lsit_cc = lpMPC_LINESEARCH_BACKTRACKING_COMBINED(lpMPC_z, lpMPC_v, lpMPC_l, lpMPC_s, lpMPC_dz_cc, lpMPC_dv_cc, lpMPC_dl_cc, lpMPC_ds_cc, &info->step_cc, &info->mu);
if( info->lsit_cc == lpMPC_NOPROGRESS ){
exitcode = lpMPC_NOPROGRESS; break;
}
info->it++;
}
output->z1[0] = lpMPC_z00[0];
output->z1[1] = lpMPC_z00[1];
output->z1[2] = lpMPC_z00[2];
output->z1[3] = lpMPC_z00[3];
output->z2[0] = lpMPC_z01[0];
output->z2[1] = lpMPC_z01[1];
output->z2[2] = lpMPC_z01[2];
output->z2[3] = lpMPC_z01[3];
output->z3[0] = lpMPC_z02[0];
output->z3[1] = lpMPC_z02[1];
output->z3[2] = lpMPC_z02[2];
output->z3[3] = lpMPC_z02[3];
output->z4[0] = lpMPC_z03[0];
output->z4[1] = lpMPC_z03[1];
output->z4[2] = lpMPC_z03[2];
output->z4[3] = lpMPC_z03[3];
output->z5[0] = lpMPC_z04[0];
output->z5[1] = lpMPC_z04[1];
output->z5[2] = lpMPC_z04[2];
output->z5[3] = lpMPC_z04[3];
output->z6[0] = lpMPC_z05[0];
output->z6[1] = lpMPC_z05[1];
output->z6[2] = lpMPC_z05[2];
output->z6[3] = lpMPC_z05[3];
output->z7[0] = lpMPC_z06[0];
output->z7[1] = lpMPC_z06[1];
output->z7[2] = lpMPC_z06[2];
output->z7[3] = lpMPC_z06[3];
output->z8[0] = lpMPC_z07[0];
output->z8[1] = lpMPC_z07[1];
output->z8[2] = lpMPC_z07[2];
output->z8[3] = lpMPC_z07[3];
output->z9[0] = lpMPC_z08[0];
output->z9[1] = lpMPC_z08[1];
output->z9[2] = lpMPC_z08[2];
output->z9[3] = lpMPC_z08[3];
output->z10[0] = lpMPC_z09[0];
output->z10[1] = lpMPC_z09[1];
output->z10[2] = lpMPC_z09[2];
output->z10[3] = lpMPC_z09[3];
output->z11[0] = lpMPC_z10[0];
output->z11[1] = lpMPC_z10[1];
output->z11[2] = lpMPC_z10[2];
output->z11[3] = lpMPC_z10[3];
output->z12[0] = lpMPC_z11[0];
output->z12[1] = lpMPC_z11[1];
output->z12[2] = lpMPC_z11[2];
output->z12[3] = lpMPC_z11[3];
output->z13[0] = lpMPC_z12[0];
output->z13[1] = lpMPC_z12[1];
output->z13[2] = lpMPC_z12[2];
output->z13[3] = lpMPC_z12[3];
output->z14[0] = lpMPC_z13[0];
output->z14[1] = lpMPC_z13[1];
output->z14[2] = lpMPC_z13[2];
output->z14[3] = lpMPC_z13[3];
output->z15[0] = lpMPC_z14[0];
output->z15[1] = lpMPC_z14[1];

#if lpMPC_SET_TIMING == 1
info->solvetime = lpMPC_toc(&solvertimer);
#if lpMPC_SET_PRINTLEVEL > 0 && lpMPC_SET_TIMING == 1
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
