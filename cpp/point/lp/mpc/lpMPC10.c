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
 * Initializes a vector of length 38 with a value.
 */
void lpMPC_LA_INITIALIZEVECTOR_38(lpMPC_FLOAT* vec, lpMPC_FLOAT value)
{
	int i;
	for( i=0; i<38; i++ )
	{
		vec[i] = value;
	}
}


/*
 * Initializes a vector of length 20 with a value.
 */
void lpMPC_LA_INITIALIZEVECTOR_20(lpMPC_FLOAT* vec, lpMPC_FLOAT value)
{
	int i;
	for( i=0; i<20; i++ )
	{
		vec[i] = value;
	}
}


/*
 * Initializes a vector of length 76 with a value.
 */
void lpMPC_LA_INITIALIZEVECTOR_76(lpMPC_FLOAT* vec, lpMPC_FLOAT value)
{
	int i;
	for( i=0; i<76; i++ )
	{
		vec[i] = value;
	}
}


/* 
 * Calculates a dot product and adds it to a variable: z += x'*y; 
 * This function is for vectors of length 76.
 */
void lpMPC_LA_DOTACC_76(lpMPC_FLOAT *x, lpMPC_FLOAT *y, lpMPC_FLOAT *z)
{
	int i;
	for( i=0; i<76; i++ ){
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
 * of length 38.
 */
void lpMPC_LA_VVADD3_38(lpMPC_FLOAT *u, lpMPC_FLOAT *v, lpMPC_FLOAT *w, lpMPC_FLOAT *z)
{
	int i;
	for( i=0; i<38; i++ ){
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
 * Vector subtraction z = -x - y for vectors of length 38.
 */
void lpMPC_LA_VSUB2_38(lpMPC_FLOAT *x, lpMPC_FLOAT *y, lpMPC_FLOAT *z)
{
	int i;
	for( i=0; i<38; i++){
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
        for( i=0; i<76; i++ ){
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
        if( i == 76 ){
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
    *mu_aff = mymu / (lpMPC_FLOAT)76;
    return lsIt;
}


/*
 * Vector subtraction x = u.*v - a where a is a scalar
*  and x,u,v are vectors of length 76.
 */
void lpMPC_LA_VSUB5_76(lpMPC_FLOAT *u, lpMPC_FLOAT *v, lpMPC_FLOAT a, lpMPC_FLOAT *x)
{
	int i;
	for( i=0; i<76; i++){
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
 * Vector subtraction z = x - y for vectors of length 38.
 */
void lpMPC_LA_VSUB_38(lpMPC_FLOAT *x, lpMPC_FLOAT *y, lpMPC_FLOAT *z)
{
	int i;
	for( i=0; i<38; i++){
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
 * Computes ds = -l.\(r + s.*dl) for vectors of length 76.
 */
void lpMPC_LA_VSUB7_76(lpMPC_FLOAT *l, lpMPC_FLOAT *r, lpMPC_FLOAT *s, lpMPC_FLOAT *dl, lpMPC_FLOAT *ds)
{
	int i;
	for( i=0; i<76; i++){
		ds[i] = -(r[i] + s[i]*dl[i])/l[i];
	}
}


/*
 * Vector addition x = x + y for vectors of length 38.
 */
void lpMPC_LA_VADD_38(lpMPC_FLOAT *x, lpMPC_FLOAT *y)
{
	int i;
	for( i=0; i<38; i++){
		x[i] += y[i];
	}
}


/*
 * Vector addition x = x + y for vectors of length 20.
 */
void lpMPC_LA_VADD_20(lpMPC_FLOAT *x, lpMPC_FLOAT *y)
{
	int i;
	for( i=0; i<20; i++){
		x[i] += y[i];
	}
}


/*
 * Vector addition x = x + y for vectors of length 76.
 */
void lpMPC_LA_VADD_76(lpMPC_FLOAT *x, lpMPC_FLOAT *y)
{
	int i;
	for( i=0; i<76; i++){
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
        for( i=0; i<76; i++ ){
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
        if( i == 76 ){
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
    for( i=0; i<38; i++ ){
        z[i] += a_gamma*dz[i];
    }
    
    /* equality constraint multipliers */
    for( i=0; i<20; i++ ){
        v[i] += a_gamma*dv[i];
    }
    
    /* inequality constraint multipliers & slacks, also update mu */
    *mu = 0;
    for( i=0; i<76; i++ ){
        dltemp = l[i] + a_gamma*dl[i]; l[i] = dltemp;
        dstemp = s[i] + a_gamma*ds[i]; s[i] = dstemp;
        *mu += dltemp*dstemp;
    }
    
    *a = a_gamma;
    *mu /= (lpMPC_FLOAT)76;
    return lsIt;
}




/* VARIABLE DEFINITIONS ------------------------------------------------ */
lpMPC_FLOAT lpMPC_z[38];
lpMPC_FLOAT lpMPC_v[20];
lpMPC_FLOAT lpMPC_dz_aff[38];
lpMPC_FLOAT lpMPC_dv_aff[20];
lpMPC_FLOAT lpMPC_grad_cost[38];
lpMPC_FLOAT lpMPC_grad_eq[38];
lpMPC_FLOAT lpMPC_rd[38];
lpMPC_FLOAT lpMPC_l[76];
lpMPC_FLOAT lpMPC_s[76];
lpMPC_FLOAT lpMPC_lbys[76];
lpMPC_FLOAT lpMPC_dl_aff[76];
lpMPC_FLOAT lpMPC_ds_aff[76];
lpMPC_FLOAT lpMPC_dz_cc[38];
lpMPC_FLOAT lpMPC_dv_cc[20];
lpMPC_FLOAT lpMPC_dl_cc[76];
lpMPC_FLOAT lpMPC_ds_cc[76];
lpMPC_FLOAT lpMPC_ccrhs[76];
lpMPC_FLOAT lpMPC_grad_ineq[38];
lpMPC_FLOAT lpMPC_H0[4] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
lpMPC_FLOAT* lpMPC_z0 = lpMPC_z + 0;
lpMPC_FLOAT* lpMPC_dzaff0 = lpMPC_dz_aff + 0;
lpMPC_FLOAT* lpMPC_dzcc0 = lpMPC_dz_cc + 0;
lpMPC_FLOAT* lpMPC_rd0 = lpMPC_rd + 0;
lpMPC_FLOAT lpMPC_Lbyrd0[4];
lpMPC_FLOAT* lpMPC_grad_cost0 = lpMPC_grad_cost + 0;
lpMPC_FLOAT* lpMPC_grad_eq0 = lpMPC_grad_eq + 0;
lpMPC_FLOAT* lpMPC_grad_ineq0 = lpMPC_grad_ineq + 0;
lpMPC_FLOAT lpMPC_ctv0[4];
lpMPC_FLOAT* lpMPC_v0 = lpMPC_v + 0;
lpMPC_FLOAT lpMPC_re0[2];
lpMPC_FLOAT lpMPC_beta0[2];
lpMPC_FLOAT lpMPC_betacc0[2];
lpMPC_FLOAT* lpMPC_dvaff0 = lpMPC_dv_aff + 0;
lpMPC_FLOAT* lpMPC_dvcc0 = lpMPC_dv_cc + 0;
lpMPC_FLOAT lpMPC_V0[8];
lpMPC_FLOAT lpMPC_Yd0[3];
lpMPC_FLOAT lpMPC_Ld0[3];
lpMPC_FLOAT lpMPC_yy0[2];
lpMPC_FLOAT lpMPC_bmy0[2];
int lpMPC_lbIdx0[4] = {0, 1, 2, 3};
lpMPC_FLOAT* lpMPC_llb0 = lpMPC_l + 0;
lpMPC_FLOAT* lpMPC_slb0 = lpMPC_s + 0;
lpMPC_FLOAT* lpMPC_llbbyslb0 = lpMPC_lbys + 0;
lpMPC_FLOAT lpMPC_rilb0[4];
lpMPC_FLOAT* lpMPC_dllbaff0 = lpMPC_dl_aff + 0;
lpMPC_FLOAT* lpMPC_dslbaff0 = lpMPC_ds_aff + 0;
lpMPC_FLOAT* lpMPC_dllbcc0 = lpMPC_dl_cc + 0;
lpMPC_FLOAT* lpMPC_dslbcc0 = lpMPC_ds_cc + 0;
lpMPC_FLOAT* lpMPC_ccrhsl0 = lpMPC_ccrhs + 0;
int lpMPC_ubIdx0[4] = {0, 1, 2, 3};
lpMPC_FLOAT* lpMPC_lub0 = lpMPC_l + 4;
lpMPC_FLOAT* lpMPC_sub0 = lpMPC_s + 4;
lpMPC_FLOAT* lpMPC_lubbysub0 = lpMPC_lbys + 4;
lpMPC_FLOAT lpMPC_riub0[4];
lpMPC_FLOAT* lpMPC_dlubaff0 = lpMPC_dl_aff + 4;
lpMPC_FLOAT* lpMPC_dsubaff0 = lpMPC_ds_aff + 4;
lpMPC_FLOAT* lpMPC_dlubcc0 = lpMPC_dl_cc + 4;
lpMPC_FLOAT* lpMPC_dsubcc0 = lpMPC_ds_cc + 4;
lpMPC_FLOAT* lpMPC_ccrhsub0 = lpMPC_ccrhs + 4;
lpMPC_FLOAT lpMPC_Phi0[4];
lpMPC_FLOAT lpMPC_D0[4] = {1.0000000000000000E+000, 
1.0000000000000000E+000};
lpMPC_FLOAT lpMPC_W0[4];
lpMPC_FLOAT* lpMPC_z1 = lpMPC_z + 4;
lpMPC_FLOAT* lpMPC_dzaff1 = lpMPC_dz_aff + 4;
lpMPC_FLOAT* lpMPC_dzcc1 = lpMPC_dz_cc + 4;
lpMPC_FLOAT* lpMPC_rd1 = lpMPC_rd + 4;
lpMPC_FLOAT lpMPC_Lbyrd1[4];
lpMPC_FLOAT* lpMPC_grad_cost1 = lpMPC_grad_cost + 4;
lpMPC_FLOAT* lpMPC_grad_eq1 = lpMPC_grad_eq + 4;
lpMPC_FLOAT* lpMPC_grad_ineq1 = lpMPC_grad_ineq + 4;
lpMPC_FLOAT lpMPC_ctv1[4];
lpMPC_FLOAT* lpMPC_v1 = lpMPC_v + 2;
lpMPC_FLOAT lpMPC_re1[2];
lpMPC_FLOAT lpMPC_beta1[2];
lpMPC_FLOAT lpMPC_betacc1[2];
lpMPC_FLOAT* lpMPC_dvaff1 = lpMPC_dv_aff + 2;
lpMPC_FLOAT* lpMPC_dvcc1 = lpMPC_dv_cc + 2;
lpMPC_FLOAT lpMPC_V1[8];
lpMPC_FLOAT lpMPC_Yd1[3];
lpMPC_FLOAT lpMPC_Ld1[3];
lpMPC_FLOAT lpMPC_yy1[2];
lpMPC_FLOAT lpMPC_bmy1[2];
int lpMPC_lbIdx1[4] = {0, 1, 2, 3};
lpMPC_FLOAT* lpMPC_llb1 = lpMPC_l + 8;
lpMPC_FLOAT* lpMPC_slb1 = lpMPC_s + 8;
lpMPC_FLOAT* lpMPC_llbbyslb1 = lpMPC_lbys + 8;
lpMPC_FLOAT lpMPC_rilb1[4];
lpMPC_FLOAT* lpMPC_dllbaff1 = lpMPC_dl_aff + 8;
lpMPC_FLOAT* lpMPC_dslbaff1 = lpMPC_ds_aff + 8;
lpMPC_FLOAT* lpMPC_dllbcc1 = lpMPC_dl_cc + 8;
lpMPC_FLOAT* lpMPC_dslbcc1 = lpMPC_ds_cc + 8;
lpMPC_FLOAT* lpMPC_ccrhsl1 = lpMPC_ccrhs + 8;
int lpMPC_ubIdx1[4] = {0, 1, 2, 3};
lpMPC_FLOAT* lpMPC_lub1 = lpMPC_l + 12;
lpMPC_FLOAT* lpMPC_sub1 = lpMPC_s + 12;
lpMPC_FLOAT* lpMPC_lubbysub1 = lpMPC_lbys + 12;
lpMPC_FLOAT lpMPC_riub1[4];
lpMPC_FLOAT* lpMPC_dlubaff1 = lpMPC_dl_aff + 12;
lpMPC_FLOAT* lpMPC_dsubaff1 = lpMPC_ds_aff + 12;
lpMPC_FLOAT* lpMPC_dlubcc1 = lpMPC_dl_cc + 12;
lpMPC_FLOAT* lpMPC_dsubcc1 = lpMPC_ds_cc + 12;
lpMPC_FLOAT* lpMPC_ccrhsub1 = lpMPC_ccrhs + 12;
lpMPC_FLOAT lpMPC_Phi1[4];
lpMPC_FLOAT lpMPC_D1[4] = {-1.0000000000000000E+000, 
-1.0000000000000000E+000};
lpMPC_FLOAT lpMPC_W1[4];
lpMPC_FLOAT lpMPC_Ysd1[4];
lpMPC_FLOAT lpMPC_Lsd1[4];
lpMPC_FLOAT* lpMPC_z2 = lpMPC_z + 8;
lpMPC_FLOAT* lpMPC_dzaff2 = lpMPC_dz_aff + 8;
lpMPC_FLOAT* lpMPC_dzcc2 = lpMPC_dz_cc + 8;
lpMPC_FLOAT* lpMPC_rd2 = lpMPC_rd + 8;
lpMPC_FLOAT lpMPC_Lbyrd2[4];
lpMPC_FLOAT* lpMPC_grad_cost2 = lpMPC_grad_cost + 8;
lpMPC_FLOAT* lpMPC_grad_eq2 = lpMPC_grad_eq + 8;
lpMPC_FLOAT* lpMPC_grad_ineq2 = lpMPC_grad_ineq + 8;
lpMPC_FLOAT lpMPC_ctv2[4];
lpMPC_FLOAT* lpMPC_v2 = lpMPC_v + 4;
lpMPC_FLOAT lpMPC_re2[2];
lpMPC_FLOAT lpMPC_beta2[2];
lpMPC_FLOAT lpMPC_betacc2[2];
lpMPC_FLOAT* lpMPC_dvaff2 = lpMPC_dv_aff + 4;
lpMPC_FLOAT* lpMPC_dvcc2 = lpMPC_dv_cc + 4;
lpMPC_FLOAT lpMPC_V2[8];
lpMPC_FLOAT lpMPC_Yd2[3];
lpMPC_FLOAT lpMPC_Ld2[3];
lpMPC_FLOAT lpMPC_yy2[2];
lpMPC_FLOAT lpMPC_bmy2[2];
int lpMPC_lbIdx2[4] = {0, 1, 2, 3};
lpMPC_FLOAT* lpMPC_llb2 = lpMPC_l + 16;
lpMPC_FLOAT* lpMPC_slb2 = lpMPC_s + 16;
lpMPC_FLOAT* lpMPC_llbbyslb2 = lpMPC_lbys + 16;
lpMPC_FLOAT lpMPC_rilb2[4];
lpMPC_FLOAT* lpMPC_dllbaff2 = lpMPC_dl_aff + 16;
lpMPC_FLOAT* lpMPC_dslbaff2 = lpMPC_ds_aff + 16;
lpMPC_FLOAT* lpMPC_dllbcc2 = lpMPC_dl_cc + 16;
lpMPC_FLOAT* lpMPC_dslbcc2 = lpMPC_ds_cc + 16;
lpMPC_FLOAT* lpMPC_ccrhsl2 = lpMPC_ccrhs + 16;
int lpMPC_ubIdx2[4] = {0, 1, 2, 3};
lpMPC_FLOAT* lpMPC_lub2 = lpMPC_l + 20;
lpMPC_FLOAT* lpMPC_sub2 = lpMPC_s + 20;
lpMPC_FLOAT* lpMPC_lubbysub2 = lpMPC_lbys + 20;
lpMPC_FLOAT lpMPC_riub2[4];
lpMPC_FLOAT* lpMPC_dlubaff2 = lpMPC_dl_aff + 20;
lpMPC_FLOAT* lpMPC_dsubaff2 = lpMPC_ds_aff + 20;
lpMPC_FLOAT* lpMPC_dlubcc2 = lpMPC_dl_cc + 20;
lpMPC_FLOAT* lpMPC_dsubcc2 = lpMPC_ds_cc + 20;
lpMPC_FLOAT* lpMPC_ccrhsub2 = lpMPC_ccrhs + 20;
lpMPC_FLOAT lpMPC_Phi2[4];
lpMPC_FLOAT lpMPC_W2[4];
lpMPC_FLOAT lpMPC_Ysd2[4];
lpMPC_FLOAT lpMPC_Lsd2[4];
lpMPC_FLOAT* lpMPC_z3 = lpMPC_z + 12;
lpMPC_FLOAT* lpMPC_dzaff3 = lpMPC_dz_aff + 12;
lpMPC_FLOAT* lpMPC_dzcc3 = lpMPC_dz_cc + 12;
lpMPC_FLOAT* lpMPC_rd3 = lpMPC_rd + 12;
lpMPC_FLOAT lpMPC_Lbyrd3[4];
lpMPC_FLOAT* lpMPC_grad_cost3 = lpMPC_grad_cost + 12;
lpMPC_FLOAT* lpMPC_grad_eq3 = lpMPC_grad_eq + 12;
lpMPC_FLOAT* lpMPC_grad_ineq3 = lpMPC_grad_ineq + 12;
lpMPC_FLOAT lpMPC_ctv3[4];
lpMPC_FLOAT* lpMPC_v3 = lpMPC_v + 6;
lpMPC_FLOAT lpMPC_re3[2];
lpMPC_FLOAT lpMPC_beta3[2];
lpMPC_FLOAT lpMPC_betacc3[2];
lpMPC_FLOAT* lpMPC_dvaff3 = lpMPC_dv_aff + 6;
lpMPC_FLOAT* lpMPC_dvcc3 = lpMPC_dv_cc + 6;
lpMPC_FLOAT lpMPC_V3[8];
lpMPC_FLOAT lpMPC_Yd3[3];
lpMPC_FLOAT lpMPC_Ld3[3];
lpMPC_FLOAT lpMPC_yy3[2];
lpMPC_FLOAT lpMPC_bmy3[2];
int lpMPC_lbIdx3[4] = {0, 1, 2, 3};
lpMPC_FLOAT* lpMPC_llb3 = lpMPC_l + 24;
lpMPC_FLOAT* lpMPC_slb3 = lpMPC_s + 24;
lpMPC_FLOAT* lpMPC_llbbyslb3 = lpMPC_lbys + 24;
lpMPC_FLOAT lpMPC_rilb3[4];
lpMPC_FLOAT* lpMPC_dllbaff3 = lpMPC_dl_aff + 24;
lpMPC_FLOAT* lpMPC_dslbaff3 = lpMPC_ds_aff + 24;
lpMPC_FLOAT* lpMPC_dllbcc3 = lpMPC_dl_cc + 24;
lpMPC_FLOAT* lpMPC_dslbcc3 = lpMPC_ds_cc + 24;
lpMPC_FLOAT* lpMPC_ccrhsl3 = lpMPC_ccrhs + 24;
int lpMPC_ubIdx3[4] = {0, 1, 2, 3};
lpMPC_FLOAT* lpMPC_lub3 = lpMPC_l + 28;
lpMPC_FLOAT* lpMPC_sub3 = lpMPC_s + 28;
lpMPC_FLOAT* lpMPC_lubbysub3 = lpMPC_lbys + 28;
lpMPC_FLOAT lpMPC_riub3[4];
lpMPC_FLOAT* lpMPC_dlubaff3 = lpMPC_dl_aff + 28;
lpMPC_FLOAT* lpMPC_dsubaff3 = lpMPC_ds_aff + 28;
lpMPC_FLOAT* lpMPC_dlubcc3 = lpMPC_dl_cc + 28;
lpMPC_FLOAT* lpMPC_dsubcc3 = lpMPC_ds_cc + 28;
lpMPC_FLOAT* lpMPC_ccrhsub3 = lpMPC_ccrhs + 28;
lpMPC_FLOAT lpMPC_Phi3[4];
lpMPC_FLOAT lpMPC_W3[4];
lpMPC_FLOAT lpMPC_Ysd3[4];
lpMPC_FLOAT lpMPC_Lsd3[4];
lpMPC_FLOAT* lpMPC_z4 = lpMPC_z + 16;
lpMPC_FLOAT* lpMPC_dzaff4 = lpMPC_dz_aff + 16;
lpMPC_FLOAT* lpMPC_dzcc4 = lpMPC_dz_cc + 16;
lpMPC_FLOAT* lpMPC_rd4 = lpMPC_rd + 16;
lpMPC_FLOAT lpMPC_Lbyrd4[4];
lpMPC_FLOAT* lpMPC_grad_cost4 = lpMPC_grad_cost + 16;
lpMPC_FLOAT* lpMPC_grad_eq4 = lpMPC_grad_eq + 16;
lpMPC_FLOAT* lpMPC_grad_ineq4 = lpMPC_grad_ineq + 16;
lpMPC_FLOAT lpMPC_ctv4[4];
lpMPC_FLOAT* lpMPC_v4 = lpMPC_v + 8;
lpMPC_FLOAT lpMPC_re4[2];
lpMPC_FLOAT lpMPC_beta4[2];
lpMPC_FLOAT lpMPC_betacc4[2];
lpMPC_FLOAT* lpMPC_dvaff4 = lpMPC_dv_aff + 8;
lpMPC_FLOAT* lpMPC_dvcc4 = lpMPC_dv_cc + 8;
lpMPC_FLOAT lpMPC_V4[8];
lpMPC_FLOAT lpMPC_Yd4[3];
lpMPC_FLOAT lpMPC_Ld4[3];
lpMPC_FLOAT lpMPC_yy4[2];
lpMPC_FLOAT lpMPC_bmy4[2];
int lpMPC_lbIdx4[4] = {0, 1, 2, 3};
lpMPC_FLOAT* lpMPC_llb4 = lpMPC_l + 32;
lpMPC_FLOAT* lpMPC_slb4 = lpMPC_s + 32;
lpMPC_FLOAT* lpMPC_llbbyslb4 = lpMPC_lbys + 32;
lpMPC_FLOAT lpMPC_rilb4[4];
lpMPC_FLOAT* lpMPC_dllbaff4 = lpMPC_dl_aff + 32;
lpMPC_FLOAT* lpMPC_dslbaff4 = lpMPC_ds_aff + 32;
lpMPC_FLOAT* lpMPC_dllbcc4 = lpMPC_dl_cc + 32;
lpMPC_FLOAT* lpMPC_dslbcc4 = lpMPC_ds_cc + 32;
lpMPC_FLOAT* lpMPC_ccrhsl4 = lpMPC_ccrhs + 32;
int lpMPC_ubIdx4[4] = {0, 1, 2, 3};
lpMPC_FLOAT* lpMPC_lub4 = lpMPC_l + 36;
lpMPC_FLOAT* lpMPC_sub4 = lpMPC_s + 36;
lpMPC_FLOAT* lpMPC_lubbysub4 = lpMPC_lbys + 36;
lpMPC_FLOAT lpMPC_riub4[4];
lpMPC_FLOAT* lpMPC_dlubaff4 = lpMPC_dl_aff + 36;
lpMPC_FLOAT* lpMPC_dsubaff4 = lpMPC_ds_aff + 36;
lpMPC_FLOAT* lpMPC_dlubcc4 = lpMPC_dl_cc + 36;
lpMPC_FLOAT* lpMPC_dsubcc4 = lpMPC_ds_cc + 36;
lpMPC_FLOAT* lpMPC_ccrhsub4 = lpMPC_ccrhs + 36;
lpMPC_FLOAT lpMPC_Phi4[4];
lpMPC_FLOAT lpMPC_W4[4];
lpMPC_FLOAT lpMPC_Ysd4[4];
lpMPC_FLOAT lpMPC_Lsd4[4];
lpMPC_FLOAT* lpMPC_z5 = lpMPC_z + 20;
lpMPC_FLOAT* lpMPC_dzaff5 = lpMPC_dz_aff + 20;
lpMPC_FLOAT* lpMPC_dzcc5 = lpMPC_dz_cc + 20;
lpMPC_FLOAT* lpMPC_rd5 = lpMPC_rd + 20;
lpMPC_FLOAT lpMPC_Lbyrd5[4];
lpMPC_FLOAT* lpMPC_grad_cost5 = lpMPC_grad_cost + 20;
lpMPC_FLOAT* lpMPC_grad_eq5 = lpMPC_grad_eq + 20;
lpMPC_FLOAT* lpMPC_grad_ineq5 = lpMPC_grad_ineq + 20;
lpMPC_FLOAT lpMPC_ctv5[4];
lpMPC_FLOAT* lpMPC_v5 = lpMPC_v + 10;
lpMPC_FLOAT lpMPC_re5[2];
lpMPC_FLOAT lpMPC_beta5[2];
lpMPC_FLOAT lpMPC_betacc5[2];
lpMPC_FLOAT* lpMPC_dvaff5 = lpMPC_dv_aff + 10;
lpMPC_FLOAT* lpMPC_dvcc5 = lpMPC_dv_cc + 10;
lpMPC_FLOAT lpMPC_V5[8];
lpMPC_FLOAT lpMPC_Yd5[3];
lpMPC_FLOAT lpMPC_Ld5[3];
lpMPC_FLOAT lpMPC_yy5[2];
lpMPC_FLOAT lpMPC_bmy5[2];
int lpMPC_lbIdx5[4] = {0, 1, 2, 3};
lpMPC_FLOAT* lpMPC_llb5 = lpMPC_l + 40;
lpMPC_FLOAT* lpMPC_slb5 = lpMPC_s + 40;
lpMPC_FLOAT* lpMPC_llbbyslb5 = lpMPC_lbys + 40;
lpMPC_FLOAT lpMPC_rilb5[4];
lpMPC_FLOAT* lpMPC_dllbaff5 = lpMPC_dl_aff + 40;
lpMPC_FLOAT* lpMPC_dslbaff5 = lpMPC_ds_aff + 40;
lpMPC_FLOAT* lpMPC_dllbcc5 = lpMPC_dl_cc + 40;
lpMPC_FLOAT* lpMPC_dslbcc5 = lpMPC_ds_cc + 40;
lpMPC_FLOAT* lpMPC_ccrhsl5 = lpMPC_ccrhs + 40;
int lpMPC_ubIdx5[4] = {0, 1, 2, 3};
lpMPC_FLOAT* lpMPC_lub5 = lpMPC_l + 44;
lpMPC_FLOAT* lpMPC_sub5 = lpMPC_s + 44;
lpMPC_FLOAT* lpMPC_lubbysub5 = lpMPC_lbys + 44;
lpMPC_FLOAT lpMPC_riub5[4];
lpMPC_FLOAT* lpMPC_dlubaff5 = lpMPC_dl_aff + 44;
lpMPC_FLOAT* lpMPC_dsubaff5 = lpMPC_ds_aff + 44;
lpMPC_FLOAT* lpMPC_dlubcc5 = lpMPC_dl_cc + 44;
lpMPC_FLOAT* lpMPC_dsubcc5 = lpMPC_ds_cc + 44;
lpMPC_FLOAT* lpMPC_ccrhsub5 = lpMPC_ccrhs + 44;
lpMPC_FLOAT lpMPC_Phi5[4];
lpMPC_FLOAT lpMPC_W5[4];
lpMPC_FLOAT lpMPC_Ysd5[4];
lpMPC_FLOAT lpMPC_Lsd5[4];
lpMPC_FLOAT* lpMPC_z6 = lpMPC_z + 24;
lpMPC_FLOAT* lpMPC_dzaff6 = lpMPC_dz_aff + 24;
lpMPC_FLOAT* lpMPC_dzcc6 = lpMPC_dz_cc + 24;
lpMPC_FLOAT* lpMPC_rd6 = lpMPC_rd + 24;
lpMPC_FLOAT lpMPC_Lbyrd6[4];
lpMPC_FLOAT* lpMPC_grad_cost6 = lpMPC_grad_cost + 24;
lpMPC_FLOAT* lpMPC_grad_eq6 = lpMPC_grad_eq + 24;
lpMPC_FLOAT* lpMPC_grad_ineq6 = lpMPC_grad_ineq + 24;
lpMPC_FLOAT lpMPC_ctv6[4];
lpMPC_FLOAT* lpMPC_v6 = lpMPC_v + 12;
lpMPC_FLOAT lpMPC_re6[2];
lpMPC_FLOAT lpMPC_beta6[2];
lpMPC_FLOAT lpMPC_betacc6[2];
lpMPC_FLOAT* lpMPC_dvaff6 = lpMPC_dv_aff + 12;
lpMPC_FLOAT* lpMPC_dvcc6 = lpMPC_dv_cc + 12;
lpMPC_FLOAT lpMPC_V6[8];
lpMPC_FLOAT lpMPC_Yd6[3];
lpMPC_FLOAT lpMPC_Ld6[3];
lpMPC_FLOAT lpMPC_yy6[2];
lpMPC_FLOAT lpMPC_bmy6[2];
int lpMPC_lbIdx6[4] = {0, 1, 2, 3};
lpMPC_FLOAT* lpMPC_llb6 = lpMPC_l + 48;
lpMPC_FLOAT* lpMPC_slb6 = lpMPC_s + 48;
lpMPC_FLOAT* lpMPC_llbbyslb6 = lpMPC_lbys + 48;
lpMPC_FLOAT lpMPC_rilb6[4];
lpMPC_FLOAT* lpMPC_dllbaff6 = lpMPC_dl_aff + 48;
lpMPC_FLOAT* lpMPC_dslbaff6 = lpMPC_ds_aff + 48;
lpMPC_FLOAT* lpMPC_dllbcc6 = lpMPC_dl_cc + 48;
lpMPC_FLOAT* lpMPC_dslbcc6 = lpMPC_ds_cc + 48;
lpMPC_FLOAT* lpMPC_ccrhsl6 = lpMPC_ccrhs + 48;
int lpMPC_ubIdx6[4] = {0, 1, 2, 3};
lpMPC_FLOAT* lpMPC_lub6 = lpMPC_l + 52;
lpMPC_FLOAT* lpMPC_sub6 = lpMPC_s + 52;
lpMPC_FLOAT* lpMPC_lubbysub6 = lpMPC_lbys + 52;
lpMPC_FLOAT lpMPC_riub6[4];
lpMPC_FLOAT* lpMPC_dlubaff6 = lpMPC_dl_aff + 52;
lpMPC_FLOAT* lpMPC_dsubaff6 = lpMPC_ds_aff + 52;
lpMPC_FLOAT* lpMPC_dlubcc6 = lpMPC_dl_cc + 52;
lpMPC_FLOAT* lpMPC_dsubcc6 = lpMPC_ds_cc + 52;
lpMPC_FLOAT* lpMPC_ccrhsub6 = lpMPC_ccrhs + 52;
lpMPC_FLOAT lpMPC_Phi6[4];
lpMPC_FLOAT lpMPC_W6[4];
lpMPC_FLOAT lpMPC_Ysd6[4];
lpMPC_FLOAT lpMPC_Lsd6[4];
lpMPC_FLOAT* lpMPC_z7 = lpMPC_z + 28;
lpMPC_FLOAT* lpMPC_dzaff7 = lpMPC_dz_aff + 28;
lpMPC_FLOAT* lpMPC_dzcc7 = lpMPC_dz_cc + 28;
lpMPC_FLOAT* lpMPC_rd7 = lpMPC_rd + 28;
lpMPC_FLOAT lpMPC_Lbyrd7[4];
lpMPC_FLOAT* lpMPC_grad_cost7 = lpMPC_grad_cost + 28;
lpMPC_FLOAT* lpMPC_grad_eq7 = lpMPC_grad_eq + 28;
lpMPC_FLOAT* lpMPC_grad_ineq7 = lpMPC_grad_ineq + 28;
lpMPC_FLOAT lpMPC_ctv7[4];
lpMPC_FLOAT* lpMPC_v7 = lpMPC_v + 14;
lpMPC_FLOAT lpMPC_re7[2];
lpMPC_FLOAT lpMPC_beta7[2];
lpMPC_FLOAT lpMPC_betacc7[2];
lpMPC_FLOAT* lpMPC_dvaff7 = lpMPC_dv_aff + 14;
lpMPC_FLOAT* lpMPC_dvcc7 = lpMPC_dv_cc + 14;
lpMPC_FLOAT lpMPC_V7[8];
lpMPC_FLOAT lpMPC_Yd7[3];
lpMPC_FLOAT lpMPC_Ld7[3];
lpMPC_FLOAT lpMPC_yy7[2];
lpMPC_FLOAT lpMPC_bmy7[2];
int lpMPC_lbIdx7[4] = {0, 1, 2, 3};
lpMPC_FLOAT* lpMPC_llb7 = lpMPC_l + 56;
lpMPC_FLOAT* lpMPC_slb7 = lpMPC_s + 56;
lpMPC_FLOAT* lpMPC_llbbyslb7 = lpMPC_lbys + 56;
lpMPC_FLOAT lpMPC_rilb7[4];
lpMPC_FLOAT* lpMPC_dllbaff7 = lpMPC_dl_aff + 56;
lpMPC_FLOAT* lpMPC_dslbaff7 = lpMPC_ds_aff + 56;
lpMPC_FLOAT* lpMPC_dllbcc7 = lpMPC_dl_cc + 56;
lpMPC_FLOAT* lpMPC_dslbcc7 = lpMPC_ds_cc + 56;
lpMPC_FLOAT* lpMPC_ccrhsl7 = lpMPC_ccrhs + 56;
int lpMPC_ubIdx7[4] = {0, 1, 2, 3};
lpMPC_FLOAT* lpMPC_lub7 = lpMPC_l + 60;
lpMPC_FLOAT* lpMPC_sub7 = lpMPC_s + 60;
lpMPC_FLOAT* lpMPC_lubbysub7 = lpMPC_lbys + 60;
lpMPC_FLOAT lpMPC_riub7[4];
lpMPC_FLOAT* lpMPC_dlubaff7 = lpMPC_dl_aff + 60;
lpMPC_FLOAT* lpMPC_dsubaff7 = lpMPC_ds_aff + 60;
lpMPC_FLOAT* lpMPC_dlubcc7 = lpMPC_dl_cc + 60;
lpMPC_FLOAT* lpMPC_dsubcc7 = lpMPC_ds_cc + 60;
lpMPC_FLOAT* lpMPC_ccrhsub7 = lpMPC_ccrhs + 60;
lpMPC_FLOAT lpMPC_Phi7[4];
lpMPC_FLOAT lpMPC_W7[4];
lpMPC_FLOAT lpMPC_Ysd7[4];
lpMPC_FLOAT lpMPC_Lsd7[4];
lpMPC_FLOAT* lpMPC_z8 = lpMPC_z + 32;
lpMPC_FLOAT* lpMPC_dzaff8 = lpMPC_dz_aff + 32;
lpMPC_FLOAT* lpMPC_dzcc8 = lpMPC_dz_cc + 32;
lpMPC_FLOAT* lpMPC_rd8 = lpMPC_rd + 32;
lpMPC_FLOAT lpMPC_Lbyrd8[4];
lpMPC_FLOAT* lpMPC_grad_cost8 = lpMPC_grad_cost + 32;
lpMPC_FLOAT* lpMPC_grad_eq8 = lpMPC_grad_eq + 32;
lpMPC_FLOAT* lpMPC_grad_ineq8 = lpMPC_grad_ineq + 32;
lpMPC_FLOAT lpMPC_ctv8[4];
lpMPC_FLOAT* lpMPC_v8 = lpMPC_v + 16;
lpMPC_FLOAT lpMPC_re8[2];
lpMPC_FLOAT lpMPC_beta8[2];
lpMPC_FLOAT lpMPC_betacc8[2];
lpMPC_FLOAT* lpMPC_dvaff8 = lpMPC_dv_aff + 16;
lpMPC_FLOAT* lpMPC_dvcc8 = lpMPC_dv_cc + 16;
lpMPC_FLOAT lpMPC_V8[8];
lpMPC_FLOAT lpMPC_Yd8[3];
lpMPC_FLOAT lpMPC_Ld8[3];
lpMPC_FLOAT lpMPC_yy8[2];
lpMPC_FLOAT lpMPC_bmy8[2];
int lpMPC_lbIdx8[4] = {0, 1, 2, 3};
lpMPC_FLOAT* lpMPC_llb8 = lpMPC_l + 64;
lpMPC_FLOAT* lpMPC_slb8 = lpMPC_s + 64;
lpMPC_FLOAT* lpMPC_llbbyslb8 = lpMPC_lbys + 64;
lpMPC_FLOAT lpMPC_rilb8[4];
lpMPC_FLOAT* lpMPC_dllbaff8 = lpMPC_dl_aff + 64;
lpMPC_FLOAT* lpMPC_dslbaff8 = lpMPC_ds_aff + 64;
lpMPC_FLOAT* lpMPC_dllbcc8 = lpMPC_dl_cc + 64;
lpMPC_FLOAT* lpMPC_dslbcc8 = lpMPC_ds_cc + 64;
lpMPC_FLOAT* lpMPC_ccrhsl8 = lpMPC_ccrhs + 64;
int lpMPC_ubIdx8[4] = {0, 1, 2, 3};
lpMPC_FLOAT* lpMPC_lub8 = lpMPC_l + 68;
lpMPC_FLOAT* lpMPC_sub8 = lpMPC_s + 68;
lpMPC_FLOAT* lpMPC_lubbysub8 = lpMPC_lbys + 68;
lpMPC_FLOAT lpMPC_riub8[4];
lpMPC_FLOAT* lpMPC_dlubaff8 = lpMPC_dl_aff + 68;
lpMPC_FLOAT* lpMPC_dsubaff8 = lpMPC_ds_aff + 68;
lpMPC_FLOAT* lpMPC_dlubcc8 = lpMPC_dl_cc + 68;
lpMPC_FLOAT* lpMPC_dsubcc8 = lpMPC_ds_cc + 68;
lpMPC_FLOAT* lpMPC_ccrhsub8 = lpMPC_ccrhs + 68;
lpMPC_FLOAT lpMPC_Phi8[4];
lpMPC_FLOAT lpMPC_W8[4];
lpMPC_FLOAT lpMPC_Ysd8[4];
lpMPC_FLOAT lpMPC_Lsd8[4];
lpMPC_FLOAT* lpMPC_z9 = lpMPC_z + 36;
lpMPC_FLOAT* lpMPC_dzaff9 = lpMPC_dz_aff + 36;
lpMPC_FLOAT* lpMPC_dzcc9 = lpMPC_dz_cc + 36;
lpMPC_FLOAT* lpMPC_rd9 = lpMPC_rd + 36;
lpMPC_FLOAT lpMPC_Lbyrd9[2];
lpMPC_FLOAT* lpMPC_grad_cost9 = lpMPC_grad_cost + 36;
lpMPC_FLOAT* lpMPC_grad_eq9 = lpMPC_grad_eq + 36;
lpMPC_FLOAT* lpMPC_grad_ineq9 = lpMPC_grad_ineq + 36;
lpMPC_FLOAT lpMPC_ctv9[2];
lpMPC_FLOAT* lpMPC_v9 = lpMPC_v + 18;
lpMPC_FLOAT lpMPC_re9[2];
lpMPC_FLOAT lpMPC_beta9[2];
lpMPC_FLOAT lpMPC_betacc9[2];
lpMPC_FLOAT* lpMPC_dvaff9 = lpMPC_dv_aff + 18;
lpMPC_FLOAT* lpMPC_dvcc9 = lpMPC_dv_cc + 18;
lpMPC_FLOAT lpMPC_V9[4];
lpMPC_FLOAT lpMPC_Yd9[3];
lpMPC_FLOAT lpMPC_Ld9[3];
lpMPC_FLOAT lpMPC_yy9[2];
lpMPC_FLOAT lpMPC_bmy9[2];
int lpMPC_lbIdx9[2] = {0, 1};
lpMPC_FLOAT* lpMPC_llb9 = lpMPC_l + 72;
lpMPC_FLOAT* lpMPC_slb9 = lpMPC_s + 72;
lpMPC_FLOAT* lpMPC_llbbyslb9 = lpMPC_lbys + 72;
lpMPC_FLOAT lpMPC_rilb9[2];
lpMPC_FLOAT* lpMPC_dllbaff9 = lpMPC_dl_aff + 72;
lpMPC_FLOAT* lpMPC_dslbaff9 = lpMPC_ds_aff + 72;
lpMPC_FLOAT* lpMPC_dllbcc9 = lpMPC_dl_cc + 72;
lpMPC_FLOAT* lpMPC_dslbcc9 = lpMPC_ds_cc + 72;
lpMPC_FLOAT* lpMPC_ccrhsl9 = lpMPC_ccrhs + 72;
int lpMPC_ubIdx9[2] = {0, 1};
lpMPC_FLOAT* lpMPC_lub9 = lpMPC_l + 74;
lpMPC_FLOAT* lpMPC_sub9 = lpMPC_s + 74;
lpMPC_FLOAT* lpMPC_lubbysub9 = lpMPC_lbys + 74;
lpMPC_FLOAT lpMPC_riub9[2];
lpMPC_FLOAT* lpMPC_dlubaff9 = lpMPC_dl_aff + 74;
lpMPC_FLOAT* lpMPC_dsubaff9 = lpMPC_ds_aff + 74;
lpMPC_FLOAT* lpMPC_dlubcc9 = lpMPC_dl_cc + 74;
lpMPC_FLOAT* lpMPC_dsubcc9 = lpMPC_ds_cc + 74;
lpMPC_FLOAT* lpMPC_ccrhsub9 = lpMPC_ccrhs + 74;
lpMPC_FLOAT lpMPC_Phi9[2];
lpMPC_FLOAT lpMPC_D9[2] = {-1.0000000000000000E+000, 
-1.0000000000000000E+000};
lpMPC_FLOAT lpMPC_W9[2];
lpMPC_FLOAT lpMPC_Ysd9[4];
lpMPC_FLOAT lpMPC_Lsd9[4];
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
lpMPC_LA_INITIALIZEVECTOR_38(lpMPC_z, 0);
lpMPC_LA_INITIALIZEVECTOR_20(lpMPC_v, 1);
lpMPC_LA_INITIALIZEVECTOR_76(lpMPC_l, 1);
lpMPC_LA_INITIALIZEVECTOR_76(lpMPC_s, 1);
info->mu = 0;
lpMPC_LA_DOTACC_76(lpMPC_l, lpMPC_s, &info->mu);
info->mu /= 76;
while( 1 ){
info->pobj = 0;
lpMPC_LA_DIAG_QUADFCN_4(lpMPC_H0, params->f1, lpMPC_z0, lpMPC_grad_cost0, &info->pobj);
lpMPC_LA_DIAG_QUADFCN_4(lpMPC_H0, params->f2, lpMPC_z1, lpMPC_grad_cost1, &info->pobj);
lpMPC_LA_DIAG_QUADFCN_4(lpMPC_H0, params->f3, lpMPC_z2, lpMPC_grad_cost2, &info->pobj);
lpMPC_LA_DIAG_QUADFCN_4(lpMPC_H0, params->f4, lpMPC_z3, lpMPC_grad_cost3, &info->pobj);
lpMPC_LA_DIAG_QUADFCN_4(lpMPC_H0, params->f5, lpMPC_z4, lpMPC_grad_cost4, &info->pobj);
lpMPC_LA_DIAG_QUADFCN_4(lpMPC_H0, params->f6, lpMPC_z5, lpMPC_grad_cost5, &info->pobj);
lpMPC_LA_DIAG_QUADFCN_4(lpMPC_H0, params->f7, lpMPC_z6, lpMPC_grad_cost6, &info->pobj);
lpMPC_LA_DIAG_QUADFCN_4(lpMPC_H0, params->f8, lpMPC_z7, lpMPC_grad_cost7, &info->pobj);
lpMPC_LA_DIAG_QUADFCN_4(lpMPC_H0, params->f9, lpMPC_z8, lpMPC_grad_cost8, &info->pobj);
lpMPC_LA_DIAG_QUADFCN_2(lpMPC_H0, params->f10, lpMPC_z9, lpMPC_grad_cost9, &info->pobj);
info->res_eq = 0;
info->dgap = 0;
lpMPC_LA_DIAGZERO_MVMSUB6_2(lpMPC_D0, lpMPC_z0, params->e1, lpMPC_v0, lpMPC_re0, &info->dgap, &info->res_eq);
lpMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_4_4(params->C1, lpMPC_z0, lpMPC_D1, lpMPC_z1, params->e2, lpMPC_v1, lpMPC_re1, &info->dgap, &info->res_eq);
lpMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_4_4(params->C2, lpMPC_z1, lpMPC_D1, lpMPC_z2, params->e3, lpMPC_v2, lpMPC_re2, &info->dgap, &info->res_eq);
lpMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_4_4(params->C3, lpMPC_z2, lpMPC_D1, lpMPC_z3, params->e4, lpMPC_v3, lpMPC_re3, &info->dgap, &info->res_eq);
lpMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_4_4(params->C4, lpMPC_z3, lpMPC_D1, lpMPC_z4, params->e5, lpMPC_v4, lpMPC_re4, &info->dgap, &info->res_eq);
lpMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_4_4(params->C5, lpMPC_z4, lpMPC_D1, lpMPC_z5, params->e6, lpMPC_v5, lpMPC_re5, &info->dgap, &info->res_eq);
lpMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_4_4(params->C6, lpMPC_z5, lpMPC_D1, lpMPC_z6, params->e7, lpMPC_v6, lpMPC_re6, &info->dgap, &info->res_eq);
lpMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_4_4(params->C7, lpMPC_z6, lpMPC_D1, lpMPC_z7, params->e8, lpMPC_v7, lpMPC_re7, &info->dgap, &info->res_eq);
lpMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_4_4(params->C8, lpMPC_z7, lpMPC_D1, lpMPC_z8, params->e9, lpMPC_v8, lpMPC_re8, &info->dgap, &info->res_eq);
lpMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_4_2(params->C9, lpMPC_z8, lpMPC_D9, lpMPC_z9, params->e10, lpMPC_v9, lpMPC_re9, &info->dgap, &info->res_eq);
lpMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C1, lpMPC_v1, lpMPC_D0, lpMPC_v0, lpMPC_grad_eq0);
lpMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C2, lpMPC_v2, lpMPC_D1, lpMPC_v1, lpMPC_grad_eq1);
lpMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C3, lpMPC_v3, lpMPC_D1, lpMPC_v2, lpMPC_grad_eq2);
lpMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C4, lpMPC_v4, lpMPC_D1, lpMPC_v3, lpMPC_grad_eq3);
lpMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C5, lpMPC_v5, lpMPC_D1, lpMPC_v4, lpMPC_grad_eq4);
lpMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C6, lpMPC_v6, lpMPC_D1, lpMPC_v5, lpMPC_grad_eq5);
lpMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C7, lpMPC_v7, lpMPC_D1, lpMPC_v6, lpMPC_grad_eq6);
lpMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C8, lpMPC_v8, lpMPC_D1, lpMPC_v7, lpMPC_grad_eq7);
lpMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C9, lpMPC_v9, lpMPC_D1, lpMPC_v8, lpMPC_grad_eq8);
lpMPC_LA_DIAGZERO_MTVM_2_2(lpMPC_D9, lpMPC_v9, lpMPC_grad_eq9);
info->res_ineq = 0;
lpMPC_LA_VSUBADD3_4(params->lb1, lpMPC_z0, lpMPC_lbIdx0, lpMPC_llb0, lpMPC_slb0, lpMPC_rilb0, &info->dgap, &info->res_ineq);
lpMPC_LA_VSUBADD2_4(lpMPC_z0, lpMPC_ubIdx0, params->ub1, lpMPC_lub0, lpMPC_sub0, lpMPC_riub0, &info->dgap, &info->res_ineq);
lpMPC_LA_VSUBADD3_4(params->lb2, lpMPC_z1, lpMPC_lbIdx1, lpMPC_llb1, lpMPC_slb1, lpMPC_rilb1, &info->dgap, &info->res_ineq);
lpMPC_LA_VSUBADD2_4(lpMPC_z1, lpMPC_ubIdx1, params->ub2, lpMPC_lub1, lpMPC_sub1, lpMPC_riub1, &info->dgap, &info->res_ineq);
lpMPC_LA_VSUBADD3_4(params->lb3, lpMPC_z2, lpMPC_lbIdx2, lpMPC_llb2, lpMPC_slb2, lpMPC_rilb2, &info->dgap, &info->res_ineq);
lpMPC_LA_VSUBADD2_4(lpMPC_z2, lpMPC_ubIdx2, params->ub3, lpMPC_lub2, lpMPC_sub2, lpMPC_riub2, &info->dgap, &info->res_ineq);
lpMPC_LA_VSUBADD3_4(params->lb4, lpMPC_z3, lpMPC_lbIdx3, lpMPC_llb3, lpMPC_slb3, lpMPC_rilb3, &info->dgap, &info->res_ineq);
lpMPC_LA_VSUBADD2_4(lpMPC_z3, lpMPC_ubIdx3, params->ub4, lpMPC_lub3, lpMPC_sub3, lpMPC_riub3, &info->dgap, &info->res_ineq);
lpMPC_LA_VSUBADD3_4(params->lb5, lpMPC_z4, lpMPC_lbIdx4, lpMPC_llb4, lpMPC_slb4, lpMPC_rilb4, &info->dgap, &info->res_ineq);
lpMPC_LA_VSUBADD2_4(lpMPC_z4, lpMPC_ubIdx4, params->ub5, lpMPC_lub4, lpMPC_sub4, lpMPC_riub4, &info->dgap, &info->res_ineq);
lpMPC_LA_VSUBADD3_4(params->lb6, lpMPC_z5, lpMPC_lbIdx5, lpMPC_llb5, lpMPC_slb5, lpMPC_rilb5, &info->dgap, &info->res_ineq);
lpMPC_LA_VSUBADD2_4(lpMPC_z5, lpMPC_ubIdx5, params->ub6, lpMPC_lub5, lpMPC_sub5, lpMPC_riub5, &info->dgap, &info->res_ineq);
lpMPC_LA_VSUBADD3_4(params->lb7, lpMPC_z6, lpMPC_lbIdx6, lpMPC_llb6, lpMPC_slb6, lpMPC_rilb6, &info->dgap, &info->res_ineq);
lpMPC_LA_VSUBADD2_4(lpMPC_z6, lpMPC_ubIdx6, params->ub7, lpMPC_lub6, lpMPC_sub6, lpMPC_riub6, &info->dgap, &info->res_ineq);
lpMPC_LA_VSUBADD3_4(params->lb8, lpMPC_z7, lpMPC_lbIdx7, lpMPC_llb7, lpMPC_slb7, lpMPC_rilb7, &info->dgap, &info->res_ineq);
lpMPC_LA_VSUBADD2_4(lpMPC_z7, lpMPC_ubIdx7, params->ub8, lpMPC_lub7, lpMPC_sub7, lpMPC_riub7, &info->dgap, &info->res_ineq);
lpMPC_LA_VSUBADD3_4(params->lb9, lpMPC_z8, lpMPC_lbIdx8, lpMPC_llb8, lpMPC_slb8, lpMPC_rilb8, &info->dgap, &info->res_ineq);
lpMPC_LA_VSUBADD2_4(lpMPC_z8, lpMPC_ubIdx8, params->ub9, lpMPC_lub8, lpMPC_sub8, lpMPC_riub8, &info->dgap, &info->res_ineq);
lpMPC_LA_VSUBADD3_2(params->lb10, lpMPC_z9, lpMPC_lbIdx9, lpMPC_llb9, lpMPC_slb9, lpMPC_rilb9, &info->dgap, &info->res_ineq);
lpMPC_LA_VSUBADD2_2(lpMPC_z9, lpMPC_ubIdx9, params->ub10, lpMPC_lub9, lpMPC_sub9, lpMPC_riub9, &info->dgap, &info->res_ineq);
lpMPC_LA_INEQ_B_GRAD_4_4_4(lpMPC_lub0, lpMPC_sub0, lpMPC_riub0, lpMPC_llb0, lpMPC_slb0, lpMPC_rilb0, lpMPC_lbIdx0, lpMPC_ubIdx0, lpMPC_grad_ineq0, lpMPC_lubbysub0, lpMPC_llbbyslb0);
lpMPC_LA_INEQ_B_GRAD_4_4_4(lpMPC_lub1, lpMPC_sub1, lpMPC_riub1, lpMPC_llb1, lpMPC_slb1, lpMPC_rilb1, lpMPC_lbIdx1, lpMPC_ubIdx1, lpMPC_grad_ineq1, lpMPC_lubbysub1, lpMPC_llbbyslb1);
lpMPC_LA_INEQ_B_GRAD_4_4_4(lpMPC_lub2, lpMPC_sub2, lpMPC_riub2, lpMPC_llb2, lpMPC_slb2, lpMPC_rilb2, lpMPC_lbIdx2, lpMPC_ubIdx2, lpMPC_grad_ineq2, lpMPC_lubbysub2, lpMPC_llbbyslb2);
lpMPC_LA_INEQ_B_GRAD_4_4_4(lpMPC_lub3, lpMPC_sub3, lpMPC_riub3, lpMPC_llb3, lpMPC_slb3, lpMPC_rilb3, lpMPC_lbIdx3, lpMPC_ubIdx3, lpMPC_grad_ineq3, lpMPC_lubbysub3, lpMPC_llbbyslb3);
lpMPC_LA_INEQ_B_GRAD_4_4_4(lpMPC_lub4, lpMPC_sub4, lpMPC_riub4, lpMPC_llb4, lpMPC_slb4, lpMPC_rilb4, lpMPC_lbIdx4, lpMPC_ubIdx4, lpMPC_grad_ineq4, lpMPC_lubbysub4, lpMPC_llbbyslb4);
lpMPC_LA_INEQ_B_GRAD_4_4_4(lpMPC_lub5, lpMPC_sub5, lpMPC_riub5, lpMPC_llb5, lpMPC_slb5, lpMPC_rilb5, lpMPC_lbIdx5, lpMPC_ubIdx5, lpMPC_grad_ineq5, lpMPC_lubbysub5, lpMPC_llbbyslb5);
lpMPC_LA_INEQ_B_GRAD_4_4_4(lpMPC_lub6, lpMPC_sub6, lpMPC_riub6, lpMPC_llb6, lpMPC_slb6, lpMPC_rilb6, lpMPC_lbIdx6, lpMPC_ubIdx6, lpMPC_grad_ineq6, lpMPC_lubbysub6, lpMPC_llbbyslb6);
lpMPC_LA_INEQ_B_GRAD_4_4_4(lpMPC_lub7, lpMPC_sub7, lpMPC_riub7, lpMPC_llb7, lpMPC_slb7, lpMPC_rilb7, lpMPC_lbIdx7, lpMPC_ubIdx7, lpMPC_grad_ineq7, lpMPC_lubbysub7, lpMPC_llbbyslb7);
lpMPC_LA_INEQ_B_GRAD_4_4_4(lpMPC_lub8, lpMPC_sub8, lpMPC_riub8, lpMPC_llb8, lpMPC_slb8, lpMPC_rilb8, lpMPC_lbIdx8, lpMPC_ubIdx8, lpMPC_grad_ineq8, lpMPC_lubbysub8, lpMPC_llbbyslb8);
lpMPC_LA_INEQ_B_GRAD_2_2_2(lpMPC_lub9, lpMPC_sub9, lpMPC_riub9, lpMPC_llb9, lpMPC_slb9, lpMPC_rilb9, lpMPC_lbIdx9, lpMPC_ubIdx9, lpMPC_grad_ineq9, lpMPC_lubbysub9, lpMPC_llbbyslb9);
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
lpMPC_LA_VVADD3_38(lpMPC_grad_cost, lpMPC_grad_eq, lpMPC_grad_ineq, lpMPC_rd);
lpMPC_LA_DIAG_CHOL_ONELOOP_LBUB_4_4_4(lpMPC_H0, lpMPC_llbbyslb0, lpMPC_lbIdx0, lpMPC_lubbysub0, lpMPC_ubIdx0, lpMPC_Phi0);
lpMPC_LA_DIAG_MATRIXFORWARDSUB_2_4(lpMPC_Phi0, params->C1, lpMPC_V0);
lpMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_4(lpMPC_Phi0, lpMPC_D0, lpMPC_W0);
lpMPC_LA_DENSE_DIAGZERO_MMTM_2_4_2(lpMPC_W0, lpMPC_V0, lpMPC_Ysd1);
lpMPC_LA_DIAG_FORWARDSUB_4(lpMPC_Phi0, lpMPC_rd0, lpMPC_Lbyrd0);
lpMPC_LA_DIAG_CHOL_ONELOOP_LBUB_4_4_4(lpMPC_H0, lpMPC_llbbyslb1, lpMPC_lbIdx1, lpMPC_lubbysub1, lpMPC_ubIdx1, lpMPC_Phi1);
lpMPC_LA_DIAG_MATRIXFORWARDSUB_2_4(lpMPC_Phi1, params->C2, lpMPC_V1);
lpMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_4(lpMPC_Phi1, lpMPC_D1, lpMPC_W1);
lpMPC_LA_DENSE_DIAGZERO_MMTM_2_4_2(lpMPC_W1, lpMPC_V1, lpMPC_Ysd2);
lpMPC_LA_DIAG_FORWARDSUB_4(lpMPC_Phi1, lpMPC_rd1, lpMPC_Lbyrd1);
lpMPC_LA_DIAG_CHOL_ONELOOP_LBUB_4_4_4(lpMPC_H0, lpMPC_llbbyslb2, lpMPC_lbIdx2, lpMPC_lubbysub2, lpMPC_ubIdx2, lpMPC_Phi2);
lpMPC_LA_DIAG_MATRIXFORWARDSUB_2_4(lpMPC_Phi2, params->C3, lpMPC_V2);
lpMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_4(lpMPC_Phi2, lpMPC_D1, lpMPC_W2);
lpMPC_LA_DENSE_DIAGZERO_MMTM_2_4_2(lpMPC_W2, lpMPC_V2, lpMPC_Ysd3);
lpMPC_LA_DIAG_FORWARDSUB_4(lpMPC_Phi2, lpMPC_rd2, lpMPC_Lbyrd2);
lpMPC_LA_DIAG_CHOL_ONELOOP_LBUB_4_4_4(lpMPC_H0, lpMPC_llbbyslb3, lpMPC_lbIdx3, lpMPC_lubbysub3, lpMPC_ubIdx3, lpMPC_Phi3);
lpMPC_LA_DIAG_MATRIXFORWARDSUB_2_4(lpMPC_Phi3, params->C4, lpMPC_V3);
lpMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_4(lpMPC_Phi3, lpMPC_D1, lpMPC_W3);
lpMPC_LA_DENSE_DIAGZERO_MMTM_2_4_2(lpMPC_W3, lpMPC_V3, lpMPC_Ysd4);
lpMPC_LA_DIAG_FORWARDSUB_4(lpMPC_Phi3, lpMPC_rd3, lpMPC_Lbyrd3);
lpMPC_LA_DIAG_CHOL_ONELOOP_LBUB_4_4_4(lpMPC_H0, lpMPC_llbbyslb4, lpMPC_lbIdx4, lpMPC_lubbysub4, lpMPC_ubIdx4, lpMPC_Phi4);
lpMPC_LA_DIAG_MATRIXFORWARDSUB_2_4(lpMPC_Phi4, params->C5, lpMPC_V4);
lpMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_4(lpMPC_Phi4, lpMPC_D1, lpMPC_W4);
lpMPC_LA_DENSE_DIAGZERO_MMTM_2_4_2(lpMPC_W4, lpMPC_V4, lpMPC_Ysd5);
lpMPC_LA_DIAG_FORWARDSUB_4(lpMPC_Phi4, lpMPC_rd4, lpMPC_Lbyrd4);
lpMPC_LA_DIAG_CHOL_ONELOOP_LBUB_4_4_4(lpMPC_H0, lpMPC_llbbyslb5, lpMPC_lbIdx5, lpMPC_lubbysub5, lpMPC_ubIdx5, lpMPC_Phi5);
lpMPC_LA_DIAG_MATRIXFORWARDSUB_2_4(lpMPC_Phi5, params->C6, lpMPC_V5);
lpMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_4(lpMPC_Phi5, lpMPC_D1, lpMPC_W5);
lpMPC_LA_DENSE_DIAGZERO_MMTM_2_4_2(lpMPC_W5, lpMPC_V5, lpMPC_Ysd6);
lpMPC_LA_DIAG_FORWARDSUB_4(lpMPC_Phi5, lpMPC_rd5, lpMPC_Lbyrd5);
lpMPC_LA_DIAG_CHOL_ONELOOP_LBUB_4_4_4(lpMPC_H0, lpMPC_llbbyslb6, lpMPC_lbIdx6, lpMPC_lubbysub6, lpMPC_ubIdx6, lpMPC_Phi6);
lpMPC_LA_DIAG_MATRIXFORWARDSUB_2_4(lpMPC_Phi6, params->C7, lpMPC_V6);
lpMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_4(lpMPC_Phi6, lpMPC_D1, lpMPC_W6);
lpMPC_LA_DENSE_DIAGZERO_MMTM_2_4_2(lpMPC_W6, lpMPC_V6, lpMPC_Ysd7);
lpMPC_LA_DIAG_FORWARDSUB_4(lpMPC_Phi6, lpMPC_rd6, lpMPC_Lbyrd6);
lpMPC_LA_DIAG_CHOL_ONELOOP_LBUB_4_4_4(lpMPC_H0, lpMPC_llbbyslb7, lpMPC_lbIdx7, lpMPC_lubbysub7, lpMPC_ubIdx7, lpMPC_Phi7);
lpMPC_LA_DIAG_MATRIXFORWARDSUB_2_4(lpMPC_Phi7, params->C8, lpMPC_V7);
lpMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_4(lpMPC_Phi7, lpMPC_D1, lpMPC_W7);
lpMPC_LA_DENSE_DIAGZERO_MMTM_2_4_2(lpMPC_W7, lpMPC_V7, lpMPC_Ysd8);
lpMPC_LA_DIAG_FORWARDSUB_4(lpMPC_Phi7, lpMPC_rd7, lpMPC_Lbyrd7);
lpMPC_LA_DIAG_CHOL_ONELOOP_LBUB_4_4_4(lpMPC_H0, lpMPC_llbbyslb8, lpMPC_lbIdx8, lpMPC_lubbysub8, lpMPC_ubIdx8, lpMPC_Phi8);
lpMPC_LA_DIAG_MATRIXFORWARDSUB_2_4(lpMPC_Phi8, params->C9, lpMPC_V8);
lpMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_4(lpMPC_Phi8, lpMPC_D1, lpMPC_W8);
lpMPC_LA_DENSE_DIAGZERO_MMTM_2_4_2(lpMPC_W8, lpMPC_V8, lpMPC_Ysd9);
lpMPC_LA_DIAG_FORWARDSUB_4(lpMPC_Phi8, lpMPC_rd8, lpMPC_Lbyrd8);
lpMPC_LA_DIAG_CHOL_ONELOOP_LBUB_2_2_2(lpMPC_H0, lpMPC_llbbyslb9, lpMPC_lbIdx9, lpMPC_lubbysub9, lpMPC_ubIdx9, lpMPC_Phi9);
lpMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_2(lpMPC_Phi9, lpMPC_D9, lpMPC_W9);
lpMPC_LA_DIAG_FORWARDSUB_2(lpMPC_Phi9, lpMPC_rd9, lpMPC_Lbyrd9);
lpMPC_LA_DIAGZERO_MMT_2(lpMPC_W0, lpMPC_Yd0);
lpMPC_LA_DIAGZERO_MVMSUB7_2(lpMPC_W0, lpMPC_Lbyrd0, lpMPC_re0, lpMPC_beta0);
lpMPC_LA_DENSE_DIAGZERO_MMT2_2_4_4(lpMPC_V0, lpMPC_W1, lpMPC_Yd1);
lpMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_4_4(lpMPC_V0, lpMPC_Lbyrd0, lpMPC_W1, lpMPC_Lbyrd1, lpMPC_re1, lpMPC_beta1);
lpMPC_LA_DENSE_DIAGZERO_MMT2_2_4_4(lpMPC_V1, lpMPC_W2, lpMPC_Yd2);
lpMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_4_4(lpMPC_V1, lpMPC_Lbyrd1, lpMPC_W2, lpMPC_Lbyrd2, lpMPC_re2, lpMPC_beta2);
lpMPC_LA_DENSE_DIAGZERO_MMT2_2_4_4(lpMPC_V2, lpMPC_W3, lpMPC_Yd3);
lpMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_4_4(lpMPC_V2, lpMPC_Lbyrd2, lpMPC_W3, lpMPC_Lbyrd3, lpMPC_re3, lpMPC_beta3);
lpMPC_LA_DENSE_DIAGZERO_MMT2_2_4_4(lpMPC_V3, lpMPC_W4, lpMPC_Yd4);
lpMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_4_4(lpMPC_V3, lpMPC_Lbyrd3, lpMPC_W4, lpMPC_Lbyrd4, lpMPC_re4, lpMPC_beta4);
lpMPC_LA_DENSE_DIAGZERO_MMT2_2_4_4(lpMPC_V4, lpMPC_W5, lpMPC_Yd5);
lpMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_4_4(lpMPC_V4, lpMPC_Lbyrd4, lpMPC_W5, lpMPC_Lbyrd5, lpMPC_re5, lpMPC_beta5);
lpMPC_LA_DENSE_DIAGZERO_MMT2_2_4_4(lpMPC_V5, lpMPC_W6, lpMPC_Yd6);
lpMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_4_4(lpMPC_V5, lpMPC_Lbyrd5, lpMPC_W6, lpMPC_Lbyrd6, lpMPC_re6, lpMPC_beta6);
lpMPC_LA_DENSE_DIAGZERO_MMT2_2_4_4(lpMPC_V6, lpMPC_W7, lpMPC_Yd7);
lpMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_4_4(lpMPC_V6, lpMPC_Lbyrd6, lpMPC_W7, lpMPC_Lbyrd7, lpMPC_re7, lpMPC_beta7);
lpMPC_LA_DENSE_DIAGZERO_MMT2_2_4_4(lpMPC_V7, lpMPC_W8, lpMPC_Yd8);
lpMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_4_4(lpMPC_V7, lpMPC_Lbyrd7, lpMPC_W8, lpMPC_Lbyrd8, lpMPC_re8, lpMPC_beta8);
lpMPC_LA_DENSE_DIAGZERO_MMT2_2_4_2(lpMPC_V8, lpMPC_W9, lpMPC_Yd9);
lpMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_4_2(lpMPC_V8, lpMPC_Lbyrd8, lpMPC_W9, lpMPC_Lbyrd9, lpMPC_re9, lpMPC_beta9);
lpMPC_LA_DENSE_CHOL_2(lpMPC_Yd0, lpMPC_Ld0);
lpMPC_LA_DENSE_FORWARDSUB_2(lpMPC_Ld0, lpMPC_beta0, lpMPC_yy0);
lpMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(lpMPC_Ld0, lpMPC_Ysd1, lpMPC_Lsd1);
lpMPC_LA_DENSE_MMTSUB_2_2(lpMPC_Lsd1, lpMPC_Yd1);
lpMPC_LA_DENSE_CHOL_2(lpMPC_Yd1, lpMPC_Ld1);
lpMPC_LA_DENSE_MVMSUB1_2_2(lpMPC_Lsd1, lpMPC_yy0, lpMPC_beta1, lpMPC_bmy1);
lpMPC_LA_DENSE_FORWARDSUB_2(lpMPC_Ld1, lpMPC_bmy1, lpMPC_yy1);
lpMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(lpMPC_Ld1, lpMPC_Ysd2, lpMPC_Lsd2);
lpMPC_LA_DENSE_MMTSUB_2_2(lpMPC_Lsd2, lpMPC_Yd2);
lpMPC_LA_DENSE_CHOL_2(lpMPC_Yd2, lpMPC_Ld2);
lpMPC_LA_DENSE_MVMSUB1_2_2(lpMPC_Lsd2, lpMPC_yy1, lpMPC_beta2, lpMPC_bmy2);
lpMPC_LA_DENSE_FORWARDSUB_2(lpMPC_Ld2, lpMPC_bmy2, lpMPC_yy2);
lpMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(lpMPC_Ld2, lpMPC_Ysd3, lpMPC_Lsd3);
lpMPC_LA_DENSE_MMTSUB_2_2(lpMPC_Lsd3, lpMPC_Yd3);
lpMPC_LA_DENSE_CHOL_2(lpMPC_Yd3, lpMPC_Ld3);
lpMPC_LA_DENSE_MVMSUB1_2_2(lpMPC_Lsd3, lpMPC_yy2, lpMPC_beta3, lpMPC_bmy3);
lpMPC_LA_DENSE_FORWARDSUB_2(lpMPC_Ld3, lpMPC_bmy3, lpMPC_yy3);
lpMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(lpMPC_Ld3, lpMPC_Ysd4, lpMPC_Lsd4);
lpMPC_LA_DENSE_MMTSUB_2_2(lpMPC_Lsd4, lpMPC_Yd4);
lpMPC_LA_DENSE_CHOL_2(lpMPC_Yd4, lpMPC_Ld4);
lpMPC_LA_DENSE_MVMSUB1_2_2(lpMPC_Lsd4, lpMPC_yy3, lpMPC_beta4, lpMPC_bmy4);
lpMPC_LA_DENSE_FORWARDSUB_2(lpMPC_Ld4, lpMPC_bmy4, lpMPC_yy4);
lpMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(lpMPC_Ld4, lpMPC_Ysd5, lpMPC_Lsd5);
lpMPC_LA_DENSE_MMTSUB_2_2(lpMPC_Lsd5, lpMPC_Yd5);
lpMPC_LA_DENSE_CHOL_2(lpMPC_Yd5, lpMPC_Ld5);
lpMPC_LA_DENSE_MVMSUB1_2_2(lpMPC_Lsd5, lpMPC_yy4, lpMPC_beta5, lpMPC_bmy5);
lpMPC_LA_DENSE_FORWARDSUB_2(lpMPC_Ld5, lpMPC_bmy5, lpMPC_yy5);
lpMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(lpMPC_Ld5, lpMPC_Ysd6, lpMPC_Lsd6);
lpMPC_LA_DENSE_MMTSUB_2_2(lpMPC_Lsd6, lpMPC_Yd6);
lpMPC_LA_DENSE_CHOL_2(lpMPC_Yd6, lpMPC_Ld6);
lpMPC_LA_DENSE_MVMSUB1_2_2(lpMPC_Lsd6, lpMPC_yy5, lpMPC_beta6, lpMPC_bmy6);
lpMPC_LA_DENSE_FORWARDSUB_2(lpMPC_Ld6, lpMPC_bmy6, lpMPC_yy6);
lpMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(lpMPC_Ld6, lpMPC_Ysd7, lpMPC_Lsd7);
lpMPC_LA_DENSE_MMTSUB_2_2(lpMPC_Lsd7, lpMPC_Yd7);
lpMPC_LA_DENSE_CHOL_2(lpMPC_Yd7, lpMPC_Ld7);
lpMPC_LA_DENSE_MVMSUB1_2_2(lpMPC_Lsd7, lpMPC_yy6, lpMPC_beta7, lpMPC_bmy7);
lpMPC_LA_DENSE_FORWARDSUB_2(lpMPC_Ld7, lpMPC_bmy7, lpMPC_yy7);
lpMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(lpMPC_Ld7, lpMPC_Ysd8, lpMPC_Lsd8);
lpMPC_LA_DENSE_MMTSUB_2_2(lpMPC_Lsd8, lpMPC_Yd8);
lpMPC_LA_DENSE_CHOL_2(lpMPC_Yd8, lpMPC_Ld8);
lpMPC_LA_DENSE_MVMSUB1_2_2(lpMPC_Lsd8, lpMPC_yy7, lpMPC_beta8, lpMPC_bmy8);
lpMPC_LA_DENSE_FORWARDSUB_2(lpMPC_Ld8, lpMPC_bmy8, lpMPC_yy8);
lpMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(lpMPC_Ld8, lpMPC_Ysd9, lpMPC_Lsd9);
lpMPC_LA_DENSE_MMTSUB_2_2(lpMPC_Lsd9, lpMPC_Yd9);
lpMPC_LA_DENSE_CHOL_2(lpMPC_Yd9, lpMPC_Ld9);
lpMPC_LA_DENSE_MVMSUB1_2_2(lpMPC_Lsd9, lpMPC_yy8, lpMPC_beta9, lpMPC_bmy9);
lpMPC_LA_DENSE_FORWARDSUB_2(lpMPC_Ld9, lpMPC_bmy9, lpMPC_yy9);
lpMPC_LA_DENSE_BACKWARDSUB_2(lpMPC_Ld9, lpMPC_yy9, lpMPC_dvaff9);
lpMPC_LA_DENSE_MTVMSUB_2_2(lpMPC_Lsd9, lpMPC_dvaff9, lpMPC_yy8, lpMPC_bmy8);
lpMPC_LA_DENSE_BACKWARDSUB_2(lpMPC_Ld8, lpMPC_bmy8, lpMPC_dvaff8);
lpMPC_LA_DENSE_MTVMSUB_2_2(lpMPC_Lsd8, lpMPC_dvaff8, lpMPC_yy7, lpMPC_bmy7);
lpMPC_LA_DENSE_BACKWARDSUB_2(lpMPC_Ld7, lpMPC_bmy7, lpMPC_dvaff7);
lpMPC_LA_DENSE_MTVMSUB_2_2(lpMPC_Lsd7, lpMPC_dvaff7, lpMPC_yy6, lpMPC_bmy6);
lpMPC_LA_DENSE_BACKWARDSUB_2(lpMPC_Ld6, lpMPC_bmy6, lpMPC_dvaff6);
lpMPC_LA_DENSE_MTVMSUB_2_2(lpMPC_Lsd6, lpMPC_dvaff6, lpMPC_yy5, lpMPC_bmy5);
lpMPC_LA_DENSE_BACKWARDSUB_2(lpMPC_Ld5, lpMPC_bmy5, lpMPC_dvaff5);
lpMPC_LA_DENSE_MTVMSUB_2_2(lpMPC_Lsd5, lpMPC_dvaff5, lpMPC_yy4, lpMPC_bmy4);
lpMPC_LA_DENSE_BACKWARDSUB_2(lpMPC_Ld4, lpMPC_bmy4, lpMPC_dvaff4);
lpMPC_LA_DENSE_MTVMSUB_2_2(lpMPC_Lsd4, lpMPC_dvaff4, lpMPC_yy3, lpMPC_bmy3);
lpMPC_LA_DENSE_BACKWARDSUB_2(lpMPC_Ld3, lpMPC_bmy3, lpMPC_dvaff3);
lpMPC_LA_DENSE_MTVMSUB_2_2(lpMPC_Lsd3, lpMPC_dvaff3, lpMPC_yy2, lpMPC_bmy2);
lpMPC_LA_DENSE_BACKWARDSUB_2(lpMPC_Ld2, lpMPC_bmy2, lpMPC_dvaff2);
lpMPC_LA_DENSE_MTVMSUB_2_2(lpMPC_Lsd2, lpMPC_dvaff2, lpMPC_yy1, lpMPC_bmy1);
lpMPC_LA_DENSE_BACKWARDSUB_2(lpMPC_Ld1, lpMPC_bmy1, lpMPC_dvaff1);
lpMPC_LA_DENSE_MTVMSUB_2_2(lpMPC_Lsd1, lpMPC_dvaff1, lpMPC_yy0, lpMPC_bmy0);
lpMPC_LA_DENSE_BACKWARDSUB_2(lpMPC_Ld0, lpMPC_bmy0, lpMPC_dvaff0);
lpMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C1, lpMPC_dvaff1, lpMPC_D0, lpMPC_dvaff0, lpMPC_grad_eq0);
lpMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C2, lpMPC_dvaff2, lpMPC_D1, lpMPC_dvaff1, lpMPC_grad_eq1);
lpMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C3, lpMPC_dvaff3, lpMPC_D1, lpMPC_dvaff2, lpMPC_grad_eq2);
lpMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C4, lpMPC_dvaff4, lpMPC_D1, lpMPC_dvaff3, lpMPC_grad_eq3);
lpMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C5, lpMPC_dvaff5, lpMPC_D1, lpMPC_dvaff4, lpMPC_grad_eq4);
lpMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C6, lpMPC_dvaff6, lpMPC_D1, lpMPC_dvaff5, lpMPC_grad_eq5);
lpMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C7, lpMPC_dvaff7, lpMPC_D1, lpMPC_dvaff6, lpMPC_grad_eq6);
lpMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C8, lpMPC_dvaff8, lpMPC_D1, lpMPC_dvaff7, lpMPC_grad_eq7);
lpMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C9, lpMPC_dvaff9, lpMPC_D1, lpMPC_dvaff8, lpMPC_grad_eq8);
lpMPC_LA_DIAGZERO_MTVM_2_2(lpMPC_D9, lpMPC_dvaff9, lpMPC_grad_eq9);
lpMPC_LA_VSUB2_38(lpMPC_rd, lpMPC_grad_eq, lpMPC_rd);
lpMPC_LA_DIAG_FORWARDBACKWARDSUB_4(lpMPC_Phi0, lpMPC_rd0, lpMPC_dzaff0);
lpMPC_LA_DIAG_FORWARDBACKWARDSUB_4(lpMPC_Phi1, lpMPC_rd1, lpMPC_dzaff1);
lpMPC_LA_DIAG_FORWARDBACKWARDSUB_4(lpMPC_Phi2, lpMPC_rd2, lpMPC_dzaff2);
lpMPC_LA_DIAG_FORWARDBACKWARDSUB_4(lpMPC_Phi3, lpMPC_rd3, lpMPC_dzaff3);
lpMPC_LA_DIAG_FORWARDBACKWARDSUB_4(lpMPC_Phi4, lpMPC_rd4, lpMPC_dzaff4);
lpMPC_LA_DIAG_FORWARDBACKWARDSUB_4(lpMPC_Phi5, lpMPC_rd5, lpMPC_dzaff5);
lpMPC_LA_DIAG_FORWARDBACKWARDSUB_4(lpMPC_Phi6, lpMPC_rd6, lpMPC_dzaff6);
lpMPC_LA_DIAG_FORWARDBACKWARDSUB_4(lpMPC_Phi7, lpMPC_rd7, lpMPC_dzaff7);
lpMPC_LA_DIAG_FORWARDBACKWARDSUB_4(lpMPC_Phi8, lpMPC_rd8, lpMPC_dzaff8);
lpMPC_LA_DIAG_FORWARDBACKWARDSUB_2(lpMPC_Phi9, lpMPC_rd9, lpMPC_dzaff9);
lpMPC_LA_VSUB_INDEXED_4(lpMPC_dzaff0, lpMPC_lbIdx0, lpMPC_rilb0, lpMPC_dslbaff0);
lpMPC_LA_VSUB3_4(lpMPC_llbbyslb0, lpMPC_dslbaff0, lpMPC_llb0, lpMPC_dllbaff0);
lpMPC_LA_VSUB2_INDEXED_4(lpMPC_riub0, lpMPC_dzaff0, lpMPC_ubIdx0, lpMPC_dsubaff0);
lpMPC_LA_VSUB3_4(lpMPC_lubbysub0, lpMPC_dsubaff0, lpMPC_lub0, lpMPC_dlubaff0);
lpMPC_LA_VSUB_INDEXED_4(lpMPC_dzaff1, lpMPC_lbIdx1, lpMPC_rilb1, lpMPC_dslbaff1);
lpMPC_LA_VSUB3_4(lpMPC_llbbyslb1, lpMPC_dslbaff1, lpMPC_llb1, lpMPC_dllbaff1);
lpMPC_LA_VSUB2_INDEXED_4(lpMPC_riub1, lpMPC_dzaff1, lpMPC_ubIdx1, lpMPC_dsubaff1);
lpMPC_LA_VSUB3_4(lpMPC_lubbysub1, lpMPC_dsubaff1, lpMPC_lub1, lpMPC_dlubaff1);
lpMPC_LA_VSUB_INDEXED_4(lpMPC_dzaff2, lpMPC_lbIdx2, lpMPC_rilb2, lpMPC_dslbaff2);
lpMPC_LA_VSUB3_4(lpMPC_llbbyslb2, lpMPC_dslbaff2, lpMPC_llb2, lpMPC_dllbaff2);
lpMPC_LA_VSUB2_INDEXED_4(lpMPC_riub2, lpMPC_dzaff2, lpMPC_ubIdx2, lpMPC_dsubaff2);
lpMPC_LA_VSUB3_4(lpMPC_lubbysub2, lpMPC_dsubaff2, lpMPC_lub2, lpMPC_dlubaff2);
lpMPC_LA_VSUB_INDEXED_4(lpMPC_dzaff3, lpMPC_lbIdx3, lpMPC_rilb3, lpMPC_dslbaff3);
lpMPC_LA_VSUB3_4(lpMPC_llbbyslb3, lpMPC_dslbaff3, lpMPC_llb3, lpMPC_dllbaff3);
lpMPC_LA_VSUB2_INDEXED_4(lpMPC_riub3, lpMPC_dzaff3, lpMPC_ubIdx3, lpMPC_dsubaff3);
lpMPC_LA_VSUB3_4(lpMPC_lubbysub3, lpMPC_dsubaff3, lpMPC_lub3, lpMPC_dlubaff3);
lpMPC_LA_VSUB_INDEXED_4(lpMPC_dzaff4, lpMPC_lbIdx4, lpMPC_rilb4, lpMPC_dslbaff4);
lpMPC_LA_VSUB3_4(lpMPC_llbbyslb4, lpMPC_dslbaff4, lpMPC_llb4, lpMPC_dllbaff4);
lpMPC_LA_VSUB2_INDEXED_4(lpMPC_riub4, lpMPC_dzaff4, lpMPC_ubIdx4, lpMPC_dsubaff4);
lpMPC_LA_VSUB3_4(lpMPC_lubbysub4, lpMPC_dsubaff4, lpMPC_lub4, lpMPC_dlubaff4);
lpMPC_LA_VSUB_INDEXED_4(lpMPC_dzaff5, lpMPC_lbIdx5, lpMPC_rilb5, lpMPC_dslbaff5);
lpMPC_LA_VSUB3_4(lpMPC_llbbyslb5, lpMPC_dslbaff5, lpMPC_llb5, lpMPC_dllbaff5);
lpMPC_LA_VSUB2_INDEXED_4(lpMPC_riub5, lpMPC_dzaff5, lpMPC_ubIdx5, lpMPC_dsubaff5);
lpMPC_LA_VSUB3_4(lpMPC_lubbysub5, lpMPC_dsubaff5, lpMPC_lub5, lpMPC_dlubaff5);
lpMPC_LA_VSUB_INDEXED_4(lpMPC_dzaff6, lpMPC_lbIdx6, lpMPC_rilb6, lpMPC_dslbaff6);
lpMPC_LA_VSUB3_4(lpMPC_llbbyslb6, lpMPC_dslbaff6, lpMPC_llb6, lpMPC_dllbaff6);
lpMPC_LA_VSUB2_INDEXED_4(lpMPC_riub6, lpMPC_dzaff6, lpMPC_ubIdx6, lpMPC_dsubaff6);
lpMPC_LA_VSUB3_4(lpMPC_lubbysub6, lpMPC_dsubaff6, lpMPC_lub6, lpMPC_dlubaff6);
lpMPC_LA_VSUB_INDEXED_4(lpMPC_dzaff7, lpMPC_lbIdx7, lpMPC_rilb7, lpMPC_dslbaff7);
lpMPC_LA_VSUB3_4(lpMPC_llbbyslb7, lpMPC_dslbaff7, lpMPC_llb7, lpMPC_dllbaff7);
lpMPC_LA_VSUB2_INDEXED_4(lpMPC_riub7, lpMPC_dzaff7, lpMPC_ubIdx7, lpMPC_dsubaff7);
lpMPC_LA_VSUB3_4(lpMPC_lubbysub7, lpMPC_dsubaff7, lpMPC_lub7, lpMPC_dlubaff7);
lpMPC_LA_VSUB_INDEXED_4(lpMPC_dzaff8, lpMPC_lbIdx8, lpMPC_rilb8, lpMPC_dslbaff8);
lpMPC_LA_VSUB3_4(lpMPC_llbbyslb8, lpMPC_dslbaff8, lpMPC_llb8, lpMPC_dllbaff8);
lpMPC_LA_VSUB2_INDEXED_4(lpMPC_riub8, lpMPC_dzaff8, lpMPC_ubIdx8, lpMPC_dsubaff8);
lpMPC_LA_VSUB3_4(lpMPC_lubbysub8, lpMPC_dsubaff8, lpMPC_lub8, lpMPC_dlubaff8);
lpMPC_LA_VSUB_INDEXED_2(lpMPC_dzaff9, lpMPC_lbIdx9, lpMPC_rilb9, lpMPC_dslbaff9);
lpMPC_LA_VSUB3_2(lpMPC_llbbyslb9, lpMPC_dslbaff9, lpMPC_llb9, lpMPC_dllbaff9);
lpMPC_LA_VSUB2_INDEXED_2(lpMPC_riub9, lpMPC_dzaff9, lpMPC_ubIdx9, lpMPC_dsubaff9);
lpMPC_LA_VSUB3_2(lpMPC_lubbysub9, lpMPC_dsubaff9, lpMPC_lub9, lpMPC_dlubaff9);
info->lsit_aff = lpMPC_LINESEARCH_BACKTRACKING_AFFINE(lpMPC_l, lpMPC_s, lpMPC_dl_aff, lpMPC_ds_aff, &info->step_aff, &info->mu_aff);
if( info->lsit_aff == lpMPC_NOPROGRESS ){
exitcode = lpMPC_NOPROGRESS; break;
}
sigma_3rdroot = info->mu_aff / info->mu;
info->sigma = sigma_3rdroot*sigma_3rdroot*sigma_3rdroot;
musigma = info->mu * info->sigma;
lpMPC_LA_VSUB5_76(lpMPC_ds_aff, lpMPC_dl_aff, musigma, lpMPC_ccrhs);
lpMPC_LA_VSUB6_INDEXED_4_4_4(lpMPC_ccrhsub0, lpMPC_sub0, lpMPC_ubIdx0, lpMPC_ccrhsl0, lpMPC_slb0, lpMPC_lbIdx0, lpMPC_rd0);
lpMPC_LA_VSUB6_INDEXED_4_4_4(lpMPC_ccrhsub1, lpMPC_sub1, lpMPC_ubIdx1, lpMPC_ccrhsl1, lpMPC_slb1, lpMPC_lbIdx1, lpMPC_rd1);
lpMPC_LA_DIAG_FORWARDSUB_4(lpMPC_Phi0, lpMPC_rd0, lpMPC_Lbyrd0);
lpMPC_LA_DIAG_FORWARDSUB_4(lpMPC_Phi1, lpMPC_rd1, lpMPC_Lbyrd1);
lpMPC_LA_DIAGZERO_MVM_2(lpMPC_W0, lpMPC_Lbyrd0, lpMPC_beta0);
lpMPC_LA_DENSE_FORWARDSUB_2(lpMPC_Ld0, lpMPC_beta0, lpMPC_yy0);
lpMPC_LA_DENSE_DIAGZERO_2MVMADD_2_4_4(lpMPC_V0, lpMPC_Lbyrd0, lpMPC_W1, lpMPC_Lbyrd1, lpMPC_beta1);
lpMPC_LA_DENSE_MVMSUB1_2_2(lpMPC_Lsd1, lpMPC_yy0, lpMPC_beta1, lpMPC_bmy1);
lpMPC_LA_DENSE_FORWARDSUB_2(lpMPC_Ld1, lpMPC_bmy1, lpMPC_yy1);
lpMPC_LA_VSUB6_INDEXED_4_4_4(lpMPC_ccrhsub2, lpMPC_sub2, lpMPC_ubIdx2, lpMPC_ccrhsl2, lpMPC_slb2, lpMPC_lbIdx2, lpMPC_rd2);
lpMPC_LA_DIAG_FORWARDSUB_4(lpMPC_Phi2, lpMPC_rd2, lpMPC_Lbyrd2);
lpMPC_LA_DENSE_DIAGZERO_2MVMADD_2_4_4(lpMPC_V1, lpMPC_Lbyrd1, lpMPC_W2, lpMPC_Lbyrd2, lpMPC_beta2);
lpMPC_LA_DENSE_MVMSUB1_2_2(lpMPC_Lsd2, lpMPC_yy1, lpMPC_beta2, lpMPC_bmy2);
lpMPC_LA_DENSE_FORWARDSUB_2(lpMPC_Ld2, lpMPC_bmy2, lpMPC_yy2);
lpMPC_LA_VSUB6_INDEXED_4_4_4(lpMPC_ccrhsub3, lpMPC_sub3, lpMPC_ubIdx3, lpMPC_ccrhsl3, lpMPC_slb3, lpMPC_lbIdx3, lpMPC_rd3);
lpMPC_LA_DIAG_FORWARDSUB_4(lpMPC_Phi3, lpMPC_rd3, lpMPC_Lbyrd3);
lpMPC_LA_DENSE_DIAGZERO_2MVMADD_2_4_4(lpMPC_V2, lpMPC_Lbyrd2, lpMPC_W3, lpMPC_Lbyrd3, lpMPC_beta3);
lpMPC_LA_DENSE_MVMSUB1_2_2(lpMPC_Lsd3, lpMPC_yy2, lpMPC_beta3, lpMPC_bmy3);
lpMPC_LA_DENSE_FORWARDSUB_2(lpMPC_Ld3, lpMPC_bmy3, lpMPC_yy3);
lpMPC_LA_VSUB6_INDEXED_4_4_4(lpMPC_ccrhsub4, lpMPC_sub4, lpMPC_ubIdx4, lpMPC_ccrhsl4, lpMPC_slb4, lpMPC_lbIdx4, lpMPC_rd4);
lpMPC_LA_DIAG_FORWARDSUB_4(lpMPC_Phi4, lpMPC_rd4, lpMPC_Lbyrd4);
lpMPC_LA_DENSE_DIAGZERO_2MVMADD_2_4_4(lpMPC_V3, lpMPC_Lbyrd3, lpMPC_W4, lpMPC_Lbyrd4, lpMPC_beta4);
lpMPC_LA_DENSE_MVMSUB1_2_2(lpMPC_Lsd4, lpMPC_yy3, lpMPC_beta4, lpMPC_bmy4);
lpMPC_LA_DENSE_FORWARDSUB_2(lpMPC_Ld4, lpMPC_bmy4, lpMPC_yy4);
lpMPC_LA_VSUB6_INDEXED_4_4_4(lpMPC_ccrhsub5, lpMPC_sub5, lpMPC_ubIdx5, lpMPC_ccrhsl5, lpMPC_slb5, lpMPC_lbIdx5, lpMPC_rd5);
lpMPC_LA_DIAG_FORWARDSUB_4(lpMPC_Phi5, lpMPC_rd5, lpMPC_Lbyrd5);
lpMPC_LA_DENSE_DIAGZERO_2MVMADD_2_4_4(lpMPC_V4, lpMPC_Lbyrd4, lpMPC_W5, lpMPC_Lbyrd5, lpMPC_beta5);
lpMPC_LA_DENSE_MVMSUB1_2_2(lpMPC_Lsd5, lpMPC_yy4, lpMPC_beta5, lpMPC_bmy5);
lpMPC_LA_DENSE_FORWARDSUB_2(lpMPC_Ld5, lpMPC_bmy5, lpMPC_yy5);
lpMPC_LA_VSUB6_INDEXED_4_4_4(lpMPC_ccrhsub6, lpMPC_sub6, lpMPC_ubIdx6, lpMPC_ccrhsl6, lpMPC_slb6, lpMPC_lbIdx6, lpMPC_rd6);
lpMPC_LA_DIAG_FORWARDSUB_4(lpMPC_Phi6, lpMPC_rd6, lpMPC_Lbyrd6);
lpMPC_LA_DENSE_DIAGZERO_2MVMADD_2_4_4(lpMPC_V5, lpMPC_Lbyrd5, lpMPC_W6, lpMPC_Lbyrd6, lpMPC_beta6);
lpMPC_LA_DENSE_MVMSUB1_2_2(lpMPC_Lsd6, lpMPC_yy5, lpMPC_beta6, lpMPC_bmy6);
lpMPC_LA_DENSE_FORWARDSUB_2(lpMPC_Ld6, lpMPC_bmy6, lpMPC_yy6);
lpMPC_LA_VSUB6_INDEXED_4_4_4(lpMPC_ccrhsub7, lpMPC_sub7, lpMPC_ubIdx7, lpMPC_ccrhsl7, lpMPC_slb7, lpMPC_lbIdx7, lpMPC_rd7);
lpMPC_LA_DIAG_FORWARDSUB_4(lpMPC_Phi7, lpMPC_rd7, lpMPC_Lbyrd7);
lpMPC_LA_DENSE_DIAGZERO_2MVMADD_2_4_4(lpMPC_V6, lpMPC_Lbyrd6, lpMPC_W7, lpMPC_Lbyrd7, lpMPC_beta7);
lpMPC_LA_DENSE_MVMSUB1_2_2(lpMPC_Lsd7, lpMPC_yy6, lpMPC_beta7, lpMPC_bmy7);
lpMPC_LA_DENSE_FORWARDSUB_2(lpMPC_Ld7, lpMPC_bmy7, lpMPC_yy7);
lpMPC_LA_VSUB6_INDEXED_4_4_4(lpMPC_ccrhsub8, lpMPC_sub8, lpMPC_ubIdx8, lpMPC_ccrhsl8, lpMPC_slb8, lpMPC_lbIdx8, lpMPC_rd8);
lpMPC_LA_DIAG_FORWARDSUB_4(lpMPC_Phi8, lpMPC_rd8, lpMPC_Lbyrd8);
lpMPC_LA_DENSE_DIAGZERO_2MVMADD_2_4_4(lpMPC_V7, lpMPC_Lbyrd7, lpMPC_W8, lpMPC_Lbyrd8, lpMPC_beta8);
lpMPC_LA_DENSE_MVMSUB1_2_2(lpMPC_Lsd8, lpMPC_yy7, lpMPC_beta8, lpMPC_bmy8);
lpMPC_LA_DENSE_FORWARDSUB_2(lpMPC_Ld8, lpMPC_bmy8, lpMPC_yy8);
lpMPC_LA_VSUB6_INDEXED_2_2_2(lpMPC_ccrhsub9, lpMPC_sub9, lpMPC_ubIdx9, lpMPC_ccrhsl9, lpMPC_slb9, lpMPC_lbIdx9, lpMPC_rd9);
lpMPC_LA_DIAG_FORWARDSUB_2(lpMPC_Phi9, lpMPC_rd9, lpMPC_Lbyrd9);
lpMPC_LA_DENSE_DIAGZERO_2MVMADD_2_4_2(lpMPC_V8, lpMPC_Lbyrd8, lpMPC_W9, lpMPC_Lbyrd9, lpMPC_beta9);
lpMPC_LA_DENSE_MVMSUB1_2_2(lpMPC_Lsd9, lpMPC_yy8, lpMPC_beta9, lpMPC_bmy9);
lpMPC_LA_DENSE_FORWARDSUB_2(lpMPC_Ld9, lpMPC_bmy9, lpMPC_yy9);
lpMPC_LA_DENSE_BACKWARDSUB_2(lpMPC_Ld9, lpMPC_yy9, lpMPC_dvcc9);
lpMPC_LA_DENSE_MTVMSUB_2_2(lpMPC_Lsd9, lpMPC_dvcc9, lpMPC_yy8, lpMPC_bmy8);
lpMPC_LA_DENSE_BACKWARDSUB_2(lpMPC_Ld8, lpMPC_bmy8, lpMPC_dvcc8);
lpMPC_LA_DENSE_MTVMSUB_2_2(lpMPC_Lsd8, lpMPC_dvcc8, lpMPC_yy7, lpMPC_bmy7);
lpMPC_LA_DENSE_BACKWARDSUB_2(lpMPC_Ld7, lpMPC_bmy7, lpMPC_dvcc7);
lpMPC_LA_DENSE_MTVMSUB_2_2(lpMPC_Lsd7, lpMPC_dvcc7, lpMPC_yy6, lpMPC_bmy6);
lpMPC_LA_DENSE_BACKWARDSUB_2(lpMPC_Ld6, lpMPC_bmy6, lpMPC_dvcc6);
lpMPC_LA_DENSE_MTVMSUB_2_2(lpMPC_Lsd6, lpMPC_dvcc6, lpMPC_yy5, lpMPC_bmy5);
lpMPC_LA_DENSE_BACKWARDSUB_2(lpMPC_Ld5, lpMPC_bmy5, lpMPC_dvcc5);
lpMPC_LA_DENSE_MTVMSUB_2_2(lpMPC_Lsd5, lpMPC_dvcc5, lpMPC_yy4, lpMPC_bmy4);
lpMPC_LA_DENSE_BACKWARDSUB_2(lpMPC_Ld4, lpMPC_bmy4, lpMPC_dvcc4);
lpMPC_LA_DENSE_MTVMSUB_2_2(lpMPC_Lsd4, lpMPC_dvcc4, lpMPC_yy3, lpMPC_bmy3);
lpMPC_LA_DENSE_BACKWARDSUB_2(lpMPC_Ld3, lpMPC_bmy3, lpMPC_dvcc3);
lpMPC_LA_DENSE_MTVMSUB_2_2(lpMPC_Lsd3, lpMPC_dvcc3, lpMPC_yy2, lpMPC_bmy2);
lpMPC_LA_DENSE_BACKWARDSUB_2(lpMPC_Ld2, lpMPC_bmy2, lpMPC_dvcc2);
lpMPC_LA_DENSE_MTVMSUB_2_2(lpMPC_Lsd2, lpMPC_dvcc2, lpMPC_yy1, lpMPC_bmy1);
lpMPC_LA_DENSE_BACKWARDSUB_2(lpMPC_Ld1, lpMPC_bmy1, lpMPC_dvcc1);
lpMPC_LA_DENSE_MTVMSUB_2_2(lpMPC_Lsd1, lpMPC_dvcc1, lpMPC_yy0, lpMPC_bmy0);
lpMPC_LA_DENSE_BACKWARDSUB_2(lpMPC_Ld0, lpMPC_bmy0, lpMPC_dvcc0);
lpMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C1, lpMPC_dvcc1, lpMPC_D0, lpMPC_dvcc0, lpMPC_grad_eq0);
lpMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C2, lpMPC_dvcc2, lpMPC_D1, lpMPC_dvcc1, lpMPC_grad_eq1);
lpMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C3, lpMPC_dvcc3, lpMPC_D1, lpMPC_dvcc2, lpMPC_grad_eq2);
lpMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C4, lpMPC_dvcc4, lpMPC_D1, lpMPC_dvcc3, lpMPC_grad_eq3);
lpMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C5, lpMPC_dvcc5, lpMPC_D1, lpMPC_dvcc4, lpMPC_grad_eq4);
lpMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C6, lpMPC_dvcc6, lpMPC_D1, lpMPC_dvcc5, lpMPC_grad_eq5);
lpMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C7, lpMPC_dvcc7, lpMPC_D1, lpMPC_dvcc6, lpMPC_grad_eq6);
lpMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C8, lpMPC_dvcc8, lpMPC_D1, lpMPC_dvcc7, lpMPC_grad_eq7);
lpMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C9, lpMPC_dvcc9, lpMPC_D1, lpMPC_dvcc8, lpMPC_grad_eq8);
lpMPC_LA_DIAGZERO_MTVM_2_2(lpMPC_D9, lpMPC_dvcc9, lpMPC_grad_eq9);
lpMPC_LA_VSUB_38(lpMPC_rd, lpMPC_grad_eq, lpMPC_rd);
lpMPC_LA_DIAG_FORWARDBACKWARDSUB_4(lpMPC_Phi0, lpMPC_rd0, lpMPC_dzcc0);
lpMPC_LA_DIAG_FORWARDBACKWARDSUB_4(lpMPC_Phi1, lpMPC_rd1, lpMPC_dzcc1);
lpMPC_LA_DIAG_FORWARDBACKWARDSUB_4(lpMPC_Phi2, lpMPC_rd2, lpMPC_dzcc2);
lpMPC_LA_DIAG_FORWARDBACKWARDSUB_4(lpMPC_Phi3, lpMPC_rd3, lpMPC_dzcc3);
lpMPC_LA_DIAG_FORWARDBACKWARDSUB_4(lpMPC_Phi4, lpMPC_rd4, lpMPC_dzcc4);
lpMPC_LA_DIAG_FORWARDBACKWARDSUB_4(lpMPC_Phi5, lpMPC_rd5, lpMPC_dzcc5);
lpMPC_LA_DIAG_FORWARDBACKWARDSUB_4(lpMPC_Phi6, lpMPC_rd6, lpMPC_dzcc6);
lpMPC_LA_DIAG_FORWARDBACKWARDSUB_4(lpMPC_Phi7, lpMPC_rd7, lpMPC_dzcc7);
lpMPC_LA_DIAG_FORWARDBACKWARDSUB_4(lpMPC_Phi8, lpMPC_rd8, lpMPC_dzcc8);
lpMPC_LA_DIAG_FORWARDBACKWARDSUB_2(lpMPC_Phi9, lpMPC_rd9, lpMPC_dzcc9);
lpMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_4(lpMPC_ccrhsl0, lpMPC_slb0, lpMPC_llbbyslb0, lpMPC_dzcc0, lpMPC_lbIdx0, lpMPC_dllbcc0);
lpMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_4(lpMPC_ccrhsub0, lpMPC_sub0, lpMPC_lubbysub0, lpMPC_dzcc0, lpMPC_ubIdx0, lpMPC_dlubcc0);
lpMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_4(lpMPC_ccrhsl1, lpMPC_slb1, lpMPC_llbbyslb1, lpMPC_dzcc1, lpMPC_lbIdx1, lpMPC_dllbcc1);
lpMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_4(lpMPC_ccrhsub1, lpMPC_sub1, lpMPC_lubbysub1, lpMPC_dzcc1, lpMPC_ubIdx1, lpMPC_dlubcc1);
lpMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_4(lpMPC_ccrhsl2, lpMPC_slb2, lpMPC_llbbyslb2, lpMPC_dzcc2, lpMPC_lbIdx2, lpMPC_dllbcc2);
lpMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_4(lpMPC_ccrhsub2, lpMPC_sub2, lpMPC_lubbysub2, lpMPC_dzcc2, lpMPC_ubIdx2, lpMPC_dlubcc2);
lpMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_4(lpMPC_ccrhsl3, lpMPC_slb3, lpMPC_llbbyslb3, lpMPC_dzcc3, lpMPC_lbIdx3, lpMPC_dllbcc3);
lpMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_4(lpMPC_ccrhsub3, lpMPC_sub3, lpMPC_lubbysub3, lpMPC_dzcc3, lpMPC_ubIdx3, lpMPC_dlubcc3);
lpMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_4(lpMPC_ccrhsl4, lpMPC_slb4, lpMPC_llbbyslb4, lpMPC_dzcc4, lpMPC_lbIdx4, lpMPC_dllbcc4);
lpMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_4(lpMPC_ccrhsub4, lpMPC_sub4, lpMPC_lubbysub4, lpMPC_dzcc4, lpMPC_ubIdx4, lpMPC_dlubcc4);
lpMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_4(lpMPC_ccrhsl5, lpMPC_slb5, lpMPC_llbbyslb5, lpMPC_dzcc5, lpMPC_lbIdx5, lpMPC_dllbcc5);
lpMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_4(lpMPC_ccrhsub5, lpMPC_sub5, lpMPC_lubbysub5, lpMPC_dzcc5, lpMPC_ubIdx5, lpMPC_dlubcc5);
lpMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_4(lpMPC_ccrhsl6, lpMPC_slb6, lpMPC_llbbyslb6, lpMPC_dzcc6, lpMPC_lbIdx6, lpMPC_dllbcc6);
lpMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_4(lpMPC_ccrhsub6, lpMPC_sub6, lpMPC_lubbysub6, lpMPC_dzcc6, lpMPC_ubIdx6, lpMPC_dlubcc6);
lpMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_4(lpMPC_ccrhsl7, lpMPC_slb7, lpMPC_llbbyslb7, lpMPC_dzcc7, lpMPC_lbIdx7, lpMPC_dllbcc7);
lpMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_4(lpMPC_ccrhsub7, lpMPC_sub7, lpMPC_lubbysub7, lpMPC_dzcc7, lpMPC_ubIdx7, lpMPC_dlubcc7);
lpMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_4(lpMPC_ccrhsl8, lpMPC_slb8, lpMPC_llbbyslb8, lpMPC_dzcc8, lpMPC_lbIdx8, lpMPC_dllbcc8);
lpMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_4(lpMPC_ccrhsub8, lpMPC_sub8, lpMPC_lubbysub8, lpMPC_dzcc8, lpMPC_ubIdx8, lpMPC_dlubcc8);
lpMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_2(lpMPC_ccrhsl9, lpMPC_slb9, lpMPC_llbbyslb9, lpMPC_dzcc9, lpMPC_lbIdx9, lpMPC_dllbcc9);
lpMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_2(lpMPC_ccrhsub9, lpMPC_sub9, lpMPC_lubbysub9, lpMPC_dzcc9, lpMPC_ubIdx9, lpMPC_dlubcc9);
lpMPC_LA_VSUB7_76(lpMPC_l, lpMPC_ccrhs, lpMPC_s, lpMPC_dl_cc, lpMPC_ds_cc);
lpMPC_LA_VADD_38(lpMPC_dz_cc, lpMPC_dz_aff);
lpMPC_LA_VADD_20(lpMPC_dv_cc, lpMPC_dv_aff);
lpMPC_LA_VADD_76(lpMPC_dl_cc, lpMPC_dl_aff);
lpMPC_LA_VADD_76(lpMPC_ds_cc, lpMPC_ds_aff);
info->lsit_cc = lpMPC_LINESEARCH_BACKTRACKING_COMBINED(lpMPC_z, lpMPC_v, lpMPC_l, lpMPC_s, lpMPC_dz_cc, lpMPC_dv_cc, lpMPC_dl_cc, lpMPC_ds_cc, &info->step_cc, &info->mu);
if( info->lsit_cc == lpMPC_NOPROGRESS ){
exitcode = lpMPC_NOPROGRESS; break;
}
info->it++;
}
output->z1[0] = lpMPC_z0[0];
output->z1[1] = lpMPC_z0[1];
output->z1[2] = lpMPC_z0[2];
output->z1[3] = lpMPC_z0[3];
output->z2[0] = lpMPC_z1[0];
output->z2[1] = lpMPC_z1[1];
output->z2[2] = lpMPC_z1[2];
output->z2[3] = lpMPC_z1[3];
output->z3[0] = lpMPC_z2[0];
output->z3[1] = lpMPC_z2[1];
output->z3[2] = lpMPC_z2[2];
output->z3[3] = lpMPC_z2[3];
output->z4[0] = lpMPC_z3[0];
output->z4[1] = lpMPC_z3[1];
output->z4[2] = lpMPC_z3[2];
output->z4[3] = lpMPC_z3[3];
output->z5[0] = lpMPC_z4[0];
output->z5[1] = lpMPC_z4[1];
output->z5[2] = lpMPC_z4[2];
output->z5[3] = lpMPC_z4[3];
output->z6[0] = lpMPC_z5[0];
output->z6[1] = lpMPC_z5[1];
output->z6[2] = lpMPC_z5[2];
output->z6[3] = lpMPC_z5[3];
output->z7[0] = lpMPC_z6[0];
output->z7[1] = lpMPC_z6[1];
output->z7[2] = lpMPC_z6[2];
output->z7[3] = lpMPC_z6[3];
output->z8[0] = lpMPC_z7[0];
output->z8[1] = lpMPC_z7[1];
output->z8[2] = lpMPC_z7[2];
output->z8[3] = lpMPC_z7[3];
output->z9[0] = lpMPC_z8[0];
output->z9[1] = lpMPC_z8[1];
output->z9[2] = lpMPC_z8[2];
output->z9[3] = lpMPC_z8[3];
output->z10[0] = lpMPC_z9[0];
output->z10[1] = lpMPC_z9[1];

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
