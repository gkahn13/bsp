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

#include "../include/stateMPC.h"

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
void stateMPC_LA_INITIALIZEVECTOR_58(stateMPC_FLOAT* vec, stateMPC_FLOAT value)
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
void stateMPC_LA_INITIALIZEVECTOR_30(stateMPC_FLOAT* vec, stateMPC_FLOAT value)
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
void stateMPC_LA_INITIALIZEVECTOR_116(stateMPC_FLOAT* vec, stateMPC_FLOAT value)
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
void stateMPC_LA_DOTACC_116(stateMPC_FLOAT *x, stateMPC_FLOAT *y, stateMPC_FLOAT *z)
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
void stateMPC_LA_DIAG_QUADFCN_4(stateMPC_FLOAT* H, stateMPC_FLOAT* f, stateMPC_FLOAT* z, stateMPC_FLOAT* grad, stateMPC_FLOAT* value)
{
	int i;
	stateMPC_FLOAT hz;	
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
void stateMPC_LA_DIAG_QUADFCN_2(stateMPC_FLOAT* H, stateMPC_FLOAT* f, stateMPC_FLOAT* z, stateMPC_FLOAT* grad, stateMPC_FLOAT* value)
{
	int i;
	stateMPC_FLOAT hz;	
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
void stateMPC_LA_DIAGZERO_MVMSUB6_2(stateMPC_FLOAT *B, stateMPC_FLOAT *u, stateMPC_FLOAT *b, stateMPC_FLOAT *l, stateMPC_FLOAT *r, stateMPC_FLOAT *z, stateMPC_FLOAT *y)
{
	int i;
	stateMPC_FLOAT Bu[2];
	stateMPC_FLOAT norm = *y;
	stateMPC_FLOAT lr = 0;

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
void stateMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_4_4(stateMPC_FLOAT *A, stateMPC_FLOAT *x, stateMPC_FLOAT *B, stateMPC_FLOAT *u, stateMPC_FLOAT *b, stateMPC_FLOAT *l, stateMPC_FLOAT *r, stateMPC_FLOAT *z, stateMPC_FLOAT *y)
{
	int i;
	int j;
	int k = 0;
	stateMPC_FLOAT AxBu[2];
	stateMPC_FLOAT norm = *y;
	stateMPC_FLOAT lr = 0;

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
void stateMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_4_2(stateMPC_FLOAT *A, stateMPC_FLOAT *x, stateMPC_FLOAT *B, stateMPC_FLOAT *u, stateMPC_FLOAT *b, stateMPC_FLOAT *l, stateMPC_FLOAT *r, stateMPC_FLOAT *z, stateMPC_FLOAT *y)
{
	int i;
	int j;
	int k = 0;
	stateMPC_FLOAT AxBu[2];
	stateMPC_FLOAT norm = *y;
	stateMPC_FLOAT lr = 0;

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
void stateMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(stateMPC_FLOAT *A, stateMPC_FLOAT *x, stateMPC_FLOAT *B, stateMPC_FLOAT *y, stateMPC_FLOAT *z)
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
void stateMPC_LA_DIAGZERO_MTVM_2_2(stateMPC_FLOAT *M, stateMPC_FLOAT *x, stateMPC_FLOAT *y)
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
void stateMPC_LA_VSUBADD3_4(stateMPC_FLOAT* t, stateMPC_FLOAT* u, int* uidx, stateMPC_FLOAT* v, stateMPC_FLOAT* w, stateMPC_FLOAT* y, stateMPC_FLOAT* z, stateMPC_FLOAT* r)
{
	int i;
	stateMPC_FLOAT norm = *r;
	stateMPC_FLOAT vx = 0;
	stateMPC_FLOAT x;
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
void stateMPC_LA_VSUBADD2_4(stateMPC_FLOAT* t, int* tidx, stateMPC_FLOAT* u, stateMPC_FLOAT* v, stateMPC_FLOAT* w, stateMPC_FLOAT* y, stateMPC_FLOAT* z, stateMPC_FLOAT* r)
{
	int i;
	stateMPC_FLOAT norm = *r;
	stateMPC_FLOAT vx = 0;
	stateMPC_FLOAT x;
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
void stateMPC_LA_VSUBADD3_2(stateMPC_FLOAT* t, stateMPC_FLOAT* u, int* uidx, stateMPC_FLOAT* v, stateMPC_FLOAT* w, stateMPC_FLOAT* y, stateMPC_FLOAT* z, stateMPC_FLOAT* r)
{
	int i;
	stateMPC_FLOAT norm = *r;
	stateMPC_FLOAT vx = 0;
	stateMPC_FLOAT x;
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
void stateMPC_LA_VSUBADD2_2(stateMPC_FLOAT* t, int* tidx, stateMPC_FLOAT* u, stateMPC_FLOAT* v, stateMPC_FLOAT* w, stateMPC_FLOAT* y, stateMPC_FLOAT* z, stateMPC_FLOAT* r)
{
	int i;
	stateMPC_FLOAT norm = *r;
	stateMPC_FLOAT vx = 0;
	stateMPC_FLOAT x;
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
void stateMPC_LA_INEQ_B_GRAD_4_4_4(stateMPC_FLOAT *lu, stateMPC_FLOAT *su, stateMPC_FLOAT *ru, stateMPC_FLOAT *ll, stateMPC_FLOAT *sl, stateMPC_FLOAT *rl, int* lbIdx, int* ubIdx, stateMPC_FLOAT *grad, stateMPC_FLOAT *lubysu, stateMPC_FLOAT *llbysl)
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
void stateMPC_LA_INEQ_B_GRAD_2_2_2(stateMPC_FLOAT *lu, stateMPC_FLOAT *su, stateMPC_FLOAT *ru, stateMPC_FLOAT *ll, stateMPC_FLOAT *sl, stateMPC_FLOAT *rl, int* lbIdx, int* ubIdx, stateMPC_FLOAT *grad, stateMPC_FLOAT *lubysu, stateMPC_FLOAT *llbysl)
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
void stateMPC_LA_VVADD3_58(stateMPC_FLOAT *u, stateMPC_FLOAT *v, stateMPC_FLOAT *w, stateMPC_FLOAT *z)
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
void stateMPC_LA_DIAG_CHOL_ONELOOP_LBUB_4_4_4(stateMPC_FLOAT *H, stateMPC_FLOAT *llbysl, int* lbIdx, stateMPC_FLOAT *lubysu, int* ubIdx, stateMPC_FLOAT *Phi)


{
	int i;
	
	/* compute cholesky */
	for( i=0; i<4; i++ ){
		Phi[i] = H[i] + llbysl[i] + lubysu[i];

#if stateMPC_SET_PRINTLEVEL > 0 && defined PRINTNUMERICALWARNINGS
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
void stateMPC_LA_DIAG_MATRIXFORWARDSUB_2_4(stateMPC_FLOAT *L, stateMPC_FLOAT *B, stateMPC_FLOAT *A)
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
void stateMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_4(stateMPC_FLOAT *L, stateMPC_FLOAT *B, stateMPC_FLOAT *A)
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
void stateMPC_LA_DENSE_DIAGZERO_MMTM_2_4_2(stateMPC_FLOAT *A, stateMPC_FLOAT *B, stateMPC_FLOAT *C)
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
void stateMPC_LA_DIAG_FORWARDSUB_4(stateMPC_FLOAT *L, stateMPC_FLOAT *b, stateMPC_FLOAT *y)
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
void stateMPC_LA_DIAG_CHOL_ONELOOP_LBUB_2_2_2(stateMPC_FLOAT *H, stateMPC_FLOAT *llbysl, int* lbIdx, stateMPC_FLOAT *lubysu, int* ubIdx, stateMPC_FLOAT *Phi)


{
	int i;
	
	/* compute cholesky */
	for( i=0; i<2; i++ ){
		Phi[i] = H[i] + llbysl[i] + lubysu[i];

#if stateMPC_SET_PRINTLEVEL > 0 && defined PRINTNUMERICALWARNINGS
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
void stateMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_2(stateMPC_FLOAT *L, stateMPC_FLOAT *B, stateMPC_FLOAT *A)
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
void stateMPC_LA_DIAG_FORWARDSUB_2(stateMPC_FLOAT *L, stateMPC_FLOAT *b, stateMPC_FLOAT *y)
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
void stateMPC_LA_DIAGZERO_MMT_2(stateMPC_FLOAT *B, stateMPC_FLOAT *L)
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
void stateMPC_LA_DIAGZERO_MVMSUB7_2(stateMPC_FLOAT *B, stateMPC_FLOAT *u, stateMPC_FLOAT *b, stateMPC_FLOAT *r)
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
void stateMPC_LA_DENSE_DIAGZERO_MMT2_2_4_4(stateMPC_FLOAT *A, stateMPC_FLOAT *B, stateMPC_FLOAT *L)
{
    int i, j, k, ii, di;
    stateMPC_FLOAT ltemp;
    
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
void stateMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_4_4(stateMPC_FLOAT *A, stateMPC_FLOAT *x, stateMPC_FLOAT *B, stateMPC_FLOAT *u, stateMPC_FLOAT *b, stateMPC_FLOAT *r)
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
void stateMPC_LA_DENSE_DIAGZERO_MMT2_2_4_2(stateMPC_FLOAT *A, stateMPC_FLOAT *B, stateMPC_FLOAT *L)
{
    int i, j, k, ii, di;
    stateMPC_FLOAT ltemp;
    
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
void stateMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_4_2(stateMPC_FLOAT *A, stateMPC_FLOAT *x, stateMPC_FLOAT *B, stateMPC_FLOAT *u, stateMPC_FLOAT *b, stateMPC_FLOAT *r)
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
void stateMPC_LA_DENSE_CHOL_2(stateMPC_FLOAT *A, stateMPC_FLOAT *L)
{
    int i, j, k, di, dj;
	 int ii, jj;

    stateMPC_FLOAT l;
    stateMPC_FLOAT Mii;

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
        
#if stateMPC_SET_PRINTLEVEL > 0 && defined PRINTNUMERICALWARNINGS
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
void stateMPC_LA_DENSE_FORWARDSUB_2(stateMPC_FLOAT *L, stateMPC_FLOAT *b, stateMPC_FLOAT *y)
{
    int i,j,ii,di;
    stateMPC_FLOAT yel;
            
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
void stateMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(stateMPC_FLOAT *L, stateMPC_FLOAT *B, stateMPC_FLOAT *A)
{
    int i,j,k,ii,di;
    stateMPC_FLOAT a;
    
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
void stateMPC_LA_DENSE_MMTSUB_2_2(stateMPC_FLOAT *A, stateMPC_FLOAT *L)
{
    int i, j, k, ii, di;
    stateMPC_FLOAT ltemp;
    
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
void stateMPC_LA_DENSE_MVMSUB1_2_2(stateMPC_FLOAT *A, stateMPC_FLOAT *x, stateMPC_FLOAT *b, stateMPC_FLOAT *r)
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
void stateMPC_LA_DENSE_BACKWARDSUB_2(stateMPC_FLOAT *L, stateMPC_FLOAT *y, stateMPC_FLOAT *x)
{
    int i, ii, di, j, jj, dj;
    stateMPC_FLOAT xel;    
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
void stateMPC_LA_DENSE_MTVMSUB_2_2(stateMPC_FLOAT *A, stateMPC_FLOAT *x, stateMPC_FLOAT *b, stateMPC_FLOAT *r)
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
void stateMPC_LA_VSUB2_58(stateMPC_FLOAT *x, stateMPC_FLOAT *y, stateMPC_FLOAT *z)
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
void stateMPC_LA_DIAG_FORWARDBACKWARDSUB_4(stateMPC_FLOAT *L, stateMPC_FLOAT *b, stateMPC_FLOAT *x)
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
void stateMPC_LA_DIAG_FORWARDBACKWARDSUB_2(stateMPC_FLOAT *L, stateMPC_FLOAT *b, stateMPC_FLOAT *x)
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
void stateMPC_LA_VSUB_INDEXED_4(stateMPC_FLOAT *x, int* xidx, stateMPC_FLOAT *y, stateMPC_FLOAT *z)
{
	int i;
	for( i=0; i<4; i++){
		z[i] = x[xidx[i]] - y[i];
	}
}


/*
 * Vector subtraction x = -u.*v - w for vectors of length 4.
 */
void stateMPC_LA_VSUB3_4(stateMPC_FLOAT *u, stateMPC_FLOAT *v, stateMPC_FLOAT *w, stateMPC_FLOAT *x)
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
void stateMPC_LA_VSUB2_INDEXED_4(stateMPC_FLOAT *x, stateMPC_FLOAT *y, int* yidx, stateMPC_FLOAT *z)
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
void stateMPC_LA_VSUB_INDEXED_2(stateMPC_FLOAT *x, int* xidx, stateMPC_FLOAT *y, stateMPC_FLOAT *z)
{
	int i;
	for( i=0; i<2; i++){
		z[i] = x[xidx[i]] - y[i];
	}
}


/*
 * Vector subtraction x = -u.*v - w for vectors of length 2.
 */
void stateMPC_LA_VSUB3_2(stateMPC_FLOAT *u, stateMPC_FLOAT *v, stateMPC_FLOAT *w, stateMPC_FLOAT *x)
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
void stateMPC_LA_VSUB2_INDEXED_2(stateMPC_FLOAT *x, stateMPC_FLOAT *y, int* yidx, stateMPC_FLOAT *z)
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
 * stateMPC_NOPROGRESS (should be negative).
 */
int stateMPC_LINESEARCH_BACKTRACKING_AFFINE(stateMPC_FLOAT *l, stateMPC_FLOAT *s, stateMPC_FLOAT *dl, stateMPC_FLOAT *ds, stateMPC_FLOAT *a, stateMPC_FLOAT *mu_aff)
{
    int i;
	int lsIt=1;    
    stateMPC_FLOAT dltemp;
    stateMPC_FLOAT dstemp;
    stateMPC_FLOAT mya = 1.0;
    stateMPC_FLOAT mymu;
        
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
            mya *= stateMPC_SET_LS_SCALE_AFF;
            if( mya < stateMPC_SET_LS_MINSTEP ){
                return stateMPC_NOPROGRESS;
            }
        }
    }
    
    /* return new values and iteration counter */
    *a = mya;
    *mu_aff = mymu / (stateMPC_FLOAT)116;
    return lsIt;
}


/*
 * Vector subtraction x = u.*v - a where a is a scalar
*  and x,u,v are vectors of length 116.
 */
void stateMPC_LA_VSUB5_116(stateMPC_FLOAT *u, stateMPC_FLOAT *v, stateMPC_FLOAT a, stateMPC_FLOAT *x)
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
void stateMPC_LA_VSUB6_INDEXED_4_4_4(stateMPC_FLOAT *u, stateMPC_FLOAT *su, int* uidx, stateMPC_FLOAT *v, stateMPC_FLOAT *sv, int* vidx, stateMPC_FLOAT *x)
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
void stateMPC_LA_DIAGZERO_MVM_2(stateMPC_FLOAT *B, stateMPC_FLOAT *u, stateMPC_FLOAT *r)
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
void stateMPC_LA_DENSE_DIAGZERO_2MVMADD_2_4_4(stateMPC_FLOAT *A, stateMPC_FLOAT *x, stateMPC_FLOAT *B, stateMPC_FLOAT *u, stateMPC_FLOAT *r)
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
void stateMPC_LA_VSUB6_INDEXED_2_2_2(stateMPC_FLOAT *u, stateMPC_FLOAT *su, int* uidx, stateMPC_FLOAT *v, stateMPC_FLOAT *sv, int* vidx, stateMPC_FLOAT *x)
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
void stateMPC_LA_DENSE_DIAGZERO_2MVMADD_2_4_2(stateMPC_FLOAT *A, stateMPC_FLOAT *x, stateMPC_FLOAT *B, stateMPC_FLOAT *u, stateMPC_FLOAT *r)
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
void stateMPC_LA_VSUB_58(stateMPC_FLOAT *x, stateMPC_FLOAT *y, stateMPC_FLOAT *z)
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
void stateMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_4(stateMPC_FLOAT *r, stateMPC_FLOAT *s, stateMPC_FLOAT *u, stateMPC_FLOAT *y, int* yidx, stateMPC_FLOAT *z)
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
void stateMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_4(stateMPC_FLOAT *r, stateMPC_FLOAT *s, stateMPC_FLOAT *u, stateMPC_FLOAT *y, int* yidx, stateMPC_FLOAT *z)
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
void stateMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_2(stateMPC_FLOAT *r, stateMPC_FLOAT *s, stateMPC_FLOAT *u, stateMPC_FLOAT *y, int* yidx, stateMPC_FLOAT *z)
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
void stateMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_2(stateMPC_FLOAT *r, stateMPC_FLOAT *s, stateMPC_FLOAT *u, stateMPC_FLOAT *y, int* yidx, stateMPC_FLOAT *z)
{
	int i;
	for( i=0; i<2; i++ ){
		z[i] = -r[i]/s[i] + u[i]*y[yidx[i]];
	}
}


/*
 * Computes ds = -l.\(r + s.*dl) for vectors of length 116.
 */
void stateMPC_LA_VSUB7_116(stateMPC_FLOAT *l, stateMPC_FLOAT *r, stateMPC_FLOAT *s, stateMPC_FLOAT *dl, stateMPC_FLOAT *ds)
{
	int i;
	for( i=0; i<116; i++){
		ds[i] = -(r[i] + s[i]*dl[i])/l[i];
	}
}


/*
 * Vector addition x = x + y for vectors of length 58.
 */
void stateMPC_LA_VADD_58(stateMPC_FLOAT *x, stateMPC_FLOAT *y)
{
	int i;
	for( i=0; i<58; i++){
		x[i] += y[i];
	}
}


/*
 * Vector addition x = x + y for vectors of length 30.
 */
void stateMPC_LA_VADD_30(stateMPC_FLOAT *x, stateMPC_FLOAT *y)
{
	int i;
	for( i=0; i<30; i++){
		x[i] += y[i];
	}
}


/*
 * Vector addition x = x + y for vectors of length 116.
 */
void stateMPC_LA_VADD_116(stateMPC_FLOAT *x, stateMPC_FLOAT *y)
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
int stateMPC_LINESEARCH_BACKTRACKING_COMBINED(stateMPC_FLOAT *z, stateMPC_FLOAT *v, stateMPC_FLOAT *l, stateMPC_FLOAT *s, stateMPC_FLOAT *dz, stateMPC_FLOAT *dv, stateMPC_FLOAT *dl, stateMPC_FLOAT *ds, stateMPC_FLOAT *a, stateMPC_FLOAT *mu)
{
    int i, lsIt=1;       
    stateMPC_FLOAT dltemp;
    stateMPC_FLOAT dstemp;    
    stateMPC_FLOAT a_gamma;
            
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
            *a *= stateMPC_SET_LS_SCALE;
            if( *a < stateMPC_SET_LS_MINSTEP ){
                return stateMPC_NOPROGRESS;
            }
        }
    }
    
    /* update variables with safety margin */
    a_gamma = (*a)*stateMPC_SET_LS_MAXSTEP;
    
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
    *mu /= (stateMPC_FLOAT)116;
    return lsIt;
}




/* VARIABLE DEFINITIONS ------------------------------------------------ */
stateMPC_FLOAT stateMPC_z[58];
stateMPC_FLOAT stateMPC_v[30];
stateMPC_FLOAT stateMPC_dz_aff[58];
stateMPC_FLOAT stateMPC_dv_aff[30];
stateMPC_FLOAT stateMPC_grad_cost[58];
stateMPC_FLOAT stateMPC_grad_eq[58];
stateMPC_FLOAT stateMPC_rd[58];
stateMPC_FLOAT stateMPC_l[116];
stateMPC_FLOAT stateMPC_s[116];
stateMPC_FLOAT stateMPC_lbys[116];
stateMPC_FLOAT stateMPC_dl_aff[116];
stateMPC_FLOAT stateMPC_ds_aff[116];
stateMPC_FLOAT stateMPC_dz_cc[58];
stateMPC_FLOAT stateMPC_dv_cc[30];
stateMPC_FLOAT stateMPC_dl_cc[116];
stateMPC_FLOAT stateMPC_ds_cc[116];
stateMPC_FLOAT stateMPC_ccrhs[116];
stateMPC_FLOAT stateMPC_grad_ineq[58];
stateMPC_FLOAT* stateMPC_z00 = stateMPC_z + 0;
stateMPC_FLOAT* stateMPC_dzaff00 = stateMPC_dz_aff + 0;
stateMPC_FLOAT* stateMPC_dzcc00 = stateMPC_dz_cc + 0;
stateMPC_FLOAT* stateMPC_rd00 = stateMPC_rd + 0;
stateMPC_FLOAT stateMPC_Lbyrd00[4];
stateMPC_FLOAT* stateMPC_grad_cost00 = stateMPC_grad_cost + 0;
stateMPC_FLOAT* stateMPC_grad_eq00 = stateMPC_grad_eq + 0;
stateMPC_FLOAT* stateMPC_grad_ineq00 = stateMPC_grad_ineq + 0;
stateMPC_FLOAT stateMPC_ctv00[4];
stateMPC_FLOAT* stateMPC_v00 = stateMPC_v + 0;
stateMPC_FLOAT stateMPC_re00[2];
stateMPC_FLOAT stateMPC_beta00[2];
stateMPC_FLOAT stateMPC_betacc00[2];
stateMPC_FLOAT* stateMPC_dvaff00 = stateMPC_dv_aff + 0;
stateMPC_FLOAT* stateMPC_dvcc00 = stateMPC_dv_cc + 0;
stateMPC_FLOAT stateMPC_V00[8];
stateMPC_FLOAT stateMPC_Yd00[3];
stateMPC_FLOAT stateMPC_Ld00[3];
stateMPC_FLOAT stateMPC_yy00[2];
stateMPC_FLOAT stateMPC_bmy00[2];
int stateMPC_lbIdx00[4] = {0, 1, 2, 3};
stateMPC_FLOAT* stateMPC_llb00 = stateMPC_l + 0;
stateMPC_FLOAT* stateMPC_slb00 = stateMPC_s + 0;
stateMPC_FLOAT* stateMPC_llbbyslb00 = stateMPC_lbys + 0;
stateMPC_FLOAT stateMPC_rilb00[4];
stateMPC_FLOAT* stateMPC_dllbaff00 = stateMPC_dl_aff + 0;
stateMPC_FLOAT* stateMPC_dslbaff00 = stateMPC_ds_aff + 0;
stateMPC_FLOAT* stateMPC_dllbcc00 = stateMPC_dl_cc + 0;
stateMPC_FLOAT* stateMPC_dslbcc00 = stateMPC_ds_cc + 0;
stateMPC_FLOAT* stateMPC_ccrhsl00 = stateMPC_ccrhs + 0;
int stateMPC_ubIdx00[4] = {0, 1, 2, 3};
stateMPC_FLOAT* stateMPC_lub00 = stateMPC_l + 4;
stateMPC_FLOAT* stateMPC_sub00 = stateMPC_s + 4;
stateMPC_FLOAT* stateMPC_lubbysub00 = stateMPC_lbys + 4;
stateMPC_FLOAT stateMPC_riub00[4];
stateMPC_FLOAT* stateMPC_dlubaff00 = stateMPC_dl_aff + 4;
stateMPC_FLOAT* stateMPC_dsubaff00 = stateMPC_ds_aff + 4;
stateMPC_FLOAT* stateMPC_dlubcc00 = stateMPC_dl_cc + 4;
stateMPC_FLOAT* stateMPC_dsubcc00 = stateMPC_ds_cc + 4;
stateMPC_FLOAT* stateMPC_ccrhsub00 = stateMPC_ccrhs + 4;
stateMPC_FLOAT stateMPC_Phi00[4];
stateMPC_FLOAT stateMPC_D00[4] = {1.0000000000000000E+000, 
1.0000000000000000E+000};
stateMPC_FLOAT stateMPC_W00[4];
stateMPC_FLOAT* stateMPC_z01 = stateMPC_z + 4;
stateMPC_FLOAT* stateMPC_dzaff01 = stateMPC_dz_aff + 4;
stateMPC_FLOAT* stateMPC_dzcc01 = stateMPC_dz_cc + 4;
stateMPC_FLOAT* stateMPC_rd01 = stateMPC_rd + 4;
stateMPC_FLOAT stateMPC_Lbyrd01[4];
stateMPC_FLOAT* stateMPC_grad_cost01 = stateMPC_grad_cost + 4;
stateMPC_FLOAT* stateMPC_grad_eq01 = stateMPC_grad_eq + 4;
stateMPC_FLOAT* stateMPC_grad_ineq01 = stateMPC_grad_ineq + 4;
stateMPC_FLOAT stateMPC_ctv01[4];
stateMPC_FLOAT* stateMPC_v01 = stateMPC_v + 2;
stateMPC_FLOAT stateMPC_re01[2];
stateMPC_FLOAT stateMPC_beta01[2];
stateMPC_FLOAT stateMPC_betacc01[2];
stateMPC_FLOAT* stateMPC_dvaff01 = stateMPC_dv_aff + 2;
stateMPC_FLOAT* stateMPC_dvcc01 = stateMPC_dv_cc + 2;
stateMPC_FLOAT stateMPC_V01[8];
stateMPC_FLOAT stateMPC_Yd01[3];
stateMPC_FLOAT stateMPC_Ld01[3];
stateMPC_FLOAT stateMPC_yy01[2];
stateMPC_FLOAT stateMPC_bmy01[2];
int stateMPC_lbIdx01[4] = {0, 1, 2, 3};
stateMPC_FLOAT* stateMPC_llb01 = stateMPC_l + 8;
stateMPC_FLOAT* stateMPC_slb01 = stateMPC_s + 8;
stateMPC_FLOAT* stateMPC_llbbyslb01 = stateMPC_lbys + 8;
stateMPC_FLOAT stateMPC_rilb01[4];
stateMPC_FLOAT* stateMPC_dllbaff01 = stateMPC_dl_aff + 8;
stateMPC_FLOAT* stateMPC_dslbaff01 = stateMPC_ds_aff + 8;
stateMPC_FLOAT* stateMPC_dllbcc01 = stateMPC_dl_cc + 8;
stateMPC_FLOAT* stateMPC_dslbcc01 = stateMPC_ds_cc + 8;
stateMPC_FLOAT* stateMPC_ccrhsl01 = stateMPC_ccrhs + 8;
int stateMPC_ubIdx01[4] = {0, 1, 2, 3};
stateMPC_FLOAT* stateMPC_lub01 = stateMPC_l + 12;
stateMPC_FLOAT* stateMPC_sub01 = stateMPC_s + 12;
stateMPC_FLOAT* stateMPC_lubbysub01 = stateMPC_lbys + 12;
stateMPC_FLOAT stateMPC_riub01[4];
stateMPC_FLOAT* stateMPC_dlubaff01 = stateMPC_dl_aff + 12;
stateMPC_FLOAT* stateMPC_dsubaff01 = stateMPC_ds_aff + 12;
stateMPC_FLOAT* stateMPC_dlubcc01 = stateMPC_dl_cc + 12;
stateMPC_FLOAT* stateMPC_dsubcc01 = stateMPC_ds_cc + 12;
stateMPC_FLOAT* stateMPC_ccrhsub01 = stateMPC_ccrhs + 12;
stateMPC_FLOAT stateMPC_Phi01[4];
stateMPC_FLOAT stateMPC_D01[4] = {-1.0000000000000000E+000, 
-1.0000000000000000E+000};
stateMPC_FLOAT stateMPC_W01[4];
stateMPC_FLOAT stateMPC_Ysd01[4];
stateMPC_FLOAT stateMPC_Lsd01[4];
stateMPC_FLOAT* stateMPC_z02 = stateMPC_z + 8;
stateMPC_FLOAT* stateMPC_dzaff02 = stateMPC_dz_aff + 8;
stateMPC_FLOAT* stateMPC_dzcc02 = stateMPC_dz_cc + 8;
stateMPC_FLOAT* stateMPC_rd02 = stateMPC_rd + 8;
stateMPC_FLOAT stateMPC_Lbyrd02[4];
stateMPC_FLOAT* stateMPC_grad_cost02 = stateMPC_grad_cost + 8;
stateMPC_FLOAT* stateMPC_grad_eq02 = stateMPC_grad_eq + 8;
stateMPC_FLOAT* stateMPC_grad_ineq02 = stateMPC_grad_ineq + 8;
stateMPC_FLOAT stateMPC_ctv02[4];
stateMPC_FLOAT* stateMPC_v02 = stateMPC_v + 4;
stateMPC_FLOAT stateMPC_re02[2];
stateMPC_FLOAT stateMPC_beta02[2];
stateMPC_FLOAT stateMPC_betacc02[2];
stateMPC_FLOAT* stateMPC_dvaff02 = stateMPC_dv_aff + 4;
stateMPC_FLOAT* stateMPC_dvcc02 = stateMPC_dv_cc + 4;
stateMPC_FLOAT stateMPC_V02[8];
stateMPC_FLOAT stateMPC_Yd02[3];
stateMPC_FLOAT stateMPC_Ld02[3];
stateMPC_FLOAT stateMPC_yy02[2];
stateMPC_FLOAT stateMPC_bmy02[2];
int stateMPC_lbIdx02[4] = {0, 1, 2, 3};
stateMPC_FLOAT* stateMPC_llb02 = stateMPC_l + 16;
stateMPC_FLOAT* stateMPC_slb02 = stateMPC_s + 16;
stateMPC_FLOAT* stateMPC_llbbyslb02 = stateMPC_lbys + 16;
stateMPC_FLOAT stateMPC_rilb02[4];
stateMPC_FLOAT* stateMPC_dllbaff02 = stateMPC_dl_aff + 16;
stateMPC_FLOAT* stateMPC_dslbaff02 = stateMPC_ds_aff + 16;
stateMPC_FLOAT* stateMPC_dllbcc02 = stateMPC_dl_cc + 16;
stateMPC_FLOAT* stateMPC_dslbcc02 = stateMPC_ds_cc + 16;
stateMPC_FLOAT* stateMPC_ccrhsl02 = stateMPC_ccrhs + 16;
int stateMPC_ubIdx02[4] = {0, 1, 2, 3};
stateMPC_FLOAT* stateMPC_lub02 = stateMPC_l + 20;
stateMPC_FLOAT* stateMPC_sub02 = stateMPC_s + 20;
stateMPC_FLOAT* stateMPC_lubbysub02 = stateMPC_lbys + 20;
stateMPC_FLOAT stateMPC_riub02[4];
stateMPC_FLOAT* stateMPC_dlubaff02 = stateMPC_dl_aff + 20;
stateMPC_FLOAT* stateMPC_dsubaff02 = stateMPC_ds_aff + 20;
stateMPC_FLOAT* stateMPC_dlubcc02 = stateMPC_dl_cc + 20;
stateMPC_FLOAT* stateMPC_dsubcc02 = stateMPC_ds_cc + 20;
stateMPC_FLOAT* stateMPC_ccrhsub02 = stateMPC_ccrhs + 20;
stateMPC_FLOAT stateMPC_Phi02[4];
stateMPC_FLOAT stateMPC_W02[4];
stateMPC_FLOAT stateMPC_Ysd02[4];
stateMPC_FLOAT stateMPC_Lsd02[4];
stateMPC_FLOAT* stateMPC_z03 = stateMPC_z + 12;
stateMPC_FLOAT* stateMPC_dzaff03 = stateMPC_dz_aff + 12;
stateMPC_FLOAT* stateMPC_dzcc03 = stateMPC_dz_cc + 12;
stateMPC_FLOAT* stateMPC_rd03 = stateMPC_rd + 12;
stateMPC_FLOAT stateMPC_Lbyrd03[4];
stateMPC_FLOAT* stateMPC_grad_cost03 = stateMPC_grad_cost + 12;
stateMPC_FLOAT* stateMPC_grad_eq03 = stateMPC_grad_eq + 12;
stateMPC_FLOAT* stateMPC_grad_ineq03 = stateMPC_grad_ineq + 12;
stateMPC_FLOAT stateMPC_ctv03[4];
stateMPC_FLOAT* stateMPC_v03 = stateMPC_v + 6;
stateMPC_FLOAT stateMPC_re03[2];
stateMPC_FLOAT stateMPC_beta03[2];
stateMPC_FLOAT stateMPC_betacc03[2];
stateMPC_FLOAT* stateMPC_dvaff03 = stateMPC_dv_aff + 6;
stateMPC_FLOAT* stateMPC_dvcc03 = stateMPC_dv_cc + 6;
stateMPC_FLOAT stateMPC_V03[8];
stateMPC_FLOAT stateMPC_Yd03[3];
stateMPC_FLOAT stateMPC_Ld03[3];
stateMPC_FLOAT stateMPC_yy03[2];
stateMPC_FLOAT stateMPC_bmy03[2];
int stateMPC_lbIdx03[4] = {0, 1, 2, 3};
stateMPC_FLOAT* stateMPC_llb03 = stateMPC_l + 24;
stateMPC_FLOAT* stateMPC_slb03 = stateMPC_s + 24;
stateMPC_FLOAT* stateMPC_llbbyslb03 = stateMPC_lbys + 24;
stateMPC_FLOAT stateMPC_rilb03[4];
stateMPC_FLOAT* stateMPC_dllbaff03 = stateMPC_dl_aff + 24;
stateMPC_FLOAT* stateMPC_dslbaff03 = stateMPC_ds_aff + 24;
stateMPC_FLOAT* stateMPC_dllbcc03 = stateMPC_dl_cc + 24;
stateMPC_FLOAT* stateMPC_dslbcc03 = stateMPC_ds_cc + 24;
stateMPC_FLOAT* stateMPC_ccrhsl03 = stateMPC_ccrhs + 24;
int stateMPC_ubIdx03[4] = {0, 1, 2, 3};
stateMPC_FLOAT* stateMPC_lub03 = stateMPC_l + 28;
stateMPC_FLOAT* stateMPC_sub03 = stateMPC_s + 28;
stateMPC_FLOAT* stateMPC_lubbysub03 = stateMPC_lbys + 28;
stateMPC_FLOAT stateMPC_riub03[4];
stateMPC_FLOAT* stateMPC_dlubaff03 = stateMPC_dl_aff + 28;
stateMPC_FLOAT* stateMPC_dsubaff03 = stateMPC_ds_aff + 28;
stateMPC_FLOAT* stateMPC_dlubcc03 = stateMPC_dl_cc + 28;
stateMPC_FLOAT* stateMPC_dsubcc03 = stateMPC_ds_cc + 28;
stateMPC_FLOAT* stateMPC_ccrhsub03 = stateMPC_ccrhs + 28;
stateMPC_FLOAT stateMPC_Phi03[4];
stateMPC_FLOAT stateMPC_W03[4];
stateMPC_FLOAT stateMPC_Ysd03[4];
stateMPC_FLOAT stateMPC_Lsd03[4];
stateMPC_FLOAT* stateMPC_z04 = stateMPC_z + 16;
stateMPC_FLOAT* stateMPC_dzaff04 = stateMPC_dz_aff + 16;
stateMPC_FLOAT* stateMPC_dzcc04 = stateMPC_dz_cc + 16;
stateMPC_FLOAT* stateMPC_rd04 = stateMPC_rd + 16;
stateMPC_FLOAT stateMPC_Lbyrd04[4];
stateMPC_FLOAT* stateMPC_grad_cost04 = stateMPC_grad_cost + 16;
stateMPC_FLOAT* stateMPC_grad_eq04 = stateMPC_grad_eq + 16;
stateMPC_FLOAT* stateMPC_grad_ineq04 = stateMPC_grad_ineq + 16;
stateMPC_FLOAT stateMPC_ctv04[4];
stateMPC_FLOAT* stateMPC_v04 = stateMPC_v + 8;
stateMPC_FLOAT stateMPC_re04[2];
stateMPC_FLOAT stateMPC_beta04[2];
stateMPC_FLOAT stateMPC_betacc04[2];
stateMPC_FLOAT* stateMPC_dvaff04 = stateMPC_dv_aff + 8;
stateMPC_FLOAT* stateMPC_dvcc04 = stateMPC_dv_cc + 8;
stateMPC_FLOAT stateMPC_V04[8];
stateMPC_FLOAT stateMPC_Yd04[3];
stateMPC_FLOAT stateMPC_Ld04[3];
stateMPC_FLOAT stateMPC_yy04[2];
stateMPC_FLOAT stateMPC_bmy04[2];
int stateMPC_lbIdx04[4] = {0, 1, 2, 3};
stateMPC_FLOAT* stateMPC_llb04 = stateMPC_l + 32;
stateMPC_FLOAT* stateMPC_slb04 = stateMPC_s + 32;
stateMPC_FLOAT* stateMPC_llbbyslb04 = stateMPC_lbys + 32;
stateMPC_FLOAT stateMPC_rilb04[4];
stateMPC_FLOAT* stateMPC_dllbaff04 = stateMPC_dl_aff + 32;
stateMPC_FLOAT* stateMPC_dslbaff04 = stateMPC_ds_aff + 32;
stateMPC_FLOAT* stateMPC_dllbcc04 = stateMPC_dl_cc + 32;
stateMPC_FLOAT* stateMPC_dslbcc04 = stateMPC_ds_cc + 32;
stateMPC_FLOAT* stateMPC_ccrhsl04 = stateMPC_ccrhs + 32;
int stateMPC_ubIdx04[4] = {0, 1, 2, 3};
stateMPC_FLOAT* stateMPC_lub04 = stateMPC_l + 36;
stateMPC_FLOAT* stateMPC_sub04 = stateMPC_s + 36;
stateMPC_FLOAT* stateMPC_lubbysub04 = stateMPC_lbys + 36;
stateMPC_FLOAT stateMPC_riub04[4];
stateMPC_FLOAT* stateMPC_dlubaff04 = stateMPC_dl_aff + 36;
stateMPC_FLOAT* stateMPC_dsubaff04 = stateMPC_ds_aff + 36;
stateMPC_FLOAT* stateMPC_dlubcc04 = stateMPC_dl_cc + 36;
stateMPC_FLOAT* stateMPC_dsubcc04 = stateMPC_ds_cc + 36;
stateMPC_FLOAT* stateMPC_ccrhsub04 = stateMPC_ccrhs + 36;
stateMPC_FLOAT stateMPC_Phi04[4];
stateMPC_FLOAT stateMPC_W04[4];
stateMPC_FLOAT stateMPC_Ysd04[4];
stateMPC_FLOAT stateMPC_Lsd04[4];
stateMPC_FLOAT* stateMPC_z05 = stateMPC_z + 20;
stateMPC_FLOAT* stateMPC_dzaff05 = stateMPC_dz_aff + 20;
stateMPC_FLOAT* stateMPC_dzcc05 = stateMPC_dz_cc + 20;
stateMPC_FLOAT* stateMPC_rd05 = stateMPC_rd + 20;
stateMPC_FLOAT stateMPC_Lbyrd05[4];
stateMPC_FLOAT* stateMPC_grad_cost05 = stateMPC_grad_cost + 20;
stateMPC_FLOAT* stateMPC_grad_eq05 = stateMPC_grad_eq + 20;
stateMPC_FLOAT* stateMPC_grad_ineq05 = stateMPC_grad_ineq + 20;
stateMPC_FLOAT stateMPC_ctv05[4];
stateMPC_FLOAT* stateMPC_v05 = stateMPC_v + 10;
stateMPC_FLOAT stateMPC_re05[2];
stateMPC_FLOAT stateMPC_beta05[2];
stateMPC_FLOAT stateMPC_betacc05[2];
stateMPC_FLOAT* stateMPC_dvaff05 = stateMPC_dv_aff + 10;
stateMPC_FLOAT* stateMPC_dvcc05 = stateMPC_dv_cc + 10;
stateMPC_FLOAT stateMPC_V05[8];
stateMPC_FLOAT stateMPC_Yd05[3];
stateMPC_FLOAT stateMPC_Ld05[3];
stateMPC_FLOAT stateMPC_yy05[2];
stateMPC_FLOAT stateMPC_bmy05[2];
int stateMPC_lbIdx05[4] = {0, 1, 2, 3};
stateMPC_FLOAT* stateMPC_llb05 = stateMPC_l + 40;
stateMPC_FLOAT* stateMPC_slb05 = stateMPC_s + 40;
stateMPC_FLOAT* stateMPC_llbbyslb05 = stateMPC_lbys + 40;
stateMPC_FLOAT stateMPC_rilb05[4];
stateMPC_FLOAT* stateMPC_dllbaff05 = stateMPC_dl_aff + 40;
stateMPC_FLOAT* stateMPC_dslbaff05 = stateMPC_ds_aff + 40;
stateMPC_FLOAT* stateMPC_dllbcc05 = stateMPC_dl_cc + 40;
stateMPC_FLOAT* stateMPC_dslbcc05 = stateMPC_ds_cc + 40;
stateMPC_FLOAT* stateMPC_ccrhsl05 = stateMPC_ccrhs + 40;
int stateMPC_ubIdx05[4] = {0, 1, 2, 3};
stateMPC_FLOAT* stateMPC_lub05 = stateMPC_l + 44;
stateMPC_FLOAT* stateMPC_sub05 = stateMPC_s + 44;
stateMPC_FLOAT* stateMPC_lubbysub05 = stateMPC_lbys + 44;
stateMPC_FLOAT stateMPC_riub05[4];
stateMPC_FLOAT* stateMPC_dlubaff05 = stateMPC_dl_aff + 44;
stateMPC_FLOAT* stateMPC_dsubaff05 = stateMPC_ds_aff + 44;
stateMPC_FLOAT* stateMPC_dlubcc05 = stateMPC_dl_cc + 44;
stateMPC_FLOAT* stateMPC_dsubcc05 = stateMPC_ds_cc + 44;
stateMPC_FLOAT* stateMPC_ccrhsub05 = stateMPC_ccrhs + 44;
stateMPC_FLOAT stateMPC_Phi05[4];
stateMPC_FLOAT stateMPC_W05[4];
stateMPC_FLOAT stateMPC_Ysd05[4];
stateMPC_FLOAT stateMPC_Lsd05[4];
stateMPC_FLOAT* stateMPC_z06 = stateMPC_z + 24;
stateMPC_FLOAT* stateMPC_dzaff06 = stateMPC_dz_aff + 24;
stateMPC_FLOAT* stateMPC_dzcc06 = stateMPC_dz_cc + 24;
stateMPC_FLOAT* stateMPC_rd06 = stateMPC_rd + 24;
stateMPC_FLOAT stateMPC_Lbyrd06[4];
stateMPC_FLOAT* stateMPC_grad_cost06 = stateMPC_grad_cost + 24;
stateMPC_FLOAT* stateMPC_grad_eq06 = stateMPC_grad_eq + 24;
stateMPC_FLOAT* stateMPC_grad_ineq06 = stateMPC_grad_ineq + 24;
stateMPC_FLOAT stateMPC_ctv06[4];
stateMPC_FLOAT* stateMPC_v06 = stateMPC_v + 12;
stateMPC_FLOAT stateMPC_re06[2];
stateMPC_FLOAT stateMPC_beta06[2];
stateMPC_FLOAT stateMPC_betacc06[2];
stateMPC_FLOAT* stateMPC_dvaff06 = stateMPC_dv_aff + 12;
stateMPC_FLOAT* stateMPC_dvcc06 = stateMPC_dv_cc + 12;
stateMPC_FLOAT stateMPC_V06[8];
stateMPC_FLOAT stateMPC_Yd06[3];
stateMPC_FLOAT stateMPC_Ld06[3];
stateMPC_FLOAT stateMPC_yy06[2];
stateMPC_FLOAT stateMPC_bmy06[2];
int stateMPC_lbIdx06[4] = {0, 1, 2, 3};
stateMPC_FLOAT* stateMPC_llb06 = stateMPC_l + 48;
stateMPC_FLOAT* stateMPC_slb06 = stateMPC_s + 48;
stateMPC_FLOAT* stateMPC_llbbyslb06 = stateMPC_lbys + 48;
stateMPC_FLOAT stateMPC_rilb06[4];
stateMPC_FLOAT* stateMPC_dllbaff06 = stateMPC_dl_aff + 48;
stateMPC_FLOAT* stateMPC_dslbaff06 = stateMPC_ds_aff + 48;
stateMPC_FLOAT* stateMPC_dllbcc06 = stateMPC_dl_cc + 48;
stateMPC_FLOAT* stateMPC_dslbcc06 = stateMPC_ds_cc + 48;
stateMPC_FLOAT* stateMPC_ccrhsl06 = stateMPC_ccrhs + 48;
int stateMPC_ubIdx06[4] = {0, 1, 2, 3};
stateMPC_FLOAT* stateMPC_lub06 = stateMPC_l + 52;
stateMPC_FLOAT* stateMPC_sub06 = stateMPC_s + 52;
stateMPC_FLOAT* stateMPC_lubbysub06 = stateMPC_lbys + 52;
stateMPC_FLOAT stateMPC_riub06[4];
stateMPC_FLOAT* stateMPC_dlubaff06 = stateMPC_dl_aff + 52;
stateMPC_FLOAT* stateMPC_dsubaff06 = stateMPC_ds_aff + 52;
stateMPC_FLOAT* stateMPC_dlubcc06 = stateMPC_dl_cc + 52;
stateMPC_FLOAT* stateMPC_dsubcc06 = stateMPC_ds_cc + 52;
stateMPC_FLOAT* stateMPC_ccrhsub06 = stateMPC_ccrhs + 52;
stateMPC_FLOAT stateMPC_Phi06[4];
stateMPC_FLOAT stateMPC_W06[4];
stateMPC_FLOAT stateMPC_Ysd06[4];
stateMPC_FLOAT stateMPC_Lsd06[4];
stateMPC_FLOAT* stateMPC_z07 = stateMPC_z + 28;
stateMPC_FLOAT* stateMPC_dzaff07 = stateMPC_dz_aff + 28;
stateMPC_FLOAT* stateMPC_dzcc07 = stateMPC_dz_cc + 28;
stateMPC_FLOAT* stateMPC_rd07 = stateMPC_rd + 28;
stateMPC_FLOAT stateMPC_Lbyrd07[4];
stateMPC_FLOAT* stateMPC_grad_cost07 = stateMPC_grad_cost + 28;
stateMPC_FLOAT* stateMPC_grad_eq07 = stateMPC_grad_eq + 28;
stateMPC_FLOAT* stateMPC_grad_ineq07 = stateMPC_grad_ineq + 28;
stateMPC_FLOAT stateMPC_ctv07[4];
stateMPC_FLOAT* stateMPC_v07 = stateMPC_v + 14;
stateMPC_FLOAT stateMPC_re07[2];
stateMPC_FLOAT stateMPC_beta07[2];
stateMPC_FLOAT stateMPC_betacc07[2];
stateMPC_FLOAT* stateMPC_dvaff07 = stateMPC_dv_aff + 14;
stateMPC_FLOAT* stateMPC_dvcc07 = stateMPC_dv_cc + 14;
stateMPC_FLOAT stateMPC_V07[8];
stateMPC_FLOAT stateMPC_Yd07[3];
stateMPC_FLOAT stateMPC_Ld07[3];
stateMPC_FLOAT stateMPC_yy07[2];
stateMPC_FLOAT stateMPC_bmy07[2];
int stateMPC_lbIdx07[4] = {0, 1, 2, 3};
stateMPC_FLOAT* stateMPC_llb07 = stateMPC_l + 56;
stateMPC_FLOAT* stateMPC_slb07 = stateMPC_s + 56;
stateMPC_FLOAT* stateMPC_llbbyslb07 = stateMPC_lbys + 56;
stateMPC_FLOAT stateMPC_rilb07[4];
stateMPC_FLOAT* stateMPC_dllbaff07 = stateMPC_dl_aff + 56;
stateMPC_FLOAT* stateMPC_dslbaff07 = stateMPC_ds_aff + 56;
stateMPC_FLOAT* stateMPC_dllbcc07 = stateMPC_dl_cc + 56;
stateMPC_FLOAT* stateMPC_dslbcc07 = stateMPC_ds_cc + 56;
stateMPC_FLOAT* stateMPC_ccrhsl07 = stateMPC_ccrhs + 56;
int stateMPC_ubIdx07[4] = {0, 1, 2, 3};
stateMPC_FLOAT* stateMPC_lub07 = stateMPC_l + 60;
stateMPC_FLOAT* stateMPC_sub07 = stateMPC_s + 60;
stateMPC_FLOAT* stateMPC_lubbysub07 = stateMPC_lbys + 60;
stateMPC_FLOAT stateMPC_riub07[4];
stateMPC_FLOAT* stateMPC_dlubaff07 = stateMPC_dl_aff + 60;
stateMPC_FLOAT* stateMPC_dsubaff07 = stateMPC_ds_aff + 60;
stateMPC_FLOAT* stateMPC_dlubcc07 = stateMPC_dl_cc + 60;
stateMPC_FLOAT* stateMPC_dsubcc07 = stateMPC_ds_cc + 60;
stateMPC_FLOAT* stateMPC_ccrhsub07 = stateMPC_ccrhs + 60;
stateMPC_FLOAT stateMPC_Phi07[4];
stateMPC_FLOAT stateMPC_W07[4];
stateMPC_FLOAT stateMPC_Ysd07[4];
stateMPC_FLOAT stateMPC_Lsd07[4];
stateMPC_FLOAT* stateMPC_z08 = stateMPC_z + 32;
stateMPC_FLOAT* stateMPC_dzaff08 = stateMPC_dz_aff + 32;
stateMPC_FLOAT* stateMPC_dzcc08 = stateMPC_dz_cc + 32;
stateMPC_FLOAT* stateMPC_rd08 = stateMPC_rd + 32;
stateMPC_FLOAT stateMPC_Lbyrd08[4];
stateMPC_FLOAT* stateMPC_grad_cost08 = stateMPC_grad_cost + 32;
stateMPC_FLOAT* stateMPC_grad_eq08 = stateMPC_grad_eq + 32;
stateMPC_FLOAT* stateMPC_grad_ineq08 = stateMPC_grad_ineq + 32;
stateMPC_FLOAT stateMPC_ctv08[4];
stateMPC_FLOAT* stateMPC_v08 = stateMPC_v + 16;
stateMPC_FLOAT stateMPC_re08[2];
stateMPC_FLOAT stateMPC_beta08[2];
stateMPC_FLOAT stateMPC_betacc08[2];
stateMPC_FLOAT* stateMPC_dvaff08 = stateMPC_dv_aff + 16;
stateMPC_FLOAT* stateMPC_dvcc08 = stateMPC_dv_cc + 16;
stateMPC_FLOAT stateMPC_V08[8];
stateMPC_FLOAT stateMPC_Yd08[3];
stateMPC_FLOAT stateMPC_Ld08[3];
stateMPC_FLOAT stateMPC_yy08[2];
stateMPC_FLOAT stateMPC_bmy08[2];
int stateMPC_lbIdx08[4] = {0, 1, 2, 3};
stateMPC_FLOAT* stateMPC_llb08 = stateMPC_l + 64;
stateMPC_FLOAT* stateMPC_slb08 = stateMPC_s + 64;
stateMPC_FLOAT* stateMPC_llbbyslb08 = stateMPC_lbys + 64;
stateMPC_FLOAT stateMPC_rilb08[4];
stateMPC_FLOAT* stateMPC_dllbaff08 = stateMPC_dl_aff + 64;
stateMPC_FLOAT* stateMPC_dslbaff08 = stateMPC_ds_aff + 64;
stateMPC_FLOAT* stateMPC_dllbcc08 = stateMPC_dl_cc + 64;
stateMPC_FLOAT* stateMPC_dslbcc08 = stateMPC_ds_cc + 64;
stateMPC_FLOAT* stateMPC_ccrhsl08 = stateMPC_ccrhs + 64;
int stateMPC_ubIdx08[4] = {0, 1, 2, 3};
stateMPC_FLOAT* stateMPC_lub08 = stateMPC_l + 68;
stateMPC_FLOAT* stateMPC_sub08 = stateMPC_s + 68;
stateMPC_FLOAT* stateMPC_lubbysub08 = stateMPC_lbys + 68;
stateMPC_FLOAT stateMPC_riub08[4];
stateMPC_FLOAT* stateMPC_dlubaff08 = stateMPC_dl_aff + 68;
stateMPC_FLOAT* stateMPC_dsubaff08 = stateMPC_ds_aff + 68;
stateMPC_FLOAT* stateMPC_dlubcc08 = stateMPC_dl_cc + 68;
stateMPC_FLOAT* stateMPC_dsubcc08 = stateMPC_ds_cc + 68;
stateMPC_FLOAT* stateMPC_ccrhsub08 = stateMPC_ccrhs + 68;
stateMPC_FLOAT stateMPC_Phi08[4];
stateMPC_FLOAT stateMPC_W08[4];
stateMPC_FLOAT stateMPC_Ysd08[4];
stateMPC_FLOAT stateMPC_Lsd08[4];
stateMPC_FLOAT* stateMPC_z09 = stateMPC_z + 36;
stateMPC_FLOAT* stateMPC_dzaff09 = stateMPC_dz_aff + 36;
stateMPC_FLOAT* stateMPC_dzcc09 = stateMPC_dz_cc + 36;
stateMPC_FLOAT* stateMPC_rd09 = stateMPC_rd + 36;
stateMPC_FLOAT stateMPC_Lbyrd09[4];
stateMPC_FLOAT* stateMPC_grad_cost09 = stateMPC_grad_cost + 36;
stateMPC_FLOAT* stateMPC_grad_eq09 = stateMPC_grad_eq + 36;
stateMPC_FLOAT* stateMPC_grad_ineq09 = stateMPC_grad_ineq + 36;
stateMPC_FLOAT stateMPC_ctv09[4];
stateMPC_FLOAT* stateMPC_v09 = stateMPC_v + 18;
stateMPC_FLOAT stateMPC_re09[2];
stateMPC_FLOAT stateMPC_beta09[2];
stateMPC_FLOAT stateMPC_betacc09[2];
stateMPC_FLOAT* stateMPC_dvaff09 = stateMPC_dv_aff + 18;
stateMPC_FLOAT* stateMPC_dvcc09 = stateMPC_dv_cc + 18;
stateMPC_FLOAT stateMPC_V09[8];
stateMPC_FLOAT stateMPC_Yd09[3];
stateMPC_FLOAT stateMPC_Ld09[3];
stateMPC_FLOAT stateMPC_yy09[2];
stateMPC_FLOAT stateMPC_bmy09[2];
int stateMPC_lbIdx09[4] = {0, 1, 2, 3};
stateMPC_FLOAT* stateMPC_llb09 = stateMPC_l + 72;
stateMPC_FLOAT* stateMPC_slb09 = stateMPC_s + 72;
stateMPC_FLOAT* stateMPC_llbbyslb09 = stateMPC_lbys + 72;
stateMPC_FLOAT stateMPC_rilb09[4];
stateMPC_FLOAT* stateMPC_dllbaff09 = stateMPC_dl_aff + 72;
stateMPC_FLOAT* stateMPC_dslbaff09 = stateMPC_ds_aff + 72;
stateMPC_FLOAT* stateMPC_dllbcc09 = stateMPC_dl_cc + 72;
stateMPC_FLOAT* stateMPC_dslbcc09 = stateMPC_ds_cc + 72;
stateMPC_FLOAT* stateMPC_ccrhsl09 = stateMPC_ccrhs + 72;
int stateMPC_ubIdx09[4] = {0, 1, 2, 3};
stateMPC_FLOAT* stateMPC_lub09 = stateMPC_l + 76;
stateMPC_FLOAT* stateMPC_sub09 = stateMPC_s + 76;
stateMPC_FLOAT* stateMPC_lubbysub09 = stateMPC_lbys + 76;
stateMPC_FLOAT stateMPC_riub09[4];
stateMPC_FLOAT* stateMPC_dlubaff09 = stateMPC_dl_aff + 76;
stateMPC_FLOAT* stateMPC_dsubaff09 = stateMPC_ds_aff + 76;
stateMPC_FLOAT* stateMPC_dlubcc09 = stateMPC_dl_cc + 76;
stateMPC_FLOAT* stateMPC_dsubcc09 = stateMPC_ds_cc + 76;
stateMPC_FLOAT* stateMPC_ccrhsub09 = stateMPC_ccrhs + 76;
stateMPC_FLOAT stateMPC_Phi09[4];
stateMPC_FLOAT stateMPC_W09[4];
stateMPC_FLOAT stateMPC_Ysd09[4];
stateMPC_FLOAT stateMPC_Lsd09[4];
stateMPC_FLOAT* stateMPC_z10 = stateMPC_z + 40;
stateMPC_FLOAT* stateMPC_dzaff10 = stateMPC_dz_aff + 40;
stateMPC_FLOAT* stateMPC_dzcc10 = stateMPC_dz_cc + 40;
stateMPC_FLOAT* stateMPC_rd10 = stateMPC_rd + 40;
stateMPC_FLOAT stateMPC_Lbyrd10[4];
stateMPC_FLOAT* stateMPC_grad_cost10 = stateMPC_grad_cost + 40;
stateMPC_FLOAT* stateMPC_grad_eq10 = stateMPC_grad_eq + 40;
stateMPC_FLOAT* stateMPC_grad_ineq10 = stateMPC_grad_ineq + 40;
stateMPC_FLOAT stateMPC_ctv10[4];
stateMPC_FLOAT* stateMPC_v10 = stateMPC_v + 20;
stateMPC_FLOAT stateMPC_re10[2];
stateMPC_FLOAT stateMPC_beta10[2];
stateMPC_FLOAT stateMPC_betacc10[2];
stateMPC_FLOAT* stateMPC_dvaff10 = stateMPC_dv_aff + 20;
stateMPC_FLOAT* stateMPC_dvcc10 = stateMPC_dv_cc + 20;
stateMPC_FLOAT stateMPC_V10[8];
stateMPC_FLOAT stateMPC_Yd10[3];
stateMPC_FLOAT stateMPC_Ld10[3];
stateMPC_FLOAT stateMPC_yy10[2];
stateMPC_FLOAT stateMPC_bmy10[2];
int stateMPC_lbIdx10[4] = {0, 1, 2, 3};
stateMPC_FLOAT* stateMPC_llb10 = stateMPC_l + 80;
stateMPC_FLOAT* stateMPC_slb10 = stateMPC_s + 80;
stateMPC_FLOAT* stateMPC_llbbyslb10 = stateMPC_lbys + 80;
stateMPC_FLOAT stateMPC_rilb10[4];
stateMPC_FLOAT* stateMPC_dllbaff10 = stateMPC_dl_aff + 80;
stateMPC_FLOAT* stateMPC_dslbaff10 = stateMPC_ds_aff + 80;
stateMPC_FLOAT* stateMPC_dllbcc10 = stateMPC_dl_cc + 80;
stateMPC_FLOAT* stateMPC_dslbcc10 = stateMPC_ds_cc + 80;
stateMPC_FLOAT* stateMPC_ccrhsl10 = stateMPC_ccrhs + 80;
int stateMPC_ubIdx10[4] = {0, 1, 2, 3};
stateMPC_FLOAT* stateMPC_lub10 = stateMPC_l + 84;
stateMPC_FLOAT* stateMPC_sub10 = stateMPC_s + 84;
stateMPC_FLOAT* stateMPC_lubbysub10 = stateMPC_lbys + 84;
stateMPC_FLOAT stateMPC_riub10[4];
stateMPC_FLOAT* stateMPC_dlubaff10 = stateMPC_dl_aff + 84;
stateMPC_FLOAT* stateMPC_dsubaff10 = stateMPC_ds_aff + 84;
stateMPC_FLOAT* stateMPC_dlubcc10 = stateMPC_dl_cc + 84;
stateMPC_FLOAT* stateMPC_dsubcc10 = stateMPC_ds_cc + 84;
stateMPC_FLOAT* stateMPC_ccrhsub10 = stateMPC_ccrhs + 84;
stateMPC_FLOAT stateMPC_Phi10[4];
stateMPC_FLOAT stateMPC_W10[4];
stateMPC_FLOAT stateMPC_Ysd10[4];
stateMPC_FLOAT stateMPC_Lsd10[4];
stateMPC_FLOAT* stateMPC_z11 = stateMPC_z + 44;
stateMPC_FLOAT* stateMPC_dzaff11 = stateMPC_dz_aff + 44;
stateMPC_FLOAT* stateMPC_dzcc11 = stateMPC_dz_cc + 44;
stateMPC_FLOAT* stateMPC_rd11 = stateMPC_rd + 44;
stateMPC_FLOAT stateMPC_Lbyrd11[4];
stateMPC_FLOAT* stateMPC_grad_cost11 = stateMPC_grad_cost + 44;
stateMPC_FLOAT* stateMPC_grad_eq11 = stateMPC_grad_eq + 44;
stateMPC_FLOAT* stateMPC_grad_ineq11 = stateMPC_grad_ineq + 44;
stateMPC_FLOAT stateMPC_ctv11[4];
stateMPC_FLOAT* stateMPC_v11 = stateMPC_v + 22;
stateMPC_FLOAT stateMPC_re11[2];
stateMPC_FLOAT stateMPC_beta11[2];
stateMPC_FLOAT stateMPC_betacc11[2];
stateMPC_FLOAT* stateMPC_dvaff11 = stateMPC_dv_aff + 22;
stateMPC_FLOAT* stateMPC_dvcc11 = stateMPC_dv_cc + 22;
stateMPC_FLOAT stateMPC_V11[8];
stateMPC_FLOAT stateMPC_Yd11[3];
stateMPC_FLOAT stateMPC_Ld11[3];
stateMPC_FLOAT stateMPC_yy11[2];
stateMPC_FLOAT stateMPC_bmy11[2];
int stateMPC_lbIdx11[4] = {0, 1, 2, 3};
stateMPC_FLOAT* stateMPC_llb11 = stateMPC_l + 88;
stateMPC_FLOAT* stateMPC_slb11 = stateMPC_s + 88;
stateMPC_FLOAT* stateMPC_llbbyslb11 = stateMPC_lbys + 88;
stateMPC_FLOAT stateMPC_rilb11[4];
stateMPC_FLOAT* stateMPC_dllbaff11 = stateMPC_dl_aff + 88;
stateMPC_FLOAT* stateMPC_dslbaff11 = stateMPC_ds_aff + 88;
stateMPC_FLOAT* stateMPC_dllbcc11 = stateMPC_dl_cc + 88;
stateMPC_FLOAT* stateMPC_dslbcc11 = stateMPC_ds_cc + 88;
stateMPC_FLOAT* stateMPC_ccrhsl11 = stateMPC_ccrhs + 88;
int stateMPC_ubIdx11[4] = {0, 1, 2, 3};
stateMPC_FLOAT* stateMPC_lub11 = stateMPC_l + 92;
stateMPC_FLOAT* stateMPC_sub11 = stateMPC_s + 92;
stateMPC_FLOAT* stateMPC_lubbysub11 = stateMPC_lbys + 92;
stateMPC_FLOAT stateMPC_riub11[4];
stateMPC_FLOAT* stateMPC_dlubaff11 = stateMPC_dl_aff + 92;
stateMPC_FLOAT* stateMPC_dsubaff11 = stateMPC_ds_aff + 92;
stateMPC_FLOAT* stateMPC_dlubcc11 = stateMPC_dl_cc + 92;
stateMPC_FLOAT* stateMPC_dsubcc11 = stateMPC_ds_cc + 92;
stateMPC_FLOAT* stateMPC_ccrhsub11 = stateMPC_ccrhs + 92;
stateMPC_FLOAT stateMPC_Phi11[4];
stateMPC_FLOAT stateMPC_W11[4];
stateMPC_FLOAT stateMPC_Ysd11[4];
stateMPC_FLOAT stateMPC_Lsd11[4];
stateMPC_FLOAT* stateMPC_z12 = stateMPC_z + 48;
stateMPC_FLOAT* stateMPC_dzaff12 = stateMPC_dz_aff + 48;
stateMPC_FLOAT* stateMPC_dzcc12 = stateMPC_dz_cc + 48;
stateMPC_FLOAT* stateMPC_rd12 = stateMPC_rd + 48;
stateMPC_FLOAT stateMPC_Lbyrd12[4];
stateMPC_FLOAT* stateMPC_grad_cost12 = stateMPC_grad_cost + 48;
stateMPC_FLOAT* stateMPC_grad_eq12 = stateMPC_grad_eq + 48;
stateMPC_FLOAT* stateMPC_grad_ineq12 = stateMPC_grad_ineq + 48;
stateMPC_FLOAT stateMPC_ctv12[4];
stateMPC_FLOAT* stateMPC_v12 = stateMPC_v + 24;
stateMPC_FLOAT stateMPC_re12[2];
stateMPC_FLOAT stateMPC_beta12[2];
stateMPC_FLOAT stateMPC_betacc12[2];
stateMPC_FLOAT* stateMPC_dvaff12 = stateMPC_dv_aff + 24;
stateMPC_FLOAT* stateMPC_dvcc12 = stateMPC_dv_cc + 24;
stateMPC_FLOAT stateMPC_V12[8];
stateMPC_FLOAT stateMPC_Yd12[3];
stateMPC_FLOAT stateMPC_Ld12[3];
stateMPC_FLOAT stateMPC_yy12[2];
stateMPC_FLOAT stateMPC_bmy12[2];
int stateMPC_lbIdx12[4] = {0, 1, 2, 3};
stateMPC_FLOAT* stateMPC_llb12 = stateMPC_l + 96;
stateMPC_FLOAT* stateMPC_slb12 = stateMPC_s + 96;
stateMPC_FLOAT* stateMPC_llbbyslb12 = stateMPC_lbys + 96;
stateMPC_FLOAT stateMPC_rilb12[4];
stateMPC_FLOAT* stateMPC_dllbaff12 = stateMPC_dl_aff + 96;
stateMPC_FLOAT* stateMPC_dslbaff12 = stateMPC_ds_aff + 96;
stateMPC_FLOAT* stateMPC_dllbcc12 = stateMPC_dl_cc + 96;
stateMPC_FLOAT* stateMPC_dslbcc12 = stateMPC_ds_cc + 96;
stateMPC_FLOAT* stateMPC_ccrhsl12 = stateMPC_ccrhs + 96;
int stateMPC_ubIdx12[4] = {0, 1, 2, 3};
stateMPC_FLOAT* stateMPC_lub12 = stateMPC_l + 100;
stateMPC_FLOAT* stateMPC_sub12 = stateMPC_s + 100;
stateMPC_FLOAT* stateMPC_lubbysub12 = stateMPC_lbys + 100;
stateMPC_FLOAT stateMPC_riub12[4];
stateMPC_FLOAT* stateMPC_dlubaff12 = stateMPC_dl_aff + 100;
stateMPC_FLOAT* stateMPC_dsubaff12 = stateMPC_ds_aff + 100;
stateMPC_FLOAT* stateMPC_dlubcc12 = stateMPC_dl_cc + 100;
stateMPC_FLOAT* stateMPC_dsubcc12 = stateMPC_ds_cc + 100;
stateMPC_FLOAT* stateMPC_ccrhsub12 = stateMPC_ccrhs + 100;
stateMPC_FLOAT stateMPC_Phi12[4];
stateMPC_FLOAT stateMPC_W12[4];
stateMPC_FLOAT stateMPC_Ysd12[4];
stateMPC_FLOAT stateMPC_Lsd12[4];
stateMPC_FLOAT* stateMPC_z13 = stateMPC_z + 52;
stateMPC_FLOAT* stateMPC_dzaff13 = stateMPC_dz_aff + 52;
stateMPC_FLOAT* stateMPC_dzcc13 = stateMPC_dz_cc + 52;
stateMPC_FLOAT* stateMPC_rd13 = stateMPC_rd + 52;
stateMPC_FLOAT stateMPC_Lbyrd13[4];
stateMPC_FLOAT* stateMPC_grad_cost13 = stateMPC_grad_cost + 52;
stateMPC_FLOAT* stateMPC_grad_eq13 = stateMPC_grad_eq + 52;
stateMPC_FLOAT* stateMPC_grad_ineq13 = stateMPC_grad_ineq + 52;
stateMPC_FLOAT stateMPC_ctv13[4];
stateMPC_FLOAT* stateMPC_v13 = stateMPC_v + 26;
stateMPC_FLOAT stateMPC_re13[2];
stateMPC_FLOAT stateMPC_beta13[2];
stateMPC_FLOAT stateMPC_betacc13[2];
stateMPC_FLOAT* stateMPC_dvaff13 = stateMPC_dv_aff + 26;
stateMPC_FLOAT* stateMPC_dvcc13 = stateMPC_dv_cc + 26;
stateMPC_FLOAT stateMPC_V13[8];
stateMPC_FLOAT stateMPC_Yd13[3];
stateMPC_FLOAT stateMPC_Ld13[3];
stateMPC_FLOAT stateMPC_yy13[2];
stateMPC_FLOAT stateMPC_bmy13[2];
int stateMPC_lbIdx13[4] = {0, 1, 2, 3};
stateMPC_FLOAT* stateMPC_llb13 = stateMPC_l + 104;
stateMPC_FLOAT* stateMPC_slb13 = stateMPC_s + 104;
stateMPC_FLOAT* stateMPC_llbbyslb13 = stateMPC_lbys + 104;
stateMPC_FLOAT stateMPC_rilb13[4];
stateMPC_FLOAT* stateMPC_dllbaff13 = stateMPC_dl_aff + 104;
stateMPC_FLOAT* stateMPC_dslbaff13 = stateMPC_ds_aff + 104;
stateMPC_FLOAT* stateMPC_dllbcc13 = stateMPC_dl_cc + 104;
stateMPC_FLOAT* stateMPC_dslbcc13 = stateMPC_ds_cc + 104;
stateMPC_FLOAT* stateMPC_ccrhsl13 = stateMPC_ccrhs + 104;
int stateMPC_ubIdx13[4] = {0, 1, 2, 3};
stateMPC_FLOAT* stateMPC_lub13 = stateMPC_l + 108;
stateMPC_FLOAT* stateMPC_sub13 = stateMPC_s + 108;
stateMPC_FLOAT* stateMPC_lubbysub13 = stateMPC_lbys + 108;
stateMPC_FLOAT stateMPC_riub13[4];
stateMPC_FLOAT* stateMPC_dlubaff13 = stateMPC_dl_aff + 108;
stateMPC_FLOAT* stateMPC_dsubaff13 = stateMPC_ds_aff + 108;
stateMPC_FLOAT* stateMPC_dlubcc13 = stateMPC_dl_cc + 108;
stateMPC_FLOAT* stateMPC_dsubcc13 = stateMPC_ds_cc + 108;
stateMPC_FLOAT* stateMPC_ccrhsub13 = stateMPC_ccrhs + 108;
stateMPC_FLOAT stateMPC_Phi13[4];
stateMPC_FLOAT stateMPC_W13[4];
stateMPC_FLOAT stateMPC_Ysd13[4];
stateMPC_FLOAT stateMPC_Lsd13[4];
stateMPC_FLOAT* stateMPC_z14 = stateMPC_z + 56;
stateMPC_FLOAT* stateMPC_dzaff14 = stateMPC_dz_aff + 56;
stateMPC_FLOAT* stateMPC_dzcc14 = stateMPC_dz_cc + 56;
stateMPC_FLOAT* stateMPC_rd14 = stateMPC_rd + 56;
stateMPC_FLOAT stateMPC_Lbyrd14[2];
stateMPC_FLOAT* stateMPC_grad_cost14 = stateMPC_grad_cost + 56;
stateMPC_FLOAT* stateMPC_grad_eq14 = stateMPC_grad_eq + 56;
stateMPC_FLOAT* stateMPC_grad_ineq14 = stateMPC_grad_ineq + 56;
stateMPC_FLOAT stateMPC_ctv14[2];
stateMPC_FLOAT* stateMPC_v14 = stateMPC_v + 28;
stateMPC_FLOAT stateMPC_re14[2];
stateMPC_FLOAT stateMPC_beta14[2];
stateMPC_FLOAT stateMPC_betacc14[2];
stateMPC_FLOAT* stateMPC_dvaff14 = stateMPC_dv_aff + 28;
stateMPC_FLOAT* stateMPC_dvcc14 = stateMPC_dv_cc + 28;
stateMPC_FLOAT stateMPC_V14[4];
stateMPC_FLOAT stateMPC_Yd14[3];
stateMPC_FLOAT stateMPC_Ld14[3];
stateMPC_FLOAT stateMPC_yy14[2];
stateMPC_FLOAT stateMPC_bmy14[2];
int stateMPC_lbIdx14[2] = {0, 1};
stateMPC_FLOAT* stateMPC_llb14 = stateMPC_l + 112;
stateMPC_FLOAT* stateMPC_slb14 = stateMPC_s + 112;
stateMPC_FLOAT* stateMPC_llbbyslb14 = stateMPC_lbys + 112;
stateMPC_FLOAT stateMPC_rilb14[2];
stateMPC_FLOAT* stateMPC_dllbaff14 = stateMPC_dl_aff + 112;
stateMPC_FLOAT* stateMPC_dslbaff14 = stateMPC_ds_aff + 112;
stateMPC_FLOAT* stateMPC_dllbcc14 = stateMPC_dl_cc + 112;
stateMPC_FLOAT* stateMPC_dslbcc14 = stateMPC_ds_cc + 112;
stateMPC_FLOAT* stateMPC_ccrhsl14 = stateMPC_ccrhs + 112;
int stateMPC_ubIdx14[2] = {0, 1};
stateMPC_FLOAT* stateMPC_lub14 = stateMPC_l + 114;
stateMPC_FLOAT* stateMPC_sub14 = stateMPC_s + 114;
stateMPC_FLOAT* stateMPC_lubbysub14 = stateMPC_lbys + 114;
stateMPC_FLOAT stateMPC_riub14[2];
stateMPC_FLOAT* stateMPC_dlubaff14 = stateMPC_dl_aff + 114;
stateMPC_FLOAT* stateMPC_dsubaff14 = stateMPC_ds_aff + 114;
stateMPC_FLOAT* stateMPC_dlubcc14 = stateMPC_dl_cc + 114;
stateMPC_FLOAT* stateMPC_dsubcc14 = stateMPC_ds_cc + 114;
stateMPC_FLOAT* stateMPC_ccrhsub14 = stateMPC_ccrhs + 114;
stateMPC_FLOAT stateMPC_Phi14[2];
stateMPC_FLOAT stateMPC_D14[2] = {-1.0000000000000000E+000, 
-1.0000000000000000E+000};
stateMPC_FLOAT stateMPC_W14[2];
stateMPC_FLOAT stateMPC_Ysd14[4];
stateMPC_FLOAT stateMPC_Lsd14[4];
stateMPC_FLOAT musigma;
stateMPC_FLOAT sigma_3rdroot;
stateMPC_FLOAT stateMPC_Diag1_0[4];
stateMPC_FLOAT stateMPC_Diag2_0[4];
stateMPC_FLOAT stateMPC_L_0[6];




/* SOLVER CODE --------------------------------------------------------- */
int stateMPC_solve(stateMPC_params* params, stateMPC_output* output, stateMPC_info* info)
{	
int exitcode;

#if stateMPC_SET_TIMING == 1
	stateMPC_timer solvertimer;
	stateMPC_tic(&solvertimer);
#endif
/* FUNCTION CALLS INTO LA LIBRARY -------------------------------------- */
info->it = 0;
stateMPC_LA_INITIALIZEVECTOR_58(stateMPC_z, 0);
stateMPC_LA_INITIALIZEVECTOR_30(stateMPC_v, 1);
stateMPC_LA_INITIALIZEVECTOR_116(stateMPC_l, 1);
stateMPC_LA_INITIALIZEVECTOR_116(stateMPC_s, 1);
info->mu = 0;
stateMPC_LA_DOTACC_116(stateMPC_l, stateMPC_s, &info->mu);
info->mu /= 116;
while( 1 ){
info->pobj = 0;
stateMPC_LA_DIAG_QUADFCN_4(params->H1, params->f1, stateMPC_z00, stateMPC_grad_cost00, &info->pobj);
stateMPC_LA_DIAG_QUADFCN_4(params->H2, params->f2, stateMPC_z01, stateMPC_grad_cost01, &info->pobj);
stateMPC_LA_DIAG_QUADFCN_4(params->H3, params->f3, stateMPC_z02, stateMPC_grad_cost02, &info->pobj);
stateMPC_LA_DIAG_QUADFCN_4(params->H4, params->f4, stateMPC_z03, stateMPC_grad_cost03, &info->pobj);
stateMPC_LA_DIAG_QUADFCN_4(params->H5, params->f5, stateMPC_z04, stateMPC_grad_cost04, &info->pobj);
stateMPC_LA_DIAG_QUADFCN_4(params->H6, params->f6, stateMPC_z05, stateMPC_grad_cost05, &info->pobj);
stateMPC_LA_DIAG_QUADFCN_4(params->H7, params->f7, stateMPC_z06, stateMPC_grad_cost06, &info->pobj);
stateMPC_LA_DIAG_QUADFCN_4(params->H8, params->f8, stateMPC_z07, stateMPC_grad_cost07, &info->pobj);
stateMPC_LA_DIAG_QUADFCN_4(params->H9, params->f9, stateMPC_z08, stateMPC_grad_cost08, &info->pobj);
stateMPC_LA_DIAG_QUADFCN_4(params->H10, params->f10, stateMPC_z09, stateMPC_grad_cost09, &info->pobj);
stateMPC_LA_DIAG_QUADFCN_4(params->H11, params->f11, stateMPC_z10, stateMPC_grad_cost10, &info->pobj);
stateMPC_LA_DIAG_QUADFCN_4(params->H12, params->f12, stateMPC_z11, stateMPC_grad_cost11, &info->pobj);
stateMPC_LA_DIAG_QUADFCN_4(params->H13, params->f13, stateMPC_z12, stateMPC_grad_cost12, &info->pobj);
stateMPC_LA_DIAG_QUADFCN_4(params->H14, params->f14, stateMPC_z13, stateMPC_grad_cost13, &info->pobj);
stateMPC_LA_DIAG_QUADFCN_2(params->H15, params->f15, stateMPC_z14, stateMPC_grad_cost14, &info->pobj);
info->res_eq = 0;
info->dgap = 0;
stateMPC_LA_DIAGZERO_MVMSUB6_2(stateMPC_D00, stateMPC_z00, params->e1, stateMPC_v00, stateMPC_re00, &info->dgap, &info->res_eq);
stateMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_4_4(params->C1, stateMPC_z00, stateMPC_D01, stateMPC_z01, params->e2, stateMPC_v01, stateMPC_re01, &info->dgap, &info->res_eq);
stateMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_4_4(params->C2, stateMPC_z01, stateMPC_D01, stateMPC_z02, params->e3, stateMPC_v02, stateMPC_re02, &info->dgap, &info->res_eq);
stateMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_4_4(params->C3, stateMPC_z02, stateMPC_D01, stateMPC_z03, params->e4, stateMPC_v03, stateMPC_re03, &info->dgap, &info->res_eq);
stateMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_4_4(params->C4, stateMPC_z03, stateMPC_D01, stateMPC_z04, params->e5, stateMPC_v04, stateMPC_re04, &info->dgap, &info->res_eq);
stateMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_4_4(params->C5, stateMPC_z04, stateMPC_D01, stateMPC_z05, params->e6, stateMPC_v05, stateMPC_re05, &info->dgap, &info->res_eq);
stateMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_4_4(params->C6, stateMPC_z05, stateMPC_D01, stateMPC_z06, params->e7, stateMPC_v06, stateMPC_re06, &info->dgap, &info->res_eq);
stateMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_4_4(params->C7, stateMPC_z06, stateMPC_D01, stateMPC_z07, params->e8, stateMPC_v07, stateMPC_re07, &info->dgap, &info->res_eq);
stateMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_4_4(params->C8, stateMPC_z07, stateMPC_D01, stateMPC_z08, params->e9, stateMPC_v08, stateMPC_re08, &info->dgap, &info->res_eq);
stateMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_4_4(params->C9, stateMPC_z08, stateMPC_D01, stateMPC_z09, params->e10, stateMPC_v09, stateMPC_re09, &info->dgap, &info->res_eq);
stateMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_4_4(params->C10, stateMPC_z09, stateMPC_D01, stateMPC_z10, params->e11, stateMPC_v10, stateMPC_re10, &info->dgap, &info->res_eq);
stateMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_4_4(params->C11, stateMPC_z10, stateMPC_D01, stateMPC_z11, params->e12, stateMPC_v11, stateMPC_re11, &info->dgap, &info->res_eq);
stateMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_4_4(params->C12, stateMPC_z11, stateMPC_D01, stateMPC_z12, params->e13, stateMPC_v12, stateMPC_re12, &info->dgap, &info->res_eq);
stateMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_4_4(params->C13, stateMPC_z12, stateMPC_D01, stateMPC_z13, params->e14, stateMPC_v13, stateMPC_re13, &info->dgap, &info->res_eq);
stateMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_4_2(params->C14, stateMPC_z13, stateMPC_D14, stateMPC_z14, params->e15, stateMPC_v14, stateMPC_re14, &info->dgap, &info->res_eq);
stateMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C1, stateMPC_v01, stateMPC_D00, stateMPC_v00, stateMPC_grad_eq00);
stateMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C2, stateMPC_v02, stateMPC_D01, stateMPC_v01, stateMPC_grad_eq01);
stateMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C3, stateMPC_v03, stateMPC_D01, stateMPC_v02, stateMPC_grad_eq02);
stateMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C4, stateMPC_v04, stateMPC_D01, stateMPC_v03, stateMPC_grad_eq03);
stateMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C5, stateMPC_v05, stateMPC_D01, stateMPC_v04, stateMPC_grad_eq04);
stateMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C6, stateMPC_v06, stateMPC_D01, stateMPC_v05, stateMPC_grad_eq05);
stateMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C7, stateMPC_v07, stateMPC_D01, stateMPC_v06, stateMPC_grad_eq06);
stateMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C8, stateMPC_v08, stateMPC_D01, stateMPC_v07, stateMPC_grad_eq07);
stateMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C9, stateMPC_v09, stateMPC_D01, stateMPC_v08, stateMPC_grad_eq08);
stateMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C10, stateMPC_v10, stateMPC_D01, stateMPC_v09, stateMPC_grad_eq09);
stateMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C11, stateMPC_v11, stateMPC_D01, stateMPC_v10, stateMPC_grad_eq10);
stateMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C12, stateMPC_v12, stateMPC_D01, stateMPC_v11, stateMPC_grad_eq11);
stateMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C13, stateMPC_v13, stateMPC_D01, stateMPC_v12, stateMPC_grad_eq12);
stateMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C14, stateMPC_v14, stateMPC_D01, stateMPC_v13, stateMPC_grad_eq13);
stateMPC_LA_DIAGZERO_MTVM_2_2(stateMPC_D14, stateMPC_v14, stateMPC_grad_eq14);
info->res_ineq = 0;
stateMPC_LA_VSUBADD3_4(params->lb1, stateMPC_z00, stateMPC_lbIdx00, stateMPC_llb00, stateMPC_slb00, stateMPC_rilb00, &info->dgap, &info->res_ineq);
stateMPC_LA_VSUBADD2_4(stateMPC_z00, stateMPC_ubIdx00, params->ub1, stateMPC_lub00, stateMPC_sub00, stateMPC_riub00, &info->dgap, &info->res_ineq);
stateMPC_LA_VSUBADD3_4(params->lb2, stateMPC_z01, stateMPC_lbIdx01, stateMPC_llb01, stateMPC_slb01, stateMPC_rilb01, &info->dgap, &info->res_ineq);
stateMPC_LA_VSUBADD2_4(stateMPC_z01, stateMPC_ubIdx01, params->ub2, stateMPC_lub01, stateMPC_sub01, stateMPC_riub01, &info->dgap, &info->res_ineq);
stateMPC_LA_VSUBADD3_4(params->lb3, stateMPC_z02, stateMPC_lbIdx02, stateMPC_llb02, stateMPC_slb02, stateMPC_rilb02, &info->dgap, &info->res_ineq);
stateMPC_LA_VSUBADD2_4(stateMPC_z02, stateMPC_ubIdx02, params->ub3, stateMPC_lub02, stateMPC_sub02, stateMPC_riub02, &info->dgap, &info->res_ineq);
stateMPC_LA_VSUBADD3_4(params->lb4, stateMPC_z03, stateMPC_lbIdx03, stateMPC_llb03, stateMPC_slb03, stateMPC_rilb03, &info->dgap, &info->res_ineq);
stateMPC_LA_VSUBADD2_4(stateMPC_z03, stateMPC_ubIdx03, params->ub4, stateMPC_lub03, stateMPC_sub03, stateMPC_riub03, &info->dgap, &info->res_ineq);
stateMPC_LA_VSUBADD3_4(params->lb5, stateMPC_z04, stateMPC_lbIdx04, stateMPC_llb04, stateMPC_slb04, stateMPC_rilb04, &info->dgap, &info->res_ineq);
stateMPC_LA_VSUBADD2_4(stateMPC_z04, stateMPC_ubIdx04, params->ub5, stateMPC_lub04, stateMPC_sub04, stateMPC_riub04, &info->dgap, &info->res_ineq);
stateMPC_LA_VSUBADD3_4(params->lb6, stateMPC_z05, stateMPC_lbIdx05, stateMPC_llb05, stateMPC_slb05, stateMPC_rilb05, &info->dgap, &info->res_ineq);
stateMPC_LA_VSUBADD2_4(stateMPC_z05, stateMPC_ubIdx05, params->ub6, stateMPC_lub05, stateMPC_sub05, stateMPC_riub05, &info->dgap, &info->res_ineq);
stateMPC_LA_VSUBADD3_4(params->lb7, stateMPC_z06, stateMPC_lbIdx06, stateMPC_llb06, stateMPC_slb06, stateMPC_rilb06, &info->dgap, &info->res_ineq);
stateMPC_LA_VSUBADD2_4(stateMPC_z06, stateMPC_ubIdx06, params->ub7, stateMPC_lub06, stateMPC_sub06, stateMPC_riub06, &info->dgap, &info->res_ineq);
stateMPC_LA_VSUBADD3_4(params->lb8, stateMPC_z07, stateMPC_lbIdx07, stateMPC_llb07, stateMPC_slb07, stateMPC_rilb07, &info->dgap, &info->res_ineq);
stateMPC_LA_VSUBADD2_4(stateMPC_z07, stateMPC_ubIdx07, params->ub8, stateMPC_lub07, stateMPC_sub07, stateMPC_riub07, &info->dgap, &info->res_ineq);
stateMPC_LA_VSUBADD3_4(params->lb9, stateMPC_z08, stateMPC_lbIdx08, stateMPC_llb08, stateMPC_slb08, stateMPC_rilb08, &info->dgap, &info->res_ineq);
stateMPC_LA_VSUBADD2_4(stateMPC_z08, stateMPC_ubIdx08, params->ub9, stateMPC_lub08, stateMPC_sub08, stateMPC_riub08, &info->dgap, &info->res_ineq);
stateMPC_LA_VSUBADD3_4(params->lb10, stateMPC_z09, stateMPC_lbIdx09, stateMPC_llb09, stateMPC_slb09, stateMPC_rilb09, &info->dgap, &info->res_ineq);
stateMPC_LA_VSUBADD2_4(stateMPC_z09, stateMPC_ubIdx09, params->ub10, stateMPC_lub09, stateMPC_sub09, stateMPC_riub09, &info->dgap, &info->res_ineq);
stateMPC_LA_VSUBADD3_4(params->lb11, stateMPC_z10, stateMPC_lbIdx10, stateMPC_llb10, stateMPC_slb10, stateMPC_rilb10, &info->dgap, &info->res_ineq);
stateMPC_LA_VSUBADD2_4(stateMPC_z10, stateMPC_ubIdx10, params->ub11, stateMPC_lub10, stateMPC_sub10, stateMPC_riub10, &info->dgap, &info->res_ineq);
stateMPC_LA_VSUBADD3_4(params->lb12, stateMPC_z11, stateMPC_lbIdx11, stateMPC_llb11, stateMPC_slb11, stateMPC_rilb11, &info->dgap, &info->res_ineq);
stateMPC_LA_VSUBADD2_4(stateMPC_z11, stateMPC_ubIdx11, params->ub12, stateMPC_lub11, stateMPC_sub11, stateMPC_riub11, &info->dgap, &info->res_ineq);
stateMPC_LA_VSUBADD3_4(params->lb13, stateMPC_z12, stateMPC_lbIdx12, stateMPC_llb12, stateMPC_slb12, stateMPC_rilb12, &info->dgap, &info->res_ineq);
stateMPC_LA_VSUBADD2_4(stateMPC_z12, stateMPC_ubIdx12, params->ub13, stateMPC_lub12, stateMPC_sub12, stateMPC_riub12, &info->dgap, &info->res_ineq);
stateMPC_LA_VSUBADD3_4(params->lb14, stateMPC_z13, stateMPC_lbIdx13, stateMPC_llb13, stateMPC_slb13, stateMPC_rilb13, &info->dgap, &info->res_ineq);
stateMPC_LA_VSUBADD2_4(stateMPC_z13, stateMPC_ubIdx13, params->ub14, stateMPC_lub13, stateMPC_sub13, stateMPC_riub13, &info->dgap, &info->res_ineq);
stateMPC_LA_VSUBADD3_2(params->lb15, stateMPC_z14, stateMPC_lbIdx14, stateMPC_llb14, stateMPC_slb14, stateMPC_rilb14, &info->dgap, &info->res_ineq);
stateMPC_LA_VSUBADD2_2(stateMPC_z14, stateMPC_ubIdx14, params->ub15, stateMPC_lub14, stateMPC_sub14, stateMPC_riub14, &info->dgap, &info->res_ineq);
stateMPC_LA_INEQ_B_GRAD_4_4_4(stateMPC_lub00, stateMPC_sub00, stateMPC_riub00, stateMPC_llb00, stateMPC_slb00, stateMPC_rilb00, stateMPC_lbIdx00, stateMPC_ubIdx00, stateMPC_grad_ineq00, stateMPC_lubbysub00, stateMPC_llbbyslb00);
stateMPC_LA_INEQ_B_GRAD_4_4_4(stateMPC_lub01, stateMPC_sub01, stateMPC_riub01, stateMPC_llb01, stateMPC_slb01, stateMPC_rilb01, stateMPC_lbIdx01, stateMPC_ubIdx01, stateMPC_grad_ineq01, stateMPC_lubbysub01, stateMPC_llbbyslb01);
stateMPC_LA_INEQ_B_GRAD_4_4_4(stateMPC_lub02, stateMPC_sub02, stateMPC_riub02, stateMPC_llb02, stateMPC_slb02, stateMPC_rilb02, stateMPC_lbIdx02, stateMPC_ubIdx02, stateMPC_grad_ineq02, stateMPC_lubbysub02, stateMPC_llbbyslb02);
stateMPC_LA_INEQ_B_GRAD_4_4_4(stateMPC_lub03, stateMPC_sub03, stateMPC_riub03, stateMPC_llb03, stateMPC_slb03, stateMPC_rilb03, stateMPC_lbIdx03, stateMPC_ubIdx03, stateMPC_grad_ineq03, stateMPC_lubbysub03, stateMPC_llbbyslb03);
stateMPC_LA_INEQ_B_GRAD_4_4_4(stateMPC_lub04, stateMPC_sub04, stateMPC_riub04, stateMPC_llb04, stateMPC_slb04, stateMPC_rilb04, stateMPC_lbIdx04, stateMPC_ubIdx04, stateMPC_grad_ineq04, stateMPC_lubbysub04, stateMPC_llbbyslb04);
stateMPC_LA_INEQ_B_GRAD_4_4_4(stateMPC_lub05, stateMPC_sub05, stateMPC_riub05, stateMPC_llb05, stateMPC_slb05, stateMPC_rilb05, stateMPC_lbIdx05, stateMPC_ubIdx05, stateMPC_grad_ineq05, stateMPC_lubbysub05, stateMPC_llbbyslb05);
stateMPC_LA_INEQ_B_GRAD_4_4_4(stateMPC_lub06, stateMPC_sub06, stateMPC_riub06, stateMPC_llb06, stateMPC_slb06, stateMPC_rilb06, stateMPC_lbIdx06, stateMPC_ubIdx06, stateMPC_grad_ineq06, stateMPC_lubbysub06, stateMPC_llbbyslb06);
stateMPC_LA_INEQ_B_GRAD_4_4_4(stateMPC_lub07, stateMPC_sub07, stateMPC_riub07, stateMPC_llb07, stateMPC_slb07, stateMPC_rilb07, stateMPC_lbIdx07, stateMPC_ubIdx07, stateMPC_grad_ineq07, stateMPC_lubbysub07, stateMPC_llbbyslb07);
stateMPC_LA_INEQ_B_GRAD_4_4_4(stateMPC_lub08, stateMPC_sub08, stateMPC_riub08, stateMPC_llb08, stateMPC_slb08, stateMPC_rilb08, stateMPC_lbIdx08, stateMPC_ubIdx08, stateMPC_grad_ineq08, stateMPC_lubbysub08, stateMPC_llbbyslb08);
stateMPC_LA_INEQ_B_GRAD_4_4_4(stateMPC_lub09, stateMPC_sub09, stateMPC_riub09, stateMPC_llb09, stateMPC_slb09, stateMPC_rilb09, stateMPC_lbIdx09, stateMPC_ubIdx09, stateMPC_grad_ineq09, stateMPC_lubbysub09, stateMPC_llbbyslb09);
stateMPC_LA_INEQ_B_GRAD_4_4_4(stateMPC_lub10, stateMPC_sub10, stateMPC_riub10, stateMPC_llb10, stateMPC_slb10, stateMPC_rilb10, stateMPC_lbIdx10, stateMPC_ubIdx10, stateMPC_grad_ineq10, stateMPC_lubbysub10, stateMPC_llbbyslb10);
stateMPC_LA_INEQ_B_GRAD_4_4_4(stateMPC_lub11, stateMPC_sub11, stateMPC_riub11, stateMPC_llb11, stateMPC_slb11, stateMPC_rilb11, stateMPC_lbIdx11, stateMPC_ubIdx11, stateMPC_grad_ineq11, stateMPC_lubbysub11, stateMPC_llbbyslb11);
stateMPC_LA_INEQ_B_GRAD_4_4_4(stateMPC_lub12, stateMPC_sub12, stateMPC_riub12, stateMPC_llb12, stateMPC_slb12, stateMPC_rilb12, stateMPC_lbIdx12, stateMPC_ubIdx12, stateMPC_grad_ineq12, stateMPC_lubbysub12, stateMPC_llbbyslb12);
stateMPC_LA_INEQ_B_GRAD_4_4_4(stateMPC_lub13, stateMPC_sub13, stateMPC_riub13, stateMPC_llb13, stateMPC_slb13, stateMPC_rilb13, stateMPC_lbIdx13, stateMPC_ubIdx13, stateMPC_grad_ineq13, stateMPC_lubbysub13, stateMPC_llbbyslb13);
stateMPC_LA_INEQ_B_GRAD_2_2_2(stateMPC_lub14, stateMPC_sub14, stateMPC_riub14, stateMPC_llb14, stateMPC_slb14, stateMPC_rilb14, stateMPC_lbIdx14, stateMPC_ubIdx14, stateMPC_grad_ineq14, stateMPC_lubbysub14, stateMPC_llbbyslb14);
info->dobj = info->pobj - info->dgap;
info->rdgap = info->pobj ? info->dgap / info->pobj : 1e6;
if( info->rdgap < 0 ) info->rdgap = -info->rdgap;
if( info->mu < stateMPC_SET_ACC_KKTCOMPL
    && (info->rdgap < stateMPC_SET_ACC_RDGAP || info->dgap < stateMPC_SET_ACC_KKTCOMPL)
    && info->res_eq < stateMPC_SET_ACC_RESEQ
    && info->res_ineq < stateMPC_SET_ACC_RESINEQ ){
exitcode = stateMPC_OPTIMAL; break; }
if( info->it == stateMPC_SET_MAXIT ){
exitcode = stateMPC_MAXITREACHED; break; }
stateMPC_LA_VVADD3_58(stateMPC_grad_cost, stateMPC_grad_eq, stateMPC_grad_ineq, stateMPC_rd);
stateMPC_LA_DIAG_CHOL_ONELOOP_LBUB_4_4_4(params->H1, stateMPC_llbbyslb00, stateMPC_lbIdx00, stateMPC_lubbysub00, stateMPC_ubIdx00, stateMPC_Phi00);
stateMPC_LA_DIAG_MATRIXFORWARDSUB_2_4(stateMPC_Phi00, params->C1, stateMPC_V00);
stateMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_4(stateMPC_Phi00, stateMPC_D00, stateMPC_W00);
stateMPC_LA_DENSE_DIAGZERO_MMTM_2_4_2(stateMPC_W00, stateMPC_V00, stateMPC_Ysd01);
stateMPC_LA_DIAG_FORWARDSUB_4(stateMPC_Phi00, stateMPC_rd00, stateMPC_Lbyrd00);
stateMPC_LA_DIAG_CHOL_ONELOOP_LBUB_4_4_4(params->H2, stateMPC_llbbyslb01, stateMPC_lbIdx01, stateMPC_lubbysub01, stateMPC_ubIdx01, stateMPC_Phi01);
stateMPC_LA_DIAG_MATRIXFORWARDSUB_2_4(stateMPC_Phi01, params->C2, stateMPC_V01);
stateMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_4(stateMPC_Phi01, stateMPC_D01, stateMPC_W01);
stateMPC_LA_DENSE_DIAGZERO_MMTM_2_4_2(stateMPC_W01, stateMPC_V01, stateMPC_Ysd02);
stateMPC_LA_DIAG_FORWARDSUB_4(stateMPC_Phi01, stateMPC_rd01, stateMPC_Lbyrd01);
stateMPC_LA_DIAG_CHOL_ONELOOP_LBUB_4_4_4(params->H3, stateMPC_llbbyslb02, stateMPC_lbIdx02, stateMPC_lubbysub02, stateMPC_ubIdx02, stateMPC_Phi02);
stateMPC_LA_DIAG_MATRIXFORWARDSUB_2_4(stateMPC_Phi02, params->C3, stateMPC_V02);
stateMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_4(stateMPC_Phi02, stateMPC_D01, stateMPC_W02);
stateMPC_LA_DENSE_DIAGZERO_MMTM_2_4_2(stateMPC_W02, stateMPC_V02, stateMPC_Ysd03);
stateMPC_LA_DIAG_FORWARDSUB_4(stateMPC_Phi02, stateMPC_rd02, stateMPC_Lbyrd02);
stateMPC_LA_DIAG_CHOL_ONELOOP_LBUB_4_4_4(params->H4, stateMPC_llbbyslb03, stateMPC_lbIdx03, stateMPC_lubbysub03, stateMPC_ubIdx03, stateMPC_Phi03);
stateMPC_LA_DIAG_MATRIXFORWARDSUB_2_4(stateMPC_Phi03, params->C4, stateMPC_V03);
stateMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_4(stateMPC_Phi03, stateMPC_D01, stateMPC_W03);
stateMPC_LA_DENSE_DIAGZERO_MMTM_2_4_2(stateMPC_W03, stateMPC_V03, stateMPC_Ysd04);
stateMPC_LA_DIAG_FORWARDSUB_4(stateMPC_Phi03, stateMPC_rd03, stateMPC_Lbyrd03);
stateMPC_LA_DIAG_CHOL_ONELOOP_LBUB_4_4_4(params->H5, stateMPC_llbbyslb04, stateMPC_lbIdx04, stateMPC_lubbysub04, stateMPC_ubIdx04, stateMPC_Phi04);
stateMPC_LA_DIAG_MATRIXFORWARDSUB_2_4(stateMPC_Phi04, params->C5, stateMPC_V04);
stateMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_4(stateMPC_Phi04, stateMPC_D01, stateMPC_W04);
stateMPC_LA_DENSE_DIAGZERO_MMTM_2_4_2(stateMPC_W04, stateMPC_V04, stateMPC_Ysd05);
stateMPC_LA_DIAG_FORWARDSUB_4(stateMPC_Phi04, stateMPC_rd04, stateMPC_Lbyrd04);
stateMPC_LA_DIAG_CHOL_ONELOOP_LBUB_4_4_4(params->H6, stateMPC_llbbyslb05, stateMPC_lbIdx05, stateMPC_lubbysub05, stateMPC_ubIdx05, stateMPC_Phi05);
stateMPC_LA_DIAG_MATRIXFORWARDSUB_2_4(stateMPC_Phi05, params->C6, stateMPC_V05);
stateMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_4(stateMPC_Phi05, stateMPC_D01, stateMPC_W05);
stateMPC_LA_DENSE_DIAGZERO_MMTM_2_4_2(stateMPC_W05, stateMPC_V05, stateMPC_Ysd06);
stateMPC_LA_DIAG_FORWARDSUB_4(stateMPC_Phi05, stateMPC_rd05, stateMPC_Lbyrd05);
stateMPC_LA_DIAG_CHOL_ONELOOP_LBUB_4_4_4(params->H7, stateMPC_llbbyslb06, stateMPC_lbIdx06, stateMPC_lubbysub06, stateMPC_ubIdx06, stateMPC_Phi06);
stateMPC_LA_DIAG_MATRIXFORWARDSUB_2_4(stateMPC_Phi06, params->C7, stateMPC_V06);
stateMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_4(stateMPC_Phi06, stateMPC_D01, stateMPC_W06);
stateMPC_LA_DENSE_DIAGZERO_MMTM_2_4_2(stateMPC_W06, stateMPC_V06, stateMPC_Ysd07);
stateMPC_LA_DIAG_FORWARDSUB_4(stateMPC_Phi06, stateMPC_rd06, stateMPC_Lbyrd06);
stateMPC_LA_DIAG_CHOL_ONELOOP_LBUB_4_4_4(params->H8, stateMPC_llbbyslb07, stateMPC_lbIdx07, stateMPC_lubbysub07, stateMPC_ubIdx07, stateMPC_Phi07);
stateMPC_LA_DIAG_MATRIXFORWARDSUB_2_4(stateMPC_Phi07, params->C8, stateMPC_V07);
stateMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_4(stateMPC_Phi07, stateMPC_D01, stateMPC_W07);
stateMPC_LA_DENSE_DIAGZERO_MMTM_2_4_2(stateMPC_W07, stateMPC_V07, stateMPC_Ysd08);
stateMPC_LA_DIAG_FORWARDSUB_4(stateMPC_Phi07, stateMPC_rd07, stateMPC_Lbyrd07);
stateMPC_LA_DIAG_CHOL_ONELOOP_LBUB_4_4_4(params->H9, stateMPC_llbbyslb08, stateMPC_lbIdx08, stateMPC_lubbysub08, stateMPC_ubIdx08, stateMPC_Phi08);
stateMPC_LA_DIAG_MATRIXFORWARDSUB_2_4(stateMPC_Phi08, params->C9, stateMPC_V08);
stateMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_4(stateMPC_Phi08, stateMPC_D01, stateMPC_W08);
stateMPC_LA_DENSE_DIAGZERO_MMTM_2_4_2(stateMPC_W08, stateMPC_V08, stateMPC_Ysd09);
stateMPC_LA_DIAG_FORWARDSUB_4(stateMPC_Phi08, stateMPC_rd08, stateMPC_Lbyrd08);
stateMPC_LA_DIAG_CHOL_ONELOOP_LBUB_4_4_4(params->H10, stateMPC_llbbyslb09, stateMPC_lbIdx09, stateMPC_lubbysub09, stateMPC_ubIdx09, stateMPC_Phi09);
stateMPC_LA_DIAG_MATRIXFORWARDSUB_2_4(stateMPC_Phi09, params->C10, stateMPC_V09);
stateMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_4(stateMPC_Phi09, stateMPC_D01, stateMPC_W09);
stateMPC_LA_DENSE_DIAGZERO_MMTM_2_4_2(stateMPC_W09, stateMPC_V09, stateMPC_Ysd10);
stateMPC_LA_DIAG_FORWARDSUB_4(stateMPC_Phi09, stateMPC_rd09, stateMPC_Lbyrd09);
stateMPC_LA_DIAG_CHOL_ONELOOP_LBUB_4_4_4(params->H11, stateMPC_llbbyslb10, stateMPC_lbIdx10, stateMPC_lubbysub10, stateMPC_ubIdx10, stateMPC_Phi10);
stateMPC_LA_DIAG_MATRIXFORWARDSUB_2_4(stateMPC_Phi10, params->C11, stateMPC_V10);
stateMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_4(stateMPC_Phi10, stateMPC_D01, stateMPC_W10);
stateMPC_LA_DENSE_DIAGZERO_MMTM_2_4_2(stateMPC_W10, stateMPC_V10, stateMPC_Ysd11);
stateMPC_LA_DIAG_FORWARDSUB_4(stateMPC_Phi10, stateMPC_rd10, stateMPC_Lbyrd10);
stateMPC_LA_DIAG_CHOL_ONELOOP_LBUB_4_4_4(params->H12, stateMPC_llbbyslb11, stateMPC_lbIdx11, stateMPC_lubbysub11, stateMPC_ubIdx11, stateMPC_Phi11);
stateMPC_LA_DIAG_MATRIXFORWARDSUB_2_4(stateMPC_Phi11, params->C12, stateMPC_V11);
stateMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_4(stateMPC_Phi11, stateMPC_D01, stateMPC_W11);
stateMPC_LA_DENSE_DIAGZERO_MMTM_2_4_2(stateMPC_W11, stateMPC_V11, stateMPC_Ysd12);
stateMPC_LA_DIAG_FORWARDSUB_4(stateMPC_Phi11, stateMPC_rd11, stateMPC_Lbyrd11);
stateMPC_LA_DIAG_CHOL_ONELOOP_LBUB_4_4_4(params->H13, stateMPC_llbbyslb12, stateMPC_lbIdx12, stateMPC_lubbysub12, stateMPC_ubIdx12, stateMPC_Phi12);
stateMPC_LA_DIAG_MATRIXFORWARDSUB_2_4(stateMPC_Phi12, params->C13, stateMPC_V12);
stateMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_4(stateMPC_Phi12, stateMPC_D01, stateMPC_W12);
stateMPC_LA_DENSE_DIAGZERO_MMTM_2_4_2(stateMPC_W12, stateMPC_V12, stateMPC_Ysd13);
stateMPC_LA_DIAG_FORWARDSUB_4(stateMPC_Phi12, stateMPC_rd12, stateMPC_Lbyrd12);
stateMPC_LA_DIAG_CHOL_ONELOOP_LBUB_4_4_4(params->H14, stateMPC_llbbyslb13, stateMPC_lbIdx13, stateMPC_lubbysub13, stateMPC_ubIdx13, stateMPC_Phi13);
stateMPC_LA_DIAG_MATRIXFORWARDSUB_2_4(stateMPC_Phi13, params->C14, stateMPC_V13);
stateMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_4(stateMPC_Phi13, stateMPC_D01, stateMPC_W13);
stateMPC_LA_DENSE_DIAGZERO_MMTM_2_4_2(stateMPC_W13, stateMPC_V13, stateMPC_Ysd14);
stateMPC_LA_DIAG_FORWARDSUB_4(stateMPC_Phi13, stateMPC_rd13, stateMPC_Lbyrd13);
stateMPC_LA_DIAG_CHOL_ONELOOP_LBUB_2_2_2(params->H15, stateMPC_llbbyslb14, stateMPC_lbIdx14, stateMPC_lubbysub14, stateMPC_ubIdx14, stateMPC_Phi14);
stateMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_2(stateMPC_Phi14, stateMPC_D14, stateMPC_W14);
stateMPC_LA_DIAG_FORWARDSUB_2(stateMPC_Phi14, stateMPC_rd14, stateMPC_Lbyrd14);
stateMPC_LA_DIAGZERO_MMT_2(stateMPC_W00, stateMPC_Yd00);
stateMPC_LA_DIAGZERO_MVMSUB7_2(stateMPC_W00, stateMPC_Lbyrd00, stateMPC_re00, stateMPC_beta00);
stateMPC_LA_DENSE_DIAGZERO_MMT2_2_4_4(stateMPC_V00, stateMPC_W01, stateMPC_Yd01);
stateMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_4_4(stateMPC_V00, stateMPC_Lbyrd00, stateMPC_W01, stateMPC_Lbyrd01, stateMPC_re01, stateMPC_beta01);
stateMPC_LA_DENSE_DIAGZERO_MMT2_2_4_4(stateMPC_V01, stateMPC_W02, stateMPC_Yd02);
stateMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_4_4(stateMPC_V01, stateMPC_Lbyrd01, stateMPC_W02, stateMPC_Lbyrd02, stateMPC_re02, stateMPC_beta02);
stateMPC_LA_DENSE_DIAGZERO_MMT2_2_4_4(stateMPC_V02, stateMPC_W03, stateMPC_Yd03);
stateMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_4_4(stateMPC_V02, stateMPC_Lbyrd02, stateMPC_W03, stateMPC_Lbyrd03, stateMPC_re03, stateMPC_beta03);
stateMPC_LA_DENSE_DIAGZERO_MMT2_2_4_4(stateMPC_V03, stateMPC_W04, stateMPC_Yd04);
stateMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_4_4(stateMPC_V03, stateMPC_Lbyrd03, stateMPC_W04, stateMPC_Lbyrd04, stateMPC_re04, stateMPC_beta04);
stateMPC_LA_DENSE_DIAGZERO_MMT2_2_4_4(stateMPC_V04, stateMPC_W05, stateMPC_Yd05);
stateMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_4_4(stateMPC_V04, stateMPC_Lbyrd04, stateMPC_W05, stateMPC_Lbyrd05, stateMPC_re05, stateMPC_beta05);
stateMPC_LA_DENSE_DIAGZERO_MMT2_2_4_4(stateMPC_V05, stateMPC_W06, stateMPC_Yd06);
stateMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_4_4(stateMPC_V05, stateMPC_Lbyrd05, stateMPC_W06, stateMPC_Lbyrd06, stateMPC_re06, stateMPC_beta06);
stateMPC_LA_DENSE_DIAGZERO_MMT2_2_4_4(stateMPC_V06, stateMPC_W07, stateMPC_Yd07);
stateMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_4_4(stateMPC_V06, stateMPC_Lbyrd06, stateMPC_W07, stateMPC_Lbyrd07, stateMPC_re07, stateMPC_beta07);
stateMPC_LA_DENSE_DIAGZERO_MMT2_2_4_4(stateMPC_V07, stateMPC_W08, stateMPC_Yd08);
stateMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_4_4(stateMPC_V07, stateMPC_Lbyrd07, stateMPC_W08, stateMPC_Lbyrd08, stateMPC_re08, stateMPC_beta08);
stateMPC_LA_DENSE_DIAGZERO_MMT2_2_4_4(stateMPC_V08, stateMPC_W09, stateMPC_Yd09);
stateMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_4_4(stateMPC_V08, stateMPC_Lbyrd08, stateMPC_W09, stateMPC_Lbyrd09, stateMPC_re09, stateMPC_beta09);
stateMPC_LA_DENSE_DIAGZERO_MMT2_2_4_4(stateMPC_V09, stateMPC_W10, stateMPC_Yd10);
stateMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_4_4(stateMPC_V09, stateMPC_Lbyrd09, stateMPC_W10, stateMPC_Lbyrd10, stateMPC_re10, stateMPC_beta10);
stateMPC_LA_DENSE_DIAGZERO_MMT2_2_4_4(stateMPC_V10, stateMPC_W11, stateMPC_Yd11);
stateMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_4_4(stateMPC_V10, stateMPC_Lbyrd10, stateMPC_W11, stateMPC_Lbyrd11, stateMPC_re11, stateMPC_beta11);
stateMPC_LA_DENSE_DIAGZERO_MMT2_2_4_4(stateMPC_V11, stateMPC_W12, stateMPC_Yd12);
stateMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_4_4(stateMPC_V11, stateMPC_Lbyrd11, stateMPC_W12, stateMPC_Lbyrd12, stateMPC_re12, stateMPC_beta12);
stateMPC_LA_DENSE_DIAGZERO_MMT2_2_4_4(stateMPC_V12, stateMPC_W13, stateMPC_Yd13);
stateMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_4_4(stateMPC_V12, stateMPC_Lbyrd12, stateMPC_W13, stateMPC_Lbyrd13, stateMPC_re13, stateMPC_beta13);
stateMPC_LA_DENSE_DIAGZERO_MMT2_2_4_2(stateMPC_V13, stateMPC_W14, stateMPC_Yd14);
stateMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_4_2(stateMPC_V13, stateMPC_Lbyrd13, stateMPC_W14, stateMPC_Lbyrd14, stateMPC_re14, stateMPC_beta14);
stateMPC_LA_DENSE_CHOL_2(stateMPC_Yd00, stateMPC_Ld00);
stateMPC_LA_DENSE_FORWARDSUB_2(stateMPC_Ld00, stateMPC_beta00, stateMPC_yy00);
stateMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(stateMPC_Ld00, stateMPC_Ysd01, stateMPC_Lsd01);
stateMPC_LA_DENSE_MMTSUB_2_2(stateMPC_Lsd01, stateMPC_Yd01);
stateMPC_LA_DENSE_CHOL_2(stateMPC_Yd01, stateMPC_Ld01);
stateMPC_LA_DENSE_MVMSUB1_2_2(stateMPC_Lsd01, stateMPC_yy00, stateMPC_beta01, stateMPC_bmy01);
stateMPC_LA_DENSE_FORWARDSUB_2(stateMPC_Ld01, stateMPC_bmy01, stateMPC_yy01);
stateMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(stateMPC_Ld01, stateMPC_Ysd02, stateMPC_Lsd02);
stateMPC_LA_DENSE_MMTSUB_2_2(stateMPC_Lsd02, stateMPC_Yd02);
stateMPC_LA_DENSE_CHOL_2(stateMPC_Yd02, stateMPC_Ld02);
stateMPC_LA_DENSE_MVMSUB1_2_2(stateMPC_Lsd02, stateMPC_yy01, stateMPC_beta02, stateMPC_bmy02);
stateMPC_LA_DENSE_FORWARDSUB_2(stateMPC_Ld02, stateMPC_bmy02, stateMPC_yy02);
stateMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(stateMPC_Ld02, stateMPC_Ysd03, stateMPC_Lsd03);
stateMPC_LA_DENSE_MMTSUB_2_2(stateMPC_Lsd03, stateMPC_Yd03);
stateMPC_LA_DENSE_CHOL_2(stateMPC_Yd03, stateMPC_Ld03);
stateMPC_LA_DENSE_MVMSUB1_2_2(stateMPC_Lsd03, stateMPC_yy02, stateMPC_beta03, stateMPC_bmy03);
stateMPC_LA_DENSE_FORWARDSUB_2(stateMPC_Ld03, stateMPC_bmy03, stateMPC_yy03);
stateMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(stateMPC_Ld03, stateMPC_Ysd04, stateMPC_Lsd04);
stateMPC_LA_DENSE_MMTSUB_2_2(stateMPC_Lsd04, stateMPC_Yd04);
stateMPC_LA_DENSE_CHOL_2(stateMPC_Yd04, stateMPC_Ld04);
stateMPC_LA_DENSE_MVMSUB1_2_2(stateMPC_Lsd04, stateMPC_yy03, stateMPC_beta04, stateMPC_bmy04);
stateMPC_LA_DENSE_FORWARDSUB_2(stateMPC_Ld04, stateMPC_bmy04, stateMPC_yy04);
stateMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(stateMPC_Ld04, stateMPC_Ysd05, stateMPC_Lsd05);
stateMPC_LA_DENSE_MMTSUB_2_2(stateMPC_Lsd05, stateMPC_Yd05);
stateMPC_LA_DENSE_CHOL_2(stateMPC_Yd05, stateMPC_Ld05);
stateMPC_LA_DENSE_MVMSUB1_2_2(stateMPC_Lsd05, stateMPC_yy04, stateMPC_beta05, stateMPC_bmy05);
stateMPC_LA_DENSE_FORWARDSUB_2(stateMPC_Ld05, stateMPC_bmy05, stateMPC_yy05);
stateMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(stateMPC_Ld05, stateMPC_Ysd06, stateMPC_Lsd06);
stateMPC_LA_DENSE_MMTSUB_2_2(stateMPC_Lsd06, stateMPC_Yd06);
stateMPC_LA_DENSE_CHOL_2(stateMPC_Yd06, stateMPC_Ld06);
stateMPC_LA_DENSE_MVMSUB1_2_2(stateMPC_Lsd06, stateMPC_yy05, stateMPC_beta06, stateMPC_bmy06);
stateMPC_LA_DENSE_FORWARDSUB_2(stateMPC_Ld06, stateMPC_bmy06, stateMPC_yy06);
stateMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(stateMPC_Ld06, stateMPC_Ysd07, stateMPC_Lsd07);
stateMPC_LA_DENSE_MMTSUB_2_2(stateMPC_Lsd07, stateMPC_Yd07);
stateMPC_LA_DENSE_CHOL_2(stateMPC_Yd07, stateMPC_Ld07);
stateMPC_LA_DENSE_MVMSUB1_2_2(stateMPC_Lsd07, stateMPC_yy06, stateMPC_beta07, stateMPC_bmy07);
stateMPC_LA_DENSE_FORWARDSUB_2(stateMPC_Ld07, stateMPC_bmy07, stateMPC_yy07);
stateMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(stateMPC_Ld07, stateMPC_Ysd08, stateMPC_Lsd08);
stateMPC_LA_DENSE_MMTSUB_2_2(stateMPC_Lsd08, stateMPC_Yd08);
stateMPC_LA_DENSE_CHOL_2(stateMPC_Yd08, stateMPC_Ld08);
stateMPC_LA_DENSE_MVMSUB1_2_2(stateMPC_Lsd08, stateMPC_yy07, stateMPC_beta08, stateMPC_bmy08);
stateMPC_LA_DENSE_FORWARDSUB_2(stateMPC_Ld08, stateMPC_bmy08, stateMPC_yy08);
stateMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(stateMPC_Ld08, stateMPC_Ysd09, stateMPC_Lsd09);
stateMPC_LA_DENSE_MMTSUB_2_2(stateMPC_Lsd09, stateMPC_Yd09);
stateMPC_LA_DENSE_CHOL_2(stateMPC_Yd09, stateMPC_Ld09);
stateMPC_LA_DENSE_MVMSUB1_2_2(stateMPC_Lsd09, stateMPC_yy08, stateMPC_beta09, stateMPC_bmy09);
stateMPC_LA_DENSE_FORWARDSUB_2(stateMPC_Ld09, stateMPC_bmy09, stateMPC_yy09);
stateMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(stateMPC_Ld09, stateMPC_Ysd10, stateMPC_Lsd10);
stateMPC_LA_DENSE_MMTSUB_2_2(stateMPC_Lsd10, stateMPC_Yd10);
stateMPC_LA_DENSE_CHOL_2(stateMPC_Yd10, stateMPC_Ld10);
stateMPC_LA_DENSE_MVMSUB1_2_2(stateMPC_Lsd10, stateMPC_yy09, stateMPC_beta10, stateMPC_bmy10);
stateMPC_LA_DENSE_FORWARDSUB_2(stateMPC_Ld10, stateMPC_bmy10, stateMPC_yy10);
stateMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(stateMPC_Ld10, stateMPC_Ysd11, stateMPC_Lsd11);
stateMPC_LA_DENSE_MMTSUB_2_2(stateMPC_Lsd11, stateMPC_Yd11);
stateMPC_LA_DENSE_CHOL_2(stateMPC_Yd11, stateMPC_Ld11);
stateMPC_LA_DENSE_MVMSUB1_2_2(stateMPC_Lsd11, stateMPC_yy10, stateMPC_beta11, stateMPC_bmy11);
stateMPC_LA_DENSE_FORWARDSUB_2(stateMPC_Ld11, stateMPC_bmy11, stateMPC_yy11);
stateMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(stateMPC_Ld11, stateMPC_Ysd12, stateMPC_Lsd12);
stateMPC_LA_DENSE_MMTSUB_2_2(stateMPC_Lsd12, stateMPC_Yd12);
stateMPC_LA_DENSE_CHOL_2(stateMPC_Yd12, stateMPC_Ld12);
stateMPC_LA_DENSE_MVMSUB1_2_2(stateMPC_Lsd12, stateMPC_yy11, stateMPC_beta12, stateMPC_bmy12);
stateMPC_LA_DENSE_FORWARDSUB_2(stateMPC_Ld12, stateMPC_bmy12, stateMPC_yy12);
stateMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(stateMPC_Ld12, stateMPC_Ysd13, stateMPC_Lsd13);
stateMPC_LA_DENSE_MMTSUB_2_2(stateMPC_Lsd13, stateMPC_Yd13);
stateMPC_LA_DENSE_CHOL_2(stateMPC_Yd13, stateMPC_Ld13);
stateMPC_LA_DENSE_MVMSUB1_2_2(stateMPC_Lsd13, stateMPC_yy12, stateMPC_beta13, stateMPC_bmy13);
stateMPC_LA_DENSE_FORWARDSUB_2(stateMPC_Ld13, stateMPC_bmy13, stateMPC_yy13);
stateMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(stateMPC_Ld13, stateMPC_Ysd14, stateMPC_Lsd14);
stateMPC_LA_DENSE_MMTSUB_2_2(stateMPC_Lsd14, stateMPC_Yd14);
stateMPC_LA_DENSE_CHOL_2(stateMPC_Yd14, stateMPC_Ld14);
stateMPC_LA_DENSE_MVMSUB1_2_2(stateMPC_Lsd14, stateMPC_yy13, stateMPC_beta14, stateMPC_bmy14);
stateMPC_LA_DENSE_FORWARDSUB_2(stateMPC_Ld14, stateMPC_bmy14, stateMPC_yy14);
stateMPC_LA_DENSE_BACKWARDSUB_2(stateMPC_Ld14, stateMPC_yy14, stateMPC_dvaff14);
stateMPC_LA_DENSE_MTVMSUB_2_2(stateMPC_Lsd14, stateMPC_dvaff14, stateMPC_yy13, stateMPC_bmy13);
stateMPC_LA_DENSE_BACKWARDSUB_2(stateMPC_Ld13, stateMPC_bmy13, stateMPC_dvaff13);
stateMPC_LA_DENSE_MTVMSUB_2_2(stateMPC_Lsd13, stateMPC_dvaff13, stateMPC_yy12, stateMPC_bmy12);
stateMPC_LA_DENSE_BACKWARDSUB_2(stateMPC_Ld12, stateMPC_bmy12, stateMPC_dvaff12);
stateMPC_LA_DENSE_MTVMSUB_2_2(stateMPC_Lsd12, stateMPC_dvaff12, stateMPC_yy11, stateMPC_bmy11);
stateMPC_LA_DENSE_BACKWARDSUB_2(stateMPC_Ld11, stateMPC_bmy11, stateMPC_dvaff11);
stateMPC_LA_DENSE_MTVMSUB_2_2(stateMPC_Lsd11, stateMPC_dvaff11, stateMPC_yy10, stateMPC_bmy10);
stateMPC_LA_DENSE_BACKWARDSUB_2(stateMPC_Ld10, stateMPC_bmy10, stateMPC_dvaff10);
stateMPC_LA_DENSE_MTVMSUB_2_2(stateMPC_Lsd10, stateMPC_dvaff10, stateMPC_yy09, stateMPC_bmy09);
stateMPC_LA_DENSE_BACKWARDSUB_2(stateMPC_Ld09, stateMPC_bmy09, stateMPC_dvaff09);
stateMPC_LA_DENSE_MTVMSUB_2_2(stateMPC_Lsd09, stateMPC_dvaff09, stateMPC_yy08, stateMPC_bmy08);
stateMPC_LA_DENSE_BACKWARDSUB_2(stateMPC_Ld08, stateMPC_bmy08, stateMPC_dvaff08);
stateMPC_LA_DENSE_MTVMSUB_2_2(stateMPC_Lsd08, stateMPC_dvaff08, stateMPC_yy07, stateMPC_bmy07);
stateMPC_LA_DENSE_BACKWARDSUB_2(stateMPC_Ld07, stateMPC_bmy07, stateMPC_dvaff07);
stateMPC_LA_DENSE_MTVMSUB_2_2(stateMPC_Lsd07, stateMPC_dvaff07, stateMPC_yy06, stateMPC_bmy06);
stateMPC_LA_DENSE_BACKWARDSUB_2(stateMPC_Ld06, stateMPC_bmy06, stateMPC_dvaff06);
stateMPC_LA_DENSE_MTVMSUB_2_2(stateMPC_Lsd06, stateMPC_dvaff06, stateMPC_yy05, stateMPC_bmy05);
stateMPC_LA_DENSE_BACKWARDSUB_2(stateMPC_Ld05, stateMPC_bmy05, stateMPC_dvaff05);
stateMPC_LA_DENSE_MTVMSUB_2_2(stateMPC_Lsd05, stateMPC_dvaff05, stateMPC_yy04, stateMPC_bmy04);
stateMPC_LA_DENSE_BACKWARDSUB_2(stateMPC_Ld04, stateMPC_bmy04, stateMPC_dvaff04);
stateMPC_LA_DENSE_MTVMSUB_2_2(stateMPC_Lsd04, stateMPC_dvaff04, stateMPC_yy03, stateMPC_bmy03);
stateMPC_LA_DENSE_BACKWARDSUB_2(stateMPC_Ld03, stateMPC_bmy03, stateMPC_dvaff03);
stateMPC_LA_DENSE_MTVMSUB_2_2(stateMPC_Lsd03, stateMPC_dvaff03, stateMPC_yy02, stateMPC_bmy02);
stateMPC_LA_DENSE_BACKWARDSUB_2(stateMPC_Ld02, stateMPC_bmy02, stateMPC_dvaff02);
stateMPC_LA_DENSE_MTVMSUB_2_2(stateMPC_Lsd02, stateMPC_dvaff02, stateMPC_yy01, stateMPC_bmy01);
stateMPC_LA_DENSE_BACKWARDSUB_2(stateMPC_Ld01, stateMPC_bmy01, stateMPC_dvaff01);
stateMPC_LA_DENSE_MTVMSUB_2_2(stateMPC_Lsd01, stateMPC_dvaff01, stateMPC_yy00, stateMPC_bmy00);
stateMPC_LA_DENSE_BACKWARDSUB_2(stateMPC_Ld00, stateMPC_bmy00, stateMPC_dvaff00);
stateMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C1, stateMPC_dvaff01, stateMPC_D00, stateMPC_dvaff00, stateMPC_grad_eq00);
stateMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C2, stateMPC_dvaff02, stateMPC_D01, stateMPC_dvaff01, stateMPC_grad_eq01);
stateMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C3, stateMPC_dvaff03, stateMPC_D01, stateMPC_dvaff02, stateMPC_grad_eq02);
stateMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C4, stateMPC_dvaff04, stateMPC_D01, stateMPC_dvaff03, stateMPC_grad_eq03);
stateMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C5, stateMPC_dvaff05, stateMPC_D01, stateMPC_dvaff04, stateMPC_grad_eq04);
stateMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C6, stateMPC_dvaff06, stateMPC_D01, stateMPC_dvaff05, stateMPC_grad_eq05);
stateMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C7, stateMPC_dvaff07, stateMPC_D01, stateMPC_dvaff06, stateMPC_grad_eq06);
stateMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C8, stateMPC_dvaff08, stateMPC_D01, stateMPC_dvaff07, stateMPC_grad_eq07);
stateMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C9, stateMPC_dvaff09, stateMPC_D01, stateMPC_dvaff08, stateMPC_grad_eq08);
stateMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C10, stateMPC_dvaff10, stateMPC_D01, stateMPC_dvaff09, stateMPC_grad_eq09);
stateMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C11, stateMPC_dvaff11, stateMPC_D01, stateMPC_dvaff10, stateMPC_grad_eq10);
stateMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C12, stateMPC_dvaff12, stateMPC_D01, stateMPC_dvaff11, stateMPC_grad_eq11);
stateMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C13, stateMPC_dvaff13, stateMPC_D01, stateMPC_dvaff12, stateMPC_grad_eq12);
stateMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C14, stateMPC_dvaff14, stateMPC_D01, stateMPC_dvaff13, stateMPC_grad_eq13);
stateMPC_LA_DIAGZERO_MTVM_2_2(stateMPC_D14, stateMPC_dvaff14, stateMPC_grad_eq14);
stateMPC_LA_VSUB2_58(stateMPC_rd, stateMPC_grad_eq, stateMPC_rd);
stateMPC_LA_DIAG_FORWARDBACKWARDSUB_4(stateMPC_Phi00, stateMPC_rd00, stateMPC_dzaff00);
stateMPC_LA_DIAG_FORWARDBACKWARDSUB_4(stateMPC_Phi01, stateMPC_rd01, stateMPC_dzaff01);
stateMPC_LA_DIAG_FORWARDBACKWARDSUB_4(stateMPC_Phi02, stateMPC_rd02, stateMPC_dzaff02);
stateMPC_LA_DIAG_FORWARDBACKWARDSUB_4(stateMPC_Phi03, stateMPC_rd03, stateMPC_dzaff03);
stateMPC_LA_DIAG_FORWARDBACKWARDSUB_4(stateMPC_Phi04, stateMPC_rd04, stateMPC_dzaff04);
stateMPC_LA_DIAG_FORWARDBACKWARDSUB_4(stateMPC_Phi05, stateMPC_rd05, stateMPC_dzaff05);
stateMPC_LA_DIAG_FORWARDBACKWARDSUB_4(stateMPC_Phi06, stateMPC_rd06, stateMPC_dzaff06);
stateMPC_LA_DIAG_FORWARDBACKWARDSUB_4(stateMPC_Phi07, stateMPC_rd07, stateMPC_dzaff07);
stateMPC_LA_DIAG_FORWARDBACKWARDSUB_4(stateMPC_Phi08, stateMPC_rd08, stateMPC_dzaff08);
stateMPC_LA_DIAG_FORWARDBACKWARDSUB_4(stateMPC_Phi09, stateMPC_rd09, stateMPC_dzaff09);
stateMPC_LA_DIAG_FORWARDBACKWARDSUB_4(stateMPC_Phi10, stateMPC_rd10, stateMPC_dzaff10);
stateMPC_LA_DIAG_FORWARDBACKWARDSUB_4(stateMPC_Phi11, stateMPC_rd11, stateMPC_dzaff11);
stateMPC_LA_DIAG_FORWARDBACKWARDSUB_4(stateMPC_Phi12, stateMPC_rd12, stateMPC_dzaff12);
stateMPC_LA_DIAG_FORWARDBACKWARDSUB_4(stateMPC_Phi13, stateMPC_rd13, stateMPC_dzaff13);
stateMPC_LA_DIAG_FORWARDBACKWARDSUB_2(stateMPC_Phi14, stateMPC_rd14, stateMPC_dzaff14);
stateMPC_LA_VSUB_INDEXED_4(stateMPC_dzaff00, stateMPC_lbIdx00, stateMPC_rilb00, stateMPC_dslbaff00);
stateMPC_LA_VSUB3_4(stateMPC_llbbyslb00, stateMPC_dslbaff00, stateMPC_llb00, stateMPC_dllbaff00);
stateMPC_LA_VSUB2_INDEXED_4(stateMPC_riub00, stateMPC_dzaff00, stateMPC_ubIdx00, stateMPC_dsubaff00);
stateMPC_LA_VSUB3_4(stateMPC_lubbysub00, stateMPC_dsubaff00, stateMPC_lub00, stateMPC_dlubaff00);
stateMPC_LA_VSUB_INDEXED_4(stateMPC_dzaff01, stateMPC_lbIdx01, stateMPC_rilb01, stateMPC_dslbaff01);
stateMPC_LA_VSUB3_4(stateMPC_llbbyslb01, stateMPC_dslbaff01, stateMPC_llb01, stateMPC_dllbaff01);
stateMPC_LA_VSUB2_INDEXED_4(stateMPC_riub01, stateMPC_dzaff01, stateMPC_ubIdx01, stateMPC_dsubaff01);
stateMPC_LA_VSUB3_4(stateMPC_lubbysub01, stateMPC_dsubaff01, stateMPC_lub01, stateMPC_dlubaff01);
stateMPC_LA_VSUB_INDEXED_4(stateMPC_dzaff02, stateMPC_lbIdx02, stateMPC_rilb02, stateMPC_dslbaff02);
stateMPC_LA_VSUB3_4(stateMPC_llbbyslb02, stateMPC_dslbaff02, stateMPC_llb02, stateMPC_dllbaff02);
stateMPC_LA_VSUB2_INDEXED_4(stateMPC_riub02, stateMPC_dzaff02, stateMPC_ubIdx02, stateMPC_dsubaff02);
stateMPC_LA_VSUB3_4(stateMPC_lubbysub02, stateMPC_dsubaff02, stateMPC_lub02, stateMPC_dlubaff02);
stateMPC_LA_VSUB_INDEXED_4(stateMPC_dzaff03, stateMPC_lbIdx03, stateMPC_rilb03, stateMPC_dslbaff03);
stateMPC_LA_VSUB3_4(stateMPC_llbbyslb03, stateMPC_dslbaff03, stateMPC_llb03, stateMPC_dllbaff03);
stateMPC_LA_VSUB2_INDEXED_4(stateMPC_riub03, stateMPC_dzaff03, stateMPC_ubIdx03, stateMPC_dsubaff03);
stateMPC_LA_VSUB3_4(stateMPC_lubbysub03, stateMPC_dsubaff03, stateMPC_lub03, stateMPC_dlubaff03);
stateMPC_LA_VSUB_INDEXED_4(stateMPC_dzaff04, stateMPC_lbIdx04, stateMPC_rilb04, stateMPC_dslbaff04);
stateMPC_LA_VSUB3_4(stateMPC_llbbyslb04, stateMPC_dslbaff04, stateMPC_llb04, stateMPC_dllbaff04);
stateMPC_LA_VSUB2_INDEXED_4(stateMPC_riub04, stateMPC_dzaff04, stateMPC_ubIdx04, stateMPC_dsubaff04);
stateMPC_LA_VSUB3_4(stateMPC_lubbysub04, stateMPC_dsubaff04, stateMPC_lub04, stateMPC_dlubaff04);
stateMPC_LA_VSUB_INDEXED_4(stateMPC_dzaff05, stateMPC_lbIdx05, stateMPC_rilb05, stateMPC_dslbaff05);
stateMPC_LA_VSUB3_4(stateMPC_llbbyslb05, stateMPC_dslbaff05, stateMPC_llb05, stateMPC_dllbaff05);
stateMPC_LA_VSUB2_INDEXED_4(stateMPC_riub05, stateMPC_dzaff05, stateMPC_ubIdx05, stateMPC_dsubaff05);
stateMPC_LA_VSUB3_4(stateMPC_lubbysub05, stateMPC_dsubaff05, stateMPC_lub05, stateMPC_dlubaff05);
stateMPC_LA_VSUB_INDEXED_4(stateMPC_dzaff06, stateMPC_lbIdx06, stateMPC_rilb06, stateMPC_dslbaff06);
stateMPC_LA_VSUB3_4(stateMPC_llbbyslb06, stateMPC_dslbaff06, stateMPC_llb06, stateMPC_dllbaff06);
stateMPC_LA_VSUB2_INDEXED_4(stateMPC_riub06, stateMPC_dzaff06, stateMPC_ubIdx06, stateMPC_dsubaff06);
stateMPC_LA_VSUB3_4(stateMPC_lubbysub06, stateMPC_dsubaff06, stateMPC_lub06, stateMPC_dlubaff06);
stateMPC_LA_VSUB_INDEXED_4(stateMPC_dzaff07, stateMPC_lbIdx07, stateMPC_rilb07, stateMPC_dslbaff07);
stateMPC_LA_VSUB3_4(stateMPC_llbbyslb07, stateMPC_dslbaff07, stateMPC_llb07, stateMPC_dllbaff07);
stateMPC_LA_VSUB2_INDEXED_4(stateMPC_riub07, stateMPC_dzaff07, stateMPC_ubIdx07, stateMPC_dsubaff07);
stateMPC_LA_VSUB3_4(stateMPC_lubbysub07, stateMPC_dsubaff07, stateMPC_lub07, stateMPC_dlubaff07);
stateMPC_LA_VSUB_INDEXED_4(stateMPC_dzaff08, stateMPC_lbIdx08, stateMPC_rilb08, stateMPC_dslbaff08);
stateMPC_LA_VSUB3_4(stateMPC_llbbyslb08, stateMPC_dslbaff08, stateMPC_llb08, stateMPC_dllbaff08);
stateMPC_LA_VSUB2_INDEXED_4(stateMPC_riub08, stateMPC_dzaff08, stateMPC_ubIdx08, stateMPC_dsubaff08);
stateMPC_LA_VSUB3_4(stateMPC_lubbysub08, stateMPC_dsubaff08, stateMPC_lub08, stateMPC_dlubaff08);
stateMPC_LA_VSUB_INDEXED_4(stateMPC_dzaff09, stateMPC_lbIdx09, stateMPC_rilb09, stateMPC_dslbaff09);
stateMPC_LA_VSUB3_4(stateMPC_llbbyslb09, stateMPC_dslbaff09, stateMPC_llb09, stateMPC_dllbaff09);
stateMPC_LA_VSUB2_INDEXED_4(stateMPC_riub09, stateMPC_dzaff09, stateMPC_ubIdx09, stateMPC_dsubaff09);
stateMPC_LA_VSUB3_4(stateMPC_lubbysub09, stateMPC_dsubaff09, stateMPC_lub09, stateMPC_dlubaff09);
stateMPC_LA_VSUB_INDEXED_4(stateMPC_dzaff10, stateMPC_lbIdx10, stateMPC_rilb10, stateMPC_dslbaff10);
stateMPC_LA_VSUB3_4(stateMPC_llbbyslb10, stateMPC_dslbaff10, stateMPC_llb10, stateMPC_dllbaff10);
stateMPC_LA_VSUB2_INDEXED_4(stateMPC_riub10, stateMPC_dzaff10, stateMPC_ubIdx10, stateMPC_dsubaff10);
stateMPC_LA_VSUB3_4(stateMPC_lubbysub10, stateMPC_dsubaff10, stateMPC_lub10, stateMPC_dlubaff10);
stateMPC_LA_VSUB_INDEXED_4(stateMPC_dzaff11, stateMPC_lbIdx11, stateMPC_rilb11, stateMPC_dslbaff11);
stateMPC_LA_VSUB3_4(stateMPC_llbbyslb11, stateMPC_dslbaff11, stateMPC_llb11, stateMPC_dllbaff11);
stateMPC_LA_VSUB2_INDEXED_4(stateMPC_riub11, stateMPC_dzaff11, stateMPC_ubIdx11, stateMPC_dsubaff11);
stateMPC_LA_VSUB3_4(stateMPC_lubbysub11, stateMPC_dsubaff11, stateMPC_lub11, stateMPC_dlubaff11);
stateMPC_LA_VSUB_INDEXED_4(stateMPC_dzaff12, stateMPC_lbIdx12, stateMPC_rilb12, stateMPC_dslbaff12);
stateMPC_LA_VSUB3_4(stateMPC_llbbyslb12, stateMPC_dslbaff12, stateMPC_llb12, stateMPC_dllbaff12);
stateMPC_LA_VSUB2_INDEXED_4(stateMPC_riub12, stateMPC_dzaff12, stateMPC_ubIdx12, stateMPC_dsubaff12);
stateMPC_LA_VSUB3_4(stateMPC_lubbysub12, stateMPC_dsubaff12, stateMPC_lub12, stateMPC_dlubaff12);
stateMPC_LA_VSUB_INDEXED_4(stateMPC_dzaff13, stateMPC_lbIdx13, stateMPC_rilb13, stateMPC_dslbaff13);
stateMPC_LA_VSUB3_4(stateMPC_llbbyslb13, stateMPC_dslbaff13, stateMPC_llb13, stateMPC_dllbaff13);
stateMPC_LA_VSUB2_INDEXED_4(stateMPC_riub13, stateMPC_dzaff13, stateMPC_ubIdx13, stateMPC_dsubaff13);
stateMPC_LA_VSUB3_4(stateMPC_lubbysub13, stateMPC_dsubaff13, stateMPC_lub13, stateMPC_dlubaff13);
stateMPC_LA_VSUB_INDEXED_2(stateMPC_dzaff14, stateMPC_lbIdx14, stateMPC_rilb14, stateMPC_dslbaff14);
stateMPC_LA_VSUB3_2(stateMPC_llbbyslb14, stateMPC_dslbaff14, stateMPC_llb14, stateMPC_dllbaff14);
stateMPC_LA_VSUB2_INDEXED_2(stateMPC_riub14, stateMPC_dzaff14, stateMPC_ubIdx14, stateMPC_dsubaff14);
stateMPC_LA_VSUB3_2(stateMPC_lubbysub14, stateMPC_dsubaff14, stateMPC_lub14, stateMPC_dlubaff14);
info->lsit_aff = stateMPC_LINESEARCH_BACKTRACKING_AFFINE(stateMPC_l, stateMPC_s, stateMPC_dl_aff, stateMPC_ds_aff, &info->step_aff, &info->mu_aff);
if( info->lsit_aff == stateMPC_NOPROGRESS ){
exitcode = stateMPC_NOPROGRESS; break;
}
sigma_3rdroot = info->mu_aff / info->mu;
info->sigma = sigma_3rdroot*sigma_3rdroot*sigma_3rdroot;
musigma = info->mu * info->sigma;
stateMPC_LA_VSUB5_116(stateMPC_ds_aff, stateMPC_dl_aff, musigma, stateMPC_ccrhs);
stateMPC_LA_VSUB6_INDEXED_4_4_4(stateMPC_ccrhsub00, stateMPC_sub00, stateMPC_ubIdx00, stateMPC_ccrhsl00, stateMPC_slb00, stateMPC_lbIdx00, stateMPC_rd00);
stateMPC_LA_VSUB6_INDEXED_4_4_4(stateMPC_ccrhsub01, stateMPC_sub01, stateMPC_ubIdx01, stateMPC_ccrhsl01, stateMPC_slb01, stateMPC_lbIdx01, stateMPC_rd01);
stateMPC_LA_DIAG_FORWARDSUB_4(stateMPC_Phi00, stateMPC_rd00, stateMPC_Lbyrd00);
stateMPC_LA_DIAG_FORWARDSUB_4(stateMPC_Phi01, stateMPC_rd01, stateMPC_Lbyrd01);
stateMPC_LA_DIAGZERO_MVM_2(stateMPC_W00, stateMPC_Lbyrd00, stateMPC_beta00);
stateMPC_LA_DENSE_FORWARDSUB_2(stateMPC_Ld00, stateMPC_beta00, stateMPC_yy00);
stateMPC_LA_DENSE_DIAGZERO_2MVMADD_2_4_4(stateMPC_V00, stateMPC_Lbyrd00, stateMPC_W01, stateMPC_Lbyrd01, stateMPC_beta01);
stateMPC_LA_DENSE_MVMSUB1_2_2(stateMPC_Lsd01, stateMPC_yy00, stateMPC_beta01, stateMPC_bmy01);
stateMPC_LA_DENSE_FORWARDSUB_2(stateMPC_Ld01, stateMPC_bmy01, stateMPC_yy01);
stateMPC_LA_VSUB6_INDEXED_4_4_4(stateMPC_ccrhsub02, stateMPC_sub02, stateMPC_ubIdx02, stateMPC_ccrhsl02, stateMPC_slb02, stateMPC_lbIdx02, stateMPC_rd02);
stateMPC_LA_DIAG_FORWARDSUB_4(stateMPC_Phi02, stateMPC_rd02, stateMPC_Lbyrd02);
stateMPC_LA_DENSE_DIAGZERO_2MVMADD_2_4_4(stateMPC_V01, stateMPC_Lbyrd01, stateMPC_W02, stateMPC_Lbyrd02, stateMPC_beta02);
stateMPC_LA_DENSE_MVMSUB1_2_2(stateMPC_Lsd02, stateMPC_yy01, stateMPC_beta02, stateMPC_bmy02);
stateMPC_LA_DENSE_FORWARDSUB_2(stateMPC_Ld02, stateMPC_bmy02, stateMPC_yy02);
stateMPC_LA_VSUB6_INDEXED_4_4_4(stateMPC_ccrhsub03, stateMPC_sub03, stateMPC_ubIdx03, stateMPC_ccrhsl03, stateMPC_slb03, stateMPC_lbIdx03, stateMPC_rd03);
stateMPC_LA_DIAG_FORWARDSUB_4(stateMPC_Phi03, stateMPC_rd03, stateMPC_Lbyrd03);
stateMPC_LA_DENSE_DIAGZERO_2MVMADD_2_4_4(stateMPC_V02, stateMPC_Lbyrd02, stateMPC_W03, stateMPC_Lbyrd03, stateMPC_beta03);
stateMPC_LA_DENSE_MVMSUB1_2_2(stateMPC_Lsd03, stateMPC_yy02, stateMPC_beta03, stateMPC_bmy03);
stateMPC_LA_DENSE_FORWARDSUB_2(stateMPC_Ld03, stateMPC_bmy03, stateMPC_yy03);
stateMPC_LA_VSUB6_INDEXED_4_4_4(stateMPC_ccrhsub04, stateMPC_sub04, stateMPC_ubIdx04, stateMPC_ccrhsl04, stateMPC_slb04, stateMPC_lbIdx04, stateMPC_rd04);
stateMPC_LA_DIAG_FORWARDSUB_4(stateMPC_Phi04, stateMPC_rd04, stateMPC_Lbyrd04);
stateMPC_LA_DENSE_DIAGZERO_2MVMADD_2_4_4(stateMPC_V03, stateMPC_Lbyrd03, stateMPC_W04, stateMPC_Lbyrd04, stateMPC_beta04);
stateMPC_LA_DENSE_MVMSUB1_2_2(stateMPC_Lsd04, stateMPC_yy03, stateMPC_beta04, stateMPC_bmy04);
stateMPC_LA_DENSE_FORWARDSUB_2(stateMPC_Ld04, stateMPC_bmy04, stateMPC_yy04);
stateMPC_LA_VSUB6_INDEXED_4_4_4(stateMPC_ccrhsub05, stateMPC_sub05, stateMPC_ubIdx05, stateMPC_ccrhsl05, stateMPC_slb05, stateMPC_lbIdx05, stateMPC_rd05);
stateMPC_LA_DIAG_FORWARDSUB_4(stateMPC_Phi05, stateMPC_rd05, stateMPC_Lbyrd05);
stateMPC_LA_DENSE_DIAGZERO_2MVMADD_2_4_4(stateMPC_V04, stateMPC_Lbyrd04, stateMPC_W05, stateMPC_Lbyrd05, stateMPC_beta05);
stateMPC_LA_DENSE_MVMSUB1_2_2(stateMPC_Lsd05, stateMPC_yy04, stateMPC_beta05, stateMPC_bmy05);
stateMPC_LA_DENSE_FORWARDSUB_2(stateMPC_Ld05, stateMPC_bmy05, stateMPC_yy05);
stateMPC_LA_VSUB6_INDEXED_4_4_4(stateMPC_ccrhsub06, stateMPC_sub06, stateMPC_ubIdx06, stateMPC_ccrhsl06, stateMPC_slb06, stateMPC_lbIdx06, stateMPC_rd06);
stateMPC_LA_DIAG_FORWARDSUB_4(stateMPC_Phi06, stateMPC_rd06, stateMPC_Lbyrd06);
stateMPC_LA_DENSE_DIAGZERO_2MVMADD_2_4_4(stateMPC_V05, stateMPC_Lbyrd05, stateMPC_W06, stateMPC_Lbyrd06, stateMPC_beta06);
stateMPC_LA_DENSE_MVMSUB1_2_2(stateMPC_Lsd06, stateMPC_yy05, stateMPC_beta06, stateMPC_bmy06);
stateMPC_LA_DENSE_FORWARDSUB_2(stateMPC_Ld06, stateMPC_bmy06, stateMPC_yy06);
stateMPC_LA_VSUB6_INDEXED_4_4_4(stateMPC_ccrhsub07, stateMPC_sub07, stateMPC_ubIdx07, stateMPC_ccrhsl07, stateMPC_slb07, stateMPC_lbIdx07, stateMPC_rd07);
stateMPC_LA_DIAG_FORWARDSUB_4(stateMPC_Phi07, stateMPC_rd07, stateMPC_Lbyrd07);
stateMPC_LA_DENSE_DIAGZERO_2MVMADD_2_4_4(stateMPC_V06, stateMPC_Lbyrd06, stateMPC_W07, stateMPC_Lbyrd07, stateMPC_beta07);
stateMPC_LA_DENSE_MVMSUB1_2_2(stateMPC_Lsd07, stateMPC_yy06, stateMPC_beta07, stateMPC_bmy07);
stateMPC_LA_DENSE_FORWARDSUB_2(stateMPC_Ld07, stateMPC_bmy07, stateMPC_yy07);
stateMPC_LA_VSUB6_INDEXED_4_4_4(stateMPC_ccrhsub08, stateMPC_sub08, stateMPC_ubIdx08, stateMPC_ccrhsl08, stateMPC_slb08, stateMPC_lbIdx08, stateMPC_rd08);
stateMPC_LA_DIAG_FORWARDSUB_4(stateMPC_Phi08, stateMPC_rd08, stateMPC_Lbyrd08);
stateMPC_LA_DENSE_DIAGZERO_2MVMADD_2_4_4(stateMPC_V07, stateMPC_Lbyrd07, stateMPC_W08, stateMPC_Lbyrd08, stateMPC_beta08);
stateMPC_LA_DENSE_MVMSUB1_2_2(stateMPC_Lsd08, stateMPC_yy07, stateMPC_beta08, stateMPC_bmy08);
stateMPC_LA_DENSE_FORWARDSUB_2(stateMPC_Ld08, stateMPC_bmy08, stateMPC_yy08);
stateMPC_LA_VSUB6_INDEXED_4_4_4(stateMPC_ccrhsub09, stateMPC_sub09, stateMPC_ubIdx09, stateMPC_ccrhsl09, stateMPC_slb09, stateMPC_lbIdx09, stateMPC_rd09);
stateMPC_LA_DIAG_FORWARDSUB_4(stateMPC_Phi09, stateMPC_rd09, stateMPC_Lbyrd09);
stateMPC_LA_DENSE_DIAGZERO_2MVMADD_2_4_4(stateMPC_V08, stateMPC_Lbyrd08, stateMPC_W09, stateMPC_Lbyrd09, stateMPC_beta09);
stateMPC_LA_DENSE_MVMSUB1_2_2(stateMPC_Lsd09, stateMPC_yy08, stateMPC_beta09, stateMPC_bmy09);
stateMPC_LA_DENSE_FORWARDSUB_2(stateMPC_Ld09, stateMPC_bmy09, stateMPC_yy09);
stateMPC_LA_VSUB6_INDEXED_4_4_4(stateMPC_ccrhsub10, stateMPC_sub10, stateMPC_ubIdx10, stateMPC_ccrhsl10, stateMPC_slb10, stateMPC_lbIdx10, stateMPC_rd10);
stateMPC_LA_DIAG_FORWARDSUB_4(stateMPC_Phi10, stateMPC_rd10, stateMPC_Lbyrd10);
stateMPC_LA_DENSE_DIAGZERO_2MVMADD_2_4_4(stateMPC_V09, stateMPC_Lbyrd09, stateMPC_W10, stateMPC_Lbyrd10, stateMPC_beta10);
stateMPC_LA_DENSE_MVMSUB1_2_2(stateMPC_Lsd10, stateMPC_yy09, stateMPC_beta10, stateMPC_bmy10);
stateMPC_LA_DENSE_FORWARDSUB_2(stateMPC_Ld10, stateMPC_bmy10, stateMPC_yy10);
stateMPC_LA_VSUB6_INDEXED_4_4_4(stateMPC_ccrhsub11, stateMPC_sub11, stateMPC_ubIdx11, stateMPC_ccrhsl11, stateMPC_slb11, stateMPC_lbIdx11, stateMPC_rd11);
stateMPC_LA_DIAG_FORWARDSUB_4(stateMPC_Phi11, stateMPC_rd11, stateMPC_Lbyrd11);
stateMPC_LA_DENSE_DIAGZERO_2MVMADD_2_4_4(stateMPC_V10, stateMPC_Lbyrd10, stateMPC_W11, stateMPC_Lbyrd11, stateMPC_beta11);
stateMPC_LA_DENSE_MVMSUB1_2_2(stateMPC_Lsd11, stateMPC_yy10, stateMPC_beta11, stateMPC_bmy11);
stateMPC_LA_DENSE_FORWARDSUB_2(stateMPC_Ld11, stateMPC_bmy11, stateMPC_yy11);
stateMPC_LA_VSUB6_INDEXED_4_4_4(stateMPC_ccrhsub12, stateMPC_sub12, stateMPC_ubIdx12, stateMPC_ccrhsl12, stateMPC_slb12, stateMPC_lbIdx12, stateMPC_rd12);
stateMPC_LA_DIAG_FORWARDSUB_4(stateMPC_Phi12, stateMPC_rd12, stateMPC_Lbyrd12);
stateMPC_LA_DENSE_DIAGZERO_2MVMADD_2_4_4(stateMPC_V11, stateMPC_Lbyrd11, stateMPC_W12, stateMPC_Lbyrd12, stateMPC_beta12);
stateMPC_LA_DENSE_MVMSUB1_2_2(stateMPC_Lsd12, stateMPC_yy11, stateMPC_beta12, stateMPC_bmy12);
stateMPC_LA_DENSE_FORWARDSUB_2(stateMPC_Ld12, stateMPC_bmy12, stateMPC_yy12);
stateMPC_LA_VSUB6_INDEXED_4_4_4(stateMPC_ccrhsub13, stateMPC_sub13, stateMPC_ubIdx13, stateMPC_ccrhsl13, stateMPC_slb13, stateMPC_lbIdx13, stateMPC_rd13);
stateMPC_LA_DIAG_FORWARDSUB_4(stateMPC_Phi13, stateMPC_rd13, stateMPC_Lbyrd13);
stateMPC_LA_DENSE_DIAGZERO_2MVMADD_2_4_4(stateMPC_V12, stateMPC_Lbyrd12, stateMPC_W13, stateMPC_Lbyrd13, stateMPC_beta13);
stateMPC_LA_DENSE_MVMSUB1_2_2(stateMPC_Lsd13, stateMPC_yy12, stateMPC_beta13, stateMPC_bmy13);
stateMPC_LA_DENSE_FORWARDSUB_2(stateMPC_Ld13, stateMPC_bmy13, stateMPC_yy13);
stateMPC_LA_VSUB6_INDEXED_2_2_2(stateMPC_ccrhsub14, stateMPC_sub14, stateMPC_ubIdx14, stateMPC_ccrhsl14, stateMPC_slb14, stateMPC_lbIdx14, stateMPC_rd14);
stateMPC_LA_DIAG_FORWARDSUB_2(stateMPC_Phi14, stateMPC_rd14, stateMPC_Lbyrd14);
stateMPC_LA_DENSE_DIAGZERO_2MVMADD_2_4_2(stateMPC_V13, stateMPC_Lbyrd13, stateMPC_W14, stateMPC_Lbyrd14, stateMPC_beta14);
stateMPC_LA_DENSE_MVMSUB1_2_2(stateMPC_Lsd14, stateMPC_yy13, stateMPC_beta14, stateMPC_bmy14);
stateMPC_LA_DENSE_FORWARDSUB_2(stateMPC_Ld14, stateMPC_bmy14, stateMPC_yy14);
stateMPC_LA_DENSE_BACKWARDSUB_2(stateMPC_Ld14, stateMPC_yy14, stateMPC_dvcc14);
stateMPC_LA_DENSE_MTVMSUB_2_2(stateMPC_Lsd14, stateMPC_dvcc14, stateMPC_yy13, stateMPC_bmy13);
stateMPC_LA_DENSE_BACKWARDSUB_2(stateMPC_Ld13, stateMPC_bmy13, stateMPC_dvcc13);
stateMPC_LA_DENSE_MTVMSUB_2_2(stateMPC_Lsd13, stateMPC_dvcc13, stateMPC_yy12, stateMPC_bmy12);
stateMPC_LA_DENSE_BACKWARDSUB_2(stateMPC_Ld12, stateMPC_bmy12, stateMPC_dvcc12);
stateMPC_LA_DENSE_MTVMSUB_2_2(stateMPC_Lsd12, stateMPC_dvcc12, stateMPC_yy11, stateMPC_bmy11);
stateMPC_LA_DENSE_BACKWARDSUB_2(stateMPC_Ld11, stateMPC_bmy11, stateMPC_dvcc11);
stateMPC_LA_DENSE_MTVMSUB_2_2(stateMPC_Lsd11, stateMPC_dvcc11, stateMPC_yy10, stateMPC_bmy10);
stateMPC_LA_DENSE_BACKWARDSUB_2(stateMPC_Ld10, stateMPC_bmy10, stateMPC_dvcc10);
stateMPC_LA_DENSE_MTVMSUB_2_2(stateMPC_Lsd10, stateMPC_dvcc10, stateMPC_yy09, stateMPC_bmy09);
stateMPC_LA_DENSE_BACKWARDSUB_2(stateMPC_Ld09, stateMPC_bmy09, stateMPC_dvcc09);
stateMPC_LA_DENSE_MTVMSUB_2_2(stateMPC_Lsd09, stateMPC_dvcc09, stateMPC_yy08, stateMPC_bmy08);
stateMPC_LA_DENSE_BACKWARDSUB_2(stateMPC_Ld08, stateMPC_bmy08, stateMPC_dvcc08);
stateMPC_LA_DENSE_MTVMSUB_2_2(stateMPC_Lsd08, stateMPC_dvcc08, stateMPC_yy07, stateMPC_bmy07);
stateMPC_LA_DENSE_BACKWARDSUB_2(stateMPC_Ld07, stateMPC_bmy07, stateMPC_dvcc07);
stateMPC_LA_DENSE_MTVMSUB_2_2(stateMPC_Lsd07, stateMPC_dvcc07, stateMPC_yy06, stateMPC_bmy06);
stateMPC_LA_DENSE_BACKWARDSUB_2(stateMPC_Ld06, stateMPC_bmy06, stateMPC_dvcc06);
stateMPC_LA_DENSE_MTVMSUB_2_2(stateMPC_Lsd06, stateMPC_dvcc06, stateMPC_yy05, stateMPC_bmy05);
stateMPC_LA_DENSE_BACKWARDSUB_2(stateMPC_Ld05, stateMPC_bmy05, stateMPC_dvcc05);
stateMPC_LA_DENSE_MTVMSUB_2_2(stateMPC_Lsd05, stateMPC_dvcc05, stateMPC_yy04, stateMPC_bmy04);
stateMPC_LA_DENSE_BACKWARDSUB_2(stateMPC_Ld04, stateMPC_bmy04, stateMPC_dvcc04);
stateMPC_LA_DENSE_MTVMSUB_2_2(stateMPC_Lsd04, stateMPC_dvcc04, stateMPC_yy03, stateMPC_bmy03);
stateMPC_LA_DENSE_BACKWARDSUB_2(stateMPC_Ld03, stateMPC_bmy03, stateMPC_dvcc03);
stateMPC_LA_DENSE_MTVMSUB_2_2(stateMPC_Lsd03, stateMPC_dvcc03, stateMPC_yy02, stateMPC_bmy02);
stateMPC_LA_DENSE_BACKWARDSUB_2(stateMPC_Ld02, stateMPC_bmy02, stateMPC_dvcc02);
stateMPC_LA_DENSE_MTVMSUB_2_2(stateMPC_Lsd02, stateMPC_dvcc02, stateMPC_yy01, stateMPC_bmy01);
stateMPC_LA_DENSE_BACKWARDSUB_2(stateMPC_Ld01, stateMPC_bmy01, stateMPC_dvcc01);
stateMPC_LA_DENSE_MTVMSUB_2_2(stateMPC_Lsd01, stateMPC_dvcc01, stateMPC_yy00, stateMPC_bmy00);
stateMPC_LA_DENSE_BACKWARDSUB_2(stateMPC_Ld00, stateMPC_bmy00, stateMPC_dvcc00);
stateMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C1, stateMPC_dvcc01, stateMPC_D00, stateMPC_dvcc00, stateMPC_grad_eq00);
stateMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C2, stateMPC_dvcc02, stateMPC_D01, stateMPC_dvcc01, stateMPC_grad_eq01);
stateMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C3, stateMPC_dvcc03, stateMPC_D01, stateMPC_dvcc02, stateMPC_grad_eq02);
stateMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C4, stateMPC_dvcc04, stateMPC_D01, stateMPC_dvcc03, stateMPC_grad_eq03);
stateMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C5, stateMPC_dvcc05, stateMPC_D01, stateMPC_dvcc04, stateMPC_grad_eq04);
stateMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C6, stateMPC_dvcc06, stateMPC_D01, stateMPC_dvcc05, stateMPC_grad_eq05);
stateMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C7, stateMPC_dvcc07, stateMPC_D01, stateMPC_dvcc06, stateMPC_grad_eq06);
stateMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C8, stateMPC_dvcc08, stateMPC_D01, stateMPC_dvcc07, stateMPC_grad_eq07);
stateMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C9, stateMPC_dvcc09, stateMPC_D01, stateMPC_dvcc08, stateMPC_grad_eq08);
stateMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C10, stateMPC_dvcc10, stateMPC_D01, stateMPC_dvcc09, stateMPC_grad_eq09);
stateMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C11, stateMPC_dvcc11, stateMPC_D01, stateMPC_dvcc10, stateMPC_grad_eq10);
stateMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C12, stateMPC_dvcc12, stateMPC_D01, stateMPC_dvcc11, stateMPC_grad_eq11);
stateMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C13, stateMPC_dvcc13, stateMPC_D01, stateMPC_dvcc12, stateMPC_grad_eq12);
stateMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C14, stateMPC_dvcc14, stateMPC_D01, stateMPC_dvcc13, stateMPC_grad_eq13);
stateMPC_LA_DIAGZERO_MTVM_2_2(stateMPC_D14, stateMPC_dvcc14, stateMPC_grad_eq14);
stateMPC_LA_VSUB_58(stateMPC_rd, stateMPC_grad_eq, stateMPC_rd);
stateMPC_LA_DIAG_FORWARDBACKWARDSUB_4(stateMPC_Phi00, stateMPC_rd00, stateMPC_dzcc00);
stateMPC_LA_DIAG_FORWARDBACKWARDSUB_4(stateMPC_Phi01, stateMPC_rd01, stateMPC_dzcc01);
stateMPC_LA_DIAG_FORWARDBACKWARDSUB_4(stateMPC_Phi02, stateMPC_rd02, stateMPC_dzcc02);
stateMPC_LA_DIAG_FORWARDBACKWARDSUB_4(stateMPC_Phi03, stateMPC_rd03, stateMPC_dzcc03);
stateMPC_LA_DIAG_FORWARDBACKWARDSUB_4(stateMPC_Phi04, stateMPC_rd04, stateMPC_dzcc04);
stateMPC_LA_DIAG_FORWARDBACKWARDSUB_4(stateMPC_Phi05, stateMPC_rd05, stateMPC_dzcc05);
stateMPC_LA_DIAG_FORWARDBACKWARDSUB_4(stateMPC_Phi06, stateMPC_rd06, stateMPC_dzcc06);
stateMPC_LA_DIAG_FORWARDBACKWARDSUB_4(stateMPC_Phi07, stateMPC_rd07, stateMPC_dzcc07);
stateMPC_LA_DIAG_FORWARDBACKWARDSUB_4(stateMPC_Phi08, stateMPC_rd08, stateMPC_dzcc08);
stateMPC_LA_DIAG_FORWARDBACKWARDSUB_4(stateMPC_Phi09, stateMPC_rd09, stateMPC_dzcc09);
stateMPC_LA_DIAG_FORWARDBACKWARDSUB_4(stateMPC_Phi10, stateMPC_rd10, stateMPC_dzcc10);
stateMPC_LA_DIAG_FORWARDBACKWARDSUB_4(stateMPC_Phi11, stateMPC_rd11, stateMPC_dzcc11);
stateMPC_LA_DIAG_FORWARDBACKWARDSUB_4(stateMPC_Phi12, stateMPC_rd12, stateMPC_dzcc12);
stateMPC_LA_DIAG_FORWARDBACKWARDSUB_4(stateMPC_Phi13, stateMPC_rd13, stateMPC_dzcc13);
stateMPC_LA_DIAG_FORWARDBACKWARDSUB_2(stateMPC_Phi14, stateMPC_rd14, stateMPC_dzcc14);
stateMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_4(stateMPC_ccrhsl00, stateMPC_slb00, stateMPC_llbbyslb00, stateMPC_dzcc00, stateMPC_lbIdx00, stateMPC_dllbcc00);
stateMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_4(stateMPC_ccrhsub00, stateMPC_sub00, stateMPC_lubbysub00, stateMPC_dzcc00, stateMPC_ubIdx00, stateMPC_dlubcc00);
stateMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_4(stateMPC_ccrhsl01, stateMPC_slb01, stateMPC_llbbyslb01, stateMPC_dzcc01, stateMPC_lbIdx01, stateMPC_dllbcc01);
stateMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_4(stateMPC_ccrhsub01, stateMPC_sub01, stateMPC_lubbysub01, stateMPC_dzcc01, stateMPC_ubIdx01, stateMPC_dlubcc01);
stateMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_4(stateMPC_ccrhsl02, stateMPC_slb02, stateMPC_llbbyslb02, stateMPC_dzcc02, stateMPC_lbIdx02, stateMPC_dllbcc02);
stateMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_4(stateMPC_ccrhsub02, stateMPC_sub02, stateMPC_lubbysub02, stateMPC_dzcc02, stateMPC_ubIdx02, stateMPC_dlubcc02);
stateMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_4(stateMPC_ccrhsl03, stateMPC_slb03, stateMPC_llbbyslb03, stateMPC_dzcc03, stateMPC_lbIdx03, stateMPC_dllbcc03);
stateMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_4(stateMPC_ccrhsub03, stateMPC_sub03, stateMPC_lubbysub03, stateMPC_dzcc03, stateMPC_ubIdx03, stateMPC_dlubcc03);
stateMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_4(stateMPC_ccrhsl04, stateMPC_slb04, stateMPC_llbbyslb04, stateMPC_dzcc04, stateMPC_lbIdx04, stateMPC_dllbcc04);
stateMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_4(stateMPC_ccrhsub04, stateMPC_sub04, stateMPC_lubbysub04, stateMPC_dzcc04, stateMPC_ubIdx04, stateMPC_dlubcc04);
stateMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_4(stateMPC_ccrhsl05, stateMPC_slb05, stateMPC_llbbyslb05, stateMPC_dzcc05, stateMPC_lbIdx05, stateMPC_dllbcc05);
stateMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_4(stateMPC_ccrhsub05, stateMPC_sub05, stateMPC_lubbysub05, stateMPC_dzcc05, stateMPC_ubIdx05, stateMPC_dlubcc05);
stateMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_4(stateMPC_ccrhsl06, stateMPC_slb06, stateMPC_llbbyslb06, stateMPC_dzcc06, stateMPC_lbIdx06, stateMPC_dllbcc06);
stateMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_4(stateMPC_ccrhsub06, stateMPC_sub06, stateMPC_lubbysub06, stateMPC_dzcc06, stateMPC_ubIdx06, stateMPC_dlubcc06);
stateMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_4(stateMPC_ccrhsl07, stateMPC_slb07, stateMPC_llbbyslb07, stateMPC_dzcc07, stateMPC_lbIdx07, stateMPC_dllbcc07);
stateMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_4(stateMPC_ccrhsub07, stateMPC_sub07, stateMPC_lubbysub07, stateMPC_dzcc07, stateMPC_ubIdx07, stateMPC_dlubcc07);
stateMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_4(stateMPC_ccrhsl08, stateMPC_slb08, stateMPC_llbbyslb08, stateMPC_dzcc08, stateMPC_lbIdx08, stateMPC_dllbcc08);
stateMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_4(stateMPC_ccrhsub08, stateMPC_sub08, stateMPC_lubbysub08, stateMPC_dzcc08, stateMPC_ubIdx08, stateMPC_dlubcc08);
stateMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_4(stateMPC_ccrhsl09, stateMPC_slb09, stateMPC_llbbyslb09, stateMPC_dzcc09, stateMPC_lbIdx09, stateMPC_dllbcc09);
stateMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_4(stateMPC_ccrhsub09, stateMPC_sub09, stateMPC_lubbysub09, stateMPC_dzcc09, stateMPC_ubIdx09, stateMPC_dlubcc09);
stateMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_4(stateMPC_ccrhsl10, stateMPC_slb10, stateMPC_llbbyslb10, stateMPC_dzcc10, stateMPC_lbIdx10, stateMPC_dllbcc10);
stateMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_4(stateMPC_ccrhsub10, stateMPC_sub10, stateMPC_lubbysub10, stateMPC_dzcc10, stateMPC_ubIdx10, stateMPC_dlubcc10);
stateMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_4(stateMPC_ccrhsl11, stateMPC_slb11, stateMPC_llbbyslb11, stateMPC_dzcc11, stateMPC_lbIdx11, stateMPC_dllbcc11);
stateMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_4(stateMPC_ccrhsub11, stateMPC_sub11, stateMPC_lubbysub11, stateMPC_dzcc11, stateMPC_ubIdx11, stateMPC_dlubcc11);
stateMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_4(stateMPC_ccrhsl12, stateMPC_slb12, stateMPC_llbbyslb12, stateMPC_dzcc12, stateMPC_lbIdx12, stateMPC_dllbcc12);
stateMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_4(stateMPC_ccrhsub12, stateMPC_sub12, stateMPC_lubbysub12, stateMPC_dzcc12, stateMPC_ubIdx12, stateMPC_dlubcc12);
stateMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_4(stateMPC_ccrhsl13, stateMPC_slb13, stateMPC_llbbyslb13, stateMPC_dzcc13, stateMPC_lbIdx13, stateMPC_dllbcc13);
stateMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_4(stateMPC_ccrhsub13, stateMPC_sub13, stateMPC_lubbysub13, stateMPC_dzcc13, stateMPC_ubIdx13, stateMPC_dlubcc13);
stateMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_2(stateMPC_ccrhsl14, stateMPC_slb14, stateMPC_llbbyslb14, stateMPC_dzcc14, stateMPC_lbIdx14, stateMPC_dllbcc14);
stateMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_2(stateMPC_ccrhsub14, stateMPC_sub14, stateMPC_lubbysub14, stateMPC_dzcc14, stateMPC_ubIdx14, stateMPC_dlubcc14);
stateMPC_LA_VSUB7_116(stateMPC_l, stateMPC_ccrhs, stateMPC_s, stateMPC_dl_cc, stateMPC_ds_cc);
stateMPC_LA_VADD_58(stateMPC_dz_cc, stateMPC_dz_aff);
stateMPC_LA_VADD_30(stateMPC_dv_cc, stateMPC_dv_aff);
stateMPC_LA_VADD_116(stateMPC_dl_cc, stateMPC_dl_aff);
stateMPC_LA_VADD_116(stateMPC_ds_cc, stateMPC_ds_aff);
info->lsit_cc = stateMPC_LINESEARCH_BACKTRACKING_COMBINED(stateMPC_z, stateMPC_v, stateMPC_l, stateMPC_s, stateMPC_dz_cc, stateMPC_dv_cc, stateMPC_dl_cc, stateMPC_ds_cc, &info->step_cc, &info->mu);
if( info->lsit_cc == stateMPC_NOPROGRESS ){
exitcode = stateMPC_NOPROGRESS; break;
}
info->it++;
}
output->z1[0] = stateMPC_z00[0];
output->z1[1] = stateMPC_z00[1];
output->z1[2] = stateMPC_z00[2];
output->z1[3] = stateMPC_z00[3];
output->z2[0] = stateMPC_z01[0];
output->z2[1] = stateMPC_z01[1];
output->z2[2] = stateMPC_z01[2];
output->z2[3] = stateMPC_z01[3];
output->z3[0] = stateMPC_z02[0];
output->z3[1] = stateMPC_z02[1];
output->z3[2] = stateMPC_z02[2];
output->z3[3] = stateMPC_z02[3];
output->z4[0] = stateMPC_z03[0];
output->z4[1] = stateMPC_z03[1];
output->z4[2] = stateMPC_z03[2];
output->z4[3] = stateMPC_z03[3];
output->z5[0] = stateMPC_z04[0];
output->z5[1] = stateMPC_z04[1];
output->z5[2] = stateMPC_z04[2];
output->z5[3] = stateMPC_z04[3];
output->z6[0] = stateMPC_z05[0];
output->z6[1] = stateMPC_z05[1];
output->z6[2] = stateMPC_z05[2];
output->z6[3] = stateMPC_z05[3];
output->z7[0] = stateMPC_z06[0];
output->z7[1] = stateMPC_z06[1];
output->z7[2] = stateMPC_z06[2];
output->z7[3] = stateMPC_z06[3];
output->z8[0] = stateMPC_z07[0];
output->z8[1] = stateMPC_z07[1];
output->z8[2] = stateMPC_z07[2];
output->z8[3] = stateMPC_z07[3];
output->z9[0] = stateMPC_z08[0];
output->z9[1] = stateMPC_z08[1];
output->z9[2] = stateMPC_z08[2];
output->z9[3] = stateMPC_z08[3];
output->z10[0] = stateMPC_z09[0];
output->z10[1] = stateMPC_z09[1];
output->z10[2] = stateMPC_z09[2];
output->z10[3] = stateMPC_z09[3];
output->z11[0] = stateMPC_z10[0];
output->z11[1] = stateMPC_z10[1];
output->z11[2] = stateMPC_z10[2];
output->z11[3] = stateMPC_z10[3];
output->z12[0] = stateMPC_z11[0];
output->z12[1] = stateMPC_z11[1];
output->z12[2] = stateMPC_z11[2];
output->z12[3] = stateMPC_z11[3];
output->z13[0] = stateMPC_z12[0];
output->z13[1] = stateMPC_z12[1];
output->z13[2] = stateMPC_z12[2];
output->z13[3] = stateMPC_z12[3];
output->z14[0] = stateMPC_z13[0];
output->z14[1] = stateMPC_z13[1];
output->z14[2] = stateMPC_z13[2];
output->z14[3] = stateMPC_z13[3];
output->z15[0] = stateMPC_z14[0];
output->z15[1] = stateMPC_z14[1];

#if stateMPC_SET_TIMING == 1
info->solvetime = stateMPC_toc(&solvertimer);
#if stateMPC_SET_PRINTLEVEL > 0 && stateMPC_SET_TIMING == 1
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
