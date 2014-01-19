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

#include "beliefMPC.h"

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
 * Initializes a vector of length 68 with a value.
 */
void beliefMPC_LA_INITIALIZEVECTOR_68(beliefMPC_FLOAT* vec, beliefMPC_FLOAT value)
{
	int i;
	for( i=0; i<68; i++ )
	{
		vec[i] = value;
	}
}


/*
 * Initializes a vector of length 50 with a value.
 */
void beliefMPC_LA_INITIALIZEVECTOR_50(beliefMPC_FLOAT* vec, beliefMPC_FLOAT value)
{
	int i;
	for( i=0; i<50; i++ )
	{
		vec[i] = value;
	}
}


/*
 * Initializes a vector of length 136 with a value.
 */
void beliefMPC_LA_INITIALIZEVECTOR_136(beliefMPC_FLOAT* vec, beliefMPC_FLOAT value)
{
	int i;
	for( i=0; i<136; i++ )
	{
		vec[i] = value;
	}
}


/* 
 * Calculates a dot product and adds it to a variable: z += x'*y; 
 * This function is for vectors of length 136.
 */
void beliefMPC_LA_DOTACC_136(beliefMPC_FLOAT *x, beliefMPC_FLOAT *y, beliefMPC_FLOAT *z)
{
	int i;
	for( i=0; i<136; i++ ){
		*z += x[i]*y[i];
	}
}


/*
 * Calculates the gradient and the value for a quadratic function 0.5*z'*H*z + f'*z
 *
 * INPUTS:     H  - Symmetric Hessian, diag matrix of size [7 x 7]
 *             f  - column vector of size 7
 *             z  - column vector of size 7
 *
 * OUTPUTS: grad  - gradient at z (= H*z + f), column vector of size 7
 *          value <-- value + 0.5*z'*H*z + f'*z (value will be modified)
 */
void beliefMPC_LA_DIAG_QUADFCN_7(beliefMPC_FLOAT* H, beliefMPC_FLOAT* f, beliefMPC_FLOAT* z, beliefMPC_FLOAT* grad, beliefMPC_FLOAT* value)
{
	int i;
	beliefMPC_FLOAT hz;	
	for( i=0; i<7; i++){
		hz = H[i]*z[i];
		grad[i] = hz + f[i];
		*value += 0.5*hz*z[i] + f[i]*z[i];
	}
}


/*
 * Calculates the gradient and the value for a quadratic function 0.5*z'*H*z + f'*z
 *
 * INPUTS:     H  - Symmetric Hessian, diag matrix of size [5 x 5]
 *             f  - column vector of size 5
 *             z  - column vector of size 5
 *
 * OUTPUTS: grad  - gradient at z (= H*z + f), column vector of size 5
 *          value <-- value + 0.5*z'*H*z + f'*z (value will be modified)
 */
void beliefMPC_LA_DIAG_QUADFCN_5(beliefMPC_FLOAT* H, beliefMPC_FLOAT* f, beliefMPC_FLOAT* z, beliefMPC_FLOAT* grad, beliefMPC_FLOAT* value)
{
	int i;
	beliefMPC_FLOAT hz;	
	for( i=0; i<5; i++){
		hz = H[i]*z[i];
		grad[i] = hz + f[i];
		*value += 0.5*hz*z[i] + f[i]*z[i];
	}
}


/* 
 * Computes r = A*x + B*u - b
 * and      y = max([norm(r,inf), y])
 * and      z -= l'*r
 * where A is stored in column major format
 */
void beliefMPC_LA_DENSE_MVMSUB3_10_7_7(beliefMPC_FLOAT *A, beliefMPC_FLOAT *x, beliefMPC_FLOAT *B, beliefMPC_FLOAT *u, beliefMPC_FLOAT *b, beliefMPC_FLOAT *l, beliefMPC_FLOAT *r, beliefMPC_FLOAT *z, beliefMPC_FLOAT *y)
{
	int i;
	int j;
	int k = 0;
	int m = 0;
	int n;
	beliefMPC_FLOAT AxBu[10];
	beliefMPC_FLOAT norm = *y;
	beliefMPC_FLOAT lr = 0;

	/* do A*x + B*u first */
	for( i=0; i<10; i++ ){
		AxBu[i] = A[k++]*x[0] + B[m++]*u[0];
	}	
	for( j=1; j<7; j++ ){		
		for( i=0; i<10; i++ ){
			AxBu[i] += A[k++]*x[j];
		}
	}
	
	for( n=1; n<7; n++ ){
		for( i=0; i<10; i++ ){
			AxBu[i] += B[m++]*u[n];
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
void beliefMPC_LA_DENSE_DIAGZERO_MVMSUB3_5_7_7(beliefMPC_FLOAT *A, beliefMPC_FLOAT *x, beliefMPC_FLOAT *B, beliefMPC_FLOAT *u, beliefMPC_FLOAT *b, beliefMPC_FLOAT *l, beliefMPC_FLOAT *r, beliefMPC_FLOAT *z, beliefMPC_FLOAT *y)
{
	int i;
	int j;
	int k = 0;
	beliefMPC_FLOAT AxBu[5];
	beliefMPC_FLOAT norm = *y;
	beliefMPC_FLOAT lr = 0;

	/* do A*x + B*u first */
	for( i=0; i<5; i++ ){
		AxBu[i] = A[k++]*x[0] + B[i]*u[i];
	}	

	for( j=1; j<7; j++ ){		
		for( i=0; i<5; i++ ){
			AxBu[i] += A[k++]*x[j];
		}
	}

	for( i=0; i<5; i++ ){
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
void beliefMPC_LA_DENSE_DIAGZERO_MVMSUB3_5_7_5(beliefMPC_FLOAT *A, beliefMPC_FLOAT *x, beliefMPC_FLOAT *B, beliefMPC_FLOAT *u, beliefMPC_FLOAT *b, beliefMPC_FLOAT *l, beliefMPC_FLOAT *r, beliefMPC_FLOAT *z, beliefMPC_FLOAT *y)
{
	int i;
	int j;
	int k = 0;
	beliefMPC_FLOAT AxBu[5];
	beliefMPC_FLOAT norm = *y;
	beliefMPC_FLOAT lr = 0;

	/* do A*x + B*u first */
	for( i=0; i<5; i++ ){
		AxBu[i] = A[k++]*x[0] + B[i]*u[i];
	}	

	for( j=1; j<7; j++ ){		
		for( i=0; i<5; i++ ){
			AxBu[i] += A[k++]*x[j];
		}
	}

	for( i=0; i<5; i++ ){
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
 * Matrix vector multiplication y = M'*x where M is of size [10 x 7]
 * and stored in column major format. Note the transpose of M!
 */
void beliefMPC_LA_DENSE_MTVM_10_7(beliefMPC_FLOAT *M, beliefMPC_FLOAT *x, beliefMPC_FLOAT *y)
{
	int i;
	int j;
	int k = 0; 
	for( i=0; i<7; i++ ){
		y[i] = 0;
		for( j=0; j<10; j++ ){
			y[i] += M[k++]*x[j];
		}
	}
}


/*
 * Matrix vector multiplication z = A'*x + B'*y 
 * where A is of size [5 x 7]
 * and B is of size [10 x 7]
 * and stored in column major format. Note the transposes of A and B!
 */
void beliefMPC_LA_DENSE_MTVM2_5_7_10(beliefMPC_FLOAT *A, beliefMPC_FLOAT *x, beliefMPC_FLOAT *B, beliefMPC_FLOAT *y, beliefMPC_FLOAT *z)
{
	int i;
	int j;
	int k = 0;
	int n;
	int m = 0;
	for( i=0; i<7; i++ ){
		z[i] = 0;
		for( j=0; j<5; j++ ){
			z[i] += A[k++]*x[j];
		}
		for( n=0; n<10; n++ ){
			z[i] += B[m++]*y[n];
		}
	}
}


/*
 * Matrix vector multiplication z = A'*x + B'*y 
 * where A is of size [5 x 7] and stored in column major format.
 * and B is of size [5 x 7] and stored in diagzero format
 * Note the transposes of A and B!
 */
void beliefMPC_LA_DENSE_DIAGZERO_MTVM2_5_7_5(beliefMPC_FLOAT *A, beliefMPC_FLOAT *x, beliefMPC_FLOAT *B, beliefMPC_FLOAT *y, beliefMPC_FLOAT *z)
{
	int i;
	int j;
	int k = 0;
	for( i=0; i<5; i++ ){
		z[i] = 0;
		for( j=0; j<5; j++ ){
			z[i] += A[k++]*x[j];
		}
		z[i] += B[i]*y[i];
	}
	for( i=5 ;i<7; i++ ){
		z[i] = 0;
		for( j=0; j<5; j++ ){
			z[i] += A[k++]*x[j];
		}
	}
}


/*
 * Matrix vector multiplication y = M'*x where M is of size [5 x 5]
 * and stored in diagzero format. Note the transpose of M!
 */
void beliefMPC_LA_DIAGZERO_MTVM_5_5(beliefMPC_FLOAT *M, beliefMPC_FLOAT *x, beliefMPC_FLOAT *y)
{
	int i;
	for( i=0; i<5; i++ ){
		y[i] = M[i]*x[i];
	}
}


/*
 * Vector subtraction and addition.
 *	 Input: five vectors t, tidx, u, v, w and two scalars z and r
 *	 Output: y = t(tidx) - u + w
 *           z = z - v'*x;
 *           r = max([norm(y,inf), z]);
 * for vectors of length 7. Output z is of course scalar.
 */
void beliefMPC_LA_VSUBADD3_7(beliefMPC_FLOAT* t, beliefMPC_FLOAT* u, int* uidx, beliefMPC_FLOAT* v, beliefMPC_FLOAT* w, beliefMPC_FLOAT* y, beliefMPC_FLOAT* z, beliefMPC_FLOAT* r)
{
	int i;
	beliefMPC_FLOAT norm = *r;
	beliefMPC_FLOAT vx = 0;
	beliefMPC_FLOAT x;
	for( i=0; i<7; i++){
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
 * for vectors of length 7. Output z is of course scalar.
 */
void beliefMPC_LA_VSUBADD2_7(beliefMPC_FLOAT* t, int* tidx, beliefMPC_FLOAT* u, beliefMPC_FLOAT* v, beliefMPC_FLOAT* w, beliefMPC_FLOAT* y, beliefMPC_FLOAT* z, beliefMPC_FLOAT* r)
{
	int i;
	beliefMPC_FLOAT norm = *r;
	beliefMPC_FLOAT vx = 0;
	beliefMPC_FLOAT x;
	for( i=0; i<7; i++){
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
 * for vectors of length 5. Output z is of course scalar.
 */
void beliefMPC_LA_VSUBADD3_5(beliefMPC_FLOAT* t, beliefMPC_FLOAT* u, int* uidx, beliefMPC_FLOAT* v, beliefMPC_FLOAT* w, beliefMPC_FLOAT* y, beliefMPC_FLOAT* z, beliefMPC_FLOAT* r)
{
	int i;
	beliefMPC_FLOAT norm = *r;
	beliefMPC_FLOAT vx = 0;
	beliefMPC_FLOAT x;
	for( i=0; i<5; i++){
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
void beliefMPC_LA_VSUBADD2_5(beliefMPC_FLOAT* t, int* tidx, beliefMPC_FLOAT* u, beliefMPC_FLOAT* v, beliefMPC_FLOAT* w, beliefMPC_FLOAT* y, beliefMPC_FLOAT* z, beliefMPC_FLOAT* r)
{
	int i;
	beliefMPC_FLOAT norm = *r;
	beliefMPC_FLOAT vx = 0;
	beliefMPC_FLOAT x;
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
 * Computes inequality constraints gradient-
 * Special function for box constraints of length 7
 * Returns also L/S, a value that is often used elsewhere.
 */
void beliefMPC_LA_INEQ_B_GRAD_7_7_7(beliefMPC_FLOAT *lu, beliefMPC_FLOAT *su, beliefMPC_FLOAT *ru, beliefMPC_FLOAT *ll, beliefMPC_FLOAT *sl, beliefMPC_FLOAT *rl, int* lbIdx, int* ubIdx, beliefMPC_FLOAT *grad, beliefMPC_FLOAT *lubysu, beliefMPC_FLOAT *llbysl)
{
	int i;
	for( i=0; i<7; i++ ){
		grad[i] = 0;
	}
	for( i=0; i<7; i++ ){		
		llbysl[i] = ll[i] / sl[i];
		grad[lbIdx[i]] -= llbysl[i]*rl[i];
	}
	for( i=0; i<7; i++ ){
		lubysu[i] = lu[i] / su[i];
		grad[ubIdx[i]] += lubysu[i]*ru[i];
	}
}


/*
 * Computes inequality constraints gradient-
 * Special function for box constraints of length 5
 * Returns also L/S, a value that is often used elsewhere.
 */
void beliefMPC_LA_INEQ_B_GRAD_5_5_5(beliefMPC_FLOAT *lu, beliefMPC_FLOAT *su, beliefMPC_FLOAT *ru, beliefMPC_FLOAT *ll, beliefMPC_FLOAT *sl, beliefMPC_FLOAT *rl, int* lbIdx, int* ubIdx, beliefMPC_FLOAT *grad, beliefMPC_FLOAT *lubysu, beliefMPC_FLOAT *llbysl)
{
	int i;
	for( i=0; i<5; i++ ){
		grad[i] = 0;
	}
	for( i=0; i<5; i++ ){		
		llbysl[i] = ll[i] / sl[i];
		grad[lbIdx[i]] -= llbysl[i]*rl[i];
	}
	for( i=0; i<5; i++ ){
		lubysu[i] = lu[i] / su[i];
		grad[ubIdx[i]] += lubysu[i]*ru[i];
	}
}


/*
 * Addition of three vectors  z = u + w + v
 * of length 68.
 */
void beliefMPC_LA_VVADD3_68(beliefMPC_FLOAT *u, beliefMPC_FLOAT *v, beliefMPC_FLOAT *w, beliefMPC_FLOAT *z)
{
	int i;
	for( i=0; i<68; i++ ){
		z[i] = u[i] + v[i] + w[i];
	}
}


/*
 * Special function to compute the diagonal cholesky factorization of the 
 * positive definite augmented Hessian for block size 7.
 *
 * Inputs: - H = diagonal cost Hessian in diagonal storage format
 *         - llbysl = L / S of lower bounds
 *         - lubysu = L / S of upper bounds
 *
 * Output: Phi = sqrt(H + diag(llbysl) + diag(lubysu))
 * where Phi is stored in diagonal storage format
 */
void beliefMPC_LA_DIAG_CHOL_ONELOOP_LBUB_7_7_7(beliefMPC_FLOAT *H, beliefMPC_FLOAT *llbysl, int* lbIdx, beliefMPC_FLOAT *lubysu, int* ubIdx, beliefMPC_FLOAT *Phi)


{
	int i;
	
	/* compute cholesky */
	for( i=0; i<7; i++ ){
		Phi[i] = H[i] + llbysl[i] + lubysu[i];

#if beliefMPC_SET_PRINTLEVEL > 0 && defined PRINTNUMERICALWARNINGS
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
 * where A is to be computed and is of size [10 x 7],
 * B is given and of size [10 x 7], L is a diagonal
 * matrix of size 10 stored in diagonal matrix 
 * storage format. Note the transpose of L has no impact!
 *
 * Result: A in column major storage format.
 *
 */
void beliefMPC_LA_DIAG_MATRIXFORWARDSUB_10_7(beliefMPC_FLOAT *L, beliefMPC_FLOAT *B, beliefMPC_FLOAT *A)
{
    int i,j;
	 int k = 0;

	for( j=0; j<7; j++){
		for( i=0; i<10; i++){
			A[k] = B[k]/L[j];
			k++;
		}
	}

}


/**
 * Forward substitution to solve L*y = b where L is a
 * diagonal matrix in vector storage format.
 * 
 * The dimensions involved are 7.
 */
void beliefMPC_LA_DIAG_FORWARDSUB_7(beliefMPC_FLOAT *L, beliefMPC_FLOAT *b, beliefMPC_FLOAT *y)
{
    int i;

    for( i=0; i<7; i++ ){
		y[i] = b[i]/L[i];
    }
}


/**
 * Forward substitution for the matrix equation A*L' = B
 * where A is to be computed and is of size [5 x 7],
 * B is given and of size [5 x 7], L is a diagonal
 * matrix of size 5 stored in diagonal matrix 
 * storage format. Note the transpose of L has no impact!
 *
 * Result: A in column major storage format.
 *
 */
void beliefMPC_LA_DIAG_MATRIXFORWARDSUB_5_7(beliefMPC_FLOAT *L, beliefMPC_FLOAT *B, beliefMPC_FLOAT *A)
{
    int i,j;
	 int k = 0;

	for( j=0; j<7; j++){
		for( i=0; i<5; i++){
			A[k] = B[k]/L[j];
			k++;
		}
	}

}


/**
 * Compute C = A*B' where 
 *
 *	size(A) = [10 x 7]
 *  size(B) = [5 x 7]
 * 
 * and all matrices are stored in column major format.
 *
 * THIS ONE HAS THE WORST ACCES PATTERN POSSIBLE.  
 * 
 */
void beliefMPC_LA_DENSE_MMTM_10_7_5(beliefMPC_FLOAT *A, beliefMPC_FLOAT *B, beliefMPC_FLOAT *C)
{
    int i, j, k;
    beliefMPC_FLOAT temp;
    
    for( i=0; i<10; i++ ){        
        for( j=0; j<5; j++ ){
            temp = 0; 
            for( k=0; k<7; k++ ){
                temp += A[k*10+i]*B[k*5+j];
            }						
            C[j*10+i] = temp;
        }
    }
}


/**
 * Forward substitution for the matrix equation A*L' = B
 * where A is to be computed and is of size [5 x 7],
 * B is given and of size [5 x 7], L is a diagonal
 *  matrix of size 7 stored in diagonal 
 * storage format. Note the transpose of L!
 *
 * Result: A in diagonalzero storage format.
 *
 */
void beliefMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_5_7(beliefMPC_FLOAT *L, beliefMPC_FLOAT *B, beliefMPC_FLOAT *A)
{
	int j;
    for( j=0; j<7; j++ ){   
		A[j] = B[j]/L[j];
     }
}


/**
 * Compute C = A*B' where 
 *
 *	size(A) = [5 x 7]
 *  size(B) = [5 x 7] in diagzero format
 * 
 * A and C matrices are stored in column major format.
 * 
 * 
 */
void beliefMPC_LA_DENSE_DIAGZERO_MMTM_5_7_5(beliefMPC_FLOAT *A, beliefMPC_FLOAT *B, beliefMPC_FLOAT *C)
{
    int i, j;
	
	for( i=0; i<5; i++ ){
		for( j=0; j<5; j++){
			C[j*5+i] = B[i*5+j]*A[i];
		}
	}

}


/*
 * Special function to compute the diagonal cholesky factorization of the 
 * positive definite augmented Hessian for block size 5.
 *
 * Inputs: - H = diagonal cost Hessian in diagonal storage format
 *         - llbysl = L / S of lower bounds
 *         - lubysu = L / S of upper bounds
 *
 * Output: Phi = sqrt(H + diag(llbysl) + diag(lubysu))
 * where Phi is stored in diagonal storage format
 */
void beliefMPC_LA_DIAG_CHOL_ONELOOP_LBUB_5_5_5(beliefMPC_FLOAT *H, beliefMPC_FLOAT *llbysl, int* lbIdx, beliefMPC_FLOAT *lubysu, int* ubIdx, beliefMPC_FLOAT *Phi)


{
	int i;
	
	/* compute cholesky */
	for( i=0; i<5; i++ ){
		Phi[i] = H[i] + llbysl[i] + lubysu[i];

#if beliefMPC_SET_PRINTLEVEL > 0 && defined PRINTNUMERICALWARNINGS
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
 * where A is to be computed and is of size [5 x 5],
 * B is given and of size [5 x 5], L is a diagonal
 *  matrix of size 5 stored in diagonal 
 * storage format. Note the transpose of L!
 *
 * Result: A in diagonalzero storage format.
 *
 */
void beliefMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_5_5(beliefMPC_FLOAT *L, beliefMPC_FLOAT *B, beliefMPC_FLOAT *A)
{
	int j;
    for( j=0; j<5; j++ ){   
		A[j] = B[j]/L[j];
     }
}


/**
 * Forward substitution to solve L*y = b where L is a
 * diagonal matrix in vector storage format.
 * 
 * The dimensions involved are 5.
 */
void beliefMPC_LA_DIAG_FORWARDSUB_5(beliefMPC_FLOAT *L, beliefMPC_FLOAT *b, beliefMPC_FLOAT *y)
{
    int i;

    for( i=0; i<5; i++ ){
		y[i] = b[i]/L[i];
    }
}


/**
 * Compute L = A*A' + B*B', where L is lower triangular of size NXp1
 * and A is a dense matrix of size [10 x 7] in column
 * storage format, and B is of size [10 x 7] also in column
 * storage format.
 * 
 * THIS ONE HAS THE WORST ACCES PATTERN POSSIBLE. 
 * POSSIBKE FIX: PUT A AND B INTO ROW MAJOR FORMAT FIRST.
 * 
 */
void beliefMPC_LA_DENSE_MMT2_10_7_7(beliefMPC_FLOAT *A, beliefMPC_FLOAT *B, beliefMPC_FLOAT *L)
{
    int i, j, k, ii, di;
    beliefMPC_FLOAT ltemp;
    
    ii = 0; di = 0;
    for( i=0; i<10; i++ ){        
        for( j=0; j<=i; j++ ){
            ltemp = 0; 
            for( k=0; k<7; k++ ){
                ltemp += A[k*10+i]*A[k*10+j];
            }			
			for( k=0; k<7; k++ ){
                ltemp += B[k*10+i]*B[k*10+j];
            }
            L[ii+j] = ltemp;
        }
        ii += ++di;
    }
}


/* 
 * Computes r = b - A*x - B*u
 * where A an B are stored in column major format
 */
void beliefMPC_LA_DENSE_MVMSUB2_10_7_7(beliefMPC_FLOAT *A, beliefMPC_FLOAT *x, beliefMPC_FLOAT *B, beliefMPC_FLOAT *u, beliefMPC_FLOAT *b, beliefMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;
	int m = 0;
	int n;

	for( i=0; i<10; i++ ){
		r[i] = b[i] - A[k++]*x[0] - B[m++]*u[0];
	}	
	for( j=1; j<7; j++ ){		
		for( i=0; i<10; i++ ){
			r[i] -= A[k++]*x[j];
		}
	}
	
	for( n=1; n<7; n++ ){
		for( i=0; i<10; i++ ){
			r[i] -= B[m++]*u[n];
		}		
	}
}


/**
 * Compute L = A*A' + B*B', where L is lower triangular of size NXp1
 * and A is a dense matrix of size [5 x 7] in column
 * storage format, and B is of size [5 x 7] diagonalzero
 * storage format.
 * 
 * THIS ONE HAS THE WORST ACCES PATTERN POSSIBLE. 
 * POSSIBKE FIX: PUT A AND B INTO ROW MAJOR FORMAT FIRST.
 * 
 */
void beliefMPC_LA_DENSE_DIAGZERO_MMT2_5_7_7(beliefMPC_FLOAT *A, beliefMPC_FLOAT *B, beliefMPC_FLOAT *L)
{
    int i, j, k, ii, di;
    beliefMPC_FLOAT ltemp;
    
    ii = 0; di = 0;
    for( i=0; i<5; i++ ){        
        for( j=0; j<=i; j++ ){
            ltemp = 0; 
            for( k=0; k<7; k++ ){
                ltemp += A[k*5+i]*A[k*5+j];
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
void beliefMPC_LA_DENSE_DIAGZERO_2MVMSUB2_5_7_7(beliefMPC_FLOAT *A, beliefMPC_FLOAT *x, beliefMPC_FLOAT *B, beliefMPC_FLOAT *u, beliefMPC_FLOAT *b, beliefMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;

	for( i=0; i<5; i++ ){
		r[i] = b[i] - A[k++]*x[0] - B[i]*u[i];
	}	

	for( j=1; j<7; j++ ){		
		for( i=0; i<5; i++ ){
			r[i] -= A[k++]*x[j];
		}
	}
	
}


/**
 * Compute L = A*A' + B*B', where L is lower triangular of size NXp1
 * and A is a dense matrix of size [5 x 7] in column
 * storage format, and B is of size [5 x 5] diagonalzero
 * storage format.
 * 
 * THIS ONE HAS THE WORST ACCES PATTERN POSSIBLE. 
 * POSSIBKE FIX: PUT A AND B INTO ROW MAJOR FORMAT FIRST.
 * 
 */
void beliefMPC_LA_DENSE_DIAGZERO_MMT2_5_7_5(beliefMPC_FLOAT *A, beliefMPC_FLOAT *B, beliefMPC_FLOAT *L)
{
    int i, j, k, ii, di;
    beliefMPC_FLOAT ltemp;
    
    ii = 0; di = 0;
    for( i=0; i<5; i++ ){        
        for( j=0; j<=i; j++ ){
            ltemp = 0; 
            for( k=0; k<7; k++ ){
                ltemp += A[k*5+i]*A[k*5+j];
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
void beliefMPC_LA_DENSE_DIAGZERO_2MVMSUB2_5_7_5(beliefMPC_FLOAT *A, beliefMPC_FLOAT *x, beliefMPC_FLOAT *B, beliefMPC_FLOAT *u, beliefMPC_FLOAT *b, beliefMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;

	for( i=0; i<5; i++ ){
		r[i] = b[i] - A[k++]*x[0] - B[i]*u[i];
	}	

	for( j=1; j<7; j++ ){		
		for( i=0; i<5; i++ ){
			r[i] -= A[k++]*x[j];
		}
	}
	
}


/**
 * Cholesky factorization as above, but working on a matrix in 
 * lower triangular storage format of size 10 and outputting
 * the Cholesky factor to matrix L in lower triangular format.
 */
void beliefMPC_LA_DENSE_CHOL_10(beliefMPC_FLOAT *A, beliefMPC_FLOAT *L)
{
    int i, j, k, di, dj;
	 int ii, jj;

    beliefMPC_FLOAT l;
    beliefMPC_FLOAT Mii;

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
        
#if beliefMPC_SET_PRINTLEVEL > 0 && defined PRINTNUMERICALWARNINGS
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
void beliefMPC_LA_DENSE_FORWARDSUB_10(beliefMPC_FLOAT *L, beliefMPC_FLOAT *b, beliefMPC_FLOAT *y)
{
    int i,j,ii,di;
    beliefMPC_FLOAT yel;
            
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
 * where A is to be computed and is of size [5 x 10],
 * B is given and of size [5 x 10], L is a lower tri-
 * angular matrix of size 10 stored in lower triangular 
 * storage format. Note the transpose of L AND B!
 *
 * Result: A in column major storage format.
 *
 */
void beliefMPC_LA_DENSE_MATRIXTFORWARDSUB_5_10(beliefMPC_FLOAT *L, beliefMPC_FLOAT *B, beliefMPC_FLOAT *A)
{
    int i,j,k,ii,di;
    beliefMPC_FLOAT a;
    
    ii=0; di=0;
    for( j=0; j<10; j++ ){        
        for( i=0; i<5; i++ ){
            a = B[i*10+j];
            for( k=0; k<j; k++ ){
                a -= A[k*5+i]*L[ii+k];
            }    

			/* saturate for numerical stability */
			a = MIN(a, BIGM);
			a = MAX(a, -BIGM); 

			A[j*5+i] = a/L[ii+j];			
        }
        ii += ++di;
    }
}


/**
 * Compute L = L - A*A', where L is lower triangular of size 5
 * and A is a dense matrix of size [5 x 10] in column
 * storage format.
 * 
 * THIS ONE HAS THE WORST ACCES PATTERN POSSIBLE. 
 * POSSIBKE FIX: PUT A INTO ROW MAJOR FORMAT FIRST.
 * 
 */
void beliefMPC_LA_DENSE_MMTSUB_5_10(beliefMPC_FLOAT *A, beliefMPC_FLOAT *L)
{
    int i, j, k, ii, di;
    beliefMPC_FLOAT ltemp;
    
    ii = 0; di = 0;
    for( i=0; i<5; i++ ){        
        for( j=0; j<=i; j++ ){
            ltemp = 0; 
            for( k=0; k<10; k++ ){
                ltemp += A[k*5+i]*A[k*5+j];
            }						
            L[ii+j] -= ltemp;
        }
        ii += ++di;
    }
}


/**
 * Cholesky factorization as above, but working on a matrix in 
 * lower triangular storage format of size 5 and outputting
 * the Cholesky factor to matrix L in lower triangular format.
 */
void beliefMPC_LA_DENSE_CHOL_5(beliefMPC_FLOAT *A, beliefMPC_FLOAT *L)
{
    int i, j, k, di, dj;
	 int ii, jj;

    beliefMPC_FLOAT l;
    beliefMPC_FLOAT Mii;

	/* copy A to L first and then operate on L */
	/* COULD BE OPTIMIZED */
	ii=0; di=0;
	for( i=0; i<5; i++ ){
		for( j=0; j<=i; j++ ){
			L[ii+j] = A[ii+j];
		}
		ii += ++di;
	}    
	
	/* factor L */
	ii=0; di=0;
    for( i=0; i<5; i++ ){
        l = 0;
        for( k=0; k<i; k++ ){
            l += L[ii+k]*L[ii+k];
        }        
        
        Mii = L[ii+i] - l;
        
#if beliefMPC_SET_PRINTLEVEL > 0 && defined PRINTNUMERICALWARNINGS
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
        for( j=i+1; j<5; j++ ){
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


/* 
 * Computes r = b - A*x
 * where A is stored in column major format
 */
void beliefMPC_LA_DENSE_MVMSUB1_5_10(beliefMPC_FLOAT *A, beliefMPC_FLOAT *x, beliefMPC_FLOAT *b, beliefMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;

	for( i=0; i<5; i++ ){
		r[i] = b[i] - A[k++]*x[0];
	}	
	for( j=1; j<10; j++ ){		
		for( i=0; i<5; i++ ){
			r[i] -= A[k++]*x[j];
		}
	}
}


/**
 * Forward substitution to solve L*y = b where L is a
 * lower triangular matrix in triangular storage format.
 * 
 * The dimensions involved are 5.
 */
void beliefMPC_LA_DENSE_FORWARDSUB_5(beliefMPC_FLOAT *L, beliefMPC_FLOAT *b, beliefMPC_FLOAT *y)
{
    int i,j,ii,di;
    beliefMPC_FLOAT yel;
            
    ii = 0; di = 0;
    for( i=0; i<5; i++ ){
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
 * where A is to be computed and is of size [5 x 5],
 * B is given and of size [5 x 5], L is a lower tri-
 * angular matrix of size 5 stored in lower triangular 
 * storage format. Note the transpose of L AND B!
 *
 * Result: A in column major storage format.
 *
 */
void beliefMPC_LA_DENSE_MATRIXTFORWARDSUB_5_5(beliefMPC_FLOAT *L, beliefMPC_FLOAT *B, beliefMPC_FLOAT *A)
{
    int i,j,k,ii,di;
    beliefMPC_FLOAT a;
    
    ii=0; di=0;
    for( j=0; j<5; j++ ){        
        for( i=0; i<5; i++ ){
            a = B[i*5+j];
            for( k=0; k<j; k++ ){
                a -= A[k*5+i]*L[ii+k];
            }    

			/* saturate for numerical stability */
			a = MIN(a, BIGM);
			a = MAX(a, -BIGM); 

			A[j*5+i] = a/L[ii+j];			
        }
        ii += ++di;
    }
}


/**
 * Compute L = L - A*A', where L is lower triangular of size 5
 * and A is a dense matrix of size [5 x 5] in column
 * storage format.
 * 
 * THIS ONE HAS THE WORST ACCES PATTERN POSSIBLE. 
 * POSSIBKE FIX: PUT A INTO ROW MAJOR FORMAT FIRST.
 * 
 */
void beliefMPC_LA_DENSE_MMTSUB_5_5(beliefMPC_FLOAT *A, beliefMPC_FLOAT *L)
{
    int i, j, k, ii, di;
    beliefMPC_FLOAT ltemp;
    
    ii = 0; di = 0;
    for( i=0; i<5; i++ ){        
        for( j=0; j<=i; j++ ){
            ltemp = 0; 
            for( k=0; k<5; k++ ){
                ltemp += A[k*5+i]*A[k*5+j];
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
void beliefMPC_LA_DENSE_MVMSUB1_5_5(beliefMPC_FLOAT *A, beliefMPC_FLOAT *x, beliefMPC_FLOAT *b, beliefMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;

	for( i=0; i<5; i++ ){
		r[i] = b[i] - A[k++]*x[0];
	}	
	for( j=1; j<5; j++ ){		
		for( i=0; i<5; i++ ){
			r[i] -= A[k++]*x[j];
		}
	}
}


/**
 * Backward Substitution to solve L^T*x = y where L is a
 * lower triangular matrix in triangular storage format.
 * 
 * All involved dimensions are 5.
 */
void beliefMPC_LA_DENSE_BACKWARDSUB_5(beliefMPC_FLOAT *L, beliefMPC_FLOAT *y, beliefMPC_FLOAT *x)
{
    int i, ii, di, j, jj, dj;
    beliefMPC_FLOAT xel;    
	int start = 10;
    
    /* now solve L^T*x = y by backward substitution */
    ii = start; di = 4;
    for( i=4; i>=0; i-- ){        
        xel = y[i];        
        jj = start; dj = 4;
        for( j=4; j>i; j-- ){
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
 * Matrix vector multiplication y = b - M'*x where M is of size [5 x 5]
 * and stored in column major format. Note the transpose of M!
 */
void beliefMPC_LA_DENSE_MTVMSUB_5_5(beliefMPC_FLOAT *A, beliefMPC_FLOAT *x, beliefMPC_FLOAT *b, beliefMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0; 
	for( i=0; i<5; i++ ){
		r[i] = b[i];
		for( j=0; j<5; j++ ){
			r[i] -= A[k++]*x[j];
		}
	}
}


/*
 * Matrix vector multiplication y = b - M'*x where M is of size [5 x 10]
 * and stored in column major format. Note the transpose of M!
 */
void beliefMPC_LA_DENSE_MTVMSUB_5_10(beliefMPC_FLOAT *A, beliefMPC_FLOAT *x, beliefMPC_FLOAT *b, beliefMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0; 
	for( i=0; i<10; i++ ){
		r[i] = b[i];
		for( j=0; j<5; j++ ){
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
void beliefMPC_LA_DENSE_BACKWARDSUB_10(beliefMPC_FLOAT *L, beliefMPC_FLOAT *y, beliefMPC_FLOAT *x)
{
    int i, ii, di, j, jj, dj;
    beliefMPC_FLOAT xel;    
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
 * Vector subtraction z = -x - y for vectors of length 68.
 */
void beliefMPC_LA_VSUB2_68(beliefMPC_FLOAT *x, beliefMPC_FLOAT *y, beliefMPC_FLOAT *z)
{
	int i;
	for( i=0; i<68; i++){
		z[i] = -x[i] - y[i];
	}
}


/**
 * Forward-Backward-Substitution to solve L*L^T*x = b where L is a
 * diagonal matrix of size 7 in vector
 * storage format.
 */
void beliefMPC_LA_DIAG_FORWARDBACKWARDSUB_7(beliefMPC_FLOAT *L, beliefMPC_FLOAT *b, beliefMPC_FLOAT *x)
{
    int i;
            
    /* solve Ly = b by forward and backward substitution */
    for( i=0; i<7; i++ ){
		x[i] = b[i]/(L[i]*L[i]);
    }
    
}


/**
 * Forward-Backward-Substitution to solve L*L^T*x = b where L is a
 * diagonal matrix of size 5 in vector
 * storage format.
 */
void beliefMPC_LA_DIAG_FORWARDBACKWARDSUB_5(beliefMPC_FLOAT *L, beliefMPC_FLOAT *b, beliefMPC_FLOAT *x)
{
    int i;
            
    /* solve Ly = b by forward and backward substitution */
    for( i=0; i<5; i++ ){
		x[i] = b[i]/(L[i]*L[i]);
    }
    
}


/*
 * Vector subtraction z = x(xidx) - y where y, z and xidx are of length 7,
 * and x has length 7 and is indexed through yidx.
 */
void beliefMPC_LA_VSUB_INDEXED_7(beliefMPC_FLOAT *x, int* xidx, beliefMPC_FLOAT *y, beliefMPC_FLOAT *z)
{
	int i;
	for( i=0; i<7; i++){
		z[i] = x[xidx[i]] - y[i];
	}
}


/*
 * Vector subtraction x = -u.*v - w for vectors of length 7.
 */
void beliefMPC_LA_VSUB3_7(beliefMPC_FLOAT *u, beliefMPC_FLOAT *v, beliefMPC_FLOAT *w, beliefMPC_FLOAT *x)
{
	int i;
	for( i=0; i<7; i++){
		x[i] = -u[i]*v[i] - w[i];
	}
}


/*
 * Vector subtraction z = -x - y(yidx) where y is of length 7
 * and z, x and yidx are of length 7.
 */
void beliefMPC_LA_VSUB2_INDEXED_7(beliefMPC_FLOAT *x, beliefMPC_FLOAT *y, int* yidx, beliefMPC_FLOAT *z)
{
	int i;
	for( i=0; i<7; i++){
		z[i] = -x[i] - y[yidx[i]];
	}
}


/*
 * Vector subtraction z = x(xidx) - y where y, z and xidx are of length 5,
 * and x has length 5 and is indexed through yidx.
 */
void beliefMPC_LA_VSUB_INDEXED_5(beliefMPC_FLOAT *x, int* xidx, beliefMPC_FLOAT *y, beliefMPC_FLOAT *z)
{
	int i;
	for( i=0; i<5; i++){
		z[i] = x[xidx[i]] - y[i];
	}
}


/*
 * Vector subtraction x = -u.*v - w for vectors of length 5.
 */
void beliefMPC_LA_VSUB3_5(beliefMPC_FLOAT *u, beliefMPC_FLOAT *v, beliefMPC_FLOAT *w, beliefMPC_FLOAT *x)
{
	int i;
	for( i=0; i<5; i++){
		x[i] = -u[i]*v[i] - w[i];
	}
}


/*
 * Vector subtraction z = -x - y(yidx) where y is of length 5
 * and z, x and yidx are of length 5.
 */
void beliefMPC_LA_VSUB2_INDEXED_5(beliefMPC_FLOAT *x, beliefMPC_FLOAT *y, int* yidx, beliefMPC_FLOAT *z)
{
	int i;
	for( i=0; i<5; i++){
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
 * beliefMPC_NOPROGRESS (should be negative).
 */
int beliefMPC_LINESEARCH_BACKTRACKING_AFFINE(beliefMPC_FLOAT *l, beliefMPC_FLOAT *s, beliefMPC_FLOAT *dl, beliefMPC_FLOAT *ds, beliefMPC_FLOAT *a, beliefMPC_FLOAT *mu_aff)
{
    int i;
	int lsIt=1;    
    beliefMPC_FLOAT dltemp;
    beliefMPC_FLOAT dstemp;
    beliefMPC_FLOAT mya = 1.0;
    beliefMPC_FLOAT mymu;
        
    while( 1 ){                        

        /* 
         * Compute both snew and wnew together.
         * We compute also mu_affine along the way here, as the
         * values might be in registers, so it should be cheaper.
         */
        mymu = 0;
        for( i=0; i<136; i++ ){
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
        if( i == 136 ){
            break;
        } else {
            mya *= beliefMPC_SET_LS_SCALE_AFF;
            if( mya < beliefMPC_SET_LS_MINSTEP ){
                return beliefMPC_NOPROGRESS;
            }
        }
    }
    
    /* return new values and iteration counter */
    *a = mya;
    *mu_aff = mymu / (beliefMPC_FLOAT)136;
    return lsIt;
}


/*
 * Vector subtraction x = u.*v - a where a is a scalar
*  and x,u,v are vectors of length 136.
 */
void beliefMPC_LA_VSUB5_136(beliefMPC_FLOAT *u, beliefMPC_FLOAT *v, beliefMPC_FLOAT a, beliefMPC_FLOAT *x)
{
	int i;
	for( i=0; i<136; i++){
		x[i] = u[i]*v[i] - a;
	}
}


/*
 * Computes x=0; x(uidx) += u/su; x(vidx) -= v/sv where x is of length 7,
 * u, su, uidx are of length 7 and v, sv, vidx are of length 7.
 */
void beliefMPC_LA_VSUB6_INDEXED_7_7_7(beliefMPC_FLOAT *u, beliefMPC_FLOAT *su, int* uidx, beliefMPC_FLOAT *v, beliefMPC_FLOAT *sv, int* vidx, beliefMPC_FLOAT *x)
{
	int i;
	for( i=0; i<7; i++ ){
		x[i] = 0;
	}
	for( i=0; i<7; i++){
		x[uidx[i]] += u[i]/su[i];
	}
	for( i=0; i<7; i++){
		x[vidx[i]] -= v[i]/sv[i];
	}
}


/* 
 * Computes r = A*x + B*u
 * where A an B are stored in column major format
 */
void beliefMPC_LA_DENSE_2MVMADD_10_7_7(beliefMPC_FLOAT *A, beliefMPC_FLOAT *x, beliefMPC_FLOAT *B, beliefMPC_FLOAT *u, beliefMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;
	int m = 0;
	int n;

	for( i=0; i<10; i++ ){
		r[i] = A[k++]*x[0] + B[m++]*u[0];
	}	

	for( j=1; j<7; j++ ){		
		for( i=0; i<10; i++ ){
			r[i] += A[k++]*x[j];
		}
	}
	
	for( n=1; n<7; n++ ){
		for( i=0; i<10; i++ ){
			r[i] += B[m++]*u[n];
		}		
	}
}


/* 
 * Computes r = A*x + B*u
 * where A is stored in column major format
 * and B is stored in diagzero format
 */
void beliefMPC_LA_DENSE_DIAGZERO_2MVMADD_5_7_7(beliefMPC_FLOAT *A, beliefMPC_FLOAT *x, beliefMPC_FLOAT *B, beliefMPC_FLOAT *u, beliefMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;

	for( i=0; i<5; i++ ){
		r[i] = A[k++]*x[0] + B[i]*u[i];
	}	

	for( j=1; j<7; j++ ){		
		for( i=0; i<5; i++ ){
			r[i] += A[k++]*x[j];
		}
	}
	
}


/*
 * Computes x=0; x(uidx) += u/su; x(vidx) -= v/sv where x is of length 5,
 * u, su, uidx are of length 5 and v, sv, vidx are of length 5.
 */
void beliefMPC_LA_VSUB6_INDEXED_5_5_5(beliefMPC_FLOAT *u, beliefMPC_FLOAT *su, int* uidx, beliefMPC_FLOAT *v, beliefMPC_FLOAT *sv, int* vidx, beliefMPC_FLOAT *x)
{
	int i;
	for( i=0; i<5; i++ ){
		x[i] = 0;
	}
	for( i=0; i<5; i++){
		x[uidx[i]] += u[i]/su[i];
	}
	for( i=0; i<5; i++){
		x[vidx[i]] -= v[i]/sv[i];
	}
}


/* 
 * Computes r = A*x + B*u
 * where A is stored in column major format
 * and B is stored in diagzero format
 */
void beliefMPC_LA_DENSE_DIAGZERO_2MVMADD_5_7_5(beliefMPC_FLOAT *A, beliefMPC_FLOAT *x, beliefMPC_FLOAT *B, beliefMPC_FLOAT *u, beliefMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;

	for( i=0; i<5; i++ ){
		r[i] = A[k++]*x[0] + B[i]*u[i];
	}	

	for( j=1; j<7; j++ ){		
		for( i=0; i<5; i++ ){
			r[i] += A[k++]*x[j];
		}
	}
	
}


/*
 * Vector subtraction z = x - y for vectors of length 68.
 */
void beliefMPC_LA_VSUB_68(beliefMPC_FLOAT *x, beliefMPC_FLOAT *y, beliefMPC_FLOAT *z)
{
	int i;
	for( i=0; i<68; i++){
		z[i] = x[i] - y[i];
	}
}


/** 
 * Computes z = -r./s - u.*y(y)
 * where all vectors except of y are of length 7 (length of y >= 7).
 */
void beliefMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_7(beliefMPC_FLOAT *r, beliefMPC_FLOAT *s, beliefMPC_FLOAT *u, beliefMPC_FLOAT *y, int* yidx, beliefMPC_FLOAT *z)
{
	int i;
	for( i=0; i<7; i++ ){
		z[i] = -r[i]/s[i] - u[i]*y[yidx[i]];
	}
}


/** 
 * Computes z = -r./s + u.*y(y)
 * where all vectors except of y are of length 7 (length of y >= 7).
 */
void beliefMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_7(beliefMPC_FLOAT *r, beliefMPC_FLOAT *s, beliefMPC_FLOAT *u, beliefMPC_FLOAT *y, int* yidx, beliefMPC_FLOAT *z)
{
	int i;
	for( i=0; i<7; i++ ){
		z[i] = -r[i]/s[i] + u[i]*y[yidx[i]];
	}
}


/** 
 * Computes z = -r./s - u.*y(y)
 * where all vectors except of y are of length 5 (length of y >= 5).
 */
void beliefMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_5(beliefMPC_FLOAT *r, beliefMPC_FLOAT *s, beliefMPC_FLOAT *u, beliefMPC_FLOAT *y, int* yidx, beliefMPC_FLOAT *z)
{
	int i;
	for( i=0; i<5; i++ ){
		z[i] = -r[i]/s[i] - u[i]*y[yidx[i]];
	}
}


/** 
 * Computes z = -r./s + u.*y(y)
 * where all vectors except of y are of length 5 (length of y >= 5).
 */
void beliefMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(beliefMPC_FLOAT *r, beliefMPC_FLOAT *s, beliefMPC_FLOAT *u, beliefMPC_FLOAT *y, int* yidx, beliefMPC_FLOAT *z)
{
	int i;
	for( i=0; i<5; i++ ){
		z[i] = -r[i]/s[i] + u[i]*y[yidx[i]];
	}
}


/*
 * Computes ds = -l.\(r + s.*dl) for vectors of length 136.
 */
void beliefMPC_LA_VSUB7_136(beliefMPC_FLOAT *l, beliefMPC_FLOAT *r, beliefMPC_FLOAT *s, beliefMPC_FLOAT *dl, beliefMPC_FLOAT *ds)
{
	int i;
	for( i=0; i<136; i++){
		ds[i] = -(r[i] + s[i]*dl[i])/l[i];
	}
}


/*
 * Vector addition x = x + y for vectors of length 68.
 */
void beliefMPC_LA_VADD_68(beliefMPC_FLOAT *x, beliefMPC_FLOAT *y)
{
	int i;
	for( i=0; i<68; i++){
		x[i] += y[i];
	}
}


/*
 * Vector addition x = x + y for vectors of length 50.
 */
void beliefMPC_LA_VADD_50(beliefMPC_FLOAT *x, beliefMPC_FLOAT *y)
{
	int i;
	for( i=0; i<50; i++){
		x[i] += y[i];
	}
}


/*
 * Vector addition x = x + y for vectors of length 136.
 */
void beliefMPC_LA_VADD_136(beliefMPC_FLOAT *x, beliefMPC_FLOAT *y)
{
	int i;
	for( i=0; i<136; i++){
		x[i] += y[i];
	}
}


/**
 * Backtracking line search for combined predictor/corrector step.
 * Update on variables with safety factor gamma (to keep us away from
 * boundary).
 */
int beliefMPC_LINESEARCH_BACKTRACKING_COMBINED(beliefMPC_FLOAT *z, beliefMPC_FLOAT *v, beliefMPC_FLOAT *l, beliefMPC_FLOAT *s, beliefMPC_FLOAT *dz, beliefMPC_FLOAT *dv, beliefMPC_FLOAT *dl, beliefMPC_FLOAT *ds, beliefMPC_FLOAT *a, beliefMPC_FLOAT *mu)
{
    int i, lsIt=1;       
    beliefMPC_FLOAT dltemp;
    beliefMPC_FLOAT dstemp;    
    beliefMPC_FLOAT a_gamma;
            
    *a = 1.0;
    while( 1 ){                        

        /* check whether search criterion is fulfilled */
        for( i=0; i<136; i++ ){
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
        if( i == 136 ){
            break;
        } else {
            *a *= beliefMPC_SET_LS_SCALE;
            if( *a < beliefMPC_SET_LS_MINSTEP ){
                return beliefMPC_NOPROGRESS;
            }
        }
    }
    
    /* update variables with safety margin */
    a_gamma = (*a)*beliefMPC_SET_LS_MAXSTEP;
    
    /* primal variables */
    for( i=0; i<68; i++ ){
        z[i] += a_gamma*dz[i];
    }
    
    /* equality constraint multipliers */
    for( i=0; i<50; i++ ){
        v[i] += a_gamma*dv[i];
    }
    
    /* inequality constraint multipliers & slacks, also update mu */
    *mu = 0;
    for( i=0; i<136; i++ ){
        dltemp = l[i] + a_gamma*dl[i]; l[i] = dltemp;
        dstemp = s[i] + a_gamma*ds[i]; s[i] = dstemp;
        *mu += dltemp*dstemp;
    }
    
    *a = a_gamma;
    *mu /= (beliefMPC_FLOAT)136;
    return lsIt;
}




/* VARIABLE DEFINITIONS ------------------------------------------------ */
beliefMPC_FLOAT beliefMPC_z[68];
beliefMPC_FLOAT beliefMPC_v[50];
beliefMPC_FLOAT beliefMPC_dz_aff[68];
beliefMPC_FLOAT beliefMPC_dv_aff[50];
beliefMPC_FLOAT beliefMPC_grad_cost[68];
beliefMPC_FLOAT beliefMPC_grad_eq[68];
beliefMPC_FLOAT beliefMPC_rd[68];
beliefMPC_FLOAT beliefMPC_l[136];
beliefMPC_FLOAT beliefMPC_s[136];
beliefMPC_FLOAT beliefMPC_lbys[136];
beliefMPC_FLOAT beliefMPC_dl_aff[136];
beliefMPC_FLOAT beliefMPC_ds_aff[136];
beliefMPC_FLOAT beliefMPC_dz_cc[68];
beliefMPC_FLOAT beliefMPC_dv_cc[50];
beliefMPC_FLOAT beliefMPC_dl_cc[136];
beliefMPC_FLOAT beliefMPC_ds_cc[136];
beliefMPC_FLOAT beliefMPC_ccrhs[136];
beliefMPC_FLOAT beliefMPC_grad_ineq[68];
beliefMPC_FLOAT beliefMPC_H0[7] = {0.0000000000000000E+000, 0.0000000000000000E+000, 2.0000000000000000E+001, 2.0000000000000000E+001, 2.0000000000000000E+001, 2.0000000000000000E+000, 2.0000000000000000E+000};
beliefMPC_FLOAT beliefMPC_f0[7] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
beliefMPC_FLOAT* beliefMPC_z0 = beliefMPC_z + 0;
beliefMPC_FLOAT* beliefMPC_dzaff0 = beliefMPC_dz_aff + 0;
beliefMPC_FLOAT* beliefMPC_dzcc0 = beliefMPC_dz_cc + 0;
beliefMPC_FLOAT* beliefMPC_rd0 = beliefMPC_rd + 0;
beliefMPC_FLOAT beliefMPC_Lbyrd0[7];
beliefMPC_FLOAT* beliefMPC_grad_cost0 = beliefMPC_grad_cost + 0;
beliefMPC_FLOAT* beliefMPC_grad_eq0 = beliefMPC_grad_eq + 0;
beliefMPC_FLOAT* beliefMPC_grad_ineq0 = beliefMPC_grad_ineq + 0;
beliefMPC_FLOAT beliefMPC_ctv0[7];
beliefMPC_FLOAT* beliefMPC_v0 = beliefMPC_v + 0;
beliefMPC_FLOAT beliefMPC_re0[10];
beliefMPC_FLOAT beliefMPC_beta0[10];
beliefMPC_FLOAT beliefMPC_betacc0[10];
beliefMPC_FLOAT* beliefMPC_dvaff0 = beliefMPC_dv_aff + 0;
beliefMPC_FLOAT* beliefMPC_dvcc0 = beliefMPC_dv_cc + 0;
beliefMPC_FLOAT beliefMPC_V0[70];
beliefMPC_FLOAT beliefMPC_Yd0[55];
beliefMPC_FLOAT beliefMPC_Ld0[55];
beliefMPC_FLOAT beliefMPC_yy0[10];
beliefMPC_FLOAT beliefMPC_bmy0[10];
int beliefMPC_lbIdx0[7] = {0, 1, 2, 3, 4, 5, 6};
beliefMPC_FLOAT* beliefMPC_llb0 = beliefMPC_l + 0;
beliefMPC_FLOAT* beliefMPC_slb0 = beliefMPC_s + 0;
beliefMPC_FLOAT* beliefMPC_llbbyslb0 = beliefMPC_lbys + 0;
beliefMPC_FLOAT beliefMPC_rilb0[7];
beliefMPC_FLOAT* beliefMPC_dllbaff0 = beliefMPC_dl_aff + 0;
beliefMPC_FLOAT* beliefMPC_dslbaff0 = beliefMPC_ds_aff + 0;
beliefMPC_FLOAT* beliefMPC_dllbcc0 = beliefMPC_dl_cc + 0;
beliefMPC_FLOAT* beliefMPC_dslbcc0 = beliefMPC_ds_cc + 0;
beliefMPC_FLOAT* beliefMPC_ccrhsl0 = beliefMPC_ccrhs + 0;
int beliefMPC_ubIdx0[7] = {0, 1, 2, 3, 4, 5, 6};
beliefMPC_FLOAT* beliefMPC_lub0 = beliefMPC_l + 7;
beliefMPC_FLOAT* beliefMPC_sub0 = beliefMPC_s + 7;
beliefMPC_FLOAT* beliefMPC_lubbysub0 = beliefMPC_lbys + 7;
beliefMPC_FLOAT beliefMPC_riub0[7];
beliefMPC_FLOAT* beliefMPC_dlubaff0 = beliefMPC_dl_aff + 7;
beliefMPC_FLOAT* beliefMPC_dsubaff0 = beliefMPC_ds_aff + 7;
beliefMPC_FLOAT* beliefMPC_dlubcc0 = beliefMPC_dl_cc + 7;
beliefMPC_FLOAT* beliefMPC_dsubcc0 = beliefMPC_ds_cc + 7;
beliefMPC_FLOAT* beliefMPC_ccrhsub0 = beliefMPC_ccrhs + 7;
beliefMPC_FLOAT beliefMPC_Phi0[7];
beliefMPC_FLOAT beliefMPC_f1[7] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
beliefMPC_FLOAT* beliefMPC_z1 = beliefMPC_z + 7;
beliefMPC_FLOAT* beliefMPC_dzaff1 = beliefMPC_dz_aff + 7;
beliefMPC_FLOAT* beliefMPC_dzcc1 = beliefMPC_dz_cc + 7;
beliefMPC_FLOAT* beliefMPC_rd1 = beliefMPC_rd + 7;
beliefMPC_FLOAT beliefMPC_Lbyrd1[7];
beliefMPC_FLOAT* beliefMPC_grad_cost1 = beliefMPC_grad_cost + 7;
beliefMPC_FLOAT* beliefMPC_grad_eq1 = beliefMPC_grad_eq + 7;
beliefMPC_FLOAT* beliefMPC_grad_ineq1 = beliefMPC_grad_ineq + 7;
beliefMPC_FLOAT beliefMPC_ctv1[7];
beliefMPC_FLOAT* beliefMPC_v1 = beliefMPC_v + 10;
beliefMPC_FLOAT beliefMPC_re1[5];
beliefMPC_FLOAT beliefMPC_beta1[5];
beliefMPC_FLOAT beliefMPC_betacc1[5];
beliefMPC_FLOAT* beliefMPC_dvaff1 = beliefMPC_dv_aff + 10;
beliefMPC_FLOAT* beliefMPC_dvcc1 = beliefMPC_dv_cc + 10;
beliefMPC_FLOAT beliefMPC_V1[35];
beliefMPC_FLOAT beliefMPC_Yd1[15];
beliefMPC_FLOAT beliefMPC_Ld1[15];
beliefMPC_FLOAT beliefMPC_yy1[5];
beliefMPC_FLOAT beliefMPC_bmy1[5];
int beliefMPC_lbIdx1[7] = {0, 1, 2, 3, 4, 5, 6};
beliefMPC_FLOAT* beliefMPC_llb1 = beliefMPC_l + 14;
beliefMPC_FLOAT* beliefMPC_slb1 = beliefMPC_s + 14;
beliefMPC_FLOAT* beliefMPC_llbbyslb1 = beliefMPC_lbys + 14;
beliefMPC_FLOAT beliefMPC_rilb1[7];
beliefMPC_FLOAT* beliefMPC_dllbaff1 = beliefMPC_dl_aff + 14;
beliefMPC_FLOAT* beliefMPC_dslbaff1 = beliefMPC_ds_aff + 14;
beliefMPC_FLOAT* beliefMPC_dllbcc1 = beliefMPC_dl_cc + 14;
beliefMPC_FLOAT* beliefMPC_dslbcc1 = beliefMPC_ds_cc + 14;
beliefMPC_FLOAT* beliefMPC_ccrhsl1 = beliefMPC_ccrhs + 14;
int beliefMPC_ubIdx1[7] = {0, 1, 2, 3, 4, 5, 6};
beliefMPC_FLOAT* beliefMPC_lub1 = beliefMPC_l + 21;
beliefMPC_FLOAT* beliefMPC_sub1 = beliefMPC_s + 21;
beliefMPC_FLOAT* beliefMPC_lubbysub1 = beliefMPC_lbys + 21;
beliefMPC_FLOAT beliefMPC_riub1[7];
beliefMPC_FLOAT* beliefMPC_dlubaff1 = beliefMPC_dl_aff + 21;
beliefMPC_FLOAT* beliefMPC_dsubaff1 = beliefMPC_ds_aff + 21;
beliefMPC_FLOAT* beliefMPC_dlubcc1 = beliefMPC_dl_cc + 21;
beliefMPC_FLOAT* beliefMPC_dsubcc1 = beliefMPC_ds_cc + 21;
beliefMPC_FLOAT* beliefMPC_ccrhsub1 = beliefMPC_ccrhs + 21;
beliefMPC_FLOAT beliefMPC_Phi1[7];
beliefMPC_FLOAT beliefMPC_D1[70] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, -1.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 
0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, -1.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 
0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, -1.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 
0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, -1.0000000000000000E+000, 0.0000000000000000E+000, 
0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, -1.0000000000000000E+000, 
0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 
0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
beliefMPC_FLOAT beliefMPC_W1[70];
beliefMPC_FLOAT beliefMPC_Ysd1[50];
beliefMPC_FLOAT beliefMPC_Lsd1[50];
beliefMPC_FLOAT beliefMPC_f2[7] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
beliefMPC_FLOAT* beliefMPC_z2 = beliefMPC_z + 14;
beliefMPC_FLOAT* beliefMPC_dzaff2 = beliefMPC_dz_aff + 14;
beliefMPC_FLOAT* beliefMPC_dzcc2 = beliefMPC_dz_cc + 14;
beliefMPC_FLOAT* beliefMPC_rd2 = beliefMPC_rd + 14;
beliefMPC_FLOAT beliefMPC_Lbyrd2[7];
beliefMPC_FLOAT* beliefMPC_grad_cost2 = beliefMPC_grad_cost + 14;
beliefMPC_FLOAT* beliefMPC_grad_eq2 = beliefMPC_grad_eq + 14;
beliefMPC_FLOAT* beliefMPC_grad_ineq2 = beliefMPC_grad_ineq + 14;
beliefMPC_FLOAT beliefMPC_ctv2[7];
beliefMPC_FLOAT* beliefMPC_v2 = beliefMPC_v + 15;
beliefMPC_FLOAT beliefMPC_re2[5];
beliefMPC_FLOAT beliefMPC_beta2[5];
beliefMPC_FLOAT beliefMPC_betacc2[5];
beliefMPC_FLOAT* beliefMPC_dvaff2 = beliefMPC_dv_aff + 15;
beliefMPC_FLOAT* beliefMPC_dvcc2 = beliefMPC_dv_cc + 15;
beliefMPC_FLOAT beliefMPC_V2[35];
beliefMPC_FLOAT beliefMPC_Yd2[15];
beliefMPC_FLOAT beliefMPC_Ld2[15];
beliefMPC_FLOAT beliefMPC_yy2[5];
beliefMPC_FLOAT beliefMPC_bmy2[5];
int beliefMPC_lbIdx2[7] = {0, 1, 2, 3, 4, 5, 6};
beliefMPC_FLOAT* beliefMPC_llb2 = beliefMPC_l + 28;
beliefMPC_FLOAT* beliefMPC_slb2 = beliefMPC_s + 28;
beliefMPC_FLOAT* beliefMPC_llbbyslb2 = beliefMPC_lbys + 28;
beliefMPC_FLOAT beliefMPC_rilb2[7];
beliefMPC_FLOAT* beliefMPC_dllbaff2 = beliefMPC_dl_aff + 28;
beliefMPC_FLOAT* beliefMPC_dslbaff2 = beliefMPC_ds_aff + 28;
beliefMPC_FLOAT* beliefMPC_dllbcc2 = beliefMPC_dl_cc + 28;
beliefMPC_FLOAT* beliefMPC_dslbcc2 = beliefMPC_ds_cc + 28;
beliefMPC_FLOAT* beliefMPC_ccrhsl2 = beliefMPC_ccrhs + 28;
int beliefMPC_ubIdx2[7] = {0, 1, 2, 3, 4, 5, 6};
beliefMPC_FLOAT* beliefMPC_lub2 = beliefMPC_l + 35;
beliefMPC_FLOAT* beliefMPC_sub2 = beliefMPC_s + 35;
beliefMPC_FLOAT* beliefMPC_lubbysub2 = beliefMPC_lbys + 35;
beliefMPC_FLOAT beliefMPC_riub2[7];
beliefMPC_FLOAT* beliefMPC_dlubaff2 = beliefMPC_dl_aff + 35;
beliefMPC_FLOAT* beliefMPC_dsubaff2 = beliefMPC_ds_aff + 35;
beliefMPC_FLOAT* beliefMPC_dlubcc2 = beliefMPC_dl_cc + 35;
beliefMPC_FLOAT* beliefMPC_dsubcc2 = beliefMPC_ds_cc + 35;
beliefMPC_FLOAT* beliefMPC_ccrhsub2 = beliefMPC_ccrhs + 35;
beliefMPC_FLOAT beliefMPC_Phi2[7];
beliefMPC_FLOAT beliefMPC_D2[7] = {-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000};
beliefMPC_FLOAT beliefMPC_W2[7];
beliefMPC_FLOAT beliefMPC_Ysd2[25];
beliefMPC_FLOAT beliefMPC_Lsd2[25];
beliefMPC_FLOAT beliefMPC_f3[7] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
beliefMPC_FLOAT* beliefMPC_z3 = beliefMPC_z + 21;
beliefMPC_FLOAT* beliefMPC_dzaff3 = beliefMPC_dz_aff + 21;
beliefMPC_FLOAT* beliefMPC_dzcc3 = beliefMPC_dz_cc + 21;
beliefMPC_FLOAT* beliefMPC_rd3 = beliefMPC_rd + 21;
beliefMPC_FLOAT beliefMPC_Lbyrd3[7];
beliefMPC_FLOAT* beliefMPC_grad_cost3 = beliefMPC_grad_cost + 21;
beliefMPC_FLOAT* beliefMPC_grad_eq3 = beliefMPC_grad_eq + 21;
beliefMPC_FLOAT* beliefMPC_grad_ineq3 = beliefMPC_grad_ineq + 21;
beliefMPC_FLOAT beliefMPC_ctv3[7];
beliefMPC_FLOAT* beliefMPC_v3 = beliefMPC_v + 20;
beliefMPC_FLOAT beliefMPC_re3[5];
beliefMPC_FLOAT beliefMPC_beta3[5];
beliefMPC_FLOAT beliefMPC_betacc3[5];
beliefMPC_FLOAT* beliefMPC_dvaff3 = beliefMPC_dv_aff + 20;
beliefMPC_FLOAT* beliefMPC_dvcc3 = beliefMPC_dv_cc + 20;
beliefMPC_FLOAT beliefMPC_V3[35];
beliefMPC_FLOAT beliefMPC_Yd3[15];
beliefMPC_FLOAT beliefMPC_Ld3[15];
beliefMPC_FLOAT beliefMPC_yy3[5];
beliefMPC_FLOAT beliefMPC_bmy3[5];
int beliefMPC_lbIdx3[7] = {0, 1, 2, 3, 4, 5, 6};
beliefMPC_FLOAT* beliefMPC_llb3 = beliefMPC_l + 42;
beliefMPC_FLOAT* beliefMPC_slb3 = beliefMPC_s + 42;
beliefMPC_FLOAT* beliefMPC_llbbyslb3 = beliefMPC_lbys + 42;
beliefMPC_FLOAT beliefMPC_rilb3[7];
beliefMPC_FLOAT* beliefMPC_dllbaff3 = beliefMPC_dl_aff + 42;
beliefMPC_FLOAT* beliefMPC_dslbaff3 = beliefMPC_ds_aff + 42;
beliefMPC_FLOAT* beliefMPC_dllbcc3 = beliefMPC_dl_cc + 42;
beliefMPC_FLOAT* beliefMPC_dslbcc3 = beliefMPC_ds_cc + 42;
beliefMPC_FLOAT* beliefMPC_ccrhsl3 = beliefMPC_ccrhs + 42;
int beliefMPC_ubIdx3[7] = {0, 1, 2, 3, 4, 5, 6};
beliefMPC_FLOAT* beliefMPC_lub3 = beliefMPC_l + 49;
beliefMPC_FLOAT* beliefMPC_sub3 = beliefMPC_s + 49;
beliefMPC_FLOAT* beliefMPC_lubbysub3 = beliefMPC_lbys + 49;
beliefMPC_FLOAT beliefMPC_riub3[7];
beliefMPC_FLOAT* beliefMPC_dlubaff3 = beliefMPC_dl_aff + 49;
beliefMPC_FLOAT* beliefMPC_dsubaff3 = beliefMPC_ds_aff + 49;
beliefMPC_FLOAT* beliefMPC_dlubcc3 = beliefMPC_dl_cc + 49;
beliefMPC_FLOAT* beliefMPC_dsubcc3 = beliefMPC_ds_cc + 49;
beliefMPC_FLOAT* beliefMPC_ccrhsub3 = beliefMPC_ccrhs + 49;
beliefMPC_FLOAT beliefMPC_Phi3[7];
beliefMPC_FLOAT beliefMPC_W3[7];
beliefMPC_FLOAT beliefMPC_Ysd3[25];
beliefMPC_FLOAT beliefMPC_Lsd3[25];
beliefMPC_FLOAT beliefMPC_f4[7] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
beliefMPC_FLOAT* beliefMPC_z4 = beliefMPC_z + 28;
beliefMPC_FLOAT* beliefMPC_dzaff4 = beliefMPC_dz_aff + 28;
beliefMPC_FLOAT* beliefMPC_dzcc4 = beliefMPC_dz_cc + 28;
beliefMPC_FLOAT* beliefMPC_rd4 = beliefMPC_rd + 28;
beliefMPC_FLOAT beliefMPC_Lbyrd4[7];
beliefMPC_FLOAT* beliefMPC_grad_cost4 = beliefMPC_grad_cost + 28;
beliefMPC_FLOAT* beliefMPC_grad_eq4 = beliefMPC_grad_eq + 28;
beliefMPC_FLOAT* beliefMPC_grad_ineq4 = beliefMPC_grad_ineq + 28;
beliefMPC_FLOAT beliefMPC_ctv4[7];
beliefMPC_FLOAT* beliefMPC_v4 = beliefMPC_v + 25;
beliefMPC_FLOAT beliefMPC_re4[5];
beliefMPC_FLOAT beliefMPC_beta4[5];
beliefMPC_FLOAT beliefMPC_betacc4[5];
beliefMPC_FLOAT* beliefMPC_dvaff4 = beliefMPC_dv_aff + 25;
beliefMPC_FLOAT* beliefMPC_dvcc4 = beliefMPC_dv_cc + 25;
beliefMPC_FLOAT beliefMPC_V4[35];
beliefMPC_FLOAT beliefMPC_Yd4[15];
beliefMPC_FLOAT beliefMPC_Ld4[15];
beliefMPC_FLOAT beliefMPC_yy4[5];
beliefMPC_FLOAT beliefMPC_bmy4[5];
int beliefMPC_lbIdx4[7] = {0, 1, 2, 3, 4, 5, 6};
beliefMPC_FLOAT* beliefMPC_llb4 = beliefMPC_l + 56;
beliefMPC_FLOAT* beliefMPC_slb4 = beliefMPC_s + 56;
beliefMPC_FLOAT* beliefMPC_llbbyslb4 = beliefMPC_lbys + 56;
beliefMPC_FLOAT beliefMPC_rilb4[7];
beliefMPC_FLOAT* beliefMPC_dllbaff4 = beliefMPC_dl_aff + 56;
beliefMPC_FLOAT* beliefMPC_dslbaff4 = beliefMPC_ds_aff + 56;
beliefMPC_FLOAT* beliefMPC_dllbcc4 = beliefMPC_dl_cc + 56;
beliefMPC_FLOAT* beliefMPC_dslbcc4 = beliefMPC_ds_cc + 56;
beliefMPC_FLOAT* beliefMPC_ccrhsl4 = beliefMPC_ccrhs + 56;
int beliefMPC_ubIdx4[7] = {0, 1, 2, 3, 4, 5, 6};
beliefMPC_FLOAT* beliefMPC_lub4 = beliefMPC_l + 63;
beliefMPC_FLOAT* beliefMPC_sub4 = beliefMPC_s + 63;
beliefMPC_FLOAT* beliefMPC_lubbysub4 = beliefMPC_lbys + 63;
beliefMPC_FLOAT beliefMPC_riub4[7];
beliefMPC_FLOAT* beliefMPC_dlubaff4 = beliefMPC_dl_aff + 63;
beliefMPC_FLOAT* beliefMPC_dsubaff4 = beliefMPC_ds_aff + 63;
beliefMPC_FLOAT* beliefMPC_dlubcc4 = beliefMPC_dl_cc + 63;
beliefMPC_FLOAT* beliefMPC_dsubcc4 = beliefMPC_ds_cc + 63;
beliefMPC_FLOAT* beliefMPC_ccrhsub4 = beliefMPC_ccrhs + 63;
beliefMPC_FLOAT beliefMPC_Phi4[7];
beliefMPC_FLOAT beliefMPC_W4[7];
beliefMPC_FLOAT beliefMPC_Ysd4[25];
beliefMPC_FLOAT beliefMPC_Lsd4[25];
beliefMPC_FLOAT beliefMPC_f5[7] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
beliefMPC_FLOAT* beliefMPC_z5 = beliefMPC_z + 35;
beliefMPC_FLOAT* beliefMPC_dzaff5 = beliefMPC_dz_aff + 35;
beliefMPC_FLOAT* beliefMPC_dzcc5 = beliefMPC_dz_cc + 35;
beliefMPC_FLOAT* beliefMPC_rd5 = beliefMPC_rd + 35;
beliefMPC_FLOAT beliefMPC_Lbyrd5[7];
beliefMPC_FLOAT* beliefMPC_grad_cost5 = beliefMPC_grad_cost + 35;
beliefMPC_FLOAT* beliefMPC_grad_eq5 = beliefMPC_grad_eq + 35;
beliefMPC_FLOAT* beliefMPC_grad_ineq5 = beliefMPC_grad_ineq + 35;
beliefMPC_FLOAT beliefMPC_ctv5[7];
beliefMPC_FLOAT* beliefMPC_v5 = beliefMPC_v + 30;
beliefMPC_FLOAT beliefMPC_re5[5];
beliefMPC_FLOAT beliefMPC_beta5[5];
beliefMPC_FLOAT beliefMPC_betacc5[5];
beliefMPC_FLOAT* beliefMPC_dvaff5 = beliefMPC_dv_aff + 30;
beliefMPC_FLOAT* beliefMPC_dvcc5 = beliefMPC_dv_cc + 30;
beliefMPC_FLOAT beliefMPC_V5[35];
beliefMPC_FLOAT beliefMPC_Yd5[15];
beliefMPC_FLOAT beliefMPC_Ld5[15];
beliefMPC_FLOAT beliefMPC_yy5[5];
beliefMPC_FLOAT beliefMPC_bmy5[5];
int beliefMPC_lbIdx5[7] = {0, 1, 2, 3, 4, 5, 6};
beliefMPC_FLOAT* beliefMPC_llb5 = beliefMPC_l + 70;
beliefMPC_FLOAT* beliefMPC_slb5 = beliefMPC_s + 70;
beliefMPC_FLOAT* beliefMPC_llbbyslb5 = beliefMPC_lbys + 70;
beliefMPC_FLOAT beliefMPC_rilb5[7];
beliefMPC_FLOAT* beliefMPC_dllbaff5 = beliefMPC_dl_aff + 70;
beliefMPC_FLOAT* beliefMPC_dslbaff5 = beliefMPC_ds_aff + 70;
beliefMPC_FLOAT* beliefMPC_dllbcc5 = beliefMPC_dl_cc + 70;
beliefMPC_FLOAT* beliefMPC_dslbcc5 = beliefMPC_ds_cc + 70;
beliefMPC_FLOAT* beliefMPC_ccrhsl5 = beliefMPC_ccrhs + 70;
int beliefMPC_ubIdx5[7] = {0, 1, 2, 3, 4, 5, 6};
beliefMPC_FLOAT* beliefMPC_lub5 = beliefMPC_l + 77;
beliefMPC_FLOAT* beliefMPC_sub5 = beliefMPC_s + 77;
beliefMPC_FLOAT* beliefMPC_lubbysub5 = beliefMPC_lbys + 77;
beliefMPC_FLOAT beliefMPC_riub5[7];
beliefMPC_FLOAT* beliefMPC_dlubaff5 = beliefMPC_dl_aff + 77;
beliefMPC_FLOAT* beliefMPC_dsubaff5 = beliefMPC_ds_aff + 77;
beliefMPC_FLOAT* beliefMPC_dlubcc5 = beliefMPC_dl_cc + 77;
beliefMPC_FLOAT* beliefMPC_dsubcc5 = beliefMPC_ds_cc + 77;
beliefMPC_FLOAT* beliefMPC_ccrhsub5 = beliefMPC_ccrhs + 77;
beliefMPC_FLOAT beliefMPC_Phi5[7];
beliefMPC_FLOAT beliefMPC_W5[7];
beliefMPC_FLOAT beliefMPC_Ysd5[25];
beliefMPC_FLOAT beliefMPC_Lsd5[25];
beliefMPC_FLOAT beliefMPC_f6[7] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
beliefMPC_FLOAT* beliefMPC_z6 = beliefMPC_z + 42;
beliefMPC_FLOAT* beliefMPC_dzaff6 = beliefMPC_dz_aff + 42;
beliefMPC_FLOAT* beliefMPC_dzcc6 = beliefMPC_dz_cc + 42;
beliefMPC_FLOAT* beliefMPC_rd6 = beliefMPC_rd + 42;
beliefMPC_FLOAT beliefMPC_Lbyrd6[7];
beliefMPC_FLOAT* beliefMPC_grad_cost6 = beliefMPC_grad_cost + 42;
beliefMPC_FLOAT* beliefMPC_grad_eq6 = beliefMPC_grad_eq + 42;
beliefMPC_FLOAT* beliefMPC_grad_ineq6 = beliefMPC_grad_ineq + 42;
beliefMPC_FLOAT beliefMPC_ctv6[7];
beliefMPC_FLOAT* beliefMPC_v6 = beliefMPC_v + 35;
beliefMPC_FLOAT beliefMPC_re6[5];
beliefMPC_FLOAT beliefMPC_beta6[5];
beliefMPC_FLOAT beliefMPC_betacc6[5];
beliefMPC_FLOAT* beliefMPC_dvaff6 = beliefMPC_dv_aff + 35;
beliefMPC_FLOAT* beliefMPC_dvcc6 = beliefMPC_dv_cc + 35;
beliefMPC_FLOAT beliefMPC_V6[35];
beliefMPC_FLOAT beliefMPC_Yd6[15];
beliefMPC_FLOAT beliefMPC_Ld6[15];
beliefMPC_FLOAT beliefMPC_yy6[5];
beliefMPC_FLOAT beliefMPC_bmy6[5];
int beliefMPC_lbIdx6[7] = {0, 1, 2, 3, 4, 5, 6};
beliefMPC_FLOAT* beliefMPC_llb6 = beliefMPC_l + 84;
beliefMPC_FLOAT* beliefMPC_slb6 = beliefMPC_s + 84;
beliefMPC_FLOAT* beliefMPC_llbbyslb6 = beliefMPC_lbys + 84;
beliefMPC_FLOAT beliefMPC_rilb6[7];
beliefMPC_FLOAT* beliefMPC_dllbaff6 = beliefMPC_dl_aff + 84;
beliefMPC_FLOAT* beliefMPC_dslbaff6 = beliefMPC_ds_aff + 84;
beliefMPC_FLOAT* beliefMPC_dllbcc6 = beliefMPC_dl_cc + 84;
beliefMPC_FLOAT* beliefMPC_dslbcc6 = beliefMPC_ds_cc + 84;
beliefMPC_FLOAT* beliefMPC_ccrhsl6 = beliefMPC_ccrhs + 84;
int beliefMPC_ubIdx6[7] = {0, 1, 2, 3, 4, 5, 6};
beliefMPC_FLOAT* beliefMPC_lub6 = beliefMPC_l + 91;
beliefMPC_FLOAT* beliefMPC_sub6 = beliefMPC_s + 91;
beliefMPC_FLOAT* beliefMPC_lubbysub6 = beliefMPC_lbys + 91;
beliefMPC_FLOAT beliefMPC_riub6[7];
beliefMPC_FLOAT* beliefMPC_dlubaff6 = beliefMPC_dl_aff + 91;
beliefMPC_FLOAT* beliefMPC_dsubaff6 = beliefMPC_ds_aff + 91;
beliefMPC_FLOAT* beliefMPC_dlubcc6 = beliefMPC_dl_cc + 91;
beliefMPC_FLOAT* beliefMPC_dsubcc6 = beliefMPC_ds_cc + 91;
beliefMPC_FLOAT* beliefMPC_ccrhsub6 = beliefMPC_ccrhs + 91;
beliefMPC_FLOAT beliefMPC_Phi6[7];
beliefMPC_FLOAT beliefMPC_W6[7];
beliefMPC_FLOAT beliefMPC_Ysd6[25];
beliefMPC_FLOAT beliefMPC_Lsd6[25];
beliefMPC_FLOAT beliefMPC_f7[7] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
beliefMPC_FLOAT* beliefMPC_z7 = beliefMPC_z + 49;
beliefMPC_FLOAT* beliefMPC_dzaff7 = beliefMPC_dz_aff + 49;
beliefMPC_FLOAT* beliefMPC_dzcc7 = beliefMPC_dz_cc + 49;
beliefMPC_FLOAT* beliefMPC_rd7 = beliefMPC_rd + 49;
beliefMPC_FLOAT beliefMPC_Lbyrd7[7];
beliefMPC_FLOAT* beliefMPC_grad_cost7 = beliefMPC_grad_cost + 49;
beliefMPC_FLOAT* beliefMPC_grad_eq7 = beliefMPC_grad_eq + 49;
beliefMPC_FLOAT* beliefMPC_grad_ineq7 = beliefMPC_grad_ineq + 49;
beliefMPC_FLOAT beliefMPC_ctv7[7];
beliefMPC_FLOAT* beliefMPC_v7 = beliefMPC_v + 40;
beliefMPC_FLOAT beliefMPC_re7[5];
beliefMPC_FLOAT beliefMPC_beta7[5];
beliefMPC_FLOAT beliefMPC_betacc7[5];
beliefMPC_FLOAT* beliefMPC_dvaff7 = beliefMPC_dv_aff + 40;
beliefMPC_FLOAT* beliefMPC_dvcc7 = beliefMPC_dv_cc + 40;
beliefMPC_FLOAT beliefMPC_V7[35];
beliefMPC_FLOAT beliefMPC_Yd7[15];
beliefMPC_FLOAT beliefMPC_Ld7[15];
beliefMPC_FLOAT beliefMPC_yy7[5];
beliefMPC_FLOAT beliefMPC_bmy7[5];
int beliefMPC_lbIdx7[7] = {0, 1, 2, 3, 4, 5, 6};
beliefMPC_FLOAT* beliefMPC_llb7 = beliefMPC_l + 98;
beliefMPC_FLOAT* beliefMPC_slb7 = beliefMPC_s + 98;
beliefMPC_FLOAT* beliefMPC_llbbyslb7 = beliefMPC_lbys + 98;
beliefMPC_FLOAT beliefMPC_rilb7[7];
beliefMPC_FLOAT* beliefMPC_dllbaff7 = beliefMPC_dl_aff + 98;
beliefMPC_FLOAT* beliefMPC_dslbaff7 = beliefMPC_ds_aff + 98;
beliefMPC_FLOAT* beliefMPC_dllbcc7 = beliefMPC_dl_cc + 98;
beliefMPC_FLOAT* beliefMPC_dslbcc7 = beliefMPC_ds_cc + 98;
beliefMPC_FLOAT* beliefMPC_ccrhsl7 = beliefMPC_ccrhs + 98;
int beliefMPC_ubIdx7[7] = {0, 1, 2, 3, 4, 5, 6};
beliefMPC_FLOAT* beliefMPC_lub7 = beliefMPC_l + 105;
beliefMPC_FLOAT* beliefMPC_sub7 = beliefMPC_s + 105;
beliefMPC_FLOAT* beliefMPC_lubbysub7 = beliefMPC_lbys + 105;
beliefMPC_FLOAT beliefMPC_riub7[7];
beliefMPC_FLOAT* beliefMPC_dlubaff7 = beliefMPC_dl_aff + 105;
beliefMPC_FLOAT* beliefMPC_dsubaff7 = beliefMPC_ds_aff + 105;
beliefMPC_FLOAT* beliefMPC_dlubcc7 = beliefMPC_dl_cc + 105;
beliefMPC_FLOAT* beliefMPC_dsubcc7 = beliefMPC_ds_cc + 105;
beliefMPC_FLOAT* beliefMPC_ccrhsub7 = beliefMPC_ccrhs + 105;
beliefMPC_FLOAT beliefMPC_Phi7[7];
beliefMPC_FLOAT beliefMPC_W7[7];
beliefMPC_FLOAT beliefMPC_Ysd7[25];
beliefMPC_FLOAT beliefMPC_Lsd7[25];
beliefMPC_FLOAT beliefMPC_f8[7] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
beliefMPC_FLOAT* beliefMPC_z8 = beliefMPC_z + 56;
beliefMPC_FLOAT* beliefMPC_dzaff8 = beliefMPC_dz_aff + 56;
beliefMPC_FLOAT* beliefMPC_dzcc8 = beliefMPC_dz_cc + 56;
beliefMPC_FLOAT* beliefMPC_rd8 = beliefMPC_rd + 56;
beliefMPC_FLOAT beliefMPC_Lbyrd8[7];
beliefMPC_FLOAT* beliefMPC_grad_cost8 = beliefMPC_grad_cost + 56;
beliefMPC_FLOAT* beliefMPC_grad_eq8 = beliefMPC_grad_eq + 56;
beliefMPC_FLOAT* beliefMPC_grad_ineq8 = beliefMPC_grad_ineq + 56;
beliefMPC_FLOAT beliefMPC_ctv8[7];
beliefMPC_FLOAT* beliefMPC_v8 = beliefMPC_v + 45;
beliefMPC_FLOAT beliefMPC_re8[5];
beliefMPC_FLOAT beliefMPC_beta8[5];
beliefMPC_FLOAT beliefMPC_betacc8[5];
beliefMPC_FLOAT* beliefMPC_dvaff8 = beliefMPC_dv_aff + 45;
beliefMPC_FLOAT* beliefMPC_dvcc8 = beliefMPC_dv_cc + 45;
beliefMPC_FLOAT beliefMPC_V8[35];
beliefMPC_FLOAT beliefMPC_Yd8[15];
beliefMPC_FLOAT beliefMPC_Ld8[15];
beliefMPC_FLOAT beliefMPC_yy8[5];
beliefMPC_FLOAT beliefMPC_bmy8[5];
int beliefMPC_lbIdx8[7] = {0, 1, 2, 3, 4, 5, 6};
beliefMPC_FLOAT* beliefMPC_llb8 = beliefMPC_l + 112;
beliefMPC_FLOAT* beliefMPC_slb8 = beliefMPC_s + 112;
beliefMPC_FLOAT* beliefMPC_llbbyslb8 = beliefMPC_lbys + 112;
beliefMPC_FLOAT beliefMPC_rilb8[7];
beliefMPC_FLOAT* beliefMPC_dllbaff8 = beliefMPC_dl_aff + 112;
beliefMPC_FLOAT* beliefMPC_dslbaff8 = beliefMPC_ds_aff + 112;
beliefMPC_FLOAT* beliefMPC_dllbcc8 = beliefMPC_dl_cc + 112;
beliefMPC_FLOAT* beliefMPC_dslbcc8 = beliefMPC_ds_cc + 112;
beliefMPC_FLOAT* beliefMPC_ccrhsl8 = beliefMPC_ccrhs + 112;
int beliefMPC_ubIdx8[7] = {0, 1, 2, 3, 4, 5, 6};
beliefMPC_FLOAT* beliefMPC_lub8 = beliefMPC_l + 119;
beliefMPC_FLOAT* beliefMPC_sub8 = beliefMPC_s + 119;
beliefMPC_FLOAT* beliefMPC_lubbysub8 = beliefMPC_lbys + 119;
beliefMPC_FLOAT beliefMPC_riub8[7];
beliefMPC_FLOAT* beliefMPC_dlubaff8 = beliefMPC_dl_aff + 119;
beliefMPC_FLOAT* beliefMPC_dsubaff8 = beliefMPC_ds_aff + 119;
beliefMPC_FLOAT* beliefMPC_dlubcc8 = beliefMPC_dl_cc + 119;
beliefMPC_FLOAT* beliefMPC_dsubcc8 = beliefMPC_ds_cc + 119;
beliefMPC_FLOAT* beliefMPC_ccrhsub8 = beliefMPC_ccrhs + 119;
beliefMPC_FLOAT beliefMPC_Phi8[7];
beliefMPC_FLOAT beliefMPC_W8[7];
beliefMPC_FLOAT beliefMPC_Ysd8[25];
beliefMPC_FLOAT beliefMPC_Lsd8[25];
beliefMPC_FLOAT beliefMPC_H9[5] = {0.0000000000000000E+000, 0.0000000000000000E+000, 2.0000000000000000E+001, 2.0000000000000000E+001, 2.0000000000000000E+001};
beliefMPC_FLOAT beliefMPC_f9[5] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
beliefMPC_FLOAT* beliefMPC_z9 = beliefMPC_z + 63;
beliefMPC_FLOAT* beliefMPC_dzaff9 = beliefMPC_dz_aff + 63;
beliefMPC_FLOAT* beliefMPC_dzcc9 = beliefMPC_dz_cc + 63;
beliefMPC_FLOAT* beliefMPC_rd9 = beliefMPC_rd + 63;
beliefMPC_FLOAT beliefMPC_Lbyrd9[5];
beliefMPC_FLOAT* beliefMPC_grad_cost9 = beliefMPC_grad_cost + 63;
beliefMPC_FLOAT* beliefMPC_grad_eq9 = beliefMPC_grad_eq + 63;
beliefMPC_FLOAT* beliefMPC_grad_ineq9 = beliefMPC_grad_ineq + 63;
beliefMPC_FLOAT beliefMPC_ctv9[5];
int beliefMPC_lbIdx9[5] = {0, 1, 2, 3, 4};
beliefMPC_FLOAT* beliefMPC_llb9 = beliefMPC_l + 126;
beliefMPC_FLOAT* beliefMPC_slb9 = beliefMPC_s + 126;
beliefMPC_FLOAT* beliefMPC_llbbyslb9 = beliefMPC_lbys + 126;
beliefMPC_FLOAT beliefMPC_rilb9[5];
beliefMPC_FLOAT* beliefMPC_dllbaff9 = beliefMPC_dl_aff + 126;
beliefMPC_FLOAT* beliefMPC_dslbaff9 = beliefMPC_ds_aff + 126;
beliefMPC_FLOAT* beliefMPC_dllbcc9 = beliefMPC_dl_cc + 126;
beliefMPC_FLOAT* beliefMPC_dslbcc9 = beliefMPC_ds_cc + 126;
beliefMPC_FLOAT* beliefMPC_ccrhsl9 = beliefMPC_ccrhs + 126;
int beliefMPC_ubIdx9[5] = {0, 1, 2, 3, 4};
beliefMPC_FLOAT* beliefMPC_lub9 = beliefMPC_l + 131;
beliefMPC_FLOAT* beliefMPC_sub9 = beliefMPC_s + 131;
beliefMPC_FLOAT* beliefMPC_lubbysub9 = beliefMPC_lbys + 131;
beliefMPC_FLOAT beliefMPC_riub9[5];
beliefMPC_FLOAT* beliefMPC_dlubaff9 = beliefMPC_dl_aff + 131;
beliefMPC_FLOAT* beliefMPC_dsubaff9 = beliefMPC_ds_aff + 131;
beliefMPC_FLOAT* beliefMPC_dlubcc9 = beliefMPC_dl_cc + 131;
beliefMPC_FLOAT* beliefMPC_dsubcc9 = beliefMPC_ds_cc + 131;
beliefMPC_FLOAT* beliefMPC_ccrhsub9 = beliefMPC_ccrhs + 131;
beliefMPC_FLOAT beliefMPC_Phi9[5];
beliefMPC_FLOAT beliefMPC_D9[5] = {-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000};
beliefMPC_FLOAT beliefMPC_W9[5];
beliefMPC_FLOAT musigma;
beliefMPC_FLOAT sigma_3rdroot;
beliefMPC_FLOAT beliefMPC_Diag1_0[7];
beliefMPC_FLOAT beliefMPC_Diag2_0[7];
beliefMPC_FLOAT beliefMPC_L_0[21];




/* SOLVER CODE --------------------------------------------------------- */
int beliefMPC_solve(beliefMPC_params* params, beliefMPC_output* output, beliefMPC_info* info)
{	
int exitcode;

#if beliefMPC_SET_TIMING == 1
	beliefMPC_timer solvertimer;
	beliefMPC_tic(&solvertimer);
#endif
/* FUNCTION CALLS INTO LA LIBRARY -------------------------------------- */
info->it = 0;
beliefMPC_LA_INITIALIZEVECTOR_68(beliefMPC_z, 0);
beliefMPC_LA_INITIALIZEVECTOR_50(beliefMPC_v, 1);
beliefMPC_LA_INITIALIZEVECTOR_136(beliefMPC_l, 1);
beliefMPC_LA_INITIALIZEVECTOR_136(beliefMPC_s, 1);
info->mu = 0;
beliefMPC_LA_DOTACC_136(beliefMPC_l, beliefMPC_s, &info->mu);
info->mu /= 136;
while( 1 ){
info->pobj = 0;
beliefMPC_LA_DIAG_QUADFCN_7(beliefMPC_H0, beliefMPC_f0, beliefMPC_z0, beliefMPC_grad_cost0, &info->pobj);
beliefMPC_LA_DIAG_QUADFCN_7(beliefMPC_H0, beliefMPC_f1, beliefMPC_z1, beliefMPC_grad_cost1, &info->pobj);
beliefMPC_LA_DIAG_QUADFCN_7(beliefMPC_H0, beliefMPC_f2, beliefMPC_z2, beliefMPC_grad_cost2, &info->pobj);
beliefMPC_LA_DIAG_QUADFCN_7(beliefMPC_H0, beliefMPC_f3, beliefMPC_z3, beliefMPC_grad_cost3, &info->pobj);
beliefMPC_LA_DIAG_QUADFCN_7(beliefMPC_H0, beliefMPC_f4, beliefMPC_z4, beliefMPC_grad_cost4, &info->pobj);
beliefMPC_LA_DIAG_QUADFCN_7(beliefMPC_H0, beliefMPC_f5, beliefMPC_z5, beliefMPC_grad_cost5, &info->pobj);
beliefMPC_LA_DIAG_QUADFCN_7(beliefMPC_H0, beliefMPC_f6, beliefMPC_z6, beliefMPC_grad_cost6, &info->pobj);
beliefMPC_LA_DIAG_QUADFCN_7(beliefMPC_H0, beliefMPC_f7, beliefMPC_z7, beliefMPC_grad_cost7, &info->pobj);
beliefMPC_LA_DIAG_QUADFCN_7(beliefMPC_H0, beliefMPC_f8, beliefMPC_z8, beliefMPC_grad_cost8, &info->pobj);
beliefMPC_LA_DIAG_QUADFCN_5(beliefMPC_H9, beliefMPC_f9, beliefMPC_z9, beliefMPC_grad_cost9, &info->pobj);
info->res_eq = 0;
info->dgap = 0;
beliefMPC_LA_DENSE_MVMSUB3_10_7_7(params->C1, beliefMPC_z0, beliefMPC_D1, beliefMPC_z1, params->e1, beliefMPC_v0, beliefMPC_re0, &info->dgap, &info->res_eq);
beliefMPC_LA_DENSE_DIAGZERO_MVMSUB3_5_7_7(params->C2, beliefMPC_z1, beliefMPC_D2, beliefMPC_z2, params->e2, beliefMPC_v1, beliefMPC_re1, &info->dgap, &info->res_eq);
beliefMPC_LA_DENSE_DIAGZERO_MVMSUB3_5_7_7(params->C3, beliefMPC_z2, beliefMPC_D2, beliefMPC_z3, params->e3, beliefMPC_v2, beliefMPC_re2, &info->dgap, &info->res_eq);
beliefMPC_LA_DENSE_DIAGZERO_MVMSUB3_5_7_7(params->C4, beliefMPC_z3, beliefMPC_D2, beliefMPC_z4, params->e4, beliefMPC_v3, beliefMPC_re3, &info->dgap, &info->res_eq);
beliefMPC_LA_DENSE_DIAGZERO_MVMSUB3_5_7_7(params->C5, beliefMPC_z4, beliefMPC_D2, beliefMPC_z5, params->e5, beliefMPC_v4, beliefMPC_re4, &info->dgap, &info->res_eq);
beliefMPC_LA_DENSE_DIAGZERO_MVMSUB3_5_7_7(params->C6, beliefMPC_z5, beliefMPC_D2, beliefMPC_z6, params->e6, beliefMPC_v5, beliefMPC_re5, &info->dgap, &info->res_eq);
beliefMPC_LA_DENSE_DIAGZERO_MVMSUB3_5_7_7(params->C7, beliefMPC_z6, beliefMPC_D2, beliefMPC_z7, params->e7, beliefMPC_v6, beliefMPC_re6, &info->dgap, &info->res_eq);
beliefMPC_LA_DENSE_DIAGZERO_MVMSUB3_5_7_7(params->C8, beliefMPC_z7, beliefMPC_D2, beliefMPC_z8, params->e8, beliefMPC_v7, beliefMPC_re7, &info->dgap, &info->res_eq);
beliefMPC_LA_DENSE_DIAGZERO_MVMSUB3_5_7_5(params->C9, beliefMPC_z8, beliefMPC_D9, beliefMPC_z9, params->e9, beliefMPC_v8, beliefMPC_re8, &info->dgap, &info->res_eq);
beliefMPC_LA_DENSE_MTVM_10_7(params->C1, beliefMPC_v0, beliefMPC_grad_eq0);
beliefMPC_LA_DENSE_MTVM2_5_7_10(params->C2, beliefMPC_v1, beliefMPC_D1, beliefMPC_v0, beliefMPC_grad_eq1);
beliefMPC_LA_DENSE_DIAGZERO_MTVM2_5_7_5(params->C3, beliefMPC_v2, beliefMPC_D2, beliefMPC_v1, beliefMPC_grad_eq2);
beliefMPC_LA_DENSE_DIAGZERO_MTVM2_5_7_5(params->C4, beliefMPC_v3, beliefMPC_D2, beliefMPC_v2, beliefMPC_grad_eq3);
beliefMPC_LA_DENSE_DIAGZERO_MTVM2_5_7_5(params->C5, beliefMPC_v4, beliefMPC_D2, beliefMPC_v3, beliefMPC_grad_eq4);
beliefMPC_LA_DENSE_DIAGZERO_MTVM2_5_7_5(params->C6, beliefMPC_v5, beliefMPC_D2, beliefMPC_v4, beliefMPC_grad_eq5);
beliefMPC_LA_DENSE_DIAGZERO_MTVM2_5_7_5(params->C7, beliefMPC_v6, beliefMPC_D2, beliefMPC_v5, beliefMPC_grad_eq6);
beliefMPC_LA_DENSE_DIAGZERO_MTVM2_5_7_5(params->C8, beliefMPC_v7, beliefMPC_D2, beliefMPC_v6, beliefMPC_grad_eq7);
beliefMPC_LA_DENSE_DIAGZERO_MTVM2_5_7_5(params->C9, beliefMPC_v8, beliefMPC_D2, beliefMPC_v7, beliefMPC_grad_eq8);
beliefMPC_LA_DIAGZERO_MTVM_5_5(beliefMPC_D9, beliefMPC_v8, beliefMPC_grad_eq9);
info->res_ineq = 0;
beliefMPC_LA_VSUBADD3_7(params->lb1, beliefMPC_z0, beliefMPC_lbIdx0, beliefMPC_llb0, beliefMPC_slb0, beliefMPC_rilb0, &info->dgap, &info->res_ineq);
beliefMPC_LA_VSUBADD2_7(beliefMPC_z0, beliefMPC_ubIdx0, params->ub1, beliefMPC_lub0, beliefMPC_sub0, beliefMPC_riub0, &info->dgap, &info->res_ineq);
beliefMPC_LA_VSUBADD3_7(params->lb2, beliefMPC_z1, beliefMPC_lbIdx1, beliefMPC_llb1, beliefMPC_slb1, beliefMPC_rilb1, &info->dgap, &info->res_ineq);
beliefMPC_LA_VSUBADD2_7(beliefMPC_z1, beliefMPC_ubIdx1, params->ub2, beliefMPC_lub1, beliefMPC_sub1, beliefMPC_riub1, &info->dgap, &info->res_ineq);
beliefMPC_LA_VSUBADD3_7(params->lb3, beliefMPC_z2, beliefMPC_lbIdx2, beliefMPC_llb2, beliefMPC_slb2, beliefMPC_rilb2, &info->dgap, &info->res_ineq);
beliefMPC_LA_VSUBADD2_7(beliefMPC_z2, beliefMPC_ubIdx2, params->ub3, beliefMPC_lub2, beliefMPC_sub2, beliefMPC_riub2, &info->dgap, &info->res_ineq);
beliefMPC_LA_VSUBADD3_7(params->lb4, beliefMPC_z3, beliefMPC_lbIdx3, beliefMPC_llb3, beliefMPC_slb3, beliefMPC_rilb3, &info->dgap, &info->res_ineq);
beliefMPC_LA_VSUBADD2_7(beliefMPC_z3, beliefMPC_ubIdx3, params->ub4, beliefMPC_lub3, beliefMPC_sub3, beliefMPC_riub3, &info->dgap, &info->res_ineq);
beliefMPC_LA_VSUBADD3_7(params->lb5, beliefMPC_z4, beliefMPC_lbIdx4, beliefMPC_llb4, beliefMPC_slb4, beliefMPC_rilb4, &info->dgap, &info->res_ineq);
beliefMPC_LA_VSUBADD2_7(beliefMPC_z4, beliefMPC_ubIdx4, params->ub5, beliefMPC_lub4, beliefMPC_sub4, beliefMPC_riub4, &info->dgap, &info->res_ineq);
beliefMPC_LA_VSUBADD3_7(params->lb6, beliefMPC_z5, beliefMPC_lbIdx5, beliefMPC_llb5, beliefMPC_slb5, beliefMPC_rilb5, &info->dgap, &info->res_ineq);
beliefMPC_LA_VSUBADD2_7(beliefMPC_z5, beliefMPC_ubIdx5, params->ub6, beliefMPC_lub5, beliefMPC_sub5, beliefMPC_riub5, &info->dgap, &info->res_ineq);
beliefMPC_LA_VSUBADD3_7(params->lb7, beliefMPC_z6, beliefMPC_lbIdx6, beliefMPC_llb6, beliefMPC_slb6, beliefMPC_rilb6, &info->dgap, &info->res_ineq);
beliefMPC_LA_VSUBADD2_7(beliefMPC_z6, beliefMPC_ubIdx6, params->ub7, beliefMPC_lub6, beliefMPC_sub6, beliefMPC_riub6, &info->dgap, &info->res_ineq);
beliefMPC_LA_VSUBADD3_7(params->lb8, beliefMPC_z7, beliefMPC_lbIdx7, beliefMPC_llb7, beliefMPC_slb7, beliefMPC_rilb7, &info->dgap, &info->res_ineq);
beliefMPC_LA_VSUBADD2_7(beliefMPC_z7, beliefMPC_ubIdx7, params->ub8, beliefMPC_lub7, beliefMPC_sub7, beliefMPC_riub7, &info->dgap, &info->res_ineq);
beliefMPC_LA_VSUBADD3_7(params->lb9, beliefMPC_z8, beliefMPC_lbIdx8, beliefMPC_llb8, beliefMPC_slb8, beliefMPC_rilb8, &info->dgap, &info->res_ineq);
beliefMPC_LA_VSUBADD2_7(beliefMPC_z8, beliefMPC_ubIdx8, params->ub9, beliefMPC_lub8, beliefMPC_sub8, beliefMPC_riub8, &info->dgap, &info->res_ineq);
beliefMPC_LA_VSUBADD3_5(params->lb10, beliefMPC_z9, beliefMPC_lbIdx9, beliefMPC_llb9, beliefMPC_slb9, beliefMPC_rilb9, &info->dgap, &info->res_ineq);
beliefMPC_LA_VSUBADD2_5(beliefMPC_z9, beliefMPC_ubIdx9, params->ub10, beliefMPC_lub9, beliefMPC_sub9, beliefMPC_riub9, &info->dgap, &info->res_ineq);
beliefMPC_LA_INEQ_B_GRAD_7_7_7(beliefMPC_lub0, beliefMPC_sub0, beliefMPC_riub0, beliefMPC_llb0, beliefMPC_slb0, beliefMPC_rilb0, beliefMPC_lbIdx0, beliefMPC_ubIdx0, beliefMPC_grad_ineq0, beliefMPC_lubbysub0, beliefMPC_llbbyslb0);
beliefMPC_LA_INEQ_B_GRAD_7_7_7(beliefMPC_lub1, beliefMPC_sub1, beliefMPC_riub1, beliefMPC_llb1, beliefMPC_slb1, beliefMPC_rilb1, beliefMPC_lbIdx1, beliefMPC_ubIdx1, beliefMPC_grad_ineq1, beliefMPC_lubbysub1, beliefMPC_llbbyslb1);
beliefMPC_LA_INEQ_B_GRAD_7_7_7(beliefMPC_lub2, beliefMPC_sub2, beliefMPC_riub2, beliefMPC_llb2, beliefMPC_slb2, beliefMPC_rilb2, beliefMPC_lbIdx2, beliefMPC_ubIdx2, beliefMPC_grad_ineq2, beliefMPC_lubbysub2, beliefMPC_llbbyslb2);
beliefMPC_LA_INEQ_B_GRAD_7_7_7(beliefMPC_lub3, beliefMPC_sub3, beliefMPC_riub3, beliefMPC_llb3, beliefMPC_slb3, beliefMPC_rilb3, beliefMPC_lbIdx3, beliefMPC_ubIdx3, beliefMPC_grad_ineq3, beliefMPC_lubbysub3, beliefMPC_llbbyslb3);
beliefMPC_LA_INEQ_B_GRAD_7_7_7(beliefMPC_lub4, beliefMPC_sub4, beliefMPC_riub4, beliefMPC_llb4, beliefMPC_slb4, beliefMPC_rilb4, beliefMPC_lbIdx4, beliefMPC_ubIdx4, beliefMPC_grad_ineq4, beliefMPC_lubbysub4, beliefMPC_llbbyslb4);
beliefMPC_LA_INEQ_B_GRAD_7_7_7(beliefMPC_lub5, beliefMPC_sub5, beliefMPC_riub5, beliefMPC_llb5, beliefMPC_slb5, beliefMPC_rilb5, beliefMPC_lbIdx5, beliefMPC_ubIdx5, beliefMPC_grad_ineq5, beliefMPC_lubbysub5, beliefMPC_llbbyslb5);
beliefMPC_LA_INEQ_B_GRAD_7_7_7(beliefMPC_lub6, beliefMPC_sub6, beliefMPC_riub6, beliefMPC_llb6, beliefMPC_slb6, beliefMPC_rilb6, beliefMPC_lbIdx6, beliefMPC_ubIdx6, beliefMPC_grad_ineq6, beliefMPC_lubbysub6, beliefMPC_llbbyslb6);
beliefMPC_LA_INEQ_B_GRAD_7_7_7(beliefMPC_lub7, beliefMPC_sub7, beliefMPC_riub7, beliefMPC_llb7, beliefMPC_slb7, beliefMPC_rilb7, beliefMPC_lbIdx7, beliefMPC_ubIdx7, beliefMPC_grad_ineq7, beliefMPC_lubbysub7, beliefMPC_llbbyslb7);
beliefMPC_LA_INEQ_B_GRAD_7_7_7(beliefMPC_lub8, beliefMPC_sub8, beliefMPC_riub8, beliefMPC_llb8, beliefMPC_slb8, beliefMPC_rilb8, beliefMPC_lbIdx8, beliefMPC_ubIdx8, beliefMPC_grad_ineq8, beliefMPC_lubbysub8, beliefMPC_llbbyslb8);
beliefMPC_LA_INEQ_B_GRAD_5_5_5(beliefMPC_lub9, beliefMPC_sub9, beliefMPC_riub9, beliefMPC_llb9, beliefMPC_slb9, beliefMPC_rilb9, beliefMPC_lbIdx9, beliefMPC_ubIdx9, beliefMPC_grad_ineq9, beliefMPC_lubbysub9, beliefMPC_llbbyslb9);
info->dobj = info->pobj - info->dgap;
info->rdgap = info->pobj ? info->dgap / info->pobj : 1e6;
if( info->rdgap < 0 ) info->rdgap = -info->rdgap;
if( info->mu < beliefMPC_SET_ACC_KKTCOMPL
    && (info->rdgap < beliefMPC_SET_ACC_RDGAP || info->dgap < beliefMPC_SET_ACC_KKTCOMPL)
    && info->res_eq < beliefMPC_SET_ACC_RESEQ
    && info->res_ineq < beliefMPC_SET_ACC_RESINEQ ){
exitcode = beliefMPC_OPTIMAL; break; }
if( info->it == beliefMPC_SET_MAXIT ){
exitcode = beliefMPC_MAXITREACHED; break; }
beliefMPC_LA_VVADD3_68(beliefMPC_grad_cost, beliefMPC_grad_eq, beliefMPC_grad_ineq, beliefMPC_rd);
beliefMPC_LA_DIAG_CHOL_ONELOOP_LBUB_7_7_7(beliefMPC_H0, beliefMPC_llbbyslb0, beliefMPC_lbIdx0, beliefMPC_lubbysub0, beliefMPC_ubIdx0, beliefMPC_Phi0);
beliefMPC_LA_DIAG_MATRIXFORWARDSUB_10_7(beliefMPC_Phi0, params->C1, beliefMPC_V0);
beliefMPC_LA_DIAG_FORWARDSUB_7(beliefMPC_Phi0, beliefMPC_rd0, beliefMPC_Lbyrd0);
beliefMPC_LA_DIAG_CHOL_ONELOOP_LBUB_7_7_7(beliefMPC_H0, beliefMPC_llbbyslb1, beliefMPC_lbIdx1, beliefMPC_lubbysub1, beliefMPC_ubIdx1, beliefMPC_Phi1);
beliefMPC_LA_DIAG_MATRIXFORWARDSUB_5_7(beliefMPC_Phi1, params->C2, beliefMPC_V1);
beliefMPC_LA_DIAG_MATRIXFORWARDSUB_10_7(beliefMPC_Phi1, beliefMPC_D1, beliefMPC_W1);
beliefMPC_LA_DENSE_MMTM_10_7_5(beliefMPC_W1, beliefMPC_V1, beliefMPC_Ysd1);
beliefMPC_LA_DIAG_FORWARDSUB_7(beliefMPC_Phi1, beliefMPC_rd1, beliefMPC_Lbyrd1);
beliefMPC_LA_DIAG_CHOL_ONELOOP_LBUB_7_7_7(beliefMPC_H0, beliefMPC_llbbyslb2, beliefMPC_lbIdx2, beliefMPC_lubbysub2, beliefMPC_ubIdx2, beliefMPC_Phi2);
beliefMPC_LA_DIAG_MATRIXFORWARDSUB_5_7(beliefMPC_Phi2, params->C3, beliefMPC_V2);
beliefMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_5_7(beliefMPC_Phi2, beliefMPC_D2, beliefMPC_W2);
beliefMPC_LA_DENSE_DIAGZERO_MMTM_5_7_5(beliefMPC_W2, beliefMPC_V2, beliefMPC_Ysd2);
beliefMPC_LA_DIAG_FORWARDSUB_7(beliefMPC_Phi2, beliefMPC_rd2, beliefMPC_Lbyrd2);
beliefMPC_LA_DIAG_CHOL_ONELOOP_LBUB_7_7_7(beliefMPC_H0, beliefMPC_llbbyslb3, beliefMPC_lbIdx3, beliefMPC_lubbysub3, beliefMPC_ubIdx3, beliefMPC_Phi3);
beliefMPC_LA_DIAG_MATRIXFORWARDSUB_5_7(beliefMPC_Phi3, params->C4, beliefMPC_V3);
beliefMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_5_7(beliefMPC_Phi3, beliefMPC_D2, beliefMPC_W3);
beliefMPC_LA_DENSE_DIAGZERO_MMTM_5_7_5(beliefMPC_W3, beliefMPC_V3, beliefMPC_Ysd3);
beliefMPC_LA_DIAG_FORWARDSUB_7(beliefMPC_Phi3, beliefMPC_rd3, beliefMPC_Lbyrd3);
beliefMPC_LA_DIAG_CHOL_ONELOOP_LBUB_7_7_7(beliefMPC_H0, beliefMPC_llbbyslb4, beliefMPC_lbIdx4, beliefMPC_lubbysub4, beliefMPC_ubIdx4, beliefMPC_Phi4);
beliefMPC_LA_DIAG_MATRIXFORWARDSUB_5_7(beliefMPC_Phi4, params->C5, beliefMPC_V4);
beliefMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_5_7(beliefMPC_Phi4, beliefMPC_D2, beliefMPC_W4);
beliefMPC_LA_DENSE_DIAGZERO_MMTM_5_7_5(beliefMPC_W4, beliefMPC_V4, beliefMPC_Ysd4);
beliefMPC_LA_DIAG_FORWARDSUB_7(beliefMPC_Phi4, beliefMPC_rd4, beliefMPC_Lbyrd4);
beliefMPC_LA_DIAG_CHOL_ONELOOP_LBUB_7_7_7(beliefMPC_H0, beliefMPC_llbbyslb5, beliefMPC_lbIdx5, beliefMPC_lubbysub5, beliefMPC_ubIdx5, beliefMPC_Phi5);
beliefMPC_LA_DIAG_MATRIXFORWARDSUB_5_7(beliefMPC_Phi5, params->C6, beliefMPC_V5);
beliefMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_5_7(beliefMPC_Phi5, beliefMPC_D2, beliefMPC_W5);
beliefMPC_LA_DENSE_DIAGZERO_MMTM_5_7_5(beliefMPC_W5, beliefMPC_V5, beliefMPC_Ysd5);
beliefMPC_LA_DIAG_FORWARDSUB_7(beliefMPC_Phi5, beliefMPC_rd5, beliefMPC_Lbyrd5);
beliefMPC_LA_DIAG_CHOL_ONELOOP_LBUB_7_7_7(beliefMPC_H0, beliefMPC_llbbyslb6, beliefMPC_lbIdx6, beliefMPC_lubbysub6, beliefMPC_ubIdx6, beliefMPC_Phi6);
beliefMPC_LA_DIAG_MATRIXFORWARDSUB_5_7(beliefMPC_Phi6, params->C7, beliefMPC_V6);
beliefMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_5_7(beliefMPC_Phi6, beliefMPC_D2, beliefMPC_W6);
beliefMPC_LA_DENSE_DIAGZERO_MMTM_5_7_5(beliefMPC_W6, beliefMPC_V6, beliefMPC_Ysd6);
beliefMPC_LA_DIAG_FORWARDSUB_7(beliefMPC_Phi6, beliefMPC_rd6, beliefMPC_Lbyrd6);
beliefMPC_LA_DIAG_CHOL_ONELOOP_LBUB_7_7_7(beliefMPC_H0, beliefMPC_llbbyslb7, beliefMPC_lbIdx7, beliefMPC_lubbysub7, beliefMPC_ubIdx7, beliefMPC_Phi7);
beliefMPC_LA_DIAG_MATRIXFORWARDSUB_5_7(beliefMPC_Phi7, params->C8, beliefMPC_V7);
beliefMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_5_7(beliefMPC_Phi7, beliefMPC_D2, beliefMPC_W7);
beliefMPC_LA_DENSE_DIAGZERO_MMTM_5_7_5(beliefMPC_W7, beliefMPC_V7, beliefMPC_Ysd7);
beliefMPC_LA_DIAG_FORWARDSUB_7(beliefMPC_Phi7, beliefMPC_rd7, beliefMPC_Lbyrd7);
beliefMPC_LA_DIAG_CHOL_ONELOOP_LBUB_7_7_7(beliefMPC_H0, beliefMPC_llbbyslb8, beliefMPC_lbIdx8, beliefMPC_lubbysub8, beliefMPC_ubIdx8, beliefMPC_Phi8);
beliefMPC_LA_DIAG_MATRIXFORWARDSUB_5_7(beliefMPC_Phi8, params->C9, beliefMPC_V8);
beliefMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_5_7(beliefMPC_Phi8, beliefMPC_D2, beliefMPC_W8);
beliefMPC_LA_DENSE_DIAGZERO_MMTM_5_7_5(beliefMPC_W8, beliefMPC_V8, beliefMPC_Ysd8);
beliefMPC_LA_DIAG_FORWARDSUB_7(beliefMPC_Phi8, beliefMPC_rd8, beliefMPC_Lbyrd8);
beliefMPC_LA_DIAG_CHOL_ONELOOP_LBUB_5_5_5(beliefMPC_H9, beliefMPC_llbbyslb9, beliefMPC_lbIdx9, beliefMPC_lubbysub9, beliefMPC_ubIdx9, beliefMPC_Phi9);
beliefMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_5_5(beliefMPC_Phi9, beliefMPC_D9, beliefMPC_W9);
beliefMPC_LA_DIAG_FORWARDSUB_5(beliefMPC_Phi9, beliefMPC_rd9, beliefMPC_Lbyrd9);
beliefMPC_LA_DENSE_MMT2_10_7_7(beliefMPC_V0, beliefMPC_W1, beliefMPC_Yd0);
beliefMPC_LA_DENSE_MVMSUB2_10_7_7(beliefMPC_V0, beliefMPC_Lbyrd0, beliefMPC_W1, beliefMPC_Lbyrd1, beliefMPC_re0, beliefMPC_beta0);
beliefMPC_LA_DENSE_DIAGZERO_MMT2_5_7_7(beliefMPC_V1, beliefMPC_W2, beliefMPC_Yd1);
beliefMPC_LA_DENSE_DIAGZERO_2MVMSUB2_5_7_7(beliefMPC_V1, beliefMPC_Lbyrd1, beliefMPC_W2, beliefMPC_Lbyrd2, beliefMPC_re1, beliefMPC_beta1);
beliefMPC_LA_DENSE_DIAGZERO_MMT2_5_7_7(beliefMPC_V2, beliefMPC_W3, beliefMPC_Yd2);
beliefMPC_LA_DENSE_DIAGZERO_2MVMSUB2_5_7_7(beliefMPC_V2, beliefMPC_Lbyrd2, beliefMPC_W3, beliefMPC_Lbyrd3, beliefMPC_re2, beliefMPC_beta2);
beliefMPC_LA_DENSE_DIAGZERO_MMT2_5_7_7(beliefMPC_V3, beliefMPC_W4, beliefMPC_Yd3);
beliefMPC_LA_DENSE_DIAGZERO_2MVMSUB2_5_7_7(beliefMPC_V3, beliefMPC_Lbyrd3, beliefMPC_W4, beliefMPC_Lbyrd4, beliefMPC_re3, beliefMPC_beta3);
beliefMPC_LA_DENSE_DIAGZERO_MMT2_5_7_7(beliefMPC_V4, beliefMPC_W5, beliefMPC_Yd4);
beliefMPC_LA_DENSE_DIAGZERO_2MVMSUB2_5_7_7(beliefMPC_V4, beliefMPC_Lbyrd4, beliefMPC_W5, beliefMPC_Lbyrd5, beliefMPC_re4, beliefMPC_beta4);
beliefMPC_LA_DENSE_DIAGZERO_MMT2_5_7_7(beliefMPC_V5, beliefMPC_W6, beliefMPC_Yd5);
beliefMPC_LA_DENSE_DIAGZERO_2MVMSUB2_5_7_7(beliefMPC_V5, beliefMPC_Lbyrd5, beliefMPC_W6, beliefMPC_Lbyrd6, beliefMPC_re5, beliefMPC_beta5);
beliefMPC_LA_DENSE_DIAGZERO_MMT2_5_7_7(beliefMPC_V6, beliefMPC_W7, beliefMPC_Yd6);
beliefMPC_LA_DENSE_DIAGZERO_2MVMSUB2_5_7_7(beliefMPC_V6, beliefMPC_Lbyrd6, beliefMPC_W7, beliefMPC_Lbyrd7, beliefMPC_re6, beliefMPC_beta6);
beliefMPC_LA_DENSE_DIAGZERO_MMT2_5_7_7(beliefMPC_V7, beliefMPC_W8, beliefMPC_Yd7);
beliefMPC_LA_DENSE_DIAGZERO_2MVMSUB2_5_7_7(beliefMPC_V7, beliefMPC_Lbyrd7, beliefMPC_W8, beliefMPC_Lbyrd8, beliefMPC_re7, beliefMPC_beta7);
beliefMPC_LA_DENSE_DIAGZERO_MMT2_5_7_5(beliefMPC_V8, beliefMPC_W9, beliefMPC_Yd8);
beliefMPC_LA_DENSE_DIAGZERO_2MVMSUB2_5_7_5(beliefMPC_V8, beliefMPC_Lbyrd8, beliefMPC_W9, beliefMPC_Lbyrd9, beliefMPC_re8, beliefMPC_beta8);
beliefMPC_LA_DENSE_CHOL_10(beliefMPC_Yd0, beliefMPC_Ld0);
beliefMPC_LA_DENSE_FORWARDSUB_10(beliefMPC_Ld0, beliefMPC_beta0, beliefMPC_yy0);
beliefMPC_LA_DENSE_MATRIXTFORWARDSUB_5_10(beliefMPC_Ld0, beliefMPC_Ysd1, beliefMPC_Lsd1);
beliefMPC_LA_DENSE_MMTSUB_5_10(beliefMPC_Lsd1, beliefMPC_Yd1);
beliefMPC_LA_DENSE_CHOL_5(beliefMPC_Yd1, beliefMPC_Ld1);
beliefMPC_LA_DENSE_MVMSUB1_5_10(beliefMPC_Lsd1, beliefMPC_yy0, beliefMPC_beta1, beliefMPC_bmy1);
beliefMPC_LA_DENSE_FORWARDSUB_5(beliefMPC_Ld1, beliefMPC_bmy1, beliefMPC_yy1);
beliefMPC_LA_DENSE_MATRIXTFORWARDSUB_5_5(beliefMPC_Ld1, beliefMPC_Ysd2, beliefMPC_Lsd2);
beliefMPC_LA_DENSE_MMTSUB_5_5(beliefMPC_Lsd2, beliefMPC_Yd2);
beliefMPC_LA_DENSE_CHOL_5(beliefMPC_Yd2, beliefMPC_Ld2);
beliefMPC_LA_DENSE_MVMSUB1_5_5(beliefMPC_Lsd2, beliefMPC_yy1, beliefMPC_beta2, beliefMPC_bmy2);
beliefMPC_LA_DENSE_FORWARDSUB_5(beliefMPC_Ld2, beliefMPC_bmy2, beliefMPC_yy2);
beliefMPC_LA_DENSE_MATRIXTFORWARDSUB_5_5(beliefMPC_Ld2, beliefMPC_Ysd3, beliefMPC_Lsd3);
beliefMPC_LA_DENSE_MMTSUB_5_5(beliefMPC_Lsd3, beliefMPC_Yd3);
beliefMPC_LA_DENSE_CHOL_5(beliefMPC_Yd3, beliefMPC_Ld3);
beliefMPC_LA_DENSE_MVMSUB1_5_5(beliefMPC_Lsd3, beliefMPC_yy2, beliefMPC_beta3, beliefMPC_bmy3);
beliefMPC_LA_DENSE_FORWARDSUB_5(beliefMPC_Ld3, beliefMPC_bmy3, beliefMPC_yy3);
beliefMPC_LA_DENSE_MATRIXTFORWARDSUB_5_5(beliefMPC_Ld3, beliefMPC_Ysd4, beliefMPC_Lsd4);
beliefMPC_LA_DENSE_MMTSUB_5_5(beliefMPC_Lsd4, beliefMPC_Yd4);
beliefMPC_LA_DENSE_CHOL_5(beliefMPC_Yd4, beliefMPC_Ld4);
beliefMPC_LA_DENSE_MVMSUB1_5_5(beliefMPC_Lsd4, beliefMPC_yy3, beliefMPC_beta4, beliefMPC_bmy4);
beliefMPC_LA_DENSE_FORWARDSUB_5(beliefMPC_Ld4, beliefMPC_bmy4, beliefMPC_yy4);
beliefMPC_LA_DENSE_MATRIXTFORWARDSUB_5_5(beliefMPC_Ld4, beliefMPC_Ysd5, beliefMPC_Lsd5);
beliefMPC_LA_DENSE_MMTSUB_5_5(beliefMPC_Lsd5, beliefMPC_Yd5);
beliefMPC_LA_DENSE_CHOL_5(beliefMPC_Yd5, beliefMPC_Ld5);
beliefMPC_LA_DENSE_MVMSUB1_5_5(beliefMPC_Lsd5, beliefMPC_yy4, beliefMPC_beta5, beliefMPC_bmy5);
beliefMPC_LA_DENSE_FORWARDSUB_5(beliefMPC_Ld5, beliefMPC_bmy5, beliefMPC_yy5);
beliefMPC_LA_DENSE_MATRIXTFORWARDSUB_5_5(beliefMPC_Ld5, beliefMPC_Ysd6, beliefMPC_Lsd6);
beliefMPC_LA_DENSE_MMTSUB_5_5(beliefMPC_Lsd6, beliefMPC_Yd6);
beliefMPC_LA_DENSE_CHOL_5(beliefMPC_Yd6, beliefMPC_Ld6);
beliefMPC_LA_DENSE_MVMSUB1_5_5(beliefMPC_Lsd6, beliefMPC_yy5, beliefMPC_beta6, beliefMPC_bmy6);
beliefMPC_LA_DENSE_FORWARDSUB_5(beliefMPC_Ld6, beliefMPC_bmy6, beliefMPC_yy6);
beliefMPC_LA_DENSE_MATRIXTFORWARDSUB_5_5(beliefMPC_Ld6, beliefMPC_Ysd7, beliefMPC_Lsd7);
beliefMPC_LA_DENSE_MMTSUB_5_5(beliefMPC_Lsd7, beliefMPC_Yd7);
beliefMPC_LA_DENSE_CHOL_5(beliefMPC_Yd7, beliefMPC_Ld7);
beliefMPC_LA_DENSE_MVMSUB1_5_5(beliefMPC_Lsd7, beliefMPC_yy6, beliefMPC_beta7, beliefMPC_bmy7);
beliefMPC_LA_DENSE_FORWARDSUB_5(beliefMPC_Ld7, beliefMPC_bmy7, beliefMPC_yy7);
beliefMPC_LA_DENSE_MATRIXTFORWARDSUB_5_5(beliefMPC_Ld7, beliefMPC_Ysd8, beliefMPC_Lsd8);
beliefMPC_LA_DENSE_MMTSUB_5_5(beliefMPC_Lsd8, beliefMPC_Yd8);
beliefMPC_LA_DENSE_CHOL_5(beliefMPC_Yd8, beliefMPC_Ld8);
beliefMPC_LA_DENSE_MVMSUB1_5_5(beliefMPC_Lsd8, beliefMPC_yy7, beliefMPC_beta8, beliefMPC_bmy8);
beliefMPC_LA_DENSE_FORWARDSUB_5(beliefMPC_Ld8, beliefMPC_bmy8, beliefMPC_yy8);
beliefMPC_LA_DENSE_BACKWARDSUB_5(beliefMPC_Ld8, beliefMPC_yy8, beliefMPC_dvaff8);
beliefMPC_LA_DENSE_MTVMSUB_5_5(beliefMPC_Lsd8, beliefMPC_dvaff8, beliefMPC_yy7, beliefMPC_bmy7);
beliefMPC_LA_DENSE_BACKWARDSUB_5(beliefMPC_Ld7, beliefMPC_bmy7, beliefMPC_dvaff7);
beliefMPC_LA_DENSE_MTVMSUB_5_5(beliefMPC_Lsd7, beliefMPC_dvaff7, beliefMPC_yy6, beliefMPC_bmy6);
beliefMPC_LA_DENSE_BACKWARDSUB_5(beliefMPC_Ld6, beliefMPC_bmy6, beliefMPC_dvaff6);
beliefMPC_LA_DENSE_MTVMSUB_5_5(beliefMPC_Lsd6, beliefMPC_dvaff6, beliefMPC_yy5, beliefMPC_bmy5);
beliefMPC_LA_DENSE_BACKWARDSUB_5(beliefMPC_Ld5, beliefMPC_bmy5, beliefMPC_dvaff5);
beliefMPC_LA_DENSE_MTVMSUB_5_5(beliefMPC_Lsd5, beliefMPC_dvaff5, beliefMPC_yy4, beliefMPC_bmy4);
beliefMPC_LA_DENSE_BACKWARDSUB_5(beliefMPC_Ld4, beliefMPC_bmy4, beliefMPC_dvaff4);
beliefMPC_LA_DENSE_MTVMSUB_5_5(beliefMPC_Lsd4, beliefMPC_dvaff4, beliefMPC_yy3, beliefMPC_bmy3);
beliefMPC_LA_DENSE_BACKWARDSUB_5(beliefMPC_Ld3, beliefMPC_bmy3, beliefMPC_dvaff3);
beliefMPC_LA_DENSE_MTVMSUB_5_5(beliefMPC_Lsd3, beliefMPC_dvaff3, beliefMPC_yy2, beliefMPC_bmy2);
beliefMPC_LA_DENSE_BACKWARDSUB_5(beliefMPC_Ld2, beliefMPC_bmy2, beliefMPC_dvaff2);
beliefMPC_LA_DENSE_MTVMSUB_5_5(beliefMPC_Lsd2, beliefMPC_dvaff2, beliefMPC_yy1, beliefMPC_bmy1);
beliefMPC_LA_DENSE_BACKWARDSUB_5(beliefMPC_Ld1, beliefMPC_bmy1, beliefMPC_dvaff1);
beliefMPC_LA_DENSE_MTVMSUB_5_10(beliefMPC_Lsd1, beliefMPC_dvaff1, beliefMPC_yy0, beliefMPC_bmy0);
beliefMPC_LA_DENSE_BACKWARDSUB_10(beliefMPC_Ld0, beliefMPC_bmy0, beliefMPC_dvaff0);
beliefMPC_LA_DENSE_MTVM_10_7(params->C1, beliefMPC_dvaff0, beliefMPC_grad_eq0);
beliefMPC_LA_DENSE_MTVM2_5_7_10(params->C2, beliefMPC_dvaff1, beliefMPC_D1, beliefMPC_dvaff0, beliefMPC_grad_eq1);
beliefMPC_LA_DENSE_DIAGZERO_MTVM2_5_7_5(params->C3, beliefMPC_dvaff2, beliefMPC_D2, beliefMPC_dvaff1, beliefMPC_grad_eq2);
beliefMPC_LA_DENSE_DIAGZERO_MTVM2_5_7_5(params->C4, beliefMPC_dvaff3, beliefMPC_D2, beliefMPC_dvaff2, beliefMPC_grad_eq3);
beliefMPC_LA_DENSE_DIAGZERO_MTVM2_5_7_5(params->C5, beliefMPC_dvaff4, beliefMPC_D2, beliefMPC_dvaff3, beliefMPC_grad_eq4);
beliefMPC_LA_DENSE_DIAGZERO_MTVM2_5_7_5(params->C6, beliefMPC_dvaff5, beliefMPC_D2, beliefMPC_dvaff4, beliefMPC_grad_eq5);
beliefMPC_LA_DENSE_DIAGZERO_MTVM2_5_7_5(params->C7, beliefMPC_dvaff6, beliefMPC_D2, beliefMPC_dvaff5, beliefMPC_grad_eq6);
beliefMPC_LA_DENSE_DIAGZERO_MTVM2_5_7_5(params->C8, beliefMPC_dvaff7, beliefMPC_D2, beliefMPC_dvaff6, beliefMPC_grad_eq7);
beliefMPC_LA_DENSE_DIAGZERO_MTVM2_5_7_5(params->C9, beliefMPC_dvaff8, beliefMPC_D2, beliefMPC_dvaff7, beliefMPC_grad_eq8);
beliefMPC_LA_DIAGZERO_MTVM_5_5(beliefMPC_D9, beliefMPC_dvaff8, beliefMPC_grad_eq9);
beliefMPC_LA_VSUB2_68(beliefMPC_rd, beliefMPC_grad_eq, beliefMPC_rd);
beliefMPC_LA_DIAG_FORWARDBACKWARDSUB_7(beliefMPC_Phi0, beliefMPC_rd0, beliefMPC_dzaff0);
beliefMPC_LA_DIAG_FORWARDBACKWARDSUB_7(beliefMPC_Phi1, beliefMPC_rd1, beliefMPC_dzaff1);
beliefMPC_LA_DIAG_FORWARDBACKWARDSUB_7(beliefMPC_Phi2, beliefMPC_rd2, beliefMPC_dzaff2);
beliefMPC_LA_DIAG_FORWARDBACKWARDSUB_7(beliefMPC_Phi3, beliefMPC_rd3, beliefMPC_dzaff3);
beliefMPC_LA_DIAG_FORWARDBACKWARDSUB_7(beliefMPC_Phi4, beliefMPC_rd4, beliefMPC_dzaff4);
beliefMPC_LA_DIAG_FORWARDBACKWARDSUB_7(beliefMPC_Phi5, beliefMPC_rd5, beliefMPC_dzaff5);
beliefMPC_LA_DIAG_FORWARDBACKWARDSUB_7(beliefMPC_Phi6, beliefMPC_rd6, beliefMPC_dzaff6);
beliefMPC_LA_DIAG_FORWARDBACKWARDSUB_7(beliefMPC_Phi7, beliefMPC_rd7, beliefMPC_dzaff7);
beliefMPC_LA_DIAG_FORWARDBACKWARDSUB_7(beliefMPC_Phi8, beliefMPC_rd8, beliefMPC_dzaff8);
beliefMPC_LA_DIAG_FORWARDBACKWARDSUB_5(beliefMPC_Phi9, beliefMPC_rd9, beliefMPC_dzaff9);
beliefMPC_LA_VSUB_INDEXED_7(beliefMPC_dzaff0, beliefMPC_lbIdx0, beliefMPC_rilb0, beliefMPC_dslbaff0);
beliefMPC_LA_VSUB3_7(beliefMPC_llbbyslb0, beliefMPC_dslbaff0, beliefMPC_llb0, beliefMPC_dllbaff0);
beliefMPC_LA_VSUB2_INDEXED_7(beliefMPC_riub0, beliefMPC_dzaff0, beliefMPC_ubIdx0, beliefMPC_dsubaff0);
beliefMPC_LA_VSUB3_7(beliefMPC_lubbysub0, beliefMPC_dsubaff0, beliefMPC_lub0, beliefMPC_dlubaff0);
beliefMPC_LA_VSUB_INDEXED_7(beliefMPC_dzaff1, beliefMPC_lbIdx1, beliefMPC_rilb1, beliefMPC_dslbaff1);
beliefMPC_LA_VSUB3_7(beliefMPC_llbbyslb1, beliefMPC_dslbaff1, beliefMPC_llb1, beliefMPC_dllbaff1);
beliefMPC_LA_VSUB2_INDEXED_7(beliefMPC_riub1, beliefMPC_dzaff1, beliefMPC_ubIdx1, beliefMPC_dsubaff1);
beliefMPC_LA_VSUB3_7(beliefMPC_lubbysub1, beliefMPC_dsubaff1, beliefMPC_lub1, beliefMPC_dlubaff1);
beliefMPC_LA_VSUB_INDEXED_7(beliefMPC_dzaff2, beliefMPC_lbIdx2, beliefMPC_rilb2, beliefMPC_dslbaff2);
beliefMPC_LA_VSUB3_7(beliefMPC_llbbyslb2, beliefMPC_dslbaff2, beliefMPC_llb2, beliefMPC_dllbaff2);
beliefMPC_LA_VSUB2_INDEXED_7(beliefMPC_riub2, beliefMPC_dzaff2, beliefMPC_ubIdx2, beliefMPC_dsubaff2);
beliefMPC_LA_VSUB3_7(beliefMPC_lubbysub2, beliefMPC_dsubaff2, beliefMPC_lub2, beliefMPC_dlubaff2);
beliefMPC_LA_VSUB_INDEXED_7(beliefMPC_dzaff3, beliefMPC_lbIdx3, beliefMPC_rilb3, beliefMPC_dslbaff3);
beliefMPC_LA_VSUB3_7(beliefMPC_llbbyslb3, beliefMPC_dslbaff3, beliefMPC_llb3, beliefMPC_dllbaff3);
beliefMPC_LA_VSUB2_INDEXED_7(beliefMPC_riub3, beliefMPC_dzaff3, beliefMPC_ubIdx3, beliefMPC_dsubaff3);
beliefMPC_LA_VSUB3_7(beliefMPC_lubbysub3, beliefMPC_dsubaff3, beliefMPC_lub3, beliefMPC_dlubaff3);
beliefMPC_LA_VSUB_INDEXED_7(beliefMPC_dzaff4, beliefMPC_lbIdx4, beliefMPC_rilb4, beliefMPC_dslbaff4);
beliefMPC_LA_VSUB3_7(beliefMPC_llbbyslb4, beliefMPC_dslbaff4, beliefMPC_llb4, beliefMPC_dllbaff4);
beliefMPC_LA_VSUB2_INDEXED_7(beliefMPC_riub4, beliefMPC_dzaff4, beliefMPC_ubIdx4, beliefMPC_dsubaff4);
beliefMPC_LA_VSUB3_7(beliefMPC_lubbysub4, beliefMPC_dsubaff4, beliefMPC_lub4, beliefMPC_dlubaff4);
beliefMPC_LA_VSUB_INDEXED_7(beliefMPC_dzaff5, beliefMPC_lbIdx5, beliefMPC_rilb5, beliefMPC_dslbaff5);
beliefMPC_LA_VSUB3_7(beliefMPC_llbbyslb5, beliefMPC_dslbaff5, beliefMPC_llb5, beliefMPC_dllbaff5);
beliefMPC_LA_VSUB2_INDEXED_7(beliefMPC_riub5, beliefMPC_dzaff5, beliefMPC_ubIdx5, beliefMPC_dsubaff5);
beliefMPC_LA_VSUB3_7(beliefMPC_lubbysub5, beliefMPC_dsubaff5, beliefMPC_lub5, beliefMPC_dlubaff5);
beliefMPC_LA_VSUB_INDEXED_7(beliefMPC_dzaff6, beliefMPC_lbIdx6, beliefMPC_rilb6, beliefMPC_dslbaff6);
beliefMPC_LA_VSUB3_7(beliefMPC_llbbyslb6, beliefMPC_dslbaff6, beliefMPC_llb6, beliefMPC_dllbaff6);
beliefMPC_LA_VSUB2_INDEXED_7(beliefMPC_riub6, beliefMPC_dzaff6, beliefMPC_ubIdx6, beliefMPC_dsubaff6);
beliefMPC_LA_VSUB3_7(beliefMPC_lubbysub6, beliefMPC_dsubaff6, beliefMPC_lub6, beliefMPC_dlubaff6);
beliefMPC_LA_VSUB_INDEXED_7(beliefMPC_dzaff7, beliefMPC_lbIdx7, beliefMPC_rilb7, beliefMPC_dslbaff7);
beliefMPC_LA_VSUB3_7(beliefMPC_llbbyslb7, beliefMPC_dslbaff7, beliefMPC_llb7, beliefMPC_dllbaff7);
beliefMPC_LA_VSUB2_INDEXED_7(beliefMPC_riub7, beliefMPC_dzaff7, beliefMPC_ubIdx7, beliefMPC_dsubaff7);
beliefMPC_LA_VSUB3_7(beliefMPC_lubbysub7, beliefMPC_dsubaff7, beliefMPC_lub7, beliefMPC_dlubaff7);
beliefMPC_LA_VSUB_INDEXED_7(beliefMPC_dzaff8, beliefMPC_lbIdx8, beliefMPC_rilb8, beliefMPC_dslbaff8);
beliefMPC_LA_VSUB3_7(beliefMPC_llbbyslb8, beliefMPC_dslbaff8, beliefMPC_llb8, beliefMPC_dllbaff8);
beliefMPC_LA_VSUB2_INDEXED_7(beliefMPC_riub8, beliefMPC_dzaff8, beliefMPC_ubIdx8, beliefMPC_dsubaff8);
beliefMPC_LA_VSUB3_7(beliefMPC_lubbysub8, beliefMPC_dsubaff8, beliefMPC_lub8, beliefMPC_dlubaff8);
beliefMPC_LA_VSUB_INDEXED_5(beliefMPC_dzaff9, beliefMPC_lbIdx9, beliefMPC_rilb9, beliefMPC_dslbaff9);
beliefMPC_LA_VSUB3_5(beliefMPC_llbbyslb9, beliefMPC_dslbaff9, beliefMPC_llb9, beliefMPC_dllbaff9);
beliefMPC_LA_VSUB2_INDEXED_5(beliefMPC_riub9, beliefMPC_dzaff9, beliefMPC_ubIdx9, beliefMPC_dsubaff9);
beliefMPC_LA_VSUB3_5(beliefMPC_lubbysub9, beliefMPC_dsubaff9, beliefMPC_lub9, beliefMPC_dlubaff9);
info->lsit_aff = beliefMPC_LINESEARCH_BACKTRACKING_AFFINE(beliefMPC_l, beliefMPC_s, beliefMPC_dl_aff, beliefMPC_ds_aff, &info->step_aff, &info->mu_aff);
if( info->lsit_aff == beliefMPC_NOPROGRESS ){
exitcode = beliefMPC_NOPROGRESS; break;
}
sigma_3rdroot = info->mu_aff / info->mu;
info->sigma = sigma_3rdroot*sigma_3rdroot*sigma_3rdroot;
musigma = info->mu * info->sigma;
beliefMPC_LA_VSUB5_136(beliefMPC_ds_aff, beliefMPC_dl_aff, musigma, beliefMPC_ccrhs);
beliefMPC_LA_VSUB6_INDEXED_7_7_7(beliefMPC_ccrhsub0, beliefMPC_sub0, beliefMPC_ubIdx0, beliefMPC_ccrhsl0, beliefMPC_slb0, beliefMPC_lbIdx0, beliefMPC_rd0);
beliefMPC_LA_VSUB6_INDEXED_7_7_7(beliefMPC_ccrhsub1, beliefMPC_sub1, beliefMPC_ubIdx1, beliefMPC_ccrhsl1, beliefMPC_slb1, beliefMPC_lbIdx1, beliefMPC_rd1);
beliefMPC_LA_DIAG_FORWARDSUB_7(beliefMPC_Phi0, beliefMPC_rd0, beliefMPC_Lbyrd0);
beliefMPC_LA_DIAG_FORWARDSUB_7(beliefMPC_Phi1, beliefMPC_rd1, beliefMPC_Lbyrd1);
beliefMPC_LA_DENSE_2MVMADD_10_7_7(beliefMPC_V0, beliefMPC_Lbyrd0, beliefMPC_W1, beliefMPC_Lbyrd1, beliefMPC_beta0);
beliefMPC_LA_DENSE_FORWARDSUB_10(beliefMPC_Ld0, beliefMPC_beta0, beliefMPC_yy0);
beliefMPC_LA_VSUB6_INDEXED_7_7_7(beliefMPC_ccrhsub2, beliefMPC_sub2, beliefMPC_ubIdx2, beliefMPC_ccrhsl2, beliefMPC_slb2, beliefMPC_lbIdx2, beliefMPC_rd2);
beliefMPC_LA_DIAG_FORWARDSUB_7(beliefMPC_Phi2, beliefMPC_rd2, beliefMPC_Lbyrd2);
beliefMPC_LA_DENSE_DIAGZERO_2MVMADD_5_7_7(beliefMPC_V1, beliefMPC_Lbyrd1, beliefMPC_W2, beliefMPC_Lbyrd2, beliefMPC_beta1);
beliefMPC_LA_DENSE_MVMSUB1_5_10(beliefMPC_Lsd1, beliefMPC_yy0, beliefMPC_beta1, beliefMPC_bmy1);
beliefMPC_LA_DENSE_FORWARDSUB_5(beliefMPC_Ld1, beliefMPC_bmy1, beliefMPC_yy1);
beliefMPC_LA_VSUB6_INDEXED_7_7_7(beliefMPC_ccrhsub3, beliefMPC_sub3, beliefMPC_ubIdx3, beliefMPC_ccrhsl3, beliefMPC_slb3, beliefMPC_lbIdx3, beliefMPC_rd3);
beliefMPC_LA_DIAG_FORWARDSUB_7(beliefMPC_Phi3, beliefMPC_rd3, beliefMPC_Lbyrd3);
beliefMPC_LA_DENSE_DIAGZERO_2MVMADD_5_7_7(beliefMPC_V2, beliefMPC_Lbyrd2, beliefMPC_W3, beliefMPC_Lbyrd3, beliefMPC_beta2);
beliefMPC_LA_DENSE_MVMSUB1_5_5(beliefMPC_Lsd2, beliefMPC_yy1, beliefMPC_beta2, beliefMPC_bmy2);
beliefMPC_LA_DENSE_FORWARDSUB_5(beliefMPC_Ld2, beliefMPC_bmy2, beliefMPC_yy2);
beliefMPC_LA_VSUB6_INDEXED_7_7_7(beliefMPC_ccrhsub4, beliefMPC_sub4, beliefMPC_ubIdx4, beliefMPC_ccrhsl4, beliefMPC_slb4, beliefMPC_lbIdx4, beliefMPC_rd4);
beliefMPC_LA_DIAG_FORWARDSUB_7(beliefMPC_Phi4, beliefMPC_rd4, beliefMPC_Lbyrd4);
beliefMPC_LA_DENSE_DIAGZERO_2MVMADD_5_7_7(beliefMPC_V3, beliefMPC_Lbyrd3, beliefMPC_W4, beliefMPC_Lbyrd4, beliefMPC_beta3);
beliefMPC_LA_DENSE_MVMSUB1_5_5(beliefMPC_Lsd3, beliefMPC_yy2, beliefMPC_beta3, beliefMPC_bmy3);
beliefMPC_LA_DENSE_FORWARDSUB_5(beliefMPC_Ld3, beliefMPC_bmy3, beliefMPC_yy3);
beliefMPC_LA_VSUB6_INDEXED_7_7_7(beliefMPC_ccrhsub5, beliefMPC_sub5, beliefMPC_ubIdx5, beliefMPC_ccrhsl5, beliefMPC_slb5, beliefMPC_lbIdx5, beliefMPC_rd5);
beliefMPC_LA_DIAG_FORWARDSUB_7(beliefMPC_Phi5, beliefMPC_rd5, beliefMPC_Lbyrd5);
beliefMPC_LA_DENSE_DIAGZERO_2MVMADD_5_7_7(beliefMPC_V4, beliefMPC_Lbyrd4, beliefMPC_W5, beliefMPC_Lbyrd5, beliefMPC_beta4);
beliefMPC_LA_DENSE_MVMSUB1_5_5(beliefMPC_Lsd4, beliefMPC_yy3, beliefMPC_beta4, beliefMPC_bmy4);
beliefMPC_LA_DENSE_FORWARDSUB_5(beliefMPC_Ld4, beliefMPC_bmy4, beliefMPC_yy4);
beliefMPC_LA_VSUB6_INDEXED_7_7_7(beliefMPC_ccrhsub6, beliefMPC_sub6, beliefMPC_ubIdx6, beliefMPC_ccrhsl6, beliefMPC_slb6, beliefMPC_lbIdx6, beliefMPC_rd6);
beliefMPC_LA_DIAG_FORWARDSUB_7(beliefMPC_Phi6, beliefMPC_rd6, beliefMPC_Lbyrd6);
beliefMPC_LA_DENSE_DIAGZERO_2MVMADD_5_7_7(beliefMPC_V5, beliefMPC_Lbyrd5, beliefMPC_W6, beliefMPC_Lbyrd6, beliefMPC_beta5);
beliefMPC_LA_DENSE_MVMSUB1_5_5(beliefMPC_Lsd5, beliefMPC_yy4, beliefMPC_beta5, beliefMPC_bmy5);
beliefMPC_LA_DENSE_FORWARDSUB_5(beliefMPC_Ld5, beliefMPC_bmy5, beliefMPC_yy5);
beliefMPC_LA_VSUB6_INDEXED_7_7_7(beliefMPC_ccrhsub7, beliefMPC_sub7, beliefMPC_ubIdx7, beliefMPC_ccrhsl7, beliefMPC_slb7, beliefMPC_lbIdx7, beliefMPC_rd7);
beliefMPC_LA_DIAG_FORWARDSUB_7(beliefMPC_Phi7, beliefMPC_rd7, beliefMPC_Lbyrd7);
beliefMPC_LA_DENSE_DIAGZERO_2MVMADD_5_7_7(beliefMPC_V6, beliefMPC_Lbyrd6, beliefMPC_W7, beliefMPC_Lbyrd7, beliefMPC_beta6);
beliefMPC_LA_DENSE_MVMSUB1_5_5(beliefMPC_Lsd6, beliefMPC_yy5, beliefMPC_beta6, beliefMPC_bmy6);
beliefMPC_LA_DENSE_FORWARDSUB_5(beliefMPC_Ld6, beliefMPC_bmy6, beliefMPC_yy6);
beliefMPC_LA_VSUB6_INDEXED_7_7_7(beliefMPC_ccrhsub8, beliefMPC_sub8, beliefMPC_ubIdx8, beliefMPC_ccrhsl8, beliefMPC_slb8, beliefMPC_lbIdx8, beliefMPC_rd8);
beliefMPC_LA_DIAG_FORWARDSUB_7(beliefMPC_Phi8, beliefMPC_rd8, beliefMPC_Lbyrd8);
beliefMPC_LA_DENSE_DIAGZERO_2MVMADD_5_7_7(beliefMPC_V7, beliefMPC_Lbyrd7, beliefMPC_W8, beliefMPC_Lbyrd8, beliefMPC_beta7);
beliefMPC_LA_DENSE_MVMSUB1_5_5(beliefMPC_Lsd7, beliefMPC_yy6, beliefMPC_beta7, beliefMPC_bmy7);
beliefMPC_LA_DENSE_FORWARDSUB_5(beliefMPC_Ld7, beliefMPC_bmy7, beliefMPC_yy7);
beliefMPC_LA_VSUB6_INDEXED_5_5_5(beliefMPC_ccrhsub9, beliefMPC_sub9, beliefMPC_ubIdx9, beliefMPC_ccrhsl9, beliefMPC_slb9, beliefMPC_lbIdx9, beliefMPC_rd9);
beliefMPC_LA_DIAG_FORWARDSUB_5(beliefMPC_Phi9, beliefMPC_rd9, beliefMPC_Lbyrd9);
beliefMPC_LA_DENSE_DIAGZERO_2MVMADD_5_7_5(beliefMPC_V8, beliefMPC_Lbyrd8, beliefMPC_W9, beliefMPC_Lbyrd9, beliefMPC_beta8);
beliefMPC_LA_DENSE_MVMSUB1_5_5(beliefMPC_Lsd8, beliefMPC_yy7, beliefMPC_beta8, beliefMPC_bmy8);
beliefMPC_LA_DENSE_FORWARDSUB_5(beliefMPC_Ld8, beliefMPC_bmy8, beliefMPC_yy8);
beliefMPC_LA_DENSE_BACKWARDSUB_5(beliefMPC_Ld8, beliefMPC_yy8, beliefMPC_dvcc8);
beliefMPC_LA_DENSE_MTVMSUB_5_5(beliefMPC_Lsd8, beliefMPC_dvcc8, beliefMPC_yy7, beliefMPC_bmy7);
beliefMPC_LA_DENSE_BACKWARDSUB_5(beliefMPC_Ld7, beliefMPC_bmy7, beliefMPC_dvcc7);
beliefMPC_LA_DENSE_MTVMSUB_5_5(beliefMPC_Lsd7, beliefMPC_dvcc7, beliefMPC_yy6, beliefMPC_bmy6);
beliefMPC_LA_DENSE_BACKWARDSUB_5(beliefMPC_Ld6, beliefMPC_bmy6, beliefMPC_dvcc6);
beliefMPC_LA_DENSE_MTVMSUB_5_5(beliefMPC_Lsd6, beliefMPC_dvcc6, beliefMPC_yy5, beliefMPC_bmy5);
beliefMPC_LA_DENSE_BACKWARDSUB_5(beliefMPC_Ld5, beliefMPC_bmy5, beliefMPC_dvcc5);
beliefMPC_LA_DENSE_MTVMSUB_5_5(beliefMPC_Lsd5, beliefMPC_dvcc5, beliefMPC_yy4, beliefMPC_bmy4);
beliefMPC_LA_DENSE_BACKWARDSUB_5(beliefMPC_Ld4, beliefMPC_bmy4, beliefMPC_dvcc4);
beliefMPC_LA_DENSE_MTVMSUB_5_5(beliefMPC_Lsd4, beliefMPC_dvcc4, beliefMPC_yy3, beliefMPC_bmy3);
beliefMPC_LA_DENSE_BACKWARDSUB_5(beliefMPC_Ld3, beliefMPC_bmy3, beliefMPC_dvcc3);
beliefMPC_LA_DENSE_MTVMSUB_5_5(beliefMPC_Lsd3, beliefMPC_dvcc3, beliefMPC_yy2, beliefMPC_bmy2);
beliefMPC_LA_DENSE_BACKWARDSUB_5(beliefMPC_Ld2, beliefMPC_bmy2, beliefMPC_dvcc2);
beliefMPC_LA_DENSE_MTVMSUB_5_5(beliefMPC_Lsd2, beliefMPC_dvcc2, beliefMPC_yy1, beliefMPC_bmy1);
beliefMPC_LA_DENSE_BACKWARDSUB_5(beliefMPC_Ld1, beliefMPC_bmy1, beliefMPC_dvcc1);
beliefMPC_LA_DENSE_MTVMSUB_5_10(beliefMPC_Lsd1, beliefMPC_dvcc1, beliefMPC_yy0, beliefMPC_bmy0);
beliefMPC_LA_DENSE_BACKWARDSUB_10(beliefMPC_Ld0, beliefMPC_bmy0, beliefMPC_dvcc0);
beliefMPC_LA_DENSE_MTVM_10_7(params->C1, beliefMPC_dvcc0, beliefMPC_grad_eq0);
beliefMPC_LA_DENSE_MTVM2_5_7_10(params->C2, beliefMPC_dvcc1, beliefMPC_D1, beliefMPC_dvcc0, beliefMPC_grad_eq1);
beliefMPC_LA_DENSE_DIAGZERO_MTVM2_5_7_5(params->C3, beliefMPC_dvcc2, beliefMPC_D2, beliefMPC_dvcc1, beliefMPC_grad_eq2);
beliefMPC_LA_DENSE_DIAGZERO_MTVM2_5_7_5(params->C4, beliefMPC_dvcc3, beliefMPC_D2, beliefMPC_dvcc2, beliefMPC_grad_eq3);
beliefMPC_LA_DENSE_DIAGZERO_MTVM2_5_7_5(params->C5, beliefMPC_dvcc4, beliefMPC_D2, beliefMPC_dvcc3, beliefMPC_grad_eq4);
beliefMPC_LA_DENSE_DIAGZERO_MTVM2_5_7_5(params->C6, beliefMPC_dvcc5, beliefMPC_D2, beliefMPC_dvcc4, beliefMPC_grad_eq5);
beliefMPC_LA_DENSE_DIAGZERO_MTVM2_5_7_5(params->C7, beliefMPC_dvcc6, beliefMPC_D2, beliefMPC_dvcc5, beliefMPC_grad_eq6);
beliefMPC_LA_DENSE_DIAGZERO_MTVM2_5_7_5(params->C8, beliefMPC_dvcc7, beliefMPC_D2, beliefMPC_dvcc6, beliefMPC_grad_eq7);
beliefMPC_LA_DENSE_DIAGZERO_MTVM2_5_7_5(params->C9, beliefMPC_dvcc8, beliefMPC_D2, beliefMPC_dvcc7, beliefMPC_grad_eq8);
beliefMPC_LA_DIAGZERO_MTVM_5_5(beliefMPC_D9, beliefMPC_dvcc8, beliefMPC_grad_eq9);
beliefMPC_LA_VSUB_68(beliefMPC_rd, beliefMPC_grad_eq, beliefMPC_rd);
beliefMPC_LA_DIAG_FORWARDBACKWARDSUB_7(beliefMPC_Phi0, beliefMPC_rd0, beliefMPC_dzcc0);
beliefMPC_LA_DIAG_FORWARDBACKWARDSUB_7(beliefMPC_Phi1, beliefMPC_rd1, beliefMPC_dzcc1);
beliefMPC_LA_DIAG_FORWARDBACKWARDSUB_7(beliefMPC_Phi2, beliefMPC_rd2, beliefMPC_dzcc2);
beliefMPC_LA_DIAG_FORWARDBACKWARDSUB_7(beliefMPC_Phi3, beliefMPC_rd3, beliefMPC_dzcc3);
beliefMPC_LA_DIAG_FORWARDBACKWARDSUB_7(beliefMPC_Phi4, beliefMPC_rd4, beliefMPC_dzcc4);
beliefMPC_LA_DIAG_FORWARDBACKWARDSUB_7(beliefMPC_Phi5, beliefMPC_rd5, beliefMPC_dzcc5);
beliefMPC_LA_DIAG_FORWARDBACKWARDSUB_7(beliefMPC_Phi6, beliefMPC_rd6, beliefMPC_dzcc6);
beliefMPC_LA_DIAG_FORWARDBACKWARDSUB_7(beliefMPC_Phi7, beliefMPC_rd7, beliefMPC_dzcc7);
beliefMPC_LA_DIAG_FORWARDBACKWARDSUB_7(beliefMPC_Phi8, beliefMPC_rd8, beliefMPC_dzcc8);
beliefMPC_LA_DIAG_FORWARDBACKWARDSUB_5(beliefMPC_Phi9, beliefMPC_rd9, beliefMPC_dzcc9);
beliefMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_7(beliefMPC_ccrhsl0, beliefMPC_slb0, beliefMPC_llbbyslb0, beliefMPC_dzcc0, beliefMPC_lbIdx0, beliefMPC_dllbcc0);
beliefMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_7(beliefMPC_ccrhsub0, beliefMPC_sub0, beliefMPC_lubbysub0, beliefMPC_dzcc0, beliefMPC_ubIdx0, beliefMPC_dlubcc0);
beliefMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_7(beliefMPC_ccrhsl1, beliefMPC_slb1, beliefMPC_llbbyslb1, beliefMPC_dzcc1, beliefMPC_lbIdx1, beliefMPC_dllbcc1);
beliefMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_7(beliefMPC_ccrhsub1, beliefMPC_sub1, beliefMPC_lubbysub1, beliefMPC_dzcc1, beliefMPC_ubIdx1, beliefMPC_dlubcc1);
beliefMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_7(beliefMPC_ccrhsl2, beliefMPC_slb2, beliefMPC_llbbyslb2, beliefMPC_dzcc2, beliefMPC_lbIdx2, beliefMPC_dllbcc2);
beliefMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_7(beliefMPC_ccrhsub2, beliefMPC_sub2, beliefMPC_lubbysub2, beliefMPC_dzcc2, beliefMPC_ubIdx2, beliefMPC_dlubcc2);
beliefMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_7(beliefMPC_ccrhsl3, beliefMPC_slb3, beliefMPC_llbbyslb3, beliefMPC_dzcc3, beliefMPC_lbIdx3, beliefMPC_dllbcc3);
beliefMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_7(beliefMPC_ccrhsub3, beliefMPC_sub3, beliefMPC_lubbysub3, beliefMPC_dzcc3, beliefMPC_ubIdx3, beliefMPC_dlubcc3);
beliefMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_7(beliefMPC_ccrhsl4, beliefMPC_slb4, beliefMPC_llbbyslb4, beliefMPC_dzcc4, beliefMPC_lbIdx4, beliefMPC_dllbcc4);
beliefMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_7(beliefMPC_ccrhsub4, beliefMPC_sub4, beliefMPC_lubbysub4, beliefMPC_dzcc4, beliefMPC_ubIdx4, beliefMPC_dlubcc4);
beliefMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_7(beliefMPC_ccrhsl5, beliefMPC_slb5, beliefMPC_llbbyslb5, beliefMPC_dzcc5, beliefMPC_lbIdx5, beliefMPC_dllbcc5);
beliefMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_7(beliefMPC_ccrhsub5, beliefMPC_sub5, beliefMPC_lubbysub5, beliefMPC_dzcc5, beliefMPC_ubIdx5, beliefMPC_dlubcc5);
beliefMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_7(beliefMPC_ccrhsl6, beliefMPC_slb6, beliefMPC_llbbyslb6, beliefMPC_dzcc6, beliefMPC_lbIdx6, beliefMPC_dllbcc6);
beliefMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_7(beliefMPC_ccrhsub6, beliefMPC_sub6, beliefMPC_lubbysub6, beliefMPC_dzcc6, beliefMPC_ubIdx6, beliefMPC_dlubcc6);
beliefMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_7(beliefMPC_ccrhsl7, beliefMPC_slb7, beliefMPC_llbbyslb7, beliefMPC_dzcc7, beliefMPC_lbIdx7, beliefMPC_dllbcc7);
beliefMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_7(beliefMPC_ccrhsub7, beliefMPC_sub7, beliefMPC_lubbysub7, beliefMPC_dzcc7, beliefMPC_ubIdx7, beliefMPC_dlubcc7);
beliefMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_7(beliefMPC_ccrhsl8, beliefMPC_slb8, beliefMPC_llbbyslb8, beliefMPC_dzcc8, beliefMPC_lbIdx8, beliefMPC_dllbcc8);
beliefMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_7(beliefMPC_ccrhsub8, beliefMPC_sub8, beliefMPC_lubbysub8, beliefMPC_dzcc8, beliefMPC_ubIdx8, beliefMPC_dlubcc8);
beliefMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_5(beliefMPC_ccrhsl9, beliefMPC_slb9, beliefMPC_llbbyslb9, beliefMPC_dzcc9, beliefMPC_lbIdx9, beliefMPC_dllbcc9);
beliefMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(beliefMPC_ccrhsub9, beliefMPC_sub9, beliefMPC_lubbysub9, beliefMPC_dzcc9, beliefMPC_ubIdx9, beliefMPC_dlubcc9);
beliefMPC_LA_VSUB7_136(beliefMPC_l, beliefMPC_ccrhs, beliefMPC_s, beliefMPC_dl_cc, beliefMPC_ds_cc);
beliefMPC_LA_VADD_68(beliefMPC_dz_cc, beliefMPC_dz_aff);
beliefMPC_LA_VADD_50(beliefMPC_dv_cc, beliefMPC_dv_aff);
beliefMPC_LA_VADD_136(beliefMPC_dl_cc, beliefMPC_dl_aff);
beliefMPC_LA_VADD_136(beliefMPC_ds_cc, beliefMPC_ds_aff);
info->lsit_cc = beliefMPC_LINESEARCH_BACKTRACKING_COMBINED(beliefMPC_z, beliefMPC_v, beliefMPC_l, beliefMPC_s, beliefMPC_dz_cc, beliefMPC_dv_cc, beliefMPC_dl_cc, beliefMPC_ds_cc, &info->step_cc, &info->mu);
if( info->lsit_cc == beliefMPC_NOPROGRESS ){
exitcode = beliefMPC_NOPROGRESS; break;
}
info->it++;
}
output->z1[0] = beliefMPC_z0[0];
output->z1[1] = beliefMPC_z0[1];
output->z1[2] = beliefMPC_z0[2];
output->z1[3] = beliefMPC_z0[3];
output->z1[4] = beliefMPC_z0[4];
output->z1[5] = beliefMPC_z0[5];
output->z1[6] = beliefMPC_z0[6];
output->z2[0] = beliefMPC_z1[0];
output->z2[1] = beliefMPC_z1[1];
output->z2[2] = beliefMPC_z1[2];
output->z2[3] = beliefMPC_z1[3];
output->z2[4] = beliefMPC_z1[4];
output->z2[5] = beliefMPC_z1[5];
output->z2[6] = beliefMPC_z1[6];
output->z3[0] = beliefMPC_z2[0];
output->z3[1] = beliefMPC_z2[1];
output->z3[2] = beliefMPC_z2[2];
output->z3[3] = beliefMPC_z2[3];
output->z3[4] = beliefMPC_z2[4];
output->z3[5] = beliefMPC_z2[5];
output->z3[6] = beliefMPC_z2[6];
output->z4[0] = beliefMPC_z3[0];
output->z4[1] = beliefMPC_z3[1];
output->z4[2] = beliefMPC_z3[2];
output->z4[3] = beliefMPC_z3[3];
output->z4[4] = beliefMPC_z3[4];
output->z4[5] = beliefMPC_z3[5];
output->z4[6] = beliefMPC_z3[6];
output->z5[0] = beliefMPC_z4[0];
output->z5[1] = beliefMPC_z4[1];
output->z5[2] = beliefMPC_z4[2];
output->z5[3] = beliefMPC_z4[3];
output->z5[4] = beliefMPC_z4[4];
output->z5[5] = beliefMPC_z4[5];
output->z5[6] = beliefMPC_z4[6];
output->z6[0] = beliefMPC_z5[0];
output->z6[1] = beliefMPC_z5[1];
output->z6[2] = beliefMPC_z5[2];
output->z6[3] = beliefMPC_z5[3];
output->z6[4] = beliefMPC_z5[4];
output->z6[5] = beliefMPC_z5[5];
output->z6[6] = beliefMPC_z5[6];
output->z7[0] = beliefMPC_z6[0];
output->z7[1] = beliefMPC_z6[1];
output->z7[2] = beliefMPC_z6[2];
output->z7[3] = beliefMPC_z6[3];
output->z7[4] = beliefMPC_z6[4];
output->z7[5] = beliefMPC_z6[5];
output->z7[6] = beliefMPC_z6[6];
output->z8[0] = beliefMPC_z7[0];
output->z8[1] = beliefMPC_z7[1];
output->z8[2] = beliefMPC_z7[2];
output->z8[3] = beliefMPC_z7[3];
output->z8[4] = beliefMPC_z7[4];
output->z8[5] = beliefMPC_z7[5];
output->z8[6] = beliefMPC_z7[6];
output->z9[0] = beliefMPC_z8[0];
output->z9[1] = beliefMPC_z8[1];
output->z9[2] = beliefMPC_z8[2];
output->z9[3] = beliefMPC_z8[3];
output->z9[4] = beliefMPC_z8[4];
output->z9[5] = beliefMPC_z8[5];
output->z9[6] = beliefMPC_z8[6];
output->z10[0] = beliefMPC_z9[0];
output->z10[1] = beliefMPC_z9[1];
output->z10[2] = beliefMPC_z9[2];
output->z10[3] = beliefMPC_z9[3];
output->z10[4] = beliefMPC_z9[4];

#if beliefMPC_SET_TIMING == 1
info->solvetime = beliefMPC_toc(&solvertimer);
#if beliefMPC_SET_PRINTLEVEL > 0 && beliefMPC_SET_TIMING == 1
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
