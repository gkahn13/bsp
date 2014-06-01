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
 * Initializes a vector of length 103 with a value.
 */
void beliefMPC_LA_INITIALIZEVECTOR_103(beliefMPC_FLOAT* vec, beliefMPC_FLOAT value)
{
	int i;
	for( i=0; i<103; i++ )
	{
		vec[i] = value;
	}
}


/*
 * Initializes a vector of length 75 with a value.
 */
void beliefMPC_LA_INITIALIZEVECTOR_75(beliefMPC_FLOAT* vec, beliefMPC_FLOAT value)
{
	int i;
	for( i=0; i<75; i++ )
	{
		vec[i] = value;
	}
}


/*
 * Initializes a vector of length 206 with a value.
 */
void beliefMPC_LA_INITIALIZEVECTOR_206(beliefMPC_FLOAT* vec, beliefMPC_FLOAT value)
{
	int i;
	for( i=0; i<206; i++ )
	{
		vec[i] = value;
	}
}


/* 
 * Calculates a dot product and adds it to a variable: z += x'*y; 
 * This function is for vectors of length 206.
 */
void beliefMPC_LA_DOTACC_206(beliefMPC_FLOAT *x, beliefMPC_FLOAT *y, beliefMPC_FLOAT *z)
{
	int i;
	for( i=0; i<206; i++ ){
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
 * of length 103.
 */
void beliefMPC_LA_VVADD3_103(beliefMPC_FLOAT *u, beliefMPC_FLOAT *v, beliefMPC_FLOAT *w, beliefMPC_FLOAT *z)
{
	int i;
	for( i=0; i<103; i++ ){
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
 * Vector subtraction z = -x - y for vectors of length 103.
 */
void beliefMPC_LA_VSUB2_103(beliefMPC_FLOAT *x, beliefMPC_FLOAT *y, beliefMPC_FLOAT *z)
{
	int i;
	for( i=0; i<103; i++){
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
        for( i=0; i<206; i++ ){
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
        if( i == 206 ){
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
    *mu_aff = mymu / (beliefMPC_FLOAT)206;
    return lsIt;
}


/*
 * Vector subtraction x = u.*v - a where a is a scalar
*  and x,u,v are vectors of length 206.
 */
void beliefMPC_LA_VSUB5_206(beliefMPC_FLOAT *u, beliefMPC_FLOAT *v, beliefMPC_FLOAT a, beliefMPC_FLOAT *x)
{
	int i;
	for( i=0; i<206; i++){
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
 * Vector subtraction z = x - y for vectors of length 103.
 */
void beliefMPC_LA_VSUB_103(beliefMPC_FLOAT *x, beliefMPC_FLOAT *y, beliefMPC_FLOAT *z)
{
	int i;
	for( i=0; i<103; i++){
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
 * Computes ds = -l.\(r + s.*dl) for vectors of length 206.
 */
void beliefMPC_LA_VSUB7_206(beliefMPC_FLOAT *l, beliefMPC_FLOAT *r, beliefMPC_FLOAT *s, beliefMPC_FLOAT *dl, beliefMPC_FLOAT *ds)
{
	int i;
	for( i=0; i<206; i++){
		ds[i] = -(r[i] + s[i]*dl[i])/l[i];
	}
}


/*
 * Vector addition x = x + y for vectors of length 103.
 */
void beliefMPC_LA_VADD_103(beliefMPC_FLOAT *x, beliefMPC_FLOAT *y)
{
	int i;
	for( i=0; i<103; i++){
		x[i] += y[i];
	}
}


/*
 * Vector addition x = x + y for vectors of length 75.
 */
void beliefMPC_LA_VADD_75(beliefMPC_FLOAT *x, beliefMPC_FLOAT *y)
{
	int i;
	for( i=0; i<75; i++){
		x[i] += y[i];
	}
}


/*
 * Vector addition x = x + y for vectors of length 206.
 */
void beliefMPC_LA_VADD_206(beliefMPC_FLOAT *x, beliefMPC_FLOAT *y)
{
	int i;
	for( i=0; i<206; i++){
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
        for( i=0; i<206; i++ ){
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
        if( i == 206 ){
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
    for( i=0; i<103; i++ ){
        z[i] += a_gamma*dz[i];
    }
    
    /* equality constraint multipliers */
    for( i=0; i<75; i++ ){
        v[i] += a_gamma*dv[i];
    }
    
    /* inequality constraint multipliers & slacks, also update mu */
    *mu = 0;
    for( i=0; i<206; i++ ){
        dltemp = l[i] + a_gamma*dl[i]; l[i] = dltemp;
        dstemp = s[i] + a_gamma*ds[i]; s[i] = dstemp;
        *mu += dltemp*dstemp;
    }
    
    *a = a_gamma;
    *mu /= (beliefMPC_FLOAT)206;
    return lsIt;
}




/* VARIABLE DEFINITIONS ------------------------------------------------ */
beliefMPC_FLOAT beliefMPC_z[103];
beliefMPC_FLOAT beliefMPC_v[75];
beliefMPC_FLOAT beliefMPC_dz_aff[103];
beliefMPC_FLOAT beliefMPC_dv_aff[75];
beliefMPC_FLOAT beliefMPC_grad_cost[103];
beliefMPC_FLOAT beliefMPC_grad_eq[103];
beliefMPC_FLOAT beliefMPC_rd[103];
beliefMPC_FLOAT beliefMPC_l[206];
beliefMPC_FLOAT beliefMPC_s[206];
beliefMPC_FLOAT beliefMPC_lbys[206];
beliefMPC_FLOAT beliefMPC_dl_aff[206];
beliefMPC_FLOAT beliefMPC_ds_aff[206];
beliefMPC_FLOAT beliefMPC_dz_cc[103];
beliefMPC_FLOAT beliefMPC_dv_cc[75];
beliefMPC_FLOAT beliefMPC_dl_cc[206];
beliefMPC_FLOAT beliefMPC_ds_cc[206];
beliefMPC_FLOAT beliefMPC_ccrhs[206];
beliefMPC_FLOAT beliefMPC_grad_ineq[103];
beliefMPC_FLOAT beliefMPC_H00[7] = {0.0000000000000000E+000, 0.0000000000000000E+000, 2.0000000000000000E+001, 2.0000000000000000E+001, 2.0000000000000000E+001, 2.0000000000000000E+000, 2.0000000000000000E+000};
beliefMPC_FLOAT beliefMPC_f00[7] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
beliefMPC_FLOAT* beliefMPC_z00 = beliefMPC_z + 0;
beliefMPC_FLOAT* beliefMPC_dzaff00 = beliefMPC_dz_aff + 0;
beliefMPC_FLOAT* beliefMPC_dzcc00 = beliefMPC_dz_cc + 0;
beliefMPC_FLOAT* beliefMPC_rd00 = beliefMPC_rd + 0;
beliefMPC_FLOAT beliefMPC_Lbyrd00[7];
beliefMPC_FLOAT* beliefMPC_grad_cost00 = beliefMPC_grad_cost + 0;
beliefMPC_FLOAT* beliefMPC_grad_eq00 = beliefMPC_grad_eq + 0;
beliefMPC_FLOAT* beliefMPC_grad_ineq00 = beliefMPC_grad_ineq + 0;
beliefMPC_FLOAT beliefMPC_ctv00[7];
beliefMPC_FLOAT* beliefMPC_v00 = beliefMPC_v + 0;
beliefMPC_FLOAT beliefMPC_re00[10];
beliefMPC_FLOAT beliefMPC_beta00[10];
beliefMPC_FLOAT beliefMPC_betacc00[10];
beliefMPC_FLOAT* beliefMPC_dvaff00 = beliefMPC_dv_aff + 0;
beliefMPC_FLOAT* beliefMPC_dvcc00 = beliefMPC_dv_cc + 0;
beliefMPC_FLOAT beliefMPC_V00[70];
beliefMPC_FLOAT beliefMPC_Yd00[55];
beliefMPC_FLOAT beliefMPC_Ld00[55];
beliefMPC_FLOAT beliefMPC_yy00[10];
beliefMPC_FLOAT beliefMPC_bmy00[10];
int beliefMPC_lbIdx00[7] = {0, 1, 2, 3, 4, 5, 6};
beliefMPC_FLOAT* beliefMPC_llb00 = beliefMPC_l + 0;
beliefMPC_FLOAT* beliefMPC_slb00 = beliefMPC_s + 0;
beliefMPC_FLOAT* beliefMPC_llbbyslb00 = beliefMPC_lbys + 0;
beliefMPC_FLOAT beliefMPC_rilb00[7];
beliefMPC_FLOAT* beliefMPC_dllbaff00 = beliefMPC_dl_aff + 0;
beliefMPC_FLOAT* beliefMPC_dslbaff00 = beliefMPC_ds_aff + 0;
beliefMPC_FLOAT* beliefMPC_dllbcc00 = beliefMPC_dl_cc + 0;
beliefMPC_FLOAT* beliefMPC_dslbcc00 = beliefMPC_ds_cc + 0;
beliefMPC_FLOAT* beliefMPC_ccrhsl00 = beliefMPC_ccrhs + 0;
int beliefMPC_ubIdx00[7] = {0, 1, 2, 3, 4, 5, 6};
beliefMPC_FLOAT* beliefMPC_lub00 = beliefMPC_l + 7;
beliefMPC_FLOAT* beliefMPC_sub00 = beliefMPC_s + 7;
beliefMPC_FLOAT* beliefMPC_lubbysub00 = beliefMPC_lbys + 7;
beliefMPC_FLOAT beliefMPC_riub00[7];
beliefMPC_FLOAT* beliefMPC_dlubaff00 = beliefMPC_dl_aff + 7;
beliefMPC_FLOAT* beliefMPC_dsubaff00 = beliefMPC_ds_aff + 7;
beliefMPC_FLOAT* beliefMPC_dlubcc00 = beliefMPC_dl_cc + 7;
beliefMPC_FLOAT* beliefMPC_dsubcc00 = beliefMPC_ds_cc + 7;
beliefMPC_FLOAT* beliefMPC_ccrhsub00 = beliefMPC_ccrhs + 7;
beliefMPC_FLOAT beliefMPC_Phi00[7];
beliefMPC_FLOAT beliefMPC_f01[7] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
beliefMPC_FLOAT* beliefMPC_z01 = beliefMPC_z + 7;
beliefMPC_FLOAT* beliefMPC_dzaff01 = beliefMPC_dz_aff + 7;
beliefMPC_FLOAT* beliefMPC_dzcc01 = beliefMPC_dz_cc + 7;
beliefMPC_FLOAT* beliefMPC_rd01 = beliefMPC_rd + 7;
beliefMPC_FLOAT beliefMPC_Lbyrd01[7];
beliefMPC_FLOAT* beliefMPC_grad_cost01 = beliefMPC_grad_cost + 7;
beliefMPC_FLOAT* beliefMPC_grad_eq01 = beliefMPC_grad_eq + 7;
beliefMPC_FLOAT* beliefMPC_grad_ineq01 = beliefMPC_grad_ineq + 7;
beliefMPC_FLOAT beliefMPC_ctv01[7];
beliefMPC_FLOAT* beliefMPC_v01 = beliefMPC_v + 10;
beliefMPC_FLOAT beliefMPC_re01[5];
beliefMPC_FLOAT beliefMPC_beta01[5];
beliefMPC_FLOAT beliefMPC_betacc01[5];
beliefMPC_FLOAT* beliefMPC_dvaff01 = beliefMPC_dv_aff + 10;
beliefMPC_FLOAT* beliefMPC_dvcc01 = beliefMPC_dv_cc + 10;
beliefMPC_FLOAT beliefMPC_V01[35];
beliefMPC_FLOAT beliefMPC_Yd01[15];
beliefMPC_FLOAT beliefMPC_Ld01[15];
beliefMPC_FLOAT beliefMPC_yy01[5];
beliefMPC_FLOAT beliefMPC_bmy01[5];
int beliefMPC_lbIdx01[7] = {0, 1, 2, 3, 4, 5, 6};
beliefMPC_FLOAT* beliefMPC_llb01 = beliefMPC_l + 14;
beliefMPC_FLOAT* beliefMPC_slb01 = beliefMPC_s + 14;
beliefMPC_FLOAT* beliefMPC_llbbyslb01 = beliefMPC_lbys + 14;
beliefMPC_FLOAT beliefMPC_rilb01[7];
beliefMPC_FLOAT* beliefMPC_dllbaff01 = beliefMPC_dl_aff + 14;
beliefMPC_FLOAT* beliefMPC_dslbaff01 = beliefMPC_ds_aff + 14;
beliefMPC_FLOAT* beliefMPC_dllbcc01 = beliefMPC_dl_cc + 14;
beliefMPC_FLOAT* beliefMPC_dslbcc01 = beliefMPC_ds_cc + 14;
beliefMPC_FLOAT* beliefMPC_ccrhsl01 = beliefMPC_ccrhs + 14;
int beliefMPC_ubIdx01[7] = {0, 1, 2, 3, 4, 5, 6};
beliefMPC_FLOAT* beliefMPC_lub01 = beliefMPC_l + 21;
beliefMPC_FLOAT* beliefMPC_sub01 = beliefMPC_s + 21;
beliefMPC_FLOAT* beliefMPC_lubbysub01 = beliefMPC_lbys + 21;
beliefMPC_FLOAT beliefMPC_riub01[7];
beliefMPC_FLOAT* beliefMPC_dlubaff01 = beliefMPC_dl_aff + 21;
beliefMPC_FLOAT* beliefMPC_dsubaff01 = beliefMPC_ds_aff + 21;
beliefMPC_FLOAT* beliefMPC_dlubcc01 = beliefMPC_dl_cc + 21;
beliefMPC_FLOAT* beliefMPC_dsubcc01 = beliefMPC_ds_cc + 21;
beliefMPC_FLOAT* beliefMPC_ccrhsub01 = beliefMPC_ccrhs + 21;
beliefMPC_FLOAT beliefMPC_Phi01[7];
beliefMPC_FLOAT beliefMPC_D01[70] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, -1.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 
0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, -1.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 
0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, -1.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 
0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, -1.0000000000000000E+000, 0.0000000000000000E+000, 
0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, -1.0000000000000000E+000, 
0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 
0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
beliefMPC_FLOAT beliefMPC_W01[70];
beliefMPC_FLOAT beliefMPC_Ysd01[50];
beliefMPC_FLOAT beliefMPC_Lsd01[50];
beliefMPC_FLOAT beliefMPC_f02[7] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
beliefMPC_FLOAT* beliefMPC_z02 = beliefMPC_z + 14;
beliefMPC_FLOAT* beliefMPC_dzaff02 = beliefMPC_dz_aff + 14;
beliefMPC_FLOAT* beliefMPC_dzcc02 = beliefMPC_dz_cc + 14;
beliefMPC_FLOAT* beliefMPC_rd02 = beliefMPC_rd + 14;
beliefMPC_FLOAT beliefMPC_Lbyrd02[7];
beliefMPC_FLOAT* beliefMPC_grad_cost02 = beliefMPC_grad_cost + 14;
beliefMPC_FLOAT* beliefMPC_grad_eq02 = beliefMPC_grad_eq + 14;
beliefMPC_FLOAT* beliefMPC_grad_ineq02 = beliefMPC_grad_ineq + 14;
beliefMPC_FLOAT beliefMPC_ctv02[7];
beliefMPC_FLOAT* beliefMPC_v02 = beliefMPC_v + 15;
beliefMPC_FLOAT beliefMPC_re02[5];
beliefMPC_FLOAT beliefMPC_beta02[5];
beliefMPC_FLOAT beliefMPC_betacc02[5];
beliefMPC_FLOAT* beliefMPC_dvaff02 = beliefMPC_dv_aff + 15;
beliefMPC_FLOAT* beliefMPC_dvcc02 = beliefMPC_dv_cc + 15;
beliefMPC_FLOAT beliefMPC_V02[35];
beliefMPC_FLOAT beliefMPC_Yd02[15];
beliefMPC_FLOAT beliefMPC_Ld02[15];
beliefMPC_FLOAT beliefMPC_yy02[5];
beliefMPC_FLOAT beliefMPC_bmy02[5];
int beliefMPC_lbIdx02[7] = {0, 1, 2, 3, 4, 5, 6};
beliefMPC_FLOAT* beliefMPC_llb02 = beliefMPC_l + 28;
beliefMPC_FLOAT* beliefMPC_slb02 = beliefMPC_s + 28;
beliefMPC_FLOAT* beliefMPC_llbbyslb02 = beliefMPC_lbys + 28;
beliefMPC_FLOAT beliefMPC_rilb02[7];
beliefMPC_FLOAT* beliefMPC_dllbaff02 = beliefMPC_dl_aff + 28;
beliefMPC_FLOAT* beliefMPC_dslbaff02 = beliefMPC_ds_aff + 28;
beliefMPC_FLOAT* beliefMPC_dllbcc02 = beliefMPC_dl_cc + 28;
beliefMPC_FLOAT* beliefMPC_dslbcc02 = beliefMPC_ds_cc + 28;
beliefMPC_FLOAT* beliefMPC_ccrhsl02 = beliefMPC_ccrhs + 28;
int beliefMPC_ubIdx02[7] = {0, 1, 2, 3, 4, 5, 6};
beliefMPC_FLOAT* beliefMPC_lub02 = beliefMPC_l + 35;
beliefMPC_FLOAT* beliefMPC_sub02 = beliefMPC_s + 35;
beliefMPC_FLOAT* beliefMPC_lubbysub02 = beliefMPC_lbys + 35;
beliefMPC_FLOAT beliefMPC_riub02[7];
beliefMPC_FLOAT* beliefMPC_dlubaff02 = beliefMPC_dl_aff + 35;
beliefMPC_FLOAT* beliefMPC_dsubaff02 = beliefMPC_ds_aff + 35;
beliefMPC_FLOAT* beliefMPC_dlubcc02 = beliefMPC_dl_cc + 35;
beliefMPC_FLOAT* beliefMPC_dsubcc02 = beliefMPC_ds_cc + 35;
beliefMPC_FLOAT* beliefMPC_ccrhsub02 = beliefMPC_ccrhs + 35;
beliefMPC_FLOAT beliefMPC_Phi02[7];
beliefMPC_FLOAT beliefMPC_D02[7] = {-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000};
beliefMPC_FLOAT beliefMPC_W02[7];
beliefMPC_FLOAT beliefMPC_Ysd02[25];
beliefMPC_FLOAT beliefMPC_Lsd02[25];
beliefMPC_FLOAT beliefMPC_f03[7] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
beliefMPC_FLOAT* beliefMPC_z03 = beliefMPC_z + 21;
beliefMPC_FLOAT* beliefMPC_dzaff03 = beliefMPC_dz_aff + 21;
beliefMPC_FLOAT* beliefMPC_dzcc03 = beliefMPC_dz_cc + 21;
beliefMPC_FLOAT* beliefMPC_rd03 = beliefMPC_rd + 21;
beliefMPC_FLOAT beliefMPC_Lbyrd03[7];
beliefMPC_FLOAT* beliefMPC_grad_cost03 = beliefMPC_grad_cost + 21;
beliefMPC_FLOAT* beliefMPC_grad_eq03 = beliefMPC_grad_eq + 21;
beliefMPC_FLOAT* beliefMPC_grad_ineq03 = beliefMPC_grad_ineq + 21;
beliefMPC_FLOAT beliefMPC_ctv03[7];
beliefMPC_FLOAT* beliefMPC_v03 = beliefMPC_v + 20;
beliefMPC_FLOAT beliefMPC_re03[5];
beliefMPC_FLOAT beliefMPC_beta03[5];
beliefMPC_FLOAT beliefMPC_betacc03[5];
beliefMPC_FLOAT* beliefMPC_dvaff03 = beliefMPC_dv_aff + 20;
beliefMPC_FLOAT* beliefMPC_dvcc03 = beliefMPC_dv_cc + 20;
beliefMPC_FLOAT beliefMPC_V03[35];
beliefMPC_FLOAT beliefMPC_Yd03[15];
beliefMPC_FLOAT beliefMPC_Ld03[15];
beliefMPC_FLOAT beliefMPC_yy03[5];
beliefMPC_FLOAT beliefMPC_bmy03[5];
int beliefMPC_lbIdx03[7] = {0, 1, 2, 3, 4, 5, 6};
beliefMPC_FLOAT* beliefMPC_llb03 = beliefMPC_l + 42;
beliefMPC_FLOAT* beliefMPC_slb03 = beliefMPC_s + 42;
beliefMPC_FLOAT* beliefMPC_llbbyslb03 = beliefMPC_lbys + 42;
beliefMPC_FLOAT beliefMPC_rilb03[7];
beliefMPC_FLOAT* beliefMPC_dllbaff03 = beliefMPC_dl_aff + 42;
beliefMPC_FLOAT* beliefMPC_dslbaff03 = beliefMPC_ds_aff + 42;
beliefMPC_FLOAT* beliefMPC_dllbcc03 = beliefMPC_dl_cc + 42;
beliefMPC_FLOAT* beliefMPC_dslbcc03 = beliefMPC_ds_cc + 42;
beliefMPC_FLOAT* beliefMPC_ccrhsl03 = beliefMPC_ccrhs + 42;
int beliefMPC_ubIdx03[7] = {0, 1, 2, 3, 4, 5, 6};
beliefMPC_FLOAT* beliefMPC_lub03 = beliefMPC_l + 49;
beliefMPC_FLOAT* beliefMPC_sub03 = beliefMPC_s + 49;
beliefMPC_FLOAT* beliefMPC_lubbysub03 = beliefMPC_lbys + 49;
beliefMPC_FLOAT beliefMPC_riub03[7];
beliefMPC_FLOAT* beliefMPC_dlubaff03 = beliefMPC_dl_aff + 49;
beliefMPC_FLOAT* beliefMPC_dsubaff03 = beliefMPC_ds_aff + 49;
beliefMPC_FLOAT* beliefMPC_dlubcc03 = beliefMPC_dl_cc + 49;
beliefMPC_FLOAT* beliefMPC_dsubcc03 = beliefMPC_ds_cc + 49;
beliefMPC_FLOAT* beliefMPC_ccrhsub03 = beliefMPC_ccrhs + 49;
beliefMPC_FLOAT beliefMPC_Phi03[7];
beliefMPC_FLOAT beliefMPC_W03[7];
beliefMPC_FLOAT beliefMPC_Ysd03[25];
beliefMPC_FLOAT beliefMPC_Lsd03[25];
beliefMPC_FLOAT beliefMPC_f04[7] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
beliefMPC_FLOAT* beliefMPC_z04 = beliefMPC_z + 28;
beliefMPC_FLOAT* beliefMPC_dzaff04 = beliefMPC_dz_aff + 28;
beliefMPC_FLOAT* beliefMPC_dzcc04 = beliefMPC_dz_cc + 28;
beliefMPC_FLOAT* beliefMPC_rd04 = beliefMPC_rd + 28;
beliefMPC_FLOAT beliefMPC_Lbyrd04[7];
beliefMPC_FLOAT* beliefMPC_grad_cost04 = beliefMPC_grad_cost + 28;
beliefMPC_FLOAT* beliefMPC_grad_eq04 = beliefMPC_grad_eq + 28;
beliefMPC_FLOAT* beliefMPC_grad_ineq04 = beliefMPC_grad_ineq + 28;
beliefMPC_FLOAT beliefMPC_ctv04[7];
beliefMPC_FLOAT* beliefMPC_v04 = beliefMPC_v + 25;
beliefMPC_FLOAT beliefMPC_re04[5];
beliefMPC_FLOAT beliefMPC_beta04[5];
beliefMPC_FLOAT beliefMPC_betacc04[5];
beliefMPC_FLOAT* beliefMPC_dvaff04 = beliefMPC_dv_aff + 25;
beliefMPC_FLOAT* beliefMPC_dvcc04 = beliefMPC_dv_cc + 25;
beliefMPC_FLOAT beliefMPC_V04[35];
beliefMPC_FLOAT beliefMPC_Yd04[15];
beliefMPC_FLOAT beliefMPC_Ld04[15];
beliefMPC_FLOAT beliefMPC_yy04[5];
beliefMPC_FLOAT beliefMPC_bmy04[5];
int beliefMPC_lbIdx04[7] = {0, 1, 2, 3, 4, 5, 6};
beliefMPC_FLOAT* beliefMPC_llb04 = beliefMPC_l + 56;
beliefMPC_FLOAT* beliefMPC_slb04 = beliefMPC_s + 56;
beliefMPC_FLOAT* beliefMPC_llbbyslb04 = beliefMPC_lbys + 56;
beliefMPC_FLOAT beliefMPC_rilb04[7];
beliefMPC_FLOAT* beliefMPC_dllbaff04 = beliefMPC_dl_aff + 56;
beliefMPC_FLOAT* beliefMPC_dslbaff04 = beliefMPC_ds_aff + 56;
beliefMPC_FLOAT* beliefMPC_dllbcc04 = beliefMPC_dl_cc + 56;
beliefMPC_FLOAT* beliefMPC_dslbcc04 = beliefMPC_ds_cc + 56;
beliefMPC_FLOAT* beliefMPC_ccrhsl04 = beliefMPC_ccrhs + 56;
int beliefMPC_ubIdx04[7] = {0, 1, 2, 3, 4, 5, 6};
beliefMPC_FLOAT* beliefMPC_lub04 = beliefMPC_l + 63;
beliefMPC_FLOAT* beliefMPC_sub04 = beliefMPC_s + 63;
beliefMPC_FLOAT* beliefMPC_lubbysub04 = beliefMPC_lbys + 63;
beliefMPC_FLOAT beliefMPC_riub04[7];
beliefMPC_FLOAT* beliefMPC_dlubaff04 = beliefMPC_dl_aff + 63;
beliefMPC_FLOAT* beliefMPC_dsubaff04 = beliefMPC_ds_aff + 63;
beliefMPC_FLOAT* beliefMPC_dlubcc04 = beliefMPC_dl_cc + 63;
beliefMPC_FLOAT* beliefMPC_dsubcc04 = beliefMPC_ds_cc + 63;
beliefMPC_FLOAT* beliefMPC_ccrhsub04 = beliefMPC_ccrhs + 63;
beliefMPC_FLOAT beliefMPC_Phi04[7];
beliefMPC_FLOAT beliefMPC_W04[7];
beliefMPC_FLOAT beliefMPC_Ysd04[25];
beliefMPC_FLOAT beliefMPC_Lsd04[25];
beliefMPC_FLOAT beliefMPC_f05[7] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
beliefMPC_FLOAT* beliefMPC_z05 = beliefMPC_z + 35;
beliefMPC_FLOAT* beliefMPC_dzaff05 = beliefMPC_dz_aff + 35;
beliefMPC_FLOAT* beliefMPC_dzcc05 = beliefMPC_dz_cc + 35;
beliefMPC_FLOAT* beliefMPC_rd05 = beliefMPC_rd + 35;
beliefMPC_FLOAT beliefMPC_Lbyrd05[7];
beliefMPC_FLOAT* beliefMPC_grad_cost05 = beliefMPC_grad_cost + 35;
beliefMPC_FLOAT* beliefMPC_grad_eq05 = beliefMPC_grad_eq + 35;
beliefMPC_FLOAT* beliefMPC_grad_ineq05 = beliefMPC_grad_ineq + 35;
beliefMPC_FLOAT beliefMPC_ctv05[7];
beliefMPC_FLOAT* beliefMPC_v05 = beliefMPC_v + 30;
beliefMPC_FLOAT beliefMPC_re05[5];
beliefMPC_FLOAT beliefMPC_beta05[5];
beliefMPC_FLOAT beliefMPC_betacc05[5];
beliefMPC_FLOAT* beliefMPC_dvaff05 = beliefMPC_dv_aff + 30;
beliefMPC_FLOAT* beliefMPC_dvcc05 = beliefMPC_dv_cc + 30;
beliefMPC_FLOAT beliefMPC_V05[35];
beliefMPC_FLOAT beliefMPC_Yd05[15];
beliefMPC_FLOAT beliefMPC_Ld05[15];
beliefMPC_FLOAT beliefMPC_yy05[5];
beliefMPC_FLOAT beliefMPC_bmy05[5];
int beliefMPC_lbIdx05[7] = {0, 1, 2, 3, 4, 5, 6};
beliefMPC_FLOAT* beliefMPC_llb05 = beliefMPC_l + 70;
beliefMPC_FLOAT* beliefMPC_slb05 = beliefMPC_s + 70;
beliefMPC_FLOAT* beliefMPC_llbbyslb05 = beliefMPC_lbys + 70;
beliefMPC_FLOAT beliefMPC_rilb05[7];
beliefMPC_FLOAT* beliefMPC_dllbaff05 = beliefMPC_dl_aff + 70;
beliefMPC_FLOAT* beliefMPC_dslbaff05 = beliefMPC_ds_aff + 70;
beliefMPC_FLOAT* beliefMPC_dllbcc05 = beliefMPC_dl_cc + 70;
beliefMPC_FLOAT* beliefMPC_dslbcc05 = beliefMPC_ds_cc + 70;
beliefMPC_FLOAT* beliefMPC_ccrhsl05 = beliefMPC_ccrhs + 70;
int beliefMPC_ubIdx05[7] = {0, 1, 2, 3, 4, 5, 6};
beliefMPC_FLOAT* beliefMPC_lub05 = beliefMPC_l + 77;
beliefMPC_FLOAT* beliefMPC_sub05 = beliefMPC_s + 77;
beliefMPC_FLOAT* beliefMPC_lubbysub05 = beliefMPC_lbys + 77;
beliefMPC_FLOAT beliefMPC_riub05[7];
beliefMPC_FLOAT* beliefMPC_dlubaff05 = beliefMPC_dl_aff + 77;
beliefMPC_FLOAT* beliefMPC_dsubaff05 = beliefMPC_ds_aff + 77;
beliefMPC_FLOAT* beliefMPC_dlubcc05 = beliefMPC_dl_cc + 77;
beliefMPC_FLOAT* beliefMPC_dsubcc05 = beliefMPC_ds_cc + 77;
beliefMPC_FLOAT* beliefMPC_ccrhsub05 = beliefMPC_ccrhs + 77;
beliefMPC_FLOAT beliefMPC_Phi05[7];
beliefMPC_FLOAT beliefMPC_W05[7];
beliefMPC_FLOAT beliefMPC_Ysd05[25];
beliefMPC_FLOAT beliefMPC_Lsd05[25];
beliefMPC_FLOAT beliefMPC_f06[7] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
beliefMPC_FLOAT* beliefMPC_z06 = beliefMPC_z + 42;
beliefMPC_FLOAT* beliefMPC_dzaff06 = beliefMPC_dz_aff + 42;
beliefMPC_FLOAT* beliefMPC_dzcc06 = beliefMPC_dz_cc + 42;
beliefMPC_FLOAT* beliefMPC_rd06 = beliefMPC_rd + 42;
beliefMPC_FLOAT beliefMPC_Lbyrd06[7];
beliefMPC_FLOAT* beliefMPC_grad_cost06 = beliefMPC_grad_cost + 42;
beliefMPC_FLOAT* beliefMPC_grad_eq06 = beliefMPC_grad_eq + 42;
beliefMPC_FLOAT* beliefMPC_grad_ineq06 = beliefMPC_grad_ineq + 42;
beliefMPC_FLOAT beliefMPC_ctv06[7];
beliefMPC_FLOAT* beliefMPC_v06 = beliefMPC_v + 35;
beliefMPC_FLOAT beliefMPC_re06[5];
beliefMPC_FLOAT beliefMPC_beta06[5];
beliefMPC_FLOAT beliefMPC_betacc06[5];
beliefMPC_FLOAT* beliefMPC_dvaff06 = beliefMPC_dv_aff + 35;
beliefMPC_FLOAT* beliefMPC_dvcc06 = beliefMPC_dv_cc + 35;
beliefMPC_FLOAT beliefMPC_V06[35];
beliefMPC_FLOAT beliefMPC_Yd06[15];
beliefMPC_FLOAT beliefMPC_Ld06[15];
beliefMPC_FLOAT beliefMPC_yy06[5];
beliefMPC_FLOAT beliefMPC_bmy06[5];
int beliefMPC_lbIdx06[7] = {0, 1, 2, 3, 4, 5, 6};
beliefMPC_FLOAT* beliefMPC_llb06 = beliefMPC_l + 84;
beliefMPC_FLOAT* beliefMPC_slb06 = beliefMPC_s + 84;
beliefMPC_FLOAT* beliefMPC_llbbyslb06 = beliefMPC_lbys + 84;
beliefMPC_FLOAT beliefMPC_rilb06[7];
beliefMPC_FLOAT* beliefMPC_dllbaff06 = beliefMPC_dl_aff + 84;
beliefMPC_FLOAT* beliefMPC_dslbaff06 = beliefMPC_ds_aff + 84;
beliefMPC_FLOAT* beliefMPC_dllbcc06 = beliefMPC_dl_cc + 84;
beliefMPC_FLOAT* beliefMPC_dslbcc06 = beliefMPC_ds_cc + 84;
beliefMPC_FLOAT* beliefMPC_ccrhsl06 = beliefMPC_ccrhs + 84;
int beliefMPC_ubIdx06[7] = {0, 1, 2, 3, 4, 5, 6};
beliefMPC_FLOAT* beliefMPC_lub06 = beliefMPC_l + 91;
beliefMPC_FLOAT* beliefMPC_sub06 = beliefMPC_s + 91;
beliefMPC_FLOAT* beliefMPC_lubbysub06 = beliefMPC_lbys + 91;
beliefMPC_FLOAT beliefMPC_riub06[7];
beliefMPC_FLOAT* beliefMPC_dlubaff06 = beliefMPC_dl_aff + 91;
beliefMPC_FLOAT* beliefMPC_dsubaff06 = beliefMPC_ds_aff + 91;
beliefMPC_FLOAT* beliefMPC_dlubcc06 = beliefMPC_dl_cc + 91;
beliefMPC_FLOAT* beliefMPC_dsubcc06 = beliefMPC_ds_cc + 91;
beliefMPC_FLOAT* beliefMPC_ccrhsub06 = beliefMPC_ccrhs + 91;
beliefMPC_FLOAT beliefMPC_Phi06[7];
beliefMPC_FLOAT beliefMPC_W06[7];
beliefMPC_FLOAT beliefMPC_Ysd06[25];
beliefMPC_FLOAT beliefMPC_Lsd06[25];
beliefMPC_FLOAT beliefMPC_f07[7] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
beliefMPC_FLOAT* beliefMPC_z07 = beliefMPC_z + 49;
beliefMPC_FLOAT* beliefMPC_dzaff07 = beliefMPC_dz_aff + 49;
beliefMPC_FLOAT* beliefMPC_dzcc07 = beliefMPC_dz_cc + 49;
beliefMPC_FLOAT* beliefMPC_rd07 = beliefMPC_rd + 49;
beliefMPC_FLOAT beliefMPC_Lbyrd07[7];
beliefMPC_FLOAT* beliefMPC_grad_cost07 = beliefMPC_grad_cost + 49;
beliefMPC_FLOAT* beliefMPC_grad_eq07 = beliefMPC_grad_eq + 49;
beliefMPC_FLOAT* beliefMPC_grad_ineq07 = beliefMPC_grad_ineq + 49;
beliefMPC_FLOAT beliefMPC_ctv07[7];
beliefMPC_FLOAT* beliefMPC_v07 = beliefMPC_v + 40;
beliefMPC_FLOAT beliefMPC_re07[5];
beliefMPC_FLOAT beliefMPC_beta07[5];
beliefMPC_FLOAT beliefMPC_betacc07[5];
beliefMPC_FLOAT* beliefMPC_dvaff07 = beliefMPC_dv_aff + 40;
beliefMPC_FLOAT* beliefMPC_dvcc07 = beliefMPC_dv_cc + 40;
beliefMPC_FLOAT beliefMPC_V07[35];
beliefMPC_FLOAT beliefMPC_Yd07[15];
beliefMPC_FLOAT beliefMPC_Ld07[15];
beliefMPC_FLOAT beliefMPC_yy07[5];
beliefMPC_FLOAT beliefMPC_bmy07[5];
int beliefMPC_lbIdx07[7] = {0, 1, 2, 3, 4, 5, 6};
beliefMPC_FLOAT* beliefMPC_llb07 = beliefMPC_l + 98;
beliefMPC_FLOAT* beliefMPC_slb07 = beliefMPC_s + 98;
beliefMPC_FLOAT* beliefMPC_llbbyslb07 = beliefMPC_lbys + 98;
beliefMPC_FLOAT beliefMPC_rilb07[7];
beliefMPC_FLOAT* beliefMPC_dllbaff07 = beliefMPC_dl_aff + 98;
beliefMPC_FLOAT* beliefMPC_dslbaff07 = beliefMPC_ds_aff + 98;
beliefMPC_FLOAT* beliefMPC_dllbcc07 = beliefMPC_dl_cc + 98;
beliefMPC_FLOAT* beliefMPC_dslbcc07 = beliefMPC_ds_cc + 98;
beliefMPC_FLOAT* beliefMPC_ccrhsl07 = beliefMPC_ccrhs + 98;
int beliefMPC_ubIdx07[7] = {0, 1, 2, 3, 4, 5, 6};
beliefMPC_FLOAT* beliefMPC_lub07 = beliefMPC_l + 105;
beliefMPC_FLOAT* beliefMPC_sub07 = beliefMPC_s + 105;
beliefMPC_FLOAT* beliefMPC_lubbysub07 = beliefMPC_lbys + 105;
beliefMPC_FLOAT beliefMPC_riub07[7];
beliefMPC_FLOAT* beliefMPC_dlubaff07 = beliefMPC_dl_aff + 105;
beliefMPC_FLOAT* beliefMPC_dsubaff07 = beliefMPC_ds_aff + 105;
beliefMPC_FLOAT* beliefMPC_dlubcc07 = beliefMPC_dl_cc + 105;
beliefMPC_FLOAT* beliefMPC_dsubcc07 = beliefMPC_ds_cc + 105;
beliefMPC_FLOAT* beliefMPC_ccrhsub07 = beliefMPC_ccrhs + 105;
beliefMPC_FLOAT beliefMPC_Phi07[7];
beliefMPC_FLOAT beliefMPC_W07[7];
beliefMPC_FLOAT beliefMPC_Ysd07[25];
beliefMPC_FLOAT beliefMPC_Lsd07[25];
beliefMPC_FLOAT beliefMPC_f08[7] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
beliefMPC_FLOAT* beliefMPC_z08 = beliefMPC_z + 56;
beliefMPC_FLOAT* beliefMPC_dzaff08 = beliefMPC_dz_aff + 56;
beliefMPC_FLOAT* beliefMPC_dzcc08 = beliefMPC_dz_cc + 56;
beliefMPC_FLOAT* beliefMPC_rd08 = beliefMPC_rd + 56;
beliefMPC_FLOAT beliefMPC_Lbyrd08[7];
beliefMPC_FLOAT* beliefMPC_grad_cost08 = beliefMPC_grad_cost + 56;
beliefMPC_FLOAT* beliefMPC_grad_eq08 = beliefMPC_grad_eq + 56;
beliefMPC_FLOAT* beliefMPC_grad_ineq08 = beliefMPC_grad_ineq + 56;
beliefMPC_FLOAT beliefMPC_ctv08[7];
beliefMPC_FLOAT* beliefMPC_v08 = beliefMPC_v + 45;
beliefMPC_FLOAT beliefMPC_re08[5];
beliefMPC_FLOAT beliefMPC_beta08[5];
beliefMPC_FLOAT beliefMPC_betacc08[5];
beliefMPC_FLOAT* beliefMPC_dvaff08 = beliefMPC_dv_aff + 45;
beliefMPC_FLOAT* beliefMPC_dvcc08 = beliefMPC_dv_cc + 45;
beliefMPC_FLOAT beliefMPC_V08[35];
beliefMPC_FLOAT beliefMPC_Yd08[15];
beliefMPC_FLOAT beliefMPC_Ld08[15];
beliefMPC_FLOAT beliefMPC_yy08[5];
beliefMPC_FLOAT beliefMPC_bmy08[5];
int beliefMPC_lbIdx08[7] = {0, 1, 2, 3, 4, 5, 6};
beliefMPC_FLOAT* beliefMPC_llb08 = beliefMPC_l + 112;
beliefMPC_FLOAT* beliefMPC_slb08 = beliefMPC_s + 112;
beliefMPC_FLOAT* beliefMPC_llbbyslb08 = beliefMPC_lbys + 112;
beliefMPC_FLOAT beliefMPC_rilb08[7];
beliefMPC_FLOAT* beliefMPC_dllbaff08 = beliefMPC_dl_aff + 112;
beliefMPC_FLOAT* beliefMPC_dslbaff08 = beliefMPC_ds_aff + 112;
beliefMPC_FLOAT* beliefMPC_dllbcc08 = beliefMPC_dl_cc + 112;
beliefMPC_FLOAT* beliefMPC_dslbcc08 = beliefMPC_ds_cc + 112;
beliefMPC_FLOAT* beliefMPC_ccrhsl08 = beliefMPC_ccrhs + 112;
int beliefMPC_ubIdx08[7] = {0, 1, 2, 3, 4, 5, 6};
beliefMPC_FLOAT* beliefMPC_lub08 = beliefMPC_l + 119;
beliefMPC_FLOAT* beliefMPC_sub08 = beliefMPC_s + 119;
beliefMPC_FLOAT* beliefMPC_lubbysub08 = beliefMPC_lbys + 119;
beliefMPC_FLOAT beliefMPC_riub08[7];
beliefMPC_FLOAT* beliefMPC_dlubaff08 = beliefMPC_dl_aff + 119;
beliefMPC_FLOAT* beliefMPC_dsubaff08 = beliefMPC_ds_aff + 119;
beliefMPC_FLOAT* beliefMPC_dlubcc08 = beliefMPC_dl_cc + 119;
beliefMPC_FLOAT* beliefMPC_dsubcc08 = beliefMPC_ds_cc + 119;
beliefMPC_FLOAT* beliefMPC_ccrhsub08 = beliefMPC_ccrhs + 119;
beliefMPC_FLOAT beliefMPC_Phi08[7];
beliefMPC_FLOAT beliefMPC_W08[7];
beliefMPC_FLOAT beliefMPC_Ysd08[25];
beliefMPC_FLOAT beliefMPC_Lsd08[25];
beliefMPC_FLOAT beliefMPC_f09[7] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
beliefMPC_FLOAT* beliefMPC_z09 = beliefMPC_z + 63;
beliefMPC_FLOAT* beliefMPC_dzaff09 = beliefMPC_dz_aff + 63;
beliefMPC_FLOAT* beliefMPC_dzcc09 = beliefMPC_dz_cc + 63;
beliefMPC_FLOAT* beliefMPC_rd09 = beliefMPC_rd + 63;
beliefMPC_FLOAT beliefMPC_Lbyrd09[7];
beliefMPC_FLOAT* beliefMPC_grad_cost09 = beliefMPC_grad_cost + 63;
beliefMPC_FLOAT* beliefMPC_grad_eq09 = beliefMPC_grad_eq + 63;
beliefMPC_FLOAT* beliefMPC_grad_ineq09 = beliefMPC_grad_ineq + 63;
beliefMPC_FLOAT beliefMPC_ctv09[7];
beliefMPC_FLOAT* beliefMPC_v09 = beliefMPC_v + 50;
beliefMPC_FLOAT beliefMPC_re09[5];
beliefMPC_FLOAT beliefMPC_beta09[5];
beliefMPC_FLOAT beliefMPC_betacc09[5];
beliefMPC_FLOAT* beliefMPC_dvaff09 = beliefMPC_dv_aff + 50;
beliefMPC_FLOAT* beliefMPC_dvcc09 = beliefMPC_dv_cc + 50;
beliefMPC_FLOAT beliefMPC_V09[35];
beliefMPC_FLOAT beliefMPC_Yd09[15];
beliefMPC_FLOAT beliefMPC_Ld09[15];
beliefMPC_FLOAT beliefMPC_yy09[5];
beliefMPC_FLOAT beliefMPC_bmy09[5];
int beliefMPC_lbIdx09[7] = {0, 1, 2, 3, 4, 5, 6};
beliefMPC_FLOAT* beliefMPC_llb09 = beliefMPC_l + 126;
beliefMPC_FLOAT* beliefMPC_slb09 = beliefMPC_s + 126;
beliefMPC_FLOAT* beliefMPC_llbbyslb09 = beliefMPC_lbys + 126;
beliefMPC_FLOAT beliefMPC_rilb09[7];
beliefMPC_FLOAT* beliefMPC_dllbaff09 = beliefMPC_dl_aff + 126;
beliefMPC_FLOAT* beliefMPC_dslbaff09 = beliefMPC_ds_aff + 126;
beliefMPC_FLOAT* beliefMPC_dllbcc09 = beliefMPC_dl_cc + 126;
beliefMPC_FLOAT* beliefMPC_dslbcc09 = beliefMPC_ds_cc + 126;
beliefMPC_FLOAT* beliefMPC_ccrhsl09 = beliefMPC_ccrhs + 126;
int beliefMPC_ubIdx09[7] = {0, 1, 2, 3, 4, 5, 6};
beliefMPC_FLOAT* beliefMPC_lub09 = beliefMPC_l + 133;
beliefMPC_FLOAT* beliefMPC_sub09 = beliefMPC_s + 133;
beliefMPC_FLOAT* beliefMPC_lubbysub09 = beliefMPC_lbys + 133;
beliefMPC_FLOAT beliefMPC_riub09[7];
beliefMPC_FLOAT* beliefMPC_dlubaff09 = beliefMPC_dl_aff + 133;
beliefMPC_FLOAT* beliefMPC_dsubaff09 = beliefMPC_ds_aff + 133;
beliefMPC_FLOAT* beliefMPC_dlubcc09 = beliefMPC_dl_cc + 133;
beliefMPC_FLOAT* beliefMPC_dsubcc09 = beliefMPC_ds_cc + 133;
beliefMPC_FLOAT* beliefMPC_ccrhsub09 = beliefMPC_ccrhs + 133;
beliefMPC_FLOAT beliefMPC_Phi09[7];
beliefMPC_FLOAT beliefMPC_W09[7];
beliefMPC_FLOAT beliefMPC_Ysd09[25];
beliefMPC_FLOAT beliefMPC_Lsd09[25];
beliefMPC_FLOAT beliefMPC_f10[7] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
beliefMPC_FLOAT* beliefMPC_z10 = beliefMPC_z + 70;
beliefMPC_FLOAT* beliefMPC_dzaff10 = beliefMPC_dz_aff + 70;
beliefMPC_FLOAT* beliefMPC_dzcc10 = beliefMPC_dz_cc + 70;
beliefMPC_FLOAT* beliefMPC_rd10 = beliefMPC_rd + 70;
beliefMPC_FLOAT beliefMPC_Lbyrd10[7];
beliefMPC_FLOAT* beliefMPC_grad_cost10 = beliefMPC_grad_cost + 70;
beliefMPC_FLOAT* beliefMPC_grad_eq10 = beliefMPC_grad_eq + 70;
beliefMPC_FLOAT* beliefMPC_grad_ineq10 = beliefMPC_grad_ineq + 70;
beliefMPC_FLOAT beliefMPC_ctv10[7];
beliefMPC_FLOAT* beliefMPC_v10 = beliefMPC_v + 55;
beliefMPC_FLOAT beliefMPC_re10[5];
beliefMPC_FLOAT beliefMPC_beta10[5];
beliefMPC_FLOAT beliefMPC_betacc10[5];
beliefMPC_FLOAT* beliefMPC_dvaff10 = beliefMPC_dv_aff + 55;
beliefMPC_FLOAT* beliefMPC_dvcc10 = beliefMPC_dv_cc + 55;
beliefMPC_FLOAT beliefMPC_V10[35];
beliefMPC_FLOAT beliefMPC_Yd10[15];
beliefMPC_FLOAT beliefMPC_Ld10[15];
beliefMPC_FLOAT beliefMPC_yy10[5];
beliefMPC_FLOAT beliefMPC_bmy10[5];
int beliefMPC_lbIdx10[7] = {0, 1, 2, 3, 4, 5, 6};
beliefMPC_FLOAT* beliefMPC_llb10 = beliefMPC_l + 140;
beliefMPC_FLOAT* beliefMPC_slb10 = beliefMPC_s + 140;
beliefMPC_FLOAT* beliefMPC_llbbyslb10 = beliefMPC_lbys + 140;
beliefMPC_FLOAT beliefMPC_rilb10[7];
beliefMPC_FLOAT* beliefMPC_dllbaff10 = beliefMPC_dl_aff + 140;
beliefMPC_FLOAT* beliefMPC_dslbaff10 = beliefMPC_ds_aff + 140;
beliefMPC_FLOAT* beliefMPC_dllbcc10 = beliefMPC_dl_cc + 140;
beliefMPC_FLOAT* beliefMPC_dslbcc10 = beliefMPC_ds_cc + 140;
beliefMPC_FLOAT* beliefMPC_ccrhsl10 = beliefMPC_ccrhs + 140;
int beliefMPC_ubIdx10[7] = {0, 1, 2, 3, 4, 5, 6};
beliefMPC_FLOAT* beliefMPC_lub10 = beliefMPC_l + 147;
beliefMPC_FLOAT* beliefMPC_sub10 = beliefMPC_s + 147;
beliefMPC_FLOAT* beliefMPC_lubbysub10 = beliefMPC_lbys + 147;
beliefMPC_FLOAT beliefMPC_riub10[7];
beliefMPC_FLOAT* beliefMPC_dlubaff10 = beliefMPC_dl_aff + 147;
beliefMPC_FLOAT* beliefMPC_dsubaff10 = beliefMPC_ds_aff + 147;
beliefMPC_FLOAT* beliefMPC_dlubcc10 = beliefMPC_dl_cc + 147;
beliefMPC_FLOAT* beliefMPC_dsubcc10 = beliefMPC_ds_cc + 147;
beliefMPC_FLOAT* beliefMPC_ccrhsub10 = beliefMPC_ccrhs + 147;
beliefMPC_FLOAT beliefMPC_Phi10[7];
beliefMPC_FLOAT beliefMPC_W10[7];
beliefMPC_FLOAT beliefMPC_Ysd10[25];
beliefMPC_FLOAT beliefMPC_Lsd10[25];
beliefMPC_FLOAT beliefMPC_f11[7] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
beliefMPC_FLOAT* beliefMPC_z11 = beliefMPC_z + 77;
beliefMPC_FLOAT* beliefMPC_dzaff11 = beliefMPC_dz_aff + 77;
beliefMPC_FLOAT* beliefMPC_dzcc11 = beliefMPC_dz_cc + 77;
beliefMPC_FLOAT* beliefMPC_rd11 = beliefMPC_rd + 77;
beliefMPC_FLOAT beliefMPC_Lbyrd11[7];
beliefMPC_FLOAT* beliefMPC_grad_cost11 = beliefMPC_grad_cost + 77;
beliefMPC_FLOAT* beliefMPC_grad_eq11 = beliefMPC_grad_eq + 77;
beliefMPC_FLOAT* beliefMPC_grad_ineq11 = beliefMPC_grad_ineq + 77;
beliefMPC_FLOAT beliefMPC_ctv11[7];
beliefMPC_FLOAT* beliefMPC_v11 = beliefMPC_v + 60;
beliefMPC_FLOAT beliefMPC_re11[5];
beliefMPC_FLOAT beliefMPC_beta11[5];
beliefMPC_FLOAT beliefMPC_betacc11[5];
beliefMPC_FLOAT* beliefMPC_dvaff11 = beliefMPC_dv_aff + 60;
beliefMPC_FLOAT* beliefMPC_dvcc11 = beliefMPC_dv_cc + 60;
beliefMPC_FLOAT beliefMPC_V11[35];
beliefMPC_FLOAT beliefMPC_Yd11[15];
beliefMPC_FLOAT beliefMPC_Ld11[15];
beliefMPC_FLOAT beliefMPC_yy11[5];
beliefMPC_FLOAT beliefMPC_bmy11[5];
int beliefMPC_lbIdx11[7] = {0, 1, 2, 3, 4, 5, 6};
beliefMPC_FLOAT* beliefMPC_llb11 = beliefMPC_l + 154;
beliefMPC_FLOAT* beliefMPC_slb11 = beliefMPC_s + 154;
beliefMPC_FLOAT* beliefMPC_llbbyslb11 = beliefMPC_lbys + 154;
beliefMPC_FLOAT beliefMPC_rilb11[7];
beliefMPC_FLOAT* beliefMPC_dllbaff11 = beliefMPC_dl_aff + 154;
beliefMPC_FLOAT* beliefMPC_dslbaff11 = beliefMPC_ds_aff + 154;
beliefMPC_FLOAT* beliefMPC_dllbcc11 = beliefMPC_dl_cc + 154;
beliefMPC_FLOAT* beliefMPC_dslbcc11 = beliefMPC_ds_cc + 154;
beliefMPC_FLOAT* beliefMPC_ccrhsl11 = beliefMPC_ccrhs + 154;
int beliefMPC_ubIdx11[7] = {0, 1, 2, 3, 4, 5, 6};
beliefMPC_FLOAT* beliefMPC_lub11 = beliefMPC_l + 161;
beliefMPC_FLOAT* beliefMPC_sub11 = beliefMPC_s + 161;
beliefMPC_FLOAT* beliefMPC_lubbysub11 = beliefMPC_lbys + 161;
beliefMPC_FLOAT beliefMPC_riub11[7];
beliefMPC_FLOAT* beliefMPC_dlubaff11 = beliefMPC_dl_aff + 161;
beliefMPC_FLOAT* beliefMPC_dsubaff11 = beliefMPC_ds_aff + 161;
beliefMPC_FLOAT* beliefMPC_dlubcc11 = beliefMPC_dl_cc + 161;
beliefMPC_FLOAT* beliefMPC_dsubcc11 = beliefMPC_ds_cc + 161;
beliefMPC_FLOAT* beliefMPC_ccrhsub11 = beliefMPC_ccrhs + 161;
beliefMPC_FLOAT beliefMPC_Phi11[7];
beliefMPC_FLOAT beliefMPC_W11[7];
beliefMPC_FLOAT beliefMPC_Ysd11[25];
beliefMPC_FLOAT beliefMPC_Lsd11[25];
beliefMPC_FLOAT beliefMPC_f12[7] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
beliefMPC_FLOAT* beliefMPC_z12 = beliefMPC_z + 84;
beliefMPC_FLOAT* beliefMPC_dzaff12 = beliefMPC_dz_aff + 84;
beliefMPC_FLOAT* beliefMPC_dzcc12 = beliefMPC_dz_cc + 84;
beliefMPC_FLOAT* beliefMPC_rd12 = beliefMPC_rd + 84;
beliefMPC_FLOAT beliefMPC_Lbyrd12[7];
beliefMPC_FLOAT* beliefMPC_grad_cost12 = beliefMPC_grad_cost + 84;
beliefMPC_FLOAT* beliefMPC_grad_eq12 = beliefMPC_grad_eq + 84;
beliefMPC_FLOAT* beliefMPC_grad_ineq12 = beliefMPC_grad_ineq + 84;
beliefMPC_FLOAT beliefMPC_ctv12[7];
beliefMPC_FLOAT* beliefMPC_v12 = beliefMPC_v + 65;
beliefMPC_FLOAT beliefMPC_re12[5];
beliefMPC_FLOAT beliefMPC_beta12[5];
beliefMPC_FLOAT beliefMPC_betacc12[5];
beliefMPC_FLOAT* beliefMPC_dvaff12 = beliefMPC_dv_aff + 65;
beliefMPC_FLOAT* beliefMPC_dvcc12 = beliefMPC_dv_cc + 65;
beliefMPC_FLOAT beliefMPC_V12[35];
beliefMPC_FLOAT beliefMPC_Yd12[15];
beliefMPC_FLOAT beliefMPC_Ld12[15];
beliefMPC_FLOAT beliefMPC_yy12[5];
beliefMPC_FLOAT beliefMPC_bmy12[5];
int beliefMPC_lbIdx12[7] = {0, 1, 2, 3, 4, 5, 6};
beliefMPC_FLOAT* beliefMPC_llb12 = beliefMPC_l + 168;
beliefMPC_FLOAT* beliefMPC_slb12 = beliefMPC_s + 168;
beliefMPC_FLOAT* beliefMPC_llbbyslb12 = beliefMPC_lbys + 168;
beliefMPC_FLOAT beliefMPC_rilb12[7];
beliefMPC_FLOAT* beliefMPC_dllbaff12 = beliefMPC_dl_aff + 168;
beliefMPC_FLOAT* beliefMPC_dslbaff12 = beliefMPC_ds_aff + 168;
beliefMPC_FLOAT* beliefMPC_dllbcc12 = beliefMPC_dl_cc + 168;
beliefMPC_FLOAT* beliefMPC_dslbcc12 = beliefMPC_ds_cc + 168;
beliefMPC_FLOAT* beliefMPC_ccrhsl12 = beliefMPC_ccrhs + 168;
int beliefMPC_ubIdx12[7] = {0, 1, 2, 3, 4, 5, 6};
beliefMPC_FLOAT* beliefMPC_lub12 = beliefMPC_l + 175;
beliefMPC_FLOAT* beliefMPC_sub12 = beliefMPC_s + 175;
beliefMPC_FLOAT* beliefMPC_lubbysub12 = beliefMPC_lbys + 175;
beliefMPC_FLOAT beliefMPC_riub12[7];
beliefMPC_FLOAT* beliefMPC_dlubaff12 = beliefMPC_dl_aff + 175;
beliefMPC_FLOAT* beliefMPC_dsubaff12 = beliefMPC_ds_aff + 175;
beliefMPC_FLOAT* beliefMPC_dlubcc12 = beliefMPC_dl_cc + 175;
beliefMPC_FLOAT* beliefMPC_dsubcc12 = beliefMPC_ds_cc + 175;
beliefMPC_FLOAT* beliefMPC_ccrhsub12 = beliefMPC_ccrhs + 175;
beliefMPC_FLOAT beliefMPC_Phi12[7];
beliefMPC_FLOAT beliefMPC_W12[7];
beliefMPC_FLOAT beliefMPC_Ysd12[25];
beliefMPC_FLOAT beliefMPC_Lsd12[25];
beliefMPC_FLOAT beliefMPC_f13[7] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
beliefMPC_FLOAT* beliefMPC_z13 = beliefMPC_z + 91;
beliefMPC_FLOAT* beliefMPC_dzaff13 = beliefMPC_dz_aff + 91;
beliefMPC_FLOAT* beliefMPC_dzcc13 = beliefMPC_dz_cc + 91;
beliefMPC_FLOAT* beliefMPC_rd13 = beliefMPC_rd + 91;
beliefMPC_FLOAT beliefMPC_Lbyrd13[7];
beliefMPC_FLOAT* beliefMPC_grad_cost13 = beliefMPC_grad_cost + 91;
beliefMPC_FLOAT* beliefMPC_grad_eq13 = beliefMPC_grad_eq + 91;
beliefMPC_FLOAT* beliefMPC_grad_ineq13 = beliefMPC_grad_ineq + 91;
beliefMPC_FLOAT beliefMPC_ctv13[7];
beliefMPC_FLOAT* beliefMPC_v13 = beliefMPC_v + 70;
beliefMPC_FLOAT beliefMPC_re13[5];
beliefMPC_FLOAT beliefMPC_beta13[5];
beliefMPC_FLOAT beliefMPC_betacc13[5];
beliefMPC_FLOAT* beliefMPC_dvaff13 = beliefMPC_dv_aff + 70;
beliefMPC_FLOAT* beliefMPC_dvcc13 = beliefMPC_dv_cc + 70;
beliefMPC_FLOAT beliefMPC_V13[35];
beliefMPC_FLOAT beliefMPC_Yd13[15];
beliefMPC_FLOAT beliefMPC_Ld13[15];
beliefMPC_FLOAT beliefMPC_yy13[5];
beliefMPC_FLOAT beliefMPC_bmy13[5];
int beliefMPC_lbIdx13[7] = {0, 1, 2, 3, 4, 5, 6};
beliefMPC_FLOAT* beliefMPC_llb13 = beliefMPC_l + 182;
beliefMPC_FLOAT* beliefMPC_slb13 = beliefMPC_s + 182;
beliefMPC_FLOAT* beliefMPC_llbbyslb13 = beliefMPC_lbys + 182;
beliefMPC_FLOAT beliefMPC_rilb13[7];
beliefMPC_FLOAT* beliefMPC_dllbaff13 = beliefMPC_dl_aff + 182;
beliefMPC_FLOAT* beliefMPC_dslbaff13 = beliefMPC_ds_aff + 182;
beliefMPC_FLOAT* beliefMPC_dllbcc13 = beliefMPC_dl_cc + 182;
beliefMPC_FLOAT* beliefMPC_dslbcc13 = beliefMPC_ds_cc + 182;
beliefMPC_FLOAT* beliefMPC_ccrhsl13 = beliefMPC_ccrhs + 182;
int beliefMPC_ubIdx13[7] = {0, 1, 2, 3, 4, 5, 6};
beliefMPC_FLOAT* beliefMPC_lub13 = beliefMPC_l + 189;
beliefMPC_FLOAT* beliefMPC_sub13 = beliefMPC_s + 189;
beliefMPC_FLOAT* beliefMPC_lubbysub13 = beliefMPC_lbys + 189;
beliefMPC_FLOAT beliefMPC_riub13[7];
beliefMPC_FLOAT* beliefMPC_dlubaff13 = beliefMPC_dl_aff + 189;
beliefMPC_FLOAT* beliefMPC_dsubaff13 = beliefMPC_ds_aff + 189;
beliefMPC_FLOAT* beliefMPC_dlubcc13 = beliefMPC_dl_cc + 189;
beliefMPC_FLOAT* beliefMPC_dsubcc13 = beliefMPC_ds_cc + 189;
beliefMPC_FLOAT* beliefMPC_ccrhsub13 = beliefMPC_ccrhs + 189;
beliefMPC_FLOAT beliefMPC_Phi13[7];
beliefMPC_FLOAT beliefMPC_W13[7];
beliefMPC_FLOAT beliefMPC_Ysd13[25];
beliefMPC_FLOAT beliefMPC_Lsd13[25];
beliefMPC_FLOAT beliefMPC_H14[5] = {0.0000000000000000E+000, 0.0000000000000000E+000, 2.0000000000000000E+001, 2.0000000000000000E+001, 2.0000000000000000E+001};
beliefMPC_FLOAT beliefMPC_f14[5] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
beliefMPC_FLOAT* beliefMPC_z14 = beliefMPC_z + 98;
beliefMPC_FLOAT* beliefMPC_dzaff14 = beliefMPC_dz_aff + 98;
beliefMPC_FLOAT* beliefMPC_dzcc14 = beliefMPC_dz_cc + 98;
beliefMPC_FLOAT* beliefMPC_rd14 = beliefMPC_rd + 98;
beliefMPC_FLOAT beliefMPC_Lbyrd14[5];
beliefMPC_FLOAT* beliefMPC_grad_cost14 = beliefMPC_grad_cost + 98;
beliefMPC_FLOAT* beliefMPC_grad_eq14 = beliefMPC_grad_eq + 98;
beliefMPC_FLOAT* beliefMPC_grad_ineq14 = beliefMPC_grad_ineq + 98;
beliefMPC_FLOAT beliefMPC_ctv14[5];
int beliefMPC_lbIdx14[5] = {0, 1, 2, 3, 4};
beliefMPC_FLOAT* beliefMPC_llb14 = beliefMPC_l + 196;
beliefMPC_FLOAT* beliefMPC_slb14 = beliefMPC_s + 196;
beliefMPC_FLOAT* beliefMPC_llbbyslb14 = beliefMPC_lbys + 196;
beliefMPC_FLOAT beliefMPC_rilb14[5];
beliefMPC_FLOAT* beliefMPC_dllbaff14 = beliefMPC_dl_aff + 196;
beliefMPC_FLOAT* beliefMPC_dslbaff14 = beliefMPC_ds_aff + 196;
beliefMPC_FLOAT* beliefMPC_dllbcc14 = beliefMPC_dl_cc + 196;
beliefMPC_FLOAT* beliefMPC_dslbcc14 = beliefMPC_ds_cc + 196;
beliefMPC_FLOAT* beliefMPC_ccrhsl14 = beliefMPC_ccrhs + 196;
int beliefMPC_ubIdx14[5] = {0, 1, 2, 3, 4};
beliefMPC_FLOAT* beliefMPC_lub14 = beliefMPC_l + 201;
beliefMPC_FLOAT* beliefMPC_sub14 = beliefMPC_s + 201;
beliefMPC_FLOAT* beliefMPC_lubbysub14 = beliefMPC_lbys + 201;
beliefMPC_FLOAT beliefMPC_riub14[5];
beliefMPC_FLOAT* beliefMPC_dlubaff14 = beliefMPC_dl_aff + 201;
beliefMPC_FLOAT* beliefMPC_dsubaff14 = beliefMPC_ds_aff + 201;
beliefMPC_FLOAT* beliefMPC_dlubcc14 = beliefMPC_dl_cc + 201;
beliefMPC_FLOAT* beliefMPC_dsubcc14 = beliefMPC_ds_cc + 201;
beliefMPC_FLOAT* beliefMPC_ccrhsub14 = beliefMPC_ccrhs + 201;
beliefMPC_FLOAT beliefMPC_Phi14[5];
beliefMPC_FLOAT beliefMPC_D14[5] = {-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000};
beliefMPC_FLOAT beliefMPC_W14[5];
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
beliefMPC_LA_INITIALIZEVECTOR_103(beliefMPC_z, 0);
beliefMPC_LA_INITIALIZEVECTOR_75(beliefMPC_v, 1);
beliefMPC_LA_INITIALIZEVECTOR_206(beliefMPC_l, 1);
beliefMPC_LA_INITIALIZEVECTOR_206(beliefMPC_s, 1);
info->mu = 0;
beliefMPC_LA_DOTACC_206(beliefMPC_l, beliefMPC_s, &info->mu);
info->mu /= 206;
while( 1 ){
info->pobj = 0;
beliefMPC_LA_DIAG_QUADFCN_7(beliefMPC_H00, beliefMPC_f00, beliefMPC_z00, beliefMPC_grad_cost00, &info->pobj);
beliefMPC_LA_DIAG_QUADFCN_7(beliefMPC_H00, beliefMPC_f01, beliefMPC_z01, beliefMPC_grad_cost01, &info->pobj);
beliefMPC_LA_DIAG_QUADFCN_7(beliefMPC_H00, beliefMPC_f02, beliefMPC_z02, beliefMPC_grad_cost02, &info->pobj);
beliefMPC_LA_DIAG_QUADFCN_7(beliefMPC_H00, beliefMPC_f03, beliefMPC_z03, beliefMPC_grad_cost03, &info->pobj);
beliefMPC_LA_DIAG_QUADFCN_7(beliefMPC_H00, beliefMPC_f04, beliefMPC_z04, beliefMPC_grad_cost04, &info->pobj);
beliefMPC_LA_DIAG_QUADFCN_7(beliefMPC_H00, beliefMPC_f05, beliefMPC_z05, beliefMPC_grad_cost05, &info->pobj);
beliefMPC_LA_DIAG_QUADFCN_7(beliefMPC_H00, beliefMPC_f06, beliefMPC_z06, beliefMPC_grad_cost06, &info->pobj);
beliefMPC_LA_DIAG_QUADFCN_7(beliefMPC_H00, beliefMPC_f07, beliefMPC_z07, beliefMPC_grad_cost07, &info->pobj);
beliefMPC_LA_DIAG_QUADFCN_7(beliefMPC_H00, beliefMPC_f08, beliefMPC_z08, beliefMPC_grad_cost08, &info->pobj);
beliefMPC_LA_DIAG_QUADFCN_7(beliefMPC_H00, beliefMPC_f09, beliefMPC_z09, beliefMPC_grad_cost09, &info->pobj);
beliefMPC_LA_DIAG_QUADFCN_7(beliefMPC_H00, beliefMPC_f10, beliefMPC_z10, beliefMPC_grad_cost10, &info->pobj);
beliefMPC_LA_DIAG_QUADFCN_7(beliefMPC_H00, beliefMPC_f11, beliefMPC_z11, beliefMPC_grad_cost11, &info->pobj);
beliefMPC_LA_DIAG_QUADFCN_7(beliefMPC_H00, beliefMPC_f12, beliefMPC_z12, beliefMPC_grad_cost12, &info->pobj);
beliefMPC_LA_DIAG_QUADFCN_7(beliefMPC_H00, beliefMPC_f13, beliefMPC_z13, beliefMPC_grad_cost13, &info->pobj);
beliefMPC_LA_DIAG_QUADFCN_5(beliefMPC_H14, beliefMPC_f14, beliefMPC_z14, beliefMPC_grad_cost14, &info->pobj);
info->res_eq = 0;
info->dgap = 0;
beliefMPC_LA_DENSE_MVMSUB3_10_7_7(params->C1, beliefMPC_z00, beliefMPC_D01, beliefMPC_z01, params->e1, beliefMPC_v00, beliefMPC_re00, &info->dgap, &info->res_eq);
beliefMPC_LA_DENSE_DIAGZERO_MVMSUB3_5_7_7(params->C2, beliefMPC_z01, beliefMPC_D02, beliefMPC_z02, params->e2, beliefMPC_v01, beliefMPC_re01, &info->dgap, &info->res_eq);
beliefMPC_LA_DENSE_DIAGZERO_MVMSUB3_5_7_7(params->C3, beliefMPC_z02, beliefMPC_D02, beliefMPC_z03, params->e3, beliefMPC_v02, beliefMPC_re02, &info->dgap, &info->res_eq);
beliefMPC_LA_DENSE_DIAGZERO_MVMSUB3_5_7_7(params->C4, beliefMPC_z03, beliefMPC_D02, beliefMPC_z04, params->e4, beliefMPC_v03, beliefMPC_re03, &info->dgap, &info->res_eq);
beliefMPC_LA_DENSE_DIAGZERO_MVMSUB3_5_7_7(params->C5, beliefMPC_z04, beliefMPC_D02, beliefMPC_z05, params->e5, beliefMPC_v04, beliefMPC_re04, &info->dgap, &info->res_eq);
beliefMPC_LA_DENSE_DIAGZERO_MVMSUB3_5_7_7(params->C6, beliefMPC_z05, beliefMPC_D02, beliefMPC_z06, params->e6, beliefMPC_v05, beliefMPC_re05, &info->dgap, &info->res_eq);
beliefMPC_LA_DENSE_DIAGZERO_MVMSUB3_5_7_7(params->C7, beliefMPC_z06, beliefMPC_D02, beliefMPC_z07, params->e7, beliefMPC_v06, beliefMPC_re06, &info->dgap, &info->res_eq);
beliefMPC_LA_DENSE_DIAGZERO_MVMSUB3_5_7_7(params->C8, beliefMPC_z07, beliefMPC_D02, beliefMPC_z08, params->e8, beliefMPC_v07, beliefMPC_re07, &info->dgap, &info->res_eq);
beliefMPC_LA_DENSE_DIAGZERO_MVMSUB3_5_7_7(params->C9, beliefMPC_z08, beliefMPC_D02, beliefMPC_z09, params->e9, beliefMPC_v08, beliefMPC_re08, &info->dgap, &info->res_eq);
beliefMPC_LA_DENSE_DIAGZERO_MVMSUB3_5_7_7(params->C10, beliefMPC_z09, beliefMPC_D02, beliefMPC_z10, params->e10, beliefMPC_v09, beliefMPC_re09, &info->dgap, &info->res_eq);
beliefMPC_LA_DENSE_DIAGZERO_MVMSUB3_5_7_7(params->C11, beliefMPC_z10, beliefMPC_D02, beliefMPC_z11, params->e11, beliefMPC_v10, beliefMPC_re10, &info->dgap, &info->res_eq);
beliefMPC_LA_DENSE_DIAGZERO_MVMSUB3_5_7_7(params->C12, beliefMPC_z11, beliefMPC_D02, beliefMPC_z12, params->e12, beliefMPC_v11, beliefMPC_re11, &info->dgap, &info->res_eq);
beliefMPC_LA_DENSE_DIAGZERO_MVMSUB3_5_7_7(params->C13, beliefMPC_z12, beliefMPC_D02, beliefMPC_z13, params->e13, beliefMPC_v12, beliefMPC_re12, &info->dgap, &info->res_eq);
beliefMPC_LA_DENSE_DIAGZERO_MVMSUB3_5_7_5(params->C14, beliefMPC_z13, beliefMPC_D14, beliefMPC_z14, params->e14, beliefMPC_v13, beliefMPC_re13, &info->dgap, &info->res_eq);
beliefMPC_LA_DENSE_MTVM_10_7(params->C1, beliefMPC_v00, beliefMPC_grad_eq00);
beliefMPC_LA_DENSE_MTVM2_5_7_10(params->C2, beliefMPC_v01, beliefMPC_D01, beliefMPC_v00, beliefMPC_grad_eq01);
beliefMPC_LA_DENSE_DIAGZERO_MTVM2_5_7_5(params->C3, beliefMPC_v02, beliefMPC_D02, beliefMPC_v01, beliefMPC_grad_eq02);
beliefMPC_LA_DENSE_DIAGZERO_MTVM2_5_7_5(params->C4, beliefMPC_v03, beliefMPC_D02, beliefMPC_v02, beliefMPC_grad_eq03);
beliefMPC_LA_DENSE_DIAGZERO_MTVM2_5_7_5(params->C5, beliefMPC_v04, beliefMPC_D02, beliefMPC_v03, beliefMPC_grad_eq04);
beliefMPC_LA_DENSE_DIAGZERO_MTVM2_5_7_5(params->C6, beliefMPC_v05, beliefMPC_D02, beliefMPC_v04, beliefMPC_grad_eq05);
beliefMPC_LA_DENSE_DIAGZERO_MTVM2_5_7_5(params->C7, beliefMPC_v06, beliefMPC_D02, beliefMPC_v05, beliefMPC_grad_eq06);
beliefMPC_LA_DENSE_DIAGZERO_MTVM2_5_7_5(params->C8, beliefMPC_v07, beliefMPC_D02, beliefMPC_v06, beliefMPC_grad_eq07);
beliefMPC_LA_DENSE_DIAGZERO_MTVM2_5_7_5(params->C9, beliefMPC_v08, beliefMPC_D02, beliefMPC_v07, beliefMPC_grad_eq08);
beliefMPC_LA_DENSE_DIAGZERO_MTVM2_5_7_5(params->C10, beliefMPC_v09, beliefMPC_D02, beliefMPC_v08, beliefMPC_grad_eq09);
beliefMPC_LA_DENSE_DIAGZERO_MTVM2_5_7_5(params->C11, beliefMPC_v10, beliefMPC_D02, beliefMPC_v09, beliefMPC_grad_eq10);
beliefMPC_LA_DENSE_DIAGZERO_MTVM2_5_7_5(params->C12, beliefMPC_v11, beliefMPC_D02, beliefMPC_v10, beliefMPC_grad_eq11);
beliefMPC_LA_DENSE_DIAGZERO_MTVM2_5_7_5(params->C13, beliefMPC_v12, beliefMPC_D02, beliefMPC_v11, beliefMPC_grad_eq12);
beliefMPC_LA_DENSE_DIAGZERO_MTVM2_5_7_5(params->C14, beliefMPC_v13, beliefMPC_D02, beliefMPC_v12, beliefMPC_grad_eq13);
beliefMPC_LA_DIAGZERO_MTVM_5_5(beliefMPC_D14, beliefMPC_v13, beliefMPC_grad_eq14);
info->res_ineq = 0;
beliefMPC_LA_VSUBADD3_7(params->lb1, beliefMPC_z00, beliefMPC_lbIdx00, beliefMPC_llb00, beliefMPC_slb00, beliefMPC_rilb00, &info->dgap, &info->res_ineq);
beliefMPC_LA_VSUBADD2_7(beliefMPC_z00, beliefMPC_ubIdx00, params->ub1, beliefMPC_lub00, beliefMPC_sub00, beliefMPC_riub00, &info->dgap, &info->res_ineq);
beliefMPC_LA_VSUBADD3_7(params->lb2, beliefMPC_z01, beliefMPC_lbIdx01, beliefMPC_llb01, beliefMPC_slb01, beliefMPC_rilb01, &info->dgap, &info->res_ineq);
beliefMPC_LA_VSUBADD2_7(beliefMPC_z01, beliefMPC_ubIdx01, params->ub2, beliefMPC_lub01, beliefMPC_sub01, beliefMPC_riub01, &info->dgap, &info->res_ineq);
beliefMPC_LA_VSUBADD3_7(params->lb3, beliefMPC_z02, beliefMPC_lbIdx02, beliefMPC_llb02, beliefMPC_slb02, beliefMPC_rilb02, &info->dgap, &info->res_ineq);
beliefMPC_LA_VSUBADD2_7(beliefMPC_z02, beliefMPC_ubIdx02, params->ub3, beliefMPC_lub02, beliefMPC_sub02, beliefMPC_riub02, &info->dgap, &info->res_ineq);
beliefMPC_LA_VSUBADD3_7(params->lb4, beliefMPC_z03, beliefMPC_lbIdx03, beliefMPC_llb03, beliefMPC_slb03, beliefMPC_rilb03, &info->dgap, &info->res_ineq);
beliefMPC_LA_VSUBADD2_7(beliefMPC_z03, beliefMPC_ubIdx03, params->ub4, beliefMPC_lub03, beliefMPC_sub03, beliefMPC_riub03, &info->dgap, &info->res_ineq);
beliefMPC_LA_VSUBADD3_7(params->lb5, beliefMPC_z04, beliefMPC_lbIdx04, beliefMPC_llb04, beliefMPC_slb04, beliefMPC_rilb04, &info->dgap, &info->res_ineq);
beliefMPC_LA_VSUBADD2_7(beliefMPC_z04, beliefMPC_ubIdx04, params->ub5, beliefMPC_lub04, beliefMPC_sub04, beliefMPC_riub04, &info->dgap, &info->res_ineq);
beliefMPC_LA_VSUBADD3_7(params->lb6, beliefMPC_z05, beliefMPC_lbIdx05, beliefMPC_llb05, beliefMPC_slb05, beliefMPC_rilb05, &info->dgap, &info->res_ineq);
beliefMPC_LA_VSUBADD2_7(beliefMPC_z05, beliefMPC_ubIdx05, params->ub6, beliefMPC_lub05, beliefMPC_sub05, beliefMPC_riub05, &info->dgap, &info->res_ineq);
beliefMPC_LA_VSUBADD3_7(params->lb7, beliefMPC_z06, beliefMPC_lbIdx06, beliefMPC_llb06, beliefMPC_slb06, beliefMPC_rilb06, &info->dgap, &info->res_ineq);
beliefMPC_LA_VSUBADD2_7(beliefMPC_z06, beliefMPC_ubIdx06, params->ub7, beliefMPC_lub06, beliefMPC_sub06, beliefMPC_riub06, &info->dgap, &info->res_ineq);
beliefMPC_LA_VSUBADD3_7(params->lb8, beliefMPC_z07, beliefMPC_lbIdx07, beliefMPC_llb07, beliefMPC_slb07, beliefMPC_rilb07, &info->dgap, &info->res_ineq);
beliefMPC_LA_VSUBADD2_7(beliefMPC_z07, beliefMPC_ubIdx07, params->ub8, beliefMPC_lub07, beliefMPC_sub07, beliefMPC_riub07, &info->dgap, &info->res_ineq);
beliefMPC_LA_VSUBADD3_7(params->lb9, beliefMPC_z08, beliefMPC_lbIdx08, beliefMPC_llb08, beliefMPC_slb08, beliefMPC_rilb08, &info->dgap, &info->res_ineq);
beliefMPC_LA_VSUBADD2_7(beliefMPC_z08, beliefMPC_ubIdx08, params->ub9, beliefMPC_lub08, beliefMPC_sub08, beliefMPC_riub08, &info->dgap, &info->res_ineq);
beliefMPC_LA_VSUBADD3_7(params->lb10, beliefMPC_z09, beliefMPC_lbIdx09, beliefMPC_llb09, beliefMPC_slb09, beliefMPC_rilb09, &info->dgap, &info->res_ineq);
beliefMPC_LA_VSUBADD2_7(beliefMPC_z09, beliefMPC_ubIdx09, params->ub10, beliefMPC_lub09, beliefMPC_sub09, beliefMPC_riub09, &info->dgap, &info->res_ineq);
beliefMPC_LA_VSUBADD3_7(params->lb11, beliefMPC_z10, beliefMPC_lbIdx10, beliefMPC_llb10, beliefMPC_slb10, beliefMPC_rilb10, &info->dgap, &info->res_ineq);
beliefMPC_LA_VSUBADD2_7(beliefMPC_z10, beliefMPC_ubIdx10, params->ub11, beliefMPC_lub10, beliefMPC_sub10, beliefMPC_riub10, &info->dgap, &info->res_ineq);
beliefMPC_LA_VSUBADD3_7(params->lb12, beliefMPC_z11, beliefMPC_lbIdx11, beliefMPC_llb11, beliefMPC_slb11, beliefMPC_rilb11, &info->dgap, &info->res_ineq);
beliefMPC_LA_VSUBADD2_7(beliefMPC_z11, beliefMPC_ubIdx11, params->ub12, beliefMPC_lub11, beliefMPC_sub11, beliefMPC_riub11, &info->dgap, &info->res_ineq);
beliefMPC_LA_VSUBADD3_7(params->lb13, beliefMPC_z12, beliefMPC_lbIdx12, beliefMPC_llb12, beliefMPC_slb12, beliefMPC_rilb12, &info->dgap, &info->res_ineq);
beliefMPC_LA_VSUBADD2_7(beliefMPC_z12, beliefMPC_ubIdx12, params->ub13, beliefMPC_lub12, beliefMPC_sub12, beliefMPC_riub12, &info->dgap, &info->res_ineq);
beliefMPC_LA_VSUBADD3_7(params->lb14, beliefMPC_z13, beliefMPC_lbIdx13, beliefMPC_llb13, beliefMPC_slb13, beliefMPC_rilb13, &info->dgap, &info->res_ineq);
beliefMPC_LA_VSUBADD2_7(beliefMPC_z13, beliefMPC_ubIdx13, params->ub14, beliefMPC_lub13, beliefMPC_sub13, beliefMPC_riub13, &info->dgap, &info->res_ineq);
beliefMPC_LA_VSUBADD3_5(params->lb15, beliefMPC_z14, beliefMPC_lbIdx14, beliefMPC_llb14, beliefMPC_slb14, beliefMPC_rilb14, &info->dgap, &info->res_ineq);
beliefMPC_LA_VSUBADD2_5(beliefMPC_z14, beliefMPC_ubIdx14, params->ub15, beliefMPC_lub14, beliefMPC_sub14, beliefMPC_riub14, &info->dgap, &info->res_ineq);
beliefMPC_LA_INEQ_B_GRAD_7_7_7(beliefMPC_lub00, beliefMPC_sub00, beliefMPC_riub00, beliefMPC_llb00, beliefMPC_slb00, beliefMPC_rilb00, beliefMPC_lbIdx00, beliefMPC_ubIdx00, beliefMPC_grad_ineq00, beliefMPC_lubbysub00, beliefMPC_llbbyslb00);
beliefMPC_LA_INEQ_B_GRAD_7_7_7(beliefMPC_lub01, beliefMPC_sub01, beliefMPC_riub01, beliefMPC_llb01, beliefMPC_slb01, beliefMPC_rilb01, beliefMPC_lbIdx01, beliefMPC_ubIdx01, beliefMPC_grad_ineq01, beliefMPC_lubbysub01, beliefMPC_llbbyslb01);
beliefMPC_LA_INEQ_B_GRAD_7_7_7(beliefMPC_lub02, beliefMPC_sub02, beliefMPC_riub02, beliefMPC_llb02, beliefMPC_slb02, beliefMPC_rilb02, beliefMPC_lbIdx02, beliefMPC_ubIdx02, beliefMPC_grad_ineq02, beliefMPC_lubbysub02, beliefMPC_llbbyslb02);
beliefMPC_LA_INEQ_B_GRAD_7_7_7(beliefMPC_lub03, beliefMPC_sub03, beliefMPC_riub03, beliefMPC_llb03, beliefMPC_slb03, beliefMPC_rilb03, beliefMPC_lbIdx03, beliefMPC_ubIdx03, beliefMPC_grad_ineq03, beliefMPC_lubbysub03, beliefMPC_llbbyslb03);
beliefMPC_LA_INEQ_B_GRAD_7_7_7(beliefMPC_lub04, beliefMPC_sub04, beliefMPC_riub04, beliefMPC_llb04, beliefMPC_slb04, beliefMPC_rilb04, beliefMPC_lbIdx04, beliefMPC_ubIdx04, beliefMPC_grad_ineq04, beliefMPC_lubbysub04, beliefMPC_llbbyslb04);
beliefMPC_LA_INEQ_B_GRAD_7_7_7(beliefMPC_lub05, beliefMPC_sub05, beliefMPC_riub05, beliefMPC_llb05, beliefMPC_slb05, beliefMPC_rilb05, beliefMPC_lbIdx05, beliefMPC_ubIdx05, beliefMPC_grad_ineq05, beliefMPC_lubbysub05, beliefMPC_llbbyslb05);
beliefMPC_LA_INEQ_B_GRAD_7_7_7(beliefMPC_lub06, beliefMPC_sub06, beliefMPC_riub06, beliefMPC_llb06, beliefMPC_slb06, beliefMPC_rilb06, beliefMPC_lbIdx06, beliefMPC_ubIdx06, beliefMPC_grad_ineq06, beliefMPC_lubbysub06, beliefMPC_llbbyslb06);
beliefMPC_LA_INEQ_B_GRAD_7_7_7(beliefMPC_lub07, beliefMPC_sub07, beliefMPC_riub07, beliefMPC_llb07, beliefMPC_slb07, beliefMPC_rilb07, beliefMPC_lbIdx07, beliefMPC_ubIdx07, beliefMPC_grad_ineq07, beliefMPC_lubbysub07, beliefMPC_llbbyslb07);
beliefMPC_LA_INEQ_B_GRAD_7_7_7(beliefMPC_lub08, beliefMPC_sub08, beliefMPC_riub08, beliefMPC_llb08, beliefMPC_slb08, beliefMPC_rilb08, beliefMPC_lbIdx08, beliefMPC_ubIdx08, beliefMPC_grad_ineq08, beliefMPC_lubbysub08, beliefMPC_llbbyslb08);
beliefMPC_LA_INEQ_B_GRAD_7_7_7(beliefMPC_lub09, beliefMPC_sub09, beliefMPC_riub09, beliefMPC_llb09, beliefMPC_slb09, beliefMPC_rilb09, beliefMPC_lbIdx09, beliefMPC_ubIdx09, beliefMPC_grad_ineq09, beliefMPC_lubbysub09, beliefMPC_llbbyslb09);
beliefMPC_LA_INEQ_B_GRAD_7_7_7(beliefMPC_lub10, beliefMPC_sub10, beliefMPC_riub10, beliefMPC_llb10, beliefMPC_slb10, beliefMPC_rilb10, beliefMPC_lbIdx10, beliefMPC_ubIdx10, beliefMPC_grad_ineq10, beliefMPC_lubbysub10, beliefMPC_llbbyslb10);
beliefMPC_LA_INEQ_B_GRAD_7_7_7(beliefMPC_lub11, beliefMPC_sub11, beliefMPC_riub11, beliefMPC_llb11, beliefMPC_slb11, beliefMPC_rilb11, beliefMPC_lbIdx11, beliefMPC_ubIdx11, beliefMPC_grad_ineq11, beliefMPC_lubbysub11, beliefMPC_llbbyslb11);
beliefMPC_LA_INEQ_B_GRAD_7_7_7(beliefMPC_lub12, beliefMPC_sub12, beliefMPC_riub12, beliefMPC_llb12, beliefMPC_slb12, beliefMPC_rilb12, beliefMPC_lbIdx12, beliefMPC_ubIdx12, beliefMPC_grad_ineq12, beliefMPC_lubbysub12, beliefMPC_llbbyslb12);
beliefMPC_LA_INEQ_B_GRAD_7_7_7(beliefMPC_lub13, beliefMPC_sub13, beliefMPC_riub13, beliefMPC_llb13, beliefMPC_slb13, beliefMPC_rilb13, beliefMPC_lbIdx13, beliefMPC_ubIdx13, beliefMPC_grad_ineq13, beliefMPC_lubbysub13, beliefMPC_llbbyslb13);
beliefMPC_LA_INEQ_B_GRAD_5_5_5(beliefMPC_lub14, beliefMPC_sub14, beliefMPC_riub14, beliefMPC_llb14, beliefMPC_slb14, beliefMPC_rilb14, beliefMPC_lbIdx14, beliefMPC_ubIdx14, beliefMPC_grad_ineq14, beliefMPC_lubbysub14, beliefMPC_llbbyslb14);
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
beliefMPC_LA_VVADD3_103(beliefMPC_grad_cost, beliefMPC_grad_eq, beliefMPC_grad_ineq, beliefMPC_rd);
beliefMPC_LA_DIAG_CHOL_ONELOOP_LBUB_7_7_7(beliefMPC_H00, beliefMPC_llbbyslb00, beliefMPC_lbIdx00, beliefMPC_lubbysub00, beliefMPC_ubIdx00, beliefMPC_Phi00);
beliefMPC_LA_DIAG_MATRIXFORWARDSUB_10_7(beliefMPC_Phi00, params->C1, beliefMPC_V00);
beliefMPC_LA_DIAG_FORWARDSUB_7(beliefMPC_Phi00, beliefMPC_rd00, beliefMPC_Lbyrd00);
beliefMPC_LA_DIAG_CHOL_ONELOOP_LBUB_7_7_7(beliefMPC_H00, beliefMPC_llbbyslb01, beliefMPC_lbIdx01, beliefMPC_lubbysub01, beliefMPC_ubIdx01, beliefMPC_Phi01);
beliefMPC_LA_DIAG_MATRIXFORWARDSUB_5_7(beliefMPC_Phi01, params->C2, beliefMPC_V01);
beliefMPC_LA_DIAG_MATRIXFORWARDSUB_10_7(beliefMPC_Phi01, beliefMPC_D01, beliefMPC_W01);
beliefMPC_LA_DENSE_MMTM_10_7_5(beliefMPC_W01, beliefMPC_V01, beliefMPC_Ysd01);
beliefMPC_LA_DIAG_FORWARDSUB_7(beliefMPC_Phi01, beliefMPC_rd01, beliefMPC_Lbyrd01);
beliefMPC_LA_DIAG_CHOL_ONELOOP_LBUB_7_7_7(beliefMPC_H00, beliefMPC_llbbyslb02, beliefMPC_lbIdx02, beliefMPC_lubbysub02, beliefMPC_ubIdx02, beliefMPC_Phi02);
beliefMPC_LA_DIAG_MATRIXFORWARDSUB_5_7(beliefMPC_Phi02, params->C3, beliefMPC_V02);
beliefMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_5_7(beliefMPC_Phi02, beliefMPC_D02, beliefMPC_W02);
beliefMPC_LA_DENSE_DIAGZERO_MMTM_5_7_5(beliefMPC_W02, beliefMPC_V02, beliefMPC_Ysd02);
beliefMPC_LA_DIAG_FORWARDSUB_7(beliefMPC_Phi02, beliefMPC_rd02, beliefMPC_Lbyrd02);
beliefMPC_LA_DIAG_CHOL_ONELOOP_LBUB_7_7_7(beliefMPC_H00, beliefMPC_llbbyslb03, beliefMPC_lbIdx03, beliefMPC_lubbysub03, beliefMPC_ubIdx03, beliefMPC_Phi03);
beliefMPC_LA_DIAG_MATRIXFORWARDSUB_5_7(beliefMPC_Phi03, params->C4, beliefMPC_V03);
beliefMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_5_7(beliefMPC_Phi03, beliefMPC_D02, beliefMPC_W03);
beliefMPC_LA_DENSE_DIAGZERO_MMTM_5_7_5(beliefMPC_W03, beliefMPC_V03, beliefMPC_Ysd03);
beliefMPC_LA_DIAG_FORWARDSUB_7(beliefMPC_Phi03, beliefMPC_rd03, beliefMPC_Lbyrd03);
beliefMPC_LA_DIAG_CHOL_ONELOOP_LBUB_7_7_7(beliefMPC_H00, beliefMPC_llbbyslb04, beliefMPC_lbIdx04, beliefMPC_lubbysub04, beliefMPC_ubIdx04, beliefMPC_Phi04);
beliefMPC_LA_DIAG_MATRIXFORWARDSUB_5_7(beliefMPC_Phi04, params->C5, beliefMPC_V04);
beliefMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_5_7(beliefMPC_Phi04, beliefMPC_D02, beliefMPC_W04);
beliefMPC_LA_DENSE_DIAGZERO_MMTM_5_7_5(beliefMPC_W04, beliefMPC_V04, beliefMPC_Ysd04);
beliefMPC_LA_DIAG_FORWARDSUB_7(beliefMPC_Phi04, beliefMPC_rd04, beliefMPC_Lbyrd04);
beliefMPC_LA_DIAG_CHOL_ONELOOP_LBUB_7_7_7(beliefMPC_H00, beliefMPC_llbbyslb05, beliefMPC_lbIdx05, beliefMPC_lubbysub05, beliefMPC_ubIdx05, beliefMPC_Phi05);
beliefMPC_LA_DIAG_MATRIXFORWARDSUB_5_7(beliefMPC_Phi05, params->C6, beliefMPC_V05);
beliefMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_5_7(beliefMPC_Phi05, beliefMPC_D02, beliefMPC_W05);
beliefMPC_LA_DENSE_DIAGZERO_MMTM_5_7_5(beliefMPC_W05, beliefMPC_V05, beliefMPC_Ysd05);
beliefMPC_LA_DIAG_FORWARDSUB_7(beliefMPC_Phi05, beliefMPC_rd05, beliefMPC_Lbyrd05);
beliefMPC_LA_DIAG_CHOL_ONELOOP_LBUB_7_7_7(beliefMPC_H00, beliefMPC_llbbyslb06, beliefMPC_lbIdx06, beliefMPC_lubbysub06, beliefMPC_ubIdx06, beliefMPC_Phi06);
beliefMPC_LA_DIAG_MATRIXFORWARDSUB_5_7(beliefMPC_Phi06, params->C7, beliefMPC_V06);
beliefMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_5_7(beliefMPC_Phi06, beliefMPC_D02, beliefMPC_W06);
beliefMPC_LA_DENSE_DIAGZERO_MMTM_5_7_5(beliefMPC_W06, beliefMPC_V06, beliefMPC_Ysd06);
beliefMPC_LA_DIAG_FORWARDSUB_7(beliefMPC_Phi06, beliefMPC_rd06, beliefMPC_Lbyrd06);
beliefMPC_LA_DIAG_CHOL_ONELOOP_LBUB_7_7_7(beliefMPC_H00, beliefMPC_llbbyslb07, beliefMPC_lbIdx07, beliefMPC_lubbysub07, beliefMPC_ubIdx07, beliefMPC_Phi07);
beliefMPC_LA_DIAG_MATRIXFORWARDSUB_5_7(beliefMPC_Phi07, params->C8, beliefMPC_V07);
beliefMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_5_7(beliefMPC_Phi07, beliefMPC_D02, beliefMPC_W07);
beliefMPC_LA_DENSE_DIAGZERO_MMTM_5_7_5(beliefMPC_W07, beliefMPC_V07, beliefMPC_Ysd07);
beliefMPC_LA_DIAG_FORWARDSUB_7(beliefMPC_Phi07, beliefMPC_rd07, beliefMPC_Lbyrd07);
beliefMPC_LA_DIAG_CHOL_ONELOOP_LBUB_7_7_7(beliefMPC_H00, beliefMPC_llbbyslb08, beliefMPC_lbIdx08, beliefMPC_lubbysub08, beliefMPC_ubIdx08, beliefMPC_Phi08);
beliefMPC_LA_DIAG_MATRIXFORWARDSUB_5_7(beliefMPC_Phi08, params->C9, beliefMPC_V08);
beliefMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_5_7(beliefMPC_Phi08, beliefMPC_D02, beliefMPC_W08);
beliefMPC_LA_DENSE_DIAGZERO_MMTM_5_7_5(beliefMPC_W08, beliefMPC_V08, beliefMPC_Ysd08);
beliefMPC_LA_DIAG_FORWARDSUB_7(beliefMPC_Phi08, beliefMPC_rd08, beliefMPC_Lbyrd08);
beliefMPC_LA_DIAG_CHOL_ONELOOP_LBUB_7_7_7(beliefMPC_H00, beliefMPC_llbbyslb09, beliefMPC_lbIdx09, beliefMPC_lubbysub09, beliefMPC_ubIdx09, beliefMPC_Phi09);
beliefMPC_LA_DIAG_MATRIXFORWARDSUB_5_7(beliefMPC_Phi09, params->C10, beliefMPC_V09);
beliefMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_5_7(beliefMPC_Phi09, beliefMPC_D02, beliefMPC_W09);
beliefMPC_LA_DENSE_DIAGZERO_MMTM_5_7_5(beliefMPC_W09, beliefMPC_V09, beliefMPC_Ysd09);
beliefMPC_LA_DIAG_FORWARDSUB_7(beliefMPC_Phi09, beliefMPC_rd09, beliefMPC_Lbyrd09);
beliefMPC_LA_DIAG_CHOL_ONELOOP_LBUB_7_7_7(beliefMPC_H00, beliefMPC_llbbyslb10, beliefMPC_lbIdx10, beliefMPC_lubbysub10, beliefMPC_ubIdx10, beliefMPC_Phi10);
beliefMPC_LA_DIAG_MATRIXFORWARDSUB_5_7(beliefMPC_Phi10, params->C11, beliefMPC_V10);
beliefMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_5_7(beliefMPC_Phi10, beliefMPC_D02, beliefMPC_W10);
beliefMPC_LA_DENSE_DIAGZERO_MMTM_5_7_5(beliefMPC_W10, beliefMPC_V10, beliefMPC_Ysd10);
beliefMPC_LA_DIAG_FORWARDSUB_7(beliefMPC_Phi10, beliefMPC_rd10, beliefMPC_Lbyrd10);
beliefMPC_LA_DIAG_CHOL_ONELOOP_LBUB_7_7_7(beliefMPC_H00, beliefMPC_llbbyslb11, beliefMPC_lbIdx11, beliefMPC_lubbysub11, beliefMPC_ubIdx11, beliefMPC_Phi11);
beliefMPC_LA_DIAG_MATRIXFORWARDSUB_5_7(beliefMPC_Phi11, params->C12, beliefMPC_V11);
beliefMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_5_7(beliefMPC_Phi11, beliefMPC_D02, beliefMPC_W11);
beliefMPC_LA_DENSE_DIAGZERO_MMTM_5_7_5(beliefMPC_W11, beliefMPC_V11, beliefMPC_Ysd11);
beliefMPC_LA_DIAG_FORWARDSUB_7(beliefMPC_Phi11, beliefMPC_rd11, beliefMPC_Lbyrd11);
beliefMPC_LA_DIAG_CHOL_ONELOOP_LBUB_7_7_7(beliefMPC_H00, beliefMPC_llbbyslb12, beliefMPC_lbIdx12, beliefMPC_lubbysub12, beliefMPC_ubIdx12, beliefMPC_Phi12);
beliefMPC_LA_DIAG_MATRIXFORWARDSUB_5_7(beliefMPC_Phi12, params->C13, beliefMPC_V12);
beliefMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_5_7(beliefMPC_Phi12, beliefMPC_D02, beliefMPC_W12);
beliefMPC_LA_DENSE_DIAGZERO_MMTM_5_7_5(beliefMPC_W12, beliefMPC_V12, beliefMPC_Ysd12);
beliefMPC_LA_DIAG_FORWARDSUB_7(beliefMPC_Phi12, beliefMPC_rd12, beliefMPC_Lbyrd12);
beliefMPC_LA_DIAG_CHOL_ONELOOP_LBUB_7_7_7(beliefMPC_H00, beliefMPC_llbbyslb13, beliefMPC_lbIdx13, beliefMPC_lubbysub13, beliefMPC_ubIdx13, beliefMPC_Phi13);
beliefMPC_LA_DIAG_MATRIXFORWARDSUB_5_7(beliefMPC_Phi13, params->C14, beliefMPC_V13);
beliefMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_5_7(beliefMPC_Phi13, beliefMPC_D02, beliefMPC_W13);
beliefMPC_LA_DENSE_DIAGZERO_MMTM_5_7_5(beliefMPC_W13, beliefMPC_V13, beliefMPC_Ysd13);
beliefMPC_LA_DIAG_FORWARDSUB_7(beliefMPC_Phi13, beliefMPC_rd13, beliefMPC_Lbyrd13);
beliefMPC_LA_DIAG_CHOL_ONELOOP_LBUB_5_5_5(beliefMPC_H14, beliefMPC_llbbyslb14, beliefMPC_lbIdx14, beliefMPC_lubbysub14, beliefMPC_ubIdx14, beliefMPC_Phi14);
beliefMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_5_5(beliefMPC_Phi14, beliefMPC_D14, beliefMPC_W14);
beliefMPC_LA_DIAG_FORWARDSUB_5(beliefMPC_Phi14, beliefMPC_rd14, beliefMPC_Lbyrd14);
beliefMPC_LA_DENSE_MMT2_10_7_7(beliefMPC_V00, beliefMPC_W01, beliefMPC_Yd00);
beliefMPC_LA_DENSE_MVMSUB2_10_7_7(beliefMPC_V00, beliefMPC_Lbyrd00, beliefMPC_W01, beliefMPC_Lbyrd01, beliefMPC_re00, beliefMPC_beta00);
beliefMPC_LA_DENSE_DIAGZERO_MMT2_5_7_7(beliefMPC_V01, beliefMPC_W02, beliefMPC_Yd01);
beliefMPC_LA_DENSE_DIAGZERO_2MVMSUB2_5_7_7(beliefMPC_V01, beliefMPC_Lbyrd01, beliefMPC_W02, beliefMPC_Lbyrd02, beliefMPC_re01, beliefMPC_beta01);
beliefMPC_LA_DENSE_DIAGZERO_MMT2_5_7_7(beliefMPC_V02, beliefMPC_W03, beliefMPC_Yd02);
beliefMPC_LA_DENSE_DIAGZERO_2MVMSUB2_5_7_7(beliefMPC_V02, beliefMPC_Lbyrd02, beliefMPC_W03, beliefMPC_Lbyrd03, beliefMPC_re02, beliefMPC_beta02);
beliefMPC_LA_DENSE_DIAGZERO_MMT2_5_7_7(beliefMPC_V03, beliefMPC_W04, beliefMPC_Yd03);
beliefMPC_LA_DENSE_DIAGZERO_2MVMSUB2_5_7_7(beliefMPC_V03, beliefMPC_Lbyrd03, beliefMPC_W04, beliefMPC_Lbyrd04, beliefMPC_re03, beliefMPC_beta03);
beliefMPC_LA_DENSE_DIAGZERO_MMT2_5_7_7(beliefMPC_V04, beliefMPC_W05, beliefMPC_Yd04);
beliefMPC_LA_DENSE_DIAGZERO_2MVMSUB2_5_7_7(beliefMPC_V04, beliefMPC_Lbyrd04, beliefMPC_W05, beliefMPC_Lbyrd05, beliefMPC_re04, beliefMPC_beta04);
beliefMPC_LA_DENSE_DIAGZERO_MMT2_5_7_7(beliefMPC_V05, beliefMPC_W06, beliefMPC_Yd05);
beliefMPC_LA_DENSE_DIAGZERO_2MVMSUB2_5_7_7(beliefMPC_V05, beliefMPC_Lbyrd05, beliefMPC_W06, beliefMPC_Lbyrd06, beliefMPC_re05, beliefMPC_beta05);
beliefMPC_LA_DENSE_DIAGZERO_MMT2_5_7_7(beliefMPC_V06, beliefMPC_W07, beliefMPC_Yd06);
beliefMPC_LA_DENSE_DIAGZERO_2MVMSUB2_5_7_7(beliefMPC_V06, beliefMPC_Lbyrd06, beliefMPC_W07, beliefMPC_Lbyrd07, beliefMPC_re06, beliefMPC_beta06);
beliefMPC_LA_DENSE_DIAGZERO_MMT2_5_7_7(beliefMPC_V07, beliefMPC_W08, beliefMPC_Yd07);
beliefMPC_LA_DENSE_DIAGZERO_2MVMSUB2_5_7_7(beliefMPC_V07, beliefMPC_Lbyrd07, beliefMPC_W08, beliefMPC_Lbyrd08, beliefMPC_re07, beliefMPC_beta07);
beliefMPC_LA_DENSE_DIAGZERO_MMT2_5_7_7(beliefMPC_V08, beliefMPC_W09, beliefMPC_Yd08);
beliefMPC_LA_DENSE_DIAGZERO_2MVMSUB2_5_7_7(beliefMPC_V08, beliefMPC_Lbyrd08, beliefMPC_W09, beliefMPC_Lbyrd09, beliefMPC_re08, beliefMPC_beta08);
beliefMPC_LA_DENSE_DIAGZERO_MMT2_5_7_7(beliefMPC_V09, beliefMPC_W10, beliefMPC_Yd09);
beliefMPC_LA_DENSE_DIAGZERO_2MVMSUB2_5_7_7(beliefMPC_V09, beliefMPC_Lbyrd09, beliefMPC_W10, beliefMPC_Lbyrd10, beliefMPC_re09, beliefMPC_beta09);
beliefMPC_LA_DENSE_DIAGZERO_MMT2_5_7_7(beliefMPC_V10, beliefMPC_W11, beliefMPC_Yd10);
beliefMPC_LA_DENSE_DIAGZERO_2MVMSUB2_5_7_7(beliefMPC_V10, beliefMPC_Lbyrd10, beliefMPC_W11, beliefMPC_Lbyrd11, beliefMPC_re10, beliefMPC_beta10);
beliefMPC_LA_DENSE_DIAGZERO_MMT2_5_7_7(beliefMPC_V11, beliefMPC_W12, beliefMPC_Yd11);
beliefMPC_LA_DENSE_DIAGZERO_2MVMSUB2_5_7_7(beliefMPC_V11, beliefMPC_Lbyrd11, beliefMPC_W12, beliefMPC_Lbyrd12, beliefMPC_re11, beliefMPC_beta11);
beliefMPC_LA_DENSE_DIAGZERO_MMT2_5_7_7(beliefMPC_V12, beliefMPC_W13, beliefMPC_Yd12);
beliefMPC_LA_DENSE_DIAGZERO_2MVMSUB2_5_7_7(beliefMPC_V12, beliefMPC_Lbyrd12, beliefMPC_W13, beliefMPC_Lbyrd13, beliefMPC_re12, beliefMPC_beta12);
beliefMPC_LA_DENSE_DIAGZERO_MMT2_5_7_5(beliefMPC_V13, beliefMPC_W14, beliefMPC_Yd13);
beliefMPC_LA_DENSE_DIAGZERO_2MVMSUB2_5_7_5(beliefMPC_V13, beliefMPC_Lbyrd13, beliefMPC_W14, beliefMPC_Lbyrd14, beliefMPC_re13, beliefMPC_beta13);
beliefMPC_LA_DENSE_CHOL_10(beliefMPC_Yd00, beliefMPC_Ld00);
beliefMPC_LA_DENSE_FORWARDSUB_10(beliefMPC_Ld00, beliefMPC_beta00, beliefMPC_yy00);
beliefMPC_LA_DENSE_MATRIXTFORWARDSUB_5_10(beliefMPC_Ld00, beliefMPC_Ysd01, beliefMPC_Lsd01);
beliefMPC_LA_DENSE_MMTSUB_5_10(beliefMPC_Lsd01, beliefMPC_Yd01);
beliefMPC_LA_DENSE_CHOL_5(beliefMPC_Yd01, beliefMPC_Ld01);
beliefMPC_LA_DENSE_MVMSUB1_5_10(beliefMPC_Lsd01, beliefMPC_yy00, beliefMPC_beta01, beliefMPC_bmy01);
beliefMPC_LA_DENSE_FORWARDSUB_5(beliefMPC_Ld01, beliefMPC_bmy01, beliefMPC_yy01);
beliefMPC_LA_DENSE_MATRIXTFORWARDSUB_5_5(beliefMPC_Ld01, beliefMPC_Ysd02, beliefMPC_Lsd02);
beliefMPC_LA_DENSE_MMTSUB_5_5(beliefMPC_Lsd02, beliefMPC_Yd02);
beliefMPC_LA_DENSE_CHOL_5(beliefMPC_Yd02, beliefMPC_Ld02);
beliefMPC_LA_DENSE_MVMSUB1_5_5(beliefMPC_Lsd02, beliefMPC_yy01, beliefMPC_beta02, beliefMPC_bmy02);
beliefMPC_LA_DENSE_FORWARDSUB_5(beliefMPC_Ld02, beliefMPC_bmy02, beliefMPC_yy02);
beliefMPC_LA_DENSE_MATRIXTFORWARDSUB_5_5(beliefMPC_Ld02, beliefMPC_Ysd03, beliefMPC_Lsd03);
beliefMPC_LA_DENSE_MMTSUB_5_5(beliefMPC_Lsd03, beliefMPC_Yd03);
beliefMPC_LA_DENSE_CHOL_5(beliefMPC_Yd03, beliefMPC_Ld03);
beliefMPC_LA_DENSE_MVMSUB1_5_5(beliefMPC_Lsd03, beliefMPC_yy02, beliefMPC_beta03, beliefMPC_bmy03);
beliefMPC_LA_DENSE_FORWARDSUB_5(beliefMPC_Ld03, beliefMPC_bmy03, beliefMPC_yy03);
beliefMPC_LA_DENSE_MATRIXTFORWARDSUB_5_5(beliefMPC_Ld03, beliefMPC_Ysd04, beliefMPC_Lsd04);
beliefMPC_LA_DENSE_MMTSUB_5_5(beliefMPC_Lsd04, beliefMPC_Yd04);
beliefMPC_LA_DENSE_CHOL_5(beliefMPC_Yd04, beliefMPC_Ld04);
beliefMPC_LA_DENSE_MVMSUB1_5_5(beliefMPC_Lsd04, beliefMPC_yy03, beliefMPC_beta04, beliefMPC_bmy04);
beliefMPC_LA_DENSE_FORWARDSUB_5(beliefMPC_Ld04, beliefMPC_bmy04, beliefMPC_yy04);
beliefMPC_LA_DENSE_MATRIXTFORWARDSUB_5_5(beliefMPC_Ld04, beliefMPC_Ysd05, beliefMPC_Lsd05);
beliefMPC_LA_DENSE_MMTSUB_5_5(beliefMPC_Lsd05, beliefMPC_Yd05);
beliefMPC_LA_DENSE_CHOL_5(beliefMPC_Yd05, beliefMPC_Ld05);
beliefMPC_LA_DENSE_MVMSUB1_5_5(beliefMPC_Lsd05, beliefMPC_yy04, beliefMPC_beta05, beliefMPC_bmy05);
beliefMPC_LA_DENSE_FORWARDSUB_5(beliefMPC_Ld05, beliefMPC_bmy05, beliefMPC_yy05);
beliefMPC_LA_DENSE_MATRIXTFORWARDSUB_5_5(beliefMPC_Ld05, beliefMPC_Ysd06, beliefMPC_Lsd06);
beliefMPC_LA_DENSE_MMTSUB_5_5(beliefMPC_Lsd06, beliefMPC_Yd06);
beliefMPC_LA_DENSE_CHOL_5(beliefMPC_Yd06, beliefMPC_Ld06);
beliefMPC_LA_DENSE_MVMSUB1_5_5(beliefMPC_Lsd06, beliefMPC_yy05, beliefMPC_beta06, beliefMPC_bmy06);
beliefMPC_LA_DENSE_FORWARDSUB_5(beliefMPC_Ld06, beliefMPC_bmy06, beliefMPC_yy06);
beliefMPC_LA_DENSE_MATRIXTFORWARDSUB_5_5(beliefMPC_Ld06, beliefMPC_Ysd07, beliefMPC_Lsd07);
beliefMPC_LA_DENSE_MMTSUB_5_5(beliefMPC_Lsd07, beliefMPC_Yd07);
beliefMPC_LA_DENSE_CHOL_5(beliefMPC_Yd07, beliefMPC_Ld07);
beliefMPC_LA_DENSE_MVMSUB1_5_5(beliefMPC_Lsd07, beliefMPC_yy06, beliefMPC_beta07, beliefMPC_bmy07);
beliefMPC_LA_DENSE_FORWARDSUB_5(beliefMPC_Ld07, beliefMPC_bmy07, beliefMPC_yy07);
beliefMPC_LA_DENSE_MATRIXTFORWARDSUB_5_5(beliefMPC_Ld07, beliefMPC_Ysd08, beliefMPC_Lsd08);
beliefMPC_LA_DENSE_MMTSUB_5_5(beliefMPC_Lsd08, beliefMPC_Yd08);
beliefMPC_LA_DENSE_CHOL_5(beliefMPC_Yd08, beliefMPC_Ld08);
beliefMPC_LA_DENSE_MVMSUB1_5_5(beliefMPC_Lsd08, beliefMPC_yy07, beliefMPC_beta08, beliefMPC_bmy08);
beliefMPC_LA_DENSE_FORWARDSUB_5(beliefMPC_Ld08, beliefMPC_bmy08, beliefMPC_yy08);
beliefMPC_LA_DENSE_MATRIXTFORWARDSUB_5_5(beliefMPC_Ld08, beliefMPC_Ysd09, beliefMPC_Lsd09);
beliefMPC_LA_DENSE_MMTSUB_5_5(beliefMPC_Lsd09, beliefMPC_Yd09);
beliefMPC_LA_DENSE_CHOL_5(beliefMPC_Yd09, beliefMPC_Ld09);
beliefMPC_LA_DENSE_MVMSUB1_5_5(beliefMPC_Lsd09, beliefMPC_yy08, beliefMPC_beta09, beliefMPC_bmy09);
beliefMPC_LA_DENSE_FORWARDSUB_5(beliefMPC_Ld09, beliefMPC_bmy09, beliefMPC_yy09);
beliefMPC_LA_DENSE_MATRIXTFORWARDSUB_5_5(beliefMPC_Ld09, beliefMPC_Ysd10, beliefMPC_Lsd10);
beliefMPC_LA_DENSE_MMTSUB_5_5(beliefMPC_Lsd10, beliefMPC_Yd10);
beliefMPC_LA_DENSE_CHOL_5(beliefMPC_Yd10, beliefMPC_Ld10);
beliefMPC_LA_DENSE_MVMSUB1_5_5(beliefMPC_Lsd10, beliefMPC_yy09, beliefMPC_beta10, beliefMPC_bmy10);
beliefMPC_LA_DENSE_FORWARDSUB_5(beliefMPC_Ld10, beliefMPC_bmy10, beliefMPC_yy10);
beliefMPC_LA_DENSE_MATRIXTFORWARDSUB_5_5(beliefMPC_Ld10, beliefMPC_Ysd11, beliefMPC_Lsd11);
beliefMPC_LA_DENSE_MMTSUB_5_5(beliefMPC_Lsd11, beliefMPC_Yd11);
beliefMPC_LA_DENSE_CHOL_5(beliefMPC_Yd11, beliefMPC_Ld11);
beliefMPC_LA_DENSE_MVMSUB1_5_5(beliefMPC_Lsd11, beliefMPC_yy10, beliefMPC_beta11, beliefMPC_bmy11);
beliefMPC_LA_DENSE_FORWARDSUB_5(beliefMPC_Ld11, beliefMPC_bmy11, beliefMPC_yy11);
beliefMPC_LA_DENSE_MATRIXTFORWARDSUB_5_5(beliefMPC_Ld11, beliefMPC_Ysd12, beliefMPC_Lsd12);
beliefMPC_LA_DENSE_MMTSUB_5_5(beliefMPC_Lsd12, beliefMPC_Yd12);
beliefMPC_LA_DENSE_CHOL_5(beliefMPC_Yd12, beliefMPC_Ld12);
beliefMPC_LA_DENSE_MVMSUB1_5_5(beliefMPC_Lsd12, beliefMPC_yy11, beliefMPC_beta12, beliefMPC_bmy12);
beliefMPC_LA_DENSE_FORWARDSUB_5(beliefMPC_Ld12, beliefMPC_bmy12, beliefMPC_yy12);
beliefMPC_LA_DENSE_MATRIXTFORWARDSUB_5_5(beliefMPC_Ld12, beliefMPC_Ysd13, beliefMPC_Lsd13);
beliefMPC_LA_DENSE_MMTSUB_5_5(beliefMPC_Lsd13, beliefMPC_Yd13);
beliefMPC_LA_DENSE_CHOL_5(beliefMPC_Yd13, beliefMPC_Ld13);
beliefMPC_LA_DENSE_MVMSUB1_5_5(beliefMPC_Lsd13, beliefMPC_yy12, beliefMPC_beta13, beliefMPC_bmy13);
beliefMPC_LA_DENSE_FORWARDSUB_5(beliefMPC_Ld13, beliefMPC_bmy13, beliefMPC_yy13);
beliefMPC_LA_DENSE_BACKWARDSUB_5(beliefMPC_Ld13, beliefMPC_yy13, beliefMPC_dvaff13);
beliefMPC_LA_DENSE_MTVMSUB_5_5(beliefMPC_Lsd13, beliefMPC_dvaff13, beliefMPC_yy12, beliefMPC_bmy12);
beliefMPC_LA_DENSE_BACKWARDSUB_5(beliefMPC_Ld12, beliefMPC_bmy12, beliefMPC_dvaff12);
beliefMPC_LA_DENSE_MTVMSUB_5_5(beliefMPC_Lsd12, beliefMPC_dvaff12, beliefMPC_yy11, beliefMPC_bmy11);
beliefMPC_LA_DENSE_BACKWARDSUB_5(beliefMPC_Ld11, beliefMPC_bmy11, beliefMPC_dvaff11);
beliefMPC_LA_DENSE_MTVMSUB_5_5(beliefMPC_Lsd11, beliefMPC_dvaff11, beliefMPC_yy10, beliefMPC_bmy10);
beliefMPC_LA_DENSE_BACKWARDSUB_5(beliefMPC_Ld10, beliefMPC_bmy10, beliefMPC_dvaff10);
beliefMPC_LA_DENSE_MTVMSUB_5_5(beliefMPC_Lsd10, beliefMPC_dvaff10, beliefMPC_yy09, beliefMPC_bmy09);
beliefMPC_LA_DENSE_BACKWARDSUB_5(beliefMPC_Ld09, beliefMPC_bmy09, beliefMPC_dvaff09);
beliefMPC_LA_DENSE_MTVMSUB_5_5(beliefMPC_Lsd09, beliefMPC_dvaff09, beliefMPC_yy08, beliefMPC_bmy08);
beliefMPC_LA_DENSE_BACKWARDSUB_5(beliefMPC_Ld08, beliefMPC_bmy08, beliefMPC_dvaff08);
beliefMPC_LA_DENSE_MTVMSUB_5_5(beliefMPC_Lsd08, beliefMPC_dvaff08, beliefMPC_yy07, beliefMPC_bmy07);
beliefMPC_LA_DENSE_BACKWARDSUB_5(beliefMPC_Ld07, beliefMPC_bmy07, beliefMPC_dvaff07);
beliefMPC_LA_DENSE_MTVMSUB_5_5(beliefMPC_Lsd07, beliefMPC_dvaff07, beliefMPC_yy06, beliefMPC_bmy06);
beliefMPC_LA_DENSE_BACKWARDSUB_5(beliefMPC_Ld06, beliefMPC_bmy06, beliefMPC_dvaff06);
beliefMPC_LA_DENSE_MTVMSUB_5_5(beliefMPC_Lsd06, beliefMPC_dvaff06, beliefMPC_yy05, beliefMPC_bmy05);
beliefMPC_LA_DENSE_BACKWARDSUB_5(beliefMPC_Ld05, beliefMPC_bmy05, beliefMPC_dvaff05);
beliefMPC_LA_DENSE_MTVMSUB_5_5(beliefMPC_Lsd05, beliefMPC_dvaff05, beliefMPC_yy04, beliefMPC_bmy04);
beliefMPC_LA_DENSE_BACKWARDSUB_5(beliefMPC_Ld04, beliefMPC_bmy04, beliefMPC_dvaff04);
beliefMPC_LA_DENSE_MTVMSUB_5_5(beliefMPC_Lsd04, beliefMPC_dvaff04, beliefMPC_yy03, beliefMPC_bmy03);
beliefMPC_LA_DENSE_BACKWARDSUB_5(beliefMPC_Ld03, beliefMPC_bmy03, beliefMPC_dvaff03);
beliefMPC_LA_DENSE_MTVMSUB_5_5(beliefMPC_Lsd03, beliefMPC_dvaff03, beliefMPC_yy02, beliefMPC_bmy02);
beliefMPC_LA_DENSE_BACKWARDSUB_5(beliefMPC_Ld02, beliefMPC_bmy02, beliefMPC_dvaff02);
beliefMPC_LA_DENSE_MTVMSUB_5_5(beliefMPC_Lsd02, beliefMPC_dvaff02, beliefMPC_yy01, beliefMPC_bmy01);
beliefMPC_LA_DENSE_BACKWARDSUB_5(beliefMPC_Ld01, beliefMPC_bmy01, beliefMPC_dvaff01);
beliefMPC_LA_DENSE_MTVMSUB_5_10(beliefMPC_Lsd01, beliefMPC_dvaff01, beliefMPC_yy00, beliefMPC_bmy00);
beliefMPC_LA_DENSE_BACKWARDSUB_10(beliefMPC_Ld00, beliefMPC_bmy00, beliefMPC_dvaff00);
beliefMPC_LA_DENSE_MTVM_10_7(params->C1, beliefMPC_dvaff00, beliefMPC_grad_eq00);
beliefMPC_LA_DENSE_MTVM2_5_7_10(params->C2, beliefMPC_dvaff01, beliefMPC_D01, beliefMPC_dvaff00, beliefMPC_grad_eq01);
beliefMPC_LA_DENSE_DIAGZERO_MTVM2_5_7_5(params->C3, beliefMPC_dvaff02, beliefMPC_D02, beliefMPC_dvaff01, beliefMPC_grad_eq02);
beliefMPC_LA_DENSE_DIAGZERO_MTVM2_5_7_5(params->C4, beliefMPC_dvaff03, beliefMPC_D02, beliefMPC_dvaff02, beliefMPC_grad_eq03);
beliefMPC_LA_DENSE_DIAGZERO_MTVM2_5_7_5(params->C5, beliefMPC_dvaff04, beliefMPC_D02, beliefMPC_dvaff03, beliefMPC_grad_eq04);
beliefMPC_LA_DENSE_DIAGZERO_MTVM2_5_7_5(params->C6, beliefMPC_dvaff05, beliefMPC_D02, beliefMPC_dvaff04, beliefMPC_grad_eq05);
beliefMPC_LA_DENSE_DIAGZERO_MTVM2_5_7_5(params->C7, beliefMPC_dvaff06, beliefMPC_D02, beliefMPC_dvaff05, beliefMPC_grad_eq06);
beliefMPC_LA_DENSE_DIAGZERO_MTVM2_5_7_5(params->C8, beliefMPC_dvaff07, beliefMPC_D02, beliefMPC_dvaff06, beliefMPC_grad_eq07);
beliefMPC_LA_DENSE_DIAGZERO_MTVM2_5_7_5(params->C9, beliefMPC_dvaff08, beliefMPC_D02, beliefMPC_dvaff07, beliefMPC_grad_eq08);
beliefMPC_LA_DENSE_DIAGZERO_MTVM2_5_7_5(params->C10, beliefMPC_dvaff09, beliefMPC_D02, beliefMPC_dvaff08, beliefMPC_grad_eq09);
beliefMPC_LA_DENSE_DIAGZERO_MTVM2_5_7_5(params->C11, beliefMPC_dvaff10, beliefMPC_D02, beliefMPC_dvaff09, beliefMPC_grad_eq10);
beliefMPC_LA_DENSE_DIAGZERO_MTVM2_5_7_5(params->C12, beliefMPC_dvaff11, beliefMPC_D02, beliefMPC_dvaff10, beliefMPC_grad_eq11);
beliefMPC_LA_DENSE_DIAGZERO_MTVM2_5_7_5(params->C13, beliefMPC_dvaff12, beliefMPC_D02, beliefMPC_dvaff11, beliefMPC_grad_eq12);
beliefMPC_LA_DENSE_DIAGZERO_MTVM2_5_7_5(params->C14, beliefMPC_dvaff13, beliefMPC_D02, beliefMPC_dvaff12, beliefMPC_grad_eq13);
beliefMPC_LA_DIAGZERO_MTVM_5_5(beliefMPC_D14, beliefMPC_dvaff13, beliefMPC_grad_eq14);
beliefMPC_LA_VSUB2_103(beliefMPC_rd, beliefMPC_grad_eq, beliefMPC_rd);
beliefMPC_LA_DIAG_FORWARDBACKWARDSUB_7(beliefMPC_Phi00, beliefMPC_rd00, beliefMPC_dzaff00);
beliefMPC_LA_DIAG_FORWARDBACKWARDSUB_7(beliefMPC_Phi01, beliefMPC_rd01, beliefMPC_dzaff01);
beliefMPC_LA_DIAG_FORWARDBACKWARDSUB_7(beliefMPC_Phi02, beliefMPC_rd02, beliefMPC_dzaff02);
beliefMPC_LA_DIAG_FORWARDBACKWARDSUB_7(beliefMPC_Phi03, beliefMPC_rd03, beliefMPC_dzaff03);
beliefMPC_LA_DIAG_FORWARDBACKWARDSUB_7(beliefMPC_Phi04, beliefMPC_rd04, beliefMPC_dzaff04);
beliefMPC_LA_DIAG_FORWARDBACKWARDSUB_7(beliefMPC_Phi05, beliefMPC_rd05, beliefMPC_dzaff05);
beliefMPC_LA_DIAG_FORWARDBACKWARDSUB_7(beliefMPC_Phi06, beliefMPC_rd06, beliefMPC_dzaff06);
beliefMPC_LA_DIAG_FORWARDBACKWARDSUB_7(beliefMPC_Phi07, beliefMPC_rd07, beliefMPC_dzaff07);
beliefMPC_LA_DIAG_FORWARDBACKWARDSUB_7(beliefMPC_Phi08, beliefMPC_rd08, beliefMPC_dzaff08);
beliefMPC_LA_DIAG_FORWARDBACKWARDSUB_7(beliefMPC_Phi09, beliefMPC_rd09, beliefMPC_dzaff09);
beliefMPC_LA_DIAG_FORWARDBACKWARDSUB_7(beliefMPC_Phi10, beliefMPC_rd10, beliefMPC_dzaff10);
beliefMPC_LA_DIAG_FORWARDBACKWARDSUB_7(beliefMPC_Phi11, beliefMPC_rd11, beliefMPC_dzaff11);
beliefMPC_LA_DIAG_FORWARDBACKWARDSUB_7(beliefMPC_Phi12, beliefMPC_rd12, beliefMPC_dzaff12);
beliefMPC_LA_DIAG_FORWARDBACKWARDSUB_7(beliefMPC_Phi13, beliefMPC_rd13, beliefMPC_dzaff13);
beliefMPC_LA_DIAG_FORWARDBACKWARDSUB_5(beliefMPC_Phi14, beliefMPC_rd14, beliefMPC_dzaff14);
beliefMPC_LA_VSUB_INDEXED_7(beliefMPC_dzaff00, beliefMPC_lbIdx00, beliefMPC_rilb00, beliefMPC_dslbaff00);
beliefMPC_LA_VSUB3_7(beliefMPC_llbbyslb00, beliefMPC_dslbaff00, beliefMPC_llb00, beliefMPC_dllbaff00);
beliefMPC_LA_VSUB2_INDEXED_7(beliefMPC_riub00, beliefMPC_dzaff00, beliefMPC_ubIdx00, beliefMPC_dsubaff00);
beliefMPC_LA_VSUB3_7(beliefMPC_lubbysub00, beliefMPC_dsubaff00, beliefMPC_lub00, beliefMPC_dlubaff00);
beliefMPC_LA_VSUB_INDEXED_7(beliefMPC_dzaff01, beliefMPC_lbIdx01, beliefMPC_rilb01, beliefMPC_dslbaff01);
beliefMPC_LA_VSUB3_7(beliefMPC_llbbyslb01, beliefMPC_dslbaff01, beliefMPC_llb01, beliefMPC_dllbaff01);
beliefMPC_LA_VSUB2_INDEXED_7(beliefMPC_riub01, beliefMPC_dzaff01, beliefMPC_ubIdx01, beliefMPC_dsubaff01);
beliefMPC_LA_VSUB3_7(beliefMPC_lubbysub01, beliefMPC_dsubaff01, beliefMPC_lub01, beliefMPC_dlubaff01);
beliefMPC_LA_VSUB_INDEXED_7(beliefMPC_dzaff02, beliefMPC_lbIdx02, beliefMPC_rilb02, beliefMPC_dslbaff02);
beliefMPC_LA_VSUB3_7(beliefMPC_llbbyslb02, beliefMPC_dslbaff02, beliefMPC_llb02, beliefMPC_dllbaff02);
beliefMPC_LA_VSUB2_INDEXED_7(beliefMPC_riub02, beliefMPC_dzaff02, beliefMPC_ubIdx02, beliefMPC_dsubaff02);
beliefMPC_LA_VSUB3_7(beliefMPC_lubbysub02, beliefMPC_dsubaff02, beliefMPC_lub02, beliefMPC_dlubaff02);
beliefMPC_LA_VSUB_INDEXED_7(beliefMPC_dzaff03, beliefMPC_lbIdx03, beliefMPC_rilb03, beliefMPC_dslbaff03);
beliefMPC_LA_VSUB3_7(beliefMPC_llbbyslb03, beliefMPC_dslbaff03, beliefMPC_llb03, beliefMPC_dllbaff03);
beliefMPC_LA_VSUB2_INDEXED_7(beliefMPC_riub03, beliefMPC_dzaff03, beliefMPC_ubIdx03, beliefMPC_dsubaff03);
beliefMPC_LA_VSUB3_7(beliefMPC_lubbysub03, beliefMPC_dsubaff03, beliefMPC_lub03, beliefMPC_dlubaff03);
beliefMPC_LA_VSUB_INDEXED_7(beliefMPC_dzaff04, beliefMPC_lbIdx04, beliefMPC_rilb04, beliefMPC_dslbaff04);
beliefMPC_LA_VSUB3_7(beliefMPC_llbbyslb04, beliefMPC_dslbaff04, beliefMPC_llb04, beliefMPC_dllbaff04);
beliefMPC_LA_VSUB2_INDEXED_7(beliefMPC_riub04, beliefMPC_dzaff04, beliefMPC_ubIdx04, beliefMPC_dsubaff04);
beliefMPC_LA_VSUB3_7(beliefMPC_lubbysub04, beliefMPC_dsubaff04, beliefMPC_lub04, beliefMPC_dlubaff04);
beliefMPC_LA_VSUB_INDEXED_7(beliefMPC_dzaff05, beliefMPC_lbIdx05, beliefMPC_rilb05, beliefMPC_dslbaff05);
beliefMPC_LA_VSUB3_7(beliefMPC_llbbyslb05, beliefMPC_dslbaff05, beliefMPC_llb05, beliefMPC_dllbaff05);
beliefMPC_LA_VSUB2_INDEXED_7(beliefMPC_riub05, beliefMPC_dzaff05, beliefMPC_ubIdx05, beliefMPC_dsubaff05);
beliefMPC_LA_VSUB3_7(beliefMPC_lubbysub05, beliefMPC_dsubaff05, beliefMPC_lub05, beliefMPC_dlubaff05);
beliefMPC_LA_VSUB_INDEXED_7(beliefMPC_dzaff06, beliefMPC_lbIdx06, beliefMPC_rilb06, beliefMPC_dslbaff06);
beliefMPC_LA_VSUB3_7(beliefMPC_llbbyslb06, beliefMPC_dslbaff06, beliefMPC_llb06, beliefMPC_dllbaff06);
beliefMPC_LA_VSUB2_INDEXED_7(beliefMPC_riub06, beliefMPC_dzaff06, beliefMPC_ubIdx06, beliefMPC_dsubaff06);
beliefMPC_LA_VSUB3_7(beliefMPC_lubbysub06, beliefMPC_dsubaff06, beliefMPC_lub06, beliefMPC_dlubaff06);
beliefMPC_LA_VSUB_INDEXED_7(beliefMPC_dzaff07, beliefMPC_lbIdx07, beliefMPC_rilb07, beliefMPC_dslbaff07);
beliefMPC_LA_VSUB3_7(beliefMPC_llbbyslb07, beliefMPC_dslbaff07, beliefMPC_llb07, beliefMPC_dllbaff07);
beliefMPC_LA_VSUB2_INDEXED_7(beliefMPC_riub07, beliefMPC_dzaff07, beliefMPC_ubIdx07, beliefMPC_dsubaff07);
beliefMPC_LA_VSUB3_7(beliefMPC_lubbysub07, beliefMPC_dsubaff07, beliefMPC_lub07, beliefMPC_dlubaff07);
beliefMPC_LA_VSUB_INDEXED_7(beliefMPC_dzaff08, beliefMPC_lbIdx08, beliefMPC_rilb08, beliefMPC_dslbaff08);
beliefMPC_LA_VSUB3_7(beliefMPC_llbbyslb08, beliefMPC_dslbaff08, beliefMPC_llb08, beliefMPC_dllbaff08);
beliefMPC_LA_VSUB2_INDEXED_7(beliefMPC_riub08, beliefMPC_dzaff08, beliefMPC_ubIdx08, beliefMPC_dsubaff08);
beliefMPC_LA_VSUB3_7(beliefMPC_lubbysub08, beliefMPC_dsubaff08, beliefMPC_lub08, beliefMPC_dlubaff08);
beliefMPC_LA_VSUB_INDEXED_7(beliefMPC_dzaff09, beliefMPC_lbIdx09, beliefMPC_rilb09, beliefMPC_dslbaff09);
beliefMPC_LA_VSUB3_7(beliefMPC_llbbyslb09, beliefMPC_dslbaff09, beliefMPC_llb09, beliefMPC_dllbaff09);
beliefMPC_LA_VSUB2_INDEXED_7(beliefMPC_riub09, beliefMPC_dzaff09, beliefMPC_ubIdx09, beliefMPC_dsubaff09);
beliefMPC_LA_VSUB3_7(beliefMPC_lubbysub09, beliefMPC_dsubaff09, beliefMPC_lub09, beliefMPC_dlubaff09);
beliefMPC_LA_VSUB_INDEXED_7(beliefMPC_dzaff10, beliefMPC_lbIdx10, beliefMPC_rilb10, beliefMPC_dslbaff10);
beliefMPC_LA_VSUB3_7(beliefMPC_llbbyslb10, beliefMPC_dslbaff10, beliefMPC_llb10, beliefMPC_dllbaff10);
beliefMPC_LA_VSUB2_INDEXED_7(beliefMPC_riub10, beliefMPC_dzaff10, beliefMPC_ubIdx10, beliefMPC_dsubaff10);
beliefMPC_LA_VSUB3_7(beliefMPC_lubbysub10, beliefMPC_dsubaff10, beliefMPC_lub10, beliefMPC_dlubaff10);
beliefMPC_LA_VSUB_INDEXED_7(beliefMPC_dzaff11, beliefMPC_lbIdx11, beliefMPC_rilb11, beliefMPC_dslbaff11);
beliefMPC_LA_VSUB3_7(beliefMPC_llbbyslb11, beliefMPC_dslbaff11, beliefMPC_llb11, beliefMPC_dllbaff11);
beliefMPC_LA_VSUB2_INDEXED_7(beliefMPC_riub11, beliefMPC_dzaff11, beliefMPC_ubIdx11, beliefMPC_dsubaff11);
beliefMPC_LA_VSUB3_7(beliefMPC_lubbysub11, beliefMPC_dsubaff11, beliefMPC_lub11, beliefMPC_dlubaff11);
beliefMPC_LA_VSUB_INDEXED_7(beliefMPC_dzaff12, beliefMPC_lbIdx12, beliefMPC_rilb12, beliefMPC_dslbaff12);
beliefMPC_LA_VSUB3_7(beliefMPC_llbbyslb12, beliefMPC_dslbaff12, beliefMPC_llb12, beliefMPC_dllbaff12);
beliefMPC_LA_VSUB2_INDEXED_7(beliefMPC_riub12, beliefMPC_dzaff12, beliefMPC_ubIdx12, beliefMPC_dsubaff12);
beliefMPC_LA_VSUB3_7(beliefMPC_lubbysub12, beliefMPC_dsubaff12, beliefMPC_lub12, beliefMPC_dlubaff12);
beliefMPC_LA_VSUB_INDEXED_7(beliefMPC_dzaff13, beliefMPC_lbIdx13, beliefMPC_rilb13, beliefMPC_dslbaff13);
beliefMPC_LA_VSUB3_7(beliefMPC_llbbyslb13, beliefMPC_dslbaff13, beliefMPC_llb13, beliefMPC_dllbaff13);
beliefMPC_LA_VSUB2_INDEXED_7(beliefMPC_riub13, beliefMPC_dzaff13, beliefMPC_ubIdx13, beliefMPC_dsubaff13);
beliefMPC_LA_VSUB3_7(beliefMPC_lubbysub13, beliefMPC_dsubaff13, beliefMPC_lub13, beliefMPC_dlubaff13);
beliefMPC_LA_VSUB_INDEXED_5(beliefMPC_dzaff14, beliefMPC_lbIdx14, beliefMPC_rilb14, beliefMPC_dslbaff14);
beliefMPC_LA_VSUB3_5(beliefMPC_llbbyslb14, beliefMPC_dslbaff14, beliefMPC_llb14, beliefMPC_dllbaff14);
beliefMPC_LA_VSUB2_INDEXED_5(beliefMPC_riub14, beliefMPC_dzaff14, beliefMPC_ubIdx14, beliefMPC_dsubaff14);
beliefMPC_LA_VSUB3_5(beliefMPC_lubbysub14, beliefMPC_dsubaff14, beliefMPC_lub14, beliefMPC_dlubaff14);
info->lsit_aff = beliefMPC_LINESEARCH_BACKTRACKING_AFFINE(beliefMPC_l, beliefMPC_s, beliefMPC_dl_aff, beliefMPC_ds_aff, &info->step_aff, &info->mu_aff);
if( info->lsit_aff == beliefMPC_NOPROGRESS ){
exitcode = beliefMPC_NOPROGRESS; break;
}
sigma_3rdroot = info->mu_aff / info->mu;
info->sigma = sigma_3rdroot*sigma_3rdroot*sigma_3rdroot;
musigma = info->mu * info->sigma;
beliefMPC_LA_VSUB5_206(beliefMPC_ds_aff, beliefMPC_dl_aff, musigma, beliefMPC_ccrhs);
beliefMPC_LA_VSUB6_INDEXED_7_7_7(beliefMPC_ccrhsub00, beliefMPC_sub00, beliefMPC_ubIdx00, beliefMPC_ccrhsl00, beliefMPC_slb00, beliefMPC_lbIdx00, beliefMPC_rd00);
beliefMPC_LA_VSUB6_INDEXED_7_7_7(beliefMPC_ccrhsub01, beliefMPC_sub01, beliefMPC_ubIdx01, beliefMPC_ccrhsl01, beliefMPC_slb01, beliefMPC_lbIdx01, beliefMPC_rd01);
beliefMPC_LA_DIAG_FORWARDSUB_7(beliefMPC_Phi00, beliefMPC_rd00, beliefMPC_Lbyrd00);
beliefMPC_LA_DIAG_FORWARDSUB_7(beliefMPC_Phi01, beliefMPC_rd01, beliefMPC_Lbyrd01);
beliefMPC_LA_DENSE_2MVMADD_10_7_7(beliefMPC_V00, beliefMPC_Lbyrd00, beliefMPC_W01, beliefMPC_Lbyrd01, beliefMPC_beta00);
beliefMPC_LA_DENSE_FORWARDSUB_10(beliefMPC_Ld00, beliefMPC_beta00, beliefMPC_yy00);
beliefMPC_LA_VSUB6_INDEXED_7_7_7(beliefMPC_ccrhsub02, beliefMPC_sub02, beliefMPC_ubIdx02, beliefMPC_ccrhsl02, beliefMPC_slb02, beliefMPC_lbIdx02, beliefMPC_rd02);
beliefMPC_LA_DIAG_FORWARDSUB_7(beliefMPC_Phi02, beliefMPC_rd02, beliefMPC_Lbyrd02);
beliefMPC_LA_DENSE_DIAGZERO_2MVMADD_5_7_7(beliefMPC_V01, beliefMPC_Lbyrd01, beliefMPC_W02, beliefMPC_Lbyrd02, beliefMPC_beta01);
beliefMPC_LA_DENSE_MVMSUB1_5_10(beliefMPC_Lsd01, beliefMPC_yy00, beliefMPC_beta01, beliefMPC_bmy01);
beliefMPC_LA_DENSE_FORWARDSUB_5(beliefMPC_Ld01, beliefMPC_bmy01, beliefMPC_yy01);
beliefMPC_LA_VSUB6_INDEXED_7_7_7(beliefMPC_ccrhsub03, beliefMPC_sub03, beliefMPC_ubIdx03, beliefMPC_ccrhsl03, beliefMPC_slb03, beliefMPC_lbIdx03, beliefMPC_rd03);
beliefMPC_LA_DIAG_FORWARDSUB_7(beliefMPC_Phi03, beliefMPC_rd03, beliefMPC_Lbyrd03);
beliefMPC_LA_DENSE_DIAGZERO_2MVMADD_5_7_7(beliefMPC_V02, beliefMPC_Lbyrd02, beliefMPC_W03, beliefMPC_Lbyrd03, beliefMPC_beta02);
beliefMPC_LA_DENSE_MVMSUB1_5_5(beliefMPC_Lsd02, beliefMPC_yy01, beliefMPC_beta02, beliefMPC_bmy02);
beliefMPC_LA_DENSE_FORWARDSUB_5(beliefMPC_Ld02, beliefMPC_bmy02, beliefMPC_yy02);
beliefMPC_LA_VSUB6_INDEXED_7_7_7(beliefMPC_ccrhsub04, beliefMPC_sub04, beliefMPC_ubIdx04, beliefMPC_ccrhsl04, beliefMPC_slb04, beliefMPC_lbIdx04, beliefMPC_rd04);
beliefMPC_LA_DIAG_FORWARDSUB_7(beliefMPC_Phi04, beliefMPC_rd04, beliefMPC_Lbyrd04);
beliefMPC_LA_DENSE_DIAGZERO_2MVMADD_5_7_7(beliefMPC_V03, beliefMPC_Lbyrd03, beliefMPC_W04, beliefMPC_Lbyrd04, beliefMPC_beta03);
beliefMPC_LA_DENSE_MVMSUB1_5_5(beliefMPC_Lsd03, beliefMPC_yy02, beliefMPC_beta03, beliefMPC_bmy03);
beliefMPC_LA_DENSE_FORWARDSUB_5(beliefMPC_Ld03, beliefMPC_bmy03, beliefMPC_yy03);
beliefMPC_LA_VSUB6_INDEXED_7_7_7(beliefMPC_ccrhsub05, beliefMPC_sub05, beliefMPC_ubIdx05, beliefMPC_ccrhsl05, beliefMPC_slb05, beliefMPC_lbIdx05, beliefMPC_rd05);
beliefMPC_LA_DIAG_FORWARDSUB_7(beliefMPC_Phi05, beliefMPC_rd05, beliefMPC_Lbyrd05);
beliefMPC_LA_DENSE_DIAGZERO_2MVMADD_5_7_7(beliefMPC_V04, beliefMPC_Lbyrd04, beliefMPC_W05, beliefMPC_Lbyrd05, beliefMPC_beta04);
beliefMPC_LA_DENSE_MVMSUB1_5_5(beliefMPC_Lsd04, beliefMPC_yy03, beliefMPC_beta04, beliefMPC_bmy04);
beliefMPC_LA_DENSE_FORWARDSUB_5(beliefMPC_Ld04, beliefMPC_bmy04, beliefMPC_yy04);
beliefMPC_LA_VSUB6_INDEXED_7_7_7(beliefMPC_ccrhsub06, beliefMPC_sub06, beliefMPC_ubIdx06, beliefMPC_ccrhsl06, beliefMPC_slb06, beliefMPC_lbIdx06, beliefMPC_rd06);
beliefMPC_LA_DIAG_FORWARDSUB_7(beliefMPC_Phi06, beliefMPC_rd06, beliefMPC_Lbyrd06);
beliefMPC_LA_DENSE_DIAGZERO_2MVMADD_5_7_7(beliefMPC_V05, beliefMPC_Lbyrd05, beliefMPC_W06, beliefMPC_Lbyrd06, beliefMPC_beta05);
beliefMPC_LA_DENSE_MVMSUB1_5_5(beliefMPC_Lsd05, beliefMPC_yy04, beliefMPC_beta05, beliefMPC_bmy05);
beliefMPC_LA_DENSE_FORWARDSUB_5(beliefMPC_Ld05, beliefMPC_bmy05, beliefMPC_yy05);
beliefMPC_LA_VSUB6_INDEXED_7_7_7(beliefMPC_ccrhsub07, beliefMPC_sub07, beliefMPC_ubIdx07, beliefMPC_ccrhsl07, beliefMPC_slb07, beliefMPC_lbIdx07, beliefMPC_rd07);
beliefMPC_LA_DIAG_FORWARDSUB_7(beliefMPC_Phi07, beliefMPC_rd07, beliefMPC_Lbyrd07);
beliefMPC_LA_DENSE_DIAGZERO_2MVMADD_5_7_7(beliefMPC_V06, beliefMPC_Lbyrd06, beliefMPC_W07, beliefMPC_Lbyrd07, beliefMPC_beta06);
beliefMPC_LA_DENSE_MVMSUB1_5_5(beliefMPC_Lsd06, beliefMPC_yy05, beliefMPC_beta06, beliefMPC_bmy06);
beliefMPC_LA_DENSE_FORWARDSUB_5(beliefMPC_Ld06, beliefMPC_bmy06, beliefMPC_yy06);
beliefMPC_LA_VSUB6_INDEXED_7_7_7(beliefMPC_ccrhsub08, beliefMPC_sub08, beliefMPC_ubIdx08, beliefMPC_ccrhsl08, beliefMPC_slb08, beliefMPC_lbIdx08, beliefMPC_rd08);
beliefMPC_LA_DIAG_FORWARDSUB_7(beliefMPC_Phi08, beliefMPC_rd08, beliefMPC_Lbyrd08);
beliefMPC_LA_DENSE_DIAGZERO_2MVMADD_5_7_7(beliefMPC_V07, beliefMPC_Lbyrd07, beliefMPC_W08, beliefMPC_Lbyrd08, beliefMPC_beta07);
beliefMPC_LA_DENSE_MVMSUB1_5_5(beliefMPC_Lsd07, beliefMPC_yy06, beliefMPC_beta07, beliefMPC_bmy07);
beliefMPC_LA_DENSE_FORWARDSUB_5(beliefMPC_Ld07, beliefMPC_bmy07, beliefMPC_yy07);
beliefMPC_LA_VSUB6_INDEXED_7_7_7(beliefMPC_ccrhsub09, beliefMPC_sub09, beliefMPC_ubIdx09, beliefMPC_ccrhsl09, beliefMPC_slb09, beliefMPC_lbIdx09, beliefMPC_rd09);
beliefMPC_LA_DIAG_FORWARDSUB_7(beliefMPC_Phi09, beliefMPC_rd09, beliefMPC_Lbyrd09);
beliefMPC_LA_DENSE_DIAGZERO_2MVMADD_5_7_7(beliefMPC_V08, beliefMPC_Lbyrd08, beliefMPC_W09, beliefMPC_Lbyrd09, beliefMPC_beta08);
beliefMPC_LA_DENSE_MVMSUB1_5_5(beliefMPC_Lsd08, beliefMPC_yy07, beliefMPC_beta08, beliefMPC_bmy08);
beliefMPC_LA_DENSE_FORWARDSUB_5(beliefMPC_Ld08, beliefMPC_bmy08, beliefMPC_yy08);
beliefMPC_LA_VSUB6_INDEXED_7_7_7(beliefMPC_ccrhsub10, beliefMPC_sub10, beliefMPC_ubIdx10, beliefMPC_ccrhsl10, beliefMPC_slb10, beliefMPC_lbIdx10, beliefMPC_rd10);
beliefMPC_LA_DIAG_FORWARDSUB_7(beliefMPC_Phi10, beliefMPC_rd10, beliefMPC_Lbyrd10);
beliefMPC_LA_DENSE_DIAGZERO_2MVMADD_5_7_7(beliefMPC_V09, beliefMPC_Lbyrd09, beliefMPC_W10, beliefMPC_Lbyrd10, beliefMPC_beta09);
beliefMPC_LA_DENSE_MVMSUB1_5_5(beliefMPC_Lsd09, beliefMPC_yy08, beliefMPC_beta09, beliefMPC_bmy09);
beliefMPC_LA_DENSE_FORWARDSUB_5(beliefMPC_Ld09, beliefMPC_bmy09, beliefMPC_yy09);
beliefMPC_LA_VSUB6_INDEXED_7_7_7(beliefMPC_ccrhsub11, beliefMPC_sub11, beliefMPC_ubIdx11, beliefMPC_ccrhsl11, beliefMPC_slb11, beliefMPC_lbIdx11, beliefMPC_rd11);
beliefMPC_LA_DIAG_FORWARDSUB_7(beliefMPC_Phi11, beliefMPC_rd11, beliefMPC_Lbyrd11);
beliefMPC_LA_DENSE_DIAGZERO_2MVMADD_5_7_7(beliefMPC_V10, beliefMPC_Lbyrd10, beliefMPC_W11, beliefMPC_Lbyrd11, beliefMPC_beta10);
beliefMPC_LA_DENSE_MVMSUB1_5_5(beliefMPC_Lsd10, beliefMPC_yy09, beliefMPC_beta10, beliefMPC_bmy10);
beliefMPC_LA_DENSE_FORWARDSUB_5(beliefMPC_Ld10, beliefMPC_bmy10, beliefMPC_yy10);
beliefMPC_LA_VSUB6_INDEXED_7_7_7(beliefMPC_ccrhsub12, beliefMPC_sub12, beliefMPC_ubIdx12, beliefMPC_ccrhsl12, beliefMPC_slb12, beliefMPC_lbIdx12, beliefMPC_rd12);
beliefMPC_LA_DIAG_FORWARDSUB_7(beliefMPC_Phi12, beliefMPC_rd12, beliefMPC_Lbyrd12);
beliefMPC_LA_DENSE_DIAGZERO_2MVMADD_5_7_7(beliefMPC_V11, beliefMPC_Lbyrd11, beliefMPC_W12, beliefMPC_Lbyrd12, beliefMPC_beta11);
beliefMPC_LA_DENSE_MVMSUB1_5_5(beliefMPC_Lsd11, beliefMPC_yy10, beliefMPC_beta11, beliefMPC_bmy11);
beliefMPC_LA_DENSE_FORWARDSUB_5(beliefMPC_Ld11, beliefMPC_bmy11, beliefMPC_yy11);
beliefMPC_LA_VSUB6_INDEXED_7_7_7(beliefMPC_ccrhsub13, beliefMPC_sub13, beliefMPC_ubIdx13, beliefMPC_ccrhsl13, beliefMPC_slb13, beliefMPC_lbIdx13, beliefMPC_rd13);
beliefMPC_LA_DIAG_FORWARDSUB_7(beliefMPC_Phi13, beliefMPC_rd13, beliefMPC_Lbyrd13);
beliefMPC_LA_DENSE_DIAGZERO_2MVMADD_5_7_7(beliefMPC_V12, beliefMPC_Lbyrd12, beliefMPC_W13, beliefMPC_Lbyrd13, beliefMPC_beta12);
beliefMPC_LA_DENSE_MVMSUB1_5_5(beliefMPC_Lsd12, beliefMPC_yy11, beliefMPC_beta12, beliefMPC_bmy12);
beliefMPC_LA_DENSE_FORWARDSUB_5(beliefMPC_Ld12, beliefMPC_bmy12, beliefMPC_yy12);
beliefMPC_LA_VSUB6_INDEXED_5_5_5(beliefMPC_ccrhsub14, beliefMPC_sub14, beliefMPC_ubIdx14, beliefMPC_ccrhsl14, beliefMPC_slb14, beliefMPC_lbIdx14, beliefMPC_rd14);
beliefMPC_LA_DIAG_FORWARDSUB_5(beliefMPC_Phi14, beliefMPC_rd14, beliefMPC_Lbyrd14);
beliefMPC_LA_DENSE_DIAGZERO_2MVMADD_5_7_5(beliefMPC_V13, beliefMPC_Lbyrd13, beliefMPC_W14, beliefMPC_Lbyrd14, beliefMPC_beta13);
beliefMPC_LA_DENSE_MVMSUB1_5_5(beliefMPC_Lsd13, beliefMPC_yy12, beliefMPC_beta13, beliefMPC_bmy13);
beliefMPC_LA_DENSE_FORWARDSUB_5(beliefMPC_Ld13, beliefMPC_bmy13, beliefMPC_yy13);
beliefMPC_LA_DENSE_BACKWARDSUB_5(beliefMPC_Ld13, beliefMPC_yy13, beliefMPC_dvcc13);
beliefMPC_LA_DENSE_MTVMSUB_5_5(beliefMPC_Lsd13, beliefMPC_dvcc13, beliefMPC_yy12, beliefMPC_bmy12);
beliefMPC_LA_DENSE_BACKWARDSUB_5(beliefMPC_Ld12, beliefMPC_bmy12, beliefMPC_dvcc12);
beliefMPC_LA_DENSE_MTVMSUB_5_5(beliefMPC_Lsd12, beliefMPC_dvcc12, beliefMPC_yy11, beliefMPC_bmy11);
beliefMPC_LA_DENSE_BACKWARDSUB_5(beliefMPC_Ld11, beliefMPC_bmy11, beliefMPC_dvcc11);
beliefMPC_LA_DENSE_MTVMSUB_5_5(beliefMPC_Lsd11, beliefMPC_dvcc11, beliefMPC_yy10, beliefMPC_bmy10);
beliefMPC_LA_DENSE_BACKWARDSUB_5(beliefMPC_Ld10, beliefMPC_bmy10, beliefMPC_dvcc10);
beliefMPC_LA_DENSE_MTVMSUB_5_5(beliefMPC_Lsd10, beliefMPC_dvcc10, beliefMPC_yy09, beliefMPC_bmy09);
beliefMPC_LA_DENSE_BACKWARDSUB_5(beliefMPC_Ld09, beliefMPC_bmy09, beliefMPC_dvcc09);
beliefMPC_LA_DENSE_MTVMSUB_5_5(beliefMPC_Lsd09, beliefMPC_dvcc09, beliefMPC_yy08, beliefMPC_bmy08);
beliefMPC_LA_DENSE_BACKWARDSUB_5(beliefMPC_Ld08, beliefMPC_bmy08, beliefMPC_dvcc08);
beliefMPC_LA_DENSE_MTVMSUB_5_5(beliefMPC_Lsd08, beliefMPC_dvcc08, beliefMPC_yy07, beliefMPC_bmy07);
beliefMPC_LA_DENSE_BACKWARDSUB_5(beliefMPC_Ld07, beliefMPC_bmy07, beliefMPC_dvcc07);
beliefMPC_LA_DENSE_MTVMSUB_5_5(beliefMPC_Lsd07, beliefMPC_dvcc07, beliefMPC_yy06, beliefMPC_bmy06);
beliefMPC_LA_DENSE_BACKWARDSUB_5(beliefMPC_Ld06, beliefMPC_bmy06, beliefMPC_dvcc06);
beliefMPC_LA_DENSE_MTVMSUB_5_5(beliefMPC_Lsd06, beliefMPC_dvcc06, beliefMPC_yy05, beliefMPC_bmy05);
beliefMPC_LA_DENSE_BACKWARDSUB_5(beliefMPC_Ld05, beliefMPC_bmy05, beliefMPC_dvcc05);
beliefMPC_LA_DENSE_MTVMSUB_5_5(beliefMPC_Lsd05, beliefMPC_dvcc05, beliefMPC_yy04, beliefMPC_bmy04);
beliefMPC_LA_DENSE_BACKWARDSUB_5(beliefMPC_Ld04, beliefMPC_bmy04, beliefMPC_dvcc04);
beliefMPC_LA_DENSE_MTVMSUB_5_5(beliefMPC_Lsd04, beliefMPC_dvcc04, beliefMPC_yy03, beliefMPC_bmy03);
beliefMPC_LA_DENSE_BACKWARDSUB_5(beliefMPC_Ld03, beliefMPC_bmy03, beliefMPC_dvcc03);
beliefMPC_LA_DENSE_MTVMSUB_5_5(beliefMPC_Lsd03, beliefMPC_dvcc03, beliefMPC_yy02, beliefMPC_bmy02);
beliefMPC_LA_DENSE_BACKWARDSUB_5(beliefMPC_Ld02, beliefMPC_bmy02, beliefMPC_dvcc02);
beliefMPC_LA_DENSE_MTVMSUB_5_5(beliefMPC_Lsd02, beliefMPC_dvcc02, beliefMPC_yy01, beliefMPC_bmy01);
beliefMPC_LA_DENSE_BACKWARDSUB_5(beliefMPC_Ld01, beliefMPC_bmy01, beliefMPC_dvcc01);
beliefMPC_LA_DENSE_MTVMSUB_5_10(beliefMPC_Lsd01, beliefMPC_dvcc01, beliefMPC_yy00, beliefMPC_bmy00);
beliefMPC_LA_DENSE_BACKWARDSUB_10(beliefMPC_Ld00, beliefMPC_bmy00, beliefMPC_dvcc00);
beliefMPC_LA_DENSE_MTVM_10_7(params->C1, beliefMPC_dvcc00, beliefMPC_grad_eq00);
beliefMPC_LA_DENSE_MTVM2_5_7_10(params->C2, beliefMPC_dvcc01, beliefMPC_D01, beliefMPC_dvcc00, beliefMPC_grad_eq01);
beliefMPC_LA_DENSE_DIAGZERO_MTVM2_5_7_5(params->C3, beliefMPC_dvcc02, beliefMPC_D02, beliefMPC_dvcc01, beliefMPC_grad_eq02);
beliefMPC_LA_DENSE_DIAGZERO_MTVM2_5_7_5(params->C4, beliefMPC_dvcc03, beliefMPC_D02, beliefMPC_dvcc02, beliefMPC_grad_eq03);
beliefMPC_LA_DENSE_DIAGZERO_MTVM2_5_7_5(params->C5, beliefMPC_dvcc04, beliefMPC_D02, beliefMPC_dvcc03, beliefMPC_grad_eq04);
beliefMPC_LA_DENSE_DIAGZERO_MTVM2_5_7_5(params->C6, beliefMPC_dvcc05, beliefMPC_D02, beliefMPC_dvcc04, beliefMPC_grad_eq05);
beliefMPC_LA_DENSE_DIAGZERO_MTVM2_5_7_5(params->C7, beliefMPC_dvcc06, beliefMPC_D02, beliefMPC_dvcc05, beliefMPC_grad_eq06);
beliefMPC_LA_DENSE_DIAGZERO_MTVM2_5_7_5(params->C8, beliefMPC_dvcc07, beliefMPC_D02, beliefMPC_dvcc06, beliefMPC_grad_eq07);
beliefMPC_LA_DENSE_DIAGZERO_MTVM2_5_7_5(params->C9, beliefMPC_dvcc08, beliefMPC_D02, beliefMPC_dvcc07, beliefMPC_grad_eq08);
beliefMPC_LA_DENSE_DIAGZERO_MTVM2_5_7_5(params->C10, beliefMPC_dvcc09, beliefMPC_D02, beliefMPC_dvcc08, beliefMPC_grad_eq09);
beliefMPC_LA_DENSE_DIAGZERO_MTVM2_5_7_5(params->C11, beliefMPC_dvcc10, beliefMPC_D02, beliefMPC_dvcc09, beliefMPC_grad_eq10);
beliefMPC_LA_DENSE_DIAGZERO_MTVM2_5_7_5(params->C12, beliefMPC_dvcc11, beliefMPC_D02, beliefMPC_dvcc10, beliefMPC_grad_eq11);
beliefMPC_LA_DENSE_DIAGZERO_MTVM2_5_7_5(params->C13, beliefMPC_dvcc12, beliefMPC_D02, beliefMPC_dvcc11, beliefMPC_grad_eq12);
beliefMPC_LA_DENSE_DIAGZERO_MTVM2_5_7_5(params->C14, beliefMPC_dvcc13, beliefMPC_D02, beliefMPC_dvcc12, beliefMPC_grad_eq13);
beliefMPC_LA_DIAGZERO_MTVM_5_5(beliefMPC_D14, beliefMPC_dvcc13, beliefMPC_grad_eq14);
beliefMPC_LA_VSUB_103(beliefMPC_rd, beliefMPC_grad_eq, beliefMPC_rd);
beliefMPC_LA_DIAG_FORWARDBACKWARDSUB_7(beliefMPC_Phi00, beliefMPC_rd00, beliefMPC_dzcc00);
beliefMPC_LA_DIAG_FORWARDBACKWARDSUB_7(beliefMPC_Phi01, beliefMPC_rd01, beliefMPC_dzcc01);
beliefMPC_LA_DIAG_FORWARDBACKWARDSUB_7(beliefMPC_Phi02, beliefMPC_rd02, beliefMPC_dzcc02);
beliefMPC_LA_DIAG_FORWARDBACKWARDSUB_7(beliefMPC_Phi03, beliefMPC_rd03, beliefMPC_dzcc03);
beliefMPC_LA_DIAG_FORWARDBACKWARDSUB_7(beliefMPC_Phi04, beliefMPC_rd04, beliefMPC_dzcc04);
beliefMPC_LA_DIAG_FORWARDBACKWARDSUB_7(beliefMPC_Phi05, beliefMPC_rd05, beliefMPC_dzcc05);
beliefMPC_LA_DIAG_FORWARDBACKWARDSUB_7(beliefMPC_Phi06, beliefMPC_rd06, beliefMPC_dzcc06);
beliefMPC_LA_DIAG_FORWARDBACKWARDSUB_7(beliefMPC_Phi07, beliefMPC_rd07, beliefMPC_dzcc07);
beliefMPC_LA_DIAG_FORWARDBACKWARDSUB_7(beliefMPC_Phi08, beliefMPC_rd08, beliefMPC_dzcc08);
beliefMPC_LA_DIAG_FORWARDBACKWARDSUB_7(beliefMPC_Phi09, beliefMPC_rd09, beliefMPC_dzcc09);
beliefMPC_LA_DIAG_FORWARDBACKWARDSUB_7(beliefMPC_Phi10, beliefMPC_rd10, beliefMPC_dzcc10);
beliefMPC_LA_DIAG_FORWARDBACKWARDSUB_7(beliefMPC_Phi11, beliefMPC_rd11, beliefMPC_dzcc11);
beliefMPC_LA_DIAG_FORWARDBACKWARDSUB_7(beliefMPC_Phi12, beliefMPC_rd12, beliefMPC_dzcc12);
beliefMPC_LA_DIAG_FORWARDBACKWARDSUB_7(beliefMPC_Phi13, beliefMPC_rd13, beliefMPC_dzcc13);
beliefMPC_LA_DIAG_FORWARDBACKWARDSUB_5(beliefMPC_Phi14, beliefMPC_rd14, beliefMPC_dzcc14);
beliefMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_7(beliefMPC_ccrhsl00, beliefMPC_slb00, beliefMPC_llbbyslb00, beliefMPC_dzcc00, beliefMPC_lbIdx00, beliefMPC_dllbcc00);
beliefMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_7(beliefMPC_ccrhsub00, beliefMPC_sub00, beliefMPC_lubbysub00, beliefMPC_dzcc00, beliefMPC_ubIdx00, beliefMPC_dlubcc00);
beliefMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_7(beliefMPC_ccrhsl01, beliefMPC_slb01, beliefMPC_llbbyslb01, beliefMPC_dzcc01, beliefMPC_lbIdx01, beliefMPC_dllbcc01);
beliefMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_7(beliefMPC_ccrhsub01, beliefMPC_sub01, beliefMPC_lubbysub01, beliefMPC_dzcc01, beliefMPC_ubIdx01, beliefMPC_dlubcc01);
beliefMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_7(beliefMPC_ccrhsl02, beliefMPC_slb02, beliefMPC_llbbyslb02, beliefMPC_dzcc02, beliefMPC_lbIdx02, beliefMPC_dllbcc02);
beliefMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_7(beliefMPC_ccrhsub02, beliefMPC_sub02, beliefMPC_lubbysub02, beliefMPC_dzcc02, beliefMPC_ubIdx02, beliefMPC_dlubcc02);
beliefMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_7(beliefMPC_ccrhsl03, beliefMPC_slb03, beliefMPC_llbbyslb03, beliefMPC_dzcc03, beliefMPC_lbIdx03, beliefMPC_dllbcc03);
beliefMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_7(beliefMPC_ccrhsub03, beliefMPC_sub03, beliefMPC_lubbysub03, beliefMPC_dzcc03, beliefMPC_ubIdx03, beliefMPC_dlubcc03);
beliefMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_7(beliefMPC_ccrhsl04, beliefMPC_slb04, beliefMPC_llbbyslb04, beliefMPC_dzcc04, beliefMPC_lbIdx04, beliefMPC_dllbcc04);
beliefMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_7(beliefMPC_ccrhsub04, beliefMPC_sub04, beliefMPC_lubbysub04, beliefMPC_dzcc04, beliefMPC_ubIdx04, beliefMPC_dlubcc04);
beliefMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_7(beliefMPC_ccrhsl05, beliefMPC_slb05, beliefMPC_llbbyslb05, beliefMPC_dzcc05, beliefMPC_lbIdx05, beliefMPC_dllbcc05);
beliefMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_7(beliefMPC_ccrhsub05, beliefMPC_sub05, beliefMPC_lubbysub05, beliefMPC_dzcc05, beliefMPC_ubIdx05, beliefMPC_dlubcc05);
beliefMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_7(beliefMPC_ccrhsl06, beliefMPC_slb06, beliefMPC_llbbyslb06, beliefMPC_dzcc06, beliefMPC_lbIdx06, beliefMPC_dllbcc06);
beliefMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_7(beliefMPC_ccrhsub06, beliefMPC_sub06, beliefMPC_lubbysub06, beliefMPC_dzcc06, beliefMPC_ubIdx06, beliefMPC_dlubcc06);
beliefMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_7(beliefMPC_ccrhsl07, beliefMPC_slb07, beliefMPC_llbbyslb07, beliefMPC_dzcc07, beliefMPC_lbIdx07, beliefMPC_dllbcc07);
beliefMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_7(beliefMPC_ccrhsub07, beliefMPC_sub07, beliefMPC_lubbysub07, beliefMPC_dzcc07, beliefMPC_ubIdx07, beliefMPC_dlubcc07);
beliefMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_7(beliefMPC_ccrhsl08, beliefMPC_slb08, beliefMPC_llbbyslb08, beliefMPC_dzcc08, beliefMPC_lbIdx08, beliefMPC_dllbcc08);
beliefMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_7(beliefMPC_ccrhsub08, beliefMPC_sub08, beliefMPC_lubbysub08, beliefMPC_dzcc08, beliefMPC_ubIdx08, beliefMPC_dlubcc08);
beliefMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_7(beliefMPC_ccrhsl09, beliefMPC_slb09, beliefMPC_llbbyslb09, beliefMPC_dzcc09, beliefMPC_lbIdx09, beliefMPC_dllbcc09);
beliefMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_7(beliefMPC_ccrhsub09, beliefMPC_sub09, beliefMPC_lubbysub09, beliefMPC_dzcc09, beliefMPC_ubIdx09, beliefMPC_dlubcc09);
beliefMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_7(beliefMPC_ccrhsl10, beliefMPC_slb10, beliefMPC_llbbyslb10, beliefMPC_dzcc10, beliefMPC_lbIdx10, beliefMPC_dllbcc10);
beliefMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_7(beliefMPC_ccrhsub10, beliefMPC_sub10, beliefMPC_lubbysub10, beliefMPC_dzcc10, beliefMPC_ubIdx10, beliefMPC_dlubcc10);
beliefMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_7(beliefMPC_ccrhsl11, beliefMPC_slb11, beliefMPC_llbbyslb11, beliefMPC_dzcc11, beliefMPC_lbIdx11, beliefMPC_dllbcc11);
beliefMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_7(beliefMPC_ccrhsub11, beliefMPC_sub11, beliefMPC_lubbysub11, beliefMPC_dzcc11, beliefMPC_ubIdx11, beliefMPC_dlubcc11);
beliefMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_7(beliefMPC_ccrhsl12, beliefMPC_slb12, beliefMPC_llbbyslb12, beliefMPC_dzcc12, beliefMPC_lbIdx12, beliefMPC_dllbcc12);
beliefMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_7(beliefMPC_ccrhsub12, beliefMPC_sub12, beliefMPC_lubbysub12, beliefMPC_dzcc12, beliefMPC_ubIdx12, beliefMPC_dlubcc12);
beliefMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_7(beliefMPC_ccrhsl13, beliefMPC_slb13, beliefMPC_llbbyslb13, beliefMPC_dzcc13, beliefMPC_lbIdx13, beliefMPC_dllbcc13);
beliefMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_7(beliefMPC_ccrhsub13, beliefMPC_sub13, beliefMPC_lubbysub13, beliefMPC_dzcc13, beliefMPC_ubIdx13, beliefMPC_dlubcc13);
beliefMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_5(beliefMPC_ccrhsl14, beliefMPC_slb14, beliefMPC_llbbyslb14, beliefMPC_dzcc14, beliefMPC_lbIdx14, beliefMPC_dllbcc14);
beliefMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(beliefMPC_ccrhsub14, beliefMPC_sub14, beliefMPC_lubbysub14, beliefMPC_dzcc14, beliefMPC_ubIdx14, beliefMPC_dlubcc14);
beliefMPC_LA_VSUB7_206(beliefMPC_l, beliefMPC_ccrhs, beliefMPC_s, beliefMPC_dl_cc, beliefMPC_ds_cc);
beliefMPC_LA_VADD_103(beliefMPC_dz_cc, beliefMPC_dz_aff);
beliefMPC_LA_VADD_75(beliefMPC_dv_cc, beliefMPC_dv_aff);
beliefMPC_LA_VADD_206(beliefMPC_dl_cc, beliefMPC_dl_aff);
beliefMPC_LA_VADD_206(beliefMPC_ds_cc, beliefMPC_ds_aff);
info->lsit_cc = beliefMPC_LINESEARCH_BACKTRACKING_COMBINED(beliefMPC_z, beliefMPC_v, beliefMPC_l, beliefMPC_s, beliefMPC_dz_cc, beliefMPC_dv_cc, beliefMPC_dl_cc, beliefMPC_ds_cc, &info->step_cc, &info->mu);
if( info->lsit_cc == beliefMPC_NOPROGRESS ){
exitcode = beliefMPC_NOPROGRESS; break;
}
info->it++;
}
output->z1[0] = beliefMPC_z00[0];
output->z1[1] = beliefMPC_z00[1];
output->z1[2] = beliefMPC_z00[2];
output->z1[3] = beliefMPC_z00[3];
output->z1[4] = beliefMPC_z00[4];
output->z1[5] = beliefMPC_z00[5];
output->z1[6] = beliefMPC_z00[6];
output->z2[0] = beliefMPC_z01[0];
output->z2[1] = beliefMPC_z01[1];
output->z2[2] = beliefMPC_z01[2];
output->z2[3] = beliefMPC_z01[3];
output->z2[4] = beliefMPC_z01[4];
output->z2[5] = beliefMPC_z01[5];
output->z2[6] = beliefMPC_z01[6];
output->z3[0] = beliefMPC_z02[0];
output->z3[1] = beliefMPC_z02[1];
output->z3[2] = beliefMPC_z02[2];
output->z3[3] = beliefMPC_z02[3];
output->z3[4] = beliefMPC_z02[4];
output->z3[5] = beliefMPC_z02[5];
output->z3[6] = beliefMPC_z02[6];
output->z4[0] = beliefMPC_z03[0];
output->z4[1] = beliefMPC_z03[1];
output->z4[2] = beliefMPC_z03[2];
output->z4[3] = beliefMPC_z03[3];
output->z4[4] = beliefMPC_z03[4];
output->z4[5] = beliefMPC_z03[5];
output->z4[6] = beliefMPC_z03[6];
output->z5[0] = beliefMPC_z04[0];
output->z5[1] = beliefMPC_z04[1];
output->z5[2] = beliefMPC_z04[2];
output->z5[3] = beliefMPC_z04[3];
output->z5[4] = beliefMPC_z04[4];
output->z5[5] = beliefMPC_z04[5];
output->z5[6] = beliefMPC_z04[6];
output->z6[0] = beliefMPC_z05[0];
output->z6[1] = beliefMPC_z05[1];
output->z6[2] = beliefMPC_z05[2];
output->z6[3] = beliefMPC_z05[3];
output->z6[4] = beliefMPC_z05[4];
output->z6[5] = beliefMPC_z05[5];
output->z6[6] = beliefMPC_z05[6];
output->z7[0] = beliefMPC_z06[0];
output->z7[1] = beliefMPC_z06[1];
output->z7[2] = beliefMPC_z06[2];
output->z7[3] = beliefMPC_z06[3];
output->z7[4] = beliefMPC_z06[4];
output->z7[5] = beliefMPC_z06[5];
output->z7[6] = beliefMPC_z06[6];
output->z8[0] = beliefMPC_z07[0];
output->z8[1] = beliefMPC_z07[1];
output->z8[2] = beliefMPC_z07[2];
output->z8[3] = beliefMPC_z07[3];
output->z8[4] = beliefMPC_z07[4];
output->z8[5] = beliefMPC_z07[5];
output->z8[6] = beliefMPC_z07[6];
output->z9[0] = beliefMPC_z08[0];
output->z9[1] = beliefMPC_z08[1];
output->z9[2] = beliefMPC_z08[2];
output->z9[3] = beliefMPC_z08[3];
output->z9[4] = beliefMPC_z08[4];
output->z9[5] = beliefMPC_z08[5];
output->z9[6] = beliefMPC_z08[6];
output->z10[0] = beliefMPC_z09[0];
output->z10[1] = beliefMPC_z09[1];
output->z10[2] = beliefMPC_z09[2];
output->z10[3] = beliefMPC_z09[3];
output->z10[4] = beliefMPC_z09[4];
output->z10[5] = beliefMPC_z09[5];
output->z10[6] = beliefMPC_z09[6];
output->z11[0] = beliefMPC_z10[0];
output->z11[1] = beliefMPC_z10[1];
output->z11[2] = beliefMPC_z10[2];
output->z11[3] = beliefMPC_z10[3];
output->z11[4] = beliefMPC_z10[4];
output->z11[5] = beliefMPC_z10[5];
output->z11[6] = beliefMPC_z10[6];
output->z12[0] = beliefMPC_z11[0];
output->z12[1] = beliefMPC_z11[1];
output->z12[2] = beliefMPC_z11[2];
output->z12[3] = beliefMPC_z11[3];
output->z12[4] = beliefMPC_z11[4];
output->z12[5] = beliefMPC_z11[5];
output->z12[6] = beliefMPC_z11[6];
output->z13[0] = beliefMPC_z12[0];
output->z13[1] = beliefMPC_z12[1];
output->z13[2] = beliefMPC_z12[2];
output->z13[3] = beliefMPC_z12[3];
output->z13[4] = beliefMPC_z12[4];
output->z13[5] = beliefMPC_z12[5];
output->z13[6] = beliefMPC_z12[6];
output->z14[0] = beliefMPC_z13[0];
output->z14[1] = beliefMPC_z13[1];
output->z14[2] = beliefMPC_z13[2];
output->z14[3] = beliefMPC_z13[3];
output->z14[4] = beliefMPC_z13[4];
output->z14[5] = beliefMPC_z13[5];
output->z14[6] = beliefMPC_z13[6];
output->z15[0] = beliefMPC_z14[0];
output->z15[1] = beliefMPC_z14[1];
output->z15[2] = beliefMPC_z14[2];
output->z15[3] = beliefMPC_z14[3];
output->z15[4] = beliefMPC_z14[4];

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
