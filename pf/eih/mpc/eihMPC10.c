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

#include "eihMPC.h"

/* for square root */
#include <math.h> 

/* SAFE DIVISION ------------------------------------------------------- */
#define MAX(X,Y)  ((X) < (Y) ? (Y) : (X))
#define MIN(X,Y)  ((X) < (Y) ? (X) : (Y))
/*#define SAFEDIV_POS(X,Y)  ( (Y) < EPS ? ((X)/EPS) : (X)/(Y) ) 
#define EPS (1.0000E-013) */
#define BIGM (1E24)
#define BIGMM (1E36)

/* includes for parallel computation if necessary */


/* SYSTEM INCLUDES FOR PRINTING ---------------------------------------- */




/* LINEAR ALGEBRA LIBRARY ---------------------------------------------- */
/*
 * Initializes a vector of length 133 with a value.
 */
void eihMPC_LA_INITIALIZEVECTOR_133(eihMPC_FLOAT* vec, eihMPC_FLOAT value)
{
	int i;
	for( i=0; i<133; i++ )
	{
		vec[i] = value;
	}
}


/*
 * Initializes a vector of length 70 with a value.
 */
void eihMPC_LA_INITIALIZEVECTOR_70(eihMPC_FLOAT* vec, eihMPC_FLOAT value)
{
	int i;
	for( i=0; i<70; i++ )
	{
		vec[i] = value;
	}
}


/*
 * Initializes a vector of length 266 with a value.
 */
void eihMPC_LA_INITIALIZEVECTOR_266(eihMPC_FLOAT* vec, eihMPC_FLOAT value)
{
	int i;
	for( i=0; i<266; i++ )
	{
		vec[i] = value;
	}
}


/* 
 * Calculates a dot product and adds it to a variable: z += x'*y; 
 * This function is for vectors of length 266.
 */
void eihMPC_LA_DOTACC_266(eihMPC_FLOAT *x, eihMPC_FLOAT *y, eihMPC_FLOAT *z)
{
	int i;
	for( i=0; i<266; i++ ){
		*z += x[i]*y[i];
	}
}


/*
 * Calculates the gradient and the value for a quadratic function 0.5*z'*H*z + f'*z
 *
 * INPUTS:     H  - Symmetric Hessian, diag matrix of size [14 x 14]
 *             f  - column vector of size 14
 *             z  - column vector of size 14
 *
 * OUTPUTS: grad  - gradient at z (= H*z + f), column vector of size 14
 *          value <-- value + 0.5*z'*H*z + f'*z (value will be modified)
 */
void eihMPC_LA_DIAG_QUADFCN_14(eihMPC_FLOAT* H, eihMPC_FLOAT* f, eihMPC_FLOAT* z, eihMPC_FLOAT* grad, eihMPC_FLOAT* value)
{
	int i;
	eihMPC_FLOAT hz;	
	for( i=0; i<14; i++){
		hz = H[i]*z[i];
		grad[i] = hz + f[i];
		*value += 0.5*hz*z[i] + f[i]*z[i];
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
void eihMPC_LA_DIAG_QUADFCN_7(eihMPC_FLOAT* H, eihMPC_FLOAT* f, eihMPC_FLOAT* z, eihMPC_FLOAT* grad, eihMPC_FLOAT* value)
{
	int i;
	eihMPC_FLOAT hz;	
	for( i=0; i<7; i++){
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
void eihMPC_LA_DIAGZERO_MVMSUB6_7(eihMPC_FLOAT *B, eihMPC_FLOAT *u, eihMPC_FLOAT *b, eihMPC_FLOAT *l, eihMPC_FLOAT *r, eihMPC_FLOAT *z, eihMPC_FLOAT *y)
{
	int i;
	eihMPC_FLOAT Bu[7];
	eihMPC_FLOAT norm = *y;
	eihMPC_FLOAT lr = 0;

	/* do A*x + B*u first */
	for( i=0; i<7; i++ ){
		Bu[i] = B[i]*u[i];
	}	

	for( i=0; i<7; i++ ){
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
void eihMPC_LA_DENSE_DIAGZERO_MVMSUB3_7_14_14(eihMPC_FLOAT *A, eihMPC_FLOAT *x, eihMPC_FLOAT *B, eihMPC_FLOAT *u, eihMPC_FLOAT *b, eihMPC_FLOAT *l, eihMPC_FLOAT *r, eihMPC_FLOAT *z, eihMPC_FLOAT *y)
{
	int i;
	int j;
	int k = 0;
	eihMPC_FLOAT AxBu[7];
	eihMPC_FLOAT norm = *y;
	eihMPC_FLOAT lr = 0;

	/* do A*x + B*u first */
	for( i=0; i<7; i++ ){
		AxBu[i] = A[k++]*x[0] + B[i]*u[i];
	}	

	for( j=1; j<14; j++ ){		
		for( i=0; i<7; i++ ){
			AxBu[i] += A[k++]*x[j];
		}
	}

	for( i=0; i<7; i++ ){
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
void eihMPC_LA_DENSE_DIAGZERO_MVMSUB3_7_14_7(eihMPC_FLOAT *A, eihMPC_FLOAT *x, eihMPC_FLOAT *B, eihMPC_FLOAT *u, eihMPC_FLOAT *b, eihMPC_FLOAT *l, eihMPC_FLOAT *r, eihMPC_FLOAT *z, eihMPC_FLOAT *y)
{
	int i;
	int j;
	int k = 0;
	eihMPC_FLOAT AxBu[7];
	eihMPC_FLOAT norm = *y;
	eihMPC_FLOAT lr = 0;

	/* do A*x + B*u first */
	for( i=0; i<7; i++ ){
		AxBu[i] = A[k++]*x[0] + B[i]*u[i];
	}	

	for( j=1; j<14; j++ ){		
		for( i=0; i<7; i++ ){
			AxBu[i] += A[k++]*x[j];
		}
	}

	for( i=0; i<7; i++ ){
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
 * where A is of size [7 x 14] and stored in column major format.
 * and B is of size [7 x 14] and stored in diagzero format
 * Note the transposes of A and B!
 */
void eihMPC_LA_DENSE_DIAGZERO_MTVM2_7_14_7(eihMPC_FLOAT *A, eihMPC_FLOAT *x, eihMPC_FLOAT *B, eihMPC_FLOAT *y, eihMPC_FLOAT *z)
{
	int i;
	int j;
	int k = 0;
	for( i=0; i<7; i++ ){
		z[i] = 0;
		for( j=0; j<7; j++ ){
			z[i] += A[k++]*x[j];
		}
		z[i] += B[i]*y[i];
	}
	for( i=7 ;i<14; i++ ){
		z[i] = 0;
		for( j=0; j<7; j++ ){
			z[i] += A[k++]*x[j];
		}
	}
}


/*
 * Matrix vector multiplication y = M'*x where M is of size [7 x 7]
 * and stored in diagzero format. Note the transpose of M!
 */
void eihMPC_LA_DIAGZERO_MTVM_7_7(eihMPC_FLOAT *M, eihMPC_FLOAT *x, eihMPC_FLOAT *y)
{
	int i;
	for( i=0; i<7; i++ ){
		y[i] = M[i]*x[i];
	}
}


/*
 * Vector subtraction and addition.
 *	 Input: five vectors t, tidx, u, v, w and two scalars z and r
 *	 Output: y = t(tidx) - u + w
 *           z = z - v'*x;
 *           r = max([norm(y,inf), z]);
 * for vectors of length 14. Output z is of course scalar.
 */
void eihMPC_LA_VSUBADD3_14(eihMPC_FLOAT* t, eihMPC_FLOAT* u, int* uidx, eihMPC_FLOAT* v, eihMPC_FLOAT* w, eihMPC_FLOAT* y, eihMPC_FLOAT* z, eihMPC_FLOAT* r)
{
	int i;
	eihMPC_FLOAT norm = *r;
	eihMPC_FLOAT vx = 0;
	eihMPC_FLOAT x;
	for( i=0; i<14; i++){
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
 * for vectors of length 14. Output z is of course scalar.
 */
void eihMPC_LA_VSUBADD2_14(eihMPC_FLOAT* t, int* tidx, eihMPC_FLOAT* u, eihMPC_FLOAT* v, eihMPC_FLOAT* w, eihMPC_FLOAT* y, eihMPC_FLOAT* z, eihMPC_FLOAT* r)
{
	int i;
	eihMPC_FLOAT norm = *r;
	eihMPC_FLOAT vx = 0;
	eihMPC_FLOAT x;
	for( i=0; i<14; i++){
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
 * for vectors of length 7. Output z is of course scalar.
 */
void eihMPC_LA_VSUBADD3_7(eihMPC_FLOAT* t, eihMPC_FLOAT* u, int* uidx, eihMPC_FLOAT* v, eihMPC_FLOAT* w, eihMPC_FLOAT* y, eihMPC_FLOAT* z, eihMPC_FLOAT* r)
{
	int i;
	eihMPC_FLOAT norm = *r;
	eihMPC_FLOAT vx = 0;
	eihMPC_FLOAT x;
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
void eihMPC_LA_VSUBADD2_7(eihMPC_FLOAT* t, int* tidx, eihMPC_FLOAT* u, eihMPC_FLOAT* v, eihMPC_FLOAT* w, eihMPC_FLOAT* y, eihMPC_FLOAT* z, eihMPC_FLOAT* r)
{
	int i;
	eihMPC_FLOAT norm = *r;
	eihMPC_FLOAT vx = 0;
	eihMPC_FLOAT x;
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
 * Computes inequality constraints gradient-
 * Special function for box constraints of length 14
 * Returns also L/S, a value that is often used elsewhere.
 */
void eihMPC_LA_INEQ_B_GRAD_14_14_14(eihMPC_FLOAT *lu, eihMPC_FLOAT *su, eihMPC_FLOAT *ru, eihMPC_FLOAT *ll, eihMPC_FLOAT *sl, eihMPC_FLOAT *rl, int* lbIdx, int* ubIdx, eihMPC_FLOAT *grad, eihMPC_FLOAT *lubysu, eihMPC_FLOAT *llbysl)
{
	int i;
	for( i=0; i<14; i++ ){
		grad[i] = 0;
	}
	for( i=0; i<14; i++ ){		
		llbysl[i] = ll[i] / sl[i];
		grad[lbIdx[i]] -= llbysl[i]*rl[i];
	}
	for( i=0; i<14; i++ ){
		lubysu[i] = lu[i] / su[i];
		grad[ubIdx[i]] += lubysu[i]*ru[i];
	}
}


/*
 * Computes inequality constraints gradient-
 * Special function for box constraints of length 7
 * Returns also L/S, a value that is often used elsewhere.
 */
void eihMPC_LA_INEQ_B_GRAD_7_7_7(eihMPC_FLOAT *lu, eihMPC_FLOAT *su, eihMPC_FLOAT *ru, eihMPC_FLOAT *ll, eihMPC_FLOAT *sl, eihMPC_FLOAT *rl, int* lbIdx, int* ubIdx, eihMPC_FLOAT *grad, eihMPC_FLOAT *lubysu, eihMPC_FLOAT *llbysl)
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
 * Addition of three vectors  z = u + w + v
 * of length 133.
 */
void eihMPC_LA_VVADD3_133(eihMPC_FLOAT *u, eihMPC_FLOAT *v, eihMPC_FLOAT *w, eihMPC_FLOAT *z)
{
	int i;
	for( i=0; i<133; i++ ){
		z[i] = u[i] + v[i] + w[i];
	}
}


/*
 * Special function to compute the diagonal cholesky factorization of the 
 * positive definite augmented Hessian for block size 14.
 *
 * Inputs: - H = diagonal cost Hessian in diagonal storage format
 *         - llbysl = L / S of lower bounds
 *         - lubysu = L / S of upper bounds
 *
 * Output: Phi = sqrt(H + diag(llbysl) + diag(lubysu))
 * where Phi is stored in diagonal storage format
 */
void eihMPC_LA_DIAG_CHOL_ONELOOP_LBUB_14_14_14(eihMPC_FLOAT *H, eihMPC_FLOAT *llbysl, int* lbIdx, eihMPC_FLOAT *lubysu, int* ubIdx, eihMPC_FLOAT *Phi)


{
	int i;
	
	/* compute cholesky */
	for( i=0; i<14; i++ ){
		Phi[i] = H[i] + llbysl[i] + lubysu[i];

#if eihMPC_SET_PRINTLEVEL > 0 && defined PRINTNUMERICALWARNINGS
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
 * where A is to be computed and is of size [7 x 14],
 * B is given and of size [7 x 14], L is a diagonal
 * matrix of size 7 stored in diagonal matrix 
 * storage format. Note the transpose of L has no impact!
 *
 * Result: A in column major storage format.
 *
 */
void eihMPC_LA_DIAG_MATRIXFORWARDSUB_7_14(eihMPC_FLOAT *L, eihMPC_FLOAT *B, eihMPC_FLOAT *A)
{
    int i,j;
	 int k = 0;

	for( j=0; j<14; j++){
		for( i=0; i<7; i++){
			A[k] = B[k]/L[j];
			k++;
		}
	}

}


/**
 * Forward substitution for the matrix equation A*L' = B
 * where A is to be computed and is of size [7 x 14],
 * B is given and of size [7 x 14], L is a diagonal
 *  matrix of size 14 stored in diagonal 
 * storage format. Note the transpose of L!
 *
 * Result: A in diagonalzero storage format.
 *
 */
void eihMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_7_14(eihMPC_FLOAT *L, eihMPC_FLOAT *B, eihMPC_FLOAT *A)
{
	int j;
    for( j=0; j<14; j++ ){   
		A[j] = B[j]/L[j];
     }
}


/**
 * Compute C = A*B' where 
 *
 *	size(A) = [7 x 14]
 *  size(B) = [7 x 14] in diagzero format
 * 
 * A and C matrices are stored in column major format.
 * 
 * 
 */
void eihMPC_LA_DENSE_DIAGZERO_MMTM_7_14_7(eihMPC_FLOAT *A, eihMPC_FLOAT *B, eihMPC_FLOAT *C)
{
    int i, j;
	
	for( i=0; i<7; i++ ){
		for( j=0; j<7; j++){
			C[j*7+i] = B[i*7+j]*A[i];
		}
	}

}


/**
 * Forward substitution to solve L*y = b where L is a
 * diagonal matrix in vector storage format.
 * 
 * The dimensions involved are 14.
 */
void eihMPC_LA_DIAG_FORWARDSUB_14(eihMPC_FLOAT *L, eihMPC_FLOAT *b, eihMPC_FLOAT *y)
{
    int i;

    for( i=0; i<14; i++ ){
		y[i] = b[i]/L[i];
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
void eihMPC_LA_DIAG_CHOL_ONELOOP_LBUB_7_7_7(eihMPC_FLOAT *H, eihMPC_FLOAT *llbysl, int* lbIdx, eihMPC_FLOAT *lubysu, int* ubIdx, eihMPC_FLOAT *Phi)


{
	int i;
	
	/* compute cholesky */
	for( i=0; i<7; i++ ){
		Phi[i] = H[i] + llbysl[i] + lubysu[i];

#if eihMPC_SET_PRINTLEVEL > 0 && defined PRINTNUMERICALWARNINGS
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
 * where A is to be computed and is of size [7 x 7],
 * B is given and of size [7 x 7], L is a diagonal
 *  matrix of size 7 stored in diagonal 
 * storage format. Note the transpose of L!
 *
 * Result: A in diagonalzero storage format.
 *
 */
void eihMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_7_7(eihMPC_FLOAT *L, eihMPC_FLOAT *B, eihMPC_FLOAT *A)
{
	int j;
    for( j=0; j<7; j++ ){   
		A[j] = B[j]/L[j];
     }
}


/**
 * Forward substitution to solve L*y = b where L is a
 * diagonal matrix in vector storage format.
 * 
 * The dimensions involved are 7.
 */
void eihMPC_LA_DIAG_FORWARDSUB_7(eihMPC_FLOAT *L, eihMPC_FLOAT *b, eihMPC_FLOAT *y)
{
    int i;

    for( i=0; i<7; i++ ){
		y[i] = b[i]/L[i];
    }
}


/**
 * Compute L = B*B', where L is lower triangular of size NXp1
 * and B is a diagzero matrix of size [7 x 14] in column
 * storage format.
 * 
 */
void eihMPC_LA_DIAGZERO_MMT_7(eihMPC_FLOAT *B, eihMPC_FLOAT *L)
{
    int i, ii, di;
    
    ii = 0; di = 0;
    for( i=0; i<7; i++ ){        
		L[ii+i] = B[i]*B[i];
        ii += ++di;
    }
}


/* 
 * Computes r = b - B*u
 * B is stored in diagzero format
 */
void eihMPC_LA_DIAGZERO_MVMSUB7_7(eihMPC_FLOAT *B, eihMPC_FLOAT *u, eihMPC_FLOAT *b, eihMPC_FLOAT *r)
{
	int i;

	for( i=0; i<7; i++ ){
		r[i] = b[i] - B[i]*u[i];
	}	
	
}


/**
 * Compute L = A*A' + B*B', where L is lower triangular of size NXp1
 * and A is a dense matrix of size [7 x 14] in column
 * storage format, and B is of size [7 x 14] diagonalzero
 * storage format.
 * 
 * THIS ONE HAS THE WORST ACCES PATTERN POSSIBLE. 
 * POSSIBKE FIX: PUT A AND B INTO ROW MAJOR FORMAT FIRST.
 * 
 */
void eihMPC_LA_DENSE_DIAGZERO_MMT2_7_14_14(eihMPC_FLOAT *A, eihMPC_FLOAT *B, eihMPC_FLOAT *L)
{
    int i, j, k, ii, di;
    eihMPC_FLOAT ltemp;
    
    ii = 0; di = 0;
    for( i=0; i<7; i++ ){        
        for( j=0; j<=i; j++ ){
            ltemp = 0; 
            for( k=0; k<14; k++ ){
                ltemp += A[k*7+i]*A[k*7+j];
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
void eihMPC_LA_DENSE_DIAGZERO_2MVMSUB2_7_14_14(eihMPC_FLOAT *A, eihMPC_FLOAT *x, eihMPC_FLOAT *B, eihMPC_FLOAT *u, eihMPC_FLOAT *b, eihMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;

	for( i=0; i<7; i++ ){
		r[i] = b[i] - A[k++]*x[0] - B[i]*u[i];
	}	

	for( j=1; j<14; j++ ){		
		for( i=0; i<7; i++ ){
			r[i] -= A[k++]*x[j];
		}
	}
	
}


/**
 * Compute L = A*A' + B*B', where L is lower triangular of size NXp1
 * and A is a dense matrix of size [7 x 14] in column
 * storage format, and B is of size [7 x 7] diagonalzero
 * storage format.
 * 
 * THIS ONE HAS THE WORST ACCES PATTERN POSSIBLE. 
 * POSSIBKE FIX: PUT A AND B INTO ROW MAJOR FORMAT FIRST.
 * 
 */
void eihMPC_LA_DENSE_DIAGZERO_MMT2_7_14_7(eihMPC_FLOAT *A, eihMPC_FLOAT *B, eihMPC_FLOAT *L)
{
    int i, j, k, ii, di;
    eihMPC_FLOAT ltemp;
    
    ii = 0; di = 0;
    for( i=0; i<7; i++ ){        
        for( j=0; j<=i; j++ ){
            ltemp = 0; 
            for( k=0; k<14; k++ ){
                ltemp += A[k*7+i]*A[k*7+j];
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
void eihMPC_LA_DENSE_DIAGZERO_2MVMSUB2_7_14_7(eihMPC_FLOAT *A, eihMPC_FLOAT *x, eihMPC_FLOAT *B, eihMPC_FLOAT *u, eihMPC_FLOAT *b, eihMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;

	for( i=0; i<7; i++ ){
		r[i] = b[i] - A[k++]*x[0] - B[i]*u[i];
	}	

	for( j=1; j<14; j++ ){		
		for( i=0; i<7; i++ ){
			r[i] -= A[k++]*x[j];
		}
	}
	
}


/**
 * Cholesky factorization as above, but working on a matrix in 
 * lower triangular storage format of size 7 and outputting
 * the Cholesky factor to matrix L in lower triangular format.
 */
void eihMPC_LA_DENSE_CHOL_7(eihMPC_FLOAT *A, eihMPC_FLOAT *L)
{
    int i, j, k, di, dj;
	 int ii, jj;

    eihMPC_FLOAT l;
    eihMPC_FLOAT Mii;

	/* copy A to L first and then operate on L */
	/* COULD BE OPTIMIZED */
	ii=0; di=0;
	for( i=0; i<7; i++ ){
		for( j=0; j<=i; j++ ){
			L[ii+j] = A[ii+j];
		}
		ii += ++di;
	}    
	
	/* factor L */
	ii=0; di=0;
    for( i=0; i<7; i++ ){
        l = 0;
        for( k=0; k<i; k++ ){
            l += L[ii+k]*L[ii+k];
        }        
        
        Mii = L[ii+i] - l;
        
#if eihMPC_SET_PRINTLEVEL > 0 && defined PRINTNUMERICALWARNINGS
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
        for( j=i+1; j<7; j++ ){
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
 * The dimensions involved are 7.
 */
void eihMPC_LA_DENSE_FORWARDSUB_7(eihMPC_FLOAT *L, eihMPC_FLOAT *b, eihMPC_FLOAT *y)
{
    int i,j,ii,di;
    eihMPC_FLOAT yel;
            
    ii = 0; di = 0;
    for( i=0; i<7; i++ ){
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
 * where A is to be computed and is of size [7 x 7],
 * B is given and of size [7 x 7], L is a lower tri-
 * angular matrix of size 7 stored in lower triangular 
 * storage format. Note the transpose of L AND B!
 *
 * Result: A in column major storage format.
 *
 */
void eihMPC_LA_DENSE_MATRIXTFORWARDSUB_7_7(eihMPC_FLOAT *L, eihMPC_FLOAT *B, eihMPC_FLOAT *A)
{
    int i,j,k,ii,di;
    eihMPC_FLOAT a;
    
    ii=0; di=0;
    for( j=0; j<7; j++ ){        
        for( i=0; i<7; i++ ){
            a = B[i*7+j];
            for( k=0; k<j; k++ ){
                a -= A[k*7+i]*L[ii+k];
            }    

			/* saturate for numerical stability */
			a = MIN(a, BIGM);
			a = MAX(a, -BIGM); 

			A[j*7+i] = a/L[ii+j];			
        }
        ii += ++di;
    }
}


/**
 * Compute L = L - A*A', where L is lower triangular of size 7
 * and A is a dense matrix of size [7 x 7] in column
 * storage format.
 * 
 * THIS ONE HAS THE WORST ACCES PATTERN POSSIBLE. 
 * POSSIBKE FIX: PUT A INTO ROW MAJOR FORMAT FIRST.
 * 
 */
void eihMPC_LA_DENSE_MMTSUB_7_7(eihMPC_FLOAT *A, eihMPC_FLOAT *L)
{
    int i, j, k, ii, di;
    eihMPC_FLOAT ltemp;
    
    ii = 0; di = 0;
    for( i=0; i<7; i++ ){        
        for( j=0; j<=i; j++ ){
            ltemp = 0; 
            for( k=0; k<7; k++ ){
                ltemp += A[k*7+i]*A[k*7+j];
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
void eihMPC_LA_DENSE_MVMSUB1_7_7(eihMPC_FLOAT *A, eihMPC_FLOAT *x, eihMPC_FLOAT *b, eihMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;

	for( i=0; i<7; i++ ){
		r[i] = b[i] - A[k++]*x[0];
	}	
	for( j=1; j<7; j++ ){		
		for( i=0; i<7; i++ ){
			r[i] -= A[k++]*x[j];
		}
	}
}


/**
 * Backward Substitution to solve L^T*x = y where L is a
 * lower triangular matrix in triangular storage format.
 * 
 * All involved dimensions are 7.
 */
void eihMPC_LA_DENSE_BACKWARDSUB_7(eihMPC_FLOAT *L, eihMPC_FLOAT *y, eihMPC_FLOAT *x)
{
    int i, ii, di, j, jj, dj;
    eihMPC_FLOAT xel;    
	int start = 21;
    
    /* now solve L^T*x = y by backward substitution */
    ii = start; di = 6;
    for( i=6; i>=0; i-- ){        
        xel = y[i];        
        jj = start; dj = 6;
        for( j=6; j>i; j-- ){
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
 * Matrix vector multiplication y = b - M'*x where M is of size [7 x 7]
 * and stored in column major format. Note the transpose of M!
 */
void eihMPC_LA_DENSE_MTVMSUB_7_7(eihMPC_FLOAT *A, eihMPC_FLOAT *x, eihMPC_FLOAT *b, eihMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0; 
	for( i=0; i<7; i++ ){
		r[i] = b[i];
		for( j=0; j<7; j++ ){
			r[i] -= A[k++]*x[j];
		}
	}
}


/*
 * Vector subtraction z = -x - y for vectors of length 133.
 */
void eihMPC_LA_VSUB2_133(eihMPC_FLOAT *x, eihMPC_FLOAT *y, eihMPC_FLOAT *z)
{
	int i;
	for( i=0; i<133; i++){
		z[i] = -x[i] - y[i];
	}
}


/**
 * Forward-Backward-Substitution to solve L*L^T*x = b where L is a
 * diagonal matrix of size 14 in vector
 * storage format.
 */
void eihMPC_LA_DIAG_FORWARDBACKWARDSUB_14(eihMPC_FLOAT *L, eihMPC_FLOAT *b, eihMPC_FLOAT *x)
{
    int i;
            
    /* solve Ly = b by forward and backward substitution */
    for( i=0; i<14; i++ ){
		x[i] = b[i]/(L[i]*L[i]);
    }
    
}


/**
 * Forward-Backward-Substitution to solve L*L^T*x = b where L is a
 * diagonal matrix of size 7 in vector
 * storage format.
 */
void eihMPC_LA_DIAG_FORWARDBACKWARDSUB_7(eihMPC_FLOAT *L, eihMPC_FLOAT *b, eihMPC_FLOAT *x)
{
    int i;
            
    /* solve Ly = b by forward and backward substitution */
    for( i=0; i<7; i++ ){
		x[i] = b[i]/(L[i]*L[i]);
    }
    
}


/*
 * Vector subtraction z = x(xidx) - y where y, z and xidx are of length 14,
 * and x has length 14 and is indexed through yidx.
 */
void eihMPC_LA_VSUB_INDEXED_14(eihMPC_FLOAT *x, int* xidx, eihMPC_FLOAT *y, eihMPC_FLOAT *z)
{
	int i;
	for( i=0; i<14; i++){
		z[i] = x[xidx[i]] - y[i];
	}
}


/*
 * Vector subtraction x = -u.*v - w for vectors of length 14.
 */
void eihMPC_LA_VSUB3_14(eihMPC_FLOAT *u, eihMPC_FLOAT *v, eihMPC_FLOAT *w, eihMPC_FLOAT *x)
{
	int i;
	for( i=0; i<14; i++){
		x[i] = -u[i]*v[i] - w[i];
	}
}


/*
 * Vector subtraction z = -x - y(yidx) where y is of length 14
 * and z, x and yidx are of length 14.
 */
void eihMPC_LA_VSUB2_INDEXED_14(eihMPC_FLOAT *x, eihMPC_FLOAT *y, int* yidx, eihMPC_FLOAT *z)
{
	int i;
	for( i=0; i<14; i++){
		z[i] = -x[i] - y[yidx[i]];
	}
}


/*
 * Vector subtraction z = x(xidx) - y where y, z and xidx are of length 7,
 * and x has length 7 and is indexed through yidx.
 */
void eihMPC_LA_VSUB_INDEXED_7(eihMPC_FLOAT *x, int* xidx, eihMPC_FLOAT *y, eihMPC_FLOAT *z)
{
	int i;
	for( i=0; i<7; i++){
		z[i] = x[xidx[i]] - y[i];
	}
}


/*
 * Vector subtraction x = -u.*v - w for vectors of length 7.
 */
void eihMPC_LA_VSUB3_7(eihMPC_FLOAT *u, eihMPC_FLOAT *v, eihMPC_FLOAT *w, eihMPC_FLOAT *x)
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
void eihMPC_LA_VSUB2_INDEXED_7(eihMPC_FLOAT *x, eihMPC_FLOAT *y, int* yidx, eihMPC_FLOAT *z)
{
	int i;
	for( i=0; i<7; i++){
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
 * eihMPC_NOPROGRESS (should be negative).
 */
int eihMPC_LINESEARCH_BACKTRACKING_AFFINE(eihMPC_FLOAT *l, eihMPC_FLOAT *s, eihMPC_FLOAT *dl, eihMPC_FLOAT *ds, eihMPC_FLOAT *a, eihMPC_FLOAT *mu_aff)
{
    int i;
	int lsIt=1;    
    eihMPC_FLOAT dltemp;
    eihMPC_FLOAT dstemp;
    eihMPC_FLOAT mya = 1.0;
    eihMPC_FLOAT mymu;
        
    while( 1 ){                        

        /* 
         * Compute both snew and wnew together.
         * We compute also mu_affine along the way here, as the
         * values might be in registers, so it should be cheaper.
         */
        mymu = 0;
        for( i=0; i<266; i++ ){
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
        if( i == 266 ){
            break;
        } else {
            mya *= eihMPC_SET_LS_SCALE_AFF;
            if( mya < eihMPC_SET_LS_MINSTEP ){
                return eihMPC_NOPROGRESS;
            }
        }
    }
    
    /* return new values and iteration counter */
    *a = mya;
    *mu_aff = mymu / (eihMPC_FLOAT)266;
    return lsIt;
}


/*
 * Vector subtraction x = u.*v - a where a is a scalar
*  and x,u,v are vectors of length 266.
 */
void eihMPC_LA_VSUB5_266(eihMPC_FLOAT *u, eihMPC_FLOAT *v, eihMPC_FLOAT a, eihMPC_FLOAT *x)
{
	int i;
	for( i=0; i<266; i++){
		x[i] = u[i]*v[i] - a;
	}
}


/*
 * Computes x=0; x(uidx) += u/su; x(vidx) -= v/sv where x is of length 14,
 * u, su, uidx are of length 14 and v, sv, vidx are of length 14.
 */
void eihMPC_LA_VSUB6_INDEXED_14_14_14(eihMPC_FLOAT *u, eihMPC_FLOAT *su, int* uidx, eihMPC_FLOAT *v, eihMPC_FLOAT *sv, int* vidx, eihMPC_FLOAT *x)
{
	int i;
	for( i=0; i<14; i++ ){
		x[i] = 0;
	}
	for( i=0; i<14; i++){
		x[uidx[i]] += u[i]/su[i];
	}
	for( i=0; i<14; i++){
		x[vidx[i]] -= v[i]/sv[i];
	}
}


/* 
 * Computes r =  B*u
 * where B is stored in diagzero format
 */
void eihMPC_LA_DIAGZERO_MVM_7(eihMPC_FLOAT *B, eihMPC_FLOAT *u, eihMPC_FLOAT *r)
{
	int i;

	for( i=0; i<7; i++ ){
		r[i] = B[i]*u[i];
	}	
	
}


/* 
 * Computes r = A*x + B*u
 * where A is stored in column major format
 * and B is stored in diagzero format
 */
void eihMPC_LA_DENSE_DIAGZERO_2MVMADD_7_14_14(eihMPC_FLOAT *A, eihMPC_FLOAT *x, eihMPC_FLOAT *B, eihMPC_FLOAT *u, eihMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;

	for( i=0; i<7; i++ ){
		r[i] = A[k++]*x[0] + B[i]*u[i];
	}	

	for( j=1; j<14; j++ ){		
		for( i=0; i<7; i++ ){
			r[i] += A[k++]*x[j];
		}
	}
	
}


/*
 * Computes x=0; x(uidx) += u/su; x(vidx) -= v/sv where x is of length 7,
 * u, su, uidx are of length 7 and v, sv, vidx are of length 7.
 */
void eihMPC_LA_VSUB6_INDEXED_7_7_7(eihMPC_FLOAT *u, eihMPC_FLOAT *su, int* uidx, eihMPC_FLOAT *v, eihMPC_FLOAT *sv, int* vidx, eihMPC_FLOAT *x)
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
 * where A is stored in column major format
 * and B is stored in diagzero format
 */
void eihMPC_LA_DENSE_DIAGZERO_2MVMADD_7_14_7(eihMPC_FLOAT *A, eihMPC_FLOAT *x, eihMPC_FLOAT *B, eihMPC_FLOAT *u, eihMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;

	for( i=0; i<7; i++ ){
		r[i] = A[k++]*x[0] + B[i]*u[i];
	}	

	for( j=1; j<14; j++ ){		
		for( i=0; i<7; i++ ){
			r[i] += A[k++]*x[j];
		}
	}
	
}


/*
 * Vector subtraction z = x - y for vectors of length 133.
 */
void eihMPC_LA_VSUB_133(eihMPC_FLOAT *x, eihMPC_FLOAT *y, eihMPC_FLOAT *z)
{
	int i;
	for( i=0; i<133; i++){
		z[i] = x[i] - y[i];
	}
}


/** 
 * Computes z = -r./s - u.*y(y)
 * where all vectors except of y are of length 14 (length of y >= 14).
 */
void eihMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_14(eihMPC_FLOAT *r, eihMPC_FLOAT *s, eihMPC_FLOAT *u, eihMPC_FLOAT *y, int* yidx, eihMPC_FLOAT *z)
{
	int i;
	for( i=0; i<14; i++ ){
		z[i] = -r[i]/s[i] - u[i]*y[yidx[i]];
	}
}


/** 
 * Computes z = -r./s + u.*y(y)
 * where all vectors except of y are of length 14 (length of y >= 14).
 */
void eihMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_14(eihMPC_FLOAT *r, eihMPC_FLOAT *s, eihMPC_FLOAT *u, eihMPC_FLOAT *y, int* yidx, eihMPC_FLOAT *z)
{
	int i;
	for( i=0; i<14; i++ ){
		z[i] = -r[i]/s[i] + u[i]*y[yidx[i]];
	}
}


/** 
 * Computes z = -r./s - u.*y(y)
 * where all vectors except of y are of length 7 (length of y >= 7).
 */
void eihMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_7(eihMPC_FLOAT *r, eihMPC_FLOAT *s, eihMPC_FLOAT *u, eihMPC_FLOAT *y, int* yidx, eihMPC_FLOAT *z)
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
void eihMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_7(eihMPC_FLOAT *r, eihMPC_FLOAT *s, eihMPC_FLOAT *u, eihMPC_FLOAT *y, int* yidx, eihMPC_FLOAT *z)
{
	int i;
	for( i=0; i<7; i++ ){
		z[i] = -r[i]/s[i] + u[i]*y[yidx[i]];
	}
}


/*
 * Computes ds = -l.\(r + s.*dl) for vectors of length 266.
 */
void eihMPC_LA_VSUB7_266(eihMPC_FLOAT *l, eihMPC_FLOAT *r, eihMPC_FLOAT *s, eihMPC_FLOAT *dl, eihMPC_FLOAT *ds)
{
	int i;
	for( i=0; i<266; i++){
		ds[i] = -(r[i] + s[i]*dl[i])/l[i];
	}
}


/*
 * Vector addition x = x + y for vectors of length 133.
 */
void eihMPC_LA_VADD_133(eihMPC_FLOAT *x, eihMPC_FLOAT *y)
{
	int i;
	for( i=0; i<133; i++){
		x[i] += y[i];
	}
}


/*
 * Vector addition x = x + y for vectors of length 70.
 */
void eihMPC_LA_VADD_70(eihMPC_FLOAT *x, eihMPC_FLOAT *y)
{
	int i;
	for( i=0; i<70; i++){
		x[i] += y[i];
	}
}


/*
 * Vector addition x = x + y for vectors of length 266.
 */
void eihMPC_LA_VADD_266(eihMPC_FLOAT *x, eihMPC_FLOAT *y)
{
	int i;
	for( i=0; i<266; i++){
		x[i] += y[i];
	}
}


/**
 * Backtracking line search for combined predictor/corrector step.
 * Update on variables with safety factor gamma (to keep us away from
 * boundary).
 */
int eihMPC_LINESEARCH_BACKTRACKING_COMBINED(eihMPC_FLOAT *z, eihMPC_FLOAT *v, eihMPC_FLOAT *l, eihMPC_FLOAT *s, eihMPC_FLOAT *dz, eihMPC_FLOAT *dv, eihMPC_FLOAT *dl, eihMPC_FLOAT *ds, eihMPC_FLOAT *a, eihMPC_FLOAT *mu)
{
    int i, lsIt=1;       
    eihMPC_FLOAT dltemp;
    eihMPC_FLOAT dstemp;    
    eihMPC_FLOAT a_gamma;
            
    *a = 1.0;
    while( 1 ){                        

        /* check whether search criterion is fulfilled */
        for( i=0; i<266; i++ ){
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
        if( i == 266 ){
            break;
        } else {
            *a *= eihMPC_SET_LS_SCALE;
            if( *a < eihMPC_SET_LS_MINSTEP ){
                return eihMPC_NOPROGRESS;
            }
        }
    }
    
    /* update variables with safety margin */
    a_gamma = (*a)*eihMPC_SET_LS_MAXSTEP;
    
    /* primal variables */
    for( i=0; i<133; i++ ){
        z[i] += a_gamma*dz[i];
    }
    
    /* equality constraint multipliers */
    for( i=0; i<70; i++ ){
        v[i] += a_gamma*dv[i];
    }
    
    /* inequality constraint multipliers & slacks, also update mu */
    *mu = 0;
    for( i=0; i<266; i++ ){
        dltemp = l[i] + a_gamma*dl[i]; l[i] = dltemp;
        dstemp = s[i] + a_gamma*ds[i]; s[i] = dstemp;
        *mu += dltemp*dstemp;
    }
    
    *a = a_gamma;
    *mu /= (eihMPC_FLOAT)266;
    return lsIt;
}




/* VARIABLE DEFINITIONS ------------------------------------------------ */
eihMPC_FLOAT eihMPC_z[133];
eihMPC_FLOAT eihMPC_v[70];
eihMPC_FLOAT eihMPC_dz_aff[133];
eihMPC_FLOAT eihMPC_dv_aff[70];
eihMPC_FLOAT eihMPC_grad_cost[133];
eihMPC_FLOAT eihMPC_grad_eq[133];
eihMPC_FLOAT eihMPC_rd[133];
eihMPC_FLOAT eihMPC_l[266];
eihMPC_FLOAT eihMPC_s[266];
eihMPC_FLOAT eihMPC_lbys[266];
eihMPC_FLOAT eihMPC_dl_aff[266];
eihMPC_FLOAT eihMPC_ds_aff[266];
eihMPC_FLOAT eihMPC_dz_cc[133];
eihMPC_FLOAT eihMPC_dv_cc[70];
eihMPC_FLOAT eihMPC_dl_cc[266];
eihMPC_FLOAT eihMPC_ds_cc[266];
eihMPC_FLOAT eihMPC_ccrhs[266];
eihMPC_FLOAT eihMPC_grad_ineq[133];
eihMPC_FLOAT* eihMPC_z0 = eihMPC_z + 0;
eihMPC_FLOAT* eihMPC_dzaff0 = eihMPC_dz_aff + 0;
eihMPC_FLOAT* eihMPC_dzcc0 = eihMPC_dz_cc + 0;
eihMPC_FLOAT* eihMPC_rd0 = eihMPC_rd + 0;
eihMPC_FLOAT eihMPC_Lbyrd0[14];
eihMPC_FLOAT* eihMPC_grad_cost0 = eihMPC_grad_cost + 0;
eihMPC_FLOAT* eihMPC_grad_eq0 = eihMPC_grad_eq + 0;
eihMPC_FLOAT* eihMPC_grad_ineq0 = eihMPC_grad_ineq + 0;
eihMPC_FLOAT eihMPC_ctv0[14];
eihMPC_FLOAT eihMPC_C0[98] = {1.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 
0.0000000000000000E+000, 1.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 
0.0000000000000000E+000, 0.0000000000000000E+000, 1.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 
0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 1.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 
0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 1.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 
0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 1.0000000000000000E+000, 0.0000000000000000E+000, 
0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 1.0000000000000000E+000, 
1.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 
0.0000000000000000E+000, 1.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 
0.0000000000000000E+000, 0.0000000000000000E+000, 1.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 
0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 1.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 
0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 1.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 
0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 1.0000000000000000E+000, 0.0000000000000000E+000, 
0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 1.0000000000000000E+000};
eihMPC_FLOAT* eihMPC_v0 = eihMPC_v + 0;
eihMPC_FLOAT eihMPC_re0[7];
eihMPC_FLOAT eihMPC_beta0[7];
eihMPC_FLOAT eihMPC_betacc0[7];
eihMPC_FLOAT* eihMPC_dvaff0 = eihMPC_dv_aff + 0;
eihMPC_FLOAT* eihMPC_dvcc0 = eihMPC_dv_cc + 0;
eihMPC_FLOAT eihMPC_V0[98];
eihMPC_FLOAT eihMPC_Yd0[28];
eihMPC_FLOAT eihMPC_Ld0[28];
eihMPC_FLOAT eihMPC_yy0[7];
eihMPC_FLOAT eihMPC_bmy0[7];
int eihMPC_lbIdx0[14] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13};
eihMPC_FLOAT* eihMPC_llb0 = eihMPC_l + 0;
eihMPC_FLOAT* eihMPC_slb0 = eihMPC_s + 0;
eihMPC_FLOAT* eihMPC_llbbyslb0 = eihMPC_lbys + 0;
eihMPC_FLOAT eihMPC_rilb0[14];
eihMPC_FLOAT* eihMPC_dllbaff0 = eihMPC_dl_aff + 0;
eihMPC_FLOAT* eihMPC_dslbaff0 = eihMPC_ds_aff + 0;
eihMPC_FLOAT* eihMPC_dllbcc0 = eihMPC_dl_cc + 0;
eihMPC_FLOAT* eihMPC_dslbcc0 = eihMPC_ds_cc + 0;
eihMPC_FLOAT* eihMPC_ccrhsl0 = eihMPC_ccrhs + 0;
int eihMPC_ubIdx0[14] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13};
eihMPC_FLOAT* eihMPC_lub0 = eihMPC_l + 14;
eihMPC_FLOAT* eihMPC_sub0 = eihMPC_s + 14;
eihMPC_FLOAT* eihMPC_lubbysub0 = eihMPC_lbys + 14;
eihMPC_FLOAT eihMPC_riub0[14];
eihMPC_FLOAT* eihMPC_dlubaff0 = eihMPC_dl_aff + 14;
eihMPC_FLOAT* eihMPC_dsubaff0 = eihMPC_ds_aff + 14;
eihMPC_FLOAT* eihMPC_dlubcc0 = eihMPC_dl_cc + 14;
eihMPC_FLOAT* eihMPC_dsubcc0 = eihMPC_ds_cc + 14;
eihMPC_FLOAT* eihMPC_ccrhsub0 = eihMPC_ccrhs + 14;
eihMPC_FLOAT eihMPC_Phi0[14];
eihMPC_FLOAT eihMPC_D0[14] = {1.0000000000000000E+000, 
1.0000000000000000E+000, 
1.0000000000000000E+000, 
1.0000000000000000E+000, 
1.0000000000000000E+000, 
1.0000000000000000E+000, 
1.0000000000000000E+000};
eihMPC_FLOAT eihMPC_W0[14];
eihMPC_FLOAT* eihMPC_z1 = eihMPC_z + 14;
eihMPC_FLOAT* eihMPC_dzaff1 = eihMPC_dz_aff + 14;
eihMPC_FLOAT* eihMPC_dzcc1 = eihMPC_dz_cc + 14;
eihMPC_FLOAT* eihMPC_rd1 = eihMPC_rd + 14;
eihMPC_FLOAT eihMPC_Lbyrd1[14];
eihMPC_FLOAT* eihMPC_grad_cost1 = eihMPC_grad_cost + 14;
eihMPC_FLOAT* eihMPC_grad_eq1 = eihMPC_grad_eq + 14;
eihMPC_FLOAT* eihMPC_grad_ineq1 = eihMPC_grad_ineq + 14;
eihMPC_FLOAT eihMPC_ctv1[14];
eihMPC_FLOAT* eihMPC_v1 = eihMPC_v + 7;
eihMPC_FLOAT eihMPC_re1[7];
eihMPC_FLOAT eihMPC_beta1[7];
eihMPC_FLOAT eihMPC_betacc1[7];
eihMPC_FLOAT* eihMPC_dvaff1 = eihMPC_dv_aff + 7;
eihMPC_FLOAT* eihMPC_dvcc1 = eihMPC_dv_cc + 7;
eihMPC_FLOAT eihMPC_V1[98];
eihMPC_FLOAT eihMPC_Yd1[28];
eihMPC_FLOAT eihMPC_Ld1[28];
eihMPC_FLOAT eihMPC_yy1[7];
eihMPC_FLOAT eihMPC_bmy1[7];
eihMPC_FLOAT eihMPC_c1[7] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
int eihMPC_lbIdx1[14] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13};
eihMPC_FLOAT* eihMPC_llb1 = eihMPC_l + 28;
eihMPC_FLOAT* eihMPC_slb1 = eihMPC_s + 28;
eihMPC_FLOAT* eihMPC_llbbyslb1 = eihMPC_lbys + 28;
eihMPC_FLOAT eihMPC_rilb1[14];
eihMPC_FLOAT* eihMPC_dllbaff1 = eihMPC_dl_aff + 28;
eihMPC_FLOAT* eihMPC_dslbaff1 = eihMPC_ds_aff + 28;
eihMPC_FLOAT* eihMPC_dllbcc1 = eihMPC_dl_cc + 28;
eihMPC_FLOAT* eihMPC_dslbcc1 = eihMPC_ds_cc + 28;
eihMPC_FLOAT* eihMPC_ccrhsl1 = eihMPC_ccrhs + 28;
int eihMPC_ubIdx1[14] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13};
eihMPC_FLOAT* eihMPC_lub1 = eihMPC_l + 42;
eihMPC_FLOAT* eihMPC_sub1 = eihMPC_s + 42;
eihMPC_FLOAT* eihMPC_lubbysub1 = eihMPC_lbys + 42;
eihMPC_FLOAT eihMPC_riub1[14];
eihMPC_FLOAT* eihMPC_dlubaff1 = eihMPC_dl_aff + 42;
eihMPC_FLOAT* eihMPC_dsubaff1 = eihMPC_ds_aff + 42;
eihMPC_FLOAT* eihMPC_dlubcc1 = eihMPC_dl_cc + 42;
eihMPC_FLOAT* eihMPC_dsubcc1 = eihMPC_ds_cc + 42;
eihMPC_FLOAT* eihMPC_ccrhsub1 = eihMPC_ccrhs + 42;
eihMPC_FLOAT eihMPC_Phi1[14];
eihMPC_FLOAT eihMPC_D1[14] = {-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000};
eihMPC_FLOAT eihMPC_W1[14];
eihMPC_FLOAT eihMPC_Ysd1[49];
eihMPC_FLOAT eihMPC_Lsd1[49];
eihMPC_FLOAT* eihMPC_z2 = eihMPC_z + 28;
eihMPC_FLOAT* eihMPC_dzaff2 = eihMPC_dz_aff + 28;
eihMPC_FLOAT* eihMPC_dzcc2 = eihMPC_dz_cc + 28;
eihMPC_FLOAT* eihMPC_rd2 = eihMPC_rd + 28;
eihMPC_FLOAT eihMPC_Lbyrd2[14];
eihMPC_FLOAT* eihMPC_grad_cost2 = eihMPC_grad_cost + 28;
eihMPC_FLOAT* eihMPC_grad_eq2 = eihMPC_grad_eq + 28;
eihMPC_FLOAT* eihMPC_grad_ineq2 = eihMPC_grad_ineq + 28;
eihMPC_FLOAT eihMPC_ctv2[14];
eihMPC_FLOAT* eihMPC_v2 = eihMPC_v + 14;
eihMPC_FLOAT eihMPC_re2[7];
eihMPC_FLOAT eihMPC_beta2[7];
eihMPC_FLOAT eihMPC_betacc2[7];
eihMPC_FLOAT* eihMPC_dvaff2 = eihMPC_dv_aff + 14;
eihMPC_FLOAT* eihMPC_dvcc2 = eihMPC_dv_cc + 14;
eihMPC_FLOAT eihMPC_V2[98];
eihMPC_FLOAT eihMPC_Yd2[28];
eihMPC_FLOAT eihMPC_Ld2[28];
eihMPC_FLOAT eihMPC_yy2[7];
eihMPC_FLOAT eihMPC_bmy2[7];
eihMPC_FLOAT eihMPC_c2[7] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
int eihMPC_lbIdx2[14] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13};
eihMPC_FLOAT* eihMPC_llb2 = eihMPC_l + 56;
eihMPC_FLOAT* eihMPC_slb2 = eihMPC_s + 56;
eihMPC_FLOAT* eihMPC_llbbyslb2 = eihMPC_lbys + 56;
eihMPC_FLOAT eihMPC_rilb2[14];
eihMPC_FLOAT* eihMPC_dllbaff2 = eihMPC_dl_aff + 56;
eihMPC_FLOAT* eihMPC_dslbaff2 = eihMPC_ds_aff + 56;
eihMPC_FLOAT* eihMPC_dllbcc2 = eihMPC_dl_cc + 56;
eihMPC_FLOAT* eihMPC_dslbcc2 = eihMPC_ds_cc + 56;
eihMPC_FLOAT* eihMPC_ccrhsl2 = eihMPC_ccrhs + 56;
int eihMPC_ubIdx2[14] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13};
eihMPC_FLOAT* eihMPC_lub2 = eihMPC_l + 70;
eihMPC_FLOAT* eihMPC_sub2 = eihMPC_s + 70;
eihMPC_FLOAT* eihMPC_lubbysub2 = eihMPC_lbys + 70;
eihMPC_FLOAT eihMPC_riub2[14];
eihMPC_FLOAT* eihMPC_dlubaff2 = eihMPC_dl_aff + 70;
eihMPC_FLOAT* eihMPC_dsubaff2 = eihMPC_ds_aff + 70;
eihMPC_FLOAT* eihMPC_dlubcc2 = eihMPC_dl_cc + 70;
eihMPC_FLOAT* eihMPC_dsubcc2 = eihMPC_ds_cc + 70;
eihMPC_FLOAT* eihMPC_ccrhsub2 = eihMPC_ccrhs + 70;
eihMPC_FLOAT eihMPC_Phi2[14];
eihMPC_FLOAT eihMPC_W2[14];
eihMPC_FLOAT eihMPC_Ysd2[49];
eihMPC_FLOAT eihMPC_Lsd2[49];
eihMPC_FLOAT* eihMPC_z3 = eihMPC_z + 42;
eihMPC_FLOAT* eihMPC_dzaff3 = eihMPC_dz_aff + 42;
eihMPC_FLOAT* eihMPC_dzcc3 = eihMPC_dz_cc + 42;
eihMPC_FLOAT* eihMPC_rd3 = eihMPC_rd + 42;
eihMPC_FLOAT eihMPC_Lbyrd3[14];
eihMPC_FLOAT* eihMPC_grad_cost3 = eihMPC_grad_cost + 42;
eihMPC_FLOAT* eihMPC_grad_eq3 = eihMPC_grad_eq + 42;
eihMPC_FLOAT* eihMPC_grad_ineq3 = eihMPC_grad_ineq + 42;
eihMPC_FLOAT eihMPC_ctv3[14];
eihMPC_FLOAT* eihMPC_v3 = eihMPC_v + 21;
eihMPC_FLOAT eihMPC_re3[7];
eihMPC_FLOAT eihMPC_beta3[7];
eihMPC_FLOAT eihMPC_betacc3[7];
eihMPC_FLOAT* eihMPC_dvaff3 = eihMPC_dv_aff + 21;
eihMPC_FLOAT* eihMPC_dvcc3 = eihMPC_dv_cc + 21;
eihMPC_FLOAT eihMPC_V3[98];
eihMPC_FLOAT eihMPC_Yd3[28];
eihMPC_FLOAT eihMPC_Ld3[28];
eihMPC_FLOAT eihMPC_yy3[7];
eihMPC_FLOAT eihMPC_bmy3[7];
eihMPC_FLOAT eihMPC_c3[7] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
int eihMPC_lbIdx3[14] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13};
eihMPC_FLOAT* eihMPC_llb3 = eihMPC_l + 84;
eihMPC_FLOAT* eihMPC_slb3 = eihMPC_s + 84;
eihMPC_FLOAT* eihMPC_llbbyslb3 = eihMPC_lbys + 84;
eihMPC_FLOAT eihMPC_rilb3[14];
eihMPC_FLOAT* eihMPC_dllbaff3 = eihMPC_dl_aff + 84;
eihMPC_FLOAT* eihMPC_dslbaff3 = eihMPC_ds_aff + 84;
eihMPC_FLOAT* eihMPC_dllbcc3 = eihMPC_dl_cc + 84;
eihMPC_FLOAT* eihMPC_dslbcc3 = eihMPC_ds_cc + 84;
eihMPC_FLOAT* eihMPC_ccrhsl3 = eihMPC_ccrhs + 84;
int eihMPC_ubIdx3[14] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13};
eihMPC_FLOAT* eihMPC_lub3 = eihMPC_l + 98;
eihMPC_FLOAT* eihMPC_sub3 = eihMPC_s + 98;
eihMPC_FLOAT* eihMPC_lubbysub3 = eihMPC_lbys + 98;
eihMPC_FLOAT eihMPC_riub3[14];
eihMPC_FLOAT* eihMPC_dlubaff3 = eihMPC_dl_aff + 98;
eihMPC_FLOAT* eihMPC_dsubaff3 = eihMPC_ds_aff + 98;
eihMPC_FLOAT* eihMPC_dlubcc3 = eihMPC_dl_cc + 98;
eihMPC_FLOAT* eihMPC_dsubcc3 = eihMPC_ds_cc + 98;
eihMPC_FLOAT* eihMPC_ccrhsub3 = eihMPC_ccrhs + 98;
eihMPC_FLOAT eihMPC_Phi3[14];
eihMPC_FLOAT eihMPC_W3[14];
eihMPC_FLOAT eihMPC_Ysd3[49];
eihMPC_FLOAT eihMPC_Lsd3[49];
eihMPC_FLOAT* eihMPC_z4 = eihMPC_z + 56;
eihMPC_FLOAT* eihMPC_dzaff4 = eihMPC_dz_aff + 56;
eihMPC_FLOAT* eihMPC_dzcc4 = eihMPC_dz_cc + 56;
eihMPC_FLOAT* eihMPC_rd4 = eihMPC_rd + 56;
eihMPC_FLOAT eihMPC_Lbyrd4[14];
eihMPC_FLOAT* eihMPC_grad_cost4 = eihMPC_grad_cost + 56;
eihMPC_FLOAT* eihMPC_grad_eq4 = eihMPC_grad_eq + 56;
eihMPC_FLOAT* eihMPC_grad_ineq4 = eihMPC_grad_ineq + 56;
eihMPC_FLOAT eihMPC_ctv4[14];
eihMPC_FLOAT* eihMPC_v4 = eihMPC_v + 28;
eihMPC_FLOAT eihMPC_re4[7];
eihMPC_FLOAT eihMPC_beta4[7];
eihMPC_FLOAT eihMPC_betacc4[7];
eihMPC_FLOAT* eihMPC_dvaff4 = eihMPC_dv_aff + 28;
eihMPC_FLOAT* eihMPC_dvcc4 = eihMPC_dv_cc + 28;
eihMPC_FLOAT eihMPC_V4[98];
eihMPC_FLOAT eihMPC_Yd4[28];
eihMPC_FLOAT eihMPC_Ld4[28];
eihMPC_FLOAT eihMPC_yy4[7];
eihMPC_FLOAT eihMPC_bmy4[7];
eihMPC_FLOAT eihMPC_c4[7] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
int eihMPC_lbIdx4[14] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13};
eihMPC_FLOAT* eihMPC_llb4 = eihMPC_l + 112;
eihMPC_FLOAT* eihMPC_slb4 = eihMPC_s + 112;
eihMPC_FLOAT* eihMPC_llbbyslb4 = eihMPC_lbys + 112;
eihMPC_FLOAT eihMPC_rilb4[14];
eihMPC_FLOAT* eihMPC_dllbaff4 = eihMPC_dl_aff + 112;
eihMPC_FLOAT* eihMPC_dslbaff4 = eihMPC_ds_aff + 112;
eihMPC_FLOAT* eihMPC_dllbcc4 = eihMPC_dl_cc + 112;
eihMPC_FLOAT* eihMPC_dslbcc4 = eihMPC_ds_cc + 112;
eihMPC_FLOAT* eihMPC_ccrhsl4 = eihMPC_ccrhs + 112;
int eihMPC_ubIdx4[14] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13};
eihMPC_FLOAT* eihMPC_lub4 = eihMPC_l + 126;
eihMPC_FLOAT* eihMPC_sub4 = eihMPC_s + 126;
eihMPC_FLOAT* eihMPC_lubbysub4 = eihMPC_lbys + 126;
eihMPC_FLOAT eihMPC_riub4[14];
eihMPC_FLOAT* eihMPC_dlubaff4 = eihMPC_dl_aff + 126;
eihMPC_FLOAT* eihMPC_dsubaff4 = eihMPC_ds_aff + 126;
eihMPC_FLOAT* eihMPC_dlubcc4 = eihMPC_dl_cc + 126;
eihMPC_FLOAT* eihMPC_dsubcc4 = eihMPC_ds_cc + 126;
eihMPC_FLOAT* eihMPC_ccrhsub4 = eihMPC_ccrhs + 126;
eihMPC_FLOAT eihMPC_Phi4[14];
eihMPC_FLOAT eihMPC_W4[14];
eihMPC_FLOAT eihMPC_Ysd4[49];
eihMPC_FLOAT eihMPC_Lsd4[49];
eihMPC_FLOAT* eihMPC_z5 = eihMPC_z + 70;
eihMPC_FLOAT* eihMPC_dzaff5 = eihMPC_dz_aff + 70;
eihMPC_FLOAT* eihMPC_dzcc5 = eihMPC_dz_cc + 70;
eihMPC_FLOAT* eihMPC_rd5 = eihMPC_rd + 70;
eihMPC_FLOAT eihMPC_Lbyrd5[14];
eihMPC_FLOAT* eihMPC_grad_cost5 = eihMPC_grad_cost + 70;
eihMPC_FLOAT* eihMPC_grad_eq5 = eihMPC_grad_eq + 70;
eihMPC_FLOAT* eihMPC_grad_ineq5 = eihMPC_grad_ineq + 70;
eihMPC_FLOAT eihMPC_ctv5[14];
eihMPC_FLOAT* eihMPC_v5 = eihMPC_v + 35;
eihMPC_FLOAT eihMPC_re5[7];
eihMPC_FLOAT eihMPC_beta5[7];
eihMPC_FLOAT eihMPC_betacc5[7];
eihMPC_FLOAT* eihMPC_dvaff5 = eihMPC_dv_aff + 35;
eihMPC_FLOAT* eihMPC_dvcc5 = eihMPC_dv_cc + 35;
eihMPC_FLOAT eihMPC_V5[98];
eihMPC_FLOAT eihMPC_Yd5[28];
eihMPC_FLOAT eihMPC_Ld5[28];
eihMPC_FLOAT eihMPC_yy5[7];
eihMPC_FLOAT eihMPC_bmy5[7];
eihMPC_FLOAT eihMPC_c5[7] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
int eihMPC_lbIdx5[14] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13};
eihMPC_FLOAT* eihMPC_llb5 = eihMPC_l + 140;
eihMPC_FLOAT* eihMPC_slb5 = eihMPC_s + 140;
eihMPC_FLOAT* eihMPC_llbbyslb5 = eihMPC_lbys + 140;
eihMPC_FLOAT eihMPC_rilb5[14];
eihMPC_FLOAT* eihMPC_dllbaff5 = eihMPC_dl_aff + 140;
eihMPC_FLOAT* eihMPC_dslbaff5 = eihMPC_ds_aff + 140;
eihMPC_FLOAT* eihMPC_dllbcc5 = eihMPC_dl_cc + 140;
eihMPC_FLOAT* eihMPC_dslbcc5 = eihMPC_ds_cc + 140;
eihMPC_FLOAT* eihMPC_ccrhsl5 = eihMPC_ccrhs + 140;
int eihMPC_ubIdx5[14] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13};
eihMPC_FLOAT* eihMPC_lub5 = eihMPC_l + 154;
eihMPC_FLOAT* eihMPC_sub5 = eihMPC_s + 154;
eihMPC_FLOAT* eihMPC_lubbysub5 = eihMPC_lbys + 154;
eihMPC_FLOAT eihMPC_riub5[14];
eihMPC_FLOAT* eihMPC_dlubaff5 = eihMPC_dl_aff + 154;
eihMPC_FLOAT* eihMPC_dsubaff5 = eihMPC_ds_aff + 154;
eihMPC_FLOAT* eihMPC_dlubcc5 = eihMPC_dl_cc + 154;
eihMPC_FLOAT* eihMPC_dsubcc5 = eihMPC_ds_cc + 154;
eihMPC_FLOAT* eihMPC_ccrhsub5 = eihMPC_ccrhs + 154;
eihMPC_FLOAT eihMPC_Phi5[14];
eihMPC_FLOAT eihMPC_W5[14];
eihMPC_FLOAT eihMPC_Ysd5[49];
eihMPC_FLOAT eihMPC_Lsd5[49];
eihMPC_FLOAT* eihMPC_z6 = eihMPC_z + 84;
eihMPC_FLOAT* eihMPC_dzaff6 = eihMPC_dz_aff + 84;
eihMPC_FLOAT* eihMPC_dzcc6 = eihMPC_dz_cc + 84;
eihMPC_FLOAT* eihMPC_rd6 = eihMPC_rd + 84;
eihMPC_FLOAT eihMPC_Lbyrd6[14];
eihMPC_FLOAT* eihMPC_grad_cost6 = eihMPC_grad_cost + 84;
eihMPC_FLOAT* eihMPC_grad_eq6 = eihMPC_grad_eq + 84;
eihMPC_FLOAT* eihMPC_grad_ineq6 = eihMPC_grad_ineq + 84;
eihMPC_FLOAT eihMPC_ctv6[14];
eihMPC_FLOAT* eihMPC_v6 = eihMPC_v + 42;
eihMPC_FLOAT eihMPC_re6[7];
eihMPC_FLOAT eihMPC_beta6[7];
eihMPC_FLOAT eihMPC_betacc6[7];
eihMPC_FLOAT* eihMPC_dvaff6 = eihMPC_dv_aff + 42;
eihMPC_FLOAT* eihMPC_dvcc6 = eihMPC_dv_cc + 42;
eihMPC_FLOAT eihMPC_V6[98];
eihMPC_FLOAT eihMPC_Yd6[28];
eihMPC_FLOAT eihMPC_Ld6[28];
eihMPC_FLOAT eihMPC_yy6[7];
eihMPC_FLOAT eihMPC_bmy6[7];
eihMPC_FLOAT eihMPC_c6[7] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
int eihMPC_lbIdx6[14] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13};
eihMPC_FLOAT* eihMPC_llb6 = eihMPC_l + 168;
eihMPC_FLOAT* eihMPC_slb6 = eihMPC_s + 168;
eihMPC_FLOAT* eihMPC_llbbyslb6 = eihMPC_lbys + 168;
eihMPC_FLOAT eihMPC_rilb6[14];
eihMPC_FLOAT* eihMPC_dllbaff6 = eihMPC_dl_aff + 168;
eihMPC_FLOAT* eihMPC_dslbaff6 = eihMPC_ds_aff + 168;
eihMPC_FLOAT* eihMPC_dllbcc6 = eihMPC_dl_cc + 168;
eihMPC_FLOAT* eihMPC_dslbcc6 = eihMPC_ds_cc + 168;
eihMPC_FLOAT* eihMPC_ccrhsl6 = eihMPC_ccrhs + 168;
int eihMPC_ubIdx6[14] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13};
eihMPC_FLOAT* eihMPC_lub6 = eihMPC_l + 182;
eihMPC_FLOAT* eihMPC_sub6 = eihMPC_s + 182;
eihMPC_FLOAT* eihMPC_lubbysub6 = eihMPC_lbys + 182;
eihMPC_FLOAT eihMPC_riub6[14];
eihMPC_FLOAT* eihMPC_dlubaff6 = eihMPC_dl_aff + 182;
eihMPC_FLOAT* eihMPC_dsubaff6 = eihMPC_ds_aff + 182;
eihMPC_FLOAT* eihMPC_dlubcc6 = eihMPC_dl_cc + 182;
eihMPC_FLOAT* eihMPC_dsubcc6 = eihMPC_ds_cc + 182;
eihMPC_FLOAT* eihMPC_ccrhsub6 = eihMPC_ccrhs + 182;
eihMPC_FLOAT eihMPC_Phi6[14];
eihMPC_FLOAT eihMPC_W6[14];
eihMPC_FLOAT eihMPC_Ysd6[49];
eihMPC_FLOAT eihMPC_Lsd6[49];
eihMPC_FLOAT* eihMPC_z7 = eihMPC_z + 98;
eihMPC_FLOAT* eihMPC_dzaff7 = eihMPC_dz_aff + 98;
eihMPC_FLOAT* eihMPC_dzcc7 = eihMPC_dz_cc + 98;
eihMPC_FLOAT* eihMPC_rd7 = eihMPC_rd + 98;
eihMPC_FLOAT eihMPC_Lbyrd7[14];
eihMPC_FLOAT* eihMPC_grad_cost7 = eihMPC_grad_cost + 98;
eihMPC_FLOAT* eihMPC_grad_eq7 = eihMPC_grad_eq + 98;
eihMPC_FLOAT* eihMPC_grad_ineq7 = eihMPC_grad_ineq + 98;
eihMPC_FLOAT eihMPC_ctv7[14];
eihMPC_FLOAT* eihMPC_v7 = eihMPC_v + 49;
eihMPC_FLOAT eihMPC_re7[7];
eihMPC_FLOAT eihMPC_beta7[7];
eihMPC_FLOAT eihMPC_betacc7[7];
eihMPC_FLOAT* eihMPC_dvaff7 = eihMPC_dv_aff + 49;
eihMPC_FLOAT* eihMPC_dvcc7 = eihMPC_dv_cc + 49;
eihMPC_FLOAT eihMPC_V7[98];
eihMPC_FLOAT eihMPC_Yd7[28];
eihMPC_FLOAT eihMPC_Ld7[28];
eihMPC_FLOAT eihMPC_yy7[7];
eihMPC_FLOAT eihMPC_bmy7[7];
eihMPC_FLOAT eihMPC_c7[7] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
int eihMPC_lbIdx7[14] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13};
eihMPC_FLOAT* eihMPC_llb7 = eihMPC_l + 196;
eihMPC_FLOAT* eihMPC_slb7 = eihMPC_s + 196;
eihMPC_FLOAT* eihMPC_llbbyslb7 = eihMPC_lbys + 196;
eihMPC_FLOAT eihMPC_rilb7[14];
eihMPC_FLOAT* eihMPC_dllbaff7 = eihMPC_dl_aff + 196;
eihMPC_FLOAT* eihMPC_dslbaff7 = eihMPC_ds_aff + 196;
eihMPC_FLOAT* eihMPC_dllbcc7 = eihMPC_dl_cc + 196;
eihMPC_FLOAT* eihMPC_dslbcc7 = eihMPC_ds_cc + 196;
eihMPC_FLOAT* eihMPC_ccrhsl7 = eihMPC_ccrhs + 196;
int eihMPC_ubIdx7[14] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13};
eihMPC_FLOAT* eihMPC_lub7 = eihMPC_l + 210;
eihMPC_FLOAT* eihMPC_sub7 = eihMPC_s + 210;
eihMPC_FLOAT* eihMPC_lubbysub7 = eihMPC_lbys + 210;
eihMPC_FLOAT eihMPC_riub7[14];
eihMPC_FLOAT* eihMPC_dlubaff7 = eihMPC_dl_aff + 210;
eihMPC_FLOAT* eihMPC_dsubaff7 = eihMPC_ds_aff + 210;
eihMPC_FLOAT* eihMPC_dlubcc7 = eihMPC_dl_cc + 210;
eihMPC_FLOAT* eihMPC_dsubcc7 = eihMPC_ds_cc + 210;
eihMPC_FLOAT* eihMPC_ccrhsub7 = eihMPC_ccrhs + 210;
eihMPC_FLOAT eihMPC_Phi7[14];
eihMPC_FLOAT eihMPC_W7[14];
eihMPC_FLOAT eihMPC_Ysd7[49];
eihMPC_FLOAT eihMPC_Lsd7[49];
eihMPC_FLOAT* eihMPC_z8 = eihMPC_z + 112;
eihMPC_FLOAT* eihMPC_dzaff8 = eihMPC_dz_aff + 112;
eihMPC_FLOAT* eihMPC_dzcc8 = eihMPC_dz_cc + 112;
eihMPC_FLOAT* eihMPC_rd8 = eihMPC_rd + 112;
eihMPC_FLOAT eihMPC_Lbyrd8[14];
eihMPC_FLOAT* eihMPC_grad_cost8 = eihMPC_grad_cost + 112;
eihMPC_FLOAT* eihMPC_grad_eq8 = eihMPC_grad_eq + 112;
eihMPC_FLOAT* eihMPC_grad_ineq8 = eihMPC_grad_ineq + 112;
eihMPC_FLOAT eihMPC_ctv8[14];
eihMPC_FLOAT* eihMPC_v8 = eihMPC_v + 56;
eihMPC_FLOAT eihMPC_re8[7];
eihMPC_FLOAT eihMPC_beta8[7];
eihMPC_FLOAT eihMPC_betacc8[7];
eihMPC_FLOAT* eihMPC_dvaff8 = eihMPC_dv_aff + 56;
eihMPC_FLOAT* eihMPC_dvcc8 = eihMPC_dv_cc + 56;
eihMPC_FLOAT eihMPC_V8[98];
eihMPC_FLOAT eihMPC_Yd8[28];
eihMPC_FLOAT eihMPC_Ld8[28];
eihMPC_FLOAT eihMPC_yy8[7];
eihMPC_FLOAT eihMPC_bmy8[7];
eihMPC_FLOAT eihMPC_c8[7] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
int eihMPC_lbIdx8[14] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13};
eihMPC_FLOAT* eihMPC_llb8 = eihMPC_l + 224;
eihMPC_FLOAT* eihMPC_slb8 = eihMPC_s + 224;
eihMPC_FLOAT* eihMPC_llbbyslb8 = eihMPC_lbys + 224;
eihMPC_FLOAT eihMPC_rilb8[14];
eihMPC_FLOAT* eihMPC_dllbaff8 = eihMPC_dl_aff + 224;
eihMPC_FLOAT* eihMPC_dslbaff8 = eihMPC_ds_aff + 224;
eihMPC_FLOAT* eihMPC_dllbcc8 = eihMPC_dl_cc + 224;
eihMPC_FLOAT* eihMPC_dslbcc8 = eihMPC_ds_cc + 224;
eihMPC_FLOAT* eihMPC_ccrhsl8 = eihMPC_ccrhs + 224;
int eihMPC_ubIdx8[14] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13};
eihMPC_FLOAT* eihMPC_lub8 = eihMPC_l + 238;
eihMPC_FLOAT* eihMPC_sub8 = eihMPC_s + 238;
eihMPC_FLOAT* eihMPC_lubbysub8 = eihMPC_lbys + 238;
eihMPC_FLOAT eihMPC_riub8[14];
eihMPC_FLOAT* eihMPC_dlubaff8 = eihMPC_dl_aff + 238;
eihMPC_FLOAT* eihMPC_dsubaff8 = eihMPC_ds_aff + 238;
eihMPC_FLOAT* eihMPC_dlubcc8 = eihMPC_dl_cc + 238;
eihMPC_FLOAT* eihMPC_dsubcc8 = eihMPC_ds_cc + 238;
eihMPC_FLOAT* eihMPC_ccrhsub8 = eihMPC_ccrhs + 238;
eihMPC_FLOAT eihMPC_Phi8[14];
eihMPC_FLOAT eihMPC_W8[14];
eihMPC_FLOAT eihMPC_Ysd8[49];
eihMPC_FLOAT eihMPC_Lsd8[49];
eihMPC_FLOAT* eihMPC_z9 = eihMPC_z + 126;
eihMPC_FLOAT* eihMPC_dzaff9 = eihMPC_dz_aff + 126;
eihMPC_FLOAT* eihMPC_dzcc9 = eihMPC_dz_cc + 126;
eihMPC_FLOAT* eihMPC_rd9 = eihMPC_rd + 126;
eihMPC_FLOAT eihMPC_Lbyrd9[7];
eihMPC_FLOAT* eihMPC_grad_cost9 = eihMPC_grad_cost + 126;
eihMPC_FLOAT* eihMPC_grad_eq9 = eihMPC_grad_eq + 126;
eihMPC_FLOAT* eihMPC_grad_ineq9 = eihMPC_grad_ineq + 126;
eihMPC_FLOAT eihMPC_ctv9[7];
eihMPC_FLOAT* eihMPC_v9 = eihMPC_v + 63;
eihMPC_FLOAT eihMPC_re9[7];
eihMPC_FLOAT eihMPC_beta9[7];
eihMPC_FLOAT eihMPC_betacc9[7];
eihMPC_FLOAT* eihMPC_dvaff9 = eihMPC_dv_aff + 63;
eihMPC_FLOAT* eihMPC_dvcc9 = eihMPC_dv_cc + 63;
eihMPC_FLOAT eihMPC_V9[49];
eihMPC_FLOAT eihMPC_Yd9[28];
eihMPC_FLOAT eihMPC_Ld9[28];
eihMPC_FLOAT eihMPC_yy9[7];
eihMPC_FLOAT eihMPC_bmy9[7];
eihMPC_FLOAT eihMPC_c9[7] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
int eihMPC_lbIdx9[7] = {0, 1, 2, 3, 4, 5, 6};
eihMPC_FLOAT* eihMPC_llb9 = eihMPC_l + 252;
eihMPC_FLOAT* eihMPC_slb9 = eihMPC_s + 252;
eihMPC_FLOAT* eihMPC_llbbyslb9 = eihMPC_lbys + 252;
eihMPC_FLOAT eihMPC_rilb9[7];
eihMPC_FLOAT* eihMPC_dllbaff9 = eihMPC_dl_aff + 252;
eihMPC_FLOAT* eihMPC_dslbaff9 = eihMPC_ds_aff + 252;
eihMPC_FLOAT* eihMPC_dllbcc9 = eihMPC_dl_cc + 252;
eihMPC_FLOAT* eihMPC_dslbcc9 = eihMPC_ds_cc + 252;
eihMPC_FLOAT* eihMPC_ccrhsl9 = eihMPC_ccrhs + 252;
int eihMPC_ubIdx9[7] = {0, 1, 2, 3, 4, 5, 6};
eihMPC_FLOAT* eihMPC_lub9 = eihMPC_l + 259;
eihMPC_FLOAT* eihMPC_sub9 = eihMPC_s + 259;
eihMPC_FLOAT* eihMPC_lubbysub9 = eihMPC_lbys + 259;
eihMPC_FLOAT eihMPC_riub9[7];
eihMPC_FLOAT* eihMPC_dlubaff9 = eihMPC_dl_aff + 259;
eihMPC_FLOAT* eihMPC_dsubaff9 = eihMPC_ds_aff + 259;
eihMPC_FLOAT* eihMPC_dlubcc9 = eihMPC_dl_cc + 259;
eihMPC_FLOAT* eihMPC_dsubcc9 = eihMPC_ds_cc + 259;
eihMPC_FLOAT* eihMPC_ccrhsub9 = eihMPC_ccrhs + 259;
eihMPC_FLOAT eihMPC_Phi9[7];
eihMPC_FLOAT eihMPC_D9[7] = {-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000};
eihMPC_FLOAT eihMPC_W9[7];
eihMPC_FLOAT eihMPC_Ysd9[49];
eihMPC_FLOAT eihMPC_Lsd9[49];
eihMPC_FLOAT musigma;
eihMPC_FLOAT sigma_3rdroot;
eihMPC_FLOAT eihMPC_Diag1_0[14];
eihMPC_FLOAT eihMPC_Diag2_0[14];
eihMPC_FLOAT eihMPC_L_0[91];




/* SOLVER CODE --------------------------------------------------------- */
int eihMPC_solve(eihMPC_params* params, eihMPC_output* output, eihMPC_info* info)
{	
int exitcode;

#if eihMPC_SET_TIMING == 1
	eihMPC_timer solvertimer;
	eihMPC_tic(&solvertimer);
#endif
/* FUNCTION CALLS INTO LA LIBRARY -------------------------------------- */
info->it = 0;
eihMPC_LA_INITIALIZEVECTOR_133(eihMPC_z, 0);
eihMPC_LA_INITIALIZEVECTOR_70(eihMPC_v, 1);
eihMPC_LA_INITIALIZEVECTOR_266(eihMPC_l, 1);
eihMPC_LA_INITIALIZEVECTOR_266(eihMPC_s, 1);
info->mu = 0;
eihMPC_LA_DOTACC_266(eihMPC_l, eihMPC_s, &info->mu);
info->mu /= 266;
while( 1 ){
info->pobj = 0;
eihMPC_LA_DIAG_QUADFCN_14(params->H1, params->f1, eihMPC_z0, eihMPC_grad_cost0, &info->pobj);
eihMPC_LA_DIAG_QUADFCN_14(params->H2, params->f2, eihMPC_z1, eihMPC_grad_cost1, &info->pobj);
eihMPC_LA_DIAG_QUADFCN_14(params->H3, params->f3, eihMPC_z2, eihMPC_grad_cost2, &info->pobj);
eihMPC_LA_DIAG_QUADFCN_14(params->H4, params->f4, eihMPC_z3, eihMPC_grad_cost3, &info->pobj);
eihMPC_LA_DIAG_QUADFCN_14(params->H5, params->f5, eihMPC_z4, eihMPC_grad_cost4, &info->pobj);
eihMPC_LA_DIAG_QUADFCN_14(params->H6, params->f6, eihMPC_z5, eihMPC_grad_cost5, &info->pobj);
eihMPC_LA_DIAG_QUADFCN_14(params->H7, params->f7, eihMPC_z6, eihMPC_grad_cost6, &info->pobj);
eihMPC_LA_DIAG_QUADFCN_14(params->H8, params->f8, eihMPC_z7, eihMPC_grad_cost7, &info->pobj);
eihMPC_LA_DIAG_QUADFCN_14(params->H9, params->f9, eihMPC_z8, eihMPC_grad_cost8, &info->pobj);
eihMPC_LA_DIAG_QUADFCN_7(params->H10, params->f10, eihMPC_z9, eihMPC_grad_cost9, &info->pobj);
info->res_eq = 0;
info->dgap = 0;
eihMPC_LA_DIAGZERO_MVMSUB6_7(eihMPC_D0, eihMPC_z0, params->c1, eihMPC_v0, eihMPC_re0, &info->dgap, &info->res_eq);
eihMPC_LA_DENSE_DIAGZERO_MVMSUB3_7_14_14(eihMPC_C0, eihMPC_z0, eihMPC_D1, eihMPC_z1, eihMPC_c1, eihMPC_v1, eihMPC_re1, &info->dgap, &info->res_eq);
eihMPC_LA_DENSE_DIAGZERO_MVMSUB3_7_14_14(eihMPC_C0, eihMPC_z1, eihMPC_D1, eihMPC_z2, eihMPC_c2, eihMPC_v2, eihMPC_re2, &info->dgap, &info->res_eq);
eihMPC_LA_DENSE_DIAGZERO_MVMSUB3_7_14_14(eihMPC_C0, eihMPC_z2, eihMPC_D1, eihMPC_z3, eihMPC_c3, eihMPC_v3, eihMPC_re3, &info->dgap, &info->res_eq);
eihMPC_LA_DENSE_DIAGZERO_MVMSUB3_7_14_14(eihMPC_C0, eihMPC_z3, eihMPC_D1, eihMPC_z4, eihMPC_c4, eihMPC_v4, eihMPC_re4, &info->dgap, &info->res_eq);
eihMPC_LA_DENSE_DIAGZERO_MVMSUB3_7_14_14(eihMPC_C0, eihMPC_z4, eihMPC_D1, eihMPC_z5, eihMPC_c5, eihMPC_v5, eihMPC_re5, &info->dgap, &info->res_eq);
eihMPC_LA_DENSE_DIAGZERO_MVMSUB3_7_14_14(eihMPC_C0, eihMPC_z5, eihMPC_D1, eihMPC_z6, eihMPC_c6, eihMPC_v6, eihMPC_re6, &info->dgap, &info->res_eq);
eihMPC_LA_DENSE_DIAGZERO_MVMSUB3_7_14_14(eihMPC_C0, eihMPC_z6, eihMPC_D1, eihMPC_z7, eihMPC_c7, eihMPC_v7, eihMPC_re7, &info->dgap, &info->res_eq);
eihMPC_LA_DENSE_DIAGZERO_MVMSUB3_7_14_14(eihMPC_C0, eihMPC_z7, eihMPC_D1, eihMPC_z8, eihMPC_c8, eihMPC_v8, eihMPC_re8, &info->dgap, &info->res_eq);
eihMPC_LA_DENSE_DIAGZERO_MVMSUB3_7_14_7(eihMPC_C0, eihMPC_z8, eihMPC_D9, eihMPC_z9, eihMPC_c9, eihMPC_v9, eihMPC_re9, &info->dgap, &info->res_eq);
eihMPC_LA_DENSE_DIAGZERO_MTVM2_7_14_7(eihMPC_C0, eihMPC_v1, eihMPC_D0, eihMPC_v0, eihMPC_grad_eq0);
eihMPC_LA_DENSE_DIAGZERO_MTVM2_7_14_7(eihMPC_C0, eihMPC_v2, eihMPC_D1, eihMPC_v1, eihMPC_grad_eq1);
eihMPC_LA_DENSE_DIAGZERO_MTVM2_7_14_7(eihMPC_C0, eihMPC_v3, eihMPC_D1, eihMPC_v2, eihMPC_grad_eq2);
eihMPC_LA_DENSE_DIAGZERO_MTVM2_7_14_7(eihMPC_C0, eihMPC_v4, eihMPC_D1, eihMPC_v3, eihMPC_grad_eq3);
eihMPC_LA_DENSE_DIAGZERO_MTVM2_7_14_7(eihMPC_C0, eihMPC_v5, eihMPC_D1, eihMPC_v4, eihMPC_grad_eq4);
eihMPC_LA_DENSE_DIAGZERO_MTVM2_7_14_7(eihMPC_C0, eihMPC_v6, eihMPC_D1, eihMPC_v5, eihMPC_grad_eq5);
eihMPC_LA_DENSE_DIAGZERO_MTVM2_7_14_7(eihMPC_C0, eihMPC_v7, eihMPC_D1, eihMPC_v6, eihMPC_grad_eq6);
eihMPC_LA_DENSE_DIAGZERO_MTVM2_7_14_7(eihMPC_C0, eihMPC_v8, eihMPC_D1, eihMPC_v7, eihMPC_grad_eq7);
eihMPC_LA_DENSE_DIAGZERO_MTVM2_7_14_7(eihMPC_C0, eihMPC_v9, eihMPC_D1, eihMPC_v8, eihMPC_grad_eq8);
eihMPC_LA_DIAGZERO_MTVM_7_7(eihMPC_D9, eihMPC_v9, eihMPC_grad_eq9);
info->res_ineq = 0;
eihMPC_LA_VSUBADD3_14(params->lb1, eihMPC_z0, eihMPC_lbIdx0, eihMPC_llb0, eihMPC_slb0, eihMPC_rilb0, &info->dgap, &info->res_ineq);
eihMPC_LA_VSUBADD2_14(eihMPC_z0, eihMPC_ubIdx0, params->ub1, eihMPC_lub0, eihMPC_sub0, eihMPC_riub0, &info->dgap, &info->res_ineq);
eihMPC_LA_VSUBADD3_14(params->lb2, eihMPC_z1, eihMPC_lbIdx1, eihMPC_llb1, eihMPC_slb1, eihMPC_rilb1, &info->dgap, &info->res_ineq);
eihMPC_LA_VSUBADD2_14(eihMPC_z1, eihMPC_ubIdx1, params->ub2, eihMPC_lub1, eihMPC_sub1, eihMPC_riub1, &info->dgap, &info->res_ineq);
eihMPC_LA_VSUBADD3_14(params->lb3, eihMPC_z2, eihMPC_lbIdx2, eihMPC_llb2, eihMPC_slb2, eihMPC_rilb2, &info->dgap, &info->res_ineq);
eihMPC_LA_VSUBADD2_14(eihMPC_z2, eihMPC_ubIdx2, params->ub3, eihMPC_lub2, eihMPC_sub2, eihMPC_riub2, &info->dgap, &info->res_ineq);
eihMPC_LA_VSUBADD3_14(params->lb4, eihMPC_z3, eihMPC_lbIdx3, eihMPC_llb3, eihMPC_slb3, eihMPC_rilb3, &info->dgap, &info->res_ineq);
eihMPC_LA_VSUBADD2_14(eihMPC_z3, eihMPC_ubIdx3, params->ub4, eihMPC_lub3, eihMPC_sub3, eihMPC_riub3, &info->dgap, &info->res_ineq);
eihMPC_LA_VSUBADD3_14(params->lb5, eihMPC_z4, eihMPC_lbIdx4, eihMPC_llb4, eihMPC_slb4, eihMPC_rilb4, &info->dgap, &info->res_ineq);
eihMPC_LA_VSUBADD2_14(eihMPC_z4, eihMPC_ubIdx4, params->ub5, eihMPC_lub4, eihMPC_sub4, eihMPC_riub4, &info->dgap, &info->res_ineq);
eihMPC_LA_VSUBADD3_14(params->lb6, eihMPC_z5, eihMPC_lbIdx5, eihMPC_llb5, eihMPC_slb5, eihMPC_rilb5, &info->dgap, &info->res_ineq);
eihMPC_LA_VSUBADD2_14(eihMPC_z5, eihMPC_ubIdx5, params->ub6, eihMPC_lub5, eihMPC_sub5, eihMPC_riub5, &info->dgap, &info->res_ineq);
eihMPC_LA_VSUBADD3_14(params->lb7, eihMPC_z6, eihMPC_lbIdx6, eihMPC_llb6, eihMPC_slb6, eihMPC_rilb6, &info->dgap, &info->res_ineq);
eihMPC_LA_VSUBADD2_14(eihMPC_z6, eihMPC_ubIdx6, params->ub7, eihMPC_lub6, eihMPC_sub6, eihMPC_riub6, &info->dgap, &info->res_ineq);
eihMPC_LA_VSUBADD3_14(params->lb8, eihMPC_z7, eihMPC_lbIdx7, eihMPC_llb7, eihMPC_slb7, eihMPC_rilb7, &info->dgap, &info->res_ineq);
eihMPC_LA_VSUBADD2_14(eihMPC_z7, eihMPC_ubIdx7, params->ub8, eihMPC_lub7, eihMPC_sub7, eihMPC_riub7, &info->dgap, &info->res_ineq);
eihMPC_LA_VSUBADD3_14(params->lb9, eihMPC_z8, eihMPC_lbIdx8, eihMPC_llb8, eihMPC_slb8, eihMPC_rilb8, &info->dgap, &info->res_ineq);
eihMPC_LA_VSUBADD2_14(eihMPC_z8, eihMPC_ubIdx8, params->ub9, eihMPC_lub8, eihMPC_sub8, eihMPC_riub8, &info->dgap, &info->res_ineq);
eihMPC_LA_VSUBADD3_7(params->lb10, eihMPC_z9, eihMPC_lbIdx9, eihMPC_llb9, eihMPC_slb9, eihMPC_rilb9, &info->dgap, &info->res_ineq);
eihMPC_LA_VSUBADD2_7(eihMPC_z9, eihMPC_ubIdx9, params->ub10, eihMPC_lub9, eihMPC_sub9, eihMPC_riub9, &info->dgap, &info->res_ineq);
eihMPC_LA_INEQ_B_GRAD_14_14_14(eihMPC_lub0, eihMPC_sub0, eihMPC_riub0, eihMPC_llb0, eihMPC_slb0, eihMPC_rilb0, eihMPC_lbIdx0, eihMPC_ubIdx0, eihMPC_grad_ineq0, eihMPC_lubbysub0, eihMPC_llbbyslb0);
eihMPC_LA_INEQ_B_GRAD_14_14_14(eihMPC_lub1, eihMPC_sub1, eihMPC_riub1, eihMPC_llb1, eihMPC_slb1, eihMPC_rilb1, eihMPC_lbIdx1, eihMPC_ubIdx1, eihMPC_grad_ineq1, eihMPC_lubbysub1, eihMPC_llbbyslb1);
eihMPC_LA_INEQ_B_GRAD_14_14_14(eihMPC_lub2, eihMPC_sub2, eihMPC_riub2, eihMPC_llb2, eihMPC_slb2, eihMPC_rilb2, eihMPC_lbIdx2, eihMPC_ubIdx2, eihMPC_grad_ineq2, eihMPC_lubbysub2, eihMPC_llbbyslb2);
eihMPC_LA_INEQ_B_GRAD_14_14_14(eihMPC_lub3, eihMPC_sub3, eihMPC_riub3, eihMPC_llb3, eihMPC_slb3, eihMPC_rilb3, eihMPC_lbIdx3, eihMPC_ubIdx3, eihMPC_grad_ineq3, eihMPC_lubbysub3, eihMPC_llbbyslb3);
eihMPC_LA_INEQ_B_GRAD_14_14_14(eihMPC_lub4, eihMPC_sub4, eihMPC_riub4, eihMPC_llb4, eihMPC_slb4, eihMPC_rilb4, eihMPC_lbIdx4, eihMPC_ubIdx4, eihMPC_grad_ineq4, eihMPC_lubbysub4, eihMPC_llbbyslb4);
eihMPC_LA_INEQ_B_GRAD_14_14_14(eihMPC_lub5, eihMPC_sub5, eihMPC_riub5, eihMPC_llb5, eihMPC_slb5, eihMPC_rilb5, eihMPC_lbIdx5, eihMPC_ubIdx5, eihMPC_grad_ineq5, eihMPC_lubbysub5, eihMPC_llbbyslb5);
eihMPC_LA_INEQ_B_GRAD_14_14_14(eihMPC_lub6, eihMPC_sub6, eihMPC_riub6, eihMPC_llb6, eihMPC_slb6, eihMPC_rilb6, eihMPC_lbIdx6, eihMPC_ubIdx6, eihMPC_grad_ineq6, eihMPC_lubbysub6, eihMPC_llbbyslb6);
eihMPC_LA_INEQ_B_GRAD_14_14_14(eihMPC_lub7, eihMPC_sub7, eihMPC_riub7, eihMPC_llb7, eihMPC_slb7, eihMPC_rilb7, eihMPC_lbIdx7, eihMPC_ubIdx7, eihMPC_grad_ineq7, eihMPC_lubbysub7, eihMPC_llbbyslb7);
eihMPC_LA_INEQ_B_GRAD_14_14_14(eihMPC_lub8, eihMPC_sub8, eihMPC_riub8, eihMPC_llb8, eihMPC_slb8, eihMPC_rilb8, eihMPC_lbIdx8, eihMPC_ubIdx8, eihMPC_grad_ineq8, eihMPC_lubbysub8, eihMPC_llbbyslb8);
eihMPC_LA_INEQ_B_GRAD_7_7_7(eihMPC_lub9, eihMPC_sub9, eihMPC_riub9, eihMPC_llb9, eihMPC_slb9, eihMPC_rilb9, eihMPC_lbIdx9, eihMPC_ubIdx9, eihMPC_grad_ineq9, eihMPC_lubbysub9, eihMPC_llbbyslb9);
info->dobj = info->pobj - info->dgap;
info->rdgap = info->pobj ? info->dgap / info->pobj : 1e6;
if( info->rdgap < 0 ) info->rdgap = -info->rdgap;
if( info->mu < eihMPC_SET_ACC_KKTCOMPL
    && (info->rdgap < eihMPC_SET_ACC_RDGAP || info->dgap < eihMPC_SET_ACC_KKTCOMPL)
    && info->res_eq < eihMPC_SET_ACC_RESEQ
    && info->res_ineq < eihMPC_SET_ACC_RESINEQ ){
exitcode = eihMPC_OPTIMAL; break; }
if( info->it == eihMPC_SET_MAXIT ){
exitcode = eihMPC_MAXITREACHED; break; }
eihMPC_LA_VVADD3_133(eihMPC_grad_cost, eihMPC_grad_eq, eihMPC_grad_ineq, eihMPC_rd);
eihMPC_LA_DIAG_CHOL_ONELOOP_LBUB_14_14_14(params->H1, eihMPC_llbbyslb0, eihMPC_lbIdx0, eihMPC_lubbysub0, eihMPC_ubIdx0, eihMPC_Phi0);
eihMPC_LA_DIAG_MATRIXFORWARDSUB_7_14(eihMPC_Phi0, eihMPC_C0, eihMPC_V0);
eihMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_7_14(eihMPC_Phi0, eihMPC_D0, eihMPC_W0);
eihMPC_LA_DENSE_DIAGZERO_MMTM_7_14_7(eihMPC_W0, eihMPC_V0, eihMPC_Ysd1);
eihMPC_LA_DIAG_FORWARDSUB_14(eihMPC_Phi0, eihMPC_rd0, eihMPC_Lbyrd0);
eihMPC_LA_DIAG_CHOL_ONELOOP_LBUB_14_14_14(params->H2, eihMPC_llbbyslb1, eihMPC_lbIdx1, eihMPC_lubbysub1, eihMPC_ubIdx1, eihMPC_Phi1);
eihMPC_LA_DIAG_MATRIXFORWARDSUB_7_14(eihMPC_Phi1, eihMPC_C0, eihMPC_V1);
eihMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_7_14(eihMPC_Phi1, eihMPC_D1, eihMPC_W1);
eihMPC_LA_DENSE_DIAGZERO_MMTM_7_14_7(eihMPC_W1, eihMPC_V1, eihMPC_Ysd2);
eihMPC_LA_DIAG_FORWARDSUB_14(eihMPC_Phi1, eihMPC_rd1, eihMPC_Lbyrd1);
eihMPC_LA_DIAG_CHOL_ONELOOP_LBUB_14_14_14(params->H3, eihMPC_llbbyslb2, eihMPC_lbIdx2, eihMPC_lubbysub2, eihMPC_ubIdx2, eihMPC_Phi2);
eihMPC_LA_DIAG_MATRIXFORWARDSUB_7_14(eihMPC_Phi2, eihMPC_C0, eihMPC_V2);
eihMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_7_14(eihMPC_Phi2, eihMPC_D1, eihMPC_W2);
eihMPC_LA_DENSE_DIAGZERO_MMTM_7_14_7(eihMPC_W2, eihMPC_V2, eihMPC_Ysd3);
eihMPC_LA_DIAG_FORWARDSUB_14(eihMPC_Phi2, eihMPC_rd2, eihMPC_Lbyrd2);
eihMPC_LA_DIAG_CHOL_ONELOOP_LBUB_14_14_14(params->H4, eihMPC_llbbyslb3, eihMPC_lbIdx3, eihMPC_lubbysub3, eihMPC_ubIdx3, eihMPC_Phi3);
eihMPC_LA_DIAG_MATRIXFORWARDSUB_7_14(eihMPC_Phi3, eihMPC_C0, eihMPC_V3);
eihMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_7_14(eihMPC_Phi3, eihMPC_D1, eihMPC_W3);
eihMPC_LA_DENSE_DIAGZERO_MMTM_7_14_7(eihMPC_W3, eihMPC_V3, eihMPC_Ysd4);
eihMPC_LA_DIAG_FORWARDSUB_14(eihMPC_Phi3, eihMPC_rd3, eihMPC_Lbyrd3);
eihMPC_LA_DIAG_CHOL_ONELOOP_LBUB_14_14_14(params->H5, eihMPC_llbbyslb4, eihMPC_lbIdx4, eihMPC_lubbysub4, eihMPC_ubIdx4, eihMPC_Phi4);
eihMPC_LA_DIAG_MATRIXFORWARDSUB_7_14(eihMPC_Phi4, eihMPC_C0, eihMPC_V4);
eihMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_7_14(eihMPC_Phi4, eihMPC_D1, eihMPC_W4);
eihMPC_LA_DENSE_DIAGZERO_MMTM_7_14_7(eihMPC_W4, eihMPC_V4, eihMPC_Ysd5);
eihMPC_LA_DIAG_FORWARDSUB_14(eihMPC_Phi4, eihMPC_rd4, eihMPC_Lbyrd4);
eihMPC_LA_DIAG_CHOL_ONELOOP_LBUB_14_14_14(params->H6, eihMPC_llbbyslb5, eihMPC_lbIdx5, eihMPC_lubbysub5, eihMPC_ubIdx5, eihMPC_Phi5);
eihMPC_LA_DIAG_MATRIXFORWARDSUB_7_14(eihMPC_Phi5, eihMPC_C0, eihMPC_V5);
eihMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_7_14(eihMPC_Phi5, eihMPC_D1, eihMPC_W5);
eihMPC_LA_DENSE_DIAGZERO_MMTM_7_14_7(eihMPC_W5, eihMPC_V5, eihMPC_Ysd6);
eihMPC_LA_DIAG_FORWARDSUB_14(eihMPC_Phi5, eihMPC_rd5, eihMPC_Lbyrd5);
eihMPC_LA_DIAG_CHOL_ONELOOP_LBUB_14_14_14(params->H7, eihMPC_llbbyslb6, eihMPC_lbIdx6, eihMPC_lubbysub6, eihMPC_ubIdx6, eihMPC_Phi6);
eihMPC_LA_DIAG_MATRIXFORWARDSUB_7_14(eihMPC_Phi6, eihMPC_C0, eihMPC_V6);
eihMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_7_14(eihMPC_Phi6, eihMPC_D1, eihMPC_W6);
eihMPC_LA_DENSE_DIAGZERO_MMTM_7_14_7(eihMPC_W6, eihMPC_V6, eihMPC_Ysd7);
eihMPC_LA_DIAG_FORWARDSUB_14(eihMPC_Phi6, eihMPC_rd6, eihMPC_Lbyrd6);
eihMPC_LA_DIAG_CHOL_ONELOOP_LBUB_14_14_14(params->H8, eihMPC_llbbyslb7, eihMPC_lbIdx7, eihMPC_lubbysub7, eihMPC_ubIdx7, eihMPC_Phi7);
eihMPC_LA_DIAG_MATRIXFORWARDSUB_7_14(eihMPC_Phi7, eihMPC_C0, eihMPC_V7);
eihMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_7_14(eihMPC_Phi7, eihMPC_D1, eihMPC_W7);
eihMPC_LA_DENSE_DIAGZERO_MMTM_7_14_7(eihMPC_W7, eihMPC_V7, eihMPC_Ysd8);
eihMPC_LA_DIAG_FORWARDSUB_14(eihMPC_Phi7, eihMPC_rd7, eihMPC_Lbyrd7);
eihMPC_LA_DIAG_CHOL_ONELOOP_LBUB_14_14_14(params->H9, eihMPC_llbbyslb8, eihMPC_lbIdx8, eihMPC_lubbysub8, eihMPC_ubIdx8, eihMPC_Phi8);
eihMPC_LA_DIAG_MATRIXFORWARDSUB_7_14(eihMPC_Phi8, eihMPC_C0, eihMPC_V8);
eihMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_7_14(eihMPC_Phi8, eihMPC_D1, eihMPC_W8);
eihMPC_LA_DENSE_DIAGZERO_MMTM_7_14_7(eihMPC_W8, eihMPC_V8, eihMPC_Ysd9);
eihMPC_LA_DIAG_FORWARDSUB_14(eihMPC_Phi8, eihMPC_rd8, eihMPC_Lbyrd8);
eihMPC_LA_DIAG_CHOL_ONELOOP_LBUB_7_7_7(params->H10, eihMPC_llbbyslb9, eihMPC_lbIdx9, eihMPC_lubbysub9, eihMPC_ubIdx9, eihMPC_Phi9);
eihMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_7_7(eihMPC_Phi9, eihMPC_D9, eihMPC_W9);
eihMPC_LA_DIAG_FORWARDSUB_7(eihMPC_Phi9, eihMPC_rd9, eihMPC_Lbyrd9);
eihMPC_LA_DIAGZERO_MMT_7(eihMPC_W0, eihMPC_Yd0);
eihMPC_LA_DIAGZERO_MVMSUB7_7(eihMPC_W0, eihMPC_Lbyrd0, eihMPC_re0, eihMPC_beta0);
eihMPC_LA_DENSE_DIAGZERO_MMT2_7_14_14(eihMPC_V0, eihMPC_W1, eihMPC_Yd1);
eihMPC_LA_DENSE_DIAGZERO_2MVMSUB2_7_14_14(eihMPC_V0, eihMPC_Lbyrd0, eihMPC_W1, eihMPC_Lbyrd1, eihMPC_re1, eihMPC_beta1);
eihMPC_LA_DENSE_DIAGZERO_MMT2_7_14_14(eihMPC_V1, eihMPC_W2, eihMPC_Yd2);
eihMPC_LA_DENSE_DIAGZERO_2MVMSUB2_7_14_14(eihMPC_V1, eihMPC_Lbyrd1, eihMPC_W2, eihMPC_Lbyrd2, eihMPC_re2, eihMPC_beta2);
eihMPC_LA_DENSE_DIAGZERO_MMT2_7_14_14(eihMPC_V2, eihMPC_W3, eihMPC_Yd3);
eihMPC_LA_DENSE_DIAGZERO_2MVMSUB2_7_14_14(eihMPC_V2, eihMPC_Lbyrd2, eihMPC_W3, eihMPC_Lbyrd3, eihMPC_re3, eihMPC_beta3);
eihMPC_LA_DENSE_DIAGZERO_MMT2_7_14_14(eihMPC_V3, eihMPC_W4, eihMPC_Yd4);
eihMPC_LA_DENSE_DIAGZERO_2MVMSUB2_7_14_14(eihMPC_V3, eihMPC_Lbyrd3, eihMPC_W4, eihMPC_Lbyrd4, eihMPC_re4, eihMPC_beta4);
eihMPC_LA_DENSE_DIAGZERO_MMT2_7_14_14(eihMPC_V4, eihMPC_W5, eihMPC_Yd5);
eihMPC_LA_DENSE_DIAGZERO_2MVMSUB2_7_14_14(eihMPC_V4, eihMPC_Lbyrd4, eihMPC_W5, eihMPC_Lbyrd5, eihMPC_re5, eihMPC_beta5);
eihMPC_LA_DENSE_DIAGZERO_MMT2_7_14_14(eihMPC_V5, eihMPC_W6, eihMPC_Yd6);
eihMPC_LA_DENSE_DIAGZERO_2MVMSUB2_7_14_14(eihMPC_V5, eihMPC_Lbyrd5, eihMPC_W6, eihMPC_Lbyrd6, eihMPC_re6, eihMPC_beta6);
eihMPC_LA_DENSE_DIAGZERO_MMT2_7_14_14(eihMPC_V6, eihMPC_W7, eihMPC_Yd7);
eihMPC_LA_DENSE_DIAGZERO_2MVMSUB2_7_14_14(eihMPC_V6, eihMPC_Lbyrd6, eihMPC_W7, eihMPC_Lbyrd7, eihMPC_re7, eihMPC_beta7);
eihMPC_LA_DENSE_DIAGZERO_MMT2_7_14_14(eihMPC_V7, eihMPC_W8, eihMPC_Yd8);
eihMPC_LA_DENSE_DIAGZERO_2MVMSUB2_7_14_14(eihMPC_V7, eihMPC_Lbyrd7, eihMPC_W8, eihMPC_Lbyrd8, eihMPC_re8, eihMPC_beta8);
eihMPC_LA_DENSE_DIAGZERO_MMT2_7_14_7(eihMPC_V8, eihMPC_W9, eihMPC_Yd9);
eihMPC_LA_DENSE_DIAGZERO_2MVMSUB2_7_14_7(eihMPC_V8, eihMPC_Lbyrd8, eihMPC_W9, eihMPC_Lbyrd9, eihMPC_re9, eihMPC_beta9);
eihMPC_LA_DENSE_CHOL_7(eihMPC_Yd0, eihMPC_Ld0);
eihMPC_LA_DENSE_FORWARDSUB_7(eihMPC_Ld0, eihMPC_beta0, eihMPC_yy0);
eihMPC_LA_DENSE_MATRIXTFORWARDSUB_7_7(eihMPC_Ld0, eihMPC_Ysd1, eihMPC_Lsd1);
eihMPC_LA_DENSE_MMTSUB_7_7(eihMPC_Lsd1, eihMPC_Yd1);
eihMPC_LA_DENSE_CHOL_7(eihMPC_Yd1, eihMPC_Ld1);
eihMPC_LA_DENSE_MVMSUB1_7_7(eihMPC_Lsd1, eihMPC_yy0, eihMPC_beta1, eihMPC_bmy1);
eihMPC_LA_DENSE_FORWARDSUB_7(eihMPC_Ld1, eihMPC_bmy1, eihMPC_yy1);
eihMPC_LA_DENSE_MATRIXTFORWARDSUB_7_7(eihMPC_Ld1, eihMPC_Ysd2, eihMPC_Lsd2);
eihMPC_LA_DENSE_MMTSUB_7_7(eihMPC_Lsd2, eihMPC_Yd2);
eihMPC_LA_DENSE_CHOL_7(eihMPC_Yd2, eihMPC_Ld2);
eihMPC_LA_DENSE_MVMSUB1_7_7(eihMPC_Lsd2, eihMPC_yy1, eihMPC_beta2, eihMPC_bmy2);
eihMPC_LA_DENSE_FORWARDSUB_7(eihMPC_Ld2, eihMPC_bmy2, eihMPC_yy2);
eihMPC_LA_DENSE_MATRIXTFORWARDSUB_7_7(eihMPC_Ld2, eihMPC_Ysd3, eihMPC_Lsd3);
eihMPC_LA_DENSE_MMTSUB_7_7(eihMPC_Lsd3, eihMPC_Yd3);
eihMPC_LA_DENSE_CHOL_7(eihMPC_Yd3, eihMPC_Ld3);
eihMPC_LA_DENSE_MVMSUB1_7_7(eihMPC_Lsd3, eihMPC_yy2, eihMPC_beta3, eihMPC_bmy3);
eihMPC_LA_DENSE_FORWARDSUB_7(eihMPC_Ld3, eihMPC_bmy3, eihMPC_yy3);
eihMPC_LA_DENSE_MATRIXTFORWARDSUB_7_7(eihMPC_Ld3, eihMPC_Ysd4, eihMPC_Lsd4);
eihMPC_LA_DENSE_MMTSUB_7_7(eihMPC_Lsd4, eihMPC_Yd4);
eihMPC_LA_DENSE_CHOL_7(eihMPC_Yd4, eihMPC_Ld4);
eihMPC_LA_DENSE_MVMSUB1_7_7(eihMPC_Lsd4, eihMPC_yy3, eihMPC_beta4, eihMPC_bmy4);
eihMPC_LA_DENSE_FORWARDSUB_7(eihMPC_Ld4, eihMPC_bmy4, eihMPC_yy4);
eihMPC_LA_DENSE_MATRIXTFORWARDSUB_7_7(eihMPC_Ld4, eihMPC_Ysd5, eihMPC_Lsd5);
eihMPC_LA_DENSE_MMTSUB_7_7(eihMPC_Lsd5, eihMPC_Yd5);
eihMPC_LA_DENSE_CHOL_7(eihMPC_Yd5, eihMPC_Ld5);
eihMPC_LA_DENSE_MVMSUB1_7_7(eihMPC_Lsd5, eihMPC_yy4, eihMPC_beta5, eihMPC_bmy5);
eihMPC_LA_DENSE_FORWARDSUB_7(eihMPC_Ld5, eihMPC_bmy5, eihMPC_yy5);
eihMPC_LA_DENSE_MATRIXTFORWARDSUB_7_7(eihMPC_Ld5, eihMPC_Ysd6, eihMPC_Lsd6);
eihMPC_LA_DENSE_MMTSUB_7_7(eihMPC_Lsd6, eihMPC_Yd6);
eihMPC_LA_DENSE_CHOL_7(eihMPC_Yd6, eihMPC_Ld6);
eihMPC_LA_DENSE_MVMSUB1_7_7(eihMPC_Lsd6, eihMPC_yy5, eihMPC_beta6, eihMPC_bmy6);
eihMPC_LA_DENSE_FORWARDSUB_7(eihMPC_Ld6, eihMPC_bmy6, eihMPC_yy6);
eihMPC_LA_DENSE_MATRIXTFORWARDSUB_7_7(eihMPC_Ld6, eihMPC_Ysd7, eihMPC_Lsd7);
eihMPC_LA_DENSE_MMTSUB_7_7(eihMPC_Lsd7, eihMPC_Yd7);
eihMPC_LA_DENSE_CHOL_7(eihMPC_Yd7, eihMPC_Ld7);
eihMPC_LA_DENSE_MVMSUB1_7_7(eihMPC_Lsd7, eihMPC_yy6, eihMPC_beta7, eihMPC_bmy7);
eihMPC_LA_DENSE_FORWARDSUB_7(eihMPC_Ld7, eihMPC_bmy7, eihMPC_yy7);
eihMPC_LA_DENSE_MATRIXTFORWARDSUB_7_7(eihMPC_Ld7, eihMPC_Ysd8, eihMPC_Lsd8);
eihMPC_LA_DENSE_MMTSUB_7_7(eihMPC_Lsd8, eihMPC_Yd8);
eihMPC_LA_DENSE_CHOL_7(eihMPC_Yd8, eihMPC_Ld8);
eihMPC_LA_DENSE_MVMSUB1_7_7(eihMPC_Lsd8, eihMPC_yy7, eihMPC_beta8, eihMPC_bmy8);
eihMPC_LA_DENSE_FORWARDSUB_7(eihMPC_Ld8, eihMPC_bmy8, eihMPC_yy8);
eihMPC_LA_DENSE_MATRIXTFORWARDSUB_7_7(eihMPC_Ld8, eihMPC_Ysd9, eihMPC_Lsd9);
eihMPC_LA_DENSE_MMTSUB_7_7(eihMPC_Lsd9, eihMPC_Yd9);
eihMPC_LA_DENSE_CHOL_7(eihMPC_Yd9, eihMPC_Ld9);
eihMPC_LA_DENSE_MVMSUB1_7_7(eihMPC_Lsd9, eihMPC_yy8, eihMPC_beta9, eihMPC_bmy9);
eihMPC_LA_DENSE_FORWARDSUB_7(eihMPC_Ld9, eihMPC_bmy9, eihMPC_yy9);
eihMPC_LA_DENSE_BACKWARDSUB_7(eihMPC_Ld9, eihMPC_yy9, eihMPC_dvaff9);
eihMPC_LA_DENSE_MTVMSUB_7_7(eihMPC_Lsd9, eihMPC_dvaff9, eihMPC_yy8, eihMPC_bmy8);
eihMPC_LA_DENSE_BACKWARDSUB_7(eihMPC_Ld8, eihMPC_bmy8, eihMPC_dvaff8);
eihMPC_LA_DENSE_MTVMSUB_7_7(eihMPC_Lsd8, eihMPC_dvaff8, eihMPC_yy7, eihMPC_bmy7);
eihMPC_LA_DENSE_BACKWARDSUB_7(eihMPC_Ld7, eihMPC_bmy7, eihMPC_dvaff7);
eihMPC_LA_DENSE_MTVMSUB_7_7(eihMPC_Lsd7, eihMPC_dvaff7, eihMPC_yy6, eihMPC_bmy6);
eihMPC_LA_DENSE_BACKWARDSUB_7(eihMPC_Ld6, eihMPC_bmy6, eihMPC_dvaff6);
eihMPC_LA_DENSE_MTVMSUB_7_7(eihMPC_Lsd6, eihMPC_dvaff6, eihMPC_yy5, eihMPC_bmy5);
eihMPC_LA_DENSE_BACKWARDSUB_7(eihMPC_Ld5, eihMPC_bmy5, eihMPC_dvaff5);
eihMPC_LA_DENSE_MTVMSUB_7_7(eihMPC_Lsd5, eihMPC_dvaff5, eihMPC_yy4, eihMPC_bmy4);
eihMPC_LA_DENSE_BACKWARDSUB_7(eihMPC_Ld4, eihMPC_bmy4, eihMPC_dvaff4);
eihMPC_LA_DENSE_MTVMSUB_7_7(eihMPC_Lsd4, eihMPC_dvaff4, eihMPC_yy3, eihMPC_bmy3);
eihMPC_LA_DENSE_BACKWARDSUB_7(eihMPC_Ld3, eihMPC_bmy3, eihMPC_dvaff3);
eihMPC_LA_DENSE_MTVMSUB_7_7(eihMPC_Lsd3, eihMPC_dvaff3, eihMPC_yy2, eihMPC_bmy2);
eihMPC_LA_DENSE_BACKWARDSUB_7(eihMPC_Ld2, eihMPC_bmy2, eihMPC_dvaff2);
eihMPC_LA_DENSE_MTVMSUB_7_7(eihMPC_Lsd2, eihMPC_dvaff2, eihMPC_yy1, eihMPC_bmy1);
eihMPC_LA_DENSE_BACKWARDSUB_7(eihMPC_Ld1, eihMPC_bmy1, eihMPC_dvaff1);
eihMPC_LA_DENSE_MTVMSUB_7_7(eihMPC_Lsd1, eihMPC_dvaff1, eihMPC_yy0, eihMPC_bmy0);
eihMPC_LA_DENSE_BACKWARDSUB_7(eihMPC_Ld0, eihMPC_bmy0, eihMPC_dvaff0);
eihMPC_LA_DENSE_DIAGZERO_MTVM2_7_14_7(eihMPC_C0, eihMPC_dvaff1, eihMPC_D0, eihMPC_dvaff0, eihMPC_grad_eq0);
eihMPC_LA_DENSE_DIAGZERO_MTVM2_7_14_7(eihMPC_C0, eihMPC_dvaff2, eihMPC_D1, eihMPC_dvaff1, eihMPC_grad_eq1);
eihMPC_LA_DENSE_DIAGZERO_MTVM2_7_14_7(eihMPC_C0, eihMPC_dvaff3, eihMPC_D1, eihMPC_dvaff2, eihMPC_grad_eq2);
eihMPC_LA_DENSE_DIAGZERO_MTVM2_7_14_7(eihMPC_C0, eihMPC_dvaff4, eihMPC_D1, eihMPC_dvaff3, eihMPC_grad_eq3);
eihMPC_LA_DENSE_DIAGZERO_MTVM2_7_14_7(eihMPC_C0, eihMPC_dvaff5, eihMPC_D1, eihMPC_dvaff4, eihMPC_grad_eq4);
eihMPC_LA_DENSE_DIAGZERO_MTVM2_7_14_7(eihMPC_C0, eihMPC_dvaff6, eihMPC_D1, eihMPC_dvaff5, eihMPC_grad_eq5);
eihMPC_LA_DENSE_DIAGZERO_MTVM2_7_14_7(eihMPC_C0, eihMPC_dvaff7, eihMPC_D1, eihMPC_dvaff6, eihMPC_grad_eq6);
eihMPC_LA_DENSE_DIAGZERO_MTVM2_7_14_7(eihMPC_C0, eihMPC_dvaff8, eihMPC_D1, eihMPC_dvaff7, eihMPC_grad_eq7);
eihMPC_LA_DENSE_DIAGZERO_MTVM2_7_14_7(eihMPC_C0, eihMPC_dvaff9, eihMPC_D1, eihMPC_dvaff8, eihMPC_grad_eq8);
eihMPC_LA_DIAGZERO_MTVM_7_7(eihMPC_D9, eihMPC_dvaff9, eihMPC_grad_eq9);
eihMPC_LA_VSUB2_133(eihMPC_rd, eihMPC_grad_eq, eihMPC_rd);
eihMPC_LA_DIAG_FORWARDBACKWARDSUB_14(eihMPC_Phi0, eihMPC_rd0, eihMPC_dzaff0);
eihMPC_LA_DIAG_FORWARDBACKWARDSUB_14(eihMPC_Phi1, eihMPC_rd1, eihMPC_dzaff1);
eihMPC_LA_DIAG_FORWARDBACKWARDSUB_14(eihMPC_Phi2, eihMPC_rd2, eihMPC_dzaff2);
eihMPC_LA_DIAG_FORWARDBACKWARDSUB_14(eihMPC_Phi3, eihMPC_rd3, eihMPC_dzaff3);
eihMPC_LA_DIAG_FORWARDBACKWARDSUB_14(eihMPC_Phi4, eihMPC_rd4, eihMPC_dzaff4);
eihMPC_LA_DIAG_FORWARDBACKWARDSUB_14(eihMPC_Phi5, eihMPC_rd5, eihMPC_dzaff5);
eihMPC_LA_DIAG_FORWARDBACKWARDSUB_14(eihMPC_Phi6, eihMPC_rd6, eihMPC_dzaff6);
eihMPC_LA_DIAG_FORWARDBACKWARDSUB_14(eihMPC_Phi7, eihMPC_rd7, eihMPC_dzaff7);
eihMPC_LA_DIAG_FORWARDBACKWARDSUB_14(eihMPC_Phi8, eihMPC_rd8, eihMPC_dzaff8);
eihMPC_LA_DIAG_FORWARDBACKWARDSUB_7(eihMPC_Phi9, eihMPC_rd9, eihMPC_dzaff9);
eihMPC_LA_VSUB_INDEXED_14(eihMPC_dzaff0, eihMPC_lbIdx0, eihMPC_rilb0, eihMPC_dslbaff0);
eihMPC_LA_VSUB3_14(eihMPC_llbbyslb0, eihMPC_dslbaff0, eihMPC_llb0, eihMPC_dllbaff0);
eihMPC_LA_VSUB2_INDEXED_14(eihMPC_riub0, eihMPC_dzaff0, eihMPC_ubIdx0, eihMPC_dsubaff0);
eihMPC_LA_VSUB3_14(eihMPC_lubbysub0, eihMPC_dsubaff0, eihMPC_lub0, eihMPC_dlubaff0);
eihMPC_LA_VSUB_INDEXED_14(eihMPC_dzaff1, eihMPC_lbIdx1, eihMPC_rilb1, eihMPC_dslbaff1);
eihMPC_LA_VSUB3_14(eihMPC_llbbyslb1, eihMPC_dslbaff1, eihMPC_llb1, eihMPC_dllbaff1);
eihMPC_LA_VSUB2_INDEXED_14(eihMPC_riub1, eihMPC_dzaff1, eihMPC_ubIdx1, eihMPC_dsubaff1);
eihMPC_LA_VSUB3_14(eihMPC_lubbysub1, eihMPC_dsubaff1, eihMPC_lub1, eihMPC_dlubaff1);
eihMPC_LA_VSUB_INDEXED_14(eihMPC_dzaff2, eihMPC_lbIdx2, eihMPC_rilb2, eihMPC_dslbaff2);
eihMPC_LA_VSUB3_14(eihMPC_llbbyslb2, eihMPC_dslbaff2, eihMPC_llb2, eihMPC_dllbaff2);
eihMPC_LA_VSUB2_INDEXED_14(eihMPC_riub2, eihMPC_dzaff2, eihMPC_ubIdx2, eihMPC_dsubaff2);
eihMPC_LA_VSUB3_14(eihMPC_lubbysub2, eihMPC_dsubaff2, eihMPC_lub2, eihMPC_dlubaff2);
eihMPC_LA_VSUB_INDEXED_14(eihMPC_dzaff3, eihMPC_lbIdx3, eihMPC_rilb3, eihMPC_dslbaff3);
eihMPC_LA_VSUB3_14(eihMPC_llbbyslb3, eihMPC_dslbaff3, eihMPC_llb3, eihMPC_dllbaff3);
eihMPC_LA_VSUB2_INDEXED_14(eihMPC_riub3, eihMPC_dzaff3, eihMPC_ubIdx3, eihMPC_dsubaff3);
eihMPC_LA_VSUB3_14(eihMPC_lubbysub3, eihMPC_dsubaff3, eihMPC_lub3, eihMPC_dlubaff3);
eihMPC_LA_VSUB_INDEXED_14(eihMPC_dzaff4, eihMPC_lbIdx4, eihMPC_rilb4, eihMPC_dslbaff4);
eihMPC_LA_VSUB3_14(eihMPC_llbbyslb4, eihMPC_dslbaff4, eihMPC_llb4, eihMPC_dllbaff4);
eihMPC_LA_VSUB2_INDEXED_14(eihMPC_riub4, eihMPC_dzaff4, eihMPC_ubIdx4, eihMPC_dsubaff4);
eihMPC_LA_VSUB3_14(eihMPC_lubbysub4, eihMPC_dsubaff4, eihMPC_lub4, eihMPC_dlubaff4);
eihMPC_LA_VSUB_INDEXED_14(eihMPC_dzaff5, eihMPC_lbIdx5, eihMPC_rilb5, eihMPC_dslbaff5);
eihMPC_LA_VSUB3_14(eihMPC_llbbyslb5, eihMPC_dslbaff5, eihMPC_llb5, eihMPC_dllbaff5);
eihMPC_LA_VSUB2_INDEXED_14(eihMPC_riub5, eihMPC_dzaff5, eihMPC_ubIdx5, eihMPC_dsubaff5);
eihMPC_LA_VSUB3_14(eihMPC_lubbysub5, eihMPC_dsubaff5, eihMPC_lub5, eihMPC_dlubaff5);
eihMPC_LA_VSUB_INDEXED_14(eihMPC_dzaff6, eihMPC_lbIdx6, eihMPC_rilb6, eihMPC_dslbaff6);
eihMPC_LA_VSUB3_14(eihMPC_llbbyslb6, eihMPC_dslbaff6, eihMPC_llb6, eihMPC_dllbaff6);
eihMPC_LA_VSUB2_INDEXED_14(eihMPC_riub6, eihMPC_dzaff6, eihMPC_ubIdx6, eihMPC_dsubaff6);
eihMPC_LA_VSUB3_14(eihMPC_lubbysub6, eihMPC_dsubaff6, eihMPC_lub6, eihMPC_dlubaff6);
eihMPC_LA_VSUB_INDEXED_14(eihMPC_dzaff7, eihMPC_lbIdx7, eihMPC_rilb7, eihMPC_dslbaff7);
eihMPC_LA_VSUB3_14(eihMPC_llbbyslb7, eihMPC_dslbaff7, eihMPC_llb7, eihMPC_dllbaff7);
eihMPC_LA_VSUB2_INDEXED_14(eihMPC_riub7, eihMPC_dzaff7, eihMPC_ubIdx7, eihMPC_dsubaff7);
eihMPC_LA_VSUB3_14(eihMPC_lubbysub7, eihMPC_dsubaff7, eihMPC_lub7, eihMPC_dlubaff7);
eihMPC_LA_VSUB_INDEXED_14(eihMPC_dzaff8, eihMPC_lbIdx8, eihMPC_rilb8, eihMPC_dslbaff8);
eihMPC_LA_VSUB3_14(eihMPC_llbbyslb8, eihMPC_dslbaff8, eihMPC_llb8, eihMPC_dllbaff8);
eihMPC_LA_VSUB2_INDEXED_14(eihMPC_riub8, eihMPC_dzaff8, eihMPC_ubIdx8, eihMPC_dsubaff8);
eihMPC_LA_VSUB3_14(eihMPC_lubbysub8, eihMPC_dsubaff8, eihMPC_lub8, eihMPC_dlubaff8);
eihMPC_LA_VSUB_INDEXED_7(eihMPC_dzaff9, eihMPC_lbIdx9, eihMPC_rilb9, eihMPC_dslbaff9);
eihMPC_LA_VSUB3_7(eihMPC_llbbyslb9, eihMPC_dslbaff9, eihMPC_llb9, eihMPC_dllbaff9);
eihMPC_LA_VSUB2_INDEXED_7(eihMPC_riub9, eihMPC_dzaff9, eihMPC_ubIdx9, eihMPC_dsubaff9);
eihMPC_LA_VSUB3_7(eihMPC_lubbysub9, eihMPC_dsubaff9, eihMPC_lub9, eihMPC_dlubaff9);
info->lsit_aff = eihMPC_LINESEARCH_BACKTRACKING_AFFINE(eihMPC_l, eihMPC_s, eihMPC_dl_aff, eihMPC_ds_aff, &info->step_aff, &info->mu_aff);
if( info->lsit_aff == eihMPC_NOPROGRESS ){
exitcode = eihMPC_NOPROGRESS; break;
}
sigma_3rdroot = info->mu_aff / info->mu;
info->sigma = sigma_3rdroot*sigma_3rdroot*sigma_3rdroot;
musigma = info->mu * info->sigma;
eihMPC_LA_VSUB5_266(eihMPC_ds_aff, eihMPC_dl_aff, musigma, eihMPC_ccrhs);
eihMPC_LA_VSUB6_INDEXED_14_14_14(eihMPC_ccrhsub0, eihMPC_sub0, eihMPC_ubIdx0, eihMPC_ccrhsl0, eihMPC_slb0, eihMPC_lbIdx0, eihMPC_rd0);
eihMPC_LA_VSUB6_INDEXED_14_14_14(eihMPC_ccrhsub1, eihMPC_sub1, eihMPC_ubIdx1, eihMPC_ccrhsl1, eihMPC_slb1, eihMPC_lbIdx1, eihMPC_rd1);
eihMPC_LA_DIAG_FORWARDSUB_14(eihMPC_Phi0, eihMPC_rd0, eihMPC_Lbyrd0);
eihMPC_LA_DIAG_FORWARDSUB_14(eihMPC_Phi1, eihMPC_rd1, eihMPC_Lbyrd1);
eihMPC_LA_DIAGZERO_MVM_7(eihMPC_W0, eihMPC_Lbyrd0, eihMPC_beta0);
eihMPC_LA_DENSE_FORWARDSUB_7(eihMPC_Ld0, eihMPC_beta0, eihMPC_yy0);
eihMPC_LA_DENSE_DIAGZERO_2MVMADD_7_14_14(eihMPC_V0, eihMPC_Lbyrd0, eihMPC_W1, eihMPC_Lbyrd1, eihMPC_beta1);
eihMPC_LA_DENSE_MVMSUB1_7_7(eihMPC_Lsd1, eihMPC_yy0, eihMPC_beta1, eihMPC_bmy1);
eihMPC_LA_DENSE_FORWARDSUB_7(eihMPC_Ld1, eihMPC_bmy1, eihMPC_yy1);
eihMPC_LA_VSUB6_INDEXED_14_14_14(eihMPC_ccrhsub2, eihMPC_sub2, eihMPC_ubIdx2, eihMPC_ccrhsl2, eihMPC_slb2, eihMPC_lbIdx2, eihMPC_rd2);
eihMPC_LA_DIAG_FORWARDSUB_14(eihMPC_Phi2, eihMPC_rd2, eihMPC_Lbyrd2);
eihMPC_LA_DENSE_DIAGZERO_2MVMADD_7_14_14(eihMPC_V1, eihMPC_Lbyrd1, eihMPC_W2, eihMPC_Lbyrd2, eihMPC_beta2);
eihMPC_LA_DENSE_MVMSUB1_7_7(eihMPC_Lsd2, eihMPC_yy1, eihMPC_beta2, eihMPC_bmy2);
eihMPC_LA_DENSE_FORWARDSUB_7(eihMPC_Ld2, eihMPC_bmy2, eihMPC_yy2);
eihMPC_LA_VSUB6_INDEXED_14_14_14(eihMPC_ccrhsub3, eihMPC_sub3, eihMPC_ubIdx3, eihMPC_ccrhsl3, eihMPC_slb3, eihMPC_lbIdx3, eihMPC_rd3);
eihMPC_LA_DIAG_FORWARDSUB_14(eihMPC_Phi3, eihMPC_rd3, eihMPC_Lbyrd3);
eihMPC_LA_DENSE_DIAGZERO_2MVMADD_7_14_14(eihMPC_V2, eihMPC_Lbyrd2, eihMPC_W3, eihMPC_Lbyrd3, eihMPC_beta3);
eihMPC_LA_DENSE_MVMSUB1_7_7(eihMPC_Lsd3, eihMPC_yy2, eihMPC_beta3, eihMPC_bmy3);
eihMPC_LA_DENSE_FORWARDSUB_7(eihMPC_Ld3, eihMPC_bmy3, eihMPC_yy3);
eihMPC_LA_VSUB6_INDEXED_14_14_14(eihMPC_ccrhsub4, eihMPC_sub4, eihMPC_ubIdx4, eihMPC_ccrhsl4, eihMPC_slb4, eihMPC_lbIdx4, eihMPC_rd4);
eihMPC_LA_DIAG_FORWARDSUB_14(eihMPC_Phi4, eihMPC_rd4, eihMPC_Lbyrd4);
eihMPC_LA_DENSE_DIAGZERO_2MVMADD_7_14_14(eihMPC_V3, eihMPC_Lbyrd3, eihMPC_W4, eihMPC_Lbyrd4, eihMPC_beta4);
eihMPC_LA_DENSE_MVMSUB1_7_7(eihMPC_Lsd4, eihMPC_yy3, eihMPC_beta4, eihMPC_bmy4);
eihMPC_LA_DENSE_FORWARDSUB_7(eihMPC_Ld4, eihMPC_bmy4, eihMPC_yy4);
eihMPC_LA_VSUB6_INDEXED_14_14_14(eihMPC_ccrhsub5, eihMPC_sub5, eihMPC_ubIdx5, eihMPC_ccrhsl5, eihMPC_slb5, eihMPC_lbIdx5, eihMPC_rd5);
eihMPC_LA_DIAG_FORWARDSUB_14(eihMPC_Phi5, eihMPC_rd5, eihMPC_Lbyrd5);
eihMPC_LA_DENSE_DIAGZERO_2MVMADD_7_14_14(eihMPC_V4, eihMPC_Lbyrd4, eihMPC_W5, eihMPC_Lbyrd5, eihMPC_beta5);
eihMPC_LA_DENSE_MVMSUB1_7_7(eihMPC_Lsd5, eihMPC_yy4, eihMPC_beta5, eihMPC_bmy5);
eihMPC_LA_DENSE_FORWARDSUB_7(eihMPC_Ld5, eihMPC_bmy5, eihMPC_yy5);
eihMPC_LA_VSUB6_INDEXED_14_14_14(eihMPC_ccrhsub6, eihMPC_sub6, eihMPC_ubIdx6, eihMPC_ccrhsl6, eihMPC_slb6, eihMPC_lbIdx6, eihMPC_rd6);
eihMPC_LA_DIAG_FORWARDSUB_14(eihMPC_Phi6, eihMPC_rd6, eihMPC_Lbyrd6);
eihMPC_LA_DENSE_DIAGZERO_2MVMADD_7_14_14(eihMPC_V5, eihMPC_Lbyrd5, eihMPC_W6, eihMPC_Lbyrd6, eihMPC_beta6);
eihMPC_LA_DENSE_MVMSUB1_7_7(eihMPC_Lsd6, eihMPC_yy5, eihMPC_beta6, eihMPC_bmy6);
eihMPC_LA_DENSE_FORWARDSUB_7(eihMPC_Ld6, eihMPC_bmy6, eihMPC_yy6);
eihMPC_LA_VSUB6_INDEXED_14_14_14(eihMPC_ccrhsub7, eihMPC_sub7, eihMPC_ubIdx7, eihMPC_ccrhsl7, eihMPC_slb7, eihMPC_lbIdx7, eihMPC_rd7);
eihMPC_LA_DIAG_FORWARDSUB_14(eihMPC_Phi7, eihMPC_rd7, eihMPC_Lbyrd7);
eihMPC_LA_DENSE_DIAGZERO_2MVMADD_7_14_14(eihMPC_V6, eihMPC_Lbyrd6, eihMPC_W7, eihMPC_Lbyrd7, eihMPC_beta7);
eihMPC_LA_DENSE_MVMSUB1_7_7(eihMPC_Lsd7, eihMPC_yy6, eihMPC_beta7, eihMPC_bmy7);
eihMPC_LA_DENSE_FORWARDSUB_7(eihMPC_Ld7, eihMPC_bmy7, eihMPC_yy7);
eihMPC_LA_VSUB6_INDEXED_14_14_14(eihMPC_ccrhsub8, eihMPC_sub8, eihMPC_ubIdx8, eihMPC_ccrhsl8, eihMPC_slb8, eihMPC_lbIdx8, eihMPC_rd8);
eihMPC_LA_DIAG_FORWARDSUB_14(eihMPC_Phi8, eihMPC_rd8, eihMPC_Lbyrd8);
eihMPC_LA_DENSE_DIAGZERO_2MVMADD_7_14_14(eihMPC_V7, eihMPC_Lbyrd7, eihMPC_W8, eihMPC_Lbyrd8, eihMPC_beta8);
eihMPC_LA_DENSE_MVMSUB1_7_7(eihMPC_Lsd8, eihMPC_yy7, eihMPC_beta8, eihMPC_bmy8);
eihMPC_LA_DENSE_FORWARDSUB_7(eihMPC_Ld8, eihMPC_bmy8, eihMPC_yy8);
eihMPC_LA_VSUB6_INDEXED_7_7_7(eihMPC_ccrhsub9, eihMPC_sub9, eihMPC_ubIdx9, eihMPC_ccrhsl9, eihMPC_slb9, eihMPC_lbIdx9, eihMPC_rd9);
eihMPC_LA_DIAG_FORWARDSUB_7(eihMPC_Phi9, eihMPC_rd9, eihMPC_Lbyrd9);
eihMPC_LA_DENSE_DIAGZERO_2MVMADD_7_14_7(eihMPC_V8, eihMPC_Lbyrd8, eihMPC_W9, eihMPC_Lbyrd9, eihMPC_beta9);
eihMPC_LA_DENSE_MVMSUB1_7_7(eihMPC_Lsd9, eihMPC_yy8, eihMPC_beta9, eihMPC_bmy9);
eihMPC_LA_DENSE_FORWARDSUB_7(eihMPC_Ld9, eihMPC_bmy9, eihMPC_yy9);
eihMPC_LA_DENSE_BACKWARDSUB_7(eihMPC_Ld9, eihMPC_yy9, eihMPC_dvcc9);
eihMPC_LA_DENSE_MTVMSUB_7_7(eihMPC_Lsd9, eihMPC_dvcc9, eihMPC_yy8, eihMPC_bmy8);
eihMPC_LA_DENSE_BACKWARDSUB_7(eihMPC_Ld8, eihMPC_bmy8, eihMPC_dvcc8);
eihMPC_LA_DENSE_MTVMSUB_7_7(eihMPC_Lsd8, eihMPC_dvcc8, eihMPC_yy7, eihMPC_bmy7);
eihMPC_LA_DENSE_BACKWARDSUB_7(eihMPC_Ld7, eihMPC_bmy7, eihMPC_dvcc7);
eihMPC_LA_DENSE_MTVMSUB_7_7(eihMPC_Lsd7, eihMPC_dvcc7, eihMPC_yy6, eihMPC_bmy6);
eihMPC_LA_DENSE_BACKWARDSUB_7(eihMPC_Ld6, eihMPC_bmy6, eihMPC_dvcc6);
eihMPC_LA_DENSE_MTVMSUB_7_7(eihMPC_Lsd6, eihMPC_dvcc6, eihMPC_yy5, eihMPC_bmy5);
eihMPC_LA_DENSE_BACKWARDSUB_7(eihMPC_Ld5, eihMPC_bmy5, eihMPC_dvcc5);
eihMPC_LA_DENSE_MTVMSUB_7_7(eihMPC_Lsd5, eihMPC_dvcc5, eihMPC_yy4, eihMPC_bmy4);
eihMPC_LA_DENSE_BACKWARDSUB_7(eihMPC_Ld4, eihMPC_bmy4, eihMPC_dvcc4);
eihMPC_LA_DENSE_MTVMSUB_7_7(eihMPC_Lsd4, eihMPC_dvcc4, eihMPC_yy3, eihMPC_bmy3);
eihMPC_LA_DENSE_BACKWARDSUB_7(eihMPC_Ld3, eihMPC_bmy3, eihMPC_dvcc3);
eihMPC_LA_DENSE_MTVMSUB_7_7(eihMPC_Lsd3, eihMPC_dvcc3, eihMPC_yy2, eihMPC_bmy2);
eihMPC_LA_DENSE_BACKWARDSUB_7(eihMPC_Ld2, eihMPC_bmy2, eihMPC_dvcc2);
eihMPC_LA_DENSE_MTVMSUB_7_7(eihMPC_Lsd2, eihMPC_dvcc2, eihMPC_yy1, eihMPC_bmy1);
eihMPC_LA_DENSE_BACKWARDSUB_7(eihMPC_Ld1, eihMPC_bmy1, eihMPC_dvcc1);
eihMPC_LA_DENSE_MTVMSUB_7_7(eihMPC_Lsd1, eihMPC_dvcc1, eihMPC_yy0, eihMPC_bmy0);
eihMPC_LA_DENSE_BACKWARDSUB_7(eihMPC_Ld0, eihMPC_bmy0, eihMPC_dvcc0);
eihMPC_LA_DENSE_DIAGZERO_MTVM2_7_14_7(eihMPC_C0, eihMPC_dvcc1, eihMPC_D0, eihMPC_dvcc0, eihMPC_grad_eq0);
eihMPC_LA_DENSE_DIAGZERO_MTVM2_7_14_7(eihMPC_C0, eihMPC_dvcc2, eihMPC_D1, eihMPC_dvcc1, eihMPC_grad_eq1);
eihMPC_LA_DENSE_DIAGZERO_MTVM2_7_14_7(eihMPC_C0, eihMPC_dvcc3, eihMPC_D1, eihMPC_dvcc2, eihMPC_grad_eq2);
eihMPC_LA_DENSE_DIAGZERO_MTVM2_7_14_7(eihMPC_C0, eihMPC_dvcc4, eihMPC_D1, eihMPC_dvcc3, eihMPC_grad_eq3);
eihMPC_LA_DENSE_DIAGZERO_MTVM2_7_14_7(eihMPC_C0, eihMPC_dvcc5, eihMPC_D1, eihMPC_dvcc4, eihMPC_grad_eq4);
eihMPC_LA_DENSE_DIAGZERO_MTVM2_7_14_7(eihMPC_C0, eihMPC_dvcc6, eihMPC_D1, eihMPC_dvcc5, eihMPC_grad_eq5);
eihMPC_LA_DENSE_DIAGZERO_MTVM2_7_14_7(eihMPC_C0, eihMPC_dvcc7, eihMPC_D1, eihMPC_dvcc6, eihMPC_grad_eq6);
eihMPC_LA_DENSE_DIAGZERO_MTVM2_7_14_7(eihMPC_C0, eihMPC_dvcc8, eihMPC_D1, eihMPC_dvcc7, eihMPC_grad_eq7);
eihMPC_LA_DENSE_DIAGZERO_MTVM2_7_14_7(eihMPC_C0, eihMPC_dvcc9, eihMPC_D1, eihMPC_dvcc8, eihMPC_grad_eq8);
eihMPC_LA_DIAGZERO_MTVM_7_7(eihMPC_D9, eihMPC_dvcc9, eihMPC_grad_eq9);
eihMPC_LA_VSUB_133(eihMPC_rd, eihMPC_grad_eq, eihMPC_rd);
eihMPC_LA_DIAG_FORWARDBACKWARDSUB_14(eihMPC_Phi0, eihMPC_rd0, eihMPC_dzcc0);
eihMPC_LA_DIAG_FORWARDBACKWARDSUB_14(eihMPC_Phi1, eihMPC_rd1, eihMPC_dzcc1);
eihMPC_LA_DIAG_FORWARDBACKWARDSUB_14(eihMPC_Phi2, eihMPC_rd2, eihMPC_dzcc2);
eihMPC_LA_DIAG_FORWARDBACKWARDSUB_14(eihMPC_Phi3, eihMPC_rd3, eihMPC_dzcc3);
eihMPC_LA_DIAG_FORWARDBACKWARDSUB_14(eihMPC_Phi4, eihMPC_rd4, eihMPC_dzcc4);
eihMPC_LA_DIAG_FORWARDBACKWARDSUB_14(eihMPC_Phi5, eihMPC_rd5, eihMPC_dzcc5);
eihMPC_LA_DIAG_FORWARDBACKWARDSUB_14(eihMPC_Phi6, eihMPC_rd6, eihMPC_dzcc6);
eihMPC_LA_DIAG_FORWARDBACKWARDSUB_14(eihMPC_Phi7, eihMPC_rd7, eihMPC_dzcc7);
eihMPC_LA_DIAG_FORWARDBACKWARDSUB_14(eihMPC_Phi8, eihMPC_rd8, eihMPC_dzcc8);
eihMPC_LA_DIAG_FORWARDBACKWARDSUB_7(eihMPC_Phi9, eihMPC_rd9, eihMPC_dzcc9);
eihMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_14(eihMPC_ccrhsl0, eihMPC_slb0, eihMPC_llbbyslb0, eihMPC_dzcc0, eihMPC_lbIdx0, eihMPC_dllbcc0);
eihMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_14(eihMPC_ccrhsub0, eihMPC_sub0, eihMPC_lubbysub0, eihMPC_dzcc0, eihMPC_ubIdx0, eihMPC_dlubcc0);
eihMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_14(eihMPC_ccrhsl1, eihMPC_slb1, eihMPC_llbbyslb1, eihMPC_dzcc1, eihMPC_lbIdx1, eihMPC_dllbcc1);
eihMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_14(eihMPC_ccrhsub1, eihMPC_sub1, eihMPC_lubbysub1, eihMPC_dzcc1, eihMPC_ubIdx1, eihMPC_dlubcc1);
eihMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_14(eihMPC_ccrhsl2, eihMPC_slb2, eihMPC_llbbyslb2, eihMPC_dzcc2, eihMPC_lbIdx2, eihMPC_dllbcc2);
eihMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_14(eihMPC_ccrhsub2, eihMPC_sub2, eihMPC_lubbysub2, eihMPC_dzcc2, eihMPC_ubIdx2, eihMPC_dlubcc2);
eihMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_14(eihMPC_ccrhsl3, eihMPC_slb3, eihMPC_llbbyslb3, eihMPC_dzcc3, eihMPC_lbIdx3, eihMPC_dllbcc3);
eihMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_14(eihMPC_ccrhsub3, eihMPC_sub3, eihMPC_lubbysub3, eihMPC_dzcc3, eihMPC_ubIdx3, eihMPC_dlubcc3);
eihMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_14(eihMPC_ccrhsl4, eihMPC_slb4, eihMPC_llbbyslb4, eihMPC_dzcc4, eihMPC_lbIdx4, eihMPC_dllbcc4);
eihMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_14(eihMPC_ccrhsub4, eihMPC_sub4, eihMPC_lubbysub4, eihMPC_dzcc4, eihMPC_ubIdx4, eihMPC_dlubcc4);
eihMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_14(eihMPC_ccrhsl5, eihMPC_slb5, eihMPC_llbbyslb5, eihMPC_dzcc5, eihMPC_lbIdx5, eihMPC_dllbcc5);
eihMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_14(eihMPC_ccrhsub5, eihMPC_sub5, eihMPC_lubbysub5, eihMPC_dzcc5, eihMPC_ubIdx5, eihMPC_dlubcc5);
eihMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_14(eihMPC_ccrhsl6, eihMPC_slb6, eihMPC_llbbyslb6, eihMPC_dzcc6, eihMPC_lbIdx6, eihMPC_dllbcc6);
eihMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_14(eihMPC_ccrhsub6, eihMPC_sub6, eihMPC_lubbysub6, eihMPC_dzcc6, eihMPC_ubIdx6, eihMPC_dlubcc6);
eihMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_14(eihMPC_ccrhsl7, eihMPC_slb7, eihMPC_llbbyslb7, eihMPC_dzcc7, eihMPC_lbIdx7, eihMPC_dllbcc7);
eihMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_14(eihMPC_ccrhsub7, eihMPC_sub7, eihMPC_lubbysub7, eihMPC_dzcc7, eihMPC_ubIdx7, eihMPC_dlubcc7);
eihMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_14(eihMPC_ccrhsl8, eihMPC_slb8, eihMPC_llbbyslb8, eihMPC_dzcc8, eihMPC_lbIdx8, eihMPC_dllbcc8);
eihMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_14(eihMPC_ccrhsub8, eihMPC_sub8, eihMPC_lubbysub8, eihMPC_dzcc8, eihMPC_ubIdx8, eihMPC_dlubcc8);
eihMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_7(eihMPC_ccrhsl9, eihMPC_slb9, eihMPC_llbbyslb9, eihMPC_dzcc9, eihMPC_lbIdx9, eihMPC_dllbcc9);
eihMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_7(eihMPC_ccrhsub9, eihMPC_sub9, eihMPC_lubbysub9, eihMPC_dzcc9, eihMPC_ubIdx9, eihMPC_dlubcc9);
eihMPC_LA_VSUB7_266(eihMPC_l, eihMPC_ccrhs, eihMPC_s, eihMPC_dl_cc, eihMPC_ds_cc);
eihMPC_LA_VADD_133(eihMPC_dz_cc, eihMPC_dz_aff);
eihMPC_LA_VADD_70(eihMPC_dv_cc, eihMPC_dv_aff);
eihMPC_LA_VADD_266(eihMPC_dl_cc, eihMPC_dl_aff);
eihMPC_LA_VADD_266(eihMPC_ds_cc, eihMPC_ds_aff);
info->lsit_cc = eihMPC_LINESEARCH_BACKTRACKING_COMBINED(eihMPC_z, eihMPC_v, eihMPC_l, eihMPC_s, eihMPC_dz_cc, eihMPC_dv_cc, eihMPC_dl_cc, eihMPC_ds_cc, &info->step_cc, &info->mu);
if( info->lsit_cc == eihMPC_NOPROGRESS ){
exitcode = eihMPC_NOPROGRESS; break;
}
info->it++;
}
output->z1[0] = eihMPC_z0[0];
output->z1[1] = eihMPC_z0[1];
output->z1[2] = eihMPC_z0[2];
output->z1[3] = eihMPC_z0[3];
output->z1[4] = eihMPC_z0[4];
output->z1[5] = eihMPC_z0[5];
output->z1[6] = eihMPC_z0[6];
output->z1[7] = eihMPC_z0[7];
output->z1[8] = eihMPC_z0[8];
output->z1[9] = eihMPC_z0[9];
output->z1[10] = eihMPC_z0[10];
output->z1[11] = eihMPC_z0[11];
output->z1[12] = eihMPC_z0[12];
output->z1[13] = eihMPC_z0[13];
output->z2[0] = eihMPC_z1[0];
output->z2[1] = eihMPC_z1[1];
output->z2[2] = eihMPC_z1[2];
output->z2[3] = eihMPC_z1[3];
output->z2[4] = eihMPC_z1[4];
output->z2[5] = eihMPC_z1[5];
output->z2[6] = eihMPC_z1[6];
output->z2[7] = eihMPC_z1[7];
output->z2[8] = eihMPC_z1[8];
output->z2[9] = eihMPC_z1[9];
output->z2[10] = eihMPC_z1[10];
output->z2[11] = eihMPC_z1[11];
output->z2[12] = eihMPC_z1[12];
output->z2[13] = eihMPC_z1[13];
output->z3[0] = eihMPC_z2[0];
output->z3[1] = eihMPC_z2[1];
output->z3[2] = eihMPC_z2[2];
output->z3[3] = eihMPC_z2[3];
output->z3[4] = eihMPC_z2[4];
output->z3[5] = eihMPC_z2[5];
output->z3[6] = eihMPC_z2[6];
output->z3[7] = eihMPC_z2[7];
output->z3[8] = eihMPC_z2[8];
output->z3[9] = eihMPC_z2[9];
output->z3[10] = eihMPC_z2[10];
output->z3[11] = eihMPC_z2[11];
output->z3[12] = eihMPC_z2[12];
output->z3[13] = eihMPC_z2[13];
output->z4[0] = eihMPC_z3[0];
output->z4[1] = eihMPC_z3[1];
output->z4[2] = eihMPC_z3[2];
output->z4[3] = eihMPC_z3[3];
output->z4[4] = eihMPC_z3[4];
output->z4[5] = eihMPC_z3[5];
output->z4[6] = eihMPC_z3[6];
output->z4[7] = eihMPC_z3[7];
output->z4[8] = eihMPC_z3[8];
output->z4[9] = eihMPC_z3[9];
output->z4[10] = eihMPC_z3[10];
output->z4[11] = eihMPC_z3[11];
output->z4[12] = eihMPC_z3[12];
output->z4[13] = eihMPC_z3[13];
output->z5[0] = eihMPC_z4[0];
output->z5[1] = eihMPC_z4[1];
output->z5[2] = eihMPC_z4[2];
output->z5[3] = eihMPC_z4[3];
output->z5[4] = eihMPC_z4[4];
output->z5[5] = eihMPC_z4[5];
output->z5[6] = eihMPC_z4[6];
output->z5[7] = eihMPC_z4[7];
output->z5[8] = eihMPC_z4[8];
output->z5[9] = eihMPC_z4[9];
output->z5[10] = eihMPC_z4[10];
output->z5[11] = eihMPC_z4[11];
output->z5[12] = eihMPC_z4[12];
output->z5[13] = eihMPC_z4[13];
output->z6[0] = eihMPC_z5[0];
output->z6[1] = eihMPC_z5[1];
output->z6[2] = eihMPC_z5[2];
output->z6[3] = eihMPC_z5[3];
output->z6[4] = eihMPC_z5[4];
output->z6[5] = eihMPC_z5[5];
output->z6[6] = eihMPC_z5[6];
output->z6[7] = eihMPC_z5[7];
output->z6[8] = eihMPC_z5[8];
output->z6[9] = eihMPC_z5[9];
output->z6[10] = eihMPC_z5[10];
output->z6[11] = eihMPC_z5[11];
output->z6[12] = eihMPC_z5[12];
output->z6[13] = eihMPC_z5[13];
output->z7[0] = eihMPC_z6[0];
output->z7[1] = eihMPC_z6[1];
output->z7[2] = eihMPC_z6[2];
output->z7[3] = eihMPC_z6[3];
output->z7[4] = eihMPC_z6[4];
output->z7[5] = eihMPC_z6[5];
output->z7[6] = eihMPC_z6[6];
output->z7[7] = eihMPC_z6[7];
output->z7[8] = eihMPC_z6[8];
output->z7[9] = eihMPC_z6[9];
output->z7[10] = eihMPC_z6[10];
output->z7[11] = eihMPC_z6[11];
output->z7[12] = eihMPC_z6[12];
output->z7[13] = eihMPC_z6[13];
output->z8[0] = eihMPC_z7[0];
output->z8[1] = eihMPC_z7[1];
output->z8[2] = eihMPC_z7[2];
output->z8[3] = eihMPC_z7[3];
output->z8[4] = eihMPC_z7[4];
output->z8[5] = eihMPC_z7[5];
output->z8[6] = eihMPC_z7[6];
output->z8[7] = eihMPC_z7[7];
output->z8[8] = eihMPC_z7[8];
output->z8[9] = eihMPC_z7[9];
output->z8[10] = eihMPC_z7[10];
output->z8[11] = eihMPC_z7[11];
output->z8[12] = eihMPC_z7[12];
output->z8[13] = eihMPC_z7[13];
output->z9[0] = eihMPC_z8[0];
output->z9[1] = eihMPC_z8[1];
output->z9[2] = eihMPC_z8[2];
output->z9[3] = eihMPC_z8[3];
output->z9[4] = eihMPC_z8[4];
output->z9[5] = eihMPC_z8[5];
output->z9[6] = eihMPC_z8[6];
output->z9[7] = eihMPC_z8[7];
output->z9[8] = eihMPC_z8[8];
output->z9[9] = eihMPC_z8[9];
output->z9[10] = eihMPC_z8[10];
output->z9[11] = eihMPC_z8[11];
output->z9[12] = eihMPC_z8[12];
output->z9[13] = eihMPC_z8[13];
output->z10[0] = eihMPC_z9[0];
output->z10[1] = eihMPC_z9[1];
output->z10[2] = eihMPC_z9[2];
output->z10[3] = eihMPC_z9[3];
output->z10[4] = eihMPC_z9[4];
output->z10[5] = eihMPC_z9[5];
output->z10[6] = eihMPC_z9[6];

#if eihMPC_SET_TIMING == 1
info->solvetime = eihMPC_toc(&solvertimer);
#if eihMPC_SET_PRINTLEVEL > 0 && eihMPC_SET_TIMING == 1
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
