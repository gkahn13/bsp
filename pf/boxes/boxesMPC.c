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

#include "boxesMPC.h"

/* for square root */
#include <math.h> 

/* SAFE DIVISION ------------------------------------------------------- */
#define MAX(X,Y)  ((X) < (Y) ? (Y) : (X))
#define MIN(X,Y)  ((X) < (Y) ? (X) : (Y))
/*#define SAFEDIV_POS(X,Y)  ( (Y) < EPS ? ((X)/EPS) : (X)/(Y) ) 
#define EPS (1.0000E-013) */
#define BIGM (1E8)
#define BIGMM (1E16)

/* includes for parallel computation if necessary */


/* SYSTEM INCLUDES FOR PRINTING ---------------------------------------- */




/* LINEAR ALGEBRA LIBRARY ---------------------------------------------- */
/*
 * Initializes a vector of length 38 with a value.
 */
void boxesMPC_LA_INITIALIZEVECTOR_38(boxesMPC_FLOAT* vec, boxesMPC_FLOAT value)
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
void boxesMPC_LA_INITIALIZEVECTOR_20(boxesMPC_FLOAT* vec, boxesMPC_FLOAT value)
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
void boxesMPC_LA_INITIALIZEVECTOR_76(boxesMPC_FLOAT* vec, boxesMPC_FLOAT value)
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
void boxesMPC_LA_DOTACC_76(boxesMPC_FLOAT *x, boxesMPC_FLOAT *y, boxesMPC_FLOAT *z)
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
void boxesMPC_LA_DIAG_QUADFCN_4(boxesMPC_FLOAT* H, boxesMPC_FLOAT* f, boxesMPC_FLOAT* z, boxesMPC_FLOAT* grad, boxesMPC_FLOAT* value)
{
	int i;
	boxesMPC_FLOAT hz;	
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
void boxesMPC_LA_DIAG_QUADFCN_2(boxesMPC_FLOAT* H, boxesMPC_FLOAT* f, boxesMPC_FLOAT* z, boxesMPC_FLOAT* grad, boxesMPC_FLOAT* value)
{
	int i;
	boxesMPC_FLOAT hz;	
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
void boxesMPC_LA_DIAGZERO_MVMSUB6_2(boxesMPC_FLOAT *B, boxesMPC_FLOAT *u, boxesMPC_FLOAT *b, boxesMPC_FLOAT *l, boxesMPC_FLOAT *r, boxesMPC_FLOAT *z, boxesMPC_FLOAT *y)
{
	int i;
	boxesMPC_FLOAT Bu[2];
	boxesMPC_FLOAT norm = *y;
	boxesMPC_FLOAT lr = 0;

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
void boxesMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_4_4(boxesMPC_FLOAT *A, boxesMPC_FLOAT *x, boxesMPC_FLOAT *B, boxesMPC_FLOAT *u, boxesMPC_FLOAT *b, boxesMPC_FLOAT *l, boxesMPC_FLOAT *r, boxesMPC_FLOAT *z, boxesMPC_FLOAT *y)
{
	int i;
	int j;
	int k = 0;
	boxesMPC_FLOAT AxBu[2];
	boxesMPC_FLOAT norm = *y;
	boxesMPC_FLOAT lr = 0;

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
void boxesMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_4_2(boxesMPC_FLOAT *A, boxesMPC_FLOAT *x, boxesMPC_FLOAT *B, boxesMPC_FLOAT *u, boxesMPC_FLOAT *b, boxesMPC_FLOAT *l, boxesMPC_FLOAT *r, boxesMPC_FLOAT *z, boxesMPC_FLOAT *y)
{
	int i;
	int j;
	int k = 0;
	boxesMPC_FLOAT AxBu[2];
	boxesMPC_FLOAT norm = *y;
	boxesMPC_FLOAT lr = 0;

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
void boxesMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(boxesMPC_FLOAT *A, boxesMPC_FLOAT *x, boxesMPC_FLOAT *B, boxesMPC_FLOAT *y, boxesMPC_FLOAT *z)
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
void boxesMPC_LA_DIAGZERO_MTVM_2_2(boxesMPC_FLOAT *M, boxesMPC_FLOAT *x, boxesMPC_FLOAT *y)
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
void boxesMPC_LA_VSUBADD3_4(boxesMPC_FLOAT* t, boxesMPC_FLOAT* u, int* uidx, boxesMPC_FLOAT* v, boxesMPC_FLOAT* w, boxesMPC_FLOAT* y, boxesMPC_FLOAT* z, boxesMPC_FLOAT* r)
{
	int i;
	boxesMPC_FLOAT norm = *r;
	boxesMPC_FLOAT vx = 0;
	boxesMPC_FLOAT x;
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
void boxesMPC_LA_VSUBADD2_4(boxesMPC_FLOAT* t, int* tidx, boxesMPC_FLOAT* u, boxesMPC_FLOAT* v, boxesMPC_FLOAT* w, boxesMPC_FLOAT* y, boxesMPC_FLOAT* z, boxesMPC_FLOAT* r)
{
	int i;
	boxesMPC_FLOAT norm = *r;
	boxesMPC_FLOAT vx = 0;
	boxesMPC_FLOAT x;
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
void boxesMPC_LA_VSUBADD3_2(boxesMPC_FLOAT* t, boxesMPC_FLOAT* u, int* uidx, boxesMPC_FLOAT* v, boxesMPC_FLOAT* w, boxesMPC_FLOAT* y, boxesMPC_FLOAT* z, boxesMPC_FLOAT* r)
{
	int i;
	boxesMPC_FLOAT norm = *r;
	boxesMPC_FLOAT vx = 0;
	boxesMPC_FLOAT x;
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
void boxesMPC_LA_VSUBADD2_2(boxesMPC_FLOAT* t, int* tidx, boxesMPC_FLOAT* u, boxesMPC_FLOAT* v, boxesMPC_FLOAT* w, boxesMPC_FLOAT* y, boxesMPC_FLOAT* z, boxesMPC_FLOAT* r)
{
	int i;
	boxesMPC_FLOAT norm = *r;
	boxesMPC_FLOAT vx = 0;
	boxesMPC_FLOAT x;
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
void boxesMPC_LA_INEQ_B_GRAD_4_4_4(boxesMPC_FLOAT *lu, boxesMPC_FLOAT *su, boxesMPC_FLOAT *ru, boxesMPC_FLOAT *ll, boxesMPC_FLOAT *sl, boxesMPC_FLOAT *rl, int* lbIdx, int* ubIdx, boxesMPC_FLOAT *grad, boxesMPC_FLOAT *lubysu, boxesMPC_FLOAT *llbysl)
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
void boxesMPC_LA_INEQ_B_GRAD_2_2_2(boxesMPC_FLOAT *lu, boxesMPC_FLOAT *su, boxesMPC_FLOAT *ru, boxesMPC_FLOAT *ll, boxesMPC_FLOAT *sl, boxesMPC_FLOAT *rl, int* lbIdx, int* ubIdx, boxesMPC_FLOAT *grad, boxesMPC_FLOAT *lubysu, boxesMPC_FLOAT *llbysl)
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
void boxesMPC_LA_VVADD3_38(boxesMPC_FLOAT *u, boxesMPC_FLOAT *v, boxesMPC_FLOAT *w, boxesMPC_FLOAT *z)
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
void boxesMPC_LA_DIAG_CHOL_ONELOOP_LBUB_4_4_4(boxesMPC_FLOAT *H, boxesMPC_FLOAT *llbysl, int* lbIdx, boxesMPC_FLOAT *lubysu, int* ubIdx, boxesMPC_FLOAT *Phi)


{
	int i;
	
	/* compute cholesky */
	for( i=0; i<4; i++ ){
		Phi[i] = H[i] + llbysl[i] + lubysu[i];

#if boxesMPC_SET_PRINTLEVEL > 0 && defined PRINTNUMERICALWARNINGS
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
void boxesMPC_LA_DIAG_MATRIXFORWARDSUB_2_4(boxesMPC_FLOAT *L, boxesMPC_FLOAT *B, boxesMPC_FLOAT *A)
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
void boxesMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_4(boxesMPC_FLOAT *L, boxesMPC_FLOAT *B, boxesMPC_FLOAT *A)
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
void boxesMPC_LA_DENSE_DIAGZERO_MMTM_2_4_2(boxesMPC_FLOAT *A, boxesMPC_FLOAT *B, boxesMPC_FLOAT *C)
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
void boxesMPC_LA_DIAG_FORWARDSUB_4(boxesMPC_FLOAT *L, boxesMPC_FLOAT *b, boxesMPC_FLOAT *y)
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
void boxesMPC_LA_DIAG_CHOL_ONELOOP_LBUB_2_2_2(boxesMPC_FLOAT *H, boxesMPC_FLOAT *llbysl, int* lbIdx, boxesMPC_FLOAT *lubysu, int* ubIdx, boxesMPC_FLOAT *Phi)


{
	int i;
	
	/* compute cholesky */
	for( i=0; i<2; i++ ){
		Phi[i] = H[i] + llbysl[i] + lubysu[i];

#if boxesMPC_SET_PRINTLEVEL > 0 && defined PRINTNUMERICALWARNINGS
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
void boxesMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_2(boxesMPC_FLOAT *L, boxesMPC_FLOAT *B, boxesMPC_FLOAT *A)
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
void boxesMPC_LA_DIAG_FORWARDSUB_2(boxesMPC_FLOAT *L, boxesMPC_FLOAT *b, boxesMPC_FLOAT *y)
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
void boxesMPC_LA_DIAGZERO_MMT_2(boxesMPC_FLOAT *B, boxesMPC_FLOAT *L)
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
void boxesMPC_LA_DIAGZERO_MVMSUB7_2(boxesMPC_FLOAT *B, boxesMPC_FLOAT *u, boxesMPC_FLOAT *b, boxesMPC_FLOAT *r)
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
void boxesMPC_LA_DENSE_DIAGZERO_MMT2_2_4_4(boxesMPC_FLOAT *A, boxesMPC_FLOAT *B, boxesMPC_FLOAT *L)
{
    int i, j, k, ii, di;
    boxesMPC_FLOAT ltemp;
    
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
void boxesMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_4_4(boxesMPC_FLOAT *A, boxesMPC_FLOAT *x, boxesMPC_FLOAT *B, boxesMPC_FLOAT *u, boxesMPC_FLOAT *b, boxesMPC_FLOAT *r)
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
void boxesMPC_LA_DENSE_DIAGZERO_MMT2_2_4_2(boxesMPC_FLOAT *A, boxesMPC_FLOAT *B, boxesMPC_FLOAT *L)
{
    int i, j, k, ii, di;
    boxesMPC_FLOAT ltemp;
    
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
void boxesMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_4_2(boxesMPC_FLOAT *A, boxesMPC_FLOAT *x, boxesMPC_FLOAT *B, boxesMPC_FLOAT *u, boxesMPC_FLOAT *b, boxesMPC_FLOAT *r)
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
void boxesMPC_LA_DENSE_CHOL_2(boxesMPC_FLOAT *A, boxesMPC_FLOAT *L)
{
    int i, j, k, di, dj;
	 int ii, jj;

    boxesMPC_FLOAT l;
    boxesMPC_FLOAT Mii;

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
        
#if boxesMPC_SET_PRINTLEVEL > 0 && defined PRINTNUMERICALWARNINGS
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
void boxesMPC_LA_DENSE_FORWARDSUB_2(boxesMPC_FLOAT *L, boxesMPC_FLOAT *b, boxesMPC_FLOAT *y)
{
    int i,j,ii,di;
    boxesMPC_FLOAT yel;
            
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
void boxesMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(boxesMPC_FLOAT *L, boxesMPC_FLOAT *B, boxesMPC_FLOAT *A)
{
    int i,j,k,ii,di;
    boxesMPC_FLOAT a;
    
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
void boxesMPC_LA_DENSE_MMTSUB_2_2(boxesMPC_FLOAT *A, boxesMPC_FLOAT *L)
{
    int i, j, k, ii, di;
    boxesMPC_FLOAT ltemp;
    
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
void boxesMPC_LA_DENSE_MVMSUB1_2_2(boxesMPC_FLOAT *A, boxesMPC_FLOAT *x, boxesMPC_FLOAT *b, boxesMPC_FLOAT *r)
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
void boxesMPC_LA_DENSE_BACKWARDSUB_2(boxesMPC_FLOAT *L, boxesMPC_FLOAT *y, boxesMPC_FLOAT *x)
{
    int i, ii, di, j, jj, dj;
    boxesMPC_FLOAT xel;    
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
void boxesMPC_LA_DENSE_MTVMSUB_2_2(boxesMPC_FLOAT *A, boxesMPC_FLOAT *x, boxesMPC_FLOAT *b, boxesMPC_FLOAT *r)
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
void boxesMPC_LA_VSUB2_38(boxesMPC_FLOAT *x, boxesMPC_FLOAT *y, boxesMPC_FLOAT *z)
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
void boxesMPC_LA_DIAG_FORWARDBACKWARDSUB_4(boxesMPC_FLOAT *L, boxesMPC_FLOAT *b, boxesMPC_FLOAT *x)
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
void boxesMPC_LA_DIAG_FORWARDBACKWARDSUB_2(boxesMPC_FLOAT *L, boxesMPC_FLOAT *b, boxesMPC_FLOAT *x)
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
void boxesMPC_LA_VSUB_INDEXED_4(boxesMPC_FLOAT *x, int* xidx, boxesMPC_FLOAT *y, boxesMPC_FLOAT *z)
{
	int i;
	for( i=0; i<4; i++){
		z[i] = x[xidx[i]] - y[i];
	}
}


/*
 * Vector subtraction x = -u.*v - w for vectors of length 4.
 */
void boxesMPC_LA_VSUB3_4(boxesMPC_FLOAT *u, boxesMPC_FLOAT *v, boxesMPC_FLOAT *w, boxesMPC_FLOAT *x)
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
void boxesMPC_LA_VSUB2_INDEXED_4(boxesMPC_FLOAT *x, boxesMPC_FLOAT *y, int* yidx, boxesMPC_FLOAT *z)
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
void boxesMPC_LA_VSUB_INDEXED_2(boxesMPC_FLOAT *x, int* xidx, boxesMPC_FLOAT *y, boxesMPC_FLOAT *z)
{
	int i;
	for( i=0; i<2; i++){
		z[i] = x[xidx[i]] - y[i];
	}
}


/*
 * Vector subtraction x = -u.*v - w for vectors of length 2.
 */
void boxesMPC_LA_VSUB3_2(boxesMPC_FLOAT *u, boxesMPC_FLOAT *v, boxesMPC_FLOAT *w, boxesMPC_FLOAT *x)
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
void boxesMPC_LA_VSUB2_INDEXED_2(boxesMPC_FLOAT *x, boxesMPC_FLOAT *y, int* yidx, boxesMPC_FLOAT *z)
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
 * boxesMPC_NOPROGRESS (should be negative).
 */
int boxesMPC_LINESEARCH_BACKTRACKING_AFFINE(boxesMPC_FLOAT *l, boxesMPC_FLOAT *s, boxesMPC_FLOAT *dl, boxesMPC_FLOAT *ds, boxesMPC_FLOAT *a, boxesMPC_FLOAT *mu_aff)
{
    int i;
	int lsIt=1;    
    boxesMPC_FLOAT dltemp;
    boxesMPC_FLOAT dstemp;
    boxesMPC_FLOAT mya = 1.0;
    boxesMPC_FLOAT mymu;
        
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
            mya *= boxesMPC_SET_LS_SCALE_AFF;
            if( mya < boxesMPC_SET_LS_MINSTEP ){
                return boxesMPC_NOPROGRESS;
            }
        }
    }
    
    /* return new values and iteration counter */
    *a = mya;
    *mu_aff = mymu / (boxesMPC_FLOAT)76;
    return lsIt;
}


/*
 * Vector subtraction x = u.*v - a where a is a scalar
*  and x,u,v are vectors of length 76.
 */
void boxesMPC_LA_VSUB5_76(boxesMPC_FLOAT *u, boxesMPC_FLOAT *v, boxesMPC_FLOAT a, boxesMPC_FLOAT *x)
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
void boxesMPC_LA_VSUB6_INDEXED_4_4_4(boxesMPC_FLOAT *u, boxesMPC_FLOAT *su, int* uidx, boxesMPC_FLOAT *v, boxesMPC_FLOAT *sv, int* vidx, boxesMPC_FLOAT *x)
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
void boxesMPC_LA_DIAGZERO_MVM_2(boxesMPC_FLOAT *B, boxesMPC_FLOAT *u, boxesMPC_FLOAT *r)
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
void boxesMPC_LA_DENSE_DIAGZERO_2MVMADD_2_4_4(boxesMPC_FLOAT *A, boxesMPC_FLOAT *x, boxesMPC_FLOAT *B, boxesMPC_FLOAT *u, boxesMPC_FLOAT *r)
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
void boxesMPC_LA_VSUB6_INDEXED_2_2_2(boxesMPC_FLOAT *u, boxesMPC_FLOAT *su, int* uidx, boxesMPC_FLOAT *v, boxesMPC_FLOAT *sv, int* vidx, boxesMPC_FLOAT *x)
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
void boxesMPC_LA_DENSE_DIAGZERO_2MVMADD_2_4_2(boxesMPC_FLOAT *A, boxesMPC_FLOAT *x, boxesMPC_FLOAT *B, boxesMPC_FLOAT *u, boxesMPC_FLOAT *r)
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
void boxesMPC_LA_VSUB_38(boxesMPC_FLOAT *x, boxesMPC_FLOAT *y, boxesMPC_FLOAT *z)
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
void boxesMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_4(boxesMPC_FLOAT *r, boxesMPC_FLOAT *s, boxesMPC_FLOAT *u, boxesMPC_FLOAT *y, int* yidx, boxesMPC_FLOAT *z)
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
void boxesMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_4(boxesMPC_FLOAT *r, boxesMPC_FLOAT *s, boxesMPC_FLOAT *u, boxesMPC_FLOAT *y, int* yidx, boxesMPC_FLOAT *z)
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
void boxesMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_2(boxesMPC_FLOAT *r, boxesMPC_FLOAT *s, boxesMPC_FLOAT *u, boxesMPC_FLOAT *y, int* yidx, boxesMPC_FLOAT *z)
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
void boxesMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_2(boxesMPC_FLOAT *r, boxesMPC_FLOAT *s, boxesMPC_FLOAT *u, boxesMPC_FLOAT *y, int* yidx, boxesMPC_FLOAT *z)
{
	int i;
	for( i=0; i<2; i++ ){
		z[i] = -r[i]/s[i] + u[i]*y[yidx[i]];
	}
}


/*
 * Computes ds = -l.\(r + s.*dl) for vectors of length 76.
 */
void boxesMPC_LA_VSUB7_76(boxesMPC_FLOAT *l, boxesMPC_FLOAT *r, boxesMPC_FLOAT *s, boxesMPC_FLOAT *dl, boxesMPC_FLOAT *ds)
{
	int i;
	for( i=0; i<76; i++){
		ds[i] = -(r[i] + s[i]*dl[i])/l[i];
	}
}


/*
 * Vector addition x = x + y for vectors of length 38.
 */
void boxesMPC_LA_VADD_38(boxesMPC_FLOAT *x, boxesMPC_FLOAT *y)
{
	int i;
	for( i=0; i<38; i++){
		x[i] += y[i];
	}
}


/*
 * Vector addition x = x + y for vectors of length 20.
 */
void boxesMPC_LA_VADD_20(boxesMPC_FLOAT *x, boxesMPC_FLOAT *y)
{
	int i;
	for( i=0; i<20; i++){
		x[i] += y[i];
	}
}


/*
 * Vector addition x = x + y for vectors of length 76.
 */
void boxesMPC_LA_VADD_76(boxesMPC_FLOAT *x, boxesMPC_FLOAT *y)
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
int boxesMPC_LINESEARCH_BACKTRACKING_COMBINED(boxesMPC_FLOAT *z, boxesMPC_FLOAT *v, boxesMPC_FLOAT *l, boxesMPC_FLOAT *s, boxesMPC_FLOAT *dz, boxesMPC_FLOAT *dv, boxesMPC_FLOAT *dl, boxesMPC_FLOAT *ds, boxesMPC_FLOAT *a, boxesMPC_FLOAT *mu)
{
    int i, lsIt=1;       
    boxesMPC_FLOAT dltemp;
    boxesMPC_FLOAT dstemp;    
    boxesMPC_FLOAT a_gamma;
            
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
            *a *= boxesMPC_SET_LS_SCALE;
            if( *a < boxesMPC_SET_LS_MINSTEP ){
                return boxesMPC_NOPROGRESS;
            }
        }
    }
    
    /* update variables with safety margin */
    a_gamma = (*a)*boxesMPC_SET_LS_MAXSTEP;
    
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
    *mu /= (boxesMPC_FLOAT)76;
    return lsIt;
}




/* VARIABLE DEFINITIONS ------------------------------------------------ */
boxesMPC_FLOAT boxesMPC_z[38];
boxesMPC_FLOAT boxesMPC_v[20];
boxesMPC_FLOAT boxesMPC_dz_aff[38];
boxesMPC_FLOAT boxesMPC_dv_aff[20];
boxesMPC_FLOAT boxesMPC_grad_cost[38];
boxesMPC_FLOAT boxesMPC_grad_eq[38];
boxesMPC_FLOAT boxesMPC_rd[38];
boxesMPC_FLOAT boxesMPC_l[76];
boxesMPC_FLOAT boxesMPC_s[76];
boxesMPC_FLOAT boxesMPC_lbys[76];
boxesMPC_FLOAT boxesMPC_dl_aff[76];
boxesMPC_FLOAT boxesMPC_ds_aff[76];
boxesMPC_FLOAT boxesMPC_dz_cc[38];
boxesMPC_FLOAT boxesMPC_dv_cc[20];
boxesMPC_FLOAT boxesMPC_dl_cc[76];
boxesMPC_FLOAT boxesMPC_ds_cc[76];
boxesMPC_FLOAT boxesMPC_ccrhs[76];
boxesMPC_FLOAT boxesMPC_grad_ineq[38];
boxesMPC_FLOAT* boxesMPC_z0 = boxesMPC_z + 0;
boxesMPC_FLOAT* boxesMPC_dzaff0 = boxesMPC_dz_aff + 0;
boxesMPC_FLOAT* boxesMPC_dzcc0 = boxesMPC_dz_cc + 0;
boxesMPC_FLOAT* boxesMPC_rd0 = boxesMPC_rd + 0;
boxesMPC_FLOAT boxesMPC_Lbyrd0[4];
boxesMPC_FLOAT* boxesMPC_grad_cost0 = boxesMPC_grad_cost + 0;
boxesMPC_FLOAT* boxesMPC_grad_eq0 = boxesMPC_grad_eq + 0;
boxesMPC_FLOAT* boxesMPC_grad_ineq0 = boxesMPC_grad_ineq + 0;
boxesMPC_FLOAT boxesMPC_ctv0[4];
boxesMPC_FLOAT boxesMPC_C0[8] = {1.0000000000000000E+000, 0.0000000000000000E+000, 
0.0000000000000000E+000, 1.0000000000000000E+000, 
1.0000000000000000E+000, 0.0000000000000000E+000, 
0.0000000000000000E+000, 1.0000000000000000E+000};
boxesMPC_FLOAT* boxesMPC_v0 = boxesMPC_v + 0;
boxesMPC_FLOAT boxesMPC_re0[2];
boxesMPC_FLOAT boxesMPC_beta0[2];
boxesMPC_FLOAT boxesMPC_betacc0[2];
boxesMPC_FLOAT* boxesMPC_dvaff0 = boxesMPC_dv_aff + 0;
boxesMPC_FLOAT* boxesMPC_dvcc0 = boxesMPC_dv_cc + 0;
boxesMPC_FLOAT boxesMPC_V0[8];
boxesMPC_FLOAT boxesMPC_Yd0[3];
boxesMPC_FLOAT boxesMPC_Ld0[3];
boxesMPC_FLOAT boxesMPC_yy0[2];
boxesMPC_FLOAT boxesMPC_bmy0[2];
int boxesMPC_lbIdx0[4] = {0, 1, 2, 3};
boxesMPC_FLOAT* boxesMPC_llb0 = boxesMPC_l + 0;
boxesMPC_FLOAT* boxesMPC_slb0 = boxesMPC_s + 0;
boxesMPC_FLOAT* boxesMPC_llbbyslb0 = boxesMPC_lbys + 0;
boxesMPC_FLOAT boxesMPC_rilb0[4];
boxesMPC_FLOAT* boxesMPC_dllbaff0 = boxesMPC_dl_aff + 0;
boxesMPC_FLOAT* boxesMPC_dslbaff0 = boxesMPC_ds_aff + 0;
boxesMPC_FLOAT* boxesMPC_dllbcc0 = boxesMPC_dl_cc + 0;
boxesMPC_FLOAT* boxesMPC_dslbcc0 = boxesMPC_ds_cc + 0;
boxesMPC_FLOAT* boxesMPC_ccrhsl0 = boxesMPC_ccrhs + 0;
int boxesMPC_ubIdx0[4] = {0, 1, 2, 3};
boxesMPC_FLOAT* boxesMPC_lub0 = boxesMPC_l + 4;
boxesMPC_FLOAT* boxesMPC_sub0 = boxesMPC_s + 4;
boxesMPC_FLOAT* boxesMPC_lubbysub0 = boxesMPC_lbys + 4;
boxesMPC_FLOAT boxesMPC_riub0[4];
boxesMPC_FLOAT* boxesMPC_dlubaff0 = boxesMPC_dl_aff + 4;
boxesMPC_FLOAT* boxesMPC_dsubaff0 = boxesMPC_ds_aff + 4;
boxesMPC_FLOAT* boxesMPC_dlubcc0 = boxesMPC_dl_cc + 4;
boxesMPC_FLOAT* boxesMPC_dsubcc0 = boxesMPC_ds_cc + 4;
boxesMPC_FLOAT* boxesMPC_ccrhsub0 = boxesMPC_ccrhs + 4;
boxesMPC_FLOAT boxesMPC_Phi0[4];
boxesMPC_FLOAT boxesMPC_D0[4] = {1.0000000000000000E+000, 
1.0000000000000000E+000};
boxesMPC_FLOAT boxesMPC_W0[4];
boxesMPC_FLOAT* boxesMPC_z1 = boxesMPC_z + 4;
boxesMPC_FLOAT* boxesMPC_dzaff1 = boxesMPC_dz_aff + 4;
boxesMPC_FLOAT* boxesMPC_dzcc1 = boxesMPC_dz_cc + 4;
boxesMPC_FLOAT* boxesMPC_rd1 = boxesMPC_rd + 4;
boxesMPC_FLOAT boxesMPC_Lbyrd1[4];
boxesMPC_FLOAT* boxesMPC_grad_cost1 = boxesMPC_grad_cost + 4;
boxesMPC_FLOAT* boxesMPC_grad_eq1 = boxesMPC_grad_eq + 4;
boxesMPC_FLOAT* boxesMPC_grad_ineq1 = boxesMPC_grad_ineq + 4;
boxesMPC_FLOAT boxesMPC_ctv1[4];
boxesMPC_FLOAT* boxesMPC_v1 = boxesMPC_v + 2;
boxesMPC_FLOAT boxesMPC_re1[2];
boxesMPC_FLOAT boxesMPC_beta1[2];
boxesMPC_FLOAT boxesMPC_betacc1[2];
boxesMPC_FLOAT* boxesMPC_dvaff1 = boxesMPC_dv_aff + 2;
boxesMPC_FLOAT* boxesMPC_dvcc1 = boxesMPC_dv_cc + 2;
boxesMPC_FLOAT boxesMPC_V1[8];
boxesMPC_FLOAT boxesMPC_Yd1[3];
boxesMPC_FLOAT boxesMPC_Ld1[3];
boxesMPC_FLOAT boxesMPC_yy1[2];
boxesMPC_FLOAT boxesMPC_bmy1[2];
boxesMPC_FLOAT boxesMPC_c1[2] = {0.0000000000000000E+000, 0.0000000000000000E+000};
int boxesMPC_lbIdx1[4] = {0, 1, 2, 3};
boxesMPC_FLOAT* boxesMPC_llb1 = boxesMPC_l + 8;
boxesMPC_FLOAT* boxesMPC_slb1 = boxesMPC_s + 8;
boxesMPC_FLOAT* boxesMPC_llbbyslb1 = boxesMPC_lbys + 8;
boxesMPC_FLOAT boxesMPC_rilb1[4];
boxesMPC_FLOAT* boxesMPC_dllbaff1 = boxesMPC_dl_aff + 8;
boxesMPC_FLOAT* boxesMPC_dslbaff1 = boxesMPC_ds_aff + 8;
boxesMPC_FLOAT* boxesMPC_dllbcc1 = boxesMPC_dl_cc + 8;
boxesMPC_FLOAT* boxesMPC_dslbcc1 = boxesMPC_ds_cc + 8;
boxesMPC_FLOAT* boxesMPC_ccrhsl1 = boxesMPC_ccrhs + 8;
int boxesMPC_ubIdx1[4] = {0, 1, 2, 3};
boxesMPC_FLOAT* boxesMPC_lub1 = boxesMPC_l + 12;
boxesMPC_FLOAT* boxesMPC_sub1 = boxesMPC_s + 12;
boxesMPC_FLOAT* boxesMPC_lubbysub1 = boxesMPC_lbys + 12;
boxesMPC_FLOAT boxesMPC_riub1[4];
boxesMPC_FLOAT* boxesMPC_dlubaff1 = boxesMPC_dl_aff + 12;
boxesMPC_FLOAT* boxesMPC_dsubaff1 = boxesMPC_ds_aff + 12;
boxesMPC_FLOAT* boxesMPC_dlubcc1 = boxesMPC_dl_cc + 12;
boxesMPC_FLOAT* boxesMPC_dsubcc1 = boxesMPC_ds_cc + 12;
boxesMPC_FLOAT* boxesMPC_ccrhsub1 = boxesMPC_ccrhs + 12;
boxesMPC_FLOAT boxesMPC_Phi1[4];
boxesMPC_FLOAT boxesMPC_D1[4] = {-1.0000000000000000E+000, 
-1.0000000000000000E+000};
boxesMPC_FLOAT boxesMPC_W1[4];
boxesMPC_FLOAT boxesMPC_Ysd1[4];
boxesMPC_FLOAT boxesMPC_Lsd1[4];
boxesMPC_FLOAT* boxesMPC_z2 = boxesMPC_z + 8;
boxesMPC_FLOAT* boxesMPC_dzaff2 = boxesMPC_dz_aff + 8;
boxesMPC_FLOAT* boxesMPC_dzcc2 = boxesMPC_dz_cc + 8;
boxesMPC_FLOAT* boxesMPC_rd2 = boxesMPC_rd + 8;
boxesMPC_FLOAT boxesMPC_Lbyrd2[4];
boxesMPC_FLOAT* boxesMPC_grad_cost2 = boxesMPC_grad_cost + 8;
boxesMPC_FLOAT* boxesMPC_grad_eq2 = boxesMPC_grad_eq + 8;
boxesMPC_FLOAT* boxesMPC_grad_ineq2 = boxesMPC_grad_ineq + 8;
boxesMPC_FLOAT boxesMPC_ctv2[4];
boxesMPC_FLOAT* boxesMPC_v2 = boxesMPC_v + 4;
boxesMPC_FLOAT boxesMPC_re2[2];
boxesMPC_FLOAT boxesMPC_beta2[2];
boxesMPC_FLOAT boxesMPC_betacc2[2];
boxesMPC_FLOAT* boxesMPC_dvaff2 = boxesMPC_dv_aff + 4;
boxesMPC_FLOAT* boxesMPC_dvcc2 = boxesMPC_dv_cc + 4;
boxesMPC_FLOAT boxesMPC_V2[8];
boxesMPC_FLOAT boxesMPC_Yd2[3];
boxesMPC_FLOAT boxesMPC_Ld2[3];
boxesMPC_FLOAT boxesMPC_yy2[2];
boxesMPC_FLOAT boxesMPC_bmy2[2];
boxesMPC_FLOAT boxesMPC_c2[2] = {0.0000000000000000E+000, 0.0000000000000000E+000};
int boxesMPC_lbIdx2[4] = {0, 1, 2, 3};
boxesMPC_FLOAT* boxesMPC_llb2 = boxesMPC_l + 16;
boxesMPC_FLOAT* boxesMPC_slb2 = boxesMPC_s + 16;
boxesMPC_FLOAT* boxesMPC_llbbyslb2 = boxesMPC_lbys + 16;
boxesMPC_FLOAT boxesMPC_rilb2[4];
boxesMPC_FLOAT* boxesMPC_dllbaff2 = boxesMPC_dl_aff + 16;
boxesMPC_FLOAT* boxesMPC_dslbaff2 = boxesMPC_ds_aff + 16;
boxesMPC_FLOAT* boxesMPC_dllbcc2 = boxesMPC_dl_cc + 16;
boxesMPC_FLOAT* boxesMPC_dslbcc2 = boxesMPC_ds_cc + 16;
boxesMPC_FLOAT* boxesMPC_ccrhsl2 = boxesMPC_ccrhs + 16;
int boxesMPC_ubIdx2[4] = {0, 1, 2, 3};
boxesMPC_FLOAT* boxesMPC_lub2 = boxesMPC_l + 20;
boxesMPC_FLOAT* boxesMPC_sub2 = boxesMPC_s + 20;
boxesMPC_FLOAT* boxesMPC_lubbysub2 = boxesMPC_lbys + 20;
boxesMPC_FLOAT boxesMPC_riub2[4];
boxesMPC_FLOAT* boxesMPC_dlubaff2 = boxesMPC_dl_aff + 20;
boxesMPC_FLOAT* boxesMPC_dsubaff2 = boxesMPC_ds_aff + 20;
boxesMPC_FLOAT* boxesMPC_dlubcc2 = boxesMPC_dl_cc + 20;
boxesMPC_FLOAT* boxesMPC_dsubcc2 = boxesMPC_ds_cc + 20;
boxesMPC_FLOAT* boxesMPC_ccrhsub2 = boxesMPC_ccrhs + 20;
boxesMPC_FLOAT boxesMPC_Phi2[4];
boxesMPC_FLOAT boxesMPC_W2[4];
boxesMPC_FLOAT boxesMPC_Ysd2[4];
boxesMPC_FLOAT boxesMPC_Lsd2[4];
boxesMPC_FLOAT* boxesMPC_z3 = boxesMPC_z + 12;
boxesMPC_FLOAT* boxesMPC_dzaff3 = boxesMPC_dz_aff + 12;
boxesMPC_FLOAT* boxesMPC_dzcc3 = boxesMPC_dz_cc + 12;
boxesMPC_FLOAT* boxesMPC_rd3 = boxesMPC_rd + 12;
boxesMPC_FLOAT boxesMPC_Lbyrd3[4];
boxesMPC_FLOAT* boxesMPC_grad_cost3 = boxesMPC_grad_cost + 12;
boxesMPC_FLOAT* boxesMPC_grad_eq3 = boxesMPC_grad_eq + 12;
boxesMPC_FLOAT* boxesMPC_grad_ineq3 = boxesMPC_grad_ineq + 12;
boxesMPC_FLOAT boxesMPC_ctv3[4];
boxesMPC_FLOAT* boxesMPC_v3 = boxesMPC_v + 6;
boxesMPC_FLOAT boxesMPC_re3[2];
boxesMPC_FLOAT boxesMPC_beta3[2];
boxesMPC_FLOAT boxesMPC_betacc3[2];
boxesMPC_FLOAT* boxesMPC_dvaff3 = boxesMPC_dv_aff + 6;
boxesMPC_FLOAT* boxesMPC_dvcc3 = boxesMPC_dv_cc + 6;
boxesMPC_FLOAT boxesMPC_V3[8];
boxesMPC_FLOAT boxesMPC_Yd3[3];
boxesMPC_FLOAT boxesMPC_Ld3[3];
boxesMPC_FLOAT boxesMPC_yy3[2];
boxesMPC_FLOAT boxesMPC_bmy3[2];
boxesMPC_FLOAT boxesMPC_c3[2] = {0.0000000000000000E+000, 0.0000000000000000E+000};
int boxesMPC_lbIdx3[4] = {0, 1, 2, 3};
boxesMPC_FLOAT* boxesMPC_llb3 = boxesMPC_l + 24;
boxesMPC_FLOAT* boxesMPC_slb3 = boxesMPC_s + 24;
boxesMPC_FLOAT* boxesMPC_llbbyslb3 = boxesMPC_lbys + 24;
boxesMPC_FLOAT boxesMPC_rilb3[4];
boxesMPC_FLOAT* boxesMPC_dllbaff3 = boxesMPC_dl_aff + 24;
boxesMPC_FLOAT* boxesMPC_dslbaff3 = boxesMPC_ds_aff + 24;
boxesMPC_FLOAT* boxesMPC_dllbcc3 = boxesMPC_dl_cc + 24;
boxesMPC_FLOAT* boxesMPC_dslbcc3 = boxesMPC_ds_cc + 24;
boxesMPC_FLOAT* boxesMPC_ccrhsl3 = boxesMPC_ccrhs + 24;
int boxesMPC_ubIdx3[4] = {0, 1, 2, 3};
boxesMPC_FLOAT* boxesMPC_lub3 = boxesMPC_l + 28;
boxesMPC_FLOAT* boxesMPC_sub3 = boxesMPC_s + 28;
boxesMPC_FLOAT* boxesMPC_lubbysub3 = boxesMPC_lbys + 28;
boxesMPC_FLOAT boxesMPC_riub3[4];
boxesMPC_FLOAT* boxesMPC_dlubaff3 = boxesMPC_dl_aff + 28;
boxesMPC_FLOAT* boxesMPC_dsubaff3 = boxesMPC_ds_aff + 28;
boxesMPC_FLOAT* boxesMPC_dlubcc3 = boxesMPC_dl_cc + 28;
boxesMPC_FLOAT* boxesMPC_dsubcc3 = boxesMPC_ds_cc + 28;
boxesMPC_FLOAT* boxesMPC_ccrhsub3 = boxesMPC_ccrhs + 28;
boxesMPC_FLOAT boxesMPC_Phi3[4];
boxesMPC_FLOAT boxesMPC_W3[4];
boxesMPC_FLOAT boxesMPC_Ysd3[4];
boxesMPC_FLOAT boxesMPC_Lsd3[4];
boxesMPC_FLOAT* boxesMPC_z4 = boxesMPC_z + 16;
boxesMPC_FLOAT* boxesMPC_dzaff4 = boxesMPC_dz_aff + 16;
boxesMPC_FLOAT* boxesMPC_dzcc4 = boxesMPC_dz_cc + 16;
boxesMPC_FLOAT* boxesMPC_rd4 = boxesMPC_rd + 16;
boxesMPC_FLOAT boxesMPC_Lbyrd4[4];
boxesMPC_FLOAT* boxesMPC_grad_cost4 = boxesMPC_grad_cost + 16;
boxesMPC_FLOAT* boxesMPC_grad_eq4 = boxesMPC_grad_eq + 16;
boxesMPC_FLOAT* boxesMPC_grad_ineq4 = boxesMPC_grad_ineq + 16;
boxesMPC_FLOAT boxesMPC_ctv4[4];
boxesMPC_FLOAT* boxesMPC_v4 = boxesMPC_v + 8;
boxesMPC_FLOAT boxesMPC_re4[2];
boxesMPC_FLOAT boxesMPC_beta4[2];
boxesMPC_FLOAT boxesMPC_betacc4[2];
boxesMPC_FLOAT* boxesMPC_dvaff4 = boxesMPC_dv_aff + 8;
boxesMPC_FLOAT* boxesMPC_dvcc4 = boxesMPC_dv_cc + 8;
boxesMPC_FLOAT boxesMPC_V4[8];
boxesMPC_FLOAT boxesMPC_Yd4[3];
boxesMPC_FLOAT boxesMPC_Ld4[3];
boxesMPC_FLOAT boxesMPC_yy4[2];
boxesMPC_FLOAT boxesMPC_bmy4[2];
boxesMPC_FLOAT boxesMPC_c4[2] = {0.0000000000000000E+000, 0.0000000000000000E+000};
int boxesMPC_lbIdx4[4] = {0, 1, 2, 3};
boxesMPC_FLOAT* boxesMPC_llb4 = boxesMPC_l + 32;
boxesMPC_FLOAT* boxesMPC_slb4 = boxesMPC_s + 32;
boxesMPC_FLOAT* boxesMPC_llbbyslb4 = boxesMPC_lbys + 32;
boxesMPC_FLOAT boxesMPC_rilb4[4];
boxesMPC_FLOAT* boxesMPC_dllbaff4 = boxesMPC_dl_aff + 32;
boxesMPC_FLOAT* boxesMPC_dslbaff4 = boxesMPC_ds_aff + 32;
boxesMPC_FLOAT* boxesMPC_dllbcc4 = boxesMPC_dl_cc + 32;
boxesMPC_FLOAT* boxesMPC_dslbcc4 = boxesMPC_ds_cc + 32;
boxesMPC_FLOAT* boxesMPC_ccrhsl4 = boxesMPC_ccrhs + 32;
int boxesMPC_ubIdx4[4] = {0, 1, 2, 3};
boxesMPC_FLOAT* boxesMPC_lub4 = boxesMPC_l + 36;
boxesMPC_FLOAT* boxesMPC_sub4 = boxesMPC_s + 36;
boxesMPC_FLOAT* boxesMPC_lubbysub4 = boxesMPC_lbys + 36;
boxesMPC_FLOAT boxesMPC_riub4[4];
boxesMPC_FLOAT* boxesMPC_dlubaff4 = boxesMPC_dl_aff + 36;
boxesMPC_FLOAT* boxesMPC_dsubaff4 = boxesMPC_ds_aff + 36;
boxesMPC_FLOAT* boxesMPC_dlubcc4 = boxesMPC_dl_cc + 36;
boxesMPC_FLOAT* boxesMPC_dsubcc4 = boxesMPC_ds_cc + 36;
boxesMPC_FLOAT* boxesMPC_ccrhsub4 = boxesMPC_ccrhs + 36;
boxesMPC_FLOAT boxesMPC_Phi4[4];
boxesMPC_FLOAT boxesMPC_W4[4];
boxesMPC_FLOAT boxesMPC_Ysd4[4];
boxesMPC_FLOAT boxesMPC_Lsd4[4];
boxesMPC_FLOAT* boxesMPC_z5 = boxesMPC_z + 20;
boxesMPC_FLOAT* boxesMPC_dzaff5 = boxesMPC_dz_aff + 20;
boxesMPC_FLOAT* boxesMPC_dzcc5 = boxesMPC_dz_cc + 20;
boxesMPC_FLOAT* boxesMPC_rd5 = boxesMPC_rd + 20;
boxesMPC_FLOAT boxesMPC_Lbyrd5[4];
boxesMPC_FLOAT* boxesMPC_grad_cost5 = boxesMPC_grad_cost + 20;
boxesMPC_FLOAT* boxesMPC_grad_eq5 = boxesMPC_grad_eq + 20;
boxesMPC_FLOAT* boxesMPC_grad_ineq5 = boxesMPC_grad_ineq + 20;
boxesMPC_FLOAT boxesMPC_ctv5[4];
boxesMPC_FLOAT* boxesMPC_v5 = boxesMPC_v + 10;
boxesMPC_FLOAT boxesMPC_re5[2];
boxesMPC_FLOAT boxesMPC_beta5[2];
boxesMPC_FLOAT boxesMPC_betacc5[2];
boxesMPC_FLOAT* boxesMPC_dvaff5 = boxesMPC_dv_aff + 10;
boxesMPC_FLOAT* boxesMPC_dvcc5 = boxesMPC_dv_cc + 10;
boxesMPC_FLOAT boxesMPC_V5[8];
boxesMPC_FLOAT boxesMPC_Yd5[3];
boxesMPC_FLOAT boxesMPC_Ld5[3];
boxesMPC_FLOAT boxesMPC_yy5[2];
boxesMPC_FLOAT boxesMPC_bmy5[2];
boxesMPC_FLOAT boxesMPC_c5[2] = {0.0000000000000000E+000, 0.0000000000000000E+000};
int boxesMPC_lbIdx5[4] = {0, 1, 2, 3};
boxesMPC_FLOAT* boxesMPC_llb5 = boxesMPC_l + 40;
boxesMPC_FLOAT* boxesMPC_slb5 = boxesMPC_s + 40;
boxesMPC_FLOAT* boxesMPC_llbbyslb5 = boxesMPC_lbys + 40;
boxesMPC_FLOAT boxesMPC_rilb5[4];
boxesMPC_FLOAT* boxesMPC_dllbaff5 = boxesMPC_dl_aff + 40;
boxesMPC_FLOAT* boxesMPC_dslbaff5 = boxesMPC_ds_aff + 40;
boxesMPC_FLOAT* boxesMPC_dllbcc5 = boxesMPC_dl_cc + 40;
boxesMPC_FLOAT* boxesMPC_dslbcc5 = boxesMPC_ds_cc + 40;
boxesMPC_FLOAT* boxesMPC_ccrhsl5 = boxesMPC_ccrhs + 40;
int boxesMPC_ubIdx5[4] = {0, 1, 2, 3};
boxesMPC_FLOAT* boxesMPC_lub5 = boxesMPC_l + 44;
boxesMPC_FLOAT* boxesMPC_sub5 = boxesMPC_s + 44;
boxesMPC_FLOAT* boxesMPC_lubbysub5 = boxesMPC_lbys + 44;
boxesMPC_FLOAT boxesMPC_riub5[4];
boxesMPC_FLOAT* boxesMPC_dlubaff5 = boxesMPC_dl_aff + 44;
boxesMPC_FLOAT* boxesMPC_dsubaff5 = boxesMPC_ds_aff + 44;
boxesMPC_FLOAT* boxesMPC_dlubcc5 = boxesMPC_dl_cc + 44;
boxesMPC_FLOAT* boxesMPC_dsubcc5 = boxesMPC_ds_cc + 44;
boxesMPC_FLOAT* boxesMPC_ccrhsub5 = boxesMPC_ccrhs + 44;
boxesMPC_FLOAT boxesMPC_Phi5[4];
boxesMPC_FLOAT boxesMPC_W5[4];
boxesMPC_FLOAT boxesMPC_Ysd5[4];
boxesMPC_FLOAT boxesMPC_Lsd5[4];
boxesMPC_FLOAT* boxesMPC_z6 = boxesMPC_z + 24;
boxesMPC_FLOAT* boxesMPC_dzaff6 = boxesMPC_dz_aff + 24;
boxesMPC_FLOAT* boxesMPC_dzcc6 = boxesMPC_dz_cc + 24;
boxesMPC_FLOAT* boxesMPC_rd6 = boxesMPC_rd + 24;
boxesMPC_FLOAT boxesMPC_Lbyrd6[4];
boxesMPC_FLOAT* boxesMPC_grad_cost6 = boxesMPC_grad_cost + 24;
boxesMPC_FLOAT* boxesMPC_grad_eq6 = boxesMPC_grad_eq + 24;
boxesMPC_FLOAT* boxesMPC_grad_ineq6 = boxesMPC_grad_ineq + 24;
boxesMPC_FLOAT boxesMPC_ctv6[4];
boxesMPC_FLOAT* boxesMPC_v6 = boxesMPC_v + 12;
boxesMPC_FLOAT boxesMPC_re6[2];
boxesMPC_FLOAT boxesMPC_beta6[2];
boxesMPC_FLOAT boxesMPC_betacc6[2];
boxesMPC_FLOAT* boxesMPC_dvaff6 = boxesMPC_dv_aff + 12;
boxesMPC_FLOAT* boxesMPC_dvcc6 = boxesMPC_dv_cc + 12;
boxesMPC_FLOAT boxesMPC_V6[8];
boxesMPC_FLOAT boxesMPC_Yd6[3];
boxesMPC_FLOAT boxesMPC_Ld6[3];
boxesMPC_FLOAT boxesMPC_yy6[2];
boxesMPC_FLOAT boxesMPC_bmy6[2];
boxesMPC_FLOAT boxesMPC_c6[2] = {0.0000000000000000E+000, 0.0000000000000000E+000};
int boxesMPC_lbIdx6[4] = {0, 1, 2, 3};
boxesMPC_FLOAT* boxesMPC_llb6 = boxesMPC_l + 48;
boxesMPC_FLOAT* boxesMPC_slb6 = boxesMPC_s + 48;
boxesMPC_FLOAT* boxesMPC_llbbyslb6 = boxesMPC_lbys + 48;
boxesMPC_FLOAT boxesMPC_rilb6[4];
boxesMPC_FLOAT* boxesMPC_dllbaff6 = boxesMPC_dl_aff + 48;
boxesMPC_FLOAT* boxesMPC_dslbaff6 = boxesMPC_ds_aff + 48;
boxesMPC_FLOAT* boxesMPC_dllbcc6 = boxesMPC_dl_cc + 48;
boxesMPC_FLOAT* boxesMPC_dslbcc6 = boxesMPC_ds_cc + 48;
boxesMPC_FLOAT* boxesMPC_ccrhsl6 = boxesMPC_ccrhs + 48;
int boxesMPC_ubIdx6[4] = {0, 1, 2, 3};
boxesMPC_FLOAT* boxesMPC_lub6 = boxesMPC_l + 52;
boxesMPC_FLOAT* boxesMPC_sub6 = boxesMPC_s + 52;
boxesMPC_FLOAT* boxesMPC_lubbysub6 = boxesMPC_lbys + 52;
boxesMPC_FLOAT boxesMPC_riub6[4];
boxesMPC_FLOAT* boxesMPC_dlubaff6 = boxesMPC_dl_aff + 52;
boxesMPC_FLOAT* boxesMPC_dsubaff6 = boxesMPC_ds_aff + 52;
boxesMPC_FLOAT* boxesMPC_dlubcc6 = boxesMPC_dl_cc + 52;
boxesMPC_FLOAT* boxesMPC_dsubcc6 = boxesMPC_ds_cc + 52;
boxesMPC_FLOAT* boxesMPC_ccrhsub6 = boxesMPC_ccrhs + 52;
boxesMPC_FLOAT boxesMPC_Phi6[4];
boxesMPC_FLOAT boxesMPC_W6[4];
boxesMPC_FLOAT boxesMPC_Ysd6[4];
boxesMPC_FLOAT boxesMPC_Lsd6[4];
boxesMPC_FLOAT* boxesMPC_z7 = boxesMPC_z + 28;
boxesMPC_FLOAT* boxesMPC_dzaff7 = boxesMPC_dz_aff + 28;
boxesMPC_FLOAT* boxesMPC_dzcc7 = boxesMPC_dz_cc + 28;
boxesMPC_FLOAT* boxesMPC_rd7 = boxesMPC_rd + 28;
boxesMPC_FLOAT boxesMPC_Lbyrd7[4];
boxesMPC_FLOAT* boxesMPC_grad_cost7 = boxesMPC_grad_cost + 28;
boxesMPC_FLOAT* boxesMPC_grad_eq7 = boxesMPC_grad_eq + 28;
boxesMPC_FLOAT* boxesMPC_grad_ineq7 = boxesMPC_grad_ineq + 28;
boxesMPC_FLOAT boxesMPC_ctv7[4];
boxesMPC_FLOAT* boxesMPC_v7 = boxesMPC_v + 14;
boxesMPC_FLOAT boxesMPC_re7[2];
boxesMPC_FLOAT boxesMPC_beta7[2];
boxesMPC_FLOAT boxesMPC_betacc7[2];
boxesMPC_FLOAT* boxesMPC_dvaff7 = boxesMPC_dv_aff + 14;
boxesMPC_FLOAT* boxesMPC_dvcc7 = boxesMPC_dv_cc + 14;
boxesMPC_FLOAT boxesMPC_V7[8];
boxesMPC_FLOAT boxesMPC_Yd7[3];
boxesMPC_FLOAT boxesMPC_Ld7[3];
boxesMPC_FLOAT boxesMPC_yy7[2];
boxesMPC_FLOAT boxesMPC_bmy7[2];
boxesMPC_FLOAT boxesMPC_c7[2] = {0.0000000000000000E+000, 0.0000000000000000E+000};
int boxesMPC_lbIdx7[4] = {0, 1, 2, 3};
boxesMPC_FLOAT* boxesMPC_llb7 = boxesMPC_l + 56;
boxesMPC_FLOAT* boxesMPC_slb7 = boxesMPC_s + 56;
boxesMPC_FLOAT* boxesMPC_llbbyslb7 = boxesMPC_lbys + 56;
boxesMPC_FLOAT boxesMPC_rilb7[4];
boxesMPC_FLOAT* boxesMPC_dllbaff7 = boxesMPC_dl_aff + 56;
boxesMPC_FLOAT* boxesMPC_dslbaff7 = boxesMPC_ds_aff + 56;
boxesMPC_FLOAT* boxesMPC_dllbcc7 = boxesMPC_dl_cc + 56;
boxesMPC_FLOAT* boxesMPC_dslbcc7 = boxesMPC_ds_cc + 56;
boxesMPC_FLOAT* boxesMPC_ccrhsl7 = boxesMPC_ccrhs + 56;
int boxesMPC_ubIdx7[4] = {0, 1, 2, 3};
boxesMPC_FLOAT* boxesMPC_lub7 = boxesMPC_l + 60;
boxesMPC_FLOAT* boxesMPC_sub7 = boxesMPC_s + 60;
boxesMPC_FLOAT* boxesMPC_lubbysub7 = boxesMPC_lbys + 60;
boxesMPC_FLOAT boxesMPC_riub7[4];
boxesMPC_FLOAT* boxesMPC_dlubaff7 = boxesMPC_dl_aff + 60;
boxesMPC_FLOAT* boxesMPC_dsubaff7 = boxesMPC_ds_aff + 60;
boxesMPC_FLOAT* boxesMPC_dlubcc7 = boxesMPC_dl_cc + 60;
boxesMPC_FLOAT* boxesMPC_dsubcc7 = boxesMPC_ds_cc + 60;
boxesMPC_FLOAT* boxesMPC_ccrhsub7 = boxesMPC_ccrhs + 60;
boxesMPC_FLOAT boxesMPC_Phi7[4];
boxesMPC_FLOAT boxesMPC_W7[4];
boxesMPC_FLOAT boxesMPC_Ysd7[4];
boxesMPC_FLOAT boxesMPC_Lsd7[4];
boxesMPC_FLOAT* boxesMPC_z8 = boxesMPC_z + 32;
boxesMPC_FLOAT* boxesMPC_dzaff8 = boxesMPC_dz_aff + 32;
boxesMPC_FLOAT* boxesMPC_dzcc8 = boxesMPC_dz_cc + 32;
boxesMPC_FLOAT* boxesMPC_rd8 = boxesMPC_rd + 32;
boxesMPC_FLOAT boxesMPC_Lbyrd8[4];
boxesMPC_FLOAT* boxesMPC_grad_cost8 = boxesMPC_grad_cost + 32;
boxesMPC_FLOAT* boxesMPC_grad_eq8 = boxesMPC_grad_eq + 32;
boxesMPC_FLOAT* boxesMPC_grad_ineq8 = boxesMPC_grad_ineq + 32;
boxesMPC_FLOAT boxesMPC_ctv8[4];
boxesMPC_FLOAT* boxesMPC_v8 = boxesMPC_v + 16;
boxesMPC_FLOAT boxesMPC_re8[2];
boxesMPC_FLOAT boxesMPC_beta8[2];
boxesMPC_FLOAT boxesMPC_betacc8[2];
boxesMPC_FLOAT* boxesMPC_dvaff8 = boxesMPC_dv_aff + 16;
boxesMPC_FLOAT* boxesMPC_dvcc8 = boxesMPC_dv_cc + 16;
boxesMPC_FLOAT boxesMPC_V8[8];
boxesMPC_FLOAT boxesMPC_Yd8[3];
boxesMPC_FLOAT boxesMPC_Ld8[3];
boxesMPC_FLOAT boxesMPC_yy8[2];
boxesMPC_FLOAT boxesMPC_bmy8[2];
boxesMPC_FLOAT boxesMPC_c8[2] = {0.0000000000000000E+000, 0.0000000000000000E+000};
int boxesMPC_lbIdx8[4] = {0, 1, 2, 3};
boxesMPC_FLOAT* boxesMPC_llb8 = boxesMPC_l + 64;
boxesMPC_FLOAT* boxesMPC_slb8 = boxesMPC_s + 64;
boxesMPC_FLOAT* boxesMPC_llbbyslb8 = boxesMPC_lbys + 64;
boxesMPC_FLOAT boxesMPC_rilb8[4];
boxesMPC_FLOAT* boxesMPC_dllbaff8 = boxesMPC_dl_aff + 64;
boxesMPC_FLOAT* boxesMPC_dslbaff8 = boxesMPC_ds_aff + 64;
boxesMPC_FLOAT* boxesMPC_dllbcc8 = boxesMPC_dl_cc + 64;
boxesMPC_FLOAT* boxesMPC_dslbcc8 = boxesMPC_ds_cc + 64;
boxesMPC_FLOAT* boxesMPC_ccrhsl8 = boxesMPC_ccrhs + 64;
int boxesMPC_ubIdx8[4] = {0, 1, 2, 3};
boxesMPC_FLOAT* boxesMPC_lub8 = boxesMPC_l + 68;
boxesMPC_FLOAT* boxesMPC_sub8 = boxesMPC_s + 68;
boxesMPC_FLOAT* boxesMPC_lubbysub8 = boxesMPC_lbys + 68;
boxesMPC_FLOAT boxesMPC_riub8[4];
boxesMPC_FLOAT* boxesMPC_dlubaff8 = boxesMPC_dl_aff + 68;
boxesMPC_FLOAT* boxesMPC_dsubaff8 = boxesMPC_ds_aff + 68;
boxesMPC_FLOAT* boxesMPC_dlubcc8 = boxesMPC_dl_cc + 68;
boxesMPC_FLOAT* boxesMPC_dsubcc8 = boxesMPC_ds_cc + 68;
boxesMPC_FLOAT* boxesMPC_ccrhsub8 = boxesMPC_ccrhs + 68;
boxesMPC_FLOAT boxesMPC_Phi8[4];
boxesMPC_FLOAT boxesMPC_W8[4];
boxesMPC_FLOAT boxesMPC_Ysd8[4];
boxesMPC_FLOAT boxesMPC_Lsd8[4];
boxesMPC_FLOAT* boxesMPC_z9 = boxesMPC_z + 36;
boxesMPC_FLOAT* boxesMPC_dzaff9 = boxesMPC_dz_aff + 36;
boxesMPC_FLOAT* boxesMPC_dzcc9 = boxesMPC_dz_cc + 36;
boxesMPC_FLOAT* boxesMPC_rd9 = boxesMPC_rd + 36;
boxesMPC_FLOAT boxesMPC_Lbyrd9[2];
boxesMPC_FLOAT* boxesMPC_grad_cost9 = boxesMPC_grad_cost + 36;
boxesMPC_FLOAT* boxesMPC_grad_eq9 = boxesMPC_grad_eq + 36;
boxesMPC_FLOAT* boxesMPC_grad_ineq9 = boxesMPC_grad_ineq + 36;
boxesMPC_FLOAT boxesMPC_ctv9[2];
boxesMPC_FLOAT* boxesMPC_v9 = boxesMPC_v + 18;
boxesMPC_FLOAT boxesMPC_re9[2];
boxesMPC_FLOAT boxesMPC_beta9[2];
boxesMPC_FLOAT boxesMPC_betacc9[2];
boxesMPC_FLOAT* boxesMPC_dvaff9 = boxesMPC_dv_aff + 18;
boxesMPC_FLOAT* boxesMPC_dvcc9 = boxesMPC_dv_cc + 18;
boxesMPC_FLOAT boxesMPC_V9[4];
boxesMPC_FLOAT boxesMPC_Yd9[3];
boxesMPC_FLOAT boxesMPC_Ld9[3];
boxesMPC_FLOAT boxesMPC_yy9[2];
boxesMPC_FLOAT boxesMPC_bmy9[2];
boxesMPC_FLOAT boxesMPC_c9[2] = {0.0000000000000000E+000, 0.0000000000000000E+000};
int boxesMPC_lbIdx9[2] = {0, 1};
boxesMPC_FLOAT* boxesMPC_llb9 = boxesMPC_l + 72;
boxesMPC_FLOAT* boxesMPC_slb9 = boxesMPC_s + 72;
boxesMPC_FLOAT* boxesMPC_llbbyslb9 = boxesMPC_lbys + 72;
boxesMPC_FLOAT boxesMPC_rilb9[2];
boxesMPC_FLOAT* boxesMPC_dllbaff9 = boxesMPC_dl_aff + 72;
boxesMPC_FLOAT* boxesMPC_dslbaff9 = boxesMPC_ds_aff + 72;
boxesMPC_FLOAT* boxesMPC_dllbcc9 = boxesMPC_dl_cc + 72;
boxesMPC_FLOAT* boxesMPC_dslbcc9 = boxesMPC_ds_cc + 72;
boxesMPC_FLOAT* boxesMPC_ccrhsl9 = boxesMPC_ccrhs + 72;
int boxesMPC_ubIdx9[2] = {0, 1};
boxesMPC_FLOAT* boxesMPC_lub9 = boxesMPC_l + 74;
boxesMPC_FLOAT* boxesMPC_sub9 = boxesMPC_s + 74;
boxesMPC_FLOAT* boxesMPC_lubbysub9 = boxesMPC_lbys + 74;
boxesMPC_FLOAT boxesMPC_riub9[2];
boxesMPC_FLOAT* boxesMPC_dlubaff9 = boxesMPC_dl_aff + 74;
boxesMPC_FLOAT* boxesMPC_dsubaff9 = boxesMPC_ds_aff + 74;
boxesMPC_FLOAT* boxesMPC_dlubcc9 = boxesMPC_dl_cc + 74;
boxesMPC_FLOAT* boxesMPC_dsubcc9 = boxesMPC_ds_cc + 74;
boxesMPC_FLOAT* boxesMPC_ccrhsub9 = boxesMPC_ccrhs + 74;
boxesMPC_FLOAT boxesMPC_Phi9[2];
boxesMPC_FLOAT boxesMPC_D9[2] = {-1.0000000000000000E+000, 
-1.0000000000000000E+000};
boxesMPC_FLOAT boxesMPC_W9[2];
boxesMPC_FLOAT boxesMPC_Ysd9[4];
boxesMPC_FLOAT boxesMPC_Lsd9[4];
boxesMPC_FLOAT musigma;
boxesMPC_FLOAT sigma_3rdroot;
boxesMPC_FLOAT boxesMPC_Diag1_0[4];
boxesMPC_FLOAT boxesMPC_Diag2_0[4];
boxesMPC_FLOAT boxesMPC_L_0[6];




/* SOLVER CODE --------------------------------------------------------- */
int boxesMPC_solve(boxesMPC_params* params, boxesMPC_output* output, boxesMPC_info* info)
{	
int exitcode;

#if boxesMPC_SET_TIMING == 1
	boxesMPC_timer solvertimer;
	boxesMPC_tic(&solvertimer);
#endif
/* FUNCTION CALLS INTO LA LIBRARY -------------------------------------- */
info->it = 0;
boxesMPC_LA_INITIALIZEVECTOR_38(boxesMPC_z, 0);
boxesMPC_LA_INITIALIZEVECTOR_20(boxesMPC_v, 1);
boxesMPC_LA_INITIALIZEVECTOR_76(boxesMPC_l, 1);
boxesMPC_LA_INITIALIZEVECTOR_76(boxesMPC_s, 1);
info->mu = 0;
boxesMPC_LA_DOTACC_76(boxesMPC_l, boxesMPC_s, &info->mu);
info->mu /= 76;
while( 1 ){
info->pobj = 0;
boxesMPC_LA_DIAG_QUADFCN_4(params->H1, params->f1, boxesMPC_z0, boxesMPC_grad_cost0, &info->pobj);
boxesMPC_LA_DIAG_QUADFCN_4(params->H2, params->f2, boxesMPC_z1, boxesMPC_grad_cost1, &info->pobj);
boxesMPC_LA_DIAG_QUADFCN_4(params->H3, params->f3, boxesMPC_z2, boxesMPC_grad_cost2, &info->pobj);
boxesMPC_LA_DIAG_QUADFCN_4(params->H4, params->f4, boxesMPC_z3, boxesMPC_grad_cost3, &info->pobj);
boxesMPC_LA_DIAG_QUADFCN_4(params->H5, params->f5, boxesMPC_z4, boxesMPC_grad_cost4, &info->pobj);
boxesMPC_LA_DIAG_QUADFCN_4(params->H6, params->f6, boxesMPC_z5, boxesMPC_grad_cost5, &info->pobj);
boxesMPC_LA_DIAG_QUADFCN_4(params->H7, params->f7, boxesMPC_z6, boxesMPC_grad_cost6, &info->pobj);
boxesMPC_LA_DIAG_QUADFCN_4(params->H8, params->f8, boxesMPC_z7, boxesMPC_grad_cost7, &info->pobj);
boxesMPC_LA_DIAG_QUADFCN_4(params->H9, params->f9, boxesMPC_z8, boxesMPC_grad_cost8, &info->pobj);
boxesMPC_LA_DIAG_QUADFCN_2(params->H10, params->f10, boxesMPC_z9, boxesMPC_grad_cost9, &info->pobj);
info->res_eq = 0;
info->dgap = 0;
boxesMPC_LA_DIAGZERO_MVMSUB6_2(boxesMPC_D0, boxesMPC_z0, params->c1, boxesMPC_v0, boxesMPC_re0, &info->dgap, &info->res_eq);
boxesMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_4_4(boxesMPC_C0, boxesMPC_z0, boxesMPC_D1, boxesMPC_z1, boxesMPC_c1, boxesMPC_v1, boxesMPC_re1, &info->dgap, &info->res_eq);
boxesMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_4_4(boxesMPC_C0, boxesMPC_z1, boxesMPC_D1, boxesMPC_z2, boxesMPC_c2, boxesMPC_v2, boxesMPC_re2, &info->dgap, &info->res_eq);
boxesMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_4_4(boxesMPC_C0, boxesMPC_z2, boxesMPC_D1, boxesMPC_z3, boxesMPC_c3, boxesMPC_v3, boxesMPC_re3, &info->dgap, &info->res_eq);
boxesMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_4_4(boxesMPC_C0, boxesMPC_z3, boxesMPC_D1, boxesMPC_z4, boxesMPC_c4, boxesMPC_v4, boxesMPC_re4, &info->dgap, &info->res_eq);
boxesMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_4_4(boxesMPC_C0, boxesMPC_z4, boxesMPC_D1, boxesMPC_z5, boxesMPC_c5, boxesMPC_v5, boxesMPC_re5, &info->dgap, &info->res_eq);
boxesMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_4_4(boxesMPC_C0, boxesMPC_z5, boxesMPC_D1, boxesMPC_z6, boxesMPC_c6, boxesMPC_v6, boxesMPC_re6, &info->dgap, &info->res_eq);
boxesMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_4_4(boxesMPC_C0, boxesMPC_z6, boxesMPC_D1, boxesMPC_z7, boxesMPC_c7, boxesMPC_v7, boxesMPC_re7, &info->dgap, &info->res_eq);
boxesMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_4_4(boxesMPC_C0, boxesMPC_z7, boxesMPC_D1, boxesMPC_z8, boxesMPC_c8, boxesMPC_v8, boxesMPC_re8, &info->dgap, &info->res_eq);
boxesMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_4_2(boxesMPC_C0, boxesMPC_z8, boxesMPC_D9, boxesMPC_z9, boxesMPC_c9, boxesMPC_v9, boxesMPC_re9, &info->dgap, &info->res_eq);
boxesMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(boxesMPC_C0, boxesMPC_v1, boxesMPC_D0, boxesMPC_v0, boxesMPC_grad_eq0);
boxesMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(boxesMPC_C0, boxesMPC_v2, boxesMPC_D1, boxesMPC_v1, boxesMPC_grad_eq1);
boxesMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(boxesMPC_C0, boxesMPC_v3, boxesMPC_D1, boxesMPC_v2, boxesMPC_grad_eq2);
boxesMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(boxesMPC_C0, boxesMPC_v4, boxesMPC_D1, boxesMPC_v3, boxesMPC_grad_eq3);
boxesMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(boxesMPC_C0, boxesMPC_v5, boxesMPC_D1, boxesMPC_v4, boxesMPC_grad_eq4);
boxesMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(boxesMPC_C0, boxesMPC_v6, boxesMPC_D1, boxesMPC_v5, boxesMPC_grad_eq5);
boxesMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(boxesMPC_C0, boxesMPC_v7, boxesMPC_D1, boxesMPC_v6, boxesMPC_grad_eq6);
boxesMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(boxesMPC_C0, boxesMPC_v8, boxesMPC_D1, boxesMPC_v7, boxesMPC_grad_eq7);
boxesMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(boxesMPC_C0, boxesMPC_v9, boxesMPC_D1, boxesMPC_v8, boxesMPC_grad_eq8);
boxesMPC_LA_DIAGZERO_MTVM_2_2(boxesMPC_D9, boxesMPC_v9, boxesMPC_grad_eq9);
info->res_ineq = 0;
boxesMPC_LA_VSUBADD3_4(params->lb1, boxesMPC_z0, boxesMPC_lbIdx0, boxesMPC_llb0, boxesMPC_slb0, boxesMPC_rilb0, &info->dgap, &info->res_ineq);
boxesMPC_LA_VSUBADD2_4(boxesMPC_z0, boxesMPC_ubIdx0, params->ub1, boxesMPC_lub0, boxesMPC_sub0, boxesMPC_riub0, &info->dgap, &info->res_ineq);
boxesMPC_LA_VSUBADD3_4(params->lb2, boxesMPC_z1, boxesMPC_lbIdx1, boxesMPC_llb1, boxesMPC_slb1, boxesMPC_rilb1, &info->dgap, &info->res_ineq);
boxesMPC_LA_VSUBADD2_4(boxesMPC_z1, boxesMPC_ubIdx1, params->ub2, boxesMPC_lub1, boxesMPC_sub1, boxesMPC_riub1, &info->dgap, &info->res_ineq);
boxesMPC_LA_VSUBADD3_4(params->lb3, boxesMPC_z2, boxesMPC_lbIdx2, boxesMPC_llb2, boxesMPC_slb2, boxesMPC_rilb2, &info->dgap, &info->res_ineq);
boxesMPC_LA_VSUBADD2_4(boxesMPC_z2, boxesMPC_ubIdx2, params->ub3, boxesMPC_lub2, boxesMPC_sub2, boxesMPC_riub2, &info->dgap, &info->res_ineq);
boxesMPC_LA_VSUBADD3_4(params->lb4, boxesMPC_z3, boxesMPC_lbIdx3, boxesMPC_llb3, boxesMPC_slb3, boxesMPC_rilb3, &info->dgap, &info->res_ineq);
boxesMPC_LA_VSUBADD2_4(boxesMPC_z3, boxesMPC_ubIdx3, params->ub4, boxesMPC_lub3, boxesMPC_sub3, boxesMPC_riub3, &info->dgap, &info->res_ineq);
boxesMPC_LA_VSUBADD3_4(params->lb5, boxesMPC_z4, boxesMPC_lbIdx4, boxesMPC_llb4, boxesMPC_slb4, boxesMPC_rilb4, &info->dgap, &info->res_ineq);
boxesMPC_LA_VSUBADD2_4(boxesMPC_z4, boxesMPC_ubIdx4, params->ub5, boxesMPC_lub4, boxesMPC_sub4, boxesMPC_riub4, &info->dgap, &info->res_ineq);
boxesMPC_LA_VSUBADD3_4(params->lb6, boxesMPC_z5, boxesMPC_lbIdx5, boxesMPC_llb5, boxesMPC_slb5, boxesMPC_rilb5, &info->dgap, &info->res_ineq);
boxesMPC_LA_VSUBADD2_4(boxesMPC_z5, boxesMPC_ubIdx5, params->ub6, boxesMPC_lub5, boxesMPC_sub5, boxesMPC_riub5, &info->dgap, &info->res_ineq);
boxesMPC_LA_VSUBADD3_4(params->lb7, boxesMPC_z6, boxesMPC_lbIdx6, boxesMPC_llb6, boxesMPC_slb6, boxesMPC_rilb6, &info->dgap, &info->res_ineq);
boxesMPC_LA_VSUBADD2_4(boxesMPC_z6, boxesMPC_ubIdx6, params->ub7, boxesMPC_lub6, boxesMPC_sub6, boxesMPC_riub6, &info->dgap, &info->res_ineq);
boxesMPC_LA_VSUBADD3_4(params->lb8, boxesMPC_z7, boxesMPC_lbIdx7, boxesMPC_llb7, boxesMPC_slb7, boxesMPC_rilb7, &info->dgap, &info->res_ineq);
boxesMPC_LA_VSUBADD2_4(boxesMPC_z7, boxesMPC_ubIdx7, params->ub8, boxesMPC_lub7, boxesMPC_sub7, boxesMPC_riub7, &info->dgap, &info->res_ineq);
boxesMPC_LA_VSUBADD3_4(params->lb9, boxesMPC_z8, boxesMPC_lbIdx8, boxesMPC_llb8, boxesMPC_slb8, boxesMPC_rilb8, &info->dgap, &info->res_ineq);
boxesMPC_LA_VSUBADD2_4(boxesMPC_z8, boxesMPC_ubIdx8, params->ub9, boxesMPC_lub8, boxesMPC_sub8, boxesMPC_riub8, &info->dgap, &info->res_ineq);
boxesMPC_LA_VSUBADD3_2(params->lb10, boxesMPC_z9, boxesMPC_lbIdx9, boxesMPC_llb9, boxesMPC_slb9, boxesMPC_rilb9, &info->dgap, &info->res_ineq);
boxesMPC_LA_VSUBADD2_2(boxesMPC_z9, boxesMPC_ubIdx9, params->ub10, boxesMPC_lub9, boxesMPC_sub9, boxesMPC_riub9, &info->dgap, &info->res_ineq);
boxesMPC_LA_INEQ_B_GRAD_4_4_4(boxesMPC_lub0, boxesMPC_sub0, boxesMPC_riub0, boxesMPC_llb0, boxesMPC_slb0, boxesMPC_rilb0, boxesMPC_lbIdx0, boxesMPC_ubIdx0, boxesMPC_grad_ineq0, boxesMPC_lubbysub0, boxesMPC_llbbyslb0);
boxesMPC_LA_INEQ_B_GRAD_4_4_4(boxesMPC_lub1, boxesMPC_sub1, boxesMPC_riub1, boxesMPC_llb1, boxesMPC_slb1, boxesMPC_rilb1, boxesMPC_lbIdx1, boxesMPC_ubIdx1, boxesMPC_grad_ineq1, boxesMPC_lubbysub1, boxesMPC_llbbyslb1);
boxesMPC_LA_INEQ_B_GRAD_4_4_4(boxesMPC_lub2, boxesMPC_sub2, boxesMPC_riub2, boxesMPC_llb2, boxesMPC_slb2, boxesMPC_rilb2, boxesMPC_lbIdx2, boxesMPC_ubIdx2, boxesMPC_grad_ineq2, boxesMPC_lubbysub2, boxesMPC_llbbyslb2);
boxesMPC_LA_INEQ_B_GRAD_4_4_4(boxesMPC_lub3, boxesMPC_sub3, boxesMPC_riub3, boxesMPC_llb3, boxesMPC_slb3, boxesMPC_rilb3, boxesMPC_lbIdx3, boxesMPC_ubIdx3, boxesMPC_grad_ineq3, boxesMPC_lubbysub3, boxesMPC_llbbyslb3);
boxesMPC_LA_INEQ_B_GRAD_4_4_4(boxesMPC_lub4, boxesMPC_sub4, boxesMPC_riub4, boxesMPC_llb4, boxesMPC_slb4, boxesMPC_rilb4, boxesMPC_lbIdx4, boxesMPC_ubIdx4, boxesMPC_grad_ineq4, boxesMPC_lubbysub4, boxesMPC_llbbyslb4);
boxesMPC_LA_INEQ_B_GRAD_4_4_4(boxesMPC_lub5, boxesMPC_sub5, boxesMPC_riub5, boxesMPC_llb5, boxesMPC_slb5, boxesMPC_rilb5, boxesMPC_lbIdx5, boxesMPC_ubIdx5, boxesMPC_grad_ineq5, boxesMPC_lubbysub5, boxesMPC_llbbyslb5);
boxesMPC_LA_INEQ_B_GRAD_4_4_4(boxesMPC_lub6, boxesMPC_sub6, boxesMPC_riub6, boxesMPC_llb6, boxesMPC_slb6, boxesMPC_rilb6, boxesMPC_lbIdx6, boxesMPC_ubIdx6, boxesMPC_grad_ineq6, boxesMPC_lubbysub6, boxesMPC_llbbyslb6);
boxesMPC_LA_INEQ_B_GRAD_4_4_4(boxesMPC_lub7, boxesMPC_sub7, boxesMPC_riub7, boxesMPC_llb7, boxesMPC_slb7, boxesMPC_rilb7, boxesMPC_lbIdx7, boxesMPC_ubIdx7, boxesMPC_grad_ineq7, boxesMPC_lubbysub7, boxesMPC_llbbyslb7);
boxesMPC_LA_INEQ_B_GRAD_4_4_4(boxesMPC_lub8, boxesMPC_sub8, boxesMPC_riub8, boxesMPC_llb8, boxesMPC_slb8, boxesMPC_rilb8, boxesMPC_lbIdx8, boxesMPC_ubIdx8, boxesMPC_grad_ineq8, boxesMPC_lubbysub8, boxesMPC_llbbyslb8);
boxesMPC_LA_INEQ_B_GRAD_2_2_2(boxesMPC_lub9, boxesMPC_sub9, boxesMPC_riub9, boxesMPC_llb9, boxesMPC_slb9, boxesMPC_rilb9, boxesMPC_lbIdx9, boxesMPC_ubIdx9, boxesMPC_grad_ineq9, boxesMPC_lubbysub9, boxesMPC_llbbyslb9);
info->dobj = info->pobj - info->dgap;
info->rdgap = info->pobj ? info->dgap / info->pobj : 1e6;
if( info->rdgap < 0 ) info->rdgap = -info->rdgap;
if( info->mu < boxesMPC_SET_ACC_KKTCOMPL
    && (info->rdgap < boxesMPC_SET_ACC_RDGAP || info->dgap < boxesMPC_SET_ACC_KKTCOMPL)
    && info->res_eq < boxesMPC_SET_ACC_RESEQ
    && info->res_ineq < boxesMPC_SET_ACC_RESINEQ ){
exitcode = boxesMPC_OPTIMAL; break; }
if( info->it == boxesMPC_SET_MAXIT ){
exitcode = boxesMPC_MAXITREACHED; break; }
boxesMPC_LA_VVADD3_38(boxesMPC_grad_cost, boxesMPC_grad_eq, boxesMPC_grad_ineq, boxesMPC_rd);
boxesMPC_LA_DIAG_CHOL_ONELOOP_LBUB_4_4_4(params->H1, boxesMPC_llbbyslb0, boxesMPC_lbIdx0, boxesMPC_lubbysub0, boxesMPC_ubIdx0, boxesMPC_Phi0);
boxesMPC_LA_DIAG_MATRIXFORWARDSUB_2_4(boxesMPC_Phi0, boxesMPC_C0, boxesMPC_V0);
boxesMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_4(boxesMPC_Phi0, boxesMPC_D0, boxesMPC_W0);
boxesMPC_LA_DENSE_DIAGZERO_MMTM_2_4_2(boxesMPC_W0, boxesMPC_V0, boxesMPC_Ysd1);
boxesMPC_LA_DIAG_FORWARDSUB_4(boxesMPC_Phi0, boxesMPC_rd0, boxesMPC_Lbyrd0);
boxesMPC_LA_DIAG_CHOL_ONELOOP_LBUB_4_4_4(params->H2, boxesMPC_llbbyslb1, boxesMPC_lbIdx1, boxesMPC_lubbysub1, boxesMPC_ubIdx1, boxesMPC_Phi1);
boxesMPC_LA_DIAG_MATRIXFORWARDSUB_2_4(boxesMPC_Phi1, boxesMPC_C0, boxesMPC_V1);
boxesMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_4(boxesMPC_Phi1, boxesMPC_D1, boxesMPC_W1);
boxesMPC_LA_DENSE_DIAGZERO_MMTM_2_4_2(boxesMPC_W1, boxesMPC_V1, boxesMPC_Ysd2);
boxesMPC_LA_DIAG_FORWARDSUB_4(boxesMPC_Phi1, boxesMPC_rd1, boxesMPC_Lbyrd1);
boxesMPC_LA_DIAG_CHOL_ONELOOP_LBUB_4_4_4(params->H3, boxesMPC_llbbyslb2, boxesMPC_lbIdx2, boxesMPC_lubbysub2, boxesMPC_ubIdx2, boxesMPC_Phi2);
boxesMPC_LA_DIAG_MATRIXFORWARDSUB_2_4(boxesMPC_Phi2, boxesMPC_C0, boxesMPC_V2);
boxesMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_4(boxesMPC_Phi2, boxesMPC_D1, boxesMPC_W2);
boxesMPC_LA_DENSE_DIAGZERO_MMTM_2_4_2(boxesMPC_W2, boxesMPC_V2, boxesMPC_Ysd3);
boxesMPC_LA_DIAG_FORWARDSUB_4(boxesMPC_Phi2, boxesMPC_rd2, boxesMPC_Lbyrd2);
boxesMPC_LA_DIAG_CHOL_ONELOOP_LBUB_4_4_4(params->H4, boxesMPC_llbbyslb3, boxesMPC_lbIdx3, boxesMPC_lubbysub3, boxesMPC_ubIdx3, boxesMPC_Phi3);
boxesMPC_LA_DIAG_MATRIXFORWARDSUB_2_4(boxesMPC_Phi3, boxesMPC_C0, boxesMPC_V3);
boxesMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_4(boxesMPC_Phi3, boxesMPC_D1, boxesMPC_W3);
boxesMPC_LA_DENSE_DIAGZERO_MMTM_2_4_2(boxesMPC_W3, boxesMPC_V3, boxesMPC_Ysd4);
boxesMPC_LA_DIAG_FORWARDSUB_4(boxesMPC_Phi3, boxesMPC_rd3, boxesMPC_Lbyrd3);
boxesMPC_LA_DIAG_CHOL_ONELOOP_LBUB_4_4_4(params->H5, boxesMPC_llbbyslb4, boxesMPC_lbIdx4, boxesMPC_lubbysub4, boxesMPC_ubIdx4, boxesMPC_Phi4);
boxesMPC_LA_DIAG_MATRIXFORWARDSUB_2_4(boxesMPC_Phi4, boxesMPC_C0, boxesMPC_V4);
boxesMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_4(boxesMPC_Phi4, boxesMPC_D1, boxesMPC_W4);
boxesMPC_LA_DENSE_DIAGZERO_MMTM_2_4_2(boxesMPC_W4, boxesMPC_V4, boxesMPC_Ysd5);
boxesMPC_LA_DIAG_FORWARDSUB_4(boxesMPC_Phi4, boxesMPC_rd4, boxesMPC_Lbyrd4);
boxesMPC_LA_DIAG_CHOL_ONELOOP_LBUB_4_4_4(params->H6, boxesMPC_llbbyslb5, boxesMPC_lbIdx5, boxesMPC_lubbysub5, boxesMPC_ubIdx5, boxesMPC_Phi5);
boxesMPC_LA_DIAG_MATRIXFORWARDSUB_2_4(boxesMPC_Phi5, boxesMPC_C0, boxesMPC_V5);
boxesMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_4(boxesMPC_Phi5, boxesMPC_D1, boxesMPC_W5);
boxesMPC_LA_DENSE_DIAGZERO_MMTM_2_4_2(boxesMPC_W5, boxesMPC_V5, boxesMPC_Ysd6);
boxesMPC_LA_DIAG_FORWARDSUB_4(boxesMPC_Phi5, boxesMPC_rd5, boxesMPC_Lbyrd5);
boxesMPC_LA_DIAG_CHOL_ONELOOP_LBUB_4_4_4(params->H7, boxesMPC_llbbyslb6, boxesMPC_lbIdx6, boxesMPC_lubbysub6, boxesMPC_ubIdx6, boxesMPC_Phi6);
boxesMPC_LA_DIAG_MATRIXFORWARDSUB_2_4(boxesMPC_Phi6, boxesMPC_C0, boxesMPC_V6);
boxesMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_4(boxesMPC_Phi6, boxesMPC_D1, boxesMPC_W6);
boxesMPC_LA_DENSE_DIAGZERO_MMTM_2_4_2(boxesMPC_W6, boxesMPC_V6, boxesMPC_Ysd7);
boxesMPC_LA_DIAG_FORWARDSUB_4(boxesMPC_Phi6, boxesMPC_rd6, boxesMPC_Lbyrd6);
boxesMPC_LA_DIAG_CHOL_ONELOOP_LBUB_4_4_4(params->H8, boxesMPC_llbbyslb7, boxesMPC_lbIdx7, boxesMPC_lubbysub7, boxesMPC_ubIdx7, boxesMPC_Phi7);
boxesMPC_LA_DIAG_MATRIXFORWARDSUB_2_4(boxesMPC_Phi7, boxesMPC_C0, boxesMPC_V7);
boxesMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_4(boxesMPC_Phi7, boxesMPC_D1, boxesMPC_W7);
boxesMPC_LA_DENSE_DIAGZERO_MMTM_2_4_2(boxesMPC_W7, boxesMPC_V7, boxesMPC_Ysd8);
boxesMPC_LA_DIAG_FORWARDSUB_4(boxesMPC_Phi7, boxesMPC_rd7, boxesMPC_Lbyrd7);
boxesMPC_LA_DIAG_CHOL_ONELOOP_LBUB_4_4_4(params->H9, boxesMPC_llbbyslb8, boxesMPC_lbIdx8, boxesMPC_lubbysub8, boxesMPC_ubIdx8, boxesMPC_Phi8);
boxesMPC_LA_DIAG_MATRIXFORWARDSUB_2_4(boxesMPC_Phi8, boxesMPC_C0, boxesMPC_V8);
boxesMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_4(boxesMPC_Phi8, boxesMPC_D1, boxesMPC_W8);
boxesMPC_LA_DENSE_DIAGZERO_MMTM_2_4_2(boxesMPC_W8, boxesMPC_V8, boxesMPC_Ysd9);
boxesMPC_LA_DIAG_FORWARDSUB_4(boxesMPC_Phi8, boxesMPC_rd8, boxesMPC_Lbyrd8);
boxesMPC_LA_DIAG_CHOL_ONELOOP_LBUB_2_2_2(params->H10, boxesMPC_llbbyslb9, boxesMPC_lbIdx9, boxesMPC_lubbysub9, boxesMPC_ubIdx9, boxesMPC_Phi9);
boxesMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_2(boxesMPC_Phi9, boxesMPC_D9, boxesMPC_W9);
boxesMPC_LA_DIAG_FORWARDSUB_2(boxesMPC_Phi9, boxesMPC_rd9, boxesMPC_Lbyrd9);
boxesMPC_LA_DIAGZERO_MMT_2(boxesMPC_W0, boxesMPC_Yd0);
boxesMPC_LA_DIAGZERO_MVMSUB7_2(boxesMPC_W0, boxesMPC_Lbyrd0, boxesMPC_re0, boxesMPC_beta0);
boxesMPC_LA_DENSE_DIAGZERO_MMT2_2_4_4(boxesMPC_V0, boxesMPC_W1, boxesMPC_Yd1);
boxesMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_4_4(boxesMPC_V0, boxesMPC_Lbyrd0, boxesMPC_W1, boxesMPC_Lbyrd1, boxesMPC_re1, boxesMPC_beta1);
boxesMPC_LA_DENSE_DIAGZERO_MMT2_2_4_4(boxesMPC_V1, boxesMPC_W2, boxesMPC_Yd2);
boxesMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_4_4(boxesMPC_V1, boxesMPC_Lbyrd1, boxesMPC_W2, boxesMPC_Lbyrd2, boxesMPC_re2, boxesMPC_beta2);
boxesMPC_LA_DENSE_DIAGZERO_MMT2_2_4_4(boxesMPC_V2, boxesMPC_W3, boxesMPC_Yd3);
boxesMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_4_4(boxesMPC_V2, boxesMPC_Lbyrd2, boxesMPC_W3, boxesMPC_Lbyrd3, boxesMPC_re3, boxesMPC_beta3);
boxesMPC_LA_DENSE_DIAGZERO_MMT2_2_4_4(boxesMPC_V3, boxesMPC_W4, boxesMPC_Yd4);
boxesMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_4_4(boxesMPC_V3, boxesMPC_Lbyrd3, boxesMPC_W4, boxesMPC_Lbyrd4, boxesMPC_re4, boxesMPC_beta4);
boxesMPC_LA_DENSE_DIAGZERO_MMT2_2_4_4(boxesMPC_V4, boxesMPC_W5, boxesMPC_Yd5);
boxesMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_4_4(boxesMPC_V4, boxesMPC_Lbyrd4, boxesMPC_W5, boxesMPC_Lbyrd5, boxesMPC_re5, boxesMPC_beta5);
boxesMPC_LA_DENSE_DIAGZERO_MMT2_2_4_4(boxesMPC_V5, boxesMPC_W6, boxesMPC_Yd6);
boxesMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_4_4(boxesMPC_V5, boxesMPC_Lbyrd5, boxesMPC_W6, boxesMPC_Lbyrd6, boxesMPC_re6, boxesMPC_beta6);
boxesMPC_LA_DENSE_DIAGZERO_MMT2_2_4_4(boxesMPC_V6, boxesMPC_W7, boxesMPC_Yd7);
boxesMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_4_4(boxesMPC_V6, boxesMPC_Lbyrd6, boxesMPC_W7, boxesMPC_Lbyrd7, boxesMPC_re7, boxesMPC_beta7);
boxesMPC_LA_DENSE_DIAGZERO_MMT2_2_4_4(boxesMPC_V7, boxesMPC_W8, boxesMPC_Yd8);
boxesMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_4_4(boxesMPC_V7, boxesMPC_Lbyrd7, boxesMPC_W8, boxesMPC_Lbyrd8, boxesMPC_re8, boxesMPC_beta8);
boxesMPC_LA_DENSE_DIAGZERO_MMT2_2_4_2(boxesMPC_V8, boxesMPC_W9, boxesMPC_Yd9);
boxesMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_4_2(boxesMPC_V8, boxesMPC_Lbyrd8, boxesMPC_W9, boxesMPC_Lbyrd9, boxesMPC_re9, boxesMPC_beta9);
boxesMPC_LA_DENSE_CHOL_2(boxesMPC_Yd0, boxesMPC_Ld0);
boxesMPC_LA_DENSE_FORWARDSUB_2(boxesMPC_Ld0, boxesMPC_beta0, boxesMPC_yy0);
boxesMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(boxesMPC_Ld0, boxesMPC_Ysd1, boxesMPC_Lsd1);
boxesMPC_LA_DENSE_MMTSUB_2_2(boxesMPC_Lsd1, boxesMPC_Yd1);
boxesMPC_LA_DENSE_CHOL_2(boxesMPC_Yd1, boxesMPC_Ld1);
boxesMPC_LA_DENSE_MVMSUB1_2_2(boxesMPC_Lsd1, boxesMPC_yy0, boxesMPC_beta1, boxesMPC_bmy1);
boxesMPC_LA_DENSE_FORWARDSUB_2(boxesMPC_Ld1, boxesMPC_bmy1, boxesMPC_yy1);
boxesMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(boxesMPC_Ld1, boxesMPC_Ysd2, boxesMPC_Lsd2);
boxesMPC_LA_DENSE_MMTSUB_2_2(boxesMPC_Lsd2, boxesMPC_Yd2);
boxesMPC_LA_DENSE_CHOL_2(boxesMPC_Yd2, boxesMPC_Ld2);
boxesMPC_LA_DENSE_MVMSUB1_2_2(boxesMPC_Lsd2, boxesMPC_yy1, boxesMPC_beta2, boxesMPC_bmy2);
boxesMPC_LA_DENSE_FORWARDSUB_2(boxesMPC_Ld2, boxesMPC_bmy2, boxesMPC_yy2);
boxesMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(boxesMPC_Ld2, boxesMPC_Ysd3, boxesMPC_Lsd3);
boxesMPC_LA_DENSE_MMTSUB_2_2(boxesMPC_Lsd3, boxesMPC_Yd3);
boxesMPC_LA_DENSE_CHOL_2(boxesMPC_Yd3, boxesMPC_Ld3);
boxesMPC_LA_DENSE_MVMSUB1_2_2(boxesMPC_Lsd3, boxesMPC_yy2, boxesMPC_beta3, boxesMPC_bmy3);
boxesMPC_LA_DENSE_FORWARDSUB_2(boxesMPC_Ld3, boxesMPC_bmy3, boxesMPC_yy3);
boxesMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(boxesMPC_Ld3, boxesMPC_Ysd4, boxesMPC_Lsd4);
boxesMPC_LA_DENSE_MMTSUB_2_2(boxesMPC_Lsd4, boxesMPC_Yd4);
boxesMPC_LA_DENSE_CHOL_2(boxesMPC_Yd4, boxesMPC_Ld4);
boxesMPC_LA_DENSE_MVMSUB1_2_2(boxesMPC_Lsd4, boxesMPC_yy3, boxesMPC_beta4, boxesMPC_bmy4);
boxesMPC_LA_DENSE_FORWARDSUB_2(boxesMPC_Ld4, boxesMPC_bmy4, boxesMPC_yy4);
boxesMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(boxesMPC_Ld4, boxesMPC_Ysd5, boxesMPC_Lsd5);
boxesMPC_LA_DENSE_MMTSUB_2_2(boxesMPC_Lsd5, boxesMPC_Yd5);
boxesMPC_LA_DENSE_CHOL_2(boxesMPC_Yd5, boxesMPC_Ld5);
boxesMPC_LA_DENSE_MVMSUB1_2_2(boxesMPC_Lsd5, boxesMPC_yy4, boxesMPC_beta5, boxesMPC_bmy5);
boxesMPC_LA_DENSE_FORWARDSUB_2(boxesMPC_Ld5, boxesMPC_bmy5, boxesMPC_yy5);
boxesMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(boxesMPC_Ld5, boxesMPC_Ysd6, boxesMPC_Lsd6);
boxesMPC_LA_DENSE_MMTSUB_2_2(boxesMPC_Lsd6, boxesMPC_Yd6);
boxesMPC_LA_DENSE_CHOL_2(boxesMPC_Yd6, boxesMPC_Ld6);
boxesMPC_LA_DENSE_MVMSUB1_2_2(boxesMPC_Lsd6, boxesMPC_yy5, boxesMPC_beta6, boxesMPC_bmy6);
boxesMPC_LA_DENSE_FORWARDSUB_2(boxesMPC_Ld6, boxesMPC_bmy6, boxesMPC_yy6);
boxesMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(boxesMPC_Ld6, boxesMPC_Ysd7, boxesMPC_Lsd7);
boxesMPC_LA_DENSE_MMTSUB_2_2(boxesMPC_Lsd7, boxesMPC_Yd7);
boxesMPC_LA_DENSE_CHOL_2(boxesMPC_Yd7, boxesMPC_Ld7);
boxesMPC_LA_DENSE_MVMSUB1_2_2(boxesMPC_Lsd7, boxesMPC_yy6, boxesMPC_beta7, boxesMPC_bmy7);
boxesMPC_LA_DENSE_FORWARDSUB_2(boxesMPC_Ld7, boxesMPC_bmy7, boxesMPC_yy7);
boxesMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(boxesMPC_Ld7, boxesMPC_Ysd8, boxesMPC_Lsd8);
boxesMPC_LA_DENSE_MMTSUB_2_2(boxesMPC_Lsd8, boxesMPC_Yd8);
boxesMPC_LA_DENSE_CHOL_2(boxesMPC_Yd8, boxesMPC_Ld8);
boxesMPC_LA_DENSE_MVMSUB1_2_2(boxesMPC_Lsd8, boxesMPC_yy7, boxesMPC_beta8, boxesMPC_bmy8);
boxesMPC_LA_DENSE_FORWARDSUB_2(boxesMPC_Ld8, boxesMPC_bmy8, boxesMPC_yy8);
boxesMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(boxesMPC_Ld8, boxesMPC_Ysd9, boxesMPC_Lsd9);
boxesMPC_LA_DENSE_MMTSUB_2_2(boxesMPC_Lsd9, boxesMPC_Yd9);
boxesMPC_LA_DENSE_CHOL_2(boxesMPC_Yd9, boxesMPC_Ld9);
boxesMPC_LA_DENSE_MVMSUB1_2_2(boxesMPC_Lsd9, boxesMPC_yy8, boxesMPC_beta9, boxesMPC_bmy9);
boxesMPC_LA_DENSE_FORWARDSUB_2(boxesMPC_Ld9, boxesMPC_bmy9, boxesMPC_yy9);
boxesMPC_LA_DENSE_BACKWARDSUB_2(boxesMPC_Ld9, boxesMPC_yy9, boxesMPC_dvaff9);
boxesMPC_LA_DENSE_MTVMSUB_2_2(boxesMPC_Lsd9, boxesMPC_dvaff9, boxesMPC_yy8, boxesMPC_bmy8);
boxesMPC_LA_DENSE_BACKWARDSUB_2(boxesMPC_Ld8, boxesMPC_bmy8, boxesMPC_dvaff8);
boxesMPC_LA_DENSE_MTVMSUB_2_2(boxesMPC_Lsd8, boxesMPC_dvaff8, boxesMPC_yy7, boxesMPC_bmy7);
boxesMPC_LA_DENSE_BACKWARDSUB_2(boxesMPC_Ld7, boxesMPC_bmy7, boxesMPC_dvaff7);
boxesMPC_LA_DENSE_MTVMSUB_2_2(boxesMPC_Lsd7, boxesMPC_dvaff7, boxesMPC_yy6, boxesMPC_bmy6);
boxesMPC_LA_DENSE_BACKWARDSUB_2(boxesMPC_Ld6, boxesMPC_bmy6, boxesMPC_dvaff6);
boxesMPC_LA_DENSE_MTVMSUB_2_2(boxesMPC_Lsd6, boxesMPC_dvaff6, boxesMPC_yy5, boxesMPC_bmy5);
boxesMPC_LA_DENSE_BACKWARDSUB_2(boxesMPC_Ld5, boxesMPC_bmy5, boxesMPC_dvaff5);
boxesMPC_LA_DENSE_MTVMSUB_2_2(boxesMPC_Lsd5, boxesMPC_dvaff5, boxesMPC_yy4, boxesMPC_bmy4);
boxesMPC_LA_DENSE_BACKWARDSUB_2(boxesMPC_Ld4, boxesMPC_bmy4, boxesMPC_dvaff4);
boxesMPC_LA_DENSE_MTVMSUB_2_2(boxesMPC_Lsd4, boxesMPC_dvaff4, boxesMPC_yy3, boxesMPC_bmy3);
boxesMPC_LA_DENSE_BACKWARDSUB_2(boxesMPC_Ld3, boxesMPC_bmy3, boxesMPC_dvaff3);
boxesMPC_LA_DENSE_MTVMSUB_2_2(boxesMPC_Lsd3, boxesMPC_dvaff3, boxesMPC_yy2, boxesMPC_bmy2);
boxesMPC_LA_DENSE_BACKWARDSUB_2(boxesMPC_Ld2, boxesMPC_bmy2, boxesMPC_dvaff2);
boxesMPC_LA_DENSE_MTVMSUB_2_2(boxesMPC_Lsd2, boxesMPC_dvaff2, boxesMPC_yy1, boxesMPC_bmy1);
boxesMPC_LA_DENSE_BACKWARDSUB_2(boxesMPC_Ld1, boxesMPC_bmy1, boxesMPC_dvaff1);
boxesMPC_LA_DENSE_MTVMSUB_2_2(boxesMPC_Lsd1, boxesMPC_dvaff1, boxesMPC_yy0, boxesMPC_bmy0);
boxesMPC_LA_DENSE_BACKWARDSUB_2(boxesMPC_Ld0, boxesMPC_bmy0, boxesMPC_dvaff0);
boxesMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(boxesMPC_C0, boxesMPC_dvaff1, boxesMPC_D0, boxesMPC_dvaff0, boxesMPC_grad_eq0);
boxesMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(boxesMPC_C0, boxesMPC_dvaff2, boxesMPC_D1, boxesMPC_dvaff1, boxesMPC_grad_eq1);
boxesMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(boxesMPC_C0, boxesMPC_dvaff3, boxesMPC_D1, boxesMPC_dvaff2, boxesMPC_grad_eq2);
boxesMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(boxesMPC_C0, boxesMPC_dvaff4, boxesMPC_D1, boxesMPC_dvaff3, boxesMPC_grad_eq3);
boxesMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(boxesMPC_C0, boxesMPC_dvaff5, boxesMPC_D1, boxesMPC_dvaff4, boxesMPC_grad_eq4);
boxesMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(boxesMPC_C0, boxesMPC_dvaff6, boxesMPC_D1, boxesMPC_dvaff5, boxesMPC_grad_eq5);
boxesMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(boxesMPC_C0, boxesMPC_dvaff7, boxesMPC_D1, boxesMPC_dvaff6, boxesMPC_grad_eq6);
boxesMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(boxesMPC_C0, boxesMPC_dvaff8, boxesMPC_D1, boxesMPC_dvaff7, boxesMPC_grad_eq7);
boxesMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(boxesMPC_C0, boxesMPC_dvaff9, boxesMPC_D1, boxesMPC_dvaff8, boxesMPC_grad_eq8);
boxesMPC_LA_DIAGZERO_MTVM_2_2(boxesMPC_D9, boxesMPC_dvaff9, boxesMPC_grad_eq9);
boxesMPC_LA_VSUB2_38(boxesMPC_rd, boxesMPC_grad_eq, boxesMPC_rd);
boxesMPC_LA_DIAG_FORWARDBACKWARDSUB_4(boxesMPC_Phi0, boxesMPC_rd0, boxesMPC_dzaff0);
boxesMPC_LA_DIAG_FORWARDBACKWARDSUB_4(boxesMPC_Phi1, boxesMPC_rd1, boxesMPC_dzaff1);
boxesMPC_LA_DIAG_FORWARDBACKWARDSUB_4(boxesMPC_Phi2, boxesMPC_rd2, boxesMPC_dzaff2);
boxesMPC_LA_DIAG_FORWARDBACKWARDSUB_4(boxesMPC_Phi3, boxesMPC_rd3, boxesMPC_dzaff3);
boxesMPC_LA_DIAG_FORWARDBACKWARDSUB_4(boxesMPC_Phi4, boxesMPC_rd4, boxesMPC_dzaff4);
boxesMPC_LA_DIAG_FORWARDBACKWARDSUB_4(boxesMPC_Phi5, boxesMPC_rd5, boxesMPC_dzaff5);
boxesMPC_LA_DIAG_FORWARDBACKWARDSUB_4(boxesMPC_Phi6, boxesMPC_rd6, boxesMPC_dzaff6);
boxesMPC_LA_DIAG_FORWARDBACKWARDSUB_4(boxesMPC_Phi7, boxesMPC_rd7, boxesMPC_dzaff7);
boxesMPC_LA_DIAG_FORWARDBACKWARDSUB_4(boxesMPC_Phi8, boxesMPC_rd8, boxesMPC_dzaff8);
boxesMPC_LA_DIAG_FORWARDBACKWARDSUB_2(boxesMPC_Phi9, boxesMPC_rd9, boxesMPC_dzaff9);
boxesMPC_LA_VSUB_INDEXED_4(boxesMPC_dzaff0, boxesMPC_lbIdx0, boxesMPC_rilb0, boxesMPC_dslbaff0);
boxesMPC_LA_VSUB3_4(boxesMPC_llbbyslb0, boxesMPC_dslbaff0, boxesMPC_llb0, boxesMPC_dllbaff0);
boxesMPC_LA_VSUB2_INDEXED_4(boxesMPC_riub0, boxesMPC_dzaff0, boxesMPC_ubIdx0, boxesMPC_dsubaff0);
boxesMPC_LA_VSUB3_4(boxesMPC_lubbysub0, boxesMPC_dsubaff0, boxesMPC_lub0, boxesMPC_dlubaff0);
boxesMPC_LA_VSUB_INDEXED_4(boxesMPC_dzaff1, boxesMPC_lbIdx1, boxesMPC_rilb1, boxesMPC_dslbaff1);
boxesMPC_LA_VSUB3_4(boxesMPC_llbbyslb1, boxesMPC_dslbaff1, boxesMPC_llb1, boxesMPC_dllbaff1);
boxesMPC_LA_VSUB2_INDEXED_4(boxesMPC_riub1, boxesMPC_dzaff1, boxesMPC_ubIdx1, boxesMPC_dsubaff1);
boxesMPC_LA_VSUB3_4(boxesMPC_lubbysub1, boxesMPC_dsubaff1, boxesMPC_lub1, boxesMPC_dlubaff1);
boxesMPC_LA_VSUB_INDEXED_4(boxesMPC_dzaff2, boxesMPC_lbIdx2, boxesMPC_rilb2, boxesMPC_dslbaff2);
boxesMPC_LA_VSUB3_4(boxesMPC_llbbyslb2, boxesMPC_dslbaff2, boxesMPC_llb2, boxesMPC_dllbaff2);
boxesMPC_LA_VSUB2_INDEXED_4(boxesMPC_riub2, boxesMPC_dzaff2, boxesMPC_ubIdx2, boxesMPC_dsubaff2);
boxesMPC_LA_VSUB3_4(boxesMPC_lubbysub2, boxesMPC_dsubaff2, boxesMPC_lub2, boxesMPC_dlubaff2);
boxesMPC_LA_VSUB_INDEXED_4(boxesMPC_dzaff3, boxesMPC_lbIdx3, boxesMPC_rilb3, boxesMPC_dslbaff3);
boxesMPC_LA_VSUB3_4(boxesMPC_llbbyslb3, boxesMPC_dslbaff3, boxesMPC_llb3, boxesMPC_dllbaff3);
boxesMPC_LA_VSUB2_INDEXED_4(boxesMPC_riub3, boxesMPC_dzaff3, boxesMPC_ubIdx3, boxesMPC_dsubaff3);
boxesMPC_LA_VSUB3_4(boxesMPC_lubbysub3, boxesMPC_dsubaff3, boxesMPC_lub3, boxesMPC_dlubaff3);
boxesMPC_LA_VSUB_INDEXED_4(boxesMPC_dzaff4, boxesMPC_lbIdx4, boxesMPC_rilb4, boxesMPC_dslbaff4);
boxesMPC_LA_VSUB3_4(boxesMPC_llbbyslb4, boxesMPC_dslbaff4, boxesMPC_llb4, boxesMPC_dllbaff4);
boxesMPC_LA_VSUB2_INDEXED_4(boxesMPC_riub4, boxesMPC_dzaff4, boxesMPC_ubIdx4, boxesMPC_dsubaff4);
boxesMPC_LA_VSUB3_4(boxesMPC_lubbysub4, boxesMPC_dsubaff4, boxesMPC_lub4, boxesMPC_dlubaff4);
boxesMPC_LA_VSUB_INDEXED_4(boxesMPC_dzaff5, boxesMPC_lbIdx5, boxesMPC_rilb5, boxesMPC_dslbaff5);
boxesMPC_LA_VSUB3_4(boxesMPC_llbbyslb5, boxesMPC_dslbaff5, boxesMPC_llb5, boxesMPC_dllbaff5);
boxesMPC_LA_VSUB2_INDEXED_4(boxesMPC_riub5, boxesMPC_dzaff5, boxesMPC_ubIdx5, boxesMPC_dsubaff5);
boxesMPC_LA_VSUB3_4(boxesMPC_lubbysub5, boxesMPC_dsubaff5, boxesMPC_lub5, boxesMPC_dlubaff5);
boxesMPC_LA_VSUB_INDEXED_4(boxesMPC_dzaff6, boxesMPC_lbIdx6, boxesMPC_rilb6, boxesMPC_dslbaff6);
boxesMPC_LA_VSUB3_4(boxesMPC_llbbyslb6, boxesMPC_dslbaff6, boxesMPC_llb6, boxesMPC_dllbaff6);
boxesMPC_LA_VSUB2_INDEXED_4(boxesMPC_riub6, boxesMPC_dzaff6, boxesMPC_ubIdx6, boxesMPC_dsubaff6);
boxesMPC_LA_VSUB3_4(boxesMPC_lubbysub6, boxesMPC_dsubaff6, boxesMPC_lub6, boxesMPC_dlubaff6);
boxesMPC_LA_VSUB_INDEXED_4(boxesMPC_dzaff7, boxesMPC_lbIdx7, boxesMPC_rilb7, boxesMPC_dslbaff7);
boxesMPC_LA_VSUB3_4(boxesMPC_llbbyslb7, boxesMPC_dslbaff7, boxesMPC_llb7, boxesMPC_dllbaff7);
boxesMPC_LA_VSUB2_INDEXED_4(boxesMPC_riub7, boxesMPC_dzaff7, boxesMPC_ubIdx7, boxesMPC_dsubaff7);
boxesMPC_LA_VSUB3_4(boxesMPC_lubbysub7, boxesMPC_dsubaff7, boxesMPC_lub7, boxesMPC_dlubaff7);
boxesMPC_LA_VSUB_INDEXED_4(boxesMPC_dzaff8, boxesMPC_lbIdx8, boxesMPC_rilb8, boxesMPC_dslbaff8);
boxesMPC_LA_VSUB3_4(boxesMPC_llbbyslb8, boxesMPC_dslbaff8, boxesMPC_llb8, boxesMPC_dllbaff8);
boxesMPC_LA_VSUB2_INDEXED_4(boxesMPC_riub8, boxesMPC_dzaff8, boxesMPC_ubIdx8, boxesMPC_dsubaff8);
boxesMPC_LA_VSUB3_4(boxesMPC_lubbysub8, boxesMPC_dsubaff8, boxesMPC_lub8, boxesMPC_dlubaff8);
boxesMPC_LA_VSUB_INDEXED_2(boxesMPC_dzaff9, boxesMPC_lbIdx9, boxesMPC_rilb9, boxesMPC_dslbaff9);
boxesMPC_LA_VSUB3_2(boxesMPC_llbbyslb9, boxesMPC_dslbaff9, boxesMPC_llb9, boxesMPC_dllbaff9);
boxesMPC_LA_VSUB2_INDEXED_2(boxesMPC_riub9, boxesMPC_dzaff9, boxesMPC_ubIdx9, boxesMPC_dsubaff9);
boxesMPC_LA_VSUB3_2(boxesMPC_lubbysub9, boxesMPC_dsubaff9, boxesMPC_lub9, boxesMPC_dlubaff9);
info->lsit_aff = boxesMPC_LINESEARCH_BACKTRACKING_AFFINE(boxesMPC_l, boxesMPC_s, boxesMPC_dl_aff, boxesMPC_ds_aff, &info->step_aff, &info->mu_aff);
if( info->lsit_aff == boxesMPC_NOPROGRESS ){
exitcode = boxesMPC_NOPROGRESS; break;
}
sigma_3rdroot = info->mu_aff / info->mu;
info->sigma = sigma_3rdroot*sigma_3rdroot*sigma_3rdroot;
musigma = info->mu * info->sigma;
boxesMPC_LA_VSUB5_76(boxesMPC_ds_aff, boxesMPC_dl_aff, musigma, boxesMPC_ccrhs);
boxesMPC_LA_VSUB6_INDEXED_4_4_4(boxesMPC_ccrhsub0, boxesMPC_sub0, boxesMPC_ubIdx0, boxesMPC_ccrhsl0, boxesMPC_slb0, boxesMPC_lbIdx0, boxesMPC_rd0);
boxesMPC_LA_VSUB6_INDEXED_4_4_4(boxesMPC_ccrhsub1, boxesMPC_sub1, boxesMPC_ubIdx1, boxesMPC_ccrhsl1, boxesMPC_slb1, boxesMPC_lbIdx1, boxesMPC_rd1);
boxesMPC_LA_DIAG_FORWARDSUB_4(boxesMPC_Phi0, boxesMPC_rd0, boxesMPC_Lbyrd0);
boxesMPC_LA_DIAG_FORWARDSUB_4(boxesMPC_Phi1, boxesMPC_rd1, boxesMPC_Lbyrd1);
boxesMPC_LA_DIAGZERO_MVM_2(boxesMPC_W0, boxesMPC_Lbyrd0, boxesMPC_beta0);
boxesMPC_LA_DENSE_FORWARDSUB_2(boxesMPC_Ld0, boxesMPC_beta0, boxesMPC_yy0);
boxesMPC_LA_DENSE_DIAGZERO_2MVMADD_2_4_4(boxesMPC_V0, boxesMPC_Lbyrd0, boxesMPC_W1, boxesMPC_Lbyrd1, boxesMPC_beta1);
boxesMPC_LA_DENSE_MVMSUB1_2_2(boxesMPC_Lsd1, boxesMPC_yy0, boxesMPC_beta1, boxesMPC_bmy1);
boxesMPC_LA_DENSE_FORWARDSUB_2(boxesMPC_Ld1, boxesMPC_bmy1, boxesMPC_yy1);
boxesMPC_LA_VSUB6_INDEXED_4_4_4(boxesMPC_ccrhsub2, boxesMPC_sub2, boxesMPC_ubIdx2, boxesMPC_ccrhsl2, boxesMPC_slb2, boxesMPC_lbIdx2, boxesMPC_rd2);
boxesMPC_LA_DIAG_FORWARDSUB_4(boxesMPC_Phi2, boxesMPC_rd2, boxesMPC_Lbyrd2);
boxesMPC_LA_DENSE_DIAGZERO_2MVMADD_2_4_4(boxesMPC_V1, boxesMPC_Lbyrd1, boxesMPC_W2, boxesMPC_Lbyrd2, boxesMPC_beta2);
boxesMPC_LA_DENSE_MVMSUB1_2_2(boxesMPC_Lsd2, boxesMPC_yy1, boxesMPC_beta2, boxesMPC_bmy2);
boxesMPC_LA_DENSE_FORWARDSUB_2(boxesMPC_Ld2, boxesMPC_bmy2, boxesMPC_yy2);
boxesMPC_LA_VSUB6_INDEXED_4_4_4(boxesMPC_ccrhsub3, boxesMPC_sub3, boxesMPC_ubIdx3, boxesMPC_ccrhsl3, boxesMPC_slb3, boxesMPC_lbIdx3, boxesMPC_rd3);
boxesMPC_LA_DIAG_FORWARDSUB_4(boxesMPC_Phi3, boxesMPC_rd3, boxesMPC_Lbyrd3);
boxesMPC_LA_DENSE_DIAGZERO_2MVMADD_2_4_4(boxesMPC_V2, boxesMPC_Lbyrd2, boxesMPC_W3, boxesMPC_Lbyrd3, boxesMPC_beta3);
boxesMPC_LA_DENSE_MVMSUB1_2_2(boxesMPC_Lsd3, boxesMPC_yy2, boxesMPC_beta3, boxesMPC_bmy3);
boxesMPC_LA_DENSE_FORWARDSUB_2(boxesMPC_Ld3, boxesMPC_bmy3, boxesMPC_yy3);
boxesMPC_LA_VSUB6_INDEXED_4_4_4(boxesMPC_ccrhsub4, boxesMPC_sub4, boxesMPC_ubIdx4, boxesMPC_ccrhsl4, boxesMPC_slb4, boxesMPC_lbIdx4, boxesMPC_rd4);
boxesMPC_LA_DIAG_FORWARDSUB_4(boxesMPC_Phi4, boxesMPC_rd4, boxesMPC_Lbyrd4);
boxesMPC_LA_DENSE_DIAGZERO_2MVMADD_2_4_4(boxesMPC_V3, boxesMPC_Lbyrd3, boxesMPC_W4, boxesMPC_Lbyrd4, boxesMPC_beta4);
boxesMPC_LA_DENSE_MVMSUB1_2_2(boxesMPC_Lsd4, boxesMPC_yy3, boxesMPC_beta4, boxesMPC_bmy4);
boxesMPC_LA_DENSE_FORWARDSUB_2(boxesMPC_Ld4, boxesMPC_bmy4, boxesMPC_yy4);
boxesMPC_LA_VSUB6_INDEXED_4_4_4(boxesMPC_ccrhsub5, boxesMPC_sub5, boxesMPC_ubIdx5, boxesMPC_ccrhsl5, boxesMPC_slb5, boxesMPC_lbIdx5, boxesMPC_rd5);
boxesMPC_LA_DIAG_FORWARDSUB_4(boxesMPC_Phi5, boxesMPC_rd5, boxesMPC_Lbyrd5);
boxesMPC_LA_DENSE_DIAGZERO_2MVMADD_2_4_4(boxesMPC_V4, boxesMPC_Lbyrd4, boxesMPC_W5, boxesMPC_Lbyrd5, boxesMPC_beta5);
boxesMPC_LA_DENSE_MVMSUB1_2_2(boxesMPC_Lsd5, boxesMPC_yy4, boxesMPC_beta5, boxesMPC_bmy5);
boxesMPC_LA_DENSE_FORWARDSUB_2(boxesMPC_Ld5, boxesMPC_bmy5, boxesMPC_yy5);
boxesMPC_LA_VSUB6_INDEXED_4_4_4(boxesMPC_ccrhsub6, boxesMPC_sub6, boxesMPC_ubIdx6, boxesMPC_ccrhsl6, boxesMPC_slb6, boxesMPC_lbIdx6, boxesMPC_rd6);
boxesMPC_LA_DIAG_FORWARDSUB_4(boxesMPC_Phi6, boxesMPC_rd6, boxesMPC_Lbyrd6);
boxesMPC_LA_DENSE_DIAGZERO_2MVMADD_2_4_4(boxesMPC_V5, boxesMPC_Lbyrd5, boxesMPC_W6, boxesMPC_Lbyrd6, boxesMPC_beta6);
boxesMPC_LA_DENSE_MVMSUB1_2_2(boxesMPC_Lsd6, boxesMPC_yy5, boxesMPC_beta6, boxesMPC_bmy6);
boxesMPC_LA_DENSE_FORWARDSUB_2(boxesMPC_Ld6, boxesMPC_bmy6, boxesMPC_yy6);
boxesMPC_LA_VSUB6_INDEXED_4_4_4(boxesMPC_ccrhsub7, boxesMPC_sub7, boxesMPC_ubIdx7, boxesMPC_ccrhsl7, boxesMPC_slb7, boxesMPC_lbIdx7, boxesMPC_rd7);
boxesMPC_LA_DIAG_FORWARDSUB_4(boxesMPC_Phi7, boxesMPC_rd7, boxesMPC_Lbyrd7);
boxesMPC_LA_DENSE_DIAGZERO_2MVMADD_2_4_4(boxesMPC_V6, boxesMPC_Lbyrd6, boxesMPC_W7, boxesMPC_Lbyrd7, boxesMPC_beta7);
boxesMPC_LA_DENSE_MVMSUB1_2_2(boxesMPC_Lsd7, boxesMPC_yy6, boxesMPC_beta7, boxesMPC_bmy7);
boxesMPC_LA_DENSE_FORWARDSUB_2(boxesMPC_Ld7, boxesMPC_bmy7, boxesMPC_yy7);
boxesMPC_LA_VSUB6_INDEXED_4_4_4(boxesMPC_ccrhsub8, boxesMPC_sub8, boxesMPC_ubIdx8, boxesMPC_ccrhsl8, boxesMPC_slb8, boxesMPC_lbIdx8, boxesMPC_rd8);
boxesMPC_LA_DIAG_FORWARDSUB_4(boxesMPC_Phi8, boxesMPC_rd8, boxesMPC_Lbyrd8);
boxesMPC_LA_DENSE_DIAGZERO_2MVMADD_2_4_4(boxesMPC_V7, boxesMPC_Lbyrd7, boxesMPC_W8, boxesMPC_Lbyrd8, boxesMPC_beta8);
boxesMPC_LA_DENSE_MVMSUB1_2_2(boxesMPC_Lsd8, boxesMPC_yy7, boxesMPC_beta8, boxesMPC_bmy8);
boxesMPC_LA_DENSE_FORWARDSUB_2(boxesMPC_Ld8, boxesMPC_bmy8, boxesMPC_yy8);
boxesMPC_LA_VSUB6_INDEXED_2_2_2(boxesMPC_ccrhsub9, boxesMPC_sub9, boxesMPC_ubIdx9, boxesMPC_ccrhsl9, boxesMPC_slb9, boxesMPC_lbIdx9, boxesMPC_rd9);
boxesMPC_LA_DIAG_FORWARDSUB_2(boxesMPC_Phi9, boxesMPC_rd9, boxesMPC_Lbyrd9);
boxesMPC_LA_DENSE_DIAGZERO_2MVMADD_2_4_2(boxesMPC_V8, boxesMPC_Lbyrd8, boxesMPC_W9, boxesMPC_Lbyrd9, boxesMPC_beta9);
boxesMPC_LA_DENSE_MVMSUB1_2_2(boxesMPC_Lsd9, boxesMPC_yy8, boxesMPC_beta9, boxesMPC_bmy9);
boxesMPC_LA_DENSE_FORWARDSUB_2(boxesMPC_Ld9, boxesMPC_bmy9, boxesMPC_yy9);
boxesMPC_LA_DENSE_BACKWARDSUB_2(boxesMPC_Ld9, boxesMPC_yy9, boxesMPC_dvcc9);
boxesMPC_LA_DENSE_MTVMSUB_2_2(boxesMPC_Lsd9, boxesMPC_dvcc9, boxesMPC_yy8, boxesMPC_bmy8);
boxesMPC_LA_DENSE_BACKWARDSUB_2(boxesMPC_Ld8, boxesMPC_bmy8, boxesMPC_dvcc8);
boxesMPC_LA_DENSE_MTVMSUB_2_2(boxesMPC_Lsd8, boxesMPC_dvcc8, boxesMPC_yy7, boxesMPC_bmy7);
boxesMPC_LA_DENSE_BACKWARDSUB_2(boxesMPC_Ld7, boxesMPC_bmy7, boxesMPC_dvcc7);
boxesMPC_LA_DENSE_MTVMSUB_2_2(boxesMPC_Lsd7, boxesMPC_dvcc7, boxesMPC_yy6, boxesMPC_bmy6);
boxesMPC_LA_DENSE_BACKWARDSUB_2(boxesMPC_Ld6, boxesMPC_bmy6, boxesMPC_dvcc6);
boxesMPC_LA_DENSE_MTVMSUB_2_2(boxesMPC_Lsd6, boxesMPC_dvcc6, boxesMPC_yy5, boxesMPC_bmy5);
boxesMPC_LA_DENSE_BACKWARDSUB_2(boxesMPC_Ld5, boxesMPC_bmy5, boxesMPC_dvcc5);
boxesMPC_LA_DENSE_MTVMSUB_2_2(boxesMPC_Lsd5, boxesMPC_dvcc5, boxesMPC_yy4, boxesMPC_bmy4);
boxesMPC_LA_DENSE_BACKWARDSUB_2(boxesMPC_Ld4, boxesMPC_bmy4, boxesMPC_dvcc4);
boxesMPC_LA_DENSE_MTVMSUB_2_2(boxesMPC_Lsd4, boxesMPC_dvcc4, boxesMPC_yy3, boxesMPC_bmy3);
boxesMPC_LA_DENSE_BACKWARDSUB_2(boxesMPC_Ld3, boxesMPC_bmy3, boxesMPC_dvcc3);
boxesMPC_LA_DENSE_MTVMSUB_2_2(boxesMPC_Lsd3, boxesMPC_dvcc3, boxesMPC_yy2, boxesMPC_bmy2);
boxesMPC_LA_DENSE_BACKWARDSUB_2(boxesMPC_Ld2, boxesMPC_bmy2, boxesMPC_dvcc2);
boxesMPC_LA_DENSE_MTVMSUB_2_2(boxesMPC_Lsd2, boxesMPC_dvcc2, boxesMPC_yy1, boxesMPC_bmy1);
boxesMPC_LA_DENSE_BACKWARDSUB_2(boxesMPC_Ld1, boxesMPC_bmy1, boxesMPC_dvcc1);
boxesMPC_LA_DENSE_MTVMSUB_2_2(boxesMPC_Lsd1, boxesMPC_dvcc1, boxesMPC_yy0, boxesMPC_bmy0);
boxesMPC_LA_DENSE_BACKWARDSUB_2(boxesMPC_Ld0, boxesMPC_bmy0, boxesMPC_dvcc0);
boxesMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(boxesMPC_C0, boxesMPC_dvcc1, boxesMPC_D0, boxesMPC_dvcc0, boxesMPC_grad_eq0);
boxesMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(boxesMPC_C0, boxesMPC_dvcc2, boxesMPC_D1, boxesMPC_dvcc1, boxesMPC_grad_eq1);
boxesMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(boxesMPC_C0, boxesMPC_dvcc3, boxesMPC_D1, boxesMPC_dvcc2, boxesMPC_grad_eq2);
boxesMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(boxesMPC_C0, boxesMPC_dvcc4, boxesMPC_D1, boxesMPC_dvcc3, boxesMPC_grad_eq3);
boxesMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(boxesMPC_C0, boxesMPC_dvcc5, boxesMPC_D1, boxesMPC_dvcc4, boxesMPC_grad_eq4);
boxesMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(boxesMPC_C0, boxesMPC_dvcc6, boxesMPC_D1, boxesMPC_dvcc5, boxesMPC_grad_eq5);
boxesMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(boxesMPC_C0, boxesMPC_dvcc7, boxesMPC_D1, boxesMPC_dvcc6, boxesMPC_grad_eq6);
boxesMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(boxesMPC_C0, boxesMPC_dvcc8, boxesMPC_D1, boxesMPC_dvcc7, boxesMPC_grad_eq7);
boxesMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(boxesMPC_C0, boxesMPC_dvcc9, boxesMPC_D1, boxesMPC_dvcc8, boxesMPC_grad_eq8);
boxesMPC_LA_DIAGZERO_MTVM_2_2(boxesMPC_D9, boxesMPC_dvcc9, boxesMPC_grad_eq9);
boxesMPC_LA_VSUB_38(boxesMPC_rd, boxesMPC_grad_eq, boxesMPC_rd);
boxesMPC_LA_DIAG_FORWARDBACKWARDSUB_4(boxesMPC_Phi0, boxesMPC_rd0, boxesMPC_dzcc0);
boxesMPC_LA_DIAG_FORWARDBACKWARDSUB_4(boxesMPC_Phi1, boxesMPC_rd1, boxesMPC_dzcc1);
boxesMPC_LA_DIAG_FORWARDBACKWARDSUB_4(boxesMPC_Phi2, boxesMPC_rd2, boxesMPC_dzcc2);
boxesMPC_LA_DIAG_FORWARDBACKWARDSUB_4(boxesMPC_Phi3, boxesMPC_rd3, boxesMPC_dzcc3);
boxesMPC_LA_DIAG_FORWARDBACKWARDSUB_4(boxesMPC_Phi4, boxesMPC_rd4, boxesMPC_dzcc4);
boxesMPC_LA_DIAG_FORWARDBACKWARDSUB_4(boxesMPC_Phi5, boxesMPC_rd5, boxesMPC_dzcc5);
boxesMPC_LA_DIAG_FORWARDBACKWARDSUB_4(boxesMPC_Phi6, boxesMPC_rd6, boxesMPC_dzcc6);
boxesMPC_LA_DIAG_FORWARDBACKWARDSUB_4(boxesMPC_Phi7, boxesMPC_rd7, boxesMPC_dzcc7);
boxesMPC_LA_DIAG_FORWARDBACKWARDSUB_4(boxesMPC_Phi8, boxesMPC_rd8, boxesMPC_dzcc8);
boxesMPC_LA_DIAG_FORWARDBACKWARDSUB_2(boxesMPC_Phi9, boxesMPC_rd9, boxesMPC_dzcc9);
boxesMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_4(boxesMPC_ccrhsl0, boxesMPC_slb0, boxesMPC_llbbyslb0, boxesMPC_dzcc0, boxesMPC_lbIdx0, boxesMPC_dllbcc0);
boxesMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_4(boxesMPC_ccrhsub0, boxesMPC_sub0, boxesMPC_lubbysub0, boxesMPC_dzcc0, boxesMPC_ubIdx0, boxesMPC_dlubcc0);
boxesMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_4(boxesMPC_ccrhsl1, boxesMPC_slb1, boxesMPC_llbbyslb1, boxesMPC_dzcc1, boxesMPC_lbIdx1, boxesMPC_dllbcc1);
boxesMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_4(boxesMPC_ccrhsub1, boxesMPC_sub1, boxesMPC_lubbysub1, boxesMPC_dzcc1, boxesMPC_ubIdx1, boxesMPC_dlubcc1);
boxesMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_4(boxesMPC_ccrhsl2, boxesMPC_slb2, boxesMPC_llbbyslb2, boxesMPC_dzcc2, boxesMPC_lbIdx2, boxesMPC_dllbcc2);
boxesMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_4(boxesMPC_ccrhsub2, boxesMPC_sub2, boxesMPC_lubbysub2, boxesMPC_dzcc2, boxesMPC_ubIdx2, boxesMPC_dlubcc2);
boxesMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_4(boxesMPC_ccrhsl3, boxesMPC_slb3, boxesMPC_llbbyslb3, boxesMPC_dzcc3, boxesMPC_lbIdx3, boxesMPC_dllbcc3);
boxesMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_4(boxesMPC_ccrhsub3, boxesMPC_sub3, boxesMPC_lubbysub3, boxesMPC_dzcc3, boxesMPC_ubIdx3, boxesMPC_dlubcc3);
boxesMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_4(boxesMPC_ccrhsl4, boxesMPC_slb4, boxesMPC_llbbyslb4, boxesMPC_dzcc4, boxesMPC_lbIdx4, boxesMPC_dllbcc4);
boxesMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_4(boxesMPC_ccrhsub4, boxesMPC_sub4, boxesMPC_lubbysub4, boxesMPC_dzcc4, boxesMPC_ubIdx4, boxesMPC_dlubcc4);
boxesMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_4(boxesMPC_ccrhsl5, boxesMPC_slb5, boxesMPC_llbbyslb5, boxesMPC_dzcc5, boxesMPC_lbIdx5, boxesMPC_dllbcc5);
boxesMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_4(boxesMPC_ccrhsub5, boxesMPC_sub5, boxesMPC_lubbysub5, boxesMPC_dzcc5, boxesMPC_ubIdx5, boxesMPC_dlubcc5);
boxesMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_4(boxesMPC_ccrhsl6, boxesMPC_slb6, boxesMPC_llbbyslb6, boxesMPC_dzcc6, boxesMPC_lbIdx6, boxesMPC_dllbcc6);
boxesMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_4(boxesMPC_ccrhsub6, boxesMPC_sub6, boxesMPC_lubbysub6, boxesMPC_dzcc6, boxesMPC_ubIdx6, boxesMPC_dlubcc6);
boxesMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_4(boxesMPC_ccrhsl7, boxesMPC_slb7, boxesMPC_llbbyslb7, boxesMPC_dzcc7, boxesMPC_lbIdx7, boxesMPC_dllbcc7);
boxesMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_4(boxesMPC_ccrhsub7, boxesMPC_sub7, boxesMPC_lubbysub7, boxesMPC_dzcc7, boxesMPC_ubIdx7, boxesMPC_dlubcc7);
boxesMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_4(boxesMPC_ccrhsl8, boxesMPC_slb8, boxesMPC_llbbyslb8, boxesMPC_dzcc8, boxesMPC_lbIdx8, boxesMPC_dllbcc8);
boxesMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_4(boxesMPC_ccrhsub8, boxesMPC_sub8, boxesMPC_lubbysub8, boxesMPC_dzcc8, boxesMPC_ubIdx8, boxesMPC_dlubcc8);
boxesMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_2(boxesMPC_ccrhsl9, boxesMPC_slb9, boxesMPC_llbbyslb9, boxesMPC_dzcc9, boxesMPC_lbIdx9, boxesMPC_dllbcc9);
boxesMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_2(boxesMPC_ccrhsub9, boxesMPC_sub9, boxesMPC_lubbysub9, boxesMPC_dzcc9, boxesMPC_ubIdx9, boxesMPC_dlubcc9);
boxesMPC_LA_VSUB7_76(boxesMPC_l, boxesMPC_ccrhs, boxesMPC_s, boxesMPC_dl_cc, boxesMPC_ds_cc);
boxesMPC_LA_VADD_38(boxesMPC_dz_cc, boxesMPC_dz_aff);
boxesMPC_LA_VADD_20(boxesMPC_dv_cc, boxesMPC_dv_aff);
boxesMPC_LA_VADD_76(boxesMPC_dl_cc, boxesMPC_dl_aff);
boxesMPC_LA_VADD_76(boxesMPC_ds_cc, boxesMPC_ds_aff);
info->lsit_cc = boxesMPC_LINESEARCH_BACKTRACKING_COMBINED(boxesMPC_z, boxesMPC_v, boxesMPC_l, boxesMPC_s, boxesMPC_dz_cc, boxesMPC_dv_cc, boxesMPC_dl_cc, boxesMPC_ds_cc, &info->step_cc, &info->mu);
if( info->lsit_cc == boxesMPC_NOPROGRESS ){
exitcode = boxesMPC_NOPROGRESS; break;
}
info->it++;
}
output->z1[0] = boxesMPC_z0[0];
output->z1[1] = boxesMPC_z0[1];
output->z1[2] = boxesMPC_z0[2];
output->z1[3] = boxesMPC_z0[3];
output->z2[0] = boxesMPC_z1[0];
output->z2[1] = boxesMPC_z1[1];
output->z2[2] = boxesMPC_z1[2];
output->z2[3] = boxesMPC_z1[3];
output->z3[0] = boxesMPC_z2[0];
output->z3[1] = boxesMPC_z2[1];
output->z3[2] = boxesMPC_z2[2];
output->z3[3] = boxesMPC_z2[3];
output->z4[0] = boxesMPC_z3[0];
output->z4[1] = boxesMPC_z3[1];
output->z4[2] = boxesMPC_z3[2];
output->z4[3] = boxesMPC_z3[3];
output->z5[0] = boxesMPC_z4[0];
output->z5[1] = boxesMPC_z4[1];
output->z5[2] = boxesMPC_z4[2];
output->z5[3] = boxesMPC_z4[3];
output->z6[0] = boxesMPC_z5[0];
output->z6[1] = boxesMPC_z5[1];
output->z6[2] = boxesMPC_z5[2];
output->z6[3] = boxesMPC_z5[3];
output->z7[0] = boxesMPC_z6[0];
output->z7[1] = boxesMPC_z6[1];
output->z7[2] = boxesMPC_z6[2];
output->z7[3] = boxesMPC_z6[3];
output->z8[0] = boxesMPC_z7[0];
output->z8[1] = boxesMPC_z7[1];
output->z8[2] = boxesMPC_z7[2];
output->z8[3] = boxesMPC_z7[3];
output->z9[0] = boxesMPC_z8[0];
output->z9[1] = boxesMPC_z8[1];
output->z9[2] = boxesMPC_z8[2];
output->z9[3] = boxesMPC_z8[3];
output->z10[0] = boxesMPC_z9[0];
output->z10[1] = boxesMPC_z9[1];

#if boxesMPC_SET_TIMING == 1
info->solvetime = boxesMPC_toc(&solvertimer);
#if boxesMPC_SET_PRINTLEVEL > 0 && boxesMPC_SET_TIMING == 1
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
