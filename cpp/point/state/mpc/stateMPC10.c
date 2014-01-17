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

#include "stateMPC.h"

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
void stateMPC_LA_INITIALIZEVECTOR_38(stateMPC_FLOAT* vec, stateMPC_FLOAT value)
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
void stateMPC_LA_INITIALIZEVECTOR_20(stateMPC_FLOAT* vec, stateMPC_FLOAT value)
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
void stateMPC_LA_INITIALIZEVECTOR_76(stateMPC_FLOAT* vec, stateMPC_FLOAT value)
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
void stateMPC_LA_DOTACC_76(stateMPC_FLOAT *x, stateMPC_FLOAT *y, stateMPC_FLOAT *z)
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
 * of length 38.
 */
void stateMPC_LA_VVADD3_38(stateMPC_FLOAT *u, stateMPC_FLOAT *v, stateMPC_FLOAT *w, stateMPC_FLOAT *z)
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
 * Vector subtraction z = -x - y for vectors of length 38.
 */
void stateMPC_LA_VSUB2_38(stateMPC_FLOAT *x, stateMPC_FLOAT *y, stateMPC_FLOAT *z)
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
            mya *= stateMPC_SET_LS_SCALE_AFF;
            if( mya < stateMPC_SET_LS_MINSTEP ){
                return stateMPC_NOPROGRESS;
            }
        }
    }
    
    /* return new values and iteration counter */
    *a = mya;
    *mu_aff = mymu / (stateMPC_FLOAT)76;
    return lsIt;
}


/*
 * Vector subtraction x = u.*v - a where a is a scalar
*  and x,u,v are vectors of length 76.
 */
void stateMPC_LA_VSUB5_76(stateMPC_FLOAT *u, stateMPC_FLOAT *v, stateMPC_FLOAT a, stateMPC_FLOAT *x)
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
 * Vector subtraction z = x - y for vectors of length 38.
 */
void stateMPC_LA_VSUB_38(stateMPC_FLOAT *x, stateMPC_FLOAT *y, stateMPC_FLOAT *z)
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
 * Computes ds = -l.\(r + s.*dl) for vectors of length 76.
 */
void stateMPC_LA_VSUB7_76(stateMPC_FLOAT *l, stateMPC_FLOAT *r, stateMPC_FLOAT *s, stateMPC_FLOAT *dl, stateMPC_FLOAT *ds)
{
	int i;
	for( i=0; i<76; i++){
		ds[i] = -(r[i] + s[i]*dl[i])/l[i];
	}
}


/*
 * Vector addition x = x + y for vectors of length 38.
 */
void stateMPC_LA_VADD_38(stateMPC_FLOAT *x, stateMPC_FLOAT *y)
{
	int i;
	for( i=0; i<38; i++){
		x[i] += y[i];
	}
}


/*
 * Vector addition x = x + y for vectors of length 20.
 */
void stateMPC_LA_VADD_20(stateMPC_FLOAT *x, stateMPC_FLOAT *y)
{
	int i;
	for( i=0; i<20; i++){
		x[i] += y[i];
	}
}


/*
 * Vector addition x = x + y for vectors of length 76.
 */
void stateMPC_LA_VADD_76(stateMPC_FLOAT *x, stateMPC_FLOAT *y)
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
int stateMPC_LINESEARCH_BACKTRACKING_COMBINED(stateMPC_FLOAT *z, stateMPC_FLOAT *v, stateMPC_FLOAT *l, stateMPC_FLOAT *s, stateMPC_FLOAT *dz, stateMPC_FLOAT *dv, stateMPC_FLOAT *dl, stateMPC_FLOAT *ds, stateMPC_FLOAT *a, stateMPC_FLOAT *mu)
{
    int i, lsIt=1;       
    stateMPC_FLOAT dltemp;
    stateMPC_FLOAT dstemp;    
    stateMPC_FLOAT a_gamma;
            
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
            *a *= stateMPC_SET_LS_SCALE;
            if( *a < stateMPC_SET_LS_MINSTEP ){
                return stateMPC_NOPROGRESS;
            }
        }
    }
    
    /* update variables with safety margin */
    a_gamma = (*a)*stateMPC_SET_LS_MAXSTEP;
    
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
    *mu /= (stateMPC_FLOAT)76;
    return lsIt;
}




/* VARIABLE DEFINITIONS ------------------------------------------------ */
stateMPC_FLOAT stateMPC_z[38];
stateMPC_FLOAT stateMPC_v[20];
stateMPC_FLOAT stateMPC_dz_aff[38];
stateMPC_FLOAT stateMPC_dv_aff[20];
stateMPC_FLOAT stateMPC_grad_cost[38];
stateMPC_FLOAT stateMPC_grad_eq[38];
stateMPC_FLOAT stateMPC_rd[38];
stateMPC_FLOAT stateMPC_l[76];
stateMPC_FLOAT stateMPC_s[76];
stateMPC_FLOAT stateMPC_lbys[76];
stateMPC_FLOAT stateMPC_dl_aff[76];
stateMPC_FLOAT stateMPC_ds_aff[76];
stateMPC_FLOAT stateMPC_dz_cc[38];
stateMPC_FLOAT stateMPC_dv_cc[20];
stateMPC_FLOAT stateMPC_dl_cc[76];
stateMPC_FLOAT stateMPC_ds_cc[76];
stateMPC_FLOAT stateMPC_ccrhs[76];
stateMPC_FLOAT stateMPC_grad_ineq[38];
stateMPC_FLOAT* stateMPC_z0 = stateMPC_z + 0;
stateMPC_FLOAT* stateMPC_dzaff0 = stateMPC_dz_aff + 0;
stateMPC_FLOAT* stateMPC_dzcc0 = stateMPC_dz_cc + 0;
stateMPC_FLOAT* stateMPC_rd0 = stateMPC_rd + 0;
stateMPC_FLOAT stateMPC_Lbyrd0[4];
stateMPC_FLOAT* stateMPC_grad_cost0 = stateMPC_grad_cost + 0;
stateMPC_FLOAT* stateMPC_grad_eq0 = stateMPC_grad_eq + 0;
stateMPC_FLOAT* stateMPC_grad_ineq0 = stateMPC_grad_ineq + 0;
stateMPC_FLOAT stateMPC_ctv0[4];
stateMPC_FLOAT* stateMPC_v0 = stateMPC_v + 0;
stateMPC_FLOAT stateMPC_re0[2];
stateMPC_FLOAT stateMPC_beta0[2];
stateMPC_FLOAT stateMPC_betacc0[2];
stateMPC_FLOAT* stateMPC_dvaff0 = stateMPC_dv_aff + 0;
stateMPC_FLOAT* stateMPC_dvcc0 = stateMPC_dv_cc + 0;
stateMPC_FLOAT stateMPC_V0[8];
stateMPC_FLOAT stateMPC_Yd0[3];
stateMPC_FLOAT stateMPC_Ld0[3];
stateMPC_FLOAT stateMPC_yy0[2];
stateMPC_FLOAT stateMPC_bmy0[2];
int stateMPC_lbIdx0[4] = {0, 1, 2, 3};
stateMPC_FLOAT* stateMPC_llb0 = stateMPC_l + 0;
stateMPC_FLOAT* stateMPC_slb0 = stateMPC_s + 0;
stateMPC_FLOAT* stateMPC_llbbyslb0 = stateMPC_lbys + 0;
stateMPC_FLOAT stateMPC_rilb0[4];
stateMPC_FLOAT* stateMPC_dllbaff0 = stateMPC_dl_aff + 0;
stateMPC_FLOAT* stateMPC_dslbaff0 = stateMPC_ds_aff + 0;
stateMPC_FLOAT* stateMPC_dllbcc0 = stateMPC_dl_cc + 0;
stateMPC_FLOAT* stateMPC_dslbcc0 = stateMPC_ds_cc + 0;
stateMPC_FLOAT* stateMPC_ccrhsl0 = stateMPC_ccrhs + 0;
int stateMPC_ubIdx0[4] = {0, 1, 2, 3};
stateMPC_FLOAT* stateMPC_lub0 = stateMPC_l + 4;
stateMPC_FLOAT* stateMPC_sub0 = stateMPC_s + 4;
stateMPC_FLOAT* stateMPC_lubbysub0 = stateMPC_lbys + 4;
stateMPC_FLOAT stateMPC_riub0[4];
stateMPC_FLOAT* stateMPC_dlubaff0 = stateMPC_dl_aff + 4;
stateMPC_FLOAT* stateMPC_dsubaff0 = stateMPC_ds_aff + 4;
stateMPC_FLOAT* stateMPC_dlubcc0 = stateMPC_dl_cc + 4;
stateMPC_FLOAT* stateMPC_dsubcc0 = stateMPC_ds_cc + 4;
stateMPC_FLOAT* stateMPC_ccrhsub0 = stateMPC_ccrhs + 4;
stateMPC_FLOAT stateMPC_Phi0[4];
stateMPC_FLOAT stateMPC_D0[4] = {1.0000000000000000E+000, 
1.0000000000000000E+000};
stateMPC_FLOAT stateMPC_W0[4];
stateMPC_FLOAT* stateMPC_z1 = stateMPC_z + 4;
stateMPC_FLOAT* stateMPC_dzaff1 = stateMPC_dz_aff + 4;
stateMPC_FLOAT* stateMPC_dzcc1 = stateMPC_dz_cc + 4;
stateMPC_FLOAT* stateMPC_rd1 = stateMPC_rd + 4;
stateMPC_FLOAT stateMPC_Lbyrd1[4];
stateMPC_FLOAT* stateMPC_grad_cost1 = stateMPC_grad_cost + 4;
stateMPC_FLOAT* stateMPC_grad_eq1 = stateMPC_grad_eq + 4;
stateMPC_FLOAT* stateMPC_grad_ineq1 = stateMPC_grad_ineq + 4;
stateMPC_FLOAT stateMPC_ctv1[4];
stateMPC_FLOAT* stateMPC_v1 = stateMPC_v + 2;
stateMPC_FLOAT stateMPC_re1[2];
stateMPC_FLOAT stateMPC_beta1[2];
stateMPC_FLOAT stateMPC_betacc1[2];
stateMPC_FLOAT* stateMPC_dvaff1 = stateMPC_dv_aff + 2;
stateMPC_FLOAT* stateMPC_dvcc1 = stateMPC_dv_cc + 2;
stateMPC_FLOAT stateMPC_V1[8];
stateMPC_FLOAT stateMPC_Yd1[3];
stateMPC_FLOAT stateMPC_Ld1[3];
stateMPC_FLOAT stateMPC_yy1[2];
stateMPC_FLOAT stateMPC_bmy1[2];
int stateMPC_lbIdx1[4] = {0, 1, 2, 3};
stateMPC_FLOAT* stateMPC_llb1 = stateMPC_l + 8;
stateMPC_FLOAT* stateMPC_slb1 = stateMPC_s + 8;
stateMPC_FLOAT* stateMPC_llbbyslb1 = stateMPC_lbys + 8;
stateMPC_FLOAT stateMPC_rilb1[4];
stateMPC_FLOAT* stateMPC_dllbaff1 = stateMPC_dl_aff + 8;
stateMPC_FLOAT* stateMPC_dslbaff1 = stateMPC_ds_aff + 8;
stateMPC_FLOAT* stateMPC_dllbcc1 = stateMPC_dl_cc + 8;
stateMPC_FLOAT* stateMPC_dslbcc1 = stateMPC_ds_cc + 8;
stateMPC_FLOAT* stateMPC_ccrhsl1 = stateMPC_ccrhs + 8;
int stateMPC_ubIdx1[4] = {0, 1, 2, 3};
stateMPC_FLOAT* stateMPC_lub1 = stateMPC_l + 12;
stateMPC_FLOAT* stateMPC_sub1 = stateMPC_s + 12;
stateMPC_FLOAT* stateMPC_lubbysub1 = stateMPC_lbys + 12;
stateMPC_FLOAT stateMPC_riub1[4];
stateMPC_FLOAT* stateMPC_dlubaff1 = stateMPC_dl_aff + 12;
stateMPC_FLOAT* stateMPC_dsubaff1 = stateMPC_ds_aff + 12;
stateMPC_FLOAT* stateMPC_dlubcc1 = stateMPC_dl_cc + 12;
stateMPC_FLOAT* stateMPC_dsubcc1 = stateMPC_ds_cc + 12;
stateMPC_FLOAT* stateMPC_ccrhsub1 = stateMPC_ccrhs + 12;
stateMPC_FLOAT stateMPC_Phi1[4];
stateMPC_FLOAT stateMPC_D1[4] = {-1.0000000000000000E+000, 
-1.0000000000000000E+000};
stateMPC_FLOAT stateMPC_W1[4];
stateMPC_FLOAT stateMPC_Ysd1[4];
stateMPC_FLOAT stateMPC_Lsd1[4];
stateMPC_FLOAT* stateMPC_z2 = stateMPC_z + 8;
stateMPC_FLOAT* stateMPC_dzaff2 = stateMPC_dz_aff + 8;
stateMPC_FLOAT* stateMPC_dzcc2 = stateMPC_dz_cc + 8;
stateMPC_FLOAT* stateMPC_rd2 = stateMPC_rd + 8;
stateMPC_FLOAT stateMPC_Lbyrd2[4];
stateMPC_FLOAT* stateMPC_grad_cost2 = stateMPC_grad_cost + 8;
stateMPC_FLOAT* stateMPC_grad_eq2 = stateMPC_grad_eq + 8;
stateMPC_FLOAT* stateMPC_grad_ineq2 = stateMPC_grad_ineq + 8;
stateMPC_FLOAT stateMPC_ctv2[4];
stateMPC_FLOAT* stateMPC_v2 = stateMPC_v + 4;
stateMPC_FLOAT stateMPC_re2[2];
stateMPC_FLOAT stateMPC_beta2[2];
stateMPC_FLOAT stateMPC_betacc2[2];
stateMPC_FLOAT* stateMPC_dvaff2 = stateMPC_dv_aff + 4;
stateMPC_FLOAT* stateMPC_dvcc2 = stateMPC_dv_cc + 4;
stateMPC_FLOAT stateMPC_V2[8];
stateMPC_FLOAT stateMPC_Yd2[3];
stateMPC_FLOAT stateMPC_Ld2[3];
stateMPC_FLOAT stateMPC_yy2[2];
stateMPC_FLOAT stateMPC_bmy2[2];
int stateMPC_lbIdx2[4] = {0, 1, 2, 3};
stateMPC_FLOAT* stateMPC_llb2 = stateMPC_l + 16;
stateMPC_FLOAT* stateMPC_slb2 = stateMPC_s + 16;
stateMPC_FLOAT* stateMPC_llbbyslb2 = stateMPC_lbys + 16;
stateMPC_FLOAT stateMPC_rilb2[4];
stateMPC_FLOAT* stateMPC_dllbaff2 = stateMPC_dl_aff + 16;
stateMPC_FLOAT* stateMPC_dslbaff2 = stateMPC_ds_aff + 16;
stateMPC_FLOAT* stateMPC_dllbcc2 = stateMPC_dl_cc + 16;
stateMPC_FLOAT* stateMPC_dslbcc2 = stateMPC_ds_cc + 16;
stateMPC_FLOAT* stateMPC_ccrhsl2 = stateMPC_ccrhs + 16;
int stateMPC_ubIdx2[4] = {0, 1, 2, 3};
stateMPC_FLOAT* stateMPC_lub2 = stateMPC_l + 20;
stateMPC_FLOAT* stateMPC_sub2 = stateMPC_s + 20;
stateMPC_FLOAT* stateMPC_lubbysub2 = stateMPC_lbys + 20;
stateMPC_FLOAT stateMPC_riub2[4];
stateMPC_FLOAT* stateMPC_dlubaff2 = stateMPC_dl_aff + 20;
stateMPC_FLOAT* stateMPC_dsubaff2 = stateMPC_ds_aff + 20;
stateMPC_FLOAT* stateMPC_dlubcc2 = stateMPC_dl_cc + 20;
stateMPC_FLOAT* stateMPC_dsubcc2 = stateMPC_ds_cc + 20;
stateMPC_FLOAT* stateMPC_ccrhsub2 = stateMPC_ccrhs + 20;
stateMPC_FLOAT stateMPC_Phi2[4];
stateMPC_FLOAT stateMPC_W2[4];
stateMPC_FLOAT stateMPC_Ysd2[4];
stateMPC_FLOAT stateMPC_Lsd2[4];
stateMPC_FLOAT* stateMPC_z3 = stateMPC_z + 12;
stateMPC_FLOAT* stateMPC_dzaff3 = stateMPC_dz_aff + 12;
stateMPC_FLOAT* stateMPC_dzcc3 = stateMPC_dz_cc + 12;
stateMPC_FLOAT* stateMPC_rd3 = stateMPC_rd + 12;
stateMPC_FLOAT stateMPC_Lbyrd3[4];
stateMPC_FLOAT* stateMPC_grad_cost3 = stateMPC_grad_cost + 12;
stateMPC_FLOAT* stateMPC_grad_eq3 = stateMPC_grad_eq + 12;
stateMPC_FLOAT* stateMPC_grad_ineq3 = stateMPC_grad_ineq + 12;
stateMPC_FLOAT stateMPC_ctv3[4];
stateMPC_FLOAT* stateMPC_v3 = stateMPC_v + 6;
stateMPC_FLOAT stateMPC_re3[2];
stateMPC_FLOAT stateMPC_beta3[2];
stateMPC_FLOAT stateMPC_betacc3[2];
stateMPC_FLOAT* stateMPC_dvaff3 = stateMPC_dv_aff + 6;
stateMPC_FLOAT* stateMPC_dvcc3 = stateMPC_dv_cc + 6;
stateMPC_FLOAT stateMPC_V3[8];
stateMPC_FLOAT stateMPC_Yd3[3];
stateMPC_FLOAT stateMPC_Ld3[3];
stateMPC_FLOAT stateMPC_yy3[2];
stateMPC_FLOAT stateMPC_bmy3[2];
int stateMPC_lbIdx3[4] = {0, 1, 2, 3};
stateMPC_FLOAT* stateMPC_llb3 = stateMPC_l + 24;
stateMPC_FLOAT* stateMPC_slb3 = stateMPC_s + 24;
stateMPC_FLOAT* stateMPC_llbbyslb3 = stateMPC_lbys + 24;
stateMPC_FLOAT stateMPC_rilb3[4];
stateMPC_FLOAT* stateMPC_dllbaff3 = stateMPC_dl_aff + 24;
stateMPC_FLOAT* stateMPC_dslbaff3 = stateMPC_ds_aff + 24;
stateMPC_FLOAT* stateMPC_dllbcc3 = stateMPC_dl_cc + 24;
stateMPC_FLOAT* stateMPC_dslbcc3 = stateMPC_ds_cc + 24;
stateMPC_FLOAT* stateMPC_ccrhsl3 = stateMPC_ccrhs + 24;
int stateMPC_ubIdx3[4] = {0, 1, 2, 3};
stateMPC_FLOAT* stateMPC_lub3 = stateMPC_l + 28;
stateMPC_FLOAT* stateMPC_sub3 = stateMPC_s + 28;
stateMPC_FLOAT* stateMPC_lubbysub3 = stateMPC_lbys + 28;
stateMPC_FLOAT stateMPC_riub3[4];
stateMPC_FLOAT* stateMPC_dlubaff3 = stateMPC_dl_aff + 28;
stateMPC_FLOAT* stateMPC_dsubaff3 = stateMPC_ds_aff + 28;
stateMPC_FLOAT* stateMPC_dlubcc3 = stateMPC_dl_cc + 28;
stateMPC_FLOAT* stateMPC_dsubcc3 = stateMPC_ds_cc + 28;
stateMPC_FLOAT* stateMPC_ccrhsub3 = stateMPC_ccrhs + 28;
stateMPC_FLOAT stateMPC_Phi3[4];
stateMPC_FLOAT stateMPC_W3[4];
stateMPC_FLOAT stateMPC_Ysd3[4];
stateMPC_FLOAT stateMPC_Lsd3[4];
stateMPC_FLOAT* stateMPC_z4 = stateMPC_z + 16;
stateMPC_FLOAT* stateMPC_dzaff4 = stateMPC_dz_aff + 16;
stateMPC_FLOAT* stateMPC_dzcc4 = stateMPC_dz_cc + 16;
stateMPC_FLOAT* stateMPC_rd4 = stateMPC_rd + 16;
stateMPC_FLOAT stateMPC_Lbyrd4[4];
stateMPC_FLOAT* stateMPC_grad_cost4 = stateMPC_grad_cost + 16;
stateMPC_FLOAT* stateMPC_grad_eq4 = stateMPC_grad_eq + 16;
stateMPC_FLOAT* stateMPC_grad_ineq4 = stateMPC_grad_ineq + 16;
stateMPC_FLOAT stateMPC_ctv4[4];
stateMPC_FLOAT* stateMPC_v4 = stateMPC_v + 8;
stateMPC_FLOAT stateMPC_re4[2];
stateMPC_FLOAT stateMPC_beta4[2];
stateMPC_FLOAT stateMPC_betacc4[2];
stateMPC_FLOAT* stateMPC_dvaff4 = stateMPC_dv_aff + 8;
stateMPC_FLOAT* stateMPC_dvcc4 = stateMPC_dv_cc + 8;
stateMPC_FLOAT stateMPC_V4[8];
stateMPC_FLOAT stateMPC_Yd4[3];
stateMPC_FLOAT stateMPC_Ld4[3];
stateMPC_FLOAT stateMPC_yy4[2];
stateMPC_FLOAT stateMPC_bmy4[2];
int stateMPC_lbIdx4[4] = {0, 1, 2, 3};
stateMPC_FLOAT* stateMPC_llb4 = stateMPC_l + 32;
stateMPC_FLOAT* stateMPC_slb4 = stateMPC_s + 32;
stateMPC_FLOAT* stateMPC_llbbyslb4 = stateMPC_lbys + 32;
stateMPC_FLOAT stateMPC_rilb4[4];
stateMPC_FLOAT* stateMPC_dllbaff4 = stateMPC_dl_aff + 32;
stateMPC_FLOAT* stateMPC_dslbaff4 = stateMPC_ds_aff + 32;
stateMPC_FLOAT* stateMPC_dllbcc4 = stateMPC_dl_cc + 32;
stateMPC_FLOAT* stateMPC_dslbcc4 = stateMPC_ds_cc + 32;
stateMPC_FLOAT* stateMPC_ccrhsl4 = stateMPC_ccrhs + 32;
int stateMPC_ubIdx4[4] = {0, 1, 2, 3};
stateMPC_FLOAT* stateMPC_lub4 = stateMPC_l + 36;
stateMPC_FLOAT* stateMPC_sub4 = stateMPC_s + 36;
stateMPC_FLOAT* stateMPC_lubbysub4 = stateMPC_lbys + 36;
stateMPC_FLOAT stateMPC_riub4[4];
stateMPC_FLOAT* stateMPC_dlubaff4 = stateMPC_dl_aff + 36;
stateMPC_FLOAT* stateMPC_dsubaff4 = stateMPC_ds_aff + 36;
stateMPC_FLOAT* stateMPC_dlubcc4 = stateMPC_dl_cc + 36;
stateMPC_FLOAT* stateMPC_dsubcc4 = stateMPC_ds_cc + 36;
stateMPC_FLOAT* stateMPC_ccrhsub4 = stateMPC_ccrhs + 36;
stateMPC_FLOAT stateMPC_Phi4[4];
stateMPC_FLOAT stateMPC_W4[4];
stateMPC_FLOAT stateMPC_Ysd4[4];
stateMPC_FLOAT stateMPC_Lsd4[4];
stateMPC_FLOAT* stateMPC_z5 = stateMPC_z + 20;
stateMPC_FLOAT* stateMPC_dzaff5 = stateMPC_dz_aff + 20;
stateMPC_FLOAT* stateMPC_dzcc5 = stateMPC_dz_cc + 20;
stateMPC_FLOAT* stateMPC_rd5 = stateMPC_rd + 20;
stateMPC_FLOAT stateMPC_Lbyrd5[4];
stateMPC_FLOAT* stateMPC_grad_cost5 = stateMPC_grad_cost + 20;
stateMPC_FLOAT* stateMPC_grad_eq5 = stateMPC_grad_eq + 20;
stateMPC_FLOAT* stateMPC_grad_ineq5 = stateMPC_grad_ineq + 20;
stateMPC_FLOAT stateMPC_ctv5[4];
stateMPC_FLOAT* stateMPC_v5 = stateMPC_v + 10;
stateMPC_FLOAT stateMPC_re5[2];
stateMPC_FLOAT stateMPC_beta5[2];
stateMPC_FLOAT stateMPC_betacc5[2];
stateMPC_FLOAT* stateMPC_dvaff5 = stateMPC_dv_aff + 10;
stateMPC_FLOAT* stateMPC_dvcc5 = stateMPC_dv_cc + 10;
stateMPC_FLOAT stateMPC_V5[8];
stateMPC_FLOAT stateMPC_Yd5[3];
stateMPC_FLOAT stateMPC_Ld5[3];
stateMPC_FLOAT stateMPC_yy5[2];
stateMPC_FLOAT stateMPC_bmy5[2];
int stateMPC_lbIdx5[4] = {0, 1, 2, 3};
stateMPC_FLOAT* stateMPC_llb5 = stateMPC_l + 40;
stateMPC_FLOAT* stateMPC_slb5 = stateMPC_s + 40;
stateMPC_FLOAT* stateMPC_llbbyslb5 = stateMPC_lbys + 40;
stateMPC_FLOAT stateMPC_rilb5[4];
stateMPC_FLOAT* stateMPC_dllbaff5 = stateMPC_dl_aff + 40;
stateMPC_FLOAT* stateMPC_dslbaff5 = stateMPC_ds_aff + 40;
stateMPC_FLOAT* stateMPC_dllbcc5 = stateMPC_dl_cc + 40;
stateMPC_FLOAT* stateMPC_dslbcc5 = stateMPC_ds_cc + 40;
stateMPC_FLOAT* stateMPC_ccrhsl5 = stateMPC_ccrhs + 40;
int stateMPC_ubIdx5[4] = {0, 1, 2, 3};
stateMPC_FLOAT* stateMPC_lub5 = stateMPC_l + 44;
stateMPC_FLOAT* stateMPC_sub5 = stateMPC_s + 44;
stateMPC_FLOAT* stateMPC_lubbysub5 = stateMPC_lbys + 44;
stateMPC_FLOAT stateMPC_riub5[4];
stateMPC_FLOAT* stateMPC_dlubaff5 = stateMPC_dl_aff + 44;
stateMPC_FLOAT* stateMPC_dsubaff5 = stateMPC_ds_aff + 44;
stateMPC_FLOAT* stateMPC_dlubcc5 = stateMPC_dl_cc + 44;
stateMPC_FLOAT* stateMPC_dsubcc5 = stateMPC_ds_cc + 44;
stateMPC_FLOAT* stateMPC_ccrhsub5 = stateMPC_ccrhs + 44;
stateMPC_FLOAT stateMPC_Phi5[4];
stateMPC_FLOAT stateMPC_W5[4];
stateMPC_FLOAT stateMPC_Ysd5[4];
stateMPC_FLOAT stateMPC_Lsd5[4];
stateMPC_FLOAT* stateMPC_z6 = stateMPC_z + 24;
stateMPC_FLOAT* stateMPC_dzaff6 = stateMPC_dz_aff + 24;
stateMPC_FLOAT* stateMPC_dzcc6 = stateMPC_dz_cc + 24;
stateMPC_FLOAT* stateMPC_rd6 = stateMPC_rd + 24;
stateMPC_FLOAT stateMPC_Lbyrd6[4];
stateMPC_FLOAT* stateMPC_grad_cost6 = stateMPC_grad_cost + 24;
stateMPC_FLOAT* stateMPC_grad_eq6 = stateMPC_grad_eq + 24;
stateMPC_FLOAT* stateMPC_grad_ineq6 = stateMPC_grad_ineq + 24;
stateMPC_FLOAT stateMPC_ctv6[4];
stateMPC_FLOAT* stateMPC_v6 = stateMPC_v + 12;
stateMPC_FLOAT stateMPC_re6[2];
stateMPC_FLOAT stateMPC_beta6[2];
stateMPC_FLOAT stateMPC_betacc6[2];
stateMPC_FLOAT* stateMPC_dvaff6 = stateMPC_dv_aff + 12;
stateMPC_FLOAT* stateMPC_dvcc6 = stateMPC_dv_cc + 12;
stateMPC_FLOAT stateMPC_V6[8];
stateMPC_FLOAT stateMPC_Yd6[3];
stateMPC_FLOAT stateMPC_Ld6[3];
stateMPC_FLOAT stateMPC_yy6[2];
stateMPC_FLOAT stateMPC_bmy6[2];
int stateMPC_lbIdx6[4] = {0, 1, 2, 3};
stateMPC_FLOAT* stateMPC_llb6 = stateMPC_l + 48;
stateMPC_FLOAT* stateMPC_slb6 = stateMPC_s + 48;
stateMPC_FLOAT* stateMPC_llbbyslb6 = stateMPC_lbys + 48;
stateMPC_FLOAT stateMPC_rilb6[4];
stateMPC_FLOAT* stateMPC_dllbaff6 = stateMPC_dl_aff + 48;
stateMPC_FLOAT* stateMPC_dslbaff6 = stateMPC_ds_aff + 48;
stateMPC_FLOAT* stateMPC_dllbcc6 = stateMPC_dl_cc + 48;
stateMPC_FLOAT* stateMPC_dslbcc6 = stateMPC_ds_cc + 48;
stateMPC_FLOAT* stateMPC_ccrhsl6 = stateMPC_ccrhs + 48;
int stateMPC_ubIdx6[4] = {0, 1, 2, 3};
stateMPC_FLOAT* stateMPC_lub6 = stateMPC_l + 52;
stateMPC_FLOAT* stateMPC_sub6 = stateMPC_s + 52;
stateMPC_FLOAT* stateMPC_lubbysub6 = stateMPC_lbys + 52;
stateMPC_FLOAT stateMPC_riub6[4];
stateMPC_FLOAT* stateMPC_dlubaff6 = stateMPC_dl_aff + 52;
stateMPC_FLOAT* stateMPC_dsubaff6 = stateMPC_ds_aff + 52;
stateMPC_FLOAT* stateMPC_dlubcc6 = stateMPC_dl_cc + 52;
stateMPC_FLOAT* stateMPC_dsubcc6 = stateMPC_ds_cc + 52;
stateMPC_FLOAT* stateMPC_ccrhsub6 = stateMPC_ccrhs + 52;
stateMPC_FLOAT stateMPC_Phi6[4];
stateMPC_FLOAT stateMPC_W6[4];
stateMPC_FLOAT stateMPC_Ysd6[4];
stateMPC_FLOAT stateMPC_Lsd6[4];
stateMPC_FLOAT* stateMPC_z7 = stateMPC_z + 28;
stateMPC_FLOAT* stateMPC_dzaff7 = stateMPC_dz_aff + 28;
stateMPC_FLOAT* stateMPC_dzcc7 = stateMPC_dz_cc + 28;
stateMPC_FLOAT* stateMPC_rd7 = stateMPC_rd + 28;
stateMPC_FLOAT stateMPC_Lbyrd7[4];
stateMPC_FLOAT* stateMPC_grad_cost7 = stateMPC_grad_cost + 28;
stateMPC_FLOAT* stateMPC_grad_eq7 = stateMPC_grad_eq + 28;
stateMPC_FLOAT* stateMPC_grad_ineq7 = stateMPC_grad_ineq + 28;
stateMPC_FLOAT stateMPC_ctv7[4];
stateMPC_FLOAT* stateMPC_v7 = stateMPC_v + 14;
stateMPC_FLOAT stateMPC_re7[2];
stateMPC_FLOAT stateMPC_beta7[2];
stateMPC_FLOAT stateMPC_betacc7[2];
stateMPC_FLOAT* stateMPC_dvaff7 = stateMPC_dv_aff + 14;
stateMPC_FLOAT* stateMPC_dvcc7 = stateMPC_dv_cc + 14;
stateMPC_FLOAT stateMPC_V7[8];
stateMPC_FLOAT stateMPC_Yd7[3];
stateMPC_FLOAT stateMPC_Ld7[3];
stateMPC_FLOAT stateMPC_yy7[2];
stateMPC_FLOAT stateMPC_bmy7[2];
int stateMPC_lbIdx7[4] = {0, 1, 2, 3};
stateMPC_FLOAT* stateMPC_llb7 = stateMPC_l + 56;
stateMPC_FLOAT* stateMPC_slb7 = stateMPC_s + 56;
stateMPC_FLOAT* stateMPC_llbbyslb7 = stateMPC_lbys + 56;
stateMPC_FLOAT stateMPC_rilb7[4];
stateMPC_FLOAT* stateMPC_dllbaff7 = stateMPC_dl_aff + 56;
stateMPC_FLOAT* stateMPC_dslbaff7 = stateMPC_ds_aff + 56;
stateMPC_FLOAT* stateMPC_dllbcc7 = stateMPC_dl_cc + 56;
stateMPC_FLOAT* stateMPC_dslbcc7 = stateMPC_ds_cc + 56;
stateMPC_FLOAT* stateMPC_ccrhsl7 = stateMPC_ccrhs + 56;
int stateMPC_ubIdx7[4] = {0, 1, 2, 3};
stateMPC_FLOAT* stateMPC_lub7 = stateMPC_l + 60;
stateMPC_FLOAT* stateMPC_sub7 = stateMPC_s + 60;
stateMPC_FLOAT* stateMPC_lubbysub7 = stateMPC_lbys + 60;
stateMPC_FLOAT stateMPC_riub7[4];
stateMPC_FLOAT* stateMPC_dlubaff7 = stateMPC_dl_aff + 60;
stateMPC_FLOAT* stateMPC_dsubaff7 = stateMPC_ds_aff + 60;
stateMPC_FLOAT* stateMPC_dlubcc7 = stateMPC_dl_cc + 60;
stateMPC_FLOAT* stateMPC_dsubcc7 = stateMPC_ds_cc + 60;
stateMPC_FLOAT* stateMPC_ccrhsub7 = stateMPC_ccrhs + 60;
stateMPC_FLOAT stateMPC_Phi7[4];
stateMPC_FLOAT stateMPC_W7[4];
stateMPC_FLOAT stateMPC_Ysd7[4];
stateMPC_FLOAT stateMPC_Lsd7[4];
stateMPC_FLOAT* stateMPC_z8 = stateMPC_z + 32;
stateMPC_FLOAT* stateMPC_dzaff8 = stateMPC_dz_aff + 32;
stateMPC_FLOAT* stateMPC_dzcc8 = stateMPC_dz_cc + 32;
stateMPC_FLOAT* stateMPC_rd8 = stateMPC_rd + 32;
stateMPC_FLOAT stateMPC_Lbyrd8[4];
stateMPC_FLOAT* stateMPC_grad_cost8 = stateMPC_grad_cost + 32;
stateMPC_FLOAT* stateMPC_grad_eq8 = stateMPC_grad_eq + 32;
stateMPC_FLOAT* stateMPC_grad_ineq8 = stateMPC_grad_ineq + 32;
stateMPC_FLOAT stateMPC_ctv8[4];
stateMPC_FLOAT* stateMPC_v8 = stateMPC_v + 16;
stateMPC_FLOAT stateMPC_re8[2];
stateMPC_FLOAT stateMPC_beta8[2];
stateMPC_FLOAT stateMPC_betacc8[2];
stateMPC_FLOAT* stateMPC_dvaff8 = stateMPC_dv_aff + 16;
stateMPC_FLOAT* stateMPC_dvcc8 = stateMPC_dv_cc + 16;
stateMPC_FLOAT stateMPC_V8[8];
stateMPC_FLOAT stateMPC_Yd8[3];
stateMPC_FLOAT stateMPC_Ld8[3];
stateMPC_FLOAT stateMPC_yy8[2];
stateMPC_FLOAT stateMPC_bmy8[2];
int stateMPC_lbIdx8[4] = {0, 1, 2, 3};
stateMPC_FLOAT* stateMPC_llb8 = stateMPC_l + 64;
stateMPC_FLOAT* stateMPC_slb8 = stateMPC_s + 64;
stateMPC_FLOAT* stateMPC_llbbyslb8 = stateMPC_lbys + 64;
stateMPC_FLOAT stateMPC_rilb8[4];
stateMPC_FLOAT* stateMPC_dllbaff8 = stateMPC_dl_aff + 64;
stateMPC_FLOAT* stateMPC_dslbaff8 = stateMPC_ds_aff + 64;
stateMPC_FLOAT* stateMPC_dllbcc8 = stateMPC_dl_cc + 64;
stateMPC_FLOAT* stateMPC_dslbcc8 = stateMPC_ds_cc + 64;
stateMPC_FLOAT* stateMPC_ccrhsl8 = stateMPC_ccrhs + 64;
int stateMPC_ubIdx8[4] = {0, 1, 2, 3};
stateMPC_FLOAT* stateMPC_lub8 = stateMPC_l + 68;
stateMPC_FLOAT* stateMPC_sub8 = stateMPC_s + 68;
stateMPC_FLOAT* stateMPC_lubbysub8 = stateMPC_lbys + 68;
stateMPC_FLOAT stateMPC_riub8[4];
stateMPC_FLOAT* stateMPC_dlubaff8 = stateMPC_dl_aff + 68;
stateMPC_FLOAT* stateMPC_dsubaff8 = stateMPC_ds_aff + 68;
stateMPC_FLOAT* stateMPC_dlubcc8 = stateMPC_dl_cc + 68;
stateMPC_FLOAT* stateMPC_dsubcc8 = stateMPC_ds_cc + 68;
stateMPC_FLOAT* stateMPC_ccrhsub8 = stateMPC_ccrhs + 68;
stateMPC_FLOAT stateMPC_Phi8[4];
stateMPC_FLOAT stateMPC_W8[4];
stateMPC_FLOAT stateMPC_Ysd8[4];
stateMPC_FLOAT stateMPC_Lsd8[4];
stateMPC_FLOAT* stateMPC_z9 = stateMPC_z + 36;
stateMPC_FLOAT* stateMPC_dzaff9 = stateMPC_dz_aff + 36;
stateMPC_FLOAT* stateMPC_dzcc9 = stateMPC_dz_cc + 36;
stateMPC_FLOAT* stateMPC_rd9 = stateMPC_rd + 36;
stateMPC_FLOAT stateMPC_Lbyrd9[2];
stateMPC_FLOAT* stateMPC_grad_cost9 = stateMPC_grad_cost + 36;
stateMPC_FLOAT* stateMPC_grad_eq9 = stateMPC_grad_eq + 36;
stateMPC_FLOAT* stateMPC_grad_ineq9 = stateMPC_grad_ineq + 36;
stateMPC_FLOAT stateMPC_ctv9[2];
stateMPC_FLOAT* stateMPC_v9 = stateMPC_v + 18;
stateMPC_FLOAT stateMPC_re9[2];
stateMPC_FLOAT stateMPC_beta9[2];
stateMPC_FLOAT stateMPC_betacc9[2];
stateMPC_FLOAT* stateMPC_dvaff9 = stateMPC_dv_aff + 18;
stateMPC_FLOAT* stateMPC_dvcc9 = stateMPC_dv_cc + 18;
stateMPC_FLOAT stateMPC_V9[4];
stateMPC_FLOAT stateMPC_Yd9[3];
stateMPC_FLOAT stateMPC_Ld9[3];
stateMPC_FLOAT stateMPC_yy9[2];
stateMPC_FLOAT stateMPC_bmy9[2];
int stateMPC_lbIdx9[2] = {0, 1};
stateMPC_FLOAT* stateMPC_llb9 = stateMPC_l + 72;
stateMPC_FLOAT* stateMPC_slb9 = stateMPC_s + 72;
stateMPC_FLOAT* stateMPC_llbbyslb9 = stateMPC_lbys + 72;
stateMPC_FLOAT stateMPC_rilb9[2];
stateMPC_FLOAT* stateMPC_dllbaff9 = stateMPC_dl_aff + 72;
stateMPC_FLOAT* stateMPC_dslbaff9 = stateMPC_ds_aff + 72;
stateMPC_FLOAT* stateMPC_dllbcc9 = stateMPC_dl_cc + 72;
stateMPC_FLOAT* stateMPC_dslbcc9 = stateMPC_ds_cc + 72;
stateMPC_FLOAT* stateMPC_ccrhsl9 = stateMPC_ccrhs + 72;
int stateMPC_ubIdx9[2] = {0, 1};
stateMPC_FLOAT* stateMPC_lub9 = stateMPC_l + 74;
stateMPC_FLOAT* stateMPC_sub9 = stateMPC_s + 74;
stateMPC_FLOAT* stateMPC_lubbysub9 = stateMPC_lbys + 74;
stateMPC_FLOAT stateMPC_riub9[2];
stateMPC_FLOAT* stateMPC_dlubaff9 = stateMPC_dl_aff + 74;
stateMPC_FLOAT* stateMPC_dsubaff9 = stateMPC_ds_aff + 74;
stateMPC_FLOAT* stateMPC_dlubcc9 = stateMPC_dl_cc + 74;
stateMPC_FLOAT* stateMPC_dsubcc9 = stateMPC_ds_cc + 74;
stateMPC_FLOAT* stateMPC_ccrhsub9 = stateMPC_ccrhs + 74;
stateMPC_FLOAT stateMPC_Phi9[2];
stateMPC_FLOAT stateMPC_D9[2] = {-1.0000000000000000E+000, 
-1.0000000000000000E+000};
stateMPC_FLOAT stateMPC_W9[2];
stateMPC_FLOAT stateMPC_Ysd9[4];
stateMPC_FLOAT stateMPC_Lsd9[4];
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
stateMPC_LA_INITIALIZEVECTOR_38(stateMPC_z, 0);
stateMPC_LA_INITIALIZEVECTOR_20(stateMPC_v, 1);
stateMPC_LA_INITIALIZEVECTOR_76(stateMPC_l, 1);
stateMPC_LA_INITIALIZEVECTOR_76(stateMPC_s, 1);
info->mu = 0;
stateMPC_LA_DOTACC_76(stateMPC_l, stateMPC_s, &info->mu);
info->mu /= 76;
while( 1 ){
info->pobj = 0;
stateMPC_LA_DIAG_QUADFCN_4(params->H1, params->f1, stateMPC_z0, stateMPC_grad_cost0, &info->pobj);
stateMPC_LA_DIAG_QUADFCN_4(params->H2, params->f2, stateMPC_z1, stateMPC_grad_cost1, &info->pobj);
stateMPC_LA_DIAG_QUADFCN_4(params->H3, params->f3, stateMPC_z2, stateMPC_grad_cost2, &info->pobj);
stateMPC_LA_DIAG_QUADFCN_4(params->H4, params->f4, stateMPC_z3, stateMPC_grad_cost3, &info->pobj);
stateMPC_LA_DIAG_QUADFCN_4(params->H5, params->f5, stateMPC_z4, stateMPC_grad_cost4, &info->pobj);
stateMPC_LA_DIAG_QUADFCN_4(params->H6, params->f6, stateMPC_z5, stateMPC_grad_cost5, &info->pobj);
stateMPC_LA_DIAG_QUADFCN_4(params->H7, params->f7, stateMPC_z6, stateMPC_grad_cost6, &info->pobj);
stateMPC_LA_DIAG_QUADFCN_4(params->H8, params->f8, stateMPC_z7, stateMPC_grad_cost7, &info->pobj);
stateMPC_LA_DIAG_QUADFCN_4(params->H9, params->f9, stateMPC_z8, stateMPC_grad_cost8, &info->pobj);
stateMPC_LA_DIAG_QUADFCN_2(params->H10, params->f10, stateMPC_z9, stateMPC_grad_cost9, &info->pobj);
info->res_eq = 0;
info->dgap = 0;
stateMPC_LA_DIAGZERO_MVMSUB6_2(stateMPC_D0, stateMPC_z0, params->e1, stateMPC_v0, stateMPC_re0, &info->dgap, &info->res_eq);
stateMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_4_4(params->C1, stateMPC_z0, stateMPC_D1, stateMPC_z1, params->e2, stateMPC_v1, stateMPC_re1, &info->dgap, &info->res_eq);
stateMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_4_4(params->C2, stateMPC_z1, stateMPC_D1, stateMPC_z2, params->e3, stateMPC_v2, stateMPC_re2, &info->dgap, &info->res_eq);
stateMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_4_4(params->C3, stateMPC_z2, stateMPC_D1, stateMPC_z3, params->e4, stateMPC_v3, stateMPC_re3, &info->dgap, &info->res_eq);
stateMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_4_4(params->C4, stateMPC_z3, stateMPC_D1, stateMPC_z4, params->e5, stateMPC_v4, stateMPC_re4, &info->dgap, &info->res_eq);
stateMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_4_4(params->C5, stateMPC_z4, stateMPC_D1, stateMPC_z5, params->e6, stateMPC_v5, stateMPC_re5, &info->dgap, &info->res_eq);
stateMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_4_4(params->C6, stateMPC_z5, stateMPC_D1, stateMPC_z6, params->e7, stateMPC_v6, stateMPC_re6, &info->dgap, &info->res_eq);
stateMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_4_4(params->C7, stateMPC_z6, stateMPC_D1, stateMPC_z7, params->e8, stateMPC_v7, stateMPC_re7, &info->dgap, &info->res_eq);
stateMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_4_4(params->C8, stateMPC_z7, stateMPC_D1, stateMPC_z8, params->e9, stateMPC_v8, stateMPC_re8, &info->dgap, &info->res_eq);
stateMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_4_2(params->C9, stateMPC_z8, stateMPC_D9, stateMPC_z9, params->e10, stateMPC_v9, stateMPC_re9, &info->dgap, &info->res_eq);
stateMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C1, stateMPC_v1, stateMPC_D0, stateMPC_v0, stateMPC_grad_eq0);
stateMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C2, stateMPC_v2, stateMPC_D1, stateMPC_v1, stateMPC_grad_eq1);
stateMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C3, stateMPC_v3, stateMPC_D1, stateMPC_v2, stateMPC_grad_eq2);
stateMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C4, stateMPC_v4, stateMPC_D1, stateMPC_v3, stateMPC_grad_eq3);
stateMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C5, stateMPC_v5, stateMPC_D1, stateMPC_v4, stateMPC_grad_eq4);
stateMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C6, stateMPC_v6, stateMPC_D1, stateMPC_v5, stateMPC_grad_eq5);
stateMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C7, stateMPC_v7, stateMPC_D1, stateMPC_v6, stateMPC_grad_eq6);
stateMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C8, stateMPC_v8, stateMPC_D1, stateMPC_v7, stateMPC_grad_eq7);
stateMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C9, stateMPC_v9, stateMPC_D1, stateMPC_v8, stateMPC_grad_eq8);
stateMPC_LA_DIAGZERO_MTVM_2_2(stateMPC_D9, stateMPC_v9, stateMPC_grad_eq9);
info->res_ineq = 0;
stateMPC_LA_VSUBADD3_4(params->lb1, stateMPC_z0, stateMPC_lbIdx0, stateMPC_llb0, stateMPC_slb0, stateMPC_rilb0, &info->dgap, &info->res_ineq);
stateMPC_LA_VSUBADD2_4(stateMPC_z0, stateMPC_ubIdx0, params->ub1, stateMPC_lub0, stateMPC_sub0, stateMPC_riub0, &info->dgap, &info->res_ineq);
stateMPC_LA_VSUBADD3_4(params->lb2, stateMPC_z1, stateMPC_lbIdx1, stateMPC_llb1, stateMPC_slb1, stateMPC_rilb1, &info->dgap, &info->res_ineq);
stateMPC_LA_VSUBADD2_4(stateMPC_z1, stateMPC_ubIdx1, params->ub2, stateMPC_lub1, stateMPC_sub1, stateMPC_riub1, &info->dgap, &info->res_ineq);
stateMPC_LA_VSUBADD3_4(params->lb3, stateMPC_z2, stateMPC_lbIdx2, stateMPC_llb2, stateMPC_slb2, stateMPC_rilb2, &info->dgap, &info->res_ineq);
stateMPC_LA_VSUBADD2_4(stateMPC_z2, stateMPC_ubIdx2, params->ub3, stateMPC_lub2, stateMPC_sub2, stateMPC_riub2, &info->dgap, &info->res_ineq);
stateMPC_LA_VSUBADD3_4(params->lb4, stateMPC_z3, stateMPC_lbIdx3, stateMPC_llb3, stateMPC_slb3, stateMPC_rilb3, &info->dgap, &info->res_ineq);
stateMPC_LA_VSUBADD2_4(stateMPC_z3, stateMPC_ubIdx3, params->ub4, stateMPC_lub3, stateMPC_sub3, stateMPC_riub3, &info->dgap, &info->res_ineq);
stateMPC_LA_VSUBADD3_4(params->lb5, stateMPC_z4, stateMPC_lbIdx4, stateMPC_llb4, stateMPC_slb4, stateMPC_rilb4, &info->dgap, &info->res_ineq);
stateMPC_LA_VSUBADD2_4(stateMPC_z4, stateMPC_ubIdx4, params->ub5, stateMPC_lub4, stateMPC_sub4, stateMPC_riub4, &info->dgap, &info->res_ineq);
stateMPC_LA_VSUBADD3_4(params->lb6, stateMPC_z5, stateMPC_lbIdx5, stateMPC_llb5, stateMPC_slb5, stateMPC_rilb5, &info->dgap, &info->res_ineq);
stateMPC_LA_VSUBADD2_4(stateMPC_z5, stateMPC_ubIdx5, params->ub6, stateMPC_lub5, stateMPC_sub5, stateMPC_riub5, &info->dgap, &info->res_ineq);
stateMPC_LA_VSUBADD3_4(params->lb7, stateMPC_z6, stateMPC_lbIdx6, stateMPC_llb6, stateMPC_slb6, stateMPC_rilb6, &info->dgap, &info->res_ineq);
stateMPC_LA_VSUBADD2_4(stateMPC_z6, stateMPC_ubIdx6, params->ub7, stateMPC_lub6, stateMPC_sub6, stateMPC_riub6, &info->dgap, &info->res_ineq);
stateMPC_LA_VSUBADD3_4(params->lb8, stateMPC_z7, stateMPC_lbIdx7, stateMPC_llb7, stateMPC_slb7, stateMPC_rilb7, &info->dgap, &info->res_ineq);
stateMPC_LA_VSUBADD2_4(stateMPC_z7, stateMPC_ubIdx7, params->ub8, stateMPC_lub7, stateMPC_sub7, stateMPC_riub7, &info->dgap, &info->res_ineq);
stateMPC_LA_VSUBADD3_4(params->lb9, stateMPC_z8, stateMPC_lbIdx8, stateMPC_llb8, stateMPC_slb8, stateMPC_rilb8, &info->dgap, &info->res_ineq);
stateMPC_LA_VSUBADD2_4(stateMPC_z8, stateMPC_ubIdx8, params->ub9, stateMPC_lub8, stateMPC_sub8, stateMPC_riub8, &info->dgap, &info->res_ineq);
stateMPC_LA_VSUBADD3_2(params->lb10, stateMPC_z9, stateMPC_lbIdx9, stateMPC_llb9, stateMPC_slb9, stateMPC_rilb9, &info->dgap, &info->res_ineq);
stateMPC_LA_VSUBADD2_2(stateMPC_z9, stateMPC_ubIdx9, params->ub10, stateMPC_lub9, stateMPC_sub9, stateMPC_riub9, &info->dgap, &info->res_ineq);
stateMPC_LA_INEQ_B_GRAD_4_4_4(stateMPC_lub0, stateMPC_sub0, stateMPC_riub0, stateMPC_llb0, stateMPC_slb0, stateMPC_rilb0, stateMPC_lbIdx0, stateMPC_ubIdx0, stateMPC_grad_ineq0, stateMPC_lubbysub0, stateMPC_llbbyslb0);
stateMPC_LA_INEQ_B_GRAD_4_4_4(stateMPC_lub1, stateMPC_sub1, stateMPC_riub1, stateMPC_llb1, stateMPC_slb1, stateMPC_rilb1, stateMPC_lbIdx1, stateMPC_ubIdx1, stateMPC_grad_ineq1, stateMPC_lubbysub1, stateMPC_llbbyslb1);
stateMPC_LA_INEQ_B_GRAD_4_4_4(stateMPC_lub2, stateMPC_sub2, stateMPC_riub2, stateMPC_llb2, stateMPC_slb2, stateMPC_rilb2, stateMPC_lbIdx2, stateMPC_ubIdx2, stateMPC_grad_ineq2, stateMPC_lubbysub2, stateMPC_llbbyslb2);
stateMPC_LA_INEQ_B_GRAD_4_4_4(stateMPC_lub3, stateMPC_sub3, stateMPC_riub3, stateMPC_llb3, stateMPC_slb3, stateMPC_rilb3, stateMPC_lbIdx3, stateMPC_ubIdx3, stateMPC_grad_ineq3, stateMPC_lubbysub3, stateMPC_llbbyslb3);
stateMPC_LA_INEQ_B_GRAD_4_4_4(stateMPC_lub4, stateMPC_sub4, stateMPC_riub4, stateMPC_llb4, stateMPC_slb4, stateMPC_rilb4, stateMPC_lbIdx4, stateMPC_ubIdx4, stateMPC_grad_ineq4, stateMPC_lubbysub4, stateMPC_llbbyslb4);
stateMPC_LA_INEQ_B_GRAD_4_4_4(stateMPC_lub5, stateMPC_sub5, stateMPC_riub5, stateMPC_llb5, stateMPC_slb5, stateMPC_rilb5, stateMPC_lbIdx5, stateMPC_ubIdx5, stateMPC_grad_ineq5, stateMPC_lubbysub5, stateMPC_llbbyslb5);
stateMPC_LA_INEQ_B_GRAD_4_4_4(stateMPC_lub6, stateMPC_sub6, stateMPC_riub6, stateMPC_llb6, stateMPC_slb6, stateMPC_rilb6, stateMPC_lbIdx6, stateMPC_ubIdx6, stateMPC_grad_ineq6, stateMPC_lubbysub6, stateMPC_llbbyslb6);
stateMPC_LA_INEQ_B_GRAD_4_4_4(stateMPC_lub7, stateMPC_sub7, stateMPC_riub7, stateMPC_llb7, stateMPC_slb7, stateMPC_rilb7, stateMPC_lbIdx7, stateMPC_ubIdx7, stateMPC_grad_ineq7, stateMPC_lubbysub7, stateMPC_llbbyslb7);
stateMPC_LA_INEQ_B_GRAD_4_4_4(stateMPC_lub8, stateMPC_sub8, stateMPC_riub8, stateMPC_llb8, stateMPC_slb8, stateMPC_rilb8, stateMPC_lbIdx8, stateMPC_ubIdx8, stateMPC_grad_ineq8, stateMPC_lubbysub8, stateMPC_llbbyslb8);
stateMPC_LA_INEQ_B_GRAD_2_2_2(stateMPC_lub9, stateMPC_sub9, stateMPC_riub9, stateMPC_llb9, stateMPC_slb9, stateMPC_rilb9, stateMPC_lbIdx9, stateMPC_ubIdx9, stateMPC_grad_ineq9, stateMPC_lubbysub9, stateMPC_llbbyslb9);
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
stateMPC_LA_VVADD3_38(stateMPC_grad_cost, stateMPC_grad_eq, stateMPC_grad_ineq, stateMPC_rd);
stateMPC_LA_DIAG_CHOL_ONELOOP_LBUB_4_4_4(params->H1, stateMPC_llbbyslb0, stateMPC_lbIdx0, stateMPC_lubbysub0, stateMPC_ubIdx0, stateMPC_Phi0);
stateMPC_LA_DIAG_MATRIXFORWARDSUB_2_4(stateMPC_Phi0, params->C1, stateMPC_V0);
stateMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_4(stateMPC_Phi0, stateMPC_D0, stateMPC_W0);
stateMPC_LA_DENSE_DIAGZERO_MMTM_2_4_2(stateMPC_W0, stateMPC_V0, stateMPC_Ysd1);
stateMPC_LA_DIAG_FORWARDSUB_4(stateMPC_Phi0, stateMPC_rd0, stateMPC_Lbyrd0);
stateMPC_LA_DIAG_CHOL_ONELOOP_LBUB_4_4_4(params->H2, stateMPC_llbbyslb1, stateMPC_lbIdx1, stateMPC_lubbysub1, stateMPC_ubIdx1, stateMPC_Phi1);
stateMPC_LA_DIAG_MATRIXFORWARDSUB_2_4(stateMPC_Phi1, params->C2, stateMPC_V1);
stateMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_4(stateMPC_Phi1, stateMPC_D1, stateMPC_W1);
stateMPC_LA_DENSE_DIAGZERO_MMTM_2_4_2(stateMPC_W1, stateMPC_V1, stateMPC_Ysd2);
stateMPC_LA_DIAG_FORWARDSUB_4(stateMPC_Phi1, stateMPC_rd1, stateMPC_Lbyrd1);
stateMPC_LA_DIAG_CHOL_ONELOOP_LBUB_4_4_4(params->H3, stateMPC_llbbyslb2, stateMPC_lbIdx2, stateMPC_lubbysub2, stateMPC_ubIdx2, stateMPC_Phi2);
stateMPC_LA_DIAG_MATRIXFORWARDSUB_2_4(stateMPC_Phi2, params->C3, stateMPC_V2);
stateMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_4(stateMPC_Phi2, stateMPC_D1, stateMPC_W2);
stateMPC_LA_DENSE_DIAGZERO_MMTM_2_4_2(stateMPC_W2, stateMPC_V2, stateMPC_Ysd3);
stateMPC_LA_DIAG_FORWARDSUB_4(stateMPC_Phi2, stateMPC_rd2, stateMPC_Lbyrd2);
stateMPC_LA_DIAG_CHOL_ONELOOP_LBUB_4_4_4(params->H4, stateMPC_llbbyslb3, stateMPC_lbIdx3, stateMPC_lubbysub3, stateMPC_ubIdx3, stateMPC_Phi3);
stateMPC_LA_DIAG_MATRIXFORWARDSUB_2_4(stateMPC_Phi3, params->C4, stateMPC_V3);
stateMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_4(stateMPC_Phi3, stateMPC_D1, stateMPC_W3);
stateMPC_LA_DENSE_DIAGZERO_MMTM_2_4_2(stateMPC_W3, stateMPC_V3, stateMPC_Ysd4);
stateMPC_LA_DIAG_FORWARDSUB_4(stateMPC_Phi3, stateMPC_rd3, stateMPC_Lbyrd3);
stateMPC_LA_DIAG_CHOL_ONELOOP_LBUB_4_4_4(params->H5, stateMPC_llbbyslb4, stateMPC_lbIdx4, stateMPC_lubbysub4, stateMPC_ubIdx4, stateMPC_Phi4);
stateMPC_LA_DIAG_MATRIXFORWARDSUB_2_4(stateMPC_Phi4, params->C5, stateMPC_V4);
stateMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_4(stateMPC_Phi4, stateMPC_D1, stateMPC_W4);
stateMPC_LA_DENSE_DIAGZERO_MMTM_2_4_2(stateMPC_W4, stateMPC_V4, stateMPC_Ysd5);
stateMPC_LA_DIAG_FORWARDSUB_4(stateMPC_Phi4, stateMPC_rd4, stateMPC_Lbyrd4);
stateMPC_LA_DIAG_CHOL_ONELOOP_LBUB_4_4_4(params->H6, stateMPC_llbbyslb5, stateMPC_lbIdx5, stateMPC_lubbysub5, stateMPC_ubIdx5, stateMPC_Phi5);
stateMPC_LA_DIAG_MATRIXFORWARDSUB_2_4(stateMPC_Phi5, params->C6, stateMPC_V5);
stateMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_4(stateMPC_Phi5, stateMPC_D1, stateMPC_W5);
stateMPC_LA_DENSE_DIAGZERO_MMTM_2_4_2(stateMPC_W5, stateMPC_V5, stateMPC_Ysd6);
stateMPC_LA_DIAG_FORWARDSUB_4(stateMPC_Phi5, stateMPC_rd5, stateMPC_Lbyrd5);
stateMPC_LA_DIAG_CHOL_ONELOOP_LBUB_4_4_4(params->H7, stateMPC_llbbyslb6, stateMPC_lbIdx6, stateMPC_lubbysub6, stateMPC_ubIdx6, stateMPC_Phi6);
stateMPC_LA_DIAG_MATRIXFORWARDSUB_2_4(stateMPC_Phi6, params->C7, stateMPC_V6);
stateMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_4(stateMPC_Phi6, stateMPC_D1, stateMPC_W6);
stateMPC_LA_DENSE_DIAGZERO_MMTM_2_4_2(stateMPC_W6, stateMPC_V6, stateMPC_Ysd7);
stateMPC_LA_DIAG_FORWARDSUB_4(stateMPC_Phi6, stateMPC_rd6, stateMPC_Lbyrd6);
stateMPC_LA_DIAG_CHOL_ONELOOP_LBUB_4_4_4(params->H8, stateMPC_llbbyslb7, stateMPC_lbIdx7, stateMPC_lubbysub7, stateMPC_ubIdx7, stateMPC_Phi7);
stateMPC_LA_DIAG_MATRIXFORWARDSUB_2_4(stateMPC_Phi7, params->C8, stateMPC_V7);
stateMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_4(stateMPC_Phi7, stateMPC_D1, stateMPC_W7);
stateMPC_LA_DENSE_DIAGZERO_MMTM_2_4_2(stateMPC_W7, stateMPC_V7, stateMPC_Ysd8);
stateMPC_LA_DIAG_FORWARDSUB_4(stateMPC_Phi7, stateMPC_rd7, stateMPC_Lbyrd7);
stateMPC_LA_DIAG_CHOL_ONELOOP_LBUB_4_4_4(params->H9, stateMPC_llbbyslb8, stateMPC_lbIdx8, stateMPC_lubbysub8, stateMPC_ubIdx8, stateMPC_Phi8);
stateMPC_LA_DIAG_MATRIXFORWARDSUB_2_4(stateMPC_Phi8, params->C9, stateMPC_V8);
stateMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_4(stateMPC_Phi8, stateMPC_D1, stateMPC_W8);
stateMPC_LA_DENSE_DIAGZERO_MMTM_2_4_2(stateMPC_W8, stateMPC_V8, stateMPC_Ysd9);
stateMPC_LA_DIAG_FORWARDSUB_4(stateMPC_Phi8, stateMPC_rd8, stateMPC_Lbyrd8);
stateMPC_LA_DIAG_CHOL_ONELOOP_LBUB_2_2_2(params->H10, stateMPC_llbbyslb9, stateMPC_lbIdx9, stateMPC_lubbysub9, stateMPC_ubIdx9, stateMPC_Phi9);
stateMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_2(stateMPC_Phi9, stateMPC_D9, stateMPC_W9);
stateMPC_LA_DIAG_FORWARDSUB_2(stateMPC_Phi9, stateMPC_rd9, stateMPC_Lbyrd9);
stateMPC_LA_DIAGZERO_MMT_2(stateMPC_W0, stateMPC_Yd0);
stateMPC_LA_DIAGZERO_MVMSUB7_2(stateMPC_W0, stateMPC_Lbyrd0, stateMPC_re0, stateMPC_beta0);
stateMPC_LA_DENSE_DIAGZERO_MMT2_2_4_4(stateMPC_V0, stateMPC_W1, stateMPC_Yd1);
stateMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_4_4(stateMPC_V0, stateMPC_Lbyrd0, stateMPC_W1, stateMPC_Lbyrd1, stateMPC_re1, stateMPC_beta1);
stateMPC_LA_DENSE_DIAGZERO_MMT2_2_4_4(stateMPC_V1, stateMPC_W2, stateMPC_Yd2);
stateMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_4_4(stateMPC_V1, stateMPC_Lbyrd1, stateMPC_W2, stateMPC_Lbyrd2, stateMPC_re2, stateMPC_beta2);
stateMPC_LA_DENSE_DIAGZERO_MMT2_2_4_4(stateMPC_V2, stateMPC_W3, stateMPC_Yd3);
stateMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_4_4(stateMPC_V2, stateMPC_Lbyrd2, stateMPC_W3, stateMPC_Lbyrd3, stateMPC_re3, stateMPC_beta3);
stateMPC_LA_DENSE_DIAGZERO_MMT2_2_4_4(stateMPC_V3, stateMPC_W4, stateMPC_Yd4);
stateMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_4_4(stateMPC_V3, stateMPC_Lbyrd3, stateMPC_W4, stateMPC_Lbyrd4, stateMPC_re4, stateMPC_beta4);
stateMPC_LA_DENSE_DIAGZERO_MMT2_2_4_4(stateMPC_V4, stateMPC_W5, stateMPC_Yd5);
stateMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_4_4(stateMPC_V4, stateMPC_Lbyrd4, stateMPC_W5, stateMPC_Lbyrd5, stateMPC_re5, stateMPC_beta5);
stateMPC_LA_DENSE_DIAGZERO_MMT2_2_4_4(stateMPC_V5, stateMPC_W6, stateMPC_Yd6);
stateMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_4_4(stateMPC_V5, stateMPC_Lbyrd5, stateMPC_W6, stateMPC_Lbyrd6, stateMPC_re6, stateMPC_beta6);
stateMPC_LA_DENSE_DIAGZERO_MMT2_2_4_4(stateMPC_V6, stateMPC_W7, stateMPC_Yd7);
stateMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_4_4(stateMPC_V6, stateMPC_Lbyrd6, stateMPC_W7, stateMPC_Lbyrd7, stateMPC_re7, stateMPC_beta7);
stateMPC_LA_DENSE_DIAGZERO_MMT2_2_4_4(stateMPC_V7, stateMPC_W8, stateMPC_Yd8);
stateMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_4_4(stateMPC_V7, stateMPC_Lbyrd7, stateMPC_W8, stateMPC_Lbyrd8, stateMPC_re8, stateMPC_beta8);
stateMPC_LA_DENSE_DIAGZERO_MMT2_2_4_2(stateMPC_V8, stateMPC_W9, stateMPC_Yd9);
stateMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_4_2(stateMPC_V8, stateMPC_Lbyrd8, stateMPC_W9, stateMPC_Lbyrd9, stateMPC_re9, stateMPC_beta9);
stateMPC_LA_DENSE_CHOL_2(stateMPC_Yd0, stateMPC_Ld0);
stateMPC_LA_DENSE_FORWARDSUB_2(stateMPC_Ld0, stateMPC_beta0, stateMPC_yy0);
stateMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(stateMPC_Ld0, stateMPC_Ysd1, stateMPC_Lsd1);
stateMPC_LA_DENSE_MMTSUB_2_2(stateMPC_Lsd1, stateMPC_Yd1);
stateMPC_LA_DENSE_CHOL_2(stateMPC_Yd1, stateMPC_Ld1);
stateMPC_LA_DENSE_MVMSUB1_2_2(stateMPC_Lsd1, stateMPC_yy0, stateMPC_beta1, stateMPC_bmy1);
stateMPC_LA_DENSE_FORWARDSUB_2(stateMPC_Ld1, stateMPC_bmy1, stateMPC_yy1);
stateMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(stateMPC_Ld1, stateMPC_Ysd2, stateMPC_Lsd2);
stateMPC_LA_DENSE_MMTSUB_2_2(stateMPC_Lsd2, stateMPC_Yd2);
stateMPC_LA_DENSE_CHOL_2(stateMPC_Yd2, stateMPC_Ld2);
stateMPC_LA_DENSE_MVMSUB1_2_2(stateMPC_Lsd2, stateMPC_yy1, stateMPC_beta2, stateMPC_bmy2);
stateMPC_LA_DENSE_FORWARDSUB_2(stateMPC_Ld2, stateMPC_bmy2, stateMPC_yy2);
stateMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(stateMPC_Ld2, stateMPC_Ysd3, stateMPC_Lsd3);
stateMPC_LA_DENSE_MMTSUB_2_2(stateMPC_Lsd3, stateMPC_Yd3);
stateMPC_LA_DENSE_CHOL_2(stateMPC_Yd3, stateMPC_Ld3);
stateMPC_LA_DENSE_MVMSUB1_2_2(stateMPC_Lsd3, stateMPC_yy2, stateMPC_beta3, stateMPC_bmy3);
stateMPC_LA_DENSE_FORWARDSUB_2(stateMPC_Ld3, stateMPC_bmy3, stateMPC_yy3);
stateMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(stateMPC_Ld3, stateMPC_Ysd4, stateMPC_Lsd4);
stateMPC_LA_DENSE_MMTSUB_2_2(stateMPC_Lsd4, stateMPC_Yd4);
stateMPC_LA_DENSE_CHOL_2(stateMPC_Yd4, stateMPC_Ld4);
stateMPC_LA_DENSE_MVMSUB1_2_2(stateMPC_Lsd4, stateMPC_yy3, stateMPC_beta4, stateMPC_bmy4);
stateMPC_LA_DENSE_FORWARDSUB_2(stateMPC_Ld4, stateMPC_bmy4, stateMPC_yy4);
stateMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(stateMPC_Ld4, stateMPC_Ysd5, stateMPC_Lsd5);
stateMPC_LA_DENSE_MMTSUB_2_2(stateMPC_Lsd5, stateMPC_Yd5);
stateMPC_LA_DENSE_CHOL_2(stateMPC_Yd5, stateMPC_Ld5);
stateMPC_LA_DENSE_MVMSUB1_2_2(stateMPC_Lsd5, stateMPC_yy4, stateMPC_beta5, stateMPC_bmy5);
stateMPC_LA_DENSE_FORWARDSUB_2(stateMPC_Ld5, stateMPC_bmy5, stateMPC_yy5);
stateMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(stateMPC_Ld5, stateMPC_Ysd6, stateMPC_Lsd6);
stateMPC_LA_DENSE_MMTSUB_2_2(stateMPC_Lsd6, stateMPC_Yd6);
stateMPC_LA_DENSE_CHOL_2(stateMPC_Yd6, stateMPC_Ld6);
stateMPC_LA_DENSE_MVMSUB1_2_2(stateMPC_Lsd6, stateMPC_yy5, stateMPC_beta6, stateMPC_bmy6);
stateMPC_LA_DENSE_FORWARDSUB_2(stateMPC_Ld6, stateMPC_bmy6, stateMPC_yy6);
stateMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(stateMPC_Ld6, stateMPC_Ysd7, stateMPC_Lsd7);
stateMPC_LA_DENSE_MMTSUB_2_2(stateMPC_Lsd7, stateMPC_Yd7);
stateMPC_LA_DENSE_CHOL_2(stateMPC_Yd7, stateMPC_Ld7);
stateMPC_LA_DENSE_MVMSUB1_2_2(stateMPC_Lsd7, stateMPC_yy6, stateMPC_beta7, stateMPC_bmy7);
stateMPC_LA_DENSE_FORWARDSUB_2(stateMPC_Ld7, stateMPC_bmy7, stateMPC_yy7);
stateMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(stateMPC_Ld7, stateMPC_Ysd8, stateMPC_Lsd8);
stateMPC_LA_DENSE_MMTSUB_2_2(stateMPC_Lsd8, stateMPC_Yd8);
stateMPC_LA_DENSE_CHOL_2(stateMPC_Yd8, stateMPC_Ld8);
stateMPC_LA_DENSE_MVMSUB1_2_2(stateMPC_Lsd8, stateMPC_yy7, stateMPC_beta8, stateMPC_bmy8);
stateMPC_LA_DENSE_FORWARDSUB_2(stateMPC_Ld8, stateMPC_bmy8, stateMPC_yy8);
stateMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(stateMPC_Ld8, stateMPC_Ysd9, stateMPC_Lsd9);
stateMPC_LA_DENSE_MMTSUB_2_2(stateMPC_Lsd9, stateMPC_Yd9);
stateMPC_LA_DENSE_CHOL_2(stateMPC_Yd9, stateMPC_Ld9);
stateMPC_LA_DENSE_MVMSUB1_2_2(stateMPC_Lsd9, stateMPC_yy8, stateMPC_beta9, stateMPC_bmy9);
stateMPC_LA_DENSE_FORWARDSUB_2(stateMPC_Ld9, stateMPC_bmy9, stateMPC_yy9);
stateMPC_LA_DENSE_BACKWARDSUB_2(stateMPC_Ld9, stateMPC_yy9, stateMPC_dvaff9);
stateMPC_LA_DENSE_MTVMSUB_2_2(stateMPC_Lsd9, stateMPC_dvaff9, stateMPC_yy8, stateMPC_bmy8);
stateMPC_LA_DENSE_BACKWARDSUB_2(stateMPC_Ld8, stateMPC_bmy8, stateMPC_dvaff8);
stateMPC_LA_DENSE_MTVMSUB_2_2(stateMPC_Lsd8, stateMPC_dvaff8, stateMPC_yy7, stateMPC_bmy7);
stateMPC_LA_DENSE_BACKWARDSUB_2(stateMPC_Ld7, stateMPC_bmy7, stateMPC_dvaff7);
stateMPC_LA_DENSE_MTVMSUB_2_2(stateMPC_Lsd7, stateMPC_dvaff7, stateMPC_yy6, stateMPC_bmy6);
stateMPC_LA_DENSE_BACKWARDSUB_2(stateMPC_Ld6, stateMPC_bmy6, stateMPC_dvaff6);
stateMPC_LA_DENSE_MTVMSUB_2_2(stateMPC_Lsd6, stateMPC_dvaff6, stateMPC_yy5, stateMPC_bmy5);
stateMPC_LA_DENSE_BACKWARDSUB_2(stateMPC_Ld5, stateMPC_bmy5, stateMPC_dvaff5);
stateMPC_LA_DENSE_MTVMSUB_2_2(stateMPC_Lsd5, stateMPC_dvaff5, stateMPC_yy4, stateMPC_bmy4);
stateMPC_LA_DENSE_BACKWARDSUB_2(stateMPC_Ld4, stateMPC_bmy4, stateMPC_dvaff4);
stateMPC_LA_DENSE_MTVMSUB_2_2(stateMPC_Lsd4, stateMPC_dvaff4, stateMPC_yy3, stateMPC_bmy3);
stateMPC_LA_DENSE_BACKWARDSUB_2(stateMPC_Ld3, stateMPC_bmy3, stateMPC_dvaff3);
stateMPC_LA_DENSE_MTVMSUB_2_2(stateMPC_Lsd3, stateMPC_dvaff3, stateMPC_yy2, stateMPC_bmy2);
stateMPC_LA_DENSE_BACKWARDSUB_2(stateMPC_Ld2, stateMPC_bmy2, stateMPC_dvaff2);
stateMPC_LA_DENSE_MTVMSUB_2_2(stateMPC_Lsd2, stateMPC_dvaff2, stateMPC_yy1, stateMPC_bmy1);
stateMPC_LA_DENSE_BACKWARDSUB_2(stateMPC_Ld1, stateMPC_bmy1, stateMPC_dvaff1);
stateMPC_LA_DENSE_MTVMSUB_2_2(stateMPC_Lsd1, stateMPC_dvaff1, stateMPC_yy0, stateMPC_bmy0);
stateMPC_LA_DENSE_BACKWARDSUB_2(stateMPC_Ld0, stateMPC_bmy0, stateMPC_dvaff0);
stateMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C1, stateMPC_dvaff1, stateMPC_D0, stateMPC_dvaff0, stateMPC_grad_eq0);
stateMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C2, stateMPC_dvaff2, stateMPC_D1, stateMPC_dvaff1, stateMPC_grad_eq1);
stateMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C3, stateMPC_dvaff3, stateMPC_D1, stateMPC_dvaff2, stateMPC_grad_eq2);
stateMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C4, stateMPC_dvaff4, stateMPC_D1, stateMPC_dvaff3, stateMPC_grad_eq3);
stateMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C5, stateMPC_dvaff5, stateMPC_D1, stateMPC_dvaff4, stateMPC_grad_eq4);
stateMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C6, stateMPC_dvaff6, stateMPC_D1, stateMPC_dvaff5, stateMPC_grad_eq5);
stateMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C7, stateMPC_dvaff7, stateMPC_D1, stateMPC_dvaff6, stateMPC_grad_eq6);
stateMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C8, stateMPC_dvaff8, stateMPC_D1, stateMPC_dvaff7, stateMPC_grad_eq7);
stateMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C9, stateMPC_dvaff9, stateMPC_D1, stateMPC_dvaff8, stateMPC_grad_eq8);
stateMPC_LA_DIAGZERO_MTVM_2_2(stateMPC_D9, stateMPC_dvaff9, stateMPC_grad_eq9);
stateMPC_LA_VSUB2_38(stateMPC_rd, stateMPC_grad_eq, stateMPC_rd);
stateMPC_LA_DIAG_FORWARDBACKWARDSUB_4(stateMPC_Phi0, stateMPC_rd0, stateMPC_dzaff0);
stateMPC_LA_DIAG_FORWARDBACKWARDSUB_4(stateMPC_Phi1, stateMPC_rd1, stateMPC_dzaff1);
stateMPC_LA_DIAG_FORWARDBACKWARDSUB_4(stateMPC_Phi2, stateMPC_rd2, stateMPC_dzaff2);
stateMPC_LA_DIAG_FORWARDBACKWARDSUB_4(stateMPC_Phi3, stateMPC_rd3, stateMPC_dzaff3);
stateMPC_LA_DIAG_FORWARDBACKWARDSUB_4(stateMPC_Phi4, stateMPC_rd4, stateMPC_dzaff4);
stateMPC_LA_DIAG_FORWARDBACKWARDSUB_4(stateMPC_Phi5, stateMPC_rd5, stateMPC_dzaff5);
stateMPC_LA_DIAG_FORWARDBACKWARDSUB_4(stateMPC_Phi6, stateMPC_rd6, stateMPC_dzaff6);
stateMPC_LA_DIAG_FORWARDBACKWARDSUB_4(stateMPC_Phi7, stateMPC_rd7, stateMPC_dzaff7);
stateMPC_LA_DIAG_FORWARDBACKWARDSUB_4(stateMPC_Phi8, stateMPC_rd8, stateMPC_dzaff8);
stateMPC_LA_DIAG_FORWARDBACKWARDSUB_2(stateMPC_Phi9, stateMPC_rd9, stateMPC_dzaff9);
stateMPC_LA_VSUB_INDEXED_4(stateMPC_dzaff0, stateMPC_lbIdx0, stateMPC_rilb0, stateMPC_dslbaff0);
stateMPC_LA_VSUB3_4(stateMPC_llbbyslb0, stateMPC_dslbaff0, stateMPC_llb0, stateMPC_dllbaff0);
stateMPC_LA_VSUB2_INDEXED_4(stateMPC_riub0, stateMPC_dzaff0, stateMPC_ubIdx0, stateMPC_dsubaff0);
stateMPC_LA_VSUB3_4(stateMPC_lubbysub0, stateMPC_dsubaff0, stateMPC_lub0, stateMPC_dlubaff0);
stateMPC_LA_VSUB_INDEXED_4(stateMPC_dzaff1, stateMPC_lbIdx1, stateMPC_rilb1, stateMPC_dslbaff1);
stateMPC_LA_VSUB3_4(stateMPC_llbbyslb1, stateMPC_dslbaff1, stateMPC_llb1, stateMPC_dllbaff1);
stateMPC_LA_VSUB2_INDEXED_4(stateMPC_riub1, stateMPC_dzaff1, stateMPC_ubIdx1, stateMPC_dsubaff1);
stateMPC_LA_VSUB3_4(stateMPC_lubbysub1, stateMPC_dsubaff1, stateMPC_lub1, stateMPC_dlubaff1);
stateMPC_LA_VSUB_INDEXED_4(stateMPC_dzaff2, stateMPC_lbIdx2, stateMPC_rilb2, stateMPC_dslbaff2);
stateMPC_LA_VSUB3_4(stateMPC_llbbyslb2, stateMPC_dslbaff2, stateMPC_llb2, stateMPC_dllbaff2);
stateMPC_LA_VSUB2_INDEXED_4(stateMPC_riub2, stateMPC_dzaff2, stateMPC_ubIdx2, stateMPC_dsubaff2);
stateMPC_LA_VSUB3_4(stateMPC_lubbysub2, stateMPC_dsubaff2, stateMPC_lub2, stateMPC_dlubaff2);
stateMPC_LA_VSUB_INDEXED_4(stateMPC_dzaff3, stateMPC_lbIdx3, stateMPC_rilb3, stateMPC_dslbaff3);
stateMPC_LA_VSUB3_4(stateMPC_llbbyslb3, stateMPC_dslbaff3, stateMPC_llb3, stateMPC_dllbaff3);
stateMPC_LA_VSUB2_INDEXED_4(stateMPC_riub3, stateMPC_dzaff3, stateMPC_ubIdx3, stateMPC_dsubaff3);
stateMPC_LA_VSUB3_4(stateMPC_lubbysub3, stateMPC_dsubaff3, stateMPC_lub3, stateMPC_dlubaff3);
stateMPC_LA_VSUB_INDEXED_4(stateMPC_dzaff4, stateMPC_lbIdx4, stateMPC_rilb4, stateMPC_dslbaff4);
stateMPC_LA_VSUB3_4(stateMPC_llbbyslb4, stateMPC_dslbaff4, stateMPC_llb4, stateMPC_dllbaff4);
stateMPC_LA_VSUB2_INDEXED_4(stateMPC_riub4, stateMPC_dzaff4, stateMPC_ubIdx4, stateMPC_dsubaff4);
stateMPC_LA_VSUB3_4(stateMPC_lubbysub4, stateMPC_dsubaff4, stateMPC_lub4, stateMPC_dlubaff4);
stateMPC_LA_VSUB_INDEXED_4(stateMPC_dzaff5, stateMPC_lbIdx5, stateMPC_rilb5, stateMPC_dslbaff5);
stateMPC_LA_VSUB3_4(stateMPC_llbbyslb5, stateMPC_dslbaff5, stateMPC_llb5, stateMPC_dllbaff5);
stateMPC_LA_VSUB2_INDEXED_4(stateMPC_riub5, stateMPC_dzaff5, stateMPC_ubIdx5, stateMPC_dsubaff5);
stateMPC_LA_VSUB3_4(stateMPC_lubbysub5, stateMPC_dsubaff5, stateMPC_lub5, stateMPC_dlubaff5);
stateMPC_LA_VSUB_INDEXED_4(stateMPC_dzaff6, stateMPC_lbIdx6, stateMPC_rilb6, stateMPC_dslbaff6);
stateMPC_LA_VSUB3_4(stateMPC_llbbyslb6, stateMPC_dslbaff6, stateMPC_llb6, stateMPC_dllbaff6);
stateMPC_LA_VSUB2_INDEXED_4(stateMPC_riub6, stateMPC_dzaff6, stateMPC_ubIdx6, stateMPC_dsubaff6);
stateMPC_LA_VSUB3_4(stateMPC_lubbysub6, stateMPC_dsubaff6, stateMPC_lub6, stateMPC_dlubaff6);
stateMPC_LA_VSUB_INDEXED_4(stateMPC_dzaff7, stateMPC_lbIdx7, stateMPC_rilb7, stateMPC_dslbaff7);
stateMPC_LA_VSUB3_4(stateMPC_llbbyslb7, stateMPC_dslbaff7, stateMPC_llb7, stateMPC_dllbaff7);
stateMPC_LA_VSUB2_INDEXED_4(stateMPC_riub7, stateMPC_dzaff7, stateMPC_ubIdx7, stateMPC_dsubaff7);
stateMPC_LA_VSUB3_4(stateMPC_lubbysub7, stateMPC_dsubaff7, stateMPC_lub7, stateMPC_dlubaff7);
stateMPC_LA_VSUB_INDEXED_4(stateMPC_dzaff8, stateMPC_lbIdx8, stateMPC_rilb8, stateMPC_dslbaff8);
stateMPC_LA_VSUB3_4(stateMPC_llbbyslb8, stateMPC_dslbaff8, stateMPC_llb8, stateMPC_dllbaff8);
stateMPC_LA_VSUB2_INDEXED_4(stateMPC_riub8, stateMPC_dzaff8, stateMPC_ubIdx8, stateMPC_dsubaff8);
stateMPC_LA_VSUB3_4(stateMPC_lubbysub8, stateMPC_dsubaff8, stateMPC_lub8, stateMPC_dlubaff8);
stateMPC_LA_VSUB_INDEXED_2(stateMPC_dzaff9, stateMPC_lbIdx9, stateMPC_rilb9, stateMPC_dslbaff9);
stateMPC_LA_VSUB3_2(stateMPC_llbbyslb9, stateMPC_dslbaff9, stateMPC_llb9, stateMPC_dllbaff9);
stateMPC_LA_VSUB2_INDEXED_2(stateMPC_riub9, stateMPC_dzaff9, stateMPC_ubIdx9, stateMPC_dsubaff9);
stateMPC_LA_VSUB3_2(stateMPC_lubbysub9, stateMPC_dsubaff9, stateMPC_lub9, stateMPC_dlubaff9);
info->lsit_aff = stateMPC_LINESEARCH_BACKTRACKING_AFFINE(stateMPC_l, stateMPC_s, stateMPC_dl_aff, stateMPC_ds_aff, &info->step_aff, &info->mu_aff);
if( info->lsit_aff == stateMPC_NOPROGRESS ){
exitcode = stateMPC_NOPROGRESS; break;
}
sigma_3rdroot = info->mu_aff / info->mu;
info->sigma = sigma_3rdroot*sigma_3rdroot*sigma_3rdroot;
musigma = info->mu * info->sigma;
stateMPC_LA_VSUB5_76(stateMPC_ds_aff, stateMPC_dl_aff, musigma, stateMPC_ccrhs);
stateMPC_LA_VSUB6_INDEXED_4_4_4(stateMPC_ccrhsub0, stateMPC_sub0, stateMPC_ubIdx0, stateMPC_ccrhsl0, stateMPC_slb0, stateMPC_lbIdx0, stateMPC_rd0);
stateMPC_LA_VSUB6_INDEXED_4_4_4(stateMPC_ccrhsub1, stateMPC_sub1, stateMPC_ubIdx1, stateMPC_ccrhsl1, stateMPC_slb1, stateMPC_lbIdx1, stateMPC_rd1);
stateMPC_LA_DIAG_FORWARDSUB_4(stateMPC_Phi0, stateMPC_rd0, stateMPC_Lbyrd0);
stateMPC_LA_DIAG_FORWARDSUB_4(stateMPC_Phi1, stateMPC_rd1, stateMPC_Lbyrd1);
stateMPC_LA_DIAGZERO_MVM_2(stateMPC_W0, stateMPC_Lbyrd0, stateMPC_beta0);
stateMPC_LA_DENSE_FORWARDSUB_2(stateMPC_Ld0, stateMPC_beta0, stateMPC_yy0);
stateMPC_LA_DENSE_DIAGZERO_2MVMADD_2_4_4(stateMPC_V0, stateMPC_Lbyrd0, stateMPC_W1, stateMPC_Lbyrd1, stateMPC_beta1);
stateMPC_LA_DENSE_MVMSUB1_2_2(stateMPC_Lsd1, stateMPC_yy0, stateMPC_beta1, stateMPC_bmy1);
stateMPC_LA_DENSE_FORWARDSUB_2(stateMPC_Ld1, stateMPC_bmy1, stateMPC_yy1);
stateMPC_LA_VSUB6_INDEXED_4_4_4(stateMPC_ccrhsub2, stateMPC_sub2, stateMPC_ubIdx2, stateMPC_ccrhsl2, stateMPC_slb2, stateMPC_lbIdx2, stateMPC_rd2);
stateMPC_LA_DIAG_FORWARDSUB_4(stateMPC_Phi2, stateMPC_rd2, stateMPC_Lbyrd2);
stateMPC_LA_DENSE_DIAGZERO_2MVMADD_2_4_4(stateMPC_V1, stateMPC_Lbyrd1, stateMPC_W2, stateMPC_Lbyrd2, stateMPC_beta2);
stateMPC_LA_DENSE_MVMSUB1_2_2(stateMPC_Lsd2, stateMPC_yy1, stateMPC_beta2, stateMPC_bmy2);
stateMPC_LA_DENSE_FORWARDSUB_2(stateMPC_Ld2, stateMPC_bmy2, stateMPC_yy2);
stateMPC_LA_VSUB6_INDEXED_4_4_4(stateMPC_ccrhsub3, stateMPC_sub3, stateMPC_ubIdx3, stateMPC_ccrhsl3, stateMPC_slb3, stateMPC_lbIdx3, stateMPC_rd3);
stateMPC_LA_DIAG_FORWARDSUB_4(stateMPC_Phi3, stateMPC_rd3, stateMPC_Lbyrd3);
stateMPC_LA_DENSE_DIAGZERO_2MVMADD_2_4_4(stateMPC_V2, stateMPC_Lbyrd2, stateMPC_W3, stateMPC_Lbyrd3, stateMPC_beta3);
stateMPC_LA_DENSE_MVMSUB1_2_2(stateMPC_Lsd3, stateMPC_yy2, stateMPC_beta3, stateMPC_bmy3);
stateMPC_LA_DENSE_FORWARDSUB_2(stateMPC_Ld3, stateMPC_bmy3, stateMPC_yy3);
stateMPC_LA_VSUB6_INDEXED_4_4_4(stateMPC_ccrhsub4, stateMPC_sub4, stateMPC_ubIdx4, stateMPC_ccrhsl4, stateMPC_slb4, stateMPC_lbIdx4, stateMPC_rd4);
stateMPC_LA_DIAG_FORWARDSUB_4(stateMPC_Phi4, stateMPC_rd4, stateMPC_Lbyrd4);
stateMPC_LA_DENSE_DIAGZERO_2MVMADD_2_4_4(stateMPC_V3, stateMPC_Lbyrd3, stateMPC_W4, stateMPC_Lbyrd4, stateMPC_beta4);
stateMPC_LA_DENSE_MVMSUB1_2_2(stateMPC_Lsd4, stateMPC_yy3, stateMPC_beta4, stateMPC_bmy4);
stateMPC_LA_DENSE_FORWARDSUB_2(stateMPC_Ld4, stateMPC_bmy4, stateMPC_yy4);
stateMPC_LA_VSUB6_INDEXED_4_4_4(stateMPC_ccrhsub5, stateMPC_sub5, stateMPC_ubIdx5, stateMPC_ccrhsl5, stateMPC_slb5, stateMPC_lbIdx5, stateMPC_rd5);
stateMPC_LA_DIAG_FORWARDSUB_4(stateMPC_Phi5, stateMPC_rd5, stateMPC_Lbyrd5);
stateMPC_LA_DENSE_DIAGZERO_2MVMADD_2_4_4(stateMPC_V4, stateMPC_Lbyrd4, stateMPC_W5, stateMPC_Lbyrd5, stateMPC_beta5);
stateMPC_LA_DENSE_MVMSUB1_2_2(stateMPC_Lsd5, stateMPC_yy4, stateMPC_beta5, stateMPC_bmy5);
stateMPC_LA_DENSE_FORWARDSUB_2(stateMPC_Ld5, stateMPC_bmy5, stateMPC_yy5);
stateMPC_LA_VSUB6_INDEXED_4_4_4(stateMPC_ccrhsub6, stateMPC_sub6, stateMPC_ubIdx6, stateMPC_ccrhsl6, stateMPC_slb6, stateMPC_lbIdx6, stateMPC_rd6);
stateMPC_LA_DIAG_FORWARDSUB_4(stateMPC_Phi6, stateMPC_rd6, stateMPC_Lbyrd6);
stateMPC_LA_DENSE_DIAGZERO_2MVMADD_2_4_4(stateMPC_V5, stateMPC_Lbyrd5, stateMPC_W6, stateMPC_Lbyrd6, stateMPC_beta6);
stateMPC_LA_DENSE_MVMSUB1_2_2(stateMPC_Lsd6, stateMPC_yy5, stateMPC_beta6, stateMPC_bmy6);
stateMPC_LA_DENSE_FORWARDSUB_2(stateMPC_Ld6, stateMPC_bmy6, stateMPC_yy6);
stateMPC_LA_VSUB6_INDEXED_4_4_4(stateMPC_ccrhsub7, stateMPC_sub7, stateMPC_ubIdx7, stateMPC_ccrhsl7, stateMPC_slb7, stateMPC_lbIdx7, stateMPC_rd7);
stateMPC_LA_DIAG_FORWARDSUB_4(stateMPC_Phi7, stateMPC_rd7, stateMPC_Lbyrd7);
stateMPC_LA_DENSE_DIAGZERO_2MVMADD_2_4_4(stateMPC_V6, stateMPC_Lbyrd6, stateMPC_W7, stateMPC_Lbyrd7, stateMPC_beta7);
stateMPC_LA_DENSE_MVMSUB1_2_2(stateMPC_Lsd7, stateMPC_yy6, stateMPC_beta7, stateMPC_bmy7);
stateMPC_LA_DENSE_FORWARDSUB_2(stateMPC_Ld7, stateMPC_bmy7, stateMPC_yy7);
stateMPC_LA_VSUB6_INDEXED_4_4_4(stateMPC_ccrhsub8, stateMPC_sub8, stateMPC_ubIdx8, stateMPC_ccrhsl8, stateMPC_slb8, stateMPC_lbIdx8, stateMPC_rd8);
stateMPC_LA_DIAG_FORWARDSUB_4(stateMPC_Phi8, stateMPC_rd8, stateMPC_Lbyrd8);
stateMPC_LA_DENSE_DIAGZERO_2MVMADD_2_4_4(stateMPC_V7, stateMPC_Lbyrd7, stateMPC_W8, stateMPC_Lbyrd8, stateMPC_beta8);
stateMPC_LA_DENSE_MVMSUB1_2_2(stateMPC_Lsd8, stateMPC_yy7, stateMPC_beta8, stateMPC_bmy8);
stateMPC_LA_DENSE_FORWARDSUB_2(stateMPC_Ld8, stateMPC_bmy8, stateMPC_yy8);
stateMPC_LA_VSUB6_INDEXED_2_2_2(stateMPC_ccrhsub9, stateMPC_sub9, stateMPC_ubIdx9, stateMPC_ccrhsl9, stateMPC_slb9, stateMPC_lbIdx9, stateMPC_rd9);
stateMPC_LA_DIAG_FORWARDSUB_2(stateMPC_Phi9, stateMPC_rd9, stateMPC_Lbyrd9);
stateMPC_LA_DENSE_DIAGZERO_2MVMADD_2_4_2(stateMPC_V8, stateMPC_Lbyrd8, stateMPC_W9, stateMPC_Lbyrd9, stateMPC_beta9);
stateMPC_LA_DENSE_MVMSUB1_2_2(stateMPC_Lsd9, stateMPC_yy8, stateMPC_beta9, stateMPC_bmy9);
stateMPC_LA_DENSE_FORWARDSUB_2(stateMPC_Ld9, stateMPC_bmy9, stateMPC_yy9);
stateMPC_LA_DENSE_BACKWARDSUB_2(stateMPC_Ld9, stateMPC_yy9, stateMPC_dvcc9);
stateMPC_LA_DENSE_MTVMSUB_2_2(stateMPC_Lsd9, stateMPC_dvcc9, stateMPC_yy8, stateMPC_bmy8);
stateMPC_LA_DENSE_BACKWARDSUB_2(stateMPC_Ld8, stateMPC_bmy8, stateMPC_dvcc8);
stateMPC_LA_DENSE_MTVMSUB_2_2(stateMPC_Lsd8, stateMPC_dvcc8, stateMPC_yy7, stateMPC_bmy7);
stateMPC_LA_DENSE_BACKWARDSUB_2(stateMPC_Ld7, stateMPC_bmy7, stateMPC_dvcc7);
stateMPC_LA_DENSE_MTVMSUB_2_2(stateMPC_Lsd7, stateMPC_dvcc7, stateMPC_yy6, stateMPC_bmy6);
stateMPC_LA_DENSE_BACKWARDSUB_2(stateMPC_Ld6, stateMPC_bmy6, stateMPC_dvcc6);
stateMPC_LA_DENSE_MTVMSUB_2_2(stateMPC_Lsd6, stateMPC_dvcc6, stateMPC_yy5, stateMPC_bmy5);
stateMPC_LA_DENSE_BACKWARDSUB_2(stateMPC_Ld5, stateMPC_bmy5, stateMPC_dvcc5);
stateMPC_LA_DENSE_MTVMSUB_2_2(stateMPC_Lsd5, stateMPC_dvcc5, stateMPC_yy4, stateMPC_bmy4);
stateMPC_LA_DENSE_BACKWARDSUB_2(stateMPC_Ld4, stateMPC_bmy4, stateMPC_dvcc4);
stateMPC_LA_DENSE_MTVMSUB_2_2(stateMPC_Lsd4, stateMPC_dvcc4, stateMPC_yy3, stateMPC_bmy3);
stateMPC_LA_DENSE_BACKWARDSUB_2(stateMPC_Ld3, stateMPC_bmy3, stateMPC_dvcc3);
stateMPC_LA_DENSE_MTVMSUB_2_2(stateMPC_Lsd3, stateMPC_dvcc3, stateMPC_yy2, stateMPC_bmy2);
stateMPC_LA_DENSE_BACKWARDSUB_2(stateMPC_Ld2, stateMPC_bmy2, stateMPC_dvcc2);
stateMPC_LA_DENSE_MTVMSUB_2_2(stateMPC_Lsd2, stateMPC_dvcc2, stateMPC_yy1, stateMPC_bmy1);
stateMPC_LA_DENSE_BACKWARDSUB_2(stateMPC_Ld1, stateMPC_bmy1, stateMPC_dvcc1);
stateMPC_LA_DENSE_MTVMSUB_2_2(stateMPC_Lsd1, stateMPC_dvcc1, stateMPC_yy0, stateMPC_bmy0);
stateMPC_LA_DENSE_BACKWARDSUB_2(stateMPC_Ld0, stateMPC_bmy0, stateMPC_dvcc0);
stateMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C1, stateMPC_dvcc1, stateMPC_D0, stateMPC_dvcc0, stateMPC_grad_eq0);
stateMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C2, stateMPC_dvcc2, stateMPC_D1, stateMPC_dvcc1, stateMPC_grad_eq1);
stateMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C3, stateMPC_dvcc3, stateMPC_D1, stateMPC_dvcc2, stateMPC_grad_eq2);
stateMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C4, stateMPC_dvcc4, stateMPC_D1, stateMPC_dvcc3, stateMPC_grad_eq3);
stateMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C5, stateMPC_dvcc5, stateMPC_D1, stateMPC_dvcc4, stateMPC_grad_eq4);
stateMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C6, stateMPC_dvcc6, stateMPC_D1, stateMPC_dvcc5, stateMPC_grad_eq5);
stateMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C7, stateMPC_dvcc7, stateMPC_D1, stateMPC_dvcc6, stateMPC_grad_eq6);
stateMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C8, stateMPC_dvcc8, stateMPC_D1, stateMPC_dvcc7, stateMPC_grad_eq7);
stateMPC_LA_DENSE_DIAGZERO_MTVM2_2_4_2(params->C9, stateMPC_dvcc9, stateMPC_D1, stateMPC_dvcc8, stateMPC_grad_eq8);
stateMPC_LA_DIAGZERO_MTVM_2_2(stateMPC_D9, stateMPC_dvcc9, stateMPC_grad_eq9);
stateMPC_LA_VSUB_38(stateMPC_rd, stateMPC_grad_eq, stateMPC_rd);
stateMPC_LA_DIAG_FORWARDBACKWARDSUB_4(stateMPC_Phi0, stateMPC_rd0, stateMPC_dzcc0);
stateMPC_LA_DIAG_FORWARDBACKWARDSUB_4(stateMPC_Phi1, stateMPC_rd1, stateMPC_dzcc1);
stateMPC_LA_DIAG_FORWARDBACKWARDSUB_4(stateMPC_Phi2, stateMPC_rd2, stateMPC_dzcc2);
stateMPC_LA_DIAG_FORWARDBACKWARDSUB_4(stateMPC_Phi3, stateMPC_rd3, stateMPC_dzcc3);
stateMPC_LA_DIAG_FORWARDBACKWARDSUB_4(stateMPC_Phi4, stateMPC_rd4, stateMPC_dzcc4);
stateMPC_LA_DIAG_FORWARDBACKWARDSUB_4(stateMPC_Phi5, stateMPC_rd5, stateMPC_dzcc5);
stateMPC_LA_DIAG_FORWARDBACKWARDSUB_4(stateMPC_Phi6, stateMPC_rd6, stateMPC_dzcc6);
stateMPC_LA_DIAG_FORWARDBACKWARDSUB_4(stateMPC_Phi7, stateMPC_rd7, stateMPC_dzcc7);
stateMPC_LA_DIAG_FORWARDBACKWARDSUB_4(stateMPC_Phi8, stateMPC_rd8, stateMPC_dzcc8);
stateMPC_LA_DIAG_FORWARDBACKWARDSUB_2(stateMPC_Phi9, stateMPC_rd9, stateMPC_dzcc9);
stateMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_4(stateMPC_ccrhsl0, stateMPC_slb0, stateMPC_llbbyslb0, stateMPC_dzcc0, stateMPC_lbIdx0, stateMPC_dllbcc0);
stateMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_4(stateMPC_ccrhsub0, stateMPC_sub0, stateMPC_lubbysub0, stateMPC_dzcc0, stateMPC_ubIdx0, stateMPC_dlubcc0);
stateMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_4(stateMPC_ccrhsl1, stateMPC_slb1, stateMPC_llbbyslb1, stateMPC_dzcc1, stateMPC_lbIdx1, stateMPC_dllbcc1);
stateMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_4(stateMPC_ccrhsub1, stateMPC_sub1, stateMPC_lubbysub1, stateMPC_dzcc1, stateMPC_ubIdx1, stateMPC_dlubcc1);
stateMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_4(stateMPC_ccrhsl2, stateMPC_slb2, stateMPC_llbbyslb2, stateMPC_dzcc2, stateMPC_lbIdx2, stateMPC_dllbcc2);
stateMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_4(stateMPC_ccrhsub2, stateMPC_sub2, stateMPC_lubbysub2, stateMPC_dzcc2, stateMPC_ubIdx2, stateMPC_dlubcc2);
stateMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_4(stateMPC_ccrhsl3, stateMPC_slb3, stateMPC_llbbyslb3, stateMPC_dzcc3, stateMPC_lbIdx3, stateMPC_dllbcc3);
stateMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_4(stateMPC_ccrhsub3, stateMPC_sub3, stateMPC_lubbysub3, stateMPC_dzcc3, stateMPC_ubIdx3, stateMPC_dlubcc3);
stateMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_4(stateMPC_ccrhsl4, stateMPC_slb4, stateMPC_llbbyslb4, stateMPC_dzcc4, stateMPC_lbIdx4, stateMPC_dllbcc4);
stateMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_4(stateMPC_ccrhsub4, stateMPC_sub4, stateMPC_lubbysub4, stateMPC_dzcc4, stateMPC_ubIdx4, stateMPC_dlubcc4);
stateMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_4(stateMPC_ccrhsl5, stateMPC_slb5, stateMPC_llbbyslb5, stateMPC_dzcc5, stateMPC_lbIdx5, stateMPC_dllbcc5);
stateMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_4(stateMPC_ccrhsub5, stateMPC_sub5, stateMPC_lubbysub5, stateMPC_dzcc5, stateMPC_ubIdx5, stateMPC_dlubcc5);
stateMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_4(stateMPC_ccrhsl6, stateMPC_slb6, stateMPC_llbbyslb6, stateMPC_dzcc6, stateMPC_lbIdx6, stateMPC_dllbcc6);
stateMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_4(stateMPC_ccrhsub6, stateMPC_sub6, stateMPC_lubbysub6, stateMPC_dzcc6, stateMPC_ubIdx6, stateMPC_dlubcc6);
stateMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_4(stateMPC_ccrhsl7, stateMPC_slb7, stateMPC_llbbyslb7, stateMPC_dzcc7, stateMPC_lbIdx7, stateMPC_dllbcc7);
stateMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_4(stateMPC_ccrhsub7, stateMPC_sub7, stateMPC_lubbysub7, stateMPC_dzcc7, stateMPC_ubIdx7, stateMPC_dlubcc7);
stateMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_4(stateMPC_ccrhsl8, stateMPC_slb8, stateMPC_llbbyslb8, stateMPC_dzcc8, stateMPC_lbIdx8, stateMPC_dllbcc8);
stateMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_4(stateMPC_ccrhsub8, stateMPC_sub8, stateMPC_lubbysub8, stateMPC_dzcc8, stateMPC_ubIdx8, stateMPC_dlubcc8);
stateMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_2(stateMPC_ccrhsl9, stateMPC_slb9, stateMPC_llbbyslb9, stateMPC_dzcc9, stateMPC_lbIdx9, stateMPC_dllbcc9);
stateMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_2(stateMPC_ccrhsub9, stateMPC_sub9, stateMPC_lubbysub9, stateMPC_dzcc9, stateMPC_ubIdx9, stateMPC_dlubcc9);
stateMPC_LA_VSUB7_76(stateMPC_l, stateMPC_ccrhs, stateMPC_s, stateMPC_dl_cc, stateMPC_ds_cc);
stateMPC_LA_VADD_38(stateMPC_dz_cc, stateMPC_dz_aff);
stateMPC_LA_VADD_20(stateMPC_dv_cc, stateMPC_dv_aff);
stateMPC_LA_VADD_76(stateMPC_dl_cc, stateMPC_dl_aff);
stateMPC_LA_VADD_76(stateMPC_ds_cc, stateMPC_ds_aff);
info->lsit_cc = stateMPC_LINESEARCH_BACKTRACKING_COMBINED(stateMPC_z, stateMPC_v, stateMPC_l, stateMPC_s, stateMPC_dz_cc, stateMPC_dv_cc, stateMPC_dl_cc, stateMPC_ds_cc, &info->step_cc, &info->mu);
if( info->lsit_cc == stateMPC_NOPROGRESS ){
exitcode = stateMPC_NOPROGRESS; break;
}
info->it++;
}
output->z1[0] = stateMPC_z0[0];
output->z1[1] = stateMPC_z0[1];
output->z1[2] = stateMPC_z0[2];
output->z1[3] = stateMPC_z0[3];
output->z2[0] = stateMPC_z1[0];
output->z2[1] = stateMPC_z1[1];
output->z2[2] = stateMPC_z1[2];
output->z2[3] = stateMPC_z1[3];
output->z3[0] = stateMPC_z2[0];
output->z3[1] = stateMPC_z2[1];
output->z3[2] = stateMPC_z2[2];
output->z3[3] = stateMPC_z2[3];
output->z4[0] = stateMPC_z3[0];
output->z4[1] = stateMPC_z3[1];
output->z4[2] = stateMPC_z3[2];
output->z4[3] = stateMPC_z3[3];
output->z5[0] = stateMPC_z4[0];
output->z5[1] = stateMPC_z4[1];
output->z5[2] = stateMPC_z4[2];
output->z5[3] = stateMPC_z4[3];
output->z6[0] = stateMPC_z5[0];
output->z6[1] = stateMPC_z5[1];
output->z6[2] = stateMPC_z5[2];
output->z6[3] = stateMPC_z5[3];
output->z7[0] = stateMPC_z6[0];
output->z7[1] = stateMPC_z6[1];
output->z7[2] = stateMPC_z6[2];
output->z7[3] = stateMPC_z6[3];
output->z8[0] = stateMPC_z7[0];
output->z8[1] = stateMPC_z7[1];
output->z8[2] = stateMPC_z7[2];
output->z8[3] = stateMPC_z7[3];
output->z9[0] = stateMPC_z8[0];
output->z9[1] = stateMPC_z8[1];
output->z9[2] = stateMPC_z8[2];
output->z9[3] = stateMPC_z8[3];
output->z10[0] = stateMPC_z9[0];
output->z10[1] = stateMPC_z9[1];

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
