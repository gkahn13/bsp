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

#include "planarMPC.h"

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
 * Initializes a vector of length 63 with a value.
 */
void planarMPC_LA_INITIALIZEVECTOR_63(planarMPC_FLOAT* vec, planarMPC_FLOAT value)
{
	int i;
	for( i=0; i<63; i++ )
	{
		vec[i] = value;
	}
}


/*
 * Initializes a vector of length 35 with a value.
 */
void planarMPC_LA_INITIALIZEVECTOR_35(planarMPC_FLOAT* vec, planarMPC_FLOAT value)
{
	int i;
	for( i=0; i<35; i++ )
	{
		vec[i] = value;
	}
}


/*
 * Initializes a vector of length 126 with a value.
 */
void planarMPC_LA_INITIALIZEVECTOR_126(planarMPC_FLOAT* vec, planarMPC_FLOAT value)
{
	int i;
	for( i=0; i<126; i++ )
	{
		vec[i] = value;
	}
}


/* 
 * Calculates a dot product and adds it to a variable: z += x'*y; 
 * This function is for vectors of length 126.
 */
void planarMPC_LA_DOTACC_126(planarMPC_FLOAT *x, planarMPC_FLOAT *y, planarMPC_FLOAT *z)
{
	int i;
	for( i=0; i<126; i++ ){
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
void planarMPC_LA_DIAG_QUADFCN_14(planarMPC_FLOAT* H, planarMPC_FLOAT* f, planarMPC_FLOAT* z, planarMPC_FLOAT* grad, planarMPC_FLOAT* value)
{
	int i;
	planarMPC_FLOAT hz;	
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
void planarMPC_LA_DIAG_QUADFCN_7(planarMPC_FLOAT* H, planarMPC_FLOAT* f, planarMPC_FLOAT* z, planarMPC_FLOAT* grad, planarMPC_FLOAT* value)
{
	int i;
	planarMPC_FLOAT hz;	
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
void planarMPC_LA_DIAGZERO_MVMSUB6_7(planarMPC_FLOAT *B, planarMPC_FLOAT *u, planarMPC_FLOAT *b, planarMPC_FLOAT *l, planarMPC_FLOAT *r, planarMPC_FLOAT *z, planarMPC_FLOAT *y)
{
	int i;
	planarMPC_FLOAT Bu[7];
	planarMPC_FLOAT norm = *y;
	planarMPC_FLOAT lr = 0;

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
void planarMPC_LA_DENSE_DIAGZERO_MVMSUB3_7_14_14(planarMPC_FLOAT *A, planarMPC_FLOAT *x, planarMPC_FLOAT *B, planarMPC_FLOAT *u, planarMPC_FLOAT *b, planarMPC_FLOAT *l, planarMPC_FLOAT *r, planarMPC_FLOAT *z, planarMPC_FLOAT *y)
{
	int i;
	int j;
	int k = 0;
	planarMPC_FLOAT AxBu[7];
	planarMPC_FLOAT norm = *y;
	planarMPC_FLOAT lr = 0;

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
void planarMPC_LA_DENSE_DIAGZERO_MVMSUB3_7_14_7(planarMPC_FLOAT *A, planarMPC_FLOAT *x, planarMPC_FLOAT *B, planarMPC_FLOAT *u, planarMPC_FLOAT *b, planarMPC_FLOAT *l, planarMPC_FLOAT *r, planarMPC_FLOAT *z, planarMPC_FLOAT *y)
{
	int i;
	int j;
	int k = 0;
	planarMPC_FLOAT AxBu[7];
	planarMPC_FLOAT norm = *y;
	planarMPC_FLOAT lr = 0;

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
void planarMPC_LA_DENSE_DIAGZERO_MTVM2_7_14_7(planarMPC_FLOAT *A, planarMPC_FLOAT *x, planarMPC_FLOAT *B, planarMPC_FLOAT *y, planarMPC_FLOAT *z)
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
void planarMPC_LA_DIAGZERO_MTVM_7_7(planarMPC_FLOAT *M, planarMPC_FLOAT *x, planarMPC_FLOAT *y)
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
void planarMPC_LA_VSUBADD3_14(planarMPC_FLOAT* t, planarMPC_FLOAT* u, int* uidx, planarMPC_FLOAT* v, planarMPC_FLOAT* w, planarMPC_FLOAT* y, planarMPC_FLOAT* z, planarMPC_FLOAT* r)
{
	int i;
	planarMPC_FLOAT norm = *r;
	planarMPC_FLOAT vx = 0;
	planarMPC_FLOAT x;
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
void planarMPC_LA_VSUBADD2_14(planarMPC_FLOAT* t, int* tidx, planarMPC_FLOAT* u, planarMPC_FLOAT* v, planarMPC_FLOAT* w, planarMPC_FLOAT* y, planarMPC_FLOAT* z, planarMPC_FLOAT* r)
{
	int i;
	planarMPC_FLOAT norm = *r;
	planarMPC_FLOAT vx = 0;
	planarMPC_FLOAT x;
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
void planarMPC_LA_VSUBADD3_7(planarMPC_FLOAT* t, planarMPC_FLOAT* u, int* uidx, planarMPC_FLOAT* v, planarMPC_FLOAT* w, planarMPC_FLOAT* y, planarMPC_FLOAT* z, planarMPC_FLOAT* r)
{
	int i;
	planarMPC_FLOAT norm = *r;
	planarMPC_FLOAT vx = 0;
	planarMPC_FLOAT x;
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
void planarMPC_LA_VSUBADD2_7(planarMPC_FLOAT* t, int* tidx, planarMPC_FLOAT* u, planarMPC_FLOAT* v, planarMPC_FLOAT* w, planarMPC_FLOAT* y, planarMPC_FLOAT* z, planarMPC_FLOAT* r)
{
	int i;
	planarMPC_FLOAT norm = *r;
	planarMPC_FLOAT vx = 0;
	planarMPC_FLOAT x;
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
void planarMPC_LA_INEQ_B_GRAD_14_14_14(planarMPC_FLOAT *lu, planarMPC_FLOAT *su, planarMPC_FLOAT *ru, planarMPC_FLOAT *ll, planarMPC_FLOAT *sl, planarMPC_FLOAT *rl, int* lbIdx, int* ubIdx, planarMPC_FLOAT *grad, planarMPC_FLOAT *lubysu, planarMPC_FLOAT *llbysl)
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
void planarMPC_LA_INEQ_B_GRAD_7_7_7(planarMPC_FLOAT *lu, planarMPC_FLOAT *su, planarMPC_FLOAT *ru, planarMPC_FLOAT *ll, planarMPC_FLOAT *sl, planarMPC_FLOAT *rl, int* lbIdx, int* ubIdx, planarMPC_FLOAT *grad, planarMPC_FLOAT *lubysu, planarMPC_FLOAT *llbysl)
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
 * of length 63.
 */
void planarMPC_LA_VVADD3_63(planarMPC_FLOAT *u, planarMPC_FLOAT *v, planarMPC_FLOAT *w, planarMPC_FLOAT *z)
{
	int i;
	for( i=0; i<63; i++ ){
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
void planarMPC_LA_DIAG_CHOL_ONELOOP_LBUB_14_14_14(planarMPC_FLOAT *H, planarMPC_FLOAT *llbysl, int* lbIdx, planarMPC_FLOAT *lubysu, int* ubIdx, planarMPC_FLOAT *Phi)


{
	int i;
	
	/* compute cholesky */
	for( i=0; i<14; i++ ){
		Phi[i] = H[i] + llbysl[i] + lubysu[i];

#if planarMPC_SET_PRINTLEVEL > 0 && defined PRINTNUMERICALWARNINGS
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
void planarMPC_LA_DIAG_MATRIXFORWARDSUB_7_14(planarMPC_FLOAT *L, planarMPC_FLOAT *B, planarMPC_FLOAT *A)
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
void planarMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_7_14(planarMPC_FLOAT *L, planarMPC_FLOAT *B, planarMPC_FLOAT *A)
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
void planarMPC_LA_DENSE_DIAGZERO_MMTM_7_14_7(planarMPC_FLOAT *A, planarMPC_FLOAT *B, planarMPC_FLOAT *C)
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
void planarMPC_LA_DIAG_FORWARDSUB_14(planarMPC_FLOAT *L, planarMPC_FLOAT *b, planarMPC_FLOAT *y)
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
void planarMPC_LA_DIAG_CHOL_ONELOOP_LBUB_7_7_7(planarMPC_FLOAT *H, planarMPC_FLOAT *llbysl, int* lbIdx, planarMPC_FLOAT *lubysu, int* ubIdx, planarMPC_FLOAT *Phi)


{
	int i;
	
	/* compute cholesky */
	for( i=0; i<7; i++ ){
		Phi[i] = H[i] + llbysl[i] + lubysu[i];

#if planarMPC_SET_PRINTLEVEL > 0 && defined PRINTNUMERICALWARNINGS
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
void planarMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_7_7(planarMPC_FLOAT *L, planarMPC_FLOAT *B, planarMPC_FLOAT *A)
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
void planarMPC_LA_DIAG_FORWARDSUB_7(planarMPC_FLOAT *L, planarMPC_FLOAT *b, planarMPC_FLOAT *y)
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
void planarMPC_LA_DIAGZERO_MMT_7(planarMPC_FLOAT *B, planarMPC_FLOAT *L)
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
void planarMPC_LA_DIAGZERO_MVMSUB7_7(planarMPC_FLOAT *B, planarMPC_FLOAT *u, planarMPC_FLOAT *b, planarMPC_FLOAT *r)
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
void planarMPC_LA_DENSE_DIAGZERO_MMT2_7_14_14(planarMPC_FLOAT *A, planarMPC_FLOAT *B, planarMPC_FLOAT *L)
{
    int i, j, k, ii, di;
    planarMPC_FLOAT ltemp;
    
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
void planarMPC_LA_DENSE_DIAGZERO_2MVMSUB2_7_14_14(planarMPC_FLOAT *A, planarMPC_FLOAT *x, planarMPC_FLOAT *B, planarMPC_FLOAT *u, planarMPC_FLOAT *b, planarMPC_FLOAT *r)
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
void planarMPC_LA_DENSE_DIAGZERO_MMT2_7_14_7(planarMPC_FLOAT *A, planarMPC_FLOAT *B, planarMPC_FLOAT *L)
{
    int i, j, k, ii, di;
    planarMPC_FLOAT ltemp;
    
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
void planarMPC_LA_DENSE_DIAGZERO_2MVMSUB2_7_14_7(planarMPC_FLOAT *A, planarMPC_FLOAT *x, planarMPC_FLOAT *B, planarMPC_FLOAT *u, planarMPC_FLOAT *b, planarMPC_FLOAT *r)
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
void planarMPC_LA_DENSE_CHOL_7(planarMPC_FLOAT *A, planarMPC_FLOAT *L)
{
    int i, j, k, di, dj;
	 int ii, jj;

    planarMPC_FLOAT l;
    planarMPC_FLOAT Mii;

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
        
#if planarMPC_SET_PRINTLEVEL > 0 && defined PRINTNUMERICALWARNINGS
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
void planarMPC_LA_DENSE_FORWARDSUB_7(planarMPC_FLOAT *L, planarMPC_FLOAT *b, planarMPC_FLOAT *y)
{
    int i,j,ii,di;
    planarMPC_FLOAT yel;
            
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
void planarMPC_LA_DENSE_MATRIXTFORWARDSUB_7_7(planarMPC_FLOAT *L, planarMPC_FLOAT *B, planarMPC_FLOAT *A)
{
    int i,j,k,ii,di;
    planarMPC_FLOAT a;
    
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
void planarMPC_LA_DENSE_MMTSUB_7_7(planarMPC_FLOAT *A, planarMPC_FLOAT *L)
{
    int i, j, k, ii, di;
    planarMPC_FLOAT ltemp;
    
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
void planarMPC_LA_DENSE_MVMSUB1_7_7(planarMPC_FLOAT *A, planarMPC_FLOAT *x, planarMPC_FLOAT *b, planarMPC_FLOAT *r)
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
void planarMPC_LA_DENSE_BACKWARDSUB_7(planarMPC_FLOAT *L, planarMPC_FLOAT *y, planarMPC_FLOAT *x)
{
    int i, ii, di, j, jj, dj;
    planarMPC_FLOAT xel;    
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
void planarMPC_LA_DENSE_MTVMSUB_7_7(planarMPC_FLOAT *A, planarMPC_FLOAT *x, planarMPC_FLOAT *b, planarMPC_FLOAT *r)
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
 * Vector subtraction z = -x - y for vectors of length 63.
 */
void planarMPC_LA_VSUB2_63(planarMPC_FLOAT *x, planarMPC_FLOAT *y, planarMPC_FLOAT *z)
{
	int i;
	for( i=0; i<63; i++){
		z[i] = -x[i] - y[i];
	}
}


/**
 * Forward-Backward-Substitution to solve L*L^T*x = b where L is a
 * diagonal matrix of size 14 in vector
 * storage format.
 */
void planarMPC_LA_DIAG_FORWARDBACKWARDSUB_14(planarMPC_FLOAT *L, planarMPC_FLOAT *b, planarMPC_FLOAT *x)
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
void planarMPC_LA_DIAG_FORWARDBACKWARDSUB_7(planarMPC_FLOAT *L, planarMPC_FLOAT *b, planarMPC_FLOAT *x)
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
void planarMPC_LA_VSUB_INDEXED_14(planarMPC_FLOAT *x, int* xidx, planarMPC_FLOAT *y, planarMPC_FLOAT *z)
{
	int i;
	for( i=0; i<14; i++){
		z[i] = x[xidx[i]] - y[i];
	}
}


/*
 * Vector subtraction x = -u.*v - w for vectors of length 14.
 */
void planarMPC_LA_VSUB3_14(planarMPC_FLOAT *u, planarMPC_FLOAT *v, planarMPC_FLOAT *w, planarMPC_FLOAT *x)
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
void planarMPC_LA_VSUB2_INDEXED_14(planarMPC_FLOAT *x, planarMPC_FLOAT *y, int* yidx, planarMPC_FLOAT *z)
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
void planarMPC_LA_VSUB_INDEXED_7(planarMPC_FLOAT *x, int* xidx, planarMPC_FLOAT *y, planarMPC_FLOAT *z)
{
	int i;
	for( i=0; i<7; i++){
		z[i] = x[xidx[i]] - y[i];
	}
}


/*
 * Vector subtraction x = -u.*v - w for vectors of length 7.
 */
void planarMPC_LA_VSUB3_7(planarMPC_FLOAT *u, planarMPC_FLOAT *v, planarMPC_FLOAT *w, planarMPC_FLOAT *x)
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
void planarMPC_LA_VSUB2_INDEXED_7(planarMPC_FLOAT *x, planarMPC_FLOAT *y, int* yidx, planarMPC_FLOAT *z)
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
 * planarMPC_NOPROGRESS (should be negative).
 */
int planarMPC_LINESEARCH_BACKTRACKING_AFFINE(planarMPC_FLOAT *l, planarMPC_FLOAT *s, planarMPC_FLOAT *dl, planarMPC_FLOAT *ds, planarMPC_FLOAT *a, planarMPC_FLOAT *mu_aff)
{
    int i;
	int lsIt=1;    
    planarMPC_FLOAT dltemp;
    planarMPC_FLOAT dstemp;
    planarMPC_FLOAT mya = 1.0;
    planarMPC_FLOAT mymu;
        
    while( 1 ){                        

        /* 
         * Compute both snew and wnew together.
         * We compute also mu_affine along the way here, as the
         * values might be in registers, so it should be cheaper.
         */
        mymu = 0;
        for( i=0; i<126; i++ ){
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
        if( i == 126 ){
            break;
        } else {
            mya *= planarMPC_SET_LS_SCALE_AFF;
            if( mya < planarMPC_SET_LS_MINSTEP ){
                return planarMPC_NOPROGRESS;
            }
        }
    }
    
    /* return new values and iteration counter */
    *a = mya;
    *mu_aff = mymu / (planarMPC_FLOAT)126;
    return lsIt;
}


/*
 * Vector subtraction x = (u.*v - mu)*sigma where a is a scalar
*  and x,u,v are vectors of length 126.
 */
void planarMPC_LA_VSUB5_126(planarMPC_FLOAT *u, planarMPC_FLOAT *v, planarMPC_FLOAT mu,  planarMPC_FLOAT sigma, planarMPC_FLOAT *x)
{
	int i;
	for( i=0; i<126; i++){
		x[i] = u[i]*v[i] - mu;
		x[i] *= sigma;
	}
}


/*
 * Computes x=0; x(uidx) += u/su; x(vidx) -= v/sv where x is of length 14,
 * u, su, uidx are of length 14 and v, sv, vidx are of length 14.
 */
void planarMPC_LA_VSUB6_INDEXED_14_14_14(planarMPC_FLOAT *u, planarMPC_FLOAT *su, int* uidx, planarMPC_FLOAT *v, planarMPC_FLOAT *sv, int* vidx, planarMPC_FLOAT *x)
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
void planarMPC_LA_DIAGZERO_MVM_7(planarMPC_FLOAT *B, planarMPC_FLOAT *u, planarMPC_FLOAT *r)
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
void planarMPC_LA_DENSE_DIAGZERO_2MVMADD_7_14_14(planarMPC_FLOAT *A, planarMPC_FLOAT *x, planarMPC_FLOAT *B, planarMPC_FLOAT *u, planarMPC_FLOAT *r)
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
void planarMPC_LA_VSUB6_INDEXED_7_7_7(planarMPC_FLOAT *u, planarMPC_FLOAT *su, int* uidx, planarMPC_FLOAT *v, planarMPC_FLOAT *sv, int* vidx, planarMPC_FLOAT *x)
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
void planarMPC_LA_DENSE_DIAGZERO_2MVMADD_7_14_7(planarMPC_FLOAT *A, planarMPC_FLOAT *x, planarMPC_FLOAT *B, planarMPC_FLOAT *u, planarMPC_FLOAT *r)
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
 * Vector subtraction z = x - y for vectors of length 63.
 */
void planarMPC_LA_VSUB_63(planarMPC_FLOAT *x, planarMPC_FLOAT *y, planarMPC_FLOAT *z)
{
	int i;
	for( i=0; i<63; i++){
		z[i] = x[i] - y[i];
	}
}


/** 
 * Computes z = -r./s - u.*y(y)
 * where all vectors except of y are of length 14 (length of y >= 14).
 */
void planarMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_14(planarMPC_FLOAT *r, planarMPC_FLOAT *s, planarMPC_FLOAT *u, planarMPC_FLOAT *y, int* yidx, planarMPC_FLOAT *z)
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
void planarMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_14(planarMPC_FLOAT *r, planarMPC_FLOAT *s, planarMPC_FLOAT *u, planarMPC_FLOAT *y, int* yidx, planarMPC_FLOAT *z)
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
void planarMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_7(planarMPC_FLOAT *r, planarMPC_FLOAT *s, planarMPC_FLOAT *u, planarMPC_FLOAT *y, int* yidx, planarMPC_FLOAT *z)
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
void planarMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_7(planarMPC_FLOAT *r, planarMPC_FLOAT *s, planarMPC_FLOAT *u, planarMPC_FLOAT *y, int* yidx, planarMPC_FLOAT *z)
{
	int i;
	for( i=0; i<7; i++ ){
		z[i] = -r[i]/s[i] + u[i]*y[yidx[i]];
	}
}


/*
 * Computes ds = -l.\(r + s.*dl) for vectors of length 126.
 */
void planarMPC_LA_VSUB7_126(planarMPC_FLOAT *l, planarMPC_FLOAT *r, planarMPC_FLOAT *s, planarMPC_FLOAT *dl, planarMPC_FLOAT *ds)
{
	int i;
	for( i=0; i<126; i++){
		ds[i] = -(r[i] + s[i]*dl[i])/l[i];
	}
}


/*
 * Vector addition x = x + y for vectors of length 63.
 */
void planarMPC_LA_VADD_63(planarMPC_FLOAT *x, planarMPC_FLOAT *y)
{
	int i;
	for( i=0; i<63; i++){
		x[i] += y[i];
	}
}


/*
 * Vector addition x = x + y for vectors of length 35.
 */
void planarMPC_LA_VADD_35(planarMPC_FLOAT *x, planarMPC_FLOAT *y)
{
	int i;
	for( i=0; i<35; i++){
		x[i] += y[i];
	}
}


/*
 * Vector addition x = x + y for vectors of length 126.
 */
void planarMPC_LA_VADD_126(planarMPC_FLOAT *x, planarMPC_FLOAT *y)
{
	int i;
	for( i=0; i<126; i++){
		x[i] += y[i];
	}
}


/**
 * Backtracking line search for combined predictor/corrector step.
 * Update on variables with safety factor gamma (to keep us away from
 * boundary).
 */
int planarMPC_LINESEARCH_BACKTRACKING_COMBINED(planarMPC_FLOAT *z, planarMPC_FLOAT *v, planarMPC_FLOAT *l, planarMPC_FLOAT *s, planarMPC_FLOAT *dz, planarMPC_FLOAT *dv, planarMPC_FLOAT *dl, planarMPC_FLOAT *ds, planarMPC_FLOAT *a, planarMPC_FLOAT *mu)
{
    int i, lsIt=1;       
    planarMPC_FLOAT dltemp;
    planarMPC_FLOAT dstemp;    
    planarMPC_FLOAT a_gamma;
            
    *a = 1.0;
    while( 1 ){                        

        /* check whether search criterion is fulfilled */
        for( i=0; i<126; i++ ){
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
        if( i == 126 ){
            break;
        } else {
            *a *= planarMPC_SET_LS_SCALE;
            if( *a < planarMPC_SET_LS_MINSTEP ){
                return planarMPC_NOPROGRESS;
            }
        }
    }
    
    /* update variables with safety margin */
    a_gamma = (*a)*planarMPC_SET_LS_MAXSTEP;
    
    /* primal variables */
    for( i=0; i<63; i++ ){
        z[i] += a_gamma*dz[i];
    }
    
    /* equality constraint multipliers */
    for( i=0; i<35; i++ ){
        v[i] += a_gamma*dv[i];
    }
    
    /* inequality constraint multipliers & slacks, also update mu */
    *mu = 0;
    for( i=0; i<126; i++ ){
        dltemp = l[i] + a_gamma*dl[i]; l[i] = dltemp;
        dstemp = s[i] + a_gamma*ds[i]; s[i] = dstemp;
        *mu += dltemp*dstemp;
    }
    
    *a = a_gamma;
    *mu /= (planarMPC_FLOAT)126;
    return lsIt;
}




/* VARIABLE DEFINITIONS ------------------------------------------------ */
planarMPC_FLOAT planarMPC_z[63];
planarMPC_FLOAT planarMPC_v[35];
planarMPC_FLOAT planarMPC_dz_aff[63];
planarMPC_FLOAT planarMPC_dv_aff[35];
planarMPC_FLOAT planarMPC_grad_cost[63];
planarMPC_FLOAT planarMPC_grad_eq[63];
planarMPC_FLOAT planarMPC_rd[63];
planarMPC_FLOAT planarMPC_l[126];
planarMPC_FLOAT planarMPC_s[126];
planarMPC_FLOAT planarMPC_lbys[126];
planarMPC_FLOAT planarMPC_dl_aff[126];
planarMPC_FLOAT planarMPC_ds_aff[126];
planarMPC_FLOAT planarMPC_dz_cc[63];
planarMPC_FLOAT planarMPC_dv_cc[35];
planarMPC_FLOAT planarMPC_dl_cc[126];
planarMPC_FLOAT planarMPC_ds_cc[126];
planarMPC_FLOAT planarMPC_ccrhs[126];
planarMPC_FLOAT planarMPC_grad_ineq[63];
planarMPC_FLOAT* planarMPC_z0 = planarMPC_z + 0;
planarMPC_FLOAT* planarMPC_dzaff0 = planarMPC_dz_aff + 0;
planarMPC_FLOAT* planarMPC_dzcc0 = planarMPC_dz_cc + 0;
planarMPC_FLOAT* planarMPC_rd0 = planarMPC_rd + 0;
planarMPC_FLOAT planarMPC_Lbyrd0[14];
planarMPC_FLOAT* planarMPC_grad_cost0 = planarMPC_grad_cost + 0;
planarMPC_FLOAT* planarMPC_grad_eq0 = planarMPC_grad_eq + 0;
planarMPC_FLOAT* planarMPC_grad_ineq0 = planarMPC_grad_ineq + 0;
planarMPC_FLOAT planarMPC_ctv0[14];
planarMPC_FLOAT planarMPC_C0[98] = {1.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 
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
planarMPC_FLOAT* planarMPC_v0 = planarMPC_v + 0;
planarMPC_FLOAT planarMPC_re0[7];
planarMPC_FLOAT planarMPC_beta0[7];
planarMPC_FLOAT planarMPC_betacc0[7];
planarMPC_FLOAT* planarMPC_dvaff0 = planarMPC_dv_aff + 0;
planarMPC_FLOAT* planarMPC_dvcc0 = planarMPC_dv_cc + 0;
planarMPC_FLOAT planarMPC_V0[98];
planarMPC_FLOAT planarMPC_Yd0[28];
planarMPC_FLOAT planarMPC_Ld0[28];
planarMPC_FLOAT planarMPC_yy0[7];
planarMPC_FLOAT planarMPC_bmy0[7];
int planarMPC_lbIdx0[14] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13};
planarMPC_FLOAT* planarMPC_llb0 = planarMPC_l + 0;
planarMPC_FLOAT* planarMPC_slb0 = planarMPC_s + 0;
planarMPC_FLOAT* planarMPC_llbbyslb0 = planarMPC_lbys + 0;
planarMPC_FLOAT planarMPC_rilb0[14];
planarMPC_FLOAT* planarMPC_dllbaff0 = planarMPC_dl_aff + 0;
planarMPC_FLOAT* planarMPC_dslbaff0 = planarMPC_ds_aff + 0;
planarMPC_FLOAT* planarMPC_dllbcc0 = planarMPC_dl_cc + 0;
planarMPC_FLOAT* planarMPC_dslbcc0 = planarMPC_ds_cc + 0;
planarMPC_FLOAT* planarMPC_ccrhsl0 = planarMPC_ccrhs + 0;
int planarMPC_ubIdx0[14] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13};
planarMPC_FLOAT* planarMPC_lub0 = planarMPC_l + 14;
planarMPC_FLOAT* planarMPC_sub0 = planarMPC_s + 14;
planarMPC_FLOAT* planarMPC_lubbysub0 = planarMPC_lbys + 14;
planarMPC_FLOAT planarMPC_riub0[14];
planarMPC_FLOAT* planarMPC_dlubaff0 = planarMPC_dl_aff + 14;
planarMPC_FLOAT* planarMPC_dsubaff0 = planarMPC_ds_aff + 14;
planarMPC_FLOAT* planarMPC_dlubcc0 = planarMPC_dl_cc + 14;
planarMPC_FLOAT* planarMPC_dsubcc0 = planarMPC_ds_cc + 14;
planarMPC_FLOAT* planarMPC_ccrhsub0 = planarMPC_ccrhs + 14;
planarMPC_FLOAT planarMPC_Phi0[14];
planarMPC_FLOAT planarMPC_D0[14] = {1.0000000000000000E+000, 
1.0000000000000000E+000, 
1.0000000000000000E+000, 
1.0000000000000000E+000, 
1.0000000000000000E+000, 
1.0000000000000000E+000, 
1.0000000000000000E+000};
planarMPC_FLOAT planarMPC_W0[14];
planarMPC_FLOAT* planarMPC_z1 = planarMPC_z + 14;
planarMPC_FLOAT* planarMPC_dzaff1 = planarMPC_dz_aff + 14;
planarMPC_FLOAT* planarMPC_dzcc1 = planarMPC_dz_cc + 14;
planarMPC_FLOAT* planarMPC_rd1 = planarMPC_rd + 14;
planarMPC_FLOAT planarMPC_Lbyrd1[14];
planarMPC_FLOAT* planarMPC_grad_cost1 = planarMPC_grad_cost + 14;
planarMPC_FLOAT* planarMPC_grad_eq1 = planarMPC_grad_eq + 14;
planarMPC_FLOAT* planarMPC_grad_ineq1 = planarMPC_grad_ineq + 14;
planarMPC_FLOAT planarMPC_ctv1[14];
planarMPC_FLOAT* planarMPC_v1 = planarMPC_v + 7;
planarMPC_FLOAT planarMPC_re1[7];
planarMPC_FLOAT planarMPC_beta1[7];
planarMPC_FLOAT planarMPC_betacc1[7];
planarMPC_FLOAT* planarMPC_dvaff1 = planarMPC_dv_aff + 7;
planarMPC_FLOAT* planarMPC_dvcc1 = planarMPC_dv_cc + 7;
planarMPC_FLOAT planarMPC_V1[98];
planarMPC_FLOAT planarMPC_Yd1[28];
planarMPC_FLOAT planarMPC_Ld1[28];
planarMPC_FLOAT planarMPC_yy1[7];
planarMPC_FLOAT planarMPC_bmy1[7];
planarMPC_FLOAT planarMPC_c1[7] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
int planarMPC_lbIdx1[14] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13};
planarMPC_FLOAT* planarMPC_llb1 = planarMPC_l + 28;
planarMPC_FLOAT* planarMPC_slb1 = planarMPC_s + 28;
planarMPC_FLOAT* planarMPC_llbbyslb1 = planarMPC_lbys + 28;
planarMPC_FLOAT planarMPC_rilb1[14];
planarMPC_FLOAT* planarMPC_dllbaff1 = planarMPC_dl_aff + 28;
planarMPC_FLOAT* planarMPC_dslbaff1 = planarMPC_ds_aff + 28;
planarMPC_FLOAT* planarMPC_dllbcc1 = planarMPC_dl_cc + 28;
planarMPC_FLOAT* planarMPC_dslbcc1 = planarMPC_ds_cc + 28;
planarMPC_FLOAT* planarMPC_ccrhsl1 = planarMPC_ccrhs + 28;
int planarMPC_ubIdx1[14] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13};
planarMPC_FLOAT* planarMPC_lub1 = planarMPC_l + 42;
planarMPC_FLOAT* planarMPC_sub1 = planarMPC_s + 42;
planarMPC_FLOAT* planarMPC_lubbysub1 = planarMPC_lbys + 42;
planarMPC_FLOAT planarMPC_riub1[14];
planarMPC_FLOAT* planarMPC_dlubaff1 = planarMPC_dl_aff + 42;
planarMPC_FLOAT* planarMPC_dsubaff1 = planarMPC_ds_aff + 42;
planarMPC_FLOAT* planarMPC_dlubcc1 = planarMPC_dl_cc + 42;
planarMPC_FLOAT* planarMPC_dsubcc1 = planarMPC_ds_cc + 42;
planarMPC_FLOAT* planarMPC_ccrhsub1 = planarMPC_ccrhs + 42;
planarMPC_FLOAT planarMPC_Phi1[14];
planarMPC_FLOAT planarMPC_D1[14] = {-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000};
planarMPC_FLOAT planarMPC_W1[14];
planarMPC_FLOAT planarMPC_Ysd1[49];
planarMPC_FLOAT planarMPC_Lsd1[49];
planarMPC_FLOAT* planarMPC_z2 = planarMPC_z + 28;
planarMPC_FLOAT* planarMPC_dzaff2 = planarMPC_dz_aff + 28;
planarMPC_FLOAT* planarMPC_dzcc2 = planarMPC_dz_cc + 28;
planarMPC_FLOAT* planarMPC_rd2 = planarMPC_rd + 28;
planarMPC_FLOAT planarMPC_Lbyrd2[14];
planarMPC_FLOAT* planarMPC_grad_cost2 = planarMPC_grad_cost + 28;
planarMPC_FLOAT* planarMPC_grad_eq2 = planarMPC_grad_eq + 28;
planarMPC_FLOAT* planarMPC_grad_ineq2 = planarMPC_grad_ineq + 28;
planarMPC_FLOAT planarMPC_ctv2[14];
planarMPC_FLOAT* planarMPC_v2 = planarMPC_v + 14;
planarMPC_FLOAT planarMPC_re2[7];
planarMPC_FLOAT planarMPC_beta2[7];
planarMPC_FLOAT planarMPC_betacc2[7];
planarMPC_FLOAT* planarMPC_dvaff2 = planarMPC_dv_aff + 14;
planarMPC_FLOAT* planarMPC_dvcc2 = planarMPC_dv_cc + 14;
planarMPC_FLOAT planarMPC_V2[98];
planarMPC_FLOAT planarMPC_Yd2[28];
planarMPC_FLOAT planarMPC_Ld2[28];
planarMPC_FLOAT planarMPC_yy2[7];
planarMPC_FLOAT planarMPC_bmy2[7];
planarMPC_FLOAT planarMPC_c2[7] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
int planarMPC_lbIdx2[14] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13};
planarMPC_FLOAT* planarMPC_llb2 = planarMPC_l + 56;
planarMPC_FLOAT* planarMPC_slb2 = planarMPC_s + 56;
planarMPC_FLOAT* planarMPC_llbbyslb2 = planarMPC_lbys + 56;
planarMPC_FLOAT planarMPC_rilb2[14];
planarMPC_FLOAT* planarMPC_dllbaff2 = planarMPC_dl_aff + 56;
planarMPC_FLOAT* planarMPC_dslbaff2 = planarMPC_ds_aff + 56;
planarMPC_FLOAT* planarMPC_dllbcc2 = planarMPC_dl_cc + 56;
planarMPC_FLOAT* planarMPC_dslbcc2 = planarMPC_ds_cc + 56;
planarMPC_FLOAT* planarMPC_ccrhsl2 = planarMPC_ccrhs + 56;
int planarMPC_ubIdx2[14] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13};
planarMPC_FLOAT* planarMPC_lub2 = planarMPC_l + 70;
planarMPC_FLOAT* planarMPC_sub2 = planarMPC_s + 70;
planarMPC_FLOAT* planarMPC_lubbysub2 = planarMPC_lbys + 70;
planarMPC_FLOAT planarMPC_riub2[14];
planarMPC_FLOAT* planarMPC_dlubaff2 = planarMPC_dl_aff + 70;
planarMPC_FLOAT* planarMPC_dsubaff2 = planarMPC_ds_aff + 70;
planarMPC_FLOAT* planarMPC_dlubcc2 = planarMPC_dl_cc + 70;
planarMPC_FLOAT* planarMPC_dsubcc2 = planarMPC_ds_cc + 70;
planarMPC_FLOAT* planarMPC_ccrhsub2 = planarMPC_ccrhs + 70;
planarMPC_FLOAT planarMPC_Phi2[14];
planarMPC_FLOAT planarMPC_W2[14];
planarMPC_FLOAT planarMPC_Ysd2[49];
planarMPC_FLOAT planarMPC_Lsd2[49];
planarMPC_FLOAT* planarMPC_z3 = planarMPC_z + 42;
planarMPC_FLOAT* planarMPC_dzaff3 = planarMPC_dz_aff + 42;
planarMPC_FLOAT* planarMPC_dzcc3 = planarMPC_dz_cc + 42;
planarMPC_FLOAT* planarMPC_rd3 = planarMPC_rd + 42;
planarMPC_FLOAT planarMPC_Lbyrd3[14];
planarMPC_FLOAT* planarMPC_grad_cost3 = planarMPC_grad_cost + 42;
planarMPC_FLOAT* planarMPC_grad_eq3 = planarMPC_grad_eq + 42;
planarMPC_FLOAT* planarMPC_grad_ineq3 = planarMPC_grad_ineq + 42;
planarMPC_FLOAT planarMPC_ctv3[14];
planarMPC_FLOAT* planarMPC_v3 = planarMPC_v + 21;
planarMPC_FLOAT planarMPC_re3[7];
planarMPC_FLOAT planarMPC_beta3[7];
planarMPC_FLOAT planarMPC_betacc3[7];
planarMPC_FLOAT* planarMPC_dvaff3 = planarMPC_dv_aff + 21;
planarMPC_FLOAT* planarMPC_dvcc3 = planarMPC_dv_cc + 21;
planarMPC_FLOAT planarMPC_V3[98];
planarMPC_FLOAT planarMPC_Yd3[28];
planarMPC_FLOAT planarMPC_Ld3[28];
planarMPC_FLOAT planarMPC_yy3[7];
planarMPC_FLOAT planarMPC_bmy3[7];
planarMPC_FLOAT planarMPC_c3[7] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
int planarMPC_lbIdx3[14] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13};
planarMPC_FLOAT* planarMPC_llb3 = planarMPC_l + 84;
planarMPC_FLOAT* planarMPC_slb3 = planarMPC_s + 84;
planarMPC_FLOAT* planarMPC_llbbyslb3 = planarMPC_lbys + 84;
planarMPC_FLOAT planarMPC_rilb3[14];
planarMPC_FLOAT* planarMPC_dllbaff3 = planarMPC_dl_aff + 84;
planarMPC_FLOAT* planarMPC_dslbaff3 = planarMPC_ds_aff + 84;
planarMPC_FLOAT* planarMPC_dllbcc3 = planarMPC_dl_cc + 84;
planarMPC_FLOAT* planarMPC_dslbcc3 = planarMPC_ds_cc + 84;
planarMPC_FLOAT* planarMPC_ccrhsl3 = planarMPC_ccrhs + 84;
int planarMPC_ubIdx3[14] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13};
planarMPC_FLOAT* planarMPC_lub3 = planarMPC_l + 98;
planarMPC_FLOAT* planarMPC_sub3 = planarMPC_s + 98;
planarMPC_FLOAT* planarMPC_lubbysub3 = planarMPC_lbys + 98;
planarMPC_FLOAT planarMPC_riub3[14];
planarMPC_FLOAT* planarMPC_dlubaff3 = planarMPC_dl_aff + 98;
planarMPC_FLOAT* planarMPC_dsubaff3 = planarMPC_ds_aff + 98;
planarMPC_FLOAT* planarMPC_dlubcc3 = planarMPC_dl_cc + 98;
planarMPC_FLOAT* planarMPC_dsubcc3 = planarMPC_ds_cc + 98;
planarMPC_FLOAT* planarMPC_ccrhsub3 = planarMPC_ccrhs + 98;
planarMPC_FLOAT planarMPC_Phi3[14];
planarMPC_FLOAT planarMPC_W3[14];
planarMPC_FLOAT planarMPC_Ysd3[49];
planarMPC_FLOAT planarMPC_Lsd3[49];
planarMPC_FLOAT* planarMPC_z4 = planarMPC_z + 56;
planarMPC_FLOAT* planarMPC_dzaff4 = planarMPC_dz_aff + 56;
planarMPC_FLOAT* planarMPC_dzcc4 = planarMPC_dz_cc + 56;
planarMPC_FLOAT* planarMPC_rd4 = planarMPC_rd + 56;
planarMPC_FLOAT planarMPC_Lbyrd4[7];
planarMPC_FLOAT* planarMPC_grad_cost4 = planarMPC_grad_cost + 56;
planarMPC_FLOAT* planarMPC_grad_eq4 = planarMPC_grad_eq + 56;
planarMPC_FLOAT* planarMPC_grad_ineq4 = planarMPC_grad_ineq + 56;
planarMPC_FLOAT planarMPC_ctv4[7];
planarMPC_FLOAT* planarMPC_v4 = planarMPC_v + 28;
planarMPC_FLOAT planarMPC_re4[7];
planarMPC_FLOAT planarMPC_beta4[7];
planarMPC_FLOAT planarMPC_betacc4[7];
planarMPC_FLOAT* planarMPC_dvaff4 = planarMPC_dv_aff + 28;
planarMPC_FLOAT* planarMPC_dvcc4 = planarMPC_dv_cc + 28;
planarMPC_FLOAT planarMPC_V4[49];
planarMPC_FLOAT planarMPC_Yd4[28];
planarMPC_FLOAT planarMPC_Ld4[28];
planarMPC_FLOAT planarMPC_yy4[7];
planarMPC_FLOAT planarMPC_bmy4[7];
planarMPC_FLOAT planarMPC_c4[7] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
int planarMPC_lbIdx4[7] = {0, 1, 2, 3, 4, 5, 6};
planarMPC_FLOAT* planarMPC_llb4 = planarMPC_l + 112;
planarMPC_FLOAT* planarMPC_slb4 = planarMPC_s + 112;
planarMPC_FLOAT* planarMPC_llbbyslb4 = planarMPC_lbys + 112;
planarMPC_FLOAT planarMPC_rilb4[7];
planarMPC_FLOAT* planarMPC_dllbaff4 = planarMPC_dl_aff + 112;
planarMPC_FLOAT* planarMPC_dslbaff4 = planarMPC_ds_aff + 112;
planarMPC_FLOAT* planarMPC_dllbcc4 = planarMPC_dl_cc + 112;
planarMPC_FLOAT* planarMPC_dslbcc4 = planarMPC_ds_cc + 112;
planarMPC_FLOAT* planarMPC_ccrhsl4 = planarMPC_ccrhs + 112;
int planarMPC_ubIdx4[7] = {0, 1, 2, 3, 4, 5, 6};
planarMPC_FLOAT* planarMPC_lub4 = planarMPC_l + 119;
planarMPC_FLOAT* planarMPC_sub4 = planarMPC_s + 119;
planarMPC_FLOAT* planarMPC_lubbysub4 = planarMPC_lbys + 119;
planarMPC_FLOAT planarMPC_riub4[7];
planarMPC_FLOAT* planarMPC_dlubaff4 = planarMPC_dl_aff + 119;
planarMPC_FLOAT* planarMPC_dsubaff4 = planarMPC_ds_aff + 119;
planarMPC_FLOAT* planarMPC_dlubcc4 = planarMPC_dl_cc + 119;
planarMPC_FLOAT* planarMPC_dsubcc4 = planarMPC_ds_cc + 119;
planarMPC_FLOAT* planarMPC_ccrhsub4 = planarMPC_ccrhs + 119;
planarMPC_FLOAT planarMPC_Phi4[7];
planarMPC_FLOAT planarMPC_D4[7] = {-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000};
planarMPC_FLOAT planarMPC_W4[7];
planarMPC_FLOAT planarMPC_Ysd4[49];
planarMPC_FLOAT planarMPC_Lsd4[49];
planarMPC_FLOAT musigma;
planarMPC_FLOAT sigma_3rdroot;
planarMPC_FLOAT planarMPC_Diag1_0[14];
planarMPC_FLOAT planarMPC_Diag2_0[14];
planarMPC_FLOAT planarMPC_L_0[91];




/* SOLVER CODE --------------------------------------------------------- */
int planarMPC_solve(planarMPC_params* params, planarMPC_output* output, planarMPC_info* info)
{	
int exitcode;

#if planarMPC_SET_TIMING == 1
	planarMPC_timer solvertimer;
	planarMPC_tic(&solvertimer);
#endif
/* FUNCTION CALLS INTO LA LIBRARY -------------------------------------- */
info->it = 0;
planarMPC_LA_INITIALIZEVECTOR_63(planarMPC_z, 0);
planarMPC_LA_INITIALIZEVECTOR_35(planarMPC_v, 1);
planarMPC_LA_INITIALIZEVECTOR_126(planarMPC_l, 10);
planarMPC_LA_INITIALIZEVECTOR_126(planarMPC_s, 10);
info->mu = 0;
planarMPC_LA_DOTACC_126(planarMPC_l, planarMPC_s, &info->mu);
info->mu /= 126;
while( 1 ){
info->pobj = 0;
planarMPC_LA_DIAG_QUADFCN_14(params->H1, params->f1, planarMPC_z0, planarMPC_grad_cost0, &info->pobj);
planarMPC_LA_DIAG_QUADFCN_14(params->H2, params->f2, planarMPC_z1, planarMPC_grad_cost1, &info->pobj);
planarMPC_LA_DIAG_QUADFCN_14(params->H3, params->f3, planarMPC_z2, planarMPC_grad_cost2, &info->pobj);
planarMPC_LA_DIAG_QUADFCN_14(params->H4, params->f4, planarMPC_z3, planarMPC_grad_cost3, &info->pobj);
planarMPC_LA_DIAG_QUADFCN_7(params->H5, params->f5, planarMPC_z4, planarMPC_grad_cost4, &info->pobj);
info->res_eq = 0;
info->dgap = 0;
planarMPC_LA_DIAGZERO_MVMSUB6_7(planarMPC_D0, planarMPC_z0, params->c1, planarMPC_v0, planarMPC_re0, &info->dgap, &info->res_eq);
planarMPC_LA_DENSE_DIAGZERO_MVMSUB3_7_14_14(planarMPC_C0, planarMPC_z0, planarMPC_D1, planarMPC_z1, planarMPC_c1, planarMPC_v1, planarMPC_re1, &info->dgap, &info->res_eq);
planarMPC_LA_DENSE_DIAGZERO_MVMSUB3_7_14_14(planarMPC_C0, planarMPC_z1, planarMPC_D1, planarMPC_z2, planarMPC_c2, planarMPC_v2, planarMPC_re2, &info->dgap, &info->res_eq);
planarMPC_LA_DENSE_DIAGZERO_MVMSUB3_7_14_14(planarMPC_C0, planarMPC_z2, planarMPC_D1, planarMPC_z3, planarMPC_c3, planarMPC_v3, planarMPC_re3, &info->dgap, &info->res_eq);
planarMPC_LA_DENSE_DIAGZERO_MVMSUB3_7_14_7(planarMPC_C0, planarMPC_z3, planarMPC_D4, planarMPC_z4, planarMPC_c4, planarMPC_v4, planarMPC_re4, &info->dgap, &info->res_eq);
planarMPC_LA_DENSE_DIAGZERO_MTVM2_7_14_7(planarMPC_C0, planarMPC_v1, planarMPC_D0, planarMPC_v0, planarMPC_grad_eq0);
planarMPC_LA_DENSE_DIAGZERO_MTVM2_7_14_7(planarMPC_C0, planarMPC_v2, planarMPC_D1, planarMPC_v1, planarMPC_grad_eq1);
planarMPC_LA_DENSE_DIAGZERO_MTVM2_7_14_7(planarMPC_C0, planarMPC_v3, planarMPC_D1, planarMPC_v2, planarMPC_grad_eq2);
planarMPC_LA_DENSE_DIAGZERO_MTVM2_7_14_7(planarMPC_C0, planarMPC_v4, planarMPC_D1, planarMPC_v3, planarMPC_grad_eq3);
planarMPC_LA_DIAGZERO_MTVM_7_7(planarMPC_D4, planarMPC_v4, planarMPC_grad_eq4);
info->res_ineq = 0;
planarMPC_LA_VSUBADD3_14(params->lb1, planarMPC_z0, planarMPC_lbIdx0, planarMPC_llb0, planarMPC_slb0, planarMPC_rilb0, &info->dgap, &info->res_ineq);
planarMPC_LA_VSUBADD2_14(planarMPC_z0, planarMPC_ubIdx0, params->ub1, planarMPC_lub0, planarMPC_sub0, planarMPC_riub0, &info->dgap, &info->res_ineq);
planarMPC_LA_VSUBADD3_14(params->lb2, planarMPC_z1, planarMPC_lbIdx1, planarMPC_llb1, planarMPC_slb1, planarMPC_rilb1, &info->dgap, &info->res_ineq);
planarMPC_LA_VSUBADD2_14(planarMPC_z1, planarMPC_ubIdx1, params->ub2, planarMPC_lub1, planarMPC_sub1, planarMPC_riub1, &info->dgap, &info->res_ineq);
planarMPC_LA_VSUBADD3_14(params->lb3, planarMPC_z2, planarMPC_lbIdx2, planarMPC_llb2, planarMPC_slb2, planarMPC_rilb2, &info->dgap, &info->res_ineq);
planarMPC_LA_VSUBADD2_14(planarMPC_z2, planarMPC_ubIdx2, params->ub3, planarMPC_lub2, planarMPC_sub2, planarMPC_riub2, &info->dgap, &info->res_ineq);
planarMPC_LA_VSUBADD3_14(params->lb4, planarMPC_z3, planarMPC_lbIdx3, planarMPC_llb3, planarMPC_slb3, planarMPC_rilb3, &info->dgap, &info->res_ineq);
planarMPC_LA_VSUBADD2_14(planarMPC_z3, planarMPC_ubIdx3, params->ub4, planarMPC_lub3, planarMPC_sub3, planarMPC_riub3, &info->dgap, &info->res_ineq);
planarMPC_LA_VSUBADD3_7(params->lb5, planarMPC_z4, planarMPC_lbIdx4, planarMPC_llb4, planarMPC_slb4, planarMPC_rilb4, &info->dgap, &info->res_ineq);
planarMPC_LA_VSUBADD2_7(planarMPC_z4, planarMPC_ubIdx4, params->ub5, planarMPC_lub4, planarMPC_sub4, planarMPC_riub4, &info->dgap, &info->res_ineq);
planarMPC_LA_INEQ_B_GRAD_14_14_14(planarMPC_lub0, planarMPC_sub0, planarMPC_riub0, planarMPC_llb0, planarMPC_slb0, planarMPC_rilb0, planarMPC_lbIdx0, planarMPC_ubIdx0, planarMPC_grad_ineq0, planarMPC_lubbysub0, planarMPC_llbbyslb0);
planarMPC_LA_INEQ_B_GRAD_14_14_14(planarMPC_lub1, planarMPC_sub1, planarMPC_riub1, planarMPC_llb1, planarMPC_slb1, planarMPC_rilb1, planarMPC_lbIdx1, planarMPC_ubIdx1, planarMPC_grad_ineq1, planarMPC_lubbysub1, planarMPC_llbbyslb1);
planarMPC_LA_INEQ_B_GRAD_14_14_14(planarMPC_lub2, planarMPC_sub2, planarMPC_riub2, planarMPC_llb2, planarMPC_slb2, planarMPC_rilb2, planarMPC_lbIdx2, planarMPC_ubIdx2, planarMPC_grad_ineq2, planarMPC_lubbysub2, planarMPC_llbbyslb2);
planarMPC_LA_INEQ_B_GRAD_14_14_14(planarMPC_lub3, planarMPC_sub3, planarMPC_riub3, planarMPC_llb3, planarMPC_slb3, planarMPC_rilb3, planarMPC_lbIdx3, planarMPC_ubIdx3, planarMPC_grad_ineq3, planarMPC_lubbysub3, planarMPC_llbbyslb3);
planarMPC_LA_INEQ_B_GRAD_7_7_7(planarMPC_lub4, planarMPC_sub4, planarMPC_riub4, planarMPC_llb4, planarMPC_slb4, planarMPC_rilb4, planarMPC_lbIdx4, planarMPC_ubIdx4, planarMPC_grad_ineq4, planarMPC_lubbysub4, planarMPC_llbbyslb4);
info->dobj = info->pobj - info->dgap;
info->rdgap = info->pobj ? info->dgap / info->pobj : 1e6;
if( info->rdgap < 0 ) info->rdgap = -info->rdgap;
if( info->mu < planarMPC_SET_ACC_KKTCOMPL
    && (info->rdgap < planarMPC_SET_ACC_RDGAP || info->dgap < planarMPC_SET_ACC_KKTCOMPL)
    && info->res_eq < planarMPC_SET_ACC_RESEQ
    && info->res_ineq < planarMPC_SET_ACC_RESINEQ ){
exitcode = planarMPC_OPTIMAL; break; }
if( info->it == planarMPC_SET_MAXIT ){
exitcode = planarMPC_MAXITREACHED; break; }
planarMPC_LA_VVADD3_63(planarMPC_grad_cost, planarMPC_grad_eq, planarMPC_grad_ineq, planarMPC_rd);
planarMPC_LA_DIAG_CHOL_ONELOOP_LBUB_14_14_14(params->H1, planarMPC_llbbyslb0, planarMPC_lbIdx0, planarMPC_lubbysub0, planarMPC_ubIdx0, planarMPC_Phi0);
planarMPC_LA_DIAG_MATRIXFORWARDSUB_7_14(planarMPC_Phi0, planarMPC_C0, planarMPC_V0);
planarMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_7_14(planarMPC_Phi0, planarMPC_D0, planarMPC_W0);
planarMPC_LA_DENSE_DIAGZERO_MMTM_7_14_7(planarMPC_W0, planarMPC_V0, planarMPC_Ysd1);
planarMPC_LA_DIAG_FORWARDSUB_14(planarMPC_Phi0, planarMPC_rd0, planarMPC_Lbyrd0);
planarMPC_LA_DIAG_CHOL_ONELOOP_LBUB_14_14_14(params->H2, planarMPC_llbbyslb1, planarMPC_lbIdx1, planarMPC_lubbysub1, planarMPC_ubIdx1, planarMPC_Phi1);
planarMPC_LA_DIAG_MATRIXFORWARDSUB_7_14(planarMPC_Phi1, planarMPC_C0, planarMPC_V1);
planarMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_7_14(planarMPC_Phi1, planarMPC_D1, planarMPC_W1);
planarMPC_LA_DENSE_DIAGZERO_MMTM_7_14_7(planarMPC_W1, planarMPC_V1, planarMPC_Ysd2);
planarMPC_LA_DIAG_FORWARDSUB_14(planarMPC_Phi1, planarMPC_rd1, planarMPC_Lbyrd1);
planarMPC_LA_DIAG_CHOL_ONELOOP_LBUB_14_14_14(params->H3, planarMPC_llbbyslb2, planarMPC_lbIdx2, planarMPC_lubbysub2, planarMPC_ubIdx2, planarMPC_Phi2);
planarMPC_LA_DIAG_MATRIXFORWARDSUB_7_14(planarMPC_Phi2, planarMPC_C0, planarMPC_V2);
planarMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_7_14(planarMPC_Phi2, planarMPC_D1, planarMPC_W2);
planarMPC_LA_DENSE_DIAGZERO_MMTM_7_14_7(planarMPC_W2, planarMPC_V2, planarMPC_Ysd3);
planarMPC_LA_DIAG_FORWARDSUB_14(planarMPC_Phi2, planarMPC_rd2, planarMPC_Lbyrd2);
planarMPC_LA_DIAG_CHOL_ONELOOP_LBUB_14_14_14(params->H4, planarMPC_llbbyslb3, planarMPC_lbIdx3, planarMPC_lubbysub3, planarMPC_ubIdx3, planarMPC_Phi3);
planarMPC_LA_DIAG_MATRIXFORWARDSUB_7_14(planarMPC_Phi3, planarMPC_C0, planarMPC_V3);
planarMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_7_14(planarMPC_Phi3, planarMPC_D1, planarMPC_W3);
planarMPC_LA_DENSE_DIAGZERO_MMTM_7_14_7(planarMPC_W3, planarMPC_V3, planarMPC_Ysd4);
planarMPC_LA_DIAG_FORWARDSUB_14(planarMPC_Phi3, planarMPC_rd3, planarMPC_Lbyrd3);
planarMPC_LA_DIAG_CHOL_ONELOOP_LBUB_7_7_7(params->H5, planarMPC_llbbyslb4, planarMPC_lbIdx4, planarMPC_lubbysub4, planarMPC_ubIdx4, planarMPC_Phi4);
planarMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_7_7(planarMPC_Phi4, planarMPC_D4, planarMPC_W4);
planarMPC_LA_DIAG_FORWARDSUB_7(planarMPC_Phi4, planarMPC_rd4, planarMPC_Lbyrd4);
planarMPC_LA_DIAGZERO_MMT_7(planarMPC_W0, planarMPC_Yd0);
planarMPC_LA_DIAGZERO_MVMSUB7_7(planarMPC_W0, planarMPC_Lbyrd0, planarMPC_re0, planarMPC_beta0);
planarMPC_LA_DENSE_DIAGZERO_MMT2_7_14_14(planarMPC_V0, planarMPC_W1, planarMPC_Yd1);
planarMPC_LA_DENSE_DIAGZERO_2MVMSUB2_7_14_14(planarMPC_V0, planarMPC_Lbyrd0, planarMPC_W1, planarMPC_Lbyrd1, planarMPC_re1, planarMPC_beta1);
planarMPC_LA_DENSE_DIAGZERO_MMT2_7_14_14(planarMPC_V1, planarMPC_W2, planarMPC_Yd2);
planarMPC_LA_DENSE_DIAGZERO_2MVMSUB2_7_14_14(planarMPC_V1, planarMPC_Lbyrd1, planarMPC_W2, planarMPC_Lbyrd2, planarMPC_re2, planarMPC_beta2);
planarMPC_LA_DENSE_DIAGZERO_MMT2_7_14_14(planarMPC_V2, planarMPC_W3, planarMPC_Yd3);
planarMPC_LA_DENSE_DIAGZERO_2MVMSUB2_7_14_14(planarMPC_V2, planarMPC_Lbyrd2, planarMPC_W3, planarMPC_Lbyrd3, planarMPC_re3, planarMPC_beta3);
planarMPC_LA_DENSE_DIAGZERO_MMT2_7_14_7(planarMPC_V3, planarMPC_W4, planarMPC_Yd4);
planarMPC_LA_DENSE_DIAGZERO_2MVMSUB2_7_14_7(planarMPC_V3, planarMPC_Lbyrd3, planarMPC_W4, planarMPC_Lbyrd4, planarMPC_re4, planarMPC_beta4);
planarMPC_LA_DENSE_CHOL_7(planarMPC_Yd0, planarMPC_Ld0);
planarMPC_LA_DENSE_FORWARDSUB_7(planarMPC_Ld0, planarMPC_beta0, planarMPC_yy0);
planarMPC_LA_DENSE_MATRIXTFORWARDSUB_7_7(planarMPC_Ld0, planarMPC_Ysd1, planarMPC_Lsd1);
planarMPC_LA_DENSE_MMTSUB_7_7(planarMPC_Lsd1, planarMPC_Yd1);
planarMPC_LA_DENSE_CHOL_7(planarMPC_Yd1, planarMPC_Ld1);
planarMPC_LA_DENSE_MVMSUB1_7_7(planarMPC_Lsd1, planarMPC_yy0, planarMPC_beta1, planarMPC_bmy1);
planarMPC_LA_DENSE_FORWARDSUB_7(planarMPC_Ld1, planarMPC_bmy1, planarMPC_yy1);
planarMPC_LA_DENSE_MATRIXTFORWARDSUB_7_7(planarMPC_Ld1, planarMPC_Ysd2, planarMPC_Lsd2);
planarMPC_LA_DENSE_MMTSUB_7_7(planarMPC_Lsd2, planarMPC_Yd2);
planarMPC_LA_DENSE_CHOL_7(planarMPC_Yd2, planarMPC_Ld2);
planarMPC_LA_DENSE_MVMSUB1_7_7(planarMPC_Lsd2, planarMPC_yy1, planarMPC_beta2, planarMPC_bmy2);
planarMPC_LA_DENSE_FORWARDSUB_7(planarMPC_Ld2, planarMPC_bmy2, planarMPC_yy2);
planarMPC_LA_DENSE_MATRIXTFORWARDSUB_7_7(planarMPC_Ld2, planarMPC_Ysd3, planarMPC_Lsd3);
planarMPC_LA_DENSE_MMTSUB_7_7(planarMPC_Lsd3, planarMPC_Yd3);
planarMPC_LA_DENSE_CHOL_7(planarMPC_Yd3, planarMPC_Ld3);
planarMPC_LA_DENSE_MVMSUB1_7_7(planarMPC_Lsd3, planarMPC_yy2, planarMPC_beta3, planarMPC_bmy3);
planarMPC_LA_DENSE_FORWARDSUB_7(planarMPC_Ld3, planarMPC_bmy3, planarMPC_yy3);
planarMPC_LA_DENSE_MATRIXTFORWARDSUB_7_7(planarMPC_Ld3, planarMPC_Ysd4, planarMPC_Lsd4);
planarMPC_LA_DENSE_MMTSUB_7_7(planarMPC_Lsd4, planarMPC_Yd4);
planarMPC_LA_DENSE_CHOL_7(planarMPC_Yd4, planarMPC_Ld4);
planarMPC_LA_DENSE_MVMSUB1_7_7(planarMPC_Lsd4, planarMPC_yy3, planarMPC_beta4, planarMPC_bmy4);
planarMPC_LA_DENSE_FORWARDSUB_7(planarMPC_Ld4, planarMPC_bmy4, planarMPC_yy4);
planarMPC_LA_DENSE_BACKWARDSUB_7(planarMPC_Ld4, planarMPC_yy4, planarMPC_dvaff4);
planarMPC_LA_DENSE_MTVMSUB_7_7(planarMPC_Lsd4, planarMPC_dvaff4, planarMPC_yy3, planarMPC_bmy3);
planarMPC_LA_DENSE_BACKWARDSUB_7(planarMPC_Ld3, planarMPC_bmy3, planarMPC_dvaff3);
planarMPC_LA_DENSE_MTVMSUB_7_7(planarMPC_Lsd3, planarMPC_dvaff3, planarMPC_yy2, planarMPC_bmy2);
planarMPC_LA_DENSE_BACKWARDSUB_7(planarMPC_Ld2, planarMPC_bmy2, planarMPC_dvaff2);
planarMPC_LA_DENSE_MTVMSUB_7_7(planarMPC_Lsd2, planarMPC_dvaff2, planarMPC_yy1, planarMPC_bmy1);
planarMPC_LA_DENSE_BACKWARDSUB_7(planarMPC_Ld1, planarMPC_bmy1, planarMPC_dvaff1);
planarMPC_LA_DENSE_MTVMSUB_7_7(planarMPC_Lsd1, planarMPC_dvaff1, planarMPC_yy0, planarMPC_bmy0);
planarMPC_LA_DENSE_BACKWARDSUB_7(planarMPC_Ld0, planarMPC_bmy0, planarMPC_dvaff0);
planarMPC_LA_DENSE_DIAGZERO_MTVM2_7_14_7(planarMPC_C0, planarMPC_dvaff1, planarMPC_D0, planarMPC_dvaff0, planarMPC_grad_eq0);
planarMPC_LA_DENSE_DIAGZERO_MTVM2_7_14_7(planarMPC_C0, planarMPC_dvaff2, planarMPC_D1, planarMPC_dvaff1, planarMPC_grad_eq1);
planarMPC_LA_DENSE_DIAGZERO_MTVM2_7_14_7(planarMPC_C0, planarMPC_dvaff3, planarMPC_D1, planarMPC_dvaff2, planarMPC_grad_eq2);
planarMPC_LA_DENSE_DIAGZERO_MTVM2_7_14_7(planarMPC_C0, planarMPC_dvaff4, planarMPC_D1, planarMPC_dvaff3, planarMPC_grad_eq3);
planarMPC_LA_DIAGZERO_MTVM_7_7(planarMPC_D4, planarMPC_dvaff4, planarMPC_grad_eq4);
planarMPC_LA_VSUB2_63(planarMPC_rd, planarMPC_grad_eq, planarMPC_rd);
planarMPC_LA_DIAG_FORWARDBACKWARDSUB_14(planarMPC_Phi0, planarMPC_rd0, planarMPC_dzaff0);
planarMPC_LA_DIAG_FORWARDBACKWARDSUB_14(planarMPC_Phi1, planarMPC_rd1, planarMPC_dzaff1);
planarMPC_LA_DIAG_FORWARDBACKWARDSUB_14(planarMPC_Phi2, planarMPC_rd2, planarMPC_dzaff2);
planarMPC_LA_DIAG_FORWARDBACKWARDSUB_14(planarMPC_Phi3, planarMPC_rd3, planarMPC_dzaff3);
planarMPC_LA_DIAG_FORWARDBACKWARDSUB_7(planarMPC_Phi4, planarMPC_rd4, planarMPC_dzaff4);
planarMPC_LA_VSUB_INDEXED_14(planarMPC_dzaff0, planarMPC_lbIdx0, planarMPC_rilb0, planarMPC_dslbaff0);
planarMPC_LA_VSUB3_14(planarMPC_llbbyslb0, planarMPC_dslbaff0, planarMPC_llb0, planarMPC_dllbaff0);
planarMPC_LA_VSUB2_INDEXED_14(planarMPC_riub0, planarMPC_dzaff0, planarMPC_ubIdx0, planarMPC_dsubaff0);
planarMPC_LA_VSUB3_14(planarMPC_lubbysub0, planarMPC_dsubaff0, planarMPC_lub0, planarMPC_dlubaff0);
planarMPC_LA_VSUB_INDEXED_14(planarMPC_dzaff1, planarMPC_lbIdx1, planarMPC_rilb1, planarMPC_dslbaff1);
planarMPC_LA_VSUB3_14(planarMPC_llbbyslb1, planarMPC_dslbaff1, planarMPC_llb1, planarMPC_dllbaff1);
planarMPC_LA_VSUB2_INDEXED_14(planarMPC_riub1, planarMPC_dzaff1, planarMPC_ubIdx1, planarMPC_dsubaff1);
planarMPC_LA_VSUB3_14(planarMPC_lubbysub1, planarMPC_dsubaff1, planarMPC_lub1, planarMPC_dlubaff1);
planarMPC_LA_VSUB_INDEXED_14(planarMPC_dzaff2, planarMPC_lbIdx2, planarMPC_rilb2, planarMPC_dslbaff2);
planarMPC_LA_VSUB3_14(planarMPC_llbbyslb2, planarMPC_dslbaff2, planarMPC_llb2, planarMPC_dllbaff2);
planarMPC_LA_VSUB2_INDEXED_14(planarMPC_riub2, planarMPC_dzaff2, planarMPC_ubIdx2, planarMPC_dsubaff2);
planarMPC_LA_VSUB3_14(planarMPC_lubbysub2, planarMPC_dsubaff2, planarMPC_lub2, planarMPC_dlubaff2);
planarMPC_LA_VSUB_INDEXED_14(planarMPC_dzaff3, planarMPC_lbIdx3, planarMPC_rilb3, planarMPC_dslbaff3);
planarMPC_LA_VSUB3_14(planarMPC_llbbyslb3, planarMPC_dslbaff3, planarMPC_llb3, planarMPC_dllbaff3);
planarMPC_LA_VSUB2_INDEXED_14(planarMPC_riub3, planarMPC_dzaff3, planarMPC_ubIdx3, planarMPC_dsubaff3);
planarMPC_LA_VSUB3_14(planarMPC_lubbysub3, planarMPC_dsubaff3, planarMPC_lub3, planarMPC_dlubaff3);
planarMPC_LA_VSUB_INDEXED_7(planarMPC_dzaff4, planarMPC_lbIdx4, planarMPC_rilb4, planarMPC_dslbaff4);
planarMPC_LA_VSUB3_7(planarMPC_llbbyslb4, planarMPC_dslbaff4, planarMPC_llb4, planarMPC_dllbaff4);
planarMPC_LA_VSUB2_INDEXED_7(planarMPC_riub4, planarMPC_dzaff4, planarMPC_ubIdx4, planarMPC_dsubaff4);
planarMPC_LA_VSUB3_7(planarMPC_lubbysub4, planarMPC_dsubaff4, planarMPC_lub4, planarMPC_dlubaff4);
info->lsit_aff = planarMPC_LINESEARCH_BACKTRACKING_AFFINE(planarMPC_l, planarMPC_s, planarMPC_dl_aff, planarMPC_ds_aff, &info->step_aff, &info->mu_aff);
if( info->lsit_aff == planarMPC_NOPROGRESS ){
exitcode = planarMPC_NOPROGRESS; break;
}
sigma_3rdroot = info->mu_aff / info->mu;
info->sigma = sigma_3rdroot*sigma_3rdroot*sigma_3rdroot;
musigma = info->mu * info->sigma;
planarMPC_LA_VSUB5_126(planarMPC_ds_aff, planarMPC_dl_aff, info->mu, info->sigma, planarMPC_ccrhs);
planarMPC_LA_VSUB6_INDEXED_14_14_14(planarMPC_ccrhsub0, planarMPC_sub0, planarMPC_ubIdx0, planarMPC_ccrhsl0, planarMPC_slb0, planarMPC_lbIdx0, planarMPC_rd0);
planarMPC_LA_VSUB6_INDEXED_14_14_14(planarMPC_ccrhsub1, planarMPC_sub1, planarMPC_ubIdx1, planarMPC_ccrhsl1, planarMPC_slb1, planarMPC_lbIdx1, planarMPC_rd1);
planarMPC_LA_DIAG_FORWARDSUB_14(planarMPC_Phi0, planarMPC_rd0, planarMPC_Lbyrd0);
planarMPC_LA_DIAG_FORWARDSUB_14(planarMPC_Phi1, planarMPC_rd1, planarMPC_Lbyrd1);
planarMPC_LA_DIAGZERO_MVM_7(planarMPC_W0, planarMPC_Lbyrd0, planarMPC_beta0);
planarMPC_LA_DENSE_FORWARDSUB_7(planarMPC_Ld0, planarMPC_beta0, planarMPC_yy0);
planarMPC_LA_DENSE_DIAGZERO_2MVMADD_7_14_14(planarMPC_V0, planarMPC_Lbyrd0, planarMPC_W1, planarMPC_Lbyrd1, planarMPC_beta1);
planarMPC_LA_DENSE_MVMSUB1_7_7(planarMPC_Lsd1, planarMPC_yy0, planarMPC_beta1, planarMPC_bmy1);
planarMPC_LA_DENSE_FORWARDSUB_7(planarMPC_Ld1, planarMPC_bmy1, planarMPC_yy1);
planarMPC_LA_VSUB6_INDEXED_14_14_14(planarMPC_ccrhsub2, planarMPC_sub2, planarMPC_ubIdx2, planarMPC_ccrhsl2, planarMPC_slb2, planarMPC_lbIdx2, planarMPC_rd2);
planarMPC_LA_DIAG_FORWARDSUB_14(planarMPC_Phi2, planarMPC_rd2, planarMPC_Lbyrd2);
planarMPC_LA_DENSE_DIAGZERO_2MVMADD_7_14_14(planarMPC_V1, planarMPC_Lbyrd1, planarMPC_W2, planarMPC_Lbyrd2, planarMPC_beta2);
planarMPC_LA_DENSE_MVMSUB1_7_7(planarMPC_Lsd2, planarMPC_yy1, planarMPC_beta2, planarMPC_bmy2);
planarMPC_LA_DENSE_FORWARDSUB_7(planarMPC_Ld2, planarMPC_bmy2, planarMPC_yy2);
planarMPC_LA_VSUB6_INDEXED_14_14_14(planarMPC_ccrhsub3, planarMPC_sub3, planarMPC_ubIdx3, planarMPC_ccrhsl3, planarMPC_slb3, planarMPC_lbIdx3, planarMPC_rd3);
planarMPC_LA_DIAG_FORWARDSUB_14(planarMPC_Phi3, planarMPC_rd3, planarMPC_Lbyrd3);
planarMPC_LA_DENSE_DIAGZERO_2MVMADD_7_14_14(planarMPC_V2, planarMPC_Lbyrd2, planarMPC_W3, planarMPC_Lbyrd3, planarMPC_beta3);
planarMPC_LA_DENSE_MVMSUB1_7_7(planarMPC_Lsd3, planarMPC_yy2, planarMPC_beta3, planarMPC_bmy3);
planarMPC_LA_DENSE_FORWARDSUB_7(planarMPC_Ld3, planarMPC_bmy3, planarMPC_yy3);
planarMPC_LA_VSUB6_INDEXED_7_7_7(planarMPC_ccrhsub4, planarMPC_sub4, planarMPC_ubIdx4, planarMPC_ccrhsl4, planarMPC_slb4, planarMPC_lbIdx4, planarMPC_rd4);
planarMPC_LA_DIAG_FORWARDSUB_7(planarMPC_Phi4, planarMPC_rd4, planarMPC_Lbyrd4);
planarMPC_LA_DENSE_DIAGZERO_2MVMADD_7_14_7(planarMPC_V3, planarMPC_Lbyrd3, planarMPC_W4, planarMPC_Lbyrd4, planarMPC_beta4);
planarMPC_LA_DENSE_MVMSUB1_7_7(planarMPC_Lsd4, planarMPC_yy3, planarMPC_beta4, planarMPC_bmy4);
planarMPC_LA_DENSE_FORWARDSUB_7(planarMPC_Ld4, planarMPC_bmy4, planarMPC_yy4);
planarMPC_LA_DENSE_BACKWARDSUB_7(planarMPC_Ld4, planarMPC_yy4, planarMPC_dvcc4);
planarMPC_LA_DENSE_MTVMSUB_7_7(planarMPC_Lsd4, planarMPC_dvcc4, planarMPC_yy3, planarMPC_bmy3);
planarMPC_LA_DENSE_BACKWARDSUB_7(planarMPC_Ld3, planarMPC_bmy3, planarMPC_dvcc3);
planarMPC_LA_DENSE_MTVMSUB_7_7(planarMPC_Lsd3, planarMPC_dvcc3, planarMPC_yy2, planarMPC_bmy2);
planarMPC_LA_DENSE_BACKWARDSUB_7(planarMPC_Ld2, planarMPC_bmy2, planarMPC_dvcc2);
planarMPC_LA_DENSE_MTVMSUB_7_7(planarMPC_Lsd2, planarMPC_dvcc2, planarMPC_yy1, planarMPC_bmy1);
planarMPC_LA_DENSE_BACKWARDSUB_7(planarMPC_Ld1, planarMPC_bmy1, planarMPC_dvcc1);
planarMPC_LA_DENSE_MTVMSUB_7_7(planarMPC_Lsd1, planarMPC_dvcc1, planarMPC_yy0, planarMPC_bmy0);
planarMPC_LA_DENSE_BACKWARDSUB_7(planarMPC_Ld0, planarMPC_bmy0, planarMPC_dvcc0);
planarMPC_LA_DENSE_DIAGZERO_MTVM2_7_14_7(planarMPC_C0, planarMPC_dvcc1, planarMPC_D0, planarMPC_dvcc0, planarMPC_grad_eq0);
planarMPC_LA_DENSE_DIAGZERO_MTVM2_7_14_7(planarMPC_C0, planarMPC_dvcc2, planarMPC_D1, planarMPC_dvcc1, planarMPC_grad_eq1);
planarMPC_LA_DENSE_DIAGZERO_MTVM2_7_14_7(planarMPC_C0, planarMPC_dvcc3, planarMPC_D1, planarMPC_dvcc2, planarMPC_grad_eq2);
planarMPC_LA_DENSE_DIAGZERO_MTVM2_7_14_7(planarMPC_C0, planarMPC_dvcc4, planarMPC_D1, planarMPC_dvcc3, planarMPC_grad_eq3);
planarMPC_LA_DIAGZERO_MTVM_7_7(planarMPC_D4, planarMPC_dvcc4, planarMPC_grad_eq4);
planarMPC_LA_VSUB_63(planarMPC_rd, planarMPC_grad_eq, planarMPC_rd);
planarMPC_LA_DIAG_FORWARDBACKWARDSUB_14(planarMPC_Phi0, planarMPC_rd0, planarMPC_dzcc0);
planarMPC_LA_DIAG_FORWARDBACKWARDSUB_14(planarMPC_Phi1, planarMPC_rd1, planarMPC_dzcc1);
planarMPC_LA_DIAG_FORWARDBACKWARDSUB_14(planarMPC_Phi2, planarMPC_rd2, planarMPC_dzcc2);
planarMPC_LA_DIAG_FORWARDBACKWARDSUB_14(planarMPC_Phi3, planarMPC_rd3, planarMPC_dzcc3);
planarMPC_LA_DIAG_FORWARDBACKWARDSUB_7(planarMPC_Phi4, planarMPC_rd4, planarMPC_dzcc4);
planarMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_14(planarMPC_ccrhsl0, planarMPC_slb0, planarMPC_llbbyslb0, planarMPC_dzcc0, planarMPC_lbIdx0, planarMPC_dllbcc0);
planarMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_14(planarMPC_ccrhsub0, planarMPC_sub0, planarMPC_lubbysub0, planarMPC_dzcc0, planarMPC_ubIdx0, planarMPC_dlubcc0);
planarMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_14(planarMPC_ccrhsl1, planarMPC_slb1, planarMPC_llbbyslb1, planarMPC_dzcc1, planarMPC_lbIdx1, planarMPC_dllbcc1);
planarMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_14(planarMPC_ccrhsub1, planarMPC_sub1, planarMPC_lubbysub1, planarMPC_dzcc1, planarMPC_ubIdx1, planarMPC_dlubcc1);
planarMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_14(planarMPC_ccrhsl2, planarMPC_slb2, planarMPC_llbbyslb2, planarMPC_dzcc2, planarMPC_lbIdx2, planarMPC_dllbcc2);
planarMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_14(planarMPC_ccrhsub2, planarMPC_sub2, planarMPC_lubbysub2, planarMPC_dzcc2, planarMPC_ubIdx2, planarMPC_dlubcc2);
planarMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_14(planarMPC_ccrhsl3, planarMPC_slb3, planarMPC_llbbyslb3, planarMPC_dzcc3, planarMPC_lbIdx3, planarMPC_dllbcc3);
planarMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_14(planarMPC_ccrhsub3, planarMPC_sub3, planarMPC_lubbysub3, planarMPC_dzcc3, planarMPC_ubIdx3, planarMPC_dlubcc3);
planarMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_7(planarMPC_ccrhsl4, planarMPC_slb4, planarMPC_llbbyslb4, planarMPC_dzcc4, planarMPC_lbIdx4, planarMPC_dllbcc4);
planarMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_7(planarMPC_ccrhsub4, planarMPC_sub4, planarMPC_lubbysub4, planarMPC_dzcc4, planarMPC_ubIdx4, planarMPC_dlubcc4);
planarMPC_LA_VSUB7_126(planarMPC_l, planarMPC_ccrhs, planarMPC_s, planarMPC_dl_cc, planarMPC_ds_cc);
planarMPC_LA_VADD_63(planarMPC_dz_cc, planarMPC_dz_aff);
planarMPC_LA_VADD_35(planarMPC_dv_cc, planarMPC_dv_aff);
planarMPC_LA_VADD_126(planarMPC_dl_cc, planarMPC_dl_aff);
planarMPC_LA_VADD_126(planarMPC_ds_cc, planarMPC_ds_aff);
info->lsit_cc = planarMPC_LINESEARCH_BACKTRACKING_COMBINED(planarMPC_z, planarMPC_v, planarMPC_l, planarMPC_s, planarMPC_dz_cc, planarMPC_dv_cc, planarMPC_dl_cc, planarMPC_ds_cc, &info->step_cc, &info->mu);
if( info->lsit_cc == planarMPC_NOPROGRESS ){
exitcode = planarMPC_NOPROGRESS; break;
}
info->it++;
}
output->z1[0] = planarMPC_z0[0];
output->z1[1] = planarMPC_z0[1];
output->z1[2] = planarMPC_z0[2];
output->z1[3] = planarMPC_z0[3];
output->z1[4] = planarMPC_z0[4];
output->z1[5] = planarMPC_z0[5];
output->z1[6] = planarMPC_z0[6];
output->z1[7] = planarMPC_z0[7];
output->z1[8] = planarMPC_z0[8];
output->z1[9] = planarMPC_z0[9];
output->z1[10] = planarMPC_z0[10];
output->z1[11] = planarMPC_z0[11];
output->z1[12] = planarMPC_z0[12];
output->z1[13] = planarMPC_z0[13];
output->z2[0] = planarMPC_z1[0];
output->z2[1] = planarMPC_z1[1];
output->z2[2] = planarMPC_z1[2];
output->z2[3] = planarMPC_z1[3];
output->z2[4] = planarMPC_z1[4];
output->z2[5] = planarMPC_z1[5];
output->z2[6] = planarMPC_z1[6];
output->z2[7] = planarMPC_z1[7];
output->z2[8] = planarMPC_z1[8];
output->z2[9] = planarMPC_z1[9];
output->z2[10] = planarMPC_z1[10];
output->z2[11] = planarMPC_z1[11];
output->z2[12] = planarMPC_z1[12];
output->z2[13] = planarMPC_z1[13];
output->z3[0] = planarMPC_z2[0];
output->z3[1] = planarMPC_z2[1];
output->z3[2] = planarMPC_z2[2];
output->z3[3] = planarMPC_z2[3];
output->z3[4] = planarMPC_z2[4];
output->z3[5] = planarMPC_z2[5];
output->z3[6] = planarMPC_z2[6];
output->z3[7] = planarMPC_z2[7];
output->z3[8] = planarMPC_z2[8];
output->z3[9] = planarMPC_z2[9];
output->z3[10] = planarMPC_z2[10];
output->z3[11] = planarMPC_z2[11];
output->z3[12] = planarMPC_z2[12];
output->z3[13] = planarMPC_z2[13];
output->z4[0] = planarMPC_z3[0];
output->z4[1] = planarMPC_z3[1];
output->z4[2] = planarMPC_z3[2];
output->z4[3] = planarMPC_z3[3];
output->z4[4] = planarMPC_z3[4];
output->z4[5] = planarMPC_z3[5];
output->z4[6] = planarMPC_z3[6];
output->z4[7] = planarMPC_z3[7];
output->z4[8] = planarMPC_z3[8];
output->z4[9] = planarMPC_z3[9];
output->z4[10] = planarMPC_z3[10];
output->z4[11] = planarMPC_z3[11];
output->z4[12] = planarMPC_z3[12];
output->z4[13] = planarMPC_z3[13];
output->z5[0] = planarMPC_z4[0];
output->z5[1] = planarMPC_z4[1];
output->z5[2] = planarMPC_z4[2];
output->z5[3] = planarMPC_z4[3];
output->z5[4] = planarMPC_z4[4];
output->z5[5] = planarMPC_z4[5];
output->z5[6] = planarMPC_z4[6];

#if planarMPC_SET_TIMING == 1
info->solvetime = planarMPC_toc(&solvertimer);
#if planarMPC_SET_PRINTLEVEL > 0 && planarMPC_SET_TIMING == 1
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
