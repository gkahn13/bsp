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

#include "pr2MPC.h"

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
void pr2MPC_LA_INITIALIZEVECTOR_63(pr2MPC_FLOAT* vec, pr2MPC_FLOAT value)
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
void pr2MPC_LA_INITIALIZEVECTOR_35(pr2MPC_FLOAT* vec, pr2MPC_FLOAT value)
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
void pr2MPC_LA_INITIALIZEVECTOR_126(pr2MPC_FLOAT* vec, pr2MPC_FLOAT value)
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
void pr2MPC_LA_DOTACC_126(pr2MPC_FLOAT *x, pr2MPC_FLOAT *y, pr2MPC_FLOAT *z)
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
void pr2MPC_LA_DIAG_QUADFCN_14(pr2MPC_FLOAT* H, pr2MPC_FLOAT* f, pr2MPC_FLOAT* z, pr2MPC_FLOAT* grad, pr2MPC_FLOAT* value)
{
	int i;
	pr2MPC_FLOAT hz;	
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
void pr2MPC_LA_DIAG_QUADFCN_7(pr2MPC_FLOAT* H, pr2MPC_FLOAT* f, pr2MPC_FLOAT* z, pr2MPC_FLOAT* grad, pr2MPC_FLOAT* value)
{
	int i;
	pr2MPC_FLOAT hz;	
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
void pr2MPC_LA_DIAGZERO_MVMSUB6_7(pr2MPC_FLOAT *B, pr2MPC_FLOAT *u, pr2MPC_FLOAT *b, pr2MPC_FLOAT *l, pr2MPC_FLOAT *r, pr2MPC_FLOAT *z, pr2MPC_FLOAT *y)
{
	int i;
	pr2MPC_FLOAT Bu[7];
	pr2MPC_FLOAT norm = *y;
	pr2MPC_FLOAT lr = 0;

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
void pr2MPC_LA_DENSE_DIAGZERO_MVMSUB3_7_14_14(pr2MPC_FLOAT *A, pr2MPC_FLOAT *x, pr2MPC_FLOAT *B, pr2MPC_FLOAT *u, pr2MPC_FLOAT *b, pr2MPC_FLOAT *l, pr2MPC_FLOAT *r, pr2MPC_FLOAT *z, pr2MPC_FLOAT *y)
{
	int i;
	int j;
	int k = 0;
	pr2MPC_FLOAT AxBu[7];
	pr2MPC_FLOAT norm = *y;
	pr2MPC_FLOAT lr = 0;

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
void pr2MPC_LA_DENSE_DIAGZERO_MVMSUB3_7_14_7(pr2MPC_FLOAT *A, pr2MPC_FLOAT *x, pr2MPC_FLOAT *B, pr2MPC_FLOAT *u, pr2MPC_FLOAT *b, pr2MPC_FLOAT *l, pr2MPC_FLOAT *r, pr2MPC_FLOAT *z, pr2MPC_FLOAT *y)
{
	int i;
	int j;
	int k = 0;
	pr2MPC_FLOAT AxBu[7];
	pr2MPC_FLOAT norm = *y;
	pr2MPC_FLOAT lr = 0;

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
void pr2MPC_LA_DENSE_DIAGZERO_MTVM2_7_14_7(pr2MPC_FLOAT *A, pr2MPC_FLOAT *x, pr2MPC_FLOAT *B, pr2MPC_FLOAT *y, pr2MPC_FLOAT *z)
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
void pr2MPC_LA_DIAGZERO_MTVM_7_7(pr2MPC_FLOAT *M, pr2MPC_FLOAT *x, pr2MPC_FLOAT *y)
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
void pr2MPC_LA_VSUBADD3_14(pr2MPC_FLOAT* t, pr2MPC_FLOAT* u, int* uidx, pr2MPC_FLOAT* v, pr2MPC_FLOAT* w, pr2MPC_FLOAT* y, pr2MPC_FLOAT* z, pr2MPC_FLOAT* r)
{
	int i;
	pr2MPC_FLOAT norm = *r;
	pr2MPC_FLOAT vx = 0;
	pr2MPC_FLOAT x;
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
void pr2MPC_LA_VSUBADD2_14(pr2MPC_FLOAT* t, int* tidx, pr2MPC_FLOAT* u, pr2MPC_FLOAT* v, pr2MPC_FLOAT* w, pr2MPC_FLOAT* y, pr2MPC_FLOAT* z, pr2MPC_FLOAT* r)
{
	int i;
	pr2MPC_FLOAT norm = *r;
	pr2MPC_FLOAT vx = 0;
	pr2MPC_FLOAT x;
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
void pr2MPC_LA_VSUBADD3_7(pr2MPC_FLOAT* t, pr2MPC_FLOAT* u, int* uidx, pr2MPC_FLOAT* v, pr2MPC_FLOAT* w, pr2MPC_FLOAT* y, pr2MPC_FLOAT* z, pr2MPC_FLOAT* r)
{
	int i;
	pr2MPC_FLOAT norm = *r;
	pr2MPC_FLOAT vx = 0;
	pr2MPC_FLOAT x;
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
void pr2MPC_LA_VSUBADD2_7(pr2MPC_FLOAT* t, int* tidx, pr2MPC_FLOAT* u, pr2MPC_FLOAT* v, pr2MPC_FLOAT* w, pr2MPC_FLOAT* y, pr2MPC_FLOAT* z, pr2MPC_FLOAT* r)
{
	int i;
	pr2MPC_FLOAT norm = *r;
	pr2MPC_FLOAT vx = 0;
	pr2MPC_FLOAT x;
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
void pr2MPC_LA_INEQ_B_GRAD_14_14_14(pr2MPC_FLOAT *lu, pr2MPC_FLOAT *su, pr2MPC_FLOAT *ru, pr2MPC_FLOAT *ll, pr2MPC_FLOAT *sl, pr2MPC_FLOAT *rl, int* lbIdx, int* ubIdx, pr2MPC_FLOAT *grad, pr2MPC_FLOAT *lubysu, pr2MPC_FLOAT *llbysl)
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
void pr2MPC_LA_INEQ_B_GRAD_7_7_7(pr2MPC_FLOAT *lu, pr2MPC_FLOAT *su, pr2MPC_FLOAT *ru, pr2MPC_FLOAT *ll, pr2MPC_FLOAT *sl, pr2MPC_FLOAT *rl, int* lbIdx, int* ubIdx, pr2MPC_FLOAT *grad, pr2MPC_FLOAT *lubysu, pr2MPC_FLOAT *llbysl)
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
void pr2MPC_LA_VVADD3_63(pr2MPC_FLOAT *u, pr2MPC_FLOAT *v, pr2MPC_FLOAT *w, pr2MPC_FLOAT *z)
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
void pr2MPC_LA_DIAG_CHOL_ONELOOP_LBUB_14_14_14(pr2MPC_FLOAT *H, pr2MPC_FLOAT *llbysl, int* lbIdx, pr2MPC_FLOAT *lubysu, int* ubIdx, pr2MPC_FLOAT *Phi)


{
	int i;
	
	/* compute cholesky */
	for( i=0; i<14; i++ ){
		Phi[i] = H[i] + llbysl[i] + lubysu[i];

#if pr2MPC_SET_PRINTLEVEL > 0 && defined PRINTNUMERICALWARNINGS
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
void pr2MPC_LA_DIAG_MATRIXFORWARDSUB_7_14(pr2MPC_FLOAT *L, pr2MPC_FLOAT *B, pr2MPC_FLOAT *A)
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
void pr2MPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_7_14(pr2MPC_FLOAT *L, pr2MPC_FLOAT *B, pr2MPC_FLOAT *A)
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
void pr2MPC_LA_DENSE_DIAGZERO_MMTM_7_14_7(pr2MPC_FLOAT *A, pr2MPC_FLOAT *B, pr2MPC_FLOAT *C)
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
void pr2MPC_LA_DIAG_FORWARDSUB_14(pr2MPC_FLOAT *L, pr2MPC_FLOAT *b, pr2MPC_FLOAT *y)
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
void pr2MPC_LA_DIAG_CHOL_ONELOOP_LBUB_7_7_7(pr2MPC_FLOAT *H, pr2MPC_FLOAT *llbysl, int* lbIdx, pr2MPC_FLOAT *lubysu, int* ubIdx, pr2MPC_FLOAT *Phi)


{
	int i;
	
	/* compute cholesky */
	for( i=0; i<7; i++ ){
		Phi[i] = H[i] + llbysl[i] + lubysu[i];

#if pr2MPC_SET_PRINTLEVEL > 0 && defined PRINTNUMERICALWARNINGS
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
void pr2MPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_7_7(pr2MPC_FLOAT *L, pr2MPC_FLOAT *B, pr2MPC_FLOAT *A)
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
void pr2MPC_LA_DIAG_FORWARDSUB_7(pr2MPC_FLOAT *L, pr2MPC_FLOAT *b, pr2MPC_FLOAT *y)
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
void pr2MPC_LA_DIAGZERO_MMT_7(pr2MPC_FLOAT *B, pr2MPC_FLOAT *L)
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
void pr2MPC_LA_DIAGZERO_MVMSUB7_7(pr2MPC_FLOAT *B, pr2MPC_FLOAT *u, pr2MPC_FLOAT *b, pr2MPC_FLOAT *r)
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
void pr2MPC_LA_DENSE_DIAGZERO_MMT2_7_14_14(pr2MPC_FLOAT *A, pr2MPC_FLOAT *B, pr2MPC_FLOAT *L)
{
    int i, j, k, ii, di;
    pr2MPC_FLOAT ltemp;
    
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
void pr2MPC_LA_DENSE_DIAGZERO_2MVMSUB2_7_14_14(pr2MPC_FLOAT *A, pr2MPC_FLOAT *x, pr2MPC_FLOAT *B, pr2MPC_FLOAT *u, pr2MPC_FLOAT *b, pr2MPC_FLOAT *r)
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
void pr2MPC_LA_DENSE_DIAGZERO_MMT2_7_14_7(pr2MPC_FLOAT *A, pr2MPC_FLOAT *B, pr2MPC_FLOAT *L)
{
    int i, j, k, ii, di;
    pr2MPC_FLOAT ltemp;
    
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
void pr2MPC_LA_DENSE_DIAGZERO_2MVMSUB2_7_14_7(pr2MPC_FLOAT *A, pr2MPC_FLOAT *x, pr2MPC_FLOAT *B, pr2MPC_FLOAT *u, pr2MPC_FLOAT *b, pr2MPC_FLOAT *r)
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
void pr2MPC_LA_DENSE_CHOL_7(pr2MPC_FLOAT *A, pr2MPC_FLOAT *L)
{
    int i, j, k, di, dj;
	 int ii, jj;

    pr2MPC_FLOAT l;
    pr2MPC_FLOAT Mii;

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
        
#if pr2MPC_SET_PRINTLEVEL > 0 && defined PRINTNUMERICALWARNINGS
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
void pr2MPC_LA_DENSE_FORWARDSUB_7(pr2MPC_FLOAT *L, pr2MPC_FLOAT *b, pr2MPC_FLOAT *y)
{
    int i,j,ii,di;
    pr2MPC_FLOAT yel;
            
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
void pr2MPC_LA_DENSE_MATRIXTFORWARDSUB_7_7(pr2MPC_FLOAT *L, pr2MPC_FLOAT *B, pr2MPC_FLOAT *A)
{
    int i,j,k,ii,di;
    pr2MPC_FLOAT a;
    
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
void pr2MPC_LA_DENSE_MMTSUB_7_7(pr2MPC_FLOAT *A, pr2MPC_FLOAT *L)
{
    int i, j, k, ii, di;
    pr2MPC_FLOAT ltemp;
    
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
void pr2MPC_LA_DENSE_MVMSUB1_7_7(pr2MPC_FLOAT *A, pr2MPC_FLOAT *x, pr2MPC_FLOAT *b, pr2MPC_FLOAT *r)
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
void pr2MPC_LA_DENSE_BACKWARDSUB_7(pr2MPC_FLOAT *L, pr2MPC_FLOAT *y, pr2MPC_FLOAT *x)
{
    int i, ii, di, j, jj, dj;
    pr2MPC_FLOAT xel;    
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
void pr2MPC_LA_DENSE_MTVMSUB_7_7(pr2MPC_FLOAT *A, pr2MPC_FLOAT *x, pr2MPC_FLOAT *b, pr2MPC_FLOAT *r)
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
void pr2MPC_LA_VSUB2_63(pr2MPC_FLOAT *x, pr2MPC_FLOAT *y, pr2MPC_FLOAT *z)
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
void pr2MPC_LA_DIAG_FORWARDBACKWARDSUB_14(pr2MPC_FLOAT *L, pr2MPC_FLOAT *b, pr2MPC_FLOAT *x)
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
void pr2MPC_LA_DIAG_FORWARDBACKWARDSUB_7(pr2MPC_FLOAT *L, pr2MPC_FLOAT *b, pr2MPC_FLOAT *x)
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
void pr2MPC_LA_VSUB_INDEXED_14(pr2MPC_FLOAT *x, int* xidx, pr2MPC_FLOAT *y, pr2MPC_FLOAT *z)
{
	int i;
	for( i=0; i<14; i++){
		z[i] = x[xidx[i]] - y[i];
	}
}


/*
 * Vector subtraction x = -u.*v - w for vectors of length 14.
 */
void pr2MPC_LA_VSUB3_14(pr2MPC_FLOAT *u, pr2MPC_FLOAT *v, pr2MPC_FLOAT *w, pr2MPC_FLOAT *x)
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
void pr2MPC_LA_VSUB2_INDEXED_14(pr2MPC_FLOAT *x, pr2MPC_FLOAT *y, int* yidx, pr2MPC_FLOAT *z)
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
void pr2MPC_LA_VSUB_INDEXED_7(pr2MPC_FLOAT *x, int* xidx, pr2MPC_FLOAT *y, pr2MPC_FLOAT *z)
{
	int i;
	for( i=0; i<7; i++){
		z[i] = x[xidx[i]] - y[i];
	}
}


/*
 * Vector subtraction x = -u.*v - w for vectors of length 7.
 */
void pr2MPC_LA_VSUB3_7(pr2MPC_FLOAT *u, pr2MPC_FLOAT *v, pr2MPC_FLOAT *w, pr2MPC_FLOAT *x)
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
void pr2MPC_LA_VSUB2_INDEXED_7(pr2MPC_FLOAT *x, pr2MPC_FLOAT *y, int* yidx, pr2MPC_FLOAT *z)
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
 * pr2MPC_NOPROGRESS (should be negative).
 */
int pr2MPC_LINESEARCH_BACKTRACKING_AFFINE(pr2MPC_FLOAT *l, pr2MPC_FLOAT *s, pr2MPC_FLOAT *dl, pr2MPC_FLOAT *ds, pr2MPC_FLOAT *a, pr2MPC_FLOAT *mu_aff)
{
    int i;
	int lsIt=1;    
    pr2MPC_FLOAT dltemp;
    pr2MPC_FLOAT dstemp;
    pr2MPC_FLOAT mya = 1.0;
    pr2MPC_FLOAT mymu;
        
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
            mya *= pr2MPC_SET_LS_SCALE_AFF;
            if( mya < pr2MPC_SET_LS_MINSTEP ){
                return pr2MPC_NOPROGRESS;
            }
        }
    }
    
    /* return new values and iteration counter */
    *a = mya;
    *mu_aff = mymu / (pr2MPC_FLOAT)126;
    return lsIt;
}


/*
 * Vector subtraction x = (u.*v - mu)*sigma where a is a scalar
*  and x,u,v are vectors of length 126.
 */
void pr2MPC_LA_VSUB5_126(pr2MPC_FLOAT *u, pr2MPC_FLOAT *v, pr2MPC_FLOAT mu,  pr2MPC_FLOAT sigma, pr2MPC_FLOAT *x)
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
void pr2MPC_LA_VSUB6_INDEXED_14_14_14(pr2MPC_FLOAT *u, pr2MPC_FLOAT *su, int* uidx, pr2MPC_FLOAT *v, pr2MPC_FLOAT *sv, int* vidx, pr2MPC_FLOAT *x)
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
void pr2MPC_LA_DIAGZERO_MVM_7(pr2MPC_FLOAT *B, pr2MPC_FLOAT *u, pr2MPC_FLOAT *r)
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
void pr2MPC_LA_DENSE_DIAGZERO_2MVMADD_7_14_14(pr2MPC_FLOAT *A, pr2MPC_FLOAT *x, pr2MPC_FLOAT *B, pr2MPC_FLOAT *u, pr2MPC_FLOAT *r)
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
void pr2MPC_LA_VSUB6_INDEXED_7_7_7(pr2MPC_FLOAT *u, pr2MPC_FLOAT *su, int* uidx, pr2MPC_FLOAT *v, pr2MPC_FLOAT *sv, int* vidx, pr2MPC_FLOAT *x)
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
void pr2MPC_LA_DENSE_DIAGZERO_2MVMADD_7_14_7(pr2MPC_FLOAT *A, pr2MPC_FLOAT *x, pr2MPC_FLOAT *B, pr2MPC_FLOAT *u, pr2MPC_FLOAT *r)
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
void pr2MPC_LA_VSUB_63(pr2MPC_FLOAT *x, pr2MPC_FLOAT *y, pr2MPC_FLOAT *z)
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
void pr2MPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_14(pr2MPC_FLOAT *r, pr2MPC_FLOAT *s, pr2MPC_FLOAT *u, pr2MPC_FLOAT *y, int* yidx, pr2MPC_FLOAT *z)
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
void pr2MPC_LA_VEC_DIVSUB_MULTADD_INDEXED_14(pr2MPC_FLOAT *r, pr2MPC_FLOAT *s, pr2MPC_FLOAT *u, pr2MPC_FLOAT *y, int* yidx, pr2MPC_FLOAT *z)
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
void pr2MPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_7(pr2MPC_FLOAT *r, pr2MPC_FLOAT *s, pr2MPC_FLOAT *u, pr2MPC_FLOAT *y, int* yidx, pr2MPC_FLOAT *z)
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
void pr2MPC_LA_VEC_DIVSUB_MULTADD_INDEXED_7(pr2MPC_FLOAT *r, pr2MPC_FLOAT *s, pr2MPC_FLOAT *u, pr2MPC_FLOAT *y, int* yidx, pr2MPC_FLOAT *z)
{
	int i;
	for( i=0; i<7; i++ ){
		z[i] = -r[i]/s[i] + u[i]*y[yidx[i]];
	}
}


/*
 * Computes ds = -l.\(r + s.*dl) for vectors of length 126.
 */
void pr2MPC_LA_VSUB7_126(pr2MPC_FLOAT *l, pr2MPC_FLOAT *r, pr2MPC_FLOAT *s, pr2MPC_FLOAT *dl, pr2MPC_FLOAT *ds)
{
	int i;
	for( i=0; i<126; i++){
		ds[i] = -(r[i] + s[i]*dl[i])/l[i];
	}
}


/*
 * Vector addition x = x + y for vectors of length 63.
 */
void pr2MPC_LA_VADD_63(pr2MPC_FLOAT *x, pr2MPC_FLOAT *y)
{
	int i;
	for( i=0; i<63; i++){
		x[i] += y[i];
	}
}


/*
 * Vector addition x = x + y for vectors of length 35.
 */
void pr2MPC_LA_VADD_35(pr2MPC_FLOAT *x, pr2MPC_FLOAT *y)
{
	int i;
	for( i=0; i<35; i++){
		x[i] += y[i];
	}
}


/*
 * Vector addition x = x + y for vectors of length 126.
 */
void pr2MPC_LA_VADD_126(pr2MPC_FLOAT *x, pr2MPC_FLOAT *y)
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
int pr2MPC_LINESEARCH_BACKTRACKING_COMBINED(pr2MPC_FLOAT *z, pr2MPC_FLOAT *v, pr2MPC_FLOAT *l, pr2MPC_FLOAT *s, pr2MPC_FLOAT *dz, pr2MPC_FLOAT *dv, pr2MPC_FLOAT *dl, pr2MPC_FLOAT *ds, pr2MPC_FLOAT *a, pr2MPC_FLOAT *mu)
{
    int i, lsIt=1;       
    pr2MPC_FLOAT dltemp;
    pr2MPC_FLOAT dstemp;    
    pr2MPC_FLOAT a_gamma;
            
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
            *a *= pr2MPC_SET_LS_SCALE;
            if( *a < pr2MPC_SET_LS_MINSTEP ){
                return pr2MPC_NOPROGRESS;
            }
        }
    }
    
    /* update variables with safety margin */
    a_gamma = (*a)*pr2MPC_SET_LS_MAXSTEP;
    
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
    *mu /= (pr2MPC_FLOAT)126;
    return lsIt;
}




/* VARIABLE DEFINITIONS ------------------------------------------------ */
pr2MPC_FLOAT pr2MPC_z[63];
pr2MPC_FLOAT pr2MPC_v[35];
pr2MPC_FLOAT pr2MPC_dz_aff[63];
pr2MPC_FLOAT pr2MPC_dv_aff[35];
pr2MPC_FLOAT pr2MPC_grad_cost[63];
pr2MPC_FLOAT pr2MPC_grad_eq[63];
pr2MPC_FLOAT pr2MPC_rd[63];
pr2MPC_FLOAT pr2MPC_l[126];
pr2MPC_FLOAT pr2MPC_s[126];
pr2MPC_FLOAT pr2MPC_lbys[126];
pr2MPC_FLOAT pr2MPC_dl_aff[126];
pr2MPC_FLOAT pr2MPC_ds_aff[126];
pr2MPC_FLOAT pr2MPC_dz_cc[63];
pr2MPC_FLOAT pr2MPC_dv_cc[35];
pr2MPC_FLOAT pr2MPC_dl_cc[126];
pr2MPC_FLOAT pr2MPC_ds_cc[126];
pr2MPC_FLOAT pr2MPC_ccrhs[126];
pr2MPC_FLOAT pr2MPC_grad_ineq[63];
pr2MPC_FLOAT* pr2MPC_z0 = pr2MPC_z + 0;
pr2MPC_FLOAT* pr2MPC_dzaff0 = pr2MPC_dz_aff + 0;
pr2MPC_FLOAT* pr2MPC_dzcc0 = pr2MPC_dz_cc + 0;
pr2MPC_FLOAT* pr2MPC_rd0 = pr2MPC_rd + 0;
pr2MPC_FLOAT pr2MPC_Lbyrd0[14];
pr2MPC_FLOAT* pr2MPC_grad_cost0 = pr2MPC_grad_cost + 0;
pr2MPC_FLOAT* pr2MPC_grad_eq0 = pr2MPC_grad_eq + 0;
pr2MPC_FLOAT* pr2MPC_grad_ineq0 = pr2MPC_grad_ineq + 0;
pr2MPC_FLOAT pr2MPC_ctv0[14];
pr2MPC_FLOAT pr2MPC_C0[98] = {1.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 
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
pr2MPC_FLOAT* pr2MPC_v0 = pr2MPC_v + 0;
pr2MPC_FLOAT pr2MPC_re0[7];
pr2MPC_FLOAT pr2MPC_beta0[7];
pr2MPC_FLOAT pr2MPC_betacc0[7];
pr2MPC_FLOAT* pr2MPC_dvaff0 = pr2MPC_dv_aff + 0;
pr2MPC_FLOAT* pr2MPC_dvcc0 = pr2MPC_dv_cc + 0;
pr2MPC_FLOAT pr2MPC_V0[98];
pr2MPC_FLOAT pr2MPC_Yd0[28];
pr2MPC_FLOAT pr2MPC_Ld0[28];
pr2MPC_FLOAT pr2MPC_yy0[7];
pr2MPC_FLOAT pr2MPC_bmy0[7];
int pr2MPC_lbIdx0[14] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13};
pr2MPC_FLOAT* pr2MPC_llb0 = pr2MPC_l + 0;
pr2MPC_FLOAT* pr2MPC_slb0 = pr2MPC_s + 0;
pr2MPC_FLOAT* pr2MPC_llbbyslb0 = pr2MPC_lbys + 0;
pr2MPC_FLOAT pr2MPC_rilb0[14];
pr2MPC_FLOAT* pr2MPC_dllbaff0 = pr2MPC_dl_aff + 0;
pr2MPC_FLOAT* pr2MPC_dslbaff0 = pr2MPC_ds_aff + 0;
pr2MPC_FLOAT* pr2MPC_dllbcc0 = pr2MPC_dl_cc + 0;
pr2MPC_FLOAT* pr2MPC_dslbcc0 = pr2MPC_ds_cc + 0;
pr2MPC_FLOAT* pr2MPC_ccrhsl0 = pr2MPC_ccrhs + 0;
int pr2MPC_ubIdx0[14] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13};
pr2MPC_FLOAT* pr2MPC_lub0 = pr2MPC_l + 14;
pr2MPC_FLOAT* pr2MPC_sub0 = pr2MPC_s + 14;
pr2MPC_FLOAT* pr2MPC_lubbysub0 = pr2MPC_lbys + 14;
pr2MPC_FLOAT pr2MPC_riub0[14];
pr2MPC_FLOAT* pr2MPC_dlubaff0 = pr2MPC_dl_aff + 14;
pr2MPC_FLOAT* pr2MPC_dsubaff0 = pr2MPC_ds_aff + 14;
pr2MPC_FLOAT* pr2MPC_dlubcc0 = pr2MPC_dl_cc + 14;
pr2MPC_FLOAT* pr2MPC_dsubcc0 = pr2MPC_ds_cc + 14;
pr2MPC_FLOAT* pr2MPC_ccrhsub0 = pr2MPC_ccrhs + 14;
pr2MPC_FLOAT pr2MPC_Phi0[14];
pr2MPC_FLOAT pr2MPC_D0[14] = {1.0000000000000000E+000, 
1.0000000000000000E+000, 
1.0000000000000000E+000, 
1.0000000000000000E+000, 
1.0000000000000000E+000, 
1.0000000000000000E+000, 
1.0000000000000000E+000};
pr2MPC_FLOAT pr2MPC_W0[14];
pr2MPC_FLOAT* pr2MPC_z1 = pr2MPC_z + 14;
pr2MPC_FLOAT* pr2MPC_dzaff1 = pr2MPC_dz_aff + 14;
pr2MPC_FLOAT* pr2MPC_dzcc1 = pr2MPC_dz_cc + 14;
pr2MPC_FLOAT* pr2MPC_rd1 = pr2MPC_rd + 14;
pr2MPC_FLOAT pr2MPC_Lbyrd1[14];
pr2MPC_FLOAT* pr2MPC_grad_cost1 = pr2MPC_grad_cost + 14;
pr2MPC_FLOAT* pr2MPC_grad_eq1 = pr2MPC_grad_eq + 14;
pr2MPC_FLOAT* pr2MPC_grad_ineq1 = pr2MPC_grad_ineq + 14;
pr2MPC_FLOAT pr2MPC_ctv1[14];
pr2MPC_FLOAT* pr2MPC_v1 = pr2MPC_v + 7;
pr2MPC_FLOAT pr2MPC_re1[7];
pr2MPC_FLOAT pr2MPC_beta1[7];
pr2MPC_FLOAT pr2MPC_betacc1[7];
pr2MPC_FLOAT* pr2MPC_dvaff1 = pr2MPC_dv_aff + 7;
pr2MPC_FLOAT* pr2MPC_dvcc1 = pr2MPC_dv_cc + 7;
pr2MPC_FLOAT pr2MPC_V1[98];
pr2MPC_FLOAT pr2MPC_Yd1[28];
pr2MPC_FLOAT pr2MPC_Ld1[28];
pr2MPC_FLOAT pr2MPC_yy1[7];
pr2MPC_FLOAT pr2MPC_bmy1[7];
pr2MPC_FLOAT pr2MPC_c1[7] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
int pr2MPC_lbIdx1[14] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13};
pr2MPC_FLOAT* pr2MPC_llb1 = pr2MPC_l + 28;
pr2MPC_FLOAT* pr2MPC_slb1 = pr2MPC_s + 28;
pr2MPC_FLOAT* pr2MPC_llbbyslb1 = pr2MPC_lbys + 28;
pr2MPC_FLOAT pr2MPC_rilb1[14];
pr2MPC_FLOAT* pr2MPC_dllbaff1 = pr2MPC_dl_aff + 28;
pr2MPC_FLOAT* pr2MPC_dslbaff1 = pr2MPC_ds_aff + 28;
pr2MPC_FLOAT* pr2MPC_dllbcc1 = pr2MPC_dl_cc + 28;
pr2MPC_FLOAT* pr2MPC_dslbcc1 = pr2MPC_ds_cc + 28;
pr2MPC_FLOAT* pr2MPC_ccrhsl1 = pr2MPC_ccrhs + 28;
int pr2MPC_ubIdx1[14] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13};
pr2MPC_FLOAT* pr2MPC_lub1 = pr2MPC_l + 42;
pr2MPC_FLOAT* pr2MPC_sub1 = pr2MPC_s + 42;
pr2MPC_FLOAT* pr2MPC_lubbysub1 = pr2MPC_lbys + 42;
pr2MPC_FLOAT pr2MPC_riub1[14];
pr2MPC_FLOAT* pr2MPC_dlubaff1 = pr2MPC_dl_aff + 42;
pr2MPC_FLOAT* pr2MPC_dsubaff1 = pr2MPC_ds_aff + 42;
pr2MPC_FLOAT* pr2MPC_dlubcc1 = pr2MPC_dl_cc + 42;
pr2MPC_FLOAT* pr2MPC_dsubcc1 = pr2MPC_ds_cc + 42;
pr2MPC_FLOAT* pr2MPC_ccrhsub1 = pr2MPC_ccrhs + 42;
pr2MPC_FLOAT pr2MPC_Phi1[14];
pr2MPC_FLOAT pr2MPC_D1[14] = {-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000};
pr2MPC_FLOAT pr2MPC_W1[14];
pr2MPC_FLOAT pr2MPC_Ysd1[49];
pr2MPC_FLOAT pr2MPC_Lsd1[49];
pr2MPC_FLOAT* pr2MPC_z2 = pr2MPC_z + 28;
pr2MPC_FLOAT* pr2MPC_dzaff2 = pr2MPC_dz_aff + 28;
pr2MPC_FLOAT* pr2MPC_dzcc2 = pr2MPC_dz_cc + 28;
pr2MPC_FLOAT* pr2MPC_rd2 = pr2MPC_rd + 28;
pr2MPC_FLOAT pr2MPC_Lbyrd2[14];
pr2MPC_FLOAT* pr2MPC_grad_cost2 = pr2MPC_grad_cost + 28;
pr2MPC_FLOAT* pr2MPC_grad_eq2 = pr2MPC_grad_eq + 28;
pr2MPC_FLOAT* pr2MPC_grad_ineq2 = pr2MPC_grad_ineq + 28;
pr2MPC_FLOAT pr2MPC_ctv2[14];
pr2MPC_FLOAT* pr2MPC_v2 = pr2MPC_v + 14;
pr2MPC_FLOAT pr2MPC_re2[7];
pr2MPC_FLOAT pr2MPC_beta2[7];
pr2MPC_FLOAT pr2MPC_betacc2[7];
pr2MPC_FLOAT* pr2MPC_dvaff2 = pr2MPC_dv_aff + 14;
pr2MPC_FLOAT* pr2MPC_dvcc2 = pr2MPC_dv_cc + 14;
pr2MPC_FLOAT pr2MPC_V2[98];
pr2MPC_FLOAT pr2MPC_Yd2[28];
pr2MPC_FLOAT pr2MPC_Ld2[28];
pr2MPC_FLOAT pr2MPC_yy2[7];
pr2MPC_FLOAT pr2MPC_bmy2[7];
pr2MPC_FLOAT pr2MPC_c2[7] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
int pr2MPC_lbIdx2[14] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13};
pr2MPC_FLOAT* pr2MPC_llb2 = pr2MPC_l + 56;
pr2MPC_FLOAT* pr2MPC_slb2 = pr2MPC_s + 56;
pr2MPC_FLOAT* pr2MPC_llbbyslb2 = pr2MPC_lbys + 56;
pr2MPC_FLOAT pr2MPC_rilb2[14];
pr2MPC_FLOAT* pr2MPC_dllbaff2 = pr2MPC_dl_aff + 56;
pr2MPC_FLOAT* pr2MPC_dslbaff2 = pr2MPC_ds_aff + 56;
pr2MPC_FLOAT* pr2MPC_dllbcc2 = pr2MPC_dl_cc + 56;
pr2MPC_FLOAT* pr2MPC_dslbcc2 = pr2MPC_ds_cc + 56;
pr2MPC_FLOAT* pr2MPC_ccrhsl2 = pr2MPC_ccrhs + 56;
int pr2MPC_ubIdx2[14] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13};
pr2MPC_FLOAT* pr2MPC_lub2 = pr2MPC_l + 70;
pr2MPC_FLOAT* pr2MPC_sub2 = pr2MPC_s + 70;
pr2MPC_FLOAT* pr2MPC_lubbysub2 = pr2MPC_lbys + 70;
pr2MPC_FLOAT pr2MPC_riub2[14];
pr2MPC_FLOAT* pr2MPC_dlubaff2 = pr2MPC_dl_aff + 70;
pr2MPC_FLOAT* pr2MPC_dsubaff2 = pr2MPC_ds_aff + 70;
pr2MPC_FLOAT* pr2MPC_dlubcc2 = pr2MPC_dl_cc + 70;
pr2MPC_FLOAT* pr2MPC_dsubcc2 = pr2MPC_ds_cc + 70;
pr2MPC_FLOAT* pr2MPC_ccrhsub2 = pr2MPC_ccrhs + 70;
pr2MPC_FLOAT pr2MPC_Phi2[14];
pr2MPC_FLOAT pr2MPC_W2[14];
pr2MPC_FLOAT pr2MPC_Ysd2[49];
pr2MPC_FLOAT pr2MPC_Lsd2[49];
pr2MPC_FLOAT* pr2MPC_z3 = pr2MPC_z + 42;
pr2MPC_FLOAT* pr2MPC_dzaff3 = pr2MPC_dz_aff + 42;
pr2MPC_FLOAT* pr2MPC_dzcc3 = pr2MPC_dz_cc + 42;
pr2MPC_FLOAT* pr2MPC_rd3 = pr2MPC_rd + 42;
pr2MPC_FLOAT pr2MPC_Lbyrd3[14];
pr2MPC_FLOAT* pr2MPC_grad_cost3 = pr2MPC_grad_cost + 42;
pr2MPC_FLOAT* pr2MPC_grad_eq3 = pr2MPC_grad_eq + 42;
pr2MPC_FLOAT* pr2MPC_grad_ineq3 = pr2MPC_grad_ineq + 42;
pr2MPC_FLOAT pr2MPC_ctv3[14];
pr2MPC_FLOAT* pr2MPC_v3 = pr2MPC_v + 21;
pr2MPC_FLOAT pr2MPC_re3[7];
pr2MPC_FLOAT pr2MPC_beta3[7];
pr2MPC_FLOAT pr2MPC_betacc3[7];
pr2MPC_FLOAT* pr2MPC_dvaff3 = pr2MPC_dv_aff + 21;
pr2MPC_FLOAT* pr2MPC_dvcc3 = pr2MPC_dv_cc + 21;
pr2MPC_FLOAT pr2MPC_V3[98];
pr2MPC_FLOAT pr2MPC_Yd3[28];
pr2MPC_FLOAT pr2MPC_Ld3[28];
pr2MPC_FLOAT pr2MPC_yy3[7];
pr2MPC_FLOAT pr2MPC_bmy3[7];
pr2MPC_FLOAT pr2MPC_c3[7] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
int pr2MPC_lbIdx3[14] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13};
pr2MPC_FLOAT* pr2MPC_llb3 = pr2MPC_l + 84;
pr2MPC_FLOAT* pr2MPC_slb3 = pr2MPC_s + 84;
pr2MPC_FLOAT* pr2MPC_llbbyslb3 = pr2MPC_lbys + 84;
pr2MPC_FLOAT pr2MPC_rilb3[14];
pr2MPC_FLOAT* pr2MPC_dllbaff3 = pr2MPC_dl_aff + 84;
pr2MPC_FLOAT* pr2MPC_dslbaff3 = pr2MPC_ds_aff + 84;
pr2MPC_FLOAT* pr2MPC_dllbcc3 = pr2MPC_dl_cc + 84;
pr2MPC_FLOAT* pr2MPC_dslbcc3 = pr2MPC_ds_cc + 84;
pr2MPC_FLOAT* pr2MPC_ccrhsl3 = pr2MPC_ccrhs + 84;
int pr2MPC_ubIdx3[14] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13};
pr2MPC_FLOAT* pr2MPC_lub3 = pr2MPC_l + 98;
pr2MPC_FLOAT* pr2MPC_sub3 = pr2MPC_s + 98;
pr2MPC_FLOAT* pr2MPC_lubbysub3 = pr2MPC_lbys + 98;
pr2MPC_FLOAT pr2MPC_riub3[14];
pr2MPC_FLOAT* pr2MPC_dlubaff3 = pr2MPC_dl_aff + 98;
pr2MPC_FLOAT* pr2MPC_dsubaff3 = pr2MPC_ds_aff + 98;
pr2MPC_FLOAT* pr2MPC_dlubcc3 = pr2MPC_dl_cc + 98;
pr2MPC_FLOAT* pr2MPC_dsubcc3 = pr2MPC_ds_cc + 98;
pr2MPC_FLOAT* pr2MPC_ccrhsub3 = pr2MPC_ccrhs + 98;
pr2MPC_FLOAT pr2MPC_Phi3[14];
pr2MPC_FLOAT pr2MPC_W3[14];
pr2MPC_FLOAT pr2MPC_Ysd3[49];
pr2MPC_FLOAT pr2MPC_Lsd3[49];
pr2MPC_FLOAT* pr2MPC_z4 = pr2MPC_z + 56;
pr2MPC_FLOAT* pr2MPC_dzaff4 = pr2MPC_dz_aff + 56;
pr2MPC_FLOAT* pr2MPC_dzcc4 = pr2MPC_dz_cc + 56;
pr2MPC_FLOAT* pr2MPC_rd4 = pr2MPC_rd + 56;
pr2MPC_FLOAT pr2MPC_Lbyrd4[7];
pr2MPC_FLOAT* pr2MPC_grad_cost4 = pr2MPC_grad_cost + 56;
pr2MPC_FLOAT* pr2MPC_grad_eq4 = pr2MPC_grad_eq + 56;
pr2MPC_FLOAT* pr2MPC_grad_ineq4 = pr2MPC_grad_ineq + 56;
pr2MPC_FLOAT pr2MPC_ctv4[7];
pr2MPC_FLOAT* pr2MPC_v4 = pr2MPC_v + 28;
pr2MPC_FLOAT pr2MPC_re4[7];
pr2MPC_FLOAT pr2MPC_beta4[7];
pr2MPC_FLOAT pr2MPC_betacc4[7];
pr2MPC_FLOAT* pr2MPC_dvaff4 = pr2MPC_dv_aff + 28;
pr2MPC_FLOAT* pr2MPC_dvcc4 = pr2MPC_dv_cc + 28;
pr2MPC_FLOAT pr2MPC_V4[49];
pr2MPC_FLOAT pr2MPC_Yd4[28];
pr2MPC_FLOAT pr2MPC_Ld4[28];
pr2MPC_FLOAT pr2MPC_yy4[7];
pr2MPC_FLOAT pr2MPC_bmy4[7];
pr2MPC_FLOAT pr2MPC_c4[7] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
int pr2MPC_lbIdx4[7] = {0, 1, 2, 3, 4, 5, 6};
pr2MPC_FLOAT* pr2MPC_llb4 = pr2MPC_l + 112;
pr2MPC_FLOAT* pr2MPC_slb4 = pr2MPC_s + 112;
pr2MPC_FLOAT* pr2MPC_llbbyslb4 = pr2MPC_lbys + 112;
pr2MPC_FLOAT pr2MPC_rilb4[7];
pr2MPC_FLOAT* pr2MPC_dllbaff4 = pr2MPC_dl_aff + 112;
pr2MPC_FLOAT* pr2MPC_dslbaff4 = pr2MPC_ds_aff + 112;
pr2MPC_FLOAT* pr2MPC_dllbcc4 = pr2MPC_dl_cc + 112;
pr2MPC_FLOAT* pr2MPC_dslbcc4 = pr2MPC_ds_cc + 112;
pr2MPC_FLOAT* pr2MPC_ccrhsl4 = pr2MPC_ccrhs + 112;
int pr2MPC_ubIdx4[7] = {0, 1, 2, 3, 4, 5, 6};
pr2MPC_FLOAT* pr2MPC_lub4 = pr2MPC_l + 119;
pr2MPC_FLOAT* pr2MPC_sub4 = pr2MPC_s + 119;
pr2MPC_FLOAT* pr2MPC_lubbysub4 = pr2MPC_lbys + 119;
pr2MPC_FLOAT pr2MPC_riub4[7];
pr2MPC_FLOAT* pr2MPC_dlubaff4 = pr2MPC_dl_aff + 119;
pr2MPC_FLOAT* pr2MPC_dsubaff4 = pr2MPC_ds_aff + 119;
pr2MPC_FLOAT* pr2MPC_dlubcc4 = pr2MPC_dl_cc + 119;
pr2MPC_FLOAT* pr2MPC_dsubcc4 = pr2MPC_ds_cc + 119;
pr2MPC_FLOAT* pr2MPC_ccrhsub4 = pr2MPC_ccrhs + 119;
pr2MPC_FLOAT pr2MPC_Phi4[7];
pr2MPC_FLOAT pr2MPC_D4[7] = {-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000};
pr2MPC_FLOAT pr2MPC_W4[7];
pr2MPC_FLOAT pr2MPC_Ysd4[49];
pr2MPC_FLOAT pr2MPC_Lsd4[49];
pr2MPC_FLOAT musigma;
pr2MPC_FLOAT sigma_3rdroot;
pr2MPC_FLOAT pr2MPC_Diag1_0[14];
pr2MPC_FLOAT pr2MPC_Diag2_0[14];
pr2MPC_FLOAT pr2MPC_L_0[91];




/* SOLVER CODE --------------------------------------------------------- */
int pr2MPC_solve(pr2MPC_params* params, pr2MPC_output* output, pr2MPC_info* info)
{	
int exitcode;

#if pr2MPC_SET_TIMING == 1
	pr2MPC_timer solvertimer;
	pr2MPC_tic(&solvertimer);
#endif
/* FUNCTION CALLS INTO LA LIBRARY -------------------------------------- */
info->it = 0;
pr2MPC_LA_INITIALIZEVECTOR_63(pr2MPC_z, 0);
pr2MPC_LA_INITIALIZEVECTOR_35(pr2MPC_v, 1);
pr2MPC_LA_INITIALIZEVECTOR_126(pr2MPC_l, 10);
pr2MPC_LA_INITIALIZEVECTOR_126(pr2MPC_s, 10);
info->mu = 0;
pr2MPC_LA_DOTACC_126(pr2MPC_l, pr2MPC_s, &info->mu);
info->mu /= 126;
while( 1 ){
info->pobj = 0;
pr2MPC_LA_DIAG_QUADFCN_14(params->H1, params->f1, pr2MPC_z0, pr2MPC_grad_cost0, &info->pobj);
pr2MPC_LA_DIAG_QUADFCN_14(params->H2, params->f2, pr2MPC_z1, pr2MPC_grad_cost1, &info->pobj);
pr2MPC_LA_DIAG_QUADFCN_14(params->H3, params->f3, pr2MPC_z2, pr2MPC_grad_cost2, &info->pobj);
pr2MPC_LA_DIAG_QUADFCN_14(params->H4, params->f4, pr2MPC_z3, pr2MPC_grad_cost3, &info->pobj);
pr2MPC_LA_DIAG_QUADFCN_7(params->H5, params->f5, pr2MPC_z4, pr2MPC_grad_cost4, &info->pobj);
info->res_eq = 0;
info->dgap = 0;
pr2MPC_LA_DIAGZERO_MVMSUB6_7(pr2MPC_D0, pr2MPC_z0, params->c1, pr2MPC_v0, pr2MPC_re0, &info->dgap, &info->res_eq);
pr2MPC_LA_DENSE_DIAGZERO_MVMSUB3_7_14_14(pr2MPC_C0, pr2MPC_z0, pr2MPC_D1, pr2MPC_z1, pr2MPC_c1, pr2MPC_v1, pr2MPC_re1, &info->dgap, &info->res_eq);
pr2MPC_LA_DENSE_DIAGZERO_MVMSUB3_7_14_14(pr2MPC_C0, pr2MPC_z1, pr2MPC_D1, pr2MPC_z2, pr2MPC_c2, pr2MPC_v2, pr2MPC_re2, &info->dgap, &info->res_eq);
pr2MPC_LA_DENSE_DIAGZERO_MVMSUB3_7_14_14(pr2MPC_C0, pr2MPC_z2, pr2MPC_D1, pr2MPC_z3, pr2MPC_c3, pr2MPC_v3, pr2MPC_re3, &info->dgap, &info->res_eq);
pr2MPC_LA_DENSE_DIAGZERO_MVMSUB3_7_14_7(pr2MPC_C0, pr2MPC_z3, pr2MPC_D4, pr2MPC_z4, pr2MPC_c4, pr2MPC_v4, pr2MPC_re4, &info->dgap, &info->res_eq);
pr2MPC_LA_DENSE_DIAGZERO_MTVM2_7_14_7(pr2MPC_C0, pr2MPC_v1, pr2MPC_D0, pr2MPC_v0, pr2MPC_grad_eq0);
pr2MPC_LA_DENSE_DIAGZERO_MTVM2_7_14_7(pr2MPC_C0, pr2MPC_v2, pr2MPC_D1, pr2MPC_v1, pr2MPC_grad_eq1);
pr2MPC_LA_DENSE_DIAGZERO_MTVM2_7_14_7(pr2MPC_C0, pr2MPC_v3, pr2MPC_D1, pr2MPC_v2, pr2MPC_grad_eq2);
pr2MPC_LA_DENSE_DIAGZERO_MTVM2_7_14_7(pr2MPC_C0, pr2MPC_v4, pr2MPC_D1, pr2MPC_v3, pr2MPC_grad_eq3);
pr2MPC_LA_DIAGZERO_MTVM_7_7(pr2MPC_D4, pr2MPC_v4, pr2MPC_grad_eq4);
info->res_ineq = 0;
pr2MPC_LA_VSUBADD3_14(params->lb1, pr2MPC_z0, pr2MPC_lbIdx0, pr2MPC_llb0, pr2MPC_slb0, pr2MPC_rilb0, &info->dgap, &info->res_ineq);
pr2MPC_LA_VSUBADD2_14(pr2MPC_z0, pr2MPC_ubIdx0, params->ub1, pr2MPC_lub0, pr2MPC_sub0, pr2MPC_riub0, &info->dgap, &info->res_ineq);
pr2MPC_LA_VSUBADD3_14(params->lb2, pr2MPC_z1, pr2MPC_lbIdx1, pr2MPC_llb1, pr2MPC_slb1, pr2MPC_rilb1, &info->dgap, &info->res_ineq);
pr2MPC_LA_VSUBADD2_14(pr2MPC_z1, pr2MPC_ubIdx1, params->ub2, pr2MPC_lub1, pr2MPC_sub1, pr2MPC_riub1, &info->dgap, &info->res_ineq);
pr2MPC_LA_VSUBADD3_14(params->lb3, pr2MPC_z2, pr2MPC_lbIdx2, pr2MPC_llb2, pr2MPC_slb2, pr2MPC_rilb2, &info->dgap, &info->res_ineq);
pr2MPC_LA_VSUBADD2_14(pr2MPC_z2, pr2MPC_ubIdx2, params->ub3, pr2MPC_lub2, pr2MPC_sub2, pr2MPC_riub2, &info->dgap, &info->res_ineq);
pr2MPC_LA_VSUBADD3_14(params->lb4, pr2MPC_z3, pr2MPC_lbIdx3, pr2MPC_llb3, pr2MPC_slb3, pr2MPC_rilb3, &info->dgap, &info->res_ineq);
pr2MPC_LA_VSUBADD2_14(pr2MPC_z3, pr2MPC_ubIdx3, params->ub4, pr2MPC_lub3, pr2MPC_sub3, pr2MPC_riub3, &info->dgap, &info->res_ineq);
pr2MPC_LA_VSUBADD3_7(params->lb5, pr2MPC_z4, pr2MPC_lbIdx4, pr2MPC_llb4, pr2MPC_slb4, pr2MPC_rilb4, &info->dgap, &info->res_ineq);
pr2MPC_LA_VSUBADD2_7(pr2MPC_z4, pr2MPC_ubIdx4, params->ub5, pr2MPC_lub4, pr2MPC_sub4, pr2MPC_riub4, &info->dgap, &info->res_ineq);
pr2MPC_LA_INEQ_B_GRAD_14_14_14(pr2MPC_lub0, pr2MPC_sub0, pr2MPC_riub0, pr2MPC_llb0, pr2MPC_slb0, pr2MPC_rilb0, pr2MPC_lbIdx0, pr2MPC_ubIdx0, pr2MPC_grad_ineq0, pr2MPC_lubbysub0, pr2MPC_llbbyslb0);
pr2MPC_LA_INEQ_B_GRAD_14_14_14(pr2MPC_lub1, pr2MPC_sub1, pr2MPC_riub1, pr2MPC_llb1, pr2MPC_slb1, pr2MPC_rilb1, pr2MPC_lbIdx1, pr2MPC_ubIdx1, pr2MPC_grad_ineq1, pr2MPC_lubbysub1, pr2MPC_llbbyslb1);
pr2MPC_LA_INEQ_B_GRAD_14_14_14(pr2MPC_lub2, pr2MPC_sub2, pr2MPC_riub2, pr2MPC_llb2, pr2MPC_slb2, pr2MPC_rilb2, pr2MPC_lbIdx2, pr2MPC_ubIdx2, pr2MPC_grad_ineq2, pr2MPC_lubbysub2, pr2MPC_llbbyslb2);
pr2MPC_LA_INEQ_B_GRAD_14_14_14(pr2MPC_lub3, pr2MPC_sub3, pr2MPC_riub3, pr2MPC_llb3, pr2MPC_slb3, pr2MPC_rilb3, pr2MPC_lbIdx3, pr2MPC_ubIdx3, pr2MPC_grad_ineq3, pr2MPC_lubbysub3, pr2MPC_llbbyslb3);
pr2MPC_LA_INEQ_B_GRAD_7_7_7(pr2MPC_lub4, pr2MPC_sub4, pr2MPC_riub4, pr2MPC_llb4, pr2MPC_slb4, pr2MPC_rilb4, pr2MPC_lbIdx4, pr2MPC_ubIdx4, pr2MPC_grad_ineq4, pr2MPC_lubbysub4, pr2MPC_llbbyslb4);
info->dobj = info->pobj - info->dgap;
info->rdgap = info->pobj ? info->dgap / info->pobj : 1e6;
if( info->rdgap < 0 ) info->rdgap = -info->rdgap;
if( info->mu < pr2MPC_SET_ACC_KKTCOMPL
    && (info->rdgap < pr2MPC_SET_ACC_RDGAP || info->dgap < pr2MPC_SET_ACC_KKTCOMPL)
    && info->res_eq < pr2MPC_SET_ACC_RESEQ
    && info->res_ineq < pr2MPC_SET_ACC_RESINEQ ){
exitcode = pr2MPC_OPTIMAL; break; }
if( info->it == pr2MPC_SET_MAXIT ){
exitcode = pr2MPC_MAXITREACHED; break; }
pr2MPC_LA_VVADD3_63(pr2MPC_grad_cost, pr2MPC_grad_eq, pr2MPC_grad_ineq, pr2MPC_rd);
pr2MPC_LA_DIAG_CHOL_ONELOOP_LBUB_14_14_14(params->H1, pr2MPC_llbbyslb0, pr2MPC_lbIdx0, pr2MPC_lubbysub0, pr2MPC_ubIdx0, pr2MPC_Phi0);
pr2MPC_LA_DIAG_MATRIXFORWARDSUB_7_14(pr2MPC_Phi0, pr2MPC_C0, pr2MPC_V0);
pr2MPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_7_14(pr2MPC_Phi0, pr2MPC_D0, pr2MPC_W0);
pr2MPC_LA_DENSE_DIAGZERO_MMTM_7_14_7(pr2MPC_W0, pr2MPC_V0, pr2MPC_Ysd1);
pr2MPC_LA_DIAG_FORWARDSUB_14(pr2MPC_Phi0, pr2MPC_rd0, pr2MPC_Lbyrd0);
pr2MPC_LA_DIAG_CHOL_ONELOOP_LBUB_14_14_14(params->H2, pr2MPC_llbbyslb1, pr2MPC_lbIdx1, pr2MPC_lubbysub1, pr2MPC_ubIdx1, pr2MPC_Phi1);
pr2MPC_LA_DIAG_MATRIXFORWARDSUB_7_14(pr2MPC_Phi1, pr2MPC_C0, pr2MPC_V1);
pr2MPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_7_14(pr2MPC_Phi1, pr2MPC_D1, pr2MPC_W1);
pr2MPC_LA_DENSE_DIAGZERO_MMTM_7_14_7(pr2MPC_W1, pr2MPC_V1, pr2MPC_Ysd2);
pr2MPC_LA_DIAG_FORWARDSUB_14(pr2MPC_Phi1, pr2MPC_rd1, pr2MPC_Lbyrd1);
pr2MPC_LA_DIAG_CHOL_ONELOOP_LBUB_14_14_14(params->H3, pr2MPC_llbbyslb2, pr2MPC_lbIdx2, pr2MPC_lubbysub2, pr2MPC_ubIdx2, pr2MPC_Phi2);
pr2MPC_LA_DIAG_MATRIXFORWARDSUB_7_14(pr2MPC_Phi2, pr2MPC_C0, pr2MPC_V2);
pr2MPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_7_14(pr2MPC_Phi2, pr2MPC_D1, pr2MPC_W2);
pr2MPC_LA_DENSE_DIAGZERO_MMTM_7_14_7(pr2MPC_W2, pr2MPC_V2, pr2MPC_Ysd3);
pr2MPC_LA_DIAG_FORWARDSUB_14(pr2MPC_Phi2, pr2MPC_rd2, pr2MPC_Lbyrd2);
pr2MPC_LA_DIAG_CHOL_ONELOOP_LBUB_14_14_14(params->H4, pr2MPC_llbbyslb3, pr2MPC_lbIdx3, pr2MPC_lubbysub3, pr2MPC_ubIdx3, pr2MPC_Phi3);
pr2MPC_LA_DIAG_MATRIXFORWARDSUB_7_14(pr2MPC_Phi3, pr2MPC_C0, pr2MPC_V3);
pr2MPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_7_14(pr2MPC_Phi3, pr2MPC_D1, pr2MPC_W3);
pr2MPC_LA_DENSE_DIAGZERO_MMTM_7_14_7(pr2MPC_W3, pr2MPC_V3, pr2MPC_Ysd4);
pr2MPC_LA_DIAG_FORWARDSUB_14(pr2MPC_Phi3, pr2MPC_rd3, pr2MPC_Lbyrd3);
pr2MPC_LA_DIAG_CHOL_ONELOOP_LBUB_7_7_7(params->H5, pr2MPC_llbbyslb4, pr2MPC_lbIdx4, pr2MPC_lubbysub4, pr2MPC_ubIdx4, pr2MPC_Phi4);
pr2MPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_7_7(pr2MPC_Phi4, pr2MPC_D4, pr2MPC_W4);
pr2MPC_LA_DIAG_FORWARDSUB_7(pr2MPC_Phi4, pr2MPC_rd4, pr2MPC_Lbyrd4);
pr2MPC_LA_DIAGZERO_MMT_7(pr2MPC_W0, pr2MPC_Yd0);
pr2MPC_LA_DIAGZERO_MVMSUB7_7(pr2MPC_W0, pr2MPC_Lbyrd0, pr2MPC_re0, pr2MPC_beta0);
pr2MPC_LA_DENSE_DIAGZERO_MMT2_7_14_14(pr2MPC_V0, pr2MPC_W1, pr2MPC_Yd1);
pr2MPC_LA_DENSE_DIAGZERO_2MVMSUB2_7_14_14(pr2MPC_V0, pr2MPC_Lbyrd0, pr2MPC_W1, pr2MPC_Lbyrd1, pr2MPC_re1, pr2MPC_beta1);
pr2MPC_LA_DENSE_DIAGZERO_MMT2_7_14_14(pr2MPC_V1, pr2MPC_W2, pr2MPC_Yd2);
pr2MPC_LA_DENSE_DIAGZERO_2MVMSUB2_7_14_14(pr2MPC_V1, pr2MPC_Lbyrd1, pr2MPC_W2, pr2MPC_Lbyrd2, pr2MPC_re2, pr2MPC_beta2);
pr2MPC_LA_DENSE_DIAGZERO_MMT2_7_14_14(pr2MPC_V2, pr2MPC_W3, pr2MPC_Yd3);
pr2MPC_LA_DENSE_DIAGZERO_2MVMSUB2_7_14_14(pr2MPC_V2, pr2MPC_Lbyrd2, pr2MPC_W3, pr2MPC_Lbyrd3, pr2MPC_re3, pr2MPC_beta3);
pr2MPC_LA_DENSE_DIAGZERO_MMT2_7_14_7(pr2MPC_V3, pr2MPC_W4, pr2MPC_Yd4);
pr2MPC_LA_DENSE_DIAGZERO_2MVMSUB2_7_14_7(pr2MPC_V3, pr2MPC_Lbyrd3, pr2MPC_W4, pr2MPC_Lbyrd4, pr2MPC_re4, pr2MPC_beta4);
pr2MPC_LA_DENSE_CHOL_7(pr2MPC_Yd0, pr2MPC_Ld0);
pr2MPC_LA_DENSE_FORWARDSUB_7(pr2MPC_Ld0, pr2MPC_beta0, pr2MPC_yy0);
pr2MPC_LA_DENSE_MATRIXTFORWARDSUB_7_7(pr2MPC_Ld0, pr2MPC_Ysd1, pr2MPC_Lsd1);
pr2MPC_LA_DENSE_MMTSUB_7_7(pr2MPC_Lsd1, pr2MPC_Yd1);
pr2MPC_LA_DENSE_CHOL_7(pr2MPC_Yd1, pr2MPC_Ld1);
pr2MPC_LA_DENSE_MVMSUB1_7_7(pr2MPC_Lsd1, pr2MPC_yy0, pr2MPC_beta1, pr2MPC_bmy1);
pr2MPC_LA_DENSE_FORWARDSUB_7(pr2MPC_Ld1, pr2MPC_bmy1, pr2MPC_yy1);
pr2MPC_LA_DENSE_MATRIXTFORWARDSUB_7_7(pr2MPC_Ld1, pr2MPC_Ysd2, pr2MPC_Lsd2);
pr2MPC_LA_DENSE_MMTSUB_7_7(pr2MPC_Lsd2, pr2MPC_Yd2);
pr2MPC_LA_DENSE_CHOL_7(pr2MPC_Yd2, pr2MPC_Ld2);
pr2MPC_LA_DENSE_MVMSUB1_7_7(pr2MPC_Lsd2, pr2MPC_yy1, pr2MPC_beta2, pr2MPC_bmy2);
pr2MPC_LA_DENSE_FORWARDSUB_7(pr2MPC_Ld2, pr2MPC_bmy2, pr2MPC_yy2);
pr2MPC_LA_DENSE_MATRIXTFORWARDSUB_7_7(pr2MPC_Ld2, pr2MPC_Ysd3, pr2MPC_Lsd3);
pr2MPC_LA_DENSE_MMTSUB_7_7(pr2MPC_Lsd3, pr2MPC_Yd3);
pr2MPC_LA_DENSE_CHOL_7(pr2MPC_Yd3, pr2MPC_Ld3);
pr2MPC_LA_DENSE_MVMSUB1_7_7(pr2MPC_Lsd3, pr2MPC_yy2, pr2MPC_beta3, pr2MPC_bmy3);
pr2MPC_LA_DENSE_FORWARDSUB_7(pr2MPC_Ld3, pr2MPC_bmy3, pr2MPC_yy3);
pr2MPC_LA_DENSE_MATRIXTFORWARDSUB_7_7(pr2MPC_Ld3, pr2MPC_Ysd4, pr2MPC_Lsd4);
pr2MPC_LA_DENSE_MMTSUB_7_7(pr2MPC_Lsd4, pr2MPC_Yd4);
pr2MPC_LA_DENSE_CHOL_7(pr2MPC_Yd4, pr2MPC_Ld4);
pr2MPC_LA_DENSE_MVMSUB1_7_7(pr2MPC_Lsd4, pr2MPC_yy3, pr2MPC_beta4, pr2MPC_bmy4);
pr2MPC_LA_DENSE_FORWARDSUB_7(pr2MPC_Ld4, pr2MPC_bmy4, pr2MPC_yy4);
pr2MPC_LA_DENSE_BACKWARDSUB_7(pr2MPC_Ld4, pr2MPC_yy4, pr2MPC_dvaff4);
pr2MPC_LA_DENSE_MTVMSUB_7_7(pr2MPC_Lsd4, pr2MPC_dvaff4, pr2MPC_yy3, pr2MPC_bmy3);
pr2MPC_LA_DENSE_BACKWARDSUB_7(pr2MPC_Ld3, pr2MPC_bmy3, pr2MPC_dvaff3);
pr2MPC_LA_DENSE_MTVMSUB_7_7(pr2MPC_Lsd3, pr2MPC_dvaff3, pr2MPC_yy2, pr2MPC_bmy2);
pr2MPC_LA_DENSE_BACKWARDSUB_7(pr2MPC_Ld2, pr2MPC_bmy2, pr2MPC_dvaff2);
pr2MPC_LA_DENSE_MTVMSUB_7_7(pr2MPC_Lsd2, pr2MPC_dvaff2, pr2MPC_yy1, pr2MPC_bmy1);
pr2MPC_LA_DENSE_BACKWARDSUB_7(pr2MPC_Ld1, pr2MPC_bmy1, pr2MPC_dvaff1);
pr2MPC_LA_DENSE_MTVMSUB_7_7(pr2MPC_Lsd1, pr2MPC_dvaff1, pr2MPC_yy0, pr2MPC_bmy0);
pr2MPC_LA_DENSE_BACKWARDSUB_7(pr2MPC_Ld0, pr2MPC_bmy0, pr2MPC_dvaff0);
pr2MPC_LA_DENSE_DIAGZERO_MTVM2_7_14_7(pr2MPC_C0, pr2MPC_dvaff1, pr2MPC_D0, pr2MPC_dvaff0, pr2MPC_grad_eq0);
pr2MPC_LA_DENSE_DIAGZERO_MTVM2_7_14_7(pr2MPC_C0, pr2MPC_dvaff2, pr2MPC_D1, pr2MPC_dvaff1, pr2MPC_grad_eq1);
pr2MPC_LA_DENSE_DIAGZERO_MTVM2_7_14_7(pr2MPC_C0, pr2MPC_dvaff3, pr2MPC_D1, pr2MPC_dvaff2, pr2MPC_grad_eq2);
pr2MPC_LA_DENSE_DIAGZERO_MTVM2_7_14_7(pr2MPC_C0, pr2MPC_dvaff4, pr2MPC_D1, pr2MPC_dvaff3, pr2MPC_grad_eq3);
pr2MPC_LA_DIAGZERO_MTVM_7_7(pr2MPC_D4, pr2MPC_dvaff4, pr2MPC_grad_eq4);
pr2MPC_LA_VSUB2_63(pr2MPC_rd, pr2MPC_grad_eq, pr2MPC_rd);
pr2MPC_LA_DIAG_FORWARDBACKWARDSUB_14(pr2MPC_Phi0, pr2MPC_rd0, pr2MPC_dzaff0);
pr2MPC_LA_DIAG_FORWARDBACKWARDSUB_14(pr2MPC_Phi1, pr2MPC_rd1, pr2MPC_dzaff1);
pr2MPC_LA_DIAG_FORWARDBACKWARDSUB_14(pr2MPC_Phi2, pr2MPC_rd2, pr2MPC_dzaff2);
pr2MPC_LA_DIAG_FORWARDBACKWARDSUB_14(pr2MPC_Phi3, pr2MPC_rd3, pr2MPC_dzaff3);
pr2MPC_LA_DIAG_FORWARDBACKWARDSUB_7(pr2MPC_Phi4, pr2MPC_rd4, pr2MPC_dzaff4);
pr2MPC_LA_VSUB_INDEXED_14(pr2MPC_dzaff0, pr2MPC_lbIdx0, pr2MPC_rilb0, pr2MPC_dslbaff0);
pr2MPC_LA_VSUB3_14(pr2MPC_llbbyslb0, pr2MPC_dslbaff0, pr2MPC_llb0, pr2MPC_dllbaff0);
pr2MPC_LA_VSUB2_INDEXED_14(pr2MPC_riub0, pr2MPC_dzaff0, pr2MPC_ubIdx0, pr2MPC_dsubaff0);
pr2MPC_LA_VSUB3_14(pr2MPC_lubbysub0, pr2MPC_dsubaff0, pr2MPC_lub0, pr2MPC_dlubaff0);
pr2MPC_LA_VSUB_INDEXED_14(pr2MPC_dzaff1, pr2MPC_lbIdx1, pr2MPC_rilb1, pr2MPC_dslbaff1);
pr2MPC_LA_VSUB3_14(pr2MPC_llbbyslb1, pr2MPC_dslbaff1, pr2MPC_llb1, pr2MPC_dllbaff1);
pr2MPC_LA_VSUB2_INDEXED_14(pr2MPC_riub1, pr2MPC_dzaff1, pr2MPC_ubIdx1, pr2MPC_dsubaff1);
pr2MPC_LA_VSUB3_14(pr2MPC_lubbysub1, pr2MPC_dsubaff1, pr2MPC_lub1, pr2MPC_dlubaff1);
pr2MPC_LA_VSUB_INDEXED_14(pr2MPC_dzaff2, pr2MPC_lbIdx2, pr2MPC_rilb2, pr2MPC_dslbaff2);
pr2MPC_LA_VSUB3_14(pr2MPC_llbbyslb2, pr2MPC_dslbaff2, pr2MPC_llb2, pr2MPC_dllbaff2);
pr2MPC_LA_VSUB2_INDEXED_14(pr2MPC_riub2, pr2MPC_dzaff2, pr2MPC_ubIdx2, pr2MPC_dsubaff2);
pr2MPC_LA_VSUB3_14(pr2MPC_lubbysub2, pr2MPC_dsubaff2, pr2MPC_lub2, pr2MPC_dlubaff2);
pr2MPC_LA_VSUB_INDEXED_14(pr2MPC_dzaff3, pr2MPC_lbIdx3, pr2MPC_rilb3, pr2MPC_dslbaff3);
pr2MPC_LA_VSUB3_14(pr2MPC_llbbyslb3, pr2MPC_dslbaff3, pr2MPC_llb3, pr2MPC_dllbaff3);
pr2MPC_LA_VSUB2_INDEXED_14(pr2MPC_riub3, pr2MPC_dzaff3, pr2MPC_ubIdx3, pr2MPC_dsubaff3);
pr2MPC_LA_VSUB3_14(pr2MPC_lubbysub3, pr2MPC_dsubaff3, pr2MPC_lub3, pr2MPC_dlubaff3);
pr2MPC_LA_VSUB_INDEXED_7(pr2MPC_dzaff4, pr2MPC_lbIdx4, pr2MPC_rilb4, pr2MPC_dslbaff4);
pr2MPC_LA_VSUB3_7(pr2MPC_llbbyslb4, pr2MPC_dslbaff4, pr2MPC_llb4, pr2MPC_dllbaff4);
pr2MPC_LA_VSUB2_INDEXED_7(pr2MPC_riub4, pr2MPC_dzaff4, pr2MPC_ubIdx4, pr2MPC_dsubaff4);
pr2MPC_LA_VSUB3_7(pr2MPC_lubbysub4, pr2MPC_dsubaff4, pr2MPC_lub4, pr2MPC_dlubaff4);
info->lsit_aff = pr2MPC_LINESEARCH_BACKTRACKING_AFFINE(pr2MPC_l, pr2MPC_s, pr2MPC_dl_aff, pr2MPC_ds_aff, &info->step_aff, &info->mu_aff);
if( info->lsit_aff == pr2MPC_NOPROGRESS ){
exitcode = pr2MPC_NOPROGRESS; break;
}
sigma_3rdroot = info->mu_aff / info->mu;
info->sigma = sigma_3rdroot*sigma_3rdroot*sigma_3rdroot;
musigma = info->mu * info->sigma;
pr2MPC_LA_VSUB5_126(pr2MPC_ds_aff, pr2MPC_dl_aff, info->mu, info->sigma, pr2MPC_ccrhs);
pr2MPC_LA_VSUB6_INDEXED_14_14_14(pr2MPC_ccrhsub0, pr2MPC_sub0, pr2MPC_ubIdx0, pr2MPC_ccrhsl0, pr2MPC_slb0, pr2MPC_lbIdx0, pr2MPC_rd0);
pr2MPC_LA_VSUB6_INDEXED_14_14_14(pr2MPC_ccrhsub1, pr2MPC_sub1, pr2MPC_ubIdx1, pr2MPC_ccrhsl1, pr2MPC_slb1, pr2MPC_lbIdx1, pr2MPC_rd1);
pr2MPC_LA_DIAG_FORWARDSUB_14(pr2MPC_Phi0, pr2MPC_rd0, pr2MPC_Lbyrd0);
pr2MPC_LA_DIAG_FORWARDSUB_14(pr2MPC_Phi1, pr2MPC_rd1, pr2MPC_Lbyrd1);
pr2MPC_LA_DIAGZERO_MVM_7(pr2MPC_W0, pr2MPC_Lbyrd0, pr2MPC_beta0);
pr2MPC_LA_DENSE_FORWARDSUB_7(pr2MPC_Ld0, pr2MPC_beta0, pr2MPC_yy0);
pr2MPC_LA_DENSE_DIAGZERO_2MVMADD_7_14_14(pr2MPC_V0, pr2MPC_Lbyrd0, pr2MPC_W1, pr2MPC_Lbyrd1, pr2MPC_beta1);
pr2MPC_LA_DENSE_MVMSUB1_7_7(pr2MPC_Lsd1, pr2MPC_yy0, pr2MPC_beta1, pr2MPC_bmy1);
pr2MPC_LA_DENSE_FORWARDSUB_7(pr2MPC_Ld1, pr2MPC_bmy1, pr2MPC_yy1);
pr2MPC_LA_VSUB6_INDEXED_14_14_14(pr2MPC_ccrhsub2, pr2MPC_sub2, pr2MPC_ubIdx2, pr2MPC_ccrhsl2, pr2MPC_slb2, pr2MPC_lbIdx2, pr2MPC_rd2);
pr2MPC_LA_DIAG_FORWARDSUB_14(pr2MPC_Phi2, pr2MPC_rd2, pr2MPC_Lbyrd2);
pr2MPC_LA_DENSE_DIAGZERO_2MVMADD_7_14_14(pr2MPC_V1, pr2MPC_Lbyrd1, pr2MPC_W2, pr2MPC_Lbyrd2, pr2MPC_beta2);
pr2MPC_LA_DENSE_MVMSUB1_7_7(pr2MPC_Lsd2, pr2MPC_yy1, pr2MPC_beta2, pr2MPC_bmy2);
pr2MPC_LA_DENSE_FORWARDSUB_7(pr2MPC_Ld2, pr2MPC_bmy2, pr2MPC_yy2);
pr2MPC_LA_VSUB6_INDEXED_14_14_14(pr2MPC_ccrhsub3, pr2MPC_sub3, pr2MPC_ubIdx3, pr2MPC_ccrhsl3, pr2MPC_slb3, pr2MPC_lbIdx3, pr2MPC_rd3);
pr2MPC_LA_DIAG_FORWARDSUB_14(pr2MPC_Phi3, pr2MPC_rd3, pr2MPC_Lbyrd3);
pr2MPC_LA_DENSE_DIAGZERO_2MVMADD_7_14_14(pr2MPC_V2, pr2MPC_Lbyrd2, pr2MPC_W3, pr2MPC_Lbyrd3, pr2MPC_beta3);
pr2MPC_LA_DENSE_MVMSUB1_7_7(pr2MPC_Lsd3, pr2MPC_yy2, pr2MPC_beta3, pr2MPC_bmy3);
pr2MPC_LA_DENSE_FORWARDSUB_7(pr2MPC_Ld3, pr2MPC_bmy3, pr2MPC_yy3);
pr2MPC_LA_VSUB6_INDEXED_7_7_7(pr2MPC_ccrhsub4, pr2MPC_sub4, pr2MPC_ubIdx4, pr2MPC_ccrhsl4, pr2MPC_slb4, pr2MPC_lbIdx4, pr2MPC_rd4);
pr2MPC_LA_DIAG_FORWARDSUB_7(pr2MPC_Phi4, pr2MPC_rd4, pr2MPC_Lbyrd4);
pr2MPC_LA_DENSE_DIAGZERO_2MVMADD_7_14_7(pr2MPC_V3, pr2MPC_Lbyrd3, pr2MPC_W4, pr2MPC_Lbyrd4, pr2MPC_beta4);
pr2MPC_LA_DENSE_MVMSUB1_7_7(pr2MPC_Lsd4, pr2MPC_yy3, pr2MPC_beta4, pr2MPC_bmy4);
pr2MPC_LA_DENSE_FORWARDSUB_7(pr2MPC_Ld4, pr2MPC_bmy4, pr2MPC_yy4);
pr2MPC_LA_DENSE_BACKWARDSUB_7(pr2MPC_Ld4, pr2MPC_yy4, pr2MPC_dvcc4);
pr2MPC_LA_DENSE_MTVMSUB_7_7(pr2MPC_Lsd4, pr2MPC_dvcc4, pr2MPC_yy3, pr2MPC_bmy3);
pr2MPC_LA_DENSE_BACKWARDSUB_7(pr2MPC_Ld3, pr2MPC_bmy3, pr2MPC_dvcc3);
pr2MPC_LA_DENSE_MTVMSUB_7_7(pr2MPC_Lsd3, pr2MPC_dvcc3, pr2MPC_yy2, pr2MPC_bmy2);
pr2MPC_LA_DENSE_BACKWARDSUB_7(pr2MPC_Ld2, pr2MPC_bmy2, pr2MPC_dvcc2);
pr2MPC_LA_DENSE_MTVMSUB_7_7(pr2MPC_Lsd2, pr2MPC_dvcc2, pr2MPC_yy1, pr2MPC_bmy1);
pr2MPC_LA_DENSE_BACKWARDSUB_7(pr2MPC_Ld1, pr2MPC_bmy1, pr2MPC_dvcc1);
pr2MPC_LA_DENSE_MTVMSUB_7_7(pr2MPC_Lsd1, pr2MPC_dvcc1, pr2MPC_yy0, pr2MPC_bmy0);
pr2MPC_LA_DENSE_BACKWARDSUB_7(pr2MPC_Ld0, pr2MPC_bmy0, pr2MPC_dvcc0);
pr2MPC_LA_DENSE_DIAGZERO_MTVM2_7_14_7(pr2MPC_C0, pr2MPC_dvcc1, pr2MPC_D0, pr2MPC_dvcc0, pr2MPC_grad_eq0);
pr2MPC_LA_DENSE_DIAGZERO_MTVM2_7_14_7(pr2MPC_C0, pr2MPC_dvcc2, pr2MPC_D1, pr2MPC_dvcc1, pr2MPC_grad_eq1);
pr2MPC_LA_DENSE_DIAGZERO_MTVM2_7_14_7(pr2MPC_C0, pr2MPC_dvcc3, pr2MPC_D1, pr2MPC_dvcc2, pr2MPC_grad_eq2);
pr2MPC_LA_DENSE_DIAGZERO_MTVM2_7_14_7(pr2MPC_C0, pr2MPC_dvcc4, pr2MPC_D1, pr2MPC_dvcc3, pr2MPC_grad_eq3);
pr2MPC_LA_DIAGZERO_MTVM_7_7(pr2MPC_D4, pr2MPC_dvcc4, pr2MPC_grad_eq4);
pr2MPC_LA_VSUB_63(pr2MPC_rd, pr2MPC_grad_eq, pr2MPC_rd);
pr2MPC_LA_DIAG_FORWARDBACKWARDSUB_14(pr2MPC_Phi0, pr2MPC_rd0, pr2MPC_dzcc0);
pr2MPC_LA_DIAG_FORWARDBACKWARDSUB_14(pr2MPC_Phi1, pr2MPC_rd1, pr2MPC_dzcc1);
pr2MPC_LA_DIAG_FORWARDBACKWARDSUB_14(pr2MPC_Phi2, pr2MPC_rd2, pr2MPC_dzcc2);
pr2MPC_LA_DIAG_FORWARDBACKWARDSUB_14(pr2MPC_Phi3, pr2MPC_rd3, pr2MPC_dzcc3);
pr2MPC_LA_DIAG_FORWARDBACKWARDSUB_7(pr2MPC_Phi4, pr2MPC_rd4, pr2MPC_dzcc4);
pr2MPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_14(pr2MPC_ccrhsl0, pr2MPC_slb0, pr2MPC_llbbyslb0, pr2MPC_dzcc0, pr2MPC_lbIdx0, pr2MPC_dllbcc0);
pr2MPC_LA_VEC_DIVSUB_MULTADD_INDEXED_14(pr2MPC_ccrhsub0, pr2MPC_sub0, pr2MPC_lubbysub0, pr2MPC_dzcc0, pr2MPC_ubIdx0, pr2MPC_dlubcc0);
pr2MPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_14(pr2MPC_ccrhsl1, pr2MPC_slb1, pr2MPC_llbbyslb1, pr2MPC_dzcc1, pr2MPC_lbIdx1, pr2MPC_dllbcc1);
pr2MPC_LA_VEC_DIVSUB_MULTADD_INDEXED_14(pr2MPC_ccrhsub1, pr2MPC_sub1, pr2MPC_lubbysub1, pr2MPC_dzcc1, pr2MPC_ubIdx1, pr2MPC_dlubcc1);
pr2MPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_14(pr2MPC_ccrhsl2, pr2MPC_slb2, pr2MPC_llbbyslb2, pr2MPC_dzcc2, pr2MPC_lbIdx2, pr2MPC_dllbcc2);
pr2MPC_LA_VEC_DIVSUB_MULTADD_INDEXED_14(pr2MPC_ccrhsub2, pr2MPC_sub2, pr2MPC_lubbysub2, pr2MPC_dzcc2, pr2MPC_ubIdx2, pr2MPC_dlubcc2);
pr2MPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_14(pr2MPC_ccrhsl3, pr2MPC_slb3, pr2MPC_llbbyslb3, pr2MPC_dzcc3, pr2MPC_lbIdx3, pr2MPC_dllbcc3);
pr2MPC_LA_VEC_DIVSUB_MULTADD_INDEXED_14(pr2MPC_ccrhsub3, pr2MPC_sub3, pr2MPC_lubbysub3, pr2MPC_dzcc3, pr2MPC_ubIdx3, pr2MPC_dlubcc3);
pr2MPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_7(pr2MPC_ccrhsl4, pr2MPC_slb4, pr2MPC_llbbyslb4, pr2MPC_dzcc4, pr2MPC_lbIdx4, pr2MPC_dllbcc4);
pr2MPC_LA_VEC_DIVSUB_MULTADD_INDEXED_7(pr2MPC_ccrhsub4, pr2MPC_sub4, pr2MPC_lubbysub4, pr2MPC_dzcc4, pr2MPC_ubIdx4, pr2MPC_dlubcc4);
pr2MPC_LA_VSUB7_126(pr2MPC_l, pr2MPC_ccrhs, pr2MPC_s, pr2MPC_dl_cc, pr2MPC_ds_cc);
pr2MPC_LA_VADD_63(pr2MPC_dz_cc, pr2MPC_dz_aff);
pr2MPC_LA_VADD_35(pr2MPC_dv_cc, pr2MPC_dv_aff);
pr2MPC_LA_VADD_126(pr2MPC_dl_cc, pr2MPC_dl_aff);
pr2MPC_LA_VADD_126(pr2MPC_ds_cc, pr2MPC_ds_aff);
info->lsit_cc = pr2MPC_LINESEARCH_BACKTRACKING_COMBINED(pr2MPC_z, pr2MPC_v, pr2MPC_l, pr2MPC_s, pr2MPC_dz_cc, pr2MPC_dv_cc, pr2MPC_dl_cc, pr2MPC_ds_cc, &info->step_cc, &info->mu);
if( info->lsit_cc == pr2MPC_NOPROGRESS ){
exitcode = pr2MPC_NOPROGRESS; break;
}
info->it++;
}
output->z1[0] = pr2MPC_z0[0];
output->z1[1] = pr2MPC_z0[1];
output->z1[2] = pr2MPC_z0[2];
output->z1[3] = pr2MPC_z0[3];
output->z1[4] = pr2MPC_z0[4];
output->z1[5] = pr2MPC_z0[5];
output->z1[6] = pr2MPC_z0[6];
output->z1[7] = pr2MPC_z0[7];
output->z1[8] = pr2MPC_z0[8];
output->z1[9] = pr2MPC_z0[9];
output->z1[10] = pr2MPC_z0[10];
output->z1[11] = pr2MPC_z0[11];
output->z1[12] = pr2MPC_z0[12];
output->z1[13] = pr2MPC_z0[13];
output->z2[0] = pr2MPC_z1[0];
output->z2[1] = pr2MPC_z1[1];
output->z2[2] = pr2MPC_z1[2];
output->z2[3] = pr2MPC_z1[3];
output->z2[4] = pr2MPC_z1[4];
output->z2[5] = pr2MPC_z1[5];
output->z2[6] = pr2MPC_z1[6];
output->z2[7] = pr2MPC_z1[7];
output->z2[8] = pr2MPC_z1[8];
output->z2[9] = pr2MPC_z1[9];
output->z2[10] = pr2MPC_z1[10];
output->z2[11] = pr2MPC_z1[11];
output->z2[12] = pr2MPC_z1[12];
output->z2[13] = pr2MPC_z1[13];
output->z3[0] = pr2MPC_z2[0];
output->z3[1] = pr2MPC_z2[1];
output->z3[2] = pr2MPC_z2[2];
output->z3[3] = pr2MPC_z2[3];
output->z3[4] = pr2MPC_z2[4];
output->z3[5] = pr2MPC_z2[5];
output->z3[6] = pr2MPC_z2[6];
output->z3[7] = pr2MPC_z2[7];
output->z3[8] = pr2MPC_z2[8];
output->z3[9] = pr2MPC_z2[9];
output->z3[10] = pr2MPC_z2[10];
output->z3[11] = pr2MPC_z2[11];
output->z3[12] = pr2MPC_z2[12];
output->z3[13] = pr2MPC_z2[13];
output->z4[0] = pr2MPC_z3[0];
output->z4[1] = pr2MPC_z3[1];
output->z4[2] = pr2MPC_z3[2];
output->z4[3] = pr2MPC_z3[3];
output->z4[4] = pr2MPC_z3[4];
output->z4[5] = pr2MPC_z3[5];
output->z4[6] = pr2MPC_z3[6];
output->z4[7] = pr2MPC_z3[7];
output->z4[8] = pr2MPC_z3[8];
output->z4[9] = pr2MPC_z3[9];
output->z4[10] = pr2MPC_z3[10];
output->z4[11] = pr2MPC_z3[11];
output->z4[12] = pr2MPC_z3[12];
output->z4[13] = pr2MPC_z3[13];
output->z5[0] = pr2MPC_z4[0];
output->z5[1] = pr2MPC_z4[1];
output->z5[2] = pr2MPC_z4[2];
output->z5[3] = pr2MPC_z4[3];
output->z5[4] = pr2MPC_z4[4];
output->z5[5] = pr2MPC_z4[5];
output->z5[6] = pr2MPC_z4[6];

#if pr2MPC_SET_TIMING == 1
info->solvetime = pr2MPC_toc(&solvertimer);
#if pr2MPC_SET_PRINTLEVEL > 0 && pr2MPC_SET_TIMING == 1
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
