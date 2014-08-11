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

#include "pr2eihMPC.h"

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
 * Initializes a vector of length 21 with a value.
 */
void pr2eihMPC_LA_INITIALIZEVECTOR_21(pr2eihMPC_FLOAT* vec, pr2eihMPC_FLOAT value)
{
	int i;
	for( i=0; i<21; i++ )
	{
		vec[i] = value;
	}
}


/*
 * Initializes a vector of length 14 with a value.
 */
void pr2eihMPC_LA_INITIALIZEVECTOR_14(pr2eihMPC_FLOAT* vec, pr2eihMPC_FLOAT value)
{
	int i;
	for( i=0; i<14; i++ )
	{
		vec[i] = value;
	}
}


/*
 * Initializes a vector of length 42 with a value.
 */
void pr2eihMPC_LA_INITIALIZEVECTOR_42(pr2eihMPC_FLOAT* vec, pr2eihMPC_FLOAT value)
{
	int i;
	for( i=0; i<42; i++ )
	{
		vec[i] = value;
	}
}


/* 
 * Calculates a dot product and adds it to a variable: z += x'*y; 
 * This function is for vectors of length 42.
 */
void pr2eihMPC_LA_DOTACC_42(pr2eihMPC_FLOAT *x, pr2eihMPC_FLOAT *y, pr2eihMPC_FLOAT *z)
{
	int i;
	for( i=0; i<42; i++ ){
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
void pr2eihMPC_LA_DIAG_QUADFCN_14(pr2eihMPC_FLOAT* H, pr2eihMPC_FLOAT* f, pr2eihMPC_FLOAT* z, pr2eihMPC_FLOAT* grad, pr2eihMPC_FLOAT* value)
{
	int i;
	pr2eihMPC_FLOAT hz;	
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
void pr2eihMPC_LA_DIAG_QUADFCN_7(pr2eihMPC_FLOAT* H, pr2eihMPC_FLOAT* f, pr2eihMPC_FLOAT* z, pr2eihMPC_FLOAT* grad, pr2eihMPC_FLOAT* value)
{
	int i;
	pr2eihMPC_FLOAT hz;	
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
void pr2eihMPC_LA_DIAGZERO_MVMSUB6_7(pr2eihMPC_FLOAT *B, pr2eihMPC_FLOAT *u, pr2eihMPC_FLOAT *b, pr2eihMPC_FLOAT *l, pr2eihMPC_FLOAT *r, pr2eihMPC_FLOAT *z, pr2eihMPC_FLOAT *y)
{
	int i;
	pr2eihMPC_FLOAT Bu[7];
	pr2eihMPC_FLOAT norm = *y;
	pr2eihMPC_FLOAT lr = 0;

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
void pr2eihMPC_LA_DENSE_DIAGZERO_MVMSUB3_7_14_7(pr2eihMPC_FLOAT *A, pr2eihMPC_FLOAT *x, pr2eihMPC_FLOAT *B, pr2eihMPC_FLOAT *u, pr2eihMPC_FLOAT *b, pr2eihMPC_FLOAT *l, pr2eihMPC_FLOAT *r, pr2eihMPC_FLOAT *z, pr2eihMPC_FLOAT *y)
{
	int i;
	int j;
	int k = 0;
	pr2eihMPC_FLOAT AxBu[7];
	pr2eihMPC_FLOAT norm = *y;
	pr2eihMPC_FLOAT lr = 0;

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
void pr2eihMPC_LA_DENSE_DIAGZERO_MTVM2_7_14_7(pr2eihMPC_FLOAT *A, pr2eihMPC_FLOAT *x, pr2eihMPC_FLOAT *B, pr2eihMPC_FLOAT *y, pr2eihMPC_FLOAT *z)
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
void pr2eihMPC_LA_DIAGZERO_MTVM_7_7(pr2eihMPC_FLOAT *M, pr2eihMPC_FLOAT *x, pr2eihMPC_FLOAT *y)
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
void pr2eihMPC_LA_VSUBADD3_14(pr2eihMPC_FLOAT* t, pr2eihMPC_FLOAT* u, int* uidx, pr2eihMPC_FLOAT* v, pr2eihMPC_FLOAT* w, pr2eihMPC_FLOAT* y, pr2eihMPC_FLOAT* z, pr2eihMPC_FLOAT* r)
{
	int i;
	pr2eihMPC_FLOAT norm = *r;
	pr2eihMPC_FLOAT vx = 0;
	pr2eihMPC_FLOAT x;
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
void pr2eihMPC_LA_VSUBADD2_14(pr2eihMPC_FLOAT* t, int* tidx, pr2eihMPC_FLOAT* u, pr2eihMPC_FLOAT* v, pr2eihMPC_FLOAT* w, pr2eihMPC_FLOAT* y, pr2eihMPC_FLOAT* z, pr2eihMPC_FLOAT* r)
{
	int i;
	pr2eihMPC_FLOAT norm = *r;
	pr2eihMPC_FLOAT vx = 0;
	pr2eihMPC_FLOAT x;
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
void pr2eihMPC_LA_VSUBADD3_7(pr2eihMPC_FLOAT* t, pr2eihMPC_FLOAT* u, int* uidx, pr2eihMPC_FLOAT* v, pr2eihMPC_FLOAT* w, pr2eihMPC_FLOAT* y, pr2eihMPC_FLOAT* z, pr2eihMPC_FLOAT* r)
{
	int i;
	pr2eihMPC_FLOAT norm = *r;
	pr2eihMPC_FLOAT vx = 0;
	pr2eihMPC_FLOAT x;
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
void pr2eihMPC_LA_VSUBADD2_7(pr2eihMPC_FLOAT* t, int* tidx, pr2eihMPC_FLOAT* u, pr2eihMPC_FLOAT* v, pr2eihMPC_FLOAT* w, pr2eihMPC_FLOAT* y, pr2eihMPC_FLOAT* z, pr2eihMPC_FLOAT* r)
{
	int i;
	pr2eihMPC_FLOAT norm = *r;
	pr2eihMPC_FLOAT vx = 0;
	pr2eihMPC_FLOAT x;
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
void pr2eihMPC_LA_INEQ_B_GRAD_14_14_14(pr2eihMPC_FLOAT *lu, pr2eihMPC_FLOAT *su, pr2eihMPC_FLOAT *ru, pr2eihMPC_FLOAT *ll, pr2eihMPC_FLOAT *sl, pr2eihMPC_FLOAT *rl, int* lbIdx, int* ubIdx, pr2eihMPC_FLOAT *grad, pr2eihMPC_FLOAT *lubysu, pr2eihMPC_FLOAT *llbysl)
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
void pr2eihMPC_LA_INEQ_B_GRAD_7_7_7(pr2eihMPC_FLOAT *lu, pr2eihMPC_FLOAT *su, pr2eihMPC_FLOAT *ru, pr2eihMPC_FLOAT *ll, pr2eihMPC_FLOAT *sl, pr2eihMPC_FLOAT *rl, int* lbIdx, int* ubIdx, pr2eihMPC_FLOAT *grad, pr2eihMPC_FLOAT *lubysu, pr2eihMPC_FLOAT *llbysl)
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
 * of length 21.
 */
void pr2eihMPC_LA_VVADD3_21(pr2eihMPC_FLOAT *u, pr2eihMPC_FLOAT *v, pr2eihMPC_FLOAT *w, pr2eihMPC_FLOAT *z)
{
	int i;
	for( i=0; i<21; i++ ){
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
void pr2eihMPC_LA_DIAG_CHOL_ONELOOP_LBUB_14_14_14(pr2eihMPC_FLOAT *H, pr2eihMPC_FLOAT *llbysl, int* lbIdx, pr2eihMPC_FLOAT *lubysu, int* ubIdx, pr2eihMPC_FLOAT *Phi)


{
	int i;
	
	/* compute cholesky */
	for( i=0; i<14; i++ ){
		Phi[i] = H[i] + llbysl[i] + lubysu[i];

#if pr2eihMPC_SET_PRINTLEVEL > 0 && defined PRINTNUMERICALWARNINGS
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
void pr2eihMPC_LA_DIAG_MATRIXFORWARDSUB_7_14(pr2eihMPC_FLOAT *L, pr2eihMPC_FLOAT *B, pr2eihMPC_FLOAT *A)
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
void pr2eihMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_7_14(pr2eihMPC_FLOAT *L, pr2eihMPC_FLOAT *B, pr2eihMPC_FLOAT *A)
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
void pr2eihMPC_LA_DENSE_DIAGZERO_MMTM_7_14_7(pr2eihMPC_FLOAT *A, pr2eihMPC_FLOAT *B, pr2eihMPC_FLOAT *C)
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
void pr2eihMPC_LA_DIAG_FORWARDSUB_14(pr2eihMPC_FLOAT *L, pr2eihMPC_FLOAT *b, pr2eihMPC_FLOAT *y)
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
void pr2eihMPC_LA_DIAG_CHOL_ONELOOP_LBUB_7_7_7(pr2eihMPC_FLOAT *H, pr2eihMPC_FLOAT *llbysl, int* lbIdx, pr2eihMPC_FLOAT *lubysu, int* ubIdx, pr2eihMPC_FLOAT *Phi)


{
	int i;
	
	/* compute cholesky */
	for( i=0; i<7; i++ ){
		Phi[i] = H[i] + llbysl[i] + lubysu[i];

#if pr2eihMPC_SET_PRINTLEVEL > 0 && defined PRINTNUMERICALWARNINGS
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
void pr2eihMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_7_7(pr2eihMPC_FLOAT *L, pr2eihMPC_FLOAT *B, pr2eihMPC_FLOAT *A)
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
void pr2eihMPC_LA_DIAG_FORWARDSUB_7(pr2eihMPC_FLOAT *L, pr2eihMPC_FLOAT *b, pr2eihMPC_FLOAT *y)
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
void pr2eihMPC_LA_DIAGZERO_MMT_7(pr2eihMPC_FLOAT *B, pr2eihMPC_FLOAT *L)
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
void pr2eihMPC_LA_DIAGZERO_MVMSUB7_7(pr2eihMPC_FLOAT *B, pr2eihMPC_FLOAT *u, pr2eihMPC_FLOAT *b, pr2eihMPC_FLOAT *r)
{
	int i;

	for( i=0; i<7; i++ ){
		r[i] = b[i] - B[i]*u[i];
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
void pr2eihMPC_LA_DENSE_DIAGZERO_MMT2_7_14_7(pr2eihMPC_FLOAT *A, pr2eihMPC_FLOAT *B, pr2eihMPC_FLOAT *L)
{
    int i, j, k, ii, di;
    pr2eihMPC_FLOAT ltemp;
    
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
void pr2eihMPC_LA_DENSE_DIAGZERO_2MVMSUB2_7_14_7(pr2eihMPC_FLOAT *A, pr2eihMPC_FLOAT *x, pr2eihMPC_FLOAT *B, pr2eihMPC_FLOAT *u, pr2eihMPC_FLOAT *b, pr2eihMPC_FLOAT *r)
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
void pr2eihMPC_LA_DENSE_CHOL_7(pr2eihMPC_FLOAT *A, pr2eihMPC_FLOAT *L)
{
    int i, j, k, di, dj;
	 int ii, jj;

    pr2eihMPC_FLOAT l;
    pr2eihMPC_FLOAT Mii;

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
        
#if pr2eihMPC_SET_PRINTLEVEL > 0 && defined PRINTNUMERICALWARNINGS
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
void pr2eihMPC_LA_DENSE_FORWARDSUB_7(pr2eihMPC_FLOAT *L, pr2eihMPC_FLOAT *b, pr2eihMPC_FLOAT *y)
{
    int i,j,ii,di;
    pr2eihMPC_FLOAT yel;
            
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
void pr2eihMPC_LA_DENSE_MATRIXTFORWARDSUB_7_7(pr2eihMPC_FLOAT *L, pr2eihMPC_FLOAT *B, pr2eihMPC_FLOAT *A)
{
    int i,j,k,ii,di;
    pr2eihMPC_FLOAT a;
    
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
void pr2eihMPC_LA_DENSE_MMTSUB_7_7(pr2eihMPC_FLOAT *A, pr2eihMPC_FLOAT *L)
{
    int i, j, k, ii, di;
    pr2eihMPC_FLOAT ltemp;
    
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
void pr2eihMPC_LA_DENSE_MVMSUB1_7_7(pr2eihMPC_FLOAT *A, pr2eihMPC_FLOAT *x, pr2eihMPC_FLOAT *b, pr2eihMPC_FLOAT *r)
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
void pr2eihMPC_LA_DENSE_BACKWARDSUB_7(pr2eihMPC_FLOAT *L, pr2eihMPC_FLOAT *y, pr2eihMPC_FLOAT *x)
{
    int i, ii, di, j, jj, dj;
    pr2eihMPC_FLOAT xel;    
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
void pr2eihMPC_LA_DENSE_MTVMSUB_7_7(pr2eihMPC_FLOAT *A, pr2eihMPC_FLOAT *x, pr2eihMPC_FLOAT *b, pr2eihMPC_FLOAT *r)
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
 * Vector subtraction z = -x - y for vectors of length 21.
 */
void pr2eihMPC_LA_VSUB2_21(pr2eihMPC_FLOAT *x, pr2eihMPC_FLOAT *y, pr2eihMPC_FLOAT *z)
{
	int i;
	for( i=0; i<21; i++){
		z[i] = -x[i] - y[i];
	}
}


/**
 * Forward-Backward-Substitution to solve L*L^T*x = b where L is a
 * diagonal matrix of size 14 in vector
 * storage format.
 */
void pr2eihMPC_LA_DIAG_FORWARDBACKWARDSUB_14(pr2eihMPC_FLOAT *L, pr2eihMPC_FLOAT *b, pr2eihMPC_FLOAT *x)
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
void pr2eihMPC_LA_DIAG_FORWARDBACKWARDSUB_7(pr2eihMPC_FLOAT *L, pr2eihMPC_FLOAT *b, pr2eihMPC_FLOAT *x)
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
void pr2eihMPC_LA_VSUB_INDEXED_14(pr2eihMPC_FLOAT *x, int* xidx, pr2eihMPC_FLOAT *y, pr2eihMPC_FLOAT *z)
{
	int i;
	for( i=0; i<14; i++){
		z[i] = x[xidx[i]] - y[i];
	}
}


/*
 * Vector subtraction x = -u.*v - w for vectors of length 14.
 */
void pr2eihMPC_LA_VSUB3_14(pr2eihMPC_FLOAT *u, pr2eihMPC_FLOAT *v, pr2eihMPC_FLOAT *w, pr2eihMPC_FLOAT *x)
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
void pr2eihMPC_LA_VSUB2_INDEXED_14(pr2eihMPC_FLOAT *x, pr2eihMPC_FLOAT *y, int* yidx, pr2eihMPC_FLOAT *z)
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
void pr2eihMPC_LA_VSUB_INDEXED_7(pr2eihMPC_FLOAT *x, int* xidx, pr2eihMPC_FLOAT *y, pr2eihMPC_FLOAT *z)
{
	int i;
	for( i=0; i<7; i++){
		z[i] = x[xidx[i]] - y[i];
	}
}


/*
 * Vector subtraction x = -u.*v - w for vectors of length 7.
 */
void pr2eihMPC_LA_VSUB3_7(pr2eihMPC_FLOAT *u, pr2eihMPC_FLOAT *v, pr2eihMPC_FLOAT *w, pr2eihMPC_FLOAT *x)
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
void pr2eihMPC_LA_VSUB2_INDEXED_7(pr2eihMPC_FLOAT *x, pr2eihMPC_FLOAT *y, int* yidx, pr2eihMPC_FLOAT *z)
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
 * pr2eihMPC_NOPROGRESS (should be negative).
 */
int pr2eihMPC_LINESEARCH_BACKTRACKING_AFFINE(pr2eihMPC_FLOAT *l, pr2eihMPC_FLOAT *s, pr2eihMPC_FLOAT *dl, pr2eihMPC_FLOAT *ds, pr2eihMPC_FLOAT *a, pr2eihMPC_FLOAT *mu_aff)
{
    int i;
	int lsIt=1;    
    pr2eihMPC_FLOAT dltemp;
    pr2eihMPC_FLOAT dstemp;
    pr2eihMPC_FLOAT mya = 1.0;
    pr2eihMPC_FLOAT mymu;
        
    while( 1 ){                        

        /* 
         * Compute both snew and wnew together.
         * We compute also mu_affine along the way here, as the
         * values might be in registers, so it should be cheaper.
         */
        mymu = 0;
        for( i=0; i<42; i++ ){
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
        if( i == 42 ){
            break;
        } else {
            mya *= pr2eihMPC_SET_LS_SCALE_AFF;
            if( mya < pr2eihMPC_SET_LS_MINSTEP ){
                return pr2eihMPC_NOPROGRESS;
            }
        }
    }
    
    /* return new values and iteration counter */
    *a = mya;
    *mu_aff = mymu / (pr2eihMPC_FLOAT)42;
    return lsIt;
}


/*
 * Vector subtraction x = (u.*v - mu)*sigma where a is a scalar
*  and x,u,v are vectors of length 42.
 */
void pr2eihMPC_LA_VSUB5_42(pr2eihMPC_FLOAT *u, pr2eihMPC_FLOAT *v, pr2eihMPC_FLOAT mu,  pr2eihMPC_FLOAT sigma, pr2eihMPC_FLOAT *x)
{
	int i;
	for( i=0; i<42; i++){
		x[i] = u[i]*v[i] - mu;
		x[i] *= sigma;
	}
}


/*
 * Computes x=0; x(uidx) += u/su; x(vidx) -= v/sv where x is of length 14,
 * u, su, uidx are of length 14 and v, sv, vidx are of length 14.
 */
void pr2eihMPC_LA_VSUB6_INDEXED_14_14_14(pr2eihMPC_FLOAT *u, pr2eihMPC_FLOAT *su, int* uidx, pr2eihMPC_FLOAT *v, pr2eihMPC_FLOAT *sv, int* vidx, pr2eihMPC_FLOAT *x)
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
 * Computes x=0; x(uidx) += u/su; x(vidx) -= v/sv where x is of length 7,
 * u, su, uidx are of length 7 and v, sv, vidx are of length 7.
 */
void pr2eihMPC_LA_VSUB6_INDEXED_7_7_7(pr2eihMPC_FLOAT *u, pr2eihMPC_FLOAT *su, int* uidx, pr2eihMPC_FLOAT *v, pr2eihMPC_FLOAT *sv, int* vidx, pr2eihMPC_FLOAT *x)
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
 * Computes r =  B*u
 * where B is stored in diagzero format
 */
void pr2eihMPC_LA_DIAGZERO_MVM_7(pr2eihMPC_FLOAT *B, pr2eihMPC_FLOAT *u, pr2eihMPC_FLOAT *r)
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
void pr2eihMPC_LA_DENSE_DIAGZERO_2MVMADD_7_14_7(pr2eihMPC_FLOAT *A, pr2eihMPC_FLOAT *x, pr2eihMPC_FLOAT *B, pr2eihMPC_FLOAT *u, pr2eihMPC_FLOAT *r)
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
 * Vector subtraction z = x - y for vectors of length 21.
 */
void pr2eihMPC_LA_VSUB_21(pr2eihMPC_FLOAT *x, pr2eihMPC_FLOAT *y, pr2eihMPC_FLOAT *z)
{
	int i;
	for( i=0; i<21; i++){
		z[i] = x[i] - y[i];
	}
}


/** 
 * Computes z = -r./s - u.*y(y)
 * where all vectors except of y are of length 14 (length of y >= 14).
 */
void pr2eihMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_14(pr2eihMPC_FLOAT *r, pr2eihMPC_FLOAT *s, pr2eihMPC_FLOAT *u, pr2eihMPC_FLOAT *y, int* yidx, pr2eihMPC_FLOAT *z)
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
void pr2eihMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_14(pr2eihMPC_FLOAT *r, pr2eihMPC_FLOAT *s, pr2eihMPC_FLOAT *u, pr2eihMPC_FLOAT *y, int* yidx, pr2eihMPC_FLOAT *z)
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
void pr2eihMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_7(pr2eihMPC_FLOAT *r, pr2eihMPC_FLOAT *s, pr2eihMPC_FLOAT *u, pr2eihMPC_FLOAT *y, int* yidx, pr2eihMPC_FLOAT *z)
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
void pr2eihMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_7(pr2eihMPC_FLOAT *r, pr2eihMPC_FLOAT *s, pr2eihMPC_FLOAT *u, pr2eihMPC_FLOAT *y, int* yidx, pr2eihMPC_FLOAT *z)
{
	int i;
	for( i=0; i<7; i++ ){
		z[i] = -r[i]/s[i] + u[i]*y[yidx[i]];
	}
}


/*
 * Computes ds = -l.\(r + s.*dl) for vectors of length 42.
 */
void pr2eihMPC_LA_VSUB7_42(pr2eihMPC_FLOAT *l, pr2eihMPC_FLOAT *r, pr2eihMPC_FLOAT *s, pr2eihMPC_FLOAT *dl, pr2eihMPC_FLOAT *ds)
{
	int i;
	for( i=0; i<42; i++){
		ds[i] = -(r[i] + s[i]*dl[i])/l[i];
	}
}


/*
 * Vector addition x = x + y for vectors of length 21.
 */
void pr2eihMPC_LA_VADD_21(pr2eihMPC_FLOAT *x, pr2eihMPC_FLOAT *y)
{
	int i;
	for( i=0; i<21; i++){
		x[i] += y[i];
	}
}


/*
 * Vector addition x = x + y for vectors of length 14.
 */
void pr2eihMPC_LA_VADD_14(pr2eihMPC_FLOAT *x, pr2eihMPC_FLOAT *y)
{
	int i;
	for( i=0; i<14; i++){
		x[i] += y[i];
	}
}


/*
 * Vector addition x = x + y for vectors of length 42.
 */
void pr2eihMPC_LA_VADD_42(pr2eihMPC_FLOAT *x, pr2eihMPC_FLOAT *y)
{
	int i;
	for( i=0; i<42; i++){
		x[i] += y[i];
	}
}


/**
 * Backtracking line search for combined predictor/corrector step.
 * Update on variables with safety factor gamma (to keep us away from
 * boundary).
 */
int pr2eihMPC_LINESEARCH_BACKTRACKING_COMBINED(pr2eihMPC_FLOAT *z, pr2eihMPC_FLOAT *v, pr2eihMPC_FLOAT *l, pr2eihMPC_FLOAT *s, pr2eihMPC_FLOAT *dz, pr2eihMPC_FLOAT *dv, pr2eihMPC_FLOAT *dl, pr2eihMPC_FLOAT *ds, pr2eihMPC_FLOAT *a, pr2eihMPC_FLOAT *mu)
{
    int i, lsIt=1;       
    pr2eihMPC_FLOAT dltemp;
    pr2eihMPC_FLOAT dstemp;    
    pr2eihMPC_FLOAT a_gamma;
            
    *a = 1.0;
    while( 1 ){                        

        /* check whether search criterion is fulfilled */
        for( i=0; i<42; i++ ){
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
        if( i == 42 ){
            break;
        } else {
            *a *= pr2eihMPC_SET_LS_SCALE;
            if( *a < pr2eihMPC_SET_LS_MINSTEP ){
                return pr2eihMPC_NOPROGRESS;
            }
        }
    }
    
    /* update variables with safety margin */
    a_gamma = (*a)*pr2eihMPC_SET_LS_MAXSTEP;
    
    /* primal variables */
    for( i=0; i<21; i++ ){
        z[i] += a_gamma*dz[i];
    }
    
    /* equality constraint multipliers */
    for( i=0; i<14; i++ ){
        v[i] += a_gamma*dv[i];
    }
    
    /* inequality constraint multipliers & slacks, also update mu */
    *mu = 0;
    for( i=0; i<42; i++ ){
        dltemp = l[i] + a_gamma*dl[i]; l[i] = dltemp;
        dstemp = s[i] + a_gamma*ds[i]; s[i] = dstemp;
        *mu += dltemp*dstemp;
    }
    
    *a = a_gamma;
    *mu /= (pr2eihMPC_FLOAT)42;
    return lsIt;
}




/* VARIABLE DEFINITIONS ------------------------------------------------ */
pr2eihMPC_FLOAT pr2eihMPC_z[21];
pr2eihMPC_FLOAT pr2eihMPC_v[14];
pr2eihMPC_FLOAT pr2eihMPC_dz_aff[21];
pr2eihMPC_FLOAT pr2eihMPC_dv_aff[14];
pr2eihMPC_FLOAT pr2eihMPC_grad_cost[21];
pr2eihMPC_FLOAT pr2eihMPC_grad_eq[21];
pr2eihMPC_FLOAT pr2eihMPC_rd[21];
pr2eihMPC_FLOAT pr2eihMPC_l[42];
pr2eihMPC_FLOAT pr2eihMPC_s[42];
pr2eihMPC_FLOAT pr2eihMPC_lbys[42];
pr2eihMPC_FLOAT pr2eihMPC_dl_aff[42];
pr2eihMPC_FLOAT pr2eihMPC_ds_aff[42];
pr2eihMPC_FLOAT pr2eihMPC_dz_cc[21];
pr2eihMPC_FLOAT pr2eihMPC_dv_cc[14];
pr2eihMPC_FLOAT pr2eihMPC_dl_cc[42];
pr2eihMPC_FLOAT pr2eihMPC_ds_cc[42];
pr2eihMPC_FLOAT pr2eihMPC_ccrhs[42];
pr2eihMPC_FLOAT pr2eihMPC_grad_ineq[21];
pr2eihMPC_FLOAT* pr2eihMPC_z0 = pr2eihMPC_z + 0;
pr2eihMPC_FLOAT* pr2eihMPC_dzaff0 = pr2eihMPC_dz_aff + 0;
pr2eihMPC_FLOAT* pr2eihMPC_dzcc0 = pr2eihMPC_dz_cc + 0;
pr2eihMPC_FLOAT* pr2eihMPC_rd0 = pr2eihMPC_rd + 0;
pr2eihMPC_FLOAT pr2eihMPC_Lbyrd0[14];
pr2eihMPC_FLOAT* pr2eihMPC_grad_cost0 = pr2eihMPC_grad_cost + 0;
pr2eihMPC_FLOAT* pr2eihMPC_grad_eq0 = pr2eihMPC_grad_eq + 0;
pr2eihMPC_FLOAT* pr2eihMPC_grad_ineq0 = pr2eihMPC_grad_ineq + 0;
pr2eihMPC_FLOAT pr2eihMPC_ctv0[14];
pr2eihMPC_FLOAT pr2eihMPC_C0[98] = {1.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 
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
pr2eihMPC_FLOAT* pr2eihMPC_v0 = pr2eihMPC_v + 0;
pr2eihMPC_FLOAT pr2eihMPC_re0[7];
pr2eihMPC_FLOAT pr2eihMPC_beta0[7];
pr2eihMPC_FLOAT pr2eihMPC_betacc0[7];
pr2eihMPC_FLOAT* pr2eihMPC_dvaff0 = pr2eihMPC_dv_aff + 0;
pr2eihMPC_FLOAT* pr2eihMPC_dvcc0 = pr2eihMPC_dv_cc + 0;
pr2eihMPC_FLOAT pr2eihMPC_V0[98];
pr2eihMPC_FLOAT pr2eihMPC_Yd0[28];
pr2eihMPC_FLOAT pr2eihMPC_Ld0[28];
pr2eihMPC_FLOAT pr2eihMPC_yy0[7];
pr2eihMPC_FLOAT pr2eihMPC_bmy0[7];
int pr2eihMPC_lbIdx0[14] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13};
pr2eihMPC_FLOAT* pr2eihMPC_llb0 = pr2eihMPC_l + 0;
pr2eihMPC_FLOAT* pr2eihMPC_slb0 = pr2eihMPC_s + 0;
pr2eihMPC_FLOAT* pr2eihMPC_llbbyslb0 = pr2eihMPC_lbys + 0;
pr2eihMPC_FLOAT pr2eihMPC_rilb0[14];
pr2eihMPC_FLOAT* pr2eihMPC_dllbaff0 = pr2eihMPC_dl_aff + 0;
pr2eihMPC_FLOAT* pr2eihMPC_dslbaff0 = pr2eihMPC_ds_aff + 0;
pr2eihMPC_FLOAT* pr2eihMPC_dllbcc0 = pr2eihMPC_dl_cc + 0;
pr2eihMPC_FLOAT* pr2eihMPC_dslbcc0 = pr2eihMPC_ds_cc + 0;
pr2eihMPC_FLOAT* pr2eihMPC_ccrhsl0 = pr2eihMPC_ccrhs + 0;
int pr2eihMPC_ubIdx0[14] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13};
pr2eihMPC_FLOAT* pr2eihMPC_lub0 = pr2eihMPC_l + 14;
pr2eihMPC_FLOAT* pr2eihMPC_sub0 = pr2eihMPC_s + 14;
pr2eihMPC_FLOAT* pr2eihMPC_lubbysub0 = pr2eihMPC_lbys + 14;
pr2eihMPC_FLOAT pr2eihMPC_riub0[14];
pr2eihMPC_FLOAT* pr2eihMPC_dlubaff0 = pr2eihMPC_dl_aff + 14;
pr2eihMPC_FLOAT* pr2eihMPC_dsubaff0 = pr2eihMPC_ds_aff + 14;
pr2eihMPC_FLOAT* pr2eihMPC_dlubcc0 = pr2eihMPC_dl_cc + 14;
pr2eihMPC_FLOAT* pr2eihMPC_dsubcc0 = pr2eihMPC_ds_cc + 14;
pr2eihMPC_FLOAT* pr2eihMPC_ccrhsub0 = pr2eihMPC_ccrhs + 14;
pr2eihMPC_FLOAT pr2eihMPC_Phi0[14];
pr2eihMPC_FLOAT pr2eihMPC_D0[14] = {1.0000000000000000E+000, 
1.0000000000000000E+000, 
1.0000000000000000E+000, 
1.0000000000000000E+000, 
1.0000000000000000E+000, 
1.0000000000000000E+000, 
1.0000000000000000E+000};
pr2eihMPC_FLOAT pr2eihMPC_W0[14];
pr2eihMPC_FLOAT* pr2eihMPC_z1 = pr2eihMPC_z + 14;
pr2eihMPC_FLOAT* pr2eihMPC_dzaff1 = pr2eihMPC_dz_aff + 14;
pr2eihMPC_FLOAT* pr2eihMPC_dzcc1 = pr2eihMPC_dz_cc + 14;
pr2eihMPC_FLOAT* pr2eihMPC_rd1 = pr2eihMPC_rd + 14;
pr2eihMPC_FLOAT pr2eihMPC_Lbyrd1[7];
pr2eihMPC_FLOAT* pr2eihMPC_grad_cost1 = pr2eihMPC_grad_cost + 14;
pr2eihMPC_FLOAT* pr2eihMPC_grad_eq1 = pr2eihMPC_grad_eq + 14;
pr2eihMPC_FLOAT* pr2eihMPC_grad_ineq1 = pr2eihMPC_grad_ineq + 14;
pr2eihMPC_FLOAT pr2eihMPC_ctv1[7];
pr2eihMPC_FLOAT* pr2eihMPC_v1 = pr2eihMPC_v + 7;
pr2eihMPC_FLOAT pr2eihMPC_re1[7];
pr2eihMPC_FLOAT pr2eihMPC_beta1[7];
pr2eihMPC_FLOAT pr2eihMPC_betacc1[7];
pr2eihMPC_FLOAT* pr2eihMPC_dvaff1 = pr2eihMPC_dv_aff + 7;
pr2eihMPC_FLOAT* pr2eihMPC_dvcc1 = pr2eihMPC_dv_cc + 7;
pr2eihMPC_FLOAT pr2eihMPC_V1[49];
pr2eihMPC_FLOAT pr2eihMPC_Yd1[28];
pr2eihMPC_FLOAT pr2eihMPC_Ld1[28];
pr2eihMPC_FLOAT pr2eihMPC_yy1[7];
pr2eihMPC_FLOAT pr2eihMPC_bmy1[7];
pr2eihMPC_FLOAT pr2eihMPC_c1[7] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
int pr2eihMPC_lbIdx1[7] = {0, 1, 2, 3, 4, 5, 6};
pr2eihMPC_FLOAT* pr2eihMPC_llb1 = pr2eihMPC_l + 28;
pr2eihMPC_FLOAT* pr2eihMPC_slb1 = pr2eihMPC_s + 28;
pr2eihMPC_FLOAT* pr2eihMPC_llbbyslb1 = pr2eihMPC_lbys + 28;
pr2eihMPC_FLOAT pr2eihMPC_rilb1[7];
pr2eihMPC_FLOAT* pr2eihMPC_dllbaff1 = pr2eihMPC_dl_aff + 28;
pr2eihMPC_FLOAT* pr2eihMPC_dslbaff1 = pr2eihMPC_ds_aff + 28;
pr2eihMPC_FLOAT* pr2eihMPC_dllbcc1 = pr2eihMPC_dl_cc + 28;
pr2eihMPC_FLOAT* pr2eihMPC_dslbcc1 = pr2eihMPC_ds_cc + 28;
pr2eihMPC_FLOAT* pr2eihMPC_ccrhsl1 = pr2eihMPC_ccrhs + 28;
int pr2eihMPC_ubIdx1[7] = {0, 1, 2, 3, 4, 5, 6};
pr2eihMPC_FLOAT* pr2eihMPC_lub1 = pr2eihMPC_l + 35;
pr2eihMPC_FLOAT* pr2eihMPC_sub1 = pr2eihMPC_s + 35;
pr2eihMPC_FLOAT* pr2eihMPC_lubbysub1 = pr2eihMPC_lbys + 35;
pr2eihMPC_FLOAT pr2eihMPC_riub1[7];
pr2eihMPC_FLOAT* pr2eihMPC_dlubaff1 = pr2eihMPC_dl_aff + 35;
pr2eihMPC_FLOAT* pr2eihMPC_dsubaff1 = pr2eihMPC_ds_aff + 35;
pr2eihMPC_FLOAT* pr2eihMPC_dlubcc1 = pr2eihMPC_dl_cc + 35;
pr2eihMPC_FLOAT* pr2eihMPC_dsubcc1 = pr2eihMPC_ds_cc + 35;
pr2eihMPC_FLOAT* pr2eihMPC_ccrhsub1 = pr2eihMPC_ccrhs + 35;
pr2eihMPC_FLOAT pr2eihMPC_Phi1[7];
pr2eihMPC_FLOAT pr2eihMPC_D1[7] = {-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000};
pr2eihMPC_FLOAT pr2eihMPC_W1[7];
pr2eihMPC_FLOAT pr2eihMPC_Ysd1[49];
pr2eihMPC_FLOAT pr2eihMPC_Lsd1[49];
pr2eihMPC_FLOAT musigma;
pr2eihMPC_FLOAT sigma_3rdroot;
pr2eihMPC_FLOAT pr2eihMPC_Diag1_0[14];
pr2eihMPC_FLOAT pr2eihMPC_Diag2_0[14];
pr2eihMPC_FLOAT pr2eihMPC_L_0[91];




/* SOLVER CODE --------------------------------------------------------- */
int pr2eihMPC_solve(pr2eihMPC_params* params, pr2eihMPC_output* output, pr2eihMPC_info* info)
{	
int exitcode;

#if pr2eihMPC_SET_TIMING == 1
	pr2eihMPC_timer solvertimer;
	pr2eihMPC_tic(&solvertimer);
#endif
/* FUNCTION CALLS INTO LA LIBRARY -------------------------------------- */
info->it = 0;
pr2eihMPC_LA_INITIALIZEVECTOR_21(pr2eihMPC_z, 0);
pr2eihMPC_LA_INITIALIZEVECTOR_14(pr2eihMPC_v, 1);
pr2eihMPC_LA_INITIALIZEVECTOR_42(pr2eihMPC_l, 10);
pr2eihMPC_LA_INITIALIZEVECTOR_42(pr2eihMPC_s, 10);
info->mu = 0;
pr2eihMPC_LA_DOTACC_42(pr2eihMPC_l, pr2eihMPC_s, &info->mu);
info->mu /= 42;
while( 1 ){
info->pobj = 0;
pr2eihMPC_LA_DIAG_QUADFCN_14(params->H1, params->f1, pr2eihMPC_z0, pr2eihMPC_grad_cost0, &info->pobj);
pr2eihMPC_LA_DIAG_QUADFCN_7(params->H2, params->f2, pr2eihMPC_z1, pr2eihMPC_grad_cost1, &info->pobj);
info->res_eq = 0;
info->dgap = 0;
pr2eihMPC_LA_DIAGZERO_MVMSUB6_7(pr2eihMPC_D0, pr2eihMPC_z0, params->c1, pr2eihMPC_v0, pr2eihMPC_re0, &info->dgap, &info->res_eq);
pr2eihMPC_LA_DENSE_DIAGZERO_MVMSUB3_7_14_7(pr2eihMPC_C0, pr2eihMPC_z0, pr2eihMPC_D1, pr2eihMPC_z1, pr2eihMPC_c1, pr2eihMPC_v1, pr2eihMPC_re1, &info->dgap, &info->res_eq);
pr2eihMPC_LA_DENSE_DIAGZERO_MTVM2_7_14_7(pr2eihMPC_C0, pr2eihMPC_v1, pr2eihMPC_D0, pr2eihMPC_v0, pr2eihMPC_grad_eq0);
pr2eihMPC_LA_DIAGZERO_MTVM_7_7(pr2eihMPC_D1, pr2eihMPC_v1, pr2eihMPC_grad_eq1);
info->res_ineq = 0;
pr2eihMPC_LA_VSUBADD3_14(params->lb1, pr2eihMPC_z0, pr2eihMPC_lbIdx0, pr2eihMPC_llb0, pr2eihMPC_slb0, pr2eihMPC_rilb0, &info->dgap, &info->res_ineq);
pr2eihMPC_LA_VSUBADD2_14(pr2eihMPC_z0, pr2eihMPC_ubIdx0, params->ub1, pr2eihMPC_lub0, pr2eihMPC_sub0, pr2eihMPC_riub0, &info->dgap, &info->res_ineq);
pr2eihMPC_LA_VSUBADD3_7(params->lb2, pr2eihMPC_z1, pr2eihMPC_lbIdx1, pr2eihMPC_llb1, pr2eihMPC_slb1, pr2eihMPC_rilb1, &info->dgap, &info->res_ineq);
pr2eihMPC_LA_VSUBADD2_7(pr2eihMPC_z1, pr2eihMPC_ubIdx1, params->ub2, pr2eihMPC_lub1, pr2eihMPC_sub1, pr2eihMPC_riub1, &info->dgap, &info->res_ineq);
pr2eihMPC_LA_INEQ_B_GRAD_14_14_14(pr2eihMPC_lub0, pr2eihMPC_sub0, pr2eihMPC_riub0, pr2eihMPC_llb0, pr2eihMPC_slb0, pr2eihMPC_rilb0, pr2eihMPC_lbIdx0, pr2eihMPC_ubIdx0, pr2eihMPC_grad_ineq0, pr2eihMPC_lubbysub0, pr2eihMPC_llbbyslb0);
pr2eihMPC_LA_INEQ_B_GRAD_7_7_7(pr2eihMPC_lub1, pr2eihMPC_sub1, pr2eihMPC_riub1, pr2eihMPC_llb1, pr2eihMPC_slb1, pr2eihMPC_rilb1, pr2eihMPC_lbIdx1, pr2eihMPC_ubIdx1, pr2eihMPC_grad_ineq1, pr2eihMPC_lubbysub1, pr2eihMPC_llbbyslb1);
info->dobj = info->pobj - info->dgap;
info->rdgap = info->pobj ? info->dgap / info->pobj : 1e6;
if( info->rdgap < 0 ) info->rdgap = -info->rdgap;
if( info->mu < pr2eihMPC_SET_ACC_KKTCOMPL
    && (info->rdgap < pr2eihMPC_SET_ACC_RDGAP || info->dgap < pr2eihMPC_SET_ACC_KKTCOMPL)
    && info->res_eq < pr2eihMPC_SET_ACC_RESEQ
    && info->res_ineq < pr2eihMPC_SET_ACC_RESINEQ ){
exitcode = pr2eihMPC_OPTIMAL; break; }
if( info->it == pr2eihMPC_SET_MAXIT ){
exitcode = pr2eihMPC_MAXITREACHED; break; }
pr2eihMPC_LA_VVADD3_21(pr2eihMPC_grad_cost, pr2eihMPC_grad_eq, pr2eihMPC_grad_ineq, pr2eihMPC_rd);
pr2eihMPC_LA_DIAG_CHOL_ONELOOP_LBUB_14_14_14(params->H1, pr2eihMPC_llbbyslb0, pr2eihMPC_lbIdx0, pr2eihMPC_lubbysub0, pr2eihMPC_ubIdx0, pr2eihMPC_Phi0);
pr2eihMPC_LA_DIAG_MATRIXFORWARDSUB_7_14(pr2eihMPC_Phi0, pr2eihMPC_C0, pr2eihMPC_V0);
pr2eihMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_7_14(pr2eihMPC_Phi0, pr2eihMPC_D0, pr2eihMPC_W0);
pr2eihMPC_LA_DENSE_DIAGZERO_MMTM_7_14_7(pr2eihMPC_W0, pr2eihMPC_V0, pr2eihMPC_Ysd1);
pr2eihMPC_LA_DIAG_FORWARDSUB_14(pr2eihMPC_Phi0, pr2eihMPC_rd0, pr2eihMPC_Lbyrd0);
pr2eihMPC_LA_DIAG_CHOL_ONELOOP_LBUB_7_7_7(params->H2, pr2eihMPC_llbbyslb1, pr2eihMPC_lbIdx1, pr2eihMPC_lubbysub1, pr2eihMPC_ubIdx1, pr2eihMPC_Phi1);
pr2eihMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_7_7(pr2eihMPC_Phi1, pr2eihMPC_D1, pr2eihMPC_W1);
pr2eihMPC_LA_DIAG_FORWARDSUB_7(pr2eihMPC_Phi1, pr2eihMPC_rd1, pr2eihMPC_Lbyrd1);
pr2eihMPC_LA_DIAGZERO_MMT_7(pr2eihMPC_W0, pr2eihMPC_Yd0);
pr2eihMPC_LA_DIAGZERO_MVMSUB7_7(pr2eihMPC_W0, pr2eihMPC_Lbyrd0, pr2eihMPC_re0, pr2eihMPC_beta0);
pr2eihMPC_LA_DENSE_DIAGZERO_MMT2_7_14_7(pr2eihMPC_V0, pr2eihMPC_W1, pr2eihMPC_Yd1);
pr2eihMPC_LA_DENSE_DIAGZERO_2MVMSUB2_7_14_7(pr2eihMPC_V0, pr2eihMPC_Lbyrd0, pr2eihMPC_W1, pr2eihMPC_Lbyrd1, pr2eihMPC_re1, pr2eihMPC_beta1);
pr2eihMPC_LA_DENSE_CHOL_7(pr2eihMPC_Yd0, pr2eihMPC_Ld0);
pr2eihMPC_LA_DENSE_FORWARDSUB_7(pr2eihMPC_Ld0, pr2eihMPC_beta0, pr2eihMPC_yy0);
pr2eihMPC_LA_DENSE_MATRIXTFORWARDSUB_7_7(pr2eihMPC_Ld0, pr2eihMPC_Ysd1, pr2eihMPC_Lsd1);
pr2eihMPC_LA_DENSE_MMTSUB_7_7(pr2eihMPC_Lsd1, pr2eihMPC_Yd1);
pr2eihMPC_LA_DENSE_CHOL_7(pr2eihMPC_Yd1, pr2eihMPC_Ld1);
pr2eihMPC_LA_DENSE_MVMSUB1_7_7(pr2eihMPC_Lsd1, pr2eihMPC_yy0, pr2eihMPC_beta1, pr2eihMPC_bmy1);
pr2eihMPC_LA_DENSE_FORWARDSUB_7(pr2eihMPC_Ld1, pr2eihMPC_bmy1, pr2eihMPC_yy1);
pr2eihMPC_LA_DENSE_BACKWARDSUB_7(pr2eihMPC_Ld1, pr2eihMPC_yy1, pr2eihMPC_dvaff1);
pr2eihMPC_LA_DENSE_MTVMSUB_7_7(pr2eihMPC_Lsd1, pr2eihMPC_dvaff1, pr2eihMPC_yy0, pr2eihMPC_bmy0);
pr2eihMPC_LA_DENSE_BACKWARDSUB_7(pr2eihMPC_Ld0, pr2eihMPC_bmy0, pr2eihMPC_dvaff0);
pr2eihMPC_LA_DENSE_DIAGZERO_MTVM2_7_14_7(pr2eihMPC_C0, pr2eihMPC_dvaff1, pr2eihMPC_D0, pr2eihMPC_dvaff0, pr2eihMPC_grad_eq0);
pr2eihMPC_LA_DIAGZERO_MTVM_7_7(pr2eihMPC_D1, pr2eihMPC_dvaff1, pr2eihMPC_grad_eq1);
pr2eihMPC_LA_VSUB2_21(pr2eihMPC_rd, pr2eihMPC_grad_eq, pr2eihMPC_rd);
pr2eihMPC_LA_DIAG_FORWARDBACKWARDSUB_14(pr2eihMPC_Phi0, pr2eihMPC_rd0, pr2eihMPC_dzaff0);
pr2eihMPC_LA_DIAG_FORWARDBACKWARDSUB_7(pr2eihMPC_Phi1, pr2eihMPC_rd1, pr2eihMPC_dzaff1);
pr2eihMPC_LA_VSUB_INDEXED_14(pr2eihMPC_dzaff0, pr2eihMPC_lbIdx0, pr2eihMPC_rilb0, pr2eihMPC_dslbaff0);
pr2eihMPC_LA_VSUB3_14(pr2eihMPC_llbbyslb0, pr2eihMPC_dslbaff0, pr2eihMPC_llb0, pr2eihMPC_dllbaff0);
pr2eihMPC_LA_VSUB2_INDEXED_14(pr2eihMPC_riub0, pr2eihMPC_dzaff0, pr2eihMPC_ubIdx0, pr2eihMPC_dsubaff0);
pr2eihMPC_LA_VSUB3_14(pr2eihMPC_lubbysub0, pr2eihMPC_dsubaff0, pr2eihMPC_lub0, pr2eihMPC_dlubaff0);
pr2eihMPC_LA_VSUB_INDEXED_7(pr2eihMPC_dzaff1, pr2eihMPC_lbIdx1, pr2eihMPC_rilb1, pr2eihMPC_dslbaff1);
pr2eihMPC_LA_VSUB3_7(pr2eihMPC_llbbyslb1, pr2eihMPC_dslbaff1, pr2eihMPC_llb1, pr2eihMPC_dllbaff1);
pr2eihMPC_LA_VSUB2_INDEXED_7(pr2eihMPC_riub1, pr2eihMPC_dzaff1, pr2eihMPC_ubIdx1, pr2eihMPC_dsubaff1);
pr2eihMPC_LA_VSUB3_7(pr2eihMPC_lubbysub1, pr2eihMPC_dsubaff1, pr2eihMPC_lub1, pr2eihMPC_dlubaff1);
info->lsit_aff = pr2eihMPC_LINESEARCH_BACKTRACKING_AFFINE(pr2eihMPC_l, pr2eihMPC_s, pr2eihMPC_dl_aff, pr2eihMPC_ds_aff, &info->step_aff, &info->mu_aff);
if( info->lsit_aff == pr2eihMPC_NOPROGRESS ){
exitcode = pr2eihMPC_NOPROGRESS; break;
}
sigma_3rdroot = info->mu_aff / info->mu;
info->sigma = sigma_3rdroot*sigma_3rdroot*sigma_3rdroot;
musigma = info->mu * info->sigma;
pr2eihMPC_LA_VSUB5_42(pr2eihMPC_ds_aff, pr2eihMPC_dl_aff, info->mu, info->sigma, pr2eihMPC_ccrhs);
pr2eihMPC_LA_VSUB6_INDEXED_14_14_14(pr2eihMPC_ccrhsub0, pr2eihMPC_sub0, pr2eihMPC_ubIdx0, pr2eihMPC_ccrhsl0, pr2eihMPC_slb0, pr2eihMPC_lbIdx0, pr2eihMPC_rd0);
pr2eihMPC_LA_VSUB6_INDEXED_7_7_7(pr2eihMPC_ccrhsub1, pr2eihMPC_sub1, pr2eihMPC_ubIdx1, pr2eihMPC_ccrhsl1, pr2eihMPC_slb1, pr2eihMPC_lbIdx1, pr2eihMPC_rd1);
pr2eihMPC_LA_DIAG_FORWARDSUB_14(pr2eihMPC_Phi0, pr2eihMPC_rd0, pr2eihMPC_Lbyrd0);
pr2eihMPC_LA_DIAG_FORWARDSUB_7(pr2eihMPC_Phi1, pr2eihMPC_rd1, pr2eihMPC_Lbyrd1);
pr2eihMPC_LA_DIAGZERO_MVM_7(pr2eihMPC_W0, pr2eihMPC_Lbyrd0, pr2eihMPC_beta0);
pr2eihMPC_LA_DENSE_FORWARDSUB_7(pr2eihMPC_Ld0, pr2eihMPC_beta0, pr2eihMPC_yy0);
pr2eihMPC_LA_DENSE_DIAGZERO_2MVMADD_7_14_7(pr2eihMPC_V0, pr2eihMPC_Lbyrd0, pr2eihMPC_W1, pr2eihMPC_Lbyrd1, pr2eihMPC_beta1);
pr2eihMPC_LA_DENSE_MVMSUB1_7_7(pr2eihMPC_Lsd1, pr2eihMPC_yy0, pr2eihMPC_beta1, pr2eihMPC_bmy1);
pr2eihMPC_LA_DENSE_FORWARDSUB_7(pr2eihMPC_Ld1, pr2eihMPC_bmy1, pr2eihMPC_yy1);
pr2eihMPC_LA_DENSE_BACKWARDSUB_7(pr2eihMPC_Ld1, pr2eihMPC_yy1, pr2eihMPC_dvcc1);
pr2eihMPC_LA_DENSE_MTVMSUB_7_7(pr2eihMPC_Lsd1, pr2eihMPC_dvcc1, pr2eihMPC_yy0, pr2eihMPC_bmy0);
pr2eihMPC_LA_DENSE_BACKWARDSUB_7(pr2eihMPC_Ld0, pr2eihMPC_bmy0, pr2eihMPC_dvcc0);
pr2eihMPC_LA_DENSE_DIAGZERO_MTVM2_7_14_7(pr2eihMPC_C0, pr2eihMPC_dvcc1, pr2eihMPC_D0, pr2eihMPC_dvcc0, pr2eihMPC_grad_eq0);
pr2eihMPC_LA_DIAGZERO_MTVM_7_7(pr2eihMPC_D1, pr2eihMPC_dvcc1, pr2eihMPC_grad_eq1);
pr2eihMPC_LA_VSUB_21(pr2eihMPC_rd, pr2eihMPC_grad_eq, pr2eihMPC_rd);
pr2eihMPC_LA_DIAG_FORWARDBACKWARDSUB_14(pr2eihMPC_Phi0, pr2eihMPC_rd0, pr2eihMPC_dzcc0);
pr2eihMPC_LA_DIAG_FORWARDBACKWARDSUB_7(pr2eihMPC_Phi1, pr2eihMPC_rd1, pr2eihMPC_dzcc1);
pr2eihMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_14(pr2eihMPC_ccrhsl0, pr2eihMPC_slb0, pr2eihMPC_llbbyslb0, pr2eihMPC_dzcc0, pr2eihMPC_lbIdx0, pr2eihMPC_dllbcc0);
pr2eihMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_14(pr2eihMPC_ccrhsub0, pr2eihMPC_sub0, pr2eihMPC_lubbysub0, pr2eihMPC_dzcc0, pr2eihMPC_ubIdx0, pr2eihMPC_dlubcc0);
pr2eihMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_7(pr2eihMPC_ccrhsl1, pr2eihMPC_slb1, pr2eihMPC_llbbyslb1, pr2eihMPC_dzcc1, pr2eihMPC_lbIdx1, pr2eihMPC_dllbcc1);
pr2eihMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_7(pr2eihMPC_ccrhsub1, pr2eihMPC_sub1, pr2eihMPC_lubbysub1, pr2eihMPC_dzcc1, pr2eihMPC_ubIdx1, pr2eihMPC_dlubcc1);
pr2eihMPC_LA_VSUB7_42(pr2eihMPC_l, pr2eihMPC_ccrhs, pr2eihMPC_s, pr2eihMPC_dl_cc, pr2eihMPC_ds_cc);
pr2eihMPC_LA_VADD_21(pr2eihMPC_dz_cc, pr2eihMPC_dz_aff);
pr2eihMPC_LA_VADD_14(pr2eihMPC_dv_cc, pr2eihMPC_dv_aff);
pr2eihMPC_LA_VADD_42(pr2eihMPC_dl_cc, pr2eihMPC_dl_aff);
pr2eihMPC_LA_VADD_42(pr2eihMPC_ds_cc, pr2eihMPC_ds_aff);
info->lsit_cc = pr2eihMPC_LINESEARCH_BACKTRACKING_COMBINED(pr2eihMPC_z, pr2eihMPC_v, pr2eihMPC_l, pr2eihMPC_s, pr2eihMPC_dz_cc, pr2eihMPC_dv_cc, pr2eihMPC_dl_cc, pr2eihMPC_ds_cc, &info->step_cc, &info->mu);
if( info->lsit_cc == pr2eihMPC_NOPROGRESS ){
exitcode = pr2eihMPC_NOPROGRESS; break;
}
info->it++;
}
output->z1[0] = pr2eihMPC_z0[0];
output->z1[1] = pr2eihMPC_z0[1];
output->z1[2] = pr2eihMPC_z0[2];
output->z1[3] = pr2eihMPC_z0[3];
output->z1[4] = pr2eihMPC_z0[4];
output->z1[5] = pr2eihMPC_z0[5];
output->z1[6] = pr2eihMPC_z0[6];
output->z1[7] = pr2eihMPC_z0[7];
output->z1[8] = pr2eihMPC_z0[8];
output->z1[9] = pr2eihMPC_z0[9];
output->z1[10] = pr2eihMPC_z0[10];
output->z1[11] = pr2eihMPC_z0[11];
output->z1[12] = pr2eihMPC_z0[12];
output->z1[13] = pr2eihMPC_z0[13];
output->z2[0] = pr2eihMPC_z1[0];
output->z2[1] = pr2eihMPC_z1[1];
output->z2[2] = pr2eihMPC_z1[2];
output->z2[3] = pr2eihMPC_z1[3];
output->z2[4] = pr2eihMPC_z1[4];
output->z2[5] = pr2eihMPC_z1[5];
output->z2[6] = pr2eihMPC_z1[6];

#if pr2eihMPC_SET_TIMING == 1
info->solvetime = pr2eihMPC_toc(&solvertimer);
#if pr2eihMPC_SET_PRINTLEVEL > 0 && pr2eihMPC_SET_TIMING == 1
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
