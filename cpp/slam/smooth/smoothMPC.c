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

#include "smoothMPC.h"

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
#ifndef USEMEXPRINTS
#include <stdio.h>
#define PRINTTEXT printf
#else
#include "mex.h"
#define PRINTTEXT mexPrintf
#endif



/* LINEAR ALGEBRA LIBRARY ---------------------------------------------- */
/*
 * Initializes a vector of length 322 with a value.
 */
void smoothMPC_LA_INITIALIZEVECTOR_322(smoothMPC_FLOAT* vec, smoothMPC_FLOAT value)
{
	int i;
	for( i=0; i<322; i++ )
	{
		vec[i] = value;
	}
}


/*
 * Initializes a vector of length 90 with a value.
 */
void smoothMPC_LA_INITIALIZEVECTOR_90(smoothMPC_FLOAT* vec, smoothMPC_FLOAT value)
{
	int i;
	for( i=0; i<90; i++ )
	{
		vec[i] = value;
	}
}


/*
 * Initializes a vector of length 470 with a value.
 */
void smoothMPC_LA_INITIALIZEVECTOR_470(smoothMPC_FLOAT* vec, smoothMPC_FLOAT value)
{
	int i;
	for( i=0; i<470; i++ )
	{
		vec[i] = value;
	}
}


/* 
 * Calculates a dot product and adds it to a variable: z += x'*y; 
 * This function is for vectors of length 470.
 */
void smoothMPC_LA_DOTACC_470(smoothMPC_FLOAT *x, smoothMPC_FLOAT *y, smoothMPC_FLOAT *z)
{
	int i;
	for( i=0; i<470; i++ ){
		*z += x[i]*y[i];
	}
}


/*
 * Calculates the gradient and the value for a quadratic function 0.5*z'*H*z + f'*z
 *
 * INPUTS:     H  - Symmetric Hessian, diag matrix of size [11 x 11]
 *             f  - column vector of size 11
 *             z  - column vector of size 11
 *
 * OUTPUTS: grad  - gradient at z (= H*z + f), column vector of size 11
 *          value <-- value + 0.5*z'*H*z + f'*z (value will be modified)
 */
void smoothMPC_LA_DIAG_QUADFCN_11(smoothMPC_FLOAT* H, smoothMPC_FLOAT* f, smoothMPC_FLOAT* z, smoothMPC_FLOAT* grad, smoothMPC_FLOAT* value)
{
	int i;
	smoothMPC_FLOAT hz;	
	for( i=0; i<11; i++){
		hz = H[i]*z[i];
		grad[i] = hz + f[i];
		*value += 0.5*hz*z[i] + f[i]*z[i];
	}
}


/*
 * Calculates the gradient and the value for a quadratic function 0.5*z'*H*z + f'*z
 *
 * INPUTS:     H  - Symmetric Hessian, diag matrix of size [3 x 3]
 *             f  - column vector of size 3
 *             z  - column vector of size 3
 *
 * OUTPUTS: grad  - gradient at z (= H*z + f), column vector of size 3
 *          value <-- value + 0.5*z'*H*z + f'*z (value will be modified)
 */
void smoothMPC_LA_DIAG_QUADFCN_3(smoothMPC_FLOAT* H, smoothMPC_FLOAT* f, smoothMPC_FLOAT* z, smoothMPC_FLOAT* grad, smoothMPC_FLOAT* value)
{
	int i;
	smoothMPC_FLOAT hz;	
	for( i=0; i<3; i++){
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
void smoothMPC_LA_DIAGZERO_MVMSUB6_3(smoothMPC_FLOAT *B, smoothMPC_FLOAT *u, smoothMPC_FLOAT *b, smoothMPC_FLOAT *l, smoothMPC_FLOAT *r, smoothMPC_FLOAT *z, smoothMPC_FLOAT *y)
{
	int i;
	smoothMPC_FLOAT Bu[3];
	smoothMPC_FLOAT norm = *y;
	smoothMPC_FLOAT lr = 0;

	/* do A*x + B*u first */
	for( i=0; i<3; i++ ){
		Bu[i] = B[i]*u[i];
	}	

	for( i=0; i<3; i++ ){
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
void smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(smoothMPC_FLOAT *A, smoothMPC_FLOAT *x, smoothMPC_FLOAT *B, smoothMPC_FLOAT *u, smoothMPC_FLOAT *b, smoothMPC_FLOAT *l, smoothMPC_FLOAT *r, smoothMPC_FLOAT *z, smoothMPC_FLOAT *y)
{
	int i;
	int j;
	int k = 0;
	smoothMPC_FLOAT AxBu[3];
	smoothMPC_FLOAT norm = *y;
	smoothMPC_FLOAT lr = 0;

	/* do A*x + B*u first */
	for( i=0; i<3; i++ ){
		AxBu[i] = A[k++]*x[0] + B[i]*u[i];
	}	

	for( j=1; j<11; j++ ){		
		for( i=0; i<3; i++ ){
			AxBu[i] += A[k++]*x[j];
		}
	}

	for( i=0; i<3; i++ ){
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
void smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_3(smoothMPC_FLOAT *A, smoothMPC_FLOAT *x, smoothMPC_FLOAT *B, smoothMPC_FLOAT *u, smoothMPC_FLOAT *b, smoothMPC_FLOAT *l, smoothMPC_FLOAT *r, smoothMPC_FLOAT *z, smoothMPC_FLOAT *y)
{
	int i;
	int j;
	int k = 0;
	smoothMPC_FLOAT AxBu[3];
	smoothMPC_FLOAT norm = *y;
	smoothMPC_FLOAT lr = 0;

	/* do A*x + B*u first */
	for( i=0; i<3; i++ ){
		AxBu[i] = A[k++]*x[0] + B[i]*u[i];
	}	

	for( j=1; j<11; j++ ){		
		for( i=0; i<3; i++ ){
			AxBu[i] += A[k++]*x[j];
		}
	}

	for( i=0; i<3; i++ ){
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
 * where A is of size [3 x 11] and stored in column major format.
 * and B is of size [3 x 11] and stored in diagzero format
 * Note the transposes of A and B!
 */
void smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(smoothMPC_FLOAT *A, smoothMPC_FLOAT *x, smoothMPC_FLOAT *B, smoothMPC_FLOAT *y, smoothMPC_FLOAT *z)
{
	int i;
	int j;
	int k = 0;
	for( i=0; i<3; i++ ){
		z[i] = 0;
		for( j=0; j<3; j++ ){
			z[i] += A[k++]*x[j];
		}
		z[i] += B[i]*y[i];
	}
	for( i=3 ;i<11; i++ ){
		z[i] = 0;
		for( j=0; j<3; j++ ){
			z[i] += A[k++]*x[j];
		}
	}
}


/*
 * Matrix vector multiplication y = M'*x where M is of size [3 x 3]
 * and stored in diagzero format. Note the transpose of M!
 */
void smoothMPC_LA_DIAGZERO_MTVM_3_3(smoothMPC_FLOAT *M, smoothMPC_FLOAT *x, smoothMPC_FLOAT *y)
{
	int i;
	for( i=0; i<3; i++ ){
		y[i] = M[i]*x[i];
	}
}


/*
 * Vector subtraction and addition.
 *	 Input: five vectors t, tidx, u, v, w and two scalars z and r
 *	 Output: y = t(tidx) - u + w
 *           z = z - v'*x;
 *           r = max([norm(y,inf), z]);
 * for vectors of length 11. Output z is of course scalar.
 */
void smoothMPC_LA_VSUBADD3_11(smoothMPC_FLOAT* t, smoothMPC_FLOAT* u, int* uidx, smoothMPC_FLOAT* v, smoothMPC_FLOAT* w, smoothMPC_FLOAT* y, smoothMPC_FLOAT* z, smoothMPC_FLOAT* r)
{
	int i;
	smoothMPC_FLOAT norm = *r;
	smoothMPC_FLOAT vx = 0;
	smoothMPC_FLOAT x;
	for( i=0; i<11; i++){
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
void smoothMPC_LA_VSUBADD2_5(smoothMPC_FLOAT* t, int* tidx, smoothMPC_FLOAT* u, smoothMPC_FLOAT* v, smoothMPC_FLOAT* w, smoothMPC_FLOAT* y, smoothMPC_FLOAT* z, smoothMPC_FLOAT* r)
{
	int i;
	smoothMPC_FLOAT norm = *r;
	smoothMPC_FLOAT vx = 0;
	smoothMPC_FLOAT x;
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
 * Vector subtraction and addition.
 *	 Input: five vectors t, tidx, u, v, w and two scalars z and r
 *	 Output: y = t(tidx) - u + w
 *           z = z - v'*x;
 *           r = max([norm(y,inf), z]);
 * for vectors of length 3. Output z is of course scalar.
 */
void smoothMPC_LA_VSUBADD3_3(smoothMPC_FLOAT* t, smoothMPC_FLOAT* u, int* uidx, smoothMPC_FLOAT* v, smoothMPC_FLOAT* w, smoothMPC_FLOAT* y, smoothMPC_FLOAT* z, smoothMPC_FLOAT* r)
{
	int i;
	smoothMPC_FLOAT norm = *r;
	smoothMPC_FLOAT vx = 0;
	smoothMPC_FLOAT x;
	for( i=0; i<3; i++){
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
 * for vectors of length 3. Output z is of course scalar.
 */
void smoothMPC_LA_VSUBADD2_3(smoothMPC_FLOAT* t, int* tidx, smoothMPC_FLOAT* u, smoothMPC_FLOAT* v, smoothMPC_FLOAT* w, smoothMPC_FLOAT* y, smoothMPC_FLOAT* z, smoothMPC_FLOAT* r)
{
	int i;
	smoothMPC_FLOAT norm = *r;
	smoothMPC_FLOAT vx = 0;
	smoothMPC_FLOAT x;
	for( i=0; i<3; i++){
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
 * Special function for box constraints of length 11
 * Returns also L/S, a value that is often used elsewhere.
 */
void smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_FLOAT *lu, smoothMPC_FLOAT *su, smoothMPC_FLOAT *ru, smoothMPC_FLOAT *ll, smoothMPC_FLOAT *sl, smoothMPC_FLOAT *rl, int* lbIdx, int* ubIdx, smoothMPC_FLOAT *grad, smoothMPC_FLOAT *lubysu, smoothMPC_FLOAT *llbysl)
{
	int i;
	for( i=0; i<11; i++ ){
		grad[i] = 0;
	}
	for( i=0; i<11; i++ ){		
		llbysl[i] = ll[i] / sl[i];
		grad[lbIdx[i]] -= llbysl[i]*rl[i];
	}
	for( i=0; i<5; i++ ){
		lubysu[i] = lu[i] / su[i];
		grad[ubIdx[i]] += lubysu[i]*ru[i];
	}
}


/*
 * Computes inequality constraints gradient-
 * Special function for box constraints of length 3
 * Returns also L/S, a value that is often used elsewhere.
 */
void smoothMPC_LA_INEQ_B_GRAD_3_3_3(smoothMPC_FLOAT *lu, smoothMPC_FLOAT *su, smoothMPC_FLOAT *ru, smoothMPC_FLOAT *ll, smoothMPC_FLOAT *sl, smoothMPC_FLOAT *rl, int* lbIdx, int* ubIdx, smoothMPC_FLOAT *grad, smoothMPC_FLOAT *lubysu, smoothMPC_FLOAT *llbysl)
{
	int i;
	for( i=0; i<3; i++ ){
		grad[i] = 0;
	}
	for( i=0; i<3; i++ ){		
		llbysl[i] = ll[i] / sl[i];
		grad[lbIdx[i]] -= llbysl[i]*rl[i];
	}
	for( i=0; i<3; i++ ){
		lubysu[i] = lu[i] / su[i];
		grad[ubIdx[i]] += lubysu[i]*ru[i];
	}
}


/*
 * Addition of three vectors  z = u + w + v
 * of length 322.
 */
void smoothMPC_LA_VVADD3_322(smoothMPC_FLOAT *u, smoothMPC_FLOAT *v, smoothMPC_FLOAT *w, smoothMPC_FLOAT *z)
{
	int i;
	for( i=0; i<322; i++ ){
		z[i] = u[i] + v[i] + w[i];
	}
}


/*
 * Special function to compute the diagonal cholesky factorization of the 
 * positive definite augmented Hessian for block size 11.
 *
 * Inputs: - H = diagonal cost Hessian in diagonal storage format
 *         - llbysl = L / S of lower bounds
 *         - lubysu = L / S of upper bounds
 *
 * Output: Phi = sqrt(H + diag(llbysl) + diag(lubysu))
 * where Phi is stored in diagonal storage format
 */
void smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(smoothMPC_FLOAT *H, smoothMPC_FLOAT *llbysl, int* lbIdx, smoothMPC_FLOAT *lubysu, int* ubIdx, smoothMPC_FLOAT *Phi)


{
	int i;
	
	/* copy  H into PHI */
	for( i=0; i<11; i++ ){
		Phi[i] = H[i];
	}

	/* add llbysl onto Phi where necessary */
	for( i=0; i<11; i++ ){
		Phi[lbIdx[i]] += llbysl[i];
	}

	/* add lubysu onto Phi where necessary */
	for( i=0; i<5; i++){
		Phi[ubIdx[i]] +=  lubysu[i];
	}
	
	/* compute cholesky */
	for(i=0; i<11; i++)
	{
#if smoothMPC_SET_PRINTLEVEL > 0 && defined PRINTNUMERICALWARNINGS
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
 * where A is to be computed and is of size [3 x 11],
 * B is given and of size [3 x 11], L is a diagonal
 * matrix of size 3 stored in diagonal matrix 
 * storage format. Note the transpose of L has no impact!
 *
 * Result: A in column major storage format.
 *
 */
void smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_FLOAT *L, smoothMPC_FLOAT *B, smoothMPC_FLOAT *A)
{
    int i,j;
	 int k = 0;

	for( j=0; j<11; j++){
		for( i=0; i<3; i++){
			A[k] = B[k]/L[j];
			k++;
		}
	}

}


/**
 * Forward substitution for the matrix equation A*L' = B
 * where A is to be computed and is of size [3 x 11],
 * B is given and of size [3 x 11], L is a diagonal
 *  matrix of size 11 stored in diagonal 
 * storage format. Note the transpose of L!
 *
 * Result: A in diagonalzero storage format.
 *
 */
void smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_FLOAT *L, smoothMPC_FLOAT *B, smoothMPC_FLOAT *A)
{
	int j;
    for( j=0; j<11; j++ ){   
		A[j] = B[j]/L[j];
     }
}


/**
 * Compute C = A*B' where 
 *
 *	size(A) = [3 x 11]
 *  size(B) = [3 x 11] in diagzero format
 * 
 * A and C matrices are stored in column major format.
 * 
 * 
 */
void smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_FLOAT *A, smoothMPC_FLOAT *B, smoothMPC_FLOAT *C)
{
    int i, j;
	
	for( i=0; i<3; i++ ){
		for( j=0; j<3; j++){
			C[j*3+i] = B[i*3+j]*A[i];
		}
	}

}


/**
 * Forward substitution to solve L*y = b where L is a
 * diagonal matrix in vector storage format.
 * 
 * The dimensions involved are 11.
 */
void smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_FLOAT *L, smoothMPC_FLOAT *b, smoothMPC_FLOAT *y)
{
    int i;

    for( i=0; i<11; i++ ){
		y[i] = b[i]/L[i];
    }
}


/*
 * Special function to compute the diagonal cholesky factorization of the 
 * positive definite augmented Hessian for block size 3.
 *
 * Inputs: - H = diagonal cost Hessian in diagonal storage format
 *         - llbysl = L / S of lower bounds
 *         - lubysu = L / S of upper bounds
 *
 * Output: Phi = sqrt(H + diag(llbysl) + diag(lubysu))
 * where Phi is stored in diagonal storage format
 */
void smoothMPC_LA_DIAG_CHOL_ONELOOP_LBUB_3_3_3(smoothMPC_FLOAT *H, smoothMPC_FLOAT *llbysl, int* lbIdx, smoothMPC_FLOAT *lubysu, int* ubIdx, smoothMPC_FLOAT *Phi)


{
	int i;
	
	/* compute cholesky */
	for( i=0; i<3; i++ ){
		Phi[i] = H[i] + llbysl[i] + lubysu[i];

#if smoothMPC_SET_PRINTLEVEL > 0 && defined PRINTNUMERICALWARNINGS
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
 * where A is to be computed and is of size [3 x 3],
 * B is given and of size [3 x 3], L is a diagonal
 *  matrix of size 3 stored in diagonal 
 * storage format. Note the transpose of L!
 *
 * Result: A in diagonalzero storage format.
 *
 */
void smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_3(smoothMPC_FLOAT *L, smoothMPC_FLOAT *B, smoothMPC_FLOAT *A)
{
	int j;
    for( j=0; j<3; j++ ){   
		A[j] = B[j]/L[j];
     }
}


/**
 * Forward substitution to solve L*y = b where L is a
 * diagonal matrix in vector storage format.
 * 
 * The dimensions involved are 3.
 */
void smoothMPC_LA_DIAG_FORWARDSUB_3(smoothMPC_FLOAT *L, smoothMPC_FLOAT *b, smoothMPC_FLOAT *y)
{
    int i;

    for( i=0; i<3; i++ ){
		y[i] = b[i]/L[i];
    }
}


/**
 * Compute L = B*B', where L is lower triangular of size NXp1
 * and B is a diagzero matrix of size [3 x 11] in column
 * storage format.
 * 
 */
void smoothMPC_LA_DIAGZERO_MMT_3(smoothMPC_FLOAT *B, smoothMPC_FLOAT *L)
{
    int i, ii, di;
    
    ii = 0; di = 0;
    for( i=0; i<3; i++ ){        
		L[ii+i] = B[i]*B[i];
        ii += ++di;
    }
}


/* 
 * Computes r = b - B*u
 * B is stored in diagzero format
 */
void smoothMPC_LA_DIAGZERO_MVMSUB7_3(smoothMPC_FLOAT *B, smoothMPC_FLOAT *u, smoothMPC_FLOAT *b, smoothMPC_FLOAT *r)
{
	int i;

	for( i=0; i<3; i++ ){
		r[i] = b[i] - B[i]*u[i];
	}	
	
}


/**
 * Compute L = A*A' + B*B', where L is lower triangular of size NXp1
 * and A is a dense matrix of size [3 x 11] in column
 * storage format, and B is of size [3 x 11] diagonalzero
 * storage format.
 * 
 * THIS ONE HAS THE WORST ACCES PATTERN POSSIBLE. 
 * POSSIBKE FIX: PUT A AND B INTO ROW MAJOR FORMAT FIRST.
 * 
 */
void smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_FLOAT *A, smoothMPC_FLOAT *B, smoothMPC_FLOAT *L)
{
    int i, j, k, ii, di;
    smoothMPC_FLOAT ltemp;
    
    ii = 0; di = 0;
    for( i=0; i<3; i++ ){        
        for( j=0; j<=i; j++ ){
            ltemp = 0; 
            for( k=0; k<11; k++ ){
                ltemp += A[k*3+i]*A[k*3+j];
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
void smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_FLOAT *A, smoothMPC_FLOAT *x, smoothMPC_FLOAT *B, smoothMPC_FLOAT *u, smoothMPC_FLOAT *b, smoothMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;

	for( i=0; i<3; i++ ){
		r[i] = b[i] - A[k++]*x[0] - B[i]*u[i];
	}	

	for( j=1; j<11; j++ ){		
		for( i=0; i<3; i++ ){
			r[i] -= A[k++]*x[j];
		}
	}
	
}


/**
 * Compute L = A*A' + B*B', where L is lower triangular of size NXp1
 * and A is a dense matrix of size [3 x 11] in column
 * storage format, and B is of size [3 x 3] diagonalzero
 * storage format.
 * 
 * THIS ONE HAS THE WORST ACCES PATTERN POSSIBLE. 
 * POSSIBKE FIX: PUT A AND B INTO ROW MAJOR FORMAT FIRST.
 * 
 */
void smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_3(smoothMPC_FLOAT *A, smoothMPC_FLOAT *B, smoothMPC_FLOAT *L)
{
    int i, j, k, ii, di;
    smoothMPC_FLOAT ltemp;
    
    ii = 0; di = 0;
    for( i=0; i<3; i++ ){        
        for( j=0; j<=i; j++ ){
            ltemp = 0; 
            for( k=0; k<11; k++ ){
                ltemp += A[k*3+i]*A[k*3+j];
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
void smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_3(smoothMPC_FLOAT *A, smoothMPC_FLOAT *x, smoothMPC_FLOAT *B, smoothMPC_FLOAT *u, smoothMPC_FLOAT *b, smoothMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;

	for( i=0; i<3; i++ ){
		r[i] = b[i] - A[k++]*x[0] - B[i]*u[i];
	}	

	for( j=1; j<11; j++ ){		
		for( i=0; i<3; i++ ){
			r[i] -= A[k++]*x[j];
		}
	}
	
}


/**
 * Cholesky factorization as above, but working on a matrix in 
 * lower triangular storage format of size 3 and outputting
 * the Cholesky factor to matrix L in lower triangular format.
 */
void smoothMPC_LA_DENSE_CHOL_3(smoothMPC_FLOAT *A, smoothMPC_FLOAT *L)
{
    int i, j, k, di, dj;
	 int ii, jj;

    smoothMPC_FLOAT l;
    smoothMPC_FLOAT Mii;

	/* copy A to L first and then operate on L */
	/* COULD BE OPTIMIZED */
	ii=0; di=0;
	for( i=0; i<3; i++ ){
		for( j=0; j<=i; j++ ){
			L[ii+j] = A[ii+j];
		}
		ii += ++di;
	}    
	
	/* factor L */
	ii=0; di=0;
    for( i=0; i<3; i++ ){
        l = 0;
        for( k=0; k<i; k++ ){
            l += L[ii+k]*L[ii+k];
        }        
        
        Mii = L[ii+i] - l;
        
#if smoothMPC_SET_PRINTLEVEL > 0 && defined PRINTNUMERICALWARNINGS
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
        for( j=i+1; j<3; j++ ){
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
 * The dimensions involved are 3.
 */
void smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_FLOAT *L, smoothMPC_FLOAT *b, smoothMPC_FLOAT *y)
{
    int i,j,ii,di;
    smoothMPC_FLOAT yel;
            
    ii = 0; di = 0;
    for( i=0; i<3; i++ ){
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
 * where A is to be computed and is of size [3 x 3],
 * B is given and of size [3 x 3], L is a lower tri-
 * angular matrix of size 3 stored in lower triangular 
 * storage format. Note the transpose of L AND B!
 *
 * Result: A in column major storage format.
 *
 */
void smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_FLOAT *L, smoothMPC_FLOAT *B, smoothMPC_FLOAT *A)
{
    int i,j,k,ii,di;
    smoothMPC_FLOAT a;
    
    ii=0; di=0;
    for( j=0; j<3; j++ ){        
        for( i=0; i<3; i++ ){
            a = B[i*3+j];
            for( k=0; k<j; k++ ){
                a -= A[k*3+i]*L[ii+k];
            }    

			/* saturate for numerical stability */
			a = MIN(a, BIGM);
			a = MAX(a, -BIGM); 

			A[j*3+i] = a/L[ii+j];			
        }
        ii += ++di;
    }
}


/**
 * Compute L = L - A*A', where L is lower triangular of size 3
 * and A is a dense matrix of size [3 x 3] in column
 * storage format.
 * 
 * THIS ONE HAS THE WORST ACCES PATTERN POSSIBLE. 
 * POSSIBKE FIX: PUT A INTO ROW MAJOR FORMAT FIRST.
 * 
 */
void smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_FLOAT *A, smoothMPC_FLOAT *L)
{
    int i, j, k, ii, di;
    smoothMPC_FLOAT ltemp;
    
    ii = 0; di = 0;
    for( i=0; i<3; i++ ){        
        for( j=0; j<=i; j++ ){
            ltemp = 0; 
            for( k=0; k<3; k++ ){
                ltemp += A[k*3+i]*A[k*3+j];
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
void smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_FLOAT *A, smoothMPC_FLOAT *x, smoothMPC_FLOAT *b, smoothMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;

	for( i=0; i<3; i++ ){
		r[i] = b[i] - A[k++]*x[0];
	}	
	for( j=1; j<3; j++ ){		
		for( i=0; i<3; i++ ){
			r[i] -= A[k++]*x[j];
		}
	}
}


/**
 * Backward Substitution to solve L^T*x = y where L is a
 * lower triangular matrix in triangular storage format.
 * 
 * All involved dimensions are 3.
 */
void smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_FLOAT *L, smoothMPC_FLOAT *y, smoothMPC_FLOAT *x)
{
    int i, ii, di, j, jj, dj;
    smoothMPC_FLOAT xel;    
	int start = 3;
    
    /* now solve L^T*x = y by backward substitution */
    ii = start; di = 2;
    for( i=2; i>=0; i-- ){        
        xel = y[i];        
        jj = start; dj = 2;
        for( j=2; j>i; j-- ){
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
 * Matrix vector multiplication y = b - M'*x where M is of size [3 x 3]
 * and stored in column major format. Note the transpose of M!
 */
void smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_FLOAT *A, smoothMPC_FLOAT *x, smoothMPC_FLOAT *b, smoothMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0; 
	for( i=0; i<3; i++ ){
		r[i] = b[i];
		for( j=0; j<3; j++ ){
			r[i] -= A[k++]*x[j];
		}
	}
}


/*
 * Vector subtraction z = -x - y for vectors of length 322.
 */
void smoothMPC_LA_VSUB2_322(smoothMPC_FLOAT *x, smoothMPC_FLOAT *y, smoothMPC_FLOAT *z)
{
	int i;
	for( i=0; i<322; i++){
		z[i] = -x[i] - y[i];
	}
}


/**
 * Forward-Backward-Substitution to solve L*L^T*x = b where L is a
 * diagonal matrix of size 11 in vector
 * storage format.
 */
void smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_FLOAT *L, smoothMPC_FLOAT *b, smoothMPC_FLOAT *x)
{
    int i;
            
    /* solve Ly = b by forward and backward substitution */
    for( i=0; i<11; i++ ){
		x[i] = b[i]/(L[i]*L[i]);
    }
    
}


/**
 * Forward-Backward-Substitution to solve L*L^T*x = b where L is a
 * diagonal matrix of size 3 in vector
 * storage format.
 */
void smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_3(smoothMPC_FLOAT *L, smoothMPC_FLOAT *b, smoothMPC_FLOAT *x)
{
    int i;
            
    /* solve Ly = b by forward and backward substitution */
    for( i=0; i<3; i++ ){
		x[i] = b[i]/(L[i]*L[i]);
    }
    
}


/*
 * Vector subtraction z = x(xidx) - y where y, z and xidx are of length 11,
 * and x has length 11 and is indexed through yidx.
 */
void smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_FLOAT *x, int* xidx, smoothMPC_FLOAT *y, smoothMPC_FLOAT *z)
{
	int i;
	for( i=0; i<11; i++){
		z[i] = x[xidx[i]] - y[i];
	}
}


/*
 * Vector subtraction x = -u.*v - w for vectors of length 11.
 */
void smoothMPC_LA_VSUB3_11(smoothMPC_FLOAT *u, smoothMPC_FLOAT *v, smoothMPC_FLOAT *w, smoothMPC_FLOAT *x)
{
	int i;
	for( i=0; i<11; i++){
		x[i] = -u[i]*v[i] - w[i];
	}
}


/*
 * Vector subtraction z = -x - y(yidx) where y is of length 11
 * and z, x and yidx are of length 5.
 */
void smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_FLOAT *x, smoothMPC_FLOAT *y, int* yidx, smoothMPC_FLOAT *z)
{
	int i;
	for( i=0; i<5; i++){
		z[i] = -x[i] - y[yidx[i]];
	}
}


/*
 * Vector subtraction x = -u.*v - w for vectors of length 5.
 */
void smoothMPC_LA_VSUB3_5(smoothMPC_FLOAT *u, smoothMPC_FLOAT *v, smoothMPC_FLOAT *w, smoothMPC_FLOAT *x)
{
	int i;
	for( i=0; i<5; i++){
		x[i] = -u[i]*v[i] - w[i];
	}
}


/*
 * Vector subtraction z = x(xidx) - y where y, z and xidx are of length 3,
 * and x has length 3 and is indexed through yidx.
 */
void smoothMPC_LA_VSUB_INDEXED_3(smoothMPC_FLOAT *x, int* xidx, smoothMPC_FLOAT *y, smoothMPC_FLOAT *z)
{
	int i;
	for( i=0; i<3; i++){
		z[i] = x[xidx[i]] - y[i];
	}
}


/*
 * Vector subtraction x = -u.*v - w for vectors of length 3.
 */
void smoothMPC_LA_VSUB3_3(smoothMPC_FLOAT *u, smoothMPC_FLOAT *v, smoothMPC_FLOAT *w, smoothMPC_FLOAT *x)
{
	int i;
	for( i=0; i<3; i++){
		x[i] = -u[i]*v[i] - w[i];
	}
}


/*
 * Vector subtraction z = -x - y(yidx) where y is of length 3
 * and z, x and yidx are of length 3.
 */
void smoothMPC_LA_VSUB2_INDEXED_3(smoothMPC_FLOAT *x, smoothMPC_FLOAT *y, int* yidx, smoothMPC_FLOAT *z)
{
	int i;
	for( i=0; i<3; i++){
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
 * smoothMPC_NOPROGRESS (should be negative).
 */
int smoothMPC_LINESEARCH_BACKTRACKING_AFFINE(smoothMPC_FLOAT *l, smoothMPC_FLOAT *s, smoothMPC_FLOAT *dl, smoothMPC_FLOAT *ds, smoothMPC_FLOAT *a, smoothMPC_FLOAT *mu_aff)
{
    int i;
	int lsIt=1;    
    smoothMPC_FLOAT dltemp;
    smoothMPC_FLOAT dstemp;
    smoothMPC_FLOAT mya = 1.0;
    smoothMPC_FLOAT mymu;
        
    while( 1 ){                        

        /* 
         * Compute both snew and wnew together.
         * We compute also mu_affine along the way here, as the
         * values might be in registers, so it should be cheaper.
         */
        mymu = 0;
        for( i=0; i<470; i++ ){
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
        if( i == 470 ){
            break;
        } else {
            mya *= smoothMPC_SET_LS_SCALE_AFF;
            if( mya < smoothMPC_SET_LS_MINSTEP ){
                return smoothMPC_NOPROGRESS;
            }
        }
    }
    
    /* return new values and iteration counter */
    *a = mya;
    *mu_aff = mymu / (smoothMPC_FLOAT)470;
    return lsIt;
}


/*
 * Vector subtraction x = u.*v - a where a is a scalar
*  and x,u,v are vectors of length 470.
 */
void smoothMPC_LA_VSUB5_470(smoothMPC_FLOAT *u, smoothMPC_FLOAT *v, smoothMPC_FLOAT a, smoothMPC_FLOAT *x)
{
	int i;
	for( i=0; i<470; i++){
		x[i] = u[i]*v[i] - a;
	}
}


/*
 * Computes x=0; x(uidx) += u/su; x(vidx) -= v/sv where x is of length 11,
 * u, su, uidx are of length 5 and v, sv, vidx are of length 11.
 */
void smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_FLOAT *u, smoothMPC_FLOAT *su, int* uidx, smoothMPC_FLOAT *v, smoothMPC_FLOAT *sv, int* vidx, smoothMPC_FLOAT *x)
{
	int i;
	for( i=0; i<11; i++ ){
		x[i] = 0;
	}
	for( i=0; i<5; i++){
		x[uidx[i]] += u[i]/su[i];
	}
	for( i=0; i<11; i++){
		x[vidx[i]] -= v[i]/sv[i];
	}
}


/* 
 * Computes r =  B*u
 * where B is stored in diagzero format
 */
void smoothMPC_LA_DIAGZERO_MVM_3(smoothMPC_FLOAT *B, smoothMPC_FLOAT *u, smoothMPC_FLOAT *r)
{
	int i;

	for( i=0; i<3; i++ ){
		r[i] = B[i]*u[i];
	}	
	
}


/* 
 * Computes r = A*x + B*u
 * where A is stored in column major format
 * and B is stored in diagzero format
 */
void smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_FLOAT *A, smoothMPC_FLOAT *x, smoothMPC_FLOAT *B, smoothMPC_FLOAT *u, smoothMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;

	for( i=0; i<3; i++ ){
		r[i] = A[k++]*x[0] + B[i]*u[i];
	}	

	for( j=1; j<11; j++ ){		
		for( i=0; i<3; i++ ){
			r[i] += A[k++]*x[j];
		}
	}
	
}


/*
 * Computes x=0; x(uidx) += u/su; x(vidx) -= v/sv where x is of length 3,
 * u, su, uidx are of length 3 and v, sv, vidx are of length 3.
 */
void smoothMPC_LA_VSUB6_INDEXED_3_3_3(smoothMPC_FLOAT *u, smoothMPC_FLOAT *su, int* uidx, smoothMPC_FLOAT *v, smoothMPC_FLOAT *sv, int* vidx, smoothMPC_FLOAT *x)
{
	int i;
	for( i=0; i<3; i++ ){
		x[i] = 0;
	}
	for( i=0; i<3; i++){
		x[uidx[i]] += u[i]/su[i];
	}
	for( i=0; i<3; i++){
		x[vidx[i]] -= v[i]/sv[i];
	}
}


/* 
 * Computes r = A*x + B*u
 * where A is stored in column major format
 * and B is stored in diagzero format
 */
void smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_3(smoothMPC_FLOAT *A, smoothMPC_FLOAT *x, smoothMPC_FLOAT *B, smoothMPC_FLOAT *u, smoothMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;

	for( i=0; i<3; i++ ){
		r[i] = A[k++]*x[0] + B[i]*u[i];
	}	

	for( j=1; j<11; j++ ){		
		for( i=0; i<3; i++ ){
			r[i] += A[k++]*x[j];
		}
	}
	
}


/*
 * Vector subtraction z = x - y for vectors of length 322.
 */
void smoothMPC_LA_VSUB_322(smoothMPC_FLOAT *x, smoothMPC_FLOAT *y, smoothMPC_FLOAT *z)
{
	int i;
	for( i=0; i<322; i++){
		z[i] = x[i] - y[i];
	}
}


/** 
 * Computes z = -r./s - u.*y(y)
 * where all vectors except of y are of length 11 (length of y >= 11).
 */
void smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_FLOAT *r, smoothMPC_FLOAT *s, smoothMPC_FLOAT *u, smoothMPC_FLOAT *y, int* yidx, smoothMPC_FLOAT *z)
{
	int i;
	for( i=0; i<11; i++ ){
		z[i] = -r[i]/s[i] - u[i]*y[yidx[i]];
	}
}


/** 
 * Computes z = -r./s + u.*y(y)
 * where all vectors except of y are of length 5 (length of y >= 5).
 */
void smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_FLOAT *r, smoothMPC_FLOAT *s, smoothMPC_FLOAT *u, smoothMPC_FLOAT *y, int* yidx, smoothMPC_FLOAT *z)
{
	int i;
	for( i=0; i<5; i++ ){
		z[i] = -r[i]/s[i] + u[i]*y[yidx[i]];
	}
}


/** 
 * Computes z = -r./s - u.*y(y)
 * where all vectors except of y are of length 3 (length of y >= 3).
 */
void smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_3(smoothMPC_FLOAT *r, smoothMPC_FLOAT *s, smoothMPC_FLOAT *u, smoothMPC_FLOAT *y, int* yidx, smoothMPC_FLOAT *z)
{
	int i;
	for( i=0; i<3; i++ ){
		z[i] = -r[i]/s[i] - u[i]*y[yidx[i]];
	}
}


/** 
 * Computes z = -r./s + u.*y(y)
 * where all vectors except of y are of length 3 (length of y >= 3).
 */
void smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_3(smoothMPC_FLOAT *r, smoothMPC_FLOAT *s, smoothMPC_FLOAT *u, smoothMPC_FLOAT *y, int* yidx, smoothMPC_FLOAT *z)
{
	int i;
	for( i=0; i<3; i++ ){
		z[i] = -r[i]/s[i] + u[i]*y[yidx[i]];
	}
}


/*
 * Computes ds = -l.\(r + s.*dl) for vectors of length 470.
 */
void smoothMPC_LA_VSUB7_470(smoothMPC_FLOAT *l, smoothMPC_FLOAT *r, smoothMPC_FLOAT *s, smoothMPC_FLOAT *dl, smoothMPC_FLOAT *ds)
{
	int i;
	for( i=0; i<470; i++){
		ds[i] = -(r[i] + s[i]*dl[i])/l[i];
	}
}


/*
 * Vector addition x = x + y for vectors of length 322.
 */
void smoothMPC_LA_VADD_322(smoothMPC_FLOAT *x, smoothMPC_FLOAT *y)
{
	int i;
	for( i=0; i<322; i++){
		x[i] += y[i];
	}
}


/*
 * Vector addition x = x + y for vectors of length 90.
 */
void smoothMPC_LA_VADD_90(smoothMPC_FLOAT *x, smoothMPC_FLOAT *y)
{
	int i;
	for( i=0; i<90; i++){
		x[i] += y[i];
	}
}


/*
 * Vector addition x = x + y for vectors of length 470.
 */
void smoothMPC_LA_VADD_470(smoothMPC_FLOAT *x, smoothMPC_FLOAT *y)
{
	int i;
	for( i=0; i<470; i++){
		x[i] += y[i];
	}
}


/**
 * Backtracking line search for combined predictor/corrector step.
 * Update on variables with safety factor gamma (to keep us away from
 * boundary).
 */
int smoothMPC_LINESEARCH_BACKTRACKING_COMBINED(smoothMPC_FLOAT *z, smoothMPC_FLOAT *v, smoothMPC_FLOAT *l, smoothMPC_FLOAT *s, smoothMPC_FLOAT *dz, smoothMPC_FLOAT *dv, smoothMPC_FLOAT *dl, smoothMPC_FLOAT *ds, smoothMPC_FLOAT *a, smoothMPC_FLOAT *mu)
{
    int i, lsIt=1;       
    smoothMPC_FLOAT dltemp;
    smoothMPC_FLOAT dstemp;    
    smoothMPC_FLOAT a_gamma;
            
    *a = 1.0;
    while( 1 ){                        

        /* check whether search criterion is fulfilled */
        for( i=0; i<470; i++ ){
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
        if( i == 470 ){
            break;
        } else {
            *a *= smoothMPC_SET_LS_SCALE;
            if( *a < smoothMPC_SET_LS_MINSTEP ){
                return smoothMPC_NOPROGRESS;
            }
        }
    }
    
    /* update variables with safety margin */
    a_gamma = (*a)*smoothMPC_SET_LS_MAXSTEP;
    
    /* primal variables */
    for( i=0; i<322; i++ ){
        z[i] += a_gamma*dz[i];
    }
    
    /* equality constraint multipliers */
    for( i=0; i<90; i++ ){
        v[i] += a_gamma*dv[i];
    }
    
    /* inequality constraint multipliers & slacks, also update mu */
    *mu = 0;
    for( i=0; i<470; i++ ){
        dltemp = l[i] + a_gamma*dl[i]; l[i] = dltemp;
        dstemp = s[i] + a_gamma*ds[i]; s[i] = dstemp;
        *mu += dltemp*dstemp;
    }
    
    *a = a_gamma;
    *mu /= (smoothMPC_FLOAT)470;
    return lsIt;
}




/* VARIABLE DEFINITIONS ------------------------------------------------ */
smoothMPC_FLOAT smoothMPC_z[322];
smoothMPC_FLOAT smoothMPC_v[90];
smoothMPC_FLOAT smoothMPC_dz_aff[322];
smoothMPC_FLOAT smoothMPC_dv_aff[90];
smoothMPC_FLOAT smoothMPC_grad_cost[322];
smoothMPC_FLOAT smoothMPC_grad_eq[322];
smoothMPC_FLOAT smoothMPC_rd[322];
smoothMPC_FLOAT smoothMPC_l[470];
smoothMPC_FLOAT smoothMPC_s[470];
smoothMPC_FLOAT smoothMPC_lbys[470];
smoothMPC_FLOAT smoothMPC_dl_aff[470];
smoothMPC_FLOAT smoothMPC_ds_aff[470];
smoothMPC_FLOAT smoothMPC_dz_cc[322];
smoothMPC_FLOAT smoothMPC_dv_cc[90];
smoothMPC_FLOAT smoothMPC_dl_cc[470];
smoothMPC_FLOAT smoothMPC_ds_cc[470];
smoothMPC_FLOAT smoothMPC_ccrhs[470];
smoothMPC_FLOAT smoothMPC_grad_ineq[322];
smoothMPC_FLOAT* smoothMPC_z00 = smoothMPC_z + 0;
smoothMPC_FLOAT* smoothMPC_dzaff00 = smoothMPC_dz_aff + 0;
smoothMPC_FLOAT* smoothMPC_dzcc00 = smoothMPC_dz_cc + 0;
smoothMPC_FLOAT* smoothMPC_rd00 = smoothMPC_rd + 0;
smoothMPC_FLOAT smoothMPC_Lbyrd00[11];
smoothMPC_FLOAT* smoothMPC_grad_cost00 = smoothMPC_grad_cost + 0;
smoothMPC_FLOAT* smoothMPC_grad_eq00 = smoothMPC_grad_eq + 0;
smoothMPC_FLOAT* smoothMPC_grad_ineq00 = smoothMPC_grad_ineq + 0;
smoothMPC_FLOAT smoothMPC_ctv00[11];
smoothMPC_FLOAT* smoothMPC_v00 = smoothMPC_v + 0;
smoothMPC_FLOAT smoothMPC_re00[3];
smoothMPC_FLOAT smoothMPC_beta00[3];
smoothMPC_FLOAT smoothMPC_betacc00[3];
smoothMPC_FLOAT* smoothMPC_dvaff00 = smoothMPC_dv_aff + 0;
smoothMPC_FLOAT* smoothMPC_dvcc00 = smoothMPC_dv_cc + 0;
smoothMPC_FLOAT smoothMPC_V00[33];
smoothMPC_FLOAT smoothMPC_Yd00[6];
smoothMPC_FLOAT smoothMPC_Ld00[6];
smoothMPC_FLOAT smoothMPC_yy00[3];
smoothMPC_FLOAT smoothMPC_bmy00[3];
int smoothMPC_lbIdx00[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb00 = smoothMPC_l + 0;
smoothMPC_FLOAT* smoothMPC_slb00 = smoothMPC_s + 0;
smoothMPC_FLOAT* smoothMPC_llbbyslb00 = smoothMPC_lbys + 0;
smoothMPC_FLOAT smoothMPC_rilb00[11];
smoothMPC_FLOAT* smoothMPC_dllbaff00 = smoothMPC_dl_aff + 0;
smoothMPC_FLOAT* smoothMPC_dslbaff00 = smoothMPC_ds_aff + 0;
smoothMPC_FLOAT* smoothMPC_dllbcc00 = smoothMPC_dl_cc + 0;
smoothMPC_FLOAT* smoothMPC_dslbcc00 = smoothMPC_ds_cc + 0;
smoothMPC_FLOAT* smoothMPC_ccrhsl00 = smoothMPC_ccrhs + 0;
int smoothMPC_ubIdx00[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub00 = smoothMPC_l + 11;
smoothMPC_FLOAT* smoothMPC_sub00 = smoothMPC_s + 11;
smoothMPC_FLOAT* smoothMPC_lubbysub00 = smoothMPC_lbys + 11;
smoothMPC_FLOAT smoothMPC_riub00[5];
smoothMPC_FLOAT* smoothMPC_dlubaff00 = smoothMPC_dl_aff + 11;
smoothMPC_FLOAT* smoothMPC_dsubaff00 = smoothMPC_ds_aff + 11;
smoothMPC_FLOAT* smoothMPC_dlubcc00 = smoothMPC_dl_cc + 11;
smoothMPC_FLOAT* smoothMPC_dsubcc00 = smoothMPC_ds_cc + 11;
smoothMPC_FLOAT* smoothMPC_ccrhsub00 = smoothMPC_ccrhs + 11;
smoothMPC_FLOAT smoothMPC_Phi00[11];
smoothMPC_FLOAT smoothMPC_D00[11] = {1.0000000000000000E+000, 
1.0000000000000000E+000, 
1.0000000000000000E+000};
smoothMPC_FLOAT smoothMPC_W00[11];
smoothMPC_FLOAT* smoothMPC_z01 = smoothMPC_z + 11;
smoothMPC_FLOAT* smoothMPC_dzaff01 = smoothMPC_dz_aff + 11;
smoothMPC_FLOAT* smoothMPC_dzcc01 = smoothMPC_dz_cc + 11;
smoothMPC_FLOAT* smoothMPC_rd01 = smoothMPC_rd + 11;
smoothMPC_FLOAT smoothMPC_Lbyrd01[11];
smoothMPC_FLOAT* smoothMPC_grad_cost01 = smoothMPC_grad_cost + 11;
smoothMPC_FLOAT* smoothMPC_grad_eq01 = smoothMPC_grad_eq + 11;
smoothMPC_FLOAT* smoothMPC_grad_ineq01 = smoothMPC_grad_ineq + 11;
smoothMPC_FLOAT smoothMPC_ctv01[11];
smoothMPC_FLOAT* smoothMPC_v01 = smoothMPC_v + 3;
smoothMPC_FLOAT smoothMPC_re01[3];
smoothMPC_FLOAT smoothMPC_beta01[3];
smoothMPC_FLOAT smoothMPC_betacc01[3];
smoothMPC_FLOAT* smoothMPC_dvaff01 = smoothMPC_dv_aff + 3;
smoothMPC_FLOAT* smoothMPC_dvcc01 = smoothMPC_dv_cc + 3;
smoothMPC_FLOAT smoothMPC_V01[33];
smoothMPC_FLOAT smoothMPC_Yd01[6];
smoothMPC_FLOAT smoothMPC_Ld01[6];
smoothMPC_FLOAT smoothMPC_yy01[3];
smoothMPC_FLOAT smoothMPC_bmy01[3];
int smoothMPC_lbIdx01[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb01 = smoothMPC_l + 16;
smoothMPC_FLOAT* smoothMPC_slb01 = smoothMPC_s + 16;
smoothMPC_FLOAT* smoothMPC_llbbyslb01 = smoothMPC_lbys + 16;
smoothMPC_FLOAT smoothMPC_rilb01[11];
smoothMPC_FLOAT* smoothMPC_dllbaff01 = smoothMPC_dl_aff + 16;
smoothMPC_FLOAT* smoothMPC_dslbaff01 = smoothMPC_ds_aff + 16;
smoothMPC_FLOAT* smoothMPC_dllbcc01 = smoothMPC_dl_cc + 16;
smoothMPC_FLOAT* smoothMPC_dslbcc01 = smoothMPC_ds_cc + 16;
smoothMPC_FLOAT* smoothMPC_ccrhsl01 = smoothMPC_ccrhs + 16;
int smoothMPC_ubIdx01[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub01 = smoothMPC_l + 27;
smoothMPC_FLOAT* smoothMPC_sub01 = smoothMPC_s + 27;
smoothMPC_FLOAT* smoothMPC_lubbysub01 = smoothMPC_lbys + 27;
smoothMPC_FLOAT smoothMPC_riub01[5];
smoothMPC_FLOAT* smoothMPC_dlubaff01 = smoothMPC_dl_aff + 27;
smoothMPC_FLOAT* smoothMPC_dsubaff01 = smoothMPC_ds_aff + 27;
smoothMPC_FLOAT* smoothMPC_dlubcc01 = smoothMPC_dl_cc + 27;
smoothMPC_FLOAT* smoothMPC_dsubcc01 = smoothMPC_ds_cc + 27;
smoothMPC_FLOAT* smoothMPC_ccrhsub01 = smoothMPC_ccrhs + 27;
smoothMPC_FLOAT smoothMPC_Phi01[11];
smoothMPC_FLOAT smoothMPC_D01[11] = {-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000};
smoothMPC_FLOAT smoothMPC_W01[11];
smoothMPC_FLOAT smoothMPC_Ysd01[9];
smoothMPC_FLOAT smoothMPC_Lsd01[9];
smoothMPC_FLOAT* smoothMPC_z02 = smoothMPC_z + 22;
smoothMPC_FLOAT* smoothMPC_dzaff02 = smoothMPC_dz_aff + 22;
smoothMPC_FLOAT* smoothMPC_dzcc02 = smoothMPC_dz_cc + 22;
smoothMPC_FLOAT* smoothMPC_rd02 = smoothMPC_rd + 22;
smoothMPC_FLOAT smoothMPC_Lbyrd02[11];
smoothMPC_FLOAT* smoothMPC_grad_cost02 = smoothMPC_grad_cost + 22;
smoothMPC_FLOAT* smoothMPC_grad_eq02 = smoothMPC_grad_eq + 22;
smoothMPC_FLOAT* smoothMPC_grad_ineq02 = smoothMPC_grad_ineq + 22;
smoothMPC_FLOAT smoothMPC_ctv02[11];
smoothMPC_FLOAT* smoothMPC_v02 = smoothMPC_v + 6;
smoothMPC_FLOAT smoothMPC_re02[3];
smoothMPC_FLOAT smoothMPC_beta02[3];
smoothMPC_FLOAT smoothMPC_betacc02[3];
smoothMPC_FLOAT* smoothMPC_dvaff02 = smoothMPC_dv_aff + 6;
smoothMPC_FLOAT* smoothMPC_dvcc02 = smoothMPC_dv_cc + 6;
smoothMPC_FLOAT smoothMPC_V02[33];
smoothMPC_FLOAT smoothMPC_Yd02[6];
smoothMPC_FLOAT smoothMPC_Ld02[6];
smoothMPC_FLOAT smoothMPC_yy02[3];
smoothMPC_FLOAT smoothMPC_bmy02[3];
int smoothMPC_lbIdx02[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb02 = smoothMPC_l + 32;
smoothMPC_FLOAT* smoothMPC_slb02 = smoothMPC_s + 32;
smoothMPC_FLOAT* smoothMPC_llbbyslb02 = smoothMPC_lbys + 32;
smoothMPC_FLOAT smoothMPC_rilb02[11];
smoothMPC_FLOAT* smoothMPC_dllbaff02 = smoothMPC_dl_aff + 32;
smoothMPC_FLOAT* smoothMPC_dslbaff02 = smoothMPC_ds_aff + 32;
smoothMPC_FLOAT* smoothMPC_dllbcc02 = smoothMPC_dl_cc + 32;
smoothMPC_FLOAT* smoothMPC_dslbcc02 = smoothMPC_ds_cc + 32;
smoothMPC_FLOAT* smoothMPC_ccrhsl02 = smoothMPC_ccrhs + 32;
int smoothMPC_ubIdx02[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub02 = smoothMPC_l + 43;
smoothMPC_FLOAT* smoothMPC_sub02 = smoothMPC_s + 43;
smoothMPC_FLOAT* smoothMPC_lubbysub02 = smoothMPC_lbys + 43;
smoothMPC_FLOAT smoothMPC_riub02[5];
smoothMPC_FLOAT* smoothMPC_dlubaff02 = smoothMPC_dl_aff + 43;
smoothMPC_FLOAT* smoothMPC_dsubaff02 = smoothMPC_ds_aff + 43;
smoothMPC_FLOAT* smoothMPC_dlubcc02 = smoothMPC_dl_cc + 43;
smoothMPC_FLOAT* smoothMPC_dsubcc02 = smoothMPC_ds_cc + 43;
smoothMPC_FLOAT* smoothMPC_ccrhsub02 = smoothMPC_ccrhs + 43;
smoothMPC_FLOAT smoothMPC_Phi02[11];
smoothMPC_FLOAT smoothMPC_W02[11];
smoothMPC_FLOAT smoothMPC_Ysd02[9];
smoothMPC_FLOAT smoothMPC_Lsd02[9];
smoothMPC_FLOAT* smoothMPC_z03 = smoothMPC_z + 33;
smoothMPC_FLOAT* smoothMPC_dzaff03 = smoothMPC_dz_aff + 33;
smoothMPC_FLOAT* smoothMPC_dzcc03 = smoothMPC_dz_cc + 33;
smoothMPC_FLOAT* smoothMPC_rd03 = smoothMPC_rd + 33;
smoothMPC_FLOAT smoothMPC_Lbyrd03[11];
smoothMPC_FLOAT* smoothMPC_grad_cost03 = smoothMPC_grad_cost + 33;
smoothMPC_FLOAT* smoothMPC_grad_eq03 = smoothMPC_grad_eq + 33;
smoothMPC_FLOAT* smoothMPC_grad_ineq03 = smoothMPC_grad_ineq + 33;
smoothMPC_FLOAT smoothMPC_ctv03[11];
smoothMPC_FLOAT* smoothMPC_v03 = smoothMPC_v + 9;
smoothMPC_FLOAT smoothMPC_re03[3];
smoothMPC_FLOAT smoothMPC_beta03[3];
smoothMPC_FLOAT smoothMPC_betacc03[3];
smoothMPC_FLOAT* smoothMPC_dvaff03 = smoothMPC_dv_aff + 9;
smoothMPC_FLOAT* smoothMPC_dvcc03 = smoothMPC_dv_cc + 9;
smoothMPC_FLOAT smoothMPC_V03[33];
smoothMPC_FLOAT smoothMPC_Yd03[6];
smoothMPC_FLOAT smoothMPC_Ld03[6];
smoothMPC_FLOAT smoothMPC_yy03[3];
smoothMPC_FLOAT smoothMPC_bmy03[3];
int smoothMPC_lbIdx03[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb03 = smoothMPC_l + 48;
smoothMPC_FLOAT* smoothMPC_slb03 = smoothMPC_s + 48;
smoothMPC_FLOAT* smoothMPC_llbbyslb03 = smoothMPC_lbys + 48;
smoothMPC_FLOAT smoothMPC_rilb03[11];
smoothMPC_FLOAT* smoothMPC_dllbaff03 = smoothMPC_dl_aff + 48;
smoothMPC_FLOAT* smoothMPC_dslbaff03 = smoothMPC_ds_aff + 48;
smoothMPC_FLOAT* smoothMPC_dllbcc03 = smoothMPC_dl_cc + 48;
smoothMPC_FLOAT* smoothMPC_dslbcc03 = smoothMPC_ds_cc + 48;
smoothMPC_FLOAT* smoothMPC_ccrhsl03 = smoothMPC_ccrhs + 48;
int smoothMPC_ubIdx03[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub03 = smoothMPC_l + 59;
smoothMPC_FLOAT* smoothMPC_sub03 = smoothMPC_s + 59;
smoothMPC_FLOAT* smoothMPC_lubbysub03 = smoothMPC_lbys + 59;
smoothMPC_FLOAT smoothMPC_riub03[5];
smoothMPC_FLOAT* smoothMPC_dlubaff03 = smoothMPC_dl_aff + 59;
smoothMPC_FLOAT* smoothMPC_dsubaff03 = smoothMPC_ds_aff + 59;
smoothMPC_FLOAT* smoothMPC_dlubcc03 = smoothMPC_dl_cc + 59;
smoothMPC_FLOAT* smoothMPC_dsubcc03 = smoothMPC_ds_cc + 59;
smoothMPC_FLOAT* smoothMPC_ccrhsub03 = smoothMPC_ccrhs + 59;
smoothMPC_FLOAT smoothMPC_Phi03[11];
smoothMPC_FLOAT smoothMPC_W03[11];
smoothMPC_FLOAT smoothMPC_Ysd03[9];
smoothMPC_FLOAT smoothMPC_Lsd03[9];
smoothMPC_FLOAT* smoothMPC_z04 = smoothMPC_z + 44;
smoothMPC_FLOAT* smoothMPC_dzaff04 = smoothMPC_dz_aff + 44;
smoothMPC_FLOAT* smoothMPC_dzcc04 = smoothMPC_dz_cc + 44;
smoothMPC_FLOAT* smoothMPC_rd04 = smoothMPC_rd + 44;
smoothMPC_FLOAT smoothMPC_Lbyrd04[11];
smoothMPC_FLOAT* smoothMPC_grad_cost04 = smoothMPC_grad_cost + 44;
smoothMPC_FLOAT* smoothMPC_grad_eq04 = smoothMPC_grad_eq + 44;
smoothMPC_FLOAT* smoothMPC_grad_ineq04 = smoothMPC_grad_ineq + 44;
smoothMPC_FLOAT smoothMPC_ctv04[11];
smoothMPC_FLOAT* smoothMPC_v04 = smoothMPC_v + 12;
smoothMPC_FLOAT smoothMPC_re04[3];
smoothMPC_FLOAT smoothMPC_beta04[3];
smoothMPC_FLOAT smoothMPC_betacc04[3];
smoothMPC_FLOAT* smoothMPC_dvaff04 = smoothMPC_dv_aff + 12;
smoothMPC_FLOAT* smoothMPC_dvcc04 = smoothMPC_dv_cc + 12;
smoothMPC_FLOAT smoothMPC_V04[33];
smoothMPC_FLOAT smoothMPC_Yd04[6];
smoothMPC_FLOAT smoothMPC_Ld04[6];
smoothMPC_FLOAT smoothMPC_yy04[3];
smoothMPC_FLOAT smoothMPC_bmy04[3];
int smoothMPC_lbIdx04[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb04 = smoothMPC_l + 64;
smoothMPC_FLOAT* smoothMPC_slb04 = smoothMPC_s + 64;
smoothMPC_FLOAT* smoothMPC_llbbyslb04 = smoothMPC_lbys + 64;
smoothMPC_FLOAT smoothMPC_rilb04[11];
smoothMPC_FLOAT* smoothMPC_dllbaff04 = smoothMPC_dl_aff + 64;
smoothMPC_FLOAT* smoothMPC_dslbaff04 = smoothMPC_ds_aff + 64;
smoothMPC_FLOAT* smoothMPC_dllbcc04 = smoothMPC_dl_cc + 64;
smoothMPC_FLOAT* smoothMPC_dslbcc04 = smoothMPC_ds_cc + 64;
smoothMPC_FLOAT* smoothMPC_ccrhsl04 = smoothMPC_ccrhs + 64;
int smoothMPC_ubIdx04[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub04 = smoothMPC_l + 75;
smoothMPC_FLOAT* smoothMPC_sub04 = smoothMPC_s + 75;
smoothMPC_FLOAT* smoothMPC_lubbysub04 = smoothMPC_lbys + 75;
smoothMPC_FLOAT smoothMPC_riub04[5];
smoothMPC_FLOAT* smoothMPC_dlubaff04 = smoothMPC_dl_aff + 75;
smoothMPC_FLOAT* smoothMPC_dsubaff04 = smoothMPC_ds_aff + 75;
smoothMPC_FLOAT* smoothMPC_dlubcc04 = smoothMPC_dl_cc + 75;
smoothMPC_FLOAT* smoothMPC_dsubcc04 = smoothMPC_ds_cc + 75;
smoothMPC_FLOAT* smoothMPC_ccrhsub04 = smoothMPC_ccrhs + 75;
smoothMPC_FLOAT smoothMPC_Phi04[11];
smoothMPC_FLOAT smoothMPC_W04[11];
smoothMPC_FLOAT smoothMPC_Ysd04[9];
smoothMPC_FLOAT smoothMPC_Lsd04[9];
smoothMPC_FLOAT* smoothMPC_z05 = smoothMPC_z + 55;
smoothMPC_FLOAT* smoothMPC_dzaff05 = smoothMPC_dz_aff + 55;
smoothMPC_FLOAT* smoothMPC_dzcc05 = smoothMPC_dz_cc + 55;
smoothMPC_FLOAT* smoothMPC_rd05 = smoothMPC_rd + 55;
smoothMPC_FLOAT smoothMPC_Lbyrd05[11];
smoothMPC_FLOAT* smoothMPC_grad_cost05 = smoothMPC_grad_cost + 55;
smoothMPC_FLOAT* smoothMPC_grad_eq05 = smoothMPC_grad_eq + 55;
smoothMPC_FLOAT* smoothMPC_grad_ineq05 = smoothMPC_grad_ineq + 55;
smoothMPC_FLOAT smoothMPC_ctv05[11];
smoothMPC_FLOAT* smoothMPC_v05 = smoothMPC_v + 15;
smoothMPC_FLOAT smoothMPC_re05[3];
smoothMPC_FLOAT smoothMPC_beta05[3];
smoothMPC_FLOAT smoothMPC_betacc05[3];
smoothMPC_FLOAT* smoothMPC_dvaff05 = smoothMPC_dv_aff + 15;
smoothMPC_FLOAT* smoothMPC_dvcc05 = smoothMPC_dv_cc + 15;
smoothMPC_FLOAT smoothMPC_V05[33];
smoothMPC_FLOAT smoothMPC_Yd05[6];
smoothMPC_FLOAT smoothMPC_Ld05[6];
smoothMPC_FLOAT smoothMPC_yy05[3];
smoothMPC_FLOAT smoothMPC_bmy05[3];
int smoothMPC_lbIdx05[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb05 = smoothMPC_l + 80;
smoothMPC_FLOAT* smoothMPC_slb05 = smoothMPC_s + 80;
smoothMPC_FLOAT* smoothMPC_llbbyslb05 = smoothMPC_lbys + 80;
smoothMPC_FLOAT smoothMPC_rilb05[11];
smoothMPC_FLOAT* smoothMPC_dllbaff05 = smoothMPC_dl_aff + 80;
smoothMPC_FLOAT* smoothMPC_dslbaff05 = smoothMPC_ds_aff + 80;
smoothMPC_FLOAT* smoothMPC_dllbcc05 = smoothMPC_dl_cc + 80;
smoothMPC_FLOAT* smoothMPC_dslbcc05 = smoothMPC_ds_cc + 80;
smoothMPC_FLOAT* smoothMPC_ccrhsl05 = smoothMPC_ccrhs + 80;
int smoothMPC_ubIdx05[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub05 = smoothMPC_l + 91;
smoothMPC_FLOAT* smoothMPC_sub05 = smoothMPC_s + 91;
smoothMPC_FLOAT* smoothMPC_lubbysub05 = smoothMPC_lbys + 91;
smoothMPC_FLOAT smoothMPC_riub05[5];
smoothMPC_FLOAT* smoothMPC_dlubaff05 = smoothMPC_dl_aff + 91;
smoothMPC_FLOAT* smoothMPC_dsubaff05 = smoothMPC_ds_aff + 91;
smoothMPC_FLOAT* smoothMPC_dlubcc05 = smoothMPC_dl_cc + 91;
smoothMPC_FLOAT* smoothMPC_dsubcc05 = smoothMPC_ds_cc + 91;
smoothMPC_FLOAT* smoothMPC_ccrhsub05 = smoothMPC_ccrhs + 91;
smoothMPC_FLOAT smoothMPC_Phi05[11];
smoothMPC_FLOAT smoothMPC_W05[11];
smoothMPC_FLOAT smoothMPC_Ysd05[9];
smoothMPC_FLOAT smoothMPC_Lsd05[9];
smoothMPC_FLOAT* smoothMPC_z06 = smoothMPC_z + 66;
smoothMPC_FLOAT* smoothMPC_dzaff06 = smoothMPC_dz_aff + 66;
smoothMPC_FLOAT* smoothMPC_dzcc06 = smoothMPC_dz_cc + 66;
smoothMPC_FLOAT* smoothMPC_rd06 = smoothMPC_rd + 66;
smoothMPC_FLOAT smoothMPC_Lbyrd06[11];
smoothMPC_FLOAT* smoothMPC_grad_cost06 = smoothMPC_grad_cost + 66;
smoothMPC_FLOAT* smoothMPC_grad_eq06 = smoothMPC_grad_eq + 66;
smoothMPC_FLOAT* smoothMPC_grad_ineq06 = smoothMPC_grad_ineq + 66;
smoothMPC_FLOAT smoothMPC_ctv06[11];
smoothMPC_FLOAT* smoothMPC_v06 = smoothMPC_v + 18;
smoothMPC_FLOAT smoothMPC_re06[3];
smoothMPC_FLOAT smoothMPC_beta06[3];
smoothMPC_FLOAT smoothMPC_betacc06[3];
smoothMPC_FLOAT* smoothMPC_dvaff06 = smoothMPC_dv_aff + 18;
smoothMPC_FLOAT* smoothMPC_dvcc06 = smoothMPC_dv_cc + 18;
smoothMPC_FLOAT smoothMPC_V06[33];
smoothMPC_FLOAT smoothMPC_Yd06[6];
smoothMPC_FLOAT smoothMPC_Ld06[6];
smoothMPC_FLOAT smoothMPC_yy06[3];
smoothMPC_FLOAT smoothMPC_bmy06[3];
int smoothMPC_lbIdx06[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb06 = smoothMPC_l + 96;
smoothMPC_FLOAT* smoothMPC_slb06 = smoothMPC_s + 96;
smoothMPC_FLOAT* smoothMPC_llbbyslb06 = smoothMPC_lbys + 96;
smoothMPC_FLOAT smoothMPC_rilb06[11];
smoothMPC_FLOAT* smoothMPC_dllbaff06 = smoothMPC_dl_aff + 96;
smoothMPC_FLOAT* smoothMPC_dslbaff06 = smoothMPC_ds_aff + 96;
smoothMPC_FLOAT* smoothMPC_dllbcc06 = smoothMPC_dl_cc + 96;
smoothMPC_FLOAT* smoothMPC_dslbcc06 = smoothMPC_ds_cc + 96;
smoothMPC_FLOAT* smoothMPC_ccrhsl06 = smoothMPC_ccrhs + 96;
int smoothMPC_ubIdx06[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub06 = smoothMPC_l + 107;
smoothMPC_FLOAT* smoothMPC_sub06 = smoothMPC_s + 107;
smoothMPC_FLOAT* smoothMPC_lubbysub06 = smoothMPC_lbys + 107;
smoothMPC_FLOAT smoothMPC_riub06[5];
smoothMPC_FLOAT* smoothMPC_dlubaff06 = smoothMPC_dl_aff + 107;
smoothMPC_FLOAT* smoothMPC_dsubaff06 = smoothMPC_ds_aff + 107;
smoothMPC_FLOAT* smoothMPC_dlubcc06 = smoothMPC_dl_cc + 107;
smoothMPC_FLOAT* smoothMPC_dsubcc06 = smoothMPC_ds_cc + 107;
smoothMPC_FLOAT* smoothMPC_ccrhsub06 = smoothMPC_ccrhs + 107;
smoothMPC_FLOAT smoothMPC_Phi06[11];
smoothMPC_FLOAT smoothMPC_W06[11];
smoothMPC_FLOAT smoothMPC_Ysd06[9];
smoothMPC_FLOAT smoothMPC_Lsd06[9];
smoothMPC_FLOAT* smoothMPC_z07 = smoothMPC_z + 77;
smoothMPC_FLOAT* smoothMPC_dzaff07 = smoothMPC_dz_aff + 77;
smoothMPC_FLOAT* smoothMPC_dzcc07 = smoothMPC_dz_cc + 77;
smoothMPC_FLOAT* smoothMPC_rd07 = smoothMPC_rd + 77;
smoothMPC_FLOAT smoothMPC_Lbyrd07[11];
smoothMPC_FLOAT* smoothMPC_grad_cost07 = smoothMPC_grad_cost + 77;
smoothMPC_FLOAT* smoothMPC_grad_eq07 = smoothMPC_grad_eq + 77;
smoothMPC_FLOAT* smoothMPC_grad_ineq07 = smoothMPC_grad_ineq + 77;
smoothMPC_FLOAT smoothMPC_ctv07[11];
smoothMPC_FLOAT* smoothMPC_v07 = smoothMPC_v + 21;
smoothMPC_FLOAT smoothMPC_re07[3];
smoothMPC_FLOAT smoothMPC_beta07[3];
smoothMPC_FLOAT smoothMPC_betacc07[3];
smoothMPC_FLOAT* smoothMPC_dvaff07 = smoothMPC_dv_aff + 21;
smoothMPC_FLOAT* smoothMPC_dvcc07 = smoothMPC_dv_cc + 21;
smoothMPC_FLOAT smoothMPC_V07[33];
smoothMPC_FLOAT smoothMPC_Yd07[6];
smoothMPC_FLOAT smoothMPC_Ld07[6];
smoothMPC_FLOAT smoothMPC_yy07[3];
smoothMPC_FLOAT smoothMPC_bmy07[3];
int smoothMPC_lbIdx07[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb07 = smoothMPC_l + 112;
smoothMPC_FLOAT* smoothMPC_slb07 = smoothMPC_s + 112;
smoothMPC_FLOAT* smoothMPC_llbbyslb07 = smoothMPC_lbys + 112;
smoothMPC_FLOAT smoothMPC_rilb07[11];
smoothMPC_FLOAT* smoothMPC_dllbaff07 = smoothMPC_dl_aff + 112;
smoothMPC_FLOAT* smoothMPC_dslbaff07 = smoothMPC_ds_aff + 112;
smoothMPC_FLOAT* smoothMPC_dllbcc07 = smoothMPC_dl_cc + 112;
smoothMPC_FLOAT* smoothMPC_dslbcc07 = smoothMPC_ds_cc + 112;
smoothMPC_FLOAT* smoothMPC_ccrhsl07 = smoothMPC_ccrhs + 112;
int smoothMPC_ubIdx07[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub07 = smoothMPC_l + 123;
smoothMPC_FLOAT* smoothMPC_sub07 = smoothMPC_s + 123;
smoothMPC_FLOAT* smoothMPC_lubbysub07 = smoothMPC_lbys + 123;
smoothMPC_FLOAT smoothMPC_riub07[5];
smoothMPC_FLOAT* smoothMPC_dlubaff07 = smoothMPC_dl_aff + 123;
smoothMPC_FLOAT* smoothMPC_dsubaff07 = smoothMPC_ds_aff + 123;
smoothMPC_FLOAT* smoothMPC_dlubcc07 = smoothMPC_dl_cc + 123;
smoothMPC_FLOAT* smoothMPC_dsubcc07 = smoothMPC_ds_cc + 123;
smoothMPC_FLOAT* smoothMPC_ccrhsub07 = smoothMPC_ccrhs + 123;
smoothMPC_FLOAT smoothMPC_Phi07[11];
smoothMPC_FLOAT smoothMPC_W07[11];
smoothMPC_FLOAT smoothMPC_Ysd07[9];
smoothMPC_FLOAT smoothMPC_Lsd07[9];
smoothMPC_FLOAT* smoothMPC_z08 = smoothMPC_z + 88;
smoothMPC_FLOAT* smoothMPC_dzaff08 = smoothMPC_dz_aff + 88;
smoothMPC_FLOAT* smoothMPC_dzcc08 = smoothMPC_dz_cc + 88;
smoothMPC_FLOAT* smoothMPC_rd08 = smoothMPC_rd + 88;
smoothMPC_FLOAT smoothMPC_Lbyrd08[11];
smoothMPC_FLOAT* smoothMPC_grad_cost08 = smoothMPC_grad_cost + 88;
smoothMPC_FLOAT* smoothMPC_grad_eq08 = smoothMPC_grad_eq + 88;
smoothMPC_FLOAT* smoothMPC_grad_ineq08 = smoothMPC_grad_ineq + 88;
smoothMPC_FLOAT smoothMPC_ctv08[11];
smoothMPC_FLOAT* smoothMPC_v08 = smoothMPC_v + 24;
smoothMPC_FLOAT smoothMPC_re08[3];
smoothMPC_FLOAT smoothMPC_beta08[3];
smoothMPC_FLOAT smoothMPC_betacc08[3];
smoothMPC_FLOAT* smoothMPC_dvaff08 = smoothMPC_dv_aff + 24;
smoothMPC_FLOAT* smoothMPC_dvcc08 = smoothMPC_dv_cc + 24;
smoothMPC_FLOAT smoothMPC_V08[33];
smoothMPC_FLOAT smoothMPC_Yd08[6];
smoothMPC_FLOAT smoothMPC_Ld08[6];
smoothMPC_FLOAT smoothMPC_yy08[3];
smoothMPC_FLOAT smoothMPC_bmy08[3];
int smoothMPC_lbIdx08[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb08 = smoothMPC_l + 128;
smoothMPC_FLOAT* smoothMPC_slb08 = smoothMPC_s + 128;
smoothMPC_FLOAT* smoothMPC_llbbyslb08 = smoothMPC_lbys + 128;
smoothMPC_FLOAT smoothMPC_rilb08[11];
smoothMPC_FLOAT* smoothMPC_dllbaff08 = smoothMPC_dl_aff + 128;
smoothMPC_FLOAT* smoothMPC_dslbaff08 = smoothMPC_ds_aff + 128;
smoothMPC_FLOAT* smoothMPC_dllbcc08 = smoothMPC_dl_cc + 128;
smoothMPC_FLOAT* smoothMPC_dslbcc08 = smoothMPC_ds_cc + 128;
smoothMPC_FLOAT* smoothMPC_ccrhsl08 = smoothMPC_ccrhs + 128;
int smoothMPC_ubIdx08[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub08 = smoothMPC_l + 139;
smoothMPC_FLOAT* smoothMPC_sub08 = smoothMPC_s + 139;
smoothMPC_FLOAT* smoothMPC_lubbysub08 = smoothMPC_lbys + 139;
smoothMPC_FLOAT smoothMPC_riub08[5];
smoothMPC_FLOAT* smoothMPC_dlubaff08 = smoothMPC_dl_aff + 139;
smoothMPC_FLOAT* smoothMPC_dsubaff08 = smoothMPC_ds_aff + 139;
smoothMPC_FLOAT* smoothMPC_dlubcc08 = smoothMPC_dl_cc + 139;
smoothMPC_FLOAT* smoothMPC_dsubcc08 = smoothMPC_ds_cc + 139;
smoothMPC_FLOAT* smoothMPC_ccrhsub08 = smoothMPC_ccrhs + 139;
smoothMPC_FLOAT smoothMPC_Phi08[11];
smoothMPC_FLOAT smoothMPC_W08[11];
smoothMPC_FLOAT smoothMPC_Ysd08[9];
smoothMPC_FLOAT smoothMPC_Lsd08[9];
smoothMPC_FLOAT* smoothMPC_z09 = smoothMPC_z + 99;
smoothMPC_FLOAT* smoothMPC_dzaff09 = smoothMPC_dz_aff + 99;
smoothMPC_FLOAT* smoothMPC_dzcc09 = smoothMPC_dz_cc + 99;
smoothMPC_FLOAT* smoothMPC_rd09 = smoothMPC_rd + 99;
smoothMPC_FLOAT smoothMPC_Lbyrd09[11];
smoothMPC_FLOAT* smoothMPC_grad_cost09 = smoothMPC_grad_cost + 99;
smoothMPC_FLOAT* smoothMPC_grad_eq09 = smoothMPC_grad_eq + 99;
smoothMPC_FLOAT* smoothMPC_grad_ineq09 = smoothMPC_grad_ineq + 99;
smoothMPC_FLOAT smoothMPC_ctv09[11];
smoothMPC_FLOAT* smoothMPC_v09 = smoothMPC_v + 27;
smoothMPC_FLOAT smoothMPC_re09[3];
smoothMPC_FLOAT smoothMPC_beta09[3];
smoothMPC_FLOAT smoothMPC_betacc09[3];
smoothMPC_FLOAT* smoothMPC_dvaff09 = smoothMPC_dv_aff + 27;
smoothMPC_FLOAT* smoothMPC_dvcc09 = smoothMPC_dv_cc + 27;
smoothMPC_FLOAT smoothMPC_V09[33];
smoothMPC_FLOAT smoothMPC_Yd09[6];
smoothMPC_FLOAT smoothMPC_Ld09[6];
smoothMPC_FLOAT smoothMPC_yy09[3];
smoothMPC_FLOAT smoothMPC_bmy09[3];
int smoothMPC_lbIdx09[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb09 = smoothMPC_l + 144;
smoothMPC_FLOAT* smoothMPC_slb09 = smoothMPC_s + 144;
smoothMPC_FLOAT* smoothMPC_llbbyslb09 = smoothMPC_lbys + 144;
smoothMPC_FLOAT smoothMPC_rilb09[11];
smoothMPC_FLOAT* smoothMPC_dllbaff09 = smoothMPC_dl_aff + 144;
smoothMPC_FLOAT* smoothMPC_dslbaff09 = smoothMPC_ds_aff + 144;
smoothMPC_FLOAT* smoothMPC_dllbcc09 = smoothMPC_dl_cc + 144;
smoothMPC_FLOAT* smoothMPC_dslbcc09 = smoothMPC_ds_cc + 144;
smoothMPC_FLOAT* smoothMPC_ccrhsl09 = smoothMPC_ccrhs + 144;
int smoothMPC_ubIdx09[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub09 = smoothMPC_l + 155;
smoothMPC_FLOAT* smoothMPC_sub09 = smoothMPC_s + 155;
smoothMPC_FLOAT* smoothMPC_lubbysub09 = smoothMPC_lbys + 155;
smoothMPC_FLOAT smoothMPC_riub09[5];
smoothMPC_FLOAT* smoothMPC_dlubaff09 = smoothMPC_dl_aff + 155;
smoothMPC_FLOAT* smoothMPC_dsubaff09 = smoothMPC_ds_aff + 155;
smoothMPC_FLOAT* smoothMPC_dlubcc09 = smoothMPC_dl_cc + 155;
smoothMPC_FLOAT* smoothMPC_dsubcc09 = smoothMPC_ds_cc + 155;
smoothMPC_FLOAT* smoothMPC_ccrhsub09 = smoothMPC_ccrhs + 155;
smoothMPC_FLOAT smoothMPC_Phi09[11];
smoothMPC_FLOAT smoothMPC_W09[11];
smoothMPC_FLOAT smoothMPC_Ysd09[9];
smoothMPC_FLOAT smoothMPC_Lsd09[9];
smoothMPC_FLOAT* smoothMPC_z10 = smoothMPC_z + 110;
smoothMPC_FLOAT* smoothMPC_dzaff10 = smoothMPC_dz_aff + 110;
smoothMPC_FLOAT* smoothMPC_dzcc10 = smoothMPC_dz_cc + 110;
smoothMPC_FLOAT* smoothMPC_rd10 = smoothMPC_rd + 110;
smoothMPC_FLOAT smoothMPC_Lbyrd10[11];
smoothMPC_FLOAT* smoothMPC_grad_cost10 = smoothMPC_grad_cost + 110;
smoothMPC_FLOAT* smoothMPC_grad_eq10 = smoothMPC_grad_eq + 110;
smoothMPC_FLOAT* smoothMPC_grad_ineq10 = smoothMPC_grad_ineq + 110;
smoothMPC_FLOAT smoothMPC_ctv10[11];
smoothMPC_FLOAT* smoothMPC_v10 = smoothMPC_v + 30;
smoothMPC_FLOAT smoothMPC_re10[3];
smoothMPC_FLOAT smoothMPC_beta10[3];
smoothMPC_FLOAT smoothMPC_betacc10[3];
smoothMPC_FLOAT* smoothMPC_dvaff10 = smoothMPC_dv_aff + 30;
smoothMPC_FLOAT* smoothMPC_dvcc10 = smoothMPC_dv_cc + 30;
smoothMPC_FLOAT smoothMPC_V10[33];
smoothMPC_FLOAT smoothMPC_Yd10[6];
smoothMPC_FLOAT smoothMPC_Ld10[6];
smoothMPC_FLOAT smoothMPC_yy10[3];
smoothMPC_FLOAT smoothMPC_bmy10[3];
int smoothMPC_lbIdx10[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb10 = smoothMPC_l + 160;
smoothMPC_FLOAT* smoothMPC_slb10 = smoothMPC_s + 160;
smoothMPC_FLOAT* smoothMPC_llbbyslb10 = smoothMPC_lbys + 160;
smoothMPC_FLOAT smoothMPC_rilb10[11];
smoothMPC_FLOAT* smoothMPC_dllbaff10 = smoothMPC_dl_aff + 160;
smoothMPC_FLOAT* smoothMPC_dslbaff10 = smoothMPC_ds_aff + 160;
smoothMPC_FLOAT* smoothMPC_dllbcc10 = smoothMPC_dl_cc + 160;
smoothMPC_FLOAT* smoothMPC_dslbcc10 = smoothMPC_ds_cc + 160;
smoothMPC_FLOAT* smoothMPC_ccrhsl10 = smoothMPC_ccrhs + 160;
int smoothMPC_ubIdx10[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub10 = smoothMPC_l + 171;
smoothMPC_FLOAT* smoothMPC_sub10 = smoothMPC_s + 171;
smoothMPC_FLOAT* smoothMPC_lubbysub10 = smoothMPC_lbys + 171;
smoothMPC_FLOAT smoothMPC_riub10[5];
smoothMPC_FLOAT* smoothMPC_dlubaff10 = smoothMPC_dl_aff + 171;
smoothMPC_FLOAT* smoothMPC_dsubaff10 = smoothMPC_ds_aff + 171;
smoothMPC_FLOAT* smoothMPC_dlubcc10 = smoothMPC_dl_cc + 171;
smoothMPC_FLOAT* smoothMPC_dsubcc10 = smoothMPC_ds_cc + 171;
smoothMPC_FLOAT* smoothMPC_ccrhsub10 = smoothMPC_ccrhs + 171;
smoothMPC_FLOAT smoothMPC_Phi10[11];
smoothMPC_FLOAT smoothMPC_W10[11];
smoothMPC_FLOAT smoothMPC_Ysd10[9];
smoothMPC_FLOAT smoothMPC_Lsd10[9];
smoothMPC_FLOAT* smoothMPC_z11 = smoothMPC_z + 121;
smoothMPC_FLOAT* smoothMPC_dzaff11 = smoothMPC_dz_aff + 121;
smoothMPC_FLOAT* smoothMPC_dzcc11 = smoothMPC_dz_cc + 121;
smoothMPC_FLOAT* smoothMPC_rd11 = smoothMPC_rd + 121;
smoothMPC_FLOAT smoothMPC_Lbyrd11[11];
smoothMPC_FLOAT* smoothMPC_grad_cost11 = smoothMPC_grad_cost + 121;
smoothMPC_FLOAT* smoothMPC_grad_eq11 = smoothMPC_grad_eq + 121;
smoothMPC_FLOAT* smoothMPC_grad_ineq11 = smoothMPC_grad_ineq + 121;
smoothMPC_FLOAT smoothMPC_ctv11[11];
smoothMPC_FLOAT* smoothMPC_v11 = smoothMPC_v + 33;
smoothMPC_FLOAT smoothMPC_re11[3];
smoothMPC_FLOAT smoothMPC_beta11[3];
smoothMPC_FLOAT smoothMPC_betacc11[3];
smoothMPC_FLOAT* smoothMPC_dvaff11 = smoothMPC_dv_aff + 33;
smoothMPC_FLOAT* smoothMPC_dvcc11 = smoothMPC_dv_cc + 33;
smoothMPC_FLOAT smoothMPC_V11[33];
smoothMPC_FLOAT smoothMPC_Yd11[6];
smoothMPC_FLOAT smoothMPC_Ld11[6];
smoothMPC_FLOAT smoothMPC_yy11[3];
smoothMPC_FLOAT smoothMPC_bmy11[3];
int smoothMPC_lbIdx11[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb11 = smoothMPC_l + 176;
smoothMPC_FLOAT* smoothMPC_slb11 = smoothMPC_s + 176;
smoothMPC_FLOAT* smoothMPC_llbbyslb11 = smoothMPC_lbys + 176;
smoothMPC_FLOAT smoothMPC_rilb11[11];
smoothMPC_FLOAT* smoothMPC_dllbaff11 = smoothMPC_dl_aff + 176;
smoothMPC_FLOAT* smoothMPC_dslbaff11 = smoothMPC_ds_aff + 176;
smoothMPC_FLOAT* smoothMPC_dllbcc11 = smoothMPC_dl_cc + 176;
smoothMPC_FLOAT* smoothMPC_dslbcc11 = smoothMPC_ds_cc + 176;
smoothMPC_FLOAT* smoothMPC_ccrhsl11 = smoothMPC_ccrhs + 176;
int smoothMPC_ubIdx11[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub11 = smoothMPC_l + 187;
smoothMPC_FLOAT* smoothMPC_sub11 = smoothMPC_s + 187;
smoothMPC_FLOAT* smoothMPC_lubbysub11 = smoothMPC_lbys + 187;
smoothMPC_FLOAT smoothMPC_riub11[5];
smoothMPC_FLOAT* smoothMPC_dlubaff11 = smoothMPC_dl_aff + 187;
smoothMPC_FLOAT* smoothMPC_dsubaff11 = smoothMPC_ds_aff + 187;
smoothMPC_FLOAT* smoothMPC_dlubcc11 = smoothMPC_dl_cc + 187;
smoothMPC_FLOAT* smoothMPC_dsubcc11 = smoothMPC_ds_cc + 187;
smoothMPC_FLOAT* smoothMPC_ccrhsub11 = smoothMPC_ccrhs + 187;
smoothMPC_FLOAT smoothMPC_Phi11[11];
smoothMPC_FLOAT smoothMPC_W11[11];
smoothMPC_FLOAT smoothMPC_Ysd11[9];
smoothMPC_FLOAT smoothMPC_Lsd11[9];
smoothMPC_FLOAT* smoothMPC_z12 = smoothMPC_z + 132;
smoothMPC_FLOAT* smoothMPC_dzaff12 = smoothMPC_dz_aff + 132;
smoothMPC_FLOAT* smoothMPC_dzcc12 = smoothMPC_dz_cc + 132;
smoothMPC_FLOAT* smoothMPC_rd12 = smoothMPC_rd + 132;
smoothMPC_FLOAT smoothMPC_Lbyrd12[11];
smoothMPC_FLOAT* smoothMPC_grad_cost12 = smoothMPC_grad_cost + 132;
smoothMPC_FLOAT* smoothMPC_grad_eq12 = smoothMPC_grad_eq + 132;
smoothMPC_FLOAT* smoothMPC_grad_ineq12 = smoothMPC_grad_ineq + 132;
smoothMPC_FLOAT smoothMPC_ctv12[11];
smoothMPC_FLOAT* smoothMPC_v12 = smoothMPC_v + 36;
smoothMPC_FLOAT smoothMPC_re12[3];
smoothMPC_FLOAT smoothMPC_beta12[3];
smoothMPC_FLOAT smoothMPC_betacc12[3];
smoothMPC_FLOAT* smoothMPC_dvaff12 = smoothMPC_dv_aff + 36;
smoothMPC_FLOAT* smoothMPC_dvcc12 = smoothMPC_dv_cc + 36;
smoothMPC_FLOAT smoothMPC_V12[33];
smoothMPC_FLOAT smoothMPC_Yd12[6];
smoothMPC_FLOAT smoothMPC_Ld12[6];
smoothMPC_FLOAT smoothMPC_yy12[3];
smoothMPC_FLOAT smoothMPC_bmy12[3];
int smoothMPC_lbIdx12[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb12 = smoothMPC_l + 192;
smoothMPC_FLOAT* smoothMPC_slb12 = smoothMPC_s + 192;
smoothMPC_FLOAT* smoothMPC_llbbyslb12 = smoothMPC_lbys + 192;
smoothMPC_FLOAT smoothMPC_rilb12[11];
smoothMPC_FLOAT* smoothMPC_dllbaff12 = smoothMPC_dl_aff + 192;
smoothMPC_FLOAT* smoothMPC_dslbaff12 = smoothMPC_ds_aff + 192;
smoothMPC_FLOAT* smoothMPC_dllbcc12 = smoothMPC_dl_cc + 192;
smoothMPC_FLOAT* smoothMPC_dslbcc12 = smoothMPC_ds_cc + 192;
smoothMPC_FLOAT* smoothMPC_ccrhsl12 = smoothMPC_ccrhs + 192;
int smoothMPC_ubIdx12[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub12 = smoothMPC_l + 203;
smoothMPC_FLOAT* smoothMPC_sub12 = smoothMPC_s + 203;
smoothMPC_FLOAT* smoothMPC_lubbysub12 = smoothMPC_lbys + 203;
smoothMPC_FLOAT smoothMPC_riub12[5];
smoothMPC_FLOAT* smoothMPC_dlubaff12 = smoothMPC_dl_aff + 203;
smoothMPC_FLOAT* smoothMPC_dsubaff12 = smoothMPC_ds_aff + 203;
smoothMPC_FLOAT* smoothMPC_dlubcc12 = smoothMPC_dl_cc + 203;
smoothMPC_FLOAT* smoothMPC_dsubcc12 = smoothMPC_ds_cc + 203;
smoothMPC_FLOAT* smoothMPC_ccrhsub12 = smoothMPC_ccrhs + 203;
smoothMPC_FLOAT smoothMPC_Phi12[11];
smoothMPC_FLOAT smoothMPC_W12[11];
smoothMPC_FLOAT smoothMPC_Ysd12[9];
smoothMPC_FLOAT smoothMPC_Lsd12[9];
smoothMPC_FLOAT* smoothMPC_z13 = smoothMPC_z + 143;
smoothMPC_FLOAT* smoothMPC_dzaff13 = smoothMPC_dz_aff + 143;
smoothMPC_FLOAT* smoothMPC_dzcc13 = smoothMPC_dz_cc + 143;
smoothMPC_FLOAT* smoothMPC_rd13 = smoothMPC_rd + 143;
smoothMPC_FLOAT smoothMPC_Lbyrd13[11];
smoothMPC_FLOAT* smoothMPC_grad_cost13 = smoothMPC_grad_cost + 143;
smoothMPC_FLOAT* smoothMPC_grad_eq13 = smoothMPC_grad_eq + 143;
smoothMPC_FLOAT* smoothMPC_grad_ineq13 = smoothMPC_grad_ineq + 143;
smoothMPC_FLOAT smoothMPC_ctv13[11];
smoothMPC_FLOAT* smoothMPC_v13 = smoothMPC_v + 39;
smoothMPC_FLOAT smoothMPC_re13[3];
smoothMPC_FLOAT smoothMPC_beta13[3];
smoothMPC_FLOAT smoothMPC_betacc13[3];
smoothMPC_FLOAT* smoothMPC_dvaff13 = smoothMPC_dv_aff + 39;
smoothMPC_FLOAT* smoothMPC_dvcc13 = smoothMPC_dv_cc + 39;
smoothMPC_FLOAT smoothMPC_V13[33];
smoothMPC_FLOAT smoothMPC_Yd13[6];
smoothMPC_FLOAT smoothMPC_Ld13[6];
smoothMPC_FLOAT smoothMPC_yy13[3];
smoothMPC_FLOAT smoothMPC_bmy13[3];
int smoothMPC_lbIdx13[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb13 = smoothMPC_l + 208;
smoothMPC_FLOAT* smoothMPC_slb13 = smoothMPC_s + 208;
smoothMPC_FLOAT* smoothMPC_llbbyslb13 = smoothMPC_lbys + 208;
smoothMPC_FLOAT smoothMPC_rilb13[11];
smoothMPC_FLOAT* smoothMPC_dllbaff13 = smoothMPC_dl_aff + 208;
smoothMPC_FLOAT* smoothMPC_dslbaff13 = smoothMPC_ds_aff + 208;
smoothMPC_FLOAT* smoothMPC_dllbcc13 = smoothMPC_dl_cc + 208;
smoothMPC_FLOAT* smoothMPC_dslbcc13 = smoothMPC_ds_cc + 208;
smoothMPC_FLOAT* smoothMPC_ccrhsl13 = smoothMPC_ccrhs + 208;
int smoothMPC_ubIdx13[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub13 = smoothMPC_l + 219;
smoothMPC_FLOAT* smoothMPC_sub13 = smoothMPC_s + 219;
smoothMPC_FLOAT* smoothMPC_lubbysub13 = smoothMPC_lbys + 219;
smoothMPC_FLOAT smoothMPC_riub13[5];
smoothMPC_FLOAT* smoothMPC_dlubaff13 = smoothMPC_dl_aff + 219;
smoothMPC_FLOAT* smoothMPC_dsubaff13 = smoothMPC_ds_aff + 219;
smoothMPC_FLOAT* smoothMPC_dlubcc13 = smoothMPC_dl_cc + 219;
smoothMPC_FLOAT* smoothMPC_dsubcc13 = smoothMPC_ds_cc + 219;
smoothMPC_FLOAT* smoothMPC_ccrhsub13 = smoothMPC_ccrhs + 219;
smoothMPC_FLOAT smoothMPC_Phi13[11];
smoothMPC_FLOAT smoothMPC_W13[11];
smoothMPC_FLOAT smoothMPC_Ysd13[9];
smoothMPC_FLOAT smoothMPC_Lsd13[9];
smoothMPC_FLOAT* smoothMPC_z14 = smoothMPC_z + 154;
smoothMPC_FLOAT* smoothMPC_dzaff14 = smoothMPC_dz_aff + 154;
smoothMPC_FLOAT* smoothMPC_dzcc14 = smoothMPC_dz_cc + 154;
smoothMPC_FLOAT* smoothMPC_rd14 = smoothMPC_rd + 154;
smoothMPC_FLOAT smoothMPC_Lbyrd14[11];
smoothMPC_FLOAT* smoothMPC_grad_cost14 = smoothMPC_grad_cost + 154;
smoothMPC_FLOAT* smoothMPC_grad_eq14 = smoothMPC_grad_eq + 154;
smoothMPC_FLOAT* smoothMPC_grad_ineq14 = smoothMPC_grad_ineq + 154;
smoothMPC_FLOAT smoothMPC_ctv14[11];
smoothMPC_FLOAT* smoothMPC_v14 = smoothMPC_v + 42;
smoothMPC_FLOAT smoothMPC_re14[3];
smoothMPC_FLOAT smoothMPC_beta14[3];
smoothMPC_FLOAT smoothMPC_betacc14[3];
smoothMPC_FLOAT* smoothMPC_dvaff14 = smoothMPC_dv_aff + 42;
smoothMPC_FLOAT* smoothMPC_dvcc14 = smoothMPC_dv_cc + 42;
smoothMPC_FLOAT smoothMPC_V14[33];
smoothMPC_FLOAT smoothMPC_Yd14[6];
smoothMPC_FLOAT smoothMPC_Ld14[6];
smoothMPC_FLOAT smoothMPC_yy14[3];
smoothMPC_FLOAT smoothMPC_bmy14[3];
int smoothMPC_lbIdx14[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb14 = smoothMPC_l + 224;
smoothMPC_FLOAT* smoothMPC_slb14 = smoothMPC_s + 224;
smoothMPC_FLOAT* smoothMPC_llbbyslb14 = smoothMPC_lbys + 224;
smoothMPC_FLOAT smoothMPC_rilb14[11];
smoothMPC_FLOAT* smoothMPC_dllbaff14 = smoothMPC_dl_aff + 224;
smoothMPC_FLOAT* smoothMPC_dslbaff14 = smoothMPC_ds_aff + 224;
smoothMPC_FLOAT* smoothMPC_dllbcc14 = smoothMPC_dl_cc + 224;
smoothMPC_FLOAT* smoothMPC_dslbcc14 = smoothMPC_ds_cc + 224;
smoothMPC_FLOAT* smoothMPC_ccrhsl14 = smoothMPC_ccrhs + 224;
int smoothMPC_ubIdx14[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub14 = smoothMPC_l + 235;
smoothMPC_FLOAT* smoothMPC_sub14 = smoothMPC_s + 235;
smoothMPC_FLOAT* smoothMPC_lubbysub14 = smoothMPC_lbys + 235;
smoothMPC_FLOAT smoothMPC_riub14[5];
smoothMPC_FLOAT* smoothMPC_dlubaff14 = smoothMPC_dl_aff + 235;
smoothMPC_FLOAT* smoothMPC_dsubaff14 = smoothMPC_ds_aff + 235;
smoothMPC_FLOAT* smoothMPC_dlubcc14 = smoothMPC_dl_cc + 235;
smoothMPC_FLOAT* smoothMPC_dsubcc14 = smoothMPC_ds_cc + 235;
smoothMPC_FLOAT* smoothMPC_ccrhsub14 = smoothMPC_ccrhs + 235;
smoothMPC_FLOAT smoothMPC_Phi14[11];
smoothMPC_FLOAT smoothMPC_W14[11];
smoothMPC_FLOAT smoothMPC_Ysd14[9];
smoothMPC_FLOAT smoothMPC_Lsd14[9];
smoothMPC_FLOAT* smoothMPC_z15 = smoothMPC_z + 165;
smoothMPC_FLOAT* smoothMPC_dzaff15 = smoothMPC_dz_aff + 165;
smoothMPC_FLOAT* smoothMPC_dzcc15 = smoothMPC_dz_cc + 165;
smoothMPC_FLOAT* smoothMPC_rd15 = smoothMPC_rd + 165;
smoothMPC_FLOAT smoothMPC_Lbyrd15[11];
smoothMPC_FLOAT* smoothMPC_grad_cost15 = smoothMPC_grad_cost + 165;
smoothMPC_FLOAT* smoothMPC_grad_eq15 = smoothMPC_grad_eq + 165;
smoothMPC_FLOAT* smoothMPC_grad_ineq15 = smoothMPC_grad_ineq + 165;
smoothMPC_FLOAT smoothMPC_ctv15[11];
smoothMPC_FLOAT* smoothMPC_v15 = smoothMPC_v + 45;
smoothMPC_FLOAT smoothMPC_re15[3];
smoothMPC_FLOAT smoothMPC_beta15[3];
smoothMPC_FLOAT smoothMPC_betacc15[3];
smoothMPC_FLOAT* smoothMPC_dvaff15 = smoothMPC_dv_aff + 45;
smoothMPC_FLOAT* smoothMPC_dvcc15 = smoothMPC_dv_cc + 45;
smoothMPC_FLOAT smoothMPC_V15[33];
smoothMPC_FLOAT smoothMPC_Yd15[6];
smoothMPC_FLOAT smoothMPC_Ld15[6];
smoothMPC_FLOAT smoothMPC_yy15[3];
smoothMPC_FLOAT smoothMPC_bmy15[3];
int smoothMPC_lbIdx15[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb15 = smoothMPC_l + 240;
smoothMPC_FLOAT* smoothMPC_slb15 = smoothMPC_s + 240;
smoothMPC_FLOAT* smoothMPC_llbbyslb15 = smoothMPC_lbys + 240;
smoothMPC_FLOAT smoothMPC_rilb15[11];
smoothMPC_FLOAT* smoothMPC_dllbaff15 = smoothMPC_dl_aff + 240;
smoothMPC_FLOAT* smoothMPC_dslbaff15 = smoothMPC_ds_aff + 240;
smoothMPC_FLOAT* smoothMPC_dllbcc15 = smoothMPC_dl_cc + 240;
smoothMPC_FLOAT* smoothMPC_dslbcc15 = smoothMPC_ds_cc + 240;
smoothMPC_FLOAT* smoothMPC_ccrhsl15 = smoothMPC_ccrhs + 240;
int smoothMPC_ubIdx15[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub15 = smoothMPC_l + 251;
smoothMPC_FLOAT* smoothMPC_sub15 = smoothMPC_s + 251;
smoothMPC_FLOAT* smoothMPC_lubbysub15 = smoothMPC_lbys + 251;
smoothMPC_FLOAT smoothMPC_riub15[5];
smoothMPC_FLOAT* smoothMPC_dlubaff15 = smoothMPC_dl_aff + 251;
smoothMPC_FLOAT* smoothMPC_dsubaff15 = smoothMPC_ds_aff + 251;
smoothMPC_FLOAT* smoothMPC_dlubcc15 = smoothMPC_dl_cc + 251;
smoothMPC_FLOAT* smoothMPC_dsubcc15 = smoothMPC_ds_cc + 251;
smoothMPC_FLOAT* smoothMPC_ccrhsub15 = smoothMPC_ccrhs + 251;
smoothMPC_FLOAT smoothMPC_Phi15[11];
smoothMPC_FLOAT smoothMPC_W15[11];
smoothMPC_FLOAT smoothMPC_Ysd15[9];
smoothMPC_FLOAT smoothMPC_Lsd15[9];
smoothMPC_FLOAT* smoothMPC_z16 = smoothMPC_z + 176;
smoothMPC_FLOAT* smoothMPC_dzaff16 = smoothMPC_dz_aff + 176;
smoothMPC_FLOAT* smoothMPC_dzcc16 = smoothMPC_dz_cc + 176;
smoothMPC_FLOAT* smoothMPC_rd16 = smoothMPC_rd + 176;
smoothMPC_FLOAT smoothMPC_Lbyrd16[11];
smoothMPC_FLOAT* smoothMPC_grad_cost16 = smoothMPC_grad_cost + 176;
smoothMPC_FLOAT* smoothMPC_grad_eq16 = smoothMPC_grad_eq + 176;
smoothMPC_FLOAT* smoothMPC_grad_ineq16 = smoothMPC_grad_ineq + 176;
smoothMPC_FLOAT smoothMPC_ctv16[11];
smoothMPC_FLOAT* smoothMPC_v16 = smoothMPC_v + 48;
smoothMPC_FLOAT smoothMPC_re16[3];
smoothMPC_FLOAT smoothMPC_beta16[3];
smoothMPC_FLOAT smoothMPC_betacc16[3];
smoothMPC_FLOAT* smoothMPC_dvaff16 = smoothMPC_dv_aff + 48;
smoothMPC_FLOAT* smoothMPC_dvcc16 = smoothMPC_dv_cc + 48;
smoothMPC_FLOAT smoothMPC_V16[33];
smoothMPC_FLOAT smoothMPC_Yd16[6];
smoothMPC_FLOAT smoothMPC_Ld16[6];
smoothMPC_FLOAT smoothMPC_yy16[3];
smoothMPC_FLOAT smoothMPC_bmy16[3];
int smoothMPC_lbIdx16[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb16 = smoothMPC_l + 256;
smoothMPC_FLOAT* smoothMPC_slb16 = smoothMPC_s + 256;
smoothMPC_FLOAT* smoothMPC_llbbyslb16 = smoothMPC_lbys + 256;
smoothMPC_FLOAT smoothMPC_rilb16[11];
smoothMPC_FLOAT* smoothMPC_dllbaff16 = smoothMPC_dl_aff + 256;
smoothMPC_FLOAT* smoothMPC_dslbaff16 = smoothMPC_ds_aff + 256;
smoothMPC_FLOAT* smoothMPC_dllbcc16 = smoothMPC_dl_cc + 256;
smoothMPC_FLOAT* smoothMPC_dslbcc16 = smoothMPC_ds_cc + 256;
smoothMPC_FLOAT* smoothMPC_ccrhsl16 = smoothMPC_ccrhs + 256;
int smoothMPC_ubIdx16[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub16 = smoothMPC_l + 267;
smoothMPC_FLOAT* smoothMPC_sub16 = smoothMPC_s + 267;
smoothMPC_FLOAT* smoothMPC_lubbysub16 = smoothMPC_lbys + 267;
smoothMPC_FLOAT smoothMPC_riub16[5];
smoothMPC_FLOAT* smoothMPC_dlubaff16 = smoothMPC_dl_aff + 267;
smoothMPC_FLOAT* smoothMPC_dsubaff16 = smoothMPC_ds_aff + 267;
smoothMPC_FLOAT* smoothMPC_dlubcc16 = smoothMPC_dl_cc + 267;
smoothMPC_FLOAT* smoothMPC_dsubcc16 = smoothMPC_ds_cc + 267;
smoothMPC_FLOAT* smoothMPC_ccrhsub16 = smoothMPC_ccrhs + 267;
smoothMPC_FLOAT smoothMPC_Phi16[11];
smoothMPC_FLOAT smoothMPC_W16[11];
smoothMPC_FLOAT smoothMPC_Ysd16[9];
smoothMPC_FLOAT smoothMPC_Lsd16[9];
smoothMPC_FLOAT* smoothMPC_z17 = smoothMPC_z + 187;
smoothMPC_FLOAT* smoothMPC_dzaff17 = smoothMPC_dz_aff + 187;
smoothMPC_FLOAT* smoothMPC_dzcc17 = smoothMPC_dz_cc + 187;
smoothMPC_FLOAT* smoothMPC_rd17 = smoothMPC_rd + 187;
smoothMPC_FLOAT smoothMPC_Lbyrd17[11];
smoothMPC_FLOAT* smoothMPC_grad_cost17 = smoothMPC_grad_cost + 187;
smoothMPC_FLOAT* smoothMPC_grad_eq17 = smoothMPC_grad_eq + 187;
smoothMPC_FLOAT* smoothMPC_grad_ineq17 = smoothMPC_grad_ineq + 187;
smoothMPC_FLOAT smoothMPC_ctv17[11];
smoothMPC_FLOAT* smoothMPC_v17 = smoothMPC_v + 51;
smoothMPC_FLOAT smoothMPC_re17[3];
smoothMPC_FLOAT smoothMPC_beta17[3];
smoothMPC_FLOAT smoothMPC_betacc17[3];
smoothMPC_FLOAT* smoothMPC_dvaff17 = smoothMPC_dv_aff + 51;
smoothMPC_FLOAT* smoothMPC_dvcc17 = smoothMPC_dv_cc + 51;
smoothMPC_FLOAT smoothMPC_V17[33];
smoothMPC_FLOAT smoothMPC_Yd17[6];
smoothMPC_FLOAT smoothMPC_Ld17[6];
smoothMPC_FLOAT smoothMPC_yy17[3];
smoothMPC_FLOAT smoothMPC_bmy17[3];
int smoothMPC_lbIdx17[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb17 = smoothMPC_l + 272;
smoothMPC_FLOAT* smoothMPC_slb17 = smoothMPC_s + 272;
smoothMPC_FLOAT* smoothMPC_llbbyslb17 = smoothMPC_lbys + 272;
smoothMPC_FLOAT smoothMPC_rilb17[11];
smoothMPC_FLOAT* smoothMPC_dllbaff17 = smoothMPC_dl_aff + 272;
smoothMPC_FLOAT* smoothMPC_dslbaff17 = smoothMPC_ds_aff + 272;
smoothMPC_FLOAT* smoothMPC_dllbcc17 = smoothMPC_dl_cc + 272;
smoothMPC_FLOAT* smoothMPC_dslbcc17 = smoothMPC_ds_cc + 272;
smoothMPC_FLOAT* smoothMPC_ccrhsl17 = smoothMPC_ccrhs + 272;
int smoothMPC_ubIdx17[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub17 = smoothMPC_l + 283;
smoothMPC_FLOAT* smoothMPC_sub17 = smoothMPC_s + 283;
smoothMPC_FLOAT* smoothMPC_lubbysub17 = smoothMPC_lbys + 283;
smoothMPC_FLOAT smoothMPC_riub17[5];
smoothMPC_FLOAT* smoothMPC_dlubaff17 = smoothMPC_dl_aff + 283;
smoothMPC_FLOAT* smoothMPC_dsubaff17 = smoothMPC_ds_aff + 283;
smoothMPC_FLOAT* smoothMPC_dlubcc17 = smoothMPC_dl_cc + 283;
smoothMPC_FLOAT* smoothMPC_dsubcc17 = smoothMPC_ds_cc + 283;
smoothMPC_FLOAT* smoothMPC_ccrhsub17 = smoothMPC_ccrhs + 283;
smoothMPC_FLOAT smoothMPC_Phi17[11];
smoothMPC_FLOAT smoothMPC_W17[11];
smoothMPC_FLOAT smoothMPC_Ysd17[9];
smoothMPC_FLOAT smoothMPC_Lsd17[9];
smoothMPC_FLOAT* smoothMPC_z18 = smoothMPC_z + 198;
smoothMPC_FLOAT* smoothMPC_dzaff18 = smoothMPC_dz_aff + 198;
smoothMPC_FLOAT* smoothMPC_dzcc18 = smoothMPC_dz_cc + 198;
smoothMPC_FLOAT* smoothMPC_rd18 = smoothMPC_rd + 198;
smoothMPC_FLOAT smoothMPC_Lbyrd18[11];
smoothMPC_FLOAT* smoothMPC_grad_cost18 = smoothMPC_grad_cost + 198;
smoothMPC_FLOAT* smoothMPC_grad_eq18 = smoothMPC_grad_eq + 198;
smoothMPC_FLOAT* smoothMPC_grad_ineq18 = smoothMPC_grad_ineq + 198;
smoothMPC_FLOAT smoothMPC_ctv18[11];
smoothMPC_FLOAT* smoothMPC_v18 = smoothMPC_v + 54;
smoothMPC_FLOAT smoothMPC_re18[3];
smoothMPC_FLOAT smoothMPC_beta18[3];
smoothMPC_FLOAT smoothMPC_betacc18[3];
smoothMPC_FLOAT* smoothMPC_dvaff18 = smoothMPC_dv_aff + 54;
smoothMPC_FLOAT* smoothMPC_dvcc18 = smoothMPC_dv_cc + 54;
smoothMPC_FLOAT smoothMPC_V18[33];
smoothMPC_FLOAT smoothMPC_Yd18[6];
smoothMPC_FLOAT smoothMPC_Ld18[6];
smoothMPC_FLOAT smoothMPC_yy18[3];
smoothMPC_FLOAT smoothMPC_bmy18[3];
int smoothMPC_lbIdx18[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb18 = smoothMPC_l + 288;
smoothMPC_FLOAT* smoothMPC_slb18 = smoothMPC_s + 288;
smoothMPC_FLOAT* smoothMPC_llbbyslb18 = smoothMPC_lbys + 288;
smoothMPC_FLOAT smoothMPC_rilb18[11];
smoothMPC_FLOAT* smoothMPC_dllbaff18 = smoothMPC_dl_aff + 288;
smoothMPC_FLOAT* smoothMPC_dslbaff18 = smoothMPC_ds_aff + 288;
smoothMPC_FLOAT* smoothMPC_dllbcc18 = smoothMPC_dl_cc + 288;
smoothMPC_FLOAT* smoothMPC_dslbcc18 = smoothMPC_ds_cc + 288;
smoothMPC_FLOAT* smoothMPC_ccrhsl18 = smoothMPC_ccrhs + 288;
int smoothMPC_ubIdx18[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub18 = smoothMPC_l + 299;
smoothMPC_FLOAT* smoothMPC_sub18 = smoothMPC_s + 299;
smoothMPC_FLOAT* smoothMPC_lubbysub18 = smoothMPC_lbys + 299;
smoothMPC_FLOAT smoothMPC_riub18[5];
smoothMPC_FLOAT* smoothMPC_dlubaff18 = smoothMPC_dl_aff + 299;
smoothMPC_FLOAT* smoothMPC_dsubaff18 = smoothMPC_ds_aff + 299;
smoothMPC_FLOAT* smoothMPC_dlubcc18 = smoothMPC_dl_cc + 299;
smoothMPC_FLOAT* smoothMPC_dsubcc18 = smoothMPC_ds_cc + 299;
smoothMPC_FLOAT* smoothMPC_ccrhsub18 = smoothMPC_ccrhs + 299;
smoothMPC_FLOAT smoothMPC_Phi18[11];
smoothMPC_FLOAT smoothMPC_W18[11];
smoothMPC_FLOAT smoothMPC_Ysd18[9];
smoothMPC_FLOAT smoothMPC_Lsd18[9];
smoothMPC_FLOAT* smoothMPC_z19 = smoothMPC_z + 209;
smoothMPC_FLOAT* smoothMPC_dzaff19 = smoothMPC_dz_aff + 209;
smoothMPC_FLOAT* smoothMPC_dzcc19 = smoothMPC_dz_cc + 209;
smoothMPC_FLOAT* smoothMPC_rd19 = smoothMPC_rd + 209;
smoothMPC_FLOAT smoothMPC_Lbyrd19[11];
smoothMPC_FLOAT* smoothMPC_grad_cost19 = smoothMPC_grad_cost + 209;
smoothMPC_FLOAT* smoothMPC_grad_eq19 = smoothMPC_grad_eq + 209;
smoothMPC_FLOAT* smoothMPC_grad_ineq19 = smoothMPC_grad_ineq + 209;
smoothMPC_FLOAT smoothMPC_ctv19[11];
smoothMPC_FLOAT* smoothMPC_v19 = smoothMPC_v + 57;
smoothMPC_FLOAT smoothMPC_re19[3];
smoothMPC_FLOAT smoothMPC_beta19[3];
smoothMPC_FLOAT smoothMPC_betacc19[3];
smoothMPC_FLOAT* smoothMPC_dvaff19 = smoothMPC_dv_aff + 57;
smoothMPC_FLOAT* smoothMPC_dvcc19 = smoothMPC_dv_cc + 57;
smoothMPC_FLOAT smoothMPC_V19[33];
smoothMPC_FLOAT smoothMPC_Yd19[6];
smoothMPC_FLOAT smoothMPC_Ld19[6];
smoothMPC_FLOAT smoothMPC_yy19[3];
smoothMPC_FLOAT smoothMPC_bmy19[3];
int smoothMPC_lbIdx19[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb19 = smoothMPC_l + 304;
smoothMPC_FLOAT* smoothMPC_slb19 = smoothMPC_s + 304;
smoothMPC_FLOAT* smoothMPC_llbbyslb19 = smoothMPC_lbys + 304;
smoothMPC_FLOAT smoothMPC_rilb19[11];
smoothMPC_FLOAT* smoothMPC_dllbaff19 = smoothMPC_dl_aff + 304;
smoothMPC_FLOAT* smoothMPC_dslbaff19 = smoothMPC_ds_aff + 304;
smoothMPC_FLOAT* smoothMPC_dllbcc19 = smoothMPC_dl_cc + 304;
smoothMPC_FLOAT* smoothMPC_dslbcc19 = smoothMPC_ds_cc + 304;
smoothMPC_FLOAT* smoothMPC_ccrhsl19 = smoothMPC_ccrhs + 304;
int smoothMPC_ubIdx19[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub19 = smoothMPC_l + 315;
smoothMPC_FLOAT* smoothMPC_sub19 = smoothMPC_s + 315;
smoothMPC_FLOAT* smoothMPC_lubbysub19 = smoothMPC_lbys + 315;
smoothMPC_FLOAT smoothMPC_riub19[5];
smoothMPC_FLOAT* smoothMPC_dlubaff19 = smoothMPC_dl_aff + 315;
smoothMPC_FLOAT* smoothMPC_dsubaff19 = smoothMPC_ds_aff + 315;
smoothMPC_FLOAT* smoothMPC_dlubcc19 = smoothMPC_dl_cc + 315;
smoothMPC_FLOAT* smoothMPC_dsubcc19 = smoothMPC_ds_cc + 315;
smoothMPC_FLOAT* smoothMPC_ccrhsub19 = smoothMPC_ccrhs + 315;
smoothMPC_FLOAT smoothMPC_Phi19[11];
smoothMPC_FLOAT smoothMPC_W19[11];
smoothMPC_FLOAT smoothMPC_Ysd19[9];
smoothMPC_FLOAT smoothMPC_Lsd19[9];
smoothMPC_FLOAT* smoothMPC_z20 = smoothMPC_z + 220;
smoothMPC_FLOAT* smoothMPC_dzaff20 = smoothMPC_dz_aff + 220;
smoothMPC_FLOAT* smoothMPC_dzcc20 = smoothMPC_dz_cc + 220;
smoothMPC_FLOAT* smoothMPC_rd20 = smoothMPC_rd + 220;
smoothMPC_FLOAT smoothMPC_Lbyrd20[11];
smoothMPC_FLOAT* smoothMPC_grad_cost20 = smoothMPC_grad_cost + 220;
smoothMPC_FLOAT* smoothMPC_grad_eq20 = smoothMPC_grad_eq + 220;
smoothMPC_FLOAT* smoothMPC_grad_ineq20 = smoothMPC_grad_ineq + 220;
smoothMPC_FLOAT smoothMPC_ctv20[11];
smoothMPC_FLOAT* smoothMPC_v20 = smoothMPC_v + 60;
smoothMPC_FLOAT smoothMPC_re20[3];
smoothMPC_FLOAT smoothMPC_beta20[3];
smoothMPC_FLOAT smoothMPC_betacc20[3];
smoothMPC_FLOAT* smoothMPC_dvaff20 = smoothMPC_dv_aff + 60;
smoothMPC_FLOAT* smoothMPC_dvcc20 = smoothMPC_dv_cc + 60;
smoothMPC_FLOAT smoothMPC_V20[33];
smoothMPC_FLOAT smoothMPC_Yd20[6];
smoothMPC_FLOAT smoothMPC_Ld20[6];
smoothMPC_FLOAT smoothMPC_yy20[3];
smoothMPC_FLOAT smoothMPC_bmy20[3];
int smoothMPC_lbIdx20[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb20 = smoothMPC_l + 320;
smoothMPC_FLOAT* smoothMPC_slb20 = smoothMPC_s + 320;
smoothMPC_FLOAT* smoothMPC_llbbyslb20 = smoothMPC_lbys + 320;
smoothMPC_FLOAT smoothMPC_rilb20[11];
smoothMPC_FLOAT* smoothMPC_dllbaff20 = smoothMPC_dl_aff + 320;
smoothMPC_FLOAT* smoothMPC_dslbaff20 = smoothMPC_ds_aff + 320;
smoothMPC_FLOAT* smoothMPC_dllbcc20 = smoothMPC_dl_cc + 320;
smoothMPC_FLOAT* smoothMPC_dslbcc20 = smoothMPC_ds_cc + 320;
smoothMPC_FLOAT* smoothMPC_ccrhsl20 = smoothMPC_ccrhs + 320;
int smoothMPC_ubIdx20[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub20 = smoothMPC_l + 331;
smoothMPC_FLOAT* smoothMPC_sub20 = smoothMPC_s + 331;
smoothMPC_FLOAT* smoothMPC_lubbysub20 = smoothMPC_lbys + 331;
smoothMPC_FLOAT smoothMPC_riub20[5];
smoothMPC_FLOAT* smoothMPC_dlubaff20 = smoothMPC_dl_aff + 331;
smoothMPC_FLOAT* smoothMPC_dsubaff20 = smoothMPC_ds_aff + 331;
smoothMPC_FLOAT* smoothMPC_dlubcc20 = smoothMPC_dl_cc + 331;
smoothMPC_FLOAT* smoothMPC_dsubcc20 = smoothMPC_ds_cc + 331;
smoothMPC_FLOAT* smoothMPC_ccrhsub20 = smoothMPC_ccrhs + 331;
smoothMPC_FLOAT smoothMPC_Phi20[11];
smoothMPC_FLOAT smoothMPC_W20[11];
smoothMPC_FLOAT smoothMPC_Ysd20[9];
smoothMPC_FLOAT smoothMPC_Lsd20[9];
smoothMPC_FLOAT* smoothMPC_z21 = smoothMPC_z + 231;
smoothMPC_FLOAT* smoothMPC_dzaff21 = smoothMPC_dz_aff + 231;
smoothMPC_FLOAT* smoothMPC_dzcc21 = smoothMPC_dz_cc + 231;
smoothMPC_FLOAT* smoothMPC_rd21 = smoothMPC_rd + 231;
smoothMPC_FLOAT smoothMPC_Lbyrd21[11];
smoothMPC_FLOAT* smoothMPC_grad_cost21 = smoothMPC_grad_cost + 231;
smoothMPC_FLOAT* smoothMPC_grad_eq21 = smoothMPC_grad_eq + 231;
smoothMPC_FLOAT* smoothMPC_grad_ineq21 = smoothMPC_grad_ineq + 231;
smoothMPC_FLOAT smoothMPC_ctv21[11];
smoothMPC_FLOAT* smoothMPC_v21 = smoothMPC_v + 63;
smoothMPC_FLOAT smoothMPC_re21[3];
smoothMPC_FLOAT smoothMPC_beta21[3];
smoothMPC_FLOAT smoothMPC_betacc21[3];
smoothMPC_FLOAT* smoothMPC_dvaff21 = smoothMPC_dv_aff + 63;
smoothMPC_FLOAT* smoothMPC_dvcc21 = smoothMPC_dv_cc + 63;
smoothMPC_FLOAT smoothMPC_V21[33];
smoothMPC_FLOAT smoothMPC_Yd21[6];
smoothMPC_FLOAT smoothMPC_Ld21[6];
smoothMPC_FLOAT smoothMPC_yy21[3];
smoothMPC_FLOAT smoothMPC_bmy21[3];
int smoothMPC_lbIdx21[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb21 = smoothMPC_l + 336;
smoothMPC_FLOAT* smoothMPC_slb21 = smoothMPC_s + 336;
smoothMPC_FLOAT* smoothMPC_llbbyslb21 = smoothMPC_lbys + 336;
smoothMPC_FLOAT smoothMPC_rilb21[11];
smoothMPC_FLOAT* smoothMPC_dllbaff21 = smoothMPC_dl_aff + 336;
smoothMPC_FLOAT* smoothMPC_dslbaff21 = smoothMPC_ds_aff + 336;
smoothMPC_FLOAT* smoothMPC_dllbcc21 = smoothMPC_dl_cc + 336;
smoothMPC_FLOAT* smoothMPC_dslbcc21 = smoothMPC_ds_cc + 336;
smoothMPC_FLOAT* smoothMPC_ccrhsl21 = smoothMPC_ccrhs + 336;
int smoothMPC_ubIdx21[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub21 = smoothMPC_l + 347;
smoothMPC_FLOAT* smoothMPC_sub21 = smoothMPC_s + 347;
smoothMPC_FLOAT* smoothMPC_lubbysub21 = smoothMPC_lbys + 347;
smoothMPC_FLOAT smoothMPC_riub21[5];
smoothMPC_FLOAT* smoothMPC_dlubaff21 = smoothMPC_dl_aff + 347;
smoothMPC_FLOAT* smoothMPC_dsubaff21 = smoothMPC_ds_aff + 347;
smoothMPC_FLOAT* smoothMPC_dlubcc21 = smoothMPC_dl_cc + 347;
smoothMPC_FLOAT* smoothMPC_dsubcc21 = smoothMPC_ds_cc + 347;
smoothMPC_FLOAT* smoothMPC_ccrhsub21 = smoothMPC_ccrhs + 347;
smoothMPC_FLOAT smoothMPC_Phi21[11];
smoothMPC_FLOAT smoothMPC_W21[11];
smoothMPC_FLOAT smoothMPC_Ysd21[9];
smoothMPC_FLOAT smoothMPC_Lsd21[9];
smoothMPC_FLOAT* smoothMPC_z22 = smoothMPC_z + 242;
smoothMPC_FLOAT* smoothMPC_dzaff22 = smoothMPC_dz_aff + 242;
smoothMPC_FLOAT* smoothMPC_dzcc22 = smoothMPC_dz_cc + 242;
smoothMPC_FLOAT* smoothMPC_rd22 = smoothMPC_rd + 242;
smoothMPC_FLOAT smoothMPC_Lbyrd22[11];
smoothMPC_FLOAT* smoothMPC_grad_cost22 = smoothMPC_grad_cost + 242;
smoothMPC_FLOAT* smoothMPC_grad_eq22 = smoothMPC_grad_eq + 242;
smoothMPC_FLOAT* smoothMPC_grad_ineq22 = smoothMPC_grad_ineq + 242;
smoothMPC_FLOAT smoothMPC_ctv22[11];
smoothMPC_FLOAT* smoothMPC_v22 = smoothMPC_v + 66;
smoothMPC_FLOAT smoothMPC_re22[3];
smoothMPC_FLOAT smoothMPC_beta22[3];
smoothMPC_FLOAT smoothMPC_betacc22[3];
smoothMPC_FLOAT* smoothMPC_dvaff22 = smoothMPC_dv_aff + 66;
smoothMPC_FLOAT* smoothMPC_dvcc22 = smoothMPC_dv_cc + 66;
smoothMPC_FLOAT smoothMPC_V22[33];
smoothMPC_FLOAT smoothMPC_Yd22[6];
smoothMPC_FLOAT smoothMPC_Ld22[6];
smoothMPC_FLOAT smoothMPC_yy22[3];
smoothMPC_FLOAT smoothMPC_bmy22[3];
int smoothMPC_lbIdx22[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb22 = smoothMPC_l + 352;
smoothMPC_FLOAT* smoothMPC_slb22 = smoothMPC_s + 352;
smoothMPC_FLOAT* smoothMPC_llbbyslb22 = smoothMPC_lbys + 352;
smoothMPC_FLOAT smoothMPC_rilb22[11];
smoothMPC_FLOAT* smoothMPC_dllbaff22 = smoothMPC_dl_aff + 352;
smoothMPC_FLOAT* smoothMPC_dslbaff22 = smoothMPC_ds_aff + 352;
smoothMPC_FLOAT* smoothMPC_dllbcc22 = smoothMPC_dl_cc + 352;
smoothMPC_FLOAT* smoothMPC_dslbcc22 = smoothMPC_ds_cc + 352;
smoothMPC_FLOAT* smoothMPC_ccrhsl22 = smoothMPC_ccrhs + 352;
int smoothMPC_ubIdx22[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub22 = smoothMPC_l + 363;
smoothMPC_FLOAT* smoothMPC_sub22 = smoothMPC_s + 363;
smoothMPC_FLOAT* smoothMPC_lubbysub22 = smoothMPC_lbys + 363;
smoothMPC_FLOAT smoothMPC_riub22[5];
smoothMPC_FLOAT* smoothMPC_dlubaff22 = smoothMPC_dl_aff + 363;
smoothMPC_FLOAT* smoothMPC_dsubaff22 = smoothMPC_ds_aff + 363;
smoothMPC_FLOAT* smoothMPC_dlubcc22 = smoothMPC_dl_cc + 363;
smoothMPC_FLOAT* smoothMPC_dsubcc22 = smoothMPC_ds_cc + 363;
smoothMPC_FLOAT* smoothMPC_ccrhsub22 = smoothMPC_ccrhs + 363;
smoothMPC_FLOAT smoothMPC_Phi22[11];
smoothMPC_FLOAT smoothMPC_W22[11];
smoothMPC_FLOAT smoothMPC_Ysd22[9];
smoothMPC_FLOAT smoothMPC_Lsd22[9];
smoothMPC_FLOAT* smoothMPC_z23 = smoothMPC_z + 253;
smoothMPC_FLOAT* smoothMPC_dzaff23 = smoothMPC_dz_aff + 253;
smoothMPC_FLOAT* smoothMPC_dzcc23 = smoothMPC_dz_cc + 253;
smoothMPC_FLOAT* smoothMPC_rd23 = smoothMPC_rd + 253;
smoothMPC_FLOAT smoothMPC_Lbyrd23[11];
smoothMPC_FLOAT* smoothMPC_grad_cost23 = smoothMPC_grad_cost + 253;
smoothMPC_FLOAT* smoothMPC_grad_eq23 = smoothMPC_grad_eq + 253;
smoothMPC_FLOAT* smoothMPC_grad_ineq23 = smoothMPC_grad_ineq + 253;
smoothMPC_FLOAT smoothMPC_ctv23[11];
smoothMPC_FLOAT* smoothMPC_v23 = smoothMPC_v + 69;
smoothMPC_FLOAT smoothMPC_re23[3];
smoothMPC_FLOAT smoothMPC_beta23[3];
smoothMPC_FLOAT smoothMPC_betacc23[3];
smoothMPC_FLOAT* smoothMPC_dvaff23 = smoothMPC_dv_aff + 69;
smoothMPC_FLOAT* smoothMPC_dvcc23 = smoothMPC_dv_cc + 69;
smoothMPC_FLOAT smoothMPC_V23[33];
smoothMPC_FLOAT smoothMPC_Yd23[6];
smoothMPC_FLOAT smoothMPC_Ld23[6];
smoothMPC_FLOAT smoothMPC_yy23[3];
smoothMPC_FLOAT smoothMPC_bmy23[3];
int smoothMPC_lbIdx23[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb23 = smoothMPC_l + 368;
smoothMPC_FLOAT* smoothMPC_slb23 = smoothMPC_s + 368;
smoothMPC_FLOAT* smoothMPC_llbbyslb23 = smoothMPC_lbys + 368;
smoothMPC_FLOAT smoothMPC_rilb23[11];
smoothMPC_FLOAT* smoothMPC_dllbaff23 = smoothMPC_dl_aff + 368;
smoothMPC_FLOAT* smoothMPC_dslbaff23 = smoothMPC_ds_aff + 368;
smoothMPC_FLOAT* smoothMPC_dllbcc23 = smoothMPC_dl_cc + 368;
smoothMPC_FLOAT* smoothMPC_dslbcc23 = smoothMPC_ds_cc + 368;
smoothMPC_FLOAT* smoothMPC_ccrhsl23 = smoothMPC_ccrhs + 368;
int smoothMPC_ubIdx23[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub23 = smoothMPC_l + 379;
smoothMPC_FLOAT* smoothMPC_sub23 = smoothMPC_s + 379;
smoothMPC_FLOAT* smoothMPC_lubbysub23 = smoothMPC_lbys + 379;
smoothMPC_FLOAT smoothMPC_riub23[5];
smoothMPC_FLOAT* smoothMPC_dlubaff23 = smoothMPC_dl_aff + 379;
smoothMPC_FLOAT* smoothMPC_dsubaff23 = smoothMPC_ds_aff + 379;
smoothMPC_FLOAT* smoothMPC_dlubcc23 = smoothMPC_dl_cc + 379;
smoothMPC_FLOAT* smoothMPC_dsubcc23 = smoothMPC_ds_cc + 379;
smoothMPC_FLOAT* smoothMPC_ccrhsub23 = smoothMPC_ccrhs + 379;
smoothMPC_FLOAT smoothMPC_Phi23[11];
smoothMPC_FLOAT smoothMPC_W23[11];
smoothMPC_FLOAT smoothMPC_Ysd23[9];
smoothMPC_FLOAT smoothMPC_Lsd23[9];
smoothMPC_FLOAT* smoothMPC_z24 = smoothMPC_z + 264;
smoothMPC_FLOAT* smoothMPC_dzaff24 = smoothMPC_dz_aff + 264;
smoothMPC_FLOAT* smoothMPC_dzcc24 = smoothMPC_dz_cc + 264;
smoothMPC_FLOAT* smoothMPC_rd24 = smoothMPC_rd + 264;
smoothMPC_FLOAT smoothMPC_Lbyrd24[11];
smoothMPC_FLOAT* smoothMPC_grad_cost24 = smoothMPC_grad_cost + 264;
smoothMPC_FLOAT* smoothMPC_grad_eq24 = smoothMPC_grad_eq + 264;
smoothMPC_FLOAT* smoothMPC_grad_ineq24 = smoothMPC_grad_ineq + 264;
smoothMPC_FLOAT smoothMPC_ctv24[11];
smoothMPC_FLOAT* smoothMPC_v24 = smoothMPC_v + 72;
smoothMPC_FLOAT smoothMPC_re24[3];
smoothMPC_FLOAT smoothMPC_beta24[3];
smoothMPC_FLOAT smoothMPC_betacc24[3];
smoothMPC_FLOAT* smoothMPC_dvaff24 = smoothMPC_dv_aff + 72;
smoothMPC_FLOAT* smoothMPC_dvcc24 = smoothMPC_dv_cc + 72;
smoothMPC_FLOAT smoothMPC_V24[33];
smoothMPC_FLOAT smoothMPC_Yd24[6];
smoothMPC_FLOAT smoothMPC_Ld24[6];
smoothMPC_FLOAT smoothMPC_yy24[3];
smoothMPC_FLOAT smoothMPC_bmy24[3];
int smoothMPC_lbIdx24[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb24 = smoothMPC_l + 384;
smoothMPC_FLOAT* smoothMPC_slb24 = smoothMPC_s + 384;
smoothMPC_FLOAT* smoothMPC_llbbyslb24 = smoothMPC_lbys + 384;
smoothMPC_FLOAT smoothMPC_rilb24[11];
smoothMPC_FLOAT* smoothMPC_dllbaff24 = smoothMPC_dl_aff + 384;
smoothMPC_FLOAT* smoothMPC_dslbaff24 = smoothMPC_ds_aff + 384;
smoothMPC_FLOAT* smoothMPC_dllbcc24 = smoothMPC_dl_cc + 384;
smoothMPC_FLOAT* smoothMPC_dslbcc24 = smoothMPC_ds_cc + 384;
smoothMPC_FLOAT* smoothMPC_ccrhsl24 = smoothMPC_ccrhs + 384;
int smoothMPC_ubIdx24[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub24 = smoothMPC_l + 395;
smoothMPC_FLOAT* smoothMPC_sub24 = smoothMPC_s + 395;
smoothMPC_FLOAT* smoothMPC_lubbysub24 = smoothMPC_lbys + 395;
smoothMPC_FLOAT smoothMPC_riub24[5];
smoothMPC_FLOAT* smoothMPC_dlubaff24 = smoothMPC_dl_aff + 395;
smoothMPC_FLOAT* smoothMPC_dsubaff24 = smoothMPC_ds_aff + 395;
smoothMPC_FLOAT* smoothMPC_dlubcc24 = smoothMPC_dl_cc + 395;
smoothMPC_FLOAT* smoothMPC_dsubcc24 = smoothMPC_ds_cc + 395;
smoothMPC_FLOAT* smoothMPC_ccrhsub24 = smoothMPC_ccrhs + 395;
smoothMPC_FLOAT smoothMPC_Phi24[11];
smoothMPC_FLOAT smoothMPC_W24[11];
smoothMPC_FLOAT smoothMPC_Ysd24[9];
smoothMPC_FLOAT smoothMPC_Lsd24[9];
smoothMPC_FLOAT* smoothMPC_z25 = smoothMPC_z + 275;
smoothMPC_FLOAT* smoothMPC_dzaff25 = smoothMPC_dz_aff + 275;
smoothMPC_FLOAT* smoothMPC_dzcc25 = smoothMPC_dz_cc + 275;
smoothMPC_FLOAT* smoothMPC_rd25 = smoothMPC_rd + 275;
smoothMPC_FLOAT smoothMPC_Lbyrd25[11];
smoothMPC_FLOAT* smoothMPC_grad_cost25 = smoothMPC_grad_cost + 275;
smoothMPC_FLOAT* smoothMPC_grad_eq25 = smoothMPC_grad_eq + 275;
smoothMPC_FLOAT* smoothMPC_grad_ineq25 = smoothMPC_grad_ineq + 275;
smoothMPC_FLOAT smoothMPC_ctv25[11];
smoothMPC_FLOAT* smoothMPC_v25 = smoothMPC_v + 75;
smoothMPC_FLOAT smoothMPC_re25[3];
smoothMPC_FLOAT smoothMPC_beta25[3];
smoothMPC_FLOAT smoothMPC_betacc25[3];
smoothMPC_FLOAT* smoothMPC_dvaff25 = smoothMPC_dv_aff + 75;
smoothMPC_FLOAT* smoothMPC_dvcc25 = smoothMPC_dv_cc + 75;
smoothMPC_FLOAT smoothMPC_V25[33];
smoothMPC_FLOAT smoothMPC_Yd25[6];
smoothMPC_FLOAT smoothMPC_Ld25[6];
smoothMPC_FLOAT smoothMPC_yy25[3];
smoothMPC_FLOAT smoothMPC_bmy25[3];
int smoothMPC_lbIdx25[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb25 = smoothMPC_l + 400;
smoothMPC_FLOAT* smoothMPC_slb25 = smoothMPC_s + 400;
smoothMPC_FLOAT* smoothMPC_llbbyslb25 = smoothMPC_lbys + 400;
smoothMPC_FLOAT smoothMPC_rilb25[11];
smoothMPC_FLOAT* smoothMPC_dllbaff25 = smoothMPC_dl_aff + 400;
smoothMPC_FLOAT* smoothMPC_dslbaff25 = smoothMPC_ds_aff + 400;
smoothMPC_FLOAT* smoothMPC_dllbcc25 = smoothMPC_dl_cc + 400;
smoothMPC_FLOAT* smoothMPC_dslbcc25 = smoothMPC_ds_cc + 400;
smoothMPC_FLOAT* smoothMPC_ccrhsl25 = smoothMPC_ccrhs + 400;
int smoothMPC_ubIdx25[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub25 = smoothMPC_l + 411;
smoothMPC_FLOAT* smoothMPC_sub25 = smoothMPC_s + 411;
smoothMPC_FLOAT* smoothMPC_lubbysub25 = smoothMPC_lbys + 411;
smoothMPC_FLOAT smoothMPC_riub25[5];
smoothMPC_FLOAT* smoothMPC_dlubaff25 = smoothMPC_dl_aff + 411;
smoothMPC_FLOAT* smoothMPC_dsubaff25 = smoothMPC_ds_aff + 411;
smoothMPC_FLOAT* smoothMPC_dlubcc25 = smoothMPC_dl_cc + 411;
smoothMPC_FLOAT* smoothMPC_dsubcc25 = smoothMPC_ds_cc + 411;
smoothMPC_FLOAT* smoothMPC_ccrhsub25 = smoothMPC_ccrhs + 411;
smoothMPC_FLOAT smoothMPC_Phi25[11];
smoothMPC_FLOAT smoothMPC_W25[11];
smoothMPC_FLOAT smoothMPC_Ysd25[9];
smoothMPC_FLOAT smoothMPC_Lsd25[9];
smoothMPC_FLOAT* smoothMPC_z26 = smoothMPC_z + 286;
smoothMPC_FLOAT* smoothMPC_dzaff26 = smoothMPC_dz_aff + 286;
smoothMPC_FLOAT* smoothMPC_dzcc26 = smoothMPC_dz_cc + 286;
smoothMPC_FLOAT* smoothMPC_rd26 = smoothMPC_rd + 286;
smoothMPC_FLOAT smoothMPC_Lbyrd26[11];
smoothMPC_FLOAT* smoothMPC_grad_cost26 = smoothMPC_grad_cost + 286;
smoothMPC_FLOAT* smoothMPC_grad_eq26 = smoothMPC_grad_eq + 286;
smoothMPC_FLOAT* smoothMPC_grad_ineq26 = smoothMPC_grad_ineq + 286;
smoothMPC_FLOAT smoothMPC_ctv26[11];
smoothMPC_FLOAT* smoothMPC_v26 = smoothMPC_v + 78;
smoothMPC_FLOAT smoothMPC_re26[3];
smoothMPC_FLOAT smoothMPC_beta26[3];
smoothMPC_FLOAT smoothMPC_betacc26[3];
smoothMPC_FLOAT* smoothMPC_dvaff26 = smoothMPC_dv_aff + 78;
smoothMPC_FLOAT* smoothMPC_dvcc26 = smoothMPC_dv_cc + 78;
smoothMPC_FLOAT smoothMPC_V26[33];
smoothMPC_FLOAT smoothMPC_Yd26[6];
smoothMPC_FLOAT smoothMPC_Ld26[6];
smoothMPC_FLOAT smoothMPC_yy26[3];
smoothMPC_FLOAT smoothMPC_bmy26[3];
int smoothMPC_lbIdx26[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb26 = smoothMPC_l + 416;
smoothMPC_FLOAT* smoothMPC_slb26 = smoothMPC_s + 416;
smoothMPC_FLOAT* smoothMPC_llbbyslb26 = smoothMPC_lbys + 416;
smoothMPC_FLOAT smoothMPC_rilb26[11];
smoothMPC_FLOAT* smoothMPC_dllbaff26 = smoothMPC_dl_aff + 416;
smoothMPC_FLOAT* smoothMPC_dslbaff26 = smoothMPC_ds_aff + 416;
smoothMPC_FLOAT* smoothMPC_dllbcc26 = smoothMPC_dl_cc + 416;
smoothMPC_FLOAT* smoothMPC_dslbcc26 = smoothMPC_ds_cc + 416;
smoothMPC_FLOAT* smoothMPC_ccrhsl26 = smoothMPC_ccrhs + 416;
int smoothMPC_ubIdx26[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub26 = smoothMPC_l + 427;
smoothMPC_FLOAT* smoothMPC_sub26 = smoothMPC_s + 427;
smoothMPC_FLOAT* smoothMPC_lubbysub26 = smoothMPC_lbys + 427;
smoothMPC_FLOAT smoothMPC_riub26[5];
smoothMPC_FLOAT* smoothMPC_dlubaff26 = smoothMPC_dl_aff + 427;
smoothMPC_FLOAT* smoothMPC_dsubaff26 = smoothMPC_ds_aff + 427;
smoothMPC_FLOAT* smoothMPC_dlubcc26 = smoothMPC_dl_cc + 427;
smoothMPC_FLOAT* smoothMPC_dsubcc26 = smoothMPC_ds_cc + 427;
smoothMPC_FLOAT* smoothMPC_ccrhsub26 = smoothMPC_ccrhs + 427;
smoothMPC_FLOAT smoothMPC_Phi26[11];
smoothMPC_FLOAT smoothMPC_W26[11];
smoothMPC_FLOAT smoothMPC_Ysd26[9];
smoothMPC_FLOAT smoothMPC_Lsd26[9];
smoothMPC_FLOAT* smoothMPC_z27 = smoothMPC_z + 297;
smoothMPC_FLOAT* smoothMPC_dzaff27 = smoothMPC_dz_aff + 297;
smoothMPC_FLOAT* smoothMPC_dzcc27 = smoothMPC_dz_cc + 297;
smoothMPC_FLOAT* smoothMPC_rd27 = smoothMPC_rd + 297;
smoothMPC_FLOAT smoothMPC_Lbyrd27[11];
smoothMPC_FLOAT* smoothMPC_grad_cost27 = smoothMPC_grad_cost + 297;
smoothMPC_FLOAT* smoothMPC_grad_eq27 = smoothMPC_grad_eq + 297;
smoothMPC_FLOAT* smoothMPC_grad_ineq27 = smoothMPC_grad_ineq + 297;
smoothMPC_FLOAT smoothMPC_ctv27[11];
smoothMPC_FLOAT* smoothMPC_v27 = smoothMPC_v + 81;
smoothMPC_FLOAT smoothMPC_re27[3];
smoothMPC_FLOAT smoothMPC_beta27[3];
smoothMPC_FLOAT smoothMPC_betacc27[3];
smoothMPC_FLOAT* smoothMPC_dvaff27 = smoothMPC_dv_aff + 81;
smoothMPC_FLOAT* smoothMPC_dvcc27 = smoothMPC_dv_cc + 81;
smoothMPC_FLOAT smoothMPC_V27[33];
smoothMPC_FLOAT smoothMPC_Yd27[6];
smoothMPC_FLOAT smoothMPC_Ld27[6];
smoothMPC_FLOAT smoothMPC_yy27[3];
smoothMPC_FLOAT smoothMPC_bmy27[3];
int smoothMPC_lbIdx27[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb27 = smoothMPC_l + 432;
smoothMPC_FLOAT* smoothMPC_slb27 = smoothMPC_s + 432;
smoothMPC_FLOAT* smoothMPC_llbbyslb27 = smoothMPC_lbys + 432;
smoothMPC_FLOAT smoothMPC_rilb27[11];
smoothMPC_FLOAT* smoothMPC_dllbaff27 = smoothMPC_dl_aff + 432;
smoothMPC_FLOAT* smoothMPC_dslbaff27 = smoothMPC_ds_aff + 432;
smoothMPC_FLOAT* smoothMPC_dllbcc27 = smoothMPC_dl_cc + 432;
smoothMPC_FLOAT* smoothMPC_dslbcc27 = smoothMPC_ds_cc + 432;
smoothMPC_FLOAT* smoothMPC_ccrhsl27 = smoothMPC_ccrhs + 432;
int smoothMPC_ubIdx27[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub27 = smoothMPC_l + 443;
smoothMPC_FLOAT* smoothMPC_sub27 = smoothMPC_s + 443;
smoothMPC_FLOAT* smoothMPC_lubbysub27 = smoothMPC_lbys + 443;
smoothMPC_FLOAT smoothMPC_riub27[5];
smoothMPC_FLOAT* smoothMPC_dlubaff27 = smoothMPC_dl_aff + 443;
smoothMPC_FLOAT* smoothMPC_dsubaff27 = smoothMPC_ds_aff + 443;
smoothMPC_FLOAT* smoothMPC_dlubcc27 = smoothMPC_dl_cc + 443;
smoothMPC_FLOAT* smoothMPC_dsubcc27 = smoothMPC_ds_cc + 443;
smoothMPC_FLOAT* smoothMPC_ccrhsub27 = smoothMPC_ccrhs + 443;
smoothMPC_FLOAT smoothMPC_Phi27[11];
smoothMPC_FLOAT smoothMPC_W27[11];
smoothMPC_FLOAT smoothMPC_Ysd27[9];
smoothMPC_FLOAT smoothMPC_Lsd27[9];
smoothMPC_FLOAT* smoothMPC_z28 = smoothMPC_z + 308;
smoothMPC_FLOAT* smoothMPC_dzaff28 = smoothMPC_dz_aff + 308;
smoothMPC_FLOAT* smoothMPC_dzcc28 = smoothMPC_dz_cc + 308;
smoothMPC_FLOAT* smoothMPC_rd28 = smoothMPC_rd + 308;
smoothMPC_FLOAT smoothMPC_Lbyrd28[11];
smoothMPC_FLOAT* smoothMPC_grad_cost28 = smoothMPC_grad_cost + 308;
smoothMPC_FLOAT* smoothMPC_grad_eq28 = smoothMPC_grad_eq + 308;
smoothMPC_FLOAT* smoothMPC_grad_ineq28 = smoothMPC_grad_ineq + 308;
smoothMPC_FLOAT smoothMPC_ctv28[11];
smoothMPC_FLOAT* smoothMPC_v28 = smoothMPC_v + 84;
smoothMPC_FLOAT smoothMPC_re28[3];
smoothMPC_FLOAT smoothMPC_beta28[3];
smoothMPC_FLOAT smoothMPC_betacc28[3];
smoothMPC_FLOAT* smoothMPC_dvaff28 = smoothMPC_dv_aff + 84;
smoothMPC_FLOAT* smoothMPC_dvcc28 = smoothMPC_dv_cc + 84;
smoothMPC_FLOAT smoothMPC_V28[33];
smoothMPC_FLOAT smoothMPC_Yd28[6];
smoothMPC_FLOAT smoothMPC_Ld28[6];
smoothMPC_FLOAT smoothMPC_yy28[3];
smoothMPC_FLOAT smoothMPC_bmy28[3];
int smoothMPC_lbIdx28[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
smoothMPC_FLOAT* smoothMPC_llb28 = smoothMPC_l + 448;
smoothMPC_FLOAT* smoothMPC_slb28 = smoothMPC_s + 448;
smoothMPC_FLOAT* smoothMPC_llbbyslb28 = smoothMPC_lbys + 448;
smoothMPC_FLOAT smoothMPC_rilb28[11];
smoothMPC_FLOAT* smoothMPC_dllbaff28 = smoothMPC_dl_aff + 448;
smoothMPC_FLOAT* smoothMPC_dslbaff28 = smoothMPC_ds_aff + 448;
smoothMPC_FLOAT* smoothMPC_dllbcc28 = smoothMPC_dl_cc + 448;
smoothMPC_FLOAT* smoothMPC_dslbcc28 = smoothMPC_ds_cc + 448;
smoothMPC_FLOAT* smoothMPC_ccrhsl28 = smoothMPC_ccrhs + 448;
int smoothMPC_ubIdx28[5] = {0, 1, 2, 3, 4};
smoothMPC_FLOAT* smoothMPC_lub28 = smoothMPC_l + 459;
smoothMPC_FLOAT* smoothMPC_sub28 = smoothMPC_s + 459;
smoothMPC_FLOAT* smoothMPC_lubbysub28 = smoothMPC_lbys + 459;
smoothMPC_FLOAT smoothMPC_riub28[5];
smoothMPC_FLOAT* smoothMPC_dlubaff28 = smoothMPC_dl_aff + 459;
smoothMPC_FLOAT* smoothMPC_dsubaff28 = smoothMPC_ds_aff + 459;
smoothMPC_FLOAT* smoothMPC_dlubcc28 = smoothMPC_dl_cc + 459;
smoothMPC_FLOAT* smoothMPC_dsubcc28 = smoothMPC_ds_cc + 459;
smoothMPC_FLOAT* smoothMPC_ccrhsub28 = smoothMPC_ccrhs + 459;
smoothMPC_FLOAT smoothMPC_Phi28[11];
smoothMPC_FLOAT smoothMPC_W28[11];
smoothMPC_FLOAT smoothMPC_Ysd28[9];
smoothMPC_FLOAT smoothMPC_Lsd28[9];
smoothMPC_FLOAT* smoothMPC_z29 = smoothMPC_z + 319;
smoothMPC_FLOAT* smoothMPC_dzaff29 = smoothMPC_dz_aff + 319;
smoothMPC_FLOAT* smoothMPC_dzcc29 = smoothMPC_dz_cc + 319;
smoothMPC_FLOAT* smoothMPC_rd29 = smoothMPC_rd + 319;
smoothMPC_FLOAT smoothMPC_Lbyrd29[3];
smoothMPC_FLOAT* smoothMPC_grad_cost29 = smoothMPC_grad_cost + 319;
smoothMPC_FLOAT* smoothMPC_grad_eq29 = smoothMPC_grad_eq + 319;
smoothMPC_FLOAT* smoothMPC_grad_ineq29 = smoothMPC_grad_ineq + 319;
smoothMPC_FLOAT smoothMPC_ctv29[3];
smoothMPC_FLOAT* smoothMPC_v29 = smoothMPC_v + 87;
smoothMPC_FLOAT smoothMPC_re29[3];
smoothMPC_FLOAT smoothMPC_beta29[3];
smoothMPC_FLOAT smoothMPC_betacc29[3];
smoothMPC_FLOAT* smoothMPC_dvaff29 = smoothMPC_dv_aff + 87;
smoothMPC_FLOAT* smoothMPC_dvcc29 = smoothMPC_dv_cc + 87;
smoothMPC_FLOAT smoothMPC_V29[9];
smoothMPC_FLOAT smoothMPC_Yd29[6];
smoothMPC_FLOAT smoothMPC_Ld29[6];
smoothMPC_FLOAT smoothMPC_yy29[3];
smoothMPC_FLOAT smoothMPC_bmy29[3];
int smoothMPC_lbIdx29[3] = {0, 1, 2};
smoothMPC_FLOAT* smoothMPC_llb29 = smoothMPC_l + 464;
smoothMPC_FLOAT* smoothMPC_slb29 = smoothMPC_s + 464;
smoothMPC_FLOAT* smoothMPC_llbbyslb29 = smoothMPC_lbys + 464;
smoothMPC_FLOAT smoothMPC_rilb29[3];
smoothMPC_FLOAT* smoothMPC_dllbaff29 = smoothMPC_dl_aff + 464;
smoothMPC_FLOAT* smoothMPC_dslbaff29 = smoothMPC_ds_aff + 464;
smoothMPC_FLOAT* smoothMPC_dllbcc29 = smoothMPC_dl_cc + 464;
smoothMPC_FLOAT* smoothMPC_dslbcc29 = smoothMPC_ds_cc + 464;
smoothMPC_FLOAT* smoothMPC_ccrhsl29 = smoothMPC_ccrhs + 464;
int smoothMPC_ubIdx29[3] = {0, 1, 2};
smoothMPC_FLOAT* smoothMPC_lub29 = smoothMPC_l + 467;
smoothMPC_FLOAT* smoothMPC_sub29 = smoothMPC_s + 467;
smoothMPC_FLOAT* smoothMPC_lubbysub29 = smoothMPC_lbys + 467;
smoothMPC_FLOAT smoothMPC_riub29[3];
smoothMPC_FLOAT* smoothMPC_dlubaff29 = smoothMPC_dl_aff + 467;
smoothMPC_FLOAT* smoothMPC_dsubaff29 = smoothMPC_ds_aff + 467;
smoothMPC_FLOAT* smoothMPC_dlubcc29 = smoothMPC_dl_cc + 467;
smoothMPC_FLOAT* smoothMPC_dsubcc29 = smoothMPC_ds_cc + 467;
smoothMPC_FLOAT* smoothMPC_ccrhsub29 = smoothMPC_ccrhs + 467;
smoothMPC_FLOAT smoothMPC_Phi29[3];
smoothMPC_FLOAT smoothMPC_D29[3] = {-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000};
smoothMPC_FLOAT smoothMPC_W29[3];
smoothMPC_FLOAT smoothMPC_Ysd29[9];
smoothMPC_FLOAT smoothMPC_Lsd29[9];
smoothMPC_FLOAT musigma;
smoothMPC_FLOAT sigma_3rdroot;
smoothMPC_FLOAT smoothMPC_Diag1_0[11];
smoothMPC_FLOAT smoothMPC_Diag2_0[11];
smoothMPC_FLOAT smoothMPC_L_0[55];




/* SOLVER CODE --------------------------------------------------------- */
int smoothMPC_solve(smoothMPC_params* params, smoothMPC_output* output, smoothMPC_info* info)
{	
int exitcode;

#if smoothMPC_SET_TIMING == 1
	smoothMPC_timer solvertimer;
	smoothMPC_tic(&solvertimer);
#endif
/* FUNCTION CALLS INTO LA LIBRARY -------------------------------------- */
info->it = 0;
smoothMPC_LA_INITIALIZEVECTOR_322(smoothMPC_z, 0);
smoothMPC_LA_INITIALIZEVECTOR_90(smoothMPC_v, 1);
smoothMPC_LA_INITIALIZEVECTOR_470(smoothMPC_l, 1);
smoothMPC_LA_INITIALIZEVECTOR_470(smoothMPC_s, 1);
info->mu = 0;
smoothMPC_LA_DOTACC_470(smoothMPC_l, smoothMPC_s, &info->mu);
info->mu /= 470;
PRINTTEXT("This is smoothMPC, a solver generated by FORCES (forces.ethz.ch).\n");
PRINTTEXT("(c) Alexander Domahidi, Automatic Control Laboratory, ETH Zurich, 2011-2014.\n");
PRINTTEXT("\n  #it  res_eq   res_ineq     pobj         dobj       dgap     rdgap     mu\n");
PRINTTEXT("  ---------------------------------------------------------------------------\n");
while( 1 ){
info->pobj = 0;
smoothMPC_LA_DIAG_QUADFCN_11(params->H1, params->f1, smoothMPC_z00, smoothMPC_grad_cost00, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H2, params->f2, smoothMPC_z01, smoothMPC_grad_cost01, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H3, params->f3, smoothMPC_z02, smoothMPC_grad_cost02, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H4, params->f4, smoothMPC_z03, smoothMPC_grad_cost03, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H5, params->f5, smoothMPC_z04, smoothMPC_grad_cost04, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H6, params->f6, smoothMPC_z05, smoothMPC_grad_cost05, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H7, params->f7, smoothMPC_z06, smoothMPC_grad_cost06, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H8, params->f8, smoothMPC_z07, smoothMPC_grad_cost07, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H9, params->f9, smoothMPC_z08, smoothMPC_grad_cost08, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H10, params->f10, smoothMPC_z09, smoothMPC_grad_cost09, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H11, params->f11, smoothMPC_z10, smoothMPC_grad_cost10, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H12, params->f12, smoothMPC_z11, smoothMPC_grad_cost11, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H13, params->f13, smoothMPC_z12, smoothMPC_grad_cost12, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H14, params->f14, smoothMPC_z13, smoothMPC_grad_cost13, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H15, params->f15, smoothMPC_z14, smoothMPC_grad_cost14, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H16, params->f16, smoothMPC_z15, smoothMPC_grad_cost15, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H17, params->f17, smoothMPC_z16, smoothMPC_grad_cost16, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H18, params->f18, smoothMPC_z17, smoothMPC_grad_cost17, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H19, params->f19, smoothMPC_z18, smoothMPC_grad_cost18, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H20, params->f20, smoothMPC_z19, smoothMPC_grad_cost19, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H21, params->f21, smoothMPC_z20, smoothMPC_grad_cost20, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H22, params->f22, smoothMPC_z21, smoothMPC_grad_cost21, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H23, params->f23, smoothMPC_z22, smoothMPC_grad_cost22, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H24, params->f24, smoothMPC_z23, smoothMPC_grad_cost23, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H25, params->f25, smoothMPC_z24, smoothMPC_grad_cost24, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H26, params->f26, smoothMPC_z25, smoothMPC_grad_cost25, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H27, params->f27, smoothMPC_z26, smoothMPC_grad_cost26, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H28, params->f28, smoothMPC_z27, smoothMPC_grad_cost27, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_11(params->H29, params->f29, smoothMPC_z28, smoothMPC_grad_cost28, &info->pobj);
smoothMPC_LA_DIAG_QUADFCN_3(params->H30, params->f30, smoothMPC_z29, smoothMPC_grad_cost29, &info->pobj);
info->res_eq = 0;
info->dgap = 0;
smoothMPC_LA_DIAGZERO_MVMSUB6_3(smoothMPC_D00, smoothMPC_z00, params->e1, smoothMPC_v00, smoothMPC_re00, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C1, smoothMPC_z00, smoothMPC_D01, smoothMPC_z01, params->e2, smoothMPC_v01, smoothMPC_re01, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C2, smoothMPC_z01, smoothMPC_D01, smoothMPC_z02, params->e3, smoothMPC_v02, smoothMPC_re02, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C3, smoothMPC_z02, smoothMPC_D01, smoothMPC_z03, params->e4, smoothMPC_v03, smoothMPC_re03, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C4, smoothMPC_z03, smoothMPC_D01, smoothMPC_z04, params->e5, smoothMPC_v04, smoothMPC_re04, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C5, smoothMPC_z04, smoothMPC_D01, smoothMPC_z05, params->e6, smoothMPC_v05, smoothMPC_re05, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C6, smoothMPC_z05, smoothMPC_D01, smoothMPC_z06, params->e7, smoothMPC_v06, smoothMPC_re06, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C7, smoothMPC_z06, smoothMPC_D01, smoothMPC_z07, params->e8, smoothMPC_v07, smoothMPC_re07, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C8, smoothMPC_z07, smoothMPC_D01, smoothMPC_z08, params->e9, smoothMPC_v08, smoothMPC_re08, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C9, smoothMPC_z08, smoothMPC_D01, smoothMPC_z09, params->e10, smoothMPC_v09, smoothMPC_re09, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C10, smoothMPC_z09, smoothMPC_D01, smoothMPC_z10, params->e11, smoothMPC_v10, smoothMPC_re10, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C11, smoothMPC_z10, smoothMPC_D01, smoothMPC_z11, params->e12, smoothMPC_v11, smoothMPC_re11, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C12, smoothMPC_z11, smoothMPC_D01, smoothMPC_z12, params->e13, smoothMPC_v12, smoothMPC_re12, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C13, smoothMPC_z12, smoothMPC_D01, smoothMPC_z13, params->e14, smoothMPC_v13, smoothMPC_re13, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C14, smoothMPC_z13, smoothMPC_D01, smoothMPC_z14, params->e15, smoothMPC_v14, smoothMPC_re14, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C15, smoothMPC_z14, smoothMPC_D01, smoothMPC_z15, params->e16, smoothMPC_v15, smoothMPC_re15, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C16, smoothMPC_z15, smoothMPC_D01, smoothMPC_z16, params->e17, smoothMPC_v16, smoothMPC_re16, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C17, smoothMPC_z16, smoothMPC_D01, smoothMPC_z17, params->e18, smoothMPC_v17, smoothMPC_re17, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C18, smoothMPC_z17, smoothMPC_D01, smoothMPC_z18, params->e19, smoothMPC_v18, smoothMPC_re18, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C19, smoothMPC_z18, smoothMPC_D01, smoothMPC_z19, params->e20, smoothMPC_v19, smoothMPC_re19, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C20, smoothMPC_z19, smoothMPC_D01, smoothMPC_z20, params->e21, smoothMPC_v20, smoothMPC_re20, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C21, smoothMPC_z20, smoothMPC_D01, smoothMPC_z21, params->e22, smoothMPC_v21, smoothMPC_re21, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C22, smoothMPC_z21, smoothMPC_D01, smoothMPC_z22, params->e23, smoothMPC_v22, smoothMPC_re22, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C23, smoothMPC_z22, smoothMPC_D01, smoothMPC_z23, params->e24, smoothMPC_v23, smoothMPC_re23, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C24, smoothMPC_z23, smoothMPC_D01, smoothMPC_z24, params->e25, smoothMPC_v24, smoothMPC_re24, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C25, smoothMPC_z24, smoothMPC_D01, smoothMPC_z25, params->e26, smoothMPC_v25, smoothMPC_re25, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C26, smoothMPC_z25, smoothMPC_D01, smoothMPC_z26, params->e27, smoothMPC_v26, smoothMPC_re26, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C27, smoothMPC_z26, smoothMPC_D01, smoothMPC_z27, params->e28, smoothMPC_v27, smoothMPC_re27, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C28, smoothMPC_z27, smoothMPC_D01, smoothMPC_z28, params->e29, smoothMPC_v28, smoothMPC_re28, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_3(params->C29, smoothMPC_z28, smoothMPC_D29, smoothMPC_z29, params->e30, smoothMPC_v29, smoothMPC_re29, &info->dgap, &info->res_eq);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C1, smoothMPC_v01, smoothMPC_D00, smoothMPC_v00, smoothMPC_grad_eq00);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C2, smoothMPC_v02, smoothMPC_D01, smoothMPC_v01, smoothMPC_grad_eq01);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C3, smoothMPC_v03, smoothMPC_D01, smoothMPC_v02, smoothMPC_grad_eq02);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C4, smoothMPC_v04, smoothMPC_D01, smoothMPC_v03, smoothMPC_grad_eq03);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C5, smoothMPC_v05, smoothMPC_D01, smoothMPC_v04, smoothMPC_grad_eq04);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C6, smoothMPC_v06, smoothMPC_D01, smoothMPC_v05, smoothMPC_grad_eq05);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C7, smoothMPC_v07, smoothMPC_D01, smoothMPC_v06, smoothMPC_grad_eq06);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C8, smoothMPC_v08, smoothMPC_D01, smoothMPC_v07, smoothMPC_grad_eq07);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C9, smoothMPC_v09, smoothMPC_D01, smoothMPC_v08, smoothMPC_grad_eq08);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C10, smoothMPC_v10, smoothMPC_D01, smoothMPC_v09, smoothMPC_grad_eq09);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C11, smoothMPC_v11, smoothMPC_D01, smoothMPC_v10, smoothMPC_grad_eq10);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C12, smoothMPC_v12, smoothMPC_D01, smoothMPC_v11, smoothMPC_grad_eq11);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C13, smoothMPC_v13, smoothMPC_D01, smoothMPC_v12, smoothMPC_grad_eq12);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C14, smoothMPC_v14, smoothMPC_D01, smoothMPC_v13, smoothMPC_grad_eq13);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C15, smoothMPC_v15, smoothMPC_D01, smoothMPC_v14, smoothMPC_grad_eq14);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C16, smoothMPC_v16, smoothMPC_D01, smoothMPC_v15, smoothMPC_grad_eq15);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C17, smoothMPC_v17, smoothMPC_D01, smoothMPC_v16, smoothMPC_grad_eq16);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C18, smoothMPC_v18, smoothMPC_D01, smoothMPC_v17, smoothMPC_grad_eq17);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C19, smoothMPC_v19, smoothMPC_D01, smoothMPC_v18, smoothMPC_grad_eq18);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C20, smoothMPC_v20, smoothMPC_D01, smoothMPC_v19, smoothMPC_grad_eq19);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C21, smoothMPC_v21, smoothMPC_D01, smoothMPC_v20, smoothMPC_grad_eq20);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C22, smoothMPC_v22, smoothMPC_D01, smoothMPC_v21, smoothMPC_grad_eq21);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C23, smoothMPC_v23, smoothMPC_D01, smoothMPC_v22, smoothMPC_grad_eq22);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C24, smoothMPC_v24, smoothMPC_D01, smoothMPC_v23, smoothMPC_grad_eq23);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C25, smoothMPC_v25, smoothMPC_D01, smoothMPC_v24, smoothMPC_grad_eq24);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C26, smoothMPC_v26, smoothMPC_D01, smoothMPC_v25, smoothMPC_grad_eq25);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C27, smoothMPC_v27, smoothMPC_D01, smoothMPC_v26, smoothMPC_grad_eq26);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C28, smoothMPC_v28, smoothMPC_D01, smoothMPC_v27, smoothMPC_grad_eq27);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C29, smoothMPC_v29, smoothMPC_D01, smoothMPC_v28, smoothMPC_grad_eq28);
smoothMPC_LA_DIAGZERO_MTVM_3_3(smoothMPC_D29, smoothMPC_v29, smoothMPC_grad_eq29);
info->res_ineq = 0;
smoothMPC_LA_VSUBADD3_11(params->lb1, smoothMPC_z00, smoothMPC_lbIdx00, smoothMPC_llb00, smoothMPC_slb00, smoothMPC_rilb00, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z00, smoothMPC_ubIdx00, params->ub1, smoothMPC_lub00, smoothMPC_sub00, smoothMPC_riub00, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb2, smoothMPC_z01, smoothMPC_lbIdx01, smoothMPC_llb01, smoothMPC_slb01, smoothMPC_rilb01, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z01, smoothMPC_ubIdx01, params->ub2, smoothMPC_lub01, smoothMPC_sub01, smoothMPC_riub01, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb3, smoothMPC_z02, smoothMPC_lbIdx02, smoothMPC_llb02, smoothMPC_slb02, smoothMPC_rilb02, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z02, smoothMPC_ubIdx02, params->ub3, smoothMPC_lub02, smoothMPC_sub02, smoothMPC_riub02, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb4, smoothMPC_z03, smoothMPC_lbIdx03, smoothMPC_llb03, smoothMPC_slb03, smoothMPC_rilb03, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z03, smoothMPC_ubIdx03, params->ub4, smoothMPC_lub03, smoothMPC_sub03, smoothMPC_riub03, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb5, smoothMPC_z04, smoothMPC_lbIdx04, smoothMPC_llb04, smoothMPC_slb04, smoothMPC_rilb04, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z04, smoothMPC_ubIdx04, params->ub5, smoothMPC_lub04, smoothMPC_sub04, smoothMPC_riub04, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb6, smoothMPC_z05, smoothMPC_lbIdx05, smoothMPC_llb05, smoothMPC_slb05, smoothMPC_rilb05, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z05, smoothMPC_ubIdx05, params->ub6, smoothMPC_lub05, smoothMPC_sub05, smoothMPC_riub05, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb7, smoothMPC_z06, smoothMPC_lbIdx06, smoothMPC_llb06, smoothMPC_slb06, smoothMPC_rilb06, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z06, smoothMPC_ubIdx06, params->ub7, smoothMPC_lub06, smoothMPC_sub06, smoothMPC_riub06, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb8, smoothMPC_z07, smoothMPC_lbIdx07, smoothMPC_llb07, smoothMPC_slb07, smoothMPC_rilb07, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z07, smoothMPC_ubIdx07, params->ub8, smoothMPC_lub07, smoothMPC_sub07, smoothMPC_riub07, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb9, smoothMPC_z08, smoothMPC_lbIdx08, smoothMPC_llb08, smoothMPC_slb08, smoothMPC_rilb08, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z08, smoothMPC_ubIdx08, params->ub9, smoothMPC_lub08, smoothMPC_sub08, smoothMPC_riub08, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb10, smoothMPC_z09, smoothMPC_lbIdx09, smoothMPC_llb09, smoothMPC_slb09, smoothMPC_rilb09, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z09, smoothMPC_ubIdx09, params->ub10, smoothMPC_lub09, smoothMPC_sub09, smoothMPC_riub09, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb11, smoothMPC_z10, smoothMPC_lbIdx10, smoothMPC_llb10, smoothMPC_slb10, smoothMPC_rilb10, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z10, smoothMPC_ubIdx10, params->ub11, smoothMPC_lub10, smoothMPC_sub10, smoothMPC_riub10, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb12, smoothMPC_z11, smoothMPC_lbIdx11, smoothMPC_llb11, smoothMPC_slb11, smoothMPC_rilb11, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z11, smoothMPC_ubIdx11, params->ub12, smoothMPC_lub11, smoothMPC_sub11, smoothMPC_riub11, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb13, smoothMPC_z12, smoothMPC_lbIdx12, smoothMPC_llb12, smoothMPC_slb12, smoothMPC_rilb12, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z12, smoothMPC_ubIdx12, params->ub13, smoothMPC_lub12, smoothMPC_sub12, smoothMPC_riub12, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb14, smoothMPC_z13, smoothMPC_lbIdx13, smoothMPC_llb13, smoothMPC_slb13, smoothMPC_rilb13, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z13, smoothMPC_ubIdx13, params->ub14, smoothMPC_lub13, smoothMPC_sub13, smoothMPC_riub13, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb15, smoothMPC_z14, smoothMPC_lbIdx14, smoothMPC_llb14, smoothMPC_slb14, smoothMPC_rilb14, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z14, smoothMPC_ubIdx14, params->ub15, smoothMPC_lub14, smoothMPC_sub14, smoothMPC_riub14, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb16, smoothMPC_z15, smoothMPC_lbIdx15, smoothMPC_llb15, smoothMPC_slb15, smoothMPC_rilb15, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z15, smoothMPC_ubIdx15, params->ub16, smoothMPC_lub15, smoothMPC_sub15, smoothMPC_riub15, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb17, smoothMPC_z16, smoothMPC_lbIdx16, smoothMPC_llb16, smoothMPC_slb16, smoothMPC_rilb16, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z16, smoothMPC_ubIdx16, params->ub17, smoothMPC_lub16, smoothMPC_sub16, smoothMPC_riub16, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb18, smoothMPC_z17, smoothMPC_lbIdx17, smoothMPC_llb17, smoothMPC_slb17, smoothMPC_rilb17, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z17, smoothMPC_ubIdx17, params->ub18, smoothMPC_lub17, smoothMPC_sub17, smoothMPC_riub17, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb19, smoothMPC_z18, smoothMPC_lbIdx18, smoothMPC_llb18, smoothMPC_slb18, smoothMPC_rilb18, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z18, smoothMPC_ubIdx18, params->ub19, smoothMPC_lub18, smoothMPC_sub18, smoothMPC_riub18, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb20, smoothMPC_z19, smoothMPC_lbIdx19, smoothMPC_llb19, smoothMPC_slb19, smoothMPC_rilb19, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z19, smoothMPC_ubIdx19, params->ub20, smoothMPC_lub19, smoothMPC_sub19, smoothMPC_riub19, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb21, smoothMPC_z20, smoothMPC_lbIdx20, smoothMPC_llb20, smoothMPC_slb20, smoothMPC_rilb20, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z20, smoothMPC_ubIdx20, params->ub21, smoothMPC_lub20, smoothMPC_sub20, smoothMPC_riub20, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb22, smoothMPC_z21, smoothMPC_lbIdx21, smoothMPC_llb21, smoothMPC_slb21, smoothMPC_rilb21, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z21, smoothMPC_ubIdx21, params->ub22, smoothMPC_lub21, smoothMPC_sub21, smoothMPC_riub21, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb23, smoothMPC_z22, smoothMPC_lbIdx22, smoothMPC_llb22, smoothMPC_slb22, smoothMPC_rilb22, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z22, smoothMPC_ubIdx22, params->ub23, smoothMPC_lub22, smoothMPC_sub22, smoothMPC_riub22, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb24, smoothMPC_z23, smoothMPC_lbIdx23, smoothMPC_llb23, smoothMPC_slb23, smoothMPC_rilb23, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z23, smoothMPC_ubIdx23, params->ub24, smoothMPC_lub23, smoothMPC_sub23, smoothMPC_riub23, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb25, smoothMPC_z24, smoothMPC_lbIdx24, smoothMPC_llb24, smoothMPC_slb24, smoothMPC_rilb24, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z24, smoothMPC_ubIdx24, params->ub25, smoothMPC_lub24, smoothMPC_sub24, smoothMPC_riub24, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb26, smoothMPC_z25, smoothMPC_lbIdx25, smoothMPC_llb25, smoothMPC_slb25, smoothMPC_rilb25, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z25, smoothMPC_ubIdx25, params->ub26, smoothMPC_lub25, smoothMPC_sub25, smoothMPC_riub25, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb27, smoothMPC_z26, smoothMPC_lbIdx26, smoothMPC_llb26, smoothMPC_slb26, smoothMPC_rilb26, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z26, smoothMPC_ubIdx26, params->ub27, smoothMPC_lub26, smoothMPC_sub26, smoothMPC_riub26, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb28, smoothMPC_z27, smoothMPC_lbIdx27, smoothMPC_llb27, smoothMPC_slb27, smoothMPC_rilb27, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z27, smoothMPC_ubIdx27, params->ub28, smoothMPC_lub27, smoothMPC_sub27, smoothMPC_riub27, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_11(params->lb29, smoothMPC_z28, smoothMPC_lbIdx28, smoothMPC_llb28, smoothMPC_slb28, smoothMPC_rilb28, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_5(smoothMPC_z28, smoothMPC_ubIdx28, params->ub29, smoothMPC_lub28, smoothMPC_sub28, smoothMPC_riub28, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD3_3(params->lb30, smoothMPC_z29, smoothMPC_lbIdx29, smoothMPC_llb29, smoothMPC_slb29, smoothMPC_rilb29, &info->dgap, &info->res_ineq);
smoothMPC_LA_VSUBADD2_3(smoothMPC_z29, smoothMPC_ubIdx29, params->ub30, smoothMPC_lub29, smoothMPC_sub29, smoothMPC_riub29, &info->dgap, &info->res_ineq);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub00, smoothMPC_sub00, smoothMPC_riub00, smoothMPC_llb00, smoothMPC_slb00, smoothMPC_rilb00, smoothMPC_lbIdx00, smoothMPC_ubIdx00, smoothMPC_grad_ineq00, smoothMPC_lubbysub00, smoothMPC_llbbyslb00);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub01, smoothMPC_sub01, smoothMPC_riub01, smoothMPC_llb01, smoothMPC_slb01, smoothMPC_rilb01, smoothMPC_lbIdx01, smoothMPC_ubIdx01, smoothMPC_grad_ineq01, smoothMPC_lubbysub01, smoothMPC_llbbyslb01);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub02, smoothMPC_sub02, smoothMPC_riub02, smoothMPC_llb02, smoothMPC_slb02, smoothMPC_rilb02, smoothMPC_lbIdx02, smoothMPC_ubIdx02, smoothMPC_grad_ineq02, smoothMPC_lubbysub02, smoothMPC_llbbyslb02);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub03, smoothMPC_sub03, smoothMPC_riub03, smoothMPC_llb03, smoothMPC_slb03, smoothMPC_rilb03, smoothMPC_lbIdx03, smoothMPC_ubIdx03, smoothMPC_grad_ineq03, smoothMPC_lubbysub03, smoothMPC_llbbyslb03);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub04, smoothMPC_sub04, smoothMPC_riub04, smoothMPC_llb04, smoothMPC_slb04, smoothMPC_rilb04, smoothMPC_lbIdx04, smoothMPC_ubIdx04, smoothMPC_grad_ineq04, smoothMPC_lubbysub04, smoothMPC_llbbyslb04);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub05, smoothMPC_sub05, smoothMPC_riub05, smoothMPC_llb05, smoothMPC_slb05, smoothMPC_rilb05, smoothMPC_lbIdx05, smoothMPC_ubIdx05, smoothMPC_grad_ineq05, smoothMPC_lubbysub05, smoothMPC_llbbyslb05);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub06, smoothMPC_sub06, smoothMPC_riub06, smoothMPC_llb06, smoothMPC_slb06, smoothMPC_rilb06, smoothMPC_lbIdx06, smoothMPC_ubIdx06, smoothMPC_grad_ineq06, smoothMPC_lubbysub06, smoothMPC_llbbyslb06);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub07, smoothMPC_sub07, smoothMPC_riub07, smoothMPC_llb07, smoothMPC_slb07, smoothMPC_rilb07, smoothMPC_lbIdx07, smoothMPC_ubIdx07, smoothMPC_grad_ineq07, smoothMPC_lubbysub07, smoothMPC_llbbyslb07);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub08, smoothMPC_sub08, smoothMPC_riub08, smoothMPC_llb08, smoothMPC_slb08, smoothMPC_rilb08, smoothMPC_lbIdx08, smoothMPC_ubIdx08, smoothMPC_grad_ineq08, smoothMPC_lubbysub08, smoothMPC_llbbyslb08);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub09, smoothMPC_sub09, smoothMPC_riub09, smoothMPC_llb09, smoothMPC_slb09, smoothMPC_rilb09, smoothMPC_lbIdx09, smoothMPC_ubIdx09, smoothMPC_grad_ineq09, smoothMPC_lubbysub09, smoothMPC_llbbyslb09);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub10, smoothMPC_sub10, smoothMPC_riub10, smoothMPC_llb10, smoothMPC_slb10, smoothMPC_rilb10, smoothMPC_lbIdx10, smoothMPC_ubIdx10, smoothMPC_grad_ineq10, smoothMPC_lubbysub10, smoothMPC_llbbyslb10);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub11, smoothMPC_sub11, smoothMPC_riub11, smoothMPC_llb11, smoothMPC_slb11, smoothMPC_rilb11, smoothMPC_lbIdx11, smoothMPC_ubIdx11, smoothMPC_grad_ineq11, smoothMPC_lubbysub11, smoothMPC_llbbyslb11);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub12, smoothMPC_sub12, smoothMPC_riub12, smoothMPC_llb12, smoothMPC_slb12, smoothMPC_rilb12, smoothMPC_lbIdx12, smoothMPC_ubIdx12, smoothMPC_grad_ineq12, smoothMPC_lubbysub12, smoothMPC_llbbyslb12);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub13, smoothMPC_sub13, smoothMPC_riub13, smoothMPC_llb13, smoothMPC_slb13, smoothMPC_rilb13, smoothMPC_lbIdx13, smoothMPC_ubIdx13, smoothMPC_grad_ineq13, smoothMPC_lubbysub13, smoothMPC_llbbyslb13);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub14, smoothMPC_sub14, smoothMPC_riub14, smoothMPC_llb14, smoothMPC_slb14, smoothMPC_rilb14, smoothMPC_lbIdx14, smoothMPC_ubIdx14, smoothMPC_grad_ineq14, smoothMPC_lubbysub14, smoothMPC_llbbyslb14);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub15, smoothMPC_sub15, smoothMPC_riub15, smoothMPC_llb15, smoothMPC_slb15, smoothMPC_rilb15, smoothMPC_lbIdx15, smoothMPC_ubIdx15, smoothMPC_grad_ineq15, smoothMPC_lubbysub15, smoothMPC_llbbyslb15);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub16, smoothMPC_sub16, smoothMPC_riub16, smoothMPC_llb16, smoothMPC_slb16, smoothMPC_rilb16, smoothMPC_lbIdx16, smoothMPC_ubIdx16, smoothMPC_grad_ineq16, smoothMPC_lubbysub16, smoothMPC_llbbyslb16);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub17, smoothMPC_sub17, smoothMPC_riub17, smoothMPC_llb17, smoothMPC_slb17, smoothMPC_rilb17, smoothMPC_lbIdx17, smoothMPC_ubIdx17, smoothMPC_grad_ineq17, smoothMPC_lubbysub17, smoothMPC_llbbyslb17);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub18, smoothMPC_sub18, smoothMPC_riub18, smoothMPC_llb18, smoothMPC_slb18, smoothMPC_rilb18, smoothMPC_lbIdx18, smoothMPC_ubIdx18, smoothMPC_grad_ineq18, smoothMPC_lubbysub18, smoothMPC_llbbyslb18);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub19, smoothMPC_sub19, smoothMPC_riub19, smoothMPC_llb19, smoothMPC_slb19, smoothMPC_rilb19, smoothMPC_lbIdx19, smoothMPC_ubIdx19, smoothMPC_grad_ineq19, smoothMPC_lubbysub19, smoothMPC_llbbyslb19);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub20, smoothMPC_sub20, smoothMPC_riub20, smoothMPC_llb20, smoothMPC_slb20, smoothMPC_rilb20, smoothMPC_lbIdx20, smoothMPC_ubIdx20, smoothMPC_grad_ineq20, smoothMPC_lubbysub20, smoothMPC_llbbyslb20);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub21, smoothMPC_sub21, smoothMPC_riub21, smoothMPC_llb21, smoothMPC_slb21, smoothMPC_rilb21, smoothMPC_lbIdx21, smoothMPC_ubIdx21, smoothMPC_grad_ineq21, smoothMPC_lubbysub21, smoothMPC_llbbyslb21);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub22, smoothMPC_sub22, smoothMPC_riub22, smoothMPC_llb22, smoothMPC_slb22, smoothMPC_rilb22, smoothMPC_lbIdx22, smoothMPC_ubIdx22, smoothMPC_grad_ineq22, smoothMPC_lubbysub22, smoothMPC_llbbyslb22);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub23, smoothMPC_sub23, smoothMPC_riub23, smoothMPC_llb23, smoothMPC_slb23, smoothMPC_rilb23, smoothMPC_lbIdx23, smoothMPC_ubIdx23, smoothMPC_grad_ineq23, smoothMPC_lubbysub23, smoothMPC_llbbyslb23);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub24, smoothMPC_sub24, smoothMPC_riub24, smoothMPC_llb24, smoothMPC_slb24, smoothMPC_rilb24, smoothMPC_lbIdx24, smoothMPC_ubIdx24, smoothMPC_grad_ineq24, smoothMPC_lubbysub24, smoothMPC_llbbyslb24);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub25, smoothMPC_sub25, smoothMPC_riub25, smoothMPC_llb25, smoothMPC_slb25, smoothMPC_rilb25, smoothMPC_lbIdx25, smoothMPC_ubIdx25, smoothMPC_grad_ineq25, smoothMPC_lubbysub25, smoothMPC_llbbyslb25);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub26, smoothMPC_sub26, smoothMPC_riub26, smoothMPC_llb26, smoothMPC_slb26, smoothMPC_rilb26, smoothMPC_lbIdx26, smoothMPC_ubIdx26, smoothMPC_grad_ineq26, smoothMPC_lubbysub26, smoothMPC_llbbyslb26);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub27, smoothMPC_sub27, smoothMPC_riub27, smoothMPC_llb27, smoothMPC_slb27, smoothMPC_rilb27, smoothMPC_lbIdx27, smoothMPC_ubIdx27, smoothMPC_grad_ineq27, smoothMPC_lubbysub27, smoothMPC_llbbyslb27);
smoothMPC_LA_INEQ_B_GRAD_11_11_5(smoothMPC_lub28, smoothMPC_sub28, smoothMPC_riub28, smoothMPC_llb28, smoothMPC_slb28, smoothMPC_rilb28, smoothMPC_lbIdx28, smoothMPC_ubIdx28, smoothMPC_grad_ineq28, smoothMPC_lubbysub28, smoothMPC_llbbyslb28);
smoothMPC_LA_INEQ_B_GRAD_3_3_3(smoothMPC_lub29, smoothMPC_sub29, smoothMPC_riub29, smoothMPC_llb29, smoothMPC_slb29, smoothMPC_rilb29, smoothMPC_lbIdx29, smoothMPC_ubIdx29, smoothMPC_grad_ineq29, smoothMPC_lubbysub29, smoothMPC_llbbyslb29);
info->dobj = info->pobj - info->dgap;
info->rdgap = info->pobj ? info->dgap / info->pobj : 1e6;
if( info->rdgap < 0 ) info->rdgap = -info->rdgap;
PRINTTEXT("  %3d  %3.1e  %3.1e  %+6.4e  %+6.4e  %+3.1e  %3.1e  %3.1e\n",info->it, info->res_eq, info->res_ineq, info->pobj, info->dobj, info->dgap, info->rdgap, info->mu);
if( info->mu < smoothMPC_SET_ACC_KKTCOMPL
    && (info->rdgap < smoothMPC_SET_ACC_RDGAP || info->dgap < smoothMPC_SET_ACC_KKTCOMPL)
    && info->res_eq < smoothMPC_SET_ACC_RESEQ
    && info->res_ineq < smoothMPC_SET_ACC_RESINEQ ){
PRINTTEXT("OPTIMAL (within RESEQ=%2.1e, RESINEQ=%2.1e, (R)DGAP=(%2.1e)%2.1e, MU=%2.1e).\n",smoothMPC_SET_ACC_RESEQ, smoothMPC_SET_ACC_RESINEQ,smoothMPC_SET_ACC_KKTCOMPL,smoothMPC_SET_ACC_RDGAP,smoothMPC_SET_ACC_KKTCOMPL);
exitcode = smoothMPC_OPTIMAL; break; }
if( info->it == smoothMPC_SET_MAXIT ){
PRINTTEXT("Maximum number of iterations reached, exiting.\n");
exitcode = smoothMPC_MAXITREACHED; break; }
smoothMPC_LA_VVADD3_322(smoothMPC_grad_cost, smoothMPC_grad_eq, smoothMPC_grad_ineq, smoothMPC_rd);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H1, smoothMPC_llbbyslb00, smoothMPC_lbIdx00, smoothMPC_lubbysub00, smoothMPC_ubIdx00, smoothMPC_Phi00);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi00, params->C1, smoothMPC_V00);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi00, smoothMPC_D00, smoothMPC_W00);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W00, smoothMPC_V00, smoothMPC_Ysd01);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi00, smoothMPC_rd00, smoothMPC_Lbyrd00);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H2, smoothMPC_llbbyslb01, smoothMPC_lbIdx01, smoothMPC_lubbysub01, smoothMPC_ubIdx01, smoothMPC_Phi01);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi01, params->C2, smoothMPC_V01);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi01, smoothMPC_D01, smoothMPC_W01);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W01, smoothMPC_V01, smoothMPC_Ysd02);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi01, smoothMPC_rd01, smoothMPC_Lbyrd01);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H3, smoothMPC_llbbyslb02, smoothMPC_lbIdx02, smoothMPC_lubbysub02, smoothMPC_ubIdx02, smoothMPC_Phi02);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi02, params->C3, smoothMPC_V02);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi02, smoothMPC_D01, smoothMPC_W02);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W02, smoothMPC_V02, smoothMPC_Ysd03);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi02, smoothMPC_rd02, smoothMPC_Lbyrd02);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H4, smoothMPC_llbbyslb03, smoothMPC_lbIdx03, smoothMPC_lubbysub03, smoothMPC_ubIdx03, smoothMPC_Phi03);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi03, params->C4, smoothMPC_V03);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi03, smoothMPC_D01, smoothMPC_W03);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W03, smoothMPC_V03, smoothMPC_Ysd04);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi03, smoothMPC_rd03, smoothMPC_Lbyrd03);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H5, smoothMPC_llbbyslb04, smoothMPC_lbIdx04, smoothMPC_lubbysub04, smoothMPC_ubIdx04, smoothMPC_Phi04);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi04, params->C5, smoothMPC_V04);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi04, smoothMPC_D01, smoothMPC_W04);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W04, smoothMPC_V04, smoothMPC_Ysd05);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi04, smoothMPC_rd04, smoothMPC_Lbyrd04);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H6, smoothMPC_llbbyslb05, smoothMPC_lbIdx05, smoothMPC_lubbysub05, smoothMPC_ubIdx05, smoothMPC_Phi05);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi05, params->C6, smoothMPC_V05);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi05, smoothMPC_D01, smoothMPC_W05);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W05, smoothMPC_V05, smoothMPC_Ysd06);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi05, smoothMPC_rd05, smoothMPC_Lbyrd05);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H7, smoothMPC_llbbyslb06, smoothMPC_lbIdx06, smoothMPC_lubbysub06, smoothMPC_ubIdx06, smoothMPC_Phi06);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi06, params->C7, smoothMPC_V06);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi06, smoothMPC_D01, smoothMPC_W06);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W06, smoothMPC_V06, smoothMPC_Ysd07);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi06, smoothMPC_rd06, smoothMPC_Lbyrd06);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H8, smoothMPC_llbbyslb07, smoothMPC_lbIdx07, smoothMPC_lubbysub07, smoothMPC_ubIdx07, smoothMPC_Phi07);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi07, params->C8, smoothMPC_V07);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi07, smoothMPC_D01, smoothMPC_W07);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W07, smoothMPC_V07, smoothMPC_Ysd08);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi07, smoothMPC_rd07, smoothMPC_Lbyrd07);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H9, smoothMPC_llbbyslb08, smoothMPC_lbIdx08, smoothMPC_lubbysub08, smoothMPC_ubIdx08, smoothMPC_Phi08);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi08, params->C9, smoothMPC_V08);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi08, smoothMPC_D01, smoothMPC_W08);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W08, smoothMPC_V08, smoothMPC_Ysd09);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi08, smoothMPC_rd08, smoothMPC_Lbyrd08);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H10, smoothMPC_llbbyslb09, smoothMPC_lbIdx09, smoothMPC_lubbysub09, smoothMPC_ubIdx09, smoothMPC_Phi09);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi09, params->C10, smoothMPC_V09);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi09, smoothMPC_D01, smoothMPC_W09);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W09, smoothMPC_V09, smoothMPC_Ysd10);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi09, smoothMPC_rd09, smoothMPC_Lbyrd09);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H11, smoothMPC_llbbyslb10, smoothMPC_lbIdx10, smoothMPC_lubbysub10, smoothMPC_ubIdx10, smoothMPC_Phi10);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi10, params->C11, smoothMPC_V10);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi10, smoothMPC_D01, smoothMPC_W10);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W10, smoothMPC_V10, smoothMPC_Ysd11);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi10, smoothMPC_rd10, smoothMPC_Lbyrd10);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H12, smoothMPC_llbbyslb11, smoothMPC_lbIdx11, smoothMPC_lubbysub11, smoothMPC_ubIdx11, smoothMPC_Phi11);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi11, params->C12, smoothMPC_V11);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi11, smoothMPC_D01, smoothMPC_W11);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W11, smoothMPC_V11, smoothMPC_Ysd12);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi11, smoothMPC_rd11, smoothMPC_Lbyrd11);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H13, smoothMPC_llbbyslb12, smoothMPC_lbIdx12, smoothMPC_lubbysub12, smoothMPC_ubIdx12, smoothMPC_Phi12);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi12, params->C13, smoothMPC_V12);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi12, smoothMPC_D01, smoothMPC_W12);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W12, smoothMPC_V12, smoothMPC_Ysd13);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi12, smoothMPC_rd12, smoothMPC_Lbyrd12);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H14, smoothMPC_llbbyslb13, smoothMPC_lbIdx13, smoothMPC_lubbysub13, smoothMPC_ubIdx13, smoothMPC_Phi13);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi13, params->C14, smoothMPC_V13);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi13, smoothMPC_D01, smoothMPC_W13);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W13, smoothMPC_V13, smoothMPC_Ysd14);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi13, smoothMPC_rd13, smoothMPC_Lbyrd13);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H15, smoothMPC_llbbyslb14, smoothMPC_lbIdx14, smoothMPC_lubbysub14, smoothMPC_ubIdx14, smoothMPC_Phi14);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi14, params->C15, smoothMPC_V14);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi14, smoothMPC_D01, smoothMPC_W14);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W14, smoothMPC_V14, smoothMPC_Ysd15);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi14, smoothMPC_rd14, smoothMPC_Lbyrd14);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H16, smoothMPC_llbbyslb15, smoothMPC_lbIdx15, smoothMPC_lubbysub15, smoothMPC_ubIdx15, smoothMPC_Phi15);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi15, params->C16, smoothMPC_V15);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi15, smoothMPC_D01, smoothMPC_W15);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W15, smoothMPC_V15, smoothMPC_Ysd16);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi15, smoothMPC_rd15, smoothMPC_Lbyrd15);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H17, smoothMPC_llbbyslb16, smoothMPC_lbIdx16, smoothMPC_lubbysub16, smoothMPC_ubIdx16, smoothMPC_Phi16);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi16, params->C17, smoothMPC_V16);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi16, smoothMPC_D01, smoothMPC_W16);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W16, smoothMPC_V16, smoothMPC_Ysd17);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi16, smoothMPC_rd16, smoothMPC_Lbyrd16);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H18, smoothMPC_llbbyslb17, smoothMPC_lbIdx17, smoothMPC_lubbysub17, smoothMPC_ubIdx17, smoothMPC_Phi17);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi17, params->C18, smoothMPC_V17);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi17, smoothMPC_D01, smoothMPC_W17);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W17, smoothMPC_V17, smoothMPC_Ysd18);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi17, smoothMPC_rd17, smoothMPC_Lbyrd17);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H19, smoothMPC_llbbyslb18, smoothMPC_lbIdx18, smoothMPC_lubbysub18, smoothMPC_ubIdx18, smoothMPC_Phi18);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi18, params->C19, smoothMPC_V18);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi18, smoothMPC_D01, smoothMPC_W18);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W18, smoothMPC_V18, smoothMPC_Ysd19);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi18, smoothMPC_rd18, smoothMPC_Lbyrd18);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H20, smoothMPC_llbbyslb19, smoothMPC_lbIdx19, smoothMPC_lubbysub19, smoothMPC_ubIdx19, smoothMPC_Phi19);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi19, params->C20, smoothMPC_V19);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi19, smoothMPC_D01, smoothMPC_W19);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W19, smoothMPC_V19, smoothMPC_Ysd20);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi19, smoothMPC_rd19, smoothMPC_Lbyrd19);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H21, smoothMPC_llbbyslb20, smoothMPC_lbIdx20, smoothMPC_lubbysub20, smoothMPC_ubIdx20, smoothMPC_Phi20);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi20, params->C21, smoothMPC_V20);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi20, smoothMPC_D01, smoothMPC_W20);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W20, smoothMPC_V20, smoothMPC_Ysd21);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi20, smoothMPC_rd20, smoothMPC_Lbyrd20);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H22, smoothMPC_llbbyslb21, smoothMPC_lbIdx21, smoothMPC_lubbysub21, smoothMPC_ubIdx21, smoothMPC_Phi21);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi21, params->C22, smoothMPC_V21);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi21, smoothMPC_D01, smoothMPC_W21);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W21, smoothMPC_V21, smoothMPC_Ysd22);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi21, smoothMPC_rd21, smoothMPC_Lbyrd21);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H23, smoothMPC_llbbyslb22, smoothMPC_lbIdx22, smoothMPC_lubbysub22, smoothMPC_ubIdx22, smoothMPC_Phi22);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi22, params->C23, smoothMPC_V22);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi22, smoothMPC_D01, smoothMPC_W22);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W22, smoothMPC_V22, smoothMPC_Ysd23);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi22, smoothMPC_rd22, smoothMPC_Lbyrd22);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H24, smoothMPC_llbbyslb23, smoothMPC_lbIdx23, smoothMPC_lubbysub23, smoothMPC_ubIdx23, smoothMPC_Phi23);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi23, params->C24, smoothMPC_V23);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi23, smoothMPC_D01, smoothMPC_W23);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W23, smoothMPC_V23, smoothMPC_Ysd24);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi23, smoothMPC_rd23, smoothMPC_Lbyrd23);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H25, smoothMPC_llbbyslb24, smoothMPC_lbIdx24, smoothMPC_lubbysub24, smoothMPC_ubIdx24, smoothMPC_Phi24);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi24, params->C25, smoothMPC_V24);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi24, smoothMPC_D01, smoothMPC_W24);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W24, smoothMPC_V24, smoothMPC_Ysd25);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi24, smoothMPC_rd24, smoothMPC_Lbyrd24);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H26, smoothMPC_llbbyslb25, smoothMPC_lbIdx25, smoothMPC_lubbysub25, smoothMPC_ubIdx25, smoothMPC_Phi25);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi25, params->C26, smoothMPC_V25);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi25, smoothMPC_D01, smoothMPC_W25);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W25, smoothMPC_V25, smoothMPC_Ysd26);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi25, smoothMPC_rd25, smoothMPC_Lbyrd25);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H27, smoothMPC_llbbyslb26, smoothMPC_lbIdx26, smoothMPC_lubbysub26, smoothMPC_ubIdx26, smoothMPC_Phi26);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi26, params->C27, smoothMPC_V26);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi26, smoothMPC_D01, smoothMPC_W26);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W26, smoothMPC_V26, smoothMPC_Ysd27);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi26, smoothMPC_rd26, smoothMPC_Lbyrd26);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H28, smoothMPC_llbbyslb27, smoothMPC_lbIdx27, smoothMPC_lubbysub27, smoothMPC_ubIdx27, smoothMPC_Phi27);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi27, params->C28, smoothMPC_V27);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi27, smoothMPC_D01, smoothMPC_W27);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W27, smoothMPC_V27, smoothMPC_Ysd28);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi27, smoothMPC_rd27, smoothMPC_Lbyrd27);
smoothMPC_LA_DIAG_CHOL_LBUB_11_11_5(params->H29, smoothMPC_llbbyslb28, smoothMPC_lbIdx28, smoothMPC_lubbysub28, smoothMPC_ubIdx28, smoothMPC_Phi28);
smoothMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(smoothMPC_Phi28, params->C29, smoothMPC_V28);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(smoothMPC_Phi28, smoothMPC_D01, smoothMPC_W28);
smoothMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(smoothMPC_W28, smoothMPC_V28, smoothMPC_Ysd29);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi28, smoothMPC_rd28, smoothMPC_Lbyrd28);
smoothMPC_LA_DIAG_CHOL_ONELOOP_LBUB_3_3_3(params->H30, smoothMPC_llbbyslb29, smoothMPC_lbIdx29, smoothMPC_lubbysub29, smoothMPC_ubIdx29, smoothMPC_Phi29);
smoothMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_3(smoothMPC_Phi29, smoothMPC_D29, smoothMPC_W29);
smoothMPC_LA_DIAG_FORWARDSUB_3(smoothMPC_Phi29, smoothMPC_rd29, smoothMPC_Lbyrd29);
smoothMPC_LA_DIAGZERO_MMT_3(smoothMPC_W00, smoothMPC_Yd00);
smoothMPC_LA_DIAGZERO_MVMSUB7_3(smoothMPC_W00, smoothMPC_Lbyrd00, smoothMPC_re00, smoothMPC_beta00);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V00, smoothMPC_W01, smoothMPC_Yd01);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V00, smoothMPC_Lbyrd00, smoothMPC_W01, smoothMPC_Lbyrd01, smoothMPC_re01, smoothMPC_beta01);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V01, smoothMPC_W02, smoothMPC_Yd02);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V01, smoothMPC_Lbyrd01, smoothMPC_W02, smoothMPC_Lbyrd02, smoothMPC_re02, smoothMPC_beta02);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V02, smoothMPC_W03, smoothMPC_Yd03);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V02, smoothMPC_Lbyrd02, smoothMPC_W03, smoothMPC_Lbyrd03, smoothMPC_re03, smoothMPC_beta03);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V03, smoothMPC_W04, smoothMPC_Yd04);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V03, smoothMPC_Lbyrd03, smoothMPC_W04, smoothMPC_Lbyrd04, smoothMPC_re04, smoothMPC_beta04);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V04, smoothMPC_W05, smoothMPC_Yd05);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V04, smoothMPC_Lbyrd04, smoothMPC_W05, smoothMPC_Lbyrd05, smoothMPC_re05, smoothMPC_beta05);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V05, smoothMPC_W06, smoothMPC_Yd06);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V05, smoothMPC_Lbyrd05, smoothMPC_W06, smoothMPC_Lbyrd06, smoothMPC_re06, smoothMPC_beta06);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V06, smoothMPC_W07, smoothMPC_Yd07);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V06, smoothMPC_Lbyrd06, smoothMPC_W07, smoothMPC_Lbyrd07, smoothMPC_re07, smoothMPC_beta07);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V07, smoothMPC_W08, smoothMPC_Yd08);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V07, smoothMPC_Lbyrd07, smoothMPC_W08, smoothMPC_Lbyrd08, smoothMPC_re08, smoothMPC_beta08);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V08, smoothMPC_W09, smoothMPC_Yd09);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V08, smoothMPC_Lbyrd08, smoothMPC_W09, smoothMPC_Lbyrd09, smoothMPC_re09, smoothMPC_beta09);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V09, smoothMPC_W10, smoothMPC_Yd10);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V09, smoothMPC_Lbyrd09, smoothMPC_W10, smoothMPC_Lbyrd10, smoothMPC_re10, smoothMPC_beta10);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V10, smoothMPC_W11, smoothMPC_Yd11);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V10, smoothMPC_Lbyrd10, smoothMPC_W11, smoothMPC_Lbyrd11, smoothMPC_re11, smoothMPC_beta11);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V11, smoothMPC_W12, smoothMPC_Yd12);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V11, smoothMPC_Lbyrd11, smoothMPC_W12, smoothMPC_Lbyrd12, smoothMPC_re12, smoothMPC_beta12);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V12, smoothMPC_W13, smoothMPC_Yd13);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V12, smoothMPC_Lbyrd12, smoothMPC_W13, smoothMPC_Lbyrd13, smoothMPC_re13, smoothMPC_beta13);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V13, smoothMPC_W14, smoothMPC_Yd14);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V13, smoothMPC_Lbyrd13, smoothMPC_W14, smoothMPC_Lbyrd14, smoothMPC_re14, smoothMPC_beta14);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V14, smoothMPC_W15, smoothMPC_Yd15);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V14, smoothMPC_Lbyrd14, smoothMPC_W15, smoothMPC_Lbyrd15, smoothMPC_re15, smoothMPC_beta15);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V15, smoothMPC_W16, smoothMPC_Yd16);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V15, smoothMPC_Lbyrd15, smoothMPC_W16, smoothMPC_Lbyrd16, smoothMPC_re16, smoothMPC_beta16);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V16, smoothMPC_W17, smoothMPC_Yd17);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V16, smoothMPC_Lbyrd16, smoothMPC_W17, smoothMPC_Lbyrd17, smoothMPC_re17, smoothMPC_beta17);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V17, smoothMPC_W18, smoothMPC_Yd18);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V17, smoothMPC_Lbyrd17, smoothMPC_W18, smoothMPC_Lbyrd18, smoothMPC_re18, smoothMPC_beta18);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V18, smoothMPC_W19, smoothMPC_Yd19);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V18, smoothMPC_Lbyrd18, smoothMPC_W19, smoothMPC_Lbyrd19, smoothMPC_re19, smoothMPC_beta19);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V19, smoothMPC_W20, smoothMPC_Yd20);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V19, smoothMPC_Lbyrd19, smoothMPC_W20, smoothMPC_Lbyrd20, smoothMPC_re20, smoothMPC_beta20);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V20, smoothMPC_W21, smoothMPC_Yd21);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V20, smoothMPC_Lbyrd20, smoothMPC_W21, smoothMPC_Lbyrd21, smoothMPC_re21, smoothMPC_beta21);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V21, smoothMPC_W22, smoothMPC_Yd22);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V21, smoothMPC_Lbyrd21, smoothMPC_W22, smoothMPC_Lbyrd22, smoothMPC_re22, smoothMPC_beta22);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V22, smoothMPC_W23, smoothMPC_Yd23);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V22, smoothMPC_Lbyrd22, smoothMPC_W23, smoothMPC_Lbyrd23, smoothMPC_re23, smoothMPC_beta23);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V23, smoothMPC_W24, smoothMPC_Yd24);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V23, smoothMPC_Lbyrd23, smoothMPC_W24, smoothMPC_Lbyrd24, smoothMPC_re24, smoothMPC_beta24);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V24, smoothMPC_W25, smoothMPC_Yd25);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V24, smoothMPC_Lbyrd24, smoothMPC_W25, smoothMPC_Lbyrd25, smoothMPC_re25, smoothMPC_beta25);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V25, smoothMPC_W26, smoothMPC_Yd26);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V25, smoothMPC_Lbyrd25, smoothMPC_W26, smoothMPC_Lbyrd26, smoothMPC_re26, smoothMPC_beta26);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V26, smoothMPC_W27, smoothMPC_Yd27);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V26, smoothMPC_Lbyrd26, smoothMPC_W27, smoothMPC_Lbyrd27, smoothMPC_re27, smoothMPC_beta27);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(smoothMPC_V27, smoothMPC_W28, smoothMPC_Yd28);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(smoothMPC_V27, smoothMPC_Lbyrd27, smoothMPC_W28, smoothMPC_Lbyrd28, smoothMPC_re28, smoothMPC_beta28);
smoothMPC_LA_DENSE_DIAGZERO_MMT2_3_11_3(smoothMPC_V28, smoothMPC_W29, smoothMPC_Yd29);
smoothMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_3(smoothMPC_V28, smoothMPC_Lbyrd28, smoothMPC_W29, smoothMPC_Lbyrd29, smoothMPC_re29, smoothMPC_beta29);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd00, smoothMPC_Ld00);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld00, smoothMPC_beta00, smoothMPC_yy00);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld00, smoothMPC_Ysd01, smoothMPC_Lsd01);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd01, smoothMPC_Yd01);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd01, smoothMPC_Ld01);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd01, smoothMPC_yy00, smoothMPC_beta01, smoothMPC_bmy01);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld01, smoothMPC_bmy01, smoothMPC_yy01);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld01, smoothMPC_Ysd02, smoothMPC_Lsd02);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd02, smoothMPC_Yd02);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd02, smoothMPC_Ld02);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd02, smoothMPC_yy01, smoothMPC_beta02, smoothMPC_bmy02);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld02, smoothMPC_bmy02, smoothMPC_yy02);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld02, smoothMPC_Ysd03, smoothMPC_Lsd03);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd03, smoothMPC_Yd03);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd03, smoothMPC_Ld03);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd03, smoothMPC_yy02, smoothMPC_beta03, smoothMPC_bmy03);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld03, smoothMPC_bmy03, smoothMPC_yy03);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld03, smoothMPC_Ysd04, smoothMPC_Lsd04);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd04, smoothMPC_Yd04);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd04, smoothMPC_Ld04);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd04, smoothMPC_yy03, smoothMPC_beta04, smoothMPC_bmy04);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld04, smoothMPC_bmy04, smoothMPC_yy04);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld04, smoothMPC_Ysd05, smoothMPC_Lsd05);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd05, smoothMPC_Yd05);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd05, smoothMPC_Ld05);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd05, smoothMPC_yy04, smoothMPC_beta05, smoothMPC_bmy05);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld05, smoothMPC_bmy05, smoothMPC_yy05);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld05, smoothMPC_Ysd06, smoothMPC_Lsd06);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd06, smoothMPC_Yd06);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd06, smoothMPC_Ld06);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd06, smoothMPC_yy05, smoothMPC_beta06, smoothMPC_bmy06);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld06, smoothMPC_bmy06, smoothMPC_yy06);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld06, smoothMPC_Ysd07, smoothMPC_Lsd07);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd07, smoothMPC_Yd07);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd07, smoothMPC_Ld07);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd07, smoothMPC_yy06, smoothMPC_beta07, smoothMPC_bmy07);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld07, smoothMPC_bmy07, smoothMPC_yy07);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld07, smoothMPC_Ysd08, smoothMPC_Lsd08);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd08, smoothMPC_Yd08);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd08, smoothMPC_Ld08);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd08, smoothMPC_yy07, smoothMPC_beta08, smoothMPC_bmy08);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld08, smoothMPC_bmy08, smoothMPC_yy08);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld08, smoothMPC_Ysd09, smoothMPC_Lsd09);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd09, smoothMPC_Yd09);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd09, smoothMPC_Ld09);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd09, smoothMPC_yy08, smoothMPC_beta09, smoothMPC_bmy09);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld09, smoothMPC_bmy09, smoothMPC_yy09);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld09, smoothMPC_Ysd10, smoothMPC_Lsd10);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd10, smoothMPC_Yd10);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd10, smoothMPC_Ld10);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd10, smoothMPC_yy09, smoothMPC_beta10, smoothMPC_bmy10);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld10, smoothMPC_bmy10, smoothMPC_yy10);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld10, smoothMPC_Ysd11, smoothMPC_Lsd11);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd11, smoothMPC_Yd11);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd11, smoothMPC_Ld11);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd11, smoothMPC_yy10, smoothMPC_beta11, smoothMPC_bmy11);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld11, smoothMPC_bmy11, smoothMPC_yy11);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld11, smoothMPC_Ysd12, smoothMPC_Lsd12);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd12, smoothMPC_Yd12);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd12, smoothMPC_Ld12);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd12, smoothMPC_yy11, smoothMPC_beta12, smoothMPC_bmy12);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld12, smoothMPC_bmy12, smoothMPC_yy12);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld12, smoothMPC_Ysd13, smoothMPC_Lsd13);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd13, smoothMPC_Yd13);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd13, smoothMPC_Ld13);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd13, smoothMPC_yy12, smoothMPC_beta13, smoothMPC_bmy13);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld13, smoothMPC_bmy13, smoothMPC_yy13);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld13, smoothMPC_Ysd14, smoothMPC_Lsd14);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd14, smoothMPC_Yd14);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd14, smoothMPC_Ld14);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd14, smoothMPC_yy13, smoothMPC_beta14, smoothMPC_bmy14);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld14, smoothMPC_bmy14, smoothMPC_yy14);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld14, smoothMPC_Ysd15, smoothMPC_Lsd15);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd15, smoothMPC_Yd15);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd15, smoothMPC_Ld15);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd15, smoothMPC_yy14, smoothMPC_beta15, smoothMPC_bmy15);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld15, smoothMPC_bmy15, smoothMPC_yy15);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld15, smoothMPC_Ysd16, smoothMPC_Lsd16);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd16, smoothMPC_Yd16);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd16, smoothMPC_Ld16);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd16, smoothMPC_yy15, smoothMPC_beta16, smoothMPC_bmy16);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld16, smoothMPC_bmy16, smoothMPC_yy16);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld16, smoothMPC_Ysd17, smoothMPC_Lsd17);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd17, smoothMPC_Yd17);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd17, smoothMPC_Ld17);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd17, smoothMPC_yy16, smoothMPC_beta17, smoothMPC_bmy17);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld17, smoothMPC_bmy17, smoothMPC_yy17);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld17, smoothMPC_Ysd18, smoothMPC_Lsd18);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd18, smoothMPC_Yd18);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd18, smoothMPC_Ld18);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd18, smoothMPC_yy17, smoothMPC_beta18, smoothMPC_bmy18);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld18, smoothMPC_bmy18, smoothMPC_yy18);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld18, smoothMPC_Ysd19, smoothMPC_Lsd19);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd19, smoothMPC_Yd19);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd19, smoothMPC_Ld19);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd19, smoothMPC_yy18, smoothMPC_beta19, smoothMPC_bmy19);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld19, smoothMPC_bmy19, smoothMPC_yy19);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld19, smoothMPC_Ysd20, smoothMPC_Lsd20);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd20, smoothMPC_Yd20);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd20, smoothMPC_Ld20);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd20, smoothMPC_yy19, smoothMPC_beta20, smoothMPC_bmy20);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld20, smoothMPC_bmy20, smoothMPC_yy20);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld20, smoothMPC_Ysd21, smoothMPC_Lsd21);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd21, smoothMPC_Yd21);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd21, smoothMPC_Ld21);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd21, smoothMPC_yy20, smoothMPC_beta21, smoothMPC_bmy21);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld21, smoothMPC_bmy21, smoothMPC_yy21);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld21, smoothMPC_Ysd22, smoothMPC_Lsd22);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd22, smoothMPC_Yd22);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd22, smoothMPC_Ld22);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd22, smoothMPC_yy21, smoothMPC_beta22, smoothMPC_bmy22);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld22, smoothMPC_bmy22, smoothMPC_yy22);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld22, smoothMPC_Ysd23, smoothMPC_Lsd23);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd23, smoothMPC_Yd23);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd23, smoothMPC_Ld23);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd23, smoothMPC_yy22, smoothMPC_beta23, smoothMPC_bmy23);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld23, smoothMPC_bmy23, smoothMPC_yy23);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld23, smoothMPC_Ysd24, smoothMPC_Lsd24);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd24, smoothMPC_Yd24);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd24, smoothMPC_Ld24);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd24, smoothMPC_yy23, smoothMPC_beta24, smoothMPC_bmy24);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld24, smoothMPC_bmy24, smoothMPC_yy24);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld24, smoothMPC_Ysd25, smoothMPC_Lsd25);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd25, smoothMPC_Yd25);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd25, smoothMPC_Ld25);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd25, smoothMPC_yy24, smoothMPC_beta25, smoothMPC_bmy25);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld25, smoothMPC_bmy25, smoothMPC_yy25);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld25, smoothMPC_Ysd26, smoothMPC_Lsd26);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd26, smoothMPC_Yd26);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd26, smoothMPC_Ld26);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd26, smoothMPC_yy25, smoothMPC_beta26, smoothMPC_bmy26);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld26, smoothMPC_bmy26, smoothMPC_yy26);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld26, smoothMPC_Ysd27, smoothMPC_Lsd27);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd27, smoothMPC_Yd27);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd27, smoothMPC_Ld27);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd27, smoothMPC_yy26, smoothMPC_beta27, smoothMPC_bmy27);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld27, smoothMPC_bmy27, smoothMPC_yy27);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld27, smoothMPC_Ysd28, smoothMPC_Lsd28);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd28, smoothMPC_Yd28);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd28, smoothMPC_Ld28);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd28, smoothMPC_yy27, smoothMPC_beta28, smoothMPC_bmy28);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld28, smoothMPC_bmy28, smoothMPC_yy28);
smoothMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(smoothMPC_Ld28, smoothMPC_Ysd29, smoothMPC_Lsd29);
smoothMPC_LA_DENSE_MMTSUB_3_3(smoothMPC_Lsd29, smoothMPC_Yd29);
smoothMPC_LA_DENSE_CHOL_3(smoothMPC_Yd29, smoothMPC_Ld29);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd29, smoothMPC_yy28, smoothMPC_beta29, smoothMPC_bmy29);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld29, smoothMPC_bmy29, smoothMPC_yy29);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld29, smoothMPC_yy29, smoothMPC_dvaff29);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd29, smoothMPC_dvaff29, smoothMPC_yy28, smoothMPC_bmy28);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld28, smoothMPC_bmy28, smoothMPC_dvaff28);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd28, smoothMPC_dvaff28, smoothMPC_yy27, smoothMPC_bmy27);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld27, smoothMPC_bmy27, smoothMPC_dvaff27);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd27, smoothMPC_dvaff27, smoothMPC_yy26, smoothMPC_bmy26);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld26, smoothMPC_bmy26, smoothMPC_dvaff26);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd26, smoothMPC_dvaff26, smoothMPC_yy25, smoothMPC_bmy25);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld25, smoothMPC_bmy25, smoothMPC_dvaff25);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd25, smoothMPC_dvaff25, smoothMPC_yy24, smoothMPC_bmy24);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld24, smoothMPC_bmy24, smoothMPC_dvaff24);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd24, smoothMPC_dvaff24, smoothMPC_yy23, smoothMPC_bmy23);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld23, smoothMPC_bmy23, smoothMPC_dvaff23);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd23, smoothMPC_dvaff23, smoothMPC_yy22, smoothMPC_bmy22);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld22, smoothMPC_bmy22, smoothMPC_dvaff22);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd22, smoothMPC_dvaff22, smoothMPC_yy21, smoothMPC_bmy21);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld21, smoothMPC_bmy21, smoothMPC_dvaff21);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd21, smoothMPC_dvaff21, smoothMPC_yy20, smoothMPC_bmy20);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld20, smoothMPC_bmy20, smoothMPC_dvaff20);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd20, smoothMPC_dvaff20, smoothMPC_yy19, smoothMPC_bmy19);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld19, smoothMPC_bmy19, smoothMPC_dvaff19);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd19, smoothMPC_dvaff19, smoothMPC_yy18, smoothMPC_bmy18);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld18, smoothMPC_bmy18, smoothMPC_dvaff18);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd18, smoothMPC_dvaff18, smoothMPC_yy17, smoothMPC_bmy17);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld17, smoothMPC_bmy17, smoothMPC_dvaff17);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd17, smoothMPC_dvaff17, smoothMPC_yy16, smoothMPC_bmy16);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld16, smoothMPC_bmy16, smoothMPC_dvaff16);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd16, smoothMPC_dvaff16, smoothMPC_yy15, smoothMPC_bmy15);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld15, smoothMPC_bmy15, smoothMPC_dvaff15);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd15, smoothMPC_dvaff15, smoothMPC_yy14, smoothMPC_bmy14);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld14, smoothMPC_bmy14, smoothMPC_dvaff14);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd14, smoothMPC_dvaff14, smoothMPC_yy13, smoothMPC_bmy13);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld13, smoothMPC_bmy13, smoothMPC_dvaff13);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd13, smoothMPC_dvaff13, smoothMPC_yy12, smoothMPC_bmy12);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld12, smoothMPC_bmy12, smoothMPC_dvaff12);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd12, smoothMPC_dvaff12, smoothMPC_yy11, smoothMPC_bmy11);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld11, smoothMPC_bmy11, smoothMPC_dvaff11);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd11, smoothMPC_dvaff11, smoothMPC_yy10, smoothMPC_bmy10);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld10, smoothMPC_bmy10, smoothMPC_dvaff10);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd10, smoothMPC_dvaff10, smoothMPC_yy09, smoothMPC_bmy09);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld09, smoothMPC_bmy09, smoothMPC_dvaff09);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd09, smoothMPC_dvaff09, smoothMPC_yy08, smoothMPC_bmy08);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld08, smoothMPC_bmy08, smoothMPC_dvaff08);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd08, smoothMPC_dvaff08, smoothMPC_yy07, smoothMPC_bmy07);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld07, smoothMPC_bmy07, smoothMPC_dvaff07);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd07, smoothMPC_dvaff07, smoothMPC_yy06, smoothMPC_bmy06);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld06, smoothMPC_bmy06, smoothMPC_dvaff06);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd06, smoothMPC_dvaff06, smoothMPC_yy05, smoothMPC_bmy05);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld05, smoothMPC_bmy05, smoothMPC_dvaff05);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd05, smoothMPC_dvaff05, smoothMPC_yy04, smoothMPC_bmy04);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld04, smoothMPC_bmy04, smoothMPC_dvaff04);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd04, smoothMPC_dvaff04, smoothMPC_yy03, smoothMPC_bmy03);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld03, smoothMPC_bmy03, smoothMPC_dvaff03);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd03, smoothMPC_dvaff03, smoothMPC_yy02, smoothMPC_bmy02);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld02, smoothMPC_bmy02, smoothMPC_dvaff02);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd02, smoothMPC_dvaff02, smoothMPC_yy01, smoothMPC_bmy01);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld01, smoothMPC_bmy01, smoothMPC_dvaff01);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd01, smoothMPC_dvaff01, smoothMPC_yy00, smoothMPC_bmy00);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld00, smoothMPC_bmy00, smoothMPC_dvaff00);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C1, smoothMPC_dvaff01, smoothMPC_D00, smoothMPC_dvaff00, smoothMPC_grad_eq00);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C2, smoothMPC_dvaff02, smoothMPC_D01, smoothMPC_dvaff01, smoothMPC_grad_eq01);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C3, smoothMPC_dvaff03, smoothMPC_D01, smoothMPC_dvaff02, smoothMPC_grad_eq02);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C4, smoothMPC_dvaff04, smoothMPC_D01, smoothMPC_dvaff03, smoothMPC_grad_eq03);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C5, smoothMPC_dvaff05, smoothMPC_D01, smoothMPC_dvaff04, smoothMPC_grad_eq04);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C6, smoothMPC_dvaff06, smoothMPC_D01, smoothMPC_dvaff05, smoothMPC_grad_eq05);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C7, smoothMPC_dvaff07, smoothMPC_D01, smoothMPC_dvaff06, smoothMPC_grad_eq06);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C8, smoothMPC_dvaff08, smoothMPC_D01, smoothMPC_dvaff07, smoothMPC_grad_eq07);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C9, smoothMPC_dvaff09, smoothMPC_D01, smoothMPC_dvaff08, smoothMPC_grad_eq08);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C10, smoothMPC_dvaff10, smoothMPC_D01, smoothMPC_dvaff09, smoothMPC_grad_eq09);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C11, smoothMPC_dvaff11, smoothMPC_D01, smoothMPC_dvaff10, smoothMPC_grad_eq10);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C12, smoothMPC_dvaff12, smoothMPC_D01, smoothMPC_dvaff11, smoothMPC_grad_eq11);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C13, smoothMPC_dvaff13, smoothMPC_D01, smoothMPC_dvaff12, smoothMPC_grad_eq12);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C14, smoothMPC_dvaff14, smoothMPC_D01, smoothMPC_dvaff13, smoothMPC_grad_eq13);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C15, smoothMPC_dvaff15, smoothMPC_D01, smoothMPC_dvaff14, smoothMPC_grad_eq14);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C16, smoothMPC_dvaff16, smoothMPC_D01, smoothMPC_dvaff15, smoothMPC_grad_eq15);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C17, smoothMPC_dvaff17, smoothMPC_D01, smoothMPC_dvaff16, smoothMPC_grad_eq16);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C18, smoothMPC_dvaff18, smoothMPC_D01, smoothMPC_dvaff17, smoothMPC_grad_eq17);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C19, smoothMPC_dvaff19, smoothMPC_D01, smoothMPC_dvaff18, smoothMPC_grad_eq18);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C20, smoothMPC_dvaff20, smoothMPC_D01, smoothMPC_dvaff19, smoothMPC_grad_eq19);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C21, smoothMPC_dvaff21, smoothMPC_D01, smoothMPC_dvaff20, smoothMPC_grad_eq20);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C22, smoothMPC_dvaff22, smoothMPC_D01, smoothMPC_dvaff21, smoothMPC_grad_eq21);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C23, smoothMPC_dvaff23, smoothMPC_D01, smoothMPC_dvaff22, smoothMPC_grad_eq22);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C24, smoothMPC_dvaff24, smoothMPC_D01, smoothMPC_dvaff23, smoothMPC_grad_eq23);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C25, smoothMPC_dvaff25, smoothMPC_D01, smoothMPC_dvaff24, smoothMPC_grad_eq24);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C26, smoothMPC_dvaff26, smoothMPC_D01, smoothMPC_dvaff25, smoothMPC_grad_eq25);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C27, smoothMPC_dvaff27, smoothMPC_D01, smoothMPC_dvaff26, smoothMPC_grad_eq26);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C28, smoothMPC_dvaff28, smoothMPC_D01, smoothMPC_dvaff27, smoothMPC_grad_eq27);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C29, smoothMPC_dvaff29, smoothMPC_D01, smoothMPC_dvaff28, smoothMPC_grad_eq28);
smoothMPC_LA_DIAGZERO_MTVM_3_3(smoothMPC_D29, smoothMPC_dvaff29, smoothMPC_grad_eq29);
smoothMPC_LA_VSUB2_322(smoothMPC_rd, smoothMPC_grad_eq, smoothMPC_rd);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi00, smoothMPC_rd00, smoothMPC_dzaff00);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi01, smoothMPC_rd01, smoothMPC_dzaff01);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi02, smoothMPC_rd02, smoothMPC_dzaff02);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi03, smoothMPC_rd03, smoothMPC_dzaff03);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi04, smoothMPC_rd04, smoothMPC_dzaff04);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi05, smoothMPC_rd05, smoothMPC_dzaff05);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi06, smoothMPC_rd06, smoothMPC_dzaff06);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi07, smoothMPC_rd07, smoothMPC_dzaff07);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi08, smoothMPC_rd08, smoothMPC_dzaff08);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi09, smoothMPC_rd09, smoothMPC_dzaff09);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi10, smoothMPC_rd10, smoothMPC_dzaff10);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi11, smoothMPC_rd11, smoothMPC_dzaff11);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi12, smoothMPC_rd12, smoothMPC_dzaff12);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi13, smoothMPC_rd13, smoothMPC_dzaff13);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi14, smoothMPC_rd14, smoothMPC_dzaff14);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi15, smoothMPC_rd15, smoothMPC_dzaff15);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi16, smoothMPC_rd16, smoothMPC_dzaff16);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi17, smoothMPC_rd17, smoothMPC_dzaff17);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi18, smoothMPC_rd18, smoothMPC_dzaff18);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi19, smoothMPC_rd19, smoothMPC_dzaff19);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi20, smoothMPC_rd20, smoothMPC_dzaff20);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi21, smoothMPC_rd21, smoothMPC_dzaff21);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi22, smoothMPC_rd22, smoothMPC_dzaff22);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi23, smoothMPC_rd23, smoothMPC_dzaff23);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi24, smoothMPC_rd24, smoothMPC_dzaff24);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi25, smoothMPC_rd25, smoothMPC_dzaff25);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi26, smoothMPC_rd26, smoothMPC_dzaff26);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi27, smoothMPC_rd27, smoothMPC_dzaff27);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi28, smoothMPC_rd28, smoothMPC_dzaff28);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_3(smoothMPC_Phi29, smoothMPC_rd29, smoothMPC_dzaff29);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff00, smoothMPC_lbIdx00, smoothMPC_rilb00, smoothMPC_dslbaff00);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb00, smoothMPC_dslbaff00, smoothMPC_llb00, smoothMPC_dllbaff00);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub00, smoothMPC_dzaff00, smoothMPC_ubIdx00, smoothMPC_dsubaff00);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub00, smoothMPC_dsubaff00, smoothMPC_lub00, smoothMPC_dlubaff00);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff01, smoothMPC_lbIdx01, smoothMPC_rilb01, smoothMPC_dslbaff01);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb01, smoothMPC_dslbaff01, smoothMPC_llb01, smoothMPC_dllbaff01);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub01, smoothMPC_dzaff01, smoothMPC_ubIdx01, smoothMPC_dsubaff01);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub01, smoothMPC_dsubaff01, smoothMPC_lub01, smoothMPC_dlubaff01);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff02, smoothMPC_lbIdx02, smoothMPC_rilb02, smoothMPC_dslbaff02);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb02, smoothMPC_dslbaff02, smoothMPC_llb02, smoothMPC_dllbaff02);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub02, smoothMPC_dzaff02, smoothMPC_ubIdx02, smoothMPC_dsubaff02);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub02, smoothMPC_dsubaff02, smoothMPC_lub02, smoothMPC_dlubaff02);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff03, smoothMPC_lbIdx03, smoothMPC_rilb03, smoothMPC_dslbaff03);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb03, smoothMPC_dslbaff03, smoothMPC_llb03, smoothMPC_dllbaff03);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub03, smoothMPC_dzaff03, smoothMPC_ubIdx03, smoothMPC_dsubaff03);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub03, smoothMPC_dsubaff03, smoothMPC_lub03, smoothMPC_dlubaff03);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff04, smoothMPC_lbIdx04, smoothMPC_rilb04, smoothMPC_dslbaff04);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb04, smoothMPC_dslbaff04, smoothMPC_llb04, smoothMPC_dllbaff04);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub04, smoothMPC_dzaff04, smoothMPC_ubIdx04, smoothMPC_dsubaff04);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub04, smoothMPC_dsubaff04, smoothMPC_lub04, smoothMPC_dlubaff04);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff05, smoothMPC_lbIdx05, smoothMPC_rilb05, smoothMPC_dslbaff05);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb05, smoothMPC_dslbaff05, smoothMPC_llb05, smoothMPC_dllbaff05);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub05, smoothMPC_dzaff05, smoothMPC_ubIdx05, smoothMPC_dsubaff05);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub05, smoothMPC_dsubaff05, smoothMPC_lub05, smoothMPC_dlubaff05);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff06, smoothMPC_lbIdx06, smoothMPC_rilb06, smoothMPC_dslbaff06);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb06, smoothMPC_dslbaff06, smoothMPC_llb06, smoothMPC_dllbaff06);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub06, smoothMPC_dzaff06, smoothMPC_ubIdx06, smoothMPC_dsubaff06);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub06, smoothMPC_dsubaff06, smoothMPC_lub06, smoothMPC_dlubaff06);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff07, smoothMPC_lbIdx07, smoothMPC_rilb07, smoothMPC_dslbaff07);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb07, smoothMPC_dslbaff07, smoothMPC_llb07, smoothMPC_dllbaff07);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub07, smoothMPC_dzaff07, smoothMPC_ubIdx07, smoothMPC_dsubaff07);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub07, smoothMPC_dsubaff07, smoothMPC_lub07, smoothMPC_dlubaff07);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff08, smoothMPC_lbIdx08, smoothMPC_rilb08, smoothMPC_dslbaff08);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb08, smoothMPC_dslbaff08, smoothMPC_llb08, smoothMPC_dllbaff08);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub08, smoothMPC_dzaff08, smoothMPC_ubIdx08, smoothMPC_dsubaff08);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub08, smoothMPC_dsubaff08, smoothMPC_lub08, smoothMPC_dlubaff08);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff09, smoothMPC_lbIdx09, smoothMPC_rilb09, smoothMPC_dslbaff09);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb09, smoothMPC_dslbaff09, smoothMPC_llb09, smoothMPC_dllbaff09);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub09, smoothMPC_dzaff09, smoothMPC_ubIdx09, smoothMPC_dsubaff09);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub09, smoothMPC_dsubaff09, smoothMPC_lub09, smoothMPC_dlubaff09);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff10, smoothMPC_lbIdx10, smoothMPC_rilb10, smoothMPC_dslbaff10);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb10, smoothMPC_dslbaff10, smoothMPC_llb10, smoothMPC_dllbaff10);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub10, smoothMPC_dzaff10, smoothMPC_ubIdx10, smoothMPC_dsubaff10);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub10, smoothMPC_dsubaff10, smoothMPC_lub10, smoothMPC_dlubaff10);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff11, smoothMPC_lbIdx11, smoothMPC_rilb11, smoothMPC_dslbaff11);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb11, smoothMPC_dslbaff11, smoothMPC_llb11, smoothMPC_dllbaff11);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub11, smoothMPC_dzaff11, smoothMPC_ubIdx11, smoothMPC_dsubaff11);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub11, smoothMPC_dsubaff11, smoothMPC_lub11, smoothMPC_dlubaff11);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff12, smoothMPC_lbIdx12, smoothMPC_rilb12, smoothMPC_dslbaff12);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb12, smoothMPC_dslbaff12, smoothMPC_llb12, smoothMPC_dllbaff12);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub12, smoothMPC_dzaff12, smoothMPC_ubIdx12, smoothMPC_dsubaff12);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub12, smoothMPC_dsubaff12, smoothMPC_lub12, smoothMPC_dlubaff12);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff13, smoothMPC_lbIdx13, smoothMPC_rilb13, smoothMPC_dslbaff13);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb13, smoothMPC_dslbaff13, smoothMPC_llb13, smoothMPC_dllbaff13);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub13, smoothMPC_dzaff13, smoothMPC_ubIdx13, smoothMPC_dsubaff13);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub13, smoothMPC_dsubaff13, smoothMPC_lub13, smoothMPC_dlubaff13);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff14, smoothMPC_lbIdx14, smoothMPC_rilb14, smoothMPC_dslbaff14);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb14, smoothMPC_dslbaff14, smoothMPC_llb14, smoothMPC_dllbaff14);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub14, smoothMPC_dzaff14, smoothMPC_ubIdx14, smoothMPC_dsubaff14);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub14, smoothMPC_dsubaff14, smoothMPC_lub14, smoothMPC_dlubaff14);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff15, smoothMPC_lbIdx15, smoothMPC_rilb15, smoothMPC_dslbaff15);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb15, smoothMPC_dslbaff15, smoothMPC_llb15, smoothMPC_dllbaff15);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub15, smoothMPC_dzaff15, smoothMPC_ubIdx15, smoothMPC_dsubaff15);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub15, smoothMPC_dsubaff15, smoothMPC_lub15, smoothMPC_dlubaff15);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff16, smoothMPC_lbIdx16, smoothMPC_rilb16, smoothMPC_dslbaff16);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb16, smoothMPC_dslbaff16, smoothMPC_llb16, smoothMPC_dllbaff16);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub16, smoothMPC_dzaff16, smoothMPC_ubIdx16, smoothMPC_dsubaff16);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub16, smoothMPC_dsubaff16, smoothMPC_lub16, smoothMPC_dlubaff16);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff17, smoothMPC_lbIdx17, smoothMPC_rilb17, smoothMPC_dslbaff17);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb17, smoothMPC_dslbaff17, smoothMPC_llb17, smoothMPC_dllbaff17);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub17, smoothMPC_dzaff17, smoothMPC_ubIdx17, smoothMPC_dsubaff17);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub17, smoothMPC_dsubaff17, smoothMPC_lub17, smoothMPC_dlubaff17);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff18, smoothMPC_lbIdx18, smoothMPC_rilb18, smoothMPC_dslbaff18);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb18, smoothMPC_dslbaff18, smoothMPC_llb18, smoothMPC_dllbaff18);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub18, smoothMPC_dzaff18, smoothMPC_ubIdx18, smoothMPC_dsubaff18);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub18, smoothMPC_dsubaff18, smoothMPC_lub18, smoothMPC_dlubaff18);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff19, smoothMPC_lbIdx19, smoothMPC_rilb19, smoothMPC_dslbaff19);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb19, smoothMPC_dslbaff19, smoothMPC_llb19, smoothMPC_dllbaff19);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub19, smoothMPC_dzaff19, smoothMPC_ubIdx19, smoothMPC_dsubaff19);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub19, smoothMPC_dsubaff19, smoothMPC_lub19, smoothMPC_dlubaff19);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff20, smoothMPC_lbIdx20, smoothMPC_rilb20, smoothMPC_dslbaff20);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb20, smoothMPC_dslbaff20, smoothMPC_llb20, smoothMPC_dllbaff20);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub20, smoothMPC_dzaff20, smoothMPC_ubIdx20, smoothMPC_dsubaff20);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub20, smoothMPC_dsubaff20, smoothMPC_lub20, smoothMPC_dlubaff20);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff21, smoothMPC_lbIdx21, smoothMPC_rilb21, smoothMPC_dslbaff21);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb21, smoothMPC_dslbaff21, smoothMPC_llb21, smoothMPC_dllbaff21);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub21, smoothMPC_dzaff21, smoothMPC_ubIdx21, smoothMPC_dsubaff21);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub21, smoothMPC_dsubaff21, smoothMPC_lub21, smoothMPC_dlubaff21);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff22, smoothMPC_lbIdx22, smoothMPC_rilb22, smoothMPC_dslbaff22);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb22, smoothMPC_dslbaff22, smoothMPC_llb22, smoothMPC_dllbaff22);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub22, smoothMPC_dzaff22, smoothMPC_ubIdx22, smoothMPC_dsubaff22);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub22, smoothMPC_dsubaff22, smoothMPC_lub22, smoothMPC_dlubaff22);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff23, smoothMPC_lbIdx23, smoothMPC_rilb23, smoothMPC_dslbaff23);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb23, smoothMPC_dslbaff23, smoothMPC_llb23, smoothMPC_dllbaff23);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub23, smoothMPC_dzaff23, smoothMPC_ubIdx23, smoothMPC_dsubaff23);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub23, smoothMPC_dsubaff23, smoothMPC_lub23, smoothMPC_dlubaff23);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff24, smoothMPC_lbIdx24, smoothMPC_rilb24, smoothMPC_dslbaff24);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb24, smoothMPC_dslbaff24, smoothMPC_llb24, smoothMPC_dllbaff24);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub24, smoothMPC_dzaff24, smoothMPC_ubIdx24, smoothMPC_dsubaff24);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub24, smoothMPC_dsubaff24, smoothMPC_lub24, smoothMPC_dlubaff24);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff25, smoothMPC_lbIdx25, smoothMPC_rilb25, smoothMPC_dslbaff25);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb25, smoothMPC_dslbaff25, smoothMPC_llb25, smoothMPC_dllbaff25);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub25, smoothMPC_dzaff25, smoothMPC_ubIdx25, smoothMPC_dsubaff25);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub25, smoothMPC_dsubaff25, smoothMPC_lub25, smoothMPC_dlubaff25);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff26, smoothMPC_lbIdx26, smoothMPC_rilb26, smoothMPC_dslbaff26);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb26, smoothMPC_dslbaff26, smoothMPC_llb26, smoothMPC_dllbaff26);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub26, smoothMPC_dzaff26, smoothMPC_ubIdx26, smoothMPC_dsubaff26);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub26, smoothMPC_dsubaff26, smoothMPC_lub26, smoothMPC_dlubaff26);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff27, smoothMPC_lbIdx27, smoothMPC_rilb27, smoothMPC_dslbaff27);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb27, smoothMPC_dslbaff27, smoothMPC_llb27, smoothMPC_dllbaff27);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub27, smoothMPC_dzaff27, smoothMPC_ubIdx27, smoothMPC_dsubaff27);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub27, smoothMPC_dsubaff27, smoothMPC_lub27, smoothMPC_dlubaff27);
smoothMPC_LA_VSUB_INDEXED_11(smoothMPC_dzaff28, smoothMPC_lbIdx28, smoothMPC_rilb28, smoothMPC_dslbaff28);
smoothMPC_LA_VSUB3_11(smoothMPC_llbbyslb28, smoothMPC_dslbaff28, smoothMPC_llb28, smoothMPC_dllbaff28);
smoothMPC_LA_VSUB2_INDEXED_5(smoothMPC_riub28, smoothMPC_dzaff28, smoothMPC_ubIdx28, smoothMPC_dsubaff28);
smoothMPC_LA_VSUB3_5(smoothMPC_lubbysub28, smoothMPC_dsubaff28, smoothMPC_lub28, smoothMPC_dlubaff28);
smoothMPC_LA_VSUB_INDEXED_3(smoothMPC_dzaff29, smoothMPC_lbIdx29, smoothMPC_rilb29, smoothMPC_dslbaff29);
smoothMPC_LA_VSUB3_3(smoothMPC_llbbyslb29, smoothMPC_dslbaff29, smoothMPC_llb29, smoothMPC_dllbaff29);
smoothMPC_LA_VSUB2_INDEXED_3(smoothMPC_riub29, smoothMPC_dzaff29, smoothMPC_ubIdx29, smoothMPC_dsubaff29);
smoothMPC_LA_VSUB3_3(smoothMPC_lubbysub29, smoothMPC_dsubaff29, smoothMPC_lub29, smoothMPC_dlubaff29);
info->lsit_aff = smoothMPC_LINESEARCH_BACKTRACKING_AFFINE(smoothMPC_l, smoothMPC_s, smoothMPC_dl_aff, smoothMPC_ds_aff, &info->step_aff, &info->mu_aff);
if( info->lsit_aff == smoothMPC_NOPROGRESS ){
PRINTTEXT("Affine line search could not proceed at iteration %d.\nThe problem might be infeasible -- exiting.\n",info->it+1);
exitcode = smoothMPC_NOPROGRESS; break;
}
sigma_3rdroot = info->mu_aff / info->mu;
info->sigma = sigma_3rdroot*sigma_3rdroot*sigma_3rdroot;
musigma = info->mu * info->sigma;
smoothMPC_LA_VSUB5_470(smoothMPC_ds_aff, smoothMPC_dl_aff, musigma, smoothMPC_ccrhs);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub00, smoothMPC_sub00, smoothMPC_ubIdx00, smoothMPC_ccrhsl00, smoothMPC_slb00, smoothMPC_lbIdx00, smoothMPC_rd00);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub01, smoothMPC_sub01, smoothMPC_ubIdx01, smoothMPC_ccrhsl01, smoothMPC_slb01, smoothMPC_lbIdx01, smoothMPC_rd01);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi00, smoothMPC_rd00, smoothMPC_Lbyrd00);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi01, smoothMPC_rd01, smoothMPC_Lbyrd01);
smoothMPC_LA_DIAGZERO_MVM_3(smoothMPC_W00, smoothMPC_Lbyrd00, smoothMPC_beta00);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld00, smoothMPC_beta00, smoothMPC_yy00);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V00, smoothMPC_Lbyrd00, smoothMPC_W01, smoothMPC_Lbyrd01, smoothMPC_beta01);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd01, smoothMPC_yy00, smoothMPC_beta01, smoothMPC_bmy01);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld01, smoothMPC_bmy01, smoothMPC_yy01);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub02, smoothMPC_sub02, smoothMPC_ubIdx02, smoothMPC_ccrhsl02, smoothMPC_slb02, smoothMPC_lbIdx02, smoothMPC_rd02);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi02, smoothMPC_rd02, smoothMPC_Lbyrd02);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V01, smoothMPC_Lbyrd01, smoothMPC_W02, smoothMPC_Lbyrd02, smoothMPC_beta02);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd02, smoothMPC_yy01, smoothMPC_beta02, smoothMPC_bmy02);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld02, smoothMPC_bmy02, smoothMPC_yy02);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub03, smoothMPC_sub03, smoothMPC_ubIdx03, smoothMPC_ccrhsl03, smoothMPC_slb03, smoothMPC_lbIdx03, smoothMPC_rd03);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi03, smoothMPC_rd03, smoothMPC_Lbyrd03);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V02, smoothMPC_Lbyrd02, smoothMPC_W03, smoothMPC_Lbyrd03, smoothMPC_beta03);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd03, smoothMPC_yy02, smoothMPC_beta03, smoothMPC_bmy03);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld03, smoothMPC_bmy03, smoothMPC_yy03);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub04, smoothMPC_sub04, smoothMPC_ubIdx04, smoothMPC_ccrhsl04, smoothMPC_slb04, smoothMPC_lbIdx04, smoothMPC_rd04);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi04, smoothMPC_rd04, smoothMPC_Lbyrd04);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V03, smoothMPC_Lbyrd03, smoothMPC_W04, smoothMPC_Lbyrd04, smoothMPC_beta04);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd04, smoothMPC_yy03, smoothMPC_beta04, smoothMPC_bmy04);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld04, smoothMPC_bmy04, smoothMPC_yy04);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub05, smoothMPC_sub05, smoothMPC_ubIdx05, smoothMPC_ccrhsl05, smoothMPC_slb05, smoothMPC_lbIdx05, smoothMPC_rd05);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi05, smoothMPC_rd05, smoothMPC_Lbyrd05);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V04, smoothMPC_Lbyrd04, smoothMPC_W05, smoothMPC_Lbyrd05, smoothMPC_beta05);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd05, smoothMPC_yy04, smoothMPC_beta05, smoothMPC_bmy05);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld05, smoothMPC_bmy05, smoothMPC_yy05);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub06, smoothMPC_sub06, smoothMPC_ubIdx06, smoothMPC_ccrhsl06, smoothMPC_slb06, smoothMPC_lbIdx06, smoothMPC_rd06);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi06, smoothMPC_rd06, smoothMPC_Lbyrd06);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V05, smoothMPC_Lbyrd05, smoothMPC_W06, smoothMPC_Lbyrd06, smoothMPC_beta06);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd06, smoothMPC_yy05, smoothMPC_beta06, smoothMPC_bmy06);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld06, smoothMPC_bmy06, smoothMPC_yy06);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub07, smoothMPC_sub07, smoothMPC_ubIdx07, smoothMPC_ccrhsl07, smoothMPC_slb07, smoothMPC_lbIdx07, smoothMPC_rd07);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi07, smoothMPC_rd07, smoothMPC_Lbyrd07);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V06, smoothMPC_Lbyrd06, smoothMPC_W07, smoothMPC_Lbyrd07, smoothMPC_beta07);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd07, smoothMPC_yy06, smoothMPC_beta07, smoothMPC_bmy07);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld07, smoothMPC_bmy07, smoothMPC_yy07);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub08, smoothMPC_sub08, smoothMPC_ubIdx08, smoothMPC_ccrhsl08, smoothMPC_slb08, smoothMPC_lbIdx08, smoothMPC_rd08);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi08, smoothMPC_rd08, smoothMPC_Lbyrd08);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V07, smoothMPC_Lbyrd07, smoothMPC_W08, smoothMPC_Lbyrd08, smoothMPC_beta08);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd08, smoothMPC_yy07, smoothMPC_beta08, smoothMPC_bmy08);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld08, smoothMPC_bmy08, smoothMPC_yy08);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub09, smoothMPC_sub09, smoothMPC_ubIdx09, smoothMPC_ccrhsl09, smoothMPC_slb09, smoothMPC_lbIdx09, smoothMPC_rd09);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi09, smoothMPC_rd09, smoothMPC_Lbyrd09);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V08, smoothMPC_Lbyrd08, smoothMPC_W09, smoothMPC_Lbyrd09, smoothMPC_beta09);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd09, smoothMPC_yy08, smoothMPC_beta09, smoothMPC_bmy09);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld09, smoothMPC_bmy09, smoothMPC_yy09);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub10, smoothMPC_sub10, smoothMPC_ubIdx10, smoothMPC_ccrhsl10, smoothMPC_slb10, smoothMPC_lbIdx10, smoothMPC_rd10);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi10, smoothMPC_rd10, smoothMPC_Lbyrd10);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V09, smoothMPC_Lbyrd09, smoothMPC_W10, smoothMPC_Lbyrd10, smoothMPC_beta10);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd10, smoothMPC_yy09, smoothMPC_beta10, smoothMPC_bmy10);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld10, smoothMPC_bmy10, smoothMPC_yy10);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub11, smoothMPC_sub11, smoothMPC_ubIdx11, smoothMPC_ccrhsl11, smoothMPC_slb11, smoothMPC_lbIdx11, smoothMPC_rd11);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi11, smoothMPC_rd11, smoothMPC_Lbyrd11);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V10, smoothMPC_Lbyrd10, smoothMPC_W11, smoothMPC_Lbyrd11, smoothMPC_beta11);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd11, smoothMPC_yy10, smoothMPC_beta11, smoothMPC_bmy11);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld11, smoothMPC_bmy11, smoothMPC_yy11);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub12, smoothMPC_sub12, smoothMPC_ubIdx12, smoothMPC_ccrhsl12, smoothMPC_slb12, smoothMPC_lbIdx12, smoothMPC_rd12);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi12, smoothMPC_rd12, smoothMPC_Lbyrd12);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V11, smoothMPC_Lbyrd11, smoothMPC_W12, smoothMPC_Lbyrd12, smoothMPC_beta12);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd12, smoothMPC_yy11, smoothMPC_beta12, smoothMPC_bmy12);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld12, smoothMPC_bmy12, smoothMPC_yy12);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub13, smoothMPC_sub13, smoothMPC_ubIdx13, smoothMPC_ccrhsl13, smoothMPC_slb13, smoothMPC_lbIdx13, smoothMPC_rd13);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi13, smoothMPC_rd13, smoothMPC_Lbyrd13);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V12, smoothMPC_Lbyrd12, smoothMPC_W13, smoothMPC_Lbyrd13, smoothMPC_beta13);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd13, smoothMPC_yy12, smoothMPC_beta13, smoothMPC_bmy13);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld13, smoothMPC_bmy13, smoothMPC_yy13);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub14, smoothMPC_sub14, smoothMPC_ubIdx14, smoothMPC_ccrhsl14, smoothMPC_slb14, smoothMPC_lbIdx14, smoothMPC_rd14);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi14, smoothMPC_rd14, smoothMPC_Lbyrd14);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V13, smoothMPC_Lbyrd13, smoothMPC_W14, smoothMPC_Lbyrd14, smoothMPC_beta14);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd14, smoothMPC_yy13, smoothMPC_beta14, smoothMPC_bmy14);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld14, smoothMPC_bmy14, smoothMPC_yy14);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub15, smoothMPC_sub15, smoothMPC_ubIdx15, smoothMPC_ccrhsl15, smoothMPC_slb15, smoothMPC_lbIdx15, smoothMPC_rd15);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi15, smoothMPC_rd15, smoothMPC_Lbyrd15);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V14, smoothMPC_Lbyrd14, smoothMPC_W15, smoothMPC_Lbyrd15, smoothMPC_beta15);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd15, smoothMPC_yy14, smoothMPC_beta15, smoothMPC_bmy15);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld15, smoothMPC_bmy15, smoothMPC_yy15);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub16, smoothMPC_sub16, smoothMPC_ubIdx16, smoothMPC_ccrhsl16, smoothMPC_slb16, smoothMPC_lbIdx16, smoothMPC_rd16);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi16, smoothMPC_rd16, smoothMPC_Lbyrd16);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V15, smoothMPC_Lbyrd15, smoothMPC_W16, smoothMPC_Lbyrd16, smoothMPC_beta16);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd16, smoothMPC_yy15, smoothMPC_beta16, smoothMPC_bmy16);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld16, smoothMPC_bmy16, smoothMPC_yy16);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub17, smoothMPC_sub17, smoothMPC_ubIdx17, smoothMPC_ccrhsl17, smoothMPC_slb17, smoothMPC_lbIdx17, smoothMPC_rd17);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi17, smoothMPC_rd17, smoothMPC_Lbyrd17);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V16, smoothMPC_Lbyrd16, smoothMPC_W17, smoothMPC_Lbyrd17, smoothMPC_beta17);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd17, smoothMPC_yy16, smoothMPC_beta17, smoothMPC_bmy17);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld17, smoothMPC_bmy17, smoothMPC_yy17);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub18, smoothMPC_sub18, smoothMPC_ubIdx18, smoothMPC_ccrhsl18, smoothMPC_slb18, smoothMPC_lbIdx18, smoothMPC_rd18);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi18, smoothMPC_rd18, smoothMPC_Lbyrd18);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V17, smoothMPC_Lbyrd17, smoothMPC_W18, smoothMPC_Lbyrd18, smoothMPC_beta18);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd18, smoothMPC_yy17, smoothMPC_beta18, smoothMPC_bmy18);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld18, smoothMPC_bmy18, smoothMPC_yy18);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub19, smoothMPC_sub19, smoothMPC_ubIdx19, smoothMPC_ccrhsl19, smoothMPC_slb19, smoothMPC_lbIdx19, smoothMPC_rd19);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi19, smoothMPC_rd19, smoothMPC_Lbyrd19);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V18, smoothMPC_Lbyrd18, smoothMPC_W19, smoothMPC_Lbyrd19, smoothMPC_beta19);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd19, smoothMPC_yy18, smoothMPC_beta19, smoothMPC_bmy19);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld19, smoothMPC_bmy19, smoothMPC_yy19);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub20, smoothMPC_sub20, smoothMPC_ubIdx20, smoothMPC_ccrhsl20, smoothMPC_slb20, smoothMPC_lbIdx20, smoothMPC_rd20);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi20, smoothMPC_rd20, smoothMPC_Lbyrd20);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V19, smoothMPC_Lbyrd19, smoothMPC_W20, smoothMPC_Lbyrd20, smoothMPC_beta20);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd20, smoothMPC_yy19, smoothMPC_beta20, smoothMPC_bmy20);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld20, smoothMPC_bmy20, smoothMPC_yy20);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub21, smoothMPC_sub21, smoothMPC_ubIdx21, smoothMPC_ccrhsl21, smoothMPC_slb21, smoothMPC_lbIdx21, smoothMPC_rd21);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi21, smoothMPC_rd21, smoothMPC_Lbyrd21);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V20, smoothMPC_Lbyrd20, smoothMPC_W21, smoothMPC_Lbyrd21, smoothMPC_beta21);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd21, smoothMPC_yy20, smoothMPC_beta21, smoothMPC_bmy21);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld21, smoothMPC_bmy21, smoothMPC_yy21);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub22, smoothMPC_sub22, smoothMPC_ubIdx22, smoothMPC_ccrhsl22, smoothMPC_slb22, smoothMPC_lbIdx22, smoothMPC_rd22);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi22, smoothMPC_rd22, smoothMPC_Lbyrd22);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V21, smoothMPC_Lbyrd21, smoothMPC_W22, smoothMPC_Lbyrd22, smoothMPC_beta22);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd22, smoothMPC_yy21, smoothMPC_beta22, smoothMPC_bmy22);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld22, smoothMPC_bmy22, smoothMPC_yy22);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub23, smoothMPC_sub23, smoothMPC_ubIdx23, smoothMPC_ccrhsl23, smoothMPC_slb23, smoothMPC_lbIdx23, smoothMPC_rd23);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi23, smoothMPC_rd23, smoothMPC_Lbyrd23);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V22, smoothMPC_Lbyrd22, smoothMPC_W23, smoothMPC_Lbyrd23, smoothMPC_beta23);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd23, smoothMPC_yy22, smoothMPC_beta23, smoothMPC_bmy23);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld23, smoothMPC_bmy23, smoothMPC_yy23);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub24, smoothMPC_sub24, smoothMPC_ubIdx24, smoothMPC_ccrhsl24, smoothMPC_slb24, smoothMPC_lbIdx24, smoothMPC_rd24);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi24, smoothMPC_rd24, smoothMPC_Lbyrd24);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V23, smoothMPC_Lbyrd23, smoothMPC_W24, smoothMPC_Lbyrd24, smoothMPC_beta24);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd24, smoothMPC_yy23, smoothMPC_beta24, smoothMPC_bmy24);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld24, smoothMPC_bmy24, smoothMPC_yy24);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub25, smoothMPC_sub25, smoothMPC_ubIdx25, smoothMPC_ccrhsl25, smoothMPC_slb25, smoothMPC_lbIdx25, smoothMPC_rd25);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi25, smoothMPC_rd25, smoothMPC_Lbyrd25);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V24, smoothMPC_Lbyrd24, smoothMPC_W25, smoothMPC_Lbyrd25, smoothMPC_beta25);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd25, smoothMPC_yy24, smoothMPC_beta25, smoothMPC_bmy25);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld25, smoothMPC_bmy25, smoothMPC_yy25);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub26, smoothMPC_sub26, smoothMPC_ubIdx26, smoothMPC_ccrhsl26, smoothMPC_slb26, smoothMPC_lbIdx26, smoothMPC_rd26);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi26, smoothMPC_rd26, smoothMPC_Lbyrd26);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V25, smoothMPC_Lbyrd25, smoothMPC_W26, smoothMPC_Lbyrd26, smoothMPC_beta26);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd26, smoothMPC_yy25, smoothMPC_beta26, smoothMPC_bmy26);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld26, smoothMPC_bmy26, smoothMPC_yy26);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub27, smoothMPC_sub27, smoothMPC_ubIdx27, smoothMPC_ccrhsl27, smoothMPC_slb27, smoothMPC_lbIdx27, smoothMPC_rd27);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi27, smoothMPC_rd27, smoothMPC_Lbyrd27);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V26, smoothMPC_Lbyrd26, smoothMPC_W27, smoothMPC_Lbyrd27, smoothMPC_beta27);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd27, smoothMPC_yy26, smoothMPC_beta27, smoothMPC_bmy27);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld27, smoothMPC_bmy27, smoothMPC_yy27);
smoothMPC_LA_VSUB6_INDEXED_11_5_11(smoothMPC_ccrhsub28, smoothMPC_sub28, smoothMPC_ubIdx28, smoothMPC_ccrhsl28, smoothMPC_slb28, smoothMPC_lbIdx28, smoothMPC_rd28);
smoothMPC_LA_DIAG_FORWARDSUB_11(smoothMPC_Phi28, smoothMPC_rd28, smoothMPC_Lbyrd28);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(smoothMPC_V27, smoothMPC_Lbyrd27, smoothMPC_W28, smoothMPC_Lbyrd28, smoothMPC_beta28);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd28, smoothMPC_yy27, smoothMPC_beta28, smoothMPC_bmy28);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld28, smoothMPC_bmy28, smoothMPC_yy28);
smoothMPC_LA_VSUB6_INDEXED_3_3_3(smoothMPC_ccrhsub29, smoothMPC_sub29, smoothMPC_ubIdx29, smoothMPC_ccrhsl29, smoothMPC_slb29, smoothMPC_lbIdx29, smoothMPC_rd29);
smoothMPC_LA_DIAG_FORWARDSUB_3(smoothMPC_Phi29, smoothMPC_rd29, smoothMPC_Lbyrd29);
smoothMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_3(smoothMPC_V28, smoothMPC_Lbyrd28, smoothMPC_W29, smoothMPC_Lbyrd29, smoothMPC_beta29);
smoothMPC_LA_DENSE_MVMSUB1_3_3(smoothMPC_Lsd29, smoothMPC_yy28, smoothMPC_beta29, smoothMPC_bmy29);
smoothMPC_LA_DENSE_FORWARDSUB_3(smoothMPC_Ld29, smoothMPC_bmy29, smoothMPC_yy29);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld29, smoothMPC_yy29, smoothMPC_dvcc29);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd29, smoothMPC_dvcc29, smoothMPC_yy28, smoothMPC_bmy28);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld28, smoothMPC_bmy28, smoothMPC_dvcc28);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd28, smoothMPC_dvcc28, smoothMPC_yy27, smoothMPC_bmy27);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld27, smoothMPC_bmy27, smoothMPC_dvcc27);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd27, smoothMPC_dvcc27, smoothMPC_yy26, smoothMPC_bmy26);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld26, smoothMPC_bmy26, smoothMPC_dvcc26);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd26, smoothMPC_dvcc26, smoothMPC_yy25, smoothMPC_bmy25);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld25, smoothMPC_bmy25, smoothMPC_dvcc25);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd25, smoothMPC_dvcc25, smoothMPC_yy24, smoothMPC_bmy24);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld24, smoothMPC_bmy24, smoothMPC_dvcc24);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd24, smoothMPC_dvcc24, smoothMPC_yy23, smoothMPC_bmy23);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld23, smoothMPC_bmy23, smoothMPC_dvcc23);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd23, smoothMPC_dvcc23, smoothMPC_yy22, smoothMPC_bmy22);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld22, smoothMPC_bmy22, smoothMPC_dvcc22);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd22, smoothMPC_dvcc22, smoothMPC_yy21, smoothMPC_bmy21);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld21, smoothMPC_bmy21, smoothMPC_dvcc21);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd21, smoothMPC_dvcc21, smoothMPC_yy20, smoothMPC_bmy20);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld20, smoothMPC_bmy20, smoothMPC_dvcc20);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd20, smoothMPC_dvcc20, smoothMPC_yy19, smoothMPC_bmy19);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld19, smoothMPC_bmy19, smoothMPC_dvcc19);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd19, smoothMPC_dvcc19, smoothMPC_yy18, smoothMPC_bmy18);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld18, smoothMPC_bmy18, smoothMPC_dvcc18);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd18, smoothMPC_dvcc18, smoothMPC_yy17, smoothMPC_bmy17);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld17, smoothMPC_bmy17, smoothMPC_dvcc17);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd17, smoothMPC_dvcc17, smoothMPC_yy16, smoothMPC_bmy16);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld16, smoothMPC_bmy16, smoothMPC_dvcc16);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd16, smoothMPC_dvcc16, smoothMPC_yy15, smoothMPC_bmy15);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld15, smoothMPC_bmy15, smoothMPC_dvcc15);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd15, smoothMPC_dvcc15, smoothMPC_yy14, smoothMPC_bmy14);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld14, smoothMPC_bmy14, smoothMPC_dvcc14);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd14, smoothMPC_dvcc14, smoothMPC_yy13, smoothMPC_bmy13);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld13, smoothMPC_bmy13, smoothMPC_dvcc13);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd13, smoothMPC_dvcc13, smoothMPC_yy12, smoothMPC_bmy12);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld12, smoothMPC_bmy12, smoothMPC_dvcc12);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd12, smoothMPC_dvcc12, smoothMPC_yy11, smoothMPC_bmy11);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld11, smoothMPC_bmy11, smoothMPC_dvcc11);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd11, smoothMPC_dvcc11, smoothMPC_yy10, smoothMPC_bmy10);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld10, smoothMPC_bmy10, smoothMPC_dvcc10);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd10, smoothMPC_dvcc10, smoothMPC_yy09, smoothMPC_bmy09);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld09, smoothMPC_bmy09, smoothMPC_dvcc09);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd09, smoothMPC_dvcc09, smoothMPC_yy08, smoothMPC_bmy08);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld08, smoothMPC_bmy08, smoothMPC_dvcc08);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd08, smoothMPC_dvcc08, smoothMPC_yy07, smoothMPC_bmy07);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld07, smoothMPC_bmy07, smoothMPC_dvcc07);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd07, smoothMPC_dvcc07, smoothMPC_yy06, smoothMPC_bmy06);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld06, smoothMPC_bmy06, smoothMPC_dvcc06);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd06, smoothMPC_dvcc06, smoothMPC_yy05, smoothMPC_bmy05);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld05, smoothMPC_bmy05, smoothMPC_dvcc05);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd05, smoothMPC_dvcc05, smoothMPC_yy04, smoothMPC_bmy04);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld04, smoothMPC_bmy04, smoothMPC_dvcc04);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd04, smoothMPC_dvcc04, smoothMPC_yy03, smoothMPC_bmy03);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld03, smoothMPC_bmy03, smoothMPC_dvcc03);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd03, smoothMPC_dvcc03, smoothMPC_yy02, smoothMPC_bmy02);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld02, smoothMPC_bmy02, smoothMPC_dvcc02);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd02, smoothMPC_dvcc02, smoothMPC_yy01, smoothMPC_bmy01);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld01, smoothMPC_bmy01, smoothMPC_dvcc01);
smoothMPC_LA_DENSE_MTVMSUB_3_3(smoothMPC_Lsd01, smoothMPC_dvcc01, smoothMPC_yy00, smoothMPC_bmy00);
smoothMPC_LA_DENSE_BACKWARDSUB_3(smoothMPC_Ld00, smoothMPC_bmy00, smoothMPC_dvcc00);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C1, smoothMPC_dvcc01, smoothMPC_D00, smoothMPC_dvcc00, smoothMPC_grad_eq00);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C2, smoothMPC_dvcc02, smoothMPC_D01, smoothMPC_dvcc01, smoothMPC_grad_eq01);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C3, smoothMPC_dvcc03, smoothMPC_D01, smoothMPC_dvcc02, smoothMPC_grad_eq02);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C4, smoothMPC_dvcc04, smoothMPC_D01, smoothMPC_dvcc03, smoothMPC_grad_eq03);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C5, smoothMPC_dvcc05, smoothMPC_D01, smoothMPC_dvcc04, smoothMPC_grad_eq04);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C6, smoothMPC_dvcc06, smoothMPC_D01, smoothMPC_dvcc05, smoothMPC_grad_eq05);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C7, smoothMPC_dvcc07, smoothMPC_D01, smoothMPC_dvcc06, smoothMPC_grad_eq06);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C8, smoothMPC_dvcc08, smoothMPC_D01, smoothMPC_dvcc07, smoothMPC_grad_eq07);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C9, smoothMPC_dvcc09, smoothMPC_D01, smoothMPC_dvcc08, smoothMPC_grad_eq08);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C10, smoothMPC_dvcc10, smoothMPC_D01, smoothMPC_dvcc09, smoothMPC_grad_eq09);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C11, smoothMPC_dvcc11, smoothMPC_D01, smoothMPC_dvcc10, smoothMPC_grad_eq10);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C12, smoothMPC_dvcc12, smoothMPC_D01, smoothMPC_dvcc11, smoothMPC_grad_eq11);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C13, smoothMPC_dvcc13, smoothMPC_D01, smoothMPC_dvcc12, smoothMPC_grad_eq12);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C14, smoothMPC_dvcc14, smoothMPC_D01, smoothMPC_dvcc13, smoothMPC_grad_eq13);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C15, smoothMPC_dvcc15, smoothMPC_D01, smoothMPC_dvcc14, smoothMPC_grad_eq14);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C16, smoothMPC_dvcc16, smoothMPC_D01, smoothMPC_dvcc15, smoothMPC_grad_eq15);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C17, smoothMPC_dvcc17, smoothMPC_D01, smoothMPC_dvcc16, smoothMPC_grad_eq16);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C18, smoothMPC_dvcc18, smoothMPC_D01, smoothMPC_dvcc17, smoothMPC_grad_eq17);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C19, smoothMPC_dvcc19, smoothMPC_D01, smoothMPC_dvcc18, smoothMPC_grad_eq18);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C20, smoothMPC_dvcc20, smoothMPC_D01, smoothMPC_dvcc19, smoothMPC_grad_eq19);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C21, smoothMPC_dvcc21, smoothMPC_D01, smoothMPC_dvcc20, smoothMPC_grad_eq20);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C22, smoothMPC_dvcc22, smoothMPC_D01, smoothMPC_dvcc21, smoothMPC_grad_eq21);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C23, smoothMPC_dvcc23, smoothMPC_D01, smoothMPC_dvcc22, smoothMPC_grad_eq22);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C24, smoothMPC_dvcc24, smoothMPC_D01, smoothMPC_dvcc23, smoothMPC_grad_eq23);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C25, smoothMPC_dvcc25, smoothMPC_D01, smoothMPC_dvcc24, smoothMPC_grad_eq24);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C26, smoothMPC_dvcc26, smoothMPC_D01, smoothMPC_dvcc25, smoothMPC_grad_eq25);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C27, smoothMPC_dvcc27, smoothMPC_D01, smoothMPC_dvcc26, smoothMPC_grad_eq26);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C28, smoothMPC_dvcc28, smoothMPC_D01, smoothMPC_dvcc27, smoothMPC_grad_eq27);
smoothMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C29, smoothMPC_dvcc29, smoothMPC_D01, smoothMPC_dvcc28, smoothMPC_grad_eq28);
smoothMPC_LA_DIAGZERO_MTVM_3_3(smoothMPC_D29, smoothMPC_dvcc29, smoothMPC_grad_eq29);
smoothMPC_LA_VSUB_322(smoothMPC_rd, smoothMPC_grad_eq, smoothMPC_rd);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi00, smoothMPC_rd00, smoothMPC_dzcc00);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi01, smoothMPC_rd01, smoothMPC_dzcc01);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi02, smoothMPC_rd02, smoothMPC_dzcc02);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi03, smoothMPC_rd03, smoothMPC_dzcc03);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi04, smoothMPC_rd04, smoothMPC_dzcc04);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi05, smoothMPC_rd05, smoothMPC_dzcc05);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi06, smoothMPC_rd06, smoothMPC_dzcc06);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi07, smoothMPC_rd07, smoothMPC_dzcc07);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi08, smoothMPC_rd08, smoothMPC_dzcc08);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi09, smoothMPC_rd09, smoothMPC_dzcc09);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi10, smoothMPC_rd10, smoothMPC_dzcc10);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi11, smoothMPC_rd11, smoothMPC_dzcc11);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi12, smoothMPC_rd12, smoothMPC_dzcc12);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi13, smoothMPC_rd13, smoothMPC_dzcc13);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi14, smoothMPC_rd14, smoothMPC_dzcc14);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi15, smoothMPC_rd15, smoothMPC_dzcc15);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi16, smoothMPC_rd16, smoothMPC_dzcc16);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi17, smoothMPC_rd17, smoothMPC_dzcc17);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi18, smoothMPC_rd18, smoothMPC_dzcc18);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi19, smoothMPC_rd19, smoothMPC_dzcc19);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi20, smoothMPC_rd20, smoothMPC_dzcc20);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi21, smoothMPC_rd21, smoothMPC_dzcc21);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi22, smoothMPC_rd22, smoothMPC_dzcc22);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi23, smoothMPC_rd23, smoothMPC_dzcc23);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi24, smoothMPC_rd24, smoothMPC_dzcc24);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi25, smoothMPC_rd25, smoothMPC_dzcc25);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi26, smoothMPC_rd26, smoothMPC_dzcc26);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi27, smoothMPC_rd27, smoothMPC_dzcc27);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_11(smoothMPC_Phi28, smoothMPC_rd28, smoothMPC_dzcc28);
smoothMPC_LA_DIAG_FORWARDBACKWARDSUB_3(smoothMPC_Phi29, smoothMPC_rd29, smoothMPC_dzcc29);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl00, smoothMPC_slb00, smoothMPC_llbbyslb00, smoothMPC_dzcc00, smoothMPC_lbIdx00, smoothMPC_dllbcc00);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub00, smoothMPC_sub00, smoothMPC_lubbysub00, smoothMPC_dzcc00, smoothMPC_ubIdx00, smoothMPC_dlubcc00);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl01, smoothMPC_slb01, smoothMPC_llbbyslb01, smoothMPC_dzcc01, smoothMPC_lbIdx01, smoothMPC_dllbcc01);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub01, smoothMPC_sub01, smoothMPC_lubbysub01, smoothMPC_dzcc01, smoothMPC_ubIdx01, smoothMPC_dlubcc01);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl02, smoothMPC_slb02, smoothMPC_llbbyslb02, smoothMPC_dzcc02, smoothMPC_lbIdx02, smoothMPC_dllbcc02);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub02, smoothMPC_sub02, smoothMPC_lubbysub02, smoothMPC_dzcc02, smoothMPC_ubIdx02, smoothMPC_dlubcc02);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl03, smoothMPC_slb03, smoothMPC_llbbyslb03, smoothMPC_dzcc03, smoothMPC_lbIdx03, smoothMPC_dllbcc03);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub03, smoothMPC_sub03, smoothMPC_lubbysub03, smoothMPC_dzcc03, smoothMPC_ubIdx03, smoothMPC_dlubcc03);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl04, smoothMPC_slb04, smoothMPC_llbbyslb04, smoothMPC_dzcc04, smoothMPC_lbIdx04, smoothMPC_dllbcc04);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub04, smoothMPC_sub04, smoothMPC_lubbysub04, smoothMPC_dzcc04, smoothMPC_ubIdx04, smoothMPC_dlubcc04);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl05, smoothMPC_slb05, smoothMPC_llbbyslb05, smoothMPC_dzcc05, smoothMPC_lbIdx05, smoothMPC_dllbcc05);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub05, smoothMPC_sub05, smoothMPC_lubbysub05, smoothMPC_dzcc05, smoothMPC_ubIdx05, smoothMPC_dlubcc05);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl06, smoothMPC_slb06, smoothMPC_llbbyslb06, smoothMPC_dzcc06, smoothMPC_lbIdx06, smoothMPC_dllbcc06);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub06, smoothMPC_sub06, smoothMPC_lubbysub06, smoothMPC_dzcc06, smoothMPC_ubIdx06, smoothMPC_dlubcc06);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl07, smoothMPC_slb07, smoothMPC_llbbyslb07, smoothMPC_dzcc07, smoothMPC_lbIdx07, smoothMPC_dllbcc07);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub07, smoothMPC_sub07, smoothMPC_lubbysub07, smoothMPC_dzcc07, smoothMPC_ubIdx07, smoothMPC_dlubcc07);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl08, smoothMPC_slb08, smoothMPC_llbbyslb08, smoothMPC_dzcc08, smoothMPC_lbIdx08, smoothMPC_dllbcc08);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub08, smoothMPC_sub08, smoothMPC_lubbysub08, smoothMPC_dzcc08, smoothMPC_ubIdx08, smoothMPC_dlubcc08);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl09, smoothMPC_slb09, smoothMPC_llbbyslb09, smoothMPC_dzcc09, smoothMPC_lbIdx09, smoothMPC_dllbcc09);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub09, smoothMPC_sub09, smoothMPC_lubbysub09, smoothMPC_dzcc09, smoothMPC_ubIdx09, smoothMPC_dlubcc09);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl10, smoothMPC_slb10, smoothMPC_llbbyslb10, smoothMPC_dzcc10, smoothMPC_lbIdx10, smoothMPC_dllbcc10);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub10, smoothMPC_sub10, smoothMPC_lubbysub10, smoothMPC_dzcc10, smoothMPC_ubIdx10, smoothMPC_dlubcc10);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl11, smoothMPC_slb11, smoothMPC_llbbyslb11, smoothMPC_dzcc11, smoothMPC_lbIdx11, smoothMPC_dllbcc11);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub11, smoothMPC_sub11, smoothMPC_lubbysub11, smoothMPC_dzcc11, smoothMPC_ubIdx11, smoothMPC_dlubcc11);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl12, smoothMPC_slb12, smoothMPC_llbbyslb12, smoothMPC_dzcc12, smoothMPC_lbIdx12, smoothMPC_dllbcc12);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub12, smoothMPC_sub12, smoothMPC_lubbysub12, smoothMPC_dzcc12, smoothMPC_ubIdx12, smoothMPC_dlubcc12);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl13, smoothMPC_slb13, smoothMPC_llbbyslb13, smoothMPC_dzcc13, smoothMPC_lbIdx13, smoothMPC_dllbcc13);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub13, smoothMPC_sub13, smoothMPC_lubbysub13, smoothMPC_dzcc13, smoothMPC_ubIdx13, smoothMPC_dlubcc13);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl14, smoothMPC_slb14, smoothMPC_llbbyslb14, smoothMPC_dzcc14, smoothMPC_lbIdx14, smoothMPC_dllbcc14);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub14, smoothMPC_sub14, smoothMPC_lubbysub14, smoothMPC_dzcc14, smoothMPC_ubIdx14, smoothMPC_dlubcc14);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl15, smoothMPC_slb15, smoothMPC_llbbyslb15, smoothMPC_dzcc15, smoothMPC_lbIdx15, smoothMPC_dllbcc15);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub15, smoothMPC_sub15, smoothMPC_lubbysub15, smoothMPC_dzcc15, smoothMPC_ubIdx15, smoothMPC_dlubcc15);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl16, smoothMPC_slb16, smoothMPC_llbbyslb16, smoothMPC_dzcc16, smoothMPC_lbIdx16, smoothMPC_dllbcc16);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub16, smoothMPC_sub16, smoothMPC_lubbysub16, smoothMPC_dzcc16, smoothMPC_ubIdx16, smoothMPC_dlubcc16);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl17, smoothMPC_slb17, smoothMPC_llbbyslb17, smoothMPC_dzcc17, smoothMPC_lbIdx17, smoothMPC_dllbcc17);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub17, smoothMPC_sub17, smoothMPC_lubbysub17, smoothMPC_dzcc17, smoothMPC_ubIdx17, smoothMPC_dlubcc17);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl18, smoothMPC_slb18, smoothMPC_llbbyslb18, smoothMPC_dzcc18, smoothMPC_lbIdx18, smoothMPC_dllbcc18);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub18, smoothMPC_sub18, smoothMPC_lubbysub18, smoothMPC_dzcc18, smoothMPC_ubIdx18, smoothMPC_dlubcc18);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl19, smoothMPC_slb19, smoothMPC_llbbyslb19, smoothMPC_dzcc19, smoothMPC_lbIdx19, smoothMPC_dllbcc19);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub19, smoothMPC_sub19, smoothMPC_lubbysub19, smoothMPC_dzcc19, smoothMPC_ubIdx19, smoothMPC_dlubcc19);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl20, smoothMPC_slb20, smoothMPC_llbbyslb20, smoothMPC_dzcc20, smoothMPC_lbIdx20, smoothMPC_dllbcc20);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub20, smoothMPC_sub20, smoothMPC_lubbysub20, smoothMPC_dzcc20, smoothMPC_ubIdx20, smoothMPC_dlubcc20);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl21, smoothMPC_slb21, smoothMPC_llbbyslb21, smoothMPC_dzcc21, smoothMPC_lbIdx21, smoothMPC_dllbcc21);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub21, smoothMPC_sub21, smoothMPC_lubbysub21, smoothMPC_dzcc21, smoothMPC_ubIdx21, smoothMPC_dlubcc21);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl22, smoothMPC_slb22, smoothMPC_llbbyslb22, smoothMPC_dzcc22, smoothMPC_lbIdx22, smoothMPC_dllbcc22);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub22, smoothMPC_sub22, smoothMPC_lubbysub22, smoothMPC_dzcc22, smoothMPC_ubIdx22, smoothMPC_dlubcc22);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl23, smoothMPC_slb23, smoothMPC_llbbyslb23, smoothMPC_dzcc23, smoothMPC_lbIdx23, smoothMPC_dllbcc23);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub23, smoothMPC_sub23, smoothMPC_lubbysub23, smoothMPC_dzcc23, smoothMPC_ubIdx23, smoothMPC_dlubcc23);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl24, smoothMPC_slb24, smoothMPC_llbbyslb24, smoothMPC_dzcc24, smoothMPC_lbIdx24, smoothMPC_dllbcc24);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub24, smoothMPC_sub24, smoothMPC_lubbysub24, smoothMPC_dzcc24, smoothMPC_ubIdx24, smoothMPC_dlubcc24);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl25, smoothMPC_slb25, smoothMPC_llbbyslb25, smoothMPC_dzcc25, smoothMPC_lbIdx25, smoothMPC_dllbcc25);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub25, smoothMPC_sub25, smoothMPC_lubbysub25, smoothMPC_dzcc25, smoothMPC_ubIdx25, smoothMPC_dlubcc25);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl26, smoothMPC_slb26, smoothMPC_llbbyslb26, smoothMPC_dzcc26, smoothMPC_lbIdx26, smoothMPC_dllbcc26);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub26, smoothMPC_sub26, smoothMPC_lubbysub26, smoothMPC_dzcc26, smoothMPC_ubIdx26, smoothMPC_dlubcc26);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl27, smoothMPC_slb27, smoothMPC_llbbyslb27, smoothMPC_dzcc27, smoothMPC_lbIdx27, smoothMPC_dllbcc27);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub27, smoothMPC_sub27, smoothMPC_lubbysub27, smoothMPC_dzcc27, smoothMPC_ubIdx27, smoothMPC_dlubcc27);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(smoothMPC_ccrhsl28, smoothMPC_slb28, smoothMPC_llbbyslb28, smoothMPC_dzcc28, smoothMPC_lbIdx28, smoothMPC_dllbcc28);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(smoothMPC_ccrhsub28, smoothMPC_sub28, smoothMPC_lubbysub28, smoothMPC_dzcc28, smoothMPC_ubIdx28, smoothMPC_dlubcc28);
smoothMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_3(smoothMPC_ccrhsl29, smoothMPC_slb29, smoothMPC_llbbyslb29, smoothMPC_dzcc29, smoothMPC_lbIdx29, smoothMPC_dllbcc29);
smoothMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_3(smoothMPC_ccrhsub29, smoothMPC_sub29, smoothMPC_lubbysub29, smoothMPC_dzcc29, smoothMPC_ubIdx29, smoothMPC_dlubcc29);
smoothMPC_LA_VSUB7_470(smoothMPC_l, smoothMPC_ccrhs, smoothMPC_s, smoothMPC_dl_cc, smoothMPC_ds_cc);
smoothMPC_LA_VADD_322(smoothMPC_dz_cc, smoothMPC_dz_aff);
smoothMPC_LA_VADD_90(smoothMPC_dv_cc, smoothMPC_dv_aff);
smoothMPC_LA_VADD_470(smoothMPC_dl_cc, smoothMPC_dl_aff);
smoothMPC_LA_VADD_470(smoothMPC_ds_cc, smoothMPC_ds_aff);
info->lsit_cc = smoothMPC_LINESEARCH_BACKTRACKING_COMBINED(smoothMPC_z, smoothMPC_v, smoothMPC_l, smoothMPC_s, smoothMPC_dz_cc, smoothMPC_dv_cc, smoothMPC_dl_cc, smoothMPC_ds_cc, &info->step_cc, &info->mu);
if( info->lsit_cc == smoothMPC_NOPROGRESS ){
PRINTTEXT("Line search could not proceed at iteration %d, exiting.\n",info->it+1);
exitcode = smoothMPC_NOPROGRESS; break;
}
info->it++;
}
output->z1[0] = smoothMPC_z00[0];
output->z1[1] = smoothMPC_z00[1];
output->z1[2] = smoothMPC_z00[2];
output->z1[3] = smoothMPC_z00[3];
output->z1[4] = smoothMPC_z00[4];
output->z2[0] = smoothMPC_z01[0];
output->z2[1] = smoothMPC_z01[1];
output->z2[2] = smoothMPC_z01[2];
output->z2[3] = smoothMPC_z01[3];
output->z2[4] = smoothMPC_z01[4];
output->z3[0] = smoothMPC_z02[0];
output->z3[1] = smoothMPC_z02[1];
output->z3[2] = smoothMPC_z02[2];
output->z3[3] = smoothMPC_z02[3];
output->z3[4] = smoothMPC_z02[4];
output->z4[0] = smoothMPC_z03[0];
output->z4[1] = smoothMPC_z03[1];
output->z4[2] = smoothMPC_z03[2];
output->z4[3] = smoothMPC_z03[3];
output->z4[4] = smoothMPC_z03[4];
output->z5[0] = smoothMPC_z04[0];
output->z5[1] = smoothMPC_z04[1];
output->z5[2] = smoothMPC_z04[2];
output->z5[3] = smoothMPC_z04[3];
output->z5[4] = smoothMPC_z04[4];
output->z6[0] = smoothMPC_z05[0];
output->z6[1] = smoothMPC_z05[1];
output->z6[2] = smoothMPC_z05[2];
output->z6[3] = smoothMPC_z05[3];
output->z6[4] = smoothMPC_z05[4];
output->z7[0] = smoothMPC_z06[0];
output->z7[1] = smoothMPC_z06[1];
output->z7[2] = smoothMPC_z06[2];
output->z7[3] = smoothMPC_z06[3];
output->z7[4] = smoothMPC_z06[4];
output->z8[0] = smoothMPC_z07[0];
output->z8[1] = smoothMPC_z07[1];
output->z8[2] = smoothMPC_z07[2];
output->z8[3] = smoothMPC_z07[3];
output->z8[4] = smoothMPC_z07[4];
output->z9[0] = smoothMPC_z08[0];
output->z9[1] = smoothMPC_z08[1];
output->z9[2] = smoothMPC_z08[2];
output->z9[3] = smoothMPC_z08[3];
output->z9[4] = smoothMPC_z08[4];
output->z10[0] = smoothMPC_z09[0];
output->z10[1] = smoothMPC_z09[1];
output->z10[2] = smoothMPC_z09[2];
output->z10[3] = smoothMPC_z09[3];
output->z10[4] = smoothMPC_z09[4];
output->z11[0] = smoothMPC_z10[0];
output->z11[1] = smoothMPC_z10[1];
output->z11[2] = smoothMPC_z10[2];
output->z11[3] = smoothMPC_z10[3];
output->z11[4] = smoothMPC_z10[4];
output->z12[0] = smoothMPC_z11[0];
output->z12[1] = smoothMPC_z11[1];
output->z12[2] = smoothMPC_z11[2];
output->z12[3] = smoothMPC_z11[3];
output->z12[4] = smoothMPC_z11[4];
output->z13[0] = smoothMPC_z12[0];
output->z13[1] = smoothMPC_z12[1];
output->z13[2] = smoothMPC_z12[2];
output->z13[3] = smoothMPC_z12[3];
output->z13[4] = smoothMPC_z12[4];
output->z14[0] = smoothMPC_z13[0];
output->z14[1] = smoothMPC_z13[1];
output->z14[2] = smoothMPC_z13[2];
output->z14[3] = smoothMPC_z13[3];
output->z14[4] = smoothMPC_z13[4];
output->z15[0] = smoothMPC_z14[0];
output->z15[1] = smoothMPC_z14[1];
output->z15[2] = smoothMPC_z14[2];
output->z15[3] = smoothMPC_z14[3];
output->z15[4] = smoothMPC_z14[4];
output->z16[0] = smoothMPC_z15[0];
output->z16[1] = smoothMPC_z15[1];
output->z16[2] = smoothMPC_z15[2];
output->z16[3] = smoothMPC_z15[3];
output->z16[4] = smoothMPC_z15[4];
output->z17[0] = smoothMPC_z16[0];
output->z17[1] = smoothMPC_z16[1];
output->z17[2] = smoothMPC_z16[2];
output->z17[3] = smoothMPC_z16[3];
output->z17[4] = smoothMPC_z16[4];
output->z18[0] = smoothMPC_z17[0];
output->z18[1] = smoothMPC_z17[1];
output->z18[2] = smoothMPC_z17[2];
output->z18[3] = smoothMPC_z17[3];
output->z18[4] = smoothMPC_z17[4];
output->z19[0] = smoothMPC_z18[0];
output->z19[1] = smoothMPC_z18[1];
output->z19[2] = smoothMPC_z18[2];
output->z19[3] = smoothMPC_z18[3];
output->z19[4] = smoothMPC_z18[4];
output->z20[0] = smoothMPC_z19[0];
output->z20[1] = smoothMPC_z19[1];
output->z20[2] = smoothMPC_z19[2];
output->z20[3] = smoothMPC_z19[3];
output->z20[4] = smoothMPC_z19[4];
output->z21[0] = smoothMPC_z20[0];
output->z21[1] = smoothMPC_z20[1];
output->z21[2] = smoothMPC_z20[2];
output->z21[3] = smoothMPC_z20[3];
output->z21[4] = smoothMPC_z20[4];
output->z22[0] = smoothMPC_z21[0];
output->z22[1] = smoothMPC_z21[1];
output->z22[2] = smoothMPC_z21[2];
output->z22[3] = smoothMPC_z21[3];
output->z22[4] = smoothMPC_z21[4];
output->z23[0] = smoothMPC_z22[0];
output->z23[1] = smoothMPC_z22[1];
output->z23[2] = smoothMPC_z22[2];
output->z23[3] = smoothMPC_z22[3];
output->z23[4] = smoothMPC_z22[4];
output->z24[0] = smoothMPC_z23[0];
output->z24[1] = smoothMPC_z23[1];
output->z24[2] = smoothMPC_z23[2];
output->z24[3] = smoothMPC_z23[3];
output->z24[4] = smoothMPC_z23[4];
output->z25[0] = smoothMPC_z24[0];
output->z25[1] = smoothMPC_z24[1];
output->z25[2] = smoothMPC_z24[2];
output->z25[3] = smoothMPC_z24[3];
output->z25[4] = smoothMPC_z24[4];
output->z26[0] = smoothMPC_z25[0];
output->z26[1] = smoothMPC_z25[1];
output->z26[2] = smoothMPC_z25[2];
output->z26[3] = smoothMPC_z25[3];
output->z26[4] = smoothMPC_z25[4];
output->z27[0] = smoothMPC_z26[0];
output->z27[1] = smoothMPC_z26[1];
output->z27[2] = smoothMPC_z26[2];
output->z27[3] = smoothMPC_z26[3];
output->z27[4] = smoothMPC_z26[4];
output->z28[0] = smoothMPC_z27[0];
output->z28[1] = smoothMPC_z27[1];
output->z28[2] = smoothMPC_z27[2];
output->z28[3] = smoothMPC_z27[3];
output->z28[4] = smoothMPC_z27[4];
output->z29[0] = smoothMPC_z28[0];
output->z29[1] = smoothMPC_z28[1];
output->z29[2] = smoothMPC_z28[2];
output->z29[3] = smoothMPC_z28[3];
output->z29[4] = smoothMPC_z28[4];
output->z30[0] = smoothMPC_z29[0];
output->z30[1] = smoothMPC_z29[1];
output->z30[2] = smoothMPC_z29[2];

#if smoothMPC_SET_TIMING == 1
info->solvetime = smoothMPC_toc(&solvertimer);
#if smoothMPC_SET_PRINTLEVEL > 0 && smoothMPC_SET_TIMING == 1
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
