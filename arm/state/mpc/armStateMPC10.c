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

#include "armStateMPC.h"

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
 * Initializes a vector of length 114 with a value.
 */
void armStateMPC_LA_INITIALIZEVECTOR_114(armStateMPC_FLOAT* vec, armStateMPC_FLOAT value)
{
	int i;
	for( i=0; i<114; i++ )
	{
		vec[i] = value;
	}
}


/*
 * Initializes a vector of length 60 with a value.
 */
void armStateMPC_LA_INITIALIZEVECTOR_60(armStateMPC_FLOAT* vec, armStateMPC_FLOAT value)
{
	int i;
	for( i=0; i<60; i++ )
	{
		vec[i] = value;
	}
}


/*
 * Initializes a vector of length 228 with a value.
 */
void armStateMPC_LA_INITIALIZEVECTOR_228(armStateMPC_FLOAT* vec, armStateMPC_FLOAT value)
{
	int i;
	for( i=0; i<228; i++ )
	{
		vec[i] = value;
	}
}


/* 
 * Calculates a dot product and adds it to a variable: z += x'*y; 
 * This function is for vectors of length 228.
 */
void armStateMPC_LA_DOTACC_228(armStateMPC_FLOAT *x, armStateMPC_FLOAT *y, armStateMPC_FLOAT *z)
{
	int i;
	for( i=0; i<228; i++ ){
		*z += x[i]*y[i];
	}
}


/*
 * Calculates the gradient and the value for a quadratic function 0.5*z'*H*z + f'*z
 *
 * INPUTS:     H  - Symmetric Hessian, diag matrix of size [12 x 12]
 *             f  - column vector of size 12
 *             z  - column vector of size 12
 *
 * OUTPUTS: grad  - gradient at z (= H*z + f), column vector of size 12
 *          value <-- value + 0.5*z'*H*z + f'*z (value will be modified)
 */
void armStateMPC_LA_DIAG_QUADFCN_12(armStateMPC_FLOAT* H, armStateMPC_FLOAT* f, armStateMPC_FLOAT* z, armStateMPC_FLOAT* grad, armStateMPC_FLOAT* value)
{
	int i;
	armStateMPC_FLOAT hz;	
	for( i=0; i<12; i++){
		hz = H[i]*z[i];
		grad[i] = hz + f[i];
		*value += 0.5*hz*z[i] + f[i]*z[i];
	}
}


/*
 * Calculates the gradient and the value for a quadratic function 0.5*z'*H*z + f'*z
 *
 * INPUTS:     H  - Symmetric Hessian, diag matrix of size [6 x 6]
 *             f  - column vector of size 6
 *             z  - column vector of size 6
 *
 * OUTPUTS: grad  - gradient at z (= H*z + f), column vector of size 6
 *          value <-- value + 0.5*z'*H*z + f'*z (value will be modified)
 */
void armStateMPC_LA_DIAG_QUADFCN_6(armStateMPC_FLOAT* H, armStateMPC_FLOAT* f, armStateMPC_FLOAT* z, armStateMPC_FLOAT* grad, armStateMPC_FLOAT* value)
{
	int i;
	armStateMPC_FLOAT hz;	
	for( i=0; i<6; i++){
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
void armStateMPC_LA_DIAGZERO_MVMSUB6_6(armStateMPC_FLOAT *B, armStateMPC_FLOAT *u, armStateMPC_FLOAT *b, armStateMPC_FLOAT *l, armStateMPC_FLOAT *r, armStateMPC_FLOAT *z, armStateMPC_FLOAT *y)
{
	int i;
	armStateMPC_FLOAT Bu[6];
	armStateMPC_FLOAT norm = *y;
	armStateMPC_FLOAT lr = 0;

	/* do A*x + B*u first */
	for( i=0; i<6; i++ ){
		Bu[i] = B[i]*u[i];
	}	

	for( i=0; i<6; i++ ){
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
void armStateMPC_LA_DENSE_DIAGZERO_MVMSUB3_6_12_12(armStateMPC_FLOAT *A, armStateMPC_FLOAT *x, armStateMPC_FLOAT *B, armStateMPC_FLOAT *u, armStateMPC_FLOAT *b, armStateMPC_FLOAT *l, armStateMPC_FLOAT *r, armStateMPC_FLOAT *z, armStateMPC_FLOAT *y)
{
	int i;
	int j;
	int k = 0;
	armStateMPC_FLOAT AxBu[6];
	armStateMPC_FLOAT norm = *y;
	armStateMPC_FLOAT lr = 0;

	/* do A*x + B*u first */
	for( i=0; i<6; i++ ){
		AxBu[i] = A[k++]*x[0] + B[i]*u[i];
	}	

	for( j=1; j<12; j++ ){		
		for( i=0; i<6; i++ ){
			AxBu[i] += A[k++]*x[j];
		}
	}

	for( i=0; i<6; i++ ){
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
void armStateMPC_LA_DENSE_DIAGZERO_MVMSUB3_6_12_6(armStateMPC_FLOAT *A, armStateMPC_FLOAT *x, armStateMPC_FLOAT *B, armStateMPC_FLOAT *u, armStateMPC_FLOAT *b, armStateMPC_FLOAT *l, armStateMPC_FLOAT *r, armStateMPC_FLOAT *z, armStateMPC_FLOAT *y)
{
	int i;
	int j;
	int k = 0;
	armStateMPC_FLOAT AxBu[6];
	armStateMPC_FLOAT norm = *y;
	armStateMPC_FLOAT lr = 0;

	/* do A*x + B*u first */
	for( i=0; i<6; i++ ){
		AxBu[i] = A[k++]*x[0] + B[i]*u[i];
	}	

	for( j=1; j<12; j++ ){		
		for( i=0; i<6; i++ ){
			AxBu[i] += A[k++]*x[j];
		}
	}

	for( i=0; i<6; i++ ){
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
 * where A is of size [6 x 12] and stored in column major format.
 * and B is of size [6 x 12] and stored in diagzero format
 * Note the transposes of A and B!
 */
void armStateMPC_LA_DENSE_DIAGZERO_MTVM2_6_12_6(armStateMPC_FLOAT *A, armStateMPC_FLOAT *x, armStateMPC_FLOAT *B, armStateMPC_FLOAT *y, armStateMPC_FLOAT *z)
{
	int i;
	int j;
	int k = 0;
	for( i=0; i<6; i++ ){
		z[i] = 0;
		for( j=0; j<6; j++ ){
			z[i] += A[k++]*x[j];
		}
		z[i] += B[i]*y[i];
	}
	for( i=6 ;i<12; i++ ){
		z[i] = 0;
		for( j=0; j<6; j++ ){
			z[i] += A[k++]*x[j];
		}
	}
}


/*
 * Matrix vector multiplication y = M'*x where M is of size [6 x 6]
 * and stored in diagzero format. Note the transpose of M!
 */
void armStateMPC_LA_DIAGZERO_MTVM_6_6(armStateMPC_FLOAT *M, armStateMPC_FLOAT *x, armStateMPC_FLOAT *y)
{
	int i;
	for( i=0; i<6; i++ ){
		y[i] = M[i]*x[i];
	}
}


/*
 * Vector subtraction and addition.
 *	 Input: five vectors t, tidx, u, v, w and two scalars z and r
 *	 Output: y = t(tidx) - u + w
 *           z = z - v'*x;
 *           r = max([norm(y,inf), z]);
 * for vectors of length 12. Output z is of course scalar.
 */
void armStateMPC_LA_VSUBADD3_12(armStateMPC_FLOAT* t, armStateMPC_FLOAT* u, int* uidx, armStateMPC_FLOAT* v, armStateMPC_FLOAT* w, armStateMPC_FLOAT* y, armStateMPC_FLOAT* z, armStateMPC_FLOAT* r)
{
	int i;
	armStateMPC_FLOAT norm = *r;
	armStateMPC_FLOAT vx = 0;
	armStateMPC_FLOAT x;
	for( i=0; i<12; i++){
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
 * for vectors of length 12. Output z is of course scalar.
 */
void armStateMPC_LA_VSUBADD2_12(armStateMPC_FLOAT* t, int* tidx, armStateMPC_FLOAT* u, armStateMPC_FLOAT* v, armStateMPC_FLOAT* w, armStateMPC_FLOAT* y, armStateMPC_FLOAT* z, armStateMPC_FLOAT* r)
{
	int i;
	armStateMPC_FLOAT norm = *r;
	armStateMPC_FLOAT vx = 0;
	armStateMPC_FLOAT x;
	for( i=0; i<12; i++){
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
 * for vectors of length 6. Output z is of course scalar.
 */
void armStateMPC_LA_VSUBADD3_6(armStateMPC_FLOAT* t, armStateMPC_FLOAT* u, int* uidx, armStateMPC_FLOAT* v, armStateMPC_FLOAT* w, armStateMPC_FLOAT* y, armStateMPC_FLOAT* z, armStateMPC_FLOAT* r)
{
	int i;
	armStateMPC_FLOAT norm = *r;
	armStateMPC_FLOAT vx = 0;
	armStateMPC_FLOAT x;
	for( i=0; i<6; i++){
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
 * for vectors of length 6. Output z is of course scalar.
 */
void armStateMPC_LA_VSUBADD2_6(armStateMPC_FLOAT* t, int* tidx, armStateMPC_FLOAT* u, armStateMPC_FLOAT* v, armStateMPC_FLOAT* w, armStateMPC_FLOAT* y, armStateMPC_FLOAT* z, armStateMPC_FLOAT* r)
{
	int i;
	armStateMPC_FLOAT norm = *r;
	armStateMPC_FLOAT vx = 0;
	armStateMPC_FLOAT x;
	for( i=0; i<6; i++){
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
 * Special function for box constraints of length 12
 * Returns also L/S, a value that is often used elsewhere.
 */
void armStateMPC_LA_INEQ_B_GRAD_12_12_12(armStateMPC_FLOAT *lu, armStateMPC_FLOAT *su, armStateMPC_FLOAT *ru, armStateMPC_FLOAT *ll, armStateMPC_FLOAT *sl, armStateMPC_FLOAT *rl, int* lbIdx, int* ubIdx, armStateMPC_FLOAT *grad, armStateMPC_FLOAT *lubysu, armStateMPC_FLOAT *llbysl)
{
	int i;
	for( i=0; i<12; i++ ){
		grad[i] = 0;
	}
	for( i=0; i<12; i++ ){		
		llbysl[i] = ll[i] / sl[i];
		grad[lbIdx[i]] -= llbysl[i]*rl[i];
	}
	for( i=0; i<12; i++ ){
		lubysu[i] = lu[i] / su[i];
		grad[ubIdx[i]] += lubysu[i]*ru[i];
	}
}


/*
 * Computes inequality constraints gradient-
 * Special function for box constraints of length 6
 * Returns also L/S, a value that is often used elsewhere.
 */
void armStateMPC_LA_INEQ_B_GRAD_6_6_6(armStateMPC_FLOAT *lu, armStateMPC_FLOAT *su, armStateMPC_FLOAT *ru, armStateMPC_FLOAT *ll, armStateMPC_FLOAT *sl, armStateMPC_FLOAT *rl, int* lbIdx, int* ubIdx, armStateMPC_FLOAT *grad, armStateMPC_FLOAT *lubysu, armStateMPC_FLOAT *llbysl)
{
	int i;
	for( i=0; i<6; i++ ){
		grad[i] = 0;
	}
	for( i=0; i<6; i++ ){		
		llbysl[i] = ll[i] / sl[i];
		grad[lbIdx[i]] -= llbysl[i]*rl[i];
	}
	for( i=0; i<6; i++ ){
		lubysu[i] = lu[i] / su[i];
		grad[ubIdx[i]] += lubysu[i]*ru[i];
	}
}


/*
 * Addition of three vectors  z = u + w + v
 * of length 114.
 */
void armStateMPC_LA_VVADD3_114(armStateMPC_FLOAT *u, armStateMPC_FLOAT *v, armStateMPC_FLOAT *w, armStateMPC_FLOAT *z)
{
	int i;
	for( i=0; i<114; i++ ){
		z[i] = u[i] + v[i] + w[i];
	}
}


/*
 * Special function to compute the diagonal cholesky factorization of the 
 * positive definite augmented Hessian for block size 12.
 *
 * Inputs: - H = diagonal cost Hessian in diagonal storage format
 *         - llbysl = L / S of lower bounds
 *         - lubysu = L / S of upper bounds
 *
 * Output: Phi = sqrt(H + diag(llbysl) + diag(lubysu))
 * where Phi is stored in diagonal storage format
 */
void armStateMPC_LA_DIAG_CHOL_ONELOOP_LBUB_12_12_12(armStateMPC_FLOAT *H, armStateMPC_FLOAT *llbysl, int* lbIdx, armStateMPC_FLOAT *lubysu, int* ubIdx, armStateMPC_FLOAT *Phi)


{
	int i;
	
	/* compute cholesky */
	for( i=0; i<12; i++ ){
		Phi[i] = H[i] + llbysl[i] + lubysu[i];

#if armStateMPC_SET_PRINTLEVEL > 0 && defined PRINTNUMERICALWARNINGS
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
 * where A is to be computed and is of size [6 x 12],
 * B is given and of size [6 x 12], L is a diagonal
 * matrix of size 6 stored in diagonal matrix 
 * storage format. Note the transpose of L has no impact!
 *
 * Result: A in column major storage format.
 *
 */
void armStateMPC_LA_DIAG_MATRIXFORWARDSUB_6_12(armStateMPC_FLOAT *L, armStateMPC_FLOAT *B, armStateMPC_FLOAT *A)
{
    int i,j;
	 int k = 0;

	for( j=0; j<12; j++){
		for( i=0; i<6; i++){
			A[k] = B[k]/L[j];
			k++;
		}
	}

}


/**
 * Forward substitution for the matrix equation A*L' = B
 * where A is to be computed and is of size [6 x 12],
 * B is given and of size [6 x 12], L is a diagonal
 *  matrix of size 12 stored in diagonal 
 * storage format. Note the transpose of L!
 *
 * Result: A in diagonalzero storage format.
 *
 */
void armStateMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_6_12(armStateMPC_FLOAT *L, armStateMPC_FLOAT *B, armStateMPC_FLOAT *A)
{
	int j;
    for( j=0; j<12; j++ ){   
		A[j] = B[j]/L[j];
     }
}


/**
 * Compute C = A*B' where 
 *
 *	size(A) = [6 x 12]
 *  size(B) = [6 x 12] in diagzero format
 * 
 * A and C matrices are stored in column major format.
 * 
 * 
 */
void armStateMPC_LA_DENSE_DIAGZERO_MMTM_6_12_6(armStateMPC_FLOAT *A, armStateMPC_FLOAT *B, armStateMPC_FLOAT *C)
{
    int i, j;
	
	for( i=0; i<6; i++ ){
		for( j=0; j<6; j++){
			C[j*6+i] = B[i*6+j]*A[i];
		}
	}

}


/**
 * Forward substitution to solve L*y = b where L is a
 * diagonal matrix in vector storage format.
 * 
 * The dimensions involved are 12.
 */
void armStateMPC_LA_DIAG_FORWARDSUB_12(armStateMPC_FLOAT *L, armStateMPC_FLOAT *b, armStateMPC_FLOAT *y)
{
    int i;

    for( i=0; i<12; i++ ){
		y[i] = b[i]/L[i];
    }
}


/*
 * Special function to compute the diagonal cholesky factorization of the 
 * positive definite augmented Hessian for block size 6.
 *
 * Inputs: - H = diagonal cost Hessian in diagonal storage format
 *         - llbysl = L / S of lower bounds
 *         - lubysu = L / S of upper bounds
 *
 * Output: Phi = sqrt(H + diag(llbysl) + diag(lubysu))
 * where Phi is stored in diagonal storage format
 */
void armStateMPC_LA_DIAG_CHOL_ONELOOP_LBUB_6_6_6(armStateMPC_FLOAT *H, armStateMPC_FLOAT *llbysl, int* lbIdx, armStateMPC_FLOAT *lubysu, int* ubIdx, armStateMPC_FLOAT *Phi)


{
	int i;
	
	/* compute cholesky */
	for( i=0; i<6; i++ ){
		Phi[i] = H[i] + llbysl[i] + lubysu[i];

#if armStateMPC_SET_PRINTLEVEL > 0 && defined PRINTNUMERICALWARNINGS
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
 * where A is to be computed and is of size [6 x 6],
 * B is given and of size [6 x 6], L is a diagonal
 *  matrix of size 6 stored in diagonal 
 * storage format. Note the transpose of L!
 *
 * Result: A in diagonalzero storage format.
 *
 */
void armStateMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_6_6(armStateMPC_FLOAT *L, armStateMPC_FLOAT *B, armStateMPC_FLOAT *A)
{
	int j;
    for( j=0; j<6; j++ ){   
		A[j] = B[j]/L[j];
     }
}


/**
 * Forward substitution to solve L*y = b where L is a
 * diagonal matrix in vector storage format.
 * 
 * The dimensions involved are 6.
 */
void armStateMPC_LA_DIAG_FORWARDSUB_6(armStateMPC_FLOAT *L, armStateMPC_FLOAT *b, armStateMPC_FLOAT *y)
{
    int i;

    for( i=0; i<6; i++ ){
		y[i] = b[i]/L[i];
    }
}


/**
 * Compute L = B*B', where L is lower triangular of size NXp1
 * and B is a diagzero matrix of size [6 x 12] in column
 * storage format.
 * 
 */
void armStateMPC_LA_DIAGZERO_MMT_6(armStateMPC_FLOAT *B, armStateMPC_FLOAT *L)
{
    int i, ii, di;
    
    ii = 0; di = 0;
    for( i=0; i<6; i++ ){        
		L[ii+i] = B[i]*B[i];
        ii += ++di;
    }
}


/* 
 * Computes r = b - B*u
 * B is stored in diagzero format
 */
void armStateMPC_LA_DIAGZERO_MVMSUB7_6(armStateMPC_FLOAT *B, armStateMPC_FLOAT *u, armStateMPC_FLOAT *b, armStateMPC_FLOAT *r)
{
	int i;

	for( i=0; i<6; i++ ){
		r[i] = b[i] - B[i]*u[i];
	}	
	
}


/**
 * Compute L = A*A' + B*B', where L is lower triangular of size NXp1
 * and A is a dense matrix of size [6 x 12] in column
 * storage format, and B is of size [6 x 12] diagonalzero
 * storage format.
 * 
 * THIS ONE HAS THE WORST ACCES PATTERN POSSIBLE. 
 * POSSIBKE FIX: PUT A AND B INTO ROW MAJOR FORMAT FIRST.
 * 
 */
void armStateMPC_LA_DENSE_DIAGZERO_MMT2_6_12_12(armStateMPC_FLOAT *A, armStateMPC_FLOAT *B, armStateMPC_FLOAT *L)
{
    int i, j, k, ii, di;
    armStateMPC_FLOAT ltemp;
    
    ii = 0; di = 0;
    for( i=0; i<6; i++ ){        
        for( j=0; j<=i; j++ ){
            ltemp = 0; 
            for( k=0; k<12; k++ ){
                ltemp += A[k*6+i]*A[k*6+j];
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
void armStateMPC_LA_DENSE_DIAGZERO_2MVMSUB2_6_12_12(armStateMPC_FLOAT *A, armStateMPC_FLOAT *x, armStateMPC_FLOAT *B, armStateMPC_FLOAT *u, armStateMPC_FLOAT *b, armStateMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;

	for( i=0; i<6; i++ ){
		r[i] = b[i] - A[k++]*x[0] - B[i]*u[i];
	}	

	for( j=1; j<12; j++ ){		
		for( i=0; i<6; i++ ){
			r[i] -= A[k++]*x[j];
		}
	}
	
}


/**
 * Compute L = A*A' + B*B', where L is lower triangular of size NXp1
 * and A is a dense matrix of size [6 x 12] in column
 * storage format, and B is of size [6 x 6] diagonalzero
 * storage format.
 * 
 * THIS ONE HAS THE WORST ACCES PATTERN POSSIBLE. 
 * POSSIBKE FIX: PUT A AND B INTO ROW MAJOR FORMAT FIRST.
 * 
 */
void armStateMPC_LA_DENSE_DIAGZERO_MMT2_6_12_6(armStateMPC_FLOAT *A, armStateMPC_FLOAT *B, armStateMPC_FLOAT *L)
{
    int i, j, k, ii, di;
    armStateMPC_FLOAT ltemp;
    
    ii = 0; di = 0;
    for( i=0; i<6; i++ ){        
        for( j=0; j<=i; j++ ){
            ltemp = 0; 
            for( k=0; k<12; k++ ){
                ltemp += A[k*6+i]*A[k*6+j];
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
void armStateMPC_LA_DENSE_DIAGZERO_2MVMSUB2_6_12_6(armStateMPC_FLOAT *A, armStateMPC_FLOAT *x, armStateMPC_FLOAT *B, armStateMPC_FLOAT *u, armStateMPC_FLOAT *b, armStateMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;

	for( i=0; i<6; i++ ){
		r[i] = b[i] - A[k++]*x[0] - B[i]*u[i];
	}	

	for( j=1; j<12; j++ ){		
		for( i=0; i<6; i++ ){
			r[i] -= A[k++]*x[j];
		}
	}
	
}


/**
 * Cholesky factorization as above, but working on a matrix in 
 * lower triangular storage format of size 6 and outputting
 * the Cholesky factor to matrix L in lower triangular format.
 */
void armStateMPC_LA_DENSE_CHOL_6(armStateMPC_FLOAT *A, armStateMPC_FLOAT *L)
{
    int i, j, k, di, dj;
	 int ii, jj;

    armStateMPC_FLOAT l;
    armStateMPC_FLOAT Mii;

	/* copy A to L first and then operate on L */
	/* COULD BE OPTIMIZED */
	ii=0; di=0;
	for( i=0; i<6; i++ ){
		for( j=0; j<=i; j++ ){
			L[ii+j] = A[ii+j];
		}
		ii += ++di;
	}    
	
	/* factor L */
	ii=0; di=0;
    for( i=0; i<6; i++ ){
        l = 0;
        for( k=0; k<i; k++ ){
            l += L[ii+k]*L[ii+k];
        }        
        
        Mii = L[ii+i] - l;
        
#if armStateMPC_SET_PRINTLEVEL > 0 && defined PRINTNUMERICALWARNINGS
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
        for( j=i+1; j<6; j++ ){
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
 * The dimensions involved are 6.
 */
void armStateMPC_LA_DENSE_FORWARDSUB_6(armStateMPC_FLOAT *L, armStateMPC_FLOAT *b, armStateMPC_FLOAT *y)
{
    int i,j,ii,di;
    armStateMPC_FLOAT yel;
            
    ii = 0; di = 0;
    for( i=0; i<6; i++ ){
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
 * where A is to be computed and is of size [6 x 6],
 * B is given and of size [6 x 6], L is a lower tri-
 * angular matrix of size 6 stored in lower triangular 
 * storage format. Note the transpose of L AND B!
 *
 * Result: A in column major storage format.
 *
 */
void armStateMPC_LA_DENSE_MATRIXTFORWARDSUB_6_6(armStateMPC_FLOAT *L, armStateMPC_FLOAT *B, armStateMPC_FLOAT *A)
{
    int i,j,k,ii,di;
    armStateMPC_FLOAT a;
    
    ii=0; di=0;
    for( j=0; j<6; j++ ){        
        for( i=0; i<6; i++ ){
            a = B[i*6+j];
            for( k=0; k<j; k++ ){
                a -= A[k*6+i]*L[ii+k];
            }    

			/* saturate for numerical stability */
			a = MIN(a, BIGM);
			a = MAX(a, -BIGM); 

			A[j*6+i] = a/L[ii+j];			
        }
        ii += ++di;
    }
}


/**
 * Compute L = L - A*A', where L is lower triangular of size 6
 * and A is a dense matrix of size [6 x 6] in column
 * storage format.
 * 
 * THIS ONE HAS THE WORST ACCES PATTERN POSSIBLE. 
 * POSSIBKE FIX: PUT A INTO ROW MAJOR FORMAT FIRST.
 * 
 */
void armStateMPC_LA_DENSE_MMTSUB_6_6(armStateMPC_FLOAT *A, armStateMPC_FLOAT *L)
{
    int i, j, k, ii, di;
    armStateMPC_FLOAT ltemp;
    
    ii = 0; di = 0;
    for( i=0; i<6; i++ ){        
        for( j=0; j<=i; j++ ){
            ltemp = 0; 
            for( k=0; k<6; k++ ){
                ltemp += A[k*6+i]*A[k*6+j];
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
void armStateMPC_LA_DENSE_MVMSUB1_6_6(armStateMPC_FLOAT *A, armStateMPC_FLOAT *x, armStateMPC_FLOAT *b, armStateMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;

	for( i=0; i<6; i++ ){
		r[i] = b[i] - A[k++]*x[0];
	}	
	for( j=1; j<6; j++ ){		
		for( i=0; i<6; i++ ){
			r[i] -= A[k++]*x[j];
		}
	}
}


/**
 * Backward Substitution to solve L^T*x = y where L is a
 * lower triangular matrix in triangular storage format.
 * 
 * All involved dimensions are 6.
 */
void armStateMPC_LA_DENSE_BACKWARDSUB_6(armStateMPC_FLOAT *L, armStateMPC_FLOAT *y, armStateMPC_FLOAT *x)
{
    int i, ii, di, j, jj, dj;
    armStateMPC_FLOAT xel;    
	int start = 15;
    
    /* now solve L^T*x = y by backward substitution */
    ii = start; di = 5;
    for( i=5; i>=0; i-- ){        
        xel = y[i];        
        jj = start; dj = 5;
        for( j=5; j>i; j-- ){
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
 * Matrix vector multiplication y = b - M'*x where M is of size [6 x 6]
 * and stored in column major format. Note the transpose of M!
 */
void armStateMPC_LA_DENSE_MTVMSUB_6_6(armStateMPC_FLOAT *A, armStateMPC_FLOAT *x, armStateMPC_FLOAT *b, armStateMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0; 
	for( i=0; i<6; i++ ){
		r[i] = b[i];
		for( j=0; j<6; j++ ){
			r[i] -= A[k++]*x[j];
		}
	}
}


/*
 * Vector subtraction z = -x - y for vectors of length 114.
 */
void armStateMPC_LA_VSUB2_114(armStateMPC_FLOAT *x, armStateMPC_FLOAT *y, armStateMPC_FLOAT *z)
{
	int i;
	for( i=0; i<114; i++){
		z[i] = -x[i] - y[i];
	}
}


/**
 * Forward-Backward-Substitution to solve L*L^T*x = b where L is a
 * diagonal matrix of size 12 in vector
 * storage format.
 */
void armStateMPC_LA_DIAG_FORWARDBACKWARDSUB_12(armStateMPC_FLOAT *L, armStateMPC_FLOAT *b, armStateMPC_FLOAT *x)
{
    int i;
            
    /* solve Ly = b by forward and backward substitution */
    for( i=0; i<12; i++ ){
		x[i] = b[i]/(L[i]*L[i]);
    }
    
}


/**
 * Forward-Backward-Substitution to solve L*L^T*x = b where L is a
 * diagonal matrix of size 6 in vector
 * storage format.
 */
void armStateMPC_LA_DIAG_FORWARDBACKWARDSUB_6(armStateMPC_FLOAT *L, armStateMPC_FLOAT *b, armStateMPC_FLOAT *x)
{
    int i;
            
    /* solve Ly = b by forward and backward substitution */
    for( i=0; i<6; i++ ){
		x[i] = b[i]/(L[i]*L[i]);
    }
    
}


/*
 * Vector subtraction z = x(xidx) - y where y, z and xidx are of length 12,
 * and x has length 12 and is indexed through yidx.
 */
void armStateMPC_LA_VSUB_INDEXED_12(armStateMPC_FLOAT *x, int* xidx, armStateMPC_FLOAT *y, armStateMPC_FLOAT *z)
{
	int i;
	for( i=0; i<12; i++){
		z[i] = x[xidx[i]] - y[i];
	}
}


/*
 * Vector subtraction x = -u.*v - w for vectors of length 12.
 */
void armStateMPC_LA_VSUB3_12(armStateMPC_FLOAT *u, armStateMPC_FLOAT *v, armStateMPC_FLOAT *w, armStateMPC_FLOAT *x)
{
	int i;
	for( i=0; i<12; i++){
		x[i] = -u[i]*v[i] - w[i];
	}
}


/*
 * Vector subtraction z = -x - y(yidx) where y is of length 12
 * and z, x and yidx are of length 12.
 */
void armStateMPC_LA_VSUB2_INDEXED_12(armStateMPC_FLOAT *x, armStateMPC_FLOAT *y, int* yidx, armStateMPC_FLOAT *z)
{
	int i;
	for( i=0; i<12; i++){
		z[i] = -x[i] - y[yidx[i]];
	}
}


/*
 * Vector subtraction z = x(xidx) - y where y, z and xidx are of length 6,
 * and x has length 6 and is indexed through yidx.
 */
void armStateMPC_LA_VSUB_INDEXED_6(armStateMPC_FLOAT *x, int* xidx, armStateMPC_FLOAT *y, armStateMPC_FLOAT *z)
{
	int i;
	for( i=0; i<6; i++){
		z[i] = x[xidx[i]] - y[i];
	}
}


/*
 * Vector subtraction x = -u.*v - w for vectors of length 6.
 */
void armStateMPC_LA_VSUB3_6(armStateMPC_FLOAT *u, armStateMPC_FLOAT *v, armStateMPC_FLOAT *w, armStateMPC_FLOAT *x)
{
	int i;
	for( i=0; i<6; i++){
		x[i] = -u[i]*v[i] - w[i];
	}
}


/*
 * Vector subtraction z = -x - y(yidx) where y is of length 6
 * and z, x and yidx are of length 6.
 */
void armStateMPC_LA_VSUB2_INDEXED_6(armStateMPC_FLOAT *x, armStateMPC_FLOAT *y, int* yidx, armStateMPC_FLOAT *z)
{
	int i;
	for( i=0; i<6; i++){
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
 * armStateMPC_NOPROGRESS (should be negative).
 */
int armStateMPC_LINESEARCH_BACKTRACKING_AFFINE(armStateMPC_FLOAT *l, armStateMPC_FLOAT *s, armStateMPC_FLOAT *dl, armStateMPC_FLOAT *ds, armStateMPC_FLOAT *a, armStateMPC_FLOAT *mu_aff)
{
    int i;
	int lsIt=1;    
    armStateMPC_FLOAT dltemp;
    armStateMPC_FLOAT dstemp;
    armStateMPC_FLOAT mya = 1.0;
    armStateMPC_FLOAT mymu;
        
    while( 1 ){                        

        /* 
         * Compute both snew and wnew together.
         * We compute also mu_affine along the way here, as the
         * values might be in registers, so it should be cheaper.
         */
        mymu = 0;
        for( i=0; i<228; i++ ){
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
        if( i == 228 ){
            break;
        } else {
            mya *= armStateMPC_SET_LS_SCALE_AFF;
            if( mya < armStateMPC_SET_LS_MINSTEP ){
                return armStateMPC_NOPROGRESS;
            }
        }
    }
    
    /* return new values and iteration counter */
    *a = mya;
    *mu_aff = mymu / (armStateMPC_FLOAT)228;
    return lsIt;
}


/*
 * Vector subtraction x = u.*v - a where a is a scalar
*  and x,u,v are vectors of length 228.
 */
void armStateMPC_LA_VSUB5_228(armStateMPC_FLOAT *u, armStateMPC_FLOAT *v, armStateMPC_FLOAT a, armStateMPC_FLOAT *x)
{
	int i;
	for( i=0; i<228; i++){
		x[i] = u[i]*v[i] - a;
	}
}


/*
 * Computes x=0; x(uidx) += u/su; x(vidx) -= v/sv where x is of length 12,
 * u, su, uidx are of length 12 and v, sv, vidx are of length 12.
 */
void armStateMPC_LA_VSUB6_INDEXED_12_12_12(armStateMPC_FLOAT *u, armStateMPC_FLOAT *su, int* uidx, armStateMPC_FLOAT *v, armStateMPC_FLOAT *sv, int* vidx, armStateMPC_FLOAT *x)
{
	int i;
	for( i=0; i<12; i++ ){
		x[i] = 0;
	}
	for( i=0; i<12; i++){
		x[uidx[i]] += u[i]/su[i];
	}
	for( i=0; i<12; i++){
		x[vidx[i]] -= v[i]/sv[i];
	}
}


/* 
 * Computes r =  B*u
 * where B is stored in diagzero format
 */
void armStateMPC_LA_DIAGZERO_MVM_6(armStateMPC_FLOAT *B, armStateMPC_FLOAT *u, armStateMPC_FLOAT *r)
{
	int i;

	for( i=0; i<6; i++ ){
		r[i] = B[i]*u[i];
	}	
	
}


/* 
 * Computes r = A*x + B*u
 * where A is stored in column major format
 * and B is stored in diagzero format
 */
void armStateMPC_LA_DENSE_DIAGZERO_2MVMADD_6_12_12(armStateMPC_FLOAT *A, armStateMPC_FLOAT *x, armStateMPC_FLOAT *B, armStateMPC_FLOAT *u, armStateMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;

	for( i=0; i<6; i++ ){
		r[i] = A[k++]*x[0] + B[i]*u[i];
	}	

	for( j=1; j<12; j++ ){		
		for( i=0; i<6; i++ ){
			r[i] += A[k++]*x[j];
		}
	}
	
}


/*
 * Computes x=0; x(uidx) += u/su; x(vidx) -= v/sv where x is of length 6,
 * u, su, uidx are of length 6 and v, sv, vidx are of length 6.
 */
void armStateMPC_LA_VSUB6_INDEXED_6_6_6(armStateMPC_FLOAT *u, armStateMPC_FLOAT *su, int* uidx, armStateMPC_FLOAT *v, armStateMPC_FLOAT *sv, int* vidx, armStateMPC_FLOAT *x)
{
	int i;
	for( i=0; i<6; i++ ){
		x[i] = 0;
	}
	for( i=0; i<6; i++){
		x[uidx[i]] += u[i]/su[i];
	}
	for( i=0; i<6; i++){
		x[vidx[i]] -= v[i]/sv[i];
	}
}


/* 
 * Computes r = A*x + B*u
 * where A is stored in column major format
 * and B is stored in diagzero format
 */
void armStateMPC_LA_DENSE_DIAGZERO_2MVMADD_6_12_6(armStateMPC_FLOAT *A, armStateMPC_FLOAT *x, armStateMPC_FLOAT *B, armStateMPC_FLOAT *u, armStateMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;

	for( i=0; i<6; i++ ){
		r[i] = A[k++]*x[0] + B[i]*u[i];
	}	

	for( j=1; j<12; j++ ){		
		for( i=0; i<6; i++ ){
			r[i] += A[k++]*x[j];
		}
	}
	
}


/*
 * Vector subtraction z = x - y for vectors of length 114.
 */
void armStateMPC_LA_VSUB_114(armStateMPC_FLOAT *x, armStateMPC_FLOAT *y, armStateMPC_FLOAT *z)
{
	int i;
	for( i=0; i<114; i++){
		z[i] = x[i] - y[i];
	}
}


/** 
 * Computes z = -r./s - u.*y(y)
 * where all vectors except of y are of length 12 (length of y >= 12).
 */
void armStateMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_12(armStateMPC_FLOAT *r, armStateMPC_FLOAT *s, armStateMPC_FLOAT *u, armStateMPC_FLOAT *y, int* yidx, armStateMPC_FLOAT *z)
{
	int i;
	for( i=0; i<12; i++ ){
		z[i] = -r[i]/s[i] - u[i]*y[yidx[i]];
	}
}


/** 
 * Computes z = -r./s + u.*y(y)
 * where all vectors except of y are of length 12 (length of y >= 12).
 */
void armStateMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_12(armStateMPC_FLOAT *r, armStateMPC_FLOAT *s, armStateMPC_FLOAT *u, armStateMPC_FLOAT *y, int* yidx, armStateMPC_FLOAT *z)
{
	int i;
	for( i=0; i<12; i++ ){
		z[i] = -r[i]/s[i] + u[i]*y[yidx[i]];
	}
}


/** 
 * Computes z = -r./s - u.*y(y)
 * where all vectors except of y are of length 6 (length of y >= 6).
 */
void armStateMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_6(armStateMPC_FLOAT *r, armStateMPC_FLOAT *s, armStateMPC_FLOAT *u, armStateMPC_FLOAT *y, int* yidx, armStateMPC_FLOAT *z)
{
	int i;
	for( i=0; i<6; i++ ){
		z[i] = -r[i]/s[i] - u[i]*y[yidx[i]];
	}
}


/** 
 * Computes z = -r./s + u.*y(y)
 * where all vectors except of y are of length 6 (length of y >= 6).
 */
void armStateMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_6(armStateMPC_FLOAT *r, armStateMPC_FLOAT *s, armStateMPC_FLOAT *u, armStateMPC_FLOAT *y, int* yidx, armStateMPC_FLOAT *z)
{
	int i;
	for( i=0; i<6; i++ ){
		z[i] = -r[i]/s[i] + u[i]*y[yidx[i]];
	}
}


/*
 * Computes ds = -l.\(r + s.*dl) for vectors of length 228.
 */
void armStateMPC_LA_VSUB7_228(armStateMPC_FLOAT *l, armStateMPC_FLOAT *r, armStateMPC_FLOAT *s, armStateMPC_FLOAT *dl, armStateMPC_FLOAT *ds)
{
	int i;
	for( i=0; i<228; i++){
		ds[i] = -(r[i] + s[i]*dl[i])/l[i];
	}
}


/*
 * Vector addition x = x + y for vectors of length 114.
 */
void armStateMPC_LA_VADD_114(armStateMPC_FLOAT *x, armStateMPC_FLOAT *y)
{
	int i;
	for( i=0; i<114; i++){
		x[i] += y[i];
	}
}


/*
 * Vector addition x = x + y for vectors of length 60.
 */
void armStateMPC_LA_VADD_60(armStateMPC_FLOAT *x, armStateMPC_FLOAT *y)
{
	int i;
	for( i=0; i<60; i++){
		x[i] += y[i];
	}
}


/*
 * Vector addition x = x + y for vectors of length 228.
 */
void armStateMPC_LA_VADD_228(armStateMPC_FLOAT *x, armStateMPC_FLOAT *y)
{
	int i;
	for( i=0; i<228; i++){
		x[i] += y[i];
	}
}


/**
 * Backtracking line search for combined predictor/corrector step.
 * Update on variables with safety factor gamma (to keep us away from
 * boundary).
 */
int armStateMPC_LINESEARCH_BACKTRACKING_COMBINED(armStateMPC_FLOAT *z, armStateMPC_FLOAT *v, armStateMPC_FLOAT *l, armStateMPC_FLOAT *s, armStateMPC_FLOAT *dz, armStateMPC_FLOAT *dv, armStateMPC_FLOAT *dl, armStateMPC_FLOAT *ds, armStateMPC_FLOAT *a, armStateMPC_FLOAT *mu)
{
    int i, lsIt=1;       
    armStateMPC_FLOAT dltemp;
    armStateMPC_FLOAT dstemp;    
    armStateMPC_FLOAT a_gamma;
            
    *a = 1.0;
    while( 1 ){                        

        /* check whether search criterion is fulfilled */
        for( i=0; i<228; i++ ){
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
        if( i == 228 ){
            break;
        } else {
            *a *= armStateMPC_SET_LS_SCALE;
            if( *a < armStateMPC_SET_LS_MINSTEP ){
                return armStateMPC_NOPROGRESS;
            }
        }
    }
    
    /* update variables with safety margin */
    a_gamma = (*a)*armStateMPC_SET_LS_MAXSTEP;
    
    /* primal variables */
    for( i=0; i<114; i++ ){
        z[i] += a_gamma*dz[i];
    }
    
    /* equality constraint multipliers */
    for( i=0; i<60; i++ ){
        v[i] += a_gamma*dv[i];
    }
    
    /* inequality constraint multipliers & slacks, also update mu */
    *mu = 0;
    for( i=0; i<228; i++ ){
        dltemp = l[i] + a_gamma*dl[i]; l[i] = dltemp;
        dstemp = s[i] + a_gamma*ds[i]; s[i] = dstemp;
        *mu += dltemp*dstemp;
    }
    
    *a = a_gamma;
    *mu /= (armStateMPC_FLOAT)228;
    return lsIt;
}




/* VARIABLE DEFINITIONS ------------------------------------------------ */
armStateMPC_FLOAT armStateMPC_z[114];
armStateMPC_FLOAT armStateMPC_v[60];
armStateMPC_FLOAT armStateMPC_dz_aff[114];
armStateMPC_FLOAT armStateMPC_dv_aff[60];
armStateMPC_FLOAT armStateMPC_grad_cost[114];
armStateMPC_FLOAT armStateMPC_grad_eq[114];
armStateMPC_FLOAT armStateMPC_rd[114];
armStateMPC_FLOAT armStateMPC_l[228];
armStateMPC_FLOAT armStateMPC_s[228];
armStateMPC_FLOAT armStateMPC_lbys[228];
armStateMPC_FLOAT armStateMPC_dl_aff[228];
armStateMPC_FLOAT armStateMPC_ds_aff[228];
armStateMPC_FLOAT armStateMPC_dz_cc[114];
armStateMPC_FLOAT armStateMPC_dv_cc[60];
armStateMPC_FLOAT armStateMPC_dl_cc[228];
armStateMPC_FLOAT armStateMPC_ds_cc[228];
armStateMPC_FLOAT armStateMPC_ccrhs[228];
armStateMPC_FLOAT armStateMPC_grad_ineq[114];
armStateMPC_FLOAT* armStateMPC_z0 = armStateMPC_z + 0;
armStateMPC_FLOAT* armStateMPC_dzaff0 = armStateMPC_dz_aff + 0;
armStateMPC_FLOAT* armStateMPC_dzcc0 = armStateMPC_dz_cc + 0;
armStateMPC_FLOAT* armStateMPC_rd0 = armStateMPC_rd + 0;
armStateMPC_FLOAT armStateMPC_Lbyrd0[12];
armStateMPC_FLOAT* armStateMPC_grad_cost0 = armStateMPC_grad_cost + 0;
armStateMPC_FLOAT* armStateMPC_grad_eq0 = armStateMPC_grad_eq + 0;
armStateMPC_FLOAT* armStateMPC_grad_ineq0 = armStateMPC_grad_ineq + 0;
armStateMPC_FLOAT armStateMPC_ctv0[12];
armStateMPC_FLOAT* armStateMPC_v0 = armStateMPC_v + 0;
armStateMPC_FLOAT armStateMPC_re0[6];
armStateMPC_FLOAT armStateMPC_beta0[6];
armStateMPC_FLOAT armStateMPC_betacc0[6];
armStateMPC_FLOAT* armStateMPC_dvaff0 = armStateMPC_dv_aff + 0;
armStateMPC_FLOAT* armStateMPC_dvcc0 = armStateMPC_dv_cc + 0;
armStateMPC_FLOAT armStateMPC_V0[72];
armStateMPC_FLOAT armStateMPC_Yd0[21];
armStateMPC_FLOAT armStateMPC_Ld0[21];
armStateMPC_FLOAT armStateMPC_yy0[6];
armStateMPC_FLOAT armStateMPC_bmy0[6];
int armStateMPC_lbIdx0[12] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
armStateMPC_FLOAT* armStateMPC_llb0 = armStateMPC_l + 0;
armStateMPC_FLOAT* armStateMPC_slb0 = armStateMPC_s + 0;
armStateMPC_FLOAT* armStateMPC_llbbyslb0 = armStateMPC_lbys + 0;
armStateMPC_FLOAT armStateMPC_rilb0[12];
armStateMPC_FLOAT* armStateMPC_dllbaff0 = armStateMPC_dl_aff + 0;
armStateMPC_FLOAT* armStateMPC_dslbaff0 = armStateMPC_ds_aff + 0;
armStateMPC_FLOAT* armStateMPC_dllbcc0 = armStateMPC_dl_cc + 0;
armStateMPC_FLOAT* armStateMPC_dslbcc0 = armStateMPC_ds_cc + 0;
armStateMPC_FLOAT* armStateMPC_ccrhsl0 = armStateMPC_ccrhs + 0;
int armStateMPC_ubIdx0[12] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
armStateMPC_FLOAT* armStateMPC_lub0 = armStateMPC_l + 12;
armStateMPC_FLOAT* armStateMPC_sub0 = armStateMPC_s + 12;
armStateMPC_FLOAT* armStateMPC_lubbysub0 = armStateMPC_lbys + 12;
armStateMPC_FLOAT armStateMPC_riub0[12];
armStateMPC_FLOAT* armStateMPC_dlubaff0 = armStateMPC_dl_aff + 12;
armStateMPC_FLOAT* armStateMPC_dsubaff0 = armStateMPC_ds_aff + 12;
armStateMPC_FLOAT* armStateMPC_dlubcc0 = armStateMPC_dl_cc + 12;
armStateMPC_FLOAT* armStateMPC_dsubcc0 = armStateMPC_ds_cc + 12;
armStateMPC_FLOAT* armStateMPC_ccrhsub0 = armStateMPC_ccrhs + 12;
armStateMPC_FLOAT armStateMPC_Phi0[12];
armStateMPC_FLOAT armStateMPC_D0[12] = {1.0000000000000000E+000, 
1.0000000000000000E+000, 
1.0000000000000000E+000, 
1.0000000000000000E+000, 
1.0000000000000000E+000, 
1.0000000000000000E+000};
armStateMPC_FLOAT armStateMPC_W0[12];
armStateMPC_FLOAT* armStateMPC_z1 = armStateMPC_z + 12;
armStateMPC_FLOAT* armStateMPC_dzaff1 = armStateMPC_dz_aff + 12;
armStateMPC_FLOAT* armStateMPC_dzcc1 = armStateMPC_dz_cc + 12;
armStateMPC_FLOAT* armStateMPC_rd1 = armStateMPC_rd + 12;
armStateMPC_FLOAT armStateMPC_Lbyrd1[12];
armStateMPC_FLOAT* armStateMPC_grad_cost1 = armStateMPC_grad_cost + 12;
armStateMPC_FLOAT* armStateMPC_grad_eq1 = armStateMPC_grad_eq + 12;
armStateMPC_FLOAT* armStateMPC_grad_ineq1 = armStateMPC_grad_ineq + 12;
armStateMPC_FLOAT armStateMPC_ctv1[12];
armStateMPC_FLOAT* armStateMPC_v1 = armStateMPC_v + 6;
armStateMPC_FLOAT armStateMPC_re1[6];
armStateMPC_FLOAT armStateMPC_beta1[6];
armStateMPC_FLOAT armStateMPC_betacc1[6];
armStateMPC_FLOAT* armStateMPC_dvaff1 = armStateMPC_dv_aff + 6;
armStateMPC_FLOAT* armStateMPC_dvcc1 = armStateMPC_dv_cc + 6;
armStateMPC_FLOAT armStateMPC_V1[72];
armStateMPC_FLOAT armStateMPC_Yd1[21];
armStateMPC_FLOAT armStateMPC_Ld1[21];
armStateMPC_FLOAT armStateMPC_yy1[6];
armStateMPC_FLOAT armStateMPC_bmy1[6];
int armStateMPC_lbIdx1[12] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
armStateMPC_FLOAT* armStateMPC_llb1 = armStateMPC_l + 24;
armStateMPC_FLOAT* armStateMPC_slb1 = armStateMPC_s + 24;
armStateMPC_FLOAT* armStateMPC_llbbyslb1 = armStateMPC_lbys + 24;
armStateMPC_FLOAT armStateMPC_rilb1[12];
armStateMPC_FLOAT* armStateMPC_dllbaff1 = armStateMPC_dl_aff + 24;
armStateMPC_FLOAT* armStateMPC_dslbaff1 = armStateMPC_ds_aff + 24;
armStateMPC_FLOAT* armStateMPC_dllbcc1 = armStateMPC_dl_cc + 24;
armStateMPC_FLOAT* armStateMPC_dslbcc1 = armStateMPC_ds_cc + 24;
armStateMPC_FLOAT* armStateMPC_ccrhsl1 = armStateMPC_ccrhs + 24;
int armStateMPC_ubIdx1[12] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
armStateMPC_FLOAT* armStateMPC_lub1 = armStateMPC_l + 36;
armStateMPC_FLOAT* armStateMPC_sub1 = armStateMPC_s + 36;
armStateMPC_FLOAT* armStateMPC_lubbysub1 = armStateMPC_lbys + 36;
armStateMPC_FLOAT armStateMPC_riub1[12];
armStateMPC_FLOAT* armStateMPC_dlubaff1 = armStateMPC_dl_aff + 36;
armStateMPC_FLOAT* armStateMPC_dsubaff1 = armStateMPC_ds_aff + 36;
armStateMPC_FLOAT* armStateMPC_dlubcc1 = armStateMPC_dl_cc + 36;
armStateMPC_FLOAT* armStateMPC_dsubcc1 = armStateMPC_ds_cc + 36;
armStateMPC_FLOAT* armStateMPC_ccrhsub1 = armStateMPC_ccrhs + 36;
armStateMPC_FLOAT armStateMPC_Phi1[12];
armStateMPC_FLOAT armStateMPC_D1[12] = {-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000};
armStateMPC_FLOAT armStateMPC_W1[12];
armStateMPC_FLOAT armStateMPC_Ysd1[36];
armStateMPC_FLOAT armStateMPC_Lsd1[36];
armStateMPC_FLOAT* armStateMPC_z2 = armStateMPC_z + 24;
armStateMPC_FLOAT* armStateMPC_dzaff2 = armStateMPC_dz_aff + 24;
armStateMPC_FLOAT* armStateMPC_dzcc2 = armStateMPC_dz_cc + 24;
armStateMPC_FLOAT* armStateMPC_rd2 = armStateMPC_rd + 24;
armStateMPC_FLOAT armStateMPC_Lbyrd2[12];
armStateMPC_FLOAT* armStateMPC_grad_cost2 = armStateMPC_grad_cost + 24;
armStateMPC_FLOAT* armStateMPC_grad_eq2 = armStateMPC_grad_eq + 24;
armStateMPC_FLOAT* armStateMPC_grad_ineq2 = armStateMPC_grad_ineq + 24;
armStateMPC_FLOAT armStateMPC_ctv2[12];
armStateMPC_FLOAT* armStateMPC_v2 = armStateMPC_v + 12;
armStateMPC_FLOAT armStateMPC_re2[6];
armStateMPC_FLOAT armStateMPC_beta2[6];
armStateMPC_FLOAT armStateMPC_betacc2[6];
armStateMPC_FLOAT* armStateMPC_dvaff2 = armStateMPC_dv_aff + 12;
armStateMPC_FLOAT* armStateMPC_dvcc2 = armStateMPC_dv_cc + 12;
armStateMPC_FLOAT armStateMPC_V2[72];
armStateMPC_FLOAT armStateMPC_Yd2[21];
armStateMPC_FLOAT armStateMPC_Ld2[21];
armStateMPC_FLOAT armStateMPC_yy2[6];
armStateMPC_FLOAT armStateMPC_bmy2[6];
int armStateMPC_lbIdx2[12] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
armStateMPC_FLOAT* armStateMPC_llb2 = armStateMPC_l + 48;
armStateMPC_FLOAT* armStateMPC_slb2 = armStateMPC_s + 48;
armStateMPC_FLOAT* armStateMPC_llbbyslb2 = armStateMPC_lbys + 48;
armStateMPC_FLOAT armStateMPC_rilb2[12];
armStateMPC_FLOAT* armStateMPC_dllbaff2 = armStateMPC_dl_aff + 48;
armStateMPC_FLOAT* armStateMPC_dslbaff2 = armStateMPC_ds_aff + 48;
armStateMPC_FLOAT* armStateMPC_dllbcc2 = armStateMPC_dl_cc + 48;
armStateMPC_FLOAT* armStateMPC_dslbcc2 = armStateMPC_ds_cc + 48;
armStateMPC_FLOAT* armStateMPC_ccrhsl2 = armStateMPC_ccrhs + 48;
int armStateMPC_ubIdx2[12] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
armStateMPC_FLOAT* armStateMPC_lub2 = armStateMPC_l + 60;
armStateMPC_FLOAT* armStateMPC_sub2 = armStateMPC_s + 60;
armStateMPC_FLOAT* armStateMPC_lubbysub2 = armStateMPC_lbys + 60;
armStateMPC_FLOAT armStateMPC_riub2[12];
armStateMPC_FLOAT* armStateMPC_dlubaff2 = armStateMPC_dl_aff + 60;
armStateMPC_FLOAT* armStateMPC_dsubaff2 = armStateMPC_ds_aff + 60;
armStateMPC_FLOAT* armStateMPC_dlubcc2 = armStateMPC_dl_cc + 60;
armStateMPC_FLOAT* armStateMPC_dsubcc2 = armStateMPC_ds_cc + 60;
armStateMPC_FLOAT* armStateMPC_ccrhsub2 = armStateMPC_ccrhs + 60;
armStateMPC_FLOAT armStateMPC_Phi2[12];
armStateMPC_FLOAT armStateMPC_W2[12];
armStateMPC_FLOAT armStateMPC_Ysd2[36];
armStateMPC_FLOAT armStateMPC_Lsd2[36];
armStateMPC_FLOAT* armStateMPC_z3 = armStateMPC_z + 36;
armStateMPC_FLOAT* armStateMPC_dzaff3 = armStateMPC_dz_aff + 36;
armStateMPC_FLOAT* armStateMPC_dzcc3 = armStateMPC_dz_cc + 36;
armStateMPC_FLOAT* armStateMPC_rd3 = armStateMPC_rd + 36;
armStateMPC_FLOAT armStateMPC_Lbyrd3[12];
armStateMPC_FLOAT* armStateMPC_grad_cost3 = armStateMPC_grad_cost + 36;
armStateMPC_FLOAT* armStateMPC_grad_eq3 = armStateMPC_grad_eq + 36;
armStateMPC_FLOAT* armStateMPC_grad_ineq3 = armStateMPC_grad_ineq + 36;
armStateMPC_FLOAT armStateMPC_ctv3[12];
armStateMPC_FLOAT* armStateMPC_v3 = armStateMPC_v + 18;
armStateMPC_FLOAT armStateMPC_re3[6];
armStateMPC_FLOAT armStateMPC_beta3[6];
armStateMPC_FLOAT armStateMPC_betacc3[6];
armStateMPC_FLOAT* armStateMPC_dvaff3 = armStateMPC_dv_aff + 18;
armStateMPC_FLOAT* armStateMPC_dvcc3 = armStateMPC_dv_cc + 18;
armStateMPC_FLOAT armStateMPC_V3[72];
armStateMPC_FLOAT armStateMPC_Yd3[21];
armStateMPC_FLOAT armStateMPC_Ld3[21];
armStateMPC_FLOAT armStateMPC_yy3[6];
armStateMPC_FLOAT armStateMPC_bmy3[6];
int armStateMPC_lbIdx3[12] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
armStateMPC_FLOAT* armStateMPC_llb3 = armStateMPC_l + 72;
armStateMPC_FLOAT* armStateMPC_slb3 = armStateMPC_s + 72;
armStateMPC_FLOAT* armStateMPC_llbbyslb3 = armStateMPC_lbys + 72;
armStateMPC_FLOAT armStateMPC_rilb3[12];
armStateMPC_FLOAT* armStateMPC_dllbaff3 = armStateMPC_dl_aff + 72;
armStateMPC_FLOAT* armStateMPC_dslbaff3 = armStateMPC_ds_aff + 72;
armStateMPC_FLOAT* armStateMPC_dllbcc3 = armStateMPC_dl_cc + 72;
armStateMPC_FLOAT* armStateMPC_dslbcc3 = armStateMPC_ds_cc + 72;
armStateMPC_FLOAT* armStateMPC_ccrhsl3 = armStateMPC_ccrhs + 72;
int armStateMPC_ubIdx3[12] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
armStateMPC_FLOAT* armStateMPC_lub3 = armStateMPC_l + 84;
armStateMPC_FLOAT* armStateMPC_sub3 = armStateMPC_s + 84;
armStateMPC_FLOAT* armStateMPC_lubbysub3 = armStateMPC_lbys + 84;
armStateMPC_FLOAT armStateMPC_riub3[12];
armStateMPC_FLOAT* armStateMPC_dlubaff3 = armStateMPC_dl_aff + 84;
armStateMPC_FLOAT* armStateMPC_dsubaff3 = armStateMPC_ds_aff + 84;
armStateMPC_FLOAT* armStateMPC_dlubcc3 = armStateMPC_dl_cc + 84;
armStateMPC_FLOAT* armStateMPC_dsubcc3 = armStateMPC_ds_cc + 84;
armStateMPC_FLOAT* armStateMPC_ccrhsub3 = armStateMPC_ccrhs + 84;
armStateMPC_FLOAT armStateMPC_Phi3[12];
armStateMPC_FLOAT armStateMPC_W3[12];
armStateMPC_FLOAT armStateMPC_Ysd3[36];
armStateMPC_FLOAT armStateMPC_Lsd3[36];
armStateMPC_FLOAT* armStateMPC_z4 = armStateMPC_z + 48;
armStateMPC_FLOAT* armStateMPC_dzaff4 = armStateMPC_dz_aff + 48;
armStateMPC_FLOAT* armStateMPC_dzcc4 = armStateMPC_dz_cc + 48;
armStateMPC_FLOAT* armStateMPC_rd4 = armStateMPC_rd + 48;
armStateMPC_FLOAT armStateMPC_Lbyrd4[12];
armStateMPC_FLOAT* armStateMPC_grad_cost4 = armStateMPC_grad_cost + 48;
armStateMPC_FLOAT* armStateMPC_grad_eq4 = armStateMPC_grad_eq + 48;
armStateMPC_FLOAT* armStateMPC_grad_ineq4 = armStateMPC_grad_ineq + 48;
armStateMPC_FLOAT armStateMPC_ctv4[12];
armStateMPC_FLOAT* armStateMPC_v4 = armStateMPC_v + 24;
armStateMPC_FLOAT armStateMPC_re4[6];
armStateMPC_FLOAT armStateMPC_beta4[6];
armStateMPC_FLOAT armStateMPC_betacc4[6];
armStateMPC_FLOAT* armStateMPC_dvaff4 = armStateMPC_dv_aff + 24;
armStateMPC_FLOAT* armStateMPC_dvcc4 = armStateMPC_dv_cc + 24;
armStateMPC_FLOAT armStateMPC_V4[72];
armStateMPC_FLOAT armStateMPC_Yd4[21];
armStateMPC_FLOAT armStateMPC_Ld4[21];
armStateMPC_FLOAT armStateMPC_yy4[6];
armStateMPC_FLOAT armStateMPC_bmy4[6];
int armStateMPC_lbIdx4[12] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
armStateMPC_FLOAT* armStateMPC_llb4 = armStateMPC_l + 96;
armStateMPC_FLOAT* armStateMPC_slb4 = armStateMPC_s + 96;
armStateMPC_FLOAT* armStateMPC_llbbyslb4 = armStateMPC_lbys + 96;
armStateMPC_FLOAT armStateMPC_rilb4[12];
armStateMPC_FLOAT* armStateMPC_dllbaff4 = armStateMPC_dl_aff + 96;
armStateMPC_FLOAT* armStateMPC_dslbaff4 = armStateMPC_ds_aff + 96;
armStateMPC_FLOAT* armStateMPC_dllbcc4 = armStateMPC_dl_cc + 96;
armStateMPC_FLOAT* armStateMPC_dslbcc4 = armStateMPC_ds_cc + 96;
armStateMPC_FLOAT* armStateMPC_ccrhsl4 = armStateMPC_ccrhs + 96;
int armStateMPC_ubIdx4[12] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
armStateMPC_FLOAT* armStateMPC_lub4 = armStateMPC_l + 108;
armStateMPC_FLOAT* armStateMPC_sub4 = armStateMPC_s + 108;
armStateMPC_FLOAT* armStateMPC_lubbysub4 = armStateMPC_lbys + 108;
armStateMPC_FLOAT armStateMPC_riub4[12];
armStateMPC_FLOAT* armStateMPC_dlubaff4 = armStateMPC_dl_aff + 108;
armStateMPC_FLOAT* armStateMPC_dsubaff4 = armStateMPC_ds_aff + 108;
armStateMPC_FLOAT* armStateMPC_dlubcc4 = armStateMPC_dl_cc + 108;
armStateMPC_FLOAT* armStateMPC_dsubcc4 = armStateMPC_ds_cc + 108;
armStateMPC_FLOAT* armStateMPC_ccrhsub4 = armStateMPC_ccrhs + 108;
armStateMPC_FLOAT armStateMPC_Phi4[12];
armStateMPC_FLOAT armStateMPC_W4[12];
armStateMPC_FLOAT armStateMPC_Ysd4[36];
armStateMPC_FLOAT armStateMPC_Lsd4[36];
armStateMPC_FLOAT* armStateMPC_z5 = armStateMPC_z + 60;
armStateMPC_FLOAT* armStateMPC_dzaff5 = armStateMPC_dz_aff + 60;
armStateMPC_FLOAT* armStateMPC_dzcc5 = armStateMPC_dz_cc + 60;
armStateMPC_FLOAT* armStateMPC_rd5 = armStateMPC_rd + 60;
armStateMPC_FLOAT armStateMPC_Lbyrd5[12];
armStateMPC_FLOAT* armStateMPC_grad_cost5 = armStateMPC_grad_cost + 60;
armStateMPC_FLOAT* armStateMPC_grad_eq5 = armStateMPC_grad_eq + 60;
armStateMPC_FLOAT* armStateMPC_grad_ineq5 = armStateMPC_grad_ineq + 60;
armStateMPC_FLOAT armStateMPC_ctv5[12];
armStateMPC_FLOAT* armStateMPC_v5 = armStateMPC_v + 30;
armStateMPC_FLOAT armStateMPC_re5[6];
armStateMPC_FLOAT armStateMPC_beta5[6];
armStateMPC_FLOAT armStateMPC_betacc5[6];
armStateMPC_FLOAT* armStateMPC_dvaff5 = armStateMPC_dv_aff + 30;
armStateMPC_FLOAT* armStateMPC_dvcc5 = armStateMPC_dv_cc + 30;
armStateMPC_FLOAT armStateMPC_V5[72];
armStateMPC_FLOAT armStateMPC_Yd5[21];
armStateMPC_FLOAT armStateMPC_Ld5[21];
armStateMPC_FLOAT armStateMPC_yy5[6];
armStateMPC_FLOAT armStateMPC_bmy5[6];
int armStateMPC_lbIdx5[12] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
armStateMPC_FLOAT* armStateMPC_llb5 = armStateMPC_l + 120;
armStateMPC_FLOAT* armStateMPC_slb5 = armStateMPC_s + 120;
armStateMPC_FLOAT* armStateMPC_llbbyslb5 = armStateMPC_lbys + 120;
armStateMPC_FLOAT armStateMPC_rilb5[12];
armStateMPC_FLOAT* armStateMPC_dllbaff5 = armStateMPC_dl_aff + 120;
armStateMPC_FLOAT* armStateMPC_dslbaff5 = armStateMPC_ds_aff + 120;
armStateMPC_FLOAT* armStateMPC_dllbcc5 = armStateMPC_dl_cc + 120;
armStateMPC_FLOAT* armStateMPC_dslbcc5 = armStateMPC_ds_cc + 120;
armStateMPC_FLOAT* armStateMPC_ccrhsl5 = armStateMPC_ccrhs + 120;
int armStateMPC_ubIdx5[12] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
armStateMPC_FLOAT* armStateMPC_lub5 = armStateMPC_l + 132;
armStateMPC_FLOAT* armStateMPC_sub5 = armStateMPC_s + 132;
armStateMPC_FLOAT* armStateMPC_lubbysub5 = armStateMPC_lbys + 132;
armStateMPC_FLOAT armStateMPC_riub5[12];
armStateMPC_FLOAT* armStateMPC_dlubaff5 = armStateMPC_dl_aff + 132;
armStateMPC_FLOAT* armStateMPC_dsubaff5 = armStateMPC_ds_aff + 132;
armStateMPC_FLOAT* armStateMPC_dlubcc5 = armStateMPC_dl_cc + 132;
armStateMPC_FLOAT* armStateMPC_dsubcc5 = armStateMPC_ds_cc + 132;
armStateMPC_FLOAT* armStateMPC_ccrhsub5 = armStateMPC_ccrhs + 132;
armStateMPC_FLOAT armStateMPC_Phi5[12];
armStateMPC_FLOAT armStateMPC_W5[12];
armStateMPC_FLOAT armStateMPC_Ysd5[36];
armStateMPC_FLOAT armStateMPC_Lsd5[36];
armStateMPC_FLOAT* armStateMPC_z6 = armStateMPC_z + 72;
armStateMPC_FLOAT* armStateMPC_dzaff6 = armStateMPC_dz_aff + 72;
armStateMPC_FLOAT* armStateMPC_dzcc6 = armStateMPC_dz_cc + 72;
armStateMPC_FLOAT* armStateMPC_rd6 = armStateMPC_rd + 72;
armStateMPC_FLOAT armStateMPC_Lbyrd6[12];
armStateMPC_FLOAT* armStateMPC_grad_cost6 = armStateMPC_grad_cost + 72;
armStateMPC_FLOAT* armStateMPC_grad_eq6 = armStateMPC_grad_eq + 72;
armStateMPC_FLOAT* armStateMPC_grad_ineq6 = armStateMPC_grad_ineq + 72;
armStateMPC_FLOAT armStateMPC_ctv6[12];
armStateMPC_FLOAT* armStateMPC_v6 = armStateMPC_v + 36;
armStateMPC_FLOAT armStateMPC_re6[6];
armStateMPC_FLOAT armStateMPC_beta6[6];
armStateMPC_FLOAT armStateMPC_betacc6[6];
armStateMPC_FLOAT* armStateMPC_dvaff6 = armStateMPC_dv_aff + 36;
armStateMPC_FLOAT* armStateMPC_dvcc6 = armStateMPC_dv_cc + 36;
armStateMPC_FLOAT armStateMPC_V6[72];
armStateMPC_FLOAT armStateMPC_Yd6[21];
armStateMPC_FLOAT armStateMPC_Ld6[21];
armStateMPC_FLOAT armStateMPC_yy6[6];
armStateMPC_FLOAT armStateMPC_bmy6[6];
int armStateMPC_lbIdx6[12] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
armStateMPC_FLOAT* armStateMPC_llb6 = armStateMPC_l + 144;
armStateMPC_FLOAT* armStateMPC_slb6 = armStateMPC_s + 144;
armStateMPC_FLOAT* armStateMPC_llbbyslb6 = armStateMPC_lbys + 144;
armStateMPC_FLOAT armStateMPC_rilb6[12];
armStateMPC_FLOAT* armStateMPC_dllbaff6 = armStateMPC_dl_aff + 144;
armStateMPC_FLOAT* armStateMPC_dslbaff6 = armStateMPC_ds_aff + 144;
armStateMPC_FLOAT* armStateMPC_dllbcc6 = armStateMPC_dl_cc + 144;
armStateMPC_FLOAT* armStateMPC_dslbcc6 = armStateMPC_ds_cc + 144;
armStateMPC_FLOAT* armStateMPC_ccrhsl6 = armStateMPC_ccrhs + 144;
int armStateMPC_ubIdx6[12] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
armStateMPC_FLOAT* armStateMPC_lub6 = armStateMPC_l + 156;
armStateMPC_FLOAT* armStateMPC_sub6 = armStateMPC_s + 156;
armStateMPC_FLOAT* armStateMPC_lubbysub6 = armStateMPC_lbys + 156;
armStateMPC_FLOAT armStateMPC_riub6[12];
armStateMPC_FLOAT* armStateMPC_dlubaff6 = armStateMPC_dl_aff + 156;
armStateMPC_FLOAT* armStateMPC_dsubaff6 = armStateMPC_ds_aff + 156;
armStateMPC_FLOAT* armStateMPC_dlubcc6 = armStateMPC_dl_cc + 156;
armStateMPC_FLOAT* armStateMPC_dsubcc6 = armStateMPC_ds_cc + 156;
armStateMPC_FLOAT* armStateMPC_ccrhsub6 = armStateMPC_ccrhs + 156;
armStateMPC_FLOAT armStateMPC_Phi6[12];
armStateMPC_FLOAT armStateMPC_W6[12];
armStateMPC_FLOAT armStateMPC_Ysd6[36];
armStateMPC_FLOAT armStateMPC_Lsd6[36];
armStateMPC_FLOAT* armStateMPC_z7 = armStateMPC_z + 84;
armStateMPC_FLOAT* armStateMPC_dzaff7 = armStateMPC_dz_aff + 84;
armStateMPC_FLOAT* armStateMPC_dzcc7 = armStateMPC_dz_cc + 84;
armStateMPC_FLOAT* armStateMPC_rd7 = armStateMPC_rd + 84;
armStateMPC_FLOAT armStateMPC_Lbyrd7[12];
armStateMPC_FLOAT* armStateMPC_grad_cost7 = armStateMPC_grad_cost + 84;
armStateMPC_FLOAT* armStateMPC_grad_eq7 = armStateMPC_grad_eq + 84;
armStateMPC_FLOAT* armStateMPC_grad_ineq7 = armStateMPC_grad_ineq + 84;
armStateMPC_FLOAT armStateMPC_ctv7[12];
armStateMPC_FLOAT* armStateMPC_v7 = armStateMPC_v + 42;
armStateMPC_FLOAT armStateMPC_re7[6];
armStateMPC_FLOAT armStateMPC_beta7[6];
armStateMPC_FLOAT armStateMPC_betacc7[6];
armStateMPC_FLOAT* armStateMPC_dvaff7 = armStateMPC_dv_aff + 42;
armStateMPC_FLOAT* armStateMPC_dvcc7 = armStateMPC_dv_cc + 42;
armStateMPC_FLOAT armStateMPC_V7[72];
armStateMPC_FLOAT armStateMPC_Yd7[21];
armStateMPC_FLOAT armStateMPC_Ld7[21];
armStateMPC_FLOAT armStateMPC_yy7[6];
armStateMPC_FLOAT armStateMPC_bmy7[6];
int armStateMPC_lbIdx7[12] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
armStateMPC_FLOAT* armStateMPC_llb7 = armStateMPC_l + 168;
armStateMPC_FLOAT* armStateMPC_slb7 = armStateMPC_s + 168;
armStateMPC_FLOAT* armStateMPC_llbbyslb7 = armStateMPC_lbys + 168;
armStateMPC_FLOAT armStateMPC_rilb7[12];
armStateMPC_FLOAT* armStateMPC_dllbaff7 = armStateMPC_dl_aff + 168;
armStateMPC_FLOAT* armStateMPC_dslbaff7 = armStateMPC_ds_aff + 168;
armStateMPC_FLOAT* armStateMPC_dllbcc7 = armStateMPC_dl_cc + 168;
armStateMPC_FLOAT* armStateMPC_dslbcc7 = armStateMPC_ds_cc + 168;
armStateMPC_FLOAT* armStateMPC_ccrhsl7 = armStateMPC_ccrhs + 168;
int armStateMPC_ubIdx7[12] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
armStateMPC_FLOAT* armStateMPC_lub7 = armStateMPC_l + 180;
armStateMPC_FLOAT* armStateMPC_sub7 = armStateMPC_s + 180;
armStateMPC_FLOAT* armStateMPC_lubbysub7 = armStateMPC_lbys + 180;
armStateMPC_FLOAT armStateMPC_riub7[12];
armStateMPC_FLOAT* armStateMPC_dlubaff7 = armStateMPC_dl_aff + 180;
armStateMPC_FLOAT* armStateMPC_dsubaff7 = armStateMPC_ds_aff + 180;
armStateMPC_FLOAT* armStateMPC_dlubcc7 = armStateMPC_dl_cc + 180;
armStateMPC_FLOAT* armStateMPC_dsubcc7 = armStateMPC_ds_cc + 180;
armStateMPC_FLOAT* armStateMPC_ccrhsub7 = armStateMPC_ccrhs + 180;
armStateMPC_FLOAT armStateMPC_Phi7[12];
armStateMPC_FLOAT armStateMPC_W7[12];
armStateMPC_FLOAT armStateMPC_Ysd7[36];
armStateMPC_FLOAT armStateMPC_Lsd7[36];
armStateMPC_FLOAT* armStateMPC_z8 = armStateMPC_z + 96;
armStateMPC_FLOAT* armStateMPC_dzaff8 = armStateMPC_dz_aff + 96;
armStateMPC_FLOAT* armStateMPC_dzcc8 = armStateMPC_dz_cc + 96;
armStateMPC_FLOAT* armStateMPC_rd8 = armStateMPC_rd + 96;
armStateMPC_FLOAT armStateMPC_Lbyrd8[12];
armStateMPC_FLOAT* armStateMPC_grad_cost8 = armStateMPC_grad_cost + 96;
armStateMPC_FLOAT* armStateMPC_grad_eq8 = armStateMPC_grad_eq + 96;
armStateMPC_FLOAT* armStateMPC_grad_ineq8 = armStateMPC_grad_ineq + 96;
armStateMPC_FLOAT armStateMPC_ctv8[12];
armStateMPC_FLOAT* armStateMPC_v8 = armStateMPC_v + 48;
armStateMPC_FLOAT armStateMPC_re8[6];
armStateMPC_FLOAT armStateMPC_beta8[6];
armStateMPC_FLOAT armStateMPC_betacc8[6];
armStateMPC_FLOAT* armStateMPC_dvaff8 = armStateMPC_dv_aff + 48;
armStateMPC_FLOAT* armStateMPC_dvcc8 = armStateMPC_dv_cc + 48;
armStateMPC_FLOAT armStateMPC_V8[72];
armStateMPC_FLOAT armStateMPC_Yd8[21];
armStateMPC_FLOAT armStateMPC_Ld8[21];
armStateMPC_FLOAT armStateMPC_yy8[6];
armStateMPC_FLOAT armStateMPC_bmy8[6];
int armStateMPC_lbIdx8[12] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
armStateMPC_FLOAT* armStateMPC_llb8 = armStateMPC_l + 192;
armStateMPC_FLOAT* armStateMPC_slb8 = armStateMPC_s + 192;
armStateMPC_FLOAT* armStateMPC_llbbyslb8 = armStateMPC_lbys + 192;
armStateMPC_FLOAT armStateMPC_rilb8[12];
armStateMPC_FLOAT* armStateMPC_dllbaff8 = armStateMPC_dl_aff + 192;
armStateMPC_FLOAT* armStateMPC_dslbaff8 = armStateMPC_ds_aff + 192;
armStateMPC_FLOAT* armStateMPC_dllbcc8 = armStateMPC_dl_cc + 192;
armStateMPC_FLOAT* armStateMPC_dslbcc8 = armStateMPC_ds_cc + 192;
armStateMPC_FLOAT* armStateMPC_ccrhsl8 = armStateMPC_ccrhs + 192;
int armStateMPC_ubIdx8[12] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
armStateMPC_FLOAT* armStateMPC_lub8 = armStateMPC_l + 204;
armStateMPC_FLOAT* armStateMPC_sub8 = armStateMPC_s + 204;
armStateMPC_FLOAT* armStateMPC_lubbysub8 = armStateMPC_lbys + 204;
armStateMPC_FLOAT armStateMPC_riub8[12];
armStateMPC_FLOAT* armStateMPC_dlubaff8 = armStateMPC_dl_aff + 204;
armStateMPC_FLOAT* armStateMPC_dsubaff8 = armStateMPC_ds_aff + 204;
armStateMPC_FLOAT* armStateMPC_dlubcc8 = armStateMPC_dl_cc + 204;
armStateMPC_FLOAT* armStateMPC_dsubcc8 = armStateMPC_ds_cc + 204;
armStateMPC_FLOAT* armStateMPC_ccrhsub8 = armStateMPC_ccrhs + 204;
armStateMPC_FLOAT armStateMPC_Phi8[12];
armStateMPC_FLOAT armStateMPC_W8[12];
armStateMPC_FLOAT armStateMPC_Ysd8[36];
armStateMPC_FLOAT armStateMPC_Lsd8[36];
armStateMPC_FLOAT* armStateMPC_z9 = armStateMPC_z + 108;
armStateMPC_FLOAT* armStateMPC_dzaff9 = armStateMPC_dz_aff + 108;
armStateMPC_FLOAT* armStateMPC_dzcc9 = armStateMPC_dz_cc + 108;
armStateMPC_FLOAT* armStateMPC_rd9 = armStateMPC_rd + 108;
armStateMPC_FLOAT armStateMPC_Lbyrd9[6];
armStateMPC_FLOAT* armStateMPC_grad_cost9 = armStateMPC_grad_cost + 108;
armStateMPC_FLOAT* armStateMPC_grad_eq9 = armStateMPC_grad_eq + 108;
armStateMPC_FLOAT* armStateMPC_grad_ineq9 = armStateMPC_grad_ineq + 108;
armStateMPC_FLOAT armStateMPC_ctv9[6];
armStateMPC_FLOAT* armStateMPC_v9 = armStateMPC_v + 54;
armStateMPC_FLOAT armStateMPC_re9[6];
armStateMPC_FLOAT armStateMPC_beta9[6];
armStateMPC_FLOAT armStateMPC_betacc9[6];
armStateMPC_FLOAT* armStateMPC_dvaff9 = armStateMPC_dv_aff + 54;
armStateMPC_FLOAT* armStateMPC_dvcc9 = armStateMPC_dv_cc + 54;
armStateMPC_FLOAT armStateMPC_V9[36];
armStateMPC_FLOAT armStateMPC_Yd9[21];
armStateMPC_FLOAT armStateMPC_Ld9[21];
armStateMPC_FLOAT armStateMPC_yy9[6];
armStateMPC_FLOAT armStateMPC_bmy9[6];
int armStateMPC_lbIdx9[6] = {0, 1, 2, 3, 4, 5};
armStateMPC_FLOAT* armStateMPC_llb9 = armStateMPC_l + 216;
armStateMPC_FLOAT* armStateMPC_slb9 = armStateMPC_s + 216;
armStateMPC_FLOAT* armStateMPC_llbbyslb9 = armStateMPC_lbys + 216;
armStateMPC_FLOAT armStateMPC_rilb9[6];
armStateMPC_FLOAT* armStateMPC_dllbaff9 = armStateMPC_dl_aff + 216;
armStateMPC_FLOAT* armStateMPC_dslbaff9 = armStateMPC_ds_aff + 216;
armStateMPC_FLOAT* armStateMPC_dllbcc9 = armStateMPC_dl_cc + 216;
armStateMPC_FLOAT* armStateMPC_dslbcc9 = armStateMPC_ds_cc + 216;
armStateMPC_FLOAT* armStateMPC_ccrhsl9 = armStateMPC_ccrhs + 216;
int armStateMPC_ubIdx9[6] = {0, 1, 2, 3, 4, 5};
armStateMPC_FLOAT* armStateMPC_lub9 = armStateMPC_l + 222;
armStateMPC_FLOAT* armStateMPC_sub9 = armStateMPC_s + 222;
armStateMPC_FLOAT* armStateMPC_lubbysub9 = armStateMPC_lbys + 222;
armStateMPC_FLOAT armStateMPC_riub9[6];
armStateMPC_FLOAT* armStateMPC_dlubaff9 = armStateMPC_dl_aff + 222;
armStateMPC_FLOAT* armStateMPC_dsubaff9 = armStateMPC_ds_aff + 222;
armStateMPC_FLOAT* armStateMPC_dlubcc9 = armStateMPC_dl_cc + 222;
armStateMPC_FLOAT* armStateMPC_dsubcc9 = armStateMPC_ds_cc + 222;
armStateMPC_FLOAT* armStateMPC_ccrhsub9 = armStateMPC_ccrhs + 222;
armStateMPC_FLOAT armStateMPC_Phi9[6];
armStateMPC_FLOAT armStateMPC_D9[6] = {-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000};
armStateMPC_FLOAT armStateMPC_W9[6];
armStateMPC_FLOAT armStateMPC_Ysd9[36];
armStateMPC_FLOAT armStateMPC_Lsd9[36];
armStateMPC_FLOAT musigma;
armStateMPC_FLOAT sigma_3rdroot;
armStateMPC_FLOAT armStateMPC_Diag1_0[12];
armStateMPC_FLOAT armStateMPC_Diag2_0[12];
armStateMPC_FLOAT armStateMPC_L_0[66];




/* SOLVER CODE --------------------------------------------------------- */
int armStateMPC_solve(armStateMPC_params* params, armStateMPC_output* output, armStateMPC_info* info)
{	
int exitcode;

#if armStateMPC_SET_TIMING == 1
	armStateMPC_timer solvertimer;
	armStateMPC_tic(&solvertimer);
#endif
/* FUNCTION CALLS INTO LA LIBRARY -------------------------------------- */
info->it = 0;
armStateMPC_LA_INITIALIZEVECTOR_114(armStateMPC_z, 0);
armStateMPC_LA_INITIALIZEVECTOR_60(armStateMPC_v, 1);
armStateMPC_LA_INITIALIZEVECTOR_228(armStateMPC_l, 1);
armStateMPC_LA_INITIALIZEVECTOR_228(armStateMPC_s, 1);
info->mu = 0;
armStateMPC_LA_DOTACC_228(armStateMPC_l, armStateMPC_s, &info->mu);
info->mu /= 228;
while( 1 ){
info->pobj = 0;
armStateMPC_LA_DIAG_QUADFCN_12(params->H1, params->f1, armStateMPC_z0, armStateMPC_grad_cost0, &info->pobj);
armStateMPC_LA_DIAG_QUADFCN_12(params->H2, params->f2, armStateMPC_z1, armStateMPC_grad_cost1, &info->pobj);
armStateMPC_LA_DIAG_QUADFCN_12(params->H3, params->f3, armStateMPC_z2, armStateMPC_grad_cost2, &info->pobj);
armStateMPC_LA_DIAG_QUADFCN_12(params->H4, params->f4, armStateMPC_z3, armStateMPC_grad_cost3, &info->pobj);
armStateMPC_LA_DIAG_QUADFCN_12(params->H5, params->f5, armStateMPC_z4, armStateMPC_grad_cost4, &info->pobj);
armStateMPC_LA_DIAG_QUADFCN_12(params->H6, params->f6, armStateMPC_z5, armStateMPC_grad_cost5, &info->pobj);
armStateMPC_LA_DIAG_QUADFCN_12(params->H7, params->f7, armStateMPC_z6, armStateMPC_grad_cost6, &info->pobj);
armStateMPC_LA_DIAG_QUADFCN_12(params->H8, params->f8, armStateMPC_z7, armStateMPC_grad_cost7, &info->pobj);
armStateMPC_LA_DIAG_QUADFCN_12(params->H9, params->f9, armStateMPC_z8, armStateMPC_grad_cost8, &info->pobj);
armStateMPC_LA_DIAG_QUADFCN_6(params->H10, params->f10, armStateMPC_z9, armStateMPC_grad_cost9, &info->pobj);
info->res_eq = 0;
info->dgap = 0;
armStateMPC_LA_DIAGZERO_MVMSUB6_6(armStateMPC_D0, armStateMPC_z0, params->e1, armStateMPC_v0, armStateMPC_re0, &info->dgap, &info->res_eq);
armStateMPC_LA_DENSE_DIAGZERO_MVMSUB3_6_12_12(params->C1, armStateMPC_z0, armStateMPC_D1, armStateMPC_z1, params->e2, armStateMPC_v1, armStateMPC_re1, &info->dgap, &info->res_eq);
armStateMPC_LA_DENSE_DIAGZERO_MVMSUB3_6_12_12(params->C2, armStateMPC_z1, armStateMPC_D1, armStateMPC_z2, params->e3, armStateMPC_v2, armStateMPC_re2, &info->dgap, &info->res_eq);
armStateMPC_LA_DENSE_DIAGZERO_MVMSUB3_6_12_12(params->C3, armStateMPC_z2, armStateMPC_D1, armStateMPC_z3, params->e4, armStateMPC_v3, armStateMPC_re3, &info->dgap, &info->res_eq);
armStateMPC_LA_DENSE_DIAGZERO_MVMSUB3_6_12_12(params->C4, armStateMPC_z3, armStateMPC_D1, armStateMPC_z4, params->e5, armStateMPC_v4, armStateMPC_re4, &info->dgap, &info->res_eq);
armStateMPC_LA_DENSE_DIAGZERO_MVMSUB3_6_12_12(params->C5, armStateMPC_z4, armStateMPC_D1, armStateMPC_z5, params->e6, armStateMPC_v5, armStateMPC_re5, &info->dgap, &info->res_eq);
armStateMPC_LA_DENSE_DIAGZERO_MVMSUB3_6_12_12(params->C6, armStateMPC_z5, armStateMPC_D1, armStateMPC_z6, params->e7, armStateMPC_v6, armStateMPC_re6, &info->dgap, &info->res_eq);
armStateMPC_LA_DENSE_DIAGZERO_MVMSUB3_6_12_12(params->C7, armStateMPC_z6, armStateMPC_D1, armStateMPC_z7, params->e8, armStateMPC_v7, armStateMPC_re7, &info->dgap, &info->res_eq);
armStateMPC_LA_DENSE_DIAGZERO_MVMSUB3_6_12_12(params->C8, armStateMPC_z7, armStateMPC_D1, armStateMPC_z8, params->e9, armStateMPC_v8, armStateMPC_re8, &info->dgap, &info->res_eq);
armStateMPC_LA_DENSE_DIAGZERO_MVMSUB3_6_12_6(params->C9, armStateMPC_z8, armStateMPC_D9, armStateMPC_z9, params->e10, armStateMPC_v9, armStateMPC_re9, &info->dgap, &info->res_eq);
armStateMPC_LA_DENSE_DIAGZERO_MTVM2_6_12_6(params->C1, armStateMPC_v1, armStateMPC_D0, armStateMPC_v0, armStateMPC_grad_eq0);
armStateMPC_LA_DENSE_DIAGZERO_MTVM2_6_12_6(params->C2, armStateMPC_v2, armStateMPC_D1, armStateMPC_v1, armStateMPC_grad_eq1);
armStateMPC_LA_DENSE_DIAGZERO_MTVM2_6_12_6(params->C3, armStateMPC_v3, armStateMPC_D1, armStateMPC_v2, armStateMPC_grad_eq2);
armStateMPC_LA_DENSE_DIAGZERO_MTVM2_6_12_6(params->C4, armStateMPC_v4, armStateMPC_D1, armStateMPC_v3, armStateMPC_grad_eq3);
armStateMPC_LA_DENSE_DIAGZERO_MTVM2_6_12_6(params->C5, armStateMPC_v5, armStateMPC_D1, armStateMPC_v4, armStateMPC_grad_eq4);
armStateMPC_LA_DENSE_DIAGZERO_MTVM2_6_12_6(params->C6, armStateMPC_v6, armStateMPC_D1, armStateMPC_v5, armStateMPC_grad_eq5);
armStateMPC_LA_DENSE_DIAGZERO_MTVM2_6_12_6(params->C7, armStateMPC_v7, armStateMPC_D1, armStateMPC_v6, armStateMPC_grad_eq6);
armStateMPC_LA_DENSE_DIAGZERO_MTVM2_6_12_6(params->C8, armStateMPC_v8, armStateMPC_D1, armStateMPC_v7, armStateMPC_grad_eq7);
armStateMPC_LA_DENSE_DIAGZERO_MTVM2_6_12_6(params->C9, armStateMPC_v9, armStateMPC_D1, armStateMPC_v8, armStateMPC_grad_eq8);
armStateMPC_LA_DIAGZERO_MTVM_6_6(armStateMPC_D9, armStateMPC_v9, armStateMPC_grad_eq9);
info->res_ineq = 0;
armStateMPC_LA_VSUBADD3_12(params->lb1, armStateMPC_z0, armStateMPC_lbIdx0, armStateMPC_llb0, armStateMPC_slb0, armStateMPC_rilb0, &info->dgap, &info->res_ineq);
armStateMPC_LA_VSUBADD2_12(armStateMPC_z0, armStateMPC_ubIdx0, params->ub1, armStateMPC_lub0, armStateMPC_sub0, armStateMPC_riub0, &info->dgap, &info->res_ineq);
armStateMPC_LA_VSUBADD3_12(params->lb2, armStateMPC_z1, armStateMPC_lbIdx1, armStateMPC_llb1, armStateMPC_slb1, armStateMPC_rilb1, &info->dgap, &info->res_ineq);
armStateMPC_LA_VSUBADD2_12(armStateMPC_z1, armStateMPC_ubIdx1, params->ub2, armStateMPC_lub1, armStateMPC_sub1, armStateMPC_riub1, &info->dgap, &info->res_ineq);
armStateMPC_LA_VSUBADD3_12(params->lb3, armStateMPC_z2, armStateMPC_lbIdx2, armStateMPC_llb2, armStateMPC_slb2, armStateMPC_rilb2, &info->dgap, &info->res_ineq);
armStateMPC_LA_VSUBADD2_12(armStateMPC_z2, armStateMPC_ubIdx2, params->ub3, armStateMPC_lub2, armStateMPC_sub2, armStateMPC_riub2, &info->dgap, &info->res_ineq);
armStateMPC_LA_VSUBADD3_12(params->lb4, armStateMPC_z3, armStateMPC_lbIdx3, armStateMPC_llb3, armStateMPC_slb3, armStateMPC_rilb3, &info->dgap, &info->res_ineq);
armStateMPC_LA_VSUBADD2_12(armStateMPC_z3, armStateMPC_ubIdx3, params->ub4, armStateMPC_lub3, armStateMPC_sub3, armStateMPC_riub3, &info->dgap, &info->res_ineq);
armStateMPC_LA_VSUBADD3_12(params->lb5, armStateMPC_z4, armStateMPC_lbIdx4, armStateMPC_llb4, armStateMPC_slb4, armStateMPC_rilb4, &info->dgap, &info->res_ineq);
armStateMPC_LA_VSUBADD2_12(armStateMPC_z4, armStateMPC_ubIdx4, params->ub5, armStateMPC_lub4, armStateMPC_sub4, armStateMPC_riub4, &info->dgap, &info->res_ineq);
armStateMPC_LA_VSUBADD3_12(params->lb6, armStateMPC_z5, armStateMPC_lbIdx5, armStateMPC_llb5, armStateMPC_slb5, armStateMPC_rilb5, &info->dgap, &info->res_ineq);
armStateMPC_LA_VSUBADD2_12(armStateMPC_z5, armStateMPC_ubIdx5, params->ub6, armStateMPC_lub5, armStateMPC_sub5, armStateMPC_riub5, &info->dgap, &info->res_ineq);
armStateMPC_LA_VSUBADD3_12(params->lb7, armStateMPC_z6, armStateMPC_lbIdx6, armStateMPC_llb6, armStateMPC_slb6, armStateMPC_rilb6, &info->dgap, &info->res_ineq);
armStateMPC_LA_VSUBADD2_12(armStateMPC_z6, armStateMPC_ubIdx6, params->ub7, armStateMPC_lub6, armStateMPC_sub6, armStateMPC_riub6, &info->dgap, &info->res_ineq);
armStateMPC_LA_VSUBADD3_12(params->lb8, armStateMPC_z7, armStateMPC_lbIdx7, armStateMPC_llb7, armStateMPC_slb7, armStateMPC_rilb7, &info->dgap, &info->res_ineq);
armStateMPC_LA_VSUBADD2_12(armStateMPC_z7, armStateMPC_ubIdx7, params->ub8, armStateMPC_lub7, armStateMPC_sub7, armStateMPC_riub7, &info->dgap, &info->res_ineq);
armStateMPC_LA_VSUBADD3_12(params->lb9, armStateMPC_z8, armStateMPC_lbIdx8, armStateMPC_llb8, armStateMPC_slb8, armStateMPC_rilb8, &info->dgap, &info->res_ineq);
armStateMPC_LA_VSUBADD2_12(armStateMPC_z8, armStateMPC_ubIdx8, params->ub9, armStateMPC_lub8, armStateMPC_sub8, armStateMPC_riub8, &info->dgap, &info->res_ineq);
armStateMPC_LA_VSUBADD3_6(params->lb10, armStateMPC_z9, armStateMPC_lbIdx9, armStateMPC_llb9, armStateMPC_slb9, armStateMPC_rilb9, &info->dgap, &info->res_ineq);
armStateMPC_LA_VSUBADD2_6(armStateMPC_z9, armStateMPC_ubIdx9, params->ub10, armStateMPC_lub9, armStateMPC_sub9, armStateMPC_riub9, &info->dgap, &info->res_ineq);
armStateMPC_LA_INEQ_B_GRAD_12_12_12(armStateMPC_lub0, armStateMPC_sub0, armStateMPC_riub0, armStateMPC_llb0, armStateMPC_slb0, armStateMPC_rilb0, armStateMPC_lbIdx0, armStateMPC_ubIdx0, armStateMPC_grad_ineq0, armStateMPC_lubbysub0, armStateMPC_llbbyslb0);
armStateMPC_LA_INEQ_B_GRAD_12_12_12(armStateMPC_lub1, armStateMPC_sub1, armStateMPC_riub1, armStateMPC_llb1, armStateMPC_slb1, armStateMPC_rilb1, armStateMPC_lbIdx1, armStateMPC_ubIdx1, armStateMPC_grad_ineq1, armStateMPC_lubbysub1, armStateMPC_llbbyslb1);
armStateMPC_LA_INEQ_B_GRAD_12_12_12(armStateMPC_lub2, armStateMPC_sub2, armStateMPC_riub2, armStateMPC_llb2, armStateMPC_slb2, armStateMPC_rilb2, armStateMPC_lbIdx2, armStateMPC_ubIdx2, armStateMPC_grad_ineq2, armStateMPC_lubbysub2, armStateMPC_llbbyslb2);
armStateMPC_LA_INEQ_B_GRAD_12_12_12(armStateMPC_lub3, armStateMPC_sub3, armStateMPC_riub3, armStateMPC_llb3, armStateMPC_slb3, armStateMPC_rilb3, armStateMPC_lbIdx3, armStateMPC_ubIdx3, armStateMPC_grad_ineq3, armStateMPC_lubbysub3, armStateMPC_llbbyslb3);
armStateMPC_LA_INEQ_B_GRAD_12_12_12(armStateMPC_lub4, armStateMPC_sub4, armStateMPC_riub4, armStateMPC_llb4, armStateMPC_slb4, armStateMPC_rilb4, armStateMPC_lbIdx4, armStateMPC_ubIdx4, armStateMPC_grad_ineq4, armStateMPC_lubbysub4, armStateMPC_llbbyslb4);
armStateMPC_LA_INEQ_B_GRAD_12_12_12(armStateMPC_lub5, armStateMPC_sub5, armStateMPC_riub5, armStateMPC_llb5, armStateMPC_slb5, armStateMPC_rilb5, armStateMPC_lbIdx5, armStateMPC_ubIdx5, armStateMPC_grad_ineq5, armStateMPC_lubbysub5, armStateMPC_llbbyslb5);
armStateMPC_LA_INEQ_B_GRAD_12_12_12(armStateMPC_lub6, armStateMPC_sub6, armStateMPC_riub6, armStateMPC_llb6, armStateMPC_slb6, armStateMPC_rilb6, armStateMPC_lbIdx6, armStateMPC_ubIdx6, armStateMPC_grad_ineq6, armStateMPC_lubbysub6, armStateMPC_llbbyslb6);
armStateMPC_LA_INEQ_B_GRAD_12_12_12(armStateMPC_lub7, armStateMPC_sub7, armStateMPC_riub7, armStateMPC_llb7, armStateMPC_slb7, armStateMPC_rilb7, armStateMPC_lbIdx7, armStateMPC_ubIdx7, armStateMPC_grad_ineq7, armStateMPC_lubbysub7, armStateMPC_llbbyslb7);
armStateMPC_LA_INEQ_B_GRAD_12_12_12(armStateMPC_lub8, armStateMPC_sub8, armStateMPC_riub8, armStateMPC_llb8, armStateMPC_slb8, armStateMPC_rilb8, armStateMPC_lbIdx8, armStateMPC_ubIdx8, armStateMPC_grad_ineq8, armStateMPC_lubbysub8, armStateMPC_llbbyslb8);
armStateMPC_LA_INEQ_B_GRAD_6_6_6(armStateMPC_lub9, armStateMPC_sub9, armStateMPC_riub9, armStateMPC_llb9, armStateMPC_slb9, armStateMPC_rilb9, armStateMPC_lbIdx9, armStateMPC_ubIdx9, armStateMPC_grad_ineq9, armStateMPC_lubbysub9, armStateMPC_llbbyslb9);
info->dobj = info->pobj - info->dgap;
info->rdgap = info->pobj ? info->dgap / info->pobj : 1e6;
if( info->rdgap < 0 ) info->rdgap = -info->rdgap;
if( info->mu < armStateMPC_SET_ACC_KKTCOMPL
    && (info->rdgap < armStateMPC_SET_ACC_RDGAP || info->dgap < armStateMPC_SET_ACC_KKTCOMPL)
    && info->res_eq < armStateMPC_SET_ACC_RESEQ
    && info->res_ineq < armStateMPC_SET_ACC_RESINEQ ){
exitcode = armStateMPC_OPTIMAL; break; }
if( info->it == armStateMPC_SET_MAXIT ){
exitcode = armStateMPC_MAXITREACHED; break; }
armStateMPC_LA_VVADD3_114(armStateMPC_grad_cost, armStateMPC_grad_eq, armStateMPC_grad_ineq, armStateMPC_rd);
armStateMPC_LA_DIAG_CHOL_ONELOOP_LBUB_12_12_12(params->H1, armStateMPC_llbbyslb0, armStateMPC_lbIdx0, armStateMPC_lubbysub0, armStateMPC_ubIdx0, armStateMPC_Phi0);
armStateMPC_LA_DIAG_MATRIXFORWARDSUB_6_12(armStateMPC_Phi0, params->C1, armStateMPC_V0);
armStateMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_6_12(armStateMPC_Phi0, armStateMPC_D0, armStateMPC_W0);
armStateMPC_LA_DENSE_DIAGZERO_MMTM_6_12_6(armStateMPC_W0, armStateMPC_V0, armStateMPC_Ysd1);
armStateMPC_LA_DIAG_FORWARDSUB_12(armStateMPC_Phi0, armStateMPC_rd0, armStateMPC_Lbyrd0);
armStateMPC_LA_DIAG_CHOL_ONELOOP_LBUB_12_12_12(params->H2, armStateMPC_llbbyslb1, armStateMPC_lbIdx1, armStateMPC_lubbysub1, armStateMPC_ubIdx1, armStateMPC_Phi1);
armStateMPC_LA_DIAG_MATRIXFORWARDSUB_6_12(armStateMPC_Phi1, params->C2, armStateMPC_V1);
armStateMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_6_12(armStateMPC_Phi1, armStateMPC_D1, armStateMPC_W1);
armStateMPC_LA_DENSE_DIAGZERO_MMTM_6_12_6(armStateMPC_W1, armStateMPC_V1, armStateMPC_Ysd2);
armStateMPC_LA_DIAG_FORWARDSUB_12(armStateMPC_Phi1, armStateMPC_rd1, armStateMPC_Lbyrd1);
armStateMPC_LA_DIAG_CHOL_ONELOOP_LBUB_12_12_12(params->H3, armStateMPC_llbbyslb2, armStateMPC_lbIdx2, armStateMPC_lubbysub2, armStateMPC_ubIdx2, armStateMPC_Phi2);
armStateMPC_LA_DIAG_MATRIXFORWARDSUB_6_12(armStateMPC_Phi2, params->C3, armStateMPC_V2);
armStateMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_6_12(armStateMPC_Phi2, armStateMPC_D1, armStateMPC_W2);
armStateMPC_LA_DENSE_DIAGZERO_MMTM_6_12_6(armStateMPC_W2, armStateMPC_V2, armStateMPC_Ysd3);
armStateMPC_LA_DIAG_FORWARDSUB_12(armStateMPC_Phi2, armStateMPC_rd2, armStateMPC_Lbyrd2);
armStateMPC_LA_DIAG_CHOL_ONELOOP_LBUB_12_12_12(params->H4, armStateMPC_llbbyslb3, armStateMPC_lbIdx3, armStateMPC_lubbysub3, armStateMPC_ubIdx3, armStateMPC_Phi3);
armStateMPC_LA_DIAG_MATRIXFORWARDSUB_6_12(armStateMPC_Phi3, params->C4, armStateMPC_V3);
armStateMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_6_12(armStateMPC_Phi3, armStateMPC_D1, armStateMPC_W3);
armStateMPC_LA_DENSE_DIAGZERO_MMTM_6_12_6(armStateMPC_W3, armStateMPC_V3, armStateMPC_Ysd4);
armStateMPC_LA_DIAG_FORWARDSUB_12(armStateMPC_Phi3, armStateMPC_rd3, armStateMPC_Lbyrd3);
armStateMPC_LA_DIAG_CHOL_ONELOOP_LBUB_12_12_12(params->H5, armStateMPC_llbbyslb4, armStateMPC_lbIdx4, armStateMPC_lubbysub4, armStateMPC_ubIdx4, armStateMPC_Phi4);
armStateMPC_LA_DIAG_MATRIXFORWARDSUB_6_12(armStateMPC_Phi4, params->C5, armStateMPC_V4);
armStateMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_6_12(armStateMPC_Phi4, armStateMPC_D1, armStateMPC_W4);
armStateMPC_LA_DENSE_DIAGZERO_MMTM_6_12_6(armStateMPC_W4, armStateMPC_V4, armStateMPC_Ysd5);
armStateMPC_LA_DIAG_FORWARDSUB_12(armStateMPC_Phi4, armStateMPC_rd4, armStateMPC_Lbyrd4);
armStateMPC_LA_DIAG_CHOL_ONELOOP_LBUB_12_12_12(params->H6, armStateMPC_llbbyslb5, armStateMPC_lbIdx5, armStateMPC_lubbysub5, armStateMPC_ubIdx5, armStateMPC_Phi5);
armStateMPC_LA_DIAG_MATRIXFORWARDSUB_6_12(armStateMPC_Phi5, params->C6, armStateMPC_V5);
armStateMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_6_12(armStateMPC_Phi5, armStateMPC_D1, armStateMPC_W5);
armStateMPC_LA_DENSE_DIAGZERO_MMTM_6_12_6(armStateMPC_W5, armStateMPC_V5, armStateMPC_Ysd6);
armStateMPC_LA_DIAG_FORWARDSUB_12(armStateMPC_Phi5, armStateMPC_rd5, armStateMPC_Lbyrd5);
armStateMPC_LA_DIAG_CHOL_ONELOOP_LBUB_12_12_12(params->H7, armStateMPC_llbbyslb6, armStateMPC_lbIdx6, armStateMPC_lubbysub6, armStateMPC_ubIdx6, armStateMPC_Phi6);
armStateMPC_LA_DIAG_MATRIXFORWARDSUB_6_12(armStateMPC_Phi6, params->C7, armStateMPC_V6);
armStateMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_6_12(armStateMPC_Phi6, armStateMPC_D1, armStateMPC_W6);
armStateMPC_LA_DENSE_DIAGZERO_MMTM_6_12_6(armStateMPC_W6, armStateMPC_V6, armStateMPC_Ysd7);
armStateMPC_LA_DIAG_FORWARDSUB_12(armStateMPC_Phi6, armStateMPC_rd6, armStateMPC_Lbyrd6);
armStateMPC_LA_DIAG_CHOL_ONELOOP_LBUB_12_12_12(params->H8, armStateMPC_llbbyslb7, armStateMPC_lbIdx7, armStateMPC_lubbysub7, armStateMPC_ubIdx7, armStateMPC_Phi7);
armStateMPC_LA_DIAG_MATRIXFORWARDSUB_6_12(armStateMPC_Phi7, params->C8, armStateMPC_V7);
armStateMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_6_12(armStateMPC_Phi7, armStateMPC_D1, armStateMPC_W7);
armStateMPC_LA_DENSE_DIAGZERO_MMTM_6_12_6(armStateMPC_W7, armStateMPC_V7, armStateMPC_Ysd8);
armStateMPC_LA_DIAG_FORWARDSUB_12(armStateMPC_Phi7, armStateMPC_rd7, armStateMPC_Lbyrd7);
armStateMPC_LA_DIAG_CHOL_ONELOOP_LBUB_12_12_12(params->H9, armStateMPC_llbbyslb8, armStateMPC_lbIdx8, armStateMPC_lubbysub8, armStateMPC_ubIdx8, armStateMPC_Phi8);
armStateMPC_LA_DIAG_MATRIXFORWARDSUB_6_12(armStateMPC_Phi8, params->C9, armStateMPC_V8);
armStateMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_6_12(armStateMPC_Phi8, armStateMPC_D1, armStateMPC_W8);
armStateMPC_LA_DENSE_DIAGZERO_MMTM_6_12_6(armStateMPC_W8, armStateMPC_V8, armStateMPC_Ysd9);
armStateMPC_LA_DIAG_FORWARDSUB_12(armStateMPC_Phi8, armStateMPC_rd8, armStateMPC_Lbyrd8);
armStateMPC_LA_DIAG_CHOL_ONELOOP_LBUB_6_6_6(params->H10, armStateMPC_llbbyslb9, armStateMPC_lbIdx9, armStateMPC_lubbysub9, armStateMPC_ubIdx9, armStateMPC_Phi9);
armStateMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_6_6(armStateMPC_Phi9, armStateMPC_D9, armStateMPC_W9);
armStateMPC_LA_DIAG_FORWARDSUB_6(armStateMPC_Phi9, armStateMPC_rd9, armStateMPC_Lbyrd9);
armStateMPC_LA_DIAGZERO_MMT_6(armStateMPC_W0, armStateMPC_Yd0);
armStateMPC_LA_DIAGZERO_MVMSUB7_6(armStateMPC_W0, armStateMPC_Lbyrd0, armStateMPC_re0, armStateMPC_beta0);
armStateMPC_LA_DENSE_DIAGZERO_MMT2_6_12_12(armStateMPC_V0, armStateMPC_W1, armStateMPC_Yd1);
armStateMPC_LA_DENSE_DIAGZERO_2MVMSUB2_6_12_12(armStateMPC_V0, armStateMPC_Lbyrd0, armStateMPC_W1, armStateMPC_Lbyrd1, armStateMPC_re1, armStateMPC_beta1);
armStateMPC_LA_DENSE_DIAGZERO_MMT2_6_12_12(armStateMPC_V1, armStateMPC_W2, armStateMPC_Yd2);
armStateMPC_LA_DENSE_DIAGZERO_2MVMSUB2_6_12_12(armStateMPC_V1, armStateMPC_Lbyrd1, armStateMPC_W2, armStateMPC_Lbyrd2, armStateMPC_re2, armStateMPC_beta2);
armStateMPC_LA_DENSE_DIAGZERO_MMT2_6_12_12(armStateMPC_V2, armStateMPC_W3, armStateMPC_Yd3);
armStateMPC_LA_DENSE_DIAGZERO_2MVMSUB2_6_12_12(armStateMPC_V2, armStateMPC_Lbyrd2, armStateMPC_W3, armStateMPC_Lbyrd3, armStateMPC_re3, armStateMPC_beta3);
armStateMPC_LA_DENSE_DIAGZERO_MMT2_6_12_12(armStateMPC_V3, armStateMPC_W4, armStateMPC_Yd4);
armStateMPC_LA_DENSE_DIAGZERO_2MVMSUB2_6_12_12(armStateMPC_V3, armStateMPC_Lbyrd3, armStateMPC_W4, armStateMPC_Lbyrd4, armStateMPC_re4, armStateMPC_beta4);
armStateMPC_LA_DENSE_DIAGZERO_MMT2_6_12_12(armStateMPC_V4, armStateMPC_W5, armStateMPC_Yd5);
armStateMPC_LA_DENSE_DIAGZERO_2MVMSUB2_6_12_12(armStateMPC_V4, armStateMPC_Lbyrd4, armStateMPC_W5, armStateMPC_Lbyrd5, armStateMPC_re5, armStateMPC_beta5);
armStateMPC_LA_DENSE_DIAGZERO_MMT2_6_12_12(armStateMPC_V5, armStateMPC_W6, armStateMPC_Yd6);
armStateMPC_LA_DENSE_DIAGZERO_2MVMSUB2_6_12_12(armStateMPC_V5, armStateMPC_Lbyrd5, armStateMPC_W6, armStateMPC_Lbyrd6, armStateMPC_re6, armStateMPC_beta6);
armStateMPC_LA_DENSE_DIAGZERO_MMT2_6_12_12(armStateMPC_V6, armStateMPC_W7, armStateMPC_Yd7);
armStateMPC_LA_DENSE_DIAGZERO_2MVMSUB2_6_12_12(armStateMPC_V6, armStateMPC_Lbyrd6, armStateMPC_W7, armStateMPC_Lbyrd7, armStateMPC_re7, armStateMPC_beta7);
armStateMPC_LA_DENSE_DIAGZERO_MMT2_6_12_12(armStateMPC_V7, armStateMPC_W8, armStateMPC_Yd8);
armStateMPC_LA_DENSE_DIAGZERO_2MVMSUB2_6_12_12(armStateMPC_V7, armStateMPC_Lbyrd7, armStateMPC_W8, armStateMPC_Lbyrd8, armStateMPC_re8, armStateMPC_beta8);
armStateMPC_LA_DENSE_DIAGZERO_MMT2_6_12_6(armStateMPC_V8, armStateMPC_W9, armStateMPC_Yd9);
armStateMPC_LA_DENSE_DIAGZERO_2MVMSUB2_6_12_6(armStateMPC_V8, armStateMPC_Lbyrd8, armStateMPC_W9, armStateMPC_Lbyrd9, armStateMPC_re9, armStateMPC_beta9);
armStateMPC_LA_DENSE_CHOL_6(armStateMPC_Yd0, armStateMPC_Ld0);
armStateMPC_LA_DENSE_FORWARDSUB_6(armStateMPC_Ld0, armStateMPC_beta0, armStateMPC_yy0);
armStateMPC_LA_DENSE_MATRIXTFORWARDSUB_6_6(armStateMPC_Ld0, armStateMPC_Ysd1, armStateMPC_Lsd1);
armStateMPC_LA_DENSE_MMTSUB_6_6(armStateMPC_Lsd1, armStateMPC_Yd1);
armStateMPC_LA_DENSE_CHOL_6(armStateMPC_Yd1, armStateMPC_Ld1);
armStateMPC_LA_DENSE_MVMSUB1_6_6(armStateMPC_Lsd1, armStateMPC_yy0, armStateMPC_beta1, armStateMPC_bmy1);
armStateMPC_LA_DENSE_FORWARDSUB_6(armStateMPC_Ld1, armStateMPC_bmy1, armStateMPC_yy1);
armStateMPC_LA_DENSE_MATRIXTFORWARDSUB_6_6(armStateMPC_Ld1, armStateMPC_Ysd2, armStateMPC_Lsd2);
armStateMPC_LA_DENSE_MMTSUB_6_6(armStateMPC_Lsd2, armStateMPC_Yd2);
armStateMPC_LA_DENSE_CHOL_6(armStateMPC_Yd2, armStateMPC_Ld2);
armStateMPC_LA_DENSE_MVMSUB1_6_6(armStateMPC_Lsd2, armStateMPC_yy1, armStateMPC_beta2, armStateMPC_bmy2);
armStateMPC_LA_DENSE_FORWARDSUB_6(armStateMPC_Ld2, armStateMPC_bmy2, armStateMPC_yy2);
armStateMPC_LA_DENSE_MATRIXTFORWARDSUB_6_6(armStateMPC_Ld2, armStateMPC_Ysd3, armStateMPC_Lsd3);
armStateMPC_LA_DENSE_MMTSUB_6_6(armStateMPC_Lsd3, armStateMPC_Yd3);
armStateMPC_LA_DENSE_CHOL_6(armStateMPC_Yd3, armStateMPC_Ld3);
armStateMPC_LA_DENSE_MVMSUB1_6_6(armStateMPC_Lsd3, armStateMPC_yy2, armStateMPC_beta3, armStateMPC_bmy3);
armStateMPC_LA_DENSE_FORWARDSUB_6(armStateMPC_Ld3, armStateMPC_bmy3, armStateMPC_yy3);
armStateMPC_LA_DENSE_MATRIXTFORWARDSUB_6_6(armStateMPC_Ld3, armStateMPC_Ysd4, armStateMPC_Lsd4);
armStateMPC_LA_DENSE_MMTSUB_6_6(armStateMPC_Lsd4, armStateMPC_Yd4);
armStateMPC_LA_DENSE_CHOL_6(armStateMPC_Yd4, armStateMPC_Ld4);
armStateMPC_LA_DENSE_MVMSUB1_6_6(armStateMPC_Lsd4, armStateMPC_yy3, armStateMPC_beta4, armStateMPC_bmy4);
armStateMPC_LA_DENSE_FORWARDSUB_6(armStateMPC_Ld4, armStateMPC_bmy4, armStateMPC_yy4);
armStateMPC_LA_DENSE_MATRIXTFORWARDSUB_6_6(armStateMPC_Ld4, armStateMPC_Ysd5, armStateMPC_Lsd5);
armStateMPC_LA_DENSE_MMTSUB_6_6(armStateMPC_Lsd5, armStateMPC_Yd5);
armStateMPC_LA_DENSE_CHOL_6(armStateMPC_Yd5, armStateMPC_Ld5);
armStateMPC_LA_DENSE_MVMSUB1_6_6(armStateMPC_Lsd5, armStateMPC_yy4, armStateMPC_beta5, armStateMPC_bmy5);
armStateMPC_LA_DENSE_FORWARDSUB_6(armStateMPC_Ld5, armStateMPC_bmy5, armStateMPC_yy5);
armStateMPC_LA_DENSE_MATRIXTFORWARDSUB_6_6(armStateMPC_Ld5, armStateMPC_Ysd6, armStateMPC_Lsd6);
armStateMPC_LA_DENSE_MMTSUB_6_6(armStateMPC_Lsd6, armStateMPC_Yd6);
armStateMPC_LA_DENSE_CHOL_6(armStateMPC_Yd6, armStateMPC_Ld6);
armStateMPC_LA_DENSE_MVMSUB1_6_6(armStateMPC_Lsd6, armStateMPC_yy5, armStateMPC_beta6, armStateMPC_bmy6);
armStateMPC_LA_DENSE_FORWARDSUB_6(armStateMPC_Ld6, armStateMPC_bmy6, armStateMPC_yy6);
armStateMPC_LA_DENSE_MATRIXTFORWARDSUB_6_6(armStateMPC_Ld6, armStateMPC_Ysd7, armStateMPC_Lsd7);
armStateMPC_LA_DENSE_MMTSUB_6_6(armStateMPC_Lsd7, armStateMPC_Yd7);
armStateMPC_LA_DENSE_CHOL_6(armStateMPC_Yd7, armStateMPC_Ld7);
armStateMPC_LA_DENSE_MVMSUB1_6_6(armStateMPC_Lsd7, armStateMPC_yy6, armStateMPC_beta7, armStateMPC_bmy7);
armStateMPC_LA_DENSE_FORWARDSUB_6(armStateMPC_Ld7, armStateMPC_bmy7, armStateMPC_yy7);
armStateMPC_LA_DENSE_MATRIXTFORWARDSUB_6_6(armStateMPC_Ld7, armStateMPC_Ysd8, armStateMPC_Lsd8);
armStateMPC_LA_DENSE_MMTSUB_6_6(armStateMPC_Lsd8, armStateMPC_Yd8);
armStateMPC_LA_DENSE_CHOL_6(armStateMPC_Yd8, armStateMPC_Ld8);
armStateMPC_LA_DENSE_MVMSUB1_6_6(armStateMPC_Lsd8, armStateMPC_yy7, armStateMPC_beta8, armStateMPC_bmy8);
armStateMPC_LA_DENSE_FORWARDSUB_6(armStateMPC_Ld8, armStateMPC_bmy8, armStateMPC_yy8);
armStateMPC_LA_DENSE_MATRIXTFORWARDSUB_6_6(armStateMPC_Ld8, armStateMPC_Ysd9, armStateMPC_Lsd9);
armStateMPC_LA_DENSE_MMTSUB_6_6(armStateMPC_Lsd9, armStateMPC_Yd9);
armStateMPC_LA_DENSE_CHOL_6(armStateMPC_Yd9, armStateMPC_Ld9);
armStateMPC_LA_DENSE_MVMSUB1_6_6(armStateMPC_Lsd9, armStateMPC_yy8, armStateMPC_beta9, armStateMPC_bmy9);
armStateMPC_LA_DENSE_FORWARDSUB_6(armStateMPC_Ld9, armStateMPC_bmy9, armStateMPC_yy9);
armStateMPC_LA_DENSE_BACKWARDSUB_6(armStateMPC_Ld9, armStateMPC_yy9, armStateMPC_dvaff9);
armStateMPC_LA_DENSE_MTVMSUB_6_6(armStateMPC_Lsd9, armStateMPC_dvaff9, armStateMPC_yy8, armStateMPC_bmy8);
armStateMPC_LA_DENSE_BACKWARDSUB_6(armStateMPC_Ld8, armStateMPC_bmy8, armStateMPC_dvaff8);
armStateMPC_LA_DENSE_MTVMSUB_6_6(armStateMPC_Lsd8, armStateMPC_dvaff8, armStateMPC_yy7, armStateMPC_bmy7);
armStateMPC_LA_DENSE_BACKWARDSUB_6(armStateMPC_Ld7, armStateMPC_bmy7, armStateMPC_dvaff7);
armStateMPC_LA_DENSE_MTVMSUB_6_6(armStateMPC_Lsd7, armStateMPC_dvaff7, armStateMPC_yy6, armStateMPC_bmy6);
armStateMPC_LA_DENSE_BACKWARDSUB_6(armStateMPC_Ld6, armStateMPC_bmy6, armStateMPC_dvaff6);
armStateMPC_LA_DENSE_MTVMSUB_6_6(armStateMPC_Lsd6, armStateMPC_dvaff6, armStateMPC_yy5, armStateMPC_bmy5);
armStateMPC_LA_DENSE_BACKWARDSUB_6(armStateMPC_Ld5, armStateMPC_bmy5, armStateMPC_dvaff5);
armStateMPC_LA_DENSE_MTVMSUB_6_6(armStateMPC_Lsd5, armStateMPC_dvaff5, armStateMPC_yy4, armStateMPC_bmy4);
armStateMPC_LA_DENSE_BACKWARDSUB_6(armStateMPC_Ld4, armStateMPC_bmy4, armStateMPC_dvaff4);
armStateMPC_LA_DENSE_MTVMSUB_6_6(armStateMPC_Lsd4, armStateMPC_dvaff4, armStateMPC_yy3, armStateMPC_bmy3);
armStateMPC_LA_DENSE_BACKWARDSUB_6(armStateMPC_Ld3, armStateMPC_bmy3, armStateMPC_dvaff3);
armStateMPC_LA_DENSE_MTVMSUB_6_6(armStateMPC_Lsd3, armStateMPC_dvaff3, armStateMPC_yy2, armStateMPC_bmy2);
armStateMPC_LA_DENSE_BACKWARDSUB_6(armStateMPC_Ld2, armStateMPC_bmy2, armStateMPC_dvaff2);
armStateMPC_LA_DENSE_MTVMSUB_6_6(armStateMPC_Lsd2, armStateMPC_dvaff2, armStateMPC_yy1, armStateMPC_bmy1);
armStateMPC_LA_DENSE_BACKWARDSUB_6(armStateMPC_Ld1, armStateMPC_bmy1, armStateMPC_dvaff1);
armStateMPC_LA_DENSE_MTVMSUB_6_6(armStateMPC_Lsd1, armStateMPC_dvaff1, armStateMPC_yy0, armStateMPC_bmy0);
armStateMPC_LA_DENSE_BACKWARDSUB_6(armStateMPC_Ld0, armStateMPC_bmy0, armStateMPC_dvaff0);
armStateMPC_LA_DENSE_DIAGZERO_MTVM2_6_12_6(params->C1, armStateMPC_dvaff1, armStateMPC_D0, armStateMPC_dvaff0, armStateMPC_grad_eq0);
armStateMPC_LA_DENSE_DIAGZERO_MTVM2_6_12_6(params->C2, armStateMPC_dvaff2, armStateMPC_D1, armStateMPC_dvaff1, armStateMPC_grad_eq1);
armStateMPC_LA_DENSE_DIAGZERO_MTVM2_6_12_6(params->C3, armStateMPC_dvaff3, armStateMPC_D1, armStateMPC_dvaff2, armStateMPC_grad_eq2);
armStateMPC_LA_DENSE_DIAGZERO_MTVM2_6_12_6(params->C4, armStateMPC_dvaff4, armStateMPC_D1, armStateMPC_dvaff3, armStateMPC_grad_eq3);
armStateMPC_LA_DENSE_DIAGZERO_MTVM2_6_12_6(params->C5, armStateMPC_dvaff5, armStateMPC_D1, armStateMPC_dvaff4, armStateMPC_grad_eq4);
armStateMPC_LA_DENSE_DIAGZERO_MTVM2_6_12_6(params->C6, armStateMPC_dvaff6, armStateMPC_D1, armStateMPC_dvaff5, armStateMPC_grad_eq5);
armStateMPC_LA_DENSE_DIAGZERO_MTVM2_6_12_6(params->C7, armStateMPC_dvaff7, armStateMPC_D1, armStateMPC_dvaff6, armStateMPC_grad_eq6);
armStateMPC_LA_DENSE_DIAGZERO_MTVM2_6_12_6(params->C8, armStateMPC_dvaff8, armStateMPC_D1, armStateMPC_dvaff7, armStateMPC_grad_eq7);
armStateMPC_LA_DENSE_DIAGZERO_MTVM2_6_12_6(params->C9, armStateMPC_dvaff9, armStateMPC_D1, armStateMPC_dvaff8, armStateMPC_grad_eq8);
armStateMPC_LA_DIAGZERO_MTVM_6_6(armStateMPC_D9, armStateMPC_dvaff9, armStateMPC_grad_eq9);
armStateMPC_LA_VSUB2_114(armStateMPC_rd, armStateMPC_grad_eq, armStateMPC_rd);
armStateMPC_LA_DIAG_FORWARDBACKWARDSUB_12(armStateMPC_Phi0, armStateMPC_rd0, armStateMPC_dzaff0);
armStateMPC_LA_DIAG_FORWARDBACKWARDSUB_12(armStateMPC_Phi1, armStateMPC_rd1, armStateMPC_dzaff1);
armStateMPC_LA_DIAG_FORWARDBACKWARDSUB_12(armStateMPC_Phi2, armStateMPC_rd2, armStateMPC_dzaff2);
armStateMPC_LA_DIAG_FORWARDBACKWARDSUB_12(armStateMPC_Phi3, armStateMPC_rd3, armStateMPC_dzaff3);
armStateMPC_LA_DIAG_FORWARDBACKWARDSUB_12(armStateMPC_Phi4, armStateMPC_rd4, armStateMPC_dzaff4);
armStateMPC_LA_DIAG_FORWARDBACKWARDSUB_12(armStateMPC_Phi5, armStateMPC_rd5, armStateMPC_dzaff5);
armStateMPC_LA_DIAG_FORWARDBACKWARDSUB_12(armStateMPC_Phi6, armStateMPC_rd6, armStateMPC_dzaff6);
armStateMPC_LA_DIAG_FORWARDBACKWARDSUB_12(armStateMPC_Phi7, armStateMPC_rd7, armStateMPC_dzaff7);
armStateMPC_LA_DIAG_FORWARDBACKWARDSUB_12(armStateMPC_Phi8, armStateMPC_rd8, armStateMPC_dzaff8);
armStateMPC_LA_DIAG_FORWARDBACKWARDSUB_6(armStateMPC_Phi9, armStateMPC_rd9, armStateMPC_dzaff9);
armStateMPC_LA_VSUB_INDEXED_12(armStateMPC_dzaff0, armStateMPC_lbIdx0, armStateMPC_rilb0, armStateMPC_dslbaff0);
armStateMPC_LA_VSUB3_12(armStateMPC_llbbyslb0, armStateMPC_dslbaff0, armStateMPC_llb0, armStateMPC_dllbaff0);
armStateMPC_LA_VSUB2_INDEXED_12(armStateMPC_riub0, armStateMPC_dzaff0, armStateMPC_ubIdx0, armStateMPC_dsubaff0);
armStateMPC_LA_VSUB3_12(armStateMPC_lubbysub0, armStateMPC_dsubaff0, armStateMPC_lub0, armStateMPC_dlubaff0);
armStateMPC_LA_VSUB_INDEXED_12(armStateMPC_dzaff1, armStateMPC_lbIdx1, armStateMPC_rilb1, armStateMPC_dslbaff1);
armStateMPC_LA_VSUB3_12(armStateMPC_llbbyslb1, armStateMPC_dslbaff1, armStateMPC_llb1, armStateMPC_dllbaff1);
armStateMPC_LA_VSUB2_INDEXED_12(armStateMPC_riub1, armStateMPC_dzaff1, armStateMPC_ubIdx1, armStateMPC_dsubaff1);
armStateMPC_LA_VSUB3_12(armStateMPC_lubbysub1, armStateMPC_dsubaff1, armStateMPC_lub1, armStateMPC_dlubaff1);
armStateMPC_LA_VSUB_INDEXED_12(armStateMPC_dzaff2, armStateMPC_lbIdx2, armStateMPC_rilb2, armStateMPC_dslbaff2);
armStateMPC_LA_VSUB3_12(armStateMPC_llbbyslb2, armStateMPC_dslbaff2, armStateMPC_llb2, armStateMPC_dllbaff2);
armStateMPC_LA_VSUB2_INDEXED_12(armStateMPC_riub2, armStateMPC_dzaff2, armStateMPC_ubIdx2, armStateMPC_dsubaff2);
armStateMPC_LA_VSUB3_12(armStateMPC_lubbysub2, armStateMPC_dsubaff2, armStateMPC_lub2, armStateMPC_dlubaff2);
armStateMPC_LA_VSUB_INDEXED_12(armStateMPC_dzaff3, armStateMPC_lbIdx3, armStateMPC_rilb3, armStateMPC_dslbaff3);
armStateMPC_LA_VSUB3_12(armStateMPC_llbbyslb3, armStateMPC_dslbaff3, armStateMPC_llb3, armStateMPC_dllbaff3);
armStateMPC_LA_VSUB2_INDEXED_12(armStateMPC_riub3, armStateMPC_dzaff3, armStateMPC_ubIdx3, armStateMPC_dsubaff3);
armStateMPC_LA_VSUB3_12(armStateMPC_lubbysub3, armStateMPC_dsubaff3, armStateMPC_lub3, armStateMPC_dlubaff3);
armStateMPC_LA_VSUB_INDEXED_12(armStateMPC_dzaff4, armStateMPC_lbIdx4, armStateMPC_rilb4, armStateMPC_dslbaff4);
armStateMPC_LA_VSUB3_12(armStateMPC_llbbyslb4, armStateMPC_dslbaff4, armStateMPC_llb4, armStateMPC_dllbaff4);
armStateMPC_LA_VSUB2_INDEXED_12(armStateMPC_riub4, armStateMPC_dzaff4, armStateMPC_ubIdx4, armStateMPC_dsubaff4);
armStateMPC_LA_VSUB3_12(armStateMPC_lubbysub4, armStateMPC_dsubaff4, armStateMPC_lub4, armStateMPC_dlubaff4);
armStateMPC_LA_VSUB_INDEXED_12(armStateMPC_dzaff5, armStateMPC_lbIdx5, armStateMPC_rilb5, armStateMPC_dslbaff5);
armStateMPC_LA_VSUB3_12(armStateMPC_llbbyslb5, armStateMPC_dslbaff5, armStateMPC_llb5, armStateMPC_dllbaff5);
armStateMPC_LA_VSUB2_INDEXED_12(armStateMPC_riub5, armStateMPC_dzaff5, armStateMPC_ubIdx5, armStateMPC_dsubaff5);
armStateMPC_LA_VSUB3_12(armStateMPC_lubbysub5, armStateMPC_dsubaff5, armStateMPC_lub5, armStateMPC_dlubaff5);
armStateMPC_LA_VSUB_INDEXED_12(armStateMPC_dzaff6, armStateMPC_lbIdx6, armStateMPC_rilb6, armStateMPC_dslbaff6);
armStateMPC_LA_VSUB3_12(armStateMPC_llbbyslb6, armStateMPC_dslbaff6, armStateMPC_llb6, armStateMPC_dllbaff6);
armStateMPC_LA_VSUB2_INDEXED_12(armStateMPC_riub6, armStateMPC_dzaff6, armStateMPC_ubIdx6, armStateMPC_dsubaff6);
armStateMPC_LA_VSUB3_12(armStateMPC_lubbysub6, armStateMPC_dsubaff6, armStateMPC_lub6, armStateMPC_dlubaff6);
armStateMPC_LA_VSUB_INDEXED_12(armStateMPC_dzaff7, armStateMPC_lbIdx7, armStateMPC_rilb7, armStateMPC_dslbaff7);
armStateMPC_LA_VSUB3_12(armStateMPC_llbbyslb7, armStateMPC_dslbaff7, armStateMPC_llb7, armStateMPC_dllbaff7);
armStateMPC_LA_VSUB2_INDEXED_12(armStateMPC_riub7, armStateMPC_dzaff7, armStateMPC_ubIdx7, armStateMPC_dsubaff7);
armStateMPC_LA_VSUB3_12(armStateMPC_lubbysub7, armStateMPC_dsubaff7, armStateMPC_lub7, armStateMPC_dlubaff7);
armStateMPC_LA_VSUB_INDEXED_12(armStateMPC_dzaff8, armStateMPC_lbIdx8, armStateMPC_rilb8, armStateMPC_dslbaff8);
armStateMPC_LA_VSUB3_12(armStateMPC_llbbyslb8, armStateMPC_dslbaff8, armStateMPC_llb8, armStateMPC_dllbaff8);
armStateMPC_LA_VSUB2_INDEXED_12(armStateMPC_riub8, armStateMPC_dzaff8, armStateMPC_ubIdx8, armStateMPC_dsubaff8);
armStateMPC_LA_VSUB3_12(armStateMPC_lubbysub8, armStateMPC_dsubaff8, armStateMPC_lub8, armStateMPC_dlubaff8);
armStateMPC_LA_VSUB_INDEXED_6(armStateMPC_dzaff9, armStateMPC_lbIdx9, armStateMPC_rilb9, armStateMPC_dslbaff9);
armStateMPC_LA_VSUB3_6(armStateMPC_llbbyslb9, armStateMPC_dslbaff9, armStateMPC_llb9, armStateMPC_dllbaff9);
armStateMPC_LA_VSUB2_INDEXED_6(armStateMPC_riub9, armStateMPC_dzaff9, armStateMPC_ubIdx9, armStateMPC_dsubaff9);
armStateMPC_LA_VSUB3_6(armStateMPC_lubbysub9, armStateMPC_dsubaff9, armStateMPC_lub9, armStateMPC_dlubaff9);
info->lsit_aff = armStateMPC_LINESEARCH_BACKTRACKING_AFFINE(armStateMPC_l, armStateMPC_s, armStateMPC_dl_aff, armStateMPC_ds_aff, &info->step_aff, &info->mu_aff);
if( info->lsit_aff == armStateMPC_NOPROGRESS ){
exitcode = armStateMPC_NOPROGRESS; break;
}
sigma_3rdroot = info->mu_aff / info->mu;
info->sigma = sigma_3rdroot*sigma_3rdroot*sigma_3rdroot;
musigma = info->mu * info->sigma;
armStateMPC_LA_VSUB5_228(armStateMPC_ds_aff, armStateMPC_dl_aff, musigma, armStateMPC_ccrhs);
armStateMPC_LA_VSUB6_INDEXED_12_12_12(armStateMPC_ccrhsub0, armStateMPC_sub0, armStateMPC_ubIdx0, armStateMPC_ccrhsl0, armStateMPC_slb0, armStateMPC_lbIdx0, armStateMPC_rd0);
armStateMPC_LA_VSUB6_INDEXED_12_12_12(armStateMPC_ccrhsub1, armStateMPC_sub1, armStateMPC_ubIdx1, armStateMPC_ccrhsl1, armStateMPC_slb1, armStateMPC_lbIdx1, armStateMPC_rd1);
armStateMPC_LA_DIAG_FORWARDSUB_12(armStateMPC_Phi0, armStateMPC_rd0, armStateMPC_Lbyrd0);
armStateMPC_LA_DIAG_FORWARDSUB_12(armStateMPC_Phi1, armStateMPC_rd1, armStateMPC_Lbyrd1);
armStateMPC_LA_DIAGZERO_MVM_6(armStateMPC_W0, armStateMPC_Lbyrd0, armStateMPC_beta0);
armStateMPC_LA_DENSE_FORWARDSUB_6(armStateMPC_Ld0, armStateMPC_beta0, armStateMPC_yy0);
armStateMPC_LA_DENSE_DIAGZERO_2MVMADD_6_12_12(armStateMPC_V0, armStateMPC_Lbyrd0, armStateMPC_W1, armStateMPC_Lbyrd1, armStateMPC_beta1);
armStateMPC_LA_DENSE_MVMSUB1_6_6(armStateMPC_Lsd1, armStateMPC_yy0, armStateMPC_beta1, armStateMPC_bmy1);
armStateMPC_LA_DENSE_FORWARDSUB_6(armStateMPC_Ld1, armStateMPC_bmy1, armStateMPC_yy1);
armStateMPC_LA_VSUB6_INDEXED_12_12_12(armStateMPC_ccrhsub2, armStateMPC_sub2, armStateMPC_ubIdx2, armStateMPC_ccrhsl2, armStateMPC_slb2, armStateMPC_lbIdx2, armStateMPC_rd2);
armStateMPC_LA_DIAG_FORWARDSUB_12(armStateMPC_Phi2, armStateMPC_rd2, armStateMPC_Lbyrd2);
armStateMPC_LA_DENSE_DIAGZERO_2MVMADD_6_12_12(armStateMPC_V1, armStateMPC_Lbyrd1, armStateMPC_W2, armStateMPC_Lbyrd2, armStateMPC_beta2);
armStateMPC_LA_DENSE_MVMSUB1_6_6(armStateMPC_Lsd2, armStateMPC_yy1, armStateMPC_beta2, armStateMPC_bmy2);
armStateMPC_LA_DENSE_FORWARDSUB_6(armStateMPC_Ld2, armStateMPC_bmy2, armStateMPC_yy2);
armStateMPC_LA_VSUB6_INDEXED_12_12_12(armStateMPC_ccrhsub3, armStateMPC_sub3, armStateMPC_ubIdx3, armStateMPC_ccrhsl3, armStateMPC_slb3, armStateMPC_lbIdx3, armStateMPC_rd3);
armStateMPC_LA_DIAG_FORWARDSUB_12(armStateMPC_Phi3, armStateMPC_rd3, armStateMPC_Lbyrd3);
armStateMPC_LA_DENSE_DIAGZERO_2MVMADD_6_12_12(armStateMPC_V2, armStateMPC_Lbyrd2, armStateMPC_W3, armStateMPC_Lbyrd3, armStateMPC_beta3);
armStateMPC_LA_DENSE_MVMSUB1_6_6(armStateMPC_Lsd3, armStateMPC_yy2, armStateMPC_beta3, armStateMPC_bmy3);
armStateMPC_LA_DENSE_FORWARDSUB_6(armStateMPC_Ld3, armStateMPC_bmy3, armStateMPC_yy3);
armStateMPC_LA_VSUB6_INDEXED_12_12_12(armStateMPC_ccrhsub4, armStateMPC_sub4, armStateMPC_ubIdx4, armStateMPC_ccrhsl4, armStateMPC_slb4, armStateMPC_lbIdx4, armStateMPC_rd4);
armStateMPC_LA_DIAG_FORWARDSUB_12(armStateMPC_Phi4, armStateMPC_rd4, armStateMPC_Lbyrd4);
armStateMPC_LA_DENSE_DIAGZERO_2MVMADD_6_12_12(armStateMPC_V3, armStateMPC_Lbyrd3, armStateMPC_W4, armStateMPC_Lbyrd4, armStateMPC_beta4);
armStateMPC_LA_DENSE_MVMSUB1_6_6(armStateMPC_Lsd4, armStateMPC_yy3, armStateMPC_beta4, armStateMPC_bmy4);
armStateMPC_LA_DENSE_FORWARDSUB_6(armStateMPC_Ld4, armStateMPC_bmy4, armStateMPC_yy4);
armStateMPC_LA_VSUB6_INDEXED_12_12_12(armStateMPC_ccrhsub5, armStateMPC_sub5, armStateMPC_ubIdx5, armStateMPC_ccrhsl5, armStateMPC_slb5, armStateMPC_lbIdx5, armStateMPC_rd5);
armStateMPC_LA_DIAG_FORWARDSUB_12(armStateMPC_Phi5, armStateMPC_rd5, armStateMPC_Lbyrd5);
armStateMPC_LA_DENSE_DIAGZERO_2MVMADD_6_12_12(armStateMPC_V4, armStateMPC_Lbyrd4, armStateMPC_W5, armStateMPC_Lbyrd5, armStateMPC_beta5);
armStateMPC_LA_DENSE_MVMSUB1_6_6(armStateMPC_Lsd5, armStateMPC_yy4, armStateMPC_beta5, armStateMPC_bmy5);
armStateMPC_LA_DENSE_FORWARDSUB_6(armStateMPC_Ld5, armStateMPC_bmy5, armStateMPC_yy5);
armStateMPC_LA_VSUB6_INDEXED_12_12_12(armStateMPC_ccrhsub6, armStateMPC_sub6, armStateMPC_ubIdx6, armStateMPC_ccrhsl6, armStateMPC_slb6, armStateMPC_lbIdx6, armStateMPC_rd6);
armStateMPC_LA_DIAG_FORWARDSUB_12(armStateMPC_Phi6, armStateMPC_rd6, armStateMPC_Lbyrd6);
armStateMPC_LA_DENSE_DIAGZERO_2MVMADD_6_12_12(armStateMPC_V5, armStateMPC_Lbyrd5, armStateMPC_W6, armStateMPC_Lbyrd6, armStateMPC_beta6);
armStateMPC_LA_DENSE_MVMSUB1_6_6(armStateMPC_Lsd6, armStateMPC_yy5, armStateMPC_beta6, armStateMPC_bmy6);
armStateMPC_LA_DENSE_FORWARDSUB_6(armStateMPC_Ld6, armStateMPC_bmy6, armStateMPC_yy6);
armStateMPC_LA_VSUB6_INDEXED_12_12_12(armStateMPC_ccrhsub7, armStateMPC_sub7, armStateMPC_ubIdx7, armStateMPC_ccrhsl7, armStateMPC_slb7, armStateMPC_lbIdx7, armStateMPC_rd7);
armStateMPC_LA_DIAG_FORWARDSUB_12(armStateMPC_Phi7, armStateMPC_rd7, armStateMPC_Lbyrd7);
armStateMPC_LA_DENSE_DIAGZERO_2MVMADD_6_12_12(armStateMPC_V6, armStateMPC_Lbyrd6, armStateMPC_W7, armStateMPC_Lbyrd7, armStateMPC_beta7);
armStateMPC_LA_DENSE_MVMSUB1_6_6(armStateMPC_Lsd7, armStateMPC_yy6, armStateMPC_beta7, armStateMPC_bmy7);
armStateMPC_LA_DENSE_FORWARDSUB_6(armStateMPC_Ld7, armStateMPC_bmy7, armStateMPC_yy7);
armStateMPC_LA_VSUB6_INDEXED_12_12_12(armStateMPC_ccrhsub8, armStateMPC_sub8, armStateMPC_ubIdx8, armStateMPC_ccrhsl8, armStateMPC_slb8, armStateMPC_lbIdx8, armStateMPC_rd8);
armStateMPC_LA_DIAG_FORWARDSUB_12(armStateMPC_Phi8, armStateMPC_rd8, armStateMPC_Lbyrd8);
armStateMPC_LA_DENSE_DIAGZERO_2MVMADD_6_12_12(armStateMPC_V7, armStateMPC_Lbyrd7, armStateMPC_W8, armStateMPC_Lbyrd8, armStateMPC_beta8);
armStateMPC_LA_DENSE_MVMSUB1_6_6(armStateMPC_Lsd8, armStateMPC_yy7, armStateMPC_beta8, armStateMPC_bmy8);
armStateMPC_LA_DENSE_FORWARDSUB_6(armStateMPC_Ld8, armStateMPC_bmy8, armStateMPC_yy8);
armStateMPC_LA_VSUB6_INDEXED_6_6_6(armStateMPC_ccrhsub9, armStateMPC_sub9, armStateMPC_ubIdx9, armStateMPC_ccrhsl9, armStateMPC_slb9, armStateMPC_lbIdx9, armStateMPC_rd9);
armStateMPC_LA_DIAG_FORWARDSUB_6(armStateMPC_Phi9, armStateMPC_rd9, armStateMPC_Lbyrd9);
armStateMPC_LA_DENSE_DIAGZERO_2MVMADD_6_12_6(armStateMPC_V8, armStateMPC_Lbyrd8, armStateMPC_W9, armStateMPC_Lbyrd9, armStateMPC_beta9);
armStateMPC_LA_DENSE_MVMSUB1_6_6(armStateMPC_Lsd9, armStateMPC_yy8, armStateMPC_beta9, armStateMPC_bmy9);
armStateMPC_LA_DENSE_FORWARDSUB_6(armStateMPC_Ld9, armStateMPC_bmy9, armStateMPC_yy9);
armStateMPC_LA_DENSE_BACKWARDSUB_6(armStateMPC_Ld9, armStateMPC_yy9, armStateMPC_dvcc9);
armStateMPC_LA_DENSE_MTVMSUB_6_6(armStateMPC_Lsd9, armStateMPC_dvcc9, armStateMPC_yy8, armStateMPC_bmy8);
armStateMPC_LA_DENSE_BACKWARDSUB_6(armStateMPC_Ld8, armStateMPC_bmy8, armStateMPC_dvcc8);
armStateMPC_LA_DENSE_MTVMSUB_6_6(armStateMPC_Lsd8, armStateMPC_dvcc8, armStateMPC_yy7, armStateMPC_bmy7);
armStateMPC_LA_DENSE_BACKWARDSUB_6(armStateMPC_Ld7, armStateMPC_bmy7, armStateMPC_dvcc7);
armStateMPC_LA_DENSE_MTVMSUB_6_6(armStateMPC_Lsd7, armStateMPC_dvcc7, armStateMPC_yy6, armStateMPC_bmy6);
armStateMPC_LA_DENSE_BACKWARDSUB_6(armStateMPC_Ld6, armStateMPC_bmy6, armStateMPC_dvcc6);
armStateMPC_LA_DENSE_MTVMSUB_6_6(armStateMPC_Lsd6, armStateMPC_dvcc6, armStateMPC_yy5, armStateMPC_bmy5);
armStateMPC_LA_DENSE_BACKWARDSUB_6(armStateMPC_Ld5, armStateMPC_bmy5, armStateMPC_dvcc5);
armStateMPC_LA_DENSE_MTVMSUB_6_6(armStateMPC_Lsd5, armStateMPC_dvcc5, armStateMPC_yy4, armStateMPC_bmy4);
armStateMPC_LA_DENSE_BACKWARDSUB_6(armStateMPC_Ld4, armStateMPC_bmy4, armStateMPC_dvcc4);
armStateMPC_LA_DENSE_MTVMSUB_6_6(armStateMPC_Lsd4, armStateMPC_dvcc4, armStateMPC_yy3, armStateMPC_bmy3);
armStateMPC_LA_DENSE_BACKWARDSUB_6(armStateMPC_Ld3, armStateMPC_bmy3, armStateMPC_dvcc3);
armStateMPC_LA_DENSE_MTVMSUB_6_6(armStateMPC_Lsd3, armStateMPC_dvcc3, armStateMPC_yy2, armStateMPC_bmy2);
armStateMPC_LA_DENSE_BACKWARDSUB_6(armStateMPC_Ld2, armStateMPC_bmy2, armStateMPC_dvcc2);
armStateMPC_LA_DENSE_MTVMSUB_6_6(armStateMPC_Lsd2, armStateMPC_dvcc2, armStateMPC_yy1, armStateMPC_bmy1);
armStateMPC_LA_DENSE_BACKWARDSUB_6(armStateMPC_Ld1, armStateMPC_bmy1, armStateMPC_dvcc1);
armStateMPC_LA_DENSE_MTVMSUB_6_6(armStateMPC_Lsd1, armStateMPC_dvcc1, armStateMPC_yy0, armStateMPC_bmy0);
armStateMPC_LA_DENSE_BACKWARDSUB_6(armStateMPC_Ld0, armStateMPC_bmy0, armStateMPC_dvcc0);
armStateMPC_LA_DENSE_DIAGZERO_MTVM2_6_12_6(params->C1, armStateMPC_dvcc1, armStateMPC_D0, armStateMPC_dvcc0, armStateMPC_grad_eq0);
armStateMPC_LA_DENSE_DIAGZERO_MTVM2_6_12_6(params->C2, armStateMPC_dvcc2, armStateMPC_D1, armStateMPC_dvcc1, armStateMPC_grad_eq1);
armStateMPC_LA_DENSE_DIAGZERO_MTVM2_6_12_6(params->C3, armStateMPC_dvcc3, armStateMPC_D1, armStateMPC_dvcc2, armStateMPC_grad_eq2);
armStateMPC_LA_DENSE_DIAGZERO_MTVM2_6_12_6(params->C4, armStateMPC_dvcc4, armStateMPC_D1, armStateMPC_dvcc3, armStateMPC_grad_eq3);
armStateMPC_LA_DENSE_DIAGZERO_MTVM2_6_12_6(params->C5, armStateMPC_dvcc5, armStateMPC_D1, armStateMPC_dvcc4, armStateMPC_grad_eq4);
armStateMPC_LA_DENSE_DIAGZERO_MTVM2_6_12_6(params->C6, armStateMPC_dvcc6, armStateMPC_D1, armStateMPC_dvcc5, armStateMPC_grad_eq5);
armStateMPC_LA_DENSE_DIAGZERO_MTVM2_6_12_6(params->C7, armStateMPC_dvcc7, armStateMPC_D1, armStateMPC_dvcc6, armStateMPC_grad_eq6);
armStateMPC_LA_DENSE_DIAGZERO_MTVM2_6_12_6(params->C8, armStateMPC_dvcc8, armStateMPC_D1, armStateMPC_dvcc7, armStateMPC_grad_eq7);
armStateMPC_LA_DENSE_DIAGZERO_MTVM2_6_12_6(params->C9, armStateMPC_dvcc9, armStateMPC_D1, armStateMPC_dvcc8, armStateMPC_grad_eq8);
armStateMPC_LA_DIAGZERO_MTVM_6_6(armStateMPC_D9, armStateMPC_dvcc9, armStateMPC_grad_eq9);
armStateMPC_LA_VSUB_114(armStateMPC_rd, armStateMPC_grad_eq, armStateMPC_rd);
armStateMPC_LA_DIAG_FORWARDBACKWARDSUB_12(armStateMPC_Phi0, armStateMPC_rd0, armStateMPC_dzcc0);
armStateMPC_LA_DIAG_FORWARDBACKWARDSUB_12(armStateMPC_Phi1, armStateMPC_rd1, armStateMPC_dzcc1);
armStateMPC_LA_DIAG_FORWARDBACKWARDSUB_12(armStateMPC_Phi2, armStateMPC_rd2, armStateMPC_dzcc2);
armStateMPC_LA_DIAG_FORWARDBACKWARDSUB_12(armStateMPC_Phi3, armStateMPC_rd3, armStateMPC_dzcc3);
armStateMPC_LA_DIAG_FORWARDBACKWARDSUB_12(armStateMPC_Phi4, armStateMPC_rd4, armStateMPC_dzcc4);
armStateMPC_LA_DIAG_FORWARDBACKWARDSUB_12(armStateMPC_Phi5, armStateMPC_rd5, armStateMPC_dzcc5);
armStateMPC_LA_DIAG_FORWARDBACKWARDSUB_12(armStateMPC_Phi6, armStateMPC_rd6, armStateMPC_dzcc6);
armStateMPC_LA_DIAG_FORWARDBACKWARDSUB_12(armStateMPC_Phi7, armStateMPC_rd7, armStateMPC_dzcc7);
armStateMPC_LA_DIAG_FORWARDBACKWARDSUB_12(armStateMPC_Phi8, armStateMPC_rd8, armStateMPC_dzcc8);
armStateMPC_LA_DIAG_FORWARDBACKWARDSUB_6(armStateMPC_Phi9, armStateMPC_rd9, armStateMPC_dzcc9);
armStateMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_12(armStateMPC_ccrhsl0, armStateMPC_slb0, armStateMPC_llbbyslb0, armStateMPC_dzcc0, armStateMPC_lbIdx0, armStateMPC_dllbcc0);
armStateMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_12(armStateMPC_ccrhsub0, armStateMPC_sub0, armStateMPC_lubbysub0, armStateMPC_dzcc0, armStateMPC_ubIdx0, armStateMPC_dlubcc0);
armStateMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_12(armStateMPC_ccrhsl1, armStateMPC_slb1, armStateMPC_llbbyslb1, armStateMPC_dzcc1, armStateMPC_lbIdx1, armStateMPC_dllbcc1);
armStateMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_12(armStateMPC_ccrhsub1, armStateMPC_sub1, armStateMPC_lubbysub1, armStateMPC_dzcc1, armStateMPC_ubIdx1, armStateMPC_dlubcc1);
armStateMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_12(armStateMPC_ccrhsl2, armStateMPC_slb2, armStateMPC_llbbyslb2, armStateMPC_dzcc2, armStateMPC_lbIdx2, armStateMPC_dllbcc2);
armStateMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_12(armStateMPC_ccrhsub2, armStateMPC_sub2, armStateMPC_lubbysub2, armStateMPC_dzcc2, armStateMPC_ubIdx2, armStateMPC_dlubcc2);
armStateMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_12(armStateMPC_ccrhsl3, armStateMPC_slb3, armStateMPC_llbbyslb3, armStateMPC_dzcc3, armStateMPC_lbIdx3, armStateMPC_dllbcc3);
armStateMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_12(armStateMPC_ccrhsub3, armStateMPC_sub3, armStateMPC_lubbysub3, armStateMPC_dzcc3, armStateMPC_ubIdx3, armStateMPC_dlubcc3);
armStateMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_12(armStateMPC_ccrhsl4, armStateMPC_slb4, armStateMPC_llbbyslb4, armStateMPC_dzcc4, armStateMPC_lbIdx4, armStateMPC_dllbcc4);
armStateMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_12(armStateMPC_ccrhsub4, armStateMPC_sub4, armStateMPC_lubbysub4, armStateMPC_dzcc4, armStateMPC_ubIdx4, armStateMPC_dlubcc4);
armStateMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_12(armStateMPC_ccrhsl5, armStateMPC_slb5, armStateMPC_llbbyslb5, armStateMPC_dzcc5, armStateMPC_lbIdx5, armStateMPC_dllbcc5);
armStateMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_12(armStateMPC_ccrhsub5, armStateMPC_sub5, armStateMPC_lubbysub5, armStateMPC_dzcc5, armStateMPC_ubIdx5, armStateMPC_dlubcc5);
armStateMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_12(armStateMPC_ccrhsl6, armStateMPC_slb6, armStateMPC_llbbyslb6, armStateMPC_dzcc6, armStateMPC_lbIdx6, armStateMPC_dllbcc6);
armStateMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_12(armStateMPC_ccrhsub6, armStateMPC_sub6, armStateMPC_lubbysub6, armStateMPC_dzcc6, armStateMPC_ubIdx6, armStateMPC_dlubcc6);
armStateMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_12(armStateMPC_ccrhsl7, armStateMPC_slb7, armStateMPC_llbbyslb7, armStateMPC_dzcc7, armStateMPC_lbIdx7, armStateMPC_dllbcc7);
armStateMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_12(armStateMPC_ccrhsub7, armStateMPC_sub7, armStateMPC_lubbysub7, armStateMPC_dzcc7, armStateMPC_ubIdx7, armStateMPC_dlubcc7);
armStateMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_12(armStateMPC_ccrhsl8, armStateMPC_slb8, armStateMPC_llbbyslb8, armStateMPC_dzcc8, armStateMPC_lbIdx8, armStateMPC_dllbcc8);
armStateMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_12(armStateMPC_ccrhsub8, armStateMPC_sub8, armStateMPC_lubbysub8, armStateMPC_dzcc8, armStateMPC_ubIdx8, armStateMPC_dlubcc8);
armStateMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_6(armStateMPC_ccrhsl9, armStateMPC_slb9, armStateMPC_llbbyslb9, armStateMPC_dzcc9, armStateMPC_lbIdx9, armStateMPC_dllbcc9);
armStateMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_6(armStateMPC_ccrhsub9, armStateMPC_sub9, armStateMPC_lubbysub9, armStateMPC_dzcc9, armStateMPC_ubIdx9, armStateMPC_dlubcc9);
armStateMPC_LA_VSUB7_228(armStateMPC_l, armStateMPC_ccrhs, armStateMPC_s, armStateMPC_dl_cc, armStateMPC_ds_cc);
armStateMPC_LA_VADD_114(armStateMPC_dz_cc, armStateMPC_dz_aff);
armStateMPC_LA_VADD_60(armStateMPC_dv_cc, armStateMPC_dv_aff);
armStateMPC_LA_VADD_228(armStateMPC_dl_cc, armStateMPC_dl_aff);
armStateMPC_LA_VADD_228(armStateMPC_ds_cc, armStateMPC_ds_aff);
info->lsit_cc = armStateMPC_LINESEARCH_BACKTRACKING_COMBINED(armStateMPC_z, armStateMPC_v, armStateMPC_l, armStateMPC_s, armStateMPC_dz_cc, armStateMPC_dv_cc, armStateMPC_dl_cc, armStateMPC_ds_cc, &info->step_cc, &info->mu);
if( info->lsit_cc == armStateMPC_NOPROGRESS ){
exitcode = armStateMPC_NOPROGRESS; break;
}
info->it++;
}
output->z1[0] = armStateMPC_z0[0];
output->z1[1] = armStateMPC_z0[1];
output->z1[2] = armStateMPC_z0[2];
output->z1[3] = armStateMPC_z0[3];
output->z1[4] = armStateMPC_z0[4];
output->z1[5] = armStateMPC_z0[5];
output->z1[6] = armStateMPC_z0[6];
output->z1[7] = armStateMPC_z0[7];
output->z1[8] = armStateMPC_z0[8];
output->z1[9] = armStateMPC_z0[9];
output->z1[10] = armStateMPC_z0[10];
output->z1[11] = armStateMPC_z0[11];
output->z2[0] = armStateMPC_z1[0];
output->z2[1] = armStateMPC_z1[1];
output->z2[2] = armStateMPC_z1[2];
output->z2[3] = armStateMPC_z1[3];
output->z2[4] = armStateMPC_z1[4];
output->z2[5] = armStateMPC_z1[5];
output->z2[6] = armStateMPC_z1[6];
output->z2[7] = armStateMPC_z1[7];
output->z2[8] = armStateMPC_z1[8];
output->z2[9] = armStateMPC_z1[9];
output->z2[10] = armStateMPC_z1[10];
output->z2[11] = armStateMPC_z1[11];
output->z3[0] = armStateMPC_z2[0];
output->z3[1] = armStateMPC_z2[1];
output->z3[2] = armStateMPC_z2[2];
output->z3[3] = armStateMPC_z2[3];
output->z3[4] = armStateMPC_z2[4];
output->z3[5] = armStateMPC_z2[5];
output->z3[6] = armStateMPC_z2[6];
output->z3[7] = armStateMPC_z2[7];
output->z3[8] = armStateMPC_z2[8];
output->z3[9] = armStateMPC_z2[9];
output->z3[10] = armStateMPC_z2[10];
output->z3[11] = armStateMPC_z2[11];
output->z4[0] = armStateMPC_z3[0];
output->z4[1] = armStateMPC_z3[1];
output->z4[2] = armStateMPC_z3[2];
output->z4[3] = armStateMPC_z3[3];
output->z4[4] = armStateMPC_z3[4];
output->z4[5] = armStateMPC_z3[5];
output->z4[6] = armStateMPC_z3[6];
output->z4[7] = armStateMPC_z3[7];
output->z4[8] = armStateMPC_z3[8];
output->z4[9] = armStateMPC_z3[9];
output->z4[10] = armStateMPC_z3[10];
output->z4[11] = armStateMPC_z3[11];
output->z5[0] = armStateMPC_z4[0];
output->z5[1] = armStateMPC_z4[1];
output->z5[2] = armStateMPC_z4[2];
output->z5[3] = armStateMPC_z4[3];
output->z5[4] = armStateMPC_z4[4];
output->z5[5] = armStateMPC_z4[5];
output->z5[6] = armStateMPC_z4[6];
output->z5[7] = armStateMPC_z4[7];
output->z5[8] = armStateMPC_z4[8];
output->z5[9] = armStateMPC_z4[9];
output->z5[10] = armStateMPC_z4[10];
output->z5[11] = armStateMPC_z4[11];
output->z6[0] = armStateMPC_z5[0];
output->z6[1] = armStateMPC_z5[1];
output->z6[2] = armStateMPC_z5[2];
output->z6[3] = armStateMPC_z5[3];
output->z6[4] = armStateMPC_z5[4];
output->z6[5] = armStateMPC_z5[5];
output->z6[6] = armStateMPC_z5[6];
output->z6[7] = armStateMPC_z5[7];
output->z6[8] = armStateMPC_z5[8];
output->z6[9] = armStateMPC_z5[9];
output->z6[10] = armStateMPC_z5[10];
output->z6[11] = armStateMPC_z5[11];
output->z7[0] = armStateMPC_z6[0];
output->z7[1] = armStateMPC_z6[1];
output->z7[2] = armStateMPC_z6[2];
output->z7[3] = armStateMPC_z6[3];
output->z7[4] = armStateMPC_z6[4];
output->z7[5] = armStateMPC_z6[5];
output->z7[6] = armStateMPC_z6[6];
output->z7[7] = armStateMPC_z6[7];
output->z7[8] = armStateMPC_z6[8];
output->z7[9] = armStateMPC_z6[9];
output->z7[10] = armStateMPC_z6[10];
output->z7[11] = armStateMPC_z6[11];
output->z8[0] = armStateMPC_z7[0];
output->z8[1] = armStateMPC_z7[1];
output->z8[2] = armStateMPC_z7[2];
output->z8[3] = armStateMPC_z7[3];
output->z8[4] = armStateMPC_z7[4];
output->z8[5] = armStateMPC_z7[5];
output->z8[6] = armStateMPC_z7[6];
output->z8[7] = armStateMPC_z7[7];
output->z8[8] = armStateMPC_z7[8];
output->z8[9] = armStateMPC_z7[9];
output->z8[10] = armStateMPC_z7[10];
output->z8[11] = armStateMPC_z7[11];
output->z9[0] = armStateMPC_z8[0];
output->z9[1] = armStateMPC_z8[1];
output->z9[2] = armStateMPC_z8[2];
output->z9[3] = armStateMPC_z8[3];
output->z9[4] = armStateMPC_z8[4];
output->z9[5] = armStateMPC_z8[5];
output->z9[6] = armStateMPC_z8[6];
output->z9[7] = armStateMPC_z8[7];
output->z9[8] = armStateMPC_z8[8];
output->z9[9] = armStateMPC_z8[9];
output->z9[10] = armStateMPC_z8[10];
output->z9[11] = armStateMPC_z8[11];
output->z10[0] = armStateMPC_z9[0];
output->z10[1] = armStateMPC_z9[1];
output->z10[2] = armStateMPC_z9[2];
output->z10[3] = armStateMPC_z9[3];
output->z10[4] = armStateMPC_z9[4];
output->z10[5] = armStateMPC_z9[5];

#if armStateMPC_SET_TIMING == 1
info->solvetime = armStateMPC_toc(&solvertimer);
#if armStateMPC_SET_PRINTLEVEL > 0 && armStateMPC_SET_TIMING == 1
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
