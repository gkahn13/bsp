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
 * Initializes a vector of length 96 with a value.
 */
void planarMPC_LA_INITIALIZEVECTOR_96(planarMPC_FLOAT* vec, planarMPC_FLOAT value)
{
	int i;
	for( i=0; i<96; i++ )
	{
		vec[i] = value;
	}
}


/*
 * Initializes a vector of length 60 with a value.
 */
void planarMPC_LA_INITIALIZEVECTOR_60(planarMPC_FLOAT* vec, planarMPC_FLOAT value)
{
	int i;
	for( i=0; i<60; i++ )
	{
		vec[i] = value;
	}
}


/*
 * Initializes a vector of length 192 with a value.
 */
void planarMPC_LA_INITIALIZEVECTOR_192(planarMPC_FLOAT* vec, planarMPC_FLOAT value)
{
	int i;
	for( i=0; i<192; i++ )
	{
		vec[i] = value;
	}
}


/* 
 * Calculates a dot product and adds it to a variable: z += x'*y; 
 * This function is for vectors of length 192.
 */
void planarMPC_LA_DOTACC_192(planarMPC_FLOAT *x, planarMPC_FLOAT *y, planarMPC_FLOAT *z)
{
	int i;
	for( i=0; i<192; i++ ){
		*z += x[i]*y[i];
	}
}


/*
 * Calculates the gradient and the value for a quadratic function 0.5*z'*H*z + f'*z
 *
 * INPUTS:     H  - Symmetric Hessian, diag matrix of size [10 x 10]
 *             f  - column vector of size 10
 *             z  - column vector of size 10
 *
 * OUTPUTS: grad  - gradient at z (= H*z + f), column vector of size 10
 *          value <-- value + 0.5*z'*H*z + f'*z (value will be modified)
 */
void planarMPC_LA_DIAG_QUADFCN_10(planarMPC_FLOAT* H, planarMPC_FLOAT* f, planarMPC_FLOAT* z, planarMPC_FLOAT* grad, planarMPC_FLOAT* value)
{
	int i;
	planarMPC_FLOAT hz;	
	for( i=0; i<10; i++){
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
void planarMPC_LA_DIAG_QUADFCN_6(planarMPC_FLOAT* H, planarMPC_FLOAT* f, planarMPC_FLOAT* z, planarMPC_FLOAT* grad, planarMPC_FLOAT* value)
{
	int i;
	planarMPC_FLOAT hz;	
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
void planarMPC_LA_DIAGZERO_MVMSUB6_6(planarMPC_FLOAT *B, planarMPC_FLOAT *u, planarMPC_FLOAT *b, planarMPC_FLOAT *l, planarMPC_FLOAT *r, planarMPC_FLOAT *z, planarMPC_FLOAT *y)
{
	int i;
	planarMPC_FLOAT Bu[6];
	planarMPC_FLOAT norm = *y;
	planarMPC_FLOAT lr = 0;

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
void planarMPC_LA_DENSE_DIAGZERO_MVMSUB3_6_10_10(planarMPC_FLOAT *A, planarMPC_FLOAT *x, planarMPC_FLOAT *B, planarMPC_FLOAT *u, planarMPC_FLOAT *b, planarMPC_FLOAT *l, planarMPC_FLOAT *r, planarMPC_FLOAT *z, planarMPC_FLOAT *y)
{
	int i;
	int j;
	int k = 0;
	planarMPC_FLOAT AxBu[6];
	planarMPC_FLOAT norm = *y;
	planarMPC_FLOAT lr = 0;

	/* do A*x + B*u first */
	for( i=0; i<6; i++ ){
		AxBu[i] = A[k++]*x[0] + B[i]*u[i];
	}	

	for( j=1; j<10; j++ ){		
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
void planarMPC_LA_DENSE_DIAGZERO_MVMSUB3_6_10_6(planarMPC_FLOAT *A, planarMPC_FLOAT *x, planarMPC_FLOAT *B, planarMPC_FLOAT *u, planarMPC_FLOAT *b, planarMPC_FLOAT *l, planarMPC_FLOAT *r, planarMPC_FLOAT *z, planarMPC_FLOAT *y)
{
	int i;
	int j;
	int k = 0;
	planarMPC_FLOAT AxBu[6];
	planarMPC_FLOAT norm = *y;
	planarMPC_FLOAT lr = 0;

	/* do A*x + B*u first */
	for( i=0; i<6; i++ ){
		AxBu[i] = A[k++]*x[0] + B[i]*u[i];
	}	

	for( j=1; j<10; j++ ){		
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
 * where A is of size [6 x 10] and stored in column major format.
 * and B is of size [6 x 10] and stored in diagzero format
 * Note the transposes of A and B!
 */
void planarMPC_LA_DENSE_DIAGZERO_MTVM2_6_10_6(planarMPC_FLOAT *A, planarMPC_FLOAT *x, planarMPC_FLOAT *B, planarMPC_FLOAT *y, planarMPC_FLOAT *z)
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
	for( i=6 ;i<10; i++ ){
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
void planarMPC_LA_DIAGZERO_MTVM_6_6(planarMPC_FLOAT *M, planarMPC_FLOAT *x, planarMPC_FLOAT *y)
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
 * for vectors of length 10. Output z is of course scalar.
 */
void planarMPC_LA_VSUBADD3_10(planarMPC_FLOAT* t, planarMPC_FLOAT* u, int* uidx, planarMPC_FLOAT* v, planarMPC_FLOAT* w, planarMPC_FLOAT* y, planarMPC_FLOAT* z, planarMPC_FLOAT* r)
{
	int i;
	planarMPC_FLOAT norm = *r;
	planarMPC_FLOAT vx = 0;
	planarMPC_FLOAT x;
	for( i=0; i<10; i++){
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
 * for vectors of length 10. Output z is of course scalar.
 */
void planarMPC_LA_VSUBADD2_10(planarMPC_FLOAT* t, int* tidx, planarMPC_FLOAT* u, planarMPC_FLOAT* v, planarMPC_FLOAT* w, planarMPC_FLOAT* y, planarMPC_FLOAT* z, planarMPC_FLOAT* r)
{
	int i;
	planarMPC_FLOAT norm = *r;
	planarMPC_FLOAT vx = 0;
	planarMPC_FLOAT x;
	for( i=0; i<10; i++){
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
void planarMPC_LA_VSUBADD3_6(planarMPC_FLOAT* t, planarMPC_FLOAT* u, int* uidx, planarMPC_FLOAT* v, planarMPC_FLOAT* w, planarMPC_FLOAT* y, planarMPC_FLOAT* z, planarMPC_FLOAT* r)
{
	int i;
	planarMPC_FLOAT norm = *r;
	planarMPC_FLOAT vx = 0;
	planarMPC_FLOAT x;
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
void planarMPC_LA_VSUBADD2_6(planarMPC_FLOAT* t, int* tidx, planarMPC_FLOAT* u, planarMPC_FLOAT* v, planarMPC_FLOAT* w, planarMPC_FLOAT* y, planarMPC_FLOAT* z, planarMPC_FLOAT* r)
{
	int i;
	planarMPC_FLOAT norm = *r;
	planarMPC_FLOAT vx = 0;
	planarMPC_FLOAT x;
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
 * Special function for box constraints of length 10
 * Returns also L/S, a value that is often used elsewhere.
 */
void planarMPC_LA_INEQ_B_GRAD_10_10_10(planarMPC_FLOAT *lu, planarMPC_FLOAT *su, planarMPC_FLOAT *ru, planarMPC_FLOAT *ll, planarMPC_FLOAT *sl, planarMPC_FLOAT *rl, int* lbIdx, int* ubIdx, planarMPC_FLOAT *grad, planarMPC_FLOAT *lubysu, planarMPC_FLOAT *llbysl)
{
	int i;
	for( i=0; i<10; i++ ){
		grad[i] = 0;
	}
	for( i=0; i<10; i++ ){		
		llbysl[i] = ll[i] / sl[i];
		grad[lbIdx[i]] -= llbysl[i]*rl[i];
	}
	for( i=0; i<10; i++ ){
		lubysu[i] = lu[i] / su[i];
		grad[ubIdx[i]] += lubysu[i]*ru[i];
	}
}


/*
 * Computes inequality constraints gradient-
 * Special function for box constraints of length 6
 * Returns also L/S, a value that is often used elsewhere.
 */
void planarMPC_LA_INEQ_B_GRAD_6_6_6(planarMPC_FLOAT *lu, planarMPC_FLOAT *su, planarMPC_FLOAT *ru, planarMPC_FLOAT *ll, planarMPC_FLOAT *sl, planarMPC_FLOAT *rl, int* lbIdx, int* ubIdx, planarMPC_FLOAT *grad, planarMPC_FLOAT *lubysu, planarMPC_FLOAT *llbysl)
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
 * of length 96.
 */
void planarMPC_LA_VVADD3_96(planarMPC_FLOAT *u, planarMPC_FLOAT *v, planarMPC_FLOAT *w, planarMPC_FLOAT *z)
{
	int i;
	for( i=0; i<96; i++ ){
		z[i] = u[i] + v[i] + w[i];
	}
}


/*
 * Special function to compute the diagonal cholesky factorization of the 
 * positive definite augmented Hessian for block size 10.
 *
 * Inputs: - H = diagonal cost Hessian in diagonal storage format
 *         - llbysl = L / S of lower bounds
 *         - lubysu = L / S of upper bounds
 *
 * Output: Phi = sqrt(H + diag(llbysl) + diag(lubysu))
 * where Phi is stored in diagonal storage format
 */
void planarMPC_LA_DIAG_CHOL_ONELOOP_LBUB_10_10_10(planarMPC_FLOAT *H, planarMPC_FLOAT *llbysl, int* lbIdx, planarMPC_FLOAT *lubysu, int* ubIdx, planarMPC_FLOAT *Phi)


{
	int i;
	
	/* compute cholesky */
	for( i=0; i<10; i++ ){
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
 * where A is to be computed and is of size [6 x 10],
 * B is given and of size [6 x 10], L is a diagonal
 * matrix of size 6 stored in diagonal matrix 
 * storage format. Note the transpose of L has no impact!
 *
 * Result: A in column major storage format.
 *
 */
void planarMPC_LA_DIAG_MATRIXFORWARDSUB_6_10(planarMPC_FLOAT *L, planarMPC_FLOAT *B, planarMPC_FLOAT *A)
{
    int i,j;
	 int k = 0;

	for( j=0; j<10; j++){
		for( i=0; i<6; i++){
			A[k] = B[k]/L[j];
			k++;
		}
	}

}


/**
 * Forward substitution for the matrix equation A*L' = B
 * where A is to be computed and is of size [6 x 10],
 * B is given and of size [6 x 10], L is a diagonal
 *  matrix of size 10 stored in diagonal 
 * storage format. Note the transpose of L!
 *
 * Result: A in diagonalzero storage format.
 *
 */
void planarMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_6_10(planarMPC_FLOAT *L, planarMPC_FLOAT *B, planarMPC_FLOAT *A)
{
	int j;
    for( j=0; j<10; j++ ){   
		A[j] = B[j]/L[j];
     }
}


/**
 * Compute C = A*B' where 
 *
 *	size(A) = [6 x 10]
 *  size(B) = [6 x 10] in diagzero format
 * 
 * A and C matrices are stored in column major format.
 * 
 * 
 */
void planarMPC_LA_DENSE_DIAGZERO_MMTM_6_10_6(planarMPC_FLOAT *A, planarMPC_FLOAT *B, planarMPC_FLOAT *C)
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
 * The dimensions involved are 10.
 */
void planarMPC_LA_DIAG_FORWARDSUB_10(planarMPC_FLOAT *L, planarMPC_FLOAT *b, planarMPC_FLOAT *y)
{
    int i;

    for( i=0; i<10; i++ ){
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
void planarMPC_LA_DIAG_CHOL_ONELOOP_LBUB_6_6_6(planarMPC_FLOAT *H, planarMPC_FLOAT *llbysl, int* lbIdx, planarMPC_FLOAT *lubysu, int* ubIdx, planarMPC_FLOAT *Phi)


{
	int i;
	
	/* compute cholesky */
	for( i=0; i<6; i++ ){
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
 * where A is to be computed and is of size [6 x 6],
 * B is given and of size [6 x 6], L is a diagonal
 *  matrix of size 6 stored in diagonal 
 * storage format. Note the transpose of L!
 *
 * Result: A in diagonalzero storage format.
 *
 */
void planarMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_6_6(planarMPC_FLOAT *L, planarMPC_FLOAT *B, planarMPC_FLOAT *A)
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
void planarMPC_LA_DIAG_FORWARDSUB_6(planarMPC_FLOAT *L, planarMPC_FLOAT *b, planarMPC_FLOAT *y)
{
    int i;

    for( i=0; i<6; i++ ){
		y[i] = b[i]/L[i];
    }
}


/**
 * Compute L = B*B', where L is lower triangular of size NXp1
 * and B is a diagzero matrix of size [6 x 10] in column
 * storage format.
 * 
 */
void planarMPC_LA_DIAGZERO_MMT_6(planarMPC_FLOAT *B, planarMPC_FLOAT *L)
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
void planarMPC_LA_DIAGZERO_MVMSUB7_6(planarMPC_FLOAT *B, planarMPC_FLOAT *u, planarMPC_FLOAT *b, planarMPC_FLOAT *r)
{
	int i;

	for( i=0; i<6; i++ ){
		r[i] = b[i] - B[i]*u[i];
	}	
	
}


/**
 * Compute L = A*A' + B*B', where L is lower triangular of size NXp1
 * and A is a dense matrix of size [6 x 10] in column
 * storage format, and B is of size [6 x 10] diagonalzero
 * storage format.
 * 
 * THIS ONE HAS THE WORST ACCES PATTERN POSSIBLE. 
 * POSSIBKE FIX: PUT A AND B INTO ROW MAJOR FORMAT FIRST.
 * 
 */
void planarMPC_LA_DENSE_DIAGZERO_MMT2_6_10_10(planarMPC_FLOAT *A, planarMPC_FLOAT *B, planarMPC_FLOAT *L)
{
    int i, j, k, ii, di;
    planarMPC_FLOAT ltemp;
    
    ii = 0; di = 0;
    for( i=0; i<6; i++ ){        
        for( j=0; j<=i; j++ ){
            ltemp = 0; 
            for( k=0; k<10; k++ ){
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
void planarMPC_LA_DENSE_DIAGZERO_2MVMSUB2_6_10_10(planarMPC_FLOAT *A, planarMPC_FLOAT *x, planarMPC_FLOAT *B, planarMPC_FLOAT *u, planarMPC_FLOAT *b, planarMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;

	for( i=0; i<6; i++ ){
		r[i] = b[i] - A[k++]*x[0] - B[i]*u[i];
	}	

	for( j=1; j<10; j++ ){		
		for( i=0; i<6; i++ ){
			r[i] -= A[k++]*x[j];
		}
	}
	
}


/**
 * Compute L = A*A' + B*B', where L is lower triangular of size NXp1
 * and A is a dense matrix of size [6 x 10] in column
 * storage format, and B is of size [6 x 6] diagonalzero
 * storage format.
 * 
 * THIS ONE HAS THE WORST ACCES PATTERN POSSIBLE. 
 * POSSIBKE FIX: PUT A AND B INTO ROW MAJOR FORMAT FIRST.
 * 
 */
void planarMPC_LA_DENSE_DIAGZERO_MMT2_6_10_6(planarMPC_FLOAT *A, planarMPC_FLOAT *B, planarMPC_FLOAT *L)
{
    int i, j, k, ii, di;
    planarMPC_FLOAT ltemp;
    
    ii = 0; di = 0;
    for( i=0; i<6; i++ ){        
        for( j=0; j<=i; j++ ){
            ltemp = 0; 
            for( k=0; k<10; k++ ){
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
void planarMPC_LA_DENSE_DIAGZERO_2MVMSUB2_6_10_6(planarMPC_FLOAT *A, planarMPC_FLOAT *x, planarMPC_FLOAT *B, planarMPC_FLOAT *u, planarMPC_FLOAT *b, planarMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;

	for( i=0; i<6; i++ ){
		r[i] = b[i] - A[k++]*x[0] - B[i]*u[i];
	}	

	for( j=1; j<10; j++ ){		
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
void planarMPC_LA_DENSE_CHOL_6(planarMPC_FLOAT *A, planarMPC_FLOAT *L)
{
    int i, j, k, di, dj;
	 int ii, jj;

    planarMPC_FLOAT l;
    planarMPC_FLOAT Mii;

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
void planarMPC_LA_DENSE_FORWARDSUB_6(planarMPC_FLOAT *L, planarMPC_FLOAT *b, planarMPC_FLOAT *y)
{
    int i,j,ii,di;
    planarMPC_FLOAT yel;
            
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
void planarMPC_LA_DENSE_MATRIXTFORWARDSUB_6_6(planarMPC_FLOAT *L, planarMPC_FLOAT *B, planarMPC_FLOAT *A)
{
    int i,j,k,ii,di;
    planarMPC_FLOAT a;
    
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
void planarMPC_LA_DENSE_MMTSUB_6_6(planarMPC_FLOAT *A, planarMPC_FLOAT *L)
{
    int i, j, k, ii, di;
    planarMPC_FLOAT ltemp;
    
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
void planarMPC_LA_DENSE_MVMSUB1_6_6(planarMPC_FLOAT *A, planarMPC_FLOAT *x, planarMPC_FLOAT *b, planarMPC_FLOAT *r)
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
void planarMPC_LA_DENSE_BACKWARDSUB_6(planarMPC_FLOAT *L, planarMPC_FLOAT *y, planarMPC_FLOAT *x)
{
    int i, ii, di, j, jj, dj;
    planarMPC_FLOAT xel;    
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
void planarMPC_LA_DENSE_MTVMSUB_6_6(planarMPC_FLOAT *A, planarMPC_FLOAT *x, planarMPC_FLOAT *b, planarMPC_FLOAT *r)
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
 * Vector subtraction z = -x - y for vectors of length 96.
 */
void planarMPC_LA_VSUB2_96(planarMPC_FLOAT *x, planarMPC_FLOAT *y, planarMPC_FLOAT *z)
{
	int i;
	for( i=0; i<96; i++){
		z[i] = -x[i] - y[i];
	}
}


/**
 * Forward-Backward-Substitution to solve L*L^T*x = b where L is a
 * diagonal matrix of size 10 in vector
 * storage format.
 */
void planarMPC_LA_DIAG_FORWARDBACKWARDSUB_10(planarMPC_FLOAT *L, planarMPC_FLOAT *b, planarMPC_FLOAT *x)
{
    int i;
            
    /* solve Ly = b by forward and backward substitution */
    for( i=0; i<10; i++ ){
		x[i] = b[i]/(L[i]*L[i]);
    }
    
}


/**
 * Forward-Backward-Substitution to solve L*L^T*x = b where L is a
 * diagonal matrix of size 6 in vector
 * storage format.
 */
void planarMPC_LA_DIAG_FORWARDBACKWARDSUB_6(planarMPC_FLOAT *L, planarMPC_FLOAT *b, planarMPC_FLOAT *x)
{
    int i;
            
    /* solve Ly = b by forward and backward substitution */
    for( i=0; i<6; i++ ){
		x[i] = b[i]/(L[i]*L[i]);
    }
    
}


/*
 * Vector subtraction z = x(xidx) - y where y, z and xidx are of length 10,
 * and x has length 10 and is indexed through yidx.
 */
void planarMPC_LA_VSUB_INDEXED_10(planarMPC_FLOAT *x, int* xidx, planarMPC_FLOAT *y, planarMPC_FLOAT *z)
{
	int i;
	for( i=0; i<10; i++){
		z[i] = x[xidx[i]] - y[i];
	}
}


/*
 * Vector subtraction x = -u.*v - w for vectors of length 10.
 */
void planarMPC_LA_VSUB3_10(planarMPC_FLOAT *u, planarMPC_FLOAT *v, planarMPC_FLOAT *w, planarMPC_FLOAT *x)
{
	int i;
	for( i=0; i<10; i++){
		x[i] = -u[i]*v[i] - w[i];
	}
}


/*
 * Vector subtraction z = -x - y(yidx) where y is of length 10
 * and z, x and yidx are of length 10.
 */
void planarMPC_LA_VSUB2_INDEXED_10(planarMPC_FLOAT *x, planarMPC_FLOAT *y, int* yidx, planarMPC_FLOAT *z)
{
	int i;
	for( i=0; i<10; i++){
		z[i] = -x[i] - y[yidx[i]];
	}
}


/*
 * Vector subtraction z = x(xidx) - y where y, z and xidx are of length 6,
 * and x has length 6 and is indexed through yidx.
 */
void planarMPC_LA_VSUB_INDEXED_6(planarMPC_FLOAT *x, int* xidx, planarMPC_FLOAT *y, planarMPC_FLOAT *z)
{
	int i;
	for( i=0; i<6; i++){
		z[i] = x[xidx[i]] - y[i];
	}
}


/*
 * Vector subtraction x = -u.*v - w for vectors of length 6.
 */
void planarMPC_LA_VSUB3_6(planarMPC_FLOAT *u, planarMPC_FLOAT *v, planarMPC_FLOAT *w, planarMPC_FLOAT *x)
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
void planarMPC_LA_VSUB2_INDEXED_6(planarMPC_FLOAT *x, planarMPC_FLOAT *y, int* yidx, planarMPC_FLOAT *z)
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
        for( i=0; i<192; i++ ){
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
        if( i == 192 ){
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
    *mu_aff = mymu / (planarMPC_FLOAT)192;
    return lsIt;
}


/*
 * Vector subtraction x = u.*v - a where a is a scalar
*  and x,u,v are vectors of length 192.
 */
void planarMPC_LA_VSUB5_192(planarMPC_FLOAT *u, planarMPC_FLOAT *v, planarMPC_FLOAT a, planarMPC_FLOAT *x)
{
	int i;
	for( i=0; i<192; i++){
		x[i] = u[i]*v[i] - a;
	}
}


/*
 * Computes x=0; x(uidx) += u/su; x(vidx) -= v/sv where x is of length 10,
 * u, su, uidx are of length 10 and v, sv, vidx are of length 10.
 */
void planarMPC_LA_VSUB6_INDEXED_10_10_10(planarMPC_FLOAT *u, planarMPC_FLOAT *su, int* uidx, planarMPC_FLOAT *v, planarMPC_FLOAT *sv, int* vidx, planarMPC_FLOAT *x)
{
	int i;
	for( i=0; i<10; i++ ){
		x[i] = 0;
	}
	for( i=0; i<10; i++){
		x[uidx[i]] += u[i]/su[i];
	}
	for( i=0; i<10; i++){
		x[vidx[i]] -= v[i]/sv[i];
	}
}


/* 
 * Computes r =  B*u
 * where B is stored in diagzero format
 */
void planarMPC_LA_DIAGZERO_MVM_6(planarMPC_FLOAT *B, planarMPC_FLOAT *u, planarMPC_FLOAT *r)
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
void planarMPC_LA_DENSE_DIAGZERO_2MVMADD_6_10_10(planarMPC_FLOAT *A, planarMPC_FLOAT *x, planarMPC_FLOAT *B, planarMPC_FLOAT *u, planarMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;

	for( i=0; i<6; i++ ){
		r[i] = A[k++]*x[0] + B[i]*u[i];
	}	

	for( j=1; j<10; j++ ){		
		for( i=0; i<6; i++ ){
			r[i] += A[k++]*x[j];
		}
	}
	
}


/*
 * Computes x=0; x(uidx) += u/su; x(vidx) -= v/sv where x is of length 6,
 * u, su, uidx are of length 6 and v, sv, vidx are of length 6.
 */
void planarMPC_LA_VSUB6_INDEXED_6_6_6(planarMPC_FLOAT *u, planarMPC_FLOAT *su, int* uidx, planarMPC_FLOAT *v, planarMPC_FLOAT *sv, int* vidx, planarMPC_FLOAT *x)
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
void planarMPC_LA_DENSE_DIAGZERO_2MVMADD_6_10_6(planarMPC_FLOAT *A, planarMPC_FLOAT *x, planarMPC_FLOAT *B, planarMPC_FLOAT *u, planarMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;

	for( i=0; i<6; i++ ){
		r[i] = A[k++]*x[0] + B[i]*u[i];
	}	

	for( j=1; j<10; j++ ){		
		for( i=0; i<6; i++ ){
			r[i] += A[k++]*x[j];
		}
	}
	
}


/*
 * Vector subtraction z = x - y for vectors of length 96.
 */
void planarMPC_LA_VSUB_96(planarMPC_FLOAT *x, planarMPC_FLOAT *y, planarMPC_FLOAT *z)
{
	int i;
	for( i=0; i<96; i++){
		z[i] = x[i] - y[i];
	}
}


/** 
 * Computes z = -r./s - u.*y(y)
 * where all vectors except of y are of length 10 (length of y >= 10).
 */
void planarMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_10(planarMPC_FLOAT *r, planarMPC_FLOAT *s, planarMPC_FLOAT *u, planarMPC_FLOAT *y, int* yidx, planarMPC_FLOAT *z)
{
	int i;
	for( i=0; i<10; i++ ){
		z[i] = -r[i]/s[i] - u[i]*y[yidx[i]];
	}
}


/** 
 * Computes z = -r./s + u.*y(y)
 * where all vectors except of y are of length 10 (length of y >= 10).
 */
void planarMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_10(planarMPC_FLOAT *r, planarMPC_FLOAT *s, planarMPC_FLOAT *u, planarMPC_FLOAT *y, int* yidx, planarMPC_FLOAT *z)
{
	int i;
	for( i=0; i<10; i++ ){
		z[i] = -r[i]/s[i] + u[i]*y[yidx[i]];
	}
}


/** 
 * Computes z = -r./s - u.*y(y)
 * where all vectors except of y are of length 6 (length of y >= 6).
 */
void planarMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_6(planarMPC_FLOAT *r, planarMPC_FLOAT *s, planarMPC_FLOAT *u, planarMPC_FLOAT *y, int* yidx, planarMPC_FLOAT *z)
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
void planarMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_6(planarMPC_FLOAT *r, planarMPC_FLOAT *s, planarMPC_FLOAT *u, planarMPC_FLOAT *y, int* yidx, planarMPC_FLOAT *z)
{
	int i;
	for( i=0; i<6; i++ ){
		z[i] = -r[i]/s[i] + u[i]*y[yidx[i]];
	}
}


/*
 * Computes ds = -l.\(r + s.*dl) for vectors of length 192.
 */
void planarMPC_LA_VSUB7_192(planarMPC_FLOAT *l, planarMPC_FLOAT *r, planarMPC_FLOAT *s, planarMPC_FLOAT *dl, planarMPC_FLOAT *ds)
{
	int i;
	for( i=0; i<192; i++){
		ds[i] = -(r[i] + s[i]*dl[i])/l[i];
	}
}


/*
 * Vector addition x = x + y for vectors of length 96.
 */
void planarMPC_LA_VADD_96(planarMPC_FLOAT *x, planarMPC_FLOAT *y)
{
	int i;
	for( i=0; i<96; i++){
		x[i] += y[i];
	}
}


/*
 * Vector addition x = x + y for vectors of length 60.
 */
void planarMPC_LA_VADD_60(planarMPC_FLOAT *x, planarMPC_FLOAT *y)
{
	int i;
	for( i=0; i<60; i++){
		x[i] += y[i];
	}
}


/*
 * Vector addition x = x + y for vectors of length 192.
 */
void planarMPC_LA_VADD_192(planarMPC_FLOAT *x, planarMPC_FLOAT *y)
{
	int i;
	for( i=0; i<192; i++){
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
        for( i=0; i<192; i++ ){
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
        if( i == 192 ){
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
    for( i=0; i<96; i++ ){
        z[i] += a_gamma*dz[i];
    }
    
    /* equality constraint multipliers */
    for( i=0; i<60; i++ ){
        v[i] += a_gamma*dv[i];
    }
    
    /* inequality constraint multipliers & slacks, also update mu */
    *mu = 0;
    for( i=0; i<192; i++ ){
        dltemp = l[i] + a_gamma*dl[i]; l[i] = dltemp;
        dstemp = s[i] + a_gamma*ds[i]; s[i] = dstemp;
        *mu += dltemp*dstemp;
    }
    
    *a = a_gamma;
    *mu /= (planarMPC_FLOAT)192;
    return lsIt;
}




/* VARIABLE DEFINITIONS ------------------------------------------------ */
planarMPC_FLOAT planarMPC_z[96];
planarMPC_FLOAT planarMPC_v[60];
planarMPC_FLOAT planarMPC_dz_aff[96];
planarMPC_FLOAT planarMPC_dv_aff[60];
planarMPC_FLOAT planarMPC_grad_cost[96];
planarMPC_FLOAT planarMPC_grad_eq[96];
planarMPC_FLOAT planarMPC_rd[96];
planarMPC_FLOAT planarMPC_l[192];
planarMPC_FLOAT planarMPC_s[192];
planarMPC_FLOAT planarMPC_lbys[192];
planarMPC_FLOAT planarMPC_dl_aff[192];
planarMPC_FLOAT planarMPC_ds_aff[192];
planarMPC_FLOAT planarMPC_dz_cc[96];
planarMPC_FLOAT planarMPC_dv_cc[60];
planarMPC_FLOAT planarMPC_dl_cc[192];
planarMPC_FLOAT planarMPC_ds_cc[192];
planarMPC_FLOAT planarMPC_ccrhs[192];
planarMPC_FLOAT planarMPC_grad_ineq[96];
planarMPC_FLOAT* planarMPC_z0 = planarMPC_z + 0;
planarMPC_FLOAT* planarMPC_dzaff0 = planarMPC_dz_aff + 0;
planarMPC_FLOAT* planarMPC_dzcc0 = planarMPC_dz_cc + 0;
planarMPC_FLOAT* planarMPC_rd0 = planarMPC_rd + 0;
planarMPC_FLOAT planarMPC_Lbyrd0[10];
planarMPC_FLOAT* planarMPC_grad_cost0 = planarMPC_grad_cost + 0;
planarMPC_FLOAT* planarMPC_grad_eq0 = planarMPC_grad_eq + 0;
planarMPC_FLOAT* planarMPC_grad_ineq0 = planarMPC_grad_ineq + 0;
planarMPC_FLOAT planarMPC_ctv0[10];
planarMPC_FLOAT planarMPC_C0[60] = {1.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 
0.0000000000000000E+000, 1.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 
0.0000000000000000E+000, 0.0000000000000000E+000, 1.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 
0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 1.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 
0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 1.0000000000000000E+000, 0.0000000000000000E+000, 
0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 1.0000000000000000E+000, 
1.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 
0.0000000000000000E+000, 1.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 
0.0000000000000000E+000, 0.0000000000000000E+000, 1.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 
0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 1.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
planarMPC_FLOAT* planarMPC_v0 = planarMPC_v + 0;
planarMPC_FLOAT planarMPC_re0[6];
planarMPC_FLOAT planarMPC_beta0[6];
planarMPC_FLOAT planarMPC_betacc0[6];
planarMPC_FLOAT* planarMPC_dvaff0 = planarMPC_dv_aff + 0;
planarMPC_FLOAT* planarMPC_dvcc0 = planarMPC_dv_cc + 0;
planarMPC_FLOAT planarMPC_V0[60];
planarMPC_FLOAT planarMPC_Yd0[21];
planarMPC_FLOAT planarMPC_Ld0[21];
planarMPC_FLOAT planarMPC_yy0[6];
planarMPC_FLOAT planarMPC_bmy0[6];
int planarMPC_lbIdx0[10] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
planarMPC_FLOAT* planarMPC_llb0 = planarMPC_l + 0;
planarMPC_FLOAT* planarMPC_slb0 = planarMPC_s + 0;
planarMPC_FLOAT* planarMPC_llbbyslb0 = planarMPC_lbys + 0;
planarMPC_FLOAT planarMPC_rilb0[10];
planarMPC_FLOAT* planarMPC_dllbaff0 = planarMPC_dl_aff + 0;
planarMPC_FLOAT* planarMPC_dslbaff0 = planarMPC_ds_aff + 0;
planarMPC_FLOAT* planarMPC_dllbcc0 = planarMPC_dl_cc + 0;
planarMPC_FLOAT* planarMPC_dslbcc0 = planarMPC_ds_cc + 0;
planarMPC_FLOAT* planarMPC_ccrhsl0 = planarMPC_ccrhs + 0;
int planarMPC_ubIdx0[10] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
planarMPC_FLOAT* planarMPC_lub0 = planarMPC_l + 10;
planarMPC_FLOAT* planarMPC_sub0 = planarMPC_s + 10;
planarMPC_FLOAT* planarMPC_lubbysub0 = planarMPC_lbys + 10;
planarMPC_FLOAT planarMPC_riub0[10];
planarMPC_FLOAT* planarMPC_dlubaff0 = planarMPC_dl_aff + 10;
planarMPC_FLOAT* planarMPC_dsubaff0 = planarMPC_ds_aff + 10;
planarMPC_FLOAT* planarMPC_dlubcc0 = planarMPC_dl_cc + 10;
planarMPC_FLOAT* planarMPC_dsubcc0 = planarMPC_ds_cc + 10;
planarMPC_FLOAT* planarMPC_ccrhsub0 = planarMPC_ccrhs + 10;
planarMPC_FLOAT planarMPC_Phi0[10];
planarMPC_FLOAT planarMPC_D0[10] = {1.0000000000000000E+000, 
1.0000000000000000E+000, 
1.0000000000000000E+000, 
1.0000000000000000E+000, 
1.0000000000000000E+000, 
1.0000000000000000E+000};
planarMPC_FLOAT planarMPC_W0[10];
planarMPC_FLOAT* planarMPC_z1 = planarMPC_z + 10;
planarMPC_FLOAT* planarMPC_dzaff1 = planarMPC_dz_aff + 10;
planarMPC_FLOAT* planarMPC_dzcc1 = planarMPC_dz_cc + 10;
planarMPC_FLOAT* planarMPC_rd1 = planarMPC_rd + 10;
planarMPC_FLOAT planarMPC_Lbyrd1[10];
planarMPC_FLOAT* planarMPC_grad_cost1 = planarMPC_grad_cost + 10;
planarMPC_FLOAT* planarMPC_grad_eq1 = planarMPC_grad_eq + 10;
planarMPC_FLOAT* planarMPC_grad_ineq1 = planarMPC_grad_ineq + 10;
planarMPC_FLOAT planarMPC_ctv1[10];
planarMPC_FLOAT* planarMPC_v1 = planarMPC_v + 6;
planarMPC_FLOAT planarMPC_re1[6];
planarMPC_FLOAT planarMPC_beta1[6];
planarMPC_FLOAT planarMPC_betacc1[6];
planarMPC_FLOAT* planarMPC_dvaff1 = planarMPC_dv_aff + 6;
planarMPC_FLOAT* planarMPC_dvcc1 = planarMPC_dv_cc + 6;
planarMPC_FLOAT planarMPC_V1[60];
planarMPC_FLOAT planarMPC_Yd1[21];
planarMPC_FLOAT planarMPC_Ld1[21];
planarMPC_FLOAT planarMPC_yy1[6];
planarMPC_FLOAT planarMPC_bmy1[6];
planarMPC_FLOAT planarMPC_c1[6] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
int planarMPC_lbIdx1[10] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
planarMPC_FLOAT* planarMPC_llb1 = planarMPC_l + 20;
planarMPC_FLOAT* planarMPC_slb1 = planarMPC_s + 20;
planarMPC_FLOAT* planarMPC_llbbyslb1 = planarMPC_lbys + 20;
planarMPC_FLOAT planarMPC_rilb1[10];
planarMPC_FLOAT* planarMPC_dllbaff1 = planarMPC_dl_aff + 20;
planarMPC_FLOAT* planarMPC_dslbaff1 = planarMPC_ds_aff + 20;
planarMPC_FLOAT* planarMPC_dllbcc1 = planarMPC_dl_cc + 20;
planarMPC_FLOAT* planarMPC_dslbcc1 = planarMPC_ds_cc + 20;
planarMPC_FLOAT* planarMPC_ccrhsl1 = planarMPC_ccrhs + 20;
int planarMPC_ubIdx1[10] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
planarMPC_FLOAT* planarMPC_lub1 = planarMPC_l + 30;
planarMPC_FLOAT* planarMPC_sub1 = planarMPC_s + 30;
planarMPC_FLOAT* planarMPC_lubbysub1 = planarMPC_lbys + 30;
planarMPC_FLOAT planarMPC_riub1[10];
planarMPC_FLOAT* planarMPC_dlubaff1 = planarMPC_dl_aff + 30;
planarMPC_FLOAT* planarMPC_dsubaff1 = planarMPC_ds_aff + 30;
planarMPC_FLOAT* planarMPC_dlubcc1 = planarMPC_dl_cc + 30;
planarMPC_FLOAT* planarMPC_dsubcc1 = planarMPC_ds_cc + 30;
planarMPC_FLOAT* planarMPC_ccrhsub1 = planarMPC_ccrhs + 30;
planarMPC_FLOAT planarMPC_Phi1[10];
planarMPC_FLOAT planarMPC_D1[10] = {-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000};
planarMPC_FLOAT planarMPC_W1[10];
planarMPC_FLOAT planarMPC_Ysd1[36];
planarMPC_FLOAT planarMPC_Lsd1[36];
planarMPC_FLOAT* planarMPC_z2 = planarMPC_z + 20;
planarMPC_FLOAT* planarMPC_dzaff2 = planarMPC_dz_aff + 20;
planarMPC_FLOAT* planarMPC_dzcc2 = planarMPC_dz_cc + 20;
planarMPC_FLOAT* planarMPC_rd2 = planarMPC_rd + 20;
planarMPC_FLOAT planarMPC_Lbyrd2[10];
planarMPC_FLOAT* planarMPC_grad_cost2 = planarMPC_grad_cost + 20;
planarMPC_FLOAT* planarMPC_grad_eq2 = planarMPC_grad_eq + 20;
planarMPC_FLOAT* planarMPC_grad_ineq2 = planarMPC_grad_ineq + 20;
planarMPC_FLOAT planarMPC_ctv2[10];
planarMPC_FLOAT* planarMPC_v2 = planarMPC_v + 12;
planarMPC_FLOAT planarMPC_re2[6];
planarMPC_FLOAT planarMPC_beta2[6];
planarMPC_FLOAT planarMPC_betacc2[6];
planarMPC_FLOAT* planarMPC_dvaff2 = planarMPC_dv_aff + 12;
planarMPC_FLOAT* planarMPC_dvcc2 = planarMPC_dv_cc + 12;
planarMPC_FLOAT planarMPC_V2[60];
planarMPC_FLOAT planarMPC_Yd2[21];
planarMPC_FLOAT planarMPC_Ld2[21];
planarMPC_FLOAT planarMPC_yy2[6];
planarMPC_FLOAT planarMPC_bmy2[6];
planarMPC_FLOAT planarMPC_c2[6] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
int planarMPC_lbIdx2[10] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
planarMPC_FLOAT* planarMPC_llb2 = planarMPC_l + 40;
planarMPC_FLOAT* planarMPC_slb2 = planarMPC_s + 40;
planarMPC_FLOAT* planarMPC_llbbyslb2 = planarMPC_lbys + 40;
planarMPC_FLOAT planarMPC_rilb2[10];
planarMPC_FLOAT* planarMPC_dllbaff2 = planarMPC_dl_aff + 40;
planarMPC_FLOAT* planarMPC_dslbaff2 = planarMPC_ds_aff + 40;
planarMPC_FLOAT* planarMPC_dllbcc2 = planarMPC_dl_cc + 40;
planarMPC_FLOAT* planarMPC_dslbcc2 = planarMPC_ds_cc + 40;
planarMPC_FLOAT* planarMPC_ccrhsl2 = planarMPC_ccrhs + 40;
int planarMPC_ubIdx2[10] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
planarMPC_FLOAT* planarMPC_lub2 = planarMPC_l + 50;
planarMPC_FLOAT* planarMPC_sub2 = planarMPC_s + 50;
planarMPC_FLOAT* planarMPC_lubbysub2 = planarMPC_lbys + 50;
planarMPC_FLOAT planarMPC_riub2[10];
planarMPC_FLOAT* planarMPC_dlubaff2 = planarMPC_dl_aff + 50;
planarMPC_FLOAT* planarMPC_dsubaff2 = planarMPC_ds_aff + 50;
planarMPC_FLOAT* planarMPC_dlubcc2 = planarMPC_dl_cc + 50;
planarMPC_FLOAT* planarMPC_dsubcc2 = planarMPC_ds_cc + 50;
planarMPC_FLOAT* planarMPC_ccrhsub2 = planarMPC_ccrhs + 50;
planarMPC_FLOAT planarMPC_Phi2[10];
planarMPC_FLOAT planarMPC_W2[10];
planarMPC_FLOAT planarMPC_Ysd2[36];
planarMPC_FLOAT planarMPC_Lsd2[36];
planarMPC_FLOAT* planarMPC_z3 = planarMPC_z + 30;
planarMPC_FLOAT* planarMPC_dzaff3 = planarMPC_dz_aff + 30;
planarMPC_FLOAT* planarMPC_dzcc3 = planarMPC_dz_cc + 30;
planarMPC_FLOAT* planarMPC_rd3 = planarMPC_rd + 30;
planarMPC_FLOAT planarMPC_Lbyrd3[10];
planarMPC_FLOAT* planarMPC_grad_cost3 = planarMPC_grad_cost + 30;
planarMPC_FLOAT* planarMPC_grad_eq3 = planarMPC_grad_eq + 30;
planarMPC_FLOAT* planarMPC_grad_ineq3 = planarMPC_grad_ineq + 30;
planarMPC_FLOAT planarMPC_ctv3[10];
planarMPC_FLOAT* planarMPC_v3 = planarMPC_v + 18;
planarMPC_FLOAT planarMPC_re3[6];
planarMPC_FLOAT planarMPC_beta3[6];
planarMPC_FLOAT planarMPC_betacc3[6];
planarMPC_FLOAT* planarMPC_dvaff3 = planarMPC_dv_aff + 18;
planarMPC_FLOAT* planarMPC_dvcc3 = planarMPC_dv_cc + 18;
planarMPC_FLOAT planarMPC_V3[60];
planarMPC_FLOAT planarMPC_Yd3[21];
planarMPC_FLOAT planarMPC_Ld3[21];
planarMPC_FLOAT planarMPC_yy3[6];
planarMPC_FLOAT planarMPC_bmy3[6];
planarMPC_FLOAT planarMPC_c3[6] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
int planarMPC_lbIdx3[10] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
planarMPC_FLOAT* planarMPC_llb3 = planarMPC_l + 60;
planarMPC_FLOAT* planarMPC_slb3 = planarMPC_s + 60;
planarMPC_FLOAT* planarMPC_llbbyslb3 = planarMPC_lbys + 60;
planarMPC_FLOAT planarMPC_rilb3[10];
planarMPC_FLOAT* planarMPC_dllbaff3 = planarMPC_dl_aff + 60;
planarMPC_FLOAT* planarMPC_dslbaff3 = planarMPC_ds_aff + 60;
planarMPC_FLOAT* planarMPC_dllbcc3 = planarMPC_dl_cc + 60;
planarMPC_FLOAT* planarMPC_dslbcc3 = planarMPC_ds_cc + 60;
planarMPC_FLOAT* planarMPC_ccrhsl3 = planarMPC_ccrhs + 60;
int planarMPC_ubIdx3[10] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
planarMPC_FLOAT* planarMPC_lub3 = planarMPC_l + 70;
planarMPC_FLOAT* planarMPC_sub3 = planarMPC_s + 70;
planarMPC_FLOAT* planarMPC_lubbysub3 = planarMPC_lbys + 70;
planarMPC_FLOAT planarMPC_riub3[10];
planarMPC_FLOAT* planarMPC_dlubaff3 = planarMPC_dl_aff + 70;
planarMPC_FLOAT* planarMPC_dsubaff3 = planarMPC_ds_aff + 70;
planarMPC_FLOAT* planarMPC_dlubcc3 = planarMPC_dl_cc + 70;
planarMPC_FLOAT* planarMPC_dsubcc3 = planarMPC_ds_cc + 70;
planarMPC_FLOAT* planarMPC_ccrhsub3 = planarMPC_ccrhs + 70;
planarMPC_FLOAT planarMPC_Phi3[10];
planarMPC_FLOAT planarMPC_W3[10];
planarMPC_FLOAT planarMPC_Ysd3[36];
planarMPC_FLOAT planarMPC_Lsd3[36];
planarMPC_FLOAT* planarMPC_z4 = planarMPC_z + 40;
planarMPC_FLOAT* planarMPC_dzaff4 = planarMPC_dz_aff + 40;
planarMPC_FLOAT* planarMPC_dzcc4 = planarMPC_dz_cc + 40;
planarMPC_FLOAT* planarMPC_rd4 = planarMPC_rd + 40;
planarMPC_FLOAT planarMPC_Lbyrd4[10];
planarMPC_FLOAT* planarMPC_grad_cost4 = planarMPC_grad_cost + 40;
planarMPC_FLOAT* planarMPC_grad_eq4 = planarMPC_grad_eq + 40;
planarMPC_FLOAT* planarMPC_grad_ineq4 = planarMPC_grad_ineq + 40;
planarMPC_FLOAT planarMPC_ctv4[10];
planarMPC_FLOAT* planarMPC_v4 = planarMPC_v + 24;
planarMPC_FLOAT planarMPC_re4[6];
planarMPC_FLOAT planarMPC_beta4[6];
planarMPC_FLOAT planarMPC_betacc4[6];
planarMPC_FLOAT* planarMPC_dvaff4 = planarMPC_dv_aff + 24;
planarMPC_FLOAT* planarMPC_dvcc4 = planarMPC_dv_cc + 24;
planarMPC_FLOAT planarMPC_V4[60];
planarMPC_FLOAT planarMPC_Yd4[21];
planarMPC_FLOAT planarMPC_Ld4[21];
planarMPC_FLOAT planarMPC_yy4[6];
planarMPC_FLOAT planarMPC_bmy4[6];
planarMPC_FLOAT planarMPC_c4[6] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
int planarMPC_lbIdx4[10] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
planarMPC_FLOAT* planarMPC_llb4 = planarMPC_l + 80;
planarMPC_FLOAT* planarMPC_slb4 = planarMPC_s + 80;
planarMPC_FLOAT* planarMPC_llbbyslb4 = planarMPC_lbys + 80;
planarMPC_FLOAT planarMPC_rilb4[10];
planarMPC_FLOAT* planarMPC_dllbaff4 = planarMPC_dl_aff + 80;
planarMPC_FLOAT* planarMPC_dslbaff4 = planarMPC_ds_aff + 80;
planarMPC_FLOAT* planarMPC_dllbcc4 = planarMPC_dl_cc + 80;
planarMPC_FLOAT* planarMPC_dslbcc4 = planarMPC_ds_cc + 80;
planarMPC_FLOAT* planarMPC_ccrhsl4 = planarMPC_ccrhs + 80;
int planarMPC_ubIdx4[10] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
planarMPC_FLOAT* planarMPC_lub4 = planarMPC_l + 90;
planarMPC_FLOAT* planarMPC_sub4 = planarMPC_s + 90;
planarMPC_FLOAT* planarMPC_lubbysub4 = planarMPC_lbys + 90;
planarMPC_FLOAT planarMPC_riub4[10];
planarMPC_FLOAT* planarMPC_dlubaff4 = planarMPC_dl_aff + 90;
planarMPC_FLOAT* planarMPC_dsubaff4 = planarMPC_ds_aff + 90;
planarMPC_FLOAT* planarMPC_dlubcc4 = planarMPC_dl_cc + 90;
planarMPC_FLOAT* planarMPC_dsubcc4 = planarMPC_ds_cc + 90;
planarMPC_FLOAT* planarMPC_ccrhsub4 = planarMPC_ccrhs + 90;
planarMPC_FLOAT planarMPC_Phi4[10];
planarMPC_FLOAT planarMPC_W4[10];
planarMPC_FLOAT planarMPC_Ysd4[36];
planarMPC_FLOAT planarMPC_Lsd4[36];
planarMPC_FLOAT* planarMPC_z5 = planarMPC_z + 50;
planarMPC_FLOAT* planarMPC_dzaff5 = planarMPC_dz_aff + 50;
planarMPC_FLOAT* planarMPC_dzcc5 = planarMPC_dz_cc + 50;
planarMPC_FLOAT* planarMPC_rd5 = planarMPC_rd + 50;
planarMPC_FLOAT planarMPC_Lbyrd5[10];
planarMPC_FLOAT* planarMPC_grad_cost5 = planarMPC_grad_cost + 50;
planarMPC_FLOAT* planarMPC_grad_eq5 = planarMPC_grad_eq + 50;
planarMPC_FLOAT* planarMPC_grad_ineq5 = planarMPC_grad_ineq + 50;
planarMPC_FLOAT planarMPC_ctv5[10];
planarMPC_FLOAT* planarMPC_v5 = planarMPC_v + 30;
planarMPC_FLOAT planarMPC_re5[6];
planarMPC_FLOAT planarMPC_beta5[6];
planarMPC_FLOAT planarMPC_betacc5[6];
planarMPC_FLOAT* planarMPC_dvaff5 = planarMPC_dv_aff + 30;
planarMPC_FLOAT* planarMPC_dvcc5 = planarMPC_dv_cc + 30;
planarMPC_FLOAT planarMPC_V5[60];
planarMPC_FLOAT planarMPC_Yd5[21];
planarMPC_FLOAT planarMPC_Ld5[21];
planarMPC_FLOAT planarMPC_yy5[6];
planarMPC_FLOAT planarMPC_bmy5[6];
planarMPC_FLOAT planarMPC_c5[6] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
int planarMPC_lbIdx5[10] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
planarMPC_FLOAT* planarMPC_llb5 = planarMPC_l + 100;
planarMPC_FLOAT* planarMPC_slb5 = planarMPC_s + 100;
planarMPC_FLOAT* planarMPC_llbbyslb5 = planarMPC_lbys + 100;
planarMPC_FLOAT planarMPC_rilb5[10];
planarMPC_FLOAT* planarMPC_dllbaff5 = planarMPC_dl_aff + 100;
planarMPC_FLOAT* planarMPC_dslbaff5 = planarMPC_ds_aff + 100;
planarMPC_FLOAT* planarMPC_dllbcc5 = planarMPC_dl_cc + 100;
planarMPC_FLOAT* planarMPC_dslbcc5 = planarMPC_ds_cc + 100;
planarMPC_FLOAT* planarMPC_ccrhsl5 = planarMPC_ccrhs + 100;
int planarMPC_ubIdx5[10] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
planarMPC_FLOAT* planarMPC_lub5 = planarMPC_l + 110;
planarMPC_FLOAT* planarMPC_sub5 = planarMPC_s + 110;
planarMPC_FLOAT* planarMPC_lubbysub5 = planarMPC_lbys + 110;
planarMPC_FLOAT planarMPC_riub5[10];
planarMPC_FLOAT* planarMPC_dlubaff5 = planarMPC_dl_aff + 110;
planarMPC_FLOAT* planarMPC_dsubaff5 = planarMPC_ds_aff + 110;
planarMPC_FLOAT* planarMPC_dlubcc5 = planarMPC_dl_cc + 110;
planarMPC_FLOAT* planarMPC_dsubcc5 = planarMPC_ds_cc + 110;
planarMPC_FLOAT* planarMPC_ccrhsub5 = planarMPC_ccrhs + 110;
planarMPC_FLOAT planarMPC_Phi5[10];
planarMPC_FLOAT planarMPC_W5[10];
planarMPC_FLOAT planarMPC_Ysd5[36];
planarMPC_FLOAT planarMPC_Lsd5[36];
planarMPC_FLOAT* planarMPC_z6 = planarMPC_z + 60;
planarMPC_FLOAT* planarMPC_dzaff6 = planarMPC_dz_aff + 60;
planarMPC_FLOAT* planarMPC_dzcc6 = planarMPC_dz_cc + 60;
planarMPC_FLOAT* planarMPC_rd6 = planarMPC_rd + 60;
planarMPC_FLOAT planarMPC_Lbyrd6[10];
planarMPC_FLOAT* planarMPC_grad_cost6 = planarMPC_grad_cost + 60;
planarMPC_FLOAT* planarMPC_grad_eq6 = planarMPC_grad_eq + 60;
planarMPC_FLOAT* planarMPC_grad_ineq6 = planarMPC_grad_ineq + 60;
planarMPC_FLOAT planarMPC_ctv6[10];
planarMPC_FLOAT* planarMPC_v6 = planarMPC_v + 36;
planarMPC_FLOAT planarMPC_re6[6];
planarMPC_FLOAT planarMPC_beta6[6];
planarMPC_FLOAT planarMPC_betacc6[6];
planarMPC_FLOAT* planarMPC_dvaff6 = planarMPC_dv_aff + 36;
planarMPC_FLOAT* planarMPC_dvcc6 = planarMPC_dv_cc + 36;
planarMPC_FLOAT planarMPC_V6[60];
planarMPC_FLOAT planarMPC_Yd6[21];
planarMPC_FLOAT planarMPC_Ld6[21];
planarMPC_FLOAT planarMPC_yy6[6];
planarMPC_FLOAT planarMPC_bmy6[6];
planarMPC_FLOAT planarMPC_c6[6] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
int planarMPC_lbIdx6[10] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
planarMPC_FLOAT* planarMPC_llb6 = planarMPC_l + 120;
planarMPC_FLOAT* planarMPC_slb6 = planarMPC_s + 120;
planarMPC_FLOAT* planarMPC_llbbyslb6 = planarMPC_lbys + 120;
planarMPC_FLOAT planarMPC_rilb6[10];
planarMPC_FLOAT* planarMPC_dllbaff6 = planarMPC_dl_aff + 120;
planarMPC_FLOAT* planarMPC_dslbaff6 = planarMPC_ds_aff + 120;
planarMPC_FLOAT* planarMPC_dllbcc6 = planarMPC_dl_cc + 120;
planarMPC_FLOAT* planarMPC_dslbcc6 = planarMPC_ds_cc + 120;
planarMPC_FLOAT* planarMPC_ccrhsl6 = planarMPC_ccrhs + 120;
int planarMPC_ubIdx6[10] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
planarMPC_FLOAT* planarMPC_lub6 = planarMPC_l + 130;
planarMPC_FLOAT* planarMPC_sub6 = planarMPC_s + 130;
planarMPC_FLOAT* planarMPC_lubbysub6 = planarMPC_lbys + 130;
planarMPC_FLOAT planarMPC_riub6[10];
planarMPC_FLOAT* planarMPC_dlubaff6 = planarMPC_dl_aff + 130;
planarMPC_FLOAT* planarMPC_dsubaff6 = planarMPC_ds_aff + 130;
planarMPC_FLOAT* planarMPC_dlubcc6 = planarMPC_dl_cc + 130;
planarMPC_FLOAT* planarMPC_dsubcc6 = planarMPC_ds_cc + 130;
planarMPC_FLOAT* planarMPC_ccrhsub6 = planarMPC_ccrhs + 130;
planarMPC_FLOAT planarMPC_Phi6[10];
planarMPC_FLOAT planarMPC_W6[10];
planarMPC_FLOAT planarMPC_Ysd6[36];
planarMPC_FLOAT planarMPC_Lsd6[36];
planarMPC_FLOAT* planarMPC_z7 = planarMPC_z + 70;
planarMPC_FLOAT* planarMPC_dzaff7 = planarMPC_dz_aff + 70;
planarMPC_FLOAT* planarMPC_dzcc7 = planarMPC_dz_cc + 70;
planarMPC_FLOAT* planarMPC_rd7 = planarMPC_rd + 70;
planarMPC_FLOAT planarMPC_Lbyrd7[10];
planarMPC_FLOAT* planarMPC_grad_cost7 = planarMPC_grad_cost + 70;
planarMPC_FLOAT* planarMPC_grad_eq7 = planarMPC_grad_eq + 70;
planarMPC_FLOAT* planarMPC_grad_ineq7 = planarMPC_grad_ineq + 70;
planarMPC_FLOAT planarMPC_ctv7[10];
planarMPC_FLOAT* planarMPC_v7 = planarMPC_v + 42;
planarMPC_FLOAT planarMPC_re7[6];
planarMPC_FLOAT planarMPC_beta7[6];
planarMPC_FLOAT planarMPC_betacc7[6];
planarMPC_FLOAT* planarMPC_dvaff7 = planarMPC_dv_aff + 42;
planarMPC_FLOAT* planarMPC_dvcc7 = planarMPC_dv_cc + 42;
planarMPC_FLOAT planarMPC_V7[60];
planarMPC_FLOAT planarMPC_Yd7[21];
planarMPC_FLOAT planarMPC_Ld7[21];
planarMPC_FLOAT planarMPC_yy7[6];
planarMPC_FLOAT planarMPC_bmy7[6];
planarMPC_FLOAT planarMPC_c7[6] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
int planarMPC_lbIdx7[10] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
planarMPC_FLOAT* planarMPC_llb7 = planarMPC_l + 140;
planarMPC_FLOAT* planarMPC_slb7 = planarMPC_s + 140;
planarMPC_FLOAT* planarMPC_llbbyslb7 = planarMPC_lbys + 140;
planarMPC_FLOAT planarMPC_rilb7[10];
planarMPC_FLOAT* planarMPC_dllbaff7 = planarMPC_dl_aff + 140;
planarMPC_FLOAT* planarMPC_dslbaff7 = planarMPC_ds_aff + 140;
planarMPC_FLOAT* planarMPC_dllbcc7 = planarMPC_dl_cc + 140;
planarMPC_FLOAT* planarMPC_dslbcc7 = planarMPC_ds_cc + 140;
planarMPC_FLOAT* planarMPC_ccrhsl7 = planarMPC_ccrhs + 140;
int planarMPC_ubIdx7[10] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
planarMPC_FLOAT* planarMPC_lub7 = planarMPC_l + 150;
planarMPC_FLOAT* planarMPC_sub7 = planarMPC_s + 150;
planarMPC_FLOAT* planarMPC_lubbysub7 = planarMPC_lbys + 150;
planarMPC_FLOAT planarMPC_riub7[10];
planarMPC_FLOAT* planarMPC_dlubaff7 = planarMPC_dl_aff + 150;
planarMPC_FLOAT* planarMPC_dsubaff7 = planarMPC_ds_aff + 150;
planarMPC_FLOAT* planarMPC_dlubcc7 = planarMPC_dl_cc + 150;
planarMPC_FLOAT* planarMPC_dsubcc7 = planarMPC_ds_cc + 150;
planarMPC_FLOAT* planarMPC_ccrhsub7 = planarMPC_ccrhs + 150;
planarMPC_FLOAT planarMPC_Phi7[10];
planarMPC_FLOAT planarMPC_W7[10];
planarMPC_FLOAT planarMPC_Ysd7[36];
planarMPC_FLOAT planarMPC_Lsd7[36];
planarMPC_FLOAT* planarMPC_z8 = planarMPC_z + 80;
planarMPC_FLOAT* planarMPC_dzaff8 = planarMPC_dz_aff + 80;
planarMPC_FLOAT* planarMPC_dzcc8 = planarMPC_dz_cc + 80;
planarMPC_FLOAT* planarMPC_rd8 = planarMPC_rd + 80;
planarMPC_FLOAT planarMPC_Lbyrd8[10];
planarMPC_FLOAT* planarMPC_grad_cost8 = planarMPC_grad_cost + 80;
planarMPC_FLOAT* planarMPC_grad_eq8 = planarMPC_grad_eq + 80;
planarMPC_FLOAT* planarMPC_grad_ineq8 = planarMPC_grad_ineq + 80;
planarMPC_FLOAT planarMPC_ctv8[10];
planarMPC_FLOAT* planarMPC_v8 = planarMPC_v + 48;
planarMPC_FLOAT planarMPC_re8[6];
planarMPC_FLOAT planarMPC_beta8[6];
planarMPC_FLOAT planarMPC_betacc8[6];
planarMPC_FLOAT* planarMPC_dvaff8 = planarMPC_dv_aff + 48;
planarMPC_FLOAT* planarMPC_dvcc8 = planarMPC_dv_cc + 48;
planarMPC_FLOAT planarMPC_V8[60];
planarMPC_FLOAT planarMPC_Yd8[21];
planarMPC_FLOAT planarMPC_Ld8[21];
planarMPC_FLOAT planarMPC_yy8[6];
planarMPC_FLOAT planarMPC_bmy8[6];
planarMPC_FLOAT planarMPC_c8[6] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
int planarMPC_lbIdx8[10] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
planarMPC_FLOAT* planarMPC_llb8 = planarMPC_l + 160;
planarMPC_FLOAT* planarMPC_slb8 = planarMPC_s + 160;
planarMPC_FLOAT* planarMPC_llbbyslb8 = planarMPC_lbys + 160;
planarMPC_FLOAT planarMPC_rilb8[10];
planarMPC_FLOAT* planarMPC_dllbaff8 = planarMPC_dl_aff + 160;
planarMPC_FLOAT* planarMPC_dslbaff8 = planarMPC_ds_aff + 160;
planarMPC_FLOAT* planarMPC_dllbcc8 = planarMPC_dl_cc + 160;
planarMPC_FLOAT* planarMPC_dslbcc8 = planarMPC_ds_cc + 160;
planarMPC_FLOAT* planarMPC_ccrhsl8 = planarMPC_ccrhs + 160;
int planarMPC_ubIdx8[10] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
planarMPC_FLOAT* planarMPC_lub8 = planarMPC_l + 170;
planarMPC_FLOAT* planarMPC_sub8 = planarMPC_s + 170;
planarMPC_FLOAT* planarMPC_lubbysub8 = planarMPC_lbys + 170;
planarMPC_FLOAT planarMPC_riub8[10];
planarMPC_FLOAT* planarMPC_dlubaff8 = planarMPC_dl_aff + 170;
planarMPC_FLOAT* planarMPC_dsubaff8 = planarMPC_ds_aff + 170;
planarMPC_FLOAT* planarMPC_dlubcc8 = planarMPC_dl_cc + 170;
planarMPC_FLOAT* planarMPC_dsubcc8 = planarMPC_ds_cc + 170;
planarMPC_FLOAT* planarMPC_ccrhsub8 = planarMPC_ccrhs + 170;
planarMPC_FLOAT planarMPC_Phi8[10];
planarMPC_FLOAT planarMPC_W8[10];
planarMPC_FLOAT planarMPC_Ysd8[36];
planarMPC_FLOAT planarMPC_Lsd8[36];
planarMPC_FLOAT* planarMPC_z9 = planarMPC_z + 90;
planarMPC_FLOAT* planarMPC_dzaff9 = planarMPC_dz_aff + 90;
planarMPC_FLOAT* planarMPC_dzcc9 = planarMPC_dz_cc + 90;
planarMPC_FLOAT* planarMPC_rd9 = planarMPC_rd + 90;
planarMPC_FLOAT planarMPC_Lbyrd9[6];
planarMPC_FLOAT* planarMPC_grad_cost9 = planarMPC_grad_cost + 90;
planarMPC_FLOAT* planarMPC_grad_eq9 = planarMPC_grad_eq + 90;
planarMPC_FLOAT* planarMPC_grad_ineq9 = planarMPC_grad_ineq + 90;
planarMPC_FLOAT planarMPC_ctv9[6];
planarMPC_FLOAT* planarMPC_v9 = planarMPC_v + 54;
planarMPC_FLOAT planarMPC_re9[6];
planarMPC_FLOAT planarMPC_beta9[6];
planarMPC_FLOAT planarMPC_betacc9[6];
planarMPC_FLOAT* planarMPC_dvaff9 = planarMPC_dv_aff + 54;
planarMPC_FLOAT* planarMPC_dvcc9 = planarMPC_dv_cc + 54;
planarMPC_FLOAT planarMPC_V9[36];
planarMPC_FLOAT planarMPC_Yd9[21];
planarMPC_FLOAT planarMPC_Ld9[21];
planarMPC_FLOAT planarMPC_yy9[6];
planarMPC_FLOAT planarMPC_bmy9[6];
planarMPC_FLOAT planarMPC_c9[6] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
int planarMPC_lbIdx9[6] = {0, 1, 2, 3, 4, 5};
planarMPC_FLOAT* planarMPC_llb9 = planarMPC_l + 180;
planarMPC_FLOAT* planarMPC_slb9 = planarMPC_s + 180;
planarMPC_FLOAT* planarMPC_llbbyslb9 = planarMPC_lbys + 180;
planarMPC_FLOAT planarMPC_rilb9[6];
planarMPC_FLOAT* planarMPC_dllbaff9 = planarMPC_dl_aff + 180;
planarMPC_FLOAT* planarMPC_dslbaff9 = planarMPC_ds_aff + 180;
planarMPC_FLOAT* planarMPC_dllbcc9 = planarMPC_dl_cc + 180;
planarMPC_FLOAT* planarMPC_dslbcc9 = planarMPC_ds_cc + 180;
planarMPC_FLOAT* planarMPC_ccrhsl9 = planarMPC_ccrhs + 180;
int planarMPC_ubIdx9[6] = {0, 1, 2, 3, 4, 5};
planarMPC_FLOAT* planarMPC_lub9 = planarMPC_l + 186;
planarMPC_FLOAT* planarMPC_sub9 = planarMPC_s + 186;
planarMPC_FLOAT* planarMPC_lubbysub9 = planarMPC_lbys + 186;
planarMPC_FLOAT planarMPC_riub9[6];
planarMPC_FLOAT* planarMPC_dlubaff9 = planarMPC_dl_aff + 186;
planarMPC_FLOAT* planarMPC_dsubaff9 = planarMPC_ds_aff + 186;
planarMPC_FLOAT* planarMPC_dlubcc9 = planarMPC_dl_cc + 186;
planarMPC_FLOAT* planarMPC_dsubcc9 = planarMPC_ds_cc + 186;
planarMPC_FLOAT* planarMPC_ccrhsub9 = planarMPC_ccrhs + 186;
planarMPC_FLOAT planarMPC_Phi9[6];
planarMPC_FLOAT planarMPC_D9[6] = {-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000};
planarMPC_FLOAT planarMPC_W9[6];
planarMPC_FLOAT planarMPC_Ysd9[36];
planarMPC_FLOAT planarMPC_Lsd9[36];
planarMPC_FLOAT musigma;
planarMPC_FLOAT sigma_3rdroot;
planarMPC_FLOAT planarMPC_Diag1_0[10];
planarMPC_FLOAT planarMPC_Diag2_0[10];
planarMPC_FLOAT planarMPC_L_0[45];




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
planarMPC_LA_INITIALIZEVECTOR_96(planarMPC_z, 0);
planarMPC_LA_INITIALIZEVECTOR_60(planarMPC_v, 1);
planarMPC_LA_INITIALIZEVECTOR_192(planarMPC_l, 1);
planarMPC_LA_INITIALIZEVECTOR_192(planarMPC_s, 1);
info->mu = 0;
planarMPC_LA_DOTACC_192(planarMPC_l, planarMPC_s, &info->mu);
info->mu /= 192;
while( 1 ){
info->pobj = 0;
planarMPC_LA_DIAG_QUADFCN_10(params->H1, params->f1, planarMPC_z0, planarMPC_grad_cost0, &info->pobj);
planarMPC_LA_DIAG_QUADFCN_10(params->H2, params->f2, planarMPC_z1, planarMPC_grad_cost1, &info->pobj);
planarMPC_LA_DIAG_QUADFCN_10(params->H3, params->f3, planarMPC_z2, planarMPC_grad_cost2, &info->pobj);
planarMPC_LA_DIAG_QUADFCN_10(params->H4, params->f4, planarMPC_z3, planarMPC_grad_cost3, &info->pobj);
planarMPC_LA_DIAG_QUADFCN_10(params->H5, params->f5, planarMPC_z4, planarMPC_grad_cost4, &info->pobj);
planarMPC_LA_DIAG_QUADFCN_10(params->H6, params->f6, planarMPC_z5, planarMPC_grad_cost5, &info->pobj);
planarMPC_LA_DIAG_QUADFCN_10(params->H7, params->f7, planarMPC_z6, planarMPC_grad_cost6, &info->pobj);
planarMPC_LA_DIAG_QUADFCN_10(params->H8, params->f8, planarMPC_z7, planarMPC_grad_cost7, &info->pobj);
planarMPC_LA_DIAG_QUADFCN_10(params->H9, params->f9, planarMPC_z8, planarMPC_grad_cost8, &info->pobj);
planarMPC_LA_DIAG_QUADFCN_6(params->H10, params->f10, planarMPC_z9, planarMPC_grad_cost9, &info->pobj);
info->res_eq = 0;
info->dgap = 0;
planarMPC_LA_DIAGZERO_MVMSUB6_6(planarMPC_D0, planarMPC_z0, params->c1, planarMPC_v0, planarMPC_re0, &info->dgap, &info->res_eq);
planarMPC_LA_DENSE_DIAGZERO_MVMSUB3_6_10_10(planarMPC_C0, planarMPC_z0, planarMPC_D1, planarMPC_z1, planarMPC_c1, planarMPC_v1, planarMPC_re1, &info->dgap, &info->res_eq);
planarMPC_LA_DENSE_DIAGZERO_MVMSUB3_6_10_10(planarMPC_C0, planarMPC_z1, planarMPC_D1, planarMPC_z2, planarMPC_c2, planarMPC_v2, planarMPC_re2, &info->dgap, &info->res_eq);
planarMPC_LA_DENSE_DIAGZERO_MVMSUB3_6_10_10(planarMPC_C0, planarMPC_z2, planarMPC_D1, planarMPC_z3, planarMPC_c3, planarMPC_v3, planarMPC_re3, &info->dgap, &info->res_eq);
planarMPC_LA_DENSE_DIAGZERO_MVMSUB3_6_10_10(planarMPC_C0, planarMPC_z3, planarMPC_D1, planarMPC_z4, planarMPC_c4, planarMPC_v4, planarMPC_re4, &info->dgap, &info->res_eq);
planarMPC_LA_DENSE_DIAGZERO_MVMSUB3_6_10_10(planarMPC_C0, planarMPC_z4, planarMPC_D1, planarMPC_z5, planarMPC_c5, planarMPC_v5, planarMPC_re5, &info->dgap, &info->res_eq);
planarMPC_LA_DENSE_DIAGZERO_MVMSUB3_6_10_10(planarMPC_C0, planarMPC_z5, planarMPC_D1, planarMPC_z6, planarMPC_c6, planarMPC_v6, planarMPC_re6, &info->dgap, &info->res_eq);
planarMPC_LA_DENSE_DIAGZERO_MVMSUB3_6_10_10(planarMPC_C0, planarMPC_z6, planarMPC_D1, planarMPC_z7, planarMPC_c7, planarMPC_v7, planarMPC_re7, &info->dgap, &info->res_eq);
planarMPC_LA_DENSE_DIAGZERO_MVMSUB3_6_10_10(planarMPC_C0, planarMPC_z7, planarMPC_D1, planarMPC_z8, planarMPC_c8, planarMPC_v8, planarMPC_re8, &info->dgap, &info->res_eq);
planarMPC_LA_DENSE_DIAGZERO_MVMSUB3_6_10_6(planarMPC_C0, planarMPC_z8, planarMPC_D9, planarMPC_z9, planarMPC_c9, planarMPC_v9, planarMPC_re9, &info->dgap, &info->res_eq);
planarMPC_LA_DENSE_DIAGZERO_MTVM2_6_10_6(planarMPC_C0, planarMPC_v1, planarMPC_D0, planarMPC_v0, planarMPC_grad_eq0);
planarMPC_LA_DENSE_DIAGZERO_MTVM2_6_10_6(planarMPC_C0, planarMPC_v2, planarMPC_D1, planarMPC_v1, planarMPC_grad_eq1);
planarMPC_LA_DENSE_DIAGZERO_MTVM2_6_10_6(planarMPC_C0, planarMPC_v3, planarMPC_D1, planarMPC_v2, planarMPC_grad_eq2);
planarMPC_LA_DENSE_DIAGZERO_MTVM2_6_10_6(planarMPC_C0, planarMPC_v4, planarMPC_D1, planarMPC_v3, planarMPC_grad_eq3);
planarMPC_LA_DENSE_DIAGZERO_MTVM2_6_10_6(planarMPC_C0, planarMPC_v5, planarMPC_D1, planarMPC_v4, planarMPC_grad_eq4);
planarMPC_LA_DENSE_DIAGZERO_MTVM2_6_10_6(planarMPC_C0, planarMPC_v6, planarMPC_D1, planarMPC_v5, planarMPC_grad_eq5);
planarMPC_LA_DENSE_DIAGZERO_MTVM2_6_10_6(planarMPC_C0, planarMPC_v7, planarMPC_D1, planarMPC_v6, planarMPC_grad_eq6);
planarMPC_LA_DENSE_DIAGZERO_MTVM2_6_10_6(planarMPC_C0, planarMPC_v8, planarMPC_D1, planarMPC_v7, planarMPC_grad_eq7);
planarMPC_LA_DENSE_DIAGZERO_MTVM2_6_10_6(planarMPC_C0, planarMPC_v9, planarMPC_D1, planarMPC_v8, planarMPC_grad_eq8);
planarMPC_LA_DIAGZERO_MTVM_6_6(planarMPC_D9, planarMPC_v9, planarMPC_grad_eq9);
info->res_ineq = 0;
planarMPC_LA_VSUBADD3_10(params->lb1, planarMPC_z0, planarMPC_lbIdx0, planarMPC_llb0, planarMPC_slb0, planarMPC_rilb0, &info->dgap, &info->res_ineq);
planarMPC_LA_VSUBADD2_10(planarMPC_z0, planarMPC_ubIdx0, params->ub1, planarMPC_lub0, planarMPC_sub0, planarMPC_riub0, &info->dgap, &info->res_ineq);
planarMPC_LA_VSUBADD3_10(params->lb2, planarMPC_z1, planarMPC_lbIdx1, planarMPC_llb1, planarMPC_slb1, planarMPC_rilb1, &info->dgap, &info->res_ineq);
planarMPC_LA_VSUBADD2_10(planarMPC_z1, planarMPC_ubIdx1, params->ub2, planarMPC_lub1, planarMPC_sub1, planarMPC_riub1, &info->dgap, &info->res_ineq);
planarMPC_LA_VSUBADD3_10(params->lb3, planarMPC_z2, planarMPC_lbIdx2, planarMPC_llb2, planarMPC_slb2, planarMPC_rilb2, &info->dgap, &info->res_ineq);
planarMPC_LA_VSUBADD2_10(planarMPC_z2, planarMPC_ubIdx2, params->ub3, planarMPC_lub2, planarMPC_sub2, planarMPC_riub2, &info->dgap, &info->res_ineq);
planarMPC_LA_VSUBADD3_10(params->lb4, planarMPC_z3, planarMPC_lbIdx3, planarMPC_llb3, planarMPC_slb3, planarMPC_rilb3, &info->dgap, &info->res_ineq);
planarMPC_LA_VSUBADD2_10(planarMPC_z3, planarMPC_ubIdx3, params->ub4, planarMPC_lub3, planarMPC_sub3, planarMPC_riub3, &info->dgap, &info->res_ineq);
planarMPC_LA_VSUBADD3_10(params->lb5, planarMPC_z4, planarMPC_lbIdx4, planarMPC_llb4, planarMPC_slb4, planarMPC_rilb4, &info->dgap, &info->res_ineq);
planarMPC_LA_VSUBADD2_10(planarMPC_z4, planarMPC_ubIdx4, params->ub5, planarMPC_lub4, planarMPC_sub4, planarMPC_riub4, &info->dgap, &info->res_ineq);
planarMPC_LA_VSUBADD3_10(params->lb6, planarMPC_z5, planarMPC_lbIdx5, planarMPC_llb5, planarMPC_slb5, planarMPC_rilb5, &info->dgap, &info->res_ineq);
planarMPC_LA_VSUBADD2_10(planarMPC_z5, planarMPC_ubIdx5, params->ub6, planarMPC_lub5, planarMPC_sub5, planarMPC_riub5, &info->dgap, &info->res_ineq);
planarMPC_LA_VSUBADD3_10(params->lb7, planarMPC_z6, planarMPC_lbIdx6, planarMPC_llb6, planarMPC_slb6, planarMPC_rilb6, &info->dgap, &info->res_ineq);
planarMPC_LA_VSUBADD2_10(planarMPC_z6, planarMPC_ubIdx6, params->ub7, planarMPC_lub6, planarMPC_sub6, planarMPC_riub6, &info->dgap, &info->res_ineq);
planarMPC_LA_VSUBADD3_10(params->lb8, planarMPC_z7, planarMPC_lbIdx7, planarMPC_llb7, planarMPC_slb7, planarMPC_rilb7, &info->dgap, &info->res_ineq);
planarMPC_LA_VSUBADD2_10(planarMPC_z7, planarMPC_ubIdx7, params->ub8, planarMPC_lub7, planarMPC_sub7, planarMPC_riub7, &info->dgap, &info->res_ineq);
planarMPC_LA_VSUBADD3_10(params->lb9, planarMPC_z8, planarMPC_lbIdx8, planarMPC_llb8, planarMPC_slb8, planarMPC_rilb8, &info->dgap, &info->res_ineq);
planarMPC_LA_VSUBADD2_10(planarMPC_z8, planarMPC_ubIdx8, params->ub9, planarMPC_lub8, planarMPC_sub8, planarMPC_riub8, &info->dgap, &info->res_ineq);
planarMPC_LA_VSUBADD3_6(params->lb10, planarMPC_z9, planarMPC_lbIdx9, planarMPC_llb9, planarMPC_slb9, planarMPC_rilb9, &info->dgap, &info->res_ineq);
planarMPC_LA_VSUBADD2_6(planarMPC_z9, planarMPC_ubIdx9, params->ub10, planarMPC_lub9, planarMPC_sub9, planarMPC_riub9, &info->dgap, &info->res_ineq);
planarMPC_LA_INEQ_B_GRAD_10_10_10(planarMPC_lub0, planarMPC_sub0, planarMPC_riub0, planarMPC_llb0, planarMPC_slb0, planarMPC_rilb0, planarMPC_lbIdx0, planarMPC_ubIdx0, planarMPC_grad_ineq0, planarMPC_lubbysub0, planarMPC_llbbyslb0);
planarMPC_LA_INEQ_B_GRAD_10_10_10(planarMPC_lub1, planarMPC_sub1, planarMPC_riub1, planarMPC_llb1, planarMPC_slb1, planarMPC_rilb1, planarMPC_lbIdx1, planarMPC_ubIdx1, planarMPC_grad_ineq1, planarMPC_lubbysub1, planarMPC_llbbyslb1);
planarMPC_LA_INEQ_B_GRAD_10_10_10(planarMPC_lub2, planarMPC_sub2, planarMPC_riub2, planarMPC_llb2, planarMPC_slb2, planarMPC_rilb2, planarMPC_lbIdx2, planarMPC_ubIdx2, planarMPC_grad_ineq2, planarMPC_lubbysub2, planarMPC_llbbyslb2);
planarMPC_LA_INEQ_B_GRAD_10_10_10(planarMPC_lub3, planarMPC_sub3, planarMPC_riub3, planarMPC_llb3, planarMPC_slb3, planarMPC_rilb3, planarMPC_lbIdx3, planarMPC_ubIdx3, planarMPC_grad_ineq3, planarMPC_lubbysub3, planarMPC_llbbyslb3);
planarMPC_LA_INEQ_B_GRAD_10_10_10(planarMPC_lub4, planarMPC_sub4, planarMPC_riub4, planarMPC_llb4, planarMPC_slb4, planarMPC_rilb4, planarMPC_lbIdx4, planarMPC_ubIdx4, planarMPC_grad_ineq4, planarMPC_lubbysub4, planarMPC_llbbyslb4);
planarMPC_LA_INEQ_B_GRAD_10_10_10(planarMPC_lub5, planarMPC_sub5, planarMPC_riub5, planarMPC_llb5, planarMPC_slb5, planarMPC_rilb5, planarMPC_lbIdx5, planarMPC_ubIdx5, planarMPC_grad_ineq5, planarMPC_lubbysub5, planarMPC_llbbyslb5);
planarMPC_LA_INEQ_B_GRAD_10_10_10(planarMPC_lub6, planarMPC_sub6, planarMPC_riub6, planarMPC_llb6, planarMPC_slb6, planarMPC_rilb6, planarMPC_lbIdx6, planarMPC_ubIdx6, planarMPC_grad_ineq6, planarMPC_lubbysub6, planarMPC_llbbyslb6);
planarMPC_LA_INEQ_B_GRAD_10_10_10(planarMPC_lub7, planarMPC_sub7, planarMPC_riub7, planarMPC_llb7, planarMPC_slb7, planarMPC_rilb7, planarMPC_lbIdx7, planarMPC_ubIdx7, planarMPC_grad_ineq7, planarMPC_lubbysub7, planarMPC_llbbyslb7);
planarMPC_LA_INEQ_B_GRAD_10_10_10(planarMPC_lub8, planarMPC_sub8, planarMPC_riub8, planarMPC_llb8, planarMPC_slb8, planarMPC_rilb8, planarMPC_lbIdx8, planarMPC_ubIdx8, planarMPC_grad_ineq8, planarMPC_lubbysub8, planarMPC_llbbyslb8);
planarMPC_LA_INEQ_B_GRAD_6_6_6(planarMPC_lub9, planarMPC_sub9, planarMPC_riub9, planarMPC_llb9, planarMPC_slb9, planarMPC_rilb9, planarMPC_lbIdx9, planarMPC_ubIdx9, planarMPC_grad_ineq9, planarMPC_lubbysub9, planarMPC_llbbyslb9);
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
planarMPC_LA_VVADD3_96(planarMPC_grad_cost, planarMPC_grad_eq, planarMPC_grad_ineq, planarMPC_rd);
planarMPC_LA_DIAG_CHOL_ONELOOP_LBUB_10_10_10(params->H1, planarMPC_llbbyslb0, planarMPC_lbIdx0, planarMPC_lubbysub0, planarMPC_ubIdx0, planarMPC_Phi0);
planarMPC_LA_DIAG_MATRIXFORWARDSUB_6_10(planarMPC_Phi0, planarMPC_C0, planarMPC_V0);
planarMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_6_10(planarMPC_Phi0, planarMPC_D0, planarMPC_W0);
planarMPC_LA_DENSE_DIAGZERO_MMTM_6_10_6(planarMPC_W0, planarMPC_V0, planarMPC_Ysd1);
planarMPC_LA_DIAG_FORWARDSUB_10(planarMPC_Phi0, planarMPC_rd0, planarMPC_Lbyrd0);
planarMPC_LA_DIAG_CHOL_ONELOOP_LBUB_10_10_10(params->H2, planarMPC_llbbyslb1, planarMPC_lbIdx1, planarMPC_lubbysub1, planarMPC_ubIdx1, planarMPC_Phi1);
planarMPC_LA_DIAG_MATRIXFORWARDSUB_6_10(planarMPC_Phi1, planarMPC_C0, planarMPC_V1);
planarMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_6_10(planarMPC_Phi1, planarMPC_D1, planarMPC_W1);
planarMPC_LA_DENSE_DIAGZERO_MMTM_6_10_6(planarMPC_W1, planarMPC_V1, planarMPC_Ysd2);
planarMPC_LA_DIAG_FORWARDSUB_10(planarMPC_Phi1, planarMPC_rd1, planarMPC_Lbyrd1);
planarMPC_LA_DIAG_CHOL_ONELOOP_LBUB_10_10_10(params->H3, planarMPC_llbbyslb2, planarMPC_lbIdx2, planarMPC_lubbysub2, planarMPC_ubIdx2, planarMPC_Phi2);
planarMPC_LA_DIAG_MATRIXFORWARDSUB_6_10(planarMPC_Phi2, planarMPC_C0, planarMPC_V2);
planarMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_6_10(planarMPC_Phi2, planarMPC_D1, planarMPC_W2);
planarMPC_LA_DENSE_DIAGZERO_MMTM_6_10_6(planarMPC_W2, planarMPC_V2, planarMPC_Ysd3);
planarMPC_LA_DIAG_FORWARDSUB_10(planarMPC_Phi2, planarMPC_rd2, planarMPC_Lbyrd2);
planarMPC_LA_DIAG_CHOL_ONELOOP_LBUB_10_10_10(params->H4, planarMPC_llbbyslb3, planarMPC_lbIdx3, planarMPC_lubbysub3, planarMPC_ubIdx3, planarMPC_Phi3);
planarMPC_LA_DIAG_MATRIXFORWARDSUB_6_10(planarMPC_Phi3, planarMPC_C0, planarMPC_V3);
planarMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_6_10(planarMPC_Phi3, planarMPC_D1, planarMPC_W3);
planarMPC_LA_DENSE_DIAGZERO_MMTM_6_10_6(planarMPC_W3, planarMPC_V3, planarMPC_Ysd4);
planarMPC_LA_DIAG_FORWARDSUB_10(planarMPC_Phi3, planarMPC_rd3, planarMPC_Lbyrd3);
planarMPC_LA_DIAG_CHOL_ONELOOP_LBUB_10_10_10(params->H5, planarMPC_llbbyslb4, planarMPC_lbIdx4, planarMPC_lubbysub4, planarMPC_ubIdx4, planarMPC_Phi4);
planarMPC_LA_DIAG_MATRIXFORWARDSUB_6_10(planarMPC_Phi4, planarMPC_C0, planarMPC_V4);
planarMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_6_10(planarMPC_Phi4, planarMPC_D1, planarMPC_W4);
planarMPC_LA_DENSE_DIAGZERO_MMTM_6_10_6(planarMPC_W4, planarMPC_V4, planarMPC_Ysd5);
planarMPC_LA_DIAG_FORWARDSUB_10(planarMPC_Phi4, planarMPC_rd4, planarMPC_Lbyrd4);
planarMPC_LA_DIAG_CHOL_ONELOOP_LBUB_10_10_10(params->H6, planarMPC_llbbyslb5, planarMPC_lbIdx5, planarMPC_lubbysub5, planarMPC_ubIdx5, planarMPC_Phi5);
planarMPC_LA_DIAG_MATRIXFORWARDSUB_6_10(planarMPC_Phi5, planarMPC_C0, planarMPC_V5);
planarMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_6_10(planarMPC_Phi5, planarMPC_D1, planarMPC_W5);
planarMPC_LA_DENSE_DIAGZERO_MMTM_6_10_6(planarMPC_W5, planarMPC_V5, planarMPC_Ysd6);
planarMPC_LA_DIAG_FORWARDSUB_10(planarMPC_Phi5, planarMPC_rd5, planarMPC_Lbyrd5);
planarMPC_LA_DIAG_CHOL_ONELOOP_LBUB_10_10_10(params->H7, planarMPC_llbbyslb6, planarMPC_lbIdx6, planarMPC_lubbysub6, planarMPC_ubIdx6, planarMPC_Phi6);
planarMPC_LA_DIAG_MATRIXFORWARDSUB_6_10(planarMPC_Phi6, planarMPC_C0, planarMPC_V6);
planarMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_6_10(planarMPC_Phi6, planarMPC_D1, planarMPC_W6);
planarMPC_LA_DENSE_DIAGZERO_MMTM_6_10_6(planarMPC_W6, planarMPC_V6, planarMPC_Ysd7);
planarMPC_LA_DIAG_FORWARDSUB_10(planarMPC_Phi6, planarMPC_rd6, planarMPC_Lbyrd6);
planarMPC_LA_DIAG_CHOL_ONELOOP_LBUB_10_10_10(params->H8, planarMPC_llbbyslb7, planarMPC_lbIdx7, planarMPC_lubbysub7, planarMPC_ubIdx7, planarMPC_Phi7);
planarMPC_LA_DIAG_MATRIXFORWARDSUB_6_10(planarMPC_Phi7, planarMPC_C0, planarMPC_V7);
planarMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_6_10(planarMPC_Phi7, planarMPC_D1, planarMPC_W7);
planarMPC_LA_DENSE_DIAGZERO_MMTM_6_10_6(planarMPC_W7, planarMPC_V7, planarMPC_Ysd8);
planarMPC_LA_DIAG_FORWARDSUB_10(planarMPC_Phi7, planarMPC_rd7, planarMPC_Lbyrd7);
planarMPC_LA_DIAG_CHOL_ONELOOP_LBUB_10_10_10(params->H9, planarMPC_llbbyslb8, planarMPC_lbIdx8, planarMPC_lubbysub8, planarMPC_ubIdx8, planarMPC_Phi8);
planarMPC_LA_DIAG_MATRIXFORWARDSUB_6_10(planarMPC_Phi8, planarMPC_C0, planarMPC_V8);
planarMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_6_10(planarMPC_Phi8, planarMPC_D1, planarMPC_W8);
planarMPC_LA_DENSE_DIAGZERO_MMTM_6_10_6(planarMPC_W8, planarMPC_V8, planarMPC_Ysd9);
planarMPC_LA_DIAG_FORWARDSUB_10(planarMPC_Phi8, planarMPC_rd8, planarMPC_Lbyrd8);
planarMPC_LA_DIAG_CHOL_ONELOOP_LBUB_6_6_6(params->H10, planarMPC_llbbyslb9, planarMPC_lbIdx9, planarMPC_lubbysub9, planarMPC_ubIdx9, planarMPC_Phi9);
planarMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_6_6(planarMPC_Phi9, planarMPC_D9, planarMPC_W9);
planarMPC_LA_DIAG_FORWARDSUB_6(planarMPC_Phi9, planarMPC_rd9, planarMPC_Lbyrd9);
planarMPC_LA_DIAGZERO_MMT_6(planarMPC_W0, planarMPC_Yd0);
planarMPC_LA_DIAGZERO_MVMSUB7_6(planarMPC_W0, planarMPC_Lbyrd0, planarMPC_re0, planarMPC_beta0);
planarMPC_LA_DENSE_DIAGZERO_MMT2_6_10_10(planarMPC_V0, planarMPC_W1, planarMPC_Yd1);
planarMPC_LA_DENSE_DIAGZERO_2MVMSUB2_6_10_10(planarMPC_V0, planarMPC_Lbyrd0, planarMPC_W1, planarMPC_Lbyrd1, planarMPC_re1, planarMPC_beta1);
planarMPC_LA_DENSE_DIAGZERO_MMT2_6_10_10(planarMPC_V1, planarMPC_W2, planarMPC_Yd2);
planarMPC_LA_DENSE_DIAGZERO_2MVMSUB2_6_10_10(planarMPC_V1, planarMPC_Lbyrd1, planarMPC_W2, planarMPC_Lbyrd2, planarMPC_re2, planarMPC_beta2);
planarMPC_LA_DENSE_DIAGZERO_MMT2_6_10_10(planarMPC_V2, planarMPC_W3, planarMPC_Yd3);
planarMPC_LA_DENSE_DIAGZERO_2MVMSUB2_6_10_10(planarMPC_V2, planarMPC_Lbyrd2, planarMPC_W3, planarMPC_Lbyrd3, planarMPC_re3, planarMPC_beta3);
planarMPC_LA_DENSE_DIAGZERO_MMT2_6_10_10(planarMPC_V3, planarMPC_W4, planarMPC_Yd4);
planarMPC_LA_DENSE_DIAGZERO_2MVMSUB2_6_10_10(planarMPC_V3, planarMPC_Lbyrd3, planarMPC_W4, planarMPC_Lbyrd4, planarMPC_re4, planarMPC_beta4);
planarMPC_LA_DENSE_DIAGZERO_MMT2_6_10_10(planarMPC_V4, planarMPC_W5, planarMPC_Yd5);
planarMPC_LA_DENSE_DIAGZERO_2MVMSUB2_6_10_10(planarMPC_V4, planarMPC_Lbyrd4, planarMPC_W5, planarMPC_Lbyrd5, planarMPC_re5, planarMPC_beta5);
planarMPC_LA_DENSE_DIAGZERO_MMT2_6_10_10(planarMPC_V5, planarMPC_W6, planarMPC_Yd6);
planarMPC_LA_DENSE_DIAGZERO_2MVMSUB2_6_10_10(planarMPC_V5, planarMPC_Lbyrd5, planarMPC_W6, planarMPC_Lbyrd6, planarMPC_re6, planarMPC_beta6);
planarMPC_LA_DENSE_DIAGZERO_MMT2_6_10_10(planarMPC_V6, planarMPC_W7, planarMPC_Yd7);
planarMPC_LA_DENSE_DIAGZERO_2MVMSUB2_6_10_10(planarMPC_V6, planarMPC_Lbyrd6, planarMPC_W7, planarMPC_Lbyrd7, planarMPC_re7, planarMPC_beta7);
planarMPC_LA_DENSE_DIAGZERO_MMT2_6_10_10(planarMPC_V7, planarMPC_W8, planarMPC_Yd8);
planarMPC_LA_DENSE_DIAGZERO_2MVMSUB2_6_10_10(planarMPC_V7, planarMPC_Lbyrd7, planarMPC_W8, planarMPC_Lbyrd8, planarMPC_re8, planarMPC_beta8);
planarMPC_LA_DENSE_DIAGZERO_MMT2_6_10_6(planarMPC_V8, planarMPC_W9, planarMPC_Yd9);
planarMPC_LA_DENSE_DIAGZERO_2MVMSUB2_6_10_6(planarMPC_V8, planarMPC_Lbyrd8, planarMPC_W9, planarMPC_Lbyrd9, planarMPC_re9, planarMPC_beta9);
planarMPC_LA_DENSE_CHOL_6(planarMPC_Yd0, planarMPC_Ld0);
planarMPC_LA_DENSE_FORWARDSUB_6(planarMPC_Ld0, planarMPC_beta0, planarMPC_yy0);
planarMPC_LA_DENSE_MATRIXTFORWARDSUB_6_6(planarMPC_Ld0, planarMPC_Ysd1, planarMPC_Lsd1);
planarMPC_LA_DENSE_MMTSUB_6_6(planarMPC_Lsd1, planarMPC_Yd1);
planarMPC_LA_DENSE_CHOL_6(planarMPC_Yd1, planarMPC_Ld1);
planarMPC_LA_DENSE_MVMSUB1_6_6(planarMPC_Lsd1, planarMPC_yy0, planarMPC_beta1, planarMPC_bmy1);
planarMPC_LA_DENSE_FORWARDSUB_6(planarMPC_Ld1, planarMPC_bmy1, planarMPC_yy1);
planarMPC_LA_DENSE_MATRIXTFORWARDSUB_6_6(planarMPC_Ld1, planarMPC_Ysd2, planarMPC_Lsd2);
planarMPC_LA_DENSE_MMTSUB_6_6(planarMPC_Lsd2, planarMPC_Yd2);
planarMPC_LA_DENSE_CHOL_6(planarMPC_Yd2, planarMPC_Ld2);
planarMPC_LA_DENSE_MVMSUB1_6_6(planarMPC_Lsd2, planarMPC_yy1, planarMPC_beta2, planarMPC_bmy2);
planarMPC_LA_DENSE_FORWARDSUB_6(planarMPC_Ld2, planarMPC_bmy2, planarMPC_yy2);
planarMPC_LA_DENSE_MATRIXTFORWARDSUB_6_6(planarMPC_Ld2, planarMPC_Ysd3, planarMPC_Lsd3);
planarMPC_LA_DENSE_MMTSUB_6_6(planarMPC_Lsd3, planarMPC_Yd3);
planarMPC_LA_DENSE_CHOL_6(planarMPC_Yd3, planarMPC_Ld3);
planarMPC_LA_DENSE_MVMSUB1_6_6(planarMPC_Lsd3, planarMPC_yy2, planarMPC_beta3, planarMPC_bmy3);
planarMPC_LA_DENSE_FORWARDSUB_6(planarMPC_Ld3, planarMPC_bmy3, planarMPC_yy3);
planarMPC_LA_DENSE_MATRIXTFORWARDSUB_6_6(planarMPC_Ld3, planarMPC_Ysd4, planarMPC_Lsd4);
planarMPC_LA_DENSE_MMTSUB_6_6(planarMPC_Lsd4, planarMPC_Yd4);
planarMPC_LA_DENSE_CHOL_6(planarMPC_Yd4, planarMPC_Ld4);
planarMPC_LA_DENSE_MVMSUB1_6_6(planarMPC_Lsd4, planarMPC_yy3, planarMPC_beta4, planarMPC_bmy4);
planarMPC_LA_DENSE_FORWARDSUB_6(planarMPC_Ld4, planarMPC_bmy4, planarMPC_yy4);
planarMPC_LA_DENSE_MATRIXTFORWARDSUB_6_6(planarMPC_Ld4, planarMPC_Ysd5, planarMPC_Lsd5);
planarMPC_LA_DENSE_MMTSUB_6_6(planarMPC_Lsd5, planarMPC_Yd5);
planarMPC_LA_DENSE_CHOL_6(planarMPC_Yd5, planarMPC_Ld5);
planarMPC_LA_DENSE_MVMSUB1_6_6(planarMPC_Lsd5, planarMPC_yy4, planarMPC_beta5, planarMPC_bmy5);
planarMPC_LA_DENSE_FORWARDSUB_6(planarMPC_Ld5, planarMPC_bmy5, planarMPC_yy5);
planarMPC_LA_DENSE_MATRIXTFORWARDSUB_6_6(planarMPC_Ld5, planarMPC_Ysd6, planarMPC_Lsd6);
planarMPC_LA_DENSE_MMTSUB_6_6(planarMPC_Lsd6, planarMPC_Yd6);
planarMPC_LA_DENSE_CHOL_6(planarMPC_Yd6, planarMPC_Ld6);
planarMPC_LA_DENSE_MVMSUB1_6_6(planarMPC_Lsd6, planarMPC_yy5, planarMPC_beta6, planarMPC_bmy6);
planarMPC_LA_DENSE_FORWARDSUB_6(planarMPC_Ld6, planarMPC_bmy6, planarMPC_yy6);
planarMPC_LA_DENSE_MATRIXTFORWARDSUB_6_6(planarMPC_Ld6, planarMPC_Ysd7, planarMPC_Lsd7);
planarMPC_LA_DENSE_MMTSUB_6_6(planarMPC_Lsd7, planarMPC_Yd7);
planarMPC_LA_DENSE_CHOL_6(planarMPC_Yd7, planarMPC_Ld7);
planarMPC_LA_DENSE_MVMSUB1_6_6(planarMPC_Lsd7, planarMPC_yy6, planarMPC_beta7, planarMPC_bmy7);
planarMPC_LA_DENSE_FORWARDSUB_6(planarMPC_Ld7, planarMPC_bmy7, planarMPC_yy7);
planarMPC_LA_DENSE_MATRIXTFORWARDSUB_6_6(planarMPC_Ld7, planarMPC_Ysd8, planarMPC_Lsd8);
planarMPC_LA_DENSE_MMTSUB_6_6(planarMPC_Lsd8, planarMPC_Yd8);
planarMPC_LA_DENSE_CHOL_6(planarMPC_Yd8, planarMPC_Ld8);
planarMPC_LA_DENSE_MVMSUB1_6_6(planarMPC_Lsd8, planarMPC_yy7, planarMPC_beta8, planarMPC_bmy8);
planarMPC_LA_DENSE_FORWARDSUB_6(planarMPC_Ld8, planarMPC_bmy8, planarMPC_yy8);
planarMPC_LA_DENSE_MATRIXTFORWARDSUB_6_6(planarMPC_Ld8, planarMPC_Ysd9, planarMPC_Lsd9);
planarMPC_LA_DENSE_MMTSUB_6_6(planarMPC_Lsd9, planarMPC_Yd9);
planarMPC_LA_DENSE_CHOL_6(planarMPC_Yd9, planarMPC_Ld9);
planarMPC_LA_DENSE_MVMSUB1_6_6(planarMPC_Lsd9, planarMPC_yy8, planarMPC_beta9, planarMPC_bmy9);
planarMPC_LA_DENSE_FORWARDSUB_6(planarMPC_Ld9, planarMPC_bmy9, planarMPC_yy9);
planarMPC_LA_DENSE_BACKWARDSUB_6(planarMPC_Ld9, planarMPC_yy9, planarMPC_dvaff9);
planarMPC_LA_DENSE_MTVMSUB_6_6(planarMPC_Lsd9, planarMPC_dvaff9, planarMPC_yy8, planarMPC_bmy8);
planarMPC_LA_DENSE_BACKWARDSUB_6(planarMPC_Ld8, planarMPC_bmy8, planarMPC_dvaff8);
planarMPC_LA_DENSE_MTVMSUB_6_6(planarMPC_Lsd8, planarMPC_dvaff8, planarMPC_yy7, planarMPC_bmy7);
planarMPC_LA_DENSE_BACKWARDSUB_6(planarMPC_Ld7, planarMPC_bmy7, planarMPC_dvaff7);
planarMPC_LA_DENSE_MTVMSUB_6_6(planarMPC_Lsd7, planarMPC_dvaff7, planarMPC_yy6, planarMPC_bmy6);
planarMPC_LA_DENSE_BACKWARDSUB_6(planarMPC_Ld6, planarMPC_bmy6, planarMPC_dvaff6);
planarMPC_LA_DENSE_MTVMSUB_6_6(planarMPC_Lsd6, planarMPC_dvaff6, planarMPC_yy5, planarMPC_bmy5);
planarMPC_LA_DENSE_BACKWARDSUB_6(planarMPC_Ld5, planarMPC_bmy5, planarMPC_dvaff5);
planarMPC_LA_DENSE_MTVMSUB_6_6(planarMPC_Lsd5, planarMPC_dvaff5, planarMPC_yy4, planarMPC_bmy4);
planarMPC_LA_DENSE_BACKWARDSUB_6(planarMPC_Ld4, planarMPC_bmy4, planarMPC_dvaff4);
planarMPC_LA_DENSE_MTVMSUB_6_6(planarMPC_Lsd4, planarMPC_dvaff4, planarMPC_yy3, planarMPC_bmy3);
planarMPC_LA_DENSE_BACKWARDSUB_6(planarMPC_Ld3, planarMPC_bmy3, planarMPC_dvaff3);
planarMPC_LA_DENSE_MTVMSUB_6_6(planarMPC_Lsd3, planarMPC_dvaff3, planarMPC_yy2, planarMPC_bmy2);
planarMPC_LA_DENSE_BACKWARDSUB_6(planarMPC_Ld2, planarMPC_bmy2, planarMPC_dvaff2);
planarMPC_LA_DENSE_MTVMSUB_6_6(planarMPC_Lsd2, planarMPC_dvaff2, planarMPC_yy1, planarMPC_bmy1);
planarMPC_LA_DENSE_BACKWARDSUB_6(planarMPC_Ld1, planarMPC_bmy1, planarMPC_dvaff1);
planarMPC_LA_DENSE_MTVMSUB_6_6(planarMPC_Lsd1, planarMPC_dvaff1, planarMPC_yy0, planarMPC_bmy0);
planarMPC_LA_DENSE_BACKWARDSUB_6(planarMPC_Ld0, planarMPC_bmy0, planarMPC_dvaff0);
planarMPC_LA_DENSE_DIAGZERO_MTVM2_6_10_6(planarMPC_C0, planarMPC_dvaff1, planarMPC_D0, planarMPC_dvaff0, planarMPC_grad_eq0);
planarMPC_LA_DENSE_DIAGZERO_MTVM2_6_10_6(planarMPC_C0, planarMPC_dvaff2, planarMPC_D1, planarMPC_dvaff1, planarMPC_grad_eq1);
planarMPC_LA_DENSE_DIAGZERO_MTVM2_6_10_6(planarMPC_C0, planarMPC_dvaff3, planarMPC_D1, planarMPC_dvaff2, planarMPC_grad_eq2);
planarMPC_LA_DENSE_DIAGZERO_MTVM2_6_10_6(planarMPC_C0, planarMPC_dvaff4, planarMPC_D1, planarMPC_dvaff3, planarMPC_grad_eq3);
planarMPC_LA_DENSE_DIAGZERO_MTVM2_6_10_6(planarMPC_C0, planarMPC_dvaff5, planarMPC_D1, planarMPC_dvaff4, planarMPC_grad_eq4);
planarMPC_LA_DENSE_DIAGZERO_MTVM2_6_10_6(planarMPC_C0, planarMPC_dvaff6, planarMPC_D1, planarMPC_dvaff5, planarMPC_grad_eq5);
planarMPC_LA_DENSE_DIAGZERO_MTVM2_6_10_6(planarMPC_C0, planarMPC_dvaff7, planarMPC_D1, planarMPC_dvaff6, planarMPC_grad_eq6);
planarMPC_LA_DENSE_DIAGZERO_MTVM2_6_10_6(planarMPC_C0, planarMPC_dvaff8, planarMPC_D1, planarMPC_dvaff7, planarMPC_grad_eq7);
planarMPC_LA_DENSE_DIAGZERO_MTVM2_6_10_6(planarMPC_C0, planarMPC_dvaff9, planarMPC_D1, planarMPC_dvaff8, planarMPC_grad_eq8);
planarMPC_LA_DIAGZERO_MTVM_6_6(planarMPC_D9, planarMPC_dvaff9, planarMPC_grad_eq9);
planarMPC_LA_VSUB2_96(planarMPC_rd, planarMPC_grad_eq, planarMPC_rd);
planarMPC_LA_DIAG_FORWARDBACKWARDSUB_10(planarMPC_Phi0, planarMPC_rd0, planarMPC_dzaff0);
planarMPC_LA_DIAG_FORWARDBACKWARDSUB_10(planarMPC_Phi1, planarMPC_rd1, planarMPC_dzaff1);
planarMPC_LA_DIAG_FORWARDBACKWARDSUB_10(planarMPC_Phi2, planarMPC_rd2, planarMPC_dzaff2);
planarMPC_LA_DIAG_FORWARDBACKWARDSUB_10(planarMPC_Phi3, planarMPC_rd3, planarMPC_dzaff3);
planarMPC_LA_DIAG_FORWARDBACKWARDSUB_10(planarMPC_Phi4, planarMPC_rd4, planarMPC_dzaff4);
planarMPC_LA_DIAG_FORWARDBACKWARDSUB_10(planarMPC_Phi5, planarMPC_rd5, planarMPC_dzaff5);
planarMPC_LA_DIAG_FORWARDBACKWARDSUB_10(planarMPC_Phi6, planarMPC_rd6, planarMPC_dzaff6);
planarMPC_LA_DIAG_FORWARDBACKWARDSUB_10(planarMPC_Phi7, planarMPC_rd7, planarMPC_dzaff7);
planarMPC_LA_DIAG_FORWARDBACKWARDSUB_10(planarMPC_Phi8, planarMPC_rd8, planarMPC_dzaff8);
planarMPC_LA_DIAG_FORWARDBACKWARDSUB_6(planarMPC_Phi9, planarMPC_rd9, planarMPC_dzaff9);
planarMPC_LA_VSUB_INDEXED_10(planarMPC_dzaff0, planarMPC_lbIdx0, planarMPC_rilb0, planarMPC_dslbaff0);
planarMPC_LA_VSUB3_10(planarMPC_llbbyslb0, planarMPC_dslbaff0, planarMPC_llb0, planarMPC_dllbaff0);
planarMPC_LA_VSUB2_INDEXED_10(planarMPC_riub0, planarMPC_dzaff0, planarMPC_ubIdx0, planarMPC_dsubaff0);
planarMPC_LA_VSUB3_10(planarMPC_lubbysub0, planarMPC_dsubaff0, planarMPC_lub0, planarMPC_dlubaff0);
planarMPC_LA_VSUB_INDEXED_10(planarMPC_dzaff1, planarMPC_lbIdx1, planarMPC_rilb1, planarMPC_dslbaff1);
planarMPC_LA_VSUB3_10(planarMPC_llbbyslb1, planarMPC_dslbaff1, planarMPC_llb1, planarMPC_dllbaff1);
planarMPC_LA_VSUB2_INDEXED_10(planarMPC_riub1, planarMPC_dzaff1, planarMPC_ubIdx1, planarMPC_dsubaff1);
planarMPC_LA_VSUB3_10(planarMPC_lubbysub1, planarMPC_dsubaff1, planarMPC_lub1, planarMPC_dlubaff1);
planarMPC_LA_VSUB_INDEXED_10(planarMPC_dzaff2, planarMPC_lbIdx2, planarMPC_rilb2, planarMPC_dslbaff2);
planarMPC_LA_VSUB3_10(planarMPC_llbbyslb2, planarMPC_dslbaff2, planarMPC_llb2, planarMPC_dllbaff2);
planarMPC_LA_VSUB2_INDEXED_10(planarMPC_riub2, planarMPC_dzaff2, planarMPC_ubIdx2, planarMPC_dsubaff2);
planarMPC_LA_VSUB3_10(planarMPC_lubbysub2, planarMPC_dsubaff2, planarMPC_lub2, planarMPC_dlubaff2);
planarMPC_LA_VSUB_INDEXED_10(planarMPC_dzaff3, planarMPC_lbIdx3, planarMPC_rilb3, planarMPC_dslbaff3);
planarMPC_LA_VSUB3_10(planarMPC_llbbyslb3, planarMPC_dslbaff3, planarMPC_llb3, planarMPC_dllbaff3);
planarMPC_LA_VSUB2_INDEXED_10(planarMPC_riub3, planarMPC_dzaff3, planarMPC_ubIdx3, planarMPC_dsubaff3);
planarMPC_LA_VSUB3_10(planarMPC_lubbysub3, planarMPC_dsubaff3, planarMPC_lub3, planarMPC_dlubaff3);
planarMPC_LA_VSUB_INDEXED_10(planarMPC_dzaff4, planarMPC_lbIdx4, planarMPC_rilb4, planarMPC_dslbaff4);
planarMPC_LA_VSUB3_10(planarMPC_llbbyslb4, planarMPC_dslbaff4, planarMPC_llb4, planarMPC_dllbaff4);
planarMPC_LA_VSUB2_INDEXED_10(planarMPC_riub4, planarMPC_dzaff4, planarMPC_ubIdx4, planarMPC_dsubaff4);
planarMPC_LA_VSUB3_10(planarMPC_lubbysub4, planarMPC_dsubaff4, planarMPC_lub4, planarMPC_dlubaff4);
planarMPC_LA_VSUB_INDEXED_10(planarMPC_dzaff5, planarMPC_lbIdx5, planarMPC_rilb5, planarMPC_dslbaff5);
planarMPC_LA_VSUB3_10(planarMPC_llbbyslb5, planarMPC_dslbaff5, planarMPC_llb5, planarMPC_dllbaff5);
planarMPC_LA_VSUB2_INDEXED_10(planarMPC_riub5, planarMPC_dzaff5, planarMPC_ubIdx5, planarMPC_dsubaff5);
planarMPC_LA_VSUB3_10(planarMPC_lubbysub5, planarMPC_dsubaff5, planarMPC_lub5, planarMPC_dlubaff5);
planarMPC_LA_VSUB_INDEXED_10(planarMPC_dzaff6, planarMPC_lbIdx6, planarMPC_rilb6, planarMPC_dslbaff6);
planarMPC_LA_VSUB3_10(planarMPC_llbbyslb6, planarMPC_dslbaff6, planarMPC_llb6, planarMPC_dllbaff6);
planarMPC_LA_VSUB2_INDEXED_10(planarMPC_riub6, planarMPC_dzaff6, planarMPC_ubIdx6, planarMPC_dsubaff6);
planarMPC_LA_VSUB3_10(planarMPC_lubbysub6, planarMPC_dsubaff6, planarMPC_lub6, planarMPC_dlubaff6);
planarMPC_LA_VSUB_INDEXED_10(planarMPC_dzaff7, planarMPC_lbIdx7, planarMPC_rilb7, planarMPC_dslbaff7);
planarMPC_LA_VSUB3_10(planarMPC_llbbyslb7, planarMPC_dslbaff7, planarMPC_llb7, planarMPC_dllbaff7);
planarMPC_LA_VSUB2_INDEXED_10(planarMPC_riub7, planarMPC_dzaff7, planarMPC_ubIdx7, planarMPC_dsubaff7);
planarMPC_LA_VSUB3_10(planarMPC_lubbysub7, planarMPC_dsubaff7, planarMPC_lub7, planarMPC_dlubaff7);
planarMPC_LA_VSUB_INDEXED_10(planarMPC_dzaff8, planarMPC_lbIdx8, planarMPC_rilb8, planarMPC_dslbaff8);
planarMPC_LA_VSUB3_10(planarMPC_llbbyslb8, planarMPC_dslbaff8, planarMPC_llb8, planarMPC_dllbaff8);
planarMPC_LA_VSUB2_INDEXED_10(planarMPC_riub8, planarMPC_dzaff8, planarMPC_ubIdx8, planarMPC_dsubaff8);
planarMPC_LA_VSUB3_10(planarMPC_lubbysub8, planarMPC_dsubaff8, planarMPC_lub8, planarMPC_dlubaff8);
planarMPC_LA_VSUB_INDEXED_6(planarMPC_dzaff9, planarMPC_lbIdx9, planarMPC_rilb9, planarMPC_dslbaff9);
planarMPC_LA_VSUB3_6(planarMPC_llbbyslb9, planarMPC_dslbaff9, planarMPC_llb9, planarMPC_dllbaff9);
planarMPC_LA_VSUB2_INDEXED_6(planarMPC_riub9, planarMPC_dzaff9, planarMPC_ubIdx9, planarMPC_dsubaff9);
planarMPC_LA_VSUB3_6(planarMPC_lubbysub9, planarMPC_dsubaff9, planarMPC_lub9, planarMPC_dlubaff9);
info->lsit_aff = planarMPC_LINESEARCH_BACKTRACKING_AFFINE(planarMPC_l, planarMPC_s, planarMPC_dl_aff, planarMPC_ds_aff, &info->step_aff, &info->mu_aff);
if( info->lsit_aff == planarMPC_NOPROGRESS ){
exitcode = planarMPC_NOPROGRESS; break;
}
sigma_3rdroot = info->mu_aff / info->mu;
info->sigma = sigma_3rdroot*sigma_3rdroot*sigma_3rdroot;
musigma = info->mu * info->sigma;
planarMPC_LA_VSUB5_192(planarMPC_ds_aff, planarMPC_dl_aff, musigma, planarMPC_ccrhs);
planarMPC_LA_VSUB6_INDEXED_10_10_10(planarMPC_ccrhsub0, planarMPC_sub0, planarMPC_ubIdx0, planarMPC_ccrhsl0, planarMPC_slb0, planarMPC_lbIdx0, planarMPC_rd0);
planarMPC_LA_VSUB6_INDEXED_10_10_10(planarMPC_ccrhsub1, planarMPC_sub1, planarMPC_ubIdx1, planarMPC_ccrhsl1, planarMPC_slb1, planarMPC_lbIdx1, planarMPC_rd1);
planarMPC_LA_DIAG_FORWARDSUB_10(planarMPC_Phi0, planarMPC_rd0, planarMPC_Lbyrd0);
planarMPC_LA_DIAG_FORWARDSUB_10(planarMPC_Phi1, planarMPC_rd1, planarMPC_Lbyrd1);
planarMPC_LA_DIAGZERO_MVM_6(planarMPC_W0, planarMPC_Lbyrd0, planarMPC_beta0);
planarMPC_LA_DENSE_FORWARDSUB_6(planarMPC_Ld0, planarMPC_beta0, planarMPC_yy0);
planarMPC_LA_DENSE_DIAGZERO_2MVMADD_6_10_10(planarMPC_V0, planarMPC_Lbyrd0, planarMPC_W1, planarMPC_Lbyrd1, planarMPC_beta1);
planarMPC_LA_DENSE_MVMSUB1_6_6(planarMPC_Lsd1, planarMPC_yy0, planarMPC_beta1, planarMPC_bmy1);
planarMPC_LA_DENSE_FORWARDSUB_6(planarMPC_Ld1, planarMPC_bmy1, planarMPC_yy1);
planarMPC_LA_VSUB6_INDEXED_10_10_10(planarMPC_ccrhsub2, planarMPC_sub2, planarMPC_ubIdx2, planarMPC_ccrhsl2, planarMPC_slb2, planarMPC_lbIdx2, planarMPC_rd2);
planarMPC_LA_DIAG_FORWARDSUB_10(planarMPC_Phi2, planarMPC_rd2, planarMPC_Lbyrd2);
planarMPC_LA_DENSE_DIAGZERO_2MVMADD_6_10_10(planarMPC_V1, planarMPC_Lbyrd1, planarMPC_W2, planarMPC_Lbyrd2, planarMPC_beta2);
planarMPC_LA_DENSE_MVMSUB1_6_6(planarMPC_Lsd2, planarMPC_yy1, planarMPC_beta2, planarMPC_bmy2);
planarMPC_LA_DENSE_FORWARDSUB_6(planarMPC_Ld2, planarMPC_bmy2, planarMPC_yy2);
planarMPC_LA_VSUB6_INDEXED_10_10_10(planarMPC_ccrhsub3, planarMPC_sub3, planarMPC_ubIdx3, planarMPC_ccrhsl3, planarMPC_slb3, planarMPC_lbIdx3, planarMPC_rd3);
planarMPC_LA_DIAG_FORWARDSUB_10(planarMPC_Phi3, planarMPC_rd3, planarMPC_Lbyrd3);
planarMPC_LA_DENSE_DIAGZERO_2MVMADD_6_10_10(planarMPC_V2, planarMPC_Lbyrd2, planarMPC_W3, planarMPC_Lbyrd3, planarMPC_beta3);
planarMPC_LA_DENSE_MVMSUB1_6_6(planarMPC_Lsd3, planarMPC_yy2, planarMPC_beta3, planarMPC_bmy3);
planarMPC_LA_DENSE_FORWARDSUB_6(planarMPC_Ld3, planarMPC_bmy3, planarMPC_yy3);
planarMPC_LA_VSUB6_INDEXED_10_10_10(planarMPC_ccrhsub4, planarMPC_sub4, planarMPC_ubIdx4, planarMPC_ccrhsl4, planarMPC_slb4, planarMPC_lbIdx4, planarMPC_rd4);
planarMPC_LA_DIAG_FORWARDSUB_10(planarMPC_Phi4, planarMPC_rd4, planarMPC_Lbyrd4);
planarMPC_LA_DENSE_DIAGZERO_2MVMADD_6_10_10(planarMPC_V3, planarMPC_Lbyrd3, planarMPC_W4, planarMPC_Lbyrd4, planarMPC_beta4);
planarMPC_LA_DENSE_MVMSUB1_6_6(planarMPC_Lsd4, planarMPC_yy3, planarMPC_beta4, planarMPC_bmy4);
planarMPC_LA_DENSE_FORWARDSUB_6(planarMPC_Ld4, planarMPC_bmy4, planarMPC_yy4);
planarMPC_LA_VSUB6_INDEXED_10_10_10(planarMPC_ccrhsub5, planarMPC_sub5, planarMPC_ubIdx5, planarMPC_ccrhsl5, planarMPC_slb5, planarMPC_lbIdx5, planarMPC_rd5);
planarMPC_LA_DIAG_FORWARDSUB_10(planarMPC_Phi5, planarMPC_rd5, planarMPC_Lbyrd5);
planarMPC_LA_DENSE_DIAGZERO_2MVMADD_6_10_10(planarMPC_V4, planarMPC_Lbyrd4, planarMPC_W5, planarMPC_Lbyrd5, planarMPC_beta5);
planarMPC_LA_DENSE_MVMSUB1_6_6(planarMPC_Lsd5, planarMPC_yy4, planarMPC_beta5, planarMPC_bmy5);
planarMPC_LA_DENSE_FORWARDSUB_6(planarMPC_Ld5, planarMPC_bmy5, planarMPC_yy5);
planarMPC_LA_VSUB6_INDEXED_10_10_10(planarMPC_ccrhsub6, planarMPC_sub6, planarMPC_ubIdx6, planarMPC_ccrhsl6, planarMPC_slb6, planarMPC_lbIdx6, planarMPC_rd6);
planarMPC_LA_DIAG_FORWARDSUB_10(planarMPC_Phi6, planarMPC_rd6, planarMPC_Lbyrd6);
planarMPC_LA_DENSE_DIAGZERO_2MVMADD_6_10_10(planarMPC_V5, planarMPC_Lbyrd5, planarMPC_W6, planarMPC_Lbyrd6, planarMPC_beta6);
planarMPC_LA_DENSE_MVMSUB1_6_6(planarMPC_Lsd6, planarMPC_yy5, planarMPC_beta6, planarMPC_bmy6);
planarMPC_LA_DENSE_FORWARDSUB_6(planarMPC_Ld6, planarMPC_bmy6, planarMPC_yy6);
planarMPC_LA_VSUB6_INDEXED_10_10_10(planarMPC_ccrhsub7, planarMPC_sub7, planarMPC_ubIdx7, planarMPC_ccrhsl7, planarMPC_slb7, planarMPC_lbIdx7, planarMPC_rd7);
planarMPC_LA_DIAG_FORWARDSUB_10(planarMPC_Phi7, planarMPC_rd7, planarMPC_Lbyrd7);
planarMPC_LA_DENSE_DIAGZERO_2MVMADD_6_10_10(planarMPC_V6, planarMPC_Lbyrd6, planarMPC_W7, planarMPC_Lbyrd7, planarMPC_beta7);
planarMPC_LA_DENSE_MVMSUB1_6_6(planarMPC_Lsd7, planarMPC_yy6, planarMPC_beta7, planarMPC_bmy7);
planarMPC_LA_DENSE_FORWARDSUB_6(planarMPC_Ld7, planarMPC_bmy7, planarMPC_yy7);
planarMPC_LA_VSUB6_INDEXED_10_10_10(planarMPC_ccrhsub8, planarMPC_sub8, planarMPC_ubIdx8, planarMPC_ccrhsl8, planarMPC_slb8, planarMPC_lbIdx8, planarMPC_rd8);
planarMPC_LA_DIAG_FORWARDSUB_10(planarMPC_Phi8, planarMPC_rd8, planarMPC_Lbyrd8);
planarMPC_LA_DENSE_DIAGZERO_2MVMADD_6_10_10(planarMPC_V7, planarMPC_Lbyrd7, planarMPC_W8, planarMPC_Lbyrd8, planarMPC_beta8);
planarMPC_LA_DENSE_MVMSUB1_6_6(planarMPC_Lsd8, planarMPC_yy7, planarMPC_beta8, planarMPC_bmy8);
planarMPC_LA_DENSE_FORWARDSUB_6(planarMPC_Ld8, planarMPC_bmy8, planarMPC_yy8);
planarMPC_LA_VSUB6_INDEXED_6_6_6(planarMPC_ccrhsub9, planarMPC_sub9, planarMPC_ubIdx9, planarMPC_ccrhsl9, planarMPC_slb9, planarMPC_lbIdx9, planarMPC_rd9);
planarMPC_LA_DIAG_FORWARDSUB_6(planarMPC_Phi9, planarMPC_rd9, planarMPC_Lbyrd9);
planarMPC_LA_DENSE_DIAGZERO_2MVMADD_6_10_6(planarMPC_V8, planarMPC_Lbyrd8, planarMPC_W9, planarMPC_Lbyrd9, planarMPC_beta9);
planarMPC_LA_DENSE_MVMSUB1_6_6(planarMPC_Lsd9, planarMPC_yy8, planarMPC_beta9, planarMPC_bmy9);
planarMPC_LA_DENSE_FORWARDSUB_6(planarMPC_Ld9, planarMPC_bmy9, planarMPC_yy9);
planarMPC_LA_DENSE_BACKWARDSUB_6(planarMPC_Ld9, planarMPC_yy9, planarMPC_dvcc9);
planarMPC_LA_DENSE_MTVMSUB_6_6(planarMPC_Lsd9, planarMPC_dvcc9, planarMPC_yy8, planarMPC_bmy8);
planarMPC_LA_DENSE_BACKWARDSUB_6(planarMPC_Ld8, planarMPC_bmy8, planarMPC_dvcc8);
planarMPC_LA_DENSE_MTVMSUB_6_6(planarMPC_Lsd8, planarMPC_dvcc8, planarMPC_yy7, planarMPC_bmy7);
planarMPC_LA_DENSE_BACKWARDSUB_6(planarMPC_Ld7, planarMPC_bmy7, planarMPC_dvcc7);
planarMPC_LA_DENSE_MTVMSUB_6_6(planarMPC_Lsd7, planarMPC_dvcc7, planarMPC_yy6, planarMPC_bmy6);
planarMPC_LA_DENSE_BACKWARDSUB_6(planarMPC_Ld6, planarMPC_bmy6, planarMPC_dvcc6);
planarMPC_LA_DENSE_MTVMSUB_6_6(planarMPC_Lsd6, planarMPC_dvcc6, planarMPC_yy5, planarMPC_bmy5);
planarMPC_LA_DENSE_BACKWARDSUB_6(planarMPC_Ld5, planarMPC_bmy5, planarMPC_dvcc5);
planarMPC_LA_DENSE_MTVMSUB_6_6(planarMPC_Lsd5, planarMPC_dvcc5, planarMPC_yy4, planarMPC_bmy4);
planarMPC_LA_DENSE_BACKWARDSUB_6(planarMPC_Ld4, planarMPC_bmy4, planarMPC_dvcc4);
planarMPC_LA_DENSE_MTVMSUB_6_6(planarMPC_Lsd4, planarMPC_dvcc4, planarMPC_yy3, planarMPC_bmy3);
planarMPC_LA_DENSE_BACKWARDSUB_6(planarMPC_Ld3, planarMPC_bmy3, planarMPC_dvcc3);
planarMPC_LA_DENSE_MTVMSUB_6_6(planarMPC_Lsd3, planarMPC_dvcc3, planarMPC_yy2, planarMPC_bmy2);
planarMPC_LA_DENSE_BACKWARDSUB_6(planarMPC_Ld2, planarMPC_bmy2, planarMPC_dvcc2);
planarMPC_LA_DENSE_MTVMSUB_6_6(planarMPC_Lsd2, planarMPC_dvcc2, planarMPC_yy1, planarMPC_bmy1);
planarMPC_LA_DENSE_BACKWARDSUB_6(planarMPC_Ld1, planarMPC_bmy1, planarMPC_dvcc1);
planarMPC_LA_DENSE_MTVMSUB_6_6(planarMPC_Lsd1, planarMPC_dvcc1, planarMPC_yy0, planarMPC_bmy0);
planarMPC_LA_DENSE_BACKWARDSUB_6(planarMPC_Ld0, planarMPC_bmy0, planarMPC_dvcc0);
planarMPC_LA_DENSE_DIAGZERO_MTVM2_6_10_6(planarMPC_C0, planarMPC_dvcc1, planarMPC_D0, planarMPC_dvcc0, planarMPC_grad_eq0);
planarMPC_LA_DENSE_DIAGZERO_MTVM2_6_10_6(planarMPC_C0, planarMPC_dvcc2, planarMPC_D1, planarMPC_dvcc1, planarMPC_grad_eq1);
planarMPC_LA_DENSE_DIAGZERO_MTVM2_6_10_6(planarMPC_C0, planarMPC_dvcc3, planarMPC_D1, planarMPC_dvcc2, planarMPC_grad_eq2);
planarMPC_LA_DENSE_DIAGZERO_MTVM2_6_10_6(planarMPC_C0, planarMPC_dvcc4, planarMPC_D1, planarMPC_dvcc3, planarMPC_grad_eq3);
planarMPC_LA_DENSE_DIAGZERO_MTVM2_6_10_6(planarMPC_C0, planarMPC_dvcc5, planarMPC_D1, planarMPC_dvcc4, planarMPC_grad_eq4);
planarMPC_LA_DENSE_DIAGZERO_MTVM2_6_10_6(planarMPC_C0, planarMPC_dvcc6, planarMPC_D1, planarMPC_dvcc5, planarMPC_grad_eq5);
planarMPC_LA_DENSE_DIAGZERO_MTVM2_6_10_6(planarMPC_C0, planarMPC_dvcc7, planarMPC_D1, planarMPC_dvcc6, planarMPC_grad_eq6);
planarMPC_LA_DENSE_DIAGZERO_MTVM2_6_10_6(planarMPC_C0, planarMPC_dvcc8, planarMPC_D1, planarMPC_dvcc7, planarMPC_grad_eq7);
planarMPC_LA_DENSE_DIAGZERO_MTVM2_6_10_6(planarMPC_C0, planarMPC_dvcc9, planarMPC_D1, planarMPC_dvcc8, planarMPC_grad_eq8);
planarMPC_LA_DIAGZERO_MTVM_6_6(planarMPC_D9, planarMPC_dvcc9, planarMPC_grad_eq9);
planarMPC_LA_VSUB_96(planarMPC_rd, planarMPC_grad_eq, planarMPC_rd);
planarMPC_LA_DIAG_FORWARDBACKWARDSUB_10(planarMPC_Phi0, planarMPC_rd0, planarMPC_dzcc0);
planarMPC_LA_DIAG_FORWARDBACKWARDSUB_10(planarMPC_Phi1, planarMPC_rd1, planarMPC_dzcc1);
planarMPC_LA_DIAG_FORWARDBACKWARDSUB_10(planarMPC_Phi2, planarMPC_rd2, planarMPC_dzcc2);
planarMPC_LA_DIAG_FORWARDBACKWARDSUB_10(planarMPC_Phi3, planarMPC_rd3, planarMPC_dzcc3);
planarMPC_LA_DIAG_FORWARDBACKWARDSUB_10(planarMPC_Phi4, planarMPC_rd4, planarMPC_dzcc4);
planarMPC_LA_DIAG_FORWARDBACKWARDSUB_10(planarMPC_Phi5, planarMPC_rd5, planarMPC_dzcc5);
planarMPC_LA_DIAG_FORWARDBACKWARDSUB_10(planarMPC_Phi6, planarMPC_rd6, planarMPC_dzcc6);
planarMPC_LA_DIAG_FORWARDBACKWARDSUB_10(planarMPC_Phi7, planarMPC_rd7, planarMPC_dzcc7);
planarMPC_LA_DIAG_FORWARDBACKWARDSUB_10(planarMPC_Phi8, planarMPC_rd8, planarMPC_dzcc8);
planarMPC_LA_DIAG_FORWARDBACKWARDSUB_6(planarMPC_Phi9, planarMPC_rd9, planarMPC_dzcc9);
planarMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_10(planarMPC_ccrhsl0, planarMPC_slb0, planarMPC_llbbyslb0, planarMPC_dzcc0, planarMPC_lbIdx0, planarMPC_dllbcc0);
planarMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_10(planarMPC_ccrhsub0, planarMPC_sub0, planarMPC_lubbysub0, planarMPC_dzcc0, planarMPC_ubIdx0, planarMPC_dlubcc0);
planarMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_10(planarMPC_ccrhsl1, planarMPC_slb1, planarMPC_llbbyslb1, planarMPC_dzcc1, planarMPC_lbIdx1, planarMPC_dllbcc1);
planarMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_10(planarMPC_ccrhsub1, planarMPC_sub1, planarMPC_lubbysub1, planarMPC_dzcc1, planarMPC_ubIdx1, planarMPC_dlubcc1);
planarMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_10(planarMPC_ccrhsl2, planarMPC_slb2, planarMPC_llbbyslb2, planarMPC_dzcc2, planarMPC_lbIdx2, planarMPC_dllbcc2);
planarMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_10(planarMPC_ccrhsub2, planarMPC_sub2, planarMPC_lubbysub2, planarMPC_dzcc2, planarMPC_ubIdx2, planarMPC_dlubcc2);
planarMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_10(planarMPC_ccrhsl3, planarMPC_slb3, planarMPC_llbbyslb3, planarMPC_dzcc3, planarMPC_lbIdx3, planarMPC_dllbcc3);
planarMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_10(planarMPC_ccrhsub3, planarMPC_sub3, planarMPC_lubbysub3, planarMPC_dzcc3, planarMPC_ubIdx3, planarMPC_dlubcc3);
planarMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_10(planarMPC_ccrhsl4, planarMPC_slb4, planarMPC_llbbyslb4, planarMPC_dzcc4, planarMPC_lbIdx4, planarMPC_dllbcc4);
planarMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_10(planarMPC_ccrhsub4, planarMPC_sub4, planarMPC_lubbysub4, planarMPC_dzcc4, planarMPC_ubIdx4, planarMPC_dlubcc4);
planarMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_10(planarMPC_ccrhsl5, planarMPC_slb5, planarMPC_llbbyslb5, planarMPC_dzcc5, planarMPC_lbIdx5, planarMPC_dllbcc5);
planarMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_10(planarMPC_ccrhsub5, planarMPC_sub5, planarMPC_lubbysub5, planarMPC_dzcc5, planarMPC_ubIdx5, planarMPC_dlubcc5);
planarMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_10(planarMPC_ccrhsl6, planarMPC_slb6, planarMPC_llbbyslb6, planarMPC_dzcc6, planarMPC_lbIdx6, planarMPC_dllbcc6);
planarMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_10(planarMPC_ccrhsub6, planarMPC_sub6, planarMPC_lubbysub6, planarMPC_dzcc6, planarMPC_ubIdx6, planarMPC_dlubcc6);
planarMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_10(planarMPC_ccrhsl7, planarMPC_slb7, planarMPC_llbbyslb7, planarMPC_dzcc7, planarMPC_lbIdx7, planarMPC_dllbcc7);
planarMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_10(planarMPC_ccrhsub7, planarMPC_sub7, planarMPC_lubbysub7, planarMPC_dzcc7, planarMPC_ubIdx7, planarMPC_dlubcc7);
planarMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_10(planarMPC_ccrhsl8, planarMPC_slb8, planarMPC_llbbyslb8, planarMPC_dzcc8, planarMPC_lbIdx8, planarMPC_dllbcc8);
planarMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_10(planarMPC_ccrhsub8, planarMPC_sub8, planarMPC_lubbysub8, planarMPC_dzcc8, planarMPC_ubIdx8, planarMPC_dlubcc8);
planarMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_6(planarMPC_ccrhsl9, planarMPC_slb9, planarMPC_llbbyslb9, planarMPC_dzcc9, planarMPC_lbIdx9, planarMPC_dllbcc9);
planarMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_6(planarMPC_ccrhsub9, planarMPC_sub9, planarMPC_lubbysub9, planarMPC_dzcc9, planarMPC_ubIdx9, planarMPC_dlubcc9);
planarMPC_LA_VSUB7_192(planarMPC_l, planarMPC_ccrhs, planarMPC_s, planarMPC_dl_cc, planarMPC_ds_cc);
planarMPC_LA_VADD_96(planarMPC_dz_cc, planarMPC_dz_aff);
planarMPC_LA_VADD_60(planarMPC_dv_cc, planarMPC_dv_aff);
planarMPC_LA_VADD_192(planarMPC_dl_cc, planarMPC_dl_aff);
planarMPC_LA_VADD_192(planarMPC_ds_cc, planarMPC_ds_aff);
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
output->z5[0] = planarMPC_z4[0];
output->z5[1] = planarMPC_z4[1];
output->z5[2] = planarMPC_z4[2];
output->z5[3] = planarMPC_z4[3];
output->z5[4] = planarMPC_z4[4];
output->z5[5] = planarMPC_z4[5];
output->z5[6] = planarMPC_z4[6];
output->z5[7] = planarMPC_z4[7];
output->z5[8] = planarMPC_z4[8];
output->z5[9] = planarMPC_z4[9];
output->z6[0] = planarMPC_z5[0];
output->z6[1] = planarMPC_z5[1];
output->z6[2] = planarMPC_z5[2];
output->z6[3] = planarMPC_z5[3];
output->z6[4] = planarMPC_z5[4];
output->z6[5] = planarMPC_z5[5];
output->z6[6] = planarMPC_z5[6];
output->z6[7] = planarMPC_z5[7];
output->z6[8] = planarMPC_z5[8];
output->z6[9] = planarMPC_z5[9];
output->z7[0] = planarMPC_z6[0];
output->z7[1] = planarMPC_z6[1];
output->z7[2] = planarMPC_z6[2];
output->z7[3] = planarMPC_z6[3];
output->z7[4] = planarMPC_z6[4];
output->z7[5] = planarMPC_z6[5];
output->z7[6] = planarMPC_z6[6];
output->z7[7] = planarMPC_z6[7];
output->z7[8] = planarMPC_z6[8];
output->z7[9] = planarMPC_z6[9];
output->z8[0] = planarMPC_z7[0];
output->z8[1] = planarMPC_z7[1];
output->z8[2] = planarMPC_z7[2];
output->z8[3] = planarMPC_z7[3];
output->z8[4] = planarMPC_z7[4];
output->z8[5] = planarMPC_z7[5];
output->z8[6] = planarMPC_z7[6];
output->z8[7] = planarMPC_z7[7];
output->z8[8] = planarMPC_z7[8];
output->z8[9] = planarMPC_z7[9];
output->z9[0] = planarMPC_z8[0];
output->z9[1] = planarMPC_z8[1];
output->z9[2] = planarMPC_z8[2];
output->z9[3] = planarMPC_z8[3];
output->z9[4] = planarMPC_z8[4];
output->z9[5] = planarMPC_z8[5];
output->z9[6] = planarMPC_z8[6];
output->z9[7] = planarMPC_z8[7];
output->z9[8] = planarMPC_z8[8];
output->z9[9] = planarMPC_z8[9];
output->z10[0] = planarMPC_z9[0];
output->z10[1] = planarMPC_z9[1];
output->z10[2] = planarMPC_z9[2];
output->z10[3] = planarMPC_z9[3];
output->z10[4] = planarMPC_z9[4];
output->z10[5] = planarMPC_z9[5];

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
