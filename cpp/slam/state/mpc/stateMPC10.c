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
 * Initializes a vector of length 270 with a value.
 */
void stateMPC_LA_INITIALIZEVECTOR_270(stateMPC_FLOAT* vec, stateMPC_FLOAT value)
{
	int i;
	for( i=0; i<270; i++ )
	{
		vec[i] = value;
	}
}


/*
 * Initializes a vector of length 90 with a value.
 */
void stateMPC_LA_INITIALIZEVECTOR_90(stateMPC_FLOAT* vec, stateMPC_FLOAT value)
{
	int i;
	for( i=0; i<90; i++ )
	{
		vec[i] = value;
	}
}


/*
 * Initializes a vector of length 378 with a value.
 */
void stateMPC_LA_INITIALIZEVECTOR_378(stateMPC_FLOAT* vec, stateMPC_FLOAT value)
{
	int i;
	for( i=0; i<378; i++ )
	{
		vec[i] = value;
	}
}


/* 
 * Calculates a dot product and adds it to a variable: z += x'*y; 
 * This function is for vectors of length 378.
 */
void stateMPC_LA_DOTACC_378(stateMPC_FLOAT *x, stateMPC_FLOAT *y, stateMPC_FLOAT *z)
{
	int i;
	for( i=0; i<378; i++ ){
		*z += x[i]*y[i];
	}
}


/*
 * Calculates the gradient and the value for a quadratic function 0.5*z'*H*z + f'*z
 *
 * INPUTS:     H  - Symmetric Hessian, diag matrix of size [29 x 29]
 *             f  - column vector of size 29
 *             z  - column vector of size 29
 *
 * OUTPUTS: grad  - gradient at z (= H*z + f), column vector of size 29
 *          value <-- value + 0.5*z'*H*z + f'*z (value will be modified)
 */
void stateMPC_LA_DIAG_QUADFCN_29(stateMPC_FLOAT* H, stateMPC_FLOAT* f, stateMPC_FLOAT* z, stateMPC_FLOAT* grad, stateMPC_FLOAT* value)
{
	int i;
	stateMPC_FLOAT hz;	
	for( i=0; i<29; i++){
		hz = H[i]*z[i];
		grad[i] = hz + f[i];
		*value += 0.5*hz*z[i] + f[i]*z[i];
	}
}


/*
 * Calculates the gradient and the value for a quadratic function 0.5*z'*H*z + f'*z
 *
 * INPUTS:     H  - Symmetric Hessian, diag matrix of size [9 x 9]
 *             f  - column vector of size 9
 *             z  - column vector of size 9
 *
 * OUTPUTS: grad  - gradient at z (= H*z + f), column vector of size 9
 *          value <-- value + 0.5*z'*H*z + f'*z (value will be modified)
 */
void stateMPC_LA_DIAG_QUADFCN_9(stateMPC_FLOAT* H, stateMPC_FLOAT* f, stateMPC_FLOAT* z, stateMPC_FLOAT* grad, stateMPC_FLOAT* value)
{
	int i;
	stateMPC_FLOAT hz;	
	for( i=0; i<9; i++){
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
void stateMPC_LA_DIAGZERO_MVMSUB6_9(stateMPC_FLOAT *B, stateMPC_FLOAT *u, stateMPC_FLOAT *b, stateMPC_FLOAT *l, stateMPC_FLOAT *r, stateMPC_FLOAT *z, stateMPC_FLOAT *y)
{
	int i;
	stateMPC_FLOAT Bu[9];
	stateMPC_FLOAT norm = *y;
	stateMPC_FLOAT lr = 0;

	/* do A*x + B*u first */
	for( i=0; i<9; i++ ){
		Bu[i] = B[i]*u[i];
	}	

	for( i=0; i<9; i++ ){
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
void stateMPC_LA_DENSE_DIAGZERO_MVMSUB3_9_29_29(stateMPC_FLOAT *A, stateMPC_FLOAT *x, stateMPC_FLOAT *B, stateMPC_FLOAT *u, stateMPC_FLOAT *b, stateMPC_FLOAT *l, stateMPC_FLOAT *r, stateMPC_FLOAT *z, stateMPC_FLOAT *y)
{
	int i;
	int j;
	int k = 0;
	stateMPC_FLOAT AxBu[9];
	stateMPC_FLOAT norm = *y;
	stateMPC_FLOAT lr = 0;

	/* do A*x + B*u first */
	for( i=0; i<9; i++ ){
		AxBu[i] = A[k++]*x[0] + B[i]*u[i];
	}	

	for( j=1; j<29; j++ ){		
		for( i=0; i<9; i++ ){
			AxBu[i] += A[k++]*x[j];
		}
	}

	for( i=0; i<9; i++ ){
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
void stateMPC_LA_DENSE_DIAGZERO_MVMSUB3_9_29_9(stateMPC_FLOAT *A, stateMPC_FLOAT *x, stateMPC_FLOAT *B, stateMPC_FLOAT *u, stateMPC_FLOAT *b, stateMPC_FLOAT *l, stateMPC_FLOAT *r, stateMPC_FLOAT *z, stateMPC_FLOAT *y)
{
	int i;
	int j;
	int k = 0;
	stateMPC_FLOAT AxBu[9];
	stateMPC_FLOAT norm = *y;
	stateMPC_FLOAT lr = 0;

	/* do A*x + B*u first */
	for( i=0; i<9; i++ ){
		AxBu[i] = A[k++]*x[0] + B[i]*u[i];
	}	

	for( j=1; j<29; j++ ){		
		for( i=0; i<9; i++ ){
			AxBu[i] += A[k++]*x[j];
		}
	}

	for( i=0; i<9; i++ ){
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
 * where A is of size [9 x 29] and stored in column major format.
 * and B is of size [9 x 29] and stored in diagzero format
 * Note the transposes of A and B!
 */
void stateMPC_LA_DENSE_DIAGZERO_MTVM2_9_29_9(stateMPC_FLOAT *A, stateMPC_FLOAT *x, stateMPC_FLOAT *B, stateMPC_FLOAT *y, stateMPC_FLOAT *z)
{
	int i;
	int j;
	int k = 0;
	for( i=0; i<9; i++ ){
		z[i] = 0;
		for( j=0; j<9; j++ ){
			z[i] += A[k++]*x[j];
		}
		z[i] += B[i]*y[i];
	}
	for( i=9 ;i<29; i++ ){
		z[i] = 0;
		for( j=0; j<9; j++ ){
			z[i] += A[k++]*x[j];
		}
	}
}


/*
 * Matrix vector multiplication y = M'*x where M is of size [9 x 9]
 * and stored in diagzero format. Note the transpose of M!
 */
void stateMPC_LA_DIAGZERO_MTVM_9_9(stateMPC_FLOAT *M, stateMPC_FLOAT *x, stateMPC_FLOAT *y)
{
	int i;
	for( i=0; i<9; i++ ){
		y[i] = M[i]*x[i];
	}
}


/*
 * Vector subtraction and addition.
 *	 Input: five vectors t, tidx, u, v, w and two scalars z and r
 *	 Output: y = t(tidx) - u + w
 *           z = z - v'*x;
 *           r = max([norm(y,inf), z]);
 * for vectors of length 29. Output z is of course scalar.
 */
void stateMPC_LA_VSUBADD3_29(stateMPC_FLOAT* t, stateMPC_FLOAT* u, int* uidx, stateMPC_FLOAT* v, stateMPC_FLOAT* w, stateMPC_FLOAT* y, stateMPC_FLOAT* z, stateMPC_FLOAT* r)
{
	int i;
	stateMPC_FLOAT norm = *r;
	stateMPC_FLOAT vx = 0;
	stateMPC_FLOAT x;
	for( i=0; i<29; i++){
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
 * for vectors of length 11. Output z is of course scalar.
 */
void stateMPC_LA_VSUBADD2_11(stateMPC_FLOAT* t, int* tidx, stateMPC_FLOAT* u, stateMPC_FLOAT* v, stateMPC_FLOAT* w, stateMPC_FLOAT* y, stateMPC_FLOAT* z, stateMPC_FLOAT* r)
{
	int i;
	stateMPC_FLOAT norm = *r;
	stateMPC_FLOAT vx = 0;
	stateMPC_FLOAT x;
	for( i=0; i<11; i++){
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
 * for vectors of length 9. Output z is of course scalar.
 */
void stateMPC_LA_VSUBADD3_9(stateMPC_FLOAT* t, stateMPC_FLOAT* u, int* uidx, stateMPC_FLOAT* v, stateMPC_FLOAT* w, stateMPC_FLOAT* y, stateMPC_FLOAT* z, stateMPC_FLOAT* r)
{
	int i;
	stateMPC_FLOAT norm = *r;
	stateMPC_FLOAT vx = 0;
	stateMPC_FLOAT x;
	for( i=0; i<9; i++){
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
 * for vectors of length 9. Output z is of course scalar.
 */
void stateMPC_LA_VSUBADD2_9(stateMPC_FLOAT* t, int* tidx, stateMPC_FLOAT* u, stateMPC_FLOAT* v, stateMPC_FLOAT* w, stateMPC_FLOAT* y, stateMPC_FLOAT* z, stateMPC_FLOAT* r)
{
	int i;
	stateMPC_FLOAT norm = *r;
	stateMPC_FLOAT vx = 0;
	stateMPC_FLOAT x;
	for( i=0; i<9; i++){
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
 * Special function for box constraints of length 29
 * Returns also L/S, a value that is often used elsewhere.
 */
void stateMPC_LA_INEQ_B_GRAD_29_29_11(stateMPC_FLOAT *lu, stateMPC_FLOAT *su, stateMPC_FLOAT *ru, stateMPC_FLOAT *ll, stateMPC_FLOAT *sl, stateMPC_FLOAT *rl, int* lbIdx, int* ubIdx, stateMPC_FLOAT *grad, stateMPC_FLOAT *lubysu, stateMPC_FLOAT *llbysl)
{
	int i;
	for( i=0; i<29; i++ ){
		grad[i] = 0;
	}
	for( i=0; i<29; i++ ){		
		llbysl[i] = ll[i] / sl[i];
		grad[lbIdx[i]] -= llbysl[i]*rl[i];
	}
	for( i=0; i<11; i++ ){
		lubysu[i] = lu[i] / su[i];
		grad[ubIdx[i]] += lubysu[i]*ru[i];
	}
}


/*
 * Computes inequality constraints gradient-
 * Special function for box constraints of length 9
 * Returns also L/S, a value that is often used elsewhere.
 */
void stateMPC_LA_INEQ_B_GRAD_9_9_9(stateMPC_FLOAT *lu, stateMPC_FLOAT *su, stateMPC_FLOAT *ru, stateMPC_FLOAT *ll, stateMPC_FLOAT *sl, stateMPC_FLOAT *rl, int* lbIdx, int* ubIdx, stateMPC_FLOAT *grad, stateMPC_FLOAT *lubysu, stateMPC_FLOAT *llbysl)
{
	int i;
	for( i=0; i<9; i++ ){
		grad[i] = 0;
	}
	for( i=0; i<9; i++ ){		
		llbysl[i] = ll[i] / sl[i];
		grad[lbIdx[i]] -= llbysl[i]*rl[i];
	}
	for( i=0; i<9; i++ ){
		lubysu[i] = lu[i] / su[i];
		grad[ubIdx[i]] += lubysu[i]*ru[i];
	}
}


/*
 * Addition of three vectors  z = u + w + v
 * of length 270.
 */
void stateMPC_LA_VVADD3_270(stateMPC_FLOAT *u, stateMPC_FLOAT *v, stateMPC_FLOAT *w, stateMPC_FLOAT *z)
{
	int i;
	for( i=0; i<270; i++ ){
		z[i] = u[i] + v[i] + w[i];
	}
}


/*
 * Special function to compute the diagonal cholesky factorization of the 
 * positive definite augmented Hessian for block size 29.
 *
 * Inputs: - H = diagonal cost Hessian in diagonal storage format
 *         - llbysl = L / S of lower bounds
 *         - lubysu = L / S of upper bounds
 *
 * Output: Phi = sqrt(H + diag(llbysl) + diag(lubysu))
 * where Phi is stored in diagonal storage format
 */
void stateMPC_LA_DIAG_CHOL_LBUB_29_29_11(stateMPC_FLOAT *H, stateMPC_FLOAT *llbysl, int* lbIdx, stateMPC_FLOAT *lubysu, int* ubIdx, stateMPC_FLOAT *Phi)


{
	int i;
	
	/* copy  H into PHI */
	for( i=0; i<29; i++ ){
		Phi[i] = H[i];
	}

	/* add llbysl onto Phi where necessary */
	for( i=0; i<29; i++ ){
		Phi[lbIdx[i]] += llbysl[i];
	}

	/* add lubysu onto Phi where necessary */
	for( i=0; i<11; i++){
		Phi[ubIdx[i]] +=  lubysu[i];
	}
	
	/* compute cholesky */
	for(i=0; i<29; i++)
	{
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
 * where A is to be computed and is of size [9 x 29],
 * B is given and of size [9 x 29], L is a diagonal
 * matrix of size 9 stored in diagonal matrix 
 * storage format. Note the transpose of L has no impact!
 *
 * Result: A in column major storage format.
 *
 */
void stateMPC_LA_DIAG_MATRIXFORWARDSUB_9_29(stateMPC_FLOAT *L, stateMPC_FLOAT *B, stateMPC_FLOAT *A)
{
    int i,j;
	 int k = 0;

	for( j=0; j<29; j++){
		for( i=0; i<9; i++){
			A[k] = B[k]/L[j];
			k++;
		}
	}

}


/**
 * Forward substitution for the matrix equation A*L' = B
 * where A is to be computed and is of size [9 x 29],
 * B is given and of size [9 x 29], L is a diagonal
 *  matrix of size 29 stored in diagonal 
 * storage format. Note the transpose of L!
 *
 * Result: A in diagonalzero storage format.
 *
 */
void stateMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_9_29(stateMPC_FLOAT *L, stateMPC_FLOAT *B, stateMPC_FLOAT *A)
{
	int j;
    for( j=0; j<29; j++ ){   
		A[j] = B[j]/L[j];
     }
}


/**
 * Compute C = A*B' where 
 *
 *	size(A) = [9 x 29]
 *  size(B) = [9 x 29] in diagzero format
 * 
 * A and C matrices are stored in column major format.
 * 
 * 
 */
void stateMPC_LA_DENSE_DIAGZERO_MMTM_9_29_9(stateMPC_FLOAT *A, stateMPC_FLOAT *B, stateMPC_FLOAT *C)
{
    int i, j;
	
	for( i=0; i<9; i++ ){
		for( j=0; j<9; j++){
			C[j*9+i] = B[i*9+j]*A[i];
		}
	}

}


/**
 * Forward substitution to solve L*y = b where L is a
 * diagonal matrix in vector storage format.
 * 
 * The dimensions involved are 29.
 */
void stateMPC_LA_DIAG_FORWARDSUB_29(stateMPC_FLOAT *L, stateMPC_FLOAT *b, stateMPC_FLOAT *y)
{
    int i;

    for( i=0; i<29; i++ ){
		y[i] = b[i]/L[i];
    }
}


/*
 * Special function to compute the diagonal cholesky factorization of the 
 * positive definite augmented Hessian for block size 9.
 *
 * Inputs: - H = diagonal cost Hessian in diagonal storage format
 *         - llbysl = L / S of lower bounds
 *         - lubysu = L / S of upper bounds
 *
 * Output: Phi = sqrt(H + diag(llbysl) + diag(lubysu))
 * where Phi is stored in diagonal storage format
 */
void stateMPC_LA_DIAG_CHOL_ONELOOP_LBUB_9_9_9(stateMPC_FLOAT *H, stateMPC_FLOAT *llbysl, int* lbIdx, stateMPC_FLOAT *lubysu, int* ubIdx, stateMPC_FLOAT *Phi)


{
	int i;
	
	/* compute cholesky */
	for( i=0; i<9; i++ ){
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
 * where A is to be computed and is of size [9 x 9],
 * B is given and of size [9 x 9], L is a diagonal
 *  matrix of size 9 stored in diagonal 
 * storage format. Note the transpose of L!
 *
 * Result: A in diagonalzero storage format.
 *
 */
void stateMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_9_9(stateMPC_FLOAT *L, stateMPC_FLOAT *B, stateMPC_FLOAT *A)
{
	int j;
    for( j=0; j<9; j++ ){   
		A[j] = B[j]/L[j];
     }
}


/**
 * Forward substitution to solve L*y = b where L is a
 * diagonal matrix in vector storage format.
 * 
 * The dimensions involved are 9.
 */
void stateMPC_LA_DIAG_FORWARDSUB_9(stateMPC_FLOAT *L, stateMPC_FLOAT *b, stateMPC_FLOAT *y)
{
    int i;

    for( i=0; i<9; i++ ){
		y[i] = b[i]/L[i];
    }
}


/**
 * Compute L = B*B', where L is lower triangular of size NXp1
 * and B is a diagzero matrix of size [9 x 29] in column
 * storage format.
 * 
 */
void stateMPC_LA_DIAGZERO_MMT_9(stateMPC_FLOAT *B, stateMPC_FLOAT *L)
{
    int i, ii, di;
    
    ii = 0; di = 0;
    for( i=0; i<9; i++ ){        
		L[ii+i] = B[i]*B[i];
        ii += ++di;
    }
}


/* 
 * Computes r = b - B*u
 * B is stored in diagzero format
 */
void stateMPC_LA_DIAGZERO_MVMSUB7_9(stateMPC_FLOAT *B, stateMPC_FLOAT *u, stateMPC_FLOAT *b, stateMPC_FLOAT *r)
{
	int i;

	for( i=0; i<9; i++ ){
		r[i] = b[i] - B[i]*u[i];
	}	
	
}


/**
 * Compute L = A*A' + B*B', where L is lower triangular of size NXp1
 * and A is a dense matrix of size [9 x 29] in column
 * storage format, and B is of size [9 x 29] diagonalzero
 * storage format.
 * 
 * THIS ONE HAS THE WORST ACCES PATTERN POSSIBLE. 
 * POSSIBKE FIX: PUT A AND B INTO ROW MAJOR FORMAT FIRST.
 * 
 */
void stateMPC_LA_DENSE_DIAGZERO_MMT2_9_29_29(stateMPC_FLOAT *A, stateMPC_FLOAT *B, stateMPC_FLOAT *L)
{
    int i, j, k, ii, di;
    stateMPC_FLOAT ltemp;
    
    ii = 0; di = 0;
    for( i=0; i<9; i++ ){        
        for( j=0; j<=i; j++ ){
            ltemp = 0; 
            for( k=0; k<29; k++ ){
                ltemp += A[k*9+i]*A[k*9+j];
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
void stateMPC_LA_DENSE_DIAGZERO_2MVMSUB2_9_29_29(stateMPC_FLOAT *A, stateMPC_FLOAT *x, stateMPC_FLOAT *B, stateMPC_FLOAT *u, stateMPC_FLOAT *b, stateMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;

	for( i=0; i<9; i++ ){
		r[i] = b[i] - A[k++]*x[0] - B[i]*u[i];
	}	

	for( j=1; j<29; j++ ){		
		for( i=0; i<9; i++ ){
			r[i] -= A[k++]*x[j];
		}
	}
	
}


/**
 * Compute L = A*A' + B*B', where L is lower triangular of size NXp1
 * and A is a dense matrix of size [9 x 29] in column
 * storage format, and B is of size [9 x 9] diagonalzero
 * storage format.
 * 
 * THIS ONE HAS THE WORST ACCES PATTERN POSSIBLE. 
 * POSSIBKE FIX: PUT A AND B INTO ROW MAJOR FORMAT FIRST.
 * 
 */
void stateMPC_LA_DENSE_DIAGZERO_MMT2_9_29_9(stateMPC_FLOAT *A, stateMPC_FLOAT *B, stateMPC_FLOAT *L)
{
    int i, j, k, ii, di;
    stateMPC_FLOAT ltemp;
    
    ii = 0; di = 0;
    for( i=0; i<9; i++ ){        
        for( j=0; j<=i; j++ ){
            ltemp = 0; 
            for( k=0; k<29; k++ ){
                ltemp += A[k*9+i]*A[k*9+j];
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
void stateMPC_LA_DENSE_DIAGZERO_2MVMSUB2_9_29_9(stateMPC_FLOAT *A, stateMPC_FLOAT *x, stateMPC_FLOAT *B, stateMPC_FLOAT *u, stateMPC_FLOAT *b, stateMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;

	for( i=0; i<9; i++ ){
		r[i] = b[i] - A[k++]*x[0] - B[i]*u[i];
	}	

	for( j=1; j<29; j++ ){		
		for( i=0; i<9; i++ ){
			r[i] -= A[k++]*x[j];
		}
	}
	
}


/**
 * Cholesky factorization as above, but working on a matrix in 
 * lower triangular storage format of size 9 and outputting
 * the Cholesky factor to matrix L in lower triangular format.
 */
void stateMPC_LA_DENSE_CHOL_9(stateMPC_FLOAT *A, stateMPC_FLOAT *L)
{
    int i, j, k, di, dj;
	 int ii, jj;

    stateMPC_FLOAT l;
    stateMPC_FLOAT Mii;

	/* copy A to L first and then operate on L */
	/* COULD BE OPTIMIZED */
	ii=0; di=0;
	for( i=0; i<9; i++ ){
		for( j=0; j<=i; j++ ){
			L[ii+j] = A[ii+j];
		}
		ii += ++di;
	}    
	
	/* factor L */
	ii=0; di=0;
    for( i=0; i<9; i++ ){
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
        for( j=i+1; j<9; j++ ){
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
 * The dimensions involved are 9.
 */
void stateMPC_LA_DENSE_FORWARDSUB_9(stateMPC_FLOAT *L, stateMPC_FLOAT *b, stateMPC_FLOAT *y)
{
    int i,j,ii,di;
    stateMPC_FLOAT yel;
            
    ii = 0; di = 0;
    for( i=0; i<9; i++ ){
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
 * where A is to be computed and is of size [9 x 9],
 * B is given and of size [9 x 9], L is a lower tri-
 * angular matrix of size 9 stored in lower triangular 
 * storage format. Note the transpose of L AND B!
 *
 * Result: A in column major storage format.
 *
 */
void stateMPC_LA_DENSE_MATRIXTFORWARDSUB_9_9(stateMPC_FLOAT *L, stateMPC_FLOAT *B, stateMPC_FLOAT *A)
{
    int i,j,k,ii,di;
    stateMPC_FLOAT a;
    
    ii=0; di=0;
    for( j=0; j<9; j++ ){        
        for( i=0; i<9; i++ ){
            a = B[i*9+j];
            for( k=0; k<j; k++ ){
                a -= A[k*9+i]*L[ii+k];
            }    

			/* saturate for numerical stability */
			a = MIN(a, BIGM);
			a = MAX(a, -BIGM); 

			A[j*9+i] = a/L[ii+j];			
        }
        ii += ++di;
    }
}


/**
 * Compute L = L - A*A', where L is lower triangular of size 9
 * and A is a dense matrix of size [9 x 9] in column
 * storage format.
 * 
 * THIS ONE HAS THE WORST ACCES PATTERN POSSIBLE. 
 * POSSIBKE FIX: PUT A INTO ROW MAJOR FORMAT FIRST.
 * 
 */
void stateMPC_LA_DENSE_MMTSUB_9_9(stateMPC_FLOAT *A, stateMPC_FLOAT *L)
{
    int i, j, k, ii, di;
    stateMPC_FLOAT ltemp;
    
    ii = 0; di = 0;
    for( i=0; i<9; i++ ){        
        for( j=0; j<=i; j++ ){
            ltemp = 0; 
            for( k=0; k<9; k++ ){
                ltemp += A[k*9+i]*A[k*9+j];
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
void stateMPC_LA_DENSE_MVMSUB1_9_9(stateMPC_FLOAT *A, stateMPC_FLOAT *x, stateMPC_FLOAT *b, stateMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;

	for( i=0; i<9; i++ ){
		r[i] = b[i] - A[k++]*x[0];
	}	
	for( j=1; j<9; j++ ){		
		for( i=0; i<9; i++ ){
			r[i] -= A[k++]*x[j];
		}
	}
}


/**
 * Backward Substitution to solve L^T*x = y where L is a
 * lower triangular matrix in triangular storage format.
 * 
 * All involved dimensions are 9.
 */
void stateMPC_LA_DENSE_BACKWARDSUB_9(stateMPC_FLOAT *L, stateMPC_FLOAT *y, stateMPC_FLOAT *x)
{
    int i, ii, di, j, jj, dj;
    stateMPC_FLOAT xel;    
	int start = 36;
    
    /* now solve L^T*x = y by backward substitution */
    ii = start; di = 8;
    for( i=8; i>=0; i-- ){        
        xel = y[i];        
        jj = start; dj = 8;
        for( j=8; j>i; j-- ){
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
 * Matrix vector multiplication y = b - M'*x where M is of size [9 x 9]
 * and stored in column major format. Note the transpose of M!
 */
void stateMPC_LA_DENSE_MTVMSUB_9_9(stateMPC_FLOAT *A, stateMPC_FLOAT *x, stateMPC_FLOAT *b, stateMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0; 
	for( i=0; i<9; i++ ){
		r[i] = b[i];
		for( j=0; j<9; j++ ){
			r[i] -= A[k++]*x[j];
		}
	}
}


/*
 * Vector subtraction z = -x - y for vectors of length 270.
 */
void stateMPC_LA_VSUB2_270(stateMPC_FLOAT *x, stateMPC_FLOAT *y, stateMPC_FLOAT *z)
{
	int i;
	for( i=0; i<270; i++){
		z[i] = -x[i] - y[i];
	}
}


/**
 * Forward-Backward-Substitution to solve L*L^T*x = b where L is a
 * diagonal matrix of size 29 in vector
 * storage format.
 */
void stateMPC_LA_DIAG_FORWARDBACKWARDSUB_29(stateMPC_FLOAT *L, stateMPC_FLOAT *b, stateMPC_FLOAT *x)
{
    int i;
            
    /* solve Ly = b by forward and backward substitution */
    for( i=0; i<29; i++ ){
		x[i] = b[i]/(L[i]*L[i]);
    }
    
}


/**
 * Forward-Backward-Substitution to solve L*L^T*x = b where L is a
 * diagonal matrix of size 9 in vector
 * storage format.
 */
void stateMPC_LA_DIAG_FORWARDBACKWARDSUB_9(stateMPC_FLOAT *L, stateMPC_FLOAT *b, stateMPC_FLOAT *x)
{
    int i;
            
    /* solve Ly = b by forward and backward substitution */
    for( i=0; i<9; i++ ){
		x[i] = b[i]/(L[i]*L[i]);
    }
    
}


/*
 * Vector subtraction z = x(xidx) - y where y, z and xidx are of length 29,
 * and x has length 29 and is indexed through yidx.
 */
void stateMPC_LA_VSUB_INDEXED_29(stateMPC_FLOAT *x, int* xidx, stateMPC_FLOAT *y, stateMPC_FLOAT *z)
{
	int i;
	for( i=0; i<29; i++){
		z[i] = x[xidx[i]] - y[i];
	}
}


/*
 * Vector subtraction x = -u.*v - w for vectors of length 29.
 */
void stateMPC_LA_VSUB3_29(stateMPC_FLOAT *u, stateMPC_FLOAT *v, stateMPC_FLOAT *w, stateMPC_FLOAT *x)
{
	int i;
	for( i=0; i<29; i++){
		x[i] = -u[i]*v[i] - w[i];
	}
}


/*
 * Vector subtraction z = -x - y(yidx) where y is of length 29
 * and z, x and yidx are of length 11.
 */
void stateMPC_LA_VSUB2_INDEXED_11(stateMPC_FLOAT *x, stateMPC_FLOAT *y, int* yidx, stateMPC_FLOAT *z)
{
	int i;
	for( i=0; i<11; i++){
		z[i] = -x[i] - y[yidx[i]];
	}
}


/*
 * Vector subtraction x = -u.*v - w for vectors of length 11.
 */
void stateMPC_LA_VSUB3_11(stateMPC_FLOAT *u, stateMPC_FLOAT *v, stateMPC_FLOAT *w, stateMPC_FLOAT *x)
{
	int i;
	for( i=0; i<11; i++){
		x[i] = -u[i]*v[i] - w[i];
	}
}


/*
 * Vector subtraction z = x(xidx) - y where y, z and xidx are of length 9,
 * and x has length 9 and is indexed through yidx.
 */
void stateMPC_LA_VSUB_INDEXED_9(stateMPC_FLOAT *x, int* xidx, stateMPC_FLOAT *y, stateMPC_FLOAT *z)
{
	int i;
	for( i=0; i<9; i++){
		z[i] = x[xidx[i]] - y[i];
	}
}


/*
 * Vector subtraction x = -u.*v - w for vectors of length 9.
 */
void stateMPC_LA_VSUB3_9(stateMPC_FLOAT *u, stateMPC_FLOAT *v, stateMPC_FLOAT *w, stateMPC_FLOAT *x)
{
	int i;
	for( i=0; i<9; i++){
		x[i] = -u[i]*v[i] - w[i];
	}
}


/*
 * Vector subtraction z = -x - y(yidx) where y is of length 9
 * and z, x and yidx are of length 9.
 */
void stateMPC_LA_VSUB2_INDEXED_9(stateMPC_FLOAT *x, stateMPC_FLOAT *y, int* yidx, stateMPC_FLOAT *z)
{
	int i;
	for( i=0; i<9; i++){
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
        for( i=0; i<378; i++ ){
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
        if( i == 378 ){
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
    *mu_aff = mymu / (stateMPC_FLOAT)378;
    return lsIt;
}


/*
 * Vector subtraction x = u.*v - a where a is a scalar
*  and x,u,v are vectors of length 378.
 */
void stateMPC_LA_VSUB5_378(stateMPC_FLOAT *u, stateMPC_FLOAT *v, stateMPC_FLOAT a, stateMPC_FLOAT *x)
{
	int i;
	for( i=0; i<378; i++){
		x[i] = u[i]*v[i] - a;
	}
}


/*
 * Computes x=0; x(uidx) += u/su; x(vidx) -= v/sv where x is of length 29,
 * u, su, uidx are of length 11 and v, sv, vidx are of length 29.
 */
void stateMPC_LA_VSUB6_INDEXED_29_11_29(stateMPC_FLOAT *u, stateMPC_FLOAT *su, int* uidx, stateMPC_FLOAT *v, stateMPC_FLOAT *sv, int* vidx, stateMPC_FLOAT *x)
{
	int i;
	for( i=0; i<29; i++ ){
		x[i] = 0;
	}
	for( i=0; i<11; i++){
		x[uidx[i]] += u[i]/su[i];
	}
	for( i=0; i<29; i++){
		x[vidx[i]] -= v[i]/sv[i];
	}
}


/* 
 * Computes r =  B*u
 * where B is stored in diagzero format
 */
void stateMPC_LA_DIAGZERO_MVM_9(stateMPC_FLOAT *B, stateMPC_FLOAT *u, stateMPC_FLOAT *r)
{
	int i;

	for( i=0; i<9; i++ ){
		r[i] = B[i]*u[i];
	}	
	
}


/* 
 * Computes r = A*x + B*u
 * where A is stored in column major format
 * and B is stored in diagzero format
 */
void stateMPC_LA_DENSE_DIAGZERO_2MVMADD_9_29_29(stateMPC_FLOAT *A, stateMPC_FLOAT *x, stateMPC_FLOAT *B, stateMPC_FLOAT *u, stateMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;

	for( i=0; i<9; i++ ){
		r[i] = A[k++]*x[0] + B[i]*u[i];
	}	

	for( j=1; j<29; j++ ){		
		for( i=0; i<9; i++ ){
			r[i] += A[k++]*x[j];
		}
	}
	
}


/*
 * Computes x=0; x(uidx) += u/su; x(vidx) -= v/sv where x is of length 9,
 * u, su, uidx are of length 9 and v, sv, vidx are of length 9.
 */
void stateMPC_LA_VSUB6_INDEXED_9_9_9(stateMPC_FLOAT *u, stateMPC_FLOAT *su, int* uidx, stateMPC_FLOAT *v, stateMPC_FLOAT *sv, int* vidx, stateMPC_FLOAT *x)
{
	int i;
	for( i=0; i<9; i++ ){
		x[i] = 0;
	}
	for( i=0; i<9; i++){
		x[uidx[i]] += u[i]/su[i];
	}
	for( i=0; i<9; i++){
		x[vidx[i]] -= v[i]/sv[i];
	}
}


/* 
 * Computes r = A*x + B*u
 * where A is stored in column major format
 * and B is stored in diagzero format
 */
void stateMPC_LA_DENSE_DIAGZERO_2MVMADD_9_29_9(stateMPC_FLOAT *A, stateMPC_FLOAT *x, stateMPC_FLOAT *B, stateMPC_FLOAT *u, stateMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;

	for( i=0; i<9; i++ ){
		r[i] = A[k++]*x[0] + B[i]*u[i];
	}	

	for( j=1; j<29; j++ ){		
		for( i=0; i<9; i++ ){
			r[i] += A[k++]*x[j];
		}
	}
	
}


/*
 * Vector subtraction z = x - y for vectors of length 270.
 */
void stateMPC_LA_VSUB_270(stateMPC_FLOAT *x, stateMPC_FLOAT *y, stateMPC_FLOAT *z)
{
	int i;
	for( i=0; i<270; i++){
		z[i] = x[i] - y[i];
	}
}


/** 
 * Computes z = -r./s - u.*y(y)
 * where all vectors except of y are of length 29 (length of y >= 29).
 */
void stateMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_29(stateMPC_FLOAT *r, stateMPC_FLOAT *s, stateMPC_FLOAT *u, stateMPC_FLOAT *y, int* yidx, stateMPC_FLOAT *z)
{
	int i;
	for( i=0; i<29; i++ ){
		z[i] = -r[i]/s[i] - u[i]*y[yidx[i]];
	}
}


/** 
 * Computes z = -r./s + u.*y(y)
 * where all vectors except of y are of length 11 (length of y >= 11).
 */
void stateMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_11(stateMPC_FLOAT *r, stateMPC_FLOAT *s, stateMPC_FLOAT *u, stateMPC_FLOAT *y, int* yidx, stateMPC_FLOAT *z)
{
	int i;
	for( i=0; i<11; i++ ){
		z[i] = -r[i]/s[i] + u[i]*y[yidx[i]];
	}
}


/** 
 * Computes z = -r./s - u.*y(y)
 * where all vectors except of y are of length 9 (length of y >= 9).
 */
void stateMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_9(stateMPC_FLOAT *r, stateMPC_FLOAT *s, stateMPC_FLOAT *u, stateMPC_FLOAT *y, int* yidx, stateMPC_FLOAT *z)
{
	int i;
	for( i=0; i<9; i++ ){
		z[i] = -r[i]/s[i] - u[i]*y[yidx[i]];
	}
}


/** 
 * Computes z = -r./s + u.*y(y)
 * where all vectors except of y are of length 9 (length of y >= 9).
 */
void stateMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_9(stateMPC_FLOAT *r, stateMPC_FLOAT *s, stateMPC_FLOAT *u, stateMPC_FLOAT *y, int* yidx, stateMPC_FLOAT *z)
{
	int i;
	for( i=0; i<9; i++ ){
		z[i] = -r[i]/s[i] + u[i]*y[yidx[i]];
	}
}


/*
 * Computes ds = -l.\(r + s.*dl) for vectors of length 378.
 */
void stateMPC_LA_VSUB7_378(stateMPC_FLOAT *l, stateMPC_FLOAT *r, stateMPC_FLOAT *s, stateMPC_FLOAT *dl, stateMPC_FLOAT *ds)
{
	int i;
	for( i=0; i<378; i++){
		ds[i] = -(r[i] + s[i]*dl[i])/l[i];
	}
}


/*
 * Vector addition x = x + y for vectors of length 270.
 */
void stateMPC_LA_VADD_270(stateMPC_FLOAT *x, stateMPC_FLOAT *y)
{
	int i;
	for( i=0; i<270; i++){
		x[i] += y[i];
	}
}


/*
 * Vector addition x = x + y for vectors of length 90.
 */
void stateMPC_LA_VADD_90(stateMPC_FLOAT *x, stateMPC_FLOAT *y)
{
	int i;
	for( i=0; i<90; i++){
		x[i] += y[i];
	}
}


/*
 * Vector addition x = x + y for vectors of length 378.
 */
void stateMPC_LA_VADD_378(stateMPC_FLOAT *x, stateMPC_FLOAT *y)
{
	int i;
	for( i=0; i<378; i++){
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
        for( i=0; i<378; i++ ){
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
        if( i == 378 ){
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
    for( i=0; i<270; i++ ){
        z[i] += a_gamma*dz[i];
    }
    
    /* equality constraint multipliers */
    for( i=0; i<90; i++ ){
        v[i] += a_gamma*dv[i];
    }
    
    /* inequality constraint multipliers & slacks, also update mu */
    *mu = 0;
    for( i=0; i<378; i++ ){
        dltemp = l[i] + a_gamma*dl[i]; l[i] = dltemp;
        dstemp = s[i] + a_gamma*ds[i]; s[i] = dstemp;
        *mu += dltemp*dstemp;
    }
    
    *a = a_gamma;
    *mu /= (stateMPC_FLOAT)378;
    return lsIt;
}




/* VARIABLE DEFINITIONS ------------------------------------------------ */
stateMPC_FLOAT stateMPC_z[270];
stateMPC_FLOAT stateMPC_v[90];
stateMPC_FLOAT stateMPC_dz_aff[270];
stateMPC_FLOAT stateMPC_dv_aff[90];
stateMPC_FLOAT stateMPC_grad_cost[270];
stateMPC_FLOAT stateMPC_grad_eq[270];
stateMPC_FLOAT stateMPC_rd[270];
stateMPC_FLOAT stateMPC_l[378];
stateMPC_FLOAT stateMPC_s[378];
stateMPC_FLOAT stateMPC_lbys[378];
stateMPC_FLOAT stateMPC_dl_aff[378];
stateMPC_FLOAT stateMPC_ds_aff[378];
stateMPC_FLOAT stateMPC_dz_cc[270];
stateMPC_FLOAT stateMPC_dv_cc[90];
stateMPC_FLOAT stateMPC_dl_cc[378];
stateMPC_FLOAT stateMPC_ds_cc[378];
stateMPC_FLOAT stateMPC_ccrhs[378];
stateMPC_FLOAT stateMPC_grad_ineq[270];
stateMPC_FLOAT* stateMPC_z0 = stateMPC_z + 0;
stateMPC_FLOAT* stateMPC_dzaff0 = stateMPC_dz_aff + 0;
stateMPC_FLOAT* stateMPC_dzcc0 = stateMPC_dz_cc + 0;
stateMPC_FLOAT* stateMPC_rd0 = stateMPC_rd + 0;
stateMPC_FLOAT stateMPC_Lbyrd0[29];
stateMPC_FLOAT* stateMPC_grad_cost0 = stateMPC_grad_cost + 0;
stateMPC_FLOAT* stateMPC_grad_eq0 = stateMPC_grad_eq + 0;
stateMPC_FLOAT* stateMPC_grad_ineq0 = stateMPC_grad_ineq + 0;
stateMPC_FLOAT stateMPC_ctv0[29];
stateMPC_FLOAT* stateMPC_v0 = stateMPC_v + 0;
stateMPC_FLOAT stateMPC_re0[9];
stateMPC_FLOAT stateMPC_beta0[9];
stateMPC_FLOAT stateMPC_betacc0[9];
stateMPC_FLOAT* stateMPC_dvaff0 = stateMPC_dv_aff + 0;
stateMPC_FLOAT* stateMPC_dvcc0 = stateMPC_dv_cc + 0;
stateMPC_FLOAT stateMPC_V0[261];
stateMPC_FLOAT stateMPC_Yd0[45];
stateMPC_FLOAT stateMPC_Ld0[45];
stateMPC_FLOAT stateMPC_yy0[9];
stateMPC_FLOAT stateMPC_bmy0[9];
int stateMPC_lbIdx0[29] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28};
stateMPC_FLOAT* stateMPC_llb0 = stateMPC_l + 0;
stateMPC_FLOAT* stateMPC_slb0 = stateMPC_s + 0;
stateMPC_FLOAT* stateMPC_llbbyslb0 = stateMPC_lbys + 0;
stateMPC_FLOAT stateMPC_rilb0[29];
stateMPC_FLOAT* stateMPC_dllbaff0 = stateMPC_dl_aff + 0;
stateMPC_FLOAT* stateMPC_dslbaff0 = stateMPC_ds_aff + 0;
stateMPC_FLOAT* stateMPC_dllbcc0 = stateMPC_dl_cc + 0;
stateMPC_FLOAT* stateMPC_dslbcc0 = stateMPC_ds_cc + 0;
stateMPC_FLOAT* stateMPC_ccrhsl0 = stateMPC_ccrhs + 0;
int stateMPC_ubIdx0[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
stateMPC_FLOAT* stateMPC_lub0 = stateMPC_l + 29;
stateMPC_FLOAT* stateMPC_sub0 = stateMPC_s + 29;
stateMPC_FLOAT* stateMPC_lubbysub0 = stateMPC_lbys + 29;
stateMPC_FLOAT stateMPC_riub0[11];
stateMPC_FLOAT* stateMPC_dlubaff0 = stateMPC_dl_aff + 29;
stateMPC_FLOAT* stateMPC_dsubaff0 = stateMPC_ds_aff + 29;
stateMPC_FLOAT* stateMPC_dlubcc0 = stateMPC_dl_cc + 29;
stateMPC_FLOAT* stateMPC_dsubcc0 = stateMPC_ds_cc + 29;
stateMPC_FLOAT* stateMPC_ccrhsub0 = stateMPC_ccrhs + 29;
stateMPC_FLOAT stateMPC_Phi0[29];
stateMPC_FLOAT stateMPC_D0[29] = {1.0000000000000000E+000, 
1.0000000000000000E+000, 
1.0000000000000000E+000, 
1.0000000000000000E+000, 
1.0000000000000000E+000, 
1.0000000000000000E+000, 
1.0000000000000000E+000, 
1.0000000000000000E+000, 
1.0000000000000000E+000};
stateMPC_FLOAT stateMPC_W0[29];
stateMPC_FLOAT* stateMPC_z1 = stateMPC_z + 29;
stateMPC_FLOAT* stateMPC_dzaff1 = stateMPC_dz_aff + 29;
stateMPC_FLOAT* stateMPC_dzcc1 = stateMPC_dz_cc + 29;
stateMPC_FLOAT* stateMPC_rd1 = stateMPC_rd + 29;
stateMPC_FLOAT stateMPC_Lbyrd1[29];
stateMPC_FLOAT* stateMPC_grad_cost1 = stateMPC_grad_cost + 29;
stateMPC_FLOAT* stateMPC_grad_eq1 = stateMPC_grad_eq + 29;
stateMPC_FLOAT* stateMPC_grad_ineq1 = stateMPC_grad_ineq + 29;
stateMPC_FLOAT stateMPC_ctv1[29];
stateMPC_FLOAT* stateMPC_v1 = stateMPC_v + 9;
stateMPC_FLOAT stateMPC_re1[9];
stateMPC_FLOAT stateMPC_beta1[9];
stateMPC_FLOAT stateMPC_betacc1[9];
stateMPC_FLOAT* stateMPC_dvaff1 = stateMPC_dv_aff + 9;
stateMPC_FLOAT* stateMPC_dvcc1 = stateMPC_dv_cc + 9;
stateMPC_FLOAT stateMPC_V1[261];
stateMPC_FLOAT stateMPC_Yd1[45];
stateMPC_FLOAT stateMPC_Ld1[45];
stateMPC_FLOAT stateMPC_yy1[9];
stateMPC_FLOAT stateMPC_bmy1[9];
int stateMPC_lbIdx1[29] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28};
stateMPC_FLOAT* stateMPC_llb1 = stateMPC_l + 40;
stateMPC_FLOAT* stateMPC_slb1 = stateMPC_s + 40;
stateMPC_FLOAT* stateMPC_llbbyslb1 = stateMPC_lbys + 40;
stateMPC_FLOAT stateMPC_rilb1[29];
stateMPC_FLOAT* stateMPC_dllbaff1 = stateMPC_dl_aff + 40;
stateMPC_FLOAT* stateMPC_dslbaff1 = stateMPC_ds_aff + 40;
stateMPC_FLOAT* stateMPC_dllbcc1 = stateMPC_dl_cc + 40;
stateMPC_FLOAT* stateMPC_dslbcc1 = stateMPC_ds_cc + 40;
stateMPC_FLOAT* stateMPC_ccrhsl1 = stateMPC_ccrhs + 40;
int stateMPC_ubIdx1[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
stateMPC_FLOAT* stateMPC_lub1 = stateMPC_l + 69;
stateMPC_FLOAT* stateMPC_sub1 = stateMPC_s + 69;
stateMPC_FLOAT* stateMPC_lubbysub1 = stateMPC_lbys + 69;
stateMPC_FLOAT stateMPC_riub1[11];
stateMPC_FLOAT* stateMPC_dlubaff1 = stateMPC_dl_aff + 69;
stateMPC_FLOAT* stateMPC_dsubaff1 = stateMPC_ds_aff + 69;
stateMPC_FLOAT* stateMPC_dlubcc1 = stateMPC_dl_cc + 69;
stateMPC_FLOAT* stateMPC_dsubcc1 = stateMPC_ds_cc + 69;
stateMPC_FLOAT* stateMPC_ccrhsub1 = stateMPC_ccrhs + 69;
stateMPC_FLOAT stateMPC_Phi1[29];
stateMPC_FLOAT stateMPC_D1[29] = {-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000};
stateMPC_FLOAT stateMPC_W1[29];
stateMPC_FLOAT stateMPC_Ysd1[81];
stateMPC_FLOAT stateMPC_Lsd1[81];
stateMPC_FLOAT* stateMPC_z2 = stateMPC_z + 58;
stateMPC_FLOAT* stateMPC_dzaff2 = stateMPC_dz_aff + 58;
stateMPC_FLOAT* stateMPC_dzcc2 = stateMPC_dz_cc + 58;
stateMPC_FLOAT* stateMPC_rd2 = stateMPC_rd + 58;
stateMPC_FLOAT stateMPC_Lbyrd2[29];
stateMPC_FLOAT* stateMPC_grad_cost2 = stateMPC_grad_cost + 58;
stateMPC_FLOAT* stateMPC_grad_eq2 = stateMPC_grad_eq + 58;
stateMPC_FLOAT* stateMPC_grad_ineq2 = stateMPC_grad_ineq + 58;
stateMPC_FLOAT stateMPC_ctv2[29];
stateMPC_FLOAT* stateMPC_v2 = stateMPC_v + 18;
stateMPC_FLOAT stateMPC_re2[9];
stateMPC_FLOAT stateMPC_beta2[9];
stateMPC_FLOAT stateMPC_betacc2[9];
stateMPC_FLOAT* stateMPC_dvaff2 = stateMPC_dv_aff + 18;
stateMPC_FLOAT* stateMPC_dvcc2 = stateMPC_dv_cc + 18;
stateMPC_FLOAT stateMPC_V2[261];
stateMPC_FLOAT stateMPC_Yd2[45];
stateMPC_FLOAT stateMPC_Ld2[45];
stateMPC_FLOAT stateMPC_yy2[9];
stateMPC_FLOAT stateMPC_bmy2[9];
int stateMPC_lbIdx2[29] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28};
stateMPC_FLOAT* stateMPC_llb2 = stateMPC_l + 80;
stateMPC_FLOAT* stateMPC_slb2 = stateMPC_s + 80;
stateMPC_FLOAT* stateMPC_llbbyslb2 = stateMPC_lbys + 80;
stateMPC_FLOAT stateMPC_rilb2[29];
stateMPC_FLOAT* stateMPC_dllbaff2 = stateMPC_dl_aff + 80;
stateMPC_FLOAT* stateMPC_dslbaff2 = stateMPC_ds_aff + 80;
stateMPC_FLOAT* stateMPC_dllbcc2 = stateMPC_dl_cc + 80;
stateMPC_FLOAT* stateMPC_dslbcc2 = stateMPC_ds_cc + 80;
stateMPC_FLOAT* stateMPC_ccrhsl2 = stateMPC_ccrhs + 80;
int stateMPC_ubIdx2[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
stateMPC_FLOAT* stateMPC_lub2 = stateMPC_l + 109;
stateMPC_FLOAT* stateMPC_sub2 = stateMPC_s + 109;
stateMPC_FLOAT* stateMPC_lubbysub2 = stateMPC_lbys + 109;
stateMPC_FLOAT stateMPC_riub2[11];
stateMPC_FLOAT* stateMPC_dlubaff2 = stateMPC_dl_aff + 109;
stateMPC_FLOAT* stateMPC_dsubaff2 = stateMPC_ds_aff + 109;
stateMPC_FLOAT* stateMPC_dlubcc2 = stateMPC_dl_cc + 109;
stateMPC_FLOAT* stateMPC_dsubcc2 = stateMPC_ds_cc + 109;
stateMPC_FLOAT* stateMPC_ccrhsub2 = stateMPC_ccrhs + 109;
stateMPC_FLOAT stateMPC_Phi2[29];
stateMPC_FLOAT stateMPC_W2[29];
stateMPC_FLOAT stateMPC_Ysd2[81];
stateMPC_FLOAT stateMPC_Lsd2[81];
stateMPC_FLOAT* stateMPC_z3 = stateMPC_z + 87;
stateMPC_FLOAT* stateMPC_dzaff3 = stateMPC_dz_aff + 87;
stateMPC_FLOAT* stateMPC_dzcc3 = stateMPC_dz_cc + 87;
stateMPC_FLOAT* stateMPC_rd3 = stateMPC_rd + 87;
stateMPC_FLOAT stateMPC_Lbyrd3[29];
stateMPC_FLOAT* stateMPC_grad_cost3 = stateMPC_grad_cost + 87;
stateMPC_FLOAT* stateMPC_grad_eq3 = stateMPC_grad_eq + 87;
stateMPC_FLOAT* stateMPC_grad_ineq3 = stateMPC_grad_ineq + 87;
stateMPC_FLOAT stateMPC_ctv3[29];
stateMPC_FLOAT* stateMPC_v3 = stateMPC_v + 27;
stateMPC_FLOAT stateMPC_re3[9];
stateMPC_FLOAT stateMPC_beta3[9];
stateMPC_FLOAT stateMPC_betacc3[9];
stateMPC_FLOAT* stateMPC_dvaff3 = stateMPC_dv_aff + 27;
stateMPC_FLOAT* stateMPC_dvcc3 = stateMPC_dv_cc + 27;
stateMPC_FLOAT stateMPC_V3[261];
stateMPC_FLOAT stateMPC_Yd3[45];
stateMPC_FLOAT stateMPC_Ld3[45];
stateMPC_FLOAT stateMPC_yy3[9];
stateMPC_FLOAT stateMPC_bmy3[9];
int stateMPC_lbIdx3[29] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28};
stateMPC_FLOAT* stateMPC_llb3 = stateMPC_l + 120;
stateMPC_FLOAT* stateMPC_slb3 = stateMPC_s + 120;
stateMPC_FLOAT* stateMPC_llbbyslb3 = stateMPC_lbys + 120;
stateMPC_FLOAT stateMPC_rilb3[29];
stateMPC_FLOAT* stateMPC_dllbaff3 = stateMPC_dl_aff + 120;
stateMPC_FLOAT* stateMPC_dslbaff3 = stateMPC_ds_aff + 120;
stateMPC_FLOAT* stateMPC_dllbcc3 = stateMPC_dl_cc + 120;
stateMPC_FLOAT* stateMPC_dslbcc3 = stateMPC_ds_cc + 120;
stateMPC_FLOAT* stateMPC_ccrhsl3 = stateMPC_ccrhs + 120;
int stateMPC_ubIdx3[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
stateMPC_FLOAT* stateMPC_lub3 = stateMPC_l + 149;
stateMPC_FLOAT* stateMPC_sub3 = stateMPC_s + 149;
stateMPC_FLOAT* stateMPC_lubbysub3 = stateMPC_lbys + 149;
stateMPC_FLOAT stateMPC_riub3[11];
stateMPC_FLOAT* stateMPC_dlubaff3 = stateMPC_dl_aff + 149;
stateMPC_FLOAT* stateMPC_dsubaff3 = stateMPC_ds_aff + 149;
stateMPC_FLOAT* stateMPC_dlubcc3 = stateMPC_dl_cc + 149;
stateMPC_FLOAT* stateMPC_dsubcc3 = stateMPC_ds_cc + 149;
stateMPC_FLOAT* stateMPC_ccrhsub3 = stateMPC_ccrhs + 149;
stateMPC_FLOAT stateMPC_Phi3[29];
stateMPC_FLOAT stateMPC_W3[29];
stateMPC_FLOAT stateMPC_Ysd3[81];
stateMPC_FLOAT stateMPC_Lsd3[81];
stateMPC_FLOAT* stateMPC_z4 = stateMPC_z + 116;
stateMPC_FLOAT* stateMPC_dzaff4 = stateMPC_dz_aff + 116;
stateMPC_FLOAT* stateMPC_dzcc4 = stateMPC_dz_cc + 116;
stateMPC_FLOAT* stateMPC_rd4 = stateMPC_rd + 116;
stateMPC_FLOAT stateMPC_Lbyrd4[29];
stateMPC_FLOAT* stateMPC_grad_cost4 = stateMPC_grad_cost + 116;
stateMPC_FLOAT* stateMPC_grad_eq4 = stateMPC_grad_eq + 116;
stateMPC_FLOAT* stateMPC_grad_ineq4 = stateMPC_grad_ineq + 116;
stateMPC_FLOAT stateMPC_ctv4[29];
stateMPC_FLOAT* stateMPC_v4 = stateMPC_v + 36;
stateMPC_FLOAT stateMPC_re4[9];
stateMPC_FLOAT stateMPC_beta4[9];
stateMPC_FLOAT stateMPC_betacc4[9];
stateMPC_FLOAT* stateMPC_dvaff4 = stateMPC_dv_aff + 36;
stateMPC_FLOAT* stateMPC_dvcc4 = stateMPC_dv_cc + 36;
stateMPC_FLOAT stateMPC_V4[261];
stateMPC_FLOAT stateMPC_Yd4[45];
stateMPC_FLOAT stateMPC_Ld4[45];
stateMPC_FLOAT stateMPC_yy4[9];
stateMPC_FLOAT stateMPC_bmy4[9];
int stateMPC_lbIdx4[29] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28};
stateMPC_FLOAT* stateMPC_llb4 = stateMPC_l + 160;
stateMPC_FLOAT* stateMPC_slb4 = stateMPC_s + 160;
stateMPC_FLOAT* stateMPC_llbbyslb4 = stateMPC_lbys + 160;
stateMPC_FLOAT stateMPC_rilb4[29];
stateMPC_FLOAT* stateMPC_dllbaff4 = stateMPC_dl_aff + 160;
stateMPC_FLOAT* stateMPC_dslbaff4 = stateMPC_ds_aff + 160;
stateMPC_FLOAT* stateMPC_dllbcc4 = stateMPC_dl_cc + 160;
stateMPC_FLOAT* stateMPC_dslbcc4 = stateMPC_ds_cc + 160;
stateMPC_FLOAT* stateMPC_ccrhsl4 = stateMPC_ccrhs + 160;
int stateMPC_ubIdx4[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
stateMPC_FLOAT* stateMPC_lub4 = stateMPC_l + 189;
stateMPC_FLOAT* stateMPC_sub4 = stateMPC_s + 189;
stateMPC_FLOAT* stateMPC_lubbysub4 = stateMPC_lbys + 189;
stateMPC_FLOAT stateMPC_riub4[11];
stateMPC_FLOAT* stateMPC_dlubaff4 = stateMPC_dl_aff + 189;
stateMPC_FLOAT* stateMPC_dsubaff4 = stateMPC_ds_aff + 189;
stateMPC_FLOAT* stateMPC_dlubcc4 = stateMPC_dl_cc + 189;
stateMPC_FLOAT* stateMPC_dsubcc4 = stateMPC_ds_cc + 189;
stateMPC_FLOAT* stateMPC_ccrhsub4 = stateMPC_ccrhs + 189;
stateMPC_FLOAT stateMPC_Phi4[29];
stateMPC_FLOAT stateMPC_W4[29];
stateMPC_FLOAT stateMPC_Ysd4[81];
stateMPC_FLOAT stateMPC_Lsd4[81];
stateMPC_FLOAT* stateMPC_z5 = stateMPC_z + 145;
stateMPC_FLOAT* stateMPC_dzaff5 = stateMPC_dz_aff + 145;
stateMPC_FLOAT* stateMPC_dzcc5 = stateMPC_dz_cc + 145;
stateMPC_FLOAT* stateMPC_rd5 = stateMPC_rd + 145;
stateMPC_FLOAT stateMPC_Lbyrd5[29];
stateMPC_FLOAT* stateMPC_grad_cost5 = stateMPC_grad_cost + 145;
stateMPC_FLOAT* stateMPC_grad_eq5 = stateMPC_grad_eq + 145;
stateMPC_FLOAT* stateMPC_grad_ineq5 = stateMPC_grad_ineq + 145;
stateMPC_FLOAT stateMPC_ctv5[29];
stateMPC_FLOAT* stateMPC_v5 = stateMPC_v + 45;
stateMPC_FLOAT stateMPC_re5[9];
stateMPC_FLOAT stateMPC_beta5[9];
stateMPC_FLOAT stateMPC_betacc5[9];
stateMPC_FLOAT* stateMPC_dvaff5 = stateMPC_dv_aff + 45;
stateMPC_FLOAT* stateMPC_dvcc5 = stateMPC_dv_cc + 45;
stateMPC_FLOAT stateMPC_V5[261];
stateMPC_FLOAT stateMPC_Yd5[45];
stateMPC_FLOAT stateMPC_Ld5[45];
stateMPC_FLOAT stateMPC_yy5[9];
stateMPC_FLOAT stateMPC_bmy5[9];
int stateMPC_lbIdx5[29] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28};
stateMPC_FLOAT* stateMPC_llb5 = stateMPC_l + 200;
stateMPC_FLOAT* stateMPC_slb5 = stateMPC_s + 200;
stateMPC_FLOAT* stateMPC_llbbyslb5 = stateMPC_lbys + 200;
stateMPC_FLOAT stateMPC_rilb5[29];
stateMPC_FLOAT* stateMPC_dllbaff5 = stateMPC_dl_aff + 200;
stateMPC_FLOAT* stateMPC_dslbaff5 = stateMPC_ds_aff + 200;
stateMPC_FLOAT* stateMPC_dllbcc5 = stateMPC_dl_cc + 200;
stateMPC_FLOAT* stateMPC_dslbcc5 = stateMPC_ds_cc + 200;
stateMPC_FLOAT* stateMPC_ccrhsl5 = stateMPC_ccrhs + 200;
int stateMPC_ubIdx5[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
stateMPC_FLOAT* stateMPC_lub5 = stateMPC_l + 229;
stateMPC_FLOAT* stateMPC_sub5 = stateMPC_s + 229;
stateMPC_FLOAT* stateMPC_lubbysub5 = stateMPC_lbys + 229;
stateMPC_FLOAT stateMPC_riub5[11];
stateMPC_FLOAT* stateMPC_dlubaff5 = stateMPC_dl_aff + 229;
stateMPC_FLOAT* stateMPC_dsubaff5 = stateMPC_ds_aff + 229;
stateMPC_FLOAT* stateMPC_dlubcc5 = stateMPC_dl_cc + 229;
stateMPC_FLOAT* stateMPC_dsubcc5 = stateMPC_ds_cc + 229;
stateMPC_FLOAT* stateMPC_ccrhsub5 = stateMPC_ccrhs + 229;
stateMPC_FLOAT stateMPC_Phi5[29];
stateMPC_FLOAT stateMPC_W5[29];
stateMPC_FLOAT stateMPC_Ysd5[81];
stateMPC_FLOAT stateMPC_Lsd5[81];
stateMPC_FLOAT* stateMPC_z6 = stateMPC_z + 174;
stateMPC_FLOAT* stateMPC_dzaff6 = stateMPC_dz_aff + 174;
stateMPC_FLOAT* stateMPC_dzcc6 = stateMPC_dz_cc + 174;
stateMPC_FLOAT* stateMPC_rd6 = stateMPC_rd + 174;
stateMPC_FLOAT stateMPC_Lbyrd6[29];
stateMPC_FLOAT* stateMPC_grad_cost6 = stateMPC_grad_cost + 174;
stateMPC_FLOAT* stateMPC_grad_eq6 = stateMPC_grad_eq + 174;
stateMPC_FLOAT* stateMPC_grad_ineq6 = stateMPC_grad_ineq + 174;
stateMPC_FLOAT stateMPC_ctv6[29];
stateMPC_FLOAT* stateMPC_v6 = stateMPC_v + 54;
stateMPC_FLOAT stateMPC_re6[9];
stateMPC_FLOAT stateMPC_beta6[9];
stateMPC_FLOAT stateMPC_betacc6[9];
stateMPC_FLOAT* stateMPC_dvaff6 = stateMPC_dv_aff + 54;
stateMPC_FLOAT* stateMPC_dvcc6 = stateMPC_dv_cc + 54;
stateMPC_FLOAT stateMPC_V6[261];
stateMPC_FLOAT stateMPC_Yd6[45];
stateMPC_FLOAT stateMPC_Ld6[45];
stateMPC_FLOAT stateMPC_yy6[9];
stateMPC_FLOAT stateMPC_bmy6[9];
int stateMPC_lbIdx6[29] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28};
stateMPC_FLOAT* stateMPC_llb6 = stateMPC_l + 240;
stateMPC_FLOAT* stateMPC_slb6 = stateMPC_s + 240;
stateMPC_FLOAT* stateMPC_llbbyslb6 = stateMPC_lbys + 240;
stateMPC_FLOAT stateMPC_rilb6[29];
stateMPC_FLOAT* stateMPC_dllbaff6 = stateMPC_dl_aff + 240;
stateMPC_FLOAT* stateMPC_dslbaff6 = stateMPC_ds_aff + 240;
stateMPC_FLOAT* stateMPC_dllbcc6 = stateMPC_dl_cc + 240;
stateMPC_FLOAT* stateMPC_dslbcc6 = stateMPC_ds_cc + 240;
stateMPC_FLOAT* stateMPC_ccrhsl6 = stateMPC_ccrhs + 240;
int stateMPC_ubIdx6[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
stateMPC_FLOAT* stateMPC_lub6 = stateMPC_l + 269;
stateMPC_FLOAT* stateMPC_sub6 = stateMPC_s + 269;
stateMPC_FLOAT* stateMPC_lubbysub6 = stateMPC_lbys + 269;
stateMPC_FLOAT stateMPC_riub6[11];
stateMPC_FLOAT* stateMPC_dlubaff6 = stateMPC_dl_aff + 269;
stateMPC_FLOAT* stateMPC_dsubaff6 = stateMPC_ds_aff + 269;
stateMPC_FLOAT* stateMPC_dlubcc6 = stateMPC_dl_cc + 269;
stateMPC_FLOAT* stateMPC_dsubcc6 = stateMPC_ds_cc + 269;
stateMPC_FLOAT* stateMPC_ccrhsub6 = stateMPC_ccrhs + 269;
stateMPC_FLOAT stateMPC_Phi6[29];
stateMPC_FLOAT stateMPC_W6[29];
stateMPC_FLOAT stateMPC_Ysd6[81];
stateMPC_FLOAT stateMPC_Lsd6[81];
stateMPC_FLOAT* stateMPC_z7 = stateMPC_z + 203;
stateMPC_FLOAT* stateMPC_dzaff7 = stateMPC_dz_aff + 203;
stateMPC_FLOAT* stateMPC_dzcc7 = stateMPC_dz_cc + 203;
stateMPC_FLOAT* stateMPC_rd7 = stateMPC_rd + 203;
stateMPC_FLOAT stateMPC_Lbyrd7[29];
stateMPC_FLOAT* stateMPC_grad_cost7 = stateMPC_grad_cost + 203;
stateMPC_FLOAT* stateMPC_grad_eq7 = stateMPC_grad_eq + 203;
stateMPC_FLOAT* stateMPC_grad_ineq7 = stateMPC_grad_ineq + 203;
stateMPC_FLOAT stateMPC_ctv7[29];
stateMPC_FLOAT* stateMPC_v7 = stateMPC_v + 63;
stateMPC_FLOAT stateMPC_re7[9];
stateMPC_FLOAT stateMPC_beta7[9];
stateMPC_FLOAT stateMPC_betacc7[9];
stateMPC_FLOAT* stateMPC_dvaff7 = stateMPC_dv_aff + 63;
stateMPC_FLOAT* stateMPC_dvcc7 = stateMPC_dv_cc + 63;
stateMPC_FLOAT stateMPC_V7[261];
stateMPC_FLOAT stateMPC_Yd7[45];
stateMPC_FLOAT stateMPC_Ld7[45];
stateMPC_FLOAT stateMPC_yy7[9];
stateMPC_FLOAT stateMPC_bmy7[9];
int stateMPC_lbIdx7[29] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28};
stateMPC_FLOAT* stateMPC_llb7 = stateMPC_l + 280;
stateMPC_FLOAT* stateMPC_slb7 = stateMPC_s + 280;
stateMPC_FLOAT* stateMPC_llbbyslb7 = stateMPC_lbys + 280;
stateMPC_FLOAT stateMPC_rilb7[29];
stateMPC_FLOAT* stateMPC_dllbaff7 = stateMPC_dl_aff + 280;
stateMPC_FLOAT* stateMPC_dslbaff7 = stateMPC_ds_aff + 280;
stateMPC_FLOAT* stateMPC_dllbcc7 = stateMPC_dl_cc + 280;
stateMPC_FLOAT* stateMPC_dslbcc7 = stateMPC_ds_cc + 280;
stateMPC_FLOAT* stateMPC_ccrhsl7 = stateMPC_ccrhs + 280;
int stateMPC_ubIdx7[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
stateMPC_FLOAT* stateMPC_lub7 = stateMPC_l + 309;
stateMPC_FLOAT* stateMPC_sub7 = stateMPC_s + 309;
stateMPC_FLOAT* stateMPC_lubbysub7 = stateMPC_lbys + 309;
stateMPC_FLOAT stateMPC_riub7[11];
stateMPC_FLOAT* stateMPC_dlubaff7 = stateMPC_dl_aff + 309;
stateMPC_FLOAT* stateMPC_dsubaff7 = stateMPC_ds_aff + 309;
stateMPC_FLOAT* stateMPC_dlubcc7 = stateMPC_dl_cc + 309;
stateMPC_FLOAT* stateMPC_dsubcc7 = stateMPC_ds_cc + 309;
stateMPC_FLOAT* stateMPC_ccrhsub7 = stateMPC_ccrhs + 309;
stateMPC_FLOAT stateMPC_Phi7[29];
stateMPC_FLOAT stateMPC_W7[29];
stateMPC_FLOAT stateMPC_Ysd7[81];
stateMPC_FLOAT stateMPC_Lsd7[81];
stateMPC_FLOAT* stateMPC_z8 = stateMPC_z + 232;
stateMPC_FLOAT* stateMPC_dzaff8 = stateMPC_dz_aff + 232;
stateMPC_FLOAT* stateMPC_dzcc8 = stateMPC_dz_cc + 232;
stateMPC_FLOAT* stateMPC_rd8 = stateMPC_rd + 232;
stateMPC_FLOAT stateMPC_Lbyrd8[29];
stateMPC_FLOAT* stateMPC_grad_cost8 = stateMPC_grad_cost + 232;
stateMPC_FLOAT* stateMPC_grad_eq8 = stateMPC_grad_eq + 232;
stateMPC_FLOAT* stateMPC_grad_ineq8 = stateMPC_grad_ineq + 232;
stateMPC_FLOAT stateMPC_ctv8[29];
stateMPC_FLOAT* stateMPC_v8 = stateMPC_v + 72;
stateMPC_FLOAT stateMPC_re8[9];
stateMPC_FLOAT stateMPC_beta8[9];
stateMPC_FLOAT stateMPC_betacc8[9];
stateMPC_FLOAT* stateMPC_dvaff8 = stateMPC_dv_aff + 72;
stateMPC_FLOAT* stateMPC_dvcc8 = stateMPC_dv_cc + 72;
stateMPC_FLOAT stateMPC_V8[261];
stateMPC_FLOAT stateMPC_Yd8[45];
stateMPC_FLOAT stateMPC_Ld8[45];
stateMPC_FLOAT stateMPC_yy8[9];
stateMPC_FLOAT stateMPC_bmy8[9];
int stateMPC_lbIdx8[29] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28};
stateMPC_FLOAT* stateMPC_llb8 = stateMPC_l + 320;
stateMPC_FLOAT* stateMPC_slb8 = stateMPC_s + 320;
stateMPC_FLOAT* stateMPC_llbbyslb8 = stateMPC_lbys + 320;
stateMPC_FLOAT stateMPC_rilb8[29];
stateMPC_FLOAT* stateMPC_dllbaff8 = stateMPC_dl_aff + 320;
stateMPC_FLOAT* stateMPC_dslbaff8 = stateMPC_ds_aff + 320;
stateMPC_FLOAT* stateMPC_dllbcc8 = stateMPC_dl_cc + 320;
stateMPC_FLOAT* stateMPC_dslbcc8 = stateMPC_ds_cc + 320;
stateMPC_FLOAT* stateMPC_ccrhsl8 = stateMPC_ccrhs + 320;
int stateMPC_ubIdx8[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
stateMPC_FLOAT* stateMPC_lub8 = stateMPC_l + 349;
stateMPC_FLOAT* stateMPC_sub8 = stateMPC_s + 349;
stateMPC_FLOAT* stateMPC_lubbysub8 = stateMPC_lbys + 349;
stateMPC_FLOAT stateMPC_riub8[11];
stateMPC_FLOAT* stateMPC_dlubaff8 = stateMPC_dl_aff + 349;
stateMPC_FLOAT* stateMPC_dsubaff8 = stateMPC_ds_aff + 349;
stateMPC_FLOAT* stateMPC_dlubcc8 = stateMPC_dl_cc + 349;
stateMPC_FLOAT* stateMPC_dsubcc8 = stateMPC_ds_cc + 349;
stateMPC_FLOAT* stateMPC_ccrhsub8 = stateMPC_ccrhs + 349;
stateMPC_FLOAT stateMPC_Phi8[29];
stateMPC_FLOAT stateMPC_W8[29];
stateMPC_FLOAT stateMPC_Ysd8[81];
stateMPC_FLOAT stateMPC_Lsd8[81];
stateMPC_FLOAT* stateMPC_z9 = stateMPC_z + 261;
stateMPC_FLOAT* stateMPC_dzaff9 = stateMPC_dz_aff + 261;
stateMPC_FLOAT* stateMPC_dzcc9 = stateMPC_dz_cc + 261;
stateMPC_FLOAT* stateMPC_rd9 = stateMPC_rd + 261;
stateMPC_FLOAT stateMPC_Lbyrd9[9];
stateMPC_FLOAT* stateMPC_grad_cost9 = stateMPC_grad_cost + 261;
stateMPC_FLOAT* stateMPC_grad_eq9 = stateMPC_grad_eq + 261;
stateMPC_FLOAT* stateMPC_grad_ineq9 = stateMPC_grad_ineq + 261;
stateMPC_FLOAT stateMPC_ctv9[9];
stateMPC_FLOAT* stateMPC_v9 = stateMPC_v + 81;
stateMPC_FLOAT stateMPC_re9[9];
stateMPC_FLOAT stateMPC_beta9[9];
stateMPC_FLOAT stateMPC_betacc9[9];
stateMPC_FLOAT* stateMPC_dvaff9 = stateMPC_dv_aff + 81;
stateMPC_FLOAT* stateMPC_dvcc9 = stateMPC_dv_cc + 81;
stateMPC_FLOAT stateMPC_V9[81];
stateMPC_FLOAT stateMPC_Yd9[45];
stateMPC_FLOAT stateMPC_Ld9[45];
stateMPC_FLOAT stateMPC_yy9[9];
stateMPC_FLOAT stateMPC_bmy9[9];
int stateMPC_lbIdx9[9] = {0, 1, 2, 3, 4, 5, 6, 7, 8};
stateMPC_FLOAT* stateMPC_llb9 = stateMPC_l + 360;
stateMPC_FLOAT* stateMPC_slb9 = stateMPC_s + 360;
stateMPC_FLOAT* stateMPC_llbbyslb9 = stateMPC_lbys + 360;
stateMPC_FLOAT stateMPC_rilb9[9];
stateMPC_FLOAT* stateMPC_dllbaff9 = stateMPC_dl_aff + 360;
stateMPC_FLOAT* stateMPC_dslbaff9 = stateMPC_ds_aff + 360;
stateMPC_FLOAT* stateMPC_dllbcc9 = stateMPC_dl_cc + 360;
stateMPC_FLOAT* stateMPC_dslbcc9 = stateMPC_ds_cc + 360;
stateMPC_FLOAT* stateMPC_ccrhsl9 = stateMPC_ccrhs + 360;
int stateMPC_ubIdx9[9] = {0, 1, 2, 3, 4, 5, 6, 7, 8};
stateMPC_FLOAT* stateMPC_lub9 = stateMPC_l + 369;
stateMPC_FLOAT* stateMPC_sub9 = stateMPC_s + 369;
stateMPC_FLOAT* stateMPC_lubbysub9 = stateMPC_lbys + 369;
stateMPC_FLOAT stateMPC_riub9[9];
stateMPC_FLOAT* stateMPC_dlubaff9 = stateMPC_dl_aff + 369;
stateMPC_FLOAT* stateMPC_dsubaff9 = stateMPC_ds_aff + 369;
stateMPC_FLOAT* stateMPC_dlubcc9 = stateMPC_dl_cc + 369;
stateMPC_FLOAT* stateMPC_dsubcc9 = stateMPC_ds_cc + 369;
stateMPC_FLOAT* stateMPC_ccrhsub9 = stateMPC_ccrhs + 369;
stateMPC_FLOAT stateMPC_Phi9[9];
stateMPC_FLOAT stateMPC_D9[9] = {-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000};
stateMPC_FLOAT stateMPC_W9[9];
stateMPC_FLOAT stateMPC_Ysd9[81];
stateMPC_FLOAT stateMPC_Lsd9[81];
stateMPC_FLOAT musigma;
stateMPC_FLOAT sigma_3rdroot;
stateMPC_FLOAT stateMPC_Diag1_0[29];
stateMPC_FLOAT stateMPC_Diag2_0[29];
stateMPC_FLOAT stateMPC_L_0[406];




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
stateMPC_LA_INITIALIZEVECTOR_270(stateMPC_z, 0);
stateMPC_LA_INITIALIZEVECTOR_90(stateMPC_v, 1);
stateMPC_LA_INITIALIZEVECTOR_378(stateMPC_l, 1);
stateMPC_LA_INITIALIZEVECTOR_378(stateMPC_s, 1);
info->mu = 0;
stateMPC_LA_DOTACC_378(stateMPC_l, stateMPC_s, &info->mu);
info->mu /= 378;
while( 1 ){
info->pobj = 0;
stateMPC_LA_DIAG_QUADFCN_29(params->H1, params->f1, stateMPC_z0, stateMPC_grad_cost0, &info->pobj);
stateMPC_LA_DIAG_QUADFCN_29(params->H2, params->f2, stateMPC_z1, stateMPC_grad_cost1, &info->pobj);
stateMPC_LA_DIAG_QUADFCN_29(params->H3, params->f3, stateMPC_z2, stateMPC_grad_cost2, &info->pobj);
stateMPC_LA_DIAG_QUADFCN_29(params->H4, params->f4, stateMPC_z3, stateMPC_grad_cost3, &info->pobj);
stateMPC_LA_DIAG_QUADFCN_29(params->H5, params->f5, stateMPC_z4, stateMPC_grad_cost4, &info->pobj);
stateMPC_LA_DIAG_QUADFCN_29(params->H6, params->f6, stateMPC_z5, stateMPC_grad_cost5, &info->pobj);
stateMPC_LA_DIAG_QUADFCN_29(params->H7, params->f7, stateMPC_z6, stateMPC_grad_cost6, &info->pobj);
stateMPC_LA_DIAG_QUADFCN_29(params->H8, params->f8, stateMPC_z7, stateMPC_grad_cost7, &info->pobj);
stateMPC_LA_DIAG_QUADFCN_29(params->H9, params->f9, stateMPC_z8, stateMPC_grad_cost8, &info->pobj);
stateMPC_LA_DIAG_QUADFCN_9(params->H10, params->f10, stateMPC_z9, stateMPC_grad_cost9, &info->pobj);
info->res_eq = 0;
info->dgap = 0;
stateMPC_LA_DIAGZERO_MVMSUB6_9(stateMPC_D0, stateMPC_z0, params->e1, stateMPC_v0, stateMPC_re0, &info->dgap, &info->res_eq);
stateMPC_LA_DENSE_DIAGZERO_MVMSUB3_9_29_29(params->C1, stateMPC_z0, stateMPC_D1, stateMPC_z1, params->e2, stateMPC_v1, stateMPC_re1, &info->dgap, &info->res_eq);
stateMPC_LA_DENSE_DIAGZERO_MVMSUB3_9_29_29(params->C2, stateMPC_z1, stateMPC_D1, stateMPC_z2, params->e3, stateMPC_v2, stateMPC_re2, &info->dgap, &info->res_eq);
stateMPC_LA_DENSE_DIAGZERO_MVMSUB3_9_29_29(params->C3, stateMPC_z2, stateMPC_D1, stateMPC_z3, params->e4, stateMPC_v3, stateMPC_re3, &info->dgap, &info->res_eq);
stateMPC_LA_DENSE_DIAGZERO_MVMSUB3_9_29_29(params->C4, stateMPC_z3, stateMPC_D1, stateMPC_z4, params->e5, stateMPC_v4, stateMPC_re4, &info->dgap, &info->res_eq);
stateMPC_LA_DENSE_DIAGZERO_MVMSUB3_9_29_29(params->C5, stateMPC_z4, stateMPC_D1, stateMPC_z5, params->e6, stateMPC_v5, stateMPC_re5, &info->dgap, &info->res_eq);
stateMPC_LA_DENSE_DIAGZERO_MVMSUB3_9_29_29(params->C6, stateMPC_z5, stateMPC_D1, stateMPC_z6, params->e7, stateMPC_v6, stateMPC_re6, &info->dgap, &info->res_eq);
stateMPC_LA_DENSE_DIAGZERO_MVMSUB3_9_29_29(params->C7, stateMPC_z6, stateMPC_D1, stateMPC_z7, params->e8, stateMPC_v7, stateMPC_re7, &info->dgap, &info->res_eq);
stateMPC_LA_DENSE_DIAGZERO_MVMSUB3_9_29_29(params->C8, stateMPC_z7, stateMPC_D1, stateMPC_z8, params->e9, stateMPC_v8, stateMPC_re8, &info->dgap, &info->res_eq);
stateMPC_LA_DENSE_DIAGZERO_MVMSUB3_9_29_9(params->C9, stateMPC_z8, stateMPC_D9, stateMPC_z9, params->e10, stateMPC_v9, stateMPC_re9, &info->dgap, &info->res_eq);
stateMPC_LA_DENSE_DIAGZERO_MTVM2_9_29_9(params->C1, stateMPC_v1, stateMPC_D0, stateMPC_v0, stateMPC_grad_eq0);
stateMPC_LA_DENSE_DIAGZERO_MTVM2_9_29_9(params->C2, stateMPC_v2, stateMPC_D1, stateMPC_v1, stateMPC_grad_eq1);
stateMPC_LA_DENSE_DIAGZERO_MTVM2_9_29_9(params->C3, stateMPC_v3, stateMPC_D1, stateMPC_v2, stateMPC_grad_eq2);
stateMPC_LA_DENSE_DIAGZERO_MTVM2_9_29_9(params->C4, stateMPC_v4, stateMPC_D1, stateMPC_v3, stateMPC_grad_eq3);
stateMPC_LA_DENSE_DIAGZERO_MTVM2_9_29_9(params->C5, stateMPC_v5, stateMPC_D1, stateMPC_v4, stateMPC_grad_eq4);
stateMPC_LA_DENSE_DIAGZERO_MTVM2_9_29_9(params->C6, stateMPC_v6, stateMPC_D1, stateMPC_v5, stateMPC_grad_eq5);
stateMPC_LA_DENSE_DIAGZERO_MTVM2_9_29_9(params->C7, stateMPC_v7, stateMPC_D1, stateMPC_v6, stateMPC_grad_eq6);
stateMPC_LA_DENSE_DIAGZERO_MTVM2_9_29_9(params->C8, stateMPC_v8, stateMPC_D1, stateMPC_v7, stateMPC_grad_eq7);
stateMPC_LA_DENSE_DIAGZERO_MTVM2_9_29_9(params->C9, stateMPC_v9, stateMPC_D1, stateMPC_v8, stateMPC_grad_eq8);
stateMPC_LA_DIAGZERO_MTVM_9_9(stateMPC_D9, stateMPC_v9, stateMPC_grad_eq9);
info->res_ineq = 0;
stateMPC_LA_VSUBADD3_29(params->lb1, stateMPC_z0, stateMPC_lbIdx0, stateMPC_llb0, stateMPC_slb0, stateMPC_rilb0, &info->dgap, &info->res_ineq);
stateMPC_LA_VSUBADD2_11(stateMPC_z0, stateMPC_ubIdx0, params->ub1, stateMPC_lub0, stateMPC_sub0, stateMPC_riub0, &info->dgap, &info->res_ineq);
stateMPC_LA_VSUBADD3_29(params->lb2, stateMPC_z1, stateMPC_lbIdx1, stateMPC_llb1, stateMPC_slb1, stateMPC_rilb1, &info->dgap, &info->res_ineq);
stateMPC_LA_VSUBADD2_11(stateMPC_z1, stateMPC_ubIdx1, params->ub2, stateMPC_lub1, stateMPC_sub1, stateMPC_riub1, &info->dgap, &info->res_ineq);
stateMPC_LA_VSUBADD3_29(params->lb3, stateMPC_z2, stateMPC_lbIdx2, stateMPC_llb2, stateMPC_slb2, stateMPC_rilb2, &info->dgap, &info->res_ineq);
stateMPC_LA_VSUBADD2_11(stateMPC_z2, stateMPC_ubIdx2, params->ub3, stateMPC_lub2, stateMPC_sub2, stateMPC_riub2, &info->dgap, &info->res_ineq);
stateMPC_LA_VSUBADD3_29(params->lb4, stateMPC_z3, stateMPC_lbIdx3, stateMPC_llb3, stateMPC_slb3, stateMPC_rilb3, &info->dgap, &info->res_ineq);
stateMPC_LA_VSUBADD2_11(stateMPC_z3, stateMPC_ubIdx3, params->ub4, stateMPC_lub3, stateMPC_sub3, stateMPC_riub3, &info->dgap, &info->res_ineq);
stateMPC_LA_VSUBADD3_29(params->lb5, stateMPC_z4, stateMPC_lbIdx4, stateMPC_llb4, stateMPC_slb4, stateMPC_rilb4, &info->dgap, &info->res_ineq);
stateMPC_LA_VSUBADD2_11(stateMPC_z4, stateMPC_ubIdx4, params->ub5, stateMPC_lub4, stateMPC_sub4, stateMPC_riub4, &info->dgap, &info->res_ineq);
stateMPC_LA_VSUBADD3_29(params->lb6, stateMPC_z5, stateMPC_lbIdx5, stateMPC_llb5, stateMPC_slb5, stateMPC_rilb5, &info->dgap, &info->res_ineq);
stateMPC_LA_VSUBADD2_11(stateMPC_z5, stateMPC_ubIdx5, params->ub6, stateMPC_lub5, stateMPC_sub5, stateMPC_riub5, &info->dgap, &info->res_ineq);
stateMPC_LA_VSUBADD3_29(params->lb7, stateMPC_z6, stateMPC_lbIdx6, stateMPC_llb6, stateMPC_slb6, stateMPC_rilb6, &info->dgap, &info->res_ineq);
stateMPC_LA_VSUBADD2_11(stateMPC_z6, stateMPC_ubIdx6, params->ub7, stateMPC_lub6, stateMPC_sub6, stateMPC_riub6, &info->dgap, &info->res_ineq);
stateMPC_LA_VSUBADD3_29(params->lb8, stateMPC_z7, stateMPC_lbIdx7, stateMPC_llb7, stateMPC_slb7, stateMPC_rilb7, &info->dgap, &info->res_ineq);
stateMPC_LA_VSUBADD2_11(stateMPC_z7, stateMPC_ubIdx7, params->ub8, stateMPC_lub7, stateMPC_sub7, stateMPC_riub7, &info->dgap, &info->res_ineq);
stateMPC_LA_VSUBADD3_29(params->lb9, stateMPC_z8, stateMPC_lbIdx8, stateMPC_llb8, stateMPC_slb8, stateMPC_rilb8, &info->dgap, &info->res_ineq);
stateMPC_LA_VSUBADD2_11(stateMPC_z8, stateMPC_ubIdx8, params->ub9, stateMPC_lub8, stateMPC_sub8, stateMPC_riub8, &info->dgap, &info->res_ineq);
stateMPC_LA_VSUBADD3_9(params->lb10, stateMPC_z9, stateMPC_lbIdx9, stateMPC_llb9, stateMPC_slb9, stateMPC_rilb9, &info->dgap, &info->res_ineq);
stateMPC_LA_VSUBADD2_9(stateMPC_z9, stateMPC_ubIdx9, params->ub10, stateMPC_lub9, stateMPC_sub9, stateMPC_riub9, &info->dgap, &info->res_ineq);
stateMPC_LA_INEQ_B_GRAD_29_29_11(stateMPC_lub0, stateMPC_sub0, stateMPC_riub0, stateMPC_llb0, stateMPC_slb0, stateMPC_rilb0, stateMPC_lbIdx0, stateMPC_ubIdx0, stateMPC_grad_ineq0, stateMPC_lubbysub0, stateMPC_llbbyslb0);
stateMPC_LA_INEQ_B_GRAD_29_29_11(stateMPC_lub1, stateMPC_sub1, stateMPC_riub1, stateMPC_llb1, stateMPC_slb1, stateMPC_rilb1, stateMPC_lbIdx1, stateMPC_ubIdx1, stateMPC_grad_ineq1, stateMPC_lubbysub1, stateMPC_llbbyslb1);
stateMPC_LA_INEQ_B_GRAD_29_29_11(stateMPC_lub2, stateMPC_sub2, stateMPC_riub2, stateMPC_llb2, stateMPC_slb2, stateMPC_rilb2, stateMPC_lbIdx2, stateMPC_ubIdx2, stateMPC_grad_ineq2, stateMPC_lubbysub2, stateMPC_llbbyslb2);
stateMPC_LA_INEQ_B_GRAD_29_29_11(stateMPC_lub3, stateMPC_sub3, stateMPC_riub3, stateMPC_llb3, stateMPC_slb3, stateMPC_rilb3, stateMPC_lbIdx3, stateMPC_ubIdx3, stateMPC_grad_ineq3, stateMPC_lubbysub3, stateMPC_llbbyslb3);
stateMPC_LA_INEQ_B_GRAD_29_29_11(stateMPC_lub4, stateMPC_sub4, stateMPC_riub4, stateMPC_llb4, stateMPC_slb4, stateMPC_rilb4, stateMPC_lbIdx4, stateMPC_ubIdx4, stateMPC_grad_ineq4, stateMPC_lubbysub4, stateMPC_llbbyslb4);
stateMPC_LA_INEQ_B_GRAD_29_29_11(stateMPC_lub5, stateMPC_sub5, stateMPC_riub5, stateMPC_llb5, stateMPC_slb5, stateMPC_rilb5, stateMPC_lbIdx5, stateMPC_ubIdx5, stateMPC_grad_ineq5, stateMPC_lubbysub5, stateMPC_llbbyslb5);
stateMPC_LA_INEQ_B_GRAD_29_29_11(stateMPC_lub6, stateMPC_sub6, stateMPC_riub6, stateMPC_llb6, stateMPC_slb6, stateMPC_rilb6, stateMPC_lbIdx6, stateMPC_ubIdx6, stateMPC_grad_ineq6, stateMPC_lubbysub6, stateMPC_llbbyslb6);
stateMPC_LA_INEQ_B_GRAD_29_29_11(stateMPC_lub7, stateMPC_sub7, stateMPC_riub7, stateMPC_llb7, stateMPC_slb7, stateMPC_rilb7, stateMPC_lbIdx7, stateMPC_ubIdx7, stateMPC_grad_ineq7, stateMPC_lubbysub7, stateMPC_llbbyslb7);
stateMPC_LA_INEQ_B_GRAD_29_29_11(stateMPC_lub8, stateMPC_sub8, stateMPC_riub8, stateMPC_llb8, stateMPC_slb8, stateMPC_rilb8, stateMPC_lbIdx8, stateMPC_ubIdx8, stateMPC_grad_ineq8, stateMPC_lubbysub8, stateMPC_llbbyslb8);
stateMPC_LA_INEQ_B_GRAD_9_9_9(stateMPC_lub9, stateMPC_sub9, stateMPC_riub9, stateMPC_llb9, stateMPC_slb9, stateMPC_rilb9, stateMPC_lbIdx9, stateMPC_ubIdx9, stateMPC_grad_ineq9, stateMPC_lubbysub9, stateMPC_llbbyslb9);
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
stateMPC_LA_VVADD3_270(stateMPC_grad_cost, stateMPC_grad_eq, stateMPC_grad_ineq, stateMPC_rd);
stateMPC_LA_DIAG_CHOL_LBUB_29_29_11(params->H1, stateMPC_llbbyslb0, stateMPC_lbIdx0, stateMPC_lubbysub0, stateMPC_ubIdx0, stateMPC_Phi0);
stateMPC_LA_DIAG_MATRIXFORWARDSUB_9_29(stateMPC_Phi0, params->C1, stateMPC_V0);
stateMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_9_29(stateMPC_Phi0, stateMPC_D0, stateMPC_W0);
stateMPC_LA_DENSE_DIAGZERO_MMTM_9_29_9(stateMPC_W0, stateMPC_V0, stateMPC_Ysd1);
stateMPC_LA_DIAG_FORWARDSUB_29(stateMPC_Phi0, stateMPC_rd0, stateMPC_Lbyrd0);
stateMPC_LA_DIAG_CHOL_LBUB_29_29_11(params->H2, stateMPC_llbbyslb1, stateMPC_lbIdx1, stateMPC_lubbysub1, stateMPC_ubIdx1, stateMPC_Phi1);
stateMPC_LA_DIAG_MATRIXFORWARDSUB_9_29(stateMPC_Phi1, params->C2, stateMPC_V1);
stateMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_9_29(stateMPC_Phi1, stateMPC_D1, stateMPC_W1);
stateMPC_LA_DENSE_DIAGZERO_MMTM_9_29_9(stateMPC_W1, stateMPC_V1, stateMPC_Ysd2);
stateMPC_LA_DIAG_FORWARDSUB_29(stateMPC_Phi1, stateMPC_rd1, stateMPC_Lbyrd1);
stateMPC_LA_DIAG_CHOL_LBUB_29_29_11(params->H3, stateMPC_llbbyslb2, stateMPC_lbIdx2, stateMPC_lubbysub2, stateMPC_ubIdx2, stateMPC_Phi2);
stateMPC_LA_DIAG_MATRIXFORWARDSUB_9_29(stateMPC_Phi2, params->C3, stateMPC_V2);
stateMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_9_29(stateMPC_Phi2, stateMPC_D1, stateMPC_W2);
stateMPC_LA_DENSE_DIAGZERO_MMTM_9_29_9(stateMPC_W2, stateMPC_V2, stateMPC_Ysd3);
stateMPC_LA_DIAG_FORWARDSUB_29(stateMPC_Phi2, stateMPC_rd2, stateMPC_Lbyrd2);
stateMPC_LA_DIAG_CHOL_LBUB_29_29_11(params->H4, stateMPC_llbbyslb3, stateMPC_lbIdx3, stateMPC_lubbysub3, stateMPC_ubIdx3, stateMPC_Phi3);
stateMPC_LA_DIAG_MATRIXFORWARDSUB_9_29(stateMPC_Phi3, params->C4, stateMPC_V3);
stateMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_9_29(stateMPC_Phi3, stateMPC_D1, stateMPC_W3);
stateMPC_LA_DENSE_DIAGZERO_MMTM_9_29_9(stateMPC_W3, stateMPC_V3, stateMPC_Ysd4);
stateMPC_LA_DIAG_FORWARDSUB_29(stateMPC_Phi3, stateMPC_rd3, stateMPC_Lbyrd3);
stateMPC_LA_DIAG_CHOL_LBUB_29_29_11(params->H5, stateMPC_llbbyslb4, stateMPC_lbIdx4, stateMPC_lubbysub4, stateMPC_ubIdx4, stateMPC_Phi4);
stateMPC_LA_DIAG_MATRIXFORWARDSUB_9_29(stateMPC_Phi4, params->C5, stateMPC_V4);
stateMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_9_29(stateMPC_Phi4, stateMPC_D1, stateMPC_W4);
stateMPC_LA_DENSE_DIAGZERO_MMTM_9_29_9(stateMPC_W4, stateMPC_V4, stateMPC_Ysd5);
stateMPC_LA_DIAG_FORWARDSUB_29(stateMPC_Phi4, stateMPC_rd4, stateMPC_Lbyrd4);
stateMPC_LA_DIAG_CHOL_LBUB_29_29_11(params->H6, stateMPC_llbbyslb5, stateMPC_lbIdx5, stateMPC_lubbysub5, stateMPC_ubIdx5, stateMPC_Phi5);
stateMPC_LA_DIAG_MATRIXFORWARDSUB_9_29(stateMPC_Phi5, params->C6, stateMPC_V5);
stateMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_9_29(stateMPC_Phi5, stateMPC_D1, stateMPC_W5);
stateMPC_LA_DENSE_DIAGZERO_MMTM_9_29_9(stateMPC_W5, stateMPC_V5, stateMPC_Ysd6);
stateMPC_LA_DIAG_FORWARDSUB_29(stateMPC_Phi5, stateMPC_rd5, stateMPC_Lbyrd5);
stateMPC_LA_DIAG_CHOL_LBUB_29_29_11(params->H7, stateMPC_llbbyslb6, stateMPC_lbIdx6, stateMPC_lubbysub6, stateMPC_ubIdx6, stateMPC_Phi6);
stateMPC_LA_DIAG_MATRIXFORWARDSUB_9_29(stateMPC_Phi6, params->C7, stateMPC_V6);
stateMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_9_29(stateMPC_Phi6, stateMPC_D1, stateMPC_W6);
stateMPC_LA_DENSE_DIAGZERO_MMTM_9_29_9(stateMPC_W6, stateMPC_V6, stateMPC_Ysd7);
stateMPC_LA_DIAG_FORWARDSUB_29(stateMPC_Phi6, stateMPC_rd6, stateMPC_Lbyrd6);
stateMPC_LA_DIAG_CHOL_LBUB_29_29_11(params->H8, stateMPC_llbbyslb7, stateMPC_lbIdx7, stateMPC_lubbysub7, stateMPC_ubIdx7, stateMPC_Phi7);
stateMPC_LA_DIAG_MATRIXFORWARDSUB_9_29(stateMPC_Phi7, params->C8, stateMPC_V7);
stateMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_9_29(stateMPC_Phi7, stateMPC_D1, stateMPC_W7);
stateMPC_LA_DENSE_DIAGZERO_MMTM_9_29_9(stateMPC_W7, stateMPC_V7, stateMPC_Ysd8);
stateMPC_LA_DIAG_FORWARDSUB_29(stateMPC_Phi7, stateMPC_rd7, stateMPC_Lbyrd7);
stateMPC_LA_DIAG_CHOL_LBUB_29_29_11(params->H9, stateMPC_llbbyslb8, stateMPC_lbIdx8, stateMPC_lubbysub8, stateMPC_ubIdx8, stateMPC_Phi8);
stateMPC_LA_DIAG_MATRIXFORWARDSUB_9_29(stateMPC_Phi8, params->C9, stateMPC_V8);
stateMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_9_29(stateMPC_Phi8, stateMPC_D1, stateMPC_W8);
stateMPC_LA_DENSE_DIAGZERO_MMTM_9_29_9(stateMPC_W8, stateMPC_V8, stateMPC_Ysd9);
stateMPC_LA_DIAG_FORWARDSUB_29(stateMPC_Phi8, stateMPC_rd8, stateMPC_Lbyrd8);
stateMPC_LA_DIAG_CHOL_ONELOOP_LBUB_9_9_9(params->H10, stateMPC_llbbyslb9, stateMPC_lbIdx9, stateMPC_lubbysub9, stateMPC_ubIdx9, stateMPC_Phi9);
stateMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_9_9(stateMPC_Phi9, stateMPC_D9, stateMPC_W9);
stateMPC_LA_DIAG_FORWARDSUB_9(stateMPC_Phi9, stateMPC_rd9, stateMPC_Lbyrd9);
stateMPC_LA_DIAGZERO_MMT_9(stateMPC_W0, stateMPC_Yd0);
stateMPC_LA_DIAGZERO_MVMSUB7_9(stateMPC_W0, stateMPC_Lbyrd0, stateMPC_re0, stateMPC_beta0);
stateMPC_LA_DENSE_DIAGZERO_MMT2_9_29_29(stateMPC_V0, stateMPC_W1, stateMPC_Yd1);
stateMPC_LA_DENSE_DIAGZERO_2MVMSUB2_9_29_29(stateMPC_V0, stateMPC_Lbyrd0, stateMPC_W1, stateMPC_Lbyrd1, stateMPC_re1, stateMPC_beta1);
stateMPC_LA_DENSE_DIAGZERO_MMT2_9_29_29(stateMPC_V1, stateMPC_W2, stateMPC_Yd2);
stateMPC_LA_DENSE_DIAGZERO_2MVMSUB2_9_29_29(stateMPC_V1, stateMPC_Lbyrd1, stateMPC_W2, stateMPC_Lbyrd2, stateMPC_re2, stateMPC_beta2);
stateMPC_LA_DENSE_DIAGZERO_MMT2_9_29_29(stateMPC_V2, stateMPC_W3, stateMPC_Yd3);
stateMPC_LA_DENSE_DIAGZERO_2MVMSUB2_9_29_29(stateMPC_V2, stateMPC_Lbyrd2, stateMPC_W3, stateMPC_Lbyrd3, stateMPC_re3, stateMPC_beta3);
stateMPC_LA_DENSE_DIAGZERO_MMT2_9_29_29(stateMPC_V3, stateMPC_W4, stateMPC_Yd4);
stateMPC_LA_DENSE_DIAGZERO_2MVMSUB2_9_29_29(stateMPC_V3, stateMPC_Lbyrd3, stateMPC_W4, stateMPC_Lbyrd4, stateMPC_re4, stateMPC_beta4);
stateMPC_LA_DENSE_DIAGZERO_MMT2_9_29_29(stateMPC_V4, stateMPC_W5, stateMPC_Yd5);
stateMPC_LA_DENSE_DIAGZERO_2MVMSUB2_9_29_29(stateMPC_V4, stateMPC_Lbyrd4, stateMPC_W5, stateMPC_Lbyrd5, stateMPC_re5, stateMPC_beta5);
stateMPC_LA_DENSE_DIAGZERO_MMT2_9_29_29(stateMPC_V5, stateMPC_W6, stateMPC_Yd6);
stateMPC_LA_DENSE_DIAGZERO_2MVMSUB2_9_29_29(stateMPC_V5, stateMPC_Lbyrd5, stateMPC_W6, stateMPC_Lbyrd6, stateMPC_re6, stateMPC_beta6);
stateMPC_LA_DENSE_DIAGZERO_MMT2_9_29_29(stateMPC_V6, stateMPC_W7, stateMPC_Yd7);
stateMPC_LA_DENSE_DIAGZERO_2MVMSUB2_9_29_29(stateMPC_V6, stateMPC_Lbyrd6, stateMPC_W7, stateMPC_Lbyrd7, stateMPC_re7, stateMPC_beta7);
stateMPC_LA_DENSE_DIAGZERO_MMT2_9_29_29(stateMPC_V7, stateMPC_W8, stateMPC_Yd8);
stateMPC_LA_DENSE_DIAGZERO_2MVMSUB2_9_29_29(stateMPC_V7, stateMPC_Lbyrd7, stateMPC_W8, stateMPC_Lbyrd8, stateMPC_re8, stateMPC_beta8);
stateMPC_LA_DENSE_DIAGZERO_MMT2_9_29_9(stateMPC_V8, stateMPC_W9, stateMPC_Yd9);
stateMPC_LA_DENSE_DIAGZERO_2MVMSUB2_9_29_9(stateMPC_V8, stateMPC_Lbyrd8, stateMPC_W9, stateMPC_Lbyrd9, stateMPC_re9, stateMPC_beta9);
stateMPC_LA_DENSE_CHOL_9(stateMPC_Yd0, stateMPC_Ld0);
stateMPC_LA_DENSE_FORWARDSUB_9(stateMPC_Ld0, stateMPC_beta0, stateMPC_yy0);
stateMPC_LA_DENSE_MATRIXTFORWARDSUB_9_9(stateMPC_Ld0, stateMPC_Ysd1, stateMPC_Lsd1);
stateMPC_LA_DENSE_MMTSUB_9_9(stateMPC_Lsd1, stateMPC_Yd1);
stateMPC_LA_DENSE_CHOL_9(stateMPC_Yd1, stateMPC_Ld1);
stateMPC_LA_DENSE_MVMSUB1_9_9(stateMPC_Lsd1, stateMPC_yy0, stateMPC_beta1, stateMPC_bmy1);
stateMPC_LA_DENSE_FORWARDSUB_9(stateMPC_Ld1, stateMPC_bmy1, stateMPC_yy1);
stateMPC_LA_DENSE_MATRIXTFORWARDSUB_9_9(stateMPC_Ld1, stateMPC_Ysd2, stateMPC_Lsd2);
stateMPC_LA_DENSE_MMTSUB_9_9(stateMPC_Lsd2, stateMPC_Yd2);
stateMPC_LA_DENSE_CHOL_9(stateMPC_Yd2, stateMPC_Ld2);
stateMPC_LA_DENSE_MVMSUB1_9_9(stateMPC_Lsd2, stateMPC_yy1, stateMPC_beta2, stateMPC_bmy2);
stateMPC_LA_DENSE_FORWARDSUB_9(stateMPC_Ld2, stateMPC_bmy2, stateMPC_yy2);
stateMPC_LA_DENSE_MATRIXTFORWARDSUB_9_9(stateMPC_Ld2, stateMPC_Ysd3, stateMPC_Lsd3);
stateMPC_LA_DENSE_MMTSUB_9_9(stateMPC_Lsd3, stateMPC_Yd3);
stateMPC_LA_DENSE_CHOL_9(stateMPC_Yd3, stateMPC_Ld3);
stateMPC_LA_DENSE_MVMSUB1_9_9(stateMPC_Lsd3, stateMPC_yy2, stateMPC_beta3, stateMPC_bmy3);
stateMPC_LA_DENSE_FORWARDSUB_9(stateMPC_Ld3, stateMPC_bmy3, stateMPC_yy3);
stateMPC_LA_DENSE_MATRIXTFORWARDSUB_9_9(stateMPC_Ld3, stateMPC_Ysd4, stateMPC_Lsd4);
stateMPC_LA_DENSE_MMTSUB_9_9(stateMPC_Lsd4, stateMPC_Yd4);
stateMPC_LA_DENSE_CHOL_9(stateMPC_Yd4, stateMPC_Ld4);
stateMPC_LA_DENSE_MVMSUB1_9_9(stateMPC_Lsd4, stateMPC_yy3, stateMPC_beta4, stateMPC_bmy4);
stateMPC_LA_DENSE_FORWARDSUB_9(stateMPC_Ld4, stateMPC_bmy4, stateMPC_yy4);
stateMPC_LA_DENSE_MATRIXTFORWARDSUB_9_9(stateMPC_Ld4, stateMPC_Ysd5, stateMPC_Lsd5);
stateMPC_LA_DENSE_MMTSUB_9_9(stateMPC_Lsd5, stateMPC_Yd5);
stateMPC_LA_DENSE_CHOL_9(stateMPC_Yd5, stateMPC_Ld5);
stateMPC_LA_DENSE_MVMSUB1_9_9(stateMPC_Lsd5, stateMPC_yy4, stateMPC_beta5, stateMPC_bmy5);
stateMPC_LA_DENSE_FORWARDSUB_9(stateMPC_Ld5, stateMPC_bmy5, stateMPC_yy5);
stateMPC_LA_DENSE_MATRIXTFORWARDSUB_9_9(stateMPC_Ld5, stateMPC_Ysd6, stateMPC_Lsd6);
stateMPC_LA_DENSE_MMTSUB_9_9(stateMPC_Lsd6, stateMPC_Yd6);
stateMPC_LA_DENSE_CHOL_9(stateMPC_Yd6, stateMPC_Ld6);
stateMPC_LA_DENSE_MVMSUB1_9_9(stateMPC_Lsd6, stateMPC_yy5, stateMPC_beta6, stateMPC_bmy6);
stateMPC_LA_DENSE_FORWARDSUB_9(stateMPC_Ld6, stateMPC_bmy6, stateMPC_yy6);
stateMPC_LA_DENSE_MATRIXTFORWARDSUB_9_9(stateMPC_Ld6, stateMPC_Ysd7, stateMPC_Lsd7);
stateMPC_LA_DENSE_MMTSUB_9_9(stateMPC_Lsd7, stateMPC_Yd7);
stateMPC_LA_DENSE_CHOL_9(stateMPC_Yd7, stateMPC_Ld7);
stateMPC_LA_DENSE_MVMSUB1_9_9(stateMPC_Lsd7, stateMPC_yy6, stateMPC_beta7, stateMPC_bmy7);
stateMPC_LA_DENSE_FORWARDSUB_9(stateMPC_Ld7, stateMPC_bmy7, stateMPC_yy7);
stateMPC_LA_DENSE_MATRIXTFORWARDSUB_9_9(stateMPC_Ld7, stateMPC_Ysd8, stateMPC_Lsd8);
stateMPC_LA_DENSE_MMTSUB_9_9(stateMPC_Lsd8, stateMPC_Yd8);
stateMPC_LA_DENSE_CHOL_9(stateMPC_Yd8, stateMPC_Ld8);
stateMPC_LA_DENSE_MVMSUB1_9_9(stateMPC_Lsd8, stateMPC_yy7, stateMPC_beta8, stateMPC_bmy8);
stateMPC_LA_DENSE_FORWARDSUB_9(stateMPC_Ld8, stateMPC_bmy8, stateMPC_yy8);
stateMPC_LA_DENSE_MATRIXTFORWARDSUB_9_9(stateMPC_Ld8, stateMPC_Ysd9, stateMPC_Lsd9);
stateMPC_LA_DENSE_MMTSUB_9_9(stateMPC_Lsd9, stateMPC_Yd9);
stateMPC_LA_DENSE_CHOL_9(stateMPC_Yd9, stateMPC_Ld9);
stateMPC_LA_DENSE_MVMSUB1_9_9(stateMPC_Lsd9, stateMPC_yy8, stateMPC_beta9, stateMPC_bmy9);
stateMPC_LA_DENSE_FORWARDSUB_9(stateMPC_Ld9, stateMPC_bmy9, stateMPC_yy9);
stateMPC_LA_DENSE_BACKWARDSUB_9(stateMPC_Ld9, stateMPC_yy9, stateMPC_dvaff9);
stateMPC_LA_DENSE_MTVMSUB_9_9(stateMPC_Lsd9, stateMPC_dvaff9, stateMPC_yy8, stateMPC_bmy8);
stateMPC_LA_DENSE_BACKWARDSUB_9(stateMPC_Ld8, stateMPC_bmy8, stateMPC_dvaff8);
stateMPC_LA_DENSE_MTVMSUB_9_9(stateMPC_Lsd8, stateMPC_dvaff8, stateMPC_yy7, stateMPC_bmy7);
stateMPC_LA_DENSE_BACKWARDSUB_9(stateMPC_Ld7, stateMPC_bmy7, stateMPC_dvaff7);
stateMPC_LA_DENSE_MTVMSUB_9_9(stateMPC_Lsd7, stateMPC_dvaff7, stateMPC_yy6, stateMPC_bmy6);
stateMPC_LA_DENSE_BACKWARDSUB_9(stateMPC_Ld6, stateMPC_bmy6, stateMPC_dvaff6);
stateMPC_LA_DENSE_MTVMSUB_9_9(stateMPC_Lsd6, stateMPC_dvaff6, stateMPC_yy5, stateMPC_bmy5);
stateMPC_LA_DENSE_BACKWARDSUB_9(stateMPC_Ld5, stateMPC_bmy5, stateMPC_dvaff5);
stateMPC_LA_DENSE_MTVMSUB_9_9(stateMPC_Lsd5, stateMPC_dvaff5, stateMPC_yy4, stateMPC_bmy4);
stateMPC_LA_DENSE_BACKWARDSUB_9(stateMPC_Ld4, stateMPC_bmy4, stateMPC_dvaff4);
stateMPC_LA_DENSE_MTVMSUB_9_9(stateMPC_Lsd4, stateMPC_dvaff4, stateMPC_yy3, stateMPC_bmy3);
stateMPC_LA_DENSE_BACKWARDSUB_9(stateMPC_Ld3, stateMPC_bmy3, stateMPC_dvaff3);
stateMPC_LA_DENSE_MTVMSUB_9_9(stateMPC_Lsd3, stateMPC_dvaff3, stateMPC_yy2, stateMPC_bmy2);
stateMPC_LA_DENSE_BACKWARDSUB_9(stateMPC_Ld2, stateMPC_bmy2, stateMPC_dvaff2);
stateMPC_LA_DENSE_MTVMSUB_9_9(stateMPC_Lsd2, stateMPC_dvaff2, stateMPC_yy1, stateMPC_bmy1);
stateMPC_LA_DENSE_BACKWARDSUB_9(stateMPC_Ld1, stateMPC_bmy1, stateMPC_dvaff1);
stateMPC_LA_DENSE_MTVMSUB_9_9(stateMPC_Lsd1, stateMPC_dvaff1, stateMPC_yy0, stateMPC_bmy0);
stateMPC_LA_DENSE_BACKWARDSUB_9(stateMPC_Ld0, stateMPC_bmy0, stateMPC_dvaff0);
stateMPC_LA_DENSE_DIAGZERO_MTVM2_9_29_9(params->C1, stateMPC_dvaff1, stateMPC_D0, stateMPC_dvaff0, stateMPC_grad_eq0);
stateMPC_LA_DENSE_DIAGZERO_MTVM2_9_29_9(params->C2, stateMPC_dvaff2, stateMPC_D1, stateMPC_dvaff1, stateMPC_grad_eq1);
stateMPC_LA_DENSE_DIAGZERO_MTVM2_9_29_9(params->C3, stateMPC_dvaff3, stateMPC_D1, stateMPC_dvaff2, stateMPC_grad_eq2);
stateMPC_LA_DENSE_DIAGZERO_MTVM2_9_29_9(params->C4, stateMPC_dvaff4, stateMPC_D1, stateMPC_dvaff3, stateMPC_grad_eq3);
stateMPC_LA_DENSE_DIAGZERO_MTVM2_9_29_9(params->C5, stateMPC_dvaff5, stateMPC_D1, stateMPC_dvaff4, stateMPC_grad_eq4);
stateMPC_LA_DENSE_DIAGZERO_MTVM2_9_29_9(params->C6, stateMPC_dvaff6, stateMPC_D1, stateMPC_dvaff5, stateMPC_grad_eq5);
stateMPC_LA_DENSE_DIAGZERO_MTVM2_9_29_9(params->C7, stateMPC_dvaff7, stateMPC_D1, stateMPC_dvaff6, stateMPC_grad_eq6);
stateMPC_LA_DENSE_DIAGZERO_MTVM2_9_29_9(params->C8, stateMPC_dvaff8, stateMPC_D1, stateMPC_dvaff7, stateMPC_grad_eq7);
stateMPC_LA_DENSE_DIAGZERO_MTVM2_9_29_9(params->C9, stateMPC_dvaff9, stateMPC_D1, stateMPC_dvaff8, stateMPC_grad_eq8);
stateMPC_LA_DIAGZERO_MTVM_9_9(stateMPC_D9, stateMPC_dvaff9, stateMPC_grad_eq9);
stateMPC_LA_VSUB2_270(stateMPC_rd, stateMPC_grad_eq, stateMPC_rd);
stateMPC_LA_DIAG_FORWARDBACKWARDSUB_29(stateMPC_Phi0, stateMPC_rd0, stateMPC_dzaff0);
stateMPC_LA_DIAG_FORWARDBACKWARDSUB_29(stateMPC_Phi1, stateMPC_rd1, stateMPC_dzaff1);
stateMPC_LA_DIAG_FORWARDBACKWARDSUB_29(stateMPC_Phi2, stateMPC_rd2, stateMPC_dzaff2);
stateMPC_LA_DIAG_FORWARDBACKWARDSUB_29(stateMPC_Phi3, stateMPC_rd3, stateMPC_dzaff3);
stateMPC_LA_DIAG_FORWARDBACKWARDSUB_29(stateMPC_Phi4, stateMPC_rd4, stateMPC_dzaff4);
stateMPC_LA_DIAG_FORWARDBACKWARDSUB_29(stateMPC_Phi5, stateMPC_rd5, stateMPC_dzaff5);
stateMPC_LA_DIAG_FORWARDBACKWARDSUB_29(stateMPC_Phi6, stateMPC_rd6, stateMPC_dzaff6);
stateMPC_LA_DIAG_FORWARDBACKWARDSUB_29(stateMPC_Phi7, stateMPC_rd7, stateMPC_dzaff7);
stateMPC_LA_DIAG_FORWARDBACKWARDSUB_29(stateMPC_Phi8, stateMPC_rd8, stateMPC_dzaff8);
stateMPC_LA_DIAG_FORWARDBACKWARDSUB_9(stateMPC_Phi9, stateMPC_rd9, stateMPC_dzaff9);
stateMPC_LA_VSUB_INDEXED_29(stateMPC_dzaff0, stateMPC_lbIdx0, stateMPC_rilb0, stateMPC_dslbaff0);
stateMPC_LA_VSUB3_29(stateMPC_llbbyslb0, stateMPC_dslbaff0, stateMPC_llb0, stateMPC_dllbaff0);
stateMPC_LA_VSUB2_INDEXED_11(stateMPC_riub0, stateMPC_dzaff0, stateMPC_ubIdx0, stateMPC_dsubaff0);
stateMPC_LA_VSUB3_11(stateMPC_lubbysub0, stateMPC_dsubaff0, stateMPC_lub0, stateMPC_dlubaff0);
stateMPC_LA_VSUB_INDEXED_29(stateMPC_dzaff1, stateMPC_lbIdx1, stateMPC_rilb1, stateMPC_dslbaff1);
stateMPC_LA_VSUB3_29(stateMPC_llbbyslb1, stateMPC_dslbaff1, stateMPC_llb1, stateMPC_dllbaff1);
stateMPC_LA_VSUB2_INDEXED_11(stateMPC_riub1, stateMPC_dzaff1, stateMPC_ubIdx1, stateMPC_dsubaff1);
stateMPC_LA_VSUB3_11(stateMPC_lubbysub1, stateMPC_dsubaff1, stateMPC_lub1, stateMPC_dlubaff1);
stateMPC_LA_VSUB_INDEXED_29(stateMPC_dzaff2, stateMPC_lbIdx2, stateMPC_rilb2, stateMPC_dslbaff2);
stateMPC_LA_VSUB3_29(stateMPC_llbbyslb2, stateMPC_dslbaff2, stateMPC_llb2, stateMPC_dllbaff2);
stateMPC_LA_VSUB2_INDEXED_11(stateMPC_riub2, stateMPC_dzaff2, stateMPC_ubIdx2, stateMPC_dsubaff2);
stateMPC_LA_VSUB3_11(stateMPC_lubbysub2, stateMPC_dsubaff2, stateMPC_lub2, stateMPC_dlubaff2);
stateMPC_LA_VSUB_INDEXED_29(stateMPC_dzaff3, stateMPC_lbIdx3, stateMPC_rilb3, stateMPC_dslbaff3);
stateMPC_LA_VSUB3_29(stateMPC_llbbyslb3, stateMPC_dslbaff3, stateMPC_llb3, stateMPC_dllbaff3);
stateMPC_LA_VSUB2_INDEXED_11(stateMPC_riub3, stateMPC_dzaff3, stateMPC_ubIdx3, stateMPC_dsubaff3);
stateMPC_LA_VSUB3_11(stateMPC_lubbysub3, stateMPC_dsubaff3, stateMPC_lub3, stateMPC_dlubaff3);
stateMPC_LA_VSUB_INDEXED_29(stateMPC_dzaff4, stateMPC_lbIdx4, stateMPC_rilb4, stateMPC_dslbaff4);
stateMPC_LA_VSUB3_29(stateMPC_llbbyslb4, stateMPC_dslbaff4, stateMPC_llb4, stateMPC_dllbaff4);
stateMPC_LA_VSUB2_INDEXED_11(stateMPC_riub4, stateMPC_dzaff4, stateMPC_ubIdx4, stateMPC_dsubaff4);
stateMPC_LA_VSUB3_11(stateMPC_lubbysub4, stateMPC_dsubaff4, stateMPC_lub4, stateMPC_dlubaff4);
stateMPC_LA_VSUB_INDEXED_29(stateMPC_dzaff5, stateMPC_lbIdx5, stateMPC_rilb5, stateMPC_dslbaff5);
stateMPC_LA_VSUB3_29(stateMPC_llbbyslb5, stateMPC_dslbaff5, stateMPC_llb5, stateMPC_dllbaff5);
stateMPC_LA_VSUB2_INDEXED_11(stateMPC_riub5, stateMPC_dzaff5, stateMPC_ubIdx5, stateMPC_dsubaff5);
stateMPC_LA_VSUB3_11(stateMPC_lubbysub5, stateMPC_dsubaff5, stateMPC_lub5, stateMPC_dlubaff5);
stateMPC_LA_VSUB_INDEXED_29(stateMPC_dzaff6, stateMPC_lbIdx6, stateMPC_rilb6, stateMPC_dslbaff6);
stateMPC_LA_VSUB3_29(stateMPC_llbbyslb6, stateMPC_dslbaff6, stateMPC_llb6, stateMPC_dllbaff6);
stateMPC_LA_VSUB2_INDEXED_11(stateMPC_riub6, stateMPC_dzaff6, stateMPC_ubIdx6, stateMPC_dsubaff6);
stateMPC_LA_VSUB3_11(stateMPC_lubbysub6, stateMPC_dsubaff6, stateMPC_lub6, stateMPC_dlubaff6);
stateMPC_LA_VSUB_INDEXED_29(stateMPC_dzaff7, stateMPC_lbIdx7, stateMPC_rilb7, stateMPC_dslbaff7);
stateMPC_LA_VSUB3_29(stateMPC_llbbyslb7, stateMPC_dslbaff7, stateMPC_llb7, stateMPC_dllbaff7);
stateMPC_LA_VSUB2_INDEXED_11(stateMPC_riub7, stateMPC_dzaff7, stateMPC_ubIdx7, stateMPC_dsubaff7);
stateMPC_LA_VSUB3_11(stateMPC_lubbysub7, stateMPC_dsubaff7, stateMPC_lub7, stateMPC_dlubaff7);
stateMPC_LA_VSUB_INDEXED_29(stateMPC_dzaff8, stateMPC_lbIdx8, stateMPC_rilb8, stateMPC_dslbaff8);
stateMPC_LA_VSUB3_29(stateMPC_llbbyslb8, stateMPC_dslbaff8, stateMPC_llb8, stateMPC_dllbaff8);
stateMPC_LA_VSUB2_INDEXED_11(stateMPC_riub8, stateMPC_dzaff8, stateMPC_ubIdx8, stateMPC_dsubaff8);
stateMPC_LA_VSUB3_11(stateMPC_lubbysub8, stateMPC_dsubaff8, stateMPC_lub8, stateMPC_dlubaff8);
stateMPC_LA_VSUB_INDEXED_9(stateMPC_dzaff9, stateMPC_lbIdx9, stateMPC_rilb9, stateMPC_dslbaff9);
stateMPC_LA_VSUB3_9(stateMPC_llbbyslb9, stateMPC_dslbaff9, stateMPC_llb9, stateMPC_dllbaff9);
stateMPC_LA_VSUB2_INDEXED_9(stateMPC_riub9, stateMPC_dzaff9, stateMPC_ubIdx9, stateMPC_dsubaff9);
stateMPC_LA_VSUB3_9(stateMPC_lubbysub9, stateMPC_dsubaff9, stateMPC_lub9, stateMPC_dlubaff9);
info->lsit_aff = stateMPC_LINESEARCH_BACKTRACKING_AFFINE(stateMPC_l, stateMPC_s, stateMPC_dl_aff, stateMPC_ds_aff, &info->step_aff, &info->mu_aff);
if( info->lsit_aff == stateMPC_NOPROGRESS ){
exitcode = stateMPC_NOPROGRESS; break;
}
sigma_3rdroot = info->mu_aff / info->mu;
info->sigma = sigma_3rdroot*sigma_3rdroot*sigma_3rdroot;
musigma = info->mu * info->sigma;
stateMPC_LA_VSUB5_378(stateMPC_ds_aff, stateMPC_dl_aff, musigma, stateMPC_ccrhs);
stateMPC_LA_VSUB6_INDEXED_29_11_29(stateMPC_ccrhsub0, stateMPC_sub0, stateMPC_ubIdx0, stateMPC_ccrhsl0, stateMPC_slb0, stateMPC_lbIdx0, stateMPC_rd0);
stateMPC_LA_VSUB6_INDEXED_29_11_29(stateMPC_ccrhsub1, stateMPC_sub1, stateMPC_ubIdx1, stateMPC_ccrhsl1, stateMPC_slb1, stateMPC_lbIdx1, stateMPC_rd1);
stateMPC_LA_DIAG_FORWARDSUB_29(stateMPC_Phi0, stateMPC_rd0, stateMPC_Lbyrd0);
stateMPC_LA_DIAG_FORWARDSUB_29(stateMPC_Phi1, stateMPC_rd1, stateMPC_Lbyrd1);
stateMPC_LA_DIAGZERO_MVM_9(stateMPC_W0, stateMPC_Lbyrd0, stateMPC_beta0);
stateMPC_LA_DENSE_FORWARDSUB_9(stateMPC_Ld0, stateMPC_beta0, stateMPC_yy0);
stateMPC_LA_DENSE_DIAGZERO_2MVMADD_9_29_29(stateMPC_V0, stateMPC_Lbyrd0, stateMPC_W1, stateMPC_Lbyrd1, stateMPC_beta1);
stateMPC_LA_DENSE_MVMSUB1_9_9(stateMPC_Lsd1, stateMPC_yy0, stateMPC_beta1, stateMPC_bmy1);
stateMPC_LA_DENSE_FORWARDSUB_9(stateMPC_Ld1, stateMPC_bmy1, stateMPC_yy1);
stateMPC_LA_VSUB6_INDEXED_29_11_29(stateMPC_ccrhsub2, stateMPC_sub2, stateMPC_ubIdx2, stateMPC_ccrhsl2, stateMPC_slb2, stateMPC_lbIdx2, stateMPC_rd2);
stateMPC_LA_DIAG_FORWARDSUB_29(stateMPC_Phi2, stateMPC_rd2, stateMPC_Lbyrd2);
stateMPC_LA_DENSE_DIAGZERO_2MVMADD_9_29_29(stateMPC_V1, stateMPC_Lbyrd1, stateMPC_W2, stateMPC_Lbyrd2, stateMPC_beta2);
stateMPC_LA_DENSE_MVMSUB1_9_9(stateMPC_Lsd2, stateMPC_yy1, stateMPC_beta2, stateMPC_bmy2);
stateMPC_LA_DENSE_FORWARDSUB_9(stateMPC_Ld2, stateMPC_bmy2, stateMPC_yy2);
stateMPC_LA_VSUB6_INDEXED_29_11_29(stateMPC_ccrhsub3, stateMPC_sub3, stateMPC_ubIdx3, stateMPC_ccrhsl3, stateMPC_slb3, stateMPC_lbIdx3, stateMPC_rd3);
stateMPC_LA_DIAG_FORWARDSUB_29(stateMPC_Phi3, stateMPC_rd3, stateMPC_Lbyrd3);
stateMPC_LA_DENSE_DIAGZERO_2MVMADD_9_29_29(stateMPC_V2, stateMPC_Lbyrd2, stateMPC_W3, stateMPC_Lbyrd3, stateMPC_beta3);
stateMPC_LA_DENSE_MVMSUB1_9_9(stateMPC_Lsd3, stateMPC_yy2, stateMPC_beta3, stateMPC_bmy3);
stateMPC_LA_DENSE_FORWARDSUB_9(stateMPC_Ld3, stateMPC_bmy3, stateMPC_yy3);
stateMPC_LA_VSUB6_INDEXED_29_11_29(stateMPC_ccrhsub4, stateMPC_sub4, stateMPC_ubIdx4, stateMPC_ccrhsl4, stateMPC_slb4, stateMPC_lbIdx4, stateMPC_rd4);
stateMPC_LA_DIAG_FORWARDSUB_29(stateMPC_Phi4, stateMPC_rd4, stateMPC_Lbyrd4);
stateMPC_LA_DENSE_DIAGZERO_2MVMADD_9_29_29(stateMPC_V3, stateMPC_Lbyrd3, stateMPC_W4, stateMPC_Lbyrd4, stateMPC_beta4);
stateMPC_LA_DENSE_MVMSUB1_9_9(stateMPC_Lsd4, stateMPC_yy3, stateMPC_beta4, stateMPC_bmy4);
stateMPC_LA_DENSE_FORWARDSUB_9(stateMPC_Ld4, stateMPC_bmy4, stateMPC_yy4);
stateMPC_LA_VSUB6_INDEXED_29_11_29(stateMPC_ccrhsub5, stateMPC_sub5, stateMPC_ubIdx5, stateMPC_ccrhsl5, stateMPC_slb5, stateMPC_lbIdx5, stateMPC_rd5);
stateMPC_LA_DIAG_FORWARDSUB_29(stateMPC_Phi5, stateMPC_rd5, stateMPC_Lbyrd5);
stateMPC_LA_DENSE_DIAGZERO_2MVMADD_9_29_29(stateMPC_V4, stateMPC_Lbyrd4, stateMPC_W5, stateMPC_Lbyrd5, stateMPC_beta5);
stateMPC_LA_DENSE_MVMSUB1_9_9(stateMPC_Lsd5, stateMPC_yy4, stateMPC_beta5, stateMPC_bmy5);
stateMPC_LA_DENSE_FORWARDSUB_9(stateMPC_Ld5, stateMPC_bmy5, stateMPC_yy5);
stateMPC_LA_VSUB6_INDEXED_29_11_29(stateMPC_ccrhsub6, stateMPC_sub6, stateMPC_ubIdx6, stateMPC_ccrhsl6, stateMPC_slb6, stateMPC_lbIdx6, stateMPC_rd6);
stateMPC_LA_DIAG_FORWARDSUB_29(stateMPC_Phi6, stateMPC_rd6, stateMPC_Lbyrd6);
stateMPC_LA_DENSE_DIAGZERO_2MVMADD_9_29_29(stateMPC_V5, stateMPC_Lbyrd5, stateMPC_W6, stateMPC_Lbyrd6, stateMPC_beta6);
stateMPC_LA_DENSE_MVMSUB1_9_9(stateMPC_Lsd6, stateMPC_yy5, stateMPC_beta6, stateMPC_bmy6);
stateMPC_LA_DENSE_FORWARDSUB_9(stateMPC_Ld6, stateMPC_bmy6, stateMPC_yy6);
stateMPC_LA_VSUB6_INDEXED_29_11_29(stateMPC_ccrhsub7, stateMPC_sub7, stateMPC_ubIdx7, stateMPC_ccrhsl7, stateMPC_slb7, stateMPC_lbIdx7, stateMPC_rd7);
stateMPC_LA_DIAG_FORWARDSUB_29(stateMPC_Phi7, stateMPC_rd7, stateMPC_Lbyrd7);
stateMPC_LA_DENSE_DIAGZERO_2MVMADD_9_29_29(stateMPC_V6, stateMPC_Lbyrd6, stateMPC_W7, stateMPC_Lbyrd7, stateMPC_beta7);
stateMPC_LA_DENSE_MVMSUB1_9_9(stateMPC_Lsd7, stateMPC_yy6, stateMPC_beta7, stateMPC_bmy7);
stateMPC_LA_DENSE_FORWARDSUB_9(stateMPC_Ld7, stateMPC_bmy7, stateMPC_yy7);
stateMPC_LA_VSUB6_INDEXED_29_11_29(stateMPC_ccrhsub8, stateMPC_sub8, stateMPC_ubIdx8, stateMPC_ccrhsl8, stateMPC_slb8, stateMPC_lbIdx8, stateMPC_rd8);
stateMPC_LA_DIAG_FORWARDSUB_29(stateMPC_Phi8, stateMPC_rd8, stateMPC_Lbyrd8);
stateMPC_LA_DENSE_DIAGZERO_2MVMADD_9_29_29(stateMPC_V7, stateMPC_Lbyrd7, stateMPC_W8, stateMPC_Lbyrd8, stateMPC_beta8);
stateMPC_LA_DENSE_MVMSUB1_9_9(stateMPC_Lsd8, stateMPC_yy7, stateMPC_beta8, stateMPC_bmy8);
stateMPC_LA_DENSE_FORWARDSUB_9(stateMPC_Ld8, stateMPC_bmy8, stateMPC_yy8);
stateMPC_LA_VSUB6_INDEXED_9_9_9(stateMPC_ccrhsub9, stateMPC_sub9, stateMPC_ubIdx9, stateMPC_ccrhsl9, stateMPC_slb9, stateMPC_lbIdx9, stateMPC_rd9);
stateMPC_LA_DIAG_FORWARDSUB_9(stateMPC_Phi9, stateMPC_rd9, stateMPC_Lbyrd9);
stateMPC_LA_DENSE_DIAGZERO_2MVMADD_9_29_9(stateMPC_V8, stateMPC_Lbyrd8, stateMPC_W9, stateMPC_Lbyrd9, stateMPC_beta9);
stateMPC_LA_DENSE_MVMSUB1_9_9(stateMPC_Lsd9, stateMPC_yy8, stateMPC_beta9, stateMPC_bmy9);
stateMPC_LA_DENSE_FORWARDSUB_9(stateMPC_Ld9, stateMPC_bmy9, stateMPC_yy9);
stateMPC_LA_DENSE_BACKWARDSUB_9(stateMPC_Ld9, stateMPC_yy9, stateMPC_dvcc9);
stateMPC_LA_DENSE_MTVMSUB_9_9(stateMPC_Lsd9, stateMPC_dvcc9, stateMPC_yy8, stateMPC_bmy8);
stateMPC_LA_DENSE_BACKWARDSUB_9(stateMPC_Ld8, stateMPC_bmy8, stateMPC_dvcc8);
stateMPC_LA_DENSE_MTVMSUB_9_9(stateMPC_Lsd8, stateMPC_dvcc8, stateMPC_yy7, stateMPC_bmy7);
stateMPC_LA_DENSE_BACKWARDSUB_9(stateMPC_Ld7, stateMPC_bmy7, stateMPC_dvcc7);
stateMPC_LA_DENSE_MTVMSUB_9_9(stateMPC_Lsd7, stateMPC_dvcc7, stateMPC_yy6, stateMPC_bmy6);
stateMPC_LA_DENSE_BACKWARDSUB_9(stateMPC_Ld6, stateMPC_bmy6, stateMPC_dvcc6);
stateMPC_LA_DENSE_MTVMSUB_9_9(stateMPC_Lsd6, stateMPC_dvcc6, stateMPC_yy5, stateMPC_bmy5);
stateMPC_LA_DENSE_BACKWARDSUB_9(stateMPC_Ld5, stateMPC_bmy5, stateMPC_dvcc5);
stateMPC_LA_DENSE_MTVMSUB_9_9(stateMPC_Lsd5, stateMPC_dvcc5, stateMPC_yy4, stateMPC_bmy4);
stateMPC_LA_DENSE_BACKWARDSUB_9(stateMPC_Ld4, stateMPC_bmy4, stateMPC_dvcc4);
stateMPC_LA_DENSE_MTVMSUB_9_9(stateMPC_Lsd4, stateMPC_dvcc4, stateMPC_yy3, stateMPC_bmy3);
stateMPC_LA_DENSE_BACKWARDSUB_9(stateMPC_Ld3, stateMPC_bmy3, stateMPC_dvcc3);
stateMPC_LA_DENSE_MTVMSUB_9_9(stateMPC_Lsd3, stateMPC_dvcc3, stateMPC_yy2, stateMPC_bmy2);
stateMPC_LA_DENSE_BACKWARDSUB_9(stateMPC_Ld2, stateMPC_bmy2, stateMPC_dvcc2);
stateMPC_LA_DENSE_MTVMSUB_9_9(stateMPC_Lsd2, stateMPC_dvcc2, stateMPC_yy1, stateMPC_bmy1);
stateMPC_LA_DENSE_BACKWARDSUB_9(stateMPC_Ld1, stateMPC_bmy1, stateMPC_dvcc1);
stateMPC_LA_DENSE_MTVMSUB_9_9(stateMPC_Lsd1, stateMPC_dvcc1, stateMPC_yy0, stateMPC_bmy0);
stateMPC_LA_DENSE_BACKWARDSUB_9(stateMPC_Ld0, stateMPC_bmy0, stateMPC_dvcc0);
stateMPC_LA_DENSE_DIAGZERO_MTVM2_9_29_9(params->C1, stateMPC_dvcc1, stateMPC_D0, stateMPC_dvcc0, stateMPC_grad_eq0);
stateMPC_LA_DENSE_DIAGZERO_MTVM2_9_29_9(params->C2, stateMPC_dvcc2, stateMPC_D1, stateMPC_dvcc1, stateMPC_grad_eq1);
stateMPC_LA_DENSE_DIAGZERO_MTVM2_9_29_9(params->C3, stateMPC_dvcc3, stateMPC_D1, stateMPC_dvcc2, stateMPC_grad_eq2);
stateMPC_LA_DENSE_DIAGZERO_MTVM2_9_29_9(params->C4, stateMPC_dvcc4, stateMPC_D1, stateMPC_dvcc3, stateMPC_grad_eq3);
stateMPC_LA_DENSE_DIAGZERO_MTVM2_9_29_9(params->C5, stateMPC_dvcc5, stateMPC_D1, stateMPC_dvcc4, stateMPC_grad_eq4);
stateMPC_LA_DENSE_DIAGZERO_MTVM2_9_29_9(params->C6, stateMPC_dvcc6, stateMPC_D1, stateMPC_dvcc5, stateMPC_grad_eq5);
stateMPC_LA_DENSE_DIAGZERO_MTVM2_9_29_9(params->C7, stateMPC_dvcc7, stateMPC_D1, stateMPC_dvcc6, stateMPC_grad_eq6);
stateMPC_LA_DENSE_DIAGZERO_MTVM2_9_29_9(params->C8, stateMPC_dvcc8, stateMPC_D1, stateMPC_dvcc7, stateMPC_grad_eq7);
stateMPC_LA_DENSE_DIAGZERO_MTVM2_9_29_9(params->C9, stateMPC_dvcc9, stateMPC_D1, stateMPC_dvcc8, stateMPC_grad_eq8);
stateMPC_LA_DIAGZERO_MTVM_9_9(stateMPC_D9, stateMPC_dvcc9, stateMPC_grad_eq9);
stateMPC_LA_VSUB_270(stateMPC_rd, stateMPC_grad_eq, stateMPC_rd);
stateMPC_LA_DIAG_FORWARDBACKWARDSUB_29(stateMPC_Phi0, stateMPC_rd0, stateMPC_dzcc0);
stateMPC_LA_DIAG_FORWARDBACKWARDSUB_29(stateMPC_Phi1, stateMPC_rd1, stateMPC_dzcc1);
stateMPC_LA_DIAG_FORWARDBACKWARDSUB_29(stateMPC_Phi2, stateMPC_rd2, stateMPC_dzcc2);
stateMPC_LA_DIAG_FORWARDBACKWARDSUB_29(stateMPC_Phi3, stateMPC_rd3, stateMPC_dzcc3);
stateMPC_LA_DIAG_FORWARDBACKWARDSUB_29(stateMPC_Phi4, stateMPC_rd4, stateMPC_dzcc4);
stateMPC_LA_DIAG_FORWARDBACKWARDSUB_29(stateMPC_Phi5, stateMPC_rd5, stateMPC_dzcc5);
stateMPC_LA_DIAG_FORWARDBACKWARDSUB_29(stateMPC_Phi6, stateMPC_rd6, stateMPC_dzcc6);
stateMPC_LA_DIAG_FORWARDBACKWARDSUB_29(stateMPC_Phi7, stateMPC_rd7, stateMPC_dzcc7);
stateMPC_LA_DIAG_FORWARDBACKWARDSUB_29(stateMPC_Phi8, stateMPC_rd8, stateMPC_dzcc8);
stateMPC_LA_DIAG_FORWARDBACKWARDSUB_9(stateMPC_Phi9, stateMPC_rd9, stateMPC_dzcc9);
stateMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_29(stateMPC_ccrhsl0, stateMPC_slb0, stateMPC_llbbyslb0, stateMPC_dzcc0, stateMPC_lbIdx0, stateMPC_dllbcc0);
stateMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_11(stateMPC_ccrhsub0, stateMPC_sub0, stateMPC_lubbysub0, stateMPC_dzcc0, stateMPC_ubIdx0, stateMPC_dlubcc0);
stateMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_29(stateMPC_ccrhsl1, stateMPC_slb1, stateMPC_llbbyslb1, stateMPC_dzcc1, stateMPC_lbIdx1, stateMPC_dllbcc1);
stateMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_11(stateMPC_ccrhsub1, stateMPC_sub1, stateMPC_lubbysub1, stateMPC_dzcc1, stateMPC_ubIdx1, stateMPC_dlubcc1);
stateMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_29(stateMPC_ccrhsl2, stateMPC_slb2, stateMPC_llbbyslb2, stateMPC_dzcc2, stateMPC_lbIdx2, stateMPC_dllbcc2);
stateMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_11(stateMPC_ccrhsub2, stateMPC_sub2, stateMPC_lubbysub2, stateMPC_dzcc2, stateMPC_ubIdx2, stateMPC_dlubcc2);
stateMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_29(stateMPC_ccrhsl3, stateMPC_slb3, stateMPC_llbbyslb3, stateMPC_dzcc3, stateMPC_lbIdx3, stateMPC_dllbcc3);
stateMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_11(stateMPC_ccrhsub3, stateMPC_sub3, stateMPC_lubbysub3, stateMPC_dzcc3, stateMPC_ubIdx3, stateMPC_dlubcc3);
stateMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_29(stateMPC_ccrhsl4, stateMPC_slb4, stateMPC_llbbyslb4, stateMPC_dzcc4, stateMPC_lbIdx4, stateMPC_dllbcc4);
stateMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_11(stateMPC_ccrhsub4, stateMPC_sub4, stateMPC_lubbysub4, stateMPC_dzcc4, stateMPC_ubIdx4, stateMPC_dlubcc4);
stateMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_29(stateMPC_ccrhsl5, stateMPC_slb5, stateMPC_llbbyslb5, stateMPC_dzcc5, stateMPC_lbIdx5, stateMPC_dllbcc5);
stateMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_11(stateMPC_ccrhsub5, stateMPC_sub5, stateMPC_lubbysub5, stateMPC_dzcc5, stateMPC_ubIdx5, stateMPC_dlubcc5);
stateMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_29(stateMPC_ccrhsl6, stateMPC_slb6, stateMPC_llbbyslb6, stateMPC_dzcc6, stateMPC_lbIdx6, stateMPC_dllbcc6);
stateMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_11(stateMPC_ccrhsub6, stateMPC_sub6, stateMPC_lubbysub6, stateMPC_dzcc6, stateMPC_ubIdx6, stateMPC_dlubcc6);
stateMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_29(stateMPC_ccrhsl7, stateMPC_slb7, stateMPC_llbbyslb7, stateMPC_dzcc7, stateMPC_lbIdx7, stateMPC_dllbcc7);
stateMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_11(stateMPC_ccrhsub7, stateMPC_sub7, stateMPC_lubbysub7, stateMPC_dzcc7, stateMPC_ubIdx7, stateMPC_dlubcc7);
stateMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_29(stateMPC_ccrhsl8, stateMPC_slb8, stateMPC_llbbyslb8, stateMPC_dzcc8, stateMPC_lbIdx8, stateMPC_dllbcc8);
stateMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_11(stateMPC_ccrhsub8, stateMPC_sub8, stateMPC_lubbysub8, stateMPC_dzcc8, stateMPC_ubIdx8, stateMPC_dlubcc8);
stateMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_9(stateMPC_ccrhsl9, stateMPC_slb9, stateMPC_llbbyslb9, stateMPC_dzcc9, stateMPC_lbIdx9, stateMPC_dllbcc9);
stateMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_9(stateMPC_ccrhsub9, stateMPC_sub9, stateMPC_lubbysub9, stateMPC_dzcc9, stateMPC_ubIdx9, stateMPC_dlubcc9);
stateMPC_LA_VSUB7_378(stateMPC_l, stateMPC_ccrhs, stateMPC_s, stateMPC_dl_cc, stateMPC_ds_cc);
stateMPC_LA_VADD_270(stateMPC_dz_cc, stateMPC_dz_aff);
stateMPC_LA_VADD_90(stateMPC_dv_cc, stateMPC_dv_aff);
stateMPC_LA_VADD_378(stateMPC_dl_cc, stateMPC_dl_aff);
stateMPC_LA_VADD_378(stateMPC_ds_cc, stateMPC_ds_aff);
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
output->z1[4] = stateMPC_z0[4];
output->z1[5] = stateMPC_z0[5];
output->z1[6] = stateMPC_z0[6];
output->z1[7] = stateMPC_z0[7];
output->z1[8] = stateMPC_z0[8];
output->z1[9] = stateMPC_z0[9];
output->z1[10] = stateMPC_z0[10];
output->z2[0] = stateMPC_z1[0];
output->z2[1] = stateMPC_z1[1];
output->z2[2] = stateMPC_z1[2];
output->z2[3] = stateMPC_z1[3];
output->z2[4] = stateMPC_z1[4];
output->z2[5] = stateMPC_z1[5];
output->z2[6] = stateMPC_z1[6];
output->z2[7] = stateMPC_z1[7];
output->z2[8] = stateMPC_z1[8];
output->z2[9] = stateMPC_z1[9];
output->z2[10] = stateMPC_z1[10];
output->z3[0] = stateMPC_z2[0];
output->z3[1] = stateMPC_z2[1];
output->z3[2] = stateMPC_z2[2];
output->z3[3] = stateMPC_z2[3];
output->z3[4] = stateMPC_z2[4];
output->z3[5] = stateMPC_z2[5];
output->z3[6] = stateMPC_z2[6];
output->z3[7] = stateMPC_z2[7];
output->z3[8] = stateMPC_z2[8];
output->z3[9] = stateMPC_z2[9];
output->z3[10] = stateMPC_z2[10];
output->z4[0] = stateMPC_z3[0];
output->z4[1] = stateMPC_z3[1];
output->z4[2] = stateMPC_z3[2];
output->z4[3] = stateMPC_z3[3];
output->z4[4] = stateMPC_z3[4];
output->z4[5] = stateMPC_z3[5];
output->z4[6] = stateMPC_z3[6];
output->z4[7] = stateMPC_z3[7];
output->z4[8] = stateMPC_z3[8];
output->z4[9] = stateMPC_z3[9];
output->z4[10] = stateMPC_z3[10];
output->z5[0] = stateMPC_z4[0];
output->z5[1] = stateMPC_z4[1];
output->z5[2] = stateMPC_z4[2];
output->z5[3] = stateMPC_z4[3];
output->z5[4] = stateMPC_z4[4];
output->z5[5] = stateMPC_z4[5];
output->z5[6] = stateMPC_z4[6];
output->z5[7] = stateMPC_z4[7];
output->z5[8] = stateMPC_z4[8];
output->z5[9] = stateMPC_z4[9];
output->z5[10] = stateMPC_z4[10];
output->z6[0] = stateMPC_z5[0];
output->z6[1] = stateMPC_z5[1];
output->z6[2] = stateMPC_z5[2];
output->z6[3] = stateMPC_z5[3];
output->z6[4] = stateMPC_z5[4];
output->z6[5] = stateMPC_z5[5];
output->z6[6] = stateMPC_z5[6];
output->z6[7] = stateMPC_z5[7];
output->z6[8] = stateMPC_z5[8];
output->z6[9] = stateMPC_z5[9];
output->z6[10] = stateMPC_z5[10];
output->z7[0] = stateMPC_z6[0];
output->z7[1] = stateMPC_z6[1];
output->z7[2] = stateMPC_z6[2];
output->z7[3] = stateMPC_z6[3];
output->z7[4] = stateMPC_z6[4];
output->z7[5] = stateMPC_z6[5];
output->z7[6] = stateMPC_z6[6];
output->z7[7] = stateMPC_z6[7];
output->z7[8] = stateMPC_z6[8];
output->z7[9] = stateMPC_z6[9];
output->z7[10] = stateMPC_z6[10];
output->z8[0] = stateMPC_z7[0];
output->z8[1] = stateMPC_z7[1];
output->z8[2] = stateMPC_z7[2];
output->z8[3] = stateMPC_z7[3];
output->z8[4] = stateMPC_z7[4];
output->z8[5] = stateMPC_z7[5];
output->z8[6] = stateMPC_z7[6];
output->z8[7] = stateMPC_z7[7];
output->z8[8] = stateMPC_z7[8];
output->z8[9] = stateMPC_z7[9];
output->z8[10] = stateMPC_z7[10];
output->z9[0] = stateMPC_z8[0];
output->z9[1] = stateMPC_z8[1];
output->z9[2] = stateMPC_z8[2];
output->z9[3] = stateMPC_z8[3];
output->z9[4] = stateMPC_z8[4];
output->z9[5] = stateMPC_z8[5];
output->z9[6] = stateMPC_z8[6];
output->z9[7] = stateMPC_z8[7];
output->z9[8] = stateMPC_z8[8];
output->z9[9] = stateMPC_z8[9];
output->z9[10] = stateMPC_z8[10];
output->z10[0] = stateMPC_z9[0];
output->z10[1] = stateMPC_z9[1];
output->z10[2] = stateMPC_z9[2];
output->z10[3] = stateMPC_z9[3];
output->z10[4] = stateMPC_z9[4];
output->z10[5] = stateMPC_z9[5];
output->z10[6] = stateMPC_z9[6];
output->z10[7] = stateMPC_z9[7];
output->z10[8] = stateMPC_z9[8];

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
