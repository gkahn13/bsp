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

#include "exploreMPC.h"

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
 * Initializes a vector of length 20 with a value.
 */
void exploreMPC_LA_INITIALIZEVECTOR_20(exploreMPC_FLOAT* vec, exploreMPC_FLOAT value)
{
	int i;
	for( i=0; i<20; i++ )
	{
		vec[i] = value;
	}
}


/*
 * Initializes a vector of length 12 with a value.
 */
void exploreMPC_LA_INITIALIZEVECTOR_12(exploreMPC_FLOAT* vec, exploreMPC_FLOAT value)
{
	int i;
	for( i=0; i<12; i++ )
	{
		vec[i] = value;
	}
}


/*
 * Initializes a vector of length 40 with a value.
 */
void exploreMPC_LA_INITIALIZEVECTOR_40(exploreMPC_FLOAT* vec, exploreMPC_FLOAT value)
{
	int i;
	for( i=0; i<40; i++ )
	{
		vec[i] = value;
	}
}


/* 
 * Calculates a dot product and adds it to a variable: z += x'*y; 
 * This function is for vectors of length 40.
 */
void exploreMPC_LA_DOTACC_40(exploreMPC_FLOAT *x, exploreMPC_FLOAT *y, exploreMPC_FLOAT *z)
{
	int i;
	for( i=0; i<40; i++ ){
		*z += x[i]*y[i];
	}
}


/*
 * Calculates the gradient and the value for a quadratic function 0.5*z'*H*z + f'*z
 *
 * INPUTS:     H  - Symmetric Hessian, diag matrix of size [8 x 8]
 *             f  - column vector of size 8
 *             z  - column vector of size 8
 *
 * OUTPUTS: grad  - gradient at z (= H*z + f), column vector of size 8
 *          value <-- value + 0.5*z'*H*z + f'*z (value will be modified)
 */
void exploreMPC_LA_DIAG_QUADFCN_8(exploreMPC_FLOAT* H, exploreMPC_FLOAT* f, exploreMPC_FLOAT* z, exploreMPC_FLOAT* grad, exploreMPC_FLOAT* value)
{
	int i;
	exploreMPC_FLOAT hz;	
	for( i=0; i<8; i++){
		hz = H[i]*z[i];
		grad[i] = hz + f[i];
		*value += 0.5*hz*z[i] + f[i]*z[i];
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
void exploreMPC_LA_DIAG_QUADFCN_4(exploreMPC_FLOAT* H, exploreMPC_FLOAT* f, exploreMPC_FLOAT* z, exploreMPC_FLOAT* grad, exploreMPC_FLOAT* value)
{
	int i;
	exploreMPC_FLOAT hz;	
	for( i=0; i<4; i++){
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
void exploreMPC_LA_DIAGZERO_MVMSUB6_4(exploreMPC_FLOAT *B, exploreMPC_FLOAT *u, exploreMPC_FLOAT *b, exploreMPC_FLOAT *l, exploreMPC_FLOAT *r, exploreMPC_FLOAT *z, exploreMPC_FLOAT *y)
{
	int i;
	exploreMPC_FLOAT Bu[4];
	exploreMPC_FLOAT norm = *y;
	exploreMPC_FLOAT lr = 0;

	/* do A*x + B*u first */
	for( i=0; i<4; i++ ){
		Bu[i] = B[i]*u[i];
	}	

	for( i=0; i<4; i++ ){
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
void exploreMPC_LA_DENSE_DIAGZERO_MVMSUB3_4_8_8(exploreMPC_FLOAT *A, exploreMPC_FLOAT *x, exploreMPC_FLOAT *B, exploreMPC_FLOAT *u, exploreMPC_FLOAT *b, exploreMPC_FLOAT *l, exploreMPC_FLOAT *r, exploreMPC_FLOAT *z, exploreMPC_FLOAT *y)
{
	int i;
	int j;
	int k = 0;
	exploreMPC_FLOAT AxBu[4];
	exploreMPC_FLOAT norm = *y;
	exploreMPC_FLOAT lr = 0;

	/* do A*x + B*u first */
	for( i=0; i<4; i++ ){
		AxBu[i] = A[k++]*x[0] + B[i]*u[i];
	}	

	for( j=1; j<8; j++ ){		
		for( i=0; i<4; i++ ){
			AxBu[i] += A[k++]*x[j];
		}
	}

	for( i=0; i<4; i++ ){
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
void exploreMPC_LA_DENSE_DIAGZERO_MVMSUB3_4_8_4(exploreMPC_FLOAT *A, exploreMPC_FLOAT *x, exploreMPC_FLOAT *B, exploreMPC_FLOAT *u, exploreMPC_FLOAT *b, exploreMPC_FLOAT *l, exploreMPC_FLOAT *r, exploreMPC_FLOAT *z, exploreMPC_FLOAT *y)
{
	int i;
	int j;
	int k = 0;
	exploreMPC_FLOAT AxBu[4];
	exploreMPC_FLOAT norm = *y;
	exploreMPC_FLOAT lr = 0;

	/* do A*x + B*u first */
	for( i=0; i<4; i++ ){
		AxBu[i] = A[k++]*x[0] + B[i]*u[i];
	}	

	for( j=1; j<8; j++ ){		
		for( i=0; i<4; i++ ){
			AxBu[i] += A[k++]*x[j];
		}
	}

	for( i=0; i<4; i++ ){
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
 * where A is of size [4 x 8] and stored in column major format.
 * and B is of size [4 x 8] and stored in diagzero format
 * Note the transposes of A and B!
 */
void exploreMPC_LA_DENSE_DIAGZERO_MTVM2_4_8_4(exploreMPC_FLOAT *A, exploreMPC_FLOAT *x, exploreMPC_FLOAT *B, exploreMPC_FLOAT *y, exploreMPC_FLOAT *z)
{
	int i;
	int j;
	int k = 0;
	for( i=0; i<4; i++ ){
		z[i] = 0;
		for( j=0; j<4; j++ ){
			z[i] += A[k++]*x[j];
		}
		z[i] += B[i]*y[i];
	}
	for( i=4 ;i<8; i++ ){
		z[i] = 0;
		for( j=0; j<4; j++ ){
			z[i] += A[k++]*x[j];
		}
	}
}


/*
 * Matrix vector multiplication y = M'*x where M is of size [4 x 4]
 * and stored in diagzero format. Note the transpose of M!
 */
void exploreMPC_LA_DIAGZERO_MTVM_4_4(exploreMPC_FLOAT *M, exploreMPC_FLOAT *x, exploreMPC_FLOAT *y)
{
	int i;
	for( i=0; i<4; i++ ){
		y[i] = M[i]*x[i];
	}
}


/*
 * Vector subtraction and addition.
 *	 Input: five vectors t, tidx, u, v, w and two scalars z and r
 *	 Output: y = t(tidx) - u + w
 *           z = z - v'*x;
 *           r = max([norm(y,inf), z]);
 * for vectors of length 8. Output z is of course scalar.
 */
void exploreMPC_LA_VSUBADD3_8(exploreMPC_FLOAT* t, exploreMPC_FLOAT* u, int* uidx, exploreMPC_FLOAT* v, exploreMPC_FLOAT* w, exploreMPC_FLOAT* y, exploreMPC_FLOAT* z, exploreMPC_FLOAT* r)
{
	int i;
	exploreMPC_FLOAT norm = *r;
	exploreMPC_FLOAT vx = 0;
	exploreMPC_FLOAT x;
	for( i=0; i<8; i++){
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
 * for vectors of length 8. Output z is of course scalar.
 */
void exploreMPC_LA_VSUBADD2_8(exploreMPC_FLOAT* t, int* tidx, exploreMPC_FLOAT* u, exploreMPC_FLOAT* v, exploreMPC_FLOAT* w, exploreMPC_FLOAT* y, exploreMPC_FLOAT* z, exploreMPC_FLOAT* r)
{
	int i;
	exploreMPC_FLOAT norm = *r;
	exploreMPC_FLOAT vx = 0;
	exploreMPC_FLOAT x;
	for( i=0; i<8; i++){
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
 * for vectors of length 4. Output z is of course scalar.
 */
void exploreMPC_LA_VSUBADD3_4(exploreMPC_FLOAT* t, exploreMPC_FLOAT* u, int* uidx, exploreMPC_FLOAT* v, exploreMPC_FLOAT* w, exploreMPC_FLOAT* y, exploreMPC_FLOAT* z, exploreMPC_FLOAT* r)
{
	int i;
	exploreMPC_FLOAT norm = *r;
	exploreMPC_FLOAT vx = 0;
	exploreMPC_FLOAT x;
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
void exploreMPC_LA_VSUBADD2_4(exploreMPC_FLOAT* t, int* tidx, exploreMPC_FLOAT* u, exploreMPC_FLOAT* v, exploreMPC_FLOAT* w, exploreMPC_FLOAT* y, exploreMPC_FLOAT* z, exploreMPC_FLOAT* r)
{
	int i;
	exploreMPC_FLOAT norm = *r;
	exploreMPC_FLOAT vx = 0;
	exploreMPC_FLOAT x;
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
 * Computes inequality constraints gradient-
 * Special function for box constraints of length 8
 * Returns also L/S, a value that is often used elsewhere.
 */
void exploreMPC_LA_INEQ_B_GRAD_8_8_8(exploreMPC_FLOAT *lu, exploreMPC_FLOAT *su, exploreMPC_FLOAT *ru, exploreMPC_FLOAT *ll, exploreMPC_FLOAT *sl, exploreMPC_FLOAT *rl, int* lbIdx, int* ubIdx, exploreMPC_FLOAT *grad, exploreMPC_FLOAT *lubysu, exploreMPC_FLOAT *llbysl)
{
	int i;
	for( i=0; i<8; i++ ){
		grad[i] = 0;
	}
	for( i=0; i<8; i++ ){		
		llbysl[i] = ll[i] / sl[i];
		grad[lbIdx[i]] -= llbysl[i]*rl[i];
	}
	for( i=0; i<8; i++ ){
		lubysu[i] = lu[i] / su[i];
		grad[ubIdx[i]] += lubysu[i]*ru[i];
	}
}


/*
 * Computes inequality constraints gradient-
 * Special function for box constraints of length 4
 * Returns also L/S, a value that is often used elsewhere.
 */
void exploreMPC_LA_INEQ_B_GRAD_4_4_4(exploreMPC_FLOAT *lu, exploreMPC_FLOAT *su, exploreMPC_FLOAT *ru, exploreMPC_FLOAT *ll, exploreMPC_FLOAT *sl, exploreMPC_FLOAT *rl, int* lbIdx, int* ubIdx, exploreMPC_FLOAT *grad, exploreMPC_FLOAT *lubysu, exploreMPC_FLOAT *llbysl)
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
 * Addition of three vectors  z = u + w + v
 * of length 20.
 */
void exploreMPC_LA_VVADD3_20(exploreMPC_FLOAT *u, exploreMPC_FLOAT *v, exploreMPC_FLOAT *w, exploreMPC_FLOAT *z)
{
	int i;
	for( i=0; i<20; i++ ){
		z[i] = u[i] + v[i] + w[i];
	}
}


/*
 * Special function to compute the diagonal cholesky factorization of the 
 * positive definite augmented Hessian for block size 8.
 *
 * Inputs: - H = diagonal cost Hessian in diagonal storage format
 *         - llbysl = L / S of lower bounds
 *         - lubysu = L / S of upper bounds
 *
 * Output: Phi = sqrt(H + diag(llbysl) + diag(lubysu))
 * where Phi is stored in diagonal storage format
 */
void exploreMPC_LA_DIAG_CHOL_ONELOOP_LBUB_8_8_8(exploreMPC_FLOAT *H, exploreMPC_FLOAT *llbysl, int* lbIdx, exploreMPC_FLOAT *lubysu, int* ubIdx, exploreMPC_FLOAT *Phi)


{
	int i;
	
	/* compute cholesky */
	for( i=0; i<8; i++ ){
		Phi[i] = H[i] + llbysl[i] + lubysu[i];

#if exploreMPC_SET_PRINTLEVEL > 0 && defined PRINTNUMERICALWARNINGS
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
 * where A is to be computed and is of size [4 x 8],
 * B is given and of size [4 x 8], L is a diagonal
 * matrix of size 4 stored in diagonal matrix 
 * storage format. Note the transpose of L has no impact!
 *
 * Result: A in column major storage format.
 *
 */
void exploreMPC_LA_DIAG_MATRIXFORWARDSUB_4_8(exploreMPC_FLOAT *L, exploreMPC_FLOAT *B, exploreMPC_FLOAT *A)
{
    int i,j;
	 int k = 0;

	for( j=0; j<8; j++){
		for( i=0; i<4; i++){
			A[k] = B[k]/L[j];
			k++;
		}
	}

}


/**
 * Forward substitution for the matrix equation A*L' = B
 * where A is to be computed and is of size [4 x 8],
 * B is given and of size [4 x 8], L is a diagonal
 *  matrix of size 8 stored in diagonal 
 * storage format. Note the transpose of L!
 *
 * Result: A in diagonalzero storage format.
 *
 */
void exploreMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_4_8(exploreMPC_FLOAT *L, exploreMPC_FLOAT *B, exploreMPC_FLOAT *A)
{
	int j;
    for( j=0; j<8; j++ ){   
		A[j] = B[j]/L[j];
     }
}


/**
 * Compute C = A*B' where 
 *
 *	size(A) = [4 x 8]
 *  size(B) = [4 x 8] in diagzero format
 * 
 * A and C matrices are stored in column major format.
 * 
 * 
 */
void exploreMPC_LA_DENSE_DIAGZERO_MMTM_4_8_4(exploreMPC_FLOAT *A, exploreMPC_FLOAT *B, exploreMPC_FLOAT *C)
{
    int i, j;
	
	for( i=0; i<4; i++ ){
		for( j=0; j<4; j++){
			C[j*4+i] = B[i*4+j]*A[i];
		}
	}

}


/**
 * Forward substitution to solve L*y = b where L is a
 * diagonal matrix in vector storage format.
 * 
 * The dimensions involved are 8.
 */
void exploreMPC_LA_DIAG_FORWARDSUB_8(exploreMPC_FLOAT *L, exploreMPC_FLOAT *b, exploreMPC_FLOAT *y)
{
    int i;

    for( i=0; i<8; i++ ){
		y[i] = b[i]/L[i];
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
void exploreMPC_LA_DIAG_CHOL_ONELOOP_LBUB_4_4_4(exploreMPC_FLOAT *H, exploreMPC_FLOAT *llbysl, int* lbIdx, exploreMPC_FLOAT *lubysu, int* ubIdx, exploreMPC_FLOAT *Phi)


{
	int i;
	
	/* compute cholesky */
	for( i=0; i<4; i++ ){
		Phi[i] = H[i] + llbysl[i] + lubysu[i];

#if exploreMPC_SET_PRINTLEVEL > 0 && defined PRINTNUMERICALWARNINGS
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
 * where A is to be computed and is of size [4 x 4],
 * B is given and of size [4 x 4], L is a diagonal
 *  matrix of size 4 stored in diagonal 
 * storage format. Note the transpose of L!
 *
 * Result: A in diagonalzero storage format.
 *
 */
void exploreMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_4_4(exploreMPC_FLOAT *L, exploreMPC_FLOAT *B, exploreMPC_FLOAT *A)
{
	int j;
    for( j=0; j<4; j++ ){   
		A[j] = B[j]/L[j];
     }
}


/**
 * Forward substitution to solve L*y = b where L is a
 * diagonal matrix in vector storage format.
 * 
 * The dimensions involved are 4.
 */
void exploreMPC_LA_DIAG_FORWARDSUB_4(exploreMPC_FLOAT *L, exploreMPC_FLOAT *b, exploreMPC_FLOAT *y)
{
    int i;

    for( i=0; i<4; i++ ){
		y[i] = b[i]/L[i];
    }
}


/**
 * Compute L = B*B', where L is lower triangular of size NXp1
 * and B is a diagzero matrix of size [4 x 8] in column
 * storage format.
 * 
 */
void exploreMPC_LA_DIAGZERO_MMT_4(exploreMPC_FLOAT *B, exploreMPC_FLOAT *L)
{
    int i, ii, di;
    
    ii = 0; di = 0;
    for( i=0; i<4; i++ ){        
		L[ii+i] = B[i]*B[i];
        ii += ++di;
    }
}


/* 
 * Computes r = b - B*u
 * B is stored in diagzero format
 */
void exploreMPC_LA_DIAGZERO_MVMSUB7_4(exploreMPC_FLOAT *B, exploreMPC_FLOAT *u, exploreMPC_FLOAT *b, exploreMPC_FLOAT *r)
{
	int i;

	for( i=0; i<4; i++ ){
		r[i] = b[i] - B[i]*u[i];
	}	
	
}


/**
 * Compute L = A*A' + B*B', where L is lower triangular of size NXp1
 * and A is a dense matrix of size [4 x 8] in column
 * storage format, and B is of size [4 x 8] diagonalzero
 * storage format.
 * 
 * THIS ONE HAS THE WORST ACCES PATTERN POSSIBLE. 
 * POSSIBKE FIX: PUT A AND B INTO ROW MAJOR FORMAT FIRST.
 * 
 */
void exploreMPC_LA_DENSE_DIAGZERO_MMT2_4_8_8(exploreMPC_FLOAT *A, exploreMPC_FLOAT *B, exploreMPC_FLOAT *L)
{
    int i, j, k, ii, di;
    exploreMPC_FLOAT ltemp;
    
    ii = 0; di = 0;
    for( i=0; i<4; i++ ){        
        for( j=0; j<=i; j++ ){
            ltemp = 0; 
            for( k=0; k<8; k++ ){
                ltemp += A[k*4+i]*A[k*4+j];
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
void exploreMPC_LA_DENSE_DIAGZERO_2MVMSUB2_4_8_8(exploreMPC_FLOAT *A, exploreMPC_FLOAT *x, exploreMPC_FLOAT *B, exploreMPC_FLOAT *u, exploreMPC_FLOAT *b, exploreMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;

	for( i=0; i<4; i++ ){
		r[i] = b[i] - A[k++]*x[0] - B[i]*u[i];
	}	

	for( j=1; j<8; j++ ){		
		for( i=0; i<4; i++ ){
			r[i] -= A[k++]*x[j];
		}
	}
	
}


/**
 * Compute L = A*A' + B*B', where L is lower triangular of size NXp1
 * and A is a dense matrix of size [4 x 8] in column
 * storage format, and B is of size [4 x 4] diagonalzero
 * storage format.
 * 
 * THIS ONE HAS THE WORST ACCES PATTERN POSSIBLE. 
 * POSSIBKE FIX: PUT A AND B INTO ROW MAJOR FORMAT FIRST.
 * 
 */
void exploreMPC_LA_DENSE_DIAGZERO_MMT2_4_8_4(exploreMPC_FLOAT *A, exploreMPC_FLOAT *B, exploreMPC_FLOAT *L)
{
    int i, j, k, ii, di;
    exploreMPC_FLOAT ltemp;
    
    ii = 0; di = 0;
    for( i=0; i<4; i++ ){        
        for( j=0; j<=i; j++ ){
            ltemp = 0; 
            for( k=0; k<8; k++ ){
                ltemp += A[k*4+i]*A[k*4+j];
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
void exploreMPC_LA_DENSE_DIAGZERO_2MVMSUB2_4_8_4(exploreMPC_FLOAT *A, exploreMPC_FLOAT *x, exploreMPC_FLOAT *B, exploreMPC_FLOAT *u, exploreMPC_FLOAT *b, exploreMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;

	for( i=0; i<4; i++ ){
		r[i] = b[i] - A[k++]*x[0] - B[i]*u[i];
	}	

	for( j=1; j<8; j++ ){		
		for( i=0; i<4; i++ ){
			r[i] -= A[k++]*x[j];
		}
	}
	
}


/**
 * Cholesky factorization as above, but working on a matrix in 
 * lower triangular storage format of size 4 and outputting
 * the Cholesky factor to matrix L in lower triangular format.
 */
void exploreMPC_LA_DENSE_CHOL_4(exploreMPC_FLOAT *A, exploreMPC_FLOAT *L)
{
    int i, j, k, di, dj;
	 int ii, jj;

    exploreMPC_FLOAT l;
    exploreMPC_FLOAT Mii;

	/* copy A to L first and then operate on L */
	/* COULD BE OPTIMIZED */
	ii=0; di=0;
	for( i=0; i<4; i++ ){
		for( j=0; j<=i; j++ ){
			L[ii+j] = A[ii+j];
		}
		ii += ++di;
	}    
	
	/* factor L */
	ii=0; di=0;
    for( i=0; i<4; i++ ){
        l = 0;
        for( k=0; k<i; k++ ){
            l += L[ii+k]*L[ii+k];
        }        
        
        Mii = L[ii+i] - l;
        
#if exploreMPC_SET_PRINTLEVEL > 0 && defined PRINTNUMERICALWARNINGS
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
        for( j=i+1; j<4; j++ ){
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
 * The dimensions involved are 4.
 */
void exploreMPC_LA_DENSE_FORWARDSUB_4(exploreMPC_FLOAT *L, exploreMPC_FLOAT *b, exploreMPC_FLOAT *y)
{
    int i,j,ii,di;
    exploreMPC_FLOAT yel;
            
    ii = 0; di = 0;
    for( i=0; i<4; i++ ){
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
 * where A is to be computed and is of size [4 x 4],
 * B is given and of size [4 x 4], L is a lower tri-
 * angular matrix of size 4 stored in lower triangular 
 * storage format. Note the transpose of L AND B!
 *
 * Result: A in column major storage format.
 *
 */
void exploreMPC_LA_DENSE_MATRIXTFORWARDSUB_4_4(exploreMPC_FLOAT *L, exploreMPC_FLOAT *B, exploreMPC_FLOAT *A)
{
    int i,j,k,ii,di;
    exploreMPC_FLOAT a;
    
    ii=0; di=0;
    for( j=0; j<4; j++ ){        
        for( i=0; i<4; i++ ){
            a = B[i*4+j];
            for( k=0; k<j; k++ ){
                a -= A[k*4+i]*L[ii+k];
            }    

			/* saturate for numerical stability */
			a = MIN(a, BIGM);
			a = MAX(a, -BIGM); 

			A[j*4+i] = a/L[ii+j];			
        }
        ii += ++di;
    }
}


/**
 * Compute L = L - A*A', where L is lower triangular of size 4
 * and A is a dense matrix of size [4 x 4] in column
 * storage format.
 * 
 * THIS ONE HAS THE WORST ACCES PATTERN POSSIBLE. 
 * POSSIBKE FIX: PUT A INTO ROW MAJOR FORMAT FIRST.
 * 
 */
void exploreMPC_LA_DENSE_MMTSUB_4_4(exploreMPC_FLOAT *A, exploreMPC_FLOAT *L)
{
    int i, j, k, ii, di;
    exploreMPC_FLOAT ltemp;
    
    ii = 0; di = 0;
    for( i=0; i<4; i++ ){        
        for( j=0; j<=i; j++ ){
            ltemp = 0; 
            for( k=0; k<4; k++ ){
                ltemp += A[k*4+i]*A[k*4+j];
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
void exploreMPC_LA_DENSE_MVMSUB1_4_4(exploreMPC_FLOAT *A, exploreMPC_FLOAT *x, exploreMPC_FLOAT *b, exploreMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;

	for( i=0; i<4; i++ ){
		r[i] = b[i] - A[k++]*x[0];
	}	
	for( j=1; j<4; j++ ){		
		for( i=0; i<4; i++ ){
			r[i] -= A[k++]*x[j];
		}
	}
}


/**
 * Backward Substitution to solve L^T*x = y where L is a
 * lower triangular matrix in triangular storage format.
 * 
 * All involved dimensions are 4.
 */
void exploreMPC_LA_DENSE_BACKWARDSUB_4(exploreMPC_FLOAT *L, exploreMPC_FLOAT *y, exploreMPC_FLOAT *x)
{
    int i, ii, di, j, jj, dj;
    exploreMPC_FLOAT xel;    
	int start = 6;
    
    /* now solve L^T*x = y by backward substitution */
    ii = start; di = 3;
    for( i=3; i>=0; i-- ){        
        xel = y[i];        
        jj = start; dj = 3;
        for( j=3; j>i; j-- ){
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
 * Matrix vector multiplication y = b - M'*x where M is of size [4 x 4]
 * and stored in column major format. Note the transpose of M!
 */
void exploreMPC_LA_DENSE_MTVMSUB_4_4(exploreMPC_FLOAT *A, exploreMPC_FLOAT *x, exploreMPC_FLOAT *b, exploreMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0; 
	for( i=0; i<4; i++ ){
		r[i] = b[i];
		for( j=0; j<4; j++ ){
			r[i] -= A[k++]*x[j];
		}
	}
}


/*
 * Vector subtraction z = -x - y for vectors of length 20.
 */
void exploreMPC_LA_VSUB2_20(exploreMPC_FLOAT *x, exploreMPC_FLOAT *y, exploreMPC_FLOAT *z)
{
	int i;
	for( i=0; i<20; i++){
		z[i] = -x[i] - y[i];
	}
}


/**
 * Forward-Backward-Substitution to solve L*L^T*x = b where L is a
 * diagonal matrix of size 8 in vector
 * storage format.
 */
void exploreMPC_LA_DIAG_FORWARDBACKWARDSUB_8(exploreMPC_FLOAT *L, exploreMPC_FLOAT *b, exploreMPC_FLOAT *x)
{
    int i;
            
    /* solve Ly = b by forward and backward substitution */
    for( i=0; i<8; i++ ){
		x[i] = b[i]/(L[i]*L[i]);
    }
    
}


/**
 * Forward-Backward-Substitution to solve L*L^T*x = b where L is a
 * diagonal matrix of size 4 in vector
 * storage format.
 */
void exploreMPC_LA_DIAG_FORWARDBACKWARDSUB_4(exploreMPC_FLOAT *L, exploreMPC_FLOAT *b, exploreMPC_FLOAT *x)
{
    int i;
            
    /* solve Ly = b by forward and backward substitution */
    for( i=0; i<4; i++ ){
		x[i] = b[i]/(L[i]*L[i]);
    }
    
}


/*
 * Vector subtraction z = x(xidx) - y where y, z and xidx are of length 8,
 * and x has length 8 and is indexed through yidx.
 */
void exploreMPC_LA_VSUB_INDEXED_8(exploreMPC_FLOAT *x, int* xidx, exploreMPC_FLOAT *y, exploreMPC_FLOAT *z)
{
	int i;
	for( i=0; i<8; i++){
		z[i] = x[xidx[i]] - y[i];
	}
}


/*
 * Vector subtraction x = -u.*v - w for vectors of length 8.
 */
void exploreMPC_LA_VSUB3_8(exploreMPC_FLOAT *u, exploreMPC_FLOAT *v, exploreMPC_FLOAT *w, exploreMPC_FLOAT *x)
{
	int i;
	for( i=0; i<8; i++){
		x[i] = -u[i]*v[i] - w[i];
	}
}


/*
 * Vector subtraction z = -x - y(yidx) where y is of length 8
 * and z, x and yidx are of length 8.
 */
void exploreMPC_LA_VSUB2_INDEXED_8(exploreMPC_FLOAT *x, exploreMPC_FLOAT *y, int* yidx, exploreMPC_FLOAT *z)
{
	int i;
	for( i=0; i<8; i++){
		z[i] = -x[i] - y[yidx[i]];
	}
}


/*
 * Vector subtraction z = x(xidx) - y where y, z and xidx are of length 4,
 * and x has length 4 and is indexed through yidx.
 */
void exploreMPC_LA_VSUB_INDEXED_4(exploreMPC_FLOAT *x, int* xidx, exploreMPC_FLOAT *y, exploreMPC_FLOAT *z)
{
	int i;
	for( i=0; i<4; i++){
		z[i] = x[xidx[i]] - y[i];
	}
}


/*
 * Vector subtraction x = -u.*v - w for vectors of length 4.
 */
void exploreMPC_LA_VSUB3_4(exploreMPC_FLOAT *u, exploreMPC_FLOAT *v, exploreMPC_FLOAT *w, exploreMPC_FLOAT *x)
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
void exploreMPC_LA_VSUB2_INDEXED_4(exploreMPC_FLOAT *x, exploreMPC_FLOAT *y, int* yidx, exploreMPC_FLOAT *z)
{
	int i;
	for( i=0; i<4; i++){
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
 * exploreMPC_NOPROGRESS (should be negative).
 */
int exploreMPC_LINESEARCH_BACKTRACKING_AFFINE(exploreMPC_FLOAT *l, exploreMPC_FLOAT *s, exploreMPC_FLOAT *dl, exploreMPC_FLOAT *ds, exploreMPC_FLOAT *a, exploreMPC_FLOAT *mu_aff)
{
    int i;
	int lsIt=1;    
    exploreMPC_FLOAT dltemp;
    exploreMPC_FLOAT dstemp;
    exploreMPC_FLOAT mya = 1.0;
    exploreMPC_FLOAT mymu;
        
    while( 1 ){                        

        /* 
         * Compute both snew and wnew together.
         * We compute also mu_affine along the way here, as the
         * values might be in registers, so it should be cheaper.
         */
        mymu = 0;
        for( i=0; i<40; i++ ){
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
        if( i == 40 ){
            break;
        } else {
            mya *= exploreMPC_SET_LS_SCALE_AFF;
            if( mya < exploreMPC_SET_LS_MINSTEP ){
                return exploreMPC_NOPROGRESS;
            }
        }
    }
    
    /* return new values and iteration counter */
    *a = mya;
    *mu_aff = mymu / (exploreMPC_FLOAT)40;
    return lsIt;
}


/*
 * Vector subtraction x = u.*v - a where a is a scalar
*  and x,u,v are vectors of length 40.
 */
void exploreMPC_LA_VSUB5_40(exploreMPC_FLOAT *u, exploreMPC_FLOAT *v, exploreMPC_FLOAT a, exploreMPC_FLOAT *x)
{
	int i;
	for( i=0; i<40; i++){
		x[i] = u[i]*v[i] - a;
	}
}


/*
 * Computes x=0; x(uidx) += u/su; x(vidx) -= v/sv where x is of length 8,
 * u, su, uidx are of length 8 and v, sv, vidx are of length 8.
 */
void exploreMPC_LA_VSUB6_INDEXED_8_8_8(exploreMPC_FLOAT *u, exploreMPC_FLOAT *su, int* uidx, exploreMPC_FLOAT *v, exploreMPC_FLOAT *sv, int* vidx, exploreMPC_FLOAT *x)
{
	int i;
	for( i=0; i<8; i++ ){
		x[i] = 0;
	}
	for( i=0; i<8; i++){
		x[uidx[i]] += u[i]/su[i];
	}
	for( i=0; i<8; i++){
		x[vidx[i]] -= v[i]/sv[i];
	}
}


/* 
 * Computes r =  B*u
 * where B is stored in diagzero format
 */
void exploreMPC_LA_DIAGZERO_MVM_4(exploreMPC_FLOAT *B, exploreMPC_FLOAT *u, exploreMPC_FLOAT *r)
{
	int i;

	for( i=0; i<4; i++ ){
		r[i] = B[i]*u[i];
	}	
	
}


/* 
 * Computes r = A*x + B*u
 * where A is stored in column major format
 * and B is stored in diagzero format
 */
void exploreMPC_LA_DENSE_DIAGZERO_2MVMADD_4_8_8(exploreMPC_FLOAT *A, exploreMPC_FLOAT *x, exploreMPC_FLOAT *B, exploreMPC_FLOAT *u, exploreMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;

	for( i=0; i<4; i++ ){
		r[i] = A[k++]*x[0] + B[i]*u[i];
	}	

	for( j=1; j<8; j++ ){		
		for( i=0; i<4; i++ ){
			r[i] += A[k++]*x[j];
		}
	}
	
}


/*
 * Computes x=0; x(uidx) += u/su; x(vidx) -= v/sv where x is of length 4,
 * u, su, uidx are of length 4 and v, sv, vidx are of length 4.
 */
void exploreMPC_LA_VSUB6_INDEXED_4_4_4(exploreMPC_FLOAT *u, exploreMPC_FLOAT *su, int* uidx, exploreMPC_FLOAT *v, exploreMPC_FLOAT *sv, int* vidx, exploreMPC_FLOAT *x)
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
 * Computes r = A*x + B*u
 * where A is stored in column major format
 * and B is stored in diagzero format
 */
void exploreMPC_LA_DENSE_DIAGZERO_2MVMADD_4_8_4(exploreMPC_FLOAT *A, exploreMPC_FLOAT *x, exploreMPC_FLOAT *B, exploreMPC_FLOAT *u, exploreMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;

	for( i=0; i<4; i++ ){
		r[i] = A[k++]*x[0] + B[i]*u[i];
	}	

	for( j=1; j<8; j++ ){		
		for( i=0; i<4; i++ ){
			r[i] += A[k++]*x[j];
		}
	}
	
}


/*
 * Vector subtraction z = x - y for vectors of length 20.
 */
void exploreMPC_LA_VSUB_20(exploreMPC_FLOAT *x, exploreMPC_FLOAT *y, exploreMPC_FLOAT *z)
{
	int i;
	for( i=0; i<20; i++){
		z[i] = x[i] - y[i];
	}
}


/** 
 * Computes z = -r./s - u.*y(y)
 * where all vectors except of y are of length 8 (length of y >= 8).
 */
void exploreMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_8(exploreMPC_FLOAT *r, exploreMPC_FLOAT *s, exploreMPC_FLOAT *u, exploreMPC_FLOAT *y, int* yidx, exploreMPC_FLOAT *z)
{
	int i;
	for( i=0; i<8; i++ ){
		z[i] = -r[i]/s[i] - u[i]*y[yidx[i]];
	}
}


/** 
 * Computes z = -r./s + u.*y(y)
 * where all vectors except of y are of length 8 (length of y >= 8).
 */
void exploreMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_8(exploreMPC_FLOAT *r, exploreMPC_FLOAT *s, exploreMPC_FLOAT *u, exploreMPC_FLOAT *y, int* yidx, exploreMPC_FLOAT *z)
{
	int i;
	for( i=0; i<8; i++ ){
		z[i] = -r[i]/s[i] + u[i]*y[yidx[i]];
	}
}


/** 
 * Computes z = -r./s - u.*y(y)
 * where all vectors except of y are of length 4 (length of y >= 4).
 */
void exploreMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_4(exploreMPC_FLOAT *r, exploreMPC_FLOAT *s, exploreMPC_FLOAT *u, exploreMPC_FLOAT *y, int* yidx, exploreMPC_FLOAT *z)
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
void exploreMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_4(exploreMPC_FLOAT *r, exploreMPC_FLOAT *s, exploreMPC_FLOAT *u, exploreMPC_FLOAT *y, int* yidx, exploreMPC_FLOAT *z)
{
	int i;
	for( i=0; i<4; i++ ){
		z[i] = -r[i]/s[i] + u[i]*y[yidx[i]];
	}
}


/*
 * Computes ds = -l.\(r + s.*dl) for vectors of length 40.
 */
void exploreMPC_LA_VSUB7_40(exploreMPC_FLOAT *l, exploreMPC_FLOAT *r, exploreMPC_FLOAT *s, exploreMPC_FLOAT *dl, exploreMPC_FLOAT *ds)
{
	int i;
	for( i=0; i<40; i++){
		ds[i] = -(r[i] + s[i]*dl[i])/l[i];
	}
}


/*
 * Vector addition x = x + y for vectors of length 20.
 */
void exploreMPC_LA_VADD_20(exploreMPC_FLOAT *x, exploreMPC_FLOAT *y)
{
	int i;
	for( i=0; i<20; i++){
		x[i] += y[i];
	}
}


/*
 * Vector addition x = x + y for vectors of length 12.
 */
void exploreMPC_LA_VADD_12(exploreMPC_FLOAT *x, exploreMPC_FLOAT *y)
{
	int i;
	for( i=0; i<12; i++){
		x[i] += y[i];
	}
}


/*
 * Vector addition x = x + y for vectors of length 40.
 */
void exploreMPC_LA_VADD_40(exploreMPC_FLOAT *x, exploreMPC_FLOAT *y)
{
	int i;
	for( i=0; i<40; i++){
		x[i] += y[i];
	}
}


/**
 * Backtracking line search for combined predictor/corrector step.
 * Update on variables with safety factor gamma (to keep us away from
 * boundary).
 */
int exploreMPC_LINESEARCH_BACKTRACKING_COMBINED(exploreMPC_FLOAT *z, exploreMPC_FLOAT *v, exploreMPC_FLOAT *l, exploreMPC_FLOAT *s, exploreMPC_FLOAT *dz, exploreMPC_FLOAT *dv, exploreMPC_FLOAT *dl, exploreMPC_FLOAT *ds, exploreMPC_FLOAT *a, exploreMPC_FLOAT *mu)
{
    int i, lsIt=1;       
    exploreMPC_FLOAT dltemp;
    exploreMPC_FLOAT dstemp;    
    exploreMPC_FLOAT a_gamma;
            
    *a = 1.0;
    while( 1 ){                        

        /* check whether search criterion is fulfilled */
        for( i=0; i<40; i++ ){
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
        if( i == 40 ){
            break;
        } else {
            *a *= exploreMPC_SET_LS_SCALE;
            if( *a < exploreMPC_SET_LS_MINSTEP ){
                return exploreMPC_NOPROGRESS;
            }
        }
    }
    
    /* update variables with safety margin */
    a_gamma = (*a)*exploreMPC_SET_LS_MAXSTEP;
    
    /* primal variables */
    for( i=0; i<20; i++ ){
        z[i] += a_gamma*dz[i];
    }
    
    /* equality constraint multipliers */
    for( i=0; i<12; i++ ){
        v[i] += a_gamma*dv[i];
    }
    
    /* inequality constraint multipliers & slacks, also update mu */
    *mu = 0;
    for( i=0; i<40; i++ ){
        dltemp = l[i] + a_gamma*dl[i]; l[i] = dltemp;
        dstemp = s[i] + a_gamma*ds[i]; s[i] = dstemp;
        *mu += dltemp*dstemp;
    }
    
    *a = a_gamma;
    *mu /= (exploreMPC_FLOAT)40;
    return lsIt;
}




/* VARIABLE DEFINITIONS ------------------------------------------------ */
exploreMPC_FLOAT exploreMPC_z[20];
exploreMPC_FLOAT exploreMPC_v[12];
exploreMPC_FLOAT exploreMPC_dz_aff[20];
exploreMPC_FLOAT exploreMPC_dv_aff[12];
exploreMPC_FLOAT exploreMPC_grad_cost[20];
exploreMPC_FLOAT exploreMPC_grad_eq[20];
exploreMPC_FLOAT exploreMPC_rd[20];
exploreMPC_FLOAT exploreMPC_l[40];
exploreMPC_FLOAT exploreMPC_s[40];
exploreMPC_FLOAT exploreMPC_lbys[40];
exploreMPC_FLOAT exploreMPC_dl_aff[40];
exploreMPC_FLOAT exploreMPC_ds_aff[40];
exploreMPC_FLOAT exploreMPC_dz_cc[20];
exploreMPC_FLOAT exploreMPC_dv_cc[12];
exploreMPC_FLOAT exploreMPC_dl_cc[40];
exploreMPC_FLOAT exploreMPC_ds_cc[40];
exploreMPC_FLOAT exploreMPC_ccrhs[40];
exploreMPC_FLOAT exploreMPC_grad_ineq[20];
exploreMPC_FLOAT* exploreMPC_z0 = exploreMPC_z + 0;
exploreMPC_FLOAT* exploreMPC_dzaff0 = exploreMPC_dz_aff + 0;
exploreMPC_FLOAT* exploreMPC_dzcc0 = exploreMPC_dz_cc + 0;
exploreMPC_FLOAT* exploreMPC_rd0 = exploreMPC_rd + 0;
exploreMPC_FLOAT exploreMPC_Lbyrd0[8];
exploreMPC_FLOAT* exploreMPC_grad_cost0 = exploreMPC_grad_cost + 0;
exploreMPC_FLOAT* exploreMPC_grad_eq0 = exploreMPC_grad_eq + 0;
exploreMPC_FLOAT* exploreMPC_grad_ineq0 = exploreMPC_grad_ineq + 0;
exploreMPC_FLOAT exploreMPC_ctv0[8];
exploreMPC_FLOAT exploreMPC_C0[32] = {1.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 
0.0000000000000000E+000, 1.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 
0.0000000000000000E+000, 0.0000000000000000E+000, 1.0000000000000000E+000, 0.0000000000000000E+000, 
0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 1.0000000000000000E+000, 
1.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 
0.0000000000000000E+000, 1.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 
0.0000000000000000E+000, 0.0000000000000000E+000, 1.0000000000000000E+000, 0.0000000000000000E+000, 
0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 1.0000000000000000E+000};
exploreMPC_FLOAT* exploreMPC_v0 = exploreMPC_v + 0;
exploreMPC_FLOAT exploreMPC_re0[4];
exploreMPC_FLOAT exploreMPC_beta0[4];
exploreMPC_FLOAT exploreMPC_betacc0[4];
exploreMPC_FLOAT* exploreMPC_dvaff0 = exploreMPC_dv_aff + 0;
exploreMPC_FLOAT* exploreMPC_dvcc0 = exploreMPC_dv_cc + 0;
exploreMPC_FLOAT exploreMPC_V0[32];
exploreMPC_FLOAT exploreMPC_Yd0[10];
exploreMPC_FLOAT exploreMPC_Ld0[10];
exploreMPC_FLOAT exploreMPC_yy0[4];
exploreMPC_FLOAT exploreMPC_bmy0[4];
int exploreMPC_lbIdx0[8] = {0, 1, 2, 3, 4, 5, 6, 7};
exploreMPC_FLOAT* exploreMPC_llb0 = exploreMPC_l + 0;
exploreMPC_FLOAT* exploreMPC_slb0 = exploreMPC_s + 0;
exploreMPC_FLOAT* exploreMPC_llbbyslb0 = exploreMPC_lbys + 0;
exploreMPC_FLOAT exploreMPC_rilb0[8];
exploreMPC_FLOAT* exploreMPC_dllbaff0 = exploreMPC_dl_aff + 0;
exploreMPC_FLOAT* exploreMPC_dslbaff0 = exploreMPC_ds_aff + 0;
exploreMPC_FLOAT* exploreMPC_dllbcc0 = exploreMPC_dl_cc + 0;
exploreMPC_FLOAT* exploreMPC_dslbcc0 = exploreMPC_ds_cc + 0;
exploreMPC_FLOAT* exploreMPC_ccrhsl0 = exploreMPC_ccrhs + 0;
int exploreMPC_ubIdx0[8] = {0, 1, 2, 3, 4, 5, 6, 7};
exploreMPC_FLOAT* exploreMPC_lub0 = exploreMPC_l + 8;
exploreMPC_FLOAT* exploreMPC_sub0 = exploreMPC_s + 8;
exploreMPC_FLOAT* exploreMPC_lubbysub0 = exploreMPC_lbys + 8;
exploreMPC_FLOAT exploreMPC_riub0[8];
exploreMPC_FLOAT* exploreMPC_dlubaff0 = exploreMPC_dl_aff + 8;
exploreMPC_FLOAT* exploreMPC_dsubaff0 = exploreMPC_ds_aff + 8;
exploreMPC_FLOAT* exploreMPC_dlubcc0 = exploreMPC_dl_cc + 8;
exploreMPC_FLOAT* exploreMPC_dsubcc0 = exploreMPC_ds_cc + 8;
exploreMPC_FLOAT* exploreMPC_ccrhsub0 = exploreMPC_ccrhs + 8;
exploreMPC_FLOAT exploreMPC_Phi0[8];
exploreMPC_FLOAT exploreMPC_D0[8] = {1.0000000000000000E+000, 
1.0000000000000000E+000, 
1.0000000000000000E+000, 
1.0000000000000000E+000};
exploreMPC_FLOAT exploreMPC_W0[8];
exploreMPC_FLOAT* exploreMPC_z1 = exploreMPC_z + 8;
exploreMPC_FLOAT* exploreMPC_dzaff1 = exploreMPC_dz_aff + 8;
exploreMPC_FLOAT* exploreMPC_dzcc1 = exploreMPC_dz_cc + 8;
exploreMPC_FLOAT* exploreMPC_rd1 = exploreMPC_rd + 8;
exploreMPC_FLOAT exploreMPC_Lbyrd1[8];
exploreMPC_FLOAT* exploreMPC_grad_cost1 = exploreMPC_grad_cost + 8;
exploreMPC_FLOAT* exploreMPC_grad_eq1 = exploreMPC_grad_eq + 8;
exploreMPC_FLOAT* exploreMPC_grad_ineq1 = exploreMPC_grad_ineq + 8;
exploreMPC_FLOAT exploreMPC_ctv1[8];
exploreMPC_FLOAT* exploreMPC_v1 = exploreMPC_v + 4;
exploreMPC_FLOAT exploreMPC_re1[4];
exploreMPC_FLOAT exploreMPC_beta1[4];
exploreMPC_FLOAT exploreMPC_betacc1[4];
exploreMPC_FLOAT* exploreMPC_dvaff1 = exploreMPC_dv_aff + 4;
exploreMPC_FLOAT* exploreMPC_dvcc1 = exploreMPC_dv_cc + 4;
exploreMPC_FLOAT exploreMPC_V1[32];
exploreMPC_FLOAT exploreMPC_Yd1[10];
exploreMPC_FLOAT exploreMPC_Ld1[10];
exploreMPC_FLOAT exploreMPC_yy1[4];
exploreMPC_FLOAT exploreMPC_bmy1[4];
exploreMPC_FLOAT exploreMPC_c1[4] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
int exploreMPC_lbIdx1[8] = {0, 1, 2, 3, 4, 5, 6, 7};
exploreMPC_FLOAT* exploreMPC_llb1 = exploreMPC_l + 16;
exploreMPC_FLOAT* exploreMPC_slb1 = exploreMPC_s + 16;
exploreMPC_FLOAT* exploreMPC_llbbyslb1 = exploreMPC_lbys + 16;
exploreMPC_FLOAT exploreMPC_rilb1[8];
exploreMPC_FLOAT* exploreMPC_dllbaff1 = exploreMPC_dl_aff + 16;
exploreMPC_FLOAT* exploreMPC_dslbaff1 = exploreMPC_ds_aff + 16;
exploreMPC_FLOAT* exploreMPC_dllbcc1 = exploreMPC_dl_cc + 16;
exploreMPC_FLOAT* exploreMPC_dslbcc1 = exploreMPC_ds_cc + 16;
exploreMPC_FLOAT* exploreMPC_ccrhsl1 = exploreMPC_ccrhs + 16;
int exploreMPC_ubIdx1[8] = {0, 1, 2, 3, 4, 5, 6, 7};
exploreMPC_FLOAT* exploreMPC_lub1 = exploreMPC_l + 24;
exploreMPC_FLOAT* exploreMPC_sub1 = exploreMPC_s + 24;
exploreMPC_FLOAT* exploreMPC_lubbysub1 = exploreMPC_lbys + 24;
exploreMPC_FLOAT exploreMPC_riub1[8];
exploreMPC_FLOAT* exploreMPC_dlubaff1 = exploreMPC_dl_aff + 24;
exploreMPC_FLOAT* exploreMPC_dsubaff1 = exploreMPC_ds_aff + 24;
exploreMPC_FLOAT* exploreMPC_dlubcc1 = exploreMPC_dl_cc + 24;
exploreMPC_FLOAT* exploreMPC_dsubcc1 = exploreMPC_ds_cc + 24;
exploreMPC_FLOAT* exploreMPC_ccrhsub1 = exploreMPC_ccrhs + 24;
exploreMPC_FLOAT exploreMPC_Phi1[8];
exploreMPC_FLOAT exploreMPC_D1[8] = {-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000};
exploreMPC_FLOAT exploreMPC_W1[8];
exploreMPC_FLOAT exploreMPC_Ysd1[16];
exploreMPC_FLOAT exploreMPC_Lsd1[16];
exploreMPC_FLOAT* exploreMPC_z2 = exploreMPC_z + 16;
exploreMPC_FLOAT* exploreMPC_dzaff2 = exploreMPC_dz_aff + 16;
exploreMPC_FLOAT* exploreMPC_dzcc2 = exploreMPC_dz_cc + 16;
exploreMPC_FLOAT* exploreMPC_rd2 = exploreMPC_rd + 16;
exploreMPC_FLOAT exploreMPC_Lbyrd2[4];
exploreMPC_FLOAT* exploreMPC_grad_cost2 = exploreMPC_grad_cost + 16;
exploreMPC_FLOAT* exploreMPC_grad_eq2 = exploreMPC_grad_eq + 16;
exploreMPC_FLOAT* exploreMPC_grad_ineq2 = exploreMPC_grad_ineq + 16;
exploreMPC_FLOAT exploreMPC_ctv2[4];
exploreMPC_FLOAT* exploreMPC_v2 = exploreMPC_v + 8;
exploreMPC_FLOAT exploreMPC_re2[4];
exploreMPC_FLOAT exploreMPC_beta2[4];
exploreMPC_FLOAT exploreMPC_betacc2[4];
exploreMPC_FLOAT* exploreMPC_dvaff2 = exploreMPC_dv_aff + 8;
exploreMPC_FLOAT* exploreMPC_dvcc2 = exploreMPC_dv_cc + 8;
exploreMPC_FLOAT exploreMPC_V2[16];
exploreMPC_FLOAT exploreMPC_Yd2[10];
exploreMPC_FLOAT exploreMPC_Ld2[10];
exploreMPC_FLOAT exploreMPC_yy2[4];
exploreMPC_FLOAT exploreMPC_bmy2[4];
exploreMPC_FLOAT exploreMPC_c2[4] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
int exploreMPC_lbIdx2[4] = {0, 1, 2, 3};
exploreMPC_FLOAT* exploreMPC_llb2 = exploreMPC_l + 32;
exploreMPC_FLOAT* exploreMPC_slb2 = exploreMPC_s + 32;
exploreMPC_FLOAT* exploreMPC_llbbyslb2 = exploreMPC_lbys + 32;
exploreMPC_FLOAT exploreMPC_rilb2[4];
exploreMPC_FLOAT* exploreMPC_dllbaff2 = exploreMPC_dl_aff + 32;
exploreMPC_FLOAT* exploreMPC_dslbaff2 = exploreMPC_ds_aff + 32;
exploreMPC_FLOAT* exploreMPC_dllbcc2 = exploreMPC_dl_cc + 32;
exploreMPC_FLOAT* exploreMPC_dslbcc2 = exploreMPC_ds_cc + 32;
exploreMPC_FLOAT* exploreMPC_ccrhsl2 = exploreMPC_ccrhs + 32;
int exploreMPC_ubIdx2[4] = {0, 1, 2, 3};
exploreMPC_FLOAT* exploreMPC_lub2 = exploreMPC_l + 36;
exploreMPC_FLOAT* exploreMPC_sub2 = exploreMPC_s + 36;
exploreMPC_FLOAT* exploreMPC_lubbysub2 = exploreMPC_lbys + 36;
exploreMPC_FLOAT exploreMPC_riub2[4];
exploreMPC_FLOAT* exploreMPC_dlubaff2 = exploreMPC_dl_aff + 36;
exploreMPC_FLOAT* exploreMPC_dsubaff2 = exploreMPC_ds_aff + 36;
exploreMPC_FLOAT* exploreMPC_dlubcc2 = exploreMPC_dl_cc + 36;
exploreMPC_FLOAT* exploreMPC_dsubcc2 = exploreMPC_ds_cc + 36;
exploreMPC_FLOAT* exploreMPC_ccrhsub2 = exploreMPC_ccrhs + 36;
exploreMPC_FLOAT exploreMPC_Phi2[4];
exploreMPC_FLOAT exploreMPC_D2[4] = {-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000};
exploreMPC_FLOAT exploreMPC_W2[4];
exploreMPC_FLOAT exploreMPC_Ysd2[16];
exploreMPC_FLOAT exploreMPC_Lsd2[16];
exploreMPC_FLOAT musigma;
exploreMPC_FLOAT sigma_3rdroot;
exploreMPC_FLOAT exploreMPC_Diag1_0[8];
exploreMPC_FLOAT exploreMPC_Diag2_0[8];
exploreMPC_FLOAT exploreMPC_L_0[28];




/* SOLVER CODE --------------------------------------------------------- */
int exploreMPC_solve(exploreMPC_params* params, exploreMPC_output* output, exploreMPC_info* info)
{	
int exitcode;

#if exploreMPC_SET_TIMING == 1
	exploreMPC_timer solvertimer;
	exploreMPC_tic(&solvertimer);
#endif
/* FUNCTION CALLS INTO LA LIBRARY -------------------------------------- */
info->it = 0;
exploreMPC_LA_INITIALIZEVECTOR_20(exploreMPC_z, 0);
exploreMPC_LA_INITIALIZEVECTOR_12(exploreMPC_v, 1);
exploreMPC_LA_INITIALIZEVECTOR_40(exploreMPC_l, 1);
exploreMPC_LA_INITIALIZEVECTOR_40(exploreMPC_s, 1);
info->mu = 0;
exploreMPC_LA_DOTACC_40(exploreMPC_l, exploreMPC_s, &info->mu);
info->mu /= 40;
while( 1 ){
info->pobj = 0;
exploreMPC_LA_DIAG_QUADFCN_8(params->H1, params->f1, exploreMPC_z0, exploreMPC_grad_cost0, &info->pobj);
exploreMPC_LA_DIAG_QUADFCN_8(params->H2, params->f2, exploreMPC_z1, exploreMPC_grad_cost1, &info->pobj);
exploreMPC_LA_DIAG_QUADFCN_4(params->H3, params->f3, exploreMPC_z2, exploreMPC_grad_cost2, &info->pobj);
info->res_eq = 0;
info->dgap = 0;
exploreMPC_LA_DIAGZERO_MVMSUB6_4(exploreMPC_D0, exploreMPC_z0, params->c1, exploreMPC_v0, exploreMPC_re0, &info->dgap, &info->res_eq);
exploreMPC_LA_DENSE_DIAGZERO_MVMSUB3_4_8_8(exploreMPC_C0, exploreMPC_z0, exploreMPC_D1, exploreMPC_z1, exploreMPC_c1, exploreMPC_v1, exploreMPC_re1, &info->dgap, &info->res_eq);
exploreMPC_LA_DENSE_DIAGZERO_MVMSUB3_4_8_4(exploreMPC_C0, exploreMPC_z1, exploreMPC_D2, exploreMPC_z2, exploreMPC_c2, exploreMPC_v2, exploreMPC_re2, &info->dgap, &info->res_eq);
exploreMPC_LA_DENSE_DIAGZERO_MTVM2_4_8_4(exploreMPC_C0, exploreMPC_v1, exploreMPC_D0, exploreMPC_v0, exploreMPC_grad_eq0);
exploreMPC_LA_DENSE_DIAGZERO_MTVM2_4_8_4(exploreMPC_C0, exploreMPC_v2, exploreMPC_D1, exploreMPC_v1, exploreMPC_grad_eq1);
exploreMPC_LA_DIAGZERO_MTVM_4_4(exploreMPC_D2, exploreMPC_v2, exploreMPC_grad_eq2);
info->res_ineq = 0;
exploreMPC_LA_VSUBADD3_8(params->lb1, exploreMPC_z0, exploreMPC_lbIdx0, exploreMPC_llb0, exploreMPC_slb0, exploreMPC_rilb0, &info->dgap, &info->res_ineq);
exploreMPC_LA_VSUBADD2_8(exploreMPC_z0, exploreMPC_ubIdx0, params->ub1, exploreMPC_lub0, exploreMPC_sub0, exploreMPC_riub0, &info->dgap, &info->res_ineq);
exploreMPC_LA_VSUBADD3_8(params->lb2, exploreMPC_z1, exploreMPC_lbIdx1, exploreMPC_llb1, exploreMPC_slb1, exploreMPC_rilb1, &info->dgap, &info->res_ineq);
exploreMPC_LA_VSUBADD2_8(exploreMPC_z1, exploreMPC_ubIdx1, params->ub2, exploreMPC_lub1, exploreMPC_sub1, exploreMPC_riub1, &info->dgap, &info->res_ineq);
exploreMPC_LA_VSUBADD3_4(params->lb3, exploreMPC_z2, exploreMPC_lbIdx2, exploreMPC_llb2, exploreMPC_slb2, exploreMPC_rilb2, &info->dgap, &info->res_ineq);
exploreMPC_LA_VSUBADD2_4(exploreMPC_z2, exploreMPC_ubIdx2, params->ub3, exploreMPC_lub2, exploreMPC_sub2, exploreMPC_riub2, &info->dgap, &info->res_ineq);
exploreMPC_LA_INEQ_B_GRAD_8_8_8(exploreMPC_lub0, exploreMPC_sub0, exploreMPC_riub0, exploreMPC_llb0, exploreMPC_slb0, exploreMPC_rilb0, exploreMPC_lbIdx0, exploreMPC_ubIdx0, exploreMPC_grad_ineq0, exploreMPC_lubbysub0, exploreMPC_llbbyslb0);
exploreMPC_LA_INEQ_B_GRAD_8_8_8(exploreMPC_lub1, exploreMPC_sub1, exploreMPC_riub1, exploreMPC_llb1, exploreMPC_slb1, exploreMPC_rilb1, exploreMPC_lbIdx1, exploreMPC_ubIdx1, exploreMPC_grad_ineq1, exploreMPC_lubbysub1, exploreMPC_llbbyslb1);
exploreMPC_LA_INEQ_B_GRAD_4_4_4(exploreMPC_lub2, exploreMPC_sub2, exploreMPC_riub2, exploreMPC_llb2, exploreMPC_slb2, exploreMPC_rilb2, exploreMPC_lbIdx2, exploreMPC_ubIdx2, exploreMPC_grad_ineq2, exploreMPC_lubbysub2, exploreMPC_llbbyslb2);
info->dobj = info->pobj - info->dgap;
info->rdgap = info->pobj ? info->dgap / info->pobj : 1e6;
if( info->rdgap < 0 ) info->rdgap = -info->rdgap;
if( info->mu < exploreMPC_SET_ACC_KKTCOMPL
    && (info->rdgap < exploreMPC_SET_ACC_RDGAP || info->dgap < exploreMPC_SET_ACC_KKTCOMPL)
    && info->res_eq < exploreMPC_SET_ACC_RESEQ
    && info->res_ineq < exploreMPC_SET_ACC_RESINEQ ){
exitcode = exploreMPC_OPTIMAL; break; }
if( info->it == exploreMPC_SET_MAXIT ){
exitcode = exploreMPC_MAXITREACHED; break; }
exploreMPC_LA_VVADD3_20(exploreMPC_grad_cost, exploreMPC_grad_eq, exploreMPC_grad_ineq, exploreMPC_rd);
exploreMPC_LA_DIAG_CHOL_ONELOOP_LBUB_8_8_8(params->H1, exploreMPC_llbbyslb0, exploreMPC_lbIdx0, exploreMPC_lubbysub0, exploreMPC_ubIdx0, exploreMPC_Phi0);
exploreMPC_LA_DIAG_MATRIXFORWARDSUB_4_8(exploreMPC_Phi0, exploreMPC_C0, exploreMPC_V0);
exploreMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_4_8(exploreMPC_Phi0, exploreMPC_D0, exploreMPC_W0);
exploreMPC_LA_DENSE_DIAGZERO_MMTM_4_8_4(exploreMPC_W0, exploreMPC_V0, exploreMPC_Ysd1);
exploreMPC_LA_DIAG_FORWARDSUB_8(exploreMPC_Phi0, exploreMPC_rd0, exploreMPC_Lbyrd0);
exploreMPC_LA_DIAG_CHOL_ONELOOP_LBUB_8_8_8(params->H2, exploreMPC_llbbyslb1, exploreMPC_lbIdx1, exploreMPC_lubbysub1, exploreMPC_ubIdx1, exploreMPC_Phi1);
exploreMPC_LA_DIAG_MATRIXFORWARDSUB_4_8(exploreMPC_Phi1, exploreMPC_C0, exploreMPC_V1);
exploreMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_4_8(exploreMPC_Phi1, exploreMPC_D1, exploreMPC_W1);
exploreMPC_LA_DENSE_DIAGZERO_MMTM_4_8_4(exploreMPC_W1, exploreMPC_V1, exploreMPC_Ysd2);
exploreMPC_LA_DIAG_FORWARDSUB_8(exploreMPC_Phi1, exploreMPC_rd1, exploreMPC_Lbyrd1);
exploreMPC_LA_DIAG_CHOL_ONELOOP_LBUB_4_4_4(params->H3, exploreMPC_llbbyslb2, exploreMPC_lbIdx2, exploreMPC_lubbysub2, exploreMPC_ubIdx2, exploreMPC_Phi2);
exploreMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_4_4(exploreMPC_Phi2, exploreMPC_D2, exploreMPC_W2);
exploreMPC_LA_DIAG_FORWARDSUB_4(exploreMPC_Phi2, exploreMPC_rd2, exploreMPC_Lbyrd2);
exploreMPC_LA_DIAGZERO_MMT_4(exploreMPC_W0, exploreMPC_Yd0);
exploreMPC_LA_DIAGZERO_MVMSUB7_4(exploreMPC_W0, exploreMPC_Lbyrd0, exploreMPC_re0, exploreMPC_beta0);
exploreMPC_LA_DENSE_DIAGZERO_MMT2_4_8_8(exploreMPC_V0, exploreMPC_W1, exploreMPC_Yd1);
exploreMPC_LA_DENSE_DIAGZERO_2MVMSUB2_4_8_8(exploreMPC_V0, exploreMPC_Lbyrd0, exploreMPC_W1, exploreMPC_Lbyrd1, exploreMPC_re1, exploreMPC_beta1);
exploreMPC_LA_DENSE_DIAGZERO_MMT2_4_8_4(exploreMPC_V1, exploreMPC_W2, exploreMPC_Yd2);
exploreMPC_LA_DENSE_DIAGZERO_2MVMSUB2_4_8_4(exploreMPC_V1, exploreMPC_Lbyrd1, exploreMPC_W2, exploreMPC_Lbyrd2, exploreMPC_re2, exploreMPC_beta2);
exploreMPC_LA_DENSE_CHOL_4(exploreMPC_Yd0, exploreMPC_Ld0);
exploreMPC_LA_DENSE_FORWARDSUB_4(exploreMPC_Ld0, exploreMPC_beta0, exploreMPC_yy0);
exploreMPC_LA_DENSE_MATRIXTFORWARDSUB_4_4(exploreMPC_Ld0, exploreMPC_Ysd1, exploreMPC_Lsd1);
exploreMPC_LA_DENSE_MMTSUB_4_4(exploreMPC_Lsd1, exploreMPC_Yd1);
exploreMPC_LA_DENSE_CHOL_4(exploreMPC_Yd1, exploreMPC_Ld1);
exploreMPC_LA_DENSE_MVMSUB1_4_4(exploreMPC_Lsd1, exploreMPC_yy0, exploreMPC_beta1, exploreMPC_bmy1);
exploreMPC_LA_DENSE_FORWARDSUB_4(exploreMPC_Ld1, exploreMPC_bmy1, exploreMPC_yy1);
exploreMPC_LA_DENSE_MATRIXTFORWARDSUB_4_4(exploreMPC_Ld1, exploreMPC_Ysd2, exploreMPC_Lsd2);
exploreMPC_LA_DENSE_MMTSUB_4_4(exploreMPC_Lsd2, exploreMPC_Yd2);
exploreMPC_LA_DENSE_CHOL_4(exploreMPC_Yd2, exploreMPC_Ld2);
exploreMPC_LA_DENSE_MVMSUB1_4_4(exploreMPC_Lsd2, exploreMPC_yy1, exploreMPC_beta2, exploreMPC_bmy2);
exploreMPC_LA_DENSE_FORWARDSUB_4(exploreMPC_Ld2, exploreMPC_bmy2, exploreMPC_yy2);
exploreMPC_LA_DENSE_BACKWARDSUB_4(exploreMPC_Ld2, exploreMPC_yy2, exploreMPC_dvaff2);
exploreMPC_LA_DENSE_MTVMSUB_4_4(exploreMPC_Lsd2, exploreMPC_dvaff2, exploreMPC_yy1, exploreMPC_bmy1);
exploreMPC_LA_DENSE_BACKWARDSUB_4(exploreMPC_Ld1, exploreMPC_bmy1, exploreMPC_dvaff1);
exploreMPC_LA_DENSE_MTVMSUB_4_4(exploreMPC_Lsd1, exploreMPC_dvaff1, exploreMPC_yy0, exploreMPC_bmy0);
exploreMPC_LA_DENSE_BACKWARDSUB_4(exploreMPC_Ld0, exploreMPC_bmy0, exploreMPC_dvaff0);
exploreMPC_LA_DENSE_DIAGZERO_MTVM2_4_8_4(exploreMPC_C0, exploreMPC_dvaff1, exploreMPC_D0, exploreMPC_dvaff0, exploreMPC_grad_eq0);
exploreMPC_LA_DENSE_DIAGZERO_MTVM2_4_8_4(exploreMPC_C0, exploreMPC_dvaff2, exploreMPC_D1, exploreMPC_dvaff1, exploreMPC_grad_eq1);
exploreMPC_LA_DIAGZERO_MTVM_4_4(exploreMPC_D2, exploreMPC_dvaff2, exploreMPC_grad_eq2);
exploreMPC_LA_VSUB2_20(exploreMPC_rd, exploreMPC_grad_eq, exploreMPC_rd);
exploreMPC_LA_DIAG_FORWARDBACKWARDSUB_8(exploreMPC_Phi0, exploreMPC_rd0, exploreMPC_dzaff0);
exploreMPC_LA_DIAG_FORWARDBACKWARDSUB_8(exploreMPC_Phi1, exploreMPC_rd1, exploreMPC_dzaff1);
exploreMPC_LA_DIAG_FORWARDBACKWARDSUB_4(exploreMPC_Phi2, exploreMPC_rd2, exploreMPC_dzaff2);
exploreMPC_LA_VSUB_INDEXED_8(exploreMPC_dzaff0, exploreMPC_lbIdx0, exploreMPC_rilb0, exploreMPC_dslbaff0);
exploreMPC_LA_VSUB3_8(exploreMPC_llbbyslb0, exploreMPC_dslbaff0, exploreMPC_llb0, exploreMPC_dllbaff0);
exploreMPC_LA_VSUB2_INDEXED_8(exploreMPC_riub0, exploreMPC_dzaff0, exploreMPC_ubIdx0, exploreMPC_dsubaff0);
exploreMPC_LA_VSUB3_8(exploreMPC_lubbysub0, exploreMPC_dsubaff0, exploreMPC_lub0, exploreMPC_dlubaff0);
exploreMPC_LA_VSUB_INDEXED_8(exploreMPC_dzaff1, exploreMPC_lbIdx1, exploreMPC_rilb1, exploreMPC_dslbaff1);
exploreMPC_LA_VSUB3_8(exploreMPC_llbbyslb1, exploreMPC_dslbaff1, exploreMPC_llb1, exploreMPC_dllbaff1);
exploreMPC_LA_VSUB2_INDEXED_8(exploreMPC_riub1, exploreMPC_dzaff1, exploreMPC_ubIdx1, exploreMPC_dsubaff1);
exploreMPC_LA_VSUB3_8(exploreMPC_lubbysub1, exploreMPC_dsubaff1, exploreMPC_lub1, exploreMPC_dlubaff1);
exploreMPC_LA_VSUB_INDEXED_4(exploreMPC_dzaff2, exploreMPC_lbIdx2, exploreMPC_rilb2, exploreMPC_dslbaff2);
exploreMPC_LA_VSUB3_4(exploreMPC_llbbyslb2, exploreMPC_dslbaff2, exploreMPC_llb2, exploreMPC_dllbaff2);
exploreMPC_LA_VSUB2_INDEXED_4(exploreMPC_riub2, exploreMPC_dzaff2, exploreMPC_ubIdx2, exploreMPC_dsubaff2);
exploreMPC_LA_VSUB3_4(exploreMPC_lubbysub2, exploreMPC_dsubaff2, exploreMPC_lub2, exploreMPC_dlubaff2);
info->lsit_aff = exploreMPC_LINESEARCH_BACKTRACKING_AFFINE(exploreMPC_l, exploreMPC_s, exploreMPC_dl_aff, exploreMPC_ds_aff, &info->step_aff, &info->mu_aff);
if( info->lsit_aff == exploreMPC_NOPROGRESS ){
exitcode = exploreMPC_NOPROGRESS; break;
}
sigma_3rdroot = info->mu_aff / info->mu;
info->sigma = sigma_3rdroot*sigma_3rdroot*sigma_3rdroot;
musigma = info->mu * info->sigma;
exploreMPC_LA_VSUB5_40(exploreMPC_ds_aff, exploreMPC_dl_aff, musigma, exploreMPC_ccrhs);
exploreMPC_LA_VSUB6_INDEXED_8_8_8(exploreMPC_ccrhsub0, exploreMPC_sub0, exploreMPC_ubIdx0, exploreMPC_ccrhsl0, exploreMPC_slb0, exploreMPC_lbIdx0, exploreMPC_rd0);
exploreMPC_LA_VSUB6_INDEXED_8_8_8(exploreMPC_ccrhsub1, exploreMPC_sub1, exploreMPC_ubIdx1, exploreMPC_ccrhsl1, exploreMPC_slb1, exploreMPC_lbIdx1, exploreMPC_rd1);
exploreMPC_LA_DIAG_FORWARDSUB_8(exploreMPC_Phi0, exploreMPC_rd0, exploreMPC_Lbyrd0);
exploreMPC_LA_DIAG_FORWARDSUB_8(exploreMPC_Phi1, exploreMPC_rd1, exploreMPC_Lbyrd1);
exploreMPC_LA_DIAGZERO_MVM_4(exploreMPC_W0, exploreMPC_Lbyrd0, exploreMPC_beta0);
exploreMPC_LA_DENSE_FORWARDSUB_4(exploreMPC_Ld0, exploreMPC_beta0, exploreMPC_yy0);
exploreMPC_LA_DENSE_DIAGZERO_2MVMADD_4_8_8(exploreMPC_V0, exploreMPC_Lbyrd0, exploreMPC_W1, exploreMPC_Lbyrd1, exploreMPC_beta1);
exploreMPC_LA_DENSE_MVMSUB1_4_4(exploreMPC_Lsd1, exploreMPC_yy0, exploreMPC_beta1, exploreMPC_bmy1);
exploreMPC_LA_DENSE_FORWARDSUB_4(exploreMPC_Ld1, exploreMPC_bmy1, exploreMPC_yy1);
exploreMPC_LA_VSUB6_INDEXED_4_4_4(exploreMPC_ccrhsub2, exploreMPC_sub2, exploreMPC_ubIdx2, exploreMPC_ccrhsl2, exploreMPC_slb2, exploreMPC_lbIdx2, exploreMPC_rd2);
exploreMPC_LA_DIAG_FORWARDSUB_4(exploreMPC_Phi2, exploreMPC_rd2, exploreMPC_Lbyrd2);
exploreMPC_LA_DENSE_DIAGZERO_2MVMADD_4_8_4(exploreMPC_V1, exploreMPC_Lbyrd1, exploreMPC_W2, exploreMPC_Lbyrd2, exploreMPC_beta2);
exploreMPC_LA_DENSE_MVMSUB1_4_4(exploreMPC_Lsd2, exploreMPC_yy1, exploreMPC_beta2, exploreMPC_bmy2);
exploreMPC_LA_DENSE_FORWARDSUB_4(exploreMPC_Ld2, exploreMPC_bmy2, exploreMPC_yy2);
exploreMPC_LA_DENSE_BACKWARDSUB_4(exploreMPC_Ld2, exploreMPC_yy2, exploreMPC_dvcc2);
exploreMPC_LA_DENSE_MTVMSUB_4_4(exploreMPC_Lsd2, exploreMPC_dvcc2, exploreMPC_yy1, exploreMPC_bmy1);
exploreMPC_LA_DENSE_BACKWARDSUB_4(exploreMPC_Ld1, exploreMPC_bmy1, exploreMPC_dvcc1);
exploreMPC_LA_DENSE_MTVMSUB_4_4(exploreMPC_Lsd1, exploreMPC_dvcc1, exploreMPC_yy0, exploreMPC_bmy0);
exploreMPC_LA_DENSE_BACKWARDSUB_4(exploreMPC_Ld0, exploreMPC_bmy0, exploreMPC_dvcc0);
exploreMPC_LA_DENSE_DIAGZERO_MTVM2_4_8_4(exploreMPC_C0, exploreMPC_dvcc1, exploreMPC_D0, exploreMPC_dvcc0, exploreMPC_grad_eq0);
exploreMPC_LA_DENSE_DIAGZERO_MTVM2_4_8_4(exploreMPC_C0, exploreMPC_dvcc2, exploreMPC_D1, exploreMPC_dvcc1, exploreMPC_grad_eq1);
exploreMPC_LA_DIAGZERO_MTVM_4_4(exploreMPC_D2, exploreMPC_dvcc2, exploreMPC_grad_eq2);
exploreMPC_LA_VSUB_20(exploreMPC_rd, exploreMPC_grad_eq, exploreMPC_rd);
exploreMPC_LA_DIAG_FORWARDBACKWARDSUB_8(exploreMPC_Phi0, exploreMPC_rd0, exploreMPC_dzcc0);
exploreMPC_LA_DIAG_FORWARDBACKWARDSUB_8(exploreMPC_Phi1, exploreMPC_rd1, exploreMPC_dzcc1);
exploreMPC_LA_DIAG_FORWARDBACKWARDSUB_4(exploreMPC_Phi2, exploreMPC_rd2, exploreMPC_dzcc2);
exploreMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_8(exploreMPC_ccrhsl0, exploreMPC_slb0, exploreMPC_llbbyslb0, exploreMPC_dzcc0, exploreMPC_lbIdx0, exploreMPC_dllbcc0);
exploreMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_8(exploreMPC_ccrhsub0, exploreMPC_sub0, exploreMPC_lubbysub0, exploreMPC_dzcc0, exploreMPC_ubIdx0, exploreMPC_dlubcc0);
exploreMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_8(exploreMPC_ccrhsl1, exploreMPC_slb1, exploreMPC_llbbyslb1, exploreMPC_dzcc1, exploreMPC_lbIdx1, exploreMPC_dllbcc1);
exploreMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_8(exploreMPC_ccrhsub1, exploreMPC_sub1, exploreMPC_lubbysub1, exploreMPC_dzcc1, exploreMPC_ubIdx1, exploreMPC_dlubcc1);
exploreMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_4(exploreMPC_ccrhsl2, exploreMPC_slb2, exploreMPC_llbbyslb2, exploreMPC_dzcc2, exploreMPC_lbIdx2, exploreMPC_dllbcc2);
exploreMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_4(exploreMPC_ccrhsub2, exploreMPC_sub2, exploreMPC_lubbysub2, exploreMPC_dzcc2, exploreMPC_ubIdx2, exploreMPC_dlubcc2);
exploreMPC_LA_VSUB7_40(exploreMPC_l, exploreMPC_ccrhs, exploreMPC_s, exploreMPC_dl_cc, exploreMPC_ds_cc);
exploreMPC_LA_VADD_20(exploreMPC_dz_cc, exploreMPC_dz_aff);
exploreMPC_LA_VADD_12(exploreMPC_dv_cc, exploreMPC_dv_aff);
exploreMPC_LA_VADD_40(exploreMPC_dl_cc, exploreMPC_dl_aff);
exploreMPC_LA_VADD_40(exploreMPC_ds_cc, exploreMPC_ds_aff);
info->lsit_cc = exploreMPC_LINESEARCH_BACKTRACKING_COMBINED(exploreMPC_z, exploreMPC_v, exploreMPC_l, exploreMPC_s, exploreMPC_dz_cc, exploreMPC_dv_cc, exploreMPC_dl_cc, exploreMPC_ds_cc, &info->step_cc, &info->mu);
if( info->lsit_cc == exploreMPC_NOPROGRESS ){
exitcode = exploreMPC_NOPROGRESS; break;
}
info->it++;
}
output->z1[0] = exploreMPC_z0[0];
output->z1[1] = exploreMPC_z0[1];
output->z1[2] = exploreMPC_z0[2];
output->z1[3] = exploreMPC_z0[3];
output->z1[4] = exploreMPC_z0[4];
output->z1[5] = exploreMPC_z0[5];
output->z1[6] = exploreMPC_z0[6];
output->z1[7] = exploreMPC_z0[7];
output->z2[0] = exploreMPC_z1[0];
output->z2[1] = exploreMPC_z1[1];
output->z2[2] = exploreMPC_z1[2];
output->z2[3] = exploreMPC_z1[3];
output->z2[4] = exploreMPC_z1[4];
output->z2[5] = exploreMPC_z1[5];
output->z2[6] = exploreMPC_z1[6];
output->z2[7] = exploreMPC_z1[7];
output->z3[0] = exploreMPC_z2[0];
output->z3[1] = exploreMPC_z2[1];
output->z3[2] = exploreMPC_z2[2];
output->z3[3] = exploreMPC_z2[3];

#if exploreMPC_SET_TIMING == 1
info->solvetime = exploreMPC_toc(&solvertimer);
#if exploreMPC_SET_PRINTLEVEL > 0 && exploreMPC_SET_TIMING == 1
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
