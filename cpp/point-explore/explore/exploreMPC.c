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
 * Initializes a vector of length 76 with a value.
 */
void exploreMPC_LA_INITIALIZEVECTOR_76(exploreMPC_FLOAT* vec, exploreMPC_FLOAT value)
{
	int i;
	for( i=0; i<76; i++ )
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
 * Initializes a vector of length 152 with a value.
 */
void exploreMPC_LA_INITIALIZEVECTOR_152(exploreMPC_FLOAT* vec, exploreMPC_FLOAT value)
{
	int i;
	for( i=0; i<152; i++ )
	{
		vec[i] = value;
	}
}


/* 
 * Calculates a dot product and adds it to a variable: z += x'*y; 
 * This function is for vectors of length 152.
 */
void exploreMPC_LA_DOTACC_152(exploreMPC_FLOAT *x, exploreMPC_FLOAT *y, exploreMPC_FLOAT *z)
{
	int i;
	for( i=0; i<152; i++ ){
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
 * of length 76.
 */
void exploreMPC_LA_VVADD3_76(exploreMPC_FLOAT *u, exploreMPC_FLOAT *v, exploreMPC_FLOAT *w, exploreMPC_FLOAT *z)
{
	int i;
	for( i=0; i<76; i++ ){
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
 * Vector subtraction z = -x - y for vectors of length 76.
 */
void exploreMPC_LA_VSUB2_76(exploreMPC_FLOAT *x, exploreMPC_FLOAT *y, exploreMPC_FLOAT *z)
{
	int i;
	for( i=0; i<76; i++){
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
        for( i=0; i<152; i++ ){
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
        if( i == 152 ){
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
    *mu_aff = mymu / (exploreMPC_FLOAT)152;
    return lsIt;
}


/*
 * Vector subtraction x = u.*v - a where a is a scalar
*  and x,u,v are vectors of length 152.
 */
void exploreMPC_LA_VSUB5_152(exploreMPC_FLOAT *u, exploreMPC_FLOAT *v, exploreMPC_FLOAT a, exploreMPC_FLOAT *x)
{
	int i;
	for( i=0; i<152; i++){
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
 * Vector subtraction z = x - y for vectors of length 76.
 */
void exploreMPC_LA_VSUB_76(exploreMPC_FLOAT *x, exploreMPC_FLOAT *y, exploreMPC_FLOAT *z)
{
	int i;
	for( i=0; i<76; i++){
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
 * Computes ds = -l.\(r + s.*dl) for vectors of length 152.
 */
void exploreMPC_LA_VSUB7_152(exploreMPC_FLOAT *l, exploreMPC_FLOAT *r, exploreMPC_FLOAT *s, exploreMPC_FLOAT *dl, exploreMPC_FLOAT *ds)
{
	int i;
	for( i=0; i<152; i++){
		ds[i] = -(r[i] + s[i]*dl[i])/l[i];
	}
}


/*
 * Vector addition x = x + y for vectors of length 76.
 */
void exploreMPC_LA_VADD_76(exploreMPC_FLOAT *x, exploreMPC_FLOAT *y)
{
	int i;
	for( i=0; i<76; i++){
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


/*
 * Vector addition x = x + y for vectors of length 152.
 */
void exploreMPC_LA_VADD_152(exploreMPC_FLOAT *x, exploreMPC_FLOAT *y)
{
	int i;
	for( i=0; i<152; i++){
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
        for( i=0; i<152; i++ ){
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
        if( i == 152 ){
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
    for( i=0; i<76; i++ ){
        z[i] += a_gamma*dz[i];
    }
    
    /* equality constraint multipliers */
    for( i=0; i<40; i++ ){
        v[i] += a_gamma*dv[i];
    }
    
    /* inequality constraint multipliers & slacks, also update mu */
    *mu = 0;
    for( i=0; i<152; i++ ){
        dltemp = l[i] + a_gamma*dl[i]; l[i] = dltemp;
        dstemp = s[i] + a_gamma*ds[i]; s[i] = dstemp;
        *mu += dltemp*dstemp;
    }
    
    *a = a_gamma;
    *mu /= (exploreMPC_FLOAT)152;
    return lsIt;
}




/* VARIABLE DEFINITIONS ------------------------------------------------ */
exploreMPC_FLOAT exploreMPC_z[76];
exploreMPC_FLOAT exploreMPC_v[40];
exploreMPC_FLOAT exploreMPC_dz_aff[76];
exploreMPC_FLOAT exploreMPC_dv_aff[40];
exploreMPC_FLOAT exploreMPC_grad_cost[76];
exploreMPC_FLOAT exploreMPC_grad_eq[76];
exploreMPC_FLOAT exploreMPC_rd[76];
exploreMPC_FLOAT exploreMPC_l[152];
exploreMPC_FLOAT exploreMPC_s[152];
exploreMPC_FLOAT exploreMPC_lbys[152];
exploreMPC_FLOAT exploreMPC_dl_aff[152];
exploreMPC_FLOAT exploreMPC_ds_aff[152];
exploreMPC_FLOAT exploreMPC_dz_cc[76];
exploreMPC_FLOAT exploreMPC_dv_cc[40];
exploreMPC_FLOAT exploreMPC_dl_cc[152];
exploreMPC_FLOAT exploreMPC_ds_cc[152];
exploreMPC_FLOAT exploreMPC_ccrhs[152];
exploreMPC_FLOAT exploreMPC_grad_ineq[76];
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
exploreMPC_FLOAT exploreMPC_Lbyrd2[8];
exploreMPC_FLOAT* exploreMPC_grad_cost2 = exploreMPC_grad_cost + 16;
exploreMPC_FLOAT* exploreMPC_grad_eq2 = exploreMPC_grad_eq + 16;
exploreMPC_FLOAT* exploreMPC_grad_ineq2 = exploreMPC_grad_ineq + 16;
exploreMPC_FLOAT exploreMPC_ctv2[8];
exploreMPC_FLOAT* exploreMPC_v2 = exploreMPC_v + 8;
exploreMPC_FLOAT exploreMPC_re2[4];
exploreMPC_FLOAT exploreMPC_beta2[4];
exploreMPC_FLOAT exploreMPC_betacc2[4];
exploreMPC_FLOAT* exploreMPC_dvaff2 = exploreMPC_dv_aff + 8;
exploreMPC_FLOAT* exploreMPC_dvcc2 = exploreMPC_dv_cc + 8;
exploreMPC_FLOAT exploreMPC_V2[32];
exploreMPC_FLOAT exploreMPC_Yd2[10];
exploreMPC_FLOAT exploreMPC_Ld2[10];
exploreMPC_FLOAT exploreMPC_yy2[4];
exploreMPC_FLOAT exploreMPC_bmy2[4];
exploreMPC_FLOAT exploreMPC_c2[4] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
int exploreMPC_lbIdx2[8] = {0, 1, 2, 3, 4, 5, 6, 7};
exploreMPC_FLOAT* exploreMPC_llb2 = exploreMPC_l + 32;
exploreMPC_FLOAT* exploreMPC_slb2 = exploreMPC_s + 32;
exploreMPC_FLOAT* exploreMPC_llbbyslb2 = exploreMPC_lbys + 32;
exploreMPC_FLOAT exploreMPC_rilb2[8];
exploreMPC_FLOAT* exploreMPC_dllbaff2 = exploreMPC_dl_aff + 32;
exploreMPC_FLOAT* exploreMPC_dslbaff2 = exploreMPC_ds_aff + 32;
exploreMPC_FLOAT* exploreMPC_dllbcc2 = exploreMPC_dl_cc + 32;
exploreMPC_FLOAT* exploreMPC_dslbcc2 = exploreMPC_ds_cc + 32;
exploreMPC_FLOAT* exploreMPC_ccrhsl2 = exploreMPC_ccrhs + 32;
int exploreMPC_ubIdx2[8] = {0, 1, 2, 3, 4, 5, 6, 7};
exploreMPC_FLOAT* exploreMPC_lub2 = exploreMPC_l + 40;
exploreMPC_FLOAT* exploreMPC_sub2 = exploreMPC_s + 40;
exploreMPC_FLOAT* exploreMPC_lubbysub2 = exploreMPC_lbys + 40;
exploreMPC_FLOAT exploreMPC_riub2[8];
exploreMPC_FLOAT* exploreMPC_dlubaff2 = exploreMPC_dl_aff + 40;
exploreMPC_FLOAT* exploreMPC_dsubaff2 = exploreMPC_ds_aff + 40;
exploreMPC_FLOAT* exploreMPC_dlubcc2 = exploreMPC_dl_cc + 40;
exploreMPC_FLOAT* exploreMPC_dsubcc2 = exploreMPC_ds_cc + 40;
exploreMPC_FLOAT* exploreMPC_ccrhsub2 = exploreMPC_ccrhs + 40;
exploreMPC_FLOAT exploreMPC_Phi2[8];
exploreMPC_FLOAT exploreMPC_W2[8];
exploreMPC_FLOAT exploreMPC_Ysd2[16];
exploreMPC_FLOAT exploreMPC_Lsd2[16];
exploreMPC_FLOAT* exploreMPC_z3 = exploreMPC_z + 24;
exploreMPC_FLOAT* exploreMPC_dzaff3 = exploreMPC_dz_aff + 24;
exploreMPC_FLOAT* exploreMPC_dzcc3 = exploreMPC_dz_cc + 24;
exploreMPC_FLOAT* exploreMPC_rd3 = exploreMPC_rd + 24;
exploreMPC_FLOAT exploreMPC_Lbyrd3[8];
exploreMPC_FLOAT* exploreMPC_grad_cost3 = exploreMPC_grad_cost + 24;
exploreMPC_FLOAT* exploreMPC_grad_eq3 = exploreMPC_grad_eq + 24;
exploreMPC_FLOAT* exploreMPC_grad_ineq3 = exploreMPC_grad_ineq + 24;
exploreMPC_FLOAT exploreMPC_ctv3[8];
exploreMPC_FLOAT* exploreMPC_v3 = exploreMPC_v + 12;
exploreMPC_FLOAT exploreMPC_re3[4];
exploreMPC_FLOAT exploreMPC_beta3[4];
exploreMPC_FLOAT exploreMPC_betacc3[4];
exploreMPC_FLOAT* exploreMPC_dvaff3 = exploreMPC_dv_aff + 12;
exploreMPC_FLOAT* exploreMPC_dvcc3 = exploreMPC_dv_cc + 12;
exploreMPC_FLOAT exploreMPC_V3[32];
exploreMPC_FLOAT exploreMPC_Yd3[10];
exploreMPC_FLOAT exploreMPC_Ld3[10];
exploreMPC_FLOAT exploreMPC_yy3[4];
exploreMPC_FLOAT exploreMPC_bmy3[4];
exploreMPC_FLOAT exploreMPC_c3[4] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
int exploreMPC_lbIdx3[8] = {0, 1, 2, 3, 4, 5, 6, 7};
exploreMPC_FLOAT* exploreMPC_llb3 = exploreMPC_l + 48;
exploreMPC_FLOAT* exploreMPC_slb3 = exploreMPC_s + 48;
exploreMPC_FLOAT* exploreMPC_llbbyslb3 = exploreMPC_lbys + 48;
exploreMPC_FLOAT exploreMPC_rilb3[8];
exploreMPC_FLOAT* exploreMPC_dllbaff3 = exploreMPC_dl_aff + 48;
exploreMPC_FLOAT* exploreMPC_dslbaff3 = exploreMPC_ds_aff + 48;
exploreMPC_FLOAT* exploreMPC_dllbcc3 = exploreMPC_dl_cc + 48;
exploreMPC_FLOAT* exploreMPC_dslbcc3 = exploreMPC_ds_cc + 48;
exploreMPC_FLOAT* exploreMPC_ccrhsl3 = exploreMPC_ccrhs + 48;
int exploreMPC_ubIdx3[8] = {0, 1, 2, 3, 4, 5, 6, 7};
exploreMPC_FLOAT* exploreMPC_lub3 = exploreMPC_l + 56;
exploreMPC_FLOAT* exploreMPC_sub3 = exploreMPC_s + 56;
exploreMPC_FLOAT* exploreMPC_lubbysub3 = exploreMPC_lbys + 56;
exploreMPC_FLOAT exploreMPC_riub3[8];
exploreMPC_FLOAT* exploreMPC_dlubaff3 = exploreMPC_dl_aff + 56;
exploreMPC_FLOAT* exploreMPC_dsubaff3 = exploreMPC_ds_aff + 56;
exploreMPC_FLOAT* exploreMPC_dlubcc3 = exploreMPC_dl_cc + 56;
exploreMPC_FLOAT* exploreMPC_dsubcc3 = exploreMPC_ds_cc + 56;
exploreMPC_FLOAT* exploreMPC_ccrhsub3 = exploreMPC_ccrhs + 56;
exploreMPC_FLOAT exploreMPC_Phi3[8];
exploreMPC_FLOAT exploreMPC_W3[8];
exploreMPC_FLOAT exploreMPC_Ysd3[16];
exploreMPC_FLOAT exploreMPC_Lsd3[16];
exploreMPC_FLOAT* exploreMPC_z4 = exploreMPC_z + 32;
exploreMPC_FLOAT* exploreMPC_dzaff4 = exploreMPC_dz_aff + 32;
exploreMPC_FLOAT* exploreMPC_dzcc4 = exploreMPC_dz_cc + 32;
exploreMPC_FLOAT* exploreMPC_rd4 = exploreMPC_rd + 32;
exploreMPC_FLOAT exploreMPC_Lbyrd4[8];
exploreMPC_FLOAT* exploreMPC_grad_cost4 = exploreMPC_grad_cost + 32;
exploreMPC_FLOAT* exploreMPC_grad_eq4 = exploreMPC_grad_eq + 32;
exploreMPC_FLOAT* exploreMPC_grad_ineq4 = exploreMPC_grad_ineq + 32;
exploreMPC_FLOAT exploreMPC_ctv4[8];
exploreMPC_FLOAT* exploreMPC_v4 = exploreMPC_v + 16;
exploreMPC_FLOAT exploreMPC_re4[4];
exploreMPC_FLOAT exploreMPC_beta4[4];
exploreMPC_FLOAT exploreMPC_betacc4[4];
exploreMPC_FLOAT* exploreMPC_dvaff4 = exploreMPC_dv_aff + 16;
exploreMPC_FLOAT* exploreMPC_dvcc4 = exploreMPC_dv_cc + 16;
exploreMPC_FLOAT exploreMPC_V4[32];
exploreMPC_FLOAT exploreMPC_Yd4[10];
exploreMPC_FLOAT exploreMPC_Ld4[10];
exploreMPC_FLOAT exploreMPC_yy4[4];
exploreMPC_FLOAT exploreMPC_bmy4[4];
exploreMPC_FLOAT exploreMPC_c4[4] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
int exploreMPC_lbIdx4[8] = {0, 1, 2, 3, 4, 5, 6, 7};
exploreMPC_FLOAT* exploreMPC_llb4 = exploreMPC_l + 64;
exploreMPC_FLOAT* exploreMPC_slb4 = exploreMPC_s + 64;
exploreMPC_FLOAT* exploreMPC_llbbyslb4 = exploreMPC_lbys + 64;
exploreMPC_FLOAT exploreMPC_rilb4[8];
exploreMPC_FLOAT* exploreMPC_dllbaff4 = exploreMPC_dl_aff + 64;
exploreMPC_FLOAT* exploreMPC_dslbaff4 = exploreMPC_ds_aff + 64;
exploreMPC_FLOAT* exploreMPC_dllbcc4 = exploreMPC_dl_cc + 64;
exploreMPC_FLOAT* exploreMPC_dslbcc4 = exploreMPC_ds_cc + 64;
exploreMPC_FLOAT* exploreMPC_ccrhsl4 = exploreMPC_ccrhs + 64;
int exploreMPC_ubIdx4[8] = {0, 1, 2, 3, 4, 5, 6, 7};
exploreMPC_FLOAT* exploreMPC_lub4 = exploreMPC_l + 72;
exploreMPC_FLOAT* exploreMPC_sub4 = exploreMPC_s + 72;
exploreMPC_FLOAT* exploreMPC_lubbysub4 = exploreMPC_lbys + 72;
exploreMPC_FLOAT exploreMPC_riub4[8];
exploreMPC_FLOAT* exploreMPC_dlubaff4 = exploreMPC_dl_aff + 72;
exploreMPC_FLOAT* exploreMPC_dsubaff4 = exploreMPC_ds_aff + 72;
exploreMPC_FLOAT* exploreMPC_dlubcc4 = exploreMPC_dl_cc + 72;
exploreMPC_FLOAT* exploreMPC_dsubcc4 = exploreMPC_ds_cc + 72;
exploreMPC_FLOAT* exploreMPC_ccrhsub4 = exploreMPC_ccrhs + 72;
exploreMPC_FLOAT exploreMPC_Phi4[8];
exploreMPC_FLOAT exploreMPC_W4[8];
exploreMPC_FLOAT exploreMPC_Ysd4[16];
exploreMPC_FLOAT exploreMPC_Lsd4[16];
exploreMPC_FLOAT* exploreMPC_z5 = exploreMPC_z + 40;
exploreMPC_FLOAT* exploreMPC_dzaff5 = exploreMPC_dz_aff + 40;
exploreMPC_FLOAT* exploreMPC_dzcc5 = exploreMPC_dz_cc + 40;
exploreMPC_FLOAT* exploreMPC_rd5 = exploreMPC_rd + 40;
exploreMPC_FLOAT exploreMPC_Lbyrd5[8];
exploreMPC_FLOAT* exploreMPC_grad_cost5 = exploreMPC_grad_cost + 40;
exploreMPC_FLOAT* exploreMPC_grad_eq5 = exploreMPC_grad_eq + 40;
exploreMPC_FLOAT* exploreMPC_grad_ineq5 = exploreMPC_grad_ineq + 40;
exploreMPC_FLOAT exploreMPC_ctv5[8];
exploreMPC_FLOAT* exploreMPC_v5 = exploreMPC_v + 20;
exploreMPC_FLOAT exploreMPC_re5[4];
exploreMPC_FLOAT exploreMPC_beta5[4];
exploreMPC_FLOAT exploreMPC_betacc5[4];
exploreMPC_FLOAT* exploreMPC_dvaff5 = exploreMPC_dv_aff + 20;
exploreMPC_FLOAT* exploreMPC_dvcc5 = exploreMPC_dv_cc + 20;
exploreMPC_FLOAT exploreMPC_V5[32];
exploreMPC_FLOAT exploreMPC_Yd5[10];
exploreMPC_FLOAT exploreMPC_Ld5[10];
exploreMPC_FLOAT exploreMPC_yy5[4];
exploreMPC_FLOAT exploreMPC_bmy5[4];
exploreMPC_FLOAT exploreMPC_c5[4] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
int exploreMPC_lbIdx5[8] = {0, 1, 2, 3, 4, 5, 6, 7};
exploreMPC_FLOAT* exploreMPC_llb5 = exploreMPC_l + 80;
exploreMPC_FLOAT* exploreMPC_slb5 = exploreMPC_s + 80;
exploreMPC_FLOAT* exploreMPC_llbbyslb5 = exploreMPC_lbys + 80;
exploreMPC_FLOAT exploreMPC_rilb5[8];
exploreMPC_FLOAT* exploreMPC_dllbaff5 = exploreMPC_dl_aff + 80;
exploreMPC_FLOAT* exploreMPC_dslbaff5 = exploreMPC_ds_aff + 80;
exploreMPC_FLOAT* exploreMPC_dllbcc5 = exploreMPC_dl_cc + 80;
exploreMPC_FLOAT* exploreMPC_dslbcc5 = exploreMPC_ds_cc + 80;
exploreMPC_FLOAT* exploreMPC_ccrhsl5 = exploreMPC_ccrhs + 80;
int exploreMPC_ubIdx5[8] = {0, 1, 2, 3, 4, 5, 6, 7};
exploreMPC_FLOAT* exploreMPC_lub5 = exploreMPC_l + 88;
exploreMPC_FLOAT* exploreMPC_sub5 = exploreMPC_s + 88;
exploreMPC_FLOAT* exploreMPC_lubbysub5 = exploreMPC_lbys + 88;
exploreMPC_FLOAT exploreMPC_riub5[8];
exploreMPC_FLOAT* exploreMPC_dlubaff5 = exploreMPC_dl_aff + 88;
exploreMPC_FLOAT* exploreMPC_dsubaff5 = exploreMPC_ds_aff + 88;
exploreMPC_FLOAT* exploreMPC_dlubcc5 = exploreMPC_dl_cc + 88;
exploreMPC_FLOAT* exploreMPC_dsubcc5 = exploreMPC_ds_cc + 88;
exploreMPC_FLOAT* exploreMPC_ccrhsub5 = exploreMPC_ccrhs + 88;
exploreMPC_FLOAT exploreMPC_Phi5[8];
exploreMPC_FLOAT exploreMPC_W5[8];
exploreMPC_FLOAT exploreMPC_Ysd5[16];
exploreMPC_FLOAT exploreMPC_Lsd5[16];
exploreMPC_FLOAT* exploreMPC_z6 = exploreMPC_z + 48;
exploreMPC_FLOAT* exploreMPC_dzaff6 = exploreMPC_dz_aff + 48;
exploreMPC_FLOAT* exploreMPC_dzcc6 = exploreMPC_dz_cc + 48;
exploreMPC_FLOAT* exploreMPC_rd6 = exploreMPC_rd + 48;
exploreMPC_FLOAT exploreMPC_Lbyrd6[8];
exploreMPC_FLOAT* exploreMPC_grad_cost6 = exploreMPC_grad_cost + 48;
exploreMPC_FLOAT* exploreMPC_grad_eq6 = exploreMPC_grad_eq + 48;
exploreMPC_FLOAT* exploreMPC_grad_ineq6 = exploreMPC_grad_ineq + 48;
exploreMPC_FLOAT exploreMPC_ctv6[8];
exploreMPC_FLOAT* exploreMPC_v6 = exploreMPC_v + 24;
exploreMPC_FLOAT exploreMPC_re6[4];
exploreMPC_FLOAT exploreMPC_beta6[4];
exploreMPC_FLOAT exploreMPC_betacc6[4];
exploreMPC_FLOAT* exploreMPC_dvaff6 = exploreMPC_dv_aff + 24;
exploreMPC_FLOAT* exploreMPC_dvcc6 = exploreMPC_dv_cc + 24;
exploreMPC_FLOAT exploreMPC_V6[32];
exploreMPC_FLOAT exploreMPC_Yd6[10];
exploreMPC_FLOAT exploreMPC_Ld6[10];
exploreMPC_FLOAT exploreMPC_yy6[4];
exploreMPC_FLOAT exploreMPC_bmy6[4];
exploreMPC_FLOAT exploreMPC_c6[4] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
int exploreMPC_lbIdx6[8] = {0, 1, 2, 3, 4, 5, 6, 7};
exploreMPC_FLOAT* exploreMPC_llb6 = exploreMPC_l + 96;
exploreMPC_FLOAT* exploreMPC_slb6 = exploreMPC_s + 96;
exploreMPC_FLOAT* exploreMPC_llbbyslb6 = exploreMPC_lbys + 96;
exploreMPC_FLOAT exploreMPC_rilb6[8];
exploreMPC_FLOAT* exploreMPC_dllbaff6 = exploreMPC_dl_aff + 96;
exploreMPC_FLOAT* exploreMPC_dslbaff6 = exploreMPC_ds_aff + 96;
exploreMPC_FLOAT* exploreMPC_dllbcc6 = exploreMPC_dl_cc + 96;
exploreMPC_FLOAT* exploreMPC_dslbcc6 = exploreMPC_ds_cc + 96;
exploreMPC_FLOAT* exploreMPC_ccrhsl6 = exploreMPC_ccrhs + 96;
int exploreMPC_ubIdx6[8] = {0, 1, 2, 3, 4, 5, 6, 7};
exploreMPC_FLOAT* exploreMPC_lub6 = exploreMPC_l + 104;
exploreMPC_FLOAT* exploreMPC_sub6 = exploreMPC_s + 104;
exploreMPC_FLOAT* exploreMPC_lubbysub6 = exploreMPC_lbys + 104;
exploreMPC_FLOAT exploreMPC_riub6[8];
exploreMPC_FLOAT* exploreMPC_dlubaff6 = exploreMPC_dl_aff + 104;
exploreMPC_FLOAT* exploreMPC_dsubaff6 = exploreMPC_ds_aff + 104;
exploreMPC_FLOAT* exploreMPC_dlubcc6 = exploreMPC_dl_cc + 104;
exploreMPC_FLOAT* exploreMPC_dsubcc6 = exploreMPC_ds_cc + 104;
exploreMPC_FLOAT* exploreMPC_ccrhsub6 = exploreMPC_ccrhs + 104;
exploreMPC_FLOAT exploreMPC_Phi6[8];
exploreMPC_FLOAT exploreMPC_W6[8];
exploreMPC_FLOAT exploreMPC_Ysd6[16];
exploreMPC_FLOAT exploreMPC_Lsd6[16];
exploreMPC_FLOAT* exploreMPC_z7 = exploreMPC_z + 56;
exploreMPC_FLOAT* exploreMPC_dzaff7 = exploreMPC_dz_aff + 56;
exploreMPC_FLOAT* exploreMPC_dzcc7 = exploreMPC_dz_cc + 56;
exploreMPC_FLOAT* exploreMPC_rd7 = exploreMPC_rd + 56;
exploreMPC_FLOAT exploreMPC_Lbyrd7[8];
exploreMPC_FLOAT* exploreMPC_grad_cost7 = exploreMPC_grad_cost + 56;
exploreMPC_FLOAT* exploreMPC_grad_eq7 = exploreMPC_grad_eq + 56;
exploreMPC_FLOAT* exploreMPC_grad_ineq7 = exploreMPC_grad_ineq + 56;
exploreMPC_FLOAT exploreMPC_ctv7[8];
exploreMPC_FLOAT* exploreMPC_v7 = exploreMPC_v + 28;
exploreMPC_FLOAT exploreMPC_re7[4];
exploreMPC_FLOAT exploreMPC_beta7[4];
exploreMPC_FLOAT exploreMPC_betacc7[4];
exploreMPC_FLOAT* exploreMPC_dvaff7 = exploreMPC_dv_aff + 28;
exploreMPC_FLOAT* exploreMPC_dvcc7 = exploreMPC_dv_cc + 28;
exploreMPC_FLOAT exploreMPC_V7[32];
exploreMPC_FLOAT exploreMPC_Yd7[10];
exploreMPC_FLOAT exploreMPC_Ld7[10];
exploreMPC_FLOAT exploreMPC_yy7[4];
exploreMPC_FLOAT exploreMPC_bmy7[4];
exploreMPC_FLOAT exploreMPC_c7[4] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
int exploreMPC_lbIdx7[8] = {0, 1, 2, 3, 4, 5, 6, 7};
exploreMPC_FLOAT* exploreMPC_llb7 = exploreMPC_l + 112;
exploreMPC_FLOAT* exploreMPC_slb7 = exploreMPC_s + 112;
exploreMPC_FLOAT* exploreMPC_llbbyslb7 = exploreMPC_lbys + 112;
exploreMPC_FLOAT exploreMPC_rilb7[8];
exploreMPC_FLOAT* exploreMPC_dllbaff7 = exploreMPC_dl_aff + 112;
exploreMPC_FLOAT* exploreMPC_dslbaff7 = exploreMPC_ds_aff + 112;
exploreMPC_FLOAT* exploreMPC_dllbcc7 = exploreMPC_dl_cc + 112;
exploreMPC_FLOAT* exploreMPC_dslbcc7 = exploreMPC_ds_cc + 112;
exploreMPC_FLOAT* exploreMPC_ccrhsl7 = exploreMPC_ccrhs + 112;
int exploreMPC_ubIdx7[8] = {0, 1, 2, 3, 4, 5, 6, 7};
exploreMPC_FLOAT* exploreMPC_lub7 = exploreMPC_l + 120;
exploreMPC_FLOAT* exploreMPC_sub7 = exploreMPC_s + 120;
exploreMPC_FLOAT* exploreMPC_lubbysub7 = exploreMPC_lbys + 120;
exploreMPC_FLOAT exploreMPC_riub7[8];
exploreMPC_FLOAT* exploreMPC_dlubaff7 = exploreMPC_dl_aff + 120;
exploreMPC_FLOAT* exploreMPC_dsubaff7 = exploreMPC_ds_aff + 120;
exploreMPC_FLOAT* exploreMPC_dlubcc7 = exploreMPC_dl_cc + 120;
exploreMPC_FLOAT* exploreMPC_dsubcc7 = exploreMPC_ds_cc + 120;
exploreMPC_FLOAT* exploreMPC_ccrhsub7 = exploreMPC_ccrhs + 120;
exploreMPC_FLOAT exploreMPC_Phi7[8];
exploreMPC_FLOAT exploreMPC_W7[8];
exploreMPC_FLOAT exploreMPC_Ysd7[16];
exploreMPC_FLOAT exploreMPC_Lsd7[16];
exploreMPC_FLOAT* exploreMPC_z8 = exploreMPC_z + 64;
exploreMPC_FLOAT* exploreMPC_dzaff8 = exploreMPC_dz_aff + 64;
exploreMPC_FLOAT* exploreMPC_dzcc8 = exploreMPC_dz_cc + 64;
exploreMPC_FLOAT* exploreMPC_rd8 = exploreMPC_rd + 64;
exploreMPC_FLOAT exploreMPC_Lbyrd8[8];
exploreMPC_FLOAT* exploreMPC_grad_cost8 = exploreMPC_grad_cost + 64;
exploreMPC_FLOAT* exploreMPC_grad_eq8 = exploreMPC_grad_eq + 64;
exploreMPC_FLOAT* exploreMPC_grad_ineq8 = exploreMPC_grad_ineq + 64;
exploreMPC_FLOAT exploreMPC_ctv8[8];
exploreMPC_FLOAT* exploreMPC_v8 = exploreMPC_v + 32;
exploreMPC_FLOAT exploreMPC_re8[4];
exploreMPC_FLOAT exploreMPC_beta8[4];
exploreMPC_FLOAT exploreMPC_betacc8[4];
exploreMPC_FLOAT* exploreMPC_dvaff8 = exploreMPC_dv_aff + 32;
exploreMPC_FLOAT* exploreMPC_dvcc8 = exploreMPC_dv_cc + 32;
exploreMPC_FLOAT exploreMPC_V8[32];
exploreMPC_FLOAT exploreMPC_Yd8[10];
exploreMPC_FLOAT exploreMPC_Ld8[10];
exploreMPC_FLOAT exploreMPC_yy8[4];
exploreMPC_FLOAT exploreMPC_bmy8[4];
exploreMPC_FLOAT exploreMPC_c8[4] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
int exploreMPC_lbIdx8[8] = {0, 1, 2, 3, 4, 5, 6, 7};
exploreMPC_FLOAT* exploreMPC_llb8 = exploreMPC_l + 128;
exploreMPC_FLOAT* exploreMPC_slb8 = exploreMPC_s + 128;
exploreMPC_FLOAT* exploreMPC_llbbyslb8 = exploreMPC_lbys + 128;
exploreMPC_FLOAT exploreMPC_rilb8[8];
exploreMPC_FLOAT* exploreMPC_dllbaff8 = exploreMPC_dl_aff + 128;
exploreMPC_FLOAT* exploreMPC_dslbaff8 = exploreMPC_ds_aff + 128;
exploreMPC_FLOAT* exploreMPC_dllbcc8 = exploreMPC_dl_cc + 128;
exploreMPC_FLOAT* exploreMPC_dslbcc8 = exploreMPC_ds_cc + 128;
exploreMPC_FLOAT* exploreMPC_ccrhsl8 = exploreMPC_ccrhs + 128;
int exploreMPC_ubIdx8[8] = {0, 1, 2, 3, 4, 5, 6, 7};
exploreMPC_FLOAT* exploreMPC_lub8 = exploreMPC_l + 136;
exploreMPC_FLOAT* exploreMPC_sub8 = exploreMPC_s + 136;
exploreMPC_FLOAT* exploreMPC_lubbysub8 = exploreMPC_lbys + 136;
exploreMPC_FLOAT exploreMPC_riub8[8];
exploreMPC_FLOAT* exploreMPC_dlubaff8 = exploreMPC_dl_aff + 136;
exploreMPC_FLOAT* exploreMPC_dsubaff8 = exploreMPC_ds_aff + 136;
exploreMPC_FLOAT* exploreMPC_dlubcc8 = exploreMPC_dl_cc + 136;
exploreMPC_FLOAT* exploreMPC_dsubcc8 = exploreMPC_ds_cc + 136;
exploreMPC_FLOAT* exploreMPC_ccrhsub8 = exploreMPC_ccrhs + 136;
exploreMPC_FLOAT exploreMPC_Phi8[8];
exploreMPC_FLOAT exploreMPC_W8[8];
exploreMPC_FLOAT exploreMPC_Ysd8[16];
exploreMPC_FLOAT exploreMPC_Lsd8[16];
exploreMPC_FLOAT* exploreMPC_z9 = exploreMPC_z + 72;
exploreMPC_FLOAT* exploreMPC_dzaff9 = exploreMPC_dz_aff + 72;
exploreMPC_FLOAT* exploreMPC_dzcc9 = exploreMPC_dz_cc + 72;
exploreMPC_FLOAT* exploreMPC_rd9 = exploreMPC_rd + 72;
exploreMPC_FLOAT exploreMPC_Lbyrd9[4];
exploreMPC_FLOAT* exploreMPC_grad_cost9 = exploreMPC_grad_cost + 72;
exploreMPC_FLOAT* exploreMPC_grad_eq9 = exploreMPC_grad_eq + 72;
exploreMPC_FLOAT* exploreMPC_grad_ineq9 = exploreMPC_grad_ineq + 72;
exploreMPC_FLOAT exploreMPC_ctv9[4];
exploreMPC_FLOAT* exploreMPC_v9 = exploreMPC_v + 36;
exploreMPC_FLOAT exploreMPC_re9[4];
exploreMPC_FLOAT exploreMPC_beta9[4];
exploreMPC_FLOAT exploreMPC_betacc9[4];
exploreMPC_FLOAT* exploreMPC_dvaff9 = exploreMPC_dv_aff + 36;
exploreMPC_FLOAT* exploreMPC_dvcc9 = exploreMPC_dv_cc + 36;
exploreMPC_FLOAT exploreMPC_V9[16];
exploreMPC_FLOAT exploreMPC_Yd9[10];
exploreMPC_FLOAT exploreMPC_Ld9[10];
exploreMPC_FLOAT exploreMPC_yy9[4];
exploreMPC_FLOAT exploreMPC_bmy9[4];
exploreMPC_FLOAT exploreMPC_c9[4] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
int exploreMPC_lbIdx9[4] = {0, 1, 2, 3};
exploreMPC_FLOAT* exploreMPC_llb9 = exploreMPC_l + 144;
exploreMPC_FLOAT* exploreMPC_slb9 = exploreMPC_s + 144;
exploreMPC_FLOAT* exploreMPC_llbbyslb9 = exploreMPC_lbys + 144;
exploreMPC_FLOAT exploreMPC_rilb9[4];
exploreMPC_FLOAT* exploreMPC_dllbaff9 = exploreMPC_dl_aff + 144;
exploreMPC_FLOAT* exploreMPC_dslbaff9 = exploreMPC_ds_aff + 144;
exploreMPC_FLOAT* exploreMPC_dllbcc9 = exploreMPC_dl_cc + 144;
exploreMPC_FLOAT* exploreMPC_dslbcc9 = exploreMPC_ds_cc + 144;
exploreMPC_FLOAT* exploreMPC_ccrhsl9 = exploreMPC_ccrhs + 144;
int exploreMPC_ubIdx9[4] = {0, 1, 2, 3};
exploreMPC_FLOAT* exploreMPC_lub9 = exploreMPC_l + 148;
exploreMPC_FLOAT* exploreMPC_sub9 = exploreMPC_s + 148;
exploreMPC_FLOAT* exploreMPC_lubbysub9 = exploreMPC_lbys + 148;
exploreMPC_FLOAT exploreMPC_riub9[4];
exploreMPC_FLOAT* exploreMPC_dlubaff9 = exploreMPC_dl_aff + 148;
exploreMPC_FLOAT* exploreMPC_dsubaff9 = exploreMPC_ds_aff + 148;
exploreMPC_FLOAT* exploreMPC_dlubcc9 = exploreMPC_dl_cc + 148;
exploreMPC_FLOAT* exploreMPC_dsubcc9 = exploreMPC_ds_cc + 148;
exploreMPC_FLOAT* exploreMPC_ccrhsub9 = exploreMPC_ccrhs + 148;
exploreMPC_FLOAT exploreMPC_Phi9[4];
exploreMPC_FLOAT exploreMPC_D9[4] = {-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000};
exploreMPC_FLOAT exploreMPC_W9[4];
exploreMPC_FLOAT exploreMPC_Ysd9[16];
exploreMPC_FLOAT exploreMPC_Lsd9[16];
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
exploreMPC_LA_INITIALIZEVECTOR_76(exploreMPC_z, 0);
exploreMPC_LA_INITIALIZEVECTOR_40(exploreMPC_v, 1);
exploreMPC_LA_INITIALIZEVECTOR_152(exploreMPC_l, 1);
exploreMPC_LA_INITIALIZEVECTOR_152(exploreMPC_s, 1);
info->mu = 0;
exploreMPC_LA_DOTACC_152(exploreMPC_l, exploreMPC_s, &info->mu);
info->mu /= 152;
while( 1 ){
info->pobj = 0;
exploreMPC_LA_DIAG_QUADFCN_8(params->H1, params->f1, exploreMPC_z0, exploreMPC_grad_cost0, &info->pobj);
exploreMPC_LA_DIAG_QUADFCN_8(params->H2, params->f2, exploreMPC_z1, exploreMPC_grad_cost1, &info->pobj);
exploreMPC_LA_DIAG_QUADFCN_8(params->H3, params->f3, exploreMPC_z2, exploreMPC_grad_cost2, &info->pobj);
exploreMPC_LA_DIAG_QUADFCN_8(params->H4, params->f4, exploreMPC_z3, exploreMPC_grad_cost3, &info->pobj);
exploreMPC_LA_DIAG_QUADFCN_8(params->H5, params->f5, exploreMPC_z4, exploreMPC_grad_cost4, &info->pobj);
exploreMPC_LA_DIAG_QUADFCN_8(params->H6, params->f6, exploreMPC_z5, exploreMPC_grad_cost5, &info->pobj);
exploreMPC_LA_DIAG_QUADFCN_8(params->H7, params->f7, exploreMPC_z6, exploreMPC_grad_cost6, &info->pobj);
exploreMPC_LA_DIAG_QUADFCN_8(params->H8, params->f8, exploreMPC_z7, exploreMPC_grad_cost7, &info->pobj);
exploreMPC_LA_DIAG_QUADFCN_8(params->H9, params->f9, exploreMPC_z8, exploreMPC_grad_cost8, &info->pobj);
exploreMPC_LA_DIAG_QUADFCN_4(params->H10, params->f10, exploreMPC_z9, exploreMPC_grad_cost9, &info->pobj);
info->res_eq = 0;
info->dgap = 0;
exploreMPC_LA_DIAGZERO_MVMSUB6_4(exploreMPC_D0, exploreMPC_z0, params->c1, exploreMPC_v0, exploreMPC_re0, &info->dgap, &info->res_eq);
exploreMPC_LA_DENSE_DIAGZERO_MVMSUB3_4_8_8(exploreMPC_C0, exploreMPC_z0, exploreMPC_D1, exploreMPC_z1, exploreMPC_c1, exploreMPC_v1, exploreMPC_re1, &info->dgap, &info->res_eq);
exploreMPC_LA_DENSE_DIAGZERO_MVMSUB3_4_8_8(exploreMPC_C0, exploreMPC_z1, exploreMPC_D1, exploreMPC_z2, exploreMPC_c2, exploreMPC_v2, exploreMPC_re2, &info->dgap, &info->res_eq);
exploreMPC_LA_DENSE_DIAGZERO_MVMSUB3_4_8_8(exploreMPC_C0, exploreMPC_z2, exploreMPC_D1, exploreMPC_z3, exploreMPC_c3, exploreMPC_v3, exploreMPC_re3, &info->dgap, &info->res_eq);
exploreMPC_LA_DENSE_DIAGZERO_MVMSUB3_4_8_8(exploreMPC_C0, exploreMPC_z3, exploreMPC_D1, exploreMPC_z4, exploreMPC_c4, exploreMPC_v4, exploreMPC_re4, &info->dgap, &info->res_eq);
exploreMPC_LA_DENSE_DIAGZERO_MVMSUB3_4_8_8(exploreMPC_C0, exploreMPC_z4, exploreMPC_D1, exploreMPC_z5, exploreMPC_c5, exploreMPC_v5, exploreMPC_re5, &info->dgap, &info->res_eq);
exploreMPC_LA_DENSE_DIAGZERO_MVMSUB3_4_8_8(exploreMPC_C0, exploreMPC_z5, exploreMPC_D1, exploreMPC_z6, exploreMPC_c6, exploreMPC_v6, exploreMPC_re6, &info->dgap, &info->res_eq);
exploreMPC_LA_DENSE_DIAGZERO_MVMSUB3_4_8_8(exploreMPC_C0, exploreMPC_z6, exploreMPC_D1, exploreMPC_z7, exploreMPC_c7, exploreMPC_v7, exploreMPC_re7, &info->dgap, &info->res_eq);
exploreMPC_LA_DENSE_DIAGZERO_MVMSUB3_4_8_8(exploreMPC_C0, exploreMPC_z7, exploreMPC_D1, exploreMPC_z8, exploreMPC_c8, exploreMPC_v8, exploreMPC_re8, &info->dgap, &info->res_eq);
exploreMPC_LA_DENSE_DIAGZERO_MVMSUB3_4_8_4(exploreMPC_C0, exploreMPC_z8, exploreMPC_D9, exploreMPC_z9, exploreMPC_c9, exploreMPC_v9, exploreMPC_re9, &info->dgap, &info->res_eq);
exploreMPC_LA_DENSE_DIAGZERO_MTVM2_4_8_4(exploreMPC_C0, exploreMPC_v1, exploreMPC_D0, exploreMPC_v0, exploreMPC_grad_eq0);
exploreMPC_LA_DENSE_DIAGZERO_MTVM2_4_8_4(exploreMPC_C0, exploreMPC_v2, exploreMPC_D1, exploreMPC_v1, exploreMPC_grad_eq1);
exploreMPC_LA_DENSE_DIAGZERO_MTVM2_4_8_4(exploreMPC_C0, exploreMPC_v3, exploreMPC_D1, exploreMPC_v2, exploreMPC_grad_eq2);
exploreMPC_LA_DENSE_DIAGZERO_MTVM2_4_8_4(exploreMPC_C0, exploreMPC_v4, exploreMPC_D1, exploreMPC_v3, exploreMPC_grad_eq3);
exploreMPC_LA_DENSE_DIAGZERO_MTVM2_4_8_4(exploreMPC_C0, exploreMPC_v5, exploreMPC_D1, exploreMPC_v4, exploreMPC_grad_eq4);
exploreMPC_LA_DENSE_DIAGZERO_MTVM2_4_8_4(exploreMPC_C0, exploreMPC_v6, exploreMPC_D1, exploreMPC_v5, exploreMPC_grad_eq5);
exploreMPC_LA_DENSE_DIAGZERO_MTVM2_4_8_4(exploreMPC_C0, exploreMPC_v7, exploreMPC_D1, exploreMPC_v6, exploreMPC_grad_eq6);
exploreMPC_LA_DENSE_DIAGZERO_MTVM2_4_8_4(exploreMPC_C0, exploreMPC_v8, exploreMPC_D1, exploreMPC_v7, exploreMPC_grad_eq7);
exploreMPC_LA_DENSE_DIAGZERO_MTVM2_4_8_4(exploreMPC_C0, exploreMPC_v9, exploreMPC_D1, exploreMPC_v8, exploreMPC_grad_eq8);
exploreMPC_LA_DIAGZERO_MTVM_4_4(exploreMPC_D9, exploreMPC_v9, exploreMPC_grad_eq9);
info->res_ineq = 0;
exploreMPC_LA_VSUBADD3_8(params->lb1, exploreMPC_z0, exploreMPC_lbIdx0, exploreMPC_llb0, exploreMPC_slb0, exploreMPC_rilb0, &info->dgap, &info->res_ineq);
exploreMPC_LA_VSUBADD2_8(exploreMPC_z0, exploreMPC_ubIdx0, params->ub1, exploreMPC_lub0, exploreMPC_sub0, exploreMPC_riub0, &info->dgap, &info->res_ineq);
exploreMPC_LA_VSUBADD3_8(params->lb2, exploreMPC_z1, exploreMPC_lbIdx1, exploreMPC_llb1, exploreMPC_slb1, exploreMPC_rilb1, &info->dgap, &info->res_ineq);
exploreMPC_LA_VSUBADD2_8(exploreMPC_z1, exploreMPC_ubIdx1, params->ub2, exploreMPC_lub1, exploreMPC_sub1, exploreMPC_riub1, &info->dgap, &info->res_ineq);
exploreMPC_LA_VSUBADD3_8(params->lb3, exploreMPC_z2, exploreMPC_lbIdx2, exploreMPC_llb2, exploreMPC_slb2, exploreMPC_rilb2, &info->dgap, &info->res_ineq);
exploreMPC_LA_VSUBADD2_8(exploreMPC_z2, exploreMPC_ubIdx2, params->ub3, exploreMPC_lub2, exploreMPC_sub2, exploreMPC_riub2, &info->dgap, &info->res_ineq);
exploreMPC_LA_VSUBADD3_8(params->lb4, exploreMPC_z3, exploreMPC_lbIdx3, exploreMPC_llb3, exploreMPC_slb3, exploreMPC_rilb3, &info->dgap, &info->res_ineq);
exploreMPC_LA_VSUBADD2_8(exploreMPC_z3, exploreMPC_ubIdx3, params->ub4, exploreMPC_lub3, exploreMPC_sub3, exploreMPC_riub3, &info->dgap, &info->res_ineq);
exploreMPC_LA_VSUBADD3_8(params->lb5, exploreMPC_z4, exploreMPC_lbIdx4, exploreMPC_llb4, exploreMPC_slb4, exploreMPC_rilb4, &info->dgap, &info->res_ineq);
exploreMPC_LA_VSUBADD2_8(exploreMPC_z4, exploreMPC_ubIdx4, params->ub5, exploreMPC_lub4, exploreMPC_sub4, exploreMPC_riub4, &info->dgap, &info->res_ineq);
exploreMPC_LA_VSUBADD3_8(params->lb6, exploreMPC_z5, exploreMPC_lbIdx5, exploreMPC_llb5, exploreMPC_slb5, exploreMPC_rilb5, &info->dgap, &info->res_ineq);
exploreMPC_LA_VSUBADD2_8(exploreMPC_z5, exploreMPC_ubIdx5, params->ub6, exploreMPC_lub5, exploreMPC_sub5, exploreMPC_riub5, &info->dgap, &info->res_ineq);
exploreMPC_LA_VSUBADD3_8(params->lb7, exploreMPC_z6, exploreMPC_lbIdx6, exploreMPC_llb6, exploreMPC_slb6, exploreMPC_rilb6, &info->dgap, &info->res_ineq);
exploreMPC_LA_VSUBADD2_8(exploreMPC_z6, exploreMPC_ubIdx6, params->ub7, exploreMPC_lub6, exploreMPC_sub6, exploreMPC_riub6, &info->dgap, &info->res_ineq);
exploreMPC_LA_VSUBADD3_8(params->lb8, exploreMPC_z7, exploreMPC_lbIdx7, exploreMPC_llb7, exploreMPC_slb7, exploreMPC_rilb7, &info->dgap, &info->res_ineq);
exploreMPC_LA_VSUBADD2_8(exploreMPC_z7, exploreMPC_ubIdx7, params->ub8, exploreMPC_lub7, exploreMPC_sub7, exploreMPC_riub7, &info->dgap, &info->res_ineq);
exploreMPC_LA_VSUBADD3_8(params->lb9, exploreMPC_z8, exploreMPC_lbIdx8, exploreMPC_llb8, exploreMPC_slb8, exploreMPC_rilb8, &info->dgap, &info->res_ineq);
exploreMPC_LA_VSUBADD2_8(exploreMPC_z8, exploreMPC_ubIdx8, params->ub9, exploreMPC_lub8, exploreMPC_sub8, exploreMPC_riub8, &info->dgap, &info->res_ineq);
exploreMPC_LA_VSUBADD3_4(params->lb10, exploreMPC_z9, exploreMPC_lbIdx9, exploreMPC_llb9, exploreMPC_slb9, exploreMPC_rilb9, &info->dgap, &info->res_ineq);
exploreMPC_LA_VSUBADD2_4(exploreMPC_z9, exploreMPC_ubIdx9, params->ub10, exploreMPC_lub9, exploreMPC_sub9, exploreMPC_riub9, &info->dgap, &info->res_ineq);
exploreMPC_LA_INEQ_B_GRAD_8_8_8(exploreMPC_lub0, exploreMPC_sub0, exploreMPC_riub0, exploreMPC_llb0, exploreMPC_slb0, exploreMPC_rilb0, exploreMPC_lbIdx0, exploreMPC_ubIdx0, exploreMPC_grad_ineq0, exploreMPC_lubbysub0, exploreMPC_llbbyslb0);
exploreMPC_LA_INEQ_B_GRAD_8_8_8(exploreMPC_lub1, exploreMPC_sub1, exploreMPC_riub1, exploreMPC_llb1, exploreMPC_slb1, exploreMPC_rilb1, exploreMPC_lbIdx1, exploreMPC_ubIdx1, exploreMPC_grad_ineq1, exploreMPC_lubbysub1, exploreMPC_llbbyslb1);
exploreMPC_LA_INEQ_B_GRAD_8_8_8(exploreMPC_lub2, exploreMPC_sub2, exploreMPC_riub2, exploreMPC_llb2, exploreMPC_slb2, exploreMPC_rilb2, exploreMPC_lbIdx2, exploreMPC_ubIdx2, exploreMPC_grad_ineq2, exploreMPC_lubbysub2, exploreMPC_llbbyslb2);
exploreMPC_LA_INEQ_B_GRAD_8_8_8(exploreMPC_lub3, exploreMPC_sub3, exploreMPC_riub3, exploreMPC_llb3, exploreMPC_slb3, exploreMPC_rilb3, exploreMPC_lbIdx3, exploreMPC_ubIdx3, exploreMPC_grad_ineq3, exploreMPC_lubbysub3, exploreMPC_llbbyslb3);
exploreMPC_LA_INEQ_B_GRAD_8_8_8(exploreMPC_lub4, exploreMPC_sub4, exploreMPC_riub4, exploreMPC_llb4, exploreMPC_slb4, exploreMPC_rilb4, exploreMPC_lbIdx4, exploreMPC_ubIdx4, exploreMPC_grad_ineq4, exploreMPC_lubbysub4, exploreMPC_llbbyslb4);
exploreMPC_LA_INEQ_B_GRAD_8_8_8(exploreMPC_lub5, exploreMPC_sub5, exploreMPC_riub5, exploreMPC_llb5, exploreMPC_slb5, exploreMPC_rilb5, exploreMPC_lbIdx5, exploreMPC_ubIdx5, exploreMPC_grad_ineq5, exploreMPC_lubbysub5, exploreMPC_llbbyslb5);
exploreMPC_LA_INEQ_B_GRAD_8_8_8(exploreMPC_lub6, exploreMPC_sub6, exploreMPC_riub6, exploreMPC_llb6, exploreMPC_slb6, exploreMPC_rilb6, exploreMPC_lbIdx6, exploreMPC_ubIdx6, exploreMPC_grad_ineq6, exploreMPC_lubbysub6, exploreMPC_llbbyslb6);
exploreMPC_LA_INEQ_B_GRAD_8_8_8(exploreMPC_lub7, exploreMPC_sub7, exploreMPC_riub7, exploreMPC_llb7, exploreMPC_slb7, exploreMPC_rilb7, exploreMPC_lbIdx7, exploreMPC_ubIdx7, exploreMPC_grad_ineq7, exploreMPC_lubbysub7, exploreMPC_llbbyslb7);
exploreMPC_LA_INEQ_B_GRAD_8_8_8(exploreMPC_lub8, exploreMPC_sub8, exploreMPC_riub8, exploreMPC_llb8, exploreMPC_slb8, exploreMPC_rilb8, exploreMPC_lbIdx8, exploreMPC_ubIdx8, exploreMPC_grad_ineq8, exploreMPC_lubbysub8, exploreMPC_llbbyslb8);
exploreMPC_LA_INEQ_B_GRAD_4_4_4(exploreMPC_lub9, exploreMPC_sub9, exploreMPC_riub9, exploreMPC_llb9, exploreMPC_slb9, exploreMPC_rilb9, exploreMPC_lbIdx9, exploreMPC_ubIdx9, exploreMPC_grad_ineq9, exploreMPC_lubbysub9, exploreMPC_llbbyslb9);
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
exploreMPC_LA_VVADD3_76(exploreMPC_grad_cost, exploreMPC_grad_eq, exploreMPC_grad_ineq, exploreMPC_rd);
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
exploreMPC_LA_DIAG_CHOL_ONELOOP_LBUB_8_8_8(params->H3, exploreMPC_llbbyslb2, exploreMPC_lbIdx2, exploreMPC_lubbysub2, exploreMPC_ubIdx2, exploreMPC_Phi2);
exploreMPC_LA_DIAG_MATRIXFORWARDSUB_4_8(exploreMPC_Phi2, exploreMPC_C0, exploreMPC_V2);
exploreMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_4_8(exploreMPC_Phi2, exploreMPC_D1, exploreMPC_W2);
exploreMPC_LA_DENSE_DIAGZERO_MMTM_4_8_4(exploreMPC_W2, exploreMPC_V2, exploreMPC_Ysd3);
exploreMPC_LA_DIAG_FORWARDSUB_8(exploreMPC_Phi2, exploreMPC_rd2, exploreMPC_Lbyrd2);
exploreMPC_LA_DIAG_CHOL_ONELOOP_LBUB_8_8_8(params->H4, exploreMPC_llbbyslb3, exploreMPC_lbIdx3, exploreMPC_lubbysub3, exploreMPC_ubIdx3, exploreMPC_Phi3);
exploreMPC_LA_DIAG_MATRIXFORWARDSUB_4_8(exploreMPC_Phi3, exploreMPC_C0, exploreMPC_V3);
exploreMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_4_8(exploreMPC_Phi3, exploreMPC_D1, exploreMPC_W3);
exploreMPC_LA_DENSE_DIAGZERO_MMTM_4_8_4(exploreMPC_W3, exploreMPC_V3, exploreMPC_Ysd4);
exploreMPC_LA_DIAG_FORWARDSUB_8(exploreMPC_Phi3, exploreMPC_rd3, exploreMPC_Lbyrd3);
exploreMPC_LA_DIAG_CHOL_ONELOOP_LBUB_8_8_8(params->H5, exploreMPC_llbbyslb4, exploreMPC_lbIdx4, exploreMPC_lubbysub4, exploreMPC_ubIdx4, exploreMPC_Phi4);
exploreMPC_LA_DIAG_MATRIXFORWARDSUB_4_8(exploreMPC_Phi4, exploreMPC_C0, exploreMPC_V4);
exploreMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_4_8(exploreMPC_Phi4, exploreMPC_D1, exploreMPC_W4);
exploreMPC_LA_DENSE_DIAGZERO_MMTM_4_8_4(exploreMPC_W4, exploreMPC_V4, exploreMPC_Ysd5);
exploreMPC_LA_DIAG_FORWARDSUB_8(exploreMPC_Phi4, exploreMPC_rd4, exploreMPC_Lbyrd4);
exploreMPC_LA_DIAG_CHOL_ONELOOP_LBUB_8_8_8(params->H6, exploreMPC_llbbyslb5, exploreMPC_lbIdx5, exploreMPC_lubbysub5, exploreMPC_ubIdx5, exploreMPC_Phi5);
exploreMPC_LA_DIAG_MATRIXFORWARDSUB_4_8(exploreMPC_Phi5, exploreMPC_C0, exploreMPC_V5);
exploreMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_4_8(exploreMPC_Phi5, exploreMPC_D1, exploreMPC_W5);
exploreMPC_LA_DENSE_DIAGZERO_MMTM_4_8_4(exploreMPC_W5, exploreMPC_V5, exploreMPC_Ysd6);
exploreMPC_LA_DIAG_FORWARDSUB_8(exploreMPC_Phi5, exploreMPC_rd5, exploreMPC_Lbyrd5);
exploreMPC_LA_DIAG_CHOL_ONELOOP_LBUB_8_8_8(params->H7, exploreMPC_llbbyslb6, exploreMPC_lbIdx6, exploreMPC_lubbysub6, exploreMPC_ubIdx6, exploreMPC_Phi6);
exploreMPC_LA_DIAG_MATRIXFORWARDSUB_4_8(exploreMPC_Phi6, exploreMPC_C0, exploreMPC_V6);
exploreMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_4_8(exploreMPC_Phi6, exploreMPC_D1, exploreMPC_W6);
exploreMPC_LA_DENSE_DIAGZERO_MMTM_4_8_4(exploreMPC_W6, exploreMPC_V6, exploreMPC_Ysd7);
exploreMPC_LA_DIAG_FORWARDSUB_8(exploreMPC_Phi6, exploreMPC_rd6, exploreMPC_Lbyrd6);
exploreMPC_LA_DIAG_CHOL_ONELOOP_LBUB_8_8_8(params->H8, exploreMPC_llbbyslb7, exploreMPC_lbIdx7, exploreMPC_lubbysub7, exploreMPC_ubIdx7, exploreMPC_Phi7);
exploreMPC_LA_DIAG_MATRIXFORWARDSUB_4_8(exploreMPC_Phi7, exploreMPC_C0, exploreMPC_V7);
exploreMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_4_8(exploreMPC_Phi7, exploreMPC_D1, exploreMPC_W7);
exploreMPC_LA_DENSE_DIAGZERO_MMTM_4_8_4(exploreMPC_W7, exploreMPC_V7, exploreMPC_Ysd8);
exploreMPC_LA_DIAG_FORWARDSUB_8(exploreMPC_Phi7, exploreMPC_rd7, exploreMPC_Lbyrd7);
exploreMPC_LA_DIAG_CHOL_ONELOOP_LBUB_8_8_8(params->H9, exploreMPC_llbbyslb8, exploreMPC_lbIdx8, exploreMPC_lubbysub8, exploreMPC_ubIdx8, exploreMPC_Phi8);
exploreMPC_LA_DIAG_MATRIXFORWARDSUB_4_8(exploreMPC_Phi8, exploreMPC_C0, exploreMPC_V8);
exploreMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_4_8(exploreMPC_Phi8, exploreMPC_D1, exploreMPC_W8);
exploreMPC_LA_DENSE_DIAGZERO_MMTM_4_8_4(exploreMPC_W8, exploreMPC_V8, exploreMPC_Ysd9);
exploreMPC_LA_DIAG_FORWARDSUB_8(exploreMPC_Phi8, exploreMPC_rd8, exploreMPC_Lbyrd8);
exploreMPC_LA_DIAG_CHOL_ONELOOP_LBUB_4_4_4(params->H10, exploreMPC_llbbyslb9, exploreMPC_lbIdx9, exploreMPC_lubbysub9, exploreMPC_ubIdx9, exploreMPC_Phi9);
exploreMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_4_4(exploreMPC_Phi9, exploreMPC_D9, exploreMPC_W9);
exploreMPC_LA_DIAG_FORWARDSUB_4(exploreMPC_Phi9, exploreMPC_rd9, exploreMPC_Lbyrd9);
exploreMPC_LA_DIAGZERO_MMT_4(exploreMPC_W0, exploreMPC_Yd0);
exploreMPC_LA_DIAGZERO_MVMSUB7_4(exploreMPC_W0, exploreMPC_Lbyrd0, exploreMPC_re0, exploreMPC_beta0);
exploreMPC_LA_DENSE_DIAGZERO_MMT2_4_8_8(exploreMPC_V0, exploreMPC_W1, exploreMPC_Yd1);
exploreMPC_LA_DENSE_DIAGZERO_2MVMSUB2_4_8_8(exploreMPC_V0, exploreMPC_Lbyrd0, exploreMPC_W1, exploreMPC_Lbyrd1, exploreMPC_re1, exploreMPC_beta1);
exploreMPC_LA_DENSE_DIAGZERO_MMT2_4_8_8(exploreMPC_V1, exploreMPC_W2, exploreMPC_Yd2);
exploreMPC_LA_DENSE_DIAGZERO_2MVMSUB2_4_8_8(exploreMPC_V1, exploreMPC_Lbyrd1, exploreMPC_W2, exploreMPC_Lbyrd2, exploreMPC_re2, exploreMPC_beta2);
exploreMPC_LA_DENSE_DIAGZERO_MMT2_4_8_8(exploreMPC_V2, exploreMPC_W3, exploreMPC_Yd3);
exploreMPC_LA_DENSE_DIAGZERO_2MVMSUB2_4_8_8(exploreMPC_V2, exploreMPC_Lbyrd2, exploreMPC_W3, exploreMPC_Lbyrd3, exploreMPC_re3, exploreMPC_beta3);
exploreMPC_LA_DENSE_DIAGZERO_MMT2_4_8_8(exploreMPC_V3, exploreMPC_W4, exploreMPC_Yd4);
exploreMPC_LA_DENSE_DIAGZERO_2MVMSUB2_4_8_8(exploreMPC_V3, exploreMPC_Lbyrd3, exploreMPC_W4, exploreMPC_Lbyrd4, exploreMPC_re4, exploreMPC_beta4);
exploreMPC_LA_DENSE_DIAGZERO_MMT2_4_8_8(exploreMPC_V4, exploreMPC_W5, exploreMPC_Yd5);
exploreMPC_LA_DENSE_DIAGZERO_2MVMSUB2_4_8_8(exploreMPC_V4, exploreMPC_Lbyrd4, exploreMPC_W5, exploreMPC_Lbyrd5, exploreMPC_re5, exploreMPC_beta5);
exploreMPC_LA_DENSE_DIAGZERO_MMT2_4_8_8(exploreMPC_V5, exploreMPC_W6, exploreMPC_Yd6);
exploreMPC_LA_DENSE_DIAGZERO_2MVMSUB2_4_8_8(exploreMPC_V5, exploreMPC_Lbyrd5, exploreMPC_W6, exploreMPC_Lbyrd6, exploreMPC_re6, exploreMPC_beta6);
exploreMPC_LA_DENSE_DIAGZERO_MMT2_4_8_8(exploreMPC_V6, exploreMPC_W7, exploreMPC_Yd7);
exploreMPC_LA_DENSE_DIAGZERO_2MVMSUB2_4_8_8(exploreMPC_V6, exploreMPC_Lbyrd6, exploreMPC_W7, exploreMPC_Lbyrd7, exploreMPC_re7, exploreMPC_beta7);
exploreMPC_LA_DENSE_DIAGZERO_MMT2_4_8_8(exploreMPC_V7, exploreMPC_W8, exploreMPC_Yd8);
exploreMPC_LA_DENSE_DIAGZERO_2MVMSUB2_4_8_8(exploreMPC_V7, exploreMPC_Lbyrd7, exploreMPC_W8, exploreMPC_Lbyrd8, exploreMPC_re8, exploreMPC_beta8);
exploreMPC_LA_DENSE_DIAGZERO_MMT2_4_8_4(exploreMPC_V8, exploreMPC_W9, exploreMPC_Yd9);
exploreMPC_LA_DENSE_DIAGZERO_2MVMSUB2_4_8_4(exploreMPC_V8, exploreMPC_Lbyrd8, exploreMPC_W9, exploreMPC_Lbyrd9, exploreMPC_re9, exploreMPC_beta9);
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
exploreMPC_LA_DENSE_MATRIXTFORWARDSUB_4_4(exploreMPC_Ld2, exploreMPC_Ysd3, exploreMPC_Lsd3);
exploreMPC_LA_DENSE_MMTSUB_4_4(exploreMPC_Lsd3, exploreMPC_Yd3);
exploreMPC_LA_DENSE_CHOL_4(exploreMPC_Yd3, exploreMPC_Ld3);
exploreMPC_LA_DENSE_MVMSUB1_4_4(exploreMPC_Lsd3, exploreMPC_yy2, exploreMPC_beta3, exploreMPC_bmy3);
exploreMPC_LA_DENSE_FORWARDSUB_4(exploreMPC_Ld3, exploreMPC_bmy3, exploreMPC_yy3);
exploreMPC_LA_DENSE_MATRIXTFORWARDSUB_4_4(exploreMPC_Ld3, exploreMPC_Ysd4, exploreMPC_Lsd4);
exploreMPC_LA_DENSE_MMTSUB_4_4(exploreMPC_Lsd4, exploreMPC_Yd4);
exploreMPC_LA_DENSE_CHOL_4(exploreMPC_Yd4, exploreMPC_Ld4);
exploreMPC_LA_DENSE_MVMSUB1_4_4(exploreMPC_Lsd4, exploreMPC_yy3, exploreMPC_beta4, exploreMPC_bmy4);
exploreMPC_LA_DENSE_FORWARDSUB_4(exploreMPC_Ld4, exploreMPC_bmy4, exploreMPC_yy4);
exploreMPC_LA_DENSE_MATRIXTFORWARDSUB_4_4(exploreMPC_Ld4, exploreMPC_Ysd5, exploreMPC_Lsd5);
exploreMPC_LA_DENSE_MMTSUB_4_4(exploreMPC_Lsd5, exploreMPC_Yd5);
exploreMPC_LA_DENSE_CHOL_4(exploreMPC_Yd5, exploreMPC_Ld5);
exploreMPC_LA_DENSE_MVMSUB1_4_4(exploreMPC_Lsd5, exploreMPC_yy4, exploreMPC_beta5, exploreMPC_bmy5);
exploreMPC_LA_DENSE_FORWARDSUB_4(exploreMPC_Ld5, exploreMPC_bmy5, exploreMPC_yy5);
exploreMPC_LA_DENSE_MATRIXTFORWARDSUB_4_4(exploreMPC_Ld5, exploreMPC_Ysd6, exploreMPC_Lsd6);
exploreMPC_LA_DENSE_MMTSUB_4_4(exploreMPC_Lsd6, exploreMPC_Yd6);
exploreMPC_LA_DENSE_CHOL_4(exploreMPC_Yd6, exploreMPC_Ld6);
exploreMPC_LA_DENSE_MVMSUB1_4_4(exploreMPC_Lsd6, exploreMPC_yy5, exploreMPC_beta6, exploreMPC_bmy6);
exploreMPC_LA_DENSE_FORWARDSUB_4(exploreMPC_Ld6, exploreMPC_bmy6, exploreMPC_yy6);
exploreMPC_LA_DENSE_MATRIXTFORWARDSUB_4_4(exploreMPC_Ld6, exploreMPC_Ysd7, exploreMPC_Lsd7);
exploreMPC_LA_DENSE_MMTSUB_4_4(exploreMPC_Lsd7, exploreMPC_Yd7);
exploreMPC_LA_DENSE_CHOL_4(exploreMPC_Yd7, exploreMPC_Ld7);
exploreMPC_LA_DENSE_MVMSUB1_4_4(exploreMPC_Lsd7, exploreMPC_yy6, exploreMPC_beta7, exploreMPC_bmy7);
exploreMPC_LA_DENSE_FORWARDSUB_4(exploreMPC_Ld7, exploreMPC_bmy7, exploreMPC_yy7);
exploreMPC_LA_DENSE_MATRIXTFORWARDSUB_4_4(exploreMPC_Ld7, exploreMPC_Ysd8, exploreMPC_Lsd8);
exploreMPC_LA_DENSE_MMTSUB_4_4(exploreMPC_Lsd8, exploreMPC_Yd8);
exploreMPC_LA_DENSE_CHOL_4(exploreMPC_Yd8, exploreMPC_Ld8);
exploreMPC_LA_DENSE_MVMSUB1_4_4(exploreMPC_Lsd8, exploreMPC_yy7, exploreMPC_beta8, exploreMPC_bmy8);
exploreMPC_LA_DENSE_FORWARDSUB_4(exploreMPC_Ld8, exploreMPC_bmy8, exploreMPC_yy8);
exploreMPC_LA_DENSE_MATRIXTFORWARDSUB_4_4(exploreMPC_Ld8, exploreMPC_Ysd9, exploreMPC_Lsd9);
exploreMPC_LA_DENSE_MMTSUB_4_4(exploreMPC_Lsd9, exploreMPC_Yd9);
exploreMPC_LA_DENSE_CHOL_4(exploreMPC_Yd9, exploreMPC_Ld9);
exploreMPC_LA_DENSE_MVMSUB1_4_4(exploreMPC_Lsd9, exploreMPC_yy8, exploreMPC_beta9, exploreMPC_bmy9);
exploreMPC_LA_DENSE_FORWARDSUB_4(exploreMPC_Ld9, exploreMPC_bmy9, exploreMPC_yy9);
exploreMPC_LA_DENSE_BACKWARDSUB_4(exploreMPC_Ld9, exploreMPC_yy9, exploreMPC_dvaff9);
exploreMPC_LA_DENSE_MTVMSUB_4_4(exploreMPC_Lsd9, exploreMPC_dvaff9, exploreMPC_yy8, exploreMPC_bmy8);
exploreMPC_LA_DENSE_BACKWARDSUB_4(exploreMPC_Ld8, exploreMPC_bmy8, exploreMPC_dvaff8);
exploreMPC_LA_DENSE_MTVMSUB_4_4(exploreMPC_Lsd8, exploreMPC_dvaff8, exploreMPC_yy7, exploreMPC_bmy7);
exploreMPC_LA_DENSE_BACKWARDSUB_4(exploreMPC_Ld7, exploreMPC_bmy7, exploreMPC_dvaff7);
exploreMPC_LA_DENSE_MTVMSUB_4_4(exploreMPC_Lsd7, exploreMPC_dvaff7, exploreMPC_yy6, exploreMPC_bmy6);
exploreMPC_LA_DENSE_BACKWARDSUB_4(exploreMPC_Ld6, exploreMPC_bmy6, exploreMPC_dvaff6);
exploreMPC_LA_DENSE_MTVMSUB_4_4(exploreMPC_Lsd6, exploreMPC_dvaff6, exploreMPC_yy5, exploreMPC_bmy5);
exploreMPC_LA_DENSE_BACKWARDSUB_4(exploreMPC_Ld5, exploreMPC_bmy5, exploreMPC_dvaff5);
exploreMPC_LA_DENSE_MTVMSUB_4_4(exploreMPC_Lsd5, exploreMPC_dvaff5, exploreMPC_yy4, exploreMPC_bmy4);
exploreMPC_LA_DENSE_BACKWARDSUB_4(exploreMPC_Ld4, exploreMPC_bmy4, exploreMPC_dvaff4);
exploreMPC_LA_DENSE_MTVMSUB_4_4(exploreMPC_Lsd4, exploreMPC_dvaff4, exploreMPC_yy3, exploreMPC_bmy3);
exploreMPC_LA_DENSE_BACKWARDSUB_4(exploreMPC_Ld3, exploreMPC_bmy3, exploreMPC_dvaff3);
exploreMPC_LA_DENSE_MTVMSUB_4_4(exploreMPC_Lsd3, exploreMPC_dvaff3, exploreMPC_yy2, exploreMPC_bmy2);
exploreMPC_LA_DENSE_BACKWARDSUB_4(exploreMPC_Ld2, exploreMPC_bmy2, exploreMPC_dvaff2);
exploreMPC_LA_DENSE_MTVMSUB_4_4(exploreMPC_Lsd2, exploreMPC_dvaff2, exploreMPC_yy1, exploreMPC_bmy1);
exploreMPC_LA_DENSE_BACKWARDSUB_4(exploreMPC_Ld1, exploreMPC_bmy1, exploreMPC_dvaff1);
exploreMPC_LA_DENSE_MTVMSUB_4_4(exploreMPC_Lsd1, exploreMPC_dvaff1, exploreMPC_yy0, exploreMPC_bmy0);
exploreMPC_LA_DENSE_BACKWARDSUB_4(exploreMPC_Ld0, exploreMPC_bmy0, exploreMPC_dvaff0);
exploreMPC_LA_DENSE_DIAGZERO_MTVM2_4_8_4(exploreMPC_C0, exploreMPC_dvaff1, exploreMPC_D0, exploreMPC_dvaff0, exploreMPC_grad_eq0);
exploreMPC_LA_DENSE_DIAGZERO_MTVM2_4_8_4(exploreMPC_C0, exploreMPC_dvaff2, exploreMPC_D1, exploreMPC_dvaff1, exploreMPC_grad_eq1);
exploreMPC_LA_DENSE_DIAGZERO_MTVM2_4_8_4(exploreMPC_C0, exploreMPC_dvaff3, exploreMPC_D1, exploreMPC_dvaff2, exploreMPC_grad_eq2);
exploreMPC_LA_DENSE_DIAGZERO_MTVM2_4_8_4(exploreMPC_C0, exploreMPC_dvaff4, exploreMPC_D1, exploreMPC_dvaff3, exploreMPC_grad_eq3);
exploreMPC_LA_DENSE_DIAGZERO_MTVM2_4_8_4(exploreMPC_C0, exploreMPC_dvaff5, exploreMPC_D1, exploreMPC_dvaff4, exploreMPC_grad_eq4);
exploreMPC_LA_DENSE_DIAGZERO_MTVM2_4_8_4(exploreMPC_C0, exploreMPC_dvaff6, exploreMPC_D1, exploreMPC_dvaff5, exploreMPC_grad_eq5);
exploreMPC_LA_DENSE_DIAGZERO_MTVM2_4_8_4(exploreMPC_C0, exploreMPC_dvaff7, exploreMPC_D1, exploreMPC_dvaff6, exploreMPC_grad_eq6);
exploreMPC_LA_DENSE_DIAGZERO_MTVM2_4_8_4(exploreMPC_C0, exploreMPC_dvaff8, exploreMPC_D1, exploreMPC_dvaff7, exploreMPC_grad_eq7);
exploreMPC_LA_DENSE_DIAGZERO_MTVM2_4_8_4(exploreMPC_C0, exploreMPC_dvaff9, exploreMPC_D1, exploreMPC_dvaff8, exploreMPC_grad_eq8);
exploreMPC_LA_DIAGZERO_MTVM_4_4(exploreMPC_D9, exploreMPC_dvaff9, exploreMPC_grad_eq9);
exploreMPC_LA_VSUB2_76(exploreMPC_rd, exploreMPC_grad_eq, exploreMPC_rd);
exploreMPC_LA_DIAG_FORWARDBACKWARDSUB_8(exploreMPC_Phi0, exploreMPC_rd0, exploreMPC_dzaff0);
exploreMPC_LA_DIAG_FORWARDBACKWARDSUB_8(exploreMPC_Phi1, exploreMPC_rd1, exploreMPC_dzaff1);
exploreMPC_LA_DIAG_FORWARDBACKWARDSUB_8(exploreMPC_Phi2, exploreMPC_rd2, exploreMPC_dzaff2);
exploreMPC_LA_DIAG_FORWARDBACKWARDSUB_8(exploreMPC_Phi3, exploreMPC_rd3, exploreMPC_dzaff3);
exploreMPC_LA_DIAG_FORWARDBACKWARDSUB_8(exploreMPC_Phi4, exploreMPC_rd4, exploreMPC_dzaff4);
exploreMPC_LA_DIAG_FORWARDBACKWARDSUB_8(exploreMPC_Phi5, exploreMPC_rd5, exploreMPC_dzaff5);
exploreMPC_LA_DIAG_FORWARDBACKWARDSUB_8(exploreMPC_Phi6, exploreMPC_rd6, exploreMPC_dzaff6);
exploreMPC_LA_DIAG_FORWARDBACKWARDSUB_8(exploreMPC_Phi7, exploreMPC_rd7, exploreMPC_dzaff7);
exploreMPC_LA_DIAG_FORWARDBACKWARDSUB_8(exploreMPC_Phi8, exploreMPC_rd8, exploreMPC_dzaff8);
exploreMPC_LA_DIAG_FORWARDBACKWARDSUB_4(exploreMPC_Phi9, exploreMPC_rd9, exploreMPC_dzaff9);
exploreMPC_LA_VSUB_INDEXED_8(exploreMPC_dzaff0, exploreMPC_lbIdx0, exploreMPC_rilb0, exploreMPC_dslbaff0);
exploreMPC_LA_VSUB3_8(exploreMPC_llbbyslb0, exploreMPC_dslbaff0, exploreMPC_llb0, exploreMPC_dllbaff0);
exploreMPC_LA_VSUB2_INDEXED_8(exploreMPC_riub0, exploreMPC_dzaff0, exploreMPC_ubIdx0, exploreMPC_dsubaff0);
exploreMPC_LA_VSUB3_8(exploreMPC_lubbysub0, exploreMPC_dsubaff0, exploreMPC_lub0, exploreMPC_dlubaff0);
exploreMPC_LA_VSUB_INDEXED_8(exploreMPC_dzaff1, exploreMPC_lbIdx1, exploreMPC_rilb1, exploreMPC_dslbaff1);
exploreMPC_LA_VSUB3_8(exploreMPC_llbbyslb1, exploreMPC_dslbaff1, exploreMPC_llb1, exploreMPC_dllbaff1);
exploreMPC_LA_VSUB2_INDEXED_8(exploreMPC_riub1, exploreMPC_dzaff1, exploreMPC_ubIdx1, exploreMPC_dsubaff1);
exploreMPC_LA_VSUB3_8(exploreMPC_lubbysub1, exploreMPC_dsubaff1, exploreMPC_lub1, exploreMPC_dlubaff1);
exploreMPC_LA_VSUB_INDEXED_8(exploreMPC_dzaff2, exploreMPC_lbIdx2, exploreMPC_rilb2, exploreMPC_dslbaff2);
exploreMPC_LA_VSUB3_8(exploreMPC_llbbyslb2, exploreMPC_dslbaff2, exploreMPC_llb2, exploreMPC_dllbaff2);
exploreMPC_LA_VSUB2_INDEXED_8(exploreMPC_riub2, exploreMPC_dzaff2, exploreMPC_ubIdx2, exploreMPC_dsubaff2);
exploreMPC_LA_VSUB3_8(exploreMPC_lubbysub2, exploreMPC_dsubaff2, exploreMPC_lub2, exploreMPC_dlubaff2);
exploreMPC_LA_VSUB_INDEXED_8(exploreMPC_dzaff3, exploreMPC_lbIdx3, exploreMPC_rilb3, exploreMPC_dslbaff3);
exploreMPC_LA_VSUB3_8(exploreMPC_llbbyslb3, exploreMPC_dslbaff3, exploreMPC_llb3, exploreMPC_dllbaff3);
exploreMPC_LA_VSUB2_INDEXED_8(exploreMPC_riub3, exploreMPC_dzaff3, exploreMPC_ubIdx3, exploreMPC_dsubaff3);
exploreMPC_LA_VSUB3_8(exploreMPC_lubbysub3, exploreMPC_dsubaff3, exploreMPC_lub3, exploreMPC_dlubaff3);
exploreMPC_LA_VSUB_INDEXED_8(exploreMPC_dzaff4, exploreMPC_lbIdx4, exploreMPC_rilb4, exploreMPC_dslbaff4);
exploreMPC_LA_VSUB3_8(exploreMPC_llbbyslb4, exploreMPC_dslbaff4, exploreMPC_llb4, exploreMPC_dllbaff4);
exploreMPC_LA_VSUB2_INDEXED_8(exploreMPC_riub4, exploreMPC_dzaff4, exploreMPC_ubIdx4, exploreMPC_dsubaff4);
exploreMPC_LA_VSUB3_8(exploreMPC_lubbysub4, exploreMPC_dsubaff4, exploreMPC_lub4, exploreMPC_dlubaff4);
exploreMPC_LA_VSUB_INDEXED_8(exploreMPC_dzaff5, exploreMPC_lbIdx5, exploreMPC_rilb5, exploreMPC_dslbaff5);
exploreMPC_LA_VSUB3_8(exploreMPC_llbbyslb5, exploreMPC_dslbaff5, exploreMPC_llb5, exploreMPC_dllbaff5);
exploreMPC_LA_VSUB2_INDEXED_8(exploreMPC_riub5, exploreMPC_dzaff5, exploreMPC_ubIdx5, exploreMPC_dsubaff5);
exploreMPC_LA_VSUB3_8(exploreMPC_lubbysub5, exploreMPC_dsubaff5, exploreMPC_lub5, exploreMPC_dlubaff5);
exploreMPC_LA_VSUB_INDEXED_8(exploreMPC_dzaff6, exploreMPC_lbIdx6, exploreMPC_rilb6, exploreMPC_dslbaff6);
exploreMPC_LA_VSUB3_8(exploreMPC_llbbyslb6, exploreMPC_dslbaff6, exploreMPC_llb6, exploreMPC_dllbaff6);
exploreMPC_LA_VSUB2_INDEXED_8(exploreMPC_riub6, exploreMPC_dzaff6, exploreMPC_ubIdx6, exploreMPC_dsubaff6);
exploreMPC_LA_VSUB3_8(exploreMPC_lubbysub6, exploreMPC_dsubaff6, exploreMPC_lub6, exploreMPC_dlubaff6);
exploreMPC_LA_VSUB_INDEXED_8(exploreMPC_dzaff7, exploreMPC_lbIdx7, exploreMPC_rilb7, exploreMPC_dslbaff7);
exploreMPC_LA_VSUB3_8(exploreMPC_llbbyslb7, exploreMPC_dslbaff7, exploreMPC_llb7, exploreMPC_dllbaff7);
exploreMPC_LA_VSUB2_INDEXED_8(exploreMPC_riub7, exploreMPC_dzaff7, exploreMPC_ubIdx7, exploreMPC_dsubaff7);
exploreMPC_LA_VSUB3_8(exploreMPC_lubbysub7, exploreMPC_dsubaff7, exploreMPC_lub7, exploreMPC_dlubaff7);
exploreMPC_LA_VSUB_INDEXED_8(exploreMPC_dzaff8, exploreMPC_lbIdx8, exploreMPC_rilb8, exploreMPC_dslbaff8);
exploreMPC_LA_VSUB3_8(exploreMPC_llbbyslb8, exploreMPC_dslbaff8, exploreMPC_llb8, exploreMPC_dllbaff8);
exploreMPC_LA_VSUB2_INDEXED_8(exploreMPC_riub8, exploreMPC_dzaff8, exploreMPC_ubIdx8, exploreMPC_dsubaff8);
exploreMPC_LA_VSUB3_8(exploreMPC_lubbysub8, exploreMPC_dsubaff8, exploreMPC_lub8, exploreMPC_dlubaff8);
exploreMPC_LA_VSUB_INDEXED_4(exploreMPC_dzaff9, exploreMPC_lbIdx9, exploreMPC_rilb9, exploreMPC_dslbaff9);
exploreMPC_LA_VSUB3_4(exploreMPC_llbbyslb9, exploreMPC_dslbaff9, exploreMPC_llb9, exploreMPC_dllbaff9);
exploreMPC_LA_VSUB2_INDEXED_4(exploreMPC_riub9, exploreMPC_dzaff9, exploreMPC_ubIdx9, exploreMPC_dsubaff9);
exploreMPC_LA_VSUB3_4(exploreMPC_lubbysub9, exploreMPC_dsubaff9, exploreMPC_lub9, exploreMPC_dlubaff9);
info->lsit_aff = exploreMPC_LINESEARCH_BACKTRACKING_AFFINE(exploreMPC_l, exploreMPC_s, exploreMPC_dl_aff, exploreMPC_ds_aff, &info->step_aff, &info->mu_aff);
if( info->lsit_aff == exploreMPC_NOPROGRESS ){
exitcode = exploreMPC_NOPROGRESS; break;
}
sigma_3rdroot = info->mu_aff / info->mu;
info->sigma = sigma_3rdroot*sigma_3rdroot*sigma_3rdroot;
musigma = info->mu * info->sigma;
exploreMPC_LA_VSUB5_152(exploreMPC_ds_aff, exploreMPC_dl_aff, musigma, exploreMPC_ccrhs);
exploreMPC_LA_VSUB6_INDEXED_8_8_8(exploreMPC_ccrhsub0, exploreMPC_sub0, exploreMPC_ubIdx0, exploreMPC_ccrhsl0, exploreMPC_slb0, exploreMPC_lbIdx0, exploreMPC_rd0);
exploreMPC_LA_VSUB6_INDEXED_8_8_8(exploreMPC_ccrhsub1, exploreMPC_sub1, exploreMPC_ubIdx1, exploreMPC_ccrhsl1, exploreMPC_slb1, exploreMPC_lbIdx1, exploreMPC_rd1);
exploreMPC_LA_DIAG_FORWARDSUB_8(exploreMPC_Phi0, exploreMPC_rd0, exploreMPC_Lbyrd0);
exploreMPC_LA_DIAG_FORWARDSUB_8(exploreMPC_Phi1, exploreMPC_rd1, exploreMPC_Lbyrd1);
exploreMPC_LA_DIAGZERO_MVM_4(exploreMPC_W0, exploreMPC_Lbyrd0, exploreMPC_beta0);
exploreMPC_LA_DENSE_FORWARDSUB_4(exploreMPC_Ld0, exploreMPC_beta0, exploreMPC_yy0);
exploreMPC_LA_DENSE_DIAGZERO_2MVMADD_4_8_8(exploreMPC_V0, exploreMPC_Lbyrd0, exploreMPC_W1, exploreMPC_Lbyrd1, exploreMPC_beta1);
exploreMPC_LA_DENSE_MVMSUB1_4_4(exploreMPC_Lsd1, exploreMPC_yy0, exploreMPC_beta1, exploreMPC_bmy1);
exploreMPC_LA_DENSE_FORWARDSUB_4(exploreMPC_Ld1, exploreMPC_bmy1, exploreMPC_yy1);
exploreMPC_LA_VSUB6_INDEXED_8_8_8(exploreMPC_ccrhsub2, exploreMPC_sub2, exploreMPC_ubIdx2, exploreMPC_ccrhsl2, exploreMPC_slb2, exploreMPC_lbIdx2, exploreMPC_rd2);
exploreMPC_LA_DIAG_FORWARDSUB_8(exploreMPC_Phi2, exploreMPC_rd2, exploreMPC_Lbyrd2);
exploreMPC_LA_DENSE_DIAGZERO_2MVMADD_4_8_8(exploreMPC_V1, exploreMPC_Lbyrd1, exploreMPC_W2, exploreMPC_Lbyrd2, exploreMPC_beta2);
exploreMPC_LA_DENSE_MVMSUB1_4_4(exploreMPC_Lsd2, exploreMPC_yy1, exploreMPC_beta2, exploreMPC_bmy2);
exploreMPC_LA_DENSE_FORWARDSUB_4(exploreMPC_Ld2, exploreMPC_bmy2, exploreMPC_yy2);
exploreMPC_LA_VSUB6_INDEXED_8_8_8(exploreMPC_ccrhsub3, exploreMPC_sub3, exploreMPC_ubIdx3, exploreMPC_ccrhsl3, exploreMPC_slb3, exploreMPC_lbIdx3, exploreMPC_rd3);
exploreMPC_LA_DIAG_FORWARDSUB_8(exploreMPC_Phi3, exploreMPC_rd3, exploreMPC_Lbyrd3);
exploreMPC_LA_DENSE_DIAGZERO_2MVMADD_4_8_8(exploreMPC_V2, exploreMPC_Lbyrd2, exploreMPC_W3, exploreMPC_Lbyrd3, exploreMPC_beta3);
exploreMPC_LA_DENSE_MVMSUB1_4_4(exploreMPC_Lsd3, exploreMPC_yy2, exploreMPC_beta3, exploreMPC_bmy3);
exploreMPC_LA_DENSE_FORWARDSUB_4(exploreMPC_Ld3, exploreMPC_bmy3, exploreMPC_yy3);
exploreMPC_LA_VSUB6_INDEXED_8_8_8(exploreMPC_ccrhsub4, exploreMPC_sub4, exploreMPC_ubIdx4, exploreMPC_ccrhsl4, exploreMPC_slb4, exploreMPC_lbIdx4, exploreMPC_rd4);
exploreMPC_LA_DIAG_FORWARDSUB_8(exploreMPC_Phi4, exploreMPC_rd4, exploreMPC_Lbyrd4);
exploreMPC_LA_DENSE_DIAGZERO_2MVMADD_4_8_8(exploreMPC_V3, exploreMPC_Lbyrd3, exploreMPC_W4, exploreMPC_Lbyrd4, exploreMPC_beta4);
exploreMPC_LA_DENSE_MVMSUB1_4_4(exploreMPC_Lsd4, exploreMPC_yy3, exploreMPC_beta4, exploreMPC_bmy4);
exploreMPC_LA_DENSE_FORWARDSUB_4(exploreMPC_Ld4, exploreMPC_bmy4, exploreMPC_yy4);
exploreMPC_LA_VSUB6_INDEXED_8_8_8(exploreMPC_ccrhsub5, exploreMPC_sub5, exploreMPC_ubIdx5, exploreMPC_ccrhsl5, exploreMPC_slb5, exploreMPC_lbIdx5, exploreMPC_rd5);
exploreMPC_LA_DIAG_FORWARDSUB_8(exploreMPC_Phi5, exploreMPC_rd5, exploreMPC_Lbyrd5);
exploreMPC_LA_DENSE_DIAGZERO_2MVMADD_4_8_8(exploreMPC_V4, exploreMPC_Lbyrd4, exploreMPC_W5, exploreMPC_Lbyrd5, exploreMPC_beta5);
exploreMPC_LA_DENSE_MVMSUB1_4_4(exploreMPC_Lsd5, exploreMPC_yy4, exploreMPC_beta5, exploreMPC_bmy5);
exploreMPC_LA_DENSE_FORWARDSUB_4(exploreMPC_Ld5, exploreMPC_bmy5, exploreMPC_yy5);
exploreMPC_LA_VSUB6_INDEXED_8_8_8(exploreMPC_ccrhsub6, exploreMPC_sub6, exploreMPC_ubIdx6, exploreMPC_ccrhsl6, exploreMPC_slb6, exploreMPC_lbIdx6, exploreMPC_rd6);
exploreMPC_LA_DIAG_FORWARDSUB_8(exploreMPC_Phi6, exploreMPC_rd6, exploreMPC_Lbyrd6);
exploreMPC_LA_DENSE_DIAGZERO_2MVMADD_4_8_8(exploreMPC_V5, exploreMPC_Lbyrd5, exploreMPC_W6, exploreMPC_Lbyrd6, exploreMPC_beta6);
exploreMPC_LA_DENSE_MVMSUB1_4_4(exploreMPC_Lsd6, exploreMPC_yy5, exploreMPC_beta6, exploreMPC_bmy6);
exploreMPC_LA_DENSE_FORWARDSUB_4(exploreMPC_Ld6, exploreMPC_bmy6, exploreMPC_yy6);
exploreMPC_LA_VSUB6_INDEXED_8_8_8(exploreMPC_ccrhsub7, exploreMPC_sub7, exploreMPC_ubIdx7, exploreMPC_ccrhsl7, exploreMPC_slb7, exploreMPC_lbIdx7, exploreMPC_rd7);
exploreMPC_LA_DIAG_FORWARDSUB_8(exploreMPC_Phi7, exploreMPC_rd7, exploreMPC_Lbyrd7);
exploreMPC_LA_DENSE_DIAGZERO_2MVMADD_4_8_8(exploreMPC_V6, exploreMPC_Lbyrd6, exploreMPC_W7, exploreMPC_Lbyrd7, exploreMPC_beta7);
exploreMPC_LA_DENSE_MVMSUB1_4_4(exploreMPC_Lsd7, exploreMPC_yy6, exploreMPC_beta7, exploreMPC_bmy7);
exploreMPC_LA_DENSE_FORWARDSUB_4(exploreMPC_Ld7, exploreMPC_bmy7, exploreMPC_yy7);
exploreMPC_LA_VSUB6_INDEXED_8_8_8(exploreMPC_ccrhsub8, exploreMPC_sub8, exploreMPC_ubIdx8, exploreMPC_ccrhsl8, exploreMPC_slb8, exploreMPC_lbIdx8, exploreMPC_rd8);
exploreMPC_LA_DIAG_FORWARDSUB_8(exploreMPC_Phi8, exploreMPC_rd8, exploreMPC_Lbyrd8);
exploreMPC_LA_DENSE_DIAGZERO_2MVMADD_4_8_8(exploreMPC_V7, exploreMPC_Lbyrd7, exploreMPC_W8, exploreMPC_Lbyrd8, exploreMPC_beta8);
exploreMPC_LA_DENSE_MVMSUB1_4_4(exploreMPC_Lsd8, exploreMPC_yy7, exploreMPC_beta8, exploreMPC_bmy8);
exploreMPC_LA_DENSE_FORWARDSUB_4(exploreMPC_Ld8, exploreMPC_bmy8, exploreMPC_yy8);
exploreMPC_LA_VSUB6_INDEXED_4_4_4(exploreMPC_ccrhsub9, exploreMPC_sub9, exploreMPC_ubIdx9, exploreMPC_ccrhsl9, exploreMPC_slb9, exploreMPC_lbIdx9, exploreMPC_rd9);
exploreMPC_LA_DIAG_FORWARDSUB_4(exploreMPC_Phi9, exploreMPC_rd9, exploreMPC_Lbyrd9);
exploreMPC_LA_DENSE_DIAGZERO_2MVMADD_4_8_4(exploreMPC_V8, exploreMPC_Lbyrd8, exploreMPC_W9, exploreMPC_Lbyrd9, exploreMPC_beta9);
exploreMPC_LA_DENSE_MVMSUB1_4_4(exploreMPC_Lsd9, exploreMPC_yy8, exploreMPC_beta9, exploreMPC_bmy9);
exploreMPC_LA_DENSE_FORWARDSUB_4(exploreMPC_Ld9, exploreMPC_bmy9, exploreMPC_yy9);
exploreMPC_LA_DENSE_BACKWARDSUB_4(exploreMPC_Ld9, exploreMPC_yy9, exploreMPC_dvcc9);
exploreMPC_LA_DENSE_MTVMSUB_4_4(exploreMPC_Lsd9, exploreMPC_dvcc9, exploreMPC_yy8, exploreMPC_bmy8);
exploreMPC_LA_DENSE_BACKWARDSUB_4(exploreMPC_Ld8, exploreMPC_bmy8, exploreMPC_dvcc8);
exploreMPC_LA_DENSE_MTVMSUB_4_4(exploreMPC_Lsd8, exploreMPC_dvcc8, exploreMPC_yy7, exploreMPC_bmy7);
exploreMPC_LA_DENSE_BACKWARDSUB_4(exploreMPC_Ld7, exploreMPC_bmy7, exploreMPC_dvcc7);
exploreMPC_LA_DENSE_MTVMSUB_4_4(exploreMPC_Lsd7, exploreMPC_dvcc7, exploreMPC_yy6, exploreMPC_bmy6);
exploreMPC_LA_DENSE_BACKWARDSUB_4(exploreMPC_Ld6, exploreMPC_bmy6, exploreMPC_dvcc6);
exploreMPC_LA_DENSE_MTVMSUB_4_4(exploreMPC_Lsd6, exploreMPC_dvcc6, exploreMPC_yy5, exploreMPC_bmy5);
exploreMPC_LA_DENSE_BACKWARDSUB_4(exploreMPC_Ld5, exploreMPC_bmy5, exploreMPC_dvcc5);
exploreMPC_LA_DENSE_MTVMSUB_4_4(exploreMPC_Lsd5, exploreMPC_dvcc5, exploreMPC_yy4, exploreMPC_bmy4);
exploreMPC_LA_DENSE_BACKWARDSUB_4(exploreMPC_Ld4, exploreMPC_bmy4, exploreMPC_dvcc4);
exploreMPC_LA_DENSE_MTVMSUB_4_4(exploreMPC_Lsd4, exploreMPC_dvcc4, exploreMPC_yy3, exploreMPC_bmy3);
exploreMPC_LA_DENSE_BACKWARDSUB_4(exploreMPC_Ld3, exploreMPC_bmy3, exploreMPC_dvcc3);
exploreMPC_LA_DENSE_MTVMSUB_4_4(exploreMPC_Lsd3, exploreMPC_dvcc3, exploreMPC_yy2, exploreMPC_bmy2);
exploreMPC_LA_DENSE_BACKWARDSUB_4(exploreMPC_Ld2, exploreMPC_bmy2, exploreMPC_dvcc2);
exploreMPC_LA_DENSE_MTVMSUB_4_4(exploreMPC_Lsd2, exploreMPC_dvcc2, exploreMPC_yy1, exploreMPC_bmy1);
exploreMPC_LA_DENSE_BACKWARDSUB_4(exploreMPC_Ld1, exploreMPC_bmy1, exploreMPC_dvcc1);
exploreMPC_LA_DENSE_MTVMSUB_4_4(exploreMPC_Lsd1, exploreMPC_dvcc1, exploreMPC_yy0, exploreMPC_bmy0);
exploreMPC_LA_DENSE_BACKWARDSUB_4(exploreMPC_Ld0, exploreMPC_bmy0, exploreMPC_dvcc0);
exploreMPC_LA_DENSE_DIAGZERO_MTVM2_4_8_4(exploreMPC_C0, exploreMPC_dvcc1, exploreMPC_D0, exploreMPC_dvcc0, exploreMPC_grad_eq0);
exploreMPC_LA_DENSE_DIAGZERO_MTVM2_4_8_4(exploreMPC_C0, exploreMPC_dvcc2, exploreMPC_D1, exploreMPC_dvcc1, exploreMPC_grad_eq1);
exploreMPC_LA_DENSE_DIAGZERO_MTVM2_4_8_4(exploreMPC_C0, exploreMPC_dvcc3, exploreMPC_D1, exploreMPC_dvcc2, exploreMPC_grad_eq2);
exploreMPC_LA_DENSE_DIAGZERO_MTVM2_4_8_4(exploreMPC_C0, exploreMPC_dvcc4, exploreMPC_D1, exploreMPC_dvcc3, exploreMPC_grad_eq3);
exploreMPC_LA_DENSE_DIAGZERO_MTVM2_4_8_4(exploreMPC_C0, exploreMPC_dvcc5, exploreMPC_D1, exploreMPC_dvcc4, exploreMPC_grad_eq4);
exploreMPC_LA_DENSE_DIAGZERO_MTVM2_4_8_4(exploreMPC_C0, exploreMPC_dvcc6, exploreMPC_D1, exploreMPC_dvcc5, exploreMPC_grad_eq5);
exploreMPC_LA_DENSE_DIAGZERO_MTVM2_4_8_4(exploreMPC_C0, exploreMPC_dvcc7, exploreMPC_D1, exploreMPC_dvcc6, exploreMPC_grad_eq6);
exploreMPC_LA_DENSE_DIAGZERO_MTVM2_4_8_4(exploreMPC_C0, exploreMPC_dvcc8, exploreMPC_D1, exploreMPC_dvcc7, exploreMPC_grad_eq7);
exploreMPC_LA_DENSE_DIAGZERO_MTVM2_4_8_4(exploreMPC_C0, exploreMPC_dvcc9, exploreMPC_D1, exploreMPC_dvcc8, exploreMPC_grad_eq8);
exploreMPC_LA_DIAGZERO_MTVM_4_4(exploreMPC_D9, exploreMPC_dvcc9, exploreMPC_grad_eq9);
exploreMPC_LA_VSUB_76(exploreMPC_rd, exploreMPC_grad_eq, exploreMPC_rd);
exploreMPC_LA_DIAG_FORWARDBACKWARDSUB_8(exploreMPC_Phi0, exploreMPC_rd0, exploreMPC_dzcc0);
exploreMPC_LA_DIAG_FORWARDBACKWARDSUB_8(exploreMPC_Phi1, exploreMPC_rd1, exploreMPC_dzcc1);
exploreMPC_LA_DIAG_FORWARDBACKWARDSUB_8(exploreMPC_Phi2, exploreMPC_rd2, exploreMPC_dzcc2);
exploreMPC_LA_DIAG_FORWARDBACKWARDSUB_8(exploreMPC_Phi3, exploreMPC_rd3, exploreMPC_dzcc3);
exploreMPC_LA_DIAG_FORWARDBACKWARDSUB_8(exploreMPC_Phi4, exploreMPC_rd4, exploreMPC_dzcc4);
exploreMPC_LA_DIAG_FORWARDBACKWARDSUB_8(exploreMPC_Phi5, exploreMPC_rd5, exploreMPC_dzcc5);
exploreMPC_LA_DIAG_FORWARDBACKWARDSUB_8(exploreMPC_Phi6, exploreMPC_rd6, exploreMPC_dzcc6);
exploreMPC_LA_DIAG_FORWARDBACKWARDSUB_8(exploreMPC_Phi7, exploreMPC_rd7, exploreMPC_dzcc7);
exploreMPC_LA_DIAG_FORWARDBACKWARDSUB_8(exploreMPC_Phi8, exploreMPC_rd8, exploreMPC_dzcc8);
exploreMPC_LA_DIAG_FORWARDBACKWARDSUB_4(exploreMPC_Phi9, exploreMPC_rd9, exploreMPC_dzcc9);
exploreMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_8(exploreMPC_ccrhsl0, exploreMPC_slb0, exploreMPC_llbbyslb0, exploreMPC_dzcc0, exploreMPC_lbIdx0, exploreMPC_dllbcc0);
exploreMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_8(exploreMPC_ccrhsub0, exploreMPC_sub0, exploreMPC_lubbysub0, exploreMPC_dzcc0, exploreMPC_ubIdx0, exploreMPC_dlubcc0);
exploreMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_8(exploreMPC_ccrhsl1, exploreMPC_slb1, exploreMPC_llbbyslb1, exploreMPC_dzcc1, exploreMPC_lbIdx1, exploreMPC_dllbcc1);
exploreMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_8(exploreMPC_ccrhsub1, exploreMPC_sub1, exploreMPC_lubbysub1, exploreMPC_dzcc1, exploreMPC_ubIdx1, exploreMPC_dlubcc1);
exploreMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_8(exploreMPC_ccrhsl2, exploreMPC_slb2, exploreMPC_llbbyslb2, exploreMPC_dzcc2, exploreMPC_lbIdx2, exploreMPC_dllbcc2);
exploreMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_8(exploreMPC_ccrhsub2, exploreMPC_sub2, exploreMPC_lubbysub2, exploreMPC_dzcc2, exploreMPC_ubIdx2, exploreMPC_dlubcc2);
exploreMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_8(exploreMPC_ccrhsl3, exploreMPC_slb3, exploreMPC_llbbyslb3, exploreMPC_dzcc3, exploreMPC_lbIdx3, exploreMPC_dllbcc3);
exploreMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_8(exploreMPC_ccrhsub3, exploreMPC_sub3, exploreMPC_lubbysub3, exploreMPC_dzcc3, exploreMPC_ubIdx3, exploreMPC_dlubcc3);
exploreMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_8(exploreMPC_ccrhsl4, exploreMPC_slb4, exploreMPC_llbbyslb4, exploreMPC_dzcc4, exploreMPC_lbIdx4, exploreMPC_dllbcc4);
exploreMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_8(exploreMPC_ccrhsub4, exploreMPC_sub4, exploreMPC_lubbysub4, exploreMPC_dzcc4, exploreMPC_ubIdx4, exploreMPC_dlubcc4);
exploreMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_8(exploreMPC_ccrhsl5, exploreMPC_slb5, exploreMPC_llbbyslb5, exploreMPC_dzcc5, exploreMPC_lbIdx5, exploreMPC_dllbcc5);
exploreMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_8(exploreMPC_ccrhsub5, exploreMPC_sub5, exploreMPC_lubbysub5, exploreMPC_dzcc5, exploreMPC_ubIdx5, exploreMPC_dlubcc5);
exploreMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_8(exploreMPC_ccrhsl6, exploreMPC_slb6, exploreMPC_llbbyslb6, exploreMPC_dzcc6, exploreMPC_lbIdx6, exploreMPC_dllbcc6);
exploreMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_8(exploreMPC_ccrhsub6, exploreMPC_sub6, exploreMPC_lubbysub6, exploreMPC_dzcc6, exploreMPC_ubIdx6, exploreMPC_dlubcc6);
exploreMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_8(exploreMPC_ccrhsl7, exploreMPC_slb7, exploreMPC_llbbyslb7, exploreMPC_dzcc7, exploreMPC_lbIdx7, exploreMPC_dllbcc7);
exploreMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_8(exploreMPC_ccrhsub7, exploreMPC_sub7, exploreMPC_lubbysub7, exploreMPC_dzcc7, exploreMPC_ubIdx7, exploreMPC_dlubcc7);
exploreMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_8(exploreMPC_ccrhsl8, exploreMPC_slb8, exploreMPC_llbbyslb8, exploreMPC_dzcc8, exploreMPC_lbIdx8, exploreMPC_dllbcc8);
exploreMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_8(exploreMPC_ccrhsub8, exploreMPC_sub8, exploreMPC_lubbysub8, exploreMPC_dzcc8, exploreMPC_ubIdx8, exploreMPC_dlubcc8);
exploreMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_4(exploreMPC_ccrhsl9, exploreMPC_slb9, exploreMPC_llbbyslb9, exploreMPC_dzcc9, exploreMPC_lbIdx9, exploreMPC_dllbcc9);
exploreMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_4(exploreMPC_ccrhsub9, exploreMPC_sub9, exploreMPC_lubbysub9, exploreMPC_dzcc9, exploreMPC_ubIdx9, exploreMPC_dlubcc9);
exploreMPC_LA_VSUB7_152(exploreMPC_l, exploreMPC_ccrhs, exploreMPC_s, exploreMPC_dl_cc, exploreMPC_ds_cc);
exploreMPC_LA_VADD_76(exploreMPC_dz_cc, exploreMPC_dz_aff);
exploreMPC_LA_VADD_40(exploreMPC_dv_cc, exploreMPC_dv_aff);
exploreMPC_LA_VADD_152(exploreMPC_dl_cc, exploreMPC_dl_aff);
exploreMPC_LA_VADD_152(exploreMPC_ds_cc, exploreMPC_ds_aff);
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
output->z3[4] = exploreMPC_z2[4];
output->z3[5] = exploreMPC_z2[5];
output->z3[6] = exploreMPC_z2[6];
output->z3[7] = exploreMPC_z2[7];
output->z4[0] = exploreMPC_z3[0];
output->z4[1] = exploreMPC_z3[1];
output->z4[2] = exploreMPC_z3[2];
output->z4[3] = exploreMPC_z3[3];
output->z4[4] = exploreMPC_z3[4];
output->z4[5] = exploreMPC_z3[5];
output->z4[6] = exploreMPC_z3[6];
output->z4[7] = exploreMPC_z3[7];
output->z5[0] = exploreMPC_z4[0];
output->z5[1] = exploreMPC_z4[1];
output->z5[2] = exploreMPC_z4[2];
output->z5[3] = exploreMPC_z4[3];
output->z5[4] = exploreMPC_z4[4];
output->z5[5] = exploreMPC_z4[5];
output->z5[6] = exploreMPC_z4[6];
output->z5[7] = exploreMPC_z4[7];
output->z6[0] = exploreMPC_z5[0];
output->z6[1] = exploreMPC_z5[1];
output->z6[2] = exploreMPC_z5[2];
output->z6[3] = exploreMPC_z5[3];
output->z6[4] = exploreMPC_z5[4];
output->z6[5] = exploreMPC_z5[5];
output->z6[6] = exploreMPC_z5[6];
output->z6[7] = exploreMPC_z5[7];
output->z7[0] = exploreMPC_z6[0];
output->z7[1] = exploreMPC_z6[1];
output->z7[2] = exploreMPC_z6[2];
output->z7[3] = exploreMPC_z6[3];
output->z7[4] = exploreMPC_z6[4];
output->z7[5] = exploreMPC_z6[5];
output->z7[6] = exploreMPC_z6[6];
output->z7[7] = exploreMPC_z6[7];
output->z8[0] = exploreMPC_z7[0];
output->z8[1] = exploreMPC_z7[1];
output->z8[2] = exploreMPC_z7[2];
output->z8[3] = exploreMPC_z7[3];
output->z8[4] = exploreMPC_z7[4];
output->z8[5] = exploreMPC_z7[5];
output->z8[6] = exploreMPC_z7[6];
output->z8[7] = exploreMPC_z7[7];
output->z9[0] = exploreMPC_z8[0];
output->z9[1] = exploreMPC_z8[1];
output->z9[2] = exploreMPC_z8[2];
output->z9[3] = exploreMPC_z8[3];
output->z9[4] = exploreMPC_z8[4];
output->z9[5] = exploreMPC_z8[5];
output->z9[6] = exploreMPC_z8[6];
output->z9[7] = exploreMPC_z8[7];
output->z10[0] = exploreMPC_z9[0];
output->z10[1] = exploreMPC_z9[1];
output->z10[2] = exploreMPC_z9[2];
output->z10[3] = exploreMPC_z9[3];

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
