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

#include "plattMPC.h"

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
 * Initializes a vector of length 118 with a value.
 */
void plattMPC_LA_INITIALIZEVECTOR_118(plattMPC_FLOAT* vec, plattMPC_FLOAT value)
{
	int i;
	for( i=0; i<118; i++ )
	{
		vec[i] = value;
	}
}


/*
 * Initializes a vector of length 90 with a value.
 */
void plattMPC_LA_INITIALIZEVECTOR_90(plattMPC_FLOAT* vec, plattMPC_FLOAT value)
{
	int i;
	for( i=0; i<90; i++ )
	{
		vec[i] = value;
	}
}


/*
 * Initializes a vector of length 236 with a value.
 */
void plattMPC_LA_INITIALIZEVECTOR_236(plattMPC_FLOAT* vec, plattMPC_FLOAT value)
{
	int i;
	for( i=0; i<236; i++ )
	{
		vec[i] = value;
	}
}


/* 
 * Calculates a dot product and adds it to a variable: z += x'*y; 
 * This function is for vectors of length 236.
 */
void plattMPC_LA_DOTACC_236(plattMPC_FLOAT *x, plattMPC_FLOAT *y, plattMPC_FLOAT *z)
{
	int i;
	for( i=0; i<236; i++ ){
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
void plattMPC_LA_DIAG_QUADFCN_8(plattMPC_FLOAT* H, plattMPC_FLOAT* f, plattMPC_FLOAT* z, plattMPC_FLOAT* grad, plattMPC_FLOAT* value)
{
	int i;
	plattMPC_FLOAT hz;	
	for( i=0; i<8; i++){
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
void plattMPC_LA_DIAG_QUADFCN_6(plattMPC_FLOAT* H, plattMPC_FLOAT* f, plattMPC_FLOAT* z, plattMPC_FLOAT* grad, plattMPC_FLOAT* value)
{
	int i;
	plattMPC_FLOAT hz;	
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
void plattMPC_LA_DIAGZERO_MVMSUB6_6(plattMPC_FLOAT *B, plattMPC_FLOAT *u, plattMPC_FLOAT *b, plattMPC_FLOAT *l, plattMPC_FLOAT *r, plattMPC_FLOAT *z, plattMPC_FLOAT *y)
{
	int i;
	plattMPC_FLOAT Bu[6];
	plattMPC_FLOAT norm = *y;
	plattMPC_FLOAT lr = 0;

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
void plattMPC_LA_DENSE_DIAGZERO_MVMSUB3_6_8_8(plattMPC_FLOAT *A, plattMPC_FLOAT *x, plattMPC_FLOAT *B, plattMPC_FLOAT *u, plattMPC_FLOAT *b, plattMPC_FLOAT *l, plattMPC_FLOAT *r, plattMPC_FLOAT *z, plattMPC_FLOAT *y)
{
	int i;
	int j;
	int k = 0;
	plattMPC_FLOAT AxBu[6];
	plattMPC_FLOAT norm = *y;
	plattMPC_FLOAT lr = 0;

	/* do A*x + B*u first */
	for( i=0; i<6; i++ ){
		AxBu[i] = A[k++]*x[0] + B[i]*u[i];
	}	

	for( j=1; j<8; j++ ){		
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
void plattMPC_LA_DENSE_DIAGZERO_MVMSUB3_6_8_6(plattMPC_FLOAT *A, plattMPC_FLOAT *x, plattMPC_FLOAT *B, plattMPC_FLOAT *u, plattMPC_FLOAT *b, plattMPC_FLOAT *l, plattMPC_FLOAT *r, plattMPC_FLOAT *z, plattMPC_FLOAT *y)
{
	int i;
	int j;
	int k = 0;
	plattMPC_FLOAT AxBu[6];
	plattMPC_FLOAT norm = *y;
	plattMPC_FLOAT lr = 0;

	/* do A*x + B*u first */
	for( i=0; i<6; i++ ){
		AxBu[i] = A[k++]*x[0] + B[i]*u[i];
	}	

	for( j=1; j<8; j++ ){		
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
 * where A is of size [6 x 8] and stored in column major format.
 * and B is of size [6 x 8] and stored in diagzero format
 * Note the transposes of A and B!
 */
void plattMPC_LA_DENSE_DIAGZERO_MTVM2_6_8_6(plattMPC_FLOAT *A, plattMPC_FLOAT *x, plattMPC_FLOAT *B, plattMPC_FLOAT *y, plattMPC_FLOAT *z)
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
	for( i=6 ;i<8; i++ ){
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
void plattMPC_LA_DIAGZERO_MTVM_6_6(plattMPC_FLOAT *M, plattMPC_FLOAT *x, plattMPC_FLOAT *y)
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
 * for vectors of length 8. Output z is of course scalar.
 */
void plattMPC_LA_VSUBADD3_8(plattMPC_FLOAT* t, plattMPC_FLOAT* u, int* uidx, plattMPC_FLOAT* v, plattMPC_FLOAT* w, plattMPC_FLOAT* y, plattMPC_FLOAT* z, plattMPC_FLOAT* r)
{
	int i;
	plattMPC_FLOAT norm = *r;
	plattMPC_FLOAT vx = 0;
	plattMPC_FLOAT x;
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
void plattMPC_LA_VSUBADD2_8(plattMPC_FLOAT* t, int* tidx, plattMPC_FLOAT* u, plattMPC_FLOAT* v, plattMPC_FLOAT* w, plattMPC_FLOAT* y, plattMPC_FLOAT* z, plattMPC_FLOAT* r)
{
	int i;
	plattMPC_FLOAT norm = *r;
	plattMPC_FLOAT vx = 0;
	plattMPC_FLOAT x;
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
 * for vectors of length 6. Output z is of course scalar.
 */
void plattMPC_LA_VSUBADD3_6(plattMPC_FLOAT* t, plattMPC_FLOAT* u, int* uidx, plattMPC_FLOAT* v, plattMPC_FLOAT* w, plattMPC_FLOAT* y, plattMPC_FLOAT* z, plattMPC_FLOAT* r)
{
	int i;
	plattMPC_FLOAT norm = *r;
	plattMPC_FLOAT vx = 0;
	plattMPC_FLOAT x;
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
void plattMPC_LA_VSUBADD2_6(plattMPC_FLOAT* t, int* tidx, plattMPC_FLOAT* u, plattMPC_FLOAT* v, plattMPC_FLOAT* w, plattMPC_FLOAT* y, plattMPC_FLOAT* z, plattMPC_FLOAT* r)
{
	int i;
	plattMPC_FLOAT norm = *r;
	plattMPC_FLOAT vx = 0;
	plattMPC_FLOAT x;
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
 * Special function for box constraints of length 8
 * Returns also L/S, a value that is often used elsewhere.
 */
void plattMPC_LA_INEQ_B_GRAD_8_8_8(plattMPC_FLOAT *lu, plattMPC_FLOAT *su, plattMPC_FLOAT *ru, plattMPC_FLOAT *ll, plattMPC_FLOAT *sl, plattMPC_FLOAT *rl, int* lbIdx, int* ubIdx, plattMPC_FLOAT *grad, plattMPC_FLOAT *lubysu, plattMPC_FLOAT *llbysl)
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
 * Special function for box constraints of length 6
 * Returns also L/S, a value that is often used elsewhere.
 */
void plattMPC_LA_INEQ_B_GRAD_6_6_6(plattMPC_FLOAT *lu, plattMPC_FLOAT *su, plattMPC_FLOAT *ru, plattMPC_FLOAT *ll, plattMPC_FLOAT *sl, plattMPC_FLOAT *rl, int* lbIdx, int* ubIdx, plattMPC_FLOAT *grad, plattMPC_FLOAT *lubysu, plattMPC_FLOAT *llbysl)
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
 * of length 118.
 */
void plattMPC_LA_VVADD3_118(plattMPC_FLOAT *u, plattMPC_FLOAT *v, plattMPC_FLOAT *w, plattMPC_FLOAT *z)
{
	int i;
	for( i=0; i<118; i++ ){
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
void plattMPC_LA_DIAG_CHOL_ONELOOP_LBUB_8_8_8(plattMPC_FLOAT *H, plattMPC_FLOAT *llbysl, int* lbIdx, plattMPC_FLOAT *lubysu, int* ubIdx, plattMPC_FLOAT *Phi)


{
	int i;
	
	/* compute cholesky */
	for( i=0; i<8; i++ ){
		Phi[i] = H[i] + llbysl[i] + lubysu[i];

#if plattMPC_SET_PRINTLEVEL > 0 && defined PRINTNUMERICALWARNINGS
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
 * where A is to be computed and is of size [6 x 8],
 * B is given and of size [6 x 8], L is a diagonal
 * matrix of size 6 stored in diagonal matrix 
 * storage format. Note the transpose of L has no impact!
 *
 * Result: A in column major storage format.
 *
 */
void plattMPC_LA_DIAG_MATRIXFORWARDSUB_6_8(plattMPC_FLOAT *L, plattMPC_FLOAT *B, plattMPC_FLOAT *A)
{
    int i,j;
	 int k = 0;

	for( j=0; j<8; j++){
		for( i=0; i<6; i++){
			A[k] = B[k]/L[j];
			k++;
		}
	}

}


/**
 * Forward substitution for the matrix equation A*L' = B
 * where A is to be computed and is of size [6 x 8],
 * B is given and of size [6 x 8], L is a diagonal
 *  matrix of size 8 stored in diagonal 
 * storage format. Note the transpose of L!
 *
 * Result: A in diagonalzero storage format.
 *
 */
void plattMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_6_8(plattMPC_FLOAT *L, plattMPC_FLOAT *B, plattMPC_FLOAT *A)
{
	int j;
    for( j=0; j<8; j++ ){   
		A[j] = B[j]/L[j];
     }
}


/**
 * Compute C = A*B' where 
 *
 *	size(A) = [6 x 8]
 *  size(B) = [6 x 8] in diagzero format
 * 
 * A and C matrices are stored in column major format.
 * 
 * 
 */
void plattMPC_LA_DENSE_DIAGZERO_MMTM_6_8_6(plattMPC_FLOAT *A, plattMPC_FLOAT *B, plattMPC_FLOAT *C)
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
 * The dimensions involved are 8.
 */
void plattMPC_LA_DIAG_FORWARDSUB_8(plattMPC_FLOAT *L, plattMPC_FLOAT *b, plattMPC_FLOAT *y)
{
    int i;

    for( i=0; i<8; i++ ){
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
void plattMPC_LA_DIAG_CHOL_ONELOOP_LBUB_6_6_6(plattMPC_FLOAT *H, plattMPC_FLOAT *llbysl, int* lbIdx, plattMPC_FLOAT *lubysu, int* ubIdx, plattMPC_FLOAT *Phi)


{
	int i;
	
	/* compute cholesky */
	for( i=0; i<6; i++ ){
		Phi[i] = H[i] + llbysl[i] + lubysu[i];

#if plattMPC_SET_PRINTLEVEL > 0 && defined PRINTNUMERICALWARNINGS
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
void plattMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_6_6(plattMPC_FLOAT *L, plattMPC_FLOAT *B, plattMPC_FLOAT *A)
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
void plattMPC_LA_DIAG_FORWARDSUB_6(plattMPC_FLOAT *L, plattMPC_FLOAT *b, plattMPC_FLOAT *y)
{
    int i;

    for( i=0; i<6; i++ ){
		y[i] = b[i]/L[i];
    }
}


/**
 * Compute L = B*B', where L is lower triangular of size NXp1
 * and B is a diagzero matrix of size [6 x 8] in column
 * storage format.
 * 
 */
void plattMPC_LA_DIAGZERO_MMT_6(plattMPC_FLOAT *B, plattMPC_FLOAT *L)
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
void plattMPC_LA_DIAGZERO_MVMSUB7_6(plattMPC_FLOAT *B, plattMPC_FLOAT *u, plattMPC_FLOAT *b, plattMPC_FLOAT *r)
{
	int i;

	for( i=0; i<6; i++ ){
		r[i] = b[i] - B[i]*u[i];
	}	
	
}


/**
 * Compute L = A*A' + B*B', where L is lower triangular of size NXp1
 * and A is a dense matrix of size [6 x 8] in column
 * storage format, and B is of size [6 x 8] diagonalzero
 * storage format.
 * 
 * THIS ONE HAS THE WORST ACCES PATTERN POSSIBLE. 
 * POSSIBKE FIX: PUT A AND B INTO ROW MAJOR FORMAT FIRST.
 * 
 */
void plattMPC_LA_DENSE_DIAGZERO_MMT2_6_8_8(plattMPC_FLOAT *A, plattMPC_FLOAT *B, plattMPC_FLOAT *L)
{
    int i, j, k, ii, di;
    plattMPC_FLOAT ltemp;
    
    ii = 0; di = 0;
    for( i=0; i<6; i++ ){        
        for( j=0; j<=i; j++ ){
            ltemp = 0; 
            for( k=0; k<8; k++ ){
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
void plattMPC_LA_DENSE_DIAGZERO_2MVMSUB2_6_8_8(plattMPC_FLOAT *A, plattMPC_FLOAT *x, plattMPC_FLOAT *B, plattMPC_FLOAT *u, plattMPC_FLOAT *b, plattMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;

	for( i=0; i<6; i++ ){
		r[i] = b[i] - A[k++]*x[0] - B[i]*u[i];
	}	

	for( j=1; j<8; j++ ){		
		for( i=0; i<6; i++ ){
			r[i] -= A[k++]*x[j];
		}
	}
	
}


/**
 * Compute L = A*A' + B*B', where L is lower triangular of size NXp1
 * and A is a dense matrix of size [6 x 8] in column
 * storage format, and B is of size [6 x 6] diagonalzero
 * storage format.
 * 
 * THIS ONE HAS THE WORST ACCES PATTERN POSSIBLE. 
 * POSSIBKE FIX: PUT A AND B INTO ROW MAJOR FORMAT FIRST.
 * 
 */
void plattMPC_LA_DENSE_DIAGZERO_MMT2_6_8_6(plattMPC_FLOAT *A, plattMPC_FLOAT *B, plattMPC_FLOAT *L)
{
    int i, j, k, ii, di;
    plattMPC_FLOAT ltemp;
    
    ii = 0; di = 0;
    for( i=0; i<6; i++ ){        
        for( j=0; j<=i; j++ ){
            ltemp = 0; 
            for( k=0; k<8; k++ ){
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
void plattMPC_LA_DENSE_DIAGZERO_2MVMSUB2_6_8_6(plattMPC_FLOAT *A, plattMPC_FLOAT *x, plattMPC_FLOAT *B, plattMPC_FLOAT *u, plattMPC_FLOAT *b, plattMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;

	for( i=0; i<6; i++ ){
		r[i] = b[i] - A[k++]*x[0] - B[i]*u[i];
	}	

	for( j=1; j<8; j++ ){		
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
void plattMPC_LA_DENSE_CHOL_6(plattMPC_FLOAT *A, plattMPC_FLOAT *L)
{
    int i, j, k, di, dj;
	 int ii, jj;

    plattMPC_FLOAT l;
    plattMPC_FLOAT Mii;

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
        
#if plattMPC_SET_PRINTLEVEL > 0 && defined PRINTNUMERICALWARNINGS
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
void plattMPC_LA_DENSE_FORWARDSUB_6(plattMPC_FLOAT *L, plattMPC_FLOAT *b, plattMPC_FLOAT *y)
{
    int i,j,ii,di;
    plattMPC_FLOAT yel;
            
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
void plattMPC_LA_DENSE_MATRIXTFORWARDSUB_6_6(plattMPC_FLOAT *L, plattMPC_FLOAT *B, plattMPC_FLOAT *A)
{
    int i,j,k,ii,di;
    plattMPC_FLOAT a;
    
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
void plattMPC_LA_DENSE_MMTSUB_6_6(plattMPC_FLOAT *A, plattMPC_FLOAT *L)
{
    int i, j, k, ii, di;
    plattMPC_FLOAT ltemp;
    
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
void plattMPC_LA_DENSE_MVMSUB1_6_6(plattMPC_FLOAT *A, plattMPC_FLOAT *x, plattMPC_FLOAT *b, plattMPC_FLOAT *r)
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
void plattMPC_LA_DENSE_BACKWARDSUB_6(plattMPC_FLOAT *L, plattMPC_FLOAT *y, plattMPC_FLOAT *x)
{
    int i, ii, di, j, jj, dj;
    plattMPC_FLOAT xel;    
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
void plattMPC_LA_DENSE_MTVMSUB_6_6(plattMPC_FLOAT *A, plattMPC_FLOAT *x, plattMPC_FLOAT *b, plattMPC_FLOAT *r)
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
 * Vector subtraction z = -x - y for vectors of length 118.
 */
void plattMPC_LA_VSUB2_118(plattMPC_FLOAT *x, plattMPC_FLOAT *y, plattMPC_FLOAT *z)
{
	int i;
	for( i=0; i<118; i++){
		z[i] = -x[i] - y[i];
	}
}


/**
 * Forward-Backward-Substitution to solve L*L^T*x = b where L is a
 * diagonal matrix of size 8 in vector
 * storage format.
 */
void plattMPC_LA_DIAG_FORWARDBACKWARDSUB_8(plattMPC_FLOAT *L, plattMPC_FLOAT *b, plattMPC_FLOAT *x)
{
    int i;
            
    /* solve Ly = b by forward and backward substitution */
    for( i=0; i<8; i++ ){
		x[i] = b[i]/(L[i]*L[i]);
    }
    
}


/**
 * Forward-Backward-Substitution to solve L*L^T*x = b where L is a
 * diagonal matrix of size 6 in vector
 * storage format.
 */
void plattMPC_LA_DIAG_FORWARDBACKWARDSUB_6(plattMPC_FLOAT *L, plattMPC_FLOAT *b, plattMPC_FLOAT *x)
{
    int i;
            
    /* solve Ly = b by forward and backward substitution */
    for( i=0; i<6; i++ ){
		x[i] = b[i]/(L[i]*L[i]);
    }
    
}


/*
 * Vector subtraction z = x(xidx) - y where y, z and xidx are of length 8,
 * and x has length 8 and is indexed through yidx.
 */
void plattMPC_LA_VSUB_INDEXED_8(plattMPC_FLOAT *x, int* xidx, plattMPC_FLOAT *y, plattMPC_FLOAT *z)
{
	int i;
	for( i=0; i<8; i++){
		z[i] = x[xidx[i]] - y[i];
	}
}


/*
 * Vector subtraction x = -u.*v - w for vectors of length 8.
 */
void plattMPC_LA_VSUB3_8(plattMPC_FLOAT *u, plattMPC_FLOAT *v, plattMPC_FLOAT *w, plattMPC_FLOAT *x)
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
void plattMPC_LA_VSUB2_INDEXED_8(plattMPC_FLOAT *x, plattMPC_FLOAT *y, int* yidx, plattMPC_FLOAT *z)
{
	int i;
	for( i=0; i<8; i++){
		z[i] = -x[i] - y[yidx[i]];
	}
}


/*
 * Vector subtraction z = x(xidx) - y where y, z and xidx are of length 6,
 * and x has length 6 and is indexed through yidx.
 */
void plattMPC_LA_VSUB_INDEXED_6(plattMPC_FLOAT *x, int* xidx, plattMPC_FLOAT *y, plattMPC_FLOAT *z)
{
	int i;
	for( i=0; i<6; i++){
		z[i] = x[xidx[i]] - y[i];
	}
}


/*
 * Vector subtraction x = -u.*v - w for vectors of length 6.
 */
void plattMPC_LA_VSUB3_6(plattMPC_FLOAT *u, plattMPC_FLOAT *v, plattMPC_FLOAT *w, plattMPC_FLOAT *x)
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
void plattMPC_LA_VSUB2_INDEXED_6(plattMPC_FLOAT *x, plattMPC_FLOAT *y, int* yidx, plattMPC_FLOAT *z)
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
 * plattMPC_NOPROGRESS (should be negative).
 */
int plattMPC_LINESEARCH_BACKTRACKING_AFFINE(plattMPC_FLOAT *l, plattMPC_FLOAT *s, plattMPC_FLOAT *dl, plattMPC_FLOAT *ds, plattMPC_FLOAT *a, plattMPC_FLOAT *mu_aff)
{
    int i;
	int lsIt=1;    
    plattMPC_FLOAT dltemp;
    plattMPC_FLOAT dstemp;
    plattMPC_FLOAT mya = 1.0;
    plattMPC_FLOAT mymu;
        
    while( 1 ){                        

        /* 
         * Compute both snew and wnew together.
         * We compute also mu_affine along the way here, as the
         * values might be in registers, so it should be cheaper.
         */
        mymu = 0;
        for( i=0; i<236; i++ ){
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
        if( i == 236 ){
            break;
        } else {
            mya *= plattMPC_SET_LS_SCALE_AFF;
            if( mya < plattMPC_SET_LS_MINSTEP ){
                return plattMPC_NOPROGRESS;
            }
        }
    }
    
    /* return new values and iteration counter */
    *a = mya;
    *mu_aff = mymu / (plattMPC_FLOAT)236;
    return lsIt;
}


/*
 * Vector subtraction x = u.*v - a where a is a scalar
*  and x,u,v are vectors of length 236.
 */
void plattMPC_LA_VSUB5_236(plattMPC_FLOAT *u, plattMPC_FLOAT *v, plattMPC_FLOAT a, plattMPC_FLOAT *x)
{
	int i;
	for( i=0; i<236; i++){
		x[i] = u[i]*v[i] - a;
	}
}


/*
 * Computes x=0; x(uidx) += u/su; x(vidx) -= v/sv where x is of length 8,
 * u, su, uidx are of length 8 and v, sv, vidx are of length 8.
 */
void plattMPC_LA_VSUB6_INDEXED_8_8_8(plattMPC_FLOAT *u, plattMPC_FLOAT *su, int* uidx, plattMPC_FLOAT *v, plattMPC_FLOAT *sv, int* vidx, plattMPC_FLOAT *x)
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
void plattMPC_LA_DIAGZERO_MVM_6(plattMPC_FLOAT *B, plattMPC_FLOAT *u, plattMPC_FLOAT *r)
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
void plattMPC_LA_DENSE_DIAGZERO_2MVMADD_6_8_8(plattMPC_FLOAT *A, plattMPC_FLOAT *x, plattMPC_FLOAT *B, plattMPC_FLOAT *u, plattMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;

	for( i=0; i<6; i++ ){
		r[i] = A[k++]*x[0] + B[i]*u[i];
	}	

	for( j=1; j<8; j++ ){		
		for( i=0; i<6; i++ ){
			r[i] += A[k++]*x[j];
		}
	}
	
}


/*
 * Computes x=0; x(uidx) += u/su; x(vidx) -= v/sv where x is of length 6,
 * u, su, uidx are of length 6 and v, sv, vidx are of length 6.
 */
void plattMPC_LA_VSUB6_INDEXED_6_6_6(plattMPC_FLOAT *u, plattMPC_FLOAT *su, int* uidx, plattMPC_FLOAT *v, plattMPC_FLOAT *sv, int* vidx, plattMPC_FLOAT *x)
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
void plattMPC_LA_DENSE_DIAGZERO_2MVMADD_6_8_6(plattMPC_FLOAT *A, plattMPC_FLOAT *x, plattMPC_FLOAT *B, plattMPC_FLOAT *u, plattMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;

	for( i=0; i<6; i++ ){
		r[i] = A[k++]*x[0] + B[i]*u[i];
	}	

	for( j=1; j<8; j++ ){		
		for( i=0; i<6; i++ ){
			r[i] += A[k++]*x[j];
		}
	}
	
}


/*
 * Vector subtraction z = x - y for vectors of length 118.
 */
void plattMPC_LA_VSUB_118(plattMPC_FLOAT *x, plattMPC_FLOAT *y, plattMPC_FLOAT *z)
{
	int i;
	for( i=0; i<118; i++){
		z[i] = x[i] - y[i];
	}
}


/** 
 * Computes z = -r./s - u.*y(y)
 * where all vectors except of y are of length 8 (length of y >= 8).
 */
void plattMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_8(plattMPC_FLOAT *r, plattMPC_FLOAT *s, plattMPC_FLOAT *u, plattMPC_FLOAT *y, int* yidx, plattMPC_FLOAT *z)
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
void plattMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_8(plattMPC_FLOAT *r, plattMPC_FLOAT *s, plattMPC_FLOAT *u, plattMPC_FLOAT *y, int* yidx, plattMPC_FLOAT *z)
{
	int i;
	for( i=0; i<8; i++ ){
		z[i] = -r[i]/s[i] + u[i]*y[yidx[i]];
	}
}


/** 
 * Computes z = -r./s - u.*y(y)
 * where all vectors except of y are of length 6 (length of y >= 6).
 */
void plattMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_6(plattMPC_FLOAT *r, plattMPC_FLOAT *s, plattMPC_FLOAT *u, plattMPC_FLOAT *y, int* yidx, plattMPC_FLOAT *z)
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
void plattMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_6(plattMPC_FLOAT *r, plattMPC_FLOAT *s, plattMPC_FLOAT *u, plattMPC_FLOAT *y, int* yidx, plattMPC_FLOAT *z)
{
	int i;
	for( i=0; i<6; i++ ){
		z[i] = -r[i]/s[i] + u[i]*y[yidx[i]];
	}
}


/*
 * Computes ds = -l.\(r + s.*dl) for vectors of length 236.
 */
void plattMPC_LA_VSUB7_236(plattMPC_FLOAT *l, plattMPC_FLOAT *r, plattMPC_FLOAT *s, plattMPC_FLOAT *dl, plattMPC_FLOAT *ds)
{
	int i;
	for( i=0; i<236; i++){
		ds[i] = -(r[i] + s[i]*dl[i])/l[i];
	}
}


/*
 * Vector addition x = x + y for vectors of length 118.
 */
void plattMPC_LA_VADD_118(plattMPC_FLOAT *x, plattMPC_FLOAT *y)
{
	int i;
	for( i=0; i<118; i++){
		x[i] += y[i];
	}
}


/*
 * Vector addition x = x + y for vectors of length 90.
 */
void plattMPC_LA_VADD_90(plattMPC_FLOAT *x, plattMPC_FLOAT *y)
{
	int i;
	for( i=0; i<90; i++){
		x[i] += y[i];
	}
}


/*
 * Vector addition x = x + y for vectors of length 236.
 */
void plattMPC_LA_VADD_236(plattMPC_FLOAT *x, plattMPC_FLOAT *y)
{
	int i;
	for( i=0; i<236; i++){
		x[i] += y[i];
	}
}


/**
 * Backtracking line search for combined predictor/corrector step.
 * Update on variables with safety factor gamma (to keep us away from
 * boundary).
 */
int plattMPC_LINESEARCH_BACKTRACKING_COMBINED(plattMPC_FLOAT *z, plattMPC_FLOAT *v, plattMPC_FLOAT *l, plattMPC_FLOAT *s, plattMPC_FLOAT *dz, plattMPC_FLOAT *dv, plattMPC_FLOAT *dl, plattMPC_FLOAT *ds, plattMPC_FLOAT *a, plattMPC_FLOAT *mu)
{
    int i, lsIt=1;       
    plattMPC_FLOAT dltemp;
    plattMPC_FLOAT dstemp;    
    plattMPC_FLOAT a_gamma;
            
    *a = 1.0;
    while( 1 ){                        

        /* check whether search criterion is fulfilled */
        for( i=0; i<236; i++ ){
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
        if( i == 236 ){
            break;
        } else {
            *a *= plattMPC_SET_LS_SCALE;
            if( *a < plattMPC_SET_LS_MINSTEP ){
                return plattMPC_NOPROGRESS;
            }
        }
    }
    
    /* update variables with safety margin */
    a_gamma = (*a)*plattMPC_SET_LS_MAXSTEP;
    
    /* primal variables */
    for( i=0; i<118; i++ ){
        z[i] += a_gamma*dz[i];
    }
    
    /* equality constraint multipliers */
    for( i=0; i<90; i++ ){
        v[i] += a_gamma*dv[i];
    }
    
    /* inequality constraint multipliers & slacks, also update mu */
    *mu = 0;
    for( i=0; i<236; i++ ){
        dltemp = l[i] + a_gamma*dl[i]; l[i] = dltemp;
        dstemp = s[i] + a_gamma*ds[i]; s[i] = dstemp;
        *mu += dltemp*dstemp;
    }
    
    *a = a_gamma;
    *mu /= (plattMPC_FLOAT)236;
    return lsIt;
}




/* VARIABLE DEFINITIONS ------------------------------------------------ */
plattMPC_FLOAT plattMPC_z[118];
plattMPC_FLOAT plattMPC_v[90];
plattMPC_FLOAT plattMPC_dz_aff[118];
plattMPC_FLOAT plattMPC_dv_aff[90];
plattMPC_FLOAT plattMPC_grad_cost[118];
plattMPC_FLOAT plattMPC_grad_eq[118];
plattMPC_FLOAT plattMPC_rd[118];
plattMPC_FLOAT plattMPC_l[236];
plattMPC_FLOAT plattMPC_s[236];
plattMPC_FLOAT plattMPC_lbys[236];
plattMPC_FLOAT plattMPC_dl_aff[236];
plattMPC_FLOAT plattMPC_ds_aff[236];
plattMPC_FLOAT plattMPC_dz_cc[118];
plattMPC_FLOAT plattMPC_dv_cc[90];
plattMPC_FLOAT plattMPC_dl_cc[236];
plattMPC_FLOAT plattMPC_ds_cc[236];
plattMPC_FLOAT plattMPC_ccrhs[236];
plattMPC_FLOAT plattMPC_grad_ineq[118];
plattMPC_FLOAT* plattMPC_z00 = plattMPC_z + 0;
plattMPC_FLOAT* plattMPC_dzaff00 = plattMPC_dz_aff + 0;
plattMPC_FLOAT* plattMPC_dzcc00 = plattMPC_dz_cc + 0;
plattMPC_FLOAT* plattMPC_rd00 = plattMPC_rd + 0;
plattMPC_FLOAT plattMPC_Lbyrd00[8];
plattMPC_FLOAT* plattMPC_grad_cost00 = plattMPC_grad_cost + 0;
plattMPC_FLOAT* plattMPC_grad_eq00 = plattMPC_grad_eq + 0;
plattMPC_FLOAT* plattMPC_grad_ineq00 = plattMPC_grad_ineq + 0;
plattMPC_FLOAT plattMPC_ctv00[8];
plattMPC_FLOAT plattMPC_C00[48] = {1.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 
0.0000000000000000E+000, 1.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 
0.0000000000000000E+000, 0.0000000000000000E+000, 1.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 
0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 1.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 
0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 1.0000000000000000E+000, 0.0000000000000000E+000, 
0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 1.0000000000000000E+000, 
1.0000000000000000E+000, 0.0000000000000000E+000, 1.0000000000000000E+000, 0.0000000000000000E+000, 1.0000000000000000E+000, 0.0000000000000000E+000, 
0.0000000000000000E+000, 1.0000000000000000E+000, 0.0000000000000000E+000, 1.0000000000000000E+000, 0.0000000000000000E+000, 1.0000000000000000E+000};
plattMPC_FLOAT* plattMPC_v00 = plattMPC_v + 0;
plattMPC_FLOAT plattMPC_re00[6];
plattMPC_FLOAT plattMPC_beta00[6];
plattMPC_FLOAT plattMPC_betacc00[6];
plattMPC_FLOAT* plattMPC_dvaff00 = plattMPC_dv_aff + 0;
plattMPC_FLOAT* plattMPC_dvcc00 = plattMPC_dv_cc + 0;
plattMPC_FLOAT plattMPC_V00[48];
plattMPC_FLOAT plattMPC_Yd00[21];
plattMPC_FLOAT plattMPC_Ld00[21];
plattMPC_FLOAT plattMPC_yy00[6];
plattMPC_FLOAT plattMPC_bmy00[6];
plattMPC_FLOAT plattMPC_c00[6] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
int plattMPC_lbIdx00[8] = {0, 1, 2, 3, 4, 5, 6, 7};
plattMPC_FLOAT* plattMPC_llb00 = plattMPC_l + 0;
plattMPC_FLOAT* plattMPC_slb00 = plattMPC_s + 0;
plattMPC_FLOAT* plattMPC_llbbyslb00 = plattMPC_lbys + 0;
plattMPC_FLOAT plattMPC_rilb00[8];
plattMPC_FLOAT* plattMPC_dllbaff00 = plattMPC_dl_aff + 0;
plattMPC_FLOAT* plattMPC_dslbaff00 = plattMPC_ds_aff + 0;
plattMPC_FLOAT* plattMPC_dllbcc00 = plattMPC_dl_cc + 0;
plattMPC_FLOAT* plattMPC_dslbcc00 = plattMPC_ds_cc + 0;
plattMPC_FLOAT* plattMPC_ccrhsl00 = plattMPC_ccrhs + 0;
int plattMPC_ubIdx00[8] = {0, 1, 2, 3, 4, 5, 6, 7};
plattMPC_FLOAT* plattMPC_lub00 = plattMPC_l + 8;
plattMPC_FLOAT* plattMPC_sub00 = plattMPC_s + 8;
plattMPC_FLOAT* plattMPC_lubbysub00 = plattMPC_lbys + 8;
plattMPC_FLOAT plattMPC_riub00[8];
plattMPC_FLOAT* plattMPC_dlubaff00 = plattMPC_dl_aff + 8;
plattMPC_FLOAT* plattMPC_dsubaff00 = plattMPC_ds_aff + 8;
plattMPC_FLOAT* plattMPC_dlubcc00 = plattMPC_dl_cc + 8;
plattMPC_FLOAT* plattMPC_dsubcc00 = plattMPC_ds_cc + 8;
plattMPC_FLOAT* plattMPC_ccrhsub00 = plattMPC_ccrhs + 8;
plattMPC_FLOAT plattMPC_Phi00[8];
plattMPC_FLOAT plattMPC_D00[8] = {1.0000000000000000E+000, 
1.0000000000000000E+000, 
1.0000000000000000E+000, 
1.0000000000000000E+000, 
1.0000000000000000E+000, 
1.0000000000000000E+000};
plattMPC_FLOAT plattMPC_W00[8];
plattMPC_FLOAT* plattMPC_z01 = plattMPC_z + 8;
plattMPC_FLOAT* plattMPC_dzaff01 = plattMPC_dz_aff + 8;
plattMPC_FLOAT* plattMPC_dzcc01 = plattMPC_dz_cc + 8;
plattMPC_FLOAT* plattMPC_rd01 = plattMPC_rd + 8;
plattMPC_FLOAT plattMPC_Lbyrd01[8];
plattMPC_FLOAT* plattMPC_grad_cost01 = plattMPC_grad_cost + 8;
plattMPC_FLOAT* plattMPC_grad_eq01 = plattMPC_grad_eq + 8;
plattMPC_FLOAT* plattMPC_grad_ineq01 = plattMPC_grad_ineq + 8;
plattMPC_FLOAT plattMPC_ctv01[8];
plattMPC_FLOAT* plattMPC_v01 = plattMPC_v + 6;
plattMPC_FLOAT plattMPC_re01[6];
plattMPC_FLOAT plattMPC_beta01[6];
plattMPC_FLOAT plattMPC_betacc01[6];
plattMPC_FLOAT* plattMPC_dvaff01 = plattMPC_dv_aff + 6;
plattMPC_FLOAT* plattMPC_dvcc01 = plattMPC_dv_cc + 6;
plattMPC_FLOAT plattMPC_V01[48];
plattMPC_FLOAT plattMPC_Yd01[21];
plattMPC_FLOAT plattMPC_Ld01[21];
plattMPC_FLOAT plattMPC_yy01[6];
plattMPC_FLOAT plattMPC_bmy01[6];
plattMPC_FLOAT plattMPC_c01[6] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
int plattMPC_lbIdx01[8] = {0, 1, 2, 3, 4, 5, 6, 7};
plattMPC_FLOAT* plattMPC_llb01 = plattMPC_l + 16;
plattMPC_FLOAT* plattMPC_slb01 = plattMPC_s + 16;
plattMPC_FLOAT* plattMPC_llbbyslb01 = plattMPC_lbys + 16;
plattMPC_FLOAT plattMPC_rilb01[8];
plattMPC_FLOAT* plattMPC_dllbaff01 = plattMPC_dl_aff + 16;
plattMPC_FLOAT* plattMPC_dslbaff01 = plattMPC_ds_aff + 16;
plattMPC_FLOAT* plattMPC_dllbcc01 = plattMPC_dl_cc + 16;
plattMPC_FLOAT* plattMPC_dslbcc01 = plattMPC_ds_cc + 16;
plattMPC_FLOAT* plattMPC_ccrhsl01 = plattMPC_ccrhs + 16;
int plattMPC_ubIdx01[8] = {0, 1, 2, 3, 4, 5, 6, 7};
plattMPC_FLOAT* plattMPC_lub01 = plattMPC_l + 24;
plattMPC_FLOAT* plattMPC_sub01 = plattMPC_s + 24;
plattMPC_FLOAT* plattMPC_lubbysub01 = plattMPC_lbys + 24;
plattMPC_FLOAT plattMPC_riub01[8];
plattMPC_FLOAT* plattMPC_dlubaff01 = plattMPC_dl_aff + 24;
plattMPC_FLOAT* plattMPC_dsubaff01 = plattMPC_ds_aff + 24;
plattMPC_FLOAT* plattMPC_dlubcc01 = plattMPC_dl_cc + 24;
plattMPC_FLOAT* plattMPC_dsubcc01 = plattMPC_ds_cc + 24;
plattMPC_FLOAT* plattMPC_ccrhsub01 = plattMPC_ccrhs + 24;
plattMPC_FLOAT plattMPC_Phi01[8];
plattMPC_FLOAT plattMPC_D01[8] = {-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000};
plattMPC_FLOAT plattMPC_W01[8];
plattMPC_FLOAT plattMPC_Ysd01[36];
plattMPC_FLOAT plattMPC_Lsd01[36];
plattMPC_FLOAT* plattMPC_z02 = plattMPC_z + 16;
plattMPC_FLOAT* plattMPC_dzaff02 = plattMPC_dz_aff + 16;
plattMPC_FLOAT* plattMPC_dzcc02 = plattMPC_dz_cc + 16;
plattMPC_FLOAT* plattMPC_rd02 = plattMPC_rd + 16;
plattMPC_FLOAT plattMPC_Lbyrd02[8];
plattMPC_FLOAT* plattMPC_grad_cost02 = plattMPC_grad_cost + 16;
plattMPC_FLOAT* plattMPC_grad_eq02 = plattMPC_grad_eq + 16;
plattMPC_FLOAT* plattMPC_grad_ineq02 = plattMPC_grad_ineq + 16;
plattMPC_FLOAT plattMPC_ctv02[8];
plattMPC_FLOAT* plattMPC_v02 = plattMPC_v + 12;
plattMPC_FLOAT plattMPC_re02[6];
plattMPC_FLOAT plattMPC_beta02[6];
plattMPC_FLOAT plattMPC_betacc02[6];
plattMPC_FLOAT* plattMPC_dvaff02 = plattMPC_dv_aff + 12;
plattMPC_FLOAT* plattMPC_dvcc02 = plattMPC_dv_cc + 12;
plattMPC_FLOAT plattMPC_V02[48];
plattMPC_FLOAT plattMPC_Yd02[21];
plattMPC_FLOAT plattMPC_Ld02[21];
plattMPC_FLOAT plattMPC_yy02[6];
plattMPC_FLOAT plattMPC_bmy02[6];
plattMPC_FLOAT plattMPC_c02[6] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
int plattMPC_lbIdx02[8] = {0, 1, 2, 3, 4, 5, 6, 7};
plattMPC_FLOAT* plattMPC_llb02 = plattMPC_l + 32;
plattMPC_FLOAT* plattMPC_slb02 = plattMPC_s + 32;
plattMPC_FLOAT* plattMPC_llbbyslb02 = plattMPC_lbys + 32;
plattMPC_FLOAT plattMPC_rilb02[8];
plattMPC_FLOAT* plattMPC_dllbaff02 = plattMPC_dl_aff + 32;
plattMPC_FLOAT* plattMPC_dslbaff02 = plattMPC_ds_aff + 32;
plattMPC_FLOAT* plattMPC_dllbcc02 = plattMPC_dl_cc + 32;
plattMPC_FLOAT* plattMPC_dslbcc02 = plattMPC_ds_cc + 32;
plattMPC_FLOAT* plattMPC_ccrhsl02 = plattMPC_ccrhs + 32;
int plattMPC_ubIdx02[8] = {0, 1, 2, 3, 4, 5, 6, 7};
plattMPC_FLOAT* plattMPC_lub02 = plattMPC_l + 40;
plattMPC_FLOAT* plattMPC_sub02 = plattMPC_s + 40;
plattMPC_FLOAT* plattMPC_lubbysub02 = plattMPC_lbys + 40;
plattMPC_FLOAT plattMPC_riub02[8];
plattMPC_FLOAT* plattMPC_dlubaff02 = plattMPC_dl_aff + 40;
plattMPC_FLOAT* plattMPC_dsubaff02 = plattMPC_ds_aff + 40;
plattMPC_FLOAT* plattMPC_dlubcc02 = plattMPC_dl_cc + 40;
plattMPC_FLOAT* plattMPC_dsubcc02 = plattMPC_ds_cc + 40;
plattMPC_FLOAT* plattMPC_ccrhsub02 = plattMPC_ccrhs + 40;
plattMPC_FLOAT plattMPC_Phi02[8];
plattMPC_FLOAT plattMPC_W02[8];
plattMPC_FLOAT plattMPC_Ysd02[36];
plattMPC_FLOAT plattMPC_Lsd02[36];
plattMPC_FLOAT* plattMPC_z03 = plattMPC_z + 24;
plattMPC_FLOAT* plattMPC_dzaff03 = plattMPC_dz_aff + 24;
plattMPC_FLOAT* plattMPC_dzcc03 = plattMPC_dz_cc + 24;
plattMPC_FLOAT* plattMPC_rd03 = plattMPC_rd + 24;
plattMPC_FLOAT plattMPC_Lbyrd03[8];
plattMPC_FLOAT* plattMPC_grad_cost03 = plattMPC_grad_cost + 24;
plattMPC_FLOAT* plattMPC_grad_eq03 = plattMPC_grad_eq + 24;
plattMPC_FLOAT* plattMPC_grad_ineq03 = plattMPC_grad_ineq + 24;
plattMPC_FLOAT plattMPC_ctv03[8];
plattMPC_FLOAT* plattMPC_v03 = plattMPC_v + 18;
plattMPC_FLOAT plattMPC_re03[6];
plattMPC_FLOAT plattMPC_beta03[6];
plattMPC_FLOAT plattMPC_betacc03[6];
plattMPC_FLOAT* plattMPC_dvaff03 = plattMPC_dv_aff + 18;
plattMPC_FLOAT* plattMPC_dvcc03 = plattMPC_dv_cc + 18;
plattMPC_FLOAT plattMPC_V03[48];
plattMPC_FLOAT plattMPC_Yd03[21];
plattMPC_FLOAT plattMPC_Ld03[21];
plattMPC_FLOAT plattMPC_yy03[6];
plattMPC_FLOAT plattMPC_bmy03[6];
plattMPC_FLOAT plattMPC_c03[6] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
int plattMPC_lbIdx03[8] = {0, 1, 2, 3, 4, 5, 6, 7};
plattMPC_FLOAT* plattMPC_llb03 = plattMPC_l + 48;
plattMPC_FLOAT* plattMPC_slb03 = plattMPC_s + 48;
plattMPC_FLOAT* plattMPC_llbbyslb03 = plattMPC_lbys + 48;
plattMPC_FLOAT plattMPC_rilb03[8];
plattMPC_FLOAT* plattMPC_dllbaff03 = plattMPC_dl_aff + 48;
plattMPC_FLOAT* plattMPC_dslbaff03 = plattMPC_ds_aff + 48;
plattMPC_FLOAT* plattMPC_dllbcc03 = plattMPC_dl_cc + 48;
plattMPC_FLOAT* plattMPC_dslbcc03 = plattMPC_ds_cc + 48;
plattMPC_FLOAT* plattMPC_ccrhsl03 = plattMPC_ccrhs + 48;
int plattMPC_ubIdx03[8] = {0, 1, 2, 3, 4, 5, 6, 7};
plattMPC_FLOAT* plattMPC_lub03 = plattMPC_l + 56;
plattMPC_FLOAT* plattMPC_sub03 = plattMPC_s + 56;
plattMPC_FLOAT* plattMPC_lubbysub03 = plattMPC_lbys + 56;
plattMPC_FLOAT plattMPC_riub03[8];
plattMPC_FLOAT* plattMPC_dlubaff03 = plattMPC_dl_aff + 56;
plattMPC_FLOAT* plattMPC_dsubaff03 = plattMPC_ds_aff + 56;
plattMPC_FLOAT* plattMPC_dlubcc03 = plattMPC_dl_cc + 56;
plattMPC_FLOAT* plattMPC_dsubcc03 = plattMPC_ds_cc + 56;
plattMPC_FLOAT* plattMPC_ccrhsub03 = plattMPC_ccrhs + 56;
plattMPC_FLOAT plattMPC_Phi03[8];
plattMPC_FLOAT plattMPC_W03[8];
plattMPC_FLOAT plattMPC_Ysd03[36];
plattMPC_FLOAT plattMPC_Lsd03[36];
plattMPC_FLOAT* plattMPC_z04 = plattMPC_z + 32;
plattMPC_FLOAT* plattMPC_dzaff04 = plattMPC_dz_aff + 32;
plattMPC_FLOAT* plattMPC_dzcc04 = plattMPC_dz_cc + 32;
plattMPC_FLOAT* plattMPC_rd04 = plattMPC_rd + 32;
plattMPC_FLOAT plattMPC_Lbyrd04[8];
plattMPC_FLOAT* plattMPC_grad_cost04 = plattMPC_grad_cost + 32;
plattMPC_FLOAT* plattMPC_grad_eq04 = plattMPC_grad_eq + 32;
plattMPC_FLOAT* plattMPC_grad_ineq04 = plattMPC_grad_ineq + 32;
plattMPC_FLOAT plattMPC_ctv04[8];
plattMPC_FLOAT* plattMPC_v04 = plattMPC_v + 24;
plattMPC_FLOAT plattMPC_re04[6];
plattMPC_FLOAT plattMPC_beta04[6];
plattMPC_FLOAT plattMPC_betacc04[6];
plattMPC_FLOAT* plattMPC_dvaff04 = plattMPC_dv_aff + 24;
plattMPC_FLOAT* plattMPC_dvcc04 = plattMPC_dv_cc + 24;
plattMPC_FLOAT plattMPC_V04[48];
plattMPC_FLOAT plattMPC_Yd04[21];
plattMPC_FLOAT plattMPC_Ld04[21];
plattMPC_FLOAT plattMPC_yy04[6];
plattMPC_FLOAT plattMPC_bmy04[6];
plattMPC_FLOAT plattMPC_c04[6] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
int plattMPC_lbIdx04[8] = {0, 1, 2, 3, 4, 5, 6, 7};
plattMPC_FLOAT* plattMPC_llb04 = plattMPC_l + 64;
plattMPC_FLOAT* plattMPC_slb04 = plattMPC_s + 64;
plattMPC_FLOAT* plattMPC_llbbyslb04 = plattMPC_lbys + 64;
plattMPC_FLOAT plattMPC_rilb04[8];
plattMPC_FLOAT* plattMPC_dllbaff04 = plattMPC_dl_aff + 64;
plattMPC_FLOAT* plattMPC_dslbaff04 = plattMPC_ds_aff + 64;
plattMPC_FLOAT* plattMPC_dllbcc04 = plattMPC_dl_cc + 64;
plattMPC_FLOAT* plattMPC_dslbcc04 = plattMPC_ds_cc + 64;
plattMPC_FLOAT* plattMPC_ccrhsl04 = plattMPC_ccrhs + 64;
int plattMPC_ubIdx04[8] = {0, 1, 2, 3, 4, 5, 6, 7};
plattMPC_FLOAT* plattMPC_lub04 = plattMPC_l + 72;
plattMPC_FLOAT* plattMPC_sub04 = plattMPC_s + 72;
plattMPC_FLOAT* plattMPC_lubbysub04 = plattMPC_lbys + 72;
plattMPC_FLOAT plattMPC_riub04[8];
plattMPC_FLOAT* plattMPC_dlubaff04 = plattMPC_dl_aff + 72;
plattMPC_FLOAT* plattMPC_dsubaff04 = plattMPC_ds_aff + 72;
plattMPC_FLOAT* plattMPC_dlubcc04 = plattMPC_dl_cc + 72;
plattMPC_FLOAT* plattMPC_dsubcc04 = plattMPC_ds_cc + 72;
plattMPC_FLOAT* plattMPC_ccrhsub04 = plattMPC_ccrhs + 72;
plattMPC_FLOAT plattMPC_Phi04[8];
plattMPC_FLOAT plattMPC_W04[8];
plattMPC_FLOAT plattMPC_Ysd04[36];
plattMPC_FLOAT plattMPC_Lsd04[36];
plattMPC_FLOAT* plattMPC_z05 = plattMPC_z + 40;
plattMPC_FLOAT* plattMPC_dzaff05 = plattMPC_dz_aff + 40;
plattMPC_FLOAT* plattMPC_dzcc05 = plattMPC_dz_cc + 40;
plattMPC_FLOAT* plattMPC_rd05 = plattMPC_rd + 40;
plattMPC_FLOAT plattMPC_Lbyrd05[8];
plattMPC_FLOAT* plattMPC_grad_cost05 = plattMPC_grad_cost + 40;
plattMPC_FLOAT* plattMPC_grad_eq05 = plattMPC_grad_eq + 40;
plattMPC_FLOAT* plattMPC_grad_ineq05 = plattMPC_grad_ineq + 40;
plattMPC_FLOAT plattMPC_ctv05[8];
plattMPC_FLOAT* plattMPC_v05 = plattMPC_v + 30;
plattMPC_FLOAT plattMPC_re05[6];
plattMPC_FLOAT plattMPC_beta05[6];
plattMPC_FLOAT plattMPC_betacc05[6];
plattMPC_FLOAT* plattMPC_dvaff05 = plattMPC_dv_aff + 30;
plattMPC_FLOAT* plattMPC_dvcc05 = plattMPC_dv_cc + 30;
plattMPC_FLOAT plattMPC_V05[48];
plattMPC_FLOAT plattMPC_Yd05[21];
plattMPC_FLOAT plattMPC_Ld05[21];
plattMPC_FLOAT plattMPC_yy05[6];
plattMPC_FLOAT plattMPC_bmy05[6];
plattMPC_FLOAT plattMPC_c05[6] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
int plattMPC_lbIdx05[8] = {0, 1, 2, 3, 4, 5, 6, 7};
plattMPC_FLOAT* plattMPC_llb05 = plattMPC_l + 80;
plattMPC_FLOAT* plattMPC_slb05 = plattMPC_s + 80;
plattMPC_FLOAT* plattMPC_llbbyslb05 = plattMPC_lbys + 80;
plattMPC_FLOAT plattMPC_rilb05[8];
plattMPC_FLOAT* plattMPC_dllbaff05 = plattMPC_dl_aff + 80;
plattMPC_FLOAT* plattMPC_dslbaff05 = plattMPC_ds_aff + 80;
plattMPC_FLOAT* plattMPC_dllbcc05 = plattMPC_dl_cc + 80;
plattMPC_FLOAT* plattMPC_dslbcc05 = plattMPC_ds_cc + 80;
plattMPC_FLOAT* plattMPC_ccrhsl05 = plattMPC_ccrhs + 80;
int plattMPC_ubIdx05[8] = {0, 1, 2, 3, 4, 5, 6, 7};
plattMPC_FLOAT* plattMPC_lub05 = plattMPC_l + 88;
plattMPC_FLOAT* plattMPC_sub05 = plattMPC_s + 88;
plattMPC_FLOAT* plattMPC_lubbysub05 = plattMPC_lbys + 88;
plattMPC_FLOAT plattMPC_riub05[8];
plattMPC_FLOAT* plattMPC_dlubaff05 = plattMPC_dl_aff + 88;
plattMPC_FLOAT* plattMPC_dsubaff05 = plattMPC_ds_aff + 88;
plattMPC_FLOAT* plattMPC_dlubcc05 = plattMPC_dl_cc + 88;
plattMPC_FLOAT* plattMPC_dsubcc05 = plattMPC_ds_cc + 88;
plattMPC_FLOAT* plattMPC_ccrhsub05 = plattMPC_ccrhs + 88;
plattMPC_FLOAT plattMPC_Phi05[8];
plattMPC_FLOAT plattMPC_W05[8];
plattMPC_FLOAT plattMPC_Ysd05[36];
plattMPC_FLOAT plattMPC_Lsd05[36];
plattMPC_FLOAT* plattMPC_z06 = plattMPC_z + 48;
plattMPC_FLOAT* plattMPC_dzaff06 = plattMPC_dz_aff + 48;
plattMPC_FLOAT* plattMPC_dzcc06 = plattMPC_dz_cc + 48;
plattMPC_FLOAT* plattMPC_rd06 = plattMPC_rd + 48;
plattMPC_FLOAT plattMPC_Lbyrd06[8];
plattMPC_FLOAT* plattMPC_grad_cost06 = plattMPC_grad_cost + 48;
plattMPC_FLOAT* plattMPC_grad_eq06 = plattMPC_grad_eq + 48;
plattMPC_FLOAT* plattMPC_grad_ineq06 = plattMPC_grad_ineq + 48;
plattMPC_FLOAT plattMPC_ctv06[8];
plattMPC_FLOAT* plattMPC_v06 = plattMPC_v + 36;
plattMPC_FLOAT plattMPC_re06[6];
plattMPC_FLOAT plattMPC_beta06[6];
plattMPC_FLOAT plattMPC_betacc06[6];
plattMPC_FLOAT* plattMPC_dvaff06 = plattMPC_dv_aff + 36;
plattMPC_FLOAT* plattMPC_dvcc06 = plattMPC_dv_cc + 36;
plattMPC_FLOAT plattMPC_V06[48];
plattMPC_FLOAT plattMPC_Yd06[21];
plattMPC_FLOAT plattMPC_Ld06[21];
plattMPC_FLOAT plattMPC_yy06[6];
plattMPC_FLOAT plattMPC_bmy06[6];
plattMPC_FLOAT plattMPC_c06[6] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
int plattMPC_lbIdx06[8] = {0, 1, 2, 3, 4, 5, 6, 7};
plattMPC_FLOAT* plattMPC_llb06 = plattMPC_l + 96;
plattMPC_FLOAT* plattMPC_slb06 = plattMPC_s + 96;
plattMPC_FLOAT* plattMPC_llbbyslb06 = plattMPC_lbys + 96;
plattMPC_FLOAT plattMPC_rilb06[8];
plattMPC_FLOAT* plattMPC_dllbaff06 = plattMPC_dl_aff + 96;
plattMPC_FLOAT* plattMPC_dslbaff06 = plattMPC_ds_aff + 96;
plattMPC_FLOAT* plattMPC_dllbcc06 = plattMPC_dl_cc + 96;
plattMPC_FLOAT* plattMPC_dslbcc06 = plattMPC_ds_cc + 96;
plattMPC_FLOAT* plattMPC_ccrhsl06 = plattMPC_ccrhs + 96;
int plattMPC_ubIdx06[8] = {0, 1, 2, 3, 4, 5, 6, 7};
plattMPC_FLOAT* plattMPC_lub06 = plattMPC_l + 104;
plattMPC_FLOAT* plattMPC_sub06 = plattMPC_s + 104;
plattMPC_FLOAT* plattMPC_lubbysub06 = plattMPC_lbys + 104;
plattMPC_FLOAT plattMPC_riub06[8];
plattMPC_FLOAT* plattMPC_dlubaff06 = plattMPC_dl_aff + 104;
plattMPC_FLOAT* plattMPC_dsubaff06 = plattMPC_ds_aff + 104;
plattMPC_FLOAT* plattMPC_dlubcc06 = plattMPC_dl_cc + 104;
plattMPC_FLOAT* plattMPC_dsubcc06 = plattMPC_ds_cc + 104;
plattMPC_FLOAT* plattMPC_ccrhsub06 = plattMPC_ccrhs + 104;
plattMPC_FLOAT plattMPC_Phi06[8];
plattMPC_FLOAT plattMPC_W06[8];
plattMPC_FLOAT plattMPC_Ysd06[36];
plattMPC_FLOAT plattMPC_Lsd06[36];
plattMPC_FLOAT* plattMPC_z07 = plattMPC_z + 56;
plattMPC_FLOAT* plattMPC_dzaff07 = plattMPC_dz_aff + 56;
plattMPC_FLOAT* plattMPC_dzcc07 = plattMPC_dz_cc + 56;
plattMPC_FLOAT* plattMPC_rd07 = plattMPC_rd + 56;
plattMPC_FLOAT plattMPC_Lbyrd07[8];
plattMPC_FLOAT* plattMPC_grad_cost07 = plattMPC_grad_cost + 56;
plattMPC_FLOAT* plattMPC_grad_eq07 = plattMPC_grad_eq + 56;
plattMPC_FLOAT* plattMPC_grad_ineq07 = plattMPC_grad_ineq + 56;
plattMPC_FLOAT plattMPC_ctv07[8];
plattMPC_FLOAT* plattMPC_v07 = plattMPC_v + 42;
plattMPC_FLOAT plattMPC_re07[6];
plattMPC_FLOAT plattMPC_beta07[6];
plattMPC_FLOAT plattMPC_betacc07[6];
plattMPC_FLOAT* plattMPC_dvaff07 = plattMPC_dv_aff + 42;
plattMPC_FLOAT* plattMPC_dvcc07 = plattMPC_dv_cc + 42;
plattMPC_FLOAT plattMPC_V07[48];
plattMPC_FLOAT plattMPC_Yd07[21];
plattMPC_FLOAT plattMPC_Ld07[21];
plattMPC_FLOAT plattMPC_yy07[6];
plattMPC_FLOAT plattMPC_bmy07[6];
plattMPC_FLOAT plattMPC_c07[6] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
int plattMPC_lbIdx07[8] = {0, 1, 2, 3, 4, 5, 6, 7};
plattMPC_FLOAT* plattMPC_llb07 = plattMPC_l + 112;
plattMPC_FLOAT* plattMPC_slb07 = plattMPC_s + 112;
plattMPC_FLOAT* plattMPC_llbbyslb07 = plattMPC_lbys + 112;
plattMPC_FLOAT plattMPC_rilb07[8];
plattMPC_FLOAT* plattMPC_dllbaff07 = plattMPC_dl_aff + 112;
plattMPC_FLOAT* plattMPC_dslbaff07 = plattMPC_ds_aff + 112;
plattMPC_FLOAT* plattMPC_dllbcc07 = plattMPC_dl_cc + 112;
plattMPC_FLOAT* plattMPC_dslbcc07 = plattMPC_ds_cc + 112;
plattMPC_FLOAT* plattMPC_ccrhsl07 = plattMPC_ccrhs + 112;
int plattMPC_ubIdx07[8] = {0, 1, 2, 3, 4, 5, 6, 7};
plattMPC_FLOAT* plattMPC_lub07 = plattMPC_l + 120;
plattMPC_FLOAT* plattMPC_sub07 = plattMPC_s + 120;
plattMPC_FLOAT* plattMPC_lubbysub07 = plattMPC_lbys + 120;
plattMPC_FLOAT plattMPC_riub07[8];
plattMPC_FLOAT* plattMPC_dlubaff07 = plattMPC_dl_aff + 120;
plattMPC_FLOAT* plattMPC_dsubaff07 = plattMPC_ds_aff + 120;
plattMPC_FLOAT* plattMPC_dlubcc07 = plattMPC_dl_cc + 120;
plattMPC_FLOAT* plattMPC_dsubcc07 = plattMPC_ds_cc + 120;
plattMPC_FLOAT* plattMPC_ccrhsub07 = plattMPC_ccrhs + 120;
plattMPC_FLOAT plattMPC_Phi07[8];
plattMPC_FLOAT plattMPC_W07[8];
plattMPC_FLOAT plattMPC_Ysd07[36];
plattMPC_FLOAT plattMPC_Lsd07[36];
plattMPC_FLOAT* plattMPC_z08 = plattMPC_z + 64;
plattMPC_FLOAT* plattMPC_dzaff08 = plattMPC_dz_aff + 64;
plattMPC_FLOAT* plattMPC_dzcc08 = plattMPC_dz_cc + 64;
plattMPC_FLOAT* plattMPC_rd08 = plattMPC_rd + 64;
plattMPC_FLOAT plattMPC_Lbyrd08[8];
plattMPC_FLOAT* plattMPC_grad_cost08 = plattMPC_grad_cost + 64;
plattMPC_FLOAT* plattMPC_grad_eq08 = plattMPC_grad_eq + 64;
plattMPC_FLOAT* plattMPC_grad_ineq08 = plattMPC_grad_ineq + 64;
plattMPC_FLOAT plattMPC_ctv08[8];
plattMPC_FLOAT* plattMPC_v08 = plattMPC_v + 48;
plattMPC_FLOAT plattMPC_re08[6];
plattMPC_FLOAT plattMPC_beta08[6];
plattMPC_FLOAT plattMPC_betacc08[6];
plattMPC_FLOAT* plattMPC_dvaff08 = plattMPC_dv_aff + 48;
plattMPC_FLOAT* plattMPC_dvcc08 = plattMPC_dv_cc + 48;
plattMPC_FLOAT plattMPC_V08[48];
plattMPC_FLOAT plattMPC_Yd08[21];
plattMPC_FLOAT plattMPC_Ld08[21];
plattMPC_FLOAT plattMPC_yy08[6];
plattMPC_FLOAT plattMPC_bmy08[6];
plattMPC_FLOAT plattMPC_c08[6] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
int plattMPC_lbIdx08[8] = {0, 1, 2, 3, 4, 5, 6, 7};
plattMPC_FLOAT* plattMPC_llb08 = plattMPC_l + 128;
plattMPC_FLOAT* plattMPC_slb08 = plattMPC_s + 128;
plattMPC_FLOAT* plattMPC_llbbyslb08 = plattMPC_lbys + 128;
plattMPC_FLOAT plattMPC_rilb08[8];
plattMPC_FLOAT* plattMPC_dllbaff08 = plattMPC_dl_aff + 128;
plattMPC_FLOAT* plattMPC_dslbaff08 = plattMPC_ds_aff + 128;
plattMPC_FLOAT* plattMPC_dllbcc08 = plattMPC_dl_cc + 128;
plattMPC_FLOAT* plattMPC_dslbcc08 = plattMPC_ds_cc + 128;
plattMPC_FLOAT* plattMPC_ccrhsl08 = plattMPC_ccrhs + 128;
int plattMPC_ubIdx08[8] = {0, 1, 2, 3, 4, 5, 6, 7};
plattMPC_FLOAT* plattMPC_lub08 = plattMPC_l + 136;
plattMPC_FLOAT* plattMPC_sub08 = plattMPC_s + 136;
plattMPC_FLOAT* plattMPC_lubbysub08 = plattMPC_lbys + 136;
plattMPC_FLOAT plattMPC_riub08[8];
plattMPC_FLOAT* plattMPC_dlubaff08 = plattMPC_dl_aff + 136;
plattMPC_FLOAT* plattMPC_dsubaff08 = plattMPC_ds_aff + 136;
plattMPC_FLOAT* plattMPC_dlubcc08 = plattMPC_dl_cc + 136;
plattMPC_FLOAT* plattMPC_dsubcc08 = plattMPC_ds_cc + 136;
plattMPC_FLOAT* plattMPC_ccrhsub08 = plattMPC_ccrhs + 136;
plattMPC_FLOAT plattMPC_Phi08[8];
plattMPC_FLOAT plattMPC_W08[8];
plattMPC_FLOAT plattMPC_Ysd08[36];
plattMPC_FLOAT plattMPC_Lsd08[36];
plattMPC_FLOAT* plattMPC_z09 = plattMPC_z + 72;
plattMPC_FLOAT* plattMPC_dzaff09 = plattMPC_dz_aff + 72;
plattMPC_FLOAT* plattMPC_dzcc09 = plattMPC_dz_cc + 72;
plattMPC_FLOAT* plattMPC_rd09 = plattMPC_rd + 72;
plattMPC_FLOAT plattMPC_Lbyrd09[8];
plattMPC_FLOAT* plattMPC_grad_cost09 = plattMPC_grad_cost + 72;
plattMPC_FLOAT* plattMPC_grad_eq09 = plattMPC_grad_eq + 72;
plattMPC_FLOAT* plattMPC_grad_ineq09 = plattMPC_grad_ineq + 72;
plattMPC_FLOAT plattMPC_ctv09[8];
plattMPC_FLOAT* plattMPC_v09 = plattMPC_v + 54;
plattMPC_FLOAT plattMPC_re09[6];
plattMPC_FLOAT plattMPC_beta09[6];
plattMPC_FLOAT plattMPC_betacc09[6];
plattMPC_FLOAT* plattMPC_dvaff09 = plattMPC_dv_aff + 54;
plattMPC_FLOAT* plattMPC_dvcc09 = plattMPC_dv_cc + 54;
plattMPC_FLOAT plattMPC_V09[48];
plattMPC_FLOAT plattMPC_Yd09[21];
plattMPC_FLOAT plattMPC_Ld09[21];
plattMPC_FLOAT plattMPC_yy09[6];
plattMPC_FLOAT plattMPC_bmy09[6];
plattMPC_FLOAT plattMPC_c09[6] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
int plattMPC_lbIdx09[8] = {0, 1, 2, 3, 4, 5, 6, 7};
plattMPC_FLOAT* plattMPC_llb09 = plattMPC_l + 144;
plattMPC_FLOAT* plattMPC_slb09 = plattMPC_s + 144;
plattMPC_FLOAT* plattMPC_llbbyslb09 = plattMPC_lbys + 144;
plattMPC_FLOAT plattMPC_rilb09[8];
plattMPC_FLOAT* plattMPC_dllbaff09 = plattMPC_dl_aff + 144;
plattMPC_FLOAT* plattMPC_dslbaff09 = plattMPC_ds_aff + 144;
plattMPC_FLOAT* plattMPC_dllbcc09 = plattMPC_dl_cc + 144;
plattMPC_FLOAT* plattMPC_dslbcc09 = plattMPC_ds_cc + 144;
plattMPC_FLOAT* plattMPC_ccrhsl09 = plattMPC_ccrhs + 144;
int plattMPC_ubIdx09[8] = {0, 1, 2, 3, 4, 5, 6, 7};
plattMPC_FLOAT* plattMPC_lub09 = plattMPC_l + 152;
plattMPC_FLOAT* plattMPC_sub09 = plattMPC_s + 152;
plattMPC_FLOAT* plattMPC_lubbysub09 = plattMPC_lbys + 152;
plattMPC_FLOAT plattMPC_riub09[8];
plattMPC_FLOAT* plattMPC_dlubaff09 = plattMPC_dl_aff + 152;
plattMPC_FLOAT* plattMPC_dsubaff09 = plattMPC_ds_aff + 152;
plattMPC_FLOAT* plattMPC_dlubcc09 = plattMPC_dl_cc + 152;
plattMPC_FLOAT* plattMPC_dsubcc09 = plattMPC_ds_cc + 152;
plattMPC_FLOAT* plattMPC_ccrhsub09 = plattMPC_ccrhs + 152;
plattMPC_FLOAT plattMPC_Phi09[8];
plattMPC_FLOAT plattMPC_W09[8];
plattMPC_FLOAT plattMPC_Ysd09[36];
plattMPC_FLOAT plattMPC_Lsd09[36];
plattMPC_FLOAT* plattMPC_z10 = plattMPC_z + 80;
plattMPC_FLOAT* plattMPC_dzaff10 = plattMPC_dz_aff + 80;
plattMPC_FLOAT* plattMPC_dzcc10 = plattMPC_dz_cc + 80;
plattMPC_FLOAT* plattMPC_rd10 = plattMPC_rd + 80;
plattMPC_FLOAT plattMPC_Lbyrd10[8];
plattMPC_FLOAT* plattMPC_grad_cost10 = plattMPC_grad_cost + 80;
plattMPC_FLOAT* plattMPC_grad_eq10 = plattMPC_grad_eq + 80;
plattMPC_FLOAT* plattMPC_grad_ineq10 = plattMPC_grad_ineq + 80;
plattMPC_FLOAT plattMPC_ctv10[8];
plattMPC_FLOAT* plattMPC_v10 = plattMPC_v + 60;
plattMPC_FLOAT plattMPC_re10[6];
plattMPC_FLOAT plattMPC_beta10[6];
plattMPC_FLOAT plattMPC_betacc10[6];
plattMPC_FLOAT* plattMPC_dvaff10 = plattMPC_dv_aff + 60;
plattMPC_FLOAT* plattMPC_dvcc10 = plattMPC_dv_cc + 60;
plattMPC_FLOAT plattMPC_V10[48];
plattMPC_FLOAT plattMPC_Yd10[21];
plattMPC_FLOAT plattMPC_Ld10[21];
plattMPC_FLOAT plattMPC_yy10[6];
plattMPC_FLOAT plattMPC_bmy10[6];
plattMPC_FLOAT plattMPC_c10[6] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
int plattMPC_lbIdx10[8] = {0, 1, 2, 3, 4, 5, 6, 7};
plattMPC_FLOAT* plattMPC_llb10 = plattMPC_l + 160;
plattMPC_FLOAT* plattMPC_slb10 = plattMPC_s + 160;
plattMPC_FLOAT* plattMPC_llbbyslb10 = plattMPC_lbys + 160;
plattMPC_FLOAT plattMPC_rilb10[8];
plattMPC_FLOAT* plattMPC_dllbaff10 = plattMPC_dl_aff + 160;
plattMPC_FLOAT* plattMPC_dslbaff10 = plattMPC_ds_aff + 160;
plattMPC_FLOAT* plattMPC_dllbcc10 = plattMPC_dl_cc + 160;
plattMPC_FLOAT* plattMPC_dslbcc10 = plattMPC_ds_cc + 160;
plattMPC_FLOAT* plattMPC_ccrhsl10 = plattMPC_ccrhs + 160;
int plattMPC_ubIdx10[8] = {0, 1, 2, 3, 4, 5, 6, 7};
plattMPC_FLOAT* plattMPC_lub10 = plattMPC_l + 168;
plattMPC_FLOAT* plattMPC_sub10 = plattMPC_s + 168;
plattMPC_FLOAT* plattMPC_lubbysub10 = plattMPC_lbys + 168;
plattMPC_FLOAT plattMPC_riub10[8];
plattMPC_FLOAT* plattMPC_dlubaff10 = plattMPC_dl_aff + 168;
plattMPC_FLOAT* plattMPC_dsubaff10 = plattMPC_ds_aff + 168;
plattMPC_FLOAT* plattMPC_dlubcc10 = plattMPC_dl_cc + 168;
plattMPC_FLOAT* plattMPC_dsubcc10 = plattMPC_ds_cc + 168;
plattMPC_FLOAT* plattMPC_ccrhsub10 = plattMPC_ccrhs + 168;
plattMPC_FLOAT plattMPC_Phi10[8];
plattMPC_FLOAT plattMPC_W10[8];
plattMPC_FLOAT plattMPC_Ysd10[36];
plattMPC_FLOAT plattMPC_Lsd10[36];
plattMPC_FLOAT* plattMPC_z11 = plattMPC_z + 88;
plattMPC_FLOAT* plattMPC_dzaff11 = plattMPC_dz_aff + 88;
plattMPC_FLOAT* plattMPC_dzcc11 = plattMPC_dz_cc + 88;
plattMPC_FLOAT* plattMPC_rd11 = plattMPC_rd + 88;
plattMPC_FLOAT plattMPC_Lbyrd11[8];
plattMPC_FLOAT* plattMPC_grad_cost11 = plattMPC_grad_cost + 88;
plattMPC_FLOAT* plattMPC_grad_eq11 = plattMPC_grad_eq + 88;
plattMPC_FLOAT* plattMPC_grad_ineq11 = plattMPC_grad_ineq + 88;
plattMPC_FLOAT plattMPC_ctv11[8];
plattMPC_FLOAT* plattMPC_v11 = plattMPC_v + 66;
plattMPC_FLOAT plattMPC_re11[6];
plattMPC_FLOAT plattMPC_beta11[6];
plattMPC_FLOAT plattMPC_betacc11[6];
plattMPC_FLOAT* plattMPC_dvaff11 = plattMPC_dv_aff + 66;
plattMPC_FLOAT* plattMPC_dvcc11 = plattMPC_dv_cc + 66;
plattMPC_FLOAT plattMPC_V11[48];
plattMPC_FLOAT plattMPC_Yd11[21];
plattMPC_FLOAT plattMPC_Ld11[21];
plattMPC_FLOAT plattMPC_yy11[6];
plattMPC_FLOAT plattMPC_bmy11[6];
plattMPC_FLOAT plattMPC_c11[6] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
int plattMPC_lbIdx11[8] = {0, 1, 2, 3, 4, 5, 6, 7};
plattMPC_FLOAT* plattMPC_llb11 = plattMPC_l + 176;
plattMPC_FLOAT* plattMPC_slb11 = plattMPC_s + 176;
plattMPC_FLOAT* plattMPC_llbbyslb11 = plattMPC_lbys + 176;
plattMPC_FLOAT plattMPC_rilb11[8];
plattMPC_FLOAT* plattMPC_dllbaff11 = plattMPC_dl_aff + 176;
plattMPC_FLOAT* plattMPC_dslbaff11 = plattMPC_ds_aff + 176;
plattMPC_FLOAT* plattMPC_dllbcc11 = plattMPC_dl_cc + 176;
plattMPC_FLOAT* plattMPC_dslbcc11 = plattMPC_ds_cc + 176;
plattMPC_FLOAT* plattMPC_ccrhsl11 = plattMPC_ccrhs + 176;
int plattMPC_ubIdx11[8] = {0, 1, 2, 3, 4, 5, 6, 7};
plattMPC_FLOAT* plattMPC_lub11 = plattMPC_l + 184;
plattMPC_FLOAT* plattMPC_sub11 = plattMPC_s + 184;
plattMPC_FLOAT* plattMPC_lubbysub11 = plattMPC_lbys + 184;
plattMPC_FLOAT plattMPC_riub11[8];
plattMPC_FLOAT* plattMPC_dlubaff11 = plattMPC_dl_aff + 184;
plattMPC_FLOAT* plattMPC_dsubaff11 = plattMPC_ds_aff + 184;
plattMPC_FLOAT* plattMPC_dlubcc11 = plattMPC_dl_cc + 184;
plattMPC_FLOAT* plattMPC_dsubcc11 = plattMPC_ds_cc + 184;
plattMPC_FLOAT* plattMPC_ccrhsub11 = plattMPC_ccrhs + 184;
plattMPC_FLOAT plattMPC_Phi11[8];
plattMPC_FLOAT plattMPC_W11[8];
plattMPC_FLOAT plattMPC_Ysd11[36];
plattMPC_FLOAT plattMPC_Lsd11[36];
plattMPC_FLOAT* plattMPC_z12 = plattMPC_z + 96;
plattMPC_FLOAT* plattMPC_dzaff12 = plattMPC_dz_aff + 96;
plattMPC_FLOAT* plattMPC_dzcc12 = plattMPC_dz_cc + 96;
plattMPC_FLOAT* plattMPC_rd12 = plattMPC_rd + 96;
plattMPC_FLOAT plattMPC_Lbyrd12[8];
plattMPC_FLOAT* plattMPC_grad_cost12 = plattMPC_grad_cost + 96;
plattMPC_FLOAT* plattMPC_grad_eq12 = plattMPC_grad_eq + 96;
plattMPC_FLOAT* plattMPC_grad_ineq12 = plattMPC_grad_ineq + 96;
plattMPC_FLOAT plattMPC_ctv12[8];
plattMPC_FLOAT* plattMPC_v12 = plattMPC_v + 72;
plattMPC_FLOAT plattMPC_re12[6];
plattMPC_FLOAT plattMPC_beta12[6];
plattMPC_FLOAT plattMPC_betacc12[6];
plattMPC_FLOAT* plattMPC_dvaff12 = plattMPC_dv_aff + 72;
plattMPC_FLOAT* plattMPC_dvcc12 = plattMPC_dv_cc + 72;
plattMPC_FLOAT plattMPC_V12[48];
plattMPC_FLOAT plattMPC_Yd12[21];
plattMPC_FLOAT plattMPC_Ld12[21];
plattMPC_FLOAT plattMPC_yy12[6];
plattMPC_FLOAT plattMPC_bmy12[6];
plattMPC_FLOAT plattMPC_c12[6] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
int plattMPC_lbIdx12[8] = {0, 1, 2, 3, 4, 5, 6, 7};
plattMPC_FLOAT* plattMPC_llb12 = plattMPC_l + 192;
plattMPC_FLOAT* plattMPC_slb12 = plattMPC_s + 192;
plattMPC_FLOAT* plattMPC_llbbyslb12 = plattMPC_lbys + 192;
plattMPC_FLOAT plattMPC_rilb12[8];
plattMPC_FLOAT* plattMPC_dllbaff12 = plattMPC_dl_aff + 192;
plattMPC_FLOAT* plattMPC_dslbaff12 = plattMPC_ds_aff + 192;
plattMPC_FLOAT* plattMPC_dllbcc12 = plattMPC_dl_cc + 192;
plattMPC_FLOAT* plattMPC_dslbcc12 = plattMPC_ds_cc + 192;
plattMPC_FLOAT* plattMPC_ccrhsl12 = plattMPC_ccrhs + 192;
int plattMPC_ubIdx12[8] = {0, 1, 2, 3, 4, 5, 6, 7};
plattMPC_FLOAT* plattMPC_lub12 = plattMPC_l + 200;
plattMPC_FLOAT* plattMPC_sub12 = plattMPC_s + 200;
plattMPC_FLOAT* plattMPC_lubbysub12 = plattMPC_lbys + 200;
plattMPC_FLOAT plattMPC_riub12[8];
plattMPC_FLOAT* plattMPC_dlubaff12 = plattMPC_dl_aff + 200;
plattMPC_FLOAT* plattMPC_dsubaff12 = plattMPC_ds_aff + 200;
plattMPC_FLOAT* plattMPC_dlubcc12 = plattMPC_dl_cc + 200;
plattMPC_FLOAT* plattMPC_dsubcc12 = plattMPC_ds_cc + 200;
plattMPC_FLOAT* plattMPC_ccrhsub12 = plattMPC_ccrhs + 200;
plattMPC_FLOAT plattMPC_Phi12[8];
plattMPC_FLOAT plattMPC_W12[8];
plattMPC_FLOAT plattMPC_Ysd12[36];
plattMPC_FLOAT plattMPC_Lsd12[36];
plattMPC_FLOAT* plattMPC_z13 = plattMPC_z + 104;
plattMPC_FLOAT* plattMPC_dzaff13 = plattMPC_dz_aff + 104;
plattMPC_FLOAT* plattMPC_dzcc13 = plattMPC_dz_cc + 104;
plattMPC_FLOAT* plattMPC_rd13 = plattMPC_rd + 104;
plattMPC_FLOAT plattMPC_Lbyrd13[8];
plattMPC_FLOAT* plattMPC_grad_cost13 = plattMPC_grad_cost + 104;
plattMPC_FLOAT* plattMPC_grad_eq13 = plattMPC_grad_eq + 104;
plattMPC_FLOAT* plattMPC_grad_ineq13 = plattMPC_grad_ineq + 104;
plattMPC_FLOAT plattMPC_ctv13[8];
plattMPC_FLOAT* plattMPC_v13 = plattMPC_v + 78;
plattMPC_FLOAT plattMPC_re13[6];
plattMPC_FLOAT plattMPC_beta13[6];
plattMPC_FLOAT plattMPC_betacc13[6];
plattMPC_FLOAT* plattMPC_dvaff13 = plattMPC_dv_aff + 78;
plattMPC_FLOAT* plattMPC_dvcc13 = plattMPC_dv_cc + 78;
plattMPC_FLOAT plattMPC_V13[48];
plattMPC_FLOAT plattMPC_Yd13[21];
plattMPC_FLOAT plattMPC_Ld13[21];
plattMPC_FLOAT plattMPC_yy13[6];
plattMPC_FLOAT plattMPC_bmy13[6];
plattMPC_FLOAT plattMPC_c13[6] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
int plattMPC_lbIdx13[8] = {0, 1, 2, 3, 4, 5, 6, 7};
plattMPC_FLOAT* plattMPC_llb13 = plattMPC_l + 208;
plattMPC_FLOAT* plattMPC_slb13 = plattMPC_s + 208;
plattMPC_FLOAT* plattMPC_llbbyslb13 = plattMPC_lbys + 208;
plattMPC_FLOAT plattMPC_rilb13[8];
plattMPC_FLOAT* plattMPC_dllbaff13 = plattMPC_dl_aff + 208;
plattMPC_FLOAT* plattMPC_dslbaff13 = plattMPC_ds_aff + 208;
plattMPC_FLOAT* plattMPC_dllbcc13 = plattMPC_dl_cc + 208;
plattMPC_FLOAT* plattMPC_dslbcc13 = plattMPC_ds_cc + 208;
plattMPC_FLOAT* plattMPC_ccrhsl13 = plattMPC_ccrhs + 208;
int plattMPC_ubIdx13[8] = {0, 1, 2, 3, 4, 5, 6, 7};
plattMPC_FLOAT* plattMPC_lub13 = plattMPC_l + 216;
plattMPC_FLOAT* plattMPC_sub13 = plattMPC_s + 216;
plattMPC_FLOAT* plattMPC_lubbysub13 = plattMPC_lbys + 216;
plattMPC_FLOAT plattMPC_riub13[8];
plattMPC_FLOAT* plattMPC_dlubaff13 = plattMPC_dl_aff + 216;
plattMPC_FLOAT* plattMPC_dsubaff13 = plattMPC_ds_aff + 216;
plattMPC_FLOAT* plattMPC_dlubcc13 = plattMPC_dl_cc + 216;
plattMPC_FLOAT* plattMPC_dsubcc13 = plattMPC_ds_cc + 216;
plattMPC_FLOAT* plattMPC_ccrhsub13 = plattMPC_ccrhs + 216;
plattMPC_FLOAT plattMPC_Phi13[8];
plattMPC_FLOAT plattMPC_W13[8];
plattMPC_FLOAT plattMPC_Ysd13[36];
plattMPC_FLOAT plattMPC_Lsd13[36];
plattMPC_FLOAT* plattMPC_z14 = plattMPC_z + 112;
plattMPC_FLOAT* plattMPC_dzaff14 = plattMPC_dz_aff + 112;
plattMPC_FLOAT* plattMPC_dzcc14 = plattMPC_dz_cc + 112;
plattMPC_FLOAT* plattMPC_rd14 = plattMPC_rd + 112;
plattMPC_FLOAT plattMPC_Lbyrd14[6];
plattMPC_FLOAT* plattMPC_grad_cost14 = plattMPC_grad_cost + 112;
plattMPC_FLOAT* plattMPC_grad_eq14 = plattMPC_grad_eq + 112;
plattMPC_FLOAT* plattMPC_grad_ineq14 = plattMPC_grad_ineq + 112;
plattMPC_FLOAT plattMPC_ctv14[6];
plattMPC_FLOAT* plattMPC_v14 = plattMPC_v + 84;
plattMPC_FLOAT plattMPC_re14[6];
plattMPC_FLOAT plattMPC_beta14[6];
plattMPC_FLOAT plattMPC_betacc14[6];
plattMPC_FLOAT* plattMPC_dvaff14 = plattMPC_dv_aff + 84;
plattMPC_FLOAT* plattMPC_dvcc14 = plattMPC_dv_cc + 84;
plattMPC_FLOAT plattMPC_V14[36];
plattMPC_FLOAT plattMPC_Yd14[21];
plattMPC_FLOAT plattMPC_Ld14[21];
plattMPC_FLOAT plattMPC_yy14[6];
plattMPC_FLOAT plattMPC_bmy14[6];
plattMPC_FLOAT plattMPC_c14[6] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
int plattMPC_lbIdx14[6] = {0, 1, 2, 3, 4, 5};
plattMPC_FLOAT* plattMPC_llb14 = plattMPC_l + 224;
plattMPC_FLOAT* plattMPC_slb14 = plattMPC_s + 224;
plattMPC_FLOAT* plattMPC_llbbyslb14 = plattMPC_lbys + 224;
plattMPC_FLOAT plattMPC_rilb14[6];
plattMPC_FLOAT* plattMPC_dllbaff14 = plattMPC_dl_aff + 224;
plattMPC_FLOAT* plattMPC_dslbaff14 = plattMPC_ds_aff + 224;
plattMPC_FLOAT* plattMPC_dllbcc14 = plattMPC_dl_cc + 224;
plattMPC_FLOAT* plattMPC_dslbcc14 = plattMPC_ds_cc + 224;
plattMPC_FLOAT* plattMPC_ccrhsl14 = plattMPC_ccrhs + 224;
int plattMPC_ubIdx14[6] = {0, 1, 2, 3, 4, 5};
plattMPC_FLOAT* plattMPC_lub14 = plattMPC_l + 230;
plattMPC_FLOAT* plattMPC_sub14 = plattMPC_s + 230;
plattMPC_FLOAT* plattMPC_lubbysub14 = plattMPC_lbys + 230;
plattMPC_FLOAT plattMPC_riub14[6];
plattMPC_FLOAT* plattMPC_dlubaff14 = plattMPC_dl_aff + 230;
plattMPC_FLOAT* plattMPC_dsubaff14 = plattMPC_ds_aff + 230;
plattMPC_FLOAT* plattMPC_dlubcc14 = plattMPC_dl_cc + 230;
plattMPC_FLOAT* plattMPC_dsubcc14 = plattMPC_ds_cc + 230;
plattMPC_FLOAT* plattMPC_ccrhsub14 = plattMPC_ccrhs + 230;
plattMPC_FLOAT plattMPC_Phi14[6];
plattMPC_FLOAT plattMPC_D14[6] = {-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000};
plattMPC_FLOAT plattMPC_W14[6];
plattMPC_FLOAT plattMPC_Ysd14[36];
plattMPC_FLOAT plattMPC_Lsd14[36];
plattMPC_FLOAT musigma;
plattMPC_FLOAT sigma_3rdroot;
plattMPC_FLOAT plattMPC_Diag1_0[8];
plattMPC_FLOAT plattMPC_Diag2_0[8];
plattMPC_FLOAT plattMPC_L_0[28];




/* SOLVER CODE --------------------------------------------------------- */
int plattMPC_solve(plattMPC_params* params, plattMPC_output* output, plattMPC_info* info)
{	
int exitcode;

#if plattMPC_SET_TIMING == 1
	plattMPC_timer solvertimer;
	plattMPC_tic(&solvertimer);
#endif
/* FUNCTION CALLS INTO LA LIBRARY -------------------------------------- */
info->it = 0;
plattMPC_LA_INITIALIZEVECTOR_118(plattMPC_z, 0);
plattMPC_LA_INITIALIZEVECTOR_90(plattMPC_v, 1);
plattMPC_LA_INITIALIZEVECTOR_236(plattMPC_l, 1);
plattMPC_LA_INITIALIZEVECTOR_236(plattMPC_s, 1);
info->mu = 0;
plattMPC_LA_DOTACC_236(plattMPC_l, plattMPC_s, &info->mu);
info->mu /= 236;
while( 1 ){
info->pobj = 0;
plattMPC_LA_DIAG_QUADFCN_8(params->H1, params->f1, plattMPC_z00, plattMPC_grad_cost00, &info->pobj);
plattMPC_LA_DIAG_QUADFCN_8(params->H2, params->f2, plattMPC_z01, plattMPC_grad_cost01, &info->pobj);
plattMPC_LA_DIAG_QUADFCN_8(params->H3, params->f3, plattMPC_z02, plattMPC_grad_cost02, &info->pobj);
plattMPC_LA_DIAG_QUADFCN_8(params->H4, params->f4, plattMPC_z03, plattMPC_grad_cost03, &info->pobj);
plattMPC_LA_DIAG_QUADFCN_8(params->H5, params->f5, plattMPC_z04, plattMPC_grad_cost04, &info->pobj);
plattMPC_LA_DIAG_QUADFCN_8(params->H6, params->f6, plattMPC_z05, plattMPC_grad_cost05, &info->pobj);
plattMPC_LA_DIAG_QUADFCN_8(params->H7, params->f7, plattMPC_z06, plattMPC_grad_cost06, &info->pobj);
plattMPC_LA_DIAG_QUADFCN_8(params->H8, params->f8, plattMPC_z07, plattMPC_grad_cost07, &info->pobj);
plattMPC_LA_DIAG_QUADFCN_8(params->H9, params->f9, plattMPC_z08, plattMPC_grad_cost08, &info->pobj);
plattMPC_LA_DIAG_QUADFCN_8(params->H10, params->f10, plattMPC_z09, plattMPC_grad_cost09, &info->pobj);
plattMPC_LA_DIAG_QUADFCN_8(params->H11, params->f11, plattMPC_z10, plattMPC_grad_cost10, &info->pobj);
plattMPC_LA_DIAG_QUADFCN_8(params->H12, params->f12, plattMPC_z11, plattMPC_grad_cost11, &info->pobj);
plattMPC_LA_DIAG_QUADFCN_8(params->H13, params->f13, plattMPC_z12, plattMPC_grad_cost12, &info->pobj);
plattMPC_LA_DIAG_QUADFCN_8(params->H14, params->f14, plattMPC_z13, plattMPC_grad_cost13, &info->pobj);
plattMPC_LA_DIAG_QUADFCN_6(params->H15, params->f15, plattMPC_z14, plattMPC_grad_cost14, &info->pobj);
info->res_eq = 0;
info->dgap = 0;
plattMPC_LA_DIAGZERO_MVMSUB6_6(plattMPC_D00, plattMPC_z00, plattMPC_c00, plattMPC_v00, plattMPC_re00, &info->dgap, &info->res_eq);
plattMPC_LA_DENSE_DIAGZERO_MVMSUB3_6_8_8(plattMPC_C00, plattMPC_z00, plattMPC_D01, plattMPC_z01, plattMPC_c01, plattMPC_v01, plattMPC_re01, &info->dgap, &info->res_eq);
plattMPC_LA_DENSE_DIAGZERO_MVMSUB3_6_8_8(plattMPC_C00, plattMPC_z01, plattMPC_D01, plattMPC_z02, plattMPC_c02, plattMPC_v02, plattMPC_re02, &info->dgap, &info->res_eq);
plattMPC_LA_DENSE_DIAGZERO_MVMSUB3_6_8_8(plattMPC_C00, plattMPC_z02, plattMPC_D01, plattMPC_z03, plattMPC_c03, plattMPC_v03, plattMPC_re03, &info->dgap, &info->res_eq);
plattMPC_LA_DENSE_DIAGZERO_MVMSUB3_6_8_8(plattMPC_C00, plattMPC_z03, plattMPC_D01, plattMPC_z04, plattMPC_c04, plattMPC_v04, plattMPC_re04, &info->dgap, &info->res_eq);
plattMPC_LA_DENSE_DIAGZERO_MVMSUB3_6_8_8(plattMPC_C00, plattMPC_z04, plattMPC_D01, plattMPC_z05, plattMPC_c05, plattMPC_v05, plattMPC_re05, &info->dgap, &info->res_eq);
plattMPC_LA_DENSE_DIAGZERO_MVMSUB3_6_8_8(plattMPC_C00, plattMPC_z05, plattMPC_D01, plattMPC_z06, plattMPC_c06, plattMPC_v06, plattMPC_re06, &info->dgap, &info->res_eq);
plattMPC_LA_DENSE_DIAGZERO_MVMSUB3_6_8_8(plattMPC_C00, plattMPC_z06, plattMPC_D01, plattMPC_z07, plattMPC_c07, plattMPC_v07, plattMPC_re07, &info->dgap, &info->res_eq);
plattMPC_LA_DENSE_DIAGZERO_MVMSUB3_6_8_8(plattMPC_C00, plattMPC_z07, plattMPC_D01, plattMPC_z08, plattMPC_c08, plattMPC_v08, plattMPC_re08, &info->dgap, &info->res_eq);
plattMPC_LA_DENSE_DIAGZERO_MVMSUB3_6_8_8(plattMPC_C00, plattMPC_z08, plattMPC_D01, plattMPC_z09, plattMPC_c09, plattMPC_v09, plattMPC_re09, &info->dgap, &info->res_eq);
plattMPC_LA_DENSE_DIAGZERO_MVMSUB3_6_8_8(plattMPC_C00, plattMPC_z09, plattMPC_D01, plattMPC_z10, plattMPC_c10, plattMPC_v10, plattMPC_re10, &info->dgap, &info->res_eq);
plattMPC_LA_DENSE_DIAGZERO_MVMSUB3_6_8_8(plattMPC_C00, plattMPC_z10, plattMPC_D01, plattMPC_z11, plattMPC_c11, plattMPC_v11, plattMPC_re11, &info->dgap, &info->res_eq);
plattMPC_LA_DENSE_DIAGZERO_MVMSUB3_6_8_8(plattMPC_C00, plattMPC_z11, plattMPC_D01, plattMPC_z12, plattMPC_c12, plattMPC_v12, plattMPC_re12, &info->dgap, &info->res_eq);
plattMPC_LA_DENSE_DIAGZERO_MVMSUB3_6_8_8(plattMPC_C00, plattMPC_z12, plattMPC_D01, plattMPC_z13, plattMPC_c13, plattMPC_v13, plattMPC_re13, &info->dgap, &info->res_eq);
plattMPC_LA_DENSE_DIAGZERO_MVMSUB3_6_8_6(plattMPC_C00, plattMPC_z13, plattMPC_D14, plattMPC_z14, plattMPC_c14, plattMPC_v14, plattMPC_re14, &info->dgap, &info->res_eq);
plattMPC_LA_DENSE_DIAGZERO_MTVM2_6_8_6(plattMPC_C00, plattMPC_v01, plattMPC_D00, plattMPC_v00, plattMPC_grad_eq00);
plattMPC_LA_DENSE_DIAGZERO_MTVM2_6_8_6(plattMPC_C00, plattMPC_v02, plattMPC_D01, plattMPC_v01, plattMPC_grad_eq01);
plattMPC_LA_DENSE_DIAGZERO_MTVM2_6_8_6(plattMPC_C00, plattMPC_v03, plattMPC_D01, plattMPC_v02, plattMPC_grad_eq02);
plattMPC_LA_DENSE_DIAGZERO_MTVM2_6_8_6(plattMPC_C00, plattMPC_v04, plattMPC_D01, plattMPC_v03, plattMPC_grad_eq03);
plattMPC_LA_DENSE_DIAGZERO_MTVM2_6_8_6(plattMPC_C00, plattMPC_v05, plattMPC_D01, plattMPC_v04, plattMPC_grad_eq04);
plattMPC_LA_DENSE_DIAGZERO_MTVM2_6_8_6(plattMPC_C00, plattMPC_v06, plattMPC_D01, plattMPC_v05, plattMPC_grad_eq05);
plattMPC_LA_DENSE_DIAGZERO_MTVM2_6_8_6(plattMPC_C00, plattMPC_v07, plattMPC_D01, plattMPC_v06, plattMPC_grad_eq06);
plattMPC_LA_DENSE_DIAGZERO_MTVM2_6_8_6(plattMPC_C00, plattMPC_v08, plattMPC_D01, plattMPC_v07, plattMPC_grad_eq07);
plattMPC_LA_DENSE_DIAGZERO_MTVM2_6_8_6(plattMPC_C00, plattMPC_v09, plattMPC_D01, plattMPC_v08, plattMPC_grad_eq08);
plattMPC_LA_DENSE_DIAGZERO_MTVM2_6_8_6(plattMPC_C00, plattMPC_v10, plattMPC_D01, plattMPC_v09, plattMPC_grad_eq09);
plattMPC_LA_DENSE_DIAGZERO_MTVM2_6_8_6(plattMPC_C00, plattMPC_v11, plattMPC_D01, plattMPC_v10, plattMPC_grad_eq10);
plattMPC_LA_DENSE_DIAGZERO_MTVM2_6_8_6(plattMPC_C00, plattMPC_v12, plattMPC_D01, plattMPC_v11, plattMPC_grad_eq11);
plattMPC_LA_DENSE_DIAGZERO_MTVM2_6_8_6(plattMPC_C00, plattMPC_v13, plattMPC_D01, plattMPC_v12, plattMPC_grad_eq12);
plattMPC_LA_DENSE_DIAGZERO_MTVM2_6_8_6(plattMPC_C00, plattMPC_v14, plattMPC_D01, plattMPC_v13, plattMPC_grad_eq13);
plattMPC_LA_DIAGZERO_MTVM_6_6(plattMPC_D14, plattMPC_v14, plattMPC_grad_eq14);
info->res_ineq = 0;
plattMPC_LA_VSUBADD3_8(params->lb1, plattMPC_z00, plattMPC_lbIdx00, plattMPC_llb00, plattMPC_slb00, plattMPC_rilb00, &info->dgap, &info->res_ineq);
plattMPC_LA_VSUBADD2_8(plattMPC_z00, plattMPC_ubIdx00, params->ub1, plattMPC_lub00, plattMPC_sub00, plattMPC_riub00, &info->dgap, &info->res_ineq);
plattMPC_LA_VSUBADD3_8(params->lb2, plattMPC_z01, plattMPC_lbIdx01, plattMPC_llb01, plattMPC_slb01, plattMPC_rilb01, &info->dgap, &info->res_ineq);
plattMPC_LA_VSUBADD2_8(plattMPC_z01, plattMPC_ubIdx01, params->ub2, plattMPC_lub01, plattMPC_sub01, plattMPC_riub01, &info->dgap, &info->res_ineq);
plattMPC_LA_VSUBADD3_8(params->lb3, plattMPC_z02, plattMPC_lbIdx02, plattMPC_llb02, plattMPC_slb02, plattMPC_rilb02, &info->dgap, &info->res_ineq);
plattMPC_LA_VSUBADD2_8(plattMPC_z02, plattMPC_ubIdx02, params->ub3, plattMPC_lub02, plattMPC_sub02, plattMPC_riub02, &info->dgap, &info->res_ineq);
plattMPC_LA_VSUBADD3_8(params->lb4, plattMPC_z03, plattMPC_lbIdx03, plattMPC_llb03, plattMPC_slb03, plattMPC_rilb03, &info->dgap, &info->res_ineq);
plattMPC_LA_VSUBADD2_8(plattMPC_z03, plattMPC_ubIdx03, params->ub4, plattMPC_lub03, plattMPC_sub03, plattMPC_riub03, &info->dgap, &info->res_ineq);
plattMPC_LA_VSUBADD3_8(params->lb5, plattMPC_z04, plattMPC_lbIdx04, plattMPC_llb04, plattMPC_slb04, plattMPC_rilb04, &info->dgap, &info->res_ineq);
plattMPC_LA_VSUBADD2_8(plattMPC_z04, plattMPC_ubIdx04, params->ub5, plattMPC_lub04, plattMPC_sub04, plattMPC_riub04, &info->dgap, &info->res_ineq);
plattMPC_LA_VSUBADD3_8(params->lb6, plattMPC_z05, plattMPC_lbIdx05, plattMPC_llb05, plattMPC_slb05, plattMPC_rilb05, &info->dgap, &info->res_ineq);
plattMPC_LA_VSUBADD2_8(plattMPC_z05, plattMPC_ubIdx05, params->ub6, plattMPC_lub05, plattMPC_sub05, plattMPC_riub05, &info->dgap, &info->res_ineq);
plattMPC_LA_VSUBADD3_8(params->lb7, plattMPC_z06, plattMPC_lbIdx06, plattMPC_llb06, plattMPC_slb06, plattMPC_rilb06, &info->dgap, &info->res_ineq);
plattMPC_LA_VSUBADD2_8(plattMPC_z06, plattMPC_ubIdx06, params->ub7, plattMPC_lub06, plattMPC_sub06, plattMPC_riub06, &info->dgap, &info->res_ineq);
plattMPC_LA_VSUBADD3_8(params->lb8, plattMPC_z07, plattMPC_lbIdx07, plattMPC_llb07, plattMPC_slb07, plattMPC_rilb07, &info->dgap, &info->res_ineq);
plattMPC_LA_VSUBADD2_8(plattMPC_z07, plattMPC_ubIdx07, params->ub8, plattMPC_lub07, plattMPC_sub07, plattMPC_riub07, &info->dgap, &info->res_ineq);
plattMPC_LA_VSUBADD3_8(params->lb9, plattMPC_z08, plattMPC_lbIdx08, plattMPC_llb08, plattMPC_slb08, plattMPC_rilb08, &info->dgap, &info->res_ineq);
plattMPC_LA_VSUBADD2_8(plattMPC_z08, plattMPC_ubIdx08, params->ub9, plattMPC_lub08, plattMPC_sub08, plattMPC_riub08, &info->dgap, &info->res_ineq);
plattMPC_LA_VSUBADD3_8(params->lb10, plattMPC_z09, plattMPC_lbIdx09, plattMPC_llb09, plattMPC_slb09, plattMPC_rilb09, &info->dgap, &info->res_ineq);
plattMPC_LA_VSUBADD2_8(plattMPC_z09, plattMPC_ubIdx09, params->ub10, plattMPC_lub09, plattMPC_sub09, plattMPC_riub09, &info->dgap, &info->res_ineq);
plattMPC_LA_VSUBADD3_8(params->lb11, plattMPC_z10, plattMPC_lbIdx10, plattMPC_llb10, plattMPC_slb10, plattMPC_rilb10, &info->dgap, &info->res_ineq);
plattMPC_LA_VSUBADD2_8(plattMPC_z10, plattMPC_ubIdx10, params->ub11, plattMPC_lub10, plattMPC_sub10, plattMPC_riub10, &info->dgap, &info->res_ineq);
plattMPC_LA_VSUBADD3_8(params->lb12, plattMPC_z11, plattMPC_lbIdx11, plattMPC_llb11, plattMPC_slb11, plattMPC_rilb11, &info->dgap, &info->res_ineq);
plattMPC_LA_VSUBADD2_8(plattMPC_z11, plattMPC_ubIdx11, params->ub12, plattMPC_lub11, plattMPC_sub11, plattMPC_riub11, &info->dgap, &info->res_ineq);
plattMPC_LA_VSUBADD3_8(params->lb13, plattMPC_z12, plattMPC_lbIdx12, plattMPC_llb12, plattMPC_slb12, plattMPC_rilb12, &info->dgap, &info->res_ineq);
plattMPC_LA_VSUBADD2_8(plattMPC_z12, plattMPC_ubIdx12, params->ub13, plattMPC_lub12, plattMPC_sub12, plattMPC_riub12, &info->dgap, &info->res_ineq);
plattMPC_LA_VSUBADD3_8(params->lb14, plattMPC_z13, plattMPC_lbIdx13, plattMPC_llb13, plattMPC_slb13, plattMPC_rilb13, &info->dgap, &info->res_ineq);
plattMPC_LA_VSUBADD2_8(plattMPC_z13, plattMPC_ubIdx13, params->ub14, plattMPC_lub13, plattMPC_sub13, plattMPC_riub13, &info->dgap, &info->res_ineq);
plattMPC_LA_VSUBADD3_6(params->lb15, plattMPC_z14, plattMPC_lbIdx14, plattMPC_llb14, plattMPC_slb14, plattMPC_rilb14, &info->dgap, &info->res_ineq);
plattMPC_LA_VSUBADD2_6(plattMPC_z14, plattMPC_ubIdx14, params->ub15, plattMPC_lub14, plattMPC_sub14, plattMPC_riub14, &info->dgap, &info->res_ineq);
plattMPC_LA_INEQ_B_GRAD_8_8_8(plattMPC_lub00, plattMPC_sub00, plattMPC_riub00, plattMPC_llb00, plattMPC_slb00, plattMPC_rilb00, plattMPC_lbIdx00, plattMPC_ubIdx00, plattMPC_grad_ineq00, plattMPC_lubbysub00, plattMPC_llbbyslb00);
plattMPC_LA_INEQ_B_GRAD_8_8_8(plattMPC_lub01, plattMPC_sub01, plattMPC_riub01, plattMPC_llb01, plattMPC_slb01, plattMPC_rilb01, plattMPC_lbIdx01, plattMPC_ubIdx01, plattMPC_grad_ineq01, plattMPC_lubbysub01, plattMPC_llbbyslb01);
plattMPC_LA_INEQ_B_GRAD_8_8_8(plattMPC_lub02, plattMPC_sub02, plattMPC_riub02, plattMPC_llb02, plattMPC_slb02, plattMPC_rilb02, plattMPC_lbIdx02, plattMPC_ubIdx02, plattMPC_grad_ineq02, plattMPC_lubbysub02, plattMPC_llbbyslb02);
plattMPC_LA_INEQ_B_GRAD_8_8_8(plattMPC_lub03, plattMPC_sub03, plattMPC_riub03, plattMPC_llb03, plattMPC_slb03, plattMPC_rilb03, plattMPC_lbIdx03, plattMPC_ubIdx03, plattMPC_grad_ineq03, plattMPC_lubbysub03, plattMPC_llbbyslb03);
plattMPC_LA_INEQ_B_GRAD_8_8_8(plattMPC_lub04, plattMPC_sub04, plattMPC_riub04, plattMPC_llb04, plattMPC_slb04, plattMPC_rilb04, plattMPC_lbIdx04, plattMPC_ubIdx04, plattMPC_grad_ineq04, plattMPC_lubbysub04, plattMPC_llbbyslb04);
plattMPC_LA_INEQ_B_GRAD_8_8_8(plattMPC_lub05, plattMPC_sub05, plattMPC_riub05, plattMPC_llb05, plattMPC_slb05, plattMPC_rilb05, plattMPC_lbIdx05, plattMPC_ubIdx05, plattMPC_grad_ineq05, plattMPC_lubbysub05, plattMPC_llbbyslb05);
plattMPC_LA_INEQ_B_GRAD_8_8_8(plattMPC_lub06, plattMPC_sub06, plattMPC_riub06, plattMPC_llb06, plattMPC_slb06, plattMPC_rilb06, plattMPC_lbIdx06, plattMPC_ubIdx06, plattMPC_grad_ineq06, plattMPC_lubbysub06, plattMPC_llbbyslb06);
plattMPC_LA_INEQ_B_GRAD_8_8_8(plattMPC_lub07, plattMPC_sub07, plattMPC_riub07, plattMPC_llb07, plattMPC_slb07, plattMPC_rilb07, plattMPC_lbIdx07, plattMPC_ubIdx07, plattMPC_grad_ineq07, plattMPC_lubbysub07, plattMPC_llbbyslb07);
plattMPC_LA_INEQ_B_GRAD_8_8_8(plattMPC_lub08, plattMPC_sub08, plattMPC_riub08, plattMPC_llb08, plattMPC_slb08, plattMPC_rilb08, plattMPC_lbIdx08, plattMPC_ubIdx08, plattMPC_grad_ineq08, plattMPC_lubbysub08, plattMPC_llbbyslb08);
plattMPC_LA_INEQ_B_GRAD_8_8_8(plattMPC_lub09, plattMPC_sub09, plattMPC_riub09, plattMPC_llb09, plattMPC_slb09, plattMPC_rilb09, plattMPC_lbIdx09, plattMPC_ubIdx09, plattMPC_grad_ineq09, plattMPC_lubbysub09, plattMPC_llbbyslb09);
plattMPC_LA_INEQ_B_GRAD_8_8_8(plattMPC_lub10, plattMPC_sub10, plattMPC_riub10, plattMPC_llb10, plattMPC_slb10, plattMPC_rilb10, plattMPC_lbIdx10, plattMPC_ubIdx10, plattMPC_grad_ineq10, plattMPC_lubbysub10, plattMPC_llbbyslb10);
plattMPC_LA_INEQ_B_GRAD_8_8_8(plattMPC_lub11, plattMPC_sub11, plattMPC_riub11, plattMPC_llb11, plattMPC_slb11, plattMPC_rilb11, plattMPC_lbIdx11, plattMPC_ubIdx11, plattMPC_grad_ineq11, plattMPC_lubbysub11, plattMPC_llbbyslb11);
plattMPC_LA_INEQ_B_GRAD_8_8_8(plattMPC_lub12, plattMPC_sub12, plattMPC_riub12, plattMPC_llb12, plattMPC_slb12, plattMPC_rilb12, plattMPC_lbIdx12, plattMPC_ubIdx12, plattMPC_grad_ineq12, plattMPC_lubbysub12, plattMPC_llbbyslb12);
plattMPC_LA_INEQ_B_GRAD_8_8_8(plattMPC_lub13, plattMPC_sub13, plattMPC_riub13, plattMPC_llb13, plattMPC_slb13, plattMPC_rilb13, plattMPC_lbIdx13, plattMPC_ubIdx13, plattMPC_grad_ineq13, plattMPC_lubbysub13, plattMPC_llbbyslb13);
plattMPC_LA_INEQ_B_GRAD_6_6_6(plattMPC_lub14, plattMPC_sub14, plattMPC_riub14, plattMPC_llb14, plattMPC_slb14, plattMPC_rilb14, plattMPC_lbIdx14, plattMPC_ubIdx14, plattMPC_grad_ineq14, plattMPC_lubbysub14, plattMPC_llbbyslb14);
info->dobj = info->pobj - info->dgap;
info->rdgap = info->pobj ? info->dgap / info->pobj : 1e6;
if( info->rdgap < 0 ) info->rdgap = -info->rdgap;
if( info->mu < plattMPC_SET_ACC_KKTCOMPL
    && (info->rdgap < plattMPC_SET_ACC_RDGAP || info->dgap < plattMPC_SET_ACC_KKTCOMPL)
    && info->res_eq < plattMPC_SET_ACC_RESEQ
    && info->res_ineq < plattMPC_SET_ACC_RESINEQ ){
exitcode = plattMPC_OPTIMAL; break; }
if( info->it == plattMPC_SET_MAXIT ){
exitcode = plattMPC_MAXITREACHED; break; }
plattMPC_LA_VVADD3_118(plattMPC_grad_cost, plattMPC_grad_eq, plattMPC_grad_ineq, plattMPC_rd);
plattMPC_LA_DIAG_CHOL_ONELOOP_LBUB_8_8_8(params->H1, plattMPC_llbbyslb00, plattMPC_lbIdx00, plattMPC_lubbysub00, plattMPC_ubIdx00, plattMPC_Phi00);
plattMPC_LA_DIAG_MATRIXFORWARDSUB_6_8(plattMPC_Phi00, plattMPC_C00, plattMPC_V00);
plattMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_6_8(plattMPC_Phi00, plattMPC_D00, plattMPC_W00);
plattMPC_LA_DENSE_DIAGZERO_MMTM_6_8_6(plattMPC_W00, plattMPC_V00, plattMPC_Ysd01);
plattMPC_LA_DIAG_FORWARDSUB_8(plattMPC_Phi00, plattMPC_rd00, plattMPC_Lbyrd00);
plattMPC_LA_DIAG_CHOL_ONELOOP_LBUB_8_8_8(params->H2, plattMPC_llbbyslb01, plattMPC_lbIdx01, plattMPC_lubbysub01, plattMPC_ubIdx01, plattMPC_Phi01);
plattMPC_LA_DIAG_MATRIXFORWARDSUB_6_8(plattMPC_Phi01, plattMPC_C00, plattMPC_V01);
plattMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_6_8(plattMPC_Phi01, plattMPC_D01, plattMPC_W01);
plattMPC_LA_DENSE_DIAGZERO_MMTM_6_8_6(plattMPC_W01, plattMPC_V01, plattMPC_Ysd02);
plattMPC_LA_DIAG_FORWARDSUB_8(plattMPC_Phi01, plattMPC_rd01, plattMPC_Lbyrd01);
plattMPC_LA_DIAG_CHOL_ONELOOP_LBUB_8_8_8(params->H3, plattMPC_llbbyslb02, plattMPC_lbIdx02, plattMPC_lubbysub02, plattMPC_ubIdx02, plattMPC_Phi02);
plattMPC_LA_DIAG_MATRIXFORWARDSUB_6_8(plattMPC_Phi02, plattMPC_C00, plattMPC_V02);
plattMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_6_8(plattMPC_Phi02, plattMPC_D01, plattMPC_W02);
plattMPC_LA_DENSE_DIAGZERO_MMTM_6_8_6(plattMPC_W02, plattMPC_V02, plattMPC_Ysd03);
plattMPC_LA_DIAG_FORWARDSUB_8(plattMPC_Phi02, plattMPC_rd02, plattMPC_Lbyrd02);
plattMPC_LA_DIAG_CHOL_ONELOOP_LBUB_8_8_8(params->H4, plattMPC_llbbyslb03, plattMPC_lbIdx03, plattMPC_lubbysub03, plattMPC_ubIdx03, plattMPC_Phi03);
plattMPC_LA_DIAG_MATRIXFORWARDSUB_6_8(plattMPC_Phi03, plattMPC_C00, plattMPC_V03);
plattMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_6_8(plattMPC_Phi03, plattMPC_D01, plattMPC_W03);
plattMPC_LA_DENSE_DIAGZERO_MMTM_6_8_6(plattMPC_W03, plattMPC_V03, plattMPC_Ysd04);
plattMPC_LA_DIAG_FORWARDSUB_8(plattMPC_Phi03, plattMPC_rd03, plattMPC_Lbyrd03);
plattMPC_LA_DIAG_CHOL_ONELOOP_LBUB_8_8_8(params->H5, plattMPC_llbbyslb04, plattMPC_lbIdx04, plattMPC_lubbysub04, plattMPC_ubIdx04, plattMPC_Phi04);
plattMPC_LA_DIAG_MATRIXFORWARDSUB_6_8(plattMPC_Phi04, plattMPC_C00, plattMPC_V04);
plattMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_6_8(plattMPC_Phi04, plattMPC_D01, plattMPC_W04);
plattMPC_LA_DENSE_DIAGZERO_MMTM_6_8_6(plattMPC_W04, plattMPC_V04, plattMPC_Ysd05);
plattMPC_LA_DIAG_FORWARDSUB_8(plattMPC_Phi04, plattMPC_rd04, plattMPC_Lbyrd04);
plattMPC_LA_DIAG_CHOL_ONELOOP_LBUB_8_8_8(params->H6, plattMPC_llbbyslb05, plattMPC_lbIdx05, plattMPC_lubbysub05, plattMPC_ubIdx05, plattMPC_Phi05);
plattMPC_LA_DIAG_MATRIXFORWARDSUB_6_8(plattMPC_Phi05, plattMPC_C00, plattMPC_V05);
plattMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_6_8(plattMPC_Phi05, plattMPC_D01, plattMPC_W05);
plattMPC_LA_DENSE_DIAGZERO_MMTM_6_8_6(plattMPC_W05, plattMPC_V05, plattMPC_Ysd06);
plattMPC_LA_DIAG_FORWARDSUB_8(plattMPC_Phi05, plattMPC_rd05, plattMPC_Lbyrd05);
plattMPC_LA_DIAG_CHOL_ONELOOP_LBUB_8_8_8(params->H7, plattMPC_llbbyslb06, plattMPC_lbIdx06, plattMPC_lubbysub06, plattMPC_ubIdx06, plattMPC_Phi06);
plattMPC_LA_DIAG_MATRIXFORWARDSUB_6_8(plattMPC_Phi06, plattMPC_C00, plattMPC_V06);
plattMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_6_8(plattMPC_Phi06, plattMPC_D01, plattMPC_W06);
plattMPC_LA_DENSE_DIAGZERO_MMTM_6_8_6(plattMPC_W06, plattMPC_V06, plattMPC_Ysd07);
plattMPC_LA_DIAG_FORWARDSUB_8(plattMPC_Phi06, plattMPC_rd06, plattMPC_Lbyrd06);
plattMPC_LA_DIAG_CHOL_ONELOOP_LBUB_8_8_8(params->H8, plattMPC_llbbyslb07, plattMPC_lbIdx07, plattMPC_lubbysub07, plattMPC_ubIdx07, plattMPC_Phi07);
plattMPC_LA_DIAG_MATRIXFORWARDSUB_6_8(plattMPC_Phi07, plattMPC_C00, plattMPC_V07);
plattMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_6_8(plattMPC_Phi07, plattMPC_D01, plattMPC_W07);
plattMPC_LA_DENSE_DIAGZERO_MMTM_6_8_6(plattMPC_W07, plattMPC_V07, plattMPC_Ysd08);
plattMPC_LA_DIAG_FORWARDSUB_8(plattMPC_Phi07, plattMPC_rd07, plattMPC_Lbyrd07);
plattMPC_LA_DIAG_CHOL_ONELOOP_LBUB_8_8_8(params->H9, plattMPC_llbbyslb08, plattMPC_lbIdx08, plattMPC_lubbysub08, plattMPC_ubIdx08, plattMPC_Phi08);
plattMPC_LA_DIAG_MATRIXFORWARDSUB_6_8(plattMPC_Phi08, plattMPC_C00, plattMPC_V08);
plattMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_6_8(plattMPC_Phi08, plattMPC_D01, plattMPC_W08);
plattMPC_LA_DENSE_DIAGZERO_MMTM_6_8_6(plattMPC_W08, plattMPC_V08, plattMPC_Ysd09);
plattMPC_LA_DIAG_FORWARDSUB_8(plattMPC_Phi08, plattMPC_rd08, plattMPC_Lbyrd08);
plattMPC_LA_DIAG_CHOL_ONELOOP_LBUB_8_8_8(params->H10, plattMPC_llbbyslb09, plattMPC_lbIdx09, plattMPC_lubbysub09, plattMPC_ubIdx09, plattMPC_Phi09);
plattMPC_LA_DIAG_MATRIXFORWARDSUB_6_8(plattMPC_Phi09, plattMPC_C00, plattMPC_V09);
plattMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_6_8(plattMPC_Phi09, plattMPC_D01, plattMPC_W09);
plattMPC_LA_DENSE_DIAGZERO_MMTM_6_8_6(plattMPC_W09, plattMPC_V09, plattMPC_Ysd10);
plattMPC_LA_DIAG_FORWARDSUB_8(plattMPC_Phi09, plattMPC_rd09, plattMPC_Lbyrd09);
plattMPC_LA_DIAG_CHOL_ONELOOP_LBUB_8_8_8(params->H11, plattMPC_llbbyslb10, plattMPC_lbIdx10, plattMPC_lubbysub10, plattMPC_ubIdx10, plattMPC_Phi10);
plattMPC_LA_DIAG_MATRIXFORWARDSUB_6_8(plattMPC_Phi10, plattMPC_C00, plattMPC_V10);
plattMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_6_8(plattMPC_Phi10, plattMPC_D01, plattMPC_W10);
plattMPC_LA_DENSE_DIAGZERO_MMTM_6_8_6(plattMPC_W10, plattMPC_V10, plattMPC_Ysd11);
plattMPC_LA_DIAG_FORWARDSUB_8(plattMPC_Phi10, plattMPC_rd10, plattMPC_Lbyrd10);
plattMPC_LA_DIAG_CHOL_ONELOOP_LBUB_8_8_8(params->H12, plattMPC_llbbyslb11, plattMPC_lbIdx11, plattMPC_lubbysub11, plattMPC_ubIdx11, plattMPC_Phi11);
plattMPC_LA_DIAG_MATRIXFORWARDSUB_6_8(plattMPC_Phi11, plattMPC_C00, plattMPC_V11);
plattMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_6_8(plattMPC_Phi11, plattMPC_D01, plattMPC_W11);
plattMPC_LA_DENSE_DIAGZERO_MMTM_6_8_6(plattMPC_W11, plattMPC_V11, plattMPC_Ysd12);
plattMPC_LA_DIAG_FORWARDSUB_8(plattMPC_Phi11, plattMPC_rd11, plattMPC_Lbyrd11);
plattMPC_LA_DIAG_CHOL_ONELOOP_LBUB_8_8_8(params->H13, plattMPC_llbbyslb12, plattMPC_lbIdx12, plattMPC_lubbysub12, plattMPC_ubIdx12, plattMPC_Phi12);
plattMPC_LA_DIAG_MATRIXFORWARDSUB_6_8(plattMPC_Phi12, plattMPC_C00, plattMPC_V12);
plattMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_6_8(plattMPC_Phi12, plattMPC_D01, plattMPC_W12);
plattMPC_LA_DENSE_DIAGZERO_MMTM_6_8_6(plattMPC_W12, plattMPC_V12, plattMPC_Ysd13);
plattMPC_LA_DIAG_FORWARDSUB_8(plattMPC_Phi12, plattMPC_rd12, plattMPC_Lbyrd12);
plattMPC_LA_DIAG_CHOL_ONELOOP_LBUB_8_8_8(params->H14, plattMPC_llbbyslb13, plattMPC_lbIdx13, plattMPC_lubbysub13, plattMPC_ubIdx13, plattMPC_Phi13);
plattMPC_LA_DIAG_MATRIXFORWARDSUB_6_8(plattMPC_Phi13, plattMPC_C00, plattMPC_V13);
plattMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_6_8(plattMPC_Phi13, plattMPC_D01, plattMPC_W13);
plattMPC_LA_DENSE_DIAGZERO_MMTM_6_8_6(plattMPC_W13, plattMPC_V13, plattMPC_Ysd14);
plattMPC_LA_DIAG_FORWARDSUB_8(plattMPC_Phi13, plattMPC_rd13, plattMPC_Lbyrd13);
plattMPC_LA_DIAG_CHOL_ONELOOP_LBUB_6_6_6(params->H15, plattMPC_llbbyslb14, plattMPC_lbIdx14, plattMPC_lubbysub14, plattMPC_ubIdx14, plattMPC_Phi14);
plattMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_6_6(plattMPC_Phi14, plattMPC_D14, plattMPC_W14);
plattMPC_LA_DIAG_FORWARDSUB_6(plattMPC_Phi14, plattMPC_rd14, plattMPC_Lbyrd14);
plattMPC_LA_DIAGZERO_MMT_6(plattMPC_W00, plattMPC_Yd00);
plattMPC_LA_DIAGZERO_MVMSUB7_6(plattMPC_W00, plattMPC_Lbyrd00, plattMPC_re00, plattMPC_beta00);
plattMPC_LA_DENSE_DIAGZERO_MMT2_6_8_8(plattMPC_V00, plattMPC_W01, plattMPC_Yd01);
plattMPC_LA_DENSE_DIAGZERO_2MVMSUB2_6_8_8(plattMPC_V00, plattMPC_Lbyrd00, plattMPC_W01, plattMPC_Lbyrd01, plattMPC_re01, plattMPC_beta01);
plattMPC_LA_DENSE_DIAGZERO_MMT2_6_8_8(plattMPC_V01, plattMPC_W02, plattMPC_Yd02);
plattMPC_LA_DENSE_DIAGZERO_2MVMSUB2_6_8_8(plattMPC_V01, plattMPC_Lbyrd01, plattMPC_W02, plattMPC_Lbyrd02, plattMPC_re02, plattMPC_beta02);
plattMPC_LA_DENSE_DIAGZERO_MMT2_6_8_8(plattMPC_V02, plattMPC_W03, plattMPC_Yd03);
plattMPC_LA_DENSE_DIAGZERO_2MVMSUB2_6_8_8(plattMPC_V02, plattMPC_Lbyrd02, plattMPC_W03, plattMPC_Lbyrd03, plattMPC_re03, plattMPC_beta03);
plattMPC_LA_DENSE_DIAGZERO_MMT2_6_8_8(plattMPC_V03, plattMPC_W04, plattMPC_Yd04);
plattMPC_LA_DENSE_DIAGZERO_2MVMSUB2_6_8_8(plattMPC_V03, plattMPC_Lbyrd03, plattMPC_W04, plattMPC_Lbyrd04, plattMPC_re04, plattMPC_beta04);
plattMPC_LA_DENSE_DIAGZERO_MMT2_6_8_8(plattMPC_V04, plattMPC_W05, plattMPC_Yd05);
plattMPC_LA_DENSE_DIAGZERO_2MVMSUB2_6_8_8(plattMPC_V04, plattMPC_Lbyrd04, plattMPC_W05, plattMPC_Lbyrd05, plattMPC_re05, plattMPC_beta05);
plattMPC_LA_DENSE_DIAGZERO_MMT2_6_8_8(plattMPC_V05, plattMPC_W06, plattMPC_Yd06);
plattMPC_LA_DENSE_DIAGZERO_2MVMSUB2_6_8_8(plattMPC_V05, plattMPC_Lbyrd05, plattMPC_W06, plattMPC_Lbyrd06, plattMPC_re06, plattMPC_beta06);
plattMPC_LA_DENSE_DIAGZERO_MMT2_6_8_8(plattMPC_V06, plattMPC_W07, plattMPC_Yd07);
plattMPC_LA_DENSE_DIAGZERO_2MVMSUB2_6_8_8(plattMPC_V06, plattMPC_Lbyrd06, plattMPC_W07, plattMPC_Lbyrd07, plattMPC_re07, plattMPC_beta07);
plattMPC_LA_DENSE_DIAGZERO_MMT2_6_8_8(plattMPC_V07, plattMPC_W08, plattMPC_Yd08);
plattMPC_LA_DENSE_DIAGZERO_2MVMSUB2_6_8_8(plattMPC_V07, plattMPC_Lbyrd07, plattMPC_W08, plattMPC_Lbyrd08, plattMPC_re08, plattMPC_beta08);
plattMPC_LA_DENSE_DIAGZERO_MMT2_6_8_8(plattMPC_V08, plattMPC_W09, plattMPC_Yd09);
plattMPC_LA_DENSE_DIAGZERO_2MVMSUB2_6_8_8(plattMPC_V08, plattMPC_Lbyrd08, plattMPC_W09, plattMPC_Lbyrd09, plattMPC_re09, plattMPC_beta09);
plattMPC_LA_DENSE_DIAGZERO_MMT2_6_8_8(plattMPC_V09, plattMPC_W10, plattMPC_Yd10);
plattMPC_LA_DENSE_DIAGZERO_2MVMSUB2_6_8_8(plattMPC_V09, plattMPC_Lbyrd09, plattMPC_W10, plattMPC_Lbyrd10, plattMPC_re10, plattMPC_beta10);
plattMPC_LA_DENSE_DIAGZERO_MMT2_6_8_8(plattMPC_V10, plattMPC_W11, plattMPC_Yd11);
plattMPC_LA_DENSE_DIAGZERO_2MVMSUB2_6_8_8(plattMPC_V10, plattMPC_Lbyrd10, plattMPC_W11, plattMPC_Lbyrd11, plattMPC_re11, plattMPC_beta11);
plattMPC_LA_DENSE_DIAGZERO_MMT2_6_8_8(plattMPC_V11, plattMPC_W12, plattMPC_Yd12);
plattMPC_LA_DENSE_DIAGZERO_2MVMSUB2_6_8_8(plattMPC_V11, plattMPC_Lbyrd11, plattMPC_W12, plattMPC_Lbyrd12, plattMPC_re12, plattMPC_beta12);
plattMPC_LA_DENSE_DIAGZERO_MMT2_6_8_8(plattMPC_V12, plattMPC_W13, plattMPC_Yd13);
plattMPC_LA_DENSE_DIAGZERO_2MVMSUB2_6_8_8(plattMPC_V12, plattMPC_Lbyrd12, plattMPC_W13, plattMPC_Lbyrd13, plattMPC_re13, plattMPC_beta13);
plattMPC_LA_DENSE_DIAGZERO_MMT2_6_8_6(plattMPC_V13, plattMPC_W14, plattMPC_Yd14);
plattMPC_LA_DENSE_DIAGZERO_2MVMSUB2_6_8_6(plattMPC_V13, plattMPC_Lbyrd13, plattMPC_W14, plattMPC_Lbyrd14, plattMPC_re14, plattMPC_beta14);
plattMPC_LA_DENSE_CHOL_6(plattMPC_Yd00, plattMPC_Ld00);
plattMPC_LA_DENSE_FORWARDSUB_6(plattMPC_Ld00, plattMPC_beta00, plattMPC_yy00);
plattMPC_LA_DENSE_MATRIXTFORWARDSUB_6_6(plattMPC_Ld00, plattMPC_Ysd01, plattMPC_Lsd01);
plattMPC_LA_DENSE_MMTSUB_6_6(plattMPC_Lsd01, plattMPC_Yd01);
plattMPC_LA_DENSE_CHOL_6(plattMPC_Yd01, plattMPC_Ld01);
plattMPC_LA_DENSE_MVMSUB1_6_6(plattMPC_Lsd01, plattMPC_yy00, plattMPC_beta01, plattMPC_bmy01);
plattMPC_LA_DENSE_FORWARDSUB_6(plattMPC_Ld01, plattMPC_bmy01, plattMPC_yy01);
plattMPC_LA_DENSE_MATRIXTFORWARDSUB_6_6(plattMPC_Ld01, plattMPC_Ysd02, plattMPC_Lsd02);
plattMPC_LA_DENSE_MMTSUB_6_6(plattMPC_Lsd02, plattMPC_Yd02);
plattMPC_LA_DENSE_CHOL_6(plattMPC_Yd02, plattMPC_Ld02);
plattMPC_LA_DENSE_MVMSUB1_6_6(plattMPC_Lsd02, plattMPC_yy01, plattMPC_beta02, plattMPC_bmy02);
plattMPC_LA_DENSE_FORWARDSUB_6(plattMPC_Ld02, plattMPC_bmy02, plattMPC_yy02);
plattMPC_LA_DENSE_MATRIXTFORWARDSUB_6_6(plattMPC_Ld02, plattMPC_Ysd03, plattMPC_Lsd03);
plattMPC_LA_DENSE_MMTSUB_6_6(plattMPC_Lsd03, plattMPC_Yd03);
plattMPC_LA_DENSE_CHOL_6(plattMPC_Yd03, plattMPC_Ld03);
plattMPC_LA_DENSE_MVMSUB1_6_6(plattMPC_Lsd03, plattMPC_yy02, plattMPC_beta03, plattMPC_bmy03);
plattMPC_LA_DENSE_FORWARDSUB_6(plattMPC_Ld03, plattMPC_bmy03, plattMPC_yy03);
plattMPC_LA_DENSE_MATRIXTFORWARDSUB_6_6(plattMPC_Ld03, plattMPC_Ysd04, plattMPC_Lsd04);
plattMPC_LA_DENSE_MMTSUB_6_6(plattMPC_Lsd04, plattMPC_Yd04);
plattMPC_LA_DENSE_CHOL_6(plattMPC_Yd04, plattMPC_Ld04);
plattMPC_LA_DENSE_MVMSUB1_6_6(plattMPC_Lsd04, plattMPC_yy03, plattMPC_beta04, plattMPC_bmy04);
plattMPC_LA_DENSE_FORWARDSUB_6(plattMPC_Ld04, plattMPC_bmy04, plattMPC_yy04);
plattMPC_LA_DENSE_MATRIXTFORWARDSUB_6_6(plattMPC_Ld04, plattMPC_Ysd05, plattMPC_Lsd05);
plattMPC_LA_DENSE_MMTSUB_6_6(plattMPC_Lsd05, plattMPC_Yd05);
plattMPC_LA_DENSE_CHOL_6(plattMPC_Yd05, plattMPC_Ld05);
plattMPC_LA_DENSE_MVMSUB1_6_6(plattMPC_Lsd05, plattMPC_yy04, plattMPC_beta05, plattMPC_bmy05);
plattMPC_LA_DENSE_FORWARDSUB_6(plattMPC_Ld05, plattMPC_bmy05, plattMPC_yy05);
plattMPC_LA_DENSE_MATRIXTFORWARDSUB_6_6(plattMPC_Ld05, plattMPC_Ysd06, plattMPC_Lsd06);
plattMPC_LA_DENSE_MMTSUB_6_6(plattMPC_Lsd06, plattMPC_Yd06);
plattMPC_LA_DENSE_CHOL_6(plattMPC_Yd06, plattMPC_Ld06);
plattMPC_LA_DENSE_MVMSUB1_6_6(plattMPC_Lsd06, plattMPC_yy05, plattMPC_beta06, plattMPC_bmy06);
plattMPC_LA_DENSE_FORWARDSUB_6(plattMPC_Ld06, plattMPC_bmy06, plattMPC_yy06);
plattMPC_LA_DENSE_MATRIXTFORWARDSUB_6_6(plattMPC_Ld06, plattMPC_Ysd07, plattMPC_Lsd07);
plattMPC_LA_DENSE_MMTSUB_6_6(plattMPC_Lsd07, plattMPC_Yd07);
plattMPC_LA_DENSE_CHOL_6(plattMPC_Yd07, plattMPC_Ld07);
plattMPC_LA_DENSE_MVMSUB1_6_6(plattMPC_Lsd07, plattMPC_yy06, plattMPC_beta07, plattMPC_bmy07);
plattMPC_LA_DENSE_FORWARDSUB_6(plattMPC_Ld07, plattMPC_bmy07, plattMPC_yy07);
plattMPC_LA_DENSE_MATRIXTFORWARDSUB_6_6(plattMPC_Ld07, plattMPC_Ysd08, plattMPC_Lsd08);
plattMPC_LA_DENSE_MMTSUB_6_6(plattMPC_Lsd08, plattMPC_Yd08);
plattMPC_LA_DENSE_CHOL_6(plattMPC_Yd08, plattMPC_Ld08);
plattMPC_LA_DENSE_MVMSUB1_6_6(plattMPC_Lsd08, plattMPC_yy07, plattMPC_beta08, plattMPC_bmy08);
plattMPC_LA_DENSE_FORWARDSUB_6(plattMPC_Ld08, plattMPC_bmy08, plattMPC_yy08);
plattMPC_LA_DENSE_MATRIXTFORWARDSUB_6_6(plattMPC_Ld08, plattMPC_Ysd09, plattMPC_Lsd09);
plattMPC_LA_DENSE_MMTSUB_6_6(plattMPC_Lsd09, plattMPC_Yd09);
plattMPC_LA_DENSE_CHOL_6(plattMPC_Yd09, plattMPC_Ld09);
plattMPC_LA_DENSE_MVMSUB1_6_6(plattMPC_Lsd09, plattMPC_yy08, plattMPC_beta09, plattMPC_bmy09);
plattMPC_LA_DENSE_FORWARDSUB_6(plattMPC_Ld09, plattMPC_bmy09, plattMPC_yy09);
plattMPC_LA_DENSE_MATRIXTFORWARDSUB_6_6(plattMPC_Ld09, plattMPC_Ysd10, plattMPC_Lsd10);
plattMPC_LA_DENSE_MMTSUB_6_6(plattMPC_Lsd10, plattMPC_Yd10);
plattMPC_LA_DENSE_CHOL_6(plattMPC_Yd10, plattMPC_Ld10);
plattMPC_LA_DENSE_MVMSUB1_6_6(plattMPC_Lsd10, plattMPC_yy09, plattMPC_beta10, plattMPC_bmy10);
plattMPC_LA_DENSE_FORWARDSUB_6(plattMPC_Ld10, plattMPC_bmy10, plattMPC_yy10);
plattMPC_LA_DENSE_MATRIXTFORWARDSUB_6_6(plattMPC_Ld10, plattMPC_Ysd11, plattMPC_Lsd11);
plattMPC_LA_DENSE_MMTSUB_6_6(plattMPC_Lsd11, plattMPC_Yd11);
plattMPC_LA_DENSE_CHOL_6(plattMPC_Yd11, plattMPC_Ld11);
plattMPC_LA_DENSE_MVMSUB1_6_6(plattMPC_Lsd11, plattMPC_yy10, plattMPC_beta11, plattMPC_bmy11);
plattMPC_LA_DENSE_FORWARDSUB_6(plattMPC_Ld11, plattMPC_bmy11, plattMPC_yy11);
plattMPC_LA_DENSE_MATRIXTFORWARDSUB_6_6(plattMPC_Ld11, plattMPC_Ysd12, plattMPC_Lsd12);
plattMPC_LA_DENSE_MMTSUB_6_6(plattMPC_Lsd12, plattMPC_Yd12);
plattMPC_LA_DENSE_CHOL_6(plattMPC_Yd12, plattMPC_Ld12);
plattMPC_LA_DENSE_MVMSUB1_6_6(plattMPC_Lsd12, plattMPC_yy11, plattMPC_beta12, plattMPC_bmy12);
plattMPC_LA_DENSE_FORWARDSUB_6(plattMPC_Ld12, plattMPC_bmy12, plattMPC_yy12);
plattMPC_LA_DENSE_MATRIXTFORWARDSUB_6_6(plattMPC_Ld12, plattMPC_Ysd13, plattMPC_Lsd13);
plattMPC_LA_DENSE_MMTSUB_6_6(plattMPC_Lsd13, plattMPC_Yd13);
plattMPC_LA_DENSE_CHOL_6(plattMPC_Yd13, plattMPC_Ld13);
plattMPC_LA_DENSE_MVMSUB1_6_6(plattMPC_Lsd13, plattMPC_yy12, plattMPC_beta13, plattMPC_bmy13);
plattMPC_LA_DENSE_FORWARDSUB_6(plattMPC_Ld13, plattMPC_bmy13, plattMPC_yy13);
plattMPC_LA_DENSE_MATRIXTFORWARDSUB_6_6(plattMPC_Ld13, plattMPC_Ysd14, plattMPC_Lsd14);
plattMPC_LA_DENSE_MMTSUB_6_6(plattMPC_Lsd14, plattMPC_Yd14);
plattMPC_LA_DENSE_CHOL_6(plattMPC_Yd14, plattMPC_Ld14);
plattMPC_LA_DENSE_MVMSUB1_6_6(plattMPC_Lsd14, plattMPC_yy13, plattMPC_beta14, plattMPC_bmy14);
plattMPC_LA_DENSE_FORWARDSUB_6(plattMPC_Ld14, plattMPC_bmy14, plattMPC_yy14);
plattMPC_LA_DENSE_BACKWARDSUB_6(plattMPC_Ld14, plattMPC_yy14, plattMPC_dvaff14);
plattMPC_LA_DENSE_MTVMSUB_6_6(plattMPC_Lsd14, plattMPC_dvaff14, plattMPC_yy13, plattMPC_bmy13);
plattMPC_LA_DENSE_BACKWARDSUB_6(plattMPC_Ld13, plattMPC_bmy13, plattMPC_dvaff13);
plattMPC_LA_DENSE_MTVMSUB_6_6(plattMPC_Lsd13, plattMPC_dvaff13, plattMPC_yy12, plattMPC_bmy12);
plattMPC_LA_DENSE_BACKWARDSUB_6(plattMPC_Ld12, plattMPC_bmy12, plattMPC_dvaff12);
plattMPC_LA_DENSE_MTVMSUB_6_6(plattMPC_Lsd12, plattMPC_dvaff12, plattMPC_yy11, plattMPC_bmy11);
plattMPC_LA_DENSE_BACKWARDSUB_6(plattMPC_Ld11, plattMPC_bmy11, plattMPC_dvaff11);
plattMPC_LA_DENSE_MTVMSUB_6_6(plattMPC_Lsd11, plattMPC_dvaff11, plattMPC_yy10, plattMPC_bmy10);
plattMPC_LA_DENSE_BACKWARDSUB_6(plattMPC_Ld10, plattMPC_bmy10, plattMPC_dvaff10);
plattMPC_LA_DENSE_MTVMSUB_6_6(plattMPC_Lsd10, plattMPC_dvaff10, plattMPC_yy09, plattMPC_bmy09);
plattMPC_LA_DENSE_BACKWARDSUB_6(plattMPC_Ld09, plattMPC_bmy09, plattMPC_dvaff09);
plattMPC_LA_DENSE_MTVMSUB_6_6(plattMPC_Lsd09, plattMPC_dvaff09, plattMPC_yy08, plattMPC_bmy08);
plattMPC_LA_DENSE_BACKWARDSUB_6(plattMPC_Ld08, plattMPC_bmy08, plattMPC_dvaff08);
plattMPC_LA_DENSE_MTVMSUB_6_6(plattMPC_Lsd08, plattMPC_dvaff08, plattMPC_yy07, plattMPC_bmy07);
plattMPC_LA_DENSE_BACKWARDSUB_6(plattMPC_Ld07, plattMPC_bmy07, plattMPC_dvaff07);
plattMPC_LA_DENSE_MTVMSUB_6_6(plattMPC_Lsd07, plattMPC_dvaff07, plattMPC_yy06, plattMPC_bmy06);
plattMPC_LA_DENSE_BACKWARDSUB_6(plattMPC_Ld06, plattMPC_bmy06, plattMPC_dvaff06);
plattMPC_LA_DENSE_MTVMSUB_6_6(plattMPC_Lsd06, plattMPC_dvaff06, plattMPC_yy05, plattMPC_bmy05);
plattMPC_LA_DENSE_BACKWARDSUB_6(plattMPC_Ld05, plattMPC_bmy05, plattMPC_dvaff05);
plattMPC_LA_DENSE_MTVMSUB_6_6(plattMPC_Lsd05, plattMPC_dvaff05, plattMPC_yy04, plattMPC_bmy04);
plattMPC_LA_DENSE_BACKWARDSUB_6(plattMPC_Ld04, plattMPC_bmy04, plattMPC_dvaff04);
plattMPC_LA_DENSE_MTVMSUB_6_6(plattMPC_Lsd04, plattMPC_dvaff04, plattMPC_yy03, plattMPC_bmy03);
plattMPC_LA_DENSE_BACKWARDSUB_6(plattMPC_Ld03, plattMPC_bmy03, plattMPC_dvaff03);
plattMPC_LA_DENSE_MTVMSUB_6_6(plattMPC_Lsd03, plattMPC_dvaff03, plattMPC_yy02, plattMPC_bmy02);
plattMPC_LA_DENSE_BACKWARDSUB_6(plattMPC_Ld02, plattMPC_bmy02, plattMPC_dvaff02);
plattMPC_LA_DENSE_MTVMSUB_6_6(plattMPC_Lsd02, plattMPC_dvaff02, plattMPC_yy01, plattMPC_bmy01);
plattMPC_LA_DENSE_BACKWARDSUB_6(plattMPC_Ld01, plattMPC_bmy01, plattMPC_dvaff01);
plattMPC_LA_DENSE_MTVMSUB_6_6(plattMPC_Lsd01, plattMPC_dvaff01, plattMPC_yy00, plattMPC_bmy00);
plattMPC_LA_DENSE_BACKWARDSUB_6(plattMPC_Ld00, plattMPC_bmy00, plattMPC_dvaff00);
plattMPC_LA_DENSE_DIAGZERO_MTVM2_6_8_6(plattMPC_C00, plattMPC_dvaff01, plattMPC_D00, plattMPC_dvaff00, plattMPC_grad_eq00);
plattMPC_LA_DENSE_DIAGZERO_MTVM2_6_8_6(plattMPC_C00, plattMPC_dvaff02, plattMPC_D01, plattMPC_dvaff01, plattMPC_grad_eq01);
plattMPC_LA_DENSE_DIAGZERO_MTVM2_6_8_6(plattMPC_C00, plattMPC_dvaff03, plattMPC_D01, plattMPC_dvaff02, plattMPC_grad_eq02);
plattMPC_LA_DENSE_DIAGZERO_MTVM2_6_8_6(plattMPC_C00, plattMPC_dvaff04, plattMPC_D01, plattMPC_dvaff03, plattMPC_grad_eq03);
plattMPC_LA_DENSE_DIAGZERO_MTVM2_6_8_6(plattMPC_C00, plattMPC_dvaff05, plattMPC_D01, plattMPC_dvaff04, plattMPC_grad_eq04);
plattMPC_LA_DENSE_DIAGZERO_MTVM2_6_8_6(plattMPC_C00, plattMPC_dvaff06, plattMPC_D01, plattMPC_dvaff05, plattMPC_grad_eq05);
plattMPC_LA_DENSE_DIAGZERO_MTVM2_6_8_6(plattMPC_C00, plattMPC_dvaff07, plattMPC_D01, plattMPC_dvaff06, plattMPC_grad_eq06);
plattMPC_LA_DENSE_DIAGZERO_MTVM2_6_8_6(plattMPC_C00, plattMPC_dvaff08, plattMPC_D01, plattMPC_dvaff07, plattMPC_grad_eq07);
plattMPC_LA_DENSE_DIAGZERO_MTVM2_6_8_6(plattMPC_C00, plattMPC_dvaff09, plattMPC_D01, plattMPC_dvaff08, plattMPC_grad_eq08);
plattMPC_LA_DENSE_DIAGZERO_MTVM2_6_8_6(plattMPC_C00, plattMPC_dvaff10, plattMPC_D01, plattMPC_dvaff09, plattMPC_grad_eq09);
plattMPC_LA_DENSE_DIAGZERO_MTVM2_6_8_6(plattMPC_C00, plattMPC_dvaff11, plattMPC_D01, plattMPC_dvaff10, plattMPC_grad_eq10);
plattMPC_LA_DENSE_DIAGZERO_MTVM2_6_8_6(plattMPC_C00, plattMPC_dvaff12, plattMPC_D01, plattMPC_dvaff11, plattMPC_grad_eq11);
plattMPC_LA_DENSE_DIAGZERO_MTVM2_6_8_6(plattMPC_C00, plattMPC_dvaff13, plattMPC_D01, plattMPC_dvaff12, plattMPC_grad_eq12);
plattMPC_LA_DENSE_DIAGZERO_MTVM2_6_8_6(plattMPC_C00, plattMPC_dvaff14, plattMPC_D01, plattMPC_dvaff13, plattMPC_grad_eq13);
plattMPC_LA_DIAGZERO_MTVM_6_6(plattMPC_D14, plattMPC_dvaff14, plattMPC_grad_eq14);
plattMPC_LA_VSUB2_118(plattMPC_rd, plattMPC_grad_eq, plattMPC_rd);
plattMPC_LA_DIAG_FORWARDBACKWARDSUB_8(plattMPC_Phi00, plattMPC_rd00, plattMPC_dzaff00);
plattMPC_LA_DIAG_FORWARDBACKWARDSUB_8(plattMPC_Phi01, plattMPC_rd01, plattMPC_dzaff01);
plattMPC_LA_DIAG_FORWARDBACKWARDSUB_8(plattMPC_Phi02, plattMPC_rd02, plattMPC_dzaff02);
plattMPC_LA_DIAG_FORWARDBACKWARDSUB_8(plattMPC_Phi03, plattMPC_rd03, plattMPC_dzaff03);
plattMPC_LA_DIAG_FORWARDBACKWARDSUB_8(plattMPC_Phi04, plattMPC_rd04, plattMPC_dzaff04);
plattMPC_LA_DIAG_FORWARDBACKWARDSUB_8(plattMPC_Phi05, plattMPC_rd05, plattMPC_dzaff05);
plattMPC_LA_DIAG_FORWARDBACKWARDSUB_8(plattMPC_Phi06, plattMPC_rd06, plattMPC_dzaff06);
plattMPC_LA_DIAG_FORWARDBACKWARDSUB_8(plattMPC_Phi07, plattMPC_rd07, plattMPC_dzaff07);
plattMPC_LA_DIAG_FORWARDBACKWARDSUB_8(plattMPC_Phi08, plattMPC_rd08, plattMPC_dzaff08);
plattMPC_LA_DIAG_FORWARDBACKWARDSUB_8(plattMPC_Phi09, plattMPC_rd09, plattMPC_dzaff09);
plattMPC_LA_DIAG_FORWARDBACKWARDSUB_8(plattMPC_Phi10, plattMPC_rd10, plattMPC_dzaff10);
plattMPC_LA_DIAG_FORWARDBACKWARDSUB_8(plattMPC_Phi11, plattMPC_rd11, plattMPC_dzaff11);
plattMPC_LA_DIAG_FORWARDBACKWARDSUB_8(plattMPC_Phi12, plattMPC_rd12, plattMPC_dzaff12);
plattMPC_LA_DIAG_FORWARDBACKWARDSUB_8(plattMPC_Phi13, plattMPC_rd13, plattMPC_dzaff13);
plattMPC_LA_DIAG_FORWARDBACKWARDSUB_6(plattMPC_Phi14, plattMPC_rd14, plattMPC_dzaff14);
plattMPC_LA_VSUB_INDEXED_8(plattMPC_dzaff00, plattMPC_lbIdx00, plattMPC_rilb00, plattMPC_dslbaff00);
plattMPC_LA_VSUB3_8(plattMPC_llbbyslb00, plattMPC_dslbaff00, plattMPC_llb00, plattMPC_dllbaff00);
plattMPC_LA_VSUB2_INDEXED_8(plattMPC_riub00, plattMPC_dzaff00, plattMPC_ubIdx00, plattMPC_dsubaff00);
plattMPC_LA_VSUB3_8(plattMPC_lubbysub00, plattMPC_dsubaff00, plattMPC_lub00, plattMPC_dlubaff00);
plattMPC_LA_VSUB_INDEXED_8(plattMPC_dzaff01, plattMPC_lbIdx01, plattMPC_rilb01, plattMPC_dslbaff01);
plattMPC_LA_VSUB3_8(plattMPC_llbbyslb01, plattMPC_dslbaff01, plattMPC_llb01, plattMPC_dllbaff01);
plattMPC_LA_VSUB2_INDEXED_8(plattMPC_riub01, plattMPC_dzaff01, plattMPC_ubIdx01, plattMPC_dsubaff01);
plattMPC_LA_VSUB3_8(plattMPC_lubbysub01, plattMPC_dsubaff01, plattMPC_lub01, plattMPC_dlubaff01);
plattMPC_LA_VSUB_INDEXED_8(plattMPC_dzaff02, plattMPC_lbIdx02, plattMPC_rilb02, plattMPC_dslbaff02);
plattMPC_LA_VSUB3_8(plattMPC_llbbyslb02, plattMPC_dslbaff02, plattMPC_llb02, plattMPC_dllbaff02);
plattMPC_LA_VSUB2_INDEXED_8(plattMPC_riub02, plattMPC_dzaff02, plattMPC_ubIdx02, plattMPC_dsubaff02);
plattMPC_LA_VSUB3_8(plattMPC_lubbysub02, plattMPC_dsubaff02, plattMPC_lub02, plattMPC_dlubaff02);
plattMPC_LA_VSUB_INDEXED_8(plattMPC_dzaff03, plattMPC_lbIdx03, plattMPC_rilb03, plattMPC_dslbaff03);
plattMPC_LA_VSUB3_8(plattMPC_llbbyslb03, plattMPC_dslbaff03, plattMPC_llb03, plattMPC_dllbaff03);
plattMPC_LA_VSUB2_INDEXED_8(plattMPC_riub03, plattMPC_dzaff03, plattMPC_ubIdx03, plattMPC_dsubaff03);
plattMPC_LA_VSUB3_8(plattMPC_lubbysub03, plattMPC_dsubaff03, plattMPC_lub03, plattMPC_dlubaff03);
plattMPC_LA_VSUB_INDEXED_8(plattMPC_dzaff04, plattMPC_lbIdx04, plattMPC_rilb04, plattMPC_dslbaff04);
plattMPC_LA_VSUB3_8(plattMPC_llbbyslb04, plattMPC_dslbaff04, plattMPC_llb04, plattMPC_dllbaff04);
plattMPC_LA_VSUB2_INDEXED_8(plattMPC_riub04, plattMPC_dzaff04, plattMPC_ubIdx04, plattMPC_dsubaff04);
plattMPC_LA_VSUB3_8(plattMPC_lubbysub04, plattMPC_dsubaff04, plattMPC_lub04, plattMPC_dlubaff04);
plattMPC_LA_VSUB_INDEXED_8(plattMPC_dzaff05, plattMPC_lbIdx05, plattMPC_rilb05, plattMPC_dslbaff05);
plattMPC_LA_VSUB3_8(plattMPC_llbbyslb05, plattMPC_dslbaff05, plattMPC_llb05, plattMPC_dllbaff05);
plattMPC_LA_VSUB2_INDEXED_8(plattMPC_riub05, plattMPC_dzaff05, plattMPC_ubIdx05, plattMPC_dsubaff05);
plattMPC_LA_VSUB3_8(plattMPC_lubbysub05, plattMPC_dsubaff05, plattMPC_lub05, plattMPC_dlubaff05);
plattMPC_LA_VSUB_INDEXED_8(plattMPC_dzaff06, plattMPC_lbIdx06, plattMPC_rilb06, plattMPC_dslbaff06);
plattMPC_LA_VSUB3_8(plattMPC_llbbyslb06, plattMPC_dslbaff06, plattMPC_llb06, plattMPC_dllbaff06);
plattMPC_LA_VSUB2_INDEXED_8(plattMPC_riub06, plattMPC_dzaff06, plattMPC_ubIdx06, plattMPC_dsubaff06);
plattMPC_LA_VSUB3_8(plattMPC_lubbysub06, plattMPC_dsubaff06, plattMPC_lub06, plattMPC_dlubaff06);
plattMPC_LA_VSUB_INDEXED_8(plattMPC_dzaff07, plattMPC_lbIdx07, plattMPC_rilb07, plattMPC_dslbaff07);
plattMPC_LA_VSUB3_8(plattMPC_llbbyslb07, plattMPC_dslbaff07, plattMPC_llb07, plattMPC_dllbaff07);
plattMPC_LA_VSUB2_INDEXED_8(plattMPC_riub07, plattMPC_dzaff07, plattMPC_ubIdx07, plattMPC_dsubaff07);
plattMPC_LA_VSUB3_8(plattMPC_lubbysub07, plattMPC_dsubaff07, plattMPC_lub07, plattMPC_dlubaff07);
plattMPC_LA_VSUB_INDEXED_8(plattMPC_dzaff08, plattMPC_lbIdx08, plattMPC_rilb08, plattMPC_dslbaff08);
plattMPC_LA_VSUB3_8(plattMPC_llbbyslb08, plattMPC_dslbaff08, plattMPC_llb08, plattMPC_dllbaff08);
plattMPC_LA_VSUB2_INDEXED_8(plattMPC_riub08, plattMPC_dzaff08, plattMPC_ubIdx08, plattMPC_dsubaff08);
plattMPC_LA_VSUB3_8(plattMPC_lubbysub08, plattMPC_dsubaff08, plattMPC_lub08, plattMPC_dlubaff08);
plattMPC_LA_VSUB_INDEXED_8(plattMPC_dzaff09, plattMPC_lbIdx09, plattMPC_rilb09, plattMPC_dslbaff09);
plattMPC_LA_VSUB3_8(plattMPC_llbbyslb09, plattMPC_dslbaff09, plattMPC_llb09, plattMPC_dllbaff09);
plattMPC_LA_VSUB2_INDEXED_8(plattMPC_riub09, plattMPC_dzaff09, plattMPC_ubIdx09, plattMPC_dsubaff09);
plattMPC_LA_VSUB3_8(plattMPC_lubbysub09, plattMPC_dsubaff09, plattMPC_lub09, plattMPC_dlubaff09);
plattMPC_LA_VSUB_INDEXED_8(plattMPC_dzaff10, plattMPC_lbIdx10, plattMPC_rilb10, plattMPC_dslbaff10);
plattMPC_LA_VSUB3_8(plattMPC_llbbyslb10, plattMPC_dslbaff10, plattMPC_llb10, plattMPC_dllbaff10);
plattMPC_LA_VSUB2_INDEXED_8(plattMPC_riub10, plattMPC_dzaff10, plattMPC_ubIdx10, plattMPC_dsubaff10);
plattMPC_LA_VSUB3_8(plattMPC_lubbysub10, plattMPC_dsubaff10, plattMPC_lub10, plattMPC_dlubaff10);
plattMPC_LA_VSUB_INDEXED_8(plattMPC_dzaff11, plattMPC_lbIdx11, plattMPC_rilb11, plattMPC_dslbaff11);
plattMPC_LA_VSUB3_8(plattMPC_llbbyslb11, plattMPC_dslbaff11, plattMPC_llb11, plattMPC_dllbaff11);
plattMPC_LA_VSUB2_INDEXED_8(plattMPC_riub11, plattMPC_dzaff11, plattMPC_ubIdx11, plattMPC_dsubaff11);
plattMPC_LA_VSUB3_8(plattMPC_lubbysub11, plattMPC_dsubaff11, plattMPC_lub11, plattMPC_dlubaff11);
plattMPC_LA_VSUB_INDEXED_8(plattMPC_dzaff12, plattMPC_lbIdx12, plattMPC_rilb12, plattMPC_dslbaff12);
plattMPC_LA_VSUB3_8(plattMPC_llbbyslb12, plattMPC_dslbaff12, plattMPC_llb12, plattMPC_dllbaff12);
plattMPC_LA_VSUB2_INDEXED_8(plattMPC_riub12, plattMPC_dzaff12, plattMPC_ubIdx12, plattMPC_dsubaff12);
plattMPC_LA_VSUB3_8(plattMPC_lubbysub12, plattMPC_dsubaff12, plattMPC_lub12, plattMPC_dlubaff12);
plattMPC_LA_VSUB_INDEXED_8(plattMPC_dzaff13, plattMPC_lbIdx13, plattMPC_rilb13, plattMPC_dslbaff13);
plattMPC_LA_VSUB3_8(plattMPC_llbbyslb13, plattMPC_dslbaff13, plattMPC_llb13, plattMPC_dllbaff13);
plattMPC_LA_VSUB2_INDEXED_8(plattMPC_riub13, plattMPC_dzaff13, plattMPC_ubIdx13, plattMPC_dsubaff13);
plattMPC_LA_VSUB3_8(plattMPC_lubbysub13, plattMPC_dsubaff13, plattMPC_lub13, plattMPC_dlubaff13);
plattMPC_LA_VSUB_INDEXED_6(plattMPC_dzaff14, plattMPC_lbIdx14, plattMPC_rilb14, plattMPC_dslbaff14);
plattMPC_LA_VSUB3_6(plattMPC_llbbyslb14, plattMPC_dslbaff14, plattMPC_llb14, plattMPC_dllbaff14);
plattMPC_LA_VSUB2_INDEXED_6(plattMPC_riub14, plattMPC_dzaff14, plattMPC_ubIdx14, plattMPC_dsubaff14);
plattMPC_LA_VSUB3_6(plattMPC_lubbysub14, plattMPC_dsubaff14, plattMPC_lub14, plattMPC_dlubaff14);
info->lsit_aff = plattMPC_LINESEARCH_BACKTRACKING_AFFINE(plattMPC_l, plattMPC_s, plattMPC_dl_aff, plattMPC_ds_aff, &info->step_aff, &info->mu_aff);
if( info->lsit_aff == plattMPC_NOPROGRESS ){
exitcode = plattMPC_NOPROGRESS; break;
}
sigma_3rdroot = info->mu_aff / info->mu;
info->sigma = sigma_3rdroot*sigma_3rdroot*sigma_3rdroot;
musigma = info->mu * info->sigma;
plattMPC_LA_VSUB5_236(plattMPC_ds_aff, plattMPC_dl_aff, musigma, plattMPC_ccrhs);
plattMPC_LA_VSUB6_INDEXED_8_8_8(plattMPC_ccrhsub00, plattMPC_sub00, plattMPC_ubIdx00, plattMPC_ccrhsl00, plattMPC_slb00, plattMPC_lbIdx00, plattMPC_rd00);
plattMPC_LA_VSUB6_INDEXED_8_8_8(plattMPC_ccrhsub01, plattMPC_sub01, plattMPC_ubIdx01, plattMPC_ccrhsl01, plattMPC_slb01, plattMPC_lbIdx01, plattMPC_rd01);
plattMPC_LA_DIAG_FORWARDSUB_8(plattMPC_Phi00, plattMPC_rd00, plattMPC_Lbyrd00);
plattMPC_LA_DIAG_FORWARDSUB_8(plattMPC_Phi01, plattMPC_rd01, plattMPC_Lbyrd01);
plattMPC_LA_DIAGZERO_MVM_6(plattMPC_W00, plattMPC_Lbyrd00, plattMPC_beta00);
plattMPC_LA_DENSE_FORWARDSUB_6(plattMPC_Ld00, plattMPC_beta00, plattMPC_yy00);
plattMPC_LA_DENSE_DIAGZERO_2MVMADD_6_8_8(plattMPC_V00, plattMPC_Lbyrd00, plattMPC_W01, plattMPC_Lbyrd01, plattMPC_beta01);
plattMPC_LA_DENSE_MVMSUB1_6_6(plattMPC_Lsd01, plattMPC_yy00, plattMPC_beta01, plattMPC_bmy01);
plattMPC_LA_DENSE_FORWARDSUB_6(plattMPC_Ld01, plattMPC_bmy01, plattMPC_yy01);
plattMPC_LA_VSUB6_INDEXED_8_8_8(plattMPC_ccrhsub02, plattMPC_sub02, plattMPC_ubIdx02, plattMPC_ccrhsl02, plattMPC_slb02, plattMPC_lbIdx02, plattMPC_rd02);
plattMPC_LA_DIAG_FORWARDSUB_8(plattMPC_Phi02, plattMPC_rd02, plattMPC_Lbyrd02);
plattMPC_LA_DENSE_DIAGZERO_2MVMADD_6_8_8(plattMPC_V01, plattMPC_Lbyrd01, plattMPC_W02, plattMPC_Lbyrd02, plattMPC_beta02);
plattMPC_LA_DENSE_MVMSUB1_6_6(plattMPC_Lsd02, plattMPC_yy01, plattMPC_beta02, plattMPC_bmy02);
plattMPC_LA_DENSE_FORWARDSUB_6(plattMPC_Ld02, plattMPC_bmy02, plattMPC_yy02);
plattMPC_LA_VSUB6_INDEXED_8_8_8(plattMPC_ccrhsub03, plattMPC_sub03, plattMPC_ubIdx03, plattMPC_ccrhsl03, plattMPC_slb03, plattMPC_lbIdx03, plattMPC_rd03);
plattMPC_LA_DIAG_FORWARDSUB_8(plattMPC_Phi03, plattMPC_rd03, plattMPC_Lbyrd03);
plattMPC_LA_DENSE_DIAGZERO_2MVMADD_6_8_8(plattMPC_V02, plattMPC_Lbyrd02, plattMPC_W03, plattMPC_Lbyrd03, plattMPC_beta03);
plattMPC_LA_DENSE_MVMSUB1_6_6(plattMPC_Lsd03, plattMPC_yy02, plattMPC_beta03, plattMPC_bmy03);
plattMPC_LA_DENSE_FORWARDSUB_6(plattMPC_Ld03, plattMPC_bmy03, plattMPC_yy03);
plattMPC_LA_VSUB6_INDEXED_8_8_8(plattMPC_ccrhsub04, plattMPC_sub04, plattMPC_ubIdx04, plattMPC_ccrhsl04, plattMPC_slb04, plattMPC_lbIdx04, plattMPC_rd04);
plattMPC_LA_DIAG_FORWARDSUB_8(plattMPC_Phi04, plattMPC_rd04, plattMPC_Lbyrd04);
plattMPC_LA_DENSE_DIAGZERO_2MVMADD_6_8_8(plattMPC_V03, plattMPC_Lbyrd03, plattMPC_W04, plattMPC_Lbyrd04, plattMPC_beta04);
plattMPC_LA_DENSE_MVMSUB1_6_6(plattMPC_Lsd04, plattMPC_yy03, plattMPC_beta04, plattMPC_bmy04);
plattMPC_LA_DENSE_FORWARDSUB_6(plattMPC_Ld04, plattMPC_bmy04, plattMPC_yy04);
plattMPC_LA_VSUB6_INDEXED_8_8_8(plattMPC_ccrhsub05, plattMPC_sub05, plattMPC_ubIdx05, plattMPC_ccrhsl05, plattMPC_slb05, plattMPC_lbIdx05, plattMPC_rd05);
plattMPC_LA_DIAG_FORWARDSUB_8(plattMPC_Phi05, plattMPC_rd05, plattMPC_Lbyrd05);
plattMPC_LA_DENSE_DIAGZERO_2MVMADD_6_8_8(plattMPC_V04, plattMPC_Lbyrd04, plattMPC_W05, plattMPC_Lbyrd05, plattMPC_beta05);
plattMPC_LA_DENSE_MVMSUB1_6_6(plattMPC_Lsd05, plattMPC_yy04, plattMPC_beta05, plattMPC_bmy05);
plattMPC_LA_DENSE_FORWARDSUB_6(plattMPC_Ld05, plattMPC_bmy05, plattMPC_yy05);
plattMPC_LA_VSUB6_INDEXED_8_8_8(plattMPC_ccrhsub06, plattMPC_sub06, plattMPC_ubIdx06, plattMPC_ccrhsl06, plattMPC_slb06, plattMPC_lbIdx06, plattMPC_rd06);
plattMPC_LA_DIAG_FORWARDSUB_8(plattMPC_Phi06, plattMPC_rd06, plattMPC_Lbyrd06);
plattMPC_LA_DENSE_DIAGZERO_2MVMADD_6_8_8(plattMPC_V05, plattMPC_Lbyrd05, plattMPC_W06, plattMPC_Lbyrd06, plattMPC_beta06);
plattMPC_LA_DENSE_MVMSUB1_6_6(plattMPC_Lsd06, plattMPC_yy05, plattMPC_beta06, plattMPC_bmy06);
plattMPC_LA_DENSE_FORWARDSUB_6(plattMPC_Ld06, plattMPC_bmy06, plattMPC_yy06);
plattMPC_LA_VSUB6_INDEXED_8_8_8(plattMPC_ccrhsub07, plattMPC_sub07, plattMPC_ubIdx07, plattMPC_ccrhsl07, plattMPC_slb07, plattMPC_lbIdx07, plattMPC_rd07);
plattMPC_LA_DIAG_FORWARDSUB_8(plattMPC_Phi07, plattMPC_rd07, plattMPC_Lbyrd07);
plattMPC_LA_DENSE_DIAGZERO_2MVMADD_6_8_8(plattMPC_V06, plattMPC_Lbyrd06, plattMPC_W07, plattMPC_Lbyrd07, plattMPC_beta07);
plattMPC_LA_DENSE_MVMSUB1_6_6(plattMPC_Lsd07, plattMPC_yy06, plattMPC_beta07, plattMPC_bmy07);
plattMPC_LA_DENSE_FORWARDSUB_6(plattMPC_Ld07, plattMPC_bmy07, plattMPC_yy07);
plattMPC_LA_VSUB6_INDEXED_8_8_8(plattMPC_ccrhsub08, plattMPC_sub08, plattMPC_ubIdx08, plattMPC_ccrhsl08, plattMPC_slb08, plattMPC_lbIdx08, plattMPC_rd08);
plattMPC_LA_DIAG_FORWARDSUB_8(plattMPC_Phi08, plattMPC_rd08, plattMPC_Lbyrd08);
plattMPC_LA_DENSE_DIAGZERO_2MVMADD_6_8_8(plattMPC_V07, plattMPC_Lbyrd07, plattMPC_W08, plattMPC_Lbyrd08, plattMPC_beta08);
plattMPC_LA_DENSE_MVMSUB1_6_6(plattMPC_Lsd08, plattMPC_yy07, plattMPC_beta08, plattMPC_bmy08);
plattMPC_LA_DENSE_FORWARDSUB_6(plattMPC_Ld08, plattMPC_bmy08, plattMPC_yy08);
plattMPC_LA_VSUB6_INDEXED_8_8_8(plattMPC_ccrhsub09, plattMPC_sub09, plattMPC_ubIdx09, plattMPC_ccrhsl09, plattMPC_slb09, plattMPC_lbIdx09, plattMPC_rd09);
plattMPC_LA_DIAG_FORWARDSUB_8(plattMPC_Phi09, plattMPC_rd09, plattMPC_Lbyrd09);
plattMPC_LA_DENSE_DIAGZERO_2MVMADD_6_8_8(plattMPC_V08, plattMPC_Lbyrd08, plattMPC_W09, plattMPC_Lbyrd09, plattMPC_beta09);
plattMPC_LA_DENSE_MVMSUB1_6_6(plattMPC_Lsd09, plattMPC_yy08, plattMPC_beta09, plattMPC_bmy09);
plattMPC_LA_DENSE_FORWARDSUB_6(plattMPC_Ld09, plattMPC_bmy09, plattMPC_yy09);
plattMPC_LA_VSUB6_INDEXED_8_8_8(plattMPC_ccrhsub10, plattMPC_sub10, plattMPC_ubIdx10, plattMPC_ccrhsl10, plattMPC_slb10, plattMPC_lbIdx10, plattMPC_rd10);
plattMPC_LA_DIAG_FORWARDSUB_8(plattMPC_Phi10, plattMPC_rd10, plattMPC_Lbyrd10);
plattMPC_LA_DENSE_DIAGZERO_2MVMADD_6_8_8(plattMPC_V09, plattMPC_Lbyrd09, plattMPC_W10, plattMPC_Lbyrd10, plattMPC_beta10);
plattMPC_LA_DENSE_MVMSUB1_6_6(plattMPC_Lsd10, plattMPC_yy09, plattMPC_beta10, plattMPC_bmy10);
plattMPC_LA_DENSE_FORWARDSUB_6(plattMPC_Ld10, plattMPC_bmy10, plattMPC_yy10);
plattMPC_LA_VSUB6_INDEXED_8_8_8(plattMPC_ccrhsub11, plattMPC_sub11, plattMPC_ubIdx11, plattMPC_ccrhsl11, plattMPC_slb11, plattMPC_lbIdx11, plattMPC_rd11);
plattMPC_LA_DIAG_FORWARDSUB_8(plattMPC_Phi11, plattMPC_rd11, plattMPC_Lbyrd11);
plattMPC_LA_DENSE_DIAGZERO_2MVMADD_6_8_8(plattMPC_V10, plattMPC_Lbyrd10, plattMPC_W11, plattMPC_Lbyrd11, plattMPC_beta11);
plattMPC_LA_DENSE_MVMSUB1_6_6(plattMPC_Lsd11, plattMPC_yy10, plattMPC_beta11, plattMPC_bmy11);
plattMPC_LA_DENSE_FORWARDSUB_6(plattMPC_Ld11, plattMPC_bmy11, plattMPC_yy11);
plattMPC_LA_VSUB6_INDEXED_8_8_8(plattMPC_ccrhsub12, plattMPC_sub12, plattMPC_ubIdx12, plattMPC_ccrhsl12, plattMPC_slb12, plattMPC_lbIdx12, plattMPC_rd12);
plattMPC_LA_DIAG_FORWARDSUB_8(plattMPC_Phi12, plattMPC_rd12, plattMPC_Lbyrd12);
plattMPC_LA_DENSE_DIAGZERO_2MVMADD_6_8_8(plattMPC_V11, plattMPC_Lbyrd11, plattMPC_W12, plattMPC_Lbyrd12, plattMPC_beta12);
plattMPC_LA_DENSE_MVMSUB1_6_6(plattMPC_Lsd12, plattMPC_yy11, plattMPC_beta12, plattMPC_bmy12);
plattMPC_LA_DENSE_FORWARDSUB_6(plattMPC_Ld12, plattMPC_bmy12, plattMPC_yy12);
plattMPC_LA_VSUB6_INDEXED_8_8_8(plattMPC_ccrhsub13, plattMPC_sub13, plattMPC_ubIdx13, plattMPC_ccrhsl13, plattMPC_slb13, plattMPC_lbIdx13, plattMPC_rd13);
plattMPC_LA_DIAG_FORWARDSUB_8(plattMPC_Phi13, plattMPC_rd13, plattMPC_Lbyrd13);
plattMPC_LA_DENSE_DIAGZERO_2MVMADD_6_8_8(plattMPC_V12, plattMPC_Lbyrd12, plattMPC_W13, plattMPC_Lbyrd13, plattMPC_beta13);
plattMPC_LA_DENSE_MVMSUB1_6_6(plattMPC_Lsd13, plattMPC_yy12, plattMPC_beta13, plattMPC_bmy13);
plattMPC_LA_DENSE_FORWARDSUB_6(plattMPC_Ld13, plattMPC_bmy13, plattMPC_yy13);
plattMPC_LA_VSUB6_INDEXED_6_6_6(plattMPC_ccrhsub14, plattMPC_sub14, plattMPC_ubIdx14, plattMPC_ccrhsl14, plattMPC_slb14, plattMPC_lbIdx14, plattMPC_rd14);
plattMPC_LA_DIAG_FORWARDSUB_6(plattMPC_Phi14, plattMPC_rd14, plattMPC_Lbyrd14);
plattMPC_LA_DENSE_DIAGZERO_2MVMADD_6_8_6(plattMPC_V13, plattMPC_Lbyrd13, plattMPC_W14, plattMPC_Lbyrd14, plattMPC_beta14);
plattMPC_LA_DENSE_MVMSUB1_6_6(plattMPC_Lsd14, plattMPC_yy13, plattMPC_beta14, plattMPC_bmy14);
plattMPC_LA_DENSE_FORWARDSUB_6(plattMPC_Ld14, plattMPC_bmy14, plattMPC_yy14);
plattMPC_LA_DENSE_BACKWARDSUB_6(plattMPC_Ld14, plattMPC_yy14, plattMPC_dvcc14);
plattMPC_LA_DENSE_MTVMSUB_6_6(plattMPC_Lsd14, plattMPC_dvcc14, plattMPC_yy13, plattMPC_bmy13);
plattMPC_LA_DENSE_BACKWARDSUB_6(plattMPC_Ld13, plattMPC_bmy13, plattMPC_dvcc13);
plattMPC_LA_DENSE_MTVMSUB_6_6(plattMPC_Lsd13, plattMPC_dvcc13, plattMPC_yy12, plattMPC_bmy12);
plattMPC_LA_DENSE_BACKWARDSUB_6(plattMPC_Ld12, plattMPC_bmy12, plattMPC_dvcc12);
plattMPC_LA_DENSE_MTVMSUB_6_6(plattMPC_Lsd12, plattMPC_dvcc12, plattMPC_yy11, plattMPC_bmy11);
plattMPC_LA_DENSE_BACKWARDSUB_6(plattMPC_Ld11, plattMPC_bmy11, plattMPC_dvcc11);
plattMPC_LA_DENSE_MTVMSUB_6_6(plattMPC_Lsd11, plattMPC_dvcc11, plattMPC_yy10, plattMPC_bmy10);
plattMPC_LA_DENSE_BACKWARDSUB_6(plattMPC_Ld10, plattMPC_bmy10, plattMPC_dvcc10);
plattMPC_LA_DENSE_MTVMSUB_6_6(plattMPC_Lsd10, plattMPC_dvcc10, plattMPC_yy09, plattMPC_bmy09);
plattMPC_LA_DENSE_BACKWARDSUB_6(plattMPC_Ld09, plattMPC_bmy09, plattMPC_dvcc09);
plattMPC_LA_DENSE_MTVMSUB_6_6(plattMPC_Lsd09, plattMPC_dvcc09, plattMPC_yy08, plattMPC_bmy08);
plattMPC_LA_DENSE_BACKWARDSUB_6(plattMPC_Ld08, plattMPC_bmy08, plattMPC_dvcc08);
plattMPC_LA_DENSE_MTVMSUB_6_6(plattMPC_Lsd08, plattMPC_dvcc08, plattMPC_yy07, plattMPC_bmy07);
plattMPC_LA_DENSE_BACKWARDSUB_6(plattMPC_Ld07, plattMPC_bmy07, plattMPC_dvcc07);
plattMPC_LA_DENSE_MTVMSUB_6_6(plattMPC_Lsd07, plattMPC_dvcc07, plattMPC_yy06, plattMPC_bmy06);
plattMPC_LA_DENSE_BACKWARDSUB_6(plattMPC_Ld06, plattMPC_bmy06, plattMPC_dvcc06);
plattMPC_LA_DENSE_MTVMSUB_6_6(plattMPC_Lsd06, plattMPC_dvcc06, plattMPC_yy05, plattMPC_bmy05);
plattMPC_LA_DENSE_BACKWARDSUB_6(plattMPC_Ld05, plattMPC_bmy05, plattMPC_dvcc05);
plattMPC_LA_DENSE_MTVMSUB_6_6(plattMPC_Lsd05, plattMPC_dvcc05, plattMPC_yy04, plattMPC_bmy04);
plattMPC_LA_DENSE_BACKWARDSUB_6(plattMPC_Ld04, plattMPC_bmy04, plattMPC_dvcc04);
plattMPC_LA_DENSE_MTVMSUB_6_6(plattMPC_Lsd04, plattMPC_dvcc04, plattMPC_yy03, plattMPC_bmy03);
plattMPC_LA_DENSE_BACKWARDSUB_6(plattMPC_Ld03, plattMPC_bmy03, plattMPC_dvcc03);
plattMPC_LA_DENSE_MTVMSUB_6_6(plattMPC_Lsd03, plattMPC_dvcc03, plattMPC_yy02, plattMPC_bmy02);
plattMPC_LA_DENSE_BACKWARDSUB_6(plattMPC_Ld02, plattMPC_bmy02, plattMPC_dvcc02);
plattMPC_LA_DENSE_MTVMSUB_6_6(plattMPC_Lsd02, plattMPC_dvcc02, plattMPC_yy01, plattMPC_bmy01);
plattMPC_LA_DENSE_BACKWARDSUB_6(plattMPC_Ld01, plattMPC_bmy01, plattMPC_dvcc01);
plattMPC_LA_DENSE_MTVMSUB_6_6(plattMPC_Lsd01, plattMPC_dvcc01, plattMPC_yy00, plattMPC_bmy00);
plattMPC_LA_DENSE_BACKWARDSUB_6(plattMPC_Ld00, plattMPC_bmy00, plattMPC_dvcc00);
plattMPC_LA_DENSE_DIAGZERO_MTVM2_6_8_6(plattMPC_C00, plattMPC_dvcc01, plattMPC_D00, plattMPC_dvcc00, plattMPC_grad_eq00);
plattMPC_LA_DENSE_DIAGZERO_MTVM2_6_8_6(plattMPC_C00, plattMPC_dvcc02, plattMPC_D01, plattMPC_dvcc01, plattMPC_grad_eq01);
plattMPC_LA_DENSE_DIAGZERO_MTVM2_6_8_6(plattMPC_C00, plattMPC_dvcc03, plattMPC_D01, plattMPC_dvcc02, plattMPC_grad_eq02);
plattMPC_LA_DENSE_DIAGZERO_MTVM2_6_8_6(plattMPC_C00, plattMPC_dvcc04, plattMPC_D01, plattMPC_dvcc03, plattMPC_grad_eq03);
plattMPC_LA_DENSE_DIAGZERO_MTVM2_6_8_6(plattMPC_C00, plattMPC_dvcc05, plattMPC_D01, plattMPC_dvcc04, plattMPC_grad_eq04);
plattMPC_LA_DENSE_DIAGZERO_MTVM2_6_8_6(plattMPC_C00, plattMPC_dvcc06, plattMPC_D01, plattMPC_dvcc05, plattMPC_grad_eq05);
plattMPC_LA_DENSE_DIAGZERO_MTVM2_6_8_6(plattMPC_C00, plattMPC_dvcc07, plattMPC_D01, plattMPC_dvcc06, plattMPC_grad_eq06);
plattMPC_LA_DENSE_DIAGZERO_MTVM2_6_8_6(plattMPC_C00, plattMPC_dvcc08, plattMPC_D01, plattMPC_dvcc07, plattMPC_grad_eq07);
plattMPC_LA_DENSE_DIAGZERO_MTVM2_6_8_6(plattMPC_C00, plattMPC_dvcc09, plattMPC_D01, plattMPC_dvcc08, plattMPC_grad_eq08);
plattMPC_LA_DENSE_DIAGZERO_MTVM2_6_8_6(plattMPC_C00, plattMPC_dvcc10, plattMPC_D01, plattMPC_dvcc09, plattMPC_grad_eq09);
plattMPC_LA_DENSE_DIAGZERO_MTVM2_6_8_6(plattMPC_C00, plattMPC_dvcc11, plattMPC_D01, plattMPC_dvcc10, plattMPC_grad_eq10);
plattMPC_LA_DENSE_DIAGZERO_MTVM2_6_8_6(plattMPC_C00, plattMPC_dvcc12, plattMPC_D01, plattMPC_dvcc11, plattMPC_grad_eq11);
plattMPC_LA_DENSE_DIAGZERO_MTVM2_6_8_6(plattMPC_C00, plattMPC_dvcc13, plattMPC_D01, plattMPC_dvcc12, plattMPC_grad_eq12);
plattMPC_LA_DENSE_DIAGZERO_MTVM2_6_8_6(plattMPC_C00, plattMPC_dvcc14, plattMPC_D01, plattMPC_dvcc13, plattMPC_grad_eq13);
plattMPC_LA_DIAGZERO_MTVM_6_6(plattMPC_D14, plattMPC_dvcc14, plattMPC_grad_eq14);
plattMPC_LA_VSUB_118(plattMPC_rd, plattMPC_grad_eq, plattMPC_rd);
plattMPC_LA_DIAG_FORWARDBACKWARDSUB_8(plattMPC_Phi00, plattMPC_rd00, plattMPC_dzcc00);
plattMPC_LA_DIAG_FORWARDBACKWARDSUB_8(plattMPC_Phi01, plattMPC_rd01, plattMPC_dzcc01);
plattMPC_LA_DIAG_FORWARDBACKWARDSUB_8(plattMPC_Phi02, plattMPC_rd02, plattMPC_dzcc02);
plattMPC_LA_DIAG_FORWARDBACKWARDSUB_8(plattMPC_Phi03, plattMPC_rd03, plattMPC_dzcc03);
plattMPC_LA_DIAG_FORWARDBACKWARDSUB_8(plattMPC_Phi04, plattMPC_rd04, plattMPC_dzcc04);
plattMPC_LA_DIAG_FORWARDBACKWARDSUB_8(plattMPC_Phi05, plattMPC_rd05, plattMPC_dzcc05);
plattMPC_LA_DIAG_FORWARDBACKWARDSUB_8(plattMPC_Phi06, plattMPC_rd06, plattMPC_dzcc06);
plattMPC_LA_DIAG_FORWARDBACKWARDSUB_8(plattMPC_Phi07, plattMPC_rd07, plattMPC_dzcc07);
plattMPC_LA_DIAG_FORWARDBACKWARDSUB_8(plattMPC_Phi08, plattMPC_rd08, plattMPC_dzcc08);
plattMPC_LA_DIAG_FORWARDBACKWARDSUB_8(plattMPC_Phi09, plattMPC_rd09, plattMPC_dzcc09);
plattMPC_LA_DIAG_FORWARDBACKWARDSUB_8(plattMPC_Phi10, plattMPC_rd10, plattMPC_dzcc10);
plattMPC_LA_DIAG_FORWARDBACKWARDSUB_8(plattMPC_Phi11, plattMPC_rd11, plattMPC_dzcc11);
plattMPC_LA_DIAG_FORWARDBACKWARDSUB_8(plattMPC_Phi12, plattMPC_rd12, plattMPC_dzcc12);
plattMPC_LA_DIAG_FORWARDBACKWARDSUB_8(plattMPC_Phi13, plattMPC_rd13, plattMPC_dzcc13);
plattMPC_LA_DIAG_FORWARDBACKWARDSUB_6(plattMPC_Phi14, plattMPC_rd14, plattMPC_dzcc14);
plattMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_8(plattMPC_ccrhsl00, plattMPC_slb00, plattMPC_llbbyslb00, plattMPC_dzcc00, plattMPC_lbIdx00, plattMPC_dllbcc00);
plattMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_8(plattMPC_ccrhsub00, plattMPC_sub00, plattMPC_lubbysub00, plattMPC_dzcc00, plattMPC_ubIdx00, plattMPC_dlubcc00);
plattMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_8(plattMPC_ccrhsl01, plattMPC_slb01, plattMPC_llbbyslb01, plattMPC_dzcc01, plattMPC_lbIdx01, plattMPC_dllbcc01);
plattMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_8(plattMPC_ccrhsub01, plattMPC_sub01, plattMPC_lubbysub01, plattMPC_dzcc01, plattMPC_ubIdx01, plattMPC_dlubcc01);
plattMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_8(plattMPC_ccrhsl02, plattMPC_slb02, plattMPC_llbbyslb02, plattMPC_dzcc02, plattMPC_lbIdx02, plattMPC_dllbcc02);
plattMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_8(plattMPC_ccrhsub02, plattMPC_sub02, plattMPC_lubbysub02, plattMPC_dzcc02, plattMPC_ubIdx02, plattMPC_dlubcc02);
plattMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_8(plattMPC_ccrhsl03, plattMPC_slb03, plattMPC_llbbyslb03, plattMPC_dzcc03, plattMPC_lbIdx03, plattMPC_dllbcc03);
plattMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_8(plattMPC_ccrhsub03, plattMPC_sub03, plattMPC_lubbysub03, plattMPC_dzcc03, plattMPC_ubIdx03, plattMPC_dlubcc03);
plattMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_8(plattMPC_ccrhsl04, plattMPC_slb04, plattMPC_llbbyslb04, plattMPC_dzcc04, plattMPC_lbIdx04, plattMPC_dllbcc04);
plattMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_8(plattMPC_ccrhsub04, plattMPC_sub04, plattMPC_lubbysub04, plattMPC_dzcc04, plattMPC_ubIdx04, plattMPC_dlubcc04);
plattMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_8(plattMPC_ccrhsl05, plattMPC_slb05, plattMPC_llbbyslb05, plattMPC_dzcc05, plattMPC_lbIdx05, plattMPC_dllbcc05);
plattMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_8(plattMPC_ccrhsub05, plattMPC_sub05, plattMPC_lubbysub05, plattMPC_dzcc05, plattMPC_ubIdx05, plattMPC_dlubcc05);
plattMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_8(plattMPC_ccrhsl06, plattMPC_slb06, plattMPC_llbbyslb06, plattMPC_dzcc06, plattMPC_lbIdx06, plattMPC_dllbcc06);
plattMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_8(plattMPC_ccrhsub06, plattMPC_sub06, plattMPC_lubbysub06, plattMPC_dzcc06, plattMPC_ubIdx06, plattMPC_dlubcc06);
plattMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_8(plattMPC_ccrhsl07, plattMPC_slb07, plattMPC_llbbyslb07, plattMPC_dzcc07, plattMPC_lbIdx07, plattMPC_dllbcc07);
plattMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_8(plattMPC_ccrhsub07, plattMPC_sub07, plattMPC_lubbysub07, plattMPC_dzcc07, plattMPC_ubIdx07, plattMPC_dlubcc07);
plattMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_8(plattMPC_ccrhsl08, plattMPC_slb08, plattMPC_llbbyslb08, plattMPC_dzcc08, plattMPC_lbIdx08, plattMPC_dllbcc08);
plattMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_8(plattMPC_ccrhsub08, plattMPC_sub08, plattMPC_lubbysub08, plattMPC_dzcc08, plattMPC_ubIdx08, plattMPC_dlubcc08);
plattMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_8(plattMPC_ccrhsl09, plattMPC_slb09, plattMPC_llbbyslb09, plattMPC_dzcc09, plattMPC_lbIdx09, plattMPC_dllbcc09);
plattMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_8(plattMPC_ccrhsub09, plattMPC_sub09, plattMPC_lubbysub09, plattMPC_dzcc09, plattMPC_ubIdx09, plattMPC_dlubcc09);
plattMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_8(plattMPC_ccrhsl10, plattMPC_slb10, plattMPC_llbbyslb10, plattMPC_dzcc10, plattMPC_lbIdx10, plattMPC_dllbcc10);
plattMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_8(plattMPC_ccrhsub10, plattMPC_sub10, plattMPC_lubbysub10, plattMPC_dzcc10, plattMPC_ubIdx10, plattMPC_dlubcc10);
plattMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_8(plattMPC_ccrhsl11, plattMPC_slb11, plattMPC_llbbyslb11, plattMPC_dzcc11, plattMPC_lbIdx11, plattMPC_dllbcc11);
plattMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_8(plattMPC_ccrhsub11, plattMPC_sub11, plattMPC_lubbysub11, plattMPC_dzcc11, plattMPC_ubIdx11, plattMPC_dlubcc11);
plattMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_8(plattMPC_ccrhsl12, plattMPC_slb12, plattMPC_llbbyslb12, plattMPC_dzcc12, plattMPC_lbIdx12, plattMPC_dllbcc12);
plattMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_8(plattMPC_ccrhsub12, plattMPC_sub12, plattMPC_lubbysub12, plattMPC_dzcc12, plattMPC_ubIdx12, plattMPC_dlubcc12);
plattMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_8(plattMPC_ccrhsl13, plattMPC_slb13, plattMPC_llbbyslb13, plattMPC_dzcc13, plattMPC_lbIdx13, plattMPC_dllbcc13);
plattMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_8(plattMPC_ccrhsub13, plattMPC_sub13, plattMPC_lubbysub13, plattMPC_dzcc13, plattMPC_ubIdx13, plattMPC_dlubcc13);
plattMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_6(plattMPC_ccrhsl14, plattMPC_slb14, plattMPC_llbbyslb14, plattMPC_dzcc14, plattMPC_lbIdx14, plattMPC_dllbcc14);
plattMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_6(plattMPC_ccrhsub14, plattMPC_sub14, plattMPC_lubbysub14, plattMPC_dzcc14, plattMPC_ubIdx14, plattMPC_dlubcc14);
plattMPC_LA_VSUB7_236(plattMPC_l, plattMPC_ccrhs, plattMPC_s, plattMPC_dl_cc, plattMPC_ds_cc);
plattMPC_LA_VADD_118(plattMPC_dz_cc, plattMPC_dz_aff);
plattMPC_LA_VADD_90(plattMPC_dv_cc, plattMPC_dv_aff);
plattMPC_LA_VADD_236(plattMPC_dl_cc, plattMPC_dl_aff);
plattMPC_LA_VADD_236(plattMPC_ds_cc, plattMPC_ds_aff);
info->lsit_cc = plattMPC_LINESEARCH_BACKTRACKING_COMBINED(plattMPC_z, plattMPC_v, plattMPC_l, plattMPC_s, plattMPC_dz_cc, plattMPC_dv_cc, plattMPC_dl_cc, plattMPC_ds_cc, &info->step_cc, &info->mu);
if( info->lsit_cc == plattMPC_NOPROGRESS ){
exitcode = plattMPC_NOPROGRESS; break;
}
info->it++;
}
output->z1[0] = plattMPC_z00[0];
output->z1[1] = plattMPC_z00[1];
output->z1[2] = plattMPC_z00[2];
output->z1[3] = plattMPC_z00[3];
output->z1[4] = plattMPC_z00[4];
output->z1[5] = plattMPC_z00[5];
output->z1[6] = plattMPC_z00[6];
output->z1[7] = plattMPC_z00[7];
output->z2[0] = plattMPC_z01[0];
output->z2[1] = plattMPC_z01[1];
output->z2[2] = plattMPC_z01[2];
output->z2[3] = plattMPC_z01[3];
output->z2[4] = plattMPC_z01[4];
output->z2[5] = plattMPC_z01[5];
output->z2[6] = plattMPC_z01[6];
output->z2[7] = plattMPC_z01[7];
output->z3[0] = plattMPC_z02[0];
output->z3[1] = plattMPC_z02[1];
output->z3[2] = plattMPC_z02[2];
output->z3[3] = plattMPC_z02[3];
output->z3[4] = plattMPC_z02[4];
output->z3[5] = plattMPC_z02[5];
output->z3[6] = plattMPC_z02[6];
output->z3[7] = plattMPC_z02[7];
output->z4[0] = plattMPC_z03[0];
output->z4[1] = plattMPC_z03[1];
output->z4[2] = plattMPC_z03[2];
output->z4[3] = plattMPC_z03[3];
output->z4[4] = plattMPC_z03[4];
output->z4[5] = plattMPC_z03[5];
output->z4[6] = plattMPC_z03[6];
output->z4[7] = plattMPC_z03[7];
output->z5[0] = plattMPC_z04[0];
output->z5[1] = plattMPC_z04[1];
output->z5[2] = plattMPC_z04[2];
output->z5[3] = plattMPC_z04[3];
output->z5[4] = plattMPC_z04[4];
output->z5[5] = plattMPC_z04[5];
output->z5[6] = plattMPC_z04[6];
output->z5[7] = plattMPC_z04[7];
output->z6[0] = plattMPC_z05[0];
output->z6[1] = plattMPC_z05[1];
output->z6[2] = plattMPC_z05[2];
output->z6[3] = plattMPC_z05[3];
output->z6[4] = plattMPC_z05[4];
output->z6[5] = plattMPC_z05[5];
output->z6[6] = plattMPC_z05[6];
output->z6[7] = plattMPC_z05[7];
output->z7[0] = plattMPC_z06[0];
output->z7[1] = plattMPC_z06[1];
output->z7[2] = plattMPC_z06[2];
output->z7[3] = plattMPC_z06[3];
output->z7[4] = plattMPC_z06[4];
output->z7[5] = plattMPC_z06[5];
output->z7[6] = plattMPC_z06[6];
output->z7[7] = plattMPC_z06[7];
output->z8[0] = plattMPC_z07[0];
output->z8[1] = plattMPC_z07[1];
output->z8[2] = plattMPC_z07[2];
output->z8[3] = plattMPC_z07[3];
output->z8[4] = plattMPC_z07[4];
output->z8[5] = plattMPC_z07[5];
output->z8[6] = plattMPC_z07[6];
output->z8[7] = plattMPC_z07[7];
output->z9[0] = plattMPC_z08[0];
output->z9[1] = plattMPC_z08[1];
output->z9[2] = plattMPC_z08[2];
output->z9[3] = plattMPC_z08[3];
output->z9[4] = plattMPC_z08[4];
output->z9[5] = plattMPC_z08[5];
output->z9[6] = plattMPC_z08[6];
output->z9[7] = plattMPC_z08[7];
output->z10[0] = plattMPC_z09[0];
output->z10[1] = plattMPC_z09[1];
output->z10[2] = plattMPC_z09[2];
output->z10[3] = plattMPC_z09[3];
output->z10[4] = plattMPC_z09[4];
output->z10[5] = plattMPC_z09[5];
output->z10[6] = plattMPC_z09[6];
output->z10[7] = plattMPC_z09[7];
output->z11[0] = plattMPC_z10[0];
output->z11[1] = plattMPC_z10[1];
output->z11[2] = plattMPC_z10[2];
output->z11[3] = plattMPC_z10[3];
output->z11[4] = plattMPC_z10[4];
output->z11[5] = plattMPC_z10[5];
output->z11[6] = plattMPC_z10[6];
output->z11[7] = plattMPC_z10[7];
output->z12[0] = plattMPC_z11[0];
output->z12[1] = plattMPC_z11[1];
output->z12[2] = plattMPC_z11[2];
output->z12[3] = plattMPC_z11[3];
output->z12[4] = plattMPC_z11[4];
output->z12[5] = plattMPC_z11[5];
output->z12[6] = plattMPC_z11[6];
output->z12[7] = plattMPC_z11[7];
output->z13[0] = plattMPC_z12[0];
output->z13[1] = plattMPC_z12[1];
output->z13[2] = plattMPC_z12[2];
output->z13[3] = plattMPC_z12[3];
output->z13[4] = plattMPC_z12[4];
output->z13[5] = plattMPC_z12[5];
output->z13[6] = plattMPC_z12[6];
output->z13[7] = plattMPC_z12[7];
output->z14[0] = plattMPC_z13[0];
output->z14[1] = plattMPC_z13[1];
output->z14[2] = plattMPC_z13[2];
output->z14[3] = plattMPC_z13[3];
output->z14[4] = plattMPC_z13[4];
output->z14[5] = plattMPC_z13[5];
output->z14[6] = plattMPC_z13[6];
output->z14[7] = plattMPC_z13[7];
output->z15[0] = plattMPC_z14[0];
output->z15[1] = plattMPC_z14[1];
output->z15[2] = plattMPC_z14[2];
output->z15[3] = plattMPC_z14[3];
output->z15[4] = plattMPC_z14[4];
output->z15[5] = plattMPC_z14[5];

#if plattMPC_SET_TIMING == 1
info->solvetime = plattMPC_toc(&solvertimer);
#if plattMPC_SET_PRINTLEVEL > 0 && plattMPC_SET_TIMING == 1
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
