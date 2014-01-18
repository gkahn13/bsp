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

#include "beliefPenaltyMPC.h"

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
 * Initializes a vector of length 158 with a value.
 */
void beliefPenaltyMPC_LA_INITIALIZEVECTOR_158(beliefPenaltyMPC_FLOAT* vec, beliefPenaltyMPC_FLOAT value)
{
	int i;
	for( i=0; i<158; i++ )
	{
		vec[i] = value;
	}
}


/*
 * Initializes a vector of length 50 with a value.
 */
void beliefPenaltyMPC_LA_INITIALIZEVECTOR_50(beliefPenaltyMPC_FLOAT* vec, beliefPenaltyMPC_FLOAT value)
{
	int i;
	for( i=0; i<50; i++ )
	{
		vec[i] = value;
	}
}


/*
 * Initializes a vector of length 226 with a value.
 */
void beliefPenaltyMPC_LA_INITIALIZEVECTOR_226(beliefPenaltyMPC_FLOAT* vec, beliefPenaltyMPC_FLOAT value)
{
	int i;
	for( i=0; i<226; i++ )
	{
		vec[i] = value;
	}
}


/* 
 * Calculates a dot product and adds it to a variable: z += x'*y; 
 * This function is for vectors of length 226.
 */
void beliefPenaltyMPC_LA_DOTACC_226(beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *y, beliefPenaltyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<226; i++ ){
		*z += x[i]*y[i];
	}
}


/*
 * Calculates the gradient and the value for a quadratic function 0.5*z'*H*z + f'*z
 *
 * INPUTS:     H  - Symmetric Hessian, diag matrix of size [17 x 17]
 *             f  - column vector of size 17
 *             z  - column vector of size 17
 *
 * OUTPUTS: grad  - gradient at z (= H*z + f), column vector of size 17
 *          value <-- value + 0.5*z'*H*z + f'*z (value will be modified)
 */
void beliefPenaltyMPC_LA_DIAG_QUADFCN_17(beliefPenaltyMPC_FLOAT* H, beliefPenaltyMPC_FLOAT* f, beliefPenaltyMPC_FLOAT* z, beliefPenaltyMPC_FLOAT* grad, beliefPenaltyMPC_FLOAT* value)
{
	int i;
	beliefPenaltyMPC_FLOAT hz;	
	for( i=0; i<17; i++){
		hz = H[i]*z[i];
		grad[i] = hz + f[i];
		*value += 0.5*hz*z[i] + f[i]*z[i];
	}
}


/*
 * Calculates the gradient and the value for a quadratic function 0.5*z'*H*z + f'*z
 *
 * INPUTS:     H  - Symmetric Hessian, diag matrix of size [5 x 5]
 *             f  - column vector of size 5
 *             z  - column vector of size 5
 *
 * OUTPUTS: grad  - gradient at z (= H*z + f), column vector of size 5
 *          value <-- value + 0.5*z'*H*z + f'*z (value will be modified)
 */
void beliefPenaltyMPC_LA_DIAG_QUADFCN_5(beliefPenaltyMPC_FLOAT* H, beliefPenaltyMPC_FLOAT* f, beliefPenaltyMPC_FLOAT* z, beliefPenaltyMPC_FLOAT* grad, beliefPenaltyMPC_FLOAT* value)
{
	int i;
	beliefPenaltyMPC_FLOAT hz;	
	for( i=0; i<5; i++){
		hz = H[i]*z[i];
		grad[i] = hz + f[i];
		*value += 0.5*hz*z[i] + f[i]*z[i];
	}
}


/* 
 * Computes r = A*x + B*u - b
 * and      y = max([norm(r,inf), y])
 * and      z -= l'*r
 * where A is stored in column major format
 */
void beliefPenaltyMPC_LA_DENSE_MVMSUB3_10_17_17(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *b, beliefPenaltyMPC_FLOAT *l, beliefPenaltyMPC_FLOAT *r, beliefPenaltyMPC_FLOAT *z, beliefPenaltyMPC_FLOAT *y)
{
	int i;
	int j;
	int k = 0;
	int m = 0;
	int n;
	beliefPenaltyMPC_FLOAT AxBu[10];
	beliefPenaltyMPC_FLOAT norm = *y;
	beliefPenaltyMPC_FLOAT lr = 0;

	/* do A*x + B*u first */
	for( i=0; i<10; i++ ){
		AxBu[i] = A[k++]*x[0] + B[m++]*u[0];
	}	
	for( j=1; j<17; j++ ){		
		for( i=0; i<10; i++ ){
			AxBu[i] += A[k++]*x[j];
		}
	}
	
	for( n=1; n<17; n++ ){
		for( i=0; i<10; i++ ){
			AxBu[i] += B[m++]*u[n];
		}		
	}

	for( i=0; i<10; i++ ){
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
void beliefPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_5_17_17(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *b, beliefPenaltyMPC_FLOAT *l, beliefPenaltyMPC_FLOAT *r, beliefPenaltyMPC_FLOAT *z, beliefPenaltyMPC_FLOAT *y)
{
	int i;
	int j;
	int k = 0;
	beliefPenaltyMPC_FLOAT AxBu[5];
	beliefPenaltyMPC_FLOAT norm = *y;
	beliefPenaltyMPC_FLOAT lr = 0;

	/* do A*x + B*u first */
	for( i=0; i<5; i++ ){
		AxBu[i] = A[k++]*x[0] + B[i]*u[i];
	}	

	for( j=1; j<17; j++ ){		
		for( i=0; i<5; i++ ){
			AxBu[i] += A[k++]*x[j];
		}
	}

	for( i=0; i<5; i++ ){
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
void beliefPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_5_17_5(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *b, beliefPenaltyMPC_FLOAT *l, beliefPenaltyMPC_FLOAT *r, beliefPenaltyMPC_FLOAT *z, beliefPenaltyMPC_FLOAT *y)
{
	int i;
	int j;
	int k = 0;
	beliefPenaltyMPC_FLOAT AxBu[5];
	beliefPenaltyMPC_FLOAT norm = *y;
	beliefPenaltyMPC_FLOAT lr = 0;

	/* do A*x + B*u first */
	for( i=0; i<5; i++ ){
		AxBu[i] = A[k++]*x[0] + B[i]*u[i];
	}	

	for( j=1; j<17; j++ ){		
		for( i=0; i<5; i++ ){
			AxBu[i] += A[k++]*x[j];
		}
	}

	for( i=0; i<5; i++ ){
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
 * Matrix vector multiplication y = M'*x where M is of size [10 x 17]
 * and stored in column major format. Note the transpose of M!
 */
void beliefPenaltyMPC_LA_DENSE_MTVM_10_17(beliefPenaltyMPC_FLOAT *M, beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *y)
{
	int i;
	int j;
	int k = 0; 
	for( i=0; i<17; i++ ){
		y[i] = 0;
		for( j=0; j<10; j++ ){
			y[i] += M[k++]*x[j];
		}
	}
}


/*
 * Matrix vector multiplication z = A'*x + B'*y 
 * where A is of size [5 x 17]
 * and B is of size [10 x 17]
 * and stored in column major format. Note the transposes of A and B!
 */
void beliefPenaltyMPC_LA_DENSE_MTVM2_5_17_10(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *y, beliefPenaltyMPC_FLOAT *z)
{
	int i;
	int j;
	int k = 0;
	int n;
	int m = 0;
	for( i=0; i<17; i++ ){
		z[i] = 0;
		for( j=0; j<5; j++ ){
			z[i] += A[k++]*x[j];
		}
		for( n=0; n<10; n++ ){
			z[i] += B[m++]*y[n];
		}
	}
}


/*
 * Matrix vector multiplication z = A'*x + B'*y 
 * where A is of size [5 x 17] and stored in column major format.
 * and B is of size [5 x 17] and stored in diagzero format
 * Note the transposes of A and B!
 */
void beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_5_17_5(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *y, beliefPenaltyMPC_FLOAT *z)
{
	int i;
	int j;
	int k = 0;
	for( i=0; i<5; i++ ){
		z[i] = 0;
		for( j=0; j<5; j++ ){
			z[i] += A[k++]*x[j];
		}
		z[i] += B[i]*y[i];
	}
	for( i=5 ;i<17; i++ ){
		z[i] = 0;
		for( j=0; j<5; j++ ){
			z[i] += A[k++]*x[j];
		}
	}
}


/*
 * Matrix vector multiplication y = M'*x where M is of size [5 x 5]
 * and stored in diagzero format. Note the transpose of M!
 */
void beliefPenaltyMPC_LA_DIAGZERO_MTVM_5_5(beliefPenaltyMPC_FLOAT *M, beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *y)
{
	int i;
	for( i=0; i<5; i++ ){
		y[i] = M[i]*x[i];
	}
}


/*
 * Vector subtraction and addition.
 *	 Input: five vectors t, tidx, u, v, w and two scalars z and r
 *	 Output: y = t(tidx) - u + w
 *           z = z - v'*x;
 *           r = max([norm(y,inf), z]);
 * for vectors of length 17. Output z is of course scalar.
 */
void beliefPenaltyMPC_LA_VSUBADD3_17(beliefPenaltyMPC_FLOAT* t, beliefPenaltyMPC_FLOAT* u, int* uidx, beliefPenaltyMPC_FLOAT* v, beliefPenaltyMPC_FLOAT* w, beliefPenaltyMPC_FLOAT* y, beliefPenaltyMPC_FLOAT* z, beliefPenaltyMPC_FLOAT* r)
{
	int i;
	beliefPenaltyMPC_FLOAT norm = *r;
	beliefPenaltyMPC_FLOAT vx = 0;
	beliefPenaltyMPC_FLOAT x;
	for( i=0; i<17; i++){
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
void beliefPenaltyMPC_LA_VSUBADD2_7(beliefPenaltyMPC_FLOAT* t, int* tidx, beliefPenaltyMPC_FLOAT* u, beliefPenaltyMPC_FLOAT* v, beliefPenaltyMPC_FLOAT* w, beliefPenaltyMPC_FLOAT* y, beliefPenaltyMPC_FLOAT* z, beliefPenaltyMPC_FLOAT* r)
{
	int i;
	beliefPenaltyMPC_FLOAT norm = *r;
	beliefPenaltyMPC_FLOAT vx = 0;
	beliefPenaltyMPC_FLOAT x;
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
 * Vector subtraction and addition.
 *	 Input: five vectors t, tidx, u, v, w and two scalars z and r
 *	 Output: y = t(tidx) - u + w
 *           z = z - v'*x;
 *           r = max([norm(y,inf), z]);
 * for vectors of length 5. Output z is of course scalar.
 */
void beliefPenaltyMPC_LA_VSUBADD3_5(beliefPenaltyMPC_FLOAT* t, beliefPenaltyMPC_FLOAT* u, int* uidx, beliefPenaltyMPC_FLOAT* v, beliefPenaltyMPC_FLOAT* w, beliefPenaltyMPC_FLOAT* y, beliefPenaltyMPC_FLOAT* z, beliefPenaltyMPC_FLOAT* r)
{
	int i;
	beliefPenaltyMPC_FLOAT norm = *r;
	beliefPenaltyMPC_FLOAT vx = 0;
	beliefPenaltyMPC_FLOAT x;
	for( i=0; i<5; i++){
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
void beliefPenaltyMPC_LA_VSUBADD2_5(beliefPenaltyMPC_FLOAT* t, int* tidx, beliefPenaltyMPC_FLOAT* u, beliefPenaltyMPC_FLOAT* v, beliefPenaltyMPC_FLOAT* w, beliefPenaltyMPC_FLOAT* y, beliefPenaltyMPC_FLOAT* z, beliefPenaltyMPC_FLOAT* r)
{
	int i;
	beliefPenaltyMPC_FLOAT norm = *r;
	beliefPenaltyMPC_FLOAT vx = 0;
	beliefPenaltyMPC_FLOAT x;
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
 * Computes inequality constraints gradient-
 * Special function for box constraints of length 17
 * Returns also L/S, a value that is often used elsewhere.
 */
void beliefPenaltyMPC_LA_INEQ_B_GRAD_17_17_7(beliefPenaltyMPC_FLOAT *lu, beliefPenaltyMPC_FLOAT *su, beliefPenaltyMPC_FLOAT *ru, beliefPenaltyMPC_FLOAT *ll, beliefPenaltyMPC_FLOAT *sl, beliefPenaltyMPC_FLOAT *rl, int* lbIdx, int* ubIdx, beliefPenaltyMPC_FLOAT *grad, beliefPenaltyMPC_FLOAT *lubysu, beliefPenaltyMPC_FLOAT *llbysl)
{
	int i;
	for( i=0; i<17; i++ ){
		grad[i] = 0;
	}
	for( i=0; i<17; i++ ){		
		llbysl[i] = ll[i] / sl[i];
		grad[lbIdx[i]] -= llbysl[i]*rl[i];
	}
	for( i=0; i<7; i++ ){
		lubysu[i] = lu[i] / su[i];
		grad[ubIdx[i]] += lubysu[i]*ru[i];
	}
}


/*
 * Computes inequality constraints gradient-
 * Special function for box constraints of length 5
 * Returns also L/S, a value that is often used elsewhere.
 */
void beliefPenaltyMPC_LA_INEQ_B_GRAD_5_5_5(beliefPenaltyMPC_FLOAT *lu, beliefPenaltyMPC_FLOAT *su, beliefPenaltyMPC_FLOAT *ru, beliefPenaltyMPC_FLOAT *ll, beliefPenaltyMPC_FLOAT *sl, beliefPenaltyMPC_FLOAT *rl, int* lbIdx, int* ubIdx, beliefPenaltyMPC_FLOAT *grad, beliefPenaltyMPC_FLOAT *lubysu, beliefPenaltyMPC_FLOAT *llbysl)
{
	int i;
	for( i=0; i<5; i++ ){
		grad[i] = 0;
	}
	for( i=0; i<5; i++ ){		
		llbysl[i] = ll[i] / sl[i];
		grad[lbIdx[i]] -= llbysl[i]*rl[i];
	}
	for( i=0; i<5; i++ ){
		lubysu[i] = lu[i] / su[i];
		grad[ubIdx[i]] += lubysu[i]*ru[i];
	}
}


/*
 * Addition of three vectors  z = u + w + v
 * of length 158.
 */
void beliefPenaltyMPC_LA_VVADD3_158(beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *v, beliefPenaltyMPC_FLOAT *w, beliefPenaltyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<158; i++ ){
		z[i] = u[i] + v[i] + w[i];
	}
}


/*
 * Special function to compute the diagonal cholesky factorization of the 
 * positive definite augmented Hessian for block size 17.
 *
 * Inputs: - H = diagonal cost Hessian in diagonal storage format
 *         - llbysl = L / S of lower bounds
 *         - lubysu = L / S of upper bounds
 *
 * Output: Phi = sqrt(H + diag(llbysl) + diag(lubysu))
 * where Phi is stored in diagonal storage format
 */
void beliefPenaltyMPC_LA_DIAG_CHOL_LBUB_17_17_7(beliefPenaltyMPC_FLOAT *H, beliefPenaltyMPC_FLOAT *llbysl, int* lbIdx, beliefPenaltyMPC_FLOAT *lubysu, int* ubIdx, beliefPenaltyMPC_FLOAT *Phi)


{
	int i;
	
	/* copy  H into PHI */
	for( i=0; i<17; i++ ){
		Phi[i] = H[i];
	}

	/* add llbysl onto Phi where necessary */
	for( i=0; i<17; i++ ){
		Phi[lbIdx[i]] += llbysl[i];
	}

	/* add lubysu onto Phi where necessary */
	for( i=0; i<7; i++){
		Phi[ubIdx[i]] +=  lubysu[i];
	}
	
	/* compute cholesky */
	for(i=0; i<17; i++)
	{
#if beliefPenaltyMPC_SET_PRINTLEVEL > 0 && defined PRINTNUMERICALWARNINGS
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
 * where A is to be computed and is of size [10 x 17],
 * B is given and of size [10 x 17], L is a diagonal
 * matrix of size 10 stored in diagonal matrix 
 * storage format. Note the transpose of L has no impact!
 *
 * Result: A in column major storage format.
 *
 */
void beliefPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_10_17(beliefPenaltyMPC_FLOAT *L, beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *A)
{
    int i,j;
	 int k = 0;

	for( j=0; j<17; j++){
		for( i=0; i<10; i++){
			A[k] = B[k]/L[j];
			k++;
		}
	}

}


/**
 * Forward substitution to solve L*y = b where L is a
 * diagonal matrix in vector storage format.
 * 
 * The dimensions involved are 17.
 */
void beliefPenaltyMPC_LA_DIAG_FORWARDSUB_17(beliefPenaltyMPC_FLOAT *L, beliefPenaltyMPC_FLOAT *b, beliefPenaltyMPC_FLOAT *y)
{
    int i;

    for( i=0; i<17; i++ ){
		y[i] = b[i]/L[i];
    }
}


/**
 * Forward substitution for the matrix equation A*L' = B
 * where A is to be computed and is of size [5 x 17],
 * B is given and of size [5 x 17], L is a diagonal
 * matrix of size 5 stored in diagonal matrix 
 * storage format. Note the transpose of L has no impact!
 *
 * Result: A in column major storage format.
 *
 */
void beliefPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_5_17(beliefPenaltyMPC_FLOAT *L, beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *A)
{
    int i,j;
	 int k = 0;

	for( j=0; j<17; j++){
		for( i=0; i<5; i++){
			A[k] = B[k]/L[j];
			k++;
		}
	}

}


/**
 * Compute C = A*B' where 
 *
 *	size(A) = [10 x 17]
 *  size(B) = [5 x 17]
 * 
 * and all matrices are stored in column major format.
 *
 * THIS ONE HAS THE WORST ACCES PATTERN POSSIBLE.  
 * 
 */
void beliefPenaltyMPC_LA_DENSE_MMTM_10_17_5(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *C)
{
    int i, j, k;
    beliefPenaltyMPC_FLOAT temp;
    
    for( i=0; i<10; i++ ){        
        for( j=0; j<5; j++ ){
            temp = 0; 
            for( k=0; k<17; k++ ){
                temp += A[k*10+i]*B[k*5+j];
            }						
            C[j*10+i] = temp;
        }
    }
}


/**
 * Forward substitution for the matrix equation A*L' = B
 * where A is to be computed and is of size [5 x 17],
 * B is given and of size [5 x 17], L is a diagonal
 *  matrix of size 17 stored in diagonal 
 * storage format. Note the transpose of L!
 *
 * Result: A in diagonalzero storage format.
 *
 */
void beliefPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_5_17(beliefPenaltyMPC_FLOAT *L, beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *A)
{
	int j;
    for( j=0; j<17; j++ ){   
		A[j] = B[j]/L[j];
     }
}


/**
 * Compute C = A*B' where 
 *
 *	size(A) = [5 x 17]
 *  size(B) = [5 x 17] in diagzero format
 * 
 * A and C matrices are stored in column major format.
 * 
 * 
 */
void beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMTM_5_17_5(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *C)
{
    int i, j;
	
	for( i=0; i<5; i++ ){
		for( j=0; j<5; j++){
			C[j*5+i] = B[i*5+j]*A[i];
		}
	}

}


/*
 * Special function to compute the diagonal cholesky factorization of the 
 * positive definite augmented Hessian for block size 5.
 *
 * Inputs: - H = diagonal cost Hessian in diagonal storage format
 *         - llbysl = L / S of lower bounds
 *         - lubysu = L / S of upper bounds
 *
 * Output: Phi = sqrt(H + diag(llbysl) + diag(lubysu))
 * where Phi is stored in diagonal storage format
 */
void beliefPenaltyMPC_LA_DIAG_CHOL_ONELOOP_LBUB_5_5_5(beliefPenaltyMPC_FLOAT *H, beliefPenaltyMPC_FLOAT *llbysl, int* lbIdx, beliefPenaltyMPC_FLOAT *lubysu, int* ubIdx, beliefPenaltyMPC_FLOAT *Phi)


{
	int i;
	
	/* compute cholesky */
	for( i=0; i<5; i++ ){
		Phi[i] = H[i] + llbysl[i] + lubysu[i];

#if beliefPenaltyMPC_SET_PRINTLEVEL > 0 && defined PRINTNUMERICALWARNINGS
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
 * where A is to be computed and is of size [5 x 5],
 * B is given and of size [5 x 5], L is a diagonal
 *  matrix of size 5 stored in diagonal 
 * storage format. Note the transpose of L!
 *
 * Result: A in diagonalzero storage format.
 *
 */
void beliefPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_5_5(beliefPenaltyMPC_FLOAT *L, beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *A)
{
	int j;
    for( j=0; j<5; j++ ){   
		A[j] = B[j]/L[j];
     }
}


/**
 * Forward substitution to solve L*y = b where L is a
 * diagonal matrix in vector storage format.
 * 
 * The dimensions involved are 5.
 */
void beliefPenaltyMPC_LA_DIAG_FORWARDSUB_5(beliefPenaltyMPC_FLOAT *L, beliefPenaltyMPC_FLOAT *b, beliefPenaltyMPC_FLOAT *y)
{
    int i;

    for( i=0; i<5; i++ ){
		y[i] = b[i]/L[i];
    }
}


/**
 * Compute L = A*A' + B*B', where L is lower triangular of size NXp1
 * and A is a dense matrix of size [10 x 17] in column
 * storage format, and B is of size [10 x 17] also in column
 * storage format.
 * 
 * THIS ONE HAS THE WORST ACCES PATTERN POSSIBLE. 
 * POSSIBKE FIX: PUT A AND B INTO ROW MAJOR FORMAT FIRST.
 * 
 */
void beliefPenaltyMPC_LA_DENSE_MMT2_10_17_17(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *L)
{
    int i, j, k, ii, di;
    beliefPenaltyMPC_FLOAT ltemp;
    
    ii = 0; di = 0;
    for( i=0; i<10; i++ ){        
        for( j=0; j<=i; j++ ){
            ltemp = 0; 
            for( k=0; k<17; k++ ){
                ltemp += A[k*10+i]*A[k*10+j];
            }			
			for( k=0; k<17; k++ ){
                ltemp += B[k*10+i]*B[k*10+j];
            }
            L[ii+j] = ltemp;
        }
        ii += ++di;
    }
}


/* 
 * Computes r = b - A*x - B*u
 * where A an B are stored in column major format
 */
void beliefPenaltyMPC_LA_DENSE_MVMSUB2_10_17_17(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *b, beliefPenaltyMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;
	int m = 0;
	int n;

	for( i=0; i<10; i++ ){
		r[i] = b[i] - A[k++]*x[0] - B[m++]*u[0];
	}	
	for( j=1; j<17; j++ ){		
		for( i=0; i<10; i++ ){
			r[i] -= A[k++]*x[j];
		}
	}
	
	for( n=1; n<17; n++ ){
		for( i=0; i<10; i++ ){
			r[i] -= B[m++]*u[n];
		}		
	}
}


/**
 * Compute L = A*A' + B*B', where L is lower triangular of size NXp1
 * and A is a dense matrix of size [5 x 17] in column
 * storage format, and B is of size [5 x 17] diagonalzero
 * storage format.
 * 
 * THIS ONE HAS THE WORST ACCES PATTERN POSSIBLE. 
 * POSSIBKE FIX: PUT A AND B INTO ROW MAJOR FORMAT FIRST.
 * 
 */
void beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMT2_5_17_17(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *L)
{
    int i, j, k, ii, di;
    beliefPenaltyMPC_FLOAT ltemp;
    
    ii = 0; di = 0;
    for( i=0; i<5; i++ ){        
        for( j=0; j<=i; j++ ){
            ltemp = 0; 
            for( k=0; k<17; k++ ){
                ltemp += A[k*5+i]*A[k*5+j];
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
void beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_5_17_17(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *b, beliefPenaltyMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;

	for( i=0; i<5; i++ ){
		r[i] = b[i] - A[k++]*x[0] - B[i]*u[i];
	}	

	for( j=1; j<17; j++ ){		
		for( i=0; i<5; i++ ){
			r[i] -= A[k++]*x[j];
		}
	}
	
}


/**
 * Compute L = A*A' + B*B', where L is lower triangular of size NXp1
 * and A is a dense matrix of size [5 x 17] in column
 * storage format, and B is of size [5 x 5] diagonalzero
 * storage format.
 * 
 * THIS ONE HAS THE WORST ACCES PATTERN POSSIBLE. 
 * POSSIBKE FIX: PUT A AND B INTO ROW MAJOR FORMAT FIRST.
 * 
 */
void beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMT2_5_17_5(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *L)
{
    int i, j, k, ii, di;
    beliefPenaltyMPC_FLOAT ltemp;
    
    ii = 0; di = 0;
    for( i=0; i<5; i++ ){        
        for( j=0; j<=i; j++ ){
            ltemp = 0; 
            for( k=0; k<17; k++ ){
                ltemp += A[k*5+i]*A[k*5+j];
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
void beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_5_17_5(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *b, beliefPenaltyMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;

	for( i=0; i<5; i++ ){
		r[i] = b[i] - A[k++]*x[0] - B[i]*u[i];
	}	

	for( j=1; j<17; j++ ){		
		for( i=0; i<5; i++ ){
			r[i] -= A[k++]*x[j];
		}
	}
	
}


/**
 * Cholesky factorization as above, but working on a matrix in 
 * lower triangular storage format of size 10 and outputting
 * the Cholesky factor to matrix L in lower triangular format.
 */
void beliefPenaltyMPC_LA_DENSE_CHOL_10(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *L)
{
    int i, j, k, di, dj;
	 int ii, jj;

    beliefPenaltyMPC_FLOAT l;
    beliefPenaltyMPC_FLOAT Mii;

	/* copy A to L first and then operate on L */
	/* COULD BE OPTIMIZED */
	ii=0; di=0;
	for( i=0; i<10; i++ ){
		for( j=0; j<=i; j++ ){
			L[ii+j] = A[ii+j];
		}
		ii += ++di;
	}    
	
	/* factor L */
	ii=0; di=0;
    for( i=0; i<10; i++ ){
        l = 0;
        for( k=0; k<i; k++ ){
            l += L[ii+k]*L[ii+k];
        }        
        
        Mii = L[ii+i] - l;
        
#if beliefPenaltyMPC_SET_PRINTLEVEL > 0 && defined PRINTNUMERICALWARNINGS
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
        for( j=i+1; j<10; j++ ){
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
 * The dimensions involved are 10.
 */
void beliefPenaltyMPC_LA_DENSE_FORWARDSUB_10(beliefPenaltyMPC_FLOAT *L, beliefPenaltyMPC_FLOAT *b, beliefPenaltyMPC_FLOAT *y)
{
    int i,j,ii,di;
    beliefPenaltyMPC_FLOAT yel;
            
    ii = 0; di = 0;
    for( i=0; i<10; i++ ){
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
 * where A is to be computed and is of size [5 x 10],
 * B is given and of size [5 x 10], L is a lower tri-
 * angular matrix of size 10 stored in lower triangular 
 * storage format. Note the transpose of L AND B!
 *
 * Result: A in column major storage format.
 *
 */
void beliefPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_5_10(beliefPenaltyMPC_FLOAT *L, beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *A)
{
    int i,j,k,ii,di;
    beliefPenaltyMPC_FLOAT a;
    
    ii=0; di=0;
    for( j=0; j<10; j++ ){        
        for( i=0; i<5; i++ ){
            a = B[i*10+j];
            for( k=0; k<j; k++ ){
                a -= A[k*5+i]*L[ii+k];
            }    

			/* saturate for numerical stability */
			a = MIN(a, BIGM);
			a = MAX(a, -BIGM); 

			A[j*5+i] = a/L[ii+j];			
        }
        ii += ++di;
    }
}


/**
 * Compute L = L - A*A', where L is lower triangular of size 5
 * and A is a dense matrix of size [5 x 10] in column
 * storage format.
 * 
 * THIS ONE HAS THE WORST ACCES PATTERN POSSIBLE. 
 * POSSIBKE FIX: PUT A INTO ROW MAJOR FORMAT FIRST.
 * 
 */
void beliefPenaltyMPC_LA_DENSE_MMTSUB_5_10(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *L)
{
    int i, j, k, ii, di;
    beliefPenaltyMPC_FLOAT ltemp;
    
    ii = 0; di = 0;
    for( i=0; i<5; i++ ){        
        for( j=0; j<=i; j++ ){
            ltemp = 0; 
            for( k=0; k<10; k++ ){
                ltemp += A[k*5+i]*A[k*5+j];
            }						
            L[ii+j] -= ltemp;
        }
        ii += ++di;
    }
}


/**
 * Cholesky factorization as above, but working on a matrix in 
 * lower triangular storage format of size 5 and outputting
 * the Cholesky factor to matrix L in lower triangular format.
 */
void beliefPenaltyMPC_LA_DENSE_CHOL_5(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *L)
{
    int i, j, k, di, dj;
	 int ii, jj;

    beliefPenaltyMPC_FLOAT l;
    beliefPenaltyMPC_FLOAT Mii;

	/* copy A to L first and then operate on L */
	/* COULD BE OPTIMIZED */
	ii=0; di=0;
	for( i=0; i<5; i++ ){
		for( j=0; j<=i; j++ ){
			L[ii+j] = A[ii+j];
		}
		ii += ++di;
	}    
	
	/* factor L */
	ii=0; di=0;
    for( i=0; i<5; i++ ){
        l = 0;
        for( k=0; k<i; k++ ){
            l += L[ii+k]*L[ii+k];
        }        
        
        Mii = L[ii+i] - l;
        
#if beliefPenaltyMPC_SET_PRINTLEVEL > 0 && defined PRINTNUMERICALWARNINGS
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
        for( j=i+1; j<5; j++ ){
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


/* 
 * Computes r = b - A*x
 * where A is stored in column major format
 */
void beliefPenaltyMPC_LA_DENSE_MVMSUB1_5_10(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *b, beliefPenaltyMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;

	for( i=0; i<5; i++ ){
		r[i] = b[i] - A[k++]*x[0];
	}	
	for( j=1; j<10; j++ ){		
		for( i=0; i<5; i++ ){
			r[i] -= A[k++]*x[j];
		}
	}
}


/**
 * Forward substitution to solve L*y = b where L is a
 * lower triangular matrix in triangular storage format.
 * 
 * The dimensions involved are 5.
 */
void beliefPenaltyMPC_LA_DENSE_FORWARDSUB_5(beliefPenaltyMPC_FLOAT *L, beliefPenaltyMPC_FLOAT *b, beliefPenaltyMPC_FLOAT *y)
{
    int i,j,ii,di;
    beliefPenaltyMPC_FLOAT yel;
            
    ii = 0; di = 0;
    for( i=0; i<5; i++ ){
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
 * where A is to be computed and is of size [5 x 5],
 * B is given and of size [5 x 5], L is a lower tri-
 * angular matrix of size 5 stored in lower triangular 
 * storage format. Note the transpose of L AND B!
 *
 * Result: A in column major storage format.
 *
 */
void beliefPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_5_5(beliefPenaltyMPC_FLOAT *L, beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *A)
{
    int i,j,k,ii,di;
    beliefPenaltyMPC_FLOAT a;
    
    ii=0; di=0;
    for( j=0; j<5; j++ ){        
        for( i=0; i<5; i++ ){
            a = B[i*5+j];
            for( k=0; k<j; k++ ){
                a -= A[k*5+i]*L[ii+k];
            }    

			/* saturate for numerical stability */
			a = MIN(a, BIGM);
			a = MAX(a, -BIGM); 

			A[j*5+i] = a/L[ii+j];			
        }
        ii += ++di;
    }
}


/**
 * Compute L = L - A*A', where L is lower triangular of size 5
 * and A is a dense matrix of size [5 x 5] in column
 * storage format.
 * 
 * THIS ONE HAS THE WORST ACCES PATTERN POSSIBLE. 
 * POSSIBKE FIX: PUT A INTO ROW MAJOR FORMAT FIRST.
 * 
 */
void beliefPenaltyMPC_LA_DENSE_MMTSUB_5_5(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *L)
{
    int i, j, k, ii, di;
    beliefPenaltyMPC_FLOAT ltemp;
    
    ii = 0; di = 0;
    for( i=0; i<5; i++ ){        
        for( j=0; j<=i; j++ ){
            ltemp = 0; 
            for( k=0; k<5; k++ ){
                ltemp += A[k*5+i]*A[k*5+j];
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
void beliefPenaltyMPC_LA_DENSE_MVMSUB1_5_5(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *b, beliefPenaltyMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;

	for( i=0; i<5; i++ ){
		r[i] = b[i] - A[k++]*x[0];
	}	
	for( j=1; j<5; j++ ){		
		for( i=0; i<5; i++ ){
			r[i] -= A[k++]*x[j];
		}
	}
}


/**
 * Backward Substitution to solve L^T*x = y where L is a
 * lower triangular matrix in triangular storage format.
 * 
 * All involved dimensions are 5.
 */
void beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_5(beliefPenaltyMPC_FLOAT *L, beliefPenaltyMPC_FLOAT *y, beliefPenaltyMPC_FLOAT *x)
{
    int i, ii, di, j, jj, dj;
    beliefPenaltyMPC_FLOAT xel;    
	int start = 10;
    
    /* now solve L^T*x = y by backward substitution */
    ii = start; di = 4;
    for( i=4; i>=0; i-- ){        
        xel = y[i];        
        jj = start; dj = 4;
        for( j=4; j>i; j-- ){
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
 * Matrix vector multiplication y = b - M'*x where M is of size [5 x 5]
 * and stored in column major format. Note the transpose of M!
 */
void beliefPenaltyMPC_LA_DENSE_MTVMSUB_5_5(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *b, beliefPenaltyMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0; 
	for( i=0; i<5; i++ ){
		r[i] = b[i];
		for( j=0; j<5; j++ ){
			r[i] -= A[k++]*x[j];
		}
	}
}


/*
 * Matrix vector multiplication y = b - M'*x where M is of size [5 x 10]
 * and stored in column major format. Note the transpose of M!
 */
void beliefPenaltyMPC_LA_DENSE_MTVMSUB_5_10(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *b, beliefPenaltyMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0; 
	for( i=0; i<10; i++ ){
		r[i] = b[i];
		for( j=0; j<5; j++ ){
			r[i] -= A[k++]*x[j];
		}
	}
}


/**
 * Backward Substitution to solve L^T*x = y where L is a
 * lower triangular matrix in triangular storage format.
 * 
 * All involved dimensions are 10.
 */
void beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_10(beliefPenaltyMPC_FLOAT *L, beliefPenaltyMPC_FLOAT *y, beliefPenaltyMPC_FLOAT *x)
{
    int i, ii, di, j, jj, dj;
    beliefPenaltyMPC_FLOAT xel;    
	int start = 45;
    
    /* now solve L^T*x = y by backward substitution */
    ii = start; di = 9;
    for( i=9; i>=0; i-- ){        
        xel = y[i];        
        jj = start; dj = 9;
        for( j=9; j>i; j-- ){
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
 * Vector subtraction z = -x - y for vectors of length 158.
 */
void beliefPenaltyMPC_LA_VSUB2_158(beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *y, beliefPenaltyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<158; i++){
		z[i] = -x[i] - y[i];
	}
}


/**
 * Forward-Backward-Substitution to solve L*L^T*x = b where L is a
 * diagonal matrix of size 17 in vector
 * storage format.
 */
void beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_17(beliefPenaltyMPC_FLOAT *L, beliefPenaltyMPC_FLOAT *b, beliefPenaltyMPC_FLOAT *x)
{
    int i;
            
    /* solve Ly = b by forward and backward substitution */
    for( i=0; i<17; i++ ){
		x[i] = b[i]/(L[i]*L[i]);
    }
    
}


/**
 * Forward-Backward-Substitution to solve L*L^T*x = b where L is a
 * diagonal matrix of size 5 in vector
 * storage format.
 */
void beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_5(beliefPenaltyMPC_FLOAT *L, beliefPenaltyMPC_FLOAT *b, beliefPenaltyMPC_FLOAT *x)
{
    int i;
            
    /* solve Ly = b by forward and backward substitution */
    for( i=0; i<5; i++ ){
		x[i] = b[i]/(L[i]*L[i]);
    }
    
}


/*
 * Vector subtraction z = x(xidx) - y where y, z and xidx are of length 17,
 * and x has length 17 and is indexed through yidx.
 */
void beliefPenaltyMPC_LA_VSUB_INDEXED_17(beliefPenaltyMPC_FLOAT *x, int* xidx, beliefPenaltyMPC_FLOAT *y, beliefPenaltyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<17; i++){
		z[i] = x[xidx[i]] - y[i];
	}
}


/*
 * Vector subtraction x = -u.*v - w for vectors of length 17.
 */
void beliefPenaltyMPC_LA_VSUB3_17(beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *v, beliefPenaltyMPC_FLOAT *w, beliefPenaltyMPC_FLOAT *x)
{
	int i;
	for( i=0; i<17; i++){
		x[i] = -u[i]*v[i] - w[i];
	}
}


/*
 * Vector subtraction z = -x - y(yidx) where y is of length 17
 * and z, x and yidx are of length 7.
 */
void beliefPenaltyMPC_LA_VSUB2_INDEXED_7(beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *y, int* yidx, beliefPenaltyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<7; i++){
		z[i] = -x[i] - y[yidx[i]];
	}
}


/*
 * Vector subtraction x = -u.*v - w for vectors of length 7.
 */
void beliefPenaltyMPC_LA_VSUB3_7(beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *v, beliefPenaltyMPC_FLOAT *w, beliefPenaltyMPC_FLOAT *x)
{
	int i;
	for( i=0; i<7; i++){
		x[i] = -u[i]*v[i] - w[i];
	}
}


/*
 * Vector subtraction z = x(xidx) - y where y, z and xidx are of length 5,
 * and x has length 5 and is indexed through yidx.
 */
void beliefPenaltyMPC_LA_VSUB_INDEXED_5(beliefPenaltyMPC_FLOAT *x, int* xidx, beliefPenaltyMPC_FLOAT *y, beliefPenaltyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<5; i++){
		z[i] = x[xidx[i]] - y[i];
	}
}


/*
 * Vector subtraction x = -u.*v - w for vectors of length 5.
 */
void beliefPenaltyMPC_LA_VSUB3_5(beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *v, beliefPenaltyMPC_FLOAT *w, beliefPenaltyMPC_FLOAT *x)
{
	int i;
	for( i=0; i<5; i++){
		x[i] = -u[i]*v[i] - w[i];
	}
}


/*
 * Vector subtraction z = -x - y(yidx) where y is of length 5
 * and z, x and yidx are of length 5.
 */
void beliefPenaltyMPC_LA_VSUB2_INDEXED_5(beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *y, int* yidx, beliefPenaltyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<5; i++){
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
 * beliefPenaltyMPC_NOPROGRESS (should be negative).
 */
int beliefPenaltyMPC_LINESEARCH_BACKTRACKING_AFFINE(beliefPenaltyMPC_FLOAT *l, beliefPenaltyMPC_FLOAT *s, beliefPenaltyMPC_FLOAT *dl, beliefPenaltyMPC_FLOAT *ds, beliefPenaltyMPC_FLOAT *a, beliefPenaltyMPC_FLOAT *mu_aff)
{
    int i;
	int lsIt=1;    
    beliefPenaltyMPC_FLOAT dltemp;
    beliefPenaltyMPC_FLOAT dstemp;
    beliefPenaltyMPC_FLOAT mya = 1.0;
    beliefPenaltyMPC_FLOAT mymu;
        
    while( 1 ){                        

        /* 
         * Compute both snew and wnew together.
         * We compute also mu_affine along the way here, as the
         * values might be in registers, so it should be cheaper.
         */
        mymu = 0;
        for( i=0; i<226; i++ ){
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
        if( i == 226 ){
            break;
        } else {
            mya *= beliefPenaltyMPC_SET_LS_SCALE_AFF;
            if( mya < beliefPenaltyMPC_SET_LS_MINSTEP ){
                return beliefPenaltyMPC_NOPROGRESS;
            }
        }
    }
    
    /* return new values and iteration counter */
    *a = mya;
    *mu_aff = mymu / (beliefPenaltyMPC_FLOAT)226;
    return lsIt;
}


/*
 * Vector subtraction x = u.*v - a where a is a scalar
*  and x,u,v are vectors of length 226.
 */
void beliefPenaltyMPC_LA_VSUB5_226(beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *v, beliefPenaltyMPC_FLOAT a, beliefPenaltyMPC_FLOAT *x)
{
	int i;
	for( i=0; i<226; i++){
		x[i] = u[i]*v[i] - a;
	}
}


/*
 * Computes x=0; x(uidx) += u/su; x(vidx) -= v/sv where x is of length 17,
 * u, su, uidx are of length 7 and v, sv, vidx are of length 17.
 */
void beliefPenaltyMPC_LA_VSUB6_INDEXED_17_7_17(beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *su, int* uidx, beliefPenaltyMPC_FLOAT *v, beliefPenaltyMPC_FLOAT *sv, int* vidx, beliefPenaltyMPC_FLOAT *x)
{
	int i;
	for( i=0; i<17; i++ ){
		x[i] = 0;
	}
	for( i=0; i<7; i++){
		x[uidx[i]] += u[i]/su[i];
	}
	for( i=0; i<17; i++){
		x[vidx[i]] -= v[i]/sv[i];
	}
}


/* 
 * Computes r = A*x + B*u
 * where A an B are stored in column major format
 */
void beliefPenaltyMPC_LA_DENSE_2MVMADD_10_17_17(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;
	int m = 0;
	int n;

	for( i=0; i<10; i++ ){
		r[i] = A[k++]*x[0] + B[m++]*u[0];
	}	

	for( j=1; j<17; j++ ){		
		for( i=0; i<10; i++ ){
			r[i] += A[k++]*x[j];
		}
	}
	
	for( n=1; n<17; n++ ){
		for( i=0; i<10; i++ ){
			r[i] += B[m++]*u[n];
		}		
	}
}


/* 
 * Computes r = A*x + B*u
 * where A is stored in column major format
 * and B is stored in diagzero format
 */
void beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_5_17_17(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;

	for( i=0; i<5; i++ ){
		r[i] = A[k++]*x[0] + B[i]*u[i];
	}	

	for( j=1; j<17; j++ ){		
		for( i=0; i<5; i++ ){
			r[i] += A[k++]*x[j];
		}
	}
	
}


/*
 * Computes x=0; x(uidx) += u/su; x(vidx) -= v/sv where x is of length 5,
 * u, su, uidx are of length 5 and v, sv, vidx are of length 5.
 */
void beliefPenaltyMPC_LA_VSUB6_INDEXED_5_5_5(beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *su, int* uidx, beliefPenaltyMPC_FLOAT *v, beliefPenaltyMPC_FLOAT *sv, int* vidx, beliefPenaltyMPC_FLOAT *x)
{
	int i;
	for( i=0; i<5; i++ ){
		x[i] = 0;
	}
	for( i=0; i<5; i++){
		x[uidx[i]] += u[i]/su[i];
	}
	for( i=0; i<5; i++){
		x[vidx[i]] -= v[i]/sv[i];
	}
}


/* 
 * Computes r = A*x + B*u
 * where A is stored in column major format
 * and B is stored in diagzero format
 */
void beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_5_17_5(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;

	for( i=0; i<5; i++ ){
		r[i] = A[k++]*x[0] + B[i]*u[i];
	}	

	for( j=1; j<17; j++ ){		
		for( i=0; i<5; i++ ){
			r[i] += A[k++]*x[j];
		}
	}
	
}


/*
 * Vector subtraction z = x - y for vectors of length 158.
 */
void beliefPenaltyMPC_LA_VSUB_158(beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *y, beliefPenaltyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<158; i++){
		z[i] = x[i] - y[i];
	}
}


/** 
 * Computes z = -r./s - u.*y(y)
 * where all vectors except of y are of length 17 (length of y >= 17).
 */
void beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_17(beliefPenaltyMPC_FLOAT *r, beliefPenaltyMPC_FLOAT *s, beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *y, int* yidx, beliefPenaltyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<17; i++ ){
		z[i] = -r[i]/s[i] - u[i]*y[yidx[i]];
	}
}


/** 
 * Computes z = -r./s + u.*y(y)
 * where all vectors except of y are of length 7 (length of y >= 7).
 */
void beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_7(beliefPenaltyMPC_FLOAT *r, beliefPenaltyMPC_FLOAT *s, beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *y, int* yidx, beliefPenaltyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<7; i++ ){
		z[i] = -r[i]/s[i] + u[i]*y[yidx[i]];
	}
}


/** 
 * Computes z = -r./s - u.*y(y)
 * where all vectors except of y are of length 5 (length of y >= 5).
 */
void beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_5(beliefPenaltyMPC_FLOAT *r, beliefPenaltyMPC_FLOAT *s, beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *y, int* yidx, beliefPenaltyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<5; i++ ){
		z[i] = -r[i]/s[i] - u[i]*y[yidx[i]];
	}
}


/** 
 * Computes z = -r./s + u.*y(y)
 * where all vectors except of y are of length 5 (length of y >= 5).
 */
void beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(beliefPenaltyMPC_FLOAT *r, beliefPenaltyMPC_FLOAT *s, beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *y, int* yidx, beliefPenaltyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<5; i++ ){
		z[i] = -r[i]/s[i] + u[i]*y[yidx[i]];
	}
}


/*
 * Computes ds = -l.\(r + s.*dl) for vectors of length 226.
 */
void beliefPenaltyMPC_LA_VSUB7_226(beliefPenaltyMPC_FLOAT *l, beliefPenaltyMPC_FLOAT *r, beliefPenaltyMPC_FLOAT *s, beliefPenaltyMPC_FLOAT *dl, beliefPenaltyMPC_FLOAT *ds)
{
	int i;
	for( i=0; i<226; i++){
		ds[i] = -(r[i] + s[i]*dl[i])/l[i];
	}
}


/*
 * Vector addition x = x + y for vectors of length 158.
 */
void beliefPenaltyMPC_LA_VADD_158(beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *y)
{
	int i;
	for( i=0; i<158; i++){
		x[i] += y[i];
	}
}


/*
 * Vector addition x = x + y for vectors of length 50.
 */
void beliefPenaltyMPC_LA_VADD_50(beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *y)
{
	int i;
	for( i=0; i<50; i++){
		x[i] += y[i];
	}
}


/*
 * Vector addition x = x + y for vectors of length 226.
 */
void beliefPenaltyMPC_LA_VADD_226(beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *y)
{
	int i;
	for( i=0; i<226; i++){
		x[i] += y[i];
	}
}


/**
 * Backtracking line search for combined predictor/corrector step.
 * Update on variables with safety factor gamma (to keep us away from
 * boundary).
 */
int beliefPenaltyMPC_LINESEARCH_BACKTRACKING_COMBINED(beliefPenaltyMPC_FLOAT *z, beliefPenaltyMPC_FLOAT *v, beliefPenaltyMPC_FLOAT *l, beliefPenaltyMPC_FLOAT *s, beliefPenaltyMPC_FLOAT *dz, beliefPenaltyMPC_FLOAT *dv, beliefPenaltyMPC_FLOAT *dl, beliefPenaltyMPC_FLOAT *ds, beliefPenaltyMPC_FLOAT *a, beliefPenaltyMPC_FLOAT *mu)
{
    int i, lsIt=1;       
    beliefPenaltyMPC_FLOAT dltemp;
    beliefPenaltyMPC_FLOAT dstemp;    
    beliefPenaltyMPC_FLOAT a_gamma;
            
    *a = 1.0;
    while( 1 ){                        

        /* check whether search criterion is fulfilled */
        for( i=0; i<226; i++ ){
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
        if( i == 226 ){
            break;
        } else {
            *a *= beliefPenaltyMPC_SET_LS_SCALE;
            if( *a < beliefPenaltyMPC_SET_LS_MINSTEP ){
                return beliefPenaltyMPC_NOPROGRESS;
            }
        }
    }
    
    /* update variables with safety margin */
    a_gamma = (*a)*beliefPenaltyMPC_SET_LS_MAXSTEP;
    
    /* primal variables */
    for( i=0; i<158; i++ ){
        z[i] += a_gamma*dz[i];
    }
    
    /* equality constraint multipliers */
    for( i=0; i<50; i++ ){
        v[i] += a_gamma*dv[i];
    }
    
    /* inequality constraint multipliers & slacks, also update mu */
    *mu = 0;
    for( i=0; i<226; i++ ){
        dltemp = l[i] + a_gamma*dl[i]; l[i] = dltemp;
        dstemp = s[i] + a_gamma*ds[i]; s[i] = dstemp;
        *mu += dltemp*dstemp;
    }
    
    *a = a_gamma;
    *mu /= (beliefPenaltyMPC_FLOAT)226;
    return lsIt;
}




/* VARIABLE DEFINITIONS ------------------------------------------------ */
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_z[158];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_v[50];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_dz_aff[158];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_dv_aff[50];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_grad_cost[158];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_grad_eq[158];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_rd[158];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_l[226];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_s[226];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_lbys[226];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_dl_aff[226];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ds_aff[226];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_dz_cc[158];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_dv_cc[50];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_dl_cc[226];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ds_cc[226];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ccrhs[226];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_grad_ineq[158];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_H0[17] = {0.0000000000000000E+000, 0.0000000000000000E+000, 2.0000000000000000E+001, 2.0000000000000000E+001, 2.0000000000000000E+001, 2.0000000000000000E+000, 2.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_z0 = beliefPenaltyMPC_z + 0;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzaff0 = beliefPenaltyMPC_dz_aff + 0;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzcc0 = beliefPenaltyMPC_dz_cc + 0;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_rd0 = beliefPenaltyMPC_rd + 0;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lbyrd0[17];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_cost0 = beliefPenaltyMPC_grad_cost + 0;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_eq0 = beliefPenaltyMPC_grad_eq + 0;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_ineq0 = beliefPenaltyMPC_grad_ineq + 0;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ctv0[17];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_v0 = beliefPenaltyMPC_v + 0;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_re0[10];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_beta0[10];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_betacc0[10];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvaff0 = beliefPenaltyMPC_dv_aff + 0;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvcc0 = beliefPenaltyMPC_dv_cc + 0;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_V0[170];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Yd0[55];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ld0[55];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_yy0[10];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_bmy0[10];
int beliefPenaltyMPC_lbIdx0[17] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llb0 = beliefPenaltyMPC_l + 0;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_slb0 = beliefPenaltyMPC_s + 0;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llbbyslb0 = beliefPenaltyMPC_lbys + 0;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_rilb0[17];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbaff0 = beliefPenaltyMPC_dl_aff + 0;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbaff0 = beliefPenaltyMPC_ds_aff + 0;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbcc0 = beliefPenaltyMPC_dl_cc + 0;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbcc0 = beliefPenaltyMPC_ds_cc + 0;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsl0 = beliefPenaltyMPC_ccrhs + 0;
int beliefPenaltyMPC_ubIdx0[7] = {0, 1, 2, 3, 4, 5, 6};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lub0 = beliefPenaltyMPC_l + 17;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_sub0 = beliefPenaltyMPC_s + 17;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lubbysub0 = beliefPenaltyMPC_lbys + 17;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_riub0[7];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubaff0 = beliefPenaltyMPC_dl_aff + 17;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubaff0 = beliefPenaltyMPC_ds_aff + 17;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubcc0 = beliefPenaltyMPC_dl_cc + 17;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubcc0 = beliefPenaltyMPC_ds_cc + 17;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsub0 = beliefPenaltyMPC_ccrhs + 17;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Phi0[17];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_z1 = beliefPenaltyMPC_z + 17;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzaff1 = beliefPenaltyMPC_dz_aff + 17;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzcc1 = beliefPenaltyMPC_dz_cc + 17;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_rd1 = beliefPenaltyMPC_rd + 17;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lbyrd1[17];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_cost1 = beliefPenaltyMPC_grad_cost + 17;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_eq1 = beliefPenaltyMPC_grad_eq + 17;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_ineq1 = beliefPenaltyMPC_grad_ineq + 17;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ctv1[17];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_v1 = beliefPenaltyMPC_v + 10;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_re1[5];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_beta1[5];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_betacc1[5];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvaff1 = beliefPenaltyMPC_dv_aff + 10;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvcc1 = beliefPenaltyMPC_dv_cc + 10;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_V1[85];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Yd1[15];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ld1[15];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_yy1[5];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_bmy1[5];
int beliefPenaltyMPC_lbIdx1[17] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llb1 = beliefPenaltyMPC_l + 24;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_slb1 = beliefPenaltyMPC_s + 24;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llbbyslb1 = beliefPenaltyMPC_lbys + 24;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_rilb1[17];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbaff1 = beliefPenaltyMPC_dl_aff + 24;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbaff1 = beliefPenaltyMPC_ds_aff + 24;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbcc1 = beliefPenaltyMPC_dl_cc + 24;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbcc1 = beliefPenaltyMPC_ds_cc + 24;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsl1 = beliefPenaltyMPC_ccrhs + 24;
int beliefPenaltyMPC_ubIdx1[7] = {0, 1, 2, 3, 4, 5, 6};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lub1 = beliefPenaltyMPC_l + 41;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_sub1 = beliefPenaltyMPC_s + 41;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lubbysub1 = beliefPenaltyMPC_lbys + 41;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_riub1[7];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubaff1 = beliefPenaltyMPC_dl_aff + 41;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubaff1 = beliefPenaltyMPC_ds_aff + 41;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubcc1 = beliefPenaltyMPC_dl_cc + 41;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubcc1 = beliefPenaltyMPC_ds_cc + 41;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsub1 = beliefPenaltyMPC_ccrhs + 41;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Phi1[17];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_D1[170] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, -1.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 
0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, -1.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 
0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, -1.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 
0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, -1.0000000000000000E+000, 0.0000000000000000E+000, 
0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, -1.0000000000000000E+000, 
0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 
0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 
0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 
0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 
0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 
0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 
0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 
0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 
0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 
0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 
0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 
0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_W1[170];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ysd1[50];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lsd1[50];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_z2 = beliefPenaltyMPC_z + 34;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzaff2 = beliefPenaltyMPC_dz_aff + 34;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzcc2 = beliefPenaltyMPC_dz_cc + 34;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_rd2 = beliefPenaltyMPC_rd + 34;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lbyrd2[17];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_cost2 = beliefPenaltyMPC_grad_cost + 34;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_eq2 = beliefPenaltyMPC_grad_eq + 34;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_ineq2 = beliefPenaltyMPC_grad_ineq + 34;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ctv2[17];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_v2 = beliefPenaltyMPC_v + 15;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_re2[5];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_beta2[5];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_betacc2[5];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvaff2 = beliefPenaltyMPC_dv_aff + 15;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvcc2 = beliefPenaltyMPC_dv_cc + 15;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_V2[85];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Yd2[15];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ld2[15];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_yy2[5];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_bmy2[5];
int beliefPenaltyMPC_lbIdx2[17] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llb2 = beliefPenaltyMPC_l + 48;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_slb2 = beliefPenaltyMPC_s + 48;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llbbyslb2 = beliefPenaltyMPC_lbys + 48;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_rilb2[17];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbaff2 = beliefPenaltyMPC_dl_aff + 48;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbaff2 = beliefPenaltyMPC_ds_aff + 48;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbcc2 = beliefPenaltyMPC_dl_cc + 48;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbcc2 = beliefPenaltyMPC_ds_cc + 48;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsl2 = beliefPenaltyMPC_ccrhs + 48;
int beliefPenaltyMPC_ubIdx2[7] = {0, 1, 2, 3, 4, 5, 6};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lub2 = beliefPenaltyMPC_l + 65;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_sub2 = beliefPenaltyMPC_s + 65;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lubbysub2 = beliefPenaltyMPC_lbys + 65;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_riub2[7];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubaff2 = beliefPenaltyMPC_dl_aff + 65;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubaff2 = beliefPenaltyMPC_ds_aff + 65;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubcc2 = beliefPenaltyMPC_dl_cc + 65;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubcc2 = beliefPenaltyMPC_ds_cc + 65;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsub2 = beliefPenaltyMPC_ccrhs + 65;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Phi2[17];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_D2[17] = {-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000};
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_W2[17];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ysd2[25];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lsd2[25];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_z3 = beliefPenaltyMPC_z + 51;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzaff3 = beliefPenaltyMPC_dz_aff + 51;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzcc3 = beliefPenaltyMPC_dz_cc + 51;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_rd3 = beliefPenaltyMPC_rd + 51;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lbyrd3[17];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_cost3 = beliefPenaltyMPC_grad_cost + 51;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_eq3 = beliefPenaltyMPC_grad_eq + 51;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_ineq3 = beliefPenaltyMPC_grad_ineq + 51;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ctv3[17];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_v3 = beliefPenaltyMPC_v + 20;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_re3[5];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_beta3[5];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_betacc3[5];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvaff3 = beliefPenaltyMPC_dv_aff + 20;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvcc3 = beliefPenaltyMPC_dv_cc + 20;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_V3[85];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Yd3[15];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ld3[15];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_yy3[5];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_bmy3[5];
int beliefPenaltyMPC_lbIdx3[17] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llb3 = beliefPenaltyMPC_l + 72;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_slb3 = beliefPenaltyMPC_s + 72;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llbbyslb3 = beliefPenaltyMPC_lbys + 72;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_rilb3[17];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbaff3 = beliefPenaltyMPC_dl_aff + 72;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbaff3 = beliefPenaltyMPC_ds_aff + 72;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbcc3 = beliefPenaltyMPC_dl_cc + 72;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbcc3 = beliefPenaltyMPC_ds_cc + 72;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsl3 = beliefPenaltyMPC_ccrhs + 72;
int beliefPenaltyMPC_ubIdx3[7] = {0, 1, 2, 3, 4, 5, 6};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lub3 = beliefPenaltyMPC_l + 89;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_sub3 = beliefPenaltyMPC_s + 89;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lubbysub3 = beliefPenaltyMPC_lbys + 89;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_riub3[7];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubaff3 = beliefPenaltyMPC_dl_aff + 89;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubaff3 = beliefPenaltyMPC_ds_aff + 89;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubcc3 = beliefPenaltyMPC_dl_cc + 89;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubcc3 = beliefPenaltyMPC_ds_cc + 89;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsub3 = beliefPenaltyMPC_ccrhs + 89;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Phi3[17];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_W3[17];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ysd3[25];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lsd3[25];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_z4 = beliefPenaltyMPC_z + 68;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzaff4 = beliefPenaltyMPC_dz_aff + 68;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzcc4 = beliefPenaltyMPC_dz_cc + 68;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_rd4 = beliefPenaltyMPC_rd + 68;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lbyrd4[17];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_cost4 = beliefPenaltyMPC_grad_cost + 68;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_eq4 = beliefPenaltyMPC_grad_eq + 68;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_ineq4 = beliefPenaltyMPC_grad_ineq + 68;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ctv4[17];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_v4 = beliefPenaltyMPC_v + 25;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_re4[5];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_beta4[5];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_betacc4[5];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvaff4 = beliefPenaltyMPC_dv_aff + 25;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvcc4 = beliefPenaltyMPC_dv_cc + 25;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_V4[85];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Yd4[15];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ld4[15];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_yy4[5];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_bmy4[5];
int beliefPenaltyMPC_lbIdx4[17] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llb4 = beliefPenaltyMPC_l + 96;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_slb4 = beliefPenaltyMPC_s + 96;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llbbyslb4 = beliefPenaltyMPC_lbys + 96;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_rilb4[17];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbaff4 = beliefPenaltyMPC_dl_aff + 96;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbaff4 = beliefPenaltyMPC_ds_aff + 96;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbcc4 = beliefPenaltyMPC_dl_cc + 96;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbcc4 = beliefPenaltyMPC_ds_cc + 96;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsl4 = beliefPenaltyMPC_ccrhs + 96;
int beliefPenaltyMPC_ubIdx4[7] = {0, 1, 2, 3, 4, 5, 6};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lub4 = beliefPenaltyMPC_l + 113;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_sub4 = beliefPenaltyMPC_s + 113;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lubbysub4 = beliefPenaltyMPC_lbys + 113;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_riub4[7];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubaff4 = beliefPenaltyMPC_dl_aff + 113;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubaff4 = beliefPenaltyMPC_ds_aff + 113;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubcc4 = beliefPenaltyMPC_dl_cc + 113;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubcc4 = beliefPenaltyMPC_ds_cc + 113;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsub4 = beliefPenaltyMPC_ccrhs + 113;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Phi4[17];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_W4[17];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ysd4[25];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lsd4[25];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_z5 = beliefPenaltyMPC_z + 85;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzaff5 = beliefPenaltyMPC_dz_aff + 85;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzcc5 = beliefPenaltyMPC_dz_cc + 85;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_rd5 = beliefPenaltyMPC_rd + 85;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lbyrd5[17];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_cost5 = beliefPenaltyMPC_grad_cost + 85;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_eq5 = beliefPenaltyMPC_grad_eq + 85;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_ineq5 = beliefPenaltyMPC_grad_ineq + 85;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ctv5[17];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_v5 = beliefPenaltyMPC_v + 30;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_re5[5];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_beta5[5];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_betacc5[5];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvaff5 = beliefPenaltyMPC_dv_aff + 30;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvcc5 = beliefPenaltyMPC_dv_cc + 30;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_V5[85];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Yd5[15];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ld5[15];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_yy5[5];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_bmy5[5];
int beliefPenaltyMPC_lbIdx5[17] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llb5 = beliefPenaltyMPC_l + 120;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_slb5 = beliefPenaltyMPC_s + 120;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llbbyslb5 = beliefPenaltyMPC_lbys + 120;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_rilb5[17];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbaff5 = beliefPenaltyMPC_dl_aff + 120;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbaff5 = beliefPenaltyMPC_ds_aff + 120;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbcc5 = beliefPenaltyMPC_dl_cc + 120;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbcc5 = beliefPenaltyMPC_ds_cc + 120;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsl5 = beliefPenaltyMPC_ccrhs + 120;
int beliefPenaltyMPC_ubIdx5[7] = {0, 1, 2, 3, 4, 5, 6};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lub5 = beliefPenaltyMPC_l + 137;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_sub5 = beliefPenaltyMPC_s + 137;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lubbysub5 = beliefPenaltyMPC_lbys + 137;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_riub5[7];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubaff5 = beliefPenaltyMPC_dl_aff + 137;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubaff5 = beliefPenaltyMPC_ds_aff + 137;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubcc5 = beliefPenaltyMPC_dl_cc + 137;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubcc5 = beliefPenaltyMPC_ds_cc + 137;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsub5 = beliefPenaltyMPC_ccrhs + 137;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Phi5[17];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_W5[17];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ysd5[25];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lsd5[25];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_z6 = beliefPenaltyMPC_z + 102;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzaff6 = beliefPenaltyMPC_dz_aff + 102;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzcc6 = beliefPenaltyMPC_dz_cc + 102;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_rd6 = beliefPenaltyMPC_rd + 102;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lbyrd6[17];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_cost6 = beliefPenaltyMPC_grad_cost + 102;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_eq6 = beliefPenaltyMPC_grad_eq + 102;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_ineq6 = beliefPenaltyMPC_grad_ineq + 102;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ctv6[17];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_v6 = beliefPenaltyMPC_v + 35;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_re6[5];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_beta6[5];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_betacc6[5];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvaff6 = beliefPenaltyMPC_dv_aff + 35;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvcc6 = beliefPenaltyMPC_dv_cc + 35;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_V6[85];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Yd6[15];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ld6[15];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_yy6[5];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_bmy6[5];
int beliefPenaltyMPC_lbIdx6[17] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llb6 = beliefPenaltyMPC_l + 144;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_slb6 = beliefPenaltyMPC_s + 144;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llbbyslb6 = beliefPenaltyMPC_lbys + 144;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_rilb6[17];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbaff6 = beliefPenaltyMPC_dl_aff + 144;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbaff6 = beliefPenaltyMPC_ds_aff + 144;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbcc6 = beliefPenaltyMPC_dl_cc + 144;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbcc6 = beliefPenaltyMPC_ds_cc + 144;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsl6 = beliefPenaltyMPC_ccrhs + 144;
int beliefPenaltyMPC_ubIdx6[7] = {0, 1, 2, 3, 4, 5, 6};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lub6 = beliefPenaltyMPC_l + 161;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_sub6 = beliefPenaltyMPC_s + 161;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lubbysub6 = beliefPenaltyMPC_lbys + 161;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_riub6[7];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubaff6 = beliefPenaltyMPC_dl_aff + 161;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubaff6 = beliefPenaltyMPC_ds_aff + 161;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubcc6 = beliefPenaltyMPC_dl_cc + 161;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubcc6 = beliefPenaltyMPC_ds_cc + 161;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsub6 = beliefPenaltyMPC_ccrhs + 161;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Phi6[17];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_W6[17];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ysd6[25];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lsd6[25];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_z7 = beliefPenaltyMPC_z + 119;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzaff7 = beliefPenaltyMPC_dz_aff + 119;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzcc7 = beliefPenaltyMPC_dz_cc + 119;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_rd7 = beliefPenaltyMPC_rd + 119;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lbyrd7[17];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_cost7 = beliefPenaltyMPC_grad_cost + 119;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_eq7 = beliefPenaltyMPC_grad_eq + 119;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_ineq7 = beliefPenaltyMPC_grad_ineq + 119;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ctv7[17];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_v7 = beliefPenaltyMPC_v + 40;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_re7[5];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_beta7[5];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_betacc7[5];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvaff7 = beliefPenaltyMPC_dv_aff + 40;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvcc7 = beliefPenaltyMPC_dv_cc + 40;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_V7[85];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Yd7[15];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ld7[15];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_yy7[5];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_bmy7[5];
int beliefPenaltyMPC_lbIdx7[17] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llb7 = beliefPenaltyMPC_l + 168;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_slb7 = beliefPenaltyMPC_s + 168;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llbbyslb7 = beliefPenaltyMPC_lbys + 168;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_rilb7[17];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbaff7 = beliefPenaltyMPC_dl_aff + 168;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbaff7 = beliefPenaltyMPC_ds_aff + 168;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbcc7 = beliefPenaltyMPC_dl_cc + 168;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbcc7 = beliefPenaltyMPC_ds_cc + 168;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsl7 = beliefPenaltyMPC_ccrhs + 168;
int beliefPenaltyMPC_ubIdx7[7] = {0, 1, 2, 3, 4, 5, 6};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lub7 = beliefPenaltyMPC_l + 185;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_sub7 = beliefPenaltyMPC_s + 185;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lubbysub7 = beliefPenaltyMPC_lbys + 185;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_riub7[7];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubaff7 = beliefPenaltyMPC_dl_aff + 185;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubaff7 = beliefPenaltyMPC_ds_aff + 185;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubcc7 = beliefPenaltyMPC_dl_cc + 185;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubcc7 = beliefPenaltyMPC_ds_cc + 185;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsub7 = beliefPenaltyMPC_ccrhs + 185;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Phi7[17];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_W7[17];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ysd7[25];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lsd7[25];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_z8 = beliefPenaltyMPC_z + 136;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzaff8 = beliefPenaltyMPC_dz_aff + 136;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzcc8 = beliefPenaltyMPC_dz_cc + 136;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_rd8 = beliefPenaltyMPC_rd + 136;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lbyrd8[17];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_cost8 = beliefPenaltyMPC_grad_cost + 136;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_eq8 = beliefPenaltyMPC_grad_eq + 136;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_ineq8 = beliefPenaltyMPC_grad_ineq + 136;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ctv8[17];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_v8 = beliefPenaltyMPC_v + 45;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_re8[5];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_beta8[5];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_betacc8[5];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvaff8 = beliefPenaltyMPC_dv_aff + 45;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvcc8 = beliefPenaltyMPC_dv_cc + 45;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_V8[85];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Yd8[15];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ld8[15];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_yy8[5];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_bmy8[5];
int beliefPenaltyMPC_lbIdx8[17] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llb8 = beliefPenaltyMPC_l + 192;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_slb8 = beliefPenaltyMPC_s + 192;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llbbyslb8 = beliefPenaltyMPC_lbys + 192;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_rilb8[17];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbaff8 = beliefPenaltyMPC_dl_aff + 192;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbaff8 = beliefPenaltyMPC_ds_aff + 192;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbcc8 = beliefPenaltyMPC_dl_cc + 192;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbcc8 = beliefPenaltyMPC_ds_cc + 192;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsl8 = beliefPenaltyMPC_ccrhs + 192;
int beliefPenaltyMPC_ubIdx8[7] = {0, 1, 2, 3, 4, 5, 6};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lub8 = beliefPenaltyMPC_l + 209;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_sub8 = beliefPenaltyMPC_s + 209;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lubbysub8 = beliefPenaltyMPC_lbys + 209;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_riub8[7];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubaff8 = beliefPenaltyMPC_dl_aff + 209;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubaff8 = beliefPenaltyMPC_ds_aff + 209;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubcc8 = beliefPenaltyMPC_dl_cc + 209;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubcc8 = beliefPenaltyMPC_ds_cc + 209;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsub8 = beliefPenaltyMPC_ccrhs + 209;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Phi8[17];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_W8[17];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ysd8[25];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lsd8[25];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_H9[5] = {0.0000000000000000E+000, 0.0000000000000000E+000, 2.0000000000000000E+001, 2.0000000000000000E+001, 2.0000000000000000E+001};
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_f9[5] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_z9 = beliefPenaltyMPC_z + 153;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzaff9 = beliefPenaltyMPC_dz_aff + 153;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzcc9 = beliefPenaltyMPC_dz_cc + 153;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_rd9 = beliefPenaltyMPC_rd + 153;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lbyrd9[5];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_cost9 = beliefPenaltyMPC_grad_cost + 153;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_eq9 = beliefPenaltyMPC_grad_eq + 153;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_ineq9 = beliefPenaltyMPC_grad_ineq + 153;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ctv9[5];
int beliefPenaltyMPC_lbIdx9[5] = {0, 1, 2, 3, 4};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llb9 = beliefPenaltyMPC_l + 216;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_slb9 = beliefPenaltyMPC_s + 216;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llbbyslb9 = beliefPenaltyMPC_lbys + 216;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_rilb9[5];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbaff9 = beliefPenaltyMPC_dl_aff + 216;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbaff9 = beliefPenaltyMPC_ds_aff + 216;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbcc9 = beliefPenaltyMPC_dl_cc + 216;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbcc9 = beliefPenaltyMPC_ds_cc + 216;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsl9 = beliefPenaltyMPC_ccrhs + 216;
int beliefPenaltyMPC_ubIdx9[5] = {0, 1, 2, 3, 4};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lub9 = beliefPenaltyMPC_l + 221;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_sub9 = beliefPenaltyMPC_s + 221;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lubbysub9 = beliefPenaltyMPC_lbys + 221;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_riub9[5];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubaff9 = beliefPenaltyMPC_dl_aff + 221;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubaff9 = beliefPenaltyMPC_ds_aff + 221;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubcc9 = beliefPenaltyMPC_dl_cc + 221;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubcc9 = beliefPenaltyMPC_ds_cc + 221;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsub9 = beliefPenaltyMPC_ccrhs + 221;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Phi9[5];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_D9[5] = {-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000};
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_W9[5];
beliefPenaltyMPC_FLOAT musigma;
beliefPenaltyMPC_FLOAT sigma_3rdroot;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Diag1_0[17];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Diag2_0[17];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_L_0[136];




/* SOLVER CODE --------------------------------------------------------- */
int beliefPenaltyMPC_solve(beliefPenaltyMPC_params* params, beliefPenaltyMPC_output* output, beliefPenaltyMPC_info* info)
{	
int exitcode;

#if beliefPenaltyMPC_SET_TIMING == 1
	beliefPenaltyMPC_timer solvertimer;
	beliefPenaltyMPC_tic(&solvertimer);
#endif
/* FUNCTION CALLS INTO LA LIBRARY -------------------------------------- */
info->it = 0;
beliefPenaltyMPC_LA_INITIALIZEVECTOR_158(beliefPenaltyMPC_z, 0);
beliefPenaltyMPC_LA_INITIALIZEVECTOR_50(beliefPenaltyMPC_v, 1);
beliefPenaltyMPC_LA_INITIALIZEVECTOR_226(beliefPenaltyMPC_l, 1);
beliefPenaltyMPC_LA_INITIALIZEVECTOR_226(beliefPenaltyMPC_s, 1);
info->mu = 0;
beliefPenaltyMPC_LA_DOTACC_226(beliefPenaltyMPC_l, beliefPenaltyMPC_s, &info->mu);
info->mu /= 226;
while( 1 ){
info->pobj = 0;
beliefPenaltyMPC_LA_DIAG_QUADFCN_17(beliefPenaltyMPC_H0, params->f1, beliefPenaltyMPC_z0, beliefPenaltyMPC_grad_cost0, &info->pobj);
beliefPenaltyMPC_LA_DIAG_QUADFCN_17(beliefPenaltyMPC_H0, params->f2, beliefPenaltyMPC_z1, beliefPenaltyMPC_grad_cost1, &info->pobj);
beliefPenaltyMPC_LA_DIAG_QUADFCN_17(beliefPenaltyMPC_H0, params->f3, beliefPenaltyMPC_z2, beliefPenaltyMPC_grad_cost2, &info->pobj);
beliefPenaltyMPC_LA_DIAG_QUADFCN_17(beliefPenaltyMPC_H0, params->f4, beliefPenaltyMPC_z3, beliefPenaltyMPC_grad_cost3, &info->pobj);
beliefPenaltyMPC_LA_DIAG_QUADFCN_17(beliefPenaltyMPC_H0, params->f5, beliefPenaltyMPC_z4, beliefPenaltyMPC_grad_cost4, &info->pobj);
beliefPenaltyMPC_LA_DIAG_QUADFCN_17(beliefPenaltyMPC_H0, params->f6, beliefPenaltyMPC_z5, beliefPenaltyMPC_grad_cost5, &info->pobj);
beliefPenaltyMPC_LA_DIAG_QUADFCN_17(beliefPenaltyMPC_H0, params->f7, beliefPenaltyMPC_z6, beliefPenaltyMPC_grad_cost6, &info->pobj);
beliefPenaltyMPC_LA_DIAG_QUADFCN_17(beliefPenaltyMPC_H0, params->f8, beliefPenaltyMPC_z7, beliefPenaltyMPC_grad_cost7, &info->pobj);
beliefPenaltyMPC_LA_DIAG_QUADFCN_17(beliefPenaltyMPC_H0, params->f9, beliefPenaltyMPC_z8, beliefPenaltyMPC_grad_cost8, &info->pobj);
beliefPenaltyMPC_LA_DIAG_QUADFCN_5(beliefPenaltyMPC_H9, beliefPenaltyMPC_f9, beliefPenaltyMPC_z9, beliefPenaltyMPC_grad_cost9, &info->pobj);
info->res_eq = 0;
info->dgap = 0;
beliefPenaltyMPC_LA_DENSE_MVMSUB3_10_17_17(params->C1, beliefPenaltyMPC_z0, beliefPenaltyMPC_D1, beliefPenaltyMPC_z1, params->e1, beliefPenaltyMPC_v0, beliefPenaltyMPC_re0, &info->dgap, &info->res_eq);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_5_17_17(params->C2, beliefPenaltyMPC_z1, beliefPenaltyMPC_D2, beliefPenaltyMPC_z2, params->e2, beliefPenaltyMPC_v1, beliefPenaltyMPC_re1, &info->dgap, &info->res_eq);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_5_17_17(params->C3, beliefPenaltyMPC_z2, beliefPenaltyMPC_D2, beliefPenaltyMPC_z3, params->e3, beliefPenaltyMPC_v2, beliefPenaltyMPC_re2, &info->dgap, &info->res_eq);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_5_17_17(params->C4, beliefPenaltyMPC_z3, beliefPenaltyMPC_D2, beliefPenaltyMPC_z4, params->e4, beliefPenaltyMPC_v3, beliefPenaltyMPC_re3, &info->dgap, &info->res_eq);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_5_17_17(params->C5, beliefPenaltyMPC_z4, beliefPenaltyMPC_D2, beliefPenaltyMPC_z5, params->e5, beliefPenaltyMPC_v4, beliefPenaltyMPC_re4, &info->dgap, &info->res_eq);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_5_17_17(params->C6, beliefPenaltyMPC_z5, beliefPenaltyMPC_D2, beliefPenaltyMPC_z6, params->e6, beliefPenaltyMPC_v5, beliefPenaltyMPC_re5, &info->dgap, &info->res_eq);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_5_17_17(params->C7, beliefPenaltyMPC_z6, beliefPenaltyMPC_D2, beliefPenaltyMPC_z7, params->e7, beliefPenaltyMPC_v6, beliefPenaltyMPC_re6, &info->dgap, &info->res_eq);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_5_17_17(params->C8, beliefPenaltyMPC_z7, beliefPenaltyMPC_D2, beliefPenaltyMPC_z8, params->e8, beliefPenaltyMPC_v7, beliefPenaltyMPC_re7, &info->dgap, &info->res_eq);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_5_17_5(params->C9, beliefPenaltyMPC_z8, beliefPenaltyMPC_D9, beliefPenaltyMPC_z9, params->e9, beliefPenaltyMPC_v8, beliefPenaltyMPC_re8, &info->dgap, &info->res_eq);
beliefPenaltyMPC_LA_DENSE_MTVM_10_17(params->C1, beliefPenaltyMPC_v0, beliefPenaltyMPC_grad_eq0);
beliefPenaltyMPC_LA_DENSE_MTVM2_5_17_10(params->C2, beliefPenaltyMPC_v1, beliefPenaltyMPC_D1, beliefPenaltyMPC_v0, beliefPenaltyMPC_grad_eq1);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_5_17_5(params->C3, beliefPenaltyMPC_v2, beliefPenaltyMPC_D2, beliefPenaltyMPC_v1, beliefPenaltyMPC_grad_eq2);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_5_17_5(params->C4, beliefPenaltyMPC_v3, beliefPenaltyMPC_D2, beliefPenaltyMPC_v2, beliefPenaltyMPC_grad_eq3);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_5_17_5(params->C5, beliefPenaltyMPC_v4, beliefPenaltyMPC_D2, beliefPenaltyMPC_v3, beliefPenaltyMPC_grad_eq4);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_5_17_5(params->C6, beliefPenaltyMPC_v5, beliefPenaltyMPC_D2, beliefPenaltyMPC_v4, beliefPenaltyMPC_grad_eq5);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_5_17_5(params->C7, beliefPenaltyMPC_v6, beliefPenaltyMPC_D2, beliefPenaltyMPC_v5, beliefPenaltyMPC_grad_eq6);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_5_17_5(params->C8, beliefPenaltyMPC_v7, beliefPenaltyMPC_D2, beliefPenaltyMPC_v6, beliefPenaltyMPC_grad_eq7);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_5_17_5(params->C9, beliefPenaltyMPC_v8, beliefPenaltyMPC_D2, beliefPenaltyMPC_v7, beliefPenaltyMPC_grad_eq8);
beliefPenaltyMPC_LA_DIAGZERO_MTVM_5_5(beliefPenaltyMPC_D9, beliefPenaltyMPC_v8, beliefPenaltyMPC_grad_eq9);
info->res_ineq = 0;
beliefPenaltyMPC_LA_VSUBADD3_17(params->lb1, beliefPenaltyMPC_z0, beliefPenaltyMPC_lbIdx0, beliefPenaltyMPC_llb0, beliefPenaltyMPC_slb0, beliefPenaltyMPC_rilb0, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD2_7(beliefPenaltyMPC_z0, beliefPenaltyMPC_ubIdx0, params->ub1, beliefPenaltyMPC_lub0, beliefPenaltyMPC_sub0, beliefPenaltyMPC_riub0, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD3_17(params->lb2, beliefPenaltyMPC_z1, beliefPenaltyMPC_lbIdx1, beliefPenaltyMPC_llb1, beliefPenaltyMPC_slb1, beliefPenaltyMPC_rilb1, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD2_7(beliefPenaltyMPC_z1, beliefPenaltyMPC_ubIdx1, params->ub2, beliefPenaltyMPC_lub1, beliefPenaltyMPC_sub1, beliefPenaltyMPC_riub1, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD3_17(params->lb3, beliefPenaltyMPC_z2, beliefPenaltyMPC_lbIdx2, beliefPenaltyMPC_llb2, beliefPenaltyMPC_slb2, beliefPenaltyMPC_rilb2, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD2_7(beliefPenaltyMPC_z2, beliefPenaltyMPC_ubIdx2, params->ub3, beliefPenaltyMPC_lub2, beliefPenaltyMPC_sub2, beliefPenaltyMPC_riub2, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD3_17(params->lb4, beliefPenaltyMPC_z3, beliefPenaltyMPC_lbIdx3, beliefPenaltyMPC_llb3, beliefPenaltyMPC_slb3, beliefPenaltyMPC_rilb3, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD2_7(beliefPenaltyMPC_z3, beliefPenaltyMPC_ubIdx3, params->ub4, beliefPenaltyMPC_lub3, beliefPenaltyMPC_sub3, beliefPenaltyMPC_riub3, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD3_17(params->lb5, beliefPenaltyMPC_z4, beliefPenaltyMPC_lbIdx4, beliefPenaltyMPC_llb4, beliefPenaltyMPC_slb4, beliefPenaltyMPC_rilb4, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD2_7(beliefPenaltyMPC_z4, beliefPenaltyMPC_ubIdx4, params->ub5, beliefPenaltyMPC_lub4, beliefPenaltyMPC_sub4, beliefPenaltyMPC_riub4, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD3_17(params->lb6, beliefPenaltyMPC_z5, beliefPenaltyMPC_lbIdx5, beliefPenaltyMPC_llb5, beliefPenaltyMPC_slb5, beliefPenaltyMPC_rilb5, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD2_7(beliefPenaltyMPC_z5, beliefPenaltyMPC_ubIdx5, params->ub6, beliefPenaltyMPC_lub5, beliefPenaltyMPC_sub5, beliefPenaltyMPC_riub5, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD3_17(params->lb7, beliefPenaltyMPC_z6, beliefPenaltyMPC_lbIdx6, beliefPenaltyMPC_llb6, beliefPenaltyMPC_slb6, beliefPenaltyMPC_rilb6, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD2_7(beliefPenaltyMPC_z6, beliefPenaltyMPC_ubIdx6, params->ub7, beliefPenaltyMPC_lub6, beliefPenaltyMPC_sub6, beliefPenaltyMPC_riub6, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD3_17(params->lb8, beliefPenaltyMPC_z7, beliefPenaltyMPC_lbIdx7, beliefPenaltyMPC_llb7, beliefPenaltyMPC_slb7, beliefPenaltyMPC_rilb7, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD2_7(beliefPenaltyMPC_z7, beliefPenaltyMPC_ubIdx7, params->ub8, beliefPenaltyMPC_lub7, beliefPenaltyMPC_sub7, beliefPenaltyMPC_riub7, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD3_17(params->lb9, beliefPenaltyMPC_z8, beliefPenaltyMPC_lbIdx8, beliefPenaltyMPC_llb8, beliefPenaltyMPC_slb8, beliefPenaltyMPC_rilb8, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD2_7(beliefPenaltyMPC_z8, beliefPenaltyMPC_ubIdx8, params->ub9, beliefPenaltyMPC_lub8, beliefPenaltyMPC_sub8, beliefPenaltyMPC_riub8, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD3_5(params->lb10, beliefPenaltyMPC_z9, beliefPenaltyMPC_lbIdx9, beliefPenaltyMPC_llb9, beliefPenaltyMPC_slb9, beliefPenaltyMPC_rilb9, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD2_5(beliefPenaltyMPC_z9, beliefPenaltyMPC_ubIdx9, params->ub10, beliefPenaltyMPC_lub9, beliefPenaltyMPC_sub9, beliefPenaltyMPC_riub9, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_INEQ_B_GRAD_17_17_7(beliefPenaltyMPC_lub0, beliefPenaltyMPC_sub0, beliefPenaltyMPC_riub0, beliefPenaltyMPC_llb0, beliefPenaltyMPC_slb0, beliefPenaltyMPC_rilb0, beliefPenaltyMPC_lbIdx0, beliefPenaltyMPC_ubIdx0, beliefPenaltyMPC_grad_ineq0, beliefPenaltyMPC_lubbysub0, beliefPenaltyMPC_llbbyslb0);
beliefPenaltyMPC_LA_INEQ_B_GRAD_17_17_7(beliefPenaltyMPC_lub1, beliefPenaltyMPC_sub1, beliefPenaltyMPC_riub1, beliefPenaltyMPC_llb1, beliefPenaltyMPC_slb1, beliefPenaltyMPC_rilb1, beliefPenaltyMPC_lbIdx1, beliefPenaltyMPC_ubIdx1, beliefPenaltyMPC_grad_ineq1, beliefPenaltyMPC_lubbysub1, beliefPenaltyMPC_llbbyslb1);
beliefPenaltyMPC_LA_INEQ_B_GRAD_17_17_7(beliefPenaltyMPC_lub2, beliefPenaltyMPC_sub2, beliefPenaltyMPC_riub2, beliefPenaltyMPC_llb2, beliefPenaltyMPC_slb2, beliefPenaltyMPC_rilb2, beliefPenaltyMPC_lbIdx2, beliefPenaltyMPC_ubIdx2, beliefPenaltyMPC_grad_ineq2, beliefPenaltyMPC_lubbysub2, beliefPenaltyMPC_llbbyslb2);
beliefPenaltyMPC_LA_INEQ_B_GRAD_17_17_7(beliefPenaltyMPC_lub3, beliefPenaltyMPC_sub3, beliefPenaltyMPC_riub3, beliefPenaltyMPC_llb3, beliefPenaltyMPC_slb3, beliefPenaltyMPC_rilb3, beliefPenaltyMPC_lbIdx3, beliefPenaltyMPC_ubIdx3, beliefPenaltyMPC_grad_ineq3, beliefPenaltyMPC_lubbysub3, beliefPenaltyMPC_llbbyslb3);
beliefPenaltyMPC_LA_INEQ_B_GRAD_17_17_7(beliefPenaltyMPC_lub4, beliefPenaltyMPC_sub4, beliefPenaltyMPC_riub4, beliefPenaltyMPC_llb4, beliefPenaltyMPC_slb4, beliefPenaltyMPC_rilb4, beliefPenaltyMPC_lbIdx4, beliefPenaltyMPC_ubIdx4, beliefPenaltyMPC_grad_ineq4, beliefPenaltyMPC_lubbysub4, beliefPenaltyMPC_llbbyslb4);
beliefPenaltyMPC_LA_INEQ_B_GRAD_17_17_7(beliefPenaltyMPC_lub5, beliefPenaltyMPC_sub5, beliefPenaltyMPC_riub5, beliefPenaltyMPC_llb5, beliefPenaltyMPC_slb5, beliefPenaltyMPC_rilb5, beliefPenaltyMPC_lbIdx5, beliefPenaltyMPC_ubIdx5, beliefPenaltyMPC_grad_ineq5, beliefPenaltyMPC_lubbysub5, beliefPenaltyMPC_llbbyslb5);
beliefPenaltyMPC_LA_INEQ_B_GRAD_17_17_7(beliefPenaltyMPC_lub6, beliefPenaltyMPC_sub6, beliefPenaltyMPC_riub6, beliefPenaltyMPC_llb6, beliefPenaltyMPC_slb6, beliefPenaltyMPC_rilb6, beliefPenaltyMPC_lbIdx6, beliefPenaltyMPC_ubIdx6, beliefPenaltyMPC_grad_ineq6, beliefPenaltyMPC_lubbysub6, beliefPenaltyMPC_llbbyslb6);
beliefPenaltyMPC_LA_INEQ_B_GRAD_17_17_7(beliefPenaltyMPC_lub7, beliefPenaltyMPC_sub7, beliefPenaltyMPC_riub7, beliefPenaltyMPC_llb7, beliefPenaltyMPC_slb7, beliefPenaltyMPC_rilb7, beliefPenaltyMPC_lbIdx7, beliefPenaltyMPC_ubIdx7, beliefPenaltyMPC_grad_ineq7, beliefPenaltyMPC_lubbysub7, beliefPenaltyMPC_llbbyslb7);
beliefPenaltyMPC_LA_INEQ_B_GRAD_17_17_7(beliefPenaltyMPC_lub8, beliefPenaltyMPC_sub8, beliefPenaltyMPC_riub8, beliefPenaltyMPC_llb8, beliefPenaltyMPC_slb8, beliefPenaltyMPC_rilb8, beliefPenaltyMPC_lbIdx8, beliefPenaltyMPC_ubIdx8, beliefPenaltyMPC_grad_ineq8, beliefPenaltyMPC_lubbysub8, beliefPenaltyMPC_llbbyslb8);
beliefPenaltyMPC_LA_INEQ_B_GRAD_5_5_5(beliefPenaltyMPC_lub9, beliefPenaltyMPC_sub9, beliefPenaltyMPC_riub9, beliefPenaltyMPC_llb9, beliefPenaltyMPC_slb9, beliefPenaltyMPC_rilb9, beliefPenaltyMPC_lbIdx9, beliefPenaltyMPC_ubIdx9, beliefPenaltyMPC_grad_ineq9, beliefPenaltyMPC_lubbysub9, beliefPenaltyMPC_llbbyslb9);
info->dobj = info->pobj - info->dgap;
info->rdgap = info->pobj ? info->dgap / info->pobj : 1e6;
if( info->rdgap < 0 ) info->rdgap = -info->rdgap;
if( info->mu < beliefPenaltyMPC_SET_ACC_KKTCOMPL
    && (info->rdgap < beliefPenaltyMPC_SET_ACC_RDGAP || info->dgap < beliefPenaltyMPC_SET_ACC_KKTCOMPL)
    && info->res_eq < beliefPenaltyMPC_SET_ACC_RESEQ
    && info->res_ineq < beliefPenaltyMPC_SET_ACC_RESINEQ ){
exitcode = beliefPenaltyMPC_OPTIMAL; break; }
if( info->it == beliefPenaltyMPC_SET_MAXIT ){
exitcode = beliefPenaltyMPC_MAXITREACHED; break; }
beliefPenaltyMPC_LA_VVADD3_158(beliefPenaltyMPC_grad_cost, beliefPenaltyMPC_grad_eq, beliefPenaltyMPC_grad_ineq, beliefPenaltyMPC_rd);
beliefPenaltyMPC_LA_DIAG_CHOL_LBUB_17_17_7(beliefPenaltyMPC_H0, beliefPenaltyMPC_llbbyslb0, beliefPenaltyMPC_lbIdx0, beliefPenaltyMPC_lubbysub0, beliefPenaltyMPC_ubIdx0, beliefPenaltyMPC_Phi0);
beliefPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_10_17(beliefPenaltyMPC_Phi0, params->C1, beliefPenaltyMPC_V0);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_17(beliefPenaltyMPC_Phi0, beliefPenaltyMPC_rd0, beliefPenaltyMPC_Lbyrd0);
beliefPenaltyMPC_LA_DIAG_CHOL_LBUB_17_17_7(beliefPenaltyMPC_H0, beliefPenaltyMPC_llbbyslb1, beliefPenaltyMPC_lbIdx1, beliefPenaltyMPC_lubbysub1, beliefPenaltyMPC_ubIdx1, beliefPenaltyMPC_Phi1);
beliefPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_5_17(beliefPenaltyMPC_Phi1, params->C2, beliefPenaltyMPC_V1);
beliefPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_10_17(beliefPenaltyMPC_Phi1, beliefPenaltyMPC_D1, beliefPenaltyMPC_W1);
beliefPenaltyMPC_LA_DENSE_MMTM_10_17_5(beliefPenaltyMPC_W1, beliefPenaltyMPC_V1, beliefPenaltyMPC_Ysd1);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_17(beliefPenaltyMPC_Phi1, beliefPenaltyMPC_rd1, beliefPenaltyMPC_Lbyrd1);
beliefPenaltyMPC_LA_DIAG_CHOL_LBUB_17_17_7(beliefPenaltyMPC_H0, beliefPenaltyMPC_llbbyslb2, beliefPenaltyMPC_lbIdx2, beliefPenaltyMPC_lubbysub2, beliefPenaltyMPC_ubIdx2, beliefPenaltyMPC_Phi2);
beliefPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_5_17(beliefPenaltyMPC_Phi2, params->C3, beliefPenaltyMPC_V2);
beliefPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_5_17(beliefPenaltyMPC_Phi2, beliefPenaltyMPC_D2, beliefPenaltyMPC_W2);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMTM_5_17_5(beliefPenaltyMPC_W2, beliefPenaltyMPC_V2, beliefPenaltyMPC_Ysd2);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_17(beliefPenaltyMPC_Phi2, beliefPenaltyMPC_rd2, beliefPenaltyMPC_Lbyrd2);
beliefPenaltyMPC_LA_DIAG_CHOL_LBUB_17_17_7(beliefPenaltyMPC_H0, beliefPenaltyMPC_llbbyslb3, beliefPenaltyMPC_lbIdx3, beliefPenaltyMPC_lubbysub3, beliefPenaltyMPC_ubIdx3, beliefPenaltyMPC_Phi3);
beliefPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_5_17(beliefPenaltyMPC_Phi3, params->C4, beliefPenaltyMPC_V3);
beliefPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_5_17(beliefPenaltyMPC_Phi3, beliefPenaltyMPC_D2, beliefPenaltyMPC_W3);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMTM_5_17_5(beliefPenaltyMPC_W3, beliefPenaltyMPC_V3, beliefPenaltyMPC_Ysd3);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_17(beliefPenaltyMPC_Phi3, beliefPenaltyMPC_rd3, beliefPenaltyMPC_Lbyrd3);
beliefPenaltyMPC_LA_DIAG_CHOL_LBUB_17_17_7(beliefPenaltyMPC_H0, beliefPenaltyMPC_llbbyslb4, beliefPenaltyMPC_lbIdx4, beliefPenaltyMPC_lubbysub4, beliefPenaltyMPC_ubIdx4, beliefPenaltyMPC_Phi4);
beliefPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_5_17(beliefPenaltyMPC_Phi4, params->C5, beliefPenaltyMPC_V4);
beliefPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_5_17(beliefPenaltyMPC_Phi4, beliefPenaltyMPC_D2, beliefPenaltyMPC_W4);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMTM_5_17_5(beliefPenaltyMPC_W4, beliefPenaltyMPC_V4, beliefPenaltyMPC_Ysd4);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_17(beliefPenaltyMPC_Phi4, beliefPenaltyMPC_rd4, beliefPenaltyMPC_Lbyrd4);
beliefPenaltyMPC_LA_DIAG_CHOL_LBUB_17_17_7(beliefPenaltyMPC_H0, beliefPenaltyMPC_llbbyslb5, beliefPenaltyMPC_lbIdx5, beliefPenaltyMPC_lubbysub5, beliefPenaltyMPC_ubIdx5, beliefPenaltyMPC_Phi5);
beliefPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_5_17(beliefPenaltyMPC_Phi5, params->C6, beliefPenaltyMPC_V5);
beliefPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_5_17(beliefPenaltyMPC_Phi5, beliefPenaltyMPC_D2, beliefPenaltyMPC_W5);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMTM_5_17_5(beliefPenaltyMPC_W5, beliefPenaltyMPC_V5, beliefPenaltyMPC_Ysd5);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_17(beliefPenaltyMPC_Phi5, beliefPenaltyMPC_rd5, beliefPenaltyMPC_Lbyrd5);
beliefPenaltyMPC_LA_DIAG_CHOL_LBUB_17_17_7(beliefPenaltyMPC_H0, beliefPenaltyMPC_llbbyslb6, beliefPenaltyMPC_lbIdx6, beliefPenaltyMPC_lubbysub6, beliefPenaltyMPC_ubIdx6, beliefPenaltyMPC_Phi6);
beliefPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_5_17(beliefPenaltyMPC_Phi6, params->C7, beliefPenaltyMPC_V6);
beliefPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_5_17(beliefPenaltyMPC_Phi6, beliefPenaltyMPC_D2, beliefPenaltyMPC_W6);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMTM_5_17_5(beliefPenaltyMPC_W6, beliefPenaltyMPC_V6, beliefPenaltyMPC_Ysd6);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_17(beliefPenaltyMPC_Phi6, beliefPenaltyMPC_rd6, beliefPenaltyMPC_Lbyrd6);
beliefPenaltyMPC_LA_DIAG_CHOL_LBUB_17_17_7(beliefPenaltyMPC_H0, beliefPenaltyMPC_llbbyslb7, beliefPenaltyMPC_lbIdx7, beliefPenaltyMPC_lubbysub7, beliefPenaltyMPC_ubIdx7, beliefPenaltyMPC_Phi7);
beliefPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_5_17(beliefPenaltyMPC_Phi7, params->C8, beliefPenaltyMPC_V7);
beliefPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_5_17(beliefPenaltyMPC_Phi7, beliefPenaltyMPC_D2, beliefPenaltyMPC_W7);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMTM_5_17_5(beliefPenaltyMPC_W7, beliefPenaltyMPC_V7, beliefPenaltyMPC_Ysd7);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_17(beliefPenaltyMPC_Phi7, beliefPenaltyMPC_rd7, beliefPenaltyMPC_Lbyrd7);
beliefPenaltyMPC_LA_DIAG_CHOL_LBUB_17_17_7(beliefPenaltyMPC_H0, beliefPenaltyMPC_llbbyslb8, beliefPenaltyMPC_lbIdx8, beliefPenaltyMPC_lubbysub8, beliefPenaltyMPC_ubIdx8, beliefPenaltyMPC_Phi8);
beliefPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_5_17(beliefPenaltyMPC_Phi8, params->C9, beliefPenaltyMPC_V8);
beliefPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_5_17(beliefPenaltyMPC_Phi8, beliefPenaltyMPC_D2, beliefPenaltyMPC_W8);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMTM_5_17_5(beliefPenaltyMPC_W8, beliefPenaltyMPC_V8, beliefPenaltyMPC_Ysd8);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_17(beliefPenaltyMPC_Phi8, beliefPenaltyMPC_rd8, beliefPenaltyMPC_Lbyrd8);
beliefPenaltyMPC_LA_DIAG_CHOL_ONELOOP_LBUB_5_5_5(beliefPenaltyMPC_H9, beliefPenaltyMPC_llbbyslb9, beliefPenaltyMPC_lbIdx9, beliefPenaltyMPC_lubbysub9, beliefPenaltyMPC_ubIdx9, beliefPenaltyMPC_Phi9);
beliefPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_5_5(beliefPenaltyMPC_Phi9, beliefPenaltyMPC_D9, beliefPenaltyMPC_W9);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_5(beliefPenaltyMPC_Phi9, beliefPenaltyMPC_rd9, beliefPenaltyMPC_Lbyrd9);
beliefPenaltyMPC_LA_DENSE_MMT2_10_17_17(beliefPenaltyMPC_V0, beliefPenaltyMPC_W1, beliefPenaltyMPC_Yd0);
beliefPenaltyMPC_LA_DENSE_MVMSUB2_10_17_17(beliefPenaltyMPC_V0, beliefPenaltyMPC_Lbyrd0, beliefPenaltyMPC_W1, beliefPenaltyMPC_Lbyrd1, beliefPenaltyMPC_re0, beliefPenaltyMPC_beta0);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMT2_5_17_17(beliefPenaltyMPC_V1, beliefPenaltyMPC_W2, beliefPenaltyMPC_Yd1);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_5_17_17(beliefPenaltyMPC_V1, beliefPenaltyMPC_Lbyrd1, beliefPenaltyMPC_W2, beliefPenaltyMPC_Lbyrd2, beliefPenaltyMPC_re1, beliefPenaltyMPC_beta1);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMT2_5_17_17(beliefPenaltyMPC_V2, beliefPenaltyMPC_W3, beliefPenaltyMPC_Yd2);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_5_17_17(beliefPenaltyMPC_V2, beliefPenaltyMPC_Lbyrd2, beliefPenaltyMPC_W3, beliefPenaltyMPC_Lbyrd3, beliefPenaltyMPC_re2, beliefPenaltyMPC_beta2);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMT2_5_17_17(beliefPenaltyMPC_V3, beliefPenaltyMPC_W4, beliefPenaltyMPC_Yd3);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_5_17_17(beliefPenaltyMPC_V3, beliefPenaltyMPC_Lbyrd3, beliefPenaltyMPC_W4, beliefPenaltyMPC_Lbyrd4, beliefPenaltyMPC_re3, beliefPenaltyMPC_beta3);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMT2_5_17_17(beliefPenaltyMPC_V4, beliefPenaltyMPC_W5, beliefPenaltyMPC_Yd4);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_5_17_17(beliefPenaltyMPC_V4, beliefPenaltyMPC_Lbyrd4, beliefPenaltyMPC_W5, beliefPenaltyMPC_Lbyrd5, beliefPenaltyMPC_re4, beliefPenaltyMPC_beta4);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMT2_5_17_17(beliefPenaltyMPC_V5, beliefPenaltyMPC_W6, beliefPenaltyMPC_Yd5);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_5_17_17(beliefPenaltyMPC_V5, beliefPenaltyMPC_Lbyrd5, beliefPenaltyMPC_W6, beliefPenaltyMPC_Lbyrd6, beliefPenaltyMPC_re5, beliefPenaltyMPC_beta5);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMT2_5_17_17(beliefPenaltyMPC_V6, beliefPenaltyMPC_W7, beliefPenaltyMPC_Yd6);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_5_17_17(beliefPenaltyMPC_V6, beliefPenaltyMPC_Lbyrd6, beliefPenaltyMPC_W7, beliefPenaltyMPC_Lbyrd7, beliefPenaltyMPC_re6, beliefPenaltyMPC_beta6);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMT2_5_17_17(beliefPenaltyMPC_V7, beliefPenaltyMPC_W8, beliefPenaltyMPC_Yd7);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_5_17_17(beliefPenaltyMPC_V7, beliefPenaltyMPC_Lbyrd7, beliefPenaltyMPC_W8, beliefPenaltyMPC_Lbyrd8, beliefPenaltyMPC_re7, beliefPenaltyMPC_beta7);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMT2_5_17_5(beliefPenaltyMPC_V8, beliefPenaltyMPC_W9, beliefPenaltyMPC_Yd8);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_5_17_5(beliefPenaltyMPC_V8, beliefPenaltyMPC_Lbyrd8, beliefPenaltyMPC_W9, beliefPenaltyMPC_Lbyrd9, beliefPenaltyMPC_re8, beliefPenaltyMPC_beta8);
beliefPenaltyMPC_LA_DENSE_CHOL_10(beliefPenaltyMPC_Yd0, beliefPenaltyMPC_Ld0);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_10(beliefPenaltyMPC_Ld0, beliefPenaltyMPC_beta0, beliefPenaltyMPC_yy0);
beliefPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_5_10(beliefPenaltyMPC_Ld0, beliefPenaltyMPC_Ysd1, beliefPenaltyMPC_Lsd1);
beliefPenaltyMPC_LA_DENSE_MMTSUB_5_10(beliefPenaltyMPC_Lsd1, beliefPenaltyMPC_Yd1);
beliefPenaltyMPC_LA_DENSE_CHOL_5(beliefPenaltyMPC_Yd1, beliefPenaltyMPC_Ld1);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_5_10(beliefPenaltyMPC_Lsd1, beliefPenaltyMPC_yy0, beliefPenaltyMPC_beta1, beliefPenaltyMPC_bmy1);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_5(beliefPenaltyMPC_Ld1, beliefPenaltyMPC_bmy1, beliefPenaltyMPC_yy1);
beliefPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_5_5(beliefPenaltyMPC_Ld1, beliefPenaltyMPC_Ysd2, beliefPenaltyMPC_Lsd2);
beliefPenaltyMPC_LA_DENSE_MMTSUB_5_5(beliefPenaltyMPC_Lsd2, beliefPenaltyMPC_Yd2);
beliefPenaltyMPC_LA_DENSE_CHOL_5(beliefPenaltyMPC_Yd2, beliefPenaltyMPC_Ld2);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_5_5(beliefPenaltyMPC_Lsd2, beliefPenaltyMPC_yy1, beliefPenaltyMPC_beta2, beliefPenaltyMPC_bmy2);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_5(beliefPenaltyMPC_Ld2, beliefPenaltyMPC_bmy2, beliefPenaltyMPC_yy2);
beliefPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_5_5(beliefPenaltyMPC_Ld2, beliefPenaltyMPC_Ysd3, beliefPenaltyMPC_Lsd3);
beliefPenaltyMPC_LA_DENSE_MMTSUB_5_5(beliefPenaltyMPC_Lsd3, beliefPenaltyMPC_Yd3);
beliefPenaltyMPC_LA_DENSE_CHOL_5(beliefPenaltyMPC_Yd3, beliefPenaltyMPC_Ld3);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_5_5(beliefPenaltyMPC_Lsd3, beliefPenaltyMPC_yy2, beliefPenaltyMPC_beta3, beliefPenaltyMPC_bmy3);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_5(beliefPenaltyMPC_Ld3, beliefPenaltyMPC_bmy3, beliefPenaltyMPC_yy3);
beliefPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_5_5(beliefPenaltyMPC_Ld3, beliefPenaltyMPC_Ysd4, beliefPenaltyMPC_Lsd4);
beliefPenaltyMPC_LA_DENSE_MMTSUB_5_5(beliefPenaltyMPC_Lsd4, beliefPenaltyMPC_Yd4);
beliefPenaltyMPC_LA_DENSE_CHOL_5(beliefPenaltyMPC_Yd4, beliefPenaltyMPC_Ld4);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_5_5(beliefPenaltyMPC_Lsd4, beliefPenaltyMPC_yy3, beliefPenaltyMPC_beta4, beliefPenaltyMPC_bmy4);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_5(beliefPenaltyMPC_Ld4, beliefPenaltyMPC_bmy4, beliefPenaltyMPC_yy4);
beliefPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_5_5(beliefPenaltyMPC_Ld4, beliefPenaltyMPC_Ysd5, beliefPenaltyMPC_Lsd5);
beliefPenaltyMPC_LA_DENSE_MMTSUB_5_5(beliefPenaltyMPC_Lsd5, beliefPenaltyMPC_Yd5);
beliefPenaltyMPC_LA_DENSE_CHOL_5(beliefPenaltyMPC_Yd5, beliefPenaltyMPC_Ld5);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_5_5(beliefPenaltyMPC_Lsd5, beliefPenaltyMPC_yy4, beliefPenaltyMPC_beta5, beliefPenaltyMPC_bmy5);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_5(beliefPenaltyMPC_Ld5, beliefPenaltyMPC_bmy5, beliefPenaltyMPC_yy5);
beliefPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_5_5(beliefPenaltyMPC_Ld5, beliefPenaltyMPC_Ysd6, beliefPenaltyMPC_Lsd6);
beliefPenaltyMPC_LA_DENSE_MMTSUB_5_5(beliefPenaltyMPC_Lsd6, beliefPenaltyMPC_Yd6);
beliefPenaltyMPC_LA_DENSE_CHOL_5(beliefPenaltyMPC_Yd6, beliefPenaltyMPC_Ld6);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_5_5(beliefPenaltyMPC_Lsd6, beliefPenaltyMPC_yy5, beliefPenaltyMPC_beta6, beliefPenaltyMPC_bmy6);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_5(beliefPenaltyMPC_Ld6, beliefPenaltyMPC_bmy6, beliefPenaltyMPC_yy6);
beliefPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_5_5(beliefPenaltyMPC_Ld6, beliefPenaltyMPC_Ysd7, beliefPenaltyMPC_Lsd7);
beliefPenaltyMPC_LA_DENSE_MMTSUB_5_5(beliefPenaltyMPC_Lsd7, beliefPenaltyMPC_Yd7);
beliefPenaltyMPC_LA_DENSE_CHOL_5(beliefPenaltyMPC_Yd7, beliefPenaltyMPC_Ld7);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_5_5(beliefPenaltyMPC_Lsd7, beliefPenaltyMPC_yy6, beliefPenaltyMPC_beta7, beliefPenaltyMPC_bmy7);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_5(beliefPenaltyMPC_Ld7, beliefPenaltyMPC_bmy7, beliefPenaltyMPC_yy7);
beliefPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_5_5(beliefPenaltyMPC_Ld7, beliefPenaltyMPC_Ysd8, beliefPenaltyMPC_Lsd8);
beliefPenaltyMPC_LA_DENSE_MMTSUB_5_5(beliefPenaltyMPC_Lsd8, beliefPenaltyMPC_Yd8);
beliefPenaltyMPC_LA_DENSE_CHOL_5(beliefPenaltyMPC_Yd8, beliefPenaltyMPC_Ld8);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_5_5(beliefPenaltyMPC_Lsd8, beliefPenaltyMPC_yy7, beliefPenaltyMPC_beta8, beliefPenaltyMPC_bmy8);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_5(beliefPenaltyMPC_Ld8, beliefPenaltyMPC_bmy8, beliefPenaltyMPC_yy8);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_5(beliefPenaltyMPC_Ld8, beliefPenaltyMPC_yy8, beliefPenaltyMPC_dvaff8);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_5_5(beliefPenaltyMPC_Lsd8, beliefPenaltyMPC_dvaff8, beliefPenaltyMPC_yy7, beliefPenaltyMPC_bmy7);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_5(beliefPenaltyMPC_Ld7, beliefPenaltyMPC_bmy7, beliefPenaltyMPC_dvaff7);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_5_5(beliefPenaltyMPC_Lsd7, beliefPenaltyMPC_dvaff7, beliefPenaltyMPC_yy6, beliefPenaltyMPC_bmy6);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_5(beliefPenaltyMPC_Ld6, beliefPenaltyMPC_bmy6, beliefPenaltyMPC_dvaff6);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_5_5(beliefPenaltyMPC_Lsd6, beliefPenaltyMPC_dvaff6, beliefPenaltyMPC_yy5, beliefPenaltyMPC_bmy5);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_5(beliefPenaltyMPC_Ld5, beliefPenaltyMPC_bmy5, beliefPenaltyMPC_dvaff5);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_5_5(beliefPenaltyMPC_Lsd5, beliefPenaltyMPC_dvaff5, beliefPenaltyMPC_yy4, beliefPenaltyMPC_bmy4);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_5(beliefPenaltyMPC_Ld4, beliefPenaltyMPC_bmy4, beliefPenaltyMPC_dvaff4);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_5_5(beliefPenaltyMPC_Lsd4, beliefPenaltyMPC_dvaff4, beliefPenaltyMPC_yy3, beliefPenaltyMPC_bmy3);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_5(beliefPenaltyMPC_Ld3, beliefPenaltyMPC_bmy3, beliefPenaltyMPC_dvaff3);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_5_5(beliefPenaltyMPC_Lsd3, beliefPenaltyMPC_dvaff3, beliefPenaltyMPC_yy2, beliefPenaltyMPC_bmy2);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_5(beliefPenaltyMPC_Ld2, beliefPenaltyMPC_bmy2, beliefPenaltyMPC_dvaff2);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_5_5(beliefPenaltyMPC_Lsd2, beliefPenaltyMPC_dvaff2, beliefPenaltyMPC_yy1, beliefPenaltyMPC_bmy1);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_5(beliefPenaltyMPC_Ld1, beliefPenaltyMPC_bmy1, beliefPenaltyMPC_dvaff1);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_5_10(beliefPenaltyMPC_Lsd1, beliefPenaltyMPC_dvaff1, beliefPenaltyMPC_yy0, beliefPenaltyMPC_bmy0);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_10(beliefPenaltyMPC_Ld0, beliefPenaltyMPC_bmy0, beliefPenaltyMPC_dvaff0);
beliefPenaltyMPC_LA_DENSE_MTVM_10_17(params->C1, beliefPenaltyMPC_dvaff0, beliefPenaltyMPC_grad_eq0);
beliefPenaltyMPC_LA_DENSE_MTVM2_5_17_10(params->C2, beliefPenaltyMPC_dvaff1, beliefPenaltyMPC_D1, beliefPenaltyMPC_dvaff0, beliefPenaltyMPC_grad_eq1);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_5_17_5(params->C3, beliefPenaltyMPC_dvaff2, beliefPenaltyMPC_D2, beliefPenaltyMPC_dvaff1, beliefPenaltyMPC_grad_eq2);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_5_17_5(params->C4, beliefPenaltyMPC_dvaff3, beliefPenaltyMPC_D2, beliefPenaltyMPC_dvaff2, beliefPenaltyMPC_grad_eq3);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_5_17_5(params->C5, beliefPenaltyMPC_dvaff4, beliefPenaltyMPC_D2, beliefPenaltyMPC_dvaff3, beliefPenaltyMPC_grad_eq4);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_5_17_5(params->C6, beliefPenaltyMPC_dvaff5, beliefPenaltyMPC_D2, beliefPenaltyMPC_dvaff4, beliefPenaltyMPC_grad_eq5);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_5_17_5(params->C7, beliefPenaltyMPC_dvaff6, beliefPenaltyMPC_D2, beliefPenaltyMPC_dvaff5, beliefPenaltyMPC_grad_eq6);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_5_17_5(params->C8, beliefPenaltyMPC_dvaff7, beliefPenaltyMPC_D2, beliefPenaltyMPC_dvaff6, beliefPenaltyMPC_grad_eq7);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_5_17_5(params->C9, beliefPenaltyMPC_dvaff8, beliefPenaltyMPC_D2, beliefPenaltyMPC_dvaff7, beliefPenaltyMPC_grad_eq8);
beliefPenaltyMPC_LA_DIAGZERO_MTVM_5_5(beliefPenaltyMPC_D9, beliefPenaltyMPC_dvaff8, beliefPenaltyMPC_grad_eq9);
beliefPenaltyMPC_LA_VSUB2_158(beliefPenaltyMPC_rd, beliefPenaltyMPC_grad_eq, beliefPenaltyMPC_rd);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_17(beliefPenaltyMPC_Phi0, beliefPenaltyMPC_rd0, beliefPenaltyMPC_dzaff0);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_17(beliefPenaltyMPC_Phi1, beliefPenaltyMPC_rd1, beliefPenaltyMPC_dzaff1);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_17(beliefPenaltyMPC_Phi2, beliefPenaltyMPC_rd2, beliefPenaltyMPC_dzaff2);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_17(beliefPenaltyMPC_Phi3, beliefPenaltyMPC_rd3, beliefPenaltyMPC_dzaff3);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_17(beliefPenaltyMPC_Phi4, beliefPenaltyMPC_rd4, beliefPenaltyMPC_dzaff4);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_17(beliefPenaltyMPC_Phi5, beliefPenaltyMPC_rd5, beliefPenaltyMPC_dzaff5);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_17(beliefPenaltyMPC_Phi6, beliefPenaltyMPC_rd6, beliefPenaltyMPC_dzaff6);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_17(beliefPenaltyMPC_Phi7, beliefPenaltyMPC_rd7, beliefPenaltyMPC_dzaff7);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_17(beliefPenaltyMPC_Phi8, beliefPenaltyMPC_rd8, beliefPenaltyMPC_dzaff8);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_5(beliefPenaltyMPC_Phi9, beliefPenaltyMPC_rd9, beliefPenaltyMPC_dzaff9);
beliefPenaltyMPC_LA_VSUB_INDEXED_17(beliefPenaltyMPC_dzaff0, beliefPenaltyMPC_lbIdx0, beliefPenaltyMPC_rilb0, beliefPenaltyMPC_dslbaff0);
beliefPenaltyMPC_LA_VSUB3_17(beliefPenaltyMPC_llbbyslb0, beliefPenaltyMPC_dslbaff0, beliefPenaltyMPC_llb0, beliefPenaltyMPC_dllbaff0);
beliefPenaltyMPC_LA_VSUB2_INDEXED_7(beliefPenaltyMPC_riub0, beliefPenaltyMPC_dzaff0, beliefPenaltyMPC_ubIdx0, beliefPenaltyMPC_dsubaff0);
beliefPenaltyMPC_LA_VSUB3_7(beliefPenaltyMPC_lubbysub0, beliefPenaltyMPC_dsubaff0, beliefPenaltyMPC_lub0, beliefPenaltyMPC_dlubaff0);
beliefPenaltyMPC_LA_VSUB_INDEXED_17(beliefPenaltyMPC_dzaff1, beliefPenaltyMPC_lbIdx1, beliefPenaltyMPC_rilb1, beliefPenaltyMPC_dslbaff1);
beliefPenaltyMPC_LA_VSUB3_17(beliefPenaltyMPC_llbbyslb1, beliefPenaltyMPC_dslbaff1, beliefPenaltyMPC_llb1, beliefPenaltyMPC_dllbaff1);
beliefPenaltyMPC_LA_VSUB2_INDEXED_7(beliefPenaltyMPC_riub1, beliefPenaltyMPC_dzaff1, beliefPenaltyMPC_ubIdx1, beliefPenaltyMPC_dsubaff1);
beliefPenaltyMPC_LA_VSUB3_7(beliefPenaltyMPC_lubbysub1, beliefPenaltyMPC_dsubaff1, beliefPenaltyMPC_lub1, beliefPenaltyMPC_dlubaff1);
beliefPenaltyMPC_LA_VSUB_INDEXED_17(beliefPenaltyMPC_dzaff2, beliefPenaltyMPC_lbIdx2, beliefPenaltyMPC_rilb2, beliefPenaltyMPC_dslbaff2);
beliefPenaltyMPC_LA_VSUB3_17(beliefPenaltyMPC_llbbyslb2, beliefPenaltyMPC_dslbaff2, beliefPenaltyMPC_llb2, beliefPenaltyMPC_dllbaff2);
beliefPenaltyMPC_LA_VSUB2_INDEXED_7(beliefPenaltyMPC_riub2, beliefPenaltyMPC_dzaff2, beliefPenaltyMPC_ubIdx2, beliefPenaltyMPC_dsubaff2);
beliefPenaltyMPC_LA_VSUB3_7(beliefPenaltyMPC_lubbysub2, beliefPenaltyMPC_dsubaff2, beliefPenaltyMPC_lub2, beliefPenaltyMPC_dlubaff2);
beliefPenaltyMPC_LA_VSUB_INDEXED_17(beliefPenaltyMPC_dzaff3, beliefPenaltyMPC_lbIdx3, beliefPenaltyMPC_rilb3, beliefPenaltyMPC_dslbaff3);
beliefPenaltyMPC_LA_VSUB3_17(beliefPenaltyMPC_llbbyslb3, beliefPenaltyMPC_dslbaff3, beliefPenaltyMPC_llb3, beliefPenaltyMPC_dllbaff3);
beliefPenaltyMPC_LA_VSUB2_INDEXED_7(beliefPenaltyMPC_riub3, beliefPenaltyMPC_dzaff3, beliefPenaltyMPC_ubIdx3, beliefPenaltyMPC_dsubaff3);
beliefPenaltyMPC_LA_VSUB3_7(beliefPenaltyMPC_lubbysub3, beliefPenaltyMPC_dsubaff3, beliefPenaltyMPC_lub3, beliefPenaltyMPC_dlubaff3);
beliefPenaltyMPC_LA_VSUB_INDEXED_17(beliefPenaltyMPC_dzaff4, beliefPenaltyMPC_lbIdx4, beliefPenaltyMPC_rilb4, beliefPenaltyMPC_dslbaff4);
beliefPenaltyMPC_LA_VSUB3_17(beliefPenaltyMPC_llbbyslb4, beliefPenaltyMPC_dslbaff4, beliefPenaltyMPC_llb4, beliefPenaltyMPC_dllbaff4);
beliefPenaltyMPC_LA_VSUB2_INDEXED_7(beliefPenaltyMPC_riub4, beliefPenaltyMPC_dzaff4, beliefPenaltyMPC_ubIdx4, beliefPenaltyMPC_dsubaff4);
beliefPenaltyMPC_LA_VSUB3_7(beliefPenaltyMPC_lubbysub4, beliefPenaltyMPC_dsubaff4, beliefPenaltyMPC_lub4, beliefPenaltyMPC_dlubaff4);
beliefPenaltyMPC_LA_VSUB_INDEXED_17(beliefPenaltyMPC_dzaff5, beliefPenaltyMPC_lbIdx5, beliefPenaltyMPC_rilb5, beliefPenaltyMPC_dslbaff5);
beliefPenaltyMPC_LA_VSUB3_17(beliefPenaltyMPC_llbbyslb5, beliefPenaltyMPC_dslbaff5, beliefPenaltyMPC_llb5, beliefPenaltyMPC_dllbaff5);
beliefPenaltyMPC_LA_VSUB2_INDEXED_7(beliefPenaltyMPC_riub5, beliefPenaltyMPC_dzaff5, beliefPenaltyMPC_ubIdx5, beliefPenaltyMPC_dsubaff5);
beliefPenaltyMPC_LA_VSUB3_7(beliefPenaltyMPC_lubbysub5, beliefPenaltyMPC_dsubaff5, beliefPenaltyMPC_lub5, beliefPenaltyMPC_dlubaff5);
beliefPenaltyMPC_LA_VSUB_INDEXED_17(beliefPenaltyMPC_dzaff6, beliefPenaltyMPC_lbIdx6, beliefPenaltyMPC_rilb6, beliefPenaltyMPC_dslbaff6);
beliefPenaltyMPC_LA_VSUB3_17(beliefPenaltyMPC_llbbyslb6, beliefPenaltyMPC_dslbaff6, beliefPenaltyMPC_llb6, beliefPenaltyMPC_dllbaff6);
beliefPenaltyMPC_LA_VSUB2_INDEXED_7(beliefPenaltyMPC_riub6, beliefPenaltyMPC_dzaff6, beliefPenaltyMPC_ubIdx6, beliefPenaltyMPC_dsubaff6);
beliefPenaltyMPC_LA_VSUB3_7(beliefPenaltyMPC_lubbysub6, beliefPenaltyMPC_dsubaff6, beliefPenaltyMPC_lub6, beliefPenaltyMPC_dlubaff6);
beliefPenaltyMPC_LA_VSUB_INDEXED_17(beliefPenaltyMPC_dzaff7, beliefPenaltyMPC_lbIdx7, beliefPenaltyMPC_rilb7, beliefPenaltyMPC_dslbaff7);
beliefPenaltyMPC_LA_VSUB3_17(beliefPenaltyMPC_llbbyslb7, beliefPenaltyMPC_dslbaff7, beliefPenaltyMPC_llb7, beliefPenaltyMPC_dllbaff7);
beliefPenaltyMPC_LA_VSUB2_INDEXED_7(beliefPenaltyMPC_riub7, beliefPenaltyMPC_dzaff7, beliefPenaltyMPC_ubIdx7, beliefPenaltyMPC_dsubaff7);
beliefPenaltyMPC_LA_VSUB3_7(beliefPenaltyMPC_lubbysub7, beliefPenaltyMPC_dsubaff7, beliefPenaltyMPC_lub7, beliefPenaltyMPC_dlubaff7);
beliefPenaltyMPC_LA_VSUB_INDEXED_17(beliefPenaltyMPC_dzaff8, beliefPenaltyMPC_lbIdx8, beliefPenaltyMPC_rilb8, beliefPenaltyMPC_dslbaff8);
beliefPenaltyMPC_LA_VSUB3_17(beliefPenaltyMPC_llbbyslb8, beliefPenaltyMPC_dslbaff8, beliefPenaltyMPC_llb8, beliefPenaltyMPC_dllbaff8);
beliefPenaltyMPC_LA_VSUB2_INDEXED_7(beliefPenaltyMPC_riub8, beliefPenaltyMPC_dzaff8, beliefPenaltyMPC_ubIdx8, beliefPenaltyMPC_dsubaff8);
beliefPenaltyMPC_LA_VSUB3_7(beliefPenaltyMPC_lubbysub8, beliefPenaltyMPC_dsubaff8, beliefPenaltyMPC_lub8, beliefPenaltyMPC_dlubaff8);
beliefPenaltyMPC_LA_VSUB_INDEXED_5(beliefPenaltyMPC_dzaff9, beliefPenaltyMPC_lbIdx9, beliefPenaltyMPC_rilb9, beliefPenaltyMPC_dslbaff9);
beliefPenaltyMPC_LA_VSUB3_5(beliefPenaltyMPC_llbbyslb9, beliefPenaltyMPC_dslbaff9, beliefPenaltyMPC_llb9, beliefPenaltyMPC_dllbaff9);
beliefPenaltyMPC_LA_VSUB2_INDEXED_5(beliefPenaltyMPC_riub9, beliefPenaltyMPC_dzaff9, beliefPenaltyMPC_ubIdx9, beliefPenaltyMPC_dsubaff9);
beliefPenaltyMPC_LA_VSUB3_5(beliefPenaltyMPC_lubbysub9, beliefPenaltyMPC_dsubaff9, beliefPenaltyMPC_lub9, beliefPenaltyMPC_dlubaff9);
info->lsit_aff = beliefPenaltyMPC_LINESEARCH_BACKTRACKING_AFFINE(beliefPenaltyMPC_l, beliefPenaltyMPC_s, beliefPenaltyMPC_dl_aff, beliefPenaltyMPC_ds_aff, &info->step_aff, &info->mu_aff);
if( info->lsit_aff == beliefPenaltyMPC_NOPROGRESS ){
exitcode = beliefPenaltyMPC_NOPROGRESS; break;
}
sigma_3rdroot = info->mu_aff / info->mu;
info->sigma = sigma_3rdroot*sigma_3rdroot*sigma_3rdroot;
musigma = info->mu * info->sigma;
beliefPenaltyMPC_LA_VSUB5_226(beliefPenaltyMPC_ds_aff, beliefPenaltyMPC_dl_aff, musigma, beliefPenaltyMPC_ccrhs);
beliefPenaltyMPC_LA_VSUB6_INDEXED_17_7_17(beliefPenaltyMPC_ccrhsub0, beliefPenaltyMPC_sub0, beliefPenaltyMPC_ubIdx0, beliefPenaltyMPC_ccrhsl0, beliefPenaltyMPC_slb0, beliefPenaltyMPC_lbIdx0, beliefPenaltyMPC_rd0);
beliefPenaltyMPC_LA_VSUB6_INDEXED_17_7_17(beliefPenaltyMPC_ccrhsub1, beliefPenaltyMPC_sub1, beliefPenaltyMPC_ubIdx1, beliefPenaltyMPC_ccrhsl1, beliefPenaltyMPC_slb1, beliefPenaltyMPC_lbIdx1, beliefPenaltyMPC_rd1);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_17(beliefPenaltyMPC_Phi0, beliefPenaltyMPC_rd0, beliefPenaltyMPC_Lbyrd0);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_17(beliefPenaltyMPC_Phi1, beliefPenaltyMPC_rd1, beliefPenaltyMPC_Lbyrd1);
beliefPenaltyMPC_LA_DENSE_2MVMADD_10_17_17(beliefPenaltyMPC_V0, beliefPenaltyMPC_Lbyrd0, beliefPenaltyMPC_W1, beliefPenaltyMPC_Lbyrd1, beliefPenaltyMPC_beta0);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_10(beliefPenaltyMPC_Ld0, beliefPenaltyMPC_beta0, beliefPenaltyMPC_yy0);
beliefPenaltyMPC_LA_VSUB6_INDEXED_17_7_17(beliefPenaltyMPC_ccrhsub2, beliefPenaltyMPC_sub2, beliefPenaltyMPC_ubIdx2, beliefPenaltyMPC_ccrhsl2, beliefPenaltyMPC_slb2, beliefPenaltyMPC_lbIdx2, beliefPenaltyMPC_rd2);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_17(beliefPenaltyMPC_Phi2, beliefPenaltyMPC_rd2, beliefPenaltyMPC_Lbyrd2);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_5_17_17(beliefPenaltyMPC_V1, beliefPenaltyMPC_Lbyrd1, beliefPenaltyMPC_W2, beliefPenaltyMPC_Lbyrd2, beliefPenaltyMPC_beta1);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_5_10(beliefPenaltyMPC_Lsd1, beliefPenaltyMPC_yy0, beliefPenaltyMPC_beta1, beliefPenaltyMPC_bmy1);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_5(beliefPenaltyMPC_Ld1, beliefPenaltyMPC_bmy1, beliefPenaltyMPC_yy1);
beliefPenaltyMPC_LA_VSUB6_INDEXED_17_7_17(beliefPenaltyMPC_ccrhsub3, beliefPenaltyMPC_sub3, beliefPenaltyMPC_ubIdx3, beliefPenaltyMPC_ccrhsl3, beliefPenaltyMPC_slb3, beliefPenaltyMPC_lbIdx3, beliefPenaltyMPC_rd3);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_17(beliefPenaltyMPC_Phi3, beliefPenaltyMPC_rd3, beliefPenaltyMPC_Lbyrd3);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_5_17_17(beliefPenaltyMPC_V2, beliefPenaltyMPC_Lbyrd2, beliefPenaltyMPC_W3, beliefPenaltyMPC_Lbyrd3, beliefPenaltyMPC_beta2);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_5_5(beliefPenaltyMPC_Lsd2, beliefPenaltyMPC_yy1, beliefPenaltyMPC_beta2, beliefPenaltyMPC_bmy2);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_5(beliefPenaltyMPC_Ld2, beliefPenaltyMPC_bmy2, beliefPenaltyMPC_yy2);
beliefPenaltyMPC_LA_VSUB6_INDEXED_17_7_17(beliefPenaltyMPC_ccrhsub4, beliefPenaltyMPC_sub4, beliefPenaltyMPC_ubIdx4, beliefPenaltyMPC_ccrhsl4, beliefPenaltyMPC_slb4, beliefPenaltyMPC_lbIdx4, beliefPenaltyMPC_rd4);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_17(beliefPenaltyMPC_Phi4, beliefPenaltyMPC_rd4, beliefPenaltyMPC_Lbyrd4);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_5_17_17(beliefPenaltyMPC_V3, beliefPenaltyMPC_Lbyrd3, beliefPenaltyMPC_W4, beliefPenaltyMPC_Lbyrd4, beliefPenaltyMPC_beta3);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_5_5(beliefPenaltyMPC_Lsd3, beliefPenaltyMPC_yy2, beliefPenaltyMPC_beta3, beliefPenaltyMPC_bmy3);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_5(beliefPenaltyMPC_Ld3, beliefPenaltyMPC_bmy3, beliefPenaltyMPC_yy3);
beliefPenaltyMPC_LA_VSUB6_INDEXED_17_7_17(beliefPenaltyMPC_ccrhsub5, beliefPenaltyMPC_sub5, beliefPenaltyMPC_ubIdx5, beliefPenaltyMPC_ccrhsl5, beliefPenaltyMPC_slb5, beliefPenaltyMPC_lbIdx5, beliefPenaltyMPC_rd5);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_17(beliefPenaltyMPC_Phi5, beliefPenaltyMPC_rd5, beliefPenaltyMPC_Lbyrd5);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_5_17_17(beliefPenaltyMPC_V4, beliefPenaltyMPC_Lbyrd4, beliefPenaltyMPC_W5, beliefPenaltyMPC_Lbyrd5, beliefPenaltyMPC_beta4);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_5_5(beliefPenaltyMPC_Lsd4, beliefPenaltyMPC_yy3, beliefPenaltyMPC_beta4, beliefPenaltyMPC_bmy4);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_5(beliefPenaltyMPC_Ld4, beliefPenaltyMPC_bmy4, beliefPenaltyMPC_yy4);
beliefPenaltyMPC_LA_VSUB6_INDEXED_17_7_17(beliefPenaltyMPC_ccrhsub6, beliefPenaltyMPC_sub6, beliefPenaltyMPC_ubIdx6, beliefPenaltyMPC_ccrhsl6, beliefPenaltyMPC_slb6, beliefPenaltyMPC_lbIdx6, beliefPenaltyMPC_rd6);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_17(beliefPenaltyMPC_Phi6, beliefPenaltyMPC_rd6, beliefPenaltyMPC_Lbyrd6);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_5_17_17(beliefPenaltyMPC_V5, beliefPenaltyMPC_Lbyrd5, beliefPenaltyMPC_W6, beliefPenaltyMPC_Lbyrd6, beliefPenaltyMPC_beta5);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_5_5(beliefPenaltyMPC_Lsd5, beliefPenaltyMPC_yy4, beliefPenaltyMPC_beta5, beliefPenaltyMPC_bmy5);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_5(beliefPenaltyMPC_Ld5, beliefPenaltyMPC_bmy5, beliefPenaltyMPC_yy5);
beliefPenaltyMPC_LA_VSUB6_INDEXED_17_7_17(beliefPenaltyMPC_ccrhsub7, beliefPenaltyMPC_sub7, beliefPenaltyMPC_ubIdx7, beliefPenaltyMPC_ccrhsl7, beliefPenaltyMPC_slb7, beliefPenaltyMPC_lbIdx7, beliefPenaltyMPC_rd7);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_17(beliefPenaltyMPC_Phi7, beliefPenaltyMPC_rd7, beliefPenaltyMPC_Lbyrd7);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_5_17_17(beliefPenaltyMPC_V6, beliefPenaltyMPC_Lbyrd6, beliefPenaltyMPC_W7, beliefPenaltyMPC_Lbyrd7, beliefPenaltyMPC_beta6);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_5_5(beliefPenaltyMPC_Lsd6, beliefPenaltyMPC_yy5, beliefPenaltyMPC_beta6, beliefPenaltyMPC_bmy6);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_5(beliefPenaltyMPC_Ld6, beliefPenaltyMPC_bmy6, beliefPenaltyMPC_yy6);
beliefPenaltyMPC_LA_VSUB6_INDEXED_17_7_17(beliefPenaltyMPC_ccrhsub8, beliefPenaltyMPC_sub8, beliefPenaltyMPC_ubIdx8, beliefPenaltyMPC_ccrhsl8, beliefPenaltyMPC_slb8, beliefPenaltyMPC_lbIdx8, beliefPenaltyMPC_rd8);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_17(beliefPenaltyMPC_Phi8, beliefPenaltyMPC_rd8, beliefPenaltyMPC_Lbyrd8);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_5_17_17(beliefPenaltyMPC_V7, beliefPenaltyMPC_Lbyrd7, beliefPenaltyMPC_W8, beliefPenaltyMPC_Lbyrd8, beliefPenaltyMPC_beta7);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_5_5(beliefPenaltyMPC_Lsd7, beliefPenaltyMPC_yy6, beliefPenaltyMPC_beta7, beliefPenaltyMPC_bmy7);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_5(beliefPenaltyMPC_Ld7, beliefPenaltyMPC_bmy7, beliefPenaltyMPC_yy7);
beliefPenaltyMPC_LA_VSUB6_INDEXED_5_5_5(beliefPenaltyMPC_ccrhsub9, beliefPenaltyMPC_sub9, beliefPenaltyMPC_ubIdx9, beliefPenaltyMPC_ccrhsl9, beliefPenaltyMPC_slb9, beliefPenaltyMPC_lbIdx9, beliefPenaltyMPC_rd9);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_5(beliefPenaltyMPC_Phi9, beliefPenaltyMPC_rd9, beliefPenaltyMPC_Lbyrd9);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_5_17_5(beliefPenaltyMPC_V8, beliefPenaltyMPC_Lbyrd8, beliefPenaltyMPC_W9, beliefPenaltyMPC_Lbyrd9, beliefPenaltyMPC_beta8);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_5_5(beliefPenaltyMPC_Lsd8, beliefPenaltyMPC_yy7, beliefPenaltyMPC_beta8, beliefPenaltyMPC_bmy8);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_5(beliefPenaltyMPC_Ld8, beliefPenaltyMPC_bmy8, beliefPenaltyMPC_yy8);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_5(beliefPenaltyMPC_Ld8, beliefPenaltyMPC_yy8, beliefPenaltyMPC_dvcc8);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_5_5(beliefPenaltyMPC_Lsd8, beliefPenaltyMPC_dvcc8, beliefPenaltyMPC_yy7, beliefPenaltyMPC_bmy7);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_5(beliefPenaltyMPC_Ld7, beliefPenaltyMPC_bmy7, beliefPenaltyMPC_dvcc7);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_5_5(beliefPenaltyMPC_Lsd7, beliefPenaltyMPC_dvcc7, beliefPenaltyMPC_yy6, beliefPenaltyMPC_bmy6);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_5(beliefPenaltyMPC_Ld6, beliefPenaltyMPC_bmy6, beliefPenaltyMPC_dvcc6);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_5_5(beliefPenaltyMPC_Lsd6, beliefPenaltyMPC_dvcc6, beliefPenaltyMPC_yy5, beliefPenaltyMPC_bmy5);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_5(beliefPenaltyMPC_Ld5, beliefPenaltyMPC_bmy5, beliefPenaltyMPC_dvcc5);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_5_5(beliefPenaltyMPC_Lsd5, beliefPenaltyMPC_dvcc5, beliefPenaltyMPC_yy4, beliefPenaltyMPC_bmy4);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_5(beliefPenaltyMPC_Ld4, beliefPenaltyMPC_bmy4, beliefPenaltyMPC_dvcc4);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_5_5(beliefPenaltyMPC_Lsd4, beliefPenaltyMPC_dvcc4, beliefPenaltyMPC_yy3, beliefPenaltyMPC_bmy3);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_5(beliefPenaltyMPC_Ld3, beliefPenaltyMPC_bmy3, beliefPenaltyMPC_dvcc3);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_5_5(beliefPenaltyMPC_Lsd3, beliefPenaltyMPC_dvcc3, beliefPenaltyMPC_yy2, beliefPenaltyMPC_bmy2);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_5(beliefPenaltyMPC_Ld2, beliefPenaltyMPC_bmy2, beliefPenaltyMPC_dvcc2);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_5_5(beliefPenaltyMPC_Lsd2, beliefPenaltyMPC_dvcc2, beliefPenaltyMPC_yy1, beliefPenaltyMPC_bmy1);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_5(beliefPenaltyMPC_Ld1, beliefPenaltyMPC_bmy1, beliefPenaltyMPC_dvcc1);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_5_10(beliefPenaltyMPC_Lsd1, beliefPenaltyMPC_dvcc1, beliefPenaltyMPC_yy0, beliefPenaltyMPC_bmy0);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_10(beliefPenaltyMPC_Ld0, beliefPenaltyMPC_bmy0, beliefPenaltyMPC_dvcc0);
beliefPenaltyMPC_LA_DENSE_MTVM_10_17(params->C1, beliefPenaltyMPC_dvcc0, beliefPenaltyMPC_grad_eq0);
beliefPenaltyMPC_LA_DENSE_MTVM2_5_17_10(params->C2, beliefPenaltyMPC_dvcc1, beliefPenaltyMPC_D1, beliefPenaltyMPC_dvcc0, beliefPenaltyMPC_grad_eq1);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_5_17_5(params->C3, beliefPenaltyMPC_dvcc2, beliefPenaltyMPC_D2, beliefPenaltyMPC_dvcc1, beliefPenaltyMPC_grad_eq2);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_5_17_5(params->C4, beliefPenaltyMPC_dvcc3, beliefPenaltyMPC_D2, beliefPenaltyMPC_dvcc2, beliefPenaltyMPC_grad_eq3);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_5_17_5(params->C5, beliefPenaltyMPC_dvcc4, beliefPenaltyMPC_D2, beliefPenaltyMPC_dvcc3, beliefPenaltyMPC_grad_eq4);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_5_17_5(params->C6, beliefPenaltyMPC_dvcc5, beliefPenaltyMPC_D2, beliefPenaltyMPC_dvcc4, beliefPenaltyMPC_grad_eq5);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_5_17_5(params->C7, beliefPenaltyMPC_dvcc6, beliefPenaltyMPC_D2, beliefPenaltyMPC_dvcc5, beliefPenaltyMPC_grad_eq6);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_5_17_5(params->C8, beliefPenaltyMPC_dvcc7, beliefPenaltyMPC_D2, beliefPenaltyMPC_dvcc6, beliefPenaltyMPC_grad_eq7);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_5_17_5(params->C9, beliefPenaltyMPC_dvcc8, beliefPenaltyMPC_D2, beliefPenaltyMPC_dvcc7, beliefPenaltyMPC_grad_eq8);
beliefPenaltyMPC_LA_DIAGZERO_MTVM_5_5(beliefPenaltyMPC_D9, beliefPenaltyMPC_dvcc8, beliefPenaltyMPC_grad_eq9);
beliefPenaltyMPC_LA_VSUB_158(beliefPenaltyMPC_rd, beliefPenaltyMPC_grad_eq, beliefPenaltyMPC_rd);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_17(beliefPenaltyMPC_Phi0, beliefPenaltyMPC_rd0, beliefPenaltyMPC_dzcc0);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_17(beliefPenaltyMPC_Phi1, beliefPenaltyMPC_rd1, beliefPenaltyMPC_dzcc1);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_17(beliefPenaltyMPC_Phi2, beliefPenaltyMPC_rd2, beliefPenaltyMPC_dzcc2);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_17(beliefPenaltyMPC_Phi3, beliefPenaltyMPC_rd3, beliefPenaltyMPC_dzcc3);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_17(beliefPenaltyMPC_Phi4, beliefPenaltyMPC_rd4, beliefPenaltyMPC_dzcc4);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_17(beliefPenaltyMPC_Phi5, beliefPenaltyMPC_rd5, beliefPenaltyMPC_dzcc5);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_17(beliefPenaltyMPC_Phi6, beliefPenaltyMPC_rd6, beliefPenaltyMPC_dzcc6);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_17(beliefPenaltyMPC_Phi7, beliefPenaltyMPC_rd7, beliefPenaltyMPC_dzcc7);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_17(beliefPenaltyMPC_Phi8, beliefPenaltyMPC_rd8, beliefPenaltyMPC_dzcc8);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_5(beliefPenaltyMPC_Phi9, beliefPenaltyMPC_rd9, beliefPenaltyMPC_dzcc9);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_17(beliefPenaltyMPC_ccrhsl0, beliefPenaltyMPC_slb0, beliefPenaltyMPC_llbbyslb0, beliefPenaltyMPC_dzcc0, beliefPenaltyMPC_lbIdx0, beliefPenaltyMPC_dllbcc0);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_7(beliefPenaltyMPC_ccrhsub0, beliefPenaltyMPC_sub0, beliefPenaltyMPC_lubbysub0, beliefPenaltyMPC_dzcc0, beliefPenaltyMPC_ubIdx0, beliefPenaltyMPC_dlubcc0);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_17(beliefPenaltyMPC_ccrhsl1, beliefPenaltyMPC_slb1, beliefPenaltyMPC_llbbyslb1, beliefPenaltyMPC_dzcc1, beliefPenaltyMPC_lbIdx1, beliefPenaltyMPC_dllbcc1);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_7(beliefPenaltyMPC_ccrhsub1, beliefPenaltyMPC_sub1, beliefPenaltyMPC_lubbysub1, beliefPenaltyMPC_dzcc1, beliefPenaltyMPC_ubIdx1, beliefPenaltyMPC_dlubcc1);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_17(beliefPenaltyMPC_ccrhsl2, beliefPenaltyMPC_slb2, beliefPenaltyMPC_llbbyslb2, beliefPenaltyMPC_dzcc2, beliefPenaltyMPC_lbIdx2, beliefPenaltyMPC_dllbcc2);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_7(beliefPenaltyMPC_ccrhsub2, beliefPenaltyMPC_sub2, beliefPenaltyMPC_lubbysub2, beliefPenaltyMPC_dzcc2, beliefPenaltyMPC_ubIdx2, beliefPenaltyMPC_dlubcc2);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_17(beliefPenaltyMPC_ccrhsl3, beliefPenaltyMPC_slb3, beliefPenaltyMPC_llbbyslb3, beliefPenaltyMPC_dzcc3, beliefPenaltyMPC_lbIdx3, beliefPenaltyMPC_dllbcc3);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_7(beliefPenaltyMPC_ccrhsub3, beliefPenaltyMPC_sub3, beliefPenaltyMPC_lubbysub3, beliefPenaltyMPC_dzcc3, beliefPenaltyMPC_ubIdx3, beliefPenaltyMPC_dlubcc3);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_17(beliefPenaltyMPC_ccrhsl4, beliefPenaltyMPC_slb4, beliefPenaltyMPC_llbbyslb4, beliefPenaltyMPC_dzcc4, beliefPenaltyMPC_lbIdx4, beliefPenaltyMPC_dllbcc4);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_7(beliefPenaltyMPC_ccrhsub4, beliefPenaltyMPC_sub4, beliefPenaltyMPC_lubbysub4, beliefPenaltyMPC_dzcc4, beliefPenaltyMPC_ubIdx4, beliefPenaltyMPC_dlubcc4);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_17(beliefPenaltyMPC_ccrhsl5, beliefPenaltyMPC_slb5, beliefPenaltyMPC_llbbyslb5, beliefPenaltyMPC_dzcc5, beliefPenaltyMPC_lbIdx5, beliefPenaltyMPC_dllbcc5);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_7(beliefPenaltyMPC_ccrhsub5, beliefPenaltyMPC_sub5, beliefPenaltyMPC_lubbysub5, beliefPenaltyMPC_dzcc5, beliefPenaltyMPC_ubIdx5, beliefPenaltyMPC_dlubcc5);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_17(beliefPenaltyMPC_ccrhsl6, beliefPenaltyMPC_slb6, beliefPenaltyMPC_llbbyslb6, beliefPenaltyMPC_dzcc6, beliefPenaltyMPC_lbIdx6, beliefPenaltyMPC_dllbcc6);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_7(beliefPenaltyMPC_ccrhsub6, beliefPenaltyMPC_sub6, beliefPenaltyMPC_lubbysub6, beliefPenaltyMPC_dzcc6, beliefPenaltyMPC_ubIdx6, beliefPenaltyMPC_dlubcc6);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_17(beliefPenaltyMPC_ccrhsl7, beliefPenaltyMPC_slb7, beliefPenaltyMPC_llbbyslb7, beliefPenaltyMPC_dzcc7, beliefPenaltyMPC_lbIdx7, beliefPenaltyMPC_dllbcc7);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_7(beliefPenaltyMPC_ccrhsub7, beliefPenaltyMPC_sub7, beliefPenaltyMPC_lubbysub7, beliefPenaltyMPC_dzcc7, beliefPenaltyMPC_ubIdx7, beliefPenaltyMPC_dlubcc7);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_17(beliefPenaltyMPC_ccrhsl8, beliefPenaltyMPC_slb8, beliefPenaltyMPC_llbbyslb8, beliefPenaltyMPC_dzcc8, beliefPenaltyMPC_lbIdx8, beliefPenaltyMPC_dllbcc8);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_7(beliefPenaltyMPC_ccrhsub8, beliefPenaltyMPC_sub8, beliefPenaltyMPC_lubbysub8, beliefPenaltyMPC_dzcc8, beliefPenaltyMPC_ubIdx8, beliefPenaltyMPC_dlubcc8);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_5(beliefPenaltyMPC_ccrhsl9, beliefPenaltyMPC_slb9, beliefPenaltyMPC_llbbyslb9, beliefPenaltyMPC_dzcc9, beliefPenaltyMPC_lbIdx9, beliefPenaltyMPC_dllbcc9);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(beliefPenaltyMPC_ccrhsub9, beliefPenaltyMPC_sub9, beliefPenaltyMPC_lubbysub9, beliefPenaltyMPC_dzcc9, beliefPenaltyMPC_ubIdx9, beliefPenaltyMPC_dlubcc9);
beliefPenaltyMPC_LA_VSUB7_226(beliefPenaltyMPC_l, beliefPenaltyMPC_ccrhs, beliefPenaltyMPC_s, beliefPenaltyMPC_dl_cc, beliefPenaltyMPC_ds_cc);
beliefPenaltyMPC_LA_VADD_158(beliefPenaltyMPC_dz_cc, beliefPenaltyMPC_dz_aff);
beliefPenaltyMPC_LA_VADD_50(beliefPenaltyMPC_dv_cc, beliefPenaltyMPC_dv_aff);
beliefPenaltyMPC_LA_VADD_226(beliefPenaltyMPC_dl_cc, beliefPenaltyMPC_dl_aff);
beliefPenaltyMPC_LA_VADD_226(beliefPenaltyMPC_ds_cc, beliefPenaltyMPC_ds_aff);
info->lsit_cc = beliefPenaltyMPC_LINESEARCH_BACKTRACKING_COMBINED(beliefPenaltyMPC_z, beliefPenaltyMPC_v, beliefPenaltyMPC_l, beliefPenaltyMPC_s, beliefPenaltyMPC_dz_cc, beliefPenaltyMPC_dv_cc, beliefPenaltyMPC_dl_cc, beliefPenaltyMPC_ds_cc, &info->step_cc, &info->mu);
if( info->lsit_cc == beliefPenaltyMPC_NOPROGRESS ){
exitcode = beliefPenaltyMPC_NOPROGRESS; break;
}
info->it++;
}
output->z1[0] = beliefPenaltyMPC_z0[0];
output->z1[1] = beliefPenaltyMPC_z0[1];
output->z1[2] = beliefPenaltyMPC_z0[2];
output->z1[3] = beliefPenaltyMPC_z0[3];
output->z1[4] = beliefPenaltyMPC_z0[4];
output->z1[5] = beliefPenaltyMPC_z0[5];
output->z1[6] = beliefPenaltyMPC_z0[6];
output->z2[0] = beliefPenaltyMPC_z1[0];
output->z2[1] = beliefPenaltyMPC_z1[1];
output->z2[2] = beliefPenaltyMPC_z1[2];
output->z2[3] = beliefPenaltyMPC_z1[3];
output->z2[4] = beliefPenaltyMPC_z1[4];
output->z2[5] = beliefPenaltyMPC_z1[5];
output->z2[6] = beliefPenaltyMPC_z1[6];
output->z3[0] = beliefPenaltyMPC_z2[0];
output->z3[1] = beliefPenaltyMPC_z2[1];
output->z3[2] = beliefPenaltyMPC_z2[2];
output->z3[3] = beliefPenaltyMPC_z2[3];
output->z3[4] = beliefPenaltyMPC_z2[4];
output->z3[5] = beliefPenaltyMPC_z2[5];
output->z3[6] = beliefPenaltyMPC_z2[6];
output->z4[0] = beliefPenaltyMPC_z3[0];
output->z4[1] = beliefPenaltyMPC_z3[1];
output->z4[2] = beliefPenaltyMPC_z3[2];
output->z4[3] = beliefPenaltyMPC_z3[3];
output->z4[4] = beliefPenaltyMPC_z3[4];
output->z4[5] = beliefPenaltyMPC_z3[5];
output->z4[6] = beliefPenaltyMPC_z3[6];
output->z5[0] = beliefPenaltyMPC_z4[0];
output->z5[1] = beliefPenaltyMPC_z4[1];
output->z5[2] = beliefPenaltyMPC_z4[2];
output->z5[3] = beliefPenaltyMPC_z4[3];
output->z5[4] = beliefPenaltyMPC_z4[4];
output->z5[5] = beliefPenaltyMPC_z4[5];
output->z5[6] = beliefPenaltyMPC_z4[6];
output->z6[0] = beliefPenaltyMPC_z5[0];
output->z6[1] = beliefPenaltyMPC_z5[1];
output->z6[2] = beliefPenaltyMPC_z5[2];
output->z6[3] = beliefPenaltyMPC_z5[3];
output->z6[4] = beliefPenaltyMPC_z5[4];
output->z6[5] = beliefPenaltyMPC_z5[5];
output->z6[6] = beliefPenaltyMPC_z5[6];
output->z7[0] = beliefPenaltyMPC_z6[0];
output->z7[1] = beliefPenaltyMPC_z6[1];
output->z7[2] = beliefPenaltyMPC_z6[2];
output->z7[3] = beliefPenaltyMPC_z6[3];
output->z7[4] = beliefPenaltyMPC_z6[4];
output->z7[5] = beliefPenaltyMPC_z6[5];
output->z7[6] = beliefPenaltyMPC_z6[6];
output->z8[0] = beliefPenaltyMPC_z7[0];
output->z8[1] = beliefPenaltyMPC_z7[1];
output->z8[2] = beliefPenaltyMPC_z7[2];
output->z8[3] = beliefPenaltyMPC_z7[3];
output->z8[4] = beliefPenaltyMPC_z7[4];
output->z8[5] = beliefPenaltyMPC_z7[5];
output->z8[6] = beliefPenaltyMPC_z7[6];
output->z9[0] = beliefPenaltyMPC_z8[0];
output->z9[1] = beliefPenaltyMPC_z8[1];
output->z9[2] = beliefPenaltyMPC_z8[2];
output->z9[3] = beliefPenaltyMPC_z8[3];
output->z9[4] = beliefPenaltyMPC_z8[4];
output->z9[5] = beliefPenaltyMPC_z8[5];
output->z9[6] = beliefPenaltyMPC_z8[6];
output->z10[0] = beliefPenaltyMPC_z9[0];
output->z10[1] = beliefPenaltyMPC_z9[1];
output->z10[2] = beliefPenaltyMPC_z9[2];
output->z10[3] = beliefPenaltyMPC_z9[3];
output->z10[4] = beliefPenaltyMPC_z9[4];

#if beliefPenaltyMPC_SET_TIMING == 1
info->solvetime = beliefPenaltyMPC_toc(&solvertimer);
#if beliefPenaltyMPC_SET_PRINTLEVEL > 0 && beliefPenaltyMPC_SET_TIMING == 1
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
