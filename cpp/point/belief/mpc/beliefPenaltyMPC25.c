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
 * Initializes a vector of length 413 with a value.
 */
void beliefPenaltyMPC_LA_INITIALIZEVECTOR_413(beliefPenaltyMPC_FLOAT* vec, beliefPenaltyMPC_FLOAT value)
{
	int i;
	for( i=0; i<413; i++ )
	{
		vec[i] = value;
	}
}


/*
 * Initializes a vector of length 125 with a value.
 */
void beliefPenaltyMPC_LA_INITIALIZEVECTOR_125(beliefPenaltyMPC_FLOAT* vec, beliefPenaltyMPC_FLOAT value)
{
	int i;
	for( i=0; i<125; i++ )
	{
		vec[i] = value;
	}
}


/*
 * Initializes a vector of length 586 with a value.
 */
void beliefPenaltyMPC_LA_INITIALIZEVECTOR_586(beliefPenaltyMPC_FLOAT* vec, beliefPenaltyMPC_FLOAT value)
{
	int i;
	for( i=0; i<586; i++ )
	{
		vec[i] = value;
	}
}


/* 
 * Calculates a dot product and adds it to a variable: z += x'*y; 
 * This function is for vectors of length 586.
 */
void beliefPenaltyMPC_LA_DOTACC_586(beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *y, beliefPenaltyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<586; i++ ){
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
 * of length 413.
 */
void beliefPenaltyMPC_LA_VVADD3_413(beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *v, beliefPenaltyMPC_FLOAT *w, beliefPenaltyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<413; i++ ){
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
 * Vector subtraction z = -x - y for vectors of length 413.
 */
void beliefPenaltyMPC_LA_VSUB2_413(beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *y, beliefPenaltyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<413; i++){
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
        for( i=0; i<586; i++ ){
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
        if( i == 586 ){
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
    *mu_aff = mymu / (beliefPenaltyMPC_FLOAT)586;
    return lsIt;
}


/*
 * Vector subtraction x = u.*v - a where a is a scalar
*  and x,u,v are vectors of length 586.
 */
void beliefPenaltyMPC_LA_VSUB5_586(beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *v, beliefPenaltyMPC_FLOAT a, beliefPenaltyMPC_FLOAT *x)
{
	int i;
	for( i=0; i<586; i++){
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
 * Vector subtraction z = x - y for vectors of length 413.
 */
void beliefPenaltyMPC_LA_VSUB_413(beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *y, beliefPenaltyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<413; i++){
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
 * Computes ds = -l.\(r + s.*dl) for vectors of length 586.
 */
void beliefPenaltyMPC_LA_VSUB7_586(beliefPenaltyMPC_FLOAT *l, beliefPenaltyMPC_FLOAT *r, beliefPenaltyMPC_FLOAT *s, beliefPenaltyMPC_FLOAT *dl, beliefPenaltyMPC_FLOAT *ds)
{
	int i;
	for( i=0; i<586; i++){
		ds[i] = -(r[i] + s[i]*dl[i])/l[i];
	}
}


/*
 * Vector addition x = x + y for vectors of length 413.
 */
void beliefPenaltyMPC_LA_VADD_413(beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *y)
{
	int i;
	for( i=0; i<413; i++){
		x[i] += y[i];
	}
}


/*
 * Vector addition x = x + y for vectors of length 125.
 */
void beliefPenaltyMPC_LA_VADD_125(beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *y)
{
	int i;
	for( i=0; i<125; i++){
		x[i] += y[i];
	}
}


/*
 * Vector addition x = x + y for vectors of length 586.
 */
void beliefPenaltyMPC_LA_VADD_586(beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *y)
{
	int i;
	for( i=0; i<586; i++){
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
        for( i=0; i<586; i++ ){
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
        if( i == 586 ){
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
    for( i=0; i<413; i++ ){
        z[i] += a_gamma*dz[i];
    }
    
    /* equality constraint multipliers */
    for( i=0; i<125; i++ ){
        v[i] += a_gamma*dv[i];
    }
    
    /* inequality constraint multipliers & slacks, also update mu */
    *mu = 0;
    for( i=0; i<586; i++ ){
        dltemp = l[i] + a_gamma*dl[i]; l[i] = dltemp;
        dstemp = s[i] + a_gamma*ds[i]; s[i] = dstemp;
        *mu += dltemp*dstemp;
    }
    
    *a = a_gamma;
    *mu /= (beliefPenaltyMPC_FLOAT)586;
    return lsIt;
}




/* VARIABLE DEFINITIONS ------------------------------------------------ */
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_z[413];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_v[125];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_dz_aff[413];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_dv_aff[125];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_grad_cost[413];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_grad_eq[413];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_rd[413];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_l[586];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_s[586];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_lbys[586];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_dl_aff[586];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ds_aff[586];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_dz_cc[413];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_dv_cc[125];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_dl_cc[586];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ds_cc[586];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ccrhs[586];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_grad_ineq[413];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_H00[17] = {0.0000000000000000E+000, 0.0000000000000000E+000, 2.0000000000000000E+001, 2.0000000000000000E+001, 2.0000000000000000E+001, 2.0000000000000000E+000, 2.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_z00 = beliefPenaltyMPC_z + 0;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzaff00 = beliefPenaltyMPC_dz_aff + 0;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzcc00 = beliefPenaltyMPC_dz_cc + 0;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_rd00 = beliefPenaltyMPC_rd + 0;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lbyrd00[17];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_cost00 = beliefPenaltyMPC_grad_cost + 0;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_eq00 = beliefPenaltyMPC_grad_eq + 0;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_ineq00 = beliefPenaltyMPC_grad_ineq + 0;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ctv00[17];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_v00 = beliefPenaltyMPC_v + 0;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_re00[10];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_beta00[10];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_betacc00[10];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvaff00 = beliefPenaltyMPC_dv_aff + 0;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvcc00 = beliefPenaltyMPC_dv_cc + 0;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_V00[170];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Yd00[55];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ld00[55];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_yy00[10];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_bmy00[10];
int beliefPenaltyMPC_lbIdx00[17] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llb00 = beliefPenaltyMPC_l + 0;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_slb00 = beliefPenaltyMPC_s + 0;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llbbyslb00 = beliefPenaltyMPC_lbys + 0;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_rilb00[17];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbaff00 = beliefPenaltyMPC_dl_aff + 0;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbaff00 = beliefPenaltyMPC_ds_aff + 0;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbcc00 = beliefPenaltyMPC_dl_cc + 0;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbcc00 = beliefPenaltyMPC_ds_cc + 0;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsl00 = beliefPenaltyMPC_ccrhs + 0;
int beliefPenaltyMPC_ubIdx00[7] = {0, 1, 2, 3, 4, 5, 6};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lub00 = beliefPenaltyMPC_l + 17;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_sub00 = beliefPenaltyMPC_s + 17;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lubbysub00 = beliefPenaltyMPC_lbys + 17;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_riub00[7];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubaff00 = beliefPenaltyMPC_dl_aff + 17;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubaff00 = beliefPenaltyMPC_ds_aff + 17;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubcc00 = beliefPenaltyMPC_dl_cc + 17;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubcc00 = beliefPenaltyMPC_ds_cc + 17;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsub00 = beliefPenaltyMPC_ccrhs + 17;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Phi00[17];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_z01 = beliefPenaltyMPC_z + 17;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzaff01 = beliefPenaltyMPC_dz_aff + 17;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzcc01 = beliefPenaltyMPC_dz_cc + 17;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_rd01 = beliefPenaltyMPC_rd + 17;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lbyrd01[17];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_cost01 = beliefPenaltyMPC_grad_cost + 17;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_eq01 = beliefPenaltyMPC_grad_eq + 17;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_ineq01 = beliefPenaltyMPC_grad_ineq + 17;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ctv01[17];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_v01 = beliefPenaltyMPC_v + 10;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_re01[5];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_beta01[5];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_betacc01[5];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvaff01 = beliefPenaltyMPC_dv_aff + 10;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvcc01 = beliefPenaltyMPC_dv_cc + 10;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_V01[85];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Yd01[15];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ld01[15];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_yy01[5];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_bmy01[5];
int beliefPenaltyMPC_lbIdx01[17] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llb01 = beliefPenaltyMPC_l + 24;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_slb01 = beliefPenaltyMPC_s + 24;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llbbyslb01 = beliefPenaltyMPC_lbys + 24;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_rilb01[17];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbaff01 = beliefPenaltyMPC_dl_aff + 24;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbaff01 = beliefPenaltyMPC_ds_aff + 24;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbcc01 = beliefPenaltyMPC_dl_cc + 24;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbcc01 = beliefPenaltyMPC_ds_cc + 24;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsl01 = beliefPenaltyMPC_ccrhs + 24;
int beliefPenaltyMPC_ubIdx01[7] = {0, 1, 2, 3, 4, 5, 6};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lub01 = beliefPenaltyMPC_l + 41;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_sub01 = beliefPenaltyMPC_s + 41;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lubbysub01 = beliefPenaltyMPC_lbys + 41;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_riub01[7];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubaff01 = beliefPenaltyMPC_dl_aff + 41;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubaff01 = beliefPenaltyMPC_ds_aff + 41;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubcc01 = beliefPenaltyMPC_dl_cc + 41;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubcc01 = beliefPenaltyMPC_ds_cc + 41;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsub01 = beliefPenaltyMPC_ccrhs + 41;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Phi01[17];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_D01[170] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, -1.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 
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
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_W01[170];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ysd01[50];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lsd01[50];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_z02 = beliefPenaltyMPC_z + 34;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzaff02 = beliefPenaltyMPC_dz_aff + 34;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzcc02 = beliefPenaltyMPC_dz_cc + 34;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_rd02 = beliefPenaltyMPC_rd + 34;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lbyrd02[17];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_cost02 = beliefPenaltyMPC_grad_cost + 34;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_eq02 = beliefPenaltyMPC_grad_eq + 34;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_ineq02 = beliefPenaltyMPC_grad_ineq + 34;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ctv02[17];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_v02 = beliefPenaltyMPC_v + 15;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_re02[5];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_beta02[5];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_betacc02[5];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvaff02 = beliefPenaltyMPC_dv_aff + 15;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvcc02 = beliefPenaltyMPC_dv_cc + 15;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_V02[85];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Yd02[15];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ld02[15];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_yy02[5];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_bmy02[5];
int beliefPenaltyMPC_lbIdx02[17] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llb02 = beliefPenaltyMPC_l + 48;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_slb02 = beliefPenaltyMPC_s + 48;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llbbyslb02 = beliefPenaltyMPC_lbys + 48;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_rilb02[17];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbaff02 = beliefPenaltyMPC_dl_aff + 48;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbaff02 = beliefPenaltyMPC_ds_aff + 48;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbcc02 = beliefPenaltyMPC_dl_cc + 48;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbcc02 = beliefPenaltyMPC_ds_cc + 48;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsl02 = beliefPenaltyMPC_ccrhs + 48;
int beliefPenaltyMPC_ubIdx02[7] = {0, 1, 2, 3, 4, 5, 6};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lub02 = beliefPenaltyMPC_l + 65;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_sub02 = beliefPenaltyMPC_s + 65;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lubbysub02 = beliefPenaltyMPC_lbys + 65;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_riub02[7];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubaff02 = beliefPenaltyMPC_dl_aff + 65;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubaff02 = beliefPenaltyMPC_ds_aff + 65;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubcc02 = beliefPenaltyMPC_dl_cc + 65;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubcc02 = beliefPenaltyMPC_ds_cc + 65;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsub02 = beliefPenaltyMPC_ccrhs + 65;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Phi02[17];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_D02[17] = {-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000};
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_W02[17];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ysd02[25];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lsd02[25];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_z03 = beliefPenaltyMPC_z + 51;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzaff03 = beliefPenaltyMPC_dz_aff + 51;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzcc03 = beliefPenaltyMPC_dz_cc + 51;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_rd03 = beliefPenaltyMPC_rd + 51;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lbyrd03[17];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_cost03 = beliefPenaltyMPC_grad_cost + 51;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_eq03 = beliefPenaltyMPC_grad_eq + 51;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_ineq03 = beliefPenaltyMPC_grad_ineq + 51;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ctv03[17];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_v03 = beliefPenaltyMPC_v + 20;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_re03[5];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_beta03[5];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_betacc03[5];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvaff03 = beliefPenaltyMPC_dv_aff + 20;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvcc03 = beliefPenaltyMPC_dv_cc + 20;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_V03[85];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Yd03[15];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ld03[15];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_yy03[5];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_bmy03[5];
int beliefPenaltyMPC_lbIdx03[17] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llb03 = beliefPenaltyMPC_l + 72;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_slb03 = beliefPenaltyMPC_s + 72;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llbbyslb03 = beliefPenaltyMPC_lbys + 72;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_rilb03[17];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbaff03 = beliefPenaltyMPC_dl_aff + 72;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbaff03 = beliefPenaltyMPC_ds_aff + 72;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbcc03 = beliefPenaltyMPC_dl_cc + 72;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbcc03 = beliefPenaltyMPC_ds_cc + 72;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsl03 = beliefPenaltyMPC_ccrhs + 72;
int beliefPenaltyMPC_ubIdx03[7] = {0, 1, 2, 3, 4, 5, 6};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lub03 = beliefPenaltyMPC_l + 89;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_sub03 = beliefPenaltyMPC_s + 89;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lubbysub03 = beliefPenaltyMPC_lbys + 89;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_riub03[7];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubaff03 = beliefPenaltyMPC_dl_aff + 89;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubaff03 = beliefPenaltyMPC_ds_aff + 89;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubcc03 = beliefPenaltyMPC_dl_cc + 89;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubcc03 = beliefPenaltyMPC_ds_cc + 89;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsub03 = beliefPenaltyMPC_ccrhs + 89;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Phi03[17];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_W03[17];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ysd03[25];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lsd03[25];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_z04 = beliefPenaltyMPC_z + 68;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzaff04 = beliefPenaltyMPC_dz_aff + 68;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzcc04 = beliefPenaltyMPC_dz_cc + 68;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_rd04 = beliefPenaltyMPC_rd + 68;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lbyrd04[17];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_cost04 = beliefPenaltyMPC_grad_cost + 68;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_eq04 = beliefPenaltyMPC_grad_eq + 68;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_ineq04 = beliefPenaltyMPC_grad_ineq + 68;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ctv04[17];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_v04 = beliefPenaltyMPC_v + 25;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_re04[5];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_beta04[5];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_betacc04[5];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvaff04 = beliefPenaltyMPC_dv_aff + 25;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvcc04 = beliefPenaltyMPC_dv_cc + 25;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_V04[85];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Yd04[15];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ld04[15];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_yy04[5];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_bmy04[5];
int beliefPenaltyMPC_lbIdx04[17] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llb04 = beliefPenaltyMPC_l + 96;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_slb04 = beliefPenaltyMPC_s + 96;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llbbyslb04 = beliefPenaltyMPC_lbys + 96;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_rilb04[17];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbaff04 = beliefPenaltyMPC_dl_aff + 96;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbaff04 = beliefPenaltyMPC_ds_aff + 96;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbcc04 = beliefPenaltyMPC_dl_cc + 96;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbcc04 = beliefPenaltyMPC_ds_cc + 96;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsl04 = beliefPenaltyMPC_ccrhs + 96;
int beliefPenaltyMPC_ubIdx04[7] = {0, 1, 2, 3, 4, 5, 6};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lub04 = beliefPenaltyMPC_l + 113;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_sub04 = beliefPenaltyMPC_s + 113;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lubbysub04 = beliefPenaltyMPC_lbys + 113;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_riub04[7];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubaff04 = beliefPenaltyMPC_dl_aff + 113;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubaff04 = beliefPenaltyMPC_ds_aff + 113;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubcc04 = beliefPenaltyMPC_dl_cc + 113;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubcc04 = beliefPenaltyMPC_ds_cc + 113;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsub04 = beliefPenaltyMPC_ccrhs + 113;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Phi04[17];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_W04[17];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ysd04[25];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lsd04[25];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_z05 = beliefPenaltyMPC_z + 85;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzaff05 = beliefPenaltyMPC_dz_aff + 85;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzcc05 = beliefPenaltyMPC_dz_cc + 85;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_rd05 = beliefPenaltyMPC_rd + 85;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lbyrd05[17];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_cost05 = beliefPenaltyMPC_grad_cost + 85;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_eq05 = beliefPenaltyMPC_grad_eq + 85;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_ineq05 = beliefPenaltyMPC_grad_ineq + 85;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ctv05[17];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_v05 = beliefPenaltyMPC_v + 30;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_re05[5];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_beta05[5];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_betacc05[5];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvaff05 = beliefPenaltyMPC_dv_aff + 30;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvcc05 = beliefPenaltyMPC_dv_cc + 30;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_V05[85];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Yd05[15];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ld05[15];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_yy05[5];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_bmy05[5];
int beliefPenaltyMPC_lbIdx05[17] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llb05 = beliefPenaltyMPC_l + 120;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_slb05 = beliefPenaltyMPC_s + 120;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llbbyslb05 = beliefPenaltyMPC_lbys + 120;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_rilb05[17];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbaff05 = beliefPenaltyMPC_dl_aff + 120;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbaff05 = beliefPenaltyMPC_ds_aff + 120;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbcc05 = beliefPenaltyMPC_dl_cc + 120;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbcc05 = beliefPenaltyMPC_ds_cc + 120;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsl05 = beliefPenaltyMPC_ccrhs + 120;
int beliefPenaltyMPC_ubIdx05[7] = {0, 1, 2, 3, 4, 5, 6};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lub05 = beliefPenaltyMPC_l + 137;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_sub05 = beliefPenaltyMPC_s + 137;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lubbysub05 = beliefPenaltyMPC_lbys + 137;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_riub05[7];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubaff05 = beliefPenaltyMPC_dl_aff + 137;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubaff05 = beliefPenaltyMPC_ds_aff + 137;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubcc05 = beliefPenaltyMPC_dl_cc + 137;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubcc05 = beliefPenaltyMPC_ds_cc + 137;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsub05 = beliefPenaltyMPC_ccrhs + 137;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Phi05[17];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_W05[17];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ysd05[25];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lsd05[25];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_z06 = beliefPenaltyMPC_z + 102;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzaff06 = beliefPenaltyMPC_dz_aff + 102;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzcc06 = beliefPenaltyMPC_dz_cc + 102;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_rd06 = beliefPenaltyMPC_rd + 102;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lbyrd06[17];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_cost06 = beliefPenaltyMPC_grad_cost + 102;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_eq06 = beliefPenaltyMPC_grad_eq + 102;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_ineq06 = beliefPenaltyMPC_grad_ineq + 102;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ctv06[17];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_v06 = beliefPenaltyMPC_v + 35;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_re06[5];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_beta06[5];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_betacc06[5];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvaff06 = beliefPenaltyMPC_dv_aff + 35;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvcc06 = beliefPenaltyMPC_dv_cc + 35;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_V06[85];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Yd06[15];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ld06[15];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_yy06[5];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_bmy06[5];
int beliefPenaltyMPC_lbIdx06[17] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llb06 = beliefPenaltyMPC_l + 144;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_slb06 = beliefPenaltyMPC_s + 144;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llbbyslb06 = beliefPenaltyMPC_lbys + 144;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_rilb06[17];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbaff06 = beliefPenaltyMPC_dl_aff + 144;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbaff06 = beliefPenaltyMPC_ds_aff + 144;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbcc06 = beliefPenaltyMPC_dl_cc + 144;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbcc06 = beliefPenaltyMPC_ds_cc + 144;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsl06 = beliefPenaltyMPC_ccrhs + 144;
int beliefPenaltyMPC_ubIdx06[7] = {0, 1, 2, 3, 4, 5, 6};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lub06 = beliefPenaltyMPC_l + 161;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_sub06 = beliefPenaltyMPC_s + 161;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lubbysub06 = beliefPenaltyMPC_lbys + 161;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_riub06[7];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubaff06 = beliefPenaltyMPC_dl_aff + 161;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubaff06 = beliefPenaltyMPC_ds_aff + 161;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubcc06 = beliefPenaltyMPC_dl_cc + 161;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubcc06 = beliefPenaltyMPC_ds_cc + 161;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsub06 = beliefPenaltyMPC_ccrhs + 161;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Phi06[17];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_W06[17];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ysd06[25];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lsd06[25];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_z07 = beliefPenaltyMPC_z + 119;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzaff07 = beliefPenaltyMPC_dz_aff + 119;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzcc07 = beliefPenaltyMPC_dz_cc + 119;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_rd07 = beliefPenaltyMPC_rd + 119;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lbyrd07[17];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_cost07 = beliefPenaltyMPC_grad_cost + 119;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_eq07 = beliefPenaltyMPC_grad_eq + 119;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_ineq07 = beliefPenaltyMPC_grad_ineq + 119;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ctv07[17];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_v07 = beliefPenaltyMPC_v + 40;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_re07[5];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_beta07[5];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_betacc07[5];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvaff07 = beliefPenaltyMPC_dv_aff + 40;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvcc07 = beliefPenaltyMPC_dv_cc + 40;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_V07[85];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Yd07[15];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ld07[15];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_yy07[5];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_bmy07[5];
int beliefPenaltyMPC_lbIdx07[17] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llb07 = beliefPenaltyMPC_l + 168;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_slb07 = beliefPenaltyMPC_s + 168;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llbbyslb07 = beliefPenaltyMPC_lbys + 168;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_rilb07[17];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbaff07 = beliefPenaltyMPC_dl_aff + 168;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbaff07 = beliefPenaltyMPC_ds_aff + 168;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbcc07 = beliefPenaltyMPC_dl_cc + 168;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbcc07 = beliefPenaltyMPC_ds_cc + 168;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsl07 = beliefPenaltyMPC_ccrhs + 168;
int beliefPenaltyMPC_ubIdx07[7] = {0, 1, 2, 3, 4, 5, 6};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lub07 = beliefPenaltyMPC_l + 185;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_sub07 = beliefPenaltyMPC_s + 185;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lubbysub07 = beliefPenaltyMPC_lbys + 185;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_riub07[7];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubaff07 = beliefPenaltyMPC_dl_aff + 185;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubaff07 = beliefPenaltyMPC_ds_aff + 185;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubcc07 = beliefPenaltyMPC_dl_cc + 185;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubcc07 = beliefPenaltyMPC_ds_cc + 185;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsub07 = beliefPenaltyMPC_ccrhs + 185;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Phi07[17];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_W07[17];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ysd07[25];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lsd07[25];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_z08 = beliefPenaltyMPC_z + 136;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzaff08 = beliefPenaltyMPC_dz_aff + 136;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzcc08 = beliefPenaltyMPC_dz_cc + 136;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_rd08 = beliefPenaltyMPC_rd + 136;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lbyrd08[17];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_cost08 = beliefPenaltyMPC_grad_cost + 136;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_eq08 = beliefPenaltyMPC_grad_eq + 136;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_ineq08 = beliefPenaltyMPC_grad_ineq + 136;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ctv08[17];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_v08 = beliefPenaltyMPC_v + 45;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_re08[5];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_beta08[5];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_betacc08[5];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvaff08 = beliefPenaltyMPC_dv_aff + 45;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvcc08 = beliefPenaltyMPC_dv_cc + 45;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_V08[85];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Yd08[15];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ld08[15];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_yy08[5];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_bmy08[5];
int beliefPenaltyMPC_lbIdx08[17] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llb08 = beliefPenaltyMPC_l + 192;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_slb08 = beliefPenaltyMPC_s + 192;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llbbyslb08 = beliefPenaltyMPC_lbys + 192;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_rilb08[17];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbaff08 = beliefPenaltyMPC_dl_aff + 192;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbaff08 = beliefPenaltyMPC_ds_aff + 192;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbcc08 = beliefPenaltyMPC_dl_cc + 192;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbcc08 = beliefPenaltyMPC_ds_cc + 192;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsl08 = beliefPenaltyMPC_ccrhs + 192;
int beliefPenaltyMPC_ubIdx08[7] = {0, 1, 2, 3, 4, 5, 6};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lub08 = beliefPenaltyMPC_l + 209;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_sub08 = beliefPenaltyMPC_s + 209;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lubbysub08 = beliefPenaltyMPC_lbys + 209;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_riub08[7];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubaff08 = beliefPenaltyMPC_dl_aff + 209;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubaff08 = beliefPenaltyMPC_ds_aff + 209;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubcc08 = beliefPenaltyMPC_dl_cc + 209;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubcc08 = beliefPenaltyMPC_ds_cc + 209;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsub08 = beliefPenaltyMPC_ccrhs + 209;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Phi08[17];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_W08[17];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ysd08[25];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lsd08[25];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_z09 = beliefPenaltyMPC_z + 153;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzaff09 = beliefPenaltyMPC_dz_aff + 153;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzcc09 = beliefPenaltyMPC_dz_cc + 153;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_rd09 = beliefPenaltyMPC_rd + 153;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lbyrd09[17];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_cost09 = beliefPenaltyMPC_grad_cost + 153;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_eq09 = beliefPenaltyMPC_grad_eq + 153;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_ineq09 = beliefPenaltyMPC_grad_ineq + 153;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ctv09[17];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_v09 = beliefPenaltyMPC_v + 50;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_re09[5];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_beta09[5];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_betacc09[5];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvaff09 = beliefPenaltyMPC_dv_aff + 50;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvcc09 = beliefPenaltyMPC_dv_cc + 50;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_V09[85];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Yd09[15];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ld09[15];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_yy09[5];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_bmy09[5];
int beliefPenaltyMPC_lbIdx09[17] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llb09 = beliefPenaltyMPC_l + 216;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_slb09 = beliefPenaltyMPC_s + 216;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llbbyslb09 = beliefPenaltyMPC_lbys + 216;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_rilb09[17];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbaff09 = beliefPenaltyMPC_dl_aff + 216;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbaff09 = beliefPenaltyMPC_ds_aff + 216;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbcc09 = beliefPenaltyMPC_dl_cc + 216;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbcc09 = beliefPenaltyMPC_ds_cc + 216;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsl09 = beliefPenaltyMPC_ccrhs + 216;
int beliefPenaltyMPC_ubIdx09[7] = {0, 1, 2, 3, 4, 5, 6};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lub09 = beliefPenaltyMPC_l + 233;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_sub09 = beliefPenaltyMPC_s + 233;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lubbysub09 = beliefPenaltyMPC_lbys + 233;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_riub09[7];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubaff09 = beliefPenaltyMPC_dl_aff + 233;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubaff09 = beliefPenaltyMPC_ds_aff + 233;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubcc09 = beliefPenaltyMPC_dl_cc + 233;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubcc09 = beliefPenaltyMPC_ds_cc + 233;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsub09 = beliefPenaltyMPC_ccrhs + 233;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Phi09[17];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_W09[17];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ysd09[25];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lsd09[25];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_z10 = beliefPenaltyMPC_z + 170;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzaff10 = beliefPenaltyMPC_dz_aff + 170;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzcc10 = beliefPenaltyMPC_dz_cc + 170;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_rd10 = beliefPenaltyMPC_rd + 170;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lbyrd10[17];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_cost10 = beliefPenaltyMPC_grad_cost + 170;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_eq10 = beliefPenaltyMPC_grad_eq + 170;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_ineq10 = beliefPenaltyMPC_grad_ineq + 170;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ctv10[17];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_v10 = beliefPenaltyMPC_v + 55;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_re10[5];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_beta10[5];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_betacc10[5];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvaff10 = beliefPenaltyMPC_dv_aff + 55;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvcc10 = beliefPenaltyMPC_dv_cc + 55;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_V10[85];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Yd10[15];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ld10[15];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_yy10[5];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_bmy10[5];
int beliefPenaltyMPC_lbIdx10[17] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llb10 = beliefPenaltyMPC_l + 240;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_slb10 = beliefPenaltyMPC_s + 240;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llbbyslb10 = beliefPenaltyMPC_lbys + 240;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_rilb10[17];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbaff10 = beliefPenaltyMPC_dl_aff + 240;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbaff10 = beliefPenaltyMPC_ds_aff + 240;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbcc10 = beliefPenaltyMPC_dl_cc + 240;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbcc10 = beliefPenaltyMPC_ds_cc + 240;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsl10 = beliefPenaltyMPC_ccrhs + 240;
int beliefPenaltyMPC_ubIdx10[7] = {0, 1, 2, 3, 4, 5, 6};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lub10 = beliefPenaltyMPC_l + 257;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_sub10 = beliefPenaltyMPC_s + 257;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lubbysub10 = beliefPenaltyMPC_lbys + 257;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_riub10[7];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubaff10 = beliefPenaltyMPC_dl_aff + 257;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubaff10 = beliefPenaltyMPC_ds_aff + 257;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubcc10 = beliefPenaltyMPC_dl_cc + 257;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubcc10 = beliefPenaltyMPC_ds_cc + 257;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsub10 = beliefPenaltyMPC_ccrhs + 257;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Phi10[17];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_W10[17];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ysd10[25];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lsd10[25];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_z11 = beliefPenaltyMPC_z + 187;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzaff11 = beliefPenaltyMPC_dz_aff + 187;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzcc11 = beliefPenaltyMPC_dz_cc + 187;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_rd11 = beliefPenaltyMPC_rd + 187;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lbyrd11[17];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_cost11 = beliefPenaltyMPC_grad_cost + 187;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_eq11 = beliefPenaltyMPC_grad_eq + 187;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_ineq11 = beliefPenaltyMPC_grad_ineq + 187;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ctv11[17];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_v11 = beliefPenaltyMPC_v + 60;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_re11[5];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_beta11[5];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_betacc11[5];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvaff11 = beliefPenaltyMPC_dv_aff + 60;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvcc11 = beliefPenaltyMPC_dv_cc + 60;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_V11[85];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Yd11[15];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ld11[15];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_yy11[5];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_bmy11[5];
int beliefPenaltyMPC_lbIdx11[17] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llb11 = beliefPenaltyMPC_l + 264;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_slb11 = beliefPenaltyMPC_s + 264;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llbbyslb11 = beliefPenaltyMPC_lbys + 264;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_rilb11[17];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbaff11 = beliefPenaltyMPC_dl_aff + 264;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbaff11 = beliefPenaltyMPC_ds_aff + 264;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbcc11 = beliefPenaltyMPC_dl_cc + 264;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbcc11 = beliefPenaltyMPC_ds_cc + 264;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsl11 = beliefPenaltyMPC_ccrhs + 264;
int beliefPenaltyMPC_ubIdx11[7] = {0, 1, 2, 3, 4, 5, 6};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lub11 = beliefPenaltyMPC_l + 281;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_sub11 = beliefPenaltyMPC_s + 281;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lubbysub11 = beliefPenaltyMPC_lbys + 281;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_riub11[7];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubaff11 = beliefPenaltyMPC_dl_aff + 281;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubaff11 = beliefPenaltyMPC_ds_aff + 281;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubcc11 = beliefPenaltyMPC_dl_cc + 281;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubcc11 = beliefPenaltyMPC_ds_cc + 281;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsub11 = beliefPenaltyMPC_ccrhs + 281;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Phi11[17];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_W11[17];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ysd11[25];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lsd11[25];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_z12 = beliefPenaltyMPC_z + 204;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzaff12 = beliefPenaltyMPC_dz_aff + 204;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzcc12 = beliefPenaltyMPC_dz_cc + 204;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_rd12 = beliefPenaltyMPC_rd + 204;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lbyrd12[17];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_cost12 = beliefPenaltyMPC_grad_cost + 204;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_eq12 = beliefPenaltyMPC_grad_eq + 204;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_ineq12 = beliefPenaltyMPC_grad_ineq + 204;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ctv12[17];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_v12 = beliefPenaltyMPC_v + 65;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_re12[5];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_beta12[5];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_betacc12[5];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvaff12 = beliefPenaltyMPC_dv_aff + 65;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvcc12 = beliefPenaltyMPC_dv_cc + 65;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_V12[85];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Yd12[15];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ld12[15];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_yy12[5];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_bmy12[5];
int beliefPenaltyMPC_lbIdx12[17] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llb12 = beliefPenaltyMPC_l + 288;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_slb12 = beliefPenaltyMPC_s + 288;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llbbyslb12 = beliefPenaltyMPC_lbys + 288;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_rilb12[17];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbaff12 = beliefPenaltyMPC_dl_aff + 288;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbaff12 = beliefPenaltyMPC_ds_aff + 288;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbcc12 = beliefPenaltyMPC_dl_cc + 288;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbcc12 = beliefPenaltyMPC_ds_cc + 288;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsl12 = beliefPenaltyMPC_ccrhs + 288;
int beliefPenaltyMPC_ubIdx12[7] = {0, 1, 2, 3, 4, 5, 6};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lub12 = beliefPenaltyMPC_l + 305;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_sub12 = beliefPenaltyMPC_s + 305;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lubbysub12 = beliefPenaltyMPC_lbys + 305;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_riub12[7];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubaff12 = beliefPenaltyMPC_dl_aff + 305;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubaff12 = beliefPenaltyMPC_ds_aff + 305;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubcc12 = beliefPenaltyMPC_dl_cc + 305;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubcc12 = beliefPenaltyMPC_ds_cc + 305;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsub12 = beliefPenaltyMPC_ccrhs + 305;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Phi12[17];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_W12[17];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ysd12[25];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lsd12[25];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_z13 = beliefPenaltyMPC_z + 221;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzaff13 = beliefPenaltyMPC_dz_aff + 221;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzcc13 = beliefPenaltyMPC_dz_cc + 221;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_rd13 = beliefPenaltyMPC_rd + 221;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lbyrd13[17];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_cost13 = beliefPenaltyMPC_grad_cost + 221;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_eq13 = beliefPenaltyMPC_grad_eq + 221;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_ineq13 = beliefPenaltyMPC_grad_ineq + 221;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ctv13[17];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_v13 = beliefPenaltyMPC_v + 70;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_re13[5];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_beta13[5];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_betacc13[5];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvaff13 = beliefPenaltyMPC_dv_aff + 70;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvcc13 = beliefPenaltyMPC_dv_cc + 70;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_V13[85];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Yd13[15];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ld13[15];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_yy13[5];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_bmy13[5];
int beliefPenaltyMPC_lbIdx13[17] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llb13 = beliefPenaltyMPC_l + 312;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_slb13 = beliefPenaltyMPC_s + 312;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llbbyslb13 = beliefPenaltyMPC_lbys + 312;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_rilb13[17];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbaff13 = beliefPenaltyMPC_dl_aff + 312;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbaff13 = beliefPenaltyMPC_ds_aff + 312;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbcc13 = beliefPenaltyMPC_dl_cc + 312;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbcc13 = beliefPenaltyMPC_ds_cc + 312;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsl13 = beliefPenaltyMPC_ccrhs + 312;
int beliefPenaltyMPC_ubIdx13[7] = {0, 1, 2, 3, 4, 5, 6};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lub13 = beliefPenaltyMPC_l + 329;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_sub13 = beliefPenaltyMPC_s + 329;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lubbysub13 = beliefPenaltyMPC_lbys + 329;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_riub13[7];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubaff13 = beliefPenaltyMPC_dl_aff + 329;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubaff13 = beliefPenaltyMPC_ds_aff + 329;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubcc13 = beliefPenaltyMPC_dl_cc + 329;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubcc13 = beliefPenaltyMPC_ds_cc + 329;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsub13 = beliefPenaltyMPC_ccrhs + 329;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Phi13[17];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_W13[17];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ysd13[25];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lsd13[25];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_z14 = beliefPenaltyMPC_z + 238;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzaff14 = beliefPenaltyMPC_dz_aff + 238;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzcc14 = beliefPenaltyMPC_dz_cc + 238;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_rd14 = beliefPenaltyMPC_rd + 238;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lbyrd14[17];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_cost14 = beliefPenaltyMPC_grad_cost + 238;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_eq14 = beliefPenaltyMPC_grad_eq + 238;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_ineq14 = beliefPenaltyMPC_grad_ineq + 238;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ctv14[17];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_v14 = beliefPenaltyMPC_v + 75;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_re14[5];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_beta14[5];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_betacc14[5];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvaff14 = beliefPenaltyMPC_dv_aff + 75;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvcc14 = beliefPenaltyMPC_dv_cc + 75;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_V14[85];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Yd14[15];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ld14[15];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_yy14[5];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_bmy14[5];
int beliefPenaltyMPC_lbIdx14[17] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llb14 = beliefPenaltyMPC_l + 336;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_slb14 = beliefPenaltyMPC_s + 336;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llbbyslb14 = beliefPenaltyMPC_lbys + 336;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_rilb14[17];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbaff14 = beliefPenaltyMPC_dl_aff + 336;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbaff14 = beliefPenaltyMPC_ds_aff + 336;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbcc14 = beliefPenaltyMPC_dl_cc + 336;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbcc14 = beliefPenaltyMPC_ds_cc + 336;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsl14 = beliefPenaltyMPC_ccrhs + 336;
int beliefPenaltyMPC_ubIdx14[7] = {0, 1, 2, 3, 4, 5, 6};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lub14 = beliefPenaltyMPC_l + 353;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_sub14 = beliefPenaltyMPC_s + 353;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lubbysub14 = beliefPenaltyMPC_lbys + 353;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_riub14[7];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubaff14 = beliefPenaltyMPC_dl_aff + 353;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubaff14 = beliefPenaltyMPC_ds_aff + 353;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubcc14 = beliefPenaltyMPC_dl_cc + 353;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubcc14 = beliefPenaltyMPC_ds_cc + 353;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsub14 = beliefPenaltyMPC_ccrhs + 353;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Phi14[17];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_W14[17];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ysd14[25];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lsd14[25];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_z15 = beliefPenaltyMPC_z + 255;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzaff15 = beliefPenaltyMPC_dz_aff + 255;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzcc15 = beliefPenaltyMPC_dz_cc + 255;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_rd15 = beliefPenaltyMPC_rd + 255;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lbyrd15[17];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_cost15 = beliefPenaltyMPC_grad_cost + 255;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_eq15 = beliefPenaltyMPC_grad_eq + 255;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_ineq15 = beliefPenaltyMPC_grad_ineq + 255;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ctv15[17];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_v15 = beliefPenaltyMPC_v + 80;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_re15[5];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_beta15[5];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_betacc15[5];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvaff15 = beliefPenaltyMPC_dv_aff + 80;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvcc15 = beliefPenaltyMPC_dv_cc + 80;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_V15[85];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Yd15[15];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ld15[15];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_yy15[5];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_bmy15[5];
int beliefPenaltyMPC_lbIdx15[17] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llb15 = beliefPenaltyMPC_l + 360;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_slb15 = beliefPenaltyMPC_s + 360;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llbbyslb15 = beliefPenaltyMPC_lbys + 360;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_rilb15[17];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbaff15 = beliefPenaltyMPC_dl_aff + 360;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbaff15 = beliefPenaltyMPC_ds_aff + 360;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbcc15 = beliefPenaltyMPC_dl_cc + 360;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbcc15 = beliefPenaltyMPC_ds_cc + 360;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsl15 = beliefPenaltyMPC_ccrhs + 360;
int beliefPenaltyMPC_ubIdx15[7] = {0, 1, 2, 3, 4, 5, 6};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lub15 = beliefPenaltyMPC_l + 377;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_sub15 = beliefPenaltyMPC_s + 377;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lubbysub15 = beliefPenaltyMPC_lbys + 377;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_riub15[7];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubaff15 = beliefPenaltyMPC_dl_aff + 377;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubaff15 = beliefPenaltyMPC_ds_aff + 377;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubcc15 = beliefPenaltyMPC_dl_cc + 377;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubcc15 = beliefPenaltyMPC_ds_cc + 377;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsub15 = beliefPenaltyMPC_ccrhs + 377;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Phi15[17];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_W15[17];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ysd15[25];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lsd15[25];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_z16 = beliefPenaltyMPC_z + 272;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzaff16 = beliefPenaltyMPC_dz_aff + 272;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzcc16 = beliefPenaltyMPC_dz_cc + 272;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_rd16 = beliefPenaltyMPC_rd + 272;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lbyrd16[17];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_cost16 = beliefPenaltyMPC_grad_cost + 272;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_eq16 = beliefPenaltyMPC_grad_eq + 272;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_ineq16 = beliefPenaltyMPC_grad_ineq + 272;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ctv16[17];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_v16 = beliefPenaltyMPC_v + 85;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_re16[5];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_beta16[5];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_betacc16[5];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvaff16 = beliefPenaltyMPC_dv_aff + 85;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvcc16 = beliefPenaltyMPC_dv_cc + 85;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_V16[85];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Yd16[15];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ld16[15];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_yy16[5];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_bmy16[5];
int beliefPenaltyMPC_lbIdx16[17] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llb16 = beliefPenaltyMPC_l + 384;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_slb16 = beliefPenaltyMPC_s + 384;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llbbyslb16 = beliefPenaltyMPC_lbys + 384;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_rilb16[17];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbaff16 = beliefPenaltyMPC_dl_aff + 384;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbaff16 = beliefPenaltyMPC_ds_aff + 384;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbcc16 = beliefPenaltyMPC_dl_cc + 384;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbcc16 = beliefPenaltyMPC_ds_cc + 384;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsl16 = beliefPenaltyMPC_ccrhs + 384;
int beliefPenaltyMPC_ubIdx16[7] = {0, 1, 2, 3, 4, 5, 6};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lub16 = beliefPenaltyMPC_l + 401;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_sub16 = beliefPenaltyMPC_s + 401;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lubbysub16 = beliefPenaltyMPC_lbys + 401;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_riub16[7];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubaff16 = beliefPenaltyMPC_dl_aff + 401;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubaff16 = beliefPenaltyMPC_ds_aff + 401;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubcc16 = beliefPenaltyMPC_dl_cc + 401;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubcc16 = beliefPenaltyMPC_ds_cc + 401;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsub16 = beliefPenaltyMPC_ccrhs + 401;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Phi16[17];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_W16[17];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ysd16[25];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lsd16[25];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_z17 = beliefPenaltyMPC_z + 289;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzaff17 = beliefPenaltyMPC_dz_aff + 289;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzcc17 = beliefPenaltyMPC_dz_cc + 289;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_rd17 = beliefPenaltyMPC_rd + 289;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lbyrd17[17];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_cost17 = beliefPenaltyMPC_grad_cost + 289;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_eq17 = beliefPenaltyMPC_grad_eq + 289;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_ineq17 = beliefPenaltyMPC_grad_ineq + 289;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ctv17[17];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_v17 = beliefPenaltyMPC_v + 90;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_re17[5];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_beta17[5];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_betacc17[5];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvaff17 = beliefPenaltyMPC_dv_aff + 90;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvcc17 = beliefPenaltyMPC_dv_cc + 90;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_V17[85];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Yd17[15];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ld17[15];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_yy17[5];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_bmy17[5];
int beliefPenaltyMPC_lbIdx17[17] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llb17 = beliefPenaltyMPC_l + 408;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_slb17 = beliefPenaltyMPC_s + 408;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llbbyslb17 = beliefPenaltyMPC_lbys + 408;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_rilb17[17];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbaff17 = beliefPenaltyMPC_dl_aff + 408;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbaff17 = beliefPenaltyMPC_ds_aff + 408;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbcc17 = beliefPenaltyMPC_dl_cc + 408;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbcc17 = beliefPenaltyMPC_ds_cc + 408;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsl17 = beliefPenaltyMPC_ccrhs + 408;
int beliefPenaltyMPC_ubIdx17[7] = {0, 1, 2, 3, 4, 5, 6};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lub17 = beliefPenaltyMPC_l + 425;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_sub17 = beliefPenaltyMPC_s + 425;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lubbysub17 = beliefPenaltyMPC_lbys + 425;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_riub17[7];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubaff17 = beliefPenaltyMPC_dl_aff + 425;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubaff17 = beliefPenaltyMPC_ds_aff + 425;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubcc17 = beliefPenaltyMPC_dl_cc + 425;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubcc17 = beliefPenaltyMPC_ds_cc + 425;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsub17 = beliefPenaltyMPC_ccrhs + 425;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Phi17[17];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_W17[17];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ysd17[25];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lsd17[25];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_z18 = beliefPenaltyMPC_z + 306;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzaff18 = beliefPenaltyMPC_dz_aff + 306;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzcc18 = beliefPenaltyMPC_dz_cc + 306;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_rd18 = beliefPenaltyMPC_rd + 306;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lbyrd18[17];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_cost18 = beliefPenaltyMPC_grad_cost + 306;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_eq18 = beliefPenaltyMPC_grad_eq + 306;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_ineq18 = beliefPenaltyMPC_grad_ineq + 306;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ctv18[17];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_v18 = beliefPenaltyMPC_v + 95;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_re18[5];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_beta18[5];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_betacc18[5];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvaff18 = beliefPenaltyMPC_dv_aff + 95;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvcc18 = beliefPenaltyMPC_dv_cc + 95;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_V18[85];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Yd18[15];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ld18[15];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_yy18[5];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_bmy18[5];
int beliefPenaltyMPC_lbIdx18[17] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llb18 = beliefPenaltyMPC_l + 432;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_slb18 = beliefPenaltyMPC_s + 432;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llbbyslb18 = beliefPenaltyMPC_lbys + 432;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_rilb18[17];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbaff18 = beliefPenaltyMPC_dl_aff + 432;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbaff18 = beliefPenaltyMPC_ds_aff + 432;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbcc18 = beliefPenaltyMPC_dl_cc + 432;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbcc18 = beliefPenaltyMPC_ds_cc + 432;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsl18 = beliefPenaltyMPC_ccrhs + 432;
int beliefPenaltyMPC_ubIdx18[7] = {0, 1, 2, 3, 4, 5, 6};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lub18 = beliefPenaltyMPC_l + 449;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_sub18 = beliefPenaltyMPC_s + 449;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lubbysub18 = beliefPenaltyMPC_lbys + 449;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_riub18[7];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubaff18 = beliefPenaltyMPC_dl_aff + 449;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubaff18 = beliefPenaltyMPC_ds_aff + 449;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubcc18 = beliefPenaltyMPC_dl_cc + 449;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubcc18 = beliefPenaltyMPC_ds_cc + 449;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsub18 = beliefPenaltyMPC_ccrhs + 449;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Phi18[17];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_W18[17];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ysd18[25];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lsd18[25];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_z19 = beliefPenaltyMPC_z + 323;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzaff19 = beliefPenaltyMPC_dz_aff + 323;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzcc19 = beliefPenaltyMPC_dz_cc + 323;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_rd19 = beliefPenaltyMPC_rd + 323;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lbyrd19[17];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_cost19 = beliefPenaltyMPC_grad_cost + 323;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_eq19 = beliefPenaltyMPC_grad_eq + 323;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_ineq19 = beliefPenaltyMPC_grad_ineq + 323;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ctv19[17];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_v19 = beliefPenaltyMPC_v + 100;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_re19[5];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_beta19[5];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_betacc19[5];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvaff19 = beliefPenaltyMPC_dv_aff + 100;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvcc19 = beliefPenaltyMPC_dv_cc + 100;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_V19[85];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Yd19[15];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ld19[15];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_yy19[5];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_bmy19[5];
int beliefPenaltyMPC_lbIdx19[17] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llb19 = beliefPenaltyMPC_l + 456;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_slb19 = beliefPenaltyMPC_s + 456;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llbbyslb19 = beliefPenaltyMPC_lbys + 456;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_rilb19[17];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbaff19 = beliefPenaltyMPC_dl_aff + 456;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbaff19 = beliefPenaltyMPC_ds_aff + 456;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbcc19 = beliefPenaltyMPC_dl_cc + 456;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbcc19 = beliefPenaltyMPC_ds_cc + 456;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsl19 = beliefPenaltyMPC_ccrhs + 456;
int beliefPenaltyMPC_ubIdx19[7] = {0, 1, 2, 3, 4, 5, 6};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lub19 = beliefPenaltyMPC_l + 473;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_sub19 = beliefPenaltyMPC_s + 473;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lubbysub19 = beliefPenaltyMPC_lbys + 473;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_riub19[7];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubaff19 = beliefPenaltyMPC_dl_aff + 473;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubaff19 = beliefPenaltyMPC_ds_aff + 473;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubcc19 = beliefPenaltyMPC_dl_cc + 473;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubcc19 = beliefPenaltyMPC_ds_cc + 473;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsub19 = beliefPenaltyMPC_ccrhs + 473;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Phi19[17];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_W19[17];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ysd19[25];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lsd19[25];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_z20 = beliefPenaltyMPC_z + 340;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzaff20 = beliefPenaltyMPC_dz_aff + 340;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzcc20 = beliefPenaltyMPC_dz_cc + 340;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_rd20 = beliefPenaltyMPC_rd + 340;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lbyrd20[17];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_cost20 = beliefPenaltyMPC_grad_cost + 340;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_eq20 = beliefPenaltyMPC_grad_eq + 340;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_ineq20 = beliefPenaltyMPC_grad_ineq + 340;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ctv20[17];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_v20 = beliefPenaltyMPC_v + 105;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_re20[5];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_beta20[5];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_betacc20[5];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvaff20 = beliefPenaltyMPC_dv_aff + 105;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvcc20 = beliefPenaltyMPC_dv_cc + 105;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_V20[85];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Yd20[15];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ld20[15];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_yy20[5];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_bmy20[5];
int beliefPenaltyMPC_lbIdx20[17] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llb20 = beliefPenaltyMPC_l + 480;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_slb20 = beliefPenaltyMPC_s + 480;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llbbyslb20 = beliefPenaltyMPC_lbys + 480;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_rilb20[17];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbaff20 = beliefPenaltyMPC_dl_aff + 480;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbaff20 = beliefPenaltyMPC_ds_aff + 480;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbcc20 = beliefPenaltyMPC_dl_cc + 480;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbcc20 = beliefPenaltyMPC_ds_cc + 480;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsl20 = beliefPenaltyMPC_ccrhs + 480;
int beliefPenaltyMPC_ubIdx20[7] = {0, 1, 2, 3, 4, 5, 6};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lub20 = beliefPenaltyMPC_l + 497;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_sub20 = beliefPenaltyMPC_s + 497;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lubbysub20 = beliefPenaltyMPC_lbys + 497;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_riub20[7];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubaff20 = beliefPenaltyMPC_dl_aff + 497;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubaff20 = beliefPenaltyMPC_ds_aff + 497;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubcc20 = beliefPenaltyMPC_dl_cc + 497;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubcc20 = beliefPenaltyMPC_ds_cc + 497;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsub20 = beliefPenaltyMPC_ccrhs + 497;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Phi20[17];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_W20[17];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ysd20[25];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lsd20[25];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_z21 = beliefPenaltyMPC_z + 357;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzaff21 = beliefPenaltyMPC_dz_aff + 357;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzcc21 = beliefPenaltyMPC_dz_cc + 357;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_rd21 = beliefPenaltyMPC_rd + 357;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lbyrd21[17];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_cost21 = beliefPenaltyMPC_grad_cost + 357;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_eq21 = beliefPenaltyMPC_grad_eq + 357;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_ineq21 = beliefPenaltyMPC_grad_ineq + 357;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ctv21[17];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_v21 = beliefPenaltyMPC_v + 110;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_re21[5];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_beta21[5];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_betacc21[5];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvaff21 = beliefPenaltyMPC_dv_aff + 110;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvcc21 = beliefPenaltyMPC_dv_cc + 110;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_V21[85];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Yd21[15];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ld21[15];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_yy21[5];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_bmy21[5];
int beliefPenaltyMPC_lbIdx21[17] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llb21 = beliefPenaltyMPC_l + 504;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_slb21 = beliefPenaltyMPC_s + 504;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llbbyslb21 = beliefPenaltyMPC_lbys + 504;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_rilb21[17];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbaff21 = beliefPenaltyMPC_dl_aff + 504;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbaff21 = beliefPenaltyMPC_ds_aff + 504;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbcc21 = beliefPenaltyMPC_dl_cc + 504;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbcc21 = beliefPenaltyMPC_ds_cc + 504;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsl21 = beliefPenaltyMPC_ccrhs + 504;
int beliefPenaltyMPC_ubIdx21[7] = {0, 1, 2, 3, 4, 5, 6};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lub21 = beliefPenaltyMPC_l + 521;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_sub21 = beliefPenaltyMPC_s + 521;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lubbysub21 = beliefPenaltyMPC_lbys + 521;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_riub21[7];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubaff21 = beliefPenaltyMPC_dl_aff + 521;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubaff21 = beliefPenaltyMPC_ds_aff + 521;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubcc21 = beliefPenaltyMPC_dl_cc + 521;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubcc21 = beliefPenaltyMPC_ds_cc + 521;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsub21 = beliefPenaltyMPC_ccrhs + 521;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Phi21[17];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_W21[17];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ysd21[25];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lsd21[25];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_z22 = beliefPenaltyMPC_z + 374;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzaff22 = beliefPenaltyMPC_dz_aff + 374;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzcc22 = beliefPenaltyMPC_dz_cc + 374;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_rd22 = beliefPenaltyMPC_rd + 374;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lbyrd22[17];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_cost22 = beliefPenaltyMPC_grad_cost + 374;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_eq22 = beliefPenaltyMPC_grad_eq + 374;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_ineq22 = beliefPenaltyMPC_grad_ineq + 374;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ctv22[17];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_v22 = beliefPenaltyMPC_v + 115;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_re22[5];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_beta22[5];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_betacc22[5];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvaff22 = beliefPenaltyMPC_dv_aff + 115;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvcc22 = beliefPenaltyMPC_dv_cc + 115;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_V22[85];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Yd22[15];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ld22[15];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_yy22[5];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_bmy22[5];
int beliefPenaltyMPC_lbIdx22[17] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llb22 = beliefPenaltyMPC_l + 528;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_slb22 = beliefPenaltyMPC_s + 528;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llbbyslb22 = beliefPenaltyMPC_lbys + 528;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_rilb22[17];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbaff22 = beliefPenaltyMPC_dl_aff + 528;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbaff22 = beliefPenaltyMPC_ds_aff + 528;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbcc22 = beliefPenaltyMPC_dl_cc + 528;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbcc22 = beliefPenaltyMPC_ds_cc + 528;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsl22 = beliefPenaltyMPC_ccrhs + 528;
int beliefPenaltyMPC_ubIdx22[7] = {0, 1, 2, 3, 4, 5, 6};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lub22 = beliefPenaltyMPC_l + 545;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_sub22 = beliefPenaltyMPC_s + 545;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lubbysub22 = beliefPenaltyMPC_lbys + 545;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_riub22[7];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubaff22 = beliefPenaltyMPC_dl_aff + 545;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubaff22 = beliefPenaltyMPC_ds_aff + 545;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubcc22 = beliefPenaltyMPC_dl_cc + 545;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubcc22 = beliefPenaltyMPC_ds_cc + 545;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsub22 = beliefPenaltyMPC_ccrhs + 545;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Phi22[17];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_W22[17];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ysd22[25];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lsd22[25];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_z23 = beliefPenaltyMPC_z + 391;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzaff23 = beliefPenaltyMPC_dz_aff + 391;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzcc23 = beliefPenaltyMPC_dz_cc + 391;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_rd23 = beliefPenaltyMPC_rd + 391;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lbyrd23[17];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_cost23 = beliefPenaltyMPC_grad_cost + 391;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_eq23 = beliefPenaltyMPC_grad_eq + 391;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_ineq23 = beliefPenaltyMPC_grad_ineq + 391;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ctv23[17];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_v23 = beliefPenaltyMPC_v + 120;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_re23[5];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_beta23[5];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_betacc23[5];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvaff23 = beliefPenaltyMPC_dv_aff + 120;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvcc23 = beliefPenaltyMPC_dv_cc + 120;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_V23[85];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Yd23[15];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ld23[15];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_yy23[5];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_bmy23[5];
int beliefPenaltyMPC_lbIdx23[17] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llb23 = beliefPenaltyMPC_l + 552;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_slb23 = beliefPenaltyMPC_s + 552;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llbbyslb23 = beliefPenaltyMPC_lbys + 552;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_rilb23[17];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbaff23 = beliefPenaltyMPC_dl_aff + 552;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbaff23 = beliefPenaltyMPC_ds_aff + 552;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbcc23 = beliefPenaltyMPC_dl_cc + 552;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbcc23 = beliefPenaltyMPC_ds_cc + 552;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsl23 = beliefPenaltyMPC_ccrhs + 552;
int beliefPenaltyMPC_ubIdx23[7] = {0, 1, 2, 3, 4, 5, 6};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lub23 = beliefPenaltyMPC_l + 569;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_sub23 = beliefPenaltyMPC_s + 569;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lubbysub23 = beliefPenaltyMPC_lbys + 569;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_riub23[7];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubaff23 = beliefPenaltyMPC_dl_aff + 569;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubaff23 = beliefPenaltyMPC_ds_aff + 569;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubcc23 = beliefPenaltyMPC_dl_cc + 569;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubcc23 = beliefPenaltyMPC_ds_cc + 569;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsub23 = beliefPenaltyMPC_ccrhs + 569;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Phi23[17];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_W23[17];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ysd23[25];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lsd23[25];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_H24[5] = {0.0000000000000000E+000, 0.0000000000000000E+000, 2.0000000000000000E+001, 2.0000000000000000E+001, 2.0000000000000000E+001};
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_f24[5] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_z24 = beliefPenaltyMPC_z + 408;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzaff24 = beliefPenaltyMPC_dz_aff + 408;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzcc24 = beliefPenaltyMPC_dz_cc + 408;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_rd24 = beliefPenaltyMPC_rd + 408;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lbyrd24[5];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_cost24 = beliefPenaltyMPC_grad_cost + 408;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_eq24 = beliefPenaltyMPC_grad_eq + 408;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_ineq24 = beliefPenaltyMPC_grad_ineq + 408;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ctv24[5];
int beliefPenaltyMPC_lbIdx24[5] = {0, 1, 2, 3, 4};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llb24 = beliefPenaltyMPC_l + 576;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_slb24 = beliefPenaltyMPC_s + 576;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llbbyslb24 = beliefPenaltyMPC_lbys + 576;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_rilb24[5];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbaff24 = beliefPenaltyMPC_dl_aff + 576;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbaff24 = beliefPenaltyMPC_ds_aff + 576;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbcc24 = beliefPenaltyMPC_dl_cc + 576;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbcc24 = beliefPenaltyMPC_ds_cc + 576;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsl24 = beliefPenaltyMPC_ccrhs + 576;
int beliefPenaltyMPC_ubIdx24[5] = {0, 1, 2, 3, 4};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lub24 = beliefPenaltyMPC_l + 581;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_sub24 = beliefPenaltyMPC_s + 581;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lubbysub24 = beliefPenaltyMPC_lbys + 581;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_riub24[5];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubaff24 = beliefPenaltyMPC_dl_aff + 581;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubaff24 = beliefPenaltyMPC_ds_aff + 581;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubcc24 = beliefPenaltyMPC_dl_cc + 581;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubcc24 = beliefPenaltyMPC_ds_cc + 581;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsub24 = beliefPenaltyMPC_ccrhs + 581;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Phi24[5];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_D24[5] = {-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000};
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_W24[5];
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
beliefPenaltyMPC_LA_INITIALIZEVECTOR_413(beliefPenaltyMPC_z, 0);
beliefPenaltyMPC_LA_INITIALIZEVECTOR_125(beliefPenaltyMPC_v, 1);
beliefPenaltyMPC_LA_INITIALIZEVECTOR_586(beliefPenaltyMPC_l, 1);
beliefPenaltyMPC_LA_INITIALIZEVECTOR_586(beliefPenaltyMPC_s, 1);
info->mu = 0;
beliefPenaltyMPC_LA_DOTACC_586(beliefPenaltyMPC_l, beliefPenaltyMPC_s, &info->mu);
info->mu /= 586;
while( 1 ){
info->pobj = 0;
beliefPenaltyMPC_LA_DIAG_QUADFCN_17(beliefPenaltyMPC_H00, params->f1, beliefPenaltyMPC_z00, beliefPenaltyMPC_grad_cost00, &info->pobj);
beliefPenaltyMPC_LA_DIAG_QUADFCN_17(beliefPenaltyMPC_H00, params->f2, beliefPenaltyMPC_z01, beliefPenaltyMPC_grad_cost01, &info->pobj);
beliefPenaltyMPC_LA_DIAG_QUADFCN_17(beliefPenaltyMPC_H00, params->f3, beliefPenaltyMPC_z02, beliefPenaltyMPC_grad_cost02, &info->pobj);
beliefPenaltyMPC_LA_DIAG_QUADFCN_17(beliefPenaltyMPC_H00, params->f4, beliefPenaltyMPC_z03, beliefPenaltyMPC_grad_cost03, &info->pobj);
beliefPenaltyMPC_LA_DIAG_QUADFCN_17(beliefPenaltyMPC_H00, params->f5, beliefPenaltyMPC_z04, beliefPenaltyMPC_grad_cost04, &info->pobj);
beliefPenaltyMPC_LA_DIAG_QUADFCN_17(beliefPenaltyMPC_H00, params->f6, beliefPenaltyMPC_z05, beliefPenaltyMPC_grad_cost05, &info->pobj);
beliefPenaltyMPC_LA_DIAG_QUADFCN_17(beliefPenaltyMPC_H00, params->f7, beliefPenaltyMPC_z06, beliefPenaltyMPC_grad_cost06, &info->pobj);
beliefPenaltyMPC_LA_DIAG_QUADFCN_17(beliefPenaltyMPC_H00, params->f8, beliefPenaltyMPC_z07, beliefPenaltyMPC_grad_cost07, &info->pobj);
beliefPenaltyMPC_LA_DIAG_QUADFCN_17(beliefPenaltyMPC_H00, params->f9, beliefPenaltyMPC_z08, beliefPenaltyMPC_grad_cost08, &info->pobj);
beliefPenaltyMPC_LA_DIAG_QUADFCN_17(beliefPenaltyMPC_H00, params->f10, beliefPenaltyMPC_z09, beliefPenaltyMPC_grad_cost09, &info->pobj);
beliefPenaltyMPC_LA_DIAG_QUADFCN_17(beliefPenaltyMPC_H00, params->f11, beliefPenaltyMPC_z10, beliefPenaltyMPC_grad_cost10, &info->pobj);
beliefPenaltyMPC_LA_DIAG_QUADFCN_17(beliefPenaltyMPC_H00, params->f12, beliefPenaltyMPC_z11, beliefPenaltyMPC_grad_cost11, &info->pobj);
beliefPenaltyMPC_LA_DIAG_QUADFCN_17(beliefPenaltyMPC_H00, params->f13, beliefPenaltyMPC_z12, beliefPenaltyMPC_grad_cost12, &info->pobj);
beliefPenaltyMPC_LA_DIAG_QUADFCN_17(beliefPenaltyMPC_H00, params->f14, beliefPenaltyMPC_z13, beliefPenaltyMPC_grad_cost13, &info->pobj);
beliefPenaltyMPC_LA_DIAG_QUADFCN_17(beliefPenaltyMPC_H00, params->f15, beliefPenaltyMPC_z14, beliefPenaltyMPC_grad_cost14, &info->pobj);
beliefPenaltyMPC_LA_DIAG_QUADFCN_17(beliefPenaltyMPC_H00, params->f16, beliefPenaltyMPC_z15, beliefPenaltyMPC_grad_cost15, &info->pobj);
beliefPenaltyMPC_LA_DIAG_QUADFCN_17(beliefPenaltyMPC_H00, params->f17, beliefPenaltyMPC_z16, beliefPenaltyMPC_grad_cost16, &info->pobj);
beliefPenaltyMPC_LA_DIAG_QUADFCN_17(beliefPenaltyMPC_H00, params->f18, beliefPenaltyMPC_z17, beliefPenaltyMPC_grad_cost17, &info->pobj);
beliefPenaltyMPC_LA_DIAG_QUADFCN_17(beliefPenaltyMPC_H00, params->f19, beliefPenaltyMPC_z18, beliefPenaltyMPC_grad_cost18, &info->pobj);
beliefPenaltyMPC_LA_DIAG_QUADFCN_17(beliefPenaltyMPC_H00, params->f20, beliefPenaltyMPC_z19, beliefPenaltyMPC_grad_cost19, &info->pobj);
beliefPenaltyMPC_LA_DIAG_QUADFCN_17(beliefPenaltyMPC_H00, params->f21, beliefPenaltyMPC_z20, beliefPenaltyMPC_grad_cost20, &info->pobj);
beliefPenaltyMPC_LA_DIAG_QUADFCN_17(beliefPenaltyMPC_H00, params->f22, beliefPenaltyMPC_z21, beliefPenaltyMPC_grad_cost21, &info->pobj);
beliefPenaltyMPC_LA_DIAG_QUADFCN_17(beliefPenaltyMPC_H00, params->f23, beliefPenaltyMPC_z22, beliefPenaltyMPC_grad_cost22, &info->pobj);
beliefPenaltyMPC_LA_DIAG_QUADFCN_17(beliefPenaltyMPC_H00, params->f24, beliefPenaltyMPC_z23, beliefPenaltyMPC_grad_cost23, &info->pobj);
beliefPenaltyMPC_LA_DIAG_QUADFCN_5(beliefPenaltyMPC_H24, beliefPenaltyMPC_f24, beliefPenaltyMPC_z24, beliefPenaltyMPC_grad_cost24, &info->pobj);
info->res_eq = 0;
info->dgap = 0;
beliefPenaltyMPC_LA_DENSE_MVMSUB3_10_17_17(params->C1, beliefPenaltyMPC_z00, beliefPenaltyMPC_D01, beliefPenaltyMPC_z01, params->e1, beliefPenaltyMPC_v00, beliefPenaltyMPC_re00, &info->dgap, &info->res_eq);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_5_17_17(params->C2, beliefPenaltyMPC_z01, beliefPenaltyMPC_D02, beliefPenaltyMPC_z02, params->e2, beliefPenaltyMPC_v01, beliefPenaltyMPC_re01, &info->dgap, &info->res_eq);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_5_17_17(params->C3, beliefPenaltyMPC_z02, beliefPenaltyMPC_D02, beliefPenaltyMPC_z03, params->e3, beliefPenaltyMPC_v02, beliefPenaltyMPC_re02, &info->dgap, &info->res_eq);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_5_17_17(params->C4, beliefPenaltyMPC_z03, beliefPenaltyMPC_D02, beliefPenaltyMPC_z04, params->e4, beliefPenaltyMPC_v03, beliefPenaltyMPC_re03, &info->dgap, &info->res_eq);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_5_17_17(params->C5, beliefPenaltyMPC_z04, beliefPenaltyMPC_D02, beliefPenaltyMPC_z05, params->e5, beliefPenaltyMPC_v04, beliefPenaltyMPC_re04, &info->dgap, &info->res_eq);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_5_17_17(params->C6, beliefPenaltyMPC_z05, beliefPenaltyMPC_D02, beliefPenaltyMPC_z06, params->e6, beliefPenaltyMPC_v05, beliefPenaltyMPC_re05, &info->dgap, &info->res_eq);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_5_17_17(params->C7, beliefPenaltyMPC_z06, beliefPenaltyMPC_D02, beliefPenaltyMPC_z07, params->e7, beliefPenaltyMPC_v06, beliefPenaltyMPC_re06, &info->dgap, &info->res_eq);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_5_17_17(params->C8, beliefPenaltyMPC_z07, beliefPenaltyMPC_D02, beliefPenaltyMPC_z08, params->e8, beliefPenaltyMPC_v07, beliefPenaltyMPC_re07, &info->dgap, &info->res_eq);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_5_17_17(params->C9, beliefPenaltyMPC_z08, beliefPenaltyMPC_D02, beliefPenaltyMPC_z09, params->e9, beliefPenaltyMPC_v08, beliefPenaltyMPC_re08, &info->dgap, &info->res_eq);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_5_17_17(params->C10, beliefPenaltyMPC_z09, beliefPenaltyMPC_D02, beliefPenaltyMPC_z10, params->e10, beliefPenaltyMPC_v09, beliefPenaltyMPC_re09, &info->dgap, &info->res_eq);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_5_17_17(params->C11, beliefPenaltyMPC_z10, beliefPenaltyMPC_D02, beliefPenaltyMPC_z11, params->e11, beliefPenaltyMPC_v10, beliefPenaltyMPC_re10, &info->dgap, &info->res_eq);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_5_17_17(params->C12, beliefPenaltyMPC_z11, beliefPenaltyMPC_D02, beliefPenaltyMPC_z12, params->e12, beliefPenaltyMPC_v11, beliefPenaltyMPC_re11, &info->dgap, &info->res_eq);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_5_17_17(params->C13, beliefPenaltyMPC_z12, beliefPenaltyMPC_D02, beliefPenaltyMPC_z13, params->e13, beliefPenaltyMPC_v12, beliefPenaltyMPC_re12, &info->dgap, &info->res_eq);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_5_17_17(params->C14, beliefPenaltyMPC_z13, beliefPenaltyMPC_D02, beliefPenaltyMPC_z14, params->e14, beliefPenaltyMPC_v13, beliefPenaltyMPC_re13, &info->dgap, &info->res_eq);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_5_17_17(params->C15, beliefPenaltyMPC_z14, beliefPenaltyMPC_D02, beliefPenaltyMPC_z15, params->e15, beliefPenaltyMPC_v14, beliefPenaltyMPC_re14, &info->dgap, &info->res_eq);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_5_17_17(params->C16, beliefPenaltyMPC_z15, beliefPenaltyMPC_D02, beliefPenaltyMPC_z16, params->e16, beliefPenaltyMPC_v15, beliefPenaltyMPC_re15, &info->dgap, &info->res_eq);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_5_17_17(params->C17, beliefPenaltyMPC_z16, beliefPenaltyMPC_D02, beliefPenaltyMPC_z17, params->e17, beliefPenaltyMPC_v16, beliefPenaltyMPC_re16, &info->dgap, &info->res_eq);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_5_17_17(params->C18, beliefPenaltyMPC_z17, beliefPenaltyMPC_D02, beliefPenaltyMPC_z18, params->e18, beliefPenaltyMPC_v17, beliefPenaltyMPC_re17, &info->dgap, &info->res_eq);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_5_17_17(params->C19, beliefPenaltyMPC_z18, beliefPenaltyMPC_D02, beliefPenaltyMPC_z19, params->e19, beliefPenaltyMPC_v18, beliefPenaltyMPC_re18, &info->dgap, &info->res_eq);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_5_17_17(params->C20, beliefPenaltyMPC_z19, beliefPenaltyMPC_D02, beliefPenaltyMPC_z20, params->e20, beliefPenaltyMPC_v19, beliefPenaltyMPC_re19, &info->dgap, &info->res_eq);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_5_17_17(params->C21, beliefPenaltyMPC_z20, beliefPenaltyMPC_D02, beliefPenaltyMPC_z21, params->e21, beliefPenaltyMPC_v20, beliefPenaltyMPC_re20, &info->dgap, &info->res_eq);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_5_17_17(params->C22, beliefPenaltyMPC_z21, beliefPenaltyMPC_D02, beliefPenaltyMPC_z22, params->e22, beliefPenaltyMPC_v21, beliefPenaltyMPC_re21, &info->dgap, &info->res_eq);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_5_17_17(params->C23, beliefPenaltyMPC_z22, beliefPenaltyMPC_D02, beliefPenaltyMPC_z23, params->e23, beliefPenaltyMPC_v22, beliefPenaltyMPC_re22, &info->dgap, &info->res_eq);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_5_17_5(params->C24, beliefPenaltyMPC_z23, beliefPenaltyMPC_D24, beliefPenaltyMPC_z24, params->e24, beliefPenaltyMPC_v23, beliefPenaltyMPC_re23, &info->dgap, &info->res_eq);
beliefPenaltyMPC_LA_DENSE_MTVM_10_17(params->C1, beliefPenaltyMPC_v00, beliefPenaltyMPC_grad_eq00);
beliefPenaltyMPC_LA_DENSE_MTVM2_5_17_10(params->C2, beliefPenaltyMPC_v01, beliefPenaltyMPC_D01, beliefPenaltyMPC_v00, beliefPenaltyMPC_grad_eq01);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_5_17_5(params->C3, beliefPenaltyMPC_v02, beliefPenaltyMPC_D02, beliefPenaltyMPC_v01, beliefPenaltyMPC_grad_eq02);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_5_17_5(params->C4, beliefPenaltyMPC_v03, beliefPenaltyMPC_D02, beliefPenaltyMPC_v02, beliefPenaltyMPC_grad_eq03);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_5_17_5(params->C5, beliefPenaltyMPC_v04, beliefPenaltyMPC_D02, beliefPenaltyMPC_v03, beliefPenaltyMPC_grad_eq04);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_5_17_5(params->C6, beliefPenaltyMPC_v05, beliefPenaltyMPC_D02, beliefPenaltyMPC_v04, beliefPenaltyMPC_grad_eq05);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_5_17_5(params->C7, beliefPenaltyMPC_v06, beliefPenaltyMPC_D02, beliefPenaltyMPC_v05, beliefPenaltyMPC_grad_eq06);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_5_17_5(params->C8, beliefPenaltyMPC_v07, beliefPenaltyMPC_D02, beliefPenaltyMPC_v06, beliefPenaltyMPC_grad_eq07);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_5_17_5(params->C9, beliefPenaltyMPC_v08, beliefPenaltyMPC_D02, beliefPenaltyMPC_v07, beliefPenaltyMPC_grad_eq08);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_5_17_5(params->C10, beliefPenaltyMPC_v09, beliefPenaltyMPC_D02, beliefPenaltyMPC_v08, beliefPenaltyMPC_grad_eq09);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_5_17_5(params->C11, beliefPenaltyMPC_v10, beliefPenaltyMPC_D02, beliefPenaltyMPC_v09, beliefPenaltyMPC_grad_eq10);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_5_17_5(params->C12, beliefPenaltyMPC_v11, beliefPenaltyMPC_D02, beliefPenaltyMPC_v10, beliefPenaltyMPC_grad_eq11);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_5_17_5(params->C13, beliefPenaltyMPC_v12, beliefPenaltyMPC_D02, beliefPenaltyMPC_v11, beliefPenaltyMPC_grad_eq12);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_5_17_5(params->C14, beliefPenaltyMPC_v13, beliefPenaltyMPC_D02, beliefPenaltyMPC_v12, beliefPenaltyMPC_grad_eq13);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_5_17_5(params->C15, beliefPenaltyMPC_v14, beliefPenaltyMPC_D02, beliefPenaltyMPC_v13, beliefPenaltyMPC_grad_eq14);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_5_17_5(params->C16, beliefPenaltyMPC_v15, beliefPenaltyMPC_D02, beliefPenaltyMPC_v14, beliefPenaltyMPC_grad_eq15);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_5_17_5(params->C17, beliefPenaltyMPC_v16, beliefPenaltyMPC_D02, beliefPenaltyMPC_v15, beliefPenaltyMPC_grad_eq16);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_5_17_5(params->C18, beliefPenaltyMPC_v17, beliefPenaltyMPC_D02, beliefPenaltyMPC_v16, beliefPenaltyMPC_grad_eq17);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_5_17_5(params->C19, beliefPenaltyMPC_v18, beliefPenaltyMPC_D02, beliefPenaltyMPC_v17, beliefPenaltyMPC_grad_eq18);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_5_17_5(params->C20, beliefPenaltyMPC_v19, beliefPenaltyMPC_D02, beliefPenaltyMPC_v18, beliefPenaltyMPC_grad_eq19);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_5_17_5(params->C21, beliefPenaltyMPC_v20, beliefPenaltyMPC_D02, beliefPenaltyMPC_v19, beliefPenaltyMPC_grad_eq20);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_5_17_5(params->C22, beliefPenaltyMPC_v21, beliefPenaltyMPC_D02, beliefPenaltyMPC_v20, beliefPenaltyMPC_grad_eq21);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_5_17_5(params->C23, beliefPenaltyMPC_v22, beliefPenaltyMPC_D02, beliefPenaltyMPC_v21, beliefPenaltyMPC_grad_eq22);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_5_17_5(params->C24, beliefPenaltyMPC_v23, beliefPenaltyMPC_D02, beliefPenaltyMPC_v22, beliefPenaltyMPC_grad_eq23);
beliefPenaltyMPC_LA_DIAGZERO_MTVM_5_5(beliefPenaltyMPC_D24, beliefPenaltyMPC_v23, beliefPenaltyMPC_grad_eq24);
info->res_ineq = 0;
beliefPenaltyMPC_LA_VSUBADD3_17(params->lb1, beliefPenaltyMPC_z00, beliefPenaltyMPC_lbIdx00, beliefPenaltyMPC_llb00, beliefPenaltyMPC_slb00, beliefPenaltyMPC_rilb00, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD2_7(beliefPenaltyMPC_z00, beliefPenaltyMPC_ubIdx00, params->ub1, beliefPenaltyMPC_lub00, beliefPenaltyMPC_sub00, beliefPenaltyMPC_riub00, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD3_17(params->lb2, beliefPenaltyMPC_z01, beliefPenaltyMPC_lbIdx01, beliefPenaltyMPC_llb01, beliefPenaltyMPC_slb01, beliefPenaltyMPC_rilb01, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD2_7(beliefPenaltyMPC_z01, beliefPenaltyMPC_ubIdx01, params->ub2, beliefPenaltyMPC_lub01, beliefPenaltyMPC_sub01, beliefPenaltyMPC_riub01, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD3_17(params->lb3, beliefPenaltyMPC_z02, beliefPenaltyMPC_lbIdx02, beliefPenaltyMPC_llb02, beliefPenaltyMPC_slb02, beliefPenaltyMPC_rilb02, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD2_7(beliefPenaltyMPC_z02, beliefPenaltyMPC_ubIdx02, params->ub3, beliefPenaltyMPC_lub02, beliefPenaltyMPC_sub02, beliefPenaltyMPC_riub02, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD3_17(params->lb4, beliefPenaltyMPC_z03, beliefPenaltyMPC_lbIdx03, beliefPenaltyMPC_llb03, beliefPenaltyMPC_slb03, beliefPenaltyMPC_rilb03, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD2_7(beliefPenaltyMPC_z03, beliefPenaltyMPC_ubIdx03, params->ub4, beliefPenaltyMPC_lub03, beliefPenaltyMPC_sub03, beliefPenaltyMPC_riub03, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD3_17(params->lb5, beliefPenaltyMPC_z04, beliefPenaltyMPC_lbIdx04, beliefPenaltyMPC_llb04, beliefPenaltyMPC_slb04, beliefPenaltyMPC_rilb04, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD2_7(beliefPenaltyMPC_z04, beliefPenaltyMPC_ubIdx04, params->ub5, beliefPenaltyMPC_lub04, beliefPenaltyMPC_sub04, beliefPenaltyMPC_riub04, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD3_17(params->lb6, beliefPenaltyMPC_z05, beliefPenaltyMPC_lbIdx05, beliefPenaltyMPC_llb05, beliefPenaltyMPC_slb05, beliefPenaltyMPC_rilb05, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD2_7(beliefPenaltyMPC_z05, beliefPenaltyMPC_ubIdx05, params->ub6, beliefPenaltyMPC_lub05, beliefPenaltyMPC_sub05, beliefPenaltyMPC_riub05, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD3_17(params->lb7, beliefPenaltyMPC_z06, beliefPenaltyMPC_lbIdx06, beliefPenaltyMPC_llb06, beliefPenaltyMPC_slb06, beliefPenaltyMPC_rilb06, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD2_7(beliefPenaltyMPC_z06, beliefPenaltyMPC_ubIdx06, params->ub7, beliefPenaltyMPC_lub06, beliefPenaltyMPC_sub06, beliefPenaltyMPC_riub06, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD3_17(params->lb8, beliefPenaltyMPC_z07, beliefPenaltyMPC_lbIdx07, beliefPenaltyMPC_llb07, beliefPenaltyMPC_slb07, beliefPenaltyMPC_rilb07, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD2_7(beliefPenaltyMPC_z07, beliefPenaltyMPC_ubIdx07, params->ub8, beliefPenaltyMPC_lub07, beliefPenaltyMPC_sub07, beliefPenaltyMPC_riub07, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD3_17(params->lb9, beliefPenaltyMPC_z08, beliefPenaltyMPC_lbIdx08, beliefPenaltyMPC_llb08, beliefPenaltyMPC_slb08, beliefPenaltyMPC_rilb08, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD2_7(beliefPenaltyMPC_z08, beliefPenaltyMPC_ubIdx08, params->ub9, beliefPenaltyMPC_lub08, beliefPenaltyMPC_sub08, beliefPenaltyMPC_riub08, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD3_17(params->lb10, beliefPenaltyMPC_z09, beliefPenaltyMPC_lbIdx09, beliefPenaltyMPC_llb09, beliefPenaltyMPC_slb09, beliefPenaltyMPC_rilb09, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD2_7(beliefPenaltyMPC_z09, beliefPenaltyMPC_ubIdx09, params->ub10, beliefPenaltyMPC_lub09, beliefPenaltyMPC_sub09, beliefPenaltyMPC_riub09, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD3_17(params->lb11, beliefPenaltyMPC_z10, beliefPenaltyMPC_lbIdx10, beliefPenaltyMPC_llb10, beliefPenaltyMPC_slb10, beliefPenaltyMPC_rilb10, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD2_7(beliefPenaltyMPC_z10, beliefPenaltyMPC_ubIdx10, params->ub11, beliefPenaltyMPC_lub10, beliefPenaltyMPC_sub10, beliefPenaltyMPC_riub10, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD3_17(params->lb12, beliefPenaltyMPC_z11, beliefPenaltyMPC_lbIdx11, beliefPenaltyMPC_llb11, beliefPenaltyMPC_slb11, beliefPenaltyMPC_rilb11, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD2_7(beliefPenaltyMPC_z11, beliefPenaltyMPC_ubIdx11, params->ub12, beliefPenaltyMPC_lub11, beliefPenaltyMPC_sub11, beliefPenaltyMPC_riub11, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD3_17(params->lb13, beliefPenaltyMPC_z12, beliefPenaltyMPC_lbIdx12, beliefPenaltyMPC_llb12, beliefPenaltyMPC_slb12, beliefPenaltyMPC_rilb12, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD2_7(beliefPenaltyMPC_z12, beliefPenaltyMPC_ubIdx12, params->ub13, beliefPenaltyMPC_lub12, beliefPenaltyMPC_sub12, beliefPenaltyMPC_riub12, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD3_17(params->lb14, beliefPenaltyMPC_z13, beliefPenaltyMPC_lbIdx13, beliefPenaltyMPC_llb13, beliefPenaltyMPC_slb13, beliefPenaltyMPC_rilb13, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD2_7(beliefPenaltyMPC_z13, beliefPenaltyMPC_ubIdx13, params->ub14, beliefPenaltyMPC_lub13, beliefPenaltyMPC_sub13, beliefPenaltyMPC_riub13, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD3_17(params->lb15, beliefPenaltyMPC_z14, beliefPenaltyMPC_lbIdx14, beliefPenaltyMPC_llb14, beliefPenaltyMPC_slb14, beliefPenaltyMPC_rilb14, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD2_7(beliefPenaltyMPC_z14, beliefPenaltyMPC_ubIdx14, params->ub15, beliefPenaltyMPC_lub14, beliefPenaltyMPC_sub14, beliefPenaltyMPC_riub14, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD3_17(params->lb16, beliefPenaltyMPC_z15, beliefPenaltyMPC_lbIdx15, beliefPenaltyMPC_llb15, beliefPenaltyMPC_slb15, beliefPenaltyMPC_rilb15, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD2_7(beliefPenaltyMPC_z15, beliefPenaltyMPC_ubIdx15, params->ub16, beliefPenaltyMPC_lub15, beliefPenaltyMPC_sub15, beliefPenaltyMPC_riub15, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD3_17(params->lb17, beliefPenaltyMPC_z16, beliefPenaltyMPC_lbIdx16, beliefPenaltyMPC_llb16, beliefPenaltyMPC_slb16, beliefPenaltyMPC_rilb16, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD2_7(beliefPenaltyMPC_z16, beliefPenaltyMPC_ubIdx16, params->ub17, beliefPenaltyMPC_lub16, beliefPenaltyMPC_sub16, beliefPenaltyMPC_riub16, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD3_17(params->lb18, beliefPenaltyMPC_z17, beliefPenaltyMPC_lbIdx17, beliefPenaltyMPC_llb17, beliefPenaltyMPC_slb17, beliefPenaltyMPC_rilb17, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD2_7(beliefPenaltyMPC_z17, beliefPenaltyMPC_ubIdx17, params->ub18, beliefPenaltyMPC_lub17, beliefPenaltyMPC_sub17, beliefPenaltyMPC_riub17, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD3_17(params->lb19, beliefPenaltyMPC_z18, beliefPenaltyMPC_lbIdx18, beliefPenaltyMPC_llb18, beliefPenaltyMPC_slb18, beliefPenaltyMPC_rilb18, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD2_7(beliefPenaltyMPC_z18, beliefPenaltyMPC_ubIdx18, params->ub19, beliefPenaltyMPC_lub18, beliefPenaltyMPC_sub18, beliefPenaltyMPC_riub18, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD3_17(params->lb20, beliefPenaltyMPC_z19, beliefPenaltyMPC_lbIdx19, beliefPenaltyMPC_llb19, beliefPenaltyMPC_slb19, beliefPenaltyMPC_rilb19, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD2_7(beliefPenaltyMPC_z19, beliefPenaltyMPC_ubIdx19, params->ub20, beliefPenaltyMPC_lub19, beliefPenaltyMPC_sub19, beliefPenaltyMPC_riub19, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD3_17(params->lb21, beliefPenaltyMPC_z20, beliefPenaltyMPC_lbIdx20, beliefPenaltyMPC_llb20, beliefPenaltyMPC_slb20, beliefPenaltyMPC_rilb20, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD2_7(beliefPenaltyMPC_z20, beliefPenaltyMPC_ubIdx20, params->ub21, beliefPenaltyMPC_lub20, beliefPenaltyMPC_sub20, beliefPenaltyMPC_riub20, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD3_17(params->lb22, beliefPenaltyMPC_z21, beliefPenaltyMPC_lbIdx21, beliefPenaltyMPC_llb21, beliefPenaltyMPC_slb21, beliefPenaltyMPC_rilb21, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD2_7(beliefPenaltyMPC_z21, beliefPenaltyMPC_ubIdx21, params->ub22, beliefPenaltyMPC_lub21, beliefPenaltyMPC_sub21, beliefPenaltyMPC_riub21, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD3_17(params->lb23, beliefPenaltyMPC_z22, beliefPenaltyMPC_lbIdx22, beliefPenaltyMPC_llb22, beliefPenaltyMPC_slb22, beliefPenaltyMPC_rilb22, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD2_7(beliefPenaltyMPC_z22, beliefPenaltyMPC_ubIdx22, params->ub23, beliefPenaltyMPC_lub22, beliefPenaltyMPC_sub22, beliefPenaltyMPC_riub22, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD3_17(params->lb24, beliefPenaltyMPC_z23, beliefPenaltyMPC_lbIdx23, beliefPenaltyMPC_llb23, beliefPenaltyMPC_slb23, beliefPenaltyMPC_rilb23, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD2_7(beliefPenaltyMPC_z23, beliefPenaltyMPC_ubIdx23, params->ub24, beliefPenaltyMPC_lub23, beliefPenaltyMPC_sub23, beliefPenaltyMPC_riub23, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD3_5(params->lb25, beliefPenaltyMPC_z24, beliefPenaltyMPC_lbIdx24, beliefPenaltyMPC_llb24, beliefPenaltyMPC_slb24, beliefPenaltyMPC_rilb24, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD2_5(beliefPenaltyMPC_z24, beliefPenaltyMPC_ubIdx24, params->ub25, beliefPenaltyMPC_lub24, beliefPenaltyMPC_sub24, beliefPenaltyMPC_riub24, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_INEQ_B_GRAD_17_17_7(beliefPenaltyMPC_lub00, beliefPenaltyMPC_sub00, beliefPenaltyMPC_riub00, beliefPenaltyMPC_llb00, beliefPenaltyMPC_slb00, beliefPenaltyMPC_rilb00, beliefPenaltyMPC_lbIdx00, beliefPenaltyMPC_ubIdx00, beliefPenaltyMPC_grad_ineq00, beliefPenaltyMPC_lubbysub00, beliefPenaltyMPC_llbbyslb00);
beliefPenaltyMPC_LA_INEQ_B_GRAD_17_17_7(beliefPenaltyMPC_lub01, beliefPenaltyMPC_sub01, beliefPenaltyMPC_riub01, beliefPenaltyMPC_llb01, beliefPenaltyMPC_slb01, beliefPenaltyMPC_rilb01, beliefPenaltyMPC_lbIdx01, beliefPenaltyMPC_ubIdx01, beliefPenaltyMPC_grad_ineq01, beliefPenaltyMPC_lubbysub01, beliefPenaltyMPC_llbbyslb01);
beliefPenaltyMPC_LA_INEQ_B_GRAD_17_17_7(beliefPenaltyMPC_lub02, beliefPenaltyMPC_sub02, beliefPenaltyMPC_riub02, beliefPenaltyMPC_llb02, beliefPenaltyMPC_slb02, beliefPenaltyMPC_rilb02, beliefPenaltyMPC_lbIdx02, beliefPenaltyMPC_ubIdx02, beliefPenaltyMPC_grad_ineq02, beliefPenaltyMPC_lubbysub02, beliefPenaltyMPC_llbbyslb02);
beliefPenaltyMPC_LA_INEQ_B_GRAD_17_17_7(beliefPenaltyMPC_lub03, beliefPenaltyMPC_sub03, beliefPenaltyMPC_riub03, beliefPenaltyMPC_llb03, beliefPenaltyMPC_slb03, beliefPenaltyMPC_rilb03, beliefPenaltyMPC_lbIdx03, beliefPenaltyMPC_ubIdx03, beliefPenaltyMPC_grad_ineq03, beliefPenaltyMPC_lubbysub03, beliefPenaltyMPC_llbbyslb03);
beliefPenaltyMPC_LA_INEQ_B_GRAD_17_17_7(beliefPenaltyMPC_lub04, beliefPenaltyMPC_sub04, beliefPenaltyMPC_riub04, beliefPenaltyMPC_llb04, beliefPenaltyMPC_slb04, beliefPenaltyMPC_rilb04, beliefPenaltyMPC_lbIdx04, beliefPenaltyMPC_ubIdx04, beliefPenaltyMPC_grad_ineq04, beliefPenaltyMPC_lubbysub04, beliefPenaltyMPC_llbbyslb04);
beliefPenaltyMPC_LA_INEQ_B_GRAD_17_17_7(beliefPenaltyMPC_lub05, beliefPenaltyMPC_sub05, beliefPenaltyMPC_riub05, beliefPenaltyMPC_llb05, beliefPenaltyMPC_slb05, beliefPenaltyMPC_rilb05, beliefPenaltyMPC_lbIdx05, beliefPenaltyMPC_ubIdx05, beliefPenaltyMPC_grad_ineq05, beliefPenaltyMPC_lubbysub05, beliefPenaltyMPC_llbbyslb05);
beliefPenaltyMPC_LA_INEQ_B_GRAD_17_17_7(beliefPenaltyMPC_lub06, beliefPenaltyMPC_sub06, beliefPenaltyMPC_riub06, beliefPenaltyMPC_llb06, beliefPenaltyMPC_slb06, beliefPenaltyMPC_rilb06, beliefPenaltyMPC_lbIdx06, beliefPenaltyMPC_ubIdx06, beliefPenaltyMPC_grad_ineq06, beliefPenaltyMPC_lubbysub06, beliefPenaltyMPC_llbbyslb06);
beliefPenaltyMPC_LA_INEQ_B_GRAD_17_17_7(beliefPenaltyMPC_lub07, beliefPenaltyMPC_sub07, beliefPenaltyMPC_riub07, beliefPenaltyMPC_llb07, beliefPenaltyMPC_slb07, beliefPenaltyMPC_rilb07, beliefPenaltyMPC_lbIdx07, beliefPenaltyMPC_ubIdx07, beliefPenaltyMPC_grad_ineq07, beliefPenaltyMPC_lubbysub07, beliefPenaltyMPC_llbbyslb07);
beliefPenaltyMPC_LA_INEQ_B_GRAD_17_17_7(beliefPenaltyMPC_lub08, beliefPenaltyMPC_sub08, beliefPenaltyMPC_riub08, beliefPenaltyMPC_llb08, beliefPenaltyMPC_slb08, beliefPenaltyMPC_rilb08, beliefPenaltyMPC_lbIdx08, beliefPenaltyMPC_ubIdx08, beliefPenaltyMPC_grad_ineq08, beliefPenaltyMPC_lubbysub08, beliefPenaltyMPC_llbbyslb08);
beliefPenaltyMPC_LA_INEQ_B_GRAD_17_17_7(beliefPenaltyMPC_lub09, beliefPenaltyMPC_sub09, beliefPenaltyMPC_riub09, beliefPenaltyMPC_llb09, beliefPenaltyMPC_slb09, beliefPenaltyMPC_rilb09, beliefPenaltyMPC_lbIdx09, beliefPenaltyMPC_ubIdx09, beliefPenaltyMPC_grad_ineq09, beliefPenaltyMPC_lubbysub09, beliefPenaltyMPC_llbbyslb09);
beliefPenaltyMPC_LA_INEQ_B_GRAD_17_17_7(beliefPenaltyMPC_lub10, beliefPenaltyMPC_sub10, beliefPenaltyMPC_riub10, beliefPenaltyMPC_llb10, beliefPenaltyMPC_slb10, beliefPenaltyMPC_rilb10, beliefPenaltyMPC_lbIdx10, beliefPenaltyMPC_ubIdx10, beliefPenaltyMPC_grad_ineq10, beliefPenaltyMPC_lubbysub10, beliefPenaltyMPC_llbbyslb10);
beliefPenaltyMPC_LA_INEQ_B_GRAD_17_17_7(beliefPenaltyMPC_lub11, beliefPenaltyMPC_sub11, beliefPenaltyMPC_riub11, beliefPenaltyMPC_llb11, beliefPenaltyMPC_slb11, beliefPenaltyMPC_rilb11, beliefPenaltyMPC_lbIdx11, beliefPenaltyMPC_ubIdx11, beliefPenaltyMPC_grad_ineq11, beliefPenaltyMPC_lubbysub11, beliefPenaltyMPC_llbbyslb11);
beliefPenaltyMPC_LA_INEQ_B_GRAD_17_17_7(beliefPenaltyMPC_lub12, beliefPenaltyMPC_sub12, beliefPenaltyMPC_riub12, beliefPenaltyMPC_llb12, beliefPenaltyMPC_slb12, beliefPenaltyMPC_rilb12, beliefPenaltyMPC_lbIdx12, beliefPenaltyMPC_ubIdx12, beliefPenaltyMPC_grad_ineq12, beliefPenaltyMPC_lubbysub12, beliefPenaltyMPC_llbbyslb12);
beliefPenaltyMPC_LA_INEQ_B_GRAD_17_17_7(beliefPenaltyMPC_lub13, beliefPenaltyMPC_sub13, beliefPenaltyMPC_riub13, beliefPenaltyMPC_llb13, beliefPenaltyMPC_slb13, beliefPenaltyMPC_rilb13, beliefPenaltyMPC_lbIdx13, beliefPenaltyMPC_ubIdx13, beliefPenaltyMPC_grad_ineq13, beliefPenaltyMPC_lubbysub13, beliefPenaltyMPC_llbbyslb13);
beliefPenaltyMPC_LA_INEQ_B_GRAD_17_17_7(beliefPenaltyMPC_lub14, beliefPenaltyMPC_sub14, beliefPenaltyMPC_riub14, beliefPenaltyMPC_llb14, beliefPenaltyMPC_slb14, beliefPenaltyMPC_rilb14, beliefPenaltyMPC_lbIdx14, beliefPenaltyMPC_ubIdx14, beliefPenaltyMPC_grad_ineq14, beliefPenaltyMPC_lubbysub14, beliefPenaltyMPC_llbbyslb14);
beliefPenaltyMPC_LA_INEQ_B_GRAD_17_17_7(beliefPenaltyMPC_lub15, beliefPenaltyMPC_sub15, beliefPenaltyMPC_riub15, beliefPenaltyMPC_llb15, beliefPenaltyMPC_slb15, beliefPenaltyMPC_rilb15, beliefPenaltyMPC_lbIdx15, beliefPenaltyMPC_ubIdx15, beliefPenaltyMPC_grad_ineq15, beliefPenaltyMPC_lubbysub15, beliefPenaltyMPC_llbbyslb15);
beliefPenaltyMPC_LA_INEQ_B_GRAD_17_17_7(beliefPenaltyMPC_lub16, beliefPenaltyMPC_sub16, beliefPenaltyMPC_riub16, beliefPenaltyMPC_llb16, beliefPenaltyMPC_slb16, beliefPenaltyMPC_rilb16, beliefPenaltyMPC_lbIdx16, beliefPenaltyMPC_ubIdx16, beliefPenaltyMPC_grad_ineq16, beliefPenaltyMPC_lubbysub16, beliefPenaltyMPC_llbbyslb16);
beliefPenaltyMPC_LA_INEQ_B_GRAD_17_17_7(beliefPenaltyMPC_lub17, beliefPenaltyMPC_sub17, beliefPenaltyMPC_riub17, beliefPenaltyMPC_llb17, beliefPenaltyMPC_slb17, beliefPenaltyMPC_rilb17, beliefPenaltyMPC_lbIdx17, beliefPenaltyMPC_ubIdx17, beliefPenaltyMPC_grad_ineq17, beliefPenaltyMPC_lubbysub17, beliefPenaltyMPC_llbbyslb17);
beliefPenaltyMPC_LA_INEQ_B_GRAD_17_17_7(beliefPenaltyMPC_lub18, beliefPenaltyMPC_sub18, beliefPenaltyMPC_riub18, beliefPenaltyMPC_llb18, beliefPenaltyMPC_slb18, beliefPenaltyMPC_rilb18, beliefPenaltyMPC_lbIdx18, beliefPenaltyMPC_ubIdx18, beliefPenaltyMPC_grad_ineq18, beliefPenaltyMPC_lubbysub18, beliefPenaltyMPC_llbbyslb18);
beliefPenaltyMPC_LA_INEQ_B_GRAD_17_17_7(beliefPenaltyMPC_lub19, beliefPenaltyMPC_sub19, beliefPenaltyMPC_riub19, beliefPenaltyMPC_llb19, beliefPenaltyMPC_slb19, beliefPenaltyMPC_rilb19, beliefPenaltyMPC_lbIdx19, beliefPenaltyMPC_ubIdx19, beliefPenaltyMPC_grad_ineq19, beliefPenaltyMPC_lubbysub19, beliefPenaltyMPC_llbbyslb19);
beliefPenaltyMPC_LA_INEQ_B_GRAD_17_17_7(beliefPenaltyMPC_lub20, beliefPenaltyMPC_sub20, beliefPenaltyMPC_riub20, beliefPenaltyMPC_llb20, beliefPenaltyMPC_slb20, beliefPenaltyMPC_rilb20, beliefPenaltyMPC_lbIdx20, beliefPenaltyMPC_ubIdx20, beliefPenaltyMPC_grad_ineq20, beliefPenaltyMPC_lubbysub20, beliefPenaltyMPC_llbbyslb20);
beliefPenaltyMPC_LA_INEQ_B_GRAD_17_17_7(beliefPenaltyMPC_lub21, beliefPenaltyMPC_sub21, beliefPenaltyMPC_riub21, beliefPenaltyMPC_llb21, beliefPenaltyMPC_slb21, beliefPenaltyMPC_rilb21, beliefPenaltyMPC_lbIdx21, beliefPenaltyMPC_ubIdx21, beliefPenaltyMPC_grad_ineq21, beliefPenaltyMPC_lubbysub21, beliefPenaltyMPC_llbbyslb21);
beliefPenaltyMPC_LA_INEQ_B_GRAD_17_17_7(beliefPenaltyMPC_lub22, beliefPenaltyMPC_sub22, beliefPenaltyMPC_riub22, beliefPenaltyMPC_llb22, beliefPenaltyMPC_slb22, beliefPenaltyMPC_rilb22, beliefPenaltyMPC_lbIdx22, beliefPenaltyMPC_ubIdx22, beliefPenaltyMPC_grad_ineq22, beliefPenaltyMPC_lubbysub22, beliefPenaltyMPC_llbbyslb22);
beliefPenaltyMPC_LA_INEQ_B_GRAD_17_17_7(beliefPenaltyMPC_lub23, beliefPenaltyMPC_sub23, beliefPenaltyMPC_riub23, beliefPenaltyMPC_llb23, beliefPenaltyMPC_slb23, beliefPenaltyMPC_rilb23, beliefPenaltyMPC_lbIdx23, beliefPenaltyMPC_ubIdx23, beliefPenaltyMPC_grad_ineq23, beliefPenaltyMPC_lubbysub23, beliefPenaltyMPC_llbbyslb23);
beliefPenaltyMPC_LA_INEQ_B_GRAD_5_5_5(beliefPenaltyMPC_lub24, beliefPenaltyMPC_sub24, beliefPenaltyMPC_riub24, beliefPenaltyMPC_llb24, beliefPenaltyMPC_slb24, beliefPenaltyMPC_rilb24, beliefPenaltyMPC_lbIdx24, beliefPenaltyMPC_ubIdx24, beliefPenaltyMPC_grad_ineq24, beliefPenaltyMPC_lubbysub24, beliefPenaltyMPC_llbbyslb24);
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
beliefPenaltyMPC_LA_VVADD3_413(beliefPenaltyMPC_grad_cost, beliefPenaltyMPC_grad_eq, beliefPenaltyMPC_grad_ineq, beliefPenaltyMPC_rd);
beliefPenaltyMPC_LA_DIAG_CHOL_LBUB_17_17_7(beliefPenaltyMPC_H00, beliefPenaltyMPC_llbbyslb00, beliefPenaltyMPC_lbIdx00, beliefPenaltyMPC_lubbysub00, beliefPenaltyMPC_ubIdx00, beliefPenaltyMPC_Phi00);
beliefPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_10_17(beliefPenaltyMPC_Phi00, params->C1, beliefPenaltyMPC_V00);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_17(beliefPenaltyMPC_Phi00, beliefPenaltyMPC_rd00, beliefPenaltyMPC_Lbyrd00);
beliefPenaltyMPC_LA_DIAG_CHOL_LBUB_17_17_7(beliefPenaltyMPC_H00, beliefPenaltyMPC_llbbyslb01, beliefPenaltyMPC_lbIdx01, beliefPenaltyMPC_lubbysub01, beliefPenaltyMPC_ubIdx01, beliefPenaltyMPC_Phi01);
beliefPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_5_17(beliefPenaltyMPC_Phi01, params->C2, beliefPenaltyMPC_V01);
beliefPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_10_17(beliefPenaltyMPC_Phi01, beliefPenaltyMPC_D01, beliefPenaltyMPC_W01);
beliefPenaltyMPC_LA_DENSE_MMTM_10_17_5(beliefPenaltyMPC_W01, beliefPenaltyMPC_V01, beliefPenaltyMPC_Ysd01);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_17(beliefPenaltyMPC_Phi01, beliefPenaltyMPC_rd01, beliefPenaltyMPC_Lbyrd01);
beliefPenaltyMPC_LA_DIAG_CHOL_LBUB_17_17_7(beliefPenaltyMPC_H00, beliefPenaltyMPC_llbbyslb02, beliefPenaltyMPC_lbIdx02, beliefPenaltyMPC_lubbysub02, beliefPenaltyMPC_ubIdx02, beliefPenaltyMPC_Phi02);
beliefPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_5_17(beliefPenaltyMPC_Phi02, params->C3, beliefPenaltyMPC_V02);
beliefPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_5_17(beliefPenaltyMPC_Phi02, beliefPenaltyMPC_D02, beliefPenaltyMPC_W02);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMTM_5_17_5(beliefPenaltyMPC_W02, beliefPenaltyMPC_V02, beliefPenaltyMPC_Ysd02);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_17(beliefPenaltyMPC_Phi02, beliefPenaltyMPC_rd02, beliefPenaltyMPC_Lbyrd02);
beliefPenaltyMPC_LA_DIAG_CHOL_LBUB_17_17_7(beliefPenaltyMPC_H00, beliefPenaltyMPC_llbbyslb03, beliefPenaltyMPC_lbIdx03, beliefPenaltyMPC_lubbysub03, beliefPenaltyMPC_ubIdx03, beliefPenaltyMPC_Phi03);
beliefPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_5_17(beliefPenaltyMPC_Phi03, params->C4, beliefPenaltyMPC_V03);
beliefPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_5_17(beliefPenaltyMPC_Phi03, beliefPenaltyMPC_D02, beliefPenaltyMPC_W03);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMTM_5_17_5(beliefPenaltyMPC_W03, beliefPenaltyMPC_V03, beliefPenaltyMPC_Ysd03);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_17(beliefPenaltyMPC_Phi03, beliefPenaltyMPC_rd03, beliefPenaltyMPC_Lbyrd03);
beliefPenaltyMPC_LA_DIAG_CHOL_LBUB_17_17_7(beliefPenaltyMPC_H00, beliefPenaltyMPC_llbbyslb04, beliefPenaltyMPC_lbIdx04, beliefPenaltyMPC_lubbysub04, beliefPenaltyMPC_ubIdx04, beliefPenaltyMPC_Phi04);
beliefPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_5_17(beliefPenaltyMPC_Phi04, params->C5, beliefPenaltyMPC_V04);
beliefPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_5_17(beliefPenaltyMPC_Phi04, beliefPenaltyMPC_D02, beliefPenaltyMPC_W04);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMTM_5_17_5(beliefPenaltyMPC_W04, beliefPenaltyMPC_V04, beliefPenaltyMPC_Ysd04);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_17(beliefPenaltyMPC_Phi04, beliefPenaltyMPC_rd04, beliefPenaltyMPC_Lbyrd04);
beliefPenaltyMPC_LA_DIAG_CHOL_LBUB_17_17_7(beliefPenaltyMPC_H00, beliefPenaltyMPC_llbbyslb05, beliefPenaltyMPC_lbIdx05, beliefPenaltyMPC_lubbysub05, beliefPenaltyMPC_ubIdx05, beliefPenaltyMPC_Phi05);
beliefPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_5_17(beliefPenaltyMPC_Phi05, params->C6, beliefPenaltyMPC_V05);
beliefPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_5_17(beliefPenaltyMPC_Phi05, beliefPenaltyMPC_D02, beliefPenaltyMPC_W05);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMTM_5_17_5(beliefPenaltyMPC_W05, beliefPenaltyMPC_V05, beliefPenaltyMPC_Ysd05);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_17(beliefPenaltyMPC_Phi05, beliefPenaltyMPC_rd05, beliefPenaltyMPC_Lbyrd05);
beliefPenaltyMPC_LA_DIAG_CHOL_LBUB_17_17_7(beliefPenaltyMPC_H00, beliefPenaltyMPC_llbbyslb06, beliefPenaltyMPC_lbIdx06, beliefPenaltyMPC_lubbysub06, beliefPenaltyMPC_ubIdx06, beliefPenaltyMPC_Phi06);
beliefPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_5_17(beliefPenaltyMPC_Phi06, params->C7, beliefPenaltyMPC_V06);
beliefPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_5_17(beliefPenaltyMPC_Phi06, beliefPenaltyMPC_D02, beliefPenaltyMPC_W06);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMTM_5_17_5(beliefPenaltyMPC_W06, beliefPenaltyMPC_V06, beliefPenaltyMPC_Ysd06);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_17(beliefPenaltyMPC_Phi06, beliefPenaltyMPC_rd06, beliefPenaltyMPC_Lbyrd06);
beliefPenaltyMPC_LA_DIAG_CHOL_LBUB_17_17_7(beliefPenaltyMPC_H00, beliefPenaltyMPC_llbbyslb07, beliefPenaltyMPC_lbIdx07, beliefPenaltyMPC_lubbysub07, beliefPenaltyMPC_ubIdx07, beliefPenaltyMPC_Phi07);
beliefPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_5_17(beliefPenaltyMPC_Phi07, params->C8, beliefPenaltyMPC_V07);
beliefPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_5_17(beliefPenaltyMPC_Phi07, beliefPenaltyMPC_D02, beliefPenaltyMPC_W07);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMTM_5_17_5(beliefPenaltyMPC_W07, beliefPenaltyMPC_V07, beliefPenaltyMPC_Ysd07);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_17(beliefPenaltyMPC_Phi07, beliefPenaltyMPC_rd07, beliefPenaltyMPC_Lbyrd07);
beliefPenaltyMPC_LA_DIAG_CHOL_LBUB_17_17_7(beliefPenaltyMPC_H00, beliefPenaltyMPC_llbbyslb08, beliefPenaltyMPC_lbIdx08, beliefPenaltyMPC_lubbysub08, beliefPenaltyMPC_ubIdx08, beliefPenaltyMPC_Phi08);
beliefPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_5_17(beliefPenaltyMPC_Phi08, params->C9, beliefPenaltyMPC_V08);
beliefPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_5_17(beliefPenaltyMPC_Phi08, beliefPenaltyMPC_D02, beliefPenaltyMPC_W08);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMTM_5_17_5(beliefPenaltyMPC_W08, beliefPenaltyMPC_V08, beliefPenaltyMPC_Ysd08);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_17(beliefPenaltyMPC_Phi08, beliefPenaltyMPC_rd08, beliefPenaltyMPC_Lbyrd08);
beliefPenaltyMPC_LA_DIAG_CHOL_LBUB_17_17_7(beliefPenaltyMPC_H00, beliefPenaltyMPC_llbbyslb09, beliefPenaltyMPC_lbIdx09, beliefPenaltyMPC_lubbysub09, beliefPenaltyMPC_ubIdx09, beliefPenaltyMPC_Phi09);
beliefPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_5_17(beliefPenaltyMPC_Phi09, params->C10, beliefPenaltyMPC_V09);
beliefPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_5_17(beliefPenaltyMPC_Phi09, beliefPenaltyMPC_D02, beliefPenaltyMPC_W09);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMTM_5_17_5(beliefPenaltyMPC_W09, beliefPenaltyMPC_V09, beliefPenaltyMPC_Ysd09);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_17(beliefPenaltyMPC_Phi09, beliefPenaltyMPC_rd09, beliefPenaltyMPC_Lbyrd09);
beliefPenaltyMPC_LA_DIAG_CHOL_LBUB_17_17_7(beliefPenaltyMPC_H00, beliefPenaltyMPC_llbbyslb10, beliefPenaltyMPC_lbIdx10, beliefPenaltyMPC_lubbysub10, beliefPenaltyMPC_ubIdx10, beliefPenaltyMPC_Phi10);
beliefPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_5_17(beliefPenaltyMPC_Phi10, params->C11, beliefPenaltyMPC_V10);
beliefPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_5_17(beliefPenaltyMPC_Phi10, beliefPenaltyMPC_D02, beliefPenaltyMPC_W10);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMTM_5_17_5(beliefPenaltyMPC_W10, beliefPenaltyMPC_V10, beliefPenaltyMPC_Ysd10);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_17(beliefPenaltyMPC_Phi10, beliefPenaltyMPC_rd10, beliefPenaltyMPC_Lbyrd10);
beliefPenaltyMPC_LA_DIAG_CHOL_LBUB_17_17_7(beliefPenaltyMPC_H00, beliefPenaltyMPC_llbbyslb11, beliefPenaltyMPC_lbIdx11, beliefPenaltyMPC_lubbysub11, beliefPenaltyMPC_ubIdx11, beliefPenaltyMPC_Phi11);
beliefPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_5_17(beliefPenaltyMPC_Phi11, params->C12, beliefPenaltyMPC_V11);
beliefPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_5_17(beliefPenaltyMPC_Phi11, beliefPenaltyMPC_D02, beliefPenaltyMPC_W11);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMTM_5_17_5(beliefPenaltyMPC_W11, beliefPenaltyMPC_V11, beliefPenaltyMPC_Ysd11);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_17(beliefPenaltyMPC_Phi11, beliefPenaltyMPC_rd11, beliefPenaltyMPC_Lbyrd11);
beliefPenaltyMPC_LA_DIAG_CHOL_LBUB_17_17_7(beliefPenaltyMPC_H00, beliefPenaltyMPC_llbbyslb12, beliefPenaltyMPC_lbIdx12, beliefPenaltyMPC_lubbysub12, beliefPenaltyMPC_ubIdx12, beliefPenaltyMPC_Phi12);
beliefPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_5_17(beliefPenaltyMPC_Phi12, params->C13, beliefPenaltyMPC_V12);
beliefPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_5_17(beliefPenaltyMPC_Phi12, beliefPenaltyMPC_D02, beliefPenaltyMPC_W12);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMTM_5_17_5(beliefPenaltyMPC_W12, beliefPenaltyMPC_V12, beliefPenaltyMPC_Ysd12);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_17(beliefPenaltyMPC_Phi12, beliefPenaltyMPC_rd12, beliefPenaltyMPC_Lbyrd12);
beliefPenaltyMPC_LA_DIAG_CHOL_LBUB_17_17_7(beliefPenaltyMPC_H00, beliefPenaltyMPC_llbbyslb13, beliefPenaltyMPC_lbIdx13, beliefPenaltyMPC_lubbysub13, beliefPenaltyMPC_ubIdx13, beliefPenaltyMPC_Phi13);
beliefPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_5_17(beliefPenaltyMPC_Phi13, params->C14, beliefPenaltyMPC_V13);
beliefPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_5_17(beliefPenaltyMPC_Phi13, beliefPenaltyMPC_D02, beliefPenaltyMPC_W13);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMTM_5_17_5(beliefPenaltyMPC_W13, beliefPenaltyMPC_V13, beliefPenaltyMPC_Ysd13);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_17(beliefPenaltyMPC_Phi13, beliefPenaltyMPC_rd13, beliefPenaltyMPC_Lbyrd13);
beliefPenaltyMPC_LA_DIAG_CHOL_LBUB_17_17_7(beliefPenaltyMPC_H00, beliefPenaltyMPC_llbbyslb14, beliefPenaltyMPC_lbIdx14, beliefPenaltyMPC_lubbysub14, beliefPenaltyMPC_ubIdx14, beliefPenaltyMPC_Phi14);
beliefPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_5_17(beliefPenaltyMPC_Phi14, params->C15, beliefPenaltyMPC_V14);
beliefPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_5_17(beliefPenaltyMPC_Phi14, beliefPenaltyMPC_D02, beliefPenaltyMPC_W14);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMTM_5_17_5(beliefPenaltyMPC_W14, beliefPenaltyMPC_V14, beliefPenaltyMPC_Ysd14);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_17(beliefPenaltyMPC_Phi14, beliefPenaltyMPC_rd14, beliefPenaltyMPC_Lbyrd14);
beliefPenaltyMPC_LA_DIAG_CHOL_LBUB_17_17_7(beliefPenaltyMPC_H00, beliefPenaltyMPC_llbbyslb15, beliefPenaltyMPC_lbIdx15, beliefPenaltyMPC_lubbysub15, beliefPenaltyMPC_ubIdx15, beliefPenaltyMPC_Phi15);
beliefPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_5_17(beliefPenaltyMPC_Phi15, params->C16, beliefPenaltyMPC_V15);
beliefPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_5_17(beliefPenaltyMPC_Phi15, beliefPenaltyMPC_D02, beliefPenaltyMPC_W15);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMTM_5_17_5(beliefPenaltyMPC_W15, beliefPenaltyMPC_V15, beliefPenaltyMPC_Ysd15);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_17(beliefPenaltyMPC_Phi15, beliefPenaltyMPC_rd15, beliefPenaltyMPC_Lbyrd15);
beliefPenaltyMPC_LA_DIAG_CHOL_LBUB_17_17_7(beliefPenaltyMPC_H00, beliefPenaltyMPC_llbbyslb16, beliefPenaltyMPC_lbIdx16, beliefPenaltyMPC_lubbysub16, beliefPenaltyMPC_ubIdx16, beliefPenaltyMPC_Phi16);
beliefPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_5_17(beliefPenaltyMPC_Phi16, params->C17, beliefPenaltyMPC_V16);
beliefPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_5_17(beliefPenaltyMPC_Phi16, beliefPenaltyMPC_D02, beliefPenaltyMPC_W16);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMTM_5_17_5(beliefPenaltyMPC_W16, beliefPenaltyMPC_V16, beliefPenaltyMPC_Ysd16);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_17(beliefPenaltyMPC_Phi16, beliefPenaltyMPC_rd16, beliefPenaltyMPC_Lbyrd16);
beliefPenaltyMPC_LA_DIAG_CHOL_LBUB_17_17_7(beliefPenaltyMPC_H00, beliefPenaltyMPC_llbbyslb17, beliefPenaltyMPC_lbIdx17, beliefPenaltyMPC_lubbysub17, beliefPenaltyMPC_ubIdx17, beliefPenaltyMPC_Phi17);
beliefPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_5_17(beliefPenaltyMPC_Phi17, params->C18, beliefPenaltyMPC_V17);
beliefPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_5_17(beliefPenaltyMPC_Phi17, beliefPenaltyMPC_D02, beliefPenaltyMPC_W17);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMTM_5_17_5(beliefPenaltyMPC_W17, beliefPenaltyMPC_V17, beliefPenaltyMPC_Ysd17);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_17(beliefPenaltyMPC_Phi17, beliefPenaltyMPC_rd17, beliefPenaltyMPC_Lbyrd17);
beliefPenaltyMPC_LA_DIAG_CHOL_LBUB_17_17_7(beliefPenaltyMPC_H00, beliefPenaltyMPC_llbbyslb18, beliefPenaltyMPC_lbIdx18, beliefPenaltyMPC_lubbysub18, beliefPenaltyMPC_ubIdx18, beliefPenaltyMPC_Phi18);
beliefPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_5_17(beliefPenaltyMPC_Phi18, params->C19, beliefPenaltyMPC_V18);
beliefPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_5_17(beliefPenaltyMPC_Phi18, beliefPenaltyMPC_D02, beliefPenaltyMPC_W18);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMTM_5_17_5(beliefPenaltyMPC_W18, beliefPenaltyMPC_V18, beliefPenaltyMPC_Ysd18);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_17(beliefPenaltyMPC_Phi18, beliefPenaltyMPC_rd18, beliefPenaltyMPC_Lbyrd18);
beliefPenaltyMPC_LA_DIAG_CHOL_LBUB_17_17_7(beliefPenaltyMPC_H00, beliefPenaltyMPC_llbbyslb19, beliefPenaltyMPC_lbIdx19, beliefPenaltyMPC_lubbysub19, beliefPenaltyMPC_ubIdx19, beliefPenaltyMPC_Phi19);
beliefPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_5_17(beliefPenaltyMPC_Phi19, params->C20, beliefPenaltyMPC_V19);
beliefPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_5_17(beliefPenaltyMPC_Phi19, beliefPenaltyMPC_D02, beliefPenaltyMPC_W19);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMTM_5_17_5(beliefPenaltyMPC_W19, beliefPenaltyMPC_V19, beliefPenaltyMPC_Ysd19);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_17(beliefPenaltyMPC_Phi19, beliefPenaltyMPC_rd19, beliefPenaltyMPC_Lbyrd19);
beliefPenaltyMPC_LA_DIAG_CHOL_LBUB_17_17_7(beliefPenaltyMPC_H00, beliefPenaltyMPC_llbbyslb20, beliefPenaltyMPC_lbIdx20, beliefPenaltyMPC_lubbysub20, beliefPenaltyMPC_ubIdx20, beliefPenaltyMPC_Phi20);
beliefPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_5_17(beliefPenaltyMPC_Phi20, params->C21, beliefPenaltyMPC_V20);
beliefPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_5_17(beliefPenaltyMPC_Phi20, beliefPenaltyMPC_D02, beliefPenaltyMPC_W20);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMTM_5_17_5(beliefPenaltyMPC_W20, beliefPenaltyMPC_V20, beliefPenaltyMPC_Ysd20);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_17(beliefPenaltyMPC_Phi20, beliefPenaltyMPC_rd20, beliefPenaltyMPC_Lbyrd20);
beliefPenaltyMPC_LA_DIAG_CHOL_LBUB_17_17_7(beliefPenaltyMPC_H00, beliefPenaltyMPC_llbbyslb21, beliefPenaltyMPC_lbIdx21, beliefPenaltyMPC_lubbysub21, beliefPenaltyMPC_ubIdx21, beliefPenaltyMPC_Phi21);
beliefPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_5_17(beliefPenaltyMPC_Phi21, params->C22, beliefPenaltyMPC_V21);
beliefPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_5_17(beliefPenaltyMPC_Phi21, beliefPenaltyMPC_D02, beliefPenaltyMPC_W21);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMTM_5_17_5(beliefPenaltyMPC_W21, beliefPenaltyMPC_V21, beliefPenaltyMPC_Ysd21);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_17(beliefPenaltyMPC_Phi21, beliefPenaltyMPC_rd21, beliefPenaltyMPC_Lbyrd21);
beliefPenaltyMPC_LA_DIAG_CHOL_LBUB_17_17_7(beliefPenaltyMPC_H00, beliefPenaltyMPC_llbbyslb22, beliefPenaltyMPC_lbIdx22, beliefPenaltyMPC_lubbysub22, beliefPenaltyMPC_ubIdx22, beliefPenaltyMPC_Phi22);
beliefPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_5_17(beliefPenaltyMPC_Phi22, params->C23, beliefPenaltyMPC_V22);
beliefPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_5_17(beliefPenaltyMPC_Phi22, beliefPenaltyMPC_D02, beliefPenaltyMPC_W22);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMTM_5_17_5(beliefPenaltyMPC_W22, beliefPenaltyMPC_V22, beliefPenaltyMPC_Ysd22);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_17(beliefPenaltyMPC_Phi22, beliefPenaltyMPC_rd22, beliefPenaltyMPC_Lbyrd22);
beliefPenaltyMPC_LA_DIAG_CHOL_LBUB_17_17_7(beliefPenaltyMPC_H00, beliefPenaltyMPC_llbbyslb23, beliefPenaltyMPC_lbIdx23, beliefPenaltyMPC_lubbysub23, beliefPenaltyMPC_ubIdx23, beliefPenaltyMPC_Phi23);
beliefPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_5_17(beliefPenaltyMPC_Phi23, params->C24, beliefPenaltyMPC_V23);
beliefPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_5_17(beliefPenaltyMPC_Phi23, beliefPenaltyMPC_D02, beliefPenaltyMPC_W23);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMTM_5_17_5(beliefPenaltyMPC_W23, beliefPenaltyMPC_V23, beliefPenaltyMPC_Ysd23);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_17(beliefPenaltyMPC_Phi23, beliefPenaltyMPC_rd23, beliefPenaltyMPC_Lbyrd23);
beliefPenaltyMPC_LA_DIAG_CHOL_ONELOOP_LBUB_5_5_5(beliefPenaltyMPC_H24, beliefPenaltyMPC_llbbyslb24, beliefPenaltyMPC_lbIdx24, beliefPenaltyMPC_lubbysub24, beliefPenaltyMPC_ubIdx24, beliefPenaltyMPC_Phi24);
beliefPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_5_5(beliefPenaltyMPC_Phi24, beliefPenaltyMPC_D24, beliefPenaltyMPC_W24);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_5(beliefPenaltyMPC_Phi24, beliefPenaltyMPC_rd24, beliefPenaltyMPC_Lbyrd24);
beliefPenaltyMPC_LA_DENSE_MMT2_10_17_17(beliefPenaltyMPC_V00, beliefPenaltyMPC_W01, beliefPenaltyMPC_Yd00);
beliefPenaltyMPC_LA_DENSE_MVMSUB2_10_17_17(beliefPenaltyMPC_V00, beliefPenaltyMPC_Lbyrd00, beliefPenaltyMPC_W01, beliefPenaltyMPC_Lbyrd01, beliefPenaltyMPC_re00, beliefPenaltyMPC_beta00);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMT2_5_17_17(beliefPenaltyMPC_V01, beliefPenaltyMPC_W02, beliefPenaltyMPC_Yd01);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_5_17_17(beliefPenaltyMPC_V01, beliefPenaltyMPC_Lbyrd01, beliefPenaltyMPC_W02, beliefPenaltyMPC_Lbyrd02, beliefPenaltyMPC_re01, beliefPenaltyMPC_beta01);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMT2_5_17_17(beliefPenaltyMPC_V02, beliefPenaltyMPC_W03, beliefPenaltyMPC_Yd02);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_5_17_17(beliefPenaltyMPC_V02, beliefPenaltyMPC_Lbyrd02, beliefPenaltyMPC_W03, beliefPenaltyMPC_Lbyrd03, beliefPenaltyMPC_re02, beliefPenaltyMPC_beta02);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMT2_5_17_17(beliefPenaltyMPC_V03, beliefPenaltyMPC_W04, beliefPenaltyMPC_Yd03);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_5_17_17(beliefPenaltyMPC_V03, beliefPenaltyMPC_Lbyrd03, beliefPenaltyMPC_W04, beliefPenaltyMPC_Lbyrd04, beliefPenaltyMPC_re03, beliefPenaltyMPC_beta03);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMT2_5_17_17(beliefPenaltyMPC_V04, beliefPenaltyMPC_W05, beliefPenaltyMPC_Yd04);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_5_17_17(beliefPenaltyMPC_V04, beliefPenaltyMPC_Lbyrd04, beliefPenaltyMPC_W05, beliefPenaltyMPC_Lbyrd05, beliefPenaltyMPC_re04, beliefPenaltyMPC_beta04);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMT2_5_17_17(beliefPenaltyMPC_V05, beliefPenaltyMPC_W06, beliefPenaltyMPC_Yd05);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_5_17_17(beliefPenaltyMPC_V05, beliefPenaltyMPC_Lbyrd05, beliefPenaltyMPC_W06, beliefPenaltyMPC_Lbyrd06, beliefPenaltyMPC_re05, beliefPenaltyMPC_beta05);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMT2_5_17_17(beliefPenaltyMPC_V06, beliefPenaltyMPC_W07, beliefPenaltyMPC_Yd06);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_5_17_17(beliefPenaltyMPC_V06, beliefPenaltyMPC_Lbyrd06, beliefPenaltyMPC_W07, beliefPenaltyMPC_Lbyrd07, beliefPenaltyMPC_re06, beliefPenaltyMPC_beta06);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMT2_5_17_17(beliefPenaltyMPC_V07, beliefPenaltyMPC_W08, beliefPenaltyMPC_Yd07);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_5_17_17(beliefPenaltyMPC_V07, beliefPenaltyMPC_Lbyrd07, beliefPenaltyMPC_W08, beliefPenaltyMPC_Lbyrd08, beliefPenaltyMPC_re07, beliefPenaltyMPC_beta07);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMT2_5_17_17(beliefPenaltyMPC_V08, beliefPenaltyMPC_W09, beliefPenaltyMPC_Yd08);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_5_17_17(beliefPenaltyMPC_V08, beliefPenaltyMPC_Lbyrd08, beliefPenaltyMPC_W09, beliefPenaltyMPC_Lbyrd09, beliefPenaltyMPC_re08, beliefPenaltyMPC_beta08);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMT2_5_17_17(beliefPenaltyMPC_V09, beliefPenaltyMPC_W10, beliefPenaltyMPC_Yd09);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_5_17_17(beliefPenaltyMPC_V09, beliefPenaltyMPC_Lbyrd09, beliefPenaltyMPC_W10, beliefPenaltyMPC_Lbyrd10, beliefPenaltyMPC_re09, beliefPenaltyMPC_beta09);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMT2_5_17_17(beliefPenaltyMPC_V10, beliefPenaltyMPC_W11, beliefPenaltyMPC_Yd10);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_5_17_17(beliefPenaltyMPC_V10, beliefPenaltyMPC_Lbyrd10, beliefPenaltyMPC_W11, beliefPenaltyMPC_Lbyrd11, beliefPenaltyMPC_re10, beliefPenaltyMPC_beta10);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMT2_5_17_17(beliefPenaltyMPC_V11, beliefPenaltyMPC_W12, beliefPenaltyMPC_Yd11);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_5_17_17(beliefPenaltyMPC_V11, beliefPenaltyMPC_Lbyrd11, beliefPenaltyMPC_W12, beliefPenaltyMPC_Lbyrd12, beliefPenaltyMPC_re11, beliefPenaltyMPC_beta11);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMT2_5_17_17(beliefPenaltyMPC_V12, beliefPenaltyMPC_W13, beliefPenaltyMPC_Yd12);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_5_17_17(beliefPenaltyMPC_V12, beliefPenaltyMPC_Lbyrd12, beliefPenaltyMPC_W13, beliefPenaltyMPC_Lbyrd13, beliefPenaltyMPC_re12, beliefPenaltyMPC_beta12);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMT2_5_17_17(beliefPenaltyMPC_V13, beliefPenaltyMPC_W14, beliefPenaltyMPC_Yd13);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_5_17_17(beliefPenaltyMPC_V13, beliefPenaltyMPC_Lbyrd13, beliefPenaltyMPC_W14, beliefPenaltyMPC_Lbyrd14, beliefPenaltyMPC_re13, beliefPenaltyMPC_beta13);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMT2_5_17_17(beliefPenaltyMPC_V14, beliefPenaltyMPC_W15, beliefPenaltyMPC_Yd14);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_5_17_17(beliefPenaltyMPC_V14, beliefPenaltyMPC_Lbyrd14, beliefPenaltyMPC_W15, beliefPenaltyMPC_Lbyrd15, beliefPenaltyMPC_re14, beliefPenaltyMPC_beta14);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMT2_5_17_17(beliefPenaltyMPC_V15, beliefPenaltyMPC_W16, beliefPenaltyMPC_Yd15);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_5_17_17(beliefPenaltyMPC_V15, beliefPenaltyMPC_Lbyrd15, beliefPenaltyMPC_W16, beliefPenaltyMPC_Lbyrd16, beliefPenaltyMPC_re15, beliefPenaltyMPC_beta15);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMT2_5_17_17(beliefPenaltyMPC_V16, beliefPenaltyMPC_W17, beliefPenaltyMPC_Yd16);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_5_17_17(beliefPenaltyMPC_V16, beliefPenaltyMPC_Lbyrd16, beliefPenaltyMPC_W17, beliefPenaltyMPC_Lbyrd17, beliefPenaltyMPC_re16, beliefPenaltyMPC_beta16);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMT2_5_17_17(beliefPenaltyMPC_V17, beliefPenaltyMPC_W18, beliefPenaltyMPC_Yd17);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_5_17_17(beliefPenaltyMPC_V17, beliefPenaltyMPC_Lbyrd17, beliefPenaltyMPC_W18, beliefPenaltyMPC_Lbyrd18, beliefPenaltyMPC_re17, beliefPenaltyMPC_beta17);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMT2_5_17_17(beliefPenaltyMPC_V18, beliefPenaltyMPC_W19, beliefPenaltyMPC_Yd18);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_5_17_17(beliefPenaltyMPC_V18, beliefPenaltyMPC_Lbyrd18, beliefPenaltyMPC_W19, beliefPenaltyMPC_Lbyrd19, beliefPenaltyMPC_re18, beliefPenaltyMPC_beta18);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMT2_5_17_17(beliefPenaltyMPC_V19, beliefPenaltyMPC_W20, beliefPenaltyMPC_Yd19);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_5_17_17(beliefPenaltyMPC_V19, beliefPenaltyMPC_Lbyrd19, beliefPenaltyMPC_W20, beliefPenaltyMPC_Lbyrd20, beliefPenaltyMPC_re19, beliefPenaltyMPC_beta19);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMT2_5_17_17(beliefPenaltyMPC_V20, beliefPenaltyMPC_W21, beliefPenaltyMPC_Yd20);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_5_17_17(beliefPenaltyMPC_V20, beliefPenaltyMPC_Lbyrd20, beliefPenaltyMPC_W21, beliefPenaltyMPC_Lbyrd21, beliefPenaltyMPC_re20, beliefPenaltyMPC_beta20);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMT2_5_17_17(beliefPenaltyMPC_V21, beliefPenaltyMPC_W22, beliefPenaltyMPC_Yd21);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_5_17_17(beliefPenaltyMPC_V21, beliefPenaltyMPC_Lbyrd21, beliefPenaltyMPC_W22, beliefPenaltyMPC_Lbyrd22, beliefPenaltyMPC_re21, beliefPenaltyMPC_beta21);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMT2_5_17_17(beliefPenaltyMPC_V22, beliefPenaltyMPC_W23, beliefPenaltyMPC_Yd22);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_5_17_17(beliefPenaltyMPC_V22, beliefPenaltyMPC_Lbyrd22, beliefPenaltyMPC_W23, beliefPenaltyMPC_Lbyrd23, beliefPenaltyMPC_re22, beliefPenaltyMPC_beta22);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMT2_5_17_5(beliefPenaltyMPC_V23, beliefPenaltyMPC_W24, beliefPenaltyMPC_Yd23);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_5_17_5(beliefPenaltyMPC_V23, beliefPenaltyMPC_Lbyrd23, beliefPenaltyMPC_W24, beliefPenaltyMPC_Lbyrd24, beliefPenaltyMPC_re23, beliefPenaltyMPC_beta23);
beliefPenaltyMPC_LA_DENSE_CHOL_10(beliefPenaltyMPC_Yd00, beliefPenaltyMPC_Ld00);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_10(beliefPenaltyMPC_Ld00, beliefPenaltyMPC_beta00, beliefPenaltyMPC_yy00);
beliefPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_5_10(beliefPenaltyMPC_Ld00, beliefPenaltyMPC_Ysd01, beliefPenaltyMPC_Lsd01);
beliefPenaltyMPC_LA_DENSE_MMTSUB_5_10(beliefPenaltyMPC_Lsd01, beliefPenaltyMPC_Yd01);
beliefPenaltyMPC_LA_DENSE_CHOL_5(beliefPenaltyMPC_Yd01, beliefPenaltyMPC_Ld01);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_5_10(beliefPenaltyMPC_Lsd01, beliefPenaltyMPC_yy00, beliefPenaltyMPC_beta01, beliefPenaltyMPC_bmy01);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_5(beliefPenaltyMPC_Ld01, beliefPenaltyMPC_bmy01, beliefPenaltyMPC_yy01);
beliefPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_5_5(beliefPenaltyMPC_Ld01, beliefPenaltyMPC_Ysd02, beliefPenaltyMPC_Lsd02);
beliefPenaltyMPC_LA_DENSE_MMTSUB_5_5(beliefPenaltyMPC_Lsd02, beliefPenaltyMPC_Yd02);
beliefPenaltyMPC_LA_DENSE_CHOL_5(beliefPenaltyMPC_Yd02, beliefPenaltyMPC_Ld02);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_5_5(beliefPenaltyMPC_Lsd02, beliefPenaltyMPC_yy01, beliefPenaltyMPC_beta02, beliefPenaltyMPC_bmy02);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_5(beliefPenaltyMPC_Ld02, beliefPenaltyMPC_bmy02, beliefPenaltyMPC_yy02);
beliefPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_5_5(beliefPenaltyMPC_Ld02, beliefPenaltyMPC_Ysd03, beliefPenaltyMPC_Lsd03);
beliefPenaltyMPC_LA_DENSE_MMTSUB_5_5(beliefPenaltyMPC_Lsd03, beliefPenaltyMPC_Yd03);
beliefPenaltyMPC_LA_DENSE_CHOL_5(beliefPenaltyMPC_Yd03, beliefPenaltyMPC_Ld03);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_5_5(beliefPenaltyMPC_Lsd03, beliefPenaltyMPC_yy02, beliefPenaltyMPC_beta03, beliefPenaltyMPC_bmy03);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_5(beliefPenaltyMPC_Ld03, beliefPenaltyMPC_bmy03, beliefPenaltyMPC_yy03);
beliefPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_5_5(beliefPenaltyMPC_Ld03, beliefPenaltyMPC_Ysd04, beliefPenaltyMPC_Lsd04);
beliefPenaltyMPC_LA_DENSE_MMTSUB_5_5(beliefPenaltyMPC_Lsd04, beliefPenaltyMPC_Yd04);
beliefPenaltyMPC_LA_DENSE_CHOL_5(beliefPenaltyMPC_Yd04, beliefPenaltyMPC_Ld04);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_5_5(beliefPenaltyMPC_Lsd04, beliefPenaltyMPC_yy03, beliefPenaltyMPC_beta04, beliefPenaltyMPC_bmy04);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_5(beliefPenaltyMPC_Ld04, beliefPenaltyMPC_bmy04, beliefPenaltyMPC_yy04);
beliefPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_5_5(beliefPenaltyMPC_Ld04, beliefPenaltyMPC_Ysd05, beliefPenaltyMPC_Lsd05);
beliefPenaltyMPC_LA_DENSE_MMTSUB_5_5(beliefPenaltyMPC_Lsd05, beliefPenaltyMPC_Yd05);
beliefPenaltyMPC_LA_DENSE_CHOL_5(beliefPenaltyMPC_Yd05, beliefPenaltyMPC_Ld05);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_5_5(beliefPenaltyMPC_Lsd05, beliefPenaltyMPC_yy04, beliefPenaltyMPC_beta05, beliefPenaltyMPC_bmy05);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_5(beliefPenaltyMPC_Ld05, beliefPenaltyMPC_bmy05, beliefPenaltyMPC_yy05);
beliefPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_5_5(beliefPenaltyMPC_Ld05, beliefPenaltyMPC_Ysd06, beliefPenaltyMPC_Lsd06);
beliefPenaltyMPC_LA_DENSE_MMTSUB_5_5(beliefPenaltyMPC_Lsd06, beliefPenaltyMPC_Yd06);
beliefPenaltyMPC_LA_DENSE_CHOL_5(beliefPenaltyMPC_Yd06, beliefPenaltyMPC_Ld06);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_5_5(beliefPenaltyMPC_Lsd06, beliefPenaltyMPC_yy05, beliefPenaltyMPC_beta06, beliefPenaltyMPC_bmy06);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_5(beliefPenaltyMPC_Ld06, beliefPenaltyMPC_bmy06, beliefPenaltyMPC_yy06);
beliefPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_5_5(beliefPenaltyMPC_Ld06, beliefPenaltyMPC_Ysd07, beliefPenaltyMPC_Lsd07);
beliefPenaltyMPC_LA_DENSE_MMTSUB_5_5(beliefPenaltyMPC_Lsd07, beliefPenaltyMPC_Yd07);
beliefPenaltyMPC_LA_DENSE_CHOL_5(beliefPenaltyMPC_Yd07, beliefPenaltyMPC_Ld07);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_5_5(beliefPenaltyMPC_Lsd07, beliefPenaltyMPC_yy06, beliefPenaltyMPC_beta07, beliefPenaltyMPC_bmy07);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_5(beliefPenaltyMPC_Ld07, beliefPenaltyMPC_bmy07, beliefPenaltyMPC_yy07);
beliefPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_5_5(beliefPenaltyMPC_Ld07, beliefPenaltyMPC_Ysd08, beliefPenaltyMPC_Lsd08);
beliefPenaltyMPC_LA_DENSE_MMTSUB_5_5(beliefPenaltyMPC_Lsd08, beliefPenaltyMPC_Yd08);
beliefPenaltyMPC_LA_DENSE_CHOL_5(beliefPenaltyMPC_Yd08, beliefPenaltyMPC_Ld08);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_5_5(beliefPenaltyMPC_Lsd08, beliefPenaltyMPC_yy07, beliefPenaltyMPC_beta08, beliefPenaltyMPC_bmy08);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_5(beliefPenaltyMPC_Ld08, beliefPenaltyMPC_bmy08, beliefPenaltyMPC_yy08);
beliefPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_5_5(beliefPenaltyMPC_Ld08, beliefPenaltyMPC_Ysd09, beliefPenaltyMPC_Lsd09);
beliefPenaltyMPC_LA_DENSE_MMTSUB_5_5(beliefPenaltyMPC_Lsd09, beliefPenaltyMPC_Yd09);
beliefPenaltyMPC_LA_DENSE_CHOL_5(beliefPenaltyMPC_Yd09, beliefPenaltyMPC_Ld09);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_5_5(beliefPenaltyMPC_Lsd09, beliefPenaltyMPC_yy08, beliefPenaltyMPC_beta09, beliefPenaltyMPC_bmy09);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_5(beliefPenaltyMPC_Ld09, beliefPenaltyMPC_bmy09, beliefPenaltyMPC_yy09);
beliefPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_5_5(beliefPenaltyMPC_Ld09, beliefPenaltyMPC_Ysd10, beliefPenaltyMPC_Lsd10);
beliefPenaltyMPC_LA_DENSE_MMTSUB_5_5(beliefPenaltyMPC_Lsd10, beliefPenaltyMPC_Yd10);
beliefPenaltyMPC_LA_DENSE_CHOL_5(beliefPenaltyMPC_Yd10, beliefPenaltyMPC_Ld10);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_5_5(beliefPenaltyMPC_Lsd10, beliefPenaltyMPC_yy09, beliefPenaltyMPC_beta10, beliefPenaltyMPC_bmy10);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_5(beliefPenaltyMPC_Ld10, beliefPenaltyMPC_bmy10, beliefPenaltyMPC_yy10);
beliefPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_5_5(beliefPenaltyMPC_Ld10, beliefPenaltyMPC_Ysd11, beliefPenaltyMPC_Lsd11);
beliefPenaltyMPC_LA_DENSE_MMTSUB_5_5(beliefPenaltyMPC_Lsd11, beliefPenaltyMPC_Yd11);
beliefPenaltyMPC_LA_DENSE_CHOL_5(beliefPenaltyMPC_Yd11, beliefPenaltyMPC_Ld11);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_5_5(beliefPenaltyMPC_Lsd11, beliefPenaltyMPC_yy10, beliefPenaltyMPC_beta11, beliefPenaltyMPC_bmy11);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_5(beliefPenaltyMPC_Ld11, beliefPenaltyMPC_bmy11, beliefPenaltyMPC_yy11);
beliefPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_5_5(beliefPenaltyMPC_Ld11, beliefPenaltyMPC_Ysd12, beliefPenaltyMPC_Lsd12);
beliefPenaltyMPC_LA_DENSE_MMTSUB_5_5(beliefPenaltyMPC_Lsd12, beliefPenaltyMPC_Yd12);
beliefPenaltyMPC_LA_DENSE_CHOL_5(beliefPenaltyMPC_Yd12, beliefPenaltyMPC_Ld12);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_5_5(beliefPenaltyMPC_Lsd12, beliefPenaltyMPC_yy11, beliefPenaltyMPC_beta12, beliefPenaltyMPC_bmy12);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_5(beliefPenaltyMPC_Ld12, beliefPenaltyMPC_bmy12, beliefPenaltyMPC_yy12);
beliefPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_5_5(beliefPenaltyMPC_Ld12, beliefPenaltyMPC_Ysd13, beliefPenaltyMPC_Lsd13);
beliefPenaltyMPC_LA_DENSE_MMTSUB_5_5(beliefPenaltyMPC_Lsd13, beliefPenaltyMPC_Yd13);
beliefPenaltyMPC_LA_DENSE_CHOL_5(beliefPenaltyMPC_Yd13, beliefPenaltyMPC_Ld13);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_5_5(beliefPenaltyMPC_Lsd13, beliefPenaltyMPC_yy12, beliefPenaltyMPC_beta13, beliefPenaltyMPC_bmy13);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_5(beliefPenaltyMPC_Ld13, beliefPenaltyMPC_bmy13, beliefPenaltyMPC_yy13);
beliefPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_5_5(beliefPenaltyMPC_Ld13, beliefPenaltyMPC_Ysd14, beliefPenaltyMPC_Lsd14);
beliefPenaltyMPC_LA_DENSE_MMTSUB_5_5(beliefPenaltyMPC_Lsd14, beliefPenaltyMPC_Yd14);
beliefPenaltyMPC_LA_DENSE_CHOL_5(beliefPenaltyMPC_Yd14, beliefPenaltyMPC_Ld14);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_5_5(beliefPenaltyMPC_Lsd14, beliefPenaltyMPC_yy13, beliefPenaltyMPC_beta14, beliefPenaltyMPC_bmy14);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_5(beliefPenaltyMPC_Ld14, beliefPenaltyMPC_bmy14, beliefPenaltyMPC_yy14);
beliefPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_5_5(beliefPenaltyMPC_Ld14, beliefPenaltyMPC_Ysd15, beliefPenaltyMPC_Lsd15);
beliefPenaltyMPC_LA_DENSE_MMTSUB_5_5(beliefPenaltyMPC_Lsd15, beliefPenaltyMPC_Yd15);
beliefPenaltyMPC_LA_DENSE_CHOL_5(beliefPenaltyMPC_Yd15, beliefPenaltyMPC_Ld15);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_5_5(beliefPenaltyMPC_Lsd15, beliefPenaltyMPC_yy14, beliefPenaltyMPC_beta15, beliefPenaltyMPC_bmy15);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_5(beliefPenaltyMPC_Ld15, beliefPenaltyMPC_bmy15, beliefPenaltyMPC_yy15);
beliefPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_5_5(beliefPenaltyMPC_Ld15, beliefPenaltyMPC_Ysd16, beliefPenaltyMPC_Lsd16);
beliefPenaltyMPC_LA_DENSE_MMTSUB_5_5(beliefPenaltyMPC_Lsd16, beliefPenaltyMPC_Yd16);
beliefPenaltyMPC_LA_DENSE_CHOL_5(beliefPenaltyMPC_Yd16, beliefPenaltyMPC_Ld16);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_5_5(beliefPenaltyMPC_Lsd16, beliefPenaltyMPC_yy15, beliefPenaltyMPC_beta16, beliefPenaltyMPC_bmy16);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_5(beliefPenaltyMPC_Ld16, beliefPenaltyMPC_bmy16, beliefPenaltyMPC_yy16);
beliefPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_5_5(beliefPenaltyMPC_Ld16, beliefPenaltyMPC_Ysd17, beliefPenaltyMPC_Lsd17);
beliefPenaltyMPC_LA_DENSE_MMTSUB_5_5(beliefPenaltyMPC_Lsd17, beliefPenaltyMPC_Yd17);
beliefPenaltyMPC_LA_DENSE_CHOL_5(beliefPenaltyMPC_Yd17, beliefPenaltyMPC_Ld17);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_5_5(beliefPenaltyMPC_Lsd17, beliefPenaltyMPC_yy16, beliefPenaltyMPC_beta17, beliefPenaltyMPC_bmy17);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_5(beliefPenaltyMPC_Ld17, beliefPenaltyMPC_bmy17, beliefPenaltyMPC_yy17);
beliefPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_5_5(beliefPenaltyMPC_Ld17, beliefPenaltyMPC_Ysd18, beliefPenaltyMPC_Lsd18);
beliefPenaltyMPC_LA_DENSE_MMTSUB_5_5(beliefPenaltyMPC_Lsd18, beliefPenaltyMPC_Yd18);
beliefPenaltyMPC_LA_DENSE_CHOL_5(beliefPenaltyMPC_Yd18, beliefPenaltyMPC_Ld18);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_5_5(beliefPenaltyMPC_Lsd18, beliefPenaltyMPC_yy17, beliefPenaltyMPC_beta18, beliefPenaltyMPC_bmy18);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_5(beliefPenaltyMPC_Ld18, beliefPenaltyMPC_bmy18, beliefPenaltyMPC_yy18);
beliefPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_5_5(beliefPenaltyMPC_Ld18, beliefPenaltyMPC_Ysd19, beliefPenaltyMPC_Lsd19);
beliefPenaltyMPC_LA_DENSE_MMTSUB_5_5(beliefPenaltyMPC_Lsd19, beliefPenaltyMPC_Yd19);
beliefPenaltyMPC_LA_DENSE_CHOL_5(beliefPenaltyMPC_Yd19, beliefPenaltyMPC_Ld19);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_5_5(beliefPenaltyMPC_Lsd19, beliefPenaltyMPC_yy18, beliefPenaltyMPC_beta19, beliefPenaltyMPC_bmy19);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_5(beliefPenaltyMPC_Ld19, beliefPenaltyMPC_bmy19, beliefPenaltyMPC_yy19);
beliefPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_5_5(beliefPenaltyMPC_Ld19, beliefPenaltyMPC_Ysd20, beliefPenaltyMPC_Lsd20);
beliefPenaltyMPC_LA_DENSE_MMTSUB_5_5(beliefPenaltyMPC_Lsd20, beliefPenaltyMPC_Yd20);
beliefPenaltyMPC_LA_DENSE_CHOL_5(beliefPenaltyMPC_Yd20, beliefPenaltyMPC_Ld20);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_5_5(beliefPenaltyMPC_Lsd20, beliefPenaltyMPC_yy19, beliefPenaltyMPC_beta20, beliefPenaltyMPC_bmy20);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_5(beliefPenaltyMPC_Ld20, beliefPenaltyMPC_bmy20, beliefPenaltyMPC_yy20);
beliefPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_5_5(beliefPenaltyMPC_Ld20, beliefPenaltyMPC_Ysd21, beliefPenaltyMPC_Lsd21);
beliefPenaltyMPC_LA_DENSE_MMTSUB_5_5(beliefPenaltyMPC_Lsd21, beliefPenaltyMPC_Yd21);
beliefPenaltyMPC_LA_DENSE_CHOL_5(beliefPenaltyMPC_Yd21, beliefPenaltyMPC_Ld21);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_5_5(beliefPenaltyMPC_Lsd21, beliefPenaltyMPC_yy20, beliefPenaltyMPC_beta21, beliefPenaltyMPC_bmy21);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_5(beliefPenaltyMPC_Ld21, beliefPenaltyMPC_bmy21, beliefPenaltyMPC_yy21);
beliefPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_5_5(beliefPenaltyMPC_Ld21, beliefPenaltyMPC_Ysd22, beliefPenaltyMPC_Lsd22);
beliefPenaltyMPC_LA_DENSE_MMTSUB_5_5(beliefPenaltyMPC_Lsd22, beliefPenaltyMPC_Yd22);
beliefPenaltyMPC_LA_DENSE_CHOL_5(beliefPenaltyMPC_Yd22, beliefPenaltyMPC_Ld22);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_5_5(beliefPenaltyMPC_Lsd22, beliefPenaltyMPC_yy21, beliefPenaltyMPC_beta22, beliefPenaltyMPC_bmy22);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_5(beliefPenaltyMPC_Ld22, beliefPenaltyMPC_bmy22, beliefPenaltyMPC_yy22);
beliefPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_5_5(beliefPenaltyMPC_Ld22, beliefPenaltyMPC_Ysd23, beliefPenaltyMPC_Lsd23);
beliefPenaltyMPC_LA_DENSE_MMTSUB_5_5(beliefPenaltyMPC_Lsd23, beliefPenaltyMPC_Yd23);
beliefPenaltyMPC_LA_DENSE_CHOL_5(beliefPenaltyMPC_Yd23, beliefPenaltyMPC_Ld23);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_5_5(beliefPenaltyMPC_Lsd23, beliefPenaltyMPC_yy22, beliefPenaltyMPC_beta23, beliefPenaltyMPC_bmy23);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_5(beliefPenaltyMPC_Ld23, beliefPenaltyMPC_bmy23, beliefPenaltyMPC_yy23);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_5(beliefPenaltyMPC_Ld23, beliefPenaltyMPC_yy23, beliefPenaltyMPC_dvaff23);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_5_5(beliefPenaltyMPC_Lsd23, beliefPenaltyMPC_dvaff23, beliefPenaltyMPC_yy22, beliefPenaltyMPC_bmy22);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_5(beliefPenaltyMPC_Ld22, beliefPenaltyMPC_bmy22, beliefPenaltyMPC_dvaff22);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_5_5(beliefPenaltyMPC_Lsd22, beliefPenaltyMPC_dvaff22, beliefPenaltyMPC_yy21, beliefPenaltyMPC_bmy21);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_5(beliefPenaltyMPC_Ld21, beliefPenaltyMPC_bmy21, beliefPenaltyMPC_dvaff21);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_5_5(beliefPenaltyMPC_Lsd21, beliefPenaltyMPC_dvaff21, beliefPenaltyMPC_yy20, beliefPenaltyMPC_bmy20);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_5(beliefPenaltyMPC_Ld20, beliefPenaltyMPC_bmy20, beliefPenaltyMPC_dvaff20);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_5_5(beliefPenaltyMPC_Lsd20, beliefPenaltyMPC_dvaff20, beliefPenaltyMPC_yy19, beliefPenaltyMPC_bmy19);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_5(beliefPenaltyMPC_Ld19, beliefPenaltyMPC_bmy19, beliefPenaltyMPC_dvaff19);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_5_5(beliefPenaltyMPC_Lsd19, beliefPenaltyMPC_dvaff19, beliefPenaltyMPC_yy18, beliefPenaltyMPC_bmy18);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_5(beliefPenaltyMPC_Ld18, beliefPenaltyMPC_bmy18, beliefPenaltyMPC_dvaff18);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_5_5(beliefPenaltyMPC_Lsd18, beliefPenaltyMPC_dvaff18, beliefPenaltyMPC_yy17, beliefPenaltyMPC_bmy17);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_5(beliefPenaltyMPC_Ld17, beliefPenaltyMPC_bmy17, beliefPenaltyMPC_dvaff17);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_5_5(beliefPenaltyMPC_Lsd17, beliefPenaltyMPC_dvaff17, beliefPenaltyMPC_yy16, beliefPenaltyMPC_bmy16);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_5(beliefPenaltyMPC_Ld16, beliefPenaltyMPC_bmy16, beliefPenaltyMPC_dvaff16);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_5_5(beliefPenaltyMPC_Lsd16, beliefPenaltyMPC_dvaff16, beliefPenaltyMPC_yy15, beliefPenaltyMPC_bmy15);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_5(beliefPenaltyMPC_Ld15, beliefPenaltyMPC_bmy15, beliefPenaltyMPC_dvaff15);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_5_5(beliefPenaltyMPC_Lsd15, beliefPenaltyMPC_dvaff15, beliefPenaltyMPC_yy14, beliefPenaltyMPC_bmy14);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_5(beliefPenaltyMPC_Ld14, beliefPenaltyMPC_bmy14, beliefPenaltyMPC_dvaff14);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_5_5(beliefPenaltyMPC_Lsd14, beliefPenaltyMPC_dvaff14, beliefPenaltyMPC_yy13, beliefPenaltyMPC_bmy13);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_5(beliefPenaltyMPC_Ld13, beliefPenaltyMPC_bmy13, beliefPenaltyMPC_dvaff13);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_5_5(beliefPenaltyMPC_Lsd13, beliefPenaltyMPC_dvaff13, beliefPenaltyMPC_yy12, beliefPenaltyMPC_bmy12);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_5(beliefPenaltyMPC_Ld12, beliefPenaltyMPC_bmy12, beliefPenaltyMPC_dvaff12);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_5_5(beliefPenaltyMPC_Lsd12, beliefPenaltyMPC_dvaff12, beliefPenaltyMPC_yy11, beliefPenaltyMPC_bmy11);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_5(beliefPenaltyMPC_Ld11, beliefPenaltyMPC_bmy11, beliefPenaltyMPC_dvaff11);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_5_5(beliefPenaltyMPC_Lsd11, beliefPenaltyMPC_dvaff11, beliefPenaltyMPC_yy10, beliefPenaltyMPC_bmy10);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_5(beliefPenaltyMPC_Ld10, beliefPenaltyMPC_bmy10, beliefPenaltyMPC_dvaff10);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_5_5(beliefPenaltyMPC_Lsd10, beliefPenaltyMPC_dvaff10, beliefPenaltyMPC_yy09, beliefPenaltyMPC_bmy09);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_5(beliefPenaltyMPC_Ld09, beliefPenaltyMPC_bmy09, beliefPenaltyMPC_dvaff09);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_5_5(beliefPenaltyMPC_Lsd09, beliefPenaltyMPC_dvaff09, beliefPenaltyMPC_yy08, beliefPenaltyMPC_bmy08);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_5(beliefPenaltyMPC_Ld08, beliefPenaltyMPC_bmy08, beliefPenaltyMPC_dvaff08);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_5_5(beliefPenaltyMPC_Lsd08, beliefPenaltyMPC_dvaff08, beliefPenaltyMPC_yy07, beliefPenaltyMPC_bmy07);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_5(beliefPenaltyMPC_Ld07, beliefPenaltyMPC_bmy07, beliefPenaltyMPC_dvaff07);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_5_5(beliefPenaltyMPC_Lsd07, beliefPenaltyMPC_dvaff07, beliefPenaltyMPC_yy06, beliefPenaltyMPC_bmy06);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_5(beliefPenaltyMPC_Ld06, beliefPenaltyMPC_bmy06, beliefPenaltyMPC_dvaff06);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_5_5(beliefPenaltyMPC_Lsd06, beliefPenaltyMPC_dvaff06, beliefPenaltyMPC_yy05, beliefPenaltyMPC_bmy05);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_5(beliefPenaltyMPC_Ld05, beliefPenaltyMPC_bmy05, beliefPenaltyMPC_dvaff05);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_5_5(beliefPenaltyMPC_Lsd05, beliefPenaltyMPC_dvaff05, beliefPenaltyMPC_yy04, beliefPenaltyMPC_bmy04);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_5(beliefPenaltyMPC_Ld04, beliefPenaltyMPC_bmy04, beliefPenaltyMPC_dvaff04);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_5_5(beliefPenaltyMPC_Lsd04, beliefPenaltyMPC_dvaff04, beliefPenaltyMPC_yy03, beliefPenaltyMPC_bmy03);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_5(beliefPenaltyMPC_Ld03, beliefPenaltyMPC_bmy03, beliefPenaltyMPC_dvaff03);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_5_5(beliefPenaltyMPC_Lsd03, beliefPenaltyMPC_dvaff03, beliefPenaltyMPC_yy02, beliefPenaltyMPC_bmy02);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_5(beliefPenaltyMPC_Ld02, beliefPenaltyMPC_bmy02, beliefPenaltyMPC_dvaff02);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_5_5(beliefPenaltyMPC_Lsd02, beliefPenaltyMPC_dvaff02, beliefPenaltyMPC_yy01, beliefPenaltyMPC_bmy01);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_5(beliefPenaltyMPC_Ld01, beliefPenaltyMPC_bmy01, beliefPenaltyMPC_dvaff01);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_5_10(beliefPenaltyMPC_Lsd01, beliefPenaltyMPC_dvaff01, beliefPenaltyMPC_yy00, beliefPenaltyMPC_bmy00);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_10(beliefPenaltyMPC_Ld00, beliefPenaltyMPC_bmy00, beliefPenaltyMPC_dvaff00);
beliefPenaltyMPC_LA_DENSE_MTVM_10_17(params->C1, beliefPenaltyMPC_dvaff00, beliefPenaltyMPC_grad_eq00);
beliefPenaltyMPC_LA_DENSE_MTVM2_5_17_10(params->C2, beliefPenaltyMPC_dvaff01, beliefPenaltyMPC_D01, beliefPenaltyMPC_dvaff00, beliefPenaltyMPC_grad_eq01);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_5_17_5(params->C3, beliefPenaltyMPC_dvaff02, beliefPenaltyMPC_D02, beliefPenaltyMPC_dvaff01, beliefPenaltyMPC_grad_eq02);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_5_17_5(params->C4, beliefPenaltyMPC_dvaff03, beliefPenaltyMPC_D02, beliefPenaltyMPC_dvaff02, beliefPenaltyMPC_grad_eq03);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_5_17_5(params->C5, beliefPenaltyMPC_dvaff04, beliefPenaltyMPC_D02, beliefPenaltyMPC_dvaff03, beliefPenaltyMPC_grad_eq04);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_5_17_5(params->C6, beliefPenaltyMPC_dvaff05, beliefPenaltyMPC_D02, beliefPenaltyMPC_dvaff04, beliefPenaltyMPC_grad_eq05);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_5_17_5(params->C7, beliefPenaltyMPC_dvaff06, beliefPenaltyMPC_D02, beliefPenaltyMPC_dvaff05, beliefPenaltyMPC_grad_eq06);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_5_17_5(params->C8, beliefPenaltyMPC_dvaff07, beliefPenaltyMPC_D02, beliefPenaltyMPC_dvaff06, beliefPenaltyMPC_grad_eq07);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_5_17_5(params->C9, beliefPenaltyMPC_dvaff08, beliefPenaltyMPC_D02, beliefPenaltyMPC_dvaff07, beliefPenaltyMPC_grad_eq08);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_5_17_5(params->C10, beliefPenaltyMPC_dvaff09, beliefPenaltyMPC_D02, beliefPenaltyMPC_dvaff08, beliefPenaltyMPC_grad_eq09);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_5_17_5(params->C11, beliefPenaltyMPC_dvaff10, beliefPenaltyMPC_D02, beliefPenaltyMPC_dvaff09, beliefPenaltyMPC_grad_eq10);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_5_17_5(params->C12, beliefPenaltyMPC_dvaff11, beliefPenaltyMPC_D02, beliefPenaltyMPC_dvaff10, beliefPenaltyMPC_grad_eq11);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_5_17_5(params->C13, beliefPenaltyMPC_dvaff12, beliefPenaltyMPC_D02, beliefPenaltyMPC_dvaff11, beliefPenaltyMPC_grad_eq12);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_5_17_5(params->C14, beliefPenaltyMPC_dvaff13, beliefPenaltyMPC_D02, beliefPenaltyMPC_dvaff12, beliefPenaltyMPC_grad_eq13);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_5_17_5(params->C15, beliefPenaltyMPC_dvaff14, beliefPenaltyMPC_D02, beliefPenaltyMPC_dvaff13, beliefPenaltyMPC_grad_eq14);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_5_17_5(params->C16, beliefPenaltyMPC_dvaff15, beliefPenaltyMPC_D02, beliefPenaltyMPC_dvaff14, beliefPenaltyMPC_grad_eq15);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_5_17_5(params->C17, beliefPenaltyMPC_dvaff16, beliefPenaltyMPC_D02, beliefPenaltyMPC_dvaff15, beliefPenaltyMPC_grad_eq16);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_5_17_5(params->C18, beliefPenaltyMPC_dvaff17, beliefPenaltyMPC_D02, beliefPenaltyMPC_dvaff16, beliefPenaltyMPC_grad_eq17);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_5_17_5(params->C19, beliefPenaltyMPC_dvaff18, beliefPenaltyMPC_D02, beliefPenaltyMPC_dvaff17, beliefPenaltyMPC_grad_eq18);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_5_17_5(params->C20, beliefPenaltyMPC_dvaff19, beliefPenaltyMPC_D02, beliefPenaltyMPC_dvaff18, beliefPenaltyMPC_grad_eq19);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_5_17_5(params->C21, beliefPenaltyMPC_dvaff20, beliefPenaltyMPC_D02, beliefPenaltyMPC_dvaff19, beliefPenaltyMPC_grad_eq20);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_5_17_5(params->C22, beliefPenaltyMPC_dvaff21, beliefPenaltyMPC_D02, beliefPenaltyMPC_dvaff20, beliefPenaltyMPC_grad_eq21);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_5_17_5(params->C23, beliefPenaltyMPC_dvaff22, beliefPenaltyMPC_D02, beliefPenaltyMPC_dvaff21, beliefPenaltyMPC_grad_eq22);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_5_17_5(params->C24, beliefPenaltyMPC_dvaff23, beliefPenaltyMPC_D02, beliefPenaltyMPC_dvaff22, beliefPenaltyMPC_grad_eq23);
beliefPenaltyMPC_LA_DIAGZERO_MTVM_5_5(beliefPenaltyMPC_D24, beliefPenaltyMPC_dvaff23, beliefPenaltyMPC_grad_eq24);
beliefPenaltyMPC_LA_VSUB2_413(beliefPenaltyMPC_rd, beliefPenaltyMPC_grad_eq, beliefPenaltyMPC_rd);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_17(beliefPenaltyMPC_Phi00, beliefPenaltyMPC_rd00, beliefPenaltyMPC_dzaff00);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_17(beliefPenaltyMPC_Phi01, beliefPenaltyMPC_rd01, beliefPenaltyMPC_dzaff01);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_17(beliefPenaltyMPC_Phi02, beliefPenaltyMPC_rd02, beliefPenaltyMPC_dzaff02);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_17(beliefPenaltyMPC_Phi03, beliefPenaltyMPC_rd03, beliefPenaltyMPC_dzaff03);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_17(beliefPenaltyMPC_Phi04, beliefPenaltyMPC_rd04, beliefPenaltyMPC_dzaff04);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_17(beliefPenaltyMPC_Phi05, beliefPenaltyMPC_rd05, beliefPenaltyMPC_dzaff05);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_17(beliefPenaltyMPC_Phi06, beliefPenaltyMPC_rd06, beliefPenaltyMPC_dzaff06);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_17(beliefPenaltyMPC_Phi07, beliefPenaltyMPC_rd07, beliefPenaltyMPC_dzaff07);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_17(beliefPenaltyMPC_Phi08, beliefPenaltyMPC_rd08, beliefPenaltyMPC_dzaff08);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_17(beliefPenaltyMPC_Phi09, beliefPenaltyMPC_rd09, beliefPenaltyMPC_dzaff09);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_17(beliefPenaltyMPC_Phi10, beliefPenaltyMPC_rd10, beliefPenaltyMPC_dzaff10);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_17(beliefPenaltyMPC_Phi11, beliefPenaltyMPC_rd11, beliefPenaltyMPC_dzaff11);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_17(beliefPenaltyMPC_Phi12, beliefPenaltyMPC_rd12, beliefPenaltyMPC_dzaff12);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_17(beliefPenaltyMPC_Phi13, beliefPenaltyMPC_rd13, beliefPenaltyMPC_dzaff13);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_17(beliefPenaltyMPC_Phi14, beliefPenaltyMPC_rd14, beliefPenaltyMPC_dzaff14);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_17(beliefPenaltyMPC_Phi15, beliefPenaltyMPC_rd15, beliefPenaltyMPC_dzaff15);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_17(beliefPenaltyMPC_Phi16, beliefPenaltyMPC_rd16, beliefPenaltyMPC_dzaff16);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_17(beliefPenaltyMPC_Phi17, beliefPenaltyMPC_rd17, beliefPenaltyMPC_dzaff17);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_17(beliefPenaltyMPC_Phi18, beliefPenaltyMPC_rd18, beliefPenaltyMPC_dzaff18);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_17(beliefPenaltyMPC_Phi19, beliefPenaltyMPC_rd19, beliefPenaltyMPC_dzaff19);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_17(beliefPenaltyMPC_Phi20, beliefPenaltyMPC_rd20, beliefPenaltyMPC_dzaff20);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_17(beliefPenaltyMPC_Phi21, beliefPenaltyMPC_rd21, beliefPenaltyMPC_dzaff21);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_17(beliefPenaltyMPC_Phi22, beliefPenaltyMPC_rd22, beliefPenaltyMPC_dzaff22);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_17(beliefPenaltyMPC_Phi23, beliefPenaltyMPC_rd23, beliefPenaltyMPC_dzaff23);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_5(beliefPenaltyMPC_Phi24, beliefPenaltyMPC_rd24, beliefPenaltyMPC_dzaff24);
beliefPenaltyMPC_LA_VSUB_INDEXED_17(beliefPenaltyMPC_dzaff00, beliefPenaltyMPC_lbIdx00, beliefPenaltyMPC_rilb00, beliefPenaltyMPC_dslbaff00);
beliefPenaltyMPC_LA_VSUB3_17(beliefPenaltyMPC_llbbyslb00, beliefPenaltyMPC_dslbaff00, beliefPenaltyMPC_llb00, beliefPenaltyMPC_dllbaff00);
beliefPenaltyMPC_LA_VSUB2_INDEXED_7(beliefPenaltyMPC_riub00, beliefPenaltyMPC_dzaff00, beliefPenaltyMPC_ubIdx00, beliefPenaltyMPC_dsubaff00);
beliefPenaltyMPC_LA_VSUB3_7(beliefPenaltyMPC_lubbysub00, beliefPenaltyMPC_dsubaff00, beliefPenaltyMPC_lub00, beliefPenaltyMPC_dlubaff00);
beliefPenaltyMPC_LA_VSUB_INDEXED_17(beliefPenaltyMPC_dzaff01, beliefPenaltyMPC_lbIdx01, beliefPenaltyMPC_rilb01, beliefPenaltyMPC_dslbaff01);
beliefPenaltyMPC_LA_VSUB3_17(beliefPenaltyMPC_llbbyslb01, beliefPenaltyMPC_dslbaff01, beliefPenaltyMPC_llb01, beliefPenaltyMPC_dllbaff01);
beliefPenaltyMPC_LA_VSUB2_INDEXED_7(beliefPenaltyMPC_riub01, beliefPenaltyMPC_dzaff01, beliefPenaltyMPC_ubIdx01, beliefPenaltyMPC_dsubaff01);
beliefPenaltyMPC_LA_VSUB3_7(beliefPenaltyMPC_lubbysub01, beliefPenaltyMPC_dsubaff01, beliefPenaltyMPC_lub01, beliefPenaltyMPC_dlubaff01);
beliefPenaltyMPC_LA_VSUB_INDEXED_17(beliefPenaltyMPC_dzaff02, beliefPenaltyMPC_lbIdx02, beliefPenaltyMPC_rilb02, beliefPenaltyMPC_dslbaff02);
beliefPenaltyMPC_LA_VSUB3_17(beliefPenaltyMPC_llbbyslb02, beliefPenaltyMPC_dslbaff02, beliefPenaltyMPC_llb02, beliefPenaltyMPC_dllbaff02);
beliefPenaltyMPC_LA_VSUB2_INDEXED_7(beliefPenaltyMPC_riub02, beliefPenaltyMPC_dzaff02, beliefPenaltyMPC_ubIdx02, beliefPenaltyMPC_dsubaff02);
beliefPenaltyMPC_LA_VSUB3_7(beliefPenaltyMPC_lubbysub02, beliefPenaltyMPC_dsubaff02, beliefPenaltyMPC_lub02, beliefPenaltyMPC_dlubaff02);
beliefPenaltyMPC_LA_VSUB_INDEXED_17(beliefPenaltyMPC_dzaff03, beliefPenaltyMPC_lbIdx03, beliefPenaltyMPC_rilb03, beliefPenaltyMPC_dslbaff03);
beliefPenaltyMPC_LA_VSUB3_17(beliefPenaltyMPC_llbbyslb03, beliefPenaltyMPC_dslbaff03, beliefPenaltyMPC_llb03, beliefPenaltyMPC_dllbaff03);
beliefPenaltyMPC_LA_VSUB2_INDEXED_7(beliefPenaltyMPC_riub03, beliefPenaltyMPC_dzaff03, beliefPenaltyMPC_ubIdx03, beliefPenaltyMPC_dsubaff03);
beliefPenaltyMPC_LA_VSUB3_7(beliefPenaltyMPC_lubbysub03, beliefPenaltyMPC_dsubaff03, beliefPenaltyMPC_lub03, beliefPenaltyMPC_dlubaff03);
beliefPenaltyMPC_LA_VSUB_INDEXED_17(beliefPenaltyMPC_dzaff04, beliefPenaltyMPC_lbIdx04, beliefPenaltyMPC_rilb04, beliefPenaltyMPC_dslbaff04);
beliefPenaltyMPC_LA_VSUB3_17(beliefPenaltyMPC_llbbyslb04, beliefPenaltyMPC_dslbaff04, beliefPenaltyMPC_llb04, beliefPenaltyMPC_dllbaff04);
beliefPenaltyMPC_LA_VSUB2_INDEXED_7(beliefPenaltyMPC_riub04, beliefPenaltyMPC_dzaff04, beliefPenaltyMPC_ubIdx04, beliefPenaltyMPC_dsubaff04);
beliefPenaltyMPC_LA_VSUB3_7(beliefPenaltyMPC_lubbysub04, beliefPenaltyMPC_dsubaff04, beliefPenaltyMPC_lub04, beliefPenaltyMPC_dlubaff04);
beliefPenaltyMPC_LA_VSUB_INDEXED_17(beliefPenaltyMPC_dzaff05, beliefPenaltyMPC_lbIdx05, beliefPenaltyMPC_rilb05, beliefPenaltyMPC_dslbaff05);
beliefPenaltyMPC_LA_VSUB3_17(beliefPenaltyMPC_llbbyslb05, beliefPenaltyMPC_dslbaff05, beliefPenaltyMPC_llb05, beliefPenaltyMPC_dllbaff05);
beliefPenaltyMPC_LA_VSUB2_INDEXED_7(beliefPenaltyMPC_riub05, beliefPenaltyMPC_dzaff05, beliefPenaltyMPC_ubIdx05, beliefPenaltyMPC_dsubaff05);
beliefPenaltyMPC_LA_VSUB3_7(beliefPenaltyMPC_lubbysub05, beliefPenaltyMPC_dsubaff05, beliefPenaltyMPC_lub05, beliefPenaltyMPC_dlubaff05);
beliefPenaltyMPC_LA_VSUB_INDEXED_17(beliefPenaltyMPC_dzaff06, beliefPenaltyMPC_lbIdx06, beliefPenaltyMPC_rilb06, beliefPenaltyMPC_dslbaff06);
beliefPenaltyMPC_LA_VSUB3_17(beliefPenaltyMPC_llbbyslb06, beliefPenaltyMPC_dslbaff06, beliefPenaltyMPC_llb06, beliefPenaltyMPC_dllbaff06);
beliefPenaltyMPC_LA_VSUB2_INDEXED_7(beliefPenaltyMPC_riub06, beliefPenaltyMPC_dzaff06, beliefPenaltyMPC_ubIdx06, beliefPenaltyMPC_dsubaff06);
beliefPenaltyMPC_LA_VSUB3_7(beliefPenaltyMPC_lubbysub06, beliefPenaltyMPC_dsubaff06, beliefPenaltyMPC_lub06, beliefPenaltyMPC_dlubaff06);
beliefPenaltyMPC_LA_VSUB_INDEXED_17(beliefPenaltyMPC_dzaff07, beliefPenaltyMPC_lbIdx07, beliefPenaltyMPC_rilb07, beliefPenaltyMPC_dslbaff07);
beliefPenaltyMPC_LA_VSUB3_17(beliefPenaltyMPC_llbbyslb07, beliefPenaltyMPC_dslbaff07, beliefPenaltyMPC_llb07, beliefPenaltyMPC_dllbaff07);
beliefPenaltyMPC_LA_VSUB2_INDEXED_7(beliefPenaltyMPC_riub07, beliefPenaltyMPC_dzaff07, beliefPenaltyMPC_ubIdx07, beliefPenaltyMPC_dsubaff07);
beliefPenaltyMPC_LA_VSUB3_7(beliefPenaltyMPC_lubbysub07, beliefPenaltyMPC_dsubaff07, beliefPenaltyMPC_lub07, beliefPenaltyMPC_dlubaff07);
beliefPenaltyMPC_LA_VSUB_INDEXED_17(beliefPenaltyMPC_dzaff08, beliefPenaltyMPC_lbIdx08, beliefPenaltyMPC_rilb08, beliefPenaltyMPC_dslbaff08);
beliefPenaltyMPC_LA_VSUB3_17(beliefPenaltyMPC_llbbyslb08, beliefPenaltyMPC_dslbaff08, beliefPenaltyMPC_llb08, beliefPenaltyMPC_dllbaff08);
beliefPenaltyMPC_LA_VSUB2_INDEXED_7(beliefPenaltyMPC_riub08, beliefPenaltyMPC_dzaff08, beliefPenaltyMPC_ubIdx08, beliefPenaltyMPC_dsubaff08);
beliefPenaltyMPC_LA_VSUB3_7(beliefPenaltyMPC_lubbysub08, beliefPenaltyMPC_dsubaff08, beliefPenaltyMPC_lub08, beliefPenaltyMPC_dlubaff08);
beliefPenaltyMPC_LA_VSUB_INDEXED_17(beliefPenaltyMPC_dzaff09, beliefPenaltyMPC_lbIdx09, beliefPenaltyMPC_rilb09, beliefPenaltyMPC_dslbaff09);
beliefPenaltyMPC_LA_VSUB3_17(beliefPenaltyMPC_llbbyslb09, beliefPenaltyMPC_dslbaff09, beliefPenaltyMPC_llb09, beliefPenaltyMPC_dllbaff09);
beliefPenaltyMPC_LA_VSUB2_INDEXED_7(beliefPenaltyMPC_riub09, beliefPenaltyMPC_dzaff09, beliefPenaltyMPC_ubIdx09, beliefPenaltyMPC_dsubaff09);
beliefPenaltyMPC_LA_VSUB3_7(beliefPenaltyMPC_lubbysub09, beliefPenaltyMPC_dsubaff09, beliefPenaltyMPC_lub09, beliefPenaltyMPC_dlubaff09);
beliefPenaltyMPC_LA_VSUB_INDEXED_17(beliefPenaltyMPC_dzaff10, beliefPenaltyMPC_lbIdx10, beliefPenaltyMPC_rilb10, beliefPenaltyMPC_dslbaff10);
beliefPenaltyMPC_LA_VSUB3_17(beliefPenaltyMPC_llbbyslb10, beliefPenaltyMPC_dslbaff10, beliefPenaltyMPC_llb10, beliefPenaltyMPC_dllbaff10);
beliefPenaltyMPC_LA_VSUB2_INDEXED_7(beliefPenaltyMPC_riub10, beliefPenaltyMPC_dzaff10, beliefPenaltyMPC_ubIdx10, beliefPenaltyMPC_dsubaff10);
beliefPenaltyMPC_LA_VSUB3_7(beliefPenaltyMPC_lubbysub10, beliefPenaltyMPC_dsubaff10, beliefPenaltyMPC_lub10, beliefPenaltyMPC_dlubaff10);
beliefPenaltyMPC_LA_VSUB_INDEXED_17(beliefPenaltyMPC_dzaff11, beliefPenaltyMPC_lbIdx11, beliefPenaltyMPC_rilb11, beliefPenaltyMPC_dslbaff11);
beliefPenaltyMPC_LA_VSUB3_17(beliefPenaltyMPC_llbbyslb11, beliefPenaltyMPC_dslbaff11, beliefPenaltyMPC_llb11, beliefPenaltyMPC_dllbaff11);
beliefPenaltyMPC_LA_VSUB2_INDEXED_7(beliefPenaltyMPC_riub11, beliefPenaltyMPC_dzaff11, beliefPenaltyMPC_ubIdx11, beliefPenaltyMPC_dsubaff11);
beliefPenaltyMPC_LA_VSUB3_7(beliefPenaltyMPC_lubbysub11, beliefPenaltyMPC_dsubaff11, beliefPenaltyMPC_lub11, beliefPenaltyMPC_dlubaff11);
beliefPenaltyMPC_LA_VSUB_INDEXED_17(beliefPenaltyMPC_dzaff12, beliefPenaltyMPC_lbIdx12, beliefPenaltyMPC_rilb12, beliefPenaltyMPC_dslbaff12);
beliefPenaltyMPC_LA_VSUB3_17(beliefPenaltyMPC_llbbyslb12, beliefPenaltyMPC_dslbaff12, beliefPenaltyMPC_llb12, beliefPenaltyMPC_dllbaff12);
beliefPenaltyMPC_LA_VSUB2_INDEXED_7(beliefPenaltyMPC_riub12, beliefPenaltyMPC_dzaff12, beliefPenaltyMPC_ubIdx12, beliefPenaltyMPC_dsubaff12);
beliefPenaltyMPC_LA_VSUB3_7(beliefPenaltyMPC_lubbysub12, beliefPenaltyMPC_dsubaff12, beliefPenaltyMPC_lub12, beliefPenaltyMPC_dlubaff12);
beliefPenaltyMPC_LA_VSUB_INDEXED_17(beliefPenaltyMPC_dzaff13, beliefPenaltyMPC_lbIdx13, beliefPenaltyMPC_rilb13, beliefPenaltyMPC_dslbaff13);
beliefPenaltyMPC_LA_VSUB3_17(beliefPenaltyMPC_llbbyslb13, beliefPenaltyMPC_dslbaff13, beliefPenaltyMPC_llb13, beliefPenaltyMPC_dllbaff13);
beliefPenaltyMPC_LA_VSUB2_INDEXED_7(beliefPenaltyMPC_riub13, beliefPenaltyMPC_dzaff13, beliefPenaltyMPC_ubIdx13, beliefPenaltyMPC_dsubaff13);
beliefPenaltyMPC_LA_VSUB3_7(beliefPenaltyMPC_lubbysub13, beliefPenaltyMPC_dsubaff13, beliefPenaltyMPC_lub13, beliefPenaltyMPC_dlubaff13);
beliefPenaltyMPC_LA_VSUB_INDEXED_17(beliefPenaltyMPC_dzaff14, beliefPenaltyMPC_lbIdx14, beliefPenaltyMPC_rilb14, beliefPenaltyMPC_dslbaff14);
beliefPenaltyMPC_LA_VSUB3_17(beliefPenaltyMPC_llbbyslb14, beliefPenaltyMPC_dslbaff14, beliefPenaltyMPC_llb14, beliefPenaltyMPC_dllbaff14);
beliefPenaltyMPC_LA_VSUB2_INDEXED_7(beliefPenaltyMPC_riub14, beliefPenaltyMPC_dzaff14, beliefPenaltyMPC_ubIdx14, beliefPenaltyMPC_dsubaff14);
beliefPenaltyMPC_LA_VSUB3_7(beliefPenaltyMPC_lubbysub14, beliefPenaltyMPC_dsubaff14, beliefPenaltyMPC_lub14, beliefPenaltyMPC_dlubaff14);
beliefPenaltyMPC_LA_VSUB_INDEXED_17(beliefPenaltyMPC_dzaff15, beliefPenaltyMPC_lbIdx15, beliefPenaltyMPC_rilb15, beliefPenaltyMPC_dslbaff15);
beliefPenaltyMPC_LA_VSUB3_17(beliefPenaltyMPC_llbbyslb15, beliefPenaltyMPC_dslbaff15, beliefPenaltyMPC_llb15, beliefPenaltyMPC_dllbaff15);
beliefPenaltyMPC_LA_VSUB2_INDEXED_7(beliefPenaltyMPC_riub15, beliefPenaltyMPC_dzaff15, beliefPenaltyMPC_ubIdx15, beliefPenaltyMPC_dsubaff15);
beliefPenaltyMPC_LA_VSUB3_7(beliefPenaltyMPC_lubbysub15, beliefPenaltyMPC_dsubaff15, beliefPenaltyMPC_lub15, beliefPenaltyMPC_dlubaff15);
beliefPenaltyMPC_LA_VSUB_INDEXED_17(beliefPenaltyMPC_dzaff16, beliefPenaltyMPC_lbIdx16, beliefPenaltyMPC_rilb16, beliefPenaltyMPC_dslbaff16);
beliefPenaltyMPC_LA_VSUB3_17(beliefPenaltyMPC_llbbyslb16, beliefPenaltyMPC_dslbaff16, beliefPenaltyMPC_llb16, beliefPenaltyMPC_dllbaff16);
beliefPenaltyMPC_LA_VSUB2_INDEXED_7(beliefPenaltyMPC_riub16, beliefPenaltyMPC_dzaff16, beliefPenaltyMPC_ubIdx16, beliefPenaltyMPC_dsubaff16);
beliefPenaltyMPC_LA_VSUB3_7(beliefPenaltyMPC_lubbysub16, beliefPenaltyMPC_dsubaff16, beliefPenaltyMPC_lub16, beliefPenaltyMPC_dlubaff16);
beliefPenaltyMPC_LA_VSUB_INDEXED_17(beliefPenaltyMPC_dzaff17, beliefPenaltyMPC_lbIdx17, beliefPenaltyMPC_rilb17, beliefPenaltyMPC_dslbaff17);
beliefPenaltyMPC_LA_VSUB3_17(beliefPenaltyMPC_llbbyslb17, beliefPenaltyMPC_dslbaff17, beliefPenaltyMPC_llb17, beliefPenaltyMPC_dllbaff17);
beliefPenaltyMPC_LA_VSUB2_INDEXED_7(beliefPenaltyMPC_riub17, beliefPenaltyMPC_dzaff17, beliefPenaltyMPC_ubIdx17, beliefPenaltyMPC_dsubaff17);
beliefPenaltyMPC_LA_VSUB3_7(beliefPenaltyMPC_lubbysub17, beliefPenaltyMPC_dsubaff17, beliefPenaltyMPC_lub17, beliefPenaltyMPC_dlubaff17);
beliefPenaltyMPC_LA_VSUB_INDEXED_17(beliefPenaltyMPC_dzaff18, beliefPenaltyMPC_lbIdx18, beliefPenaltyMPC_rilb18, beliefPenaltyMPC_dslbaff18);
beliefPenaltyMPC_LA_VSUB3_17(beliefPenaltyMPC_llbbyslb18, beliefPenaltyMPC_dslbaff18, beliefPenaltyMPC_llb18, beliefPenaltyMPC_dllbaff18);
beliefPenaltyMPC_LA_VSUB2_INDEXED_7(beliefPenaltyMPC_riub18, beliefPenaltyMPC_dzaff18, beliefPenaltyMPC_ubIdx18, beliefPenaltyMPC_dsubaff18);
beliefPenaltyMPC_LA_VSUB3_7(beliefPenaltyMPC_lubbysub18, beliefPenaltyMPC_dsubaff18, beliefPenaltyMPC_lub18, beliefPenaltyMPC_dlubaff18);
beliefPenaltyMPC_LA_VSUB_INDEXED_17(beliefPenaltyMPC_dzaff19, beliefPenaltyMPC_lbIdx19, beliefPenaltyMPC_rilb19, beliefPenaltyMPC_dslbaff19);
beliefPenaltyMPC_LA_VSUB3_17(beliefPenaltyMPC_llbbyslb19, beliefPenaltyMPC_dslbaff19, beliefPenaltyMPC_llb19, beliefPenaltyMPC_dllbaff19);
beliefPenaltyMPC_LA_VSUB2_INDEXED_7(beliefPenaltyMPC_riub19, beliefPenaltyMPC_dzaff19, beliefPenaltyMPC_ubIdx19, beliefPenaltyMPC_dsubaff19);
beliefPenaltyMPC_LA_VSUB3_7(beliefPenaltyMPC_lubbysub19, beliefPenaltyMPC_dsubaff19, beliefPenaltyMPC_lub19, beliefPenaltyMPC_dlubaff19);
beliefPenaltyMPC_LA_VSUB_INDEXED_17(beliefPenaltyMPC_dzaff20, beliefPenaltyMPC_lbIdx20, beliefPenaltyMPC_rilb20, beliefPenaltyMPC_dslbaff20);
beliefPenaltyMPC_LA_VSUB3_17(beliefPenaltyMPC_llbbyslb20, beliefPenaltyMPC_dslbaff20, beliefPenaltyMPC_llb20, beliefPenaltyMPC_dllbaff20);
beliefPenaltyMPC_LA_VSUB2_INDEXED_7(beliefPenaltyMPC_riub20, beliefPenaltyMPC_dzaff20, beliefPenaltyMPC_ubIdx20, beliefPenaltyMPC_dsubaff20);
beliefPenaltyMPC_LA_VSUB3_7(beliefPenaltyMPC_lubbysub20, beliefPenaltyMPC_dsubaff20, beliefPenaltyMPC_lub20, beliefPenaltyMPC_dlubaff20);
beliefPenaltyMPC_LA_VSUB_INDEXED_17(beliefPenaltyMPC_dzaff21, beliefPenaltyMPC_lbIdx21, beliefPenaltyMPC_rilb21, beliefPenaltyMPC_dslbaff21);
beliefPenaltyMPC_LA_VSUB3_17(beliefPenaltyMPC_llbbyslb21, beliefPenaltyMPC_dslbaff21, beliefPenaltyMPC_llb21, beliefPenaltyMPC_dllbaff21);
beliefPenaltyMPC_LA_VSUB2_INDEXED_7(beliefPenaltyMPC_riub21, beliefPenaltyMPC_dzaff21, beliefPenaltyMPC_ubIdx21, beliefPenaltyMPC_dsubaff21);
beliefPenaltyMPC_LA_VSUB3_7(beliefPenaltyMPC_lubbysub21, beliefPenaltyMPC_dsubaff21, beliefPenaltyMPC_lub21, beliefPenaltyMPC_dlubaff21);
beliefPenaltyMPC_LA_VSUB_INDEXED_17(beliefPenaltyMPC_dzaff22, beliefPenaltyMPC_lbIdx22, beliefPenaltyMPC_rilb22, beliefPenaltyMPC_dslbaff22);
beliefPenaltyMPC_LA_VSUB3_17(beliefPenaltyMPC_llbbyslb22, beliefPenaltyMPC_dslbaff22, beliefPenaltyMPC_llb22, beliefPenaltyMPC_dllbaff22);
beliefPenaltyMPC_LA_VSUB2_INDEXED_7(beliefPenaltyMPC_riub22, beliefPenaltyMPC_dzaff22, beliefPenaltyMPC_ubIdx22, beliefPenaltyMPC_dsubaff22);
beliefPenaltyMPC_LA_VSUB3_7(beliefPenaltyMPC_lubbysub22, beliefPenaltyMPC_dsubaff22, beliefPenaltyMPC_lub22, beliefPenaltyMPC_dlubaff22);
beliefPenaltyMPC_LA_VSUB_INDEXED_17(beliefPenaltyMPC_dzaff23, beliefPenaltyMPC_lbIdx23, beliefPenaltyMPC_rilb23, beliefPenaltyMPC_dslbaff23);
beliefPenaltyMPC_LA_VSUB3_17(beliefPenaltyMPC_llbbyslb23, beliefPenaltyMPC_dslbaff23, beliefPenaltyMPC_llb23, beliefPenaltyMPC_dllbaff23);
beliefPenaltyMPC_LA_VSUB2_INDEXED_7(beliefPenaltyMPC_riub23, beliefPenaltyMPC_dzaff23, beliefPenaltyMPC_ubIdx23, beliefPenaltyMPC_dsubaff23);
beliefPenaltyMPC_LA_VSUB3_7(beliefPenaltyMPC_lubbysub23, beliefPenaltyMPC_dsubaff23, beliefPenaltyMPC_lub23, beliefPenaltyMPC_dlubaff23);
beliefPenaltyMPC_LA_VSUB_INDEXED_5(beliefPenaltyMPC_dzaff24, beliefPenaltyMPC_lbIdx24, beliefPenaltyMPC_rilb24, beliefPenaltyMPC_dslbaff24);
beliefPenaltyMPC_LA_VSUB3_5(beliefPenaltyMPC_llbbyslb24, beliefPenaltyMPC_dslbaff24, beliefPenaltyMPC_llb24, beliefPenaltyMPC_dllbaff24);
beliefPenaltyMPC_LA_VSUB2_INDEXED_5(beliefPenaltyMPC_riub24, beliefPenaltyMPC_dzaff24, beliefPenaltyMPC_ubIdx24, beliefPenaltyMPC_dsubaff24);
beliefPenaltyMPC_LA_VSUB3_5(beliefPenaltyMPC_lubbysub24, beliefPenaltyMPC_dsubaff24, beliefPenaltyMPC_lub24, beliefPenaltyMPC_dlubaff24);
info->lsit_aff = beliefPenaltyMPC_LINESEARCH_BACKTRACKING_AFFINE(beliefPenaltyMPC_l, beliefPenaltyMPC_s, beliefPenaltyMPC_dl_aff, beliefPenaltyMPC_ds_aff, &info->step_aff, &info->mu_aff);
if( info->lsit_aff == beliefPenaltyMPC_NOPROGRESS ){
exitcode = beliefPenaltyMPC_NOPROGRESS; break;
}
sigma_3rdroot = info->mu_aff / info->mu;
info->sigma = sigma_3rdroot*sigma_3rdroot*sigma_3rdroot;
musigma = info->mu * info->sigma;
beliefPenaltyMPC_LA_VSUB5_586(beliefPenaltyMPC_ds_aff, beliefPenaltyMPC_dl_aff, musigma, beliefPenaltyMPC_ccrhs);
beliefPenaltyMPC_LA_VSUB6_INDEXED_17_7_17(beliefPenaltyMPC_ccrhsub00, beliefPenaltyMPC_sub00, beliefPenaltyMPC_ubIdx00, beliefPenaltyMPC_ccrhsl00, beliefPenaltyMPC_slb00, beliefPenaltyMPC_lbIdx00, beliefPenaltyMPC_rd00);
beliefPenaltyMPC_LA_VSUB6_INDEXED_17_7_17(beliefPenaltyMPC_ccrhsub01, beliefPenaltyMPC_sub01, beliefPenaltyMPC_ubIdx01, beliefPenaltyMPC_ccrhsl01, beliefPenaltyMPC_slb01, beliefPenaltyMPC_lbIdx01, beliefPenaltyMPC_rd01);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_17(beliefPenaltyMPC_Phi00, beliefPenaltyMPC_rd00, beliefPenaltyMPC_Lbyrd00);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_17(beliefPenaltyMPC_Phi01, beliefPenaltyMPC_rd01, beliefPenaltyMPC_Lbyrd01);
beliefPenaltyMPC_LA_DENSE_2MVMADD_10_17_17(beliefPenaltyMPC_V00, beliefPenaltyMPC_Lbyrd00, beliefPenaltyMPC_W01, beliefPenaltyMPC_Lbyrd01, beliefPenaltyMPC_beta00);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_10(beliefPenaltyMPC_Ld00, beliefPenaltyMPC_beta00, beliefPenaltyMPC_yy00);
beliefPenaltyMPC_LA_VSUB6_INDEXED_17_7_17(beliefPenaltyMPC_ccrhsub02, beliefPenaltyMPC_sub02, beliefPenaltyMPC_ubIdx02, beliefPenaltyMPC_ccrhsl02, beliefPenaltyMPC_slb02, beliefPenaltyMPC_lbIdx02, beliefPenaltyMPC_rd02);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_17(beliefPenaltyMPC_Phi02, beliefPenaltyMPC_rd02, beliefPenaltyMPC_Lbyrd02);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_5_17_17(beliefPenaltyMPC_V01, beliefPenaltyMPC_Lbyrd01, beliefPenaltyMPC_W02, beliefPenaltyMPC_Lbyrd02, beliefPenaltyMPC_beta01);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_5_10(beliefPenaltyMPC_Lsd01, beliefPenaltyMPC_yy00, beliefPenaltyMPC_beta01, beliefPenaltyMPC_bmy01);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_5(beliefPenaltyMPC_Ld01, beliefPenaltyMPC_bmy01, beliefPenaltyMPC_yy01);
beliefPenaltyMPC_LA_VSUB6_INDEXED_17_7_17(beliefPenaltyMPC_ccrhsub03, beliefPenaltyMPC_sub03, beliefPenaltyMPC_ubIdx03, beliefPenaltyMPC_ccrhsl03, beliefPenaltyMPC_slb03, beliefPenaltyMPC_lbIdx03, beliefPenaltyMPC_rd03);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_17(beliefPenaltyMPC_Phi03, beliefPenaltyMPC_rd03, beliefPenaltyMPC_Lbyrd03);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_5_17_17(beliefPenaltyMPC_V02, beliefPenaltyMPC_Lbyrd02, beliefPenaltyMPC_W03, beliefPenaltyMPC_Lbyrd03, beliefPenaltyMPC_beta02);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_5_5(beliefPenaltyMPC_Lsd02, beliefPenaltyMPC_yy01, beliefPenaltyMPC_beta02, beliefPenaltyMPC_bmy02);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_5(beliefPenaltyMPC_Ld02, beliefPenaltyMPC_bmy02, beliefPenaltyMPC_yy02);
beliefPenaltyMPC_LA_VSUB6_INDEXED_17_7_17(beliefPenaltyMPC_ccrhsub04, beliefPenaltyMPC_sub04, beliefPenaltyMPC_ubIdx04, beliefPenaltyMPC_ccrhsl04, beliefPenaltyMPC_slb04, beliefPenaltyMPC_lbIdx04, beliefPenaltyMPC_rd04);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_17(beliefPenaltyMPC_Phi04, beliefPenaltyMPC_rd04, beliefPenaltyMPC_Lbyrd04);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_5_17_17(beliefPenaltyMPC_V03, beliefPenaltyMPC_Lbyrd03, beliefPenaltyMPC_W04, beliefPenaltyMPC_Lbyrd04, beliefPenaltyMPC_beta03);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_5_5(beliefPenaltyMPC_Lsd03, beliefPenaltyMPC_yy02, beliefPenaltyMPC_beta03, beliefPenaltyMPC_bmy03);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_5(beliefPenaltyMPC_Ld03, beliefPenaltyMPC_bmy03, beliefPenaltyMPC_yy03);
beliefPenaltyMPC_LA_VSUB6_INDEXED_17_7_17(beliefPenaltyMPC_ccrhsub05, beliefPenaltyMPC_sub05, beliefPenaltyMPC_ubIdx05, beliefPenaltyMPC_ccrhsl05, beliefPenaltyMPC_slb05, beliefPenaltyMPC_lbIdx05, beliefPenaltyMPC_rd05);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_17(beliefPenaltyMPC_Phi05, beliefPenaltyMPC_rd05, beliefPenaltyMPC_Lbyrd05);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_5_17_17(beliefPenaltyMPC_V04, beliefPenaltyMPC_Lbyrd04, beliefPenaltyMPC_W05, beliefPenaltyMPC_Lbyrd05, beliefPenaltyMPC_beta04);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_5_5(beliefPenaltyMPC_Lsd04, beliefPenaltyMPC_yy03, beliefPenaltyMPC_beta04, beliefPenaltyMPC_bmy04);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_5(beliefPenaltyMPC_Ld04, beliefPenaltyMPC_bmy04, beliefPenaltyMPC_yy04);
beliefPenaltyMPC_LA_VSUB6_INDEXED_17_7_17(beliefPenaltyMPC_ccrhsub06, beliefPenaltyMPC_sub06, beliefPenaltyMPC_ubIdx06, beliefPenaltyMPC_ccrhsl06, beliefPenaltyMPC_slb06, beliefPenaltyMPC_lbIdx06, beliefPenaltyMPC_rd06);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_17(beliefPenaltyMPC_Phi06, beliefPenaltyMPC_rd06, beliefPenaltyMPC_Lbyrd06);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_5_17_17(beliefPenaltyMPC_V05, beliefPenaltyMPC_Lbyrd05, beliefPenaltyMPC_W06, beliefPenaltyMPC_Lbyrd06, beliefPenaltyMPC_beta05);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_5_5(beliefPenaltyMPC_Lsd05, beliefPenaltyMPC_yy04, beliefPenaltyMPC_beta05, beliefPenaltyMPC_bmy05);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_5(beliefPenaltyMPC_Ld05, beliefPenaltyMPC_bmy05, beliefPenaltyMPC_yy05);
beliefPenaltyMPC_LA_VSUB6_INDEXED_17_7_17(beliefPenaltyMPC_ccrhsub07, beliefPenaltyMPC_sub07, beliefPenaltyMPC_ubIdx07, beliefPenaltyMPC_ccrhsl07, beliefPenaltyMPC_slb07, beliefPenaltyMPC_lbIdx07, beliefPenaltyMPC_rd07);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_17(beliefPenaltyMPC_Phi07, beliefPenaltyMPC_rd07, beliefPenaltyMPC_Lbyrd07);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_5_17_17(beliefPenaltyMPC_V06, beliefPenaltyMPC_Lbyrd06, beliefPenaltyMPC_W07, beliefPenaltyMPC_Lbyrd07, beliefPenaltyMPC_beta06);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_5_5(beliefPenaltyMPC_Lsd06, beliefPenaltyMPC_yy05, beliefPenaltyMPC_beta06, beliefPenaltyMPC_bmy06);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_5(beliefPenaltyMPC_Ld06, beliefPenaltyMPC_bmy06, beliefPenaltyMPC_yy06);
beliefPenaltyMPC_LA_VSUB6_INDEXED_17_7_17(beliefPenaltyMPC_ccrhsub08, beliefPenaltyMPC_sub08, beliefPenaltyMPC_ubIdx08, beliefPenaltyMPC_ccrhsl08, beliefPenaltyMPC_slb08, beliefPenaltyMPC_lbIdx08, beliefPenaltyMPC_rd08);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_17(beliefPenaltyMPC_Phi08, beliefPenaltyMPC_rd08, beliefPenaltyMPC_Lbyrd08);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_5_17_17(beliefPenaltyMPC_V07, beliefPenaltyMPC_Lbyrd07, beliefPenaltyMPC_W08, beliefPenaltyMPC_Lbyrd08, beliefPenaltyMPC_beta07);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_5_5(beliefPenaltyMPC_Lsd07, beliefPenaltyMPC_yy06, beliefPenaltyMPC_beta07, beliefPenaltyMPC_bmy07);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_5(beliefPenaltyMPC_Ld07, beliefPenaltyMPC_bmy07, beliefPenaltyMPC_yy07);
beliefPenaltyMPC_LA_VSUB6_INDEXED_17_7_17(beliefPenaltyMPC_ccrhsub09, beliefPenaltyMPC_sub09, beliefPenaltyMPC_ubIdx09, beliefPenaltyMPC_ccrhsl09, beliefPenaltyMPC_slb09, beliefPenaltyMPC_lbIdx09, beliefPenaltyMPC_rd09);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_17(beliefPenaltyMPC_Phi09, beliefPenaltyMPC_rd09, beliefPenaltyMPC_Lbyrd09);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_5_17_17(beliefPenaltyMPC_V08, beliefPenaltyMPC_Lbyrd08, beliefPenaltyMPC_W09, beliefPenaltyMPC_Lbyrd09, beliefPenaltyMPC_beta08);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_5_5(beliefPenaltyMPC_Lsd08, beliefPenaltyMPC_yy07, beliefPenaltyMPC_beta08, beliefPenaltyMPC_bmy08);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_5(beliefPenaltyMPC_Ld08, beliefPenaltyMPC_bmy08, beliefPenaltyMPC_yy08);
beliefPenaltyMPC_LA_VSUB6_INDEXED_17_7_17(beliefPenaltyMPC_ccrhsub10, beliefPenaltyMPC_sub10, beliefPenaltyMPC_ubIdx10, beliefPenaltyMPC_ccrhsl10, beliefPenaltyMPC_slb10, beliefPenaltyMPC_lbIdx10, beliefPenaltyMPC_rd10);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_17(beliefPenaltyMPC_Phi10, beliefPenaltyMPC_rd10, beliefPenaltyMPC_Lbyrd10);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_5_17_17(beliefPenaltyMPC_V09, beliefPenaltyMPC_Lbyrd09, beliefPenaltyMPC_W10, beliefPenaltyMPC_Lbyrd10, beliefPenaltyMPC_beta09);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_5_5(beliefPenaltyMPC_Lsd09, beliefPenaltyMPC_yy08, beliefPenaltyMPC_beta09, beliefPenaltyMPC_bmy09);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_5(beliefPenaltyMPC_Ld09, beliefPenaltyMPC_bmy09, beliefPenaltyMPC_yy09);
beliefPenaltyMPC_LA_VSUB6_INDEXED_17_7_17(beliefPenaltyMPC_ccrhsub11, beliefPenaltyMPC_sub11, beliefPenaltyMPC_ubIdx11, beliefPenaltyMPC_ccrhsl11, beliefPenaltyMPC_slb11, beliefPenaltyMPC_lbIdx11, beliefPenaltyMPC_rd11);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_17(beliefPenaltyMPC_Phi11, beliefPenaltyMPC_rd11, beliefPenaltyMPC_Lbyrd11);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_5_17_17(beliefPenaltyMPC_V10, beliefPenaltyMPC_Lbyrd10, beliefPenaltyMPC_W11, beliefPenaltyMPC_Lbyrd11, beliefPenaltyMPC_beta10);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_5_5(beliefPenaltyMPC_Lsd10, beliefPenaltyMPC_yy09, beliefPenaltyMPC_beta10, beliefPenaltyMPC_bmy10);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_5(beliefPenaltyMPC_Ld10, beliefPenaltyMPC_bmy10, beliefPenaltyMPC_yy10);
beliefPenaltyMPC_LA_VSUB6_INDEXED_17_7_17(beliefPenaltyMPC_ccrhsub12, beliefPenaltyMPC_sub12, beliefPenaltyMPC_ubIdx12, beliefPenaltyMPC_ccrhsl12, beliefPenaltyMPC_slb12, beliefPenaltyMPC_lbIdx12, beliefPenaltyMPC_rd12);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_17(beliefPenaltyMPC_Phi12, beliefPenaltyMPC_rd12, beliefPenaltyMPC_Lbyrd12);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_5_17_17(beliefPenaltyMPC_V11, beliefPenaltyMPC_Lbyrd11, beliefPenaltyMPC_W12, beliefPenaltyMPC_Lbyrd12, beliefPenaltyMPC_beta11);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_5_5(beliefPenaltyMPC_Lsd11, beliefPenaltyMPC_yy10, beliefPenaltyMPC_beta11, beliefPenaltyMPC_bmy11);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_5(beliefPenaltyMPC_Ld11, beliefPenaltyMPC_bmy11, beliefPenaltyMPC_yy11);
beliefPenaltyMPC_LA_VSUB6_INDEXED_17_7_17(beliefPenaltyMPC_ccrhsub13, beliefPenaltyMPC_sub13, beliefPenaltyMPC_ubIdx13, beliefPenaltyMPC_ccrhsl13, beliefPenaltyMPC_slb13, beliefPenaltyMPC_lbIdx13, beliefPenaltyMPC_rd13);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_17(beliefPenaltyMPC_Phi13, beliefPenaltyMPC_rd13, beliefPenaltyMPC_Lbyrd13);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_5_17_17(beliefPenaltyMPC_V12, beliefPenaltyMPC_Lbyrd12, beliefPenaltyMPC_W13, beliefPenaltyMPC_Lbyrd13, beliefPenaltyMPC_beta12);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_5_5(beliefPenaltyMPC_Lsd12, beliefPenaltyMPC_yy11, beliefPenaltyMPC_beta12, beliefPenaltyMPC_bmy12);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_5(beliefPenaltyMPC_Ld12, beliefPenaltyMPC_bmy12, beliefPenaltyMPC_yy12);
beliefPenaltyMPC_LA_VSUB6_INDEXED_17_7_17(beliefPenaltyMPC_ccrhsub14, beliefPenaltyMPC_sub14, beliefPenaltyMPC_ubIdx14, beliefPenaltyMPC_ccrhsl14, beliefPenaltyMPC_slb14, beliefPenaltyMPC_lbIdx14, beliefPenaltyMPC_rd14);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_17(beliefPenaltyMPC_Phi14, beliefPenaltyMPC_rd14, beliefPenaltyMPC_Lbyrd14);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_5_17_17(beliefPenaltyMPC_V13, beliefPenaltyMPC_Lbyrd13, beliefPenaltyMPC_W14, beliefPenaltyMPC_Lbyrd14, beliefPenaltyMPC_beta13);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_5_5(beliefPenaltyMPC_Lsd13, beliefPenaltyMPC_yy12, beliefPenaltyMPC_beta13, beliefPenaltyMPC_bmy13);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_5(beliefPenaltyMPC_Ld13, beliefPenaltyMPC_bmy13, beliefPenaltyMPC_yy13);
beliefPenaltyMPC_LA_VSUB6_INDEXED_17_7_17(beliefPenaltyMPC_ccrhsub15, beliefPenaltyMPC_sub15, beliefPenaltyMPC_ubIdx15, beliefPenaltyMPC_ccrhsl15, beliefPenaltyMPC_slb15, beliefPenaltyMPC_lbIdx15, beliefPenaltyMPC_rd15);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_17(beliefPenaltyMPC_Phi15, beliefPenaltyMPC_rd15, beliefPenaltyMPC_Lbyrd15);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_5_17_17(beliefPenaltyMPC_V14, beliefPenaltyMPC_Lbyrd14, beliefPenaltyMPC_W15, beliefPenaltyMPC_Lbyrd15, beliefPenaltyMPC_beta14);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_5_5(beliefPenaltyMPC_Lsd14, beliefPenaltyMPC_yy13, beliefPenaltyMPC_beta14, beliefPenaltyMPC_bmy14);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_5(beliefPenaltyMPC_Ld14, beliefPenaltyMPC_bmy14, beliefPenaltyMPC_yy14);
beliefPenaltyMPC_LA_VSUB6_INDEXED_17_7_17(beliefPenaltyMPC_ccrhsub16, beliefPenaltyMPC_sub16, beliefPenaltyMPC_ubIdx16, beliefPenaltyMPC_ccrhsl16, beliefPenaltyMPC_slb16, beliefPenaltyMPC_lbIdx16, beliefPenaltyMPC_rd16);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_17(beliefPenaltyMPC_Phi16, beliefPenaltyMPC_rd16, beliefPenaltyMPC_Lbyrd16);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_5_17_17(beliefPenaltyMPC_V15, beliefPenaltyMPC_Lbyrd15, beliefPenaltyMPC_W16, beliefPenaltyMPC_Lbyrd16, beliefPenaltyMPC_beta15);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_5_5(beliefPenaltyMPC_Lsd15, beliefPenaltyMPC_yy14, beliefPenaltyMPC_beta15, beliefPenaltyMPC_bmy15);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_5(beliefPenaltyMPC_Ld15, beliefPenaltyMPC_bmy15, beliefPenaltyMPC_yy15);
beliefPenaltyMPC_LA_VSUB6_INDEXED_17_7_17(beliefPenaltyMPC_ccrhsub17, beliefPenaltyMPC_sub17, beliefPenaltyMPC_ubIdx17, beliefPenaltyMPC_ccrhsl17, beliefPenaltyMPC_slb17, beliefPenaltyMPC_lbIdx17, beliefPenaltyMPC_rd17);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_17(beliefPenaltyMPC_Phi17, beliefPenaltyMPC_rd17, beliefPenaltyMPC_Lbyrd17);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_5_17_17(beliefPenaltyMPC_V16, beliefPenaltyMPC_Lbyrd16, beliefPenaltyMPC_W17, beliefPenaltyMPC_Lbyrd17, beliefPenaltyMPC_beta16);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_5_5(beliefPenaltyMPC_Lsd16, beliefPenaltyMPC_yy15, beliefPenaltyMPC_beta16, beliefPenaltyMPC_bmy16);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_5(beliefPenaltyMPC_Ld16, beliefPenaltyMPC_bmy16, beliefPenaltyMPC_yy16);
beliefPenaltyMPC_LA_VSUB6_INDEXED_17_7_17(beliefPenaltyMPC_ccrhsub18, beliefPenaltyMPC_sub18, beliefPenaltyMPC_ubIdx18, beliefPenaltyMPC_ccrhsl18, beliefPenaltyMPC_slb18, beliefPenaltyMPC_lbIdx18, beliefPenaltyMPC_rd18);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_17(beliefPenaltyMPC_Phi18, beliefPenaltyMPC_rd18, beliefPenaltyMPC_Lbyrd18);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_5_17_17(beliefPenaltyMPC_V17, beliefPenaltyMPC_Lbyrd17, beliefPenaltyMPC_W18, beliefPenaltyMPC_Lbyrd18, beliefPenaltyMPC_beta17);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_5_5(beliefPenaltyMPC_Lsd17, beliefPenaltyMPC_yy16, beliefPenaltyMPC_beta17, beliefPenaltyMPC_bmy17);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_5(beliefPenaltyMPC_Ld17, beliefPenaltyMPC_bmy17, beliefPenaltyMPC_yy17);
beliefPenaltyMPC_LA_VSUB6_INDEXED_17_7_17(beliefPenaltyMPC_ccrhsub19, beliefPenaltyMPC_sub19, beliefPenaltyMPC_ubIdx19, beliefPenaltyMPC_ccrhsl19, beliefPenaltyMPC_slb19, beliefPenaltyMPC_lbIdx19, beliefPenaltyMPC_rd19);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_17(beliefPenaltyMPC_Phi19, beliefPenaltyMPC_rd19, beliefPenaltyMPC_Lbyrd19);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_5_17_17(beliefPenaltyMPC_V18, beliefPenaltyMPC_Lbyrd18, beliefPenaltyMPC_W19, beliefPenaltyMPC_Lbyrd19, beliefPenaltyMPC_beta18);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_5_5(beliefPenaltyMPC_Lsd18, beliefPenaltyMPC_yy17, beliefPenaltyMPC_beta18, beliefPenaltyMPC_bmy18);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_5(beliefPenaltyMPC_Ld18, beliefPenaltyMPC_bmy18, beliefPenaltyMPC_yy18);
beliefPenaltyMPC_LA_VSUB6_INDEXED_17_7_17(beliefPenaltyMPC_ccrhsub20, beliefPenaltyMPC_sub20, beliefPenaltyMPC_ubIdx20, beliefPenaltyMPC_ccrhsl20, beliefPenaltyMPC_slb20, beliefPenaltyMPC_lbIdx20, beliefPenaltyMPC_rd20);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_17(beliefPenaltyMPC_Phi20, beliefPenaltyMPC_rd20, beliefPenaltyMPC_Lbyrd20);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_5_17_17(beliefPenaltyMPC_V19, beliefPenaltyMPC_Lbyrd19, beliefPenaltyMPC_W20, beliefPenaltyMPC_Lbyrd20, beliefPenaltyMPC_beta19);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_5_5(beliefPenaltyMPC_Lsd19, beliefPenaltyMPC_yy18, beliefPenaltyMPC_beta19, beliefPenaltyMPC_bmy19);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_5(beliefPenaltyMPC_Ld19, beliefPenaltyMPC_bmy19, beliefPenaltyMPC_yy19);
beliefPenaltyMPC_LA_VSUB6_INDEXED_17_7_17(beliefPenaltyMPC_ccrhsub21, beliefPenaltyMPC_sub21, beliefPenaltyMPC_ubIdx21, beliefPenaltyMPC_ccrhsl21, beliefPenaltyMPC_slb21, beliefPenaltyMPC_lbIdx21, beliefPenaltyMPC_rd21);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_17(beliefPenaltyMPC_Phi21, beliefPenaltyMPC_rd21, beliefPenaltyMPC_Lbyrd21);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_5_17_17(beliefPenaltyMPC_V20, beliefPenaltyMPC_Lbyrd20, beliefPenaltyMPC_W21, beliefPenaltyMPC_Lbyrd21, beliefPenaltyMPC_beta20);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_5_5(beliefPenaltyMPC_Lsd20, beliefPenaltyMPC_yy19, beliefPenaltyMPC_beta20, beliefPenaltyMPC_bmy20);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_5(beliefPenaltyMPC_Ld20, beliefPenaltyMPC_bmy20, beliefPenaltyMPC_yy20);
beliefPenaltyMPC_LA_VSUB6_INDEXED_17_7_17(beliefPenaltyMPC_ccrhsub22, beliefPenaltyMPC_sub22, beliefPenaltyMPC_ubIdx22, beliefPenaltyMPC_ccrhsl22, beliefPenaltyMPC_slb22, beliefPenaltyMPC_lbIdx22, beliefPenaltyMPC_rd22);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_17(beliefPenaltyMPC_Phi22, beliefPenaltyMPC_rd22, beliefPenaltyMPC_Lbyrd22);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_5_17_17(beliefPenaltyMPC_V21, beliefPenaltyMPC_Lbyrd21, beliefPenaltyMPC_W22, beliefPenaltyMPC_Lbyrd22, beliefPenaltyMPC_beta21);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_5_5(beliefPenaltyMPC_Lsd21, beliefPenaltyMPC_yy20, beliefPenaltyMPC_beta21, beliefPenaltyMPC_bmy21);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_5(beliefPenaltyMPC_Ld21, beliefPenaltyMPC_bmy21, beliefPenaltyMPC_yy21);
beliefPenaltyMPC_LA_VSUB6_INDEXED_17_7_17(beliefPenaltyMPC_ccrhsub23, beliefPenaltyMPC_sub23, beliefPenaltyMPC_ubIdx23, beliefPenaltyMPC_ccrhsl23, beliefPenaltyMPC_slb23, beliefPenaltyMPC_lbIdx23, beliefPenaltyMPC_rd23);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_17(beliefPenaltyMPC_Phi23, beliefPenaltyMPC_rd23, beliefPenaltyMPC_Lbyrd23);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_5_17_17(beliefPenaltyMPC_V22, beliefPenaltyMPC_Lbyrd22, beliefPenaltyMPC_W23, beliefPenaltyMPC_Lbyrd23, beliefPenaltyMPC_beta22);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_5_5(beliefPenaltyMPC_Lsd22, beliefPenaltyMPC_yy21, beliefPenaltyMPC_beta22, beliefPenaltyMPC_bmy22);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_5(beliefPenaltyMPC_Ld22, beliefPenaltyMPC_bmy22, beliefPenaltyMPC_yy22);
beliefPenaltyMPC_LA_VSUB6_INDEXED_5_5_5(beliefPenaltyMPC_ccrhsub24, beliefPenaltyMPC_sub24, beliefPenaltyMPC_ubIdx24, beliefPenaltyMPC_ccrhsl24, beliefPenaltyMPC_slb24, beliefPenaltyMPC_lbIdx24, beliefPenaltyMPC_rd24);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_5(beliefPenaltyMPC_Phi24, beliefPenaltyMPC_rd24, beliefPenaltyMPC_Lbyrd24);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_5_17_5(beliefPenaltyMPC_V23, beliefPenaltyMPC_Lbyrd23, beliefPenaltyMPC_W24, beliefPenaltyMPC_Lbyrd24, beliefPenaltyMPC_beta23);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_5_5(beliefPenaltyMPC_Lsd23, beliefPenaltyMPC_yy22, beliefPenaltyMPC_beta23, beliefPenaltyMPC_bmy23);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_5(beliefPenaltyMPC_Ld23, beliefPenaltyMPC_bmy23, beliefPenaltyMPC_yy23);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_5(beliefPenaltyMPC_Ld23, beliefPenaltyMPC_yy23, beliefPenaltyMPC_dvcc23);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_5_5(beliefPenaltyMPC_Lsd23, beliefPenaltyMPC_dvcc23, beliefPenaltyMPC_yy22, beliefPenaltyMPC_bmy22);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_5(beliefPenaltyMPC_Ld22, beliefPenaltyMPC_bmy22, beliefPenaltyMPC_dvcc22);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_5_5(beliefPenaltyMPC_Lsd22, beliefPenaltyMPC_dvcc22, beliefPenaltyMPC_yy21, beliefPenaltyMPC_bmy21);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_5(beliefPenaltyMPC_Ld21, beliefPenaltyMPC_bmy21, beliefPenaltyMPC_dvcc21);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_5_5(beliefPenaltyMPC_Lsd21, beliefPenaltyMPC_dvcc21, beliefPenaltyMPC_yy20, beliefPenaltyMPC_bmy20);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_5(beliefPenaltyMPC_Ld20, beliefPenaltyMPC_bmy20, beliefPenaltyMPC_dvcc20);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_5_5(beliefPenaltyMPC_Lsd20, beliefPenaltyMPC_dvcc20, beliefPenaltyMPC_yy19, beliefPenaltyMPC_bmy19);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_5(beliefPenaltyMPC_Ld19, beliefPenaltyMPC_bmy19, beliefPenaltyMPC_dvcc19);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_5_5(beliefPenaltyMPC_Lsd19, beliefPenaltyMPC_dvcc19, beliefPenaltyMPC_yy18, beliefPenaltyMPC_bmy18);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_5(beliefPenaltyMPC_Ld18, beliefPenaltyMPC_bmy18, beliefPenaltyMPC_dvcc18);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_5_5(beliefPenaltyMPC_Lsd18, beliefPenaltyMPC_dvcc18, beliefPenaltyMPC_yy17, beliefPenaltyMPC_bmy17);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_5(beliefPenaltyMPC_Ld17, beliefPenaltyMPC_bmy17, beliefPenaltyMPC_dvcc17);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_5_5(beliefPenaltyMPC_Lsd17, beliefPenaltyMPC_dvcc17, beliefPenaltyMPC_yy16, beliefPenaltyMPC_bmy16);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_5(beliefPenaltyMPC_Ld16, beliefPenaltyMPC_bmy16, beliefPenaltyMPC_dvcc16);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_5_5(beliefPenaltyMPC_Lsd16, beliefPenaltyMPC_dvcc16, beliefPenaltyMPC_yy15, beliefPenaltyMPC_bmy15);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_5(beliefPenaltyMPC_Ld15, beliefPenaltyMPC_bmy15, beliefPenaltyMPC_dvcc15);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_5_5(beliefPenaltyMPC_Lsd15, beliefPenaltyMPC_dvcc15, beliefPenaltyMPC_yy14, beliefPenaltyMPC_bmy14);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_5(beliefPenaltyMPC_Ld14, beliefPenaltyMPC_bmy14, beliefPenaltyMPC_dvcc14);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_5_5(beliefPenaltyMPC_Lsd14, beliefPenaltyMPC_dvcc14, beliefPenaltyMPC_yy13, beliefPenaltyMPC_bmy13);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_5(beliefPenaltyMPC_Ld13, beliefPenaltyMPC_bmy13, beliefPenaltyMPC_dvcc13);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_5_5(beliefPenaltyMPC_Lsd13, beliefPenaltyMPC_dvcc13, beliefPenaltyMPC_yy12, beliefPenaltyMPC_bmy12);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_5(beliefPenaltyMPC_Ld12, beliefPenaltyMPC_bmy12, beliefPenaltyMPC_dvcc12);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_5_5(beliefPenaltyMPC_Lsd12, beliefPenaltyMPC_dvcc12, beliefPenaltyMPC_yy11, beliefPenaltyMPC_bmy11);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_5(beliefPenaltyMPC_Ld11, beliefPenaltyMPC_bmy11, beliefPenaltyMPC_dvcc11);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_5_5(beliefPenaltyMPC_Lsd11, beliefPenaltyMPC_dvcc11, beliefPenaltyMPC_yy10, beliefPenaltyMPC_bmy10);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_5(beliefPenaltyMPC_Ld10, beliefPenaltyMPC_bmy10, beliefPenaltyMPC_dvcc10);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_5_5(beliefPenaltyMPC_Lsd10, beliefPenaltyMPC_dvcc10, beliefPenaltyMPC_yy09, beliefPenaltyMPC_bmy09);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_5(beliefPenaltyMPC_Ld09, beliefPenaltyMPC_bmy09, beliefPenaltyMPC_dvcc09);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_5_5(beliefPenaltyMPC_Lsd09, beliefPenaltyMPC_dvcc09, beliefPenaltyMPC_yy08, beliefPenaltyMPC_bmy08);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_5(beliefPenaltyMPC_Ld08, beliefPenaltyMPC_bmy08, beliefPenaltyMPC_dvcc08);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_5_5(beliefPenaltyMPC_Lsd08, beliefPenaltyMPC_dvcc08, beliefPenaltyMPC_yy07, beliefPenaltyMPC_bmy07);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_5(beliefPenaltyMPC_Ld07, beliefPenaltyMPC_bmy07, beliefPenaltyMPC_dvcc07);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_5_5(beliefPenaltyMPC_Lsd07, beliefPenaltyMPC_dvcc07, beliefPenaltyMPC_yy06, beliefPenaltyMPC_bmy06);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_5(beliefPenaltyMPC_Ld06, beliefPenaltyMPC_bmy06, beliefPenaltyMPC_dvcc06);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_5_5(beliefPenaltyMPC_Lsd06, beliefPenaltyMPC_dvcc06, beliefPenaltyMPC_yy05, beliefPenaltyMPC_bmy05);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_5(beliefPenaltyMPC_Ld05, beliefPenaltyMPC_bmy05, beliefPenaltyMPC_dvcc05);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_5_5(beliefPenaltyMPC_Lsd05, beliefPenaltyMPC_dvcc05, beliefPenaltyMPC_yy04, beliefPenaltyMPC_bmy04);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_5(beliefPenaltyMPC_Ld04, beliefPenaltyMPC_bmy04, beliefPenaltyMPC_dvcc04);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_5_5(beliefPenaltyMPC_Lsd04, beliefPenaltyMPC_dvcc04, beliefPenaltyMPC_yy03, beliefPenaltyMPC_bmy03);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_5(beliefPenaltyMPC_Ld03, beliefPenaltyMPC_bmy03, beliefPenaltyMPC_dvcc03);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_5_5(beliefPenaltyMPC_Lsd03, beliefPenaltyMPC_dvcc03, beliefPenaltyMPC_yy02, beliefPenaltyMPC_bmy02);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_5(beliefPenaltyMPC_Ld02, beliefPenaltyMPC_bmy02, beliefPenaltyMPC_dvcc02);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_5_5(beliefPenaltyMPC_Lsd02, beliefPenaltyMPC_dvcc02, beliefPenaltyMPC_yy01, beliefPenaltyMPC_bmy01);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_5(beliefPenaltyMPC_Ld01, beliefPenaltyMPC_bmy01, beliefPenaltyMPC_dvcc01);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_5_10(beliefPenaltyMPC_Lsd01, beliefPenaltyMPC_dvcc01, beliefPenaltyMPC_yy00, beliefPenaltyMPC_bmy00);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_10(beliefPenaltyMPC_Ld00, beliefPenaltyMPC_bmy00, beliefPenaltyMPC_dvcc00);
beliefPenaltyMPC_LA_DENSE_MTVM_10_17(params->C1, beliefPenaltyMPC_dvcc00, beliefPenaltyMPC_grad_eq00);
beliefPenaltyMPC_LA_DENSE_MTVM2_5_17_10(params->C2, beliefPenaltyMPC_dvcc01, beliefPenaltyMPC_D01, beliefPenaltyMPC_dvcc00, beliefPenaltyMPC_grad_eq01);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_5_17_5(params->C3, beliefPenaltyMPC_dvcc02, beliefPenaltyMPC_D02, beliefPenaltyMPC_dvcc01, beliefPenaltyMPC_grad_eq02);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_5_17_5(params->C4, beliefPenaltyMPC_dvcc03, beliefPenaltyMPC_D02, beliefPenaltyMPC_dvcc02, beliefPenaltyMPC_grad_eq03);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_5_17_5(params->C5, beliefPenaltyMPC_dvcc04, beliefPenaltyMPC_D02, beliefPenaltyMPC_dvcc03, beliefPenaltyMPC_grad_eq04);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_5_17_5(params->C6, beliefPenaltyMPC_dvcc05, beliefPenaltyMPC_D02, beliefPenaltyMPC_dvcc04, beliefPenaltyMPC_grad_eq05);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_5_17_5(params->C7, beliefPenaltyMPC_dvcc06, beliefPenaltyMPC_D02, beliefPenaltyMPC_dvcc05, beliefPenaltyMPC_grad_eq06);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_5_17_5(params->C8, beliefPenaltyMPC_dvcc07, beliefPenaltyMPC_D02, beliefPenaltyMPC_dvcc06, beliefPenaltyMPC_grad_eq07);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_5_17_5(params->C9, beliefPenaltyMPC_dvcc08, beliefPenaltyMPC_D02, beliefPenaltyMPC_dvcc07, beliefPenaltyMPC_grad_eq08);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_5_17_5(params->C10, beliefPenaltyMPC_dvcc09, beliefPenaltyMPC_D02, beliefPenaltyMPC_dvcc08, beliefPenaltyMPC_grad_eq09);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_5_17_5(params->C11, beliefPenaltyMPC_dvcc10, beliefPenaltyMPC_D02, beliefPenaltyMPC_dvcc09, beliefPenaltyMPC_grad_eq10);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_5_17_5(params->C12, beliefPenaltyMPC_dvcc11, beliefPenaltyMPC_D02, beliefPenaltyMPC_dvcc10, beliefPenaltyMPC_grad_eq11);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_5_17_5(params->C13, beliefPenaltyMPC_dvcc12, beliefPenaltyMPC_D02, beliefPenaltyMPC_dvcc11, beliefPenaltyMPC_grad_eq12);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_5_17_5(params->C14, beliefPenaltyMPC_dvcc13, beliefPenaltyMPC_D02, beliefPenaltyMPC_dvcc12, beliefPenaltyMPC_grad_eq13);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_5_17_5(params->C15, beliefPenaltyMPC_dvcc14, beliefPenaltyMPC_D02, beliefPenaltyMPC_dvcc13, beliefPenaltyMPC_grad_eq14);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_5_17_5(params->C16, beliefPenaltyMPC_dvcc15, beliefPenaltyMPC_D02, beliefPenaltyMPC_dvcc14, beliefPenaltyMPC_grad_eq15);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_5_17_5(params->C17, beliefPenaltyMPC_dvcc16, beliefPenaltyMPC_D02, beliefPenaltyMPC_dvcc15, beliefPenaltyMPC_grad_eq16);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_5_17_5(params->C18, beliefPenaltyMPC_dvcc17, beliefPenaltyMPC_D02, beliefPenaltyMPC_dvcc16, beliefPenaltyMPC_grad_eq17);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_5_17_5(params->C19, beliefPenaltyMPC_dvcc18, beliefPenaltyMPC_D02, beliefPenaltyMPC_dvcc17, beliefPenaltyMPC_grad_eq18);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_5_17_5(params->C20, beliefPenaltyMPC_dvcc19, beliefPenaltyMPC_D02, beliefPenaltyMPC_dvcc18, beliefPenaltyMPC_grad_eq19);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_5_17_5(params->C21, beliefPenaltyMPC_dvcc20, beliefPenaltyMPC_D02, beliefPenaltyMPC_dvcc19, beliefPenaltyMPC_grad_eq20);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_5_17_5(params->C22, beliefPenaltyMPC_dvcc21, beliefPenaltyMPC_D02, beliefPenaltyMPC_dvcc20, beliefPenaltyMPC_grad_eq21);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_5_17_5(params->C23, beliefPenaltyMPC_dvcc22, beliefPenaltyMPC_D02, beliefPenaltyMPC_dvcc21, beliefPenaltyMPC_grad_eq22);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_5_17_5(params->C24, beliefPenaltyMPC_dvcc23, beliefPenaltyMPC_D02, beliefPenaltyMPC_dvcc22, beliefPenaltyMPC_grad_eq23);
beliefPenaltyMPC_LA_DIAGZERO_MTVM_5_5(beliefPenaltyMPC_D24, beliefPenaltyMPC_dvcc23, beliefPenaltyMPC_grad_eq24);
beliefPenaltyMPC_LA_VSUB_413(beliefPenaltyMPC_rd, beliefPenaltyMPC_grad_eq, beliefPenaltyMPC_rd);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_17(beliefPenaltyMPC_Phi00, beliefPenaltyMPC_rd00, beliefPenaltyMPC_dzcc00);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_17(beliefPenaltyMPC_Phi01, beliefPenaltyMPC_rd01, beliefPenaltyMPC_dzcc01);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_17(beliefPenaltyMPC_Phi02, beliefPenaltyMPC_rd02, beliefPenaltyMPC_dzcc02);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_17(beliefPenaltyMPC_Phi03, beliefPenaltyMPC_rd03, beliefPenaltyMPC_dzcc03);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_17(beliefPenaltyMPC_Phi04, beliefPenaltyMPC_rd04, beliefPenaltyMPC_dzcc04);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_17(beliefPenaltyMPC_Phi05, beliefPenaltyMPC_rd05, beliefPenaltyMPC_dzcc05);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_17(beliefPenaltyMPC_Phi06, beliefPenaltyMPC_rd06, beliefPenaltyMPC_dzcc06);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_17(beliefPenaltyMPC_Phi07, beliefPenaltyMPC_rd07, beliefPenaltyMPC_dzcc07);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_17(beliefPenaltyMPC_Phi08, beliefPenaltyMPC_rd08, beliefPenaltyMPC_dzcc08);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_17(beliefPenaltyMPC_Phi09, beliefPenaltyMPC_rd09, beliefPenaltyMPC_dzcc09);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_17(beliefPenaltyMPC_Phi10, beliefPenaltyMPC_rd10, beliefPenaltyMPC_dzcc10);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_17(beliefPenaltyMPC_Phi11, beliefPenaltyMPC_rd11, beliefPenaltyMPC_dzcc11);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_17(beliefPenaltyMPC_Phi12, beliefPenaltyMPC_rd12, beliefPenaltyMPC_dzcc12);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_17(beliefPenaltyMPC_Phi13, beliefPenaltyMPC_rd13, beliefPenaltyMPC_dzcc13);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_17(beliefPenaltyMPC_Phi14, beliefPenaltyMPC_rd14, beliefPenaltyMPC_dzcc14);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_17(beliefPenaltyMPC_Phi15, beliefPenaltyMPC_rd15, beliefPenaltyMPC_dzcc15);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_17(beliefPenaltyMPC_Phi16, beliefPenaltyMPC_rd16, beliefPenaltyMPC_dzcc16);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_17(beliefPenaltyMPC_Phi17, beliefPenaltyMPC_rd17, beliefPenaltyMPC_dzcc17);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_17(beliefPenaltyMPC_Phi18, beliefPenaltyMPC_rd18, beliefPenaltyMPC_dzcc18);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_17(beliefPenaltyMPC_Phi19, beliefPenaltyMPC_rd19, beliefPenaltyMPC_dzcc19);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_17(beliefPenaltyMPC_Phi20, beliefPenaltyMPC_rd20, beliefPenaltyMPC_dzcc20);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_17(beliefPenaltyMPC_Phi21, beliefPenaltyMPC_rd21, beliefPenaltyMPC_dzcc21);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_17(beliefPenaltyMPC_Phi22, beliefPenaltyMPC_rd22, beliefPenaltyMPC_dzcc22);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_17(beliefPenaltyMPC_Phi23, beliefPenaltyMPC_rd23, beliefPenaltyMPC_dzcc23);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_5(beliefPenaltyMPC_Phi24, beliefPenaltyMPC_rd24, beliefPenaltyMPC_dzcc24);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_17(beliefPenaltyMPC_ccrhsl00, beliefPenaltyMPC_slb00, beliefPenaltyMPC_llbbyslb00, beliefPenaltyMPC_dzcc00, beliefPenaltyMPC_lbIdx00, beliefPenaltyMPC_dllbcc00);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_7(beliefPenaltyMPC_ccrhsub00, beliefPenaltyMPC_sub00, beliefPenaltyMPC_lubbysub00, beliefPenaltyMPC_dzcc00, beliefPenaltyMPC_ubIdx00, beliefPenaltyMPC_dlubcc00);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_17(beliefPenaltyMPC_ccrhsl01, beliefPenaltyMPC_slb01, beliefPenaltyMPC_llbbyslb01, beliefPenaltyMPC_dzcc01, beliefPenaltyMPC_lbIdx01, beliefPenaltyMPC_dllbcc01);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_7(beliefPenaltyMPC_ccrhsub01, beliefPenaltyMPC_sub01, beliefPenaltyMPC_lubbysub01, beliefPenaltyMPC_dzcc01, beliefPenaltyMPC_ubIdx01, beliefPenaltyMPC_dlubcc01);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_17(beliefPenaltyMPC_ccrhsl02, beliefPenaltyMPC_slb02, beliefPenaltyMPC_llbbyslb02, beliefPenaltyMPC_dzcc02, beliefPenaltyMPC_lbIdx02, beliefPenaltyMPC_dllbcc02);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_7(beliefPenaltyMPC_ccrhsub02, beliefPenaltyMPC_sub02, beliefPenaltyMPC_lubbysub02, beliefPenaltyMPC_dzcc02, beliefPenaltyMPC_ubIdx02, beliefPenaltyMPC_dlubcc02);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_17(beliefPenaltyMPC_ccrhsl03, beliefPenaltyMPC_slb03, beliefPenaltyMPC_llbbyslb03, beliefPenaltyMPC_dzcc03, beliefPenaltyMPC_lbIdx03, beliefPenaltyMPC_dllbcc03);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_7(beliefPenaltyMPC_ccrhsub03, beliefPenaltyMPC_sub03, beliefPenaltyMPC_lubbysub03, beliefPenaltyMPC_dzcc03, beliefPenaltyMPC_ubIdx03, beliefPenaltyMPC_dlubcc03);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_17(beliefPenaltyMPC_ccrhsl04, beliefPenaltyMPC_slb04, beliefPenaltyMPC_llbbyslb04, beliefPenaltyMPC_dzcc04, beliefPenaltyMPC_lbIdx04, beliefPenaltyMPC_dllbcc04);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_7(beliefPenaltyMPC_ccrhsub04, beliefPenaltyMPC_sub04, beliefPenaltyMPC_lubbysub04, beliefPenaltyMPC_dzcc04, beliefPenaltyMPC_ubIdx04, beliefPenaltyMPC_dlubcc04);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_17(beliefPenaltyMPC_ccrhsl05, beliefPenaltyMPC_slb05, beliefPenaltyMPC_llbbyslb05, beliefPenaltyMPC_dzcc05, beliefPenaltyMPC_lbIdx05, beliefPenaltyMPC_dllbcc05);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_7(beliefPenaltyMPC_ccrhsub05, beliefPenaltyMPC_sub05, beliefPenaltyMPC_lubbysub05, beliefPenaltyMPC_dzcc05, beliefPenaltyMPC_ubIdx05, beliefPenaltyMPC_dlubcc05);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_17(beliefPenaltyMPC_ccrhsl06, beliefPenaltyMPC_slb06, beliefPenaltyMPC_llbbyslb06, beliefPenaltyMPC_dzcc06, beliefPenaltyMPC_lbIdx06, beliefPenaltyMPC_dllbcc06);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_7(beliefPenaltyMPC_ccrhsub06, beliefPenaltyMPC_sub06, beliefPenaltyMPC_lubbysub06, beliefPenaltyMPC_dzcc06, beliefPenaltyMPC_ubIdx06, beliefPenaltyMPC_dlubcc06);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_17(beliefPenaltyMPC_ccrhsl07, beliefPenaltyMPC_slb07, beliefPenaltyMPC_llbbyslb07, beliefPenaltyMPC_dzcc07, beliefPenaltyMPC_lbIdx07, beliefPenaltyMPC_dllbcc07);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_7(beliefPenaltyMPC_ccrhsub07, beliefPenaltyMPC_sub07, beliefPenaltyMPC_lubbysub07, beliefPenaltyMPC_dzcc07, beliefPenaltyMPC_ubIdx07, beliefPenaltyMPC_dlubcc07);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_17(beliefPenaltyMPC_ccrhsl08, beliefPenaltyMPC_slb08, beliefPenaltyMPC_llbbyslb08, beliefPenaltyMPC_dzcc08, beliefPenaltyMPC_lbIdx08, beliefPenaltyMPC_dllbcc08);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_7(beliefPenaltyMPC_ccrhsub08, beliefPenaltyMPC_sub08, beliefPenaltyMPC_lubbysub08, beliefPenaltyMPC_dzcc08, beliefPenaltyMPC_ubIdx08, beliefPenaltyMPC_dlubcc08);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_17(beliefPenaltyMPC_ccrhsl09, beliefPenaltyMPC_slb09, beliefPenaltyMPC_llbbyslb09, beliefPenaltyMPC_dzcc09, beliefPenaltyMPC_lbIdx09, beliefPenaltyMPC_dllbcc09);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_7(beliefPenaltyMPC_ccrhsub09, beliefPenaltyMPC_sub09, beliefPenaltyMPC_lubbysub09, beliefPenaltyMPC_dzcc09, beliefPenaltyMPC_ubIdx09, beliefPenaltyMPC_dlubcc09);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_17(beliefPenaltyMPC_ccrhsl10, beliefPenaltyMPC_slb10, beliefPenaltyMPC_llbbyslb10, beliefPenaltyMPC_dzcc10, beliefPenaltyMPC_lbIdx10, beliefPenaltyMPC_dllbcc10);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_7(beliefPenaltyMPC_ccrhsub10, beliefPenaltyMPC_sub10, beliefPenaltyMPC_lubbysub10, beliefPenaltyMPC_dzcc10, beliefPenaltyMPC_ubIdx10, beliefPenaltyMPC_dlubcc10);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_17(beliefPenaltyMPC_ccrhsl11, beliefPenaltyMPC_slb11, beliefPenaltyMPC_llbbyslb11, beliefPenaltyMPC_dzcc11, beliefPenaltyMPC_lbIdx11, beliefPenaltyMPC_dllbcc11);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_7(beliefPenaltyMPC_ccrhsub11, beliefPenaltyMPC_sub11, beliefPenaltyMPC_lubbysub11, beliefPenaltyMPC_dzcc11, beliefPenaltyMPC_ubIdx11, beliefPenaltyMPC_dlubcc11);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_17(beliefPenaltyMPC_ccrhsl12, beliefPenaltyMPC_slb12, beliefPenaltyMPC_llbbyslb12, beliefPenaltyMPC_dzcc12, beliefPenaltyMPC_lbIdx12, beliefPenaltyMPC_dllbcc12);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_7(beliefPenaltyMPC_ccrhsub12, beliefPenaltyMPC_sub12, beliefPenaltyMPC_lubbysub12, beliefPenaltyMPC_dzcc12, beliefPenaltyMPC_ubIdx12, beliefPenaltyMPC_dlubcc12);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_17(beliefPenaltyMPC_ccrhsl13, beliefPenaltyMPC_slb13, beliefPenaltyMPC_llbbyslb13, beliefPenaltyMPC_dzcc13, beliefPenaltyMPC_lbIdx13, beliefPenaltyMPC_dllbcc13);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_7(beliefPenaltyMPC_ccrhsub13, beliefPenaltyMPC_sub13, beliefPenaltyMPC_lubbysub13, beliefPenaltyMPC_dzcc13, beliefPenaltyMPC_ubIdx13, beliefPenaltyMPC_dlubcc13);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_17(beliefPenaltyMPC_ccrhsl14, beliefPenaltyMPC_slb14, beliefPenaltyMPC_llbbyslb14, beliefPenaltyMPC_dzcc14, beliefPenaltyMPC_lbIdx14, beliefPenaltyMPC_dllbcc14);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_7(beliefPenaltyMPC_ccrhsub14, beliefPenaltyMPC_sub14, beliefPenaltyMPC_lubbysub14, beliefPenaltyMPC_dzcc14, beliefPenaltyMPC_ubIdx14, beliefPenaltyMPC_dlubcc14);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_17(beliefPenaltyMPC_ccrhsl15, beliefPenaltyMPC_slb15, beliefPenaltyMPC_llbbyslb15, beliefPenaltyMPC_dzcc15, beliefPenaltyMPC_lbIdx15, beliefPenaltyMPC_dllbcc15);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_7(beliefPenaltyMPC_ccrhsub15, beliefPenaltyMPC_sub15, beliefPenaltyMPC_lubbysub15, beliefPenaltyMPC_dzcc15, beliefPenaltyMPC_ubIdx15, beliefPenaltyMPC_dlubcc15);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_17(beliefPenaltyMPC_ccrhsl16, beliefPenaltyMPC_slb16, beliefPenaltyMPC_llbbyslb16, beliefPenaltyMPC_dzcc16, beliefPenaltyMPC_lbIdx16, beliefPenaltyMPC_dllbcc16);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_7(beliefPenaltyMPC_ccrhsub16, beliefPenaltyMPC_sub16, beliefPenaltyMPC_lubbysub16, beliefPenaltyMPC_dzcc16, beliefPenaltyMPC_ubIdx16, beliefPenaltyMPC_dlubcc16);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_17(beliefPenaltyMPC_ccrhsl17, beliefPenaltyMPC_slb17, beliefPenaltyMPC_llbbyslb17, beliefPenaltyMPC_dzcc17, beliefPenaltyMPC_lbIdx17, beliefPenaltyMPC_dllbcc17);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_7(beliefPenaltyMPC_ccrhsub17, beliefPenaltyMPC_sub17, beliefPenaltyMPC_lubbysub17, beliefPenaltyMPC_dzcc17, beliefPenaltyMPC_ubIdx17, beliefPenaltyMPC_dlubcc17);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_17(beliefPenaltyMPC_ccrhsl18, beliefPenaltyMPC_slb18, beliefPenaltyMPC_llbbyslb18, beliefPenaltyMPC_dzcc18, beliefPenaltyMPC_lbIdx18, beliefPenaltyMPC_dllbcc18);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_7(beliefPenaltyMPC_ccrhsub18, beliefPenaltyMPC_sub18, beliefPenaltyMPC_lubbysub18, beliefPenaltyMPC_dzcc18, beliefPenaltyMPC_ubIdx18, beliefPenaltyMPC_dlubcc18);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_17(beliefPenaltyMPC_ccrhsl19, beliefPenaltyMPC_slb19, beliefPenaltyMPC_llbbyslb19, beliefPenaltyMPC_dzcc19, beliefPenaltyMPC_lbIdx19, beliefPenaltyMPC_dllbcc19);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_7(beliefPenaltyMPC_ccrhsub19, beliefPenaltyMPC_sub19, beliefPenaltyMPC_lubbysub19, beliefPenaltyMPC_dzcc19, beliefPenaltyMPC_ubIdx19, beliefPenaltyMPC_dlubcc19);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_17(beliefPenaltyMPC_ccrhsl20, beliefPenaltyMPC_slb20, beliefPenaltyMPC_llbbyslb20, beliefPenaltyMPC_dzcc20, beliefPenaltyMPC_lbIdx20, beliefPenaltyMPC_dllbcc20);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_7(beliefPenaltyMPC_ccrhsub20, beliefPenaltyMPC_sub20, beliefPenaltyMPC_lubbysub20, beliefPenaltyMPC_dzcc20, beliefPenaltyMPC_ubIdx20, beliefPenaltyMPC_dlubcc20);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_17(beliefPenaltyMPC_ccrhsl21, beliefPenaltyMPC_slb21, beliefPenaltyMPC_llbbyslb21, beliefPenaltyMPC_dzcc21, beliefPenaltyMPC_lbIdx21, beliefPenaltyMPC_dllbcc21);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_7(beliefPenaltyMPC_ccrhsub21, beliefPenaltyMPC_sub21, beliefPenaltyMPC_lubbysub21, beliefPenaltyMPC_dzcc21, beliefPenaltyMPC_ubIdx21, beliefPenaltyMPC_dlubcc21);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_17(beliefPenaltyMPC_ccrhsl22, beliefPenaltyMPC_slb22, beliefPenaltyMPC_llbbyslb22, beliefPenaltyMPC_dzcc22, beliefPenaltyMPC_lbIdx22, beliefPenaltyMPC_dllbcc22);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_7(beliefPenaltyMPC_ccrhsub22, beliefPenaltyMPC_sub22, beliefPenaltyMPC_lubbysub22, beliefPenaltyMPC_dzcc22, beliefPenaltyMPC_ubIdx22, beliefPenaltyMPC_dlubcc22);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_17(beliefPenaltyMPC_ccrhsl23, beliefPenaltyMPC_slb23, beliefPenaltyMPC_llbbyslb23, beliefPenaltyMPC_dzcc23, beliefPenaltyMPC_lbIdx23, beliefPenaltyMPC_dllbcc23);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_7(beliefPenaltyMPC_ccrhsub23, beliefPenaltyMPC_sub23, beliefPenaltyMPC_lubbysub23, beliefPenaltyMPC_dzcc23, beliefPenaltyMPC_ubIdx23, beliefPenaltyMPC_dlubcc23);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_5(beliefPenaltyMPC_ccrhsl24, beliefPenaltyMPC_slb24, beliefPenaltyMPC_llbbyslb24, beliefPenaltyMPC_dzcc24, beliefPenaltyMPC_lbIdx24, beliefPenaltyMPC_dllbcc24);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(beliefPenaltyMPC_ccrhsub24, beliefPenaltyMPC_sub24, beliefPenaltyMPC_lubbysub24, beliefPenaltyMPC_dzcc24, beliefPenaltyMPC_ubIdx24, beliefPenaltyMPC_dlubcc24);
beliefPenaltyMPC_LA_VSUB7_586(beliefPenaltyMPC_l, beliefPenaltyMPC_ccrhs, beliefPenaltyMPC_s, beliefPenaltyMPC_dl_cc, beliefPenaltyMPC_ds_cc);
beliefPenaltyMPC_LA_VADD_413(beliefPenaltyMPC_dz_cc, beliefPenaltyMPC_dz_aff);
beliefPenaltyMPC_LA_VADD_125(beliefPenaltyMPC_dv_cc, beliefPenaltyMPC_dv_aff);
beliefPenaltyMPC_LA_VADD_586(beliefPenaltyMPC_dl_cc, beliefPenaltyMPC_dl_aff);
beliefPenaltyMPC_LA_VADD_586(beliefPenaltyMPC_ds_cc, beliefPenaltyMPC_ds_aff);
info->lsit_cc = beliefPenaltyMPC_LINESEARCH_BACKTRACKING_COMBINED(beliefPenaltyMPC_z, beliefPenaltyMPC_v, beliefPenaltyMPC_l, beliefPenaltyMPC_s, beliefPenaltyMPC_dz_cc, beliefPenaltyMPC_dv_cc, beliefPenaltyMPC_dl_cc, beliefPenaltyMPC_ds_cc, &info->step_cc, &info->mu);
if( info->lsit_cc == beliefPenaltyMPC_NOPROGRESS ){
exitcode = beliefPenaltyMPC_NOPROGRESS; break;
}
info->it++;
}
output->z1[0] = beliefPenaltyMPC_z00[0];
output->z1[1] = beliefPenaltyMPC_z00[1];
output->z1[2] = beliefPenaltyMPC_z00[2];
output->z1[3] = beliefPenaltyMPC_z00[3];
output->z1[4] = beliefPenaltyMPC_z00[4];
output->z1[5] = beliefPenaltyMPC_z00[5];
output->z1[6] = beliefPenaltyMPC_z00[6];
output->z2[0] = beliefPenaltyMPC_z01[0];
output->z2[1] = beliefPenaltyMPC_z01[1];
output->z2[2] = beliefPenaltyMPC_z01[2];
output->z2[3] = beliefPenaltyMPC_z01[3];
output->z2[4] = beliefPenaltyMPC_z01[4];
output->z2[5] = beliefPenaltyMPC_z01[5];
output->z2[6] = beliefPenaltyMPC_z01[6];
output->z3[0] = beliefPenaltyMPC_z02[0];
output->z3[1] = beliefPenaltyMPC_z02[1];
output->z3[2] = beliefPenaltyMPC_z02[2];
output->z3[3] = beliefPenaltyMPC_z02[3];
output->z3[4] = beliefPenaltyMPC_z02[4];
output->z3[5] = beliefPenaltyMPC_z02[5];
output->z3[6] = beliefPenaltyMPC_z02[6];
output->z4[0] = beliefPenaltyMPC_z03[0];
output->z4[1] = beliefPenaltyMPC_z03[1];
output->z4[2] = beliefPenaltyMPC_z03[2];
output->z4[3] = beliefPenaltyMPC_z03[3];
output->z4[4] = beliefPenaltyMPC_z03[4];
output->z4[5] = beliefPenaltyMPC_z03[5];
output->z4[6] = beliefPenaltyMPC_z03[6];
output->z5[0] = beliefPenaltyMPC_z04[0];
output->z5[1] = beliefPenaltyMPC_z04[1];
output->z5[2] = beliefPenaltyMPC_z04[2];
output->z5[3] = beliefPenaltyMPC_z04[3];
output->z5[4] = beliefPenaltyMPC_z04[4];
output->z5[5] = beliefPenaltyMPC_z04[5];
output->z5[6] = beliefPenaltyMPC_z04[6];
output->z6[0] = beliefPenaltyMPC_z05[0];
output->z6[1] = beliefPenaltyMPC_z05[1];
output->z6[2] = beliefPenaltyMPC_z05[2];
output->z6[3] = beliefPenaltyMPC_z05[3];
output->z6[4] = beliefPenaltyMPC_z05[4];
output->z6[5] = beliefPenaltyMPC_z05[5];
output->z6[6] = beliefPenaltyMPC_z05[6];
output->z7[0] = beliefPenaltyMPC_z06[0];
output->z7[1] = beliefPenaltyMPC_z06[1];
output->z7[2] = beliefPenaltyMPC_z06[2];
output->z7[3] = beliefPenaltyMPC_z06[3];
output->z7[4] = beliefPenaltyMPC_z06[4];
output->z7[5] = beliefPenaltyMPC_z06[5];
output->z7[6] = beliefPenaltyMPC_z06[6];
output->z8[0] = beliefPenaltyMPC_z07[0];
output->z8[1] = beliefPenaltyMPC_z07[1];
output->z8[2] = beliefPenaltyMPC_z07[2];
output->z8[3] = beliefPenaltyMPC_z07[3];
output->z8[4] = beliefPenaltyMPC_z07[4];
output->z8[5] = beliefPenaltyMPC_z07[5];
output->z8[6] = beliefPenaltyMPC_z07[6];
output->z9[0] = beliefPenaltyMPC_z08[0];
output->z9[1] = beliefPenaltyMPC_z08[1];
output->z9[2] = beliefPenaltyMPC_z08[2];
output->z9[3] = beliefPenaltyMPC_z08[3];
output->z9[4] = beliefPenaltyMPC_z08[4];
output->z9[5] = beliefPenaltyMPC_z08[5];
output->z9[6] = beliefPenaltyMPC_z08[6];
output->z10[0] = beliefPenaltyMPC_z09[0];
output->z10[1] = beliefPenaltyMPC_z09[1];
output->z10[2] = beliefPenaltyMPC_z09[2];
output->z10[3] = beliefPenaltyMPC_z09[3];
output->z10[4] = beliefPenaltyMPC_z09[4];
output->z10[5] = beliefPenaltyMPC_z09[5];
output->z10[6] = beliefPenaltyMPC_z09[6];
output->z11[0] = beliefPenaltyMPC_z10[0];
output->z11[1] = beliefPenaltyMPC_z10[1];
output->z11[2] = beliefPenaltyMPC_z10[2];
output->z11[3] = beliefPenaltyMPC_z10[3];
output->z11[4] = beliefPenaltyMPC_z10[4];
output->z11[5] = beliefPenaltyMPC_z10[5];
output->z11[6] = beliefPenaltyMPC_z10[6];
output->z12[0] = beliefPenaltyMPC_z11[0];
output->z12[1] = beliefPenaltyMPC_z11[1];
output->z12[2] = beliefPenaltyMPC_z11[2];
output->z12[3] = beliefPenaltyMPC_z11[3];
output->z12[4] = beliefPenaltyMPC_z11[4];
output->z12[5] = beliefPenaltyMPC_z11[5];
output->z12[6] = beliefPenaltyMPC_z11[6];
output->z13[0] = beliefPenaltyMPC_z12[0];
output->z13[1] = beliefPenaltyMPC_z12[1];
output->z13[2] = beliefPenaltyMPC_z12[2];
output->z13[3] = beliefPenaltyMPC_z12[3];
output->z13[4] = beliefPenaltyMPC_z12[4];
output->z13[5] = beliefPenaltyMPC_z12[5];
output->z13[6] = beliefPenaltyMPC_z12[6];
output->z14[0] = beliefPenaltyMPC_z13[0];
output->z14[1] = beliefPenaltyMPC_z13[1];
output->z14[2] = beliefPenaltyMPC_z13[2];
output->z14[3] = beliefPenaltyMPC_z13[3];
output->z14[4] = beliefPenaltyMPC_z13[4];
output->z14[5] = beliefPenaltyMPC_z13[5];
output->z14[6] = beliefPenaltyMPC_z13[6];
output->z15[0] = beliefPenaltyMPC_z14[0];
output->z15[1] = beliefPenaltyMPC_z14[1];
output->z15[2] = beliefPenaltyMPC_z14[2];
output->z15[3] = beliefPenaltyMPC_z14[3];
output->z15[4] = beliefPenaltyMPC_z14[4];
output->z15[5] = beliefPenaltyMPC_z14[5];
output->z15[6] = beliefPenaltyMPC_z14[6];
output->z16[0] = beliefPenaltyMPC_z15[0];
output->z16[1] = beliefPenaltyMPC_z15[1];
output->z16[2] = beliefPenaltyMPC_z15[2];
output->z16[3] = beliefPenaltyMPC_z15[3];
output->z16[4] = beliefPenaltyMPC_z15[4];
output->z16[5] = beliefPenaltyMPC_z15[5];
output->z16[6] = beliefPenaltyMPC_z15[6];
output->z17[0] = beliefPenaltyMPC_z16[0];
output->z17[1] = beliefPenaltyMPC_z16[1];
output->z17[2] = beliefPenaltyMPC_z16[2];
output->z17[3] = beliefPenaltyMPC_z16[3];
output->z17[4] = beliefPenaltyMPC_z16[4];
output->z17[5] = beliefPenaltyMPC_z16[5];
output->z17[6] = beliefPenaltyMPC_z16[6];
output->z18[0] = beliefPenaltyMPC_z17[0];
output->z18[1] = beliefPenaltyMPC_z17[1];
output->z18[2] = beliefPenaltyMPC_z17[2];
output->z18[3] = beliefPenaltyMPC_z17[3];
output->z18[4] = beliefPenaltyMPC_z17[4];
output->z18[5] = beliefPenaltyMPC_z17[5];
output->z18[6] = beliefPenaltyMPC_z17[6];
output->z19[0] = beliefPenaltyMPC_z18[0];
output->z19[1] = beliefPenaltyMPC_z18[1];
output->z19[2] = beliefPenaltyMPC_z18[2];
output->z19[3] = beliefPenaltyMPC_z18[3];
output->z19[4] = beliefPenaltyMPC_z18[4];
output->z19[5] = beliefPenaltyMPC_z18[5];
output->z19[6] = beliefPenaltyMPC_z18[6];
output->z20[0] = beliefPenaltyMPC_z19[0];
output->z20[1] = beliefPenaltyMPC_z19[1];
output->z20[2] = beliefPenaltyMPC_z19[2];
output->z20[3] = beliefPenaltyMPC_z19[3];
output->z20[4] = beliefPenaltyMPC_z19[4];
output->z20[5] = beliefPenaltyMPC_z19[5];
output->z20[6] = beliefPenaltyMPC_z19[6];
output->z21[0] = beliefPenaltyMPC_z20[0];
output->z21[1] = beliefPenaltyMPC_z20[1];
output->z21[2] = beliefPenaltyMPC_z20[2];
output->z21[3] = beliefPenaltyMPC_z20[3];
output->z21[4] = beliefPenaltyMPC_z20[4];
output->z21[5] = beliefPenaltyMPC_z20[5];
output->z21[6] = beliefPenaltyMPC_z20[6];
output->z22[0] = beliefPenaltyMPC_z21[0];
output->z22[1] = beliefPenaltyMPC_z21[1];
output->z22[2] = beliefPenaltyMPC_z21[2];
output->z22[3] = beliefPenaltyMPC_z21[3];
output->z22[4] = beliefPenaltyMPC_z21[4];
output->z22[5] = beliefPenaltyMPC_z21[5];
output->z22[6] = beliefPenaltyMPC_z21[6];
output->z23[0] = beliefPenaltyMPC_z22[0];
output->z23[1] = beliefPenaltyMPC_z22[1];
output->z23[2] = beliefPenaltyMPC_z22[2];
output->z23[3] = beliefPenaltyMPC_z22[3];
output->z23[4] = beliefPenaltyMPC_z22[4];
output->z23[5] = beliefPenaltyMPC_z22[5];
output->z23[6] = beliefPenaltyMPC_z22[6];
output->z24[0] = beliefPenaltyMPC_z23[0];
output->z24[1] = beliefPenaltyMPC_z23[1];
output->z24[2] = beliefPenaltyMPC_z23[2];
output->z24[3] = beliefPenaltyMPC_z23[3];
output->z24[4] = beliefPenaltyMPC_z23[4];
output->z24[5] = beliefPenaltyMPC_z23[5];
output->z24[6] = beliefPenaltyMPC_z23[6];
output->z25[0] = beliefPenaltyMPC_z24[0];
output->z25[1] = beliefPenaltyMPC_z24[1];
output->z25[2] = beliefPenaltyMPC_z24[2];
output->z25[3] = beliefPenaltyMPC_z24[3];
output->z25[4] = beliefPenaltyMPC_z24[4];

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
