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
 * Initializes a vector of length 578 with a value.
 */
void beliefPenaltyMPC_LA_INITIALIZEVECTOR_578(beliefPenaltyMPC_FLOAT* vec, beliefPenaltyMPC_FLOAT value)
{
	int i;
	for( i=0; i<578; i++ )
	{
		vec[i] = value;
	}
}


/*
 * Initializes a vector of length 200 with a value.
 */
void beliefPenaltyMPC_LA_INITIALIZEVECTOR_200(beliefPenaltyMPC_FLOAT* vec, beliefPenaltyMPC_FLOAT value)
{
	int i;
	for( i=0; i<200; i++ )
	{
		vec[i] = value;
	}
}


/*
 * Initializes a vector of length 796 with a value.
 */
void beliefPenaltyMPC_LA_INITIALIZEVECTOR_796(beliefPenaltyMPC_FLOAT* vec, beliefPenaltyMPC_FLOAT value)
{
	int i;
	for( i=0; i<796; i++ )
	{
		vec[i] = value;
	}
}


/* 
 * Calculates a dot product and adds it to a variable: z += x'*y; 
 * This function is for vectors of length 796.
 */
void beliefPenaltyMPC_LA_DOTACC_796(beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *y, beliefPenaltyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<796; i++ ){
		*z += x[i]*y[i];
	}
}


/*
 * Calculates the gradient and the value for a quadratic function 0.5*z'*H*z + f'*z
 *
 * INPUTS:     H  - Symmetric Hessian, diag matrix of size [62 x 62]
 *             f  - column vector of size 62
 *             z  - column vector of size 62
 *
 * OUTPUTS: grad  - gradient at z (= H*z + f), column vector of size 62
 *          value <-- value + 0.5*z'*H*z + f'*z (value will be modified)
 */
void beliefPenaltyMPC_LA_DIAG_QUADFCN_62(beliefPenaltyMPC_FLOAT* H, beliefPenaltyMPC_FLOAT* f, beliefPenaltyMPC_FLOAT* z, beliefPenaltyMPC_FLOAT* grad, beliefPenaltyMPC_FLOAT* value)
{
	int i;
	beliefPenaltyMPC_FLOAT hz;	
	for( i=0; i<62; i++){
		hz = H[i]*z[i];
		grad[i] = hz + f[i];
		*value += 0.5*hz*z[i] + f[i]*z[i];
	}
}


/*
 * Calculates the gradient and the value for a quadratic function 0.5*z'*H*z + f'*z
 *
 * INPUTS:     H  - Symmetric Hessian, diag matrix of size [20 x 20]
 *             f  - column vector of size 20
 *             z  - column vector of size 20
 *
 * OUTPUTS: grad  - gradient at z (= H*z + f), column vector of size 20
 *          value <-- value + 0.5*z'*H*z + f'*z (value will be modified)
 */
void beliefPenaltyMPC_LA_DIAG_QUADFCN_20(beliefPenaltyMPC_FLOAT* H, beliefPenaltyMPC_FLOAT* f, beliefPenaltyMPC_FLOAT* z, beliefPenaltyMPC_FLOAT* grad, beliefPenaltyMPC_FLOAT* value)
{
	int i;
	beliefPenaltyMPC_FLOAT hz;	
	for( i=0; i<20; i++){
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
void beliefPenaltyMPC_LA_DIAGZERO_MVMSUB6_20(beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *b, beliefPenaltyMPC_FLOAT *l, beliefPenaltyMPC_FLOAT *r, beliefPenaltyMPC_FLOAT *z, beliefPenaltyMPC_FLOAT *y)
{
	int i;
	beliefPenaltyMPC_FLOAT Bu[20];
	beliefPenaltyMPC_FLOAT norm = *y;
	beliefPenaltyMPC_FLOAT lr = 0;

	/* do A*x + B*u first */
	for( i=0; i<20; i++ ){
		Bu[i] = B[i]*u[i];
	}	

	for( i=0; i<20; i++ ){
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
void beliefPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_20_62_62(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *b, beliefPenaltyMPC_FLOAT *l, beliefPenaltyMPC_FLOAT *r, beliefPenaltyMPC_FLOAT *z, beliefPenaltyMPC_FLOAT *y)
{
	int i;
	int j;
	int k = 0;
	beliefPenaltyMPC_FLOAT AxBu[20];
	beliefPenaltyMPC_FLOAT norm = *y;
	beliefPenaltyMPC_FLOAT lr = 0;

	/* do A*x + B*u first */
	for( i=0; i<20; i++ ){
		AxBu[i] = A[k++]*x[0] + B[i]*u[i];
	}	

	for( j=1; j<62; j++ ){		
		for( i=0; i<20; i++ ){
			AxBu[i] += A[k++]*x[j];
		}
	}

	for( i=0; i<20; i++ ){
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
void beliefPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_20_62_20(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *b, beliefPenaltyMPC_FLOAT *l, beliefPenaltyMPC_FLOAT *r, beliefPenaltyMPC_FLOAT *z, beliefPenaltyMPC_FLOAT *y)
{
	int i;
	int j;
	int k = 0;
	beliefPenaltyMPC_FLOAT AxBu[20];
	beliefPenaltyMPC_FLOAT norm = *y;
	beliefPenaltyMPC_FLOAT lr = 0;

	/* do A*x + B*u first */
	for( i=0; i<20; i++ ){
		AxBu[i] = A[k++]*x[0] + B[i]*u[i];
	}	

	for( j=1; j<62; j++ ){		
		for( i=0; i<20; i++ ){
			AxBu[i] += A[k++]*x[j];
		}
	}

	for( i=0; i<20; i++ ){
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
 * where A is of size [20 x 62] and stored in column major format.
 * and B is of size [20 x 62] and stored in diagzero format
 * Note the transposes of A and B!
 */
void beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_20_62_20(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *y, beliefPenaltyMPC_FLOAT *z)
{
	int i;
	int j;
	int k = 0;
	for( i=0; i<20; i++ ){
		z[i] = 0;
		for( j=0; j<20; j++ ){
			z[i] += A[k++]*x[j];
		}
		z[i] += B[i]*y[i];
	}
	for( i=20 ;i<62; i++ ){
		z[i] = 0;
		for( j=0; j<20; j++ ){
			z[i] += A[k++]*x[j];
		}
	}
}


/*
 * Matrix vector multiplication y = M'*x where M is of size [20 x 20]
 * and stored in diagzero format. Note the transpose of M!
 */
void beliefPenaltyMPC_LA_DIAGZERO_MTVM_20_20(beliefPenaltyMPC_FLOAT *M, beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *y)
{
	int i;
	for( i=0; i<20; i++ ){
		y[i] = M[i]*x[i];
	}
}


/*
 * Vector subtraction and addition.
 *	 Input: five vectors t, tidx, u, v, w and two scalars z and r
 *	 Output: y = t(tidx) - u + w
 *           z = z - v'*x;
 *           r = max([norm(y,inf), z]);
 * for vectors of length 62. Output z is of course scalar.
 */
void beliefPenaltyMPC_LA_VSUBADD3_62(beliefPenaltyMPC_FLOAT* t, beliefPenaltyMPC_FLOAT* u, int* uidx, beliefPenaltyMPC_FLOAT* v, beliefPenaltyMPC_FLOAT* w, beliefPenaltyMPC_FLOAT* y, beliefPenaltyMPC_FLOAT* z, beliefPenaltyMPC_FLOAT* r)
{
	int i;
	beliefPenaltyMPC_FLOAT norm = *r;
	beliefPenaltyMPC_FLOAT vx = 0;
	beliefPenaltyMPC_FLOAT x;
	for( i=0; i<62; i++){
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
 * for vectors of length 22. Output z is of course scalar.
 */
void beliefPenaltyMPC_LA_VSUBADD2_22(beliefPenaltyMPC_FLOAT* t, int* tidx, beliefPenaltyMPC_FLOAT* u, beliefPenaltyMPC_FLOAT* v, beliefPenaltyMPC_FLOAT* w, beliefPenaltyMPC_FLOAT* y, beliefPenaltyMPC_FLOAT* z, beliefPenaltyMPC_FLOAT* r)
{
	int i;
	beliefPenaltyMPC_FLOAT norm = *r;
	beliefPenaltyMPC_FLOAT vx = 0;
	beliefPenaltyMPC_FLOAT x;
	for( i=0; i<22; i++){
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
 * for vectors of length 20. Output z is of course scalar.
 */
void beliefPenaltyMPC_LA_VSUBADD3_20(beliefPenaltyMPC_FLOAT* t, beliefPenaltyMPC_FLOAT* u, int* uidx, beliefPenaltyMPC_FLOAT* v, beliefPenaltyMPC_FLOAT* w, beliefPenaltyMPC_FLOAT* y, beliefPenaltyMPC_FLOAT* z, beliefPenaltyMPC_FLOAT* r)
{
	int i;
	beliefPenaltyMPC_FLOAT norm = *r;
	beliefPenaltyMPC_FLOAT vx = 0;
	beliefPenaltyMPC_FLOAT x;
	for( i=0; i<20; i++){
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
 * for vectors of length 20. Output z is of course scalar.
 */
void beliefPenaltyMPC_LA_VSUBADD2_20(beliefPenaltyMPC_FLOAT* t, int* tidx, beliefPenaltyMPC_FLOAT* u, beliefPenaltyMPC_FLOAT* v, beliefPenaltyMPC_FLOAT* w, beliefPenaltyMPC_FLOAT* y, beliefPenaltyMPC_FLOAT* z, beliefPenaltyMPC_FLOAT* r)
{
	int i;
	beliefPenaltyMPC_FLOAT norm = *r;
	beliefPenaltyMPC_FLOAT vx = 0;
	beliefPenaltyMPC_FLOAT x;
	for( i=0; i<20; i++){
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
 * Special function for box constraints of length 62
 * Returns also L/S, a value that is often used elsewhere.
 */
void beliefPenaltyMPC_LA_INEQ_B_GRAD_62_62_22(beliefPenaltyMPC_FLOAT *lu, beliefPenaltyMPC_FLOAT *su, beliefPenaltyMPC_FLOAT *ru, beliefPenaltyMPC_FLOAT *ll, beliefPenaltyMPC_FLOAT *sl, beliefPenaltyMPC_FLOAT *rl, int* lbIdx, int* ubIdx, beliefPenaltyMPC_FLOAT *grad, beliefPenaltyMPC_FLOAT *lubysu, beliefPenaltyMPC_FLOAT *llbysl)
{
	int i;
	for( i=0; i<62; i++ ){
		grad[i] = 0;
	}
	for( i=0; i<62; i++ ){		
		llbysl[i] = ll[i] / sl[i];
		grad[lbIdx[i]] -= llbysl[i]*rl[i];
	}
	for( i=0; i<22; i++ ){
		lubysu[i] = lu[i] / su[i];
		grad[ubIdx[i]] += lubysu[i]*ru[i];
	}
}


/*
 * Computes inequality constraints gradient-
 * Special function for box constraints of length 20
 * Returns also L/S, a value that is often used elsewhere.
 */
void beliefPenaltyMPC_LA_INEQ_B_GRAD_20_20_20(beliefPenaltyMPC_FLOAT *lu, beliefPenaltyMPC_FLOAT *su, beliefPenaltyMPC_FLOAT *ru, beliefPenaltyMPC_FLOAT *ll, beliefPenaltyMPC_FLOAT *sl, beliefPenaltyMPC_FLOAT *rl, int* lbIdx, int* ubIdx, beliefPenaltyMPC_FLOAT *grad, beliefPenaltyMPC_FLOAT *lubysu, beliefPenaltyMPC_FLOAT *llbysl)
{
	int i;
	for( i=0; i<20; i++ ){
		grad[i] = 0;
	}
	for( i=0; i<20; i++ ){		
		llbysl[i] = ll[i] / sl[i];
		grad[lbIdx[i]] -= llbysl[i]*rl[i];
	}
	for( i=0; i<20; i++ ){
		lubysu[i] = lu[i] / su[i];
		grad[ubIdx[i]] += lubysu[i]*ru[i];
	}
}


/*
 * Addition of three vectors  z = u + w + v
 * of length 578.
 */
void beliefPenaltyMPC_LA_VVADD3_578(beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *v, beliefPenaltyMPC_FLOAT *w, beliefPenaltyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<578; i++ ){
		z[i] = u[i] + v[i] + w[i];
	}
}


/*
 * Special function to compute the diagonal cholesky factorization of the 
 * positive definite augmented Hessian for block size 62.
 *
 * Inputs: - H = diagonal cost Hessian in diagonal storage format
 *         - llbysl = L / S of lower bounds
 *         - lubysu = L / S of upper bounds
 *
 * Output: Phi = sqrt(H + diag(llbysl) + diag(lubysu))
 * where Phi is stored in diagonal storage format
 */
void beliefPenaltyMPC_LA_DIAG_CHOL_LBUB_62_62_22(beliefPenaltyMPC_FLOAT *H, beliefPenaltyMPC_FLOAT *llbysl, int* lbIdx, beliefPenaltyMPC_FLOAT *lubysu, int* ubIdx, beliefPenaltyMPC_FLOAT *Phi)


{
	int i;
	
	/* copy  H into PHI */
	for( i=0; i<62; i++ ){
		Phi[i] = H[i];
	}

	/* add llbysl onto Phi where necessary */
	for( i=0; i<62; i++ ){
		Phi[lbIdx[i]] += llbysl[i];
	}

	/* add lubysu onto Phi where necessary */
	for( i=0; i<22; i++){
		Phi[ubIdx[i]] +=  lubysu[i];
	}
	
	/* compute cholesky */
	for(i=0; i<62; i++)
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
 * where A is to be computed and is of size [20 x 62],
 * B is given and of size [20 x 62], L is a diagonal
 * matrix of size 20 stored in diagonal matrix 
 * storage format. Note the transpose of L has no impact!
 *
 * Result: A in column major storage format.
 *
 */
void beliefPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_20_62(beliefPenaltyMPC_FLOAT *L, beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *A)
{
    int i,j;
	 int k = 0;

	for( j=0; j<62; j++){
		for( i=0; i<20; i++){
			A[k] = B[k]/L[j];
			k++;
		}
	}

}


/**
 * Forward substitution for the matrix equation A*L' = B
 * where A is to be computed and is of size [20 x 62],
 * B is given and of size [20 x 62], L is a diagonal
 *  matrix of size 62 stored in diagonal 
 * storage format. Note the transpose of L!
 *
 * Result: A in diagonalzero storage format.
 *
 */
void beliefPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_20_62(beliefPenaltyMPC_FLOAT *L, beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *A)
{
	int j;
    for( j=0; j<62; j++ ){   
		A[j] = B[j]/L[j];
     }
}


/**
 * Compute C = A*B' where 
 *
 *	size(A) = [20 x 62]
 *  size(B) = [20 x 62] in diagzero format
 * 
 * A and C matrices are stored in column major format.
 * 
 * 
 */
void beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMTM_20_62_20(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *C)
{
    int i, j;
	
	for( i=0; i<20; i++ ){
		for( j=0; j<20; j++){
			C[j*20+i] = B[i*20+j]*A[i];
		}
	}

}


/**
 * Forward substitution to solve L*y = b where L is a
 * diagonal matrix in vector storage format.
 * 
 * The dimensions involved are 62.
 */
void beliefPenaltyMPC_LA_DIAG_FORWARDSUB_62(beliefPenaltyMPC_FLOAT *L, beliefPenaltyMPC_FLOAT *b, beliefPenaltyMPC_FLOAT *y)
{
    int i;

    for( i=0; i<62; i++ ){
		y[i] = b[i]/L[i];
    }
}


/*
 * Special function to compute the diagonal cholesky factorization of the 
 * positive definite augmented Hessian for block size 20.
 *
 * Inputs: - H = diagonal cost Hessian in diagonal storage format
 *         - llbysl = L / S of lower bounds
 *         - lubysu = L / S of upper bounds
 *
 * Output: Phi = sqrt(H + diag(llbysl) + diag(lubysu))
 * where Phi is stored in diagonal storage format
 */
void beliefPenaltyMPC_LA_DIAG_CHOL_ONELOOP_LBUB_20_20_20(beliefPenaltyMPC_FLOAT *H, beliefPenaltyMPC_FLOAT *llbysl, int* lbIdx, beliefPenaltyMPC_FLOAT *lubysu, int* ubIdx, beliefPenaltyMPC_FLOAT *Phi)


{
	int i;
	
	/* compute cholesky */
	for( i=0; i<20; i++ ){
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
 * where A is to be computed and is of size [20 x 20],
 * B is given and of size [20 x 20], L is a diagonal
 *  matrix of size 20 stored in diagonal 
 * storage format. Note the transpose of L!
 *
 * Result: A in diagonalzero storage format.
 *
 */
void beliefPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_20_20(beliefPenaltyMPC_FLOAT *L, beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *A)
{
	int j;
    for( j=0; j<20; j++ ){   
		A[j] = B[j]/L[j];
     }
}


/**
 * Forward substitution to solve L*y = b where L is a
 * diagonal matrix in vector storage format.
 * 
 * The dimensions involved are 20.
 */
void beliefPenaltyMPC_LA_DIAG_FORWARDSUB_20(beliefPenaltyMPC_FLOAT *L, beliefPenaltyMPC_FLOAT *b, beliefPenaltyMPC_FLOAT *y)
{
    int i;

    for( i=0; i<20; i++ ){
		y[i] = b[i]/L[i];
    }
}


/**
 * Compute L = B*B', where L is lower triangular of size NXp1
 * and B is a diagzero matrix of size [20 x 62] in column
 * storage format.
 * 
 */
void beliefPenaltyMPC_LA_DIAGZERO_MMT_20(beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *L)
{
    int i, ii, di;
    
    ii = 0; di = 0;
    for( i=0; i<20; i++ ){        
		L[ii+i] = B[i]*B[i];
        ii += ++di;
    }
}


/* 
 * Computes r = b - B*u
 * B is stored in diagzero format
 */
void beliefPenaltyMPC_LA_DIAGZERO_MVMSUB7_20(beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *b, beliefPenaltyMPC_FLOAT *r)
{
	int i;

	for( i=0; i<20; i++ ){
		r[i] = b[i] - B[i]*u[i];
	}	
	
}


/**
 * Compute L = A*A' + B*B', where L is lower triangular of size NXp1
 * and A is a dense matrix of size [20 x 62] in column
 * storage format, and B is of size [20 x 62] diagonalzero
 * storage format.
 * 
 * THIS ONE HAS THE WORST ACCES PATTERN POSSIBLE. 
 * POSSIBKE FIX: PUT A AND B INTO ROW MAJOR FORMAT FIRST.
 * 
 */
void beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMT2_20_62_62(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *L)
{
    int i, j, k, ii, di;
    beliefPenaltyMPC_FLOAT ltemp;
    
    ii = 0; di = 0;
    for( i=0; i<20; i++ ){        
        for( j=0; j<=i; j++ ){
            ltemp = 0; 
            for( k=0; k<62; k++ ){
                ltemp += A[k*20+i]*A[k*20+j];
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
void beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_20_62_62(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *b, beliefPenaltyMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;

	for( i=0; i<20; i++ ){
		r[i] = b[i] - A[k++]*x[0] - B[i]*u[i];
	}	

	for( j=1; j<62; j++ ){		
		for( i=0; i<20; i++ ){
			r[i] -= A[k++]*x[j];
		}
	}
	
}


/**
 * Compute L = A*A' + B*B', where L is lower triangular of size NXp1
 * and A is a dense matrix of size [20 x 62] in column
 * storage format, and B is of size [20 x 20] diagonalzero
 * storage format.
 * 
 * THIS ONE HAS THE WORST ACCES PATTERN POSSIBLE. 
 * POSSIBKE FIX: PUT A AND B INTO ROW MAJOR FORMAT FIRST.
 * 
 */
void beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMT2_20_62_20(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *L)
{
    int i, j, k, ii, di;
    beliefPenaltyMPC_FLOAT ltemp;
    
    ii = 0; di = 0;
    for( i=0; i<20; i++ ){        
        for( j=0; j<=i; j++ ){
            ltemp = 0; 
            for( k=0; k<62; k++ ){
                ltemp += A[k*20+i]*A[k*20+j];
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
void beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_20_62_20(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *b, beliefPenaltyMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;

	for( i=0; i<20; i++ ){
		r[i] = b[i] - A[k++]*x[0] - B[i]*u[i];
	}	

	for( j=1; j<62; j++ ){		
		for( i=0; i<20; i++ ){
			r[i] -= A[k++]*x[j];
		}
	}
	
}


/**
 * Cholesky factorization as above, but working on a matrix in 
 * lower triangular storage format of size 20 and outputting
 * the Cholesky factor to matrix L in lower triangular format.
 */
void beliefPenaltyMPC_LA_DENSE_CHOL_20(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *L)
{
    int i, j, k, di, dj;
	 int ii, jj;

    beliefPenaltyMPC_FLOAT l;
    beliefPenaltyMPC_FLOAT Mii;

	/* copy A to L first and then operate on L */
	/* COULD BE OPTIMIZED */
	ii=0; di=0;
	for( i=0; i<20; i++ ){
		for( j=0; j<=i; j++ ){
			L[ii+j] = A[ii+j];
		}
		ii += ++di;
	}    
	
	/* factor L */
	ii=0; di=0;
    for( i=0; i<20; i++ ){
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
        for( j=i+1; j<20; j++ ){
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
 * The dimensions involved are 20.
 */
void beliefPenaltyMPC_LA_DENSE_FORWARDSUB_20(beliefPenaltyMPC_FLOAT *L, beliefPenaltyMPC_FLOAT *b, beliefPenaltyMPC_FLOAT *y)
{
    int i,j,ii,di;
    beliefPenaltyMPC_FLOAT yel;
            
    ii = 0; di = 0;
    for( i=0; i<20; i++ ){
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
 * where A is to be computed and is of size [20 x 20],
 * B is given and of size [20 x 20], L is a lower tri-
 * angular matrix of size 20 stored in lower triangular 
 * storage format. Note the transpose of L AND B!
 *
 * Result: A in column major storage format.
 *
 */
void beliefPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_20_20(beliefPenaltyMPC_FLOAT *L, beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *A)
{
    int i,j,k,ii,di;
    beliefPenaltyMPC_FLOAT a;
    
    ii=0; di=0;
    for( j=0; j<20; j++ ){        
        for( i=0; i<20; i++ ){
            a = B[i*20+j];
            for( k=0; k<j; k++ ){
                a -= A[k*20+i]*L[ii+k];
            }    

			/* saturate for numerical stability */
			a = MIN(a, BIGM);
			a = MAX(a, -BIGM); 

			A[j*20+i] = a/L[ii+j];			
        }
        ii += ++di;
    }
}


/**
 * Compute L = L - A*A', where L is lower triangular of size 20
 * and A is a dense matrix of size [20 x 20] in column
 * storage format.
 * 
 * THIS ONE HAS THE WORST ACCES PATTERN POSSIBLE. 
 * POSSIBKE FIX: PUT A INTO ROW MAJOR FORMAT FIRST.
 * 
 */
void beliefPenaltyMPC_LA_DENSE_MMTSUB_20_20(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *L)
{
    int i, j, k, ii, di;
    beliefPenaltyMPC_FLOAT ltemp;
    
    ii = 0; di = 0;
    for( i=0; i<20; i++ ){        
        for( j=0; j<=i; j++ ){
            ltemp = 0; 
            for( k=0; k<20; k++ ){
                ltemp += A[k*20+i]*A[k*20+j];
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
void beliefPenaltyMPC_LA_DENSE_MVMSUB1_20_20(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *b, beliefPenaltyMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;

	for( i=0; i<20; i++ ){
		r[i] = b[i] - A[k++]*x[0];
	}	
	for( j=1; j<20; j++ ){		
		for( i=0; i<20; i++ ){
			r[i] -= A[k++]*x[j];
		}
	}
}


/**
 * Backward Substitution to solve L^T*x = y where L is a
 * lower triangular matrix in triangular storage format.
 * 
 * All involved dimensions are 20.
 */
void beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_20(beliefPenaltyMPC_FLOAT *L, beliefPenaltyMPC_FLOAT *y, beliefPenaltyMPC_FLOAT *x)
{
    int i, ii, di, j, jj, dj;
    beliefPenaltyMPC_FLOAT xel;    
	int start = 190;
    
    /* now solve L^T*x = y by backward substitution */
    ii = start; di = 19;
    for( i=19; i>=0; i-- ){        
        xel = y[i];        
        jj = start; dj = 19;
        for( j=19; j>i; j-- ){
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
 * Matrix vector multiplication y = b - M'*x where M is of size [20 x 20]
 * and stored in column major format. Note the transpose of M!
 */
void beliefPenaltyMPC_LA_DENSE_MTVMSUB_20_20(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *b, beliefPenaltyMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0; 
	for( i=0; i<20; i++ ){
		r[i] = b[i];
		for( j=0; j<20; j++ ){
			r[i] -= A[k++]*x[j];
		}
	}
}


/*
 * Vector subtraction z = -x - y for vectors of length 578.
 */
void beliefPenaltyMPC_LA_VSUB2_578(beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *y, beliefPenaltyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<578; i++){
		z[i] = -x[i] - y[i];
	}
}


/**
 * Forward-Backward-Substitution to solve L*L^T*x = b where L is a
 * diagonal matrix of size 62 in vector
 * storage format.
 */
void beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_62(beliefPenaltyMPC_FLOAT *L, beliefPenaltyMPC_FLOAT *b, beliefPenaltyMPC_FLOAT *x)
{
    int i;
            
    /* solve Ly = b by forward and backward substitution */
    for( i=0; i<62; i++ ){
		x[i] = b[i]/(L[i]*L[i]);
    }
    
}


/**
 * Forward-Backward-Substitution to solve L*L^T*x = b where L is a
 * diagonal matrix of size 20 in vector
 * storage format.
 */
void beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_20(beliefPenaltyMPC_FLOAT *L, beliefPenaltyMPC_FLOAT *b, beliefPenaltyMPC_FLOAT *x)
{
    int i;
            
    /* solve Ly = b by forward and backward substitution */
    for( i=0; i<20; i++ ){
		x[i] = b[i]/(L[i]*L[i]);
    }
    
}


/*
 * Vector subtraction z = x(xidx) - y where y, z and xidx are of length 62,
 * and x has length 62 and is indexed through yidx.
 */
void beliefPenaltyMPC_LA_VSUB_INDEXED_62(beliefPenaltyMPC_FLOAT *x, int* xidx, beliefPenaltyMPC_FLOAT *y, beliefPenaltyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<62; i++){
		z[i] = x[xidx[i]] - y[i];
	}
}


/*
 * Vector subtraction x = -u.*v - w for vectors of length 62.
 */
void beliefPenaltyMPC_LA_VSUB3_62(beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *v, beliefPenaltyMPC_FLOAT *w, beliefPenaltyMPC_FLOAT *x)
{
	int i;
	for( i=0; i<62; i++){
		x[i] = -u[i]*v[i] - w[i];
	}
}


/*
 * Vector subtraction z = -x - y(yidx) where y is of length 62
 * and z, x and yidx are of length 22.
 */
void beliefPenaltyMPC_LA_VSUB2_INDEXED_22(beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *y, int* yidx, beliefPenaltyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<22; i++){
		z[i] = -x[i] - y[yidx[i]];
	}
}


/*
 * Vector subtraction x = -u.*v - w for vectors of length 22.
 */
void beliefPenaltyMPC_LA_VSUB3_22(beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *v, beliefPenaltyMPC_FLOAT *w, beliefPenaltyMPC_FLOAT *x)
{
	int i;
	for( i=0; i<22; i++){
		x[i] = -u[i]*v[i] - w[i];
	}
}


/*
 * Vector subtraction z = x(xidx) - y where y, z and xidx are of length 20,
 * and x has length 20 and is indexed through yidx.
 */
void beliefPenaltyMPC_LA_VSUB_INDEXED_20(beliefPenaltyMPC_FLOAT *x, int* xidx, beliefPenaltyMPC_FLOAT *y, beliefPenaltyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<20; i++){
		z[i] = x[xidx[i]] - y[i];
	}
}


/*
 * Vector subtraction x = -u.*v - w for vectors of length 20.
 */
void beliefPenaltyMPC_LA_VSUB3_20(beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *v, beliefPenaltyMPC_FLOAT *w, beliefPenaltyMPC_FLOAT *x)
{
	int i;
	for( i=0; i<20; i++){
		x[i] = -u[i]*v[i] - w[i];
	}
}


/*
 * Vector subtraction z = -x - y(yidx) where y is of length 20
 * and z, x and yidx are of length 20.
 */
void beliefPenaltyMPC_LA_VSUB2_INDEXED_20(beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *y, int* yidx, beliefPenaltyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<20; i++){
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
        for( i=0; i<796; i++ ){
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
        if( i == 796 ){
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
    *mu_aff = mymu / (beliefPenaltyMPC_FLOAT)796;
    return lsIt;
}


/*
 * Vector subtraction x = u.*v - a where a is a scalar
*  and x,u,v are vectors of length 796.
 */
void beliefPenaltyMPC_LA_VSUB5_796(beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *v, beliefPenaltyMPC_FLOAT a, beliefPenaltyMPC_FLOAT *x)
{
	int i;
	for( i=0; i<796; i++){
		x[i] = u[i]*v[i] - a;
	}
}


/*
 * Computes x=0; x(uidx) += u/su; x(vidx) -= v/sv where x is of length 62,
 * u, su, uidx are of length 22 and v, sv, vidx are of length 62.
 */
void beliefPenaltyMPC_LA_VSUB6_INDEXED_62_22_62(beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *su, int* uidx, beliefPenaltyMPC_FLOAT *v, beliefPenaltyMPC_FLOAT *sv, int* vidx, beliefPenaltyMPC_FLOAT *x)
{
	int i;
	for( i=0; i<62; i++ ){
		x[i] = 0;
	}
	for( i=0; i<22; i++){
		x[uidx[i]] += u[i]/su[i];
	}
	for( i=0; i<62; i++){
		x[vidx[i]] -= v[i]/sv[i];
	}
}


/* 
 * Computes r =  B*u
 * where B is stored in diagzero format
 */
void beliefPenaltyMPC_LA_DIAGZERO_MVM_20(beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *r)
{
	int i;

	for( i=0; i<20; i++ ){
		r[i] = B[i]*u[i];
	}	
	
}


/* 
 * Computes r = A*x + B*u
 * where A is stored in column major format
 * and B is stored in diagzero format
 */
void beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_20_62_62(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;

	for( i=0; i<20; i++ ){
		r[i] = A[k++]*x[0] + B[i]*u[i];
	}	

	for( j=1; j<62; j++ ){		
		for( i=0; i<20; i++ ){
			r[i] += A[k++]*x[j];
		}
	}
	
}


/*
 * Computes x=0; x(uidx) += u/su; x(vidx) -= v/sv where x is of length 20,
 * u, su, uidx are of length 20 and v, sv, vidx are of length 20.
 */
void beliefPenaltyMPC_LA_VSUB6_INDEXED_20_20_20(beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *su, int* uidx, beliefPenaltyMPC_FLOAT *v, beliefPenaltyMPC_FLOAT *sv, int* vidx, beliefPenaltyMPC_FLOAT *x)
{
	int i;
	for( i=0; i<20; i++ ){
		x[i] = 0;
	}
	for( i=0; i<20; i++){
		x[uidx[i]] += u[i]/su[i];
	}
	for( i=0; i<20; i++){
		x[vidx[i]] -= v[i]/sv[i];
	}
}


/* 
 * Computes r = A*x + B*u
 * where A is stored in column major format
 * and B is stored in diagzero format
 */
void beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_20_62_20(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;

	for( i=0; i<20; i++ ){
		r[i] = A[k++]*x[0] + B[i]*u[i];
	}	

	for( j=1; j<62; j++ ){		
		for( i=0; i<20; i++ ){
			r[i] += A[k++]*x[j];
		}
	}
	
}


/*
 * Vector subtraction z = x - y for vectors of length 578.
 */
void beliefPenaltyMPC_LA_VSUB_578(beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *y, beliefPenaltyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<578; i++){
		z[i] = x[i] - y[i];
	}
}


/** 
 * Computes z = -r./s - u.*y(y)
 * where all vectors except of y are of length 62 (length of y >= 62).
 */
void beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_62(beliefPenaltyMPC_FLOAT *r, beliefPenaltyMPC_FLOAT *s, beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *y, int* yidx, beliefPenaltyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<62; i++ ){
		z[i] = -r[i]/s[i] - u[i]*y[yidx[i]];
	}
}


/** 
 * Computes z = -r./s + u.*y(y)
 * where all vectors except of y are of length 22 (length of y >= 22).
 */
void beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_22(beliefPenaltyMPC_FLOAT *r, beliefPenaltyMPC_FLOAT *s, beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *y, int* yidx, beliefPenaltyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<22; i++ ){
		z[i] = -r[i]/s[i] + u[i]*y[yidx[i]];
	}
}


/** 
 * Computes z = -r./s - u.*y(y)
 * where all vectors except of y are of length 20 (length of y >= 20).
 */
void beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_20(beliefPenaltyMPC_FLOAT *r, beliefPenaltyMPC_FLOAT *s, beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *y, int* yidx, beliefPenaltyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<20; i++ ){
		z[i] = -r[i]/s[i] - u[i]*y[yidx[i]];
	}
}


/** 
 * Computes z = -r./s + u.*y(y)
 * where all vectors except of y are of length 20 (length of y >= 20).
 */
void beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_20(beliefPenaltyMPC_FLOAT *r, beliefPenaltyMPC_FLOAT *s, beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *y, int* yidx, beliefPenaltyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<20; i++ ){
		z[i] = -r[i]/s[i] + u[i]*y[yidx[i]];
	}
}


/*
 * Computes ds = -l.\(r + s.*dl) for vectors of length 796.
 */
void beliefPenaltyMPC_LA_VSUB7_796(beliefPenaltyMPC_FLOAT *l, beliefPenaltyMPC_FLOAT *r, beliefPenaltyMPC_FLOAT *s, beliefPenaltyMPC_FLOAT *dl, beliefPenaltyMPC_FLOAT *ds)
{
	int i;
	for( i=0; i<796; i++){
		ds[i] = -(r[i] + s[i]*dl[i])/l[i];
	}
}


/*
 * Vector addition x = x + y for vectors of length 578.
 */
void beliefPenaltyMPC_LA_VADD_578(beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *y)
{
	int i;
	for( i=0; i<578; i++){
		x[i] += y[i];
	}
}


/*
 * Vector addition x = x + y for vectors of length 200.
 */
void beliefPenaltyMPC_LA_VADD_200(beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *y)
{
	int i;
	for( i=0; i<200; i++){
		x[i] += y[i];
	}
}


/*
 * Vector addition x = x + y for vectors of length 796.
 */
void beliefPenaltyMPC_LA_VADD_796(beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *y)
{
	int i;
	for( i=0; i<796; i++){
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
        for( i=0; i<796; i++ ){
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
        if( i == 796 ){
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
    for( i=0; i<578; i++ ){
        z[i] += a_gamma*dz[i];
    }
    
    /* equality constraint multipliers */
    for( i=0; i<200; i++ ){
        v[i] += a_gamma*dv[i];
    }
    
    /* inequality constraint multipliers & slacks, also update mu */
    *mu = 0;
    for( i=0; i<796; i++ ){
        dltemp = l[i] + a_gamma*dl[i]; l[i] = dltemp;
        dstemp = s[i] + a_gamma*ds[i]; s[i] = dstemp;
        *mu += dltemp*dstemp;
    }
    
    *a = a_gamma;
    *mu /= (beliefPenaltyMPC_FLOAT)796;
    return lsIt;
}




/* VARIABLE DEFINITIONS ------------------------------------------------ */
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_z[578];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_v[200];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_dz_aff[578];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_dv_aff[200];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_grad_cost[578];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_grad_eq[578];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_rd[578];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_l[796];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_s[796];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_lbys[796];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_dl_aff[796];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ds_aff[796];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_dz_cc[578];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_dv_cc[200];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_dl_cc[796];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ds_cc[796];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ccrhs[796];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_grad_ineq[578];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_H0[62] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 2.0000000000000000E+001, 2.0000000000000000E+001, 2.0000000000000000E+001, 2.0000000000000000E+001, 2.0000000000000000E+001, 2.0000000000000000E+001, 2.0000000000000000E+001, 2.0000000000000000E+001, 2.0000000000000000E+001, 2.0000000000000000E+001, 2.0000000000000000E+001, 2.0000000000000000E+001, 2.0000000000000000E+001, 2.0000000000000000E+001, 2.0000000000000000E+001, 2.0000000000000000E+000, 2.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_z0 = beliefPenaltyMPC_z + 0;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzaff0 = beliefPenaltyMPC_dz_aff + 0;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzcc0 = beliefPenaltyMPC_dz_cc + 0;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_rd0 = beliefPenaltyMPC_rd + 0;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lbyrd0[62];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_cost0 = beliefPenaltyMPC_grad_cost + 0;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_eq0 = beliefPenaltyMPC_grad_eq + 0;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_ineq0 = beliefPenaltyMPC_grad_ineq + 0;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ctv0[62];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_v0 = beliefPenaltyMPC_v + 0;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_re0[20];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_beta0[20];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_betacc0[20];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvaff0 = beliefPenaltyMPC_dv_aff + 0;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvcc0 = beliefPenaltyMPC_dv_cc + 0;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_V0[1240];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Yd0[210];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ld0[210];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_yy0[20];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_bmy0[20];
int beliefPenaltyMPC_lbIdx0[62] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llb0 = beliefPenaltyMPC_l + 0;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_slb0 = beliefPenaltyMPC_s + 0;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llbbyslb0 = beliefPenaltyMPC_lbys + 0;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_rilb0[62];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbaff0 = beliefPenaltyMPC_dl_aff + 0;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbaff0 = beliefPenaltyMPC_ds_aff + 0;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbcc0 = beliefPenaltyMPC_dl_cc + 0;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbcc0 = beliefPenaltyMPC_ds_cc + 0;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsl0 = beliefPenaltyMPC_ccrhs + 0;
int beliefPenaltyMPC_ubIdx0[22] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lub0 = beliefPenaltyMPC_l + 62;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_sub0 = beliefPenaltyMPC_s + 62;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lubbysub0 = beliefPenaltyMPC_lbys + 62;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_riub0[22];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubaff0 = beliefPenaltyMPC_dl_aff + 62;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubaff0 = beliefPenaltyMPC_ds_aff + 62;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubcc0 = beliefPenaltyMPC_dl_cc + 62;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubcc0 = beliefPenaltyMPC_ds_cc + 62;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsub0 = beliefPenaltyMPC_ccrhs + 62;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Phi0[62];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_D0[62] = {1.0000000000000000E+000, 
1.0000000000000000E+000, 
1.0000000000000000E+000, 
1.0000000000000000E+000, 
1.0000000000000000E+000, 
1.0000000000000000E+000, 
1.0000000000000000E+000, 
1.0000000000000000E+000, 
1.0000000000000000E+000, 
1.0000000000000000E+000, 
1.0000000000000000E+000, 
1.0000000000000000E+000, 
1.0000000000000000E+000, 
1.0000000000000000E+000, 
1.0000000000000000E+000, 
1.0000000000000000E+000, 
1.0000000000000000E+000, 
1.0000000000000000E+000, 
1.0000000000000000E+000, 
1.0000000000000000E+000};
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_W0[62];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_z1 = beliefPenaltyMPC_z + 62;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzaff1 = beliefPenaltyMPC_dz_aff + 62;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzcc1 = beliefPenaltyMPC_dz_cc + 62;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_rd1 = beliefPenaltyMPC_rd + 62;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lbyrd1[62];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_cost1 = beliefPenaltyMPC_grad_cost + 62;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_eq1 = beliefPenaltyMPC_grad_eq + 62;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_ineq1 = beliefPenaltyMPC_grad_ineq + 62;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ctv1[62];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_v1 = beliefPenaltyMPC_v + 20;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_re1[20];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_beta1[20];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_betacc1[20];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvaff1 = beliefPenaltyMPC_dv_aff + 20;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvcc1 = beliefPenaltyMPC_dv_cc + 20;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_V1[1240];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Yd1[210];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ld1[210];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_yy1[20];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_bmy1[20];
int beliefPenaltyMPC_lbIdx1[62] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llb1 = beliefPenaltyMPC_l + 84;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_slb1 = beliefPenaltyMPC_s + 84;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llbbyslb1 = beliefPenaltyMPC_lbys + 84;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_rilb1[62];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbaff1 = beliefPenaltyMPC_dl_aff + 84;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbaff1 = beliefPenaltyMPC_ds_aff + 84;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbcc1 = beliefPenaltyMPC_dl_cc + 84;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbcc1 = beliefPenaltyMPC_ds_cc + 84;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsl1 = beliefPenaltyMPC_ccrhs + 84;
int beliefPenaltyMPC_ubIdx1[22] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lub1 = beliefPenaltyMPC_l + 146;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_sub1 = beliefPenaltyMPC_s + 146;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lubbysub1 = beliefPenaltyMPC_lbys + 146;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_riub1[22];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubaff1 = beliefPenaltyMPC_dl_aff + 146;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubaff1 = beliefPenaltyMPC_ds_aff + 146;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubcc1 = beliefPenaltyMPC_dl_cc + 146;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubcc1 = beliefPenaltyMPC_ds_cc + 146;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsub1 = beliefPenaltyMPC_ccrhs + 146;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Phi1[62];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_D1[62] = {-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000};
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_W1[62];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ysd1[400];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lsd1[400];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_z2 = beliefPenaltyMPC_z + 124;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzaff2 = beliefPenaltyMPC_dz_aff + 124;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzcc2 = beliefPenaltyMPC_dz_cc + 124;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_rd2 = beliefPenaltyMPC_rd + 124;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lbyrd2[62];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_cost2 = beliefPenaltyMPC_grad_cost + 124;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_eq2 = beliefPenaltyMPC_grad_eq + 124;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_ineq2 = beliefPenaltyMPC_grad_ineq + 124;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ctv2[62];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_v2 = beliefPenaltyMPC_v + 40;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_re2[20];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_beta2[20];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_betacc2[20];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvaff2 = beliefPenaltyMPC_dv_aff + 40;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvcc2 = beliefPenaltyMPC_dv_cc + 40;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_V2[1240];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Yd2[210];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ld2[210];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_yy2[20];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_bmy2[20];
int beliefPenaltyMPC_lbIdx2[62] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llb2 = beliefPenaltyMPC_l + 168;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_slb2 = beliefPenaltyMPC_s + 168;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llbbyslb2 = beliefPenaltyMPC_lbys + 168;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_rilb2[62];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbaff2 = beliefPenaltyMPC_dl_aff + 168;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbaff2 = beliefPenaltyMPC_ds_aff + 168;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbcc2 = beliefPenaltyMPC_dl_cc + 168;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbcc2 = beliefPenaltyMPC_ds_cc + 168;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsl2 = beliefPenaltyMPC_ccrhs + 168;
int beliefPenaltyMPC_ubIdx2[22] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lub2 = beliefPenaltyMPC_l + 230;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_sub2 = beliefPenaltyMPC_s + 230;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lubbysub2 = beliefPenaltyMPC_lbys + 230;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_riub2[22];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubaff2 = beliefPenaltyMPC_dl_aff + 230;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubaff2 = beliefPenaltyMPC_ds_aff + 230;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubcc2 = beliefPenaltyMPC_dl_cc + 230;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubcc2 = beliefPenaltyMPC_ds_cc + 230;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsub2 = beliefPenaltyMPC_ccrhs + 230;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Phi2[62];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_W2[62];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ysd2[400];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lsd2[400];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_z3 = beliefPenaltyMPC_z + 186;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzaff3 = beliefPenaltyMPC_dz_aff + 186;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzcc3 = beliefPenaltyMPC_dz_cc + 186;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_rd3 = beliefPenaltyMPC_rd + 186;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lbyrd3[62];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_cost3 = beliefPenaltyMPC_grad_cost + 186;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_eq3 = beliefPenaltyMPC_grad_eq + 186;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_ineq3 = beliefPenaltyMPC_grad_ineq + 186;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ctv3[62];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_v3 = beliefPenaltyMPC_v + 60;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_re3[20];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_beta3[20];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_betacc3[20];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvaff3 = beliefPenaltyMPC_dv_aff + 60;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvcc3 = beliefPenaltyMPC_dv_cc + 60;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_V3[1240];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Yd3[210];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ld3[210];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_yy3[20];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_bmy3[20];
int beliefPenaltyMPC_lbIdx3[62] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llb3 = beliefPenaltyMPC_l + 252;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_slb3 = beliefPenaltyMPC_s + 252;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llbbyslb3 = beliefPenaltyMPC_lbys + 252;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_rilb3[62];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbaff3 = beliefPenaltyMPC_dl_aff + 252;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbaff3 = beliefPenaltyMPC_ds_aff + 252;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbcc3 = beliefPenaltyMPC_dl_cc + 252;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbcc3 = beliefPenaltyMPC_ds_cc + 252;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsl3 = beliefPenaltyMPC_ccrhs + 252;
int beliefPenaltyMPC_ubIdx3[22] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lub3 = beliefPenaltyMPC_l + 314;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_sub3 = beliefPenaltyMPC_s + 314;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lubbysub3 = beliefPenaltyMPC_lbys + 314;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_riub3[22];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubaff3 = beliefPenaltyMPC_dl_aff + 314;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubaff3 = beliefPenaltyMPC_ds_aff + 314;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubcc3 = beliefPenaltyMPC_dl_cc + 314;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubcc3 = beliefPenaltyMPC_ds_cc + 314;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsub3 = beliefPenaltyMPC_ccrhs + 314;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Phi3[62];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_W3[62];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ysd3[400];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lsd3[400];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_z4 = beliefPenaltyMPC_z + 248;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzaff4 = beliefPenaltyMPC_dz_aff + 248;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzcc4 = beliefPenaltyMPC_dz_cc + 248;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_rd4 = beliefPenaltyMPC_rd + 248;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lbyrd4[62];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_cost4 = beliefPenaltyMPC_grad_cost + 248;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_eq4 = beliefPenaltyMPC_grad_eq + 248;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_ineq4 = beliefPenaltyMPC_grad_ineq + 248;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ctv4[62];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_v4 = beliefPenaltyMPC_v + 80;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_re4[20];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_beta4[20];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_betacc4[20];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvaff4 = beliefPenaltyMPC_dv_aff + 80;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvcc4 = beliefPenaltyMPC_dv_cc + 80;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_V4[1240];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Yd4[210];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ld4[210];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_yy4[20];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_bmy4[20];
int beliefPenaltyMPC_lbIdx4[62] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llb4 = beliefPenaltyMPC_l + 336;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_slb4 = beliefPenaltyMPC_s + 336;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llbbyslb4 = beliefPenaltyMPC_lbys + 336;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_rilb4[62];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbaff4 = beliefPenaltyMPC_dl_aff + 336;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbaff4 = beliefPenaltyMPC_ds_aff + 336;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbcc4 = beliefPenaltyMPC_dl_cc + 336;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbcc4 = beliefPenaltyMPC_ds_cc + 336;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsl4 = beliefPenaltyMPC_ccrhs + 336;
int beliefPenaltyMPC_ubIdx4[22] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lub4 = beliefPenaltyMPC_l + 398;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_sub4 = beliefPenaltyMPC_s + 398;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lubbysub4 = beliefPenaltyMPC_lbys + 398;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_riub4[22];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubaff4 = beliefPenaltyMPC_dl_aff + 398;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubaff4 = beliefPenaltyMPC_ds_aff + 398;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubcc4 = beliefPenaltyMPC_dl_cc + 398;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubcc4 = beliefPenaltyMPC_ds_cc + 398;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsub4 = beliefPenaltyMPC_ccrhs + 398;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Phi4[62];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_W4[62];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ysd4[400];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lsd4[400];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_z5 = beliefPenaltyMPC_z + 310;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzaff5 = beliefPenaltyMPC_dz_aff + 310;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzcc5 = beliefPenaltyMPC_dz_cc + 310;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_rd5 = beliefPenaltyMPC_rd + 310;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lbyrd5[62];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_cost5 = beliefPenaltyMPC_grad_cost + 310;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_eq5 = beliefPenaltyMPC_grad_eq + 310;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_ineq5 = beliefPenaltyMPC_grad_ineq + 310;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ctv5[62];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_v5 = beliefPenaltyMPC_v + 100;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_re5[20];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_beta5[20];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_betacc5[20];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvaff5 = beliefPenaltyMPC_dv_aff + 100;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvcc5 = beliefPenaltyMPC_dv_cc + 100;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_V5[1240];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Yd5[210];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ld5[210];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_yy5[20];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_bmy5[20];
int beliefPenaltyMPC_lbIdx5[62] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llb5 = beliefPenaltyMPC_l + 420;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_slb5 = beliefPenaltyMPC_s + 420;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llbbyslb5 = beliefPenaltyMPC_lbys + 420;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_rilb5[62];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbaff5 = beliefPenaltyMPC_dl_aff + 420;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbaff5 = beliefPenaltyMPC_ds_aff + 420;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbcc5 = beliefPenaltyMPC_dl_cc + 420;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbcc5 = beliefPenaltyMPC_ds_cc + 420;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsl5 = beliefPenaltyMPC_ccrhs + 420;
int beliefPenaltyMPC_ubIdx5[22] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lub5 = beliefPenaltyMPC_l + 482;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_sub5 = beliefPenaltyMPC_s + 482;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lubbysub5 = beliefPenaltyMPC_lbys + 482;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_riub5[22];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubaff5 = beliefPenaltyMPC_dl_aff + 482;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubaff5 = beliefPenaltyMPC_ds_aff + 482;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubcc5 = beliefPenaltyMPC_dl_cc + 482;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubcc5 = beliefPenaltyMPC_ds_cc + 482;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsub5 = beliefPenaltyMPC_ccrhs + 482;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Phi5[62];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_W5[62];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ysd5[400];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lsd5[400];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_z6 = beliefPenaltyMPC_z + 372;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzaff6 = beliefPenaltyMPC_dz_aff + 372;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzcc6 = beliefPenaltyMPC_dz_cc + 372;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_rd6 = beliefPenaltyMPC_rd + 372;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lbyrd6[62];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_cost6 = beliefPenaltyMPC_grad_cost + 372;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_eq6 = beliefPenaltyMPC_grad_eq + 372;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_ineq6 = beliefPenaltyMPC_grad_ineq + 372;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ctv6[62];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_v6 = beliefPenaltyMPC_v + 120;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_re6[20];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_beta6[20];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_betacc6[20];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvaff6 = beliefPenaltyMPC_dv_aff + 120;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvcc6 = beliefPenaltyMPC_dv_cc + 120;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_V6[1240];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Yd6[210];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ld6[210];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_yy6[20];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_bmy6[20];
int beliefPenaltyMPC_lbIdx6[62] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llb6 = beliefPenaltyMPC_l + 504;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_slb6 = beliefPenaltyMPC_s + 504;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llbbyslb6 = beliefPenaltyMPC_lbys + 504;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_rilb6[62];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbaff6 = beliefPenaltyMPC_dl_aff + 504;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbaff6 = beliefPenaltyMPC_ds_aff + 504;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbcc6 = beliefPenaltyMPC_dl_cc + 504;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbcc6 = beliefPenaltyMPC_ds_cc + 504;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsl6 = beliefPenaltyMPC_ccrhs + 504;
int beliefPenaltyMPC_ubIdx6[22] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lub6 = beliefPenaltyMPC_l + 566;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_sub6 = beliefPenaltyMPC_s + 566;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lubbysub6 = beliefPenaltyMPC_lbys + 566;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_riub6[22];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubaff6 = beliefPenaltyMPC_dl_aff + 566;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubaff6 = beliefPenaltyMPC_ds_aff + 566;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubcc6 = beliefPenaltyMPC_dl_cc + 566;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubcc6 = beliefPenaltyMPC_ds_cc + 566;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsub6 = beliefPenaltyMPC_ccrhs + 566;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Phi6[62];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_W6[62];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ysd6[400];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lsd6[400];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_z7 = beliefPenaltyMPC_z + 434;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzaff7 = beliefPenaltyMPC_dz_aff + 434;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzcc7 = beliefPenaltyMPC_dz_cc + 434;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_rd7 = beliefPenaltyMPC_rd + 434;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lbyrd7[62];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_cost7 = beliefPenaltyMPC_grad_cost + 434;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_eq7 = beliefPenaltyMPC_grad_eq + 434;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_ineq7 = beliefPenaltyMPC_grad_ineq + 434;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ctv7[62];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_v7 = beliefPenaltyMPC_v + 140;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_re7[20];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_beta7[20];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_betacc7[20];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvaff7 = beliefPenaltyMPC_dv_aff + 140;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvcc7 = beliefPenaltyMPC_dv_cc + 140;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_V7[1240];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Yd7[210];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ld7[210];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_yy7[20];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_bmy7[20];
int beliefPenaltyMPC_lbIdx7[62] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llb7 = beliefPenaltyMPC_l + 588;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_slb7 = beliefPenaltyMPC_s + 588;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llbbyslb7 = beliefPenaltyMPC_lbys + 588;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_rilb7[62];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbaff7 = beliefPenaltyMPC_dl_aff + 588;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbaff7 = beliefPenaltyMPC_ds_aff + 588;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbcc7 = beliefPenaltyMPC_dl_cc + 588;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbcc7 = beliefPenaltyMPC_ds_cc + 588;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsl7 = beliefPenaltyMPC_ccrhs + 588;
int beliefPenaltyMPC_ubIdx7[22] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lub7 = beliefPenaltyMPC_l + 650;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_sub7 = beliefPenaltyMPC_s + 650;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lubbysub7 = beliefPenaltyMPC_lbys + 650;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_riub7[22];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubaff7 = beliefPenaltyMPC_dl_aff + 650;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubaff7 = beliefPenaltyMPC_ds_aff + 650;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubcc7 = beliefPenaltyMPC_dl_cc + 650;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubcc7 = beliefPenaltyMPC_ds_cc + 650;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsub7 = beliefPenaltyMPC_ccrhs + 650;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Phi7[62];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_W7[62];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ysd7[400];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lsd7[400];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_z8 = beliefPenaltyMPC_z + 496;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzaff8 = beliefPenaltyMPC_dz_aff + 496;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzcc8 = beliefPenaltyMPC_dz_cc + 496;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_rd8 = beliefPenaltyMPC_rd + 496;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lbyrd8[62];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_cost8 = beliefPenaltyMPC_grad_cost + 496;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_eq8 = beliefPenaltyMPC_grad_eq + 496;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_ineq8 = beliefPenaltyMPC_grad_ineq + 496;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ctv8[62];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_v8 = beliefPenaltyMPC_v + 160;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_re8[20];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_beta8[20];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_betacc8[20];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvaff8 = beliefPenaltyMPC_dv_aff + 160;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvcc8 = beliefPenaltyMPC_dv_cc + 160;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_V8[1240];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Yd8[210];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ld8[210];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_yy8[20];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_bmy8[20];
int beliefPenaltyMPC_lbIdx8[62] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llb8 = beliefPenaltyMPC_l + 672;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_slb8 = beliefPenaltyMPC_s + 672;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llbbyslb8 = beliefPenaltyMPC_lbys + 672;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_rilb8[62];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbaff8 = beliefPenaltyMPC_dl_aff + 672;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbaff8 = beliefPenaltyMPC_ds_aff + 672;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbcc8 = beliefPenaltyMPC_dl_cc + 672;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbcc8 = beliefPenaltyMPC_ds_cc + 672;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsl8 = beliefPenaltyMPC_ccrhs + 672;
int beliefPenaltyMPC_ubIdx8[22] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lub8 = beliefPenaltyMPC_l + 734;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_sub8 = beliefPenaltyMPC_s + 734;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lubbysub8 = beliefPenaltyMPC_lbys + 734;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_riub8[22];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubaff8 = beliefPenaltyMPC_dl_aff + 734;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubaff8 = beliefPenaltyMPC_ds_aff + 734;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubcc8 = beliefPenaltyMPC_dl_cc + 734;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubcc8 = beliefPenaltyMPC_ds_cc + 734;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsub8 = beliefPenaltyMPC_ccrhs + 734;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Phi8[62];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_W8[62];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ysd8[400];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lsd8[400];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_H9[20] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 2.0000000000000000E+001, 2.0000000000000000E+001, 2.0000000000000000E+001, 2.0000000000000000E+001, 2.0000000000000000E+001, 2.0000000000000000E+001, 2.0000000000000000E+001, 2.0000000000000000E+001, 2.0000000000000000E+001, 2.0000000000000000E+001, 2.0000000000000000E+001, 2.0000000000000000E+001, 2.0000000000000000E+001, 2.0000000000000000E+001, 2.0000000000000000E+001};
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_f9[20] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_z9 = beliefPenaltyMPC_z + 558;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzaff9 = beliefPenaltyMPC_dz_aff + 558;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzcc9 = beliefPenaltyMPC_dz_cc + 558;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_rd9 = beliefPenaltyMPC_rd + 558;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lbyrd9[20];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_cost9 = beliefPenaltyMPC_grad_cost + 558;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_eq9 = beliefPenaltyMPC_grad_eq + 558;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_ineq9 = beliefPenaltyMPC_grad_ineq + 558;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ctv9[20];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_v9 = beliefPenaltyMPC_v + 180;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_re9[20];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_beta9[20];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_betacc9[20];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvaff9 = beliefPenaltyMPC_dv_aff + 180;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvcc9 = beliefPenaltyMPC_dv_cc + 180;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_V9[400];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Yd9[210];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ld9[210];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_yy9[20];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_bmy9[20];
int beliefPenaltyMPC_lbIdx9[20] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llb9 = beliefPenaltyMPC_l + 756;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_slb9 = beliefPenaltyMPC_s + 756;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llbbyslb9 = beliefPenaltyMPC_lbys + 756;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_rilb9[20];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbaff9 = beliefPenaltyMPC_dl_aff + 756;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbaff9 = beliefPenaltyMPC_ds_aff + 756;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbcc9 = beliefPenaltyMPC_dl_cc + 756;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbcc9 = beliefPenaltyMPC_ds_cc + 756;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsl9 = beliefPenaltyMPC_ccrhs + 756;
int beliefPenaltyMPC_ubIdx9[20] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lub9 = beliefPenaltyMPC_l + 776;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_sub9 = beliefPenaltyMPC_s + 776;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lubbysub9 = beliefPenaltyMPC_lbys + 776;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_riub9[20];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubaff9 = beliefPenaltyMPC_dl_aff + 776;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubaff9 = beliefPenaltyMPC_ds_aff + 776;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubcc9 = beliefPenaltyMPC_dl_cc + 776;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubcc9 = beliefPenaltyMPC_ds_cc + 776;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsub9 = beliefPenaltyMPC_ccrhs + 776;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Phi9[20];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_D9[20] = {-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000};
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_W9[20];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ysd9[400];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lsd9[400];
beliefPenaltyMPC_FLOAT musigma;
beliefPenaltyMPC_FLOAT sigma_3rdroot;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Diag1_0[62];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Diag2_0[62];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_L_0[1891];




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
beliefPenaltyMPC_LA_INITIALIZEVECTOR_578(beliefPenaltyMPC_z, 0);
beliefPenaltyMPC_LA_INITIALIZEVECTOR_200(beliefPenaltyMPC_v, 1);
beliefPenaltyMPC_LA_INITIALIZEVECTOR_796(beliefPenaltyMPC_l, 1);
beliefPenaltyMPC_LA_INITIALIZEVECTOR_796(beliefPenaltyMPC_s, 1);
info->mu = 0;
beliefPenaltyMPC_LA_DOTACC_796(beliefPenaltyMPC_l, beliefPenaltyMPC_s, &info->mu);
info->mu /= 796;
while( 1 ){
info->pobj = 0;
beliefPenaltyMPC_LA_DIAG_QUADFCN_62(beliefPenaltyMPC_H0, params->f1, beliefPenaltyMPC_z0, beliefPenaltyMPC_grad_cost0, &info->pobj);
beliefPenaltyMPC_LA_DIAG_QUADFCN_62(beliefPenaltyMPC_H0, params->f2, beliefPenaltyMPC_z1, beliefPenaltyMPC_grad_cost1, &info->pobj);
beliefPenaltyMPC_LA_DIAG_QUADFCN_62(beliefPenaltyMPC_H0, params->f3, beliefPenaltyMPC_z2, beliefPenaltyMPC_grad_cost2, &info->pobj);
beliefPenaltyMPC_LA_DIAG_QUADFCN_62(beliefPenaltyMPC_H0, params->f4, beliefPenaltyMPC_z3, beliefPenaltyMPC_grad_cost3, &info->pobj);
beliefPenaltyMPC_LA_DIAG_QUADFCN_62(beliefPenaltyMPC_H0, params->f5, beliefPenaltyMPC_z4, beliefPenaltyMPC_grad_cost4, &info->pobj);
beliefPenaltyMPC_LA_DIAG_QUADFCN_62(beliefPenaltyMPC_H0, params->f6, beliefPenaltyMPC_z5, beliefPenaltyMPC_grad_cost5, &info->pobj);
beliefPenaltyMPC_LA_DIAG_QUADFCN_62(beliefPenaltyMPC_H0, params->f7, beliefPenaltyMPC_z6, beliefPenaltyMPC_grad_cost6, &info->pobj);
beliefPenaltyMPC_LA_DIAG_QUADFCN_62(beliefPenaltyMPC_H0, params->f8, beliefPenaltyMPC_z7, beliefPenaltyMPC_grad_cost7, &info->pobj);
beliefPenaltyMPC_LA_DIAG_QUADFCN_62(beliefPenaltyMPC_H0, params->f9, beliefPenaltyMPC_z8, beliefPenaltyMPC_grad_cost8, &info->pobj);
beliefPenaltyMPC_LA_DIAG_QUADFCN_20(beliefPenaltyMPC_H9, beliefPenaltyMPC_f9, beliefPenaltyMPC_z9, beliefPenaltyMPC_grad_cost9, &info->pobj);
info->res_eq = 0;
info->dgap = 0;
beliefPenaltyMPC_LA_DIAGZERO_MVMSUB6_20(beliefPenaltyMPC_D0, beliefPenaltyMPC_z0, params->e1, beliefPenaltyMPC_v0, beliefPenaltyMPC_re0, &info->dgap, &info->res_eq);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_20_62_62(params->C1, beliefPenaltyMPC_z0, beliefPenaltyMPC_D1, beliefPenaltyMPC_z1, params->e2, beliefPenaltyMPC_v1, beliefPenaltyMPC_re1, &info->dgap, &info->res_eq);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_20_62_62(params->C2, beliefPenaltyMPC_z1, beliefPenaltyMPC_D1, beliefPenaltyMPC_z2, params->e3, beliefPenaltyMPC_v2, beliefPenaltyMPC_re2, &info->dgap, &info->res_eq);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_20_62_62(params->C3, beliefPenaltyMPC_z2, beliefPenaltyMPC_D1, beliefPenaltyMPC_z3, params->e4, beliefPenaltyMPC_v3, beliefPenaltyMPC_re3, &info->dgap, &info->res_eq);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_20_62_62(params->C4, beliefPenaltyMPC_z3, beliefPenaltyMPC_D1, beliefPenaltyMPC_z4, params->e5, beliefPenaltyMPC_v4, beliefPenaltyMPC_re4, &info->dgap, &info->res_eq);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_20_62_62(params->C5, beliefPenaltyMPC_z4, beliefPenaltyMPC_D1, beliefPenaltyMPC_z5, params->e6, beliefPenaltyMPC_v5, beliefPenaltyMPC_re5, &info->dgap, &info->res_eq);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_20_62_62(params->C6, beliefPenaltyMPC_z5, beliefPenaltyMPC_D1, beliefPenaltyMPC_z6, params->e7, beliefPenaltyMPC_v6, beliefPenaltyMPC_re6, &info->dgap, &info->res_eq);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_20_62_62(params->C7, beliefPenaltyMPC_z6, beliefPenaltyMPC_D1, beliefPenaltyMPC_z7, params->e8, beliefPenaltyMPC_v7, beliefPenaltyMPC_re7, &info->dgap, &info->res_eq);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_20_62_62(params->C8, beliefPenaltyMPC_z7, beliefPenaltyMPC_D1, beliefPenaltyMPC_z8, params->e9, beliefPenaltyMPC_v8, beliefPenaltyMPC_re8, &info->dgap, &info->res_eq);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_20_62_20(params->C9, beliefPenaltyMPC_z8, beliefPenaltyMPC_D9, beliefPenaltyMPC_z9, params->e10, beliefPenaltyMPC_v9, beliefPenaltyMPC_re9, &info->dgap, &info->res_eq);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_20_62_20(params->C1, beliefPenaltyMPC_v1, beliefPenaltyMPC_D0, beliefPenaltyMPC_v0, beliefPenaltyMPC_grad_eq0);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_20_62_20(params->C2, beliefPenaltyMPC_v2, beliefPenaltyMPC_D1, beliefPenaltyMPC_v1, beliefPenaltyMPC_grad_eq1);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_20_62_20(params->C3, beliefPenaltyMPC_v3, beliefPenaltyMPC_D1, beliefPenaltyMPC_v2, beliefPenaltyMPC_grad_eq2);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_20_62_20(params->C4, beliefPenaltyMPC_v4, beliefPenaltyMPC_D1, beliefPenaltyMPC_v3, beliefPenaltyMPC_grad_eq3);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_20_62_20(params->C5, beliefPenaltyMPC_v5, beliefPenaltyMPC_D1, beliefPenaltyMPC_v4, beliefPenaltyMPC_grad_eq4);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_20_62_20(params->C6, beliefPenaltyMPC_v6, beliefPenaltyMPC_D1, beliefPenaltyMPC_v5, beliefPenaltyMPC_grad_eq5);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_20_62_20(params->C7, beliefPenaltyMPC_v7, beliefPenaltyMPC_D1, beliefPenaltyMPC_v6, beliefPenaltyMPC_grad_eq6);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_20_62_20(params->C8, beliefPenaltyMPC_v8, beliefPenaltyMPC_D1, beliefPenaltyMPC_v7, beliefPenaltyMPC_grad_eq7);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_20_62_20(params->C9, beliefPenaltyMPC_v9, beliefPenaltyMPC_D1, beliefPenaltyMPC_v8, beliefPenaltyMPC_grad_eq8);
beliefPenaltyMPC_LA_DIAGZERO_MTVM_20_20(beliefPenaltyMPC_D9, beliefPenaltyMPC_v9, beliefPenaltyMPC_grad_eq9);
info->res_ineq = 0;
beliefPenaltyMPC_LA_VSUBADD3_62(params->lb1, beliefPenaltyMPC_z0, beliefPenaltyMPC_lbIdx0, beliefPenaltyMPC_llb0, beliefPenaltyMPC_slb0, beliefPenaltyMPC_rilb0, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD2_22(beliefPenaltyMPC_z0, beliefPenaltyMPC_ubIdx0, params->ub1, beliefPenaltyMPC_lub0, beliefPenaltyMPC_sub0, beliefPenaltyMPC_riub0, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD3_62(params->lb2, beliefPenaltyMPC_z1, beliefPenaltyMPC_lbIdx1, beliefPenaltyMPC_llb1, beliefPenaltyMPC_slb1, beliefPenaltyMPC_rilb1, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD2_22(beliefPenaltyMPC_z1, beliefPenaltyMPC_ubIdx1, params->ub2, beliefPenaltyMPC_lub1, beliefPenaltyMPC_sub1, beliefPenaltyMPC_riub1, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD3_62(params->lb3, beliefPenaltyMPC_z2, beliefPenaltyMPC_lbIdx2, beliefPenaltyMPC_llb2, beliefPenaltyMPC_slb2, beliefPenaltyMPC_rilb2, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD2_22(beliefPenaltyMPC_z2, beliefPenaltyMPC_ubIdx2, params->ub3, beliefPenaltyMPC_lub2, beliefPenaltyMPC_sub2, beliefPenaltyMPC_riub2, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD3_62(params->lb4, beliefPenaltyMPC_z3, beliefPenaltyMPC_lbIdx3, beliefPenaltyMPC_llb3, beliefPenaltyMPC_slb3, beliefPenaltyMPC_rilb3, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD2_22(beliefPenaltyMPC_z3, beliefPenaltyMPC_ubIdx3, params->ub4, beliefPenaltyMPC_lub3, beliefPenaltyMPC_sub3, beliefPenaltyMPC_riub3, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD3_62(params->lb5, beliefPenaltyMPC_z4, beliefPenaltyMPC_lbIdx4, beliefPenaltyMPC_llb4, beliefPenaltyMPC_slb4, beliefPenaltyMPC_rilb4, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD2_22(beliefPenaltyMPC_z4, beliefPenaltyMPC_ubIdx4, params->ub5, beliefPenaltyMPC_lub4, beliefPenaltyMPC_sub4, beliefPenaltyMPC_riub4, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD3_62(params->lb6, beliefPenaltyMPC_z5, beliefPenaltyMPC_lbIdx5, beliefPenaltyMPC_llb5, beliefPenaltyMPC_slb5, beliefPenaltyMPC_rilb5, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD2_22(beliefPenaltyMPC_z5, beliefPenaltyMPC_ubIdx5, params->ub6, beliefPenaltyMPC_lub5, beliefPenaltyMPC_sub5, beliefPenaltyMPC_riub5, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD3_62(params->lb7, beliefPenaltyMPC_z6, beliefPenaltyMPC_lbIdx6, beliefPenaltyMPC_llb6, beliefPenaltyMPC_slb6, beliefPenaltyMPC_rilb6, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD2_22(beliefPenaltyMPC_z6, beliefPenaltyMPC_ubIdx6, params->ub7, beliefPenaltyMPC_lub6, beliefPenaltyMPC_sub6, beliefPenaltyMPC_riub6, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD3_62(params->lb8, beliefPenaltyMPC_z7, beliefPenaltyMPC_lbIdx7, beliefPenaltyMPC_llb7, beliefPenaltyMPC_slb7, beliefPenaltyMPC_rilb7, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD2_22(beliefPenaltyMPC_z7, beliefPenaltyMPC_ubIdx7, params->ub8, beliefPenaltyMPC_lub7, beliefPenaltyMPC_sub7, beliefPenaltyMPC_riub7, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD3_62(params->lb9, beliefPenaltyMPC_z8, beliefPenaltyMPC_lbIdx8, beliefPenaltyMPC_llb8, beliefPenaltyMPC_slb8, beliefPenaltyMPC_rilb8, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD2_22(beliefPenaltyMPC_z8, beliefPenaltyMPC_ubIdx8, params->ub9, beliefPenaltyMPC_lub8, beliefPenaltyMPC_sub8, beliefPenaltyMPC_riub8, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD3_20(params->lb10, beliefPenaltyMPC_z9, beliefPenaltyMPC_lbIdx9, beliefPenaltyMPC_llb9, beliefPenaltyMPC_slb9, beliefPenaltyMPC_rilb9, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD2_20(beliefPenaltyMPC_z9, beliefPenaltyMPC_ubIdx9, params->ub10, beliefPenaltyMPC_lub9, beliefPenaltyMPC_sub9, beliefPenaltyMPC_riub9, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_INEQ_B_GRAD_62_62_22(beliefPenaltyMPC_lub0, beliefPenaltyMPC_sub0, beliefPenaltyMPC_riub0, beliefPenaltyMPC_llb0, beliefPenaltyMPC_slb0, beliefPenaltyMPC_rilb0, beliefPenaltyMPC_lbIdx0, beliefPenaltyMPC_ubIdx0, beliefPenaltyMPC_grad_ineq0, beliefPenaltyMPC_lubbysub0, beliefPenaltyMPC_llbbyslb0);
beliefPenaltyMPC_LA_INEQ_B_GRAD_62_62_22(beliefPenaltyMPC_lub1, beliefPenaltyMPC_sub1, beliefPenaltyMPC_riub1, beliefPenaltyMPC_llb1, beliefPenaltyMPC_slb1, beliefPenaltyMPC_rilb1, beliefPenaltyMPC_lbIdx1, beliefPenaltyMPC_ubIdx1, beliefPenaltyMPC_grad_ineq1, beliefPenaltyMPC_lubbysub1, beliefPenaltyMPC_llbbyslb1);
beliefPenaltyMPC_LA_INEQ_B_GRAD_62_62_22(beliefPenaltyMPC_lub2, beliefPenaltyMPC_sub2, beliefPenaltyMPC_riub2, beliefPenaltyMPC_llb2, beliefPenaltyMPC_slb2, beliefPenaltyMPC_rilb2, beliefPenaltyMPC_lbIdx2, beliefPenaltyMPC_ubIdx2, beliefPenaltyMPC_grad_ineq2, beliefPenaltyMPC_lubbysub2, beliefPenaltyMPC_llbbyslb2);
beliefPenaltyMPC_LA_INEQ_B_GRAD_62_62_22(beliefPenaltyMPC_lub3, beliefPenaltyMPC_sub3, beliefPenaltyMPC_riub3, beliefPenaltyMPC_llb3, beliefPenaltyMPC_slb3, beliefPenaltyMPC_rilb3, beliefPenaltyMPC_lbIdx3, beliefPenaltyMPC_ubIdx3, beliefPenaltyMPC_grad_ineq3, beliefPenaltyMPC_lubbysub3, beliefPenaltyMPC_llbbyslb3);
beliefPenaltyMPC_LA_INEQ_B_GRAD_62_62_22(beliefPenaltyMPC_lub4, beliefPenaltyMPC_sub4, beliefPenaltyMPC_riub4, beliefPenaltyMPC_llb4, beliefPenaltyMPC_slb4, beliefPenaltyMPC_rilb4, beliefPenaltyMPC_lbIdx4, beliefPenaltyMPC_ubIdx4, beliefPenaltyMPC_grad_ineq4, beliefPenaltyMPC_lubbysub4, beliefPenaltyMPC_llbbyslb4);
beliefPenaltyMPC_LA_INEQ_B_GRAD_62_62_22(beliefPenaltyMPC_lub5, beliefPenaltyMPC_sub5, beliefPenaltyMPC_riub5, beliefPenaltyMPC_llb5, beliefPenaltyMPC_slb5, beliefPenaltyMPC_rilb5, beliefPenaltyMPC_lbIdx5, beliefPenaltyMPC_ubIdx5, beliefPenaltyMPC_grad_ineq5, beliefPenaltyMPC_lubbysub5, beliefPenaltyMPC_llbbyslb5);
beliefPenaltyMPC_LA_INEQ_B_GRAD_62_62_22(beliefPenaltyMPC_lub6, beliefPenaltyMPC_sub6, beliefPenaltyMPC_riub6, beliefPenaltyMPC_llb6, beliefPenaltyMPC_slb6, beliefPenaltyMPC_rilb6, beliefPenaltyMPC_lbIdx6, beliefPenaltyMPC_ubIdx6, beliefPenaltyMPC_grad_ineq6, beliefPenaltyMPC_lubbysub6, beliefPenaltyMPC_llbbyslb6);
beliefPenaltyMPC_LA_INEQ_B_GRAD_62_62_22(beliefPenaltyMPC_lub7, beliefPenaltyMPC_sub7, beliefPenaltyMPC_riub7, beliefPenaltyMPC_llb7, beliefPenaltyMPC_slb7, beliefPenaltyMPC_rilb7, beliefPenaltyMPC_lbIdx7, beliefPenaltyMPC_ubIdx7, beliefPenaltyMPC_grad_ineq7, beliefPenaltyMPC_lubbysub7, beliefPenaltyMPC_llbbyslb7);
beliefPenaltyMPC_LA_INEQ_B_GRAD_62_62_22(beliefPenaltyMPC_lub8, beliefPenaltyMPC_sub8, beliefPenaltyMPC_riub8, beliefPenaltyMPC_llb8, beliefPenaltyMPC_slb8, beliefPenaltyMPC_rilb8, beliefPenaltyMPC_lbIdx8, beliefPenaltyMPC_ubIdx8, beliefPenaltyMPC_grad_ineq8, beliefPenaltyMPC_lubbysub8, beliefPenaltyMPC_llbbyslb8);
beliefPenaltyMPC_LA_INEQ_B_GRAD_20_20_20(beliefPenaltyMPC_lub9, beliefPenaltyMPC_sub9, beliefPenaltyMPC_riub9, beliefPenaltyMPC_llb9, beliefPenaltyMPC_slb9, beliefPenaltyMPC_rilb9, beliefPenaltyMPC_lbIdx9, beliefPenaltyMPC_ubIdx9, beliefPenaltyMPC_grad_ineq9, beliefPenaltyMPC_lubbysub9, beliefPenaltyMPC_llbbyslb9);
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
beliefPenaltyMPC_LA_VVADD3_578(beliefPenaltyMPC_grad_cost, beliefPenaltyMPC_grad_eq, beliefPenaltyMPC_grad_ineq, beliefPenaltyMPC_rd);
beliefPenaltyMPC_LA_DIAG_CHOL_LBUB_62_62_22(beliefPenaltyMPC_H0, beliefPenaltyMPC_llbbyslb0, beliefPenaltyMPC_lbIdx0, beliefPenaltyMPC_lubbysub0, beliefPenaltyMPC_ubIdx0, beliefPenaltyMPC_Phi0);
beliefPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_20_62(beliefPenaltyMPC_Phi0, params->C1, beliefPenaltyMPC_V0);
beliefPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_20_62(beliefPenaltyMPC_Phi0, beliefPenaltyMPC_D0, beliefPenaltyMPC_W0);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMTM_20_62_20(beliefPenaltyMPC_W0, beliefPenaltyMPC_V0, beliefPenaltyMPC_Ysd1);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_62(beliefPenaltyMPC_Phi0, beliefPenaltyMPC_rd0, beliefPenaltyMPC_Lbyrd0);
beliefPenaltyMPC_LA_DIAG_CHOL_LBUB_62_62_22(beliefPenaltyMPC_H0, beliefPenaltyMPC_llbbyslb1, beliefPenaltyMPC_lbIdx1, beliefPenaltyMPC_lubbysub1, beliefPenaltyMPC_ubIdx1, beliefPenaltyMPC_Phi1);
beliefPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_20_62(beliefPenaltyMPC_Phi1, params->C2, beliefPenaltyMPC_V1);
beliefPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_20_62(beliefPenaltyMPC_Phi1, beliefPenaltyMPC_D1, beliefPenaltyMPC_W1);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMTM_20_62_20(beliefPenaltyMPC_W1, beliefPenaltyMPC_V1, beliefPenaltyMPC_Ysd2);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_62(beliefPenaltyMPC_Phi1, beliefPenaltyMPC_rd1, beliefPenaltyMPC_Lbyrd1);
beliefPenaltyMPC_LA_DIAG_CHOL_LBUB_62_62_22(beliefPenaltyMPC_H0, beliefPenaltyMPC_llbbyslb2, beliefPenaltyMPC_lbIdx2, beliefPenaltyMPC_lubbysub2, beliefPenaltyMPC_ubIdx2, beliefPenaltyMPC_Phi2);
beliefPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_20_62(beliefPenaltyMPC_Phi2, params->C3, beliefPenaltyMPC_V2);
beliefPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_20_62(beliefPenaltyMPC_Phi2, beliefPenaltyMPC_D1, beliefPenaltyMPC_W2);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMTM_20_62_20(beliefPenaltyMPC_W2, beliefPenaltyMPC_V2, beliefPenaltyMPC_Ysd3);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_62(beliefPenaltyMPC_Phi2, beliefPenaltyMPC_rd2, beliefPenaltyMPC_Lbyrd2);
beliefPenaltyMPC_LA_DIAG_CHOL_LBUB_62_62_22(beliefPenaltyMPC_H0, beliefPenaltyMPC_llbbyslb3, beliefPenaltyMPC_lbIdx3, beliefPenaltyMPC_lubbysub3, beliefPenaltyMPC_ubIdx3, beliefPenaltyMPC_Phi3);
beliefPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_20_62(beliefPenaltyMPC_Phi3, params->C4, beliefPenaltyMPC_V3);
beliefPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_20_62(beliefPenaltyMPC_Phi3, beliefPenaltyMPC_D1, beliefPenaltyMPC_W3);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMTM_20_62_20(beliefPenaltyMPC_W3, beliefPenaltyMPC_V3, beliefPenaltyMPC_Ysd4);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_62(beliefPenaltyMPC_Phi3, beliefPenaltyMPC_rd3, beliefPenaltyMPC_Lbyrd3);
beliefPenaltyMPC_LA_DIAG_CHOL_LBUB_62_62_22(beliefPenaltyMPC_H0, beliefPenaltyMPC_llbbyslb4, beliefPenaltyMPC_lbIdx4, beliefPenaltyMPC_lubbysub4, beliefPenaltyMPC_ubIdx4, beliefPenaltyMPC_Phi4);
beliefPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_20_62(beliefPenaltyMPC_Phi4, params->C5, beliefPenaltyMPC_V4);
beliefPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_20_62(beliefPenaltyMPC_Phi4, beliefPenaltyMPC_D1, beliefPenaltyMPC_W4);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMTM_20_62_20(beliefPenaltyMPC_W4, beliefPenaltyMPC_V4, beliefPenaltyMPC_Ysd5);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_62(beliefPenaltyMPC_Phi4, beliefPenaltyMPC_rd4, beliefPenaltyMPC_Lbyrd4);
beliefPenaltyMPC_LA_DIAG_CHOL_LBUB_62_62_22(beliefPenaltyMPC_H0, beliefPenaltyMPC_llbbyslb5, beliefPenaltyMPC_lbIdx5, beliefPenaltyMPC_lubbysub5, beliefPenaltyMPC_ubIdx5, beliefPenaltyMPC_Phi5);
beliefPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_20_62(beliefPenaltyMPC_Phi5, params->C6, beliefPenaltyMPC_V5);
beliefPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_20_62(beliefPenaltyMPC_Phi5, beliefPenaltyMPC_D1, beliefPenaltyMPC_W5);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMTM_20_62_20(beliefPenaltyMPC_W5, beliefPenaltyMPC_V5, beliefPenaltyMPC_Ysd6);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_62(beliefPenaltyMPC_Phi5, beliefPenaltyMPC_rd5, beliefPenaltyMPC_Lbyrd5);
beliefPenaltyMPC_LA_DIAG_CHOL_LBUB_62_62_22(beliefPenaltyMPC_H0, beliefPenaltyMPC_llbbyslb6, beliefPenaltyMPC_lbIdx6, beliefPenaltyMPC_lubbysub6, beliefPenaltyMPC_ubIdx6, beliefPenaltyMPC_Phi6);
beliefPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_20_62(beliefPenaltyMPC_Phi6, params->C7, beliefPenaltyMPC_V6);
beliefPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_20_62(beliefPenaltyMPC_Phi6, beliefPenaltyMPC_D1, beliefPenaltyMPC_W6);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMTM_20_62_20(beliefPenaltyMPC_W6, beliefPenaltyMPC_V6, beliefPenaltyMPC_Ysd7);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_62(beliefPenaltyMPC_Phi6, beliefPenaltyMPC_rd6, beliefPenaltyMPC_Lbyrd6);
beliefPenaltyMPC_LA_DIAG_CHOL_LBUB_62_62_22(beliefPenaltyMPC_H0, beliefPenaltyMPC_llbbyslb7, beliefPenaltyMPC_lbIdx7, beliefPenaltyMPC_lubbysub7, beliefPenaltyMPC_ubIdx7, beliefPenaltyMPC_Phi7);
beliefPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_20_62(beliefPenaltyMPC_Phi7, params->C8, beliefPenaltyMPC_V7);
beliefPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_20_62(beliefPenaltyMPC_Phi7, beliefPenaltyMPC_D1, beliefPenaltyMPC_W7);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMTM_20_62_20(beliefPenaltyMPC_W7, beliefPenaltyMPC_V7, beliefPenaltyMPC_Ysd8);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_62(beliefPenaltyMPC_Phi7, beliefPenaltyMPC_rd7, beliefPenaltyMPC_Lbyrd7);
beliefPenaltyMPC_LA_DIAG_CHOL_LBUB_62_62_22(beliefPenaltyMPC_H0, beliefPenaltyMPC_llbbyslb8, beliefPenaltyMPC_lbIdx8, beliefPenaltyMPC_lubbysub8, beliefPenaltyMPC_ubIdx8, beliefPenaltyMPC_Phi8);
beliefPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_20_62(beliefPenaltyMPC_Phi8, params->C9, beliefPenaltyMPC_V8);
beliefPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_20_62(beliefPenaltyMPC_Phi8, beliefPenaltyMPC_D1, beliefPenaltyMPC_W8);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMTM_20_62_20(beliefPenaltyMPC_W8, beliefPenaltyMPC_V8, beliefPenaltyMPC_Ysd9);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_62(beliefPenaltyMPC_Phi8, beliefPenaltyMPC_rd8, beliefPenaltyMPC_Lbyrd8);
beliefPenaltyMPC_LA_DIAG_CHOL_ONELOOP_LBUB_20_20_20(beliefPenaltyMPC_H9, beliefPenaltyMPC_llbbyslb9, beliefPenaltyMPC_lbIdx9, beliefPenaltyMPC_lubbysub9, beliefPenaltyMPC_ubIdx9, beliefPenaltyMPC_Phi9);
beliefPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_20_20(beliefPenaltyMPC_Phi9, beliefPenaltyMPC_D9, beliefPenaltyMPC_W9);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_20(beliefPenaltyMPC_Phi9, beliefPenaltyMPC_rd9, beliefPenaltyMPC_Lbyrd9);
beliefPenaltyMPC_LA_DIAGZERO_MMT_20(beliefPenaltyMPC_W0, beliefPenaltyMPC_Yd0);
beliefPenaltyMPC_LA_DIAGZERO_MVMSUB7_20(beliefPenaltyMPC_W0, beliefPenaltyMPC_Lbyrd0, beliefPenaltyMPC_re0, beliefPenaltyMPC_beta0);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMT2_20_62_62(beliefPenaltyMPC_V0, beliefPenaltyMPC_W1, beliefPenaltyMPC_Yd1);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_20_62_62(beliefPenaltyMPC_V0, beliefPenaltyMPC_Lbyrd0, beliefPenaltyMPC_W1, beliefPenaltyMPC_Lbyrd1, beliefPenaltyMPC_re1, beliefPenaltyMPC_beta1);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMT2_20_62_62(beliefPenaltyMPC_V1, beliefPenaltyMPC_W2, beliefPenaltyMPC_Yd2);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_20_62_62(beliefPenaltyMPC_V1, beliefPenaltyMPC_Lbyrd1, beliefPenaltyMPC_W2, beliefPenaltyMPC_Lbyrd2, beliefPenaltyMPC_re2, beliefPenaltyMPC_beta2);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMT2_20_62_62(beliefPenaltyMPC_V2, beliefPenaltyMPC_W3, beliefPenaltyMPC_Yd3);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_20_62_62(beliefPenaltyMPC_V2, beliefPenaltyMPC_Lbyrd2, beliefPenaltyMPC_W3, beliefPenaltyMPC_Lbyrd3, beliefPenaltyMPC_re3, beliefPenaltyMPC_beta3);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMT2_20_62_62(beliefPenaltyMPC_V3, beliefPenaltyMPC_W4, beliefPenaltyMPC_Yd4);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_20_62_62(beliefPenaltyMPC_V3, beliefPenaltyMPC_Lbyrd3, beliefPenaltyMPC_W4, beliefPenaltyMPC_Lbyrd4, beliefPenaltyMPC_re4, beliefPenaltyMPC_beta4);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMT2_20_62_62(beliefPenaltyMPC_V4, beliefPenaltyMPC_W5, beliefPenaltyMPC_Yd5);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_20_62_62(beliefPenaltyMPC_V4, beliefPenaltyMPC_Lbyrd4, beliefPenaltyMPC_W5, beliefPenaltyMPC_Lbyrd5, beliefPenaltyMPC_re5, beliefPenaltyMPC_beta5);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMT2_20_62_62(beliefPenaltyMPC_V5, beliefPenaltyMPC_W6, beliefPenaltyMPC_Yd6);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_20_62_62(beliefPenaltyMPC_V5, beliefPenaltyMPC_Lbyrd5, beliefPenaltyMPC_W6, beliefPenaltyMPC_Lbyrd6, beliefPenaltyMPC_re6, beliefPenaltyMPC_beta6);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMT2_20_62_62(beliefPenaltyMPC_V6, beliefPenaltyMPC_W7, beliefPenaltyMPC_Yd7);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_20_62_62(beliefPenaltyMPC_V6, beliefPenaltyMPC_Lbyrd6, beliefPenaltyMPC_W7, beliefPenaltyMPC_Lbyrd7, beliefPenaltyMPC_re7, beliefPenaltyMPC_beta7);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMT2_20_62_62(beliefPenaltyMPC_V7, beliefPenaltyMPC_W8, beliefPenaltyMPC_Yd8);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_20_62_62(beliefPenaltyMPC_V7, beliefPenaltyMPC_Lbyrd7, beliefPenaltyMPC_W8, beliefPenaltyMPC_Lbyrd8, beliefPenaltyMPC_re8, beliefPenaltyMPC_beta8);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMT2_20_62_20(beliefPenaltyMPC_V8, beliefPenaltyMPC_W9, beliefPenaltyMPC_Yd9);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_20_62_20(beliefPenaltyMPC_V8, beliefPenaltyMPC_Lbyrd8, beliefPenaltyMPC_W9, beliefPenaltyMPC_Lbyrd9, beliefPenaltyMPC_re9, beliefPenaltyMPC_beta9);
beliefPenaltyMPC_LA_DENSE_CHOL_20(beliefPenaltyMPC_Yd0, beliefPenaltyMPC_Ld0);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_20(beliefPenaltyMPC_Ld0, beliefPenaltyMPC_beta0, beliefPenaltyMPC_yy0);
beliefPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_20_20(beliefPenaltyMPC_Ld0, beliefPenaltyMPC_Ysd1, beliefPenaltyMPC_Lsd1);
beliefPenaltyMPC_LA_DENSE_MMTSUB_20_20(beliefPenaltyMPC_Lsd1, beliefPenaltyMPC_Yd1);
beliefPenaltyMPC_LA_DENSE_CHOL_20(beliefPenaltyMPC_Yd1, beliefPenaltyMPC_Ld1);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_20_20(beliefPenaltyMPC_Lsd1, beliefPenaltyMPC_yy0, beliefPenaltyMPC_beta1, beliefPenaltyMPC_bmy1);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_20(beliefPenaltyMPC_Ld1, beliefPenaltyMPC_bmy1, beliefPenaltyMPC_yy1);
beliefPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_20_20(beliefPenaltyMPC_Ld1, beliefPenaltyMPC_Ysd2, beliefPenaltyMPC_Lsd2);
beliefPenaltyMPC_LA_DENSE_MMTSUB_20_20(beliefPenaltyMPC_Lsd2, beliefPenaltyMPC_Yd2);
beliefPenaltyMPC_LA_DENSE_CHOL_20(beliefPenaltyMPC_Yd2, beliefPenaltyMPC_Ld2);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_20_20(beliefPenaltyMPC_Lsd2, beliefPenaltyMPC_yy1, beliefPenaltyMPC_beta2, beliefPenaltyMPC_bmy2);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_20(beliefPenaltyMPC_Ld2, beliefPenaltyMPC_bmy2, beliefPenaltyMPC_yy2);
beliefPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_20_20(beliefPenaltyMPC_Ld2, beliefPenaltyMPC_Ysd3, beliefPenaltyMPC_Lsd3);
beliefPenaltyMPC_LA_DENSE_MMTSUB_20_20(beliefPenaltyMPC_Lsd3, beliefPenaltyMPC_Yd3);
beliefPenaltyMPC_LA_DENSE_CHOL_20(beliefPenaltyMPC_Yd3, beliefPenaltyMPC_Ld3);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_20_20(beliefPenaltyMPC_Lsd3, beliefPenaltyMPC_yy2, beliefPenaltyMPC_beta3, beliefPenaltyMPC_bmy3);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_20(beliefPenaltyMPC_Ld3, beliefPenaltyMPC_bmy3, beliefPenaltyMPC_yy3);
beliefPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_20_20(beliefPenaltyMPC_Ld3, beliefPenaltyMPC_Ysd4, beliefPenaltyMPC_Lsd4);
beliefPenaltyMPC_LA_DENSE_MMTSUB_20_20(beliefPenaltyMPC_Lsd4, beliefPenaltyMPC_Yd4);
beliefPenaltyMPC_LA_DENSE_CHOL_20(beliefPenaltyMPC_Yd4, beliefPenaltyMPC_Ld4);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_20_20(beliefPenaltyMPC_Lsd4, beliefPenaltyMPC_yy3, beliefPenaltyMPC_beta4, beliefPenaltyMPC_bmy4);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_20(beliefPenaltyMPC_Ld4, beliefPenaltyMPC_bmy4, beliefPenaltyMPC_yy4);
beliefPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_20_20(beliefPenaltyMPC_Ld4, beliefPenaltyMPC_Ysd5, beliefPenaltyMPC_Lsd5);
beliefPenaltyMPC_LA_DENSE_MMTSUB_20_20(beliefPenaltyMPC_Lsd5, beliefPenaltyMPC_Yd5);
beliefPenaltyMPC_LA_DENSE_CHOL_20(beliefPenaltyMPC_Yd5, beliefPenaltyMPC_Ld5);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_20_20(beliefPenaltyMPC_Lsd5, beliefPenaltyMPC_yy4, beliefPenaltyMPC_beta5, beliefPenaltyMPC_bmy5);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_20(beliefPenaltyMPC_Ld5, beliefPenaltyMPC_bmy5, beliefPenaltyMPC_yy5);
beliefPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_20_20(beliefPenaltyMPC_Ld5, beliefPenaltyMPC_Ysd6, beliefPenaltyMPC_Lsd6);
beliefPenaltyMPC_LA_DENSE_MMTSUB_20_20(beliefPenaltyMPC_Lsd6, beliefPenaltyMPC_Yd6);
beliefPenaltyMPC_LA_DENSE_CHOL_20(beliefPenaltyMPC_Yd6, beliefPenaltyMPC_Ld6);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_20_20(beliefPenaltyMPC_Lsd6, beliefPenaltyMPC_yy5, beliefPenaltyMPC_beta6, beliefPenaltyMPC_bmy6);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_20(beliefPenaltyMPC_Ld6, beliefPenaltyMPC_bmy6, beliefPenaltyMPC_yy6);
beliefPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_20_20(beliefPenaltyMPC_Ld6, beliefPenaltyMPC_Ysd7, beliefPenaltyMPC_Lsd7);
beliefPenaltyMPC_LA_DENSE_MMTSUB_20_20(beliefPenaltyMPC_Lsd7, beliefPenaltyMPC_Yd7);
beliefPenaltyMPC_LA_DENSE_CHOL_20(beliefPenaltyMPC_Yd7, beliefPenaltyMPC_Ld7);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_20_20(beliefPenaltyMPC_Lsd7, beliefPenaltyMPC_yy6, beliefPenaltyMPC_beta7, beliefPenaltyMPC_bmy7);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_20(beliefPenaltyMPC_Ld7, beliefPenaltyMPC_bmy7, beliefPenaltyMPC_yy7);
beliefPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_20_20(beliefPenaltyMPC_Ld7, beliefPenaltyMPC_Ysd8, beliefPenaltyMPC_Lsd8);
beliefPenaltyMPC_LA_DENSE_MMTSUB_20_20(beliefPenaltyMPC_Lsd8, beliefPenaltyMPC_Yd8);
beliefPenaltyMPC_LA_DENSE_CHOL_20(beliefPenaltyMPC_Yd8, beliefPenaltyMPC_Ld8);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_20_20(beliefPenaltyMPC_Lsd8, beliefPenaltyMPC_yy7, beliefPenaltyMPC_beta8, beliefPenaltyMPC_bmy8);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_20(beliefPenaltyMPC_Ld8, beliefPenaltyMPC_bmy8, beliefPenaltyMPC_yy8);
beliefPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_20_20(beliefPenaltyMPC_Ld8, beliefPenaltyMPC_Ysd9, beliefPenaltyMPC_Lsd9);
beliefPenaltyMPC_LA_DENSE_MMTSUB_20_20(beliefPenaltyMPC_Lsd9, beliefPenaltyMPC_Yd9);
beliefPenaltyMPC_LA_DENSE_CHOL_20(beliefPenaltyMPC_Yd9, beliefPenaltyMPC_Ld9);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_20_20(beliefPenaltyMPC_Lsd9, beliefPenaltyMPC_yy8, beliefPenaltyMPC_beta9, beliefPenaltyMPC_bmy9);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_20(beliefPenaltyMPC_Ld9, beliefPenaltyMPC_bmy9, beliefPenaltyMPC_yy9);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_20(beliefPenaltyMPC_Ld9, beliefPenaltyMPC_yy9, beliefPenaltyMPC_dvaff9);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_20_20(beliefPenaltyMPC_Lsd9, beliefPenaltyMPC_dvaff9, beliefPenaltyMPC_yy8, beliefPenaltyMPC_bmy8);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_20(beliefPenaltyMPC_Ld8, beliefPenaltyMPC_bmy8, beliefPenaltyMPC_dvaff8);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_20_20(beliefPenaltyMPC_Lsd8, beliefPenaltyMPC_dvaff8, beliefPenaltyMPC_yy7, beliefPenaltyMPC_bmy7);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_20(beliefPenaltyMPC_Ld7, beliefPenaltyMPC_bmy7, beliefPenaltyMPC_dvaff7);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_20_20(beliefPenaltyMPC_Lsd7, beliefPenaltyMPC_dvaff7, beliefPenaltyMPC_yy6, beliefPenaltyMPC_bmy6);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_20(beliefPenaltyMPC_Ld6, beliefPenaltyMPC_bmy6, beliefPenaltyMPC_dvaff6);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_20_20(beliefPenaltyMPC_Lsd6, beliefPenaltyMPC_dvaff6, beliefPenaltyMPC_yy5, beliefPenaltyMPC_bmy5);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_20(beliefPenaltyMPC_Ld5, beliefPenaltyMPC_bmy5, beliefPenaltyMPC_dvaff5);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_20_20(beliefPenaltyMPC_Lsd5, beliefPenaltyMPC_dvaff5, beliefPenaltyMPC_yy4, beliefPenaltyMPC_bmy4);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_20(beliefPenaltyMPC_Ld4, beliefPenaltyMPC_bmy4, beliefPenaltyMPC_dvaff4);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_20_20(beliefPenaltyMPC_Lsd4, beliefPenaltyMPC_dvaff4, beliefPenaltyMPC_yy3, beliefPenaltyMPC_bmy3);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_20(beliefPenaltyMPC_Ld3, beliefPenaltyMPC_bmy3, beliefPenaltyMPC_dvaff3);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_20_20(beliefPenaltyMPC_Lsd3, beliefPenaltyMPC_dvaff3, beliefPenaltyMPC_yy2, beliefPenaltyMPC_bmy2);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_20(beliefPenaltyMPC_Ld2, beliefPenaltyMPC_bmy2, beliefPenaltyMPC_dvaff2);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_20_20(beliefPenaltyMPC_Lsd2, beliefPenaltyMPC_dvaff2, beliefPenaltyMPC_yy1, beliefPenaltyMPC_bmy1);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_20(beliefPenaltyMPC_Ld1, beliefPenaltyMPC_bmy1, beliefPenaltyMPC_dvaff1);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_20_20(beliefPenaltyMPC_Lsd1, beliefPenaltyMPC_dvaff1, beliefPenaltyMPC_yy0, beliefPenaltyMPC_bmy0);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_20(beliefPenaltyMPC_Ld0, beliefPenaltyMPC_bmy0, beliefPenaltyMPC_dvaff0);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_20_62_20(params->C1, beliefPenaltyMPC_dvaff1, beliefPenaltyMPC_D0, beliefPenaltyMPC_dvaff0, beliefPenaltyMPC_grad_eq0);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_20_62_20(params->C2, beliefPenaltyMPC_dvaff2, beliefPenaltyMPC_D1, beliefPenaltyMPC_dvaff1, beliefPenaltyMPC_grad_eq1);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_20_62_20(params->C3, beliefPenaltyMPC_dvaff3, beliefPenaltyMPC_D1, beliefPenaltyMPC_dvaff2, beliefPenaltyMPC_grad_eq2);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_20_62_20(params->C4, beliefPenaltyMPC_dvaff4, beliefPenaltyMPC_D1, beliefPenaltyMPC_dvaff3, beliefPenaltyMPC_grad_eq3);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_20_62_20(params->C5, beliefPenaltyMPC_dvaff5, beliefPenaltyMPC_D1, beliefPenaltyMPC_dvaff4, beliefPenaltyMPC_grad_eq4);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_20_62_20(params->C6, beliefPenaltyMPC_dvaff6, beliefPenaltyMPC_D1, beliefPenaltyMPC_dvaff5, beliefPenaltyMPC_grad_eq5);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_20_62_20(params->C7, beliefPenaltyMPC_dvaff7, beliefPenaltyMPC_D1, beliefPenaltyMPC_dvaff6, beliefPenaltyMPC_grad_eq6);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_20_62_20(params->C8, beliefPenaltyMPC_dvaff8, beliefPenaltyMPC_D1, beliefPenaltyMPC_dvaff7, beliefPenaltyMPC_grad_eq7);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_20_62_20(params->C9, beliefPenaltyMPC_dvaff9, beliefPenaltyMPC_D1, beliefPenaltyMPC_dvaff8, beliefPenaltyMPC_grad_eq8);
beliefPenaltyMPC_LA_DIAGZERO_MTVM_20_20(beliefPenaltyMPC_D9, beliefPenaltyMPC_dvaff9, beliefPenaltyMPC_grad_eq9);
beliefPenaltyMPC_LA_VSUB2_578(beliefPenaltyMPC_rd, beliefPenaltyMPC_grad_eq, beliefPenaltyMPC_rd);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_62(beliefPenaltyMPC_Phi0, beliefPenaltyMPC_rd0, beliefPenaltyMPC_dzaff0);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_62(beliefPenaltyMPC_Phi1, beliefPenaltyMPC_rd1, beliefPenaltyMPC_dzaff1);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_62(beliefPenaltyMPC_Phi2, beliefPenaltyMPC_rd2, beliefPenaltyMPC_dzaff2);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_62(beliefPenaltyMPC_Phi3, beliefPenaltyMPC_rd3, beliefPenaltyMPC_dzaff3);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_62(beliefPenaltyMPC_Phi4, beliefPenaltyMPC_rd4, beliefPenaltyMPC_dzaff4);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_62(beliefPenaltyMPC_Phi5, beliefPenaltyMPC_rd5, beliefPenaltyMPC_dzaff5);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_62(beliefPenaltyMPC_Phi6, beliefPenaltyMPC_rd6, beliefPenaltyMPC_dzaff6);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_62(beliefPenaltyMPC_Phi7, beliefPenaltyMPC_rd7, beliefPenaltyMPC_dzaff7);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_62(beliefPenaltyMPC_Phi8, beliefPenaltyMPC_rd8, beliefPenaltyMPC_dzaff8);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_20(beliefPenaltyMPC_Phi9, beliefPenaltyMPC_rd9, beliefPenaltyMPC_dzaff9);
beliefPenaltyMPC_LA_VSUB_INDEXED_62(beliefPenaltyMPC_dzaff0, beliefPenaltyMPC_lbIdx0, beliefPenaltyMPC_rilb0, beliefPenaltyMPC_dslbaff0);
beliefPenaltyMPC_LA_VSUB3_62(beliefPenaltyMPC_llbbyslb0, beliefPenaltyMPC_dslbaff0, beliefPenaltyMPC_llb0, beliefPenaltyMPC_dllbaff0);
beliefPenaltyMPC_LA_VSUB2_INDEXED_22(beliefPenaltyMPC_riub0, beliefPenaltyMPC_dzaff0, beliefPenaltyMPC_ubIdx0, beliefPenaltyMPC_dsubaff0);
beliefPenaltyMPC_LA_VSUB3_22(beliefPenaltyMPC_lubbysub0, beliefPenaltyMPC_dsubaff0, beliefPenaltyMPC_lub0, beliefPenaltyMPC_dlubaff0);
beliefPenaltyMPC_LA_VSUB_INDEXED_62(beliefPenaltyMPC_dzaff1, beliefPenaltyMPC_lbIdx1, beliefPenaltyMPC_rilb1, beliefPenaltyMPC_dslbaff1);
beliefPenaltyMPC_LA_VSUB3_62(beliefPenaltyMPC_llbbyslb1, beliefPenaltyMPC_dslbaff1, beliefPenaltyMPC_llb1, beliefPenaltyMPC_dllbaff1);
beliefPenaltyMPC_LA_VSUB2_INDEXED_22(beliefPenaltyMPC_riub1, beliefPenaltyMPC_dzaff1, beliefPenaltyMPC_ubIdx1, beliefPenaltyMPC_dsubaff1);
beliefPenaltyMPC_LA_VSUB3_22(beliefPenaltyMPC_lubbysub1, beliefPenaltyMPC_dsubaff1, beliefPenaltyMPC_lub1, beliefPenaltyMPC_dlubaff1);
beliefPenaltyMPC_LA_VSUB_INDEXED_62(beliefPenaltyMPC_dzaff2, beliefPenaltyMPC_lbIdx2, beliefPenaltyMPC_rilb2, beliefPenaltyMPC_dslbaff2);
beliefPenaltyMPC_LA_VSUB3_62(beliefPenaltyMPC_llbbyslb2, beliefPenaltyMPC_dslbaff2, beliefPenaltyMPC_llb2, beliefPenaltyMPC_dllbaff2);
beliefPenaltyMPC_LA_VSUB2_INDEXED_22(beliefPenaltyMPC_riub2, beliefPenaltyMPC_dzaff2, beliefPenaltyMPC_ubIdx2, beliefPenaltyMPC_dsubaff2);
beliefPenaltyMPC_LA_VSUB3_22(beliefPenaltyMPC_lubbysub2, beliefPenaltyMPC_dsubaff2, beliefPenaltyMPC_lub2, beliefPenaltyMPC_dlubaff2);
beliefPenaltyMPC_LA_VSUB_INDEXED_62(beliefPenaltyMPC_dzaff3, beliefPenaltyMPC_lbIdx3, beliefPenaltyMPC_rilb3, beliefPenaltyMPC_dslbaff3);
beliefPenaltyMPC_LA_VSUB3_62(beliefPenaltyMPC_llbbyslb3, beliefPenaltyMPC_dslbaff3, beliefPenaltyMPC_llb3, beliefPenaltyMPC_dllbaff3);
beliefPenaltyMPC_LA_VSUB2_INDEXED_22(beliefPenaltyMPC_riub3, beliefPenaltyMPC_dzaff3, beliefPenaltyMPC_ubIdx3, beliefPenaltyMPC_dsubaff3);
beliefPenaltyMPC_LA_VSUB3_22(beliefPenaltyMPC_lubbysub3, beliefPenaltyMPC_dsubaff3, beliefPenaltyMPC_lub3, beliefPenaltyMPC_dlubaff3);
beliefPenaltyMPC_LA_VSUB_INDEXED_62(beliefPenaltyMPC_dzaff4, beliefPenaltyMPC_lbIdx4, beliefPenaltyMPC_rilb4, beliefPenaltyMPC_dslbaff4);
beliefPenaltyMPC_LA_VSUB3_62(beliefPenaltyMPC_llbbyslb4, beliefPenaltyMPC_dslbaff4, beliefPenaltyMPC_llb4, beliefPenaltyMPC_dllbaff4);
beliefPenaltyMPC_LA_VSUB2_INDEXED_22(beliefPenaltyMPC_riub4, beliefPenaltyMPC_dzaff4, beliefPenaltyMPC_ubIdx4, beliefPenaltyMPC_dsubaff4);
beliefPenaltyMPC_LA_VSUB3_22(beliefPenaltyMPC_lubbysub4, beliefPenaltyMPC_dsubaff4, beliefPenaltyMPC_lub4, beliefPenaltyMPC_dlubaff4);
beliefPenaltyMPC_LA_VSUB_INDEXED_62(beliefPenaltyMPC_dzaff5, beliefPenaltyMPC_lbIdx5, beliefPenaltyMPC_rilb5, beliefPenaltyMPC_dslbaff5);
beliefPenaltyMPC_LA_VSUB3_62(beliefPenaltyMPC_llbbyslb5, beliefPenaltyMPC_dslbaff5, beliefPenaltyMPC_llb5, beliefPenaltyMPC_dllbaff5);
beliefPenaltyMPC_LA_VSUB2_INDEXED_22(beliefPenaltyMPC_riub5, beliefPenaltyMPC_dzaff5, beliefPenaltyMPC_ubIdx5, beliefPenaltyMPC_dsubaff5);
beliefPenaltyMPC_LA_VSUB3_22(beliefPenaltyMPC_lubbysub5, beliefPenaltyMPC_dsubaff5, beliefPenaltyMPC_lub5, beliefPenaltyMPC_dlubaff5);
beliefPenaltyMPC_LA_VSUB_INDEXED_62(beliefPenaltyMPC_dzaff6, beliefPenaltyMPC_lbIdx6, beliefPenaltyMPC_rilb6, beliefPenaltyMPC_dslbaff6);
beliefPenaltyMPC_LA_VSUB3_62(beliefPenaltyMPC_llbbyslb6, beliefPenaltyMPC_dslbaff6, beliefPenaltyMPC_llb6, beliefPenaltyMPC_dllbaff6);
beliefPenaltyMPC_LA_VSUB2_INDEXED_22(beliefPenaltyMPC_riub6, beliefPenaltyMPC_dzaff6, beliefPenaltyMPC_ubIdx6, beliefPenaltyMPC_dsubaff6);
beliefPenaltyMPC_LA_VSUB3_22(beliefPenaltyMPC_lubbysub6, beliefPenaltyMPC_dsubaff6, beliefPenaltyMPC_lub6, beliefPenaltyMPC_dlubaff6);
beliefPenaltyMPC_LA_VSUB_INDEXED_62(beliefPenaltyMPC_dzaff7, beliefPenaltyMPC_lbIdx7, beliefPenaltyMPC_rilb7, beliefPenaltyMPC_dslbaff7);
beliefPenaltyMPC_LA_VSUB3_62(beliefPenaltyMPC_llbbyslb7, beliefPenaltyMPC_dslbaff7, beliefPenaltyMPC_llb7, beliefPenaltyMPC_dllbaff7);
beliefPenaltyMPC_LA_VSUB2_INDEXED_22(beliefPenaltyMPC_riub7, beliefPenaltyMPC_dzaff7, beliefPenaltyMPC_ubIdx7, beliefPenaltyMPC_dsubaff7);
beliefPenaltyMPC_LA_VSUB3_22(beliefPenaltyMPC_lubbysub7, beliefPenaltyMPC_dsubaff7, beliefPenaltyMPC_lub7, beliefPenaltyMPC_dlubaff7);
beliefPenaltyMPC_LA_VSUB_INDEXED_62(beliefPenaltyMPC_dzaff8, beliefPenaltyMPC_lbIdx8, beliefPenaltyMPC_rilb8, beliefPenaltyMPC_dslbaff8);
beliefPenaltyMPC_LA_VSUB3_62(beliefPenaltyMPC_llbbyslb8, beliefPenaltyMPC_dslbaff8, beliefPenaltyMPC_llb8, beliefPenaltyMPC_dllbaff8);
beliefPenaltyMPC_LA_VSUB2_INDEXED_22(beliefPenaltyMPC_riub8, beliefPenaltyMPC_dzaff8, beliefPenaltyMPC_ubIdx8, beliefPenaltyMPC_dsubaff8);
beliefPenaltyMPC_LA_VSUB3_22(beliefPenaltyMPC_lubbysub8, beliefPenaltyMPC_dsubaff8, beliefPenaltyMPC_lub8, beliefPenaltyMPC_dlubaff8);
beliefPenaltyMPC_LA_VSUB_INDEXED_20(beliefPenaltyMPC_dzaff9, beliefPenaltyMPC_lbIdx9, beliefPenaltyMPC_rilb9, beliefPenaltyMPC_dslbaff9);
beliefPenaltyMPC_LA_VSUB3_20(beliefPenaltyMPC_llbbyslb9, beliefPenaltyMPC_dslbaff9, beliefPenaltyMPC_llb9, beliefPenaltyMPC_dllbaff9);
beliefPenaltyMPC_LA_VSUB2_INDEXED_20(beliefPenaltyMPC_riub9, beliefPenaltyMPC_dzaff9, beliefPenaltyMPC_ubIdx9, beliefPenaltyMPC_dsubaff9);
beliefPenaltyMPC_LA_VSUB3_20(beliefPenaltyMPC_lubbysub9, beliefPenaltyMPC_dsubaff9, beliefPenaltyMPC_lub9, beliefPenaltyMPC_dlubaff9);
info->lsit_aff = beliefPenaltyMPC_LINESEARCH_BACKTRACKING_AFFINE(beliefPenaltyMPC_l, beliefPenaltyMPC_s, beliefPenaltyMPC_dl_aff, beliefPenaltyMPC_ds_aff, &info->step_aff, &info->mu_aff);
if( info->lsit_aff == beliefPenaltyMPC_NOPROGRESS ){
exitcode = beliefPenaltyMPC_NOPROGRESS; break;
}
sigma_3rdroot = info->mu_aff / info->mu;
info->sigma = sigma_3rdroot*sigma_3rdroot*sigma_3rdroot;
musigma = info->mu * info->sigma;
beliefPenaltyMPC_LA_VSUB5_796(beliefPenaltyMPC_ds_aff, beliefPenaltyMPC_dl_aff, musigma, beliefPenaltyMPC_ccrhs);
beliefPenaltyMPC_LA_VSUB6_INDEXED_62_22_62(beliefPenaltyMPC_ccrhsub0, beliefPenaltyMPC_sub0, beliefPenaltyMPC_ubIdx0, beliefPenaltyMPC_ccrhsl0, beliefPenaltyMPC_slb0, beliefPenaltyMPC_lbIdx0, beliefPenaltyMPC_rd0);
beliefPenaltyMPC_LA_VSUB6_INDEXED_62_22_62(beliefPenaltyMPC_ccrhsub1, beliefPenaltyMPC_sub1, beliefPenaltyMPC_ubIdx1, beliefPenaltyMPC_ccrhsl1, beliefPenaltyMPC_slb1, beliefPenaltyMPC_lbIdx1, beliefPenaltyMPC_rd1);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_62(beliefPenaltyMPC_Phi0, beliefPenaltyMPC_rd0, beliefPenaltyMPC_Lbyrd0);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_62(beliefPenaltyMPC_Phi1, beliefPenaltyMPC_rd1, beliefPenaltyMPC_Lbyrd1);
beliefPenaltyMPC_LA_DIAGZERO_MVM_20(beliefPenaltyMPC_W0, beliefPenaltyMPC_Lbyrd0, beliefPenaltyMPC_beta0);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_20(beliefPenaltyMPC_Ld0, beliefPenaltyMPC_beta0, beliefPenaltyMPC_yy0);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_20_62_62(beliefPenaltyMPC_V0, beliefPenaltyMPC_Lbyrd0, beliefPenaltyMPC_W1, beliefPenaltyMPC_Lbyrd1, beliefPenaltyMPC_beta1);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_20_20(beliefPenaltyMPC_Lsd1, beliefPenaltyMPC_yy0, beliefPenaltyMPC_beta1, beliefPenaltyMPC_bmy1);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_20(beliefPenaltyMPC_Ld1, beliefPenaltyMPC_bmy1, beliefPenaltyMPC_yy1);
beliefPenaltyMPC_LA_VSUB6_INDEXED_62_22_62(beliefPenaltyMPC_ccrhsub2, beliefPenaltyMPC_sub2, beliefPenaltyMPC_ubIdx2, beliefPenaltyMPC_ccrhsl2, beliefPenaltyMPC_slb2, beliefPenaltyMPC_lbIdx2, beliefPenaltyMPC_rd2);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_62(beliefPenaltyMPC_Phi2, beliefPenaltyMPC_rd2, beliefPenaltyMPC_Lbyrd2);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_20_62_62(beliefPenaltyMPC_V1, beliefPenaltyMPC_Lbyrd1, beliefPenaltyMPC_W2, beliefPenaltyMPC_Lbyrd2, beliefPenaltyMPC_beta2);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_20_20(beliefPenaltyMPC_Lsd2, beliefPenaltyMPC_yy1, beliefPenaltyMPC_beta2, beliefPenaltyMPC_bmy2);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_20(beliefPenaltyMPC_Ld2, beliefPenaltyMPC_bmy2, beliefPenaltyMPC_yy2);
beliefPenaltyMPC_LA_VSUB6_INDEXED_62_22_62(beliefPenaltyMPC_ccrhsub3, beliefPenaltyMPC_sub3, beliefPenaltyMPC_ubIdx3, beliefPenaltyMPC_ccrhsl3, beliefPenaltyMPC_slb3, beliefPenaltyMPC_lbIdx3, beliefPenaltyMPC_rd3);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_62(beliefPenaltyMPC_Phi3, beliefPenaltyMPC_rd3, beliefPenaltyMPC_Lbyrd3);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_20_62_62(beliefPenaltyMPC_V2, beliefPenaltyMPC_Lbyrd2, beliefPenaltyMPC_W3, beliefPenaltyMPC_Lbyrd3, beliefPenaltyMPC_beta3);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_20_20(beliefPenaltyMPC_Lsd3, beliefPenaltyMPC_yy2, beliefPenaltyMPC_beta3, beliefPenaltyMPC_bmy3);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_20(beliefPenaltyMPC_Ld3, beliefPenaltyMPC_bmy3, beliefPenaltyMPC_yy3);
beliefPenaltyMPC_LA_VSUB6_INDEXED_62_22_62(beliefPenaltyMPC_ccrhsub4, beliefPenaltyMPC_sub4, beliefPenaltyMPC_ubIdx4, beliefPenaltyMPC_ccrhsl4, beliefPenaltyMPC_slb4, beliefPenaltyMPC_lbIdx4, beliefPenaltyMPC_rd4);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_62(beliefPenaltyMPC_Phi4, beliefPenaltyMPC_rd4, beliefPenaltyMPC_Lbyrd4);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_20_62_62(beliefPenaltyMPC_V3, beliefPenaltyMPC_Lbyrd3, beliefPenaltyMPC_W4, beliefPenaltyMPC_Lbyrd4, beliefPenaltyMPC_beta4);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_20_20(beliefPenaltyMPC_Lsd4, beliefPenaltyMPC_yy3, beliefPenaltyMPC_beta4, beliefPenaltyMPC_bmy4);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_20(beliefPenaltyMPC_Ld4, beliefPenaltyMPC_bmy4, beliefPenaltyMPC_yy4);
beliefPenaltyMPC_LA_VSUB6_INDEXED_62_22_62(beliefPenaltyMPC_ccrhsub5, beliefPenaltyMPC_sub5, beliefPenaltyMPC_ubIdx5, beliefPenaltyMPC_ccrhsl5, beliefPenaltyMPC_slb5, beliefPenaltyMPC_lbIdx5, beliefPenaltyMPC_rd5);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_62(beliefPenaltyMPC_Phi5, beliefPenaltyMPC_rd5, beliefPenaltyMPC_Lbyrd5);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_20_62_62(beliefPenaltyMPC_V4, beliefPenaltyMPC_Lbyrd4, beliefPenaltyMPC_W5, beliefPenaltyMPC_Lbyrd5, beliefPenaltyMPC_beta5);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_20_20(beliefPenaltyMPC_Lsd5, beliefPenaltyMPC_yy4, beliefPenaltyMPC_beta5, beliefPenaltyMPC_bmy5);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_20(beliefPenaltyMPC_Ld5, beliefPenaltyMPC_bmy5, beliefPenaltyMPC_yy5);
beliefPenaltyMPC_LA_VSUB6_INDEXED_62_22_62(beliefPenaltyMPC_ccrhsub6, beliefPenaltyMPC_sub6, beliefPenaltyMPC_ubIdx6, beliefPenaltyMPC_ccrhsl6, beliefPenaltyMPC_slb6, beliefPenaltyMPC_lbIdx6, beliefPenaltyMPC_rd6);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_62(beliefPenaltyMPC_Phi6, beliefPenaltyMPC_rd6, beliefPenaltyMPC_Lbyrd6);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_20_62_62(beliefPenaltyMPC_V5, beliefPenaltyMPC_Lbyrd5, beliefPenaltyMPC_W6, beliefPenaltyMPC_Lbyrd6, beliefPenaltyMPC_beta6);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_20_20(beliefPenaltyMPC_Lsd6, beliefPenaltyMPC_yy5, beliefPenaltyMPC_beta6, beliefPenaltyMPC_bmy6);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_20(beliefPenaltyMPC_Ld6, beliefPenaltyMPC_bmy6, beliefPenaltyMPC_yy6);
beliefPenaltyMPC_LA_VSUB6_INDEXED_62_22_62(beliefPenaltyMPC_ccrhsub7, beliefPenaltyMPC_sub7, beliefPenaltyMPC_ubIdx7, beliefPenaltyMPC_ccrhsl7, beliefPenaltyMPC_slb7, beliefPenaltyMPC_lbIdx7, beliefPenaltyMPC_rd7);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_62(beliefPenaltyMPC_Phi7, beliefPenaltyMPC_rd7, beliefPenaltyMPC_Lbyrd7);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_20_62_62(beliefPenaltyMPC_V6, beliefPenaltyMPC_Lbyrd6, beliefPenaltyMPC_W7, beliefPenaltyMPC_Lbyrd7, beliefPenaltyMPC_beta7);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_20_20(beliefPenaltyMPC_Lsd7, beliefPenaltyMPC_yy6, beliefPenaltyMPC_beta7, beliefPenaltyMPC_bmy7);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_20(beliefPenaltyMPC_Ld7, beliefPenaltyMPC_bmy7, beliefPenaltyMPC_yy7);
beliefPenaltyMPC_LA_VSUB6_INDEXED_62_22_62(beliefPenaltyMPC_ccrhsub8, beliefPenaltyMPC_sub8, beliefPenaltyMPC_ubIdx8, beliefPenaltyMPC_ccrhsl8, beliefPenaltyMPC_slb8, beliefPenaltyMPC_lbIdx8, beliefPenaltyMPC_rd8);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_62(beliefPenaltyMPC_Phi8, beliefPenaltyMPC_rd8, beliefPenaltyMPC_Lbyrd8);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_20_62_62(beliefPenaltyMPC_V7, beliefPenaltyMPC_Lbyrd7, beliefPenaltyMPC_W8, beliefPenaltyMPC_Lbyrd8, beliefPenaltyMPC_beta8);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_20_20(beliefPenaltyMPC_Lsd8, beliefPenaltyMPC_yy7, beliefPenaltyMPC_beta8, beliefPenaltyMPC_bmy8);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_20(beliefPenaltyMPC_Ld8, beliefPenaltyMPC_bmy8, beliefPenaltyMPC_yy8);
beliefPenaltyMPC_LA_VSUB6_INDEXED_20_20_20(beliefPenaltyMPC_ccrhsub9, beliefPenaltyMPC_sub9, beliefPenaltyMPC_ubIdx9, beliefPenaltyMPC_ccrhsl9, beliefPenaltyMPC_slb9, beliefPenaltyMPC_lbIdx9, beliefPenaltyMPC_rd9);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_20(beliefPenaltyMPC_Phi9, beliefPenaltyMPC_rd9, beliefPenaltyMPC_Lbyrd9);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_20_62_20(beliefPenaltyMPC_V8, beliefPenaltyMPC_Lbyrd8, beliefPenaltyMPC_W9, beliefPenaltyMPC_Lbyrd9, beliefPenaltyMPC_beta9);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_20_20(beliefPenaltyMPC_Lsd9, beliefPenaltyMPC_yy8, beliefPenaltyMPC_beta9, beliefPenaltyMPC_bmy9);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_20(beliefPenaltyMPC_Ld9, beliefPenaltyMPC_bmy9, beliefPenaltyMPC_yy9);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_20(beliefPenaltyMPC_Ld9, beliefPenaltyMPC_yy9, beliefPenaltyMPC_dvcc9);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_20_20(beliefPenaltyMPC_Lsd9, beliefPenaltyMPC_dvcc9, beliefPenaltyMPC_yy8, beliefPenaltyMPC_bmy8);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_20(beliefPenaltyMPC_Ld8, beliefPenaltyMPC_bmy8, beliefPenaltyMPC_dvcc8);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_20_20(beliefPenaltyMPC_Lsd8, beliefPenaltyMPC_dvcc8, beliefPenaltyMPC_yy7, beliefPenaltyMPC_bmy7);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_20(beliefPenaltyMPC_Ld7, beliefPenaltyMPC_bmy7, beliefPenaltyMPC_dvcc7);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_20_20(beliefPenaltyMPC_Lsd7, beliefPenaltyMPC_dvcc7, beliefPenaltyMPC_yy6, beliefPenaltyMPC_bmy6);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_20(beliefPenaltyMPC_Ld6, beliefPenaltyMPC_bmy6, beliefPenaltyMPC_dvcc6);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_20_20(beliefPenaltyMPC_Lsd6, beliefPenaltyMPC_dvcc6, beliefPenaltyMPC_yy5, beliefPenaltyMPC_bmy5);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_20(beliefPenaltyMPC_Ld5, beliefPenaltyMPC_bmy5, beliefPenaltyMPC_dvcc5);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_20_20(beliefPenaltyMPC_Lsd5, beliefPenaltyMPC_dvcc5, beliefPenaltyMPC_yy4, beliefPenaltyMPC_bmy4);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_20(beliefPenaltyMPC_Ld4, beliefPenaltyMPC_bmy4, beliefPenaltyMPC_dvcc4);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_20_20(beliefPenaltyMPC_Lsd4, beliefPenaltyMPC_dvcc4, beliefPenaltyMPC_yy3, beliefPenaltyMPC_bmy3);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_20(beliefPenaltyMPC_Ld3, beliefPenaltyMPC_bmy3, beliefPenaltyMPC_dvcc3);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_20_20(beliefPenaltyMPC_Lsd3, beliefPenaltyMPC_dvcc3, beliefPenaltyMPC_yy2, beliefPenaltyMPC_bmy2);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_20(beliefPenaltyMPC_Ld2, beliefPenaltyMPC_bmy2, beliefPenaltyMPC_dvcc2);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_20_20(beliefPenaltyMPC_Lsd2, beliefPenaltyMPC_dvcc2, beliefPenaltyMPC_yy1, beliefPenaltyMPC_bmy1);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_20(beliefPenaltyMPC_Ld1, beliefPenaltyMPC_bmy1, beliefPenaltyMPC_dvcc1);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_20_20(beliefPenaltyMPC_Lsd1, beliefPenaltyMPC_dvcc1, beliefPenaltyMPC_yy0, beliefPenaltyMPC_bmy0);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_20(beliefPenaltyMPC_Ld0, beliefPenaltyMPC_bmy0, beliefPenaltyMPC_dvcc0);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_20_62_20(params->C1, beliefPenaltyMPC_dvcc1, beliefPenaltyMPC_D0, beliefPenaltyMPC_dvcc0, beliefPenaltyMPC_grad_eq0);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_20_62_20(params->C2, beliefPenaltyMPC_dvcc2, beliefPenaltyMPC_D1, beliefPenaltyMPC_dvcc1, beliefPenaltyMPC_grad_eq1);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_20_62_20(params->C3, beliefPenaltyMPC_dvcc3, beliefPenaltyMPC_D1, beliefPenaltyMPC_dvcc2, beliefPenaltyMPC_grad_eq2);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_20_62_20(params->C4, beliefPenaltyMPC_dvcc4, beliefPenaltyMPC_D1, beliefPenaltyMPC_dvcc3, beliefPenaltyMPC_grad_eq3);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_20_62_20(params->C5, beliefPenaltyMPC_dvcc5, beliefPenaltyMPC_D1, beliefPenaltyMPC_dvcc4, beliefPenaltyMPC_grad_eq4);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_20_62_20(params->C6, beliefPenaltyMPC_dvcc6, beliefPenaltyMPC_D1, beliefPenaltyMPC_dvcc5, beliefPenaltyMPC_grad_eq5);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_20_62_20(params->C7, beliefPenaltyMPC_dvcc7, beliefPenaltyMPC_D1, beliefPenaltyMPC_dvcc6, beliefPenaltyMPC_grad_eq6);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_20_62_20(params->C8, beliefPenaltyMPC_dvcc8, beliefPenaltyMPC_D1, beliefPenaltyMPC_dvcc7, beliefPenaltyMPC_grad_eq7);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_20_62_20(params->C9, beliefPenaltyMPC_dvcc9, beliefPenaltyMPC_D1, beliefPenaltyMPC_dvcc8, beliefPenaltyMPC_grad_eq8);
beliefPenaltyMPC_LA_DIAGZERO_MTVM_20_20(beliefPenaltyMPC_D9, beliefPenaltyMPC_dvcc9, beliefPenaltyMPC_grad_eq9);
beliefPenaltyMPC_LA_VSUB_578(beliefPenaltyMPC_rd, beliefPenaltyMPC_grad_eq, beliefPenaltyMPC_rd);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_62(beliefPenaltyMPC_Phi0, beliefPenaltyMPC_rd0, beliefPenaltyMPC_dzcc0);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_62(beliefPenaltyMPC_Phi1, beliefPenaltyMPC_rd1, beliefPenaltyMPC_dzcc1);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_62(beliefPenaltyMPC_Phi2, beliefPenaltyMPC_rd2, beliefPenaltyMPC_dzcc2);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_62(beliefPenaltyMPC_Phi3, beliefPenaltyMPC_rd3, beliefPenaltyMPC_dzcc3);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_62(beliefPenaltyMPC_Phi4, beliefPenaltyMPC_rd4, beliefPenaltyMPC_dzcc4);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_62(beliefPenaltyMPC_Phi5, beliefPenaltyMPC_rd5, beliefPenaltyMPC_dzcc5);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_62(beliefPenaltyMPC_Phi6, beliefPenaltyMPC_rd6, beliefPenaltyMPC_dzcc6);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_62(beliefPenaltyMPC_Phi7, beliefPenaltyMPC_rd7, beliefPenaltyMPC_dzcc7);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_62(beliefPenaltyMPC_Phi8, beliefPenaltyMPC_rd8, beliefPenaltyMPC_dzcc8);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_20(beliefPenaltyMPC_Phi9, beliefPenaltyMPC_rd9, beliefPenaltyMPC_dzcc9);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_62(beliefPenaltyMPC_ccrhsl0, beliefPenaltyMPC_slb0, beliefPenaltyMPC_llbbyslb0, beliefPenaltyMPC_dzcc0, beliefPenaltyMPC_lbIdx0, beliefPenaltyMPC_dllbcc0);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_22(beliefPenaltyMPC_ccrhsub0, beliefPenaltyMPC_sub0, beliefPenaltyMPC_lubbysub0, beliefPenaltyMPC_dzcc0, beliefPenaltyMPC_ubIdx0, beliefPenaltyMPC_dlubcc0);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_62(beliefPenaltyMPC_ccrhsl1, beliefPenaltyMPC_slb1, beliefPenaltyMPC_llbbyslb1, beliefPenaltyMPC_dzcc1, beliefPenaltyMPC_lbIdx1, beliefPenaltyMPC_dllbcc1);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_22(beliefPenaltyMPC_ccrhsub1, beliefPenaltyMPC_sub1, beliefPenaltyMPC_lubbysub1, beliefPenaltyMPC_dzcc1, beliefPenaltyMPC_ubIdx1, beliefPenaltyMPC_dlubcc1);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_62(beliefPenaltyMPC_ccrhsl2, beliefPenaltyMPC_slb2, beliefPenaltyMPC_llbbyslb2, beliefPenaltyMPC_dzcc2, beliefPenaltyMPC_lbIdx2, beliefPenaltyMPC_dllbcc2);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_22(beliefPenaltyMPC_ccrhsub2, beliefPenaltyMPC_sub2, beliefPenaltyMPC_lubbysub2, beliefPenaltyMPC_dzcc2, beliefPenaltyMPC_ubIdx2, beliefPenaltyMPC_dlubcc2);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_62(beliefPenaltyMPC_ccrhsl3, beliefPenaltyMPC_slb3, beliefPenaltyMPC_llbbyslb3, beliefPenaltyMPC_dzcc3, beliefPenaltyMPC_lbIdx3, beliefPenaltyMPC_dllbcc3);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_22(beliefPenaltyMPC_ccrhsub3, beliefPenaltyMPC_sub3, beliefPenaltyMPC_lubbysub3, beliefPenaltyMPC_dzcc3, beliefPenaltyMPC_ubIdx3, beliefPenaltyMPC_dlubcc3);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_62(beliefPenaltyMPC_ccrhsl4, beliefPenaltyMPC_slb4, beliefPenaltyMPC_llbbyslb4, beliefPenaltyMPC_dzcc4, beliefPenaltyMPC_lbIdx4, beliefPenaltyMPC_dllbcc4);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_22(beliefPenaltyMPC_ccrhsub4, beliefPenaltyMPC_sub4, beliefPenaltyMPC_lubbysub4, beliefPenaltyMPC_dzcc4, beliefPenaltyMPC_ubIdx4, beliefPenaltyMPC_dlubcc4);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_62(beliefPenaltyMPC_ccrhsl5, beliefPenaltyMPC_slb5, beliefPenaltyMPC_llbbyslb5, beliefPenaltyMPC_dzcc5, beliefPenaltyMPC_lbIdx5, beliefPenaltyMPC_dllbcc5);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_22(beliefPenaltyMPC_ccrhsub5, beliefPenaltyMPC_sub5, beliefPenaltyMPC_lubbysub5, beliefPenaltyMPC_dzcc5, beliefPenaltyMPC_ubIdx5, beliefPenaltyMPC_dlubcc5);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_62(beliefPenaltyMPC_ccrhsl6, beliefPenaltyMPC_slb6, beliefPenaltyMPC_llbbyslb6, beliefPenaltyMPC_dzcc6, beliefPenaltyMPC_lbIdx6, beliefPenaltyMPC_dllbcc6);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_22(beliefPenaltyMPC_ccrhsub6, beliefPenaltyMPC_sub6, beliefPenaltyMPC_lubbysub6, beliefPenaltyMPC_dzcc6, beliefPenaltyMPC_ubIdx6, beliefPenaltyMPC_dlubcc6);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_62(beliefPenaltyMPC_ccrhsl7, beliefPenaltyMPC_slb7, beliefPenaltyMPC_llbbyslb7, beliefPenaltyMPC_dzcc7, beliefPenaltyMPC_lbIdx7, beliefPenaltyMPC_dllbcc7);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_22(beliefPenaltyMPC_ccrhsub7, beliefPenaltyMPC_sub7, beliefPenaltyMPC_lubbysub7, beliefPenaltyMPC_dzcc7, beliefPenaltyMPC_ubIdx7, beliefPenaltyMPC_dlubcc7);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_62(beliefPenaltyMPC_ccrhsl8, beliefPenaltyMPC_slb8, beliefPenaltyMPC_llbbyslb8, beliefPenaltyMPC_dzcc8, beliefPenaltyMPC_lbIdx8, beliefPenaltyMPC_dllbcc8);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_22(beliefPenaltyMPC_ccrhsub8, beliefPenaltyMPC_sub8, beliefPenaltyMPC_lubbysub8, beliefPenaltyMPC_dzcc8, beliefPenaltyMPC_ubIdx8, beliefPenaltyMPC_dlubcc8);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_20(beliefPenaltyMPC_ccrhsl9, beliefPenaltyMPC_slb9, beliefPenaltyMPC_llbbyslb9, beliefPenaltyMPC_dzcc9, beliefPenaltyMPC_lbIdx9, beliefPenaltyMPC_dllbcc9);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_20(beliefPenaltyMPC_ccrhsub9, beliefPenaltyMPC_sub9, beliefPenaltyMPC_lubbysub9, beliefPenaltyMPC_dzcc9, beliefPenaltyMPC_ubIdx9, beliefPenaltyMPC_dlubcc9);
beliefPenaltyMPC_LA_VSUB7_796(beliefPenaltyMPC_l, beliefPenaltyMPC_ccrhs, beliefPenaltyMPC_s, beliefPenaltyMPC_dl_cc, beliefPenaltyMPC_ds_cc);
beliefPenaltyMPC_LA_VADD_578(beliefPenaltyMPC_dz_cc, beliefPenaltyMPC_dz_aff);
beliefPenaltyMPC_LA_VADD_200(beliefPenaltyMPC_dv_cc, beliefPenaltyMPC_dv_aff);
beliefPenaltyMPC_LA_VADD_796(beliefPenaltyMPC_dl_cc, beliefPenaltyMPC_dl_aff);
beliefPenaltyMPC_LA_VADD_796(beliefPenaltyMPC_ds_cc, beliefPenaltyMPC_ds_aff);
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
output->z1[7] = beliefPenaltyMPC_z0[7];
output->z1[8] = beliefPenaltyMPC_z0[8];
output->z1[9] = beliefPenaltyMPC_z0[9];
output->z1[10] = beliefPenaltyMPC_z0[10];
output->z1[11] = beliefPenaltyMPC_z0[11];
output->z1[12] = beliefPenaltyMPC_z0[12];
output->z1[13] = beliefPenaltyMPC_z0[13];
output->z1[14] = beliefPenaltyMPC_z0[14];
output->z1[15] = beliefPenaltyMPC_z0[15];
output->z1[16] = beliefPenaltyMPC_z0[16];
output->z1[17] = beliefPenaltyMPC_z0[17];
output->z1[18] = beliefPenaltyMPC_z0[18];
output->z1[19] = beliefPenaltyMPC_z0[19];
output->z1[20] = beliefPenaltyMPC_z0[20];
output->z1[21] = beliefPenaltyMPC_z0[21];
output->z2[0] = beliefPenaltyMPC_z1[0];
output->z2[1] = beliefPenaltyMPC_z1[1];
output->z2[2] = beliefPenaltyMPC_z1[2];
output->z2[3] = beliefPenaltyMPC_z1[3];
output->z2[4] = beliefPenaltyMPC_z1[4];
output->z2[5] = beliefPenaltyMPC_z1[5];
output->z2[6] = beliefPenaltyMPC_z1[6];
output->z2[7] = beliefPenaltyMPC_z1[7];
output->z2[8] = beliefPenaltyMPC_z1[8];
output->z2[9] = beliefPenaltyMPC_z1[9];
output->z2[10] = beliefPenaltyMPC_z1[10];
output->z2[11] = beliefPenaltyMPC_z1[11];
output->z2[12] = beliefPenaltyMPC_z1[12];
output->z2[13] = beliefPenaltyMPC_z1[13];
output->z2[14] = beliefPenaltyMPC_z1[14];
output->z2[15] = beliefPenaltyMPC_z1[15];
output->z2[16] = beliefPenaltyMPC_z1[16];
output->z2[17] = beliefPenaltyMPC_z1[17];
output->z2[18] = beliefPenaltyMPC_z1[18];
output->z2[19] = beliefPenaltyMPC_z1[19];
output->z2[20] = beliefPenaltyMPC_z1[20];
output->z2[21] = beliefPenaltyMPC_z1[21];
output->z3[0] = beliefPenaltyMPC_z2[0];
output->z3[1] = beliefPenaltyMPC_z2[1];
output->z3[2] = beliefPenaltyMPC_z2[2];
output->z3[3] = beliefPenaltyMPC_z2[3];
output->z3[4] = beliefPenaltyMPC_z2[4];
output->z3[5] = beliefPenaltyMPC_z2[5];
output->z3[6] = beliefPenaltyMPC_z2[6];
output->z3[7] = beliefPenaltyMPC_z2[7];
output->z3[8] = beliefPenaltyMPC_z2[8];
output->z3[9] = beliefPenaltyMPC_z2[9];
output->z3[10] = beliefPenaltyMPC_z2[10];
output->z3[11] = beliefPenaltyMPC_z2[11];
output->z3[12] = beliefPenaltyMPC_z2[12];
output->z3[13] = beliefPenaltyMPC_z2[13];
output->z3[14] = beliefPenaltyMPC_z2[14];
output->z3[15] = beliefPenaltyMPC_z2[15];
output->z3[16] = beliefPenaltyMPC_z2[16];
output->z3[17] = beliefPenaltyMPC_z2[17];
output->z3[18] = beliefPenaltyMPC_z2[18];
output->z3[19] = beliefPenaltyMPC_z2[19];
output->z3[20] = beliefPenaltyMPC_z2[20];
output->z3[21] = beliefPenaltyMPC_z2[21];
output->z4[0] = beliefPenaltyMPC_z3[0];
output->z4[1] = beliefPenaltyMPC_z3[1];
output->z4[2] = beliefPenaltyMPC_z3[2];
output->z4[3] = beliefPenaltyMPC_z3[3];
output->z4[4] = beliefPenaltyMPC_z3[4];
output->z4[5] = beliefPenaltyMPC_z3[5];
output->z4[6] = beliefPenaltyMPC_z3[6];
output->z4[7] = beliefPenaltyMPC_z3[7];
output->z4[8] = beliefPenaltyMPC_z3[8];
output->z4[9] = beliefPenaltyMPC_z3[9];
output->z4[10] = beliefPenaltyMPC_z3[10];
output->z4[11] = beliefPenaltyMPC_z3[11];
output->z4[12] = beliefPenaltyMPC_z3[12];
output->z4[13] = beliefPenaltyMPC_z3[13];
output->z4[14] = beliefPenaltyMPC_z3[14];
output->z4[15] = beliefPenaltyMPC_z3[15];
output->z4[16] = beliefPenaltyMPC_z3[16];
output->z4[17] = beliefPenaltyMPC_z3[17];
output->z4[18] = beliefPenaltyMPC_z3[18];
output->z4[19] = beliefPenaltyMPC_z3[19];
output->z4[20] = beliefPenaltyMPC_z3[20];
output->z4[21] = beliefPenaltyMPC_z3[21];
output->z5[0] = beliefPenaltyMPC_z4[0];
output->z5[1] = beliefPenaltyMPC_z4[1];
output->z5[2] = beliefPenaltyMPC_z4[2];
output->z5[3] = beliefPenaltyMPC_z4[3];
output->z5[4] = beliefPenaltyMPC_z4[4];
output->z5[5] = beliefPenaltyMPC_z4[5];
output->z5[6] = beliefPenaltyMPC_z4[6];
output->z5[7] = beliefPenaltyMPC_z4[7];
output->z5[8] = beliefPenaltyMPC_z4[8];
output->z5[9] = beliefPenaltyMPC_z4[9];
output->z5[10] = beliefPenaltyMPC_z4[10];
output->z5[11] = beliefPenaltyMPC_z4[11];
output->z5[12] = beliefPenaltyMPC_z4[12];
output->z5[13] = beliefPenaltyMPC_z4[13];
output->z5[14] = beliefPenaltyMPC_z4[14];
output->z5[15] = beliefPenaltyMPC_z4[15];
output->z5[16] = beliefPenaltyMPC_z4[16];
output->z5[17] = beliefPenaltyMPC_z4[17];
output->z5[18] = beliefPenaltyMPC_z4[18];
output->z5[19] = beliefPenaltyMPC_z4[19];
output->z5[20] = beliefPenaltyMPC_z4[20];
output->z5[21] = beliefPenaltyMPC_z4[21];
output->z6[0] = beliefPenaltyMPC_z5[0];
output->z6[1] = beliefPenaltyMPC_z5[1];
output->z6[2] = beliefPenaltyMPC_z5[2];
output->z6[3] = beliefPenaltyMPC_z5[3];
output->z6[4] = beliefPenaltyMPC_z5[4];
output->z6[5] = beliefPenaltyMPC_z5[5];
output->z6[6] = beliefPenaltyMPC_z5[6];
output->z6[7] = beliefPenaltyMPC_z5[7];
output->z6[8] = beliefPenaltyMPC_z5[8];
output->z6[9] = beliefPenaltyMPC_z5[9];
output->z6[10] = beliefPenaltyMPC_z5[10];
output->z6[11] = beliefPenaltyMPC_z5[11];
output->z6[12] = beliefPenaltyMPC_z5[12];
output->z6[13] = beliefPenaltyMPC_z5[13];
output->z6[14] = beliefPenaltyMPC_z5[14];
output->z6[15] = beliefPenaltyMPC_z5[15];
output->z6[16] = beliefPenaltyMPC_z5[16];
output->z6[17] = beliefPenaltyMPC_z5[17];
output->z6[18] = beliefPenaltyMPC_z5[18];
output->z6[19] = beliefPenaltyMPC_z5[19];
output->z6[20] = beliefPenaltyMPC_z5[20];
output->z6[21] = beliefPenaltyMPC_z5[21];
output->z7[0] = beliefPenaltyMPC_z6[0];
output->z7[1] = beliefPenaltyMPC_z6[1];
output->z7[2] = beliefPenaltyMPC_z6[2];
output->z7[3] = beliefPenaltyMPC_z6[3];
output->z7[4] = beliefPenaltyMPC_z6[4];
output->z7[5] = beliefPenaltyMPC_z6[5];
output->z7[6] = beliefPenaltyMPC_z6[6];
output->z7[7] = beliefPenaltyMPC_z6[7];
output->z7[8] = beliefPenaltyMPC_z6[8];
output->z7[9] = beliefPenaltyMPC_z6[9];
output->z7[10] = beliefPenaltyMPC_z6[10];
output->z7[11] = beliefPenaltyMPC_z6[11];
output->z7[12] = beliefPenaltyMPC_z6[12];
output->z7[13] = beliefPenaltyMPC_z6[13];
output->z7[14] = beliefPenaltyMPC_z6[14];
output->z7[15] = beliefPenaltyMPC_z6[15];
output->z7[16] = beliefPenaltyMPC_z6[16];
output->z7[17] = beliefPenaltyMPC_z6[17];
output->z7[18] = beliefPenaltyMPC_z6[18];
output->z7[19] = beliefPenaltyMPC_z6[19];
output->z7[20] = beliefPenaltyMPC_z6[20];
output->z7[21] = beliefPenaltyMPC_z6[21];
output->z8[0] = beliefPenaltyMPC_z7[0];
output->z8[1] = beliefPenaltyMPC_z7[1];
output->z8[2] = beliefPenaltyMPC_z7[2];
output->z8[3] = beliefPenaltyMPC_z7[3];
output->z8[4] = beliefPenaltyMPC_z7[4];
output->z8[5] = beliefPenaltyMPC_z7[5];
output->z8[6] = beliefPenaltyMPC_z7[6];
output->z8[7] = beliefPenaltyMPC_z7[7];
output->z8[8] = beliefPenaltyMPC_z7[8];
output->z8[9] = beliefPenaltyMPC_z7[9];
output->z8[10] = beliefPenaltyMPC_z7[10];
output->z8[11] = beliefPenaltyMPC_z7[11];
output->z8[12] = beliefPenaltyMPC_z7[12];
output->z8[13] = beliefPenaltyMPC_z7[13];
output->z8[14] = beliefPenaltyMPC_z7[14];
output->z8[15] = beliefPenaltyMPC_z7[15];
output->z8[16] = beliefPenaltyMPC_z7[16];
output->z8[17] = beliefPenaltyMPC_z7[17];
output->z8[18] = beliefPenaltyMPC_z7[18];
output->z8[19] = beliefPenaltyMPC_z7[19];
output->z8[20] = beliefPenaltyMPC_z7[20];
output->z8[21] = beliefPenaltyMPC_z7[21];
output->z9[0] = beliefPenaltyMPC_z8[0];
output->z9[1] = beliefPenaltyMPC_z8[1];
output->z9[2] = beliefPenaltyMPC_z8[2];
output->z9[3] = beliefPenaltyMPC_z8[3];
output->z9[4] = beliefPenaltyMPC_z8[4];
output->z9[5] = beliefPenaltyMPC_z8[5];
output->z9[6] = beliefPenaltyMPC_z8[6];
output->z9[7] = beliefPenaltyMPC_z8[7];
output->z9[8] = beliefPenaltyMPC_z8[8];
output->z9[9] = beliefPenaltyMPC_z8[9];
output->z9[10] = beliefPenaltyMPC_z8[10];
output->z9[11] = beliefPenaltyMPC_z8[11];
output->z9[12] = beliefPenaltyMPC_z8[12];
output->z9[13] = beliefPenaltyMPC_z8[13];
output->z9[14] = beliefPenaltyMPC_z8[14];
output->z9[15] = beliefPenaltyMPC_z8[15];
output->z9[16] = beliefPenaltyMPC_z8[16];
output->z9[17] = beliefPenaltyMPC_z8[17];
output->z9[18] = beliefPenaltyMPC_z8[18];
output->z9[19] = beliefPenaltyMPC_z8[19];
output->z9[20] = beliefPenaltyMPC_z8[20];
output->z9[21] = beliefPenaltyMPC_z8[21];
output->z10[0] = beliefPenaltyMPC_z9[0];
output->z10[1] = beliefPenaltyMPC_z9[1];
output->z10[2] = beliefPenaltyMPC_z9[2];
output->z10[3] = beliefPenaltyMPC_z9[3];
output->z10[4] = beliefPenaltyMPC_z9[4];
output->z10[5] = beliefPenaltyMPC_z9[5];
output->z10[6] = beliefPenaltyMPC_z9[6];
output->z10[7] = beliefPenaltyMPC_z9[7];
output->z10[8] = beliefPenaltyMPC_z9[8];
output->z10[9] = beliefPenaltyMPC_z9[9];
output->z10[10] = beliefPenaltyMPC_z9[10];
output->z10[11] = beliefPenaltyMPC_z9[11];
output->z10[12] = beliefPenaltyMPC_z9[12];
output->z10[13] = beliefPenaltyMPC_z9[13];
output->z10[14] = beliefPenaltyMPC_z9[14];
output->z10[15] = beliefPenaltyMPC_z9[15];
output->z10[16] = beliefPenaltyMPC_z9[16];
output->z10[17] = beliefPenaltyMPC_z9[17];
output->z10[18] = beliefPenaltyMPC_z9[18];
output->z10[19] = beliefPenaltyMPC_z9[19];

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
