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
 * Initializes a vector of length 3339 with a value.
 */
void beliefPenaltyMPC_LA_INITIALIZEVECTOR_3339(beliefPenaltyMPC_FLOAT* vec, beliefPenaltyMPC_FLOAT value)
{
	int i;
	for( i=0; i<3339; i++ )
	{
		vec[i] = value;
	}
}


/*
 * Initializes a vector of length 1155 with a value.
 */
void beliefPenaltyMPC_LA_INITIALIZEVECTOR_1155(beliefPenaltyMPC_FLOAT* vec, beliefPenaltyMPC_FLOAT value)
{
	int i;
	for( i=0; i<1155; i++ )
	{
		vec[i] = value;
	}
}


/*
 * Initializes a vector of length 4522 with a value.
 */
void beliefPenaltyMPC_LA_INITIALIZEVECTOR_4522(beliefPenaltyMPC_FLOAT* vec, beliefPenaltyMPC_FLOAT value)
{
	int i;
	for( i=0; i<4522; i++ )
	{
		vec[i] = value;
	}
}


/* 
 * Calculates a dot product and adds it to a variable: z += x'*y; 
 * This function is for vectors of length 4522.
 */
void beliefPenaltyMPC_LA_DOTACC_4522(beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *y, beliefPenaltyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<4522; i++ ){
		*z += x[i]*y[i];
	}
}


/*
 * Calculates the gradient and the value for a quadratic function 0.5*z'*H*z + f'*z
 *
 * INPUTS:     H  - Symmetric Hessian, diag matrix of size [233 x 233]
 *             f  - column vector of size 233
 *             z  - column vector of size 233
 *
 * OUTPUTS: grad  - gradient at z (= H*z + f), column vector of size 233
 *          value <-- value + 0.5*z'*H*z + f'*z (value will be modified)
 */
void beliefPenaltyMPC_LA_DIAG_QUADFCN_233(beliefPenaltyMPC_FLOAT* H, beliefPenaltyMPC_FLOAT* f, beliefPenaltyMPC_FLOAT* z, beliefPenaltyMPC_FLOAT* grad, beliefPenaltyMPC_FLOAT* value)
{
	int i;
	beliefPenaltyMPC_FLOAT hz;	
	for( i=0; i<233; i++){
		hz = H[i]*z[i];
		grad[i] = hz + f[i];
		*value += 0.5*hz*z[i] + f[i]*z[i];
	}
}


/*
 * Calculates the gradient and the value for a quadratic function 0.5*z'*H*z + f'*z
 *
 * INPUTS:     H  - Symmetric Hessian, diag matrix of size [77 x 77]
 *             f  - column vector of size 77
 *             z  - column vector of size 77
 *
 * OUTPUTS: grad  - gradient at z (= H*z + f), column vector of size 77
 *          value <-- value + 0.5*z'*H*z + f'*z (value will be modified)
 */
void beliefPenaltyMPC_LA_DIAG_QUADFCN_77(beliefPenaltyMPC_FLOAT* H, beliefPenaltyMPC_FLOAT* f, beliefPenaltyMPC_FLOAT* z, beliefPenaltyMPC_FLOAT* grad, beliefPenaltyMPC_FLOAT* value)
{
	int i;
	beliefPenaltyMPC_FLOAT hz;	
	for( i=0; i<77; i++){
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
void beliefPenaltyMPC_LA_DIAGZERO_MVMSUB6_77(beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *b, beliefPenaltyMPC_FLOAT *l, beliefPenaltyMPC_FLOAT *r, beliefPenaltyMPC_FLOAT *z, beliefPenaltyMPC_FLOAT *y)
{
	int i;
	beliefPenaltyMPC_FLOAT Bu[77];
	beliefPenaltyMPC_FLOAT norm = *y;
	beliefPenaltyMPC_FLOAT lr = 0;

	/* do A*x + B*u first */
	for( i=0; i<77; i++ ){
		Bu[i] = B[i]*u[i];
	}	

	for( i=0; i<77; i++ ){
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
void beliefPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_77_233_233(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *b, beliefPenaltyMPC_FLOAT *l, beliefPenaltyMPC_FLOAT *r, beliefPenaltyMPC_FLOAT *z, beliefPenaltyMPC_FLOAT *y)
{
	int i;
	int j;
	int k = 0;
	beliefPenaltyMPC_FLOAT AxBu[77];
	beliefPenaltyMPC_FLOAT norm = *y;
	beliefPenaltyMPC_FLOAT lr = 0;

	/* do A*x + B*u first */
	for( i=0; i<77; i++ ){
		AxBu[i] = A[k++]*x[0] + B[i]*u[i];
	}	

	for( j=1; j<233; j++ ){		
		for( i=0; i<77; i++ ){
			AxBu[i] += A[k++]*x[j];
		}
	}

	for( i=0; i<77; i++ ){
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
void beliefPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_77_233_77(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *b, beliefPenaltyMPC_FLOAT *l, beliefPenaltyMPC_FLOAT *r, beliefPenaltyMPC_FLOAT *z, beliefPenaltyMPC_FLOAT *y)
{
	int i;
	int j;
	int k = 0;
	beliefPenaltyMPC_FLOAT AxBu[77];
	beliefPenaltyMPC_FLOAT norm = *y;
	beliefPenaltyMPC_FLOAT lr = 0;

	/* do A*x + B*u first */
	for( i=0; i<77; i++ ){
		AxBu[i] = A[k++]*x[0] + B[i]*u[i];
	}	

	for( j=1; j<233; j++ ){		
		for( i=0; i<77; i++ ){
			AxBu[i] += A[k++]*x[j];
		}
	}

	for( i=0; i<77; i++ ){
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
 * where A is of size [77 x 233] and stored in column major format.
 * and B is of size [77 x 233] and stored in diagzero format
 * Note the transposes of A and B!
 */
void beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_77_233_77(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *y, beliefPenaltyMPC_FLOAT *z)
{
	int i;
	int j;
	int k = 0;
	for( i=0; i<77; i++ ){
		z[i] = 0;
		for( j=0; j<77; j++ ){
			z[i] += A[k++]*x[j];
		}
		z[i] += B[i]*y[i];
	}
	for( i=77 ;i<233; i++ ){
		z[i] = 0;
		for( j=0; j<77; j++ ){
			z[i] += A[k++]*x[j];
		}
	}
}


/*
 * Matrix vector multiplication y = M'*x where M is of size [77 x 77]
 * and stored in diagzero format. Note the transpose of M!
 */
void beliefPenaltyMPC_LA_DIAGZERO_MTVM_77_77(beliefPenaltyMPC_FLOAT *M, beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *y)
{
	int i;
	for( i=0; i<77; i++ ){
		y[i] = M[i]*x[i];
	}
}


/*
 * Vector subtraction and addition.
 *	 Input: five vectors t, tidx, u, v, w and two scalars z and r
 *	 Output: y = t(tidx) - u + w
 *           z = z - v'*x;
 *           r = max([norm(y,inf), z]);
 * for vectors of length 233. Output z is of course scalar.
 */
void beliefPenaltyMPC_LA_VSUBADD3_233(beliefPenaltyMPC_FLOAT* t, beliefPenaltyMPC_FLOAT* u, int* uidx, beliefPenaltyMPC_FLOAT* v, beliefPenaltyMPC_FLOAT* w, beliefPenaltyMPC_FLOAT* y, beliefPenaltyMPC_FLOAT* z, beliefPenaltyMPC_FLOAT* r)
{
	int i;
	beliefPenaltyMPC_FLOAT norm = *r;
	beliefPenaltyMPC_FLOAT vx = 0;
	beliefPenaltyMPC_FLOAT x;
	for( i=0; i<233; i++){
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
 * for vectors of length 79. Output z is of course scalar.
 */
void beliefPenaltyMPC_LA_VSUBADD2_79(beliefPenaltyMPC_FLOAT* t, int* tidx, beliefPenaltyMPC_FLOAT* u, beliefPenaltyMPC_FLOAT* v, beliefPenaltyMPC_FLOAT* w, beliefPenaltyMPC_FLOAT* y, beliefPenaltyMPC_FLOAT* z, beliefPenaltyMPC_FLOAT* r)
{
	int i;
	beliefPenaltyMPC_FLOAT norm = *r;
	beliefPenaltyMPC_FLOAT vx = 0;
	beliefPenaltyMPC_FLOAT x;
	for( i=0; i<79; i++){
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
 * for vectors of length 77. Output z is of course scalar.
 */
void beliefPenaltyMPC_LA_VSUBADD3_77(beliefPenaltyMPC_FLOAT* t, beliefPenaltyMPC_FLOAT* u, int* uidx, beliefPenaltyMPC_FLOAT* v, beliefPenaltyMPC_FLOAT* w, beliefPenaltyMPC_FLOAT* y, beliefPenaltyMPC_FLOAT* z, beliefPenaltyMPC_FLOAT* r)
{
	int i;
	beliefPenaltyMPC_FLOAT norm = *r;
	beliefPenaltyMPC_FLOAT vx = 0;
	beliefPenaltyMPC_FLOAT x;
	for( i=0; i<77; i++){
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
 * for vectors of length 77. Output z is of course scalar.
 */
void beliefPenaltyMPC_LA_VSUBADD2_77(beliefPenaltyMPC_FLOAT* t, int* tidx, beliefPenaltyMPC_FLOAT* u, beliefPenaltyMPC_FLOAT* v, beliefPenaltyMPC_FLOAT* w, beliefPenaltyMPC_FLOAT* y, beliefPenaltyMPC_FLOAT* z, beliefPenaltyMPC_FLOAT* r)
{
	int i;
	beliefPenaltyMPC_FLOAT norm = *r;
	beliefPenaltyMPC_FLOAT vx = 0;
	beliefPenaltyMPC_FLOAT x;
	for( i=0; i<77; i++){
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
 * Special function for box constraints of length 233
 * Returns also L/S, a value that is often used elsewhere.
 */
void beliefPenaltyMPC_LA_INEQ_B_GRAD_233_233_79(beliefPenaltyMPC_FLOAT *lu, beliefPenaltyMPC_FLOAT *su, beliefPenaltyMPC_FLOAT *ru, beliefPenaltyMPC_FLOAT *ll, beliefPenaltyMPC_FLOAT *sl, beliefPenaltyMPC_FLOAT *rl, int* lbIdx, int* ubIdx, beliefPenaltyMPC_FLOAT *grad, beliefPenaltyMPC_FLOAT *lubysu, beliefPenaltyMPC_FLOAT *llbysl)
{
	int i;
	for( i=0; i<233; i++ ){
		grad[i] = 0;
	}
	for( i=0; i<233; i++ ){		
		llbysl[i] = ll[i] / sl[i];
		grad[lbIdx[i]] -= llbysl[i]*rl[i];
	}
	for( i=0; i<79; i++ ){
		lubysu[i] = lu[i] / su[i];
		grad[ubIdx[i]] += lubysu[i]*ru[i];
	}
}


/*
 * Computes inequality constraints gradient-
 * Special function for box constraints of length 77
 * Returns also L/S, a value that is often used elsewhere.
 */
void beliefPenaltyMPC_LA_INEQ_B_GRAD_77_77_77(beliefPenaltyMPC_FLOAT *lu, beliefPenaltyMPC_FLOAT *su, beliefPenaltyMPC_FLOAT *ru, beliefPenaltyMPC_FLOAT *ll, beliefPenaltyMPC_FLOAT *sl, beliefPenaltyMPC_FLOAT *rl, int* lbIdx, int* ubIdx, beliefPenaltyMPC_FLOAT *grad, beliefPenaltyMPC_FLOAT *lubysu, beliefPenaltyMPC_FLOAT *llbysl)
{
	int i;
	for( i=0; i<77; i++ ){
		grad[i] = 0;
	}
	for( i=0; i<77; i++ ){		
		llbysl[i] = ll[i] / sl[i];
		grad[lbIdx[i]] -= llbysl[i]*rl[i];
	}
	for( i=0; i<77; i++ ){
		lubysu[i] = lu[i] / su[i];
		grad[ubIdx[i]] += lubysu[i]*ru[i];
	}
}


/*
 * Addition of three vectors  z = u + w + v
 * of length 3339.
 */
void beliefPenaltyMPC_LA_VVADD3_3339(beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *v, beliefPenaltyMPC_FLOAT *w, beliefPenaltyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<3339; i++ ){
		z[i] = u[i] + v[i] + w[i];
	}
}


/*
 * Special function to compute the diagonal cholesky factorization of the 
 * positive definite augmented Hessian for block size 233.
 *
 * Inputs: - H = diagonal cost Hessian in diagonal storage format
 *         - llbysl = L / S of lower bounds
 *         - lubysu = L / S of upper bounds
 *
 * Output: Phi = sqrt(H + diag(llbysl) + diag(lubysu))
 * where Phi is stored in diagonal storage format
 */
void beliefPenaltyMPC_LA_DIAG_CHOL_LBUB_233_233_79(beliefPenaltyMPC_FLOAT *H, beliefPenaltyMPC_FLOAT *llbysl, int* lbIdx, beliefPenaltyMPC_FLOAT *lubysu, int* ubIdx, beliefPenaltyMPC_FLOAT *Phi)


{
	int i;
	
	/* copy  H into PHI */
	for( i=0; i<233; i++ ){
		Phi[i] = H[i];
	}

	/* add llbysl onto Phi where necessary */
	for( i=0; i<233; i++ ){
		Phi[lbIdx[i]] += llbysl[i];
	}

	/* add lubysu onto Phi where necessary */
	for( i=0; i<79; i++){
		Phi[ubIdx[i]] +=  lubysu[i];
	}
	
	/* compute cholesky */
	for(i=0; i<233; i++)
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
 * where A is to be computed and is of size [77 x 233],
 * B is given and of size [77 x 233], L is a diagonal
 * matrix of size 77 stored in diagonal matrix 
 * storage format. Note the transpose of L has no impact!
 *
 * Result: A in column major storage format.
 *
 */
void beliefPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_77_233(beliefPenaltyMPC_FLOAT *L, beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *A)
{
    int i,j;
	 int k = 0;

	for( j=0; j<233; j++){
		for( i=0; i<77; i++){
			A[k] = B[k]/L[j];
			k++;
		}
	}

}


/**
 * Forward substitution for the matrix equation A*L' = B
 * where A is to be computed and is of size [77 x 233],
 * B is given and of size [77 x 233], L is a diagonal
 *  matrix of size 233 stored in diagonal 
 * storage format. Note the transpose of L!
 *
 * Result: A in diagonalzero storage format.
 *
 */
void beliefPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_77_233(beliefPenaltyMPC_FLOAT *L, beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *A)
{
	int j;
    for( j=0; j<233; j++ ){   
		A[j] = B[j]/L[j];
     }
}


/**
 * Compute C = A*B' where 
 *
 *	size(A) = [77 x 233]
 *  size(B) = [77 x 233] in diagzero format
 * 
 * A and C matrices are stored in column major format.
 * 
 * 
 */
void beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMTM_77_233_77(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *C)
{
    int i, j;
	
	for( i=0; i<77; i++ ){
		for( j=0; j<77; j++){
			C[j*77+i] = B[i*77+j]*A[i];
		}
	}

}


/**
 * Forward substitution to solve L*y = b where L is a
 * diagonal matrix in vector storage format.
 * 
 * The dimensions involved are 233.
 */
void beliefPenaltyMPC_LA_DIAG_FORWARDSUB_233(beliefPenaltyMPC_FLOAT *L, beliefPenaltyMPC_FLOAT *b, beliefPenaltyMPC_FLOAT *y)
{
    int i;

    for( i=0; i<233; i++ ){
		y[i] = b[i]/L[i];
    }
}


/*
 * Special function to compute the diagonal cholesky factorization of the 
 * positive definite augmented Hessian for block size 77.
 *
 * Inputs: - H = diagonal cost Hessian in diagonal storage format
 *         - llbysl = L / S of lower bounds
 *         - lubysu = L / S of upper bounds
 *
 * Output: Phi = sqrt(H + diag(llbysl) + diag(lubysu))
 * where Phi is stored in diagonal storage format
 */
void beliefPenaltyMPC_LA_DIAG_CHOL_ONELOOP_LBUB_77_77_77(beliefPenaltyMPC_FLOAT *H, beliefPenaltyMPC_FLOAT *llbysl, int* lbIdx, beliefPenaltyMPC_FLOAT *lubysu, int* ubIdx, beliefPenaltyMPC_FLOAT *Phi)


{
	int i;
	
	/* compute cholesky */
	for( i=0; i<77; i++ ){
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
 * where A is to be computed and is of size [77 x 77],
 * B is given and of size [77 x 77], L is a diagonal
 *  matrix of size 77 stored in diagonal 
 * storage format. Note the transpose of L!
 *
 * Result: A in diagonalzero storage format.
 *
 */
void beliefPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_77_77(beliefPenaltyMPC_FLOAT *L, beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *A)
{
	int j;
    for( j=0; j<77; j++ ){   
		A[j] = B[j]/L[j];
     }
}


/**
 * Forward substitution to solve L*y = b where L is a
 * diagonal matrix in vector storage format.
 * 
 * The dimensions involved are 77.
 */
void beliefPenaltyMPC_LA_DIAG_FORWARDSUB_77(beliefPenaltyMPC_FLOAT *L, beliefPenaltyMPC_FLOAT *b, beliefPenaltyMPC_FLOAT *y)
{
    int i;

    for( i=0; i<77; i++ ){
		y[i] = b[i]/L[i];
    }
}


/**
 * Compute L = B*B', where L is lower triangular of size NXp1
 * and B is a diagzero matrix of size [77 x 233] in column
 * storage format.
 * 
 */
void beliefPenaltyMPC_LA_DIAGZERO_MMT_77(beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *L)
{
    int i, ii, di;
    
    ii = 0; di = 0;
    for( i=0; i<77; i++ ){        
		L[ii+i] = B[i]*B[i];
        ii += ++di;
    }
}


/* 
 * Computes r = b - B*u
 * B is stored in diagzero format
 */
void beliefPenaltyMPC_LA_DIAGZERO_MVMSUB7_77(beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *b, beliefPenaltyMPC_FLOAT *r)
{
	int i;

	for( i=0; i<77; i++ ){
		r[i] = b[i] - B[i]*u[i];
	}	
	
}


/**
 * Compute L = A*A' + B*B', where L is lower triangular of size NXp1
 * and A is a dense matrix of size [77 x 233] in column
 * storage format, and B is of size [77 x 233] diagonalzero
 * storage format.
 * 
 * THIS ONE HAS THE WORST ACCES PATTERN POSSIBLE. 
 * POSSIBKE FIX: PUT A AND B INTO ROW MAJOR FORMAT FIRST.
 * 
 */
void beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMT2_77_233_233(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *L)
{
    int i, j, k, ii, di;
    beliefPenaltyMPC_FLOAT ltemp;
    
    ii = 0; di = 0;
    for( i=0; i<77; i++ ){        
        for( j=0; j<=i; j++ ){
            ltemp = 0; 
            for( k=0; k<233; k++ ){
                ltemp += A[k*77+i]*A[k*77+j];
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
void beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_77_233_233(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *b, beliefPenaltyMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;

	for( i=0; i<77; i++ ){
		r[i] = b[i] - A[k++]*x[0] - B[i]*u[i];
	}	

	for( j=1; j<233; j++ ){		
		for( i=0; i<77; i++ ){
			r[i] -= A[k++]*x[j];
		}
	}
	
}


/**
 * Compute L = A*A' + B*B', where L is lower triangular of size NXp1
 * and A is a dense matrix of size [77 x 233] in column
 * storage format, and B is of size [77 x 77] diagonalzero
 * storage format.
 * 
 * THIS ONE HAS THE WORST ACCES PATTERN POSSIBLE. 
 * POSSIBKE FIX: PUT A AND B INTO ROW MAJOR FORMAT FIRST.
 * 
 */
void beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMT2_77_233_77(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *L)
{
    int i, j, k, ii, di;
    beliefPenaltyMPC_FLOAT ltemp;
    
    ii = 0; di = 0;
    for( i=0; i<77; i++ ){        
        for( j=0; j<=i; j++ ){
            ltemp = 0; 
            for( k=0; k<233; k++ ){
                ltemp += A[k*77+i]*A[k*77+j];
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
void beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_77_233_77(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *b, beliefPenaltyMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;

	for( i=0; i<77; i++ ){
		r[i] = b[i] - A[k++]*x[0] - B[i]*u[i];
	}	

	for( j=1; j<233; j++ ){		
		for( i=0; i<77; i++ ){
			r[i] -= A[k++]*x[j];
		}
	}
	
}


/**
 * Cholesky factorization as above, but working on a matrix in 
 * lower triangular storage format of size 77 and outputting
 * the Cholesky factor to matrix L in lower triangular format.
 */
void beliefPenaltyMPC_LA_DENSE_CHOL_77(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *L)
{
    int i, j, k, di, dj;
	 int ii, jj;

    beliefPenaltyMPC_FLOAT l;
    beliefPenaltyMPC_FLOAT Mii;

	/* copy A to L first and then operate on L */
	/* COULD BE OPTIMIZED */
	ii=0; di=0;
	for( i=0; i<77; i++ ){
		for( j=0; j<=i; j++ ){
			L[ii+j] = A[ii+j];
		}
		ii += ++di;
	}    
	
	/* factor L */
	ii=0; di=0;
    for( i=0; i<77; i++ ){
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
        for( j=i+1; j<77; j++ ){
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
 * The dimensions involved are 77.
 */
void beliefPenaltyMPC_LA_DENSE_FORWARDSUB_77(beliefPenaltyMPC_FLOAT *L, beliefPenaltyMPC_FLOAT *b, beliefPenaltyMPC_FLOAT *y)
{
    int i,j,ii,di;
    beliefPenaltyMPC_FLOAT yel;
            
    ii = 0; di = 0;
    for( i=0; i<77; i++ ){
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
 * where A is to be computed and is of size [77 x 77],
 * B is given and of size [77 x 77], L is a lower tri-
 * angular matrix of size 77 stored in lower triangular 
 * storage format. Note the transpose of L AND B!
 *
 * Result: A in column major storage format.
 *
 */
void beliefPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_77_77(beliefPenaltyMPC_FLOAT *L, beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *A)
{
    int i,j,k,ii,di;
    beliefPenaltyMPC_FLOAT a;
    
    ii=0; di=0;
    for( j=0; j<77; j++ ){        
        for( i=0; i<77; i++ ){
            a = B[i*77+j];
            for( k=0; k<j; k++ ){
                a -= A[k*77+i]*L[ii+k];
            }    

			/* saturate for numerical stability */
			a = MIN(a, BIGM);
			a = MAX(a, -BIGM); 

			A[j*77+i] = a/L[ii+j];			
        }
        ii += ++di;
    }
}


/**
 * Compute L = L - A*A', where L is lower triangular of size 77
 * and A is a dense matrix of size [77 x 77] in column
 * storage format.
 * 
 * THIS ONE HAS THE WORST ACCES PATTERN POSSIBLE. 
 * POSSIBKE FIX: PUT A INTO ROW MAJOR FORMAT FIRST.
 * 
 */
void beliefPenaltyMPC_LA_DENSE_MMTSUB_77_77(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *L)
{
    int i, j, k, ii, di;
    beliefPenaltyMPC_FLOAT ltemp;
    
    ii = 0; di = 0;
    for( i=0; i<77; i++ ){        
        for( j=0; j<=i; j++ ){
            ltemp = 0; 
            for( k=0; k<77; k++ ){
                ltemp += A[k*77+i]*A[k*77+j];
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
void beliefPenaltyMPC_LA_DENSE_MVMSUB1_77_77(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *b, beliefPenaltyMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;

	for( i=0; i<77; i++ ){
		r[i] = b[i] - A[k++]*x[0];
	}	
	for( j=1; j<77; j++ ){		
		for( i=0; i<77; i++ ){
			r[i] -= A[k++]*x[j];
		}
	}
}


/**
 * Backward Substitution to solve L^T*x = y where L is a
 * lower triangular matrix in triangular storage format.
 * 
 * All involved dimensions are 77.
 */
void beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_77(beliefPenaltyMPC_FLOAT *L, beliefPenaltyMPC_FLOAT *y, beliefPenaltyMPC_FLOAT *x)
{
    int i, ii, di, j, jj, dj;
    beliefPenaltyMPC_FLOAT xel;    
	int start = 2926;
    
    /* now solve L^T*x = y by backward substitution */
    ii = start; di = 76;
    for( i=76; i>=0; i-- ){        
        xel = y[i];        
        jj = start; dj = 76;
        for( j=76; j>i; j-- ){
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
 * Matrix vector multiplication y = b - M'*x where M is of size [77 x 77]
 * and stored in column major format. Note the transpose of M!
 */
void beliefPenaltyMPC_LA_DENSE_MTVMSUB_77_77(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *b, beliefPenaltyMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0; 
	for( i=0; i<77; i++ ){
		r[i] = b[i];
		for( j=0; j<77; j++ ){
			r[i] -= A[k++]*x[j];
		}
	}
}


/*
 * Vector subtraction z = -x - y for vectors of length 3339.
 */
void beliefPenaltyMPC_LA_VSUB2_3339(beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *y, beliefPenaltyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<3339; i++){
		z[i] = -x[i] - y[i];
	}
}


/**
 * Forward-Backward-Substitution to solve L*L^T*x = b where L is a
 * diagonal matrix of size 233 in vector
 * storage format.
 */
void beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_233(beliefPenaltyMPC_FLOAT *L, beliefPenaltyMPC_FLOAT *b, beliefPenaltyMPC_FLOAT *x)
{
    int i;
            
    /* solve Ly = b by forward and backward substitution */
    for( i=0; i<233; i++ ){
		x[i] = b[i]/(L[i]*L[i]);
    }
    
}


/**
 * Forward-Backward-Substitution to solve L*L^T*x = b where L is a
 * diagonal matrix of size 77 in vector
 * storage format.
 */
void beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_77(beliefPenaltyMPC_FLOAT *L, beliefPenaltyMPC_FLOAT *b, beliefPenaltyMPC_FLOAT *x)
{
    int i;
            
    /* solve Ly = b by forward and backward substitution */
    for( i=0; i<77; i++ ){
		x[i] = b[i]/(L[i]*L[i]);
    }
    
}


/*
 * Vector subtraction z = x(xidx) - y where y, z and xidx are of length 233,
 * and x has length 233 and is indexed through yidx.
 */
void beliefPenaltyMPC_LA_VSUB_INDEXED_233(beliefPenaltyMPC_FLOAT *x, int* xidx, beliefPenaltyMPC_FLOAT *y, beliefPenaltyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<233; i++){
		z[i] = x[xidx[i]] - y[i];
	}
}


/*
 * Vector subtraction x = -u.*v - w for vectors of length 233.
 */
void beliefPenaltyMPC_LA_VSUB3_233(beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *v, beliefPenaltyMPC_FLOAT *w, beliefPenaltyMPC_FLOAT *x)
{
	int i;
	for( i=0; i<233; i++){
		x[i] = -u[i]*v[i] - w[i];
	}
}


/*
 * Vector subtraction z = -x - y(yidx) where y is of length 233
 * and z, x and yidx are of length 79.
 */
void beliefPenaltyMPC_LA_VSUB2_INDEXED_79(beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *y, int* yidx, beliefPenaltyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<79; i++){
		z[i] = -x[i] - y[yidx[i]];
	}
}


/*
 * Vector subtraction x = -u.*v - w for vectors of length 79.
 */
void beliefPenaltyMPC_LA_VSUB3_79(beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *v, beliefPenaltyMPC_FLOAT *w, beliefPenaltyMPC_FLOAT *x)
{
	int i;
	for( i=0; i<79; i++){
		x[i] = -u[i]*v[i] - w[i];
	}
}


/*
 * Vector subtraction z = x(xidx) - y where y, z and xidx are of length 77,
 * and x has length 77 and is indexed through yidx.
 */
void beliefPenaltyMPC_LA_VSUB_INDEXED_77(beliefPenaltyMPC_FLOAT *x, int* xidx, beliefPenaltyMPC_FLOAT *y, beliefPenaltyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<77; i++){
		z[i] = x[xidx[i]] - y[i];
	}
}


/*
 * Vector subtraction x = -u.*v - w for vectors of length 77.
 */
void beliefPenaltyMPC_LA_VSUB3_77(beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *v, beliefPenaltyMPC_FLOAT *w, beliefPenaltyMPC_FLOAT *x)
{
	int i;
	for( i=0; i<77; i++){
		x[i] = -u[i]*v[i] - w[i];
	}
}


/*
 * Vector subtraction z = -x - y(yidx) where y is of length 77
 * and z, x and yidx are of length 77.
 */
void beliefPenaltyMPC_LA_VSUB2_INDEXED_77(beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *y, int* yidx, beliefPenaltyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<77; i++){
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
        for( i=0; i<4522; i++ ){
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
        if( i == 4522 ){
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
    *mu_aff = mymu / (beliefPenaltyMPC_FLOAT)4522;
    return lsIt;
}


/*
 * Vector subtraction x = u.*v - a where a is a scalar
*  and x,u,v are vectors of length 4522.
 */
void beliefPenaltyMPC_LA_VSUB5_4522(beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *v, beliefPenaltyMPC_FLOAT a, beliefPenaltyMPC_FLOAT *x)
{
	int i;
	for( i=0; i<4522; i++){
		x[i] = u[i]*v[i] - a;
	}
}


/*
 * Computes x=0; x(uidx) += u/su; x(vidx) -= v/sv where x is of length 233,
 * u, su, uidx are of length 79 and v, sv, vidx are of length 233.
 */
void beliefPenaltyMPC_LA_VSUB6_INDEXED_233_79_233(beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *su, int* uidx, beliefPenaltyMPC_FLOAT *v, beliefPenaltyMPC_FLOAT *sv, int* vidx, beliefPenaltyMPC_FLOAT *x)
{
	int i;
	for( i=0; i<233; i++ ){
		x[i] = 0;
	}
	for( i=0; i<79; i++){
		x[uidx[i]] += u[i]/su[i];
	}
	for( i=0; i<233; i++){
		x[vidx[i]] -= v[i]/sv[i];
	}
}


/* 
 * Computes r =  B*u
 * where B is stored in diagzero format
 */
void beliefPenaltyMPC_LA_DIAGZERO_MVM_77(beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *r)
{
	int i;

	for( i=0; i<77; i++ ){
		r[i] = B[i]*u[i];
	}	
	
}


/* 
 * Computes r = A*x + B*u
 * where A is stored in column major format
 * and B is stored in diagzero format
 */
void beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_77_233_233(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;

	for( i=0; i<77; i++ ){
		r[i] = A[k++]*x[0] + B[i]*u[i];
	}	

	for( j=1; j<233; j++ ){		
		for( i=0; i<77; i++ ){
			r[i] += A[k++]*x[j];
		}
	}
	
}


/*
 * Computes x=0; x(uidx) += u/su; x(vidx) -= v/sv where x is of length 77,
 * u, su, uidx are of length 77 and v, sv, vidx are of length 77.
 */
void beliefPenaltyMPC_LA_VSUB6_INDEXED_77_77_77(beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *su, int* uidx, beliefPenaltyMPC_FLOAT *v, beliefPenaltyMPC_FLOAT *sv, int* vidx, beliefPenaltyMPC_FLOAT *x)
{
	int i;
	for( i=0; i<77; i++ ){
		x[i] = 0;
	}
	for( i=0; i<77; i++){
		x[uidx[i]] += u[i]/su[i];
	}
	for( i=0; i<77; i++){
		x[vidx[i]] -= v[i]/sv[i];
	}
}


/* 
 * Computes r = A*x + B*u
 * where A is stored in column major format
 * and B is stored in diagzero format
 */
void beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_77_233_77(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;

	for( i=0; i<77; i++ ){
		r[i] = A[k++]*x[0] + B[i]*u[i];
	}	

	for( j=1; j<233; j++ ){		
		for( i=0; i<77; i++ ){
			r[i] += A[k++]*x[j];
		}
	}
	
}


/*
 * Vector subtraction z = x - y for vectors of length 3339.
 */
void beliefPenaltyMPC_LA_VSUB_3339(beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *y, beliefPenaltyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<3339; i++){
		z[i] = x[i] - y[i];
	}
}


/** 
 * Computes z = -r./s - u.*y(y)
 * where all vectors except of y are of length 233 (length of y >= 233).
 */
void beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_233(beliefPenaltyMPC_FLOAT *r, beliefPenaltyMPC_FLOAT *s, beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *y, int* yidx, beliefPenaltyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<233; i++ ){
		z[i] = -r[i]/s[i] - u[i]*y[yidx[i]];
	}
}


/** 
 * Computes z = -r./s + u.*y(y)
 * where all vectors except of y are of length 79 (length of y >= 79).
 */
void beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_79(beliefPenaltyMPC_FLOAT *r, beliefPenaltyMPC_FLOAT *s, beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *y, int* yidx, beliefPenaltyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<79; i++ ){
		z[i] = -r[i]/s[i] + u[i]*y[yidx[i]];
	}
}


/** 
 * Computes z = -r./s - u.*y(y)
 * where all vectors except of y are of length 77 (length of y >= 77).
 */
void beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_77(beliefPenaltyMPC_FLOAT *r, beliefPenaltyMPC_FLOAT *s, beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *y, int* yidx, beliefPenaltyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<77; i++ ){
		z[i] = -r[i]/s[i] - u[i]*y[yidx[i]];
	}
}


/** 
 * Computes z = -r./s + u.*y(y)
 * where all vectors except of y are of length 77 (length of y >= 77).
 */
void beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_77(beliefPenaltyMPC_FLOAT *r, beliefPenaltyMPC_FLOAT *s, beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *y, int* yidx, beliefPenaltyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<77; i++ ){
		z[i] = -r[i]/s[i] + u[i]*y[yidx[i]];
	}
}


/*
 * Computes ds = -l.\(r + s.*dl) for vectors of length 4522.
 */
void beliefPenaltyMPC_LA_VSUB7_4522(beliefPenaltyMPC_FLOAT *l, beliefPenaltyMPC_FLOAT *r, beliefPenaltyMPC_FLOAT *s, beliefPenaltyMPC_FLOAT *dl, beliefPenaltyMPC_FLOAT *ds)
{
	int i;
	for( i=0; i<4522; i++){
		ds[i] = -(r[i] + s[i]*dl[i])/l[i];
	}
}


/*
 * Vector addition x = x + y for vectors of length 3339.
 */
void beliefPenaltyMPC_LA_VADD_3339(beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *y)
{
	int i;
	for( i=0; i<3339; i++){
		x[i] += y[i];
	}
}


/*
 * Vector addition x = x + y for vectors of length 1155.
 */
void beliefPenaltyMPC_LA_VADD_1155(beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *y)
{
	int i;
	for( i=0; i<1155; i++){
		x[i] += y[i];
	}
}


/*
 * Vector addition x = x + y for vectors of length 4522.
 */
void beliefPenaltyMPC_LA_VADD_4522(beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *y)
{
	int i;
	for( i=0; i<4522; i++){
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
        for( i=0; i<4522; i++ ){
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
        if( i == 4522 ){
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
    for( i=0; i<3339; i++ ){
        z[i] += a_gamma*dz[i];
    }
    
    /* equality constraint multipliers */
    for( i=0; i<1155; i++ ){
        v[i] += a_gamma*dv[i];
    }
    
    /* inequality constraint multipliers & slacks, also update mu */
    *mu = 0;
    for( i=0; i<4522; i++ ){
        dltemp = l[i] + a_gamma*dl[i]; l[i] = dltemp;
        dstemp = s[i] + a_gamma*ds[i]; s[i] = dstemp;
        *mu += dltemp*dstemp;
    }
    
    *a = a_gamma;
    *mu /= (beliefPenaltyMPC_FLOAT)4522;
    return lsIt;
}




/* VARIABLE DEFINITIONS ------------------------------------------------ */
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_z[3339];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_v[1155];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_dz_aff[3339];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_dv_aff[1155];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_grad_cost[3339];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_grad_eq[3339];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_rd[3339];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_l[4522];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_s[4522];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_lbys[4522];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_dl_aff[4522];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ds_aff[4522];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_dz_cc[3339];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_dv_cc[1155];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_dl_cc[4522];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ds_cc[4522];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ccrhs[4522];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_grad_ineq[3339];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_z00 = beliefPenaltyMPC_z + 0;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzaff00 = beliefPenaltyMPC_dz_aff + 0;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzcc00 = beliefPenaltyMPC_dz_cc + 0;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_rd00 = beliefPenaltyMPC_rd + 0;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lbyrd00[233];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_cost00 = beliefPenaltyMPC_grad_cost + 0;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_eq00 = beliefPenaltyMPC_grad_eq + 0;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_ineq00 = beliefPenaltyMPC_grad_ineq + 0;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ctv00[233];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_v00 = beliefPenaltyMPC_v + 0;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_re00[77];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_beta00[77];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_betacc00[77];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvaff00 = beliefPenaltyMPC_dv_aff + 0;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvcc00 = beliefPenaltyMPC_dv_cc + 0;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_V00[17941];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Yd00[3003];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ld00[3003];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_yy00[77];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_bmy00[77];
int beliefPenaltyMPC_lbIdx00[233] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llb00 = beliefPenaltyMPC_l + 0;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_slb00 = beliefPenaltyMPC_s + 0;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llbbyslb00 = beliefPenaltyMPC_lbys + 0;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_rilb00[233];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbaff00 = beliefPenaltyMPC_dl_aff + 0;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbaff00 = beliefPenaltyMPC_ds_aff + 0;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbcc00 = beliefPenaltyMPC_dl_cc + 0;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbcc00 = beliefPenaltyMPC_ds_cc + 0;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsl00 = beliefPenaltyMPC_ccrhs + 0;
int beliefPenaltyMPC_ubIdx00[79] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lub00 = beliefPenaltyMPC_l + 233;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_sub00 = beliefPenaltyMPC_s + 233;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lubbysub00 = beliefPenaltyMPC_lbys + 233;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_riub00[79];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubaff00 = beliefPenaltyMPC_dl_aff + 233;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubaff00 = beliefPenaltyMPC_ds_aff + 233;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubcc00 = beliefPenaltyMPC_dl_cc + 233;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubcc00 = beliefPenaltyMPC_ds_cc + 233;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsub00 = beliefPenaltyMPC_ccrhs + 233;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Phi00[233];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_D00[233] = {1.0000000000000000E+000, 
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
1.0000000000000000E+000, 
1.0000000000000000E+000, 
1.0000000000000000E+000, 
1.0000000000000000E+000};
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_W00[233];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_z01 = beliefPenaltyMPC_z + 233;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzaff01 = beliefPenaltyMPC_dz_aff + 233;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzcc01 = beliefPenaltyMPC_dz_cc + 233;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_rd01 = beliefPenaltyMPC_rd + 233;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lbyrd01[233];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_cost01 = beliefPenaltyMPC_grad_cost + 233;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_eq01 = beliefPenaltyMPC_grad_eq + 233;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_ineq01 = beliefPenaltyMPC_grad_ineq + 233;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ctv01[233];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_v01 = beliefPenaltyMPC_v + 77;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_re01[77];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_beta01[77];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_betacc01[77];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvaff01 = beliefPenaltyMPC_dv_aff + 77;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvcc01 = beliefPenaltyMPC_dv_cc + 77;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_V01[17941];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Yd01[3003];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ld01[3003];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_yy01[77];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_bmy01[77];
int beliefPenaltyMPC_lbIdx01[233] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llb01 = beliefPenaltyMPC_l + 312;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_slb01 = beliefPenaltyMPC_s + 312;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llbbyslb01 = beliefPenaltyMPC_lbys + 312;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_rilb01[233];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbaff01 = beliefPenaltyMPC_dl_aff + 312;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbaff01 = beliefPenaltyMPC_ds_aff + 312;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbcc01 = beliefPenaltyMPC_dl_cc + 312;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbcc01 = beliefPenaltyMPC_ds_cc + 312;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsl01 = beliefPenaltyMPC_ccrhs + 312;
int beliefPenaltyMPC_ubIdx01[79] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lub01 = beliefPenaltyMPC_l + 545;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_sub01 = beliefPenaltyMPC_s + 545;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lubbysub01 = beliefPenaltyMPC_lbys + 545;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_riub01[79];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubaff01 = beliefPenaltyMPC_dl_aff + 545;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubaff01 = beliefPenaltyMPC_ds_aff + 545;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubcc01 = beliefPenaltyMPC_dl_cc + 545;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubcc01 = beliefPenaltyMPC_ds_cc + 545;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsub01 = beliefPenaltyMPC_ccrhs + 545;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Phi01[233];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_D01[233] = {-1.0000000000000000E+000, 
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
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000};
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_W01[233];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ysd01[5929];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lsd01[5929];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_z02 = beliefPenaltyMPC_z + 466;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzaff02 = beliefPenaltyMPC_dz_aff + 466;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzcc02 = beliefPenaltyMPC_dz_cc + 466;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_rd02 = beliefPenaltyMPC_rd + 466;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lbyrd02[233];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_cost02 = beliefPenaltyMPC_grad_cost + 466;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_eq02 = beliefPenaltyMPC_grad_eq + 466;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_ineq02 = beliefPenaltyMPC_grad_ineq + 466;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ctv02[233];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_v02 = beliefPenaltyMPC_v + 154;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_re02[77];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_beta02[77];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_betacc02[77];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvaff02 = beliefPenaltyMPC_dv_aff + 154;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvcc02 = beliefPenaltyMPC_dv_cc + 154;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_V02[17941];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Yd02[3003];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ld02[3003];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_yy02[77];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_bmy02[77];
int beliefPenaltyMPC_lbIdx02[233] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llb02 = beliefPenaltyMPC_l + 624;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_slb02 = beliefPenaltyMPC_s + 624;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llbbyslb02 = beliefPenaltyMPC_lbys + 624;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_rilb02[233];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbaff02 = beliefPenaltyMPC_dl_aff + 624;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbaff02 = beliefPenaltyMPC_ds_aff + 624;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbcc02 = beliefPenaltyMPC_dl_cc + 624;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbcc02 = beliefPenaltyMPC_ds_cc + 624;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsl02 = beliefPenaltyMPC_ccrhs + 624;
int beliefPenaltyMPC_ubIdx02[79] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lub02 = beliefPenaltyMPC_l + 857;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_sub02 = beliefPenaltyMPC_s + 857;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lubbysub02 = beliefPenaltyMPC_lbys + 857;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_riub02[79];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubaff02 = beliefPenaltyMPC_dl_aff + 857;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubaff02 = beliefPenaltyMPC_ds_aff + 857;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubcc02 = beliefPenaltyMPC_dl_cc + 857;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubcc02 = beliefPenaltyMPC_ds_cc + 857;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsub02 = beliefPenaltyMPC_ccrhs + 857;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Phi02[233];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_W02[233];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ysd02[5929];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lsd02[5929];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_z03 = beliefPenaltyMPC_z + 699;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzaff03 = beliefPenaltyMPC_dz_aff + 699;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzcc03 = beliefPenaltyMPC_dz_cc + 699;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_rd03 = beliefPenaltyMPC_rd + 699;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lbyrd03[233];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_cost03 = beliefPenaltyMPC_grad_cost + 699;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_eq03 = beliefPenaltyMPC_grad_eq + 699;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_ineq03 = beliefPenaltyMPC_grad_ineq + 699;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ctv03[233];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_v03 = beliefPenaltyMPC_v + 231;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_re03[77];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_beta03[77];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_betacc03[77];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvaff03 = beliefPenaltyMPC_dv_aff + 231;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvcc03 = beliefPenaltyMPC_dv_cc + 231;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_V03[17941];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Yd03[3003];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ld03[3003];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_yy03[77];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_bmy03[77];
int beliefPenaltyMPC_lbIdx03[233] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llb03 = beliefPenaltyMPC_l + 936;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_slb03 = beliefPenaltyMPC_s + 936;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llbbyslb03 = beliefPenaltyMPC_lbys + 936;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_rilb03[233];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbaff03 = beliefPenaltyMPC_dl_aff + 936;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbaff03 = beliefPenaltyMPC_ds_aff + 936;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbcc03 = beliefPenaltyMPC_dl_cc + 936;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbcc03 = beliefPenaltyMPC_ds_cc + 936;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsl03 = beliefPenaltyMPC_ccrhs + 936;
int beliefPenaltyMPC_ubIdx03[79] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lub03 = beliefPenaltyMPC_l + 1169;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_sub03 = beliefPenaltyMPC_s + 1169;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lubbysub03 = beliefPenaltyMPC_lbys + 1169;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_riub03[79];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubaff03 = beliefPenaltyMPC_dl_aff + 1169;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubaff03 = beliefPenaltyMPC_ds_aff + 1169;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubcc03 = beliefPenaltyMPC_dl_cc + 1169;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubcc03 = beliefPenaltyMPC_ds_cc + 1169;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsub03 = beliefPenaltyMPC_ccrhs + 1169;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Phi03[233];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_W03[233];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ysd03[5929];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lsd03[5929];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_z04 = beliefPenaltyMPC_z + 932;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzaff04 = beliefPenaltyMPC_dz_aff + 932;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzcc04 = beliefPenaltyMPC_dz_cc + 932;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_rd04 = beliefPenaltyMPC_rd + 932;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lbyrd04[233];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_cost04 = beliefPenaltyMPC_grad_cost + 932;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_eq04 = beliefPenaltyMPC_grad_eq + 932;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_ineq04 = beliefPenaltyMPC_grad_ineq + 932;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ctv04[233];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_v04 = beliefPenaltyMPC_v + 308;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_re04[77];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_beta04[77];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_betacc04[77];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvaff04 = beliefPenaltyMPC_dv_aff + 308;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvcc04 = beliefPenaltyMPC_dv_cc + 308;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_V04[17941];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Yd04[3003];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ld04[3003];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_yy04[77];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_bmy04[77];
int beliefPenaltyMPC_lbIdx04[233] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llb04 = beliefPenaltyMPC_l + 1248;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_slb04 = beliefPenaltyMPC_s + 1248;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llbbyslb04 = beliefPenaltyMPC_lbys + 1248;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_rilb04[233];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbaff04 = beliefPenaltyMPC_dl_aff + 1248;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbaff04 = beliefPenaltyMPC_ds_aff + 1248;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbcc04 = beliefPenaltyMPC_dl_cc + 1248;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbcc04 = beliefPenaltyMPC_ds_cc + 1248;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsl04 = beliefPenaltyMPC_ccrhs + 1248;
int beliefPenaltyMPC_ubIdx04[79] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lub04 = beliefPenaltyMPC_l + 1481;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_sub04 = beliefPenaltyMPC_s + 1481;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lubbysub04 = beliefPenaltyMPC_lbys + 1481;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_riub04[79];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubaff04 = beliefPenaltyMPC_dl_aff + 1481;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubaff04 = beliefPenaltyMPC_ds_aff + 1481;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubcc04 = beliefPenaltyMPC_dl_cc + 1481;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubcc04 = beliefPenaltyMPC_ds_cc + 1481;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsub04 = beliefPenaltyMPC_ccrhs + 1481;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Phi04[233];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_W04[233];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ysd04[5929];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lsd04[5929];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_z05 = beliefPenaltyMPC_z + 1165;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzaff05 = beliefPenaltyMPC_dz_aff + 1165;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzcc05 = beliefPenaltyMPC_dz_cc + 1165;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_rd05 = beliefPenaltyMPC_rd + 1165;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lbyrd05[233];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_cost05 = beliefPenaltyMPC_grad_cost + 1165;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_eq05 = beliefPenaltyMPC_grad_eq + 1165;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_ineq05 = beliefPenaltyMPC_grad_ineq + 1165;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ctv05[233];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_v05 = beliefPenaltyMPC_v + 385;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_re05[77];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_beta05[77];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_betacc05[77];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvaff05 = beliefPenaltyMPC_dv_aff + 385;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvcc05 = beliefPenaltyMPC_dv_cc + 385;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_V05[17941];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Yd05[3003];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ld05[3003];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_yy05[77];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_bmy05[77];
int beliefPenaltyMPC_lbIdx05[233] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llb05 = beliefPenaltyMPC_l + 1560;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_slb05 = beliefPenaltyMPC_s + 1560;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llbbyslb05 = beliefPenaltyMPC_lbys + 1560;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_rilb05[233];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbaff05 = beliefPenaltyMPC_dl_aff + 1560;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbaff05 = beliefPenaltyMPC_ds_aff + 1560;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbcc05 = beliefPenaltyMPC_dl_cc + 1560;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbcc05 = beliefPenaltyMPC_ds_cc + 1560;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsl05 = beliefPenaltyMPC_ccrhs + 1560;
int beliefPenaltyMPC_ubIdx05[79] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lub05 = beliefPenaltyMPC_l + 1793;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_sub05 = beliefPenaltyMPC_s + 1793;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lubbysub05 = beliefPenaltyMPC_lbys + 1793;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_riub05[79];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubaff05 = beliefPenaltyMPC_dl_aff + 1793;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubaff05 = beliefPenaltyMPC_ds_aff + 1793;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubcc05 = beliefPenaltyMPC_dl_cc + 1793;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubcc05 = beliefPenaltyMPC_ds_cc + 1793;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsub05 = beliefPenaltyMPC_ccrhs + 1793;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Phi05[233];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_W05[233];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ysd05[5929];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lsd05[5929];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_z06 = beliefPenaltyMPC_z + 1398;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzaff06 = beliefPenaltyMPC_dz_aff + 1398;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzcc06 = beliefPenaltyMPC_dz_cc + 1398;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_rd06 = beliefPenaltyMPC_rd + 1398;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lbyrd06[233];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_cost06 = beliefPenaltyMPC_grad_cost + 1398;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_eq06 = beliefPenaltyMPC_grad_eq + 1398;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_ineq06 = beliefPenaltyMPC_grad_ineq + 1398;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ctv06[233];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_v06 = beliefPenaltyMPC_v + 462;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_re06[77];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_beta06[77];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_betacc06[77];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvaff06 = beliefPenaltyMPC_dv_aff + 462;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvcc06 = beliefPenaltyMPC_dv_cc + 462;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_V06[17941];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Yd06[3003];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ld06[3003];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_yy06[77];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_bmy06[77];
int beliefPenaltyMPC_lbIdx06[233] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llb06 = beliefPenaltyMPC_l + 1872;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_slb06 = beliefPenaltyMPC_s + 1872;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llbbyslb06 = beliefPenaltyMPC_lbys + 1872;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_rilb06[233];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbaff06 = beliefPenaltyMPC_dl_aff + 1872;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbaff06 = beliefPenaltyMPC_ds_aff + 1872;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbcc06 = beliefPenaltyMPC_dl_cc + 1872;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbcc06 = beliefPenaltyMPC_ds_cc + 1872;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsl06 = beliefPenaltyMPC_ccrhs + 1872;
int beliefPenaltyMPC_ubIdx06[79] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lub06 = beliefPenaltyMPC_l + 2105;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_sub06 = beliefPenaltyMPC_s + 2105;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lubbysub06 = beliefPenaltyMPC_lbys + 2105;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_riub06[79];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubaff06 = beliefPenaltyMPC_dl_aff + 2105;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubaff06 = beliefPenaltyMPC_ds_aff + 2105;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubcc06 = beliefPenaltyMPC_dl_cc + 2105;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubcc06 = beliefPenaltyMPC_ds_cc + 2105;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsub06 = beliefPenaltyMPC_ccrhs + 2105;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Phi06[233];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_W06[233];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ysd06[5929];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lsd06[5929];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_z07 = beliefPenaltyMPC_z + 1631;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzaff07 = beliefPenaltyMPC_dz_aff + 1631;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzcc07 = beliefPenaltyMPC_dz_cc + 1631;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_rd07 = beliefPenaltyMPC_rd + 1631;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lbyrd07[233];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_cost07 = beliefPenaltyMPC_grad_cost + 1631;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_eq07 = beliefPenaltyMPC_grad_eq + 1631;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_ineq07 = beliefPenaltyMPC_grad_ineq + 1631;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ctv07[233];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_v07 = beliefPenaltyMPC_v + 539;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_re07[77];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_beta07[77];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_betacc07[77];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvaff07 = beliefPenaltyMPC_dv_aff + 539;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvcc07 = beliefPenaltyMPC_dv_cc + 539;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_V07[17941];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Yd07[3003];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ld07[3003];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_yy07[77];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_bmy07[77];
int beliefPenaltyMPC_lbIdx07[233] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llb07 = beliefPenaltyMPC_l + 2184;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_slb07 = beliefPenaltyMPC_s + 2184;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llbbyslb07 = beliefPenaltyMPC_lbys + 2184;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_rilb07[233];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbaff07 = beliefPenaltyMPC_dl_aff + 2184;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbaff07 = beliefPenaltyMPC_ds_aff + 2184;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbcc07 = beliefPenaltyMPC_dl_cc + 2184;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbcc07 = beliefPenaltyMPC_ds_cc + 2184;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsl07 = beliefPenaltyMPC_ccrhs + 2184;
int beliefPenaltyMPC_ubIdx07[79] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lub07 = beliefPenaltyMPC_l + 2417;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_sub07 = beliefPenaltyMPC_s + 2417;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lubbysub07 = beliefPenaltyMPC_lbys + 2417;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_riub07[79];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubaff07 = beliefPenaltyMPC_dl_aff + 2417;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubaff07 = beliefPenaltyMPC_ds_aff + 2417;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubcc07 = beliefPenaltyMPC_dl_cc + 2417;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubcc07 = beliefPenaltyMPC_ds_cc + 2417;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsub07 = beliefPenaltyMPC_ccrhs + 2417;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Phi07[233];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_W07[233];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ysd07[5929];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lsd07[5929];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_z08 = beliefPenaltyMPC_z + 1864;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzaff08 = beliefPenaltyMPC_dz_aff + 1864;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzcc08 = beliefPenaltyMPC_dz_cc + 1864;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_rd08 = beliefPenaltyMPC_rd + 1864;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lbyrd08[233];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_cost08 = beliefPenaltyMPC_grad_cost + 1864;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_eq08 = beliefPenaltyMPC_grad_eq + 1864;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_ineq08 = beliefPenaltyMPC_grad_ineq + 1864;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ctv08[233];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_v08 = beliefPenaltyMPC_v + 616;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_re08[77];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_beta08[77];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_betacc08[77];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvaff08 = beliefPenaltyMPC_dv_aff + 616;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvcc08 = beliefPenaltyMPC_dv_cc + 616;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_V08[17941];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Yd08[3003];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ld08[3003];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_yy08[77];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_bmy08[77];
int beliefPenaltyMPC_lbIdx08[233] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llb08 = beliefPenaltyMPC_l + 2496;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_slb08 = beliefPenaltyMPC_s + 2496;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llbbyslb08 = beliefPenaltyMPC_lbys + 2496;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_rilb08[233];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbaff08 = beliefPenaltyMPC_dl_aff + 2496;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbaff08 = beliefPenaltyMPC_ds_aff + 2496;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbcc08 = beliefPenaltyMPC_dl_cc + 2496;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbcc08 = beliefPenaltyMPC_ds_cc + 2496;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsl08 = beliefPenaltyMPC_ccrhs + 2496;
int beliefPenaltyMPC_ubIdx08[79] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lub08 = beliefPenaltyMPC_l + 2729;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_sub08 = beliefPenaltyMPC_s + 2729;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lubbysub08 = beliefPenaltyMPC_lbys + 2729;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_riub08[79];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubaff08 = beliefPenaltyMPC_dl_aff + 2729;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubaff08 = beliefPenaltyMPC_ds_aff + 2729;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubcc08 = beliefPenaltyMPC_dl_cc + 2729;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubcc08 = beliefPenaltyMPC_ds_cc + 2729;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsub08 = beliefPenaltyMPC_ccrhs + 2729;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Phi08[233];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_W08[233];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ysd08[5929];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lsd08[5929];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_z09 = beliefPenaltyMPC_z + 2097;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzaff09 = beliefPenaltyMPC_dz_aff + 2097;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzcc09 = beliefPenaltyMPC_dz_cc + 2097;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_rd09 = beliefPenaltyMPC_rd + 2097;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lbyrd09[233];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_cost09 = beliefPenaltyMPC_grad_cost + 2097;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_eq09 = beliefPenaltyMPC_grad_eq + 2097;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_ineq09 = beliefPenaltyMPC_grad_ineq + 2097;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ctv09[233];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_v09 = beliefPenaltyMPC_v + 693;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_re09[77];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_beta09[77];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_betacc09[77];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvaff09 = beliefPenaltyMPC_dv_aff + 693;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvcc09 = beliefPenaltyMPC_dv_cc + 693;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_V09[17941];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Yd09[3003];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ld09[3003];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_yy09[77];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_bmy09[77];
int beliefPenaltyMPC_lbIdx09[233] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llb09 = beliefPenaltyMPC_l + 2808;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_slb09 = beliefPenaltyMPC_s + 2808;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llbbyslb09 = beliefPenaltyMPC_lbys + 2808;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_rilb09[233];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbaff09 = beliefPenaltyMPC_dl_aff + 2808;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbaff09 = beliefPenaltyMPC_ds_aff + 2808;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbcc09 = beliefPenaltyMPC_dl_cc + 2808;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbcc09 = beliefPenaltyMPC_ds_cc + 2808;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsl09 = beliefPenaltyMPC_ccrhs + 2808;
int beliefPenaltyMPC_ubIdx09[79] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lub09 = beliefPenaltyMPC_l + 3041;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_sub09 = beliefPenaltyMPC_s + 3041;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lubbysub09 = beliefPenaltyMPC_lbys + 3041;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_riub09[79];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubaff09 = beliefPenaltyMPC_dl_aff + 3041;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubaff09 = beliefPenaltyMPC_ds_aff + 3041;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubcc09 = beliefPenaltyMPC_dl_cc + 3041;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubcc09 = beliefPenaltyMPC_ds_cc + 3041;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsub09 = beliefPenaltyMPC_ccrhs + 3041;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Phi09[233];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_W09[233];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ysd09[5929];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lsd09[5929];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_z10 = beliefPenaltyMPC_z + 2330;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzaff10 = beliefPenaltyMPC_dz_aff + 2330;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzcc10 = beliefPenaltyMPC_dz_cc + 2330;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_rd10 = beliefPenaltyMPC_rd + 2330;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lbyrd10[233];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_cost10 = beliefPenaltyMPC_grad_cost + 2330;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_eq10 = beliefPenaltyMPC_grad_eq + 2330;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_ineq10 = beliefPenaltyMPC_grad_ineq + 2330;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ctv10[233];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_v10 = beliefPenaltyMPC_v + 770;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_re10[77];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_beta10[77];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_betacc10[77];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvaff10 = beliefPenaltyMPC_dv_aff + 770;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvcc10 = beliefPenaltyMPC_dv_cc + 770;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_V10[17941];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Yd10[3003];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ld10[3003];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_yy10[77];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_bmy10[77];
int beliefPenaltyMPC_lbIdx10[233] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llb10 = beliefPenaltyMPC_l + 3120;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_slb10 = beliefPenaltyMPC_s + 3120;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llbbyslb10 = beliefPenaltyMPC_lbys + 3120;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_rilb10[233];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbaff10 = beliefPenaltyMPC_dl_aff + 3120;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbaff10 = beliefPenaltyMPC_ds_aff + 3120;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbcc10 = beliefPenaltyMPC_dl_cc + 3120;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbcc10 = beliefPenaltyMPC_ds_cc + 3120;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsl10 = beliefPenaltyMPC_ccrhs + 3120;
int beliefPenaltyMPC_ubIdx10[79] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lub10 = beliefPenaltyMPC_l + 3353;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_sub10 = beliefPenaltyMPC_s + 3353;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lubbysub10 = beliefPenaltyMPC_lbys + 3353;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_riub10[79];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubaff10 = beliefPenaltyMPC_dl_aff + 3353;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubaff10 = beliefPenaltyMPC_ds_aff + 3353;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubcc10 = beliefPenaltyMPC_dl_cc + 3353;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubcc10 = beliefPenaltyMPC_ds_cc + 3353;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsub10 = beliefPenaltyMPC_ccrhs + 3353;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Phi10[233];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_W10[233];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ysd10[5929];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lsd10[5929];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_z11 = beliefPenaltyMPC_z + 2563;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzaff11 = beliefPenaltyMPC_dz_aff + 2563;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzcc11 = beliefPenaltyMPC_dz_cc + 2563;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_rd11 = beliefPenaltyMPC_rd + 2563;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lbyrd11[233];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_cost11 = beliefPenaltyMPC_grad_cost + 2563;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_eq11 = beliefPenaltyMPC_grad_eq + 2563;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_ineq11 = beliefPenaltyMPC_grad_ineq + 2563;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ctv11[233];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_v11 = beliefPenaltyMPC_v + 847;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_re11[77];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_beta11[77];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_betacc11[77];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvaff11 = beliefPenaltyMPC_dv_aff + 847;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvcc11 = beliefPenaltyMPC_dv_cc + 847;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_V11[17941];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Yd11[3003];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ld11[3003];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_yy11[77];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_bmy11[77];
int beliefPenaltyMPC_lbIdx11[233] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llb11 = beliefPenaltyMPC_l + 3432;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_slb11 = beliefPenaltyMPC_s + 3432;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llbbyslb11 = beliefPenaltyMPC_lbys + 3432;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_rilb11[233];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbaff11 = beliefPenaltyMPC_dl_aff + 3432;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbaff11 = beliefPenaltyMPC_ds_aff + 3432;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbcc11 = beliefPenaltyMPC_dl_cc + 3432;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbcc11 = beliefPenaltyMPC_ds_cc + 3432;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsl11 = beliefPenaltyMPC_ccrhs + 3432;
int beliefPenaltyMPC_ubIdx11[79] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lub11 = beliefPenaltyMPC_l + 3665;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_sub11 = beliefPenaltyMPC_s + 3665;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lubbysub11 = beliefPenaltyMPC_lbys + 3665;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_riub11[79];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubaff11 = beliefPenaltyMPC_dl_aff + 3665;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubaff11 = beliefPenaltyMPC_ds_aff + 3665;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubcc11 = beliefPenaltyMPC_dl_cc + 3665;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubcc11 = beliefPenaltyMPC_ds_cc + 3665;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsub11 = beliefPenaltyMPC_ccrhs + 3665;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Phi11[233];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_W11[233];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ysd11[5929];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lsd11[5929];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_z12 = beliefPenaltyMPC_z + 2796;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzaff12 = beliefPenaltyMPC_dz_aff + 2796;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzcc12 = beliefPenaltyMPC_dz_cc + 2796;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_rd12 = beliefPenaltyMPC_rd + 2796;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lbyrd12[233];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_cost12 = beliefPenaltyMPC_grad_cost + 2796;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_eq12 = beliefPenaltyMPC_grad_eq + 2796;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_ineq12 = beliefPenaltyMPC_grad_ineq + 2796;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ctv12[233];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_v12 = beliefPenaltyMPC_v + 924;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_re12[77];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_beta12[77];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_betacc12[77];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvaff12 = beliefPenaltyMPC_dv_aff + 924;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvcc12 = beliefPenaltyMPC_dv_cc + 924;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_V12[17941];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Yd12[3003];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ld12[3003];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_yy12[77];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_bmy12[77];
int beliefPenaltyMPC_lbIdx12[233] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llb12 = beliefPenaltyMPC_l + 3744;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_slb12 = beliefPenaltyMPC_s + 3744;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llbbyslb12 = beliefPenaltyMPC_lbys + 3744;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_rilb12[233];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbaff12 = beliefPenaltyMPC_dl_aff + 3744;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbaff12 = beliefPenaltyMPC_ds_aff + 3744;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbcc12 = beliefPenaltyMPC_dl_cc + 3744;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbcc12 = beliefPenaltyMPC_ds_cc + 3744;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsl12 = beliefPenaltyMPC_ccrhs + 3744;
int beliefPenaltyMPC_ubIdx12[79] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lub12 = beliefPenaltyMPC_l + 3977;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_sub12 = beliefPenaltyMPC_s + 3977;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lubbysub12 = beliefPenaltyMPC_lbys + 3977;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_riub12[79];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubaff12 = beliefPenaltyMPC_dl_aff + 3977;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubaff12 = beliefPenaltyMPC_ds_aff + 3977;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubcc12 = beliefPenaltyMPC_dl_cc + 3977;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubcc12 = beliefPenaltyMPC_ds_cc + 3977;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsub12 = beliefPenaltyMPC_ccrhs + 3977;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Phi12[233];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_W12[233];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ysd12[5929];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lsd12[5929];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_z13 = beliefPenaltyMPC_z + 3029;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzaff13 = beliefPenaltyMPC_dz_aff + 3029;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzcc13 = beliefPenaltyMPC_dz_cc + 3029;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_rd13 = beliefPenaltyMPC_rd + 3029;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lbyrd13[233];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_cost13 = beliefPenaltyMPC_grad_cost + 3029;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_eq13 = beliefPenaltyMPC_grad_eq + 3029;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_ineq13 = beliefPenaltyMPC_grad_ineq + 3029;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ctv13[233];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_v13 = beliefPenaltyMPC_v + 1001;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_re13[77];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_beta13[77];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_betacc13[77];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvaff13 = beliefPenaltyMPC_dv_aff + 1001;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvcc13 = beliefPenaltyMPC_dv_cc + 1001;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_V13[17941];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Yd13[3003];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ld13[3003];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_yy13[77];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_bmy13[77];
int beliefPenaltyMPC_lbIdx13[233] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llb13 = beliefPenaltyMPC_l + 4056;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_slb13 = beliefPenaltyMPC_s + 4056;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llbbyslb13 = beliefPenaltyMPC_lbys + 4056;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_rilb13[233];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbaff13 = beliefPenaltyMPC_dl_aff + 4056;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbaff13 = beliefPenaltyMPC_ds_aff + 4056;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbcc13 = beliefPenaltyMPC_dl_cc + 4056;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbcc13 = beliefPenaltyMPC_ds_cc + 4056;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsl13 = beliefPenaltyMPC_ccrhs + 4056;
int beliefPenaltyMPC_ubIdx13[79] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lub13 = beliefPenaltyMPC_l + 4289;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_sub13 = beliefPenaltyMPC_s + 4289;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lubbysub13 = beliefPenaltyMPC_lbys + 4289;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_riub13[79];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubaff13 = beliefPenaltyMPC_dl_aff + 4289;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubaff13 = beliefPenaltyMPC_ds_aff + 4289;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubcc13 = beliefPenaltyMPC_dl_cc + 4289;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubcc13 = beliefPenaltyMPC_ds_cc + 4289;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsub13 = beliefPenaltyMPC_ccrhs + 4289;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Phi13[233];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_W13[233];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ysd13[5929];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lsd13[5929];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_f14[77] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_z14 = beliefPenaltyMPC_z + 3262;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzaff14 = beliefPenaltyMPC_dz_aff + 3262;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzcc14 = beliefPenaltyMPC_dz_cc + 3262;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_rd14 = beliefPenaltyMPC_rd + 3262;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lbyrd14[77];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_cost14 = beliefPenaltyMPC_grad_cost + 3262;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_eq14 = beliefPenaltyMPC_grad_eq + 3262;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_ineq14 = beliefPenaltyMPC_grad_ineq + 3262;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ctv14[77];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_v14 = beliefPenaltyMPC_v + 1078;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_re14[77];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_beta14[77];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_betacc14[77];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvaff14 = beliefPenaltyMPC_dv_aff + 1078;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvcc14 = beliefPenaltyMPC_dv_cc + 1078;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_V14[5929];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Yd14[3003];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ld14[3003];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_yy14[77];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_bmy14[77];
int beliefPenaltyMPC_lbIdx14[77] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llb14 = beliefPenaltyMPC_l + 4368;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_slb14 = beliefPenaltyMPC_s + 4368;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llbbyslb14 = beliefPenaltyMPC_lbys + 4368;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_rilb14[77];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbaff14 = beliefPenaltyMPC_dl_aff + 4368;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbaff14 = beliefPenaltyMPC_ds_aff + 4368;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbcc14 = beliefPenaltyMPC_dl_cc + 4368;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbcc14 = beliefPenaltyMPC_ds_cc + 4368;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsl14 = beliefPenaltyMPC_ccrhs + 4368;
int beliefPenaltyMPC_ubIdx14[77] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lub14 = beliefPenaltyMPC_l + 4445;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_sub14 = beliefPenaltyMPC_s + 4445;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lubbysub14 = beliefPenaltyMPC_lbys + 4445;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_riub14[77];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubaff14 = beliefPenaltyMPC_dl_aff + 4445;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubaff14 = beliefPenaltyMPC_ds_aff + 4445;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubcc14 = beliefPenaltyMPC_dl_cc + 4445;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubcc14 = beliefPenaltyMPC_ds_cc + 4445;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsub14 = beliefPenaltyMPC_ccrhs + 4445;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Phi14[77];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_D14[77] = {-1.0000000000000000E+000, 
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
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000};
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_W14[77];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ysd14[5929];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lsd14[5929];
beliefPenaltyMPC_FLOAT musigma;
beliefPenaltyMPC_FLOAT sigma_3rdroot;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Diag1_0[233];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Diag2_0[233];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_L_0[27028];




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
beliefPenaltyMPC_LA_INITIALIZEVECTOR_3339(beliefPenaltyMPC_z, 0);
beliefPenaltyMPC_LA_INITIALIZEVECTOR_1155(beliefPenaltyMPC_v, 1);
beliefPenaltyMPC_LA_INITIALIZEVECTOR_4522(beliefPenaltyMPC_l, 1);
beliefPenaltyMPC_LA_INITIALIZEVECTOR_4522(beliefPenaltyMPC_s, 1);
info->mu = 0;
beliefPenaltyMPC_LA_DOTACC_4522(beliefPenaltyMPC_l, beliefPenaltyMPC_s, &info->mu);
info->mu /= 4522;
while( 1 ){
info->pobj = 0;
beliefPenaltyMPC_LA_DIAG_QUADFCN_233(params->H1, params->f1, beliefPenaltyMPC_z00, beliefPenaltyMPC_grad_cost00, &info->pobj);
beliefPenaltyMPC_LA_DIAG_QUADFCN_233(params->H2, params->f2, beliefPenaltyMPC_z01, beliefPenaltyMPC_grad_cost01, &info->pobj);
beliefPenaltyMPC_LA_DIAG_QUADFCN_233(params->H3, params->f3, beliefPenaltyMPC_z02, beliefPenaltyMPC_grad_cost02, &info->pobj);
beliefPenaltyMPC_LA_DIAG_QUADFCN_233(params->H4, params->f4, beliefPenaltyMPC_z03, beliefPenaltyMPC_grad_cost03, &info->pobj);
beliefPenaltyMPC_LA_DIAG_QUADFCN_233(params->H5, params->f5, beliefPenaltyMPC_z04, beliefPenaltyMPC_grad_cost04, &info->pobj);
beliefPenaltyMPC_LA_DIAG_QUADFCN_233(params->H6, params->f6, beliefPenaltyMPC_z05, beliefPenaltyMPC_grad_cost05, &info->pobj);
beliefPenaltyMPC_LA_DIAG_QUADFCN_233(params->H7, params->f7, beliefPenaltyMPC_z06, beliefPenaltyMPC_grad_cost06, &info->pobj);
beliefPenaltyMPC_LA_DIAG_QUADFCN_233(params->H8, params->f8, beliefPenaltyMPC_z07, beliefPenaltyMPC_grad_cost07, &info->pobj);
beliefPenaltyMPC_LA_DIAG_QUADFCN_233(params->H9, params->f9, beliefPenaltyMPC_z08, beliefPenaltyMPC_grad_cost08, &info->pobj);
beliefPenaltyMPC_LA_DIAG_QUADFCN_233(params->H10, params->f10, beliefPenaltyMPC_z09, beliefPenaltyMPC_grad_cost09, &info->pobj);
beliefPenaltyMPC_LA_DIAG_QUADFCN_233(params->H11, params->f11, beliefPenaltyMPC_z10, beliefPenaltyMPC_grad_cost10, &info->pobj);
beliefPenaltyMPC_LA_DIAG_QUADFCN_233(params->H12, params->f12, beliefPenaltyMPC_z11, beliefPenaltyMPC_grad_cost11, &info->pobj);
beliefPenaltyMPC_LA_DIAG_QUADFCN_233(params->H13, params->f13, beliefPenaltyMPC_z12, beliefPenaltyMPC_grad_cost12, &info->pobj);
beliefPenaltyMPC_LA_DIAG_QUADFCN_233(params->H14, params->f14, beliefPenaltyMPC_z13, beliefPenaltyMPC_grad_cost13, &info->pobj);
beliefPenaltyMPC_LA_DIAG_QUADFCN_77(params->H15, beliefPenaltyMPC_f14, beliefPenaltyMPC_z14, beliefPenaltyMPC_grad_cost14, &info->pobj);
info->res_eq = 0;
info->dgap = 0;
beliefPenaltyMPC_LA_DIAGZERO_MVMSUB6_77(beliefPenaltyMPC_D00, beliefPenaltyMPC_z00, params->e1, beliefPenaltyMPC_v00, beliefPenaltyMPC_re00, &info->dgap, &info->res_eq);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_77_233_233(params->C1, beliefPenaltyMPC_z00, beliefPenaltyMPC_D01, beliefPenaltyMPC_z01, params->e2, beliefPenaltyMPC_v01, beliefPenaltyMPC_re01, &info->dgap, &info->res_eq);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_77_233_233(params->C2, beliefPenaltyMPC_z01, beliefPenaltyMPC_D01, beliefPenaltyMPC_z02, params->e3, beliefPenaltyMPC_v02, beliefPenaltyMPC_re02, &info->dgap, &info->res_eq);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_77_233_233(params->C3, beliefPenaltyMPC_z02, beliefPenaltyMPC_D01, beliefPenaltyMPC_z03, params->e4, beliefPenaltyMPC_v03, beliefPenaltyMPC_re03, &info->dgap, &info->res_eq);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_77_233_233(params->C4, beliefPenaltyMPC_z03, beliefPenaltyMPC_D01, beliefPenaltyMPC_z04, params->e5, beliefPenaltyMPC_v04, beliefPenaltyMPC_re04, &info->dgap, &info->res_eq);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_77_233_233(params->C5, beliefPenaltyMPC_z04, beliefPenaltyMPC_D01, beliefPenaltyMPC_z05, params->e6, beliefPenaltyMPC_v05, beliefPenaltyMPC_re05, &info->dgap, &info->res_eq);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_77_233_233(params->C6, beliefPenaltyMPC_z05, beliefPenaltyMPC_D01, beliefPenaltyMPC_z06, params->e7, beliefPenaltyMPC_v06, beliefPenaltyMPC_re06, &info->dgap, &info->res_eq);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_77_233_233(params->C7, beliefPenaltyMPC_z06, beliefPenaltyMPC_D01, beliefPenaltyMPC_z07, params->e8, beliefPenaltyMPC_v07, beliefPenaltyMPC_re07, &info->dgap, &info->res_eq);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_77_233_233(params->C8, beliefPenaltyMPC_z07, beliefPenaltyMPC_D01, beliefPenaltyMPC_z08, params->e9, beliefPenaltyMPC_v08, beliefPenaltyMPC_re08, &info->dgap, &info->res_eq);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_77_233_233(params->C9, beliefPenaltyMPC_z08, beliefPenaltyMPC_D01, beliefPenaltyMPC_z09, params->e10, beliefPenaltyMPC_v09, beliefPenaltyMPC_re09, &info->dgap, &info->res_eq);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_77_233_233(params->C10, beliefPenaltyMPC_z09, beliefPenaltyMPC_D01, beliefPenaltyMPC_z10, params->e11, beliefPenaltyMPC_v10, beliefPenaltyMPC_re10, &info->dgap, &info->res_eq);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_77_233_233(params->C11, beliefPenaltyMPC_z10, beliefPenaltyMPC_D01, beliefPenaltyMPC_z11, params->e12, beliefPenaltyMPC_v11, beliefPenaltyMPC_re11, &info->dgap, &info->res_eq);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_77_233_233(params->C12, beliefPenaltyMPC_z11, beliefPenaltyMPC_D01, beliefPenaltyMPC_z12, params->e13, beliefPenaltyMPC_v12, beliefPenaltyMPC_re12, &info->dgap, &info->res_eq);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_77_233_233(params->C13, beliefPenaltyMPC_z12, beliefPenaltyMPC_D01, beliefPenaltyMPC_z13, params->e14, beliefPenaltyMPC_v13, beliefPenaltyMPC_re13, &info->dgap, &info->res_eq);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_77_233_77(params->C14, beliefPenaltyMPC_z13, beliefPenaltyMPC_D14, beliefPenaltyMPC_z14, params->e15, beliefPenaltyMPC_v14, beliefPenaltyMPC_re14, &info->dgap, &info->res_eq);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_77_233_77(params->C1, beliefPenaltyMPC_v01, beliefPenaltyMPC_D00, beliefPenaltyMPC_v00, beliefPenaltyMPC_grad_eq00);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_77_233_77(params->C2, beliefPenaltyMPC_v02, beliefPenaltyMPC_D01, beliefPenaltyMPC_v01, beliefPenaltyMPC_grad_eq01);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_77_233_77(params->C3, beliefPenaltyMPC_v03, beliefPenaltyMPC_D01, beliefPenaltyMPC_v02, beliefPenaltyMPC_grad_eq02);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_77_233_77(params->C4, beliefPenaltyMPC_v04, beliefPenaltyMPC_D01, beliefPenaltyMPC_v03, beliefPenaltyMPC_grad_eq03);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_77_233_77(params->C5, beliefPenaltyMPC_v05, beliefPenaltyMPC_D01, beliefPenaltyMPC_v04, beliefPenaltyMPC_grad_eq04);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_77_233_77(params->C6, beliefPenaltyMPC_v06, beliefPenaltyMPC_D01, beliefPenaltyMPC_v05, beliefPenaltyMPC_grad_eq05);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_77_233_77(params->C7, beliefPenaltyMPC_v07, beliefPenaltyMPC_D01, beliefPenaltyMPC_v06, beliefPenaltyMPC_grad_eq06);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_77_233_77(params->C8, beliefPenaltyMPC_v08, beliefPenaltyMPC_D01, beliefPenaltyMPC_v07, beliefPenaltyMPC_grad_eq07);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_77_233_77(params->C9, beliefPenaltyMPC_v09, beliefPenaltyMPC_D01, beliefPenaltyMPC_v08, beliefPenaltyMPC_grad_eq08);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_77_233_77(params->C10, beliefPenaltyMPC_v10, beliefPenaltyMPC_D01, beliefPenaltyMPC_v09, beliefPenaltyMPC_grad_eq09);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_77_233_77(params->C11, beliefPenaltyMPC_v11, beliefPenaltyMPC_D01, beliefPenaltyMPC_v10, beliefPenaltyMPC_grad_eq10);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_77_233_77(params->C12, beliefPenaltyMPC_v12, beliefPenaltyMPC_D01, beliefPenaltyMPC_v11, beliefPenaltyMPC_grad_eq11);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_77_233_77(params->C13, beliefPenaltyMPC_v13, beliefPenaltyMPC_D01, beliefPenaltyMPC_v12, beliefPenaltyMPC_grad_eq12);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_77_233_77(params->C14, beliefPenaltyMPC_v14, beliefPenaltyMPC_D01, beliefPenaltyMPC_v13, beliefPenaltyMPC_grad_eq13);
beliefPenaltyMPC_LA_DIAGZERO_MTVM_77_77(beliefPenaltyMPC_D14, beliefPenaltyMPC_v14, beliefPenaltyMPC_grad_eq14);
info->res_ineq = 0;
beliefPenaltyMPC_LA_VSUBADD3_233(params->lb1, beliefPenaltyMPC_z00, beliefPenaltyMPC_lbIdx00, beliefPenaltyMPC_llb00, beliefPenaltyMPC_slb00, beliefPenaltyMPC_rilb00, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD2_79(beliefPenaltyMPC_z00, beliefPenaltyMPC_ubIdx00, params->ub1, beliefPenaltyMPC_lub00, beliefPenaltyMPC_sub00, beliefPenaltyMPC_riub00, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD3_233(params->lb2, beliefPenaltyMPC_z01, beliefPenaltyMPC_lbIdx01, beliefPenaltyMPC_llb01, beliefPenaltyMPC_slb01, beliefPenaltyMPC_rilb01, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD2_79(beliefPenaltyMPC_z01, beliefPenaltyMPC_ubIdx01, params->ub2, beliefPenaltyMPC_lub01, beliefPenaltyMPC_sub01, beliefPenaltyMPC_riub01, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD3_233(params->lb3, beliefPenaltyMPC_z02, beliefPenaltyMPC_lbIdx02, beliefPenaltyMPC_llb02, beliefPenaltyMPC_slb02, beliefPenaltyMPC_rilb02, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD2_79(beliefPenaltyMPC_z02, beliefPenaltyMPC_ubIdx02, params->ub3, beliefPenaltyMPC_lub02, beliefPenaltyMPC_sub02, beliefPenaltyMPC_riub02, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD3_233(params->lb4, beliefPenaltyMPC_z03, beliefPenaltyMPC_lbIdx03, beliefPenaltyMPC_llb03, beliefPenaltyMPC_slb03, beliefPenaltyMPC_rilb03, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD2_79(beliefPenaltyMPC_z03, beliefPenaltyMPC_ubIdx03, params->ub4, beliefPenaltyMPC_lub03, beliefPenaltyMPC_sub03, beliefPenaltyMPC_riub03, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD3_233(params->lb5, beliefPenaltyMPC_z04, beliefPenaltyMPC_lbIdx04, beliefPenaltyMPC_llb04, beliefPenaltyMPC_slb04, beliefPenaltyMPC_rilb04, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD2_79(beliefPenaltyMPC_z04, beliefPenaltyMPC_ubIdx04, params->ub5, beliefPenaltyMPC_lub04, beliefPenaltyMPC_sub04, beliefPenaltyMPC_riub04, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD3_233(params->lb6, beliefPenaltyMPC_z05, beliefPenaltyMPC_lbIdx05, beliefPenaltyMPC_llb05, beliefPenaltyMPC_slb05, beliefPenaltyMPC_rilb05, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD2_79(beliefPenaltyMPC_z05, beliefPenaltyMPC_ubIdx05, params->ub6, beliefPenaltyMPC_lub05, beliefPenaltyMPC_sub05, beliefPenaltyMPC_riub05, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD3_233(params->lb7, beliefPenaltyMPC_z06, beliefPenaltyMPC_lbIdx06, beliefPenaltyMPC_llb06, beliefPenaltyMPC_slb06, beliefPenaltyMPC_rilb06, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD2_79(beliefPenaltyMPC_z06, beliefPenaltyMPC_ubIdx06, params->ub7, beliefPenaltyMPC_lub06, beliefPenaltyMPC_sub06, beliefPenaltyMPC_riub06, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD3_233(params->lb8, beliefPenaltyMPC_z07, beliefPenaltyMPC_lbIdx07, beliefPenaltyMPC_llb07, beliefPenaltyMPC_slb07, beliefPenaltyMPC_rilb07, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD2_79(beliefPenaltyMPC_z07, beliefPenaltyMPC_ubIdx07, params->ub8, beliefPenaltyMPC_lub07, beliefPenaltyMPC_sub07, beliefPenaltyMPC_riub07, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD3_233(params->lb9, beliefPenaltyMPC_z08, beliefPenaltyMPC_lbIdx08, beliefPenaltyMPC_llb08, beliefPenaltyMPC_slb08, beliefPenaltyMPC_rilb08, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD2_79(beliefPenaltyMPC_z08, beliefPenaltyMPC_ubIdx08, params->ub9, beliefPenaltyMPC_lub08, beliefPenaltyMPC_sub08, beliefPenaltyMPC_riub08, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD3_233(params->lb10, beliefPenaltyMPC_z09, beliefPenaltyMPC_lbIdx09, beliefPenaltyMPC_llb09, beliefPenaltyMPC_slb09, beliefPenaltyMPC_rilb09, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD2_79(beliefPenaltyMPC_z09, beliefPenaltyMPC_ubIdx09, params->ub10, beliefPenaltyMPC_lub09, beliefPenaltyMPC_sub09, beliefPenaltyMPC_riub09, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD3_233(params->lb11, beliefPenaltyMPC_z10, beliefPenaltyMPC_lbIdx10, beliefPenaltyMPC_llb10, beliefPenaltyMPC_slb10, beliefPenaltyMPC_rilb10, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD2_79(beliefPenaltyMPC_z10, beliefPenaltyMPC_ubIdx10, params->ub11, beliefPenaltyMPC_lub10, beliefPenaltyMPC_sub10, beliefPenaltyMPC_riub10, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD3_233(params->lb12, beliefPenaltyMPC_z11, beliefPenaltyMPC_lbIdx11, beliefPenaltyMPC_llb11, beliefPenaltyMPC_slb11, beliefPenaltyMPC_rilb11, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD2_79(beliefPenaltyMPC_z11, beliefPenaltyMPC_ubIdx11, params->ub12, beliefPenaltyMPC_lub11, beliefPenaltyMPC_sub11, beliefPenaltyMPC_riub11, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD3_233(params->lb13, beliefPenaltyMPC_z12, beliefPenaltyMPC_lbIdx12, beliefPenaltyMPC_llb12, beliefPenaltyMPC_slb12, beliefPenaltyMPC_rilb12, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD2_79(beliefPenaltyMPC_z12, beliefPenaltyMPC_ubIdx12, params->ub13, beliefPenaltyMPC_lub12, beliefPenaltyMPC_sub12, beliefPenaltyMPC_riub12, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD3_233(params->lb14, beliefPenaltyMPC_z13, beliefPenaltyMPC_lbIdx13, beliefPenaltyMPC_llb13, beliefPenaltyMPC_slb13, beliefPenaltyMPC_rilb13, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD2_79(beliefPenaltyMPC_z13, beliefPenaltyMPC_ubIdx13, params->ub14, beliefPenaltyMPC_lub13, beliefPenaltyMPC_sub13, beliefPenaltyMPC_riub13, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD3_77(params->lb15, beliefPenaltyMPC_z14, beliefPenaltyMPC_lbIdx14, beliefPenaltyMPC_llb14, beliefPenaltyMPC_slb14, beliefPenaltyMPC_rilb14, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD2_77(beliefPenaltyMPC_z14, beliefPenaltyMPC_ubIdx14, params->ub15, beliefPenaltyMPC_lub14, beliefPenaltyMPC_sub14, beliefPenaltyMPC_riub14, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_INEQ_B_GRAD_233_233_79(beliefPenaltyMPC_lub00, beliefPenaltyMPC_sub00, beliefPenaltyMPC_riub00, beliefPenaltyMPC_llb00, beliefPenaltyMPC_slb00, beliefPenaltyMPC_rilb00, beliefPenaltyMPC_lbIdx00, beliefPenaltyMPC_ubIdx00, beliefPenaltyMPC_grad_ineq00, beliefPenaltyMPC_lubbysub00, beliefPenaltyMPC_llbbyslb00);
beliefPenaltyMPC_LA_INEQ_B_GRAD_233_233_79(beliefPenaltyMPC_lub01, beliefPenaltyMPC_sub01, beliefPenaltyMPC_riub01, beliefPenaltyMPC_llb01, beliefPenaltyMPC_slb01, beliefPenaltyMPC_rilb01, beliefPenaltyMPC_lbIdx01, beliefPenaltyMPC_ubIdx01, beliefPenaltyMPC_grad_ineq01, beliefPenaltyMPC_lubbysub01, beliefPenaltyMPC_llbbyslb01);
beliefPenaltyMPC_LA_INEQ_B_GRAD_233_233_79(beliefPenaltyMPC_lub02, beliefPenaltyMPC_sub02, beliefPenaltyMPC_riub02, beliefPenaltyMPC_llb02, beliefPenaltyMPC_slb02, beliefPenaltyMPC_rilb02, beliefPenaltyMPC_lbIdx02, beliefPenaltyMPC_ubIdx02, beliefPenaltyMPC_grad_ineq02, beliefPenaltyMPC_lubbysub02, beliefPenaltyMPC_llbbyslb02);
beliefPenaltyMPC_LA_INEQ_B_GRAD_233_233_79(beliefPenaltyMPC_lub03, beliefPenaltyMPC_sub03, beliefPenaltyMPC_riub03, beliefPenaltyMPC_llb03, beliefPenaltyMPC_slb03, beliefPenaltyMPC_rilb03, beliefPenaltyMPC_lbIdx03, beliefPenaltyMPC_ubIdx03, beliefPenaltyMPC_grad_ineq03, beliefPenaltyMPC_lubbysub03, beliefPenaltyMPC_llbbyslb03);
beliefPenaltyMPC_LA_INEQ_B_GRAD_233_233_79(beliefPenaltyMPC_lub04, beliefPenaltyMPC_sub04, beliefPenaltyMPC_riub04, beliefPenaltyMPC_llb04, beliefPenaltyMPC_slb04, beliefPenaltyMPC_rilb04, beliefPenaltyMPC_lbIdx04, beliefPenaltyMPC_ubIdx04, beliefPenaltyMPC_grad_ineq04, beliefPenaltyMPC_lubbysub04, beliefPenaltyMPC_llbbyslb04);
beliefPenaltyMPC_LA_INEQ_B_GRAD_233_233_79(beliefPenaltyMPC_lub05, beliefPenaltyMPC_sub05, beliefPenaltyMPC_riub05, beliefPenaltyMPC_llb05, beliefPenaltyMPC_slb05, beliefPenaltyMPC_rilb05, beliefPenaltyMPC_lbIdx05, beliefPenaltyMPC_ubIdx05, beliefPenaltyMPC_grad_ineq05, beliefPenaltyMPC_lubbysub05, beliefPenaltyMPC_llbbyslb05);
beliefPenaltyMPC_LA_INEQ_B_GRAD_233_233_79(beliefPenaltyMPC_lub06, beliefPenaltyMPC_sub06, beliefPenaltyMPC_riub06, beliefPenaltyMPC_llb06, beliefPenaltyMPC_slb06, beliefPenaltyMPC_rilb06, beliefPenaltyMPC_lbIdx06, beliefPenaltyMPC_ubIdx06, beliefPenaltyMPC_grad_ineq06, beliefPenaltyMPC_lubbysub06, beliefPenaltyMPC_llbbyslb06);
beliefPenaltyMPC_LA_INEQ_B_GRAD_233_233_79(beliefPenaltyMPC_lub07, beliefPenaltyMPC_sub07, beliefPenaltyMPC_riub07, beliefPenaltyMPC_llb07, beliefPenaltyMPC_slb07, beliefPenaltyMPC_rilb07, beliefPenaltyMPC_lbIdx07, beliefPenaltyMPC_ubIdx07, beliefPenaltyMPC_grad_ineq07, beliefPenaltyMPC_lubbysub07, beliefPenaltyMPC_llbbyslb07);
beliefPenaltyMPC_LA_INEQ_B_GRAD_233_233_79(beliefPenaltyMPC_lub08, beliefPenaltyMPC_sub08, beliefPenaltyMPC_riub08, beliefPenaltyMPC_llb08, beliefPenaltyMPC_slb08, beliefPenaltyMPC_rilb08, beliefPenaltyMPC_lbIdx08, beliefPenaltyMPC_ubIdx08, beliefPenaltyMPC_grad_ineq08, beliefPenaltyMPC_lubbysub08, beliefPenaltyMPC_llbbyslb08);
beliefPenaltyMPC_LA_INEQ_B_GRAD_233_233_79(beliefPenaltyMPC_lub09, beliefPenaltyMPC_sub09, beliefPenaltyMPC_riub09, beliefPenaltyMPC_llb09, beliefPenaltyMPC_slb09, beliefPenaltyMPC_rilb09, beliefPenaltyMPC_lbIdx09, beliefPenaltyMPC_ubIdx09, beliefPenaltyMPC_grad_ineq09, beliefPenaltyMPC_lubbysub09, beliefPenaltyMPC_llbbyslb09);
beliefPenaltyMPC_LA_INEQ_B_GRAD_233_233_79(beliefPenaltyMPC_lub10, beliefPenaltyMPC_sub10, beliefPenaltyMPC_riub10, beliefPenaltyMPC_llb10, beliefPenaltyMPC_slb10, beliefPenaltyMPC_rilb10, beliefPenaltyMPC_lbIdx10, beliefPenaltyMPC_ubIdx10, beliefPenaltyMPC_grad_ineq10, beliefPenaltyMPC_lubbysub10, beliefPenaltyMPC_llbbyslb10);
beliefPenaltyMPC_LA_INEQ_B_GRAD_233_233_79(beliefPenaltyMPC_lub11, beliefPenaltyMPC_sub11, beliefPenaltyMPC_riub11, beliefPenaltyMPC_llb11, beliefPenaltyMPC_slb11, beliefPenaltyMPC_rilb11, beliefPenaltyMPC_lbIdx11, beliefPenaltyMPC_ubIdx11, beliefPenaltyMPC_grad_ineq11, beliefPenaltyMPC_lubbysub11, beliefPenaltyMPC_llbbyslb11);
beliefPenaltyMPC_LA_INEQ_B_GRAD_233_233_79(beliefPenaltyMPC_lub12, beliefPenaltyMPC_sub12, beliefPenaltyMPC_riub12, beliefPenaltyMPC_llb12, beliefPenaltyMPC_slb12, beliefPenaltyMPC_rilb12, beliefPenaltyMPC_lbIdx12, beliefPenaltyMPC_ubIdx12, beliefPenaltyMPC_grad_ineq12, beliefPenaltyMPC_lubbysub12, beliefPenaltyMPC_llbbyslb12);
beliefPenaltyMPC_LA_INEQ_B_GRAD_233_233_79(beliefPenaltyMPC_lub13, beliefPenaltyMPC_sub13, beliefPenaltyMPC_riub13, beliefPenaltyMPC_llb13, beliefPenaltyMPC_slb13, beliefPenaltyMPC_rilb13, beliefPenaltyMPC_lbIdx13, beliefPenaltyMPC_ubIdx13, beliefPenaltyMPC_grad_ineq13, beliefPenaltyMPC_lubbysub13, beliefPenaltyMPC_llbbyslb13);
beliefPenaltyMPC_LA_INEQ_B_GRAD_77_77_77(beliefPenaltyMPC_lub14, beliefPenaltyMPC_sub14, beliefPenaltyMPC_riub14, beliefPenaltyMPC_llb14, beliefPenaltyMPC_slb14, beliefPenaltyMPC_rilb14, beliefPenaltyMPC_lbIdx14, beliefPenaltyMPC_ubIdx14, beliefPenaltyMPC_grad_ineq14, beliefPenaltyMPC_lubbysub14, beliefPenaltyMPC_llbbyslb14);
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
beliefPenaltyMPC_LA_VVADD3_3339(beliefPenaltyMPC_grad_cost, beliefPenaltyMPC_grad_eq, beliefPenaltyMPC_grad_ineq, beliefPenaltyMPC_rd);
beliefPenaltyMPC_LA_DIAG_CHOL_LBUB_233_233_79(params->H1, beliefPenaltyMPC_llbbyslb00, beliefPenaltyMPC_lbIdx00, beliefPenaltyMPC_lubbysub00, beliefPenaltyMPC_ubIdx00, beliefPenaltyMPC_Phi00);
beliefPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_77_233(beliefPenaltyMPC_Phi00, params->C1, beliefPenaltyMPC_V00);
beliefPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_77_233(beliefPenaltyMPC_Phi00, beliefPenaltyMPC_D00, beliefPenaltyMPC_W00);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMTM_77_233_77(beliefPenaltyMPC_W00, beliefPenaltyMPC_V00, beliefPenaltyMPC_Ysd01);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_233(beliefPenaltyMPC_Phi00, beliefPenaltyMPC_rd00, beliefPenaltyMPC_Lbyrd00);
beliefPenaltyMPC_LA_DIAG_CHOL_LBUB_233_233_79(params->H2, beliefPenaltyMPC_llbbyslb01, beliefPenaltyMPC_lbIdx01, beliefPenaltyMPC_lubbysub01, beliefPenaltyMPC_ubIdx01, beliefPenaltyMPC_Phi01);
beliefPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_77_233(beliefPenaltyMPC_Phi01, params->C2, beliefPenaltyMPC_V01);
beliefPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_77_233(beliefPenaltyMPC_Phi01, beliefPenaltyMPC_D01, beliefPenaltyMPC_W01);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMTM_77_233_77(beliefPenaltyMPC_W01, beliefPenaltyMPC_V01, beliefPenaltyMPC_Ysd02);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_233(beliefPenaltyMPC_Phi01, beliefPenaltyMPC_rd01, beliefPenaltyMPC_Lbyrd01);
beliefPenaltyMPC_LA_DIAG_CHOL_LBUB_233_233_79(params->H3, beliefPenaltyMPC_llbbyslb02, beliefPenaltyMPC_lbIdx02, beliefPenaltyMPC_lubbysub02, beliefPenaltyMPC_ubIdx02, beliefPenaltyMPC_Phi02);
beliefPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_77_233(beliefPenaltyMPC_Phi02, params->C3, beliefPenaltyMPC_V02);
beliefPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_77_233(beliefPenaltyMPC_Phi02, beliefPenaltyMPC_D01, beliefPenaltyMPC_W02);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMTM_77_233_77(beliefPenaltyMPC_W02, beliefPenaltyMPC_V02, beliefPenaltyMPC_Ysd03);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_233(beliefPenaltyMPC_Phi02, beliefPenaltyMPC_rd02, beliefPenaltyMPC_Lbyrd02);
beliefPenaltyMPC_LA_DIAG_CHOL_LBUB_233_233_79(params->H4, beliefPenaltyMPC_llbbyslb03, beliefPenaltyMPC_lbIdx03, beliefPenaltyMPC_lubbysub03, beliefPenaltyMPC_ubIdx03, beliefPenaltyMPC_Phi03);
beliefPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_77_233(beliefPenaltyMPC_Phi03, params->C4, beliefPenaltyMPC_V03);
beliefPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_77_233(beliefPenaltyMPC_Phi03, beliefPenaltyMPC_D01, beliefPenaltyMPC_W03);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMTM_77_233_77(beliefPenaltyMPC_W03, beliefPenaltyMPC_V03, beliefPenaltyMPC_Ysd04);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_233(beliefPenaltyMPC_Phi03, beliefPenaltyMPC_rd03, beliefPenaltyMPC_Lbyrd03);
beliefPenaltyMPC_LA_DIAG_CHOL_LBUB_233_233_79(params->H5, beliefPenaltyMPC_llbbyslb04, beliefPenaltyMPC_lbIdx04, beliefPenaltyMPC_lubbysub04, beliefPenaltyMPC_ubIdx04, beliefPenaltyMPC_Phi04);
beliefPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_77_233(beliefPenaltyMPC_Phi04, params->C5, beliefPenaltyMPC_V04);
beliefPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_77_233(beliefPenaltyMPC_Phi04, beliefPenaltyMPC_D01, beliefPenaltyMPC_W04);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMTM_77_233_77(beliefPenaltyMPC_W04, beliefPenaltyMPC_V04, beliefPenaltyMPC_Ysd05);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_233(beliefPenaltyMPC_Phi04, beliefPenaltyMPC_rd04, beliefPenaltyMPC_Lbyrd04);
beliefPenaltyMPC_LA_DIAG_CHOL_LBUB_233_233_79(params->H6, beliefPenaltyMPC_llbbyslb05, beliefPenaltyMPC_lbIdx05, beliefPenaltyMPC_lubbysub05, beliefPenaltyMPC_ubIdx05, beliefPenaltyMPC_Phi05);
beliefPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_77_233(beliefPenaltyMPC_Phi05, params->C6, beliefPenaltyMPC_V05);
beliefPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_77_233(beliefPenaltyMPC_Phi05, beliefPenaltyMPC_D01, beliefPenaltyMPC_W05);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMTM_77_233_77(beliefPenaltyMPC_W05, beliefPenaltyMPC_V05, beliefPenaltyMPC_Ysd06);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_233(beliefPenaltyMPC_Phi05, beliefPenaltyMPC_rd05, beliefPenaltyMPC_Lbyrd05);
beliefPenaltyMPC_LA_DIAG_CHOL_LBUB_233_233_79(params->H7, beliefPenaltyMPC_llbbyslb06, beliefPenaltyMPC_lbIdx06, beliefPenaltyMPC_lubbysub06, beliefPenaltyMPC_ubIdx06, beliefPenaltyMPC_Phi06);
beliefPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_77_233(beliefPenaltyMPC_Phi06, params->C7, beliefPenaltyMPC_V06);
beliefPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_77_233(beliefPenaltyMPC_Phi06, beliefPenaltyMPC_D01, beliefPenaltyMPC_W06);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMTM_77_233_77(beliefPenaltyMPC_W06, beliefPenaltyMPC_V06, beliefPenaltyMPC_Ysd07);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_233(beliefPenaltyMPC_Phi06, beliefPenaltyMPC_rd06, beliefPenaltyMPC_Lbyrd06);
beliefPenaltyMPC_LA_DIAG_CHOL_LBUB_233_233_79(params->H8, beliefPenaltyMPC_llbbyslb07, beliefPenaltyMPC_lbIdx07, beliefPenaltyMPC_lubbysub07, beliefPenaltyMPC_ubIdx07, beliefPenaltyMPC_Phi07);
beliefPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_77_233(beliefPenaltyMPC_Phi07, params->C8, beliefPenaltyMPC_V07);
beliefPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_77_233(beliefPenaltyMPC_Phi07, beliefPenaltyMPC_D01, beliefPenaltyMPC_W07);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMTM_77_233_77(beliefPenaltyMPC_W07, beliefPenaltyMPC_V07, beliefPenaltyMPC_Ysd08);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_233(beliefPenaltyMPC_Phi07, beliefPenaltyMPC_rd07, beliefPenaltyMPC_Lbyrd07);
beliefPenaltyMPC_LA_DIAG_CHOL_LBUB_233_233_79(params->H9, beliefPenaltyMPC_llbbyslb08, beliefPenaltyMPC_lbIdx08, beliefPenaltyMPC_lubbysub08, beliefPenaltyMPC_ubIdx08, beliefPenaltyMPC_Phi08);
beliefPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_77_233(beliefPenaltyMPC_Phi08, params->C9, beliefPenaltyMPC_V08);
beliefPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_77_233(beliefPenaltyMPC_Phi08, beliefPenaltyMPC_D01, beliefPenaltyMPC_W08);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMTM_77_233_77(beliefPenaltyMPC_W08, beliefPenaltyMPC_V08, beliefPenaltyMPC_Ysd09);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_233(beliefPenaltyMPC_Phi08, beliefPenaltyMPC_rd08, beliefPenaltyMPC_Lbyrd08);
beliefPenaltyMPC_LA_DIAG_CHOL_LBUB_233_233_79(params->H10, beliefPenaltyMPC_llbbyslb09, beliefPenaltyMPC_lbIdx09, beliefPenaltyMPC_lubbysub09, beliefPenaltyMPC_ubIdx09, beliefPenaltyMPC_Phi09);
beliefPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_77_233(beliefPenaltyMPC_Phi09, params->C10, beliefPenaltyMPC_V09);
beliefPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_77_233(beliefPenaltyMPC_Phi09, beliefPenaltyMPC_D01, beliefPenaltyMPC_W09);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMTM_77_233_77(beliefPenaltyMPC_W09, beliefPenaltyMPC_V09, beliefPenaltyMPC_Ysd10);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_233(beliefPenaltyMPC_Phi09, beliefPenaltyMPC_rd09, beliefPenaltyMPC_Lbyrd09);
beliefPenaltyMPC_LA_DIAG_CHOL_LBUB_233_233_79(params->H11, beliefPenaltyMPC_llbbyslb10, beliefPenaltyMPC_lbIdx10, beliefPenaltyMPC_lubbysub10, beliefPenaltyMPC_ubIdx10, beliefPenaltyMPC_Phi10);
beliefPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_77_233(beliefPenaltyMPC_Phi10, params->C11, beliefPenaltyMPC_V10);
beliefPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_77_233(beliefPenaltyMPC_Phi10, beliefPenaltyMPC_D01, beliefPenaltyMPC_W10);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMTM_77_233_77(beliefPenaltyMPC_W10, beliefPenaltyMPC_V10, beliefPenaltyMPC_Ysd11);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_233(beliefPenaltyMPC_Phi10, beliefPenaltyMPC_rd10, beliefPenaltyMPC_Lbyrd10);
beliefPenaltyMPC_LA_DIAG_CHOL_LBUB_233_233_79(params->H12, beliefPenaltyMPC_llbbyslb11, beliefPenaltyMPC_lbIdx11, beliefPenaltyMPC_lubbysub11, beliefPenaltyMPC_ubIdx11, beliefPenaltyMPC_Phi11);
beliefPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_77_233(beliefPenaltyMPC_Phi11, params->C12, beliefPenaltyMPC_V11);
beliefPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_77_233(beliefPenaltyMPC_Phi11, beliefPenaltyMPC_D01, beliefPenaltyMPC_W11);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMTM_77_233_77(beliefPenaltyMPC_W11, beliefPenaltyMPC_V11, beliefPenaltyMPC_Ysd12);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_233(beliefPenaltyMPC_Phi11, beliefPenaltyMPC_rd11, beliefPenaltyMPC_Lbyrd11);
beliefPenaltyMPC_LA_DIAG_CHOL_LBUB_233_233_79(params->H13, beliefPenaltyMPC_llbbyslb12, beliefPenaltyMPC_lbIdx12, beliefPenaltyMPC_lubbysub12, beliefPenaltyMPC_ubIdx12, beliefPenaltyMPC_Phi12);
beliefPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_77_233(beliefPenaltyMPC_Phi12, params->C13, beliefPenaltyMPC_V12);
beliefPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_77_233(beliefPenaltyMPC_Phi12, beliefPenaltyMPC_D01, beliefPenaltyMPC_W12);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMTM_77_233_77(beliefPenaltyMPC_W12, beliefPenaltyMPC_V12, beliefPenaltyMPC_Ysd13);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_233(beliefPenaltyMPC_Phi12, beliefPenaltyMPC_rd12, beliefPenaltyMPC_Lbyrd12);
beliefPenaltyMPC_LA_DIAG_CHOL_LBUB_233_233_79(params->H14, beliefPenaltyMPC_llbbyslb13, beliefPenaltyMPC_lbIdx13, beliefPenaltyMPC_lubbysub13, beliefPenaltyMPC_ubIdx13, beliefPenaltyMPC_Phi13);
beliefPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_77_233(beliefPenaltyMPC_Phi13, params->C14, beliefPenaltyMPC_V13);
beliefPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_77_233(beliefPenaltyMPC_Phi13, beliefPenaltyMPC_D01, beliefPenaltyMPC_W13);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMTM_77_233_77(beliefPenaltyMPC_W13, beliefPenaltyMPC_V13, beliefPenaltyMPC_Ysd14);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_233(beliefPenaltyMPC_Phi13, beliefPenaltyMPC_rd13, beliefPenaltyMPC_Lbyrd13);
beliefPenaltyMPC_LA_DIAG_CHOL_ONELOOP_LBUB_77_77_77(params->H15, beliefPenaltyMPC_llbbyslb14, beliefPenaltyMPC_lbIdx14, beliefPenaltyMPC_lubbysub14, beliefPenaltyMPC_ubIdx14, beliefPenaltyMPC_Phi14);
beliefPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_77_77(beliefPenaltyMPC_Phi14, beliefPenaltyMPC_D14, beliefPenaltyMPC_W14);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_77(beliefPenaltyMPC_Phi14, beliefPenaltyMPC_rd14, beliefPenaltyMPC_Lbyrd14);
beliefPenaltyMPC_LA_DIAGZERO_MMT_77(beliefPenaltyMPC_W00, beliefPenaltyMPC_Yd00);
beliefPenaltyMPC_LA_DIAGZERO_MVMSUB7_77(beliefPenaltyMPC_W00, beliefPenaltyMPC_Lbyrd00, beliefPenaltyMPC_re00, beliefPenaltyMPC_beta00);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMT2_77_233_233(beliefPenaltyMPC_V00, beliefPenaltyMPC_W01, beliefPenaltyMPC_Yd01);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_77_233_233(beliefPenaltyMPC_V00, beliefPenaltyMPC_Lbyrd00, beliefPenaltyMPC_W01, beliefPenaltyMPC_Lbyrd01, beliefPenaltyMPC_re01, beliefPenaltyMPC_beta01);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMT2_77_233_233(beliefPenaltyMPC_V01, beliefPenaltyMPC_W02, beliefPenaltyMPC_Yd02);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_77_233_233(beliefPenaltyMPC_V01, beliefPenaltyMPC_Lbyrd01, beliefPenaltyMPC_W02, beliefPenaltyMPC_Lbyrd02, beliefPenaltyMPC_re02, beliefPenaltyMPC_beta02);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMT2_77_233_233(beliefPenaltyMPC_V02, beliefPenaltyMPC_W03, beliefPenaltyMPC_Yd03);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_77_233_233(beliefPenaltyMPC_V02, beliefPenaltyMPC_Lbyrd02, beliefPenaltyMPC_W03, beliefPenaltyMPC_Lbyrd03, beliefPenaltyMPC_re03, beliefPenaltyMPC_beta03);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMT2_77_233_233(beliefPenaltyMPC_V03, beliefPenaltyMPC_W04, beliefPenaltyMPC_Yd04);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_77_233_233(beliefPenaltyMPC_V03, beliefPenaltyMPC_Lbyrd03, beliefPenaltyMPC_W04, beliefPenaltyMPC_Lbyrd04, beliefPenaltyMPC_re04, beliefPenaltyMPC_beta04);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMT2_77_233_233(beliefPenaltyMPC_V04, beliefPenaltyMPC_W05, beliefPenaltyMPC_Yd05);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_77_233_233(beliefPenaltyMPC_V04, beliefPenaltyMPC_Lbyrd04, beliefPenaltyMPC_W05, beliefPenaltyMPC_Lbyrd05, beliefPenaltyMPC_re05, beliefPenaltyMPC_beta05);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMT2_77_233_233(beliefPenaltyMPC_V05, beliefPenaltyMPC_W06, beliefPenaltyMPC_Yd06);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_77_233_233(beliefPenaltyMPC_V05, beliefPenaltyMPC_Lbyrd05, beliefPenaltyMPC_W06, beliefPenaltyMPC_Lbyrd06, beliefPenaltyMPC_re06, beliefPenaltyMPC_beta06);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMT2_77_233_233(beliefPenaltyMPC_V06, beliefPenaltyMPC_W07, beliefPenaltyMPC_Yd07);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_77_233_233(beliefPenaltyMPC_V06, beliefPenaltyMPC_Lbyrd06, beliefPenaltyMPC_W07, beliefPenaltyMPC_Lbyrd07, beliefPenaltyMPC_re07, beliefPenaltyMPC_beta07);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMT2_77_233_233(beliefPenaltyMPC_V07, beliefPenaltyMPC_W08, beliefPenaltyMPC_Yd08);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_77_233_233(beliefPenaltyMPC_V07, beliefPenaltyMPC_Lbyrd07, beliefPenaltyMPC_W08, beliefPenaltyMPC_Lbyrd08, beliefPenaltyMPC_re08, beliefPenaltyMPC_beta08);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMT2_77_233_233(beliefPenaltyMPC_V08, beliefPenaltyMPC_W09, beliefPenaltyMPC_Yd09);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_77_233_233(beliefPenaltyMPC_V08, beliefPenaltyMPC_Lbyrd08, beliefPenaltyMPC_W09, beliefPenaltyMPC_Lbyrd09, beliefPenaltyMPC_re09, beliefPenaltyMPC_beta09);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMT2_77_233_233(beliefPenaltyMPC_V09, beliefPenaltyMPC_W10, beliefPenaltyMPC_Yd10);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_77_233_233(beliefPenaltyMPC_V09, beliefPenaltyMPC_Lbyrd09, beliefPenaltyMPC_W10, beliefPenaltyMPC_Lbyrd10, beliefPenaltyMPC_re10, beliefPenaltyMPC_beta10);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMT2_77_233_233(beliefPenaltyMPC_V10, beliefPenaltyMPC_W11, beliefPenaltyMPC_Yd11);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_77_233_233(beliefPenaltyMPC_V10, beliefPenaltyMPC_Lbyrd10, beliefPenaltyMPC_W11, beliefPenaltyMPC_Lbyrd11, beliefPenaltyMPC_re11, beliefPenaltyMPC_beta11);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMT2_77_233_233(beliefPenaltyMPC_V11, beliefPenaltyMPC_W12, beliefPenaltyMPC_Yd12);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_77_233_233(beliefPenaltyMPC_V11, beliefPenaltyMPC_Lbyrd11, beliefPenaltyMPC_W12, beliefPenaltyMPC_Lbyrd12, beliefPenaltyMPC_re12, beliefPenaltyMPC_beta12);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMT2_77_233_233(beliefPenaltyMPC_V12, beliefPenaltyMPC_W13, beliefPenaltyMPC_Yd13);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_77_233_233(beliefPenaltyMPC_V12, beliefPenaltyMPC_Lbyrd12, beliefPenaltyMPC_W13, beliefPenaltyMPC_Lbyrd13, beliefPenaltyMPC_re13, beliefPenaltyMPC_beta13);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MMT2_77_233_77(beliefPenaltyMPC_V13, beliefPenaltyMPC_W14, beliefPenaltyMPC_Yd14);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_77_233_77(beliefPenaltyMPC_V13, beliefPenaltyMPC_Lbyrd13, beliefPenaltyMPC_W14, beliefPenaltyMPC_Lbyrd14, beliefPenaltyMPC_re14, beliefPenaltyMPC_beta14);
beliefPenaltyMPC_LA_DENSE_CHOL_77(beliefPenaltyMPC_Yd00, beliefPenaltyMPC_Ld00);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_77(beliefPenaltyMPC_Ld00, beliefPenaltyMPC_beta00, beliefPenaltyMPC_yy00);
beliefPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_77_77(beliefPenaltyMPC_Ld00, beliefPenaltyMPC_Ysd01, beliefPenaltyMPC_Lsd01);
beliefPenaltyMPC_LA_DENSE_MMTSUB_77_77(beliefPenaltyMPC_Lsd01, beliefPenaltyMPC_Yd01);
beliefPenaltyMPC_LA_DENSE_CHOL_77(beliefPenaltyMPC_Yd01, beliefPenaltyMPC_Ld01);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_77_77(beliefPenaltyMPC_Lsd01, beliefPenaltyMPC_yy00, beliefPenaltyMPC_beta01, beliefPenaltyMPC_bmy01);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_77(beliefPenaltyMPC_Ld01, beliefPenaltyMPC_bmy01, beliefPenaltyMPC_yy01);
beliefPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_77_77(beliefPenaltyMPC_Ld01, beliefPenaltyMPC_Ysd02, beliefPenaltyMPC_Lsd02);
beliefPenaltyMPC_LA_DENSE_MMTSUB_77_77(beliefPenaltyMPC_Lsd02, beliefPenaltyMPC_Yd02);
beliefPenaltyMPC_LA_DENSE_CHOL_77(beliefPenaltyMPC_Yd02, beliefPenaltyMPC_Ld02);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_77_77(beliefPenaltyMPC_Lsd02, beliefPenaltyMPC_yy01, beliefPenaltyMPC_beta02, beliefPenaltyMPC_bmy02);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_77(beliefPenaltyMPC_Ld02, beliefPenaltyMPC_bmy02, beliefPenaltyMPC_yy02);
beliefPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_77_77(beliefPenaltyMPC_Ld02, beliefPenaltyMPC_Ysd03, beliefPenaltyMPC_Lsd03);
beliefPenaltyMPC_LA_DENSE_MMTSUB_77_77(beliefPenaltyMPC_Lsd03, beliefPenaltyMPC_Yd03);
beliefPenaltyMPC_LA_DENSE_CHOL_77(beliefPenaltyMPC_Yd03, beliefPenaltyMPC_Ld03);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_77_77(beliefPenaltyMPC_Lsd03, beliefPenaltyMPC_yy02, beliefPenaltyMPC_beta03, beliefPenaltyMPC_bmy03);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_77(beliefPenaltyMPC_Ld03, beliefPenaltyMPC_bmy03, beliefPenaltyMPC_yy03);
beliefPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_77_77(beliefPenaltyMPC_Ld03, beliefPenaltyMPC_Ysd04, beliefPenaltyMPC_Lsd04);
beliefPenaltyMPC_LA_DENSE_MMTSUB_77_77(beliefPenaltyMPC_Lsd04, beliefPenaltyMPC_Yd04);
beliefPenaltyMPC_LA_DENSE_CHOL_77(beliefPenaltyMPC_Yd04, beliefPenaltyMPC_Ld04);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_77_77(beliefPenaltyMPC_Lsd04, beliefPenaltyMPC_yy03, beliefPenaltyMPC_beta04, beliefPenaltyMPC_bmy04);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_77(beliefPenaltyMPC_Ld04, beliefPenaltyMPC_bmy04, beliefPenaltyMPC_yy04);
beliefPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_77_77(beliefPenaltyMPC_Ld04, beliefPenaltyMPC_Ysd05, beliefPenaltyMPC_Lsd05);
beliefPenaltyMPC_LA_DENSE_MMTSUB_77_77(beliefPenaltyMPC_Lsd05, beliefPenaltyMPC_Yd05);
beliefPenaltyMPC_LA_DENSE_CHOL_77(beliefPenaltyMPC_Yd05, beliefPenaltyMPC_Ld05);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_77_77(beliefPenaltyMPC_Lsd05, beliefPenaltyMPC_yy04, beliefPenaltyMPC_beta05, beliefPenaltyMPC_bmy05);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_77(beliefPenaltyMPC_Ld05, beliefPenaltyMPC_bmy05, beliefPenaltyMPC_yy05);
beliefPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_77_77(beliefPenaltyMPC_Ld05, beliefPenaltyMPC_Ysd06, beliefPenaltyMPC_Lsd06);
beliefPenaltyMPC_LA_DENSE_MMTSUB_77_77(beliefPenaltyMPC_Lsd06, beliefPenaltyMPC_Yd06);
beliefPenaltyMPC_LA_DENSE_CHOL_77(beliefPenaltyMPC_Yd06, beliefPenaltyMPC_Ld06);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_77_77(beliefPenaltyMPC_Lsd06, beliefPenaltyMPC_yy05, beliefPenaltyMPC_beta06, beliefPenaltyMPC_bmy06);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_77(beliefPenaltyMPC_Ld06, beliefPenaltyMPC_bmy06, beliefPenaltyMPC_yy06);
beliefPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_77_77(beliefPenaltyMPC_Ld06, beliefPenaltyMPC_Ysd07, beliefPenaltyMPC_Lsd07);
beliefPenaltyMPC_LA_DENSE_MMTSUB_77_77(beliefPenaltyMPC_Lsd07, beliefPenaltyMPC_Yd07);
beliefPenaltyMPC_LA_DENSE_CHOL_77(beliefPenaltyMPC_Yd07, beliefPenaltyMPC_Ld07);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_77_77(beliefPenaltyMPC_Lsd07, beliefPenaltyMPC_yy06, beliefPenaltyMPC_beta07, beliefPenaltyMPC_bmy07);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_77(beliefPenaltyMPC_Ld07, beliefPenaltyMPC_bmy07, beliefPenaltyMPC_yy07);
beliefPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_77_77(beliefPenaltyMPC_Ld07, beliefPenaltyMPC_Ysd08, beliefPenaltyMPC_Lsd08);
beliefPenaltyMPC_LA_DENSE_MMTSUB_77_77(beliefPenaltyMPC_Lsd08, beliefPenaltyMPC_Yd08);
beliefPenaltyMPC_LA_DENSE_CHOL_77(beliefPenaltyMPC_Yd08, beliefPenaltyMPC_Ld08);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_77_77(beliefPenaltyMPC_Lsd08, beliefPenaltyMPC_yy07, beliefPenaltyMPC_beta08, beliefPenaltyMPC_bmy08);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_77(beliefPenaltyMPC_Ld08, beliefPenaltyMPC_bmy08, beliefPenaltyMPC_yy08);
beliefPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_77_77(beliefPenaltyMPC_Ld08, beliefPenaltyMPC_Ysd09, beliefPenaltyMPC_Lsd09);
beliefPenaltyMPC_LA_DENSE_MMTSUB_77_77(beliefPenaltyMPC_Lsd09, beliefPenaltyMPC_Yd09);
beliefPenaltyMPC_LA_DENSE_CHOL_77(beliefPenaltyMPC_Yd09, beliefPenaltyMPC_Ld09);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_77_77(beliefPenaltyMPC_Lsd09, beliefPenaltyMPC_yy08, beliefPenaltyMPC_beta09, beliefPenaltyMPC_bmy09);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_77(beliefPenaltyMPC_Ld09, beliefPenaltyMPC_bmy09, beliefPenaltyMPC_yy09);
beliefPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_77_77(beliefPenaltyMPC_Ld09, beliefPenaltyMPC_Ysd10, beliefPenaltyMPC_Lsd10);
beliefPenaltyMPC_LA_DENSE_MMTSUB_77_77(beliefPenaltyMPC_Lsd10, beliefPenaltyMPC_Yd10);
beliefPenaltyMPC_LA_DENSE_CHOL_77(beliefPenaltyMPC_Yd10, beliefPenaltyMPC_Ld10);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_77_77(beliefPenaltyMPC_Lsd10, beliefPenaltyMPC_yy09, beliefPenaltyMPC_beta10, beliefPenaltyMPC_bmy10);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_77(beliefPenaltyMPC_Ld10, beliefPenaltyMPC_bmy10, beliefPenaltyMPC_yy10);
beliefPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_77_77(beliefPenaltyMPC_Ld10, beliefPenaltyMPC_Ysd11, beliefPenaltyMPC_Lsd11);
beliefPenaltyMPC_LA_DENSE_MMTSUB_77_77(beliefPenaltyMPC_Lsd11, beliefPenaltyMPC_Yd11);
beliefPenaltyMPC_LA_DENSE_CHOL_77(beliefPenaltyMPC_Yd11, beliefPenaltyMPC_Ld11);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_77_77(beliefPenaltyMPC_Lsd11, beliefPenaltyMPC_yy10, beliefPenaltyMPC_beta11, beliefPenaltyMPC_bmy11);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_77(beliefPenaltyMPC_Ld11, beliefPenaltyMPC_bmy11, beliefPenaltyMPC_yy11);
beliefPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_77_77(beliefPenaltyMPC_Ld11, beliefPenaltyMPC_Ysd12, beliefPenaltyMPC_Lsd12);
beliefPenaltyMPC_LA_DENSE_MMTSUB_77_77(beliefPenaltyMPC_Lsd12, beliefPenaltyMPC_Yd12);
beliefPenaltyMPC_LA_DENSE_CHOL_77(beliefPenaltyMPC_Yd12, beliefPenaltyMPC_Ld12);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_77_77(beliefPenaltyMPC_Lsd12, beliefPenaltyMPC_yy11, beliefPenaltyMPC_beta12, beliefPenaltyMPC_bmy12);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_77(beliefPenaltyMPC_Ld12, beliefPenaltyMPC_bmy12, beliefPenaltyMPC_yy12);
beliefPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_77_77(beliefPenaltyMPC_Ld12, beliefPenaltyMPC_Ysd13, beliefPenaltyMPC_Lsd13);
beliefPenaltyMPC_LA_DENSE_MMTSUB_77_77(beliefPenaltyMPC_Lsd13, beliefPenaltyMPC_Yd13);
beliefPenaltyMPC_LA_DENSE_CHOL_77(beliefPenaltyMPC_Yd13, beliefPenaltyMPC_Ld13);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_77_77(beliefPenaltyMPC_Lsd13, beliefPenaltyMPC_yy12, beliefPenaltyMPC_beta13, beliefPenaltyMPC_bmy13);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_77(beliefPenaltyMPC_Ld13, beliefPenaltyMPC_bmy13, beliefPenaltyMPC_yy13);
beliefPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_77_77(beliefPenaltyMPC_Ld13, beliefPenaltyMPC_Ysd14, beliefPenaltyMPC_Lsd14);
beliefPenaltyMPC_LA_DENSE_MMTSUB_77_77(beliefPenaltyMPC_Lsd14, beliefPenaltyMPC_Yd14);
beliefPenaltyMPC_LA_DENSE_CHOL_77(beliefPenaltyMPC_Yd14, beliefPenaltyMPC_Ld14);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_77_77(beliefPenaltyMPC_Lsd14, beliefPenaltyMPC_yy13, beliefPenaltyMPC_beta14, beliefPenaltyMPC_bmy14);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_77(beliefPenaltyMPC_Ld14, beliefPenaltyMPC_bmy14, beliefPenaltyMPC_yy14);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_77(beliefPenaltyMPC_Ld14, beliefPenaltyMPC_yy14, beliefPenaltyMPC_dvaff14);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_77_77(beliefPenaltyMPC_Lsd14, beliefPenaltyMPC_dvaff14, beliefPenaltyMPC_yy13, beliefPenaltyMPC_bmy13);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_77(beliefPenaltyMPC_Ld13, beliefPenaltyMPC_bmy13, beliefPenaltyMPC_dvaff13);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_77_77(beliefPenaltyMPC_Lsd13, beliefPenaltyMPC_dvaff13, beliefPenaltyMPC_yy12, beliefPenaltyMPC_bmy12);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_77(beliefPenaltyMPC_Ld12, beliefPenaltyMPC_bmy12, beliefPenaltyMPC_dvaff12);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_77_77(beliefPenaltyMPC_Lsd12, beliefPenaltyMPC_dvaff12, beliefPenaltyMPC_yy11, beliefPenaltyMPC_bmy11);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_77(beliefPenaltyMPC_Ld11, beliefPenaltyMPC_bmy11, beliefPenaltyMPC_dvaff11);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_77_77(beliefPenaltyMPC_Lsd11, beliefPenaltyMPC_dvaff11, beliefPenaltyMPC_yy10, beliefPenaltyMPC_bmy10);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_77(beliefPenaltyMPC_Ld10, beliefPenaltyMPC_bmy10, beliefPenaltyMPC_dvaff10);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_77_77(beliefPenaltyMPC_Lsd10, beliefPenaltyMPC_dvaff10, beliefPenaltyMPC_yy09, beliefPenaltyMPC_bmy09);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_77(beliefPenaltyMPC_Ld09, beliefPenaltyMPC_bmy09, beliefPenaltyMPC_dvaff09);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_77_77(beliefPenaltyMPC_Lsd09, beliefPenaltyMPC_dvaff09, beliefPenaltyMPC_yy08, beliefPenaltyMPC_bmy08);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_77(beliefPenaltyMPC_Ld08, beliefPenaltyMPC_bmy08, beliefPenaltyMPC_dvaff08);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_77_77(beliefPenaltyMPC_Lsd08, beliefPenaltyMPC_dvaff08, beliefPenaltyMPC_yy07, beliefPenaltyMPC_bmy07);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_77(beliefPenaltyMPC_Ld07, beliefPenaltyMPC_bmy07, beliefPenaltyMPC_dvaff07);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_77_77(beliefPenaltyMPC_Lsd07, beliefPenaltyMPC_dvaff07, beliefPenaltyMPC_yy06, beliefPenaltyMPC_bmy06);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_77(beliefPenaltyMPC_Ld06, beliefPenaltyMPC_bmy06, beliefPenaltyMPC_dvaff06);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_77_77(beliefPenaltyMPC_Lsd06, beliefPenaltyMPC_dvaff06, beliefPenaltyMPC_yy05, beliefPenaltyMPC_bmy05);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_77(beliefPenaltyMPC_Ld05, beliefPenaltyMPC_bmy05, beliefPenaltyMPC_dvaff05);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_77_77(beliefPenaltyMPC_Lsd05, beliefPenaltyMPC_dvaff05, beliefPenaltyMPC_yy04, beliefPenaltyMPC_bmy04);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_77(beliefPenaltyMPC_Ld04, beliefPenaltyMPC_bmy04, beliefPenaltyMPC_dvaff04);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_77_77(beliefPenaltyMPC_Lsd04, beliefPenaltyMPC_dvaff04, beliefPenaltyMPC_yy03, beliefPenaltyMPC_bmy03);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_77(beliefPenaltyMPC_Ld03, beliefPenaltyMPC_bmy03, beliefPenaltyMPC_dvaff03);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_77_77(beliefPenaltyMPC_Lsd03, beliefPenaltyMPC_dvaff03, beliefPenaltyMPC_yy02, beliefPenaltyMPC_bmy02);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_77(beliefPenaltyMPC_Ld02, beliefPenaltyMPC_bmy02, beliefPenaltyMPC_dvaff02);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_77_77(beliefPenaltyMPC_Lsd02, beliefPenaltyMPC_dvaff02, beliefPenaltyMPC_yy01, beliefPenaltyMPC_bmy01);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_77(beliefPenaltyMPC_Ld01, beliefPenaltyMPC_bmy01, beliefPenaltyMPC_dvaff01);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_77_77(beliefPenaltyMPC_Lsd01, beliefPenaltyMPC_dvaff01, beliefPenaltyMPC_yy00, beliefPenaltyMPC_bmy00);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_77(beliefPenaltyMPC_Ld00, beliefPenaltyMPC_bmy00, beliefPenaltyMPC_dvaff00);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_77_233_77(params->C1, beliefPenaltyMPC_dvaff01, beliefPenaltyMPC_D00, beliefPenaltyMPC_dvaff00, beliefPenaltyMPC_grad_eq00);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_77_233_77(params->C2, beliefPenaltyMPC_dvaff02, beliefPenaltyMPC_D01, beliefPenaltyMPC_dvaff01, beliefPenaltyMPC_grad_eq01);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_77_233_77(params->C3, beliefPenaltyMPC_dvaff03, beliefPenaltyMPC_D01, beliefPenaltyMPC_dvaff02, beliefPenaltyMPC_grad_eq02);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_77_233_77(params->C4, beliefPenaltyMPC_dvaff04, beliefPenaltyMPC_D01, beliefPenaltyMPC_dvaff03, beliefPenaltyMPC_grad_eq03);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_77_233_77(params->C5, beliefPenaltyMPC_dvaff05, beliefPenaltyMPC_D01, beliefPenaltyMPC_dvaff04, beliefPenaltyMPC_grad_eq04);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_77_233_77(params->C6, beliefPenaltyMPC_dvaff06, beliefPenaltyMPC_D01, beliefPenaltyMPC_dvaff05, beliefPenaltyMPC_grad_eq05);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_77_233_77(params->C7, beliefPenaltyMPC_dvaff07, beliefPenaltyMPC_D01, beliefPenaltyMPC_dvaff06, beliefPenaltyMPC_grad_eq06);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_77_233_77(params->C8, beliefPenaltyMPC_dvaff08, beliefPenaltyMPC_D01, beliefPenaltyMPC_dvaff07, beliefPenaltyMPC_grad_eq07);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_77_233_77(params->C9, beliefPenaltyMPC_dvaff09, beliefPenaltyMPC_D01, beliefPenaltyMPC_dvaff08, beliefPenaltyMPC_grad_eq08);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_77_233_77(params->C10, beliefPenaltyMPC_dvaff10, beliefPenaltyMPC_D01, beliefPenaltyMPC_dvaff09, beliefPenaltyMPC_grad_eq09);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_77_233_77(params->C11, beliefPenaltyMPC_dvaff11, beliefPenaltyMPC_D01, beliefPenaltyMPC_dvaff10, beliefPenaltyMPC_grad_eq10);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_77_233_77(params->C12, beliefPenaltyMPC_dvaff12, beliefPenaltyMPC_D01, beliefPenaltyMPC_dvaff11, beliefPenaltyMPC_grad_eq11);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_77_233_77(params->C13, beliefPenaltyMPC_dvaff13, beliefPenaltyMPC_D01, beliefPenaltyMPC_dvaff12, beliefPenaltyMPC_grad_eq12);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_77_233_77(params->C14, beliefPenaltyMPC_dvaff14, beliefPenaltyMPC_D01, beliefPenaltyMPC_dvaff13, beliefPenaltyMPC_grad_eq13);
beliefPenaltyMPC_LA_DIAGZERO_MTVM_77_77(beliefPenaltyMPC_D14, beliefPenaltyMPC_dvaff14, beliefPenaltyMPC_grad_eq14);
beliefPenaltyMPC_LA_VSUB2_3339(beliefPenaltyMPC_rd, beliefPenaltyMPC_grad_eq, beliefPenaltyMPC_rd);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_233(beliefPenaltyMPC_Phi00, beliefPenaltyMPC_rd00, beliefPenaltyMPC_dzaff00);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_233(beliefPenaltyMPC_Phi01, beliefPenaltyMPC_rd01, beliefPenaltyMPC_dzaff01);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_233(beliefPenaltyMPC_Phi02, beliefPenaltyMPC_rd02, beliefPenaltyMPC_dzaff02);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_233(beliefPenaltyMPC_Phi03, beliefPenaltyMPC_rd03, beliefPenaltyMPC_dzaff03);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_233(beliefPenaltyMPC_Phi04, beliefPenaltyMPC_rd04, beliefPenaltyMPC_dzaff04);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_233(beliefPenaltyMPC_Phi05, beliefPenaltyMPC_rd05, beliefPenaltyMPC_dzaff05);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_233(beliefPenaltyMPC_Phi06, beliefPenaltyMPC_rd06, beliefPenaltyMPC_dzaff06);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_233(beliefPenaltyMPC_Phi07, beliefPenaltyMPC_rd07, beliefPenaltyMPC_dzaff07);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_233(beliefPenaltyMPC_Phi08, beliefPenaltyMPC_rd08, beliefPenaltyMPC_dzaff08);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_233(beliefPenaltyMPC_Phi09, beliefPenaltyMPC_rd09, beliefPenaltyMPC_dzaff09);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_233(beliefPenaltyMPC_Phi10, beliefPenaltyMPC_rd10, beliefPenaltyMPC_dzaff10);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_233(beliefPenaltyMPC_Phi11, beliefPenaltyMPC_rd11, beliefPenaltyMPC_dzaff11);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_233(beliefPenaltyMPC_Phi12, beliefPenaltyMPC_rd12, beliefPenaltyMPC_dzaff12);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_233(beliefPenaltyMPC_Phi13, beliefPenaltyMPC_rd13, beliefPenaltyMPC_dzaff13);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_77(beliefPenaltyMPC_Phi14, beliefPenaltyMPC_rd14, beliefPenaltyMPC_dzaff14);
beliefPenaltyMPC_LA_VSUB_INDEXED_233(beliefPenaltyMPC_dzaff00, beliefPenaltyMPC_lbIdx00, beliefPenaltyMPC_rilb00, beliefPenaltyMPC_dslbaff00);
beliefPenaltyMPC_LA_VSUB3_233(beliefPenaltyMPC_llbbyslb00, beliefPenaltyMPC_dslbaff00, beliefPenaltyMPC_llb00, beliefPenaltyMPC_dllbaff00);
beliefPenaltyMPC_LA_VSUB2_INDEXED_79(beliefPenaltyMPC_riub00, beliefPenaltyMPC_dzaff00, beliefPenaltyMPC_ubIdx00, beliefPenaltyMPC_dsubaff00);
beliefPenaltyMPC_LA_VSUB3_79(beliefPenaltyMPC_lubbysub00, beliefPenaltyMPC_dsubaff00, beliefPenaltyMPC_lub00, beliefPenaltyMPC_dlubaff00);
beliefPenaltyMPC_LA_VSUB_INDEXED_233(beliefPenaltyMPC_dzaff01, beliefPenaltyMPC_lbIdx01, beliefPenaltyMPC_rilb01, beliefPenaltyMPC_dslbaff01);
beliefPenaltyMPC_LA_VSUB3_233(beliefPenaltyMPC_llbbyslb01, beliefPenaltyMPC_dslbaff01, beliefPenaltyMPC_llb01, beliefPenaltyMPC_dllbaff01);
beliefPenaltyMPC_LA_VSUB2_INDEXED_79(beliefPenaltyMPC_riub01, beliefPenaltyMPC_dzaff01, beliefPenaltyMPC_ubIdx01, beliefPenaltyMPC_dsubaff01);
beliefPenaltyMPC_LA_VSUB3_79(beliefPenaltyMPC_lubbysub01, beliefPenaltyMPC_dsubaff01, beliefPenaltyMPC_lub01, beliefPenaltyMPC_dlubaff01);
beliefPenaltyMPC_LA_VSUB_INDEXED_233(beliefPenaltyMPC_dzaff02, beliefPenaltyMPC_lbIdx02, beliefPenaltyMPC_rilb02, beliefPenaltyMPC_dslbaff02);
beliefPenaltyMPC_LA_VSUB3_233(beliefPenaltyMPC_llbbyslb02, beliefPenaltyMPC_dslbaff02, beliefPenaltyMPC_llb02, beliefPenaltyMPC_dllbaff02);
beliefPenaltyMPC_LA_VSUB2_INDEXED_79(beliefPenaltyMPC_riub02, beliefPenaltyMPC_dzaff02, beliefPenaltyMPC_ubIdx02, beliefPenaltyMPC_dsubaff02);
beliefPenaltyMPC_LA_VSUB3_79(beliefPenaltyMPC_lubbysub02, beliefPenaltyMPC_dsubaff02, beliefPenaltyMPC_lub02, beliefPenaltyMPC_dlubaff02);
beliefPenaltyMPC_LA_VSUB_INDEXED_233(beliefPenaltyMPC_dzaff03, beliefPenaltyMPC_lbIdx03, beliefPenaltyMPC_rilb03, beliefPenaltyMPC_dslbaff03);
beliefPenaltyMPC_LA_VSUB3_233(beliefPenaltyMPC_llbbyslb03, beliefPenaltyMPC_dslbaff03, beliefPenaltyMPC_llb03, beliefPenaltyMPC_dllbaff03);
beliefPenaltyMPC_LA_VSUB2_INDEXED_79(beliefPenaltyMPC_riub03, beliefPenaltyMPC_dzaff03, beliefPenaltyMPC_ubIdx03, beliefPenaltyMPC_dsubaff03);
beliefPenaltyMPC_LA_VSUB3_79(beliefPenaltyMPC_lubbysub03, beliefPenaltyMPC_dsubaff03, beliefPenaltyMPC_lub03, beliefPenaltyMPC_dlubaff03);
beliefPenaltyMPC_LA_VSUB_INDEXED_233(beliefPenaltyMPC_dzaff04, beliefPenaltyMPC_lbIdx04, beliefPenaltyMPC_rilb04, beliefPenaltyMPC_dslbaff04);
beliefPenaltyMPC_LA_VSUB3_233(beliefPenaltyMPC_llbbyslb04, beliefPenaltyMPC_dslbaff04, beliefPenaltyMPC_llb04, beliefPenaltyMPC_dllbaff04);
beliefPenaltyMPC_LA_VSUB2_INDEXED_79(beliefPenaltyMPC_riub04, beliefPenaltyMPC_dzaff04, beliefPenaltyMPC_ubIdx04, beliefPenaltyMPC_dsubaff04);
beliefPenaltyMPC_LA_VSUB3_79(beliefPenaltyMPC_lubbysub04, beliefPenaltyMPC_dsubaff04, beliefPenaltyMPC_lub04, beliefPenaltyMPC_dlubaff04);
beliefPenaltyMPC_LA_VSUB_INDEXED_233(beliefPenaltyMPC_dzaff05, beliefPenaltyMPC_lbIdx05, beliefPenaltyMPC_rilb05, beliefPenaltyMPC_dslbaff05);
beliefPenaltyMPC_LA_VSUB3_233(beliefPenaltyMPC_llbbyslb05, beliefPenaltyMPC_dslbaff05, beliefPenaltyMPC_llb05, beliefPenaltyMPC_dllbaff05);
beliefPenaltyMPC_LA_VSUB2_INDEXED_79(beliefPenaltyMPC_riub05, beliefPenaltyMPC_dzaff05, beliefPenaltyMPC_ubIdx05, beliefPenaltyMPC_dsubaff05);
beliefPenaltyMPC_LA_VSUB3_79(beliefPenaltyMPC_lubbysub05, beliefPenaltyMPC_dsubaff05, beliefPenaltyMPC_lub05, beliefPenaltyMPC_dlubaff05);
beliefPenaltyMPC_LA_VSUB_INDEXED_233(beliefPenaltyMPC_dzaff06, beliefPenaltyMPC_lbIdx06, beliefPenaltyMPC_rilb06, beliefPenaltyMPC_dslbaff06);
beliefPenaltyMPC_LA_VSUB3_233(beliefPenaltyMPC_llbbyslb06, beliefPenaltyMPC_dslbaff06, beliefPenaltyMPC_llb06, beliefPenaltyMPC_dllbaff06);
beliefPenaltyMPC_LA_VSUB2_INDEXED_79(beliefPenaltyMPC_riub06, beliefPenaltyMPC_dzaff06, beliefPenaltyMPC_ubIdx06, beliefPenaltyMPC_dsubaff06);
beliefPenaltyMPC_LA_VSUB3_79(beliefPenaltyMPC_lubbysub06, beliefPenaltyMPC_dsubaff06, beliefPenaltyMPC_lub06, beliefPenaltyMPC_dlubaff06);
beliefPenaltyMPC_LA_VSUB_INDEXED_233(beliefPenaltyMPC_dzaff07, beliefPenaltyMPC_lbIdx07, beliefPenaltyMPC_rilb07, beliefPenaltyMPC_dslbaff07);
beliefPenaltyMPC_LA_VSUB3_233(beliefPenaltyMPC_llbbyslb07, beliefPenaltyMPC_dslbaff07, beliefPenaltyMPC_llb07, beliefPenaltyMPC_dllbaff07);
beliefPenaltyMPC_LA_VSUB2_INDEXED_79(beliefPenaltyMPC_riub07, beliefPenaltyMPC_dzaff07, beliefPenaltyMPC_ubIdx07, beliefPenaltyMPC_dsubaff07);
beliefPenaltyMPC_LA_VSUB3_79(beliefPenaltyMPC_lubbysub07, beliefPenaltyMPC_dsubaff07, beliefPenaltyMPC_lub07, beliefPenaltyMPC_dlubaff07);
beliefPenaltyMPC_LA_VSUB_INDEXED_233(beliefPenaltyMPC_dzaff08, beliefPenaltyMPC_lbIdx08, beliefPenaltyMPC_rilb08, beliefPenaltyMPC_dslbaff08);
beliefPenaltyMPC_LA_VSUB3_233(beliefPenaltyMPC_llbbyslb08, beliefPenaltyMPC_dslbaff08, beliefPenaltyMPC_llb08, beliefPenaltyMPC_dllbaff08);
beliefPenaltyMPC_LA_VSUB2_INDEXED_79(beliefPenaltyMPC_riub08, beliefPenaltyMPC_dzaff08, beliefPenaltyMPC_ubIdx08, beliefPenaltyMPC_dsubaff08);
beliefPenaltyMPC_LA_VSUB3_79(beliefPenaltyMPC_lubbysub08, beliefPenaltyMPC_dsubaff08, beliefPenaltyMPC_lub08, beliefPenaltyMPC_dlubaff08);
beliefPenaltyMPC_LA_VSUB_INDEXED_233(beliefPenaltyMPC_dzaff09, beliefPenaltyMPC_lbIdx09, beliefPenaltyMPC_rilb09, beliefPenaltyMPC_dslbaff09);
beliefPenaltyMPC_LA_VSUB3_233(beliefPenaltyMPC_llbbyslb09, beliefPenaltyMPC_dslbaff09, beliefPenaltyMPC_llb09, beliefPenaltyMPC_dllbaff09);
beliefPenaltyMPC_LA_VSUB2_INDEXED_79(beliefPenaltyMPC_riub09, beliefPenaltyMPC_dzaff09, beliefPenaltyMPC_ubIdx09, beliefPenaltyMPC_dsubaff09);
beliefPenaltyMPC_LA_VSUB3_79(beliefPenaltyMPC_lubbysub09, beliefPenaltyMPC_dsubaff09, beliefPenaltyMPC_lub09, beliefPenaltyMPC_dlubaff09);
beliefPenaltyMPC_LA_VSUB_INDEXED_233(beliefPenaltyMPC_dzaff10, beliefPenaltyMPC_lbIdx10, beliefPenaltyMPC_rilb10, beliefPenaltyMPC_dslbaff10);
beliefPenaltyMPC_LA_VSUB3_233(beliefPenaltyMPC_llbbyslb10, beliefPenaltyMPC_dslbaff10, beliefPenaltyMPC_llb10, beliefPenaltyMPC_dllbaff10);
beliefPenaltyMPC_LA_VSUB2_INDEXED_79(beliefPenaltyMPC_riub10, beliefPenaltyMPC_dzaff10, beliefPenaltyMPC_ubIdx10, beliefPenaltyMPC_dsubaff10);
beliefPenaltyMPC_LA_VSUB3_79(beliefPenaltyMPC_lubbysub10, beliefPenaltyMPC_dsubaff10, beliefPenaltyMPC_lub10, beliefPenaltyMPC_dlubaff10);
beliefPenaltyMPC_LA_VSUB_INDEXED_233(beliefPenaltyMPC_dzaff11, beliefPenaltyMPC_lbIdx11, beliefPenaltyMPC_rilb11, beliefPenaltyMPC_dslbaff11);
beliefPenaltyMPC_LA_VSUB3_233(beliefPenaltyMPC_llbbyslb11, beliefPenaltyMPC_dslbaff11, beliefPenaltyMPC_llb11, beliefPenaltyMPC_dllbaff11);
beliefPenaltyMPC_LA_VSUB2_INDEXED_79(beliefPenaltyMPC_riub11, beliefPenaltyMPC_dzaff11, beliefPenaltyMPC_ubIdx11, beliefPenaltyMPC_dsubaff11);
beliefPenaltyMPC_LA_VSUB3_79(beliefPenaltyMPC_lubbysub11, beliefPenaltyMPC_dsubaff11, beliefPenaltyMPC_lub11, beliefPenaltyMPC_dlubaff11);
beliefPenaltyMPC_LA_VSUB_INDEXED_233(beliefPenaltyMPC_dzaff12, beliefPenaltyMPC_lbIdx12, beliefPenaltyMPC_rilb12, beliefPenaltyMPC_dslbaff12);
beliefPenaltyMPC_LA_VSUB3_233(beliefPenaltyMPC_llbbyslb12, beliefPenaltyMPC_dslbaff12, beliefPenaltyMPC_llb12, beliefPenaltyMPC_dllbaff12);
beliefPenaltyMPC_LA_VSUB2_INDEXED_79(beliefPenaltyMPC_riub12, beliefPenaltyMPC_dzaff12, beliefPenaltyMPC_ubIdx12, beliefPenaltyMPC_dsubaff12);
beliefPenaltyMPC_LA_VSUB3_79(beliefPenaltyMPC_lubbysub12, beliefPenaltyMPC_dsubaff12, beliefPenaltyMPC_lub12, beliefPenaltyMPC_dlubaff12);
beliefPenaltyMPC_LA_VSUB_INDEXED_233(beliefPenaltyMPC_dzaff13, beliefPenaltyMPC_lbIdx13, beliefPenaltyMPC_rilb13, beliefPenaltyMPC_dslbaff13);
beliefPenaltyMPC_LA_VSUB3_233(beliefPenaltyMPC_llbbyslb13, beliefPenaltyMPC_dslbaff13, beliefPenaltyMPC_llb13, beliefPenaltyMPC_dllbaff13);
beliefPenaltyMPC_LA_VSUB2_INDEXED_79(beliefPenaltyMPC_riub13, beliefPenaltyMPC_dzaff13, beliefPenaltyMPC_ubIdx13, beliefPenaltyMPC_dsubaff13);
beliefPenaltyMPC_LA_VSUB3_79(beliefPenaltyMPC_lubbysub13, beliefPenaltyMPC_dsubaff13, beliefPenaltyMPC_lub13, beliefPenaltyMPC_dlubaff13);
beliefPenaltyMPC_LA_VSUB_INDEXED_77(beliefPenaltyMPC_dzaff14, beliefPenaltyMPC_lbIdx14, beliefPenaltyMPC_rilb14, beliefPenaltyMPC_dslbaff14);
beliefPenaltyMPC_LA_VSUB3_77(beliefPenaltyMPC_llbbyslb14, beliefPenaltyMPC_dslbaff14, beliefPenaltyMPC_llb14, beliefPenaltyMPC_dllbaff14);
beliefPenaltyMPC_LA_VSUB2_INDEXED_77(beliefPenaltyMPC_riub14, beliefPenaltyMPC_dzaff14, beliefPenaltyMPC_ubIdx14, beliefPenaltyMPC_dsubaff14);
beliefPenaltyMPC_LA_VSUB3_77(beliefPenaltyMPC_lubbysub14, beliefPenaltyMPC_dsubaff14, beliefPenaltyMPC_lub14, beliefPenaltyMPC_dlubaff14);
info->lsit_aff = beliefPenaltyMPC_LINESEARCH_BACKTRACKING_AFFINE(beliefPenaltyMPC_l, beliefPenaltyMPC_s, beliefPenaltyMPC_dl_aff, beliefPenaltyMPC_ds_aff, &info->step_aff, &info->mu_aff);
if( info->lsit_aff == beliefPenaltyMPC_NOPROGRESS ){
exitcode = beliefPenaltyMPC_NOPROGRESS; break;
}
sigma_3rdroot = info->mu_aff / info->mu;
info->sigma = sigma_3rdroot*sigma_3rdroot*sigma_3rdroot;
musigma = info->mu * info->sigma;
beliefPenaltyMPC_LA_VSUB5_4522(beliefPenaltyMPC_ds_aff, beliefPenaltyMPC_dl_aff, musigma, beliefPenaltyMPC_ccrhs);
beliefPenaltyMPC_LA_VSUB6_INDEXED_233_79_233(beliefPenaltyMPC_ccrhsub00, beliefPenaltyMPC_sub00, beliefPenaltyMPC_ubIdx00, beliefPenaltyMPC_ccrhsl00, beliefPenaltyMPC_slb00, beliefPenaltyMPC_lbIdx00, beliefPenaltyMPC_rd00);
beliefPenaltyMPC_LA_VSUB6_INDEXED_233_79_233(beliefPenaltyMPC_ccrhsub01, beliefPenaltyMPC_sub01, beliefPenaltyMPC_ubIdx01, beliefPenaltyMPC_ccrhsl01, beliefPenaltyMPC_slb01, beliefPenaltyMPC_lbIdx01, beliefPenaltyMPC_rd01);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_233(beliefPenaltyMPC_Phi00, beliefPenaltyMPC_rd00, beliefPenaltyMPC_Lbyrd00);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_233(beliefPenaltyMPC_Phi01, beliefPenaltyMPC_rd01, beliefPenaltyMPC_Lbyrd01);
beliefPenaltyMPC_LA_DIAGZERO_MVM_77(beliefPenaltyMPC_W00, beliefPenaltyMPC_Lbyrd00, beliefPenaltyMPC_beta00);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_77(beliefPenaltyMPC_Ld00, beliefPenaltyMPC_beta00, beliefPenaltyMPC_yy00);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_77_233_233(beliefPenaltyMPC_V00, beliefPenaltyMPC_Lbyrd00, beliefPenaltyMPC_W01, beliefPenaltyMPC_Lbyrd01, beliefPenaltyMPC_beta01);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_77_77(beliefPenaltyMPC_Lsd01, beliefPenaltyMPC_yy00, beliefPenaltyMPC_beta01, beliefPenaltyMPC_bmy01);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_77(beliefPenaltyMPC_Ld01, beliefPenaltyMPC_bmy01, beliefPenaltyMPC_yy01);
beliefPenaltyMPC_LA_VSUB6_INDEXED_233_79_233(beliefPenaltyMPC_ccrhsub02, beliefPenaltyMPC_sub02, beliefPenaltyMPC_ubIdx02, beliefPenaltyMPC_ccrhsl02, beliefPenaltyMPC_slb02, beliefPenaltyMPC_lbIdx02, beliefPenaltyMPC_rd02);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_233(beliefPenaltyMPC_Phi02, beliefPenaltyMPC_rd02, beliefPenaltyMPC_Lbyrd02);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_77_233_233(beliefPenaltyMPC_V01, beliefPenaltyMPC_Lbyrd01, beliefPenaltyMPC_W02, beliefPenaltyMPC_Lbyrd02, beliefPenaltyMPC_beta02);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_77_77(beliefPenaltyMPC_Lsd02, beliefPenaltyMPC_yy01, beliefPenaltyMPC_beta02, beliefPenaltyMPC_bmy02);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_77(beliefPenaltyMPC_Ld02, beliefPenaltyMPC_bmy02, beliefPenaltyMPC_yy02);
beliefPenaltyMPC_LA_VSUB6_INDEXED_233_79_233(beliefPenaltyMPC_ccrhsub03, beliefPenaltyMPC_sub03, beliefPenaltyMPC_ubIdx03, beliefPenaltyMPC_ccrhsl03, beliefPenaltyMPC_slb03, beliefPenaltyMPC_lbIdx03, beliefPenaltyMPC_rd03);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_233(beliefPenaltyMPC_Phi03, beliefPenaltyMPC_rd03, beliefPenaltyMPC_Lbyrd03);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_77_233_233(beliefPenaltyMPC_V02, beliefPenaltyMPC_Lbyrd02, beliefPenaltyMPC_W03, beliefPenaltyMPC_Lbyrd03, beliefPenaltyMPC_beta03);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_77_77(beliefPenaltyMPC_Lsd03, beliefPenaltyMPC_yy02, beliefPenaltyMPC_beta03, beliefPenaltyMPC_bmy03);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_77(beliefPenaltyMPC_Ld03, beliefPenaltyMPC_bmy03, beliefPenaltyMPC_yy03);
beliefPenaltyMPC_LA_VSUB6_INDEXED_233_79_233(beliefPenaltyMPC_ccrhsub04, beliefPenaltyMPC_sub04, beliefPenaltyMPC_ubIdx04, beliefPenaltyMPC_ccrhsl04, beliefPenaltyMPC_slb04, beliefPenaltyMPC_lbIdx04, beliefPenaltyMPC_rd04);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_233(beliefPenaltyMPC_Phi04, beliefPenaltyMPC_rd04, beliefPenaltyMPC_Lbyrd04);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_77_233_233(beliefPenaltyMPC_V03, beliefPenaltyMPC_Lbyrd03, beliefPenaltyMPC_W04, beliefPenaltyMPC_Lbyrd04, beliefPenaltyMPC_beta04);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_77_77(beliefPenaltyMPC_Lsd04, beliefPenaltyMPC_yy03, beliefPenaltyMPC_beta04, beliefPenaltyMPC_bmy04);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_77(beliefPenaltyMPC_Ld04, beliefPenaltyMPC_bmy04, beliefPenaltyMPC_yy04);
beliefPenaltyMPC_LA_VSUB6_INDEXED_233_79_233(beliefPenaltyMPC_ccrhsub05, beliefPenaltyMPC_sub05, beliefPenaltyMPC_ubIdx05, beliefPenaltyMPC_ccrhsl05, beliefPenaltyMPC_slb05, beliefPenaltyMPC_lbIdx05, beliefPenaltyMPC_rd05);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_233(beliefPenaltyMPC_Phi05, beliefPenaltyMPC_rd05, beliefPenaltyMPC_Lbyrd05);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_77_233_233(beliefPenaltyMPC_V04, beliefPenaltyMPC_Lbyrd04, beliefPenaltyMPC_W05, beliefPenaltyMPC_Lbyrd05, beliefPenaltyMPC_beta05);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_77_77(beliefPenaltyMPC_Lsd05, beliefPenaltyMPC_yy04, beliefPenaltyMPC_beta05, beliefPenaltyMPC_bmy05);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_77(beliefPenaltyMPC_Ld05, beliefPenaltyMPC_bmy05, beliefPenaltyMPC_yy05);
beliefPenaltyMPC_LA_VSUB6_INDEXED_233_79_233(beliefPenaltyMPC_ccrhsub06, beliefPenaltyMPC_sub06, beliefPenaltyMPC_ubIdx06, beliefPenaltyMPC_ccrhsl06, beliefPenaltyMPC_slb06, beliefPenaltyMPC_lbIdx06, beliefPenaltyMPC_rd06);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_233(beliefPenaltyMPC_Phi06, beliefPenaltyMPC_rd06, beliefPenaltyMPC_Lbyrd06);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_77_233_233(beliefPenaltyMPC_V05, beliefPenaltyMPC_Lbyrd05, beliefPenaltyMPC_W06, beliefPenaltyMPC_Lbyrd06, beliefPenaltyMPC_beta06);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_77_77(beliefPenaltyMPC_Lsd06, beliefPenaltyMPC_yy05, beliefPenaltyMPC_beta06, beliefPenaltyMPC_bmy06);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_77(beliefPenaltyMPC_Ld06, beliefPenaltyMPC_bmy06, beliefPenaltyMPC_yy06);
beliefPenaltyMPC_LA_VSUB6_INDEXED_233_79_233(beliefPenaltyMPC_ccrhsub07, beliefPenaltyMPC_sub07, beliefPenaltyMPC_ubIdx07, beliefPenaltyMPC_ccrhsl07, beliefPenaltyMPC_slb07, beliefPenaltyMPC_lbIdx07, beliefPenaltyMPC_rd07);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_233(beliefPenaltyMPC_Phi07, beliefPenaltyMPC_rd07, beliefPenaltyMPC_Lbyrd07);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_77_233_233(beliefPenaltyMPC_V06, beliefPenaltyMPC_Lbyrd06, beliefPenaltyMPC_W07, beliefPenaltyMPC_Lbyrd07, beliefPenaltyMPC_beta07);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_77_77(beliefPenaltyMPC_Lsd07, beliefPenaltyMPC_yy06, beliefPenaltyMPC_beta07, beliefPenaltyMPC_bmy07);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_77(beliefPenaltyMPC_Ld07, beliefPenaltyMPC_bmy07, beliefPenaltyMPC_yy07);
beliefPenaltyMPC_LA_VSUB6_INDEXED_233_79_233(beliefPenaltyMPC_ccrhsub08, beliefPenaltyMPC_sub08, beliefPenaltyMPC_ubIdx08, beliefPenaltyMPC_ccrhsl08, beliefPenaltyMPC_slb08, beliefPenaltyMPC_lbIdx08, beliefPenaltyMPC_rd08);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_233(beliefPenaltyMPC_Phi08, beliefPenaltyMPC_rd08, beliefPenaltyMPC_Lbyrd08);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_77_233_233(beliefPenaltyMPC_V07, beliefPenaltyMPC_Lbyrd07, beliefPenaltyMPC_W08, beliefPenaltyMPC_Lbyrd08, beliefPenaltyMPC_beta08);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_77_77(beliefPenaltyMPC_Lsd08, beliefPenaltyMPC_yy07, beliefPenaltyMPC_beta08, beliefPenaltyMPC_bmy08);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_77(beliefPenaltyMPC_Ld08, beliefPenaltyMPC_bmy08, beliefPenaltyMPC_yy08);
beliefPenaltyMPC_LA_VSUB6_INDEXED_233_79_233(beliefPenaltyMPC_ccrhsub09, beliefPenaltyMPC_sub09, beliefPenaltyMPC_ubIdx09, beliefPenaltyMPC_ccrhsl09, beliefPenaltyMPC_slb09, beliefPenaltyMPC_lbIdx09, beliefPenaltyMPC_rd09);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_233(beliefPenaltyMPC_Phi09, beliefPenaltyMPC_rd09, beliefPenaltyMPC_Lbyrd09);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_77_233_233(beliefPenaltyMPC_V08, beliefPenaltyMPC_Lbyrd08, beliefPenaltyMPC_W09, beliefPenaltyMPC_Lbyrd09, beliefPenaltyMPC_beta09);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_77_77(beliefPenaltyMPC_Lsd09, beliefPenaltyMPC_yy08, beliefPenaltyMPC_beta09, beliefPenaltyMPC_bmy09);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_77(beliefPenaltyMPC_Ld09, beliefPenaltyMPC_bmy09, beliefPenaltyMPC_yy09);
beliefPenaltyMPC_LA_VSUB6_INDEXED_233_79_233(beliefPenaltyMPC_ccrhsub10, beliefPenaltyMPC_sub10, beliefPenaltyMPC_ubIdx10, beliefPenaltyMPC_ccrhsl10, beliefPenaltyMPC_slb10, beliefPenaltyMPC_lbIdx10, beliefPenaltyMPC_rd10);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_233(beliefPenaltyMPC_Phi10, beliefPenaltyMPC_rd10, beliefPenaltyMPC_Lbyrd10);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_77_233_233(beliefPenaltyMPC_V09, beliefPenaltyMPC_Lbyrd09, beliefPenaltyMPC_W10, beliefPenaltyMPC_Lbyrd10, beliefPenaltyMPC_beta10);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_77_77(beliefPenaltyMPC_Lsd10, beliefPenaltyMPC_yy09, beliefPenaltyMPC_beta10, beliefPenaltyMPC_bmy10);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_77(beliefPenaltyMPC_Ld10, beliefPenaltyMPC_bmy10, beliefPenaltyMPC_yy10);
beliefPenaltyMPC_LA_VSUB6_INDEXED_233_79_233(beliefPenaltyMPC_ccrhsub11, beliefPenaltyMPC_sub11, beliefPenaltyMPC_ubIdx11, beliefPenaltyMPC_ccrhsl11, beliefPenaltyMPC_slb11, beliefPenaltyMPC_lbIdx11, beliefPenaltyMPC_rd11);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_233(beliefPenaltyMPC_Phi11, beliefPenaltyMPC_rd11, beliefPenaltyMPC_Lbyrd11);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_77_233_233(beliefPenaltyMPC_V10, beliefPenaltyMPC_Lbyrd10, beliefPenaltyMPC_W11, beliefPenaltyMPC_Lbyrd11, beliefPenaltyMPC_beta11);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_77_77(beliefPenaltyMPC_Lsd11, beliefPenaltyMPC_yy10, beliefPenaltyMPC_beta11, beliefPenaltyMPC_bmy11);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_77(beliefPenaltyMPC_Ld11, beliefPenaltyMPC_bmy11, beliefPenaltyMPC_yy11);
beliefPenaltyMPC_LA_VSUB6_INDEXED_233_79_233(beliefPenaltyMPC_ccrhsub12, beliefPenaltyMPC_sub12, beliefPenaltyMPC_ubIdx12, beliefPenaltyMPC_ccrhsl12, beliefPenaltyMPC_slb12, beliefPenaltyMPC_lbIdx12, beliefPenaltyMPC_rd12);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_233(beliefPenaltyMPC_Phi12, beliefPenaltyMPC_rd12, beliefPenaltyMPC_Lbyrd12);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_77_233_233(beliefPenaltyMPC_V11, beliefPenaltyMPC_Lbyrd11, beliefPenaltyMPC_W12, beliefPenaltyMPC_Lbyrd12, beliefPenaltyMPC_beta12);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_77_77(beliefPenaltyMPC_Lsd12, beliefPenaltyMPC_yy11, beliefPenaltyMPC_beta12, beliefPenaltyMPC_bmy12);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_77(beliefPenaltyMPC_Ld12, beliefPenaltyMPC_bmy12, beliefPenaltyMPC_yy12);
beliefPenaltyMPC_LA_VSUB6_INDEXED_233_79_233(beliefPenaltyMPC_ccrhsub13, beliefPenaltyMPC_sub13, beliefPenaltyMPC_ubIdx13, beliefPenaltyMPC_ccrhsl13, beliefPenaltyMPC_slb13, beliefPenaltyMPC_lbIdx13, beliefPenaltyMPC_rd13);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_233(beliefPenaltyMPC_Phi13, beliefPenaltyMPC_rd13, beliefPenaltyMPC_Lbyrd13);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_77_233_233(beliefPenaltyMPC_V12, beliefPenaltyMPC_Lbyrd12, beliefPenaltyMPC_W13, beliefPenaltyMPC_Lbyrd13, beliefPenaltyMPC_beta13);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_77_77(beliefPenaltyMPC_Lsd13, beliefPenaltyMPC_yy12, beliefPenaltyMPC_beta13, beliefPenaltyMPC_bmy13);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_77(beliefPenaltyMPC_Ld13, beliefPenaltyMPC_bmy13, beliefPenaltyMPC_yy13);
beliefPenaltyMPC_LA_VSUB6_INDEXED_77_77_77(beliefPenaltyMPC_ccrhsub14, beliefPenaltyMPC_sub14, beliefPenaltyMPC_ubIdx14, beliefPenaltyMPC_ccrhsl14, beliefPenaltyMPC_slb14, beliefPenaltyMPC_lbIdx14, beliefPenaltyMPC_rd14);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_77(beliefPenaltyMPC_Phi14, beliefPenaltyMPC_rd14, beliefPenaltyMPC_Lbyrd14);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_77_233_77(beliefPenaltyMPC_V13, beliefPenaltyMPC_Lbyrd13, beliefPenaltyMPC_W14, beliefPenaltyMPC_Lbyrd14, beliefPenaltyMPC_beta14);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_77_77(beliefPenaltyMPC_Lsd14, beliefPenaltyMPC_yy13, beliefPenaltyMPC_beta14, beliefPenaltyMPC_bmy14);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_77(beliefPenaltyMPC_Ld14, beliefPenaltyMPC_bmy14, beliefPenaltyMPC_yy14);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_77(beliefPenaltyMPC_Ld14, beliefPenaltyMPC_yy14, beliefPenaltyMPC_dvcc14);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_77_77(beliefPenaltyMPC_Lsd14, beliefPenaltyMPC_dvcc14, beliefPenaltyMPC_yy13, beliefPenaltyMPC_bmy13);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_77(beliefPenaltyMPC_Ld13, beliefPenaltyMPC_bmy13, beliefPenaltyMPC_dvcc13);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_77_77(beliefPenaltyMPC_Lsd13, beliefPenaltyMPC_dvcc13, beliefPenaltyMPC_yy12, beliefPenaltyMPC_bmy12);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_77(beliefPenaltyMPC_Ld12, beliefPenaltyMPC_bmy12, beliefPenaltyMPC_dvcc12);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_77_77(beliefPenaltyMPC_Lsd12, beliefPenaltyMPC_dvcc12, beliefPenaltyMPC_yy11, beliefPenaltyMPC_bmy11);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_77(beliefPenaltyMPC_Ld11, beliefPenaltyMPC_bmy11, beliefPenaltyMPC_dvcc11);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_77_77(beliefPenaltyMPC_Lsd11, beliefPenaltyMPC_dvcc11, beliefPenaltyMPC_yy10, beliefPenaltyMPC_bmy10);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_77(beliefPenaltyMPC_Ld10, beliefPenaltyMPC_bmy10, beliefPenaltyMPC_dvcc10);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_77_77(beliefPenaltyMPC_Lsd10, beliefPenaltyMPC_dvcc10, beliefPenaltyMPC_yy09, beliefPenaltyMPC_bmy09);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_77(beliefPenaltyMPC_Ld09, beliefPenaltyMPC_bmy09, beliefPenaltyMPC_dvcc09);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_77_77(beliefPenaltyMPC_Lsd09, beliefPenaltyMPC_dvcc09, beliefPenaltyMPC_yy08, beliefPenaltyMPC_bmy08);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_77(beliefPenaltyMPC_Ld08, beliefPenaltyMPC_bmy08, beliefPenaltyMPC_dvcc08);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_77_77(beliefPenaltyMPC_Lsd08, beliefPenaltyMPC_dvcc08, beliefPenaltyMPC_yy07, beliefPenaltyMPC_bmy07);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_77(beliefPenaltyMPC_Ld07, beliefPenaltyMPC_bmy07, beliefPenaltyMPC_dvcc07);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_77_77(beliefPenaltyMPC_Lsd07, beliefPenaltyMPC_dvcc07, beliefPenaltyMPC_yy06, beliefPenaltyMPC_bmy06);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_77(beliefPenaltyMPC_Ld06, beliefPenaltyMPC_bmy06, beliefPenaltyMPC_dvcc06);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_77_77(beliefPenaltyMPC_Lsd06, beliefPenaltyMPC_dvcc06, beliefPenaltyMPC_yy05, beliefPenaltyMPC_bmy05);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_77(beliefPenaltyMPC_Ld05, beliefPenaltyMPC_bmy05, beliefPenaltyMPC_dvcc05);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_77_77(beliefPenaltyMPC_Lsd05, beliefPenaltyMPC_dvcc05, beliefPenaltyMPC_yy04, beliefPenaltyMPC_bmy04);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_77(beliefPenaltyMPC_Ld04, beliefPenaltyMPC_bmy04, beliefPenaltyMPC_dvcc04);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_77_77(beliefPenaltyMPC_Lsd04, beliefPenaltyMPC_dvcc04, beliefPenaltyMPC_yy03, beliefPenaltyMPC_bmy03);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_77(beliefPenaltyMPC_Ld03, beliefPenaltyMPC_bmy03, beliefPenaltyMPC_dvcc03);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_77_77(beliefPenaltyMPC_Lsd03, beliefPenaltyMPC_dvcc03, beliefPenaltyMPC_yy02, beliefPenaltyMPC_bmy02);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_77(beliefPenaltyMPC_Ld02, beliefPenaltyMPC_bmy02, beliefPenaltyMPC_dvcc02);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_77_77(beliefPenaltyMPC_Lsd02, beliefPenaltyMPC_dvcc02, beliefPenaltyMPC_yy01, beliefPenaltyMPC_bmy01);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_77(beliefPenaltyMPC_Ld01, beliefPenaltyMPC_bmy01, beliefPenaltyMPC_dvcc01);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_77_77(beliefPenaltyMPC_Lsd01, beliefPenaltyMPC_dvcc01, beliefPenaltyMPC_yy00, beliefPenaltyMPC_bmy00);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_77(beliefPenaltyMPC_Ld00, beliefPenaltyMPC_bmy00, beliefPenaltyMPC_dvcc00);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_77_233_77(params->C1, beliefPenaltyMPC_dvcc01, beliefPenaltyMPC_D00, beliefPenaltyMPC_dvcc00, beliefPenaltyMPC_grad_eq00);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_77_233_77(params->C2, beliefPenaltyMPC_dvcc02, beliefPenaltyMPC_D01, beliefPenaltyMPC_dvcc01, beliefPenaltyMPC_grad_eq01);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_77_233_77(params->C3, beliefPenaltyMPC_dvcc03, beliefPenaltyMPC_D01, beliefPenaltyMPC_dvcc02, beliefPenaltyMPC_grad_eq02);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_77_233_77(params->C4, beliefPenaltyMPC_dvcc04, beliefPenaltyMPC_D01, beliefPenaltyMPC_dvcc03, beliefPenaltyMPC_grad_eq03);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_77_233_77(params->C5, beliefPenaltyMPC_dvcc05, beliefPenaltyMPC_D01, beliefPenaltyMPC_dvcc04, beliefPenaltyMPC_grad_eq04);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_77_233_77(params->C6, beliefPenaltyMPC_dvcc06, beliefPenaltyMPC_D01, beliefPenaltyMPC_dvcc05, beliefPenaltyMPC_grad_eq05);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_77_233_77(params->C7, beliefPenaltyMPC_dvcc07, beliefPenaltyMPC_D01, beliefPenaltyMPC_dvcc06, beliefPenaltyMPC_grad_eq06);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_77_233_77(params->C8, beliefPenaltyMPC_dvcc08, beliefPenaltyMPC_D01, beliefPenaltyMPC_dvcc07, beliefPenaltyMPC_grad_eq07);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_77_233_77(params->C9, beliefPenaltyMPC_dvcc09, beliefPenaltyMPC_D01, beliefPenaltyMPC_dvcc08, beliefPenaltyMPC_grad_eq08);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_77_233_77(params->C10, beliefPenaltyMPC_dvcc10, beliefPenaltyMPC_D01, beliefPenaltyMPC_dvcc09, beliefPenaltyMPC_grad_eq09);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_77_233_77(params->C11, beliefPenaltyMPC_dvcc11, beliefPenaltyMPC_D01, beliefPenaltyMPC_dvcc10, beliefPenaltyMPC_grad_eq10);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_77_233_77(params->C12, beliefPenaltyMPC_dvcc12, beliefPenaltyMPC_D01, beliefPenaltyMPC_dvcc11, beliefPenaltyMPC_grad_eq11);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_77_233_77(params->C13, beliefPenaltyMPC_dvcc13, beliefPenaltyMPC_D01, beliefPenaltyMPC_dvcc12, beliefPenaltyMPC_grad_eq12);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_77_233_77(params->C14, beliefPenaltyMPC_dvcc14, beliefPenaltyMPC_D01, beliefPenaltyMPC_dvcc13, beliefPenaltyMPC_grad_eq13);
beliefPenaltyMPC_LA_DIAGZERO_MTVM_77_77(beliefPenaltyMPC_D14, beliefPenaltyMPC_dvcc14, beliefPenaltyMPC_grad_eq14);
beliefPenaltyMPC_LA_VSUB_3339(beliefPenaltyMPC_rd, beliefPenaltyMPC_grad_eq, beliefPenaltyMPC_rd);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_233(beliefPenaltyMPC_Phi00, beliefPenaltyMPC_rd00, beliefPenaltyMPC_dzcc00);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_233(beliefPenaltyMPC_Phi01, beliefPenaltyMPC_rd01, beliefPenaltyMPC_dzcc01);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_233(beliefPenaltyMPC_Phi02, beliefPenaltyMPC_rd02, beliefPenaltyMPC_dzcc02);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_233(beliefPenaltyMPC_Phi03, beliefPenaltyMPC_rd03, beliefPenaltyMPC_dzcc03);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_233(beliefPenaltyMPC_Phi04, beliefPenaltyMPC_rd04, beliefPenaltyMPC_dzcc04);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_233(beliefPenaltyMPC_Phi05, beliefPenaltyMPC_rd05, beliefPenaltyMPC_dzcc05);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_233(beliefPenaltyMPC_Phi06, beliefPenaltyMPC_rd06, beliefPenaltyMPC_dzcc06);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_233(beliefPenaltyMPC_Phi07, beliefPenaltyMPC_rd07, beliefPenaltyMPC_dzcc07);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_233(beliefPenaltyMPC_Phi08, beliefPenaltyMPC_rd08, beliefPenaltyMPC_dzcc08);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_233(beliefPenaltyMPC_Phi09, beliefPenaltyMPC_rd09, beliefPenaltyMPC_dzcc09);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_233(beliefPenaltyMPC_Phi10, beliefPenaltyMPC_rd10, beliefPenaltyMPC_dzcc10);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_233(beliefPenaltyMPC_Phi11, beliefPenaltyMPC_rd11, beliefPenaltyMPC_dzcc11);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_233(beliefPenaltyMPC_Phi12, beliefPenaltyMPC_rd12, beliefPenaltyMPC_dzcc12);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_233(beliefPenaltyMPC_Phi13, beliefPenaltyMPC_rd13, beliefPenaltyMPC_dzcc13);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_77(beliefPenaltyMPC_Phi14, beliefPenaltyMPC_rd14, beliefPenaltyMPC_dzcc14);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_233(beliefPenaltyMPC_ccrhsl00, beliefPenaltyMPC_slb00, beliefPenaltyMPC_llbbyslb00, beliefPenaltyMPC_dzcc00, beliefPenaltyMPC_lbIdx00, beliefPenaltyMPC_dllbcc00);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_79(beliefPenaltyMPC_ccrhsub00, beliefPenaltyMPC_sub00, beliefPenaltyMPC_lubbysub00, beliefPenaltyMPC_dzcc00, beliefPenaltyMPC_ubIdx00, beliefPenaltyMPC_dlubcc00);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_233(beliefPenaltyMPC_ccrhsl01, beliefPenaltyMPC_slb01, beliefPenaltyMPC_llbbyslb01, beliefPenaltyMPC_dzcc01, beliefPenaltyMPC_lbIdx01, beliefPenaltyMPC_dllbcc01);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_79(beliefPenaltyMPC_ccrhsub01, beliefPenaltyMPC_sub01, beliefPenaltyMPC_lubbysub01, beliefPenaltyMPC_dzcc01, beliefPenaltyMPC_ubIdx01, beliefPenaltyMPC_dlubcc01);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_233(beliefPenaltyMPC_ccrhsl02, beliefPenaltyMPC_slb02, beliefPenaltyMPC_llbbyslb02, beliefPenaltyMPC_dzcc02, beliefPenaltyMPC_lbIdx02, beliefPenaltyMPC_dllbcc02);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_79(beliefPenaltyMPC_ccrhsub02, beliefPenaltyMPC_sub02, beliefPenaltyMPC_lubbysub02, beliefPenaltyMPC_dzcc02, beliefPenaltyMPC_ubIdx02, beliefPenaltyMPC_dlubcc02);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_233(beliefPenaltyMPC_ccrhsl03, beliefPenaltyMPC_slb03, beliefPenaltyMPC_llbbyslb03, beliefPenaltyMPC_dzcc03, beliefPenaltyMPC_lbIdx03, beliefPenaltyMPC_dllbcc03);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_79(beliefPenaltyMPC_ccrhsub03, beliefPenaltyMPC_sub03, beliefPenaltyMPC_lubbysub03, beliefPenaltyMPC_dzcc03, beliefPenaltyMPC_ubIdx03, beliefPenaltyMPC_dlubcc03);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_233(beliefPenaltyMPC_ccrhsl04, beliefPenaltyMPC_slb04, beliefPenaltyMPC_llbbyslb04, beliefPenaltyMPC_dzcc04, beliefPenaltyMPC_lbIdx04, beliefPenaltyMPC_dllbcc04);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_79(beliefPenaltyMPC_ccrhsub04, beliefPenaltyMPC_sub04, beliefPenaltyMPC_lubbysub04, beliefPenaltyMPC_dzcc04, beliefPenaltyMPC_ubIdx04, beliefPenaltyMPC_dlubcc04);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_233(beliefPenaltyMPC_ccrhsl05, beliefPenaltyMPC_slb05, beliefPenaltyMPC_llbbyslb05, beliefPenaltyMPC_dzcc05, beliefPenaltyMPC_lbIdx05, beliefPenaltyMPC_dllbcc05);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_79(beliefPenaltyMPC_ccrhsub05, beliefPenaltyMPC_sub05, beliefPenaltyMPC_lubbysub05, beliefPenaltyMPC_dzcc05, beliefPenaltyMPC_ubIdx05, beliefPenaltyMPC_dlubcc05);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_233(beliefPenaltyMPC_ccrhsl06, beliefPenaltyMPC_slb06, beliefPenaltyMPC_llbbyslb06, beliefPenaltyMPC_dzcc06, beliefPenaltyMPC_lbIdx06, beliefPenaltyMPC_dllbcc06);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_79(beliefPenaltyMPC_ccrhsub06, beliefPenaltyMPC_sub06, beliefPenaltyMPC_lubbysub06, beliefPenaltyMPC_dzcc06, beliefPenaltyMPC_ubIdx06, beliefPenaltyMPC_dlubcc06);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_233(beliefPenaltyMPC_ccrhsl07, beliefPenaltyMPC_slb07, beliefPenaltyMPC_llbbyslb07, beliefPenaltyMPC_dzcc07, beliefPenaltyMPC_lbIdx07, beliefPenaltyMPC_dllbcc07);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_79(beliefPenaltyMPC_ccrhsub07, beliefPenaltyMPC_sub07, beliefPenaltyMPC_lubbysub07, beliefPenaltyMPC_dzcc07, beliefPenaltyMPC_ubIdx07, beliefPenaltyMPC_dlubcc07);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_233(beliefPenaltyMPC_ccrhsl08, beliefPenaltyMPC_slb08, beliefPenaltyMPC_llbbyslb08, beliefPenaltyMPC_dzcc08, beliefPenaltyMPC_lbIdx08, beliefPenaltyMPC_dllbcc08);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_79(beliefPenaltyMPC_ccrhsub08, beliefPenaltyMPC_sub08, beliefPenaltyMPC_lubbysub08, beliefPenaltyMPC_dzcc08, beliefPenaltyMPC_ubIdx08, beliefPenaltyMPC_dlubcc08);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_233(beliefPenaltyMPC_ccrhsl09, beliefPenaltyMPC_slb09, beliefPenaltyMPC_llbbyslb09, beliefPenaltyMPC_dzcc09, beliefPenaltyMPC_lbIdx09, beliefPenaltyMPC_dllbcc09);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_79(beliefPenaltyMPC_ccrhsub09, beliefPenaltyMPC_sub09, beliefPenaltyMPC_lubbysub09, beliefPenaltyMPC_dzcc09, beliefPenaltyMPC_ubIdx09, beliefPenaltyMPC_dlubcc09);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_233(beliefPenaltyMPC_ccrhsl10, beliefPenaltyMPC_slb10, beliefPenaltyMPC_llbbyslb10, beliefPenaltyMPC_dzcc10, beliefPenaltyMPC_lbIdx10, beliefPenaltyMPC_dllbcc10);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_79(beliefPenaltyMPC_ccrhsub10, beliefPenaltyMPC_sub10, beliefPenaltyMPC_lubbysub10, beliefPenaltyMPC_dzcc10, beliefPenaltyMPC_ubIdx10, beliefPenaltyMPC_dlubcc10);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_233(beliefPenaltyMPC_ccrhsl11, beliefPenaltyMPC_slb11, beliefPenaltyMPC_llbbyslb11, beliefPenaltyMPC_dzcc11, beliefPenaltyMPC_lbIdx11, beliefPenaltyMPC_dllbcc11);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_79(beliefPenaltyMPC_ccrhsub11, beliefPenaltyMPC_sub11, beliefPenaltyMPC_lubbysub11, beliefPenaltyMPC_dzcc11, beliefPenaltyMPC_ubIdx11, beliefPenaltyMPC_dlubcc11);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_233(beliefPenaltyMPC_ccrhsl12, beliefPenaltyMPC_slb12, beliefPenaltyMPC_llbbyslb12, beliefPenaltyMPC_dzcc12, beliefPenaltyMPC_lbIdx12, beliefPenaltyMPC_dllbcc12);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_79(beliefPenaltyMPC_ccrhsub12, beliefPenaltyMPC_sub12, beliefPenaltyMPC_lubbysub12, beliefPenaltyMPC_dzcc12, beliefPenaltyMPC_ubIdx12, beliefPenaltyMPC_dlubcc12);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_233(beliefPenaltyMPC_ccrhsl13, beliefPenaltyMPC_slb13, beliefPenaltyMPC_llbbyslb13, beliefPenaltyMPC_dzcc13, beliefPenaltyMPC_lbIdx13, beliefPenaltyMPC_dllbcc13);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_79(beliefPenaltyMPC_ccrhsub13, beliefPenaltyMPC_sub13, beliefPenaltyMPC_lubbysub13, beliefPenaltyMPC_dzcc13, beliefPenaltyMPC_ubIdx13, beliefPenaltyMPC_dlubcc13);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_77(beliefPenaltyMPC_ccrhsl14, beliefPenaltyMPC_slb14, beliefPenaltyMPC_llbbyslb14, beliefPenaltyMPC_dzcc14, beliefPenaltyMPC_lbIdx14, beliefPenaltyMPC_dllbcc14);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_77(beliefPenaltyMPC_ccrhsub14, beliefPenaltyMPC_sub14, beliefPenaltyMPC_lubbysub14, beliefPenaltyMPC_dzcc14, beliefPenaltyMPC_ubIdx14, beliefPenaltyMPC_dlubcc14);
beliefPenaltyMPC_LA_VSUB7_4522(beliefPenaltyMPC_l, beliefPenaltyMPC_ccrhs, beliefPenaltyMPC_s, beliefPenaltyMPC_dl_cc, beliefPenaltyMPC_ds_cc);
beliefPenaltyMPC_LA_VADD_3339(beliefPenaltyMPC_dz_cc, beliefPenaltyMPC_dz_aff);
beliefPenaltyMPC_LA_VADD_1155(beliefPenaltyMPC_dv_cc, beliefPenaltyMPC_dv_aff);
beliefPenaltyMPC_LA_VADD_4522(beliefPenaltyMPC_dl_cc, beliefPenaltyMPC_dl_aff);
beliefPenaltyMPC_LA_VADD_4522(beliefPenaltyMPC_ds_cc, beliefPenaltyMPC_ds_aff);
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
output->z1[7] = beliefPenaltyMPC_z00[7];
output->z1[8] = beliefPenaltyMPC_z00[8];
output->z1[9] = beliefPenaltyMPC_z00[9];
output->z1[10] = beliefPenaltyMPC_z00[10];
output->z1[11] = beliefPenaltyMPC_z00[11];
output->z1[12] = beliefPenaltyMPC_z00[12];
output->z1[13] = beliefPenaltyMPC_z00[13];
output->z1[14] = beliefPenaltyMPC_z00[14];
output->z1[15] = beliefPenaltyMPC_z00[15];
output->z1[16] = beliefPenaltyMPC_z00[16];
output->z1[17] = beliefPenaltyMPC_z00[17];
output->z1[18] = beliefPenaltyMPC_z00[18];
output->z1[19] = beliefPenaltyMPC_z00[19];
output->z1[20] = beliefPenaltyMPC_z00[20];
output->z1[21] = beliefPenaltyMPC_z00[21];
output->z1[22] = beliefPenaltyMPC_z00[22];
output->z1[23] = beliefPenaltyMPC_z00[23];
output->z1[24] = beliefPenaltyMPC_z00[24];
output->z1[25] = beliefPenaltyMPC_z00[25];
output->z1[26] = beliefPenaltyMPC_z00[26];
output->z1[27] = beliefPenaltyMPC_z00[27];
output->z1[28] = beliefPenaltyMPC_z00[28];
output->z1[29] = beliefPenaltyMPC_z00[29];
output->z1[30] = beliefPenaltyMPC_z00[30];
output->z1[31] = beliefPenaltyMPC_z00[31];
output->z1[32] = beliefPenaltyMPC_z00[32];
output->z1[33] = beliefPenaltyMPC_z00[33];
output->z1[34] = beliefPenaltyMPC_z00[34];
output->z1[35] = beliefPenaltyMPC_z00[35];
output->z1[36] = beliefPenaltyMPC_z00[36];
output->z1[37] = beliefPenaltyMPC_z00[37];
output->z1[38] = beliefPenaltyMPC_z00[38];
output->z1[39] = beliefPenaltyMPC_z00[39];
output->z1[40] = beliefPenaltyMPC_z00[40];
output->z1[41] = beliefPenaltyMPC_z00[41];
output->z1[42] = beliefPenaltyMPC_z00[42];
output->z1[43] = beliefPenaltyMPC_z00[43];
output->z1[44] = beliefPenaltyMPC_z00[44];
output->z1[45] = beliefPenaltyMPC_z00[45];
output->z1[46] = beliefPenaltyMPC_z00[46];
output->z1[47] = beliefPenaltyMPC_z00[47];
output->z1[48] = beliefPenaltyMPC_z00[48];
output->z1[49] = beliefPenaltyMPC_z00[49];
output->z1[50] = beliefPenaltyMPC_z00[50];
output->z1[51] = beliefPenaltyMPC_z00[51];
output->z1[52] = beliefPenaltyMPC_z00[52];
output->z1[53] = beliefPenaltyMPC_z00[53];
output->z1[54] = beliefPenaltyMPC_z00[54];
output->z1[55] = beliefPenaltyMPC_z00[55];
output->z1[56] = beliefPenaltyMPC_z00[56];
output->z1[57] = beliefPenaltyMPC_z00[57];
output->z1[58] = beliefPenaltyMPC_z00[58];
output->z1[59] = beliefPenaltyMPC_z00[59];
output->z1[60] = beliefPenaltyMPC_z00[60];
output->z1[61] = beliefPenaltyMPC_z00[61];
output->z1[62] = beliefPenaltyMPC_z00[62];
output->z1[63] = beliefPenaltyMPC_z00[63];
output->z1[64] = beliefPenaltyMPC_z00[64];
output->z1[65] = beliefPenaltyMPC_z00[65];
output->z1[66] = beliefPenaltyMPC_z00[66];
output->z1[67] = beliefPenaltyMPC_z00[67];
output->z1[68] = beliefPenaltyMPC_z00[68];
output->z1[69] = beliefPenaltyMPC_z00[69];
output->z1[70] = beliefPenaltyMPC_z00[70];
output->z1[71] = beliefPenaltyMPC_z00[71];
output->z1[72] = beliefPenaltyMPC_z00[72];
output->z1[73] = beliefPenaltyMPC_z00[73];
output->z1[74] = beliefPenaltyMPC_z00[74];
output->z1[75] = beliefPenaltyMPC_z00[75];
output->z1[76] = beliefPenaltyMPC_z00[76];
output->z1[77] = beliefPenaltyMPC_z00[77];
output->z1[78] = beliefPenaltyMPC_z00[78];
output->z2[0] = beliefPenaltyMPC_z01[0];
output->z2[1] = beliefPenaltyMPC_z01[1];
output->z2[2] = beliefPenaltyMPC_z01[2];
output->z2[3] = beliefPenaltyMPC_z01[3];
output->z2[4] = beliefPenaltyMPC_z01[4];
output->z2[5] = beliefPenaltyMPC_z01[5];
output->z2[6] = beliefPenaltyMPC_z01[6];
output->z2[7] = beliefPenaltyMPC_z01[7];
output->z2[8] = beliefPenaltyMPC_z01[8];
output->z2[9] = beliefPenaltyMPC_z01[9];
output->z2[10] = beliefPenaltyMPC_z01[10];
output->z2[11] = beliefPenaltyMPC_z01[11];
output->z2[12] = beliefPenaltyMPC_z01[12];
output->z2[13] = beliefPenaltyMPC_z01[13];
output->z2[14] = beliefPenaltyMPC_z01[14];
output->z2[15] = beliefPenaltyMPC_z01[15];
output->z2[16] = beliefPenaltyMPC_z01[16];
output->z2[17] = beliefPenaltyMPC_z01[17];
output->z2[18] = beliefPenaltyMPC_z01[18];
output->z2[19] = beliefPenaltyMPC_z01[19];
output->z2[20] = beliefPenaltyMPC_z01[20];
output->z2[21] = beliefPenaltyMPC_z01[21];
output->z2[22] = beliefPenaltyMPC_z01[22];
output->z2[23] = beliefPenaltyMPC_z01[23];
output->z2[24] = beliefPenaltyMPC_z01[24];
output->z2[25] = beliefPenaltyMPC_z01[25];
output->z2[26] = beliefPenaltyMPC_z01[26];
output->z2[27] = beliefPenaltyMPC_z01[27];
output->z2[28] = beliefPenaltyMPC_z01[28];
output->z2[29] = beliefPenaltyMPC_z01[29];
output->z2[30] = beliefPenaltyMPC_z01[30];
output->z2[31] = beliefPenaltyMPC_z01[31];
output->z2[32] = beliefPenaltyMPC_z01[32];
output->z2[33] = beliefPenaltyMPC_z01[33];
output->z2[34] = beliefPenaltyMPC_z01[34];
output->z2[35] = beliefPenaltyMPC_z01[35];
output->z2[36] = beliefPenaltyMPC_z01[36];
output->z2[37] = beliefPenaltyMPC_z01[37];
output->z2[38] = beliefPenaltyMPC_z01[38];
output->z2[39] = beliefPenaltyMPC_z01[39];
output->z2[40] = beliefPenaltyMPC_z01[40];
output->z2[41] = beliefPenaltyMPC_z01[41];
output->z2[42] = beliefPenaltyMPC_z01[42];
output->z2[43] = beliefPenaltyMPC_z01[43];
output->z2[44] = beliefPenaltyMPC_z01[44];
output->z2[45] = beliefPenaltyMPC_z01[45];
output->z2[46] = beliefPenaltyMPC_z01[46];
output->z2[47] = beliefPenaltyMPC_z01[47];
output->z2[48] = beliefPenaltyMPC_z01[48];
output->z2[49] = beliefPenaltyMPC_z01[49];
output->z2[50] = beliefPenaltyMPC_z01[50];
output->z2[51] = beliefPenaltyMPC_z01[51];
output->z2[52] = beliefPenaltyMPC_z01[52];
output->z2[53] = beliefPenaltyMPC_z01[53];
output->z2[54] = beliefPenaltyMPC_z01[54];
output->z2[55] = beliefPenaltyMPC_z01[55];
output->z2[56] = beliefPenaltyMPC_z01[56];
output->z2[57] = beliefPenaltyMPC_z01[57];
output->z2[58] = beliefPenaltyMPC_z01[58];
output->z2[59] = beliefPenaltyMPC_z01[59];
output->z2[60] = beliefPenaltyMPC_z01[60];
output->z2[61] = beliefPenaltyMPC_z01[61];
output->z2[62] = beliefPenaltyMPC_z01[62];
output->z2[63] = beliefPenaltyMPC_z01[63];
output->z2[64] = beliefPenaltyMPC_z01[64];
output->z2[65] = beliefPenaltyMPC_z01[65];
output->z2[66] = beliefPenaltyMPC_z01[66];
output->z2[67] = beliefPenaltyMPC_z01[67];
output->z2[68] = beliefPenaltyMPC_z01[68];
output->z2[69] = beliefPenaltyMPC_z01[69];
output->z2[70] = beliefPenaltyMPC_z01[70];
output->z2[71] = beliefPenaltyMPC_z01[71];
output->z2[72] = beliefPenaltyMPC_z01[72];
output->z2[73] = beliefPenaltyMPC_z01[73];
output->z2[74] = beliefPenaltyMPC_z01[74];
output->z2[75] = beliefPenaltyMPC_z01[75];
output->z2[76] = beliefPenaltyMPC_z01[76];
output->z2[77] = beliefPenaltyMPC_z01[77];
output->z2[78] = beliefPenaltyMPC_z01[78];
output->z3[0] = beliefPenaltyMPC_z02[0];
output->z3[1] = beliefPenaltyMPC_z02[1];
output->z3[2] = beliefPenaltyMPC_z02[2];
output->z3[3] = beliefPenaltyMPC_z02[3];
output->z3[4] = beliefPenaltyMPC_z02[4];
output->z3[5] = beliefPenaltyMPC_z02[5];
output->z3[6] = beliefPenaltyMPC_z02[6];
output->z3[7] = beliefPenaltyMPC_z02[7];
output->z3[8] = beliefPenaltyMPC_z02[8];
output->z3[9] = beliefPenaltyMPC_z02[9];
output->z3[10] = beliefPenaltyMPC_z02[10];
output->z3[11] = beliefPenaltyMPC_z02[11];
output->z3[12] = beliefPenaltyMPC_z02[12];
output->z3[13] = beliefPenaltyMPC_z02[13];
output->z3[14] = beliefPenaltyMPC_z02[14];
output->z3[15] = beliefPenaltyMPC_z02[15];
output->z3[16] = beliefPenaltyMPC_z02[16];
output->z3[17] = beliefPenaltyMPC_z02[17];
output->z3[18] = beliefPenaltyMPC_z02[18];
output->z3[19] = beliefPenaltyMPC_z02[19];
output->z3[20] = beliefPenaltyMPC_z02[20];
output->z3[21] = beliefPenaltyMPC_z02[21];
output->z3[22] = beliefPenaltyMPC_z02[22];
output->z3[23] = beliefPenaltyMPC_z02[23];
output->z3[24] = beliefPenaltyMPC_z02[24];
output->z3[25] = beliefPenaltyMPC_z02[25];
output->z3[26] = beliefPenaltyMPC_z02[26];
output->z3[27] = beliefPenaltyMPC_z02[27];
output->z3[28] = beliefPenaltyMPC_z02[28];
output->z3[29] = beliefPenaltyMPC_z02[29];
output->z3[30] = beliefPenaltyMPC_z02[30];
output->z3[31] = beliefPenaltyMPC_z02[31];
output->z3[32] = beliefPenaltyMPC_z02[32];
output->z3[33] = beliefPenaltyMPC_z02[33];
output->z3[34] = beliefPenaltyMPC_z02[34];
output->z3[35] = beliefPenaltyMPC_z02[35];
output->z3[36] = beliefPenaltyMPC_z02[36];
output->z3[37] = beliefPenaltyMPC_z02[37];
output->z3[38] = beliefPenaltyMPC_z02[38];
output->z3[39] = beliefPenaltyMPC_z02[39];
output->z3[40] = beliefPenaltyMPC_z02[40];
output->z3[41] = beliefPenaltyMPC_z02[41];
output->z3[42] = beliefPenaltyMPC_z02[42];
output->z3[43] = beliefPenaltyMPC_z02[43];
output->z3[44] = beliefPenaltyMPC_z02[44];
output->z3[45] = beliefPenaltyMPC_z02[45];
output->z3[46] = beliefPenaltyMPC_z02[46];
output->z3[47] = beliefPenaltyMPC_z02[47];
output->z3[48] = beliefPenaltyMPC_z02[48];
output->z3[49] = beliefPenaltyMPC_z02[49];
output->z3[50] = beliefPenaltyMPC_z02[50];
output->z3[51] = beliefPenaltyMPC_z02[51];
output->z3[52] = beliefPenaltyMPC_z02[52];
output->z3[53] = beliefPenaltyMPC_z02[53];
output->z3[54] = beliefPenaltyMPC_z02[54];
output->z3[55] = beliefPenaltyMPC_z02[55];
output->z3[56] = beliefPenaltyMPC_z02[56];
output->z3[57] = beliefPenaltyMPC_z02[57];
output->z3[58] = beliefPenaltyMPC_z02[58];
output->z3[59] = beliefPenaltyMPC_z02[59];
output->z3[60] = beliefPenaltyMPC_z02[60];
output->z3[61] = beliefPenaltyMPC_z02[61];
output->z3[62] = beliefPenaltyMPC_z02[62];
output->z3[63] = beliefPenaltyMPC_z02[63];
output->z3[64] = beliefPenaltyMPC_z02[64];
output->z3[65] = beliefPenaltyMPC_z02[65];
output->z3[66] = beliefPenaltyMPC_z02[66];
output->z3[67] = beliefPenaltyMPC_z02[67];
output->z3[68] = beliefPenaltyMPC_z02[68];
output->z3[69] = beliefPenaltyMPC_z02[69];
output->z3[70] = beliefPenaltyMPC_z02[70];
output->z3[71] = beliefPenaltyMPC_z02[71];
output->z3[72] = beliefPenaltyMPC_z02[72];
output->z3[73] = beliefPenaltyMPC_z02[73];
output->z3[74] = beliefPenaltyMPC_z02[74];
output->z3[75] = beliefPenaltyMPC_z02[75];
output->z3[76] = beliefPenaltyMPC_z02[76];
output->z3[77] = beliefPenaltyMPC_z02[77];
output->z3[78] = beliefPenaltyMPC_z02[78];
output->z4[0] = beliefPenaltyMPC_z03[0];
output->z4[1] = beliefPenaltyMPC_z03[1];
output->z4[2] = beliefPenaltyMPC_z03[2];
output->z4[3] = beliefPenaltyMPC_z03[3];
output->z4[4] = beliefPenaltyMPC_z03[4];
output->z4[5] = beliefPenaltyMPC_z03[5];
output->z4[6] = beliefPenaltyMPC_z03[6];
output->z4[7] = beliefPenaltyMPC_z03[7];
output->z4[8] = beliefPenaltyMPC_z03[8];
output->z4[9] = beliefPenaltyMPC_z03[9];
output->z4[10] = beliefPenaltyMPC_z03[10];
output->z4[11] = beliefPenaltyMPC_z03[11];
output->z4[12] = beliefPenaltyMPC_z03[12];
output->z4[13] = beliefPenaltyMPC_z03[13];
output->z4[14] = beliefPenaltyMPC_z03[14];
output->z4[15] = beliefPenaltyMPC_z03[15];
output->z4[16] = beliefPenaltyMPC_z03[16];
output->z4[17] = beliefPenaltyMPC_z03[17];
output->z4[18] = beliefPenaltyMPC_z03[18];
output->z4[19] = beliefPenaltyMPC_z03[19];
output->z4[20] = beliefPenaltyMPC_z03[20];
output->z4[21] = beliefPenaltyMPC_z03[21];
output->z4[22] = beliefPenaltyMPC_z03[22];
output->z4[23] = beliefPenaltyMPC_z03[23];
output->z4[24] = beliefPenaltyMPC_z03[24];
output->z4[25] = beliefPenaltyMPC_z03[25];
output->z4[26] = beliefPenaltyMPC_z03[26];
output->z4[27] = beliefPenaltyMPC_z03[27];
output->z4[28] = beliefPenaltyMPC_z03[28];
output->z4[29] = beliefPenaltyMPC_z03[29];
output->z4[30] = beliefPenaltyMPC_z03[30];
output->z4[31] = beliefPenaltyMPC_z03[31];
output->z4[32] = beliefPenaltyMPC_z03[32];
output->z4[33] = beliefPenaltyMPC_z03[33];
output->z4[34] = beliefPenaltyMPC_z03[34];
output->z4[35] = beliefPenaltyMPC_z03[35];
output->z4[36] = beliefPenaltyMPC_z03[36];
output->z4[37] = beliefPenaltyMPC_z03[37];
output->z4[38] = beliefPenaltyMPC_z03[38];
output->z4[39] = beliefPenaltyMPC_z03[39];
output->z4[40] = beliefPenaltyMPC_z03[40];
output->z4[41] = beliefPenaltyMPC_z03[41];
output->z4[42] = beliefPenaltyMPC_z03[42];
output->z4[43] = beliefPenaltyMPC_z03[43];
output->z4[44] = beliefPenaltyMPC_z03[44];
output->z4[45] = beliefPenaltyMPC_z03[45];
output->z4[46] = beliefPenaltyMPC_z03[46];
output->z4[47] = beliefPenaltyMPC_z03[47];
output->z4[48] = beliefPenaltyMPC_z03[48];
output->z4[49] = beliefPenaltyMPC_z03[49];
output->z4[50] = beliefPenaltyMPC_z03[50];
output->z4[51] = beliefPenaltyMPC_z03[51];
output->z4[52] = beliefPenaltyMPC_z03[52];
output->z4[53] = beliefPenaltyMPC_z03[53];
output->z4[54] = beliefPenaltyMPC_z03[54];
output->z4[55] = beliefPenaltyMPC_z03[55];
output->z4[56] = beliefPenaltyMPC_z03[56];
output->z4[57] = beliefPenaltyMPC_z03[57];
output->z4[58] = beliefPenaltyMPC_z03[58];
output->z4[59] = beliefPenaltyMPC_z03[59];
output->z4[60] = beliefPenaltyMPC_z03[60];
output->z4[61] = beliefPenaltyMPC_z03[61];
output->z4[62] = beliefPenaltyMPC_z03[62];
output->z4[63] = beliefPenaltyMPC_z03[63];
output->z4[64] = beliefPenaltyMPC_z03[64];
output->z4[65] = beliefPenaltyMPC_z03[65];
output->z4[66] = beliefPenaltyMPC_z03[66];
output->z4[67] = beliefPenaltyMPC_z03[67];
output->z4[68] = beliefPenaltyMPC_z03[68];
output->z4[69] = beliefPenaltyMPC_z03[69];
output->z4[70] = beliefPenaltyMPC_z03[70];
output->z4[71] = beliefPenaltyMPC_z03[71];
output->z4[72] = beliefPenaltyMPC_z03[72];
output->z4[73] = beliefPenaltyMPC_z03[73];
output->z4[74] = beliefPenaltyMPC_z03[74];
output->z4[75] = beliefPenaltyMPC_z03[75];
output->z4[76] = beliefPenaltyMPC_z03[76];
output->z4[77] = beliefPenaltyMPC_z03[77];
output->z4[78] = beliefPenaltyMPC_z03[78];
output->z5[0] = beliefPenaltyMPC_z04[0];
output->z5[1] = beliefPenaltyMPC_z04[1];
output->z5[2] = beliefPenaltyMPC_z04[2];
output->z5[3] = beliefPenaltyMPC_z04[3];
output->z5[4] = beliefPenaltyMPC_z04[4];
output->z5[5] = beliefPenaltyMPC_z04[5];
output->z5[6] = beliefPenaltyMPC_z04[6];
output->z5[7] = beliefPenaltyMPC_z04[7];
output->z5[8] = beliefPenaltyMPC_z04[8];
output->z5[9] = beliefPenaltyMPC_z04[9];
output->z5[10] = beliefPenaltyMPC_z04[10];
output->z5[11] = beliefPenaltyMPC_z04[11];
output->z5[12] = beliefPenaltyMPC_z04[12];
output->z5[13] = beliefPenaltyMPC_z04[13];
output->z5[14] = beliefPenaltyMPC_z04[14];
output->z5[15] = beliefPenaltyMPC_z04[15];
output->z5[16] = beliefPenaltyMPC_z04[16];
output->z5[17] = beliefPenaltyMPC_z04[17];
output->z5[18] = beliefPenaltyMPC_z04[18];
output->z5[19] = beliefPenaltyMPC_z04[19];
output->z5[20] = beliefPenaltyMPC_z04[20];
output->z5[21] = beliefPenaltyMPC_z04[21];
output->z5[22] = beliefPenaltyMPC_z04[22];
output->z5[23] = beliefPenaltyMPC_z04[23];
output->z5[24] = beliefPenaltyMPC_z04[24];
output->z5[25] = beliefPenaltyMPC_z04[25];
output->z5[26] = beliefPenaltyMPC_z04[26];
output->z5[27] = beliefPenaltyMPC_z04[27];
output->z5[28] = beliefPenaltyMPC_z04[28];
output->z5[29] = beliefPenaltyMPC_z04[29];
output->z5[30] = beliefPenaltyMPC_z04[30];
output->z5[31] = beliefPenaltyMPC_z04[31];
output->z5[32] = beliefPenaltyMPC_z04[32];
output->z5[33] = beliefPenaltyMPC_z04[33];
output->z5[34] = beliefPenaltyMPC_z04[34];
output->z5[35] = beliefPenaltyMPC_z04[35];
output->z5[36] = beliefPenaltyMPC_z04[36];
output->z5[37] = beliefPenaltyMPC_z04[37];
output->z5[38] = beliefPenaltyMPC_z04[38];
output->z5[39] = beliefPenaltyMPC_z04[39];
output->z5[40] = beliefPenaltyMPC_z04[40];
output->z5[41] = beliefPenaltyMPC_z04[41];
output->z5[42] = beliefPenaltyMPC_z04[42];
output->z5[43] = beliefPenaltyMPC_z04[43];
output->z5[44] = beliefPenaltyMPC_z04[44];
output->z5[45] = beliefPenaltyMPC_z04[45];
output->z5[46] = beliefPenaltyMPC_z04[46];
output->z5[47] = beliefPenaltyMPC_z04[47];
output->z5[48] = beliefPenaltyMPC_z04[48];
output->z5[49] = beliefPenaltyMPC_z04[49];
output->z5[50] = beliefPenaltyMPC_z04[50];
output->z5[51] = beliefPenaltyMPC_z04[51];
output->z5[52] = beliefPenaltyMPC_z04[52];
output->z5[53] = beliefPenaltyMPC_z04[53];
output->z5[54] = beliefPenaltyMPC_z04[54];
output->z5[55] = beliefPenaltyMPC_z04[55];
output->z5[56] = beliefPenaltyMPC_z04[56];
output->z5[57] = beliefPenaltyMPC_z04[57];
output->z5[58] = beliefPenaltyMPC_z04[58];
output->z5[59] = beliefPenaltyMPC_z04[59];
output->z5[60] = beliefPenaltyMPC_z04[60];
output->z5[61] = beliefPenaltyMPC_z04[61];
output->z5[62] = beliefPenaltyMPC_z04[62];
output->z5[63] = beliefPenaltyMPC_z04[63];
output->z5[64] = beliefPenaltyMPC_z04[64];
output->z5[65] = beliefPenaltyMPC_z04[65];
output->z5[66] = beliefPenaltyMPC_z04[66];
output->z5[67] = beliefPenaltyMPC_z04[67];
output->z5[68] = beliefPenaltyMPC_z04[68];
output->z5[69] = beliefPenaltyMPC_z04[69];
output->z5[70] = beliefPenaltyMPC_z04[70];
output->z5[71] = beliefPenaltyMPC_z04[71];
output->z5[72] = beliefPenaltyMPC_z04[72];
output->z5[73] = beliefPenaltyMPC_z04[73];
output->z5[74] = beliefPenaltyMPC_z04[74];
output->z5[75] = beliefPenaltyMPC_z04[75];
output->z5[76] = beliefPenaltyMPC_z04[76];
output->z5[77] = beliefPenaltyMPC_z04[77];
output->z5[78] = beliefPenaltyMPC_z04[78];
output->z6[0] = beliefPenaltyMPC_z05[0];
output->z6[1] = beliefPenaltyMPC_z05[1];
output->z6[2] = beliefPenaltyMPC_z05[2];
output->z6[3] = beliefPenaltyMPC_z05[3];
output->z6[4] = beliefPenaltyMPC_z05[4];
output->z6[5] = beliefPenaltyMPC_z05[5];
output->z6[6] = beliefPenaltyMPC_z05[6];
output->z6[7] = beliefPenaltyMPC_z05[7];
output->z6[8] = beliefPenaltyMPC_z05[8];
output->z6[9] = beliefPenaltyMPC_z05[9];
output->z6[10] = beliefPenaltyMPC_z05[10];
output->z6[11] = beliefPenaltyMPC_z05[11];
output->z6[12] = beliefPenaltyMPC_z05[12];
output->z6[13] = beliefPenaltyMPC_z05[13];
output->z6[14] = beliefPenaltyMPC_z05[14];
output->z6[15] = beliefPenaltyMPC_z05[15];
output->z6[16] = beliefPenaltyMPC_z05[16];
output->z6[17] = beliefPenaltyMPC_z05[17];
output->z6[18] = beliefPenaltyMPC_z05[18];
output->z6[19] = beliefPenaltyMPC_z05[19];
output->z6[20] = beliefPenaltyMPC_z05[20];
output->z6[21] = beliefPenaltyMPC_z05[21];
output->z6[22] = beliefPenaltyMPC_z05[22];
output->z6[23] = beliefPenaltyMPC_z05[23];
output->z6[24] = beliefPenaltyMPC_z05[24];
output->z6[25] = beliefPenaltyMPC_z05[25];
output->z6[26] = beliefPenaltyMPC_z05[26];
output->z6[27] = beliefPenaltyMPC_z05[27];
output->z6[28] = beliefPenaltyMPC_z05[28];
output->z6[29] = beliefPenaltyMPC_z05[29];
output->z6[30] = beliefPenaltyMPC_z05[30];
output->z6[31] = beliefPenaltyMPC_z05[31];
output->z6[32] = beliefPenaltyMPC_z05[32];
output->z6[33] = beliefPenaltyMPC_z05[33];
output->z6[34] = beliefPenaltyMPC_z05[34];
output->z6[35] = beliefPenaltyMPC_z05[35];
output->z6[36] = beliefPenaltyMPC_z05[36];
output->z6[37] = beliefPenaltyMPC_z05[37];
output->z6[38] = beliefPenaltyMPC_z05[38];
output->z6[39] = beliefPenaltyMPC_z05[39];
output->z6[40] = beliefPenaltyMPC_z05[40];
output->z6[41] = beliefPenaltyMPC_z05[41];
output->z6[42] = beliefPenaltyMPC_z05[42];
output->z6[43] = beliefPenaltyMPC_z05[43];
output->z6[44] = beliefPenaltyMPC_z05[44];
output->z6[45] = beliefPenaltyMPC_z05[45];
output->z6[46] = beliefPenaltyMPC_z05[46];
output->z6[47] = beliefPenaltyMPC_z05[47];
output->z6[48] = beliefPenaltyMPC_z05[48];
output->z6[49] = beliefPenaltyMPC_z05[49];
output->z6[50] = beliefPenaltyMPC_z05[50];
output->z6[51] = beliefPenaltyMPC_z05[51];
output->z6[52] = beliefPenaltyMPC_z05[52];
output->z6[53] = beliefPenaltyMPC_z05[53];
output->z6[54] = beliefPenaltyMPC_z05[54];
output->z6[55] = beliefPenaltyMPC_z05[55];
output->z6[56] = beliefPenaltyMPC_z05[56];
output->z6[57] = beliefPenaltyMPC_z05[57];
output->z6[58] = beliefPenaltyMPC_z05[58];
output->z6[59] = beliefPenaltyMPC_z05[59];
output->z6[60] = beliefPenaltyMPC_z05[60];
output->z6[61] = beliefPenaltyMPC_z05[61];
output->z6[62] = beliefPenaltyMPC_z05[62];
output->z6[63] = beliefPenaltyMPC_z05[63];
output->z6[64] = beliefPenaltyMPC_z05[64];
output->z6[65] = beliefPenaltyMPC_z05[65];
output->z6[66] = beliefPenaltyMPC_z05[66];
output->z6[67] = beliefPenaltyMPC_z05[67];
output->z6[68] = beliefPenaltyMPC_z05[68];
output->z6[69] = beliefPenaltyMPC_z05[69];
output->z6[70] = beliefPenaltyMPC_z05[70];
output->z6[71] = beliefPenaltyMPC_z05[71];
output->z6[72] = beliefPenaltyMPC_z05[72];
output->z6[73] = beliefPenaltyMPC_z05[73];
output->z6[74] = beliefPenaltyMPC_z05[74];
output->z6[75] = beliefPenaltyMPC_z05[75];
output->z6[76] = beliefPenaltyMPC_z05[76];
output->z6[77] = beliefPenaltyMPC_z05[77];
output->z6[78] = beliefPenaltyMPC_z05[78];
output->z7[0] = beliefPenaltyMPC_z06[0];
output->z7[1] = beliefPenaltyMPC_z06[1];
output->z7[2] = beliefPenaltyMPC_z06[2];
output->z7[3] = beliefPenaltyMPC_z06[3];
output->z7[4] = beliefPenaltyMPC_z06[4];
output->z7[5] = beliefPenaltyMPC_z06[5];
output->z7[6] = beliefPenaltyMPC_z06[6];
output->z7[7] = beliefPenaltyMPC_z06[7];
output->z7[8] = beliefPenaltyMPC_z06[8];
output->z7[9] = beliefPenaltyMPC_z06[9];
output->z7[10] = beliefPenaltyMPC_z06[10];
output->z7[11] = beliefPenaltyMPC_z06[11];
output->z7[12] = beliefPenaltyMPC_z06[12];
output->z7[13] = beliefPenaltyMPC_z06[13];
output->z7[14] = beliefPenaltyMPC_z06[14];
output->z7[15] = beliefPenaltyMPC_z06[15];
output->z7[16] = beliefPenaltyMPC_z06[16];
output->z7[17] = beliefPenaltyMPC_z06[17];
output->z7[18] = beliefPenaltyMPC_z06[18];
output->z7[19] = beliefPenaltyMPC_z06[19];
output->z7[20] = beliefPenaltyMPC_z06[20];
output->z7[21] = beliefPenaltyMPC_z06[21];
output->z7[22] = beliefPenaltyMPC_z06[22];
output->z7[23] = beliefPenaltyMPC_z06[23];
output->z7[24] = beliefPenaltyMPC_z06[24];
output->z7[25] = beliefPenaltyMPC_z06[25];
output->z7[26] = beliefPenaltyMPC_z06[26];
output->z7[27] = beliefPenaltyMPC_z06[27];
output->z7[28] = beliefPenaltyMPC_z06[28];
output->z7[29] = beliefPenaltyMPC_z06[29];
output->z7[30] = beliefPenaltyMPC_z06[30];
output->z7[31] = beliefPenaltyMPC_z06[31];
output->z7[32] = beliefPenaltyMPC_z06[32];
output->z7[33] = beliefPenaltyMPC_z06[33];
output->z7[34] = beliefPenaltyMPC_z06[34];
output->z7[35] = beliefPenaltyMPC_z06[35];
output->z7[36] = beliefPenaltyMPC_z06[36];
output->z7[37] = beliefPenaltyMPC_z06[37];
output->z7[38] = beliefPenaltyMPC_z06[38];
output->z7[39] = beliefPenaltyMPC_z06[39];
output->z7[40] = beliefPenaltyMPC_z06[40];
output->z7[41] = beliefPenaltyMPC_z06[41];
output->z7[42] = beliefPenaltyMPC_z06[42];
output->z7[43] = beliefPenaltyMPC_z06[43];
output->z7[44] = beliefPenaltyMPC_z06[44];
output->z7[45] = beliefPenaltyMPC_z06[45];
output->z7[46] = beliefPenaltyMPC_z06[46];
output->z7[47] = beliefPenaltyMPC_z06[47];
output->z7[48] = beliefPenaltyMPC_z06[48];
output->z7[49] = beliefPenaltyMPC_z06[49];
output->z7[50] = beliefPenaltyMPC_z06[50];
output->z7[51] = beliefPenaltyMPC_z06[51];
output->z7[52] = beliefPenaltyMPC_z06[52];
output->z7[53] = beliefPenaltyMPC_z06[53];
output->z7[54] = beliefPenaltyMPC_z06[54];
output->z7[55] = beliefPenaltyMPC_z06[55];
output->z7[56] = beliefPenaltyMPC_z06[56];
output->z7[57] = beliefPenaltyMPC_z06[57];
output->z7[58] = beliefPenaltyMPC_z06[58];
output->z7[59] = beliefPenaltyMPC_z06[59];
output->z7[60] = beliefPenaltyMPC_z06[60];
output->z7[61] = beliefPenaltyMPC_z06[61];
output->z7[62] = beliefPenaltyMPC_z06[62];
output->z7[63] = beliefPenaltyMPC_z06[63];
output->z7[64] = beliefPenaltyMPC_z06[64];
output->z7[65] = beliefPenaltyMPC_z06[65];
output->z7[66] = beliefPenaltyMPC_z06[66];
output->z7[67] = beliefPenaltyMPC_z06[67];
output->z7[68] = beliefPenaltyMPC_z06[68];
output->z7[69] = beliefPenaltyMPC_z06[69];
output->z7[70] = beliefPenaltyMPC_z06[70];
output->z7[71] = beliefPenaltyMPC_z06[71];
output->z7[72] = beliefPenaltyMPC_z06[72];
output->z7[73] = beliefPenaltyMPC_z06[73];
output->z7[74] = beliefPenaltyMPC_z06[74];
output->z7[75] = beliefPenaltyMPC_z06[75];
output->z7[76] = beliefPenaltyMPC_z06[76];
output->z7[77] = beliefPenaltyMPC_z06[77];
output->z7[78] = beliefPenaltyMPC_z06[78];
output->z8[0] = beliefPenaltyMPC_z07[0];
output->z8[1] = beliefPenaltyMPC_z07[1];
output->z8[2] = beliefPenaltyMPC_z07[2];
output->z8[3] = beliefPenaltyMPC_z07[3];
output->z8[4] = beliefPenaltyMPC_z07[4];
output->z8[5] = beliefPenaltyMPC_z07[5];
output->z8[6] = beliefPenaltyMPC_z07[6];
output->z8[7] = beliefPenaltyMPC_z07[7];
output->z8[8] = beliefPenaltyMPC_z07[8];
output->z8[9] = beliefPenaltyMPC_z07[9];
output->z8[10] = beliefPenaltyMPC_z07[10];
output->z8[11] = beliefPenaltyMPC_z07[11];
output->z8[12] = beliefPenaltyMPC_z07[12];
output->z8[13] = beliefPenaltyMPC_z07[13];
output->z8[14] = beliefPenaltyMPC_z07[14];
output->z8[15] = beliefPenaltyMPC_z07[15];
output->z8[16] = beliefPenaltyMPC_z07[16];
output->z8[17] = beliefPenaltyMPC_z07[17];
output->z8[18] = beliefPenaltyMPC_z07[18];
output->z8[19] = beliefPenaltyMPC_z07[19];
output->z8[20] = beliefPenaltyMPC_z07[20];
output->z8[21] = beliefPenaltyMPC_z07[21];
output->z8[22] = beliefPenaltyMPC_z07[22];
output->z8[23] = beliefPenaltyMPC_z07[23];
output->z8[24] = beliefPenaltyMPC_z07[24];
output->z8[25] = beliefPenaltyMPC_z07[25];
output->z8[26] = beliefPenaltyMPC_z07[26];
output->z8[27] = beliefPenaltyMPC_z07[27];
output->z8[28] = beliefPenaltyMPC_z07[28];
output->z8[29] = beliefPenaltyMPC_z07[29];
output->z8[30] = beliefPenaltyMPC_z07[30];
output->z8[31] = beliefPenaltyMPC_z07[31];
output->z8[32] = beliefPenaltyMPC_z07[32];
output->z8[33] = beliefPenaltyMPC_z07[33];
output->z8[34] = beliefPenaltyMPC_z07[34];
output->z8[35] = beliefPenaltyMPC_z07[35];
output->z8[36] = beliefPenaltyMPC_z07[36];
output->z8[37] = beliefPenaltyMPC_z07[37];
output->z8[38] = beliefPenaltyMPC_z07[38];
output->z8[39] = beliefPenaltyMPC_z07[39];
output->z8[40] = beliefPenaltyMPC_z07[40];
output->z8[41] = beliefPenaltyMPC_z07[41];
output->z8[42] = beliefPenaltyMPC_z07[42];
output->z8[43] = beliefPenaltyMPC_z07[43];
output->z8[44] = beliefPenaltyMPC_z07[44];
output->z8[45] = beliefPenaltyMPC_z07[45];
output->z8[46] = beliefPenaltyMPC_z07[46];
output->z8[47] = beliefPenaltyMPC_z07[47];
output->z8[48] = beliefPenaltyMPC_z07[48];
output->z8[49] = beliefPenaltyMPC_z07[49];
output->z8[50] = beliefPenaltyMPC_z07[50];
output->z8[51] = beliefPenaltyMPC_z07[51];
output->z8[52] = beliefPenaltyMPC_z07[52];
output->z8[53] = beliefPenaltyMPC_z07[53];
output->z8[54] = beliefPenaltyMPC_z07[54];
output->z8[55] = beliefPenaltyMPC_z07[55];
output->z8[56] = beliefPenaltyMPC_z07[56];
output->z8[57] = beliefPenaltyMPC_z07[57];
output->z8[58] = beliefPenaltyMPC_z07[58];
output->z8[59] = beliefPenaltyMPC_z07[59];
output->z8[60] = beliefPenaltyMPC_z07[60];
output->z8[61] = beliefPenaltyMPC_z07[61];
output->z8[62] = beliefPenaltyMPC_z07[62];
output->z8[63] = beliefPenaltyMPC_z07[63];
output->z8[64] = beliefPenaltyMPC_z07[64];
output->z8[65] = beliefPenaltyMPC_z07[65];
output->z8[66] = beliefPenaltyMPC_z07[66];
output->z8[67] = beliefPenaltyMPC_z07[67];
output->z8[68] = beliefPenaltyMPC_z07[68];
output->z8[69] = beliefPenaltyMPC_z07[69];
output->z8[70] = beliefPenaltyMPC_z07[70];
output->z8[71] = beliefPenaltyMPC_z07[71];
output->z8[72] = beliefPenaltyMPC_z07[72];
output->z8[73] = beliefPenaltyMPC_z07[73];
output->z8[74] = beliefPenaltyMPC_z07[74];
output->z8[75] = beliefPenaltyMPC_z07[75];
output->z8[76] = beliefPenaltyMPC_z07[76];
output->z8[77] = beliefPenaltyMPC_z07[77];
output->z8[78] = beliefPenaltyMPC_z07[78];
output->z9[0] = beliefPenaltyMPC_z08[0];
output->z9[1] = beliefPenaltyMPC_z08[1];
output->z9[2] = beliefPenaltyMPC_z08[2];
output->z9[3] = beliefPenaltyMPC_z08[3];
output->z9[4] = beliefPenaltyMPC_z08[4];
output->z9[5] = beliefPenaltyMPC_z08[5];
output->z9[6] = beliefPenaltyMPC_z08[6];
output->z9[7] = beliefPenaltyMPC_z08[7];
output->z9[8] = beliefPenaltyMPC_z08[8];
output->z9[9] = beliefPenaltyMPC_z08[9];
output->z9[10] = beliefPenaltyMPC_z08[10];
output->z9[11] = beliefPenaltyMPC_z08[11];
output->z9[12] = beliefPenaltyMPC_z08[12];
output->z9[13] = beliefPenaltyMPC_z08[13];
output->z9[14] = beliefPenaltyMPC_z08[14];
output->z9[15] = beliefPenaltyMPC_z08[15];
output->z9[16] = beliefPenaltyMPC_z08[16];
output->z9[17] = beliefPenaltyMPC_z08[17];
output->z9[18] = beliefPenaltyMPC_z08[18];
output->z9[19] = beliefPenaltyMPC_z08[19];
output->z9[20] = beliefPenaltyMPC_z08[20];
output->z9[21] = beliefPenaltyMPC_z08[21];
output->z9[22] = beliefPenaltyMPC_z08[22];
output->z9[23] = beliefPenaltyMPC_z08[23];
output->z9[24] = beliefPenaltyMPC_z08[24];
output->z9[25] = beliefPenaltyMPC_z08[25];
output->z9[26] = beliefPenaltyMPC_z08[26];
output->z9[27] = beliefPenaltyMPC_z08[27];
output->z9[28] = beliefPenaltyMPC_z08[28];
output->z9[29] = beliefPenaltyMPC_z08[29];
output->z9[30] = beliefPenaltyMPC_z08[30];
output->z9[31] = beliefPenaltyMPC_z08[31];
output->z9[32] = beliefPenaltyMPC_z08[32];
output->z9[33] = beliefPenaltyMPC_z08[33];
output->z9[34] = beliefPenaltyMPC_z08[34];
output->z9[35] = beliefPenaltyMPC_z08[35];
output->z9[36] = beliefPenaltyMPC_z08[36];
output->z9[37] = beliefPenaltyMPC_z08[37];
output->z9[38] = beliefPenaltyMPC_z08[38];
output->z9[39] = beliefPenaltyMPC_z08[39];
output->z9[40] = beliefPenaltyMPC_z08[40];
output->z9[41] = beliefPenaltyMPC_z08[41];
output->z9[42] = beliefPenaltyMPC_z08[42];
output->z9[43] = beliefPenaltyMPC_z08[43];
output->z9[44] = beliefPenaltyMPC_z08[44];
output->z9[45] = beliefPenaltyMPC_z08[45];
output->z9[46] = beliefPenaltyMPC_z08[46];
output->z9[47] = beliefPenaltyMPC_z08[47];
output->z9[48] = beliefPenaltyMPC_z08[48];
output->z9[49] = beliefPenaltyMPC_z08[49];
output->z9[50] = beliefPenaltyMPC_z08[50];
output->z9[51] = beliefPenaltyMPC_z08[51];
output->z9[52] = beliefPenaltyMPC_z08[52];
output->z9[53] = beliefPenaltyMPC_z08[53];
output->z9[54] = beliefPenaltyMPC_z08[54];
output->z9[55] = beliefPenaltyMPC_z08[55];
output->z9[56] = beliefPenaltyMPC_z08[56];
output->z9[57] = beliefPenaltyMPC_z08[57];
output->z9[58] = beliefPenaltyMPC_z08[58];
output->z9[59] = beliefPenaltyMPC_z08[59];
output->z9[60] = beliefPenaltyMPC_z08[60];
output->z9[61] = beliefPenaltyMPC_z08[61];
output->z9[62] = beliefPenaltyMPC_z08[62];
output->z9[63] = beliefPenaltyMPC_z08[63];
output->z9[64] = beliefPenaltyMPC_z08[64];
output->z9[65] = beliefPenaltyMPC_z08[65];
output->z9[66] = beliefPenaltyMPC_z08[66];
output->z9[67] = beliefPenaltyMPC_z08[67];
output->z9[68] = beliefPenaltyMPC_z08[68];
output->z9[69] = beliefPenaltyMPC_z08[69];
output->z9[70] = beliefPenaltyMPC_z08[70];
output->z9[71] = beliefPenaltyMPC_z08[71];
output->z9[72] = beliefPenaltyMPC_z08[72];
output->z9[73] = beliefPenaltyMPC_z08[73];
output->z9[74] = beliefPenaltyMPC_z08[74];
output->z9[75] = beliefPenaltyMPC_z08[75];
output->z9[76] = beliefPenaltyMPC_z08[76];
output->z9[77] = beliefPenaltyMPC_z08[77];
output->z9[78] = beliefPenaltyMPC_z08[78];
output->z10[0] = beliefPenaltyMPC_z09[0];
output->z10[1] = beliefPenaltyMPC_z09[1];
output->z10[2] = beliefPenaltyMPC_z09[2];
output->z10[3] = beliefPenaltyMPC_z09[3];
output->z10[4] = beliefPenaltyMPC_z09[4];
output->z10[5] = beliefPenaltyMPC_z09[5];
output->z10[6] = beliefPenaltyMPC_z09[6];
output->z10[7] = beliefPenaltyMPC_z09[7];
output->z10[8] = beliefPenaltyMPC_z09[8];
output->z10[9] = beliefPenaltyMPC_z09[9];
output->z10[10] = beliefPenaltyMPC_z09[10];
output->z10[11] = beliefPenaltyMPC_z09[11];
output->z10[12] = beliefPenaltyMPC_z09[12];
output->z10[13] = beliefPenaltyMPC_z09[13];
output->z10[14] = beliefPenaltyMPC_z09[14];
output->z10[15] = beliefPenaltyMPC_z09[15];
output->z10[16] = beliefPenaltyMPC_z09[16];
output->z10[17] = beliefPenaltyMPC_z09[17];
output->z10[18] = beliefPenaltyMPC_z09[18];
output->z10[19] = beliefPenaltyMPC_z09[19];
output->z10[20] = beliefPenaltyMPC_z09[20];
output->z10[21] = beliefPenaltyMPC_z09[21];
output->z10[22] = beliefPenaltyMPC_z09[22];
output->z10[23] = beliefPenaltyMPC_z09[23];
output->z10[24] = beliefPenaltyMPC_z09[24];
output->z10[25] = beliefPenaltyMPC_z09[25];
output->z10[26] = beliefPenaltyMPC_z09[26];
output->z10[27] = beliefPenaltyMPC_z09[27];
output->z10[28] = beliefPenaltyMPC_z09[28];
output->z10[29] = beliefPenaltyMPC_z09[29];
output->z10[30] = beliefPenaltyMPC_z09[30];
output->z10[31] = beliefPenaltyMPC_z09[31];
output->z10[32] = beliefPenaltyMPC_z09[32];
output->z10[33] = beliefPenaltyMPC_z09[33];
output->z10[34] = beliefPenaltyMPC_z09[34];
output->z10[35] = beliefPenaltyMPC_z09[35];
output->z10[36] = beliefPenaltyMPC_z09[36];
output->z10[37] = beliefPenaltyMPC_z09[37];
output->z10[38] = beliefPenaltyMPC_z09[38];
output->z10[39] = beliefPenaltyMPC_z09[39];
output->z10[40] = beliefPenaltyMPC_z09[40];
output->z10[41] = beliefPenaltyMPC_z09[41];
output->z10[42] = beliefPenaltyMPC_z09[42];
output->z10[43] = beliefPenaltyMPC_z09[43];
output->z10[44] = beliefPenaltyMPC_z09[44];
output->z10[45] = beliefPenaltyMPC_z09[45];
output->z10[46] = beliefPenaltyMPC_z09[46];
output->z10[47] = beliefPenaltyMPC_z09[47];
output->z10[48] = beliefPenaltyMPC_z09[48];
output->z10[49] = beliefPenaltyMPC_z09[49];
output->z10[50] = beliefPenaltyMPC_z09[50];
output->z10[51] = beliefPenaltyMPC_z09[51];
output->z10[52] = beliefPenaltyMPC_z09[52];
output->z10[53] = beliefPenaltyMPC_z09[53];
output->z10[54] = beliefPenaltyMPC_z09[54];
output->z10[55] = beliefPenaltyMPC_z09[55];
output->z10[56] = beliefPenaltyMPC_z09[56];
output->z10[57] = beliefPenaltyMPC_z09[57];
output->z10[58] = beliefPenaltyMPC_z09[58];
output->z10[59] = beliefPenaltyMPC_z09[59];
output->z10[60] = beliefPenaltyMPC_z09[60];
output->z10[61] = beliefPenaltyMPC_z09[61];
output->z10[62] = beliefPenaltyMPC_z09[62];
output->z10[63] = beliefPenaltyMPC_z09[63];
output->z10[64] = beliefPenaltyMPC_z09[64];
output->z10[65] = beliefPenaltyMPC_z09[65];
output->z10[66] = beliefPenaltyMPC_z09[66];
output->z10[67] = beliefPenaltyMPC_z09[67];
output->z10[68] = beliefPenaltyMPC_z09[68];
output->z10[69] = beliefPenaltyMPC_z09[69];
output->z10[70] = beliefPenaltyMPC_z09[70];
output->z10[71] = beliefPenaltyMPC_z09[71];
output->z10[72] = beliefPenaltyMPC_z09[72];
output->z10[73] = beliefPenaltyMPC_z09[73];
output->z10[74] = beliefPenaltyMPC_z09[74];
output->z10[75] = beliefPenaltyMPC_z09[75];
output->z10[76] = beliefPenaltyMPC_z09[76];
output->z10[77] = beliefPenaltyMPC_z09[77];
output->z10[78] = beliefPenaltyMPC_z09[78];
output->z11[0] = beliefPenaltyMPC_z10[0];
output->z11[1] = beliefPenaltyMPC_z10[1];
output->z11[2] = beliefPenaltyMPC_z10[2];
output->z11[3] = beliefPenaltyMPC_z10[3];
output->z11[4] = beliefPenaltyMPC_z10[4];
output->z11[5] = beliefPenaltyMPC_z10[5];
output->z11[6] = beliefPenaltyMPC_z10[6];
output->z11[7] = beliefPenaltyMPC_z10[7];
output->z11[8] = beliefPenaltyMPC_z10[8];
output->z11[9] = beliefPenaltyMPC_z10[9];
output->z11[10] = beliefPenaltyMPC_z10[10];
output->z11[11] = beliefPenaltyMPC_z10[11];
output->z11[12] = beliefPenaltyMPC_z10[12];
output->z11[13] = beliefPenaltyMPC_z10[13];
output->z11[14] = beliefPenaltyMPC_z10[14];
output->z11[15] = beliefPenaltyMPC_z10[15];
output->z11[16] = beliefPenaltyMPC_z10[16];
output->z11[17] = beliefPenaltyMPC_z10[17];
output->z11[18] = beliefPenaltyMPC_z10[18];
output->z11[19] = beliefPenaltyMPC_z10[19];
output->z11[20] = beliefPenaltyMPC_z10[20];
output->z11[21] = beliefPenaltyMPC_z10[21];
output->z11[22] = beliefPenaltyMPC_z10[22];
output->z11[23] = beliefPenaltyMPC_z10[23];
output->z11[24] = beliefPenaltyMPC_z10[24];
output->z11[25] = beliefPenaltyMPC_z10[25];
output->z11[26] = beliefPenaltyMPC_z10[26];
output->z11[27] = beliefPenaltyMPC_z10[27];
output->z11[28] = beliefPenaltyMPC_z10[28];
output->z11[29] = beliefPenaltyMPC_z10[29];
output->z11[30] = beliefPenaltyMPC_z10[30];
output->z11[31] = beliefPenaltyMPC_z10[31];
output->z11[32] = beliefPenaltyMPC_z10[32];
output->z11[33] = beliefPenaltyMPC_z10[33];
output->z11[34] = beliefPenaltyMPC_z10[34];
output->z11[35] = beliefPenaltyMPC_z10[35];
output->z11[36] = beliefPenaltyMPC_z10[36];
output->z11[37] = beliefPenaltyMPC_z10[37];
output->z11[38] = beliefPenaltyMPC_z10[38];
output->z11[39] = beliefPenaltyMPC_z10[39];
output->z11[40] = beliefPenaltyMPC_z10[40];
output->z11[41] = beliefPenaltyMPC_z10[41];
output->z11[42] = beliefPenaltyMPC_z10[42];
output->z11[43] = beliefPenaltyMPC_z10[43];
output->z11[44] = beliefPenaltyMPC_z10[44];
output->z11[45] = beliefPenaltyMPC_z10[45];
output->z11[46] = beliefPenaltyMPC_z10[46];
output->z11[47] = beliefPenaltyMPC_z10[47];
output->z11[48] = beliefPenaltyMPC_z10[48];
output->z11[49] = beliefPenaltyMPC_z10[49];
output->z11[50] = beliefPenaltyMPC_z10[50];
output->z11[51] = beliefPenaltyMPC_z10[51];
output->z11[52] = beliefPenaltyMPC_z10[52];
output->z11[53] = beliefPenaltyMPC_z10[53];
output->z11[54] = beliefPenaltyMPC_z10[54];
output->z11[55] = beliefPenaltyMPC_z10[55];
output->z11[56] = beliefPenaltyMPC_z10[56];
output->z11[57] = beliefPenaltyMPC_z10[57];
output->z11[58] = beliefPenaltyMPC_z10[58];
output->z11[59] = beliefPenaltyMPC_z10[59];
output->z11[60] = beliefPenaltyMPC_z10[60];
output->z11[61] = beliefPenaltyMPC_z10[61];
output->z11[62] = beliefPenaltyMPC_z10[62];
output->z11[63] = beliefPenaltyMPC_z10[63];
output->z11[64] = beliefPenaltyMPC_z10[64];
output->z11[65] = beliefPenaltyMPC_z10[65];
output->z11[66] = beliefPenaltyMPC_z10[66];
output->z11[67] = beliefPenaltyMPC_z10[67];
output->z11[68] = beliefPenaltyMPC_z10[68];
output->z11[69] = beliefPenaltyMPC_z10[69];
output->z11[70] = beliefPenaltyMPC_z10[70];
output->z11[71] = beliefPenaltyMPC_z10[71];
output->z11[72] = beliefPenaltyMPC_z10[72];
output->z11[73] = beliefPenaltyMPC_z10[73];
output->z11[74] = beliefPenaltyMPC_z10[74];
output->z11[75] = beliefPenaltyMPC_z10[75];
output->z11[76] = beliefPenaltyMPC_z10[76];
output->z11[77] = beliefPenaltyMPC_z10[77];
output->z11[78] = beliefPenaltyMPC_z10[78];
output->z12[0] = beliefPenaltyMPC_z11[0];
output->z12[1] = beliefPenaltyMPC_z11[1];
output->z12[2] = beliefPenaltyMPC_z11[2];
output->z12[3] = beliefPenaltyMPC_z11[3];
output->z12[4] = beliefPenaltyMPC_z11[4];
output->z12[5] = beliefPenaltyMPC_z11[5];
output->z12[6] = beliefPenaltyMPC_z11[6];
output->z12[7] = beliefPenaltyMPC_z11[7];
output->z12[8] = beliefPenaltyMPC_z11[8];
output->z12[9] = beliefPenaltyMPC_z11[9];
output->z12[10] = beliefPenaltyMPC_z11[10];
output->z12[11] = beliefPenaltyMPC_z11[11];
output->z12[12] = beliefPenaltyMPC_z11[12];
output->z12[13] = beliefPenaltyMPC_z11[13];
output->z12[14] = beliefPenaltyMPC_z11[14];
output->z12[15] = beliefPenaltyMPC_z11[15];
output->z12[16] = beliefPenaltyMPC_z11[16];
output->z12[17] = beliefPenaltyMPC_z11[17];
output->z12[18] = beliefPenaltyMPC_z11[18];
output->z12[19] = beliefPenaltyMPC_z11[19];
output->z12[20] = beliefPenaltyMPC_z11[20];
output->z12[21] = beliefPenaltyMPC_z11[21];
output->z12[22] = beliefPenaltyMPC_z11[22];
output->z12[23] = beliefPenaltyMPC_z11[23];
output->z12[24] = beliefPenaltyMPC_z11[24];
output->z12[25] = beliefPenaltyMPC_z11[25];
output->z12[26] = beliefPenaltyMPC_z11[26];
output->z12[27] = beliefPenaltyMPC_z11[27];
output->z12[28] = beliefPenaltyMPC_z11[28];
output->z12[29] = beliefPenaltyMPC_z11[29];
output->z12[30] = beliefPenaltyMPC_z11[30];
output->z12[31] = beliefPenaltyMPC_z11[31];
output->z12[32] = beliefPenaltyMPC_z11[32];
output->z12[33] = beliefPenaltyMPC_z11[33];
output->z12[34] = beliefPenaltyMPC_z11[34];
output->z12[35] = beliefPenaltyMPC_z11[35];
output->z12[36] = beliefPenaltyMPC_z11[36];
output->z12[37] = beliefPenaltyMPC_z11[37];
output->z12[38] = beliefPenaltyMPC_z11[38];
output->z12[39] = beliefPenaltyMPC_z11[39];
output->z12[40] = beliefPenaltyMPC_z11[40];
output->z12[41] = beliefPenaltyMPC_z11[41];
output->z12[42] = beliefPenaltyMPC_z11[42];
output->z12[43] = beliefPenaltyMPC_z11[43];
output->z12[44] = beliefPenaltyMPC_z11[44];
output->z12[45] = beliefPenaltyMPC_z11[45];
output->z12[46] = beliefPenaltyMPC_z11[46];
output->z12[47] = beliefPenaltyMPC_z11[47];
output->z12[48] = beliefPenaltyMPC_z11[48];
output->z12[49] = beliefPenaltyMPC_z11[49];
output->z12[50] = beliefPenaltyMPC_z11[50];
output->z12[51] = beliefPenaltyMPC_z11[51];
output->z12[52] = beliefPenaltyMPC_z11[52];
output->z12[53] = beliefPenaltyMPC_z11[53];
output->z12[54] = beliefPenaltyMPC_z11[54];
output->z12[55] = beliefPenaltyMPC_z11[55];
output->z12[56] = beliefPenaltyMPC_z11[56];
output->z12[57] = beliefPenaltyMPC_z11[57];
output->z12[58] = beliefPenaltyMPC_z11[58];
output->z12[59] = beliefPenaltyMPC_z11[59];
output->z12[60] = beliefPenaltyMPC_z11[60];
output->z12[61] = beliefPenaltyMPC_z11[61];
output->z12[62] = beliefPenaltyMPC_z11[62];
output->z12[63] = beliefPenaltyMPC_z11[63];
output->z12[64] = beliefPenaltyMPC_z11[64];
output->z12[65] = beliefPenaltyMPC_z11[65];
output->z12[66] = beliefPenaltyMPC_z11[66];
output->z12[67] = beliefPenaltyMPC_z11[67];
output->z12[68] = beliefPenaltyMPC_z11[68];
output->z12[69] = beliefPenaltyMPC_z11[69];
output->z12[70] = beliefPenaltyMPC_z11[70];
output->z12[71] = beliefPenaltyMPC_z11[71];
output->z12[72] = beliefPenaltyMPC_z11[72];
output->z12[73] = beliefPenaltyMPC_z11[73];
output->z12[74] = beliefPenaltyMPC_z11[74];
output->z12[75] = beliefPenaltyMPC_z11[75];
output->z12[76] = beliefPenaltyMPC_z11[76];
output->z12[77] = beliefPenaltyMPC_z11[77];
output->z12[78] = beliefPenaltyMPC_z11[78];
output->z13[0] = beliefPenaltyMPC_z12[0];
output->z13[1] = beliefPenaltyMPC_z12[1];
output->z13[2] = beliefPenaltyMPC_z12[2];
output->z13[3] = beliefPenaltyMPC_z12[3];
output->z13[4] = beliefPenaltyMPC_z12[4];
output->z13[5] = beliefPenaltyMPC_z12[5];
output->z13[6] = beliefPenaltyMPC_z12[6];
output->z13[7] = beliefPenaltyMPC_z12[7];
output->z13[8] = beliefPenaltyMPC_z12[8];
output->z13[9] = beliefPenaltyMPC_z12[9];
output->z13[10] = beliefPenaltyMPC_z12[10];
output->z13[11] = beliefPenaltyMPC_z12[11];
output->z13[12] = beliefPenaltyMPC_z12[12];
output->z13[13] = beliefPenaltyMPC_z12[13];
output->z13[14] = beliefPenaltyMPC_z12[14];
output->z13[15] = beliefPenaltyMPC_z12[15];
output->z13[16] = beliefPenaltyMPC_z12[16];
output->z13[17] = beliefPenaltyMPC_z12[17];
output->z13[18] = beliefPenaltyMPC_z12[18];
output->z13[19] = beliefPenaltyMPC_z12[19];
output->z13[20] = beliefPenaltyMPC_z12[20];
output->z13[21] = beliefPenaltyMPC_z12[21];
output->z13[22] = beliefPenaltyMPC_z12[22];
output->z13[23] = beliefPenaltyMPC_z12[23];
output->z13[24] = beliefPenaltyMPC_z12[24];
output->z13[25] = beliefPenaltyMPC_z12[25];
output->z13[26] = beliefPenaltyMPC_z12[26];
output->z13[27] = beliefPenaltyMPC_z12[27];
output->z13[28] = beliefPenaltyMPC_z12[28];
output->z13[29] = beliefPenaltyMPC_z12[29];
output->z13[30] = beliefPenaltyMPC_z12[30];
output->z13[31] = beliefPenaltyMPC_z12[31];
output->z13[32] = beliefPenaltyMPC_z12[32];
output->z13[33] = beliefPenaltyMPC_z12[33];
output->z13[34] = beliefPenaltyMPC_z12[34];
output->z13[35] = beliefPenaltyMPC_z12[35];
output->z13[36] = beliefPenaltyMPC_z12[36];
output->z13[37] = beliefPenaltyMPC_z12[37];
output->z13[38] = beliefPenaltyMPC_z12[38];
output->z13[39] = beliefPenaltyMPC_z12[39];
output->z13[40] = beliefPenaltyMPC_z12[40];
output->z13[41] = beliefPenaltyMPC_z12[41];
output->z13[42] = beliefPenaltyMPC_z12[42];
output->z13[43] = beliefPenaltyMPC_z12[43];
output->z13[44] = beliefPenaltyMPC_z12[44];
output->z13[45] = beliefPenaltyMPC_z12[45];
output->z13[46] = beliefPenaltyMPC_z12[46];
output->z13[47] = beliefPenaltyMPC_z12[47];
output->z13[48] = beliefPenaltyMPC_z12[48];
output->z13[49] = beliefPenaltyMPC_z12[49];
output->z13[50] = beliefPenaltyMPC_z12[50];
output->z13[51] = beliefPenaltyMPC_z12[51];
output->z13[52] = beliefPenaltyMPC_z12[52];
output->z13[53] = beliefPenaltyMPC_z12[53];
output->z13[54] = beliefPenaltyMPC_z12[54];
output->z13[55] = beliefPenaltyMPC_z12[55];
output->z13[56] = beliefPenaltyMPC_z12[56];
output->z13[57] = beliefPenaltyMPC_z12[57];
output->z13[58] = beliefPenaltyMPC_z12[58];
output->z13[59] = beliefPenaltyMPC_z12[59];
output->z13[60] = beliefPenaltyMPC_z12[60];
output->z13[61] = beliefPenaltyMPC_z12[61];
output->z13[62] = beliefPenaltyMPC_z12[62];
output->z13[63] = beliefPenaltyMPC_z12[63];
output->z13[64] = beliefPenaltyMPC_z12[64];
output->z13[65] = beliefPenaltyMPC_z12[65];
output->z13[66] = beliefPenaltyMPC_z12[66];
output->z13[67] = beliefPenaltyMPC_z12[67];
output->z13[68] = beliefPenaltyMPC_z12[68];
output->z13[69] = beliefPenaltyMPC_z12[69];
output->z13[70] = beliefPenaltyMPC_z12[70];
output->z13[71] = beliefPenaltyMPC_z12[71];
output->z13[72] = beliefPenaltyMPC_z12[72];
output->z13[73] = beliefPenaltyMPC_z12[73];
output->z13[74] = beliefPenaltyMPC_z12[74];
output->z13[75] = beliefPenaltyMPC_z12[75];
output->z13[76] = beliefPenaltyMPC_z12[76];
output->z13[77] = beliefPenaltyMPC_z12[77];
output->z13[78] = beliefPenaltyMPC_z12[78];
output->z14[0] = beliefPenaltyMPC_z13[0];
output->z14[1] = beliefPenaltyMPC_z13[1];
output->z14[2] = beliefPenaltyMPC_z13[2];
output->z14[3] = beliefPenaltyMPC_z13[3];
output->z14[4] = beliefPenaltyMPC_z13[4];
output->z14[5] = beliefPenaltyMPC_z13[5];
output->z14[6] = beliefPenaltyMPC_z13[6];
output->z14[7] = beliefPenaltyMPC_z13[7];
output->z14[8] = beliefPenaltyMPC_z13[8];
output->z14[9] = beliefPenaltyMPC_z13[9];
output->z14[10] = beliefPenaltyMPC_z13[10];
output->z14[11] = beliefPenaltyMPC_z13[11];
output->z14[12] = beliefPenaltyMPC_z13[12];
output->z14[13] = beliefPenaltyMPC_z13[13];
output->z14[14] = beliefPenaltyMPC_z13[14];
output->z14[15] = beliefPenaltyMPC_z13[15];
output->z14[16] = beliefPenaltyMPC_z13[16];
output->z14[17] = beliefPenaltyMPC_z13[17];
output->z14[18] = beliefPenaltyMPC_z13[18];
output->z14[19] = beliefPenaltyMPC_z13[19];
output->z14[20] = beliefPenaltyMPC_z13[20];
output->z14[21] = beliefPenaltyMPC_z13[21];
output->z14[22] = beliefPenaltyMPC_z13[22];
output->z14[23] = beliefPenaltyMPC_z13[23];
output->z14[24] = beliefPenaltyMPC_z13[24];
output->z14[25] = beliefPenaltyMPC_z13[25];
output->z14[26] = beliefPenaltyMPC_z13[26];
output->z14[27] = beliefPenaltyMPC_z13[27];
output->z14[28] = beliefPenaltyMPC_z13[28];
output->z14[29] = beliefPenaltyMPC_z13[29];
output->z14[30] = beliefPenaltyMPC_z13[30];
output->z14[31] = beliefPenaltyMPC_z13[31];
output->z14[32] = beliefPenaltyMPC_z13[32];
output->z14[33] = beliefPenaltyMPC_z13[33];
output->z14[34] = beliefPenaltyMPC_z13[34];
output->z14[35] = beliefPenaltyMPC_z13[35];
output->z14[36] = beliefPenaltyMPC_z13[36];
output->z14[37] = beliefPenaltyMPC_z13[37];
output->z14[38] = beliefPenaltyMPC_z13[38];
output->z14[39] = beliefPenaltyMPC_z13[39];
output->z14[40] = beliefPenaltyMPC_z13[40];
output->z14[41] = beliefPenaltyMPC_z13[41];
output->z14[42] = beliefPenaltyMPC_z13[42];
output->z14[43] = beliefPenaltyMPC_z13[43];
output->z14[44] = beliefPenaltyMPC_z13[44];
output->z14[45] = beliefPenaltyMPC_z13[45];
output->z14[46] = beliefPenaltyMPC_z13[46];
output->z14[47] = beliefPenaltyMPC_z13[47];
output->z14[48] = beliefPenaltyMPC_z13[48];
output->z14[49] = beliefPenaltyMPC_z13[49];
output->z14[50] = beliefPenaltyMPC_z13[50];
output->z14[51] = beliefPenaltyMPC_z13[51];
output->z14[52] = beliefPenaltyMPC_z13[52];
output->z14[53] = beliefPenaltyMPC_z13[53];
output->z14[54] = beliefPenaltyMPC_z13[54];
output->z14[55] = beliefPenaltyMPC_z13[55];
output->z14[56] = beliefPenaltyMPC_z13[56];
output->z14[57] = beliefPenaltyMPC_z13[57];
output->z14[58] = beliefPenaltyMPC_z13[58];
output->z14[59] = beliefPenaltyMPC_z13[59];
output->z14[60] = beliefPenaltyMPC_z13[60];
output->z14[61] = beliefPenaltyMPC_z13[61];
output->z14[62] = beliefPenaltyMPC_z13[62];
output->z14[63] = beliefPenaltyMPC_z13[63];
output->z14[64] = beliefPenaltyMPC_z13[64];
output->z14[65] = beliefPenaltyMPC_z13[65];
output->z14[66] = beliefPenaltyMPC_z13[66];
output->z14[67] = beliefPenaltyMPC_z13[67];
output->z14[68] = beliefPenaltyMPC_z13[68];
output->z14[69] = beliefPenaltyMPC_z13[69];
output->z14[70] = beliefPenaltyMPC_z13[70];
output->z14[71] = beliefPenaltyMPC_z13[71];
output->z14[72] = beliefPenaltyMPC_z13[72];
output->z14[73] = beliefPenaltyMPC_z13[73];
output->z14[74] = beliefPenaltyMPC_z13[74];
output->z14[75] = beliefPenaltyMPC_z13[75];
output->z14[76] = beliefPenaltyMPC_z13[76];
output->z14[77] = beliefPenaltyMPC_z13[77];
output->z14[78] = beliefPenaltyMPC_z13[78];
output->z15[0] = beliefPenaltyMPC_z14[0];
output->z15[1] = beliefPenaltyMPC_z14[1];
output->z15[2] = beliefPenaltyMPC_z14[2];
output->z15[3] = beliefPenaltyMPC_z14[3];
output->z15[4] = beliefPenaltyMPC_z14[4];
output->z15[5] = beliefPenaltyMPC_z14[5];
output->z15[6] = beliefPenaltyMPC_z14[6];
output->z15[7] = beliefPenaltyMPC_z14[7];
output->z15[8] = beliefPenaltyMPC_z14[8];
output->z15[9] = beliefPenaltyMPC_z14[9];
output->z15[10] = beliefPenaltyMPC_z14[10];
output->z15[11] = beliefPenaltyMPC_z14[11];
output->z15[12] = beliefPenaltyMPC_z14[12];
output->z15[13] = beliefPenaltyMPC_z14[13];
output->z15[14] = beliefPenaltyMPC_z14[14];
output->z15[15] = beliefPenaltyMPC_z14[15];
output->z15[16] = beliefPenaltyMPC_z14[16];
output->z15[17] = beliefPenaltyMPC_z14[17];
output->z15[18] = beliefPenaltyMPC_z14[18];
output->z15[19] = beliefPenaltyMPC_z14[19];
output->z15[20] = beliefPenaltyMPC_z14[20];
output->z15[21] = beliefPenaltyMPC_z14[21];
output->z15[22] = beliefPenaltyMPC_z14[22];
output->z15[23] = beliefPenaltyMPC_z14[23];
output->z15[24] = beliefPenaltyMPC_z14[24];
output->z15[25] = beliefPenaltyMPC_z14[25];
output->z15[26] = beliefPenaltyMPC_z14[26];
output->z15[27] = beliefPenaltyMPC_z14[27];
output->z15[28] = beliefPenaltyMPC_z14[28];
output->z15[29] = beliefPenaltyMPC_z14[29];
output->z15[30] = beliefPenaltyMPC_z14[30];
output->z15[31] = beliefPenaltyMPC_z14[31];
output->z15[32] = beliefPenaltyMPC_z14[32];
output->z15[33] = beliefPenaltyMPC_z14[33];
output->z15[34] = beliefPenaltyMPC_z14[34];
output->z15[35] = beliefPenaltyMPC_z14[35];
output->z15[36] = beliefPenaltyMPC_z14[36];
output->z15[37] = beliefPenaltyMPC_z14[37];
output->z15[38] = beliefPenaltyMPC_z14[38];
output->z15[39] = beliefPenaltyMPC_z14[39];
output->z15[40] = beliefPenaltyMPC_z14[40];
output->z15[41] = beliefPenaltyMPC_z14[41];
output->z15[42] = beliefPenaltyMPC_z14[42];
output->z15[43] = beliefPenaltyMPC_z14[43];
output->z15[44] = beliefPenaltyMPC_z14[44];
output->z15[45] = beliefPenaltyMPC_z14[45];
output->z15[46] = beliefPenaltyMPC_z14[46];
output->z15[47] = beliefPenaltyMPC_z14[47];
output->z15[48] = beliefPenaltyMPC_z14[48];
output->z15[49] = beliefPenaltyMPC_z14[49];
output->z15[50] = beliefPenaltyMPC_z14[50];
output->z15[51] = beliefPenaltyMPC_z14[51];
output->z15[52] = beliefPenaltyMPC_z14[52];
output->z15[53] = beliefPenaltyMPC_z14[53];
output->z15[54] = beliefPenaltyMPC_z14[54];
output->z15[55] = beliefPenaltyMPC_z14[55];
output->z15[56] = beliefPenaltyMPC_z14[56];
output->z15[57] = beliefPenaltyMPC_z14[57];
output->z15[58] = beliefPenaltyMPC_z14[58];
output->z15[59] = beliefPenaltyMPC_z14[59];
output->z15[60] = beliefPenaltyMPC_z14[60];
output->z15[61] = beliefPenaltyMPC_z14[61];
output->z15[62] = beliefPenaltyMPC_z14[62];
output->z15[63] = beliefPenaltyMPC_z14[63];
output->z15[64] = beliefPenaltyMPC_z14[64];
output->z15[65] = beliefPenaltyMPC_z14[65];
output->z15[66] = beliefPenaltyMPC_z14[66];
output->z15[67] = beliefPenaltyMPC_z14[67];
output->z15[68] = beliefPenaltyMPC_z14[68];
output->z15[69] = beliefPenaltyMPC_z14[69];
output->z15[70] = beliefPenaltyMPC_z14[70];
output->z15[71] = beliefPenaltyMPC_z14[71];
output->z15[72] = beliefPenaltyMPC_z14[72];
output->z15[73] = beliefPenaltyMPC_z14[73];
output->z15[74] = beliefPenaltyMPC_z14[74];
output->z15[75] = beliefPenaltyMPC_z14[75];
output->z15[76] = beliefPenaltyMPC_z14[76];

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
