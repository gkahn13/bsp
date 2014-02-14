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

#include "trajMPC.h"

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
 * Initializes a vector of length 102 with a value.
 */
void trajMPC_LA_INITIALIZEVECTOR_102(trajMPC_FLOAT* vec, trajMPC_FLOAT value)
{
	int i;
	for( i=0; i<102; i++ )
	{
		vec[i] = value;
	}
}


/*
 * Initializes a vector of length 30 with a value.
 */
void trajMPC_LA_INITIALIZEVECTOR_30(trajMPC_FLOAT* vec, trajMPC_FLOAT value)
{
	int i;
	for( i=0; i<30; i++ )
	{
		vec[i] = value;
	}
}


/*
 * Initializes a vector of length 150 with a value.
 */
void trajMPC_LA_INITIALIZEVECTOR_150(trajMPC_FLOAT* vec, trajMPC_FLOAT value)
{
	int i;
	for( i=0; i<150; i++ )
	{
		vec[i] = value;
	}
}


/* 
 * Calculates a dot product and adds it to a variable: z += x'*y; 
 * This function is for vectors of length 150.
 */
void trajMPC_LA_DOTACC_150(trajMPC_FLOAT *x, trajMPC_FLOAT *y, trajMPC_FLOAT *z)
{
	int i;
	for( i=0; i<150; i++ ){
		*z += x[i]*y[i];
	}
}


/*
 * Calculates the gradient and the value for a quadratic function 0.5*z'*H*z + f'*z
 *
 * INPUTS:     H  - Symmetric Hessian, diag matrix of size [11 x 11]
 *             f  - column vector of size 11
 *             z  - column vector of size 11
 *
 * OUTPUTS: grad  - gradient at z (= H*z + f), column vector of size 11
 *          value <-- value + 0.5*z'*H*z + f'*z (value will be modified)
 */
void trajMPC_LA_DIAG_QUADFCN_11(trajMPC_FLOAT* H, trajMPC_FLOAT* f, trajMPC_FLOAT* z, trajMPC_FLOAT* grad, trajMPC_FLOAT* value)
{
	int i;
	trajMPC_FLOAT hz;	
	for( i=0; i<11; i++){
		hz = H[i]*z[i];
		grad[i] = hz + f[i];
		*value += 0.5*hz*z[i] + f[i]*z[i];
	}
}


/*
 * Calculates the gradient and the value for a quadratic function 0.5*z'*H*z + f'*z
 *
 * INPUTS:     H  - Symmetric Hessian, diag matrix of size [3 x 3]
 *             f  - column vector of size 3
 *             z  - column vector of size 3
 *
 * OUTPUTS: grad  - gradient at z (= H*z + f), column vector of size 3
 *          value <-- value + 0.5*z'*H*z + f'*z (value will be modified)
 */
void trajMPC_LA_DIAG_QUADFCN_3(trajMPC_FLOAT* H, trajMPC_FLOAT* f, trajMPC_FLOAT* z, trajMPC_FLOAT* grad, trajMPC_FLOAT* value)
{
	int i;
	trajMPC_FLOAT hz;	
	for( i=0; i<3; i++){
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
void trajMPC_LA_DIAGZERO_MVMSUB6_3(trajMPC_FLOAT *B, trajMPC_FLOAT *u, trajMPC_FLOAT *b, trajMPC_FLOAT *l, trajMPC_FLOAT *r, trajMPC_FLOAT *z, trajMPC_FLOAT *y)
{
	int i;
	trajMPC_FLOAT Bu[3];
	trajMPC_FLOAT norm = *y;
	trajMPC_FLOAT lr = 0;

	/* do A*x + B*u first */
	for( i=0; i<3; i++ ){
		Bu[i] = B[i]*u[i];
	}	

	for( i=0; i<3; i++ ){
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
void trajMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(trajMPC_FLOAT *A, trajMPC_FLOAT *x, trajMPC_FLOAT *B, trajMPC_FLOAT *u, trajMPC_FLOAT *b, trajMPC_FLOAT *l, trajMPC_FLOAT *r, trajMPC_FLOAT *z, trajMPC_FLOAT *y)
{
	int i;
	int j;
	int k = 0;
	trajMPC_FLOAT AxBu[3];
	trajMPC_FLOAT norm = *y;
	trajMPC_FLOAT lr = 0;

	/* do A*x + B*u first */
	for( i=0; i<3; i++ ){
		AxBu[i] = A[k++]*x[0] + B[i]*u[i];
	}	

	for( j=1; j<11; j++ ){		
		for( i=0; i<3; i++ ){
			AxBu[i] += A[k++]*x[j];
		}
	}

	for( i=0; i<3; i++ ){
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
void trajMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_3(trajMPC_FLOAT *A, trajMPC_FLOAT *x, trajMPC_FLOAT *B, trajMPC_FLOAT *u, trajMPC_FLOAT *b, trajMPC_FLOAT *l, trajMPC_FLOAT *r, trajMPC_FLOAT *z, trajMPC_FLOAT *y)
{
	int i;
	int j;
	int k = 0;
	trajMPC_FLOAT AxBu[3];
	trajMPC_FLOAT norm = *y;
	trajMPC_FLOAT lr = 0;

	/* do A*x + B*u first */
	for( i=0; i<3; i++ ){
		AxBu[i] = A[k++]*x[0] + B[i]*u[i];
	}	

	for( j=1; j<11; j++ ){		
		for( i=0; i<3; i++ ){
			AxBu[i] += A[k++]*x[j];
		}
	}

	for( i=0; i<3; i++ ){
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
 * where A is of size [3 x 11] and stored in column major format.
 * and B is of size [3 x 11] and stored in diagzero format
 * Note the transposes of A and B!
 */
void trajMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(trajMPC_FLOAT *A, trajMPC_FLOAT *x, trajMPC_FLOAT *B, trajMPC_FLOAT *y, trajMPC_FLOAT *z)
{
	int i;
	int j;
	int k = 0;
	for( i=0; i<3; i++ ){
		z[i] = 0;
		for( j=0; j<3; j++ ){
			z[i] += A[k++]*x[j];
		}
		z[i] += B[i]*y[i];
	}
	for( i=3 ;i<11; i++ ){
		z[i] = 0;
		for( j=0; j<3; j++ ){
			z[i] += A[k++]*x[j];
		}
	}
}


/*
 * Matrix vector multiplication y = M'*x where M is of size [3 x 3]
 * and stored in diagzero format. Note the transpose of M!
 */
void trajMPC_LA_DIAGZERO_MTVM_3_3(trajMPC_FLOAT *M, trajMPC_FLOAT *x, trajMPC_FLOAT *y)
{
	int i;
	for( i=0; i<3; i++ ){
		y[i] = M[i]*x[i];
	}
}


/*
 * Vector subtraction and addition.
 *	 Input: five vectors t, tidx, u, v, w and two scalars z and r
 *	 Output: y = t(tidx) - u + w
 *           z = z - v'*x;
 *           r = max([norm(y,inf), z]);
 * for vectors of length 11. Output z is of course scalar.
 */
void trajMPC_LA_VSUBADD3_11(trajMPC_FLOAT* t, trajMPC_FLOAT* u, int* uidx, trajMPC_FLOAT* v, trajMPC_FLOAT* w, trajMPC_FLOAT* y, trajMPC_FLOAT* z, trajMPC_FLOAT* r)
{
	int i;
	trajMPC_FLOAT norm = *r;
	trajMPC_FLOAT vx = 0;
	trajMPC_FLOAT x;
	for( i=0; i<11; i++){
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
void trajMPC_LA_VSUBADD2_5(trajMPC_FLOAT* t, int* tidx, trajMPC_FLOAT* u, trajMPC_FLOAT* v, trajMPC_FLOAT* w, trajMPC_FLOAT* y, trajMPC_FLOAT* z, trajMPC_FLOAT* r)
{
	int i;
	trajMPC_FLOAT norm = *r;
	trajMPC_FLOAT vx = 0;
	trajMPC_FLOAT x;
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
 * Vector subtraction and addition.
 *	 Input: five vectors t, tidx, u, v, w and two scalars z and r
 *	 Output: y = t(tidx) - u + w
 *           z = z - v'*x;
 *           r = max([norm(y,inf), z]);
 * for vectors of length 3. Output z is of course scalar.
 */
void trajMPC_LA_VSUBADD3_3(trajMPC_FLOAT* t, trajMPC_FLOAT* u, int* uidx, trajMPC_FLOAT* v, trajMPC_FLOAT* w, trajMPC_FLOAT* y, trajMPC_FLOAT* z, trajMPC_FLOAT* r)
{
	int i;
	trajMPC_FLOAT norm = *r;
	trajMPC_FLOAT vx = 0;
	trajMPC_FLOAT x;
	for( i=0; i<3; i++){
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
 * for vectors of length 3. Output z is of course scalar.
 */
void trajMPC_LA_VSUBADD2_3(trajMPC_FLOAT* t, int* tidx, trajMPC_FLOAT* u, trajMPC_FLOAT* v, trajMPC_FLOAT* w, trajMPC_FLOAT* y, trajMPC_FLOAT* z, trajMPC_FLOAT* r)
{
	int i;
	trajMPC_FLOAT norm = *r;
	trajMPC_FLOAT vx = 0;
	trajMPC_FLOAT x;
	for( i=0; i<3; i++){
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
 * Special function for box constraints of length 11
 * Returns also L/S, a value that is often used elsewhere.
 */
void trajMPC_LA_INEQ_B_GRAD_11_11_5(trajMPC_FLOAT *lu, trajMPC_FLOAT *su, trajMPC_FLOAT *ru, trajMPC_FLOAT *ll, trajMPC_FLOAT *sl, trajMPC_FLOAT *rl, int* lbIdx, int* ubIdx, trajMPC_FLOAT *grad, trajMPC_FLOAT *lubysu, trajMPC_FLOAT *llbysl)
{
	int i;
	for( i=0; i<11; i++ ){
		grad[i] = 0;
	}
	for( i=0; i<11; i++ ){		
		llbysl[i] = ll[i] / sl[i];
		grad[lbIdx[i]] -= llbysl[i]*rl[i];
	}
	for( i=0; i<5; i++ ){
		lubysu[i] = lu[i] / su[i];
		grad[ubIdx[i]] += lubysu[i]*ru[i];
	}
}


/*
 * Computes inequality constraints gradient-
 * Special function for box constraints of length 3
 * Returns also L/S, a value that is often used elsewhere.
 */
void trajMPC_LA_INEQ_B_GRAD_3_3_3(trajMPC_FLOAT *lu, trajMPC_FLOAT *su, trajMPC_FLOAT *ru, trajMPC_FLOAT *ll, trajMPC_FLOAT *sl, trajMPC_FLOAT *rl, int* lbIdx, int* ubIdx, trajMPC_FLOAT *grad, trajMPC_FLOAT *lubysu, trajMPC_FLOAT *llbysl)
{
	int i;
	for( i=0; i<3; i++ ){
		grad[i] = 0;
	}
	for( i=0; i<3; i++ ){		
		llbysl[i] = ll[i] / sl[i];
		grad[lbIdx[i]] -= llbysl[i]*rl[i];
	}
	for( i=0; i<3; i++ ){
		lubysu[i] = lu[i] / su[i];
		grad[ubIdx[i]] += lubysu[i]*ru[i];
	}
}


/*
 * Addition of three vectors  z = u + w + v
 * of length 102.
 */
void trajMPC_LA_VVADD3_102(trajMPC_FLOAT *u, trajMPC_FLOAT *v, trajMPC_FLOAT *w, trajMPC_FLOAT *z)
{
	int i;
	for( i=0; i<102; i++ ){
		z[i] = u[i] + v[i] + w[i];
	}
}


/*
 * Special function to compute the diagonal cholesky factorization of the 
 * positive definite augmented Hessian for block size 11.
 *
 * Inputs: - H = diagonal cost Hessian in diagonal storage format
 *         - llbysl = L / S of lower bounds
 *         - lubysu = L / S of upper bounds
 *
 * Output: Phi = sqrt(H + diag(llbysl) + diag(lubysu))
 * where Phi is stored in diagonal storage format
 */
void trajMPC_LA_DIAG_CHOL_LBUB_11_11_5(trajMPC_FLOAT *H, trajMPC_FLOAT *llbysl, int* lbIdx, trajMPC_FLOAT *lubysu, int* ubIdx, trajMPC_FLOAT *Phi)


{
	int i;
	
	/* copy  H into PHI */
	for( i=0; i<11; i++ ){
		Phi[i] = H[i];
	}

	/* add llbysl onto Phi where necessary */
	for( i=0; i<11; i++ ){
		Phi[lbIdx[i]] += llbysl[i];
	}

	/* add lubysu onto Phi where necessary */
	for( i=0; i<5; i++){
		Phi[ubIdx[i]] +=  lubysu[i];
	}
	
	/* compute cholesky */
	for(i=0; i<11; i++)
	{
#if trajMPC_SET_PRINTLEVEL > 0 && defined PRINTNUMERICALWARNINGS
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
 * where A is to be computed and is of size [3 x 11],
 * B is given and of size [3 x 11], L is a diagonal
 * matrix of size 3 stored in diagonal matrix 
 * storage format. Note the transpose of L has no impact!
 *
 * Result: A in column major storage format.
 *
 */
void trajMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(trajMPC_FLOAT *L, trajMPC_FLOAT *B, trajMPC_FLOAT *A)
{
    int i,j;
	 int k = 0;

	for( j=0; j<11; j++){
		for( i=0; i<3; i++){
			A[k] = B[k]/L[j];
			k++;
		}
	}

}


/**
 * Forward substitution for the matrix equation A*L' = B
 * where A is to be computed and is of size [3 x 11],
 * B is given and of size [3 x 11], L is a diagonal
 *  matrix of size 11 stored in diagonal 
 * storage format. Note the transpose of L!
 *
 * Result: A in diagonalzero storage format.
 *
 */
void trajMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(trajMPC_FLOAT *L, trajMPC_FLOAT *B, trajMPC_FLOAT *A)
{
	int j;
    for( j=0; j<11; j++ ){   
		A[j] = B[j]/L[j];
     }
}


/**
 * Compute C = A*B' where 
 *
 *	size(A) = [3 x 11]
 *  size(B) = [3 x 11] in diagzero format
 * 
 * A and C matrices are stored in column major format.
 * 
 * 
 */
void trajMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(trajMPC_FLOAT *A, trajMPC_FLOAT *B, trajMPC_FLOAT *C)
{
    int i, j;
	
	for( i=0; i<3; i++ ){
		for( j=0; j<3; j++){
			C[j*3+i] = B[i*3+j]*A[i];
		}
	}

}


/**
 * Forward substitution to solve L*y = b where L is a
 * diagonal matrix in vector storage format.
 * 
 * The dimensions involved are 11.
 */
void trajMPC_LA_DIAG_FORWARDSUB_11(trajMPC_FLOAT *L, trajMPC_FLOAT *b, trajMPC_FLOAT *y)
{
    int i;

    for( i=0; i<11; i++ ){
		y[i] = b[i]/L[i];
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
void trajMPC_LA_DIAG_CHOL_LBUB_5_11_5(trajMPC_FLOAT *H, trajMPC_FLOAT *llbysl, int* lbIdx, trajMPC_FLOAT *lubysu, int* ubIdx, trajMPC_FLOAT *Phi)


{
	int i;
	
	/* copy  H into PHI */
	for( i=0; i<5; i++ ){
		Phi[i] = H[i];
	}

	/* add llbysl onto Phi where necessary */
	for( i=0; i<11; i++ ){
		Phi[lbIdx[i]] += llbysl[i];
	}

	/* add lubysu onto Phi where necessary */
	for( i=0; i<5; i++){
		Phi[ubIdx[i]] +=  lubysu[i];
	}
	
	/* compute cholesky */
	for(i=0; i<5; i++)
	{
#if trajMPC_SET_PRINTLEVEL > 0 && defined PRINTNUMERICALWARNINGS
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


/*
 * Special function to compute the diagonal cholesky factorization of the 
 * positive definite augmented Hessian for block size 3.
 *
 * Inputs: - H = diagonal cost Hessian in diagonal storage format
 *         - llbysl = L / S of lower bounds
 *         - lubysu = L / S of upper bounds
 *
 * Output: Phi = sqrt(H + diag(llbysl) + diag(lubysu))
 * where Phi is stored in diagonal storage format
 */
void trajMPC_LA_DIAG_CHOL_ONELOOP_LBUB_3_3_3(trajMPC_FLOAT *H, trajMPC_FLOAT *llbysl, int* lbIdx, trajMPC_FLOAT *lubysu, int* ubIdx, trajMPC_FLOAT *Phi)


{
	int i;
	
	/* compute cholesky */
	for( i=0; i<3; i++ ){
		Phi[i] = H[i] + llbysl[i] + lubysu[i];

#if trajMPC_SET_PRINTLEVEL > 0 && defined PRINTNUMERICALWARNINGS
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
 * where A is to be computed and is of size [3 x 3],
 * B is given and of size [3 x 3], L is a diagonal
 *  matrix of size 3 stored in diagonal 
 * storage format. Note the transpose of L!
 *
 * Result: A in diagonalzero storage format.
 *
 */
void trajMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_3(trajMPC_FLOAT *L, trajMPC_FLOAT *B, trajMPC_FLOAT *A)
{
	int j;
    for( j=0; j<3; j++ ){   
		A[j] = B[j]/L[j];
     }
}


/**
 * Forward substitution to solve L*y = b where L is a
 * diagonal matrix in vector storage format.
 * 
 * The dimensions involved are 3.
 */
void trajMPC_LA_DIAG_FORWARDSUB_3(trajMPC_FLOAT *L, trajMPC_FLOAT *b, trajMPC_FLOAT *y)
{
    int i;

    for( i=0; i<3; i++ ){
		y[i] = b[i]/L[i];
    }
}


/**
 * Compute L = B*B', where L is lower triangular of size NXp1
 * and B is a diagzero matrix of size [3 x 11] in column
 * storage format.
 * 
 */
void trajMPC_LA_DIAGZERO_MMT_3(trajMPC_FLOAT *B, trajMPC_FLOAT *L)
{
    int i, ii, di;
    
    ii = 0; di = 0;
    for( i=0; i<3; i++ ){        
		L[ii+i] = B[i]*B[i];
        ii += ++di;
    }
}


/* 
 * Computes r = b - B*u
 * B is stored in diagzero format
 */
void trajMPC_LA_DIAGZERO_MVMSUB7_3(trajMPC_FLOAT *B, trajMPC_FLOAT *u, trajMPC_FLOAT *b, trajMPC_FLOAT *r)
{
	int i;

	for( i=0; i<3; i++ ){
		r[i] = b[i] - B[i]*u[i];
	}	
	
}


/**
 * Compute L = A*A' + B*B', where L is lower triangular of size NXp1
 * and A is a dense matrix of size [3 x 11] in column
 * storage format, and B is of size [3 x 11] diagonalzero
 * storage format.
 * 
 * THIS ONE HAS THE WORST ACCES PATTERN POSSIBLE. 
 * POSSIBKE FIX: PUT A AND B INTO ROW MAJOR FORMAT FIRST.
 * 
 */
void trajMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(trajMPC_FLOAT *A, trajMPC_FLOAT *B, trajMPC_FLOAT *L)
{
    int i, j, k, ii, di;
    trajMPC_FLOAT ltemp;
    
    ii = 0; di = 0;
    for( i=0; i<3; i++ ){        
        for( j=0; j<=i; j++ ){
            ltemp = 0; 
            for( k=0; k<11; k++ ){
                ltemp += A[k*3+i]*A[k*3+j];
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
void trajMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(trajMPC_FLOAT *A, trajMPC_FLOAT *x, trajMPC_FLOAT *B, trajMPC_FLOAT *u, trajMPC_FLOAT *b, trajMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;

	for( i=0; i<3; i++ ){
		r[i] = b[i] - A[k++]*x[0] - B[i]*u[i];
	}	

	for( j=1; j<11; j++ ){		
		for( i=0; i<3; i++ ){
			r[i] -= A[k++]*x[j];
		}
	}
	
}


/**
 * Compute L = A*A' + B*B', where L is lower triangular of size NXp1
 * and A is a dense matrix of size [3 x 11] in column
 * storage format, and B is of size [3 x 3] diagonalzero
 * storage format.
 * 
 * THIS ONE HAS THE WORST ACCES PATTERN POSSIBLE. 
 * POSSIBKE FIX: PUT A AND B INTO ROW MAJOR FORMAT FIRST.
 * 
 */
void trajMPC_LA_DENSE_DIAGZERO_MMT2_3_11_3(trajMPC_FLOAT *A, trajMPC_FLOAT *B, trajMPC_FLOAT *L)
{
    int i, j, k, ii, di;
    trajMPC_FLOAT ltemp;
    
    ii = 0; di = 0;
    for( i=0; i<3; i++ ){        
        for( j=0; j<=i; j++ ){
            ltemp = 0; 
            for( k=0; k<11; k++ ){
                ltemp += A[k*3+i]*A[k*3+j];
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
void trajMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_3(trajMPC_FLOAT *A, trajMPC_FLOAT *x, trajMPC_FLOAT *B, trajMPC_FLOAT *u, trajMPC_FLOAT *b, trajMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;

	for( i=0; i<3; i++ ){
		r[i] = b[i] - A[k++]*x[0] - B[i]*u[i];
	}	

	for( j=1; j<11; j++ ){		
		for( i=0; i<3; i++ ){
			r[i] -= A[k++]*x[j];
		}
	}
	
}


/**
 * Cholesky factorization as above, but working on a matrix in 
 * lower triangular storage format of size 3 and outputting
 * the Cholesky factor to matrix L in lower triangular format.
 */
void trajMPC_LA_DENSE_CHOL_3(trajMPC_FLOAT *A, trajMPC_FLOAT *L)
{
    int i, j, k, di, dj;
	 int ii, jj;

    trajMPC_FLOAT l;
    trajMPC_FLOAT Mii;

	/* copy A to L first and then operate on L */
	/* COULD BE OPTIMIZED */
	ii=0; di=0;
	for( i=0; i<3; i++ ){
		for( j=0; j<=i; j++ ){
			L[ii+j] = A[ii+j];
		}
		ii += ++di;
	}    
	
	/* factor L */
	ii=0; di=0;
    for( i=0; i<3; i++ ){
        l = 0;
        for( k=0; k<i; k++ ){
            l += L[ii+k]*L[ii+k];
        }        
        
        Mii = L[ii+i] - l;
        
#if trajMPC_SET_PRINTLEVEL > 0 && defined PRINTNUMERICALWARNINGS
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
        for( j=i+1; j<3; j++ ){
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
 * The dimensions involved are 3.
 */
void trajMPC_LA_DENSE_FORWARDSUB_3(trajMPC_FLOAT *L, trajMPC_FLOAT *b, trajMPC_FLOAT *y)
{
    int i,j,ii,di;
    trajMPC_FLOAT yel;
            
    ii = 0; di = 0;
    for( i=0; i<3; i++ ){
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
 * where A is to be computed and is of size [3 x 3],
 * B is given and of size [3 x 3], L is a lower tri-
 * angular matrix of size 3 stored in lower triangular 
 * storage format. Note the transpose of L AND B!
 *
 * Result: A in column major storage format.
 *
 */
void trajMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(trajMPC_FLOAT *L, trajMPC_FLOAT *B, trajMPC_FLOAT *A)
{
    int i,j,k,ii,di;
    trajMPC_FLOAT a;
    
    ii=0; di=0;
    for( j=0; j<3; j++ ){        
        for( i=0; i<3; i++ ){
            a = B[i*3+j];
            for( k=0; k<j; k++ ){
                a -= A[k*3+i]*L[ii+k];
            }    

			/* saturate for numerical stability */
			a = MIN(a, BIGM);
			a = MAX(a, -BIGM); 

			A[j*3+i] = a/L[ii+j];			
        }
        ii += ++di;
    }
}


/**
 * Compute L = L - A*A', where L is lower triangular of size 3
 * and A is a dense matrix of size [3 x 3] in column
 * storage format.
 * 
 * THIS ONE HAS THE WORST ACCES PATTERN POSSIBLE. 
 * POSSIBKE FIX: PUT A INTO ROW MAJOR FORMAT FIRST.
 * 
 */
void trajMPC_LA_DENSE_MMTSUB_3_3(trajMPC_FLOAT *A, trajMPC_FLOAT *L)
{
    int i, j, k, ii, di;
    trajMPC_FLOAT ltemp;
    
    ii = 0; di = 0;
    for( i=0; i<3; i++ ){        
        for( j=0; j<=i; j++ ){
            ltemp = 0; 
            for( k=0; k<3; k++ ){
                ltemp += A[k*3+i]*A[k*3+j];
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
void trajMPC_LA_DENSE_MVMSUB1_3_3(trajMPC_FLOAT *A, trajMPC_FLOAT *x, trajMPC_FLOAT *b, trajMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;

	for( i=0; i<3; i++ ){
		r[i] = b[i] - A[k++]*x[0];
	}	
	for( j=1; j<3; j++ ){		
		for( i=0; i<3; i++ ){
			r[i] -= A[k++]*x[j];
		}
	}
}


/**
 * Backward Substitution to solve L^T*x = y where L is a
 * lower triangular matrix in triangular storage format.
 * 
 * All involved dimensions are 3.
 */
void trajMPC_LA_DENSE_BACKWARDSUB_3(trajMPC_FLOAT *L, trajMPC_FLOAT *y, trajMPC_FLOAT *x)
{
    int i, ii, di, j, jj, dj;
    trajMPC_FLOAT xel;    
	int start = 3;
    
    /* now solve L^T*x = y by backward substitution */
    ii = start; di = 2;
    for( i=2; i>=0; i-- ){        
        xel = y[i];        
        jj = start; dj = 2;
        for( j=2; j>i; j-- ){
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
 * Matrix vector multiplication y = b - M'*x where M is of size [3 x 3]
 * and stored in column major format. Note the transpose of M!
 */
void trajMPC_LA_DENSE_MTVMSUB_3_3(trajMPC_FLOAT *A, trajMPC_FLOAT *x, trajMPC_FLOAT *b, trajMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0; 
	for( i=0; i<3; i++ ){
		r[i] = b[i];
		for( j=0; j<3; j++ ){
			r[i] -= A[k++]*x[j];
		}
	}
}


/*
 * Vector subtraction z = -x - y for vectors of length 102.
 */
void trajMPC_LA_VSUB2_102(trajMPC_FLOAT *x, trajMPC_FLOAT *y, trajMPC_FLOAT *z)
{
	int i;
	for( i=0; i<102; i++){
		z[i] = -x[i] - y[i];
	}
}


/**
 * Forward-Backward-Substitution to solve L*L^T*x = b where L is a
 * diagonal matrix of size 11 in vector
 * storage format.
 */
void trajMPC_LA_DIAG_FORWARDBACKWARDSUB_11(trajMPC_FLOAT *L, trajMPC_FLOAT *b, trajMPC_FLOAT *x)
{
    int i;
            
    /* solve Ly = b by forward and backward substitution */
    for( i=0; i<11; i++ ){
		x[i] = b[i]/(L[i]*L[i]);
    }
    
}


/**
 * Forward-Backward-Substitution to solve L*L^T*x = b where L is a
 * diagonal matrix of size 3 in vector
 * storage format.
 */
void trajMPC_LA_DIAG_FORWARDBACKWARDSUB_3(trajMPC_FLOAT *L, trajMPC_FLOAT *b, trajMPC_FLOAT *x)
{
    int i;
            
    /* solve Ly = b by forward and backward substitution */
    for( i=0; i<3; i++ ){
		x[i] = b[i]/(L[i]*L[i]);
    }
    
}


/*
 * Vector subtraction z = x(xidx) - y where y, z and xidx are of length 11,
 * and x has length 11 and is indexed through yidx.
 */
void trajMPC_LA_VSUB_INDEXED_11(trajMPC_FLOAT *x, int* xidx, trajMPC_FLOAT *y, trajMPC_FLOAT *z)
{
	int i;
	for( i=0; i<11; i++){
		z[i] = x[xidx[i]] - y[i];
	}
}


/*
 * Vector subtraction x = -u.*v - w for vectors of length 11.
 */
void trajMPC_LA_VSUB3_11(trajMPC_FLOAT *u, trajMPC_FLOAT *v, trajMPC_FLOAT *w, trajMPC_FLOAT *x)
{
	int i;
	for( i=0; i<11; i++){
		x[i] = -u[i]*v[i] - w[i];
	}
}


/*
 * Vector subtraction z = -x - y(yidx) where y is of length 11
 * and z, x and yidx are of length 5.
 */
void trajMPC_LA_VSUB2_INDEXED_5(trajMPC_FLOAT *x, trajMPC_FLOAT *y, int* yidx, trajMPC_FLOAT *z)
{
	int i;
	for( i=0; i<5; i++){
		z[i] = -x[i] - y[yidx[i]];
	}
}


/*
 * Vector subtraction x = -u.*v - w for vectors of length 5.
 */
void trajMPC_LA_VSUB3_5(trajMPC_FLOAT *u, trajMPC_FLOAT *v, trajMPC_FLOAT *w, trajMPC_FLOAT *x)
{
	int i;
	for( i=0; i<5; i++){
		x[i] = -u[i]*v[i] - w[i];
	}
}


/*
 * Vector subtraction z = x(xidx) - y where y, z and xidx are of length 3,
 * and x has length 3 and is indexed through yidx.
 */
void trajMPC_LA_VSUB_INDEXED_3(trajMPC_FLOAT *x, int* xidx, trajMPC_FLOAT *y, trajMPC_FLOAT *z)
{
	int i;
	for( i=0; i<3; i++){
		z[i] = x[xidx[i]] - y[i];
	}
}


/*
 * Vector subtraction x = -u.*v - w for vectors of length 3.
 */
void trajMPC_LA_VSUB3_3(trajMPC_FLOAT *u, trajMPC_FLOAT *v, trajMPC_FLOAT *w, trajMPC_FLOAT *x)
{
	int i;
	for( i=0; i<3; i++){
		x[i] = -u[i]*v[i] - w[i];
	}
}


/*
 * Vector subtraction z = -x - y(yidx) where y is of length 3
 * and z, x and yidx are of length 3.
 */
void trajMPC_LA_VSUB2_INDEXED_3(trajMPC_FLOAT *x, trajMPC_FLOAT *y, int* yidx, trajMPC_FLOAT *z)
{
	int i;
	for( i=0; i<3; i++){
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
 * trajMPC_NOPROGRESS (should be negative).
 */
int trajMPC_LINESEARCH_BACKTRACKING_AFFINE(trajMPC_FLOAT *l, trajMPC_FLOAT *s, trajMPC_FLOAT *dl, trajMPC_FLOAT *ds, trajMPC_FLOAT *a, trajMPC_FLOAT *mu_aff)
{
    int i;
	int lsIt=1;    
    trajMPC_FLOAT dltemp;
    trajMPC_FLOAT dstemp;
    trajMPC_FLOAT mya = 1.0;
    trajMPC_FLOAT mymu;
        
    while( 1 ){                        

        /* 
         * Compute both snew and wnew together.
         * We compute also mu_affine along the way here, as the
         * values might be in registers, so it should be cheaper.
         */
        mymu = 0;
        for( i=0; i<150; i++ ){
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
        if( i == 150 ){
            break;
        } else {
            mya *= trajMPC_SET_LS_SCALE_AFF;
            if( mya < trajMPC_SET_LS_MINSTEP ){
                return trajMPC_NOPROGRESS;
            }
        }
    }
    
    /* return new values and iteration counter */
    *a = mya;
    *mu_aff = mymu / (trajMPC_FLOAT)150;
    return lsIt;
}


/*
 * Vector subtraction x = u.*v - a where a is a scalar
*  and x,u,v are vectors of length 150.
 */
void trajMPC_LA_VSUB5_150(trajMPC_FLOAT *u, trajMPC_FLOAT *v, trajMPC_FLOAT a, trajMPC_FLOAT *x)
{
	int i;
	for( i=0; i<150; i++){
		x[i] = u[i]*v[i] - a;
	}
}


/*
 * Computes x=0; x(uidx) += u/su; x(vidx) -= v/sv where x is of length 11,
 * u, su, uidx are of length 5 and v, sv, vidx are of length 11.
 */
void trajMPC_LA_VSUB6_INDEXED_11_5_11(trajMPC_FLOAT *u, trajMPC_FLOAT *su, int* uidx, trajMPC_FLOAT *v, trajMPC_FLOAT *sv, int* vidx, trajMPC_FLOAT *x)
{
	int i;
	for( i=0; i<11; i++ ){
		x[i] = 0;
	}
	for( i=0; i<5; i++){
		x[uidx[i]] += u[i]/su[i];
	}
	for( i=0; i<11; i++){
		x[vidx[i]] -= v[i]/sv[i];
	}
}


/* 
 * Computes r =  B*u
 * where B is stored in diagzero format
 */
void trajMPC_LA_DIAGZERO_MVM_3(trajMPC_FLOAT *B, trajMPC_FLOAT *u, trajMPC_FLOAT *r)
{
	int i;

	for( i=0; i<3; i++ ){
		r[i] = B[i]*u[i];
	}	
	
}


/* 
 * Computes r = A*x + B*u
 * where A is stored in column major format
 * and B is stored in diagzero format
 */
void trajMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(trajMPC_FLOAT *A, trajMPC_FLOAT *x, trajMPC_FLOAT *B, trajMPC_FLOAT *u, trajMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;

	for( i=0; i<3; i++ ){
		r[i] = A[k++]*x[0] + B[i]*u[i];
	}	

	for( j=1; j<11; j++ ){		
		for( i=0; i<3; i++ ){
			r[i] += A[k++]*x[j];
		}
	}
	
}


/*
 * Computes x=0; x(uidx) += u/su; x(vidx) -= v/sv where x is of length 3,
 * u, su, uidx are of length 3 and v, sv, vidx are of length 3.
 */
void trajMPC_LA_VSUB6_INDEXED_3_3_3(trajMPC_FLOAT *u, trajMPC_FLOAT *su, int* uidx, trajMPC_FLOAT *v, trajMPC_FLOAT *sv, int* vidx, trajMPC_FLOAT *x)
{
	int i;
	for( i=0; i<3; i++ ){
		x[i] = 0;
	}
	for( i=0; i<3; i++){
		x[uidx[i]] += u[i]/su[i];
	}
	for( i=0; i<3; i++){
		x[vidx[i]] -= v[i]/sv[i];
	}
}


/* 
 * Computes r = A*x + B*u
 * where A is stored in column major format
 * and B is stored in diagzero format
 */
void trajMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_3(trajMPC_FLOAT *A, trajMPC_FLOAT *x, trajMPC_FLOAT *B, trajMPC_FLOAT *u, trajMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;

	for( i=0; i<3; i++ ){
		r[i] = A[k++]*x[0] + B[i]*u[i];
	}	

	for( j=1; j<11; j++ ){		
		for( i=0; i<3; i++ ){
			r[i] += A[k++]*x[j];
		}
	}
	
}


/*
 * Vector subtraction z = x - y for vectors of length 102.
 */
void trajMPC_LA_VSUB_102(trajMPC_FLOAT *x, trajMPC_FLOAT *y, trajMPC_FLOAT *z)
{
	int i;
	for( i=0; i<102; i++){
		z[i] = x[i] - y[i];
	}
}


/** 
 * Computes z = -r./s - u.*y(y)
 * where all vectors except of y are of length 11 (length of y >= 11).
 */
void trajMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(trajMPC_FLOAT *r, trajMPC_FLOAT *s, trajMPC_FLOAT *u, trajMPC_FLOAT *y, int* yidx, trajMPC_FLOAT *z)
{
	int i;
	for( i=0; i<11; i++ ){
		z[i] = -r[i]/s[i] - u[i]*y[yidx[i]];
	}
}


/** 
 * Computes z = -r./s + u.*y(y)
 * where all vectors except of y are of length 5 (length of y >= 5).
 */
void trajMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(trajMPC_FLOAT *r, trajMPC_FLOAT *s, trajMPC_FLOAT *u, trajMPC_FLOAT *y, int* yidx, trajMPC_FLOAT *z)
{
	int i;
	for( i=0; i<5; i++ ){
		z[i] = -r[i]/s[i] + u[i]*y[yidx[i]];
	}
}


/** 
 * Computes z = -r./s - u.*y(y)
 * where all vectors except of y are of length 3 (length of y >= 3).
 */
void trajMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_3(trajMPC_FLOAT *r, trajMPC_FLOAT *s, trajMPC_FLOAT *u, trajMPC_FLOAT *y, int* yidx, trajMPC_FLOAT *z)
{
	int i;
	for( i=0; i<3; i++ ){
		z[i] = -r[i]/s[i] - u[i]*y[yidx[i]];
	}
}


/** 
 * Computes z = -r./s + u.*y(y)
 * where all vectors except of y are of length 3 (length of y >= 3).
 */
void trajMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_3(trajMPC_FLOAT *r, trajMPC_FLOAT *s, trajMPC_FLOAT *u, trajMPC_FLOAT *y, int* yidx, trajMPC_FLOAT *z)
{
	int i;
	for( i=0; i<3; i++ ){
		z[i] = -r[i]/s[i] + u[i]*y[yidx[i]];
	}
}


/*
 * Computes ds = -l.\(r + s.*dl) for vectors of length 150.
 */
void trajMPC_LA_VSUB7_150(trajMPC_FLOAT *l, trajMPC_FLOAT *r, trajMPC_FLOAT *s, trajMPC_FLOAT *dl, trajMPC_FLOAT *ds)
{
	int i;
	for( i=0; i<150; i++){
		ds[i] = -(r[i] + s[i]*dl[i])/l[i];
	}
}


/*
 * Vector addition x = x + y for vectors of length 102.
 */
void trajMPC_LA_VADD_102(trajMPC_FLOAT *x, trajMPC_FLOAT *y)
{
	int i;
	for( i=0; i<102; i++){
		x[i] += y[i];
	}
}


/*
 * Vector addition x = x + y for vectors of length 30.
 */
void trajMPC_LA_VADD_30(trajMPC_FLOAT *x, trajMPC_FLOAT *y)
{
	int i;
	for( i=0; i<30; i++){
		x[i] += y[i];
	}
}


/*
 * Vector addition x = x + y for vectors of length 150.
 */
void trajMPC_LA_VADD_150(trajMPC_FLOAT *x, trajMPC_FLOAT *y)
{
	int i;
	for( i=0; i<150; i++){
		x[i] += y[i];
	}
}


/**
 * Backtracking line search for combined predictor/corrector step.
 * Update on variables with safety factor gamma (to keep us away from
 * boundary).
 */
int trajMPC_LINESEARCH_BACKTRACKING_COMBINED(trajMPC_FLOAT *z, trajMPC_FLOAT *v, trajMPC_FLOAT *l, trajMPC_FLOAT *s, trajMPC_FLOAT *dz, trajMPC_FLOAT *dv, trajMPC_FLOAT *dl, trajMPC_FLOAT *ds, trajMPC_FLOAT *a, trajMPC_FLOAT *mu)
{
    int i, lsIt=1;       
    trajMPC_FLOAT dltemp;
    trajMPC_FLOAT dstemp;    
    trajMPC_FLOAT a_gamma;
            
    *a = 1.0;
    while( 1 ){                        

        /* check whether search criterion is fulfilled */
        for( i=0; i<150; i++ ){
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
        if( i == 150 ){
            break;
        } else {
            *a *= trajMPC_SET_LS_SCALE;
            if( *a < trajMPC_SET_LS_MINSTEP ){
                return trajMPC_NOPROGRESS;
            }
        }
    }
    
    /* update variables with safety margin */
    a_gamma = (*a)*trajMPC_SET_LS_MAXSTEP;
    
    /* primal variables */
    for( i=0; i<102; i++ ){
        z[i] += a_gamma*dz[i];
    }
    
    /* equality constraint multipliers */
    for( i=0; i<30; i++ ){
        v[i] += a_gamma*dv[i];
    }
    
    /* inequality constraint multipliers & slacks, also update mu */
    *mu = 0;
    for( i=0; i<150; i++ ){
        dltemp = l[i] + a_gamma*dl[i]; l[i] = dltemp;
        dstemp = s[i] + a_gamma*ds[i]; s[i] = dstemp;
        *mu += dltemp*dstemp;
    }
    
    *a = a_gamma;
    *mu /= (trajMPC_FLOAT)150;
    return lsIt;
}




/* VARIABLE DEFINITIONS ------------------------------------------------ */
trajMPC_FLOAT trajMPC_z[102];
trajMPC_FLOAT trajMPC_v[30];
trajMPC_FLOAT trajMPC_dz_aff[102];
trajMPC_FLOAT trajMPC_dv_aff[30];
trajMPC_FLOAT trajMPC_grad_cost[102];
trajMPC_FLOAT trajMPC_grad_eq[102];
trajMPC_FLOAT trajMPC_rd[102];
trajMPC_FLOAT trajMPC_l[150];
trajMPC_FLOAT trajMPC_s[150];
trajMPC_FLOAT trajMPC_lbys[150];
trajMPC_FLOAT trajMPC_dl_aff[150];
trajMPC_FLOAT trajMPC_ds_aff[150];
trajMPC_FLOAT trajMPC_dz_cc[102];
trajMPC_FLOAT trajMPC_dv_cc[30];
trajMPC_FLOAT trajMPC_dl_cc[150];
trajMPC_FLOAT trajMPC_ds_cc[150];
trajMPC_FLOAT trajMPC_ccrhs[150];
trajMPC_FLOAT trajMPC_grad_ineq[102];
trajMPC_FLOAT trajMPC_H0[11] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 2.0000000000000000E+000, 2.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
trajMPC_FLOAT* trajMPC_z0 = trajMPC_z + 0;
trajMPC_FLOAT* trajMPC_dzaff0 = trajMPC_dz_aff + 0;
trajMPC_FLOAT* trajMPC_dzcc0 = trajMPC_dz_cc + 0;
trajMPC_FLOAT* trajMPC_rd0 = trajMPC_rd + 0;
trajMPC_FLOAT trajMPC_Lbyrd0[11];
trajMPC_FLOAT* trajMPC_grad_cost0 = trajMPC_grad_cost + 0;
trajMPC_FLOAT* trajMPC_grad_eq0 = trajMPC_grad_eq + 0;
trajMPC_FLOAT* trajMPC_grad_ineq0 = trajMPC_grad_ineq + 0;
trajMPC_FLOAT trajMPC_ctv0[11];
trajMPC_FLOAT* trajMPC_v0 = trajMPC_v + 0;
trajMPC_FLOAT trajMPC_re0[3];
trajMPC_FLOAT trajMPC_beta0[3];
trajMPC_FLOAT trajMPC_betacc0[3];
trajMPC_FLOAT* trajMPC_dvaff0 = trajMPC_dv_aff + 0;
trajMPC_FLOAT* trajMPC_dvcc0 = trajMPC_dv_cc + 0;
trajMPC_FLOAT trajMPC_V0[33];
trajMPC_FLOAT trajMPC_Yd0[6];
trajMPC_FLOAT trajMPC_Ld0[6];
trajMPC_FLOAT trajMPC_yy0[3];
trajMPC_FLOAT trajMPC_bmy0[3];
int trajMPC_lbIdx0[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
trajMPC_FLOAT* trajMPC_llb0 = trajMPC_l + 0;
trajMPC_FLOAT* trajMPC_slb0 = trajMPC_s + 0;
trajMPC_FLOAT* trajMPC_llbbyslb0 = trajMPC_lbys + 0;
trajMPC_FLOAT trajMPC_rilb0[11];
trajMPC_FLOAT* trajMPC_dllbaff0 = trajMPC_dl_aff + 0;
trajMPC_FLOAT* trajMPC_dslbaff0 = trajMPC_ds_aff + 0;
trajMPC_FLOAT* trajMPC_dllbcc0 = trajMPC_dl_cc + 0;
trajMPC_FLOAT* trajMPC_dslbcc0 = trajMPC_ds_cc + 0;
trajMPC_FLOAT* trajMPC_ccrhsl0 = trajMPC_ccrhs + 0;
int trajMPC_ubIdx0[5] = {0, 1, 2, 3, 4};
trajMPC_FLOAT* trajMPC_lub0 = trajMPC_l + 11;
trajMPC_FLOAT* trajMPC_sub0 = trajMPC_s + 11;
trajMPC_FLOAT* trajMPC_lubbysub0 = trajMPC_lbys + 11;
trajMPC_FLOAT trajMPC_riub0[5];
trajMPC_FLOAT* trajMPC_dlubaff0 = trajMPC_dl_aff + 11;
trajMPC_FLOAT* trajMPC_dsubaff0 = trajMPC_ds_aff + 11;
trajMPC_FLOAT* trajMPC_dlubcc0 = trajMPC_dl_cc + 11;
trajMPC_FLOAT* trajMPC_dsubcc0 = trajMPC_ds_cc + 11;
trajMPC_FLOAT* trajMPC_ccrhsub0 = trajMPC_ccrhs + 11;
trajMPC_FLOAT trajMPC_Phi0[11];
trajMPC_FLOAT trajMPC_D0[11] = {1.0000000000000000E+000, 
1.0000000000000000E+000, 
1.0000000000000000E+000};
trajMPC_FLOAT trajMPC_W0[11];
trajMPC_FLOAT trajMPC_H1[5] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 2.0000000000000000E+000, 2.0000000000000000E+000};
trajMPC_FLOAT* trajMPC_z1 = trajMPC_z + 11;
trajMPC_FLOAT* trajMPC_dzaff1 = trajMPC_dz_aff + 11;
trajMPC_FLOAT* trajMPC_dzcc1 = trajMPC_dz_cc + 11;
trajMPC_FLOAT* trajMPC_rd1 = trajMPC_rd + 11;
trajMPC_FLOAT trajMPC_Lbyrd1[11];
trajMPC_FLOAT* trajMPC_grad_cost1 = trajMPC_grad_cost + 11;
trajMPC_FLOAT* trajMPC_grad_eq1 = trajMPC_grad_eq + 11;
trajMPC_FLOAT* trajMPC_grad_ineq1 = trajMPC_grad_ineq + 11;
trajMPC_FLOAT trajMPC_ctv1[11];
trajMPC_FLOAT* trajMPC_v1 = trajMPC_v + 3;
trajMPC_FLOAT trajMPC_re1[3];
trajMPC_FLOAT trajMPC_beta1[3];
trajMPC_FLOAT trajMPC_betacc1[3];
trajMPC_FLOAT* trajMPC_dvaff1 = trajMPC_dv_aff + 3;
trajMPC_FLOAT* trajMPC_dvcc1 = trajMPC_dv_cc + 3;
trajMPC_FLOAT trajMPC_V1[33];
trajMPC_FLOAT trajMPC_Yd1[6];
trajMPC_FLOAT trajMPC_Ld1[6];
trajMPC_FLOAT trajMPC_yy1[3];
trajMPC_FLOAT trajMPC_bmy1[3];
int trajMPC_lbIdx1[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
trajMPC_FLOAT* trajMPC_llb1 = trajMPC_l + 16;
trajMPC_FLOAT* trajMPC_slb1 = trajMPC_s + 16;
trajMPC_FLOAT* trajMPC_llbbyslb1 = trajMPC_lbys + 16;
trajMPC_FLOAT trajMPC_rilb1[11];
trajMPC_FLOAT* trajMPC_dllbaff1 = trajMPC_dl_aff + 16;
trajMPC_FLOAT* trajMPC_dslbaff1 = trajMPC_ds_aff + 16;
trajMPC_FLOAT* trajMPC_dllbcc1 = trajMPC_dl_cc + 16;
trajMPC_FLOAT* trajMPC_dslbcc1 = trajMPC_ds_cc + 16;
trajMPC_FLOAT* trajMPC_ccrhsl1 = trajMPC_ccrhs + 16;
int trajMPC_ubIdx1[5] = {0, 1, 2, 3, 4};
trajMPC_FLOAT* trajMPC_lub1 = trajMPC_l + 27;
trajMPC_FLOAT* trajMPC_sub1 = trajMPC_s + 27;
trajMPC_FLOAT* trajMPC_lubbysub1 = trajMPC_lbys + 27;
trajMPC_FLOAT trajMPC_riub1[5];
trajMPC_FLOAT* trajMPC_dlubaff1 = trajMPC_dl_aff + 27;
trajMPC_FLOAT* trajMPC_dsubaff1 = trajMPC_ds_aff + 27;
trajMPC_FLOAT* trajMPC_dlubcc1 = trajMPC_dl_cc + 27;
trajMPC_FLOAT* trajMPC_dsubcc1 = trajMPC_ds_cc + 27;
trajMPC_FLOAT* trajMPC_ccrhsub1 = trajMPC_ccrhs + 27;
trajMPC_FLOAT trajMPC_Phi1[11];
trajMPC_FLOAT trajMPC_D1[11] = {-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000};
trajMPC_FLOAT trajMPC_W1[11];
trajMPC_FLOAT trajMPC_Ysd1[9];
trajMPC_FLOAT trajMPC_Lsd1[9];
trajMPC_FLOAT* trajMPC_z2 = trajMPC_z + 22;
trajMPC_FLOAT* trajMPC_dzaff2 = trajMPC_dz_aff + 22;
trajMPC_FLOAT* trajMPC_dzcc2 = trajMPC_dz_cc + 22;
trajMPC_FLOAT* trajMPC_rd2 = trajMPC_rd + 22;
trajMPC_FLOAT trajMPC_Lbyrd2[11];
trajMPC_FLOAT* trajMPC_grad_cost2 = trajMPC_grad_cost + 22;
trajMPC_FLOAT* trajMPC_grad_eq2 = trajMPC_grad_eq + 22;
trajMPC_FLOAT* trajMPC_grad_ineq2 = trajMPC_grad_ineq + 22;
trajMPC_FLOAT trajMPC_ctv2[11];
trajMPC_FLOAT* trajMPC_v2 = trajMPC_v + 6;
trajMPC_FLOAT trajMPC_re2[3];
trajMPC_FLOAT trajMPC_beta2[3];
trajMPC_FLOAT trajMPC_betacc2[3];
trajMPC_FLOAT* trajMPC_dvaff2 = trajMPC_dv_aff + 6;
trajMPC_FLOAT* trajMPC_dvcc2 = trajMPC_dv_cc + 6;
trajMPC_FLOAT trajMPC_V2[33];
trajMPC_FLOAT trajMPC_Yd2[6];
trajMPC_FLOAT trajMPC_Ld2[6];
trajMPC_FLOAT trajMPC_yy2[3];
trajMPC_FLOAT trajMPC_bmy2[3];
int trajMPC_lbIdx2[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
trajMPC_FLOAT* trajMPC_llb2 = trajMPC_l + 32;
trajMPC_FLOAT* trajMPC_slb2 = trajMPC_s + 32;
trajMPC_FLOAT* trajMPC_llbbyslb2 = trajMPC_lbys + 32;
trajMPC_FLOAT trajMPC_rilb2[11];
trajMPC_FLOAT* trajMPC_dllbaff2 = trajMPC_dl_aff + 32;
trajMPC_FLOAT* trajMPC_dslbaff2 = trajMPC_ds_aff + 32;
trajMPC_FLOAT* trajMPC_dllbcc2 = trajMPC_dl_cc + 32;
trajMPC_FLOAT* trajMPC_dslbcc2 = trajMPC_ds_cc + 32;
trajMPC_FLOAT* trajMPC_ccrhsl2 = trajMPC_ccrhs + 32;
int trajMPC_ubIdx2[5] = {0, 1, 2, 3, 4};
trajMPC_FLOAT* trajMPC_lub2 = trajMPC_l + 43;
trajMPC_FLOAT* trajMPC_sub2 = trajMPC_s + 43;
trajMPC_FLOAT* trajMPC_lubbysub2 = trajMPC_lbys + 43;
trajMPC_FLOAT trajMPC_riub2[5];
trajMPC_FLOAT* trajMPC_dlubaff2 = trajMPC_dl_aff + 43;
trajMPC_FLOAT* trajMPC_dsubaff2 = trajMPC_ds_aff + 43;
trajMPC_FLOAT* trajMPC_dlubcc2 = trajMPC_dl_cc + 43;
trajMPC_FLOAT* trajMPC_dsubcc2 = trajMPC_ds_cc + 43;
trajMPC_FLOAT* trajMPC_ccrhsub2 = trajMPC_ccrhs + 43;
trajMPC_FLOAT trajMPC_Phi2[11];
trajMPC_FLOAT trajMPC_W2[11];
trajMPC_FLOAT trajMPC_Ysd2[9];
trajMPC_FLOAT trajMPC_Lsd2[9];
trajMPC_FLOAT* trajMPC_z3 = trajMPC_z + 33;
trajMPC_FLOAT* trajMPC_dzaff3 = trajMPC_dz_aff + 33;
trajMPC_FLOAT* trajMPC_dzcc3 = trajMPC_dz_cc + 33;
trajMPC_FLOAT* trajMPC_rd3 = trajMPC_rd + 33;
trajMPC_FLOAT trajMPC_Lbyrd3[11];
trajMPC_FLOAT* trajMPC_grad_cost3 = trajMPC_grad_cost + 33;
trajMPC_FLOAT* trajMPC_grad_eq3 = trajMPC_grad_eq + 33;
trajMPC_FLOAT* trajMPC_grad_ineq3 = trajMPC_grad_ineq + 33;
trajMPC_FLOAT trajMPC_ctv3[11];
trajMPC_FLOAT* trajMPC_v3 = trajMPC_v + 9;
trajMPC_FLOAT trajMPC_re3[3];
trajMPC_FLOAT trajMPC_beta3[3];
trajMPC_FLOAT trajMPC_betacc3[3];
trajMPC_FLOAT* trajMPC_dvaff3 = trajMPC_dv_aff + 9;
trajMPC_FLOAT* trajMPC_dvcc3 = trajMPC_dv_cc + 9;
trajMPC_FLOAT trajMPC_V3[33];
trajMPC_FLOAT trajMPC_Yd3[6];
trajMPC_FLOAT trajMPC_Ld3[6];
trajMPC_FLOAT trajMPC_yy3[3];
trajMPC_FLOAT trajMPC_bmy3[3];
int trajMPC_lbIdx3[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
trajMPC_FLOAT* trajMPC_llb3 = trajMPC_l + 48;
trajMPC_FLOAT* trajMPC_slb3 = trajMPC_s + 48;
trajMPC_FLOAT* trajMPC_llbbyslb3 = trajMPC_lbys + 48;
trajMPC_FLOAT trajMPC_rilb3[11];
trajMPC_FLOAT* trajMPC_dllbaff3 = trajMPC_dl_aff + 48;
trajMPC_FLOAT* trajMPC_dslbaff3 = trajMPC_ds_aff + 48;
trajMPC_FLOAT* trajMPC_dllbcc3 = trajMPC_dl_cc + 48;
trajMPC_FLOAT* trajMPC_dslbcc3 = trajMPC_ds_cc + 48;
trajMPC_FLOAT* trajMPC_ccrhsl3 = trajMPC_ccrhs + 48;
int trajMPC_ubIdx3[5] = {0, 1, 2, 3, 4};
trajMPC_FLOAT* trajMPC_lub3 = trajMPC_l + 59;
trajMPC_FLOAT* trajMPC_sub3 = trajMPC_s + 59;
trajMPC_FLOAT* trajMPC_lubbysub3 = trajMPC_lbys + 59;
trajMPC_FLOAT trajMPC_riub3[5];
trajMPC_FLOAT* trajMPC_dlubaff3 = trajMPC_dl_aff + 59;
trajMPC_FLOAT* trajMPC_dsubaff3 = trajMPC_ds_aff + 59;
trajMPC_FLOAT* trajMPC_dlubcc3 = trajMPC_dl_cc + 59;
trajMPC_FLOAT* trajMPC_dsubcc3 = trajMPC_ds_cc + 59;
trajMPC_FLOAT* trajMPC_ccrhsub3 = trajMPC_ccrhs + 59;
trajMPC_FLOAT trajMPC_Phi3[11];
trajMPC_FLOAT trajMPC_W3[11];
trajMPC_FLOAT trajMPC_Ysd3[9];
trajMPC_FLOAT trajMPC_Lsd3[9];
trajMPC_FLOAT* trajMPC_z4 = trajMPC_z + 44;
trajMPC_FLOAT* trajMPC_dzaff4 = trajMPC_dz_aff + 44;
trajMPC_FLOAT* trajMPC_dzcc4 = trajMPC_dz_cc + 44;
trajMPC_FLOAT* trajMPC_rd4 = trajMPC_rd + 44;
trajMPC_FLOAT trajMPC_Lbyrd4[11];
trajMPC_FLOAT* trajMPC_grad_cost4 = trajMPC_grad_cost + 44;
trajMPC_FLOAT* trajMPC_grad_eq4 = trajMPC_grad_eq + 44;
trajMPC_FLOAT* trajMPC_grad_ineq4 = trajMPC_grad_ineq + 44;
trajMPC_FLOAT trajMPC_ctv4[11];
trajMPC_FLOAT* trajMPC_v4 = trajMPC_v + 12;
trajMPC_FLOAT trajMPC_re4[3];
trajMPC_FLOAT trajMPC_beta4[3];
trajMPC_FLOAT trajMPC_betacc4[3];
trajMPC_FLOAT* trajMPC_dvaff4 = trajMPC_dv_aff + 12;
trajMPC_FLOAT* trajMPC_dvcc4 = trajMPC_dv_cc + 12;
trajMPC_FLOAT trajMPC_V4[33];
trajMPC_FLOAT trajMPC_Yd4[6];
trajMPC_FLOAT trajMPC_Ld4[6];
trajMPC_FLOAT trajMPC_yy4[3];
trajMPC_FLOAT trajMPC_bmy4[3];
int trajMPC_lbIdx4[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
trajMPC_FLOAT* trajMPC_llb4 = trajMPC_l + 64;
trajMPC_FLOAT* trajMPC_slb4 = trajMPC_s + 64;
trajMPC_FLOAT* trajMPC_llbbyslb4 = trajMPC_lbys + 64;
trajMPC_FLOAT trajMPC_rilb4[11];
trajMPC_FLOAT* trajMPC_dllbaff4 = trajMPC_dl_aff + 64;
trajMPC_FLOAT* trajMPC_dslbaff4 = trajMPC_ds_aff + 64;
trajMPC_FLOAT* trajMPC_dllbcc4 = trajMPC_dl_cc + 64;
trajMPC_FLOAT* trajMPC_dslbcc4 = trajMPC_ds_cc + 64;
trajMPC_FLOAT* trajMPC_ccrhsl4 = trajMPC_ccrhs + 64;
int trajMPC_ubIdx4[5] = {0, 1, 2, 3, 4};
trajMPC_FLOAT* trajMPC_lub4 = trajMPC_l + 75;
trajMPC_FLOAT* trajMPC_sub4 = trajMPC_s + 75;
trajMPC_FLOAT* trajMPC_lubbysub4 = trajMPC_lbys + 75;
trajMPC_FLOAT trajMPC_riub4[5];
trajMPC_FLOAT* trajMPC_dlubaff4 = trajMPC_dl_aff + 75;
trajMPC_FLOAT* trajMPC_dsubaff4 = trajMPC_ds_aff + 75;
trajMPC_FLOAT* trajMPC_dlubcc4 = trajMPC_dl_cc + 75;
trajMPC_FLOAT* trajMPC_dsubcc4 = trajMPC_ds_cc + 75;
trajMPC_FLOAT* trajMPC_ccrhsub4 = trajMPC_ccrhs + 75;
trajMPC_FLOAT trajMPC_Phi4[11];
trajMPC_FLOAT trajMPC_W4[11];
trajMPC_FLOAT trajMPC_Ysd4[9];
trajMPC_FLOAT trajMPC_Lsd4[9];
trajMPC_FLOAT* trajMPC_z5 = trajMPC_z + 55;
trajMPC_FLOAT* trajMPC_dzaff5 = trajMPC_dz_aff + 55;
trajMPC_FLOAT* trajMPC_dzcc5 = trajMPC_dz_cc + 55;
trajMPC_FLOAT* trajMPC_rd5 = trajMPC_rd + 55;
trajMPC_FLOAT trajMPC_Lbyrd5[11];
trajMPC_FLOAT* trajMPC_grad_cost5 = trajMPC_grad_cost + 55;
trajMPC_FLOAT* trajMPC_grad_eq5 = trajMPC_grad_eq + 55;
trajMPC_FLOAT* trajMPC_grad_ineq5 = trajMPC_grad_ineq + 55;
trajMPC_FLOAT trajMPC_ctv5[11];
trajMPC_FLOAT* trajMPC_v5 = trajMPC_v + 15;
trajMPC_FLOAT trajMPC_re5[3];
trajMPC_FLOAT trajMPC_beta5[3];
trajMPC_FLOAT trajMPC_betacc5[3];
trajMPC_FLOAT* trajMPC_dvaff5 = trajMPC_dv_aff + 15;
trajMPC_FLOAT* trajMPC_dvcc5 = trajMPC_dv_cc + 15;
trajMPC_FLOAT trajMPC_V5[33];
trajMPC_FLOAT trajMPC_Yd5[6];
trajMPC_FLOAT trajMPC_Ld5[6];
trajMPC_FLOAT trajMPC_yy5[3];
trajMPC_FLOAT trajMPC_bmy5[3];
int trajMPC_lbIdx5[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
trajMPC_FLOAT* trajMPC_llb5 = trajMPC_l + 80;
trajMPC_FLOAT* trajMPC_slb5 = trajMPC_s + 80;
trajMPC_FLOAT* trajMPC_llbbyslb5 = trajMPC_lbys + 80;
trajMPC_FLOAT trajMPC_rilb5[11];
trajMPC_FLOAT* trajMPC_dllbaff5 = trajMPC_dl_aff + 80;
trajMPC_FLOAT* trajMPC_dslbaff5 = trajMPC_ds_aff + 80;
trajMPC_FLOAT* trajMPC_dllbcc5 = trajMPC_dl_cc + 80;
trajMPC_FLOAT* trajMPC_dslbcc5 = trajMPC_ds_cc + 80;
trajMPC_FLOAT* trajMPC_ccrhsl5 = trajMPC_ccrhs + 80;
int trajMPC_ubIdx5[5] = {0, 1, 2, 3, 4};
trajMPC_FLOAT* trajMPC_lub5 = trajMPC_l + 91;
trajMPC_FLOAT* trajMPC_sub5 = trajMPC_s + 91;
trajMPC_FLOAT* trajMPC_lubbysub5 = trajMPC_lbys + 91;
trajMPC_FLOAT trajMPC_riub5[5];
trajMPC_FLOAT* trajMPC_dlubaff5 = trajMPC_dl_aff + 91;
trajMPC_FLOAT* trajMPC_dsubaff5 = trajMPC_ds_aff + 91;
trajMPC_FLOAT* trajMPC_dlubcc5 = trajMPC_dl_cc + 91;
trajMPC_FLOAT* trajMPC_dsubcc5 = trajMPC_ds_cc + 91;
trajMPC_FLOAT* trajMPC_ccrhsub5 = trajMPC_ccrhs + 91;
trajMPC_FLOAT trajMPC_Phi5[11];
trajMPC_FLOAT trajMPC_W5[11];
trajMPC_FLOAT trajMPC_Ysd5[9];
trajMPC_FLOAT trajMPC_Lsd5[9];
trajMPC_FLOAT* trajMPC_z6 = trajMPC_z + 66;
trajMPC_FLOAT* trajMPC_dzaff6 = trajMPC_dz_aff + 66;
trajMPC_FLOAT* trajMPC_dzcc6 = trajMPC_dz_cc + 66;
trajMPC_FLOAT* trajMPC_rd6 = trajMPC_rd + 66;
trajMPC_FLOAT trajMPC_Lbyrd6[11];
trajMPC_FLOAT* trajMPC_grad_cost6 = trajMPC_grad_cost + 66;
trajMPC_FLOAT* trajMPC_grad_eq6 = trajMPC_grad_eq + 66;
trajMPC_FLOAT* trajMPC_grad_ineq6 = trajMPC_grad_ineq + 66;
trajMPC_FLOAT trajMPC_ctv6[11];
trajMPC_FLOAT* trajMPC_v6 = trajMPC_v + 18;
trajMPC_FLOAT trajMPC_re6[3];
trajMPC_FLOAT trajMPC_beta6[3];
trajMPC_FLOAT trajMPC_betacc6[3];
trajMPC_FLOAT* trajMPC_dvaff6 = trajMPC_dv_aff + 18;
trajMPC_FLOAT* trajMPC_dvcc6 = trajMPC_dv_cc + 18;
trajMPC_FLOAT trajMPC_V6[33];
trajMPC_FLOAT trajMPC_Yd6[6];
trajMPC_FLOAT trajMPC_Ld6[6];
trajMPC_FLOAT trajMPC_yy6[3];
trajMPC_FLOAT trajMPC_bmy6[3];
int trajMPC_lbIdx6[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
trajMPC_FLOAT* trajMPC_llb6 = trajMPC_l + 96;
trajMPC_FLOAT* trajMPC_slb6 = trajMPC_s + 96;
trajMPC_FLOAT* trajMPC_llbbyslb6 = trajMPC_lbys + 96;
trajMPC_FLOAT trajMPC_rilb6[11];
trajMPC_FLOAT* trajMPC_dllbaff6 = trajMPC_dl_aff + 96;
trajMPC_FLOAT* trajMPC_dslbaff6 = trajMPC_ds_aff + 96;
trajMPC_FLOAT* trajMPC_dllbcc6 = trajMPC_dl_cc + 96;
trajMPC_FLOAT* trajMPC_dslbcc6 = trajMPC_ds_cc + 96;
trajMPC_FLOAT* trajMPC_ccrhsl6 = trajMPC_ccrhs + 96;
int trajMPC_ubIdx6[5] = {0, 1, 2, 3, 4};
trajMPC_FLOAT* trajMPC_lub6 = trajMPC_l + 107;
trajMPC_FLOAT* trajMPC_sub6 = trajMPC_s + 107;
trajMPC_FLOAT* trajMPC_lubbysub6 = trajMPC_lbys + 107;
trajMPC_FLOAT trajMPC_riub6[5];
trajMPC_FLOAT* trajMPC_dlubaff6 = trajMPC_dl_aff + 107;
trajMPC_FLOAT* trajMPC_dsubaff6 = trajMPC_ds_aff + 107;
trajMPC_FLOAT* trajMPC_dlubcc6 = trajMPC_dl_cc + 107;
trajMPC_FLOAT* trajMPC_dsubcc6 = trajMPC_ds_cc + 107;
trajMPC_FLOAT* trajMPC_ccrhsub6 = trajMPC_ccrhs + 107;
trajMPC_FLOAT trajMPC_Phi6[11];
trajMPC_FLOAT trajMPC_W6[11];
trajMPC_FLOAT trajMPC_Ysd6[9];
trajMPC_FLOAT trajMPC_Lsd6[9];
trajMPC_FLOAT* trajMPC_z7 = trajMPC_z + 77;
trajMPC_FLOAT* trajMPC_dzaff7 = trajMPC_dz_aff + 77;
trajMPC_FLOAT* trajMPC_dzcc7 = trajMPC_dz_cc + 77;
trajMPC_FLOAT* trajMPC_rd7 = trajMPC_rd + 77;
trajMPC_FLOAT trajMPC_Lbyrd7[11];
trajMPC_FLOAT* trajMPC_grad_cost7 = trajMPC_grad_cost + 77;
trajMPC_FLOAT* trajMPC_grad_eq7 = trajMPC_grad_eq + 77;
trajMPC_FLOAT* trajMPC_grad_ineq7 = trajMPC_grad_ineq + 77;
trajMPC_FLOAT trajMPC_ctv7[11];
trajMPC_FLOAT* trajMPC_v7 = trajMPC_v + 21;
trajMPC_FLOAT trajMPC_re7[3];
trajMPC_FLOAT trajMPC_beta7[3];
trajMPC_FLOAT trajMPC_betacc7[3];
trajMPC_FLOAT* trajMPC_dvaff7 = trajMPC_dv_aff + 21;
trajMPC_FLOAT* trajMPC_dvcc7 = trajMPC_dv_cc + 21;
trajMPC_FLOAT trajMPC_V7[33];
trajMPC_FLOAT trajMPC_Yd7[6];
trajMPC_FLOAT trajMPC_Ld7[6];
trajMPC_FLOAT trajMPC_yy7[3];
trajMPC_FLOAT trajMPC_bmy7[3];
int trajMPC_lbIdx7[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
trajMPC_FLOAT* trajMPC_llb7 = trajMPC_l + 112;
trajMPC_FLOAT* trajMPC_slb7 = trajMPC_s + 112;
trajMPC_FLOAT* trajMPC_llbbyslb7 = trajMPC_lbys + 112;
trajMPC_FLOAT trajMPC_rilb7[11];
trajMPC_FLOAT* trajMPC_dllbaff7 = trajMPC_dl_aff + 112;
trajMPC_FLOAT* trajMPC_dslbaff7 = trajMPC_ds_aff + 112;
trajMPC_FLOAT* trajMPC_dllbcc7 = trajMPC_dl_cc + 112;
trajMPC_FLOAT* trajMPC_dslbcc7 = trajMPC_ds_cc + 112;
trajMPC_FLOAT* trajMPC_ccrhsl7 = trajMPC_ccrhs + 112;
int trajMPC_ubIdx7[5] = {0, 1, 2, 3, 4};
trajMPC_FLOAT* trajMPC_lub7 = trajMPC_l + 123;
trajMPC_FLOAT* trajMPC_sub7 = trajMPC_s + 123;
trajMPC_FLOAT* trajMPC_lubbysub7 = trajMPC_lbys + 123;
trajMPC_FLOAT trajMPC_riub7[5];
trajMPC_FLOAT* trajMPC_dlubaff7 = trajMPC_dl_aff + 123;
trajMPC_FLOAT* trajMPC_dsubaff7 = trajMPC_ds_aff + 123;
trajMPC_FLOAT* trajMPC_dlubcc7 = trajMPC_dl_cc + 123;
trajMPC_FLOAT* trajMPC_dsubcc7 = trajMPC_ds_cc + 123;
trajMPC_FLOAT* trajMPC_ccrhsub7 = trajMPC_ccrhs + 123;
trajMPC_FLOAT trajMPC_Phi7[11];
trajMPC_FLOAT trajMPC_W7[11];
trajMPC_FLOAT trajMPC_Ysd7[9];
trajMPC_FLOAT trajMPC_Lsd7[9];
trajMPC_FLOAT* trajMPC_z8 = trajMPC_z + 88;
trajMPC_FLOAT* trajMPC_dzaff8 = trajMPC_dz_aff + 88;
trajMPC_FLOAT* trajMPC_dzcc8 = trajMPC_dz_cc + 88;
trajMPC_FLOAT* trajMPC_rd8 = trajMPC_rd + 88;
trajMPC_FLOAT trajMPC_Lbyrd8[11];
trajMPC_FLOAT* trajMPC_grad_cost8 = trajMPC_grad_cost + 88;
trajMPC_FLOAT* trajMPC_grad_eq8 = trajMPC_grad_eq + 88;
trajMPC_FLOAT* trajMPC_grad_ineq8 = trajMPC_grad_ineq + 88;
trajMPC_FLOAT trajMPC_ctv8[11];
trajMPC_FLOAT* trajMPC_v8 = trajMPC_v + 24;
trajMPC_FLOAT trajMPC_re8[3];
trajMPC_FLOAT trajMPC_beta8[3];
trajMPC_FLOAT trajMPC_betacc8[3];
trajMPC_FLOAT* trajMPC_dvaff8 = trajMPC_dv_aff + 24;
trajMPC_FLOAT* trajMPC_dvcc8 = trajMPC_dv_cc + 24;
trajMPC_FLOAT trajMPC_V8[33];
trajMPC_FLOAT trajMPC_Yd8[6];
trajMPC_FLOAT trajMPC_Ld8[6];
trajMPC_FLOAT trajMPC_yy8[3];
trajMPC_FLOAT trajMPC_bmy8[3];
int trajMPC_lbIdx8[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
trajMPC_FLOAT* trajMPC_llb8 = trajMPC_l + 128;
trajMPC_FLOAT* trajMPC_slb8 = trajMPC_s + 128;
trajMPC_FLOAT* trajMPC_llbbyslb8 = trajMPC_lbys + 128;
trajMPC_FLOAT trajMPC_rilb8[11];
trajMPC_FLOAT* trajMPC_dllbaff8 = trajMPC_dl_aff + 128;
trajMPC_FLOAT* trajMPC_dslbaff8 = trajMPC_ds_aff + 128;
trajMPC_FLOAT* trajMPC_dllbcc8 = trajMPC_dl_cc + 128;
trajMPC_FLOAT* trajMPC_dslbcc8 = trajMPC_ds_cc + 128;
trajMPC_FLOAT* trajMPC_ccrhsl8 = trajMPC_ccrhs + 128;
int trajMPC_ubIdx8[5] = {0, 1, 2, 3, 4};
trajMPC_FLOAT* trajMPC_lub8 = trajMPC_l + 139;
trajMPC_FLOAT* trajMPC_sub8 = trajMPC_s + 139;
trajMPC_FLOAT* trajMPC_lubbysub8 = trajMPC_lbys + 139;
trajMPC_FLOAT trajMPC_riub8[5];
trajMPC_FLOAT* trajMPC_dlubaff8 = trajMPC_dl_aff + 139;
trajMPC_FLOAT* trajMPC_dsubaff8 = trajMPC_ds_aff + 139;
trajMPC_FLOAT* trajMPC_dlubcc8 = trajMPC_dl_cc + 139;
trajMPC_FLOAT* trajMPC_dsubcc8 = trajMPC_ds_cc + 139;
trajMPC_FLOAT* trajMPC_ccrhsub8 = trajMPC_ccrhs + 139;
trajMPC_FLOAT trajMPC_Phi8[11];
trajMPC_FLOAT trajMPC_W8[11];
trajMPC_FLOAT trajMPC_Ysd8[9];
trajMPC_FLOAT trajMPC_Lsd8[9];
trajMPC_FLOAT trajMPC_H9[3] = {2.0000000000000000E+001, 2.0000000000000000E+001, 0.0000000000000000E+000};
trajMPC_FLOAT* trajMPC_z9 = trajMPC_z + 99;
trajMPC_FLOAT* trajMPC_dzaff9 = trajMPC_dz_aff + 99;
trajMPC_FLOAT* trajMPC_dzcc9 = trajMPC_dz_cc + 99;
trajMPC_FLOAT* trajMPC_rd9 = trajMPC_rd + 99;
trajMPC_FLOAT trajMPC_Lbyrd9[3];
trajMPC_FLOAT* trajMPC_grad_cost9 = trajMPC_grad_cost + 99;
trajMPC_FLOAT* trajMPC_grad_eq9 = trajMPC_grad_eq + 99;
trajMPC_FLOAT* trajMPC_grad_ineq9 = trajMPC_grad_ineq + 99;
trajMPC_FLOAT trajMPC_ctv9[3];
trajMPC_FLOAT* trajMPC_v9 = trajMPC_v + 27;
trajMPC_FLOAT trajMPC_re9[3];
trajMPC_FLOAT trajMPC_beta9[3];
trajMPC_FLOAT trajMPC_betacc9[3];
trajMPC_FLOAT* trajMPC_dvaff9 = trajMPC_dv_aff + 27;
trajMPC_FLOAT* trajMPC_dvcc9 = trajMPC_dv_cc + 27;
trajMPC_FLOAT trajMPC_V9[9];
trajMPC_FLOAT trajMPC_Yd9[6];
trajMPC_FLOAT trajMPC_Ld9[6];
trajMPC_FLOAT trajMPC_yy9[3];
trajMPC_FLOAT trajMPC_bmy9[3];
int trajMPC_lbIdx9[3] = {0, 1, 2};
trajMPC_FLOAT* trajMPC_llb9 = trajMPC_l + 144;
trajMPC_FLOAT* trajMPC_slb9 = trajMPC_s + 144;
trajMPC_FLOAT* trajMPC_llbbyslb9 = trajMPC_lbys + 144;
trajMPC_FLOAT trajMPC_rilb9[3];
trajMPC_FLOAT* trajMPC_dllbaff9 = trajMPC_dl_aff + 144;
trajMPC_FLOAT* trajMPC_dslbaff9 = trajMPC_ds_aff + 144;
trajMPC_FLOAT* trajMPC_dllbcc9 = trajMPC_dl_cc + 144;
trajMPC_FLOAT* trajMPC_dslbcc9 = trajMPC_ds_cc + 144;
trajMPC_FLOAT* trajMPC_ccrhsl9 = trajMPC_ccrhs + 144;
int trajMPC_ubIdx9[3] = {0, 1, 2};
trajMPC_FLOAT* trajMPC_lub9 = trajMPC_l + 147;
trajMPC_FLOAT* trajMPC_sub9 = trajMPC_s + 147;
trajMPC_FLOAT* trajMPC_lubbysub9 = trajMPC_lbys + 147;
trajMPC_FLOAT trajMPC_riub9[3];
trajMPC_FLOAT* trajMPC_dlubaff9 = trajMPC_dl_aff + 147;
trajMPC_FLOAT* trajMPC_dsubaff9 = trajMPC_ds_aff + 147;
trajMPC_FLOAT* trajMPC_dlubcc9 = trajMPC_dl_cc + 147;
trajMPC_FLOAT* trajMPC_dsubcc9 = trajMPC_ds_cc + 147;
trajMPC_FLOAT* trajMPC_ccrhsub9 = trajMPC_ccrhs + 147;
trajMPC_FLOAT trajMPC_Phi9[3];
trajMPC_FLOAT trajMPC_D9[3] = {-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000};
trajMPC_FLOAT trajMPC_W9[3];
trajMPC_FLOAT trajMPC_Ysd9[9];
trajMPC_FLOAT trajMPC_Lsd9[9];
trajMPC_FLOAT musigma;
trajMPC_FLOAT sigma_3rdroot;
trajMPC_FLOAT trajMPC_Diag1_0[11];
trajMPC_FLOAT trajMPC_Diag2_0[11];
trajMPC_FLOAT trajMPC_L_0[55];




/* SOLVER CODE --------------------------------------------------------- */
int trajMPC_solve(trajMPC_params* params, trajMPC_output* output, trajMPC_info* info)
{	
int exitcode;

#if trajMPC_SET_TIMING == 1
	trajMPC_timer solvertimer;
	trajMPC_tic(&solvertimer);
#endif
/* FUNCTION CALLS INTO LA LIBRARY -------------------------------------- */
info->it = 0;
trajMPC_LA_INITIALIZEVECTOR_102(trajMPC_z, 0);
trajMPC_LA_INITIALIZEVECTOR_30(trajMPC_v, 1);
trajMPC_LA_INITIALIZEVECTOR_150(trajMPC_l, 1);
trajMPC_LA_INITIALIZEVECTOR_150(trajMPC_s, 1);
info->mu = 0;
trajMPC_LA_DOTACC_150(trajMPC_l, trajMPC_s, &info->mu);
info->mu /= 150;
while( 1 ){
info->pobj = 0;
trajMPC_LA_DIAG_QUADFCN_11(trajMPC_H0, params->f1, trajMPC_z0, trajMPC_grad_cost0, &info->pobj);
trajMPC_LA_DIAG_QUADFCN_11(trajMPC_H1, params->f2, trajMPC_z1, trajMPC_grad_cost1, &info->pobj);
trajMPC_LA_DIAG_QUADFCN_11(trajMPC_H1, params->f3, trajMPC_z2, trajMPC_grad_cost2, &info->pobj);
trajMPC_LA_DIAG_QUADFCN_11(trajMPC_H1, params->f4, trajMPC_z3, trajMPC_grad_cost3, &info->pobj);
trajMPC_LA_DIAG_QUADFCN_11(trajMPC_H1, params->f5, trajMPC_z4, trajMPC_grad_cost4, &info->pobj);
trajMPC_LA_DIAG_QUADFCN_11(trajMPC_H1, params->f6, trajMPC_z5, trajMPC_grad_cost5, &info->pobj);
trajMPC_LA_DIAG_QUADFCN_11(trajMPC_H1, params->f7, trajMPC_z6, trajMPC_grad_cost6, &info->pobj);
trajMPC_LA_DIAG_QUADFCN_11(trajMPC_H1, params->f8, trajMPC_z7, trajMPC_grad_cost7, &info->pobj);
trajMPC_LA_DIAG_QUADFCN_11(trajMPC_H1, params->f9, trajMPC_z8, trajMPC_grad_cost8, &info->pobj);
trajMPC_LA_DIAG_QUADFCN_3(trajMPC_H9, params->f10, trajMPC_z9, trajMPC_grad_cost9, &info->pobj);
info->res_eq = 0;
info->dgap = 0;
trajMPC_LA_DIAGZERO_MVMSUB6_3(trajMPC_D0, trajMPC_z0, params->e1, trajMPC_v0, trajMPC_re0, &info->dgap, &info->res_eq);
trajMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C1, trajMPC_z0, trajMPC_D1, trajMPC_z1, params->e2, trajMPC_v1, trajMPC_re1, &info->dgap, &info->res_eq);
trajMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C2, trajMPC_z1, trajMPC_D1, trajMPC_z2, params->e3, trajMPC_v2, trajMPC_re2, &info->dgap, &info->res_eq);
trajMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C3, trajMPC_z2, trajMPC_D1, trajMPC_z3, params->e4, trajMPC_v3, trajMPC_re3, &info->dgap, &info->res_eq);
trajMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C4, trajMPC_z3, trajMPC_D1, trajMPC_z4, params->e5, trajMPC_v4, trajMPC_re4, &info->dgap, &info->res_eq);
trajMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C5, trajMPC_z4, trajMPC_D1, trajMPC_z5, params->e6, trajMPC_v5, trajMPC_re5, &info->dgap, &info->res_eq);
trajMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C6, trajMPC_z5, trajMPC_D1, trajMPC_z6, params->e7, trajMPC_v6, trajMPC_re6, &info->dgap, &info->res_eq);
trajMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C7, trajMPC_z6, trajMPC_D1, trajMPC_z7, params->e8, trajMPC_v7, trajMPC_re7, &info->dgap, &info->res_eq);
trajMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C8, trajMPC_z7, trajMPC_D1, trajMPC_z8, params->e9, trajMPC_v8, trajMPC_re8, &info->dgap, &info->res_eq);
trajMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_3(params->C9, trajMPC_z8, trajMPC_D9, trajMPC_z9, params->e10, trajMPC_v9, trajMPC_re9, &info->dgap, &info->res_eq);
trajMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C1, trajMPC_v1, trajMPC_D0, trajMPC_v0, trajMPC_grad_eq0);
trajMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C2, trajMPC_v2, trajMPC_D1, trajMPC_v1, trajMPC_grad_eq1);
trajMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C3, trajMPC_v3, trajMPC_D1, trajMPC_v2, trajMPC_grad_eq2);
trajMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C4, trajMPC_v4, trajMPC_D1, trajMPC_v3, trajMPC_grad_eq3);
trajMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C5, trajMPC_v5, trajMPC_D1, trajMPC_v4, trajMPC_grad_eq4);
trajMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C6, trajMPC_v6, trajMPC_D1, trajMPC_v5, trajMPC_grad_eq5);
trajMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C7, trajMPC_v7, trajMPC_D1, trajMPC_v6, trajMPC_grad_eq6);
trajMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C8, trajMPC_v8, trajMPC_D1, trajMPC_v7, trajMPC_grad_eq7);
trajMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C9, trajMPC_v9, trajMPC_D1, trajMPC_v8, trajMPC_grad_eq8);
trajMPC_LA_DIAGZERO_MTVM_3_3(trajMPC_D9, trajMPC_v9, trajMPC_grad_eq9);
info->res_ineq = 0;
trajMPC_LA_VSUBADD3_11(params->lb1, trajMPC_z0, trajMPC_lbIdx0, trajMPC_llb0, trajMPC_slb0, trajMPC_rilb0, &info->dgap, &info->res_ineq);
trajMPC_LA_VSUBADD2_5(trajMPC_z0, trajMPC_ubIdx0, params->ub1, trajMPC_lub0, trajMPC_sub0, trajMPC_riub0, &info->dgap, &info->res_ineq);
trajMPC_LA_VSUBADD3_11(params->lb2, trajMPC_z1, trajMPC_lbIdx1, trajMPC_llb1, trajMPC_slb1, trajMPC_rilb1, &info->dgap, &info->res_ineq);
trajMPC_LA_VSUBADD2_5(trajMPC_z1, trajMPC_ubIdx1, params->ub2, trajMPC_lub1, trajMPC_sub1, trajMPC_riub1, &info->dgap, &info->res_ineq);
trajMPC_LA_VSUBADD3_11(params->lb3, trajMPC_z2, trajMPC_lbIdx2, trajMPC_llb2, trajMPC_slb2, trajMPC_rilb2, &info->dgap, &info->res_ineq);
trajMPC_LA_VSUBADD2_5(trajMPC_z2, trajMPC_ubIdx2, params->ub3, trajMPC_lub2, trajMPC_sub2, trajMPC_riub2, &info->dgap, &info->res_ineq);
trajMPC_LA_VSUBADD3_11(params->lb4, trajMPC_z3, trajMPC_lbIdx3, trajMPC_llb3, trajMPC_slb3, trajMPC_rilb3, &info->dgap, &info->res_ineq);
trajMPC_LA_VSUBADD2_5(trajMPC_z3, trajMPC_ubIdx3, params->ub4, trajMPC_lub3, trajMPC_sub3, trajMPC_riub3, &info->dgap, &info->res_ineq);
trajMPC_LA_VSUBADD3_11(params->lb5, trajMPC_z4, trajMPC_lbIdx4, trajMPC_llb4, trajMPC_slb4, trajMPC_rilb4, &info->dgap, &info->res_ineq);
trajMPC_LA_VSUBADD2_5(trajMPC_z4, trajMPC_ubIdx4, params->ub5, trajMPC_lub4, trajMPC_sub4, trajMPC_riub4, &info->dgap, &info->res_ineq);
trajMPC_LA_VSUBADD3_11(params->lb6, trajMPC_z5, trajMPC_lbIdx5, trajMPC_llb5, trajMPC_slb5, trajMPC_rilb5, &info->dgap, &info->res_ineq);
trajMPC_LA_VSUBADD2_5(trajMPC_z5, trajMPC_ubIdx5, params->ub6, trajMPC_lub5, trajMPC_sub5, trajMPC_riub5, &info->dgap, &info->res_ineq);
trajMPC_LA_VSUBADD3_11(params->lb7, trajMPC_z6, trajMPC_lbIdx6, trajMPC_llb6, trajMPC_slb6, trajMPC_rilb6, &info->dgap, &info->res_ineq);
trajMPC_LA_VSUBADD2_5(trajMPC_z6, trajMPC_ubIdx6, params->ub7, trajMPC_lub6, trajMPC_sub6, trajMPC_riub6, &info->dgap, &info->res_ineq);
trajMPC_LA_VSUBADD3_11(params->lb8, trajMPC_z7, trajMPC_lbIdx7, trajMPC_llb7, trajMPC_slb7, trajMPC_rilb7, &info->dgap, &info->res_ineq);
trajMPC_LA_VSUBADD2_5(trajMPC_z7, trajMPC_ubIdx7, params->ub8, trajMPC_lub7, trajMPC_sub7, trajMPC_riub7, &info->dgap, &info->res_ineq);
trajMPC_LA_VSUBADD3_11(params->lb9, trajMPC_z8, trajMPC_lbIdx8, trajMPC_llb8, trajMPC_slb8, trajMPC_rilb8, &info->dgap, &info->res_ineq);
trajMPC_LA_VSUBADD2_5(trajMPC_z8, trajMPC_ubIdx8, params->ub9, trajMPC_lub8, trajMPC_sub8, trajMPC_riub8, &info->dgap, &info->res_ineq);
trajMPC_LA_VSUBADD3_3(params->lb10, trajMPC_z9, trajMPC_lbIdx9, trajMPC_llb9, trajMPC_slb9, trajMPC_rilb9, &info->dgap, &info->res_ineq);
trajMPC_LA_VSUBADD2_3(trajMPC_z9, trajMPC_ubIdx9, params->ub10, trajMPC_lub9, trajMPC_sub9, trajMPC_riub9, &info->dgap, &info->res_ineq);
trajMPC_LA_INEQ_B_GRAD_11_11_5(trajMPC_lub0, trajMPC_sub0, trajMPC_riub0, trajMPC_llb0, trajMPC_slb0, trajMPC_rilb0, trajMPC_lbIdx0, trajMPC_ubIdx0, trajMPC_grad_ineq0, trajMPC_lubbysub0, trajMPC_llbbyslb0);
trajMPC_LA_INEQ_B_GRAD_11_11_5(trajMPC_lub1, trajMPC_sub1, trajMPC_riub1, trajMPC_llb1, trajMPC_slb1, trajMPC_rilb1, trajMPC_lbIdx1, trajMPC_ubIdx1, trajMPC_grad_ineq1, trajMPC_lubbysub1, trajMPC_llbbyslb1);
trajMPC_LA_INEQ_B_GRAD_11_11_5(trajMPC_lub2, trajMPC_sub2, trajMPC_riub2, trajMPC_llb2, trajMPC_slb2, trajMPC_rilb2, trajMPC_lbIdx2, trajMPC_ubIdx2, trajMPC_grad_ineq2, trajMPC_lubbysub2, trajMPC_llbbyslb2);
trajMPC_LA_INEQ_B_GRAD_11_11_5(trajMPC_lub3, trajMPC_sub3, trajMPC_riub3, trajMPC_llb3, trajMPC_slb3, trajMPC_rilb3, trajMPC_lbIdx3, trajMPC_ubIdx3, trajMPC_grad_ineq3, trajMPC_lubbysub3, trajMPC_llbbyslb3);
trajMPC_LA_INEQ_B_GRAD_11_11_5(trajMPC_lub4, trajMPC_sub4, trajMPC_riub4, trajMPC_llb4, trajMPC_slb4, trajMPC_rilb4, trajMPC_lbIdx4, trajMPC_ubIdx4, trajMPC_grad_ineq4, trajMPC_lubbysub4, trajMPC_llbbyslb4);
trajMPC_LA_INEQ_B_GRAD_11_11_5(trajMPC_lub5, trajMPC_sub5, trajMPC_riub5, trajMPC_llb5, trajMPC_slb5, trajMPC_rilb5, trajMPC_lbIdx5, trajMPC_ubIdx5, trajMPC_grad_ineq5, trajMPC_lubbysub5, trajMPC_llbbyslb5);
trajMPC_LA_INEQ_B_GRAD_11_11_5(trajMPC_lub6, trajMPC_sub6, trajMPC_riub6, trajMPC_llb6, trajMPC_slb6, trajMPC_rilb6, trajMPC_lbIdx6, trajMPC_ubIdx6, trajMPC_grad_ineq6, trajMPC_lubbysub6, trajMPC_llbbyslb6);
trajMPC_LA_INEQ_B_GRAD_11_11_5(trajMPC_lub7, trajMPC_sub7, trajMPC_riub7, trajMPC_llb7, trajMPC_slb7, trajMPC_rilb7, trajMPC_lbIdx7, trajMPC_ubIdx7, trajMPC_grad_ineq7, trajMPC_lubbysub7, trajMPC_llbbyslb7);
trajMPC_LA_INEQ_B_GRAD_11_11_5(trajMPC_lub8, trajMPC_sub8, trajMPC_riub8, trajMPC_llb8, trajMPC_slb8, trajMPC_rilb8, trajMPC_lbIdx8, trajMPC_ubIdx8, trajMPC_grad_ineq8, trajMPC_lubbysub8, trajMPC_llbbyslb8);
trajMPC_LA_INEQ_B_GRAD_3_3_3(trajMPC_lub9, trajMPC_sub9, trajMPC_riub9, trajMPC_llb9, trajMPC_slb9, trajMPC_rilb9, trajMPC_lbIdx9, trajMPC_ubIdx9, trajMPC_grad_ineq9, trajMPC_lubbysub9, trajMPC_llbbyslb9);
info->dobj = info->pobj - info->dgap;
info->rdgap = info->pobj ? info->dgap / info->pobj : 1e6;
if( info->rdgap < 0 ) info->rdgap = -info->rdgap;
if( info->mu < trajMPC_SET_ACC_KKTCOMPL
    && (info->rdgap < trajMPC_SET_ACC_RDGAP || info->dgap < trajMPC_SET_ACC_KKTCOMPL)
    && info->res_eq < trajMPC_SET_ACC_RESEQ
    && info->res_ineq < trajMPC_SET_ACC_RESINEQ ){
exitcode = trajMPC_OPTIMAL; break; }
if( info->it == trajMPC_SET_MAXIT ){
exitcode = trajMPC_MAXITREACHED; break; }
trajMPC_LA_VVADD3_102(trajMPC_grad_cost, trajMPC_grad_eq, trajMPC_grad_ineq, trajMPC_rd);
trajMPC_LA_DIAG_CHOL_LBUB_11_11_5(trajMPC_H0, trajMPC_llbbyslb0, trajMPC_lbIdx0, trajMPC_lubbysub0, trajMPC_ubIdx0, trajMPC_Phi0);
trajMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(trajMPC_Phi0, params->C1, trajMPC_V0);
trajMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(trajMPC_Phi0, trajMPC_D0, trajMPC_W0);
trajMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(trajMPC_W0, trajMPC_V0, trajMPC_Ysd1);
trajMPC_LA_DIAG_FORWARDSUB_11(trajMPC_Phi0, trajMPC_rd0, trajMPC_Lbyrd0);
trajMPC_LA_DIAG_CHOL_LBUB_5_11_5(trajMPC_H1, trajMPC_llbbyslb1, trajMPC_lbIdx1, trajMPC_lubbysub1, trajMPC_ubIdx1, trajMPC_Phi1);
trajMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(trajMPC_Phi1, params->C2, trajMPC_V1);
trajMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(trajMPC_Phi1, trajMPC_D1, trajMPC_W1);
trajMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(trajMPC_W1, trajMPC_V1, trajMPC_Ysd2);
trajMPC_LA_DIAG_FORWARDSUB_11(trajMPC_Phi1, trajMPC_rd1, trajMPC_Lbyrd1);
trajMPC_LA_DIAG_CHOL_LBUB_11_11_5(trajMPC_H1, trajMPC_llbbyslb2, trajMPC_lbIdx2, trajMPC_lubbysub2, trajMPC_ubIdx2, trajMPC_Phi2);
trajMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(trajMPC_Phi2, params->C3, trajMPC_V2);
trajMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(trajMPC_Phi2, trajMPC_D1, trajMPC_W2);
trajMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(trajMPC_W2, trajMPC_V2, trajMPC_Ysd3);
trajMPC_LA_DIAG_FORWARDSUB_11(trajMPC_Phi2, trajMPC_rd2, trajMPC_Lbyrd2);
trajMPC_LA_DIAG_CHOL_LBUB_11_11_5(trajMPC_H1, trajMPC_llbbyslb3, trajMPC_lbIdx3, trajMPC_lubbysub3, trajMPC_ubIdx3, trajMPC_Phi3);
trajMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(trajMPC_Phi3, params->C4, trajMPC_V3);
trajMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(trajMPC_Phi3, trajMPC_D1, trajMPC_W3);
trajMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(trajMPC_W3, trajMPC_V3, trajMPC_Ysd4);
trajMPC_LA_DIAG_FORWARDSUB_11(trajMPC_Phi3, trajMPC_rd3, trajMPC_Lbyrd3);
trajMPC_LA_DIAG_CHOL_LBUB_11_11_5(trajMPC_H1, trajMPC_llbbyslb4, trajMPC_lbIdx4, trajMPC_lubbysub4, trajMPC_ubIdx4, trajMPC_Phi4);
trajMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(trajMPC_Phi4, params->C5, trajMPC_V4);
trajMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(trajMPC_Phi4, trajMPC_D1, trajMPC_W4);
trajMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(trajMPC_W4, trajMPC_V4, trajMPC_Ysd5);
trajMPC_LA_DIAG_FORWARDSUB_11(trajMPC_Phi4, trajMPC_rd4, trajMPC_Lbyrd4);
trajMPC_LA_DIAG_CHOL_LBUB_11_11_5(trajMPC_H1, trajMPC_llbbyslb5, trajMPC_lbIdx5, trajMPC_lubbysub5, trajMPC_ubIdx5, trajMPC_Phi5);
trajMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(trajMPC_Phi5, params->C6, trajMPC_V5);
trajMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(trajMPC_Phi5, trajMPC_D1, trajMPC_W5);
trajMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(trajMPC_W5, trajMPC_V5, trajMPC_Ysd6);
trajMPC_LA_DIAG_FORWARDSUB_11(trajMPC_Phi5, trajMPC_rd5, trajMPC_Lbyrd5);
trajMPC_LA_DIAG_CHOL_LBUB_11_11_5(trajMPC_H1, trajMPC_llbbyslb6, trajMPC_lbIdx6, trajMPC_lubbysub6, trajMPC_ubIdx6, trajMPC_Phi6);
trajMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(trajMPC_Phi6, params->C7, trajMPC_V6);
trajMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(trajMPC_Phi6, trajMPC_D1, trajMPC_W6);
trajMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(trajMPC_W6, trajMPC_V6, trajMPC_Ysd7);
trajMPC_LA_DIAG_FORWARDSUB_11(trajMPC_Phi6, trajMPC_rd6, trajMPC_Lbyrd6);
trajMPC_LA_DIAG_CHOL_LBUB_11_11_5(trajMPC_H1, trajMPC_llbbyslb7, trajMPC_lbIdx7, trajMPC_lubbysub7, trajMPC_ubIdx7, trajMPC_Phi7);
trajMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(trajMPC_Phi7, params->C8, trajMPC_V7);
trajMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(trajMPC_Phi7, trajMPC_D1, trajMPC_W7);
trajMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(trajMPC_W7, trajMPC_V7, trajMPC_Ysd8);
trajMPC_LA_DIAG_FORWARDSUB_11(trajMPC_Phi7, trajMPC_rd7, trajMPC_Lbyrd7);
trajMPC_LA_DIAG_CHOL_LBUB_11_11_5(trajMPC_H1, trajMPC_llbbyslb8, trajMPC_lbIdx8, trajMPC_lubbysub8, trajMPC_ubIdx8, trajMPC_Phi8);
trajMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(trajMPC_Phi8, params->C9, trajMPC_V8);
trajMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(trajMPC_Phi8, trajMPC_D1, trajMPC_W8);
trajMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(trajMPC_W8, trajMPC_V8, trajMPC_Ysd9);
trajMPC_LA_DIAG_FORWARDSUB_11(trajMPC_Phi8, trajMPC_rd8, trajMPC_Lbyrd8);
trajMPC_LA_DIAG_CHOL_ONELOOP_LBUB_3_3_3(trajMPC_H9, trajMPC_llbbyslb9, trajMPC_lbIdx9, trajMPC_lubbysub9, trajMPC_ubIdx9, trajMPC_Phi9);
trajMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_3(trajMPC_Phi9, trajMPC_D9, trajMPC_W9);
trajMPC_LA_DIAG_FORWARDSUB_3(trajMPC_Phi9, trajMPC_rd9, trajMPC_Lbyrd9);
trajMPC_LA_DIAGZERO_MMT_3(trajMPC_W0, trajMPC_Yd0);
trajMPC_LA_DIAGZERO_MVMSUB7_3(trajMPC_W0, trajMPC_Lbyrd0, trajMPC_re0, trajMPC_beta0);
trajMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(trajMPC_V0, trajMPC_W1, trajMPC_Yd1);
trajMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(trajMPC_V0, trajMPC_Lbyrd0, trajMPC_W1, trajMPC_Lbyrd1, trajMPC_re1, trajMPC_beta1);
trajMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(trajMPC_V1, trajMPC_W2, trajMPC_Yd2);
trajMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(trajMPC_V1, trajMPC_Lbyrd1, trajMPC_W2, trajMPC_Lbyrd2, trajMPC_re2, trajMPC_beta2);
trajMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(trajMPC_V2, trajMPC_W3, trajMPC_Yd3);
trajMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(trajMPC_V2, trajMPC_Lbyrd2, trajMPC_W3, trajMPC_Lbyrd3, trajMPC_re3, trajMPC_beta3);
trajMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(trajMPC_V3, trajMPC_W4, trajMPC_Yd4);
trajMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(trajMPC_V3, trajMPC_Lbyrd3, trajMPC_W4, trajMPC_Lbyrd4, trajMPC_re4, trajMPC_beta4);
trajMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(trajMPC_V4, trajMPC_W5, trajMPC_Yd5);
trajMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(trajMPC_V4, trajMPC_Lbyrd4, trajMPC_W5, trajMPC_Lbyrd5, trajMPC_re5, trajMPC_beta5);
trajMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(trajMPC_V5, trajMPC_W6, trajMPC_Yd6);
trajMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(trajMPC_V5, trajMPC_Lbyrd5, trajMPC_W6, trajMPC_Lbyrd6, trajMPC_re6, trajMPC_beta6);
trajMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(trajMPC_V6, trajMPC_W7, trajMPC_Yd7);
trajMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(trajMPC_V6, trajMPC_Lbyrd6, trajMPC_W7, trajMPC_Lbyrd7, trajMPC_re7, trajMPC_beta7);
trajMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(trajMPC_V7, trajMPC_W8, trajMPC_Yd8);
trajMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(trajMPC_V7, trajMPC_Lbyrd7, trajMPC_W8, trajMPC_Lbyrd8, trajMPC_re8, trajMPC_beta8);
trajMPC_LA_DENSE_DIAGZERO_MMT2_3_11_3(trajMPC_V8, trajMPC_W9, trajMPC_Yd9);
trajMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_3(trajMPC_V8, trajMPC_Lbyrd8, trajMPC_W9, trajMPC_Lbyrd9, trajMPC_re9, trajMPC_beta9);
trajMPC_LA_DENSE_CHOL_3(trajMPC_Yd0, trajMPC_Ld0);
trajMPC_LA_DENSE_FORWARDSUB_3(trajMPC_Ld0, trajMPC_beta0, trajMPC_yy0);
trajMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(trajMPC_Ld0, trajMPC_Ysd1, trajMPC_Lsd1);
trajMPC_LA_DENSE_MMTSUB_3_3(trajMPC_Lsd1, trajMPC_Yd1);
trajMPC_LA_DENSE_CHOL_3(trajMPC_Yd1, trajMPC_Ld1);
trajMPC_LA_DENSE_MVMSUB1_3_3(trajMPC_Lsd1, trajMPC_yy0, trajMPC_beta1, trajMPC_bmy1);
trajMPC_LA_DENSE_FORWARDSUB_3(trajMPC_Ld1, trajMPC_bmy1, trajMPC_yy1);
trajMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(trajMPC_Ld1, trajMPC_Ysd2, trajMPC_Lsd2);
trajMPC_LA_DENSE_MMTSUB_3_3(trajMPC_Lsd2, trajMPC_Yd2);
trajMPC_LA_DENSE_CHOL_3(trajMPC_Yd2, trajMPC_Ld2);
trajMPC_LA_DENSE_MVMSUB1_3_3(trajMPC_Lsd2, trajMPC_yy1, trajMPC_beta2, trajMPC_bmy2);
trajMPC_LA_DENSE_FORWARDSUB_3(trajMPC_Ld2, trajMPC_bmy2, trajMPC_yy2);
trajMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(trajMPC_Ld2, trajMPC_Ysd3, trajMPC_Lsd3);
trajMPC_LA_DENSE_MMTSUB_3_3(trajMPC_Lsd3, trajMPC_Yd3);
trajMPC_LA_DENSE_CHOL_3(trajMPC_Yd3, trajMPC_Ld3);
trajMPC_LA_DENSE_MVMSUB1_3_3(trajMPC_Lsd3, trajMPC_yy2, trajMPC_beta3, trajMPC_bmy3);
trajMPC_LA_DENSE_FORWARDSUB_3(trajMPC_Ld3, trajMPC_bmy3, trajMPC_yy3);
trajMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(trajMPC_Ld3, trajMPC_Ysd4, trajMPC_Lsd4);
trajMPC_LA_DENSE_MMTSUB_3_3(trajMPC_Lsd4, trajMPC_Yd4);
trajMPC_LA_DENSE_CHOL_3(trajMPC_Yd4, trajMPC_Ld4);
trajMPC_LA_DENSE_MVMSUB1_3_3(trajMPC_Lsd4, trajMPC_yy3, trajMPC_beta4, trajMPC_bmy4);
trajMPC_LA_DENSE_FORWARDSUB_3(trajMPC_Ld4, trajMPC_bmy4, trajMPC_yy4);
trajMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(trajMPC_Ld4, trajMPC_Ysd5, trajMPC_Lsd5);
trajMPC_LA_DENSE_MMTSUB_3_3(trajMPC_Lsd5, trajMPC_Yd5);
trajMPC_LA_DENSE_CHOL_3(trajMPC_Yd5, trajMPC_Ld5);
trajMPC_LA_DENSE_MVMSUB1_3_3(trajMPC_Lsd5, trajMPC_yy4, trajMPC_beta5, trajMPC_bmy5);
trajMPC_LA_DENSE_FORWARDSUB_3(trajMPC_Ld5, trajMPC_bmy5, trajMPC_yy5);
trajMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(trajMPC_Ld5, trajMPC_Ysd6, trajMPC_Lsd6);
trajMPC_LA_DENSE_MMTSUB_3_3(trajMPC_Lsd6, trajMPC_Yd6);
trajMPC_LA_DENSE_CHOL_3(trajMPC_Yd6, trajMPC_Ld6);
trajMPC_LA_DENSE_MVMSUB1_3_3(trajMPC_Lsd6, trajMPC_yy5, trajMPC_beta6, trajMPC_bmy6);
trajMPC_LA_DENSE_FORWARDSUB_3(trajMPC_Ld6, trajMPC_bmy6, trajMPC_yy6);
trajMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(trajMPC_Ld6, trajMPC_Ysd7, trajMPC_Lsd7);
trajMPC_LA_DENSE_MMTSUB_3_3(trajMPC_Lsd7, trajMPC_Yd7);
trajMPC_LA_DENSE_CHOL_3(trajMPC_Yd7, trajMPC_Ld7);
trajMPC_LA_DENSE_MVMSUB1_3_3(trajMPC_Lsd7, trajMPC_yy6, trajMPC_beta7, trajMPC_bmy7);
trajMPC_LA_DENSE_FORWARDSUB_3(trajMPC_Ld7, trajMPC_bmy7, trajMPC_yy7);
trajMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(trajMPC_Ld7, trajMPC_Ysd8, trajMPC_Lsd8);
trajMPC_LA_DENSE_MMTSUB_3_3(trajMPC_Lsd8, trajMPC_Yd8);
trajMPC_LA_DENSE_CHOL_3(trajMPC_Yd8, trajMPC_Ld8);
trajMPC_LA_DENSE_MVMSUB1_3_3(trajMPC_Lsd8, trajMPC_yy7, trajMPC_beta8, trajMPC_bmy8);
trajMPC_LA_DENSE_FORWARDSUB_3(trajMPC_Ld8, trajMPC_bmy8, trajMPC_yy8);
trajMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(trajMPC_Ld8, trajMPC_Ysd9, trajMPC_Lsd9);
trajMPC_LA_DENSE_MMTSUB_3_3(trajMPC_Lsd9, trajMPC_Yd9);
trajMPC_LA_DENSE_CHOL_3(trajMPC_Yd9, trajMPC_Ld9);
trajMPC_LA_DENSE_MVMSUB1_3_3(trajMPC_Lsd9, trajMPC_yy8, trajMPC_beta9, trajMPC_bmy9);
trajMPC_LA_DENSE_FORWARDSUB_3(trajMPC_Ld9, trajMPC_bmy9, trajMPC_yy9);
trajMPC_LA_DENSE_BACKWARDSUB_3(trajMPC_Ld9, trajMPC_yy9, trajMPC_dvaff9);
trajMPC_LA_DENSE_MTVMSUB_3_3(trajMPC_Lsd9, trajMPC_dvaff9, trajMPC_yy8, trajMPC_bmy8);
trajMPC_LA_DENSE_BACKWARDSUB_3(trajMPC_Ld8, trajMPC_bmy8, trajMPC_dvaff8);
trajMPC_LA_DENSE_MTVMSUB_3_3(trajMPC_Lsd8, trajMPC_dvaff8, trajMPC_yy7, trajMPC_bmy7);
trajMPC_LA_DENSE_BACKWARDSUB_3(trajMPC_Ld7, trajMPC_bmy7, trajMPC_dvaff7);
trajMPC_LA_DENSE_MTVMSUB_3_3(trajMPC_Lsd7, trajMPC_dvaff7, trajMPC_yy6, trajMPC_bmy6);
trajMPC_LA_DENSE_BACKWARDSUB_3(trajMPC_Ld6, trajMPC_bmy6, trajMPC_dvaff6);
trajMPC_LA_DENSE_MTVMSUB_3_3(trajMPC_Lsd6, trajMPC_dvaff6, trajMPC_yy5, trajMPC_bmy5);
trajMPC_LA_DENSE_BACKWARDSUB_3(trajMPC_Ld5, trajMPC_bmy5, trajMPC_dvaff5);
trajMPC_LA_DENSE_MTVMSUB_3_3(trajMPC_Lsd5, trajMPC_dvaff5, trajMPC_yy4, trajMPC_bmy4);
trajMPC_LA_DENSE_BACKWARDSUB_3(trajMPC_Ld4, trajMPC_bmy4, trajMPC_dvaff4);
trajMPC_LA_DENSE_MTVMSUB_3_3(trajMPC_Lsd4, trajMPC_dvaff4, trajMPC_yy3, trajMPC_bmy3);
trajMPC_LA_DENSE_BACKWARDSUB_3(trajMPC_Ld3, trajMPC_bmy3, trajMPC_dvaff3);
trajMPC_LA_DENSE_MTVMSUB_3_3(trajMPC_Lsd3, trajMPC_dvaff3, trajMPC_yy2, trajMPC_bmy2);
trajMPC_LA_DENSE_BACKWARDSUB_3(trajMPC_Ld2, trajMPC_bmy2, trajMPC_dvaff2);
trajMPC_LA_DENSE_MTVMSUB_3_3(trajMPC_Lsd2, trajMPC_dvaff2, trajMPC_yy1, trajMPC_bmy1);
trajMPC_LA_DENSE_BACKWARDSUB_3(trajMPC_Ld1, trajMPC_bmy1, trajMPC_dvaff1);
trajMPC_LA_DENSE_MTVMSUB_3_3(trajMPC_Lsd1, trajMPC_dvaff1, trajMPC_yy0, trajMPC_bmy0);
trajMPC_LA_DENSE_BACKWARDSUB_3(trajMPC_Ld0, trajMPC_bmy0, trajMPC_dvaff0);
trajMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C1, trajMPC_dvaff1, trajMPC_D0, trajMPC_dvaff0, trajMPC_grad_eq0);
trajMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C2, trajMPC_dvaff2, trajMPC_D1, trajMPC_dvaff1, trajMPC_grad_eq1);
trajMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C3, trajMPC_dvaff3, trajMPC_D1, trajMPC_dvaff2, trajMPC_grad_eq2);
trajMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C4, trajMPC_dvaff4, trajMPC_D1, trajMPC_dvaff3, trajMPC_grad_eq3);
trajMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C5, trajMPC_dvaff5, trajMPC_D1, trajMPC_dvaff4, trajMPC_grad_eq4);
trajMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C6, trajMPC_dvaff6, trajMPC_D1, trajMPC_dvaff5, trajMPC_grad_eq5);
trajMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C7, trajMPC_dvaff7, trajMPC_D1, trajMPC_dvaff6, trajMPC_grad_eq6);
trajMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C8, trajMPC_dvaff8, trajMPC_D1, trajMPC_dvaff7, trajMPC_grad_eq7);
trajMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C9, trajMPC_dvaff9, trajMPC_D1, trajMPC_dvaff8, trajMPC_grad_eq8);
trajMPC_LA_DIAGZERO_MTVM_3_3(trajMPC_D9, trajMPC_dvaff9, trajMPC_grad_eq9);
trajMPC_LA_VSUB2_102(trajMPC_rd, trajMPC_grad_eq, trajMPC_rd);
trajMPC_LA_DIAG_FORWARDBACKWARDSUB_11(trajMPC_Phi0, trajMPC_rd0, trajMPC_dzaff0);
trajMPC_LA_DIAG_FORWARDBACKWARDSUB_11(trajMPC_Phi1, trajMPC_rd1, trajMPC_dzaff1);
trajMPC_LA_DIAG_FORWARDBACKWARDSUB_11(trajMPC_Phi2, trajMPC_rd2, trajMPC_dzaff2);
trajMPC_LA_DIAG_FORWARDBACKWARDSUB_11(trajMPC_Phi3, trajMPC_rd3, trajMPC_dzaff3);
trajMPC_LA_DIAG_FORWARDBACKWARDSUB_11(trajMPC_Phi4, trajMPC_rd4, trajMPC_dzaff4);
trajMPC_LA_DIAG_FORWARDBACKWARDSUB_11(trajMPC_Phi5, trajMPC_rd5, trajMPC_dzaff5);
trajMPC_LA_DIAG_FORWARDBACKWARDSUB_11(trajMPC_Phi6, trajMPC_rd6, trajMPC_dzaff6);
trajMPC_LA_DIAG_FORWARDBACKWARDSUB_11(trajMPC_Phi7, trajMPC_rd7, trajMPC_dzaff7);
trajMPC_LA_DIAG_FORWARDBACKWARDSUB_11(trajMPC_Phi8, trajMPC_rd8, trajMPC_dzaff8);
trajMPC_LA_DIAG_FORWARDBACKWARDSUB_3(trajMPC_Phi9, trajMPC_rd9, trajMPC_dzaff9);
trajMPC_LA_VSUB_INDEXED_11(trajMPC_dzaff0, trajMPC_lbIdx0, trajMPC_rilb0, trajMPC_dslbaff0);
trajMPC_LA_VSUB3_11(trajMPC_llbbyslb0, trajMPC_dslbaff0, trajMPC_llb0, trajMPC_dllbaff0);
trajMPC_LA_VSUB2_INDEXED_5(trajMPC_riub0, trajMPC_dzaff0, trajMPC_ubIdx0, trajMPC_dsubaff0);
trajMPC_LA_VSUB3_5(trajMPC_lubbysub0, trajMPC_dsubaff0, trajMPC_lub0, trajMPC_dlubaff0);
trajMPC_LA_VSUB_INDEXED_11(trajMPC_dzaff1, trajMPC_lbIdx1, trajMPC_rilb1, trajMPC_dslbaff1);
trajMPC_LA_VSUB3_11(trajMPC_llbbyslb1, trajMPC_dslbaff1, trajMPC_llb1, trajMPC_dllbaff1);
trajMPC_LA_VSUB2_INDEXED_5(trajMPC_riub1, trajMPC_dzaff1, trajMPC_ubIdx1, trajMPC_dsubaff1);
trajMPC_LA_VSUB3_5(trajMPC_lubbysub1, trajMPC_dsubaff1, trajMPC_lub1, trajMPC_dlubaff1);
trajMPC_LA_VSUB_INDEXED_11(trajMPC_dzaff2, trajMPC_lbIdx2, trajMPC_rilb2, trajMPC_dslbaff2);
trajMPC_LA_VSUB3_11(trajMPC_llbbyslb2, trajMPC_dslbaff2, trajMPC_llb2, trajMPC_dllbaff2);
trajMPC_LA_VSUB2_INDEXED_5(trajMPC_riub2, trajMPC_dzaff2, trajMPC_ubIdx2, trajMPC_dsubaff2);
trajMPC_LA_VSUB3_5(trajMPC_lubbysub2, trajMPC_dsubaff2, trajMPC_lub2, trajMPC_dlubaff2);
trajMPC_LA_VSUB_INDEXED_11(trajMPC_dzaff3, trajMPC_lbIdx3, trajMPC_rilb3, trajMPC_dslbaff3);
trajMPC_LA_VSUB3_11(trajMPC_llbbyslb3, trajMPC_dslbaff3, trajMPC_llb3, trajMPC_dllbaff3);
trajMPC_LA_VSUB2_INDEXED_5(trajMPC_riub3, trajMPC_dzaff3, trajMPC_ubIdx3, trajMPC_dsubaff3);
trajMPC_LA_VSUB3_5(trajMPC_lubbysub3, trajMPC_dsubaff3, trajMPC_lub3, trajMPC_dlubaff3);
trajMPC_LA_VSUB_INDEXED_11(trajMPC_dzaff4, trajMPC_lbIdx4, trajMPC_rilb4, trajMPC_dslbaff4);
trajMPC_LA_VSUB3_11(trajMPC_llbbyslb4, trajMPC_dslbaff4, trajMPC_llb4, trajMPC_dllbaff4);
trajMPC_LA_VSUB2_INDEXED_5(trajMPC_riub4, trajMPC_dzaff4, trajMPC_ubIdx4, trajMPC_dsubaff4);
trajMPC_LA_VSUB3_5(trajMPC_lubbysub4, trajMPC_dsubaff4, trajMPC_lub4, trajMPC_dlubaff4);
trajMPC_LA_VSUB_INDEXED_11(trajMPC_dzaff5, trajMPC_lbIdx5, trajMPC_rilb5, trajMPC_dslbaff5);
trajMPC_LA_VSUB3_11(trajMPC_llbbyslb5, trajMPC_dslbaff5, trajMPC_llb5, trajMPC_dllbaff5);
trajMPC_LA_VSUB2_INDEXED_5(trajMPC_riub5, trajMPC_dzaff5, trajMPC_ubIdx5, trajMPC_dsubaff5);
trajMPC_LA_VSUB3_5(trajMPC_lubbysub5, trajMPC_dsubaff5, trajMPC_lub5, trajMPC_dlubaff5);
trajMPC_LA_VSUB_INDEXED_11(trajMPC_dzaff6, trajMPC_lbIdx6, trajMPC_rilb6, trajMPC_dslbaff6);
trajMPC_LA_VSUB3_11(trajMPC_llbbyslb6, trajMPC_dslbaff6, trajMPC_llb6, trajMPC_dllbaff6);
trajMPC_LA_VSUB2_INDEXED_5(trajMPC_riub6, trajMPC_dzaff6, trajMPC_ubIdx6, trajMPC_dsubaff6);
trajMPC_LA_VSUB3_5(trajMPC_lubbysub6, trajMPC_dsubaff6, trajMPC_lub6, trajMPC_dlubaff6);
trajMPC_LA_VSUB_INDEXED_11(trajMPC_dzaff7, trajMPC_lbIdx7, trajMPC_rilb7, trajMPC_dslbaff7);
trajMPC_LA_VSUB3_11(trajMPC_llbbyslb7, trajMPC_dslbaff7, trajMPC_llb7, trajMPC_dllbaff7);
trajMPC_LA_VSUB2_INDEXED_5(trajMPC_riub7, trajMPC_dzaff7, trajMPC_ubIdx7, trajMPC_dsubaff7);
trajMPC_LA_VSUB3_5(trajMPC_lubbysub7, trajMPC_dsubaff7, trajMPC_lub7, trajMPC_dlubaff7);
trajMPC_LA_VSUB_INDEXED_11(trajMPC_dzaff8, trajMPC_lbIdx8, trajMPC_rilb8, trajMPC_dslbaff8);
trajMPC_LA_VSUB3_11(trajMPC_llbbyslb8, trajMPC_dslbaff8, trajMPC_llb8, trajMPC_dllbaff8);
trajMPC_LA_VSUB2_INDEXED_5(trajMPC_riub8, trajMPC_dzaff8, trajMPC_ubIdx8, trajMPC_dsubaff8);
trajMPC_LA_VSUB3_5(trajMPC_lubbysub8, trajMPC_dsubaff8, trajMPC_lub8, trajMPC_dlubaff8);
trajMPC_LA_VSUB_INDEXED_3(trajMPC_dzaff9, trajMPC_lbIdx9, trajMPC_rilb9, trajMPC_dslbaff9);
trajMPC_LA_VSUB3_3(trajMPC_llbbyslb9, trajMPC_dslbaff9, trajMPC_llb9, trajMPC_dllbaff9);
trajMPC_LA_VSUB2_INDEXED_3(trajMPC_riub9, trajMPC_dzaff9, trajMPC_ubIdx9, trajMPC_dsubaff9);
trajMPC_LA_VSUB3_3(trajMPC_lubbysub9, trajMPC_dsubaff9, trajMPC_lub9, trajMPC_dlubaff9);
info->lsit_aff = trajMPC_LINESEARCH_BACKTRACKING_AFFINE(trajMPC_l, trajMPC_s, trajMPC_dl_aff, trajMPC_ds_aff, &info->step_aff, &info->mu_aff);
if( info->lsit_aff == trajMPC_NOPROGRESS ){
exitcode = trajMPC_NOPROGRESS; break;
}
sigma_3rdroot = info->mu_aff / info->mu;
info->sigma = sigma_3rdroot*sigma_3rdroot*sigma_3rdroot;
musigma = info->mu * info->sigma;
trajMPC_LA_VSUB5_150(trajMPC_ds_aff, trajMPC_dl_aff, musigma, trajMPC_ccrhs);
trajMPC_LA_VSUB6_INDEXED_11_5_11(trajMPC_ccrhsub0, trajMPC_sub0, trajMPC_ubIdx0, trajMPC_ccrhsl0, trajMPC_slb0, trajMPC_lbIdx0, trajMPC_rd0);
trajMPC_LA_VSUB6_INDEXED_11_5_11(trajMPC_ccrhsub1, trajMPC_sub1, trajMPC_ubIdx1, trajMPC_ccrhsl1, trajMPC_slb1, trajMPC_lbIdx1, trajMPC_rd1);
trajMPC_LA_DIAG_FORWARDSUB_11(trajMPC_Phi0, trajMPC_rd0, trajMPC_Lbyrd0);
trajMPC_LA_DIAG_FORWARDSUB_11(trajMPC_Phi1, trajMPC_rd1, trajMPC_Lbyrd1);
trajMPC_LA_DIAGZERO_MVM_3(trajMPC_W0, trajMPC_Lbyrd0, trajMPC_beta0);
trajMPC_LA_DENSE_FORWARDSUB_3(trajMPC_Ld0, trajMPC_beta0, trajMPC_yy0);
trajMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(trajMPC_V0, trajMPC_Lbyrd0, trajMPC_W1, trajMPC_Lbyrd1, trajMPC_beta1);
trajMPC_LA_DENSE_MVMSUB1_3_3(trajMPC_Lsd1, trajMPC_yy0, trajMPC_beta1, trajMPC_bmy1);
trajMPC_LA_DENSE_FORWARDSUB_3(trajMPC_Ld1, trajMPC_bmy1, trajMPC_yy1);
trajMPC_LA_VSUB6_INDEXED_11_5_11(trajMPC_ccrhsub2, trajMPC_sub2, trajMPC_ubIdx2, trajMPC_ccrhsl2, trajMPC_slb2, trajMPC_lbIdx2, trajMPC_rd2);
trajMPC_LA_DIAG_FORWARDSUB_11(trajMPC_Phi2, trajMPC_rd2, trajMPC_Lbyrd2);
trajMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(trajMPC_V1, trajMPC_Lbyrd1, trajMPC_W2, trajMPC_Lbyrd2, trajMPC_beta2);
trajMPC_LA_DENSE_MVMSUB1_3_3(trajMPC_Lsd2, trajMPC_yy1, trajMPC_beta2, trajMPC_bmy2);
trajMPC_LA_DENSE_FORWARDSUB_3(trajMPC_Ld2, trajMPC_bmy2, trajMPC_yy2);
trajMPC_LA_VSUB6_INDEXED_11_5_11(trajMPC_ccrhsub3, trajMPC_sub3, trajMPC_ubIdx3, trajMPC_ccrhsl3, trajMPC_slb3, trajMPC_lbIdx3, trajMPC_rd3);
trajMPC_LA_DIAG_FORWARDSUB_11(trajMPC_Phi3, trajMPC_rd3, trajMPC_Lbyrd3);
trajMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(trajMPC_V2, trajMPC_Lbyrd2, trajMPC_W3, trajMPC_Lbyrd3, trajMPC_beta3);
trajMPC_LA_DENSE_MVMSUB1_3_3(trajMPC_Lsd3, trajMPC_yy2, trajMPC_beta3, trajMPC_bmy3);
trajMPC_LA_DENSE_FORWARDSUB_3(trajMPC_Ld3, trajMPC_bmy3, trajMPC_yy3);
trajMPC_LA_VSUB6_INDEXED_11_5_11(trajMPC_ccrhsub4, trajMPC_sub4, trajMPC_ubIdx4, trajMPC_ccrhsl4, trajMPC_slb4, trajMPC_lbIdx4, trajMPC_rd4);
trajMPC_LA_DIAG_FORWARDSUB_11(trajMPC_Phi4, trajMPC_rd4, trajMPC_Lbyrd4);
trajMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(trajMPC_V3, trajMPC_Lbyrd3, trajMPC_W4, trajMPC_Lbyrd4, trajMPC_beta4);
trajMPC_LA_DENSE_MVMSUB1_3_3(trajMPC_Lsd4, trajMPC_yy3, trajMPC_beta4, trajMPC_bmy4);
trajMPC_LA_DENSE_FORWARDSUB_3(trajMPC_Ld4, trajMPC_bmy4, trajMPC_yy4);
trajMPC_LA_VSUB6_INDEXED_11_5_11(trajMPC_ccrhsub5, trajMPC_sub5, trajMPC_ubIdx5, trajMPC_ccrhsl5, trajMPC_slb5, trajMPC_lbIdx5, trajMPC_rd5);
trajMPC_LA_DIAG_FORWARDSUB_11(trajMPC_Phi5, trajMPC_rd5, trajMPC_Lbyrd5);
trajMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(trajMPC_V4, trajMPC_Lbyrd4, trajMPC_W5, trajMPC_Lbyrd5, trajMPC_beta5);
trajMPC_LA_DENSE_MVMSUB1_3_3(trajMPC_Lsd5, trajMPC_yy4, trajMPC_beta5, trajMPC_bmy5);
trajMPC_LA_DENSE_FORWARDSUB_3(trajMPC_Ld5, trajMPC_bmy5, trajMPC_yy5);
trajMPC_LA_VSUB6_INDEXED_11_5_11(trajMPC_ccrhsub6, trajMPC_sub6, trajMPC_ubIdx6, trajMPC_ccrhsl6, trajMPC_slb6, trajMPC_lbIdx6, trajMPC_rd6);
trajMPC_LA_DIAG_FORWARDSUB_11(trajMPC_Phi6, trajMPC_rd6, trajMPC_Lbyrd6);
trajMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(trajMPC_V5, trajMPC_Lbyrd5, trajMPC_W6, trajMPC_Lbyrd6, trajMPC_beta6);
trajMPC_LA_DENSE_MVMSUB1_3_3(trajMPC_Lsd6, trajMPC_yy5, trajMPC_beta6, trajMPC_bmy6);
trajMPC_LA_DENSE_FORWARDSUB_3(trajMPC_Ld6, trajMPC_bmy6, trajMPC_yy6);
trajMPC_LA_VSUB6_INDEXED_11_5_11(trajMPC_ccrhsub7, trajMPC_sub7, trajMPC_ubIdx7, trajMPC_ccrhsl7, trajMPC_slb7, trajMPC_lbIdx7, trajMPC_rd7);
trajMPC_LA_DIAG_FORWARDSUB_11(trajMPC_Phi7, trajMPC_rd7, trajMPC_Lbyrd7);
trajMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(trajMPC_V6, trajMPC_Lbyrd6, trajMPC_W7, trajMPC_Lbyrd7, trajMPC_beta7);
trajMPC_LA_DENSE_MVMSUB1_3_3(trajMPC_Lsd7, trajMPC_yy6, trajMPC_beta7, trajMPC_bmy7);
trajMPC_LA_DENSE_FORWARDSUB_3(trajMPC_Ld7, trajMPC_bmy7, trajMPC_yy7);
trajMPC_LA_VSUB6_INDEXED_11_5_11(trajMPC_ccrhsub8, trajMPC_sub8, trajMPC_ubIdx8, trajMPC_ccrhsl8, trajMPC_slb8, trajMPC_lbIdx8, trajMPC_rd8);
trajMPC_LA_DIAG_FORWARDSUB_11(trajMPC_Phi8, trajMPC_rd8, trajMPC_Lbyrd8);
trajMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(trajMPC_V7, trajMPC_Lbyrd7, trajMPC_W8, trajMPC_Lbyrd8, trajMPC_beta8);
trajMPC_LA_DENSE_MVMSUB1_3_3(trajMPC_Lsd8, trajMPC_yy7, trajMPC_beta8, trajMPC_bmy8);
trajMPC_LA_DENSE_FORWARDSUB_3(trajMPC_Ld8, trajMPC_bmy8, trajMPC_yy8);
trajMPC_LA_VSUB6_INDEXED_3_3_3(trajMPC_ccrhsub9, trajMPC_sub9, trajMPC_ubIdx9, trajMPC_ccrhsl9, trajMPC_slb9, trajMPC_lbIdx9, trajMPC_rd9);
trajMPC_LA_DIAG_FORWARDSUB_3(trajMPC_Phi9, trajMPC_rd9, trajMPC_Lbyrd9);
trajMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_3(trajMPC_V8, trajMPC_Lbyrd8, trajMPC_W9, trajMPC_Lbyrd9, trajMPC_beta9);
trajMPC_LA_DENSE_MVMSUB1_3_3(trajMPC_Lsd9, trajMPC_yy8, trajMPC_beta9, trajMPC_bmy9);
trajMPC_LA_DENSE_FORWARDSUB_3(trajMPC_Ld9, trajMPC_bmy9, trajMPC_yy9);
trajMPC_LA_DENSE_BACKWARDSUB_3(trajMPC_Ld9, trajMPC_yy9, trajMPC_dvcc9);
trajMPC_LA_DENSE_MTVMSUB_3_3(trajMPC_Lsd9, trajMPC_dvcc9, trajMPC_yy8, trajMPC_bmy8);
trajMPC_LA_DENSE_BACKWARDSUB_3(trajMPC_Ld8, trajMPC_bmy8, trajMPC_dvcc8);
trajMPC_LA_DENSE_MTVMSUB_3_3(trajMPC_Lsd8, trajMPC_dvcc8, trajMPC_yy7, trajMPC_bmy7);
trajMPC_LA_DENSE_BACKWARDSUB_3(trajMPC_Ld7, trajMPC_bmy7, trajMPC_dvcc7);
trajMPC_LA_DENSE_MTVMSUB_3_3(trajMPC_Lsd7, trajMPC_dvcc7, trajMPC_yy6, trajMPC_bmy6);
trajMPC_LA_DENSE_BACKWARDSUB_3(trajMPC_Ld6, trajMPC_bmy6, trajMPC_dvcc6);
trajMPC_LA_DENSE_MTVMSUB_3_3(trajMPC_Lsd6, trajMPC_dvcc6, trajMPC_yy5, trajMPC_bmy5);
trajMPC_LA_DENSE_BACKWARDSUB_3(trajMPC_Ld5, trajMPC_bmy5, trajMPC_dvcc5);
trajMPC_LA_DENSE_MTVMSUB_3_3(trajMPC_Lsd5, trajMPC_dvcc5, trajMPC_yy4, trajMPC_bmy4);
trajMPC_LA_DENSE_BACKWARDSUB_3(trajMPC_Ld4, trajMPC_bmy4, trajMPC_dvcc4);
trajMPC_LA_DENSE_MTVMSUB_3_3(trajMPC_Lsd4, trajMPC_dvcc4, trajMPC_yy3, trajMPC_bmy3);
trajMPC_LA_DENSE_BACKWARDSUB_3(trajMPC_Ld3, trajMPC_bmy3, trajMPC_dvcc3);
trajMPC_LA_DENSE_MTVMSUB_3_3(trajMPC_Lsd3, trajMPC_dvcc3, trajMPC_yy2, trajMPC_bmy2);
trajMPC_LA_DENSE_BACKWARDSUB_3(trajMPC_Ld2, trajMPC_bmy2, trajMPC_dvcc2);
trajMPC_LA_DENSE_MTVMSUB_3_3(trajMPC_Lsd2, trajMPC_dvcc2, trajMPC_yy1, trajMPC_bmy1);
trajMPC_LA_DENSE_BACKWARDSUB_3(trajMPC_Ld1, trajMPC_bmy1, trajMPC_dvcc1);
trajMPC_LA_DENSE_MTVMSUB_3_3(trajMPC_Lsd1, trajMPC_dvcc1, trajMPC_yy0, trajMPC_bmy0);
trajMPC_LA_DENSE_BACKWARDSUB_3(trajMPC_Ld0, trajMPC_bmy0, trajMPC_dvcc0);
trajMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C1, trajMPC_dvcc1, trajMPC_D0, trajMPC_dvcc0, trajMPC_grad_eq0);
trajMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C2, trajMPC_dvcc2, trajMPC_D1, trajMPC_dvcc1, trajMPC_grad_eq1);
trajMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C3, trajMPC_dvcc3, trajMPC_D1, trajMPC_dvcc2, trajMPC_grad_eq2);
trajMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C4, trajMPC_dvcc4, trajMPC_D1, trajMPC_dvcc3, trajMPC_grad_eq3);
trajMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C5, trajMPC_dvcc5, trajMPC_D1, trajMPC_dvcc4, trajMPC_grad_eq4);
trajMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C6, trajMPC_dvcc6, trajMPC_D1, trajMPC_dvcc5, trajMPC_grad_eq5);
trajMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C7, trajMPC_dvcc7, trajMPC_D1, trajMPC_dvcc6, trajMPC_grad_eq6);
trajMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C8, trajMPC_dvcc8, trajMPC_D1, trajMPC_dvcc7, trajMPC_grad_eq7);
trajMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C9, trajMPC_dvcc9, trajMPC_D1, trajMPC_dvcc8, trajMPC_grad_eq8);
trajMPC_LA_DIAGZERO_MTVM_3_3(trajMPC_D9, trajMPC_dvcc9, trajMPC_grad_eq9);
trajMPC_LA_VSUB_102(trajMPC_rd, trajMPC_grad_eq, trajMPC_rd);
trajMPC_LA_DIAG_FORWARDBACKWARDSUB_11(trajMPC_Phi0, trajMPC_rd0, trajMPC_dzcc0);
trajMPC_LA_DIAG_FORWARDBACKWARDSUB_11(trajMPC_Phi1, trajMPC_rd1, trajMPC_dzcc1);
trajMPC_LA_DIAG_FORWARDBACKWARDSUB_11(trajMPC_Phi2, trajMPC_rd2, trajMPC_dzcc2);
trajMPC_LA_DIAG_FORWARDBACKWARDSUB_11(trajMPC_Phi3, trajMPC_rd3, trajMPC_dzcc3);
trajMPC_LA_DIAG_FORWARDBACKWARDSUB_11(trajMPC_Phi4, trajMPC_rd4, trajMPC_dzcc4);
trajMPC_LA_DIAG_FORWARDBACKWARDSUB_11(trajMPC_Phi5, trajMPC_rd5, trajMPC_dzcc5);
trajMPC_LA_DIAG_FORWARDBACKWARDSUB_11(trajMPC_Phi6, trajMPC_rd6, trajMPC_dzcc6);
trajMPC_LA_DIAG_FORWARDBACKWARDSUB_11(trajMPC_Phi7, trajMPC_rd7, trajMPC_dzcc7);
trajMPC_LA_DIAG_FORWARDBACKWARDSUB_11(trajMPC_Phi8, trajMPC_rd8, trajMPC_dzcc8);
trajMPC_LA_DIAG_FORWARDBACKWARDSUB_3(trajMPC_Phi9, trajMPC_rd9, trajMPC_dzcc9);
trajMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(trajMPC_ccrhsl0, trajMPC_slb0, trajMPC_llbbyslb0, trajMPC_dzcc0, trajMPC_lbIdx0, trajMPC_dllbcc0);
trajMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(trajMPC_ccrhsub0, trajMPC_sub0, trajMPC_lubbysub0, trajMPC_dzcc0, trajMPC_ubIdx0, trajMPC_dlubcc0);
trajMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(trajMPC_ccrhsl1, trajMPC_slb1, trajMPC_llbbyslb1, trajMPC_dzcc1, trajMPC_lbIdx1, trajMPC_dllbcc1);
trajMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(trajMPC_ccrhsub1, trajMPC_sub1, trajMPC_lubbysub1, trajMPC_dzcc1, trajMPC_ubIdx1, trajMPC_dlubcc1);
trajMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(trajMPC_ccrhsl2, trajMPC_slb2, trajMPC_llbbyslb2, trajMPC_dzcc2, trajMPC_lbIdx2, trajMPC_dllbcc2);
trajMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(trajMPC_ccrhsub2, trajMPC_sub2, trajMPC_lubbysub2, trajMPC_dzcc2, trajMPC_ubIdx2, trajMPC_dlubcc2);
trajMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(trajMPC_ccrhsl3, trajMPC_slb3, trajMPC_llbbyslb3, trajMPC_dzcc3, trajMPC_lbIdx3, trajMPC_dllbcc3);
trajMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(trajMPC_ccrhsub3, trajMPC_sub3, trajMPC_lubbysub3, trajMPC_dzcc3, trajMPC_ubIdx3, trajMPC_dlubcc3);
trajMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(trajMPC_ccrhsl4, trajMPC_slb4, trajMPC_llbbyslb4, trajMPC_dzcc4, trajMPC_lbIdx4, trajMPC_dllbcc4);
trajMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(trajMPC_ccrhsub4, trajMPC_sub4, trajMPC_lubbysub4, trajMPC_dzcc4, trajMPC_ubIdx4, trajMPC_dlubcc4);
trajMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(trajMPC_ccrhsl5, trajMPC_slb5, trajMPC_llbbyslb5, trajMPC_dzcc5, trajMPC_lbIdx5, trajMPC_dllbcc5);
trajMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(trajMPC_ccrhsub5, trajMPC_sub5, trajMPC_lubbysub5, trajMPC_dzcc5, trajMPC_ubIdx5, trajMPC_dlubcc5);
trajMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(trajMPC_ccrhsl6, trajMPC_slb6, trajMPC_llbbyslb6, trajMPC_dzcc6, trajMPC_lbIdx6, trajMPC_dllbcc6);
trajMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(trajMPC_ccrhsub6, trajMPC_sub6, trajMPC_lubbysub6, trajMPC_dzcc6, trajMPC_ubIdx6, trajMPC_dlubcc6);
trajMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(trajMPC_ccrhsl7, trajMPC_slb7, trajMPC_llbbyslb7, trajMPC_dzcc7, trajMPC_lbIdx7, trajMPC_dllbcc7);
trajMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(trajMPC_ccrhsub7, trajMPC_sub7, trajMPC_lubbysub7, trajMPC_dzcc7, trajMPC_ubIdx7, trajMPC_dlubcc7);
trajMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(trajMPC_ccrhsl8, trajMPC_slb8, trajMPC_llbbyslb8, trajMPC_dzcc8, trajMPC_lbIdx8, trajMPC_dllbcc8);
trajMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(trajMPC_ccrhsub8, trajMPC_sub8, trajMPC_lubbysub8, trajMPC_dzcc8, trajMPC_ubIdx8, trajMPC_dlubcc8);
trajMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_3(trajMPC_ccrhsl9, trajMPC_slb9, trajMPC_llbbyslb9, trajMPC_dzcc9, trajMPC_lbIdx9, trajMPC_dllbcc9);
trajMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_3(trajMPC_ccrhsub9, trajMPC_sub9, trajMPC_lubbysub9, trajMPC_dzcc9, trajMPC_ubIdx9, trajMPC_dlubcc9);
trajMPC_LA_VSUB7_150(trajMPC_l, trajMPC_ccrhs, trajMPC_s, trajMPC_dl_cc, trajMPC_ds_cc);
trajMPC_LA_VADD_102(trajMPC_dz_cc, trajMPC_dz_aff);
trajMPC_LA_VADD_30(trajMPC_dv_cc, trajMPC_dv_aff);
trajMPC_LA_VADD_150(trajMPC_dl_cc, trajMPC_dl_aff);
trajMPC_LA_VADD_150(trajMPC_ds_cc, trajMPC_ds_aff);
info->lsit_cc = trajMPC_LINESEARCH_BACKTRACKING_COMBINED(trajMPC_z, trajMPC_v, trajMPC_l, trajMPC_s, trajMPC_dz_cc, trajMPC_dv_cc, trajMPC_dl_cc, trajMPC_ds_cc, &info->step_cc, &info->mu);
if( info->lsit_cc == trajMPC_NOPROGRESS ){
exitcode = trajMPC_NOPROGRESS; break;
}
info->it++;
}
output->z1[0] = trajMPC_z0[0];
output->z1[1] = trajMPC_z0[1];
output->z1[2] = trajMPC_z0[2];
output->z1[3] = trajMPC_z0[3];
output->z1[4] = trajMPC_z0[4];
output->z2[0] = trajMPC_z1[0];
output->z2[1] = trajMPC_z1[1];
output->z2[2] = trajMPC_z1[2];
output->z2[3] = trajMPC_z1[3];
output->z2[4] = trajMPC_z1[4];
output->z3[0] = trajMPC_z2[0];
output->z3[1] = trajMPC_z2[1];
output->z3[2] = trajMPC_z2[2];
output->z3[3] = trajMPC_z2[3];
output->z3[4] = trajMPC_z2[4];
output->z4[0] = trajMPC_z3[0];
output->z4[1] = trajMPC_z3[1];
output->z4[2] = trajMPC_z3[2];
output->z4[3] = trajMPC_z3[3];
output->z4[4] = trajMPC_z3[4];
output->z5[0] = trajMPC_z4[0];
output->z5[1] = trajMPC_z4[1];
output->z5[2] = trajMPC_z4[2];
output->z5[3] = trajMPC_z4[3];
output->z5[4] = trajMPC_z4[4];
output->z6[0] = trajMPC_z5[0];
output->z6[1] = trajMPC_z5[1];
output->z6[2] = trajMPC_z5[2];
output->z6[3] = trajMPC_z5[3];
output->z6[4] = trajMPC_z5[4];
output->z7[0] = trajMPC_z6[0];
output->z7[1] = trajMPC_z6[1];
output->z7[2] = trajMPC_z6[2];
output->z7[3] = trajMPC_z6[3];
output->z7[4] = trajMPC_z6[4];
output->z8[0] = trajMPC_z7[0];
output->z8[1] = trajMPC_z7[1];
output->z8[2] = trajMPC_z7[2];
output->z8[3] = trajMPC_z7[3];
output->z8[4] = trajMPC_z7[4];
output->z9[0] = trajMPC_z8[0];
output->z9[1] = trajMPC_z8[1];
output->z9[2] = trajMPC_z8[2];
output->z9[3] = trajMPC_z8[3];
output->z9[4] = trajMPC_z8[4];
output->z10[0] = trajMPC_z9[0];
output->z10[1] = trajMPC_z9[1];
output->z10[2] = trajMPC_z9[2];

#if trajMPC_SET_TIMING == 1
info->solvetime = trajMPC_toc(&solvertimer);
#if trajMPC_SET_PRINTLEVEL > 0 && trajMPC_SET_TIMING == 1
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
