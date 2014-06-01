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
 * Initializes a vector of length 157 with a value.
 */
void trajMPC_LA_INITIALIZEVECTOR_157(trajMPC_FLOAT* vec, trajMPC_FLOAT value)
{
	int i;
	for( i=0; i<157; i++ )
	{
		vec[i] = value;
	}
}


/*
 * Initializes a vector of length 45 with a value.
 */
void trajMPC_LA_INITIALIZEVECTOR_45(trajMPC_FLOAT* vec, trajMPC_FLOAT value)
{
	int i;
	for( i=0; i<45; i++ )
	{
		vec[i] = value;
	}
}


/*
 * Initializes a vector of length 230 with a value.
 */
void trajMPC_LA_INITIALIZEVECTOR_230(trajMPC_FLOAT* vec, trajMPC_FLOAT value)
{
	int i;
	for( i=0; i<230; i++ )
	{
		vec[i] = value;
	}
}


/* 
 * Calculates a dot product and adds it to a variable: z += x'*y; 
 * This function is for vectors of length 230.
 */
void trajMPC_LA_DOTACC_230(trajMPC_FLOAT *x, trajMPC_FLOAT *y, trajMPC_FLOAT *z)
{
	int i;
	for( i=0; i<230; i++ ){
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
 * of length 157.
 */
void trajMPC_LA_VVADD3_157(trajMPC_FLOAT *u, trajMPC_FLOAT *v, trajMPC_FLOAT *w, trajMPC_FLOAT *z)
{
	int i;
	for( i=0; i<157; i++ ){
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
 * Vector subtraction z = -x - y for vectors of length 157.
 */
void trajMPC_LA_VSUB2_157(trajMPC_FLOAT *x, trajMPC_FLOAT *y, trajMPC_FLOAT *z)
{
	int i;
	for( i=0; i<157; i++){
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
        for( i=0; i<230; i++ ){
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
        if( i == 230 ){
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
    *mu_aff = mymu / (trajMPC_FLOAT)230;
    return lsIt;
}


/*
 * Vector subtraction x = u.*v - a where a is a scalar
*  and x,u,v are vectors of length 230.
 */
void trajMPC_LA_VSUB5_230(trajMPC_FLOAT *u, trajMPC_FLOAT *v, trajMPC_FLOAT a, trajMPC_FLOAT *x)
{
	int i;
	for( i=0; i<230; i++){
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
 * Vector subtraction z = x - y for vectors of length 157.
 */
void trajMPC_LA_VSUB_157(trajMPC_FLOAT *x, trajMPC_FLOAT *y, trajMPC_FLOAT *z)
{
	int i;
	for( i=0; i<157; i++){
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
 * Computes ds = -l.\(r + s.*dl) for vectors of length 230.
 */
void trajMPC_LA_VSUB7_230(trajMPC_FLOAT *l, trajMPC_FLOAT *r, trajMPC_FLOAT *s, trajMPC_FLOAT *dl, trajMPC_FLOAT *ds)
{
	int i;
	for( i=0; i<230; i++){
		ds[i] = -(r[i] + s[i]*dl[i])/l[i];
	}
}


/*
 * Vector addition x = x + y for vectors of length 157.
 */
void trajMPC_LA_VADD_157(trajMPC_FLOAT *x, trajMPC_FLOAT *y)
{
	int i;
	for( i=0; i<157; i++){
		x[i] += y[i];
	}
}


/*
 * Vector addition x = x + y for vectors of length 45.
 */
void trajMPC_LA_VADD_45(trajMPC_FLOAT *x, trajMPC_FLOAT *y)
{
	int i;
	for( i=0; i<45; i++){
		x[i] += y[i];
	}
}


/*
 * Vector addition x = x + y for vectors of length 230.
 */
void trajMPC_LA_VADD_230(trajMPC_FLOAT *x, trajMPC_FLOAT *y)
{
	int i;
	for( i=0; i<230; i++){
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
        for( i=0; i<230; i++ ){
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
        if( i == 230 ){
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
    for( i=0; i<157; i++ ){
        z[i] += a_gamma*dz[i];
    }
    
    /* equality constraint multipliers */
    for( i=0; i<45; i++ ){
        v[i] += a_gamma*dv[i];
    }
    
    /* inequality constraint multipliers & slacks, also update mu */
    *mu = 0;
    for( i=0; i<230; i++ ){
        dltemp = l[i] + a_gamma*dl[i]; l[i] = dltemp;
        dstemp = s[i] + a_gamma*ds[i]; s[i] = dstemp;
        *mu += dltemp*dstemp;
    }
    
    *a = a_gamma;
    *mu /= (trajMPC_FLOAT)230;
    return lsIt;
}




/* VARIABLE DEFINITIONS ------------------------------------------------ */
trajMPC_FLOAT trajMPC_z[157];
trajMPC_FLOAT trajMPC_v[45];
trajMPC_FLOAT trajMPC_dz_aff[157];
trajMPC_FLOAT trajMPC_dv_aff[45];
trajMPC_FLOAT trajMPC_grad_cost[157];
trajMPC_FLOAT trajMPC_grad_eq[157];
trajMPC_FLOAT trajMPC_rd[157];
trajMPC_FLOAT trajMPC_l[230];
trajMPC_FLOAT trajMPC_s[230];
trajMPC_FLOAT trajMPC_lbys[230];
trajMPC_FLOAT trajMPC_dl_aff[230];
trajMPC_FLOAT trajMPC_ds_aff[230];
trajMPC_FLOAT trajMPC_dz_cc[157];
trajMPC_FLOAT trajMPC_dv_cc[45];
trajMPC_FLOAT trajMPC_dl_cc[230];
trajMPC_FLOAT trajMPC_ds_cc[230];
trajMPC_FLOAT trajMPC_ccrhs[230];
trajMPC_FLOAT trajMPC_grad_ineq[157];
trajMPC_FLOAT trajMPC_H00[11] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 2.0000000000000000E+000, 2.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
trajMPC_FLOAT* trajMPC_z00 = trajMPC_z + 0;
trajMPC_FLOAT* trajMPC_dzaff00 = trajMPC_dz_aff + 0;
trajMPC_FLOAT* trajMPC_dzcc00 = trajMPC_dz_cc + 0;
trajMPC_FLOAT* trajMPC_rd00 = trajMPC_rd + 0;
trajMPC_FLOAT trajMPC_Lbyrd00[11];
trajMPC_FLOAT* trajMPC_grad_cost00 = trajMPC_grad_cost + 0;
trajMPC_FLOAT* trajMPC_grad_eq00 = trajMPC_grad_eq + 0;
trajMPC_FLOAT* trajMPC_grad_ineq00 = trajMPC_grad_ineq + 0;
trajMPC_FLOAT trajMPC_ctv00[11];
trajMPC_FLOAT* trajMPC_v00 = trajMPC_v + 0;
trajMPC_FLOAT trajMPC_re00[3];
trajMPC_FLOAT trajMPC_beta00[3];
trajMPC_FLOAT trajMPC_betacc00[3];
trajMPC_FLOAT* trajMPC_dvaff00 = trajMPC_dv_aff + 0;
trajMPC_FLOAT* trajMPC_dvcc00 = trajMPC_dv_cc + 0;
trajMPC_FLOAT trajMPC_V00[33];
trajMPC_FLOAT trajMPC_Yd00[6];
trajMPC_FLOAT trajMPC_Ld00[6];
trajMPC_FLOAT trajMPC_yy00[3];
trajMPC_FLOAT trajMPC_bmy00[3];
int trajMPC_lbIdx00[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
trajMPC_FLOAT* trajMPC_llb00 = trajMPC_l + 0;
trajMPC_FLOAT* trajMPC_slb00 = trajMPC_s + 0;
trajMPC_FLOAT* trajMPC_llbbyslb00 = trajMPC_lbys + 0;
trajMPC_FLOAT trajMPC_rilb00[11];
trajMPC_FLOAT* trajMPC_dllbaff00 = trajMPC_dl_aff + 0;
trajMPC_FLOAT* trajMPC_dslbaff00 = trajMPC_ds_aff + 0;
trajMPC_FLOAT* trajMPC_dllbcc00 = trajMPC_dl_cc + 0;
trajMPC_FLOAT* trajMPC_dslbcc00 = trajMPC_ds_cc + 0;
trajMPC_FLOAT* trajMPC_ccrhsl00 = trajMPC_ccrhs + 0;
int trajMPC_ubIdx00[5] = {0, 1, 2, 3, 4};
trajMPC_FLOAT* trajMPC_lub00 = trajMPC_l + 11;
trajMPC_FLOAT* trajMPC_sub00 = trajMPC_s + 11;
trajMPC_FLOAT* trajMPC_lubbysub00 = trajMPC_lbys + 11;
trajMPC_FLOAT trajMPC_riub00[5];
trajMPC_FLOAT* trajMPC_dlubaff00 = trajMPC_dl_aff + 11;
trajMPC_FLOAT* trajMPC_dsubaff00 = trajMPC_ds_aff + 11;
trajMPC_FLOAT* trajMPC_dlubcc00 = trajMPC_dl_cc + 11;
trajMPC_FLOAT* trajMPC_dsubcc00 = trajMPC_ds_cc + 11;
trajMPC_FLOAT* trajMPC_ccrhsub00 = trajMPC_ccrhs + 11;
trajMPC_FLOAT trajMPC_Phi00[11];
trajMPC_FLOAT trajMPC_D00[11] = {1.0000000000000000E+000, 
1.0000000000000000E+000, 
1.0000000000000000E+000};
trajMPC_FLOAT trajMPC_W00[11];
trajMPC_FLOAT trajMPC_H01[5] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 2.0000000000000000E+000, 2.0000000000000000E+000};
trajMPC_FLOAT* trajMPC_z01 = trajMPC_z + 11;
trajMPC_FLOAT* trajMPC_dzaff01 = trajMPC_dz_aff + 11;
trajMPC_FLOAT* trajMPC_dzcc01 = trajMPC_dz_cc + 11;
trajMPC_FLOAT* trajMPC_rd01 = trajMPC_rd + 11;
trajMPC_FLOAT trajMPC_Lbyrd01[11];
trajMPC_FLOAT* trajMPC_grad_cost01 = trajMPC_grad_cost + 11;
trajMPC_FLOAT* trajMPC_grad_eq01 = trajMPC_grad_eq + 11;
trajMPC_FLOAT* trajMPC_grad_ineq01 = trajMPC_grad_ineq + 11;
trajMPC_FLOAT trajMPC_ctv01[11];
trajMPC_FLOAT* trajMPC_v01 = trajMPC_v + 3;
trajMPC_FLOAT trajMPC_re01[3];
trajMPC_FLOAT trajMPC_beta01[3];
trajMPC_FLOAT trajMPC_betacc01[3];
trajMPC_FLOAT* trajMPC_dvaff01 = trajMPC_dv_aff + 3;
trajMPC_FLOAT* trajMPC_dvcc01 = trajMPC_dv_cc + 3;
trajMPC_FLOAT trajMPC_V01[33];
trajMPC_FLOAT trajMPC_Yd01[6];
trajMPC_FLOAT trajMPC_Ld01[6];
trajMPC_FLOAT trajMPC_yy01[3];
trajMPC_FLOAT trajMPC_bmy01[3];
int trajMPC_lbIdx01[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
trajMPC_FLOAT* trajMPC_llb01 = trajMPC_l + 16;
trajMPC_FLOAT* trajMPC_slb01 = trajMPC_s + 16;
trajMPC_FLOAT* trajMPC_llbbyslb01 = trajMPC_lbys + 16;
trajMPC_FLOAT trajMPC_rilb01[11];
trajMPC_FLOAT* trajMPC_dllbaff01 = trajMPC_dl_aff + 16;
trajMPC_FLOAT* trajMPC_dslbaff01 = trajMPC_ds_aff + 16;
trajMPC_FLOAT* trajMPC_dllbcc01 = trajMPC_dl_cc + 16;
trajMPC_FLOAT* trajMPC_dslbcc01 = trajMPC_ds_cc + 16;
trajMPC_FLOAT* trajMPC_ccrhsl01 = trajMPC_ccrhs + 16;
int trajMPC_ubIdx01[5] = {0, 1, 2, 3, 4};
trajMPC_FLOAT* trajMPC_lub01 = trajMPC_l + 27;
trajMPC_FLOAT* trajMPC_sub01 = trajMPC_s + 27;
trajMPC_FLOAT* trajMPC_lubbysub01 = trajMPC_lbys + 27;
trajMPC_FLOAT trajMPC_riub01[5];
trajMPC_FLOAT* trajMPC_dlubaff01 = trajMPC_dl_aff + 27;
trajMPC_FLOAT* trajMPC_dsubaff01 = trajMPC_ds_aff + 27;
trajMPC_FLOAT* trajMPC_dlubcc01 = trajMPC_dl_cc + 27;
trajMPC_FLOAT* trajMPC_dsubcc01 = trajMPC_ds_cc + 27;
trajMPC_FLOAT* trajMPC_ccrhsub01 = trajMPC_ccrhs + 27;
trajMPC_FLOAT trajMPC_Phi01[11];
trajMPC_FLOAT trajMPC_D01[11] = {-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000};
trajMPC_FLOAT trajMPC_W01[11];
trajMPC_FLOAT trajMPC_Ysd01[9];
trajMPC_FLOAT trajMPC_Lsd01[9];
trajMPC_FLOAT* trajMPC_z02 = trajMPC_z + 22;
trajMPC_FLOAT* trajMPC_dzaff02 = trajMPC_dz_aff + 22;
trajMPC_FLOAT* trajMPC_dzcc02 = trajMPC_dz_cc + 22;
trajMPC_FLOAT* trajMPC_rd02 = trajMPC_rd + 22;
trajMPC_FLOAT trajMPC_Lbyrd02[11];
trajMPC_FLOAT* trajMPC_grad_cost02 = trajMPC_grad_cost + 22;
trajMPC_FLOAT* trajMPC_grad_eq02 = trajMPC_grad_eq + 22;
trajMPC_FLOAT* trajMPC_grad_ineq02 = trajMPC_grad_ineq + 22;
trajMPC_FLOAT trajMPC_ctv02[11];
trajMPC_FLOAT* trajMPC_v02 = trajMPC_v + 6;
trajMPC_FLOAT trajMPC_re02[3];
trajMPC_FLOAT trajMPC_beta02[3];
trajMPC_FLOAT trajMPC_betacc02[3];
trajMPC_FLOAT* trajMPC_dvaff02 = trajMPC_dv_aff + 6;
trajMPC_FLOAT* trajMPC_dvcc02 = trajMPC_dv_cc + 6;
trajMPC_FLOAT trajMPC_V02[33];
trajMPC_FLOAT trajMPC_Yd02[6];
trajMPC_FLOAT trajMPC_Ld02[6];
trajMPC_FLOAT trajMPC_yy02[3];
trajMPC_FLOAT trajMPC_bmy02[3];
int trajMPC_lbIdx02[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
trajMPC_FLOAT* trajMPC_llb02 = trajMPC_l + 32;
trajMPC_FLOAT* trajMPC_slb02 = trajMPC_s + 32;
trajMPC_FLOAT* trajMPC_llbbyslb02 = trajMPC_lbys + 32;
trajMPC_FLOAT trajMPC_rilb02[11];
trajMPC_FLOAT* trajMPC_dllbaff02 = trajMPC_dl_aff + 32;
trajMPC_FLOAT* trajMPC_dslbaff02 = trajMPC_ds_aff + 32;
trajMPC_FLOAT* trajMPC_dllbcc02 = trajMPC_dl_cc + 32;
trajMPC_FLOAT* trajMPC_dslbcc02 = trajMPC_ds_cc + 32;
trajMPC_FLOAT* trajMPC_ccrhsl02 = trajMPC_ccrhs + 32;
int trajMPC_ubIdx02[5] = {0, 1, 2, 3, 4};
trajMPC_FLOAT* trajMPC_lub02 = trajMPC_l + 43;
trajMPC_FLOAT* trajMPC_sub02 = trajMPC_s + 43;
trajMPC_FLOAT* trajMPC_lubbysub02 = trajMPC_lbys + 43;
trajMPC_FLOAT trajMPC_riub02[5];
trajMPC_FLOAT* trajMPC_dlubaff02 = trajMPC_dl_aff + 43;
trajMPC_FLOAT* trajMPC_dsubaff02 = trajMPC_ds_aff + 43;
trajMPC_FLOAT* trajMPC_dlubcc02 = trajMPC_dl_cc + 43;
trajMPC_FLOAT* trajMPC_dsubcc02 = trajMPC_ds_cc + 43;
trajMPC_FLOAT* trajMPC_ccrhsub02 = trajMPC_ccrhs + 43;
trajMPC_FLOAT trajMPC_Phi02[11];
trajMPC_FLOAT trajMPC_W02[11];
trajMPC_FLOAT trajMPC_Ysd02[9];
trajMPC_FLOAT trajMPC_Lsd02[9];
trajMPC_FLOAT* trajMPC_z03 = trajMPC_z + 33;
trajMPC_FLOAT* trajMPC_dzaff03 = trajMPC_dz_aff + 33;
trajMPC_FLOAT* trajMPC_dzcc03 = trajMPC_dz_cc + 33;
trajMPC_FLOAT* trajMPC_rd03 = trajMPC_rd + 33;
trajMPC_FLOAT trajMPC_Lbyrd03[11];
trajMPC_FLOAT* trajMPC_grad_cost03 = trajMPC_grad_cost + 33;
trajMPC_FLOAT* trajMPC_grad_eq03 = trajMPC_grad_eq + 33;
trajMPC_FLOAT* trajMPC_grad_ineq03 = trajMPC_grad_ineq + 33;
trajMPC_FLOAT trajMPC_ctv03[11];
trajMPC_FLOAT* trajMPC_v03 = trajMPC_v + 9;
trajMPC_FLOAT trajMPC_re03[3];
trajMPC_FLOAT trajMPC_beta03[3];
trajMPC_FLOAT trajMPC_betacc03[3];
trajMPC_FLOAT* trajMPC_dvaff03 = trajMPC_dv_aff + 9;
trajMPC_FLOAT* trajMPC_dvcc03 = trajMPC_dv_cc + 9;
trajMPC_FLOAT trajMPC_V03[33];
trajMPC_FLOAT trajMPC_Yd03[6];
trajMPC_FLOAT trajMPC_Ld03[6];
trajMPC_FLOAT trajMPC_yy03[3];
trajMPC_FLOAT trajMPC_bmy03[3];
int trajMPC_lbIdx03[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
trajMPC_FLOAT* trajMPC_llb03 = trajMPC_l + 48;
trajMPC_FLOAT* trajMPC_slb03 = trajMPC_s + 48;
trajMPC_FLOAT* trajMPC_llbbyslb03 = trajMPC_lbys + 48;
trajMPC_FLOAT trajMPC_rilb03[11];
trajMPC_FLOAT* trajMPC_dllbaff03 = trajMPC_dl_aff + 48;
trajMPC_FLOAT* trajMPC_dslbaff03 = trajMPC_ds_aff + 48;
trajMPC_FLOAT* trajMPC_dllbcc03 = trajMPC_dl_cc + 48;
trajMPC_FLOAT* trajMPC_dslbcc03 = trajMPC_ds_cc + 48;
trajMPC_FLOAT* trajMPC_ccrhsl03 = trajMPC_ccrhs + 48;
int trajMPC_ubIdx03[5] = {0, 1, 2, 3, 4};
trajMPC_FLOAT* trajMPC_lub03 = trajMPC_l + 59;
trajMPC_FLOAT* trajMPC_sub03 = trajMPC_s + 59;
trajMPC_FLOAT* trajMPC_lubbysub03 = trajMPC_lbys + 59;
trajMPC_FLOAT trajMPC_riub03[5];
trajMPC_FLOAT* trajMPC_dlubaff03 = trajMPC_dl_aff + 59;
trajMPC_FLOAT* trajMPC_dsubaff03 = trajMPC_ds_aff + 59;
trajMPC_FLOAT* trajMPC_dlubcc03 = trajMPC_dl_cc + 59;
trajMPC_FLOAT* trajMPC_dsubcc03 = trajMPC_ds_cc + 59;
trajMPC_FLOAT* trajMPC_ccrhsub03 = trajMPC_ccrhs + 59;
trajMPC_FLOAT trajMPC_Phi03[11];
trajMPC_FLOAT trajMPC_W03[11];
trajMPC_FLOAT trajMPC_Ysd03[9];
trajMPC_FLOAT trajMPC_Lsd03[9];
trajMPC_FLOAT* trajMPC_z04 = trajMPC_z + 44;
trajMPC_FLOAT* trajMPC_dzaff04 = trajMPC_dz_aff + 44;
trajMPC_FLOAT* trajMPC_dzcc04 = trajMPC_dz_cc + 44;
trajMPC_FLOAT* trajMPC_rd04 = trajMPC_rd + 44;
trajMPC_FLOAT trajMPC_Lbyrd04[11];
trajMPC_FLOAT* trajMPC_grad_cost04 = trajMPC_grad_cost + 44;
trajMPC_FLOAT* trajMPC_grad_eq04 = trajMPC_grad_eq + 44;
trajMPC_FLOAT* trajMPC_grad_ineq04 = trajMPC_grad_ineq + 44;
trajMPC_FLOAT trajMPC_ctv04[11];
trajMPC_FLOAT* trajMPC_v04 = trajMPC_v + 12;
trajMPC_FLOAT trajMPC_re04[3];
trajMPC_FLOAT trajMPC_beta04[3];
trajMPC_FLOAT trajMPC_betacc04[3];
trajMPC_FLOAT* trajMPC_dvaff04 = trajMPC_dv_aff + 12;
trajMPC_FLOAT* trajMPC_dvcc04 = trajMPC_dv_cc + 12;
trajMPC_FLOAT trajMPC_V04[33];
trajMPC_FLOAT trajMPC_Yd04[6];
trajMPC_FLOAT trajMPC_Ld04[6];
trajMPC_FLOAT trajMPC_yy04[3];
trajMPC_FLOAT trajMPC_bmy04[3];
int trajMPC_lbIdx04[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
trajMPC_FLOAT* trajMPC_llb04 = trajMPC_l + 64;
trajMPC_FLOAT* trajMPC_slb04 = trajMPC_s + 64;
trajMPC_FLOAT* trajMPC_llbbyslb04 = trajMPC_lbys + 64;
trajMPC_FLOAT trajMPC_rilb04[11];
trajMPC_FLOAT* trajMPC_dllbaff04 = trajMPC_dl_aff + 64;
trajMPC_FLOAT* trajMPC_dslbaff04 = trajMPC_ds_aff + 64;
trajMPC_FLOAT* trajMPC_dllbcc04 = trajMPC_dl_cc + 64;
trajMPC_FLOAT* trajMPC_dslbcc04 = trajMPC_ds_cc + 64;
trajMPC_FLOAT* trajMPC_ccrhsl04 = trajMPC_ccrhs + 64;
int trajMPC_ubIdx04[5] = {0, 1, 2, 3, 4};
trajMPC_FLOAT* trajMPC_lub04 = trajMPC_l + 75;
trajMPC_FLOAT* trajMPC_sub04 = trajMPC_s + 75;
trajMPC_FLOAT* trajMPC_lubbysub04 = trajMPC_lbys + 75;
trajMPC_FLOAT trajMPC_riub04[5];
trajMPC_FLOAT* trajMPC_dlubaff04 = trajMPC_dl_aff + 75;
trajMPC_FLOAT* trajMPC_dsubaff04 = trajMPC_ds_aff + 75;
trajMPC_FLOAT* trajMPC_dlubcc04 = trajMPC_dl_cc + 75;
trajMPC_FLOAT* trajMPC_dsubcc04 = trajMPC_ds_cc + 75;
trajMPC_FLOAT* trajMPC_ccrhsub04 = trajMPC_ccrhs + 75;
trajMPC_FLOAT trajMPC_Phi04[11];
trajMPC_FLOAT trajMPC_W04[11];
trajMPC_FLOAT trajMPC_Ysd04[9];
trajMPC_FLOAT trajMPC_Lsd04[9];
trajMPC_FLOAT* trajMPC_z05 = trajMPC_z + 55;
trajMPC_FLOAT* trajMPC_dzaff05 = trajMPC_dz_aff + 55;
trajMPC_FLOAT* trajMPC_dzcc05 = trajMPC_dz_cc + 55;
trajMPC_FLOAT* trajMPC_rd05 = trajMPC_rd + 55;
trajMPC_FLOAT trajMPC_Lbyrd05[11];
trajMPC_FLOAT* trajMPC_grad_cost05 = trajMPC_grad_cost + 55;
trajMPC_FLOAT* trajMPC_grad_eq05 = trajMPC_grad_eq + 55;
trajMPC_FLOAT* trajMPC_grad_ineq05 = trajMPC_grad_ineq + 55;
trajMPC_FLOAT trajMPC_ctv05[11];
trajMPC_FLOAT* trajMPC_v05 = trajMPC_v + 15;
trajMPC_FLOAT trajMPC_re05[3];
trajMPC_FLOAT trajMPC_beta05[3];
trajMPC_FLOAT trajMPC_betacc05[3];
trajMPC_FLOAT* trajMPC_dvaff05 = trajMPC_dv_aff + 15;
trajMPC_FLOAT* trajMPC_dvcc05 = trajMPC_dv_cc + 15;
trajMPC_FLOAT trajMPC_V05[33];
trajMPC_FLOAT trajMPC_Yd05[6];
trajMPC_FLOAT trajMPC_Ld05[6];
trajMPC_FLOAT trajMPC_yy05[3];
trajMPC_FLOAT trajMPC_bmy05[3];
int trajMPC_lbIdx05[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
trajMPC_FLOAT* trajMPC_llb05 = trajMPC_l + 80;
trajMPC_FLOAT* trajMPC_slb05 = trajMPC_s + 80;
trajMPC_FLOAT* trajMPC_llbbyslb05 = trajMPC_lbys + 80;
trajMPC_FLOAT trajMPC_rilb05[11];
trajMPC_FLOAT* trajMPC_dllbaff05 = trajMPC_dl_aff + 80;
trajMPC_FLOAT* trajMPC_dslbaff05 = trajMPC_ds_aff + 80;
trajMPC_FLOAT* trajMPC_dllbcc05 = trajMPC_dl_cc + 80;
trajMPC_FLOAT* trajMPC_dslbcc05 = trajMPC_ds_cc + 80;
trajMPC_FLOAT* trajMPC_ccrhsl05 = trajMPC_ccrhs + 80;
int trajMPC_ubIdx05[5] = {0, 1, 2, 3, 4};
trajMPC_FLOAT* trajMPC_lub05 = trajMPC_l + 91;
trajMPC_FLOAT* trajMPC_sub05 = trajMPC_s + 91;
trajMPC_FLOAT* trajMPC_lubbysub05 = trajMPC_lbys + 91;
trajMPC_FLOAT trajMPC_riub05[5];
trajMPC_FLOAT* trajMPC_dlubaff05 = trajMPC_dl_aff + 91;
trajMPC_FLOAT* trajMPC_dsubaff05 = trajMPC_ds_aff + 91;
trajMPC_FLOAT* trajMPC_dlubcc05 = trajMPC_dl_cc + 91;
trajMPC_FLOAT* trajMPC_dsubcc05 = trajMPC_ds_cc + 91;
trajMPC_FLOAT* trajMPC_ccrhsub05 = trajMPC_ccrhs + 91;
trajMPC_FLOAT trajMPC_Phi05[11];
trajMPC_FLOAT trajMPC_W05[11];
trajMPC_FLOAT trajMPC_Ysd05[9];
trajMPC_FLOAT trajMPC_Lsd05[9];
trajMPC_FLOAT* trajMPC_z06 = trajMPC_z + 66;
trajMPC_FLOAT* trajMPC_dzaff06 = trajMPC_dz_aff + 66;
trajMPC_FLOAT* trajMPC_dzcc06 = trajMPC_dz_cc + 66;
trajMPC_FLOAT* trajMPC_rd06 = trajMPC_rd + 66;
trajMPC_FLOAT trajMPC_Lbyrd06[11];
trajMPC_FLOAT* trajMPC_grad_cost06 = trajMPC_grad_cost + 66;
trajMPC_FLOAT* trajMPC_grad_eq06 = trajMPC_grad_eq + 66;
trajMPC_FLOAT* trajMPC_grad_ineq06 = trajMPC_grad_ineq + 66;
trajMPC_FLOAT trajMPC_ctv06[11];
trajMPC_FLOAT* trajMPC_v06 = trajMPC_v + 18;
trajMPC_FLOAT trajMPC_re06[3];
trajMPC_FLOAT trajMPC_beta06[3];
trajMPC_FLOAT trajMPC_betacc06[3];
trajMPC_FLOAT* trajMPC_dvaff06 = trajMPC_dv_aff + 18;
trajMPC_FLOAT* trajMPC_dvcc06 = trajMPC_dv_cc + 18;
trajMPC_FLOAT trajMPC_V06[33];
trajMPC_FLOAT trajMPC_Yd06[6];
trajMPC_FLOAT trajMPC_Ld06[6];
trajMPC_FLOAT trajMPC_yy06[3];
trajMPC_FLOAT trajMPC_bmy06[3];
int trajMPC_lbIdx06[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
trajMPC_FLOAT* trajMPC_llb06 = trajMPC_l + 96;
trajMPC_FLOAT* trajMPC_slb06 = trajMPC_s + 96;
trajMPC_FLOAT* trajMPC_llbbyslb06 = trajMPC_lbys + 96;
trajMPC_FLOAT trajMPC_rilb06[11];
trajMPC_FLOAT* trajMPC_dllbaff06 = trajMPC_dl_aff + 96;
trajMPC_FLOAT* trajMPC_dslbaff06 = trajMPC_ds_aff + 96;
trajMPC_FLOAT* trajMPC_dllbcc06 = trajMPC_dl_cc + 96;
trajMPC_FLOAT* trajMPC_dslbcc06 = trajMPC_ds_cc + 96;
trajMPC_FLOAT* trajMPC_ccrhsl06 = trajMPC_ccrhs + 96;
int trajMPC_ubIdx06[5] = {0, 1, 2, 3, 4};
trajMPC_FLOAT* trajMPC_lub06 = trajMPC_l + 107;
trajMPC_FLOAT* trajMPC_sub06 = trajMPC_s + 107;
trajMPC_FLOAT* trajMPC_lubbysub06 = trajMPC_lbys + 107;
trajMPC_FLOAT trajMPC_riub06[5];
trajMPC_FLOAT* trajMPC_dlubaff06 = trajMPC_dl_aff + 107;
trajMPC_FLOAT* trajMPC_dsubaff06 = trajMPC_ds_aff + 107;
trajMPC_FLOAT* trajMPC_dlubcc06 = trajMPC_dl_cc + 107;
trajMPC_FLOAT* trajMPC_dsubcc06 = trajMPC_ds_cc + 107;
trajMPC_FLOAT* trajMPC_ccrhsub06 = trajMPC_ccrhs + 107;
trajMPC_FLOAT trajMPC_Phi06[11];
trajMPC_FLOAT trajMPC_W06[11];
trajMPC_FLOAT trajMPC_Ysd06[9];
trajMPC_FLOAT trajMPC_Lsd06[9];
trajMPC_FLOAT* trajMPC_z07 = trajMPC_z + 77;
trajMPC_FLOAT* trajMPC_dzaff07 = trajMPC_dz_aff + 77;
trajMPC_FLOAT* trajMPC_dzcc07 = trajMPC_dz_cc + 77;
trajMPC_FLOAT* trajMPC_rd07 = trajMPC_rd + 77;
trajMPC_FLOAT trajMPC_Lbyrd07[11];
trajMPC_FLOAT* trajMPC_grad_cost07 = trajMPC_grad_cost + 77;
trajMPC_FLOAT* trajMPC_grad_eq07 = trajMPC_grad_eq + 77;
trajMPC_FLOAT* trajMPC_grad_ineq07 = trajMPC_grad_ineq + 77;
trajMPC_FLOAT trajMPC_ctv07[11];
trajMPC_FLOAT* trajMPC_v07 = trajMPC_v + 21;
trajMPC_FLOAT trajMPC_re07[3];
trajMPC_FLOAT trajMPC_beta07[3];
trajMPC_FLOAT trajMPC_betacc07[3];
trajMPC_FLOAT* trajMPC_dvaff07 = trajMPC_dv_aff + 21;
trajMPC_FLOAT* trajMPC_dvcc07 = trajMPC_dv_cc + 21;
trajMPC_FLOAT trajMPC_V07[33];
trajMPC_FLOAT trajMPC_Yd07[6];
trajMPC_FLOAT trajMPC_Ld07[6];
trajMPC_FLOAT trajMPC_yy07[3];
trajMPC_FLOAT trajMPC_bmy07[3];
int trajMPC_lbIdx07[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
trajMPC_FLOAT* trajMPC_llb07 = trajMPC_l + 112;
trajMPC_FLOAT* trajMPC_slb07 = trajMPC_s + 112;
trajMPC_FLOAT* trajMPC_llbbyslb07 = trajMPC_lbys + 112;
trajMPC_FLOAT trajMPC_rilb07[11];
trajMPC_FLOAT* trajMPC_dllbaff07 = trajMPC_dl_aff + 112;
trajMPC_FLOAT* trajMPC_dslbaff07 = trajMPC_ds_aff + 112;
trajMPC_FLOAT* trajMPC_dllbcc07 = trajMPC_dl_cc + 112;
trajMPC_FLOAT* trajMPC_dslbcc07 = trajMPC_ds_cc + 112;
trajMPC_FLOAT* trajMPC_ccrhsl07 = trajMPC_ccrhs + 112;
int trajMPC_ubIdx07[5] = {0, 1, 2, 3, 4};
trajMPC_FLOAT* trajMPC_lub07 = trajMPC_l + 123;
trajMPC_FLOAT* trajMPC_sub07 = trajMPC_s + 123;
trajMPC_FLOAT* trajMPC_lubbysub07 = trajMPC_lbys + 123;
trajMPC_FLOAT trajMPC_riub07[5];
trajMPC_FLOAT* trajMPC_dlubaff07 = trajMPC_dl_aff + 123;
trajMPC_FLOAT* trajMPC_dsubaff07 = trajMPC_ds_aff + 123;
trajMPC_FLOAT* trajMPC_dlubcc07 = trajMPC_dl_cc + 123;
trajMPC_FLOAT* trajMPC_dsubcc07 = trajMPC_ds_cc + 123;
trajMPC_FLOAT* trajMPC_ccrhsub07 = trajMPC_ccrhs + 123;
trajMPC_FLOAT trajMPC_Phi07[11];
trajMPC_FLOAT trajMPC_W07[11];
trajMPC_FLOAT trajMPC_Ysd07[9];
trajMPC_FLOAT trajMPC_Lsd07[9];
trajMPC_FLOAT* trajMPC_z08 = trajMPC_z + 88;
trajMPC_FLOAT* trajMPC_dzaff08 = trajMPC_dz_aff + 88;
trajMPC_FLOAT* trajMPC_dzcc08 = trajMPC_dz_cc + 88;
trajMPC_FLOAT* trajMPC_rd08 = trajMPC_rd + 88;
trajMPC_FLOAT trajMPC_Lbyrd08[11];
trajMPC_FLOAT* trajMPC_grad_cost08 = trajMPC_grad_cost + 88;
trajMPC_FLOAT* trajMPC_grad_eq08 = trajMPC_grad_eq + 88;
trajMPC_FLOAT* trajMPC_grad_ineq08 = trajMPC_grad_ineq + 88;
trajMPC_FLOAT trajMPC_ctv08[11];
trajMPC_FLOAT* trajMPC_v08 = trajMPC_v + 24;
trajMPC_FLOAT trajMPC_re08[3];
trajMPC_FLOAT trajMPC_beta08[3];
trajMPC_FLOAT trajMPC_betacc08[3];
trajMPC_FLOAT* trajMPC_dvaff08 = trajMPC_dv_aff + 24;
trajMPC_FLOAT* trajMPC_dvcc08 = trajMPC_dv_cc + 24;
trajMPC_FLOAT trajMPC_V08[33];
trajMPC_FLOAT trajMPC_Yd08[6];
trajMPC_FLOAT trajMPC_Ld08[6];
trajMPC_FLOAT trajMPC_yy08[3];
trajMPC_FLOAT trajMPC_bmy08[3];
int trajMPC_lbIdx08[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
trajMPC_FLOAT* trajMPC_llb08 = trajMPC_l + 128;
trajMPC_FLOAT* trajMPC_slb08 = trajMPC_s + 128;
trajMPC_FLOAT* trajMPC_llbbyslb08 = trajMPC_lbys + 128;
trajMPC_FLOAT trajMPC_rilb08[11];
trajMPC_FLOAT* trajMPC_dllbaff08 = trajMPC_dl_aff + 128;
trajMPC_FLOAT* trajMPC_dslbaff08 = trajMPC_ds_aff + 128;
trajMPC_FLOAT* trajMPC_dllbcc08 = trajMPC_dl_cc + 128;
trajMPC_FLOAT* trajMPC_dslbcc08 = trajMPC_ds_cc + 128;
trajMPC_FLOAT* trajMPC_ccrhsl08 = trajMPC_ccrhs + 128;
int trajMPC_ubIdx08[5] = {0, 1, 2, 3, 4};
trajMPC_FLOAT* trajMPC_lub08 = trajMPC_l + 139;
trajMPC_FLOAT* trajMPC_sub08 = trajMPC_s + 139;
trajMPC_FLOAT* trajMPC_lubbysub08 = trajMPC_lbys + 139;
trajMPC_FLOAT trajMPC_riub08[5];
trajMPC_FLOAT* trajMPC_dlubaff08 = trajMPC_dl_aff + 139;
trajMPC_FLOAT* trajMPC_dsubaff08 = trajMPC_ds_aff + 139;
trajMPC_FLOAT* trajMPC_dlubcc08 = trajMPC_dl_cc + 139;
trajMPC_FLOAT* trajMPC_dsubcc08 = trajMPC_ds_cc + 139;
trajMPC_FLOAT* trajMPC_ccrhsub08 = trajMPC_ccrhs + 139;
trajMPC_FLOAT trajMPC_Phi08[11];
trajMPC_FLOAT trajMPC_W08[11];
trajMPC_FLOAT trajMPC_Ysd08[9];
trajMPC_FLOAT trajMPC_Lsd08[9];
trajMPC_FLOAT* trajMPC_z09 = trajMPC_z + 99;
trajMPC_FLOAT* trajMPC_dzaff09 = trajMPC_dz_aff + 99;
trajMPC_FLOAT* trajMPC_dzcc09 = trajMPC_dz_cc + 99;
trajMPC_FLOAT* trajMPC_rd09 = trajMPC_rd + 99;
trajMPC_FLOAT trajMPC_Lbyrd09[11];
trajMPC_FLOAT* trajMPC_grad_cost09 = trajMPC_grad_cost + 99;
trajMPC_FLOAT* trajMPC_grad_eq09 = trajMPC_grad_eq + 99;
trajMPC_FLOAT* trajMPC_grad_ineq09 = trajMPC_grad_ineq + 99;
trajMPC_FLOAT trajMPC_ctv09[11];
trajMPC_FLOAT* trajMPC_v09 = trajMPC_v + 27;
trajMPC_FLOAT trajMPC_re09[3];
trajMPC_FLOAT trajMPC_beta09[3];
trajMPC_FLOAT trajMPC_betacc09[3];
trajMPC_FLOAT* trajMPC_dvaff09 = trajMPC_dv_aff + 27;
trajMPC_FLOAT* trajMPC_dvcc09 = trajMPC_dv_cc + 27;
trajMPC_FLOAT trajMPC_V09[33];
trajMPC_FLOAT trajMPC_Yd09[6];
trajMPC_FLOAT trajMPC_Ld09[6];
trajMPC_FLOAT trajMPC_yy09[3];
trajMPC_FLOAT trajMPC_bmy09[3];
int trajMPC_lbIdx09[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
trajMPC_FLOAT* trajMPC_llb09 = trajMPC_l + 144;
trajMPC_FLOAT* trajMPC_slb09 = trajMPC_s + 144;
trajMPC_FLOAT* trajMPC_llbbyslb09 = trajMPC_lbys + 144;
trajMPC_FLOAT trajMPC_rilb09[11];
trajMPC_FLOAT* trajMPC_dllbaff09 = trajMPC_dl_aff + 144;
trajMPC_FLOAT* trajMPC_dslbaff09 = trajMPC_ds_aff + 144;
trajMPC_FLOAT* trajMPC_dllbcc09 = trajMPC_dl_cc + 144;
trajMPC_FLOAT* trajMPC_dslbcc09 = trajMPC_ds_cc + 144;
trajMPC_FLOAT* trajMPC_ccrhsl09 = trajMPC_ccrhs + 144;
int trajMPC_ubIdx09[5] = {0, 1, 2, 3, 4};
trajMPC_FLOAT* trajMPC_lub09 = trajMPC_l + 155;
trajMPC_FLOAT* trajMPC_sub09 = trajMPC_s + 155;
trajMPC_FLOAT* trajMPC_lubbysub09 = trajMPC_lbys + 155;
trajMPC_FLOAT trajMPC_riub09[5];
trajMPC_FLOAT* trajMPC_dlubaff09 = trajMPC_dl_aff + 155;
trajMPC_FLOAT* trajMPC_dsubaff09 = trajMPC_ds_aff + 155;
trajMPC_FLOAT* trajMPC_dlubcc09 = trajMPC_dl_cc + 155;
trajMPC_FLOAT* trajMPC_dsubcc09 = trajMPC_ds_cc + 155;
trajMPC_FLOAT* trajMPC_ccrhsub09 = trajMPC_ccrhs + 155;
trajMPC_FLOAT trajMPC_Phi09[11];
trajMPC_FLOAT trajMPC_W09[11];
trajMPC_FLOAT trajMPC_Ysd09[9];
trajMPC_FLOAT trajMPC_Lsd09[9];
trajMPC_FLOAT* trajMPC_z10 = trajMPC_z + 110;
trajMPC_FLOAT* trajMPC_dzaff10 = trajMPC_dz_aff + 110;
trajMPC_FLOAT* trajMPC_dzcc10 = trajMPC_dz_cc + 110;
trajMPC_FLOAT* trajMPC_rd10 = trajMPC_rd + 110;
trajMPC_FLOAT trajMPC_Lbyrd10[11];
trajMPC_FLOAT* trajMPC_grad_cost10 = trajMPC_grad_cost + 110;
trajMPC_FLOAT* trajMPC_grad_eq10 = trajMPC_grad_eq + 110;
trajMPC_FLOAT* trajMPC_grad_ineq10 = trajMPC_grad_ineq + 110;
trajMPC_FLOAT trajMPC_ctv10[11];
trajMPC_FLOAT* trajMPC_v10 = trajMPC_v + 30;
trajMPC_FLOAT trajMPC_re10[3];
trajMPC_FLOAT trajMPC_beta10[3];
trajMPC_FLOAT trajMPC_betacc10[3];
trajMPC_FLOAT* trajMPC_dvaff10 = trajMPC_dv_aff + 30;
trajMPC_FLOAT* trajMPC_dvcc10 = trajMPC_dv_cc + 30;
trajMPC_FLOAT trajMPC_V10[33];
trajMPC_FLOAT trajMPC_Yd10[6];
trajMPC_FLOAT trajMPC_Ld10[6];
trajMPC_FLOAT trajMPC_yy10[3];
trajMPC_FLOAT trajMPC_bmy10[3];
int trajMPC_lbIdx10[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
trajMPC_FLOAT* trajMPC_llb10 = trajMPC_l + 160;
trajMPC_FLOAT* trajMPC_slb10 = trajMPC_s + 160;
trajMPC_FLOAT* trajMPC_llbbyslb10 = trajMPC_lbys + 160;
trajMPC_FLOAT trajMPC_rilb10[11];
trajMPC_FLOAT* trajMPC_dllbaff10 = trajMPC_dl_aff + 160;
trajMPC_FLOAT* trajMPC_dslbaff10 = trajMPC_ds_aff + 160;
trajMPC_FLOAT* trajMPC_dllbcc10 = trajMPC_dl_cc + 160;
trajMPC_FLOAT* trajMPC_dslbcc10 = trajMPC_ds_cc + 160;
trajMPC_FLOAT* trajMPC_ccrhsl10 = trajMPC_ccrhs + 160;
int trajMPC_ubIdx10[5] = {0, 1, 2, 3, 4};
trajMPC_FLOAT* trajMPC_lub10 = trajMPC_l + 171;
trajMPC_FLOAT* trajMPC_sub10 = trajMPC_s + 171;
trajMPC_FLOAT* trajMPC_lubbysub10 = trajMPC_lbys + 171;
trajMPC_FLOAT trajMPC_riub10[5];
trajMPC_FLOAT* trajMPC_dlubaff10 = trajMPC_dl_aff + 171;
trajMPC_FLOAT* trajMPC_dsubaff10 = trajMPC_ds_aff + 171;
trajMPC_FLOAT* trajMPC_dlubcc10 = trajMPC_dl_cc + 171;
trajMPC_FLOAT* trajMPC_dsubcc10 = trajMPC_ds_cc + 171;
trajMPC_FLOAT* trajMPC_ccrhsub10 = trajMPC_ccrhs + 171;
trajMPC_FLOAT trajMPC_Phi10[11];
trajMPC_FLOAT trajMPC_W10[11];
trajMPC_FLOAT trajMPC_Ysd10[9];
trajMPC_FLOAT trajMPC_Lsd10[9];
trajMPC_FLOAT* trajMPC_z11 = trajMPC_z + 121;
trajMPC_FLOAT* trajMPC_dzaff11 = trajMPC_dz_aff + 121;
trajMPC_FLOAT* trajMPC_dzcc11 = trajMPC_dz_cc + 121;
trajMPC_FLOAT* trajMPC_rd11 = trajMPC_rd + 121;
trajMPC_FLOAT trajMPC_Lbyrd11[11];
trajMPC_FLOAT* trajMPC_grad_cost11 = trajMPC_grad_cost + 121;
trajMPC_FLOAT* trajMPC_grad_eq11 = trajMPC_grad_eq + 121;
trajMPC_FLOAT* trajMPC_grad_ineq11 = trajMPC_grad_ineq + 121;
trajMPC_FLOAT trajMPC_ctv11[11];
trajMPC_FLOAT* trajMPC_v11 = trajMPC_v + 33;
trajMPC_FLOAT trajMPC_re11[3];
trajMPC_FLOAT trajMPC_beta11[3];
trajMPC_FLOAT trajMPC_betacc11[3];
trajMPC_FLOAT* trajMPC_dvaff11 = trajMPC_dv_aff + 33;
trajMPC_FLOAT* trajMPC_dvcc11 = trajMPC_dv_cc + 33;
trajMPC_FLOAT trajMPC_V11[33];
trajMPC_FLOAT trajMPC_Yd11[6];
trajMPC_FLOAT trajMPC_Ld11[6];
trajMPC_FLOAT trajMPC_yy11[3];
trajMPC_FLOAT trajMPC_bmy11[3];
int trajMPC_lbIdx11[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
trajMPC_FLOAT* trajMPC_llb11 = trajMPC_l + 176;
trajMPC_FLOAT* trajMPC_slb11 = trajMPC_s + 176;
trajMPC_FLOAT* trajMPC_llbbyslb11 = trajMPC_lbys + 176;
trajMPC_FLOAT trajMPC_rilb11[11];
trajMPC_FLOAT* trajMPC_dllbaff11 = trajMPC_dl_aff + 176;
trajMPC_FLOAT* trajMPC_dslbaff11 = trajMPC_ds_aff + 176;
trajMPC_FLOAT* trajMPC_dllbcc11 = trajMPC_dl_cc + 176;
trajMPC_FLOAT* trajMPC_dslbcc11 = trajMPC_ds_cc + 176;
trajMPC_FLOAT* trajMPC_ccrhsl11 = trajMPC_ccrhs + 176;
int trajMPC_ubIdx11[5] = {0, 1, 2, 3, 4};
trajMPC_FLOAT* trajMPC_lub11 = trajMPC_l + 187;
trajMPC_FLOAT* trajMPC_sub11 = trajMPC_s + 187;
trajMPC_FLOAT* trajMPC_lubbysub11 = trajMPC_lbys + 187;
trajMPC_FLOAT trajMPC_riub11[5];
trajMPC_FLOAT* trajMPC_dlubaff11 = trajMPC_dl_aff + 187;
trajMPC_FLOAT* trajMPC_dsubaff11 = trajMPC_ds_aff + 187;
trajMPC_FLOAT* trajMPC_dlubcc11 = trajMPC_dl_cc + 187;
trajMPC_FLOAT* trajMPC_dsubcc11 = trajMPC_ds_cc + 187;
trajMPC_FLOAT* trajMPC_ccrhsub11 = trajMPC_ccrhs + 187;
trajMPC_FLOAT trajMPC_Phi11[11];
trajMPC_FLOAT trajMPC_W11[11];
trajMPC_FLOAT trajMPC_Ysd11[9];
trajMPC_FLOAT trajMPC_Lsd11[9];
trajMPC_FLOAT* trajMPC_z12 = trajMPC_z + 132;
trajMPC_FLOAT* trajMPC_dzaff12 = trajMPC_dz_aff + 132;
trajMPC_FLOAT* trajMPC_dzcc12 = trajMPC_dz_cc + 132;
trajMPC_FLOAT* trajMPC_rd12 = trajMPC_rd + 132;
trajMPC_FLOAT trajMPC_Lbyrd12[11];
trajMPC_FLOAT* trajMPC_grad_cost12 = trajMPC_grad_cost + 132;
trajMPC_FLOAT* trajMPC_grad_eq12 = trajMPC_grad_eq + 132;
trajMPC_FLOAT* trajMPC_grad_ineq12 = trajMPC_grad_ineq + 132;
trajMPC_FLOAT trajMPC_ctv12[11];
trajMPC_FLOAT* trajMPC_v12 = trajMPC_v + 36;
trajMPC_FLOAT trajMPC_re12[3];
trajMPC_FLOAT trajMPC_beta12[3];
trajMPC_FLOAT trajMPC_betacc12[3];
trajMPC_FLOAT* trajMPC_dvaff12 = trajMPC_dv_aff + 36;
trajMPC_FLOAT* trajMPC_dvcc12 = trajMPC_dv_cc + 36;
trajMPC_FLOAT trajMPC_V12[33];
trajMPC_FLOAT trajMPC_Yd12[6];
trajMPC_FLOAT trajMPC_Ld12[6];
trajMPC_FLOAT trajMPC_yy12[3];
trajMPC_FLOAT trajMPC_bmy12[3];
int trajMPC_lbIdx12[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
trajMPC_FLOAT* trajMPC_llb12 = trajMPC_l + 192;
trajMPC_FLOAT* trajMPC_slb12 = trajMPC_s + 192;
trajMPC_FLOAT* trajMPC_llbbyslb12 = trajMPC_lbys + 192;
trajMPC_FLOAT trajMPC_rilb12[11];
trajMPC_FLOAT* trajMPC_dllbaff12 = trajMPC_dl_aff + 192;
trajMPC_FLOAT* trajMPC_dslbaff12 = trajMPC_ds_aff + 192;
trajMPC_FLOAT* trajMPC_dllbcc12 = trajMPC_dl_cc + 192;
trajMPC_FLOAT* trajMPC_dslbcc12 = trajMPC_ds_cc + 192;
trajMPC_FLOAT* trajMPC_ccrhsl12 = trajMPC_ccrhs + 192;
int trajMPC_ubIdx12[5] = {0, 1, 2, 3, 4};
trajMPC_FLOAT* trajMPC_lub12 = trajMPC_l + 203;
trajMPC_FLOAT* trajMPC_sub12 = trajMPC_s + 203;
trajMPC_FLOAT* trajMPC_lubbysub12 = trajMPC_lbys + 203;
trajMPC_FLOAT trajMPC_riub12[5];
trajMPC_FLOAT* trajMPC_dlubaff12 = trajMPC_dl_aff + 203;
trajMPC_FLOAT* trajMPC_dsubaff12 = trajMPC_ds_aff + 203;
trajMPC_FLOAT* trajMPC_dlubcc12 = trajMPC_dl_cc + 203;
trajMPC_FLOAT* trajMPC_dsubcc12 = trajMPC_ds_cc + 203;
trajMPC_FLOAT* trajMPC_ccrhsub12 = trajMPC_ccrhs + 203;
trajMPC_FLOAT trajMPC_Phi12[11];
trajMPC_FLOAT trajMPC_W12[11];
trajMPC_FLOAT trajMPC_Ysd12[9];
trajMPC_FLOAT trajMPC_Lsd12[9];
trajMPC_FLOAT* trajMPC_z13 = trajMPC_z + 143;
trajMPC_FLOAT* trajMPC_dzaff13 = trajMPC_dz_aff + 143;
trajMPC_FLOAT* trajMPC_dzcc13 = trajMPC_dz_cc + 143;
trajMPC_FLOAT* trajMPC_rd13 = trajMPC_rd + 143;
trajMPC_FLOAT trajMPC_Lbyrd13[11];
trajMPC_FLOAT* trajMPC_grad_cost13 = trajMPC_grad_cost + 143;
trajMPC_FLOAT* trajMPC_grad_eq13 = trajMPC_grad_eq + 143;
trajMPC_FLOAT* trajMPC_grad_ineq13 = trajMPC_grad_ineq + 143;
trajMPC_FLOAT trajMPC_ctv13[11];
trajMPC_FLOAT* trajMPC_v13 = trajMPC_v + 39;
trajMPC_FLOAT trajMPC_re13[3];
trajMPC_FLOAT trajMPC_beta13[3];
trajMPC_FLOAT trajMPC_betacc13[3];
trajMPC_FLOAT* trajMPC_dvaff13 = trajMPC_dv_aff + 39;
trajMPC_FLOAT* trajMPC_dvcc13 = trajMPC_dv_cc + 39;
trajMPC_FLOAT trajMPC_V13[33];
trajMPC_FLOAT trajMPC_Yd13[6];
trajMPC_FLOAT trajMPC_Ld13[6];
trajMPC_FLOAT trajMPC_yy13[3];
trajMPC_FLOAT trajMPC_bmy13[3];
int trajMPC_lbIdx13[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
trajMPC_FLOAT* trajMPC_llb13 = trajMPC_l + 208;
trajMPC_FLOAT* trajMPC_slb13 = trajMPC_s + 208;
trajMPC_FLOAT* trajMPC_llbbyslb13 = trajMPC_lbys + 208;
trajMPC_FLOAT trajMPC_rilb13[11];
trajMPC_FLOAT* trajMPC_dllbaff13 = trajMPC_dl_aff + 208;
trajMPC_FLOAT* trajMPC_dslbaff13 = trajMPC_ds_aff + 208;
trajMPC_FLOAT* trajMPC_dllbcc13 = trajMPC_dl_cc + 208;
trajMPC_FLOAT* trajMPC_dslbcc13 = trajMPC_ds_cc + 208;
trajMPC_FLOAT* trajMPC_ccrhsl13 = trajMPC_ccrhs + 208;
int trajMPC_ubIdx13[5] = {0, 1, 2, 3, 4};
trajMPC_FLOAT* trajMPC_lub13 = trajMPC_l + 219;
trajMPC_FLOAT* trajMPC_sub13 = trajMPC_s + 219;
trajMPC_FLOAT* trajMPC_lubbysub13 = trajMPC_lbys + 219;
trajMPC_FLOAT trajMPC_riub13[5];
trajMPC_FLOAT* trajMPC_dlubaff13 = trajMPC_dl_aff + 219;
trajMPC_FLOAT* trajMPC_dsubaff13 = trajMPC_ds_aff + 219;
trajMPC_FLOAT* trajMPC_dlubcc13 = trajMPC_dl_cc + 219;
trajMPC_FLOAT* trajMPC_dsubcc13 = trajMPC_ds_cc + 219;
trajMPC_FLOAT* trajMPC_ccrhsub13 = trajMPC_ccrhs + 219;
trajMPC_FLOAT trajMPC_Phi13[11];
trajMPC_FLOAT trajMPC_W13[11];
trajMPC_FLOAT trajMPC_Ysd13[9];
trajMPC_FLOAT trajMPC_Lsd13[9];
trajMPC_FLOAT trajMPC_H14[3] = {2.0000000000000000E+001, 2.0000000000000000E+001, 0.0000000000000000E+000};
trajMPC_FLOAT* trajMPC_z14 = trajMPC_z + 154;
trajMPC_FLOAT* trajMPC_dzaff14 = trajMPC_dz_aff + 154;
trajMPC_FLOAT* trajMPC_dzcc14 = trajMPC_dz_cc + 154;
trajMPC_FLOAT* trajMPC_rd14 = trajMPC_rd + 154;
trajMPC_FLOAT trajMPC_Lbyrd14[3];
trajMPC_FLOAT* trajMPC_grad_cost14 = trajMPC_grad_cost + 154;
trajMPC_FLOAT* trajMPC_grad_eq14 = trajMPC_grad_eq + 154;
trajMPC_FLOAT* trajMPC_grad_ineq14 = trajMPC_grad_ineq + 154;
trajMPC_FLOAT trajMPC_ctv14[3];
trajMPC_FLOAT* trajMPC_v14 = trajMPC_v + 42;
trajMPC_FLOAT trajMPC_re14[3];
trajMPC_FLOAT trajMPC_beta14[3];
trajMPC_FLOAT trajMPC_betacc14[3];
trajMPC_FLOAT* trajMPC_dvaff14 = trajMPC_dv_aff + 42;
trajMPC_FLOAT* trajMPC_dvcc14 = trajMPC_dv_cc + 42;
trajMPC_FLOAT trajMPC_V14[9];
trajMPC_FLOAT trajMPC_Yd14[6];
trajMPC_FLOAT trajMPC_Ld14[6];
trajMPC_FLOAT trajMPC_yy14[3];
trajMPC_FLOAT trajMPC_bmy14[3];
int trajMPC_lbIdx14[3] = {0, 1, 2};
trajMPC_FLOAT* trajMPC_llb14 = trajMPC_l + 224;
trajMPC_FLOAT* trajMPC_slb14 = trajMPC_s + 224;
trajMPC_FLOAT* trajMPC_llbbyslb14 = trajMPC_lbys + 224;
trajMPC_FLOAT trajMPC_rilb14[3];
trajMPC_FLOAT* trajMPC_dllbaff14 = trajMPC_dl_aff + 224;
trajMPC_FLOAT* trajMPC_dslbaff14 = trajMPC_ds_aff + 224;
trajMPC_FLOAT* trajMPC_dllbcc14 = trajMPC_dl_cc + 224;
trajMPC_FLOAT* trajMPC_dslbcc14 = trajMPC_ds_cc + 224;
trajMPC_FLOAT* trajMPC_ccrhsl14 = trajMPC_ccrhs + 224;
int trajMPC_ubIdx14[3] = {0, 1, 2};
trajMPC_FLOAT* trajMPC_lub14 = trajMPC_l + 227;
trajMPC_FLOAT* trajMPC_sub14 = trajMPC_s + 227;
trajMPC_FLOAT* trajMPC_lubbysub14 = trajMPC_lbys + 227;
trajMPC_FLOAT trajMPC_riub14[3];
trajMPC_FLOAT* trajMPC_dlubaff14 = trajMPC_dl_aff + 227;
trajMPC_FLOAT* trajMPC_dsubaff14 = trajMPC_ds_aff + 227;
trajMPC_FLOAT* trajMPC_dlubcc14 = trajMPC_dl_cc + 227;
trajMPC_FLOAT* trajMPC_dsubcc14 = trajMPC_ds_cc + 227;
trajMPC_FLOAT* trajMPC_ccrhsub14 = trajMPC_ccrhs + 227;
trajMPC_FLOAT trajMPC_Phi14[3];
trajMPC_FLOAT trajMPC_D14[3] = {-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000};
trajMPC_FLOAT trajMPC_W14[3];
trajMPC_FLOAT trajMPC_Ysd14[9];
trajMPC_FLOAT trajMPC_Lsd14[9];
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
trajMPC_LA_INITIALIZEVECTOR_157(trajMPC_z, 0);
trajMPC_LA_INITIALIZEVECTOR_45(trajMPC_v, 1);
trajMPC_LA_INITIALIZEVECTOR_230(trajMPC_l, 1);
trajMPC_LA_INITIALIZEVECTOR_230(trajMPC_s, 1);
info->mu = 0;
trajMPC_LA_DOTACC_230(trajMPC_l, trajMPC_s, &info->mu);
info->mu /= 230;
while( 1 ){
info->pobj = 0;
trajMPC_LA_DIAG_QUADFCN_11(trajMPC_H00, params->f1, trajMPC_z00, trajMPC_grad_cost00, &info->pobj);
trajMPC_LA_DIAG_QUADFCN_11(trajMPC_H01, params->f2, trajMPC_z01, trajMPC_grad_cost01, &info->pobj);
trajMPC_LA_DIAG_QUADFCN_11(trajMPC_H01, params->f3, trajMPC_z02, trajMPC_grad_cost02, &info->pobj);
trajMPC_LA_DIAG_QUADFCN_11(trajMPC_H01, params->f4, trajMPC_z03, trajMPC_grad_cost03, &info->pobj);
trajMPC_LA_DIAG_QUADFCN_11(trajMPC_H01, params->f5, trajMPC_z04, trajMPC_grad_cost04, &info->pobj);
trajMPC_LA_DIAG_QUADFCN_11(trajMPC_H01, params->f6, trajMPC_z05, trajMPC_grad_cost05, &info->pobj);
trajMPC_LA_DIAG_QUADFCN_11(trajMPC_H01, params->f7, trajMPC_z06, trajMPC_grad_cost06, &info->pobj);
trajMPC_LA_DIAG_QUADFCN_11(trajMPC_H01, params->f8, trajMPC_z07, trajMPC_grad_cost07, &info->pobj);
trajMPC_LA_DIAG_QUADFCN_11(trajMPC_H01, params->f9, trajMPC_z08, trajMPC_grad_cost08, &info->pobj);
trajMPC_LA_DIAG_QUADFCN_11(trajMPC_H01, params->f10, trajMPC_z09, trajMPC_grad_cost09, &info->pobj);
trajMPC_LA_DIAG_QUADFCN_11(trajMPC_H01, params->f11, trajMPC_z10, trajMPC_grad_cost10, &info->pobj);
trajMPC_LA_DIAG_QUADFCN_11(trajMPC_H01, params->f12, trajMPC_z11, trajMPC_grad_cost11, &info->pobj);
trajMPC_LA_DIAG_QUADFCN_11(trajMPC_H01, params->f13, trajMPC_z12, trajMPC_grad_cost12, &info->pobj);
trajMPC_LA_DIAG_QUADFCN_11(trajMPC_H01, params->f14, trajMPC_z13, trajMPC_grad_cost13, &info->pobj);
trajMPC_LA_DIAG_QUADFCN_3(trajMPC_H14, params->f15, trajMPC_z14, trajMPC_grad_cost14, &info->pobj);
info->res_eq = 0;
info->dgap = 0;
trajMPC_LA_DIAGZERO_MVMSUB6_3(trajMPC_D00, trajMPC_z00, params->e1, trajMPC_v00, trajMPC_re00, &info->dgap, &info->res_eq);
trajMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C1, trajMPC_z00, trajMPC_D01, trajMPC_z01, params->e2, trajMPC_v01, trajMPC_re01, &info->dgap, &info->res_eq);
trajMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C2, trajMPC_z01, trajMPC_D01, trajMPC_z02, params->e3, trajMPC_v02, trajMPC_re02, &info->dgap, &info->res_eq);
trajMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C3, trajMPC_z02, trajMPC_D01, trajMPC_z03, params->e4, trajMPC_v03, trajMPC_re03, &info->dgap, &info->res_eq);
trajMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C4, trajMPC_z03, trajMPC_D01, trajMPC_z04, params->e5, trajMPC_v04, trajMPC_re04, &info->dgap, &info->res_eq);
trajMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C5, trajMPC_z04, trajMPC_D01, trajMPC_z05, params->e6, trajMPC_v05, trajMPC_re05, &info->dgap, &info->res_eq);
trajMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C6, trajMPC_z05, trajMPC_D01, trajMPC_z06, params->e7, trajMPC_v06, trajMPC_re06, &info->dgap, &info->res_eq);
trajMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C7, trajMPC_z06, trajMPC_D01, trajMPC_z07, params->e8, trajMPC_v07, trajMPC_re07, &info->dgap, &info->res_eq);
trajMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C8, trajMPC_z07, trajMPC_D01, trajMPC_z08, params->e9, trajMPC_v08, trajMPC_re08, &info->dgap, &info->res_eq);
trajMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C9, trajMPC_z08, trajMPC_D01, trajMPC_z09, params->e10, trajMPC_v09, trajMPC_re09, &info->dgap, &info->res_eq);
trajMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C10, trajMPC_z09, trajMPC_D01, trajMPC_z10, params->e11, trajMPC_v10, trajMPC_re10, &info->dgap, &info->res_eq);
trajMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C11, trajMPC_z10, trajMPC_D01, trajMPC_z11, params->e12, trajMPC_v11, trajMPC_re11, &info->dgap, &info->res_eq);
trajMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C12, trajMPC_z11, trajMPC_D01, trajMPC_z12, params->e13, trajMPC_v12, trajMPC_re12, &info->dgap, &info->res_eq);
trajMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_11(params->C13, trajMPC_z12, trajMPC_D01, trajMPC_z13, params->e14, trajMPC_v13, trajMPC_re13, &info->dgap, &info->res_eq);
trajMPC_LA_DENSE_DIAGZERO_MVMSUB3_3_11_3(params->C14, trajMPC_z13, trajMPC_D14, trajMPC_z14, params->e15, trajMPC_v14, trajMPC_re14, &info->dgap, &info->res_eq);
trajMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C1, trajMPC_v01, trajMPC_D00, trajMPC_v00, trajMPC_grad_eq00);
trajMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C2, trajMPC_v02, trajMPC_D01, trajMPC_v01, trajMPC_grad_eq01);
trajMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C3, trajMPC_v03, trajMPC_D01, trajMPC_v02, trajMPC_grad_eq02);
trajMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C4, trajMPC_v04, trajMPC_D01, trajMPC_v03, trajMPC_grad_eq03);
trajMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C5, trajMPC_v05, trajMPC_D01, trajMPC_v04, trajMPC_grad_eq04);
trajMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C6, trajMPC_v06, trajMPC_D01, trajMPC_v05, trajMPC_grad_eq05);
trajMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C7, trajMPC_v07, trajMPC_D01, trajMPC_v06, trajMPC_grad_eq06);
trajMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C8, trajMPC_v08, trajMPC_D01, trajMPC_v07, trajMPC_grad_eq07);
trajMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C9, trajMPC_v09, trajMPC_D01, trajMPC_v08, trajMPC_grad_eq08);
trajMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C10, trajMPC_v10, trajMPC_D01, trajMPC_v09, trajMPC_grad_eq09);
trajMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C11, trajMPC_v11, trajMPC_D01, trajMPC_v10, trajMPC_grad_eq10);
trajMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C12, trajMPC_v12, trajMPC_D01, trajMPC_v11, trajMPC_grad_eq11);
trajMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C13, trajMPC_v13, trajMPC_D01, trajMPC_v12, trajMPC_grad_eq12);
trajMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C14, trajMPC_v14, trajMPC_D01, trajMPC_v13, trajMPC_grad_eq13);
trajMPC_LA_DIAGZERO_MTVM_3_3(trajMPC_D14, trajMPC_v14, trajMPC_grad_eq14);
info->res_ineq = 0;
trajMPC_LA_VSUBADD3_11(params->lb1, trajMPC_z00, trajMPC_lbIdx00, trajMPC_llb00, trajMPC_slb00, trajMPC_rilb00, &info->dgap, &info->res_ineq);
trajMPC_LA_VSUBADD2_5(trajMPC_z00, trajMPC_ubIdx00, params->ub1, trajMPC_lub00, trajMPC_sub00, trajMPC_riub00, &info->dgap, &info->res_ineq);
trajMPC_LA_VSUBADD3_11(params->lb2, trajMPC_z01, trajMPC_lbIdx01, trajMPC_llb01, trajMPC_slb01, trajMPC_rilb01, &info->dgap, &info->res_ineq);
trajMPC_LA_VSUBADD2_5(trajMPC_z01, trajMPC_ubIdx01, params->ub2, trajMPC_lub01, trajMPC_sub01, trajMPC_riub01, &info->dgap, &info->res_ineq);
trajMPC_LA_VSUBADD3_11(params->lb3, trajMPC_z02, trajMPC_lbIdx02, trajMPC_llb02, trajMPC_slb02, trajMPC_rilb02, &info->dgap, &info->res_ineq);
trajMPC_LA_VSUBADD2_5(trajMPC_z02, trajMPC_ubIdx02, params->ub3, trajMPC_lub02, trajMPC_sub02, trajMPC_riub02, &info->dgap, &info->res_ineq);
trajMPC_LA_VSUBADD3_11(params->lb4, trajMPC_z03, trajMPC_lbIdx03, trajMPC_llb03, trajMPC_slb03, trajMPC_rilb03, &info->dgap, &info->res_ineq);
trajMPC_LA_VSUBADD2_5(trajMPC_z03, trajMPC_ubIdx03, params->ub4, trajMPC_lub03, trajMPC_sub03, trajMPC_riub03, &info->dgap, &info->res_ineq);
trajMPC_LA_VSUBADD3_11(params->lb5, trajMPC_z04, trajMPC_lbIdx04, trajMPC_llb04, trajMPC_slb04, trajMPC_rilb04, &info->dgap, &info->res_ineq);
trajMPC_LA_VSUBADD2_5(trajMPC_z04, trajMPC_ubIdx04, params->ub5, trajMPC_lub04, trajMPC_sub04, trajMPC_riub04, &info->dgap, &info->res_ineq);
trajMPC_LA_VSUBADD3_11(params->lb6, trajMPC_z05, trajMPC_lbIdx05, trajMPC_llb05, trajMPC_slb05, trajMPC_rilb05, &info->dgap, &info->res_ineq);
trajMPC_LA_VSUBADD2_5(trajMPC_z05, trajMPC_ubIdx05, params->ub6, trajMPC_lub05, trajMPC_sub05, trajMPC_riub05, &info->dgap, &info->res_ineq);
trajMPC_LA_VSUBADD3_11(params->lb7, trajMPC_z06, trajMPC_lbIdx06, trajMPC_llb06, trajMPC_slb06, trajMPC_rilb06, &info->dgap, &info->res_ineq);
trajMPC_LA_VSUBADD2_5(trajMPC_z06, trajMPC_ubIdx06, params->ub7, trajMPC_lub06, trajMPC_sub06, trajMPC_riub06, &info->dgap, &info->res_ineq);
trajMPC_LA_VSUBADD3_11(params->lb8, trajMPC_z07, trajMPC_lbIdx07, trajMPC_llb07, trajMPC_slb07, trajMPC_rilb07, &info->dgap, &info->res_ineq);
trajMPC_LA_VSUBADD2_5(trajMPC_z07, trajMPC_ubIdx07, params->ub8, trajMPC_lub07, trajMPC_sub07, trajMPC_riub07, &info->dgap, &info->res_ineq);
trajMPC_LA_VSUBADD3_11(params->lb9, trajMPC_z08, trajMPC_lbIdx08, trajMPC_llb08, trajMPC_slb08, trajMPC_rilb08, &info->dgap, &info->res_ineq);
trajMPC_LA_VSUBADD2_5(trajMPC_z08, trajMPC_ubIdx08, params->ub9, trajMPC_lub08, trajMPC_sub08, trajMPC_riub08, &info->dgap, &info->res_ineq);
trajMPC_LA_VSUBADD3_11(params->lb10, trajMPC_z09, trajMPC_lbIdx09, trajMPC_llb09, trajMPC_slb09, trajMPC_rilb09, &info->dgap, &info->res_ineq);
trajMPC_LA_VSUBADD2_5(trajMPC_z09, trajMPC_ubIdx09, params->ub10, trajMPC_lub09, trajMPC_sub09, trajMPC_riub09, &info->dgap, &info->res_ineq);
trajMPC_LA_VSUBADD3_11(params->lb11, trajMPC_z10, trajMPC_lbIdx10, trajMPC_llb10, trajMPC_slb10, trajMPC_rilb10, &info->dgap, &info->res_ineq);
trajMPC_LA_VSUBADD2_5(trajMPC_z10, trajMPC_ubIdx10, params->ub11, trajMPC_lub10, trajMPC_sub10, trajMPC_riub10, &info->dgap, &info->res_ineq);
trajMPC_LA_VSUBADD3_11(params->lb12, trajMPC_z11, trajMPC_lbIdx11, trajMPC_llb11, trajMPC_slb11, trajMPC_rilb11, &info->dgap, &info->res_ineq);
trajMPC_LA_VSUBADD2_5(trajMPC_z11, trajMPC_ubIdx11, params->ub12, trajMPC_lub11, trajMPC_sub11, trajMPC_riub11, &info->dgap, &info->res_ineq);
trajMPC_LA_VSUBADD3_11(params->lb13, trajMPC_z12, trajMPC_lbIdx12, trajMPC_llb12, trajMPC_slb12, trajMPC_rilb12, &info->dgap, &info->res_ineq);
trajMPC_LA_VSUBADD2_5(trajMPC_z12, trajMPC_ubIdx12, params->ub13, trajMPC_lub12, trajMPC_sub12, trajMPC_riub12, &info->dgap, &info->res_ineq);
trajMPC_LA_VSUBADD3_11(params->lb14, trajMPC_z13, trajMPC_lbIdx13, trajMPC_llb13, trajMPC_slb13, trajMPC_rilb13, &info->dgap, &info->res_ineq);
trajMPC_LA_VSUBADD2_5(trajMPC_z13, trajMPC_ubIdx13, params->ub14, trajMPC_lub13, trajMPC_sub13, trajMPC_riub13, &info->dgap, &info->res_ineq);
trajMPC_LA_VSUBADD3_3(params->lb15, trajMPC_z14, trajMPC_lbIdx14, trajMPC_llb14, trajMPC_slb14, trajMPC_rilb14, &info->dgap, &info->res_ineq);
trajMPC_LA_VSUBADD2_3(trajMPC_z14, trajMPC_ubIdx14, params->ub15, trajMPC_lub14, trajMPC_sub14, trajMPC_riub14, &info->dgap, &info->res_ineq);
trajMPC_LA_INEQ_B_GRAD_11_11_5(trajMPC_lub00, trajMPC_sub00, trajMPC_riub00, trajMPC_llb00, trajMPC_slb00, trajMPC_rilb00, trajMPC_lbIdx00, trajMPC_ubIdx00, trajMPC_grad_ineq00, trajMPC_lubbysub00, trajMPC_llbbyslb00);
trajMPC_LA_INEQ_B_GRAD_11_11_5(trajMPC_lub01, trajMPC_sub01, trajMPC_riub01, trajMPC_llb01, trajMPC_slb01, trajMPC_rilb01, trajMPC_lbIdx01, trajMPC_ubIdx01, trajMPC_grad_ineq01, trajMPC_lubbysub01, trajMPC_llbbyslb01);
trajMPC_LA_INEQ_B_GRAD_11_11_5(trajMPC_lub02, trajMPC_sub02, trajMPC_riub02, trajMPC_llb02, trajMPC_slb02, trajMPC_rilb02, trajMPC_lbIdx02, trajMPC_ubIdx02, trajMPC_grad_ineq02, trajMPC_lubbysub02, trajMPC_llbbyslb02);
trajMPC_LA_INEQ_B_GRAD_11_11_5(trajMPC_lub03, trajMPC_sub03, trajMPC_riub03, trajMPC_llb03, trajMPC_slb03, trajMPC_rilb03, trajMPC_lbIdx03, trajMPC_ubIdx03, trajMPC_grad_ineq03, trajMPC_lubbysub03, trajMPC_llbbyslb03);
trajMPC_LA_INEQ_B_GRAD_11_11_5(trajMPC_lub04, trajMPC_sub04, trajMPC_riub04, trajMPC_llb04, trajMPC_slb04, trajMPC_rilb04, trajMPC_lbIdx04, trajMPC_ubIdx04, trajMPC_grad_ineq04, trajMPC_lubbysub04, trajMPC_llbbyslb04);
trajMPC_LA_INEQ_B_GRAD_11_11_5(trajMPC_lub05, trajMPC_sub05, trajMPC_riub05, trajMPC_llb05, trajMPC_slb05, trajMPC_rilb05, trajMPC_lbIdx05, trajMPC_ubIdx05, trajMPC_grad_ineq05, trajMPC_lubbysub05, trajMPC_llbbyslb05);
trajMPC_LA_INEQ_B_GRAD_11_11_5(trajMPC_lub06, trajMPC_sub06, trajMPC_riub06, trajMPC_llb06, trajMPC_slb06, trajMPC_rilb06, trajMPC_lbIdx06, trajMPC_ubIdx06, trajMPC_grad_ineq06, trajMPC_lubbysub06, trajMPC_llbbyslb06);
trajMPC_LA_INEQ_B_GRAD_11_11_5(trajMPC_lub07, trajMPC_sub07, trajMPC_riub07, trajMPC_llb07, trajMPC_slb07, trajMPC_rilb07, trajMPC_lbIdx07, trajMPC_ubIdx07, trajMPC_grad_ineq07, trajMPC_lubbysub07, trajMPC_llbbyslb07);
trajMPC_LA_INEQ_B_GRAD_11_11_5(trajMPC_lub08, trajMPC_sub08, trajMPC_riub08, trajMPC_llb08, trajMPC_slb08, trajMPC_rilb08, trajMPC_lbIdx08, trajMPC_ubIdx08, trajMPC_grad_ineq08, trajMPC_lubbysub08, trajMPC_llbbyslb08);
trajMPC_LA_INEQ_B_GRAD_11_11_5(trajMPC_lub09, trajMPC_sub09, trajMPC_riub09, trajMPC_llb09, trajMPC_slb09, trajMPC_rilb09, trajMPC_lbIdx09, trajMPC_ubIdx09, trajMPC_grad_ineq09, trajMPC_lubbysub09, trajMPC_llbbyslb09);
trajMPC_LA_INEQ_B_GRAD_11_11_5(trajMPC_lub10, trajMPC_sub10, trajMPC_riub10, trajMPC_llb10, trajMPC_slb10, trajMPC_rilb10, trajMPC_lbIdx10, trajMPC_ubIdx10, trajMPC_grad_ineq10, trajMPC_lubbysub10, trajMPC_llbbyslb10);
trajMPC_LA_INEQ_B_GRAD_11_11_5(trajMPC_lub11, trajMPC_sub11, trajMPC_riub11, trajMPC_llb11, trajMPC_slb11, trajMPC_rilb11, trajMPC_lbIdx11, trajMPC_ubIdx11, trajMPC_grad_ineq11, trajMPC_lubbysub11, trajMPC_llbbyslb11);
trajMPC_LA_INEQ_B_GRAD_11_11_5(trajMPC_lub12, trajMPC_sub12, trajMPC_riub12, trajMPC_llb12, trajMPC_slb12, trajMPC_rilb12, trajMPC_lbIdx12, trajMPC_ubIdx12, trajMPC_grad_ineq12, trajMPC_lubbysub12, trajMPC_llbbyslb12);
trajMPC_LA_INEQ_B_GRAD_11_11_5(trajMPC_lub13, trajMPC_sub13, trajMPC_riub13, trajMPC_llb13, trajMPC_slb13, trajMPC_rilb13, trajMPC_lbIdx13, trajMPC_ubIdx13, trajMPC_grad_ineq13, trajMPC_lubbysub13, trajMPC_llbbyslb13);
trajMPC_LA_INEQ_B_GRAD_3_3_3(trajMPC_lub14, trajMPC_sub14, trajMPC_riub14, trajMPC_llb14, trajMPC_slb14, trajMPC_rilb14, trajMPC_lbIdx14, trajMPC_ubIdx14, trajMPC_grad_ineq14, trajMPC_lubbysub14, trajMPC_llbbyslb14);
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
trajMPC_LA_VVADD3_157(trajMPC_grad_cost, trajMPC_grad_eq, trajMPC_grad_ineq, trajMPC_rd);
trajMPC_LA_DIAG_CHOL_LBUB_11_11_5(trajMPC_H00, trajMPC_llbbyslb00, trajMPC_lbIdx00, trajMPC_lubbysub00, trajMPC_ubIdx00, trajMPC_Phi00);
trajMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(trajMPC_Phi00, params->C1, trajMPC_V00);
trajMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(trajMPC_Phi00, trajMPC_D00, trajMPC_W00);
trajMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(trajMPC_W00, trajMPC_V00, trajMPC_Ysd01);
trajMPC_LA_DIAG_FORWARDSUB_11(trajMPC_Phi00, trajMPC_rd00, trajMPC_Lbyrd00);
trajMPC_LA_DIAG_CHOL_LBUB_5_11_5(trajMPC_H01, trajMPC_llbbyslb01, trajMPC_lbIdx01, trajMPC_lubbysub01, trajMPC_ubIdx01, trajMPC_Phi01);
trajMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(trajMPC_Phi01, params->C2, trajMPC_V01);
trajMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(trajMPC_Phi01, trajMPC_D01, trajMPC_W01);
trajMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(trajMPC_W01, trajMPC_V01, trajMPC_Ysd02);
trajMPC_LA_DIAG_FORWARDSUB_11(trajMPC_Phi01, trajMPC_rd01, trajMPC_Lbyrd01);
trajMPC_LA_DIAG_CHOL_LBUB_11_11_5(trajMPC_H01, trajMPC_llbbyslb02, trajMPC_lbIdx02, trajMPC_lubbysub02, trajMPC_ubIdx02, trajMPC_Phi02);
trajMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(trajMPC_Phi02, params->C3, trajMPC_V02);
trajMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(trajMPC_Phi02, trajMPC_D01, trajMPC_W02);
trajMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(trajMPC_W02, trajMPC_V02, trajMPC_Ysd03);
trajMPC_LA_DIAG_FORWARDSUB_11(trajMPC_Phi02, trajMPC_rd02, trajMPC_Lbyrd02);
trajMPC_LA_DIAG_CHOL_LBUB_11_11_5(trajMPC_H01, trajMPC_llbbyslb03, trajMPC_lbIdx03, trajMPC_lubbysub03, trajMPC_ubIdx03, trajMPC_Phi03);
trajMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(trajMPC_Phi03, params->C4, trajMPC_V03);
trajMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(trajMPC_Phi03, trajMPC_D01, trajMPC_W03);
trajMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(trajMPC_W03, trajMPC_V03, trajMPC_Ysd04);
trajMPC_LA_DIAG_FORWARDSUB_11(trajMPC_Phi03, trajMPC_rd03, trajMPC_Lbyrd03);
trajMPC_LA_DIAG_CHOL_LBUB_11_11_5(trajMPC_H01, trajMPC_llbbyslb04, trajMPC_lbIdx04, trajMPC_lubbysub04, trajMPC_ubIdx04, trajMPC_Phi04);
trajMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(trajMPC_Phi04, params->C5, trajMPC_V04);
trajMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(trajMPC_Phi04, trajMPC_D01, trajMPC_W04);
trajMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(trajMPC_W04, trajMPC_V04, trajMPC_Ysd05);
trajMPC_LA_DIAG_FORWARDSUB_11(trajMPC_Phi04, trajMPC_rd04, trajMPC_Lbyrd04);
trajMPC_LA_DIAG_CHOL_LBUB_11_11_5(trajMPC_H01, trajMPC_llbbyslb05, trajMPC_lbIdx05, trajMPC_lubbysub05, trajMPC_ubIdx05, trajMPC_Phi05);
trajMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(trajMPC_Phi05, params->C6, trajMPC_V05);
trajMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(trajMPC_Phi05, trajMPC_D01, trajMPC_W05);
trajMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(trajMPC_W05, trajMPC_V05, trajMPC_Ysd06);
trajMPC_LA_DIAG_FORWARDSUB_11(trajMPC_Phi05, trajMPC_rd05, trajMPC_Lbyrd05);
trajMPC_LA_DIAG_CHOL_LBUB_11_11_5(trajMPC_H01, trajMPC_llbbyslb06, trajMPC_lbIdx06, trajMPC_lubbysub06, trajMPC_ubIdx06, trajMPC_Phi06);
trajMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(trajMPC_Phi06, params->C7, trajMPC_V06);
trajMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(trajMPC_Phi06, trajMPC_D01, trajMPC_W06);
trajMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(trajMPC_W06, trajMPC_V06, trajMPC_Ysd07);
trajMPC_LA_DIAG_FORWARDSUB_11(trajMPC_Phi06, trajMPC_rd06, trajMPC_Lbyrd06);
trajMPC_LA_DIAG_CHOL_LBUB_11_11_5(trajMPC_H01, trajMPC_llbbyslb07, trajMPC_lbIdx07, trajMPC_lubbysub07, trajMPC_ubIdx07, trajMPC_Phi07);
trajMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(trajMPC_Phi07, params->C8, trajMPC_V07);
trajMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(trajMPC_Phi07, trajMPC_D01, trajMPC_W07);
trajMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(trajMPC_W07, trajMPC_V07, trajMPC_Ysd08);
trajMPC_LA_DIAG_FORWARDSUB_11(trajMPC_Phi07, trajMPC_rd07, trajMPC_Lbyrd07);
trajMPC_LA_DIAG_CHOL_LBUB_11_11_5(trajMPC_H01, trajMPC_llbbyslb08, trajMPC_lbIdx08, trajMPC_lubbysub08, trajMPC_ubIdx08, trajMPC_Phi08);
trajMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(trajMPC_Phi08, params->C9, trajMPC_V08);
trajMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(trajMPC_Phi08, trajMPC_D01, trajMPC_W08);
trajMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(trajMPC_W08, trajMPC_V08, trajMPC_Ysd09);
trajMPC_LA_DIAG_FORWARDSUB_11(trajMPC_Phi08, trajMPC_rd08, trajMPC_Lbyrd08);
trajMPC_LA_DIAG_CHOL_LBUB_11_11_5(trajMPC_H01, trajMPC_llbbyslb09, trajMPC_lbIdx09, trajMPC_lubbysub09, trajMPC_ubIdx09, trajMPC_Phi09);
trajMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(trajMPC_Phi09, params->C10, trajMPC_V09);
trajMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(trajMPC_Phi09, trajMPC_D01, trajMPC_W09);
trajMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(trajMPC_W09, trajMPC_V09, trajMPC_Ysd10);
trajMPC_LA_DIAG_FORWARDSUB_11(trajMPC_Phi09, trajMPC_rd09, trajMPC_Lbyrd09);
trajMPC_LA_DIAG_CHOL_LBUB_11_11_5(trajMPC_H01, trajMPC_llbbyslb10, trajMPC_lbIdx10, trajMPC_lubbysub10, trajMPC_ubIdx10, trajMPC_Phi10);
trajMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(trajMPC_Phi10, params->C11, trajMPC_V10);
trajMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(trajMPC_Phi10, trajMPC_D01, trajMPC_W10);
trajMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(trajMPC_W10, trajMPC_V10, trajMPC_Ysd11);
trajMPC_LA_DIAG_FORWARDSUB_11(trajMPC_Phi10, trajMPC_rd10, trajMPC_Lbyrd10);
trajMPC_LA_DIAG_CHOL_LBUB_11_11_5(trajMPC_H01, trajMPC_llbbyslb11, trajMPC_lbIdx11, trajMPC_lubbysub11, trajMPC_ubIdx11, trajMPC_Phi11);
trajMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(trajMPC_Phi11, params->C12, trajMPC_V11);
trajMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(trajMPC_Phi11, trajMPC_D01, trajMPC_W11);
trajMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(trajMPC_W11, trajMPC_V11, trajMPC_Ysd12);
trajMPC_LA_DIAG_FORWARDSUB_11(trajMPC_Phi11, trajMPC_rd11, trajMPC_Lbyrd11);
trajMPC_LA_DIAG_CHOL_LBUB_11_11_5(trajMPC_H01, trajMPC_llbbyslb12, trajMPC_lbIdx12, trajMPC_lubbysub12, trajMPC_ubIdx12, trajMPC_Phi12);
trajMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(trajMPC_Phi12, params->C13, trajMPC_V12);
trajMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(trajMPC_Phi12, trajMPC_D01, trajMPC_W12);
trajMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(trajMPC_W12, trajMPC_V12, trajMPC_Ysd13);
trajMPC_LA_DIAG_FORWARDSUB_11(trajMPC_Phi12, trajMPC_rd12, trajMPC_Lbyrd12);
trajMPC_LA_DIAG_CHOL_LBUB_11_11_5(trajMPC_H01, trajMPC_llbbyslb13, trajMPC_lbIdx13, trajMPC_lubbysub13, trajMPC_ubIdx13, trajMPC_Phi13);
trajMPC_LA_DIAG_MATRIXFORWARDSUB_3_11(trajMPC_Phi13, params->C14, trajMPC_V13);
trajMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_11(trajMPC_Phi13, trajMPC_D01, trajMPC_W13);
trajMPC_LA_DENSE_DIAGZERO_MMTM_3_11_3(trajMPC_W13, trajMPC_V13, trajMPC_Ysd14);
trajMPC_LA_DIAG_FORWARDSUB_11(trajMPC_Phi13, trajMPC_rd13, trajMPC_Lbyrd13);
trajMPC_LA_DIAG_CHOL_ONELOOP_LBUB_3_3_3(trajMPC_H14, trajMPC_llbbyslb14, trajMPC_lbIdx14, trajMPC_lubbysub14, trajMPC_ubIdx14, trajMPC_Phi14);
trajMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_3_3(trajMPC_Phi14, trajMPC_D14, trajMPC_W14);
trajMPC_LA_DIAG_FORWARDSUB_3(trajMPC_Phi14, trajMPC_rd14, trajMPC_Lbyrd14);
trajMPC_LA_DIAGZERO_MMT_3(trajMPC_W00, trajMPC_Yd00);
trajMPC_LA_DIAGZERO_MVMSUB7_3(trajMPC_W00, trajMPC_Lbyrd00, trajMPC_re00, trajMPC_beta00);
trajMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(trajMPC_V00, trajMPC_W01, trajMPC_Yd01);
trajMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(trajMPC_V00, trajMPC_Lbyrd00, trajMPC_W01, trajMPC_Lbyrd01, trajMPC_re01, trajMPC_beta01);
trajMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(trajMPC_V01, trajMPC_W02, trajMPC_Yd02);
trajMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(trajMPC_V01, trajMPC_Lbyrd01, trajMPC_W02, trajMPC_Lbyrd02, trajMPC_re02, trajMPC_beta02);
trajMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(trajMPC_V02, trajMPC_W03, trajMPC_Yd03);
trajMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(trajMPC_V02, trajMPC_Lbyrd02, trajMPC_W03, trajMPC_Lbyrd03, trajMPC_re03, trajMPC_beta03);
trajMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(trajMPC_V03, trajMPC_W04, trajMPC_Yd04);
trajMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(trajMPC_V03, trajMPC_Lbyrd03, trajMPC_W04, trajMPC_Lbyrd04, trajMPC_re04, trajMPC_beta04);
trajMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(trajMPC_V04, trajMPC_W05, trajMPC_Yd05);
trajMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(trajMPC_V04, trajMPC_Lbyrd04, trajMPC_W05, trajMPC_Lbyrd05, trajMPC_re05, trajMPC_beta05);
trajMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(trajMPC_V05, trajMPC_W06, trajMPC_Yd06);
trajMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(trajMPC_V05, trajMPC_Lbyrd05, trajMPC_W06, trajMPC_Lbyrd06, trajMPC_re06, trajMPC_beta06);
trajMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(trajMPC_V06, trajMPC_W07, trajMPC_Yd07);
trajMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(trajMPC_V06, trajMPC_Lbyrd06, trajMPC_W07, trajMPC_Lbyrd07, trajMPC_re07, trajMPC_beta07);
trajMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(trajMPC_V07, trajMPC_W08, trajMPC_Yd08);
trajMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(trajMPC_V07, trajMPC_Lbyrd07, trajMPC_W08, trajMPC_Lbyrd08, trajMPC_re08, trajMPC_beta08);
trajMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(trajMPC_V08, trajMPC_W09, trajMPC_Yd09);
trajMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(trajMPC_V08, trajMPC_Lbyrd08, trajMPC_W09, trajMPC_Lbyrd09, trajMPC_re09, trajMPC_beta09);
trajMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(trajMPC_V09, trajMPC_W10, trajMPC_Yd10);
trajMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(trajMPC_V09, trajMPC_Lbyrd09, trajMPC_W10, trajMPC_Lbyrd10, trajMPC_re10, trajMPC_beta10);
trajMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(trajMPC_V10, trajMPC_W11, trajMPC_Yd11);
trajMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(trajMPC_V10, trajMPC_Lbyrd10, trajMPC_W11, trajMPC_Lbyrd11, trajMPC_re11, trajMPC_beta11);
trajMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(trajMPC_V11, trajMPC_W12, trajMPC_Yd12);
trajMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(trajMPC_V11, trajMPC_Lbyrd11, trajMPC_W12, trajMPC_Lbyrd12, trajMPC_re12, trajMPC_beta12);
trajMPC_LA_DENSE_DIAGZERO_MMT2_3_11_11(trajMPC_V12, trajMPC_W13, trajMPC_Yd13);
trajMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_11(trajMPC_V12, trajMPC_Lbyrd12, trajMPC_W13, trajMPC_Lbyrd13, trajMPC_re13, trajMPC_beta13);
trajMPC_LA_DENSE_DIAGZERO_MMT2_3_11_3(trajMPC_V13, trajMPC_W14, trajMPC_Yd14);
trajMPC_LA_DENSE_DIAGZERO_2MVMSUB2_3_11_3(trajMPC_V13, trajMPC_Lbyrd13, trajMPC_W14, trajMPC_Lbyrd14, trajMPC_re14, trajMPC_beta14);
trajMPC_LA_DENSE_CHOL_3(trajMPC_Yd00, trajMPC_Ld00);
trajMPC_LA_DENSE_FORWARDSUB_3(trajMPC_Ld00, trajMPC_beta00, trajMPC_yy00);
trajMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(trajMPC_Ld00, trajMPC_Ysd01, trajMPC_Lsd01);
trajMPC_LA_DENSE_MMTSUB_3_3(trajMPC_Lsd01, trajMPC_Yd01);
trajMPC_LA_DENSE_CHOL_3(trajMPC_Yd01, trajMPC_Ld01);
trajMPC_LA_DENSE_MVMSUB1_3_3(trajMPC_Lsd01, trajMPC_yy00, trajMPC_beta01, trajMPC_bmy01);
trajMPC_LA_DENSE_FORWARDSUB_3(trajMPC_Ld01, trajMPC_bmy01, trajMPC_yy01);
trajMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(trajMPC_Ld01, trajMPC_Ysd02, trajMPC_Lsd02);
trajMPC_LA_DENSE_MMTSUB_3_3(trajMPC_Lsd02, trajMPC_Yd02);
trajMPC_LA_DENSE_CHOL_3(trajMPC_Yd02, trajMPC_Ld02);
trajMPC_LA_DENSE_MVMSUB1_3_3(trajMPC_Lsd02, trajMPC_yy01, trajMPC_beta02, trajMPC_bmy02);
trajMPC_LA_DENSE_FORWARDSUB_3(trajMPC_Ld02, trajMPC_bmy02, trajMPC_yy02);
trajMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(trajMPC_Ld02, trajMPC_Ysd03, trajMPC_Lsd03);
trajMPC_LA_DENSE_MMTSUB_3_3(trajMPC_Lsd03, trajMPC_Yd03);
trajMPC_LA_DENSE_CHOL_3(trajMPC_Yd03, trajMPC_Ld03);
trajMPC_LA_DENSE_MVMSUB1_3_3(trajMPC_Lsd03, trajMPC_yy02, trajMPC_beta03, trajMPC_bmy03);
trajMPC_LA_DENSE_FORWARDSUB_3(trajMPC_Ld03, trajMPC_bmy03, trajMPC_yy03);
trajMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(trajMPC_Ld03, trajMPC_Ysd04, trajMPC_Lsd04);
trajMPC_LA_DENSE_MMTSUB_3_3(trajMPC_Lsd04, trajMPC_Yd04);
trajMPC_LA_DENSE_CHOL_3(trajMPC_Yd04, trajMPC_Ld04);
trajMPC_LA_DENSE_MVMSUB1_3_3(trajMPC_Lsd04, trajMPC_yy03, trajMPC_beta04, trajMPC_bmy04);
trajMPC_LA_DENSE_FORWARDSUB_3(trajMPC_Ld04, trajMPC_bmy04, trajMPC_yy04);
trajMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(trajMPC_Ld04, trajMPC_Ysd05, trajMPC_Lsd05);
trajMPC_LA_DENSE_MMTSUB_3_3(trajMPC_Lsd05, trajMPC_Yd05);
trajMPC_LA_DENSE_CHOL_3(trajMPC_Yd05, trajMPC_Ld05);
trajMPC_LA_DENSE_MVMSUB1_3_3(trajMPC_Lsd05, trajMPC_yy04, trajMPC_beta05, trajMPC_bmy05);
trajMPC_LA_DENSE_FORWARDSUB_3(trajMPC_Ld05, trajMPC_bmy05, trajMPC_yy05);
trajMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(trajMPC_Ld05, trajMPC_Ysd06, trajMPC_Lsd06);
trajMPC_LA_DENSE_MMTSUB_3_3(trajMPC_Lsd06, trajMPC_Yd06);
trajMPC_LA_DENSE_CHOL_3(trajMPC_Yd06, trajMPC_Ld06);
trajMPC_LA_DENSE_MVMSUB1_3_3(trajMPC_Lsd06, trajMPC_yy05, trajMPC_beta06, trajMPC_bmy06);
trajMPC_LA_DENSE_FORWARDSUB_3(trajMPC_Ld06, trajMPC_bmy06, trajMPC_yy06);
trajMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(trajMPC_Ld06, trajMPC_Ysd07, trajMPC_Lsd07);
trajMPC_LA_DENSE_MMTSUB_3_3(trajMPC_Lsd07, trajMPC_Yd07);
trajMPC_LA_DENSE_CHOL_3(trajMPC_Yd07, trajMPC_Ld07);
trajMPC_LA_DENSE_MVMSUB1_3_3(trajMPC_Lsd07, trajMPC_yy06, trajMPC_beta07, trajMPC_bmy07);
trajMPC_LA_DENSE_FORWARDSUB_3(trajMPC_Ld07, trajMPC_bmy07, trajMPC_yy07);
trajMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(trajMPC_Ld07, trajMPC_Ysd08, trajMPC_Lsd08);
trajMPC_LA_DENSE_MMTSUB_3_3(trajMPC_Lsd08, trajMPC_Yd08);
trajMPC_LA_DENSE_CHOL_3(trajMPC_Yd08, trajMPC_Ld08);
trajMPC_LA_DENSE_MVMSUB1_3_3(trajMPC_Lsd08, trajMPC_yy07, trajMPC_beta08, trajMPC_bmy08);
trajMPC_LA_DENSE_FORWARDSUB_3(trajMPC_Ld08, trajMPC_bmy08, trajMPC_yy08);
trajMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(trajMPC_Ld08, trajMPC_Ysd09, trajMPC_Lsd09);
trajMPC_LA_DENSE_MMTSUB_3_3(trajMPC_Lsd09, trajMPC_Yd09);
trajMPC_LA_DENSE_CHOL_3(trajMPC_Yd09, trajMPC_Ld09);
trajMPC_LA_DENSE_MVMSUB1_3_3(trajMPC_Lsd09, trajMPC_yy08, trajMPC_beta09, trajMPC_bmy09);
trajMPC_LA_DENSE_FORWARDSUB_3(trajMPC_Ld09, trajMPC_bmy09, trajMPC_yy09);
trajMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(trajMPC_Ld09, trajMPC_Ysd10, trajMPC_Lsd10);
trajMPC_LA_DENSE_MMTSUB_3_3(trajMPC_Lsd10, trajMPC_Yd10);
trajMPC_LA_DENSE_CHOL_3(trajMPC_Yd10, trajMPC_Ld10);
trajMPC_LA_DENSE_MVMSUB1_3_3(trajMPC_Lsd10, trajMPC_yy09, trajMPC_beta10, trajMPC_bmy10);
trajMPC_LA_DENSE_FORWARDSUB_3(trajMPC_Ld10, trajMPC_bmy10, trajMPC_yy10);
trajMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(trajMPC_Ld10, trajMPC_Ysd11, trajMPC_Lsd11);
trajMPC_LA_DENSE_MMTSUB_3_3(trajMPC_Lsd11, trajMPC_Yd11);
trajMPC_LA_DENSE_CHOL_3(trajMPC_Yd11, trajMPC_Ld11);
trajMPC_LA_DENSE_MVMSUB1_3_3(trajMPC_Lsd11, trajMPC_yy10, trajMPC_beta11, trajMPC_bmy11);
trajMPC_LA_DENSE_FORWARDSUB_3(trajMPC_Ld11, trajMPC_bmy11, trajMPC_yy11);
trajMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(trajMPC_Ld11, trajMPC_Ysd12, trajMPC_Lsd12);
trajMPC_LA_DENSE_MMTSUB_3_3(trajMPC_Lsd12, trajMPC_Yd12);
trajMPC_LA_DENSE_CHOL_3(trajMPC_Yd12, trajMPC_Ld12);
trajMPC_LA_DENSE_MVMSUB1_3_3(trajMPC_Lsd12, trajMPC_yy11, trajMPC_beta12, trajMPC_bmy12);
trajMPC_LA_DENSE_FORWARDSUB_3(trajMPC_Ld12, trajMPC_bmy12, trajMPC_yy12);
trajMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(trajMPC_Ld12, trajMPC_Ysd13, trajMPC_Lsd13);
trajMPC_LA_DENSE_MMTSUB_3_3(trajMPC_Lsd13, trajMPC_Yd13);
trajMPC_LA_DENSE_CHOL_3(trajMPC_Yd13, trajMPC_Ld13);
trajMPC_LA_DENSE_MVMSUB1_3_3(trajMPC_Lsd13, trajMPC_yy12, trajMPC_beta13, trajMPC_bmy13);
trajMPC_LA_DENSE_FORWARDSUB_3(trajMPC_Ld13, trajMPC_bmy13, trajMPC_yy13);
trajMPC_LA_DENSE_MATRIXTFORWARDSUB_3_3(trajMPC_Ld13, trajMPC_Ysd14, trajMPC_Lsd14);
trajMPC_LA_DENSE_MMTSUB_3_3(trajMPC_Lsd14, trajMPC_Yd14);
trajMPC_LA_DENSE_CHOL_3(trajMPC_Yd14, trajMPC_Ld14);
trajMPC_LA_DENSE_MVMSUB1_3_3(trajMPC_Lsd14, trajMPC_yy13, trajMPC_beta14, trajMPC_bmy14);
trajMPC_LA_DENSE_FORWARDSUB_3(trajMPC_Ld14, trajMPC_bmy14, trajMPC_yy14);
trajMPC_LA_DENSE_BACKWARDSUB_3(trajMPC_Ld14, trajMPC_yy14, trajMPC_dvaff14);
trajMPC_LA_DENSE_MTVMSUB_3_3(trajMPC_Lsd14, trajMPC_dvaff14, trajMPC_yy13, trajMPC_bmy13);
trajMPC_LA_DENSE_BACKWARDSUB_3(trajMPC_Ld13, trajMPC_bmy13, trajMPC_dvaff13);
trajMPC_LA_DENSE_MTVMSUB_3_3(trajMPC_Lsd13, trajMPC_dvaff13, trajMPC_yy12, trajMPC_bmy12);
trajMPC_LA_DENSE_BACKWARDSUB_3(trajMPC_Ld12, trajMPC_bmy12, trajMPC_dvaff12);
trajMPC_LA_DENSE_MTVMSUB_3_3(trajMPC_Lsd12, trajMPC_dvaff12, trajMPC_yy11, trajMPC_bmy11);
trajMPC_LA_DENSE_BACKWARDSUB_3(trajMPC_Ld11, trajMPC_bmy11, trajMPC_dvaff11);
trajMPC_LA_DENSE_MTVMSUB_3_3(trajMPC_Lsd11, trajMPC_dvaff11, trajMPC_yy10, trajMPC_bmy10);
trajMPC_LA_DENSE_BACKWARDSUB_3(trajMPC_Ld10, trajMPC_bmy10, trajMPC_dvaff10);
trajMPC_LA_DENSE_MTVMSUB_3_3(trajMPC_Lsd10, trajMPC_dvaff10, trajMPC_yy09, trajMPC_bmy09);
trajMPC_LA_DENSE_BACKWARDSUB_3(trajMPC_Ld09, trajMPC_bmy09, trajMPC_dvaff09);
trajMPC_LA_DENSE_MTVMSUB_3_3(trajMPC_Lsd09, trajMPC_dvaff09, trajMPC_yy08, trajMPC_bmy08);
trajMPC_LA_DENSE_BACKWARDSUB_3(trajMPC_Ld08, trajMPC_bmy08, trajMPC_dvaff08);
trajMPC_LA_DENSE_MTVMSUB_3_3(trajMPC_Lsd08, trajMPC_dvaff08, trajMPC_yy07, trajMPC_bmy07);
trajMPC_LA_DENSE_BACKWARDSUB_3(trajMPC_Ld07, trajMPC_bmy07, trajMPC_dvaff07);
trajMPC_LA_DENSE_MTVMSUB_3_3(trajMPC_Lsd07, trajMPC_dvaff07, trajMPC_yy06, trajMPC_bmy06);
trajMPC_LA_DENSE_BACKWARDSUB_3(trajMPC_Ld06, trajMPC_bmy06, trajMPC_dvaff06);
trajMPC_LA_DENSE_MTVMSUB_3_3(trajMPC_Lsd06, trajMPC_dvaff06, trajMPC_yy05, trajMPC_bmy05);
trajMPC_LA_DENSE_BACKWARDSUB_3(trajMPC_Ld05, trajMPC_bmy05, trajMPC_dvaff05);
trajMPC_LA_DENSE_MTVMSUB_3_3(trajMPC_Lsd05, trajMPC_dvaff05, trajMPC_yy04, trajMPC_bmy04);
trajMPC_LA_DENSE_BACKWARDSUB_3(trajMPC_Ld04, trajMPC_bmy04, trajMPC_dvaff04);
trajMPC_LA_DENSE_MTVMSUB_3_3(trajMPC_Lsd04, trajMPC_dvaff04, trajMPC_yy03, trajMPC_bmy03);
trajMPC_LA_DENSE_BACKWARDSUB_3(trajMPC_Ld03, trajMPC_bmy03, trajMPC_dvaff03);
trajMPC_LA_DENSE_MTVMSUB_3_3(trajMPC_Lsd03, trajMPC_dvaff03, trajMPC_yy02, trajMPC_bmy02);
trajMPC_LA_DENSE_BACKWARDSUB_3(trajMPC_Ld02, trajMPC_bmy02, trajMPC_dvaff02);
trajMPC_LA_DENSE_MTVMSUB_3_3(trajMPC_Lsd02, trajMPC_dvaff02, trajMPC_yy01, trajMPC_bmy01);
trajMPC_LA_DENSE_BACKWARDSUB_3(trajMPC_Ld01, trajMPC_bmy01, trajMPC_dvaff01);
trajMPC_LA_DENSE_MTVMSUB_3_3(trajMPC_Lsd01, trajMPC_dvaff01, trajMPC_yy00, trajMPC_bmy00);
trajMPC_LA_DENSE_BACKWARDSUB_3(trajMPC_Ld00, trajMPC_bmy00, trajMPC_dvaff00);
trajMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C1, trajMPC_dvaff01, trajMPC_D00, trajMPC_dvaff00, trajMPC_grad_eq00);
trajMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C2, trajMPC_dvaff02, trajMPC_D01, trajMPC_dvaff01, trajMPC_grad_eq01);
trajMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C3, trajMPC_dvaff03, trajMPC_D01, trajMPC_dvaff02, trajMPC_grad_eq02);
trajMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C4, trajMPC_dvaff04, trajMPC_D01, trajMPC_dvaff03, trajMPC_grad_eq03);
trajMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C5, trajMPC_dvaff05, trajMPC_D01, trajMPC_dvaff04, trajMPC_grad_eq04);
trajMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C6, trajMPC_dvaff06, trajMPC_D01, trajMPC_dvaff05, trajMPC_grad_eq05);
trajMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C7, trajMPC_dvaff07, trajMPC_D01, trajMPC_dvaff06, trajMPC_grad_eq06);
trajMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C8, trajMPC_dvaff08, trajMPC_D01, trajMPC_dvaff07, trajMPC_grad_eq07);
trajMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C9, trajMPC_dvaff09, trajMPC_D01, trajMPC_dvaff08, trajMPC_grad_eq08);
trajMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C10, trajMPC_dvaff10, trajMPC_D01, trajMPC_dvaff09, trajMPC_grad_eq09);
trajMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C11, trajMPC_dvaff11, trajMPC_D01, trajMPC_dvaff10, trajMPC_grad_eq10);
trajMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C12, trajMPC_dvaff12, trajMPC_D01, trajMPC_dvaff11, trajMPC_grad_eq11);
trajMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C13, trajMPC_dvaff13, trajMPC_D01, trajMPC_dvaff12, trajMPC_grad_eq12);
trajMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C14, trajMPC_dvaff14, trajMPC_D01, trajMPC_dvaff13, trajMPC_grad_eq13);
trajMPC_LA_DIAGZERO_MTVM_3_3(trajMPC_D14, trajMPC_dvaff14, trajMPC_grad_eq14);
trajMPC_LA_VSUB2_157(trajMPC_rd, trajMPC_grad_eq, trajMPC_rd);
trajMPC_LA_DIAG_FORWARDBACKWARDSUB_11(trajMPC_Phi00, trajMPC_rd00, trajMPC_dzaff00);
trajMPC_LA_DIAG_FORWARDBACKWARDSUB_11(trajMPC_Phi01, trajMPC_rd01, trajMPC_dzaff01);
trajMPC_LA_DIAG_FORWARDBACKWARDSUB_11(trajMPC_Phi02, trajMPC_rd02, trajMPC_dzaff02);
trajMPC_LA_DIAG_FORWARDBACKWARDSUB_11(trajMPC_Phi03, trajMPC_rd03, trajMPC_dzaff03);
trajMPC_LA_DIAG_FORWARDBACKWARDSUB_11(trajMPC_Phi04, trajMPC_rd04, trajMPC_dzaff04);
trajMPC_LA_DIAG_FORWARDBACKWARDSUB_11(trajMPC_Phi05, trajMPC_rd05, trajMPC_dzaff05);
trajMPC_LA_DIAG_FORWARDBACKWARDSUB_11(trajMPC_Phi06, trajMPC_rd06, trajMPC_dzaff06);
trajMPC_LA_DIAG_FORWARDBACKWARDSUB_11(trajMPC_Phi07, trajMPC_rd07, trajMPC_dzaff07);
trajMPC_LA_DIAG_FORWARDBACKWARDSUB_11(trajMPC_Phi08, trajMPC_rd08, trajMPC_dzaff08);
trajMPC_LA_DIAG_FORWARDBACKWARDSUB_11(trajMPC_Phi09, trajMPC_rd09, trajMPC_dzaff09);
trajMPC_LA_DIAG_FORWARDBACKWARDSUB_11(trajMPC_Phi10, trajMPC_rd10, trajMPC_dzaff10);
trajMPC_LA_DIAG_FORWARDBACKWARDSUB_11(trajMPC_Phi11, trajMPC_rd11, trajMPC_dzaff11);
trajMPC_LA_DIAG_FORWARDBACKWARDSUB_11(trajMPC_Phi12, trajMPC_rd12, trajMPC_dzaff12);
trajMPC_LA_DIAG_FORWARDBACKWARDSUB_11(trajMPC_Phi13, trajMPC_rd13, trajMPC_dzaff13);
trajMPC_LA_DIAG_FORWARDBACKWARDSUB_3(trajMPC_Phi14, trajMPC_rd14, trajMPC_dzaff14);
trajMPC_LA_VSUB_INDEXED_11(trajMPC_dzaff00, trajMPC_lbIdx00, trajMPC_rilb00, trajMPC_dslbaff00);
trajMPC_LA_VSUB3_11(trajMPC_llbbyslb00, trajMPC_dslbaff00, trajMPC_llb00, trajMPC_dllbaff00);
trajMPC_LA_VSUB2_INDEXED_5(trajMPC_riub00, trajMPC_dzaff00, trajMPC_ubIdx00, trajMPC_dsubaff00);
trajMPC_LA_VSUB3_5(trajMPC_lubbysub00, trajMPC_dsubaff00, trajMPC_lub00, trajMPC_dlubaff00);
trajMPC_LA_VSUB_INDEXED_11(trajMPC_dzaff01, trajMPC_lbIdx01, trajMPC_rilb01, trajMPC_dslbaff01);
trajMPC_LA_VSUB3_11(trajMPC_llbbyslb01, trajMPC_dslbaff01, trajMPC_llb01, trajMPC_dllbaff01);
trajMPC_LA_VSUB2_INDEXED_5(trajMPC_riub01, trajMPC_dzaff01, trajMPC_ubIdx01, trajMPC_dsubaff01);
trajMPC_LA_VSUB3_5(trajMPC_lubbysub01, trajMPC_dsubaff01, trajMPC_lub01, trajMPC_dlubaff01);
trajMPC_LA_VSUB_INDEXED_11(trajMPC_dzaff02, trajMPC_lbIdx02, trajMPC_rilb02, trajMPC_dslbaff02);
trajMPC_LA_VSUB3_11(trajMPC_llbbyslb02, trajMPC_dslbaff02, trajMPC_llb02, trajMPC_dllbaff02);
trajMPC_LA_VSUB2_INDEXED_5(trajMPC_riub02, trajMPC_dzaff02, trajMPC_ubIdx02, trajMPC_dsubaff02);
trajMPC_LA_VSUB3_5(trajMPC_lubbysub02, trajMPC_dsubaff02, trajMPC_lub02, trajMPC_dlubaff02);
trajMPC_LA_VSUB_INDEXED_11(trajMPC_dzaff03, trajMPC_lbIdx03, trajMPC_rilb03, trajMPC_dslbaff03);
trajMPC_LA_VSUB3_11(trajMPC_llbbyslb03, trajMPC_dslbaff03, trajMPC_llb03, trajMPC_dllbaff03);
trajMPC_LA_VSUB2_INDEXED_5(trajMPC_riub03, trajMPC_dzaff03, trajMPC_ubIdx03, trajMPC_dsubaff03);
trajMPC_LA_VSUB3_5(trajMPC_lubbysub03, trajMPC_dsubaff03, trajMPC_lub03, trajMPC_dlubaff03);
trajMPC_LA_VSUB_INDEXED_11(trajMPC_dzaff04, trajMPC_lbIdx04, trajMPC_rilb04, trajMPC_dslbaff04);
trajMPC_LA_VSUB3_11(trajMPC_llbbyslb04, trajMPC_dslbaff04, trajMPC_llb04, trajMPC_dllbaff04);
trajMPC_LA_VSUB2_INDEXED_5(trajMPC_riub04, trajMPC_dzaff04, trajMPC_ubIdx04, trajMPC_dsubaff04);
trajMPC_LA_VSUB3_5(trajMPC_lubbysub04, trajMPC_dsubaff04, trajMPC_lub04, trajMPC_dlubaff04);
trajMPC_LA_VSUB_INDEXED_11(trajMPC_dzaff05, trajMPC_lbIdx05, trajMPC_rilb05, trajMPC_dslbaff05);
trajMPC_LA_VSUB3_11(trajMPC_llbbyslb05, trajMPC_dslbaff05, trajMPC_llb05, trajMPC_dllbaff05);
trajMPC_LA_VSUB2_INDEXED_5(trajMPC_riub05, trajMPC_dzaff05, trajMPC_ubIdx05, trajMPC_dsubaff05);
trajMPC_LA_VSUB3_5(trajMPC_lubbysub05, trajMPC_dsubaff05, trajMPC_lub05, trajMPC_dlubaff05);
trajMPC_LA_VSUB_INDEXED_11(trajMPC_dzaff06, trajMPC_lbIdx06, trajMPC_rilb06, trajMPC_dslbaff06);
trajMPC_LA_VSUB3_11(trajMPC_llbbyslb06, trajMPC_dslbaff06, trajMPC_llb06, trajMPC_dllbaff06);
trajMPC_LA_VSUB2_INDEXED_5(trajMPC_riub06, trajMPC_dzaff06, trajMPC_ubIdx06, trajMPC_dsubaff06);
trajMPC_LA_VSUB3_5(trajMPC_lubbysub06, trajMPC_dsubaff06, trajMPC_lub06, trajMPC_dlubaff06);
trajMPC_LA_VSUB_INDEXED_11(trajMPC_dzaff07, trajMPC_lbIdx07, trajMPC_rilb07, trajMPC_dslbaff07);
trajMPC_LA_VSUB3_11(trajMPC_llbbyslb07, trajMPC_dslbaff07, trajMPC_llb07, trajMPC_dllbaff07);
trajMPC_LA_VSUB2_INDEXED_5(trajMPC_riub07, trajMPC_dzaff07, trajMPC_ubIdx07, trajMPC_dsubaff07);
trajMPC_LA_VSUB3_5(trajMPC_lubbysub07, trajMPC_dsubaff07, trajMPC_lub07, trajMPC_dlubaff07);
trajMPC_LA_VSUB_INDEXED_11(trajMPC_dzaff08, trajMPC_lbIdx08, trajMPC_rilb08, trajMPC_dslbaff08);
trajMPC_LA_VSUB3_11(trajMPC_llbbyslb08, trajMPC_dslbaff08, trajMPC_llb08, trajMPC_dllbaff08);
trajMPC_LA_VSUB2_INDEXED_5(trajMPC_riub08, trajMPC_dzaff08, trajMPC_ubIdx08, trajMPC_dsubaff08);
trajMPC_LA_VSUB3_5(trajMPC_lubbysub08, trajMPC_dsubaff08, trajMPC_lub08, trajMPC_dlubaff08);
trajMPC_LA_VSUB_INDEXED_11(trajMPC_dzaff09, trajMPC_lbIdx09, trajMPC_rilb09, trajMPC_dslbaff09);
trajMPC_LA_VSUB3_11(trajMPC_llbbyslb09, trajMPC_dslbaff09, trajMPC_llb09, trajMPC_dllbaff09);
trajMPC_LA_VSUB2_INDEXED_5(trajMPC_riub09, trajMPC_dzaff09, trajMPC_ubIdx09, trajMPC_dsubaff09);
trajMPC_LA_VSUB3_5(trajMPC_lubbysub09, trajMPC_dsubaff09, trajMPC_lub09, trajMPC_dlubaff09);
trajMPC_LA_VSUB_INDEXED_11(trajMPC_dzaff10, trajMPC_lbIdx10, trajMPC_rilb10, trajMPC_dslbaff10);
trajMPC_LA_VSUB3_11(trajMPC_llbbyslb10, trajMPC_dslbaff10, trajMPC_llb10, trajMPC_dllbaff10);
trajMPC_LA_VSUB2_INDEXED_5(trajMPC_riub10, trajMPC_dzaff10, trajMPC_ubIdx10, trajMPC_dsubaff10);
trajMPC_LA_VSUB3_5(trajMPC_lubbysub10, trajMPC_dsubaff10, trajMPC_lub10, trajMPC_dlubaff10);
trajMPC_LA_VSUB_INDEXED_11(trajMPC_dzaff11, trajMPC_lbIdx11, trajMPC_rilb11, trajMPC_dslbaff11);
trajMPC_LA_VSUB3_11(trajMPC_llbbyslb11, trajMPC_dslbaff11, trajMPC_llb11, trajMPC_dllbaff11);
trajMPC_LA_VSUB2_INDEXED_5(trajMPC_riub11, trajMPC_dzaff11, trajMPC_ubIdx11, trajMPC_dsubaff11);
trajMPC_LA_VSUB3_5(trajMPC_lubbysub11, trajMPC_dsubaff11, trajMPC_lub11, trajMPC_dlubaff11);
trajMPC_LA_VSUB_INDEXED_11(trajMPC_dzaff12, trajMPC_lbIdx12, trajMPC_rilb12, trajMPC_dslbaff12);
trajMPC_LA_VSUB3_11(trajMPC_llbbyslb12, trajMPC_dslbaff12, trajMPC_llb12, trajMPC_dllbaff12);
trajMPC_LA_VSUB2_INDEXED_5(trajMPC_riub12, trajMPC_dzaff12, trajMPC_ubIdx12, trajMPC_dsubaff12);
trajMPC_LA_VSUB3_5(trajMPC_lubbysub12, trajMPC_dsubaff12, trajMPC_lub12, trajMPC_dlubaff12);
trajMPC_LA_VSUB_INDEXED_11(trajMPC_dzaff13, trajMPC_lbIdx13, trajMPC_rilb13, trajMPC_dslbaff13);
trajMPC_LA_VSUB3_11(trajMPC_llbbyslb13, trajMPC_dslbaff13, trajMPC_llb13, trajMPC_dllbaff13);
trajMPC_LA_VSUB2_INDEXED_5(trajMPC_riub13, trajMPC_dzaff13, trajMPC_ubIdx13, trajMPC_dsubaff13);
trajMPC_LA_VSUB3_5(trajMPC_lubbysub13, trajMPC_dsubaff13, trajMPC_lub13, trajMPC_dlubaff13);
trajMPC_LA_VSUB_INDEXED_3(trajMPC_dzaff14, trajMPC_lbIdx14, trajMPC_rilb14, trajMPC_dslbaff14);
trajMPC_LA_VSUB3_3(trajMPC_llbbyslb14, trajMPC_dslbaff14, trajMPC_llb14, trajMPC_dllbaff14);
trajMPC_LA_VSUB2_INDEXED_3(trajMPC_riub14, trajMPC_dzaff14, trajMPC_ubIdx14, trajMPC_dsubaff14);
trajMPC_LA_VSUB3_3(trajMPC_lubbysub14, trajMPC_dsubaff14, trajMPC_lub14, trajMPC_dlubaff14);
info->lsit_aff = trajMPC_LINESEARCH_BACKTRACKING_AFFINE(trajMPC_l, trajMPC_s, trajMPC_dl_aff, trajMPC_ds_aff, &info->step_aff, &info->mu_aff);
if( info->lsit_aff == trajMPC_NOPROGRESS ){
exitcode = trajMPC_NOPROGRESS; break;
}
sigma_3rdroot = info->mu_aff / info->mu;
info->sigma = sigma_3rdroot*sigma_3rdroot*sigma_3rdroot;
musigma = info->mu * info->sigma;
trajMPC_LA_VSUB5_230(trajMPC_ds_aff, trajMPC_dl_aff, musigma, trajMPC_ccrhs);
trajMPC_LA_VSUB6_INDEXED_11_5_11(trajMPC_ccrhsub00, trajMPC_sub00, trajMPC_ubIdx00, trajMPC_ccrhsl00, trajMPC_slb00, trajMPC_lbIdx00, trajMPC_rd00);
trajMPC_LA_VSUB6_INDEXED_11_5_11(trajMPC_ccrhsub01, trajMPC_sub01, trajMPC_ubIdx01, trajMPC_ccrhsl01, trajMPC_slb01, trajMPC_lbIdx01, trajMPC_rd01);
trajMPC_LA_DIAG_FORWARDSUB_11(trajMPC_Phi00, trajMPC_rd00, trajMPC_Lbyrd00);
trajMPC_LA_DIAG_FORWARDSUB_11(trajMPC_Phi01, trajMPC_rd01, trajMPC_Lbyrd01);
trajMPC_LA_DIAGZERO_MVM_3(trajMPC_W00, trajMPC_Lbyrd00, trajMPC_beta00);
trajMPC_LA_DENSE_FORWARDSUB_3(trajMPC_Ld00, trajMPC_beta00, trajMPC_yy00);
trajMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(trajMPC_V00, trajMPC_Lbyrd00, trajMPC_W01, trajMPC_Lbyrd01, trajMPC_beta01);
trajMPC_LA_DENSE_MVMSUB1_3_3(trajMPC_Lsd01, trajMPC_yy00, trajMPC_beta01, trajMPC_bmy01);
trajMPC_LA_DENSE_FORWARDSUB_3(trajMPC_Ld01, trajMPC_bmy01, trajMPC_yy01);
trajMPC_LA_VSUB6_INDEXED_11_5_11(trajMPC_ccrhsub02, trajMPC_sub02, trajMPC_ubIdx02, trajMPC_ccrhsl02, trajMPC_slb02, trajMPC_lbIdx02, trajMPC_rd02);
trajMPC_LA_DIAG_FORWARDSUB_11(trajMPC_Phi02, trajMPC_rd02, trajMPC_Lbyrd02);
trajMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(trajMPC_V01, trajMPC_Lbyrd01, trajMPC_W02, trajMPC_Lbyrd02, trajMPC_beta02);
trajMPC_LA_DENSE_MVMSUB1_3_3(trajMPC_Lsd02, trajMPC_yy01, trajMPC_beta02, trajMPC_bmy02);
trajMPC_LA_DENSE_FORWARDSUB_3(trajMPC_Ld02, trajMPC_bmy02, trajMPC_yy02);
trajMPC_LA_VSUB6_INDEXED_11_5_11(trajMPC_ccrhsub03, trajMPC_sub03, trajMPC_ubIdx03, trajMPC_ccrhsl03, trajMPC_slb03, trajMPC_lbIdx03, trajMPC_rd03);
trajMPC_LA_DIAG_FORWARDSUB_11(trajMPC_Phi03, trajMPC_rd03, trajMPC_Lbyrd03);
trajMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(trajMPC_V02, trajMPC_Lbyrd02, trajMPC_W03, trajMPC_Lbyrd03, trajMPC_beta03);
trajMPC_LA_DENSE_MVMSUB1_3_3(trajMPC_Lsd03, trajMPC_yy02, trajMPC_beta03, trajMPC_bmy03);
trajMPC_LA_DENSE_FORWARDSUB_3(trajMPC_Ld03, trajMPC_bmy03, trajMPC_yy03);
trajMPC_LA_VSUB6_INDEXED_11_5_11(trajMPC_ccrhsub04, trajMPC_sub04, trajMPC_ubIdx04, trajMPC_ccrhsl04, trajMPC_slb04, trajMPC_lbIdx04, trajMPC_rd04);
trajMPC_LA_DIAG_FORWARDSUB_11(trajMPC_Phi04, trajMPC_rd04, trajMPC_Lbyrd04);
trajMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(trajMPC_V03, trajMPC_Lbyrd03, trajMPC_W04, trajMPC_Lbyrd04, trajMPC_beta04);
trajMPC_LA_DENSE_MVMSUB1_3_3(trajMPC_Lsd04, trajMPC_yy03, trajMPC_beta04, trajMPC_bmy04);
trajMPC_LA_DENSE_FORWARDSUB_3(trajMPC_Ld04, trajMPC_bmy04, trajMPC_yy04);
trajMPC_LA_VSUB6_INDEXED_11_5_11(trajMPC_ccrhsub05, trajMPC_sub05, trajMPC_ubIdx05, trajMPC_ccrhsl05, trajMPC_slb05, trajMPC_lbIdx05, trajMPC_rd05);
trajMPC_LA_DIAG_FORWARDSUB_11(trajMPC_Phi05, trajMPC_rd05, trajMPC_Lbyrd05);
trajMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(trajMPC_V04, trajMPC_Lbyrd04, trajMPC_W05, trajMPC_Lbyrd05, trajMPC_beta05);
trajMPC_LA_DENSE_MVMSUB1_3_3(trajMPC_Lsd05, trajMPC_yy04, trajMPC_beta05, trajMPC_bmy05);
trajMPC_LA_DENSE_FORWARDSUB_3(trajMPC_Ld05, trajMPC_bmy05, trajMPC_yy05);
trajMPC_LA_VSUB6_INDEXED_11_5_11(trajMPC_ccrhsub06, trajMPC_sub06, trajMPC_ubIdx06, trajMPC_ccrhsl06, trajMPC_slb06, trajMPC_lbIdx06, trajMPC_rd06);
trajMPC_LA_DIAG_FORWARDSUB_11(trajMPC_Phi06, trajMPC_rd06, trajMPC_Lbyrd06);
trajMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(trajMPC_V05, trajMPC_Lbyrd05, trajMPC_W06, trajMPC_Lbyrd06, trajMPC_beta06);
trajMPC_LA_DENSE_MVMSUB1_3_3(trajMPC_Lsd06, trajMPC_yy05, trajMPC_beta06, trajMPC_bmy06);
trajMPC_LA_DENSE_FORWARDSUB_3(trajMPC_Ld06, trajMPC_bmy06, trajMPC_yy06);
trajMPC_LA_VSUB6_INDEXED_11_5_11(trajMPC_ccrhsub07, trajMPC_sub07, trajMPC_ubIdx07, trajMPC_ccrhsl07, trajMPC_slb07, trajMPC_lbIdx07, trajMPC_rd07);
trajMPC_LA_DIAG_FORWARDSUB_11(trajMPC_Phi07, trajMPC_rd07, trajMPC_Lbyrd07);
trajMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(trajMPC_V06, trajMPC_Lbyrd06, trajMPC_W07, trajMPC_Lbyrd07, trajMPC_beta07);
trajMPC_LA_DENSE_MVMSUB1_3_3(trajMPC_Lsd07, trajMPC_yy06, trajMPC_beta07, trajMPC_bmy07);
trajMPC_LA_DENSE_FORWARDSUB_3(trajMPC_Ld07, trajMPC_bmy07, trajMPC_yy07);
trajMPC_LA_VSUB6_INDEXED_11_5_11(trajMPC_ccrhsub08, trajMPC_sub08, trajMPC_ubIdx08, trajMPC_ccrhsl08, trajMPC_slb08, trajMPC_lbIdx08, trajMPC_rd08);
trajMPC_LA_DIAG_FORWARDSUB_11(trajMPC_Phi08, trajMPC_rd08, trajMPC_Lbyrd08);
trajMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(trajMPC_V07, trajMPC_Lbyrd07, trajMPC_W08, trajMPC_Lbyrd08, trajMPC_beta08);
trajMPC_LA_DENSE_MVMSUB1_3_3(trajMPC_Lsd08, trajMPC_yy07, trajMPC_beta08, trajMPC_bmy08);
trajMPC_LA_DENSE_FORWARDSUB_3(trajMPC_Ld08, trajMPC_bmy08, trajMPC_yy08);
trajMPC_LA_VSUB6_INDEXED_11_5_11(trajMPC_ccrhsub09, trajMPC_sub09, trajMPC_ubIdx09, trajMPC_ccrhsl09, trajMPC_slb09, trajMPC_lbIdx09, trajMPC_rd09);
trajMPC_LA_DIAG_FORWARDSUB_11(trajMPC_Phi09, trajMPC_rd09, trajMPC_Lbyrd09);
trajMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(trajMPC_V08, trajMPC_Lbyrd08, trajMPC_W09, trajMPC_Lbyrd09, trajMPC_beta09);
trajMPC_LA_DENSE_MVMSUB1_3_3(trajMPC_Lsd09, trajMPC_yy08, trajMPC_beta09, trajMPC_bmy09);
trajMPC_LA_DENSE_FORWARDSUB_3(trajMPC_Ld09, trajMPC_bmy09, trajMPC_yy09);
trajMPC_LA_VSUB6_INDEXED_11_5_11(trajMPC_ccrhsub10, trajMPC_sub10, trajMPC_ubIdx10, trajMPC_ccrhsl10, trajMPC_slb10, trajMPC_lbIdx10, trajMPC_rd10);
trajMPC_LA_DIAG_FORWARDSUB_11(trajMPC_Phi10, trajMPC_rd10, trajMPC_Lbyrd10);
trajMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(trajMPC_V09, trajMPC_Lbyrd09, trajMPC_W10, trajMPC_Lbyrd10, trajMPC_beta10);
trajMPC_LA_DENSE_MVMSUB1_3_3(trajMPC_Lsd10, trajMPC_yy09, trajMPC_beta10, trajMPC_bmy10);
trajMPC_LA_DENSE_FORWARDSUB_3(trajMPC_Ld10, trajMPC_bmy10, trajMPC_yy10);
trajMPC_LA_VSUB6_INDEXED_11_5_11(trajMPC_ccrhsub11, trajMPC_sub11, trajMPC_ubIdx11, trajMPC_ccrhsl11, trajMPC_slb11, trajMPC_lbIdx11, trajMPC_rd11);
trajMPC_LA_DIAG_FORWARDSUB_11(trajMPC_Phi11, trajMPC_rd11, trajMPC_Lbyrd11);
trajMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(trajMPC_V10, trajMPC_Lbyrd10, trajMPC_W11, trajMPC_Lbyrd11, trajMPC_beta11);
trajMPC_LA_DENSE_MVMSUB1_3_3(trajMPC_Lsd11, trajMPC_yy10, trajMPC_beta11, trajMPC_bmy11);
trajMPC_LA_DENSE_FORWARDSUB_3(trajMPC_Ld11, trajMPC_bmy11, trajMPC_yy11);
trajMPC_LA_VSUB6_INDEXED_11_5_11(trajMPC_ccrhsub12, trajMPC_sub12, trajMPC_ubIdx12, trajMPC_ccrhsl12, trajMPC_slb12, trajMPC_lbIdx12, trajMPC_rd12);
trajMPC_LA_DIAG_FORWARDSUB_11(trajMPC_Phi12, trajMPC_rd12, trajMPC_Lbyrd12);
trajMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(trajMPC_V11, trajMPC_Lbyrd11, trajMPC_W12, trajMPC_Lbyrd12, trajMPC_beta12);
trajMPC_LA_DENSE_MVMSUB1_3_3(trajMPC_Lsd12, trajMPC_yy11, trajMPC_beta12, trajMPC_bmy12);
trajMPC_LA_DENSE_FORWARDSUB_3(trajMPC_Ld12, trajMPC_bmy12, trajMPC_yy12);
trajMPC_LA_VSUB6_INDEXED_11_5_11(trajMPC_ccrhsub13, trajMPC_sub13, trajMPC_ubIdx13, trajMPC_ccrhsl13, trajMPC_slb13, trajMPC_lbIdx13, trajMPC_rd13);
trajMPC_LA_DIAG_FORWARDSUB_11(trajMPC_Phi13, trajMPC_rd13, trajMPC_Lbyrd13);
trajMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_11(trajMPC_V12, trajMPC_Lbyrd12, trajMPC_W13, trajMPC_Lbyrd13, trajMPC_beta13);
trajMPC_LA_DENSE_MVMSUB1_3_3(trajMPC_Lsd13, trajMPC_yy12, trajMPC_beta13, trajMPC_bmy13);
trajMPC_LA_DENSE_FORWARDSUB_3(trajMPC_Ld13, trajMPC_bmy13, trajMPC_yy13);
trajMPC_LA_VSUB6_INDEXED_3_3_3(trajMPC_ccrhsub14, trajMPC_sub14, trajMPC_ubIdx14, trajMPC_ccrhsl14, trajMPC_slb14, trajMPC_lbIdx14, trajMPC_rd14);
trajMPC_LA_DIAG_FORWARDSUB_3(trajMPC_Phi14, trajMPC_rd14, trajMPC_Lbyrd14);
trajMPC_LA_DENSE_DIAGZERO_2MVMADD_3_11_3(trajMPC_V13, trajMPC_Lbyrd13, trajMPC_W14, trajMPC_Lbyrd14, trajMPC_beta14);
trajMPC_LA_DENSE_MVMSUB1_3_3(trajMPC_Lsd14, trajMPC_yy13, trajMPC_beta14, trajMPC_bmy14);
trajMPC_LA_DENSE_FORWARDSUB_3(trajMPC_Ld14, trajMPC_bmy14, trajMPC_yy14);
trajMPC_LA_DENSE_BACKWARDSUB_3(trajMPC_Ld14, trajMPC_yy14, trajMPC_dvcc14);
trajMPC_LA_DENSE_MTVMSUB_3_3(trajMPC_Lsd14, trajMPC_dvcc14, trajMPC_yy13, trajMPC_bmy13);
trajMPC_LA_DENSE_BACKWARDSUB_3(trajMPC_Ld13, trajMPC_bmy13, trajMPC_dvcc13);
trajMPC_LA_DENSE_MTVMSUB_3_3(trajMPC_Lsd13, trajMPC_dvcc13, trajMPC_yy12, trajMPC_bmy12);
trajMPC_LA_DENSE_BACKWARDSUB_3(trajMPC_Ld12, trajMPC_bmy12, trajMPC_dvcc12);
trajMPC_LA_DENSE_MTVMSUB_3_3(trajMPC_Lsd12, trajMPC_dvcc12, trajMPC_yy11, trajMPC_bmy11);
trajMPC_LA_DENSE_BACKWARDSUB_3(trajMPC_Ld11, trajMPC_bmy11, trajMPC_dvcc11);
trajMPC_LA_DENSE_MTVMSUB_3_3(trajMPC_Lsd11, trajMPC_dvcc11, trajMPC_yy10, trajMPC_bmy10);
trajMPC_LA_DENSE_BACKWARDSUB_3(trajMPC_Ld10, trajMPC_bmy10, trajMPC_dvcc10);
trajMPC_LA_DENSE_MTVMSUB_3_3(trajMPC_Lsd10, trajMPC_dvcc10, trajMPC_yy09, trajMPC_bmy09);
trajMPC_LA_DENSE_BACKWARDSUB_3(trajMPC_Ld09, trajMPC_bmy09, trajMPC_dvcc09);
trajMPC_LA_DENSE_MTVMSUB_3_3(trajMPC_Lsd09, trajMPC_dvcc09, trajMPC_yy08, trajMPC_bmy08);
trajMPC_LA_DENSE_BACKWARDSUB_3(trajMPC_Ld08, trajMPC_bmy08, trajMPC_dvcc08);
trajMPC_LA_DENSE_MTVMSUB_3_3(trajMPC_Lsd08, trajMPC_dvcc08, trajMPC_yy07, trajMPC_bmy07);
trajMPC_LA_DENSE_BACKWARDSUB_3(trajMPC_Ld07, trajMPC_bmy07, trajMPC_dvcc07);
trajMPC_LA_DENSE_MTVMSUB_3_3(trajMPC_Lsd07, trajMPC_dvcc07, trajMPC_yy06, trajMPC_bmy06);
trajMPC_LA_DENSE_BACKWARDSUB_3(trajMPC_Ld06, trajMPC_bmy06, trajMPC_dvcc06);
trajMPC_LA_DENSE_MTVMSUB_3_3(trajMPC_Lsd06, trajMPC_dvcc06, trajMPC_yy05, trajMPC_bmy05);
trajMPC_LA_DENSE_BACKWARDSUB_3(trajMPC_Ld05, trajMPC_bmy05, trajMPC_dvcc05);
trajMPC_LA_DENSE_MTVMSUB_3_3(trajMPC_Lsd05, trajMPC_dvcc05, trajMPC_yy04, trajMPC_bmy04);
trajMPC_LA_DENSE_BACKWARDSUB_3(trajMPC_Ld04, trajMPC_bmy04, trajMPC_dvcc04);
trajMPC_LA_DENSE_MTVMSUB_3_3(trajMPC_Lsd04, trajMPC_dvcc04, trajMPC_yy03, trajMPC_bmy03);
trajMPC_LA_DENSE_BACKWARDSUB_3(trajMPC_Ld03, trajMPC_bmy03, trajMPC_dvcc03);
trajMPC_LA_DENSE_MTVMSUB_3_3(trajMPC_Lsd03, trajMPC_dvcc03, trajMPC_yy02, trajMPC_bmy02);
trajMPC_LA_DENSE_BACKWARDSUB_3(trajMPC_Ld02, trajMPC_bmy02, trajMPC_dvcc02);
trajMPC_LA_DENSE_MTVMSUB_3_3(trajMPC_Lsd02, trajMPC_dvcc02, trajMPC_yy01, trajMPC_bmy01);
trajMPC_LA_DENSE_BACKWARDSUB_3(trajMPC_Ld01, trajMPC_bmy01, trajMPC_dvcc01);
trajMPC_LA_DENSE_MTVMSUB_3_3(trajMPC_Lsd01, trajMPC_dvcc01, trajMPC_yy00, trajMPC_bmy00);
trajMPC_LA_DENSE_BACKWARDSUB_3(trajMPC_Ld00, trajMPC_bmy00, trajMPC_dvcc00);
trajMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C1, trajMPC_dvcc01, trajMPC_D00, trajMPC_dvcc00, trajMPC_grad_eq00);
trajMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C2, trajMPC_dvcc02, trajMPC_D01, trajMPC_dvcc01, trajMPC_grad_eq01);
trajMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C3, trajMPC_dvcc03, trajMPC_D01, trajMPC_dvcc02, trajMPC_grad_eq02);
trajMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C4, trajMPC_dvcc04, trajMPC_D01, trajMPC_dvcc03, trajMPC_grad_eq03);
trajMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C5, trajMPC_dvcc05, trajMPC_D01, trajMPC_dvcc04, trajMPC_grad_eq04);
trajMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C6, trajMPC_dvcc06, trajMPC_D01, trajMPC_dvcc05, trajMPC_grad_eq05);
trajMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C7, trajMPC_dvcc07, trajMPC_D01, trajMPC_dvcc06, trajMPC_grad_eq06);
trajMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C8, trajMPC_dvcc08, trajMPC_D01, trajMPC_dvcc07, trajMPC_grad_eq07);
trajMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C9, trajMPC_dvcc09, trajMPC_D01, trajMPC_dvcc08, trajMPC_grad_eq08);
trajMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C10, trajMPC_dvcc10, trajMPC_D01, trajMPC_dvcc09, trajMPC_grad_eq09);
trajMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C11, trajMPC_dvcc11, trajMPC_D01, trajMPC_dvcc10, trajMPC_grad_eq10);
trajMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C12, trajMPC_dvcc12, trajMPC_D01, trajMPC_dvcc11, trajMPC_grad_eq11);
trajMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C13, trajMPC_dvcc13, trajMPC_D01, trajMPC_dvcc12, trajMPC_grad_eq12);
trajMPC_LA_DENSE_DIAGZERO_MTVM2_3_11_3(params->C14, trajMPC_dvcc14, trajMPC_D01, trajMPC_dvcc13, trajMPC_grad_eq13);
trajMPC_LA_DIAGZERO_MTVM_3_3(trajMPC_D14, trajMPC_dvcc14, trajMPC_grad_eq14);
trajMPC_LA_VSUB_157(trajMPC_rd, trajMPC_grad_eq, trajMPC_rd);
trajMPC_LA_DIAG_FORWARDBACKWARDSUB_11(trajMPC_Phi00, trajMPC_rd00, trajMPC_dzcc00);
trajMPC_LA_DIAG_FORWARDBACKWARDSUB_11(trajMPC_Phi01, trajMPC_rd01, trajMPC_dzcc01);
trajMPC_LA_DIAG_FORWARDBACKWARDSUB_11(trajMPC_Phi02, trajMPC_rd02, trajMPC_dzcc02);
trajMPC_LA_DIAG_FORWARDBACKWARDSUB_11(trajMPC_Phi03, trajMPC_rd03, trajMPC_dzcc03);
trajMPC_LA_DIAG_FORWARDBACKWARDSUB_11(trajMPC_Phi04, trajMPC_rd04, trajMPC_dzcc04);
trajMPC_LA_DIAG_FORWARDBACKWARDSUB_11(trajMPC_Phi05, trajMPC_rd05, trajMPC_dzcc05);
trajMPC_LA_DIAG_FORWARDBACKWARDSUB_11(trajMPC_Phi06, trajMPC_rd06, trajMPC_dzcc06);
trajMPC_LA_DIAG_FORWARDBACKWARDSUB_11(trajMPC_Phi07, trajMPC_rd07, trajMPC_dzcc07);
trajMPC_LA_DIAG_FORWARDBACKWARDSUB_11(trajMPC_Phi08, trajMPC_rd08, trajMPC_dzcc08);
trajMPC_LA_DIAG_FORWARDBACKWARDSUB_11(trajMPC_Phi09, trajMPC_rd09, trajMPC_dzcc09);
trajMPC_LA_DIAG_FORWARDBACKWARDSUB_11(trajMPC_Phi10, trajMPC_rd10, trajMPC_dzcc10);
trajMPC_LA_DIAG_FORWARDBACKWARDSUB_11(trajMPC_Phi11, trajMPC_rd11, trajMPC_dzcc11);
trajMPC_LA_DIAG_FORWARDBACKWARDSUB_11(trajMPC_Phi12, trajMPC_rd12, trajMPC_dzcc12);
trajMPC_LA_DIAG_FORWARDBACKWARDSUB_11(trajMPC_Phi13, trajMPC_rd13, trajMPC_dzcc13);
trajMPC_LA_DIAG_FORWARDBACKWARDSUB_3(trajMPC_Phi14, trajMPC_rd14, trajMPC_dzcc14);
trajMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(trajMPC_ccrhsl00, trajMPC_slb00, trajMPC_llbbyslb00, trajMPC_dzcc00, trajMPC_lbIdx00, trajMPC_dllbcc00);
trajMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(trajMPC_ccrhsub00, trajMPC_sub00, trajMPC_lubbysub00, trajMPC_dzcc00, trajMPC_ubIdx00, trajMPC_dlubcc00);
trajMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(trajMPC_ccrhsl01, trajMPC_slb01, trajMPC_llbbyslb01, trajMPC_dzcc01, trajMPC_lbIdx01, trajMPC_dllbcc01);
trajMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(trajMPC_ccrhsub01, trajMPC_sub01, trajMPC_lubbysub01, trajMPC_dzcc01, trajMPC_ubIdx01, trajMPC_dlubcc01);
trajMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(trajMPC_ccrhsl02, trajMPC_slb02, trajMPC_llbbyslb02, trajMPC_dzcc02, trajMPC_lbIdx02, trajMPC_dllbcc02);
trajMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(trajMPC_ccrhsub02, trajMPC_sub02, trajMPC_lubbysub02, trajMPC_dzcc02, trajMPC_ubIdx02, trajMPC_dlubcc02);
trajMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(trajMPC_ccrhsl03, trajMPC_slb03, trajMPC_llbbyslb03, trajMPC_dzcc03, trajMPC_lbIdx03, trajMPC_dllbcc03);
trajMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(trajMPC_ccrhsub03, trajMPC_sub03, trajMPC_lubbysub03, trajMPC_dzcc03, trajMPC_ubIdx03, trajMPC_dlubcc03);
trajMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(trajMPC_ccrhsl04, trajMPC_slb04, trajMPC_llbbyslb04, trajMPC_dzcc04, trajMPC_lbIdx04, trajMPC_dllbcc04);
trajMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(trajMPC_ccrhsub04, trajMPC_sub04, trajMPC_lubbysub04, trajMPC_dzcc04, trajMPC_ubIdx04, trajMPC_dlubcc04);
trajMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(trajMPC_ccrhsl05, trajMPC_slb05, trajMPC_llbbyslb05, trajMPC_dzcc05, trajMPC_lbIdx05, trajMPC_dllbcc05);
trajMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(trajMPC_ccrhsub05, trajMPC_sub05, trajMPC_lubbysub05, trajMPC_dzcc05, trajMPC_ubIdx05, trajMPC_dlubcc05);
trajMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(trajMPC_ccrhsl06, trajMPC_slb06, trajMPC_llbbyslb06, trajMPC_dzcc06, trajMPC_lbIdx06, trajMPC_dllbcc06);
trajMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(trajMPC_ccrhsub06, trajMPC_sub06, trajMPC_lubbysub06, trajMPC_dzcc06, trajMPC_ubIdx06, trajMPC_dlubcc06);
trajMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(trajMPC_ccrhsl07, trajMPC_slb07, trajMPC_llbbyslb07, trajMPC_dzcc07, trajMPC_lbIdx07, trajMPC_dllbcc07);
trajMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(trajMPC_ccrhsub07, trajMPC_sub07, trajMPC_lubbysub07, trajMPC_dzcc07, trajMPC_ubIdx07, trajMPC_dlubcc07);
trajMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(trajMPC_ccrhsl08, trajMPC_slb08, trajMPC_llbbyslb08, trajMPC_dzcc08, trajMPC_lbIdx08, trajMPC_dllbcc08);
trajMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(trajMPC_ccrhsub08, trajMPC_sub08, trajMPC_lubbysub08, trajMPC_dzcc08, trajMPC_ubIdx08, trajMPC_dlubcc08);
trajMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(trajMPC_ccrhsl09, trajMPC_slb09, trajMPC_llbbyslb09, trajMPC_dzcc09, trajMPC_lbIdx09, trajMPC_dllbcc09);
trajMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(trajMPC_ccrhsub09, trajMPC_sub09, trajMPC_lubbysub09, trajMPC_dzcc09, trajMPC_ubIdx09, trajMPC_dlubcc09);
trajMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(trajMPC_ccrhsl10, trajMPC_slb10, trajMPC_llbbyslb10, trajMPC_dzcc10, trajMPC_lbIdx10, trajMPC_dllbcc10);
trajMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(trajMPC_ccrhsub10, trajMPC_sub10, trajMPC_lubbysub10, trajMPC_dzcc10, trajMPC_ubIdx10, trajMPC_dlubcc10);
trajMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(trajMPC_ccrhsl11, trajMPC_slb11, trajMPC_llbbyslb11, trajMPC_dzcc11, trajMPC_lbIdx11, trajMPC_dllbcc11);
trajMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(trajMPC_ccrhsub11, trajMPC_sub11, trajMPC_lubbysub11, trajMPC_dzcc11, trajMPC_ubIdx11, trajMPC_dlubcc11);
trajMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(trajMPC_ccrhsl12, trajMPC_slb12, trajMPC_llbbyslb12, trajMPC_dzcc12, trajMPC_lbIdx12, trajMPC_dllbcc12);
trajMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(trajMPC_ccrhsub12, trajMPC_sub12, trajMPC_lubbysub12, trajMPC_dzcc12, trajMPC_ubIdx12, trajMPC_dlubcc12);
trajMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_11(trajMPC_ccrhsl13, trajMPC_slb13, trajMPC_llbbyslb13, trajMPC_dzcc13, trajMPC_lbIdx13, trajMPC_dllbcc13);
trajMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_5(trajMPC_ccrhsub13, trajMPC_sub13, trajMPC_lubbysub13, trajMPC_dzcc13, trajMPC_ubIdx13, trajMPC_dlubcc13);
trajMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_3(trajMPC_ccrhsl14, trajMPC_slb14, trajMPC_llbbyslb14, trajMPC_dzcc14, trajMPC_lbIdx14, trajMPC_dllbcc14);
trajMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_3(trajMPC_ccrhsub14, trajMPC_sub14, trajMPC_lubbysub14, trajMPC_dzcc14, trajMPC_ubIdx14, trajMPC_dlubcc14);
trajMPC_LA_VSUB7_230(trajMPC_l, trajMPC_ccrhs, trajMPC_s, trajMPC_dl_cc, trajMPC_ds_cc);
trajMPC_LA_VADD_157(trajMPC_dz_cc, trajMPC_dz_aff);
trajMPC_LA_VADD_45(trajMPC_dv_cc, trajMPC_dv_aff);
trajMPC_LA_VADD_230(trajMPC_dl_cc, trajMPC_dl_aff);
trajMPC_LA_VADD_230(trajMPC_ds_cc, trajMPC_ds_aff);
info->lsit_cc = trajMPC_LINESEARCH_BACKTRACKING_COMBINED(trajMPC_z, trajMPC_v, trajMPC_l, trajMPC_s, trajMPC_dz_cc, trajMPC_dv_cc, trajMPC_dl_cc, trajMPC_ds_cc, &info->step_cc, &info->mu);
if( info->lsit_cc == trajMPC_NOPROGRESS ){
exitcode = trajMPC_NOPROGRESS; break;
}
info->it++;
}
output->z1[0] = trajMPC_z00[0];
output->z1[1] = trajMPC_z00[1];
output->z1[2] = trajMPC_z00[2];
output->z1[3] = trajMPC_z00[3];
output->z1[4] = trajMPC_z00[4];
output->z2[0] = trajMPC_z01[0];
output->z2[1] = trajMPC_z01[1];
output->z2[2] = trajMPC_z01[2];
output->z2[3] = trajMPC_z01[3];
output->z2[4] = trajMPC_z01[4];
output->z3[0] = trajMPC_z02[0];
output->z3[1] = trajMPC_z02[1];
output->z3[2] = trajMPC_z02[2];
output->z3[3] = trajMPC_z02[3];
output->z3[4] = trajMPC_z02[4];
output->z4[0] = trajMPC_z03[0];
output->z4[1] = trajMPC_z03[1];
output->z4[2] = trajMPC_z03[2];
output->z4[3] = trajMPC_z03[3];
output->z4[4] = trajMPC_z03[4];
output->z5[0] = trajMPC_z04[0];
output->z5[1] = trajMPC_z04[1];
output->z5[2] = trajMPC_z04[2];
output->z5[3] = trajMPC_z04[3];
output->z5[4] = trajMPC_z04[4];
output->z6[0] = trajMPC_z05[0];
output->z6[1] = trajMPC_z05[1];
output->z6[2] = trajMPC_z05[2];
output->z6[3] = trajMPC_z05[3];
output->z6[4] = trajMPC_z05[4];
output->z7[0] = trajMPC_z06[0];
output->z7[1] = trajMPC_z06[1];
output->z7[2] = trajMPC_z06[2];
output->z7[3] = trajMPC_z06[3];
output->z7[4] = trajMPC_z06[4];
output->z8[0] = trajMPC_z07[0];
output->z8[1] = trajMPC_z07[1];
output->z8[2] = trajMPC_z07[2];
output->z8[3] = trajMPC_z07[3];
output->z8[4] = trajMPC_z07[4];
output->z9[0] = trajMPC_z08[0];
output->z9[1] = trajMPC_z08[1];
output->z9[2] = trajMPC_z08[2];
output->z9[3] = trajMPC_z08[3];
output->z9[4] = trajMPC_z08[4];
output->z10[0] = trajMPC_z09[0];
output->z10[1] = trajMPC_z09[1];
output->z10[2] = trajMPC_z09[2];
output->z10[3] = trajMPC_z09[3];
output->z10[4] = trajMPC_z09[4];
output->z11[0] = trajMPC_z10[0];
output->z11[1] = trajMPC_z10[1];
output->z11[2] = trajMPC_z10[2];
output->z11[3] = trajMPC_z10[3];
output->z11[4] = trajMPC_z10[4];
output->z12[0] = trajMPC_z11[0];
output->z12[1] = trajMPC_z11[1];
output->z12[2] = trajMPC_z11[2];
output->z12[3] = trajMPC_z11[3];
output->z12[4] = trajMPC_z11[4];
output->z13[0] = trajMPC_z12[0];
output->z13[1] = trajMPC_z12[1];
output->z13[2] = trajMPC_z12[2];
output->z13[3] = trajMPC_z12[3];
output->z13[4] = trajMPC_z12[4];
output->z14[0] = trajMPC_z13[0];
output->z14[1] = trajMPC_z13[1];
output->z14[2] = trajMPC_z13[2];
output->z14[3] = trajMPC_z13[3];
output->z14[4] = trajMPC_z13[4];
output->z15[0] = trajMPC_z14[0];
output->z15[1] = trajMPC_z14[1];
output->z15[2] = trajMPC_z14[2];

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
