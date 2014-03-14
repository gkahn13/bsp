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

#include "controlPenaltyMPC.h"

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
 * Initializes a vector of length 180 with a value.
 */
void controlPenaltyMPC_LA_INITIALIZEVECTOR_180(controlPenaltyMPC_FLOAT* vec, controlPenaltyMPC_FLOAT value)
{
	int i;
	for( i=0; i<180; i++ )
	{
		vec[i] = value;
	}
}


/*
 * Initializes a vector of length 90 with a value.
 */
void controlPenaltyMPC_LA_INITIALIZEVECTOR_90(controlPenaltyMPC_FLOAT* vec, controlPenaltyMPC_FLOAT value)
{
	int i;
	for( i=0; i<90; i++ )
	{
		vec[i] = value;
	}
}


/*
 * Initializes a vector of length 354 with a value.
 */
void controlPenaltyMPC_LA_INITIALIZEVECTOR_354(controlPenaltyMPC_FLOAT* vec, controlPenaltyMPC_FLOAT value)
{
	int i;
	for( i=0; i<354; i++ )
	{
		vec[i] = value;
	}
}


/* 
 * Calculates a dot product and adds it to a variable: z += x'*y; 
 * This function is for vectors of length 354.
 */
void controlPenaltyMPC_LA_DOTACC_354(controlPenaltyMPC_FLOAT *x, controlPenaltyMPC_FLOAT *y, controlPenaltyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<354; i++ ){
		*z += x[i]*y[i];
	}
}


/*
 * Calculates the gradient and the value for a quadratic function 0.5*z'*H*z + f'*z
 *
 * INPUTS:     H  - Symmetric Hessian, diag matrix of size [12 x 12]
 *             f  - column vector of size 12
 *             z  - column vector of size 12
 *
 * OUTPUTS: grad  - gradient at z (= H*z + f), column vector of size 12
 *          value <-- value + 0.5*z'*H*z + f'*z (value will be modified)
 */
void controlPenaltyMPC_LA_DIAG_QUADFCN_12(controlPenaltyMPC_FLOAT* H, controlPenaltyMPC_FLOAT* f, controlPenaltyMPC_FLOAT* z, controlPenaltyMPC_FLOAT* grad, controlPenaltyMPC_FLOAT* value)
{
	int i;
	controlPenaltyMPC_FLOAT hz;	
	for( i=0; i<12; i++){
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
void controlPenaltyMPC_LA_DIAGZERO_MVMSUB6_6(controlPenaltyMPC_FLOAT *B, controlPenaltyMPC_FLOAT *u, controlPenaltyMPC_FLOAT *b, controlPenaltyMPC_FLOAT *l, controlPenaltyMPC_FLOAT *r, controlPenaltyMPC_FLOAT *z, controlPenaltyMPC_FLOAT *y)
{
	int i;
	controlPenaltyMPC_FLOAT Bu[6];
	controlPenaltyMPC_FLOAT norm = *y;
	controlPenaltyMPC_FLOAT lr = 0;

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
void controlPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_6_12_12(controlPenaltyMPC_FLOAT *A, controlPenaltyMPC_FLOAT *x, controlPenaltyMPC_FLOAT *B, controlPenaltyMPC_FLOAT *u, controlPenaltyMPC_FLOAT *b, controlPenaltyMPC_FLOAT *l, controlPenaltyMPC_FLOAT *r, controlPenaltyMPC_FLOAT *z, controlPenaltyMPC_FLOAT *y)
{
	int i;
	int j;
	int k = 0;
	controlPenaltyMPC_FLOAT AxBu[6];
	controlPenaltyMPC_FLOAT norm = *y;
	controlPenaltyMPC_FLOAT lr = 0;

	/* do A*x + B*u first */
	for( i=0; i<6; i++ ){
		AxBu[i] = A[k++]*x[0] + B[i]*u[i];
	}	

	for( j=1; j<12; j++ ){		
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
 * where A is of size [6 x 12] and stored in column major format.
 * and B is of size [6 x 12] and stored in diagzero format
 * Note the transposes of A and B!
 */
void controlPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_6_12_6(controlPenaltyMPC_FLOAT *A, controlPenaltyMPC_FLOAT *x, controlPenaltyMPC_FLOAT *B, controlPenaltyMPC_FLOAT *y, controlPenaltyMPC_FLOAT *z)
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
	for( i=6 ;i<12; i++ ){
		z[i] = 0;
		for( j=0; j<6; j++ ){
			z[i] += A[k++]*x[j];
		}
	}
}


/*
 * Matrix vector multiplication y = M'*x where M is of size [6 x 12]
 * and stored in diagzero format. Note the transpose of M!
 */
void controlPenaltyMPC_LA_DIAGZERO_MTVM_6_12(controlPenaltyMPC_FLOAT *M, controlPenaltyMPC_FLOAT *x, controlPenaltyMPC_FLOAT *y)
{
	int i;
	for( i=0; i<12; i++ ){
		y[i] = M[i]*x[i];
	}
}


/*
 * Vector subtraction and addition.
 *	 Input: five vectors t, tidx, u, v, w and two scalars z and r
 *	 Output: y = t(tidx) - u + w
 *           z = z - v'*x;
 *           r = max([norm(y,inf), z]);
 * for vectors of length 12. Output z is of course scalar.
 */
void controlPenaltyMPC_LA_VSUBADD3_12(controlPenaltyMPC_FLOAT* t, controlPenaltyMPC_FLOAT* u, int* uidx, controlPenaltyMPC_FLOAT* v, controlPenaltyMPC_FLOAT* w, controlPenaltyMPC_FLOAT* y, controlPenaltyMPC_FLOAT* z, controlPenaltyMPC_FLOAT* r)
{
	int i;
	controlPenaltyMPC_FLOAT norm = *r;
	controlPenaltyMPC_FLOAT vx = 0;
	controlPenaltyMPC_FLOAT x;
	for( i=0; i<12; i++){
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
 * for vectors of length 12. Output z is of course scalar.
 */
void controlPenaltyMPC_LA_VSUBADD2_12(controlPenaltyMPC_FLOAT* t, int* tidx, controlPenaltyMPC_FLOAT* u, controlPenaltyMPC_FLOAT* v, controlPenaltyMPC_FLOAT* w, controlPenaltyMPC_FLOAT* y, controlPenaltyMPC_FLOAT* z, controlPenaltyMPC_FLOAT* r)
{
	int i;
	controlPenaltyMPC_FLOAT norm = *r;
	controlPenaltyMPC_FLOAT vx = 0;
	controlPenaltyMPC_FLOAT x;
	for( i=0; i<12; i++){
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
void controlPenaltyMPC_LA_VSUBADD2_6(controlPenaltyMPC_FLOAT* t, int* tidx, controlPenaltyMPC_FLOAT* u, controlPenaltyMPC_FLOAT* v, controlPenaltyMPC_FLOAT* w, controlPenaltyMPC_FLOAT* y, controlPenaltyMPC_FLOAT* z, controlPenaltyMPC_FLOAT* r)
{
	int i;
	controlPenaltyMPC_FLOAT norm = *r;
	controlPenaltyMPC_FLOAT vx = 0;
	controlPenaltyMPC_FLOAT x;
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
 * Special function for box constraints of length 12
 * Returns also L/S, a value that is often used elsewhere.
 */
void controlPenaltyMPC_LA_INEQ_B_GRAD_12_12_12(controlPenaltyMPC_FLOAT *lu, controlPenaltyMPC_FLOAT *su, controlPenaltyMPC_FLOAT *ru, controlPenaltyMPC_FLOAT *ll, controlPenaltyMPC_FLOAT *sl, controlPenaltyMPC_FLOAT *rl, int* lbIdx, int* ubIdx, controlPenaltyMPC_FLOAT *grad, controlPenaltyMPC_FLOAT *lubysu, controlPenaltyMPC_FLOAT *llbysl)
{
	int i;
	for( i=0; i<12; i++ ){
		grad[i] = 0;
	}
	for( i=0; i<12; i++ ){		
		llbysl[i] = ll[i] / sl[i];
		grad[lbIdx[i]] -= llbysl[i]*rl[i];
	}
	for( i=0; i<12; i++ ){
		lubysu[i] = lu[i] / su[i];
		grad[ubIdx[i]] += lubysu[i]*ru[i];
	}
}


/*
 * Computes inequality constraints gradient-
 * Special function for box constraints of length 12
 * Returns also L/S, a value that is often used elsewhere.
 */
void controlPenaltyMPC_LA_INEQ_B_GRAD_12_12_6(controlPenaltyMPC_FLOAT *lu, controlPenaltyMPC_FLOAT *su, controlPenaltyMPC_FLOAT *ru, controlPenaltyMPC_FLOAT *ll, controlPenaltyMPC_FLOAT *sl, controlPenaltyMPC_FLOAT *rl, int* lbIdx, int* ubIdx, controlPenaltyMPC_FLOAT *grad, controlPenaltyMPC_FLOAT *lubysu, controlPenaltyMPC_FLOAT *llbysl)
{
	int i;
	for( i=0; i<12; i++ ){
		grad[i] = 0;
	}
	for( i=0; i<12; i++ ){		
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
 * of length 180.
 */
void controlPenaltyMPC_LA_VVADD3_180(controlPenaltyMPC_FLOAT *u, controlPenaltyMPC_FLOAT *v, controlPenaltyMPC_FLOAT *w, controlPenaltyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<180; i++ ){
		z[i] = u[i] + v[i] + w[i];
	}
}


/*
 * Special function to compute the diagonal cholesky factorization of the 
 * positive definite augmented Hessian for block size 12.
 *
 * Inputs: - H = diagonal cost Hessian in diagonal storage format
 *         - llbysl = L / S of lower bounds
 *         - lubysu = L / S of upper bounds
 *
 * Output: Phi = sqrt(H + diag(llbysl) + diag(lubysu))
 * where Phi is stored in diagonal storage format
 */
void controlPenaltyMPC_LA_DIAG_CHOL_ONELOOP_LBUB_12_12_12(controlPenaltyMPC_FLOAT *H, controlPenaltyMPC_FLOAT *llbysl, int* lbIdx, controlPenaltyMPC_FLOAT *lubysu, int* ubIdx, controlPenaltyMPC_FLOAT *Phi)


{
	int i;
	
	/* compute cholesky */
	for( i=0; i<12; i++ ){
		Phi[i] = H[i] + llbysl[i] + lubysu[i];

#if controlPenaltyMPC_SET_PRINTLEVEL > 0 && defined PRINTNUMERICALWARNINGS
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
 * where A is to be computed and is of size [6 x 12],
 * B is given and of size [6 x 12], L is a diagonal
 * matrix of size 6 stored in diagonal matrix 
 * storage format. Note the transpose of L has no impact!
 *
 * Result: A in column major storage format.
 *
 */
void controlPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_6_12(controlPenaltyMPC_FLOAT *L, controlPenaltyMPC_FLOAT *B, controlPenaltyMPC_FLOAT *A)
{
    int i,j;
	 int k = 0;

	for( j=0; j<12; j++){
		for( i=0; i<6; i++){
			A[k] = B[k]/L[j];
			k++;
		}
	}

}


/**
 * Forward substitution for the matrix equation A*L' = B
 * where A is to be computed and is of size [6 x 12],
 * B is given and of size [6 x 12], L is a diagonal
 *  matrix of size 12 stored in diagonal 
 * storage format. Note the transpose of L!
 *
 * Result: A in diagonalzero storage format.
 *
 */
void controlPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_6_12(controlPenaltyMPC_FLOAT *L, controlPenaltyMPC_FLOAT *B, controlPenaltyMPC_FLOAT *A)
{
	int j;
    for( j=0; j<12; j++ ){   
		A[j] = B[j]/L[j];
     }
}


/**
 * Compute C = A*B' where 
 *
 *	size(A) = [6 x 12]
 *  size(B) = [6 x 12] in diagzero format
 * 
 * A and C matrices are stored in column major format.
 * 
 * 
 */
void controlPenaltyMPC_LA_DENSE_DIAGZERO_MMTM_6_12_6(controlPenaltyMPC_FLOAT *A, controlPenaltyMPC_FLOAT *B, controlPenaltyMPC_FLOAT *C)
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
 * The dimensions involved are 12.
 */
void controlPenaltyMPC_LA_DIAG_FORWARDSUB_12(controlPenaltyMPC_FLOAT *L, controlPenaltyMPC_FLOAT *b, controlPenaltyMPC_FLOAT *y)
{
    int i;

    for( i=0; i<12; i++ ){
		y[i] = b[i]/L[i];
    }
}


/*
 * Special function to compute the diagonal cholesky factorization of the 
 * positive definite augmented Hessian for block size 12.
 *
 * Inputs: - H = diagonal cost Hessian in diagonal storage format
 *         - llbysl = L / S of lower bounds
 *         - lubysu = L / S of upper bounds
 *
 * Output: Phi = sqrt(H + diag(llbysl) + diag(lubysu))
 * where Phi is stored in diagonal storage format
 */
void controlPenaltyMPC_LA_DIAG_CHOL_LBUB_12_12_6(controlPenaltyMPC_FLOAT *H, controlPenaltyMPC_FLOAT *llbysl, int* lbIdx, controlPenaltyMPC_FLOAT *lubysu, int* ubIdx, controlPenaltyMPC_FLOAT *Phi)


{
	int i;
	
	/* copy  H into PHI */
	for( i=0; i<12; i++ ){
		Phi[i] = H[i];
	}

	/* add llbysl onto Phi where necessary */
	for( i=0; i<12; i++ ){
		Phi[lbIdx[i]] += llbysl[i];
	}

	/* add lubysu onto Phi where necessary */
	for( i=0; i<6; i++){
		Phi[ubIdx[i]] +=  lubysu[i];
	}
	
	/* compute cholesky */
	for(i=0; i<12; i++)
	{
#if controlPenaltyMPC_SET_PRINTLEVEL > 0 && defined PRINTNUMERICALWARNINGS
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
 * Compute L = B*B', where L is lower triangular of size NXp1
 * and B is a diagzero matrix of size [6 x 12] in column
 * storage format.
 * 
 */
void controlPenaltyMPC_LA_DIAGZERO_MMT_6(controlPenaltyMPC_FLOAT *B, controlPenaltyMPC_FLOAT *L)
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
void controlPenaltyMPC_LA_DIAGZERO_MVMSUB7_6(controlPenaltyMPC_FLOAT *B, controlPenaltyMPC_FLOAT *u, controlPenaltyMPC_FLOAT *b, controlPenaltyMPC_FLOAT *r)
{
	int i;

	for( i=0; i<6; i++ ){
		r[i] = b[i] - B[i]*u[i];
	}	
	
}


/**
 * Compute L = A*A' + B*B', where L is lower triangular of size NXp1
 * and A is a dense matrix of size [6 x 12] in column
 * storage format, and B is of size [6 x 12] diagonalzero
 * storage format.
 * 
 * THIS ONE HAS THE WORST ACCES PATTERN POSSIBLE. 
 * POSSIBKE FIX: PUT A AND B INTO ROW MAJOR FORMAT FIRST.
 * 
 */
void controlPenaltyMPC_LA_DENSE_DIAGZERO_MMT2_6_12_12(controlPenaltyMPC_FLOAT *A, controlPenaltyMPC_FLOAT *B, controlPenaltyMPC_FLOAT *L)
{
    int i, j, k, ii, di;
    controlPenaltyMPC_FLOAT ltemp;
    
    ii = 0; di = 0;
    for( i=0; i<6; i++ ){        
        for( j=0; j<=i; j++ ){
            ltemp = 0; 
            for( k=0; k<12; k++ ){
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
void controlPenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_6_12_12(controlPenaltyMPC_FLOAT *A, controlPenaltyMPC_FLOAT *x, controlPenaltyMPC_FLOAT *B, controlPenaltyMPC_FLOAT *u, controlPenaltyMPC_FLOAT *b, controlPenaltyMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;

	for( i=0; i<6; i++ ){
		r[i] = b[i] - A[k++]*x[0] - B[i]*u[i];
	}	

	for( j=1; j<12; j++ ){		
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
void controlPenaltyMPC_LA_DENSE_CHOL_6(controlPenaltyMPC_FLOAT *A, controlPenaltyMPC_FLOAT *L)
{
    int i, j, k, di, dj;
	 int ii, jj;

    controlPenaltyMPC_FLOAT l;
    controlPenaltyMPC_FLOAT Mii;

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
        
#if controlPenaltyMPC_SET_PRINTLEVEL > 0 && defined PRINTNUMERICALWARNINGS
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
void controlPenaltyMPC_LA_DENSE_FORWARDSUB_6(controlPenaltyMPC_FLOAT *L, controlPenaltyMPC_FLOAT *b, controlPenaltyMPC_FLOAT *y)
{
    int i,j,ii,di;
    controlPenaltyMPC_FLOAT yel;
            
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
void controlPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_6_6(controlPenaltyMPC_FLOAT *L, controlPenaltyMPC_FLOAT *B, controlPenaltyMPC_FLOAT *A)
{
    int i,j,k,ii,di;
    controlPenaltyMPC_FLOAT a;
    
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
void controlPenaltyMPC_LA_DENSE_MMTSUB_6_6(controlPenaltyMPC_FLOAT *A, controlPenaltyMPC_FLOAT *L)
{
    int i, j, k, ii, di;
    controlPenaltyMPC_FLOAT ltemp;
    
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
void controlPenaltyMPC_LA_DENSE_MVMSUB1_6_6(controlPenaltyMPC_FLOAT *A, controlPenaltyMPC_FLOAT *x, controlPenaltyMPC_FLOAT *b, controlPenaltyMPC_FLOAT *r)
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
void controlPenaltyMPC_LA_DENSE_BACKWARDSUB_6(controlPenaltyMPC_FLOAT *L, controlPenaltyMPC_FLOAT *y, controlPenaltyMPC_FLOAT *x)
{
    int i, ii, di, j, jj, dj;
    controlPenaltyMPC_FLOAT xel;    
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
void controlPenaltyMPC_LA_DENSE_MTVMSUB_6_6(controlPenaltyMPC_FLOAT *A, controlPenaltyMPC_FLOAT *x, controlPenaltyMPC_FLOAT *b, controlPenaltyMPC_FLOAT *r)
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
 * Vector subtraction z = -x - y for vectors of length 180.
 */
void controlPenaltyMPC_LA_VSUB2_180(controlPenaltyMPC_FLOAT *x, controlPenaltyMPC_FLOAT *y, controlPenaltyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<180; i++){
		z[i] = -x[i] - y[i];
	}
}


/**
 * Forward-Backward-Substitution to solve L*L^T*x = b where L is a
 * diagonal matrix of size 12 in vector
 * storage format.
 */
void controlPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_12(controlPenaltyMPC_FLOAT *L, controlPenaltyMPC_FLOAT *b, controlPenaltyMPC_FLOAT *x)
{
    int i;
            
    /* solve Ly = b by forward and backward substitution */
    for( i=0; i<12; i++ ){
		x[i] = b[i]/(L[i]*L[i]);
    }
    
}


/*
 * Vector subtraction z = x(xidx) - y where y, z and xidx are of length 12,
 * and x has length 12 and is indexed through yidx.
 */
void controlPenaltyMPC_LA_VSUB_INDEXED_12(controlPenaltyMPC_FLOAT *x, int* xidx, controlPenaltyMPC_FLOAT *y, controlPenaltyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<12; i++){
		z[i] = x[xidx[i]] - y[i];
	}
}


/*
 * Vector subtraction x = -u.*v - w for vectors of length 12.
 */
void controlPenaltyMPC_LA_VSUB3_12(controlPenaltyMPC_FLOAT *u, controlPenaltyMPC_FLOAT *v, controlPenaltyMPC_FLOAT *w, controlPenaltyMPC_FLOAT *x)
{
	int i;
	for( i=0; i<12; i++){
		x[i] = -u[i]*v[i] - w[i];
	}
}


/*
 * Vector subtraction z = -x - y(yidx) where y is of length 12
 * and z, x and yidx are of length 12.
 */
void controlPenaltyMPC_LA_VSUB2_INDEXED_12(controlPenaltyMPC_FLOAT *x, controlPenaltyMPC_FLOAT *y, int* yidx, controlPenaltyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<12; i++){
		z[i] = -x[i] - y[yidx[i]];
	}
}


/*
 * Vector subtraction z = -x - y(yidx) where y is of length 12
 * and z, x and yidx are of length 6.
 */
void controlPenaltyMPC_LA_VSUB2_INDEXED_6(controlPenaltyMPC_FLOAT *x, controlPenaltyMPC_FLOAT *y, int* yidx, controlPenaltyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<6; i++){
		z[i] = -x[i] - y[yidx[i]];
	}
}


/*
 * Vector subtraction x = -u.*v - w for vectors of length 6.
 */
void controlPenaltyMPC_LA_VSUB3_6(controlPenaltyMPC_FLOAT *u, controlPenaltyMPC_FLOAT *v, controlPenaltyMPC_FLOAT *w, controlPenaltyMPC_FLOAT *x)
{
	int i;
	for( i=0; i<6; i++){
		x[i] = -u[i]*v[i] - w[i];
	}
}


/**
 * Backtracking line search.
 * 
 * First determine the maximum line length by a feasibility line
 * search, i.e. a ~= argmax{ a \in [0...1] s.t. l+a*dl >= 0 and s+a*ds >= 0}.
 *
 * The function returns either the number of iterations or exits the error code
 * controlPenaltyMPC_NOPROGRESS (should be negative).
 */
int controlPenaltyMPC_LINESEARCH_BACKTRACKING_AFFINE(controlPenaltyMPC_FLOAT *l, controlPenaltyMPC_FLOAT *s, controlPenaltyMPC_FLOAT *dl, controlPenaltyMPC_FLOAT *ds, controlPenaltyMPC_FLOAT *a, controlPenaltyMPC_FLOAT *mu_aff)
{
    int i;
	int lsIt=1;    
    controlPenaltyMPC_FLOAT dltemp;
    controlPenaltyMPC_FLOAT dstemp;
    controlPenaltyMPC_FLOAT mya = 1.0;
    controlPenaltyMPC_FLOAT mymu;
        
    while( 1 ){                        

        /* 
         * Compute both snew and wnew together.
         * We compute also mu_affine along the way here, as the
         * values might be in registers, so it should be cheaper.
         */
        mymu = 0;
        for( i=0; i<354; i++ ){
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
        if( i == 354 ){
            break;
        } else {
            mya *= controlPenaltyMPC_SET_LS_SCALE_AFF;
            if( mya < controlPenaltyMPC_SET_LS_MINSTEP ){
                return controlPenaltyMPC_NOPROGRESS;
            }
        }
    }
    
    /* return new values and iteration counter */
    *a = mya;
    *mu_aff = mymu / (controlPenaltyMPC_FLOAT)354;
    return lsIt;
}


/*
 * Vector subtraction x = u.*v - a where a is a scalar
*  and x,u,v are vectors of length 354.
 */
void controlPenaltyMPC_LA_VSUB5_354(controlPenaltyMPC_FLOAT *u, controlPenaltyMPC_FLOAT *v, controlPenaltyMPC_FLOAT a, controlPenaltyMPC_FLOAT *x)
{
	int i;
	for( i=0; i<354; i++){
		x[i] = u[i]*v[i] - a;
	}
}


/*
 * Computes x=0; x(uidx) += u/su; x(vidx) -= v/sv where x is of length 12,
 * u, su, uidx are of length 12 and v, sv, vidx are of length 12.
 */
void controlPenaltyMPC_LA_VSUB6_INDEXED_12_12_12(controlPenaltyMPC_FLOAT *u, controlPenaltyMPC_FLOAT *su, int* uidx, controlPenaltyMPC_FLOAT *v, controlPenaltyMPC_FLOAT *sv, int* vidx, controlPenaltyMPC_FLOAT *x)
{
	int i;
	for( i=0; i<12; i++ ){
		x[i] = 0;
	}
	for( i=0; i<12; i++){
		x[uidx[i]] += u[i]/su[i];
	}
	for( i=0; i<12; i++){
		x[vidx[i]] -= v[i]/sv[i];
	}
}


/* 
 * Computes r =  B*u
 * where B is stored in diagzero format
 */
void controlPenaltyMPC_LA_DIAGZERO_MVM_6(controlPenaltyMPC_FLOAT *B, controlPenaltyMPC_FLOAT *u, controlPenaltyMPC_FLOAT *r)
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
void controlPenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_6_12_12(controlPenaltyMPC_FLOAT *A, controlPenaltyMPC_FLOAT *x, controlPenaltyMPC_FLOAT *B, controlPenaltyMPC_FLOAT *u, controlPenaltyMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;

	for( i=0; i<6; i++ ){
		r[i] = A[k++]*x[0] + B[i]*u[i];
	}	

	for( j=1; j<12; j++ ){		
		for( i=0; i<6; i++ ){
			r[i] += A[k++]*x[j];
		}
	}
	
}


/*
 * Computes x=0; x(uidx) += u/su; x(vidx) -= v/sv where x is of length 12,
 * u, su, uidx are of length 6 and v, sv, vidx are of length 12.
 */
void controlPenaltyMPC_LA_VSUB6_INDEXED_12_6_12(controlPenaltyMPC_FLOAT *u, controlPenaltyMPC_FLOAT *su, int* uidx, controlPenaltyMPC_FLOAT *v, controlPenaltyMPC_FLOAT *sv, int* vidx, controlPenaltyMPC_FLOAT *x)
{
	int i;
	for( i=0; i<12; i++ ){
		x[i] = 0;
	}
	for( i=0; i<6; i++){
		x[uidx[i]] += u[i]/su[i];
	}
	for( i=0; i<12; i++){
		x[vidx[i]] -= v[i]/sv[i];
	}
}


/*
 * Vector subtraction z = x - y for vectors of length 180.
 */
void controlPenaltyMPC_LA_VSUB_180(controlPenaltyMPC_FLOAT *x, controlPenaltyMPC_FLOAT *y, controlPenaltyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<180; i++){
		z[i] = x[i] - y[i];
	}
}


/** 
 * Computes z = -r./s - u.*y(y)
 * where all vectors except of y are of length 12 (length of y >= 12).
 */
void controlPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_12(controlPenaltyMPC_FLOAT *r, controlPenaltyMPC_FLOAT *s, controlPenaltyMPC_FLOAT *u, controlPenaltyMPC_FLOAT *y, int* yidx, controlPenaltyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<12; i++ ){
		z[i] = -r[i]/s[i] - u[i]*y[yidx[i]];
	}
}


/** 
 * Computes z = -r./s + u.*y(y)
 * where all vectors except of y are of length 12 (length of y >= 12).
 */
void controlPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_12(controlPenaltyMPC_FLOAT *r, controlPenaltyMPC_FLOAT *s, controlPenaltyMPC_FLOAT *u, controlPenaltyMPC_FLOAT *y, int* yidx, controlPenaltyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<12; i++ ){
		z[i] = -r[i]/s[i] + u[i]*y[yidx[i]];
	}
}


/** 
 * Computes z = -r./s + u.*y(y)
 * where all vectors except of y are of length 6 (length of y >= 6).
 */
void controlPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_6(controlPenaltyMPC_FLOAT *r, controlPenaltyMPC_FLOAT *s, controlPenaltyMPC_FLOAT *u, controlPenaltyMPC_FLOAT *y, int* yidx, controlPenaltyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<6; i++ ){
		z[i] = -r[i]/s[i] + u[i]*y[yidx[i]];
	}
}


/*
 * Computes ds = -l.\(r + s.*dl) for vectors of length 354.
 */
void controlPenaltyMPC_LA_VSUB7_354(controlPenaltyMPC_FLOAT *l, controlPenaltyMPC_FLOAT *r, controlPenaltyMPC_FLOAT *s, controlPenaltyMPC_FLOAT *dl, controlPenaltyMPC_FLOAT *ds)
{
	int i;
	for( i=0; i<354; i++){
		ds[i] = -(r[i] + s[i]*dl[i])/l[i];
	}
}


/*
 * Vector addition x = x + y for vectors of length 180.
 */
void controlPenaltyMPC_LA_VADD_180(controlPenaltyMPC_FLOAT *x, controlPenaltyMPC_FLOAT *y)
{
	int i;
	for( i=0; i<180; i++){
		x[i] += y[i];
	}
}


/*
 * Vector addition x = x + y for vectors of length 90.
 */
void controlPenaltyMPC_LA_VADD_90(controlPenaltyMPC_FLOAT *x, controlPenaltyMPC_FLOAT *y)
{
	int i;
	for( i=0; i<90; i++){
		x[i] += y[i];
	}
}


/*
 * Vector addition x = x + y for vectors of length 354.
 */
void controlPenaltyMPC_LA_VADD_354(controlPenaltyMPC_FLOAT *x, controlPenaltyMPC_FLOAT *y)
{
	int i;
	for( i=0; i<354; i++){
		x[i] += y[i];
	}
}


/**
 * Backtracking line search for combined predictor/corrector step.
 * Update on variables with safety factor gamma (to keep us away from
 * boundary).
 */
int controlPenaltyMPC_LINESEARCH_BACKTRACKING_COMBINED(controlPenaltyMPC_FLOAT *z, controlPenaltyMPC_FLOAT *v, controlPenaltyMPC_FLOAT *l, controlPenaltyMPC_FLOAT *s, controlPenaltyMPC_FLOAT *dz, controlPenaltyMPC_FLOAT *dv, controlPenaltyMPC_FLOAT *dl, controlPenaltyMPC_FLOAT *ds, controlPenaltyMPC_FLOAT *a, controlPenaltyMPC_FLOAT *mu)
{
    int i, lsIt=1;       
    controlPenaltyMPC_FLOAT dltemp;
    controlPenaltyMPC_FLOAT dstemp;    
    controlPenaltyMPC_FLOAT a_gamma;
            
    *a = 1.0;
    while( 1 ){                        

        /* check whether search criterion is fulfilled */
        for( i=0; i<354; i++ ){
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
        if( i == 354 ){
            break;
        } else {
            *a *= controlPenaltyMPC_SET_LS_SCALE;
            if( *a < controlPenaltyMPC_SET_LS_MINSTEP ){
                return controlPenaltyMPC_NOPROGRESS;
            }
        }
    }
    
    /* update variables with safety margin */
    a_gamma = (*a)*controlPenaltyMPC_SET_LS_MAXSTEP;
    
    /* primal variables */
    for( i=0; i<180; i++ ){
        z[i] += a_gamma*dz[i];
    }
    
    /* equality constraint multipliers */
    for( i=0; i<90; i++ ){
        v[i] += a_gamma*dv[i];
    }
    
    /* inequality constraint multipliers & slacks, also update mu */
    *mu = 0;
    for( i=0; i<354; i++ ){
        dltemp = l[i] + a_gamma*dl[i]; l[i] = dltemp;
        dstemp = s[i] + a_gamma*ds[i]; s[i] = dstemp;
        *mu += dltemp*dstemp;
    }
    
    *a = a_gamma;
    *mu /= (controlPenaltyMPC_FLOAT)354;
    return lsIt;
}




/* VARIABLE DEFINITIONS ------------------------------------------------ */
controlPenaltyMPC_FLOAT controlPenaltyMPC_z[180];
controlPenaltyMPC_FLOAT controlPenaltyMPC_v[90];
controlPenaltyMPC_FLOAT controlPenaltyMPC_dz_aff[180];
controlPenaltyMPC_FLOAT controlPenaltyMPC_dv_aff[90];
controlPenaltyMPC_FLOAT controlPenaltyMPC_grad_cost[180];
controlPenaltyMPC_FLOAT controlPenaltyMPC_grad_eq[180];
controlPenaltyMPC_FLOAT controlPenaltyMPC_rd[180];
controlPenaltyMPC_FLOAT controlPenaltyMPC_l[354];
controlPenaltyMPC_FLOAT controlPenaltyMPC_s[354];
controlPenaltyMPC_FLOAT controlPenaltyMPC_lbys[354];
controlPenaltyMPC_FLOAT controlPenaltyMPC_dl_aff[354];
controlPenaltyMPC_FLOAT controlPenaltyMPC_ds_aff[354];
controlPenaltyMPC_FLOAT controlPenaltyMPC_dz_cc[180];
controlPenaltyMPC_FLOAT controlPenaltyMPC_dv_cc[90];
controlPenaltyMPC_FLOAT controlPenaltyMPC_dl_cc[354];
controlPenaltyMPC_FLOAT controlPenaltyMPC_ds_cc[354];
controlPenaltyMPC_FLOAT controlPenaltyMPC_ccrhs[354];
controlPenaltyMPC_FLOAT controlPenaltyMPC_grad_ineq[180];
controlPenaltyMPC_FLOAT* controlPenaltyMPC_z00 = controlPenaltyMPC_z + 0;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dzaff00 = controlPenaltyMPC_dz_aff + 0;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dzcc00 = controlPenaltyMPC_dz_cc + 0;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_rd00 = controlPenaltyMPC_rd + 0;
controlPenaltyMPC_FLOAT controlPenaltyMPC_Lbyrd00[12];
controlPenaltyMPC_FLOAT* controlPenaltyMPC_grad_cost00 = controlPenaltyMPC_grad_cost + 0;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_grad_eq00 = controlPenaltyMPC_grad_eq + 0;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_grad_ineq00 = controlPenaltyMPC_grad_ineq + 0;
controlPenaltyMPC_FLOAT controlPenaltyMPC_ctv00[12];
controlPenaltyMPC_FLOAT controlPenaltyMPC_C00[72] = {1.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 
0.0000000000000000E+000, 1.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 
0.0000000000000000E+000, 0.0000000000000000E+000, 1.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 
0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 1.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 
0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 1.0000000000000000E+000, 0.0000000000000000E+000, 
0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 1.0000000000000000E+000, 
1.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 
0.0000000000000000E+000, 1.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 
0.0000000000000000E+000, 0.0000000000000000E+000, 1.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 
0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 1.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 
0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 1.0000000000000000E+000, 0.0000000000000000E+000, 
0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 1.0000000000000000E+000};
controlPenaltyMPC_FLOAT* controlPenaltyMPC_v00 = controlPenaltyMPC_v + 0;
controlPenaltyMPC_FLOAT controlPenaltyMPC_re00[6];
controlPenaltyMPC_FLOAT controlPenaltyMPC_beta00[6];
controlPenaltyMPC_FLOAT controlPenaltyMPC_betacc00[6];
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dvaff00 = controlPenaltyMPC_dv_aff + 0;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dvcc00 = controlPenaltyMPC_dv_cc + 0;
controlPenaltyMPC_FLOAT controlPenaltyMPC_V00[72];
controlPenaltyMPC_FLOAT controlPenaltyMPC_Yd00[21];
controlPenaltyMPC_FLOAT controlPenaltyMPC_Ld00[21];
controlPenaltyMPC_FLOAT controlPenaltyMPC_yy00[6];
controlPenaltyMPC_FLOAT controlPenaltyMPC_bmy00[6];
int controlPenaltyMPC_lbIdx00[12] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
controlPenaltyMPC_FLOAT* controlPenaltyMPC_llb00 = controlPenaltyMPC_l + 0;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_slb00 = controlPenaltyMPC_s + 0;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_llbbyslb00 = controlPenaltyMPC_lbys + 0;
controlPenaltyMPC_FLOAT controlPenaltyMPC_rilb00[12];
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dllbaff00 = controlPenaltyMPC_dl_aff + 0;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dslbaff00 = controlPenaltyMPC_ds_aff + 0;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dllbcc00 = controlPenaltyMPC_dl_cc + 0;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dslbcc00 = controlPenaltyMPC_ds_cc + 0;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_ccrhsl00 = controlPenaltyMPC_ccrhs + 0;
int controlPenaltyMPC_ubIdx00[12] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
controlPenaltyMPC_FLOAT* controlPenaltyMPC_lub00 = controlPenaltyMPC_l + 12;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_sub00 = controlPenaltyMPC_s + 12;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_lubbysub00 = controlPenaltyMPC_lbys + 12;
controlPenaltyMPC_FLOAT controlPenaltyMPC_riub00[12];
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dlubaff00 = controlPenaltyMPC_dl_aff + 12;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dsubaff00 = controlPenaltyMPC_ds_aff + 12;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dlubcc00 = controlPenaltyMPC_dl_cc + 12;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dsubcc00 = controlPenaltyMPC_ds_cc + 12;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_ccrhsub00 = controlPenaltyMPC_ccrhs + 12;
controlPenaltyMPC_FLOAT controlPenaltyMPC_Phi00[12];
controlPenaltyMPC_FLOAT controlPenaltyMPC_D00[12] = {1.0000000000000000E+000, 
1.0000000000000000E+000, 
1.0000000000000000E+000, 
1.0000000000000000E+000, 
1.0000000000000000E+000, 
1.0000000000000000E+000};
controlPenaltyMPC_FLOAT controlPenaltyMPC_W00[12];
controlPenaltyMPC_FLOAT* controlPenaltyMPC_z01 = controlPenaltyMPC_z + 12;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dzaff01 = controlPenaltyMPC_dz_aff + 12;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dzcc01 = controlPenaltyMPC_dz_cc + 12;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_rd01 = controlPenaltyMPC_rd + 12;
controlPenaltyMPC_FLOAT controlPenaltyMPC_Lbyrd01[12];
controlPenaltyMPC_FLOAT* controlPenaltyMPC_grad_cost01 = controlPenaltyMPC_grad_cost + 12;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_grad_eq01 = controlPenaltyMPC_grad_eq + 12;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_grad_ineq01 = controlPenaltyMPC_grad_ineq + 12;
controlPenaltyMPC_FLOAT controlPenaltyMPC_ctv01[12];
controlPenaltyMPC_FLOAT* controlPenaltyMPC_v01 = controlPenaltyMPC_v + 6;
controlPenaltyMPC_FLOAT controlPenaltyMPC_re01[6];
controlPenaltyMPC_FLOAT controlPenaltyMPC_beta01[6];
controlPenaltyMPC_FLOAT controlPenaltyMPC_betacc01[6];
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dvaff01 = controlPenaltyMPC_dv_aff + 6;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dvcc01 = controlPenaltyMPC_dv_cc + 6;
controlPenaltyMPC_FLOAT controlPenaltyMPC_V01[72];
controlPenaltyMPC_FLOAT controlPenaltyMPC_Yd01[21];
controlPenaltyMPC_FLOAT controlPenaltyMPC_Ld01[21];
controlPenaltyMPC_FLOAT controlPenaltyMPC_yy01[6];
controlPenaltyMPC_FLOAT controlPenaltyMPC_bmy01[6];
controlPenaltyMPC_FLOAT controlPenaltyMPC_c01[6] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
int controlPenaltyMPC_lbIdx01[12] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
controlPenaltyMPC_FLOAT* controlPenaltyMPC_llb01 = controlPenaltyMPC_l + 24;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_slb01 = controlPenaltyMPC_s + 24;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_llbbyslb01 = controlPenaltyMPC_lbys + 24;
controlPenaltyMPC_FLOAT controlPenaltyMPC_rilb01[12];
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dllbaff01 = controlPenaltyMPC_dl_aff + 24;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dslbaff01 = controlPenaltyMPC_ds_aff + 24;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dllbcc01 = controlPenaltyMPC_dl_cc + 24;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dslbcc01 = controlPenaltyMPC_ds_cc + 24;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_ccrhsl01 = controlPenaltyMPC_ccrhs + 24;
int controlPenaltyMPC_ubIdx01[12] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
controlPenaltyMPC_FLOAT* controlPenaltyMPC_lub01 = controlPenaltyMPC_l + 36;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_sub01 = controlPenaltyMPC_s + 36;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_lubbysub01 = controlPenaltyMPC_lbys + 36;
controlPenaltyMPC_FLOAT controlPenaltyMPC_riub01[12];
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dlubaff01 = controlPenaltyMPC_dl_aff + 36;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dsubaff01 = controlPenaltyMPC_ds_aff + 36;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dlubcc01 = controlPenaltyMPC_dl_cc + 36;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dsubcc01 = controlPenaltyMPC_ds_cc + 36;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_ccrhsub01 = controlPenaltyMPC_ccrhs + 36;
controlPenaltyMPC_FLOAT controlPenaltyMPC_Phi01[12];
controlPenaltyMPC_FLOAT controlPenaltyMPC_D01[12] = {-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000, 
-1.0000000000000000E+000};
controlPenaltyMPC_FLOAT controlPenaltyMPC_W01[12];
controlPenaltyMPC_FLOAT controlPenaltyMPC_Ysd01[36];
controlPenaltyMPC_FLOAT controlPenaltyMPC_Lsd01[36];
controlPenaltyMPC_FLOAT* controlPenaltyMPC_z02 = controlPenaltyMPC_z + 24;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dzaff02 = controlPenaltyMPC_dz_aff + 24;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dzcc02 = controlPenaltyMPC_dz_cc + 24;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_rd02 = controlPenaltyMPC_rd + 24;
controlPenaltyMPC_FLOAT controlPenaltyMPC_Lbyrd02[12];
controlPenaltyMPC_FLOAT* controlPenaltyMPC_grad_cost02 = controlPenaltyMPC_grad_cost + 24;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_grad_eq02 = controlPenaltyMPC_grad_eq + 24;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_grad_ineq02 = controlPenaltyMPC_grad_ineq + 24;
controlPenaltyMPC_FLOAT controlPenaltyMPC_ctv02[12];
controlPenaltyMPC_FLOAT* controlPenaltyMPC_v02 = controlPenaltyMPC_v + 12;
controlPenaltyMPC_FLOAT controlPenaltyMPC_re02[6];
controlPenaltyMPC_FLOAT controlPenaltyMPC_beta02[6];
controlPenaltyMPC_FLOAT controlPenaltyMPC_betacc02[6];
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dvaff02 = controlPenaltyMPC_dv_aff + 12;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dvcc02 = controlPenaltyMPC_dv_cc + 12;
controlPenaltyMPC_FLOAT controlPenaltyMPC_V02[72];
controlPenaltyMPC_FLOAT controlPenaltyMPC_Yd02[21];
controlPenaltyMPC_FLOAT controlPenaltyMPC_Ld02[21];
controlPenaltyMPC_FLOAT controlPenaltyMPC_yy02[6];
controlPenaltyMPC_FLOAT controlPenaltyMPC_bmy02[6];
controlPenaltyMPC_FLOAT controlPenaltyMPC_c02[6] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
int controlPenaltyMPC_lbIdx02[12] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
controlPenaltyMPC_FLOAT* controlPenaltyMPC_llb02 = controlPenaltyMPC_l + 48;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_slb02 = controlPenaltyMPC_s + 48;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_llbbyslb02 = controlPenaltyMPC_lbys + 48;
controlPenaltyMPC_FLOAT controlPenaltyMPC_rilb02[12];
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dllbaff02 = controlPenaltyMPC_dl_aff + 48;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dslbaff02 = controlPenaltyMPC_ds_aff + 48;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dllbcc02 = controlPenaltyMPC_dl_cc + 48;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dslbcc02 = controlPenaltyMPC_ds_cc + 48;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_ccrhsl02 = controlPenaltyMPC_ccrhs + 48;
int controlPenaltyMPC_ubIdx02[12] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
controlPenaltyMPC_FLOAT* controlPenaltyMPC_lub02 = controlPenaltyMPC_l + 60;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_sub02 = controlPenaltyMPC_s + 60;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_lubbysub02 = controlPenaltyMPC_lbys + 60;
controlPenaltyMPC_FLOAT controlPenaltyMPC_riub02[12];
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dlubaff02 = controlPenaltyMPC_dl_aff + 60;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dsubaff02 = controlPenaltyMPC_ds_aff + 60;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dlubcc02 = controlPenaltyMPC_dl_cc + 60;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dsubcc02 = controlPenaltyMPC_ds_cc + 60;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_ccrhsub02 = controlPenaltyMPC_ccrhs + 60;
controlPenaltyMPC_FLOAT controlPenaltyMPC_Phi02[12];
controlPenaltyMPC_FLOAT controlPenaltyMPC_W02[12];
controlPenaltyMPC_FLOAT controlPenaltyMPC_Ysd02[36];
controlPenaltyMPC_FLOAT controlPenaltyMPC_Lsd02[36];
controlPenaltyMPC_FLOAT* controlPenaltyMPC_z03 = controlPenaltyMPC_z + 36;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dzaff03 = controlPenaltyMPC_dz_aff + 36;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dzcc03 = controlPenaltyMPC_dz_cc + 36;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_rd03 = controlPenaltyMPC_rd + 36;
controlPenaltyMPC_FLOAT controlPenaltyMPC_Lbyrd03[12];
controlPenaltyMPC_FLOAT* controlPenaltyMPC_grad_cost03 = controlPenaltyMPC_grad_cost + 36;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_grad_eq03 = controlPenaltyMPC_grad_eq + 36;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_grad_ineq03 = controlPenaltyMPC_grad_ineq + 36;
controlPenaltyMPC_FLOAT controlPenaltyMPC_ctv03[12];
controlPenaltyMPC_FLOAT* controlPenaltyMPC_v03 = controlPenaltyMPC_v + 18;
controlPenaltyMPC_FLOAT controlPenaltyMPC_re03[6];
controlPenaltyMPC_FLOAT controlPenaltyMPC_beta03[6];
controlPenaltyMPC_FLOAT controlPenaltyMPC_betacc03[6];
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dvaff03 = controlPenaltyMPC_dv_aff + 18;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dvcc03 = controlPenaltyMPC_dv_cc + 18;
controlPenaltyMPC_FLOAT controlPenaltyMPC_V03[72];
controlPenaltyMPC_FLOAT controlPenaltyMPC_Yd03[21];
controlPenaltyMPC_FLOAT controlPenaltyMPC_Ld03[21];
controlPenaltyMPC_FLOAT controlPenaltyMPC_yy03[6];
controlPenaltyMPC_FLOAT controlPenaltyMPC_bmy03[6];
controlPenaltyMPC_FLOAT controlPenaltyMPC_c03[6] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
int controlPenaltyMPC_lbIdx03[12] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
controlPenaltyMPC_FLOAT* controlPenaltyMPC_llb03 = controlPenaltyMPC_l + 72;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_slb03 = controlPenaltyMPC_s + 72;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_llbbyslb03 = controlPenaltyMPC_lbys + 72;
controlPenaltyMPC_FLOAT controlPenaltyMPC_rilb03[12];
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dllbaff03 = controlPenaltyMPC_dl_aff + 72;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dslbaff03 = controlPenaltyMPC_ds_aff + 72;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dllbcc03 = controlPenaltyMPC_dl_cc + 72;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dslbcc03 = controlPenaltyMPC_ds_cc + 72;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_ccrhsl03 = controlPenaltyMPC_ccrhs + 72;
int controlPenaltyMPC_ubIdx03[12] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
controlPenaltyMPC_FLOAT* controlPenaltyMPC_lub03 = controlPenaltyMPC_l + 84;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_sub03 = controlPenaltyMPC_s + 84;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_lubbysub03 = controlPenaltyMPC_lbys + 84;
controlPenaltyMPC_FLOAT controlPenaltyMPC_riub03[12];
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dlubaff03 = controlPenaltyMPC_dl_aff + 84;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dsubaff03 = controlPenaltyMPC_ds_aff + 84;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dlubcc03 = controlPenaltyMPC_dl_cc + 84;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dsubcc03 = controlPenaltyMPC_ds_cc + 84;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_ccrhsub03 = controlPenaltyMPC_ccrhs + 84;
controlPenaltyMPC_FLOAT controlPenaltyMPC_Phi03[12];
controlPenaltyMPC_FLOAT controlPenaltyMPC_W03[12];
controlPenaltyMPC_FLOAT controlPenaltyMPC_Ysd03[36];
controlPenaltyMPC_FLOAT controlPenaltyMPC_Lsd03[36];
controlPenaltyMPC_FLOAT* controlPenaltyMPC_z04 = controlPenaltyMPC_z + 48;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dzaff04 = controlPenaltyMPC_dz_aff + 48;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dzcc04 = controlPenaltyMPC_dz_cc + 48;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_rd04 = controlPenaltyMPC_rd + 48;
controlPenaltyMPC_FLOAT controlPenaltyMPC_Lbyrd04[12];
controlPenaltyMPC_FLOAT* controlPenaltyMPC_grad_cost04 = controlPenaltyMPC_grad_cost + 48;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_grad_eq04 = controlPenaltyMPC_grad_eq + 48;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_grad_ineq04 = controlPenaltyMPC_grad_ineq + 48;
controlPenaltyMPC_FLOAT controlPenaltyMPC_ctv04[12];
controlPenaltyMPC_FLOAT* controlPenaltyMPC_v04 = controlPenaltyMPC_v + 24;
controlPenaltyMPC_FLOAT controlPenaltyMPC_re04[6];
controlPenaltyMPC_FLOAT controlPenaltyMPC_beta04[6];
controlPenaltyMPC_FLOAT controlPenaltyMPC_betacc04[6];
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dvaff04 = controlPenaltyMPC_dv_aff + 24;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dvcc04 = controlPenaltyMPC_dv_cc + 24;
controlPenaltyMPC_FLOAT controlPenaltyMPC_V04[72];
controlPenaltyMPC_FLOAT controlPenaltyMPC_Yd04[21];
controlPenaltyMPC_FLOAT controlPenaltyMPC_Ld04[21];
controlPenaltyMPC_FLOAT controlPenaltyMPC_yy04[6];
controlPenaltyMPC_FLOAT controlPenaltyMPC_bmy04[6];
controlPenaltyMPC_FLOAT controlPenaltyMPC_c04[6] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
int controlPenaltyMPC_lbIdx04[12] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
controlPenaltyMPC_FLOAT* controlPenaltyMPC_llb04 = controlPenaltyMPC_l + 96;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_slb04 = controlPenaltyMPC_s + 96;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_llbbyslb04 = controlPenaltyMPC_lbys + 96;
controlPenaltyMPC_FLOAT controlPenaltyMPC_rilb04[12];
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dllbaff04 = controlPenaltyMPC_dl_aff + 96;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dslbaff04 = controlPenaltyMPC_ds_aff + 96;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dllbcc04 = controlPenaltyMPC_dl_cc + 96;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dslbcc04 = controlPenaltyMPC_ds_cc + 96;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_ccrhsl04 = controlPenaltyMPC_ccrhs + 96;
int controlPenaltyMPC_ubIdx04[12] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
controlPenaltyMPC_FLOAT* controlPenaltyMPC_lub04 = controlPenaltyMPC_l + 108;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_sub04 = controlPenaltyMPC_s + 108;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_lubbysub04 = controlPenaltyMPC_lbys + 108;
controlPenaltyMPC_FLOAT controlPenaltyMPC_riub04[12];
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dlubaff04 = controlPenaltyMPC_dl_aff + 108;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dsubaff04 = controlPenaltyMPC_ds_aff + 108;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dlubcc04 = controlPenaltyMPC_dl_cc + 108;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dsubcc04 = controlPenaltyMPC_ds_cc + 108;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_ccrhsub04 = controlPenaltyMPC_ccrhs + 108;
controlPenaltyMPC_FLOAT controlPenaltyMPC_Phi04[12];
controlPenaltyMPC_FLOAT controlPenaltyMPC_W04[12];
controlPenaltyMPC_FLOAT controlPenaltyMPC_Ysd04[36];
controlPenaltyMPC_FLOAT controlPenaltyMPC_Lsd04[36];
controlPenaltyMPC_FLOAT* controlPenaltyMPC_z05 = controlPenaltyMPC_z + 60;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dzaff05 = controlPenaltyMPC_dz_aff + 60;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dzcc05 = controlPenaltyMPC_dz_cc + 60;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_rd05 = controlPenaltyMPC_rd + 60;
controlPenaltyMPC_FLOAT controlPenaltyMPC_Lbyrd05[12];
controlPenaltyMPC_FLOAT* controlPenaltyMPC_grad_cost05 = controlPenaltyMPC_grad_cost + 60;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_grad_eq05 = controlPenaltyMPC_grad_eq + 60;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_grad_ineq05 = controlPenaltyMPC_grad_ineq + 60;
controlPenaltyMPC_FLOAT controlPenaltyMPC_ctv05[12];
controlPenaltyMPC_FLOAT* controlPenaltyMPC_v05 = controlPenaltyMPC_v + 30;
controlPenaltyMPC_FLOAT controlPenaltyMPC_re05[6];
controlPenaltyMPC_FLOAT controlPenaltyMPC_beta05[6];
controlPenaltyMPC_FLOAT controlPenaltyMPC_betacc05[6];
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dvaff05 = controlPenaltyMPC_dv_aff + 30;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dvcc05 = controlPenaltyMPC_dv_cc + 30;
controlPenaltyMPC_FLOAT controlPenaltyMPC_V05[72];
controlPenaltyMPC_FLOAT controlPenaltyMPC_Yd05[21];
controlPenaltyMPC_FLOAT controlPenaltyMPC_Ld05[21];
controlPenaltyMPC_FLOAT controlPenaltyMPC_yy05[6];
controlPenaltyMPC_FLOAT controlPenaltyMPC_bmy05[6];
controlPenaltyMPC_FLOAT controlPenaltyMPC_c05[6] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
int controlPenaltyMPC_lbIdx05[12] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
controlPenaltyMPC_FLOAT* controlPenaltyMPC_llb05 = controlPenaltyMPC_l + 120;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_slb05 = controlPenaltyMPC_s + 120;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_llbbyslb05 = controlPenaltyMPC_lbys + 120;
controlPenaltyMPC_FLOAT controlPenaltyMPC_rilb05[12];
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dllbaff05 = controlPenaltyMPC_dl_aff + 120;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dslbaff05 = controlPenaltyMPC_ds_aff + 120;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dllbcc05 = controlPenaltyMPC_dl_cc + 120;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dslbcc05 = controlPenaltyMPC_ds_cc + 120;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_ccrhsl05 = controlPenaltyMPC_ccrhs + 120;
int controlPenaltyMPC_ubIdx05[12] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
controlPenaltyMPC_FLOAT* controlPenaltyMPC_lub05 = controlPenaltyMPC_l + 132;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_sub05 = controlPenaltyMPC_s + 132;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_lubbysub05 = controlPenaltyMPC_lbys + 132;
controlPenaltyMPC_FLOAT controlPenaltyMPC_riub05[12];
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dlubaff05 = controlPenaltyMPC_dl_aff + 132;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dsubaff05 = controlPenaltyMPC_ds_aff + 132;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dlubcc05 = controlPenaltyMPC_dl_cc + 132;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dsubcc05 = controlPenaltyMPC_ds_cc + 132;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_ccrhsub05 = controlPenaltyMPC_ccrhs + 132;
controlPenaltyMPC_FLOAT controlPenaltyMPC_Phi05[12];
controlPenaltyMPC_FLOAT controlPenaltyMPC_W05[12];
controlPenaltyMPC_FLOAT controlPenaltyMPC_Ysd05[36];
controlPenaltyMPC_FLOAT controlPenaltyMPC_Lsd05[36];
controlPenaltyMPC_FLOAT* controlPenaltyMPC_z06 = controlPenaltyMPC_z + 72;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dzaff06 = controlPenaltyMPC_dz_aff + 72;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dzcc06 = controlPenaltyMPC_dz_cc + 72;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_rd06 = controlPenaltyMPC_rd + 72;
controlPenaltyMPC_FLOAT controlPenaltyMPC_Lbyrd06[12];
controlPenaltyMPC_FLOAT* controlPenaltyMPC_grad_cost06 = controlPenaltyMPC_grad_cost + 72;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_grad_eq06 = controlPenaltyMPC_grad_eq + 72;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_grad_ineq06 = controlPenaltyMPC_grad_ineq + 72;
controlPenaltyMPC_FLOAT controlPenaltyMPC_ctv06[12];
controlPenaltyMPC_FLOAT* controlPenaltyMPC_v06 = controlPenaltyMPC_v + 36;
controlPenaltyMPC_FLOAT controlPenaltyMPC_re06[6];
controlPenaltyMPC_FLOAT controlPenaltyMPC_beta06[6];
controlPenaltyMPC_FLOAT controlPenaltyMPC_betacc06[6];
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dvaff06 = controlPenaltyMPC_dv_aff + 36;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dvcc06 = controlPenaltyMPC_dv_cc + 36;
controlPenaltyMPC_FLOAT controlPenaltyMPC_V06[72];
controlPenaltyMPC_FLOAT controlPenaltyMPC_Yd06[21];
controlPenaltyMPC_FLOAT controlPenaltyMPC_Ld06[21];
controlPenaltyMPC_FLOAT controlPenaltyMPC_yy06[6];
controlPenaltyMPC_FLOAT controlPenaltyMPC_bmy06[6];
controlPenaltyMPC_FLOAT controlPenaltyMPC_c06[6] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
int controlPenaltyMPC_lbIdx06[12] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
controlPenaltyMPC_FLOAT* controlPenaltyMPC_llb06 = controlPenaltyMPC_l + 144;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_slb06 = controlPenaltyMPC_s + 144;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_llbbyslb06 = controlPenaltyMPC_lbys + 144;
controlPenaltyMPC_FLOAT controlPenaltyMPC_rilb06[12];
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dllbaff06 = controlPenaltyMPC_dl_aff + 144;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dslbaff06 = controlPenaltyMPC_ds_aff + 144;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dllbcc06 = controlPenaltyMPC_dl_cc + 144;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dslbcc06 = controlPenaltyMPC_ds_cc + 144;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_ccrhsl06 = controlPenaltyMPC_ccrhs + 144;
int controlPenaltyMPC_ubIdx06[12] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
controlPenaltyMPC_FLOAT* controlPenaltyMPC_lub06 = controlPenaltyMPC_l + 156;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_sub06 = controlPenaltyMPC_s + 156;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_lubbysub06 = controlPenaltyMPC_lbys + 156;
controlPenaltyMPC_FLOAT controlPenaltyMPC_riub06[12];
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dlubaff06 = controlPenaltyMPC_dl_aff + 156;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dsubaff06 = controlPenaltyMPC_ds_aff + 156;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dlubcc06 = controlPenaltyMPC_dl_cc + 156;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dsubcc06 = controlPenaltyMPC_ds_cc + 156;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_ccrhsub06 = controlPenaltyMPC_ccrhs + 156;
controlPenaltyMPC_FLOAT controlPenaltyMPC_Phi06[12];
controlPenaltyMPC_FLOAT controlPenaltyMPC_W06[12];
controlPenaltyMPC_FLOAT controlPenaltyMPC_Ysd06[36];
controlPenaltyMPC_FLOAT controlPenaltyMPC_Lsd06[36];
controlPenaltyMPC_FLOAT* controlPenaltyMPC_z07 = controlPenaltyMPC_z + 84;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dzaff07 = controlPenaltyMPC_dz_aff + 84;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dzcc07 = controlPenaltyMPC_dz_cc + 84;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_rd07 = controlPenaltyMPC_rd + 84;
controlPenaltyMPC_FLOAT controlPenaltyMPC_Lbyrd07[12];
controlPenaltyMPC_FLOAT* controlPenaltyMPC_grad_cost07 = controlPenaltyMPC_grad_cost + 84;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_grad_eq07 = controlPenaltyMPC_grad_eq + 84;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_grad_ineq07 = controlPenaltyMPC_grad_ineq + 84;
controlPenaltyMPC_FLOAT controlPenaltyMPC_ctv07[12];
controlPenaltyMPC_FLOAT* controlPenaltyMPC_v07 = controlPenaltyMPC_v + 42;
controlPenaltyMPC_FLOAT controlPenaltyMPC_re07[6];
controlPenaltyMPC_FLOAT controlPenaltyMPC_beta07[6];
controlPenaltyMPC_FLOAT controlPenaltyMPC_betacc07[6];
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dvaff07 = controlPenaltyMPC_dv_aff + 42;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dvcc07 = controlPenaltyMPC_dv_cc + 42;
controlPenaltyMPC_FLOAT controlPenaltyMPC_V07[72];
controlPenaltyMPC_FLOAT controlPenaltyMPC_Yd07[21];
controlPenaltyMPC_FLOAT controlPenaltyMPC_Ld07[21];
controlPenaltyMPC_FLOAT controlPenaltyMPC_yy07[6];
controlPenaltyMPC_FLOAT controlPenaltyMPC_bmy07[6];
controlPenaltyMPC_FLOAT controlPenaltyMPC_c07[6] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
int controlPenaltyMPC_lbIdx07[12] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
controlPenaltyMPC_FLOAT* controlPenaltyMPC_llb07 = controlPenaltyMPC_l + 168;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_slb07 = controlPenaltyMPC_s + 168;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_llbbyslb07 = controlPenaltyMPC_lbys + 168;
controlPenaltyMPC_FLOAT controlPenaltyMPC_rilb07[12];
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dllbaff07 = controlPenaltyMPC_dl_aff + 168;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dslbaff07 = controlPenaltyMPC_ds_aff + 168;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dllbcc07 = controlPenaltyMPC_dl_cc + 168;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dslbcc07 = controlPenaltyMPC_ds_cc + 168;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_ccrhsl07 = controlPenaltyMPC_ccrhs + 168;
int controlPenaltyMPC_ubIdx07[12] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
controlPenaltyMPC_FLOAT* controlPenaltyMPC_lub07 = controlPenaltyMPC_l + 180;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_sub07 = controlPenaltyMPC_s + 180;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_lubbysub07 = controlPenaltyMPC_lbys + 180;
controlPenaltyMPC_FLOAT controlPenaltyMPC_riub07[12];
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dlubaff07 = controlPenaltyMPC_dl_aff + 180;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dsubaff07 = controlPenaltyMPC_ds_aff + 180;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dlubcc07 = controlPenaltyMPC_dl_cc + 180;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dsubcc07 = controlPenaltyMPC_ds_cc + 180;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_ccrhsub07 = controlPenaltyMPC_ccrhs + 180;
controlPenaltyMPC_FLOAT controlPenaltyMPC_Phi07[12];
controlPenaltyMPC_FLOAT controlPenaltyMPC_W07[12];
controlPenaltyMPC_FLOAT controlPenaltyMPC_Ysd07[36];
controlPenaltyMPC_FLOAT controlPenaltyMPC_Lsd07[36];
controlPenaltyMPC_FLOAT* controlPenaltyMPC_z08 = controlPenaltyMPC_z + 96;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dzaff08 = controlPenaltyMPC_dz_aff + 96;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dzcc08 = controlPenaltyMPC_dz_cc + 96;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_rd08 = controlPenaltyMPC_rd + 96;
controlPenaltyMPC_FLOAT controlPenaltyMPC_Lbyrd08[12];
controlPenaltyMPC_FLOAT* controlPenaltyMPC_grad_cost08 = controlPenaltyMPC_grad_cost + 96;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_grad_eq08 = controlPenaltyMPC_grad_eq + 96;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_grad_ineq08 = controlPenaltyMPC_grad_ineq + 96;
controlPenaltyMPC_FLOAT controlPenaltyMPC_ctv08[12];
controlPenaltyMPC_FLOAT* controlPenaltyMPC_v08 = controlPenaltyMPC_v + 48;
controlPenaltyMPC_FLOAT controlPenaltyMPC_re08[6];
controlPenaltyMPC_FLOAT controlPenaltyMPC_beta08[6];
controlPenaltyMPC_FLOAT controlPenaltyMPC_betacc08[6];
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dvaff08 = controlPenaltyMPC_dv_aff + 48;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dvcc08 = controlPenaltyMPC_dv_cc + 48;
controlPenaltyMPC_FLOAT controlPenaltyMPC_V08[72];
controlPenaltyMPC_FLOAT controlPenaltyMPC_Yd08[21];
controlPenaltyMPC_FLOAT controlPenaltyMPC_Ld08[21];
controlPenaltyMPC_FLOAT controlPenaltyMPC_yy08[6];
controlPenaltyMPC_FLOAT controlPenaltyMPC_bmy08[6];
controlPenaltyMPC_FLOAT controlPenaltyMPC_c08[6] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
int controlPenaltyMPC_lbIdx08[12] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
controlPenaltyMPC_FLOAT* controlPenaltyMPC_llb08 = controlPenaltyMPC_l + 192;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_slb08 = controlPenaltyMPC_s + 192;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_llbbyslb08 = controlPenaltyMPC_lbys + 192;
controlPenaltyMPC_FLOAT controlPenaltyMPC_rilb08[12];
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dllbaff08 = controlPenaltyMPC_dl_aff + 192;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dslbaff08 = controlPenaltyMPC_ds_aff + 192;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dllbcc08 = controlPenaltyMPC_dl_cc + 192;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dslbcc08 = controlPenaltyMPC_ds_cc + 192;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_ccrhsl08 = controlPenaltyMPC_ccrhs + 192;
int controlPenaltyMPC_ubIdx08[12] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
controlPenaltyMPC_FLOAT* controlPenaltyMPC_lub08 = controlPenaltyMPC_l + 204;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_sub08 = controlPenaltyMPC_s + 204;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_lubbysub08 = controlPenaltyMPC_lbys + 204;
controlPenaltyMPC_FLOAT controlPenaltyMPC_riub08[12];
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dlubaff08 = controlPenaltyMPC_dl_aff + 204;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dsubaff08 = controlPenaltyMPC_ds_aff + 204;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dlubcc08 = controlPenaltyMPC_dl_cc + 204;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dsubcc08 = controlPenaltyMPC_ds_cc + 204;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_ccrhsub08 = controlPenaltyMPC_ccrhs + 204;
controlPenaltyMPC_FLOAT controlPenaltyMPC_Phi08[12];
controlPenaltyMPC_FLOAT controlPenaltyMPC_W08[12];
controlPenaltyMPC_FLOAT controlPenaltyMPC_Ysd08[36];
controlPenaltyMPC_FLOAT controlPenaltyMPC_Lsd08[36];
controlPenaltyMPC_FLOAT* controlPenaltyMPC_z09 = controlPenaltyMPC_z + 108;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dzaff09 = controlPenaltyMPC_dz_aff + 108;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dzcc09 = controlPenaltyMPC_dz_cc + 108;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_rd09 = controlPenaltyMPC_rd + 108;
controlPenaltyMPC_FLOAT controlPenaltyMPC_Lbyrd09[12];
controlPenaltyMPC_FLOAT* controlPenaltyMPC_grad_cost09 = controlPenaltyMPC_grad_cost + 108;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_grad_eq09 = controlPenaltyMPC_grad_eq + 108;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_grad_ineq09 = controlPenaltyMPC_grad_ineq + 108;
controlPenaltyMPC_FLOAT controlPenaltyMPC_ctv09[12];
controlPenaltyMPC_FLOAT* controlPenaltyMPC_v09 = controlPenaltyMPC_v + 54;
controlPenaltyMPC_FLOAT controlPenaltyMPC_re09[6];
controlPenaltyMPC_FLOAT controlPenaltyMPC_beta09[6];
controlPenaltyMPC_FLOAT controlPenaltyMPC_betacc09[6];
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dvaff09 = controlPenaltyMPC_dv_aff + 54;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dvcc09 = controlPenaltyMPC_dv_cc + 54;
controlPenaltyMPC_FLOAT controlPenaltyMPC_V09[72];
controlPenaltyMPC_FLOAT controlPenaltyMPC_Yd09[21];
controlPenaltyMPC_FLOAT controlPenaltyMPC_Ld09[21];
controlPenaltyMPC_FLOAT controlPenaltyMPC_yy09[6];
controlPenaltyMPC_FLOAT controlPenaltyMPC_bmy09[6];
controlPenaltyMPC_FLOAT controlPenaltyMPC_c09[6] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
int controlPenaltyMPC_lbIdx09[12] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
controlPenaltyMPC_FLOAT* controlPenaltyMPC_llb09 = controlPenaltyMPC_l + 216;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_slb09 = controlPenaltyMPC_s + 216;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_llbbyslb09 = controlPenaltyMPC_lbys + 216;
controlPenaltyMPC_FLOAT controlPenaltyMPC_rilb09[12];
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dllbaff09 = controlPenaltyMPC_dl_aff + 216;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dslbaff09 = controlPenaltyMPC_ds_aff + 216;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dllbcc09 = controlPenaltyMPC_dl_cc + 216;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dslbcc09 = controlPenaltyMPC_ds_cc + 216;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_ccrhsl09 = controlPenaltyMPC_ccrhs + 216;
int controlPenaltyMPC_ubIdx09[12] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
controlPenaltyMPC_FLOAT* controlPenaltyMPC_lub09 = controlPenaltyMPC_l + 228;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_sub09 = controlPenaltyMPC_s + 228;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_lubbysub09 = controlPenaltyMPC_lbys + 228;
controlPenaltyMPC_FLOAT controlPenaltyMPC_riub09[12];
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dlubaff09 = controlPenaltyMPC_dl_aff + 228;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dsubaff09 = controlPenaltyMPC_ds_aff + 228;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dlubcc09 = controlPenaltyMPC_dl_cc + 228;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dsubcc09 = controlPenaltyMPC_ds_cc + 228;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_ccrhsub09 = controlPenaltyMPC_ccrhs + 228;
controlPenaltyMPC_FLOAT controlPenaltyMPC_Phi09[12];
controlPenaltyMPC_FLOAT controlPenaltyMPC_W09[12];
controlPenaltyMPC_FLOAT controlPenaltyMPC_Ysd09[36];
controlPenaltyMPC_FLOAT controlPenaltyMPC_Lsd09[36];
controlPenaltyMPC_FLOAT* controlPenaltyMPC_z10 = controlPenaltyMPC_z + 120;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dzaff10 = controlPenaltyMPC_dz_aff + 120;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dzcc10 = controlPenaltyMPC_dz_cc + 120;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_rd10 = controlPenaltyMPC_rd + 120;
controlPenaltyMPC_FLOAT controlPenaltyMPC_Lbyrd10[12];
controlPenaltyMPC_FLOAT* controlPenaltyMPC_grad_cost10 = controlPenaltyMPC_grad_cost + 120;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_grad_eq10 = controlPenaltyMPC_grad_eq + 120;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_grad_ineq10 = controlPenaltyMPC_grad_ineq + 120;
controlPenaltyMPC_FLOAT controlPenaltyMPC_ctv10[12];
controlPenaltyMPC_FLOAT* controlPenaltyMPC_v10 = controlPenaltyMPC_v + 60;
controlPenaltyMPC_FLOAT controlPenaltyMPC_re10[6];
controlPenaltyMPC_FLOAT controlPenaltyMPC_beta10[6];
controlPenaltyMPC_FLOAT controlPenaltyMPC_betacc10[6];
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dvaff10 = controlPenaltyMPC_dv_aff + 60;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dvcc10 = controlPenaltyMPC_dv_cc + 60;
controlPenaltyMPC_FLOAT controlPenaltyMPC_V10[72];
controlPenaltyMPC_FLOAT controlPenaltyMPC_Yd10[21];
controlPenaltyMPC_FLOAT controlPenaltyMPC_Ld10[21];
controlPenaltyMPC_FLOAT controlPenaltyMPC_yy10[6];
controlPenaltyMPC_FLOAT controlPenaltyMPC_bmy10[6];
controlPenaltyMPC_FLOAT controlPenaltyMPC_c10[6] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
int controlPenaltyMPC_lbIdx10[12] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
controlPenaltyMPC_FLOAT* controlPenaltyMPC_llb10 = controlPenaltyMPC_l + 240;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_slb10 = controlPenaltyMPC_s + 240;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_llbbyslb10 = controlPenaltyMPC_lbys + 240;
controlPenaltyMPC_FLOAT controlPenaltyMPC_rilb10[12];
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dllbaff10 = controlPenaltyMPC_dl_aff + 240;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dslbaff10 = controlPenaltyMPC_ds_aff + 240;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dllbcc10 = controlPenaltyMPC_dl_cc + 240;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dslbcc10 = controlPenaltyMPC_ds_cc + 240;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_ccrhsl10 = controlPenaltyMPC_ccrhs + 240;
int controlPenaltyMPC_ubIdx10[12] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
controlPenaltyMPC_FLOAT* controlPenaltyMPC_lub10 = controlPenaltyMPC_l + 252;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_sub10 = controlPenaltyMPC_s + 252;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_lubbysub10 = controlPenaltyMPC_lbys + 252;
controlPenaltyMPC_FLOAT controlPenaltyMPC_riub10[12];
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dlubaff10 = controlPenaltyMPC_dl_aff + 252;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dsubaff10 = controlPenaltyMPC_ds_aff + 252;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dlubcc10 = controlPenaltyMPC_dl_cc + 252;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dsubcc10 = controlPenaltyMPC_ds_cc + 252;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_ccrhsub10 = controlPenaltyMPC_ccrhs + 252;
controlPenaltyMPC_FLOAT controlPenaltyMPC_Phi10[12];
controlPenaltyMPC_FLOAT controlPenaltyMPC_W10[12];
controlPenaltyMPC_FLOAT controlPenaltyMPC_Ysd10[36];
controlPenaltyMPC_FLOAT controlPenaltyMPC_Lsd10[36];
controlPenaltyMPC_FLOAT* controlPenaltyMPC_z11 = controlPenaltyMPC_z + 132;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dzaff11 = controlPenaltyMPC_dz_aff + 132;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dzcc11 = controlPenaltyMPC_dz_cc + 132;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_rd11 = controlPenaltyMPC_rd + 132;
controlPenaltyMPC_FLOAT controlPenaltyMPC_Lbyrd11[12];
controlPenaltyMPC_FLOAT* controlPenaltyMPC_grad_cost11 = controlPenaltyMPC_grad_cost + 132;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_grad_eq11 = controlPenaltyMPC_grad_eq + 132;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_grad_ineq11 = controlPenaltyMPC_grad_ineq + 132;
controlPenaltyMPC_FLOAT controlPenaltyMPC_ctv11[12];
controlPenaltyMPC_FLOAT* controlPenaltyMPC_v11 = controlPenaltyMPC_v + 66;
controlPenaltyMPC_FLOAT controlPenaltyMPC_re11[6];
controlPenaltyMPC_FLOAT controlPenaltyMPC_beta11[6];
controlPenaltyMPC_FLOAT controlPenaltyMPC_betacc11[6];
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dvaff11 = controlPenaltyMPC_dv_aff + 66;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dvcc11 = controlPenaltyMPC_dv_cc + 66;
controlPenaltyMPC_FLOAT controlPenaltyMPC_V11[72];
controlPenaltyMPC_FLOAT controlPenaltyMPC_Yd11[21];
controlPenaltyMPC_FLOAT controlPenaltyMPC_Ld11[21];
controlPenaltyMPC_FLOAT controlPenaltyMPC_yy11[6];
controlPenaltyMPC_FLOAT controlPenaltyMPC_bmy11[6];
controlPenaltyMPC_FLOAT controlPenaltyMPC_c11[6] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
int controlPenaltyMPC_lbIdx11[12] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
controlPenaltyMPC_FLOAT* controlPenaltyMPC_llb11 = controlPenaltyMPC_l + 264;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_slb11 = controlPenaltyMPC_s + 264;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_llbbyslb11 = controlPenaltyMPC_lbys + 264;
controlPenaltyMPC_FLOAT controlPenaltyMPC_rilb11[12];
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dllbaff11 = controlPenaltyMPC_dl_aff + 264;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dslbaff11 = controlPenaltyMPC_ds_aff + 264;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dllbcc11 = controlPenaltyMPC_dl_cc + 264;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dslbcc11 = controlPenaltyMPC_ds_cc + 264;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_ccrhsl11 = controlPenaltyMPC_ccrhs + 264;
int controlPenaltyMPC_ubIdx11[12] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
controlPenaltyMPC_FLOAT* controlPenaltyMPC_lub11 = controlPenaltyMPC_l + 276;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_sub11 = controlPenaltyMPC_s + 276;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_lubbysub11 = controlPenaltyMPC_lbys + 276;
controlPenaltyMPC_FLOAT controlPenaltyMPC_riub11[12];
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dlubaff11 = controlPenaltyMPC_dl_aff + 276;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dsubaff11 = controlPenaltyMPC_ds_aff + 276;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dlubcc11 = controlPenaltyMPC_dl_cc + 276;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dsubcc11 = controlPenaltyMPC_ds_cc + 276;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_ccrhsub11 = controlPenaltyMPC_ccrhs + 276;
controlPenaltyMPC_FLOAT controlPenaltyMPC_Phi11[12];
controlPenaltyMPC_FLOAT controlPenaltyMPC_W11[12];
controlPenaltyMPC_FLOAT controlPenaltyMPC_Ysd11[36];
controlPenaltyMPC_FLOAT controlPenaltyMPC_Lsd11[36];
controlPenaltyMPC_FLOAT* controlPenaltyMPC_z12 = controlPenaltyMPC_z + 144;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dzaff12 = controlPenaltyMPC_dz_aff + 144;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dzcc12 = controlPenaltyMPC_dz_cc + 144;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_rd12 = controlPenaltyMPC_rd + 144;
controlPenaltyMPC_FLOAT controlPenaltyMPC_Lbyrd12[12];
controlPenaltyMPC_FLOAT* controlPenaltyMPC_grad_cost12 = controlPenaltyMPC_grad_cost + 144;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_grad_eq12 = controlPenaltyMPC_grad_eq + 144;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_grad_ineq12 = controlPenaltyMPC_grad_ineq + 144;
controlPenaltyMPC_FLOAT controlPenaltyMPC_ctv12[12];
controlPenaltyMPC_FLOAT* controlPenaltyMPC_v12 = controlPenaltyMPC_v + 72;
controlPenaltyMPC_FLOAT controlPenaltyMPC_re12[6];
controlPenaltyMPC_FLOAT controlPenaltyMPC_beta12[6];
controlPenaltyMPC_FLOAT controlPenaltyMPC_betacc12[6];
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dvaff12 = controlPenaltyMPC_dv_aff + 72;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dvcc12 = controlPenaltyMPC_dv_cc + 72;
controlPenaltyMPC_FLOAT controlPenaltyMPC_V12[72];
controlPenaltyMPC_FLOAT controlPenaltyMPC_Yd12[21];
controlPenaltyMPC_FLOAT controlPenaltyMPC_Ld12[21];
controlPenaltyMPC_FLOAT controlPenaltyMPC_yy12[6];
controlPenaltyMPC_FLOAT controlPenaltyMPC_bmy12[6];
controlPenaltyMPC_FLOAT controlPenaltyMPC_c12[6] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
int controlPenaltyMPC_lbIdx12[12] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
controlPenaltyMPC_FLOAT* controlPenaltyMPC_llb12 = controlPenaltyMPC_l + 288;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_slb12 = controlPenaltyMPC_s + 288;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_llbbyslb12 = controlPenaltyMPC_lbys + 288;
controlPenaltyMPC_FLOAT controlPenaltyMPC_rilb12[12];
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dllbaff12 = controlPenaltyMPC_dl_aff + 288;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dslbaff12 = controlPenaltyMPC_ds_aff + 288;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dllbcc12 = controlPenaltyMPC_dl_cc + 288;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dslbcc12 = controlPenaltyMPC_ds_cc + 288;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_ccrhsl12 = controlPenaltyMPC_ccrhs + 288;
int controlPenaltyMPC_ubIdx12[12] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
controlPenaltyMPC_FLOAT* controlPenaltyMPC_lub12 = controlPenaltyMPC_l + 300;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_sub12 = controlPenaltyMPC_s + 300;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_lubbysub12 = controlPenaltyMPC_lbys + 300;
controlPenaltyMPC_FLOAT controlPenaltyMPC_riub12[12];
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dlubaff12 = controlPenaltyMPC_dl_aff + 300;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dsubaff12 = controlPenaltyMPC_ds_aff + 300;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dlubcc12 = controlPenaltyMPC_dl_cc + 300;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dsubcc12 = controlPenaltyMPC_ds_cc + 300;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_ccrhsub12 = controlPenaltyMPC_ccrhs + 300;
controlPenaltyMPC_FLOAT controlPenaltyMPC_Phi12[12];
controlPenaltyMPC_FLOAT controlPenaltyMPC_W12[12];
controlPenaltyMPC_FLOAT controlPenaltyMPC_Ysd12[36];
controlPenaltyMPC_FLOAT controlPenaltyMPC_Lsd12[36];
controlPenaltyMPC_FLOAT* controlPenaltyMPC_z13 = controlPenaltyMPC_z + 156;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dzaff13 = controlPenaltyMPC_dz_aff + 156;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dzcc13 = controlPenaltyMPC_dz_cc + 156;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_rd13 = controlPenaltyMPC_rd + 156;
controlPenaltyMPC_FLOAT controlPenaltyMPC_Lbyrd13[12];
controlPenaltyMPC_FLOAT* controlPenaltyMPC_grad_cost13 = controlPenaltyMPC_grad_cost + 156;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_grad_eq13 = controlPenaltyMPC_grad_eq + 156;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_grad_ineq13 = controlPenaltyMPC_grad_ineq + 156;
controlPenaltyMPC_FLOAT controlPenaltyMPC_ctv13[12];
controlPenaltyMPC_FLOAT* controlPenaltyMPC_v13 = controlPenaltyMPC_v + 78;
controlPenaltyMPC_FLOAT controlPenaltyMPC_re13[6];
controlPenaltyMPC_FLOAT controlPenaltyMPC_beta13[6];
controlPenaltyMPC_FLOAT controlPenaltyMPC_betacc13[6];
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dvaff13 = controlPenaltyMPC_dv_aff + 78;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dvcc13 = controlPenaltyMPC_dv_cc + 78;
controlPenaltyMPC_FLOAT controlPenaltyMPC_V13[72];
controlPenaltyMPC_FLOAT controlPenaltyMPC_Yd13[21];
controlPenaltyMPC_FLOAT controlPenaltyMPC_Ld13[21];
controlPenaltyMPC_FLOAT controlPenaltyMPC_yy13[6];
controlPenaltyMPC_FLOAT controlPenaltyMPC_bmy13[6];
controlPenaltyMPC_FLOAT controlPenaltyMPC_c13[6] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
int controlPenaltyMPC_lbIdx13[12] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
controlPenaltyMPC_FLOAT* controlPenaltyMPC_llb13 = controlPenaltyMPC_l + 312;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_slb13 = controlPenaltyMPC_s + 312;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_llbbyslb13 = controlPenaltyMPC_lbys + 312;
controlPenaltyMPC_FLOAT controlPenaltyMPC_rilb13[12];
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dllbaff13 = controlPenaltyMPC_dl_aff + 312;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dslbaff13 = controlPenaltyMPC_ds_aff + 312;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dllbcc13 = controlPenaltyMPC_dl_cc + 312;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dslbcc13 = controlPenaltyMPC_ds_cc + 312;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_ccrhsl13 = controlPenaltyMPC_ccrhs + 312;
int controlPenaltyMPC_ubIdx13[12] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
controlPenaltyMPC_FLOAT* controlPenaltyMPC_lub13 = controlPenaltyMPC_l + 324;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_sub13 = controlPenaltyMPC_s + 324;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_lubbysub13 = controlPenaltyMPC_lbys + 324;
controlPenaltyMPC_FLOAT controlPenaltyMPC_riub13[12];
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dlubaff13 = controlPenaltyMPC_dl_aff + 324;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dsubaff13 = controlPenaltyMPC_ds_aff + 324;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dlubcc13 = controlPenaltyMPC_dl_cc + 324;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dsubcc13 = controlPenaltyMPC_ds_cc + 324;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_ccrhsub13 = controlPenaltyMPC_ccrhs + 324;
controlPenaltyMPC_FLOAT controlPenaltyMPC_Phi13[12];
controlPenaltyMPC_FLOAT controlPenaltyMPC_W13[12];
controlPenaltyMPC_FLOAT controlPenaltyMPC_Ysd13[36];
controlPenaltyMPC_FLOAT controlPenaltyMPC_Lsd13[36];
controlPenaltyMPC_FLOAT* controlPenaltyMPC_z14 = controlPenaltyMPC_z + 168;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dzaff14 = controlPenaltyMPC_dz_aff + 168;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dzcc14 = controlPenaltyMPC_dz_cc + 168;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_rd14 = controlPenaltyMPC_rd + 168;
controlPenaltyMPC_FLOAT controlPenaltyMPC_Lbyrd14[12];
controlPenaltyMPC_FLOAT* controlPenaltyMPC_grad_cost14 = controlPenaltyMPC_grad_cost + 168;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_grad_eq14 = controlPenaltyMPC_grad_eq + 168;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_grad_ineq14 = controlPenaltyMPC_grad_ineq + 168;
controlPenaltyMPC_FLOAT controlPenaltyMPC_ctv14[12];
controlPenaltyMPC_FLOAT* controlPenaltyMPC_v14 = controlPenaltyMPC_v + 84;
controlPenaltyMPC_FLOAT controlPenaltyMPC_re14[6];
controlPenaltyMPC_FLOAT controlPenaltyMPC_beta14[6];
controlPenaltyMPC_FLOAT controlPenaltyMPC_betacc14[6];
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dvaff14 = controlPenaltyMPC_dv_aff + 84;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dvcc14 = controlPenaltyMPC_dv_cc + 84;
controlPenaltyMPC_FLOAT controlPenaltyMPC_V14[72];
controlPenaltyMPC_FLOAT controlPenaltyMPC_Yd14[21];
controlPenaltyMPC_FLOAT controlPenaltyMPC_Ld14[21];
controlPenaltyMPC_FLOAT controlPenaltyMPC_yy14[6];
controlPenaltyMPC_FLOAT controlPenaltyMPC_bmy14[6];
controlPenaltyMPC_FLOAT controlPenaltyMPC_c14[6] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
int controlPenaltyMPC_lbIdx14[12] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
controlPenaltyMPC_FLOAT* controlPenaltyMPC_llb14 = controlPenaltyMPC_l + 336;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_slb14 = controlPenaltyMPC_s + 336;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_llbbyslb14 = controlPenaltyMPC_lbys + 336;
controlPenaltyMPC_FLOAT controlPenaltyMPC_rilb14[12];
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dllbaff14 = controlPenaltyMPC_dl_aff + 336;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dslbaff14 = controlPenaltyMPC_ds_aff + 336;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dllbcc14 = controlPenaltyMPC_dl_cc + 336;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dslbcc14 = controlPenaltyMPC_ds_cc + 336;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_ccrhsl14 = controlPenaltyMPC_ccrhs + 336;
int controlPenaltyMPC_ubIdx14[6] = {0, 1, 2, 3, 4, 5};
controlPenaltyMPC_FLOAT* controlPenaltyMPC_lub14 = controlPenaltyMPC_l + 348;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_sub14 = controlPenaltyMPC_s + 348;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_lubbysub14 = controlPenaltyMPC_lbys + 348;
controlPenaltyMPC_FLOAT controlPenaltyMPC_riub14[6];
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dlubaff14 = controlPenaltyMPC_dl_aff + 348;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dsubaff14 = controlPenaltyMPC_ds_aff + 348;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dlubcc14 = controlPenaltyMPC_dl_cc + 348;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_dsubcc14 = controlPenaltyMPC_ds_cc + 348;
controlPenaltyMPC_FLOAT* controlPenaltyMPC_ccrhsub14 = controlPenaltyMPC_ccrhs + 348;
controlPenaltyMPC_FLOAT controlPenaltyMPC_Phi14[12];
controlPenaltyMPC_FLOAT controlPenaltyMPC_W14[12];
controlPenaltyMPC_FLOAT controlPenaltyMPC_Ysd14[36];
controlPenaltyMPC_FLOAT controlPenaltyMPC_Lsd14[36];
controlPenaltyMPC_FLOAT musigma;
controlPenaltyMPC_FLOAT sigma_3rdroot;
controlPenaltyMPC_FLOAT controlPenaltyMPC_Diag1_0[12];
controlPenaltyMPC_FLOAT controlPenaltyMPC_Diag2_0[12];
controlPenaltyMPC_FLOAT controlPenaltyMPC_L_0[66];




/* SOLVER CODE --------------------------------------------------------- */
int controlPenaltyMPC_solve(controlPenaltyMPC_params* params, controlPenaltyMPC_output* output, controlPenaltyMPC_info* info)
{	
int exitcode;

#if controlPenaltyMPC_SET_TIMING == 1
	controlPenaltyMPC_timer solvertimer;
	controlPenaltyMPC_tic(&solvertimer);
#endif
/* FUNCTION CALLS INTO LA LIBRARY -------------------------------------- */
info->it = 0;
controlPenaltyMPC_LA_INITIALIZEVECTOR_180(controlPenaltyMPC_z, 0);
controlPenaltyMPC_LA_INITIALIZEVECTOR_90(controlPenaltyMPC_v, 1);
controlPenaltyMPC_LA_INITIALIZEVECTOR_354(controlPenaltyMPC_l, 1);
controlPenaltyMPC_LA_INITIALIZEVECTOR_354(controlPenaltyMPC_s, 1);
info->mu = 0;
controlPenaltyMPC_LA_DOTACC_354(controlPenaltyMPC_l, controlPenaltyMPC_s, &info->mu);
info->mu /= 354;
while( 1 ){
info->pobj = 0;
controlPenaltyMPC_LA_DIAG_QUADFCN_12(params->Q1, params->f1, controlPenaltyMPC_z00, controlPenaltyMPC_grad_cost00, &info->pobj);
controlPenaltyMPC_LA_DIAG_QUADFCN_12(params->Q2, params->f2, controlPenaltyMPC_z01, controlPenaltyMPC_grad_cost01, &info->pobj);
controlPenaltyMPC_LA_DIAG_QUADFCN_12(params->Q3, params->f3, controlPenaltyMPC_z02, controlPenaltyMPC_grad_cost02, &info->pobj);
controlPenaltyMPC_LA_DIAG_QUADFCN_12(params->Q4, params->f4, controlPenaltyMPC_z03, controlPenaltyMPC_grad_cost03, &info->pobj);
controlPenaltyMPC_LA_DIAG_QUADFCN_12(params->Q5, params->f5, controlPenaltyMPC_z04, controlPenaltyMPC_grad_cost04, &info->pobj);
controlPenaltyMPC_LA_DIAG_QUADFCN_12(params->Q6, params->f6, controlPenaltyMPC_z05, controlPenaltyMPC_grad_cost05, &info->pobj);
controlPenaltyMPC_LA_DIAG_QUADFCN_12(params->Q7, params->f7, controlPenaltyMPC_z06, controlPenaltyMPC_grad_cost06, &info->pobj);
controlPenaltyMPC_LA_DIAG_QUADFCN_12(params->Q8, params->f8, controlPenaltyMPC_z07, controlPenaltyMPC_grad_cost07, &info->pobj);
controlPenaltyMPC_LA_DIAG_QUADFCN_12(params->Q9, params->f9, controlPenaltyMPC_z08, controlPenaltyMPC_grad_cost08, &info->pobj);
controlPenaltyMPC_LA_DIAG_QUADFCN_12(params->Q10, params->f10, controlPenaltyMPC_z09, controlPenaltyMPC_grad_cost09, &info->pobj);
controlPenaltyMPC_LA_DIAG_QUADFCN_12(params->Q11, params->f11, controlPenaltyMPC_z10, controlPenaltyMPC_grad_cost10, &info->pobj);
controlPenaltyMPC_LA_DIAG_QUADFCN_12(params->Q12, params->f12, controlPenaltyMPC_z11, controlPenaltyMPC_grad_cost11, &info->pobj);
controlPenaltyMPC_LA_DIAG_QUADFCN_12(params->Q13, params->f13, controlPenaltyMPC_z12, controlPenaltyMPC_grad_cost12, &info->pobj);
controlPenaltyMPC_LA_DIAG_QUADFCN_12(params->Q14, params->f14, controlPenaltyMPC_z13, controlPenaltyMPC_grad_cost13, &info->pobj);
controlPenaltyMPC_LA_DIAG_QUADFCN_12(params->Q15, params->f15, controlPenaltyMPC_z14, controlPenaltyMPC_grad_cost14, &info->pobj);
info->res_eq = 0;
info->dgap = 0;
controlPenaltyMPC_LA_DIAGZERO_MVMSUB6_6(controlPenaltyMPC_D00, controlPenaltyMPC_z00, params->e1, controlPenaltyMPC_v00, controlPenaltyMPC_re00, &info->dgap, &info->res_eq);
controlPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_6_12_12(controlPenaltyMPC_C00, controlPenaltyMPC_z00, controlPenaltyMPC_D01, controlPenaltyMPC_z01, controlPenaltyMPC_c01, controlPenaltyMPC_v01, controlPenaltyMPC_re01, &info->dgap, &info->res_eq);
controlPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_6_12_12(controlPenaltyMPC_C00, controlPenaltyMPC_z01, controlPenaltyMPC_D01, controlPenaltyMPC_z02, controlPenaltyMPC_c02, controlPenaltyMPC_v02, controlPenaltyMPC_re02, &info->dgap, &info->res_eq);
controlPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_6_12_12(controlPenaltyMPC_C00, controlPenaltyMPC_z02, controlPenaltyMPC_D01, controlPenaltyMPC_z03, controlPenaltyMPC_c03, controlPenaltyMPC_v03, controlPenaltyMPC_re03, &info->dgap, &info->res_eq);
controlPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_6_12_12(controlPenaltyMPC_C00, controlPenaltyMPC_z03, controlPenaltyMPC_D01, controlPenaltyMPC_z04, controlPenaltyMPC_c04, controlPenaltyMPC_v04, controlPenaltyMPC_re04, &info->dgap, &info->res_eq);
controlPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_6_12_12(controlPenaltyMPC_C00, controlPenaltyMPC_z04, controlPenaltyMPC_D01, controlPenaltyMPC_z05, controlPenaltyMPC_c05, controlPenaltyMPC_v05, controlPenaltyMPC_re05, &info->dgap, &info->res_eq);
controlPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_6_12_12(controlPenaltyMPC_C00, controlPenaltyMPC_z05, controlPenaltyMPC_D01, controlPenaltyMPC_z06, controlPenaltyMPC_c06, controlPenaltyMPC_v06, controlPenaltyMPC_re06, &info->dgap, &info->res_eq);
controlPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_6_12_12(controlPenaltyMPC_C00, controlPenaltyMPC_z06, controlPenaltyMPC_D01, controlPenaltyMPC_z07, controlPenaltyMPC_c07, controlPenaltyMPC_v07, controlPenaltyMPC_re07, &info->dgap, &info->res_eq);
controlPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_6_12_12(controlPenaltyMPC_C00, controlPenaltyMPC_z07, controlPenaltyMPC_D01, controlPenaltyMPC_z08, controlPenaltyMPC_c08, controlPenaltyMPC_v08, controlPenaltyMPC_re08, &info->dgap, &info->res_eq);
controlPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_6_12_12(controlPenaltyMPC_C00, controlPenaltyMPC_z08, controlPenaltyMPC_D01, controlPenaltyMPC_z09, controlPenaltyMPC_c09, controlPenaltyMPC_v09, controlPenaltyMPC_re09, &info->dgap, &info->res_eq);
controlPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_6_12_12(controlPenaltyMPC_C00, controlPenaltyMPC_z09, controlPenaltyMPC_D01, controlPenaltyMPC_z10, controlPenaltyMPC_c10, controlPenaltyMPC_v10, controlPenaltyMPC_re10, &info->dgap, &info->res_eq);
controlPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_6_12_12(controlPenaltyMPC_C00, controlPenaltyMPC_z10, controlPenaltyMPC_D01, controlPenaltyMPC_z11, controlPenaltyMPC_c11, controlPenaltyMPC_v11, controlPenaltyMPC_re11, &info->dgap, &info->res_eq);
controlPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_6_12_12(controlPenaltyMPC_C00, controlPenaltyMPC_z11, controlPenaltyMPC_D01, controlPenaltyMPC_z12, controlPenaltyMPC_c12, controlPenaltyMPC_v12, controlPenaltyMPC_re12, &info->dgap, &info->res_eq);
controlPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_6_12_12(controlPenaltyMPC_C00, controlPenaltyMPC_z12, controlPenaltyMPC_D01, controlPenaltyMPC_z13, controlPenaltyMPC_c13, controlPenaltyMPC_v13, controlPenaltyMPC_re13, &info->dgap, &info->res_eq);
controlPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_6_12_12(controlPenaltyMPC_C00, controlPenaltyMPC_z13, controlPenaltyMPC_D01, controlPenaltyMPC_z14, controlPenaltyMPC_c14, controlPenaltyMPC_v14, controlPenaltyMPC_re14, &info->dgap, &info->res_eq);
controlPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_6_12_6(controlPenaltyMPC_C00, controlPenaltyMPC_v01, controlPenaltyMPC_D00, controlPenaltyMPC_v00, controlPenaltyMPC_grad_eq00);
controlPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_6_12_6(controlPenaltyMPC_C00, controlPenaltyMPC_v02, controlPenaltyMPC_D01, controlPenaltyMPC_v01, controlPenaltyMPC_grad_eq01);
controlPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_6_12_6(controlPenaltyMPC_C00, controlPenaltyMPC_v03, controlPenaltyMPC_D01, controlPenaltyMPC_v02, controlPenaltyMPC_grad_eq02);
controlPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_6_12_6(controlPenaltyMPC_C00, controlPenaltyMPC_v04, controlPenaltyMPC_D01, controlPenaltyMPC_v03, controlPenaltyMPC_grad_eq03);
controlPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_6_12_6(controlPenaltyMPC_C00, controlPenaltyMPC_v05, controlPenaltyMPC_D01, controlPenaltyMPC_v04, controlPenaltyMPC_grad_eq04);
controlPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_6_12_6(controlPenaltyMPC_C00, controlPenaltyMPC_v06, controlPenaltyMPC_D01, controlPenaltyMPC_v05, controlPenaltyMPC_grad_eq05);
controlPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_6_12_6(controlPenaltyMPC_C00, controlPenaltyMPC_v07, controlPenaltyMPC_D01, controlPenaltyMPC_v06, controlPenaltyMPC_grad_eq06);
controlPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_6_12_6(controlPenaltyMPC_C00, controlPenaltyMPC_v08, controlPenaltyMPC_D01, controlPenaltyMPC_v07, controlPenaltyMPC_grad_eq07);
controlPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_6_12_6(controlPenaltyMPC_C00, controlPenaltyMPC_v09, controlPenaltyMPC_D01, controlPenaltyMPC_v08, controlPenaltyMPC_grad_eq08);
controlPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_6_12_6(controlPenaltyMPC_C00, controlPenaltyMPC_v10, controlPenaltyMPC_D01, controlPenaltyMPC_v09, controlPenaltyMPC_grad_eq09);
controlPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_6_12_6(controlPenaltyMPC_C00, controlPenaltyMPC_v11, controlPenaltyMPC_D01, controlPenaltyMPC_v10, controlPenaltyMPC_grad_eq10);
controlPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_6_12_6(controlPenaltyMPC_C00, controlPenaltyMPC_v12, controlPenaltyMPC_D01, controlPenaltyMPC_v11, controlPenaltyMPC_grad_eq11);
controlPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_6_12_6(controlPenaltyMPC_C00, controlPenaltyMPC_v13, controlPenaltyMPC_D01, controlPenaltyMPC_v12, controlPenaltyMPC_grad_eq12);
controlPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_6_12_6(controlPenaltyMPC_C00, controlPenaltyMPC_v14, controlPenaltyMPC_D01, controlPenaltyMPC_v13, controlPenaltyMPC_grad_eq13);
controlPenaltyMPC_LA_DIAGZERO_MTVM_6_12(controlPenaltyMPC_D01, controlPenaltyMPC_v14, controlPenaltyMPC_grad_eq14);
info->res_ineq = 0;
controlPenaltyMPC_LA_VSUBADD3_12(params->lb1, controlPenaltyMPC_z00, controlPenaltyMPC_lbIdx00, controlPenaltyMPC_llb00, controlPenaltyMPC_slb00, controlPenaltyMPC_rilb00, &info->dgap, &info->res_ineq);
controlPenaltyMPC_LA_VSUBADD2_12(controlPenaltyMPC_z00, controlPenaltyMPC_ubIdx00, params->ub1, controlPenaltyMPC_lub00, controlPenaltyMPC_sub00, controlPenaltyMPC_riub00, &info->dgap, &info->res_ineq);
controlPenaltyMPC_LA_VSUBADD3_12(params->lb2, controlPenaltyMPC_z01, controlPenaltyMPC_lbIdx01, controlPenaltyMPC_llb01, controlPenaltyMPC_slb01, controlPenaltyMPC_rilb01, &info->dgap, &info->res_ineq);
controlPenaltyMPC_LA_VSUBADD2_12(controlPenaltyMPC_z01, controlPenaltyMPC_ubIdx01, params->ub2, controlPenaltyMPC_lub01, controlPenaltyMPC_sub01, controlPenaltyMPC_riub01, &info->dgap, &info->res_ineq);
controlPenaltyMPC_LA_VSUBADD3_12(params->lb3, controlPenaltyMPC_z02, controlPenaltyMPC_lbIdx02, controlPenaltyMPC_llb02, controlPenaltyMPC_slb02, controlPenaltyMPC_rilb02, &info->dgap, &info->res_ineq);
controlPenaltyMPC_LA_VSUBADD2_12(controlPenaltyMPC_z02, controlPenaltyMPC_ubIdx02, params->ub3, controlPenaltyMPC_lub02, controlPenaltyMPC_sub02, controlPenaltyMPC_riub02, &info->dgap, &info->res_ineq);
controlPenaltyMPC_LA_VSUBADD3_12(params->lb4, controlPenaltyMPC_z03, controlPenaltyMPC_lbIdx03, controlPenaltyMPC_llb03, controlPenaltyMPC_slb03, controlPenaltyMPC_rilb03, &info->dgap, &info->res_ineq);
controlPenaltyMPC_LA_VSUBADD2_12(controlPenaltyMPC_z03, controlPenaltyMPC_ubIdx03, params->ub4, controlPenaltyMPC_lub03, controlPenaltyMPC_sub03, controlPenaltyMPC_riub03, &info->dgap, &info->res_ineq);
controlPenaltyMPC_LA_VSUBADD3_12(params->lb5, controlPenaltyMPC_z04, controlPenaltyMPC_lbIdx04, controlPenaltyMPC_llb04, controlPenaltyMPC_slb04, controlPenaltyMPC_rilb04, &info->dgap, &info->res_ineq);
controlPenaltyMPC_LA_VSUBADD2_12(controlPenaltyMPC_z04, controlPenaltyMPC_ubIdx04, params->ub5, controlPenaltyMPC_lub04, controlPenaltyMPC_sub04, controlPenaltyMPC_riub04, &info->dgap, &info->res_ineq);
controlPenaltyMPC_LA_VSUBADD3_12(params->lb6, controlPenaltyMPC_z05, controlPenaltyMPC_lbIdx05, controlPenaltyMPC_llb05, controlPenaltyMPC_slb05, controlPenaltyMPC_rilb05, &info->dgap, &info->res_ineq);
controlPenaltyMPC_LA_VSUBADD2_12(controlPenaltyMPC_z05, controlPenaltyMPC_ubIdx05, params->ub6, controlPenaltyMPC_lub05, controlPenaltyMPC_sub05, controlPenaltyMPC_riub05, &info->dgap, &info->res_ineq);
controlPenaltyMPC_LA_VSUBADD3_12(params->lb7, controlPenaltyMPC_z06, controlPenaltyMPC_lbIdx06, controlPenaltyMPC_llb06, controlPenaltyMPC_slb06, controlPenaltyMPC_rilb06, &info->dgap, &info->res_ineq);
controlPenaltyMPC_LA_VSUBADD2_12(controlPenaltyMPC_z06, controlPenaltyMPC_ubIdx06, params->ub7, controlPenaltyMPC_lub06, controlPenaltyMPC_sub06, controlPenaltyMPC_riub06, &info->dgap, &info->res_ineq);
controlPenaltyMPC_LA_VSUBADD3_12(params->lb8, controlPenaltyMPC_z07, controlPenaltyMPC_lbIdx07, controlPenaltyMPC_llb07, controlPenaltyMPC_slb07, controlPenaltyMPC_rilb07, &info->dgap, &info->res_ineq);
controlPenaltyMPC_LA_VSUBADD2_12(controlPenaltyMPC_z07, controlPenaltyMPC_ubIdx07, params->ub8, controlPenaltyMPC_lub07, controlPenaltyMPC_sub07, controlPenaltyMPC_riub07, &info->dgap, &info->res_ineq);
controlPenaltyMPC_LA_VSUBADD3_12(params->lb9, controlPenaltyMPC_z08, controlPenaltyMPC_lbIdx08, controlPenaltyMPC_llb08, controlPenaltyMPC_slb08, controlPenaltyMPC_rilb08, &info->dgap, &info->res_ineq);
controlPenaltyMPC_LA_VSUBADD2_12(controlPenaltyMPC_z08, controlPenaltyMPC_ubIdx08, params->ub9, controlPenaltyMPC_lub08, controlPenaltyMPC_sub08, controlPenaltyMPC_riub08, &info->dgap, &info->res_ineq);
controlPenaltyMPC_LA_VSUBADD3_12(params->lb10, controlPenaltyMPC_z09, controlPenaltyMPC_lbIdx09, controlPenaltyMPC_llb09, controlPenaltyMPC_slb09, controlPenaltyMPC_rilb09, &info->dgap, &info->res_ineq);
controlPenaltyMPC_LA_VSUBADD2_12(controlPenaltyMPC_z09, controlPenaltyMPC_ubIdx09, params->ub10, controlPenaltyMPC_lub09, controlPenaltyMPC_sub09, controlPenaltyMPC_riub09, &info->dgap, &info->res_ineq);
controlPenaltyMPC_LA_VSUBADD3_12(params->lb11, controlPenaltyMPC_z10, controlPenaltyMPC_lbIdx10, controlPenaltyMPC_llb10, controlPenaltyMPC_slb10, controlPenaltyMPC_rilb10, &info->dgap, &info->res_ineq);
controlPenaltyMPC_LA_VSUBADD2_12(controlPenaltyMPC_z10, controlPenaltyMPC_ubIdx10, params->ub11, controlPenaltyMPC_lub10, controlPenaltyMPC_sub10, controlPenaltyMPC_riub10, &info->dgap, &info->res_ineq);
controlPenaltyMPC_LA_VSUBADD3_12(params->lb12, controlPenaltyMPC_z11, controlPenaltyMPC_lbIdx11, controlPenaltyMPC_llb11, controlPenaltyMPC_slb11, controlPenaltyMPC_rilb11, &info->dgap, &info->res_ineq);
controlPenaltyMPC_LA_VSUBADD2_12(controlPenaltyMPC_z11, controlPenaltyMPC_ubIdx11, params->ub12, controlPenaltyMPC_lub11, controlPenaltyMPC_sub11, controlPenaltyMPC_riub11, &info->dgap, &info->res_ineq);
controlPenaltyMPC_LA_VSUBADD3_12(params->lb13, controlPenaltyMPC_z12, controlPenaltyMPC_lbIdx12, controlPenaltyMPC_llb12, controlPenaltyMPC_slb12, controlPenaltyMPC_rilb12, &info->dgap, &info->res_ineq);
controlPenaltyMPC_LA_VSUBADD2_12(controlPenaltyMPC_z12, controlPenaltyMPC_ubIdx12, params->ub13, controlPenaltyMPC_lub12, controlPenaltyMPC_sub12, controlPenaltyMPC_riub12, &info->dgap, &info->res_ineq);
controlPenaltyMPC_LA_VSUBADD3_12(params->lb14, controlPenaltyMPC_z13, controlPenaltyMPC_lbIdx13, controlPenaltyMPC_llb13, controlPenaltyMPC_slb13, controlPenaltyMPC_rilb13, &info->dgap, &info->res_ineq);
controlPenaltyMPC_LA_VSUBADD2_12(controlPenaltyMPC_z13, controlPenaltyMPC_ubIdx13, params->ub14, controlPenaltyMPC_lub13, controlPenaltyMPC_sub13, controlPenaltyMPC_riub13, &info->dgap, &info->res_ineq);
controlPenaltyMPC_LA_VSUBADD3_12(params->lb15, controlPenaltyMPC_z14, controlPenaltyMPC_lbIdx14, controlPenaltyMPC_llb14, controlPenaltyMPC_slb14, controlPenaltyMPC_rilb14, &info->dgap, &info->res_ineq);
controlPenaltyMPC_LA_VSUBADD2_6(controlPenaltyMPC_z14, controlPenaltyMPC_ubIdx14, params->ub15, controlPenaltyMPC_lub14, controlPenaltyMPC_sub14, controlPenaltyMPC_riub14, &info->dgap, &info->res_ineq);
controlPenaltyMPC_LA_INEQ_B_GRAD_12_12_12(controlPenaltyMPC_lub00, controlPenaltyMPC_sub00, controlPenaltyMPC_riub00, controlPenaltyMPC_llb00, controlPenaltyMPC_slb00, controlPenaltyMPC_rilb00, controlPenaltyMPC_lbIdx00, controlPenaltyMPC_ubIdx00, controlPenaltyMPC_grad_ineq00, controlPenaltyMPC_lubbysub00, controlPenaltyMPC_llbbyslb00);
controlPenaltyMPC_LA_INEQ_B_GRAD_12_12_12(controlPenaltyMPC_lub01, controlPenaltyMPC_sub01, controlPenaltyMPC_riub01, controlPenaltyMPC_llb01, controlPenaltyMPC_slb01, controlPenaltyMPC_rilb01, controlPenaltyMPC_lbIdx01, controlPenaltyMPC_ubIdx01, controlPenaltyMPC_grad_ineq01, controlPenaltyMPC_lubbysub01, controlPenaltyMPC_llbbyslb01);
controlPenaltyMPC_LA_INEQ_B_GRAD_12_12_12(controlPenaltyMPC_lub02, controlPenaltyMPC_sub02, controlPenaltyMPC_riub02, controlPenaltyMPC_llb02, controlPenaltyMPC_slb02, controlPenaltyMPC_rilb02, controlPenaltyMPC_lbIdx02, controlPenaltyMPC_ubIdx02, controlPenaltyMPC_grad_ineq02, controlPenaltyMPC_lubbysub02, controlPenaltyMPC_llbbyslb02);
controlPenaltyMPC_LA_INEQ_B_GRAD_12_12_12(controlPenaltyMPC_lub03, controlPenaltyMPC_sub03, controlPenaltyMPC_riub03, controlPenaltyMPC_llb03, controlPenaltyMPC_slb03, controlPenaltyMPC_rilb03, controlPenaltyMPC_lbIdx03, controlPenaltyMPC_ubIdx03, controlPenaltyMPC_grad_ineq03, controlPenaltyMPC_lubbysub03, controlPenaltyMPC_llbbyslb03);
controlPenaltyMPC_LA_INEQ_B_GRAD_12_12_12(controlPenaltyMPC_lub04, controlPenaltyMPC_sub04, controlPenaltyMPC_riub04, controlPenaltyMPC_llb04, controlPenaltyMPC_slb04, controlPenaltyMPC_rilb04, controlPenaltyMPC_lbIdx04, controlPenaltyMPC_ubIdx04, controlPenaltyMPC_grad_ineq04, controlPenaltyMPC_lubbysub04, controlPenaltyMPC_llbbyslb04);
controlPenaltyMPC_LA_INEQ_B_GRAD_12_12_12(controlPenaltyMPC_lub05, controlPenaltyMPC_sub05, controlPenaltyMPC_riub05, controlPenaltyMPC_llb05, controlPenaltyMPC_slb05, controlPenaltyMPC_rilb05, controlPenaltyMPC_lbIdx05, controlPenaltyMPC_ubIdx05, controlPenaltyMPC_grad_ineq05, controlPenaltyMPC_lubbysub05, controlPenaltyMPC_llbbyslb05);
controlPenaltyMPC_LA_INEQ_B_GRAD_12_12_12(controlPenaltyMPC_lub06, controlPenaltyMPC_sub06, controlPenaltyMPC_riub06, controlPenaltyMPC_llb06, controlPenaltyMPC_slb06, controlPenaltyMPC_rilb06, controlPenaltyMPC_lbIdx06, controlPenaltyMPC_ubIdx06, controlPenaltyMPC_grad_ineq06, controlPenaltyMPC_lubbysub06, controlPenaltyMPC_llbbyslb06);
controlPenaltyMPC_LA_INEQ_B_GRAD_12_12_12(controlPenaltyMPC_lub07, controlPenaltyMPC_sub07, controlPenaltyMPC_riub07, controlPenaltyMPC_llb07, controlPenaltyMPC_slb07, controlPenaltyMPC_rilb07, controlPenaltyMPC_lbIdx07, controlPenaltyMPC_ubIdx07, controlPenaltyMPC_grad_ineq07, controlPenaltyMPC_lubbysub07, controlPenaltyMPC_llbbyslb07);
controlPenaltyMPC_LA_INEQ_B_GRAD_12_12_12(controlPenaltyMPC_lub08, controlPenaltyMPC_sub08, controlPenaltyMPC_riub08, controlPenaltyMPC_llb08, controlPenaltyMPC_slb08, controlPenaltyMPC_rilb08, controlPenaltyMPC_lbIdx08, controlPenaltyMPC_ubIdx08, controlPenaltyMPC_grad_ineq08, controlPenaltyMPC_lubbysub08, controlPenaltyMPC_llbbyslb08);
controlPenaltyMPC_LA_INEQ_B_GRAD_12_12_12(controlPenaltyMPC_lub09, controlPenaltyMPC_sub09, controlPenaltyMPC_riub09, controlPenaltyMPC_llb09, controlPenaltyMPC_slb09, controlPenaltyMPC_rilb09, controlPenaltyMPC_lbIdx09, controlPenaltyMPC_ubIdx09, controlPenaltyMPC_grad_ineq09, controlPenaltyMPC_lubbysub09, controlPenaltyMPC_llbbyslb09);
controlPenaltyMPC_LA_INEQ_B_GRAD_12_12_12(controlPenaltyMPC_lub10, controlPenaltyMPC_sub10, controlPenaltyMPC_riub10, controlPenaltyMPC_llb10, controlPenaltyMPC_slb10, controlPenaltyMPC_rilb10, controlPenaltyMPC_lbIdx10, controlPenaltyMPC_ubIdx10, controlPenaltyMPC_grad_ineq10, controlPenaltyMPC_lubbysub10, controlPenaltyMPC_llbbyslb10);
controlPenaltyMPC_LA_INEQ_B_GRAD_12_12_12(controlPenaltyMPC_lub11, controlPenaltyMPC_sub11, controlPenaltyMPC_riub11, controlPenaltyMPC_llb11, controlPenaltyMPC_slb11, controlPenaltyMPC_rilb11, controlPenaltyMPC_lbIdx11, controlPenaltyMPC_ubIdx11, controlPenaltyMPC_grad_ineq11, controlPenaltyMPC_lubbysub11, controlPenaltyMPC_llbbyslb11);
controlPenaltyMPC_LA_INEQ_B_GRAD_12_12_12(controlPenaltyMPC_lub12, controlPenaltyMPC_sub12, controlPenaltyMPC_riub12, controlPenaltyMPC_llb12, controlPenaltyMPC_slb12, controlPenaltyMPC_rilb12, controlPenaltyMPC_lbIdx12, controlPenaltyMPC_ubIdx12, controlPenaltyMPC_grad_ineq12, controlPenaltyMPC_lubbysub12, controlPenaltyMPC_llbbyslb12);
controlPenaltyMPC_LA_INEQ_B_GRAD_12_12_12(controlPenaltyMPC_lub13, controlPenaltyMPC_sub13, controlPenaltyMPC_riub13, controlPenaltyMPC_llb13, controlPenaltyMPC_slb13, controlPenaltyMPC_rilb13, controlPenaltyMPC_lbIdx13, controlPenaltyMPC_ubIdx13, controlPenaltyMPC_grad_ineq13, controlPenaltyMPC_lubbysub13, controlPenaltyMPC_llbbyslb13);
controlPenaltyMPC_LA_INEQ_B_GRAD_12_12_6(controlPenaltyMPC_lub14, controlPenaltyMPC_sub14, controlPenaltyMPC_riub14, controlPenaltyMPC_llb14, controlPenaltyMPC_slb14, controlPenaltyMPC_rilb14, controlPenaltyMPC_lbIdx14, controlPenaltyMPC_ubIdx14, controlPenaltyMPC_grad_ineq14, controlPenaltyMPC_lubbysub14, controlPenaltyMPC_llbbyslb14);
info->dobj = info->pobj - info->dgap;
info->rdgap = info->pobj ? info->dgap / info->pobj : 1e6;
if( info->rdgap < 0 ) info->rdgap = -info->rdgap;
if( info->mu < controlPenaltyMPC_SET_ACC_KKTCOMPL
    && (info->rdgap < controlPenaltyMPC_SET_ACC_RDGAP || info->dgap < controlPenaltyMPC_SET_ACC_KKTCOMPL)
    && info->res_eq < controlPenaltyMPC_SET_ACC_RESEQ
    && info->res_ineq < controlPenaltyMPC_SET_ACC_RESINEQ ){
exitcode = controlPenaltyMPC_OPTIMAL; break; }
if( info->it == controlPenaltyMPC_SET_MAXIT ){
exitcode = controlPenaltyMPC_MAXITREACHED; break; }
controlPenaltyMPC_LA_VVADD3_180(controlPenaltyMPC_grad_cost, controlPenaltyMPC_grad_eq, controlPenaltyMPC_grad_ineq, controlPenaltyMPC_rd);
controlPenaltyMPC_LA_DIAG_CHOL_ONELOOP_LBUB_12_12_12(params->Q1, controlPenaltyMPC_llbbyslb00, controlPenaltyMPC_lbIdx00, controlPenaltyMPC_lubbysub00, controlPenaltyMPC_ubIdx00, controlPenaltyMPC_Phi00);
controlPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_6_12(controlPenaltyMPC_Phi00, controlPenaltyMPC_C00, controlPenaltyMPC_V00);
controlPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_6_12(controlPenaltyMPC_Phi00, controlPenaltyMPC_D00, controlPenaltyMPC_W00);
controlPenaltyMPC_LA_DENSE_DIAGZERO_MMTM_6_12_6(controlPenaltyMPC_W00, controlPenaltyMPC_V00, controlPenaltyMPC_Ysd01);
controlPenaltyMPC_LA_DIAG_FORWARDSUB_12(controlPenaltyMPC_Phi00, controlPenaltyMPC_rd00, controlPenaltyMPC_Lbyrd00);
controlPenaltyMPC_LA_DIAG_CHOL_ONELOOP_LBUB_12_12_12(params->Q2, controlPenaltyMPC_llbbyslb01, controlPenaltyMPC_lbIdx01, controlPenaltyMPC_lubbysub01, controlPenaltyMPC_ubIdx01, controlPenaltyMPC_Phi01);
controlPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_6_12(controlPenaltyMPC_Phi01, controlPenaltyMPC_C00, controlPenaltyMPC_V01);
controlPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_6_12(controlPenaltyMPC_Phi01, controlPenaltyMPC_D01, controlPenaltyMPC_W01);
controlPenaltyMPC_LA_DENSE_DIAGZERO_MMTM_6_12_6(controlPenaltyMPC_W01, controlPenaltyMPC_V01, controlPenaltyMPC_Ysd02);
controlPenaltyMPC_LA_DIAG_FORWARDSUB_12(controlPenaltyMPC_Phi01, controlPenaltyMPC_rd01, controlPenaltyMPC_Lbyrd01);
controlPenaltyMPC_LA_DIAG_CHOL_ONELOOP_LBUB_12_12_12(params->Q3, controlPenaltyMPC_llbbyslb02, controlPenaltyMPC_lbIdx02, controlPenaltyMPC_lubbysub02, controlPenaltyMPC_ubIdx02, controlPenaltyMPC_Phi02);
controlPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_6_12(controlPenaltyMPC_Phi02, controlPenaltyMPC_C00, controlPenaltyMPC_V02);
controlPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_6_12(controlPenaltyMPC_Phi02, controlPenaltyMPC_D01, controlPenaltyMPC_W02);
controlPenaltyMPC_LA_DENSE_DIAGZERO_MMTM_6_12_6(controlPenaltyMPC_W02, controlPenaltyMPC_V02, controlPenaltyMPC_Ysd03);
controlPenaltyMPC_LA_DIAG_FORWARDSUB_12(controlPenaltyMPC_Phi02, controlPenaltyMPC_rd02, controlPenaltyMPC_Lbyrd02);
controlPenaltyMPC_LA_DIAG_CHOL_ONELOOP_LBUB_12_12_12(params->Q4, controlPenaltyMPC_llbbyslb03, controlPenaltyMPC_lbIdx03, controlPenaltyMPC_lubbysub03, controlPenaltyMPC_ubIdx03, controlPenaltyMPC_Phi03);
controlPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_6_12(controlPenaltyMPC_Phi03, controlPenaltyMPC_C00, controlPenaltyMPC_V03);
controlPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_6_12(controlPenaltyMPC_Phi03, controlPenaltyMPC_D01, controlPenaltyMPC_W03);
controlPenaltyMPC_LA_DENSE_DIAGZERO_MMTM_6_12_6(controlPenaltyMPC_W03, controlPenaltyMPC_V03, controlPenaltyMPC_Ysd04);
controlPenaltyMPC_LA_DIAG_FORWARDSUB_12(controlPenaltyMPC_Phi03, controlPenaltyMPC_rd03, controlPenaltyMPC_Lbyrd03);
controlPenaltyMPC_LA_DIAG_CHOL_ONELOOP_LBUB_12_12_12(params->Q5, controlPenaltyMPC_llbbyslb04, controlPenaltyMPC_lbIdx04, controlPenaltyMPC_lubbysub04, controlPenaltyMPC_ubIdx04, controlPenaltyMPC_Phi04);
controlPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_6_12(controlPenaltyMPC_Phi04, controlPenaltyMPC_C00, controlPenaltyMPC_V04);
controlPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_6_12(controlPenaltyMPC_Phi04, controlPenaltyMPC_D01, controlPenaltyMPC_W04);
controlPenaltyMPC_LA_DENSE_DIAGZERO_MMTM_6_12_6(controlPenaltyMPC_W04, controlPenaltyMPC_V04, controlPenaltyMPC_Ysd05);
controlPenaltyMPC_LA_DIAG_FORWARDSUB_12(controlPenaltyMPC_Phi04, controlPenaltyMPC_rd04, controlPenaltyMPC_Lbyrd04);
controlPenaltyMPC_LA_DIAG_CHOL_ONELOOP_LBUB_12_12_12(params->Q6, controlPenaltyMPC_llbbyslb05, controlPenaltyMPC_lbIdx05, controlPenaltyMPC_lubbysub05, controlPenaltyMPC_ubIdx05, controlPenaltyMPC_Phi05);
controlPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_6_12(controlPenaltyMPC_Phi05, controlPenaltyMPC_C00, controlPenaltyMPC_V05);
controlPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_6_12(controlPenaltyMPC_Phi05, controlPenaltyMPC_D01, controlPenaltyMPC_W05);
controlPenaltyMPC_LA_DENSE_DIAGZERO_MMTM_6_12_6(controlPenaltyMPC_W05, controlPenaltyMPC_V05, controlPenaltyMPC_Ysd06);
controlPenaltyMPC_LA_DIAG_FORWARDSUB_12(controlPenaltyMPC_Phi05, controlPenaltyMPC_rd05, controlPenaltyMPC_Lbyrd05);
controlPenaltyMPC_LA_DIAG_CHOL_ONELOOP_LBUB_12_12_12(params->Q7, controlPenaltyMPC_llbbyslb06, controlPenaltyMPC_lbIdx06, controlPenaltyMPC_lubbysub06, controlPenaltyMPC_ubIdx06, controlPenaltyMPC_Phi06);
controlPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_6_12(controlPenaltyMPC_Phi06, controlPenaltyMPC_C00, controlPenaltyMPC_V06);
controlPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_6_12(controlPenaltyMPC_Phi06, controlPenaltyMPC_D01, controlPenaltyMPC_W06);
controlPenaltyMPC_LA_DENSE_DIAGZERO_MMTM_6_12_6(controlPenaltyMPC_W06, controlPenaltyMPC_V06, controlPenaltyMPC_Ysd07);
controlPenaltyMPC_LA_DIAG_FORWARDSUB_12(controlPenaltyMPC_Phi06, controlPenaltyMPC_rd06, controlPenaltyMPC_Lbyrd06);
controlPenaltyMPC_LA_DIAG_CHOL_ONELOOP_LBUB_12_12_12(params->Q8, controlPenaltyMPC_llbbyslb07, controlPenaltyMPC_lbIdx07, controlPenaltyMPC_lubbysub07, controlPenaltyMPC_ubIdx07, controlPenaltyMPC_Phi07);
controlPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_6_12(controlPenaltyMPC_Phi07, controlPenaltyMPC_C00, controlPenaltyMPC_V07);
controlPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_6_12(controlPenaltyMPC_Phi07, controlPenaltyMPC_D01, controlPenaltyMPC_W07);
controlPenaltyMPC_LA_DENSE_DIAGZERO_MMTM_6_12_6(controlPenaltyMPC_W07, controlPenaltyMPC_V07, controlPenaltyMPC_Ysd08);
controlPenaltyMPC_LA_DIAG_FORWARDSUB_12(controlPenaltyMPC_Phi07, controlPenaltyMPC_rd07, controlPenaltyMPC_Lbyrd07);
controlPenaltyMPC_LA_DIAG_CHOL_ONELOOP_LBUB_12_12_12(params->Q9, controlPenaltyMPC_llbbyslb08, controlPenaltyMPC_lbIdx08, controlPenaltyMPC_lubbysub08, controlPenaltyMPC_ubIdx08, controlPenaltyMPC_Phi08);
controlPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_6_12(controlPenaltyMPC_Phi08, controlPenaltyMPC_C00, controlPenaltyMPC_V08);
controlPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_6_12(controlPenaltyMPC_Phi08, controlPenaltyMPC_D01, controlPenaltyMPC_W08);
controlPenaltyMPC_LA_DENSE_DIAGZERO_MMTM_6_12_6(controlPenaltyMPC_W08, controlPenaltyMPC_V08, controlPenaltyMPC_Ysd09);
controlPenaltyMPC_LA_DIAG_FORWARDSUB_12(controlPenaltyMPC_Phi08, controlPenaltyMPC_rd08, controlPenaltyMPC_Lbyrd08);
controlPenaltyMPC_LA_DIAG_CHOL_ONELOOP_LBUB_12_12_12(params->Q10, controlPenaltyMPC_llbbyslb09, controlPenaltyMPC_lbIdx09, controlPenaltyMPC_lubbysub09, controlPenaltyMPC_ubIdx09, controlPenaltyMPC_Phi09);
controlPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_6_12(controlPenaltyMPC_Phi09, controlPenaltyMPC_C00, controlPenaltyMPC_V09);
controlPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_6_12(controlPenaltyMPC_Phi09, controlPenaltyMPC_D01, controlPenaltyMPC_W09);
controlPenaltyMPC_LA_DENSE_DIAGZERO_MMTM_6_12_6(controlPenaltyMPC_W09, controlPenaltyMPC_V09, controlPenaltyMPC_Ysd10);
controlPenaltyMPC_LA_DIAG_FORWARDSUB_12(controlPenaltyMPC_Phi09, controlPenaltyMPC_rd09, controlPenaltyMPC_Lbyrd09);
controlPenaltyMPC_LA_DIAG_CHOL_ONELOOP_LBUB_12_12_12(params->Q11, controlPenaltyMPC_llbbyslb10, controlPenaltyMPC_lbIdx10, controlPenaltyMPC_lubbysub10, controlPenaltyMPC_ubIdx10, controlPenaltyMPC_Phi10);
controlPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_6_12(controlPenaltyMPC_Phi10, controlPenaltyMPC_C00, controlPenaltyMPC_V10);
controlPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_6_12(controlPenaltyMPC_Phi10, controlPenaltyMPC_D01, controlPenaltyMPC_W10);
controlPenaltyMPC_LA_DENSE_DIAGZERO_MMTM_6_12_6(controlPenaltyMPC_W10, controlPenaltyMPC_V10, controlPenaltyMPC_Ysd11);
controlPenaltyMPC_LA_DIAG_FORWARDSUB_12(controlPenaltyMPC_Phi10, controlPenaltyMPC_rd10, controlPenaltyMPC_Lbyrd10);
controlPenaltyMPC_LA_DIAG_CHOL_ONELOOP_LBUB_12_12_12(params->Q12, controlPenaltyMPC_llbbyslb11, controlPenaltyMPC_lbIdx11, controlPenaltyMPC_lubbysub11, controlPenaltyMPC_ubIdx11, controlPenaltyMPC_Phi11);
controlPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_6_12(controlPenaltyMPC_Phi11, controlPenaltyMPC_C00, controlPenaltyMPC_V11);
controlPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_6_12(controlPenaltyMPC_Phi11, controlPenaltyMPC_D01, controlPenaltyMPC_W11);
controlPenaltyMPC_LA_DENSE_DIAGZERO_MMTM_6_12_6(controlPenaltyMPC_W11, controlPenaltyMPC_V11, controlPenaltyMPC_Ysd12);
controlPenaltyMPC_LA_DIAG_FORWARDSUB_12(controlPenaltyMPC_Phi11, controlPenaltyMPC_rd11, controlPenaltyMPC_Lbyrd11);
controlPenaltyMPC_LA_DIAG_CHOL_ONELOOP_LBUB_12_12_12(params->Q13, controlPenaltyMPC_llbbyslb12, controlPenaltyMPC_lbIdx12, controlPenaltyMPC_lubbysub12, controlPenaltyMPC_ubIdx12, controlPenaltyMPC_Phi12);
controlPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_6_12(controlPenaltyMPC_Phi12, controlPenaltyMPC_C00, controlPenaltyMPC_V12);
controlPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_6_12(controlPenaltyMPC_Phi12, controlPenaltyMPC_D01, controlPenaltyMPC_W12);
controlPenaltyMPC_LA_DENSE_DIAGZERO_MMTM_6_12_6(controlPenaltyMPC_W12, controlPenaltyMPC_V12, controlPenaltyMPC_Ysd13);
controlPenaltyMPC_LA_DIAG_FORWARDSUB_12(controlPenaltyMPC_Phi12, controlPenaltyMPC_rd12, controlPenaltyMPC_Lbyrd12);
controlPenaltyMPC_LA_DIAG_CHOL_ONELOOP_LBUB_12_12_12(params->Q14, controlPenaltyMPC_llbbyslb13, controlPenaltyMPC_lbIdx13, controlPenaltyMPC_lubbysub13, controlPenaltyMPC_ubIdx13, controlPenaltyMPC_Phi13);
controlPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_6_12(controlPenaltyMPC_Phi13, controlPenaltyMPC_C00, controlPenaltyMPC_V13);
controlPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_6_12(controlPenaltyMPC_Phi13, controlPenaltyMPC_D01, controlPenaltyMPC_W13);
controlPenaltyMPC_LA_DENSE_DIAGZERO_MMTM_6_12_6(controlPenaltyMPC_W13, controlPenaltyMPC_V13, controlPenaltyMPC_Ysd14);
controlPenaltyMPC_LA_DIAG_FORWARDSUB_12(controlPenaltyMPC_Phi13, controlPenaltyMPC_rd13, controlPenaltyMPC_Lbyrd13);
controlPenaltyMPC_LA_DIAG_CHOL_LBUB_12_12_6(params->Q15, controlPenaltyMPC_llbbyslb14, controlPenaltyMPC_lbIdx14, controlPenaltyMPC_lubbysub14, controlPenaltyMPC_ubIdx14, controlPenaltyMPC_Phi14);
controlPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_6_12(controlPenaltyMPC_Phi14, controlPenaltyMPC_D01, controlPenaltyMPC_W14);
controlPenaltyMPC_LA_DIAG_FORWARDSUB_12(controlPenaltyMPC_Phi14, controlPenaltyMPC_rd14, controlPenaltyMPC_Lbyrd14);
controlPenaltyMPC_LA_DIAGZERO_MMT_6(controlPenaltyMPC_W00, controlPenaltyMPC_Yd00);
controlPenaltyMPC_LA_DIAGZERO_MVMSUB7_6(controlPenaltyMPC_W00, controlPenaltyMPC_Lbyrd00, controlPenaltyMPC_re00, controlPenaltyMPC_beta00);
controlPenaltyMPC_LA_DENSE_DIAGZERO_MMT2_6_12_12(controlPenaltyMPC_V00, controlPenaltyMPC_W01, controlPenaltyMPC_Yd01);
controlPenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_6_12_12(controlPenaltyMPC_V00, controlPenaltyMPC_Lbyrd00, controlPenaltyMPC_W01, controlPenaltyMPC_Lbyrd01, controlPenaltyMPC_re01, controlPenaltyMPC_beta01);
controlPenaltyMPC_LA_DENSE_DIAGZERO_MMT2_6_12_12(controlPenaltyMPC_V01, controlPenaltyMPC_W02, controlPenaltyMPC_Yd02);
controlPenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_6_12_12(controlPenaltyMPC_V01, controlPenaltyMPC_Lbyrd01, controlPenaltyMPC_W02, controlPenaltyMPC_Lbyrd02, controlPenaltyMPC_re02, controlPenaltyMPC_beta02);
controlPenaltyMPC_LA_DENSE_DIAGZERO_MMT2_6_12_12(controlPenaltyMPC_V02, controlPenaltyMPC_W03, controlPenaltyMPC_Yd03);
controlPenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_6_12_12(controlPenaltyMPC_V02, controlPenaltyMPC_Lbyrd02, controlPenaltyMPC_W03, controlPenaltyMPC_Lbyrd03, controlPenaltyMPC_re03, controlPenaltyMPC_beta03);
controlPenaltyMPC_LA_DENSE_DIAGZERO_MMT2_6_12_12(controlPenaltyMPC_V03, controlPenaltyMPC_W04, controlPenaltyMPC_Yd04);
controlPenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_6_12_12(controlPenaltyMPC_V03, controlPenaltyMPC_Lbyrd03, controlPenaltyMPC_W04, controlPenaltyMPC_Lbyrd04, controlPenaltyMPC_re04, controlPenaltyMPC_beta04);
controlPenaltyMPC_LA_DENSE_DIAGZERO_MMT2_6_12_12(controlPenaltyMPC_V04, controlPenaltyMPC_W05, controlPenaltyMPC_Yd05);
controlPenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_6_12_12(controlPenaltyMPC_V04, controlPenaltyMPC_Lbyrd04, controlPenaltyMPC_W05, controlPenaltyMPC_Lbyrd05, controlPenaltyMPC_re05, controlPenaltyMPC_beta05);
controlPenaltyMPC_LA_DENSE_DIAGZERO_MMT2_6_12_12(controlPenaltyMPC_V05, controlPenaltyMPC_W06, controlPenaltyMPC_Yd06);
controlPenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_6_12_12(controlPenaltyMPC_V05, controlPenaltyMPC_Lbyrd05, controlPenaltyMPC_W06, controlPenaltyMPC_Lbyrd06, controlPenaltyMPC_re06, controlPenaltyMPC_beta06);
controlPenaltyMPC_LA_DENSE_DIAGZERO_MMT2_6_12_12(controlPenaltyMPC_V06, controlPenaltyMPC_W07, controlPenaltyMPC_Yd07);
controlPenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_6_12_12(controlPenaltyMPC_V06, controlPenaltyMPC_Lbyrd06, controlPenaltyMPC_W07, controlPenaltyMPC_Lbyrd07, controlPenaltyMPC_re07, controlPenaltyMPC_beta07);
controlPenaltyMPC_LA_DENSE_DIAGZERO_MMT2_6_12_12(controlPenaltyMPC_V07, controlPenaltyMPC_W08, controlPenaltyMPC_Yd08);
controlPenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_6_12_12(controlPenaltyMPC_V07, controlPenaltyMPC_Lbyrd07, controlPenaltyMPC_W08, controlPenaltyMPC_Lbyrd08, controlPenaltyMPC_re08, controlPenaltyMPC_beta08);
controlPenaltyMPC_LA_DENSE_DIAGZERO_MMT2_6_12_12(controlPenaltyMPC_V08, controlPenaltyMPC_W09, controlPenaltyMPC_Yd09);
controlPenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_6_12_12(controlPenaltyMPC_V08, controlPenaltyMPC_Lbyrd08, controlPenaltyMPC_W09, controlPenaltyMPC_Lbyrd09, controlPenaltyMPC_re09, controlPenaltyMPC_beta09);
controlPenaltyMPC_LA_DENSE_DIAGZERO_MMT2_6_12_12(controlPenaltyMPC_V09, controlPenaltyMPC_W10, controlPenaltyMPC_Yd10);
controlPenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_6_12_12(controlPenaltyMPC_V09, controlPenaltyMPC_Lbyrd09, controlPenaltyMPC_W10, controlPenaltyMPC_Lbyrd10, controlPenaltyMPC_re10, controlPenaltyMPC_beta10);
controlPenaltyMPC_LA_DENSE_DIAGZERO_MMT2_6_12_12(controlPenaltyMPC_V10, controlPenaltyMPC_W11, controlPenaltyMPC_Yd11);
controlPenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_6_12_12(controlPenaltyMPC_V10, controlPenaltyMPC_Lbyrd10, controlPenaltyMPC_W11, controlPenaltyMPC_Lbyrd11, controlPenaltyMPC_re11, controlPenaltyMPC_beta11);
controlPenaltyMPC_LA_DENSE_DIAGZERO_MMT2_6_12_12(controlPenaltyMPC_V11, controlPenaltyMPC_W12, controlPenaltyMPC_Yd12);
controlPenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_6_12_12(controlPenaltyMPC_V11, controlPenaltyMPC_Lbyrd11, controlPenaltyMPC_W12, controlPenaltyMPC_Lbyrd12, controlPenaltyMPC_re12, controlPenaltyMPC_beta12);
controlPenaltyMPC_LA_DENSE_DIAGZERO_MMT2_6_12_12(controlPenaltyMPC_V12, controlPenaltyMPC_W13, controlPenaltyMPC_Yd13);
controlPenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_6_12_12(controlPenaltyMPC_V12, controlPenaltyMPC_Lbyrd12, controlPenaltyMPC_W13, controlPenaltyMPC_Lbyrd13, controlPenaltyMPC_re13, controlPenaltyMPC_beta13);
controlPenaltyMPC_LA_DENSE_DIAGZERO_MMT2_6_12_12(controlPenaltyMPC_V13, controlPenaltyMPC_W14, controlPenaltyMPC_Yd14);
controlPenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_6_12_12(controlPenaltyMPC_V13, controlPenaltyMPC_Lbyrd13, controlPenaltyMPC_W14, controlPenaltyMPC_Lbyrd14, controlPenaltyMPC_re14, controlPenaltyMPC_beta14);
controlPenaltyMPC_LA_DENSE_CHOL_6(controlPenaltyMPC_Yd00, controlPenaltyMPC_Ld00);
controlPenaltyMPC_LA_DENSE_FORWARDSUB_6(controlPenaltyMPC_Ld00, controlPenaltyMPC_beta00, controlPenaltyMPC_yy00);
controlPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_6_6(controlPenaltyMPC_Ld00, controlPenaltyMPC_Ysd01, controlPenaltyMPC_Lsd01);
controlPenaltyMPC_LA_DENSE_MMTSUB_6_6(controlPenaltyMPC_Lsd01, controlPenaltyMPC_Yd01);
controlPenaltyMPC_LA_DENSE_CHOL_6(controlPenaltyMPC_Yd01, controlPenaltyMPC_Ld01);
controlPenaltyMPC_LA_DENSE_MVMSUB1_6_6(controlPenaltyMPC_Lsd01, controlPenaltyMPC_yy00, controlPenaltyMPC_beta01, controlPenaltyMPC_bmy01);
controlPenaltyMPC_LA_DENSE_FORWARDSUB_6(controlPenaltyMPC_Ld01, controlPenaltyMPC_bmy01, controlPenaltyMPC_yy01);
controlPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_6_6(controlPenaltyMPC_Ld01, controlPenaltyMPC_Ysd02, controlPenaltyMPC_Lsd02);
controlPenaltyMPC_LA_DENSE_MMTSUB_6_6(controlPenaltyMPC_Lsd02, controlPenaltyMPC_Yd02);
controlPenaltyMPC_LA_DENSE_CHOL_6(controlPenaltyMPC_Yd02, controlPenaltyMPC_Ld02);
controlPenaltyMPC_LA_DENSE_MVMSUB1_6_6(controlPenaltyMPC_Lsd02, controlPenaltyMPC_yy01, controlPenaltyMPC_beta02, controlPenaltyMPC_bmy02);
controlPenaltyMPC_LA_DENSE_FORWARDSUB_6(controlPenaltyMPC_Ld02, controlPenaltyMPC_bmy02, controlPenaltyMPC_yy02);
controlPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_6_6(controlPenaltyMPC_Ld02, controlPenaltyMPC_Ysd03, controlPenaltyMPC_Lsd03);
controlPenaltyMPC_LA_DENSE_MMTSUB_6_6(controlPenaltyMPC_Lsd03, controlPenaltyMPC_Yd03);
controlPenaltyMPC_LA_DENSE_CHOL_6(controlPenaltyMPC_Yd03, controlPenaltyMPC_Ld03);
controlPenaltyMPC_LA_DENSE_MVMSUB1_6_6(controlPenaltyMPC_Lsd03, controlPenaltyMPC_yy02, controlPenaltyMPC_beta03, controlPenaltyMPC_bmy03);
controlPenaltyMPC_LA_DENSE_FORWARDSUB_6(controlPenaltyMPC_Ld03, controlPenaltyMPC_bmy03, controlPenaltyMPC_yy03);
controlPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_6_6(controlPenaltyMPC_Ld03, controlPenaltyMPC_Ysd04, controlPenaltyMPC_Lsd04);
controlPenaltyMPC_LA_DENSE_MMTSUB_6_6(controlPenaltyMPC_Lsd04, controlPenaltyMPC_Yd04);
controlPenaltyMPC_LA_DENSE_CHOL_6(controlPenaltyMPC_Yd04, controlPenaltyMPC_Ld04);
controlPenaltyMPC_LA_DENSE_MVMSUB1_6_6(controlPenaltyMPC_Lsd04, controlPenaltyMPC_yy03, controlPenaltyMPC_beta04, controlPenaltyMPC_bmy04);
controlPenaltyMPC_LA_DENSE_FORWARDSUB_6(controlPenaltyMPC_Ld04, controlPenaltyMPC_bmy04, controlPenaltyMPC_yy04);
controlPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_6_6(controlPenaltyMPC_Ld04, controlPenaltyMPC_Ysd05, controlPenaltyMPC_Lsd05);
controlPenaltyMPC_LA_DENSE_MMTSUB_6_6(controlPenaltyMPC_Lsd05, controlPenaltyMPC_Yd05);
controlPenaltyMPC_LA_DENSE_CHOL_6(controlPenaltyMPC_Yd05, controlPenaltyMPC_Ld05);
controlPenaltyMPC_LA_DENSE_MVMSUB1_6_6(controlPenaltyMPC_Lsd05, controlPenaltyMPC_yy04, controlPenaltyMPC_beta05, controlPenaltyMPC_bmy05);
controlPenaltyMPC_LA_DENSE_FORWARDSUB_6(controlPenaltyMPC_Ld05, controlPenaltyMPC_bmy05, controlPenaltyMPC_yy05);
controlPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_6_6(controlPenaltyMPC_Ld05, controlPenaltyMPC_Ysd06, controlPenaltyMPC_Lsd06);
controlPenaltyMPC_LA_DENSE_MMTSUB_6_6(controlPenaltyMPC_Lsd06, controlPenaltyMPC_Yd06);
controlPenaltyMPC_LA_DENSE_CHOL_6(controlPenaltyMPC_Yd06, controlPenaltyMPC_Ld06);
controlPenaltyMPC_LA_DENSE_MVMSUB1_6_6(controlPenaltyMPC_Lsd06, controlPenaltyMPC_yy05, controlPenaltyMPC_beta06, controlPenaltyMPC_bmy06);
controlPenaltyMPC_LA_DENSE_FORWARDSUB_6(controlPenaltyMPC_Ld06, controlPenaltyMPC_bmy06, controlPenaltyMPC_yy06);
controlPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_6_6(controlPenaltyMPC_Ld06, controlPenaltyMPC_Ysd07, controlPenaltyMPC_Lsd07);
controlPenaltyMPC_LA_DENSE_MMTSUB_6_6(controlPenaltyMPC_Lsd07, controlPenaltyMPC_Yd07);
controlPenaltyMPC_LA_DENSE_CHOL_6(controlPenaltyMPC_Yd07, controlPenaltyMPC_Ld07);
controlPenaltyMPC_LA_DENSE_MVMSUB1_6_6(controlPenaltyMPC_Lsd07, controlPenaltyMPC_yy06, controlPenaltyMPC_beta07, controlPenaltyMPC_bmy07);
controlPenaltyMPC_LA_DENSE_FORWARDSUB_6(controlPenaltyMPC_Ld07, controlPenaltyMPC_bmy07, controlPenaltyMPC_yy07);
controlPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_6_6(controlPenaltyMPC_Ld07, controlPenaltyMPC_Ysd08, controlPenaltyMPC_Lsd08);
controlPenaltyMPC_LA_DENSE_MMTSUB_6_6(controlPenaltyMPC_Lsd08, controlPenaltyMPC_Yd08);
controlPenaltyMPC_LA_DENSE_CHOL_6(controlPenaltyMPC_Yd08, controlPenaltyMPC_Ld08);
controlPenaltyMPC_LA_DENSE_MVMSUB1_6_6(controlPenaltyMPC_Lsd08, controlPenaltyMPC_yy07, controlPenaltyMPC_beta08, controlPenaltyMPC_bmy08);
controlPenaltyMPC_LA_DENSE_FORWARDSUB_6(controlPenaltyMPC_Ld08, controlPenaltyMPC_bmy08, controlPenaltyMPC_yy08);
controlPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_6_6(controlPenaltyMPC_Ld08, controlPenaltyMPC_Ysd09, controlPenaltyMPC_Lsd09);
controlPenaltyMPC_LA_DENSE_MMTSUB_6_6(controlPenaltyMPC_Lsd09, controlPenaltyMPC_Yd09);
controlPenaltyMPC_LA_DENSE_CHOL_6(controlPenaltyMPC_Yd09, controlPenaltyMPC_Ld09);
controlPenaltyMPC_LA_DENSE_MVMSUB1_6_6(controlPenaltyMPC_Lsd09, controlPenaltyMPC_yy08, controlPenaltyMPC_beta09, controlPenaltyMPC_bmy09);
controlPenaltyMPC_LA_DENSE_FORWARDSUB_6(controlPenaltyMPC_Ld09, controlPenaltyMPC_bmy09, controlPenaltyMPC_yy09);
controlPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_6_6(controlPenaltyMPC_Ld09, controlPenaltyMPC_Ysd10, controlPenaltyMPC_Lsd10);
controlPenaltyMPC_LA_DENSE_MMTSUB_6_6(controlPenaltyMPC_Lsd10, controlPenaltyMPC_Yd10);
controlPenaltyMPC_LA_DENSE_CHOL_6(controlPenaltyMPC_Yd10, controlPenaltyMPC_Ld10);
controlPenaltyMPC_LA_DENSE_MVMSUB1_6_6(controlPenaltyMPC_Lsd10, controlPenaltyMPC_yy09, controlPenaltyMPC_beta10, controlPenaltyMPC_bmy10);
controlPenaltyMPC_LA_DENSE_FORWARDSUB_6(controlPenaltyMPC_Ld10, controlPenaltyMPC_bmy10, controlPenaltyMPC_yy10);
controlPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_6_6(controlPenaltyMPC_Ld10, controlPenaltyMPC_Ysd11, controlPenaltyMPC_Lsd11);
controlPenaltyMPC_LA_DENSE_MMTSUB_6_6(controlPenaltyMPC_Lsd11, controlPenaltyMPC_Yd11);
controlPenaltyMPC_LA_DENSE_CHOL_6(controlPenaltyMPC_Yd11, controlPenaltyMPC_Ld11);
controlPenaltyMPC_LA_DENSE_MVMSUB1_6_6(controlPenaltyMPC_Lsd11, controlPenaltyMPC_yy10, controlPenaltyMPC_beta11, controlPenaltyMPC_bmy11);
controlPenaltyMPC_LA_DENSE_FORWARDSUB_6(controlPenaltyMPC_Ld11, controlPenaltyMPC_bmy11, controlPenaltyMPC_yy11);
controlPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_6_6(controlPenaltyMPC_Ld11, controlPenaltyMPC_Ysd12, controlPenaltyMPC_Lsd12);
controlPenaltyMPC_LA_DENSE_MMTSUB_6_6(controlPenaltyMPC_Lsd12, controlPenaltyMPC_Yd12);
controlPenaltyMPC_LA_DENSE_CHOL_6(controlPenaltyMPC_Yd12, controlPenaltyMPC_Ld12);
controlPenaltyMPC_LA_DENSE_MVMSUB1_6_6(controlPenaltyMPC_Lsd12, controlPenaltyMPC_yy11, controlPenaltyMPC_beta12, controlPenaltyMPC_bmy12);
controlPenaltyMPC_LA_DENSE_FORWARDSUB_6(controlPenaltyMPC_Ld12, controlPenaltyMPC_bmy12, controlPenaltyMPC_yy12);
controlPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_6_6(controlPenaltyMPC_Ld12, controlPenaltyMPC_Ysd13, controlPenaltyMPC_Lsd13);
controlPenaltyMPC_LA_DENSE_MMTSUB_6_6(controlPenaltyMPC_Lsd13, controlPenaltyMPC_Yd13);
controlPenaltyMPC_LA_DENSE_CHOL_6(controlPenaltyMPC_Yd13, controlPenaltyMPC_Ld13);
controlPenaltyMPC_LA_DENSE_MVMSUB1_6_6(controlPenaltyMPC_Lsd13, controlPenaltyMPC_yy12, controlPenaltyMPC_beta13, controlPenaltyMPC_bmy13);
controlPenaltyMPC_LA_DENSE_FORWARDSUB_6(controlPenaltyMPC_Ld13, controlPenaltyMPC_bmy13, controlPenaltyMPC_yy13);
controlPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_6_6(controlPenaltyMPC_Ld13, controlPenaltyMPC_Ysd14, controlPenaltyMPC_Lsd14);
controlPenaltyMPC_LA_DENSE_MMTSUB_6_6(controlPenaltyMPC_Lsd14, controlPenaltyMPC_Yd14);
controlPenaltyMPC_LA_DENSE_CHOL_6(controlPenaltyMPC_Yd14, controlPenaltyMPC_Ld14);
controlPenaltyMPC_LA_DENSE_MVMSUB1_6_6(controlPenaltyMPC_Lsd14, controlPenaltyMPC_yy13, controlPenaltyMPC_beta14, controlPenaltyMPC_bmy14);
controlPenaltyMPC_LA_DENSE_FORWARDSUB_6(controlPenaltyMPC_Ld14, controlPenaltyMPC_bmy14, controlPenaltyMPC_yy14);
controlPenaltyMPC_LA_DENSE_BACKWARDSUB_6(controlPenaltyMPC_Ld14, controlPenaltyMPC_yy14, controlPenaltyMPC_dvaff14);
controlPenaltyMPC_LA_DENSE_MTVMSUB_6_6(controlPenaltyMPC_Lsd14, controlPenaltyMPC_dvaff14, controlPenaltyMPC_yy13, controlPenaltyMPC_bmy13);
controlPenaltyMPC_LA_DENSE_BACKWARDSUB_6(controlPenaltyMPC_Ld13, controlPenaltyMPC_bmy13, controlPenaltyMPC_dvaff13);
controlPenaltyMPC_LA_DENSE_MTVMSUB_6_6(controlPenaltyMPC_Lsd13, controlPenaltyMPC_dvaff13, controlPenaltyMPC_yy12, controlPenaltyMPC_bmy12);
controlPenaltyMPC_LA_DENSE_BACKWARDSUB_6(controlPenaltyMPC_Ld12, controlPenaltyMPC_bmy12, controlPenaltyMPC_dvaff12);
controlPenaltyMPC_LA_DENSE_MTVMSUB_6_6(controlPenaltyMPC_Lsd12, controlPenaltyMPC_dvaff12, controlPenaltyMPC_yy11, controlPenaltyMPC_bmy11);
controlPenaltyMPC_LA_DENSE_BACKWARDSUB_6(controlPenaltyMPC_Ld11, controlPenaltyMPC_bmy11, controlPenaltyMPC_dvaff11);
controlPenaltyMPC_LA_DENSE_MTVMSUB_6_6(controlPenaltyMPC_Lsd11, controlPenaltyMPC_dvaff11, controlPenaltyMPC_yy10, controlPenaltyMPC_bmy10);
controlPenaltyMPC_LA_DENSE_BACKWARDSUB_6(controlPenaltyMPC_Ld10, controlPenaltyMPC_bmy10, controlPenaltyMPC_dvaff10);
controlPenaltyMPC_LA_DENSE_MTVMSUB_6_6(controlPenaltyMPC_Lsd10, controlPenaltyMPC_dvaff10, controlPenaltyMPC_yy09, controlPenaltyMPC_bmy09);
controlPenaltyMPC_LA_DENSE_BACKWARDSUB_6(controlPenaltyMPC_Ld09, controlPenaltyMPC_bmy09, controlPenaltyMPC_dvaff09);
controlPenaltyMPC_LA_DENSE_MTVMSUB_6_6(controlPenaltyMPC_Lsd09, controlPenaltyMPC_dvaff09, controlPenaltyMPC_yy08, controlPenaltyMPC_bmy08);
controlPenaltyMPC_LA_DENSE_BACKWARDSUB_6(controlPenaltyMPC_Ld08, controlPenaltyMPC_bmy08, controlPenaltyMPC_dvaff08);
controlPenaltyMPC_LA_DENSE_MTVMSUB_6_6(controlPenaltyMPC_Lsd08, controlPenaltyMPC_dvaff08, controlPenaltyMPC_yy07, controlPenaltyMPC_bmy07);
controlPenaltyMPC_LA_DENSE_BACKWARDSUB_6(controlPenaltyMPC_Ld07, controlPenaltyMPC_bmy07, controlPenaltyMPC_dvaff07);
controlPenaltyMPC_LA_DENSE_MTVMSUB_6_6(controlPenaltyMPC_Lsd07, controlPenaltyMPC_dvaff07, controlPenaltyMPC_yy06, controlPenaltyMPC_bmy06);
controlPenaltyMPC_LA_DENSE_BACKWARDSUB_6(controlPenaltyMPC_Ld06, controlPenaltyMPC_bmy06, controlPenaltyMPC_dvaff06);
controlPenaltyMPC_LA_DENSE_MTVMSUB_6_6(controlPenaltyMPC_Lsd06, controlPenaltyMPC_dvaff06, controlPenaltyMPC_yy05, controlPenaltyMPC_bmy05);
controlPenaltyMPC_LA_DENSE_BACKWARDSUB_6(controlPenaltyMPC_Ld05, controlPenaltyMPC_bmy05, controlPenaltyMPC_dvaff05);
controlPenaltyMPC_LA_DENSE_MTVMSUB_6_6(controlPenaltyMPC_Lsd05, controlPenaltyMPC_dvaff05, controlPenaltyMPC_yy04, controlPenaltyMPC_bmy04);
controlPenaltyMPC_LA_DENSE_BACKWARDSUB_6(controlPenaltyMPC_Ld04, controlPenaltyMPC_bmy04, controlPenaltyMPC_dvaff04);
controlPenaltyMPC_LA_DENSE_MTVMSUB_6_6(controlPenaltyMPC_Lsd04, controlPenaltyMPC_dvaff04, controlPenaltyMPC_yy03, controlPenaltyMPC_bmy03);
controlPenaltyMPC_LA_DENSE_BACKWARDSUB_6(controlPenaltyMPC_Ld03, controlPenaltyMPC_bmy03, controlPenaltyMPC_dvaff03);
controlPenaltyMPC_LA_DENSE_MTVMSUB_6_6(controlPenaltyMPC_Lsd03, controlPenaltyMPC_dvaff03, controlPenaltyMPC_yy02, controlPenaltyMPC_bmy02);
controlPenaltyMPC_LA_DENSE_BACKWARDSUB_6(controlPenaltyMPC_Ld02, controlPenaltyMPC_bmy02, controlPenaltyMPC_dvaff02);
controlPenaltyMPC_LA_DENSE_MTVMSUB_6_6(controlPenaltyMPC_Lsd02, controlPenaltyMPC_dvaff02, controlPenaltyMPC_yy01, controlPenaltyMPC_bmy01);
controlPenaltyMPC_LA_DENSE_BACKWARDSUB_6(controlPenaltyMPC_Ld01, controlPenaltyMPC_bmy01, controlPenaltyMPC_dvaff01);
controlPenaltyMPC_LA_DENSE_MTVMSUB_6_6(controlPenaltyMPC_Lsd01, controlPenaltyMPC_dvaff01, controlPenaltyMPC_yy00, controlPenaltyMPC_bmy00);
controlPenaltyMPC_LA_DENSE_BACKWARDSUB_6(controlPenaltyMPC_Ld00, controlPenaltyMPC_bmy00, controlPenaltyMPC_dvaff00);
controlPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_6_12_6(controlPenaltyMPC_C00, controlPenaltyMPC_dvaff01, controlPenaltyMPC_D00, controlPenaltyMPC_dvaff00, controlPenaltyMPC_grad_eq00);
controlPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_6_12_6(controlPenaltyMPC_C00, controlPenaltyMPC_dvaff02, controlPenaltyMPC_D01, controlPenaltyMPC_dvaff01, controlPenaltyMPC_grad_eq01);
controlPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_6_12_6(controlPenaltyMPC_C00, controlPenaltyMPC_dvaff03, controlPenaltyMPC_D01, controlPenaltyMPC_dvaff02, controlPenaltyMPC_grad_eq02);
controlPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_6_12_6(controlPenaltyMPC_C00, controlPenaltyMPC_dvaff04, controlPenaltyMPC_D01, controlPenaltyMPC_dvaff03, controlPenaltyMPC_grad_eq03);
controlPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_6_12_6(controlPenaltyMPC_C00, controlPenaltyMPC_dvaff05, controlPenaltyMPC_D01, controlPenaltyMPC_dvaff04, controlPenaltyMPC_grad_eq04);
controlPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_6_12_6(controlPenaltyMPC_C00, controlPenaltyMPC_dvaff06, controlPenaltyMPC_D01, controlPenaltyMPC_dvaff05, controlPenaltyMPC_grad_eq05);
controlPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_6_12_6(controlPenaltyMPC_C00, controlPenaltyMPC_dvaff07, controlPenaltyMPC_D01, controlPenaltyMPC_dvaff06, controlPenaltyMPC_grad_eq06);
controlPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_6_12_6(controlPenaltyMPC_C00, controlPenaltyMPC_dvaff08, controlPenaltyMPC_D01, controlPenaltyMPC_dvaff07, controlPenaltyMPC_grad_eq07);
controlPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_6_12_6(controlPenaltyMPC_C00, controlPenaltyMPC_dvaff09, controlPenaltyMPC_D01, controlPenaltyMPC_dvaff08, controlPenaltyMPC_grad_eq08);
controlPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_6_12_6(controlPenaltyMPC_C00, controlPenaltyMPC_dvaff10, controlPenaltyMPC_D01, controlPenaltyMPC_dvaff09, controlPenaltyMPC_grad_eq09);
controlPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_6_12_6(controlPenaltyMPC_C00, controlPenaltyMPC_dvaff11, controlPenaltyMPC_D01, controlPenaltyMPC_dvaff10, controlPenaltyMPC_grad_eq10);
controlPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_6_12_6(controlPenaltyMPC_C00, controlPenaltyMPC_dvaff12, controlPenaltyMPC_D01, controlPenaltyMPC_dvaff11, controlPenaltyMPC_grad_eq11);
controlPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_6_12_6(controlPenaltyMPC_C00, controlPenaltyMPC_dvaff13, controlPenaltyMPC_D01, controlPenaltyMPC_dvaff12, controlPenaltyMPC_grad_eq12);
controlPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_6_12_6(controlPenaltyMPC_C00, controlPenaltyMPC_dvaff14, controlPenaltyMPC_D01, controlPenaltyMPC_dvaff13, controlPenaltyMPC_grad_eq13);
controlPenaltyMPC_LA_DIAGZERO_MTVM_6_12(controlPenaltyMPC_D01, controlPenaltyMPC_dvaff14, controlPenaltyMPC_grad_eq14);
controlPenaltyMPC_LA_VSUB2_180(controlPenaltyMPC_rd, controlPenaltyMPC_grad_eq, controlPenaltyMPC_rd);
controlPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_12(controlPenaltyMPC_Phi00, controlPenaltyMPC_rd00, controlPenaltyMPC_dzaff00);
controlPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_12(controlPenaltyMPC_Phi01, controlPenaltyMPC_rd01, controlPenaltyMPC_dzaff01);
controlPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_12(controlPenaltyMPC_Phi02, controlPenaltyMPC_rd02, controlPenaltyMPC_dzaff02);
controlPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_12(controlPenaltyMPC_Phi03, controlPenaltyMPC_rd03, controlPenaltyMPC_dzaff03);
controlPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_12(controlPenaltyMPC_Phi04, controlPenaltyMPC_rd04, controlPenaltyMPC_dzaff04);
controlPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_12(controlPenaltyMPC_Phi05, controlPenaltyMPC_rd05, controlPenaltyMPC_dzaff05);
controlPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_12(controlPenaltyMPC_Phi06, controlPenaltyMPC_rd06, controlPenaltyMPC_dzaff06);
controlPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_12(controlPenaltyMPC_Phi07, controlPenaltyMPC_rd07, controlPenaltyMPC_dzaff07);
controlPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_12(controlPenaltyMPC_Phi08, controlPenaltyMPC_rd08, controlPenaltyMPC_dzaff08);
controlPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_12(controlPenaltyMPC_Phi09, controlPenaltyMPC_rd09, controlPenaltyMPC_dzaff09);
controlPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_12(controlPenaltyMPC_Phi10, controlPenaltyMPC_rd10, controlPenaltyMPC_dzaff10);
controlPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_12(controlPenaltyMPC_Phi11, controlPenaltyMPC_rd11, controlPenaltyMPC_dzaff11);
controlPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_12(controlPenaltyMPC_Phi12, controlPenaltyMPC_rd12, controlPenaltyMPC_dzaff12);
controlPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_12(controlPenaltyMPC_Phi13, controlPenaltyMPC_rd13, controlPenaltyMPC_dzaff13);
controlPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_12(controlPenaltyMPC_Phi14, controlPenaltyMPC_rd14, controlPenaltyMPC_dzaff14);
controlPenaltyMPC_LA_VSUB_INDEXED_12(controlPenaltyMPC_dzaff00, controlPenaltyMPC_lbIdx00, controlPenaltyMPC_rilb00, controlPenaltyMPC_dslbaff00);
controlPenaltyMPC_LA_VSUB3_12(controlPenaltyMPC_llbbyslb00, controlPenaltyMPC_dslbaff00, controlPenaltyMPC_llb00, controlPenaltyMPC_dllbaff00);
controlPenaltyMPC_LA_VSUB2_INDEXED_12(controlPenaltyMPC_riub00, controlPenaltyMPC_dzaff00, controlPenaltyMPC_ubIdx00, controlPenaltyMPC_dsubaff00);
controlPenaltyMPC_LA_VSUB3_12(controlPenaltyMPC_lubbysub00, controlPenaltyMPC_dsubaff00, controlPenaltyMPC_lub00, controlPenaltyMPC_dlubaff00);
controlPenaltyMPC_LA_VSUB_INDEXED_12(controlPenaltyMPC_dzaff01, controlPenaltyMPC_lbIdx01, controlPenaltyMPC_rilb01, controlPenaltyMPC_dslbaff01);
controlPenaltyMPC_LA_VSUB3_12(controlPenaltyMPC_llbbyslb01, controlPenaltyMPC_dslbaff01, controlPenaltyMPC_llb01, controlPenaltyMPC_dllbaff01);
controlPenaltyMPC_LA_VSUB2_INDEXED_12(controlPenaltyMPC_riub01, controlPenaltyMPC_dzaff01, controlPenaltyMPC_ubIdx01, controlPenaltyMPC_dsubaff01);
controlPenaltyMPC_LA_VSUB3_12(controlPenaltyMPC_lubbysub01, controlPenaltyMPC_dsubaff01, controlPenaltyMPC_lub01, controlPenaltyMPC_dlubaff01);
controlPenaltyMPC_LA_VSUB_INDEXED_12(controlPenaltyMPC_dzaff02, controlPenaltyMPC_lbIdx02, controlPenaltyMPC_rilb02, controlPenaltyMPC_dslbaff02);
controlPenaltyMPC_LA_VSUB3_12(controlPenaltyMPC_llbbyslb02, controlPenaltyMPC_dslbaff02, controlPenaltyMPC_llb02, controlPenaltyMPC_dllbaff02);
controlPenaltyMPC_LA_VSUB2_INDEXED_12(controlPenaltyMPC_riub02, controlPenaltyMPC_dzaff02, controlPenaltyMPC_ubIdx02, controlPenaltyMPC_dsubaff02);
controlPenaltyMPC_LA_VSUB3_12(controlPenaltyMPC_lubbysub02, controlPenaltyMPC_dsubaff02, controlPenaltyMPC_lub02, controlPenaltyMPC_dlubaff02);
controlPenaltyMPC_LA_VSUB_INDEXED_12(controlPenaltyMPC_dzaff03, controlPenaltyMPC_lbIdx03, controlPenaltyMPC_rilb03, controlPenaltyMPC_dslbaff03);
controlPenaltyMPC_LA_VSUB3_12(controlPenaltyMPC_llbbyslb03, controlPenaltyMPC_dslbaff03, controlPenaltyMPC_llb03, controlPenaltyMPC_dllbaff03);
controlPenaltyMPC_LA_VSUB2_INDEXED_12(controlPenaltyMPC_riub03, controlPenaltyMPC_dzaff03, controlPenaltyMPC_ubIdx03, controlPenaltyMPC_dsubaff03);
controlPenaltyMPC_LA_VSUB3_12(controlPenaltyMPC_lubbysub03, controlPenaltyMPC_dsubaff03, controlPenaltyMPC_lub03, controlPenaltyMPC_dlubaff03);
controlPenaltyMPC_LA_VSUB_INDEXED_12(controlPenaltyMPC_dzaff04, controlPenaltyMPC_lbIdx04, controlPenaltyMPC_rilb04, controlPenaltyMPC_dslbaff04);
controlPenaltyMPC_LA_VSUB3_12(controlPenaltyMPC_llbbyslb04, controlPenaltyMPC_dslbaff04, controlPenaltyMPC_llb04, controlPenaltyMPC_dllbaff04);
controlPenaltyMPC_LA_VSUB2_INDEXED_12(controlPenaltyMPC_riub04, controlPenaltyMPC_dzaff04, controlPenaltyMPC_ubIdx04, controlPenaltyMPC_dsubaff04);
controlPenaltyMPC_LA_VSUB3_12(controlPenaltyMPC_lubbysub04, controlPenaltyMPC_dsubaff04, controlPenaltyMPC_lub04, controlPenaltyMPC_dlubaff04);
controlPenaltyMPC_LA_VSUB_INDEXED_12(controlPenaltyMPC_dzaff05, controlPenaltyMPC_lbIdx05, controlPenaltyMPC_rilb05, controlPenaltyMPC_dslbaff05);
controlPenaltyMPC_LA_VSUB3_12(controlPenaltyMPC_llbbyslb05, controlPenaltyMPC_dslbaff05, controlPenaltyMPC_llb05, controlPenaltyMPC_dllbaff05);
controlPenaltyMPC_LA_VSUB2_INDEXED_12(controlPenaltyMPC_riub05, controlPenaltyMPC_dzaff05, controlPenaltyMPC_ubIdx05, controlPenaltyMPC_dsubaff05);
controlPenaltyMPC_LA_VSUB3_12(controlPenaltyMPC_lubbysub05, controlPenaltyMPC_dsubaff05, controlPenaltyMPC_lub05, controlPenaltyMPC_dlubaff05);
controlPenaltyMPC_LA_VSUB_INDEXED_12(controlPenaltyMPC_dzaff06, controlPenaltyMPC_lbIdx06, controlPenaltyMPC_rilb06, controlPenaltyMPC_dslbaff06);
controlPenaltyMPC_LA_VSUB3_12(controlPenaltyMPC_llbbyslb06, controlPenaltyMPC_dslbaff06, controlPenaltyMPC_llb06, controlPenaltyMPC_dllbaff06);
controlPenaltyMPC_LA_VSUB2_INDEXED_12(controlPenaltyMPC_riub06, controlPenaltyMPC_dzaff06, controlPenaltyMPC_ubIdx06, controlPenaltyMPC_dsubaff06);
controlPenaltyMPC_LA_VSUB3_12(controlPenaltyMPC_lubbysub06, controlPenaltyMPC_dsubaff06, controlPenaltyMPC_lub06, controlPenaltyMPC_dlubaff06);
controlPenaltyMPC_LA_VSUB_INDEXED_12(controlPenaltyMPC_dzaff07, controlPenaltyMPC_lbIdx07, controlPenaltyMPC_rilb07, controlPenaltyMPC_dslbaff07);
controlPenaltyMPC_LA_VSUB3_12(controlPenaltyMPC_llbbyslb07, controlPenaltyMPC_dslbaff07, controlPenaltyMPC_llb07, controlPenaltyMPC_dllbaff07);
controlPenaltyMPC_LA_VSUB2_INDEXED_12(controlPenaltyMPC_riub07, controlPenaltyMPC_dzaff07, controlPenaltyMPC_ubIdx07, controlPenaltyMPC_dsubaff07);
controlPenaltyMPC_LA_VSUB3_12(controlPenaltyMPC_lubbysub07, controlPenaltyMPC_dsubaff07, controlPenaltyMPC_lub07, controlPenaltyMPC_dlubaff07);
controlPenaltyMPC_LA_VSUB_INDEXED_12(controlPenaltyMPC_dzaff08, controlPenaltyMPC_lbIdx08, controlPenaltyMPC_rilb08, controlPenaltyMPC_dslbaff08);
controlPenaltyMPC_LA_VSUB3_12(controlPenaltyMPC_llbbyslb08, controlPenaltyMPC_dslbaff08, controlPenaltyMPC_llb08, controlPenaltyMPC_dllbaff08);
controlPenaltyMPC_LA_VSUB2_INDEXED_12(controlPenaltyMPC_riub08, controlPenaltyMPC_dzaff08, controlPenaltyMPC_ubIdx08, controlPenaltyMPC_dsubaff08);
controlPenaltyMPC_LA_VSUB3_12(controlPenaltyMPC_lubbysub08, controlPenaltyMPC_dsubaff08, controlPenaltyMPC_lub08, controlPenaltyMPC_dlubaff08);
controlPenaltyMPC_LA_VSUB_INDEXED_12(controlPenaltyMPC_dzaff09, controlPenaltyMPC_lbIdx09, controlPenaltyMPC_rilb09, controlPenaltyMPC_dslbaff09);
controlPenaltyMPC_LA_VSUB3_12(controlPenaltyMPC_llbbyslb09, controlPenaltyMPC_dslbaff09, controlPenaltyMPC_llb09, controlPenaltyMPC_dllbaff09);
controlPenaltyMPC_LA_VSUB2_INDEXED_12(controlPenaltyMPC_riub09, controlPenaltyMPC_dzaff09, controlPenaltyMPC_ubIdx09, controlPenaltyMPC_dsubaff09);
controlPenaltyMPC_LA_VSUB3_12(controlPenaltyMPC_lubbysub09, controlPenaltyMPC_dsubaff09, controlPenaltyMPC_lub09, controlPenaltyMPC_dlubaff09);
controlPenaltyMPC_LA_VSUB_INDEXED_12(controlPenaltyMPC_dzaff10, controlPenaltyMPC_lbIdx10, controlPenaltyMPC_rilb10, controlPenaltyMPC_dslbaff10);
controlPenaltyMPC_LA_VSUB3_12(controlPenaltyMPC_llbbyslb10, controlPenaltyMPC_dslbaff10, controlPenaltyMPC_llb10, controlPenaltyMPC_dllbaff10);
controlPenaltyMPC_LA_VSUB2_INDEXED_12(controlPenaltyMPC_riub10, controlPenaltyMPC_dzaff10, controlPenaltyMPC_ubIdx10, controlPenaltyMPC_dsubaff10);
controlPenaltyMPC_LA_VSUB3_12(controlPenaltyMPC_lubbysub10, controlPenaltyMPC_dsubaff10, controlPenaltyMPC_lub10, controlPenaltyMPC_dlubaff10);
controlPenaltyMPC_LA_VSUB_INDEXED_12(controlPenaltyMPC_dzaff11, controlPenaltyMPC_lbIdx11, controlPenaltyMPC_rilb11, controlPenaltyMPC_dslbaff11);
controlPenaltyMPC_LA_VSUB3_12(controlPenaltyMPC_llbbyslb11, controlPenaltyMPC_dslbaff11, controlPenaltyMPC_llb11, controlPenaltyMPC_dllbaff11);
controlPenaltyMPC_LA_VSUB2_INDEXED_12(controlPenaltyMPC_riub11, controlPenaltyMPC_dzaff11, controlPenaltyMPC_ubIdx11, controlPenaltyMPC_dsubaff11);
controlPenaltyMPC_LA_VSUB3_12(controlPenaltyMPC_lubbysub11, controlPenaltyMPC_dsubaff11, controlPenaltyMPC_lub11, controlPenaltyMPC_dlubaff11);
controlPenaltyMPC_LA_VSUB_INDEXED_12(controlPenaltyMPC_dzaff12, controlPenaltyMPC_lbIdx12, controlPenaltyMPC_rilb12, controlPenaltyMPC_dslbaff12);
controlPenaltyMPC_LA_VSUB3_12(controlPenaltyMPC_llbbyslb12, controlPenaltyMPC_dslbaff12, controlPenaltyMPC_llb12, controlPenaltyMPC_dllbaff12);
controlPenaltyMPC_LA_VSUB2_INDEXED_12(controlPenaltyMPC_riub12, controlPenaltyMPC_dzaff12, controlPenaltyMPC_ubIdx12, controlPenaltyMPC_dsubaff12);
controlPenaltyMPC_LA_VSUB3_12(controlPenaltyMPC_lubbysub12, controlPenaltyMPC_dsubaff12, controlPenaltyMPC_lub12, controlPenaltyMPC_dlubaff12);
controlPenaltyMPC_LA_VSUB_INDEXED_12(controlPenaltyMPC_dzaff13, controlPenaltyMPC_lbIdx13, controlPenaltyMPC_rilb13, controlPenaltyMPC_dslbaff13);
controlPenaltyMPC_LA_VSUB3_12(controlPenaltyMPC_llbbyslb13, controlPenaltyMPC_dslbaff13, controlPenaltyMPC_llb13, controlPenaltyMPC_dllbaff13);
controlPenaltyMPC_LA_VSUB2_INDEXED_12(controlPenaltyMPC_riub13, controlPenaltyMPC_dzaff13, controlPenaltyMPC_ubIdx13, controlPenaltyMPC_dsubaff13);
controlPenaltyMPC_LA_VSUB3_12(controlPenaltyMPC_lubbysub13, controlPenaltyMPC_dsubaff13, controlPenaltyMPC_lub13, controlPenaltyMPC_dlubaff13);
controlPenaltyMPC_LA_VSUB_INDEXED_12(controlPenaltyMPC_dzaff14, controlPenaltyMPC_lbIdx14, controlPenaltyMPC_rilb14, controlPenaltyMPC_dslbaff14);
controlPenaltyMPC_LA_VSUB3_12(controlPenaltyMPC_llbbyslb14, controlPenaltyMPC_dslbaff14, controlPenaltyMPC_llb14, controlPenaltyMPC_dllbaff14);
controlPenaltyMPC_LA_VSUB2_INDEXED_6(controlPenaltyMPC_riub14, controlPenaltyMPC_dzaff14, controlPenaltyMPC_ubIdx14, controlPenaltyMPC_dsubaff14);
controlPenaltyMPC_LA_VSUB3_6(controlPenaltyMPC_lubbysub14, controlPenaltyMPC_dsubaff14, controlPenaltyMPC_lub14, controlPenaltyMPC_dlubaff14);
info->lsit_aff = controlPenaltyMPC_LINESEARCH_BACKTRACKING_AFFINE(controlPenaltyMPC_l, controlPenaltyMPC_s, controlPenaltyMPC_dl_aff, controlPenaltyMPC_ds_aff, &info->step_aff, &info->mu_aff);
if( info->lsit_aff == controlPenaltyMPC_NOPROGRESS ){
exitcode = controlPenaltyMPC_NOPROGRESS; break;
}
sigma_3rdroot = info->mu_aff / info->mu;
info->sigma = sigma_3rdroot*sigma_3rdroot*sigma_3rdroot;
musigma = info->mu * info->sigma;
controlPenaltyMPC_LA_VSUB5_354(controlPenaltyMPC_ds_aff, controlPenaltyMPC_dl_aff, musigma, controlPenaltyMPC_ccrhs);
controlPenaltyMPC_LA_VSUB6_INDEXED_12_12_12(controlPenaltyMPC_ccrhsub00, controlPenaltyMPC_sub00, controlPenaltyMPC_ubIdx00, controlPenaltyMPC_ccrhsl00, controlPenaltyMPC_slb00, controlPenaltyMPC_lbIdx00, controlPenaltyMPC_rd00);
controlPenaltyMPC_LA_VSUB6_INDEXED_12_12_12(controlPenaltyMPC_ccrhsub01, controlPenaltyMPC_sub01, controlPenaltyMPC_ubIdx01, controlPenaltyMPC_ccrhsl01, controlPenaltyMPC_slb01, controlPenaltyMPC_lbIdx01, controlPenaltyMPC_rd01);
controlPenaltyMPC_LA_DIAG_FORWARDSUB_12(controlPenaltyMPC_Phi00, controlPenaltyMPC_rd00, controlPenaltyMPC_Lbyrd00);
controlPenaltyMPC_LA_DIAG_FORWARDSUB_12(controlPenaltyMPC_Phi01, controlPenaltyMPC_rd01, controlPenaltyMPC_Lbyrd01);
controlPenaltyMPC_LA_DIAGZERO_MVM_6(controlPenaltyMPC_W00, controlPenaltyMPC_Lbyrd00, controlPenaltyMPC_beta00);
controlPenaltyMPC_LA_DENSE_FORWARDSUB_6(controlPenaltyMPC_Ld00, controlPenaltyMPC_beta00, controlPenaltyMPC_yy00);
controlPenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_6_12_12(controlPenaltyMPC_V00, controlPenaltyMPC_Lbyrd00, controlPenaltyMPC_W01, controlPenaltyMPC_Lbyrd01, controlPenaltyMPC_beta01);
controlPenaltyMPC_LA_DENSE_MVMSUB1_6_6(controlPenaltyMPC_Lsd01, controlPenaltyMPC_yy00, controlPenaltyMPC_beta01, controlPenaltyMPC_bmy01);
controlPenaltyMPC_LA_DENSE_FORWARDSUB_6(controlPenaltyMPC_Ld01, controlPenaltyMPC_bmy01, controlPenaltyMPC_yy01);
controlPenaltyMPC_LA_VSUB6_INDEXED_12_12_12(controlPenaltyMPC_ccrhsub02, controlPenaltyMPC_sub02, controlPenaltyMPC_ubIdx02, controlPenaltyMPC_ccrhsl02, controlPenaltyMPC_slb02, controlPenaltyMPC_lbIdx02, controlPenaltyMPC_rd02);
controlPenaltyMPC_LA_DIAG_FORWARDSUB_12(controlPenaltyMPC_Phi02, controlPenaltyMPC_rd02, controlPenaltyMPC_Lbyrd02);
controlPenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_6_12_12(controlPenaltyMPC_V01, controlPenaltyMPC_Lbyrd01, controlPenaltyMPC_W02, controlPenaltyMPC_Lbyrd02, controlPenaltyMPC_beta02);
controlPenaltyMPC_LA_DENSE_MVMSUB1_6_6(controlPenaltyMPC_Lsd02, controlPenaltyMPC_yy01, controlPenaltyMPC_beta02, controlPenaltyMPC_bmy02);
controlPenaltyMPC_LA_DENSE_FORWARDSUB_6(controlPenaltyMPC_Ld02, controlPenaltyMPC_bmy02, controlPenaltyMPC_yy02);
controlPenaltyMPC_LA_VSUB6_INDEXED_12_12_12(controlPenaltyMPC_ccrhsub03, controlPenaltyMPC_sub03, controlPenaltyMPC_ubIdx03, controlPenaltyMPC_ccrhsl03, controlPenaltyMPC_slb03, controlPenaltyMPC_lbIdx03, controlPenaltyMPC_rd03);
controlPenaltyMPC_LA_DIAG_FORWARDSUB_12(controlPenaltyMPC_Phi03, controlPenaltyMPC_rd03, controlPenaltyMPC_Lbyrd03);
controlPenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_6_12_12(controlPenaltyMPC_V02, controlPenaltyMPC_Lbyrd02, controlPenaltyMPC_W03, controlPenaltyMPC_Lbyrd03, controlPenaltyMPC_beta03);
controlPenaltyMPC_LA_DENSE_MVMSUB1_6_6(controlPenaltyMPC_Lsd03, controlPenaltyMPC_yy02, controlPenaltyMPC_beta03, controlPenaltyMPC_bmy03);
controlPenaltyMPC_LA_DENSE_FORWARDSUB_6(controlPenaltyMPC_Ld03, controlPenaltyMPC_bmy03, controlPenaltyMPC_yy03);
controlPenaltyMPC_LA_VSUB6_INDEXED_12_12_12(controlPenaltyMPC_ccrhsub04, controlPenaltyMPC_sub04, controlPenaltyMPC_ubIdx04, controlPenaltyMPC_ccrhsl04, controlPenaltyMPC_slb04, controlPenaltyMPC_lbIdx04, controlPenaltyMPC_rd04);
controlPenaltyMPC_LA_DIAG_FORWARDSUB_12(controlPenaltyMPC_Phi04, controlPenaltyMPC_rd04, controlPenaltyMPC_Lbyrd04);
controlPenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_6_12_12(controlPenaltyMPC_V03, controlPenaltyMPC_Lbyrd03, controlPenaltyMPC_W04, controlPenaltyMPC_Lbyrd04, controlPenaltyMPC_beta04);
controlPenaltyMPC_LA_DENSE_MVMSUB1_6_6(controlPenaltyMPC_Lsd04, controlPenaltyMPC_yy03, controlPenaltyMPC_beta04, controlPenaltyMPC_bmy04);
controlPenaltyMPC_LA_DENSE_FORWARDSUB_6(controlPenaltyMPC_Ld04, controlPenaltyMPC_bmy04, controlPenaltyMPC_yy04);
controlPenaltyMPC_LA_VSUB6_INDEXED_12_12_12(controlPenaltyMPC_ccrhsub05, controlPenaltyMPC_sub05, controlPenaltyMPC_ubIdx05, controlPenaltyMPC_ccrhsl05, controlPenaltyMPC_slb05, controlPenaltyMPC_lbIdx05, controlPenaltyMPC_rd05);
controlPenaltyMPC_LA_DIAG_FORWARDSUB_12(controlPenaltyMPC_Phi05, controlPenaltyMPC_rd05, controlPenaltyMPC_Lbyrd05);
controlPenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_6_12_12(controlPenaltyMPC_V04, controlPenaltyMPC_Lbyrd04, controlPenaltyMPC_W05, controlPenaltyMPC_Lbyrd05, controlPenaltyMPC_beta05);
controlPenaltyMPC_LA_DENSE_MVMSUB1_6_6(controlPenaltyMPC_Lsd05, controlPenaltyMPC_yy04, controlPenaltyMPC_beta05, controlPenaltyMPC_bmy05);
controlPenaltyMPC_LA_DENSE_FORWARDSUB_6(controlPenaltyMPC_Ld05, controlPenaltyMPC_bmy05, controlPenaltyMPC_yy05);
controlPenaltyMPC_LA_VSUB6_INDEXED_12_12_12(controlPenaltyMPC_ccrhsub06, controlPenaltyMPC_sub06, controlPenaltyMPC_ubIdx06, controlPenaltyMPC_ccrhsl06, controlPenaltyMPC_slb06, controlPenaltyMPC_lbIdx06, controlPenaltyMPC_rd06);
controlPenaltyMPC_LA_DIAG_FORWARDSUB_12(controlPenaltyMPC_Phi06, controlPenaltyMPC_rd06, controlPenaltyMPC_Lbyrd06);
controlPenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_6_12_12(controlPenaltyMPC_V05, controlPenaltyMPC_Lbyrd05, controlPenaltyMPC_W06, controlPenaltyMPC_Lbyrd06, controlPenaltyMPC_beta06);
controlPenaltyMPC_LA_DENSE_MVMSUB1_6_6(controlPenaltyMPC_Lsd06, controlPenaltyMPC_yy05, controlPenaltyMPC_beta06, controlPenaltyMPC_bmy06);
controlPenaltyMPC_LA_DENSE_FORWARDSUB_6(controlPenaltyMPC_Ld06, controlPenaltyMPC_bmy06, controlPenaltyMPC_yy06);
controlPenaltyMPC_LA_VSUB6_INDEXED_12_12_12(controlPenaltyMPC_ccrhsub07, controlPenaltyMPC_sub07, controlPenaltyMPC_ubIdx07, controlPenaltyMPC_ccrhsl07, controlPenaltyMPC_slb07, controlPenaltyMPC_lbIdx07, controlPenaltyMPC_rd07);
controlPenaltyMPC_LA_DIAG_FORWARDSUB_12(controlPenaltyMPC_Phi07, controlPenaltyMPC_rd07, controlPenaltyMPC_Lbyrd07);
controlPenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_6_12_12(controlPenaltyMPC_V06, controlPenaltyMPC_Lbyrd06, controlPenaltyMPC_W07, controlPenaltyMPC_Lbyrd07, controlPenaltyMPC_beta07);
controlPenaltyMPC_LA_DENSE_MVMSUB1_6_6(controlPenaltyMPC_Lsd07, controlPenaltyMPC_yy06, controlPenaltyMPC_beta07, controlPenaltyMPC_bmy07);
controlPenaltyMPC_LA_DENSE_FORWARDSUB_6(controlPenaltyMPC_Ld07, controlPenaltyMPC_bmy07, controlPenaltyMPC_yy07);
controlPenaltyMPC_LA_VSUB6_INDEXED_12_12_12(controlPenaltyMPC_ccrhsub08, controlPenaltyMPC_sub08, controlPenaltyMPC_ubIdx08, controlPenaltyMPC_ccrhsl08, controlPenaltyMPC_slb08, controlPenaltyMPC_lbIdx08, controlPenaltyMPC_rd08);
controlPenaltyMPC_LA_DIAG_FORWARDSUB_12(controlPenaltyMPC_Phi08, controlPenaltyMPC_rd08, controlPenaltyMPC_Lbyrd08);
controlPenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_6_12_12(controlPenaltyMPC_V07, controlPenaltyMPC_Lbyrd07, controlPenaltyMPC_W08, controlPenaltyMPC_Lbyrd08, controlPenaltyMPC_beta08);
controlPenaltyMPC_LA_DENSE_MVMSUB1_6_6(controlPenaltyMPC_Lsd08, controlPenaltyMPC_yy07, controlPenaltyMPC_beta08, controlPenaltyMPC_bmy08);
controlPenaltyMPC_LA_DENSE_FORWARDSUB_6(controlPenaltyMPC_Ld08, controlPenaltyMPC_bmy08, controlPenaltyMPC_yy08);
controlPenaltyMPC_LA_VSUB6_INDEXED_12_12_12(controlPenaltyMPC_ccrhsub09, controlPenaltyMPC_sub09, controlPenaltyMPC_ubIdx09, controlPenaltyMPC_ccrhsl09, controlPenaltyMPC_slb09, controlPenaltyMPC_lbIdx09, controlPenaltyMPC_rd09);
controlPenaltyMPC_LA_DIAG_FORWARDSUB_12(controlPenaltyMPC_Phi09, controlPenaltyMPC_rd09, controlPenaltyMPC_Lbyrd09);
controlPenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_6_12_12(controlPenaltyMPC_V08, controlPenaltyMPC_Lbyrd08, controlPenaltyMPC_W09, controlPenaltyMPC_Lbyrd09, controlPenaltyMPC_beta09);
controlPenaltyMPC_LA_DENSE_MVMSUB1_6_6(controlPenaltyMPC_Lsd09, controlPenaltyMPC_yy08, controlPenaltyMPC_beta09, controlPenaltyMPC_bmy09);
controlPenaltyMPC_LA_DENSE_FORWARDSUB_6(controlPenaltyMPC_Ld09, controlPenaltyMPC_bmy09, controlPenaltyMPC_yy09);
controlPenaltyMPC_LA_VSUB6_INDEXED_12_12_12(controlPenaltyMPC_ccrhsub10, controlPenaltyMPC_sub10, controlPenaltyMPC_ubIdx10, controlPenaltyMPC_ccrhsl10, controlPenaltyMPC_slb10, controlPenaltyMPC_lbIdx10, controlPenaltyMPC_rd10);
controlPenaltyMPC_LA_DIAG_FORWARDSUB_12(controlPenaltyMPC_Phi10, controlPenaltyMPC_rd10, controlPenaltyMPC_Lbyrd10);
controlPenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_6_12_12(controlPenaltyMPC_V09, controlPenaltyMPC_Lbyrd09, controlPenaltyMPC_W10, controlPenaltyMPC_Lbyrd10, controlPenaltyMPC_beta10);
controlPenaltyMPC_LA_DENSE_MVMSUB1_6_6(controlPenaltyMPC_Lsd10, controlPenaltyMPC_yy09, controlPenaltyMPC_beta10, controlPenaltyMPC_bmy10);
controlPenaltyMPC_LA_DENSE_FORWARDSUB_6(controlPenaltyMPC_Ld10, controlPenaltyMPC_bmy10, controlPenaltyMPC_yy10);
controlPenaltyMPC_LA_VSUB6_INDEXED_12_12_12(controlPenaltyMPC_ccrhsub11, controlPenaltyMPC_sub11, controlPenaltyMPC_ubIdx11, controlPenaltyMPC_ccrhsl11, controlPenaltyMPC_slb11, controlPenaltyMPC_lbIdx11, controlPenaltyMPC_rd11);
controlPenaltyMPC_LA_DIAG_FORWARDSUB_12(controlPenaltyMPC_Phi11, controlPenaltyMPC_rd11, controlPenaltyMPC_Lbyrd11);
controlPenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_6_12_12(controlPenaltyMPC_V10, controlPenaltyMPC_Lbyrd10, controlPenaltyMPC_W11, controlPenaltyMPC_Lbyrd11, controlPenaltyMPC_beta11);
controlPenaltyMPC_LA_DENSE_MVMSUB1_6_6(controlPenaltyMPC_Lsd11, controlPenaltyMPC_yy10, controlPenaltyMPC_beta11, controlPenaltyMPC_bmy11);
controlPenaltyMPC_LA_DENSE_FORWARDSUB_6(controlPenaltyMPC_Ld11, controlPenaltyMPC_bmy11, controlPenaltyMPC_yy11);
controlPenaltyMPC_LA_VSUB6_INDEXED_12_12_12(controlPenaltyMPC_ccrhsub12, controlPenaltyMPC_sub12, controlPenaltyMPC_ubIdx12, controlPenaltyMPC_ccrhsl12, controlPenaltyMPC_slb12, controlPenaltyMPC_lbIdx12, controlPenaltyMPC_rd12);
controlPenaltyMPC_LA_DIAG_FORWARDSUB_12(controlPenaltyMPC_Phi12, controlPenaltyMPC_rd12, controlPenaltyMPC_Lbyrd12);
controlPenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_6_12_12(controlPenaltyMPC_V11, controlPenaltyMPC_Lbyrd11, controlPenaltyMPC_W12, controlPenaltyMPC_Lbyrd12, controlPenaltyMPC_beta12);
controlPenaltyMPC_LA_DENSE_MVMSUB1_6_6(controlPenaltyMPC_Lsd12, controlPenaltyMPC_yy11, controlPenaltyMPC_beta12, controlPenaltyMPC_bmy12);
controlPenaltyMPC_LA_DENSE_FORWARDSUB_6(controlPenaltyMPC_Ld12, controlPenaltyMPC_bmy12, controlPenaltyMPC_yy12);
controlPenaltyMPC_LA_VSUB6_INDEXED_12_12_12(controlPenaltyMPC_ccrhsub13, controlPenaltyMPC_sub13, controlPenaltyMPC_ubIdx13, controlPenaltyMPC_ccrhsl13, controlPenaltyMPC_slb13, controlPenaltyMPC_lbIdx13, controlPenaltyMPC_rd13);
controlPenaltyMPC_LA_DIAG_FORWARDSUB_12(controlPenaltyMPC_Phi13, controlPenaltyMPC_rd13, controlPenaltyMPC_Lbyrd13);
controlPenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_6_12_12(controlPenaltyMPC_V12, controlPenaltyMPC_Lbyrd12, controlPenaltyMPC_W13, controlPenaltyMPC_Lbyrd13, controlPenaltyMPC_beta13);
controlPenaltyMPC_LA_DENSE_MVMSUB1_6_6(controlPenaltyMPC_Lsd13, controlPenaltyMPC_yy12, controlPenaltyMPC_beta13, controlPenaltyMPC_bmy13);
controlPenaltyMPC_LA_DENSE_FORWARDSUB_6(controlPenaltyMPC_Ld13, controlPenaltyMPC_bmy13, controlPenaltyMPC_yy13);
controlPenaltyMPC_LA_VSUB6_INDEXED_12_6_12(controlPenaltyMPC_ccrhsub14, controlPenaltyMPC_sub14, controlPenaltyMPC_ubIdx14, controlPenaltyMPC_ccrhsl14, controlPenaltyMPC_slb14, controlPenaltyMPC_lbIdx14, controlPenaltyMPC_rd14);
controlPenaltyMPC_LA_DIAG_FORWARDSUB_12(controlPenaltyMPC_Phi14, controlPenaltyMPC_rd14, controlPenaltyMPC_Lbyrd14);
controlPenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_6_12_12(controlPenaltyMPC_V13, controlPenaltyMPC_Lbyrd13, controlPenaltyMPC_W14, controlPenaltyMPC_Lbyrd14, controlPenaltyMPC_beta14);
controlPenaltyMPC_LA_DENSE_MVMSUB1_6_6(controlPenaltyMPC_Lsd14, controlPenaltyMPC_yy13, controlPenaltyMPC_beta14, controlPenaltyMPC_bmy14);
controlPenaltyMPC_LA_DENSE_FORWARDSUB_6(controlPenaltyMPC_Ld14, controlPenaltyMPC_bmy14, controlPenaltyMPC_yy14);
controlPenaltyMPC_LA_DENSE_BACKWARDSUB_6(controlPenaltyMPC_Ld14, controlPenaltyMPC_yy14, controlPenaltyMPC_dvcc14);
controlPenaltyMPC_LA_DENSE_MTVMSUB_6_6(controlPenaltyMPC_Lsd14, controlPenaltyMPC_dvcc14, controlPenaltyMPC_yy13, controlPenaltyMPC_bmy13);
controlPenaltyMPC_LA_DENSE_BACKWARDSUB_6(controlPenaltyMPC_Ld13, controlPenaltyMPC_bmy13, controlPenaltyMPC_dvcc13);
controlPenaltyMPC_LA_DENSE_MTVMSUB_6_6(controlPenaltyMPC_Lsd13, controlPenaltyMPC_dvcc13, controlPenaltyMPC_yy12, controlPenaltyMPC_bmy12);
controlPenaltyMPC_LA_DENSE_BACKWARDSUB_6(controlPenaltyMPC_Ld12, controlPenaltyMPC_bmy12, controlPenaltyMPC_dvcc12);
controlPenaltyMPC_LA_DENSE_MTVMSUB_6_6(controlPenaltyMPC_Lsd12, controlPenaltyMPC_dvcc12, controlPenaltyMPC_yy11, controlPenaltyMPC_bmy11);
controlPenaltyMPC_LA_DENSE_BACKWARDSUB_6(controlPenaltyMPC_Ld11, controlPenaltyMPC_bmy11, controlPenaltyMPC_dvcc11);
controlPenaltyMPC_LA_DENSE_MTVMSUB_6_6(controlPenaltyMPC_Lsd11, controlPenaltyMPC_dvcc11, controlPenaltyMPC_yy10, controlPenaltyMPC_bmy10);
controlPenaltyMPC_LA_DENSE_BACKWARDSUB_6(controlPenaltyMPC_Ld10, controlPenaltyMPC_bmy10, controlPenaltyMPC_dvcc10);
controlPenaltyMPC_LA_DENSE_MTVMSUB_6_6(controlPenaltyMPC_Lsd10, controlPenaltyMPC_dvcc10, controlPenaltyMPC_yy09, controlPenaltyMPC_bmy09);
controlPenaltyMPC_LA_DENSE_BACKWARDSUB_6(controlPenaltyMPC_Ld09, controlPenaltyMPC_bmy09, controlPenaltyMPC_dvcc09);
controlPenaltyMPC_LA_DENSE_MTVMSUB_6_6(controlPenaltyMPC_Lsd09, controlPenaltyMPC_dvcc09, controlPenaltyMPC_yy08, controlPenaltyMPC_bmy08);
controlPenaltyMPC_LA_DENSE_BACKWARDSUB_6(controlPenaltyMPC_Ld08, controlPenaltyMPC_bmy08, controlPenaltyMPC_dvcc08);
controlPenaltyMPC_LA_DENSE_MTVMSUB_6_6(controlPenaltyMPC_Lsd08, controlPenaltyMPC_dvcc08, controlPenaltyMPC_yy07, controlPenaltyMPC_bmy07);
controlPenaltyMPC_LA_DENSE_BACKWARDSUB_6(controlPenaltyMPC_Ld07, controlPenaltyMPC_bmy07, controlPenaltyMPC_dvcc07);
controlPenaltyMPC_LA_DENSE_MTVMSUB_6_6(controlPenaltyMPC_Lsd07, controlPenaltyMPC_dvcc07, controlPenaltyMPC_yy06, controlPenaltyMPC_bmy06);
controlPenaltyMPC_LA_DENSE_BACKWARDSUB_6(controlPenaltyMPC_Ld06, controlPenaltyMPC_bmy06, controlPenaltyMPC_dvcc06);
controlPenaltyMPC_LA_DENSE_MTVMSUB_6_6(controlPenaltyMPC_Lsd06, controlPenaltyMPC_dvcc06, controlPenaltyMPC_yy05, controlPenaltyMPC_bmy05);
controlPenaltyMPC_LA_DENSE_BACKWARDSUB_6(controlPenaltyMPC_Ld05, controlPenaltyMPC_bmy05, controlPenaltyMPC_dvcc05);
controlPenaltyMPC_LA_DENSE_MTVMSUB_6_6(controlPenaltyMPC_Lsd05, controlPenaltyMPC_dvcc05, controlPenaltyMPC_yy04, controlPenaltyMPC_bmy04);
controlPenaltyMPC_LA_DENSE_BACKWARDSUB_6(controlPenaltyMPC_Ld04, controlPenaltyMPC_bmy04, controlPenaltyMPC_dvcc04);
controlPenaltyMPC_LA_DENSE_MTVMSUB_6_6(controlPenaltyMPC_Lsd04, controlPenaltyMPC_dvcc04, controlPenaltyMPC_yy03, controlPenaltyMPC_bmy03);
controlPenaltyMPC_LA_DENSE_BACKWARDSUB_6(controlPenaltyMPC_Ld03, controlPenaltyMPC_bmy03, controlPenaltyMPC_dvcc03);
controlPenaltyMPC_LA_DENSE_MTVMSUB_6_6(controlPenaltyMPC_Lsd03, controlPenaltyMPC_dvcc03, controlPenaltyMPC_yy02, controlPenaltyMPC_bmy02);
controlPenaltyMPC_LA_DENSE_BACKWARDSUB_6(controlPenaltyMPC_Ld02, controlPenaltyMPC_bmy02, controlPenaltyMPC_dvcc02);
controlPenaltyMPC_LA_DENSE_MTVMSUB_6_6(controlPenaltyMPC_Lsd02, controlPenaltyMPC_dvcc02, controlPenaltyMPC_yy01, controlPenaltyMPC_bmy01);
controlPenaltyMPC_LA_DENSE_BACKWARDSUB_6(controlPenaltyMPC_Ld01, controlPenaltyMPC_bmy01, controlPenaltyMPC_dvcc01);
controlPenaltyMPC_LA_DENSE_MTVMSUB_6_6(controlPenaltyMPC_Lsd01, controlPenaltyMPC_dvcc01, controlPenaltyMPC_yy00, controlPenaltyMPC_bmy00);
controlPenaltyMPC_LA_DENSE_BACKWARDSUB_6(controlPenaltyMPC_Ld00, controlPenaltyMPC_bmy00, controlPenaltyMPC_dvcc00);
controlPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_6_12_6(controlPenaltyMPC_C00, controlPenaltyMPC_dvcc01, controlPenaltyMPC_D00, controlPenaltyMPC_dvcc00, controlPenaltyMPC_grad_eq00);
controlPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_6_12_6(controlPenaltyMPC_C00, controlPenaltyMPC_dvcc02, controlPenaltyMPC_D01, controlPenaltyMPC_dvcc01, controlPenaltyMPC_grad_eq01);
controlPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_6_12_6(controlPenaltyMPC_C00, controlPenaltyMPC_dvcc03, controlPenaltyMPC_D01, controlPenaltyMPC_dvcc02, controlPenaltyMPC_grad_eq02);
controlPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_6_12_6(controlPenaltyMPC_C00, controlPenaltyMPC_dvcc04, controlPenaltyMPC_D01, controlPenaltyMPC_dvcc03, controlPenaltyMPC_grad_eq03);
controlPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_6_12_6(controlPenaltyMPC_C00, controlPenaltyMPC_dvcc05, controlPenaltyMPC_D01, controlPenaltyMPC_dvcc04, controlPenaltyMPC_grad_eq04);
controlPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_6_12_6(controlPenaltyMPC_C00, controlPenaltyMPC_dvcc06, controlPenaltyMPC_D01, controlPenaltyMPC_dvcc05, controlPenaltyMPC_grad_eq05);
controlPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_6_12_6(controlPenaltyMPC_C00, controlPenaltyMPC_dvcc07, controlPenaltyMPC_D01, controlPenaltyMPC_dvcc06, controlPenaltyMPC_grad_eq06);
controlPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_6_12_6(controlPenaltyMPC_C00, controlPenaltyMPC_dvcc08, controlPenaltyMPC_D01, controlPenaltyMPC_dvcc07, controlPenaltyMPC_grad_eq07);
controlPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_6_12_6(controlPenaltyMPC_C00, controlPenaltyMPC_dvcc09, controlPenaltyMPC_D01, controlPenaltyMPC_dvcc08, controlPenaltyMPC_grad_eq08);
controlPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_6_12_6(controlPenaltyMPC_C00, controlPenaltyMPC_dvcc10, controlPenaltyMPC_D01, controlPenaltyMPC_dvcc09, controlPenaltyMPC_grad_eq09);
controlPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_6_12_6(controlPenaltyMPC_C00, controlPenaltyMPC_dvcc11, controlPenaltyMPC_D01, controlPenaltyMPC_dvcc10, controlPenaltyMPC_grad_eq10);
controlPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_6_12_6(controlPenaltyMPC_C00, controlPenaltyMPC_dvcc12, controlPenaltyMPC_D01, controlPenaltyMPC_dvcc11, controlPenaltyMPC_grad_eq11);
controlPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_6_12_6(controlPenaltyMPC_C00, controlPenaltyMPC_dvcc13, controlPenaltyMPC_D01, controlPenaltyMPC_dvcc12, controlPenaltyMPC_grad_eq12);
controlPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_6_12_6(controlPenaltyMPC_C00, controlPenaltyMPC_dvcc14, controlPenaltyMPC_D01, controlPenaltyMPC_dvcc13, controlPenaltyMPC_grad_eq13);
controlPenaltyMPC_LA_DIAGZERO_MTVM_6_12(controlPenaltyMPC_D01, controlPenaltyMPC_dvcc14, controlPenaltyMPC_grad_eq14);
controlPenaltyMPC_LA_VSUB_180(controlPenaltyMPC_rd, controlPenaltyMPC_grad_eq, controlPenaltyMPC_rd);
controlPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_12(controlPenaltyMPC_Phi00, controlPenaltyMPC_rd00, controlPenaltyMPC_dzcc00);
controlPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_12(controlPenaltyMPC_Phi01, controlPenaltyMPC_rd01, controlPenaltyMPC_dzcc01);
controlPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_12(controlPenaltyMPC_Phi02, controlPenaltyMPC_rd02, controlPenaltyMPC_dzcc02);
controlPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_12(controlPenaltyMPC_Phi03, controlPenaltyMPC_rd03, controlPenaltyMPC_dzcc03);
controlPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_12(controlPenaltyMPC_Phi04, controlPenaltyMPC_rd04, controlPenaltyMPC_dzcc04);
controlPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_12(controlPenaltyMPC_Phi05, controlPenaltyMPC_rd05, controlPenaltyMPC_dzcc05);
controlPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_12(controlPenaltyMPC_Phi06, controlPenaltyMPC_rd06, controlPenaltyMPC_dzcc06);
controlPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_12(controlPenaltyMPC_Phi07, controlPenaltyMPC_rd07, controlPenaltyMPC_dzcc07);
controlPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_12(controlPenaltyMPC_Phi08, controlPenaltyMPC_rd08, controlPenaltyMPC_dzcc08);
controlPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_12(controlPenaltyMPC_Phi09, controlPenaltyMPC_rd09, controlPenaltyMPC_dzcc09);
controlPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_12(controlPenaltyMPC_Phi10, controlPenaltyMPC_rd10, controlPenaltyMPC_dzcc10);
controlPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_12(controlPenaltyMPC_Phi11, controlPenaltyMPC_rd11, controlPenaltyMPC_dzcc11);
controlPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_12(controlPenaltyMPC_Phi12, controlPenaltyMPC_rd12, controlPenaltyMPC_dzcc12);
controlPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_12(controlPenaltyMPC_Phi13, controlPenaltyMPC_rd13, controlPenaltyMPC_dzcc13);
controlPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_12(controlPenaltyMPC_Phi14, controlPenaltyMPC_rd14, controlPenaltyMPC_dzcc14);
controlPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_12(controlPenaltyMPC_ccrhsl00, controlPenaltyMPC_slb00, controlPenaltyMPC_llbbyslb00, controlPenaltyMPC_dzcc00, controlPenaltyMPC_lbIdx00, controlPenaltyMPC_dllbcc00);
controlPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_12(controlPenaltyMPC_ccrhsub00, controlPenaltyMPC_sub00, controlPenaltyMPC_lubbysub00, controlPenaltyMPC_dzcc00, controlPenaltyMPC_ubIdx00, controlPenaltyMPC_dlubcc00);
controlPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_12(controlPenaltyMPC_ccrhsl01, controlPenaltyMPC_slb01, controlPenaltyMPC_llbbyslb01, controlPenaltyMPC_dzcc01, controlPenaltyMPC_lbIdx01, controlPenaltyMPC_dllbcc01);
controlPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_12(controlPenaltyMPC_ccrhsub01, controlPenaltyMPC_sub01, controlPenaltyMPC_lubbysub01, controlPenaltyMPC_dzcc01, controlPenaltyMPC_ubIdx01, controlPenaltyMPC_dlubcc01);
controlPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_12(controlPenaltyMPC_ccrhsl02, controlPenaltyMPC_slb02, controlPenaltyMPC_llbbyslb02, controlPenaltyMPC_dzcc02, controlPenaltyMPC_lbIdx02, controlPenaltyMPC_dllbcc02);
controlPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_12(controlPenaltyMPC_ccrhsub02, controlPenaltyMPC_sub02, controlPenaltyMPC_lubbysub02, controlPenaltyMPC_dzcc02, controlPenaltyMPC_ubIdx02, controlPenaltyMPC_dlubcc02);
controlPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_12(controlPenaltyMPC_ccrhsl03, controlPenaltyMPC_slb03, controlPenaltyMPC_llbbyslb03, controlPenaltyMPC_dzcc03, controlPenaltyMPC_lbIdx03, controlPenaltyMPC_dllbcc03);
controlPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_12(controlPenaltyMPC_ccrhsub03, controlPenaltyMPC_sub03, controlPenaltyMPC_lubbysub03, controlPenaltyMPC_dzcc03, controlPenaltyMPC_ubIdx03, controlPenaltyMPC_dlubcc03);
controlPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_12(controlPenaltyMPC_ccrhsl04, controlPenaltyMPC_slb04, controlPenaltyMPC_llbbyslb04, controlPenaltyMPC_dzcc04, controlPenaltyMPC_lbIdx04, controlPenaltyMPC_dllbcc04);
controlPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_12(controlPenaltyMPC_ccrhsub04, controlPenaltyMPC_sub04, controlPenaltyMPC_lubbysub04, controlPenaltyMPC_dzcc04, controlPenaltyMPC_ubIdx04, controlPenaltyMPC_dlubcc04);
controlPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_12(controlPenaltyMPC_ccrhsl05, controlPenaltyMPC_slb05, controlPenaltyMPC_llbbyslb05, controlPenaltyMPC_dzcc05, controlPenaltyMPC_lbIdx05, controlPenaltyMPC_dllbcc05);
controlPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_12(controlPenaltyMPC_ccrhsub05, controlPenaltyMPC_sub05, controlPenaltyMPC_lubbysub05, controlPenaltyMPC_dzcc05, controlPenaltyMPC_ubIdx05, controlPenaltyMPC_dlubcc05);
controlPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_12(controlPenaltyMPC_ccrhsl06, controlPenaltyMPC_slb06, controlPenaltyMPC_llbbyslb06, controlPenaltyMPC_dzcc06, controlPenaltyMPC_lbIdx06, controlPenaltyMPC_dllbcc06);
controlPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_12(controlPenaltyMPC_ccrhsub06, controlPenaltyMPC_sub06, controlPenaltyMPC_lubbysub06, controlPenaltyMPC_dzcc06, controlPenaltyMPC_ubIdx06, controlPenaltyMPC_dlubcc06);
controlPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_12(controlPenaltyMPC_ccrhsl07, controlPenaltyMPC_slb07, controlPenaltyMPC_llbbyslb07, controlPenaltyMPC_dzcc07, controlPenaltyMPC_lbIdx07, controlPenaltyMPC_dllbcc07);
controlPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_12(controlPenaltyMPC_ccrhsub07, controlPenaltyMPC_sub07, controlPenaltyMPC_lubbysub07, controlPenaltyMPC_dzcc07, controlPenaltyMPC_ubIdx07, controlPenaltyMPC_dlubcc07);
controlPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_12(controlPenaltyMPC_ccrhsl08, controlPenaltyMPC_slb08, controlPenaltyMPC_llbbyslb08, controlPenaltyMPC_dzcc08, controlPenaltyMPC_lbIdx08, controlPenaltyMPC_dllbcc08);
controlPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_12(controlPenaltyMPC_ccrhsub08, controlPenaltyMPC_sub08, controlPenaltyMPC_lubbysub08, controlPenaltyMPC_dzcc08, controlPenaltyMPC_ubIdx08, controlPenaltyMPC_dlubcc08);
controlPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_12(controlPenaltyMPC_ccrhsl09, controlPenaltyMPC_slb09, controlPenaltyMPC_llbbyslb09, controlPenaltyMPC_dzcc09, controlPenaltyMPC_lbIdx09, controlPenaltyMPC_dllbcc09);
controlPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_12(controlPenaltyMPC_ccrhsub09, controlPenaltyMPC_sub09, controlPenaltyMPC_lubbysub09, controlPenaltyMPC_dzcc09, controlPenaltyMPC_ubIdx09, controlPenaltyMPC_dlubcc09);
controlPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_12(controlPenaltyMPC_ccrhsl10, controlPenaltyMPC_slb10, controlPenaltyMPC_llbbyslb10, controlPenaltyMPC_dzcc10, controlPenaltyMPC_lbIdx10, controlPenaltyMPC_dllbcc10);
controlPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_12(controlPenaltyMPC_ccrhsub10, controlPenaltyMPC_sub10, controlPenaltyMPC_lubbysub10, controlPenaltyMPC_dzcc10, controlPenaltyMPC_ubIdx10, controlPenaltyMPC_dlubcc10);
controlPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_12(controlPenaltyMPC_ccrhsl11, controlPenaltyMPC_slb11, controlPenaltyMPC_llbbyslb11, controlPenaltyMPC_dzcc11, controlPenaltyMPC_lbIdx11, controlPenaltyMPC_dllbcc11);
controlPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_12(controlPenaltyMPC_ccrhsub11, controlPenaltyMPC_sub11, controlPenaltyMPC_lubbysub11, controlPenaltyMPC_dzcc11, controlPenaltyMPC_ubIdx11, controlPenaltyMPC_dlubcc11);
controlPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_12(controlPenaltyMPC_ccrhsl12, controlPenaltyMPC_slb12, controlPenaltyMPC_llbbyslb12, controlPenaltyMPC_dzcc12, controlPenaltyMPC_lbIdx12, controlPenaltyMPC_dllbcc12);
controlPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_12(controlPenaltyMPC_ccrhsub12, controlPenaltyMPC_sub12, controlPenaltyMPC_lubbysub12, controlPenaltyMPC_dzcc12, controlPenaltyMPC_ubIdx12, controlPenaltyMPC_dlubcc12);
controlPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_12(controlPenaltyMPC_ccrhsl13, controlPenaltyMPC_slb13, controlPenaltyMPC_llbbyslb13, controlPenaltyMPC_dzcc13, controlPenaltyMPC_lbIdx13, controlPenaltyMPC_dllbcc13);
controlPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_12(controlPenaltyMPC_ccrhsub13, controlPenaltyMPC_sub13, controlPenaltyMPC_lubbysub13, controlPenaltyMPC_dzcc13, controlPenaltyMPC_ubIdx13, controlPenaltyMPC_dlubcc13);
controlPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_12(controlPenaltyMPC_ccrhsl14, controlPenaltyMPC_slb14, controlPenaltyMPC_llbbyslb14, controlPenaltyMPC_dzcc14, controlPenaltyMPC_lbIdx14, controlPenaltyMPC_dllbcc14);
controlPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_6(controlPenaltyMPC_ccrhsub14, controlPenaltyMPC_sub14, controlPenaltyMPC_lubbysub14, controlPenaltyMPC_dzcc14, controlPenaltyMPC_ubIdx14, controlPenaltyMPC_dlubcc14);
controlPenaltyMPC_LA_VSUB7_354(controlPenaltyMPC_l, controlPenaltyMPC_ccrhs, controlPenaltyMPC_s, controlPenaltyMPC_dl_cc, controlPenaltyMPC_ds_cc);
controlPenaltyMPC_LA_VADD_180(controlPenaltyMPC_dz_cc, controlPenaltyMPC_dz_aff);
controlPenaltyMPC_LA_VADD_90(controlPenaltyMPC_dv_cc, controlPenaltyMPC_dv_aff);
controlPenaltyMPC_LA_VADD_354(controlPenaltyMPC_dl_cc, controlPenaltyMPC_dl_aff);
controlPenaltyMPC_LA_VADD_354(controlPenaltyMPC_ds_cc, controlPenaltyMPC_ds_aff);
info->lsit_cc = controlPenaltyMPC_LINESEARCH_BACKTRACKING_COMBINED(controlPenaltyMPC_z, controlPenaltyMPC_v, controlPenaltyMPC_l, controlPenaltyMPC_s, controlPenaltyMPC_dz_cc, controlPenaltyMPC_dv_cc, controlPenaltyMPC_dl_cc, controlPenaltyMPC_ds_cc, &info->step_cc, &info->mu);
if( info->lsit_cc == controlPenaltyMPC_NOPROGRESS ){
exitcode = controlPenaltyMPC_NOPROGRESS; break;
}
info->it++;
}
output->z1[0] = controlPenaltyMPC_z00[0];
output->z1[1] = controlPenaltyMPC_z00[1];
output->z1[2] = controlPenaltyMPC_z00[2];
output->z1[3] = controlPenaltyMPC_z00[3];
output->z1[4] = controlPenaltyMPC_z00[4];
output->z1[5] = controlPenaltyMPC_z00[5];
output->z1[6] = controlPenaltyMPC_z00[6];
output->z1[7] = controlPenaltyMPC_z00[7];
output->z1[8] = controlPenaltyMPC_z00[8];
output->z1[9] = controlPenaltyMPC_z00[9];
output->z1[10] = controlPenaltyMPC_z00[10];
output->z1[11] = controlPenaltyMPC_z00[11];
output->z2[0] = controlPenaltyMPC_z01[0];
output->z2[1] = controlPenaltyMPC_z01[1];
output->z2[2] = controlPenaltyMPC_z01[2];
output->z2[3] = controlPenaltyMPC_z01[3];
output->z2[4] = controlPenaltyMPC_z01[4];
output->z2[5] = controlPenaltyMPC_z01[5];
output->z2[6] = controlPenaltyMPC_z01[6];
output->z2[7] = controlPenaltyMPC_z01[7];
output->z2[8] = controlPenaltyMPC_z01[8];
output->z2[9] = controlPenaltyMPC_z01[9];
output->z2[10] = controlPenaltyMPC_z01[10];
output->z2[11] = controlPenaltyMPC_z01[11];
output->z3[0] = controlPenaltyMPC_z02[0];
output->z3[1] = controlPenaltyMPC_z02[1];
output->z3[2] = controlPenaltyMPC_z02[2];
output->z3[3] = controlPenaltyMPC_z02[3];
output->z3[4] = controlPenaltyMPC_z02[4];
output->z3[5] = controlPenaltyMPC_z02[5];
output->z3[6] = controlPenaltyMPC_z02[6];
output->z3[7] = controlPenaltyMPC_z02[7];
output->z3[8] = controlPenaltyMPC_z02[8];
output->z3[9] = controlPenaltyMPC_z02[9];
output->z3[10] = controlPenaltyMPC_z02[10];
output->z3[11] = controlPenaltyMPC_z02[11];
output->z4[0] = controlPenaltyMPC_z03[0];
output->z4[1] = controlPenaltyMPC_z03[1];
output->z4[2] = controlPenaltyMPC_z03[2];
output->z4[3] = controlPenaltyMPC_z03[3];
output->z4[4] = controlPenaltyMPC_z03[4];
output->z4[5] = controlPenaltyMPC_z03[5];
output->z4[6] = controlPenaltyMPC_z03[6];
output->z4[7] = controlPenaltyMPC_z03[7];
output->z4[8] = controlPenaltyMPC_z03[8];
output->z4[9] = controlPenaltyMPC_z03[9];
output->z4[10] = controlPenaltyMPC_z03[10];
output->z4[11] = controlPenaltyMPC_z03[11];
output->z5[0] = controlPenaltyMPC_z04[0];
output->z5[1] = controlPenaltyMPC_z04[1];
output->z5[2] = controlPenaltyMPC_z04[2];
output->z5[3] = controlPenaltyMPC_z04[3];
output->z5[4] = controlPenaltyMPC_z04[4];
output->z5[5] = controlPenaltyMPC_z04[5];
output->z5[6] = controlPenaltyMPC_z04[6];
output->z5[7] = controlPenaltyMPC_z04[7];
output->z5[8] = controlPenaltyMPC_z04[8];
output->z5[9] = controlPenaltyMPC_z04[9];
output->z5[10] = controlPenaltyMPC_z04[10];
output->z5[11] = controlPenaltyMPC_z04[11];
output->z6[0] = controlPenaltyMPC_z05[0];
output->z6[1] = controlPenaltyMPC_z05[1];
output->z6[2] = controlPenaltyMPC_z05[2];
output->z6[3] = controlPenaltyMPC_z05[3];
output->z6[4] = controlPenaltyMPC_z05[4];
output->z6[5] = controlPenaltyMPC_z05[5];
output->z6[6] = controlPenaltyMPC_z05[6];
output->z6[7] = controlPenaltyMPC_z05[7];
output->z6[8] = controlPenaltyMPC_z05[8];
output->z6[9] = controlPenaltyMPC_z05[9];
output->z6[10] = controlPenaltyMPC_z05[10];
output->z6[11] = controlPenaltyMPC_z05[11];
output->z7[0] = controlPenaltyMPC_z06[0];
output->z7[1] = controlPenaltyMPC_z06[1];
output->z7[2] = controlPenaltyMPC_z06[2];
output->z7[3] = controlPenaltyMPC_z06[3];
output->z7[4] = controlPenaltyMPC_z06[4];
output->z7[5] = controlPenaltyMPC_z06[5];
output->z7[6] = controlPenaltyMPC_z06[6];
output->z7[7] = controlPenaltyMPC_z06[7];
output->z7[8] = controlPenaltyMPC_z06[8];
output->z7[9] = controlPenaltyMPC_z06[9];
output->z7[10] = controlPenaltyMPC_z06[10];
output->z7[11] = controlPenaltyMPC_z06[11];
output->z8[0] = controlPenaltyMPC_z07[0];
output->z8[1] = controlPenaltyMPC_z07[1];
output->z8[2] = controlPenaltyMPC_z07[2];
output->z8[3] = controlPenaltyMPC_z07[3];
output->z8[4] = controlPenaltyMPC_z07[4];
output->z8[5] = controlPenaltyMPC_z07[5];
output->z8[6] = controlPenaltyMPC_z07[6];
output->z8[7] = controlPenaltyMPC_z07[7];
output->z8[8] = controlPenaltyMPC_z07[8];
output->z8[9] = controlPenaltyMPC_z07[9];
output->z8[10] = controlPenaltyMPC_z07[10];
output->z8[11] = controlPenaltyMPC_z07[11];
output->z9[0] = controlPenaltyMPC_z08[0];
output->z9[1] = controlPenaltyMPC_z08[1];
output->z9[2] = controlPenaltyMPC_z08[2];
output->z9[3] = controlPenaltyMPC_z08[3];
output->z9[4] = controlPenaltyMPC_z08[4];
output->z9[5] = controlPenaltyMPC_z08[5];
output->z9[6] = controlPenaltyMPC_z08[6];
output->z9[7] = controlPenaltyMPC_z08[7];
output->z9[8] = controlPenaltyMPC_z08[8];
output->z9[9] = controlPenaltyMPC_z08[9];
output->z9[10] = controlPenaltyMPC_z08[10];
output->z9[11] = controlPenaltyMPC_z08[11];
output->z10[0] = controlPenaltyMPC_z09[0];
output->z10[1] = controlPenaltyMPC_z09[1];
output->z10[2] = controlPenaltyMPC_z09[2];
output->z10[3] = controlPenaltyMPC_z09[3];
output->z10[4] = controlPenaltyMPC_z09[4];
output->z10[5] = controlPenaltyMPC_z09[5];
output->z10[6] = controlPenaltyMPC_z09[6];
output->z10[7] = controlPenaltyMPC_z09[7];
output->z10[8] = controlPenaltyMPC_z09[8];
output->z10[9] = controlPenaltyMPC_z09[9];
output->z10[10] = controlPenaltyMPC_z09[10];
output->z10[11] = controlPenaltyMPC_z09[11];
output->z11[0] = controlPenaltyMPC_z10[0];
output->z11[1] = controlPenaltyMPC_z10[1];
output->z11[2] = controlPenaltyMPC_z10[2];
output->z11[3] = controlPenaltyMPC_z10[3];
output->z11[4] = controlPenaltyMPC_z10[4];
output->z11[5] = controlPenaltyMPC_z10[5];
output->z11[6] = controlPenaltyMPC_z10[6];
output->z11[7] = controlPenaltyMPC_z10[7];
output->z11[8] = controlPenaltyMPC_z10[8];
output->z11[9] = controlPenaltyMPC_z10[9];
output->z11[10] = controlPenaltyMPC_z10[10];
output->z11[11] = controlPenaltyMPC_z10[11];
output->z12[0] = controlPenaltyMPC_z11[0];
output->z12[1] = controlPenaltyMPC_z11[1];
output->z12[2] = controlPenaltyMPC_z11[2];
output->z12[3] = controlPenaltyMPC_z11[3];
output->z12[4] = controlPenaltyMPC_z11[4];
output->z12[5] = controlPenaltyMPC_z11[5];
output->z12[6] = controlPenaltyMPC_z11[6];
output->z12[7] = controlPenaltyMPC_z11[7];
output->z12[8] = controlPenaltyMPC_z11[8];
output->z12[9] = controlPenaltyMPC_z11[9];
output->z12[10] = controlPenaltyMPC_z11[10];
output->z12[11] = controlPenaltyMPC_z11[11];
output->z13[0] = controlPenaltyMPC_z12[0];
output->z13[1] = controlPenaltyMPC_z12[1];
output->z13[2] = controlPenaltyMPC_z12[2];
output->z13[3] = controlPenaltyMPC_z12[3];
output->z13[4] = controlPenaltyMPC_z12[4];
output->z13[5] = controlPenaltyMPC_z12[5];
output->z13[6] = controlPenaltyMPC_z12[6];
output->z13[7] = controlPenaltyMPC_z12[7];
output->z13[8] = controlPenaltyMPC_z12[8];
output->z13[9] = controlPenaltyMPC_z12[9];
output->z13[10] = controlPenaltyMPC_z12[10];
output->z13[11] = controlPenaltyMPC_z12[11];
output->z14[0] = controlPenaltyMPC_z13[0];
output->z14[1] = controlPenaltyMPC_z13[1];
output->z14[2] = controlPenaltyMPC_z13[2];
output->z14[3] = controlPenaltyMPC_z13[3];
output->z14[4] = controlPenaltyMPC_z13[4];
output->z14[5] = controlPenaltyMPC_z13[5];
output->z14[6] = controlPenaltyMPC_z13[6];
output->z14[7] = controlPenaltyMPC_z13[7];
output->z14[8] = controlPenaltyMPC_z13[8];
output->z14[9] = controlPenaltyMPC_z13[9];
output->z14[10] = controlPenaltyMPC_z13[10];
output->z14[11] = controlPenaltyMPC_z13[11];
output->z15[0] = controlPenaltyMPC_z14[0];
output->z15[1] = controlPenaltyMPC_z14[1];
output->z15[2] = controlPenaltyMPC_z14[2];
output->z15[3] = controlPenaltyMPC_z14[3];
output->z15[4] = controlPenaltyMPC_z14[4];
output->z15[5] = controlPenaltyMPC_z14[5];

#if controlPenaltyMPC_SET_TIMING == 1
info->solvetime = controlPenaltyMPC_toc(&solvertimer);
#if controlPenaltyMPC_SET_PRINTLEVEL > 0 && controlPenaltyMPC_SET_TIMING == 1
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
