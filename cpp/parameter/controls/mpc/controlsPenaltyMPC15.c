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

#include "controlsPenaltyMPC.h"

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
 * Initializes a vector of length 28 with a value.
 */
void controlsPenaltyMPC_LA_INITIALIZEVECTOR_28(controlsPenaltyMPC_FLOAT* vec, controlsPenaltyMPC_FLOAT value)
{
	int i;
	for( i=0; i<28; i++ )
	{
		vec[i] = value;
	}
}


/*
 * Initializes a vector of length 56 with a value.
 */
void controlsPenaltyMPC_LA_INITIALIZEVECTOR_56(controlsPenaltyMPC_FLOAT* vec, controlsPenaltyMPC_FLOAT value)
{
	int i;
	for( i=0; i<56; i++ )
	{
		vec[i] = value;
	}
}


/* 
 * Calculates a dot product and adds it to a variable: z += x'*y; 
 * This function is for vectors of length 56.
 */
void controlsPenaltyMPC_LA_DOTACC_56(controlsPenaltyMPC_FLOAT *x, controlsPenaltyMPC_FLOAT *y, controlsPenaltyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<56; i++ ){
		*z += x[i]*y[i];
	}
}


/*
 * Calculates the gradient and the value for a quadratic function 0.5*z'*H*z + f'*z
 *
 * INPUTS:     H  - Symmetric Hessian, diag matrix of size [2 x 2]
 *             f  - column vector of size 2
 *             z  - column vector of size 2
 *
 * OUTPUTS: grad  - gradient at z (= H*z + f), column vector of size 2
 *          value <-- value + 0.5*z'*H*z + f'*z (value will be modified)
 */
void controlsPenaltyMPC_LA_DIAG_QUADFCN_2(controlsPenaltyMPC_FLOAT* H, controlsPenaltyMPC_FLOAT* f, controlsPenaltyMPC_FLOAT* z, controlsPenaltyMPC_FLOAT* grad, controlsPenaltyMPC_FLOAT* value)
{
	int i;
	controlsPenaltyMPC_FLOAT hz;	
	for( i=0; i<2; i++){
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
void controlsPenaltyMPC_LA_DIAGZERO_MVMSUB6_2(controlsPenaltyMPC_FLOAT *B, controlsPenaltyMPC_FLOAT *u, controlsPenaltyMPC_FLOAT *b, controlsPenaltyMPC_FLOAT *l, controlsPenaltyMPC_FLOAT *r, controlsPenaltyMPC_FLOAT *z, controlsPenaltyMPC_FLOAT *y)
{
	int i;
	controlsPenaltyMPC_FLOAT Bu[2];
	controlsPenaltyMPC_FLOAT norm = *y;
	controlsPenaltyMPC_FLOAT lr = 0;

	/* do A*x + B*u first */
	for( i=0; i<2; i++ ){
		Bu[i] = B[i]*u[i];
	}	

	for( i=0; i<2; i++ ){
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
void controlsPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_2_2(controlsPenaltyMPC_FLOAT *A, controlsPenaltyMPC_FLOAT *x, controlsPenaltyMPC_FLOAT *B, controlsPenaltyMPC_FLOAT *u, controlsPenaltyMPC_FLOAT *b, controlsPenaltyMPC_FLOAT *l, controlsPenaltyMPC_FLOAT *r, controlsPenaltyMPC_FLOAT *z, controlsPenaltyMPC_FLOAT *y)
{
	int i;
	int j;
	int k = 0;
	controlsPenaltyMPC_FLOAT AxBu[2];
	controlsPenaltyMPC_FLOAT norm = *y;
	controlsPenaltyMPC_FLOAT lr = 0;

	/* do A*x + B*u first */
	for( i=0; i<2; i++ ){
		AxBu[i] = A[k++]*x[0] + B[i]*u[i];
	}	

	for( j=1; j<2; j++ ){		
		for( i=0; i<2; i++ ){
			AxBu[i] += A[k++]*x[j];
		}
	}

	for( i=0; i<2; i++ ){
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
 * where A is of size [2 x 2] and stored in column major format.
 * and B is of size [2 x 2] and stored in diagzero format
 * Note the transposes of A and B!
 */
void controlsPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlsPenaltyMPC_FLOAT *A, controlsPenaltyMPC_FLOAT *x, controlsPenaltyMPC_FLOAT *B, controlsPenaltyMPC_FLOAT *y, controlsPenaltyMPC_FLOAT *z)
{
	int i;
	int j;
	int k = 0;
	for( i=0; i<2; i++ ){
		z[i] = 0;
		for( j=0; j<2; j++ ){
			z[i] += A[k++]*x[j];
		}
		z[i] += B[i]*y[i];
	}
	for( i=2 ;i<2; i++ ){
		z[i] = 0;
		for( j=0; j<2; j++ ){
			z[i] += A[k++]*x[j];
		}
	}
}


/*
 * Matrix vector multiplication y = M'*x where M is of size [2 x 2]
 * and stored in diagzero format. Note the transpose of M!
 */
void controlsPenaltyMPC_LA_DIAGZERO_MTVM_2_2(controlsPenaltyMPC_FLOAT *M, controlsPenaltyMPC_FLOAT *x, controlsPenaltyMPC_FLOAT *y)
{
	int i;
	for( i=0; i<2; i++ ){
		y[i] = M[i]*x[i];
	}
}


/*
 * Vector subtraction and addition.
 *	 Input: five vectors t, tidx, u, v, w and two scalars z and r
 *	 Output: y = t(tidx) - u + w
 *           z = z - v'*x;
 *           r = max([norm(y,inf), z]);
 * for vectors of length 2. Output z is of course scalar.
 */
void controlsPenaltyMPC_LA_VSUBADD3_2(controlsPenaltyMPC_FLOAT* t, controlsPenaltyMPC_FLOAT* u, int* uidx, controlsPenaltyMPC_FLOAT* v, controlsPenaltyMPC_FLOAT* w, controlsPenaltyMPC_FLOAT* y, controlsPenaltyMPC_FLOAT* z, controlsPenaltyMPC_FLOAT* r)
{
	int i;
	controlsPenaltyMPC_FLOAT norm = *r;
	controlsPenaltyMPC_FLOAT vx = 0;
	controlsPenaltyMPC_FLOAT x;
	for( i=0; i<2; i++){
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
 * for vectors of length 2. Output z is of course scalar.
 */
void controlsPenaltyMPC_LA_VSUBADD2_2(controlsPenaltyMPC_FLOAT* t, int* tidx, controlsPenaltyMPC_FLOAT* u, controlsPenaltyMPC_FLOAT* v, controlsPenaltyMPC_FLOAT* w, controlsPenaltyMPC_FLOAT* y, controlsPenaltyMPC_FLOAT* z, controlsPenaltyMPC_FLOAT* r)
{
	int i;
	controlsPenaltyMPC_FLOAT norm = *r;
	controlsPenaltyMPC_FLOAT vx = 0;
	controlsPenaltyMPC_FLOAT x;
	for( i=0; i<2; i++){
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
 * Special function for box constraints of length 2
 * Returns also L/S, a value that is often used elsewhere.
 */
void controlsPenaltyMPC_LA_INEQ_B_GRAD_2_2_2(controlsPenaltyMPC_FLOAT *lu, controlsPenaltyMPC_FLOAT *su, controlsPenaltyMPC_FLOAT *ru, controlsPenaltyMPC_FLOAT *ll, controlsPenaltyMPC_FLOAT *sl, controlsPenaltyMPC_FLOAT *rl, int* lbIdx, int* ubIdx, controlsPenaltyMPC_FLOAT *grad, controlsPenaltyMPC_FLOAT *lubysu, controlsPenaltyMPC_FLOAT *llbysl)
{
	int i;
	for( i=0; i<2; i++ ){
		grad[i] = 0;
	}
	for( i=0; i<2; i++ ){		
		llbysl[i] = ll[i] / sl[i];
		grad[lbIdx[i]] -= llbysl[i]*rl[i];
	}
	for( i=0; i<2; i++ ){
		lubysu[i] = lu[i] / su[i];
		grad[ubIdx[i]] += lubysu[i]*ru[i];
	}
}


/*
 * Addition of three vectors  z = u + w + v
 * of length 28.
 */
void controlsPenaltyMPC_LA_VVADD3_28(controlsPenaltyMPC_FLOAT *u, controlsPenaltyMPC_FLOAT *v, controlsPenaltyMPC_FLOAT *w, controlsPenaltyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<28; i++ ){
		z[i] = u[i] + v[i] + w[i];
	}
}


/*
 * Special function to compute the diagonal cholesky factorization of the 
 * positive definite augmented Hessian for block size 2.
 *
 * Inputs: - H = diagonal cost Hessian in diagonal storage format
 *         - llbysl = L / S of lower bounds
 *         - lubysu = L / S of upper bounds
 *
 * Output: Phi = sqrt(H + diag(llbysl) + diag(lubysu))
 * where Phi is stored in diagonal storage format
 */
void controlsPenaltyMPC_LA_DIAG_CHOL_ONELOOP_LBUB_2_2_2(controlsPenaltyMPC_FLOAT *H, controlsPenaltyMPC_FLOAT *llbysl, int* lbIdx, controlsPenaltyMPC_FLOAT *lubysu, int* ubIdx, controlsPenaltyMPC_FLOAT *Phi)


{
	int i;
	
	/* compute cholesky */
	for( i=0; i<2; i++ ){
		Phi[i] = H[i] + llbysl[i] + lubysu[i];

#if controlsPenaltyMPC_SET_PRINTLEVEL > 0 && defined PRINTNUMERICALWARNINGS
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
 * where A is to be computed and is of size [2 x 2],
 * B is given and of size [2 x 2], L is a diagonal
 * matrix of size 2 stored in diagonal matrix 
 * storage format. Note the transpose of L has no impact!
 *
 * Result: A in column major storage format.
 *
 */
void controlsPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_2_2(controlsPenaltyMPC_FLOAT *L, controlsPenaltyMPC_FLOAT *B, controlsPenaltyMPC_FLOAT *A)
{
    int i,j;
	 int k = 0;

	for( j=0; j<2; j++){
		for( i=0; i<2; i++){
			A[k] = B[k]/L[j];
			k++;
		}
	}

}


/**
 * Forward substitution for the matrix equation A*L' = B
 * where A is to be computed and is of size [2 x 2],
 * B is given and of size [2 x 2], L is a diagonal
 *  matrix of size 2 stored in diagonal 
 * storage format. Note the transpose of L!
 *
 * Result: A in diagonalzero storage format.
 *
 */
void controlsPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_2(controlsPenaltyMPC_FLOAT *L, controlsPenaltyMPC_FLOAT *B, controlsPenaltyMPC_FLOAT *A)
{
	int j;
    for( j=0; j<2; j++ ){   
		A[j] = B[j]/L[j];
     }
}


/**
 * Compute C = A*B' where 
 *
 *	size(A) = [2 x 2]
 *  size(B) = [2 x 2] in diagzero format
 * 
 * A and C matrices are stored in column major format.
 * 
 * 
 */
void controlsPenaltyMPC_LA_DENSE_DIAGZERO_MMTM_2_2_2(controlsPenaltyMPC_FLOAT *A, controlsPenaltyMPC_FLOAT *B, controlsPenaltyMPC_FLOAT *C)
{
    int i, j;
	
	for( i=0; i<2; i++ ){
		for( j=0; j<2; j++){
			C[j*2+i] = B[i*2+j]*A[i];
		}
	}

}


/**
 * Forward substitution to solve L*y = b where L is a
 * diagonal matrix in vector storage format.
 * 
 * The dimensions involved are 2.
 */
void controlsPenaltyMPC_LA_DIAG_FORWARDSUB_2(controlsPenaltyMPC_FLOAT *L, controlsPenaltyMPC_FLOAT *b, controlsPenaltyMPC_FLOAT *y)
{
    int i;

    for( i=0; i<2; i++ ){
		y[i] = b[i]/L[i];
    }
}


/**
 * Compute L = B*B', where L is lower triangular of size NXp1
 * and B is a diagzero matrix of size [2 x 2] in column
 * storage format.
 * 
 */
void controlsPenaltyMPC_LA_DIAGZERO_MMT_2(controlsPenaltyMPC_FLOAT *B, controlsPenaltyMPC_FLOAT *L)
{
    int i, ii, di;
    
    ii = 0; di = 0;
    for( i=0; i<2; i++ ){        
		L[ii+i] = B[i]*B[i];
        ii += ++di;
    }
}


/* 
 * Computes r = b - B*u
 * B is stored in diagzero format
 */
void controlsPenaltyMPC_LA_DIAGZERO_MVMSUB7_2(controlsPenaltyMPC_FLOAT *B, controlsPenaltyMPC_FLOAT *u, controlsPenaltyMPC_FLOAT *b, controlsPenaltyMPC_FLOAT *r)
{
	int i;

	for( i=0; i<2; i++ ){
		r[i] = b[i] - B[i]*u[i];
	}	
	
}


/**
 * Compute L = A*A' + B*B', where L is lower triangular of size NXp1
 * and A is a dense matrix of size [2 x 2] in column
 * storage format, and B is of size [2 x 2] diagonalzero
 * storage format.
 * 
 * THIS ONE HAS THE WORST ACCES PATTERN POSSIBLE. 
 * POSSIBKE FIX: PUT A AND B INTO ROW MAJOR FORMAT FIRST.
 * 
 */
void controlsPenaltyMPC_LA_DENSE_DIAGZERO_MMT2_2_2_2(controlsPenaltyMPC_FLOAT *A, controlsPenaltyMPC_FLOAT *B, controlsPenaltyMPC_FLOAT *L)
{
    int i, j, k, ii, di;
    controlsPenaltyMPC_FLOAT ltemp;
    
    ii = 0; di = 0;
    for( i=0; i<2; i++ ){        
        for( j=0; j<=i; j++ ){
            ltemp = 0; 
            for( k=0; k<2; k++ ){
                ltemp += A[k*2+i]*A[k*2+j];
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
void controlsPenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_2_2(controlsPenaltyMPC_FLOAT *A, controlsPenaltyMPC_FLOAT *x, controlsPenaltyMPC_FLOAT *B, controlsPenaltyMPC_FLOAT *u, controlsPenaltyMPC_FLOAT *b, controlsPenaltyMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;

	for( i=0; i<2; i++ ){
		r[i] = b[i] - A[k++]*x[0] - B[i]*u[i];
	}	

	for( j=1; j<2; j++ ){		
		for( i=0; i<2; i++ ){
			r[i] -= A[k++]*x[j];
		}
	}
	
}


/**
 * Cholesky factorization as above, but working on a matrix in 
 * lower triangular storage format of size 2 and outputting
 * the Cholesky factor to matrix L in lower triangular format.
 */
void controlsPenaltyMPC_LA_DENSE_CHOL_2(controlsPenaltyMPC_FLOAT *A, controlsPenaltyMPC_FLOAT *L)
{
    int i, j, k, di, dj;
	 int ii, jj;

    controlsPenaltyMPC_FLOAT l;
    controlsPenaltyMPC_FLOAT Mii;

	/* copy A to L first and then operate on L */
	/* COULD BE OPTIMIZED */
	ii=0; di=0;
	for( i=0; i<2; i++ ){
		for( j=0; j<=i; j++ ){
			L[ii+j] = A[ii+j];
		}
		ii += ++di;
	}    
	
	/* factor L */
	ii=0; di=0;
    for( i=0; i<2; i++ ){
        l = 0;
        for( k=0; k<i; k++ ){
            l += L[ii+k]*L[ii+k];
        }        
        
        Mii = L[ii+i] - l;
        
#if controlsPenaltyMPC_SET_PRINTLEVEL > 0 && defined PRINTNUMERICALWARNINGS
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
        for( j=i+1; j<2; j++ ){
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
 * The dimensions involved are 2.
 */
void controlsPenaltyMPC_LA_DENSE_FORWARDSUB_2(controlsPenaltyMPC_FLOAT *L, controlsPenaltyMPC_FLOAT *b, controlsPenaltyMPC_FLOAT *y)
{
    int i,j,ii,di;
    controlsPenaltyMPC_FLOAT yel;
            
    ii = 0; di = 0;
    for( i=0; i<2; i++ ){
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
 * where A is to be computed and is of size [2 x 2],
 * B is given and of size [2 x 2], L is a lower tri-
 * angular matrix of size 2 stored in lower triangular 
 * storage format. Note the transpose of L AND B!
 *
 * Result: A in column major storage format.
 *
 */
void controlsPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(controlsPenaltyMPC_FLOAT *L, controlsPenaltyMPC_FLOAT *B, controlsPenaltyMPC_FLOAT *A)
{
    int i,j,k,ii,di;
    controlsPenaltyMPC_FLOAT a;
    
    ii=0; di=0;
    for( j=0; j<2; j++ ){        
        for( i=0; i<2; i++ ){
            a = B[i*2+j];
            for( k=0; k<j; k++ ){
                a -= A[k*2+i]*L[ii+k];
            }    

			/* saturate for numerical stability */
			a = MIN(a, BIGM);
			a = MAX(a, -BIGM); 

			A[j*2+i] = a/L[ii+j];			
        }
        ii += ++di;
    }
}


/**
 * Compute L = L - A*A', where L is lower triangular of size 2
 * and A is a dense matrix of size [2 x 2] in column
 * storage format.
 * 
 * THIS ONE HAS THE WORST ACCES PATTERN POSSIBLE. 
 * POSSIBKE FIX: PUT A INTO ROW MAJOR FORMAT FIRST.
 * 
 */
void controlsPenaltyMPC_LA_DENSE_MMTSUB_2_2(controlsPenaltyMPC_FLOAT *A, controlsPenaltyMPC_FLOAT *L)
{
    int i, j, k, ii, di;
    controlsPenaltyMPC_FLOAT ltemp;
    
    ii = 0; di = 0;
    for( i=0; i<2; i++ ){        
        for( j=0; j<=i; j++ ){
            ltemp = 0; 
            for( k=0; k<2; k++ ){
                ltemp += A[k*2+i]*A[k*2+j];
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
void controlsPenaltyMPC_LA_DENSE_MVMSUB1_2_2(controlsPenaltyMPC_FLOAT *A, controlsPenaltyMPC_FLOAT *x, controlsPenaltyMPC_FLOAT *b, controlsPenaltyMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;

	for( i=0; i<2; i++ ){
		r[i] = b[i] - A[k++]*x[0];
	}	
	for( j=1; j<2; j++ ){		
		for( i=0; i<2; i++ ){
			r[i] -= A[k++]*x[j];
		}
	}
}


/**
 * Backward Substitution to solve L^T*x = y where L is a
 * lower triangular matrix in triangular storage format.
 * 
 * All involved dimensions are 2.
 */
void controlsPenaltyMPC_LA_DENSE_BACKWARDSUB_2(controlsPenaltyMPC_FLOAT *L, controlsPenaltyMPC_FLOAT *y, controlsPenaltyMPC_FLOAT *x)
{
    int i, ii, di, j, jj, dj;
    controlsPenaltyMPC_FLOAT xel;    
	int start = 1;
    
    /* now solve L^T*x = y by backward substitution */
    ii = start; di = 1;
    for( i=1; i>=0; i-- ){        
        xel = y[i];        
        jj = start; dj = 1;
        for( j=1; j>i; j-- ){
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
 * Matrix vector multiplication y = b - M'*x where M is of size [2 x 2]
 * and stored in column major format. Note the transpose of M!
 */
void controlsPenaltyMPC_LA_DENSE_MTVMSUB_2_2(controlsPenaltyMPC_FLOAT *A, controlsPenaltyMPC_FLOAT *x, controlsPenaltyMPC_FLOAT *b, controlsPenaltyMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0; 
	for( i=0; i<2; i++ ){
		r[i] = b[i];
		for( j=0; j<2; j++ ){
			r[i] -= A[k++]*x[j];
		}
	}
}


/*
 * Vector subtraction z = -x - y for vectors of length 28.
 */
void controlsPenaltyMPC_LA_VSUB2_28(controlsPenaltyMPC_FLOAT *x, controlsPenaltyMPC_FLOAT *y, controlsPenaltyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<28; i++){
		z[i] = -x[i] - y[i];
	}
}


/**
 * Forward-Backward-Substitution to solve L*L^T*x = b where L is a
 * diagonal matrix of size 2 in vector
 * storage format.
 */
void controlsPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_2(controlsPenaltyMPC_FLOAT *L, controlsPenaltyMPC_FLOAT *b, controlsPenaltyMPC_FLOAT *x)
{
    int i;
            
    /* solve Ly = b by forward and backward substitution */
    for( i=0; i<2; i++ ){
		x[i] = b[i]/(L[i]*L[i]);
    }
    
}


/*
 * Vector subtraction z = x(xidx) - y where y, z and xidx are of length 2,
 * and x has length 2 and is indexed through yidx.
 */
void controlsPenaltyMPC_LA_VSUB_INDEXED_2(controlsPenaltyMPC_FLOAT *x, int* xidx, controlsPenaltyMPC_FLOAT *y, controlsPenaltyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<2; i++){
		z[i] = x[xidx[i]] - y[i];
	}
}


/*
 * Vector subtraction x = -u.*v - w for vectors of length 2.
 */
void controlsPenaltyMPC_LA_VSUB3_2(controlsPenaltyMPC_FLOAT *u, controlsPenaltyMPC_FLOAT *v, controlsPenaltyMPC_FLOAT *w, controlsPenaltyMPC_FLOAT *x)
{
	int i;
	for( i=0; i<2; i++){
		x[i] = -u[i]*v[i] - w[i];
	}
}


/*
 * Vector subtraction z = -x - y(yidx) where y is of length 2
 * and z, x and yidx are of length 2.
 */
void controlsPenaltyMPC_LA_VSUB2_INDEXED_2(controlsPenaltyMPC_FLOAT *x, controlsPenaltyMPC_FLOAT *y, int* yidx, controlsPenaltyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<2; i++){
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
 * controlsPenaltyMPC_NOPROGRESS (should be negative).
 */
int controlsPenaltyMPC_LINESEARCH_BACKTRACKING_AFFINE(controlsPenaltyMPC_FLOAT *l, controlsPenaltyMPC_FLOAT *s, controlsPenaltyMPC_FLOAT *dl, controlsPenaltyMPC_FLOAT *ds, controlsPenaltyMPC_FLOAT *a, controlsPenaltyMPC_FLOAT *mu_aff)
{
    int i;
	int lsIt=1;    
    controlsPenaltyMPC_FLOAT dltemp;
    controlsPenaltyMPC_FLOAT dstemp;
    controlsPenaltyMPC_FLOAT mya = 1.0;
    controlsPenaltyMPC_FLOAT mymu;
        
    while( 1 ){                        

        /* 
         * Compute both snew and wnew together.
         * We compute also mu_affine along the way here, as the
         * values might be in registers, so it should be cheaper.
         */
        mymu = 0;
        for( i=0; i<56; i++ ){
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
        if( i == 56 ){
            break;
        } else {
            mya *= controlsPenaltyMPC_SET_LS_SCALE_AFF;
            if( mya < controlsPenaltyMPC_SET_LS_MINSTEP ){
                return controlsPenaltyMPC_NOPROGRESS;
            }
        }
    }
    
    /* return new values and iteration counter */
    *a = mya;
    *mu_aff = mymu / (controlsPenaltyMPC_FLOAT)56;
    return lsIt;
}


/*
 * Vector subtraction x = u.*v - a where a is a scalar
*  and x,u,v are vectors of length 56.
 */
void controlsPenaltyMPC_LA_VSUB5_56(controlsPenaltyMPC_FLOAT *u, controlsPenaltyMPC_FLOAT *v, controlsPenaltyMPC_FLOAT a, controlsPenaltyMPC_FLOAT *x)
{
	int i;
	for( i=0; i<56; i++){
		x[i] = u[i]*v[i] - a;
	}
}


/*
 * Computes x=0; x(uidx) += u/su; x(vidx) -= v/sv where x is of length 2,
 * u, su, uidx are of length 2 and v, sv, vidx are of length 2.
 */
void controlsPenaltyMPC_LA_VSUB6_INDEXED_2_2_2(controlsPenaltyMPC_FLOAT *u, controlsPenaltyMPC_FLOAT *su, int* uidx, controlsPenaltyMPC_FLOAT *v, controlsPenaltyMPC_FLOAT *sv, int* vidx, controlsPenaltyMPC_FLOAT *x)
{
	int i;
	for( i=0; i<2; i++ ){
		x[i] = 0;
	}
	for( i=0; i<2; i++){
		x[uidx[i]] += u[i]/su[i];
	}
	for( i=0; i<2; i++){
		x[vidx[i]] -= v[i]/sv[i];
	}
}


/* 
 * Computes r =  B*u
 * where B is stored in diagzero format
 */
void controlsPenaltyMPC_LA_DIAGZERO_MVM_2(controlsPenaltyMPC_FLOAT *B, controlsPenaltyMPC_FLOAT *u, controlsPenaltyMPC_FLOAT *r)
{
	int i;

	for( i=0; i<2; i++ ){
		r[i] = B[i]*u[i];
	}	
	
}


/* 
 * Computes r = A*x + B*u
 * where A is stored in column major format
 * and B is stored in diagzero format
 */
void controlsPenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_2_2_2(controlsPenaltyMPC_FLOAT *A, controlsPenaltyMPC_FLOAT *x, controlsPenaltyMPC_FLOAT *B, controlsPenaltyMPC_FLOAT *u, controlsPenaltyMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;

	for( i=0; i<2; i++ ){
		r[i] = A[k++]*x[0] + B[i]*u[i];
	}	

	for( j=1; j<2; j++ ){		
		for( i=0; i<2; i++ ){
			r[i] += A[k++]*x[j];
		}
	}
	
}


/*
 * Vector subtraction z = x - y for vectors of length 28.
 */
void controlsPenaltyMPC_LA_VSUB_28(controlsPenaltyMPC_FLOAT *x, controlsPenaltyMPC_FLOAT *y, controlsPenaltyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<28; i++){
		z[i] = x[i] - y[i];
	}
}


/** 
 * Computes z = -r./s - u.*y(y)
 * where all vectors except of y are of length 2 (length of y >= 2).
 */
void controlsPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_2(controlsPenaltyMPC_FLOAT *r, controlsPenaltyMPC_FLOAT *s, controlsPenaltyMPC_FLOAT *u, controlsPenaltyMPC_FLOAT *y, int* yidx, controlsPenaltyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<2; i++ ){
		z[i] = -r[i]/s[i] - u[i]*y[yidx[i]];
	}
}


/** 
 * Computes z = -r./s + u.*y(y)
 * where all vectors except of y are of length 2 (length of y >= 2).
 */
void controlsPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_2(controlsPenaltyMPC_FLOAT *r, controlsPenaltyMPC_FLOAT *s, controlsPenaltyMPC_FLOAT *u, controlsPenaltyMPC_FLOAT *y, int* yidx, controlsPenaltyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<2; i++ ){
		z[i] = -r[i]/s[i] + u[i]*y[yidx[i]];
	}
}


/*
 * Computes ds = -l.\(r + s.*dl) for vectors of length 56.
 */
void controlsPenaltyMPC_LA_VSUB7_56(controlsPenaltyMPC_FLOAT *l, controlsPenaltyMPC_FLOAT *r, controlsPenaltyMPC_FLOAT *s, controlsPenaltyMPC_FLOAT *dl, controlsPenaltyMPC_FLOAT *ds)
{
	int i;
	for( i=0; i<56; i++){
		ds[i] = -(r[i] + s[i]*dl[i])/l[i];
	}
}


/*
 * Vector addition x = x + y for vectors of length 28.
 */
void controlsPenaltyMPC_LA_VADD_28(controlsPenaltyMPC_FLOAT *x, controlsPenaltyMPC_FLOAT *y)
{
	int i;
	for( i=0; i<28; i++){
		x[i] += y[i];
	}
}


/*
 * Vector addition x = x + y for vectors of length 56.
 */
void controlsPenaltyMPC_LA_VADD_56(controlsPenaltyMPC_FLOAT *x, controlsPenaltyMPC_FLOAT *y)
{
	int i;
	for( i=0; i<56; i++){
		x[i] += y[i];
	}
}


/**
 * Backtracking line search for combined predictor/corrector step.
 * Update on variables with safety factor gamma (to keep us away from
 * boundary).
 */
int controlsPenaltyMPC_LINESEARCH_BACKTRACKING_COMBINED(controlsPenaltyMPC_FLOAT *z, controlsPenaltyMPC_FLOAT *v, controlsPenaltyMPC_FLOAT *l, controlsPenaltyMPC_FLOAT *s, controlsPenaltyMPC_FLOAT *dz, controlsPenaltyMPC_FLOAT *dv, controlsPenaltyMPC_FLOAT *dl, controlsPenaltyMPC_FLOAT *ds, controlsPenaltyMPC_FLOAT *a, controlsPenaltyMPC_FLOAT *mu)
{
    int i, lsIt=1;       
    controlsPenaltyMPC_FLOAT dltemp;
    controlsPenaltyMPC_FLOAT dstemp;    
    controlsPenaltyMPC_FLOAT a_gamma;
            
    *a = 1.0;
    while( 1 ){                        

        /* check whether search criterion is fulfilled */
        for( i=0; i<56; i++ ){
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
        if( i == 56 ){
            break;
        } else {
            *a *= controlsPenaltyMPC_SET_LS_SCALE;
            if( *a < controlsPenaltyMPC_SET_LS_MINSTEP ){
                return controlsPenaltyMPC_NOPROGRESS;
            }
        }
    }
    
    /* update variables with safety margin */
    a_gamma = (*a)*controlsPenaltyMPC_SET_LS_MAXSTEP;
    
    /* primal variables */
    for( i=0; i<28; i++ ){
        z[i] += a_gamma*dz[i];
    }
    
    /* equality constraint multipliers */
    for( i=0; i<28; i++ ){
        v[i] += a_gamma*dv[i];
    }
    
    /* inequality constraint multipliers & slacks, also update mu */
    *mu = 0;
    for( i=0; i<56; i++ ){
        dltemp = l[i] + a_gamma*dl[i]; l[i] = dltemp;
        dstemp = s[i] + a_gamma*ds[i]; s[i] = dstemp;
        *mu += dltemp*dstemp;
    }
    
    *a = a_gamma;
    *mu /= (controlsPenaltyMPC_FLOAT)56;
    return lsIt;
}




/* VARIABLE DEFINITIONS ------------------------------------------------ */
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_z[28];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_v[28];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_dz_aff[28];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_dv_aff[28];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_grad_cost[28];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_grad_eq[28];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_rd[28];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_l[56];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_s[56];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_lbys[56];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_dl_aff[56];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_ds_aff[56];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_dz_cc[28];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_dv_cc[28];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_dl_cc[56];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_ds_cc[56];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_ccrhs[56];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_grad_ineq[28];
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_z00 = controlsPenaltyMPC_z + 0;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dzaff00 = controlsPenaltyMPC_dz_aff + 0;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dzcc00 = controlsPenaltyMPC_dz_cc + 0;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_rd00 = controlsPenaltyMPC_rd + 0;
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_Lbyrd00[2];
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_grad_cost00 = controlsPenaltyMPC_grad_cost + 0;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_grad_eq00 = controlsPenaltyMPC_grad_eq + 0;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_grad_ineq00 = controlsPenaltyMPC_grad_ineq + 0;
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_ctv00[2];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_C00[4] = {0.0000000000000000E+000, 0.0000000000000000E+000, 
0.0000000000000000E+000, 0.0000000000000000E+000};
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_v00 = controlsPenaltyMPC_v + 0;
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_re00[2];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_beta00[2];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_betacc00[2];
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dvaff00 = controlsPenaltyMPC_dv_aff + 0;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dvcc00 = controlsPenaltyMPC_dv_cc + 0;
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_V00[4];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_Yd00[3];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_Ld00[3];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_yy00[2];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_bmy00[2];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_c00[2] = {0.0000000000000000E+000, 0.0000000000000000E+000};
int controlsPenaltyMPC_lbIdx00[2] = {0, 1};
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_llb00 = controlsPenaltyMPC_l + 0;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_slb00 = controlsPenaltyMPC_s + 0;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_llbbyslb00 = controlsPenaltyMPC_lbys + 0;
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_rilb00[2];
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dllbaff00 = controlsPenaltyMPC_dl_aff + 0;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dslbaff00 = controlsPenaltyMPC_ds_aff + 0;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dllbcc00 = controlsPenaltyMPC_dl_cc + 0;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dslbcc00 = controlsPenaltyMPC_ds_cc + 0;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_ccrhsl00 = controlsPenaltyMPC_ccrhs + 0;
int controlsPenaltyMPC_ubIdx00[2] = {0, 1};
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_lub00 = controlsPenaltyMPC_l + 2;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_sub00 = controlsPenaltyMPC_s + 2;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_lubbysub00 = controlsPenaltyMPC_lbys + 2;
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_riub00[2];
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dlubaff00 = controlsPenaltyMPC_dl_aff + 2;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dsubaff00 = controlsPenaltyMPC_ds_aff + 2;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dlubcc00 = controlsPenaltyMPC_dl_cc + 2;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dsubcc00 = controlsPenaltyMPC_ds_cc + 2;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_ccrhsub00 = controlsPenaltyMPC_ccrhs + 2;
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_Phi00[2];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_D00[2] = {0.0000000000000000E+000, 
0.0000000000000000E+000};
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_W00[2];
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_z01 = controlsPenaltyMPC_z + 2;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dzaff01 = controlsPenaltyMPC_dz_aff + 2;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dzcc01 = controlsPenaltyMPC_dz_cc + 2;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_rd01 = controlsPenaltyMPC_rd + 2;
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_Lbyrd01[2];
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_grad_cost01 = controlsPenaltyMPC_grad_cost + 2;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_grad_eq01 = controlsPenaltyMPC_grad_eq + 2;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_grad_ineq01 = controlsPenaltyMPC_grad_ineq + 2;
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_ctv01[2];
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_v01 = controlsPenaltyMPC_v + 2;
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_re01[2];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_beta01[2];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_betacc01[2];
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dvaff01 = controlsPenaltyMPC_dv_aff + 2;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dvcc01 = controlsPenaltyMPC_dv_cc + 2;
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_V01[4];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_Yd01[3];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_Ld01[3];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_yy01[2];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_bmy01[2];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_c01[2] = {0.0000000000000000E+000, 0.0000000000000000E+000};
int controlsPenaltyMPC_lbIdx01[2] = {0, 1};
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_llb01 = controlsPenaltyMPC_l + 4;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_slb01 = controlsPenaltyMPC_s + 4;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_llbbyslb01 = controlsPenaltyMPC_lbys + 4;
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_rilb01[2];
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dllbaff01 = controlsPenaltyMPC_dl_aff + 4;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dslbaff01 = controlsPenaltyMPC_ds_aff + 4;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dllbcc01 = controlsPenaltyMPC_dl_cc + 4;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dslbcc01 = controlsPenaltyMPC_ds_cc + 4;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_ccrhsl01 = controlsPenaltyMPC_ccrhs + 4;
int controlsPenaltyMPC_ubIdx01[2] = {0, 1};
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_lub01 = controlsPenaltyMPC_l + 6;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_sub01 = controlsPenaltyMPC_s + 6;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_lubbysub01 = controlsPenaltyMPC_lbys + 6;
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_riub01[2];
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dlubaff01 = controlsPenaltyMPC_dl_aff + 6;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dsubaff01 = controlsPenaltyMPC_ds_aff + 6;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dlubcc01 = controlsPenaltyMPC_dl_cc + 6;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dsubcc01 = controlsPenaltyMPC_ds_cc + 6;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_ccrhsub01 = controlsPenaltyMPC_ccrhs + 6;
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_Phi01[2];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_D01[2] = {0.0000000000000000E+000, 
0.0000000000000000E+000};
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_W01[2];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_Ysd01[4];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_Lsd01[4];
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_z02 = controlsPenaltyMPC_z + 4;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dzaff02 = controlsPenaltyMPC_dz_aff + 4;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dzcc02 = controlsPenaltyMPC_dz_cc + 4;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_rd02 = controlsPenaltyMPC_rd + 4;
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_Lbyrd02[2];
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_grad_cost02 = controlsPenaltyMPC_grad_cost + 4;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_grad_eq02 = controlsPenaltyMPC_grad_eq + 4;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_grad_ineq02 = controlsPenaltyMPC_grad_ineq + 4;
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_ctv02[2];
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_v02 = controlsPenaltyMPC_v + 4;
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_re02[2];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_beta02[2];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_betacc02[2];
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dvaff02 = controlsPenaltyMPC_dv_aff + 4;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dvcc02 = controlsPenaltyMPC_dv_cc + 4;
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_V02[4];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_Yd02[3];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_Ld02[3];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_yy02[2];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_bmy02[2];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_c02[2] = {0.0000000000000000E+000, 0.0000000000000000E+000};
int controlsPenaltyMPC_lbIdx02[2] = {0, 1};
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_llb02 = controlsPenaltyMPC_l + 8;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_slb02 = controlsPenaltyMPC_s + 8;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_llbbyslb02 = controlsPenaltyMPC_lbys + 8;
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_rilb02[2];
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dllbaff02 = controlsPenaltyMPC_dl_aff + 8;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dslbaff02 = controlsPenaltyMPC_ds_aff + 8;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dllbcc02 = controlsPenaltyMPC_dl_cc + 8;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dslbcc02 = controlsPenaltyMPC_ds_cc + 8;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_ccrhsl02 = controlsPenaltyMPC_ccrhs + 8;
int controlsPenaltyMPC_ubIdx02[2] = {0, 1};
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_lub02 = controlsPenaltyMPC_l + 10;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_sub02 = controlsPenaltyMPC_s + 10;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_lubbysub02 = controlsPenaltyMPC_lbys + 10;
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_riub02[2];
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dlubaff02 = controlsPenaltyMPC_dl_aff + 10;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dsubaff02 = controlsPenaltyMPC_ds_aff + 10;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dlubcc02 = controlsPenaltyMPC_dl_cc + 10;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dsubcc02 = controlsPenaltyMPC_ds_cc + 10;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_ccrhsub02 = controlsPenaltyMPC_ccrhs + 10;
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_Phi02[2];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_W02[2];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_Ysd02[4];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_Lsd02[4];
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_z03 = controlsPenaltyMPC_z + 6;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dzaff03 = controlsPenaltyMPC_dz_aff + 6;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dzcc03 = controlsPenaltyMPC_dz_cc + 6;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_rd03 = controlsPenaltyMPC_rd + 6;
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_Lbyrd03[2];
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_grad_cost03 = controlsPenaltyMPC_grad_cost + 6;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_grad_eq03 = controlsPenaltyMPC_grad_eq + 6;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_grad_ineq03 = controlsPenaltyMPC_grad_ineq + 6;
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_ctv03[2];
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_v03 = controlsPenaltyMPC_v + 6;
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_re03[2];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_beta03[2];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_betacc03[2];
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dvaff03 = controlsPenaltyMPC_dv_aff + 6;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dvcc03 = controlsPenaltyMPC_dv_cc + 6;
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_V03[4];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_Yd03[3];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_Ld03[3];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_yy03[2];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_bmy03[2];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_c03[2] = {0.0000000000000000E+000, 0.0000000000000000E+000};
int controlsPenaltyMPC_lbIdx03[2] = {0, 1};
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_llb03 = controlsPenaltyMPC_l + 12;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_slb03 = controlsPenaltyMPC_s + 12;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_llbbyslb03 = controlsPenaltyMPC_lbys + 12;
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_rilb03[2];
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dllbaff03 = controlsPenaltyMPC_dl_aff + 12;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dslbaff03 = controlsPenaltyMPC_ds_aff + 12;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dllbcc03 = controlsPenaltyMPC_dl_cc + 12;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dslbcc03 = controlsPenaltyMPC_ds_cc + 12;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_ccrhsl03 = controlsPenaltyMPC_ccrhs + 12;
int controlsPenaltyMPC_ubIdx03[2] = {0, 1};
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_lub03 = controlsPenaltyMPC_l + 14;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_sub03 = controlsPenaltyMPC_s + 14;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_lubbysub03 = controlsPenaltyMPC_lbys + 14;
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_riub03[2];
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dlubaff03 = controlsPenaltyMPC_dl_aff + 14;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dsubaff03 = controlsPenaltyMPC_ds_aff + 14;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dlubcc03 = controlsPenaltyMPC_dl_cc + 14;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dsubcc03 = controlsPenaltyMPC_ds_cc + 14;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_ccrhsub03 = controlsPenaltyMPC_ccrhs + 14;
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_Phi03[2];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_W03[2];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_Ysd03[4];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_Lsd03[4];
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_z04 = controlsPenaltyMPC_z + 8;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dzaff04 = controlsPenaltyMPC_dz_aff + 8;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dzcc04 = controlsPenaltyMPC_dz_cc + 8;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_rd04 = controlsPenaltyMPC_rd + 8;
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_Lbyrd04[2];
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_grad_cost04 = controlsPenaltyMPC_grad_cost + 8;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_grad_eq04 = controlsPenaltyMPC_grad_eq + 8;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_grad_ineq04 = controlsPenaltyMPC_grad_ineq + 8;
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_ctv04[2];
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_v04 = controlsPenaltyMPC_v + 8;
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_re04[2];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_beta04[2];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_betacc04[2];
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dvaff04 = controlsPenaltyMPC_dv_aff + 8;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dvcc04 = controlsPenaltyMPC_dv_cc + 8;
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_V04[4];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_Yd04[3];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_Ld04[3];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_yy04[2];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_bmy04[2];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_c04[2] = {0.0000000000000000E+000, 0.0000000000000000E+000};
int controlsPenaltyMPC_lbIdx04[2] = {0, 1};
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_llb04 = controlsPenaltyMPC_l + 16;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_slb04 = controlsPenaltyMPC_s + 16;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_llbbyslb04 = controlsPenaltyMPC_lbys + 16;
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_rilb04[2];
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dllbaff04 = controlsPenaltyMPC_dl_aff + 16;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dslbaff04 = controlsPenaltyMPC_ds_aff + 16;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dllbcc04 = controlsPenaltyMPC_dl_cc + 16;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dslbcc04 = controlsPenaltyMPC_ds_cc + 16;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_ccrhsl04 = controlsPenaltyMPC_ccrhs + 16;
int controlsPenaltyMPC_ubIdx04[2] = {0, 1};
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_lub04 = controlsPenaltyMPC_l + 18;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_sub04 = controlsPenaltyMPC_s + 18;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_lubbysub04 = controlsPenaltyMPC_lbys + 18;
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_riub04[2];
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dlubaff04 = controlsPenaltyMPC_dl_aff + 18;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dsubaff04 = controlsPenaltyMPC_ds_aff + 18;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dlubcc04 = controlsPenaltyMPC_dl_cc + 18;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dsubcc04 = controlsPenaltyMPC_ds_cc + 18;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_ccrhsub04 = controlsPenaltyMPC_ccrhs + 18;
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_Phi04[2];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_W04[2];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_Ysd04[4];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_Lsd04[4];
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_z05 = controlsPenaltyMPC_z + 10;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dzaff05 = controlsPenaltyMPC_dz_aff + 10;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dzcc05 = controlsPenaltyMPC_dz_cc + 10;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_rd05 = controlsPenaltyMPC_rd + 10;
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_Lbyrd05[2];
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_grad_cost05 = controlsPenaltyMPC_grad_cost + 10;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_grad_eq05 = controlsPenaltyMPC_grad_eq + 10;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_grad_ineq05 = controlsPenaltyMPC_grad_ineq + 10;
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_ctv05[2];
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_v05 = controlsPenaltyMPC_v + 10;
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_re05[2];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_beta05[2];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_betacc05[2];
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dvaff05 = controlsPenaltyMPC_dv_aff + 10;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dvcc05 = controlsPenaltyMPC_dv_cc + 10;
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_V05[4];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_Yd05[3];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_Ld05[3];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_yy05[2];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_bmy05[2];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_c05[2] = {0.0000000000000000E+000, 0.0000000000000000E+000};
int controlsPenaltyMPC_lbIdx05[2] = {0, 1};
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_llb05 = controlsPenaltyMPC_l + 20;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_slb05 = controlsPenaltyMPC_s + 20;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_llbbyslb05 = controlsPenaltyMPC_lbys + 20;
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_rilb05[2];
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dllbaff05 = controlsPenaltyMPC_dl_aff + 20;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dslbaff05 = controlsPenaltyMPC_ds_aff + 20;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dllbcc05 = controlsPenaltyMPC_dl_cc + 20;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dslbcc05 = controlsPenaltyMPC_ds_cc + 20;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_ccrhsl05 = controlsPenaltyMPC_ccrhs + 20;
int controlsPenaltyMPC_ubIdx05[2] = {0, 1};
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_lub05 = controlsPenaltyMPC_l + 22;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_sub05 = controlsPenaltyMPC_s + 22;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_lubbysub05 = controlsPenaltyMPC_lbys + 22;
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_riub05[2];
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dlubaff05 = controlsPenaltyMPC_dl_aff + 22;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dsubaff05 = controlsPenaltyMPC_ds_aff + 22;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dlubcc05 = controlsPenaltyMPC_dl_cc + 22;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dsubcc05 = controlsPenaltyMPC_ds_cc + 22;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_ccrhsub05 = controlsPenaltyMPC_ccrhs + 22;
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_Phi05[2];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_W05[2];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_Ysd05[4];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_Lsd05[4];
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_z06 = controlsPenaltyMPC_z + 12;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dzaff06 = controlsPenaltyMPC_dz_aff + 12;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dzcc06 = controlsPenaltyMPC_dz_cc + 12;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_rd06 = controlsPenaltyMPC_rd + 12;
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_Lbyrd06[2];
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_grad_cost06 = controlsPenaltyMPC_grad_cost + 12;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_grad_eq06 = controlsPenaltyMPC_grad_eq + 12;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_grad_ineq06 = controlsPenaltyMPC_grad_ineq + 12;
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_ctv06[2];
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_v06 = controlsPenaltyMPC_v + 12;
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_re06[2];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_beta06[2];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_betacc06[2];
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dvaff06 = controlsPenaltyMPC_dv_aff + 12;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dvcc06 = controlsPenaltyMPC_dv_cc + 12;
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_V06[4];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_Yd06[3];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_Ld06[3];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_yy06[2];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_bmy06[2];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_c06[2] = {0.0000000000000000E+000, 0.0000000000000000E+000};
int controlsPenaltyMPC_lbIdx06[2] = {0, 1};
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_llb06 = controlsPenaltyMPC_l + 24;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_slb06 = controlsPenaltyMPC_s + 24;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_llbbyslb06 = controlsPenaltyMPC_lbys + 24;
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_rilb06[2];
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dllbaff06 = controlsPenaltyMPC_dl_aff + 24;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dslbaff06 = controlsPenaltyMPC_ds_aff + 24;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dllbcc06 = controlsPenaltyMPC_dl_cc + 24;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dslbcc06 = controlsPenaltyMPC_ds_cc + 24;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_ccrhsl06 = controlsPenaltyMPC_ccrhs + 24;
int controlsPenaltyMPC_ubIdx06[2] = {0, 1};
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_lub06 = controlsPenaltyMPC_l + 26;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_sub06 = controlsPenaltyMPC_s + 26;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_lubbysub06 = controlsPenaltyMPC_lbys + 26;
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_riub06[2];
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dlubaff06 = controlsPenaltyMPC_dl_aff + 26;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dsubaff06 = controlsPenaltyMPC_ds_aff + 26;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dlubcc06 = controlsPenaltyMPC_dl_cc + 26;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dsubcc06 = controlsPenaltyMPC_ds_cc + 26;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_ccrhsub06 = controlsPenaltyMPC_ccrhs + 26;
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_Phi06[2];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_W06[2];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_Ysd06[4];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_Lsd06[4];
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_z07 = controlsPenaltyMPC_z + 14;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dzaff07 = controlsPenaltyMPC_dz_aff + 14;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dzcc07 = controlsPenaltyMPC_dz_cc + 14;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_rd07 = controlsPenaltyMPC_rd + 14;
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_Lbyrd07[2];
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_grad_cost07 = controlsPenaltyMPC_grad_cost + 14;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_grad_eq07 = controlsPenaltyMPC_grad_eq + 14;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_grad_ineq07 = controlsPenaltyMPC_grad_ineq + 14;
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_ctv07[2];
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_v07 = controlsPenaltyMPC_v + 14;
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_re07[2];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_beta07[2];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_betacc07[2];
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dvaff07 = controlsPenaltyMPC_dv_aff + 14;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dvcc07 = controlsPenaltyMPC_dv_cc + 14;
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_V07[4];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_Yd07[3];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_Ld07[3];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_yy07[2];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_bmy07[2];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_c07[2] = {0.0000000000000000E+000, 0.0000000000000000E+000};
int controlsPenaltyMPC_lbIdx07[2] = {0, 1};
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_llb07 = controlsPenaltyMPC_l + 28;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_slb07 = controlsPenaltyMPC_s + 28;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_llbbyslb07 = controlsPenaltyMPC_lbys + 28;
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_rilb07[2];
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dllbaff07 = controlsPenaltyMPC_dl_aff + 28;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dslbaff07 = controlsPenaltyMPC_ds_aff + 28;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dllbcc07 = controlsPenaltyMPC_dl_cc + 28;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dslbcc07 = controlsPenaltyMPC_ds_cc + 28;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_ccrhsl07 = controlsPenaltyMPC_ccrhs + 28;
int controlsPenaltyMPC_ubIdx07[2] = {0, 1};
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_lub07 = controlsPenaltyMPC_l + 30;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_sub07 = controlsPenaltyMPC_s + 30;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_lubbysub07 = controlsPenaltyMPC_lbys + 30;
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_riub07[2];
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dlubaff07 = controlsPenaltyMPC_dl_aff + 30;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dsubaff07 = controlsPenaltyMPC_ds_aff + 30;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dlubcc07 = controlsPenaltyMPC_dl_cc + 30;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dsubcc07 = controlsPenaltyMPC_ds_cc + 30;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_ccrhsub07 = controlsPenaltyMPC_ccrhs + 30;
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_Phi07[2];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_W07[2];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_Ysd07[4];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_Lsd07[4];
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_z08 = controlsPenaltyMPC_z + 16;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dzaff08 = controlsPenaltyMPC_dz_aff + 16;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dzcc08 = controlsPenaltyMPC_dz_cc + 16;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_rd08 = controlsPenaltyMPC_rd + 16;
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_Lbyrd08[2];
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_grad_cost08 = controlsPenaltyMPC_grad_cost + 16;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_grad_eq08 = controlsPenaltyMPC_grad_eq + 16;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_grad_ineq08 = controlsPenaltyMPC_grad_ineq + 16;
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_ctv08[2];
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_v08 = controlsPenaltyMPC_v + 16;
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_re08[2];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_beta08[2];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_betacc08[2];
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dvaff08 = controlsPenaltyMPC_dv_aff + 16;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dvcc08 = controlsPenaltyMPC_dv_cc + 16;
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_V08[4];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_Yd08[3];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_Ld08[3];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_yy08[2];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_bmy08[2];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_c08[2] = {0.0000000000000000E+000, 0.0000000000000000E+000};
int controlsPenaltyMPC_lbIdx08[2] = {0, 1};
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_llb08 = controlsPenaltyMPC_l + 32;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_slb08 = controlsPenaltyMPC_s + 32;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_llbbyslb08 = controlsPenaltyMPC_lbys + 32;
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_rilb08[2];
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dllbaff08 = controlsPenaltyMPC_dl_aff + 32;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dslbaff08 = controlsPenaltyMPC_ds_aff + 32;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dllbcc08 = controlsPenaltyMPC_dl_cc + 32;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dslbcc08 = controlsPenaltyMPC_ds_cc + 32;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_ccrhsl08 = controlsPenaltyMPC_ccrhs + 32;
int controlsPenaltyMPC_ubIdx08[2] = {0, 1};
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_lub08 = controlsPenaltyMPC_l + 34;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_sub08 = controlsPenaltyMPC_s + 34;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_lubbysub08 = controlsPenaltyMPC_lbys + 34;
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_riub08[2];
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dlubaff08 = controlsPenaltyMPC_dl_aff + 34;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dsubaff08 = controlsPenaltyMPC_ds_aff + 34;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dlubcc08 = controlsPenaltyMPC_dl_cc + 34;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dsubcc08 = controlsPenaltyMPC_ds_cc + 34;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_ccrhsub08 = controlsPenaltyMPC_ccrhs + 34;
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_Phi08[2];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_W08[2];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_Ysd08[4];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_Lsd08[4];
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_z09 = controlsPenaltyMPC_z + 18;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dzaff09 = controlsPenaltyMPC_dz_aff + 18;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dzcc09 = controlsPenaltyMPC_dz_cc + 18;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_rd09 = controlsPenaltyMPC_rd + 18;
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_Lbyrd09[2];
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_grad_cost09 = controlsPenaltyMPC_grad_cost + 18;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_grad_eq09 = controlsPenaltyMPC_grad_eq + 18;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_grad_ineq09 = controlsPenaltyMPC_grad_ineq + 18;
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_ctv09[2];
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_v09 = controlsPenaltyMPC_v + 18;
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_re09[2];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_beta09[2];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_betacc09[2];
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dvaff09 = controlsPenaltyMPC_dv_aff + 18;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dvcc09 = controlsPenaltyMPC_dv_cc + 18;
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_V09[4];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_Yd09[3];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_Ld09[3];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_yy09[2];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_bmy09[2];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_c09[2] = {0.0000000000000000E+000, 0.0000000000000000E+000};
int controlsPenaltyMPC_lbIdx09[2] = {0, 1};
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_llb09 = controlsPenaltyMPC_l + 36;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_slb09 = controlsPenaltyMPC_s + 36;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_llbbyslb09 = controlsPenaltyMPC_lbys + 36;
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_rilb09[2];
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dllbaff09 = controlsPenaltyMPC_dl_aff + 36;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dslbaff09 = controlsPenaltyMPC_ds_aff + 36;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dllbcc09 = controlsPenaltyMPC_dl_cc + 36;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dslbcc09 = controlsPenaltyMPC_ds_cc + 36;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_ccrhsl09 = controlsPenaltyMPC_ccrhs + 36;
int controlsPenaltyMPC_ubIdx09[2] = {0, 1};
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_lub09 = controlsPenaltyMPC_l + 38;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_sub09 = controlsPenaltyMPC_s + 38;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_lubbysub09 = controlsPenaltyMPC_lbys + 38;
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_riub09[2];
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dlubaff09 = controlsPenaltyMPC_dl_aff + 38;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dsubaff09 = controlsPenaltyMPC_ds_aff + 38;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dlubcc09 = controlsPenaltyMPC_dl_cc + 38;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dsubcc09 = controlsPenaltyMPC_ds_cc + 38;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_ccrhsub09 = controlsPenaltyMPC_ccrhs + 38;
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_Phi09[2];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_W09[2];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_Ysd09[4];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_Lsd09[4];
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_z10 = controlsPenaltyMPC_z + 20;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dzaff10 = controlsPenaltyMPC_dz_aff + 20;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dzcc10 = controlsPenaltyMPC_dz_cc + 20;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_rd10 = controlsPenaltyMPC_rd + 20;
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_Lbyrd10[2];
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_grad_cost10 = controlsPenaltyMPC_grad_cost + 20;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_grad_eq10 = controlsPenaltyMPC_grad_eq + 20;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_grad_ineq10 = controlsPenaltyMPC_grad_ineq + 20;
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_ctv10[2];
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_v10 = controlsPenaltyMPC_v + 20;
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_re10[2];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_beta10[2];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_betacc10[2];
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dvaff10 = controlsPenaltyMPC_dv_aff + 20;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dvcc10 = controlsPenaltyMPC_dv_cc + 20;
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_V10[4];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_Yd10[3];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_Ld10[3];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_yy10[2];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_bmy10[2];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_c10[2] = {0.0000000000000000E+000, 0.0000000000000000E+000};
int controlsPenaltyMPC_lbIdx10[2] = {0, 1};
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_llb10 = controlsPenaltyMPC_l + 40;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_slb10 = controlsPenaltyMPC_s + 40;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_llbbyslb10 = controlsPenaltyMPC_lbys + 40;
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_rilb10[2];
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dllbaff10 = controlsPenaltyMPC_dl_aff + 40;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dslbaff10 = controlsPenaltyMPC_ds_aff + 40;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dllbcc10 = controlsPenaltyMPC_dl_cc + 40;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dslbcc10 = controlsPenaltyMPC_ds_cc + 40;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_ccrhsl10 = controlsPenaltyMPC_ccrhs + 40;
int controlsPenaltyMPC_ubIdx10[2] = {0, 1};
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_lub10 = controlsPenaltyMPC_l + 42;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_sub10 = controlsPenaltyMPC_s + 42;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_lubbysub10 = controlsPenaltyMPC_lbys + 42;
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_riub10[2];
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dlubaff10 = controlsPenaltyMPC_dl_aff + 42;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dsubaff10 = controlsPenaltyMPC_ds_aff + 42;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dlubcc10 = controlsPenaltyMPC_dl_cc + 42;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dsubcc10 = controlsPenaltyMPC_ds_cc + 42;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_ccrhsub10 = controlsPenaltyMPC_ccrhs + 42;
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_Phi10[2];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_W10[2];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_Ysd10[4];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_Lsd10[4];
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_z11 = controlsPenaltyMPC_z + 22;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dzaff11 = controlsPenaltyMPC_dz_aff + 22;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dzcc11 = controlsPenaltyMPC_dz_cc + 22;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_rd11 = controlsPenaltyMPC_rd + 22;
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_Lbyrd11[2];
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_grad_cost11 = controlsPenaltyMPC_grad_cost + 22;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_grad_eq11 = controlsPenaltyMPC_grad_eq + 22;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_grad_ineq11 = controlsPenaltyMPC_grad_ineq + 22;
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_ctv11[2];
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_v11 = controlsPenaltyMPC_v + 22;
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_re11[2];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_beta11[2];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_betacc11[2];
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dvaff11 = controlsPenaltyMPC_dv_aff + 22;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dvcc11 = controlsPenaltyMPC_dv_cc + 22;
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_V11[4];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_Yd11[3];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_Ld11[3];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_yy11[2];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_bmy11[2];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_c11[2] = {0.0000000000000000E+000, 0.0000000000000000E+000};
int controlsPenaltyMPC_lbIdx11[2] = {0, 1};
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_llb11 = controlsPenaltyMPC_l + 44;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_slb11 = controlsPenaltyMPC_s + 44;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_llbbyslb11 = controlsPenaltyMPC_lbys + 44;
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_rilb11[2];
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dllbaff11 = controlsPenaltyMPC_dl_aff + 44;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dslbaff11 = controlsPenaltyMPC_ds_aff + 44;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dllbcc11 = controlsPenaltyMPC_dl_cc + 44;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dslbcc11 = controlsPenaltyMPC_ds_cc + 44;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_ccrhsl11 = controlsPenaltyMPC_ccrhs + 44;
int controlsPenaltyMPC_ubIdx11[2] = {0, 1};
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_lub11 = controlsPenaltyMPC_l + 46;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_sub11 = controlsPenaltyMPC_s + 46;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_lubbysub11 = controlsPenaltyMPC_lbys + 46;
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_riub11[2];
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dlubaff11 = controlsPenaltyMPC_dl_aff + 46;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dsubaff11 = controlsPenaltyMPC_ds_aff + 46;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dlubcc11 = controlsPenaltyMPC_dl_cc + 46;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dsubcc11 = controlsPenaltyMPC_ds_cc + 46;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_ccrhsub11 = controlsPenaltyMPC_ccrhs + 46;
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_Phi11[2];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_W11[2];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_Ysd11[4];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_Lsd11[4];
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_z12 = controlsPenaltyMPC_z + 24;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dzaff12 = controlsPenaltyMPC_dz_aff + 24;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dzcc12 = controlsPenaltyMPC_dz_cc + 24;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_rd12 = controlsPenaltyMPC_rd + 24;
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_Lbyrd12[2];
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_grad_cost12 = controlsPenaltyMPC_grad_cost + 24;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_grad_eq12 = controlsPenaltyMPC_grad_eq + 24;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_grad_ineq12 = controlsPenaltyMPC_grad_ineq + 24;
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_ctv12[2];
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_v12 = controlsPenaltyMPC_v + 24;
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_re12[2];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_beta12[2];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_betacc12[2];
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dvaff12 = controlsPenaltyMPC_dv_aff + 24;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dvcc12 = controlsPenaltyMPC_dv_cc + 24;
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_V12[4];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_Yd12[3];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_Ld12[3];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_yy12[2];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_bmy12[2];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_c12[2] = {0.0000000000000000E+000, 0.0000000000000000E+000};
int controlsPenaltyMPC_lbIdx12[2] = {0, 1};
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_llb12 = controlsPenaltyMPC_l + 48;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_slb12 = controlsPenaltyMPC_s + 48;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_llbbyslb12 = controlsPenaltyMPC_lbys + 48;
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_rilb12[2];
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dllbaff12 = controlsPenaltyMPC_dl_aff + 48;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dslbaff12 = controlsPenaltyMPC_ds_aff + 48;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dllbcc12 = controlsPenaltyMPC_dl_cc + 48;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dslbcc12 = controlsPenaltyMPC_ds_cc + 48;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_ccrhsl12 = controlsPenaltyMPC_ccrhs + 48;
int controlsPenaltyMPC_ubIdx12[2] = {0, 1};
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_lub12 = controlsPenaltyMPC_l + 50;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_sub12 = controlsPenaltyMPC_s + 50;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_lubbysub12 = controlsPenaltyMPC_lbys + 50;
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_riub12[2];
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dlubaff12 = controlsPenaltyMPC_dl_aff + 50;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dsubaff12 = controlsPenaltyMPC_ds_aff + 50;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dlubcc12 = controlsPenaltyMPC_dl_cc + 50;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dsubcc12 = controlsPenaltyMPC_ds_cc + 50;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_ccrhsub12 = controlsPenaltyMPC_ccrhs + 50;
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_Phi12[2];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_W12[2];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_Ysd12[4];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_Lsd12[4];
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_z13 = controlsPenaltyMPC_z + 26;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dzaff13 = controlsPenaltyMPC_dz_aff + 26;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dzcc13 = controlsPenaltyMPC_dz_cc + 26;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_rd13 = controlsPenaltyMPC_rd + 26;
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_Lbyrd13[2];
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_grad_cost13 = controlsPenaltyMPC_grad_cost + 26;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_grad_eq13 = controlsPenaltyMPC_grad_eq + 26;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_grad_ineq13 = controlsPenaltyMPC_grad_ineq + 26;
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_ctv13[2];
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_v13 = controlsPenaltyMPC_v + 26;
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_re13[2];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_beta13[2];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_betacc13[2];
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dvaff13 = controlsPenaltyMPC_dv_aff + 26;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dvcc13 = controlsPenaltyMPC_dv_cc + 26;
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_V13[4];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_Yd13[3];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_Ld13[3];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_yy13[2];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_bmy13[2];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_c13[2] = {0.0000000000000000E+000, 0.0000000000000000E+000};
int controlsPenaltyMPC_lbIdx13[2] = {0, 1};
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_llb13 = controlsPenaltyMPC_l + 52;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_slb13 = controlsPenaltyMPC_s + 52;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_llbbyslb13 = controlsPenaltyMPC_lbys + 52;
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_rilb13[2];
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dllbaff13 = controlsPenaltyMPC_dl_aff + 52;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dslbaff13 = controlsPenaltyMPC_ds_aff + 52;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dllbcc13 = controlsPenaltyMPC_dl_cc + 52;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dslbcc13 = controlsPenaltyMPC_ds_cc + 52;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_ccrhsl13 = controlsPenaltyMPC_ccrhs + 52;
int controlsPenaltyMPC_ubIdx13[2] = {0, 1};
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_lub13 = controlsPenaltyMPC_l + 54;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_sub13 = controlsPenaltyMPC_s + 54;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_lubbysub13 = controlsPenaltyMPC_lbys + 54;
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_riub13[2];
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dlubaff13 = controlsPenaltyMPC_dl_aff + 54;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dsubaff13 = controlsPenaltyMPC_ds_aff + 54;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dlubcc13 = controlsPenaltyMPC_dl_cc + 54;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_dsubcc13 = controlsPenaltyMPC_ds_cc + 54;
controlsPenaltyMPC_FLOAT* controlsPenaltyMPC_ccrhsub13 = controlsPenaltyMPC_ccrhs + 54;
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_Phi13[2];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_W13[2];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_Ysd13[4];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_Lsd13[4];
controlsPenaltyMPC_FLOAT musigma;
controlsPenaltyMPC_FLOAT sigma_3rdroot;
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_Diag1_0[2];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_Diag2_0[2];
controlsPenaltyMPC_FLOAT controlsPenaltyMPC_L_0[1];




/* SOLVER CODE --------------------------------------------------------- */
int controlsPenaltyMPC_solve(controlsPenaltyMPC_params* params, controlsPenaltyMPC_output* output, controlsPenaltyMPC_info* info)
{	
int exitcode;

#if controlsPenaltyMPC_SET_TIMING == 1
	controlsPenaltyMPC_timer solvertimer;
	controlsPenaltyMPC_tic(&solvertimer);
#endif
/* FUNCTION CALLS INTO LA LIBRARY -------------------------------------- */
info->it = 0;
controlsPenaltyMPC_LA_INITIALIZEVECTOR_28(controlsPenaltyMPC_z, 0);
controlsPenaltyMPC_LA_INITIALIZEVECTOR_28(controlsPenaltyMPC_v, 1);
controlsPenaltyMPC_LA_INITIALIZEVECTOR_56(controlsPenaltyMPC_l, 1);
controlsPenaltyMPC_LA_INITIALIZEVECTOR_56(controlsPenaltyMPC_s, 1);
info->mu = 0;
controlsPenaltyMPC_LA_DOTACC_56(controlsPenaltyMPC_l, controlsPenaltyMPC_s, &info->mu);
info->mu /= 56;
while( 1 ){
info->pobj = 0;
controlsPenaltyMPC_LA_DIAG_QUADFCN_2(params->H1, params->f1, controlsPenaltyMPC_z00, controlsPenaltyMPC_grad_cost00, &info->pobj);
controlsPenaltyMPC_LA_DIAG_QUADFCN_2(params->H2, params->f2, controlsPenaltyMPC_z01, controlsPenaltyMPC_grad_cost01, &info->pobj);
controlsPenaltyMPC_LA_DIAG_QUADFCN_2(params->H3, params->f3, controlsPenaltyMPC_z02, controlsPenaltyMPC_grad_cost02, &info->pobj);
controlsPenaltyMPC_LA_DIAG_QUADFCN_2(params->H4, params->f4, controlsPenaltyMPC_z03, controlsPenaltyMPC_grad_cost03, &info->pobj);
controlsPenaltyMPC_LA_DIAG_QUADFCN_2(params->H5, params->f5, controlsPenaltyMPC_z04, controlsPenaltyMPC_grad_cost04, &info->pobj);
controlsPenaltyMPC_LA_DIAG_QUADFCN_2(params->H6, params->f6, controlsPenaltyMPC_z05, controlsPenaltyMPC_grad_cost05, &info->pobj);
controlsPenaltyMPC_LA_DIAG_QUADFCN_2(params->H7, params->f7, controlsPenaltyMPC_z06, controlsPenaltyMPC_grad_cost06, &info->pobj);
controlsPenaltyMPC_LA_DIAG_QUADFCN_2(params->H8, params->f8, controlsPenaltyMPC_z07, controlsPenaltyMPC_grad_cost07, &info->pobj);
controlsPenaltyMPC_LA_DIAG_QUADFCN_2(params->H9, params->f9, controlsPenaltyMPC_z08, controlsPenaltyMPC_grad_cost08, &info->pobj);
controlsPenaltyMPC_LA_DIAG_QUADFCN_2(params->H10, params->f10, controlsPenaltyMPC_z09, controlsPenaltyMPC_grad_cost09, &info->pobj);
controlsPenaltyMPC_LA_DIAG_QUADFCN_2(params->H11, params->f11, controlsPenaltyMPC_z10, controlsPenaltyMPC_grad_cost10, &info->pobj);
controlsPenaltyMPC_LA_DIAG_QUADFCN_2(params->H12, params->f12, controlsPenaltyMPC_z11, controlsPenaltyMPC_grad_cost11, &info->pobj);
controlsPenaltyMPC_LA_DIAG_QUADFCN_2(params->H13, params->f13, controlsPenaltyMPC_z12, controlsPenaltyMPC_grad_cost12, &info->pobj);
controlsPenaltyMPC_LA_DIAG_QUADFCN_2(params->H14, params->f14, controlsPenaltyMPC_z13, controlsPenaltyMPC_grad_cost13, &info->pobj);
info->res_eq = 0;
info->dgap = 0;
controlsPenaltyMPC_LA_DIAGZERO_MVMSUB6_2(controlsPenaltyMPC_D00, controlsPenaltyMPC_z00, controlsPenaltyMPC_c00, controlsPenaltyMPC_v00, controlsPenaltyMPC_re00, &info->dgap, &info->res_eq);
controlsPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_2_2(controlsPenaltyMPC_C00, controlsPenaltyMPC_z00, controlsPenaltyMPC_D01, controlsPenaltyMPC_z01, controlsPenaltyMPC_c01, controlsPenaltyMPC_v01, controlsPenaltyMPC_re01, &info->dgap, &info->res_eq);
controlsPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_2_2(controlsPenaltyMPC_C00, controlsPenaltyMPC_z01, controlsPenaltyMPC_D01, controlsPenaltyMPC_z02, controlsPenaltyMPC_c02, controlsPenaltyMPC_v02, controlsPenaltyMPC_re02, &info->dgap, &info->res_eq);
controlsPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_2_2(controlsPenaltyMPC_C00, controlsPenaltyMPC_z02, controlsPenaltyMPC_D01, controlsPenaltyMPC_z03, controlsPenaltyMPC_c03, controlsPenaltyMPC_v03, controlsPenaltyMPC_re03, &info->dgap, &info->res_eq);
controlsPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_2_2(controlsPenaltyMPC_C00, controlsPenaltyMPC_z03, controlsPenaltyMPC_D01, controlsPenaltyMPC_z04, controlsPenaltyMPC_c04, controlsPenaltyMPC_v04, controlsPenaltyMPC_re04, &info->dgap, &info->res_eq);
controlsPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_2_2(controlsPenaltyMPC_C00, controlsPenaltyMPC_z04, controlsPenaltyMPC_D01, controlsPenaltyMPC_z05, controlsPenaltyMPC_c05, controlsPenaltyMPC_v05, controlsPenaltyMPC_re05, &info->dgap, &info->res_eq);
controlsPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_2_2(controlsPenaltyMPC_C00, controlsPenaltyMPC_z05, controlsPenaltyMPC_D01, controlsPenaltyMPC_z06, controlsPenaltyMPC_c06, controlsPenaltyMPC_v06, controlsPenaltyMPC_re06, &info->dgap, &info->res_eq);
controlsPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_2_2(controlsPenaltyMPC_C00, controlsPenaltyMPC_z06, controlsPenaltyMPC_D01, controlsPenaltyMPC_z07, controlsPenaltyMPC_c07, controlsPenaltyMPC_v07, controlsPenaltyMPC_re07, &info->dgap, &info->res_eq);
controlsPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_2_2(controlsPenaltyMPC_C00, controlsPenaltyMPC_z07, controlsPenaltyMPC_D01, controlsPenaltyMPC_z08, controlsPenaltyMPC_c08, controlsPenaltyMPC_v08, controlsPenaltyMPC_re08, &info->dgap, &info->res_eq);
controlsPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_2_2(controlsPenaltyMPC_C00, controlsPenaltyMPC_z08, controlsPenaltyMPC_D01, controlsPenaltyMPC_z09, controlsPenaltyMPC_c09, controlsPenaltyMPC_v09, controlsPenaltyMPC_re09, &info->dgap, &info->res_eq);
controlsPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_2_2(controlsPenaltyMPC_C00, controlsPenaltyMPC_z09, controlsPenaltyMPC_D01, controlsPenaltyMPC_z10, controlsPenaltyMPC_c10, controlsPenaltyMPC_v10, controlsPenaltyMPC_re10, &info->dgap, &info->res_eq);
controlsPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_2_2(controlsPenaltyMPC_C00, controlsPenaltyMPC_z10, controlsPenaltyMPC_D01, controlsPenaltyMPC_z11, controlsPenaltyMPC_c11, controlsPenaltyMPC_v11, controlsPenaltyMPC_re11, &info->dgap, &info->res_eq);
controlsPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_2_2(controlsPenaltyMPC_C00, controlsPenaltyMPC_z11, controlsPenaltyMPC_D01, controlsPenaltyMPC_z12, controlsPenaltyMPC_c12, controlsPenaltyMPC_v12, controlsPenaltyMPC_re12, &info->dgap, &info->res_eq);
controlsPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_2_2(controlsPenaltyMPC_C00, controlsPenaltyMPC_z12, controlsPenaltyMPC_D01, controlsPenaltyMPC_z13, controlsPenaltyMPC_c13, controlsPenaltyMPC_v13, controlsPenaltyMPC_re13, &info->dgap, &info->res_eq);
controlsPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlsPenaltyMPC_C00, controlsPenaltyMPC_v01, controlsPenaltyMPC_D00, controlsPenaltyMPC_v00, controlsPenaltyMPC_grad_eq00);
controlsPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlsPenaltyMPC_C00, controlsPenaltyMPC_v02, controlsPenaltyMPC_D01, controlsPenaltyMPC_v01, controlsPenaltyMPC_grad_eq01);
controlsPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlsPenaltyMPC_C00, controlsPenaltyMPC_v03, controlsPenaltyMPC_D01, controlsPenaltyMPC_v02, controlsPenaltyMPC_grad_eq02);
controlsPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlsPenaltyMPC_C00, controlsPenaltyMPC_v04, controlsPenaltyMPC_D01, controlsPenaltyMPC_v03, controlsPenaltyMPC_grad_eq03);
controlsPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlsPenaltyMPC_C00, controlsPenaltyMPC_v05, controlsPenaltyMPC_D01, controlsPenaltyMPC_v04, controlsPenaltyMPC_grad_eq04);
controlsPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlsPenaltyMPC_C00, controlsPenaltyMPC_v06, controlsPenaltyMPC_D01, controlsPenaltyMPC_v05, controlsPenaltyMPC_grad_eq05);
controlsPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlsPenaltyMPC_C00, controlsPenaltyMPC_v07, controlsPenaltyMPC_D01, controlsPenaltyMPC_v06, controlsPenaltyMPC_grad_eq06);
controlsPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlsPenaltyMPC_C00, controlsPenaltyMPC_v08, controlsPenaltyMPC_D01, controlsPenaltyMPC_v07, controlsPenaltyMPC_grad_eq07);
controlsPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlsPenaltyMPC_C00, controlsPenaltyMPC_v09, controlsPenaltyMPC_D01, controlsPenaltyMPC_v08, controlsPenaltyMPC_grad_eq08);
controlsPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlsPenaltyMPC_C00, controlsPenaltyMPC_v10, controlsPenaltyMPC_D01, controlsPenaltyMPC_v09, controlsPenaltyMPC_grad_eq09);
controlsPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlsPenaltyMPC_C00, controlsPenaltyMPC_v11, controlsPenaltyMPC_D01, controlsPenaltyMPC_v10, controlsPenaltyMPC_grad_eq10);
controlsPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlsPenaltyMPC_C00, controlsPenaltyMPC_v12, controlsPenaltyMPC_D01, controlsPenaltyMPC_v11, controlsPenaltyMPC_grad_eq11);
controlsPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlsPenaltyMPC_C00, controlsPenaltyMPC_v13, controlsPenaltyMPC_D01, controlsPenaltyMPC_v12, controlsPenaltyMPC_grad_eq12);
controlsPenaltyMPC_LA_DIAGZERO_MTVM_2_2(controlsPenaltyMPC_D01, controlsPenaltyMPC_v13, controlsPenaltyMPC_grad_eq13);
info->res_ineq = 0;
controlsPenaltyMPC_LA_VSUBADD3_2(params->lb1, controlsPenaltyMPC_z00, controlsPenaltyMPC_lbIdx00, controlsPenaltyMPC_llb00, controlsPenaltyMPC_slb00, controlsPenaltyMPC_rilb00, &info->dgap, &info->res_ineq);
controlsPenaltyMPC_LA_VSUBADD2_2(controlsPenaltyMPC_z00, controlsPenaltyMPC_ubIdx00, params->ub1, controlsPenaltyMPC_lub00, controlsPenaltyMPC_sub00, controlsPenaltyMPC_riub00, &info->dgap, &info->res_ineq);
controlsPenaltyMPC_LA_VSUBADD3_2(params->lb2, controlsPenaltyMPC_z01, controlsPenaltyMPC_lbIdx01, controlsPenaltyMPC_llb01, controlsPenaltyMPC_slb01, controlsPenaltyMPC_rilb01, &info->dgap, &info->res_ineq);
controlsPenaltyMPC_LA_VSUBADD2_2(controlsPenaltyMPC_z01, controlsPenaltyMPC_ubIdx01, params->ub2, controlsPenaltyMPC_lub01, controlsPenaltyMPC_sub01, controlsPenaltyMPC_riub01, &info->dgap, &info->res_ineq);
controlsPenaltyMPC_LA_VSUBADD3_2(params->lb3, controlsPenaltyMPC_z02, controlsPenaltyMPC_lbIdx02, controlsPenaltyMPC_llb02, controlsPenaltyMPC_slb02, controlsPenaltyMPC_rilb02, &info->dgap, &info->res_ineq);
controlsPenaltyMPC_LA_VSUBADD2_2(controlsPenaltyMPC_z02, controlsPenaltyMPC_ubIdx02, params->ub3, controlsPenaltyMPC_lub02, controlsPenaltyMPC_sub02, controlsPenaltyMPC_riub02, &info->dgap, &info->res_ineq);
controlsPenaltyMPC_LA_VSUBADD3_2(params->lb4, controlsPenaltyMPC_z03, controlsPenaltyMPC_lbIdx03, controlsPenaltyMPC_llb03, controlsPenaltyMPC_slb03, controlsPenaltyMPC_rilb03, &info->dgap, &info->res_ineq);
controlsPenaltyMPC_LA_VSUBADD2_2(controlsPenaltyMPC_z03, controlsPenaltyMPC_ubIdx03, params->ub4, controlsPenaltyMPC_lub03, controlsPenaltyMPC_sub03, controlsPenaltyMPC_riub03, &info->dgap, &info->res_ineq);
controlsPenaltyMPC_LA_VSUBADD3_2(params->lb5, controlsPenaltyMPC_z04, controlsPenaltyMPC_lbIdx04, controlsPenaltyMPC_llb04, controlsPenaltyMPC_slb04, controlsPenaltyMPC_rilb04, &info->dgap, &info->res_ineq);
controlsPenaltyMPC_LA_VSUBADD2_2(controlsPenaltyMPC_z04, controlsPenaltyMPC_ubIdx04, params->ub5, controlsPenaltyMPC_lub04, controlsPenaltyMPC_sub04, controlsPenaltyMPC_riub04, &info->dgap, &info->res_ineq);
controlsPenaltyMPC_LA_VSUBADD3_2(params->lb6, controlsPenaltyMPC_z05, controlsPenaltyMPC_lbIdx05, controlsPenaltyMPC_llb05, controlsPenaltyMPC_slb05, controlsPenaltyMPC_rilb05, &info->dgap, &info->res_ineq);
controlsPenaltyMPC_LA_VSUBADD2_2(controlsPenaltyMPC_z05, controlsPenaltyMPC_ubIdx05, params->ub6, controlsPenaltyMPC_lub05, controlsPenaltyMPC_sub05, controlsPenaltyMPC_riub05, &info->dgap, &info->res_ineq);
controlsPenaltyMPC_LA_VSUBADD3_2(params->lb7, controlsPenaltyMPC_z06, controlsPenaltyMPC_lbIdx06, controlsPenaltyMPC_llb06, controlsPenaltyMPC_slb06, controlsPenaltyMPC_rilb06, &info->dgap, &info->res_ineq);
controlsPenaltyMPC_LA_VSUBADD2_2(controlsPenaltyMPC_z06, controlsPenaltyMPC_ubIdx06, params->ub7, controlsPenaltyMPC_lub06, controlsPenaltyMPC_sub06, controlsPenaltyMPC_riub06, &info->dgap, &info->res_ineq);
controlsPenaltyMPC_LA_VSUBADD3_2(params->lb8, controlsPenaltyMPC_z07, controlsPenaltyMPC_lbIdx07, controlsPenaltyMPC_llb07, controlsPenaltyMPC_slb07, controlsPenaltyMPC_rilb07, &info->dgap, &info->res_ineq);
controlsPenaltyMPC_LA_VSUBADD2_2(controlsPenaltyMPC_z07, controlsPenaltyMPC_ubIdx07, params->ub8, controlsPenaltyMPC_lub07, controlsPenaltyMPC_sub07, controlsPenaltyMPC_riub07, &info->dgap, &info->res_ineq);
controlsPenaltyMPC_LA_VSUBADD3_2(params->lb9, controlsPenaltyMPC_z08, controlsPenaltyMPC_lbIdx08, controlsPenaltyMPC_llb08, controlsPenaltyMPC_slb08, controlsPenaltyMPC_rilb08, &info->dgap, &info->res_ineq);
controlsPenaltyMPC_LA_VSUBADD2_2(controlsPenaltyMPC_z08, controlsPenaltyMPC_ubIdx08, params->ub9, controlsPenaltyMPC_lub08, controlsPenaltyMPC_sub08, controlsPenaltyMPC_riub08, &info->dgap, &info->res_ineq);
controlsPenaltyMPC_LA_VSUBADD3_2(params->lb10, controlsPenaltyMPC_z09, controlsPenaltyMPC_lbIdx09, controlsPenaltyMPC_llb09, controlsPenaltyMPC_slb09, controlsPenaltyMPC_rilb09, &info->dgap, &info->res_ineq);
controlsPenaltyMPC_LA_VSUBADD2_2(controlsPenaltyMPC_z09, controlsPenaltyMPC_ubIdx09, params->ub10, controlsPenaltyMPC_lub09, controlsPenaltyMPC_sub09, controlsPenaltyMPC_riub09, &info->dgap, &info->res_ineq);
controlsPenaltyMPC_LA_VSUBADD3_2(params->lb11, controlsPenaltyMPC_z10, controlsPenaltyMPC_lbIdx10, controlsPenaltyMPC_llb10, controlsPenaltyMPC_slb10, controlsPenaltyMPC_rilb10, &info->dgap, &info->res_ineq);
controlsPenaltyMPC_LA_VSUBADD2_2(controlsPenaltyMPC_z10, controlsPenaltyMPC_ubIdx10, params->ub11, controlsPenaltyMPC_lub10, controlsPenaltyMPC_sub10, controlsPenaltyMPC_riub10, &info->dgap, &info->res_ineq);
controlsPenaltyMPC_LA_VSUBADD3_2(params->lb12, controlsPenaltyMPC_z11, controlsPenaltyMPC_lbIdx11, controlsPenaltyMPC_llb11, controlsPenaltyMPC_slb11, controlsPenaltyMPC_rilb11, &info->dgap, &info->res_ineq);
controlsPenaltyMPC_LA_VSUBADD2_2(controlsPenaltyMPC_z11, controlsPenaltyMPC_ubIdx11, params->ub12, controlsPenaltyMPC_lub11, controlsPenaltyMPC_sub11, controlsPenaltyMPC_riub11, &info->dgap, &info->res_ineq);
controlsPenaltyMPC_LA_VSUBADD3_2(params->lb13, controlsPenaltyMPC_z12, controlsPenaltyMPC_lbIdx12, controlsPenaltyMPC_llb12, controlsPenaltyMPC_slb12, controlsPenaltyMPC_rilb12, &info->dgap, &info->res_ineq);
controlsPenaltyMPC_LA_VSUBADD2_2(controlsPenaltyMPC_z12, controlsPenaltyMPC_ubIdx12, params->ub13, controlsPenaltyMPC_lub12, controlsPenaltyMPC_sub12, controlsPenaltyMPC_riub12, &info->dgap, &info->res_ineq);
controlsPenaltyMPC_LA_VSUBADD3_2(params->lb14, controlsPenaltyMPC_z13, controlsPenaltyMPC_lbIdx13, controlsPenaltyMPC_llb13, controlsPenaltyMPC_slb13, controlsPenaltyMPC_rilb13, &info->dgap, &info->res_ineq);
controlsPenaltyMPC_LA_VSUBADD2_2(controlsPenaltyMPC_z13, controlsPenaltyMPC_ubIdx13, params->ub14, controlsPenaltyMPC_lub13, controlsPenaltyMPC_sub13, controlsPenaltyMPC_riub13, &info->dgap, &info->res_ineq);
controlsPenaltyMPC_LA_INEQ_B_GRAD_2_2_2(controlsPenaltyMPC_lub00, controlsPenaltyMPC_sub00, controlsPenaltyMPC_riub00, controlsPenaltyMPC_llb00, controlsPenaltyMPC_slb00, controlsPenaltyMPC_rilb00, controlsPenaltyMPC_lbIdx00, controlsPenaltyMPC_ubIdx00, controlsPenaltyMPC_grad_ineq00, controlsPenaltyMPC_lubbysub00, controlsPenaltyMPC_llbbyslb00);
controlsPenaltyMPC_LA_INEQ_B_GRAD_2_2_2(controlsPenaltyMPC_lub01, controlsPenaltyMPC_sub01, controlsPenaltyMPC_riub01, controlsPenaltyMPC_llb01, controlsPenaltyMPC_slb01, controlsPenaltyMPC_rilb01, controlsPenaltyMPC_lbIdx01, controlsPenaltyMPC_ubIdx01, controlsPenaltyMPC_grad_ineq01, controlsPenaltyMPC_lubbysub01, controlsPenaltyMPC_llbbyslb01);
controlsPenaltyMPC_LA_INEQ_B_GRAD_2_2_2(controlsPenaltyMPC_lub02, controlsPenaltyMPC_sub02, controlsPenaltyMPC_riub02, controlsPenaltyMPC_llb02, controlsPenaltyMPC_slb02, controlsPenaltyMPC_rilb02, controlsPenaltyMPC_lbIdx02, controlsPenaltyMPC_ubIdx02, controlsPenaltyMPC_grad_ineq02, controlsPenaltyMPC_lubbysub02, controlsPenaltyMPC_llbbyslb02);
controlsPenaltyMPC_LA_INEQ_B_GRAD_2_2_2(controlsPenaltyMPC_lub03, controlsPenaltyMPC_sub03, controlsPenaltyMPC_riub03, controlsPenaltyMPC_llb03, controlsPenaltyMPC_slb03, controlsPenaltyMPC_rilb03, controlsPenaltyMPC_lbIdx03, controlsPenaltyMPC_ubIdx03, controlsPenaltyMPC_grad_ineq03, controlsPenaltyMPC_lubbysub03, controlsPenaltyMPC_llbbyslb03);
controlsPenaltyMPC_LA_INEQ_B_GRAD_2_2_2(controlsPenaltyMPC_lub04, controlsPenaltyMPC_sub04, controlsPenaltyMPC_riub04, controlsPenaltyMPC_llb04, controlsPenaltyMPC_slb04, controlsPenaltyMPC_rilb04, controlsPenaltyMPC_lbIdx04, controlsPenaltyMPC_ubIdx04, controlsPenaltyMPC_grad_ineq04, controlsPenaltyMPC_lubbysub04, controlsPenaltyMPC_llbbyslb04);
controlsPenaltyMPC_LA_INEQ_B_GRAD_2_2_2(controlsPenaltyMPC_lub05, controlsPenaltyMPC_sub05, controlsPenaltyMPC_riub05, controlsPenaltyMPC_llb05, controlsPenaltyMPC_slb05, controlsPenaltyMPC_rilb05, controlsPenaltyMPC_lbIdx05, controlsPenaltyMPC_ubIdx05, controlsPenaltyMPC_grad_ineq05, controlsPenaltyMPC_lubbysub05, controlsPenaltyMPC_llbbyslb05);
controlsPenaltyMPC_LA_INEQ_B_GRAD_2_2_2(controlsPenaltyMPC_lub06, controlsPenaltyMPC_sub06, controlsPenaltyMPC_riub06, controlsPenaltyMPC_llb06, controlsPenaltyMPC_slb06, controlsPenaltyMPC_rilb06, controlsPenaltyMPC_lbIdx06, controlsPenaltyMPC_ubIdx06, controlsPenaltyMPC_grad_ineq06, controlsPenaltyMPC_lubbysub06, controlsPenaltyMPC_llbbyslb06);
controlsPenaltyMPC_LA_INEQ_B_GRAD_2_2_2(controlsPenaltyMPC_lub07, controlsPenaltyMPC_sub07, controlsPenaltyMPC_riub07, controlsPenaltyMPC_llb07, controlsPenaltyMPC_slb07, controlsPenaltyMPC_rilb07, controlsPenaltyMPC_lbIdx07, controlsPenaltyMPC_ubIdx07, controlsPenaltyMPC_grad_ineq07, controlsPenaltyMPC_lubbysub07, controlsPenaltyMPC_llbbyslb07);
controlsPenaltyMPC_LA_INEQ_B_GRAD_2_2_2(controlsPenaltyMPC_lub08, controlsPenaltyMPC_sub08, controlsPenaltyMPC_riub08, controlsPenaltyMPC_llb08, controlsPenaltyMPC_slb08, controlsPenaltyMPC_rilb08, controlsPenaltyMPC_lbIdx08, controlsPenaltyMPC_ubIdx08, controlsPenaltyMPC_grad_ineq08, controlsPenaltyMPC_lubbysub08, controlsPenaltyMPC_llbbyslb08);
controlsPenaltyMPC_LA_INEQ_B_GRAD_2_2_2(controlsPenaltyMPC_lub09, controlsPenaltyMPC_sub09, controlsPenaltyMPC_riub09, controlsPenaltyMPC_llb09, controlsPenaltyMPC_slb09, controlsPenaltyMPC_rilb09, controlsPenaltyMPC_lbIdx09, controlsPenaltyMPC_ubIdx09, controlsPenaltyMPC_grad_ineq09, controlsPenaltyMPC_lubbysub09, controlsPenaltyMPC_llbbyslb09);
controlsPenaltyMPC_LA_INEQ_B_GRAD_2_2_2(controlsPenaltyMPC_lub10, controlsPenaltyMPC_sub10, controlsPenaltyMPC_riub10, controlsPenaltyMPC_llb10, controlsPenaltyMPC_slb10, controlsPenaltyMPC_rilb10, controlsPenaltyMPC_lbIdx10, controlsPenaltyMPC_ubIdx10, controlsPenaltyMPC_grad_ineq10, controlsPenaltyMPC_lubbysub10, controlsPenaltyMPC_llbbyslb10);
controlsPenaltyMPC_LA_INEQ_B_GRAD_2_2_2(controlsPenaltyMPC_lub11, controlsPenaltyMPC_sub11, controlsPenaltyMPC_riub11, controlsPenaltyMPC_llb11, controlsPenaltyMPC_slb11, controlsPenaltyMPC_rilb11, controlsPenaltyMPC_lbIdx11, controlsPenaltyMPC_ubIdx11, controlsPenaltyMPC_grad_ineq11, controlsPenaltyMPC_lubbysub11, controlsPenaltyMPC_llbbyslb11);
controlsPenaltyMPC_LA_INEQ_B_GRAD_2_2_2(controlsPenaltyMPC_lub12, controlsPenaltyMPC_sub12, controlsPenaltyMPC_riub12, controlsPenaltyMPC_llb12, controlsPenaltyMPC_slb12, controlsPenaltyMPC_rilb12, controlsPenaltyMPC_lbIdx12, controlsPenaltyMPC_ubIdx12, controlsPenaltyMPC_grad_ineq12, controlsPenaltyMPC_lubbysub12, controlsPenaltyMPC_llbbyslb12);
controlsPenaltyMPC_LA_INEQ_B_GRAD_2_2_2(controlsPenaltyMPC_lub13, controlsPenaltyMPC_sub13, controlsPenaltyMPC_riub13, controlsPenaltyMPC_llb13, controlsPenaltyMPC_slb13, controlsPenaltyMPC_rilb13, controlsPenaltyMPC_lbIdx13, controlsPenaltyMPC_ubIdx13, controlsPenaltyMPC_grad_ineq13, controlsPenaltyMPC_lubbysub13, controlsPenaltyMPC_llbbyslb13);
info->dobj = info->pobj - info->dgap;
info->rdgap = info->pobj ? info->dgap / info->pobj : 1e6;
if( info->rdgap < 0 ) info->rdgap = -info->rdgap;
if( info->mu < controlsPenaltyMPC_SET_ACC_KKTCOMPL
    && (info->rdgap < controlsPenaltyMPC_SET_ACC_RDGAP || info->dgap < controlsPenaltyMPC_SET_ACC_KKTCOMPL)
    && info->res_eq < controlsPenaltyMPC_SET_ACC_RESEQ
    && info->res_ineq < controlsPenaltyMPC_SET_ACC_RESINEQ ){
exitcode = controlsPenaltyMPC_OPTIMAL; break; }
if( info->it == controlsPenaltyMPC_SET_MAXIT ){
exitcode = controlsPenaltyMPC_MAXITREACHED; break; }
controlsPenaltyMPC_LA_VVADD3_28(controlsPenaltyMPC_grad_cost, controlsPenaltyMPC_grad_eq, controlsPenaltyMPC_grad_ineq, controlsPenaltyMPC_rd);
controlsPenaltyMPC_LA_DIAG_CHOL_ONELOOP_LBUB_2_2_2(params->H1, controlsPenaltyMPC_llbbyslb00, controlsPenaltyMPC_lbIdx00, controlsPenaltyMPC_lubbysub00, controlsPenaltyMPC_ubIdx00, controlsPenaltyMPC_Phi00);
controlsPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_2_2(controlsPenaltyMPC_Phi00, controlsPenaltyMPC_C00, controlsPenaltyMPC_V00);
controlsPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_2(controlsPenaltyMPC_Phi00, controlsPenaltyMPC_D00, controlsPenaltyMPC_W00);
controlsPenaltyMPC_LA_DENSE_DIAGZERO_MMTM_2_2_2(controlsPenaltyMPC_W00, controlsPenaltyMPC_V00, controlsPenaltyMPC_Ysd01);
controlsPenaltyMPC_LA_DIAG_FORWARDSUB_2(controlsPenaltyMPC_Phi00, controlsPenaltyMPC_rd00, controlsPenaltyMPC_Lbyrd00);
controlsPenaltyMPC_LA_DIAG_CHOL_ONELOOP_LBUB_2_2_2(params->H2, controlsPenaltyMPC_llbbyslb01, controlsPenaltyMPC_lbIdx01, controlsPenaltyMPC_lubbysub01, controlsPenaltyMPC_ubIdx01, controlsPenaltyMPC_Phi01);
controlsPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_2_2(controlsPenaltyMPC_Phi01, controlsPenaltyMPC_C00, controlsPenaltyMPC_V01);
controlsPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_2(controlsPenaltyMPC_Phi01, controlsPenaltyMPC_D01, controlsPenaltyMPC_W01);
controlsPenaltyMPC_LA_DENSE_DIAGZERO_MMTM_2_2_2(controlsPenaltyMPC_W01, controlsPenaltyMPC_V01, controlsPenaltyMPC_Ysd02);
controlsPenaltyMPC_LA_DIAG_FORWARDSUB_2(controlsPenaltyMPC_Phi01, controlsPenaltyMPC_rd01, controlsPenaltyMPC_Lbyrd01);
controlsPenaltyMPC_LA_DIAG_CHOL_ONELOOP_LBUB_2_2_2(params->H3, controlsPenaltyMPC_llbbyslb02, controlsPenaltyMPC_lbIdx02, controlsPenaltyMPC_lubbysub02, controlsPenaltyMPC_ubIdx02, controlsPenaltyMPC_Phi02);
controlsPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_2_2(controlsPenaltyMPC_Phi02, controlsPenaltyMPC_C00, controlsPenaltyMPC_V02);
controlsPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_2(controlsPenaltyMPC_Phi02, controlsPenaltyMPC_D01, controlsPenaltyMPC_W02);
controlsPenaltyMPC_LA_DENSE_DIAGZERO_MMTM_2_2_2(controlsPenaltyMPC_W02, controlsPenaltyMPC_V02, controlsPenaltyMPC_Ysd03);
controlsPenaltyMPC_LA_DIAG_FORWARDSUB_2(controlsPenaltyMPC_Phi02, controlsPenaltyMPC_rd02, controlsPenaltyMPC_Lbyrd02);
controlsPenaltyMPC_LA_DIAG_CHOL_ONELOOP_LBUB_2_2_2(params->H4, controlsPenaltyMPC_llbbyslb03, controlsPenaltyMPC_lbIdx03, controlsPenaltyMPC_lubbysub03, controlsPenaltyMPC_ubIdx03, controlsPenaltyMPC_Phi03);
controlsPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_2_2(controlsPenaltyMPC_Phi03, controlsPenaltyMPC_C00, controlsPenaltyMPC_V03);
controlsPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_2(controlsPenaltyMPC_Phi03, controlsPenaltyMPC_D01, controlsPenaltyMPC_W03);
controlsPenaltyMPC_LA_DENSE_DIAGZERO_MMTM_2_2_2(controlsPenaltyMPC_W03, controlsPenaltyMPC_V03, controlsPenaltyMPC_Ysd04);
controlsPenaltyMPC_LA_DIAG_FORWARDSUB_2(controlsPenaltyMPC_Phi03, controlsPenaltyMPC_rd03, controlsPenaltyMPC_Lbyrd03);
controlsPenaltyMPC_LA_DIAG_CHOL_ONELOOP_LBUB_2_2_2(params->H5, controlsPenaltyMPC_llbbyslb04, controlsPenaltyMPC_lbIdx04, controlsPenaltyMPC_lubbysub04, controlsPenaltyMPC_ubIdx04, controlsPenaltyMPC_Phi04);
controlsPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_2_2(controlsPenaltyMPC_Phi04, controlsPenaltyMPC_C00, controlsPenaltyMPC_V04);
controlsPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_2(controlsPenaltyMPC_Phi04, controlsPenaltyMPC_D01, controlsPenaltyMPC_W04);
controlsPenaltyMPC_LA_DENSE_DIAGZERO_MMTM_2_2_2(controlsPenaltyMPC_W04, controlsPenaltyMPC_V04, controlsPenaltyMPC_Ysd05);
controlsPenaltyMPC_LA_DIAG_FORWARDSUB_2(controlsPenaltyMPC_Phi04, controlsPenaltyMPC_rd04, controlsPenaltyMPC_Lbyrd04);
controlsPenaltyMPC_LA_DIAG_CHOL_ONELOOP_LBUB_2_2_2(params->H6, controlsPenaltyMPC_llbbyslb05, controlsPenaltyMPC_lbIdx05, controlsPenaltyMPC_lubbysub05, controlsPenaltyMPC_ubIdx05, controlsPenaltyMPC_Phi05);
controlsPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_2_2(controlsPenaltyMPC_Phi05, controlsPenaltyMPC_C00, controlsPenaltyMPC_V05);
controlsPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_2(controlsPenaltyMPC_Phi05, controlsPenaltyMPC_D01, controlsPenaltyMPC_W05);
controlsPenaltyMPC_LA_DENSE_DIAGZERO_MMTM_2_2_2(controlsPenaltyMPC_W05, controlsPenaltyMPC_V05, controlsPenaltyMPC_Ysd06);
controlsPenaltyMPC_LA_DIAG_FORWARDSUB_2(controlsPenaltyMPC_Phi05, controlsPenaltyMPC_rd05, controlsPenaltyMPC_Lbyrd05);
controlsPenaltyMPC_LA_DIAG_CHOL_ONELOOP_LBUB_2_2_2(params->H7, controlsPenaltyMPC_llbbyslb06, controlsPenaltyMPC_lbIdx06, controlsPenaltyMPC_lubbysub06, controlsPenaltyMPC_ubIdx06, controlsPenaltyMPC_Phi06);
controlsPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_2_2(controlsPenaltyMPC_Phi06, controlsPenaltyMPC_C00, controlsPenaltyMPC_V06);
controlsPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_2(controlsPenaltyMPC_Phi06, controlsPenaltyMPC_D01, controlsPenaltyMPC_W06);
controlsPenaltyMPC_LA_DENSE_DIAGZERO_MMTM_2_2_2(controlsPenaltyMPC_W06, controlsPenaltyMPC_V06, controlsPenaltyMPC_Ysd07);
controlsPenaltyMPC_LA_DIAG_FORWARDSUB_2(controlsPenaltyMPC_Phi06, controlsPenaltyMPC_rd06, controlsPenaltyMPC_Lbyrd06);
controlsPenaltyMPC_LA_DIAG_CHOL_ONELOOP_LBUB_2_2_2(params->H8, controlsPenaltyMPC_llbbyslb07, controlsPenaltyMPC_lbIdx07, controlsPenaltyMPC_lubbysub07, controlsPenaltyMPC_ubIdx07, controlsPenaltyMPC_Phi07);
controlsPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_2_2(controlsPenaltyMPC_Phi07, controlsPenaltyMPC_C00, controlsPenaltyMPC_V07);
controlsPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_2(controlsPenaltyMPC_Phi07, controlsPenaltyMPC_D01, controlsPenaltyMPC_W07);
controlsPenaltyMPC_LA_DENSE_DIAGZERO_MMTM_2_2_2(controlsPenaltyMPC_W07, controlsPenaltyMPC_V07, controlsPenaltyMPC_Ysd08);
controlsPenaltyMPC_LA_DIAG_FORWARDSUB_2(controlsPenaltyMPC_Phi07, controlsPenaltyMPC_rd07, controlsPenaltyMPC_Lbyrd07);
controlsPenaltyMPC_LA_DIAG_CHOL_ONELOOP_LBUB_2_2_2(params->H9, controlsPenaltyMPC_llbbyslb08, controlsPenaltyMPC_lbIdx08, controlsPenaltyMPC_lubbysub08, controlsPenaltyMPC_ubIdx08, controlsPenaltyMPC_Phi08);
controlsPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_2_2(controlsPenaltyMPC_Phi08, controlsPenaltyMPC_C00, controlsPenaltyMPC_V08);
controlsPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_2(controlsPenaltyMPC_Phi08, controlsPenaltyMPC_D01, controlsPenaltyMPC_W08);
controlsPenaltyMPC_LA_DENSE_DIAGZERO_MMTM_2_2_2(controlsPenaltyMPC_W08, controlsPenaltyMPC_V08, controlsPenaltyMPC_Ysd09);
controlsPenaltyMPC_LA_DIAG_FORWARDSUB_2(controlsPenaltyMPC_Phi08, controlsPenaltyMPC_rd08, controlsPenaltyMPC_Lbyrd08);
controlsPenaltyMPC_LA_DIAG_CHOL_ONELOOP_LBUB_2_2_2(params->H10, controlsPenaltyMPC_llbbyslb09, controlsPenaltyMPC_lbIdx09, controlsPenaltyMPC_lubbysub09, controlsPenaltyMPC_ubIdx09, controlsPenaltyMPC_Phi09);
controlsPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_2_2(controlsPenaltyMPC_Phi09, controlsPenaltyMPC_C00, controlsPenaltyMPC_V09);
controlsPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_2(controlsPenaltyMPC_Phi09, controlsPenaltyMPC_D01, controlsPenaltyMPC_W09);
controlsPenaltyMPC_LA_DENSE_DIAGZERO_MMTM_2_2_2(controlsPenaltyMPC_W09, controlsPenaltyMPC_V09, controlsPenaltyMPC_Ysd10);
controlsPenaltyMPC_LA_DIAG_FORWARDSUB_2(controlsPenaltyMPC_Phi09, controlsPenaltyMPC_rd09, controlsPenaltyMPC_Lbyrd09);
controlsPenaltyMPC_LA_DIAG_CHOL_ONELOOP_LBUB_2_2_2(params->H11, controlsPenaltyMPC_llbbyslb10, controlsPenaltyMPC_lbIdx10, controlsPenaltyMPC_lubbysub10, controlsPenaltyMPC_ubIdx10, controlsPenaltyMPC_Phi10);
controlsPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_2_2(controlsPenaltyMPC_Phi10, controlsPenaltyMPC_C00, controlsPenaltyMPC_V10);
controlsPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_2(controlsPenaltyMPC_Phi10, controlsPenaltyMPC_D01, controlsPenaltyMPC_W10);
controlsPenaltyMPC_LA_DENSE_DIAGZERO_MMTM_2_2_2(controlsPenaltyMPC_W10, controlsPenaltyMPC_V10, controlsPenaltyMPC_Ysd11);
controlsPenaltyMPC_LA_DIAG_FORWARDSUB_2(controlsPenaltyMPC_Phi10, controlsPenaltyMPC_rd10, controlsPenaltyMPC_Lbyrd10);
controlsPenaltyMPC_LA_DIAG_CHOL_ONELOOP_LBUB_2_2_2(params->H12, controlsPenaltyMPC_llbbyslb11, controlsPenaltyMPC_lbIdx11, controlsPenaltyMPC_lubbysub11, controlsPenaltyMPC_ubIdx11, controlsPenaltyMPC_Phi11);
controlsPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_2_2(controlsPenaltyMPC_Phi11, controlsPenaltyMPC_C00, controlsPenaltyMPC_V11);
controlsPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_2(controlsPenaltyMPC_Phi11, controlsPenaltyMPC_D01, controlsPenaltyMPC_W11);
controlsPenaltyMPC_LA_DENSE_DIAGZERO_MMTM_2_2_2(controlsPenaltyMPC_W11, controlsPenaltyMPC_V11, controlsPenaltyMPC_Ysd12);
controlsPenaltyMPC_LA_DIAG_FORWARDSUB_2(controlsPenaltyMPC_Phi11, controlsPenaltyMPC_rd11, controlsPenaltyMPC_Lbyrd11);
controlsPenaltyMPC_LA_DIAG_CHOL_ONELOOP_LBUB_2_2_2(params->H13, controlsPenaltyMPC_llbbyslb12, controlsPenaltyMPC_lbIdx12, controlsPenaltyMPC_lubbysub12, controlsPenaltyMPC_ubIdx12, controlsPenaltyMPC_Phi12);
controlsPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_2_2(controlsPenaltyMPC_Phi12, controlsPenaltyMPC_C00, controlsPenaltyMPC_V12);
controlsPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_2(controlsPenaltyMPC_Phi12, controlsPenaltyMPC_D01, controlsPenaltyMPC_W12);
controlsPenaltyMPC_LA_DENSE_DIAGZERO_MMTM_2_2_2(controlsPenaltyMPC_W12, controlsPenaltyMPC_V12, controlsPenaltyMPC_Ysd13);
controlsPenaltyMPC_LA_DIAG_FORWARDSUB_2(controlsPenaltyMPC_Phi12, controlsPenaltyMPC_rd12, controlsPenaltyMPC_Lbyrd12);
controlsPenaltyMPC_LA_DIAG_CHOL_ONELOOP_LBUB_2_2_2(params->H14, controlsPenaltyMPC_llbbyslb13, controlsPenaltyMPC_lbIdx13, controlsPenaltyMPC_lubbysub13, controlsPenaltyMPC_ubIdx13, controlsPenaltyMPC_Phi13);
controlsPenaltyMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_2(controlsPenaltyMPC_Phi13, controlsPenaltyMPC_D01, controlsPenaltyMPC_W13);
controlsPenaltyMPC_LA_DIAG_FORWARDSUB_2(controlsPenaltyMPC_Phi13, controlsPenaltyMPC_rd13, controlsPenaltyMPC_Lbyrd13);
controlsPenaltyMPC_LA_DIAGZERO_MMT_2(controlsPenaltyMPC_W00, controlsPenaltyMPC_Yd00);
controlsPenaltyMPC_LA_DIAGZERO_MVMSUB7_2(controlsPenaltyMPC_W00, controlsPenaltyMPC_Lbyrd00, controlsPenaltyMPC_re00, controlsPenaltyMPC_beta00);
controlsPenaltyMPC_LA_DENSE_DIAGZERO_MMT2_2_2_2(controlsPenaltyMPC_V00, controlsPenaltyMPC_W01, controlsPenaltyMPC_Yd01);
controlsPenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_2_2(controlsPenaltyMPC_V00, controlsPenaltyMPC_Lbyrd00, controlsPenaltyMPC_W01, controlsPenaltyMPC_Lbyrd01, controlsPenaltyMPC_re01, controlsPenaltyMPC_beta01);
controlsPenaltyMPC_LA_DENSE_DIAGZERO_MMT2_2_2_2(controlsPenaltyMPC_V01, controlsPenaltyMPC_W02, controlsPenaltyMPC_Yd02);
controlsPenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_2_2(controlsPenaltyMPC_V01, controlsPenaltyMPC_Lbyrd01, controlsPenaltyMPC_W02, controlsPenaltyMPC_Lbyrd02, controlsPenaltyMPC_re02, controlsPenaltyMPC_beta02);
controlsPenaltyMPC_LA_DENSE_DIAGZERO_MMT2_2_2_2(controlsPenaltyMPC_V02, controlsPenaltyMPC_W03, controlsPenaltyMPC_Yd03);
controlsPenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_2_2(controlsPenaltyMPC_V02, controlsPenaltyMPC_Lbyrd02, controlsPenaltyMPC_W03, controlsPenaltyMPC_Lbyrd03, controlsPenaltyMPC_re03, controlsPenaltyMPC_beta03);
controlsPenaltyMPC_LA_DENSE_DIAGZERO_MMT2_2_2_2(controlsPenaltyMPC_V03, controlsPenaltyMPC_W04, controlsPenaltyMPC_Yd04);
controlsPenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_2_2(controlsPenaltyMPC_V03, controlsPenaltyMPC_Lbyrd03, controlsPenaltyMPC_W04, controlsPenaltyMPC_Lbyrd04, controlsPenaltyMPC_re04, controlsPenaltyMPC_beta04);
controlsPenaltyMPC_LA_DENSE_DIAGZERO_MMT2_2_2_2(controlsPenaltyMPC_V04, controlsPenaltyMPC_W05, controlsPenaltyMPC_Yd05);
controlsPenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_2_2(controlsPenaltyMPC_V04, controlsPenaltyMPC_Lbyrd04, controlsPenaltyMPC_W05, controlsPenaltyMPC_Lbyrd05, controlsPenaltyMPC_re05, controlsPenaltyMPC_beta05);
controlsPenaltyMPC_LA_DENSE_DIAGZERO_MMT2_2_2_2(controlsPenaltyMPC_V05, controlsPenaltyMPC_W06, controlsPenaltyMPC_Yd06);
controlsPenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_2_2(controlsPenaltyMPC_V05, controlsPenaltyMPC_Lbyrd05, controlsPenaltyMPC_W06, controlsPenaltyMPC_Lbyrd06, controlsPenaltyMPC_re06, controlsPenaltyMPC_beta06);
controlsPenaltyMPC_LA_DENSE_DIAGZERO_MMT2_2_2_2(controlsPenaltyMPC_V06, controlsPenaltyMPC_W07, controlsPenaltyMPC_Yd07);
controlsPenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_2_2(controlsPenaltyMPC_V06, controlsPenaltyMPC_Lbyrd06, controlsPenaltyMPC_W07, controlsPenaltyMPC_Lbyrd07, controlsPenaltyMPC_re07, controlsPenaltyMPC_beta07);
controlsPenaltyMPC_LA_DENSE_DIAGZERO_MMT2_2_2_2(controlsPenaltyMPC_V07, controlsPenaltyMPC_W08, controlsPenaltyMPC_Yd08);
controlsPenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_2_2(controlsPenaltyMPC_V07, controlsPenaltyMPC_Lbyrd07, controlsPenaltyMPC_W08, controlsPenaltyMPC_Lbyrd08, controlsPenaltyMPC_re08, controlsPenaltyMPC_beta08);
controlsPenaltyMPC_LA_DENSE_DIAGZERO_MMT2_2_2_2(controlsPenaltyMPC_V08, controlsPenaltyMPC_W09, controlsPenaltyMPC_Yd09);
controlsPenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_2_2(controlsPenaltyMPC_V08, controlsPenaltyMPC_Lbyrd08, controlsPenaltyMPC_W09, controlsPenaltyMPC_Lbyrd09, controlsPenaltyMPC_re09, controlsPenaltyMPC_beta09);
controlsPenaltyMPC_LA_DENSE_DIAGZERO_MMT2_2_2_2(controlsPenaltyMPC_V09, controlsPenaltyMPC_W10, controlsPenaltyMPC_Yd10);
controlsPenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_2_2(controlsPenaltyMPC_V09, controlsPenaltyMPC_Lbyrd09, controlsPenaltyMPC_W10, controlsPenaltyMPC_Lbyrd10, controlsPenaltyMPC_re10, controlsPenaltyMPC_beta10);
controlsPenaltyMPC_LA_DENSE_DIAGZERO_MMT2_2_2_2(controlsPenaltyMPC_V10, controlsPenaltyMPC_W11, controlsPenaltyMPC_Yd11);
controlsPenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_2_2(controlsPenaltyMPC_V10, controlsPenaltyMPC_Lbyrd10, controlsPenaltyMPC_W11, controlsPenaltyMPC_Lbyrd11, controlsPenaltyMPC_re11, controlsPenaltyMPC_beta11);
controlsPenaltyMPC_LA_DENSE_DIAGZERO_MMT2_2_2_2(controlsPenaltyMPC_V11, controlsPenaltyMPC_W12, controlsPenaltyMPC_Yd12);
controlsPenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_2_2(controlsPenaltyMPC_V11, controlsPenaltyMPC_Lbyrd11, controlsPenaltyMPC_W12, controlsPenaltyMPC_Lbyrd12, controlsPenaltyMPC_re12, controlsPenaltyMPC_beta12);
controlsPenaltyMPC_LA_DENSE_DIAGZERO_MMT2_2_2_2(controlsPenaltyMPC_V12, controlsPenaltyMPC_W13, controlsPenaltyMPC_Yd13);
controlsPenaltyMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_2_2(controlsPenaltyMPC_V12, controlsPenaltyMPC_Lbyrd12, controlsPenaltyMPC_W13, controlsPenaltyMPC_Lbyrd13, controlsPenaltyMPC_re13, controlsPenaltyMPC_beta13);
controlsPenaltyMPC_LA_DENSE_CHOL_2(controlsPenaltyMPC_Yd00, controlsPenaltyMPC_Ld00);
controlsPenaltyMPC_LA_DENSE_FORWARDSUB_2(controlsPenaltyMPC_Ld00, controlsPenaltyMPC_beta00, controlsPenaltyMPC_yy00);
controlsPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(controlsPenaltyMPC_Ld00, controlsPenaltyMPC_Ysd01, controlsPenaltyMPC_Lsd01);
controlsPenaltyMPC_LA_DENSE_MMTSUB_2_2(controlsPenaltyMPC_Lsd01, controlsPenaltyMPC_Yd01);
controlsPenaltyMPC_LA_DENSE_CHOL_2(controlsPenaltyMPC_Yd01, controlsPenaltyMPC_Ld01);
controlsPenaltyMPC_LA_DENSE_MVMSUB1_2_2(controlsPenaltyMPC_Lsd01, controlsPenaltyMPC_yy00, controlsPenaltyMPC_beta01, controlsPenaltyMPC_bmy01);
controlsPenaltyMPC_LA_DENSE_FORWARDSUB_2(controlsPenaltyMPC_Ld01, controlsPenaltyMPC_bmy01, controlsPenaltyMPC_yy01);
controlsPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(controlsPenaltyMPC_Ld01, controlsPenaltyMPC_Ysd02, controlsPenaltyMPC_Lsd02);
controlsPenaltyMPC_LA_DENSE_MMTSUB_2_2(controlsPenaltyMPC_Lsd02, controlsPenaltyMPC_Yd02);
controlsPenaltyMPC_LA_DENSE_CHOL_2(controlsPenaltyMPC_Yd02, controlsPenaltyMPC_Ld02);
controlsPenaltyMPC_LA_DENSE_MVMSUB1_2_2(controlsPenaltyMPC_Lsd02, controlsPenaltyMPC_yy01, controlsPenaltyMPC_beta02, controlsPenaltyMPC_bmy02);
controlsPenaltyMPC_LA_DENSE_FORWARDSUB_2(controlsPenaltyMPC_Ld02, controlsPenaltyMPC_bmy02, controlsPenaltyMPC_yy02);
controlsPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(controlsPenaltyMPC_Ld02, controlsPenaltyMPC_Ysd03, controlsPenaltyMPC_Lsd03);
controlsPenaltyMPC_LA_DENSE_MMTSUB_2_2(controlsPenaltyMPC_Lsd03, controlsPenaltyMPC_Yd03);
controlsPenaltyMPC_LA_DENSE_CHOL_2(controlsPenaltyMPC_Yd03, controlsPenaltyMPC_Ld03);
controlsPenaltyMPC_LA_DENSE_MVMSUB1_2_2(controlsPenaltyMPC_Lsd03, controlsPenaltyMPC_yy02, controlsPenaltyMPC_beta03, controlsPenaltyMPC_bmy03);
controlsPenaltyMPC_LA_DENSE_FORWARDSUB_2(controlsPenaltyMPC_Ld03, controlsPenaltyMPC_bmy03, controlsPenaltyMPC_yy03);
controlsPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(controlsPenaltyMPC_Ld03, controlsPenaltyMPC_Ysd04, controlsPenaltyMPC_Lsd04);
controlsPenaltyMPC_LA_DENSE_MMTSUB_2_2(controlsPenaltyMPC_Lsd04, controlsPenaltyMPC_Yd04);
controlsPenaltyMPC_LA_DENSE_CHOL_2(controlsPenaltyMPC_Yd04, controlsPenaltyMPC_Ld04);
controlsPenaltyMPC_LA_DENSE_MVMSUB1_2_2(controlsPenaltyMPC_Lsd04, controlsPenaltyMPC_yy03, controlsPenaltyMPC_beta04, controlsPenaltyMPC_bmy04);
controlsPenaltyMPC_LA_DENSE_FORWARDSUB_2(controlsPenaltyMPC_Ld04, controlsPenaltyMPC_bmy04, controlsPenaltyMPC_yy04);
controlsPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(controlsPenaltyMPC_Ld04, controlsPenaltyMPC_Ysd05, controlsPenaltyMPC_Lsd05);
controlsPenaltyMPC_LA_DENSE_MMTSUB_2_2(controlsPenaltyMPC_Lsd05, controlsPenaltyMPC_Yd05);
controlsPenaltyMPC_LA_DENSE_CHOL_2(controlsPenaltyMPC_Yd05, controlsPenaltyMPC_Ld05);
controlsPenaltyMPC_LA_DENSE_MVMSUB1_2_2(controlsPenaltyMPC_Lsd05, controlsPenaltyMPC_yy04, controlsPenaltyMPC_beta05, controlsPenaltyMPC_bmy05);
controlsPenaltyMPC_LA_DENSE_FORWARDSUB_2(controlsPenaltyMPC_Ld05, controlsPenaltyMPC_bmy05, controlsPenaltyMPC_yy05);
controlsPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(controlsPenaltyMPC_Ld05, controlsPenaltyMPC_Ysd06, controlsPenaltyMPC_Lsd06);
controlsPenaltyMPC_LA_DENSE_MMTSUB_2_2(controlsPenaltyMPC_Lsd06, controlsPenaltyMPC_Yd06);
controlsPenaltyMPC_LA_DENSE_CHOL_2(controlsPenaltyMPC_Yd06, controlsPenaltyMPC_Ld06);
controlsPenaltyMPC_LA_DENSE_MVMSUB1_2_2(controlsPenaltyMPC_Lsd06, controlsPenaltyMPC_yy05, controlsPenaltyMPC_beta06, controlsPenaltyMPC_bmy06);
controlsPenaltyMPC_LA_DENSE_FORWARDSUB_2(controlsPenaltyMPC_Ld06, controlsPenaltyMPC_bmy06, controlsPenaltyMPC_yy06);
controlsPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(controlsPenaltyMPC_Ld06, controlsPenaltyMPC_Ysd07, controlsPenaltyMPC_Lsd07);
controlsPenaltyMPC_LA_DENSE_MMTSUB_2_2(controlsPenaltyMPC_Lsd07, controlsPenaltyMPC_Yd07);
controlsPenaltyMPC_LA_DENSE_CHOL_2(controlsPenaltyMPC_Yd07, controlsPenaltyMPC_Ld07);
controlsPenaltyMPC_LA_DENSE_MVMSUB1_2_2(controlsPenaltyMPC_Lsd07, controlsPenaltyMPC_yy06, controlsPenaltyMPC_beta07, controlsPenaltyMPC_bmy07);
controlsPenaltyMPC_LA_DENSE_FORWARDSUB_2(controlsPenaltyMPC_Ld07, controlsPenaltyMPC_bmy07, controlsPenaltyMPC_yy07);
controlsPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(controlsPenaltyMPC_Ld07, controlsPenaltyMPC_Ysd08, controlsPenaltyMPC_Lsd08);
controlsPenaltyMPC_LA_DENSE_MMTSUB_2_2(controlsPenaltyMPC_Lsd08, controlsPenaltyMPC_Yd08);
controlsPenaltyMPC_LA_DENSE_CHOL_2(controlsPenaltyMPC_Yd08, controlsPenaltyMPC_Ld08);
controlsPenaltyMPC_LA_DENSE_MVMSUB1_2_2(controlsPenaltyMPC_Lsd08, controlsPenaltyMPC_yy07, controlsPenaltyMPC_beta08, controlsPenaltyMPC_bmy08);
controlsPenaltyMPC_LA_DENSE_FORWARDSUB_2(controlsPenaltyMPC_Ld08, controlsPenaltyMPC_bmy08, controlsPenaltyMPC_yy08);
controlsPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(controlsPenaltyMPC_Ld08, controlsPenaltyMPC_Ysd09, controlsPenaltyMPC_Lsd09);
controlsPenaltyMPC_LA_DENSE_MMTSUB_2_2(controlsPenaltyMPC_Lsd09, controlsPenaltyMPC_Yd09);
controlsPenaltyMPC_LA_DENSE_CHOL_2(controlsPenaltyMPC_Yd09, controlsPenaltyMPC_Ld09);
controlsPenaltyMPC_LA_DENSE_MVMSUB1_2_2(controlsPenaltyMPC_Lsd09, controlsPenaltyMPC_yy08, controlsPenaltyMPC_beta09, controlsPenaltyMPC_bmy09);
controlsPenaltyMPC_LA_DENSE_FORWARDSUB_2(controlsPenaltyMPC_Ld09, controlsPenaltyMPC_bmy09, controlsPenaltyMPC_yy09);
controlsPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(controlsPenaltyMPC_Ld09, controlsPenaltyMPC_Ysd10, controlsPenaltyMPC_Lsd10);
controlsPenaltyMPC_LA_DENSE_MMTSUB_2_2(controlsPenaltyMPC_Lsd10, controlsPenaltyMPC_Yd10);
controlsPenaltyMPC_LA_DENSE_CHOL_2(controlsPenaltyMPC_Yd10, controlsPenaltyMPC_Ld10);
controlsPenaltyMPC_LA_DENSE_MVMSUB1_2_2(controlsPenaltyMPC_Lsd10, controlsPenaltyMPC_yy09, controlsPenaltyMPC_beta10, controlsPenaltyMPC_bmy10);
controlsPenaltyMPC_LA_DENSE_FORWARDSUB_2(controlsPenaltyMPC_Ld10, controlsPenaltyMPC_bmy10, controlsPenaltyMPC_yy10);
controlsPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(controlsPenaltyMPC_Ld10, controlsPenaltyMPC_Ysd11, controlsPenaltyMPC_Lsd11);
controlsPenaltyMPC_LA_DENSE_MMTSUB_2_2(controlsPenaltyMPC_Lsd11, controlsPenaltyMPC_Yd11);
controlsPenaltyMPC_LA_DENSE_CHOL_2(controlsPenaltyMPC_Yd11, controlsPenaltyMPC_Ld11);
controlsPenaltyMPC_LA_DENSE_MVMSUB1_2_2(controlsPenaltyMPC_Lsd11, controlsPenaltyMPC_yy10, controlsPenaltyMPC_beta11, controlsPenaltyMPC_bmy11);
controlsPenaltyMPC_LA_DENSE_FORWARDSUB_2(controlsPenaltyMPC_Ld11, controlsPenaltyMPC_bmy11, controlsPenaltyMPC_yy11);
controlsPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(controlsPenaltyMPC_Ld11, controlsPenaltyMPC_Ysd12, controlsPenaltyMPC_Lsd12);
controlsPenaltyMPC_LA_DENSE_MMTSUB_2_2(controlsPenaltyMPC_Lsd12, controlsPenaltyMPC_Yd12);
controlsPenaltyMPC_LA_DENSE_CHOL_2(controlsPenaltyMPC_Yd12, controlsPenaltyMPC_Ld12);
controlsPenaltyMPC_LA_DENSE_MVMSUB1_2_2(controlsPenaltyMPC_Lsd12, controlsPenaltyMPC_yy11, controlsPenaltyMPC_beta12, controlsPenaltyMPC_bmy12);
controlsPenaltyMPC_LA_DENSE_FORWARDSUB_2(controlsPenaltyMPC_Ld12, controlsPenaltyMPC_bmy12, controlsPenaltyMPC_yy12);
controlsPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(controlsPenaltyMPC_Ld12, controlsPenaltyMPC_Ysd13, controlsPenaltyMPC_Lsd13);
controlsPenaltyMPC_LA_DENSE_MMTSUB_2_2(controlsPenaltyMPC_Lsd13, controlsPenaltyMPC_Yd13);
controlsPenaltyMPC_LA_DENSE_CHOL_2(controlsPenaltyMPC_Yd13, controlsPenaltyMPC_Ld13);
controlsPenaltyMPC_LA_DENSE_MVMSUB1_2_2(controlsPenaltyMPC_Lsd13, controlsPenaltyMPC_yy12, controlsPenaltyMPC_beta13, controlsPenaltyMPC_bmy13);
controlsPenaltyMPC_LA_DENSE_FORWARDSUB_2(controlsPenaltyMPC_Ld13, controlsPenaltyMPC_bmy13, controlsPenaltyMPC_yy13);
controlsPenaltyMPC_LA_DENSE_BACKWARDSUB_2(controlsPenaltyMPC_Ld13, controlsPenaltyMPC_yy13, controlsPenaltyMPC_dvaff13);
controlsPenaltyMPC_LA_DENSE_MTVMSUB_2_2(controlsPenaltyMPC_Lsd13, controlsPenaltyMPC_dvaff13, controlsPenaltyMPC_yy12, controlsPenaltyMPC_bmy12);
controlsPenaltyMPC_LA_DENSE_BACKWARDSUB_2(controlsPenaltyMPC_Ld12, controlsPenaltyMPC_bmy12, controlsPenaltyMPC_dvaff12);
controlsPenaltyMPC_LA_DENSE_MTVMSUB_2_2(controlsPenaltyMPC_Lsd12, controlsPenaltyMPC_dvaff12, controlsPenaltyMPC_yy11, controlsPenaltyMPC_bmy11);
controlsPenaltyMPC_LA_DENSE_BACKWARDSUB_2(controlsPenaltyMPC_Ld11, controlsPenaltyMPC_bmy11, controlsPenaltyMPC_dvaff11);
controlsPenaltyMPC_LA_DENSE_MTVMSUB_2_2(controlsPenaltyMPC_Lsd11, controlsPenaltyMPC_dvaff11, controlsPenaltyMPC_yy10, controlsPenaltyMPC_bmy10);
controlsPenaltyMPC_LA_DENSE_BACKWARDSUB_2(controlsPenaltyMPC_Ld10, controlsPenaltyMPC_bmy10, controlsPenaltyMPC_dvaff10);
controlsPenaltyMPC_LA_DENSE_MTVMSUB_2_2(controlsPenaltyMPC_Lsd10, controlsPenaltyMPC_dvaff10, controlsPenaltyMPC_yy09, controlsPenaltyMPC_bmy09);
controlsPenaltyMPC_LA_DENSE_BACKWARDSUB_2(controlsPenaltyMPC_Ld09, controlsPenaltyMPC_bmy09, controlsPenaltyMPC_dvaff09);
controlsPenaltyMPC_LA_DENSE_MTVMSUB_2_2(controlsPenaltyMPC_Lsd09, controlsPenaltyMPC_dvaff09, controlsPenaltyMPC_yy08, controlsPenaltyMPC_bmy08);
controlsPenaltyMPC_LA_DENSE_BACKWARDSUB_2(controlsPenaltyMPC_Ld08, controlsPenaltyMPC_bmy08, controlsPenaltyMPC_dvaff08);
controlsPenaltyMPC_LA_DENSE_MTVMSUB_2_2(controlsPenaltyMPC_Lsd08, controlsPenaltyMPC_dvaff08, controlsPenaltyMPC_yy07, controlsPenaltyMPC_bmy07);
controlsPenaltyMPC_LA_DENSE_BACKWARDSUB_2(controlsPenaltyMPC_Ld07, controlsPenaltyMPC_bmy07, controlsPenaltyMPC_dvaff07);
controlsPenaltyMPC_LA_DENSE_MTVMSUB_2_2(controlsPenaltyMPC_Lsd07, controlsPenaltyMPC_dvaff07, controlsPenaltyMPC_yy06, controlsPenaltyMPC_bmy06);
controlsPenaltyMPC_LA_DENSE_BACKWARDSUB_2(controlsPenaltyMPC_Ld06, controlsPenaltyMPC_bmy06, controlsPenaltyMPC_dvaff06);
controlsPenaltyMPC_LA_DENSE_MTVMSUB_2_2(controlsPenaltyMPC_Lsd06, controlsPenaltyMPC_dvaff06, controlsPenaltyMPC_yy05, controlsPenaltyMPC_bmy05);
controlsPenaltyMPC_LA_DENSE_BACKWARDSUB_2(controlsPenaltyMPC_Ld05, controlsPenaltyMPC_bmy05, controlsPenaltyMPC_dvaff05);
controlsPenaltyMPC_LA_DENSE_MTVMSUB_2_2(controlsPenaltyMPC_Lsd05, controlsPenaltyMPC_dvaff05, controlsPenaltyMPC_yy04, controlsPenaltyMPC_bmy04);
controlsPenaltyMPC_LA_DENSE_BACKWARDSUB_2(controlsPenaltyMPC_Ld04, controlsPenaltyMPC_bmy04, controlsPenaltyMPC_dvaff04);
controlsPenaltyMPC_LA_DENSE_MTVMSUB_2_2(controlsPenaltyMPC_Lsd04, controlsPenaltyMPC_dvaff04, controlsPenaltyMPC_yy03, controlsPenaltyMPC_bmy03);
controlsPenaltyMPC_LA_DENSE_BACKWARDSUB_2(controlsPenaltyMPC_Ld03, controlsPenaltyMPC_bmy03, controlsPenaltyMPC_dvaff03);
controlsPenaltyMPC_LA_DENSE_MTVMSUB_2_2(controlsPenaltyMPC_Lsd03, controlsPenaltyMPC_dvaff03, controlsPenaltyMPC_yy02, controlsPenaltyMPC_bmy02);
controlsPenaltyMPC_LA_DENSE_BACKWARDSUB_2(controlsPenaltyMPC_Ld02, controlsPenaltyMPC_bmy02, controlsPenaltyMPC_dvaff02);
controlsPenaltyMPC_LA_DENSE_MTVMSUB_2_2(controlsPenaltyMPC_Lsd02, controlsPenaltyMPC_dvaff02, controlsPenaltyMPC_yy01, controlsPenaltyMPC_bmy01);
controlsPenaltyMPC_LA_DENSE_BACKWARDSUB_2(controlsPenaltyMPC_Ld01, controlsPenaltyMPC_bmy01, controlsPenaltyMPC_dvaff01);
controlsPenaltyMPC_LA_DENSE_MTVMSUB_2_2(controlsPenaltyMPC_Lsd01, controlsPenaltyMPC_dvaff01, controlsPenaltyMPC_yy00, controlsPenaltyMPC_bmy00);
controlsPenaltyMPC_LA_DENSE_BACKWARDSUB_2(controlsPenaltyMPC_Ld00, controlsPenaltyMPC_bmy00, controlsPenaltyMPC_dvaff00);
controlsPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlsPenaltyMPC_C00, controlsPenaltyMPC_dvaff01, controlsPenaltyMPC_D00, controlsPenaltyMPC_dvaff00, controlsPenaltyMPC_grad_eq00);
controlsPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlsPenaltyMPC_C00, controlsPenaltyMPC_dvaff02, controlsPenaltyMPC_D01, controlsPenaltyMPC_dvaff01, controlsPenaltyMPC_grad_eq01);
controlsPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlsPenaltyMPC_C00, controlsPenaltyMPC_dvaff03, controlsPenaltyMPC_D01, controlsPenaltyMPC_dvaff02, controlsPenaltyMPC_grad_eq02);
controlsPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlsPenaltyMPC_C00, controlsPenaltyMPC_dvaff04, controlsPenaltyMPC_D01, controlsPenaltyMPC_dvaff03, controlsPenaltyMPC_grad_eq03);
controlsPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlsPenaltyMPC_C00, controlsPenaltyMPC_dvaff05, controlsPenaltyMPC_D01, controlsPenaltyMPC_dvaff04, controlsPenaltyMPC_grad_eq04);
controlsPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlsPenaltyMPC_C00, controlsPenaltyMPC_dvaff06, controlsPenaltyMPC_D01, controlsPenaltyMPC_dvaff05, controlsPenaltyMPC_grad_eq05);
controlsPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlsPenaltyMPC_C00, controlsPenaltyMPC_dvaff07, controlsPenaltyMPC_D01, controlsPenaltyMPC_dvaff06, controlsPenaltyMPC_grad_eq06);
controlsPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlsPenaltyMPC_C00, controlsPenaltyMPC_dvaff08, controlsPenaltyMPC_D01, controlsPenaltyMPC_dvaff07, controlsPenaltyMPC_grad_eq07);
controlsPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlsPenaltyMPC_C00, controlsPenaltyMPC_dvaff09, controlsPenaltyMPC_D01, controlsPenaltyMPC_dvaff08, controlsPenaltyMPC_grad_eq08);
controlsPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlsPenaltyMPC_C00, controlsPenaltyMPC_dvaff10, controlsPenaltyMPC_D01, controlsPenaltyMPC_dvaff09, controlsPenaltyMPC_grad_eq09);
controlsPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlsPenaltyMPC_C00, controlsPenaltyMPC_dvaff11, controlsPenaltyMPC_D01, controlsPenaltyMPC_dvaff10, controlsPenaltyMPC_grad_eq10);
controlsPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlsPenaltyMPC_C00, controlsPenaltyMPC_dvaff12, controlsPenaltyMPC_D01, controlsPenaltyMPC_dvaff11, controlsPenaltyMPC_grad_eq11);
controlsPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlsPenaltyMPC_C00, controlsPenaltyMPC_dvaff13, controlsPenaltyMPC_D01, controlsPenaltyMPC_dvaff12, controlsPenaltyMPC_grad_eq12);
controlsPenaltyMPC_LA_DIAGZERO_MTVM_2_2(controlsPenaltyMPC_D01, controlsPenaltyMPC_dvaff13, controlsPenaltyMPC_grad_eq13);
controlsPenaltyMPC_LA_VSUB2_28(controlsPenaltyMPC_rd, controlsPenaltyMPC_grad_eq, controlsPenaltyMPC_rd);
controlsPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_2(controlsPenaltyMPC_Phi00, controlsPenaltyMPC_rd00, controlsPenaltyMPC_dzaff00);
controlsPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_2(controlsPenaltyMPC_Phi01, controlsPenaltyMPC_rd01, controlsPenaltyMPC_dzaff01);
controlsPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_2(controlsPenaltyMPC_Phi02, controlsPenaltyMPC_rd02, controlsPenaltyMPC_dzaff02);
controlsPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_2(controlsPenaltyMPC_Phi03, controlsPenaltyMPC_rd03, controlsPenaltyMPC_dzaff03);
controlsPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_2(controlsPenaltyMPC_Phi04, controlsPenaltyMPC_rd04, controlsPenaltyMPC_dzaff04);
controlsPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_2(controlsPenaltyMPC_Phi05, controlsPenaltyMPC_rd05, controlsPenaltyMPC_dzaff05);
controlsPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_2(controlsPenaltyMPC_Phi06, controlsPenaltyMPC_rd06, controlsPenaltyMPC_dzaff06);
controlsPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_2(controlsPenaltyMPC_Phi07, controlsPenaltyMPC_rd07, controlsPenaltyMPC_dzaff07);
controlsPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_2(controlsPenaltyMPC_Phi08, controlsPenaltyMPC_rd08, controlsPenaltyMPC_dzaff08);
controlsPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_2(controlsPenaltyMPC_Phi09, controlsPenaltyMPC_rd09, controlsPenaltyMPC_dzaff09);
controlsPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_2(controlsPenaltyMPC_Phi10, controlsPenaltyMPC_rd10, controlsPenaltyMPC_dzaff10);
controlsPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_2(controlsPenaltyMPC_Phi11, controlsPenaltyMPC_rd11, controlsPenaltyMPC_dzaff11);
controlsPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_2(controlsPenaltyMPC_Phi12, controlsPenaltyMPC_rd12, controlsPenaltyMPC_dzaff12);
controlsPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_2(controlsPenaltyMPC_Phi13, controlsPenaltyMPC_rd13, controlsPenaltyMPC_dzaff13);
controlsPenaltyMPC_LA_VSUB_INDEXED_2(controlsPenaltyMPC_dzaff00, controlsPenaltyMPC_lbIdx00, controlsPenaltyMPC_rilb00, controlsPenaltyMPC_dslbaff00);
controlsPenaltyMPC_LA_VSUB3_2(controlsPenaltyMPC_llbbyslb00, controlsPenaltyMPC_dslbaff00, controlsPenaltyMPC_llb00, controlsPenaltyMPC_dllbaff00);
controlsPenaltyMPC_LA_VSUB2_INDEXED_2(controlsPenaltyMPC_riub00, controlsPenaltyMPC_dzaff00, controlsPenaltyMPC_ubIdx00, controlsPenaltyMPC_dsubaff00);
controlsPenaltyMPC_LA_VSUB3_2(controlsPenaltyMPC_lubbysub00, controlsPenaltyMPC_dsubaff00, controlsPenaltyMPC_lub00, controlsPenaltyMPC_dlubaff00);
controlsPenaltyMPC_LA_VSUB_INDEXED_2(controlsPenaltyMPC_dzaff01, controlsPenaltyMPC_lbIdx01, controlsPenaltyMPC_rilb01, controlsPenaltyMPC_dslbaff01);
controlsPenaltyMPC_LA_VSUB3_2(controlsPenaltyMPC_llbbyslb01, controlsPenaltyMPC_dslbaff01, controlsPenaltyMPC_llb01, controlsPenaltyMPC_dllbaff01);
controlsPenaltyMPC_LA_VSUB2_INDEXED_2(controlsPenaltyMPC_riub01, controlsPenaltyMPC_dzaff01, controlsPenaltyMPC_ubIdx01, controlsPenaltyMPC_dsubaff01);
controlsPenaltyMPC_LA_VSUB3_2(controlsPenaltyMPC_lubbysub01, controlsPenaltyMPC_dsubaff01, controlsPenaltyMPC_lub01, controlsPenaltyMPC_dlubaff01);
controlsPenaltyMPC_LA_VSUB_INDEXED_2(controlsPenaltyMPC_dzaff02, controlsPenaltyMPC_lbIdx02, controlsPenaltyMPC_rilb02, controlsPenaltyMPC_dslbaff02);
controlsPenaltyMPC_LA_VSUB3_2(controlsPenaltyMPC_llbbyslb02, controlsPenaltyMPC_dslbaff02, controlsPenaltyMPC_llb02, controlsPenaltyMPC_dllbaff02);
controlsPenaltyMPC_LA_VSUB2_INDEXED_2(controlsPenaltyMPC_riub02, controlsPenaltyMPC_dzaff02, controlsPenaltyMPC_ubIdx02, controlsPenaltyMPC_dsubaff02);
controlsPenaltyMPC_LA_VSUB3_2(controlsPenaltyMPC_lubbysub02, controlsPenaltyMPC_dsubaff02, controlsPenaltyMPC_lub02, controlsPenaltyMPC_dlubaff02);
controlsPenaltyMPC_LA_VSUB_INDEXED_2(controlsPenaltyMPC_dzaff03, controlsPenaltyMPC_lbIdx03, controlsPenaltyMPC_rilb03, controlsPenaltyMPC_dslbaff03);
controlsPenaltyMPC_LA_VSUB3_2(controlsPenaltyMPC_llbbyslb03, controlsPenaltyMPC_dslbaff03, controlsPenaltyMPC_llb03, controlsPenaltyMPC_dllbaff03);
controlsPenaltyMPC_LA_VSUB2_INDEXED_2(controlsPenaltyMPC_riub03, controlsPenaltyMPC_dzaff03, controlsPenaltyMPC_ubIdx03, controlsPenaltyMPC_dsubaff03);
controlsPenaltyMPC_LA_VSUB3_2(controlsPenaltyMPC_lubbysub03, controlsPenaltyMPC_dsubaff03, controlsPenaltyMPC_lub03, controlsPenaltyMPC_dlubaff03);
controlsPenaltyMPC_LA_VSUB_INDEXED_2(controlsPenaltyMPC_dzaff04, controlsPenaltyMPC_lbIdx04, controlsPenaltyMPC_rilb04, controlsPenaltyMPC_dslbaff04);
controlsPenaltyMPC_LA_VSUB3_2(controlsPenaltyMPC_llbbyslb04, controlsPenaltyMPC_dslbaff04, controlsPenaltyMPC_llb04, controlsPenaltyMPC_dllbaff04);
controlsPenaltyMPC_LA_VSUB2_INDEXED_2(controlsPenaltyMPC_riub04, controlsPenaltyMPC_dzaff04, controlsPenaltyMPC_ubIdx04, controlsPenaltyMPC_dsubaff04);
controlsPenaltyMPC_LA_VSUB3_2(controlsPenaltyMPC_lubbysub04, controlsPenaltyMPC_dsubaff04, controlsPenaltyMPC_lub04, controlsPenaltyMPC_dlubaff04);
controlsPenaltyMPC_LA_VSUB_INDEXED_2(controlsPenaltyMPC_dzaff05, controlsPenaltyMPC_lbIdx05, controlsPenaltyMPC_rilb05, controlsPenaltyMPC_dslbaff05);
controlsPenaltyMPC_LA_VSUB3_2(controlsPenaltyMPC_llbbyslb05, controlsPenaltyMPC_dslbaff05, controlsPenaltyMPC_llb05, controlsPenaltyMPC_dllbaff05);
controlsPenaltyMPC_LA_VSUB2_INDEXED_2(controlsPenaltyMPC_riub05, controlsPenaltyMPC_dzaff05, controlsPenaltyMPC_ubIdx05, controlsPenaltyMPC_dsubaff05);
controlsPenaltyMPC_LA_VSUB3_2(controlsPenaltyMPC_lubbysub05, controlsPenaltyMPC_dsubaff05, controlsPenaltyMPC_lub05, controlsPenaltyMPC_dlubaff05);
controlsPenaltyMPC_LA_VSUB_INDEXED_2(controlsPenaltyMPC_dzaff06, controlsPenaltyMPC_lbIdx06, controlsPenaltyMPC_rilb06, controlsPenaltyMPC_dslbaff06);
controlsPenaltyMPC_LA_VSUB3_2(controlsPenaltyMPC_llbbyslb06, controlsPenaltyMPC_dslbaff06, controlsPenaltyMPC_llb06, controlsPenaltyMPC_dllbaff06);
controlsPenaltyMPC_LA_VSUB2_INDEXED_2(controlsPenaltyMPC_riub06, controlsPenaltyMPC_dzaff06, controlsPenaltyMPC_ubIdx06, controlsPenaltyMPC_dsubaff06);
controlsPenaltyMPC_LA_VSUB3_2(controlsPenaltyMPC_lubbysub06, controlsPenaltyMPC_dsubaff06, controlsPenaltyMPC_lub06, controlsPenaltyMPC_dlubaff06);
controlsPenaltyMPC_LA_VSUB_INDEXED_2(controlsPenaltyMPC_dzaff07, controlsPenaltyMPC_lbIdx07, controlsPenaltyMPC_rilb07, controlsPenaltyMPC_dslbaff07);
controlsPenaltyMPC_LA_VSUB3_2(controlsPenaltyMPC_llbbyslb07, controlsPenaltyMPC_dslbaff07, controlsPenaltyMPC_llb07, controlsPenaltyMPC_dllbaff07);
controlsPenaltyMPC_LA_VSUB2_INDEXED_2(controlsPenaltyMPC_riub07, controlsPenaltyMPC_dzaff07, controlsPenaltyMPC_ubIdx07, controlsPenaltyMPC_dsubaff07);
controlsPenaltyMPC_LA_VSUB3_2(controlsPenaltyMPC_lubbysub07, controlsPenaltyMPC_dsubaff07, controlsPenaltyMPC_lub07, controlsPenaltyMPC_dlubaff07);
controlsPenaltyMPC_LA_VSUB_INDEXED_2(controlsPenaltyMPC_dzaff08, controlsPenaltyMPC_lbIdx08, controlsPenaltyMPC_rilb08, controlsPenaltyMPC_dslbaff08);
controlsPenaltyMPC_LA_VSUB3_2(controlsPenaltyMPC_llbbyslb08, controlsPenaltyMPC_dslbaff08, controlsPenaltyMPC_llb08, controlsPenaltyMPC_dllbaff08);
controlsPenaltyMPC_LA_VSUB2_INDEXED_2(controlsPenaltyMPC_riub08, controlsPenaltyMPC_dzaff08, controlsPenaltyMPC_ubIdx08, controlsPenaltyMPC_dsubaff08);
controlsPenaltyMPC_LA_VSUB3_2(controlsPenaltyMPC_lubbysub08, controlsPenaltyMPC_dsubaff08, controlsPenaltyMPC_lub08, controlsPenaltyMPC_dlubaff08);
controlsPenaltyMPC_LA_VSUB_INDEXED_2(controlsPenaltyMPC_dzaff09, controlsPenaltyMPC_lbIdx09, controlsPenaltyMPC_rilb09, controlsPenaltyMPC_dslbaff09);
controlsPenaltyMPC_LA_VSUB3_2(controlsPenaltyMPC_llbbyslb09, controlsPenaltyMPC_dslbaff09, controlsPenaltyMPC_llb09, controlsPenaltyMPC_dllbaff09);
controlsPenaltyMPC_LA_VSUB2_INDEXED_2(controlsPenaltyMPC_riub09, controlsPenaltyMPC_dzaff09, controlsPenaltyMPC_ubIdx09, controlsPenaltyMPC_dsubaff09);
controlsPenaltyMPC_LA_VSUB3_2(controlsPenaltyMPC_lubbysub09, controlsPenaltyMPC_dsubaff09, controlsPenaltyMPC_lub09, controlsPenaltyMPC_dlubaff09);
controlsPenaltyMPC_LA_VSUB_INDEXED_2(controlsPenaltyMPC_dzaff10, controlsPenaltyMPC_lbIdx10, controlsPenaltyMPC_rilb10, controlsPenaltyMPC_dslbaff10);
controlsPenaltyMPC_LA_VSUB3_2(controlsPenaltyMPC_llbbyslb10, controlsPenaltyMPC_dslbaff10, controlsPenaltyMPC_llb10, controlsPenaltyMPC_dllbaff10);
controlsPenaltyMPC_LA_VSUB2_INDEXED_2(controlsPenaltyMPC_riub10, controlsPenaltyMPC_dzaff10, controlsPenaltyMPC_ubIdx10, controlsPenaltyMPC_dsubaff10);
controlsPenaltyMPC_LA_VSUB3_2(controlsPenaltyMPC_lubbysub10, controlsPenaltyMPC_dsubaff10, controlsPenaltyMPC_lub10, controlsPenaltyMPC_dlubaff10);
controlsPenaltyMPC_LA_VSUB_INDEXED_2(controlsPenaltyMPC_dzaff11, controlsPenaltyMPC_lbIdx11, controlsPenaltyMPC_rilb11, controlsPenaltyMPC_dslbaff11);
controlsPenaltyMPC_LA_VSUB3_2(controlsPenaltyMPC_llbbyslb11, controlsPenaltyMPC_dslbaff11, controlsPenaltyMPC_llb11, controlsPenaltyMPC_dllbaff11);
controlsPenaltyMPC_LA_VSUB2_INDEXED_2(controlsPenaltyMPC_riub11, controlsPenaltyMPC_dzaff11, controlsPenaltyMPC_ubIdx11, controlsPenaltyMPC_dsubaff11);
controlsPenaltyMPC_LA_VSUB3_2(controlsPenaltyMPC_lubbysub11, controlsPenaltyMPC_dsubaff11, controlsPenaltyMPC_lub11, controlsPenaltyMPC_dlubaff11);
controlsPenaltyMPC_LA_VSUB_INDEXED_2(controlsPenaltyMPC_dzaff12, controlsPenaltyMPC_lbIdx12, controlsPenaltyMPC_rilb12, controlsPenaltyMPC_dslbaff12);
controlsPenaltyMPC_LA_VSUB3_2(controlsPenaltyMPC_llbbyslb12, controlsPenaltyMPC_dslbaff12, controlsPenaltyMPC_llb12, controlsPenaltyMPC_dllbaff12);
controlsPenaltyMPC_LA_VSUB2_INDEXED_2(controlsPenaltyMPC_riub12, controlsPenaltyMPC_dzaff12, controlsPenaltyMPC_ubIdx12, controlsPenaltyMPC_dsubaff12);
controlsPenaltyMPC_LA_VSUB3_2(controlsPenaltyMPC_lubbysub12, controlsPenaltyMPC_dsubaff12, controlsPenaltyMPC_lub12, controlsPenaltyMPC_dlubaff12);
controlsPenaltyMPC_LA_VSUB_INDEXED_2(controlsPenaltyMPC_dzaff13, controlsPenaltyMPC_lbIdx13, controlsPenaltyMPC_rilb13, controlsPenaltyMPC_dslbaff13);
controlsPenaltyMPC_LA_VSUB3_2(controlsPenaltyMPC_llbbyslb13, controlsPenaltyMPC_dslbaff13, controlsPenaltyMPC_llb13, controlsPenaltyMPC_dllbaff13);
controlsPenaltyMPC_LA_VSUB2_INDEXED_2(controlsPenaltyMPC_riub13, controlsPenaltyMPC_dzaff13, controlsPenaltyMPC_ubIdx13, controlsPenaltyMPC_dsubaff13);
controlsPenaltyMPC_LA_VSUB3_2(controlsPenaltyMPC_lubbysub13, controlsPenaltyMPC_dsubaff13, controlsPenaltyMPC_lub13, controlsPenaltyMPC_dlubaff13);
info->lsit_aff = controlsPenaltyMPC_LINESEARCH_BACKTRACKING_AFFINE(controlsPenaltyMPC_l, controlsPenaltyMPC_s, controlsPenaltyMPC_dl_aff, controlsPenaltyMPC_ds_aff, &info->step_aff, &info->mu_aff);
if( info->lsit_aff == controlsPenaltyMPC_NOPROGRESS ){
exitcode = controlsPenaltyMPC_NOPROGRESS; break;
}
sigma_3rdroot = info->mu_aff / info->mu;
info->sigma = sigma_3rdroot*sigma_3rdroot*sigma_3rdroot;
musigma = info->mu * info->sigma;
controlsPenaltyMPC_LA_VSUB5_56(controlsPenaltyMPC_ds_aff, controlsPenaltyMPC_dl_aff, musigma, controlsPenaltyMPC_ccrhs);
controlsPenaltyMPC_LA_VSUB6_INDEXED_2_2_2(controlsPenaltyMPC_ccrhsub00, controlsPenaltyMPC_sub00, controlsPenaltyMPC_ubIdx00, controlsPenaltyMPC_ccrhsl00, controlsPenaltyMPC_slb00, controlsPenaltyMPC_lbIdx00, controlsPenaltyMPC_rd00);
controlsPenaltyMPC_LA_VSUB6_INDEXED_2_2_2(controlsPenaltyMPC_ccrhsub01, controlsPenaltyMPC_sub01, controlsPenaltyMPC_ubIdx01, controlsPenaltyMPC_ccrhsl01, controlsPenaltyMPC_slb01, controlsPenaltyMPC_lbIdx01, controlsPenaltyMPC_rd01);
controlsPenaltyMPC_LA_DIAG_FORWARDSUB_2(controlsPenaltyMPC_Phi00, controlsPenaltyMPC_rd00, controlsPenaltyMPC_Lbyrd00);
controlsPenaltyMPC_LA_DIAG_FORWARDSUB_2(controlsPenaltyMPC_Phi01, controlsPenaltyMPC_rd01, controlsPenaltyMPC_Lbyrd01);
controlsPenaltyMPC_LA_DIAGZERO_MVM_2(controlsPenaltyMPC_W00, controlsPenaltyMPC_Lbyrd00, controlsPenaltyMPC_beta00);
controlsPenaltyMPC_LA_DENSE_FORWARDSUB_2(controlsPenaltyMPC_Ld00, controlsPenaltyMPC_beta00, controlsPenaltyMPC_yy00);
controlsPenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_2_2_2(controlsPenaltyMPC_V00, controlsPenaltyMPC_Lbyrd00, controlsPenaltyMPC_W01, controlsPenaltyMPC_Lbyrd01, controlsPenaltyMPC_beta01);
controlsPenaltyMPC_LA_DENSE_MVMSUB1_2_2(controlsPenaltyMPC_Lsd01, controlsPenaltyMPC_yy00, controlsPenaltyMPC_beta01, controlsPenaltyMPC_bmy01);
controlsPenaltyMPC_LA_DENSE_FORWARDSUB_2(controlsPenaltyMPC_Ld01, controlsPenaltyMPC_bmy01, controlsPenaltyMPC_yy01);
controlsPenaltyMPC_LA_VSUB6_INDEXED_2_2_2(controlsPenaltyMPC_ccrhsub02, controlsPenaltyMPC_sub02, controlsPenaltyMPC_ubIdx02, controlsPenaltyMPC_ccrhsl02, controlsPenaltyMPC_slb02, controlsPenaltyMPC_lbIdx02, controlsPenaltyMPC_rd02);
controlsPenaltyMPC_LA_DIAG_FORWARDSUB_2(controlsPenaltyMPC_Phi02, controlsPenaltyMPC_rd02, controlsPenaltyMPC_Lbyrd02);
controlsPenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_2_2_2(controlsPenaltyMPC_V01, controlsPenaltyMPC_Lbyrd01, controlsPenaltyMPC_W02, controlsPenaltyMPC_Lbyrd02, controlsPenaltyMPC_beta02);
controlsPenaltyMPC_LA_DENSE_MVMSUB1_2_2(controlsPenaltyMPC_Lsd02, controlsPenaltyMPC_yy01, controlsPenaltyMPC_beta02, controlsPenaltyMPC_bmy02);
controlsPenaltyMPC_LA_DENSE_FORWARDSUB_2(controlsPenaltyMPC_Ld02, controlsPenaltyMPC_bmy02, controlsPenaltyMPC_yy02);
controlsPenaltyMPC_LA_VSUB6_INDEXED_2_2_2(controlsPenaltyMPC_ccrhsub03, controlsPenaltyMPC_sub03, controlsPenaltyMPC_ubIdx03, controlsPenaltyMPC_ccrhsl03, controlsPenaltyMPC_slb03, controlsPenaltyMPC_lbIdx03, controlsPenaltyMPC_rd03);
controlsPenaltyMPC_LA_DIAG_FORWARDSUB_2(controlsPenaltyMPC_Phi03, controlsPenaltyMPC_rd03, controlsPenaltyMPC_Lbyrd03);
controlsPenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_2_2_2(controlsPenaltyMPC_V02, controlsPenaltyMPC_Lbyrd02, controlsPenaltyMPC_W03, controlsPenaltyMPC_Lbyrd03, controlsPenaltyMPC_beta03);
controlsPenaltyMPC_LA_DENSE_MVMSUB1_2_2(controlsPenaltyMPC_Lsd03, controlsPenaltyMPC_yy02, controlsPenaltyMPC_beta03, controlsPenaltyMPC_bmy03);
controlsPenaltyMPC_LA_DENSE_FORWARDSUB_2(controlsPenaltyMPC_Ld03, controlsPenaltyMPC_bmy03, controlsPenaltyMPC_yy03);
controlsPenaltyMPC_LA_VSUB6_INDEXED_2_2_2(controlsPenaltyMPC_ccrhsub04, controlsPenaltyMPC_sub04, controlsPenaltyMPC_ubIdx04, controlsPenaltyMPC_ccrhsl04, controlsPenaltyMPC_slb04, controlsPenaltyMPC_lbIdx04, controlsPenaltyMPC_rd04);
controlsPenaltyMPC_LA_DIAG_FORWARDSUB_2(controlsPenaltyMPC_Phi04, controlsPenaltyMPC_rd04, controlsPenaltyMPC_Lbyrd04);
controlsPenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_2_2_2(controlsPenaltyMPC_V03, controlsPenaltyMPC_Lbyrd03, controlsPenaltyMPC_W04, controlsPenaltyMPC_Lbyrd04, controlsPenaltyMPC_beta04);
controlsPenaltyMPC_LA_DENSE_MVMSUB1_2_2(controlsPenaltyMPC_Lsd04, controlsPenaltyMPC_yy03, controlsPenaltyMPC_beta04, controlsPenaltyMPC_bmy04);
controlsPenaltyMPC_LA_DENSE_FORWARDSUB_2(controlsPenaltyMPC_Ld04, controlsPenaltyMPC_bmy04, controlsPenaltyMPC_yy04);
controlsPenaltyMPC_LA_VSUB6_INDEXED_2_2_2(controlsPenaltyMPC_ccrhsub05, controlsPenaltyMPC_sub05, controlsPenaltyMPC_ubIdx05, controlsPenaltyMPC_ccrhsl05, controlsPenaltyMPC_slb05, controlsPenaltyMPC_lbIdx05, controlsPenaltyMPC_rd05);
controlsPenaltyMPC_LA_DIAG_FORWARDSUB_2(controlsPenaltyMPC_Phi05, controlsPenaltyMPC_rd05, controlsPenaltyMPC_Lbyrd05);
controlsPenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_2_2_2(controlsPenaltyMPC_V04, controlsPenaltyMPC_Lbyrd04, controlsPenaltyMPC_W05, controlsPenaltyMPC_Lbyrd05, controlsPenaltyMPC_beta05);
controlsPenaltyMPC_LA_DENSE_MVMSUB1_2_2(controlsPenaltyMPC_Lsd05, controlsPenaltyMPC_yy04, controlsPenaltyMPC_beta05, controlsPenaltyMPC_bmy05);
controlsPenaltyMPC_LA_DENSE_FORWARDSUB_2(controlsPenaltyMPC_Ld05, controlsPenaltyMPC_bmy05, controlsPenaltyMPC_yy05);
controlsPenaltyMPC_LA_VSUB6_INDEXED_2_2_2(controlsPenaltyMPC_ccrhsub06, controlsPenaltyMPC_sub06, controlsPenaltyMPC_ubIdx06, controlsPenaltyMPC_ccrhsl06, controlsPenaltyMPC_slb06, controlsPenaltyMPC_lbIdx06, controlsPenaltyMPC_rd06);
controlsPenaltyMPC_LA_DIAG_FORWARDSUB_2(controlsPenaltyMPC_Phi06, controlsPenaltyMPC_rd06, controlsPenaltyMPC_Lbyrd06);
controlsPenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_2_2_2(controlsPenaltyMPC_V05, controlsPenaltyMPC_Lbyrd05, controlsPenaltyMPC_W06, controlsPenaltyMPC_Lbyrd06, controlsPenaltyMPC_beta06);
controlsPenaltyMPC_LA_DENSE_MVMSUB1_2_2(controlsPenaltyMPC_Lsd06, controlsPenaltyMPC_yy05, controlsPenaltyMPC_beta06, controlsPenaltyMPC_bmy06);
controlsPenaltyMPC_LA_DENSE_FORWARDSUB_2(controlsPenaltyMPC_Ld06, controlsPenaltyMPC_bmy06, controlsPenaltyMPC_yy06);
controlsPenaltyMPC_LA_VSUB6_INDEXED_2_2_2(controlsPenaltyMPC_ccrhsub07, controlsPenaltyMPC_sub07, controlsPenaltyMPC_ubIdx07, controlsPenaltyMPC_ccrhsl07, controlsPenaltyMPC_slb07, controlsPenaltyMPC_lbIdx07, controlsPenaltyMPC_rd07);
controlsPenaltyMPC_LA_DIAG_FORWARDSUB_2(controlsPenaltyMPC_Phi07, controlsPenaltyMPC_rd07, controlsPenaltyMPC_Lbyrd07);
controlsPenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_2_2_2(controlsPenaltyMPC_V06, controlsPenaltyMPC_Lbyrd06, controlsPenaltyMPC_W07, controlsPenaltyMPC_Lbyrd07, controlsPenaltyMPC_beta07);
controlsPenaltyMPC_LA_DENSE_MVMSUB1_2_2(controlsPenaltyMPC_Lsd07, controlsPenaltyMPC_yy06, controlsPenaltyMPC_beta07, controlsPenaltyMPC_bmy07);
controlsPenaltyMPC_LA_DENSE_FORWARDSUB_2(controlsPenaltyMPC_Ld07, controlsPenaltyMPC_bmy07, controlsPenaltyMPC_yy07);
controlsPenaltyMPC_LA_VSUB6_INDEXED_2_2_2(controlsPenaltyMPC_ccrhsub08, controlsPenaltyMPC_sub08, controlsPenaltyMPC_ubIdx08, controlsPenaltyMPC_ccrhsl08, controlsPenaltyMPC_slb08, controlsPenaltyMPC_lbIdx08, controlsPenaltyMPC_rd08);
controlsPenaltyMPC_LA_DIAG_FORWARDSUB_2(controlsPenaltyMPC_Phi08, controlsPenaltyMPC_rd08, controlsPenaltyMPC_Lbyrd08);
controlsPenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_2_2_2(controlsPenaltyMPC_V07, controlsPenaltyMPC_Lbyrd07, controlsPenaltyMPC_W08, controlsPenaltyMPC_Lbyrd08, controlsPenaltyMPC_beta08);
controlsPenaltyMPC_LA_DENSE_MVMSUB1_2_2(controlsPenaltyMPC_Lsd08, controlsPenaltyMPC_yy07, controlsPenaltyMPC_beta08, controlsPenaltyMPC_bmy08);
controlsPenaltyMPC_LA_DENSE_FORWARDSUB_2(controlsPenaltyMPC_Ld08, controlsPenaltyMPC_bmy08, controlsPenaltyMPC_yy08);
controlsPenaltyMPC_LA_VSUB6_INDEXED_2_2_2(controlsPenaltyMPC_ccrhsub09, controlsPenaltyMPC_sub09, controlsPenaltyMPC_ubIdx09, controlsPenaltyMPC_ccrhsl09, controlsPenaltyMPC_slb09, controlsPenaltyMPC_lbIdx09, controlsPenaltyMPC_rd09);
controlsPenaltyMPC_LA_DIAG_FORWARDSUB_2(controlsPenaltyMPC_Phi09, controlsPenaltyMPC_rd09, controlsPenaltyMPC_Lbyrd09);
controlsPenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_2_2_2(controlsPenaltyMPC_V08, controlsPenaltyMPC_Lbyrd08, controlsPenaltyMPC_W09, controlsPenaltyMPC_Lbyrd09, controlsPenaltyMPC_beta09);
controlsPenaltyMPC_LA_DENSE_MVMSUB1_2_2(controlsPenaltyMPC_Lsd09, controlsPenaltyMPC_yy08, controlsPenaltyMPC_beta09, controlsPenaltyMPC_bmy09);
controlsPenaltyMPC_LA_DENSE_FORWARDSUB_2(controlsPenaltyMPC_Ld09, controlsPenaltyMPC_bmy09, controlsPenaltyMPC_yy09);
controlsPenaltyMPC_LA_VSUB6_INDEXED_2_2_2(controlsPenaltyMPC_ccrhsub10, controlsPenaltyMPC_sub10, controlsPenaltyMPC_ubIdx10, controlsPenaltyMPC_ccrhsl10, controlsPenaltyMPC_slb10, controlsPenaltyMPC_lbIdx10, controlsPenaltyMPC_rd10);
controlsPenaltyMPC_LA_DIAG_FORWARDSUB_2(controlsPenaltyMPC_Phi10, controlsPenaltyMPC_rd10, controlsPenaltyMPC_Lbyrd10);
controlsPenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_2_2_2(controlsPenaltyMPC_V09, controlsPenaltyMPC_Lbyrd09, controlsPenaltyMPC_W10, controlsPenaltyMPC_Lbyrd10, controlsPenaltyMPC_beta10);
controlsPenaltyMPC_LA_DENSE_MVMSUB1_2_2(controlsPenaltyMPC_Lsd10, controlsPenaltyMPC_yy09, controlsPenaltyMPC_beta10, controlsPenaltyMPC_bmy10);
controlsPenaltyMPC_LA_DENSE_FORWARDSUB_2(controlsPenaltyMPC_Ld10, controlsPenaltyMPC_bmy10, controlsPenaltyMPC_yy10);
controlsPenaltyMPC_LA_VSUB6_INDEXED_2_2_2(controlsPenaltyMPC_ccrhsub11, controlsPenaltyMPC_sub11, controlsPenaltyMPC_ubIdx11, controlsPenaltyMPC_ccrhsl11, controlsPenaltyMPC_slb11, controlsPenaltyMPC_lbIdx11, controlsPenaltyMPC_rd11);
controlsPenaltyMPC_LA_DIAG_FORWARDSUB_2(controlsPenaltyMPC_Phi11, controlsPenaltyMPC_rd11, controlsPenaltyMPC_Lbyrd11);
controlsPenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_2_2_2(controlsPenaltyMPC_V10, controlsPenaltyMPC_Lbyrd10, controlsPenaltyMPC_W11, controlsPenaltyMPC_Lbyrd11, controlsPenaltyMPC_beta11);
controlsPenaltyMPC_LA_DENSE_MVMSUB1_2_2(controlsPenaltyMPC_Lsd11, controlsPenaltyMPC_yy10, controlsPenaltyMPC_beta11, controlsPenaltyMPC_bmy11);
controlsPenaltyMPC_LA_DENSE_FORWARDSUB_2(controlsPenaltyMPC_Ld11, controlsPenaltyMPC_bmy11, controlsPenaltyMPC_yy11);
controlsPenaltyMPC_LA_VSUB6_INDEXED_2_2_2(controlsPenaltyMPC_ccrhsub12, controlsPenaltyMPC_sub12, controlsPenaltyMPC_ubIdx12, controlsPenaltyMPC_ccrhsl12, controlsPenaltyMPC_slb12, controlsPenaltyMPC_lbIdx12, controlsPenaltyMPC_rd12);
controlsPenaltyMPC_LA_DIAG_FORWARDSUB_2(controlsPenaltyMPC_Phi12, controlsPenaltyMPC_rd12, controlsPenaltyMPC_Lbyrd12);
controlsPenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_2_2_2(controlsPenaltyMPC_V11, controlsPenaltyMPC_Lbyrd11, controlsPenaltyMPC_W12, controlsPenaltyMPC_Lbyrd12, controlsPenaltyMPC_beta12);
controlsPenaltyMPC_LA_DENSE_MVMSUB1_2_2(controlsPenaltyMPC_Lsd12, controlsPenaltyMPC_yy11, controlsPenaltyMPC_beta12, controlsPenaltyMPC_bmy12);
controlsPenaltyMPC_LA_DENSE_FORWARDSUB_2(controlsPenaltyMPC_Ld12, controlsPenaltyMPC_bmy12, controlsPenaltyMPC_yy12);
controlsPenaltyMPC_LA_VSUB6_INDEXED_2_2_2(controlsPenaltyMPC_ccrhsub13, controlsPenaltyMPC_sub13, controlsPenaltyMPC_ubIdx13, controlsPenaltyMPC_ccrhsl13, controlsPenaltyMPC_slb13, controlsPenaltyMPC_lbIdx13, controlsPenaltyMPC_rd13);
controlsPenaltyMPC_LA_DIAG_FORWARDSUB_2(controlsPenaltyMPC_Phi13, controlsPenaltyMPC_rd13, controlsPenaltyMPC_Lbyrd13);
controlsPenaltyMPC_LA_DENSE_DIAGZERO_2MVMADD_2_2_2(controlsPenaltyMPC_V12, controlsPenaltyMPC_Lbyrd12, controlsPenaltyMPC_W13, controlsPenaltyMPC_Lbyrd13, controlsPenaltyMPC_beta13);
controlsPenaltyMPC_LA_DENSE_MVMSUB1_2_2(controlsPenaltyMPC_Lsd13, controlsPenaltyMPC_yy12, controlsPenaltyMPC_beta13, controlsPenaltyMPC_bmy13);
controlsPenaltyMPC_LA_DENSE_FORWARDSUB_2(controlsPenaltyMPC_Ld13, controlsPenaltyMPC_bmy13, controlsPenaltyMPC_yy13);
controlsPenaltyMPC_LA_DENSE_BACKWARDSUB_2(controlsPenaltyMPC_Ld13, controlsPenaltyMPC_yy13, controlsPenaltyMPC_dvcc13);
controlsPenaltyMPC_LA_DENSE_MTVMSUB_2_2(controlsPenaltyMPC_Lsd13, controlsPenaltyMPC_dvcc13, controlsPenaltyMPC_yy12, controlsPenaltyMPC_bmy12);
controlsPenaltyMPC_LA_DENSE_BACKWARDSUB_2(controlsPenaltyMPC_Ld12, controlsPenaltyMPC_bmy12, controlsPenaltyMPC_dvcc12);
controlsPenaltyMPC_LA_DENSE_MTVMSUB_2_2(controlsPenaltyMPC_Lsd12, controlsPenaltyMPC_dvcc12, controlsPenaltyMPC_yy11, controlsPenaltyMPC_bmy11);
controlsPenaltyMPC_LA_DENSE_BACKWARDSUB_2(controlsPenaltyMPC_Ld11, controlsPenaltyMPC_bmy11, controlsPenaltyMPC_dvcc11);
controlsPenaltyMPC_LA_DENSE_MTVMSUB_2_2(controlsPenaltyMPC_Lsd11, controlsPenaltyMPC_dvcc11, controlsPenaltyMPC_yy10, controlsPenaltyMPC_bmy10);
controlsPenaltyMPC_LA_DENSE_BACKWARDSUB_2(controlsPenaltyMPC_Ld10, controlsPenaltyMPC_bmy10, controlsPenaltyMPC_dvcc10);
controlsPenaltyMPC_LA_DENSE_MTVMSUB_2_2(controlsPenaltyMPC_Lsd10, controlsPenaltyMPC_dvcc10, controlsPenaltyMPC_yy09, controlsPenaltyMPC_bmy09);
controlsPenaltyMPC_LA_DENSE_BACKWARDSUB_2(controlsPenaltyMPC_Ld09, controlsPenaltyMPC_bmy09, controlsPenaltyMPC_dvcc09);
controlsPenaltyMPC_LA_DENSE_MTVMSUB_2_2(controlsPenaltyMPC_Lsd09, controlsPenaltyMPC_dvcc09, controlsPenaltyMPC_yy08, controlsPenaltyMPC_bmy08);
controlsPenaltyMPC_LA_DENSE_BACKWARDSUB_2(controlsPenaltyMPC_Ld08, controlsPenaltyMPC_bmy08, controlsPenaltyMPC_dvcc08);
controlsPenaltyMPC_LA_DENSE_MTVMSUB_2_2(controlsPenaltyMPC_Lsd08, controlsPenaltyMPC_dvcc08, controlsPenaltyMPC_yy07, controlsPenaltyMPC_bmy07);
controlsPenaltyMPC_LA_DENSE_BACKWARDSUB_2(controlsPenaltyMPC_Ld07, controlsPenaltyMPC_bmy07, controlsPenaltyMPC_dvcc07);
controlsPenaltyMPC_LA_DENSE_MTVMSUB_2_2(controlsPenaltyMPC_Lsd07, controlsPenaltyMPC_dvcc07, controlsPenaltyMPC_yy06, controlsPenaltyMPC_bmy06);
controlsPenaltyMPC_LA_DENSE_BACKWARDSUB_2(controlsPenaltyMPC_Ld06, controlsPenaltyMPC_bmy06, controlsPenaltyMPC_dvcc06);
controlsPenaltyMPC_LA_DENSE_MTVMSUB_2_2(controlsPenaltyMPC_Lsd06, controlsPenaltyMPC_dvcc06, controlsPenaltyMPC_yy05, controlsPenaltyMPC_bmy05);
controlsPenaltyMPC_LA_DENSE_BACKWARDSUB_2(controlsPenaltyMPC_Ld05, controlsPenaltyMPC_bmy05, controlsPenaltyMPC_dvcc05);
controlsPenaltyMPC_LA_DENSE_MTVMSUB_2_2(controlsPenaltyMPC_Lsd05, controlsPenaltyMPC_dvcc05, controlsPenaltyMPC_yy04, controlsPenaltyMPC_bmy04);
controlsPenaltyMPC_LA_DENSE_BACKWARDSUB_2(controlsPenaltyMPC_Ld04, controlsPenaltyMPC_bmy04, controlsPenaltyMPC_dvcc04);
controlsPenaltyMPC_LA_DENSE_MTVMSUB_2_2(controlsPenaltyMPC_Lsd04, controlsPenaltyMPC_dvcc04, controlsPenaltyMPC_yy03, controlsPenaltyMPC_bmy03);
controlsPenaltyMPC_LA_DENSE_BACKWARDSUB_2(controlsPenaltyMPC_Ld03, controlsPenaltyMPC_bmy03, controlsPenaltyMPC_dvcc03);
controlsPenaltyMPC_LA_DENSE_MTVMSUB_2_2(controlsPenaltyMPC_Lsd03, controlsPenaltyMPC_dvcc03, controlsPenaltyMPC_yy02, controlsPenaltyMPC_bmy02);
controlsPenaltyMPC_LA_DENSE_BACKWARDSUB_2(controlsPenaltyMPC_Ld02, controlsPenaltyMPC_bmy02, controlsPenaltyMPC_dvcc02);
controlsPenaltyMPC_LA_DENSE_MTVMSUB_2_2(controlsPenaltyMPC_Lsd02, controlsPenaltyMPC_dvcc02, controlsPenaltyMPC_yy01, controlsPenaltyMPC_bmy01);
controlsPenaltyMPC_LA_DENSE_BACKWARDSUB_2(controlsPenaltyMPC_Ld01, controlsPenaltyMPC_bmy01, controlsPenaltyMPC_dvcc01);
controlsPenaltyMPC_LA_DENSE_MTVMSUB_2_2(controlsPenaltyMPC_Lsd01, controlsPenaltyMPC_dvcc01, controlsPenaltyMPC_yy00, controlsPenaltyMPC_bmy00);
controlsPenaltyMPC_LA_DENSE_BACKWARDSUB_2(controlsPenaltyMPC_Ld00, controlsPenaltyMPC_bmy00, controlsPenaltyMPC_dvcc00);
controlsPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlsPenaltyMPC_C00, controlsPenaltyMPC_dvcc01, controlsPenaltyMPC_D00, controlsPenaltyMPC_dvcc00, controlsPenaltyMPC_grad_eq00);
controlsPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlsPenaltyMPC_C00, controlsPenaltyMPC_dvcc02, controlsPenaltyMPC_D01, controlsPenaltyMPC_dvcc01, controlsPenaltyMPC_grad_eq01);
controlsPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlsPenaltyMPC_C00, controlsPenaltyMPC_dvcc03, controlsPenaltyMPC_D01, controlsPenaltyMPC_dvcc02, controlsPenaltyMPC_grad_eq02);
controlsPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlsPenaltyMPC_C00, controlsPenaltyMPC_dvcc04, controlsPenaltyMPC_D01, controlsPenaltyMPC_dvcc03, controlsPenaltyMPC_grad_eq03);
controlsPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlsPenaltyMPC_C00, controlsPenaltyMPC_dvcc05, controlsPenaltyMPC_D01, controlsPenaltyMPC_dvcc04, controlsPenaltyMPC_grad_eq04);
controlsPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlsPenaltyMPC_C00, controlsPenaltyMPC_dvcc06, controlsPenaltyMPC_D01, controlsPenaltyMPC_dvcc05, controlsPenaltyMPC_grad_eq05);
controlsPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlsPenaltyMPC_C00, controlsPenaltyMPC_dvcc07, controlsPenaltyMPC_D01, controlsPenaltyMPC_dvcc06, controlsPenaltyMPC_grad_eq06);
controlsPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlsPenaltyMPC_C00, controlsPenaltyMPC_dvcc08, controlsPenaltyMPC_D01, controlsPenaltyMPC_dvcc07, controlsPenaltyMPC_grad_eq07);
controlsPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlsPenaltyMPC_C00, controlsPenaltyMPC_dvcc09, controlsPenaltyMPC_D01, controlsPenaltyMPC_dvcc08, controlsPenaltyMPC_grad_eq08);
controlsPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlsPenaltyMPC_C00, controlsPenaltyMPC_dvcc10, controlsPenaltyMPC_D01, controlsPenaltyMPC_dvcc09, controlsPenaltyMPC_grad_eq09);
controlsPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlsPenaltyMPC_C00, controlsPenaltyMPC_dvcc11, controlsPenaltyMPC_D01, controlsPenaltyMPC_dvcc10, controlsPenaltyMPC_grad_eq10);
controlsPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlsPenaltyMPC_C00, controlsPenaltyMPC_dvcc12, controlsPenaltyMPC_D01, controlsPenaltyMPC_dvcc11, controlsPenaltyMPC_grad_eq11);
controlsPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlsPenaltyMPC_C00, controlsPenaltyMPC_dvcc13, controlsPenaltyMPC_D01, controlsPenaltyMPC_dvcc12, controlsPenaltyMPC_grad_eq12);
controlsPenaltyMPC_LA_DIAGZERO_MTVM_2_2(controlsPenaltyMPC_D01, controlsPenaltyMPC_dvcc13, controlsPenaltyMPC_grad_eq13);
controlsPenaltyMPC_LA_VSUB_28(controlsPenaltyMPC_rd, controlsPenaltyMPC_grad_eq, controlsPenaltyMPC_rd);
controlsPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_2(controlsPenaltyMPC_Phi00, controlsPenaltyMPC_rd00, controlsPenaltyMPC_dzcc00);
controlsPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_2(controlsPenaltyMPC_Phi01, controlsPenaltyMPC_rd01, controlsPenaltyMPC_dzcc01);
controlsPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_2(controlsPenaltyMPC_Phi02, controlsPenaltyMPC_rd02, controlsPenaltyMPC_dzcc02);
controlsPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_2(controlsPenaltyMPC_Phi03, controlsPenaltyMPC_rd03, controlsPenaltyMPC_dzcc03);
controlsPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_2(controlsPenaltyMPC_Phi04, controlsPenaltyMPC_rd04, controlsPenaltyMPC_dzcc04);
controlsPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_2(controlsPenaltyMPC_Phi05, controlsPenaltyMPC_rd05, controlsPenaltyMPC_dzcc05);
controlsPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_2(controlsPenaltyMPC_Phi06, controlsPenaltyMPC_rd06, controlsPenaltyMPC_dzcc06);
controlsPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_2(controlsPenaltyMPC_Phi07, controlsPenaltyMPC_rd07, controlsPenaltyMPC_dzcc07);
controlsPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_2(controlsPenaltyMPC_Phi08, controlsPenaltyMPC_rd08, controlsPenaltyMPC_dzcc08);
controlsPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_2(controlsPenaltyMPC_Phi09, controlsPenaltyMPC_rd09, controlsPenaltyMPC_dzcc09);
controlsPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_2(controlsPenaltyMPC_Phi10, controlsPenaltyMPC_rd10, controlsPenaltyMPC_dzcc10);
controlsPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_2(controlsPenaltyMPC_Phi11, controlsPenaltyMPC_rd11, controlsPenaltyMPC_dzcc11);
controlsPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_2(controlsPenaltyMPC_Phi12, controlsPenaltyMPC_rd12, controlsPenaltyMPC_dzcc12);
controlsPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_2(controlsPenaltyMPC_Phi13, controlsPenaltyMPC_rd13, controlsPenaltyMPC_dzcc13);
controlsPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_2(controlsPenaltyMPC_ccrhsl00, controlsPenaltyMPC_slb00, controlsPenaltyMPC_llbbyslb00, controlsPenaltyMPC_dzcc00, controlsPenaltyMPC_lbIdx00, controlsPenaltyMPC_dllbcc00);
controlsPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_2(controlsPenaltyMPC_ccrhsub00, controlsPenaltyMPC_sub00, controlsPenaltyMPC_lubbysub00, controlsPenaltyMPC_dzcc00, controlsPenaltyMPC_ubIdx00, controlsPenaltyMPC_dlubcc00);
controlsPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_2(controlsPenaltyMPC_ccrhsl01, controlsPenaltyMPC_slb01, controlsPenaltyMPC_llbbyslb01, controlsPenaltyMPC_dzcc01, controlsPenaltyMPC_lbIdx01, controlsPenaltyMPC_dllbcc01);
controlsPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_2(controlsPenaltyMPC_ccrhsub01, controlsPenaltyMPC_sub01, controlsPenaltyMPC_lubbysub01, controlsPenaltyMPC_dzcc01, controlsPenaltyMPC_ubIdx01, controlsPenaltyMPC_dlubcc01);
controlsPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_2(controlsPenaltyMPC_ccrhsl02, controlsPenaltyMPC_slb02, controlsPenaltyMPC_llbbyslb02, controlsPenaltyMPC_dzcc02, controlsPenaltyMPC_lbIdx02, controlsPenaltyMPC_dllbcc02);
controlsPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_2(controlsPenaltyMPC_ccrhsub02, controlsPenaltyMPC_sub02, controlsPenaltyMPC_lubbysub02, controlsPenaltyMPC_dzcc02, controlsPenaltyMPC_ubIdx02, controlsPenaltyMPC_dlubcc02);
controlsPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_2(controlsPenaltyMPC_ccrhsl03, controlsPenaltyMPC_slb03, controlsPenaltyMPC_llbbyslb03, controlsPenaltyMPC_dzcc03, controlsPenaltyMPC_lbIdx03, controlsPenaltyMPC_dllbcc03);
controlsPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_2(controlsPenaltyMPC_ccrhsub03, controlsPenaltyMPC_sub03, controlsPenaltyMPC_lubbysub03, controlsPenaltyMPC_dzcc03, controlsPenaltyMPC_ubIdx03, controlsPenaltyMPC_dlubcc03);
controlsPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_2(controlsPenaltyMPC_ccrhsl04, controlsPenaltyMPC_slb04, controlsPenaltyMPC_llbbyslb04, controlsPenaltyMPC_dzcc04, controlsPenaltyMPC_lbIdx04, controlsPenaltyMPC_dllbcc04);
controlsPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_2(controlsPenaltyMPC_ccrhsub04, controlsPenaltyMPC_sub04, controlsPenaltyMPC_lubbysub04, controlsPenaltyMPC_dzcc04, controlsPenaltyMPC_ubIdx04, controlsPenaltyMPC_dlubcc04);
controlsPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_2(controlsPenaltyMPC_ccrhsl05, controlsPenaltyMPC_slb05, controlsPenaltyMPC_llbbyslb05, controlsPenaltyMPC_dzcc05, controlsPenaltyMPC_lbIdx05, controlsPenaltyMPC_dllbcc05);
controlsPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_2(controlsPenaltyMPC_ccrhsub05, controlsPenaltyMPC_sub05, controlsPenaltyMPC_lubbysub05, controlsPenaltyMPC_dzcc05, controlsPenaltyMPC_ubIdx05, controlsPenaltyMPC_dlubcc05);
controlsPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_2(controlsPenaltyMPC_ccrhsl06, controlsPenaltyMPC_slb06, controlsPenaltyMPC_llbbyslb06, controlsPenaltyMPC_dzcc06, controlsPenaltyMPC_lbIdx06, controlsPenaltyMPC_dllbcc06);
controlsPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_2(controlsPenaltyMPC_ccrhsub06, controlsPenaltyMPC_sub06, controlsPenaltyMPC_lubbysub06, controlsPenaltyMPC_dzcc06, controlsPenaltyMPC_ubIdx06, controlsPenaltyMPC_dlubcc06);
controlsPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_2(controlsPenaltyMPC_ccrhsl07, controlsPenaltyMPC_slb07, controlsPenaltyMPC_llbbyslb07, controlsPenaltyMPC_dzcc07, controlsPenaltyMPC_lbIdx07, controlsPenaltyMPC_dllbcc07);
controlsPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_2(controlsPenaltyMPC_ccrhsub07, controlsPenaltyMPC_sub07, controlsPenaltyMPC_lubbysub07, controlsPenaltyMPC_dzcc07, controlsPenaltyMPC_ubIdx07, controlsPenaltyMPC_dlubcc07);
controlsPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_2(controlsPenaltyMPC_ccrhsl08, controlsPenaltyMPC_slb08, controlsPenaltyMPC_llbbyslb08, controlsPenaltyMPC_dzcc08, controlsPenaltyMPC_lbIdx08, controlsPenaltyMPC_dllbcc08);
controlsPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_2(controlsPenaltyMPC_ccrhsub08, controlsPenaltyMPC_sub08, controlsPenaltyMPC_lubbysub08, controlsPenaltyMPC_dzcc08, controlsPenaltyMPC_ubIdx08, controlsPenaltyMPC_dlubcc08);
controlsPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_2(controlsPenaltyMPC_ccrhsl09, controlsPenaltyMPC_slb09, controlsPenaltyMPC_llbbyslb09, controlsPenaltyMPC_dzcc09, controlsPenaltyMPC_lbIdx09, controlsPenaltyMPC_dllbcc09);
controlsPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_2(controlsPenaltyMPC_ccrhsub09, controlsPenaltyMPC_sub09, controlsPenaltyMPC_lubbysub09, controlsPenaltyMPC_dzcc09, controlsPenaltyMPC_ubIdx09, controlsPenaltyMPC_dlubcc09);
controlsPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_2(controlsPenaltyMPC_ccrhsl10, controlsPenaltyMPC_slb10, controlsPenaltyMPC_llbbyslb10, controlsPenaltyMPC_dzcc10, controlsPenaltyMPC_lbIdx10, controlsPenaltyMPC_dllbcc10);
controlsPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_2(controlsPenaltyMPC_ccrhsub10, controlsPenaltyMPC_sub10, controlsPenaltyMPC_lubbysub10, controlsPenaltyMPC_dzcc10, controlsPenaltyMPC_ubIdx10, controlsPenaltyMPC_dlubcc10);
controlsPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_2(controlsPenaltyMPC_ccrhsl11, controlsPenaltyMPC_slb11, controlsPenaltyMPC_llbbyslb11, controlsPenaltyMPC_dzcc11, controlsPenaltyMPC_lbIdx11, controlsPenaltyMPC_dllbcc11);
controlsPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_2(controlsPenaltyMPC_ccrhsub11, controlsPenaltyMPC_sub11, controlsPenaltyMPC_lubbysub11, controlsPenaltyMPC_dzcc11, controlsPenaltyMPC_ubIdx11, controlsPenaltyMPC_dlubcc11);
controlsPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_2(controlsPenaltyMPC_ccrhsl12, controlsPenaltyMPC_slb12, controlsPenaltyMPC_llbbyslb12, controlsPenaltyMPC_dzcc12, controlsPenaltyMPC_lbIdx12, controlsPenaltyMPC_dllbcc12);
controlsPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_2(controlsPenaltyMPC_ccrhsub12, controlsPenaltyMPC_sub12, controlsPenaltyMPC_lubbysub12, controlsPenaltyMPC_dzcc12, controlsPenaltyMPC_ubIdx12, controlsPenaltyMPC_dlubcc12);
controlsPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_2(controlsPenaltyMPC_ccrhsl13, controlsPenaltyMPC_slb13, controlsPenaltyMPC_llbbyslb13, controlsPenaltyMPC_dzcc13, controlsPenaltyMPC_lbIdx13, controlsPenaltyMPC_dllbcc13);
controlsPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_2(controlsPenaltyMPC_ccrhsub13, controlsPenaltyMPC_sub13, controlsPenaltyMPC_lubbysub13, controlsPenaltyMPC_dzcc13, controlsPenaltyMPC_ubIdx13, controlsPenaltyMPC_dlubcc13);
controlsPenaltyMPC_LA_VSUB7_56(controlsPenaltyMPC_l, controlsPenaltyMPC_ccrhs, controlsPenaltyMPC_s, controlsPenaltyMPC_dl_cc, controlsPenaltyMPC_ds_cc);
controlsPenaltyMPC_LA_VADD_28(controlsPenaltyMPC_dz_cc, controlsPenaltyMPC_dz_aff);
controlsPenaltyMPC_LA_VADD_28(controlsPenaltyMPC_dv_cc, controlsPenaltyMPC_dv_aff);
controlsPenaltyMPC_LA_VADD_56(controlsPenaltyMPC_dl_cc, controlsPenaltyMPC_dl_aff);
controlsPenaltyMPC_LA_VADD_56(controlsPenaltyMPC_ds_cc, controlsPenaltyMPC_ds_aff);
info->lsit_cc = controlsPenaltyMPC_LINESEARCH_BACKTRACKING_COMBINED(controlsPenaltyMPC_z, controlsPenaltyMPC_v, controlsPenaltyMPC_l, controlsPenaltyMPC_s, controlsPenaltyMPC_dz_cc, controlsPenaltyMPC_dv_cc, controlsPenaltyMPC_dl_cc, controlsPenaltyMPC_ds_cc, &info->step_cc, &info->mu);
if( info->lsit_cc == controlsPenaltyMPC_NOPROGRESS ){
exitcode = controlsPenaltyMPC_NOPROGRESS; break;
}
info->it++;
}
output->z1[0] = controlsPenaltyMPC_z00[0];
output->z1[1] = controlsPenaltyMPC_z00[1];
output->z2[0] = controlsPenaltyMPC_z01[0];
output->z2[1] = controlsPenaltyMPC_z01[1];
output->z3[0] = controlsPenaltyMPC_z02[0];
output->z3[1] = controlsPenaltyMPC_z02[1];
output->z4[0] = controlsPenaltyMPC_z03[0];
output->z4[1] = controlsPenaltyMPC_z03[1];
output->z5[0] = controlsPenaltyMPC_z04[0];
output->z5[1] = controlsPenaltyMPC_z04[1];
output->z6[0] = controlsPenaltyMPC_z05[0];
output->z6[1] = controlsPenaltyMPC_z05[1];
output->z7[0] = controlsPenaltyMPC_z06[0];
output->z7[1] = controlsPenaltyMPC_z06[1];
output->z8[0] = controlsPenaltyMPC_z07[0];
output->z8[1] = controlsPenaltyMPC_z07[1];
output->z9[0] = controlsPenaltyMPC_z08[0];
output->z9[1] = controlsPenaltyMPC_z08[1];
output->z10[0] = controlsPenaltyMPC_z09[0];
output->z10[1] = controlsPenaltyMPC_z09[1];
output->z11[0] = controlsPenaltyMPC_z10[0];
output->z11[1] = controlsPenaltyMPC_z10[1];
output->z12[0] = controlsPenaltyMPC_z11[0];
output->z12[1] = controlsPenaltyMPC_z11[1];
output->z13[0] = controlsPenaltyMPC_z12[0];
output->z13[1] = controlsPenaltyMPC_z12[1];
output->z14[0] = controlsPenaltyMPC_z13[0];
output->z14[1] = controlsPenaltyMPC_z13[1];

#if controlsPenaltyMPC_SET_TIMING == 1
info->solvetime = controlsPenaltyMPC_toc(&solvertimer);
#if controlsPenaltyMPC_SET_PRINTLEVEL > 0 && controlsPenaltyMPC_SET_TIMING == 1
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
