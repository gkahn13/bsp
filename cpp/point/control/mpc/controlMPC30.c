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

#include "../include/controlMPC.h"

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
 * Initializes a vector of length 58 with a value.
 */
void controlMPC_LA_INITIALIZEVECTOR_58(controlMPC_FLOAT* vec, controlMPC_FLOAT value)
{
	int i;
	for( i=0; i<58; i++ )
	{
		vec[i] = value;
	}
}


/*
 * Initializes a vector of length 56 with a value.
 */
void controlMPC_LA_INITIALIZEVECTOR_56(controlMPC_FLOAT* vec, controlMPC_FLOAT value)
{
	int i;
	for( i=0; i<56; i++ )
	{
		vec[i] = value;
	}
}


/*
 * Initializes a vector of length 116 with a value.
 */
void controlMPC_LA_INITIALIZEVECTOR_116(controlMPC_FLOAT* vec, controlMPC_FLOAT value)
{
	int i;
	for( i=0; i<116; i++ )
	{
		vec[i] = value;
	}
}


/* 
 * Calculates a dot product and adds it to a variable: z += x'*y; 
 * This function is for vectors of length 116.
 */
void controlMPC_LA_DOTACC_116(controlMPC_FLOAT *x, controlMPC_FLOAT *y, controlMPC_FLOAT *z)
{
	int i;
	for( i=0; i<116; i++ ){
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
void controlMPC_LA_DIAG_QUADFCN_2(controlMPC_FLOAT* H, controlMPC_FLOAT* f, controlMPC_FLOAT* z, controlMPC_FLOAT* grad, controlMPC_FLOAT* value)
{
	int i;
	controlMPC_FLOAT hz;	
	for( i=0; i<2; i++){
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
void controlMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_2_2(controlMPC_FLOAT *A, controlMPC_FLOAT *x, controlMPC_FLOAT *B, controlMPC_FLOAT *u, controlMPC_FLOAT *b, controlMPC_FLOAT *l, controlMPC_FLOAT *r, controlMPC_FLOAT *z, controlMPC_FLOAT *y)
{
	int i;
	int j;
	int k = 0;
	controlMPC_FLOAT AxBu[2];
	controlMPC_FLOAT norm = *y;
	controlMPC_FLOAT lr = 0;

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
 * Matrix vector multiplication y = M'*x where M is of size [2 x 2]
 * and stored in column major format. Note the transpose of M!
 */
void controlMPC_LA_DENSE_MTVM_2_2(controlMPC_FLOAT *M, controlMPC_FLOAT *x, controlMPC_FLOAT *y)
{
	int i;
	int j;
	int k = 0; 
	for( i=0; i<2; i++ ){
		y[i] = 0;
		for( j=0; j<2; j++ ){
			y[i] += M[k++]*x[j];
		}
	}
}


/*
 * Matrix vector multiplication z = A'*x + B'*y 
 * where A is of size [2 x 2] and stored in column major format.
 * and B is of size [2 x 2] and stored in diagzero format
 * Note the transposes of A and B!
 */
void controlMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlMPC_FLOAT *A, controlMPC_FLOAT *x, controlMPC_FLOAT *B, controlMPC_FLOAT *y, controlMPC_FLOAT *z)
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
void controlMPC_LA_DIAGZERO_MTVM_2_2(controlMPC_FLOAT *M, controlMPC_FLOAT *x, controlMPC_FLOAT *y)
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
void controlMPC_LA_VSUBADD3_2(controlMPC_FLOAT* t, controlMPC_FLOAT* u, int* uidx, controlMPC_FLOAT* v, controlMPC_FLOAT* w, controlMPC_FLOAT* y, controlMPC_FLOAT* z, controlMPC_FLOAT* r)
{
	int i;
	controlMPC_FLOAT norm = *r;
	controlMPC_FLOAT vx = 0;
	controlMPC_FLOAT x;
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
void controlMPC_LA_VSUBADD2_2(controlMPC_FLOAT* t, int* tidx, controlMPC_FLOAT* u, controlMPC_FLOAT* v, controlMPC_FLOAT* w, controlMPC_FLOAT* y, controlMPC_FLOAT* z, controlMPC_FLOAT* r)
{
	int i;
	controlMPC_FLOAT norm = *r;
	controlMPC_FLOAT vx = 0;
	controlMPC_FLOAT x;
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
void controlMPC_LA_INEQ_B_GRAD_2_2_2(controlMPC_FLOAT *lu, controlMPC_FLOAT *su, controlMPC_FLOAT *ru, controlMPC_FLOAT *ll, controlMPC_FLOAT *sl, controlMPC_FLOAT *rl, int* lbIdx, int* ubIdx, controlMPC_FLOAT *grad, controlMPC_FLOAT *lubysu, controlMPC_FLOAT *llbysl)
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
 * of length 58.
 */
void controlMPC_LA_VVADD3_58(controlMPC_FLOAT *u, controlMPC_FLOAT *v, controlMPC_FLOAT *w, controlMPC_FLOAT *z)
{
	int i;
	for( i=0; i<58; i++ ){
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
void controlMPC_LA_DIAG_CHOL_ONELOOP_LBUB_2_2_2(controlMPC_FLOAT *H, controlMPC_FLOAT *llbysl, int* lbIdx, controlMPC_FLOAT *lubysu, int* ubIdx, controlMPC_FLOAT *Phi)


{
	int i;
	
	/* compute cholesky */
	for( i=0; i<2; i++ ){
		Phi[i] = H[i] + llbysl[i] + lubysu[i];

#if controlMPC_SET_PRINTLEVEL > 0 && defined PRINTNUMERICALWARNINGS
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
void controlMPC_LA_DIAG_MATRIXFORWARDSUB_2_2(controlMPC_FLOAT *L, controlMPC_FLOAT *B, controlMPC_FLOAT *A)
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
 * Forward substitution to solve L*y = b where L is a
 * diagonal matrix in vector storage format.
 * 
 * The dimensions involved are 2.
 */
void controlMPC_LA_DIAG_FORWARDSUB_2(controlMPC_FLOAT *L, controlMPC_FLOAT *b, controlMPC_FLOAT *y)
{
    int i;

    for( i=0; i<2; i++ ){
		y[i] = b[i]/L[i];
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
void controlMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_2(controlMPC_FLOAT *L, controlMPC_FLOAT *B, controlMPC_FLOAT *A)
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
void controlMPC_LA_DENSE_DIAGZERO_MMTM_2_2_2(controlMPC_FLOAT *A, controlMPC_FLOAT *B, controlMPC_FLOAT *C)
{
    int i, j;
	
	for( i=0; i<2; i++ ){
		for( j=0; j<2; j++){
			C[j*2+i] = B[i*2+j]*A[i];
		}
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
void controlMPC_LA_DENSE_DIAGZERO_MMT2_2_2_2(controlMPC_FLOAT *A, controlMPC_FLOAT *B, controlMPC_FLOAT *L)
{
    int i, j, k, ii, di;
    controlMPC_FLOAT ltemp;
    
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
void controlMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_2_2(controlMPC_FLOAT *A, controlMPC_FLOAT *x, controlMPC_FLOAT *B, controlMPC_FLOAT *u, controlMPC_FLOAT *b, controlMPC_FLOAT *r)
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
void controlMPC_LA_DENSE_CHOL_2(controlMPC_FLOAT *A, controlMPC_FLOAT *L)
{
    int i, j, k, di, dj;
	 int ii, jj;

    controlMPC_FLOAT l;
    controlMPC_FLOAT Mii;

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
        
#if controlMPC_SET_PRINTLEVEL > 0 && defined PRINTNUMERICALWARNINGS
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
void controlMPC_LA_DENSE_FORWARDSUB_2(controlMPC_FLOAT *L, controlMPC_FLOAT *b, controlMPC_FLOAT *y)
{
    int i,j,ii,di;
    controlMPC_FLOAT yel;
            
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
void controlMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(controlMPC_FLOAT *L, controlMPC_FLOAT *B, controlMPC_FLOAT *A)
{
    int i,j,k,ii,di;
    controlMPC_FLOAT a;
    
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
void controlMPC_LA_DENSE_MMTSUB_2_2(controlMPC_FLOAT *A, controlMPC_FLOAT *L)
{
    int i, j, k, ii, di;
    controlMPC_FLOAT ltemp;
    
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
void controlMPC_LA_DENSE_MVMSUB1_2_2(controlMPC_FLOAT *A, controlMPC_FLOAT *x, controlMPC_FLOAT *b, controlMPC_FLOAT *r)
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
void controlMPC_LA_DENSE_BACKWARDSUB_2(controlMPC_FLOAT *L, controlMPC_FLOAT *y, controlMPC_FLOAT *x)
{
    int i, ii, di, j, jj, dj;
    controlMPC_FLOAT xel;    
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
void controlMPC_LA_DENSE_MTVMSUB_2_2(controlMPC_FLOAT *A, controlMPC_FLOAT *x, controlMPC_FLOAT *b, controlMPC_FLOAT *r)
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
 * Vector subtraction z = -x - y for vectors of length 58.
 */
void controlMPC_LA_VSUB2_58(controlMPC_FLOAT *x, controlMPC_FLOAT *y, controlMPC_FLOAT *z)
{
	int i;
	for( i=0; i<58; i++){
		z[i] = -x[i] - y[i];
	}
}


/**
 * Forward-Backward-Substitution to solve L*L^T*x = b where L is a
 * diagonal matrix of size 2 in vector
 * storage format.
 */
void controlMPC_LA_DIAG_FORWARDBACKWARDSUB_2(controlMPC_FLOAT *L, controlMPC_FLOAT *b, controlMPC_FLOAT *x)
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
void controlMPC_LA_VSUB_INDEXED_2(controlMPC_FLOAT *x, int* xidx, controlMPC_FLOAT *y, controlMPC_FLOAT *z)
{
	int i;
	for( i=0; i<2; i++){
		z[i] = x[xidx[i]] - y[i];
	}
}


/*
 * Vector subtraction x = -u.*v - w for vectors of length 2.
 */
void controlMPC_LA_VSUB3_2(controlMPC_FLOAT *u, controlMPC_FLOAT *v, controlMPC_FLOAT *w, controlMPC_FLOAT *x)
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
void controlMPC_LA_VSUB2_INDEXED_2(controlMPC_FLOAT *x, controlMPC_FLOAT *y, int* yidx, controlMPC_FLOAT *z)
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
 * controlMPC_NOPROGRESS (should be negative).
 */
int controlMPC_LINESEARCH_BACKTRACKING_AFFINE(controlMPC_FLOAT *l, controlMPC_FLOAT *s, controlMPC_FLOAT *dl, controlMPC_FLOAT *ds, controlMPC_FLOAT *a, controlMPC_FLOAT *mu_aff)
{
    int i;
	int lsIt=1;    
    controlMPC_FLOAT dltemp;
    controlMPC_FLOAT dstemp;
    controlMPC_FLOAT mya = 1.0;
    controlMPC_FLOAT mymu;
        
    while( 1 ){                        

        /* 
         * Compute both snew and wnew together.
         * We compute also mu_affine along the way here, as the
         * values might be in registers, so it should be cheaper.
         */
        mymu = 0;
        for( i=0; i<116; i++ ){
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
        if( i == 116 ){
            break;
        } else {
            mya *= controlMPC_SET_LS_SCALE_AFF;
            if( mya < controlMPC_SET_LS_MINSTEP ){
                return controlMPC_NOPROGRESS;
            }
        }
    }
    
    /* return new values and iteration counter */
    *a = mya;
    *mu_aff = mymu / (controlMPC_FLOAT)116;
    return lsIt;
}


/*
 * Vector subtraction x = u.*v - a where a is a scalar
*  and x,u,v are vectors of length 116.
 */
void controlMPC_LA_VSUB5_116(controlMPC_FLOAT *u, controlMPC_FLOAT *v, controlMPC_FLOAT a, controlMPC_FLOAT *x)
{
	int i;
	for( i=0; i<116; i++){
		x[i] = u[i]*v[i] - a;
	}
}


/*
 * Computes x=0; x(uidx) += u/su; x(vidx) -= v/sv where x is of length 2,
 * u, su, uidx are of length 2 and v, sv, vidx are of length 2.
 */
void controlMPC_LA_VSUB6_INDEXED_2_2_2(controlMPC_FLOAT *u, controlMPC_FLOAT *su, int* uidx, controlMPC_FLOAT *v, controlMPC_FLOAT *sv, int* vidx, controlMPC_FLOAT *x)
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
 * Computes r = A*x + B*u
 * where A is stored in column major format
 * and B is stored in diagzero format
 */
void controlMPC_LA_DENSE_DIAGZERO_2MVMADD_2_2_2(controlMPC_FLOAT *A, controlMPC_FLOAT *x, controlMPC_FLOAT *B, controlMPC_FLOAT *u, controlMPC_FLOAT *r)
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
 * Vector subtraction z = x - y for vectors of length 58.
 */
void controlMPC_LA_VSUB_58(controlMPC_FLOAT *x, controlMPC_FLOAT *y, controlMPC_FLOAT *z)
{
	int i;
	for( i=0; i<58; i++){
		z[i] = x[i] - y[i];
	}
}


/** 
 * Computes z = -r./s - u.*y(y)
 * where all vectors except of y are of length 2 (length of y >= 2).
 */
void controlMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_2(controlMPC_FLOAT *r, controlMPC_FLOAT *s, controlMPC_FLOAT *u, controlMPC_FLOAT *y, int* yidx, controlMPC_FLOAT *z)
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
void controlMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_2(controlMPC_FLOAT *r, controlMPC_FLOAT *s, controlMPC_FLOAT *u, controlMPC_FLOAT *y, int* yidx, controlMPC_FLOAT *z)
{
	int i;
	for( i=0; i<2; i++ ){
		z[i] = -r[i]/s[i] + u[i]*y[yidx[i]];
	}
}


/*
 * Computes ds = -l.\(r + s.*dl) for vectors of length 116.
 */
void controlMPC_LA_VSUB7_116(controlMPC_FLOAT *l, controlMPC_FLOAT *r, controlMPC_FLOAT *s, controlMPC_FLOAT *dl, controlMPC_FLOAT *ds)
{
	int i;
	for( i=0; i<116; i++){
		ds[i] = -(r[i] + s[i]*dl[i])/l[i];
	}
}


/*
 * Vector addition x = x + y for vectors of length 58.
 */
void controlMPC_LA_VADD_58(controlMPC_FLOAT *x, controlMPC_FLOAT *y)
{
	int i;
	for( i=0; i<58; i++){
		x[i] += y[i];
	}
}


/*
 * Vector addition x = x + y for vectors of length 56.
 */
void controlMPC_LA_VADD_56(controlMPC_FLOAT *x, controlMPC_FLOAT *y)
{
	int i;
	for( i=0; i<56; i++){
		x[i] += y[i];
	}
}


/*
 * Vector addition x = x + y for vectors of length 116.
 */
void controlMPC_LA_VADD_116(controlMPC_FLOAT *x, controlMPC_FLOAT *y)
{
	int i;
	for( i=0; i<116; i++){
		x[i] += y[i];
	}
}


/**
 * Backtracking line search for combined predictor/corrector step.
 * Update on variables with safety factor gamma (to keep us away from
 * boundary).
 */
int controlMPC_LINESEARCH_BACKTRACKING_COMBINED(controlMPC_FLOAT *z, controlMPC_FLOAT *v, controlMPC_FLOAT *l, controlMPC_FLOAT *s, controlMPC_FLOAT *dz, controlMPC_FLOAT *dv, controlMPC_FLOAT *dl, controlMPC_FLOAT *ds, controlMPC_FLOAT *a, controlMPC_FLOAT *mu)
{
    int i, lsIt=1;       
    controlMPC_FLOAT dltemp;
    controlMPC_FLOAT dstemp;    
    controlMPC_FLOAT a_gamma;
            
    *a = 1.0;
    while( 1 ){                        

        /* check whether search criterion is fulfilled */
        for( i=0; i<116; i++ ){
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
        if( i == 116 ){
            break;
        } else {
            *a *= controlMPC_SET_LS_SCALE;
            if( *a < controlMPC_SET_LS_MINSTEP ){
                return controlMPC_NOPROGRESS;
            }
        }
    }
    
    /* update variables with safety margin */
    a_gamma = (*a)*controlMPC_SET_LS_MAXSTEP;
    
    /* primal variables */
    for( i=0; i<58; i++ ){
        z[i] += a_gamma*dz[i];
    }
    
    /* equality constraint multipliers */
    for( i=0; i<56; i++ ){
        v[i] += a_gamma*dv[i];
    }
    
    /* inequality constraint multipliers & slacks, also update mu */
    *mu = 0;
    for( i=0; i<116; i++ ){
        dltemp = l[i] + a_gamma*dl[i]; l[i] = dltemp;
        dstemp = s[i] + a_gamma*ds[i]; s[i] = dstemp;
        *mu += dltemp*dstemp;
    }
    
    *a = a_gamma;
    *mu /= (controlMPC_FLOAT)116;
    return lsIt;
}




/* VARIABLE DEFINITIONS ------------------------------------------------ */
controlMPC_FLOAT controlMPC_z[58];
controlMPC_FLOAT controlMPC_v[56];
controlMPC_FLOAT controlMPC_dz_aff[58];
controlMPC_FLOAT controlMPC_dv_aff[56];
controlMPC_FLOAT controlMPC_grad_cost[58];
controlMPC_FLOAT controlMPC_grad_eq[58];
controlMPC_FLOAT controlMPC_rd[58];
controlMPC_FLOAT controlMPC_l[116];
controlMPC_FLOAT controlMPC_s[116];
controlMPC_FLOAT controlMPC_lbys[116];
controlMPC_FLOAT controlMPC_dl_aff[116];
controlMPC_FLOAT controlMPC_ds_aff[116];
controlMPC_FLOAT controlMPC_dz_cc[58];
controlMPC_FLOAT controlMPC_dv_cc[56];
controlMPC_FLOAT controlMPC_dl_cc[116];
controlMPC_FLOAT controlMPC_ds_cc[116];
controlMPC_FLOAT controlMPC_ccrhs[116];
controlMPC_FLOAT controlMPC_grad_ineq[58];
controlMPC_FLOAT* controlMPC_z00 = controlMPC_z + 0;
controlMPC_FLOAT* controlMPC_dzaff00 = controlMPC_dz_aff + 0;
controlMPC_FLOAT* controlMPC_dzcc00 = controlMPC_dz_cc + 0;
controlMPC_FLOAT* controlMPC_rd00 = controlMPC_rd + 0;
controlMPC_FLOAT controlMPC_Lbyrd00[2];
controlMPC_FLOAT* controlMPC_grad_cost00 = controlMPC_grad_cost + 0;
controlMPC_FLOAT* controlMPC_grad_eq00 = controlMPC_grad_eq + 0;
controlMPC_FLOAT* controlMPC_grad_ineq00 = controlMPC_grad_ineq + 0;
controlMPC_FLOAT controlMPC_ctv00[2];
controlMPC_FLOAT controlMPC_C00[4] = {0.0000000000000000E+000, 0.0000000000000000E+000, 
0.0000000000000000E+000, 0.0000000000000000E+000};
controlMPC_FLOAT* controlMPC_v00 = controlMPC_v + 0;
controlMPC_FLOAT controlMPC_re00[2];
controlMPC_FLOAT controlMPC_beta00[2];
controlMPC_FLOAT controlMPC_betacc00[2];
controlMPC_FLOAT* controlMPC_dvaff00 = controlMPC_dv_aff + 0;
controlMPC_FLOAT* controlMPC_dvcc00 = controlMPC_dv_cc + 0;
controlMPC_FLOAT controlMPC_V00[4];
controlMPC_FLOAT controlMPC_Yd00[3];
controlMPC_FLOAT controlMPC_Ld00[3];
controlMPC_FLOAT controlMPC_yy00[2];
controlMPC_FLOAT controlMPC_bmy00[2];
controlMPC_FLOAT controlMPC_c00[2] = {0.0000000000000000E+000, 0.0000000000000000E+000};
int controlMPC_lbIdx00[2] = {0, 1};
controlMPC_FLOAT* controlMPC_llb00 = controlMPC_l + 0;
controlMPC_FLOAT* controlMPC_slb00 = controlMPC_s + 0;
controlMPC_FLOAT* controlMPC_llbbyslb00 = controlMPC_lbys + 0;
controlMPC_FLOAT controlMPC_rilb00[2];
controlMPC_FLOAT* controlMPC_dllbaff00 = controlMPC_dl_aff + 0;
controlMPC_FLOAT* controlMPC_dslbaff00 = controlMPC_ds_aff + 0;
controlMPC_FLOAT* controlMPC_dllbcc00 = controlMPC_dl_cc + 0;
controlMPC_FLOAT* controlMPC_dslbcc00 = controlMPC_ds_cc + 0;
controlMPC_FLOAT* controlMPC_ccrhsl00 = controlMPC_ccrhs + 0;
int controlMPC_ubIdx00[2] = {0, 1};
controlMPC_FLOAT* controlMPC_lub00 = controlMPC_l + 2;
controlMPC_FLOAT* controlMPC_sub00 = controlMPC_s + 2;
controlMPC_FLOAT* controlMPC_lubbysub00 = controlMPC_lbys + 2;
controlMPC_FLOAT controlMPC_riub00[2];
controlMPC_FLOAT* controlMPC_dlubaff00 = controlMPC_dl_aff + 2;
controlMPC_FLOAT* controlMPC_dsubaff00 = controlMPC_ds_aff + 2;
controlMPC_FLOAT* controlMPC_dlubcc00 = controlMPC_dl_cc + 2;
controlMPC_FLOAT* controlMPC_dsubcc00 = controlMPC_ds_cc + 2;
controlMPC_FLOAT* controlMPC_ccrhsub00 = controlMPC_ccrhs + 2;
controlMPC_FLOAT controlMPC_Phi00[2];
controlMPC_FLOAT* controlMPC_z01 = controlMPC_z + 2;
controlMPC_FLOAT* controlMPC_dzaff01 = controlMPC_dz_aff + 2;
controlMPC_FLOAT* controlMPC_dzcc01 = controlMPC_dz_cc + 2;
controlMPC_FLOAT* controlMPC_rd01 = controlMPC_rd + 2;
controlMPC_FLOAT controlMPC_Lbyrd01[2];
controlMPC_FLOAT* controlMPC_grad_cost01 = controlMPC_grad_cost + 2;
controlMPC_FLOAT* controlMPC_grad_eq01 = controlMPC_grad_eq + 2;
controlMPC_FLOAT* controlMPC_grad_ineq01 = controlMPC_grad_ineq + 2;
controlMPC_FLOAT controlMPC_ctv01[2];
controlMPC_FLOAT* controlMPC_v01 = controlMPC_v + 2;
controlMPC_FLOAT controlMPC_re01[2];
controlMPC_FLOAT controlMPC_beta01[2];
controlMPC_FLOAT controlMPC_betacc01[2];
controlMPC_FLOAT* controlMPC_dvaff01 = controlMPC_dv_aff + 2;
controlMPC_FLOAT* controlMPC_dvcc01 = controlMPC_dv_cc + 2;
controlMPC_FLOAT controlMPC_V01[4];
controlMPC_FLOAT controlMPC_Yd01[3];
controlMPC_FLOAT controlMPC_Ld01[3];
controlMPC_FLOAT controlMPC_yy01[2];
controlMPC_FLOAT controlMPC_bmy01[2];
controlMPC_FLOAT controlMPC_c01[2] = {0.0000000000000000E+000, 0.0000000000000000E+000};
int controlMPC_lbIdx01[2] = {0, 1};
controlMPC_FLOAT* controlMPC_llb01 = controlMPC_l + 4;
controlMPC_FLOAT* controlMPC_slb01 = controlMPC_s + 4;
controlMPC_FLOAT* controlMPC_llbbyslb01 = controlMPC_lbys + 4;
controlMPC_FLOAT controlMPC_rilb01[2];
controlMPC_FLOAT* controlMPC_dllbaff01 = controlMPC_dl_aff + 4;
controlMPC_FLOAT* controlMPC_dslbaff01 = controlMPC_ds_aff + 4;
controlMPC_FLOAT* controlMPC_dllbcc01 = controlMPC_dl_cc + 4;
controlMPC_FLOAT* controlMPC_dslbcc01 = controlMPC_ds_cc + 4;
controlMPC_FLOAT* controlMPC_ccrhsl01 = controlMPC_ccrhs + 4;
int controlMPC_ubIdx01[2] = {0, 1};
controlMPC_FLOAT* controlMPC_lub01 = controlMPC_l + 6;
controlMPC_FLOAT* controlMPC_sub01 = controlMPC_s + 6;
controlMPC_FLOAT* controlMPC_lubbysub01 = controlMPC_lbys + 6;
controlMPC_FLOAT controlMPC_riub01[2];
controlMPC_FLOAT* controlMPC_dlubaff01 = controlMPC_dl_aff + 6;
controlMPC_FLOAT* controlMPC_dsubaff01 = controlMPC_ds_aff + 6;
controlMPC_FLOAT* controlMPC_dlubcc01 = controlMPC_dl_cc + 6;
controlMPC_FLOAT* controlMPC_dsubcc01 = controlMPC_ds_cc + 6;
controlMPC_FLOAT* controlMPC_ccrhsub01 = controlMPC_ccrhs + 6;
controlMPC_FLOAT controlMPC_Phi01[2];
controlMPC_FLOAT controlMPC_D01[2] = {0.0000000000000000E+000, 
0.0000000000000000E+000};
controlMPC_FLOAT controlMPC_W01[2];
controlMPC_FLOAT controlMPC_Ysd01[4];
controlMPC_FLOAT controlMPC_Lsd01[4];
controlMPC_FLOAT* controlMPC_z02 = controlMPC_z + 4;
controlMPC_FLOAT* controlMPC_dzaff02 = controlMPC_dz_aff + 4;
controlMPC_FLOAT* controlMPC_dzcc02 = controlMPC_dz_cc + 4;
controlMPC_FLOAT* controlMPC_rd02 = controlMPC_rd + 4;
controlMPC_FLOAT controlMPC_Lbyrd02[2];
controlMPC_FLOAT* controlMPC_grad_cost02 = controlMPC_grad_cost + 4;
controlMPC_FLOAT* controlMPC_grad_eq02 = controlMPC_grad_eq + 4;
controlMPC_FLOAT* controlMPC_grad_ineq02 = controlMPC_grad_ineq + 4;
controlMPC_FLOAT controlMPC_ctv02[2];
controlMPC_FLOAT* controlMPC_v02 = controlMPC_v + 4;
controlMPC_FLOAT controlMPC_re02[2];
controlMPC_FLOAT controlMPC_beta02[2];
controlMPC_FLOAT controlMPC_betacc02[2];
controlMPC_FLOAT* controlMPC_dvaff02 = controlMPC_dv_aff + 4;
controlMPC_FLOAT* controlMPC_dvcc02 = controlMPC_dv_cc + 4;
controlMPC_FLOAT controlMPC_V02[4];
controlMPC_FLOAT controlMPC_Yd02[3];
controlMPC_FLOAT controlMPC_Ld02[3];
controlMPC_FLOAT controlMPC_yy02[2];
controlMPC_FLOAT controlMPC_bmy02[2];
controlMPC_FLOAT controlMPC_c02[2] = {0.0000000000000000E+000, 0.0000000000000000E+000};
int controlMPC_lbIdx02[2] = {0, 1};
controlMPC_FLOAT* controlMPC_llb02 = controlMPC_l + 8;
controlMPC_FLOAT* controlMPC_slb02 = controlMPC_s + 8;
controlMPC_FLOAT* controlMPC_llbbyslb02 = controlMPC_lbys + 8;
controlMPC_FLOAT controlMPC_rilb02[2];
controlMPC_FLOAT* controlMPC_dllbaff02 = controlMPC_dl_aff + 8;
controlMPC_FLOAT* controlMPC_dslbaff02 = controlMPC_ds_aff + 8;
controlMPC_FLOAT* controlMPC_dllbcc02 = controlMPC_dl_cc + 8;
controlMPC_FLOAT* controlMPC_dslbcc02 = controlMPC_ds_cc + 8;
controlMPC_FLOAT* controlMPC_ccrhsl02 = controlMPC_ccrhs + 8;
int controlMPC_ubIdx02[2] = {0, 1};
controlMPC_FLOAT* controlMPC_lub02 = controlMPC_l + 10;
controlMPC_FLOAT* controlMPC_sub02 = controlMPC_s + 10;
controlMPC_FLOAT* controlMPC_lubbysub02 = controlMPC_lbys + 10;
controlMPC_FLOAT controlMPC_riub02[2];
controlMPC_FLOAT* controlMPC_dlubaff02 = controlMPC_dl_aff + 10;
controlMPC_FLOAT* controlMPC_dsubaff02 = controlMPC_ds_aff + 10;
controlMPC_FLOAT* controlMPC_dlubcc02 = controlMPC_dl_cc + 10;
controlMPC_FLOAT* controlMPC_dsubcc02 = controlMPC_ds_cc + 10;
controlMPC_FLOAT* controlMPC_ccrhsub02 = controlMPC_ccrhs + 10;
controlMPC_FLOAT controlMPC_Phi02[2];
controlMPC_FLOAT controlMPC_W02[2];
controlMPC_FLOAT controlMPC_Ysd02[4];
controlMPC_FLOAT controlMPC_Lsd02[4];
controlMPC_FLOAT* controlMPC_z03 = controlMPC_z + 6;
controlMPC_FLOAT* controlMPC_dzaff03 = controlMPC_dz_aff + 6;
controlMPC_FLOAT* controlMPC_dzcc03 = controlMPC_dz_cc + 6;
controlMPC_FLOAT* controlMPC_rd03 = controlMPC_rd + 6;
controlMPC_FLOAT controlMPC_Lbyrd03[2];
controlMPC_FLOAT* controlMPC_grad_cost03 = controlMPC_grad_cost + 6;
controlMPC_FLOAT* controlMPC_grad_eq03 = controlMPC_grad_eq + 6;
controlMPC_FLOAT* controlMPC_grad_ineq03 = controlMPC_grad_ineq + 6;
controlMPC_FLOAT controlMPC_ctv03[2];
controlMPC_FLOAT* controlMPC_v03 = controlMPC_v + 6;
controlMPC_FLOAT controlMPC_re03[2];
controlMPC_FLOAT controlMPC_beta03[2];
controlMPC_FLOAT controlMPC_betacc03[2];
controlMPC_FLOAT* controlMPC_dvaff03 = controlMPC_dv_aff + 6;
controlMPC_FLOAT* controlMPC_dvcc03 = controlMPC_dv_cc + 6;
controlMPC_FLOAT controlMPC_V03[4];
controlMPC_FLOAT controlMPC_Yd03[3];
controlMPC_FLOAT controlMPC_Ld03[3];
controlMPC_FLOAT controlMPC_yy03[2];
controlMPC_FLOAT controlMPC_bmy03[2];
controlMPC_FLOAT controlMPC_c03[2] = {0.0000000000000000E+000, 0.0000000000000000E+000};
int controlMPC_lbIdx03[2] = {0, 1};
controlMPC_FLOAT* controlMPC_llb03 = controlMPC_l + 12;
controlMPC_FLOAT* controlMPC_slb03 = controlMPC_s + 12;
controlMPC_FLOAT* controlMPC_llbbyslb03 = controlMPC_lbys + 12;
controlMPC_FLOAT controlMPC_rilb03[2];
controlMPC_FLOAT* controlMPC_dllbaff03 = controlMPC_dl_aff + 12;
controlMPC_FLOAT* controlMPC_dslbaff03 = controlMPC_ds_aff + 12;
controlMPC_FLOAT* controlMPC_dllbcc03 = controlMPC_dl_cc + 12;
controlMPC_FLOAT* controlMPC_dslbcc03 = controlMPC_ds_cc + 12;
controlMPC_FLOAT* controlMPC_ccrhsl03 = controlMPC_ccrhs + 12;
int controlMPC_ubIdx03[2] = {0, 1};
controlMPC_FLOAT* controlMPC_lub03 = controlMPC_l + 14;
controlMPC_FLOAT* controlMPC_sub03 = controlMPC_s + 14;
controlMPC_FLOAT* controlMPC_lubbysub03 = controlMPC_lbys + 14;
controlMPC_FLOAT controlMPC_riub03[2];
controlMPC_FLOAT* controlMPC_dlubaff03 = controlMPC_dl_aff + 14;
controlMPC_FLOAT* controlMPC_dsubaff03 = controlMPC_ds_aff + 14;
controlMPC_FLOAT* controlMPC_dlubcc03 = controlMPC_dl_cc + 14;
controlMPC_FLOAT* controlMPC_dsubcc03 = controlMPC_ds_cc + 14;
controlMPC_FLOAT* controlMPC_ccrhsub03 = controlMPC_ccrhs + 14;
controlMPC_FLOAT controlMPC_Phi03[2];
controlMPC_FLOAT controlMPC_W03[2];
controlMPC_FLOAT controlMPC_Ysd03[4];
controlMPC_FLOAT controlMPC_Lsd03[4];
controlMPC_FLOAT* controlMPC_z04 = controlMPC_z + 8;
controlMPC_FLOAT* controlMPC_dzaff04 = controlMPC_dz_aff + 8;
controlMPC_FLOAT* controlMPC_dzcc04 = controlMPC_dz_cc + 8;
controlMPC_FLOAT* controlMPC_rd04 = controlMPC_rd + 8;
controlMPC_FLOAT controlMPC_Lbyrd04[2];
controlMPC_FLOAT* controlMPC_grad_cost04 = controlMPC_grad_cost + 8;
controlMPC_FLOAT* controlMPC_grad_eq04 = controlMPC_grad_eq + 8;
controlMPC_FLOAT* controlMPC_grad_ineq04 = controlMPC_grad_ineq + 8;
controlMPC_FLOAT controlMPC_ctv04[2];
controlMPC_FLOAT* controlMPC_v04 = controlMPC_v + 8;
controlMPC_FLOAT controlMPC_re04[2];
controlMPC_FLOAT controlMPC_beta04[2];
controlMPC_FLOAT controlMPC_betacc04[2];
controlMPC_FLOAT* controlMPC_dvaff04 = controlMPC_dv_aff + 8;
controlMPC_FLOAT* controlMPC_dvcc04 = controlMPC_dv_cc + 8;
controlMPC_FLOAT controlMPC_V04[4];
controlMPC_FLOAT controlMPC_Yd04[3];
controlMPC_FLOAT controlMPC_Ld04[3];
controlMPC_FLOAT controlMPC_yy04[2];
controlMPC_FLOAT controlMPC_bmy04[2];
controlMPC_FLOAT controlMPC_c04[2] = {0.0000000000000000E+000, 0.0000000000000000E+000};
int controlMPC_lbIdx04[2] = {0, 1};
controlMPC_FLOAT* controlMPC_llb04 = controlMPC_l + 16;
controlMPC_FLOAT* controlMPC_slb04 = controlMPC_s + 16;
controlMPC_FLOAT* controlMPC_llbbyslb04 = controlMPC_lbys + 16;
controlMPC_FLOAT controlMPC_rilb04[2];
controlMPC_FLOAT* controlMPC_dllbaff04 = controlMPC_dl_aff + 16;
controlMPC_FLOAT* controlMPC_dslbaff04 = controlMPC_ds_aff + 16;
controlMPC_FLOAT* controlMPC_dllbcc04 = controlMPC_dl_cc + 16;
controlMPC_FLOAT* controlMPC_dslbcc04 = controlMPC_ds_cc + 16;
controlMPC_FLOAT* controlMPC_ccrhsl04 = controlMPC_ccrhs + 16;
int controlMPC_ubIdx04[2] = {0, 1};
controlMPC_FLOAT* controlMPC_lub04 = controlMPC_l + 18;
controlMPC_FLOAT* controlMPC_sub04 = controlMPC_s + 18;
controlMPC_FLOAT* controlMPC_lubbysub04 = controlMPC_lbys + 18;
controlMPC_FLOAT controlMPC_riub04[2];
controlMPC_FLOAT* controlMPC_dlubaff04 = controlMPC_dl_aff + 18;
controlMPC_FLOAT* controlMPC_dsubaff04 = controlMPC_ds_aff + 18;
controlMPC_FLOAT* controlMPC_dlubcc04 = controlMPC_dl_cc + 18;
controlMPC_FLOAT* controlMPC_dsubcc04 = controlMPC_ds_cc + 18;
controlMPC_FLOAT* controlMPC_ccrhsub04 = controlMPC_ccrhs + 18;
controlMPC_FLOAT controlMPC_Phi04[2];
controlMPC_FLOAT controlMPC_W04[2];
controlMPC_FLOAT controlMPC_Ysd04[4];
controlMPC_FLOAT controlMPC_Lsd04[4];
controlMPC_FLOAT* controlMPC_z05 = controlMPC_z + 10;
controlMPC_FLOAT* controlMPC_dzaff05 = controlMPC_dz_aff + 10;
controlMPC_FLOAT* controlMPC_dzcc05 = controlMPC_dz_cc + 10;
controlMPC_FLOAT* controlMPC_rd05 = controlMPC_rd + 10;
controlMPC_FLOAT controlMPC_Lbyrd05[2];
controlMPC_FLOAT* controlMPC_grad_cost05 = controlMPC_grad_cost + 10;
controlMPC_FLOAT* controlMPC_grad_eq05 = controlMPC_grad_eq + 10;
controlMPC_FLOAT* controlMPC_grad_ineq05 = controlMPC_grad_ineq + 10;
controlMPC_FLOAT controlMPC_ctv05[2];
controlMPC_FLOAT* controlMPC_v05 = controlMPC_v + 10;
controlMPC_FLOAT controlMPC_re05[2];
controlMPC_FLOAT controlMPC_beta05[2];
controlMPC_FLOAT controlMPC_betacc05[2];
controlMPC_FLOAT* controlMPC_dvaff05 = controlMPC_dv_aff + 10;
controlMPC_FLOAT* controlMPC_dvcc05 = controlMPC_dv_cc + 10;
controlMPC_FLOAT controlMPC_V05[4];
controlMPC_FLOAT controlMPC_Yd05[3];
controlMPC_FLOAT controlMPC_Ld05[3];
controlMPC_FLOAT controlMPC_yy05[2];
controlMPC_FLOAT controlMPC_bmy05[2];
controlMPC_FLOAT controlMPC_c05[2] = {0.0000000000000000E+000, 0.0000000000000000E+000};
int controlMPC_lbIdx05[2] = {0, 1};
controlMPC_FLOAT* controlMPC_llb05 = controlMPC_l + 20;
controlMPC_FLOAT* controlMPC_slb05 = controlMPC_s + 20;
controlMPC_FLOAT* controlMPC_llbbyslb05 = controlMPC_lbys + 20;
controlMPC_FLOAT controlMPC_rilb05[2];
controlMPC_FLOAT* controlMPC_dllbaff05 = controlMPC_dl_aff + 20;
controlMPC_FLOAT* controlMPC_dslbaff05 = controlMPC_ds_aff + 20;
controlMPC_FLOAT* controlMPC_dllbcc05 = controlMPC_dl_cc + 20;
controlMPC_FLOAT* controlMPC_dslbcc05 = controlMPC_ds_cc + 20;
controlMPC_FLOAT* controlMPC_ccrhsl05 = controlMPC_ccrhs + 20;
int controlMPC_ubIdx05[2] = {0, 1};
controlMPC_FLOAT* controlMPC_lub05 = controlMPC_l + 22;
controlMPC_FLOAT* controlMPC_sub05 = controlMPC_s + 22;
controlMPC_FLOAT* controlMPC_lubbysub05 = controlMPC_lbys + 22;
controlMPC_FLOAT controlMPC_riub05[2];
controlMPC_FLOAT* controlMPC_dlubaff05 = controlMPC_dl_aff + 22;
controlMPC_FLOAT* controlMPC_dsubaff05 = controlMPC_ds_aff + 22;
controlMPC_FLOAT* controlMPC_dlubcc05 = controlMPC_dl_cc + 22;
controlMPC_FLOAT* controlMPC_dsubcc05 = controlMPC_ds_cc + 22;
controlMPC_FLOAT* controlMPC_ccrhsub05 = controlMPC_ccrhs + 22;
controlMPC_FLOAT controlMPC_Phi05[2];
controlMPC_FLOAT controlMPC_W05[2];
controlMPC_FLOAT controlMPC_Ysd05[4];
controlMPC_FLOAT controlMPC_Lsd05[4];
controlMPC_FLOAT* controlMPC_z06 = controlMPC_z + 12;
controlMPC_FLOAT* controlMPC_dzaff06 = controlMPC_dz_aff + 12;
controlMPC_FLOAT* controlMPC_dzcc06 = controlMPC_dz_cc + 12;
controlMPC_FLOAT* controlMPC_rd06 = controlMPC_rd + 12;
controlMPC_FLOAT controlMPC_Lbyrd06[2];
controlMPC_FLOAT* controlMPC_grad_cost06 = controlMPC_grad_cost + 12;
controlMPC_FLOAT* controlMPC_grad_eq06 = controlMPC_grad_eq + 12;
controlMPC_FLOAT* controlMPC_grad_ineq06 = controlMPC_grad_ineq + 12;
controlMPC_FLOAT controlMPC_ctv06[2];
controlMPC_FLOAT* controlMPC_v06 = controlMPC_v + 12;
controlMPC_FLOAT controlMPC_re06[2];
controlMPC_FLOAT controlMPC_beta06[2];
controlMPC_FLOAT controlMPC_betacc06[2];
controlMPC_FLOAT* controlMPC_dvaff06 = controlMPC_dv_aff + 12;
controlMPC_FLOAT* controlMPC_dvcc06 = controlMPC_dv_cc + 12;
controlMPC_FLOAT controlMPC_V06[4];
controlMPC_FLOAT controlMPC_Yd06[3];
controlMPC_FLOAT controlMPC_Ld06[3];
controlMPC_FLOAT controlMPC_yy06[2];
controlMPC_FLOAT controlMPC_bmy06[2];
controlMPC_FLOAT controlMPC_c06[2] = {0.0000000000000000E+000, 0.0000000000000000E+000};
int controlMPC_lbIdx06[2] = {0, 1};
controlMPC_FLOAT* controlMPC_llb06 = controlMPC_l + 24;
controlMPC_FLOAT* controlMPC_slb06 = controlMPC_s + 24;
controlMPC_FLOAT* controlMPC_llbbyslb06 = controlMPC_lbys + 24;
controlMPC_FLOAT controlMPC_rilb06[2];
controlMPC_FLOAT* controlMPC_dllbaff06 = controlMPC_dl_aff + 24;
controlMPC_FLOAT* controlMPC_dslbaff06 = controlMPC_ds_aff + 24;
controlMPC_FLOAT* controlMPC_dllbcc06 = controlMPC_dl_cc + 24;
controlMPC_FLOAT* controlMPC_dslbcc06 = controlMPC_ds_cc + 24;
controlMPC_FLOAT* controlMPC_ccrhsl06 = controlMPC_ccrhs + 24;
int controlMPC_ubIdx06[2] = {0, 1};
controlMPC_FLOAT* controlMPC_lub06 = controlMPC_l + 26;
controlMPC_FLOAT* controlMPC_sub06 = controlMPC_s + 26;
controlMPC_FLOAT* controlMPC_lubbysub06 = controlMPC_lbys + 26;
controlMPC_FLOAT controlMPC_riub06[2];
controlMPC_FLOAT* controlMPC_dlubaff06 = controlMPC_dl_aff + 26;
controlMPC_FLOAT* controlMPC_dsubaff06 = controlMPC_ds_aff + 26;
controlMPC_FLOAT* controlMPC_dlubcc06 = controlMPC_dl_cc + 26;
controlMPC_FLOAT* controlMPC_dsubcc06 = controlMPC_ds_cc + 26;
controlMPC_FLOAT* controlMPC_ccrhsub06 = controlMPC_ccrhs + 26;
controlMPC_FLOAT controlMPC_Phi06[2];
controlMPC_FLOAT controlMPC_W06[2];
controlMPC_FLOAT controlMPC_Ysd06[4];
controlMPC_FLOAT controlMPC_Lsd06[4];
controlMPC_FLOAT* controlMPC_z07 = controlMPC_z + 14;
controlMPC_FLOAT* controlMPC_dzaff07 = controlMPC_dz_aff + 14;
controlMPC_FLOAT* controlMPC_dzcc07 = controlMPC_dz_cc + 14;
controlMPC_FLOAT* controlMPC_rd07 = controlMPC_rd + 14;
controlMPC_FLOAT controlMPC_Lbyrd07[2];
controlMPC_FLOAT* controlMPC_grad_cost07 = controlMPC_grad_cost + 14;
controlMPC_FLOAT* controlMPC_grad_eq07 = controlMPC_grad_eq + 14;
controlMPC_FLOAT* controlMPC_grad_ineq07 = controlMPC_grad_ineq + 14;
controlMPC_FLOAT controlMPC_ctv07[2];
controlMPC_FLOAT* controlMPC_v07 = controlMPC_v + 14;
controlMPC_FLOAT controlMPC_re07[2];
controlMPC_FLOAT controlMPC_beta07[2];
controlMPC_FLOAT controlMPC_betacc07[2];
controlMPC_FLOAT* controlMPC_dvaff07 = controlMPC_dv_aff + 14;
controlMPC_FLOAT* controlMPC_dvcc07 = controlMPC_dv_cc + 14;
controlMPC_FLOAT controlMPC_V07[4];
controlMPC_FLOAT controlMPC_Yd07[3];
controlMPC_FLOAT controlMPC_Ld07[3];
controlMPC_FLOAT controlMPC_yy07[2];
controlMPC_FLOAT controlMPC_bmy07[2];
controlMPC_FLOAT controlMPC_c07[2] = {0.0000000000000000E+000, 0.0000000000000000E+000};
int controlMPC_lbIdx07[2] = {0, 1};
controlMPC_FLOAT* controlMPC_llb07 = controlMPC_l + 28;
controlMPC_FLOAT* controlMPC_slb07 = controlMPC_s + 28;
controlMPC_FLOAT* controlMPC_llbbyslb07 = controlMPC_lbys + 28;
controlMPC_FLOAT controlMPC_rilb07[2];
controlMPC_FLOAT* controlMPC_dllbaff07 = controlMPC_dl_aff + 28;
controlMPC_FLOAT* controlMPC_dslbaff07 = controlMPC_ds_aff + 28;
controlMPC_FLOAT* controlMPC_dllbcc07 = controlMPC_dl_cc + 28;
controlMPC_FLOAT* controlMPC_dslbcc07 = controlMPC_ds_cc + 28;
controlMPC_FLOAT* controlMPC_ccrhsl07 = controlMPC_ccrhs + 28;
int controlMPC_ubIdx07[2] = {0, 1};
controlMPC_FLOAT* controlMPC_lub07 = controlMPC_l + 30;
controlMPC_FLOAT* controlMPC_sub07 = controlMPC_s + 30;
controlMPC_FLOAT* controlMPC_lubbysub07 = controlMPC_lbys + 30;
controlMPC_FLOAT controlMPC_riub07[2];
controlMPC_FLOAT* controlMPC_dlubaff07 = controlMPC_dl_aff + 30;
controlMPC_FLOAT* controlMPC_dsubaff07 = controlMPC_ds_aff + 30;
controlMPC_FLOAT* controlMPC_dlubcc07 = controlMPC_dl_cc + 30;
controlMPC_FLOAT* controlMPC_dsubcc07 = controlMPC_ds_cc + 30;
controlMPC_FLOAT* controlMPC_ccrhsub07 = controlMPC_ccrhs + 30;
controlMPC_FLOAT controlMPC_Phi07[2];
controlMPC_FLOAT controlMPC_W07[2];
controlMPC_FLOAT controlMPC_Ysd07[4];
controlMPC_FLOAT controlMPC_Lsd07[4];
controlMPC_FLOAT* controlMPC_z08 = controlMPC_z + 16;
controlMPC_FLOAT* controlMPC_dzaff08 = controlMPC_dz_aff + 16;
controlMPC_FLOAT* controlMPC_dzcc08 = controlMPC_dz_cc + 16;
controlMPC_FLOAT* controlMPC_rd08 = controlMPC_rd + 16;
controlMPC_FLOAT controlMPC_Lbyrd08[2];
controlMPC_FLOAT* controlMPC_grad_cost08 = controlMPC_grad_cost + 16;
controlMPC_FLOAT* controlMPC_grad_eq08 = controlMPC_grad_eq + 16;
controlMPC_FLOAT* controlMPC_grad_ineq08 = controlMPC_grad_ineq + 16;
controlMPC_FLOAT controlMPC_ctv08[2];
controlMPC_FLOAT* controlMPC_v08 = controlMPC_v + 16;
controlMPC_FLOAT controlMPC_re08[2];
controlMPC_FLOAT controlMPC_beta08[2];
controlMPC_FLOAT controlMPC_betacc08[2];
controlMPC_FLOAT* controlMPC_dvaff08 = controlMPC_dv_aff + 16;
controlMPC_FLOAT* controlMPC_dvcc08 = controlMPC_dv_cc + 16;
controlMPC_FLOAT controlMPC_V08[4];
controlMPC_FLOAT controlMPC_Yd08[3];
controlMPC_FLOAT controlMPC_Ld08[3];
controlMPC_FLOAT controlMPC_yy08[2];
controlMPC_FLOAT controlMPC_bmy08[2];
controlMPC_FLOAT controlMPC_c08[2] = {0.0000000000000000E+000, 0.0000000000000000E+000};
int controlMPC_lbIdx08[2] = {0, 1};
controlMPC_FLOAT* controlMPC_llb08 = controlMPC_l + 32;
controlMPC_FLOAT* controlMPC_slb08 = controlMPC_s + 32;
controlMPC_FLOAT* controlMPC_llbbyslb08 = controlMPC_lbys + 32;
controlMPC_FLOAT controlMPC_rilb08[2];
controlMPC_FLOAT* controlMPC_dllbaff08 = controlMPC_dl_aff + 32;
controlMPC_FLOAT* controlMPC_dslbaff08 = controlMPC_ds_aff + 32;
controlMPC_FLOAT* controlMPC_dllbcc08 = controlMPC_dl_cc + 32;
controlMPC_FLOAT* controlMPC_dslbcc08 = controlMPC_ds_cc + 32;
controlMPC_FLOAT* controlMPC_ccrhsl08 = controlMPC_ccrhs + 32;
int controlMPC_ubIdx08[2] = {0, 1};
controlMPC_FLOAT* controlMPC_lub08 = controlMPC_l + 34;
controlMPC_FLOAT* controlMPC_sub08 = controlMPC_s + 34;
controlMPC_FLOAT* controlMPC_lubbysub08 = controlMPC_lbys + 34;
controlMPC_FLOAT controlMPC_riub08[2];
controlMPC_FLOAT* controlMPC_dlubaff08 = controlMPC_dl_aff + 34;
controlMPC_FLOAT* controlMPC_dsubaff08 = controlMPC_ds_aff + 34;
controlMPC_FLOAT* controlMPC_dlubcc08 = controlMPC_dl_cc + 34;
controlMPC_FLOAT* controlMPC_dsubcc08 = controlMPC_ds_cc + 34;
controlMPC_FLOAT* controlMPC_ccrhsub08 = controlMPC_ccrhs + 34;
controlMPC_FLOAT controlMPC_Phi08[2];
controlMPC_FLOAT controlMPC_W08[2];
controlMPC_FLOAT controlMPC_Ysd08[4];
controlMPC_FLOAT controlMPC_Lsd08[4];
controlMPC_FLOAT* controlMPC_z09 = controlMPC_z + 18;
controlMPC_FLOAT* controlMPC_dzaff09 = controlMPC_dz_aff + 18;
controlMPC_FLOAT* controlMPC_dzcc09 = controlMPC_dz_cc + 18;
controlMPC_FLOAT* controlMPC_rd09 = controlMPC_rd + 18;
controlMPC_FLOAT controlMPC_Lbyrd09[2];
controlMPC_FLOAT* controlMPC_grad_cost09 = controlMPC_grad_cost + 18;
controlMPC_FLOAT* controlMPC_grad_eq09 = controlMPC_grad_eq + 18;
controlMPC_FLOAT* controlMPC_grad_ineq09 = controlMPC_grad_ineq + 18;
controlMPC_FLOAT controlMPC_ctv09[2];
controlMPC_FLOAT* controlMPC_v09 = controlMPC_v + 18;
controlMPC_FLOAT controlMPC_re09[2];
controlMPC_FLOAT controlMPC_beta09[2];
controlMPC_FLOAT controlMPC_betacc09[2];
controlMPC_FLOAT* controlMPC_dvaff09 = controlMPC_dv_aff + 18;
controlMPC_FLOAT* controlMPC_dvcc09 = controlMPC_dv_cc + 18;
controlMPC_FLOAT controlMPC_V09[4];
controlMPC_FLOAT controlMPC_Yd09[3];
controlMPC_FLOAT controlMPC_Ld09[3];
controlMPC_FLOAT controlMPC_yy09[2];
controlMPC_FLOAT controlMPC_bmy09[2];
controlMPC_FLOAT controlMPC_c09[2] = {0.0000000000000000E+000, 0.0000000000000000E+000};
int controlMPC_lbIdx09[2] = {0, 1};
controlMPC_FLOAT* controlMPC_llb09 = controlMPC_l + 36;
controlMPC_FLOAT* controlMPC_slb09 = controlMPC_s + 36;
controlMPC_FLOAT* controlMPC_llbbyslb09 = controlMPC_lbys + 36;
controlMPC_FLOAT controlMPC_rilb09[2];
controlMPC_FLOAT* controlMPC_dllbaff09 = controlMPC_dl_aff + 36;
controlMPC_FLOAT* controlMPC_dslbaff09 = controlMPC_ds_aff + 36;
controlMPC_FLOAT* controlMPC_dllbcc09 = controlMPC_dl_cc + 36;
controlMPC_FLOAT* controlMPC_dslbcc09 = controlMPC_ds_cc + 36;
controlMPC_FLOAT* controlMPC_ccrhsl09 = controlMPC_ccrhs + 36;
int controlMPC_ubIdx09[2] = {0, 1};
controlMPC_FLOAT* controlMPC_lub09 = controlMPC_l + 38;
controlMPC_FLOAT* controlMPC_sub09 = controlMPC_s + 38;
controlMPC_FLOAT* controlMPC_lubbysub09 = controlMPC_lbys + 38;
controlMPC_FLOAT controlMPC_riub09[2];
controlMPC_FLOAT* controlMPC_dlubaff09 = controlMPC_dl_aff + 38;
controlMPC_FLOAT* controlMPC_dsubaff09 = controlMPC_ds_aff + 38;
controlMPC_FLOAT* controlMPC_dlubcc09 = controlMPC_dl_cc + 38;
controlMPC_FLOAT* controlMPC_dsubcc09 = controlMPC_ds_cc + 38;
controlMPC_FLOAT* controlMPC_ccrhsub09 = controlMPC_ccrhs + 38;
controlMPC_FLOAT controlMPC_Phi09[2];
controlMPC_FLOAT controlMPC_W09[2];
controlMPC_FLOAT controlMPC_Ysd09[4];
controlMPC_FLOAT controlMPC_Lsd09[4];
controlMPC_FLOAT* controlMPC_z10 = controlMPC_z + 20;
controlMPC_FLOAT* controlMPC_dzaff10 = controlMPC_dz_aff + 20;
controlMPC_FLOAT* controlMPC_dzcc10 = controlMPC_dz_cc + 20;
controlMPC_FLOAT* controlMPC_rd10 = controlMPC_rd + 20;
controlMPC_FLOAT controlMPC_Lbyrd10[2];
controlMPC_FLOAT* controlMPC_grad_cost10 = controlMPC_grad_cost + 20;
controlMPC_FLOAT* controlMPC_grad_eq10 = controlMPC_grad_eq + 20;
controlMPC_FLOAT* controlMPC_grad_ineq10 = controlMPC_grad_ineq + 20;
controlMPC_FLOAT controlMPC_ctv10[2];
controlMPC_FLOAT* controlMPC_v10 = controlMPC_v + 20;
controlMPC_FLOAT controlMPC_re10[2];
controlMPC_FLOAT controlMPC_beta10[2];
controlMPC_FLOAT controlMPC_betacc10[2];
controlMPC_FLOAT* controlMPC_dvaff10 = controlMPC_dv_aff + 20;
controlMPC_FLOAT* controlMPC_dvcc10 = controlMPC_dv_cc + 20;
controlMPC_FLOAT controlMPC_V10[4];
controlMPC_FLOAT controlMPC_Yd10[3];
controlMPC_FLOAT controlMPC_Ld10[3];
controlMPC_FLOAT controlMPC_yy10[2];
controlMPC_FLOAT controlMPC_bmy10[2];
controlMPC_FLOAT controlMPC_c10[2] = {0.0000000000000000E+000, 0.0000000000000000E+000};
int controlMPC_lbIdx10[2] = {0, 1};
controlMPC_FLOAT* controlMPC_llb10 = controlMPC_l + 40;
controlMPC_FLOAT* controlMPC_slb10 = controlMPC_s + 40;
controlMPC_FLOAT* controlMPC_llbbyslb10 = controlMPC_lbys + 40;
controlMPC_FLOAT controlMPC_rilb10[2];
controlMPC_FLOAT* controlMPC_dllbaff10 = controlMPC_dl_aff + 40;
controlMPC_FLOAT* controlMPC_dslbaff10 = controlMPC_ds_aff + 40;
controlMPC_FLOAT* controlMPC_dllbcc10 = controlMPC_dl_cc + 40;
controlMPC_FLOAT* controlMPC_dslbcc10 = controlMPC_ds_cc + 40;
controlMPC_FLOAT* controlMPC_ccrhsl10 = controlMPC_ccrhs + 40;
int controlMPC_ubIdx10[2] = {0, 1};
controlMPC_FLOAT* controlMPC_lub10 = controlMPC_l + 42;
controlMPC_FLOAT* controlMPC_sub10 = controlMPC_s + 42;
controlMPC_FLOAT* controlMPC_lubbysub10 = controlMPC_lbys + 42;
controlMPC_FLOAT controlMPC_riub10[2];
controlMPC_FLOAT* controlMPC_dlubaff10 = controlMPC_dl_aff + 42;
controlMPC_FLOAT* controlMPC_dsubaff10 = controlMPC_ds_aff + 42;
controlMPC_FLOAT* controlMPC_dlubcc10 = controlMPC_dl_cc + 42;
controlMPC_FLOAT* controlMPC_dsubcc10 = controlMPC_ds_cc + 42;
controlMPC_FLOAT* controlMPC_ccrhsub10 = controlMPC_ccrhs + 42;
controlMPC_FLOAT controlMPC_Phi10[2];
controlMPC_FLOAT controlMPC_W10[2];
controlMPC_FLOAT controlMPC_Ysd10[4];
controlMPC_FLOAT controlMPC_Lsd10[4];
controlMPC_FLOAT* controlMPC_z11 = controlMPC_z + 22;
controlMPC_FLOAT* controlMPC_dzaff11 = controlMPC_dz_aff + 22;
controlMPC_FLOAT* controlMPC_dzcc11 = controlMPC_dz_cc + 22;
controlMPC_FLOAT* controlMPC_rd11 = controlMPC_rd + 22;
controlMPC_FLOAT controlMPC_Lbyrd11[2];
controlMPC_FLOAT* controlMPC_grad_cost11 = controlMPC_grad_cost + 22;
controlMPC_FLOAT* controlMPC_grad_eq11 = controlMPC_grad_eq + 22;
controlMPC_FLOAT* controlMPC_grad_ineq11 = controlMPC_grad_ineq + 22;
controlMPC_FLOAT controlMPC_ctv11[2];
controlMPC_FLOAT* controlMPC_v11 = controlMPC_v + 22;
controlMPC_FLOAT controlMPC_re11[2];
controlMPC_FLOAT controlMPC_beta11[2];
controlMPC_FLOAT controlMPC_betacc11[2];
controlMPC_FLOAT* controlMPC_dvaff11 = controlMPC_dv_aff + 22;
controlMPC_FLOAT* controlMPC_dvcc11 = controlMPC_dv_cc + 22;
controlMPC_FLOAT controlMPC_V11[4];
controlMPC_FLOAT controlMPC_Yd11[3];
controlMPC_FLOAT controlMPC_Ld11[3];
controlMPC_FLOAT controlMPC_yy11[2];
controlMPC_FLOAT controlMPC_bmy11[2];
controlMPC_FLOAT controlMPC_c11[2] = {0.0000000000000000E+000, 0.0000000000000000E+000};
int controlMPC_lbIdx11[2] = {0, 1};
controlMPC_FLOAT* controlMPC_llb11 = controlMPC_l + 44;
controlMPC_FLOAT* controlMPC_slb11 = controlMPC_s + 44;
controlMPC_FLOAT* controlMPC_llbbyslb11 = controlMPC_lbys + 44;
controlMPC_FLOAT controlMPC_rilb11[2];
controlMPC_FLOAT* controlMPC_dllbaff11 = controlMPC_dl_aff + 44;
controlMPC_FLOAT* controlMPC_dslbaff11 = controlMPC_ds_aff + 44;
controlMPC_FLOAT* controlMPC_dllbcc11 = controlMPC_dl_cc + 44;
controlMPC_FLOAT* controlMPC_dslbcc11 = controlMPC_ds_cc + 44;
controlMPC_FLOAT* controlMPC_ccrhsl11 = controlMPC_ccrhs + 44;
int controlMPC_ubIdx11[2] = {0, 1};
controlMPC_FLOAT* controlMPC_lub11 = controlMPC_l + 46;
controlMPC_FLOAT* controlMPC_sub11 = controlMPC_s + 46;
controlMPC_FLOAT* controlMPC_lubbysub11 = controlMPC_lbys + 46;
controlMPC_FLOAT controlMPC_riub11[2];
controlMPC_FLOAT* controlMPC_dlubaff11 = controlMPC_dl_aff + 46;
controlMPC_FLOAT* controlMPC_dsubaff11 = controlMPC_ds_aff + 46;
controlMPC_FLOAT* controlMPC_dlubcc11 = controlMPC_dl_cc + 46;
controlMPC_FLOAT* controlMPC_dsubcc11 = controlMPC_ds_cc + 46;
controlMPC_FLOAT* controlMPC_ccrhsub11 = controlMPC_ccrhs + 46;
controlMPC_FLOAT controlMPC_Phi11[2];
controlMPC_FLOAT controlMPC_W11[2];
controlMPC_FLOAT controlMPC_Ysd11[4];
controlMPC_FLOAT controlMPC_Lsd11[4];
controlMPC_FLOAT* controlMPC_z12 = controlMPC_z + 24;
controlMPC_FLOAT* controlMPC_dzaff12 = controlMPC_dz_aff + 24;
controlMPC_FLOAT* controlMPC_dzcc12 = controlMPC_dz_cc + 24;
controlMPC_FLOAT* controlMPC_rd12 = controlMPC_rd + 24;
controlMPC_FLOAT controlMPC_Lbyrd12[2];
controlMPC_FLOAT* controlMPC_grad_cost12 = controlMPC_grad_cost + 24;
controlMPC_FLOAT* controlMPC_grad_eq12 = controlMPC_grad_eq + 24;
controlMPC_FLOAT* controlMPC_grad_ineq12 = controlMPC_grad_ineq + 24;
controlMPC_FLOAT controlMPC_ctv12[2];
controlMPC_FLOAT* controlMPC_v12 = controlMPC_v + 24;
controlMPC_FLOAT controlMPC_re12[2];
controlMPC_FLOAT controlMPC_beta12[2];
controlMPC_FLOAT controlMPC_betacc12[2];
controlMPC_FLOAT* controlMPC_dvaff12 = controlMPC_dv_aff + 24;
controlMPC_FLOAT* controlMPC_dvcc12 = controlMPC_dv_cc + 24;
controlMPC_FLOAT controlMPC_V12[4];
controlMPC_FLOAT controlMPC_Yd12[3];
controlMPC_FLOAT controlMPC_Ld12[3];
controlMPC_FLOAT controlMPC_yy12[2];
controlMPC_FLOAT controlMPC_bmy12[2];
controlMPC_FLOAT controlMPC_c12[2] = {0.0000000000000000E+000, 0.0000000000000000E+000};
int controlMPC_lbIdx12[2] = {0, 1};
controlMPC_FLOAT* controlMPC_llb12 = controlMPC_l + 48;
controlMPC_FLOAT* controlMPC_slb12 = controlMPC_s + 48;
controlMPC_FLOAT* controlMPC_llbbyslb12 = controlMPC_lbys + 48;
controlMPC_FLOAT controlMPC_rilb12[2];
controlMPC_FLOAT* controlMPC_dllbaff12 = controlMPC_dl_aff + 48;
controlMPC_FLOAT* controlMPC_dslbaff12 = controlMPC_ds_aff + 48;
controlMPC_FLOAT* controlMPC_dllbcc12 = controlMPC_dl_cc + 48;
controlMPC_FLOAT* controlMPC_dslbcc12 = controlMPC_ds_cc + 48;
controlMPC_FLOAT* controlMPC_ccrhsl12 = controlMPC_ccrhs + 48;
int controlMPC_ubIdx12[2] = {0, 1};
controlMPC_FLOAT* controlMPC_lub12 = controlMPC_l + 50;
controlMPC_FLOAT* controlMPC_sub12 = controlMPC_s + 50;
controlMPC_FLOAT* controlMPC_lubbysub12 = controlMPC_lbys + 50;
controlMPC_FLOAT controlMPC_riub12[2];
controlMPC_FLOAT* controlMPC_dlubaff12 = controlMPC_dl_aff + 50;
controlMPC_FLOAT* controlMPC_dsubaff12 = controlMPC_ds_aff + 50;
controlMPC_FLOAT* controlMPC_dlubcc12 = controlMPC_dl_cc + 50;
controlMPC_FLOAT* controlMPC_dsubcc12 = controlMPC_ds_cc + 50;
controlMPC_FLOAT* controlMPC_ccrhsub12 = controlMPC_ccrhs + 50;
controlMPC_FLOAT controlMPC_Phi12[2];
controlMPC_FLOAT controlMPC_W12[2];
controlMPC_FLOAT controlMPC_Ysd12[4];
controlMPC_FLOAT controlMPC_Lsd12[4];
controlMPC_FLOAT* controlMPC_z13 = controlMPC_z + 26;
controlMPC_FLOAT* controlMPC_dzaff13 = controlMPC_dz_aff + 26;
controlMPC_FLOAT* controlMPC_dzcc13 = controlMPC_dz_cc + 26;
controlMPC_FLOAT* controlMPC_rd13 = controlMPC_rd + 26;
controlMPC_FLOAT controlMPC_Lbyrd13[2];
controlMPC_FLOAT* controlMPC_grad_cost13 = controlMPC_grad_cost + 26;
controlMPC_FLOAT* controlMPC_grad_eq13 = controlMPC_grad_eq + 26;
controlMPC_FLOAT* controlMPC_grad_ineq13 = controlMPC_grad_ineq + 26;
controlMPC_FLOAT controlMPC_ctv13[2];
controlMPC_FLOAT* controlMPC_v13 = controlMPC_v + 26;
controlMPC_FLOAT controlMPC_re13[2];
controlMPC_FLOAT controlMPC_beta13[2];
controlMPC_FLOAT controlMPC_betacc13[2];
controlMPC_FLOAT* controlMPC_dvaff13 = controlMPC_dv_aff + 26;
controlMPC_FLOAT* controlMPC_dvcc13 = controlMPC_dv_cc + 26;
controlMPC_FLOAT controlMPC_V13[4];
controlMPC_FLOAT controlMPC_Yd13[3];
controlMPC_FLOAT controlMPC_Ld13[3];
controlMPC_FLOAT controlMPC_yy13[2];
controlMPC_FLOAT controlMPC_bmy13[2];
controlMPC_FLOAT controlMPC_c13[2] = {0.0000000000000000E+000, 0.0000000000000000E+000};
int controlMPC_lbIdx13[2] = {0, 1};
controlMPC_FLOAT* controlMPC_llb13 = controlMPC_l + 52;
controlMPC_FLOAT* controlMPC_slb13 = controlMPC_s + 52;
controlMPC_FLOAT* controlMPC_llbbyslb13 = controlMPC_lbys + 52;
controlMPC_FLOAT controlMPC_rilb13[2];
controlMPC_FLOAT* controlMPC_dllbaff13 = controlMPC_dl_aff + 52;
controlMPC_FLOAT* controlMPC_dslbaff13 = controlMPC_ds_aff + 52;
controlMPC_FLOAT* controlMPC_dllbcc13 = controlMPC_dl_cc + 52;
controlMPC_FLOAT* controlMPC_dslbcc13 = controlMPC_ds_cc + 52;
controlMPC_FLOAT* controlMPC_ccrhsl13 = controlMPC_ccrhs + 52;
int controlMPC_ubIdx13[2] = {0, 1};
controlMPC_FLOAT* controlMPC_lub13 = controlMPC_l + 54;
controlMPC_FLOAT* controlMPC_sub13 = controlMPC_s + 54;
controlMPC_FLOAT* controlMPC_lubbysub13 = controlMPC_lbys + 54;
controlMPC_FLOAT controlMPC_riub13[2];
controlMPC_FLOAT* controlMPC_dlubaff13 = controlMPC_dl_aff + 54;
controlMPC_FLOAT* controlMPC_dsubaff13 = controlMPC_ds_aff + 54;
controlMPC_FLOAT* controlMPC_dlubcc13 = controlMPC_dl_cc + 54;
controlMPC_FLOAT* controlMPC_dsubcc13 = controlMPC_ds_cc + 54;
controlMPC_FLOAT* controlMPC_ccrhsub13 = controlMPC_ccrhs + 54;
controlMPC_FLOAT controlMPC_Phi13[2];
controlMPC_FLOAT controlMPC_W13[2];
controlMPC_FLOAT controlMPC_Ysd13[4];
controlMPC_FLOAT controlMPC_Lsd13[4];
controlMPC_FLOAT* controlMPC_z14 = controlMPC_z + 28;
controlMPC_FLOAT* controlMPC_dzaff14 = controlMPC_dz_aff + 28;
controlMPC_FLOAT* controlMPC_dzcc14 = controlMPC_dz_cc + 28;
controlMPC_FLOAT* controlMPC_rd14 = controlMPC_rd + 28;
controlMPC_FLOAT controlMPC_Lbyrd14[2];
controlMPC_FLOAT* controlMPC_grad_cost14 = controlMPC_grad_cost + 28;
controlMPC_FLOAT* controlMPC_grad_eq14 = controlMPC_grad_eq + 28;
controlMPC_FLOAT* controlMPC_grad_ineq14 = controlMPC_grad_ineq + 28;
controlMPC_FLOAT controlMPC_ctv14[2];
controlMPC_FLOAT* controlMPC_v14 = controlMPC_v + 28;
controlMPC_FLOAT controlMPC_re14[2];
controlMPC_FLOAT controlMPC_beta14[2];
controlMPC_FLOAT controlMPC_betacc14[2];
controlMPC_FLOAT* controlMPC_dvaff14 = controlMPC_dv_aff + 28;
controlMPC_FLOAT* controlMPC_dvcc14 = controlMPC_dv_cc + 28;
controlMPC_FLOAT controlMPC_V14[4];
controlMPC_FLOAT controlMPC_Yd14[3];
controlMPC_FLOAT controlMPC_Ld14[3];
controlMPC_FLOAT controlMPC_yy14[2];
controlMPC_FLOAT controlMPC_bmy14[2];
controlMPC_FLOAT controlMPC_c14[2] = {0.0000000000000000E+000, 0.0000000000000000E+000};
int controlMPC_lbIdx14[2] = {0, 1};
controlMPC_FLOAT* controlMPC_llb14 = controlMPC_l + 56;
controlMPC_FLOAT* controlMPC_slb14 = controlMPC_s + 56;
controlMPC_FLOAT* controlMPC_llbbyslb14 = controlMPC_lbys + 56;
controlMPC_FLOAT controlMPC_rilb14[2];
controlMPC_FLOAT* controlMPC_dllbaff14 = controlMPC_dl_aff + 56;
controlMPC_FLOAT* controlMPC_dslbaff14 = controlMPC_ds_aff + 56;
controlMPC_FLOAT* controlMPC_dllbcc14 = controlMPC_dl_cc + 56;
controlMPC_FLOAT* controlMPC_dslbcc14 = controlMPC_ds_cc + 56;
controlMPC_FLOAT* controlMPC_ccrhsl14 = controlMPC_ccrhs + 56;
int controlMPC_ubIdx14[2] = {0, 1};
controlMPC_FLOAT* controlMPC_lub14 = controlMPC_l + 58;
controlMPC_FLOAT* controlMPC_sub14 = controlMPC_s + 58;
controlMPC_FLOAT* controlMPC_lubbysub14 = controlMPC_lbys + 58;
controlMPC_FLOAT controlMPC_riub14[2];
controlMPC_FLOAT* controlMPC_dlubaff14 = controlMPC_dl_aff + 58;
controlMPC_FLOAT* controlMPC_dsubaff14 = controlMPC_ds_aff + 58;
controlMPC_FLOAT* controlMPC_dlubcc14 = controlMPC_dl_cc + 58;
controlMPC_FLOAT* controlMPC_dsubcc14 = controlMPC_ds_cc + 58;
controlMPC_FLOAT* controlMPC_ccrhsub14 = controlMPC_ccrhs + 58;
controlMPC_FLOAT controlMPC_Phi14[2];
controlMPC_FLOAT controlMPC_W14[2];
controlMPC_FLOAT controlMPC_Ysd14[4];
controlMPC_FLOAT controlMPC_Lsd14[4];
controlMPC_FLOAT* controlMPC_z15 = controlMPC_z + 30;
controlMPC_FLOAT* controlMPC_dzaff15 = controlMPC_dz_aff + 30;
controlMPC_FLOAT* controlMPC_dzcc15 = controlMPC_dz_cc + 30;
controlMPC_FLOAT* controlMPC_rd15 = controlMPC_rd + 30;
controlMPC_FLOAT controlMPC_Lbyrd15[2];
controlMPC_FLOAT* controlMPC_grad_cost15 = controlMPC_grad_cost + 30;
controlMPC_FLOAT* controlMPC_grad_eq15 = controlMPC_grad_eq + 30;
controlMPC_FLOAT* controlMPC_grad_ineq15 = controlMPC_grad_ineq + 30;
controlMPC_FLOAT controlMPC_ctv15[2];
controlMPC_FLOAT* controlMPC_v15 = controlMPC_v + 30;
controlMPC_FLOAT controlMPC_re15[2];
controlMPC_FLOAT controlMPC_beta15[2];
controlMPC_FLOAT controlMPC_betacc15[2];
controlMPC_FLOAT* controlMPC_dvaff15 = controlMPC_dv_aff + 30;
controlMPC_FLOAT* controlMPC_dvcc15 = controlMPC_dv_cc + 30;
controlMPC_FLOAT controlMPC_V15[4];
controlMPC_FLOAT controlMPC_Yd15[3];
controlMPC_FLOAT controlMPC_Ld15[3];
controlMPC_FLOAT controlMPC_yy15[2];
controlMPC_FLOAT controlMPC_bmy15[2];
controlMPC_FLOAT controlMPC_c15[2] = {0.0000000000000000E+000, 0.0000000000000000E+000};
int controlMPC_lbIdx15[2] = {0, 1};
controlMPC_FLOAT* controlMPC_llb15 = controlMPC_l + 60;
controlMPC_FLOAT* controlMPC_slb15 = controlMPC_s + 60;
controlMPC_FLOAT* controlMPC_llbbyslb15 = controlMPC_lbys + 60;
controlMPC_FLOAT controlMPC_rilb15[2];
controlMPC_FLOAT* controlMPC_dllbaff15 = controlMPC_dl_aff + 60;
controlMPC_FLOAT* controlMPC_dslbaff15 = controlMPC_ds_aff + 60;
controlMPC_FLOAT* controlMPC_dllbcc15 = controlMPC_dl_cc + 60;
controlMPC_FLOAT* controlMPC_dslbcc15 = controlMPC_ds_cc + 60;
controlMPC_FLOAT* controlMPC_ccrhsl15 = controlMPC_ccrhs + 60;
int controlMPC_ubIdx15[2] = {0, 1};
controlMPC_FLOAT* controlMPC_lub15 = controlMPC_l + 62;
controlMPC_FLOAT* controlMPC_sub15 = controlMPC_s + 62;
controlMPC_FLOAT* controlMPC_lubbysub15 = controlMPC_lbys + 62;
controlMPC_FLOAT controlMPC_riub15[2];
controlMPC_FLOAT* controlMPC_dlubaff15 = controlMPC_dl_aff + 62;
controlMPC_FLOAT* controlMPC_dsubaff15 = controlMPC_ds_aff + 62;
controlMPC_FLOAT* controlMPC_dlubcc15 = controlMPC_dl_cc + 62;
controlMPC_FLOAT* controlMPC_dsubcc15 = controlMPC_ds_cc + 62;
controlMPC_FLOAT* controlMPC_ccrhsub15 = controlMPC_ccrhs + 62;
controlMPC_FLOAT controlMPC_Phi15[2];
controlMPC_FLOAT controlMPC_W15[2];
controlMPC_FLOAT controlMPC_Ysd15[4];
controlMPC_FLOAT controlMPC_Lsd15[4];
controlMPC_FLOAT* controlMPC_z16 = controlMPC_z + 32;
controlMPC_FLOAT* controlMPC_dzaff16 = controlMPC_dz_aff + 32;
controlMPC_FLOAT* controlMPC_dzcc16 = controlMPC_dz_cc + 32;
controlMPC_FLOAT* controlMPC_rd16 = controlMPC_rd + 32;
controlMPC_FLOAT controlMPC_Lbyrd16[2];
controlMPC_FLOAT* controlMPC_grad_cost16 = controlMPC_grad_cost + 32;
controlMPC_FLOAT* controlMPC_grad_eq16 = controlMPC_grad_eq + 32;
controlMPC_FLOAT* controlMPC_grad_ineq16 = controlMPC_grad_ineq + 32;
controlMPC_FLOAT controlMPC_ctv16[2];
controlMPC_FLOAT* controlMPC_v16 = controlMPC_v + 32;
controlMPC_FLOAT controlMPC_re16[2];
controlMPC_FLOAT controlMPC_beta16[2];
controlMPC_FLOAT controlMPC_betacc16[2];
controlMPC_FLOAT* controlMPC_dvaff16 = controlMPC_dv_aff + 32;
controlMPC_FLOAT* controlMPC_dvcc16 = controlMPC_dv_cc + 32;
controlMPC_FLOAT controlMPC_V16[4];
controlMPC_FLOAT controlMPC_Yd16[3];
controlMPC_FLOAT controlMPC_Ld16[3];
controlMPC_FLOAT controlMPC_yy16[2];
controlMPC_FLOAT controlMPC_bmy16[2];
controlMPC_FLOAT controlMPC_c16[2] = {0.0000000000000000E+000, 0.0000000000000000E+000};
int controlMPC_lbIdx16[2] = {0, 1};
controlMPC_FLOAT* controlMPC_llb16 = controlMPC_l + 64;
controlMPC_FLOAT* controlMPC_slb16 = controlMPC_s + 64;
controlMPC_FLOAT* controlMPC_llbbyslb16 = controlMPC_lbys + 64;
controlMPC_FLOAT controlMPC_rilb16[2];
controlMPC_FLOAT* controlMPC_dllbaff16 = controlMPC_dl_aff + 64;
controlMPC_FLOAT* controlMPC_dslbaff16 = controlMPC_ds_aff + 64;
controlMPC_FLOAT* controlMPC_dllbcc16 = controlMPC_dl_cc + 64;
controlMPC_FLOAT* controlMPC_dslbcc16 = controlMPC_ds_cc + 64;
controlMPC_FLOAT* controlMPC_ccrhsl16 = controlMPC_ccrhs + 64;
int controlMPC_ubIdx16[2] = {0, 1};
controlMPC_FLOAT* controlMPC_lub16 = controlMPC_l + 66;
controlMPC_FLOAT* controlMPC_sub16 = controlMPC_s + 66;
controlMPC_FLOAT* controlMPC_lubbysub16 = controlMPC_lbys + 66;
controlMPC_FLOAT controlMPC_riub16[2];
controlMPC_FLOAT* controlMPC_dlubaff16 = controlMPC_dl_aff + 66;
controlMPC_FLOAT* controlMPC_dsubaff16 = controlMPC_ds_aff + 66;
controlMPC_FLOAT* controlMPC_dlubcc16 = controlMPC_dl_cc + 66;
controlMPC_FLOAT* controlMPC_dsubcc16 = controlMPC_ds_cc + 66;
controlMPC_FLOAT* controlMPC_ccrhsub16 = controlMPC_ccrhs + 66;
controlMPC_FLOAT controlMPC_Phi16[2];
controlMPC_FLOAT controlMPC_W16[2];
controlMPC_FLOAT controlMPC_Ysd16[4];
controlMPC_FLOAT controlMPC_Lsd16[4];
controlMPC_FLOAT* controlMPC_z17 = controlMPC_z + 34;
controlMPC_FLOAT* controlMPC_dzaff17 = controlMPC_dz_aff + 34;
controlMPC_FLOAT* controlMPC_dzcc17 = controlMPC_dz_cc + 34;
controlMPC_FLOAT* controlMPC_rd17 = controlMPC_rd + 34;
controlMPC_FLOAT controlMPC_Lbyrd17[2];
controlMPC_FLOAT* controlMPC_grad_cost17 = controlMPC_grad_cost + 34;
controlMPC_FLOAT* controlMPC_grad_eq17 = controlMPC_grad_eq + 34;
controlMPC_FLOAT* controlMPC_grad_ineq17 = controlMPC_grad_ineq + 34;
controlMPC_FLOAT controlMPC_ctv17[2];
controlMPC_FLOAT* controlMPC_v17 = controlMPC_v + 34;
controlMPC_FLOAT controlMPC_re17[2];
controlMPC_FLOAT controlMPC_beta17[2];
controlMPC_FLOAT controlMPC_betacc17[2];
controlMPC_FLOAT* controlMPC_dvaff17 = controlMPC_dv_aff + 34;
controlMPC_FLOAT* controlMPC_dvcc17 = controlMPC_dv_cc + 34;
controlMPC_FLOAT controlMPC_V17[4];
controlMPC_FLOAT controlMPC_Yd17[3];
controlMPC_FLOAT controlMPC_Ld17[3];
controlMPC_FLOAT controlMPC_yy17[2];
controlMPC_FLOAT controlMPC_bmy17[2];
controlMPC_FLOAT controlMPC_c17[2] = {0.0000000000000000E+000, 0.0000000000000000E+000};
int controlMPC_lbIdx17[2] = {0, 1};
controlMPC_FLOAT* controlMPC_llb17 = controlMPC_l + 68;
controlMPC_FLOAT* controlMPC_slb17 = controlMPC_s + 68;
controlMPC_FLOAT* controlMPC_llbbyslb17 = controlMPC_lbys + 68;
controlMPC_FLOAT controlMPC_rilb17[2];
controlMPC_FLOAT* controlMPC_dllbaff17 = controlMPC_dl_aff + 68;
controlMPC_FLOAT* controlMPC_dslbaff17 = controlMPC_ds_aff + 68;
controlMPC_FLOAT* controlMPC_dllbcc17 = controlMPC_dl_cc + 68;
controlMPC_FLOAT* controlMPC_dslbcc17 = controlMPC_ds_cc + 68;
controlMPC_FLOAT* controlMPC_ccrhsl17 = controlMPC_ccrhs + 68;
int controlMPC_ubIdx17[2] = {0, 1};
controlMPC_FLOAT* controlMPC_lub17 = controlMPC_l + 70;
controlMPC_FLOAT* controlMPC_sub17 = controlMPC_s + 70;
controlMPC_FLOAT* controlMPC_lubbysub17 = controlMPC_lbys + 70;
controlMPC_FLOAT controlMPC_riub17[2];
controlMPC_FLOAT* controlMPC_dlubaff17 = controlMPC_dl_aff + 70;
controlMPC_FLOAT* controlMPC_dsubaff17 = controlMPC_ds_aff + 70;
controlMPC_FLOAT* controlMPC_dlubcc17 = controlMPC_dl_cc + 70;
controlMPC_FLOAT* controlMPC_dsubcc17 = controlMPC_ds_cc + 70;
controlMPC_FLOAT* controlMPC_ccrhsub17 = controlMPC_ccrhs + 70;
controlMPC_FLOAT controlMPC_Phi17[2];
controlMPC_FLOAT controlMPC_W17[2];
controlMPC_FLOAT controlMPC_Ysd17[4];
controlMPC_FLOAT controlMPC_Lsd17[4];
controlMPC_FLOAT* controlMPC_z18 = controlMPC_z + 36;
controlMPC_FLOAT* controlMPC_dzaff18 = controlMPC_dz_aff + 36;
controlMPC_FLOAT* controlMPC_dzcc18 = controlMPC_dz_cc + 36;
controlMPC_FLOAT* controlMPC_rd18 = controlMPC_rd + 36;
controlMPC_FLOAT controlMPC_Lbyrd18[2];
controlMPC_FLOAT* controlMPC_grad_cost18 = controlMPC_grad_cost + 36;
controlMPC_FLOAT* controlMPC_grad_eq18 = controlMPC_grad_eq + 36;
controlMPC_FLOAT* controlMPC_grad_ineq18 = controlMPC_grad_ineq + 36;
controlMPC_FLOAT controlMPC_ctv18[2];
controlMPC_FLOAT* controlMPC_v18 = controlMPC_v + 36;
controlMPC_FLOAT controlMPC_re18[2];
controlMPC_FLOAT controlMPC_beta18[2];
controlMPC_FLOAT controlMPC_betacc18[2];
controlMPC_FLOAT* controlMPC_dvaff18 = controlMPC_dv_aff + 36;
controlMPC_FLOAT* controlMPC_dvcc18 = controlMPC_dv_cc + 36;
controlMPC_FLOAT controlMPC_V18[4];
controlMPC_FLOAT controlMPC_Yd18[3];
controlMPC_FLOAT controlMPC_Ld18[3];
controlMPC_FLOAT controlMPC_yy18[2];
controlMPC_FLOAT controlMPC_bmy18[2];
controlMPC_FLOAT controlMPC_c18[2] = {0.0000000000000000E+000, 0.0000000000000000E+000};
int controlMPC_lbIdx18[2] = {0, 1};
controlMPC_FLOAT* controlMPC_llb18 = controlMPC_l + 72;
controlMPC_FLOAT* controlMPC_slb18 = controlMPC_s + 72;
controlMPC_FLOAT* controlMPC_llbbyslb18 = controlMPC_lbys + 72;
controlMPC_FLOAT controlMPC_rilb18[2];
controlMPC_FLOAT* controlMPC_dllbaff18 = controlMPC_dl_aff + 72;
controlMPC_FLOAT* controlMPC_dslbaff18 = controlMPC_ds_aff + 72;
controlMPC_FLOAT* controlMPC_dllbcc18 = controlMPC_dl_cc + 72;
controlMPC_FLOAT* controlMPC_dslbcc18 = controlMPC_ds_cc + 72;
controlMPC_FLOAT* controlMPC_ccrhsl18 = controlMPC_ccrhs + 72;
int controlMPC_ubIdx18[2] = {0, 1};
controlMPC_FLOAT* controlMPC_lub18 = controlMPC_l + 74;
controlMPC_FLOAT* controlMPC_sub18 = controlMPC_s + 74;
controlMPC_FLOAT* controlMPC_lubbysub18 = controlMPC_lbys + 74;
controlMPC_FLOAT controlMPC_riub18[2];
controlMPC_FLOAT* controlMPC_dlubaff18 = controlMPC_dl_aff + 74;
controlMPC_FLOAT* controlMPC_dsubaff18 = controlMPC_ds_aff + 74;
controlMPC_FLOAT* controlMPC_dlubcc18 = controlMPC_dl_cc + 74;
controlMPC_FLOAT* controlMPC_dsubcc18 = controlMPC_ds_cc + 74;
controlMPC_FLOAT* controlMPC_ccrhsub18 = controlMPC_ccrhs + 74;
controlMPC_FLOAT controlMPC_Phi18[2];
controlMPC_FLOAT controlMPC_W18[2];
controlMPC_FLOAT controlMPC_Ysd18[4];
controlMPC_FLOAT controlMPC_Lsd18[4];
controlMPC_FLOAT* controlMPC_z19 = controlMPC_z + 38;
controlMPC_FLOAT* controlMPC_dzaff19 = controlMPC_dz_aff + 38;
controlMPC_FLOAT* controlMPC_dzcc19 = controlMPC_dz_cc + 38;
controlMPC_FLOAT* controlMPC_rd19 = controlMPC_rd + 38;
controlMPC_FLOAT controlMPC_Lbyrd19[2];
controlMPC_FLOAT* controlMPC_grad_cost19 = controlMPC_grad_cost + 38;
controlMPC_FLOAT* controlMPC_grad_eq19 = controlMPC_grad_eq + 38;
controlMPC_FLOAT* controlMPC_grad_ineq19 = controlMPC_grad_ineq + 38;
controlMPC_FLOAT controlMPC_ctv19[2];
controlMPC_FLOAT* controlMPC_v19 = controlMPC_v + 38;
controlMPC_FLOAT controlMPC_re19[2];
controlMPC_FLOAT controlMPC_beta19[2];
controlMPC_FLOAT controlMPC_betacc19[2];
controlMPC_FLOAT* controlMPC_dvaff19 = controlMPC_dv_aff + 38;
controlMPC_FLOAT* controlMPC_dvcc19 = controlMPC_dv_cc + 38;
controlMPC_FLOAT controlMPC_V19[4];
controlMPC_FLOAT controlMPC_Yd19[3];
controlMPC_FLOAT controlMPC_Ld19[3];
controlMPC_FLOAT controlMPC_yy19[2];
controlMPC_FLOAT controlMPC_bmy19[2];
controlMPC_FLOAT controlMPC_c19[2] = {0.0000000000000000E+000, 0.0000000000000000E+000};
int controlMPC_lbIdx19[2] = {0, 1};
controlMPC_FLOAT* controlMPC_llb19 = controlMPC_l + 76;
controlMPC_FLOAT* controlMPC_slb19 = controlMPC_s + 76;
controlMPC_FLOAT* controlMPC_llbbyslb19 = controlMPC_lbys + 76;
controlMPC_FLOAT controlMPC_rilb19[2];
controlMPC_FLOAT* controlMPC_dllbaff19 = controlMPC_dl_aff + 76;
controlMPC_FLOAT* controlMPC_dslbaff19 = controlMPC_ds_aff + 76;
controlMPC_FLOAT* controlMPC_dllbcc19 = controlMPC_dl_cc + 76;
controlMPC_FLOAT* controlMPC_dslbcc19 = controlMPC_ds_cc + 76;
controlMPC_FLOAT* controlMPC_ccrhsl19 = controlMPC_ccrhs + 76;
int controlMPC_ubIdx19[2] = {0, 1};
controlMPC_FLOAT* controlMPC_lub19 = controlMPC_l + 78;
controlMPC_FLOAT* controlMPC_sub19 = controlMPC_s + 78;
controlMPC_FLOAT* controlMPC_lubbysub19 = controlMPC_lbys + 78;
controlMPC_FLOAT controlMPC_riub19[2];
controlMPC_FLOAT* controlMPC_dlubaff19 = controlMPC_dl_aff + 78;
controlMPC_FLOAT* controlMPC_dsubaff19 = controlMPC_ds_aff + 78;
controlMPC_FLOAT* controlMPC_dlubcc19 = controlMPC_dl_cc + 78;
controlMPC_FLOAT* controlMPC_dsubcc19 = controlMPC_ds_cc + 78;
controlMPC_FLOAT* controlMPC_ccrhsub19 = controlMPC_ccrhs + 78;
controlMPC_FLOAT controlMPC_Phi19[2];
controlMPC_FLOAT controlMPC_W19[2];
controlMPC_FLOAT controlMPC_Ysd19[4];
controlMPC_FLOAT controlMPC_Lsd19[4];
controlMPC_FLOAT* controlMPC_z20 = controlMPC_z + 40;
controlMPC_FLOAT* controlMPC_dzaff20 = controlMPC_dz_aff + 40;
controlMPC_FLOAT* controlMPC_dzcc20 = controlMPC_dz_cc + 40;
controlMPC_FLOAT* controlMPC_rd20 = controlMPC_rd + 40;
controlMPC_FLOAT controlMPC_Lbyrd20[2];
controlMPC_FLOAT* controlMPC_grad_cost20 = controlMPC_grad_cost + 40;
controlMPC_FLOAT* controlMPC_grad_eq20 = controlMPC_grad_eq + 40;
controlMPC_FLOAT* controlMPC_grad_ineq20 = controlMPC_grad_ineq + 40;
controlMPC_FLOAT controlMPC_ctv20[2];
controlMPC_FLOAT* controlMPC_v20 = controlMPC_v + 40;
controlMPC_FLOAT controlMPC_re20[2];
controlMPC_FLOAT controlMPC_beta20[2];
controlMPC_FLOAT controlMPC_betacc20[2];
controlMPC_FLOAT* controlMPC_dvaff20 = controlMPC_dv_aff + 40;
controlMPC_FLOAT* controlMPC_dvcc20 = controlMPC_dv_cc + 40;
controlMPC_FLOAT controlMPC_V20[4];
controlMPC_FLOAT controlMPC_Yd20[3];
controlMPC_FLOAT controlMPC_Ld20[3];
controlMPC_FLOAT controlMPC_yy20[2];
controlMPC_FLOAT controlMPC_bmy20[2];
controlMPC_FLOAT controlMPC_c20[2] = {0.0000000000000000E+000, 0.0000000000000000E+000};
int controlMPC_lbIdx20[2] = {0, 1};
controlMPC_FLOAT* controlMPC_llb20 = controlMPC_l + 80;
controlMPC_FLOAT* controlMPC_slb20 = controlMPC_s + 80;
controlMPC_FLOAT* controlMPC_llbbyslb20 = controlMPC_lbys + 80;
controlMPC_FLOAT controlMPC_rilb20[2];
controlMPC_FLOAT* controlMPC_dllbaff20 = controlMPC_dl_aff + 80;
controlMPC_FLOAT* controlMPC_dslbaff20 = controlMPC_ds_aff + 80;
controlMPC_FLOAT* controlMPC_dllbcc20 = controlMPC_dl_cc + 80;
controlMPC_FLOAT* controlMPC_dslbcc20 = controlMPC_ds_cc + 80;
controlMPC_FLOAT* controlMPC_ccrhsl20 = controlMPC_ccrhs + 80;
int controlMPC_ubIdx20[2] = {0, 1};
controlMPC_FLOAT* controlMPC_lub20 = controlMPC_l + 82;
controlMPC_FLOAT* controlMPC_sub20 = controlMPC_s + 82;
controlMPC_FLOAT* controlMPC_lubbysub20 = controlMPC_lbys + 82;
controlMPC_FLOAT controlMPC_riub20[2];
controlMPC_FLOAT* controlMPC_dlubaff20 = controlMPC_dl_aff + 82;
controlMPC_FLOAT* controlMPC_dsubaff20 = controlMPC_ds_aff + 82;
controlMPC_FLOAT* controlMPC_dlubcc20 = controlMPC_dl_cc + 82;
controlMPC_FLOAT* controlMPC_dsubcc20 = controlMPC_ds_cc + 82;
controlMPC_FLOAT* controlMPC_ccrhsub20 = controlMPC_ccrhs + 82;
controlMPC_FLOAT controlMPC_Phi20[2];
controlMPC_FLOAT controlMPC_W20[2];
controlMPC_FLOAT controlMPC_Ysd20[4];
controlMPC_FLOAT controlMPC_Lsd20[4];
controlMPC_FLOAT* controlMPC_z21 = controlMPC_z + 42;
controlMPC_FLOAT* controlMPC_dzaff21 = controlMPC_dz_aff + 42;
controlMPC_FLOAT* controlMPC_dzcc21 = controlMPC_dz_cc + 42;
controlMPC_FLOAT* controlMPC_rd21 = controlMPC_rd + 42;
controlMPC_FLOAT controlMPC_Lbyrd21[2];
controlMPC_FLOAT* controlMPC_grad_cost21 = controlMPC_grad_cost + 42;
controlMPC_FLOAT* controlMPC_grad_eq21 = controlMPC_grad_eq + 42;
controlMPC_FLOAT* controlMPC_grad_ineq21 = controlMPC_grad_ineq + 42;
controlMPC_FLOAT controlMPC_ctv21[2];
controlMPC_FLOAT* controlMPC_v21 = controlMPC_v + 42;
controlMPC_FLOAT controlMPC_re21[2];
controlMPC_FLOAT controlMPC_beta21[2];
controlMPC_FLOAT controlMPC_betacc21[2];
controlMPC_FLOAT* controlMPC_dvaff21 = controlMPC_dv_aff + 42;
controlMPC_FLOAT* controlMPC_dvcc21 = controlMPC_dv_cc + 42;
controlMPC_FLOAT controlMPC_V21[4];
controlMPC_FLOAT controlMPC_Yd21[3];
controlMPC_FLOAT controlMPC_Ld21[3];
controlMPC_FLOAT controlMPC_yy21[2];
controlMPC_FLOAT controlMPC_bmy21[2];
controlMPC_FLOAT controlMPC_c21[2] = {0.0000000000000000E+000, 0.0000000000000000E+000};
int controlMPC_lbIdx21[2] = {0, 1};
controlMPC_FLOAT* controlMPC_llb21 = controlMPC_l + 84;
controlMPC_FLOAT* controlMPC_slb21 = controlMPC_s + 84;
controlMPC_FLOAT* controlMPC_llbbyslb21 = controlMPC_lbys + 84;
controlMPC_FLOAT controlMPC_rilb21[2];
controlMPC_FLOAT* controlMPC_dllbaff21 = controlMPC_dl_aff + 84;
controlMPC_FLOAT* controlMPC_dslbaff21 = controlMPC_ds_aff + 84;
controlMPC_FLOAT* controlMPC_dllbcc21 = controlMPC_dl_cc + 84;
controlMPC_FLOAT* controlMPC_dslbcc21 = controlMPC_ds_cc + 84;
controlMPC_FLOAT* controlMPC_ccrhsl21 = controlMPC_ccrhs + 84;
int controlMPC_ubIdx21[2] = {0, 1};
controlMPC_FLOAT* controlMPC_lub21 = controlMPC_l + 86;
controlMPC_FLOAT* controlMPC_sub21 = controlMPC_s + 86;
controlMPC_FLOAT* controlMPC_lubbysub21 = controlMPC_lbys + 86;
controlMPC_FLOAT controlMPC_riub21[2];
controlMPC_FLOAT* controlMPC_dlubaff21 = controlMPC_dl_aff + 86;
controlMPC_FLOAT* controlMPC_dsubaff21 = controlMPC_ds_aff + 86;
controlMPC_FLOAT* controlMPC_dlubcc21 = controlMPC_dl_cc + 86;
controlMPC_FLOAT* controlMPC_dsubcc21 = controlMPC_ds_cc + 86;
controlMPC_FLOAT* controlMPC_ccrhsub21 = controlMPC_ccrhs + 86;
controlMPC_FLOAT controlMPC_Phi21[2];
controlMPC_FLOAT controlMPC_W21[2];
controlMPC_FLOAT controlMPC_Ysd21[4];
controlMPC_FLOAT controlMPC_Lsd21[4];
controlMPC_FLOAT* controlMPC_z22 = controlMPC_z + 44;
controlMPC_FLOAT* controlMPC_dzaff22 = controlMPC_dz_aff + 44;
controlMPC_FLOAT* controlMPC_dzcc22 = controlMPC_dz_cc + 44;
controlMPC_FLOAT* controlMPC_rd22 = controlMPC_rd + 44;
controlMPC_FLOAT controlMPC_Lbyrd22[2];
controlMPC_FLOAT* controlMPC_grad_cost22 = controlMPC_grad_cost + 44;
controlMPC_FLOAT* controlMPC_grad_eq22 = controlMPC_grad_eq + 44;
controlMPC_FLOAT* controlMPC_grad_ineq22 = controlMPC_grad_ineq + 44;
controlMPC_FLOAT controlMPC_ctv22[2];
controlMPC_FLOAT* controlMPC_v22 = controlMPC_v + 44;
controlMPC_FLOAT controlMPC_re22[2];
controlMPC_FLOAT controlMPC_beta22[2];
controlMPC_FLOAT controlMPC_betacc22[2];
controlMPC_FLOAT* controlMPC_dvaff22 = controlMPC_dv_aff + 44;
controlMPC_FLOAT* controlMPC_dvcc22 = controlMPC_dv_cc + 44;
controlMPC_FLOAT controlMPC_V22[4];
controlMPC_FLOAT controlMPC_Yd22[3];
controlMPC_FLOAT controlMPC_Ld22[3];
controlMPC_FLOAT controlMPC_yy22[2];
controlMPC_FLOAT controlMPC_bmy22[2];
controlMPC_FLOAT controlMPC_c22[2] = {0.0000000000000000E+000, 0.0000000000000000E+000};
int controlMPC_lbIdx22[2] = {0, 1};
controlMPC_FLOAT* controlMPC_llb22 = controlMPC_l + 88;
controlMPC_FLOAT* controlMPC_slb22 = controlMPC_s + 88;
controlMPC_FLOAT* controlMPC_llbbyslb22 = controlMPC_lbys + 88;
controlMPC_FLOAT controlMPC_rilb22[2];
controlMPC_FLOAT* controlMPC_dllbaff22 = controlMPC_dl_aff + 88;
controlMPC_FLOAT* controlMPC_dslbaff22 = controlMPC_ds_aff + 88;
controlMPC_FLOAT* controlMPC_dllbcc22 = controlMPC_dl_cc + 88;
controlMPC_FLOAT* controlMPC_dslbcc22 = controlMPC_ds_cc + 88;
controlMPC_FLOAT* controlMPC_ccrhsl22 = controlMPC_ccrhs + 88;
int controlMPC_ubIdx22[2] = {0, 1};
controlMPC_FLOAT* controlMPC_lub22 = controlMPC_l + 90;
controlMPC_FLOAT* controlMPC_sub22 = controlMPC_s + 90;
controlMPC_FLOAT* controlMPC_lubbysub22 = controlMPC_lbys + 90;
controlMPC_FLOAT controlMPC_riub22[2];
controlMPC_FLOAT* controlMPC_dlubaff22 = controlMPC_dl_aff + 90;
controlMPC_FLOAT* controlMPC_dsubaff22 = controlMPC_ds_aff + 90;
controlMPC_FLOAT* controlMPC_dlubcc22 = controlMPC_dl_cc + 90;
controlMPC_FLOAT* controlMPC_dsubcc22 = controlMPC_ds_cc + 90;
controlMPC_FLOAT* controlMPC_ccrhsub22 = controlMPC_ccrhs + 90;
controlMPC_FLOAT controlMPC_Phi22[2];
controlMPC_FLOAT controlMPC_W22[2];
controlMPC_FLOAT controlMPC_Ysd22[4];
controlMPC_FLOAT controlMPC_Lsd22[4];
controlMPC_FLOAT* controlMPC_z23 = controlMPC_z + 46;
controlMPC_FLOAT* controlMPC_dzaff23 = controlMPC_dz_aff + 46;
controlMPC_FLOAT* controlMPC_dzcc23 = controlMPC_dz_cc + 46;
controlMPC_FLOAT* controlMPC_rd23 = controlMPC_rd + 46;
controlMPC_FLOAT controlMPC_Lbyrd23[2];
controlMPC_FLOAT* controlMPC_grad_cost23 = controlMPC_grad_cost + 46;
controlMPC_FLOAT* controlMPC_grad_eq23 = controlMPC_grad_eq + 46;
controlMPC_FLOAT* controlMPC_grad_ineq23 = controlMPC_grad_ineq + 46;
controlMPC_FLOAT controlMPC_ctv23[2];
controlMPC_FLOAT* controlMPC_v23 = controlMPC_v + 46;
controlMPC_FLOAT controlMPC_re23[2];
controlMPC_FLOAT controlMPC_beta23[2];
controlMPC_FLOAT controlMPC_betacc23[2];
controlMPC_FLOAT* controlMPC_dvaff23 = controlMPC_dv_aff + 46;
controlMPC_FLOAT* controlMPC_dvcc23 = controlMPC_dv_cc + 46;
controlMPC_FLOAT controlMPC_V23[4];
controlMPC_FLOAT controlMPC_Yd23[3];
controlMPC_FLOAT controlMPC_Ld23[3];
controlMPC_FLOAT controlMPC_yy23[2];
controlMPC_FLOAT controlMPC_bmy23[2];
controlMPC_FLOAT controlMPC_c23[2] = {0.0000000000000000E+000, 0.0000000000000000E+000};
int controlMPC_lbIdx23[2] = {0, 1};
controlMPC_FLOAT* controlMPC_llb23 = controlMPC_l + 92;
controlMPC_FLOAT* controlMPC_slb23 = controlMPC_s + 92;
controlMPC_FLOAT* controlMPC_llbbyslb23 = controlMPC_lbys + 92;
controlMPC_FLOAT controlMPC_rilb23[2];
controlMPC_FLOAT* controlMPC_dllbaff23 = controlMPC_dl_aff + 92;
controlMPC_FLOAT* controlMPC_dslbaff23 = controlMPC_ds_aff + 92;
controlMPC_FLOAT* controlMPC_dllbcc23 = controlMPC_dl_cc + 92;
controlMPC_FLOAT* controlMPC_dslbcc23 = controlMPC_ds_cc + 92;
controlMPC_FLOAT* controlMPC_ccrhsl23 = controlMPC_ccrhs + 92;
int controlMPC_ubIdx23[2] = {0, 1};
controlMPC_FLOAT* controlMPC_lub23 = controlMPC_l + 94;
controlMPC_FLOAT* controlMPC_sub23 = controlMPC_s + 94;
controlMPC_FLOAT* controlMPC_lubbysub23 = controlMPC_lbys + 94;
controlMPC_FLOAT controlMPC_riub23[2];
controlMPC_FLOAT* controlMPC_dlubaff23 = controlMPC_dl_aff + 94;
controlMPC_FLOAT* controlMPC_dsubaff23 = controlMPC_ds_aff + 94;
controlMPC_FLOAT* controlMPC_dlubcc23 = controlMPC_dl_cc + 94;
controlMPC_FLOAT* controlMPC_dsubcc23 = controlMPC_ds_cc + 94;
controlMPC_FLOAT* controlMPC_ccrhsub23 = controlMPC_ccrhs + 94;
controlMPC_FLOAT controlMPC_Phi23[2];
controlMPC_FLOAT controlMPC_W23[2];
controlMPC_FLOAT controlMPC_Ysd23[4];
controlMPC_FLOAT controlMPC_Lsd23[4];
controlMPC_FLOAT* controlMPC_z24 = controlMPC_z + 48;
controlMPC_FLOAT* controlMPC_dzaff24 = controlMPC_dz_aff + 48;
controlMPC_FLOAT* controlMPC_dzcc24 = controlMPC_dz_cc + 48;
controlMPC_FLOAT* controlMPC_rd24 = controlMPC_rd + 48;
controlMPC_FLOAT controlMPC_Lbyrd24[2];
controlMPC_FLOAT* controlMPC_grad_cost24 = controlMPC_grad_cost + 48;
controlMPC_FLOAT* controlMPC_grad_eq24 = controlMPC_grad_eq + 48;
controlMPC_FLOAT* controlMPC_grad_ineq24 = controlMPC_grad_ineq + 48;
controlMPC_FLOAT controlMPC_ctv24[2];
controlMPC_FLOAT* controlMPC_v24 = controlMPC_v + 48;
controlMPC_FLOAT controlMPC_re24[2];
controlMPC_FLOAT controlMPC_beta24[2];
controlMPC_FLOAT controlMPC_betacc24[2];
controlMPC_FLOAT* controlMPC_dvaff24 = controlMPC_dv_aff + 48;
controlMPC_FLOAT* controlMPC_dvcc24 = controlMPC_dv_cc + 48;
controlMPC_FLOAT controlMPC_V24[4];
controlMPC_FLOAT controlMPC_Yd24[3];
controlMPC_FLOAT controlMPC_Ld24[3];
controlMPC_FLOAT controlMPC_yy24[2];
controlMPC_FLOAT controlMPC_bmy24[2];
controlMPC_FLOAT controlMPC_c24[2] = {0.0000000000000000E+000, 0.0000000000000000E+000};
int controlMPC_lbIdx24[2] = {0, 1};
controlMPC_FLOAT* controlMPC_llb24 = controlMPC_l + 96;
controlMPC_FLOAT* controlMPC_slb24 = controlMPC_s + 96;
controlMPC_FLOAT* controlMPC_llbbyslb24 = controlMPC_lbys + 96;
controlMPC_FLOAT controlMPC_rilb24[2];
controlMPC_FLOAT* controlMPC_dllbaff24 = controlMPC_dl_aff + 96;
controlMPC_FLOAT* controlMPC_dslbaff24 = controlMPC_ds_aff + 96;
controlMPC_FLOAT* controlMPC_dllbcc24 = controlMPC_dl_cc + 96;
controlMPC_FLOAT* controlMPC_dslbcc24 = controlMPC_ds_cc + 96;
controlMPC_FLOAT* controlMPC_ccrhsl24 = controlMPC_ccrhs + 96;
int controlMPC_ubIdx24[2] = {0, 1};
controlMPC_FLOAT* controlMPC_lub24 = controlMPC_l + 98;
controlMPC_FLOAT* controlMPC_sub24 = controlMPC_s + 98;
controlMPC_FLOAT* controlMPC_lubbysub24 = controlMPC_lbys + 98;
controlMPC_FLOAT controlMPC_riub24[2];
controlMPC_FLOAT* controlMPC_dlubaff24 = controlMPC_dl_aff + 98;
controlMPC_FLOAT* controlMPC_dsubaff24 = controlMPC_ds_aff + 98;
controlMPC_FLOAT* controlMPC_dlubcc24 = controlMPC_dl_cc + 98;
controlMPC_FLOAT* controlMPC_dsubcc24 = controlMPC_ds_cc + 98;
controlMPC_FLOAT* controlMPC_ccrhsub24 = controlMPC_ccrhs + 98;
controlMPC_FLOAT controlMPC_Phi24[2];
controlMPC_FLOAT controlMPC_W24[2];
controlMPC_FLOAT controlMPC_Ysd24[4];
controlMPC_FLOAT controlMPC_Lsd24[4];
controlMPC_FLOAT* controlMPC_z25 = controlMPC_z + 50;
controlMPC_FLOAT* controlMPC_dzaff25 = controlMPC_dz_aff + 50;
controlMPC_FLOAT* controlMPC_dzcc25 = controlMPC_dz_cc + 50;
controlMPC_FLOAT* controlMPC_rd25 = controlMPC_rd + 50;
controlMPC_FLOAT controlMPC_Lbyrd25[2];
controlMPC_FLOAT* controlMPC_grad_cost25 = controlMPC_grad_cost + 50;
controlMPC_FLOAT* controlMPC_grad_eq25 = controlMPC_grad_eq + 50;
controlMPC_FLOAT* controlMPC_grad_ineq25 = controlMPC_grad_ineq + 50;
controlMPC_FLOAT controlMPC_ctv25[2];
controlMPC_FLOAT* controlMPC_v25 = controlMPC_v + 50;
controlMPC_FLOAT controlMPC_re25[2];
controlMPC_FLOAT controlMPC_beta25[2];
controlMPC_FLOAT controlMPC_betacc25[2];
controlMPC_FLOAT* controlMPC_dvaff25 = controlMPC_dv_aff + 50;
controlMPC_FLOAT* controlMPC_dvcc25 = controlMPC_dv_cc + 50;
controlMPC_FLOAT controlMPC_V25[4];
controlMPC_FLOAT controlMPC_Yd25[3];
controlMPC_FLOAT controlMPC_Ld25[3];
controlMPC_FLOAT controlMPC_yy25[2];
controlMPC_FLOAT controlMPC_bmy25[2];
controlMPC_FLOAT controlMPC_c25[2] = {0.0000000000000000E+000, 0.0000000000000000E+000};
int controlMPC_lbIdx25[2] = {0, 1};
controlMPC_FLOAT* controlMPC_llb25 = controlMPC_l + 100;
controlMPC_FLOAT* controlMPC_slb25 = controlMPC_s + 100;
controlMPC_FLOAT* controlMPC_llbbyslb25 = controlMPC_lbys + 100;
controlMPC_FLOAT controlMPC_rilb25[2];
controlMPC_FLOAT* controlMPC_dllbaff25 = controlMPC_dl_aff + 100;
controlMPC_FLOAT* controlMPC_dslbaff25 = controlMPC_ds_aff + 100;
controlMPC_FLOAT* controlMPC_dllbcc25 = controlMPC_dl_cc + 100;
controlMPC_FLOAT* controlMPC_dslbcc25 = controlMPC_ds_cc + 100;
controlMPC_FLOAT* controlMPC_ccrhsl25 = controlMPC_ccrhs + 100;
int controlMPC_ubIdx25[2] = {0, 1};
controlMPC_FLOAT* controlMPC_lub25 = controlMPC_l + 102;
controlMPC_FLOAT* controlMPC_sub25 = controlMPC_s + 102;
controlMPC_FLOAT* controlMPC_lubbysub25 = controlMPC_lbys + 102;
controlMPC_FLOAT controlMPC_riub25[2];
controlMPC_FLOAT* controlMPC_dlubaff25 = controlMPC_dl_aff + 102;
controlMPC_FLOAT* controlMPC_dsubaff25 = controlMPC_ds_aff + 102;
controlMPC_FLOAT* controlMPC_dlubcc25 = controlMPC_dl_cc + 102;
controlMPC_FLOAT* controlMPC_dsubcc25 = controlMPC_ds_cc + 102;
controlMPC_FLOAT* controlMPC_ccrhsub25 = controlMPC_ccrhs + 102;
controlMPC_FLOAT controlMPC_Phi25[2];
controlMPC_FLOAT controlMPC_W25[2];
controlMPC_FLOAT controlMPC_Ysd25[4];
controlMPC_FLOAT controlMPC_Lsd25[4];
controlMPC_FLOAT* controlMPC_z26 = controlMPC_z + 52;
controlMPC_FLOAT* controlMPC_dzaff26 = controlMPC_dz_aff + 52;
controlMPC_FLOAT* controlMPC_dzcc26 = controlMPC_dz_cc + 52;
controlMPC_FLOAT* controlMPC_rd26 = controlMPC_rd + 52;
controlMPC_FLOAT controlMPC_Lbyrd26[2];
controlMPC_FLOAT* controlMPC_grad_cost26 = controlMPC_grad_cost + 52;
controlMPC_FLOAT* controlMPC_grad_eq26 = controlMPC_grad_eq + 52;
controlMPC_FLOAT* controlMPC_grad_ineq26 = controlMPC_grad_ineq + 52;
controlMPC_FLOAT controlMPC_ctv26[2];
controlMPC_FLOAT* controlMPC_v26 = controlMPC_v + 52;
controlMPC_FLOAT controlMPC_re26[2];
controlMPC_FLOAT controlMPC_beta26[2];
controlMPC_FLOAT controlMPC_betacc26[2];
controlMPC_FLOAT* controlMPC_dvaff26 = controlMPC_dv_aff + 52;
controlMPC_FLOAT* controlMPC_dvcc26 = controlMPC_dv_cc + 52;
controlMPC_FLOAT controlMPC_V26[4];
controlMPC_FLOAT controlMPC_Yd26[3];
controlMPC_FLOAT controlMPC_Ld26[3];
controlMPC_FLOAT controlMPC_yy26[2];
controlMPC_FLOAT controlMPC_bmy26[2];
controlMPC_FLOAT controlMPC_c26[2] = {0.0000000000000000E+000, 0.0000000000000000E+000};
int controlMPC_lbIdx26[2] = {0, 1};
controlMPC_FLOAT* controlMPC_llb26 = controlMPC_l + 104;
controlMPC_FLOAT* controlMPC_slb26 = controlMPC_s + 104;
controlMPC_FLOAT* controlMPC_llbbyslb26 = controlMPC_lbys + 104;
controlMPC_FLOAT controlMPC_rilb26[2];
controlMPC_FLOAT* controlMPC_dllbaff26 = controlMPC_dl_aff + 104;
controlMPC_FLOAT* controlMPC_dslbaff26 = controlMPC_ds_aff + 104;
controlMPC_FLOAT* controlMPC_dllbcc26 = controlMPC_dl_cc + 104;
controlMPC_FLOAT* controlMPC_dslbcc26 = controlMPC_ds_cc + 104;
controlMPC_FLOAT* controlMPC_ccrhsl26 = controlMPC_ccrhs + 104;
int controlMPC_ubIdx26[2] = {0, 1};
controlMPC_FLOAT* controlMPC_lub26 = controlMPC_l + 106;
controlMPC_FLOAT* controlMPC_sub26 = controlMPC_s + 106;
controlMPC_FLOAT* controlMPC_lubbysub26 = controlMPC_lbys + 106;
controlMPC_FLOAT controlMPC_riub26[2];
controlMPC_FLOAT* controlMPC_dlubaff26 = controlMPC_dl_aff + 106;
controlMPC_FLOAT* controlMPC_dsubaff26 = controlMPC_ds_aff + 106;
controlMPC_FLOAT* controlMPC_dlubcc26 = controlMPC_dl_cc + 106;
controlMPC_FLOAT* controlMPC_dsubcc26 = controlMPC_ds_cc + 106;
controlMPC_FLOAT* controlMPC_ccrhsub26 = controlMPC_ccrhs + 106;
controlMPC_FLOAT controlMPC_Phi26[2];
controlMPC_FLOAT controlMPC_W26[2];
controlMPC_FLOAT controlMPC_Ysd26[4];
controlMPC_FLOAT controlMPC_Lsd26[4];
controlMPC_FLOAT* controlMPC_z27 = controlMPC_z + 54;
controlMPC_FLOAT* controlMPC_dzaff27 = controlMPC_dz_aff + 54;
controlMPC_FLOAT* controlMPC_dzcc27 = controlMPC_dz_cc + 54;
controlMPC_FLOAT* controlMPC_rd27 = controlMPC_rd + 54;
controlMPC_FLOAT controlMPC_Lbyrd27[2];
controlMPC_FLOAT* controlMPC_grad_cost27 = controlMPC_grad_cost + 54;
controlMPC_FLOAT* controlMPC_grad_eq27 = controlMPC_grad_eq + 54;
controlMPC_FLOAT* controlMPC_grad_ineq27 = controlMPC_grad_ineq + 54;
controlMPC_FLOAT controlMPC_ctv27[2];
controlMPC_FLOAT* controlMPC_v27 = controlMPC_v + 54;
controlMPC_FLOAT controlMPC_re27[2];
controlMPC_FLOAT controlMPC_beta27[2];
controlMPC_FLOAT controlMPC_betacc27[2];
controlMPC_FLOAT* controlMPC_dvaff27 = controlMPC_dv_aff + 54;
controlMPC_FLOAT* controlMPC_dvcc27 = controlMPC_dv_cc + 54;
controlMPC_FLOAT controlMPC_V27[4];
controlMPC_FLOAT controlMPC_Yd27[3];
controlMPC_FLOAT controlMPC_Ld27[3];
controlMPC_FLOAT controlMPC_yy27[2];
controlMPC_FLOAT controlMPC_bmy27[2];
controlMPC_FLOAT controlMPC_c27[2] = {0.0000000000000000E+000, 0.0000000000000000E+000};
int controlMPC_lbIdx27[2] = {0, 1};
controlMPC_FLOAT* controlMPC_llb27 = controlMPC_l + 108;
controlMPC_FLOAT* controlMPC_slb27 = controlMPC_s + 108;
controlMPC_FLOAT* controlMPC_llbbyslb27 = controlMPC_lbys + 108;
controlMPC_FLOAT controlMPC_rilb27[2];
controlMPC_FLOAT* controlMPC_dllbaff27 = controlMPC_dl_aff + 108;
controlMPC_FLOAT* controlMPC_dslbaff27 = controlMPC_ds_aff + 108;
controlMPC_FLOAT* controlMPC_dllbcc27 = controlMPC_dl_cc + 108;
controlMPC_FLOAT* controlMPC_dslbcc27 = controlMPC_ds_cc + 108;
controlMPC_FLOAT* controlMPC_ccrhsl27 = controlMPC_ccrhs + 108;
int controlMPC_ubIdx27[2] = {0, 1};
controlMPC_FLOAT* controlMPC_lub27 = controlMPC_l + 110;
controlMPC_FLOAT* controlMPC_sub27 = controlMPC_s + 110;
controlMPC_FLOAT* controlMPC_lubbysub27 = controlMPC_lbys + 110;
controlMPC_FLOAT controlMPC_riub27[2];
controlMPC_FLOAT* controlMPC_dlubaff27 = controlMPC_dl_aff + 110;
controlMPC_FLOAT* controlMPC_dsubaff27 = controlMPC_ds_aff + 110;
controlMPC_FLOAT* controlMPC_dlubcc27 = controlMPC_dl_cc + 110;
controlMPC_FLOAT* controlMPC_dsubcc27 = controlMPC_ds_cc + 110;
controlMPC_FLOAT* controlMPC_ccrhsub27 = controlMPC_ccrhs + 110;
controlMPC_FLOAT controlMPC_Phi27[2];
controlMPC_FLOAT controlMPC_W27[2];
controlMPC_FLOAT controlMPC_Ysd27[4];
controlMPC_FLOAT controlMPC_Lsd27[4];
controlMPC_FLOAT* controlMPC_z28 = controlMPC_z + 56;
controlMPC_FLOAT* controlMPC_dzaff28 = controlMPC_dz_aff + 56;
controlMPC_FLOAT* controlMPC_dzcc28 = controlMPC_dz_cc + 56;
controlMPC_FLOAT* controlMPC_rd28 = controlMPC_rd + 56;
controlMPC_FLOAT controlMPC_Lbyrd28[2];
controlMPC_FLOAT* controlMPC_grad_cost28 = controlMPC_grad_cost + 56;
controlMPC_FLOAT* controlMPC_grad_eq28 = controlMPC_grad_eq + 56;
controlMPC_FLOAT* controlMPC_grad_ineq28 = controlMPC_grad_ineq + 56;
controlMPC_FLOAT controlMPC_ctv28[2];
int controlMPC_lbIdx28[2] = {0, 1};
controlMPC_FLOAT* controlMPC_llb28 = controlMPC_l + 112;
controlMPC_FLOAT* controlMPC_slb28 = controlMPC_s + 112;
controlMPC_FLOAT* controlMPC_llbbyslb28 = controlMPC_lbys + 112;
controlMPC_FLOAT controlMPC_rilb28[2];
controlMPC_FLOAT* controlMPC_dllbaff28 = controlMPC_dl_aff + 112;
controlMPC_FLOAT* controlMPC_dslbaff28 = controlMPC_ds_aff + 112;
controlMPC_FLOAT* controlMPC_dllbcc28 = controlMPC_dl_cc + 112;
controlMPC_FLOAT* controlMPC_dslbcc28 = controlMPC_ds_cc + 112;
controlMPC_FLOAT* controlMPC_ccrhsl28 = controlMPC_ccrhs + 112;
int controlMPC_ubIdx28[2] = {0, 1};
controlMPC_FLOAT* controlMPC_lub28 = controlMPC_l + 114;
controlMPC_FLOAT* controlMPC_sub28 = controlMPC_s + 114;
controlMPC_FLOAT* controlMPC_lubbysub28 = controlMPC_lbys + 114;
controlMPC_FLOAT controlMPC_riub28[2];
controlMPC_FLOAT* controlMPC_dlubaff28 = controlMPC_dl_aff + 114;
controlMPC_FLOAT* controlMPC_dsubaff28 = controlMPC_ds_aff + 114;
controlMPC_FLOAT* controlMPC_dlubcc28 = controlMPC_dl_cc + 114;
controlMPC_FLOAT* controlMPC_dsubcc28 = controlMPC_ds_cc + 114;
controlMPC_FLOAT* controlMPC_ccrhsub28 = controlMPC_ccrhs + 114;
controlMPC_FLOAT controlMPC_Phi28[2];
controlMPC_FLOAT controlMPC_W28[2];
controlMPC_FLOAT musigma;
controlMPC_FLOAT sigma_3rdroot;
controlMPC_FLOAT controlMPC_Diag1_0[2];
controlMPC_FLOAT controlMPC_Diag2_0[2];
controlMPC_FLOAT controlMPC_L_0[1];




/* SOLVER CODE --------------------------------------------------------- */
int controlMPC_solve(controlMPC_params* params, controlMPC_output* output, controlMPC_info* info)
{	
int exitcode;

#if controlMPC_SET_TIMING == 1
	controlMPC_timer solvertimer;
	controlMPC_tic(&solvertimer);
#endif
/* FUNCTION CALLS INTO LA LIBRARY -------------------------------------- */
info->it = 0;
controlMPC_LA_INITIALIZEVECTOR_58(controlMPC_z, 0);
controlMPC_LA_INITIALIZEVECTOR_56(controlMPC_v, 1);
controlMPC_LA_INITIALIZEVECTOR_116(controlMPC_l, 1);
controlMPC_LA_INITIALIZEVECTOR_116(controlMPC_s, 1);
info->mu = 0;
controlMPC_LA_DOTACC_116(controlMPC_l, controlMPC_s, &info->mu);
info->mu /= 116;
while( 1 ){
info->pobj = 0;
controlMPC_LA_DIAG_QUADFCN_2(params->H1, params->f1, controlMPC_z00, controlMPC_grad_cost00, &info->pobj);
controlMPC_LA_DIAG_QUADFCN_2(params->H2, params->f2, controlMPC_z01, controlMPC_grad_cost01, &info->pobj);
controlMPC_LA_DIAG_QUADFCN_2(params->H3, params->f3, controlMPC_z02, controlMPC_grad_cost02, &info->pobj);
controlMPC_LA_DIAG_QUADFCN_2(params->H4, params->f4, controlMPC_z03, controlMPC_grad_cost03, &info->pobj);
controlMPC_LA_DIAG_QUADFCN_2(params->H5, params->f5, controlMPC_z04, controlMPC_grad_cost04, &info->pobj);
controlMPC_LA_DIAG_QUADFCN_2(params->H6, params->f6, controlMPC_z05, controlMPC_grad_cost05, &info->pobj);
controlMPC_LA_DIAG_QUADFCN_2(params->H7, params->f7, controlMPC_z06, controlMPC_grad_cost06, &info->pobj);
controlMPC_LA_DIAG_QUADFCN_2(params->H8, params->f8, controlMPC_z07, controlMPC_grad_cost07, &info->pobj);
controlMPC_LA_DIAG_QUADFCN_2(params->H9, params->f9, controlMPC_z08, controlMPC_grad_cost08, &info->pobj);
controlMPC_LA_DIAG_QUADFCN_2(params->H10, params->f10, controlMPC_z09, controlMPC_grad_cost09, &info->pobj);
controlMPC_LA_DIAG_QUADFCN_2(params->H11, params->f11, controlMPC_z10, controlMPC_grad_cost10, &info->pobj);
controlMPC_LA_DIAG_QUADFCN_2(params->H12, params->f12, controlMPC_z11, controlMPC_grad_cost11, &info->pobj);
controlMPC_LA_DIAG_QUADFCN_2(params->H13, params->f13, controlMPC_z12, controlMPC_grad_cost12, &info->pobj);
controlMPC_LA_DIAG_QUADFCN_2(params->H14, params->f14, controlMPC_z13, controlMPC_grad_cost13, &info->pobj);
controlMPC_LA_DIAG_QUADFCN_2(params->H15, params->f15, controlMPC_z14, controlMPC_grad_cost14, &info->pobj);
controlMPC_LA_DIAG_QUADFCN_2(params->H16, params->f16, controlMPC_z15, controlMPC_grad_cost15, &info->pobj);
controlMPC_LA_DIAG_QUADFCN_2(params->H17, params->f17, controlMPC_z16, controlMPC_grad_cost16, &info->pobj);
controlMPC_LA_DIAG_QUADFCN_2(params->H18, params->f18, controlMPC_z17, controlMPC_grad_cost17, &info->pobj);
controlMPC_LA_DIAG_QUADFCN_2(params->H19, params->f19, controlMPC_z18, controlMPC_grad_cost18, &info->pobj);
controlMPC_LA_DIAG_QUADFCN_2(params->H20, params->f20, controlMPC_z19, controlMPC_grad_cost19, &info->pobj);
controlMPC_LA_DIAG_QUADFCN_2(params->H21, params->f21, controlMPC_z20, controlMPC_grad_cost20, &info->pobj);
controlMPC_LA_DIAG_QUADFCN_2(params->H22, params->f22, controlMPC_z21, controlMPC_grad_cost21, &info->pobj);
controlMPC_LA_DIAG_QUADFCN_2(params->H23, params->f23, controlMPC_z22, controlMPC_grad_cost22, &info->pobj);
controlMPC_LA_DIAG_QUADFCN_2(params->H24, params->f24, controlMPC_z23, controlMPC_grad_cost23, &info->pobj);
controlMPC_LA_DIAG_QUADFCN_2(params->H25, params->f25, controlMPC_z24, controlMPC_grad_cost24, &info->pobj);
controlMPC_LA_DIAG_QUADFCN_2(params->H26, params->f26, controlMPC_z25, controlMPC_grad_cost25, &info->pobj);
controlMPC_LA_DIAG_QUADFCN_2(params->H27, params->f27, controlMPC_z26, controlMPC_grad_cost26, &info->pobj);
controlMPC_LA_DIAG_QUADFCN_2(params->H28, params->f28, controlMPC_z27, controlMPC_grad_cost27, &info->pobj);
controlMPC_LA_DIAG_QUADFCN_2(params->H29, params->f29, controlMPC_z28, controlMPC_grad_cost28, &info->pobj);
info->res_eq = 0;
info->dgap = 0;
controlMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_2_2(controlMPC_C00, controlMPC_z00, controlMPC_D01, controlMPC_z01, controlMPC_c00, controlMPC_v00, controlMPC_re00, &info->dgap, &info->res_eq);
controlMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_2_2(controlMPC_C00, controlMPC_z01, controlMPC_D01, controlMPC_z02, controlMPC_c01, controlMPC_v01, controlMPC_re01, &info->dgap, &info->res_eq);
controlMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_2_2(controlMPC_C00, controlMPC_z02, controlMPC_D01, controlMPC_z03, controlMPC_c02, controlMPC_v02, controlMPC_re02, &info->dgap, &info->res_eq);
controlMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_2_2(controlMPC_C00, controlMPC_z03, controlMPC_D01, controlMPC_z04, controlMPC_c03, controlMPC_v03, controlMPC_re03, &info->dgap, &info->res_eq);
controlMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_2_2(controlMPC_C00, controlMPC_z04, controlMPC_D01, controlMPC_z05, controlMPC_c04, controlMPC_v04, controlMPC_re04, &info->dgap, &info->res_eq);
controlMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_2_2(controlMPC_C00, controlMPC_z05, controlMPC_D01, controlMPC_z06, controlMPC_c05, controlMPC_v05, controlMPC_re05, &info->dgap, &info->res_eq);
controlMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_2_2(controlMPC_C00, controlMPC_z06, controlMPC_D01, controlMPC_z07, controlMPC_c06, controlMPC_v06, controlMPC_re06, &info->dgap, &info->res_eq);
controlMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_2_2(controlMPC_C00, controlMPC_z07, controlMPC_D01, controlMPC_z08, controlMPC_c07, controlMPC_v07, controlMPC_re07, &info->dgap, &info->res_eq);
controlMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_2_2(controlMPC_C00, controlMPC_z08, controlMPC_D01, controlMPC_z09, controlMPC_c08, controlMPC_v08, controlMPC_re08, &info->dgap, &info->res_eq);
controlMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_2_2(controlMPC_C00, controlMPC_z09, controlMPC_D01, controlMPC_z10, controlMPC_c09, controlMPC_v09, controlMPC_re09, &info->dgap, &info->res_eq);
controlMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_2_2(controlMPC_C00, controlMPC_z10, controlMPC_D01, controlMPC_z11, controlMPC_c10, controlMPC_v10, controlMPC_re10, &info->dgap, &info->res_eq);
controlMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_2_2(controlMPC_C00, controlMPC_z11, controlMPC_D01, controlMPC_z12, controlMPC_c11, controlMPC_v11, controlMPC_re11, &info->dgap, &info->res_eq);
controlMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_2_2(controlMPC_C00, controlMPC_z12, controlMPC_D01, controlMPC_z13, controlMPC_c12, controlMPC_v12, controlMPC_re12, &info->dgap, &info->res_eq);
controlMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_2_2(controlMPC_C00, controlMPC_z13, controlMPC_D01, controlMPC_z14, controlMPC_c13, controlMPC_v13, controlMPC_re13, &info->dgap, &info->res_eq);
controlMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_2_2(controlMPC_C00, controlMPC_z14, controlMPC_D01, controlMPC_z15, controlMPC_c14, controlMPC_v14, controlMPC_re14, &info->dgap, &info->res_eq);
controlMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_2_2(controlMPC_C00, controlMPC_z15, controlMPC_D01, controlMPC_z16, controlMPC_c15, controlMPC_v15, controlMPC_re15, &info->dgap, &info->res_eq);
controlMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_2_2(controlMPC_C00, controlMPC_z16, controlMPC_D01, controlMPC_z17, controlMPC_c16, controlMPC_v16, controlMPC_re16, &info->dgap, &info->res_eq);
controlMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_2_2(controlMPC_C00, controlMPC_z17, controlMPC_D01, controlMPC_z18, controlMPC_c17, controlMPC_v17, controlMPC_re17, &info->dgap, &info->res_eq);
controlMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_2_2(controlMPC_C00, controlMPC_z18, controlMPC_D01, controlMPC_z19, controlMPC_c18, controlMPC_v18, controlMPC_re18, &info->dgap, &info->res_eq);
controlMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_2_2(controlMPC_C00, controlMPC_z19, controlMPC_D01, controlMPC_z20, controlMPC_c19, controlMPC_v19, controlMPC_re19, &info->dgap, &info->res_eq);
controlMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_2_2(controlMPC_C00, controlMPC_z20, controlMPC_D01, controlMPC_z21, controlMPC_c20, controlMPC_v20, controlMPC_re20, &info->dgap, &info->res_eq);
controlMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_2_2(controlMPC_C00, controlMPC_z21, controlMPC_D01, controlMPC_z22, controlMPC_c21, controlMPC_v21, controlMPC_re21, &info->dgap, &info->res_eq);
controlMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_2_2(controlMPC_C00, controlMPC_z22, controlMPC_D01, controlMPC_z23, controlMPC_c22, controlMPC_v22, controlMPC_re22, &info->dgap, &info->res_eq);
controlMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_2_2(controlMPC_C00, controlMPC_z23, controlMPC_D01, controlMPC_z24, controlMPC_c23, controlMPC_v23, controlMPC_re23, &info->dgap, &info->res_eq);
controlMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_2_2(controlMPC_C00, controlMPC_z24, controlMPC_D01, controlMPC_z25, controlMPC_c24, controlMPC_v24, controlMPC_re24, &info->dgap, &info->res_eq);
controlMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_2_2(controlMPC_C00, controlMPC_z25, controlMPC_D01, controlMPC_z26, controlMPC_c25, controlMPC_v25, controlMPC_re25, &info->dgap, &info->res_eq);
controlMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_2_2(controlMPC_C00, controlMPC_z26, controlMPC_D01, controlMPC_z27, controlMPC_c26, controlMPC_v26, controlMPC_re26, &info->dgap, &info->res_eq);
controlMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_2_2(controlMPC_C00, controlMPC_z27, controlMPC_D01, controlMPC_z28, controlMPC_c27, controlMPC_v27, controlMPC_re27, &info->dgap, &info->res_eq);
controlMPC_LA_DENSE_MTVM_2_2(controlMPC_C00, controlMPC_v00, controlMPC_grad_eq00);
controlMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlMPC_C00, controlMPC_v01, controlMPC_D01, controlMPC_v00, controlMPC_grad_eq01);
controlMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlMPC_C00, controlMPC_v02, controlMPC_D01, controlMPC_v01, controlMPC_grad_eq02);
controlMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlMPC_C00, controlMPC_v03, controlMPC_D01, controlMPC_v02, controlMPC_grad_eq03);
controlMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlMPC_C00, controlMPC_v04, controlMPC_D01, controlMPC_v03, controlMPC_grad_eq04);
controlMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlMPC_C00, controlMPC_v05, controlMPC_D01, controlMPC_v04, controlMPC_grad_eq05);
controlMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlMPC_C00, controlMPC_v06, controlMPC_D01, controlMPC_v05, controlMPC_grad_eq06);
controlMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlMPC_C00, controlMPC_v07, controlMPC_D01, controlMPC_v06, controlMPC_grad_eq07);
controlMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlMPC_C00, controlMPC_v08, controlMPC_D01, controlMPC_v07, controlMPC_grad_eq08);
controlMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlMPC_C00, controlMPC_v09, controlMPC_D01, controlMPC_v08, controlMPC_grad_eq09);
controlMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlMPC_C00, controlMPC_v10, controlMPC_D01, controlMPC_v09, controlMPC_grad_eq10);
controlMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlMPC_C00, controlMPC_v11, controlMPC_D01, controlMPC_v10, controlMPC_grad_eq11);
controlMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlMPC_C00, controlMPC_v12, controlMPC_D01, controlMPC_v11, controlMPC_grad_eq12);
controlMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlMPC_C00, controlMPC_v13, controlMPC_D01, controlMPC_v12, controlMPC_grad_eq13);
controlMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlMPC_C00, controlMPC_v14, controlMPC_D01, controlMPC_v13, controlMPC_grad_eq14);
controlMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlMPC_C00, controlMPC_v15, controlMPC_D01, controlMPC_v14, controlMPC_grad_eq15);
controlMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlMPC_C00, controlMPC_v16, controlMPC_D01, controlMPC_v15, controlMPC_grad_eq16);
controlMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlMPC_C00, controlMPC_v17, controlMPC_D01, controlMPC_v16, controlMPC_grad_eq17);
controlMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlMPC_C00, controlMPC_v18, controlMPC_D01, controlMPC_v17, controlMPC_grad_eq18);
controlMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlMPC_C00, controlMPC_v19, controlMPC_D01, controlMPC_v18, controlMPC_grad_eq19);
controlMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlMPC_C00, controlMPC_v20, controlMPC_D01, controlMPC_v19, controlMPC_grad_eq20);
controlMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlMPC_C00, controlMPC_v21, controlMPC_D01, controlMPC_v20, controlMPC_grad_eq21);
controlMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlMPC_C00, controlMPC_v22, controlMPC_D01, controlMPC_v21, controlMPC_grad_eq22);
controlMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlMPC_C00, controlMPC_v23, controlMPC_D01, controlMPC_v22, controlMPC_grad_eq23);
controlMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlMPC_C00, controlMPC_v24, controlMPC_D01, controlMPC_v23, controlMPC_grad_eq24);
controlMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlMPC_C00, controlMPC_v25, controlMPC_D01, controlMPC_v24, controlMPC_grad_eq25);
controlMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlMPC_C00, controlMPC_v26, controlMPC_D01, controlMPC_v25, controlMPC_grad_eq26);
controlMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlMPC_C00, controlMPC_v27, controlMPC_D01, controlMPC_v26, controlMPC_grad_eq27);
controlMPC_LA_DIAGZERO_MTVM_2_2(controlMPC_D01, controlMPC_v27, controlMPC_grad_eq28);
info->res_ineq = 0;
controlMPC_LA_VSUBADD3_2(params->lb1, controlMPC_z00, controlMPC_lbIdx00, controlMPC_llb00, controlMPC_slb00, controlMPC_rilb00, &info->dgap, &info->res_ineq);
controlMPC_LA_VSUBADD2_2(controlMPC_z00, controlMPC_ubIdx00, params->ub1, controlMPC_lub00, controlMPC_sub00, controlMPC_riub00, &info->dgap, &info->res_ineq);
controlMPC_LA_VSUBADD3_2(params->lb2, controlMPC_z01, controlMPC_lbIdx01, controlMPC_llb01, controlMPC_slb01, controlMPC_rilb01, &info->dgap, &info->res_ineq);
controlMPC_LA_VSUBADD2_2(controlMPC_z01, controlMPC_ubIdx01, params->ub2, controlMPC_lub01, controlMPC_sub01, controlMPC_riub01, &info->dgap, &info->res_ineq);
controlMPC_LA_VSUBADD3_2(params->lb3, controlMPC_z02, controlMPC_lbIdx02, controlMPC_llb02, controlMPC_slb02, controlMPC_rilb02, &info->dgap, &info->res_ineq);
controlMPC_LA_VSUBADD2_2(controlMPC_z02, controlMPC_ubIdx02, params->ub3, controlMPC_lub02, controlMPC_sub02, controlMPC_riub02, &info->dgap, &info->res_ineq);
controlMPC_LA_VSUBADD3_2(params->lb4, controlMPC_z03, controlMPC_lbIdx03, controlMPC_llb03, controlMPC_slb03, controlMPC_rilb03, &info->dgap, &info->res_ineq);
controlMPC_LA_VSUBADD2_2(controlMPC_z03, controlMPC_ubIdx03, params->ub4, controlMPC_lub03, controlMPC_sub03, controlMPC_riub03, &info->dgap, &info->res_ineq);
controlMPC_LA_VSUBADD3_2(params->lb5, controlMPC_z04, controlMPC_lbIdx04, controlMPC_llb04, controlMPC_slb04, controlMPC_rilb04, &info->dgap, &info->res_ineq);
controlMPC_LA_VSUBADD2_2(controlMPC_z04, controlMPC_ubIdx04, params->ub5, controlMPC_lub04, controlMPC_sub04, controlMPC_riub04, &info->dgap, &info->res_ineq);
controlMPC_LA_VSUBADD3_2(params->lb6, controlMPC_z05, controlMPC_lbIdx05, controlMPC_llb05, controlMPC_slb05, controlMPC_rilb05, &info->dgap, &info->res_ineq);
controlMPC_LA_VSUBADD2_2(controlMPC_z05, controlMPC_ubIdx05, params->ub6, controlMPC_lub05, controlMPC_sub05, controlMPC_riub05, &info->dgap, &info->res_ineq);
controlMPC_LA_VSUBADD3_2(params->lb7, controlMPC_z06, controlMPC_lbIdx06, controlMPC_llb06, controlMPC_slb06, controlMPC_rilb06, &info->dgap, &info->res_ineq);
controlMPC_LA_VSUBADD2_2(controlMPC_z06, controlMPC_ubIdx06, params->ub7, controlMPC_lub06, controlMPC_sub06, controlMPC_riub06, &info->dgap, &info->res_ineq);
controlMPC_LA_VSUBADD3_2(params->lb8, controlMPC_z07, controlMPC_lbIdx07, controlMPC_llb07, controlMPC_slb07, controlMPC_rilb07, &info->dgap, &info->res_ineq);
controlMPC_LA_VSUBADD2_2(controlMPC_z07, controlMPC_ubIdx07, params->ub8, controlMPC_lub07, controlMPC_sub07, controlMPC_riub07, &info->dgap, &info->res_ineq);
controlMPC_LA_VSUBADD3_2(params->lb9, controlMPC_z08, controlMPC_lbIdx08, controlMPC_llb08, controlMPC_slb08, controlMPC_rilb08, &info->dgap, &info->res_ineq);
controlMPC_LA_VSUBADD2_2(controlMPC_z08, controlMPC_ubIdx08, params->ub9, controlMPC_lub08, controlMPC_sub08, controlMPC_riub08, &info->dgap, &info->res_ineq);
controlMPC_LA_VSUBADD3_2(params->lb10, controlMPC_z09, controlMPC_lbIdx09, controlMPC_llb09, controlMPC_slb09, controlMPC_rilb09, &info->dgap, &info->res_ineq);
controlMPC_LA_VSUBADD2_2(controlMPC_z09, controlMPC_ubIdx09, params->ub10, controlMPC_lub09, controlMPC_sub09, controlMPC_riub09, &info->dgap, &info->res_ineq);
controlMPC_LA_VSUBADD3_2(params->lb11, controlMPC_z10, controlMPC_lbIdx10, controlMPC_llb10, controlMPC_slb10, controlMPC_rilb10, &info->dgap, &info->res_ineq);
controlMPC_LA_VSUBADD2_2(controlMPC_z10, controlMPC_ubIdx10, params->ub11, controlMPC_lub10, controlMPC_sub10, controlMPC_riub10, &info->dgap, &info->res_ineq);
controlMPC_LA_VSUBADD3_2(params->lb12, controlMPC_z11, controlMPC_lbIdx11, controlMPC_llb11, controlMPC_slb11, controlMPC_rilb11, &info->dgap, &info->res_ineq);
controlMPC_LA_VSUBADD2_2(controlMPC_z11, controlMPC_ubIdx11, params->ub12, controlMPC_lub11, controlMPC_sub11, controlMPC_riub11, &info->dgap, &info->res_ineq);
controlMPC_LA_VSUBADD3_2(params->lb13, controlMPC_z12, controlMPC_lbIdx12, controlMPC_llb12, controlMPC_slb12, controlMPC_rilb12, &info->dgap, &info->res_ineq);
controlMPC_LA_VSUBADD2_2(controlMPC_z12, controlMPC_ubIdx12, params->ub13, controlMPC_lub12, controlMPC_sub12, controlMPC_riub12, &info->dgap, &info->res_ineq);
controlMPC_LA_VSUBADD3_2(params->lb14, controlMPC_z13, controlMPC_lbIdx13, controlMPC_llb13, controlMPC_slb13, controlMPC_rilb13, &info->dgap, &info->res_ineq);
controlMPC_LA_VSUBADD2_2(controlMPC_z13, controlMPC_ubIdx13, params->ub14, controlMPC_lub13, controlMPC_sub13, controlMPC_riub13, &info->dgap, &info->res_ineq);
controlMPC_LA_VSUBADD3_2(params->lb15, controlMPC_z14, controlMPC_lbIdx14, controlMPC_llb14, controlMPC_slb14, controlMPC_rilb14, &info->dgap, &info->res_ineq);
controlMPC_LA_VSUBADD2_2(controlMPC_z14, controlMPC_ubIdx14, params->ub15, controlMPC_lub14, controlMPC_sub14, controlMPC_riub14, &info->dgap, &info->res_ineq);
controlMPC_LA_VSUBADD3_2(params->lb16, controlMPC_z15, controlMPC_lbIdx15, controlMPC_llb15, controlMPC_slb15, controlMPC_rilb15, &info->dgap, &info->res_ineq);
controlMPC_LA_VSUBADD2_2(controlMPC_z15, controlMPC_ubIdx15, params->ub16, controlMPC_lub15, controlMPC_sub15, controlMPC_riub15, &info->dgap, &info->res_ineq);
controlMPC_LA_VSUBADD3_2(params->lb17, controlMPC_z16, controlMPC_lbIdx16, controlMPC_llb16, controlMPC_slb16, controlMPC_rilb16, &info->dgap, &info->res_ineq);
controlMPC_LA_VSUBADD2_2(controlMPC_z16, controlMPC_ubIdx16, params->ub17, controlMPC_lub16, controlMPC_sub16, controlMPC_riub16, &info->dgap, &info->res_ineq);
controlMPC_LA_VSUBADD3_2(params->lb18, controlMPC_z17, controlMPC_lbIdx17, controlMPC_llb17, controlMPC_slb17, controlMPC_rilb17, &info->dgap, &info->res_ineq);
controlMPC_LA_VSUBADD2_2(controlMPC_z17, controlMPC_ubIdx17, params->ub18, controlMPC_lub17, controlMPC_sub17, controlMPC_riub17, &info->dgap, &info->res_ineq);
controlMPC_LA_VSUBADD3_2(params->lb19, controlMPC_z18, controlMPC_lbIdx18, controlMPC_llb18, controlMPC_slb18, controlMPC_rilb18, &info->dgap, &info->res_ineq);
controlMPC_LA_VSUBADD2_2(controlMPC_z18, controlMPC_ubIdx18, params->ub19, controlMPC_lub18, controlMPC_sub18, controlMPC_riub18, &info->dgap, &info->res_ineq);
controlMPC_LA_VSUBADD3_2(params->lb20, controlMPC_z19, controlMPC_lbIdx19, controlMPC_llb19, controlMPC_slb19, controlMPC_rilb19, &info->dgap, &info->res_ineq);
controlMPC_LA_VSUBADD2_2(controlMPC_z19, controlMPC_ubIdx19, params->ub20, controlMPC_lub19, controlMPC_sub19, controlMPC_riub19, &info->dgap, &info->res_ineq);
controlMPC_LA_VSUBADD3_2(params->lb21, controlMPC_z20, controlMPC_lbIdx20, controlMPC_llb20, controlMPC_slb20, controlMPC_rilb20, &info->dgap, &info->res_ineq);
controlMPC_LA_VSUBADD2_2(controlMPC_z20, controlMPC_ubIdx20, params->ub21, controlMPC_lub20, controlMPC_sub20, controlMPC_riub20, &info->dgap, &info->res_ineq);
controlMPC_LA_VSUBADD3_2(params->lb22, controlMPC_z21, controlMPC_lbIdx21, controlMPC_llb21, controlMPC_slb21, controlMPC_rilb21, &info->dgap, &info->res_ineq);
controlMPC_LA_VSUBADD2_2(controlMPC_z21, controlMPC_ubIdx21, params->ub22, controlMPC_lub21, controlMPC_sub21, controlMPC_riub21, &info->dgap, &info->res_ineq);
controlMPC_LA_VSUBADD3_2(params->lb23, controlMPC_z22, controlMPC_lbIdx22, controlMPC_llb22, controlMPC_slb22, controlMPC_rilb22, &info->dgap, &info->res_ineq);
controlMPC_LA_VSUBADD2_2(controlMPC_z22, controlMPC_ubIdx22, params->ub23, controlMPC_lub22, controlMPC_sub22, controlMPC_riub22, &info->dgap, &info->res_ineq);
controlMPC_LA_VSUBADD3_2(params->lb24, controlMPC_z23, controlMPC_lbIdx23, controlMPC_llb23, controlMPC_slb23, controlMPC_rilb23, &info->dgap, &info->res_ineq);
controlMPC_LA_VSUBADD2_2(controlMPC_z23, controlMPC_ubIdx23, params->ub24, controlMPC_lub23, controlMPC_sub23, controlMPC_riub23, &info->dgap, &info->res_ineq);
controlMPC_LA_VSUBADD3_2(params->lb25, controlMPC_z24, controlMPC_lbIdx24, controlMPC_llb24, controlMPC_slb24, controlMPC_rilb24, &info->dgap, &info->res_ineq);
controlMPC_LA_VSUBADD2_2(controlMPC_z24, controlMPC_ubIdx24, params->ub25, controlMPC_lub24, controlMPC_sub24, controlMPC_riub24, &info->dgap, &info->res_ineq);
controlMPC_LA_VSUBADD3_2(params->lb26, controlMPC_z25, controlMPC_lbIdx25, controlMPC_llb25, controlMPC_slb25, controlMPC_rilb25, &info->dgap, &info->res_ineq);
controlMPC_LA_VSUBADD2_2(controlMPC_z25, controlMPC_ubIdx25, params->ub26, controlMPC_lub25, controlMPC_sub25, controlMPC_riub25, &info->dgap, &info->res_ineq);
controlMPC_LA_VSUBADD3_2(params->lb27, controlMPC_z26, controlMPC_lbIdx26, controlMPC_llb26, controlMPC_slb26, controlMPC_rilb26, &info->dgap, &info->res_ineq);
controlMPC_LA_VSUBADD2_2(controlMPC_z26, controlMPC_ubIdx26, params->ub27, controlMPC_lub26, controlMPC_sub26, controlMPC_riub26, &info->dgap, &info->res_ineq);
controlMPC_LA_VSUBADD3_2(params->lb28, controlMPC_z27, controlMPC_lbIdx27, controlMPC_llb27, controlMPC_slb27, controlMPC_rilb27, &info->dgap, &info->res_ineq);
controlMPC_LA_VSUBADD2_2(controlMPC_z27, controlMPC_ubIdx27, params->ub28, controlMPC_lub27, controlMPC_sub27, controlMPC_riub27, &info->dgap, &info->res_ineq);
controlMPC_LA_VSUBADD3_2(params->lb29, controlMPC_z28, controlMPC_lbIdx28, controlMPC_llb28, controlMPC_slb28, controlMPC_rilb28, &info->dgap, &info->res_ineq);
controlMPC_LA_VSUBADD2_2(controlMPC_z28, controlMPC_ubIdx28, params->ub29, controlMPC_lub28, controlMPC_sub28, controlMPC_riub28, &info->dgap, &info->res_ineq);
controlMPC_LA_INEQ_B_GRAD_2_2_2(controlMPC_lub00, controlMPC_sub00, controlMPC_riub00, controlMPC_llb00, controlMPC_slb00, controlMPC_rilb00, controlMPC_lbIdx00, controlMPC_ubIdx00, controlMPC_grad_ineq00, controlMPC_lubbysub00, controlMPC_llbbyslb00);
controlMPC_LA_INEQ_B_GRAD_2_2_2(controlMPC_lub01, controlMPC_sub01, controlMPC_riub01, controlMPC_llb01, controlMPC_slb01, controlMPC_rilb01, controlMPC_lbIdx01, controlMPC_ubIdx01, controlMPC_grad_ineq01, controlMPC_lubbysub01, controlMPC_llbbyslb01);
controlMPC_LA_INEQ_B_GRAD_2_2_2(controlMPC_lub02, controlMPC_sub02, controlMPC_riub02, controlMPC_llb02, controlMPC_slb02, controlMPC_rilb02, controlMPC_lbIdx02, controlMPC_ubIdx02, controlMPC_grad_ineq02, controlMPC_lubbysub02, controlMPC_llbbyslb02);
controlMPC_LA_INEQ_B_GRAD_2_2_2(controlMPC_lub03, controlMPC_sub03, controlMPC_riub03, controlMPC_llb03, controlMPC_slb03, controlMPC_rilb03, controlMPC_lbIdx03, controlMPC_ubIdx03, controlMPC_grad_ineq03, controlMPC_lubbysub03, controlMPC_llbbyslb03);
controlMPC_LA_INEQ_B_GRAD_2_2_2(controlMPC_lub04, controlMPC_sub04, controlMPC_riub04, controlMPC_llb04, controlMPC_slb04, controlMPC_rilb04, controlMPC_lbIdx04, controlMPC_ubIdx04, controlMPC_grad_ineq04, controlMPC_lubbysub04, controlMPC_llbbyslb04);
controlMPC_LA_INEQ_B_GRAD_2_2_2(controlMPC_lub05, controlMPC_sub05, controlMPC_riub05, controlMPC_llb05, controlMPC_slb05, controlMPC_rilb05, controlMPC_lbIdx05, controlMPC_ubIdx05, controlMPC_grad_ineq05, controlMPC_lubbysub05, controlMPC_llbbyslb05);
controlMPC_LA_INEQ_B_GRAD_2_2_2(controlMPC_lub06, controlMPC_sub06, controlMPC_riub06, controlMPC_llb06, controlMPC_slb06, controlMPC_rilb06, controlMPC_lbIdx06, controlMPC_ubIdx06, controlMPC_grad_ineq06, controlMPC_lubbysub06, controlMPC_llbbyslb06);
controlMPC_LA_INEQ_B_GRAD_2_2_2(controlMPC_lub07, controlMPC_sub07, controlMPC_riub07, controlMPC_llb07, controlMPC_slb07, controlMPC_rilb07, controlMPC_lbIdx07, controlMPC_ubIdx07, controlMPC_grad_ineq07, controlMPC_lubbysub07, controlMPC_llbbyslb07);
controlMPC_LA_INEQ_B_GRAD_2_2_2(controlMPC_lub08, controlMPC_sub08, controlMPC_riub08, controlMPC_llb08, controlMPC_slb08, controlMPC_rilb08, controlMPC_lbIdx08, controlMPC_ubIdx08, controlMPC_grad_ineq08, controlMPC_lubbysub08, controlMPC_llbbyslb08);
controlMPC_LA_INEQ_B_GRAD_2_2_2(controlMPC_lub09, controlMPC_sub09, controlMPC_riub09, controlMPC_llb09, controlMPC_slb09, controlMPC_rilb09, controlMPC_lbIdx09, controlMPC_ubIdx09, controlMPC_grad_ineq09, controlMPC_lubbysub09, controlMPC_llbbyslb09);
controlMPC_LA_INEQ_B_GRAD_2_2_2(controlMPC_lub10, controlMPC_sub10, controlMPC_riub10, controlMPC_llb10, controlMPC_slb10, controlMPC_rilb10, controlMPC_lbIdx10, controlMPC_ubIdx10, controlMPC_grad_ineq10, controlMPC_lubbysub10, controlMPC_llbbyslb10);
controlMPC_LA_INEQ_B_GRAD_2_2_2(controlMPC_lub11, controlMPC_sub11, controlMPC_riub11, controlMPC_llb11, controlMPC_slb11, controlMPC_rilb11, controlMPC_lbIdx11, controlMPC_ubIdx11, controlMPC_grad_ineq11, controlMPC_lubbysub11, controlMPC_llbbyslb11);
controlMPC_LA_INEQ_B_GRAD_2_2_2(controlMPC_lub12, controlMPC_sub12, controlMPC_riub12, controlMPC_llb12, controlMPC_slb12, controlMPC_rilb12, controlMPC_lbIdx12, controlMPC_ubIdx12, controlMPC_grad_ineq12, controlMPC_lubbysub12, controlMPC_llbbyslb12);
controlMPC_LA_INEQ_B_GRAD_2_2_2(controlMPC_lub13, controlMPC_sub13, controlMPC_riub13, controlMPC_llb13, controlMPC_slb13, controlMPC_rilb13, controlMPC_lbIdx13, controlMPC_ubIdx13, controlMPC_grad_ineq13, controlMPC_lubbysub13, controlMPC_llbbyslb13);
controlMPC_LA_INEQ_B_GRAD_2_2_2(controlMPC_lub14, controlMPC_sub14, controlMPC_riub14, controlMPC_llb14, controlMPC_slb14, controlMPC_rilb14, controlMPC_lbIdx14, controlMPC_ubIdx14, controlMPC_grad_ineq14, controlMPC_lubbysub14, controlMPC_llbbyslb14);
controlMPC_LA_INEQ_B_GRAD_2_2_2(controlMPC_lub15, controlMPC_sub15, controlMPC_riub15, controlMPC_llb15, controlMPC_slb15, controlMPC_rilb15, controlMPC_lbIdx15, controlMPC_ubIdx15, controlMPC_grad_ineq15, controlMPC_lubbysub15, controlMPC_llbbyslb15);
controlMPC_LA_INEQ_B_GRAD_2_2_2(controlMPC_lub16, controlMPC_sub16, controlMPC_riub16, controlMPC_llb16, controlMPC_slb16, controlMPC_rilb16, controlMPC_lbIdx16, controlMPC_ubIdx16, controlMPC_grad_ineq16, controlMPC_lubbysub16, controlMPC_llbbyslb16);
controlMPC_LA_INEQ_B_GRAD_2_2_2(controlMPC_lub17, controlMPC_sub17, controlMPC_riub17, controlMPC_llb17, controlMPC_slb17, controlMPC_rilb17, controlMPC_lbIdx17, controlMPC_ubIdx17, controlMPC_grad_ineq17, controlMPC_lubbysub17, controlMPC_llbbyslb17);
controlMPC_LA_INEQ_B_GRAD_2_2_2(controlMPC_lub18, controlMPC_sub18, controlMPC_riub18, controlMPC_llb18, controlMPC_slb18, controlMPC_rilb18, controlMPC_lbIdx18, controlMPC_ubIdx18, controlMPC_grad_ineq18, controlMPC_lubbysub18, controlMPC_llbbyslb18);
controlMPC_LA_INEQ_B_GRAD_2_2_2(controlMPC_lub19, controlMPC_sub19, controlMPC_riub19, controlMPC_llb19, controlMPC_slb19, controlMPC_rilb19, controlMPC_lbIdx19, controlMPC_ubIdx19, controlMPC_grad_ineq19, controlMPC_lubbysub19, controlMPC_llbbyslb19);
controlMPC_LA_INEQ_B_GRAD_2_2_2(controlMPC_lub20, controlMPC_sub20, controlMPC_riub20, controlMPC_llb20, controlMPC_slb20, controlMPC_rilb20, controlMPC_lbIdx20, controlMPC_ubIdx20, controlMPC_grad_ineq20, controlMPC_lubbysub20, controlMPC_llbbyslb20);
controlMPC_LA_INEQ_B_GRAD_2_2_2(controlMPC_lub21, controlMPC_sub21, controlMPC_riub21, controlMPC_llb21, controlMPC_slb21, controlMPC_rilb21, controlMPC_lbIdx21, controlMPC_ubIdx21, controlMPC_grad_ineq21, controlMPC_lubbysub21, controlMPC_llbbyslb21);
controlMPC_LA_INEQ_B_GRAD_2_2_2(controlMPC_lub22, controlMPC_sub22, controlMPC_riub22, controlMPC_llb22, controlMPC_slb22, controlMPC_rilb22, controlMPC_lbIdx22, controlMPC_ubIdx22, controlMPC_grad_ineq22, controlMPC_lubbysub22, controlMPC_llbbyslb22);
controlMPC_LA_INEQ_B_GRAD_2_2_2(controlMPC_lub23, controlMPC_sub23, controlMPC_riub23, controlMPC_llb23, controlMPC_slb23, controlMPC_rilb23, controlMPC_lbIdx23, controlMPC_ubIdx23, controlMPC_grad_ineq23, controlMPC_lubbysub23, controlMPC_llbbyslb23);
controlMPC_LA_INEQ_B_GRAD_2_2_2(controlMPC_lub24, controlMPC_sub24, controlMPC_riub24, controlMPC_llb24, controlMPC_slb24, controlMPC_rilb24, controlMPC_lbIdx24, controlMPC_ubIdx24, controlMPC_grad_ineq24, controlMPC_lubbysub24, controlMPC_llbbyslb24);
controlMPC_LA_INEQ_B_GRAD_2_2_2(controlMPC_lub25, controlMPC_sub25, controlMPC_riub25, controlMPC_llb25, controlMPC_slb25, controlMPC_rilb25, controlMPC_lbIdx25, controlMPC_ubIdx25, controlMPC_grad_ineq25, controlMPC_lubbysub25, controlMPC_llbbyslb25);
controlMPC_LA_INEQ_B_GRAD_2_2_2(controlMPC_lub26, controlMPC_sub26, controlMPC_riub26, controlMPC_llb26, controlMPC_slb26, controlMPC_rilb26, controlMPC_lbIdx26, controlMPC_ubIdx26, controlMPC_grad_ineq26, controlMPC_lubbysub26, controlMPC_llbbyslb26);
controlMPC_LA_INEQ_B_GRAD_2_2_2(controlMPC_lub27, controlMPC_sub27, controlMPC_riub27, controlMPC_llb27, controlMPC_slb27, controlMPC_rilb27, controlMPC_lbIdx27, controlMPC_ubIdx27, controlMPC_grad_ineq27, controlMPC_lubbysub27, controlMPC_llbbyslb27);
controlMPC_LA_INEQ_B_GRAD_2_2_2(controlMPC_lub28, controlMPC_sub28, controlMPC_riub28, controlMPC_llb28, controlMPC_slb28, controlMPC_rilb28, controlMPC_lbIdx28, controlMPC_ubIdx28, controlMPC_grad_ineq28, controlMPC_lubbysub28, controlMPC_llbbyslb28);
info->dobj = info->pobj - info->dgap;
info->rdgap = info->pobj ? info->dgap / info->pobj : 1e6;
if( info->rdgap < 0 ) info->rdgap = -info->rdgap;
if( info->mu < controlMPC_SET_ACC_KKTCOMPL
    && (info->rdgap < controlMPC_SET_ACC_RDGAP || info->dgap < controlMPC_SET_ACC_KKTCOMPL)
    && info->res_eq < controlMPC_SET_ACC_RESEQ
    && info->res_ineq < controlMPC_SET_ACC_RESINEQ ){
exitcode = controlMPC_OPTIMAL; break; }
if( info->it == controlMPC_SET_MAXIT ){
exitcode = controlMPC_MAXITREACHED; break; }
controlMPC_LA_VVADD3_58(controlMPC_grad_cost, controlMPC_grad_eq, controlMPC_grad_ineq, controlMPC_rd);
controlMPC_LA_DIAG_CHOL_ONELOOP_LBUB_2_2_2(params->H1, controlMPC_llbbyslb00, controlMPC_lbIdx00, controlMPC_lubbysub00, controlMPC_ubIdx00, controlMPC_Phi00);
controlMPC_LA_DIAG_MATRIXFORWARDSUB_2_2(controlMPC_Phi00, controlMPC_C00, controlMPC_V00);
controlMPC_LA_DIAG_FORWARDSUB_2(controlMPC_Phi00, controlMPC_rd00, controlMPC_Lbyrd00);
controlMPC_LA_DIAG_CHOL_ONELOOP_LBUB_2_2_2(params->H2, controlMPC_llbbyslb01, controlMPC_lbIdx01, controlMPC_lubbysub01, controlMPC_ubIdx01, controlMPC_Phi01);
controlMPC_LA_DIAG_MATRIXFORWARDSUB_2_2(controlMPC_Phi01, controlMPC_C00, controlMPC_V01);
controlMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_2(controlMPC_Phi01, controlMPC_D01, controlMPC_W01);
controlMPC_LA_DENSE_DIAGZERO_MMTM_2_2_2(controlMPC_W01, controlMPC_V01, controlMPC_Ysd01);
controlMPC_LA_DIAG_FORWARDSUB_2(controlMPC_Phi01, controlMPC_rd01, controlMPC_Lbyrd01);
controlMPC_LA_DIAG_CHOL_ONELOOP_LBUB_2_2_2(params->H3, controlMPC_llbbyslb02, controlMPC_lbIdx02, controlMPC_lubbysub02, controlMPC_ubIdx02, controlMPC_Phi02);
controlMPC_LA_DIAG_MATRIXFORWARDSUB_2_2(controlMPC_Phi02, controlMPC_C00, controlMPC_V02);
controlMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_2(controlMPC_Phi02, controlMPC_D01, controlMPC_W02);
controlMPC_LA_DENSE_DIAGZERO_MMTM_2_2_2(controlMPC_W02, controlMPC_V02, controlMPC_Ysd02);
controlMPC_LA_DIAG_FORWARDSUB_2(controlMPC_Phi02, controlMPC_rd02, controlMPC_Lbyrd02);
controlMPC_LA_DIAG_CHOL_ONELOOP_LBUB_2_2_2(params->H4, controlMPC_llbbyslb03, controlMPC_lbIdx03, controlMPC_lubbysub03, controlMPC_ubIdx03, controlMPC_Phi03);
controlMPC_LA_DIAG_MATRIXFORWARDSUB_2_2(controlMPC_Phi03, controlMPC_C00, controlMPC_V03);
controlMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_2(controlMPC_Phi03, controlMPC_D01, controlMPC_W03);
controlMPC_LA_DENSE_DIAGZERO_MMTM_2_2_2(controlMPC_W03, controlMPC_V03, controlMPC_Ysd03);
controlMPC_LA_DIAG_FORWARDSUB_2(controlMPC_Phi03, controlMPC_rd03, controlMPC_Lbyrd03);
controlMPC_LA_DIAG_CHOL_ONELOOP_LBUB_2_2_2(params->H5, controlMPC_llbbyslb04, controlMPC_lbIdx04, controlMPC_lubbysub04, controlMPC_ubIdx04, controlMPC_Phi04);
controlMPC_LA_DIAG_MATRIXFORWARDSUB_2_2(controlMPC_Phi04, controlMPC_C00, controlMPC_V04);
controlMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_2(controlMPC_Phi04, controlMPC_D01, controlMPC_W04);
controlMPC_LA_DENSE_DIAGZERO_MMTM_2_2_2(controlMPC_W04, controlMPC_V04, controlMPC_Ysd04);
controlMPC_LA_DIAG_FORWARDSUB_2(controlMPC_Phi04, controlMPC_rd04, controlMPC_Lbyrd04);
controlMPC_LA_DIAG_CHOL_ONELOOP_LBUB_2_2_2(params->H6, controlMPC_llbbyslb05, controlMPC_lbIdx05, controlMPC_lubbysub05, controlMPC_ubIdx05, controlMPC_Phi05);
controlMPC_LA_DIAG_MATRIXFORWARDSUB_2_2(controlMPC_Phi05, controlMPC_C00, controlMPC_V05);
controlMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_2(controlMPC_Phi05, controlMPC_D01, controlMPC_W05);
controlMPC_LA_DENSE_DIAGZERO_MMTM_2_2_2(controlMPC_W05, controlMPC_V05, controlMPC_Ysd05);
controlMPC_LA_DIAG_FORWARDSUB_2(controlMPC_Phi05, controlMPC_rd05, controlMPC_Lbyrd05);
controlMPC_LA_DIAG_CHOL_ONELOOP_LBUB_2_2_2(params->H7, controlMPC_llbbyslb06, controlMPC_lbIdx06, controlMPC_lubbysub06, controlMPC_ubIdx06, controlMPC_Phi06);
controlMPC_LA_DIAG_MATRIXFORWARDSUB_2_2(controlMPC_Phi06, controlMPC_C00, controlMPC_V06);
controlMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_2(controlMPC_Phi06, controlMPC_D01, controlMPC_W06);
controlMPC_LA_DENSE_DIAGZERO_MMTM_2_2_2(controlMPC_W06, controlMPC_V06, controlMPC_Ysd06);
controlMPC_LA_DIAG_FORWARDSUB_2(controlMPC_Phi06, controlMPC_rd06, controlMPC_Lbyrd06);
controlMPC_LA_DIAG_CHOL_ONELOOP_LBUB_2_2_2(params->H8, controlMPC_llbbyslb07, controlMPC_lbIdx07, controlMPC_lubbysub07, controlMPC_ubIdx07, controlMPC_Phi07);
controlMPC_LA_DIAG_MATRIXFORWARDSUB_2_2(controlMPC_Phi07, controlMPC_C00, controlMPC_V07);
controlMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_2(controlMPC_Phi07, controlMPC_D01, controlMPC_W07);
controlMPC_LA_DENSE_DIAGZERO_MMTM_2_2_2(controlMPC_W07, controlMPC_V07, controlMPC_Ysd07);
controlMPC_LA_DIAG_FORWARDSUB_2(controlMPC_Phi07, controlMPC_rd07, controlMPC_Lbyrd07);
controlMPC_LA_DIAG_CHOL_ONELOOP_LBUB_2_2_2(params->H9, controlMPC_llbbyslb08, controlMPC_lbIdx08, controlMPC_lubbysub08, controlMPC_ubIdx08, controlMPC_Phi08);
controlMPC_LA_DIAG_MATRIXFORWARDSUB_2_2(controlMPC_Phi08, controlMPC_C00, controlMPC_V08);
controlMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_2(controlMPC_Phi08, controlMPC_D01, controlMPC_W08);
controlMPC_LA_DENSE_DIAGZERO_MMTM_2_2_2(controlMPC_W08, controlMPC_V08, controlMPC_Ysd08);
controlMPC_LA_DIAG_FORWARDSUB_2(controlMPC_Phi08, controlMPC_rd08, controlMPC_Lbyrd08);
controlMPC_LA_DIAG_CHOL_ONELOOP_LBUB_2_2_2(params->H10, controlMPC_llbbyslb09, controlMPC_lbIdx09, controlMPC_lubbysub09, controlMPC_ubIdx09, controlMPC_Phi09);
controlMPC_LA_DIAG_MATRIXFORWARDSUB_2_2(controlMPC_Phi09, controlMPC_C00, controlMPC_V09);
controlMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_2(controlMPC_Phi09, controlMPC_D01, controlMPC_W09);
controlMPC_LA_DENSE_DIAGZERO_MMTM_2_2_2(controlMPC_W09, controlMPC_V09, controlMPC_Ysd09);
controlMPC_LA_DIAG_FORWARDSUB_2(controlMPC_Phi09, controlMPC_rd09, controlMPC_Lbyrd09);
controlMPC_LA_DIAG_CHOL_ONELOOP_LBUB_2_2_2(params->H11, controlMPC_llbbyslb10, controlMPC_lbIdx10, controlMPC_lubbysub10, controlMPC_ubIdx10, controlMPC_Phi10);
controlMPC_LA_DIAG_MATRIXFORWARDSUB_2_2(controlMPC_Phi10, controlMPC_C00, controlMPC_V10);
controlMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_2(controlMPC_Phi10, controlMPC_D01, controlMPC_W10);
controlMPC_LA_DENSE_DIAGZERO_MMTM_2_2_2(controlMPC_W10, controlMPC_V10, controlMPC_Ysd10);
controlMPC_LA_DIAG_FORWARDSUB_2(controlMPC_Phi10, controlMPC_rd10, controlMPC_Lbyrd10);
controlMPC_LA_DIAG_CHOL_ONELOOP_LBUB_2_2_2(params->H12, controlMPC_llbbyslb11, controlMPC_lbIdx11, controlMPC_lubbysub11, controlMPC_ubIdx11, controlMPC_Phi11);
controlMPC_LA_DIAG_MATRIXFORWARDSUB_2_2(controlMPC_Phi11, controlMPC_C00, controlMPC_V11);
controlMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_2(controlMPC_Phi11, controlMPC_D01, controlMPC_W11);
controlMPC_LA_DENSE_DIAGZERO_MMTM_2_2_2(controlMPC_W11, controlMPC_V11, controlMPC_Ysd11);
controlMPC_LA_DIAG_FORWARDSUB_2(controlMPC_Phi11, controlMPC_rd11, controlMPC_Lbyrd11);
controlMPC_LA_DIAG_CHOL_ONELOOP_LBUB_2_2_2(params->H13, controlMPC_llbbyslb12, controlMPC_lbIdx12, controlMPC_lubbysub12, controlMPC_ubIdx12, controlMPC_Phi12);
controlMPC_LA_DIAG_MATRIXFORWARDSUB_2_2(controlMPC_Phi12, controlMPC_C00, controlMPC_V12);
controlMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_2(controlMPC_Phi12, controlMPC_D01, controlMPC_W12);
controlMPC_LA_DENSE_DIAGZERO_MMTM_2_2_2(controlMPC_W12, controlMPC_V12, controlMPC_Ysd12);
controlMPC_LA_DIAG_FORWARDSUB_2(controlMPC_Phi12, controlMPC_rd12, controlMPC_Lbyrd12);
controlMPC_LA_DIAG_CHOL_ONELOOP_LBUB_2_2_2(params->H14, controlMPC_llbbyslb13, controlMPC_lbIdx13, controlMPC_lubbysub13, controlMPC_ubIdx13, controlMPC_Phi13);
controlMPC_LA_DIAG_MATRIXFORWARDSUB_2_2(controlMPC_Phi13, controlMPC_C00, controlMPC_V13);
controlMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_2(controlMPC_Phi13, controlMPC_D01, controlMPC_W13);
controlMPC_LA_DENSE_DIAGZERO_MMTM_2_2_2(controlMPC_W13, controlMPC_V13, controlMPC_Ysd13);
controlMPC_LA_DIAG_FORWARDSUB_2(controlMPC_Phi13, controlMPC_rd13, controlMPC_Lbyrd13);
controlMPC_LA_DIAG_CHOL_ONELOOP_LBUB_2_2_2(params->H15, controlMPC_llbbyslb14, controlMPC_lbIdx14, controlMPC_lubbysub14, controlMPC_ubIdx14, controlMPC_Phi14);
controlMPC_LA_DIAG_MATRIXFORWARDSUB_2_2(controlMPC_Phi14, controlMPC_C00, controlMPC_V14);
controlMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_2(controlMPC_Phi14, controlMPC_D01, controlMPC_W14);
controlMPC_LA_DENSE_DIAGZERO_MMTM_2_2_2(controlMPC_W14, controlMPC_V14, controlMPC_Ysd14);
controlMPC_LA_DIAG_FORWARDSUB_2(controlMPC_Phi14, controlMPC_rd14, controlMPC_Lbyrd14);
controlMPC_LA_DIAG_CHOL_ONELOOP_LBUB_2_2_2(params->H16, controlMPC_llbbyslb15, controlMPC_lbIdx15, controlMPC_lubbysub15, controlMPC_ubIdx15, controlMPC_Phi15);
controlMPC_LA_DIAG_MATRIXFORWARDSUB_2_2(controlMPC_Phi15, controlMPC_C00, controlMPC_V15);
controlMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_2(controlMPC_Phi15, controlMPC_D01, controlMPC_W15);
controlMPC_LA_DENSE_DIAGZERO_MMTM_2_2_2(controlMPC_W15, controlMPC_V15, controlMPC_Ysd15);
controlMPC_LA_DIAG_FORWARDSUB_2(controlMPC_Phi15, controlMPC_rd15, controlMPC_Lbyrd15);
controlMPC_LA_DIAG_CHOL_ONELOOP_LBUB_2_2_2(params->H17, controlMPC_llbbyslb16, controlMPC_lbIdx16, controlMPC_lubbysub16, controlMPC_ubIdx16, controlMPC_Phi16);
controlMPC_LA_DIAG_MATRIXFORWARDSUB_2_2(controlMPC_Phi16, controlMPC_C00, controlMPC_V16);
controlMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_2(controlMPC_Phi16, controlMPC_D01, controlMPC_W16);
controlMPC_LA_DENSE_DIAGZERO_MMTM_2_2_2(controlMPC_W16, controlMPC_V16, controlMPC_Ysd16);
controlMPC_LA_DIAG_FORWARDSUB_2(controlMPC_Phi16, controlMPC_rd16, controlMPC_Lbyrd16);
controlMPC_LA_DIAG_CHOL_ONELOOP_LBUB_2_2_2(params->H18, controlMPC_llbbyslb17, controlMPC_lbIdx17, controlMPC_lubbysub17, controlMPC_ubIdx17, controlMPC_Phi17);
controlMPC_LA_DIAG_MATRIXFORWARDSUB_2_2(controlMPC_Phi17, controlMPC_C00, controlMPC_V17);
controlMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_2(controlMPC_Phi17, controlMPC_D01, controlMPC_W17);
controlMPC_LA_DENSE_DIAGZERO_MMTM_2_2_2(controlMPC_W17, controlMPC_V17, controlMPC_Ysd17);
controlMPC_LA_DIAG_FORWARDSUB_2(controlMPC_Phi17, controlMPC_rd17, controlMPC_Lbyrd17);
controlMPC_LA_DIAG_CHOL_ONELOOP_LBUB_2_2_2(params->H19, controlMPC_llbbyslb18, controlMPC_lbIdx18, controlMPC_lubbysub18, controlMPC_ubIdx18, controlMPC_Phi18);
controlMPC_LA_DIAG_MATRIXFORWARDSUB_2_2(controlMPC_Phi18, controlMPC_C00, controlMPC_V18);
controlMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_2(controlMPC_Phi18, controlMPC_D01, controlMPC_W18);
controlMPC_LA_DENSE_DIAGZERO_MMTM_2_2_2(controlMPC_W18, controlMPC_V18, controlMPC_Ysd18);
controlMPC_LA_DIAG_FORWARDSUB_2(controlMPC_Phi18, controlMPC_rd18, controlMPC_Lbyrd18);
controlMPC_LA_DIAG_CHOL_ONELOOP_LBUB_2_2_2(params->H20, controlMPC_llbbyslb19, controlMPC_lbIdx19, controlMPC_lubbysub19, controlMPC_ubIdx19, controlMPC_Phi19);
controlMPC_LA_DIAG_MATRIXFORWARDSUB_2_2(controlMPC_Phi19, controlMPC_C00, controlMPC_V19);
controlMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_2(controlMPC_Phi19, controlMPC_D01, controlMPC_W19);
controlMPC_LA_DENSE_DIAGZERO_MMTM_2_2_2(controlMPC_W19, controlMPC_V19, controlMPC_Ysd19);
controlMPC_LA_DIAG_FORWARDSUB_2(controlMPC_Phi19, controlMPC_rd19, controlMPC_Lbyrd19);
controlMPC_LA_DIAG_CHOL_ONELOOP_LBUB_2_2_2(params->H21, controlMPC_llbbyslb20, controlMPC_lbIdx20, controlMPC_lubbysub20, controlMPC_ubIdx20, controlMPC_Phi20);
controlMPC_LA_DIAG_MATRIXFORWARDSUB_2_2(controlMPC_Phi20, controlMPC_C00, controlMPC_V20);
controlMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_2(controlMPC_Phi20, controlMPC_D01, controlMPC_W20);
controlMPC_LA_DENSE_DIAGZERO_MMTM_2_2_2(controlMPC_W20, controlMPC_V20, controlMPC_Ysd20);
controlMPC_LA_DIAG_FORWARDSUB_2(controlMPC_Phi20, controlMPC_rd20, controlMPC_Lbyrd20);
controlMPC_LA_DIAG_CHOL_ONELOOP_LBUB_2_2_2(params->H22, controlMPC_llbbyslb21, controlMPC_lbIdx21, controlMPC_lubbysub21, controlMPC_ubIdx21, controlMPC_Phi21);
controlMPC_LA_DIAG_MATRIXFORWARDSUB_2_2(controlMPC_Phi21, controlMPC_C00, controlMPC_V21);
controlMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_2(controlMPC_Phi21, controlMPC_D01, controlMPC_W21);
controlMPC_LA_DENSE_DIAGZERO_MMTM_2_2_2(controlMPC_W21, controlMPC_V21, controlMPC_Ysd21);
controlMPC_LA_DIAG_FORWARDSUB_2(controlMPC_Phi21, controlMPC_rd21, controlMPC_Lbyrd21);
controlMPC_LA_DIAG_CHOL_ONELOOP_LBUB_2_2_2(params->H23, controlMPC_llbbyslb22, controlMPC_lbIdx22, controlMPC_lubbysub22, controlMPC_ubIdx22, controlMPC_Phi22);
controlMPC_LA_DIAG_MATRIXFORWARDSUB_2_2(controlMPC_Phi22, controlMPC_C00, controlMPC_V22);
controlMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_2(controlMPC_Phi22, controlMPC_D01, controlMPC_W22);
controlMPC_LA_DENSE_DIAGZERO_MMTM_2_2_2(controlMPC_W22, controlMPC_V22, controlMPC_Ysd22);
controlMPC_LA_DIAG_FORWARDSUB_2(controlMPC_Phi22, controlMPC_rd22, controlMPC_Lbyrd22);
controlMPC_LA_DIAG_CHOL_ONELOOP_LBUB_2_2_2(params->H24, controlMPC_llbbyslb23, controlMPC_lbIdx23, controlMPC_lubbysub23, controlMPC_ubIdx23, controlMPC_Phi23);
controlMPC_LA_DIAG_MATRIXFORWARDSUB_2_2(controlMPC_Phi23, controlMPC_C00, controlMPC_V23);
controlMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_2(controlMPC_Phi23, controlMPC_D01, controlMPC_W23);
controlMPC_LA_DENSE_DIAGZERO_MMTM_2_2_2(controlMPC_W23, controlMPC_V23, controlMPC_Ysd23);
controlMPC_LA_DIAG_FORWARDSUB_2(controlMPC_Phi23, controlMPC_rd23, controlMPC_Lbyrd23);
controlMPC_LA_DIAG_CHOL_ONELOOP_LBUB_2_2_2(params->H25, controlMPC_llbbyslb24, controlMPC_lbIdx24, controlMPC_lubbysub24, controlMPC_ubIdx24, controlMPC_Phi24);
controlMPC_LA_DIAG_MATRIXFORWARDSUB_2_2(controlMPC_Phi24, controlMPC_C00, controlMPC_V24);
controlMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_2(controlMPC_Phi24, controlMPC_D01, controlMPC_W24);
controlMPC_LA_DENSE_DIAGZERO_MMTM_2_2_2(controlMPC_W24, controlMPC_V24, controlMPC_Ysd24);
controlMPC_LA_DIAG_FORWARDSUB_2(controlMPC_Phi24, controlMPC_rd24, controlMPC_Lbyrd24);
controlMPC_LA_DIAG_CHOL_ONELOOP_LBUB_2_2_2(params->H26, controlMPC_llbbyslb25, controlMPC_lbIdx25, controlMPC_lubbysub25, controlMPC_ubIdx25, controlMPC_Phi25);
controlMPC_LA_DIAG_MATRIXFORWARDSUB_2_2(controlMPC_Phi25, controlMPC_C00, controlMPC_V25);
controlMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_2(controlMPC_Phi25, controlMPC_D01, controlMPC_W25);
controlMPC_LA_DENSE_DIAGZERO_MMTM_2_2_2(controlMPC_W25, controlMPC_V25, controlMPC_Ysd25);
controlMPC_LA_DIAG_FORWARDSUB_2(controlMPC_Phi25, controlMPC_rd25, controlMPC_Lbyrd25);
controlMPC_LA_DIAG_CHOL_ONELOOP_LBUB_2_2_2(params->H27, controlMPC_llbbyslb26, controlMPC_lbIdx26, controlMPC_lubbysub26, controlMPC_ubIdx26, controlMPC_Phi26);
controlMPC_LA_DIAG_MATRIXFORWARDSUB_2_2(controlMPC_Phi26, controlMPC_C00, controlMPC_V26);
controlMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_2(controlMPC_Phi26, controlMPC_D01, controlMPC_W26);
controlMPC_LA_DENSE_DIAGZERO_MMTM_2_2_2(controlMPC_W26, controlMPC_V26, controlMPC_Ysd26);
controlMPC_LA_DIAG_FORWARDSUB_2(controlMPC_Phi26, controlMPC_rd26, controlMPC_Lbyrd26);
controlMPC_LA_DIAG_CHOL_ONELOOP_LBUB_2_2_2(params->H28, controlMPC_llbbyslb27, controlMPC_lbIdx27, controlMPC_lubbysub27, controlMPC_ubIdx27, controlMPC_Phi27);
controlMPC_LA_DIAG_MATRIXFORWARDSUB_2_2(controlMPC_Phi27, controlMPC_C00, controlMPC_V27);
controlMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_2(controlMPC_Phi27, controlMPC_D01, controlMPC_W27);
controlMPC_LA_DENSE_DIAGZERO_MMTM_2_2_2(controlMPC_W27, controlMPC_V27, controlMPC_Ysd27);
controlMPC_LA_DIAG_FORWARDSUB_2(controlMPC_Phi27, controlMPC_rd27, controlMPC_Lbyrd27);
controlMPC_LA_DIAG_CHOL_ONELOOP_LBUB_2_2_2(params->H29, controlMPC_llbbyslb28, controlMPC_lbIdx28, controlMPC_lubbysub28, controlMPC_ubIdx28, controlMPC_Phi28);
controlMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_2(controlMPC_Phi28, controlMPC_D01, controlMPC_W28);
controlMPC_LA_DIAG_FORWARDSUB_2(controlMPC_Phi28, controlMPC_rd28, controlMPC_Lbyrd28);
controlMPC_LA_DENSE_DIAGZERO_MMT2_2_2_2(controlMPC_V00, controlMPC_W01, controlMPC_Yd00);
controlMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_2_2(controlMPC_V00, controlMPC_Lbyrd00, controlMPC_W01, controlMPC_Lbyrd01, controlMPC_re00, controlMPC_beta00);
controlMPC_LA_DENSE_DIAGZERO_MMT2_2_2_2(controlMPC_V01, controlMPC_W02, controlMPC_Yd01);
controlMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_2_2(controlMPC_V01, controlMPC_Lbyrd01, controlMPC_W02, controlMPC_Lbyrd02, controlMPC_re01, controlMPC_beta01);
controlMPC_LA_DENSE_DIAGZERO_MMT2_2_2_2(controlMPC_V02, controlMPC_W03, controlMPC_Yd02);
controlMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_2_2(controlMPC_V02, controlMPC_Lbyrd02, controlMPC_W03, controlMPC_Lbyrd03, controlMPC_re02, controlMPC_beta02);
controlMPC_LA_DENSE_DIAGZERO_MMT2_2_2_2(controlMPC_V03, controlMPC_W04, controlMPC_Yd03);
controlMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_2_2(controlMPC_V03, controlMPC_Lbyrd03, controlMPC_W04, controlMPC_Lbyrd04, controlMPC_re03, controlMPC_beta03);
controlMPC_LA_DENSE_DIAGZERO_MMT2_2_2_2(controlMPC_V04, controlMPC_W05, controlMPC_Yd04);
controlMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_2_2(controlMPC_V04, controlMPC_Lbyrd04, controlMPC_W05, controlMPC_Lbyrd05, controlMPC_re04, controlMPC_beta04);
controlMPC_LA_DENSE_DIAGZERO_MMT2_2_2_2(controlMPC_V05, controlMPC_W06, controlMPC_Yd05);
controlMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_2_2(controlMPC_V05, controlMPC_Lbyrd05, controlMPC_W06, controlMPC_Lbyrd06, controlMPC_re05, controlMPC_beta05);
controlMPC_LA_DENSE_DIAGZERO_MMT2_2_2_2(controlMPC_V06, controlMPC_W07, controlMPC_Yd06);
controlMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_2_2(controlMPC_V06, controlMPC_Lbyrd06, controlMPC_W07, controlMPC_Lbyrd07, controlMPC_re06, controlMPC_beta06);
controlMPC_LA_DENSE_DIAGZERO_MMT2_2_2_2(controlMPC_V07, controlMPC_W08, controlMPC_Yd07);
controlMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_2_2(controlMPC_V07, controlMPC_Lbyrd07, controlMPC_W08, controlMPC_Lbyrd08, controlMPC_re07, controlMPC_beta07);
controlMPC_LA_DENSE_DIAGZERO_MMT2_2_2_2(controlMPC_V08, controlMPC_W09, controlMPC_Yd08);
controlMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_2_2(controlMPC_V08, controlMPC_Lbyrd08, controlMPC_W09, controlMPC_Lbyrd09, controlMPC_re08, controlMPC_beta08);
controlMPC_LA_DENSE_DIAGZERO_MMT2_2_2_2(controlMPC_V09, controlMPC_W10, controlMPC_Yd09);
controlMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_2_2(controlMPC_V09, controlMPC_Lbyrd09, controlMPC_W10, controlMPC_Lbyrd10, controlMPC_re09, controlMPC_beta09);
controlMPC_LA_DENSE_DIAGZERO_MMT2_2_2_2(controlMPC_V10, controlMPC_W11, controlMPC_Yd10);
controlMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_2_2(controlMPC_V10, controlMPC_Lbyrd10, controlMPC_W11, controlMPC_Lbyrd11, controlMPC_re10, controlMPC_beta10);
controlMPC_LA_DENSE_DIAGZERO_MMT2_2_2_2(controlMPC_V11, controlMPC_W12, controlMPC_Yd11);
controlMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_2_2(controlMPC_V11, controlMPC_Lbyrd11, controlMPC_W12, controlMPC_Lbyrd12, controlMPC_re11, controlMPC_beta11);
controlMPC_LA_DENSE_DIAGZERO_MMT2_2_2_2(controlMPC_V12, controlMPC_W13, controlMPC_Yd12);
controlMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_2_2(controlMPC_V12, controlMPC_Lbyrd12, controlMPC_W13, controlMPC_Lbyrd13, controlMPC_re12, controlMPC_beta12);
controlMPC_LA_DENSE_DIAGZERO_MMT2_2_2_2(controlMPC_V13, controlMPC_W14, controlMPC_Yd13);
controlMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_2_2(controlMPC_V13, controlMPC_Lbyrd13, controlMPC_W14, controlMPC_Lbyrd14, controlMPC_re13, controlMPC_beta13);
controlMPC_LA_DENSE_DIAGZERO_MMT2_2_2_2(controlMPC_V14, controlMPC_W15, controlMPC_Yd14);
controlMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_2_2(controlMPC_V14, controlMPC_Lbyrd14, controlMPC_W15, controlMPC_Lbyrd15, controlMPC_re14, controlMPC_beta14);
controlMPC_LA_DENSE_DIAGZERO_MMT2_2_2_2(controlMPC_V15, controlMPC_W16, controlMPC_Yd15);
controlMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_2_2(controlMPC_V15, controlMPC_Lbyrd15, controlMPC_W16, controlMPC_Lbyrd16, controlMPC_re15, controlMPC_beta15);
controlMPC_LA_DENSE_DIAGZERO_MMT2_2_2_2(controlMPC_V16, controlMPC_W17, controlMPC_Yd16);
controlMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_2_2(controlMPC_V16, controlMPC_Lbyrd16, controlMPC_W17, controlMPC_Lbyrd17, controlMPC_re16, controlMPC_beta16);
controlMPC_LA_DENSE_DIAGZERO_MMT2_2_2_2(controlMPC_V17, controlMPC_W18, controlMPC_Yd17);
controlMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_2_2(controlMPC_V17, controlMPC_Lbyrd17, controlMPC_W18, controlMPC_Lbyrd18, controlMPC_re17, controlMPC_beta17);
controlMPC_LA_DENSE_DIAGZERO_MMT2_2_2_2(controlMPC_V18, controlMPC_W19, controlMPC_Yd18);
controlMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_2_2(controlMPC_V18, controlMPC_Lbyrd18, controlMPC_W19, controlMPC_Lbyrd19, controlMPC_re18, controlMPC_beta18);
controlMPC_LA_DENSE_DIAGZERO_MMT2_2_2_2(controlMPC_V19, controlMPC_W20, controlMPC_Yd19);
controlMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_2_2(controlMPC_V19, controlMPC_Lbyrd19, controlMPC_W20, controlMPC_Lbyrd20, controlMPC_re19, controlMPC_beta19);
controlMPC_LA_DENSE_DIAGZERO_MMT2_2_2_2(controlMPC_V20, controlMPC_W21, controlMPC_Yd20);
controlMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_2_2(controlMPC_V20, controlMPC_Lbyrd20, controlMPC_W21, controlMPC_Lbyrd21, controlMPC_re20, controlMPC_beta20);
controlMPC_LA_DENSE_DIAGZERO_MMT2_2_2_2(controlMPC_V21, controlMPC_W22, controlMPC_Yd21);
controlMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_2_2(controlMPC_V21, controlMPC_Lbyrd21, controlMPC_W22, controlMPC_Lbyrd22, controlMPC_re21, controlMPC_beta21);
controlMPC_LA_DENSE_DIAGZERO_MMT2_2_2_2(controlMPC_V22, controlMPC_W23, controlMPC_Yd22);
controlMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_2_2(controlMPC_V22, controlMPC_Lbyrd22, controlMPC_W23, controlMPC_Lbyrd23, controlMPC_re22, controlMPC_beta22);
controlMPC_LA_DENSE_DIAGZERO_MMT2_2_2_2(controlMPC_V23, controlMPC_W24, controlMPC_Yd23);
controlMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_2_2(controlMPC_V23, controlMPC_Lbyrd23, controlMPC_W24, controlMPC_Lbyrd24, controlMPC_re23, controlMPC_beta23);
controlMPC_LA_DENSE_DIAGZERO_MMT2_2_2_2(controlMPC_V24, controlMPC_W25, controlMPC_Yd24);
controlMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_2_2(controlMPC_V24, controlMPC_Lbyrd24, controlMPC_W25, controlMPC_Lbyrd25, controlMPC_re24, controlMPC_beta24);
controlMPC_LA_DENSE_DIAGZERO_MMT2_2_2_2(controlMPC_V25, controlMPC_W26, controlMPC_Yd25);
controlMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_2_2(controlMPC_V25, controlMPC_Lbyrd25, controlMPC_W26, controlMPC_Lbyrd26, controlMPC_re25, controlMPC_beta25);
controlMPC_LA_DENSE_DIAGZERO_MMT2_2_2_2(controlMPC_V26, controlMPC_W27, controlMPC_Yd26);
controlMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_2_2(controlMPC_V26, controlMPC_Lbyrd26, controlMPC_W27, controlMPC_Lbyrd27, controlMPC_re26, controlMPC_beta26);
controlMPC_LA_DENSE_DIAGZERO_MMT2_2_2_2(controlMPC_V27, controlMPC_W28, controlMPC_Yd27);
controlMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_2_2(controlMPC_V27, controlMPC_Lbyrd27, controlMPC_W28, controlMPC_Lbyrd28, controlMPC_re27, controlMPC_beta27);
controlMPC_LA_DENSE_CHOL_2(controlMPC_Yd00, controlMPC_Ld00);
controlMPC_LA_DENSE_FORWARDSUB_2(controlMPC_Ld00, controlMPC_beta00, controlMPC_yy00);
controlMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(controlMPC_Ld00, controlMPC_Ysd01, controlMPC_Lsd01);
controlMPC_LA_DENSE_MMTSUB_2_2(controlMPC_Lsd01, controlMPC_Yd01);
controlMPC_LA_DENSE_CHOL_2(controlMPC_Yd01, controlMPC_Ld01);
controlMPC_LA_DENSE_MVMSUB1_2_2(controlMPC_Lsd01, controlMPC_yy00, controlMPC_beta01, controlMPC_bmy01);
controlMPC_LA_DENSE_FORWARDSUB_2(controlMPC_Ld01, controlMPC_bmy01, controlMPC_yy01);
controlMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(controlMPC_Ld01, controlMPC_Ysd02, controlMPC_Lsd02);
controlMPC_LA_DENSE_MMTSUB_2_2(controlMPC_Lsd02, controlMPC_Yd02);
controlMPC_LA_DENSE_CHOL_2(controlMPC_Yd02, controlMPC_Ld02);
controlMPC_LA_DENSE_MVMSUB1_2_2(controlMPC_Lsd02, controlMPC_yy01, controlMPC_beta02, controlMPC_bmy02);
controlMPC_LA_DENSE_FORWARDSUB_2(controlMPC_Ld02, controlMPC_bmy02, controlMPC_yy02);
controlMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(controlMPC_Ld02, controlMPC_Ysd03, controlMPC_Lsd03);
controlMPC_LA_DENSE_MMTSUB_2_2(controlMPC_Lsd03, controlMPC_Yd03);
controlMPC_LA_DENSE_CHOL_2(controlMPC_Yd03, controlMPC_Ld03);
controlMPC_LA_DENSE_MVMSUB1_2_2(controlMPC_Lsd03, controlMPC_yy02, controlMPC_beta03, controlMPC_bmy03);
controlMPC_LA_DENSE_FORWARDSUB_2(controlMPC_Ld03, controlMPC_bmy03, controlMPC_yy03);
controlMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(controlMPC_Ld03, controlMPC_Ysd04, controlMPC_Lsd04);
controlMPC_LA_DENSE_MMTSUB_2_2(controlMPC_Lsd04, controlMPC_Yd04);
controlMPC_LA_DENSE_CHOL_2(controlMPC_Yd04, controlMPC_Ld04);
controlMPC_LA_DENSE_MVMSUB1_2_2(controlMPC_Lsd04, controlMPC_yy03, controlMPC_beta04, controlMPC_bmy04);
controlMPC_LA_DENSE_FORWARDSUB_2(controlMPC_Ld04, controlMPC_bmy04, controlMPC_yy04);
controlMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(controlMPC_Ld04, controlMPC_Ysd05, controlMPC_Lsd05);
controlMPC_LA_DENSE_MMTSUB_2_2(controlMPC_Lsd05, controlMPC_Yd05);
controlMPC_LA_DENSE_CHOL_2(controlMPC_Yd05, controlMPC_Ld05);
controlMPC_LA_DENSE_MVMSUB1_2_2(controlMPC_Lsd05, controlMPC_yy04, controlMPC_beta05, controlMPC_bmy05);
controlMPC_LA_DENSE_FORWARDSUB_2(controlMPC_Ld05, controlMPC_bmy05, controlMPC_yy05);
controlMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(controlMPC_Ld05, controlMPC_Ysd06, controlMPC_Lsd06);
controlMPC_LA_DENSE_MMTSUB_2_2(controlMPC_Lsd06, controlMPC_Yd06);
controlMPC_LA_DENSE_CHOL_2(controlMPC_Yd06, controlMPC_Ld06);
controlMPC_LA_DENSE_MVMSUB1_2_2(controlMPC_Lsd06, controlMPC_yy05, controlMPC_beta06, controlMPC_bmy06);
controlMPC_LA_DENSE_FORWARDSUB_2(controlMPC_Ld06, controlMPC_bmy06, controlMPC_yy06);
controlMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(controlMPC_Ld06, controlMPC_Ysd07, controlMPC_Lsd07);
controlMPC_LA_DENSE_MMTSUB_2_2(controlMPC_Lsd07, controlMPC_Yd07);
controlMPC_LA_DENSE_CHOL_2(controlMPC_Yd07, controlMPC_Ld07);
controlMPC_LA_DENSE_MVMSUB1_2_2(controlMPC_Lsd07, controlMPC_yy06, controlMPC_beta07, controlMPC_bmy07);
controlMPC_LA_DENSE_FORWARDSUB_2(controlMPC_Ld07, controlMPC_bmy07, controlMPC_yy07);
controlMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(controlMPC_Ld07, controlMPC_Ysd08, controlMPC_Lsd08);
controlMPC_LA_DENSE_MMTSUB_2_2(controlMPC_Lsd08, controlMPC_Yd08);
controlMPC_LA_DENSE_CHOL_2(controlMPC_Yd08, controlMPC_Ld08);
controlMPC_LA_DENSE_MVMSUB1_2_2(controlMPC_Lsd08, controlMPC_yy07, controlMPC_beta08, controlMPC_bmy08);
controlMPC_LA_DENSE_FORWARDSUB_2(controlMPC_Ld08, controlMPC_bmy08, controlMPC_yy08);
controlMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(controlMPC_Ld08, controlMPC_Ysd09, controlMPC_Lsd09);
controlMPC_LA_DENSE_MMTSUB_2_2(controlMPC_Lsd09, controlMPC_Yd09);
controlMPC_LA_DENSE_CHOL_2(controlMPC_Yd09, controlMPC_Ld09);
controlMPC_LA_DENSE_MVMSUB1_2_2(controlMPC_Lsd09, controlMPC_yy08, controlMPC_beta09, controlMPC_bmy09);
controlMPC_LA_DENSE_FORWARDSUB_2(controlMPC_Ld09, controlMPC_bmy09, controlMPC_yy09);
controlMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(controlMPC_Ld09, controlMPC_Ysd10, controlMPC_Lsd10);
controlMPC_LA_DENSE_MMTSUB_2_2(controlMPC_Lsd10, controlMPC_Yd10);
controlMPC_LA_DENSE_CHOL_2(controlMPC_Yd10, controlMPC_Ld10);
controlMPC_LA_DENSE_MVMSUB1_2_2(controlMPC_Lsd10, controlMPC_yy09, controlMPC_beta10, controlMPC_bmy10);
controlMPC_LA_DENSE_FORWARDSUB_2(controlMPC_Ld10, controlMPC_bmy10, controlMPC_yy10);
controlMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(controlMPC_Ld10, controlMPC_Ysd11, controlMPC_Lsd11);
controlMPC_LA_DENSE_MMTSUB_2_2(controlMPC_Lsd11, controlMPC_Yd11);
controlMPC_LA_DENSE_CHOL_2(controlMPC_Yd11, controlMPC_Ld11);
controlMPC_LA_DENSE_MVMSUB1_2_2(controlMPC_Lsd11, controlMPC_yy10, controlMPC_beta11, controlMPC_bmy11);
controlMPC_LA_DENSE_FORWARDSUB_2(controlMPC_Ld11, controlMPC_bmy11, controlMPC_yy11);
controlMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(controlMPC_Ld11, controlMPC_Ysd12, controlMPC_Lsd12);
controlMPC_LA_DENSE_MMTSUB_2_2(controlMPC_Lsd12, controlMPC_Yd12);
controlMPC_LA_DENSE_CHOL_2(controlMPC_Yd12, controlMPC_Ld12);
controlMPC_LA_DENSE_MVMSUB1_2_2(controlMPC_Lsd12, controlMPC_yy11, controlMPC_beta12, controlMPC_bmy12);
controlMPC_LA_DENSE_FORWARDSUB_2(controlMPC_Ld12, controlMPC_bmy12, controlMPC_yy12);
controlMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(controlMPC_Ld12, controlMPC_Ysd13, controlMPC_Lsd13);
controlMPC_LA_DENSE_MMTSUB_2_2(controlMPC_Lsd13, controlMPC_Yd13);
controlMPC_LA_DENSE_CHOL_2(controlMPC_Yd13, controlMPC_Ld13);
controlMPC_LA_DENSE_MVMSUB1_2_2(controlMPC_Lsd13, controlMPC_yy12, controlMPC_beta13, controlMPC_bmy13);
controlMPC_LA_DENSE_FORWARDSUB_2(controlMPC_Ld13, controlMPC_bmy13, controlMPC_yy13);
controlMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(controlMPC_Ld13, controlMPC_Ysd14, controlMPC_Lsd14);
controlMPC_LA_DENSE_MMTSUB_2_2(controlMPC_Lsd14, controlMPC_Yd14);
controlMPC_LA_DENSE_CHOL_2(controlMPC_Yd14, controlMPC_Ld14);
controlMPC_LA_DENSE_MVMSUB1_2_2(controlMPC_Lsd14, controlMPC_yy13, controlMPC_beta14, controlMPC_bmy14);
controlMPC_LA_DENSE_FORWARDSUB_2(controlMPC_Ld14, controlMPC_bmy14, controlMPC_yy14);
controlMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(controlMPC_Ld14, controlMPC_Ysd15, controlMPC_Lsd15);
controlMPC_LA_DENSE_MMTSUB_2_2(controlMPC_Lsd15, controlMPC_Yd15);
controlMPC_LA_DENSE_CHOL_2(controlMPC_Yd15, controlMPC_Ld15);
controlMPC_LA_DENSE_MVMSUB1_2_2(controlMPC_Lsd15, controlMPC_yy14, controlMPC_beta15, controlMPC_bmy15);
controlMPC_LA_DENSE_FORWARDSUB_2(controlMPC_Ld15, controlMPC_bmy15, controlMPC_yy15);
controlMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(controlMPC_Ld15, controlMPC_Ysd16, controlMPC_Lsd16);
controlMPC_LA_DENSE_MMTSUB_2_2(controlMPC_Lsd16, controlMPC_Yd16);
controlMPC_LA_DENSE_CHOL_2(controlMPC_Yd16, controlMPC_Ld16);
controlMPC_LA_DENSE_MVMSUB1_2_2(controlMPC_Lsd16, controlMPC_yy15, controlMPC_beta16, controlMPC_bmy16);
controlMPC_LA_DENSE_FORWARDSUB_2(controlMPC_Ld16, controlMPC_bmy16, controlMPC_yy16);
controlMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(controlMPC_Ld16, controlMPC_Ysd17, controlMPC_Lsd17);
controlMPC_LA_DENSE_MMTSUB_2_2(controlMPC_Lsd17, controlMPC_Yd17);
controlMPC_LA_DENSE_CHOL_2(controlMPC_Yd17, controlMPC_Ld17);
controlMPC_LA_DENSE_MVMSUB1_2_2(controlMPC_Lsd17, controlMPC_yy16, controlMPC_beta17, controlMPC_bmy17);
controlMPC_LA_DENSE_FORWARDSUB_2(controlMPC_Ld17, controlMPC_bmy17, controlMPC_yy17);
controlMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(controlMPC_Ld17, controlMPC_Ysd18, controlMPC_Lsd18);
controlMPC_LA_DENSE_MMTSUB_2_2(controlMPC_Lsd18, controlMPC_Yd18);
controlMPC_LA_DENSE_CHOL_2(controlMPC_Yd18, controlMPC_Ld18);
controlMPC_LA_DENSE_MVMSUB1_2_2(controlMPC_Lsd18, controlMPC_yy17, controlMPC_beta18, controlMPC_bmy18);
controlMPC_LA_DENSE_FORWARDSUB_2(controlMPC_Ld18, controlMPC_bmy18, controlMPC_yy18);
controlMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(controlMPC_Ld18, controlMPC_Ysd19, controlMPC_Lsd19);
controlMPC_LA_DENSE_MMTSUB_2_2(controlMPC_Lsd19, controlMPC_Yd19);
controlMPC_LA_DENSE_CHOL_2(controlMPC_Yd19, controlMPC_Ld19);
controlMPC_LA_DENSE_MVMSUB1_2_2(controlMPC_Lsd19, controlMPC_yy18, controlMPC_beta19, controlMPC_bmy19);
controlMPC_LA_DENSE_FORWARDSUB_2(controlMPC_Ld19, controlMPC_bmy19, controlMPC_yy19);
controlMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(controlMPC_Ld19, controlMPC_Ysd20, controlMPC_Lsd20);
controlMPC_LA_DENSE_MMTSUB_2_2(controlMPC_Lsd20, controlMPC_Yd20);
controlMPC_LA_DENSE_CHOL_2(controlMPC_Yd20, controlMPC_Ld20);
controlMPC_LA_DENSE_MVMSUB1_2_2(controlMPC_Lsd20, controlMPC_yy19, controlMPC_beta20, controlMPC_bmy20);
controlMPC_LA_DENSE_FORWARDSUB_2(controlMPC_Ld20, controlMPC_bmy20, controlMPC_yy20);
controlMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(controlMPC_Ld20, controlMPC_Ysd21, controlMPC_Lsd21);
controlMPC_LA_DENSE_MMTSUB_2_2(controlMPC_Lsd21, controlMPC_Yd21);
controlMPC_LA_DENSE_CHOL_2(controlMPC_Yd21, controlMPC_Ld21);
controlMPC_LA_DENSE_MVMSUB1_2_2(controlMPC_Lsd21, controlMPC_yy20, controlMPC_beta21, controlMPC_bmy21);
controlMPC_LA_DENSE_FORWARDSUB_2(controlMPC_Ld21, controlMPC_bmy21, controlMPC_yy21);
controlMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(controlMPC_Ld21, controlMPC_Ysd22, controlMPC_Lsd22);
controlMPC_LA_DENSE_MMTSUB_2_2(controlMPC_Lsd22, controlMPC_Yd22);
controlMPC_LA_DENSE_CHOL_2(controlMPC_Yd22, controlMPC_Ld22);
controlMPC_LA_DENSE_MVMSUB1_2_2(controlMPC_Lsd22, controlMPC_yy21, controlMPC_beta22, controlMPC_bmy22);
controlMPC_LA_DENSE_FORWARDSUB_2(controlMPC_Ld22, controlMPC_bmy22, controlMPC_yy22);
controlMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(controlMPC_Ld22, controlMPC_Ysd23, controlMPC_Lsd23);
controlMPC_LA_DENSE_MMTSUB_2_2(controlMPC_Lsd23, controlMPC_Yd23);
controlMPC_LA_DENSE_CHOL_2(controlMPC_Yd23, controlMPC_Ld23);
controlMPC_LA_DENSE_MVMSUB1_2_2(controlMPC_Lsd23, controlMPC_yy22, controlMPC_beta23, controlMPC_bmy23);
controlMPC_LA_DENSE_FORWARDSUB_2(controlMPC_Ld23, controlMPC_bmy23, controlMPC_yy23);
controlMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(controlMPC_Ld23, controlMPC_Ysd24, controlMPC_Lsd24);
controlMPC_LA_DENSE_MMTSUB_2_2(controlMPC_Lsd24, controlMPC_Yd24);
controlMPC_LA_DENSE_CHOL_2(controlMPC_Yd24, controlMPC_Ld24);
controlMPC_LA_DENSE_MVMSUB1_2_2(controlMPC_Lsd24, controlMPC_yy23, controlMPC_beta24, controlMPC_bmy24);
controlMPC_LA_DENSE_FORWARDSUB_2(controlMPC_Ld24, controlMPC_bmy24, controlMPC_yy24);
controlMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(controlMPC_Ld24, controlMPC_Ysd25, controlMPC_Lsd25);
controlMPC_LA_DENSE_MMTSUB_2_2(controlMPC_Lsd25, controlMPC_Yd25);
controlMPC_LA_DENSE_CHOL_2(controlMPC_Yd25, controlMPC_Ld25);
controlMPC_LA_DENSE_MVMSUB1_2_2(controlMPC_Lsd25, controlMPC_yy24, controlMPC_beta25, controlMPC_bmy25);
controlMPC_LA_DENSE_FORWARDSUB_2(controlMPC_Ld25, controlMPC_bmy25, controlMPC_yy25);
controlMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(controlMPC_Ld25, controlMPC_Ysd26, controlMPC_Lsd26);
controlMPC_LA_DENSE_MMTSUB_2_2(controlMPC_Lsd26, controlMPC_Yd26);
controlMPC_LA_DENSE_CHOL_2(controlMPC_Yd26, controlMPC_Ld26);
controlMPC_LA_DENSE_MVMSUB1_2_2(controlMPC_Lsd26, controlMPC_yy25, controlMPC_beta26, controlMPC_bmy26);
controlMPC_LA_DENSE_FORWARDSUB_2(controlMPC_Ld26, controlMPC_bmy26, controlMPC_yy26);
controlMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(controlMPC_Ld26, controlMPC_Ysd27, controlMPC_Lsd27);
controlMPC_LA_DENSE_MMTSUB_2_2(controlMPC_Lsd27, controlMPC_Yd27);
controlMPC_LA_DENSE_CHOL_2(controlMPC_Yd27, controlMPC_Ld27);
controlMPC_LA_DENSE_MVMSUB1_2_2(controlMPC_Lsd27, controlMPC_yy26, controlMPC_beta27, controlMPC_bmy27);
controlMPC_LA_DENSE_FORWARDSUB_2(controlMPC_Ld27, controlMPC_bmy27, controlMPC_yy27);
controlMPC_LA_DENSE_BACKWARDSUB_2(controlMPC_Ld27, controlMPC_yy27, controlMPC_dvaff27);
controlMPC_LA_DENSE_MTVMSUB_2_2(controlMPC_Lsd27, controlMPC_dvaff27, controlMPC_yy26, controlMPC_bmy26);
controlMPC_LA_DENSE_BACKWARDSUB_2(controlMPC_Ld26, controlMPC_bmy26, controlMPC_dvaff26);
controlMPC_LA_DENSE_MTVMSUB_2_2(controlMPC_Lsd26, controlMPC_dvaff26, controlMPC_yy25, controlMPC_bmy25);
controlMPC_LA_DENSE_BACKWARDSUB_2(controlMPC_Ld25, controlMPC_bmy25, controlMPC_dvaff25);
controlMPC_LA_DENSE_MTVMSUB_2_2(controlMPC_Lsd25, controlMPC_dvaff25, controlMPC_yy24, controlMPC_bmy24);
controlMPC_LA_DENSE_BACKWARDSUB_2(controlMPC_Ld24, controlMPC_bmy24, controlMPC_dvaff24);
controlMPC_LA_DENSE_MTVMSUB_2_2(controlMPC_Lsd24, controlMPC_dvaff24, controlMPC_yy23, controlMPC_bmy23);
controlMPC_LA_DENSE_BACKWARDSUB_2(controlMPC_Ld23, controlMPC_bmy23, controlMPC_dvaff23);
controlMPC_LA_DENSE_MTVMSUB_2_2(controlMPC_Lsd23, controlMPC_dvaff23, controlMPC_yy22, controlMPC_bmy22);
controlMPC_LA_DENSE_BACKWARDSUB_2(controlMPC_Ld22, controlMPC_bmy22, controlMPC_dvaff22);
controlMPC_LA_DENSE_MTVMSUB_2_2(controlMPC_Lsd22, controlMPC_dvaff22, controlMPC_yy21, controlMPC_bmy21);
controlMPC_LA_DENSE_BACKWARDSUB_2(controlMPC_Ld21, controlMPC_bmy21, controlMPC_dvaff21);
controlMPC_LA_DENSE_MTVMSUB_2_2(controlMPC_Lsd21, controlMPC_dvaff21, controlMPC_yy20, controlMPC_bmy20);
controlMPC_LA_DENSE_BACKWARDSUB_2(controlMPC_Ld20, controlMPC_bmy20, controlMPC_dvaff20);
controlMPC_LA_DENSE_MTVMSUB_2_2(controlMPC_Lsd20, controlMPC_dvaff20, controlMPC_yy19, controlMPC_bmy19);
controlMPC_LA_DENSE_BACKWARDSUB_2(controlMPC_Ld19, controlMPC_bmy19, controlMPC_dvaff19);
controlMPC_LA_DENSE_MTVMSUB_2_2(controlMPC_Lsd19, controlMPC_dvaff19, controlMPC_yy18, controlMPC_bmy18);
controlMPC_LA_DENSE_BACKWARDSUB_2(controlMPC_Ld18, controlMPC_bmy18, controlMPC_dvaff18);
controlMPC_LA_DENSE_MTVMSUB_2_2(controlMPC_Lsd18, controlMPC_dvaff18, controlMPC_yy17, controlMPC_bmy17);
controlMPC_LA_DENSE_BACKWARDSUB_2(controlMPC_Ld17, controlMPC_bmy17, controlMPC_dvaff17);
controlMPC_LA_DENSE_MTVMSUB_2_2(controlMPC_Lsd17, controlMPC_dvaff17, controlMPC_yy16, controlMPC_bmy16);
controlMPC_LA_DENSE_BACKWARDSUB_2(controlMPC_Ld16, controlMPC_bmy16, controlMPC_dvaff16);
controlMPC_LA_DENSE_MTVMSUB_2_2(controlMPC_Lsd16, controlMPC_dvaff16, controlMPC_yy15, controlMPC_bmy15);
controlMPC_LA_DENSE_BACKWARDSUB_2(controlMPC_Ld15, controlMPC_bmy15, controlMPC_dvaff15);
controlMPC_LA_DENSE_MTVMSUB_2_2(controlMPC_Lsd15, controlMPC_dvaff15, controlMPC_yy14, controlMPC_bmy14);
controlMPC_LA_DENSE_BACKWARDSUB_2(controlMPC_Ld14, controlMPC_bmy14, controlMPC_dvaff14);
controlMPC_LA_DENSE_MTVMSUB_2_2(controlMPC_Lsd14, controlMPC_dvaff14, controlMPC_yy13, controlMPC_bmy13);
controlMPC_LA_DENSE_BACKWARDSUB_2(controlMPC_Ld13, controlMPC_bmy13, controlMPC_dvaff13);
controlMPC_LA_DENSE_MTVMSUB_2_2(controlMPC_Lsd13, controlMPC_dvaff13, controlMPC_yy12, controlMPC_bmy12);
controlMPC_LA_DENSE_BACKWARDSUB_2(controlMPC_Ld12, controlMPC_bmy12, controlMPC_dvaff12);
controlMPC_LA_DENSE_MTVMSUB_2_2(controlMPC_Lsd12, controlMPC_dvaff12, controlMPC_yy11, controlMPC_bmy11);
controlMPC_LA_DENSE_BACKWARDSUB_2(controlMPC_Ld11, controlMPC_bmy11, controlMPC_dvaff11);
controlMPC_LA_DENSE_MTVMSUB_2_2(controlMPC_Lsd11, controlMPC_dvaff11, controlMPC_yy10, controlMPC_bmy10);
controlMPC_LA_DENSE_BACKWARDSUB_2(controlMPC_Ld10, controlMPC_bmy10, controlMPC_dvaff10);
controlMPC_LA_DENSE_MTVMSUB_2_2(controlMPC_Lsd10, controlMPC_dvaff10, controlMPC_yy09, controlMPC_bmy09);
controlMPC_LA_DENSE_BACKWARDSUB_2(controlMPC_Ld09, controlMPC_bmy09, controlMPC_dvaff09);
controlMPC_LA_DENSE_MTVMSUB_2_2(controlMPC_Lsd09, controlMPC_dvaff09, controlMPC_yy08, controlMPC_bmy08);
controlMPC_LA_DENSE_BACKWARDSUB_2(controlMPC_Ld08, controlMPC_bmy08, controlMPC_dvaff08);
controlMPC_LA_DENSE_MTVMSUB_2_2(controlMPC_Lsd08, controlMPC_dvaff08, controlMPC_yy07, controlMPC_bmy07);
controlMPC_LA_DENSE_BACKWARDSUB_2(controlMPC_Ld07, controlMPC_bmy07, controlMPC_dvaff07);
controlMPC_LA_DENSE_MTVMSUB_2_2(controlMPC_Lsd07, controlMPC_dvaff07, controlMPC_yy06, controlMPC_bmy06);
controlMPC_LA_DENSE_BACKWARDSUB_2(controlMPC_Ld06, controlMPC_bmy06, controlMPC_dvaff06);
controlMPC_LA_DENSE_MTVMSUB_2_2(controlMPC_Lsd06, controlMPC_dvaff06, controlMPC_yy05, controlMPC_bmy05);
controlMPC_LA_DENSE_BACKWARDSUB_2(controlMPC_Ld05, controlMPC_bmy05, controlMPC_dvaff05);
controlMPC_LA_DENSE_MTVMSUB_2_2(controlMPC_Lsd05, controlMPC_dvaff05, controlMPC_yy04, controlMPC_bmy04);
controlMPC_LA_DENSE_BACKWARDSUB_2(controlMPC_Ld04, controlMPC_bmy04, controlMPC_dvaff04);
controlMPC_LA_DENSE_MTVMSUB_2_2(controlMPC_Lsd04, controlMPC_dvaff04, controlMPC_yy03, controlMPC_bmy03);
controlMPC_LA_DENSE_BACKWARDSUB_2(controlMPC_Ld03, controlMPC_bmy03, controlMPC_dvaff03);
controlMPC_LA_DENSE_MTVMSUB_2_2(controlMPC_Lsd03, controlMPC_dvaff03, controlMPC_yy02, controlMPC_bmy02);
controlMPC_LA_DENSE_BACKWARDSUB_2(controlMPC_Ld02, controlMPC_bmy02, controlMPC_dvaff02);
controlMPC_LA_DENSE_MTVMSUB_2_2(controlMPC_Lsd02, controlMPC_dvaff02, controlMPC_yy01, controlMPC_bmy01);
controlMPC_LA_DENSE_BACKWARDSUB_2(controlMPC_Ld01, controlMPC_bmy01, controlMPC_dvaff01);
controlMPC_LA_DENSE_MTVMSUB_2_2(controlMPC_Lsd01, controlMPC_dvaff01, controlMPC_yy00, controlMPC_bmy00);
controlMPC_LA_DENSE_BACKWARDSUB_2(controlMPC_Ld00, controlMPC_bmy00, controlMPC_dvaff00);
controlMPC_LA_DENSE_MTVM_2_2(controlMPC_C00, controlMPC_dvaff00, controlMPC_grad_eq00);
controlMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlMPC_C00, controlMPC_dvaff01, controlMPC_D01, controlMPC_dvaff00, controlMPC_grad_eq01);
controlMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlMPC_C00, controlMPC_dvaff02, controlMPC_D01, controlMPC_dvaff01, controlMPC_grad_eq02);
controlMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlMPC_C00, controlMPC_dvaff03, controlMPC_D01, controlMPC_dvaff02, controlMPC_grad_eq03);
controlMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlMPC_C00, controlMPC_dvaff04, controlMPC_D01, controlMPC_dvaff03, controlMPC_grad_eq04);
controlMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlMPC_C00, controlMPC_dvaff05, controlMPC_D01, controlMPC_dvaff04, controlMPC_grad_eq05);
controlMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlMPC_C00, controlMPC_dvaff06, controlMPC_D01, controlMPC_dvaff05, controlMPC_grad_eq06);
controlMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlMPC_C00, controlMPC_dvaff07, controlMPC_D01, controlMPC_dvaff06, controlMPC_grad_eq07);
controlMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlMPC_C00, controlMPC_dvaff08, controlMPC_D01, controlMPC_dvaff07, controlMPC_grad_eq08);
controlMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlMPC_C00, controlMPC_dvaff09, controlMPC_D01, controlMPC_dvaff08, controlMPC_grad_eq09);
controlMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlMPC_C00, controlMPC_dvaff10, controlMPC_D01, controlMPC_dvaff09, controlMPC_grad_eq10);
controlMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlMPC_C00, controlMPC_dvaff11, controlMPC_D01, controlMPC_dvaff10, controlMPC_grad_eq11);
controlMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlMPC_C00, controlMPC_dvaff12, controlMPC_D01, controlMPC_dvaff11, controlMPC_grad_eq12);
controlMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlMPC_C00, controlMPC_dvaff13, controlMPC_D01, controlMPC_dvaff12, controlMPC_grad_eq13);
controlMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlMPC_C00, controlMPC_dvaff14, controlMPC_D01, controlMPC_dvaff13, controlMPC_grad_eq14);
controlMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlMPC_C00, controlMPC_dvaff15, controlMPC_D01, controlMPC_dvaff14, controlMPC_grad_eq15);
controlMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlMPC_C00, controlMPC_dvaff16, controlMPC_D01, controlMPC_dvaff15, controlMPC_grad_eq16);
controlMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlMPC_C00, controlMPC_dvaff17, controlMPC_D01, controlMPC_dvaff16, controlMPC_grad_eq17);
controlMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlMPC_C00, controlMPC_dvaff18, controlMPC_D01, controlMPC_dvaff17, controlMPC_grad_eq18);
controlMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlMPC_C00, controlMPC_dvaff19, controlMPC_D01, controlMPC_dvaff18, controlMPC_grad_eq19);
controlMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlMPC_C00, controlMPC_dvaff20, controlMPC_D01, controlMPC_dvaff19, controlMPC_grad_eq20);
controlMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlMPC_C00, controlMPC_dvaff21, controlMPC_D01, controlMPC_dvaff20, controlMPC_grad_eq21);
controlMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlMPC_C00, controlMPC_dvaff22, controlMPC_D01, controlMPC_dvaff21, controlMPC_grad_eq22);
controlMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlMPC_C00, controlMPC_dvaff23, controlMPC_D01, controlMPC_dvaff22, controlMPC_grad_eq23);
controlMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlMPC_C00, controlMPC_dvaff24, controlMPC_D01, controlMPC_dvaff23, controlMPC_grad_eq24);
controlMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlMPC_C00, controlMPC_dvaff25, controlMPC_D01, controlMPC_dvaff24, controlMPC_grad_eq25);
controlMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlMPC_C00, controlMPC_dvaff26, controlMPC_D01, controlMPC_dvaff25, controlMPC_grad_eq26);
controlMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlMPC_C00, controlMPC_dvaff27, controlMPC_D01, controlMPC_dvaff26, controlMPC_grad_eq27);
controlMPC_LA_DIAGZERO_MTVM_2_2(controlMPC_D01, controlMPC_dvaff27, controlMPC_grad_eq28);
controlMPC_LA_VSUB2_58(controlMPC_rd, controlMPC_grad_eq, controlMPC_rd);
controlMPC_LA_DIAG_FORWARDBACKWARDSUB_2(controlMPC_Phi00, controlMPC_rd00, controlMPC_dzaff00);
controlMPC_LA_DIAG_FORWARDBACKWARDSUB_2(controlMPC_Phi01, controlMPC_rd01, controlMPC_dzaff01);
controlMPC_LA_DIAG_FORWARDBACKWARDSUB_2(controlMPC_Phi02, controlMPC_rd02, controlMPC_dzaff02);
controlMPC_LA_DIAG_FORWARDBACKWARDSUB_2(controlMPC_Phi03, controlMPC_rd03, controlMPC_dzaff03);
controlMPC_LA_DIAG_FORWARDBACKWARDSUB_2(controlMPC_Phi04, controlMPC_rd04, controlMPC_dzaff04);
controlMPC_LA_DIAG_FORWARDBACKWARDSUB_2(controlMPC_Phi05, controlMPC_rd05, controlMPC_dzaff05);
controlMPC_LA_DIAG_FORWARDBACKWARDSUB_2(controlMPC_Phi06, controlMPC_rd06, controlMPC_dzaff06);
controlMPC_LA_DIAG_FORWARDBACKWARDSUB_2(controlMPC_Phi07, controlMPC_rd07, controlMPC_dzaff07);
controlMPC_LA_DIAG_FORWARDBACKWARDSUB_2(controlMPC_Phi08, controlMPC_rd08, controlMPC_dzaff08);
controlMPC_LA_DIAG_FORWARDBACKWARDSUB_2(controlMPC_Phi09, controlMPC_rd09, controlMPC_dzaff09);
controlMPC_LA_DIAG_FORWARDBACKWARDSUB_2(controlMPC_Phi10, controlMPC_rd10, controlMPC_dzaff10);
controlMPC_LA_DIAG_FORWARDBACKWARDSUB_2(controlMPC_Phi11, controlMPC_rd11, controlMPC_dzaff11);
controlMPC_LA_DIAG_FORWARDBACKWARDSUB_2(controlMPC_Phi12, controlMPC_rd12, controlMPC_dzaff12);
controlMPC_LA_DIAG_FORWARDBACKWARDSUB_2(controlMPC_Phi13, controlMPC_rd13, controlMPC_dzaff13);
controlMPC_LA_DIAG_FORWARDBACKWARDSUB_2(controlMPC_Phi14, controlMPC_rd14, controlMPC_dzaff14);
controlMPC_LA_DIAG_FORWARDBACKWARDSUB_2(controlMPC_Phi15, controlMPC_rd15, controlMPC_dzaff15);
controlMPC_LA_DIAG_FORWARDBACKWARDSUB_2(controlMPC_Phi16, controlMPC_rd16, controlMPC_dzaff16);
controlMPC_LA_DIAG_FORWARDBACKWARDSUB_2(controlMPC_Phi17, controlMPC_rd17, controlMPC_dzaff17);
controlMPC_LA_DIAG_FORWARDBACKWARDSUB_2(controlMPC_Phi18, controlMPC_rd18, controlMPC_dzaff18);
controlMPC_LA_DIAG_FORWARDBACKWARDSUB_2(controlMPC_Phi19, controlMPC_rd19, controlMPC_dzaff19);
controlMPC_LA_DIAG_FORWARDBACKWARDSUB_2(controlMPC_Phi20, controlMPC_rd20, controlMPC_dzaff20);
controlMPC_LA_DIAG_FORWARDBACKWARDSUB_2(controlMPC_Phi21, controlMPC_rd21, controlMPC_dzaff21);
controlMPC_LA_DIAG_FORWARDBACKWARDSUB_2(controlMPC_Phi22, controlMPC_rd22, controlMPC_dzaff22);
controlMPC_LA_DIAG_FORWARDBACKWARDSUB_2(controlMPC_Phi23, controlMPC_rd23, controlMPC_dzaff23);
controlMPC_LA_DIAG_FORWARDBACKWARDSUB_2(controlMPC_Phi24, controlMPC_rd24, controlMPC_dzaff24);
controlMPC_LA_DIAG_FORWARDBACKWARDSUB_2(controlMPC_Phi25, controlMPC_rd25, controlMPC_dzaff25);
controlMPC_LA_DIAG_FORWARDBACKWARDSUB_2(controlMPC_Phi26, controlMPC_rd26, controlMPC_dzaff26);
controlMPC_LA_DIAG_FORWARDBACKWARDSUB_2(controlMPC_Phi27, controlMPC_rd27, controlMPC_dzaff27);
controlMPC_LA_DIAG_FORWARDBACKWARDSUB_2(controlMPC_Phi28, controlMPC_rd28, controlMPC_dzaff28);
controlMPC_LA_VSUB_INDEXED_2(controlMPC_dzaff00, controlMPC_lbIdx00, controlMPC_rilb00, controlMPC_dslbaff00);
controlMPC_LA_VSUB3_2(controlMPC_llbbyslb00, controlMPC_dslbaff00, controlMPC_llb00, controlMPC_dllbaff00);
controlMPC_LA_VSUB2_INDEXED_2(controlMPC_riub00, controlMPC_dzaff00, controlMPC_ubIdx00, controlMPC_dsubaff00);
controlMPC_LA_VSUB3_2(controlMPC_lubbysub00, controlMPC_dsubaff00, controlMPC_lub00, controlMPC_dlubaff00);
controlMPC_LA_VSUB_INDEXED_2(controlMPC_dzaff01, controlMPC_lbIdx01, controlMPC_rilb01, controlMPC_dslbaff01);
controlMPC_LA_VSUB3_2(controlMPC_llbbyslb01, controlMPC_dslbaff01, controlMPC_llb01, controlMPC_dllbaff01);
controlMPC_LA_VSUB2_INDEXED_2(controlMPC_riub01, controlMPC_dzaff01, controlMPC_ubIdx01, controlMPC_dsubaff01);
controlMPC_LA_VSUB3_2(controlMPC_lubbysub01, controlMPC_dsubaff01, controlMPC_lub01, controlMPC_dlubaff01);
controlMPC_LA_VSUB_INDEXED_2(controlMPC_dzaff02, controlMPC_lbIdx02, controlMPC_rilb02, controlMPC_dslbaff02);
controlMPC_LA_VSUB3_2(controlMPC_llbbyslb02, controlMPC_dslbaff02, controlMPC_llb02, controlMPC_dllbaff02);
controlMPC_LA_VSUB2_INDEXED_2(controlMPC_riub02, controlMPC_dzaff02, controlMPC_ubIdx02, controlMPC_dsubaff02);
controlMPC_LA_VSUB3_2(controlMPC_lubbysub02, controlMPC_dsubaff02, controlMPC_lub02, controlMPC_dlubaff02);
controlMPC_LA_VSUB_INDEXED_2(controlMPC_dzaff03, controlMPC_lbIdx03, controlMPC_rilb03, controlMPC_dslbaff03);
controlMPC_LA_VSUB3_2(controlMPC_llbbyslb03, controlMPC_dslbaff03, controlMPC_llb03, controlMPC_dllbaff03);
controlMPC_LA_VSUB2_INDEXED_2(controlMPC_riub03, controlMPC_dzaff03, controlMPC_ubIdx03, controlMPC_dsubaff03);
controlMPC_LA_VSUB3_2(controlMPC_lubbysub03, controlMPC_dsubaff03, controlMPC_lub03, controlMPC_dlubaff03);
controlMPC_LA_VSUB_INDEXED_2(controlMPC_dzaff04, controlMPC_lbIdx04, controlMPC_rilb04, controlMPC_dslbaff04);
controlMPC_LA_VSUB3_2(controlMPC_llbbyslb04, controlMPC_dslbaff04, controlMPC_llb04, controlMPC_dllbaff04);
controlMPC_LA_VSUB2_INDEXED_2(controlMPC_riub04, controlMPC_dzaff04, controlMPC_ubIdx04, controlMPC_dsubaff04);
controlMPC_LA_VSUB3_2(controlMPC_lubbysub04, controlMPC_dsubaff04, controlMPC_lub04, controlMPC_dlubaff04);
controlMPC_LA_VSUB_INDEXED_2(controlMPC_dzaff05, controlMPC_lbIdx05, controlMPC_rilb05, controlMPC_dslbaff05);
controlMPC_LA_VSUB3_2(controlMPC_llbbyslb05, controlMPC_dslbaff05, controlMPC_llb05, controlMPC_dllbaff05);
controlMPC_LA_VSUB2_INDEXED_2(controlMPC_riub05, controlMPC_dzaff05, controlMPC_ubIdx05, controlMPC_dsubaff05);
controlMPC_LA_VSUB3_2(controlMPC_lubbysub05, controlMPC_dsubaff05, controlMPC_lub05, controlMPC_dlubaff05);
controlMPC_LA_VSUB_INDEXED_2(controlMPC_dzaff06, controlMPC_lbIdx06, controlMPC_rilb06, controlMPC_dslbaff06);
controlMPC_LA_VSUB3_2(controlMPC_llbbyslb06, controlMPC_dslbaff06, controlMPC_llb06, controlMPC_dllbaff06);
controlMPC_LA_VSUB2_INDEXED_2(controlMPC_riub06, controlMPC_dzaff06, controlMPC_ubIdx06, controlMPC_dsubaff06);
controlMPC_LA_VSUB3_2(controlMPC_lubbysub06, controlMPC_dsubaff06, controlMPC_lub06, controlMPC_dlubaff06);
controlMPC_LA_VSUB_INDEXED_2(controlMPC_dzaff07, controlMPC_lbIdx07, controlMPC_rilb07, controlMPC_dslbaff07);
controlMPC_LA_VSUB3_2(controlMPC_llbbyslb07, controlMPC_dslbaff07, controlMPC_llb07, controlMPC_dllbaff07);
controlMPC_LA_VSUB2_INDEXED_2(controlMPC_riub07, controlMPC_dzaff07, controlMPC_ubIdx07, controlMPC_dsubaff07);
controlMPC_LA_VSUB3_2(controlMPC_lubbysub07, controlMPC_dsubaff07, controlMPC_lub07, controlMPC_dlubaff07);
controlMPC_LA_VSUB_INDEXED_2(controlMPC_dzaff08, controlMPC_lbIdx08, controlMPC_rilb08, controlMPC_dslbaff08);
controlMPC_LA_VSUB3_2(controlMPC_llbbyslb08, controlMPC_dslbaff08, controlMPC_llb08, controlMPC_dllbaff08);
controlMPC_LA_VSUB2_INDEXED_2(controlMPC_riub08, controlMPC_dzaff08, controlMPC_ubIdx08, controlMPC_dsubaff08);
controlMPC_LA_VSUB3_2(controlMPC_lubbysub08, controlMPC_dsubaff08, controlMPC_lub08, controlMPC_dlubaff08);
controlMPC_LA_VSUB_INDEXED_2(controlMPC_dzaff09, controlMPC_lbIdx09, controlMPC_rilb09, controlMPC_dslbaff09);
controlMPC_LA_VSUB3_2(controlMPC_llbbyslb09, controlMPC_dslbaff09, controlMPC_llb09, controlMPC_dllbaff09);
controlMPC_LA_VSUB2_INDEXED_2(controlMPC_riub09, controlMPC_dzaff09, controlMPC_ubIdx09, controlMPC_dsubaff09);
controlMPC_LA_VSUB3_2(controlMPC_lubbysub09, controlMPC_dsubaff09, controlMPC_lub09, controlMPC_dlubaff09);
controlMPC_LA_VSUB_INDEXED_2(controlMPC_dzaff10, controlMPC_lbIdx10, controlMPC_rilb10, controlMPC_dslbaff10);
controlMPC_LA_VSUB3_2(controlMPC_llbbyslb10, controlMPC_dslbaff10, controlMPC_llb10, controlMPC_dllbaff10);
controlMPC_LA_VSUB2_INDEXED_2(controlMPC_riub10, controlMPC_dzaff10, controlMPC_ubIdx10, controlMPC_dsubaff10);
controlMPC_LA_VSUB3_2(controlMPC_lubbysub10, controlMPC_dsubaff10, controlMPC_lub10, controlMPC_dlubaff10);
controlMPC_LA_VSUB_INDEXED_2(controlMPC_dzaff11, controlMPC_lbIdx11, controlMPC_rilb11, controlMPC_dslbaff11);
controlMPC_LA_VSUB3_2(controlMPC_llbbyslb11, controlMPC_dslbaff11, controlMPC_llb11, controlMPC_dllbaff11);
controlMPC_LA_VSUB2_INDEXED_2(controlMPC_riub11, controlMPC_dzaff11, controlMPC_ubIdx11, controlMPC_dsubaff11);
controlMPC_LA_VSUB3_2(controlMPC_lubbysub11, controlMPC_dsubaff11, controlMPC_lub11, controlMPC_dlubaff11);
controlMPC_LA_VSUB_INDEXED_2(controlMPC_dzaff12, controlMPC_lbIdx12, controlMPC_rilb12, controlMPC_dslbaff12);
controlMPC_LA_VSUB3_2(controlMPC_llbbyslb12, controlMPC_dslbaff12, controlMPC_llb12, controlMPC_dllbaff12);
controlMPC_LA_VSUB2_INDEXED_2(controlMPC_riub12, controlMPC_dzaff12, controlMPC_ubIdx12, controlMPC_dsubaff12);
controlMPC_LA_VSUB3_2(controlMPC_lubbysub12, controlMPC_dsubaff12, controlMPC_lub12, controlMPC_dlubaff12);
controlMPC_LA_VSUB_INDEXED_2(controlMPC_dzaff13, controlMPC_lbIdx13, controlMPC_rilb13, controlMPC_dslbaff13);
controlMPC_LA_VSUB3_2(controlMPC_llbbyslb13, controlMPC_dslbaff13, controlMPC_llb13, controlMPC_dllbaff13);
controlMPC_LA_VSUB2_INDEXED_2(controlMPC_riub13, controlMPC_dzaff13, controlMPC_ubIdx13, controlMPC_dsubaff13);
controlMPC_LA_VSUB3_2(controlMPC_lubbysub13, controlMPC_dsubaff13, controlMPC_lub13, controlMPC_dlubaff13);
controlMPC_LA_VSUB_INDEXED_2(controlMPC_dzaff14, controlMPC_lbIdx14, controlMPC_rilb14, controlMPC_dslbaff14);
controlMPC_LA_VSUB3_2(controlMPC_llbbyslb14, controlMPC_dslbaff14, controlMPC_llb14, controlMPC_dllbaff14);
controlMPC_LA_VSUB2_INDEXED_2(controlMPC_riub14, controlMPC_dzaff14, controlMPC_ubIdx14, controlMPC_dsubaff14);
controlMPC_LA_VSUB3_2(controlMPC_lubbysub14, controlMPC_dsubaff14, controlMPC_lub14, controlMPC_dlubaff14);
controlMPC_LA_VSUB_INDEXED_2(controlMPC_dzaff15, controlMPC_lbIdx15, controlMPC_rilb15, controlMPC_dslbaff15);
controlMPC_LA_VSUB3_2(controlMPC_llbbyslb15, controlMPC_dslbaff15, controlMPC_llb15, controlMPC_dllbaff15);
controlMPC_LA_VSUB2_INDEXED_2(controlMPC_riub15, controlMPC_dzaff15, controlMPC_ubIdx15, controlMPC_dsubaff15);
controlMPC_LA_VSUB3_2(controlMPC_lubbysub15, controlMPC_dsubaff15, controlMPC_lub15, controlMPC_dlubaff15);
controlMPC_LA_VSUB_INDEXED_2(controlMPC_dzaff16, controlMPC_lbIdx16, controlMPC_rilb16, controlMPC_dslbaff16);
controlMPC_LA_VSUB3_2(controlMPC_llbbyslb16, controlMPC_dslbaff16, controlMPC_llb16, controlMPC_dllbaff16);
controlMPC_LA_VSUB2_INDEXED_2(controlMPC_riub16, controlMPC_dzaff16, controlMPC_ubIdx16, controlMPC_dsubaff16);
controlMPC_LA_VSUB3_2(controlMPC_lubbysub16, controlMPC_dsubaff16, controlMPC_lub16, controlMPC_dlubaff16);
controlMPC_LA_VSUB_INDEXED_2(controlMPC_dzaff17, controlMPC_lbIdx17, controlMPC_rilb17, controlMPC_dslbaff17);
controlMPC_LA_VSUB3_2(controlMPC_llbbyslb17, controlMPC_dslbaff17, controlMPC_llb17, controlMPC_dllbaff17);
controlMPC_LA_VSUB2_INDEXED_2(controlMPC_riub17, controlMPC_dzaff17, controlMPC_ubIdx17, controlMPC_dsubaff17);
controlMPC_LA_VSUB3_2(controlMPC_lubbysub17, controlMPC_dsubaff17, controlMPC_lub17, controlMPC_dlubaff17);
controlMPC_LA_VSUB_INDEXED_2(controlMPC_dzaff18, controlMPC_lbIdx18, controlMPC_rilb18, controlMPC_dslbaff18);
controlMPC_LA_VSUB3_2(controlMPC_llbbyslb18, controlMPC_dslbaff18, controlMPC_llb18, controlMPC_dllbaff18);
controlMPC_LA_VSUB2_INDEXED_2(controlMPC_riub18, controlMPC_dzaff18, controlMPC_ubIdx18, controlMPC_dsubaff18);
controlMPC_LA_VSUB3_2(controlMPC_lubbysub18, controlMPC_dsubaff18, controlMPC_lub18, controlMPC_dlubaff18);
controlMPC_LA_VSUB_INDEXED_2(controlMPC_dzaff19, controlMPC_lbIdx19, controlMPC_rilb19, controlMPC_dslbaff19);
controlMPC_LA_VSUB3_2(controlMPC_llbbyslb19, controlMPC_dslbaff19, controlMPC_llb19, controlMPC_dllbaff19);
controlMPC_LA_VSUB2_INDEXED_2(controlMPC_riub19, controlMPC_dzaff19, controlMPC_ubIdx19, controlMPC_dsubaff19);
controlMPC_LA_VSUB3_2(controlMPC_lubbysub19, controlMPC_dsubaff19, controlMPC_lub19, controlMPC_dlubaff19);
controlMPC_LA_VSUB_INDEXED_2(controlMPC_dzaff20, controlMPC_lbIdx20, controlMPC_rilb20, controlMPC_dslbaff20);
controlMPC_LA_VSUB3_2(controlMPC_llbbyslb20, controlMPC_dslbaff20, controlMPC_llb20, controlMPC_dllbaff20);
controlMPC_LA_VSUB2_INDEXED_2(controlMPC_riub20, controlMPC_dzaff20, controlMPC_ubIdx20, controlMPC_dsubaff20);
controlMPC_LA_VSUB3_2(controlMPC_lubbysub20, controlMPC_dsubaff20, controlMPC_lub20, controlMPC_dlubaff20);
controlMPC_LA_VSUB_INDEXED_2(controlMPC_dzaff21, controlMPC_lbIdx21, controlMPC_rilb21, controlMPC_dslbaff21);
controlMPC_LA_VSUB3_2(controlMPC_llbbyslb21, controlMPC_dslbaff21, controlMPC_llb21, controlMPC_dllbaff21);
controlMPC_LA_VSUB2_INDEXED_2(controlMPC_riub21, controlMPC_dzaff21, controlMPC_ubIdx21, controlMPC_dsubaff21);
controlMPC_LA_VSUB3_2(controlMPC_lubbysub21, controlMPC_dsubaff21, controlMPC_lub21, controlMPC_dlubaff21);
controlMPC_LA_VSUB_INDEXED_2(controlMPC_dzaff22, controlMPC_lbIdx22, controlMPC_rilb22, controlMPC_dslbaff22);
controlMPC_LA_VSUB3_2(controlMPC_llbbyslb22, controlMPC_dslbaff22, controlMPC_llb22, controlMPC_dllbaff22);
controlMPC_LA_VSUB2_INDEXED_2(controlMPC_riub22, controlMPC_dzaff22, controlMPC_ubIdx22, controlMPC_dsubaff22);
controlMPC_LA_VSUB3_2(controlMPC_lubbysub22, controlMPC_dsubaff22, controlMPC_lub22, controlMPC_dlubaff22);
controlMPC_LA_VSUB_INDEXED_2(controlMPC_dzaff23, controlMPC_lbIdx23, controlMPC_rilb23, controlMPC_dslbaff23);
controlMPC_LA_VSUB3_2(controlMPC_llbbyslb23, controlMPC_dslbaff23, controlMPC_llb23, controlMPC_dllbaff23);
controlMPC_LA_VSUB2_INDEXED_2(controlMPC_riub23, controlMPC_dzaff23, controlMPC_ubIdx23, controlMPC_dsubaff23);
controlMPC_LA_VSUB3_2(controlMPC_lubbysub23, controlMPC_dsubaff23, controlMPC_lub23, controlMPC_dlubaff23);
controlMPC_LA_VSUB_INDEXED_2(controlMPC_dzaff24, controlMPC_lbIdx24, controlMPC_rilb24, controlMPC_dslbaff24);
controlMPC_LA_VSUB3_2(controlMPC_llbbyslb24, controlMPC_dslbaff24, controlMPC_llb24, controlMPC_dllbaff24);
controlMPC_LA_VSUB2_INDEXED_2(controlMPC_riub24, controlMPC_dzaff24, controlMPC_ubIdx24, controlMPC_dsubaff24);
controlMPC_LA_VSUB3_2(controlMPC_lubbysub24, controlMPC_dsubaff24, controlMPC_lub24, controlMPC_dlubaff24);
controlMPC_LA_VSUB_INDEXED_2(controlMPC_dzaff25, controlMPC_lbIdx25, controlMPC_rilb25, controlMPC_dslbaff25);
controlMPC_LA_VSUB3_2(controlMPC_llbbyslb25, controlMPC_dslbaff25, controlMPC_llb25, controlMPC_dllbaff25);
controlMPC_LA_VSUB2_INDEXED_2(controlMPC_riub25, controlMPC_dzaff25, controlMPC_ubIdx25, controlMPC_dsubaff25);
controlMPC_LA_VSUB3_2(controlMPC_lubbysub25, controlMPC_dsubaff25, controlMPC_lub25, controlMPC_dlubaff25);
controlMPC_LA_VSUB_INDEXED_2(controlMPC_dzaff26, controlMPC_lbIdx26, controlMPC_rilb26, controlMPC_dslbaff26);
controlMPC_LA_VSUB3_2(controlMPC_llbbyslb26, controlMPC_dslbaff26, controlMPC_llb26, controlMPC_dllbaff26);
controlMPC_LA_VSUB2_INDEXED_2(controlMPC_riub26, controlMPC_dzaff26, controlMPC_ubIdx26, controlMPC_dsubaff26);
controlMPC_LA_VSUB3_2(controlMPC_lubbysub26, controlMPC_dsubaff26, controlMPC_lub26, controlMPC_dlubaff26);
controlMPC_LA_VSUB_INDEXED_2(controlMPC_dzaff27, controlMPC_lbIdx27, controlMPC_rilb27, controlMPC_dslbaff27);
controlMPC_LA_VSUB3_2(controlMPC_llbbyslb27, controlMPC_dslbaff27, controlMPC_llb27, controlMPC_dllbaff27);
controlMPC_LA_VSUB2_INDEXED_2(controlMPC_riub27, controlMPC_dzaff27, controlMPC_ubIdx27, controlMPC_dsubaff27);
controlMPC_LA_VSUB3_2(controlMPC_lubbysub27, controlMPC_dsubaff27, controlMPC_lub27, controlMPC_dlubaff27);
controlMPC_LA_VSUB_INDEXED_2(controlMPC_dzaff28, controlMPC_lbIdx28, controlMPC_rilb28, controlMPC_dslbaff28);
controlMPC_LA_VSUB3_2(controlMPC_llbbyslb28, controlMPC_dslbaff28, controlMPC_llb28, controlMPC_dllbaff28);
controlMPC_LA_VSUB2_INDEXED_2(controlMPC_riub28, controlMPC_dzaff28, controlMPC_ubIdx28, controlMPC_dsubaff28);
controlMPC_LA_VSUB3_2(controlMPC_lubbysub28, controlMPC_dsubaff28, controlMPC_lub28, controlMPC_dlubaff28);
info->lsit_aff = controlMPC_LINESEARCH_BACKTRACKING_AFFINE(controlMPC_l, controlMPC_s, controlMPC_dl_aff, controlMPC_ds_aff, &info->step_aff, &info->mu_aff);
if( info->lsit_aff == controlMPC_NOPROGRESS ){
exitcode = controlMPC_NOPROGRESS; break;
}
sigma_3rdroot = info->mu_aff / info->mu;
info->sigma = sigma_3rdroot*sigma_3rdroot*sigma_3rdroot;
musigma = info->mu * info->sigma;
controlMPC_LA_VSUB5_116(controlMPC_ds_aff, controlMPC_dl_aff, musigma, controlMPC_ccrhs);
controlMPC_LA_VSUB6_INDEXED_2_2_2(controlMPC_ccrhsub00, controlMPC_sub00, controlMPC_ubIdx00, controlMPC_ccrhsl00, controlMPC_slb00, controlMPC_lbIdx00, controlMPC_rd00);
controlMPC_LA_VSUB6_INDEXED_2_2_2(controlMPC_ccrhsub01, controlMPC_sub01, controlMPC_ubIdx01, controlMPC_ccrhsl01, controlMPC_slb01, controlMPC_lbIdx01, controlMPC_rd01);
controlMPC_LA_DIAG_FORWARDSUB_2(controlMPC_Phi00, controlMPC_rd00, controlMPC_Lbyrd00);
controlMPC_LA_DIAG_FORWARDSUB_2(controlMPC_Phi01, controlMPC_rd01, controlMPC_Lbyrd01);
controlMPC_LA_DENSE_DIAGZERO_2MVMADD_2_2_2(controlMPC_V00, controlMPC_Lbyrd00, controlMPC_W01, controlMPC_Lbyrd01, controlMPC_beta00);
controlMPC_LA_DENSE_FORWARDSUB_2(controlMPC_Ld00, controlMPC_beta00, controlMPC_yy00);
controlMPC_LA_VSUB6_INDEXED_2_2_2(controlMPC_ccrhsub02, controlMPC_sub02, controlMPC_ubIdx02, controlMPC_ccrhsl02, controlMPC_slb02, controlMPC_lbIdx02, controlMPC_rd02);
controlMPC_LA_DIAG_FORWARDSUB_2(controlMPC_Phi02, controlMPC_rd02, controlMPC_Lbyrd02);
controlMPC_LA_DENSE_DIAGZERO_2MVMADD_2_2_2(controlMPC_V01, controlMPC_Lbyrd01, controlMPC_W02, controlMPC_Lbyrd02, controlMPC_beta01);
controlMPC_LA_DENSE_MVMSUB1_2_2(controlMPC_Lsd01, controlMPC_yy00, controlMPC_beta01, controlMPC_bmy01);
controlMPC_LA_DENSE_FORWARDSUB_2(controlMPC_Ld01, controlMPC_bmy01, controlMPC_yy01);
controlMPC_LA_VSUB6_INDEXED_2_2_2(controlMPC_ccrhsub03, controlMPC_sub03, controlMPC_ubIdx03, controlMPC_ccrhsl03, controlMPC_slb03, controlMPC_lbIdx03, controlMPC_rd03);
controlMPC_LA_DIAG_FORWARDSUB_2(controlMPC_Phi03, controlMPC_rd03, controlMPC_Lbyrd03);
controlMPC_LA_DENSE_DIAGZERO_2MVMADD_2_2_2(controlMPC_V02, controlMPC_Lbyrd02, controlMPC_W03, controlMPC_Lbyrd03, controlMPC_beta02);
controlMPC_LA_DENSE_MVMSUB1_2_2(controlMPC_Lsd02, controlMPC_yy01, controlMPC_beta02, controlMPC_bmy02);
controlMPC_LA_DENSE_FORWARDSUB_2(controlMPC_Ld02, controlMPC_bmy02, controlMPC_yy02);
controlMPC_LA_VSUB6_INDEXED_2_2_2(controlMPC_ccrhsub04, controlMPC_sub04, controlMPC_ubIdx04, controlMPC_ccrhsl04, controlMPC_slb04, controlMPC_lbIdx04, controlMPC_rd04);
controlMPC_LA_DIAG_FORWARDSUB_2(controlMPC_Phi04, controlMPC_rd04, controlMPC_Lbyrd04);
controlMPC_LA_DENSE_DIAGZERO_2MVMADD_2_2_2(controlMPC_V03, controlMPC_Lbyrd03, controlMPC_W04, controlMPC_Lbyrd04, controlMPC_beta03);
controlMPC_LA_DENSE_MVMSUB1_2_2(controlMPC_Lsd03, controlMPC_yy02, controlMPC_beta03, controlMPC_bmy03);
controlMPC_LA_DENSE_FORWARDSUB_2(controlMPC_Ld03, controlMPC_bmy03, controlMPC_yy03);
controlMPC_LA_VSUB6_INDEXED_2_2_2(controlMPC_ccrhsub05, controlMPC_sub05, controlMPC_ubIdx05, controlMPC_ccrhsl05, controlMPC_slb05, controlMPC_lbIdx05, controlMPC_rd05);
controlMPC_LA_DIAG_FORWARDSUB_2(controlMPC_Phi05, controlMPC_rd05, controlMPC_Lbyrd05);
controlMPC_LA_DENSE_DIAGZERO_2MVMADD_2_2_2(controlMPC_V04, controlMPC_Lbyrd04, controlMPC_W05, controlMPC_Lbyrd05, controlMPC_beta04);
controlMPC_LA_DENSE_MVMSUB1_2_2(controlMPC_Lsd04, controlMPC_yy03, controlMPC_beta04, controlMPC_bmy04);
controlMPC_LA_DENSE_FORWARDSUB_2(controlMPC_Ld04, controlMPC_bmy04, controlMPC_yy04);
controlMPC_LA_VSUB6_INDEXED_2_2_2(controlMPC_ccrhsub06, controlMPC_sub06, controlMPC_ubIdx06, controlMPC_ccrhsl06, controlMPC_slb06, controlMPC_lbIdx06, controlMPC_rd06);
controlMPC_LA_DIAG_FORWARDSUB_2(controlMPC_Phi06, controlMPC_rd06, controlMPC_Lbyrd06);
controlMPC_LA_DENSE_DIAGZERO_2MVMADD_2_2_2(controlMPC_V05, controlMPC_Lbyrd05, controlMPC_W06, controlMPC_Lbyrd06, controlMPC_beta05);
controlMPC_LA_DENSE_MVMSUB1_2_2(controlMPC_Lsd05, controlMPC_yy04, controlMPC_beta05, controlMPC_bmy05);
controlMPC_LA_DENSE_FORWARDSUB_2(controlMPC_Ld05, controlMPC_bmy05, controlMPC_yy05);
controlMPC_LA_VSUB6_INDEXED_2_2_2(controlMPC_ccrhsub07, controlMPC_sub07, controlMPC_ubIdx07, controlMPC_ccrhsl07, controlMPC_slb07, controlMPC_lbIdx07, controlMPC_rd07);
controlMPC_LA_DIAG_FORWARDSUB_2(controlMPC_Phi07, controlMPC_rd07, controlMPC_Lbyrd07);
controlMPC_LA_DENSE_DIAGZERO_2MVMADD_2_2_2(controlMPC_V06, controlMPC_Lbyrd06, controlMPC_W07, controlMPC_Lbyrd07, controlMPC_beta06);
controlMPC_LA_DENSE_MVMSUB1_2_2(controlMPC_Lsd06, controlMPC_yy05, controlMPC_beta06, controlMPC_bmy06);
controlMPC_LA_DENSE_FORWARDSUB_2(controlMPC_Ld06, controlMPC_bmy06, controlMPC_yy06);
controlMPC_LA_VSUB6_INDEXED_2_2_2(controlMPC_ccrhsub08, controlMPC_sub08, controlMPC_ubIdx08, controlMPC_ccrhsl08, controlMPC_slb08, controlMPC_lbIdx08, controlMPC_rd08);
controlMPC_LA_DIAG_FORWARDSUB_2(controlMPC_Phi08, controlMPC_rd08, controlMPC_Lbyrd08);
controlMPC_LA_DENSE_DIAGZERO_2MVMADD_2_2_2(controlMPC_V07, controlMPC_Lbyrd07, controlMPC_W08, controlMPC_Lbyrd08, controlMPC_beta07);
controlMPC_LA_DENSE_MVMSUB1_2_2(controlMPC_Lsd07, controlMPC_yy06, controlMPC_beta07, controlMPC_bmy07);
controlMPC_LA_DENSE_FORWARDSUB_2(controlMPC_Ld07, controlMPC_bmy07, controlMPC_yy07);
controlMPC_LA_VSUB6_INDEXED_2_2_2(controlMPC_ccrhsub09, controlMPC_sub09, controlMPC_ubIdx09, controlMPC_ccrhsl09, controlMPC_slb09, controlMPC_lbIdx09, controlMPC_rd09);
controlMPC_LA_DIAG_FORWARDSUB_2(controlMPC_Phi09, controlMPC_rd09, controlMPC_Lbyrd09);
controlMPC_LA_DENSE_DIAGZERO_2MVMADD_2_2_2(controlMPC_V08, controlMPC_Lbyrd08, controlMPC_W09, controlMPC_Lbyrd09, controlMPC_beta08);
controlMPC_LA_DENSE_MVMSUB1_2_2(controlMPC_Lsd08, controlMPC_yy07, controlMPC_beta08, controlMPC_bmy08);
controlMPC_LA_DENSE_FORWARDSUB_2(controlMPC_Ld08, controlMPC_bmy08, controlMPC_yy08);
controlMPC_LA_VSUB6_INDEXED_2_2_2(controlMPC_ccrhsub10, controlMPC_sub10, controlMPC_ubIdx10, controlMPC_ccrhsl10, controlMPC_slb10, controlMPC_lbIdx10, controlMPC_rd10);
controlMPC_LA_DIAG_FORWARDSUB_2(controlMPC_Phi10, controlMPC_rd10, controlMPC_Lbyrd10);
controlMPC_LA_DENSE_DIAGZERO_2MVMADD_2_2_2(controlMPC_V09, controlMPC_Lbyrd09, controlMPC_W10, controlMPC_Lbyrd10, controlMPC_beta09);
controlMPC_LA_DENSE_MVMSUB1_2_2(controlMPC_Lsd09, controlMPC_yy08, controlMPC_beta09, controlMPC_bmy09);
controlMPC_LA_DENSE_FORWARDSUB_2(controlMPC_Ld09, controlMPC_bmy09, controlMPC_yy09);
controlMPC_LA_VSUB6_INDEXED_2_2_2(controlMPC_ccrhsub11, controlMPC_sub11, controlMPC_ubIdx11, controlMPC_ccrhsl11, controlMPC_slb11, controlMPC_lbIdx11, controlMPC_rd11);
controlMPC_LA_DIAG_FORWARDSUB_2(controlMPC_Phi11, controlMPC_rd11, controlMPC_Lbyrd11);
controlMPC_LA_DENSE_DIAGZERO_2MVMADD_2_2_2(controlMPC_V10, controlMPC_Lbyrd10, controlMPC_W11, controlMPC_Lbyrd11, controlMPC_beta10);
controlMPC_LA_DENSE_MVMSUB1_2_2(controlMPC_Lsd10, controlMPC_yy09, controlMPC_beta10, controlMPC_bmy10);
controlMPC_LA_DENSE_FORWARDSUB_2(controlMPC_Ld10, controlMPC_bmy10, controlMPC_yy10);
controlMPC_LA_VSUB6_INDEXED_2_2_2(controlMPC_ccrhsub12, controlMPC_sub12, controlMPC_ubIdx12, controlMPC_ccrhsl12, controlMPC_slb12, controlMPC_lbIdx12, controlMPC_rd12);
controlMPC_LA_DIAG_FORWARDSUB_2(controlMPC_Phi12, controlMPC_rd12, controlMPC_Lbyrd12);
controlMPC_LA_DENSE_DIAGZERO_2MVMADD_2_2_2(controlMPC_V11, controlMPC_Lbyrd11, controlMPC_W12, controlMPC_Lbyrd12, controlMPC_beta11);
controlMPC_LA_DENSE_MVMSUB1_2_2(controlMPC_Lsd11, controlMPC_yy10, controlMPC_beta11, controlMPC_bmy11);
controlMPC_LA_DENSE_FORWARDSUB_2(controlMPC_Ld11, controlMPC_bmy11, controlMPC_yy11);
controlMPC_LA_VSUB6_INDEXED_2_2_2(controlMPC_ccrhsub13, controlMPC_sub13, controlMPC_ubIdx13, controlMPC_ccrhsl13, controlMPC_slb13, controlMPC_lbIdx13, controlMPC_rd13);
controlMPC_LA_DIAG_FORWARDSUB_2(controlMPC_Phi13, controlMPC_rd13, controlMPC_Lbyrd13);
controlMPC_LA_DENSE_DIAGZERO_2MVMADD_2_2_2(controlMPC_V12, controlMPC_Lbyrd12, controlMPC_W13, controlMPC_Lbyrd13, controlMPC_beta12);
controlMPC_LA_DENSE_MVMSUB1_2_2(controlMPC_Lsd12, controlMPC_yy11, controlMPC_beta12, controlMPC_bmy12);
controlMPC_LA_DENSE_FORWARDSUB_2(controlMPC_Ld12, controlMPC_bmy12, controlMPC_yy12);
controlMPC_LA_VSUB6_INDEXED_2_2_2(controlMPC_ccrhsub14, controlMPC_sub14, controlMPC_ubIdx14, controlMPC_ccrhsl14, controlMPC_slb14, controlMPC_lbIdx14, controlMPC_rd14);
controlMPC_LA_DIAG_FORWARDSUB_2(controlMPC_Phi14, controlMPC_rd14, controlMPC_Lbyrd14);
controlMPC_LA_DENSE_DIAGZERO_2MVMADD_2_2_2(controlMPC_V13, controlMPC_Lbyrd13, controlMPC_W14, controlMPC_Lbyrd14, controlMPC_beta13);
controlMPC_LA_DENSE_MVMSUB1_2_2(controlMPC_Lsd13, controlMPC_yy12, controlMPC_beta13, controlMPC_bmy13);
controlMPC_LA_DENSE_FORWARDSUB_2(controlMPC_Ld13, controlMPC_bmy13, controlMPC_yy13);
controlMPC_LA_VSUB6_INDEXED_2_2_2(controlMPC_ccrhsub15, controlMPC_sub15, controlMPC_ubIdx15, controlMPC_ccrhsl15, controlMPC_slb15, controlMPC_lbIdx15, controlMPC_rd15);
controlMPC_LA_DIAG_FORWARDSUB_2(controlMPC_Phi15, controlMPC_rd15, controlMPC_Lbyrd15);
controlMPC_LA_DENSE_DIAGZERO_2MVMADD_2_2_2(controlMPC_V14, controlMPC_Lbyrd14, controlMPC_W15, controlMPC_Lbyrd15, controlMPC_beta14);
controlMPC_LA_DENSE_MVMSUB1_2_2(controlMPC_Lsd14, controlMPC_yy13, controlMPC_beta14, controlMPC_bmy14);
controlMPC_LA_DENSE_FORWARDSUB_2(controlMPC_Ld14, controlMPC_bmy14, controlMPC_yy14);
controlMPC_LA_VSUB6_INDEXED_2_2_2(controlMPC_ccrhsub16, controlMPC_sub16, controlMPC_ubIdx16, controlMPC_ccrhsl16, controlMPC_slb16, controlMPC_lbIdx16, controlMPC_rd16);
controlMPC_LA_DIAG_FORWARDSUB_2(controlMPC_Phi16, controlMPC_rd16, controlMPC_Lbyrd16);
controlMPC_LA_DENSE_DIAGZERO_2MVMADD_2_2_2(controlMPC_V15, controlMPC_Lbyrd15, controlMPC_W16, controlMPC_Lbyrd16, controlMPC_beta15);
controlMPC_LA_DENSE_MVMSUB1_2_2(controlMPC_Lsd15, controlMPC_yy14, controlMPC_beta15, controlMPC_bmy15);
controlMPC_LA_DENSE_FORWARDSUB_2(controlMPC_Ld15, controlMPC_bmy15, controlMPC_yy15);
controlMPC_LA_VSUB6_INDEXED_2_2_2(controlMPC_ccrhsub17, controlMPC_sub17, controlMPC_ubIdx17, controlMPC_ccrhsl17, controlMPC_slb17, controlMPC_lbIdx17, controlMPC_rd17);
controlMPC_LA_DIAG_FORWARDSUB_2(controlMPC_Phi17, controlMPC_rd17, controlMPC_Lbyrd17);
controlMPC_LA_DENSE_DIAGZERO_2MVMADD_2_2_2(controlMPC_V16, controlMPC_Lbyrd16, controlMPC_W17, controlMPC_Lbyrd17, controlMPC_beta16);
controlMPC_LA_DENSE_MVMSUB1_2_2(controlMPC_Lsd16, controlMPC_yy15, controlMPC_beta16, controlMPC_bmy16);
controlMPC_LA_DENSE_FORWARDSUB_2(controlMPC_Ld16, controlMPC_bmy16, controlMPC_yy16);
controlMPC_LA_VSUB6_INDEXED_2_2_2(controlMPC_ccrhsub18, controlMPC_sub18, controlMPC_ubIdx18, controlMPC_ccrhsl18, controlMPC_slb18, controlMPC_lbIdx18, controlMPC_rd18);
controlMPC_LA_DIAG_FORWARDSUB_2(controlMPC_Phi18, controlMPC_rd18, controlMPC_Lbyrd18);
controlMPC_LA_DENSE_DIAGZERO_2MVMADD_2_2_2(controlMPC_V17, controlMPC_Lbyrd17, controlMPC_W18, controlMPC_Lbyrd18, controlMPC_beta17);
controlMPC_LA_DENSE_MVMSUB1_2_2(controlMPC_Lsd17, controlMPC_yy16, controlMPC_beta17, controlMPC_bmy17);
controlMPC_LA_DENSE_FORWARDSUB_2(controlMPC_Ld17, controlMPC_bmy17, controlMPC_yy17);
controlMPC_LA_VSUB6_INDEXED_2_2_2(controlMPC_ccrhsub19, controlMPC_sub19, controlMPC_ubIdx19, controlMPC_ccrhsl19, controlMPC_slb19, controlMPC_lbIdx19, controlMPC_rd19);
controlMPC_LA_DIAG_FORWARDSUB_2(controlMPC_Phi19, controlMPC_rd19, controlMPC_Lbyrd19);
controlMPC_LA_DENSE_DIAGZERO_2MVMADD_2_2_2(controlMPC_V18, controlMPC_Lbyrd18, controlMPC_W19, controlMPC_Lbyrd19, controlMPC_beta18);
controlMPC_LA_DENSE_MVMSUB1_2_2(controlMPC_Lsd18, controlMPC_yy17, controlMPC_beta18, controlMPC_bmy18);
controlMPC_LA_DENSE_FORWARDSUB_2(controlMPC_Ld18, controlMPC_bmy18, controlMPC_yy18);
controlMPC_LA_VSUB6_INDEXED_2_2_2(controlMPC_ccrhsub20, controlMPC_sub20, controlMPC_ubIdx20, controlMPC_ccrhsl20, controlMPC_slb20, controlMPC_lbIdx20, controlMPC_rd20);
controlMPC_LA_DIAG_FORWARDSUB_2(controlMPC_Phi20, controlMPC_rd20, controlMPC_Lbyrd20);
controlMPC_LA_DENSE_DIAGZERO_2MVMADD_2_2_2(controlMPC_V19, controlMPC_Lbyrd19, controlMPC_W20, controlMPC_Lbyrd20, controlMPC_beta19);
controlMPC_LA_DENSE_MVMSUB1_2_2(controlMPC_Lsd19, controlMPC_yy18, controlMPC_beta19, controlMPC_bmy19);
controlMPC_LA_DENSE_FORWARDSUB_2(controlMPC_Ld19, controlMPC_bmy19, controlMPC_yy19);
controlMPC_LA_VSUB6_INDEXED_2_2_2(controlMPC_ccrhsub21, controlMPC_sub21, controlMPC_ubIdx21, controlMPC_ccrhsl21, controlMPC_slb21, controlMPC_lbIdx21, controlMPC_rd21);
controlMPC_LA_DIAG_FORWARDSUB_2(controlMPC_Phi21, controlMPC_rd21, controlMPC_Lbyrd21);
controlMPC_LA_DENSE_DIAGZERO_2MVMADD_2_2_2(controlMPC_V20, controlMPC_Lbyrd20, controlMPC_W21, controlMPC_Lbyrd21, controlMPC_beta20);
controlMPC_LA_DENSE_MVMSUB1_2_2(controlMPC_Lsd20, controlMPC_yy19, controlMPC_beta20, controlMPC_bmy20);
controlMPC_LA_DENSE_FORWARDSUB_2(controlMPC_Ld20, controlMPC_bmy20, controlMPC_yy20);
controlMPC_LA_VSUB6_INDEXED_2_2_2(controlMPC_ccrhsub22, controlMPC_sub22, controlMPC_ubIdx22, controlMPC_ccrhsl22, controlMPC_slb22, controlMPC_lbIdx22, controlMPC_rd22);
controlMPC_LA_DIAG_FORWARDSUB_2(controlMPC_Phi22, controlMPC_rd22, controlMPC_Lbyrd22);
controlMPC_LA_DENSE_DIAGZERO_2MVMADD_2_2_2(controlMPC_V21, controlMPC_Lbyrd21, controlMPC_W22, controlMPC_Lbyrd22, controlMPC_beta21);
controlMPC_LA_DENSE_MVMSUB1_2_2(controlMPC_Lsd21, controlMPC_yy20, controlMPC_beta21, controlMPC_bmy21);
controlMPC_LA_DENSE_FORWARDSUB_2(controlMPC_Ld21, controlMPC_bmy21, controlMPC_yy21);
controlMPC_LA_VSUB6_INDEXED_2_2_2(controlMPC_ccrhsub23, controlMPC_sub23, controlMPC_ubIdx23, controlMPC_ccrhsl23, controlMPC_slb23, controlMPC_lbIdx23, controlMPC_rd23);
controlMPC_LA_DIAG_FORWARDSUB_2(controlMPC_Phi23, controlMPC_rd23, controlMPC_Lbyrd23);
controlMPC_LA_DENSE_DIAGZERO_2MVMADD_2_2_2(controlMPC_V22, controlMPC_Lbyrd22, controlMPC_W23, controlMPC_Lbyrd23, controlMPC_beta22);
controlMPC_LA_DENSE_MVMSUB1_2_2(controlMPC_Lsd22, controlMPC_yy21, controlMPC_beta22, controlMPC_bmy22);
controlMPC_LA_DENSE_FORWARDSUB_2(controlMPC_Ld22, controlMPC_bmy22, controlMPC_yy22);
controlMPC_LA_VSUB6_INDEXED_2_2_2(controlMPC_ccrhsub24, controlMPC_sub24, controlMPC_ubIdx24, controlMPC_ccrhsl24, controlMPC_slb24, controlMPC_lbIdx24, controlMPC_rd24);
controlMPC_LA_DIAG_FORWARDSUB_2(controlMPC_Phi24, controlMPC_rd24, controlMPC_Lbyrd24);
controlMPC_LA_DENSE_DIAGZERO_2MVMADD_2_2_2(controlMPC_V23, controlMPC_Lbyrd23, controlMPC_W24, controlMPC_Lbyrd24, controlMPC_beta23);
controlMPC_LA_DENSE_MVMSUB1_2_2(controlMPC_Lsd23, controlMPC_yy22, controlMPC_beta23, controlMPC_bmy23);
controlMPC_LA_DENSE_FORWARDSUB_2(controlMPC_Ld23, controlMPC_bmy23, controlMPC_yy23);
controlMPC_LA_VSUB6_INDEXED_2_2_2(controlMPC_ccrhsub25, controlMPC_sub25, controlMPC_ubIdx25, controlMPC_ccrhsl25, controlMPC_slb25, controlMPC_lbIdx25, controlMPC_rd25);
controlMPC_LA_DIAG_FORWARDSUB_2(controlMPC_Phi25, controlMPC_rd25, controlMPC_Lbyrd25);
controlMPC_LA_DENSE_DIAGZERO_2MVMADD_2_2_2(controlMPC_V24, controlMPC_Lbyrd24, controlMPC_W25, controlMPC_Lbyrd25, controlMPC_beta24);
controlMPC_LA_DENSE_MVMSUB1_2_2(controlMPC_Lsd24, controlMPC_yy23, controlMPC_beta24, controlMPC_bmy24);
controlMPC_LA_DENSE_FORWARDSUB_2(controlMPC_Ld24, controlMPC_bmy24, controlMPC_yy24);
controlMPC_LA_VSUB6_INDEXED_2_2_2(controlMPC_ccrhsub26, controlMPC_sub26, controlMPC_ubIdx26, controlMPC_ccrhsl26, controlMPC_slb26, controlMPC_lbIdx26, controlMPC_rd26);
controlMPC_LA_DIAG_FORWARDSUB_2(controlMPC_Phi26, controlMPC_rd26, controlMPC_Lbyrd26);
controlMPC_LA_DENSE_DIAGZERO_2MVMADD_2_2_2(controlMPC_V25, controlMPC_Lbyrd25, controlMPC_W26, controlMPC_Lbyrd26, controlMPC_beta25);
controlMPC_LA_DENSE_MVMSUB1_2_2(controlMPC_Lsd25, controlMPC_yy24, controlMPC_beta25, controlMPC_bmy25);
controlMPC_LA_DENSE_FORWARDSUB_2(controlMPC_Ld25, controlMPC_bmy25, controlMPC_yy25);
controlMPC_LA_VSUB6_INDEXED_2_2_2(controlMPC_ccrhsub27, controlMPC_sub27, controlMPC_ubIdx27, controlMPC_ccrhsl27, controlMPC_slb27, controlMPC_lbIdx27, controlMPC_rd27);
controlMPC_LA_DIAG_FORWARDSUB_2(controlMPC_Phi27, controlMPC_rd27, controlMPC_Lbyrd27);
controlMPC_LA_DENSE_DIAGZERO_2MVMADD_2_2_2(controlMPC_V26, controlMPC_Lbyrd26, controlMPC_W27, controlMPC_Lbyrd27, controlMPC_beta26);
controlMPC_LA_DENSE_MVMSUB1_2_2(controlMPC_Lsd26, controlMPC_yy25, controlMPC_beta26, controlMPC_bmy26);
controlMPC_LA_DENSE_FORWARDSUB_2(controlMPC_Ld26, controlMPC_bmy26, controlMPC_yy26);
controlMPC_LA_VSUB6_INDEXED_2_2_2(controlMPC_ccrhsub28, controlMPC_sub28, controlMPC_ubIdx28, controlMPC_ccrhsl28, controlMPC_slb28, controlMPC_lbIdx28, controlMPC_rd28);
controlMPC_LA_DIAG_FORWARDSUB_2(controlMPC_Phi28, controlMPC_rd28, controlMPC_Lbyrd28);
controlMPC_LA_DENSE_DIAGZERO_2MVMADD_2_2_2(controlMPC_V27, controlMPC_Lbyrd27, controlMPC_W28, controlMPC_Lbyrd28, controlMPC_beta27);
controlMPC_LA_DENSE_MVMSUB1_2_2(controlMPC_Lsd27, controlMPC_yy26, controlMPC_beta27, controlMPC_bmy27);
controlMPC_LA_DENSE_FORWARDSUB_2(controlMPC_Ld27, controlMPC_bmy27, controlMPC_yy27);
controlMPC_LA_DENSE_BACKWARDSUB_2(controlMPC_Ld27, controlMPC_yy27, controlMPC_dvcc27);
controlMPC_LA_DENSE_MTVMSUB_2_2(controlMPC_Lsd27, controlMPC_dvcc27, controlMPC_yy26, controlMPC_bmy26);
controlMPC_LA_DENSE_BACKWARDSUB_2(controlMPC_Ld26, controlMPC_bmy26, controlMPC_dvcc26);
controlMPC_LA_DENSE_MTVMSUB_2_2(controlMPC_Lsd26, controlMPC_dvcc26, controlMPC_yy25, controlMPC_bmy25);
controlMPC_LA_DENSE_BACKWARDSUB_2(controlMPC_Ld25, controlMPC_bmy25, controlMPC_dvcc25);
controlMPC_LA_DENSE_MTVMSUB_2_2(controlMPC_Lsd25, controlMPC_dvcc25, controlMPC_yy24, controlMPC_bmy24);
controlMPC_LA_DENSE_BACKWARDSUB_2(controlMPC_Ld24, controlMPC_bmy24, controlMPC_dvcc24);
controlMPC_LA_DENSE_MTVMSUB_2_2(controlMPC_Lsd24, controlMPC_dvcc24, controlMPC_yy23, controlMPC_bmy23);
controlMPC_LA_DENSE_BACKWARDSUB_2(controlMPC_Ld23, controlMPC_bmy23, controlMPC_dvcc23);
controlMPC_LA_DENSE_MTVMSUB_2_2(controlMPC_Lsd23, controlMPC_dvcc23, controlMPC_yy22, controlMPC_bmy22);
controlMPC_LA_DENSE_BACKWARDSUB_2(controlMPC_Ld22, controlMPC_bmy22, controlMPC_dvcc22);
controlMPC_LA_DENSE_MTVMSUB_2_2(controlMPC_Lsd22, controlMPC_dvcc22, controlMPC_yy21, controlMPC_bmy21);
controlMPC_LA_DENSE_BACKWARDSUB_2(controlMPC_Ld21, controlMPC_bmy21, controlMPC_dvcc21);
controlMPC_LA_DENSE_MTVMSUB_2_2(controlMPC_Lsd21, controlMPC_dvcc21, controlMPC_yy20, controlMPC_bmy20);
controlMPC_LA_DENSE_BACKWARDSUB_2(controlMPC_Ld20, controlMPC_bmy20, controlMPC_dvcc20);
controlMPC_LA_DENSE_MTVMSUB_2_2(controlMPC_Lsd20, controlMPC_dvcc20, controlMPC_yy19, controlMPC_bmy19);
controlMPC_LA_DENSE_BACKWARDSUB_2(controlMPC_Ld19, controlMPC_bmy19, controlMPC_dvcc19);
controlMPC_LA_DENSE_MTVMSUB_2_2(controlMPC_Lsd19, controlMPC_dvcc19, controlMPC_yy18, controlMPC_bmy18);
controlMPC_LA_DENSE_BACKWARDSUB_2(controlMPC_Ld18, controlMPC_bmy18, controlMPC_dvcc18);
controlMPC_LA_DENSE_MTVMSUB_2_2(controlMPC_Lsd18, controlMPC_dvcc18, controlMPC_yy17, controlMPC_bmy17);
controlMPC_LA_DENSE_BACKWARDSUB_2(controlMPC_Ld17, controlMPC_bmy17, controlMPC_dvcc17);
controlMPC_LA_DENSE_MTVMSUB_2_2(controlMPC_Lsd17, controlMPC_dvcc17, controlMPC_yy16, controlMPC_bmy16);
controlMPC_LA_DENSE_BACKWARDSUB_2(controlMPC_Ld16, controlMPC_bmy16, controlMPC_dvcc16);
controlMPC_LA_DENSE_MTVMSUB_2_2(controlMPC_Lsd16, controlMPC_dvcc16, controlMPC_yy15, controlMPC_bmy15);
controlMPC_LA_DENSE_BACKWARDSUB_2(controlMPC_Ld15, controlMPC_bmy15, controlMPC_dvcc15);
controlMPC_LA_DENSE_MTVMSUB_2_2(controlMPC_Lsd15, controlMPC_dvcc15, controlMPC_yy14, controlMPC_bmy14);
controlMPC_LA_DENSE_BACKWARDSUB_2(controlMPC_Ld14, controlMPC_bmy14, controlMPC_dvcc14);
controlMPC_LA_DENSE_MTVMSUB_2_2(controlMPC_Lsd14, controlMPC_dvcc14, controlMPC_yy13, controlMPC_bmy13);
controlMPC_LA_DENSE_BACKWARDSUB_2(controlMPC_Ld13, controlMPC_bmy13, controlMPC_dvcc13);
controlMPC_LA_DENSE_MTVMSUB_2_2(controlMPC_Lsd13, controlMPC_dvcc13, controlMPC_yy12, controlMPC_bmy12);
controlMPC_LA_DENSE_BACKWARDSUB_2(controlMPC_Ld12, controlMPC_bmy12, controlMPC_dvcc12);
controlMPC_LA_DENSE_MTVMSUB_2_2(controlMPC_Lsd12, controlMPC_dvcc12, controlMPC_yy11, controlMPC_bmy11);
controlMPC_LA_DENSE_BACKWARDSUB_2(controlMPC_Ld11, controlMPC_bmy11, controlMPC_dvcc11);
controlMPC_LA_DENSE_MTVMSUB_2_2(controlMPC_Lsd11, controlMPC_dvcc11, controlMPC_yy10, controlMPC_bmy10);
controlMPC_LA_DENSE_BACKWARDSUB_2(controlMPC_Ld10, controlMPC_bmy10, controlMPC_dvcc10);
controlMPC_LA_DENSE_MTVMSUB_2_2(controlMPC_Lsd10, controlMPC_dvcc10, controlMPC_yy09, controlMPC_bmy09);
controlMPC_LA_DENSE_BACKWARDSUB_2(controlMPC_Ld09, controlMPC_bmy09, controlMPC_dvcc09);
controlMPC_LA_DENSE_MTVMSUB_2_2(controlMPC_Lsd09, controlMPC_dvcc09, controlMPC_yy08, controlMPC_bmy08);
controlMPC_LA_DENSE_BACKWARDSUB_2(controlMPC_Ld08, controlMPC_bmy08, controlMPC_dvcc08);
controlMPC_LA_DENSE_MTVMSUB_2_2(controlMPC_Lsd08, controlMPC_dvcc08, controlMPC_yy07, controlMPC_bmy07);
controlMPC_LA_DENSE_BACKWARDSUB_2(controlMPC_Ld07, controlMPC_bmy07, controlMPC_dvcc07);
controlMPC_LA_DENSE_MTVMSUB_2_2(controlMPC_Lsd07, controlMPC_dvcc07, controlMPC_yy06, controlMPC_bmy06);
controlMPC_LA_DENSE_BACKWARDSUB_2(controlMPC_Ld06, controlMPC_bmy06, controlMPC_dvcc06);
controlMPC_LA_DENSE_MTVMSUB_2_2(controlMPC_Lsd06, controlMPC_dvcc06, controlMPC_yy05, controlMPC_bmy05);
controlMPC_LA_DENSE_BACKWARDSUB_2(controlMPC_Ld05, controlMPC_bmy05, controlMPC_dvcc05);
controlMPC_LA_DENSE_MTVMSUB_2_2(controlMPC_Lsd05, controlMPC_dvcc05, controlMPC_yy04, controlMPC_bmy04);
controlMPC_LA_DENSE_BACKWARDSUB_2(controlMPC_Ld04, controlMPC_bmy04, controlMPC_dvcc04);
controlMPC_LA_DENSE_MTVMSUB_2_2(controlMPC_Lsd04, controlMPC_dvcc04, controlMPC_yy03, controlMPC_bmy03);
controlMPC_LA_DENSE_BACKWARDSUB_2(controlMPC_Ld03, controlMPC_bmy03, controlMPC_dvcc03);
controlMPC_LA_DENSE_MTVMSUB_2_2(controlMPC_Lsd03, controlMPC_dvcc03, controlMPC_yy02, controlMPC_bmy02);
controlMPC_LA_DENSE_BACKWARDSUB_2(controlMPC_Ld02, controlMPC_bmy02, controlMPC_dvcc02);
controlMPC_LA_DENSE_MTVMSUB_2_2(controlMPC_Lsd02, controlMPC_dvcc02, controlMPC_yy01, controlMPC_bmy01);
controlMPC_LA_DENSE_BACKWARDSUB_2(controlMPC_Ld01, controlMPC_bmy01, controlMPC_dvcc01);
controlMPC_LA_DENSE_MTVMSUB_2_2(controlMPC_Lsd01, controlMPC_dvcc01, controlMPC_yy00, controlMPC_bmy00);
controlMPC_LA_DENSE_BACKWARDSUB_2(controlMPC_Ld00, controlMPC_bmy00, controlMPC_dvcc00);
controlMPC_LA_DENSE_MTVM_2_2(controlMPC_C00, controlMPC_dvcc00, controlMPC_grad_eq00);
controlMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlMPC_C00, controlMPC_dvcc01, controlMPC_D01, controlMPC_dvcc00, controlMPC_grad_eq01);
controlMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlMPC_C00, controlMPC_dvcc02, controlMPC_D01, controlMPC_dvcc01, controlMPC_grad_eq02);
controlMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlMPC_C00, controlMPC_dvcc03, controlMPC_D01, controlMPC_dvcc02, controlMPC_grad_eq03);
controlMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlMPC_C00, controlMPC_dvcc04, controlMPC_D01, controlMPC_dvcc03, controlMPC_grad_eq04);
controlMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlMPC_C00, controlMPC_dvcc05, controlMPC_D01, controlMPC_dvcc04, controlMPC_grad_eq05);
controlMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlMPC_C00, controlMPC_dvcc06, controlMPC_D01, controlMPC_dvcc05, controlMPC_grad_eq06);
controlMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlMPC_C00, controlMPC_dvcc07, controlMPC_D01, controlMPC_dvcc06, controlMPC_grad_eq07);
controlMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlMPC_C00, controlMPC_dvcc08, controlMPC_D01, controlMPC_dvcc07, controlMPC_grad_eq08);
controlMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlMPC_C00, controlMPC_dvcc09, controlMPC_D01, controlMPC_dvcc08, controlMPC_grad_eq09);
controlMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlMPC_C00, controlMPC_dvcc10, controlMPC_D01, controlMPC_dvcc09, controlMPC_grad_eq10);
controlMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlMPC_C00, controlMPC_dvcc11, controlMPC_D01, controlMPC_dvcc10, controlMPC_grad_eq11);
controlMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlMPC_C00, controlMPC_dvcc12, controlMPC_D01, controlMPC_dvcc11, controlMPC_grad_eq12);
controlMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlMPC_C00, controlMPC_dvcc13, controlMPC_D01, controlMPC_dvcc12, controlMPC_grad_eq13);
controlMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlMPC_C00, controlMPC_dvcc14, controlMPC_D01, controlMPC_dvcc13, controlMPC_grad_eq14);
controlMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlMPC_C00, controlMPC_dvcc15, controlMPC_D01, controlMPC_dvcc14, controlMPC_grad_eq15);
controlMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlMPC_C00, controlMPC_dvcc16, controlMPC_D01, controlMPC_dvcc15, controlMPC_grad_eq16);
controlMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlMPC_C00, controlMPC_dvcc17, controlMPC_D01, controlMPC_dvcc16, controlMPC_grad_eq17);
controlMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlMPC_C00, controlMPC_dvcc18, controlMPC_D01, controlMPC_dvcc17, controlMPC_grad_eq18);
controlMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlMPC_C00, controlMPC_dvcc19, controlMPC_D01, controlMPC_dvcc18, controlMPC_grad_eq19);
controlMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlMPC_C00, controlMPC_dvcc20, controlMPC_D01, controlMPC_dvcc19, controlMPC_grad_eq20);
controlMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlMPC_C00, controlMPC_dvcc21, controlMPC_D01, controlMPC_dvcc20, controlMPC_grad_eq21);
controlMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlMPC_C00, controlMPC_dvcc22, controlMPC_D01, controlMPC_dvcc21, controlMPC_grad_eq22);
controlMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlMPC_C00, controlMPC_dvcc23, controlMPC_D01, controlMPC_dvcc22, controlMPC_grad_eq23);
controlMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlMPC_C00, controlMPC_dvcc24, controlMPC_D01, controlMPC_dvcc23, controlMPC_grad_eq24);
controlMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlMPC_C00, controlMPC_dvcc25, controlMPC_D01, controlMPC_dvcc24, controlMPC_grad_eq25);
controlMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlMPC_C00, controlMPC_dvcc26, controlMPC_D01, controlMPC_dvcc25, controlMPC_grad_eq26);
controlMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlMPC_C00, controlMPC_dvcc27, controlMPC_D01, controlMPC_dvcc26, controlMPC_grad_eq27);
controlMPC_LA_DIAGZERO_MTVM_2_2(controlMPC_D01, controlMPC_dvcc27, controlMPC_grad_eq28);
controlMPC_LA_VSUB_58(controlMPC_rd, controlMPC_grad_eq, controlMPC_rd);
controlMPC_LA_DIAG_FORWARDBACKWARDSUB_2(controlMPC_Phi00, controlMPC_rd00, controlMPC_dzcc00);
controlMPC_LA_DIAG_FORWARDBACKWARDSUB_2(controlMPC_Phi01, controlMPC_rd01, controlMPC_dzcc01);
controlMPC_LA_DIAG_FORWARDBACKWARDSUB_2(controlMPC_Phi02, controlMPC_rd02, controlMPC_dzcc02);
controlMPC_LA_DIAG_FORWARDBACKWARDSUB_2(controlMPC_Phi03, controlMPC_rd03, controlMPC_dzcc03);
controlMPC_LA_DIAG_FORWARDBACKWARDSUB_2(controlMPC_Phi04, controlMPC_rd04, controlMPC_dzcc04);
controlMPC_LA_DIAG_FORWARDBACKWARDSUB_2(controlMPC_Phi05, controlMPC_rd05, controlMPC_dzcc05);
controlMPC_LA_DIAG_FORWARDBACKWARDSUB_2(controlMPC_Phi06, controlMPC_rd06, controlMPC_dzcc06);
controlMPC_LA_DIAG_FORWARDBACKWARDSUB_2(controlMPC_Phi07, controlMPC_rd07, controlMPC_dzcc07);
controlMPC_LA_DIAG_FORWARDBACKWARDSUB_2(controlMPC_Phi08, controlMPC_rd08, controlMPC_dzcc08);
controlMPC_LA_DIAG_FORWARDBACKWARDSUB_2(controlMPC_Phi09, controlMPC_rd09, controlMPC_dzcc09);
controlMPC_LA_DIAG_FORWARDBACKWARDSUB_2(controlMPC_Phi10, controlMPC_rd10, controlMPC_dzcc10);
controlMPC_LA_DIAG_FORWARDBACKWARDSUB_2(controlMPC_Phi11, controlMPC_rd11, controlMPC_dzcc11);
controlMPC_LA_DIAG_FORWARDBACKWARDSUB_2(controlMPC_Phi12, controlMPC_rd12, controlMPC_dzcc12);
controlMPC_LA_DIAG_FORWARDBACKWARDSUB_2(controlMPC_Phi13, controlMPC_rd13, controlMPC_dzcc13);
controlMPC_LA_DIAG_FORWARDBACKWARDSUB_2(controlMPC_Phi14, controlMPC_rd14, controlMPC_dzcc14);
controlMPC_LA_DIAG_FORWARDBACKWARDSUB_2(controlMPC_Phi15, controlMPC_rd15, controlMPC_dzcc15);
controlMPC_LA_DIAG_FORWARDBACKWARDSUB_2(controlMPC_Phi16, controlMPC_rd16, controlMPC_dzcc16);
controlMPC_LA_DIAG_FORWARDBACKWARDSUB_2(controlMPC_Phi17, controlMPC_rd17, controlMPC_dzcc17);
controlMPC_LA_DIAG_FORWARDBACKWARDSUB_2(controlMPC_Phi18, controlMPC_rd18, controlMPC_dzcc18);
controlMPC_LA_DIAG_FORWARDBACKWARDSUB_2(controlMPC_Phi19, controlMPC_rd19, controlMPC_dzcc19);
controlMPC_LA_DIAG_FORWARDBACKWARDSUB_2(controlMPC_Phi20, controlMPC_rd20, controlMPC_dzcc20);
controlMPC_LA_DIAG_FORWARDBACKWARDSUB_2(controlMPC_Phi21, controlMPC_rd21, controlMPC_dzcc21);
controlMPC_LA_DIAG_FORWARDBACKWARDSUB_2(controlMPC_Phi22, controlMPC_rd22, controlMPC_dzcc22);
controlMPC_LA_DIAG_FORWARDBACKWARDSUB_2(controlMPC_Phi23, controlMPC_rd23, controlMPC_dzcc23);
controlMPC_LA_DIAG_FORWARDBACKWARDSUB_2(controlMPC_Phi24, controlMPC_rd24, controlMPC_dzcc24);
controlMPC_LA_DIAG_FORWARDBACKWARDSUB_2(controlMPC_Phi25, controlMPC_rd25, controlMPC_dzcc25);
controlMPC_LA_DIAG_FORWARDBACKWARDSUB_2(controlMPC_Phi26, controlMPC_rd26, controlMPC_dzcc26);
controlMPC_LA_DIAG_FORWARDBACKWARDSUB_2(controlMPC_Phi27, controlMPC_rd27, controlMPC_dzcc27);
controlMPC_LA_DIAG_FORWARDBACKWARDSUB_2(controlMPC_Phi28, controlMPC_rd28, controlMPC_dzcc28);
controlMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_2(controlMPC_ccrhsl00, controlMPC_slb00, controlMPC_llbbyslb00, controlMPC_dzcc00, controlMPC_lbIdx00, controlMPC_dllbcc00);
controlMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_2(controlMPC_ccrhsub00, controlMPC_sub00, controlMPC_lubbysub00, controlMPC_dzcc00, controlMPC_ubIdx00, controlMPC_dlubcc00);
controlMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_2(controlMPC_ccrhsl01, controlMPC_slb01, controlMPC_llbbyslb01, controlMPC_dzcc01, controlMPC_lbIdx01, controlMPC_dllbcc01);
controlMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_2(controlMPC_ccrhsub01, controlMPC_sub01, controlMPC_lubbysub01, controlMPC_dzcc01, controlMPC_ubIdx01, controlMPC_dlubcc01);
controlMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_2(controlMPC_ccrhsl02, controlMPC_slb02, controlMPC_llbbyslb02, controlMPC_dzcc02, controlMPC_lbIdx02, controlMPC_dllbcc02);
controlMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_2(controlMPC_ccrhsub02, controlMPC_sub02, controlMPC_lubbysub02, controlMPC_dzcc02, controlMPC_ubIdx02, controlMPC_dlubcc02);
controlMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_2(controlMPC_ccrhsl03, controlMPC_slb03, controlMPC_llbbyslb03, controlMPC_dzcc03, controlMPC_lbIdx03, controlMPC_dllbcc03);
controlMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_2(controlMPC_ccrhsub03, controlMPC_sub03, controlMPC_lubbysub03, controlMPC_dzcc03, controlMPC_ubIdx03, controlMPC_dlubcc03);
controlMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_2(controlMPC_ccrhsl04, controlMPC_slb04, controlMPC_llbbyslb04, controlMPC_dzcc04, controlMPC_lbIdx04, controlMPC_dllbcc04);
controlMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_2(controlMPC_ccrhsub04, controlMPC_sub04, controlMPC_lubbysub04, controlMPC_dzcc04, controlMPC_ubIdx04, controlMPC_dlubcc04);
controlMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_2(controlMPC_ccrhsl05, controlMPC_slb05, controlMPC_llbbyslb05, controlMPC_dzcc05, controlMPC_lbIdx05, controlMPC_dllbcc05);
controlMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_2(controlMPC_ccrhsub05, controlMPC_sub05, controlMPC_lubbysub05, controlMPC_dzcc05, controlMPC_ubIdx05, controlMPC_dlubcc05);
controlMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_2(controlMPC_ccrhsl06, controlMPC_slb06, controlMPC_llbbyslb06, controlMPC_dzcc06, controlMPC_lbIdx06, controlMPC_dllbcc06);
controlMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_2(controlMPC_ccrhsub06, controlMPC_sub06, controlMPC_lubbysub06, controlMPC_dzcc06, controlMPC_ubIdx06, controlMPC_dlubcc06);
controlMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_2(controlMPC_ccrhsl07, controlMPC_slb07, controlMPC_llbbyslb07, controlMPC_dzcc07, controlMPC_lbIdx07, controlMPC_dllbcc07);
controlMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_2(controlMPC_ccrhsub07, controlMPC_sub07, controlMPC_lubbysub07, controlMPC_dzcc07, controlMPC_ubIdx07, controlMPC_dlubcc07);
controlMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_2(controlMPC_ccrhsl08, controlMPC_slb08, controlMPC_llbbyslb08, controlMPC_dzcc08, controlMPC_lbIdx08, controlMPC_dllbcc08);
controlMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_2(controlMPC_ccrhsub08, controlMPC_sub08, controlMPC_lubbysub08, controlMPC_dzcc08, controlMPC_ubIdx08, controlMPC_dlubcc08);
controlMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_2(controlMPC_ccrhsl09, controlMPC_slb09, controlMPC_llbbyslb09, controlMPC_dzcc09, controlMPC_lbIdx09, controlMPC_dllbcc09);
controlMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_2(controlMPC_ccrhsub09, controlMPC_sub09, controlMPC_lubbysub09, controlMPC_dzcc09, controlMPC_ubIdx09, controlMPC_dlubcc09);
controlMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_2(controlMPC_ccrhsl10, controlMPC_slb10, controlMPC_llbbyslb10, controlMPC_dzcc10, controlMPC_lbIdx10, controlMPC_dllbcc10);
controlMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_2(controlMPC_ccrhsub10, controlMPC_sub10, controlMPC_lubbysub10, controlMPC_dzcc10, controlMPC_ubIdx10, controlMPC_dlubcc10);
controlMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_2(controlMPC_ccrhsl11, controlMPC_slb11, controlMPC_llbbyslb11, controlMPC_dzcc11, controlMPC_lbIdx11, controlMPC_dllbcc11);
controlMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_2(controlMPC_ccrhsub11, controlMPC_sub11, controlMPC_lubbysub11, controlMPC_dzcc11, controlMPC_ubIdx11, controlMPC_dlubcc11);
controlMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_2(controlMPC_ccrhsl12, controlMPC_slb12, controlMPC_llbbyslb12, controlMPC_dzcc12, controlMPC_lbIdx12, controlMPC_dllbcc12);
controlMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_2(controlMPC_ccrhsub12, controlMPC_sub12, controlMPC_lubbysub12, controlMPC_dzcc12, controlMPC_ubIdx12, controlMPC_dlubcc12);
controlMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_2(controlMPC_ccrhsl13, controlMPC_slb13, controlMPC_llbbyslb13, controlMPC_dzcc13, controlMPC_lbIdx13, controlMPC_dllbcc13);
controlMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_2(controlMPC_ccrhsub13, controlMPC_sub13, controlMPC_lubbysub13, controlMPC_dzcc13, controlMPC_ubIdx13, controlMPC_dlubcc13);
controlMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_2(controlMPC_ccrhsl14, controlMPC_slb14, controlMPC_llbbyslb14, controlMPC_dzcc14, controlMPC_lbIdx14, controlMPC_dllbcc14);
controlMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_2(controlMPC_ccrhsub14, controlMPC_sub14, controlMPC_lubbysub14, controlMPC_dzcc14, controlMPC_ubIdx14, controlMPC_dlubcc14);
controlMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_2(controlMPC_ccrhsl15, controlMPC_slb15, controlMPC_llbbyslb15, controlMPC_dzcc15, controlMPC_lbIdx15, controlMPC_dllbcc15);
controlMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_2(controlMPC_ccrhsub15, controlMPC_sub15, controlMPC_lubbysub15, controlMPC_dzcc15, controlMPC_ubIdx15, controlMPC_dlubcc15);
controlMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_2(controlMPC_ccrhsl16, controlMPC_slb16, controlMPC_llbbyslb16, controlMPC_dzcc16, controlMPC_lbIdx16, controlMPC_dllbcc16);
controlMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_2(controlMPC_ccrhsub16, controlMPC_sub16, controlMPC_lubbysub16, controlMPC_dzcc16, controlMPC_ubIdx16, controlMPC_dlubcc16);
controlMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_2(controlMPC_ccrhsl17, controlMPC_slb17, controlMPC_llbbyslb17, controlMPC_dzcc17, controlMPC_lbIdx17, controlMPC_dllbcc17);
controlMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_2(controlMPC_ccrhsub17, controlMPC_sub17, controlMPC_lubbysub17, controlMPC_dzcc17, controlMPC_ubIdx17, controlMPC_dlubcc17);
controlMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_2(controlMPC_ccrhsl18, controlMPC_slb18, controlMPC_llbbyslb18, controlMPC_dzcc18, controlMPC_lbIdx18, controlMPC_dllbcc18);
controlMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_2(controlMPC_ccrhsub18, controlMPC_sub18, controlMPC_lubbysub18, controlMPC_dzcc18, controlMPC_ubIdx18, controlMPC_dlubcc18);
controlMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_2(controlMPC_ccrhsl19, controlMPC_slb19, controlMPC_llbbyslb19, controlMPC_dzcc19, controlMPC_lbIdx19, controlMPC_dllbcc19);
controlMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_2(controlMPC_ccrhsub19, controlMPC_sub19, controlMPC_lubbysub19, controlMPC_dzcc19, controlMPC_ubIdx19, controlMPC_dlubcc19);
controlMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_2(controlMPC_ccrhsl20, controlMPC_slb20, controlMPC_llbbyslb20, controlMPC_dzcc20, controlMPC_lbIdx20, controlMPC_dllbcc20);
controlMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_2(controlMPC_ccrhsub20, controlMPC_sub20, controlMPC_lubbysub20, controlMPC_dzcc20, controlMPC_ubIdx20, controlMPC_dlubcc20);
controlMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_2(controlMPC_ccrhsl21, controlMPC_slb21, controlMPC_llbbyslb21, controlMPC_dzcc21, controlMPC_lbIdx21, controlMPC_dllbcc21);
controlMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_2(controlMPC_ccrhsub21, controlMPC_sub21, controlMPC_lubbysub21, controlMPC_dzcc21, controlMPC_ubIdx21, controlMPC_dlubcc21);
controlMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_2(controlMPC_ccrhsl22, controlMPC_slb22, controlMPC_llbbyslb22, controlMPC_dzcc22, controlMPC_lbIdx22, controlMPC_dllbcc22);
controlMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_2(controlMPC_ccrhsub22, controlMPC_sub22, controlMPC_lubbysub22, controlMPC_dzcc22, controlMPC_ubIdx22, controlMPC_dlubcc22);
controlMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_2(controlMPC_ccrhsl23, controlMPC_slb23, controlMPC_llbbyslb23, controlMPC_dzcc23, controlMPC_lbIdx23, controlMPC_dllbcc23);
controlMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_2(controlMPC_ccrhsub23, controlMPC_sub23, controlMPC_lubbysub23, controlMPC_dzcc23, controlMPC_ubIdx23, controlMPC_dlubcc23);
controlMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_2(controlMPC_ccrhsl24, controlMPC_slb24, controlMPC_llbbyslb24, controlMPC_dzcc24, controlMPC_lbIdx24, controlMPC_dllbcc24);
controlMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_2(controlMPC_ccrhsub24, controlMPC_sub24, controlMPC_lubbysub24, controlMPC_dzcc24, controlMPC_ubIdx24, controlMPC_dlubcc24);
controlMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_2(controlMPC_ccrhsl25, controlMPC_slb25, controlMPC_llbbyslb25, controlMPC_dzcc25, controlMPC_lbIdx25, controlMPC_dllbcc25);
controlMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_2(controlMPC_ccrhsub25, controlMPC_sub25, controlMPC_lubbysub25, controlMPC_dzcc25, controlMPC_ubIdx25, controlMPC_dlubcc25);
controlMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_2(controlMPC_ccrhsl26, controlMPC_slb26, controlMPC_llbbyslb26, controlMPC_dzcc26, controlMPC_lbIdx26, controlMPC_dllbcc26);
controlMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_2(controlMPC_ccrhsub26, controlMPC_sub26, controlMPC_lubbysub26, controlMPC_dzcc26, controlMPC_ubIdx26, controlMPC_dlubcc26);
controlMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_2(controlMPC_ccrhsl27, controlMPC_slb27, controlMPC_llbbyslb27, controlMPC_dzcc27, controlMPC_lbIdx27, controlMPC_dllbcc27);
controlMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_2(controlMPC_ccrhsub27, controlMPC_sub27, controlMPC_lubbysub27, controlMPC_dzcc27, controlMPC_ubIdx27, controlMPC_dlubcc27);
controlMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_2(controlMPC_ccrhsl28, controlMPC_slb28, controlMPC_llbbyslb28, controlMPC_dzcc28, controlMPC_lbIdx28, controlMPC_dllbcc28);
controlMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_2(controlMPC_ccrhsub28, controlMPC_sub28, controlMPC_lubbysub28, controlMPC_dzcc28, controlMPC_ubIdx28, controlMPC_dlubcc28);
controlMPC_LA_VSUB7_116(controlMPC_l, controlMPC_ccrhs, controlMPC_s, controlMPC_dl_cc, controlMPC_ds_cc);
controlMPC_LA_VADD_58(controlMPC_dz_cc, controlMPC_dz_aff);
controlMPC_LA_VADD_56(controlMPC_dv_cc, controlMPC_dv_aff);
controlMPC_LA_VADD_116(controlMPC_dl_cc, controlMPC_dl_aff);
controlMPC_LA_VADD_116(controlMPC_ds_cc, controlMPC_ds_aff);
info->lsit_cc = controlMPC_LINESEARCH_BACKTRACKING_COMBINED(controlMPC_z, controlMPC_v, controlMPC_l, controlMPC_s, controlMPC_dz_cc, controlMPC_dv_cc, controlMPC_dl_cc, controlMPC_ds_cc, &info->step_cc, &info->mu);
if( info->lsit_cc == controlMPC_NOPROGRESS ){
exitcode = controlMPC_NOPROGRESS; break;
}
info->it++;
}
output->z1[0] = controlMPC_z00[0];
output->z1[1] = controlMPC_z00[1];
output->z2[0] = controlMPC_z01[0];
output->z2[1] = controlMPC_z01[1];
output->z3[0] = controlMPC_z02[0];
output->z3[1] = controlMPC_z02[1];
output->z4[0] = controlMPC_z03[0];
output->z4[1] = controlMPC_z03[1];
output->z5[0] = controlMPC_z04[0];
output->z5[1] = controlMPC_z04[1];
output->z6[0] = controlMPC_z05[0];
output->z6[1] = controlMPC_z05[1];
output->z7[0] = controlMPC_z06[0];
output->z7[1] = controlMPC_z06[1];
output->z8[0] = controlMPC_z07[0];
output->z8[1] = controlMPC_z07[1];
output->z9[0] = controlMPC_z08[0];
output->z9[1] = controlMPC_z08[1];
output->z10[0] = controlMPC_z09[0];
output->z10[1] = controlMPC_z09[1];
output->z11[0] = controlMPC_z10[0];
output->z11[1] = controlMPC_z10[1];
output->z12[0] = controlMPC_z11[0];
output->z12[1] = controlMPC_z11[1];
output->z13[0] = controlMPC_z12[0];
output->z13[1] = controlMPC_z12[1];
output->z14[0] = controlMPC_z13[0];
output->z14[1] = controlMPC_z13[1];
output->z15[0] = controlMPC_z14[0];
output->z15[1] = controlMPC_z14[1];
output->z16[0] = controlMPC_z15[0];
output->z16[1] = controlMPC_z15[1];
output->z17[0] = controlMPC_z16[0];
output->z17[1] = controlMPC_z16[1];
output->z18[0] = controlMPC_z17[0];
output->z18[1] = controlMPC_z17[1];
output->z19[0] = controlMPC_z18[0];
output->z19[1] = controlMPC_z18[1];
output->z20[0] = controlMPC_z19[0];
output->z20[1] = controlMPC_z19[1];
output->z21[0] = controlMPC_z20[0];
output->z21[1] = controlMPC_z20[1];
output->z22[0] = controlMPC_z21[0];
output->z22[1] = controlMPC_z21[1];
output->z23[0] = controlMPC_z22[0];
output->z23[1] = controlMPC_z22[1];
output->z24[0] = controlMPC_z23[0];
output->z24[1] = controlMPC_z23[1];
output->z25[0] = controlMPC_z24[0];
output->z25[1] = controlMPC_z24[1];
output->z26[0] = controlMPC_z25[0];
output->z26[1] = controlMPC_z25[1];
output->z27[0] = controlMPC_z26[0];
output->z27[1] = controlMPC_z26[1];
output->z28[0] = controlMPC_z27[0];
output->z28[1] = controlMPC_z27[1];
output->z29[0] = controlMPC_z28[0];
output->z29[1] = controlMPC_z28[1];

#if controlMPC_SET_TIMING == 1
info->solvetime = controlMPC_toc(&solvertimer);
#if controlMPC_SET_PRINTLEVEL > 0 && controlMPC_SET_TIMING == 1
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
