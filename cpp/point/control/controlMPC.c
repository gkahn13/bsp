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

#include "controlMPC.h"

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
 * Initializes a vector of length 18 with a value.
 */
void controlMPC_LA_INITIALIZEVECTOR_18(controlMPC_FLOAT* vec, controlMPC_FLOAT value)
{
	int i;
	for( i=0; i<18; i++ )
	{
		vec[i] = value;
	}
}


/*
 * Initializes a vector of length 16 with a value.
 */
void controlMPC_LA_INITIALIZEVECTOR_16(controlMPC_FLOAT* vec, controlMPC_FLOAT value)
{
	int i;
	for( i=0; i<16; i++ )
	{
		vec[i] = value;
	}
}


/*
 * Initializes a vector of length 36 with a value.
 */
void controlMPC_LA_INITIALIZEVECTOR_36(controlMPC_FLOAT* vec, controlMPC_FLOAT value)
{
	int i;
	for( i=0; i<36; i++ )
	{
		vec[i] = value;
	}
}


/* 
 * Calculates a dot product and adds it to a variable: z += x'*y; 
 * This function is for vectors of length 36.
 */
void controlMPC_LA_DOTACC_36(controlMPC_FLOAT *x, controlMPC_FLOAT *y, controlMPC_FLOAT *z)
{
	int i;
	for( i=0; i<36; i++ ){
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
 * of length 18.
 */
void controlMPC_LA_VVADD3_18(controlMPC_FLOAT *u, controlMPC_FLOAT *v, controlMPC_FLOAT *w, controlMPC_FLOAT *z)
{
	int i;
	for( i=0; i<18; i++ ){
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
 * Vector subtraction z = -x - y for vectors of length 18.
 */
void controlMPC_LA_VSUB2_18(controlMPC_FLOAT *x, controlMPC_FLOAT *y, controlMPC_FLOAT *z)
{
	int i;
	for( i=0; i<18; i++){
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
        for( i=0; i<36; i++ ){
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
        if( i == 36 ){
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
    *mu_aff = mymu / (controlMPC_FLOAT)36;
    return lsIt;
}


/*
 * Vector subtraction x = u.*v - a where a is a scalar
*  and x,u,v are vectors of length 36.
 */
void controlMPC_LA_VSUB5_36(controlMPC_FLOAT *u, controlMPC_FLOAT *v, controlMPC_FLOAT a, controlMPC_FLOAT *x)
{
	int i;
	for( i=0; i<36; i++){
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
 * Vector subtraction z = x - y for vectors of length 18.
 */
void controlMPC_LA_VSUB_18(controlMPC_FLOAT *x, controlMPC_FLOAT *y, controlMPC_FLOAT *z)
{
	int i;
	for( i=0; i<18; i++){
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
 * Computes ds = -l.\(r + s.*dl) for vectors of length 36.
 */
void controlMPC_LA_VSUB7_36(controlMPC_FLOAT *l, controlMPC_FLOAT *r, controlMPC_FLOAT *s, controlMPC_FLOAT *dl, controlMPC_FLOAT *ds)
{
	int i;
	for( i=0; i<36; i++){
		ds[i] = -(r[i] + s[i]*dl[i])/l[i];
	}
}


/*
 * Vector addition x = x + y for vectors of length 18.
 */
void controlMPC_LA_VADD_18(controlMPC_FLOAT *x, controlMPC_FLOAT *y)
{
	int i;
	for( i=0; i<18; i++){
		x[i] += y[i];
	}
}


/*
 * Vector addition x = x + y for vectors of length 16.
 */
void controlMPC_LA_VADD_16(controlMPC_FLOAT *x, controlMPC_FLOAT *y)
{
	int i;
	for( i=0; i<16; i++){
		x[i] += y[i];
	}
}


/*
 * Vector addition x = x + y for vectors of length 36.
 */
void controlMPC_LA_VADD_36(controlMPC_FLOAT *x, controlMPC_FLOAT *y)
{
	int i;
	for( i=0; i<36; i++){
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
        for( i=0; i<36; i++ ){
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
        if( i == 36 ){
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
    for( i=0; i<18; i++ ){
        z[i] += a_gamma*dz[i];
    }
    
    /* equality constraint multipliers */
    for( i=0; i<16; i++ ){
        v[i] += a_gamma*dv[i];
    }
    
    /* inequality constraint multipliers & slacks, also update mu */
    *mu = 0;
    for( i=0; i<36; i++ ){
        dltemp = l[i] + a_gamma*dl[i]; l[i] = dltemp;
        dstemp = s[i] + a_gamma*ds[i]; s[i] = dstemp;
        *mu += dltemp*dstemp;
    }
    
    *a = a_gamma;
    *mu /= (controlMPC_FLOAT)36;
    return lsIt;
}




/* VARIABLE DEFINITIONS ------------------------------------------------ */
controlMPC_FLOAT controlMPC_z[18];
controlMPC_FLOAT controlMPC_v[16];
controlMPC_FLOAT controlMPC_dz_aff[18];
controlMPC_FLOAT controlMPC_dv_aff[16];
controlMPC_FLOAT controlMPC_grad_cost[18];
controlMPC_FLOAT controlMPC_grad_eq[18];
controlMPC_FLOAT controlMPC_rd[18];
controlMPC_FLOAT controlMPC_l[36];
controlMPC_FLOAT controlMPC_s[36];
controlMPC_FLOAT controlMPC_lbys[36];
controlMPC_FLOAT controlMPC_dl_aff[36];
controlMPC_FLOAT controlMPC_ds_aff[36];
controlMPC_FLOAT controlMPC_dz_cc[18];
controlMPC_FLOAT controlMPC_dv_cc[16];
controlMPC_FLOAT controlMPC_dl_cc[36];
controlMPC_FLOAT controlMPC_ds_cc[36];
controlMPC_FLOAT controlMPC_ccrhs[36];
controlMPC_FLOAT controlMPC_grad_ineq[18];
controlMPC_FLOAT* controlMPC_z0 = controlMPC_z + 0;
controlMPC_FLOAT* controlMPC_dzaff0 = controlMPC_dz_aff + 0;
controlMPC_FLOAT* controlMPC_dzcc0 = controlMPC_dz_cc + 0;
controlMPC_FLOAT* controlMPC_rd0 = controlMPC_rd + 0;
controlMPC_FLOAT controlMPC_Lbyrd0[2];
controlMPC_FLOAT* controlMPC_grad_cost0 = controlMPC_grad_cost + 0;
controlMPC_FLOAT* controlMPC_grad_eq0 = controlMPC_grad_eq + 0;
controlMPC_FLOAT* controlMPC_grad_ineq0 = controlMPC_grad_ineq + 0;
controlMPC_FLOAT controlMPC_ctv0[2];
controlMPC_FLOAT controlMPC_C0[4] = {0.0000000000000000E+000, 0.0000000000000000E+000, 
0.0000000000000000E+000, 0.0000000000000000E+000};
controlMPC_FLOAT* controlMPC_v0 = controlMPC_v + 0;
controlMPC_FLOAT controlMPC_re0[2];
controlMPC_FLOAT controlMPC_beta0[2];
controlMPC_FLOAT controlMPC_betacc0[2];
controlMPC_FLOAT* controlMPC_dvaff0 = controlMPC_dv_aff + 0;
controlMPC_FLOAT* controlMPC_dvcc0 = controlMPC_dv_cc + 0;
controlMPC_FLOAT controlMPC_V0[4];
controlMPC_FLOAT controlMPC_Yd0[3];
controlMPC_FLOAT controlMPC_Ld0[3];
controlMPC_FLOAT controlMPC_yy0[2];
controlMPC_FLOAT controlMPC_bmy0[2];
controlMPC_FLOAT controlMPC_c0[2] = {0.0000000000000000E+000, 0.0000000000000000E+000};
int controlMPC_lbIdx0[2] = {0, 1};
controlMPC_FLOAT* controlMPC_llb0 = controlMPC_l + 0;
controlMPC_FLOAT* controlMPC_slb0 = controlMPC_s + 0;
controlMPC_FLOAT* controlMPC_llbbyslb0 = controlMPC_lbys + 0;
controlMPC_FLOAT controlMPC_rilb0[2];
controlMPC_FLOAT* controlMPC_dllbaff0 = controlMPC_dl_aff + 0;
controlMPC_FLOAT* controlMPC_dslbaff0 = controlMPC_ds_aff + 0;
controlMPC_FLOAT* controlMPC_dllbcc0 = controlMPC_dl_cc + 0;
controlMPC_FLOAT* controlMPC_dslbcc0 = controlMPC_ds_cc + 0;
controlMPC_FLOAT* controlMPC_ccrhsl0 = controlMPC_ccrhs + 0;
int controlMPC_ubIdx0[2] = {0, 1};
controlMPC_FLOAT* controlMPC_lub0 = controlMPC_l + 2;
controlMPC_FLOAT* controlMPC_sub0 = controlMPC_s + 2;
controlMPC_FLOAT* controlMPC_lubbysub0 = controlMPC_lbys + 2;
controlMPC_FLOAT controlMPC_riub0[2];
controlMPC_FLOAT* controlMPC_dlubaff0 = controlMPC_dl_aff + 2;
controlMPC_FLOAT* controlMPC_dsubaff0 = controlMPC_ds_aff + 2;
controlMPC_FLOAT* controlMPC_dlubcc0 = controlMPC_dl_cc + 2;
controlMPC_FLOAT* controlMPC_dsubcc0 = controlMPC_ds_cc + 2;
controlMPC_FLOAT* controlMPC_ccrhsub0 = controlMPC_ccrhs + 2;
controlMPC_FLOAT controlMPC_Phi0[2];
controlMPC_FLOAT* controlMPC_z1 = controlMPC_z + 2;
controlMPC_FLOAT* controlMPC_dzaff1 = controlMPC_dz_aff + 2;
controlMPC_FLOAT* controlMPC_dzcc1 = controlMPC_dz_cc + 2;
controlMPC_FLOAT* controlMPC_rd1 = controlMPC_rd + 2;
controlMPC_FLOAT controlMPC_Lbyrd1[2];
controlMPC_FLOAT* controlMPC_grad_cost1 = controlMPC_grad_cost + 2;
controlMPC_FLOAT* controlMPC_grad_eq1 = controlMPC_grad_eq + 2;
controlMPC_FLOAT* controlMPC_grad_ineq1 = controlMPC_grad_ineq + 2;
controlMPC_FLOAT controlMPC_ctv1[2];
controlMPC_FLOAT* controlMPC_v1 = controlMPC_v + 2;
controlMPC_FLOAT controlMPC_re1[2];
controlMPC_FLOAT controlMPC_beta1[2];
controlMPC_FLOAT controlMPC_betacc1[2];
controlMPC_FLOAT* controlMPC_dvaff1 = controlMPC_dv_aff + 2;
controlMPC_FLOAT* controlMPC_dvcc1 = controlMPC_dv_cc + 2;
controlMPC_FLOAT controlMPC_V1[4];
controlMPC_FLOAT controlMPC_Yd1[3];
controlMPC_FLOAT controlMPC_Ld1[3];
controlMPC_FLOAT controlMPC_yy1[2];
controlMPC_FLOAT controlMPC_bmy1[2];
controlMPC_FLOAT controlMPC_c1[2] = {0.0000000000000000E+000, 0.0000000000000000E+000};
int controlMPC_lbIdx1[2] = {0, 1};
controlMPC_FLOAT* controlMPC_llb1 = controlMPC_l + 4;
controlMPC_FLOAT* controlMPC_slb1 = controlMPC_s + 4;
controlMPC_FLOAT* controlMPC_llbbyslb1 = controlMPC_lbys + 4;
controlMPC_FLOAT controlMPC_rilb1[2];
controlMPC_FLOAT* controlMPC_dllbaff1 = controlMPC_dl_aff + 4;
controlMPC_FLOAT* controlMPC_dslbaff1 = controlMPC_ds_aff + 4;
controlMPC_FLOAT* controlMPC_dllbcc1 = controlMPC_dl_cc + 4;
controlMPC_FLOAT* controlMPC_dslbcc1 = controlMPC_ds_cc + 4;
controlMPC_FLOAT* controlMPC_ccrhsl1 = controlMPC_ccrhs + 4;
int controlMPC_ubIdx1[2] = {0, 1};
controlMPC_FLOAT* controlMPC_lub1 = controlMPC_l + 6;
controlMPC_FLOAT* controlMPC_sub1 = controlMPC_s + 6;
controlMPC_FLOAT* controlMPC_lubbysub1 = controlMPC_lbys + 6;
controlMPC_FLOAT controlMPC_riub1[2];
controlMPC_FLOAT* controlMPC_dlubaff1 = controlMPC_dl_aff + 6;
controlMPC_FLOAT* controlMPC_dsubaff1 = controlMPC_ds_aff + 6;
controlMPC_FLOAT* controlMPC_dlubcc1 = controlMPC_dl_cc + 6;
controlMPC_FLOAT* controlMPC_dsubcc1 = controlMPC_ds_cc + 6;
controlMPC_FLOAT* controlMPC_ccrhsub1 = controlMPC_ccrhs + 6;
controlMPC_FLOAT controlMPC_Phi1[2];
controlMPC_FLOAT controlMPC_D1[2] = {0.0000000000000000E+000, 
0.0000000000000000E+000};
controlMPC_FLOAT controlMPC_W1[2];
controlMPC_FLOAT controlMPC_Ysd1[4];
controlMPC_FLOAT controlMPC_Lsd1[4];
controlMPC_FLOAT* controlMPC_z2 = controlMPC_z + 4;
controlMPC_FLOAT* controlMPC_dzaff2 = controlMPC_dz_aff + 4;
controlMPC_FLOAT* controlMPC_dzcc2 = controlMPC_dz_cc + 4;
controlMPC_FLOAT* controlMPC_rd2 = controlMPC_rd + 4;
controlMPC_FLOAT controlMPC_Lbyrd2[2];
controlMPC_FLOAT* controlMPC_grad_cost2 = controlMPC_grad_cost + 4;
controlMPC_FLOAT* controlMPC_grad_eq2 = controlMPC_grad_eq + 4;
controlMPC_FLOAT* controlMPC_grad_ineq2 = controlMPC_grad_ineq + 4;
controlMPC_FLOAT controlMPC_ctv2[2];
controlMPC_FLOAT* controlMPC_v2 = controlMPC_v + 4;
controlMPC_FLOAT controlMPC_re2[2];
controlMPC_FLOAT controlMPC_beta2[2];
controlMPC_FLOAT controlMPC_betacc2[2];
controlMPC_FLOAT* controlMPC_dvaff2 = controlMPC_dv_aff + 4;
controlMPC_FLOAT* controlMPC_dvcc2 = controlMPC_dv_cc + 4;
controlMPC_FLOAT controlMPC_V2[4];
controlMPC_FLOAT controlMPC_Yd2[3];
controlMPC_FLOAT controlMPC_Ld2[3];
controlMPC_FLOAT controlMPC_yy2[2];
controlMPC_FLOAT controlMPC_bmy2[2];
controlMPC_FLOAT controlMPC_c2[2] = {0.0000000000000000E+000, 0.0000000000000000E+000};
int controlMPC_lbIdx2[2] = {0, 1};
controlMPC_FLOAT* controlMPC_llb2 = controlMPC_l + 8;
controlMPC_FLOAT* controlMPC_slb2 = controlMPC_s + 8;
controlMPC_FLOAT* controlMPC_llbbyslb2 = controlMPC_lbys + 8;
controlMPC_FLOAT controlMPC_rilb2[2];
controlMPC_FLOAT* controlMPC_dllbaff2 = controlMPC_dl_aff + 8;
controlMPC_FLOAT* controlMPC_dslbaff2 = controlMPC_ds_aff + 8;
controlMPC_FLOAT* controlMPC_dllbcc2 = controlMPC_dl_cc + 8;
controlMPC_FLOAT* controlMPC_dslbcc2 = controlMPC_ds_cc + 8;
controlMPC_FLOAT* controlMPC_ccrhsl2 = controlMPC_ccrhs + 8;
int controlMPC_ubIdx2[2] = {0, 1};
controlMPC_FLOAT* controlMPC_lub2 = controlMPC_l + 10;
controlMPC_FLOAT* controlMPC_sub2 = controlMPC_s + 10;
controlMPC_FLOAT* controlMPC_lubbysub2 = controlMPC_lbys + 10;
controlMPC_FLOAT controlMPC_riub2[2];
controlMPC_FLOAT* controlMPC_dlubaff2 = controlMPC_dl_aff + 10;
controlMPC_FLOAT* controlMPC_dsubaff2 = controlMPC_ds_aff + 10;
controlMPC_FLOAT* controlMPC_dlubcc2 = controlMPC_dl_cc + 10;
controlMPC_FLOAT* controlMPC_dsubcc2 = controlMPC_ds_cc + 10;
controlMPC_FLOAT* controlMPC_ccrhsub2 = controlMPC_ccrhs + 10;
controlMPC_FLOAT controlMPC_Phi2[2];
controlMPC_FLOAT controlMPC_W2[2];
controlMPC_FLOAT controlMPC_Ysd2[4];
controlMPC_FLOAT controlMPC_Lsd2[4];
controlMPC_FLOAT* controlMPC_z3 = controlMPC_z + 6;
controlMPC_FLOAT* controlMPC_dzaff3 = controlMPC_dz_aff + 6;
controlMPC_FLOAT* controlMPC_dzcc3 = controlMPC_dz_cc + 6;
controlMPC_FLOAT* controlMPC_rd3 = controlMPC_rd + 6;
controlMPC_FLOAT controlMPC_Lbyrd3[2];
controlMPC_FLOAT* controlMPC_grad_cost3 = controlMPC_grad_cost + 6;
controlMPC_FLOAT* controlMPC_grad_eq3 = controlMPC_grad_eq + 6;
controlMPC_FLOAT* controlMPC_grad_ineq3 = controlMPC_grad_ineq + 6;
controlMPC_FLOAT controlMPC_ctv3[2];
controlMPC_FLOAT* controlMPC_v3 = controlMPC_v + 6;
controlMPC_FLOAT controlMPC_re3[2];
controlMPC_FLOAT controlMPC_beta3[2];
controlMPC_FLOAT controlMPC_betacc3[2];
controlMPC_FLOAT* controlMPC_dvaff3 = controlMPC_dv_aff + 6;
controlMPC_FLOAT* controlMPC_dvcc3 = controlMPC_dv_cc + 6;
controlMPC_FLOAT controlMPC_V3[4];
controlMPC_FLOAT controlMPC_Yd3[3];
controlMPC_FLOAT controlMPC_Ld3[3];
controlMPC_FLOAT controlMPC_yy3[2];
controlMPC_FLOAT controlMPC_bmy3[2];
controlMPC_FLOAT controlMPC_c3[2] = {0.0000000000000000E+000, 0.0000000000000000E+000};
int controlMPC_lbIdx3[2] = {0, 1};
controlMPC_FLOAT* controlMPC_llb3 = controlMPC_l + 12;
controlMPC_FLOAT* controlMPC_slb3 = controlMPC_s + 12;
controlMPC_FLOAT* controlMPC_llbbyslb3 = controlMPC_lbys + 12;
controlMPC_FLOAT controlMPC_rilb3[2];
controlMPC_FLOAT* controlMPC_dllbaff3 = controlMPC_dl_aff + 12;
controlMPC_FLOAT* controlMPC_dslbaff3 = controlMPC_ds_aff + 12;
controlMPC_FLOAT* controlMPC_dllbcc3 = controlMPC_dl_cc + 12;
controlMPC_FLOAT* controlMPC_dslbcc3 = controlMPC_ds_cc + 12;
controlMPC_FLOAT* controlMPC_ccrhsl3 = controlMPC_ccrhs + 12;
int controlMPC_ubIdx3[2] = {0, 1};
controlMPC_FLOAT* controlMPC_lub3 = controlMPC_l + 14;
controlMPC_FLOAT* controlMPC_sub3 = controlMPC_s + 14;
controlMPC_FLOAT* controlMPC_lubbysub3 = controlMPC_lbys + 14;
controlMPC_FLOAT controlMPC_riub3[2];
controlMPC_FLOAT* controlMPC_dlubaff3 = controlMPC_dl_aff + 14;
controlMPC_FLOAT* controlMPC_dsubaff3 = controlMPC_ds_aff + 14;
controlMPC_FLOAT* controlMPC_dlubcc3 = controlMPC_dl_cc + 14;
controlMPC_FLOAT* controlMPC_dsubcc3 = controlMPC_ds_cc + 14;
controlMPC_FLOAT* controlMPC_ccrhsub3 = controlMPC_ccrhs + 14;
controlMPC_FLOAT controlMPC_Phi3[2];
controlMPC_FLOAT controlMPC_W3[2];
controlMPC_FLOAT controlMPC_Ysd3[4];
controlMPC_FLOAT controlMPC_Lsd3[4];
controlMPC_FLOAT* controlMPC_z4 = controlMPC_z + 8;
controlMPC_FLOAT* controlMPC_dzaff4 = controlMPC_dz_aff + 8;
controlMPC_FLOAT* controlMPC_dzcc4 = controlMPC_dz_cc + 8;
controlMPC_FLOAT* controlMPC_rd4 = controlMPC_rd + 8;
controlMPC_FLOAT controlMPC_Lbyrd4[2];
controlMPC_FLOAT* controlMPC_grad_cost4 = controlMPC_grad_cost + 8;
controlMPC_FLOAT* controlMPC_grad_eq4 = controlMPC_grad_eq + 8;
controlMPC_FLOAT* controlMPC_grad_ineq4 = controlMPC_grad_ineq + 8;
controlMPC_FLOAT controlMPC_ctv4[2];
controlMPC_FLOAT* controlMPC_v4 = controlMPC_v + 8;
controlMPC_FLOAT controlMPC_re4[2];
controlMPC_FLOAT controlMPC_beta4[2];
controlMPC_FLOAT controlMPC_betacc4[2];
controlMPC_FLOAT* controlMPC_dvaff4 = controlMPC_dv_aff + 8;
controlMPC_FLOAT* controlMPC_dvcc4 = controlMPC_dv_cc + 8;
controlMPC_FLOAT controlMPC_V4[4];
controlMPC_FLOAT controlMPC_Yd4[3];
controlMPC_FLOAT controlMPC_Ld4[3];
controlMPC_FLOAT controlMPC_yy4[2];
controlMPC_FLOAT controlMPC_bmy4[2];
controlMPC_FLOAT controlMPC_c4[2] = {0.0000000000000000E+000, 0.0000000000000000E+000};
int controlMPC_lbIdx4[2] = {0, 1};
controlMPC_FLOAT* controlMPC_llb4 = controlMPC_l + 16;
controlMPC_FLOAT* controlMPC_slb4 = controlMPC_s + 16;
controlMPC_FLOAT* controlMPC_llbbyslb4 = controlMPC_lbys + 16;
controlMPC_FLOAT controlMPC_rilb4[2];
controlMPC_FLOAT* controlMPC_dllbaff4 = controlMPC_dl_aff + 16;
controlMPC_FLOAT* controlMPC_dslbaff4 = controlMPC_ds_aff + 16;
controlMPC_FLOAT* controlMPC_dllbcc4 = controlMPC_dl_cc + 16;
controlMPC_FLOAT* controlMPC_dslbcc4 = controlMPC_ds_cc + 16;
controlMPC_FLOAT* controlMPC_ccrhsl4 = controlMPC_ccrhs + 16;
int controlMPC_ubIdx4[2] = {0, 1};
controlMPC_FLOAT* controlMPC_lub4 = controlMPC_l + 18;
controlMPC_FLOAT* controlMPC_sub4 = controlMPC_s + 18;
controlMPC_FLOAT* controlMPC_lubbysub4 = controlMPC_lbys + 18;
controlMPC_FLOAT controlMPC_riub4[2];
controlMPC_FLOAT* controlMPC_dlubaff4 = controlMPC_dl_aff + 18;
controlMPC_FLOAT* controlMPC_dsubaff4 = controlMPC_ds_aff + 18;
controlMPC_FLOAT* controlMPC_dlubcc4 = controlMPC_dl_cc + 18;
controlMPC_FLOAT* controlMPC_dsubcc4 = controlMPC_ds_cc + 18;
controlMPC_FLOAT* controlMPC_ccrhsub4 = controlMPC_ccrhs + 18;
controlMPC_FLOAT controlMPC_Phi4[2];
controlMPC_FLOAT controlMPC_W4[2];
controlMPC_FLOAT controlMPC_Ysd4[4];
controlMPC_FLOAT controlMPC_Lsd4[4];
controlMPC_FLOAT* controlMPC_z5 = controlMPC_z + 10;
controlMPC_FLOAT* controlMPC_dzaff5 = controlMPC_dz_aff + 10;
controlMPC_FLOAT* controlMPC_dzcc5 = controlMPC_dz_cc + 10;
controlMPC_FLOAT* controlMPC_rd5 = controlMPC_rd + 10;
controlMPC_FLOAT controlMPC_Lbyrd5[2];
controlMPC_FLOAT* controlMPC_grad_cost5 = controlMPC_grad_cost + 10;
controlMPC_FLOAT* controlMPC_grad_eq5 = controlMPC_grad_eq + 10;
controlMPC_FLOAT* controlMPC_grad_ineq5 = controlMPC_grad_ineq + 10;
controlMPC_FLOAT controlMPC_ctv5[2];
controlMPC_FLOAT* controlMPC_v5 = controlMPC_v + 10;
controlMPC_FLOAT controlMPC_re5[2];
controlMPC_FLOAT controlMPC_beta5[2];
controlMPC_FLOAT controlMPC_betacc5[2];
controlMPC_FLOAT* controlMPC_dvaff5 = controlMPC_dv_aff + 10;
controlMPC_FLOAT* controlMPC_dvcc5 = controlMPC_dv_cc + 10;
controlMPC_FLOAT controlMPC_V5[4];
controlMPC_FLOAT controlMPC_Yd5[3];
controlMPC_FLOAT controlMPC_Ld5[3];
controlMPC_FLOAT controlMPC_yy5[2];
controlMPC_FLOAT controlMPC_bmy5[2];
controlMPC_FLOAT controlMPC_c5[2] = {0.0000000000000000E+000, 0.0000000000000000E+000};
int controlMPC_lbIdx5[2] = {0, 1};
controlMPC_FLOAT* controlMPC_llb5 = controlMPC_l + 20;
controlMPC_FLOAT* controlMPC_slb5 = controlMPC_s + 20;
controlMPC_FLOAT* controlMPC_llbbyslb5 = controlMPC_lbys + 20;
controlMPC_FLOAT controlMPC_rilb5[2];
controlMPC_FLOAT* controlMPC_dllbaff5 = controlMPC_dl_aff + 20;
controlMPC_FLOAT* controlMPC_dslbaff5 = controlMPC_ds_aff + 20;
controlMPC_FLOAT* controlMPC_dllbcc5 = controlMPC_dl_cc + 20;
controlMPC_FLOAT* controlMPC_dslbcc5 = controlMPC_ds_cc + 20;
controlMPC_FLOAT* controlMPC_ccrhsl5 = controlMPC_ccrhs + 20;
int controlMPC_ubIdx5[2] = {0, 1};
controlMPC_FLOAT* controlMPC_lub5 = controlMPC_l + 22;
controlMPC_FLOAT* controlMPC_sub5 = controlMPC_s + 22;
controlMPC_FLOAT* controlMPC_lubbysub5 = controlMPC_lbys + 22;
controlMPC_FLOAT controlMPC_riub5[2];
controlMPC_FLOAT* controlMPC_dlubaff5 = controlMPC_dl_aff + 22;
controlMPC_FLOAT* controlMPC_dsubaff5 = controlMPC_ds_aff + 22;
controlMPC_FLOAT* controlMPC_dlubcc5 = controlMPC_dl_cc + 22;
controlMPC_FLOAT* controlMPC_dsubcc5 = controlMPC_ds_cc + 22;
controlMPC_FLOAT* controlMPC_ccrhsub5 = controlMPC_ccrhs + 22;
controlMPC_FLOAT controlMPC_Phi5[2];
controlMPC_FLOAT controlMPC_W5[2];
controlMPC_FLOAT controlMPC_Ysd5[4];
controlMPC_FLOAT controlMPC_Lsd5[4];
controlMPC_FLOAT* controlMPC_z6 = controlMPC_z + 12;
controlMPC_FLOAT* controlMPC_dzaff6 = controlMPC_dz_aff + 12;
controlMPC_FLOAT* controlMPC_dzcc6 = controlMPC_dz_cc + 12;
controlMPC_FLOAT* controlMPC_rd6 = controlMPC_rd + 12;
controlMPC_FLOAT controlMPC_Lbyrd6[2];
controlMPC_FLOAT* controlMPC_grad_cost6 = controlMPC_grad_cost + 12;
controlMPC_FLOAT* controlMPC_grad_eq6 = controlMPC_grad_eq + 12;
controlMPC_FLOAT* controlMPC_grad_ineq6 = controlMPC_grad_ineq + 12;
controlMPC_FLOAT controlMPC_ctv6[2];
controlMPC_FLOAT* controlMPC_v6 = controlMPC_v + 12;
controlMPC_FLOAT controlMPC_re6[2];
controlMPC_FLOAT controlMPC_beta6[2];
controlMPC_FLOAT controlMPC_betacc6[2];
controlMPC_FLOAT* controlMPC_dvaff6 = controlMPC_dv_aff + 12;
controlMPC_FLOAT* controlMPC_dvcc6 = controlMPC_dv_cc + 12;
controlMPC_FLOAT controlMPC_V6[4];
controlMPC_FLOAT controlMPC_Yd6[3];
controlMPC_FLOAT controlMPC_Ld6[3];
controlMPC_FLOAT controlMPC_yy6[2];
controlMPC_FLOAT controlMPC_bmy6[2];
controlMPC_FLOAT controlMPC_c6[2] = {0.0000000000000000E+000, 0.0000000000000000E+000};
int controlMPC_lbIdx6[2] = {0, 1};
controlMPC_FLOAT* controlMPC_llb6 = controlMPC_l + 24;
controlMPC_FLOAT* controlMPC_slb6 = controlMPC_s + 24;
controlMPC_FLOAT* controlMPC_llbbyslb6 = controlMPC_lbys + 24;
controlMPC_FLOAT controlMPC_rilb6[2];
controlMPC_FLOAT* controlMPC_dllbaff6 = controlMPC_dl_aff + 24;
controlMPC_FLOAT* controlMPC_dslbaff6 = controlMPC_ds_aff + 24;
controlMPC_FLOAT* controlMPC_dllbcc6 = controlMPC_dl_cc + 24;
controlMPC_FLOAT* controlMPC_dslbcc6 = controlMPC_ds_cc + 24;
controlMPC_FLOAT* controlMPC_ccrhsl6 = controlMPC_ccrhs + 24;
int controlMPC_ubIdx6[2] = {0, 1};
controlMPC_FLOAT* controlMPC_lub6 = controlMPC_l + 26;
controlMPC_FLOAT* controlMPC_sub6 = controlMPC_s + 26;
controlMPC_FLOAT* controlMPC_lubbysub6 = controlMPC_lbys + 26;
controlMPC_FLOAT controlMPC_riub6[2];
controlMPC_FLOAT* controlMPC_dlubaff6 = controlMPC_dl_aff + 26;
controlMPC_FLOAT* controlMPC_dsubaff6 = controlMPC_ds_aff + 26;
controlMPC_FLOAT* controlMPC_dlubcc6 = controlMPC_dl_cc + 26;
controlMPC_FLOAT* controlMPC_dsubcc6 = controlMPC_ds_cc + 26;
controlMPC_FLOAT* controlMPC_ccrhsub6 = controlMPC_ccrhs + 26;
controlMPC_FLOAT controlMPC_Phi6[2];
controlMPC_FLOAT controlMPC_W6[2];
controlMPC_FLOAT controlMPC_Ysd6[4];
controlMPC_FLOAT controlMPC_Lsd6[4];
controlMPC_FLOAT* controlMPC_z7 = controlMPC_z + 14;
controlMPC_FLOAT* controlMPC_dzaff7 = controlMPC_dz_aff + 14;
controlMPC_FLOAT* controlMPC_dzcc7 = controlMPC_dz_cc + 14;
controlMPC_FLOAT* controlMPC_rd7 = controlMPC_rd + 14;
controlMPC_FLOAT controlMPC_Lbyrd7[2];
controlMPC_FLOAT* controlMPC_grad_cost7 = controlMPC_grad_cost + 14;
controlMPC_FLOAT* controlMPC_grad_eq7 = controlMPC_grad_eq + 14;
controlMPC_FLOAT* controlMPC_grad_ineq7 = controlMPC_grad_ineq + 14;
controlMPC_FLOAT controlMPC_ctv7[2];
controlMPC_FLOAT* controlMPC_v7 = controlMPC_v + 14;
controlMPC_FLOAT controlMPC_re7[2];
controlMPC_FLOAT controlMPC_beta7[2];
controlMPC_FLOAT controlMPC_betacc7[2];
controlMPC_FLOAT* controlMPC_dvaff7 = controlMPC_dv_aff + 14;
controlMPC_FLOAT* controlMPC_dvcc7 = controlMPC_dv_cc + 14;
controlMPC_FLOAT controlMPC_V7[4];
controlMPC_FLOAT controlMPC_Yd7[3];
controlMPC_FLOAT controlMPC_Ld7[3];
controlMPC_FLOAT controlMPC_yy7[2];
controlMPC_FLOAT controlMPC_bmy7[2];
controlMPC_FLOAT controlMPC_c7[2] = {0.0000000000000000E+000, 0.0000000000000000E+000};
int controlMPC_lbIdx7[2] = {0, 1};
controlMPC_FLOAT* controlMPC_llb7 = controlMPC_l + 28;
controlMPC_FLOAT* controlMPC_slb7 = controlMPC_s + 28;
controlMPC_FLOAT* controlMPC_llbbyslb7 = controlMPC_lbys + 28;
controlMPC_FLOAT controlMPC_rilb7[2];
controlMPC_FLOAT* controlMPC_dllbaff7 = controlMPC_dl_aff + 28;
controlMPC_FLOAT* controlMPC_dslbaff7 = controlMPC_ds_aff + 28;
controlMPC_FLOAT* controlMPC_dllbcc7 = controlMPC_dl_cc + 28;
controlMPC_FLOAT* controlMPC_dslbcc7 = controlMPC_ds_cc + 28;
controlMPC_FLOAT* controlMPC_ccrhsl7 = controlMPC_ccrhs + 28;
int controlMPC_ubIdx7[2] = {0, 1};
controlMPC_FLOAT* controlMPC_lub7 = controlMPC_l + 30;
controlMPC_FLOAT* controlMPC_sub7 = controlMPC_s + 30;
controlMPC_FLOAT* controlMPC_lubbysub7 = controlMPC_lbys + 30;
controlMPC_FLOAT controlMPC_riub7[2];
controlMPC_FLOAT* controlMPC_dlubaff7 = controlMPC_dl_aff + 30;
controlMPC_FLOAT* controlMPC_dsubaff7 = controlMPC_ds_aff + 30;
controlMPC_FLOAT* controlMPC_dlubcc7 = controlMPC_dl_cc + 30;
controlMPC_FLOAT* controlMPC_dsubcc7 = controlMPC_ds_cc + 30;
controlMPC_FLOAT* controlMPC_ccrhsub7 = controlMPC_ccrhs + 30;
controlMPC_FLOAT controlMPC_Phi7[2];
controlMPC_FLOAT controlMPC_W7[2];
controlMPC_FLOAT controlMPC_Ysd7[4];
controlMPC_FLOAT controlMPC_Lsd7[4];
controlMPC_FLOAT* controlMPC_z8 = controlMPC_z + 16;
controlMPC_FLOAT* controlMPC_dzaff8 = controlMPC_dz_aff + 16;
controlMPC_FLOAT* controlMPC_dzcc8 = controlMPC_dz_cc + 16;
controlMPC_FLOAT* controlMPC_rd8 = controlMPC_rd + 16;
controlMPC_FLOAT controlMPC_Lbyrd8[2];
controlMPC_FLOAT* controlMPC_grad_cost8 = controlMPC_grad_cost + 16;
controlMPC_FLOAT* controlMPC_grad_eq8 = controlMPC_grad_eq + 16;
controlMPC_FLOAT* controlMPC_grad_ineq8 = controlMPC_grad_ineq + 16;
controlMPC_FLOAT controlMPC_ctv8[2];
int controlMPC_lbIdx8[2] = {0, 1};
controlMPC_FLOAT* controlMPC_llb8 = controlMPC_l + 32;
controlMPC_FLOAT* controlMPC_slb8 = controlMPC_s + 32;
controlMPC_FLOAT* controlMPC_llbbyslb8 = controlMPC_lbys + 32;
controlMPC_FLOAT controlMPC_rilb8[2];
controlMPC_FLOAT* controlMPC_dllbaff8 = controlMPC_dl_aff + 32;
controlMPC_FLOAT* controlMPC_dslbaff8 = controlMPC_ds_aff + 32;
controlMPC_FLOAT* controlMPC_dllbcc8 = controlMPC_dl_cc + 32;
controlMPC_FLOAT* controlMPC_dslbcc8 = controlMPC_ds_cc + 32;
controlMPC_FLOAT* controlMPC_ccrhsl8 = controlMPC_ccrhs + 32;
int controlMPC_ubIdx8[2] = {0, 1};
controlMPC_FLOAT* controlMPC_lub8 = controlMPC_l + 34;
controlMPC_FLOAT* controlMPC_sub8 = controlMPC_s + 34;
controlMPC_FLOAT* controlMPC_lubbysub8 = controlMPC_lbys + 34;
controlMPC_FLOAT controlMPC_riub8[2];
controlMPC_FLOAT* controlMPC_dlubaff8 = controlMPC_dl_aff + 34;
controlMPC_FLOAT* controlMPC_dsubaff8 = controlMPC_ds_aff + 34;
controlMPC_FLOAT* controlMPC_dlubcc8 = controlMPC_dl_cc + 34;
controlMPC_FLOAT* controlMPC_dsubcc8 = controlMPC_ds_cc + 34;
controlMPC_FLOAT* controlMPC_ccrhsub8 = controlMPC_ccrhs + 34;
controlMPC_FLOAT controlMPC_Phi8[2];
controlMPC_FLOAT controlMPC_W8[2];
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
controlMPC_LA_INITIALIZEVECTOR_18(controlMPC_z, 0);
controlMPC_LA_INITIALIZEVECTOR_16(controlMPC_v, 1);
controlMPC_LA_INITIALIZEVECTOR_36(controlMPC_l, 1);
controlMPC_LA_INITIALIZEVECTOR_36(controlMPC_s, 1);
info->mu = 0;
controlMPC_LA_DOTACC_36(controlMPC_l, controlMPC_s, &info->mu);
info->mu /= 36;
while( 1 ){
info->pobj = 0;
controlMPC_LA_DIAG_QUADFCN_2(params->H1, params->f1, controlMPC_z0, controlMPC_grad_cost0, &info->pobj);
controlMPC_LA_DIAG_QUADFCN_2(params->H2, params->f2, controlMPC_z1, controlMPC_grad_cost1, &info->pobj);
controlMPC_LA_DIAG_QUADFCN_2(params->H3, params->f3, controlMPC_z2, controlMPC_grad_cost2, &info->pobj);
controlMPC_LA_DIAG_QUADFCN_2(params->H4, params->f4, controlMPC_z3, controlMPC_grad_cost3, &info->pobj);
controlMPC_LA_DIAG_QUADFCN_2(params->H5, params->f5, controlMPC_z4, controlMPC_grad_cost4, &info->pobj);
controlMPC_LA_DIAG_QUADFCN_2(params->H6, params->f6, controlMPC_z5, controlMPC_grad_cost5, &info->pobj);
controlMPC_LA_DIAG_QUADFCN_2(params->H7, params->f7, controlMPC_z6, controlMPC_grad_cost6, &info->pobj);
controlMPC_LA_DIAG_QUADFCN_2(params->H8, params->f8, controlMPC_z7, controlMPC_grad_cost7, &info->pobj);
controlMPC_LA_DIAG_QUADFCN_2(params->H9, params->f9, controlMPC_z8, controlMPC_grad_cost8, &info->pobj);
info->res_eq = 0;
info->dgap = 0;
controlMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_2_2(controlMPC_C0, controlMPC_z0, controlMPC_D1, controlMPC_z1, controlMPC_c0, controlMPC_v0, controlMPC_re0, &info->dgap, &info->res_eq);
controlMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_2_2(controlMPC_C0, controlMPC_z1, controlMPC_D1, controlMPC_z2, controlMPC_c1, controlMPC_v1, controlMPC_re1, &info->dgap, &info->res_eq);
controlMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_2_2(controlMPC_C0, controlMPC_z2, controlMPC_D1, controlMPC_z3, controlMPC_c2, controlMPC_v2, controlMPC_re2, &info->dgap, &info->res_eq);
controlMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_2_2(controlMPC_C0, controlMPC_z3, controlMPC_D1, controlMPC_z4, controlMPC_c3, controlMPC_v3, controlMPC_re3, &info->dgap, &info->res_eq);
controlMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_2_2(controlMPC_C0, controlMPC_z4, controlMPC_D1, controlMPC_z5, controlMPC_c4, controlMPC_v4, controlMPC_re4, &info->dgap, &info->res_eq);
controlMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_2_2(controlMPC_C0, controlMPC_z5, controlMPC_D1, controlMPC_z6, controlMPC_c5, controlMPC_v5, controlMPC_re5, &info->dgap, &info->res_eq);
controlMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_2_2(controlMPC_C0, controlMPC_z6, controlMPC_D1, controlMPC_z7, controlMPC_c6, controlMPC_v6, controlMPC_re6, &info->dgap, &info->res_eq);
controlMPC_LA_DENSE_DIAGZERO_MVMSUB3_2_2_2(controlMPC_C0, controlMPC_z7, controlMPC_D1, controlMPC_z8, controlMPC_c7, controlMPC_v7, controlMPC_re7, &info->dgap, &info->res_eq);
controlMPC_LA_DENSE_MTVM_2_2(controlMPC_C0, controlMPC_v0, controlMPC_grad_eq0);
controlMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlMPC_C0, controlMPC_v1, controlMPC_D1, controlMPC_v0, controlMPC_grad_eq1);
controlMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlMPC_C0, controlMPC_v2, controlMPC_D1, controlMPC_v1, controlMPC_grad_eq2);
controlMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlMPC_C0, controlMPC_v3, controlMPC_D1, controlMPC_v2, controlMPC_grad_eq3);
controlMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlMPC_C0, controlMPC_v4, controlMPC_D1, controlMPC_v3, controlMPC_grad_eq4);
controlMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlMPC_C0, controlMPC_v5, controlMPC_D1, controlMPC_v4, controlMPC_grad_eq5);
controlMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlMPC_C0, controlMPC_v6, controlMPC_D1, controlMPC_v5, controlMPC_grad_eq6);
controlMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlMPC_C0, controlMPC_v7, controlMPC_D1, controlMPC_v6, controlMPC_grad_eq7);
controlMPC_LA_DIAGZERO_MTVM_2_2(controlMPC_D1, controlMPC_v7, controlMPC_grad_eq8);
info->res_ineq = 0;
controlMPC_LA_VSUBADD3_2(params->lb1, controlMPC_z0, controlMPC_lbIdx0, controlMPC_llb0, controlMPC_slb0, controlMPC_rilb0, &info->dgap, &info->res_ineq);
controlMPC_LA_VSUBADD2_2(controlMPC_z0, controlMPC_ubIdx0, params->ub1, controlMPC_lub0, controlMPC_sub0, controlMPC_riub0, &info->dgap, &info->res_ineq);
controlMPC_LA_VSUBADD3_2(params->lb2, controlMPC_z1, controlMPC_lbIdx1, controlMPC_llb1, controlMPC_slb1, controlMPC_rilb1, &info->dgap, &info->res_ineq);
controlMPC_LA_VSUBADD2_2(controlMPC_z1, controlMPC_ubIdx1, params->ub2, controlMPC_lub1, controlMPC_sub1, controlMPC_riub1, &info->dgap, &info->res_ineq);
controlMPC_LA_VSUBADD3_2(params->lb3, controlMPC_z2, controlMPC_lbIdx2, controlMPC_llb2, controlMPC_slb2, controlMPC_rilb2, &info->dgap, &info->res_ineq);
controlMPC_LA_VSUBADD2_2(controlMPC_z2, controlMPC_ubIdx2, params->ub3, controlMPC_lub2, controlMPC_sub2, controlMPC_riub2, &info->dgap, &info->res_ineq);
controlMPC_LA_VSUBADD3_2(params->lb4, controlMPC_z3, controlMPC_lbIdx3, controlMPC_llb3, controlMPC_slb3, controlMPC_rilb3, &info->dgap, &info->res_ineq);
controlMPC_LA_VSUBADD2_2(controlMPC_z3, controlMPC_ubIdx3, params->ub4, controlMPC_lub3, controlMPC_sub3, controlMPC_riub3, &info->dgap, &info->res_ineq);
controlMPC_LA_VSUBADD3_2(params->lb5, controlMPC_z4, controlMPC_lbIdx4, controlMPC_llb4, controlMPC_slb4, controlMPC_rilb4, &info->dgap, &info->res_ineq);
controlMPC_LA_VSUBADD2_2(controlMPC_z4, controlMPC_ubIdx4, params->ub5, controlMPC_lub4, controlMPC_sub4, controlMPC_riub4, &info->dgap, &info->res_ineq);
controlMPC_LA_VSUBADD3_2(params->lb6, controlMPC_z5, controlMPC_lbIdx5, controlMPC_llb5, controlMPC_slb5, controlMPC_rilb5, &info->dgap, &info->res_ineq);
controlMPC_LA_VSUBADD2_2(controlMPC_z5, controlMPC_ubIdx5, params->ub6, controlMPC_lub5, controlMPC_sub5, controlMPC_riub5, &info->dgap, &info->res_ineq);
controlMPC_LA_VSUBADD3_2(params->lb7, controlMPC_z6, controlMPC_lbIdx6, controlMPC_llb6, controlMPC_slb6, controlMPC_rilb6, &info->dgap, &info->res_ineq);
controlMPC_LA_VSUBADD2_2(controlMPC_z6, controlMPC_ubIdx6, params->ub7, controlMPC_lub6, controlMPC_sub6, controlMPC_riub6, &info->dgap, &info->res_ineq);
controlMPC_LA_VSUBADD3_2(params->lb8, controlMPC_z7, controlMPC_lbIdx7, controlMPC_llb7, controlMPC_slb7, controlMPC_rilb7, &info->dgap, &info->res_ineq);
controlMPC_LA_VSUBADD2_2(controlMPC_z7, controlMPC_ubIdx7, params->ub8, controlMPC_lub7, controlMPC_sub7, controlMPC_riub7, &info->dgap, &info->res_ineq);
controlMPC_LA_VSUBADD3_2(params->lb9, controlMPC_z8, controlMPC_lbIdx8, controlMPC_llb8, controlMPC_slb8, controlMPC_rilb8, &info->dgap, &info->res_ineq);
controlMPC_LA_VSUBADD2_2(controlMPC_z8, controlMPC_ubIdx8, params->ub9, controlMPC_lub8, controlMPC_sub8, controlMPC_riub8, &info->dgap, &info->res_ineq);
controlMPC_LA_INEQ_B_GRAD_2_2_2(controlMPC_lub0, controlMPC_sub0, controlMPC_riub0, controlMPC_llb0, controlMPC_slb0, controlMPC_rilb0, controlMPC_lbIdx0, controlMPC_ubIdx0, controlMPC_grad_ineq0, controlMPC_lubbysub0, controlMPC_llbbyslb0);
controlMPC_LA_INEQ_B_GRAD_2_2_2(controlMPC_lub1, controlMPC_sub1, controlMPC_riub1, controlMPC_llb1, controlMPC_slb1, controlMPC_rilb1, controlMPC_lbIdx1, controlMPC_ubIdx1, controlMPC_grad_ineq1, controlMPC_lubbysub1, controlMPC_llbbyslb1);
controlMPC_LA_INEQ_B_GRAD_2_2_2(controlMPC_lub2, controlMPC_sub2, controlMPC_riub2, controlMPC_llb2, controlMPC_slb2, controlMPC_rilb2, controlMPC_lbIdx2, controlMPC_ubIdx2, controlMPC_grad_ineq2, controlMPC_lubbysub2, controlMPC_llbbyslb2);
controlMPC_LA_INEQ_B_GRAD_2_2_2(controlMPC_lub3, controlMPC_sub3, controlMPC_riub3, controlMPC_llb3, controlMPC_slb3, controlMPC_rilb3, controlMPC_lbIdx3, controlMPC_ubIdx3, controlMPC_grad_ineq3, controlMPC_lubbysub3, controlMPC_llbbyslb3);
controlMPC_LA_INEQ_B_GRAD_2_2_2(controlMPC_lub4, controlMPC_sub4, controlMPC_riub4, controlMPC_llb4, controlMPC_slb4, controlMPC_rilb4, controlMPC_lbIdx4, controlMPC_ubIdx4, controlMPC_grad_ineq4, controlMPC_lubbysub4, controlMPC_llbbyslb4);
controlMPC_LA_INEQ_B_GRAD_2_2_2(controlMPC_lub5, controlMPC_sub5, controlMPC_riub5, controlMPC_llb5, controlMPC_slb5, controlMPC_rilb5, controlMPC_lbIdx5, controlMPC_ubIdx5, controlMPC_grad_ineq5, controlMPC_lubbysub5, controlMPC_llbbyslb5);
controlMPC_LA_INEQ_B_GRAD_2_2_2(controlMPC_lub6, controlMPC_sub6, controlMPC_riub6, controlMPC_llb6, controlMPC_slb6, controlMPC_rilb6, controlMPC_lbIdx6, controlMPC_ubIdx6, controlMPC_grad_ineq6, controlMPC_lubbysub6, controlMPC_llbbyslb6);
controlMPC_LA_INEQ_B_GRAD_2_2_2(controlMPC_lub7, controlMPC_sub7, controlMPC_riub7, controlMPC_llb7, controlMPC_slb7, controlMPC_rilb7, controlMPC_lbIdx7, controlMPC_ubIdx7, controlMPC_grad_ineq7, controlMPC_lubbysub7, controlMPC_llbbyslb7);
controlMPC_LA_INEQ_B_GRAD_2_2_2(controlMPC_lub8, controlMPC_sub8, controlMPC_riub8, controlMPC_llb8, controlMPC_slb8, controlMPC_rilb8, controlMPC_lbIdx8, controlMPC_ubIdx8, controlMPC_grad_ineq8, controlMPC_lubbysub8, controlMPC_llbbyslb8);
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
controlMPC_LA_VVADD3_18(controlMPC_grad_cost, controlMPC_grad_eq, controlMPC_grad_ineq, controlMPC_rd);
controlMPC_LA_DIAG_CHOL_ONELOOP_LBUB_2_2_2(params->H1, controlMPC_llbbyslb0, controlMPC_lbIdx0, controlMPC_lubbysub0, controlMPC_ubIdx0, controlMPC_Phi0);
controlMPC_LA_DIAG_MATRIXFORWARDSUB_2_2(controlMPC_Phi0, controlMPC_C0, controlMPC_V0);
controlMPC_LA_DIAG_FORWARDSUB_2(controlMPC_Phi0, controlMPC_rd0, controlMPC_Lbyrd0);
controlMPC_LA_DIAG_CHOL_ONELOOP_LBUB_2_2_2(params->H2, controlMPC_llbbyslb1, controlMPC_lbIdx1, controlMPC_lubbysub1, controlMPC_ubIdx1, controlMPC_Phi1);
controlMPC_LA_DIAG_MATRIXFORWARDSUB_2_2(controlMPC_Phi1, controlMPC_C0, controlMPC_V1);
controlMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_2(controlMPC_Phi1, controlMPC_D1, controlMPC_W1);
controlMPC_LA_DENSE_DIAGZERO_MMTM_2_2_2(controlMPC_W1, controlMPC_V1, controlMPC_Ysd1);
controlMPC_LA_DIAG_FORWARDSUB_2(controlMPC_Phi1, controlMPC_rd1, controlMPC_Lbyrd1);
controlMPC_LA_DIAG_CHOL_ONELOOP_LBUB_2_2_2(params->H3, controlMPC_llbbyslb2, controlMPC_lbIdx2, controlMPC_lubbysub2, controlMPC_ubIdx2, controlMPC_Phi2);
controlMPC_LA_DIAG_MATRIXFORWARDSUB_2_2(controlMPC_Phi2, controlMPC_C0, controlMPC_V2);
controlMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_2(controlMPC_Phi2, controlMPC_D1, controlMPC_W2);
controlMPC_LA_DENSE_DIAGZERO_MMTM_2_2_2(controlMPC_W2, controlMPC_V2, controlMPC_Ysd2);
controlMPC_LA_DIAG_FORWARDSUB_2(controlMPC_Phi2, controlMPC_rd2, controlMPC_Lbyrd2);
controlMPC_LA_DIAG_CHOL_ONELOOP_LBUB_2_2_2(params->H4, controlMPC_llbbyslb3, controlMPC_lbIdx3, controlMPC_lubbysub3, controlMPC_ubIdx3, controlMPC_Phi3);
controlMPC_LA_DIAG_MATRIXFORWARDSUB_2_2(controlMPC_Phi3, controlMPC_C0, controlMPC_V3);
controlMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_2(controlMPC_Phi3, controlMPC_D1, controlMPC_W3);
controlMPC_LA_DENSE_DIAGZERO_MMTM_2_2_2(controlMPC_W3, controlMPC_V3, controlMPC_Ysd3);
controlMPC_LA_DIAG_FORWARDSUB_2(controlMPC_Phi3, controlMPC_rd3, controlMPC_Lbyrd3);
controlMPC_LA_DIAG_CHOL_ONELOOP_LBUB_2_2_2(params->H5, controlMPC_llbbyslb4, controlMPC_lbIdx4, controlMPC_lubbysub4, controlMPC_ubIdx4, controlMPC_Phi4);
controlMPC_LA_DIAG_MATRIXFORWARDSUB_2_2(controlMPC_Phi4, controlMPC_C0, controlMPC_V4);
controlMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_2(controlMPC_Phi4, controlMPC_D1, controlMPC_W4);
controlMPC_LA_DENSE_DIAGZERO_MMTM_2_2_2(controlMPC_W4, controlMPC_V4, controlMPC_Ysd4);
controlMPC_LA_DIAG_FORWARDSUB_2(controlMPC_Phi4, controlMPC_rd4, controlMPC_Lbyrd4);
controlMPC_LA_DIAG_CHOL_ONELOOP_LBUB_2_2_2(params->H6, controlMPC_llbbyslb5, controlMPC_lbIdx5, controlMPC_lubbysub5, controlMPC_ubIdx5, controlMPC_Phi5);
controlMPC_LA_DIAG_MATRIXFORWARDSUB_2_2(controlMPC_Phi5, controlMPC_C0, controlMPC_V5);
controlMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_2(controlMPC_Phi5, controlMPC_D1, controlMPC_W5);
controlMPC_LA_DENSE_DIAGZERO_MMTM_2_2_2(controlMPC_W5, controlMPC_V5, controlMPC_Ysd5);
controlMPC_LA_DIAG_FORWARDSUB_2(controlMPC_Phi5, controlMPC_rd5, controlMPC_Lbyrd5);
controlMPC_LA_DIAG_CHOL_ONELOOP_LBUB_2_2_2(params->H7, controlMPC_llbbyslb6, controlMPC_lbIdx6, controlMPC_lubbysub6, controlMPC_ubIdx6, controlMPC_Phi6);
controlMPC_LA_DIAG_MATRIXFORWARDSUB_2_2(controlMPC_Phi6, controlMPC_C0, controlMPC_V6);
controlMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_2(controlMPC_Phi6, controlMPC_D1, controlMPC_W6);
controlMPC_LA_DENSE_DIAGZERO_MMTM_2_2_2(controlMPC_W6, controlMPC_V6, controlMPC_Ysd6);
controlMPC_LA_DIAG_FORWARDSUB_2(controlMPC_Phi6, controlMPC_rd6, controlMPC_Lbyrd6);
controlMPC_LA_DIAG_CHOL_ONELOOP_LBUB_2_2_2(params->H8, controlMPC_llbbyslb7, controlMPC_lbIdx7, controlMPC_lubbysub7, controlMPC_ubIdx7, controlMPC_Phi7);
controlMPC_LA_DIAG_MATRIXFORWARDSUB_2_2(controlMPC_Phi7, controlMPC_C0, controlMPC_V7);
controlMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_2(controlMPC_Phi7, controlMPC_D1, controlMPC_W7);
controlMPC_LA_DENSE_DIAGZERO_MMTM_2_2_2(controlMPC_W7, controlMPC_V7, controlMPC_Ysd7);
controlMPC_LA_DIAG_FORWARDSUB_2(controlMPC_Phi7, controlMPC_rd7, controlMPC_Lbyrd7);
controlMPC_LA_DIAG_CHOL_ONELOOP_LBUB_2_2_2(params->H9, controlMPC_llbbyslb8, controlMPC_lbIdx8, controlMPC_lubbysub8, controlMPC_ubIdx8, controlMPC_Phi8);
controlMPC_LA_DIAG_DIAGZERO_MATRIXTFORWARDSUB_2_2(controlMPC_Phi8, controlMPC_D1, controlMPC_W8);
controlMPC_LA_DIAG_FORWARDSUB_2(controlMPC_Phi8, controlMPC_rd8, controlMPC_Lbyrd8);
controlMPC_LA_DENSE_DIAGZERO_MMT2_2_2_2(controlMPC_V0, controlMPC_W1, controlMPC_Yd0);
controlMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_2_2(controlMPC_V0, controlMPC_Lbyrd0, controlMPC_W1, controlMPC_Lbyrd1, controlMPC_re0, controlMPC_beta0);
controlMPC_LA_DENSE_DIAGZERO_MMT2_2_2_2(controlMPC_V1, controlMPC_W2, controlMPC_Yd1);
controlMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_2_2(controlMPC_V1, controlMPC_Lbyrd1, controlMPC_W2, controlMPC_Lbyrd2, controlMPC_re1, controlMPC_beta1);
controlMPC_LA_DENSE_DIAGZERO_MMT2_2_2_2(controlMPC_V2, controlMPC_W3, controlMPC_Yd2);
controlMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_2_2(controlMPC_V2, controlMPC_Lbyrd2, controlMPC_W3, controlMPC_Lbyrd3, controlMPC_re2, controlMPC_beta2);
controlMPC_LA_DENSE_DIAGZERO_MMT2_2_2_2(controlMPC_V3, controlMPC_W4, controlMPC_Yd3);
controlMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_2_2(controlMPC_V3, controlMPC_Lbyrd3, controlMPC_W4, controlMPC_Lbyrd4, controlMPC_re3, controlMPC_beta3);
controlMPC_LA_DENSE_DIAGZERO_MMT2_2_2_2(controlMPC_V4, controlMPC_W5, controlMPC_Yd4);
controlMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_2_2(controlMPC_V4, controlMPC_Lbyrd4, controlMPC_W5, controlMPC_Lbyrd5, controlMPC_re4, controlMPC_beta4);
controlMPC_LA_DENSE_DIAGZERO_MMT2_2_2_2(controlMPC_V5, controlMPC_W6, controlMPC_Yd5);
controlMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_2_2(controlMPC_V5, controlMPC_Lbyrd5, controlMPC_W6, controlMPC_Lbyrd6, controlMPC_re5, controlMPC_beta5);
controlMPC_LA_DENSE_DIAGZERO_MMT2_2_2_2(controlMPC_V6, controlMPC_W7, controlMPC_Yd6);
controlMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_2_2(controlMPC_V6, controlMPC_Lbyrd6, controlMPC_W7, controlMPC_Lbyrd7, controlMPC_re6, controlMPC_beta6);
controlMPC_LA_DENSE_DIAGZERO_MMT2_2_2_2(controlMPC_V7, controlMPC_W8, controlMPC_Yd7);
controlMPC_LA_DENSE_DIAGZERO_2MVMSUB2_2_2_2(controlMPC_V7, controlMPC_Lbyrd7, controlMPC_W8, controlMPC_Lbyrd8, controlMPC_re7, controlMPC_beta7);
controlMPC_LA_DENSE_CHOL_2(controlMPC_Yd0, controlMPC_Ld0);
controlMPC_LA_DENSE_FORWARDSUB_2(controlMPC_Ld0, controlMPC_beta0, controlMPC_yy0);
controlMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(controlMPC_Ld0, controlMPC_Ysd1, controlMPC_Lsd1);
controlMPC_LA_DENSE_MMTSUB_2_2(controlMPC_Lsd1, controlMPC_Yd1);
controlMPC_LA_DENSE_CHOL_2(controlMPC_Yd1, controlMPC_Ld1);
controlMPC_LA_DENSE_MVMSUB1_2_2(controlMPC_Lsd1, controlMPC_yy0, controlMPC_beta1, controlMPC_bmy1);
controlMPC_LA_DENSE_FORWARDSUB_2(controlMPC_Ld1, controlMPC_bmy1, controlMPC_yy1);
controlMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(controlMPC_Ld1, controlMPC_Ysd2, controlMPC_Lsd2);
controlMPC_LA_DENSE_MMTSUB_2_2(controlMPC_Lsd2, controlMPC_Yd2);
controlMPC_LA_DENSE_CHOL_2(controlMPC_Yd2, controlMPC_Ld2);
controlMPC_LA_DENSE_MVMSUB1_2_2(controlMPC_Lsd2, controlMPC_yy1, controlMPC_beta2, controlMPC_bmy2);
controlMPC_LA_DENSE_FORWARDSUB_2(controlMPC_Ld2, controlMPC_bmy2, controlMPC_yy2);
controlMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(controlMPC_Ld2, controlMPC_Ysd3, controlMPC_Lsd3);
controlMPC_LA_DENSE_MMTSUB_2_2(controlMPC_Lsd3, controlMPC_Yd3);
controlMPC_LA_DENSE_CHOL_2(controlMPC_Yd3, controlMPC_Ld3);
controlMPC_LA_DENSE_MVMSUB1_2_2(controlMPC_Lsd3, controlMPC_yy2, controlMPC_beta3, controlMPC_bmy3);
controlMPC_LA_DENSE_FORWARDSUB_2(controlMPC_Ld3, controlMPC_bmy3, controlMPC_yy3);
controlMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(controlMPC_Ld3, controlMPC_Ysd4, controlMPC_Lsd4);
controlMPC_LA_DENSE_MMTSUB_2_2(controlMPC_Lsd4, controlMPC_Yd4);
controlMPC_LA_DENSE_CHOL_2(controlMPC_Yd4, controlMPC_Ld4);
controlMPC_LA_DENSE_MVMSUB1_2_2(controlMPC_Lsd4, controlMPC_yy3, controlMPC_beta4, controlMPC_bmy4);
controlMPC_LA_DENSE_FORWARDSUB_2(controlMPC_Ld4, controlMPC_bmy4, controlMPC_yy4);
controlMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(controlMPC_Ld4, controlMPC_Ysd5, controlMPC_Lsd5);
controlMPC_LA_DENSE_MMTSUB_2_2(controlMPC_Lsd5, controlMPC_Yd5);
controlMPC_LA_DENSE_CHOL_2(controlMPC_Yd5, controlMPC_Ld5);
controlMPC_LA_DENSE_MVMSUB1_2_2(controlMPC_Lsd5, controlMPC_yy4, controlMPC_beta5, controlMPC_bmy5);
controlMPC_LA_DENSE_FORWARDSUB_2(controlMPC_Ld5, controlMPC_bmy5, controlMPC_yy5);
controlMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(controlMPC_Ld5, controlMPC_Ysd6, controlMPC_Lsd6);
controlMPC_LA_DENSE_MMTSUB_2_2(controlMPC_Lsd6, controlMPC_Yd6);
controlMPC_LA_DENSE_CHOL_2(controlMPC_Yd6, controlMPC_Ld6);
controlMPC_LA_DENSE_MVMSUB1_2_2(controlMPC_Lsd6, controlMPC_yy5, controlMPC_beta6, controlMPC_bmy6);
controlMPC_LA_DENSE_FORWARDSUB_2(controlMPC_Ld6, controlMPC_bmy6, controlMPC_yy6);
controlMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(controlMPC_Ld6, controlMPC_Ysd7, controlMPC_Lsd7);
controlMPC_LA_DENSE_MMTSUB_2_2(controlMPC_Lsd7, controlMPC_Yd7);
controlMPC_LA_DENSE_CHOL_2(controlMPC_Yd7, controlMPC_Ld7);
controlMPC_LA_DENSE_MVMSUB1_2_2(controlMPC_Lsd7, controlMPC_yy6, controlMPC_beta7, controlMPC_bmy7);
controlMPC_LA_DENSE_FORWARDSUB_2(controlMPC_Ld7, controlMPC_bmy7, controlMPC_yy7);
controlMPC_LA_DENSE_BACKWARDSUB_2(controlMPC_Ld7, controlMPC_yy7, controlMPC_dvaff7);
controlMPC_LA_DENSE_MTVMSUB_2_2(controlMPC_Lsd7, controlMPC_dvaff7, controlMPC_yy6, controlMPC_bmy6);
controlMPC_LA_DENSE_BACKWARDSUB_2(controlMPC_Ld6, controlMPC_bmy6, controlMPC_dvaff6);
controlMPC_LA_DENSE_MTVMSUB_2_2(controlMPC_Lsd6, controlMPC_dvaff6, controlMPC_yy5, controlMPC_bmy5);
controlMPC_LA_DENSE_BACKWARDSUB_2(controlMPC_Ld5, controlMPC_bmy5, controlMPC_dvaff5);
controlMPC_LA_DENSE_MTVMSUB_2_2(controlMPC_Lsd5, controlMPC_dvaff5, controlMPC_yy4, controlMPC_bmy4);
controlMPC_LA_DENSE_BACKWARDSUB_2(controlMPC_Ld4, controlMPC_bmy4, controlMPC_dvaff4);
controlMPC_LA_DENSE_MTVMSUB_2_2(controlMPC_Lsd4, controlMPC_dvaff4, controlMPC_yy3, controlMPC_bmy3);
controlMPC_LA_DENSE_BACKWARDSUB_2(controlMPC_Ld3, controlMPC_bmy3, controlMPC_dvaff3);
controlMPC_LA_DENSE_MTVMSUB_2_2(controlMPC_Lsd3, controlMPC_dvaff3, controlMPC_yy2, controlMPC_bmy2);
controlMPC_LA_DENSE_BACKWARDSUB_2(controlMPC_Ld2, controlMPC_bmy2, controlMPC_dvaff2);
controlMPC_LA_DENSE_MTVMSUB_2_2(controlMPC_Lsd2, controlMPC_dvaff2, controlMPC_yy1, controlMPC_bmy1);
controlMPC_LA_DENSE_BACKWARDSUB_2(controlMPC_Ld1, controlMPC_bmy1, controlMPC_dvaff1);
controlMPC_LA_DENSE_MTVMSUB_2_2(controlMPC_Lsd1, controlMPC_dvaff1, controlMPC_yy0, controlMPC_bmy0);
controlMPC_LA_DENSE_BACKWARDSUB_2(controlMPC_Ld0, controlMPC_bmy0, controlMPC_dvaff0);
controlMPC_LA_DENSE_MTVM_2_2(controlMPC_C0, controlMPC_dvaff0, controlMPC_grad_eq0);
controlMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlMPC_C0, controlMPC_dvaff1, controlMPC_D1, controlMPC_dvaff0, controlMPC_grad_eq1);
controlMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlMPC_C0, controlMPC_dvaff2, controlMPC_D1, controlMPC_dvaff1, controlMPC_grad_eq2);
controlMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlMPC_C0, controlMPC_dvaff3, controlMPC_D1, controlMPC_dvaff2, controlMPC_grad_eq3);
controlMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlMPC_C0, controlMPC_dvaff4, controlMPC_D1, controlMPC_dvaff3, controlMPC_grad_eq4);
controlMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlMPC_C0, controlMPC_dvaff5, controlMPC_D1, controlMPC_dvaff4, controlMPC_grad_eq5);
controlMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlMPC_C0, controlMPC_dvaff6, controlMPC_D1, controlMPC_dvaff5, controlMPC_grad_eq6);
controlMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlMPC_C0, controlMPC_dvaff7, controlMPC_D1, controlMPC_dvaff6, controlMPC_grad_eq7);
controlMPC_LA_DIAGZERO_MTVM_2_2(controlMPC_D1, controlMPC_dvaff7, controlMPC_grad_eq8);
controlMPC_LA_VSUB2_18(controlMPC_rd, controlMPC_grad_eq, controlMPC_rd);
controlMPC_LA_DIAG_FORWARDBACKWARDSUB_2(controlMPC_Phi0, controlMPC_rd0, controlMPC_dzaff0);
controlMPC_LA_DIAG_FORWARDBACKWARDSUB_2(controlMPC_Phi1, controlMPC_rd1, controlMPC_dzaff1);
controlMPC_LA_DIAG_FORWARDBACKWARDSUB_2(controlMPC_Phi2, controlMPC_rd2, controlMPC_dzaff2);
controlMPC_LA_DIAG_FORWARDBACKWARDSUB_2(controlMPC_Phi3, controlMPC_rd3, controlMPC_dzaff3);
controlMPC_LA_DIAG_FORWARDBACKWARDSUB_2(controlMPC_Phi4, controlMPC_rd4, controlMPC_dzaff4);
controlMPC_LA_DIAG_FORWARDBACKWARDSUB_2(controlMPC_Phi5, controlMPC_rd5, controlMPC_dzaff5);
controlMPC_LA_DIAG_FORWARDBACKWARDSUB_2(controlMPC_Phi6, controlMPC_rd6, controlMPC_dzaff6);
controlMPC_LA_DIAG_FORWARDBACKWARDSUB_2(controlMPC_Phi7, controlMPC_rd7, controlMPC_dzaff7);
controlMPC_LA_DIAG_FORWARDBACKWARDSUB_2(controlMPC_Phi8, controlMPC_rd8, controlMPC_dzaff8);
controlMPC_LA_VSUB_INDEXED_2(controlMPC_dzaff0, controlMPC_lbIdx0, controlMPC_rilb0, controlMPC_dslbaff0);
controlMPC_LA_VSUB3_2(controlMPC_llbbyslb0, controlMPC_dslbaff0, controlMPC_llb0, controlMPC_dllbaff0);
controlMPC_LA_VSUB2_INDEXED_2(controlMPC_riub0, controlMPC_dzaff0, controlMPC_ubIdx0, controlMPC_dsubaff0);
controlMPC_LA_VSUB3_2(controlMPC_lubbysub0, controlMPC_dsubaff0, controlMPC_lub0, controlMPC_dlubaff0);
controlMPC_LA_VSUB_INDEXED_2(controlMPC_dzaff1, controlMPC_lbIdx1, controlMPC_rilb1, controlMPC_dslbaff1);
controlMPC_LA_VSUB3_2(controlMPC_llbbyslb1, controlMPC_dslbaff1, controlMPC_llb1, controlMPC_dllbaff1);
controlMPC_LA_VSUB2_INDEXED_2(controlMPC_riub1, controlMPC_dzaff1, controlMPC_ubIdx1, controlMPC_dsubaff1);
controlMPC_LA_VSUB3_2(controlMPC_lubbysub1, controlMPC_dsubaff1, controlMPC_lub1, controlMPC_dlubaff1);
controlMPC_LA_VSUB_INDEXED_2(controlMPC_dzaff2, controlMPC_lbIdx2, controlMPC_rilb2, controlMPC_dslbaff2);
controlMPC_LA_VSUB3_2(controlMPC_llbbyslb2, controlMPC_dslbaff2, controlMPC_llb2, controlMPC_dllbaff2);
controlMPC_LA_VSUB2_INDEXED_2(controlMPC_riub2, controlMPC_dzaff2, controlMPC_ubIdx2, controlMPC_dsubaff2);
controlMPC_LA_VSUB3_2(controlMPC_lubbysub2, controlMPC_dsubaff2, controlMPC_lub2, controlMPC_dlubaff2);
controlMPC_LA_VSUB_INDEXED_2(controlMPC_dzaff3, controlMPC_lbIdx3, controlMPC_rilb3, controlMPC_dslbaff3);
controlMPC_LA_VSUB3_2(controlMPC_llbbyslb3, controlMPC_dslbaff3, controlMPC_llb3, controlMPC_dllbaff3);
controlMPC_LA_VSUB2_INDEXED_2(controlMPC_riub3, controlMPC_dzaff3, controlMPC_ubIdx3, controlMPC_dsubaff3);
controlMPC_LA_VSUB3_2(controlMPC_lubbysub3, controlMPC_dsubaff3, controlMPC_lub3, controlMPC_dlubaff3);
controlMPC_LA_VSUB_INDEXED_2(controlMPC_dzaff4, controlMPC_lbIdx4, controlMPC_rilb4, controlMPC_dslbaff4);
controlMPC_LA_VSUB3_2(controlMPC_llbbyslb4, controlMPC_dslbaff4, controlMPC_llb4, controlMPC_dllbaff4);
controlMPC_LA_VSUB2_INDEXED_2(controlMPC_riub4, controlMPC_dzaff4, controlMPC_ubIdx4, controlMPC_dsubaff4);
controlMPC_LA_VSUB3_2(controlMPC_lubbysub4, controlMPC_dsubaff4, controlMPC_lub4, controlMPC_dlubaff4);
controlMPC_LA_VSUB_INDEXED_2(controlMPC_dzaff5, controlMPC_lbIdx5, controlMPC_rilb5, controlMPC_dslbaff5);
controlMPC_LA_VSUB3_2(controlMPC_llbbyslb5, controlMPC_dslbaff5, controlMPC_llb5, controlMPC_dllbaff5);
controlMPC_LA_VSUB2_INDEXED_2(controlMPC_riub5, controlMPC_dzaff5, controlMPC_ubIdx5, controlMPC_dsubaff5);
controlMPC_LA_VSUB3_2(controlMPC_lubbysub5, controlMPC_dsubaff5, controlMPC_lub5, controlMPC_dlubaff5);
controlMPC_LA_VSUB_INDEXED_2(controlMPC_dzaff6, controlMPC_lbIdx6, controlMPC_rilb6, controlMPC_dslbaff6);
controlMPC_LA_VSUB3_2(controlMPC_llbbyslb6, controlMPC_dslbaff6, controlMPC_llb6, controlMPC_dllbaff6);
controlMPC_LA_VSUB2_INDEXED_2(controlMPC_riub6, controlMPC_dzaff6, controlMPC_ubIdx6, controlMPC_dsubaff6);
controlMPC_LA_VSUB3_2(controlMPC_lubbysub6, controlMPC_dsubaff6, controlMPC_lub6, controlMPC_dlubaff6);
controlMPC_LA_VSUB_INDEXED_2(controlMPC_dzaff7, controlMPC_lbIdx7, controlMPC_rilb7, controlMPC_dslbaff7);
controlMPC_LA_VSUB3_2(controlMPC_llbbyslb7, controlMPC_dslbaff7, controlMPC_llb7, controlMPC_dllbaff7);
controlMPC_LA_VSUB2_INDEXED_2(controlMPC_riub7, controlMPC_dzaff7, controlMPC_ubIdx7, controlMPC_dsubaff7);
controlMPC_LA_VSUB3_2(controlMPC_lubbysub7, controlMPC_dsubaff7, controlMPC_lub7, controlMPC_dlubaff7);
controlMPC_LA_VSUB_INDEXED_2(controlMPC_dzaff8, controlMPC_lbIdx8, controlMPC_rilb8, controlMPC_dslbaff8);
controlMPC_LA_VSUB3_2(controlMPC_llbbyslb8, controlMPC_dslbaff8, controlMPC_llb8, controlMPC_dllbaff8);
controlMPC_LA_VSUB2_INDEXED_2(controlMPC_riub8, controlMPC_dzaff8, controlMPC_ubIdx8, controlMPC_dsubaff8);
controlMPC_LA_VSUB3_2(controlMPC_lubbysub8, controlMPC_dsubaff8, controlMPC_lub8, controlMPC_dlubaff8);
info->lsit_aff = controlMPC_LINESEARCH_BACKTRACKING_AFFINE(controlMPC_l, controlMPC_s, controlMPC_dl_aff, controlMPC_ds_aff, &info->step_aff, &info->mu_aff);
if( info->lsit_aff == controlMPC_NOPROGRESS ){
exitcode = controlMPC_NOPROGRESS; break;
}
sigma_3rdroot = info->mu_aff / info->mu;
info->sigma = sigma_3rdroot*sigma_3rdroot*sigma_3rdroot;
musigma = info->mu * info->sigma;
controlMPC_LA_VSUB5_36(controlMPC_ds_aff, controlMPC_dl_aff, musigma, controlMPC_ccrhs);
controlMPC_LA_VSUB6_INDEXED_2_2_2(controlMPC_ccrhsub0, controlMPC_sub0, controlMPC_ubIdx0, controlMPC_ccrhsl0, controlMPC_slb0, controlMPC_lbIdx0, controlMPC_rd0);
controlMPC_LA_VSUB6_INDEXED_2_2_2(controlMPC_ccrhsub1, controlMPC_sub1, controlMPC_ubIdx1, controlMPC_ccrhsl1, controlMPC_slb1, controlMPC_lbIdx1, controlMPC_rd1);
controlMPC_LA_DIAG_FORWARDSUB_2(controlMPC_Phi0, controlMPC_rd0, controlMPC_Lbyrd0);
controlMPC_LA_DIAG_FORWARDSUB_2(controlMPC_Phi1, controlMPC_rd1, controlMPC_Lbyrd1);
controlMPC_LA_DENSE_DIAGZERO_2MVMADD_2_2_2(controlMPC_V0, controlMPC_Lbyrd0, controlMPC_W1, controlMPC_Lbyrd1, controlMPC_beta0);
controlMPC_LA_DENSE_FORWARDSUB_2(controlMPC_Ld0, controlMPC_beta0, controlMPC_yy0);
controlMPC_LA_VSUB6_INDEXED_2_2_2(controlMPC_ccrhsub2, controlMPC_sub2, controlMPC_ubIdx2, controlMPC_ccrhsl2, controlMPC_slb2, controlMPC_lbIdx2, controlMPC_rd2);
controlMPC_LA_DIAG_FORWARDSUB_2(controlMPC_Phi2, controlMPC_rd2, controlMPC_Lbyrd2);
controlMPC_LA_DENSE_DIAGZERO_2MVMADD_2_2_2(controlMPC_V1, controlMPC_Lbyrd1, controlMPC_W2, controlMPC_Lbyrd2, controlMPC_beta1);
controlMPC_LA_DENSE_MVMSUB1_2_2(controlMPC_Lsd1, controlMPC_yy0, controlMPC_beta1, controlMPC_bmy1);
controlMPC_LA_DENSE_FORWARDSUB_2(controlMPC_Ld1, controlMPC_bmy1, controlMPC_yy1);
controlMPC_LA_VSUB6_INDEXED_2_2_2(controlMPC_ccrhsub3, controlMPC_sub3, controlMPC_ubIdx3, controlMPC_ccrhsl3, controlMPC_slb3, controlMPC_lbIdx3, controlMPC_rd3);
controlMPC_LA_DIAG_FORWARDSUB_2(controlMPC_Phi3, controlMPC_rd3, controlMPC_Lbyrd3);
controlMPC_LA_DENSE_DIAGZERO_2MVMADD_2_2_2(controlMPC_V2, controlMPC_Lbyrd2, controlMPC_W3, controlMPC_Lbyrd3, controlMPC_beta2);
controlMPC_LA_DENSE_MVMSUB1_2_2(controlMPC_Lsd2, controlMPC_yy1, controlMPC_beta2, controlMPC_bmy2);
controlMPC_LA_DENSE_FORWARDSUB_2(controlMPC_Ld2, controlMPC_bmy2, controlMPC_yy2);
controlMPC_LA_VSUB6_INDEXED_2_2_2(controlMPC_ccrhsub4, controlMPC_sub4, controlMPC_ubIdx4, controlMPC_ccrhsl4, controlMPC_slb4, controlMPC_lbIdx4, controlMPC_rd4);
controlMPC_LA_DIAG_FORWARDSUB_2(controlMPC_Phi4, controlMPC_rd4, controlMPC_Lbyrd4);
controlMPC_LA_DENSE_DIAGZERO_2MVMADD_2_2_2(controlMPC_V3, controlMPC_Lbyrd3, controlMPC_W4, controlMPC_Lbyrd4, controlMPC_beta3);
controlMPC_LA_DENSE_MVMSUB1_2_2(controlMPC_Lsd3, controlMPC_yy2, controlMPC_beta3, controlMPC_bmy3);
controlMPC_LA_DENSE_FORWARDSUB_2(controlMPC_Ld3, controlMPC_bmy3, controlMPC_yy3);
controlMPC_LA_VSUB6_INDEXED_2_2_2(controlMPC_ccrhsub5, controlMPC_sub5, controlMPC_ubIdx5, controlMPC_ccrhsl5, controlMPC_slb5, controlMPC_lbIdx5, controlMPC_rd5);
controlMPC_LA_DIAG_FORWARDSUB_2(controlMPC_Phi5, controlMPC_rd5, controlMPC_Lbyrd5);
controlMPC_LA_DENSE_DIAGZERO_2MVMADD_2_2_2(controlMPC_V4, controlMPC_Lbyrd4, controlMPC_W5, controlMPC_Lbyrd5, controlMPC_beta4);
controlMPC_LA_DENSE_MVMSUB1_2_2(controlMPC_Lsd4, controlMPC_yy3, controlMPC_beta4, controlMPC_bmy4);
controlMPC_LA_DENSE_FORWARDSUB_2(controlMPC_Ld4, controlMPC_bmy4, controlMPC_yy4);
controlMPC_LA_VSUB6_INDEXED_2_2_2(controlMPC_ccrhsub6, controlMPC_sub6, controlMPC_ubIdx6, controlMPC_ccrhsl6, controlMPC_slb6, controlMPC_lbIdx6, controlMPC_rd6);
controlMPC_LA_DIAG_FORWARDSUB_2(controlMPC_Phi6, controlMPC_rd6, controlMPC_Lbyrd6);
controlMPC_LA_DENSE_DIAGZERO_2MVMADD_2_2_2(controlMPC_V5, controlMPC_Lbyrd5, controlMPC_W6, controlMPC_Lbyrd6, controlMPC_beta5);
controlMPC_LA_DENSE_MVMSUB1_2_2(controlMPC_Lsd5, controlMPC_yy4, controlMPC_beta5, controlMPC_bmy5);
controlMPC_LA_DENSE_FORWARDSUB_2(controlMPC_Ld5, controlMPC_bmy5, controlMPC_yy5);
controlMPC_LA_VSUB6_INDEXED_2_2_2(controlMPC_ccrhsub7, controlMPC_sub7, controlMPC_ubIdx7, controlMPC_ccrhsl7, controlMPC_slb7, controlMPC_lbIdx7, controlMPC_rd7);
controlMPC_LA_DIAG_FORWARDSUB_2(controlMPC_Phi7, controlMPC_rd7, controlMPC_Lbyrd7);
controlMPC_LA_DENSE_DIAGZERO_2MVMADD_2_2_2(controlMPC_V6, controlMPC_Lbyrd6, controlMPC_W7, controlMPC_Lbyrd7, controlMPC_beta6);
controlMPC_LA_DENSE_MVMSUB1_2_2(controlMPC_Lsd6, controlMPC_yy5, controlMPC_beta6, controlMPC_bmy6);
controlMPC_LA_DENSE_FORWARDSUB_2(controlMPC_Ld6, controlMPC_bmy6, controlMPC_yy6);
controlMPC_LA_VSUB6_INDEXED_2_2_2(controlMPC_ccrhsub8, controlMPC_sub8, controlMPC_ubIdx8, controlMPC_ccrhsl8, controlMPC_slb8, controlMPC_lbIdx8, controlMPC_rd8);
controlMPC_LA_DIAG_FORWARDSUB_2(controlMPC_Phi8, controlMPC_rd8, controlMPC_Lbyrd8);
controlMPC_LA_DENSE_DIAGZERO_2MVMADD_2_2_2(controlMPC_V7, controlMPC_Lbyrd7, controlMPC_W8, controlMPC_Lbyrd8, controlMPC_beta7);
controlMPC_LA_DENSE_MVMSUB1_2_2(controlMPC_Lsd7, controlMPC_yy6, controlMPC_beta7, controlMPC_bmy7);
controlMPC_LA_DENSE_FORWARDSUB_2(controlMPC_Ld7, controlMPC_bmy7, controlMPC_yy7);
controlMPC_LA_DENSE_BACKWARDSUB_2(controlMPC_Ld7, controlMPC_yy7, controlMPC_dvcc7);
controlMPC_LA_DENSE_MTVMSUB_2_2(controlMPC_Lsd7, controlMPC_dvcc7, controlMPC_yy6, controlMPC_bmy6);
controlMPC_LA_DENSE_BACKWARDSUB_2(controlMPC_Ld6, controlMPC_bmy6, controlMPC_dvcc6);
controlMPC_LA_DENSE_MTVMSUB_2_2(controlMPC_Lsd6, controlMPC_dvcc6, controlMPC_yy5, controlMPC_bmy5);
controlMPC_LA_DENSE_BACKWARDSUB_2(controlMPC_Ld5, controlMPC_bmy5, controlMPC_dvcc5);
controlMPC_LA_DENSE_MTVMSUB_2_2(controlMPC_Lsd5, controlMPC_dvcc5, controlMPC_yy4, controlMPC_bmy4);
controlMPC_LA_DENSE_BACKWARDSUB_2(controlMPC_Ld4, controlMPC_bmy4, controlMPC_dvcc4);
controlMPC_LA_DENSE_MTVMSUB_2_2(controlMPC_Lsd4, controlMPC_dvcc4, controlMPC_yy3, controlMPC_bmy3);
controlMPC_LA_DENSE_BACKWARDSUB_2(controlMPC_Ld3, controlMPC_bmy3, controlMPC_dvcc3);
controlMPC_LA_DENSE_MTVMSUB_2_2(controlMPC_Lsd3, controlMPC_dvcc3, controlMPC_yy2, controlMPC_bmy2);
controlMPC_LA_DENSE_BACKWARDSUB_2(controlMPC_Ld2, controlMPC_bmy2, controlMPC_dvcc2);
controlMPC_LA_DENSE_MTVMSUB_2_2(controlMPC_Lsd2, controlMPC_dvcc2, controlMPC_yy1, controlMPC_bmy1);
controlMPC_LA_DENSE_BACKWARDSUB_2(controlMPC_Ld1, controlMPC_bmy1, controlMPC_dvcc1);
controlMPC_LA_DENSE_MTVMSUB_2_2(controlMPC_Lsd1, controlMPC_dvcc1, controlMPC_yy0, controlMPC_bmy0);
controlMPC_LA_DENSE_BACKWARDSUB_2(controlMPC_Ld0, controlMPC_bmy0, controlMPC_dvcc0);
controlMPC_LA_DENSE_MTVM_2_2(controlMPC_C0, controlMPC_dvcc0, controlMPC_grad_eq0);
controlMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlMPC_C0, controlMPC_dvcc1, controlMPC_D1, controlMPC_dvcc0, controlMPC_grad_eq1);
controlMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlMPC_C0, controlMPC_dvcc2, controlMPC_D1, controlMPC_dvcc1, controlMPC_grad_eq2);
controlMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlMPC_C0, controlMPC_dvcc3, controlMPC_D1, controlMPC_dvcc2, controlMPC_grad_eq3);
controlMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlMPC_C0, controlMPC_dvcc4, controlMPC_D1, controlMPC_dvcc3, controlMPC_grad_eq4);
controlMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlMPC_C0, controlMPC_dvcc5, controlMPC_D1, controlMPC_dvcc4, controlMPC_grad_eq5);
controlMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlMPC_C0, controlMPC_dvcc6, controlMPC_D1, controlMPC_dvcc5, controlMPC_grad_eq6);
controlMPC_LA_DENSE_DIAGZERO_MTVM2_2_2_2(controlMPC_C0, controlMPC_dvcc7, controlMPC_D1, controlMPC_dvcc6, controlMPC_grad_eq7);
controlMPC_LA_DIAGZERO_MTVM_2_2(controlMPC_D1, controlMPC_dvcc7, controlMPC_grad_eq8);
controlMPC_LA_VSUB_18(controlMPC_rd, controlMPC_grad_eq, controlMPC_rd);
controlMPC_LA_DIAG_FORWARDBACKWARDSUB_2(controlMPC_Phi0, controlMPC_rd0, controlMPC_dzcc0);
controlMPC_LA_DIAG_FORWARDBACKWARDSUB_2(controlMPC_Phi1, controlMPC_rd1, controlMPC_dzcc1);
controlMPC_LA_DIAG_FORWARDBACKWARDSUB_2(controlMPC_Phi2, controlMPC_rd2, controlMPC_dzcc2);
controlMPC_LA_DIAG_FORWARDBACKWARDSUB_2(controlMPC_Phi3, controlMPC_rd3, controlMPC_dzcc3);
controlMPC_LA_DIAG_FORWARDBACKWARDSUB_2(controlMPC_Phi4, controlMPC_rd4, controlMPC_dzcc4);
controlMPC_LA_DIAG_FORWARDBACKWARDSUB_2(controlMPC_Phi5, controlMPC_rd5, controlMPC_dzcc5);
controlMPC_LA_DIAG_FORWARDBACKWARDSUB_2(controlMPC_Phi6, controlMPC_rd6, controlMPC_dzcc6);
controlMPC_LA_DIAG_FORWARDBACKWARDSUB_2(controlMPC_Phi7, controlMPC_rd7, controlMPC_dzcc7);
controlMPC_LA_DIAG_FORWARDBACKWARDSUB_2(controlMPC_Phi8, controlMPC_rd8, controlMPC_dzcc8);
controlMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_2(controlMPC_ccrhsl0, controlMPC_slb0, controlMPC_llbbyslb0, controlMPC_dzcc0, controlMPC_lbIdx0, controlMPC_dllbcc0);
controlMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_2(controlMPC_ccrhsub0, controlMPC_sub0, controlMPC_lubbysub0, controlMPC_dzcc0, controlMPC_ubIdx0, controlMPC_dlubcc0);
controlMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_2(controlMPC_ccrhsl1, controlMPC_slb1, controlMPC_llbbyslb1, controlMPC_dzcc1, controlMPC_lbIdx1, controlMPC_dllbcc1);
controlMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_2(controlMPC_ccrhsub1, controlMPC_sub1, controlMPC_lubbysub1, controlMPC_dzcc1, controlMPC_ubIdx1, controlMPC_dlubcc1);
controlMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_2(controlMPC_ccrhsl2, controlMPC_slb2, controlMPC_llbbyslb2, controlMPC_dzcc2, controlMPC_lbIdx2, controlMPC_dllbcc2);
controlMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_2(controlMPC_ccrhsub2, controlMPC_sub2, controlMPC_lubbysub2, controlMPC_dzcc2, controlMPC_ubIdx2, controlMPC_dlubcc2);
controlMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_2(controlMPC_ccrhsl3, controlMPC_slb3, controlMPC_llbbyslb3, controlMPC_dzcc3, controlMPC_lbIdx3, controlMPC_dllbcc3);
controlMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_2(controlMPC_ccrhsub3, controlMPC_sub3, controlMPC_lubbysub3, controlMPC_dzcc3, controlMPC_ubIdx3, controlMPC_dlubcc3);
controlMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_2(controlMPC_ccrhsl4, controlMPC_slb4, controlMPC_llbbyslb4, controlMPC_dzcc4, controlMPC_lbIdx4, controlMPC_dllbcc4);
controlMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_2(controlMPC_ccrhsub4, controlMPC_sub4, controlMPC_lubbysub4, controlMPC_dzcc4, controlMPC_ubIdx4, controlMPC_dlubcc4);
controlMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_2(controlMPC_ccrhsl5, controlMPC_slb5, controlMPC_llbbyslb5, controlMPC_dzcc5, controlMPC_lbIdx5, controlMPC_dllbcc5);
controlMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_2(controlMPC_ccrhsub5, controlMPC_sub5, controlMPC_lubbysub5, controlMPC_dzcc5, controlMPC_ubIdx5, controlMPC_dlubcc5);
controlMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_2(controlMPC_ccrhsl6, controlMPC_slb6, controlMPC_llbbyslb6, controlMPC_dzcc6, controlMPC_lbIdx6, controlMPC_dllbcc6);
controlMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_2(controlMPC_ccrhsub6, controlMPC_sub6, controlMPC_lubbysub6, controlMPC_dzcc6, controlMPC_ubIdx6, controlMPC_dlubcc6);
controlMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_2(controlMPC_ccrhsl7, controlMPC_slb7, controlMPC_llbbyslb7, controlMPC_dzcc7, controlMPC_lbIdx7, controlMPC_dllbcc7);
controlMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_2(controlMPC_ccrhsub7, controlMPC_sub7, controlMPC_lubbysub7, controlMPC_dzcc7, controlMPC_ubIdx7, controlMPC_dlubcc7);
controlMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_2(controlMPC_ccrhsl8, controlMPC_slb8, controlMPC_llbbyslb8, controlMPC_dzcc8, controlMPC_lbIdx8, controlMPC_dllbcc8);
controlMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_2(controlMPC_ccrhsub8, controlMPC_sub8, controlMPC_lubbysub8, controlMPC_dzcc8, controlMPC_ubIdx8, controlMPC_dlubcc8);
controlMPC_LA_VSUB7_36(controlMPC_l, controlMPC_ccrhs, controlMPC_s, controlMPC_dl_cc, controlMPC_ds_cc);
controlMPC_LA_VADD_18(controlMPC_dz_cc, controlMPC_dz_aff);
controlMPC_LA_VADD_16(controlMPC_dv_cc, controlMPC_dv_aff);
controlMPC_LA_VADD_36(controlMPC_dl_cc, controlMPC_dl_aff);
controlMPC_LA_VADD_36(controlMPC_ds_cc, controlMPC_ds_aff);
info->lsit_cc = controlMPC_LINESEARCH_BACKTRACKING_COMBINED(controlMPC_z, controlMPC_v, controlMPC_l, controlMPC_s, controlMPC_dz_cc, controlMPC_dv_cc, controlMPC_dl_cc, controlMPC_ds_cc, &info->step_cc, &info->mu);
if( info->lsit_cc == controlMPC_NOPROGRESS ){
exitcode = controlMPC_NOPROGRESS; break;
}
info->it++;
}
output->z1[0] = controlMPC_z0[0];
output->z1[1] = controlMPC_z0[1];
output->z2[0] = controlMPC_z1[0];
output->z2[1] = controlMPC_z1[1];
output->z3[0] = controlMPC_z2[0];
output->z3[1] = controlMPC_z2[1];
output->z4[0] = controlMPC_z3[0];
output->z4[1] = controlMPC_z3[1];
output->z5[0] = controlMPC_z4[0];
output->z5[1] = controlMPC_z4[1];
output->z6[0] = controlMPC_z5[0];
output->z6[1] = controlMPC_z5[1];
output->z7[0] = controlMPC_z6[0];
output->z7[1] = controlMPC_z6[1];
output->z8[0] = controlMPC_z7[0];
output->z8[1] = controlMPC_z7[1];
output->z9[0] = controlMPC_z8[0];
output->z9[1] = controlMPC_z8[1];

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
