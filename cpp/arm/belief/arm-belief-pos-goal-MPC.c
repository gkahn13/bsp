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

#include "arm-belief-pos-goal-MPC.h"

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
 * Initializes a vector of length 1251 with a value.
 */
void beliefPenaltyMPC_LA_INITIALIZEVECTOR_1251(beliefPenaltyMPC_FLOAT* vec, beliefPenaltyMPC_FLOAT value)
{
	int i;
	for( i=0; i<1251; i++ )
	{
		vec[i] = value;
	}
}


/*
 * Initializes a vector of length 405 with a value.
 */
void beliefPenaltyMPC_LA_INITIALIZEVECTOR_405(beliefPenaltyMPC_FLOAT* vec, beliefPenaltyMPC_FLOAT value)
{
	int i;
	for( i=0; i<405; i++ )
	{
		vec[i] = value;
	}
}


/*
 * Initializes a vector of length 1746 with a value.
 */
void beliefPenaltyMPC_LA_INITIALIZEVECTOR_1746(beliefPenaltyMPC_FLOAT* vec, beliefPenaltyMPC_FLOAT value)
{
	int i;
	for( i=0; i<1746; i++ )
	{
		vec[i] = value;
	}
}


/* 
 * Calculates a dot product and adds it to a variable: z += x'*y; 
 * This function is for vectors of length 1746.
 */
void beliefPenaltyMPC_LA_DOTACC_1746(beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *y, beliefPenaltyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<1746; i++ ){
		*z += x[i]*y[i];
	}
}


/*
 * Calculates the gradient and the value for a quadratic function 0.5*z'*H*z + f'*z
 *
 * INPUTS:     H  - Symmetric Hessian, dense matrix of size [87 x 87]
 *             f  - column vector of size 87
 *             z  - column vector of size 87
 *
 * OUTPUTS: grad  - gradient at z (= H*z + f), column vector of size 87
 *          value <-- value + 0.5*z'*H*z + f'*z (value will be modified)
 */
void beliefPenaltyMPC_LA_DENSE_QUADFCN_87(beliefPenaltyMPC_FLOAT* H, beliefPenaltyMPC_FLOAT* f, beliefPenaltyMPC_FLOAT* z, beliefPenaltyMPC_FLOAT* grad, beliefPenaltyMPC_FLOAT* value)
{
	int i;
	int j;
	int k = 0;
	beliefPenaltyMPC_FLOAT hz;	
	for( i=0; i<87; i++){
		hz = 0;
		for( j=0; j<87; j++ )
		{
			hz += H[k++]*z[j];
		}
		grad[i] = hz + f[i];
		*value += 0.5*hz*z[i] + f[i]*z[i];
	}
}


/*
 * Calculates the gradient and the value for a quadratic function 0.5*z'*H*z + f'*z
 *
 * INPUTS:     H  - Symmetric Hessian, dense matrix of size [33 x 33]
 *             f  - column vector of size 33
 *             z  - column vector of size 33
 *
 * OUTPUTS: grad  - gradient at z (= H*z + f), column vector of size 33
 *          value <-- value + 0.5*z'*H*z + f'*z (value will be modified)
 */
void beliefPenaltyMPC_LA_DENSE_QUADFCN_33(beliefPenaltyMPC_FLOAT* H, beliefPenaltyMPC_FLOAT* f, beliefPenaltyMPC_FLOAT* z, beliefPenaltyMPC_FLOAT* grad, beliefPenaltyMPC_FLOAT* value)
{
	int i;
	int j;
	int k = 0;
	beliefPenaltyMPC_FLOAT hz;	
	for( i=0; i<33; i++){
		hz = 0;
		for( j=0; j<33; j++ )
		{
			hz += H[k++]*z[j];
		}
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
void beliefPenaltyMPC_LA_DIAGZERO_MVMSUB6_27(beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *b, beliefPenaltyMPC_FLOAT *l, beliefPenaltyMPC_FLOAT *r, beliefPenaltyMPC_FLOAT *z, beliefPenaltyMPC_FLOAT *y)
{
	int i;
	beliefPenaltyMPC_FLOAT Bu[27];
	beliefPenaltyMPC_FLOAT norm = *y;
	beliefPenaltyMPC_FLOAT lr = 0;

	/* do A*x + B*u first */
	for( i=0; i<27; i++ ){
		Bu[i] = B[i]*u[i];
	}	

	for( i=0; i<27; i++ ){
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
void beliefPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_27_87_87(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *b, beliefPenaltyMPC_FLOAT *l, beliefPenaltyMPC_FLOAT *r, beliefPenaltyMPC_FLOAT *z, beliefPenaltyMPC_FLOAT *y)
{
	int i;
	int j;
	int k = 0;
	beliefPenaltyMPC_FLOAT AxBu[27];
	beliefPenaltyMPC_FLOAT norm = *y;
	beliefPenaltyMPC_FLOAT lr = 0;

	/* do A*x + B*u first */
	for( i=0; i<27; i++ ){
		AxBu[i] = A[k++]*x[0] + B[i]*u[i];
	}	

	for( j=1; j<87; j++ ){		
		for( i=0; i<27; i++ ){
			AxBu[i] += A[k++]*x[j];
		}
	}

	for( i=0; i<27; i++ ){
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
void beliefPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_27_87_33(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *b, beliefPenaltyMPC_FLOAT *l, beliefPenaltyMPC_FLOAT *r, beliefPenaltyMPC_FLOAT *z, beliefPenaltyMPC_FLOAT *y)
{
	int i;
	int j;
	int k = 0;
	beliefPenaltyMPC_FLOAT AxBu[27];
	beliefPenaltyMPC_FLOAT norm = *y;
	beliefPenaltyMPC_FLOAT lr = 0;

	/* do A*x + B*u first */
	for( i=0; i<27; i++ ){
		AxBu[i] = A[k++]*x[0] + B[i]*u[i];
	}	

	for( j=1; j<87; j++ ){		
		for( i=0; i<27; i++ ){
			AxBu[i] += A[k++]*x[j];
		}
	}

	for( i=0; i<27; i++ ){
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
 * where A is of size [27 x 87] and stored in column major format.
 * and B is of size [27 x 87] and stored in diagzero format
 * Note the transposes of A and B!
 */
void beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_27_87_27(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *y, beliefPenaltyMPC_FLOAT *z)
{
	int i;
	int j;
	int k = 0;
	for( i=0; i<27; i++ ){
		z[i] = 0;
		for( j=0; j<27; j++ ){
			z[i] += A[k++]*x[j];
		}
		z[i] += B[i]*y[i];
	}
	for( i=27 ;i<87; i++ ){
		z[i] = 0;
		for( j=0; j<27; j++ ){
			z[i] += A[k++]*x[j];
		}
	}
}


/*
 * Matrix vector multiplication y = M'*x where M is of size [27 x 33]
 * and stored in diagzero format. Note the transpose of M!
 */
void beliefPenaltyMPC_LA_DIAGZERO_MTVM_27_33(beliefPenaltyMPC_FLOAT *M, beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *y)
{
	int i;
	for( i=0; i<33; i++ ){
		y[i] = M[i]*x[i];
	}
}


/*
 * Vector subtraction and addition.
 *	 Input: five vectors t, tidx, u, v, w and two scalars z and r
 *	 Output: y = t(tidx) - u + w
 *           z = z - v'*x;
 *           r = max([norm(y,inf), z]);
 * for vectors of length 87. Output z is of course scalar.
 */
void beliefPenaltyMPC_LA_VSUBADD3_87(beliefPenaltyMPC_FLOAT* t, beliefPenaltyMPC_FLOAT* u, int* uidx, beliefPenaltyMPC_FLOAT* v, beliefPenaltyMPC_FLOAT* w, beliefPenaltyMPC_FLOAT* y, beliefPenaltyMPC_FLOAT* z, beliefPenaltyMPC_FLOAT* r)
{
	int i;
	beliefPenaltyMPC_FLOAT norm = *r;
	beliefPenaltyMPC_FLOAT vx = 0;
	beliefPenaltyMPC_FLOAT x;
	for( i=0; i<87; i++){
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
 * for vectors of length 33. Output z is of course scalar.
 */
void beliefPenaltyMPC_LA_VSUBADD2_33(beliefPenaltyMPC_FLOAT* t, int* tidx, beliefPenaltyMPC_FLOAT* u, beliefPenaltyMPC_FLOAT* v, beliefPenaltyMPC_FLOAT* w, beliefPenaltyMPC_FLOAT* y, beliefPenaltyMPC_FLOAT* z, beliefPenaltyMPC_FLOAT* r)
{
	int i;
	beliefPenaltyMPC_FLOAT norm = *r;
	beliefPenaltyMPC_FLOAT vx = 0;
	beliefPenaltyMPC_FLOAT x;
	for( i=0; i<33; i++){
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
 * for vectors of length 33. Output z is of course scalar.
 */
void beliefPenaltyMPC_LA_VSUBADD3_33(beliefPenaltyMPC_FLOAT* t, beliefPenaltyMPC_FLOAT* u, int* uidx, beliefPenaltyMPC_FLOAT* v, beliefPenaltyMPC_FLOAT* w, beliefPenaltyMPC_FLOAT* y, beliefPenaltyMPC_FLOAT* z, beliefPenaltyMPC_FLOAT* r)
{
	int i;
	beliefPenaltyMPC_FLOAT norm = *r;
	beliefPenaltyMPC_FLOAT vx = 0;
	beliefPenaltyMPC_FLOAT x;
	for( i=0; i<33; i++){
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
 * for vectors of length 27. Output z is of course scalar.
 */
void beliefPenaltyMPC_LA_VSUBADD2_27(beliefPenaltyMPC_FLOAT* t, int* tidx, beliefPenaltyMPC_FLOAT* u, beliefPenaltyMPC_FLOAT* v, beliefPenaltyMPC_FLOAT* w, beliefPenaltyMPC_FLOAT* y, beliefPenaltyMPC_FLOAT* z, beliefPenaltyMPC_FLOAT* r)
{
	int i;
	beliefPenaltyMPC_FLOAT norm = *r;
	beliefPenaltyMPC_FLOAT vx = 0;
	beliefPenaltyMPC_FLOAT x;
	for( i=0; i<27; i++){
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
 * Computes r = A*x - b + s
 * and      y = max([norm(r,inf), y])
 * and      z -= l'*(Ax-b)
 * where A is stored in column major format
 */
void beliefPenaltyMPC_LA_MVSUBADD_6_33(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *b, beliefPenaltyMPC_FLOAT *s, beliefPenaltyMPC_FLOAT *l, beliefPenaltyMPC_FLOAT *r, beliefPenaltyMPC_FLOAT *z, beliefPenaltyMPC_FLOAT *y)
{
	int i;
	int j;
	int k = 0;
	beliefPenaltyMPC_FLOAT Ax[6];
	beliefPenaltyMPC_FLOAT Axlessb;
	beliefPenaltyMPC_FLOAT norm = *y;
	beliefPenaltyMPC_FLOAT lAxlessb = 0;

	/* do A*x first */
	for( i=0; i<6; i++ ){
		Ax[i] = A[k++]*x[0];				
	}	
	for( j=1; j<33; j++ ){		
		for( i=0; i<6; i++ ){
			Ax[i] += A[k++]*x[j];
		}
	}

	for( i=0; i<6; i++ ){
		Axlessb = Ax[i] - b[i];
		r[i] = Axlessb + s[i];
		lAxlessb += l[i]*Axlessb;
		if( r[i] > norm ){
			norm = r[i];
		}
		if( -r[i] > norm ){
			norm = -r[i];
		}
	}
	*y = norm;
	*z -= lAxlessb;
}


/*
 * Computes inequality constraints gradient-
 * Special function for box constraints of length 87
 * Returns also L/S, a value that is often used elsewhere.
 */
void beliefPenaltyMPC_LA_INEQ_B_GRAD_87_87_33(beliefPenaltyMPC_FLOAT *lu, beliefPenaltyMPC_FLOAT *su, beliefPenaltyMPC_FLOAT *ru, beliefPenaltyMPC_FLOAT *ll, beliefPenaltyMPC_FLOAT *sl, beliefPenaltyMPC_FLOAT *rl, int* lbIdx, int* ubIdx, beliefPenaltyMPC_FLOAT *grad, beliefPenaltyMPC_FLOAT *lubysu, beliefPenaltyMPC_FLOAT *llbysl)
{
	int i;
	for( i=0; i<87; i++ ){
		grad[i] = 0;
	}
	for( i=0; i<87; i++ ){		
		llbysl[i] = ll[i] / sl[i];
		grad[lbIdx[i]] -= llbysl[i]*rl[i];
	}
	for( i=0; i<33; i++ ){
		lubysu[i] = lu[i] / su[i];
		grad[ubIdx[i]] += lubysu[i]*ru[i];
	}
}


/*
 * Computes inequality constraints gradient-
 * Special function for box constraints of length 33
 * Returns also L/S, a value that is often used elsewhere.
 */
void beliefPenaltyMPC_LA_INEQ_B_GRAD_33_33_27(beliefPenaltyMPC_FLOAT *lu, beliefPenaltyMPC_FLOAT *su, beliefPenaltyMPC_FLOAT *ru, beliefPenaltyMPC_FLOAT *ll, beliefPenaltyMPC_FLOAT *sl, beliefPenaltyMPC_FLOAT *rl, int* lbIdx, int* ubIdx, beliefPenaltyMPC_FLOAT *grad, beliefPenaltyMPC_FLOAT *lubysu, beliefPenaltyMPC_FLOAT *llbysl)
{
	int i;
	for( i=0; i<33; i++ ){
		grad[i] = 0;
	}
	for( i=0; i<33; i++ ){		
		llbysl[i] = ll[i] / sl[i];
		grad[lbIdx[i]] -= llbysl[i]*rl[i];
	}
	for( i=0; i<27; i++ ){
		lubysu[i] = lu[i] / su[i];
		grad[ubIdx[i]] += lubysu[i]*ru[i];
	}
}


/*
 * Special function for gradient of inequality constraints
 * Calculates grad += A'*(L/S)*rI
 */
void beliefPenaltyMPC_LA_INEQ_P_6_33(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *lp, beliefPenaltyMPC_FLOAT *sp, beliefPenaltyMPC_FLOAT *rip, beliefPenaltyMPC_FLOAT *grad, beliefPenaltyMPC_FLOAT *lpbysp)
{
	int i;
	int j;
	int k = 0;

	beliefPenaltyMPC_FLOAT lsr[6];
	
	/* do (L/S)*ri first */
	for( j=0; j<6; j++ ){
		lpbysp[j] = lp[j] / sp[j];
		lsr[j] = lpbysp[j]*rip[j];
	}

	for( i=0; i<33; i++ ){		
		for( j=0; j<6; j++ ){
			grad[i] += A[k++]*lsr[j];
		}
	}
}


/*
 * Addition of three vectors  z = u + w + v
 * of length 1251.
 */
void beliefPenaltyMPC_LA_VVADD3_1251(beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *v, beliefPenaltyMPC_FLOAT *w, beliefPenaltyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<1251; i++ ){
		z[i] = u[i] + v[i] + w[i];
	}
}


/*
 * Special function to compute the Dense positive definite 
 * augmented Hessian for block size 87.
 *
 * Inputs: - H = dense cost Hessian in column major storage format
 *         - llbysl = L / S of lower bounds
 *         - lubysu = L / S of upper bounds
 *
 * Output: Phi = H + diag(llbysl) + diag(lubysu)
 * where Phi is stored in lower triangular row major format
 */
void beliefPenaltyMPC_LA_INEQ_DENSE_HESS_87_87_33(beliefPenaltyMPC_FLOAT *H, beliefPenaltyMPC_FLOAT *llbysl, int* lbIdx, beliefPenaltyMPC_FLOAT *lubysu, int* ubIdx, beliefPenaltyMPC_FLOAT *Phi)
{
	int i;
	int j;
	int k = 0;
	
	/* copy lower triangular part of H into PHI */
	for( i=0; i<87; i++ ){
		for( j=0; j<=i; j++ ){
			Phi[k++] = H[i*87+j];
		}		
	}

	/* add llbysl onto Phi where necessary */
	for( i=0; i<87; i++ ){
		j = lbIdx[i];
		Phi[((j+1)*(j+2))/2-1] += llbysl[i];
	}

	/* add lubysu onto Phi where necessary */
	for( i=0; i<33; i++){
		j = ubIdx[i];
		Phi[((j+1)*(j+2))/2-1] +=  lubysu[i];
	}

}


/**
 * Cholesky factorization as above, but working on a matrix in 
 * lower triangular storage format of size 87.
 */
void beliefPenaltyMPC_LA_DENSE_CHOL2_87(beliefPenaltyMPC_FLOAT *A)
{
    int i, j, k, di, dj;
	 int ii, jj;
    beliefPenaltyMPC_FLOAT l;
    beliefPenaltyMPC_FLOAT Mii;
    
	ii=0; di=0;
    for( i=0; i<87; i++ ){
        l = 0;
        for( k=0; k<i; k++ ){
            l += A[ii+k]*A[ii+k];
        }        
        
        Mii = A[ii+i] - l;
        
#if beliefPenaltyMPC_SET_PRINTLEVEL > 0 && defined PRINTNUMERICALWARNINGS
        if( Mii < 1.0000000000000000E-013 ){
             PRINTTEXT("WARNING (CHOL2): small %d-th pivot in Cholesky fact. (=%3.1e < eps=%3.1e), regularizing to %3.1e\n",i,Mii,1.0000000000000000E-013,4.0000000000000002E-004);
			 A[ii+i] = 2.0000000000000000E-002;
		} else
		{
			A[ii+i] = sqrt(Mii);
		}
#else
		A[ii+i] = Mii < 1.0000000000000000E-013 ? 2.0000000000000000E-002 : sqrt(Mii);
#endif
                    
		jj = ((i+1)*(i+2))/2; dj = i+1;
        for( j=i+1; j<87; j++ ){
            l = 0;            
            for( k=0; k<i; k++ ){
                l += A[jj+k]*A[ii+k];
            }

			/* saturate values for numerical stability */
			l = MIN(l,  BIGMM);
			l = MAX(l, -BIGMM);

            A[jj+i] = (A[jj+i] - l)/A[ii+i];            
			jj += ++dj;
        }
		ii += ++di;
    }
}


/**
 * Forward substitution for the matrix equation A*L' = B
 * where A is to be computed and is of size [27 x 87],
 * B is given and of size [27 x 87], L is a lower tri-
 * angular matrix of size 87 stored in lower triangular 
 * storage format. Note the transpose of L!
 *
 * Result: A in column major storage format.
 *
 */
void beliefPenaltyMPC_LA_DENSE_MATRIXFORWARDSUB_27_87(beliefPenaltyMPC_FLOAT *L, beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *A)
{
    int i,j,k,di;
	 int ii;
    beliefPenaltyMPC_FLOAT a;
    
    ii=0; di=0;
    for( j=0; j<87; j++ ){        
        for( i=0; i<27; i++ ){
            a = B[j*27+i];
            for( k=0; k<j; k++ ){
                a -= A[k*27+i]*L[ii+k];
            }

			/* saturate for numerical stability */
			a = MIN(a, BIGM);
			a = MAX(a, -BIGM); 

            A[j*27+i] = a/L[ii+j];
        }
        ii += ++di;
    }
}


/**
 * Forward substitution for the matrix equation A*L' = B
 * where A is to be computed and is of size [27 x 87],
 * B is given and of size [27 x 87] stored in 
 * diagzero storage format, L is a lower tri-
 * angular matrix of size 87 stored in lower triangular 
 * storage format. Note the transpose of L!
 *
 * Result: A in column major storage format.
 *
 */
void beliefPenaltyMPC_LA_DENSE_DIAGZERO_MATRIXFORWARDSUB_27_87(beliefPenaltyMPC_FLOAT *L, beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *A)
{
    int i,j,k,di;
	 int ii;
    beliefPenaltyMPC_FLOAT a;
	
	/*
	* The matrix A has the form
	*
	* d u u u r r r r r 
	* 0 d u u r r r r r 
	* 0 0 d u r r r r r 
	* 0 0 0 d r r r r r
	*
	* |Part1|| Part 2 |
	* 
	* d: diagonal
	* u: upper
	* r: right
	*/
	
	
    /* Part 1 */
    ii=0; di=0;
    for( j=0; j<27; j++ ){        
        for( i=0; i<j; i++ ){
            /* Calculate part of A which is non-zero and not diagonal "u"
             * i < j */
            a = 0;
			
            for( k=i; k<j; k++ ){
                a -= A[k*27+i]*L[ii+k];
            }
            A[j*27+i] = a/L[ii+j];
        }
        /* do the diagonal "d"
         * i = j */
        A[j*27+j] = B[i]/L[ii+j];
        
        /* fill lower triangular part with zeros "0"
         * n > i > j */
        for( i=j+1     ; i < 27; i++ ){
            A[j*27+i] = 0;
        }
        
        /* increment index of L */
        ii += ++di;	
    }
	
	/* Part 2 */ 
	for( j=27; j<87; j++ ){        
        for( i=0; i<27; i++ ){
            /* Calculate part of A which is non-zero and not diagonal "r" */
            a = 0;
			
            for( k=i; k<j; k++ ){
                a -= A[k*27+i]*L[ii+k];
            }
            A[j*27+i] = a/L[ii+j];
        }
        
        /* increment index of L */
        ii += ++di;	
    }
	
	
	
}


/**
 * Compute C = A*B' where 
 *
 *	size(A) = [27 x 87]
 *  size(B) = [27 x 87]
 * 
 * and all matrices are stored in column major format.
 *
 * THIS ONE HAS THE WORST ACCES PATTERN POSSIBLE.  
 * 
 */
void beliefPenaltyMPC_LA_DENSE_MMTM_27_87_27(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *C)
{
    int i, j, k;
    beliefPenaltyMPC_FLOAT temp;
    
    for( i=0; i<27; i++ ){        
        for( j=0; j<27; j++ ){
            temp = 0; 
            for( k=0; k<87; k++ ){
                temp += A[k*27+i]*B[k*27+j];
            }						
            C[j*27+i] = temp;
        }
    }
}


/**
 * Forward substitution to solve L*y = b where L is a
 * lower triangular matrix in triangular storage format.
 * 
 * The dimensions involved are 87.
 */
void beliefPenaltyMPC_LA_DENSE_FORWARDSUB_87(beliefPenaltyMPC_FLOAT *L, beliefPenaltyMPC_FLOAT *b, beliefPenaltyMPC_FLOAT *y)
{
    int i,j,ii,di;
    beliefPenaltyMPC_FLOAT yel;
            
    ii = 0; di = 0;
    for( i=0; i<87; i++ ){
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


/*
 * Special function to compute the Dense positive definite 
 * augmented Hessian for block size 33.
 *
 * Inputs: - H = dense cost Hessian in column major storage format
 *         - llbysl = L / S of lower bounds
 *         - lubysu = L / S of upper bounds
 *
 * Output: Phi = H + diag(llbysl) + diag(lubysu)
 * where Phi is stored in lower triangular row major format
 */
void beliefPenaltyMPC_LA_INEQ_DENSE_HESS_33_33_27(beliefPenaltyMPC_FLOAT *H, beliefPenaltyMPC_FLOAT *llbysl, int* lbIdx, beliefPenaltyMPC_FLOAT *lubysu, int* ubIdx, beliefPenaltyMPC_FLOAT *Phi)
{
	int i;
	int j;
	int k = 0;
	
	/* copy lower triangular part of H into PHI */
	for( i=0; i<33; i++ ){
		for( j=0; j<=i; j++ ){
			Phi[k++] = H[i*33+j];
		}		
	}

	/* add llbysl onto Phi where necessary */
	for( i=0; i<33; i++ ){
		j = lbIdx[i];
		Phi[((j+1)*(j+2))/2-1] += llbysl[i];
	}

	/* add lubysu onto Phi where necessary */
	for( i=0; i<27; i++){
		j = ubIdx[i];
		Phi[((j+1)*(j+2))/2-1] +=  lubysu[i];
	}

}


/**
 * Compute X = X + A'*D*A, where A is a general full matrix, D is
 * is a diagonal matrix stored in the vector d and X is a symmetric
 * positive definite matrix in lower triangular storage format. 
 * A is stored in column major format and is of size [6 x 33]
 * Phi is of size [33 x 33].
 */
void beliefPenaltyMPC_LA_DENSE_ADDMTDM_6_33(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *d, beliefPenaltyMPC_FLOAT *X)
{    
    int i,j,k,ii,di;
    beliefPenaltyMPC_FLOAT x;
    
    di = 0; ii = 0;
    for( i=0; i<33; i++ ){        
        for( j=0; j<=i; j++ ){
            x = 0;
            for( k=0; k<6; k++ ){
                x += A[i*6+k]*A[j*6+k]*d[k];
            }
            X[ii+j] += x;
        }
        ii += ++di;
    }
}


/**
 * Cholesky factorization as above, but working on a matrix in 
 * lower triangular storage format of size 33.
 */
void beliefPenaltyMPC_LA_DENSE_CHOL2_33(beliefPenaltyMPC_FLOAT *A)
{
    int i, j, k, di, dj;
	 int ii, jj;
    beliefPenaltyMPC_FLOAT l;
    beliefPenaltyMPC_FLOAT Mii;
    
	ii=0; di=0;
    for( i=0; i<33; i++ ){
        l = 0;
        for( k=0; k<i; k++ ){
            l += A[ii+k]*A[ii+k];
        }        
        
        Mii = A[ii+i] - l;
        
#if beliefPenaltyMPC_SET_PRINTLEVEL > 0 && defined PRINTNUMERICALWARNINGS
        if( Mii < 1.0000000000000000E-013 ){
             PRINTTEXT("WARNING (CHOL2): small %d-th pivot in Cholesky fact. (=%3.1e < eps=%3.1e), regularizing to %3.1e\n",i,Mii,1.0000000000000000E-013,4.0000000000000002E-004);
			 A[ii+i] = 2.0000000000000000E-002;
		} else
		{
			A[ii+i] = sqrt(Mii);
		}
#else
		A[ii+i] = Mii < 1.0000000000000000E-013 ? 2.0000000000000000E-002 : sqrt(Mii);
#endif
                    
		jj = ((i+1)*(i+2))/2; dj = i+1;
        for( j=i+1; j<33; j++ ){
            l = 0;            
            for( k=0; k<i; k++ ){
                l += A[jj+k]*A[ii+k];
            }

			/* saturate values for numerical stability */
			l = MIN(l,  BIGMM);
			l = MAX(l, -BIGMM);

            A[jj+i] = (A[jj+i] - l)/A[ii+i];            
			jj += ++dj;
        }
		ii += ++di;
    }
}


/**
 * Forward substitution for the matrix equation A*L' = B
 * where A is to be computed and is of size [27 x 33],
 * B is given and of size [27 x 33] stored in 
 * diagzero storage format, L is a lower tri-
 * angular matrix of size 33 stored in lower triangular 
 * storage format. Note the transpose of L!
 *
 * Result: A in column major storage format.
 *
 */
void beliefPenaltyMPC_LA_DENSE_DIAGZERO_MATRIXFORWARDSUB_27_33(beliefPenaltyMPC_FLOAT *L, beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *A)
{
    int i,j,k,di;
	 int ii;
    beliefPenaltyMPC_FLOAT a;
	
	/*
	* The matrix A has the form
	*
	* d u u u r r r r r 
	* 0 d u u r r r r r 
	* 0 0 d u r r r r r 
	* 0 0 0 d r r r r r
	*
	* |Part1|| Part 2 |
	* 
	* d: diagonal
	* u: upper
	* r: right
	*/
	
	
    /* Part 1 */
    ii=0; di=0;
    for( j=0; j<27; j++ ){        
        for( i=0; i<j; i++ ){
            /* Calculate part of A which is non-zero and not diagonal "u"
             * i < j */
            a = 0;
			
            for( k=i; k<j; k++ ){
                a -= A[k*27+i]*L[ii+k];
            }
            A[j*27+i] = a/L[ii+j];
        }
        /* do the diagonal "d"
         * i = j */
        A[j*27+j] = B[i]/L[ii+j];
        
        /* fill lower triangular part with zeros "0"
         * n > i > j */
        for( i=j+1     ; i < 27; i++ ){
            A[j*27+i] = 0;
        }
        
        /* increment index of L */
        ii += ++di;	
    }
	
	/* Part 2 */ 
	for( j=27; j<33; j++ ){        
        for( i=0; i<27; i++ ){
            /* Calculate part of A which is non-zero and not diagonal "r" */
            a = 0;
			
            for( k=i; k<j; k++ ){
                a -= A[k*27+i]*L[ii+k];
            }
            A[j*27+i] = a/L[ii+j];
        }
        
        /* increment index of L */
        ii += ++di;	
    }
	
	
	
}


/**
 * Forward substitution to solve L*y = b where L is a
 * lower triangular matrix in triangular storage format.
 * 
 * The dimensions involved are 33.
 */
void beliefPenaltyMPC_LA_DENSE_FORWARDSUB_33(beliefPenaltyMPC_FLOAT *L, beliefPenaltyMPC_FLOAT *b, beliefPenaltyMPC_FLOAT *y)
{
    int i,j,ii,di;
    beliefPenaltyMPC_FLOAT yel;
            
    ii = 0; di = 0;
    for( i=0; i<33; i++ ){
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
 * Compute L = B*B', where L is lower triangular of size NXp1
 * and B is a diagzero matrix of size [27 x 87] in column
 * storage format.
 * 
 * THIS ONE HAS THE WORST ACCES PATTERN POSSIBLE. 
 * POSSIBKE FIX: PUT A AND B INTO ROW MAJOR FORMAT FIRST.
 * 
 */
void beliefPenaltyMPC_LA_DENSE_MMT_27_87(beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *L)
{
    int i, j, k, ii, di;
    beliefPenaltyMPC_FLOAT ltemp;
    
    ii = 0; di = 0;
    for( i=0; i<27; i++ ){        
        for( j=0; j<=i; j++ ){
            ltemp = 0; 		
			for( k=0; k<87; k++ ){
                ltemp += B[k*27+i]*B[k*27+j];
            }
            L[ii+j] = ltemp;
        }
        ii += ++di;
    }
}


/* 
 * Computes r = b - B*u
 * where B is stored in column major format
 */
void beliefPenaltyMPC_LA_DENSE_MVMSUB7_27_87(beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *b, beliefPenaltyMPC_FLOAT *r)
{
	int i;
	int m = 0;
	int n;

	for( i=0; i<27; i++ ){
		r[i] = b[i] - B[m++]*u[0];
	}	
	
	for( n=1; n<87; n++ ){
		for( i=0; i<27; i++ ){
			r[i] -= B[m++]*u[n];
		}		
	}
}


/**
 * Compute L = A*A' + B*B', where L is lower triangular of size NXp1
 * and A is a dense matrix of size [27 x 87] in column
 * storage format, and B is of size [27 x 87] also in column
 * storage format.
 * 
 * THIS ONE HAS THE WORST ACCES PATTERN POSSIBLE. 
 * POSSIBKE FIX: PUT A AND B INTO ROW MAJOR FORMAT FIRST.
 * 
 */
void beliefPenaltyMPC_LA_DENSE_MMT2_27_87_87(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *L)
{
    int i, j, k, ii, di;
    beliefPenaltyMPC_FLOAT ltemp;
    
    ii = 0; di = 0;
    for( i=0; i<27; i++ ){        
        for( j=0; j<=i; j++ ){
            ltemp = 0; 
            for( k=0; k<87; k++ ){
                ltemp += A[k*27+i]*A[k*27+j];
            }			
			for( k=0; k<87; k++ ){
                ltemp += B[k*27+i]*B[k*27+j];
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
void beliefPenaltyMPC_LA_DENSE_MVMSUB2_27_87_87(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *b, beliefPenaltyMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;
	int m = 0;
	int n;

	for( i=0; i<27; i++ ){
		r[i] = b[i] - A[k++]*x[0] - B[m++]*u[0];
	}	
	for( j=1; j<87; j++ ){		
		for( i=0; i<27; i++ ){
			r[i] -= A[k++]*x[j];
		}
	}
	
	for( n=1; n<87; n++ ){
		for( i=0; i<27; i++ ){
			r[i] -= B[m++]*u[n];
		}		
	}
}


/**
 * Compute L = A*A' + B*B', where L is lower triangular of size NXp1
 * and A is a dense matrix of size [27 x 87] in column
 * storage format, and B is of size [27 x 33] also in column
 * storage format.
 * 
 * THIS ONE HAS THE WORST ACCES PATTERN POSSIBLE. 
 * POSSIBKE FIX: PUT A AND B INTO ROW MAJOR FORMAT FIRST.
 * 
 */
void beliefPenaltyMPC_LA_DENSE_MMT2_27_87_33(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *L)
{
    int i, j, k, ii, di;
    beliefPenaltyMPC_FLOAT ltemp;
    
    ii = 0; di = 0;
    for( i=0; i<27; i++ ){        
        for( j=0; j<=i; j++ ){
            ltemp = 0; 
            for( k=0; k<87; k++ ){
                ltemp += A[k*27+i]*A[k*27+j];
            }			
			for( k=0; k<33; k++ ){
                ltemp += B[k*27+i]*B[k*27+j];
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
void beliefPenaltyMPC_LA_DENSE_MVMSUB2_27_87_33(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *b, beliefPenaltyMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;
	int m = 0;
	int n;

	for( i=0; i<27; i++ ){
		r[i] = b[i] - A[k++]*x[0] - B[m++]*u[0];
	}	
	for( j=1; j<87; j++ ){		
		for( i=0; i<27; i++ ){
			r[i] -= A[k++]*x[j];
		}
	}
	
	for( n=1; n<33; n++ ){
		for( i=0; i<27; i++ ){
			r[i] -= B[m++]*u[n];
		}		
	}
}


/**
 * Cholesky factorization as above, but working on a matrix in 
 * lower triangular storage format of size 27 and outputting
 * the Cholesky factor to matrix L in lower triangular format.
 */
void beliefPenaltyMPC_LA_DENSE_CHOL_27(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *L)
{
    int i, j, k, di, dj;
	 int ii, jj;

    beliefPenaltyMPC_FLOAT l;
    beliefPenaltyMPC_FLOAT Mii;

	/* copy A to L first and then operate on L */
	/* COULD BE OPTIMIZED */
	ii=0; di=0;
	for( i=0; i<27; i++ ){
		for( j=0; j<=i; j++ ){
			L[ii+j] = A[ii+j];
		}
		ii += ++di;
	}    
	
	/* factor L */
	ii=0; di=0;
    for( i=0; i<27; i++ ){
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
        for( j=i+1; j<27; j++ ){
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
 * The dimensions involved are 27.
 */
void beliefPenaltyMPC_LA_DENSE_FORWARDSUB_27(beliefPenaltyMPC_FLOAT *L, beliefPenaltyMPC_FLOAT *b, beliefPenaltyMPC_FLOAT *y)
{
    int i,j,ii,di;
    beliefPenaltyMPC_FLOAT yel;
            
    ii = 0; di = 0;
    for( i=0; i<27; i++ ){
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
 * where A is to be computed and is of size [27 x 27],
 * B is given and of size [27 x 27], L is a lower tri-
 * angular matrix of size 27 stored in lower triangular 
 * storage format. Note the transpose of L AND B!
 *
 * Result: A in column major storage format.
 *
 */
void beliefPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_27_27(beliefPenaltyMPC_FLOAT *L, beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *A)
{
    int i,j,k,ii,di;
    beliefPenaltyMPC_FLOAT a;
    
    ii=0; di=0;
    for( j=0; j<27; j++ ){        
        for( i=0; i<27; i++ ){
            a = B[i*27+j];
            for( k=0; k<j; k++ ){
                a -= A[k*27+i]*L[ii+k];
            }    

			/* saturate for numerical stability */
			a = MIN(a, BIGM);
			a = MAX(a, -BIGM); 

			A[j*27+i] = a/L[ii+j];			
        }
        ii += ++di;
    }
}


/**
 * Compute L = L - A*A', where L is lower triangular of size 27
 * and A is a dense matrix of size [27 x 27] in column
 * storage format.
 * 
 * THIS ONE HAS THE WORST ACCES PATTERN POSSIBLE. 
 * POSSIBKE FIX: PUT A INTO ROW MAJOR FORMAT FIRST.
 * 
 */
void beliefPenaltyMPC_LA_DENSE_MMTSUB_27_27(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *L)
{
    int i, j, k, ii, di;
    beliefPenaltyMPC_FLOAT ltemp;
    
    ii = 0; di = 0;
    for( i=0; i<27; i++ ){        
        for( j=0; j<=i; j++ ){
            ltemp = 0; 
            for( k=0; k<27; k++ ){
                ltemp += A[k*27+i]*A[k*27+j];
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
void beliefPenaltyMPC_LA_DENSE_MVMSUB1_27_27(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *b, beliefPenaltyMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;

	for( i=0; i<27; i++ ){
		r[i] = b[i] - A[k++]*x[0];
	}	
	for( j=1; j<27; j++ ){		
		for( i=0; i<27; i++ ){
			r[i] -= A[k++]*x[j];
		}
	}
}


/**
 * Backward Substitution to solve L^T*x = y where L is a
 * lower triangular matrix in triangular storage format.
 * 
 * All involved dimensions are 27.
 */
void beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_27(beliefPenaltyMPC_FLOAT *L, beliefPenaltyMPC_FLOAT *y, beliefPenaltyMPC_FLOAT *x)
{
    int i, ii, di, j, jj, dj;
    beliefPenaltyMPC_FLOAT xel;    
	int start = 351;
    
    /* now solve L^T*x = y by backward substitution */
    ii = start; di = 26;
    for( i=26; i>=0; i-- ){        
        xel = y[i];        
        jj = start; dj = 26;
        for( j=26; j>i; j-- ){
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
 * Matrix vector multiplication y = b - M'*x where M is of size [27 x 27]
 * and stored in column major format. Note the transpose of M!
 */
void beliefPenaltyMPC_LA_DENSE_MTVMSUB_27_27(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *b, beliefPenaltyMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0; 
	for( i=0; i<27; i++ ){
		r[i] = b[i];
		for( j=0; j<27; j++ ){
			r[i] -= A[k++]*x[j];
		}
	}
}


/*
 * Vector subtraction z = -x - y for vectors of length 1251.
 */
void beliefPenaltyMPC_LA_VSUB2_1251(beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *y, beliefPenaltyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<1251; i++){
		z[i] = -x[i] - y[i];
	}
}


/**
 * Forward-Backward-Substitution to solve L*L^T*x = b where L is a
 * lower triangular matrix of size 87 in lower triangular
 * storage format.
 */
void beliefPenaltyMPC_LA_DENSE_FORWARDBACKWARDSUB_87(beliefPenaltyMPC_FLOAT *L, beliefPenaltyMPC_FLOAT *b, beliefPenaltyMPC_FLOAT *x)
{
    int i, ii, di, j, jj, dj;
    beliefPenaltyMPC_FLOAT y[87];
    beliefPenaltyMPC_FLOAT yel,xel;
	int start = 3741;
            
    /* first solve Ly = b by forward substitution */
     ii = 0; di = 0;
    for( i=0; i<87; i++ ){
        yel = b[i];        
        for( j=0; j<i; j++ ){
            yel -= y[j]*L[ii+j];
        }

		/* saturate for numerical stability */
		yel = MIN(yel, BIGM);
		yel = MAX(yel, -BIGM); 

        y[i] = yel / L[ii+i];
        ii += ++di;
    }
    
    /* now solve L^T*x = y by backward substitution */
    ii = start; di = 86;
    for( i=86; i>=0; i-- ){        
        xel = y[i];        
        jj = start; dj = 86;
        for( j=86; j>i; j-- ){
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


/**
 * Forward-Backward-Substitution to solve L*L^T*x = b where L is a
 * lower triangular matrix of size 33 in lower triangular
 * storage format.
 */
void beliefPenaltyMPC_LA_DENSE_FORWARDBACKWARDSUB_33(beliefPenaltyMPC_FLOAT *L, beliefPenaltyMPC_FLOAT *b, beliefPenaltyMPC_FLOAT *x)
{
    int i, ii, di, j, jj, dj;
    beliefPenaltyMPC_FLOAT y[33];
    beliefPenaltyMPC_FLOAT yel,xel;
	int start = 528;
            
    /* first solve Ly = b by forward substitution */
     ii = 0; di = 0;
    for( i=0; i<33; i++ ){
        yel = b[i];        
        for( j=0; j<i; j++ ){
            yel -= y[j]*L[ii+j];
        }

		/* saturate for numerical stability */
		yel = MIN(yel, BIGM);
		yel = MAX(yel, -BIGM); 

        y[i] = yel / L[ii+i];
        ii += ++di;
    }
    
    /* now solve L^T*x = y by backward substitution */
    ii = start; di = 32;
    for( i=32; i>=0; i-- ){        
        xel = y[i];        
        jj = start; dj = 32;
        for( j=32; j>i; j-- ){
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
 * Vector subtraction z = x(xidx) - y where y, z and xidx are of length 87,
 * and x has length 87 and is indexed through yidx.
 */
void beliefPenaltyMPC_LA_VSUB_INDEXED_87(beliefPenaltyMPC_FLOAT *x, int* xidx, beliefPenaltyMPC_FLOAT *y, beliefPenaltyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<87; i++){
		z[i] = x[xidx[i]] - y[i];
	}
}


/*
 * Vector subtraction x = -u.*v - w for vectors of length 87.
 */
void beliefPenaltyMPC_LA_VSUB3_87(beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *v, beliefPenaltyMPC_FLOAT *w, beliefPenaltyMPC_FLOAT *x)
{
	int i;
	for( i=0; i<87; i++){
		x[i] = -u[i]*v[i] - w[i];
	}
}


/*
 * Vector subtraction z = -x - y(yidx) where y is of length 87
 * and z, x and yidx are of length 33.
 */
void beliefPenaltyMPC_LA_VSUB2_INDEXED_33(beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *y, int* yidx, beliefPenaltyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<33; i++){
		z[i] = -x[i] - y[yidx[i]];
	}
}


/*
 * Vector subtraction x = -u.*v - w for vectors of length 33.
 */
void beliefPenaltyMPC_LA_VSUB3_33(beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *v, beliefPenaltyMPC_FLOAT *w, beliefPenaltyMPC_FLOAT *x)
{
	int i;
	for( i=0; i<33; i++){
		x[i] = -u[i]*v[i] - w[i];
	}
}


/*
 * Vector subtraction z = x(xidx) - y where y, z and xidx are of length 33,
 * and x has length 33 and is indexed through yidx.
 */
void beliefPenaltyMPC_LA_VSUB_INDEXED_33(beliefPenaltyMPC_FLOAT *x, int* xidx, beliefPenaltyMPC_FLOAT *y, beliefPenaltyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<33; i++){
		z[i] = x[xidx[i]] - y[i];
	}
}


/*
 * Vector subtraction z = -x - y(yidx) where y is of length 33
 * and z, x and yidx are of length 27.
 */
void beliefPenaltyMPC_LA_VSUB2_INDEXED_27(beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *y, int* yidx, beliefPenaltyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<27; i++){
		z[i] = -x[i] - y[yidx[i]];
	}
}


/*
 * Vector subtraction x = -u.*v - w for vectors of length 27.
 */
void beliefPenaltyMPC_LA_VSUB3_27(beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *v, beliefPenaltyMPC_FLOAT *w, beliefPenaltyMPC_FLOAT *x)
{
	int i;
	for( i=0; i<27; i++){
		x[i] = -u[i]*v[i] - w[i];
	}
}


/* 
 * Computes r = -b - A*x
 * where A is stored in column major format
 */
void beliefPenaltyMPC_LA_DENSE_MVMSUB4_6_33(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *b, beliefPenaltyMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;

	for( i=0; i<6; i++ ){
		r[i] = -b[i] - A[k++]*x[0];
	}	
	for( j=1; j<33; j++ ){		
		for( i=0; i<6; i++ ){
			r[i] -= A[k++]*x[j];
		}
	}
}


/*
 * Vector subtraction x = -u.*v - w for vectors of length 6.
 */
void beliefPenaltyMPC_LA_VSUB3_6(beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *v, beliefPenaltyMPC_FLOAT *w, beliefPenaltyMPC_FLOAT *x)
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
        for( i=0; i<1746; i++ ){
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
        if( i == 1746 ){
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
    *mu_aff = mymu / (beliefPenaltyMPC_FLOAT)1746;
    return lsIt;
}


/*
 * Vector subtraction x = u.*v - a where a is a scalar
*  and x,u,v are vectors of length 1746.
 */
void beliefPenaltyMPC_LA_VSUB5_1746(beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *v, beliefPenaltyMPC_FLOAT a, beliefPenaltyMPC_FLOAT *x)
{
	int i;
	for( i=0; i<1746; i++){
		x[i] = u[i]*v[i] - a;
	}
}


/*
 * Computes x=0; x(uidx) += u/su; x(vidx) -= v/sv where x is of length 87,
 * u, su, uidx are of length 33 and v, sv, vidx are of length 87.
 */
void beliefPenaltyMPC_LA_VSUB6_INDEXED_87_33_87(beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *su, int* uidx, beliefPenaltyMPC_FLOAT *v, beliefPenaltyMPC_FLOAT *sv, int* vidx, beliefPenaltyMPC_FLOAT *x)
{
	int i;
	for( i=0; i<87; i++ ){
		x[i] = 0;
	}
	for( i=0; i<33; i++){
		x[uidx[i]] += u[i]/su[i];
	}
	for( i=0; i<87; i++){
		x[vidx[i]] -= v[i]/sv[i];
	}
}


/* 
 * Computes r = B*u
 * where B is stored in column major format
 */
void beliefPenaltyMPC_LA_DENSE_MVM_27_87(beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *r)
{
	int i;
	int m = 0;
	int n;

	for( i=0; i<27; i++ ){
		r[i] = B[m++]*u[0];
	}	
	
	for( n=1; n<87; n++ ){
		for( i=0; i<27; i++ ){
			r[i] += B[m++]*u[n];
		}		
	}
}


/* 
 * Computes r = A*x + B*u
 * where A an B are stored in column major format
 */
void beliefPenaltyMPC_LA_DENSE_2MVMADD_27_87_87(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;
	int m = 0;
	int n;

	for( i=0; i<27; i++ ){
		r[i] = A[k++]*x[0] + B[m++]*u[0];
	}	

	for( j=1; j<87; j++ ){		
		for( i=0; i<27; i++ ){
			r[i] += A[k++]*x[j];
		}
	}
	
	for( n=1; n<87; n++ ){
		for( i=0; i<27; i++ ){
			r[i] += B[m++]*u[n];
		}		
	}
}


/*
 * Computes x=0; x(uidx) += u/su; x(vidx) -= v/sv where x is of length 33,
 * u, su, uidx are of length 27 and v, sv, vidx are of length 33.
 */
void beliefPenaltyMPC_LA_VSUB6_INDEXED_33_27_33(beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *su, int* uidx, beliefPenaltyMPC_FLOAT *v, beliefPenaltyMPC_FLOAT *sv, int* vidx, beliefPenaltyMPC_FLOAT *x)
{
	int i;
	for( i=0; i<33; i++ ){
		x[i] = 0;
	}
	for( i=0; i<27; i++){
		x[uidx[i]] += u[i]/su[i];
	}
	for( i=0; i<33; i++){
		x[vidx[i]] -= v[i]/sv[i];
	}
}


/*
 * Matrix vector multiplication z = z + A'*(x./s) where A is of size [6 x 33]
 * and stored in column major format. Note the transpose of M!
 */
void beliefPenaltyMPC_LA_DENSE_MTVMADD2_6_33(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *s, beliefPenaltyMPC_FLOAT *z)
{
	int i;
	int j;
	int k = 0; 
	beliefPenaltyMPC_FLOAT temp[6];

	for( j=0; j<6; j++ ){
		temp[j] = x[j] / s[j];
	}

	for( i=0; i<33; i++ ){
		for( j=0; j<6; j++ ){
			z[i] += A[k++]*temp[j];
		}
	}
}


/* 
 * Computes r = A*x + B*u
 * where A an B are stored in column major format
 */
void beliefPenaltyMPC_LA_DENSE_2MVMADD_27_87_33(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;
	int m = 0;
	int n;

	for( i=0; i<27; i++ ){
		r[i] = A[k++]*x[0] + B[m++]*u[0];
	}	

	for( j=1; j<87; j++ ){		
		for( i=0; i<27; i++ ){
			r[i] += A[k++]*x[j];
		}
	}
	
	for( n=1; n<33; n++ ){
		for( i=0; i<27; i++ ){
			r[i] += B[m++]*u[n];
		}		
	}
}


/*
 * Vector subtraction z = x - y for vectors of length 1251.
 */
void beliefPenaltyMPC_LA_VSUB_1251(beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *y, beliefPenaltyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<1251; i++){
		z[i] = x[i] - y[i];
	}
}


/** 
 * Computes z = -r./s - u.*y(y)
 * where all vectors except of y are of length 87 (length of y >= 87).
 */
void beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_87(beliefPenaltyMPC_FLOAT *r, beliefPenaltyMPC_FLOAT *s, beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *y, int* yidx, beliefPenaltyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<87; i++ ){
		z[i] = -r[i]/s[i] - u[i]*y[yidx[i]];
	}
}


/** 
 * Computes z = -r./s + u.*y(y)
 * where all vectors except of y are of length 33 (length of y >= 33).
 */
void beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_33(beliefPenaltyMPC_FLOAT *r, beliefPenaltyMPC_FLOAT *s, beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *y, int* yidx, beliefPenaltyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<33; i++ ){
		z[i] = -r[i]/s[i] + u[i]*y[yidx[i]];
	}
}


/** 
 * Computes z = -r./s - u.*y(y)
 * where all vectors except of y are of length 33 (length of y >= 33).
 */
void beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_33(beliefPenaltyMPC_FLOAT *r, beliefPenaltyMPC_FLOAT *s, beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *y, int* yidx, beliefPenaltyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<33; i++ ){
		z[i] = -r[i]/s[i] - u[i]*y[yidx[i]];
	}
}


/** 
 * Computes z = -r./s + u.*y(y)
 * where all vectors except of y are of length 27 (length of y >= 27).
 */
void beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_27(beliefPenaltyMPC_FLOAT *r, beliefPenaltyMPC_FLOAT *s, beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *y, int* yidx, beliefPenaltyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<27; i++ ){
		z[i] = -r[i]/s[i] + u[i]*y[yidx[i]];
	}
}


/* 
 * Computes r = (-b + l.*(A*x))./s
 * where A is stored in column major format
 */
void beliefPenaltyMPC_LA_DENSE_MVMSUB5_6_33(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *b, beliefPenaltyMPC_FLOAT *s, beliefPenaltyMPC_FLOAT *l, beliefPenaltyMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;

	beliefPenaltyMPC_FLOAT temp[6];

	
	for( i=0; i<6; i++ ){
		temp[i] = A[k++]*x[0];
	}
	

	for( j=1; j<33; j++ ){		
		for( i=0; i<6; i++ ){
			temp[i] += A[k++]*x[j];
		}
	}

	for( i=0; i<6; i++ ){
		r[i] = (-b[i] + l[i]*temp[i])/s[i]; 
	}	
	
}


/*
 * Computes ds = -l.\(r + s.*dl) for vectors of length 1746.
 */
void beliefPenaltyMPC_LA_VSUB7_1746(beliefPenaltyMPC_FLOAT *l, beliefPenaltyMPC_FLOAT *r, beliefPenaltyMPC_FLOAT *s, beliefPenaltyMPC_FLOAT *dl, beliefPenaltyMPC_FLOAT *ds)
{
	int i;
	for( i=0; i<1746; i++){
		ds[i] = -(r[i] + s[i]*dl[i])/l[i];
	}
}


/*
 * Vector addition x = x + y for vectors of length 1251.
 */
void beliefPenaltyMPC_LA_VADD_1251(beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *y)
{
	int i;
	for( i=0; i<1251; i++){
		x[i] += y[i];
	}
}


/*
 * Vector addition x = x + y for vectors of length 405.
 */
void beliefPenaltyMPC_LA_VADD_405(beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *y)
{
	int i;
	for( i=0; i<405; i++){
		x[i] += y[i];
	}
}


/*
 * Vector addition x = x + y for vectors of length 1746.
 */
void beliefPenaltyMPC_LA_VADD_1746(beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *y)
{
	int i;
	for( i=0; i<1746; i++){
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
        for( i=0; i<1746; i++ ){
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
        if( i == 1746 ){
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
    for( i=0; i<1251; i++ ){
        z[i] += a_gamma*dz[i];
    }
    
    /* equality constraint multipliers */
    for( i=0; i<405; i++ ){
        v[i] += a_gamma*dv[i];
    }
    
    /* inequality constraint multipliers & slacks, also update mu */
    *mu = 0;
    for( i=0; i<1746; i++ ){
        dltemp = l[i] + a_gamma*dl[i]; l[i] = dltemp;
        dstemp = s[i] + a_gamma*ds[i]; s[i] = dstemp;
        *mu += dltemp*dstemp;
    }
    
    *a = a_gamma;
    *mu /= (beliefPenaltyMPC_FLOAT)1746;
    return lsIt;
}




/* VARIABLE DEFINITIONS ------------------------------------------------ */
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_z[1251];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_v[405];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_dz_aff[1251];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_dv_aff[405];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_grad_cost[1251];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_grad_eq[1251];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_rd[1251];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_l[1746];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_s[1746];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_lbys[1746];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_dl_aff[1746];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ds_aff[1746];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_dz_cc[1251];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_dv_cc[405];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_dl_cc[1746];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ds_cc[1746];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ccrhs[1746];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_grad_ineq[1251];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_z00 = beliefPenaltyMPC_z + 0;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzaff00 = beliefPenaltyMPC_dz_aff + 0;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzcc00 = beliefPenaltyMPC_dz_cc + 0;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_rd00 = beliefPenaltyMPC_rd + 0;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lbyrd00[87];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_cost00 = beliefPenaltyMPC_grad_cost + 0;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_eq00 = beliefPenaltyMPC_grad_eq + 0;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_ineq00 = beliefPenaltyMPC_grad_ineq + 0;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ctv00[87];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_v00 = beliefPenaltyMPC_v + 0;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_re00[27];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_beta00[27];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_betacc00[27];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvaff00 = beliefPenaltyMPC_dv_aff + 0;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvcc00 = beliefPenaltyMPC_dv_cc + 0;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_V00[2349];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Yd00[378];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ld00[378];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_yy00[27];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_bmy00[27];
int beliefPenaltyMPC_lbIdx00[87] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llb00 = beliefPenaltyMPC_l + 0;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_slb00 = beliefPenaltyMPC_s + 0;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llbbyslb00 = beliefPenaltyMPC_lbys + 0;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_rilb00[87];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbaff00 = beliefPenaltyMPC_dl_aff + 0;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbaff00 = beliefPenaltyMPC_ds_aff + 0;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbcc00 = beliefPenaltyMPC_dl_cc + 0;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbcc00 = beliefPenaltyMPC_ds_cc + 0;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsl00 = beliefPenaltyMPC_ccrhs + 0;
int beliefPenaltyMPC_ubIdx00[33] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lub00 = beliefPenaltyMPC_l + 87;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_sub00 = beliefPenaltyMPC_s + 87;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lubbysub00 = beliefPenaltyMPC_lbys + 87;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_riub00[33];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubaff00 = beliefPenaltyMPC_dl_aff + 87;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubaff00 = beliefPenaltyMPC_ds_aff + 87;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubcc00 = beliefPenaltyMPC_dl_cc + 87;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubcc00 = beliefPenaltyMPC_ds_cc + 87;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsub00 = beliefPenaltyMPC_ccrhs + 87;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Phi00[3828];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_D00[87] = {1.0000000000000000E+000, 
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
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_W00[2349];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_z01 = beliefPenaltyMPC_z + 87;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzaff01 = beliefPenaltyMPC_dz_aff + 87;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzcc01 = beliefPenaltyMPC_dz_cc + 87;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_rd01 = beliefPenaltyMPC_rd + 87;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lbyrd01[87];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_cost01 = beliefPenaltyMPC_grad_cost + 87;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_eq01 = beliefPenaltyMPC_grad_eq + 87;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_ineq01 = beliefPenaltyMPC_grad_ineq + 87;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ctv01[87];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_v01 = beliefPenaltyMPC_v + 27;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_re01[27];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_beta01[27];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_betacc01[27];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvaff01 = beliefPenaltyMPC_dv_aff + 27;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvcc01 = beliefPenaltyMPC_dv_cc + 27;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_V01[2349];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Yd01[378];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ld01[378];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_yy01[27];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_bmy01[27];
int beliefPenaltyMPC_lbIdx01[87] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llb01 = beliefPenaltyMPC_l + 120;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_slb01 = beliefPenaltyMPC_s + 120;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llbbyslb01 = beliefPenaltyMPC_lbys + 120;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_rilb01[87];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbaff01 = beliefPenaltyMPC_dl_aff + 120;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbaff01 = beliefPenaltyMPC_ds_aff + 120;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbcc01 = beliefPenaltyMPC_dl_cc + 120;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbcc01 = beliefPenaltyMPC_ds_cc + 120;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsl01 = beliefPenaltyMPC_ccrhs + 120;
int beliefPenaltyMPC_ubIdx01[33] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lub01 = beliefPenaltyMPC_l + 207;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_sub01 = beliefPenaltyMPC_s + 207;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lubbysub01 = beliefPenaltyMPC_lbys + 207;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_riub01[33];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubaff01 = beliefPenaltyMPC_dl_aff + 207;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubaff01 = beliefPenaltyMPC_ds_aff + 207;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubcc01 = beliefPenaltyMPC_dl_cc + 207;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubcc01 = beliefPenaltyMPC_ds_cc + 207;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsub01 = beliefPenaltyMPC_ccrhs + 207;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Phi01[3828];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_D01[87] = {-1.0000000000000000E+000, 
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
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_W01[2349];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ysd01[729];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lsd01[729];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_z02 = beliefPenaltyMPC_z + 174;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzaff02 = beliefPenaltyMPC_dz_aff + 174;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzcc02 = beliefPenaltyMPC_dz_cc + 174;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_rd02 = beliefPenaltyMPC_rd + 174;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lbyrd02[87];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_cost02 = beliefPenaltyMPC_grad_cost + 174;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_eq02 = beliefPenaltyMPC_grad_eq + 174;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_ineq02 = beliefPenaltyMPC_grad_ineq + 174;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ctv02[87];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_v02 = beliefPenaltyMPC_v + 54;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_re02[27];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_beta02[27];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_betacc02[27];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvaff02 = beliefPenaltyMPC_dv_aff + 54;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvcc02 = beliefPenaltyMPC_dv_cc + 54;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_V02[2349];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Yd02[378];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ld02[378];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_yy02[27];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_bmy02[27];
int beliefPenaltyMPC_lbIdx02[87] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llb02 = beliefPenaltyMPC_l + 240;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_slb02 = beliefPenaltyMPC_s + 240;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llbbyslb02 = beliefPenaltyMPC_lbys + 240;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_rilb02[87];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbaff02 = beliefPenaltyMPC_dl_aff + 240;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbaff02 = beliefPenaltyMPC_ds_aff + 240;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbcc02 = beliefPenaltyMPC_dl_cc + 240;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbcc02 = beliefPenaltyMPC_ds_cc + 240;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsl02 = beliefPenaltyMPC_ccrhs + 240;
int beliefPenaltyMPC_ubIdx02[33] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lub02 = beliefPenaltyMPC_l + 327;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_sub02 = beliefPenaltyMPC_s + 327;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lubbysub02 = beliefPenaltyMPC_lbys + 327;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_riub02[33];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubaff02 = beliefPenaltyMPC_dl_aff + 327;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubaff02 = beliefPenaltyMPC_ds_aff + 327;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubcc02 = beliefPenaltyMPC_dl_cc + 327;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubcc02 = beliefPenaltyMPC_ds_cc + 327;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsub02 = beliefPenaltyMPC_ccrhs + 327;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Phi02[3828];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_W02[2349];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ysd02[729];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lsd02[729];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_z03 = beliefPenaltyMPC_z + 261;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzaff03 = beliefPenaltyMPC_dz_aff + 261;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzcc03 = beliefPenaltyMPC_dz_cc + 261;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_rd03 = beliefPenaltyMPC_rd + 261;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lbyrd03[87];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_cost03 = beliefPenaltyMPC_grad_cost + 261;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_eq03 = beliefPenaltyMPC_grad_eq + 261;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_ineq03 = beliefPenaltyMPC_grad_ineq + 261;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ctv03[87];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_v03 = beliefPenaltyMPC_v + 81;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_re03[27];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_beta03[27];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_betacc03[27];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvaff03 = beliefPenaltyMPC_dv_aff + 81;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvcc03 = beliefPenaltyMPC_dv_cc + 81;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_V03[2349];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Yd03[378];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ld03[378];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_yy03[27];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_bmy03[27];
int beliefPenaltyMPC_lbIdx03[87] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llb03 = beliefPenaltyMPC_l + 360;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_slb03 = beliefPenaltyMPC_s + 360;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llbbyslb03 = beliefPenaltyMPC_lbys + 360;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_rilb03[87];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbaff03 = beliefPenaltyMPC_dl_aff + 360;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbaff03 = beliefPenaltyMPC_ds_aff + 360;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbcc03 = beliefPenaltyMPC_dl_cc + 360;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbcc03 = beliefPenaltyMPC_ds_cc + 360;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsl03 = beliefPenaltyMPC_ccrhs + 360;
int beliefPenaltyMPC_ubIdx03[33] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lub03 = beliefPenaltyMPC_l + 447;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_sub03 = beliefPenaltyMPC_s + 447;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lubbysub03 = beliefPenaltyMPC_lbys + 447;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_riub03[33];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubaff03 = beliefPenaltyMPC_dl_aff + 447;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubaff03 = beliefPenaltyMPC_ds_aff + 447;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubcc03 = beliefPenaltyMPC_dl_cc + 447;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubcc03 = beliefPenaltyMPC_ds_cc + 447;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsub03 = beliefPenaltyMPC_ccrhs + 447;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Phi03[3828];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_W03[2349];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ysd03[729];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lsd03[729];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_z04 = beliefPenaltyMPC_z + 348;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzaff04 = beliefPenaltyMPC_dz_aff + 348;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzcc04 = beliefPenaltyMPC_dz_cc + 348;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_rd04 = beliefPenaltyMPC_rd + 348;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lbyrd04[87];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_cost04 = beliefPenaltyMPC_grad_cost + 348;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_eq04 = beliefPenaltyMPC_grad_eq + 348;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_ineq04 = beliefPenaltyMPC_grad_ineq + 348;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ctv04[87];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_v04 = beliefPenaltyMPC_v + 108;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_re04[27];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_beta04[27];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_betacc04[27];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvaff04 = beliefPenaltyMPC_dv_aff + 108;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvcc04 = beliefPenaltyMPC_dv_cc + 108;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_V04[2349];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Yd04[378];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ld04[378];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_yy04[27];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_bmy04[27];
int beliefPenaltyMPC_lbIdx04[87] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llb04 = beliefPenaltyMPC_l + 480;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_slb04 = beliefPenaltyMPC_s + 480;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llbbyslb04 = beliefPenaltyMPC_lbys + 480;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_rilb04[87];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbaff04 = beliefPenaltyMPC_dl_aff + 480;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbaff04 = beliefPenaltyMPC_ds_aff + 480;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbcc04 = beliefPenaltyMPC_dl_cc + 480;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbcc04 = beliefPenaltyMPC_ds_cc + 480;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsl04 = beliefPenaltyMPC_ccrhs + 480;
int beliefPenaltyMPC_ubIdx04[33] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lub04 = beliefPenaltyMPC_l + 567;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_sub04 = beliefPenaltyMPC_s + 567;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lubbysub04 = beliefPenaltyMPC_lbys + 567;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_riub04[33];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubaff04 = beliefPenaltyMPC_dl_aff + 567;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubaff04 = beliefPenaltyMPC_ds_aff + 567;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubcc04 = beliefPenaltyMPC_dl_cc + 567;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubcc04 = beliefPenaltyMPC_ds_cc + 567;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsub04 = beliefPenaltyMPC_ccrhs + 567;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Phi04[3828];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_W04[2349];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ysd04[729];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lsd04[729];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_z05 = beliefPenaltyMPC_z + 435;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzaff05 = beliefPenaltyMPC_dz_aff + 435;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzcc05 = beliefPenaltyMPC_dz_cc + 435;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_rd05 = beliefPenaltyMPC_rd + 435;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lbyrd05[87];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_cost05 = beliefPenaltyMPC_grad_cost + 435;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_eq05 = beliefPenaltyMPC_grad_eq + 435;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_ineq05 = beliefPenaltyMPC_grad_ineq + 435;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ctv05[87];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_v05 = beliefPenaltyMPC_v + 135;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_re05[27];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_beta05[27];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_betacc05[27];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvaff05 = beliefPenaltyMPC_dv_aff + 135;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvcc05 = beliefPenaltyMPC_dv_cc + 135;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_V05[2349];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Yd05[378];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ld05[378];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_yy05[27];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_bmy05[27];
int beliefPenaltyMPC_lbIdx05[87] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llb05 = beliefPenaltyMPC_l + 600;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_slb05 = beliefPenaltyMPC_s + 600;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llbbyslb05 = beliefPenaltyMPC_lbys + 600;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_rilb05[87];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbaff05 = beliefPenaltyMPC_dl_aff + 600;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbaff05 = beliefPenaltyMPC_ds_aff + 600;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbcc05 = beliefPenaltyMPC_dl_cc + 600;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbcc05 = beliefPenaltyMPC_ds_cc + 600;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsl05 = beliefPenaltyMPC_ccrhs + 600;
int beliefPenaltyMPC_ubIdx05[33] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lub05 = beliefPenaltyMPC_l + 687;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_sub05 = beliefPenaltyMPC_s + 687;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lubbysub05 = beliefPenaltyMPC_lbys + 687;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_riub05[33];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubaff05 = beliefPenaltyMPC_dl_aff + 687;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubaff05 = beliefPenaltyMPC_ds_aff + 687;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubcc05 = beliefPenaltyMPC_dl_cc + 687;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubcc05 = beliefPenaltyMPC_ds_cc + 687;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsub05 = beliefPenaltyMPC_ccrhs + 687;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Phi05[3828];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_W05[2349];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ysd05[729];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lsd05[729];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_z06 = beliefPenaltyMPC_z + 522;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzaff06 = beliefPenaltyMPC_dz_aff + 522;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzcc06 = beliefPenaltyMPC_dz_cc + 522;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_rd06 = beliefPenaltyMPC_rd + 522;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lbyrd06[87];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_cost06 = beliefPenaltyMPC_grad_cost + 522;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_eq06 = beliefPenaltyMPC_grad_eq + 522;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_ineq06 = beliefPenaltyMPC_grad_ineq + 522;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ctv06[87];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_v06 = beliefPenaltyMPC_v + 162;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_re06[27];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_beta06[27];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_betacc06[27];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvaff06 = beliefPenaltyMPC_dv_aff + 162;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvcc06 = beliefPenaltyMPC_dv_cc + 162;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_V06[2349];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Yd06[378];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ld06[378];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_yy06[27];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_bmy06[27];
int beliefPenaltyMPC_lbIdx06[87] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llb06 = beliefPenaltyMPC_l + 720;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_slb06 = beliefPenaltyMPC_s + 720;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llbbyslb06 = beliefPenaltyMPC_lbys + 720;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_rilb06[87];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbaff06 = beliefPenaltyMPC_dl_aff + 720;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbaff06 = beliefPenaltyMPC_ds_aff + 720;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbcc06 = beliefPenaltyMPC_dl_cc + 720;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbcc06 = beliefPenaltyMPC_ds_cc + 720;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsl06 = beliefPenaltyMPC_ccrhs + 720;
int beliefPenaltyMPC_ubIdx06[33] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lub06 = beliefPenaltyMPC_l + 807;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_sub06 = beliefPenaltyMPC_s + 807;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lubbysub06 = beliefPenaltyMPC_lbys + 807;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_riub06[33];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubaff06 = beliefPenaltyMPC_dl_aff + 807;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubaff06 = beliefPenaltyMPC_ds_aff + 807;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubcc06 = beliefPenaltyMPC_dl_cc + 807;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubcc06 = beliefPenaltyMPC_ds_cc + 807;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsub06 = beliefPenaltyMPC_ccrhs + 807;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Phi06[3828];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_W06[2349];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ysd06[729];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lsd06[729];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_z07 = beliefPenaltyMPC_z + 609;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzaff07 = beliefPenaltyMPC_dz_aff + 609;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzcc07 = beliefPenaltyMPC_dz_cc + 609;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_rd07 = beliefPenaltyMPC_rd + 609;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lbyrd07[87];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_cost07 = beliefPenaltyMPC_grad_cost + 609;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_eq07 = beliefPenaltyMPC_grad_eq + 609;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_ineq07 = beliefPenaltyMPC_grad_ineq + 609;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ctv07[87];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_v07 = beliefPenaltyMPC_v + 189;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_re07[27];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_beta07[27];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_betacc07[27];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvaff07 = beliefPenaltyMPC_dv_aff + 189;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvcc07 = beliefPenaltyMPC_dv_cc + 189;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_V07[2349];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Yd07[378];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ld07[378];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_yy07[27];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_bmy07[27];
int beliefPenaltyMPC_lbIdx07[87] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llb07 = beliefPenaltyMPC_l + 840;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_slb07 = beliefPenaltyMPC_s + 840;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llbbyslb07 = beliefPenaltyMPC_lbys + 840;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_rilb07[87];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbaff07 = beliefPenaltyMPC_dl_aff + 840;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbaff07 = beliefPenaltyMPC_ds_aff + 840;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbcc07 = beliefPenaltyMPC_dl_cc + 840;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbcc07 = beliefPenaltyMPC_ds_cc + 840;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsl07 = beliefPenaltyMPC_ccrhs + 840;
int beliefPenaltyMPC_ubIdx07[33] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lub07 = beliefPenaltyMPC_l + 927;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_sub07 = beliefPenaltyMPC_s + 927;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lubbysub07 = beliefPenaltyMPC_lbys + 927;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_riub07[33];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubaff07 = beliefPenaltyMPC_dl_aff + 927;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubaff07 = beliefPenaltyMPC_ds_aff + 927;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubcc07 = beliefPenaltyMPC_dl_cc + 927;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubcc07 = beliefPenaltyMPC_ds_cc + 927;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsub07 = beliefPenaltyMPC_ccrhs + 927;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Phi07[3828];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_W07[2349];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ysd07[729];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lsd07[729];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_z08 = beliefPenaltyMPC_z + 696;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzaff08 = beliefPenaltyMPC_dz_aff + 696;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzcc08 = beliefPenaltyMPC_dz_cc + 696;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_rd08 = beliefPenaltyMPC_rd + 696;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lbyrd08[87];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_cost08 = beliefPenaltyMPC_grad_cost + 696;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_eq08 = beliefPenaltyMPC_grad_eq + 696;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_ineq08 = beliefPenaltyMPC_grad_ineq + 696;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ctv08[87];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_v08 = beliefPenaltyMPC_v + 216;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_re08[27];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_beta08[27];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_betacc08[27];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvaff08 = beliefPenaltyMPC_dv_aff + 216;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvcc08 = beliefPenaltyMPC_dv_cc + 216;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_V08[2349];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Yd08[378];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ld08[378];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_yy08[27];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_bmy08[27];
int beliefPenaltyMPC_lbIdx08[87] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llb08 = beliefPenaltyMPC_l + 960;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_slb08 = beliefPenaltyMPC_s + 960;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llbbyslb08 = beliefPenaltyMPC_lbys + 960;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_rilb08[87];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbaff08 = beliefPenaltyMPC_dl_aff + 960;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbaff08 = beliefPenaltyMPC_ds_aff + 960;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbcc08 = beliefPenaltyMPC_dl_cc + 960;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbcc08 = beliefPenaltyMPC_ds_cc + 960;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsl08 = beliefPenaltyMPC_ccrhs + 960;
int beliefPenaltyMPC_ubIdx08[33] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lub08 = beliefPenaltyMPC_l + 1047;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_sub08 = beliefPenaltyMPC_s + 1047;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lubbysub08 = beliefPenaltyMPC_lbys + 1047;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_riub08[33];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubaff08 = beliefPenaltyMPC_dl_aff + 1047;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubaff08 = beliefPenaltyMPC_ds_aff + 1047;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubcc08 = beliefPenaltyMPC_dl_cc + 1047;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubcc08 = beliefPenaltyMPC_ds_cc + 1047;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsub08 = beliefPenaltyMPC_ccrhs + 1047;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Phi08[3828];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_W08[2349];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ysd08[729];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lsd08[729];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_z09 = beliefPenaltyMPC_z + 783;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzaff09 = beliefPenaltyMPC_dz_aff + 783;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzcc09 = beliefPenaltyMPC_dz_cc + 783;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_rd09 = beliefPenaltyMPC_rd + 783;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lbyrd09[87];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_cost09 = beliefPenaltyMPC_grad_cost + 783;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_eq09 = beliefPenaltyMPC_grad_eq + 783;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_ineq09 = beliefPenaltyMPC_grad_ineq + 783;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ctv09[87];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_v09 = beliefPenaltyMPC_v + 243;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_re09[27];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_beta09[27];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_betacc09[27];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvaff09 = beliefPenaltyMPC_dv_aff + 243;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvcc09 = beliefPenaltyMPC_dv_cc + 243;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_V09[2349];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Yd09[378];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ld09[378];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_yy09[27];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_bmy09[27];
int beliefPenaltyMPC_lbIdx09[87] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llb09 = beliefPenaltyMPC_l + 1080;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_slb09 = beliefPenaltyMPC_s + 1080;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llbbyslb09 = beliefPenaltyMPC_lbys + 1080;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_rilb09[87];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbaff09 = beliefPenaltyMPC_dl_aff + 1080;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbaff09 = beliefPenaltyMPC_ds_aff + 1080;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbcc09 = beliefPenaltyMPC_dl_cc + 1080;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbcc09 = beliefPenaltyMPC_ds_cc + 1080;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsl09 = beliefPenaltyMPC_ccrhs + 1080;
int beliefPenaltyMPC_ubIdx09[33] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lub09 = beliefPenaltyMPC_l + 1167;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_sub09 = beliefPenaltyMPC_s + 1167;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lubbysub09 = beliefPenaltyMPC_lbys + 1167;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_riub09[33];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubaff09 = beliefPenaltyMPC_dl_aff + 1167;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubaff09 = beliefPenaltyMPC_ds_aff + 1167;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubcc09 = beliefPenaltyMPC_dl_cc + 1167;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubcc09 = beliefPenaltyMPC_ds_cc + 1167;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsub09 = beliefPenaltyMPC_ccrhs + 1167;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Phi09[3828];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_W09[2349];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ysd09[729];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lsd09[729];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_z10 = beliefPenaltyMPC_z + 870;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzaff10 = beliefPenaltyMPC_dz_aff + 870;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzcc10 = beliefPenaltyMPC_dz_cc + 870;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_rd10 = beliefPenaltyMPC_rd + 870;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lbyrd10[87];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_cost10 = beliefPenaltyMPC_grad_cost + 870;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_eq10 = beliefPenaltyMPC_grad_eq + 870;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_ineq10 = beliefPenaltyMPC_grad_ineq + 870;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ctv10[87];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_v10 = beliefPenaltyMPC_v + 270;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_re10[27];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_beta10[27];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_betacc10[27];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvaff10 = beliefPenaltyMPC_dv_aff + 270;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvcc10 = beliefPenaltyMPC_dv_cc + 270;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_V10[2349];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Yd10[378];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ld10[378];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_yy10[27];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_bmy10[27];
int beliefPenaltyMPC_lbIdx10[87] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llb10 = beliefPenaltyMPC_l + 1200;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_slb10 = beliefPenaltyMPC_s + 1200;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llbbyslb10 = beliefPenaltyMPC_lbys + 1200;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_rilb10[87];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbaff10 = beliefPenaltyMPC_dl_aff + 1200;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbaff10 = beliefPenaltyMPC_ds_aff + 1200;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbcc10 = beliefPenaltyMPC_dl_cc + 1200;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbcc10 = beliefPenaltyMPC_ds_cc + 1200;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsl10 = beliefPenaltyMPC_ccrhs + 1200;
int beliefPenaltyMPC_ubIdx10[33] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lub10 = beliefPenaltyMPC_l + 1287;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_sub10 = beliefPenaltyMPC_s + 1287;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lubbysub10 = beliefPenaltyMPC_lbys + 1287;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_riub10[33];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubaff10 = beliefPenaltyMPC_dl_aff + 1287;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubaff10 = beliefPenaltyMPC_ds_aff + 1287;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubcc10 = beliefPenaltyMPC_dl_cc + 1287;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubcc10 = beliefPenaltyMPC_ds_cc + 1287;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsub10 = beliefPenaltyMPC_ccrhs + 1287;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Phi10[3828];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_W10[2349];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ysd10[729];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lsd10[729];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_z11 = beliefPenaltyMPC_z + 957;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzaff11 = beliefPenaltyMPC_dz_aff + 957;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzcc11 = beliefPenaltyMPC_dz_cc + 957;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_rd11 = beliefPenaltyMPC_rd + 957;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lbyrd11[87];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_cost11 = beliefPenaltyMPC_grad_cost + 957;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_eq11 = beliefPenaltyMPC_grad_eq + 957;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_ineq11 = beliefPenaltyMPC_grad_ineq + 957;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ctv11[87];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_v11 = beliefPenaltyMPC_v + 297;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_re11[27];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_beta11[27];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_betacc11[27];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvaff11 = beliefPenaltyMPC_dv_aff + 297;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvcc11 = beliefPenaltyMPC_dv_cc + 297;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_V11[2349];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Yd11[378];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ld11[378];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_yy11[27];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_bmy11[27];
int beliefPenaltyMPC_lbIdx11[87] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llb11 = beliefPenaltyMPC_l + 1320;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_slb11 = beliefPenaltyMPC_s + 1320;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llbbyslb11 = beliefPenaltyMPC_lbys + 1320;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_rilb11[87];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbaff11 = beliefPenaltyMPC_dl_aff + 1320;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbaff11 = beliefPenaltyMPC_ds_aff + 1320;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbcc11 = beliefPenaltyMPC_dl_cc + 1320;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbcc11 = beliefPenaltyMPC_ds_cc + 1320;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsl11 = beliefPenaltyMPC_ccrhs + 1320;
int beliefPenaltyMPC_ubIdx11[33] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lub11 = beliefPenaltyMPC_l + 1407;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_sub11 = beliefPenaltyMPC_s + 1407;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lubbysub11 = beliefPenaltyMPC_lbys + 1407;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_riub11[33];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubaff11 = beliefPenaltyMPC_dl_aff + 1407;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubaff11 = beliefPenaltyMPC_ds_aff + 1407;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubcc11 = beliefPenaltyMPC_dl_cc + 1407;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubcc11 = beliefPenaltyMPC_ds_cc + 1407;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsub11 = beliefPenaltyMPC_ccrhs + 1407;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Phi11[3828];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_W11[2349];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ysd11[729];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lsd11[729];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_z12 = beliefPenaltyMPC_z + 1044;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzaff12 = beliefPenaltyMPC_dz_aff + 1044;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzcc12 = beliefPenaltyMPC_dz_cc + 1044;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_rd12 = beliefPenaltyMPC_rd + 1044;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lbyrd12[87];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_cost12 = beliefPenaltyMPC_grad_cost + 1044;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_eq12 = beliefPenaltyMPC_grad_eq + 1044;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_ineq12 = beliefPenaltyMPC_grad_ineq + 1044;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ctv12[87];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_v12 = beliefPenaltyMPC_v + 324;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_re12[27];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_beta12[27];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_betacc12[27];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvaff12 = beliefPenaltyMPC_dv_aff + 324;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvcc12 = beliefPenaltyMPC_dv_cc + 324;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_V12[2349];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Yd12[378];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ld12[378];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_yy12[27];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_bmy12[27];
int beliefPenaltyMPC_lbIdx12[87] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llb12 = beliefPenaltyMPC_l + 1440;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_slb12 = beliefPenaltyMPC_s + 1440;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llbbyslb12 = beliefPenaltyMPC_lbys + 1440;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_rilb12[87];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbaff12 = beliefPenaltyMPC_dl_aff + 1440;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbaff12 = beliefPenaltyMPC_ds_aff + 1440;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbcc12 = beliefPenaltyMPC_dl_cc + 1440;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbcc12 = beliefPenaltyMPC_ds_cc + 1440;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsl12 = beliefPenaltyMPC_ccrhs + 1440;
int beliefPenaltyMPC_ubIdx12[33] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lub12 = beliefPenaltyMPC_l + 1527;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_sub12 = beliefPenaltyMPC_s + 1527;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lubbysub12 = beliefPenaltyMPC_lbys + 1527;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_riub12[33];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubaff12 = beliefPenaltyMPC_dl_aff + 1527;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubaff12 = beliefPenaltyMPC_ds_aff + 1527;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubcc12 = beliefPenaltyMPC_dl_cc + 1527;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubcc12 = beliefPenaltyMPC_ds_cc + 1527;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsub12 = beliefPenaltyMPC_ccrhs + 1527;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Phi12[3828];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_W12[2349];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ysd12[729];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lsd12[729];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_z13 = beliefPenaltyMPC_z + 1131;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzaff13 = beliefPenaltyMPC_dz_aff + 1131;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzcc13 = beliefPenaltyMPC_dz_cc + 1131;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_rd13 = beliefPenaltyMPC_rd + 1131;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lbyrd13[87];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_cost13 = beliefPenaltyMPC_grad_cost + 1131;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_eq13 = beliefPenaltyMPC_grad_eq + 1131;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_ineq13 = beliefPenaltyMPC_grad_ineq + 1131;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ctv13[87];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_v13 = beliefPenaltyMPC_v + 351;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_re13[27];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_beta13[27];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_betacc13[27];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvaff13 = beliefPenaltyMPC_dv_aff + 351;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvcc13 = beliefPenaltyMPC_dv_cc + 351;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_V13[2349];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Yd13[378];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ld13[378];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_yy13[27];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_bmy13[27];
int beliefPenaltyMPC_lbIdx13[87] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llb13 = beliefPenaltyMPC_l + 1560;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_slb13 = beliefPenaltyMPC_s + 1560;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llbbyslb13 = beliefPenaltyMPC_lbys + 1560;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_rilb13[87];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbaff13 = beliefPenaltyMPC_dl_aff + 1560;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbaff13 = beliefPenaltyMPC_ds_aff + 1560;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbcc13 = beliefPenaltyMPC_dl_cc + 1560;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbcc13 = beliefPenaltyMPC_ds_cc + 1560;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsl13 = beliefPenaltyMPC_ccrhs + 1560;
int beliefPenaltyMPC_ubIdx13[33] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lub13 = beliefPenaltyMPC_l + 1647;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_sub13 = beliefPenaltyMPC_s + 1647;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lubbysub13 = beliefPenaltyMPC_lbys + 1647;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_riub13[33];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubaff13 = beliefPenaltyMPC_dl_aff + 1647;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubaff13 = beliefPenaltyMPC_ds_aff + 1647;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubcc13 = beliefPenaltyMPC_dl_cc + 1647;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubcc13 = beliefPenaltyMPC_ds_cc + 1647;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsub13 = beliefPenaltyMPC_ccrhs + 1647;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Phi13[3828];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_W13[2349];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ysd13[729];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lsd13[729];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_z14 = beliefPenaltyMPC_z + 1218;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzaff14 = beliefPenaltyMPC_dz_aff + 1218;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzcc14 = beliefPenaltyMPC_dz_cc + 1218;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_rd14 = beliefPenaltyMPC_rd + 1218;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lbyrd14[33];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_cost14 = beliefPenaltyMPC_grad_cost + 1218;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_eq14 = beliefPenaltyMPC_grad_eq + 1218;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_ineq14 = beliefPenaltyMPC_grad_ineq + 1218;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ctv14[33];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_v14 = beliefPenaltyMPC_v + 378;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_re14[27];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_beta14[27];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_betacc14[27];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvaff14 = beliefPenaltyMPC_dv_aff + 378;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvcc14 = beliefPenaltyMPC_dv_cc + 378;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_V14[891];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Yd14[378];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ld14[378];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_yy14[27];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_bmy14[27];
int beliefPenaltyMPC_lbIdx14[33] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llb14 = beliefPenaltyMPC_l + 1680;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_slb14 = beliefPenaltyMPC_s + 1680;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llbbyslb14 = beliefPenaltyMPC_lbys + 1680;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_rilb14[33];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbaff14 = beliefPenaltyMPC_dl_aff + 1680;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbaff14 = beliefPenaltyMPC_ds_aff + 1680;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbcc14 = beliefPenaltyMPC_dl_cc + 1680;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbcc14 = beliefPenaltyMPC_ds_cc + 1680;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsl14 = beliefPenaltyMPC_ccrhs + 1680;
int beliefPenaltyMPC_ubIdx14[27] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lub14 = beliefPenaltyMPC_l + 1713;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_sub14 = beliefPenaltyMPC_s + 1713;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lubbysub14 = beliefPenaltyMPC_lbys + 1713;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_riub14[27];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubaff14 = beliefPenaltyMPC_dl_aff + 1713;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubaff14 = beliefPenaltyMPC_ds_aff + 1713;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubcc14 = beliefPenaltyMPC_dl_cc + 1713;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubcc14 = beliefPenaltyMPC_ds_cc + 1713;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsub14 = beliefPenaltyMPC_ccrhs + 1713;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_sp14 = beliefPenaltyMPC_s + 1740;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lp14 = beliefPenaltyMPC_l + 1740;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lpbysp14 = beliefPenaltyMPC_lbys + 1740;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlp_aff14 = beliefPenaltyMPC_dl_aff + 1740;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsp_aff14 = beliefPenaltyMPC_ds_aff + 1740;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlp_cc14 = beliefPenaltyMPC_dl_cc + 1740;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsp_cc14 = beliefPenaltyMPC_ds_cc + 1740;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsp14 = beliefPenaltyMPC_ccrhs + 1740;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_rip14[6];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Phi14[561];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_D14[33] = {-1.0000000000000000E+000, 
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
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_W14[891];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ysd14[729];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lsd14[729];
beliefPenaltyMPC_FLOAT musigma;
beliefPenaltyMPC_FLOAT sigma_3rdroot;




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
beliefPenaltyMPC_LA_INITIALIZEVECTOR_1251(beliefPenaltyMPC_z, 0);
beliefPenaltyMPC_LA_INITIALIZEVECTOR_405(beliefPenaltyMPC_v, 1);
beliefPenaltyMPC_LA_INITIALIZEVECTOR_1746(beliefPenaltyMPC_l, 1);
beliefPenaltyMPC_LA_INITIALIZEVECTOR_1746(beliefPenaltyMPC_s, 1);
info->mu = 0;
beliefPenaltyMPC_LA_DOTACC_1746(beliefPenaltyMPC_l, beliefPenaltyMPC_s, &info->mu);
info->mu /= 1746;
while( 1 ){
info->pobj = 0;
beliefPenaltyMPC_LA_DENSE_QUADFCN_87(params->Q1, params->f1, beliefPenaltyMPC_z00, beliefPenaltyMPC_grad_cost00, &info->pobj);
beliefPenaltyMPC_LA_DENSE_QUADFCN_87(params->Q2, params->f2, beliefPenaltyMPC_z01, beliefPenaltyMPC_grad_cost01, &info->pobj);
beliefPenaltyMPC_LA_DENSE_QUADFCN_87(params->Q3, params->f3, beliefPenaltyMPC_z02, beliefPenaltyMPC_grad_cost02, &info->pobj);
beliefPenaltyMPC_LA_DENSE_QUADFCN_87(params->Q4, params->f4, beliefPenaltyMPC_z03, beliefPenaltyMPC_grad_cost03, &info->pobj);
beliefPenaltyMPC_LA_DENSE_QUADFCN_87(params->Q5, params->f5, beliefPenaltyMPC_z04, beliefPenaltyMPC_grad_cost04, &info->pobj);
beliefPenaltyMPC_LA_DENSE_QUADFCN_87(params->Q6, params->f6, beliefPenaltyMPC_z05, beliefPenaltyMPC_grad_cost05, &info->pobj);
beliefPenaltyMPC_LA_DENSE_QUADFCN_87(params->Q7, params->f7, beliefPenaltyMPC_z06, beliefPenaltyMPC_grad_cost06, &info->pobj);
beliefPenaltyMPC_LA_DENSE_QUADFCN_87(params->Q8, params->f8, beliefPenaltyMPC_z07, beliefPenaltyMPC_grad_cost07, &info->pobj);
beliefPenaltyMPC_LA_DENSE_QUADFCN_87(params->Q9, params->f9, beliefPenaltyMPC_z08, beliefPenaltyMPC_grad_cost08, &info->pobj);
beliefPenaltyMPC_LA_DENSE_QUADFCN_87(params->Q10, params->f10, beliefPenaltyMPC_z09, beliefPenaltyMPC_grad_cost09, &info->pobj);
beliefPenaltyMPC_LA_DENSE_QUADFCN_87(params->Q11, params->f11, beliefPenaltyMPC_z10, beliefPenaltyMPC_grad_cost10, &info->pobj);
beliefPenaltyMPC_LA_DENSE_QUADFCN_87(params->Q12, params->f12, beliefPenaltyMPC_z11, beliefPenaltyMPC_grad_cost11, &info->pobj);
beliefPenaltyMPC_LA_DENSE_QUADFCN_87(params->Q13, params->f13, beliefPenaltyMPC_z12, beliefPenaltyMPC_grad_cost12, &info->pobj);
beliefPenaltyMPC_LA_DENSE_QUADFCN_87(params->Q14, params->f14, beliefPenaltyMPC_z13, beliefPenaltyMPC_grad_cost13, &info->pobj);
beliefPenaltyMPC_LA_DENSE_QUADFCN_33(params->Q15, params->f15, beliefPenaltyMPC_z14, beliefPenaltyMPC_grad_cost14, &info->pobj);
info->res_eq = 0;
info->dgap = 0;
beliefPenaltyMPC_LA_DIAGZERO_MVMSUB6_27(beliefPenaltyMPC_D00, beliefPenaltyMPC_z00, params->e1, beliefPenaltyMPC_v00, beliefPenaltyMPC_re00, &info->dgap, &info->res_eq);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_27_87_87(params->C1, beliefPenaltyMPC_z00, beliefPenaltyMPC_D01, beliefPenaltyMPC_z01, params->e2, beliefPenaltyMPC_v01, beliefPenaltyMPC_re01, &info->dgap, &info->res_eq);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_27_87_87(params->C2, beliefPenaltyMPC_z01, beliefPenaltyMPC_D01, beliefPenaltyMPC_z02, params->e3, beliefPenaltyMPC_v02, beliefPenaltyMPC_re02, &info->dgap, &info->res_eq);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_27_87_87(params->C3, beliefPenaltyMPC_z02, beliefPenaltyMPC_D01, beliefPenaltyMPC_z03, params->e4, beliefPenaltyMPC_v03, beliefPenaltyMPC_re03, &info->dgap, &info->res_eq);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_27_87_87(params->C4, beliefPenaltyMPC_z03, beliefPenaltyMPC_D01, beliefPenaltyMPC_z04, params->e5, beliefPenaltyMPC_v04, beliefPenaltyMPC_re04, &info->dgap, &info->res_eq);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_27_87_87(params->C5, beliefPenaltyMPC_z04, beliefPenaltyMPC_D01, beliefPenaltyMPC_z05, params->e6, beliefPenaltyMPC_v05, beliefPenaltyMPC_re05, &info->dgap, &info->res_eq);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_27_87_87(params->C6, beliefPenaltyMPC_z05, beliefPenaltyMPC_D01, beliefPenaltyMPC_z06, params->e7, beliefPenaltyMPC_v06, beliefPenaltyMPC_re06, &info->dgap, &info->res_eq);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_27_87_87(params->C7, beliefPenaltyMPC_z06, beliefPenaltyMPC_D01, beliefPenaltyMPC_z07, params->e8, beliefPenaltyMPC_v07, beliefPenaltyMPC_re07, &info->dgap, &info->res_eq);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_27_87_87(params->C8, beliefPenaltyMPC_z07, beliefPenaltyMPC_D01, beliefPenaltyMPC_z08, params->e9, beliefPenaltyMPC_v08, beliefPenaltyMPC_re08, &info->dgap, &info->res_eq);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_27_87_87(params->C9, beliefPenaltyMPC_z08, beliefPenaltyMPC_D01, beliefPenaltyMPC_z09, params->e10, beliefPenaltyMPC_v09, beliefPenaltyMPC_re09, &info->dgap, &info->res_eq);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_27_87_87(params->C10, beliefPenaltyMPC_z09, beliefPenaltyMPC_D01, beliefPenaltyMPC_z10, params->e11, beliefPenaltyMPC_v10, beliefPenaltyMPC_re10, &info->dgap, &info->res_eq);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_27_87_87(params->C11, beliefPenaltyMPC_z10, beliefPenaltyMPC_D01, beliefPenaltyMPC_z11, params->e12, beliefPenaltyMPC_v11, beliefPenaltyMPC_re11, &info->dgap, &info->res_eq);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_27_87_87(params->C12, beliefPenaltyMPC_z11, beliefPenaltyMPC_D01, beliefPenaltyMPC_z12, params->e13, beliefPenaltyMPC_v12, beliefPenaltyMPC_re12, &info->dgap, &info->res_eq);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_27_87_87(params->C13, beliefPenaltyMPC_z12, beliefPenaltyMPC_D01, beliefPenaltyMPC_z13, params->e14, beliefPenaltyMPC_v13, beliefPenaltyMPC_re13, &info->dgap, &info->res_eq);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MVMSUB3_27_87_33(params->C14, beliefPenaltyMPC_z13, beliefPenaltyMPC_D14, beliefPenaltyMPC_z14, params->e15, beliefPenaltyMPC_v14, beliefPenaltyMPC_re14, &info->dgap, &info->res_eq);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_27_87_27(params->C1, beliefPenaltyMPC_v01, beliefPenaltyMPC_D00, beliefPenaltyMPC_v00, beliefPenaltyMPC_grad_eq00);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_27_87_27(params->C2, beliefPenaltyMPC_v02, beliefPenaltyMPC_D01, beliefPenaltyMPC_v01, beliefPenaltyMPC_grad_eq01);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_27_87_27(params->C3, beliefPenaltyMPC_v03, beliefPenaltyMPC_D01, beliefPenaltyMPC_v02, beliefPenaltyMPC_grad_eq02);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_27_87_27(params->C4, beliefPenaltyMPC_v04, beliefPenaltyMPC_D01, beliefPenaltyMPC_v03, beliefPenaltyMPC_grad_eq03);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_27_87_27(params->C5, beliefPenaltyMPC_v05, beliefPenaltyMPC_D01, beliefPenaltyMPC_v04, beliefPenaltyMPC_grad_eq04);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_27_87_27(params->C6, beliefPenaltyMPC_v06, beliefPenaltyMPC_D01, beliefPenaltyMPC_v05, beliefPenaltyMPC_grad_eq05);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_27_87_27(params->C7, beliefPenaltyMPC_v07, beliefPenaltyMPC_D01, beliefPenaltyMPC_v06, beliefPenaltyMPC_grad_eq06);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_27_87_27(params->C8, beliefPenaltyMPC_v08, beliefPenaltyMPC_D01, beliefPenaltyMPC_v07, beliefPenaltyMPC_grad_eq07);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_27_87_27(params->C9, beliefPenaltyMPC_v09, beliefPenaltyMPC_D01, beliefPenaltyMPC_v08, beliefPenaltyMPC_grad_eq08);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_27_87_27(params->C10, beliefPenaltyMPC_v10, beliefPenaltyMPC_D01, beliefPenaltyMPC_v09, beliefPenaltyMPC_grad_eq09);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_27_87_27(params->C11, beliefPenaltyMPC_v11, beliefPenaltyMPC_D01, beliefPenaltyMPC_v10, beliefPenaltyMPC_grad_eq10);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_27_87_27(params->C12, beliefPenaltyMPC_v12, beliefPenaltyMPC_D01, beliefPenaltyMPC_v11, beliefPenaltyMPC_grad_eq11);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_27_87_27(params->C13, beliefPenaltyMPC_v13, beliefPenaltyMPC_D01, beliefPenaltyMPC_v12, beliefPenaltyMPC_grad_eq12);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_27_87_27(params->C14, beliefPenaltyMPC_v14, beliefPenaltyMPC_D01, beliefPenaltyMPC_v13, beliefPenaltyMPC_grad_eq13);
beliefPenaltyMPC_LA_DIAGZERO_MTVM_27_33(beliefPenaltyMPC_D14, beliefPenaltyMPC_v14, beliefPenaltyMPC_grad_eq14);
info->res_ineq = 0;
beliefPenaltyMPC_LA_VSUBADD3_87(params->lb1, beliefPenaltyMPC_z00, beliefPenaltyMPC_lbIdx00, beliefPenaltyMPC_llb00, beliefPenaltyMPC_slb00, beliefPenaltyMPC_rilb00, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD2_33(beliefPenaltyMPC_z00, beliefPenaltyMPC_ubIdx00, params->ub1, beliefPenaltyMPC_lub00, beliefPenaltyMPC_sub00, beliefPenaltyMPC_riub00, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD3_87(params->lb2, beliefPenaltyMPC_z01, beliefPenaltyMPC_lbIdx01, beliefPenaltyMPC_llb01, beliefPenaltyMPC_slb01, beliefPenaltyMPC_rilb01, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD2_33(beliefPenaltyMPC_z01, beliefPenaltyMPC_ubIdx01, params->ub2, beliefPenaltyMPC_lub01, beliefPenaltyMPC_sub01, beliefPenaltyMPC_riub01, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD3_87(params->lb3, beliefPenaltyMPC_z02, beliefPenaltyMPC_lbIdx02, beliefPenaltyMPC_llb02, beliefPenaltyMPC_slb02, beliefPenaltyMPC_rilb02, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD2_33(beliefPenaltyMPC_z02, beliefPenaltyMPC_ubIdx02, params->ub3, beliefPenaltyMPC_lub02, beliefPenaltyMPC_sub02, beliefPenaltyMPC_riub02, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD3_87(params->lb4, beliefPenaltyMPC_z03, beliefPenaltyMPC_lbIdx03, beliefPenaltyMPC_llb03, beliefPenaltyMPC_slb03, beliefPenaltyMPC_rilb03, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD2_33(beliefPenaltyMPC_z03, beliefPenaltyMPC_ubIdx03, params->ub4, beliefPenaltyMPC_lub03, beliefPenaltyMPC_sub03, beliefPenaltyMPC_riub03, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD3_87(params->lb5, beliefPenaltyMPC_z04, beliefPenaltyMPC_lbIdx04, beliefPenaltyMPC_llb04, beliefPenaltyMPC_slb04, beliefPenaltyMPC_rilb04, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD2_33(beliefPenaltyMPC_z04, beliefPenaltyMPC_ubIdx04, params->ub5, beliefPenaltyMPC_lub04, beliefPenaltyMPC_sub04, beliefPenaltyMPC_riub04, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD3_87(params->lb6, beliefPenaltyMPC_z05, beliefPenaltyMPC_lbIdx05, beliefPenaltyMPC_llb05, beliefPenaltyMPC_slb05, beliefPenaltyMPC_rilb05, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD2_33(beliefPenaltyMPC_z05, beliefPenaltyMPC_ubIdx05, params->ub6, beliefPenaltyMPC_lub05, beliefPenaltyMPC_sub05, beliefPenaltyMPC_riub05, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD3_87(params->lb7, beliefPenaltyMPC_z06, beliefPenaltyMPC_lbIdx06, beliefPenaltyMPC_llb06, beliefPenaltyMPC_slb06, beliefPenaltyMPC_rilb06, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD2_33(beliefPenaltyMPC_z06, beliefPenaltyMPC_ubIdx06, params->ub7, beliefPenaltyMPC_lub06, beliefPenaltyMPC_sub06, beliefPenaltyMPC_riub06, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD3_87(params->lb8, beliefPenaltyMPC_z07, beliefPenaltyMPC_lbIdx07, beliefPenaltyMPC_llb07, beliefPenaltyMPC_slb07, beliefPenaltyMPC_rilb07, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD2_33(beliefPenaltyMPC_z07, beliefPenaltyMPC_ubIdx07, params->ub8, beliefPenaltyMPC_lub07, beliefPenaltyMPC_sub07, beliefPenaltyMPC_riub07, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD3_87(params->lb9, beliefPenaltyMPC_z08, beliefPenaltyMPC_lbIdx08, beliefPenaltyMPC_llb08, beliefPenaltyMPC_slb08, beliefPenaltyMPC_rilb08, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD2_33(beliefPenaltyMPC_z08, beliefPenaltyMPC_ubIdx08, params->ub9, beliefPenaltyMPC_lub08, beliefPenaltyMPC_sub08, beliefPenaltyMPC_riub08, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD3_87(params->lb10, beliefPenaltyMPC_z09, beliefPenaltyMPC_lbIdx09, beliefPenaltyMPC_llb09, beliefPenaltyMPC_slb09, beliefPenaltyMPC_rilb09, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD2_33(beliefPenaltyMPC_z09, beliefPenaltyMPC_ubIdx09, params->ub10, beliefPenaltyMPC_lub09, beliefPenaltyMPC_sub09, beliefPenaltyMPC_riub09, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD3_87(params->lb11, beliefPenaltyMPC_z10, beliefPenaltyMPC_lbIdx10, beliefPenaltyMPC_llb10, beliefPenaltyMPC_slb10, beliefPenaltyMPC_rilb10, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD2_33(beliefPenaltyMPC_z10, beliefPenaltyMPC_ubIdx10, params->ub11, beliefPenaltyMPC_lub10, beliefPenaltyMPC_sub10, beliefPenaltyMPC_riub10, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD3_87(params->lb12, beliefPenaltyMPC_z11, beliefPenaltyMPC_lbIdx11, beliefPenaltyMPC_llb11, beliefPenaltyMPC_slb11, beliefPenaltyMPC_rilb11, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD2_33(beliefPenaltyMPC_z11, beliefPenaltyMPC_ubIdx11, params->ub12, beliefPenaltyMPC_lub11, beliefPenaltyMPC_sub11, beliefPenaltyMPC_riub11, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD3_87(params->lb13, beliefPenaltyMPC_z12, beliefPenaltyMPC_lbIdx12, beliefPenaltyMPC_llb12, beliefPenaltyMPC_slb12, beliefPenaltyMPC_rilb12, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD2_33(beliefPenaltyMPC_z12, beliefPenaltyMPC_ubIdx12, params->ub13, beliefPenaltyMPC_lub12, beliefPenaltyMPC_sub12, beliefPenaltyMPC_riub12, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD3_87(params->lb14, beliefPenaltyMPC_z13, beliefPenaltyMPC_lbIdx13, beliefPenaltyMPC_llb13, beliefPenaltyMPC_slb13, beliefPenaltyMPC_rilb13, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD2_33(beliefPenaltyMPC_z13, beliefPenaltyMPC_ubIdx13, params->ub14, beliefPenaltyMPC_lub13, beliefPenaltyMPC_sub13, beliefPenaltyMPC_riub13, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD3_33(params->lb15, beliefPenaltyMPC_z14, beliefPenaltyMPC_lbIdx14, beliefPenaltyMPC_llb14, beliefPenaltyMPC_slb14, beliefPenaltyMPC_rilb14, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD2_27(beliefPenaltyMPC_z14, beliefPenaltyMPC_ubIdx14, params->ub15, beliefPenaltyMPC_lub14, beliefPenaltyMPC_sub14, beliefPenaltyMPC_riub14, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_MVSUBADD_6_33(params->A15, beliefPenaltyMPC_z14, params->b15, beliefPenaltyMPC_sp14, beliefPenaltyMPC_lp14, beliefPenaltyMPC_rip14, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_INEQ_B_GRAD_87_87_33(beliefPenaltyMPC_lub00, beliefPenaltyMPC_sub00, beliefPenaltyMPC_riub00, beliefPenaltyMPC_llb00, beliefPenaltyMPC_slb00, beliefPenaltyMPC_rilb00, beliefPenaltyMPC_lbIdx00, beliefPenaltyMPC_ubIdx00, beliefPenaltyMPC_grad_ineq00, beliefPenaltyMPC_lubbysub00, beliefPenaltyMPC_llbbyslb00);
beliefPenaltyMPC_LA_INEQ_B_GRAD_87_87_33(beliefPenaltyMPC_lub01, beliefPenaltyMPC_sub01, beliefPenaltyMPC_riub01, beliefPenaltyMPC_llb01, beliefPenaltyMPC_slb01, beliefPenaltyMPC_rilb01, beliefPenaltyMPC_lbIdx01, beliefPenaltyMPC_ubIdx01, beliefPenaltyMPC_grad_ineq01, beliefPenaltyMPC_lubbysub01, beliefPenaltyMPC_llbbyslb01);
beliefPenaltyMPC_LA_INEQ_B_GRAD_87_87_33(beliefPenaltyMPC_lub02, beliefPenaltyMPC_sub02, beliefPenaltyMPC_riub02, beliefPenaltyMPC_llb02, beliefPenaltyMPC_slb02, beliefPenaltyMPC_rilb02, beliefPenaltyMPC_lbIdx02, beliefPenaltyMPC_ubIdx02, beliefPenaltyMPC_grad_ineq02, beliefPenaltyMPC_lubbysub02, beliefPenaltyMPC_llbbyslb02);
beliefPenaltyMPC_LA_INEQ_B_GRAD_87_87_33(beliefPenaltyMPC_lub03, beliefPenaltyMPC_sub03, beliefPenaltyMPC_riub03, beliefPenaltyMPC_llb03, beliefPenaltyMPC_slb03, beliefPenaltyMPC_rilb03, beliefPenaltyMPC_lbIdx03, beliefPenaltyMPC_ubIdx03, beliefPenaltyMPC_grad_ineq03, beliefPenaltyMPC_lubbysub03, beliefPenaltyMPC_llbbyslb03);
beliefPenaltyMPC_LA_INEQ_B_GRAD_87_87_33(beliefPenaltyMPC_lub04, beliefPenaltyMPC_sub04, beliefPenaltyMPC_riub04, beliefPenaltyMPC_llb04, beliefPenaltyMPC_slb04, beliefPenaltyMPC_rilb04, beliefPenaltyMPC_lbIdx04, beliefPenaltyMPC_ubIdx04, beliefPenaltyMPC_grad_ineq04, beliefPenaltyMPC_lubbysub04, beliefPenaltyMPC_llbbyslb04);
beliefPenaltyMPC_LA_INEQ_B_GRAD_87_87_33(beliefPenaltyMPC_lub05, beliefPenaltyMPC_sub05, beliefPenaltyMPC_riub05, beliefPenaltyMPC_llb05, beliefPenaltyMPC_slb05, beliefPenaltyMPC_rilb05, beliefPenaltyMPC_lbIdx05, beliefPenaltyMPC_ubIdx05, beliefPenaltyMPC_grad_ineq05, beliefPenaltyMPC_lubbysub05, beliefPenaltyMPC_llbbyslb05);
beliefPenaltyMPC_LA_INEQ_B_GRAD_87_87_33(beliefPenaltyMPC_lub06, beliefPenaltyMPC_sub06, beliefPenaltyMPC_riub06, beliefPenaltyMPC_llb06, beliefPenaltyMPC_slb06, beliefPenaltyMPC_rilb06, beliefPenaltyMPC_lbIdx06, beliefPenaltyMPC_ubIdx06, beliefPenaltyMPC_grad_ineq06, beliefPenaltyMPC_lubbysub06, beliefPenaltyMPC_llbbyslb06);
beliefPenaltyMPC_LA_INEQ_B_GRAD_87_87_33(beliefPenaltyMPC_lub07, beliefPenaltyMPC_sub07, beliefPenaltyMPC_riub07, beliefPenaltyMPC_llb07, beliefPenaltyMPC_slb07, beliefPenaltyMPC_rilb07, beliefPenaltyMPC_lbIdx07, beliefPenaltyMPC_ubIdx07, beliefPenaltyMPC_grad_ineq07, beliefPenaltyMPC_lubbysub07, beliefPenaltyMPC_llbbyslb07);
beliefPenaltyMPC_LA_INEQ_B_GRAD_87_87_33(beliefPenaltyMPC_lub08, beliefPenaltyMPC_sub08, beliefPenaltyMPC_riub08, beliefPenaltyMPC_llb08, beliefPenaltyMPC_slb08, beliefPenaltyMPC_rilb08, beliefPenaltyMPC_lbIdx08, beliefPenaltyMPC_ubIdx08, beliefPenaltyMPC_grad_ineq08, beliefPenaltyMPC_lubbysub08, beliefPenaltyMPC_llbbyslb08);
beliefPenaltyMPC_LA_INEQ_B_GRAD_87_87_33(beliefPenaltyMPC_lub09, beliefPenaltyMPC_sub09, beliefPenaltyMPC_riub09, beliefPenaltyMPC_llb09, beliefPenaltyMPC_slb09, beliefPenaltyMPC_rilb09, beliefPenaltyMPC_lbIdx09, beliefPenaltyMPC_ubIdx09, beliefPenaltyMPC_grad_ineq09, beliefPenaltyMPC_lubbysub09, beliefPenaltyMPC_llbbyslb09);
beliefPenaltyMPC_LA_INEQ_B_GRAD_87_87_33(beliefPenaltyMPC_lub10, beliefPenaltyMPC_sub10, beliefPenaltyMPC_riub10, beliefPenaltyMPC_llb10, beliefPenaltyMPC_slb10, beliefPenaltyMPC_rilb10, beliefPenaltyMPC_lbIdx10, beliefPenaltyMPC_ubIdx10, beliefPenaltyMPC_grad_ineq10, beliefPenaltyMPC_lubbysub10, beliefPenaltyMPC_llbbyslb10);
beliefPenaltyMPC_LA_INEQ_B_GRAD_87_87_33(beliefPenaltyMPC_lub11, beliefPenaltyMPC_sub11, beliefPenaltyMPC_riub11, beliefPenaltyMPC_llb11, beliefPenaltyMPC_slb11, beliefPenaltyMPC_rilb11, beliefPenaltyMPC_lbIdx11, beliefPenaltyMPC_ubIdx11, beliefPenaltyMPC_grad_ineq11, beliefPenaltyMPC_lubbysub11, beliefPenaltyMPC_llbbyslb11);
beliefPenaltyMPC_LA_INEQ_B_GRAD_87_87_33(beliefPenaltyMPC_lub12, beliefPenaltyMPC_sub12, beliefPenaltyMPC_riub12, beliefPenaltyMPC_llb12, beliefPenaltyMPC_slb12, beliefPenaltyMPC_rilb12, beliefPenaltyMPC_lbIdx12, beliefPenaltyMPC_ubIdx12, beliefPenaltyMPC_grad_ineq12, beliefPenaltyMPC_lubbysub12, beliefPenaltyMPC_llbbyslb12);
beliefPenaltyMPC_LA_INEQ_B_GRAD_87_87_33(beliefPenaltyMPC_lub13, beliefPenaltyMPC_sub13, beliefPenaltyMPC_riub13, beliefPenaltyMPC_llb13, beliefPenaltyMPC_slb13, beliefPenaltyMPC_rilb13, beliefPenaltyMPC_lbIdx13, beliefPenaltyMPC_ubIdx13, beliefPenaltyMPC_grad_ineq13, beliefPenaltyMPC_lubbysub13, beliefPenaltyMPC_llbbyslb13);
beliefPenaltyMPC_LA_INEQ_B_GRAD_33_33_27(beliefPenaltyMPC_lub14, beliefPenaltyMPC_sub14, beliefPenaltyMPC_riub14, beliefPenaltyMPC_llb14, beliefPenaltyMPC_slb14, beliefPenaltyMPC_rilb14, beliefPenaltyMPC_lbIdx14, beliefPenaltyMPC_ubIdx14, beliefPenaltyMPC_grad_ineq14, beliefPenaltyMPC_lubbysub14, beliefPenaltyMPC_llbbyslb14);
beliefPenaltyMPC_LA_INEQ_P_6_33(params->A15, beliefPenaltyMPC_lp14, beliefPenaltyMPC_sp14, beliefPenaltyMPC_rip14, beliefPenaltyMPC_grad_ineq14, beliefPenaltyMPC_lpbysp14);
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
beliefPenaltyMPC_LA_VVADD3_1251(beliefPenaltyMPC_grad_cost, beliefPenaltyMPC_grad_eq, beliefPenaltyMPC_grad_ineq, beliefPenaltyMPC_rd);
beliefPenaltyMPC_LA_INEQ_DENSE_HESS_87_87_33(params->Q1, beliefPenaltyMPC_llbbyslb00, beliefPenaltyMPC_lbIdx00, beliefPenaltyMPC_lubbysub00, beliefPenaltyMPC_ubIdx00, beliefPenaltyMPC_Phi00);
beliefPenaltyMPC_LA_DENSE_CHOL2_87(beliefPenaltyMPC_Phi00);
beliefPenaltyMPC_LA_DENSE_MATRIXFORWARDSUB_27_87(beliefPenaltyMPC_Phi00, params->C1, beliefPenaltyMPC_V00);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MATRIXFORWARDSUB_27_87(beliefPenaltyMPC_Phi00, beliefPenaltyMPC_D00, beliefPenaltyMPC_W00);
beliefPenaltyMPC_LA_DENSE_MMTM_27_87_27(beliefPenaltyMPC_W00, beliefPenaltyMPC_V00, beliefPenaltyMPC_Ysd01);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_87(beliefPenaltyMPC_Phi00, beliefPenaltyMPC_rd00, beliefPenaltyMPC_Lbyrd00);
beliefPenaltyMPC_LA_INEQ_DENSE_HESS_87_87_33(params->Q2, beliefPenaltyMPC_llbbyslb01, beliefPenaltyMPC_lbIdx01, beliefPenaltyMPC_lubbysub01, beliefPenaltyMPC_ubIdx01, beliefPenaltyMPC_Phi01);
beliefPenaltyMPC_LA_DENSE_CHOL2_87(beliefPenaltyMPC_Phi01);
beliefPenaltyMPC_LA_DENSE_MATRIXFORWARDSUB_27_87(beliefPenaltyMPC_Phi01, params->C2, beliefPenaltyMPC_V01);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MATRIXFORWARDSUB_27_87(beliefPenaltyMPC_Phi01, beliefPenaltyMPC_D01, beliefPenaltyMPC_W01);
beliefPenaltyMPC_LA_DENSE_MMTM_27_87_27(beliefPenaltyMPC_W01, beliefPenaltyMPC_V01, beliefPenaltyMPC_Ysd02);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_87(beliefPenaltyMPC_Phi01, beliefPenaltyMPC_rd01, beliefPenaltyMPC_Lbyrd01);
beliefPenaltyMPC_LA_INEQ_DENSE_HESS_87_87_33(params->Q3, beliefPenaltyMPC_llbbyslb02, beliefPenaltyMPC_lbIdx02, beliefPenaltyMPC_lubbysub02, beliefPenaltyMPC_ubIdx02, beliefPenaltyMPC_Phi02);
beliefPenaltyMPC_LA_DENSE_CHOL2_87(beliefPenaltyMPC_Phi02);
beliefPenaltyMPC_LA_DENSE_MATRIXFORWARDSUB_27_87(beliefPenaltyMPC_Phi02, params->C3, beliefPenaltyMPC_V02);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MATRIXFORWARDSUB_27_87(beliefPenaltyMPC_Phi02, beliefPenaltyMPC_D01, beliefPenaltyMPC_W02);
beliefPenaltyMPC_LA_DENSE_MMTM_27_87_27(beliefPenaltyMPC_W02, beliefPenaltyMPC_V02, beliefPenaltyMPC_Ysd03);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_87(beliefPenaltyMPC_Phi02, beliefPenaltyMPC_rd02, beliefPenaltyMPC_Lbyrd02);
beliefPenaltyMPC_LA_INEQ_DENSE_HESS_87_87_33(params->Q4, beliefPenaltyMPC_llbbyslb03, beliefPenaltyMPC_lbIdx03, beliefPenaltyMPC_lubbysub03, beliefPenaltyMPC_ubIdx03, beliefPenaltyMPC_Phi03);
beliefPenaltyMPC_LA_DENSE_CHOL2_87(beliefPenaltyMPC_Phi03);
beliefPenaltyMPC_LA_DENSE_MATRIXFORWARDSUB_27_87(beliefPenaltyMPC_Phi03, params->C4, beliefPenaltyMPC_V03);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MATRIXFORWARDSUB_27_87(beliefPenaltyMPC_Phi03, beliefPenaltyMPC_D01, beliefPenaltyMPC_W03);
beliefPenaltyMPC_LA_DENSE_MMTM_27_87_27(beliefPenaltyMPC_W03, beliefPenaltyMPC_V03, beliefPenaltyMPC_Ysd04);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_87(beliefPenaltyMPC_Phi03, beliefPenaltyMPC_rd03, beliefPenaltyMPC_Lbyrd03);
beliefPenaltyMPC_LA_INEQ_DENSE_HESS_87_87_33(params->Q5, beliefPenaltyMPC_llbbyslb04, beliefPenaltyMPC_lbIdx04, beliefPenaltyMPC_lubbysub04, beliefPenaltyMPC_ubIdx04, beliefPenaltyMPC_Phi04);
beliefPenaltyMPC_LA_DENSE_CHOL2_87(beliefPenaltyMPC_Phi04);
beliefPenaltyMPC_LA_DENSE_MATRIXFORWARDSUB_27_87(beliefPenaltyMPC_Phi04, params->C5, beliefPenaltyMPC_V04);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MATRIXFORWARDSUB_27_87(beliefPenaltyMPC_Phi04, beliefPenaltyMPC_D01, beliefPenaltyMPC_W04);
beliefPenaltyMPC_LA_DENSE_MMTM_27_87_27(beliefPenaltyMPC_W04, beliefPenaltyMPC_V04, beliefPenaltyMPC_Ysd05);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_87(beliefPenaltyMPC_Phi04, beliefPenaltyMPC_rd04, beliefPenaltyMPC_Lbyrd04);
beliefPenaltyMPC_LA_INEQ_DENSE_HESS_87_87_33(params->Q6, beliefPenaltyMPC_llbbyslb05, beliefPenaltyMPC_lbIdx05, beliefPenaltyMPC_lubbysub05, beliefPenaltyMPC_ubIdx05, beliefPenaltyMPC_Phi05);
beliefPenaltyMPC_LA_DENSE_CHOL2_87(beliefPenaltyMPC_Phi05);
beliefPenaltyMPC_LA_DENSE_MATRIXFORWARDSUB_27_87(beliefPenaltyMPC_Phi05, params->C6, beliefPenaltyMPC_V05);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MATRIXFORWARDSUB_27_87(beliefPenaltyMPC_Phi05, beliefPenaltyMPC_D01, beliefPenaltyMPC_W05);
beliefPenaltyMPC_LA_DENSE_MMTM_27_87_27(beliefPenaltyMPC_W05, beliefPenaltyMPC_V05, beliefPenaltyMPC_Ysd06);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_87(beliefPenaltyMPC_Phi05, beliefPenaltyMPC_rd05, beliefPenaltyMPC_Lbyrd05);
beliefPenaltyMPC_LA_INEQ_DENSE_HESS_87_87_33(params->Q7, beliefPenaltyMPC_llbbyslb06, beliefPenaltyMPC_lbIdx06, beliefPenaltyMPC_lubbysub06, beliefPenaltyMPC_ubIdx06, beliefPenaltyMPC_Phi06);
beliefPenaltyMPC_LA_DENSE_CHOL2_87(beliefPenaltyMPC_Phi06);
beliefPenaltyMPC_LA_DENSE_MATRIXFORWARDSUB_27_87(beliefPenaltyMPC_Phi06, params->C7, beliefPenaltyMPC_V06);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MATRIXFORWARDSUB_27_87(beliefPenaltyMPC_Phi06, beliefPenaltyMPC_D01, beliefPenaltyMPC_W06);
beliefPenaltyMPC_LA_DENSE_MMTM_27_87_27(beliefPenaltyMPC_W06, beliefPenaltyMPC_V06, beliefPenaltyMPC_Ysd07);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_87(beliefPenaltyMPC_Phi06, beliefPenaltyMPC_rd06, beliefPenaltyMPC_Lbyrd06);
beliefPenaltyMPC_LA_INEQ_DENSE_HESS_87_87_33(params->Q8, beliefPenaltyMPC_llbbyslb07, beliefPenaltyMPC_lbIdx07, beliefPenaltyMPC_lubbysub07, beliefPenaltyMPC_ubIdx07, beliefPenaltyMPC_Phi07);
beliefPenaltyMPC_LA_DENSE_CHOL2_87(beliefPenaltyMPC_Phi07);
beliefPenaltyMPC_LA_DENSE_MATRIXFORWARDSUB_27_87(beliefPenaltyMPC_Phi07, params->C8, beliefPenaltyMPC_V07);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MATRIXFORWARDSUB_27_87(beliefPenaltyMPC_Phi07, beliefPenaltyMPC_D01, beliefPenaltyMPC_W07);
beliefPenaltyMPC_LA_DENSE_MMTM_27_87_27(beliefPenaltyMPC_W07, beliefPenaltyMPC_V07, beliefPenaltyMPC_Ysd08);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_87(beliefPenaltyMPC_Phi07, beliefPenaltyMPC_rd07, beliefPenaltyMPC_Lbyrd07);
beliefPenaltyMPC_LA_INEQ_DENSE_HESS_87_87_33(params->Q9, beliefPenaltyMPC_llbbyslb08, beliefPenaltyMPC_lbIdx08, beliefPenaltyMPC_lubbysub08, beliefPenaltyMPC_ubIdx08, beliefPenaltyMPC_Phi08);
beliefPenaltyMPC_LA_DENSE_CHOL2_87(beliefPenaltyMPC_Phi08);
beliefPenaltyMPC_LA_DENSE_MATRIXFORWARDSUB_27_87(beliefPenaltyMPC_Phi08, params->C9, beliefPenaltyMPC_V08);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MATRIXFORWARDSUB_27_87(beliefPenaltyMPC_Phi08, beliefPenaltyMPC_D01, beliefPenaltyMPC_W08);
beliefPenaltyMPC_LA_DENSE_MMTM_27_87_27(beliefPenaltyMPC_W08, beliefPenaltyMPC_V08, beliefPenaltyMPC_Ysd09);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_87(beliefPenaltyMPC_Phi08, beliefPenaltyMPC_rd08, beliefPenaltyMPC_Lbyrd08);
beliefPenaltyMPC_LA_INEQ_DENSE_HESS_87_87_33(params->Q10, beliefPenaltyMPC_llbbyslb09, beliefPenaltyMPC_lbIdx09, beliefPenaltyMPC_lubbysub09, beliefPenaltyMPC_ubIdx09, beliefPenaltyMPC_Phi09);
beliefPenaltyMPC_LA_DENSE_CHOL2_87(beliefPenaltyMPC_Phi09);
beliefPenaltyMPC_LA_DENSE_MATRIXFORWARDSUB_27_87(beliefPenaltyMPC_Phi09, params->C10, beliefPenaltyMPC_V09);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MATRIXFORWARDSUB_27_87(beliefPenaltyMPC_Phi09, beliefPenaltyMPC_D01, beliefPenaltyMPC_W09);
beliefPenaltyMPC_LA_DENSE_MMTM_27_87_27(beliefPenaltyMPC_W09, beliefPenaltyMPC_V09, beliefPenaltyMPC_Ysd10);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_87(beliefPenaltyMPC_Phi09, beliefPenaltyMPC_rd09, beliefPenaltyMPC_Lbyrd09);
beliefPenaltyMPC_LA_INEQ_DENSE_HESS_87_87_33(params->Q11, beliefPenaltyMPC_llbbyslb10, beliefPenaltyMPC_lbIdx10, beliefPenaltyMPC_lubbysub10, beliefPenaltyMPC_ubIdx10, beliefPenaltyMPC_Phi10);
beliefPenaltyMPC_LA_DENSE_CHOL2_87(beliefPenaltyMPC_Phi10);
beliefPenaltyMPC_LA_DENSE_MATRIXFORWARDSUB_27_87(beliefPenaltyMPC_Phi10, params->C11, beliefPenaltyMPC_V10);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MATRIXFORWARDSUB_27_87(beliefPenaltyMPC_Phi10, beliefPenaltyMPC_D01, beliefPenaltyMPC_W10);
beliefPenaltyMPC_LA_DENSE_MMTM_27_87_27(beliefPenaltyMPC_W10, beliefPenaltyMPC_V10, beliefPenaltyMPC_Ysd11);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_87(beliefPenaltyMPC_Phi10, beliefPenaltyMPC_rd10, beliefPenaltyMPC_Lbyrd10);
beliefPenaltyMPC_LA_INEQ_DENSE_HESS_87_87_33(params->Q12, beliefPenaltyMPC_llbbyslb11, beliefPenaltyMPC_lbIdx11, beliefPenaltyMPC_lubbysub11, beliefPenaltyMPC_ubIdx11, beliefPenaltyMPC_Phi11);
beliefPenaltyMPC_LA_DENSE_CHOL2_87(beliefPenaltyMPC_Phi11);
beliefPenaltyMPC_LA_DENSE_MATRIXFORWARDSUB_27_87(beliefPenaltyMPC_Phi11, params->C12, beliefPenaltyMPC_V11);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MATRIXFORWARDSUB_27_87(beliefPenaltyMPC_Phi11, beliefPenaltyMPC_D01, beliefPenaltyMPC_W11);
beliefPenaltyMPC_LA_DENSE_MMTM_27_87_27(beliefPenaltyMPC_W11, beliefPenaltyMPC_V11, beliefPenaltyMPC_Ysd12);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_87(beliefPenaltyMPC_Phi11, beliefPenaltyMPC_rd11, beliefPenaltyMPC_Lbyrd11);
beliefPenaltyMPC_LA_INEQ_DENSE_HESS_87_87_33(params->Q13, beliefPenaltyMPC_llbbyslb12, beliefPenaltyMPC_lbIdx12, beliefPenaltyMPC_lubbysub12, beliefPenaltyMPC_ubIdx12, beliefPenaltyMPC_Phi12);
beliefPenaltyMPC_LA_DENSE_CHOL2_87(beliefPenaltyMPC_Phi12);
beliefPenaltyMPC_LA_DENSE_MATRIXFORWARDSUB_27_87(beliefPenaltyMPC_Phi12, params->C13, beliefPenaltyMPC_V12);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MATRIXFORWARDSUB_27_87(beliefPenaltyMPC_Phi12, beliefPenaltyMPC_D01, beliefPenaltyMPC_W12);
beliefPenaltyMPC_LA_DENSE_MMTM_27_87_27(beliefPenaltyMPC_W12, beliefPenaltyMPC_V12, beliefPenaltyMPC_Ysd13);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_87(beliefPenaltyMPC_Phi12, beliefPenaltyMPC_rd12, beliefPenaltyMPC_Lbyrd12);
beliefPenaltyMPC_LA_INEQ_DENSE_HESS_87_87_33(params->Q14, beliefPenaltyMPC_llbbyslb13, beliefPenaltyMPC_lbIdx13, beliefPenaltyMPC_lubbysub13, beliefPenaltyMPC_ubIdx13, beliefPenaltyMPC_Phi13);
beliefPenaltyMPC_LA_DENSE_CHOL2_87(beliefPenaltyMPC_Phi13);
beliefPenaltyMPC_LA_DENSE_MATRIXFORWARDSUB_27_87(beliefPenaltyMPC_Phi13, params->C14, beliefPenaltyMPC_V13);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MATRIXFORWARDSUB_27_87(beliefPenaltyMPC_Phi13, beliefPenaltyMPC_D01, beliefPenaltyMPC_W13);
beliefPenaltyMPC_LA_DENSE_MMTM_27_87_27(beliefPenaltyMPC_W13, beliefPenaltyMPC_V13, beliefPenaltyMPC_Ysd14);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_87(beliefPenaltyMPC_Phi13, beliefPenaltyMPC_rd13, beliefPenaltyMPC_Lbyrd13);
beliefPenaltyMPC_LA_INEQ_DENSE_HESS_33_33_27(params->Q15, beliefPenaltyMPC_llbbyslb14, beliefPenaltyMPC_lbIdx14, beliefPenaltyMPC_lubbysub14, beliefPenaltyMPC_ubIdx14, beliefPenaltyMPC_Phi14);
beliefPenaltyMPC_LA_DENSE_ADDMTDM_6_33(params->A15, beliefPenaltyMPC_lpbysp14, beliefPenaltyMPC_Phi14);
beliefPenaltyMPC_LA_DENSE_CHOL2_33(beliefPenaltyMPC_Phi14);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MATRIXFORWARDSUB_27_33(beliefPenaltyMPC_Phi14, beliefPenaltyMPC_D14, beliefPenaltyMPC_W14);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_33(beliefPenaltyMPC_Phi14, beliefPenaltyMPC_rd14, beliefPenaltyMPC_Lbyrd14);
beliefPenaltyMPC_LA_DENSE_MMT_27_87(beliefPenaltyMPC_W00, beliefPenaltyMPC_Yd00);
beliefPenaltyMPC_LA_DENSE_MVMSUB7_27_87(beliefPenaltyMPC_W00, beliefPenaltyMPC_Lbyrd00, beliefPenaltyMPC_re00, beliefPenaltyMPC_beta00);
beliefPenaltyMPC_LA_DENSE_MMT2_27_87_87(beliefPenaltyMPC_V00, beliefPenaltyMPC_W01, beliefPenaltyMPC_Yd01);
beliefPenaltyMPC_LA_DENSE_MVMSUB2_27_87_87(beliefPenaltyMPC_V00, beliefPenaltyMPC_Lbyrd00, beliefPenaltyMPC_W01, beliefPenaltyMPC_Lbyrd01, beliefPenaltyMPC_re01, beliefPenaltyMPC_beta01);
beliefPenaltyMPC_LA_DENSE_MMT2_27_87_87(beliefPenaltyMPC_V01, beliefPenaltyMPC_W02, beliefPenaltyMPC_Yd02);
beliefPenaltyMPC_LA_DENSE_MVMSUB2_27_87_87(beliefPenaltyMPC_V01, beliefPenaltyMPC_Lbyrd01, beliefPenaltyMPC_W02, beliefPenaltyMPC_Lbyrd02, beliefPenaltyMPC_re02, beliefPenaltyMPC_beta02);
beliefPenaltyMPC_LA_DENSE_MMT2_27_87_87(beliefPenaltyMPC_V02, beliefPenaltyMPC_W03, beliefPenaltyMPC_Yd03);
beliefPenaltyMPC_LA_DENSE_MVMSUB2_27_87_87(beliefPenaltyMPC_V02, beliefPenaltyMPC_Lbyrd02, beliefPenaltyMPC_W03, beliefPenaltyMPC_Lbyrd03, beliefPenaltyMPC_re03, beliefPenaltyMPC_beta03);
beliefPenaltyMPC_LA_DENSE_MMT2_27_87_87(beliefPenaltyMPC_V03, beliefPenaltyMPC_W04, beliefPenaltyMPC_Yd04);
beliefPenaltyMPC_LA_DENSE_MVMSUB2_27_87_87(beliefPenaltyMPC_V03, beliefPenaltyMPC_Lbyrd03, beliefPenaltyMPC_W04, beliefPenaltyMPC_Lbyrd04, beliefPenaltyMPC_re04, beliefPenaltyMPC_beta04);
beliefPenaltyMPC_LA_DENSE_MMT2_27_87_87(beliefPenaltyMPC_V04, beliefPenaltyMPC_W05, beliefPenaltyMPC_Yd05);
beliefPenaltyMPC_LA_DENSE_MVMSUB2_27_87_87(beliefPenaltyMPC_V04, beliefPenaltyMPC_Lbyrd04, beliefPenaltyMPC_W05, beliefPenaltyMPC_Lbyrd05, beliefPenaltyMPC_re05, beliefPenaltyMPC_beta05);
beliefPenaltyMPC_LA_DENSE_MMT2_27_87_87(beliefPenaltyMPC_V05, beliefPenaltyMPC_W06, beliefPenaltyMPC_Yd06);
beliefPenaltyMPC_LA_DENSE_MVMSUB2_27_87_87(beliefPenaltyMPC_V05, beliefPenaltyMPC_Lbyrd05, beliefPenaltyMPC_W06, beliefPenaltyMPC_Lbyrd06, beliefPenaltyMPC_re06, beliefPenaltyMPC_beta06);
beliefPenaltyMPC_LA_DENSE_MMT2_27_87_87(beliefPenaltyMPC_V06, beliefPenaltyMPC_W07, beliefPenaltyMPC_Yd07);
beliefPenaltyMPC_LA_DENSE_MVMSUB2_27_87_87(beliefPenaltyMPC_V06, beliefPenaltyMPC_Lbyrd06, beliefPenaltyMPC_W07, beliefPenaltyMPC_Lbyrd07, beliefPenaltyMPC_re07, beliefPenaltyMPC_beta07);
beliefPenaltyMPC_LA_DENSE_MMT2_27_87_87(beliefPenaltyMPC_V07, beliefPenaltyMPC_W08, beliefPenaltyMPC_Yd08);
beliefPenaltyMPC_LA_DENSE_MVMSUB2_27_87_87(beliefPenaltyMPC_V07, beliefPenaltyMPC_Lbyrd07, beliefPenaltyMPC_W08, beliefPenaltyMPC_Lbyrd08, beliefPenaltyMPC_re08, beliefPenaltyMPC_beta08);
beliefPenaltyMPC_LA_DENSE_MMT2_27_87_87(beliefPenaltyMPC_V08, beliefPenaltyMPC_W09, beliefPenaltyMPC_Yd09);
beliefPenaltyMPC_LA_DENSE_MVMSUB2_27_87_87(beliefPenaltyMPC_V08, beliefPenaltyMPC_Lbyrd08, beliefPenaltyMPC_W09, beliefPenaltyMPC_Lbyrd09, beliefPenaltyMPC_re09, beliefPenaltyMPC_beta09);
beliefPenaltyMPC_LA_DENSE_MMT2_27_87_87(beliefPenaltyMPC_V09, beliefPenaltyMPC_W10, beliefPenaltyMPC_Yd10);
beliefPenaltyMPC_LA_DENSE_MVMSUB2_27_87_87(beliefPenaltyMPC_V09, beliefPenaltyMPC_Lbyrd09, beliefPenaltyMPC_W10, beliefPenaltyMPC_Lbyrd10, beliefPenaltyMPC_re10, beliefPenaltyMPC_beta10);
beliefPenaltyMPC_LA_DENSE_MMT2_27_87_87(beliefPenaltyMPC_V10, beliefPenaltyMPC_W11, beliefPenaltyMPC_Yd11);
beliefPenaltyMPC_LA_DENSE_MVMSUB2_27_87_87(beliefPenaltyMPC_V10, beliefPenaltyMPC_Lbyrd10, beliefPenaltyMPC_W11, beliefPenaltyMPC_Lbyrd11, beliefPenaltyMPC_re11, beliefPenaltyMPC_beta11);
beliefPenaltyMPC_LA_DENSE_MMT2_27_87_87(beliefPenaltyMPC_V11, beliefPenaltyMPC_W12, beliefPenaltyMPC_Yd12);
beliefPenaltyMPC_LA_DENSE_MVMSUB2_27_87_87(beliefPenaltyMPC_V11, beliefPenaltyMPC_Lbyrd11, beliefPenaltyMPC_W12, beliefPenaltyMPC_Lbyrd12, beliefPenaltyMPC_re12, beliefPenaltyMPC_beta12);
beliefPenaltyMPC_LA_DENSE_MMT2_27_87_87(beliefPenaltyMPC_V12, beliefPenaltyMPC_W13, beliefPenaltyMPC_Yd13);
beliefPenaltyMPC_LA_DENSE_MVMSUB2_27_87_87(beliefPenaltyMPC_V12, beliefPenaltyMPC_Lbyrd12, beliefPenaltyMPC_W13, beliefPenaltyMPC_Lbyrd13, beliefPenaltyMPC_re13, beliefPenaltyMPC_beta13);
beliefPenaltyMPC_LA_DENSE_MMT2_27_87_33(beliefPenaltyMPC_V13, beliefPenaltyMPC_W14, beliefPenaltyMPC_Yd14);
beliefPenaltyMPC_LA_DENSE_MVMSUB2_27_87_33(beliefPenaltyMPC_V13, beliefPenaltyMPC_Lbyrd13, beliefPenaltyMPC_W14, beliefPenaltyMPC_Lbyrd14, beliefPenaltyMPC_re14, beliefPenaltyMPC_beta14);
beliefPenaltyMPC_LA_DENSE_CHOL_27(beliefPenaltyMPC_Yd00, beliefPenaltyMPC_Ld00);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_27(beliefPenaltyMPC_Ld00, beliefPenaltyMPC_beta00, beliefPenaltyMPC_yy00);
beliefPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_27_27(beliefPenaltyMPC_Ld00, beliefPenaltyMPC_Ysd01, beliefPenaltyMPC_Lsd01);
beliefPenaltyMPC_LA_DENSE_MMTSUB_27_27(beliefPenaltyMPC_Lsd01, beliefPenaltyMPC_Yd01);
beliefPenaltyMPC_LA_DENSE_CHOL_27(beliefPenaltyMPC_Yd01, beliefPenaltyMPC_Ld01);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_27_27(beliefPenaltyMPC_Lsd01, beliefPenaltyMPC_yy00, beliefPenaltyMPC_beta01, beliefPenaltyMPC_bmy01);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_27(beliefPenaltyMPC_Ld01, beliefPenaltyMPC_bmy01, beliefPenaltyMPC_yy01);
beliefPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_27_27(beliefPenaltyMPC_Ld01, beliefPenaltyMPC_Ysd02, beliefPenaltyMPC_Lsd02);
beliefPenaltyMPC_LA_DENSE_MMTSUB_27_27(beliefPenaltyMPC_Lsd02, beliefPenaltyMPC_Yd02);
beliefPenaltyMPC_LA_DENSE_CHOL_27(beliefPenaltyMPC_Yd02, beliefPenaltyMPC_Ld02);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_27_27(beliefPenaltyMPC_Lsd02, beliefPenaltyMPC_yy01, beliefPenaltyMPC_beta02, beliefPenaltyMPC_bmy02);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_27(beliefPenaltyMPC_Ld02, beliefPenaltyMPC_bmy02, beliefPenaltyMPC_yy02);
beliefPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_27_27(beliefPenaltyMPC_Ld02, beliefPenaltyMPC_Ysd03, beliefPenaltyMPC_Lsd03);
beliefPenaltyMPC_LA_DENSE_MMTSUB_27_27(beliefPenaltyMPC_Lsd03, beliefPenaltyMPC_Yd03);
beliefPenaltyMPC_LA_DENSE_CHOL_27(beliefPenaltyMPC_Yd03, beliefPenaltyMPC_Ld03);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_27_27(beliefPenaltyMPC_Lsd03, beliefPenaltyMPC_yy02, beliefPenaltyMPC_beta03, beliefPenaltyMPC_bmy03);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_27(beliefPenaltyMPC_Ld03, beliefPenaltyMPC_bmy03, beliefPenaltyMPC_yy03);
beliefPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_27_27(beliefPenaltyMPC_Ld03, beliefPenaltyMPC_Ysd04, beliefPenaltyMPC_Lsd04);
beliefPenaltyMPC_LA_DENSE_MMTSUB_27_27(beliefPenaltyMPC_Lsd04, beliefPenaltyMPC_Yd04);
beliefPenaltyMPC_LA_DENSE_CHOL_27(beliefPenaltyMPC_Yd04, beliefPenaltyMPC_Ld04);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_27_27(beliefPenaltyMPC_Lsd04, beliefPenaltyMPC_yy03, beliefPenaltyMPC_beta04, beliefPenaltyMPC_bmy04);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_27(beliefPenaltyMPC_Ld04, beliefPenaltyMPC_bmy04, beliefPenaltyMPC_yy04);
beliefPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_27_27(beliefPenaltyMPC_Ld04, beliefPenaltyMPC_Ysd05, beliefPenaltyMPC_Lsd05);
beliefPenaltyMPC_LA_DENSE_MMTSUB_27_27(beliefPenaltyMPC_Lsd05, beliefPenaltyMPC_Yd05);
beliefPenaltyMPC_LA_DENSE_CHOL_27(beliefPenaltyMPC_Yd05, beliefPenaltyMPC_Ld05);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_27_27(beliefPenaltyMPC_Lsd05, beliefPenaltyMPC_yy04, beliefPenaltyMPC_beta05, beliefPenaltyMPC_bmy05);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_27(beliefPenaltyMPC_Ld05, beliefPenaltyMPC_bmy05, beliefPenaltyMPC_yy05);
beliefPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_27_27(beliefPenaltyMPC_Ld05, beliefPenaltyMPC_Ysd06, beliefPenaltyMPC_Lsd06);
beliefPenaltyMPC_LA_DENSE_MMTSUB_27_27(beliefPenaltyMPC_Lsd06, beliefPenaltyMPC_Yd06);
beliefPenaltyMPC_LA_DENSE_CHOL_27(beliefPenaltyMPC_Yd06, beliefPenaltyMPC_Ld06);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_27_27(beliefPenaltyMPC_Lsd06, beliefPenaltyMPC_yy05, beliefPenaltyMPC_beta06, beliefPenaltyMPC_bmy06);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_27(beliefPenaltyMPC_Ld06, beliefPenaltyMPC_bmy06, beliefPenaltyMPC_yy06);
beliefPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_27_27(beliefPenaltyMPC_Ld06, beliefPenaltyMPC_Ysd07, beliefPenaltyMPC_Lsd07);
beliefPenaltyMPC_LA_DENSE_MMTSUB_27_27(beliefPenaltyMPC_Lsd07, beliefPenaltyMPC_Yd07);
beliefPenaltyMPC_LA_DENSE_CHOL_27(beliefPenaltyMPC_Yd07, beliefPenaltyMPC_Ld07);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_27_27(beliefPenaltyMPC_Lsd07, beliefPenaltyMPC_yy06, beliefPenaltyMPC_beta07, beliefPenaltyMPC_bmy07);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_27(beliefPenaltyMPC_Ld07, beliefPenaltyMPC_bmy07, beliefPenaltyMPC_yy07);
beliefPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_27_27(beliefPenaltyMPC_Ld07, beliefPenaltyMPC_Ysd08, beliefPenaltyMPC_Lsd08);
beliefPenaltyMPC_LA_DENSE_MMTSUB_27_27(beliefPenaltyMPC_Lsd08, beliefPenaltyMPC_Yd08);
beliefPenaltyMPC_LA_DENSE_CHOL_27(beliefPenaltyMPC_Yd08, beliefPenaltyMPC_Ld08);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_27_27(beliefPenaltyMPC_Lsd08, beliefPenaltyMPC_yy07, beliefPenaltyMPC_beta08, beliefPenaltyMPC_bmy08);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_27(beliefPenaltyMPC_Ld08, beliefPenaltyMPC_bmy08, beliefPenaltyMPC_yy08);
beliefPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_27_27(beliefPenaltyMPC_Ld08, beliefPenaltyMPC_Ysd09, beliefPenaltyMPC_Lsd09);
beliefPenaltyMPC_LA_DENSE_MMTSUB_27_27(beliefPenaltyMPC_Lsd09, beliefPenaltyMPC_Yd09);
beliefPenaltyMPC_LA_DENSE_CHOL_27(beliefPenaltyMPC_Yd09, beliefPenaltyMPC_Ld09);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_27_27(beliefPenaltyMPC_Lsd09, beliefPenaltyMPC_yy08, beliefPenaltyMPC_beta09, beliefPenaltyMPC_bmy09);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_27(beliefPenaltyMPC_Ld09, beliefPenaltyMPC_bmy09, beliefPenaltyMPC_yy09);
beliefPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_27_27(beliefPenaltyMPC_Ld09, beliefPenaltyMPC_Ysd10, beliefPenaltyMPC_Lsd10);
beliefPenaltyMPC_LA_DENSE_MMTSUB_27_27(beliefPenaltyMPC_Lsd10, beliefPenaltyMPC_Yd10);
beliefPenaltyMPC_LA_DENSE_CHOL_27(beliefPenaltyMPC_Yd10, beliefPenaltyMPC_Ld10);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_27_27(beliefPenaltyMPC_Lsd10, beliefPenaltyMPC_yy09, beliefPenaltyMPC_beta10, beliefPenaltyMPC_bmy10);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_27(beliefPenaltyMPC_Ld10, beliefPenaltyMPC_bmy10, beliefPenaltyMPC_yy10);
beliefPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_27_27(beliefPenaltyMPC_Ld10, beliefPenaltyMPC_Ysd11, beliefPenaltyMPC_Lsd11);
beliefPenaltyMPC_LA_DENSE_MMTSUB_27_27(beliefPenaltyMPC_Lsd11, beliefPenaltyMPC_Yd11);
beliefPenaltyMPC_LA_DENSE_CHOL_27(beliefPenaltyMPC_Yd11, beliefPenaltyMPC_Ld11);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_27_27(beliefPenaltyMPC_Lsd11, beliefPenaltyMPC_yy10, beliefPenaltyMPC_beta11, beliefPenaltyMPC_bmy11);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_27(beliefPenaltyMPC_Ld11, beliefPenaltyMPC_bmy11, beliefPenaltyMPC_yy11);
beliefPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_27_27(beliefPenaltyMPC_Ld11, beliefPenaltyMPC_Ysd12, beliefPenaltyMPC_Lsd12);
beliefPenaltyMPC_LA_DENSE_MMTSUB_27_27(beliefPenaltyMPC_Lsd12, beliefPenaltyMPC_Yd12);
beliefPenaltyMPC_LA_DENSE_CHOL_27(beliefPenaltyMPC_Yd12, beliefPenaltyMPC_Ld12);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_27_27(beliefPenaltyMPC_Lsd12, beliefPenaltyMPC_yy11, beliefPenaltyMPC_beta12, beliefPenaltyMPC_bmy12);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_27(beliefPenaltyMPC_Ld12, beliefPenaltyMPC_bmy12, beliefPenaltyMPC_yy12);
beliefPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_27_27(beliefPenaltyMPC_Ld12, beliefPenaltyMPC_Ysd13, beliefPenaltyMPC_Lsd13);
beliefPenaltyMPC_LA_DENSE_MMTSUB_27_27(beliefPenaltyMPC_Lsd13, beliefPenaltyMPC_Yd13);
beliefPenaltyMPC_LA_DENSE_CHOL_27(beliefPenaltyMPC_Yd13, beliefPenaltyMPC_Ld13);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_27_27(beliefPenaltyMPC_Lsd13, beliefPenaltyMPC_yy12, beliefPenaltyMPC_beta13, beliefPenaltyMPC_bmy13);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_27(beliefPenaltyMPC_Ld13, beliefPenaltyMPC_bmy13, beliefPenaltyMPC_yy13);
beliefPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_27_27(beliefPenaltyMPC_Ld13, beliefPenaltyMPC_Ysd14, beliefPenaltyMPC_Lsd14);
beliefPenaltyMPC_LA_DENSE_MMTSUB_27_27(beliefPenaltyMPC_Lsd14, beliefPenaltyMPC_Yd14);
beliefPenaltyMPC_LA_DENSE_CHOL_27(beliefPenaltyMPC_Yd14, beliefPenaltyMPC_Ld14);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_27_27(beliefPenaltyMPC_Lsd14, beliefPenaltyMPC_yy13, beliefPenaltyMPC_beta14, beliefPenaltyMPC_bmy14);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_27(beliefPenaltyMPC_Ld14, beliefPenaltyMPC_bmy14, beliefPenaltyMPC_yy14);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_27(beliefPenaltyMPC_Ld14, beliefPenaltyMPC_yy14, beliefPenaltyMPC_dvaff14);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_27_27(beliefPenaltyMPC_Lsd14, beliefPenaltyMPC_dvaff14, beliefPenaltyMPC_yy13, beliefPenaltyMPC_bmy13);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_27(beliefPenaltyMPC_Ld13, beliefPenaltyMPC_bmy13, beliefPenaltyMPC_dvaff13);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_27_27(beliefPenaltyMPC_Lsd13, beliefPenaltyMPC_dvaff13, beliefPenaltyMPC_yy12, beliefPenaltyMPC_bmy12);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_27(beliefPenaltyMPC_Ld12, beliefPenaltyMPC_bmy12, beliefPenaltyMPC_dvaff12);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_27_27(beliefPenaltyMPC_Lsd12, beliefPenaltyMPC_dvaff12, beliefPenaltyMPC_yy11, beliefPenaltyMPC_bmy11);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_27(beliefPenaltyMPC_Ld11, beliefPenaltyMPC_bmy11, beliefPenaltyMPC_dvaff11);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_27_27(beliefPenaltyMPC_Lsd11, beliefPenaltyMPC_dvaff11, beliefPenaltyMPC_yy10, beliefPenaltyMPC_bmy10);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_27(beliefPenaltyMPC_Ld10, beliefPenaltyMPC_bmy10, beliefPenaltyMPC_dvaff10);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_27_27(beliefPenaltyMPC_Lsd10, beliefPenaltyMPC_dvaff10, beliefPenaltyMPC_yy09, beliefPenaltyMPC_bmy09);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_27(beliefPenaltyMPC_Ld09, beliefPenaltyMPC_bmy09, beliefPenaltyMPC_dvaff09);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_27_27(beliefPenaltyMPC_Lsd09, beliefPenaltyMPC_dvaff09, beliefPenaltyMPC_yy08, beliefPenaltyMPC_bmy08);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_27(beliefPenaltyMPC_Ld08, beliefPenaltyMPC_bmy08, beliefPenaltyMPC_dvaff08);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_27_27(beliefPenaltyMPC_Lsd08, beliefPenaltyMPC_dvaff08, beliefPenaltyMPC_yy07, beliefPenaltyMPC_bmy07);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_27(beliefPenaltyMPC_Ld07, beliefPenaltyMPC_bmy07, beliefPenaltyMPC_dvaff07);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_27_27(beliefPenaltyMPC_Lsd07, beliefPenaltyMPC_dvaff07, beliefPenaltyMPC_yy06, beliefPenaltyMPC_bmy06);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_27(beliefPenaltyMPC_Ld06, beliefPenaltyMPC_bmy06, beliefPenaltyMPC_dvaff06);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_27_27(beliefPenaltyMPC_Lsd06, beliefPenaltyMPC_dvaff06, beliefPenaltyMPC_yy05, beliefPenaltyMPC_bmy05);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_27(beliefPenaltyMPC_Ld05, beliefPenaltyMPC_bmy05, beliefPenaltyMPC_dvaff05);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_27_27(beliefPenaltyMPC_Lsd05, beliefPenaltyMPC_dvaff05, beliefPenaltyMPC_yy04, beliefPenaltyMPC_bmy04);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_27(beliefPenaltyMPC_Ld04, beliefPenaltyMPC_bmy04, beliefPenaltyMPC_dvaff04);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_27_27(beliefPenaltyMPC_Lsd04, beliefPenaltyMPC_dvaff04, beliefPenaltyMPC_yy03, beliefPenaltyMPC_bmy03);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_27(beliefPenaltyMPC_Ld03, beliefPenaltyMPC_bmy03, beliefPenaltyMPC_dvaff03);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_27_27(beliefPenaltyMPC_Lsd03, beliefPenaltyMPC_dvaff03, beliefPenaltyMPC_yy02, beliefPenaltyMPC_bmy02);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_27(beliefPenaltyMPC_Ld02, beliefPenaltyMPC_bmy02, beliefPenaltyMPC_dvaff02);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_27_27(beliefPenaltyMPC_Lsd02, beliefPenaltyMPC_dvaff02, beliefPenaltyMPC_yy01, beliefPenaltyMPC_bmy01);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_27(beliefPenaltyMPC_Ld01, beliefPenaltyMPC_bmy01, beliefPenaltyMPC_dvaff01);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_27_27(beliefPenaltyMPC_Lsd01, beliefPenaltyMPC_dvaff01, beliefPenaltyMPC_yy00, beliefPenaltyMPC_bmy00);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_27(beliefPenaltyMPC_Ld00, beliefPenaltyMPC_bmy00, beliefPenaltyMPC_dvaff00);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_27_87_27(params->C1, beliefPenaltyMPC_dvaff01, beliefPenaltyMPC_D00, beliefPenaltyMPC_dvaff00, beliefPenaltyMPC_grad_eq00);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_27_87_27(params->C2, beliefPenaltyMPC_dvaff02, beliefPenaltyMPC_D01, beliefPenaltyMPC_dvaff01, beliefPenaltyMPC_grad_eq01);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_27_87_27(params->C3, beliefPenaltyMPC_dvaff03, beliefPenaltyMPC_D01, beliefPenaltyMPC_dvaff02, beliefPenaltyMPC_grad_eq02);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_27_87_27(params->C4, beliefPenaltyMPC_dvaff04, beliefPenaltyMPC_D01, beliefPenaltyMPC_dvaff03, beliefPenaltyMPC_grad_eq03);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_27_87_27(params->C5, beliefPenaltyMPC_dvaff05, beliefPenaltyMPC_D01, beliefPenaltyMPC_dvaff04, beliefPenaltyMPC_grad_eq04);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_27_87_27(params->C6, beliefPenaltyMPC_dvaff06, beliefPenaltyMPC_D01, beliefPenaltyMPC_dvaff05, beliefPenaltyMPC_grad_eq05);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_27_87_27(params->C7, beliefPenaltyMPC_dvaff07, beliefPenaltyMPC_D01, beliefPenaltyMPC_dvaff06, beliefPenaltyMPC_grad_eq06);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_27_87_27(params->C8, beliefPenaltyMPC_dvaff08, beliefPenaltyMPC_D01, beliefPenaltyMPC_dvaff07, beliefPenaltyMPC_grad_eq07);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_27_87_27(params->C9, beliefPenaltyMPC_dvaff09, beliefPenaltyMPC_D01, beliefPenaltyMPC_dvaff08, beliefPenaltyMPC_grad_eq08);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_27_87_27(params->C10, beliefPenaltyMPC_dvaff10, beliefPenaltyMPC_D01, beliefPenaltyMPC_dvaff09, beliefPenaltyMPC_grad_eq09);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_27_87_27(params->C11, beliefPenaltyMPC_dvaff11, beliefPenaltyMPC_D01, beliefPenaltyMPC_dvaff10, beliefPenaltyMPC_grad_eq10);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_27_87_27(params->C12, beliefPenaltyMPC_dvaff12, beliefPenaltyMPC_D01, beliefPenaltyMPC_dvaff11, beliefPenaltyMPC_grad_eq11);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_27_87_27(params->C13, beliefPenaltyMPC_dvaff13, beliefPenaltyMPC_D01, beliefPenaltyMPC_dvaff12, beliefPenaltyMPC_grad_eq12);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_27_87_27(params->C14, beliefPenaltyMPC_dvaff14, beliefPenaltyMPC_D01, beliefPenaltyMPC_dvaff13, beliefPenaltyMPC_grad_eq13);
beliefPenaltyMPC_LA_DIAGZERO_MTVM_27_33(beliefPenaltyMPC_D14, beliefPenaltyMPC_dvaff14, beliefPenaltyMPC_grad_eq14);
beliefPenaltyMPC_LA_VSUB2_1251(beliefPenaltyMPC_rd, beliefPenaltyMPC_grad_eq, beliefPenaltyMPC_rd);
beliefPenaltyMPC_LA_DENSE_FORWARDBACKWARDSUB_87(beliefPenaltyMPC_Phi00, beliefPenaltyMPC_rd00, beliefPenaltyMPC_dzaff00);
beliefPenaltyMPC_LA_DENSE_FORWARDBACKWARDSUB_87(beliefPenaltyMPC_Phi01, beliefPenaltyMPC_rd01, beliefPenaltyMPC_dzaff01);
beliefPenaltyMPC_LA_DENSE_FORWARDBACKWARDSUB_87(beliefPenaltyMPC_Phi02, beliefPenaltyMPC_rd02, beliefPenaltyMPC_dzaff02);
beliefPenaltyMPC_LA_DENSE_FORWARDBACKWARDSUB_87(beliefPenaltyMPC_Phi03, beliefPenaltyMPC_rd03, beliefPenaltyMPC_dzaff03);
beliefPenaltyMPC_LA_DENSE_FORWARDBACKWARDSUB_87(beliefPenaltyMPC_Phi04, beliefPenaltyMPC_rd04, beliefPenaltyMPC_dzaff04);
beliefPenaltyMPC_LA_DENSE_FORWARDBACKWARDSUB_87(beliefPenaltyMPC_Phi05, beliefPenaltyMPC_rd05, beliefPenaltyMPC_dzaff05);
beliefPenaltyMPC_LA_DENSE_FORWARDBACKWARDSUB_87(beliefPenaltyMPC_Phi06, beliefPenaltyMPC_rd06, beliefPenaltyMPC_dzaff06);
beliefPenaltyMPC_LA_DENSE_FORWARDBACKWARDSUB_87(beliefPenaltyMPC_Phi07, beliefPenaltyMPC_rd07, beliefPenaltyMPC_dzaff07);
beliefPenaltyMPC_LA_DENSE_FORWARDBACKWARDSUB_87(beliefPenaltyMPC_Phi08, beliefPenaltyMPC_rd08, beliefPenaltyMPC_dzaff08);
beliefPenaltyMPC_LA_DENSE_FORWARDBACKWARDSUB_87(beliefPenaltyMPC_Phi09, beliefPenaltyMPC_rd09, beliefPenaltyMPC_dzaff09);
beliefPenaltyMPC_LA_DENSE_FORWARDBACKWARDSUB_87(beliefPenaltyMPC_Phi10, beliefPenaltyMPC_rd10, beliefPenaltyMPC_dzaff10);
beliefPenaltyMPC_LA_DENSE_FORWARDBACKWARDSUB_87(beliefPenaltyMPC_Phi11, beliefPenaltyMPC_rd11, beliefPenaltyMPC_dzaff11);
beliefPenaltyMPC_LA_DENSE_FORWARDBACKWARDSUB_87(beliefPenaltyMPC_Phi12, beliefPenaltyMPC_rd12, beliefPenaltyMPC_dzaff12);
beliefPenaltyMPC_LA_DENSE_FORWARDBACKWARDSUB_87(beliefPenaltyMPC_Phi13, beliefPenaltyMPC_rd13, beliefPenaltyMPC_dzaff13);
beliefPenaltyMPC_LA_DENSE_FORWARDBACKWARDSUB_33(beliefPenaltyMPC_Phi14, beliefPenaltyMPC_rd14, beliefPenaltyMPC_dzaff14);
beliefPenaltyMPC_LA_VSUB_INDEXED_87(beliefPenaltyMPC_dzaff00, beliefPenaltyMPC_lbIdx00, beliefPenaltyMPC_rilb00, beliefPenaltyMPC_dslbaff00);
beliefPenaltyMPC_LA_VSUB3_87(beliefPenaltyMPC_llbbyslb00, beliefPenaltyMPC_dslbaff00, beliefPenaltyMPC_llb00, beliefPenaltyMPC_dllbaff00);
beliefPenaltyMPC_LA_VSUB2_INDEXED_33(beliefPenaltyMPC_riub00, beliefPenaltyMPC_dzaff00, beliefPenaltyMPC_ubIdx00, beliefPenaltyMPC_dsubaff00);
beliefPenaltyMPC_LA_VSUB3_33(beliefPenaltyMPC_lubbysub00, beliefPenaltyMPC_dsubaff00, beliefPenaltyMPC_lub00, beliefPenaltyMPC_dlubaff00);
beliefPenaltyMPC_LA_VSUB_INDEXED_87(beliefPenaltyMPC_dzaff01, beliefPenaltyMPC_lbIdx01, beliefPenaltyMPC_rilb01, beliefPenaltyMPC_dslbaff01);
beliefPenaltyMPC_LA_VSUB3_87(beliefPenaltyMPC_llbbyslb01, beliefPenaltyMPC_dslbaff01, beliefPenaltyMPC_llb01, beliefPenaltyMPC_dllbaff01);
beliefPenaltyMPC_LA_VSUB2_INDEXED_33(beliefPenaltyMPC_riub01, beliefPenaltyMPC_dzaff01, beliefPenaltyMPC_ubIdx01, beliefPenaltyMPC_dsubaff01);
beliefPenaltyMPC_LA_VSUB3_33(beliefPenaltyMPC_lubbysub01, beliefPenaltyMPC_dsubaff01, beliefPenaltyMPC_lub01, beliefPenaltyMPC_dlubaff01);
beliefPenaltyMPC_LA_VSUB_INDEXED_87(beliefPenaltyMPC_dzaff02, beliefPenaltyMPC_lbIdx02, beliefPenaltyMPC_rilb02, beliefPenaltyMPC_dslbaff02);
beliefPenaltyMPC_LA_VSUB3_87(beliefPenaltyMPC_llbbyslb02, beliefPenaltyMPC_dslbaff02, beliefPenaltyMPC_llb02, beliefPenaltyMPC_dllbaff02);
beliefPenaltyMPC_LA_VSUB2_INDEXED_33(beliefPenaltyMPC_riub02, beliefPenaltyMPC_dzaff02, beliefPenaltyMPC_ubIdx02, beliefPenaltyMPC_dsubaff02);
beliefPenaltyMPC_LA_VSUB3_33(beliefPenaltyMPC_lubbysub02, beliefPenaltyMPC_dsubaff02, beliefPenaltyMPC_lub02, beliefPenaltyMPC_dlubaff02);
beliefPenaltyMPC_LA_VSUB_INDEXED_87(beliefPenaltyMPC_dzaff03, beliefPenaltyMPC_lbIdx03, beliefPenaltyMPC_rilb03, beliefPenaltyMPC_dslbaff03);
beliefPenaltyMPC_LA_VSUB3_87(beliefPenaltyMPC_llbbyslb03, beliefPenaltyMPC_dslbaff03, beliefPenaltyMPC_llb03, beliefPenaltyMPC_dllbaff03);
beliefPenaltyMPC_LA_VSUB2_INDEXED_33(beliefPenaltyMPC_riub03, beliefPenaltyMPC_dzaff03, beliefPenaltyMPC_ubIdx03, beliefPenaltyMPC_dsubaff03);
beliefPenaltyMPC_LA_VSUB3_33(beliefPenaltyMPC_lubbysub03, beliefPenaltyMPC_dsubaff03, beliefPenaltyMPC_lub03, beliefPenaltyMPC_dlubaff03);
beliefPenaltyMPC_LA_VSUB_INDEXED_87(beliefPenaltyMPC_dzaff04, beliefPenaltyMPC_lbIdx04, beliefPenaltyMPC_rilb04, beliefPenaltyMPC_dslbaff04);
beliefPenaltyMPC_LA_VSUB3_87(beliefPenaltyMPC_llbbyslb04, beliefPenaltyMPC_dslbaff04, beliefPenaltyMPC_llb04, beliefPenaltyMPC_dllbaff04);
beliefPenaltyMPC_LA_VSUB2_INDEXED_33(beliefPenaltyMPC_riub04, beliefPenaltyMPC_dzaff04, beliefPenaltyMPC_ubIdx04, beliefPenaltyMPC_dsubaff04);
beliefPenaltyMPC_LA_VSUB3_33(beliefPenaltyMPC_lubbysub04, beliefPenaltyMPC_dsubaff04, beliefPenaltyMPC_lub04, beliefPenaltyMPC_dlubaff04);
beliefPenaltyMPC_LA_VSUB_INDEXED_87(beliefPenaltyMPC_dzaff05, beliefPenaltyMPC_lbIdx05, beliefPenaltyMPC_rilb05, beliefPenaltyMPC_dslbaff05);
beliefPenaltyMPC_LA_VSUB3_87(beliefPenaltyMPC_llbbyslb05, beliefPenaltyMPC_dslbaff05, beliefPenaltyMPC_llb05, beliefPenaltyMPC_dllbaff05);
beliefPenaltyMPC_LA_VSUB2_INDEXED_33(beliefPenaltyMPC_riub05, beliefPenaltyMPC_dzaff05, beliefPenaltyMPC_ubIdx05, beliefPenaltyMPC_dsubaff05);
beliefPenaltyMPC_LA_VSUB3_33(beliefPenaltyMPC_lubbysub05, beliefPenaltyMPC_dsubaff05, beliefPenaltyMPC_lub05, beliefPenaltyMPC_dlubaff05);
beliefPenaltyMPC_LA_VSUB_INDEXED_87(beliefPenaltyMPC_dzaff06, beliefPenaltyMPC_lbIdx06, beliefPenaltyMPC_rilb06, beliefPenaltyMPC_dslbaff06);
beliefPenaltyMPC_LA_VSUB3_87(beliefPenaltyMPC_llbbyslb06, beliefPenaltyMPC_dslbaff06, beliefPenaltyMPC_llb06, beliefPenaltyMPC_dllbaff06);
beliefPenaltyMPC_LA_VSUB2_INDEXED_33(beliefPenaltyMPC_riub06, beliefPenaltyMPC_dzaff06, beliefPenaltyMPC_ubIdx06, beliefPenaltyMPC_dsubaff06);
beliefPenaltyMPC_LA_VSUB3_33(beliefPenaltyMPC_lubbysub06, beliefPenaltyMPC_dsubaff06, beliefPenaltyMPC_lub06, beliefPenaltyMPC_dlubaff06);
beliefPenaltyMPC_LA_VSUB_INDEXED_87(beliefPenaltyMPC_dzaff07, beliefPenaltyMPC_lbIdx07, beliefPenaltyMPC_rilb07, beliefPenaltyMPC_dslbaff07);
beliefPenaltyMPC_LA_VSUB3_87(beliefPenaltyMPC_llbbyslb07, beliefPenaltyMPC_dslbaff07, beliefPenaltyMPC_llb07, beliefPenaltyMPC_dllbaff07);
beliefPenaltyMPC_LA_VSUB2_INDEXED_33(beliefPenaltyMPC_riub07, beliefPenaltyMPC_dzaff07, beliefPenaltyMPC_ubIdx07, beliefPenaltyMPC_dsubaff07);
beliefPenaltyMPC_LA_VSUB3_33(beliefPenaltyMPC_lubbysub07, beliefPenaltyMPC_dsubaff07, beliefPenaltyMPC_lub07, beliefPenaltyMPC_dlubaff07);
beliefPenaltyMPC_LA_VSUB_INDEXED_87(beliefPenaltyMPC_dzaff08, beliefPenaltyMPC_lbIdx08, beliefPenaltyMPC_rilb08, beliefPenaltyMPC_dslbaff08);
beliefPenaltyMPC_LA_VSUB3_87(beliefPenaltyMPC_llbbyslb08, beliefPenaltyMPC_dslbaff08, beliefPenaltyMPC_llb08, beliefPenaltyMPC_dllbaff08);
beliefPenaltyMPC_LA_VSUB2_INDEXED_33(beliefPenaltyMPC_riub08, beliefPenaltyMPC_dzaff08, beliefPenaltyMPC_ubIdx08, beliefPenaltyMPC_dsubaff08);
beliefPenaltyMPC_LA_VSUB3_33(beliefPenaltyMPC_lubbysub08, beliefPenaltyMPC_dsubaff08, beliefPenaltyMPC_lub08, beliefPenaltyMPC_dlubaff08);
beliefPenaltyMPC_LA_VSUB_INDEXED_87(beliefPenaltyMPC_dzaff09, beliefPenaltyMPC_lbIdx09, beliefPenaltyMPC_rilb09, beliefPenaltyMPC_dslbaff09);
beliefPenaltyMPC_LA_VSUB3_87(beliefPenaltyMPC_llbbyslb09, beliefPenaltyMPC_dslbaff09, beliefPenaltyMPC_llb09, beliefPenaltyMPC_dllbaff09);
beliefPenaltyMPC_LA_VSUB2_INDEXED_33(beliefPenaltyMPC_riub09, beliefPenaltyMPC_dzaff09, beliefPenaltyMPC_ubIdx09, beliefPenaltyMPC_dsubaff09);
beliefPenaltyMPC_LA_VSUB3_33(beliefPenaltyMPC_lubbysub09, beliefPenaltyMPC_dsubaff09, beliefPenaltyMPC_lub09, beliefPenaltyMPC_dlubaff09);
beliefPenaltyMPC_LA_VSUB_INDEXED_87(beliefPenaltyMPC_dzaff10, beliefPenaltyMPC_lbIdx10, beliefPenaltyMPC_rilb10, beliefPenaltyMPC_dslbaff10);
beliefPenaltyMPC_LA_VSUB3_87(beliefPenaltyMPC_llbbyslb10, beliefPenaltyMPC_dslbaff10, beliefPenaltyMPC_llb10, beliefPenaltyMPC_dllbaff10);
beliefPenaltyMPC_LA_VSUB2_INDEXED_33(beliefPenaltyMPC_riub10, beliefPenaltyMPC_dzaff10, beliefPenaltyMPC_ubIdx10, beliefPenaltyMPC_dsubaff10);
beliefPenaltyMPC_LA_VSUB3_33(beliefPenaltyMPC_lubbysub10, beliefPenaltyMPC_dsubaff10, beliefPenaltyMPC_lub10, beliefPenaltyMPC_dlubaff10);
beliefPenaltyMPC_LA_VSUB_INDEXED_87(beliefPenaltyMPC_dzaff11, beliefPenaltyMPC_lbIdx11, beliefPenaltyMPC_rilb11, beliefPenaltyMPC_dslbaff11);
beliefPenaltyMPC_LA_VSUB3_87(beliefPenaltyMPC_llbbyslb11, beliefPenaltyMPC_dslbaff11, beliefPenaltyMPC_llb11, beliefPenaltyMPC_dllbaff11);
beliefPenaltyMPC_LA_VSUB2_INDEXED_33(beliefPenaltyMPC_riub11, beliefPenaltyMPC_dzaff11, beliefPenaltyMPC_ubIdx11, beliefPenaltyMPC_dsubaff11);
beliefPenaltyMPC_LA_VSUB3_33(beliefPenaltyMPC_lubbysub11, beliefPenaltyMPC_dsubaff11, beliefPenaltyMPC_lub11, beliefPenaltyMPC_dlubaff11);
beliefPenaltyMPC_LA_VSUB_INDEXED_87(beliefPenaltyMPC_dzaff12, beliefPenaltyMPC_lbIdx12, beliefPenaltyMPC_rilb12, beliefPenaltyMPC_dslbaff12);
beliefPenaltyMPC_LA_VSUB3_87(beliefPenaltyMPC_llbbyslb12, beliefPenaltyMPC_dslbaff12, beliefPenaltyMPC_llb12, beliefPenaltyMPC_dllbaff12);
beliefPenaltyMPC_LA_VSUB2_INDEXED_33(beliefPenaltyMPC_riub12, beliefPenaltyMPC_dzaff12, beliefPenaltyMPC_ubIdx12, beliefPenaltyMPC_dsubaff12);
beliefPenaltyMPC_LA_VSUB3_33(beliefPenaltyMPC_lubbysub12, beliefPenaltyMPC_dsubaff12, beliefPenaltyMPC_lub12, beliefPenaltyMPC_dlubaff12);
beliefPenaltyMPC_LA_VSUB_INDEXED_87(beliefPenaltyMPC_dzaff13, beliefPenaltyMPC_lbIdx13, beliefPenaltyMPC_rilb13, beliefPenaltyMPC_dslbaff13);
beliefPenaltyMPC_LA_VSUB3_87(beliefPenaltyMPC_llbbyslb13, beliefPenaltyMPC_dslbaff13, beliefPenaltyMPC_llb13, beliefPenaltyMPC_dllbaff13);
beliefPenaltyMPC_LA_VSUB2_INDEXED_33(beliefPenaltyMPC_riub13, beliefPenaltyMPC_dzaff13, beliefPenaltyMPC_ubIdx13, beliefPenaltyMPC_dsubaff13);
beliefPenaltyMPC_LA_VSUB3_33(beliefPenaltyMPC_lubbysub13, beliefPenaltyMPC_dsubaff13, beliefPenaltyMPC_lub13, beliefPenaltyMPC_dlubaff13);
beliefPenaltyMPC_LA_VSUB_INDEXED_33(beliefPenaltyMPC_dzaff14, beliefPenaltyMPC_lbIdx14, beliefPenaltyMPC_rilb14, beliefPenaltyMPC_dslbaff14);
beliefPenaltyMPC_LA_VSUB3_33(beliefPenaltyMPC_llbbyslb14, beliefPenaltyMPC_dslbaff14, beliefPenaltyMPC_llb14, beliefPenaltyMPC_dllbaff14);
beliefPenaltyMPC_LA_VSUB2_INDEXED_27(beliefPenaltyMPC_riub14, beliefPenaltyMPC_dzaff14, beliefPenaltyMPC_ubIdx14, beliefPenaltyMPC_dsubaff14);
beliefPenaltyMPC_LA_VSUB3_27(beliefPenaltyMPC_lubbysub14, beliefPenaltyMPC_dsubaff14, beliefPenaltyMPC_lub14, beliefPenaltyMPC_dlubaff14);
beliefPenaltyMPC_LA_DENSE_MVMSUB4_6_33(params->A15, beliefPenaltyMPC_dzaff14, beliefPenaltyMPC_rip14, beliefPenaltyMPC_dsp_aff14);
beliefPenaltyMPC_LA_VSUB3_6(beliefPenaltyMPC_lpbysp14, beliefPenaltyMPC_dsp_aff14, beliefPenaltyMPC_lp14, beliefPenaltyMPC_dlp_aff14);
info->lsit_aff = beliefPenaltyMPC_LINESEARCH_BACKTRACKING_AFFINE(beliefPenaltyMPC_l, beliefPenaltyMPC_s, beliefPenaltyMPC_dl_aff, beliefPenaltyMPC_ds_aff, &info->step_aff, &info->mu_aff);
if( info->lsit_aff == beliefPenaltyMPC_NOPROGRESS ){
exitcode = beliefPenaltyMPC_NOPROGRESS; break;
}
sigma_3rdroot = info->mu_aff / info->mu;
info->sigma = sigma_3rdroot*sigma_3rdroot*sigma_3rdroot;
musigma = info->mu * info->sigma;
beliefPenaltyMPC_LA_VSUB5_1746(beliefPenaltyMPC_ds_aff, beliefPenaltyMPC_dl_aff, musigma, beliefPenaltyMPC_ccrhs);
beliefPenaltyMPC_LA_VSUB6_INDEXED_87_33_87(beliefPenaltyMPC_ccrhsub00, beliefPenaltyMPC_sub00, beliefPenaltyMPC_ubIdx00, beliefPenaltyMPC_ccrhsl00, beliefPenaltyMPC_slb00, beliefPenaltyMPC_lbIdx00, beliefPenaltyMPC_rd00);
beliefPenaltyMPC_LA_VSUB6_INDEXED_87_33_87(beliefPenaltyMPC_ccrhsub01, beliefPenaltyMPC_sub01, beliefPenaltyMPC_ubIdx01, beliefPenaltyMPC_ccrhsl01, beliefPenaltyMPC_slb01, beliefPenaltyMPC_lbIdx01, beliefPenaltyMPC_rd01);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_87(beliefPenaltyMPC_Phi00, beliefPenaltyMPC_rd00, beliefPenaltyMPC_Lbyrd00);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_87(beliefPenaltyMPC_Phi01, beliefPenaltyMPC_rd01, beliefPenaltyMPC_Lbyrd01);
beliefPenaltyMPC_LA_DENSE_MVM_27_87(beliefPenaltyMPC_W00, beliefPenaltyMPC_Lbyrd00, beliefPenaltyMPC_beta00);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_27(beliefPenaltyMPC_Ld00, beliefPenaltyMPC_beta00, beliefPenaltyMPC_yy00);
beliefPenaltyMPC_LA_DENSE_2MVMADD_27_87_87(beliefPenaltyMPC_V00, beliefPenaltyMPC_Lbyrd00, beliefPenaltyMPC_W01, beliefPenaltyMPC_Lbyrd01, beliefPenaltyMPC_beta01);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_27_27(beliefPenaltyMPC_Lsd01, beliefPenaltyMPC_yy00, beliefPenaltyMPC_beta01, beliefPenaltyMPC_bmy01);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_27(beliefPenaltyMPC_Ld01, beliefPenaltyMPC_bmy01, beliefPenaltyMPC_yy01);
beliefPenaltyMPC_LA_VSUB6_INDEXED_87_33_87(beliefPenaltyMPC_ccrhsub02, beliefPenaltyMPC_sub02, beliefPenaltyMPC_ubIdx02, beliefPenaltyMPC_ccrhsl02, beliefPenaltyMPC_slb02, beliefPenaltyMPC_lbIdx02, beliefPenaltyMPC_rd02);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_87(beliefPenaltyMPC_Phi02, beliefPenaltyMPC_rd02, beliefPenaltyMPC_Lbyrd02);
beliefPenaltyMPC_LA_DENSE_2MVMADD_27_87_87(beliefPenaltyMPC_V01, beliefPenaltyMPC_Lbyrd01, beliefPenaltyMPC_W02, beliefPenaltyMPC_Lbyrd02, beliefPenaltyMPC_beta02);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_27_27(beliefPenaltyMPC_Lsd02, beliefPenaltyMPC_yy01, beliefPenaltyMPC_beta02, beliefPenaltyMPC_bmy02);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_27(beliefPenaltyMPC_Ld02, beliefPenaltyMPC_bmy02, beliefPenaltyMPC_yy02);
beliefPenaltyMPC_LA_VSUB6_INDEXED_87_33_87(beliefPenaltyMPC_ccrhsub03, beliefPenaltyMPC_sub03, beliefPenaltyMPC_ubIdx03, beliefPenaltyMPC_ccrhsl03, beliefPenaltyMPC_slb03, beliefPenaltyMPC_lbIdx03, beliefPenaltyMPC_rd03);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_87(beliefPenaltyMPC_Phi03, beliefPenaltyMPC_rd03, beliefPenaltyMPC_Lbyrd03);
beliefPenaltyMPC_LA_DENSE_2MVMADD_27_87_87(beliefPenaltyMPC_V02, beliefPenaltyMPC_Lbyrd02, beliefPenaltyMPC_W03, beliefPenaltyMPC_Lbyrd03, beliefPenaltyMPC_beta03);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_27_27(beliefPenaltyMPC_Lsd03, beliefPenaltyMPC_yy02, beliefPenaltyMPC_beta03, beliefPenaltyMPC_bmy03);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_27(beliefPenaltyMPC_Ld03, beliefPenaltyMPC_bmy03, beliefPenaltyMPC_yy03);
beliefPenaltyMPC_LA_VSUB6_INDEXED_87_33_87(beliefPenaltyMPC_ccrhsub04, beliefPenaltyMPC_sub04, beliefPenaltyMPC_ubIdx04, beliefPenaltyMPC_ccrhsl04, beliefPenaltyMPC_slb04, beliefPenaltyMPC_lbIdx04, beliefPenaltyMPC_rd04);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_87(beliefPenaltyMPC_Phi04, beliefPenaltyMPC_rd04, beliefPenaltyMPC_Lbyrd04);
beliefPenaltyMPC_LA_DENSE_2MVMADD_27_87_87(beliefPenaltyMPC_V03, beliefPenaltyMPC_Lbyrd03, beliefPenaltyMPC_W04, beliefPenaltyMPC_Lbyrd04, beliefPenaltyMPC_beta04);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_27_27(beliefPenaltyMPC_Lsd04, beliefPenaltyMPC_yy03, beliefPenaltyMPC_beta04, beliefPenaltyMPC_bmy04);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_27(beliefPenaltyMPC_Ld04, beliefPenaltyMPC_bmy04, beliefPenaltyMPC_yy04);
beliefPenaltyMPC_LA_VSUB6_INDEXED_87_33_87(beliefPenaltyMPC_ccrhsub05, beliefPenaltyMPC_sub05, beliefPenaltyMPC_ubIdx05, beliefPenaltyMPC_ccrhsl05, beliefPenaltyMPC_slb05, beliefPenaltyMPC_lbIdx05, beliefPenaltyMPC_rd05);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_87(beliefPenaltyMPC_Phi05, beliefPenaltyMPC_rd05, beliefPenaltyMPC_Lbyrd05);
beliefPenaltyMPC_LA_DENSE_2MVMADD_27_87_87(beliefPenaltyMPC_V04, beliefPenaltyMPC_Lbyrd04, beliefPenaltyMPC_W05, beliefPenaltyMPC_Lbyrd05, beliefPenaltyMPC_beta05);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_27_27(beliefPenaltyMPC_Lsd05, beliefPenaltyMPC_yy04, beliefPenaltyMPC_beta05, beliefPenaltyMPC_bmy05);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_27(beliefPenaltyMPC_Ld05, beliefPenaltyMPC_bmy05, beliefPenaltyMPC_yy05);
beliefPenaltyMPC_LA_VSUB6_INDEXED_87_33_87(beliefPenaltyMPC_ccrhsub06, beliefPenaltyMPC_sub06, beliefPenaltyMPC_ubIdx06, beliefPenaltyMPC_ccrhsl06, beliefPenaltyMPC_slb06, beliefPenaltyMPC_lbIdx06, beliefPenaltyMPC_rd06);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_87(beliefPenaltyMPC_Phi06, beliefPenaltyMPC_rd06, beliefPenaltyMPC_Lbyrd06);
beliefPenaltyMPC_LA_DENSE_2MVMADD_27_87_87(beliefPenaltyMPC_V05, beliefPenaltyMPC_Lbyrd05, beliefPenaltyMPC_W06, beliefPenaltyMPC_Lbyrd06, beliefPenaltyMPC_beta06);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_27_27(beliefPenaltyMPC_Lsd06, beliefPenaltyMPC_yy05, beliefPenaltyMPC_beta06, beliefPenaltyMPC_bmy06);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_27(beliefPenaltyMPC_Ld06, beliefPenaltyMPC_bmy06, beliefPenaltyMPC_yy06);
beliefPenaltyMPC_LA_VSUB6_INDEXED_87_33_87(beliefPenaltyMPC_ccrhsub07, beliefPenaltyMPC_sub07, beliefPenaltyMPC_ubIdx07, beliefPenaltyMPC_ccrhsl07, beliefPenaltyMPC_slb07, beliefPenaltyMPC_lbIdx07, beliefPenaltyMPC_rd07);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_87(beliefPenaltyMPC_Phi07, beliefPenaltyMPC_rd07, beliefPenaltyMPC_Lbyrd07);
beliefPenaltyMPC_LA_DENSE_2MVMADD_27_87_87(beliefPenaltyMPC_V06, beliefPenaltyMPC_Lbyrd06, beliefPenaltyMPC_W07, beliefPenaltyMPC_Lbyrd07, beliefPenaltyMPC_beta07);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_27_27(beliefPenaltyMPC_Lsd07, beliefPenaltyMPC_yy06, beliefPenaltyMPC_beta07, beliefPenaltyMPC_bmy07);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_27(beliefPenaltyMPC_Ld07, beliefPenaltyMPC_bmy07, beliefPenaltyMPC_yy07);
beliefPenaltyMPC_LA_VSUB6_INDEXED_87_33_87(beliefPenaltyMPC_ccrhsub08, beliefPenaltyMPC_sub08, beliefPenaltyMPC_ubIdx08, beliefPenaltyMPC_ccrhsl08, beliefPenaltyMPC_slb08, beliefPenaltyMPC_lbIdx08, beliefPenaltyMPC_rd08);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_87(beliefPenaltyMPC_Phi08, beliefPenaltyMPC_rd08, beliefPenaltyMPC_Lbyrd08);
beliefPenaltyMPC_LA_DENSE_2MVMADD_27_87_87(beliefPenaltyMPC_V07, beliefPenaltyMPC_Lbyrd07, beliefPenaltyMPC_W08, beliefPenaltyMPC_Lbyrd08, beliefPenaltyMPC_beta08);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_27_27(beliefPenaltyMPC_Lsd08, beliefPenaltyMPC_yy07, beliefPenaltyMPC_beta08, beliefPenaltyMPC_bmy08);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_27(beliefPenaltyMPC_Ld08, beliefPenaltyMPC_bmy08, beliefPenaltyMPC_yy08);
beliefPenaltyMPC_LA_VSUB6_INDEXED_87_33_87(beliefPenaltyMPC_ccrhsub09, beliefPenaltyMPC_sub09, beliefPenaltyMPC_ubIdx09, beliefPenaltyMPC_ccrhsl09, beliefPenaltyMPC_slb09, beliefPenaltyMPC_lbIdx09, beliefPenaltyMPC_rd09);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_87(beliefPenaltyMPC_Phi09, beliefPenaltyMPC_rd09, beliefPenaltyMPC_Lbyrd09);
beliefPenaltyMPC_LA_DENSE_2MVMADD_27_87_87(beliefPenaltyMPC_V08, beliefPenaltyMPC_Lbyrd08, beliefPenaltyMPC_W09, beliefPenaltyMPC_Lbyrd09, beliefPenaltyMPC_beta09);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_27_27(beliefPenaltyMPC_Lsd09, beliefPenaltyMPC_yy08, beliefPenaltyMPC_beta09, beliefPenaltyMPC_bmy09);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_27(beliefPenaltyMPC_Ld09, beliefPenaltyMPC_bmy09, beliefPenaltyMPC_yy09);
beliefPenaltyMPC_LA_VSUB6_INDEXED_87_33_87(beliefPenaltyMPC_ccrhsub10, beliefPenaltyMPC_sub10, beliefPenaltyMPC_ubIdx10, beliefPenaltyMPC_ccrhsl10, beliefPenaltyMPC_slb10, beliefPenaltyMPC_lbIdx10, beliefPenaltyMPC_rd10);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_87(beliefPenaltyMPC_Phi10, beliefPenaltyMPC_rd10, beliefPenaltyMPC_Lbyrd10);
beliefPenaltyMPC_LA_DENSE_2MVMADD_27_87_87(beliefPenaltyMPC_V09, beliefPenaltyMPC_Lbyrd09, beliefPenaltyMPC_W10, beliefPenaltyMPC_Lbyrd10, beliefPenaltyMPC_beta10);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_27_27(beliefPenaltyMPC_Lsd10, beliefPenaltyMPC_yy09, beliefPenaltyMPC_beta10, beliefPenaltyMPC_bmy10);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_27(beliefPenaltyMPC_Ld10, beliefPenaltyMPC_bmy10, beliefPenaltyMPC_yy10);
beliefPenaltyMPC_LA_VSUB6_INDEXED_87_33_87(beliefPenaltyMPC_ccrhsub11, beliefPenaltyMPC_sub11, beliefPenaltyMPC_ubIdx11, beliefPenaltyMPC_ccrhsl11, beliefPenaltyMPC_slb11, beliefPenaltyMPC_lbIdx11, beliefPenaltyMPC_rd11);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_87(beliefPenaltyMPC_Phi11, beliefPenaltyMPC_rd11, beliefPenaltyMPC_Lbyrd11);
beliefPenaltyMPC_LA_DENSE_2MVMADD_27_87_87(beliefPenaltyMPC_V10, beliefPenaltyMPC_Lbyrd10, beliefPenaltyMPC_W11, beliefPenaltyMPC_Lbyrd11, beliefPenaltyMPC_beta11);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_27_27(beliefPenaltyMPC_Lsd11, beliefPenaltyMPC_yy10, beliefPenaltyMPC_beta11, beliefPenaltyMPC_bmy11);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_27(beliefPenaltyMPC_Ld11, beliefPenaltyMPC_bmy11, beliefPenaltyMPC_yy11);
beliefPenaltyMPC_LA_VSUB6_INDEXED_87_33_87(beliefPenaltyMPC_ccrhsub12, beliefPenaltyMPC_sub12, beliefPenaltyMPC_ubIdx12, beliefPenaltyMPC_ccrhsl12, beliefPenaltyMPC_slb12, beliefPenaltyMPC_lbIdx12, beliefPenaltyMPC_rd12);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_87(beliefPenaltyMPC_Phi12, beliefPenaltyMPC_rd12, beliefPenaltyMPC_Lbyrd12);
beliefPenaltyMPC_LA_DENSE_2MVMADD_27_87_87(beliefPenaltyMPC_V11, beliefPenaltyMPC_Lbyrd11, beliefPenaltyMPC_W12, beliefPenaltyMPC_Lbyrd12, beliefPenaltyMPC_beta12);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_27_27(beliefPenaltyMPC_Lsd12, beliefPenaltyMPC_yy11, beliefPenaltyMPC_beta12, beliefPenaltyMPC_bmy12);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_27(beliefPenaltyMPC_Ld12, beliefPenaltyMPC_bmy12, beliefPenaltyMPC_yy12);
beliefPenaltyMPC_LA_VSUB6_INDEXED_87_33_87(beliefPenaltyMPC_ccrhsub13, beliefPenaltyMPC_sub13, beliefPenaltyMPC_ubIdx13, beliefPenaltyMPC_ccrhsl13, beliefPenaltyMPC_slb13, beliefPenaltyMPC_lbIdx13, beliefPenaltyMPC_rd13);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_87(beliefPenaltyMPC_Phi13, beliefPenaltyMPC_rd13, beliefPenaltyMPC_Lbyrd13);
beliefPenaltyMPC_LA_DENSE_2MVMADD_27_87_87(beliefPenaltyMPC_V12, beliefPenaltyMPC_Lbyrd12, beliefPenaltyMPC_W13, beliefPenaltyMPC_Lbyrd13, beliefPenaltyMPC_beta13);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_27_27(beliefPenaltyMPC_Lsd13, beliefPenaltyMPC_yy12, beliefPenaltyMPC_beta13, beliefPenaltyMPC_bmy13);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_27(beliefPenaltyMPC_Ld13, beliefPenaltyMPC_bmy13, beliefPenaltyMPC_yy13);
beliefPenaltyMPC_LA_VSUB6_INDEXED_33_27_33(beliefPenaltyMPC_ccrhsub14, beliefPenaltyMPC_sub14, beliefPenaltyMPC_ubIdx14, beliefPenaltyMPC_ccrhsl14, beliefPenaltyMPC_slb14, beliefPenaltyMPC_lbIdx14, beliefPenaltyMPC_rd14);
beliefPenaltyMPC_LA_DENSE_MTVMADD2_6_33(params->A15, beliefPenaltyMPC_ccrhsp14, beliefPenaltyMPC_sp14, beliefPenaltyMPC_rd14);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_33(beliefPenaltyMPC_Phi14, beliefPenaltyMPC_rd14, beliefPenaltyMPC_Lbyrd14);
beliefPenaltyMPC_LA_DENSE_2MVMADD_27_87_33(beliefPenaltyMPC_V13, beliefPenaltyMPC_Lbyrd13, beliefPenaltyMPC_W14, beliefPenaltyMPC_Lbyrd14, beliefPenaltyMPC_beta14);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_27_27(beliefPenaltyMPC_Lsd14, beliefPenaltyMPC_yy13, beliefPenaltyMPC_beta14, beliefPenaltyMPC_bmy14);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_27(beliefPenaltyMPC_Ld14, beliefPenaltyMPC_bmy14, beliefPenaltyMPC_yy14);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_27(beliefPenaltyMPC_Ld14, beliefPenaltyMPC_yy14, beliefPenaltyMPC_dvcc14);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_27_27(beliefPenaltyMPC_Lsd14, beliefPenaltyMPC_dvcc14, beliefPenaltyMPC_yy13, beliefPenaltyMPC_bmy13);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_27(beliefPenaltyMPC_Ld13, beliefPenaltyMPC_bmy13, beliefPenaltyMPC_dvcc13);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_27_27(beliefPenaltyMPC_Lsd13, beliefPenaltyMPC_dvcc13, beliefPenaltyMPC_yy12, beliefPenaltyMPC_bmy12);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_27(beliefPenaltyMPC_Ld12, beliefPenaltyMPC_bmy12, beliefPenaltyMPC_dvcc12);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_27_27(beliefPenaltyMPC_Lsd12, beliefPenaltyMPC_dvcc12, beliefPenaltyMPC_yy11, beliefPenaltyMPC_bmy11);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_27(beliefPenaltyMPC_Ld11, beliefPenaltyMPC_bmy11, beliefPenaltyMPC_dvcc11);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_27_27(beliefPenaltyMPC_Lsd11, beliefPenaltyMPC_dvcc11, beliefPenaltyMPC_yy10, beliefPenaltyMPC_bmy10);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_27(beliefPenaltyMPC_Ld10, beliefPenaltyMPC_bmy10, beliefPenaltyMPC_dvcc10);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_27_27(beliefPenaltyMPC_Lsd10, beliefPenaltyMPC_dvcc10, beliefPenaltyMPC_yy09, beliefPenaltyMPC_bmy09);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_27(beliefPenaltyMPC_Ld09, beliefPenaltyMPC_bmy09, beliefPenaltyMPC_dvcc09);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_27_27(beliefPenaltyMPC_Lsd09, beliefPenaltyMPC_dvcc09, beliefPenaltyMPC_yy08, beliefPenaltyMPC_bmy08);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_27(beliefPenaltyMPC_Ld08, beliefPenaltyMPC_bmy08, beliefPenaltyMPC_dvcc08);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_27_27(beliefPenaltyMPC_Lsd08, beliefPenaltyMPC_dvcc08, beliefPenaltyMPC_yy07, beliefPenaltyMPC_bmy07);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_27(beliefPenaltyMPC_Ld07, beliefPenaltyMPC_bmy07, beliefPenaltyMPC_dvcc07);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_27_27(beliefPenaltyMPC_Lsd07, beliefPenaltyMPC_dvcc07, beliefPenaltyMPC_yy06, beliefPenaltyMPC_bmy06);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_27(beliefPenaltyMPC_Ld06, beliefPenaltyMPC_bmy06, beliefPenaltyMPC_dvcc06);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_27_27(beliefPenaltyMPC_Lsd06, beliefPenaltyMPC_dvcc06, beliefPenaltyMPC_yy05, beliefPenaltyMPC_bmy05);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_27(beliefPenaltyMPC_Ld05, beliefPenaltyMPC_bmy05, beliefPenaltyMPC_dvcc05);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_27_27(beliefPenaltyMPC_Lsd05, beliefPenaltyMPC_dvcc05, beliefPenaltyMPC_yy04, beliefPenaltyMPC_bmy04);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_27(beliefPenaltyMPC_Ld04, beliefPenaltyMPC_bmy04, beliefPenaltyMPC_dvcc04);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_27_27(beliefPenaltyMPC_Lsd04, beliefPenaltyMPC_dvcc04, beliefPenaltyMPC_yy03, beliefPenaltyMPC_bmy03);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_27(beliefPenaltyMPC_Ld03, beliefPenaltyMPC_bmy03, beliefPenaltyMPC_dvcc03);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_27_27(beliefPenaltyMPC_Lsd03, beliefPenaltyMPC_dvcc03, beliefPenaltyMPC_yy02, beliefPenaltyMPC_bmy02);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_27(beliefPenaltyMPC_Ld02, beliefPenaltyMPC_bmy02, beliefPenaltyMPC_dvcc02);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_27_27(beliefPenaltyMPC_Lsd02, beliefPenaltyMPC_dvcc02, beliefPenaltyMPC_yy01, beliefPenaltyMPC_bmy01);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_27(beliefPenaltyMPC_Ld01, beliefPenaltyMPC_bmy01, beliefPenaltyMPC_dvcc01);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_27_27(beliefPenaltyMPC_Lsd01, beliefPenaltyMPC_dvcc01, beliefPenaltyMPC_yy00, beliefPenaltyMPC_bmy00);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_27(beliefPenaltyMPC_Ld00, beliefPenaltyMPC_bmy00, beliefPenaltyMPC_dvcc00);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_27_87_27(params->C1, beliefPenaltyMPC_dvcc01, beliefPenaltyMPC_D00, beliefPenaltyMPC_dvcc00, beliefPenaltyMPC_grad_eq00);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_27_87_27(params->C2, beliefPenaltyMPC_dvcc02, beliefPenaltyMPC_D01, beliefPenaltyMPC_dvcc01, beliefPenaltyMPC_grad_eq01);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_27_87_27(params->C3, beliefPenaltyMPC_dvcc03, beliefPenaltyMPC_D01, beliefPenaltyMPC_dvcc02, beliefPenaltyMPC_grad_eq02);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_27_87_27(params->C4, beliefPenaltyMPC_dvcc04, beliefPenaltyMPC_D01, beliefPenaltyMPC_dvcc03, beliefPenaltyMPC_grad_eq03);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_27_87_27(params->C5, beliefPenaltyMPC_dvcc05, beliefPenaltyMPC_D01, beliefPenaltyMPC_dvcc04, beliefPenaltyMPC_grad_eq04);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_27_87_27(params->C6, beliefPenaltyMPC_dvcc06, beliefPenaltyMPC_D01, beliefPenaltyMPC_dvcc05, beliefPenaltyMPC_grad_eq05);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_27_87_27(params->C7, beliefPenaltyMPC_dvcc07, beliefPenaltyMPC_D01, beliefPenaltyMPC_dvcc06, beliefPenaltyMPC_grad_eq06);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_27_87_27(params->C8, beliefPenaltyMPC_dvcc08, beliefPenaltyMPC_D01, beliefPenaltyMPC_dvcc07, beliefPenaltyMPC_grad_eq07);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_27_87_27(params->C9, beliefPenaltyMPC_dvcc09, beliefPenaltyMPC_D01, beliefPenaltyMPC_dvcc08, beliefPenaltyMPC_grad_eq08);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_27_87_27(params->C10, beliefPenaltyMPC_dvcc10, beliefPenaltyMPC_D01, beliefPenaltyMPC_dvcc09, beliefPenaltyMPC_grad_eq09);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_27_87_27(params->C11, beliefPenaltyMPC_dvcc11, beliefPenaltyMPC_D01, beliefPenaltyMPC_dvcc10, beliefPenaltyMPC_grad_eq10);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_27_87_27(params->C12, beliefPenaltyMPC_dvcc12, beliefPenaltyMPC_D01, beliefPenaltyMPC_dvcc11, beliefPenaltyMPC_grad_eq11);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_27_87_27(params->C13, beliefPenaltyMPC_dvcc13, beliefPenaltyMPC_D01, beliefPenaltyMPC_dvcc12, beliefPenaltyMPC_grad_eq12);
beliefPenaltyMPC_LA_DENSE_DIAGZERO_MTVM2_27_87_27(params->C14, beliefPenaltyMPC_dvcc14, beliefPenaltyMPC_D01, beliefPenaltyMPC_dvcc13, beliefPenaltyMPC_grad_eq13);
beliefPenaltyMPC_LA_DIAGZERO_MTVM_27_33(beliefPenaltyMPC_D14, beliefPenaltyMPC_dvcc14, beliefPenaltyMPC_grad_eq14);
beliefPenaltyMPC_LA_VSUB_1251(beliefPenaltyMPC_rd, beliefPenaltyMPC_grad_eq, beliefPenaltyMPC_rd);
beliefPenaltyMPC_LA_DENSE_FORWARDBACKWARDSUB_87(beliefPenaltyMPC_Phi00, beliefPenaltyMPC_rd00, beliefPenaltyMPC_dzcc00);
beliefPenaltyMPC_LA_DENSE_FORWARDBACKWARDSUB_87(beliefPenaltyMPC_Phi01, beliefPenaltyMPC_rd01, beliefPenaltyMPC_dzcc01);
beliefPenaltyMPC_LA_DENSE_FORWARDBACKWARDSUB_87(beliefPenaltyMPC_Phi02, beliefPenaltyMPC_rd02, beliefPenaltyMPC_dzcc02);
beliefPenaltyMPC_LA_DENSE_FORWARDBACKWARDSUB_87(beliefPenaltyMPC_Phi03, beliefPenaltyMPC_rd03, beliefPenaltyMPC_dzcc03);
beliefPenaltyMPC_LA_DENSE_FORWARDBACKWARDSUB_87(beliefPenaltyMPC_Phi04, beliefPenaltyMPC_rd04, beliefPenaltyMPC_dzcc04);
beliefPenaltyMPC_LA_DENSE_FORWARDBACKWARDSUB_87(beliefPenaltyMPC_Phi05, beliefPenaltyMPC_rd05, beliefPenaltyMPC_dzcc05);
beliefPenaltyMPC_LA_DENSE_FORWARDBACKWARDSUB_87(beliefPenaltyMPC_Phi06, beliefPenaltyMPC_rd06, beliefPenaltyMPC_dzcc06);
beliefPenaltyMPC_LA_DENSE_FORWARDBACKWARDSUB_87(beliefPenaltyMPC_Phi07, beliefPenaltyMPC_rd07, beliefPenaltyMPC_dzcc07);
beliefPenaltyMPC_LA_DENSE_FORWARDBACKWARDSUB_87(beliefPenaltyMPC_Phi08, beliefPenaltyMPC_rd08, beliefPenaltyMPC_dzcc08);
beliefPenaltyMPC_LA_DENSE_FORWARDBACKWARDSUB_87(beliefPenaltyMPC_Phi09, beliefPenaltyMPC_rd09, beliefPenaltyMPC_dzcc09);
beliefPenaltyMPC_LA_DENSE_FORWARDBACKWARDSUB_87(beliefPenaltyMPC_Phi10, beliefPenaltyMPC_rd10, beliefPenaltyMPC_dzcc10);
beliefPenaltyMPC_LA_DENSE_FORWARDBACKWARDSUB_87(beliefPenaltyMPC_Phi11, beliefPenaltyMPC_rd11, beliefPenaltyMPC_dzcc11);
beliefPenaltyMPC_LA_DENSE_FORWARDBACKWARDSUB_87(beliefPenaltyMPC_Phi12, beliefPenaltyMPC_rd12, beliefPenaltyMPC_dzcc12);
beliefPenaltyMPC_LA_DENSE_FORWARDBACKWARDSUB_87(beliefPenaltyMPC_Phi13, beliefPenaltyMPC_rd13, beliefPenaltyMPC_dzcc13);
beliefPenaltyMPC_LA_DENSE_FORWARDBACKWARDSUB_33(beliefPenaltyMPC_Phi14, beliefPenaltyMPC_rd14, beliefPenaltyMPC_dzcc14);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_87(beliefPenaltyMPC_ccrhsl00, beliefPenaltyMPC_slb00, beliefPenaltyMPC_llbbyslb00, beliefPenaltyMPC_dzcc00, beliefPenaltyMPC_lbIdx00, beliefPenaltyMPC_dllbcc00);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_33(beliefPenaltyMPC_ccrhsub00, beliefPenaltyMPC_sub00, beliefPenaltyMPC_lubbysub00, beliefPenaltyMPC_dzcc00, beliefPenaltyMPC_ubIdx00, beliefPenaltyMPC_dlubcc00);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_87(beliefPenaltyMPC_ccrhsl01, beliefPenaltyMPC_slb01, beliefPenaltyMPC_llbbyslb01, beliefPenaltyMPC_dzcc01, beliefPenaltyMPC_lbIdx01, beliefPenaltyMPC_dllbcc01);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_33(beliefPenaltyMPC_ccrhsub01, beliefPenaltyMPC_sub01, beliefPenaltyMPC_lubbysub01, beliefPenaltyMPC_dzcc01, beliefPenaltyMPC_ubIdx01, beliefPenaltyMPC_dlubcc01);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_87(beliefPenaltyMPC_ccrhsl02, beliefPenaltyMPC_slb02, beliefPenaltyMPC_llbbyslb02, beliefPenaltyMPC_dzcc02, beliefPenaltyMPC_lbIdx02, beliefPenaltyMPC_dllbcc02);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_33(beliefPenaltyMPC_ccrhsub02, beliefPenaltyMPC_sub02, beliefPenaltyMPC_lubbysub02, beliefPenaltyMPC_dzcc02, beliefPenaltyMPC_ubIdx02, beliefPenaltyMPC_dlubcc02);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_87(beliefPenaltyMPC_ccrhsl03, beliefPenaltyMPC_slb03, beliefPenaltyMPC_llbbyslb03, beliefPenaltyMPC_dzcc03, beliefPenaltyMPC_lbIdx03, beliefPenaltyMPC_dllbcc03);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_33(beliefPenaltyMPC_ccrhsub03, beliefPenaltyMPC_sub03, beliefPenaltyMPC_lubbysub03, beliefPenaltyMPC_dzcc03, beliefPenaltyMPC_ubIdx03, beliefPenaltyMPC_dlubcc03);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_87(beliefPenaltyMPC_ccrhsl04, beliefPenaltyMPC_slb04, beliefPenaltyMPC_llbbyslb04, beliefPenaltyMPC_dzcc04, beliefPenaltyMPC_lbIdx04, beliefPenaltyMPC_dllbcc04);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_33(beliefPenaltyMPC_ccrhsub04, beliefPenaltyMPC_sub04, beliefPenaltyMPC_lubbysub04, beliefPenaltyMPC_dzcc04, beliefPenaltyMPC_ubIdx04, beliefPenaltyMPC_dlubcc04);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_87(beliefPenaltyMPC_ccrhsl05, beliefPenaltyMPC_slb05, beliefPenaltyMPC_llbbyslb05, beliefPenaltyMPC_dzcc05, beliefPenaltyMPC_lbIdx05, beliefPenaltyMPC_dllbcc05);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_33(beliefPenaltyMPC_ccrhsub05, beliefPenaltyMPC_sub05, beliefPenaltyMPC_lubbysub05, beliefPenaltyMPC_dzcc05, beliefPenaltyMPC_ubIdx05, beliefPenaltyMPC_dlubcc05);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_87(beliefPenaltyMPC_ccrhsl06, beliefPenaltyMPC_slb06, beliefPenaltyMPC_llbbyslb06, beliefPenaltyMPC_dzcc06, beliefPenaltyMPC_lbIdx06, beliefPenaltyMPC_dllbcc06);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_33(beliefPenaltyMPC_ccrhsub06, beliefPenaltyMPC_sub06, beliefPenaltyMPC_lubbysub06, beliefPenaltyMPC_dzcc06, beliefPenaltyMPC_ubIdx06, beliefPenaltyMPC_dlubcc06);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_87(beliefPenaltyMPC_ccrhsl07, beliefPenaltyMPC_slb07, beliefPenaltyMPC_llbbyslb07, beliefPenaltyMPC_dzcc07, beliefPenaltyMPC_lbIdx07, beliefPenaltyMPC_dllbcc07);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_33(beliefPenaltyMPC_ccrhsub07, beliefPenaltyMPC_sub07, beliefPenaltyMPC_lubbysub07, beliefPenaltyMPC_dzcc07, beliefPenaltyMPC_ubIdx07, beliefPenaltyMPC_dlubcc07);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_87(beliefPenaltyMPC_ccrhsl08, beliefPenaltyMPC_slb08, beliefPenaltyMPC_llbbyslb08, beliefPenaltyMPC_dzcc08, beliefPenaltyMPC_lbIdx08, beliefPenaltyMPC_dllbcc08);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_33(beliefPenaltyMPC_ccrhsub08, beliefPenaltyMPC_sub08, beliefPenaltyMPC_lubbysub08, beliefPenaltyMPC_dzcc08, beliefPenaltyMPC_ubIdx08, beliefPenaltyMPC_dlubcc08);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_87(beliefPenaltyMPC_ccrhsl09, beliefPenaltyMPC_slb09, beliefPenaltyMPC_llbbyslb09, beliefPenaltyMPC_dzcc09, beliefPenaltyMPC_lbIdx09, beliefPenaltyMPC_dllbcc09);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_33(beliefPenaltyMPC_ccrhsub09, beliefPenaltyMPC_sub09, beliefPenaltyMPC_lubbysub09, beliefPenaltyMPC_dzcc09, beliefPenaltyMPC_ubIdx09, beliefPenaltyMPC_dlubcc09);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_87(beliefPenaltyMPC_ccrhsl10, beliefPenaltyMPC_slb10, beliefPenaltyMPC_llbbyslb10, beliefPenaltyMPC_dzcc10, beliefPenaltyMPC_lbIdx10, beliefPenaltyMPC_dllbcc10);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_33(beliefPenaltyMPC_ccrhsub10, beliefPenaltyMPC_sub10, beliefPenaltyMPC_lubbysub10, beliefPenaltyMPC_dzcc10, beliefPenaltyMPC_ubIdx10, beliefPenaltyMPC_dlubcc10);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_87(beliefPenaltyMPC_ccrhsl11, beliefPenaltyMPC_slb11, beliefPenaltyMPC_llbbyslb11, beliefPenaltyMPC_dzcc11, beliefPenaltyMPC_lbIdx11, beliefPenaltyMPC_dllbcc11);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_33(beliefPenaltyMPC_ccrhsub11, beliefPenaltyMPC_sub11, beliefPenaltyMPC_lubbysub11, beliefPenaltyMPC_dzcc11, beliefPenaltyMPC_ubIdx11, beliefPenaltyMPC_dlubcc11);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_87(beliefPenaltyMPC_ccrhsl12, beliefPenaltyMPC_slb12, beliefPenaltyMPC_llbbyslb12, beliefPenaltyMPC_dzcc12, beliefPenaltyMPC_lbIdx12, beliefPenaltyMPC_dllbcc12);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_33(beliefPenaltyMPC_ccrhsub12, beliefPenaltyMPC_sub12, beliefPenaltyMPC_lubbysub12, beliefPenaltyMPC_dzcc12, beliefPenaltyMPC_ubIdx12, beliefPenaltyMPC_dlubcc12);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_87(beliefPenaltyMPC_ccrhsl13, beliefPenaltyMPC_slb13, beliefPenaltyMPC_llbbyslb13, beliefPenaltyMPC_dzcc13, beliefPenaltyMPC_lbIdx13, beliefPenaltyMPC_dllbcc13);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_33(beliefPenaltyMPC_ccrhsub13, beliefPenaltyMPC_sub13, beliefPenaltyMPC_lubbysub13, beliefPenaltyMPC_dzcc13, beliefPenaltyMPC_ubIdx13, beliefPenaltyMPC_dlubcc13);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_33(beliefPenaltyMPC_ccrhsl14, beliefPenaltyMPC_slb14, beliefPenaltyMPC_llbbyslb14, beliefPenaltyMPC_dzcc14, beliefPenaltyMPC_lbIdx14, beliefPenaltyMPC_dllbcc14);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_27(beliefPenaltyMPC_ccrhsub14, beliefPenaltyMPC_sub14, beliefPenaltyMPC_lubbysub14, beliefPenaltyMPC_dzcc14, beliefPenaltyMPC_ubIdx14, beliefPenaltyMPC_dlubcc14);
beliefPenaltyMPC_LA_DENSE_MVMSUB5_6_33(params->A15, beliefPenaltyMPC_dzcc14, beliefPenaltyMPC_ccrhsp14, beliefPenaltyMPC_sp14, beliefPenaltyMPC_lp14, beliefPenaltyMPC_dlp_cc14);
beliefPenaltyMPC_LA_VSUB7_1746(beliefPenaltyMPC_l, beliefPenaltyMPC_ccrhs, beliefPenaltyMPC_s, beliefPenaltyMPC_dl_cc, beliefPenaltyMPC_ds_cc);
beliefPenaltyMPC_LA_VADD_1251(beliefPenaltyMPC_dz_cc, beliefPenaltyMPC_dz_aff);
beliefPenaltyMPC_LA_VADD_405(beliefPenaltyMPC_dv_cc, beliefPenaltyMPC_dv_aff);
beliefPenaltyMPC_LA_VADD_1746(beliefPenaltyMPC_dl_cc, beliefPenaltyMPC_dl_aff);
beliefPenaltyMPC_LA_VADD_1746(beliefPenaltyMPC_ds_cc, beliefPenaltyMPC_ds_aff);
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
