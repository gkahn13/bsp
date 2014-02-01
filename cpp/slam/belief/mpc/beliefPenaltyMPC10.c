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
#ifndef USEMEXPRINTS
#include <stdio.h>
#define PRINTTEXT printf
#else
#include "mex.h"
#define PRINTTEXT mexPrintf
#endif

/* TIMING LIBRARY ------------------------------------------------- */

/* ARE WE ON WINDOWS? */
#if (defined WIN32 || defined _WIN64 || defined _WIN32)

/* Use Windows QueryPerformanceCounter for timing */

#include <windows.h>

typedef struct beliefPenaltyMPC_timer{
	LARGE_INTEGER tic;
	LARGE_INTEGER toc;
	LARGE_INTEGER freq;
} beliefPenaltyMPC_timer;


void beliefPenaltyMPC_tic(beliefPenaltyMPC_timer* t)
{
	QueryPerformanceFrequency(&t->freq);
	QueryPerformanceCounter(&t->tic);
}



beliefPenaltyMPC_FLOAT beliefPenaltyMPC_toc(beliefPenaltyMPC_timer* t)
{
	QueryPerformanceCounter(&t->toc);
	return ((t->toc.QuadPart - t->tic.QuadPart) / (beliefPenaltyMPC_FLOAT)t->freq.QuadPart);
}


/* WE ARE ON THE MAC */
#elif (defined __APPLE__)
#include <mach/mach_time.h>


/* Use MAC OSX  mach_time for timing */
typedef struct beliefPenaltyMPC_timer{
	uint64_t tic;
	uint64_t toc;
	mach_timebase_info_data_t tinfo;

} beliefPenaltyMPC_timer;


void beliefPenaltyMPC_tic(beliefPenaltyMPC_timer* t)
{
    /* read current clock cycles */
    t->tic = mach_absolute_time();
}



beliefPenaltyMPC_FLOAT beliefPenaltyMPC_toc(beliefPenaltyMPC_timer* t)
{
    uint64_t duration; /* elapsed time in clock cycles*/
    t->toc = mach_absolute_time();
	duration = t->toc - t->tic;

    /*conversion from clock cycles to nanoseconds*/
    mach_timebase_info(&(t->tinfo));
    duration *= t->tinfo.numer;
    duration /= t->tinfo.denom;

    return (beliefPenaltyMPC_FLOAT)duration / 1000000000;
}

/* WE ARE ON SOME TEXAS INSTRUMENTS PLATFORM */
#elif (defined __TI_COMPILER_VERSION__)

/* TimeStamps */
#include <c6x.h> /* make use of TSCL, TSCH */


typedef struct beliefPenaltyMPC_timer{
	unsigned long long tic;
	unsigned long long toc;
} beliefPenaltyMPC_timer;


void beliefPenaltyMPC_tic(beliefPenaltyMPC_timer* t)
{
	TSCL = 0;	/* Initiate CPU timer by writing any val to TSCL */
	t->tic = _itoll( TSCH, TSCL );
}



beliefPenaltyMPC_FLOAT beliefPenaltyMPC_toc(beliefPenaltyMPC_timer* t)
{
	t->toc = _itoll( TSCH, TSCL );
	unsigned long long t0;
	unsigned long long overhead;
	t0 = _itoll( TSCH, TSCL );
	overhead = _itoll( TSCH, TSCL )  - t0;

	return (beliefPenaltyMPC_FLOAT)(t->toc - t->tic - overhead) / 1000000000;
}



/* WE ARE ON SOME OTHER UNIX/LINUX SYSTEM */
#else

/* Use POSIX clocl_gettime() for timing on non-Windows machines */
#include <time.h>
typedef struct beliefPenaltyMPC_timer{
	struct timespec tic;
	struct timespec toc;
} beliefPenaltyMPC_timer;


/* read current time */
void beliefPenaltyMPC_tic(beliefPenaltyMPC_timer* t)
{
	clock_gettime(CLOCK_MONOTONIC, &t->tic);
}



/* return time passed since last call to tic on this timer */
double beliefPenaltyMPC_toc(beliefPenaltyMPC_timer* t)
{
	struct timespec temp;
	clock_gettime(CLOCK_MONOTONIC, &t->toc);	

	if ((t->toc.tv_nsec - t->tic.tv_nsec)<0) {
		temp.tv_sec = t->toc.tv_sec - t->tic.tv_sec-1;
		temp.tv_nsec = 1000000000+t->toc.tv_nsec - t->tic.tv_nsec;
	} else {
		temp.tv_sec = t->toc.tv_sec - t->tic.tv_sec;
		temp.tv_nsec = t->toc.tv_nsec - t->tic.tv_nsec;
	}

	return (beliefPenaltyMPC_FLOAT)temp.tv_sec + (beliefPenaltyMPC_FLOAT)temp.tv_nsec / 1000000000;
}


#endif

/* LINEAR ALGEBRA LIBRARY ---------------------------------------------- */
/*
 * Initializes a vector of length 3029 with a value.
 */
void beliefPenaltyMPC_LA_INITIALIZEVECTOR_3029(beliefPenaltyMPC_FLOAT* vec, beliefPenaltyMPC_FLOAT value)
{
	int i;
	for( i=0; i<3029; i++ )
	{
		vec[i] = value;
	}
}


/*
 * Initializes a vector of length 1040 with a value.
 */
void beliefPenaltyMPC_LA_INITIALIZEVECTOR_1040(beliefPenaltyMPC_FLOAT* vec, beliefPenaltyMPC_FLOAT value)
{
	int i;
	for( i=0; i<1040; i++ )
	{
		vec[i] = value;
	}
}


/*
 * Initializes a vector of length 4186 with a value.
 */
void beliefPenaltyMPC_LA_INITIALIZEVECTOR_4186(beliefPenaltyMPC_FLOAT* vec, beliefPenaltyMPC_FLOAT value)
{
	int i;
	for( i=0; i<4186; i++ )
	{
		vec[i] = value;
	}
}


/* 
 * Calculates a dot product and adds it to a variable: z += x'*y; 
 * This function is for vectors of length 4186.
 */
void beliefPenaltyMPC_LA_DOTACC_4186(beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *y, beliefPenaltyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<4186; i++ ){
		*z += x[i]*y[i];
	}
}


/*
 * Calculates the gradient and the value for a quadratic function 0.5*z'*H*z + f'*z
 *
 * INPUTS:     H  - Symmetric Hessian, diag matrix of size [325 x 325]
 *             f  - column vector of size 325
 *             z  - column vector of size 325
 *
 * OUTPUTS: grad  - gradient at z (= H*z + f), column vector of size 325
 *          value <-- value + 0.5*z'*H*z + f'*z (value will be modified)
 */
void beliefPenaltyMPC_LA_DIAG_QUADFCN_325(beliefPenaltyMPC_FLOAT* H, beliefPenaltyMPC_FLOAT* f, beliefPenaltyMPC_FLOAT* z, beliefPenaltyMPC_FLOAT* grad, beliefPenaltyMPC_FLOAT* value)
{
	int i;
	beliefPenaltyMPC_FLOAT hz;	
	for( i=0; i<325; i++){
		hz = H[i]*z[i];
		grad[i] = hz + f[i];
		*value += 0.5*hz*z[i] + f[i]*z[i];
	}
}


/*
 * Calculates the gradient and the value for a quadratic function 0.5*z'*H*z + f'*z
 *
 * INPUTS:     H  - Symmetric Hessian, diag matrix of size [104 x 104]
 *             f  - column vector of size 104
 *             z  - column vector of size 104
 *
 * OUTPUTS: grad  - gradient at z (= H*z + f), column vector of size 104
 *          value <-- value + 0.5*z'*H*z + f'*z (value will be modified)
 */
void beliefPenaltyMPC_LA_DIAG_QUADFCN_104(beliefPenaltyMPC_FLOAT* H, beliefPenaltyMPC_FLOAT* f, beliefPenaltyMPC_FLOAT* z, beliefPenaltyMPC_FLOAT* grad, beliefPenaltyMPC_FLOAT* value)
{
	int i;
	beliefPenaltyMPC_FLOAT hz;	
	for( i=0; i<104; i++){
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
void beliefPenaltyMPC_LA_DENSE_MVMSUB3_208_325_325(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *b, beliefPenaltyMPC_FLOAT *l, beliefPenaltyMPC_FLOAT *r, beliefPenaltyMPC_FLOAT *z, beliefPenaltyMPC_FLOAT *y)
{
	int i;
	int j;
	int k = 0;
	int m = 0;
	int n;
	beliefPenaltyMPC_FLOAT AxBu[208];
	beliefPenaltyMPC_FLOAT norm = *y;
	beliefPenaltyMPC_FLOAT lr = 0;

	/* do A*x + B*u first */
	for( i=0; i<208; i++ ){
		AxBu[i] = A[k++]*x[0] + B[m++]*u[0];
	}	
	for( j=1; j<325; j++ ){		
		for( i=0; i<208; i++ ){
			AxBu[i] += A[k++]*x[j];
		}
	}
	
	for( n=1; n<325; n++ ){
		for( i=0; i<208; i++ ){
			AxBu[i] += B[m++]*u[n];
		}		
	}

	for( i=0; i<208; i++ ){
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
void beliefPenaltyMPC_LA_DENSE_MVMSUB3_104_325_325(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *b, beliefPenaltyMPC_FLOAT *l, beliefPenaltyMPC_FLOAT *r, beliefPenaltyMPC_FLOAT *z, beliefPenaltyMPC_FLOAT *y)
{
	int i;
	int j;
	int k = 0;
	int m = 0;
	int n;
	beliefPenaltyMPC_FLOAT AxBu[104];
	beliefPenaltyMPC_FLOAT norm = *y;
	beliefPenaltyMPC_FLOAT lr = 0;

	/* do A*x + B*u first */
	for( i=0; i<104; i++ ){
		AxBu[i] = A[k++]*x[0] + B[m++]*u[0];
	}	
	for( j=1; j<325; j++ ){		
		for( i=0; i<104; i++ ){
			AxBu[i] += A[k++]*x[j];
		}
	}
	
	for( n=1; n<325; n++ ){
		for( i=0; i<104; i++ ){
			AxBu[i] += B[m++]*u[n];
		}		
	}

	for( i=0; i<104; i++ ){
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
void beliefPenaltyMPC_LA_DENSE_MVMSUB3_104_325_104(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *b, beliefPenaltyMPC_FLOAT *l, beliefPenaltyMPC_FLOAT *r, beliefPenaltyMPC_FLOAT *z, beliefPenaltyMPC_FLOAT *y)
{
	int i;
	int j;
	int k = 0;
	int m = 0;
	int n;
	beliefPenaltyMPC_FLOAT AxBu[104];
	beliefPenaltyMPC_FLOAT norm = *y;
	beliefPenaltyMPC_FLOAT lr = 0;

	/* do A*x + B*u first */
	for( i=0; i<104; i++ ){
		AxBu[i] = A[k++]*x[0] + B[m++]*u[0];
	}	
	for( j=1; j<325; j++ ){		
		for( i=0; i<104; i++ ){
			AxBu[i] += A[k++]*x[j];
		}
	}
	
	for( n=1; n<104; n++ ){
		for( i=0; i<104; i++ ){
			AxBu[i] += B[m++]*u[n];
		}		
	}

	for( i=0; i<104; i++ ){
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
 * Matrix vector multiplication y = M'*x where M is of size [208 x 325]
 * and stored in column major format. Note the transpose of M!
 */
void beliefPenaltyMPC_LA_DENSE_MTVM_208_325(beliefPenaltyMPC_FLOAT *M, beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *y)
{
	int i;
	int j;
	int k = 0; 
	for( i=0; i<325; i++ ){
		y[i] = 0;
		for( j=0; j<208; j++ ){
			y[i] += M[k++]*x[j];
		}
	}
}


/*
 * Matrix vector multiplication z = A'*x + B'*y 
 * where A is of size [104 x 325]
 * and B is of size [208 x 325]
 * and stored in column major format. Note the transposes of A and B!
 */
void beliefPenaltyMPC_LA_DENSE_MTVM2_104_325_208(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *y, beliefPenaltyMPC_FLOAT *z)
{
	int i;
	int j;
	int k = 0;
	int n;
	int m = 0;
	for( i=0; i<325; i++ ){
		z[i] = 0;
		for( j=0; j<104; j++ ){
			z[i] += A[k++]*x[j];
		}
		for( n=0; n<208; n++ ){
			z[i] += B[m++]*y[n];
		}
	}
}


/*
 * Matrix vector multiplication z = A'*x + B'*y 
 * where A is of size [104 x 325]
 * and B is of size [104 x 325]
 * and stored in column major format. Note the transposes of A and B!
 */
void beliefPenaltyMPC_LA_DENSE_MTVM2_104_325_104(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *y, beliefPenaltyMPC_FLOAT *z)
{
	int i;
	int j;
	int k = 0;
	int n;
	int m = 0;
	for( i=0; i<325; i++ ){
		z[i] = 0;
		for( j=0; j<104; j++ ){
			z[i] += A[k++]*x[j];
		}
		for( n=0; n<104; n++ ){
			z[i] += B[m++]*y[n];
		}
	}
}


/*
 * Matrix vector multiplication y = M'*x where M is of size [104 x 104]
 * and stored in column major format. Note the transpose of M!
 */
void beliefPenaltyMPC_LA_DENSE_MTVM_104_104(beliefPenaltyMPC_FLOAT *M, beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *y)
{
	int i;
	int j;
	int k = 0; 
	for( i=0; i<104; i++ ){
		y[i] = 0;
		for( j=0; j<104; j++ ){
			y[i] += M[k++]*x[j];
		}
	}
}


/*
 * Vector subtraction and addition.
 *	 Input: five vectors t, tidx, u, v, w and two scalars z and r
 *	 Output: y = t(tidx) - u + w
 *           z = z - v'*x;
 *           r = max([norm(y,inf), z]);
 * for vectors of length 325. Output z is of course scalar.
 */
void beliefPenaltyMPC_LA_VSUBADD3_325(beliefPenaltyMPC_FLOAT* t, beliefPenaltyMPC_FLOAT* u, int* uidx, beliefPenaltyMPC_FLOAT* v, beliefPenaltyMPC_FLOAT* w, beliefPenaltyMPC_FLOAT* y, beliefPenaltyMPC_FLOAT* z, beliefPenaltyMPC_FLOAT* r)
{
	int i;
	beliefPenaltyMPC_FLOAT norm = *r;
	beliefPenaltyMPC_FLOAT vx = 0;
	beliefPenaltyMPC_FLOAT x;
	for( i=0; i<325; i++){
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
 * for vectors of length 117. Output z is of course scalar.
 */
void beliefPenaltyMPC_LA_VSUBADD2_117(beliefPenaltyMPC_FLOAT* t, int* tidx, beliefPenaltyMPC_FLOAT* u, beliefPenaltyMPC_FLOAT* v, beliefPenaltyMPC_FLOAT* w, beliefPenaltyMPC_FLOAT* y, beliefPenaltyMPC_FLOAT* z, beliefPenaltyMPC_FLOAT* r)
{
	int i;
	beliefPenaltyMPC_FLOAT norm = *r;
	beliefPenaltyMPC_FLOAT vx = 0;
	beliefPenaltyMPC_FLOAT x;
	for( i=0; i<117; i++){
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
 * for vectors of length 104. Output z is of course scalar.
 */
void beliefPenaltyMPC_LA_VSUBADD3_104(beliefPenaltyMPC_FLOAT* t, beliefPenaltyMPC_FLOAT* u, int* uidx, beliefPenaltyMPC_FLOAT* v, beliefPenaltyMPC_FLOAT* w, beliefPenaltyMPC_FLOAT* y, beliefPenaltyMPC_FLOAT* z, beliefPenaltyMPC_FLOAT* r)
{
	int i;
	beliefPenaltyMPC_FLOAT norm = *r;
	beliefPenaltyMPC_FLOAT vx = 0;
	beliefPenaltyMPC_FLOAT x;
	for( i=0; i<104; i++){
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
 * for vectors of length 104. Output z is of course scalar.
 */
void beliefPenaltyMPC_LA_VSUBADD2_104(beliefPenaltyMPC_FLOAT* t, int* tidx, beliefPenaltyMPC_FLOAT* u, beliefPenaltyMPC_FLOAT* v, beliefPenaltyMPC_FLOAT* w, beliefPenaltyMPC_FLOAT* y, beliefPenaltyMPC_FLOAT* z, beliefPenaltyMPC_FLOAT* r)
{
	int i;
	beliefPenaltyMPC_FLOAT norm = *r;
	beliefPenaltyMPC_FLOAT vx = 0;
	beliefPenaltyMPC_FLOAT x;
	for( i=0; i<104; i++){
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
 * Special function for box constraints of length 325
 * Returns also L/S, a value that is often used elsewhere.
 */
void beliefPenaltyMPC_LA_INEQ_B_GRAD_325_325_117(beliefPenaltyMPC_FLOAT *lu, beliefPenaltyMPC_FLOAT *su, beliefPenaltyMPC_FLOAT *ru, beliefPenaltyMPC_FLOAT *ll, beliefPenaltyMPC_FLOAT *sl, beliefPenaltyMPC_FLOAT *rl, int* lbIdx, int* ubIdx, beliefPenaltyMPC_FLOAT *grad, beliefPenaltyMPC_FLOAT *lubysu, beliefPenaltyMPC_FLOAT *llbysl)
{
	int i;
	for( i=0; i<325; i++ ){
		grad[i] = 0;
	}
	for( i=0; i<325; i++ ){		
		llbysl[i] = ll[i] / sl[i];
		grad[lbIdx[i]] -= llbysl[i]*rl[i];
	}
	for( i=0; i<117; i++ ){
		lubysu[i] = lu[i] / su[i];
		grad[ubIdx[i]] += lubysu[i]*ru[i];
	}
}


/*
 * Computes inequality constraints gradient-
 * Special function for box constraints of length 104
 * Returns also L/S, a value that is often used elsewhere.
 */
void beliefPenaltyMPC_LA_INEQ_B_GRAD_104_104_104(beliefPenaltyMPC_FLOAT *lu, beliefPenaltyMPC_FLOAT *su, beliefPenaltyMPC_FLOAT *ru, beliefPenaltyMPC_FLOAT *ll, beliefPenaltyMPC_FLOAT *sl, beliefPenaltyMPC_FLOAT *rl, int* lbIdx, int* ubIdx, beliefPenaltyMPC_FLOAT *grad, beliefPenaltyMPC_FLOAT *lubysu, beliefPenaltyMPC_FLOAT *llbysl)
{
	int i;
	for( i=0; i<104; i++ ){
		grad[i] = 0;
	}
	for( i=0; i<104; i++ ){		
		llbysl[i] = ll[i] / sl[i];
		grad[lbIdx[i]] -= llbysl[i]*rl[i];
	}
	for( i=0; i<104; i++ ){
		lubysu[i] = lu[i] / su[i];
		grad[ubIdx[i]] += lubysu[i]*ru[i];
	}
}


/*
 * Addition of three vectors  z = u + w + v
 * of length 3029.
 */
void beliefPenaltyMPC_LA_VVADD3_3029(beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *v, beliefPenaltyMPC_FLOAT *w, beliefPenaltyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<3029; i++ ){
		z[i] = u[i] + v[i] + w[i];
	}
}


/*
 * Special function to compute the diagonal cholesky factorization of the 
 * positive definite augmented Hessian for block size 325.
 *
 * Inputs: - H = diagonal cost Hessian in diagonal storage format
 *         - llbysl = L / S of lower bounds
 *         - lubysu = L / S of upper bounds
 *
 * Output: Phi = sqrt(H + diag(llbysl) + diag(lubysu))
 * where Phi is stored in diagonal storage format
 */
void beliefPenaltyMPC_LA_DIAG_CHOL_LBUB_325_325_117(beliefPenaltyMPC_FLOAT *H, beliefPenaltyMPC_FLOAT *llbysl, int* lbIdx, beliefPenaltyMPC_FLOAT *lubysu, int* ubIdx, beliefPenaltyMPC_FLOAT *Phi)


{
	int i;
	
	/* copy  H into PHI */
	for( i=0; i<325; i++ ){
		Phi[i] = H[i];
	}

	/* add llbysl onto Phi where necessary */
	for( i=0; i<325; i++ ){
		Phi[lbIdx[i]] += llbysl[i];
	}

	/* add lubysu onto Phi where necessary */
	for( i=0; i<117; i++){
		Phi[ubIdx[i]] +=  lubysu[i];
	}
	
	/* compute cholesky */
	for(i=0; i<325; i++)
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
 * where A is to be computed and is of size [208 x 325],
 * B is given and of size [208 x 325], L is a diagonal
 * matrix of size 208 stored in diagonal matrix 
 * storage format. Note the transpose of L has no impact!
 *
 * Result: A in column major storage format.
 *
 */
void beliefPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_208_325(beliefPenaltyMPC_FLOAT *L, beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *A)
{
    int i,j;
	 int k = 0;

	for( j=0; j<325; j++){
		for( i=0; i<208; i++){
			A[k] = B[k]/L[j];
			k++;
		}
	}

}


/**
 * Forward substitution to solve L*y = b where L is a
 * diagonal matrix in vector storage format.
 * 
 * The dimensions involved are 325.
 */
void beliefPenaltyMPC_LA_DIAG_FORWARDSUB_325(beliefPenaltyMPC_FLOAT *L, beliefPenaltyMPC_FLOAT *b, beliefPenaltyMPC_FLOAT *y)
{
    int i;

    for( i=0; i<325; i++ ){
		y[i] = b[i]/L[i];
    }
}


/**
 * Forward substitution for the matrix equation A*L' = B
 * where A is to be computed and is of size [104 x 325],
 * B is given and of size [104 x 325], L is a diagonal
 * matrix of size 104 stored in diagonal matrix 
 * storage format. Note the transpose of L has no impact!
 *
 * Result: A in column major storage format.
 *
 */
void beliefPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_104_325(beliefPenaltyMPC_FLOAT *L, beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *A)
{
    int i,j;
	 int k = 0;

	for( j=0; j<325; j++){
		for( i=0; i<104; i++){
			A[k] = B[k]/L[j];
			k++;
		}
	}

}


/**
 * Compute C = A*B' where 
 *
 *	size(A) = [208 x 325]
 *  size(B) = [104 x 325]
 * 
 * and all matrices are stored in column major format.
 *
 * THIS ONE HAS THE WORST ACCES PATTERN POSSIBLE.  
 * 
 */
void beliefPenaltyMPC_LA_DENSE_MMTM_208_325_104(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *C)
{
    int i, j, k;
    beliefPenaltyMPC_FLOAT temp;
    
    for( i=0; i<208; i++ ){        
        for( j=0; j<104; j++ ){
            temp = 0; 
            for( k=0; k<325; k++ ){
                temp += A[k*208+i]*B[k*104+j];
            }						
            C[j*208+i] = temp;
        }
    }
}


/**
 * Compute C = A*B' where 
 *
 *	size(A) = [104 x 325]
 *  size(B) = [104 x 325]
 * 
 * and all matrices are stored in column major format.
 *
 * THIS ONE HAS THE WORST ACCES PATTERN POSSIBLE.  
 * 
 */
void beliefPenaltyMPC_LA_DENSE_MMTM_104_325_104(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *C)
{
    int i, j, k;
    beliefPenaltyMPC_FLOAT temp;
    
    for( i=0; i<104; i++ ){        
        for( j=0; j<104; j++ ){
            temp = 0; 
            for( k=0; k<325; k++ ){
                temp += A[k*104+i]*B[k*104+j];
            }						
            C[j*104+i] = temp;
        }
    }
}


/*
 * Special function to compute the diagonal cholesky factorization of the 
 * positive definite augmented Hessian for block size 104.
 *
 * Inputs: - H = diagonal cost Hessian in diagonal storage format
 *         - llbysl = L / S of lower bounds
 *         - lubysu = L / S of upper bounds
 *
 * Output: Phi = sqrt(H + diag(llbysl) + diag(lubysu))
 * where Phi is stored in diagonal storage format
 */
void beliefPenaltyMPC_LA_DIAG_CHOL_ONELOOP_LBUB_104_104_104(beliefPenaltyMPC_FLOAT *H, beliefPenaltyMPC_FLOAT *llbysl, int* lbIdx, beliefPenaltyMPC_FLOAT *lubysu, int* ubIdx, beliefPenaltyMPC_FLOAT *Phi)


{
	int i;
	
	/* compute cholesky */
	for( i=0; i<104; i++ ){
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
 * where A is to be computed and is of size [104 x 104],
 * B is given and of size [104 x 104], L is a diagonal
 * matrix of size 104 stored in diagonal matrix 
 * storage format. Note the transpose of L has no impact!
 *
 * Result: A in column major storage format.
 *
 */
void beliefPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_104_104(beliefPenaltyMPC_FLOAT *L, beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *A)
{
    int i,j;
	 int k = 0;

	for( j=0; j<104; j++){
		for( i=0; i<104; i++){
			A[k] = B[k]/L[j];
			k++;
		}
	}

}


/**
 * Forward substitution to solve L*y = b where L is a
 * diagonal matrix in vector storage format.
 * 
 * The dimensions involved are 104.
 */
void beliefPenaltyMPC_LA_DIAG_FORWARDSUB_104(beliefPenaltyMPC_FLOAT *L, beliefPenaltyMPC_FLOAT *b, beliefPenaltyMPC_FLOAT *y)
{
    int i;

    for( i=0; i<104; i++ ){
		y[i] = b[i]/L[i];
    }
}


/**
 * Compute L = A*A' + B*B', where L is lower triangular of size NXp1
 * and A is a dense matrix of size [208 x 325] in column
 * storage format, and B is of size [208 x 325] also in column
 * storage format.
 * 
 * THIS ONE HAS THE WORST ACCES PATTERN POSSIBLE. 
 * POSSIBKE FIX: PUT A AND B INTO ROW MAJOR FORMAT FIRST.
 * 
 */
void beliefPenaltyMPC_LA_DENSE_MMT2_208_325_325(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *L)
{
    int i, j, k, ii, di;
    beliefPenaltyMPC_FLOAT ltemp;
    
    ii = 0; di = 0;
    for( i=0; i<208; i++ ){        
        for( j=0; j<=i; j++ ){
            ltemp = 0; 
            for( k=0; k<325; k++ ){
                ltemp += A[k*208+i]*A[k*208+j];
            }			
			for( k=0; k<325; k++ ){
                ltemp += B[k*208+i]*B[k*208+j];
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
void beliefPenaltyMPC_LA_DENSE_MVMSUB2_208_325_325(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *b, beliefPenaltyMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;
	int m = 0;
	int n;

	for( i=0; i<208; i++ ){
		r[i] = b[i] - A[k++]*x[0] - B[m++]*u[0];
	}	
	for( j=1; j<325; j++ ){		
		for( i=0; i<208; i++ ){
			r[i] -= A[k++]*x[j];
		}
	}
	
	for( n=1; n<325; n++ ){
		for( i=0; i<208; i++ ){
			r[i] -= B[m++]*u[n];
		}		
	}
}


/**
 * Compute L = A*A' + B*B', where L is lower triangular of size NXp1
 * and A is a dense matrix of size [104 x 325] in column
 * storage format, and B is of size [104 x 325] also in column
 * storage format.
 * 
 * THIS ONE HAS THE WORST ACCES PATTERN POSSIBLE. 
 * POSSIBKE FIX: PUT A AND B INTO ROW MAJOR FORMAT FIRST.
 * 
 */
void beliefPenaltyMPC_LA_DENSE_MMT2_104_325_325(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *L)
{
    int i, j, k, ii, di;
    beliefPenaltyMPC_FLOAT ltemp;
    
    ii = 0; di = 0;
    for( i=0; i<104; i++ ){        
        for( j=0; j<=i; j++ ){
            ltemp = 0; 
            for( k=0; k<325; k++ ){
                ltemp += A[k*104+i]*A[k*104+j];
            }			
			for( k=0; k<325; k++ ){
                ltemp += B[k*104+i]*B[k*104+j];
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
void beliefPenaltyMPC_LA_DENSE_MVMSUB2_104_325_325(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *b, beliefPenaltyMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;
	int m = 0;
	int n;

	for( i=0; i<104; i++ ){
		r[i] = b[i] - A[k++]*x[0] - B[m++]*u[0];
	}	
	for( j=1; j<325; j++ ){		
		for( i=0; i<104; i++ ){
			r[i] -= A[k++]*x[j];
		}
	}
	
	for( n=1; n<325; n++ ){
		for( i=0; i<104; i++ ){
			r[i] -= B[m++]*u[n];
		}		
	}
}


/**
 * Compute L = A*A' + B*B', where L is lower triangular of size NXp1
 * and A is a dense matrix of size [104 x 325] in column
 * storage format, and B is of size [104 x 104] also in column
 * storage format.
 * 
 * THIS ONE HAS THE WORST ACCES PATTERN POSSIBLE. 
 * POSSIBKE FIX: PUT A AND B INTO ROW MAJOR FORMAT FIRST.
 * 
 */
void beliefPenaltyMPC_LA_DENSE_MMT2_104_325_104(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *L)
{
    int i, j, k, ii, di;
    beliefPenaltyMPC_FLOAT ltemp;
    
    ii = 0; di = 0;
    for( i=0; i<104; i++ ){        
        for( j=0; j<=i; j++ ){
            ltemp = 0; 
            for( k=0; k<325; k++ ){
                ltemp += A[k*104+i]*A[k*104+j];
            }			
			for( k=0; k<104; k++ ){
                ltemp += B[k*104+i]*B[k*104+j];
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
void beliefPenaltyMPC_LA_DENSE_MVMSUB2_104_325_104(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *b, beliefPenaltyMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;
	int m = 0;
	int n;

	for( i=0; i<104; i++ ){
		r[i] = b[i] - A[k++]*x[0] - B[m++]*u[0];
	}	
	for( j=1; j<325; j++ ){		
		for( i=0; i<104; i++ ){
			r[i] -= A[k++]*x[j];
		}
	}
	
	for( n=1; n<104; n++ ){
		for( i=0; i<104; i++ ){
			r[i] -= B[m++]*u[n];
		}		
	}
}


/**
 * Cholesky factorization as above, but working on a matrix in 
 * lower triangular storage format of size 208 and outputting
 * the Cholesky factor to matrix L in lower triangular format.
 */
void beliefPenaltyMPC_LA_DENSE_CHOL_208(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *L)
{
    int i, j, k, di, dj;
	 int ii, jj;

    beliefPenaltyMPC_FLOAT l;
    beliefPenaltyMPC_FLOAT Mii;

	/* copy A to L first and then operate on L */
	/* COULD BE OPTIMIZED */
	ii=0; di=0;
	for( i=0; i<208; i++ ){
		for( j=0; j<=i; j++ ){
			L[ii+j] = A[ii+j];
		}
		ii += ++di;
	}    
	
	/* factor L */
	ii=0; di=0;
    for( i=0; i<208; i++ ){
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
        for( j=i+1; j<208; j++ ){
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
 * The dimensions involved are 208.
 */
void beliefPenaltyMPC_LA_DENSE_FORWARDSUB_208(beliefPenaltyMPC_FLOAT *L, beliefPenaltyMPC_FLOAT *b, beliefPenaltyMPC_FLOAT *y)
{
    int i,j,ii,di;
    beliefPenaltyMPC_FLOAT yel;
            
    ii = 0; di = 0;
    for( i=0; i<208; i++ ){
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
 * where A is to be computed and is of size [104 x 208],
 * B is given and of size [104 x 208], L is a lower tri-
 * angular matrix of size 208 stored in lower triangular 
 * storage format. Note the transpose of L AND B!
 *
 * Result: A in column major storage format.
 *
 */
void beliefPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_104_208(beliefPenaltyMPC_FLOAT *L, beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *A)
{
    int i,j,k,ii,di;
    beliefPenaltyMPC_FLOAT a;
    
    ii=0; di=0;
    for( j=0; j<208; j++ ){        
        for( i=0; i<104; i++ ){
            a = B[i*208+j];
            for( k=0; k<j; k++ ){
                a -= A[k*104+i]*L[ii+k];
            }    

			/* saturate for numerical stability */
			a = MIN(a, BIGM);
			a = MAX(a, -BIGM); 

			A[j*104+i] = a/L[ii+j];			
        }
        ii += ++di;
    }
}


/**
 * Compute L = L - A*A', where L is lower triangular of size 104
 * and A is a dense matrix of size [104 x 208] in column
 * storage format.
 * 
 * THIS ONE HAS THE WORST ACCES PATTERN POSSIBLE. 
 * POSSIBKE FIX: PUT A INTO ROW MAJOR FORMAT FIRST.
 * 
 */
void beliefPenaltyMPC_LA_DENSE_MMTSUB_104_208(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *L)
{
    int i, j, k, ii, di;
    beliefPenaltyMPC_FLOAT ltemp;
    
    ii = 0; di = 0;
    for( i=0; i<104; i++ ){        
        for( j=0; j<=i; j++ ){
            ltemp = 0; 
            for( k=0; k<208; k++ ){
                ltemp += A[k*104+i]*A[k*104+j];
            }						
            L[ii+j] -= ltemp;
        }
        ii += ++di;
    }
}


/**
 * Cholesky factorization as above, but working on a matrix in 
 * lower triangular storage format of size 104 and outputting
 * the Cholesky factor to matrix L in lower triangular format.
 */
void beliefPenaltyMPC_LA_DENSE_CHOL_104(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *L)
{
    int i, j, k, di, dj;
	 int ii, jj;

    beliefPenaltyMPC_FLOAT l;
    beliefPenaltyMPC_FLOAT Mii;

	/* copy A to L first and then operate on L */
	/* COULD BE OPTIMIZED */
	ii=0; di=0;
	for( i=0; i<104; i++ ){
		for( j=0; j<=i; j++ ){
			L[ii+j] = A[ii+j];
		}
		ii += ++di;
	}    
	
	/* factor L */
	ii=0; di=0;
    for( i=0; i<104; i++ ){
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
        for( j=i+1; j<104; j++ ){
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
void beliefPenaltyMPC_LA_DENSE_MVMSUB1_104_208(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *b, beliefPenaltyMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;

	for( i=0; i<104; i++ ){
		r[i] = b[i] - A[k++]*x[0];
	}	
	for( j=1; j<208; j++ ){		
		for( i=0; i<104; i++ ){
			r[i] -= A[k++]*x[j];
		}
	}
}


/**
 * Forward substitution to solve L*y = b where L is a
 * lower triangular matrix in triangular storage format.
 * 
 * The dimensions involved are 104.
 */
void beliefPenaltyMPC_LA_DENSE_FORWARDSUB_104(beliefPenaltyMPC_FLOAT *L, beliefPenaltyMPC_FLOAT *b, beliefPenaltyMPC_FLOAT *y)
{
    int i,j,ii,di;
    beliefPenaltyMPC_FLOAT yel;
            
    ii = 0; di = 0;
    for( i=0; i<104; i++ ){
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
 * where A is to be computed and is of size [104 x 104],
 * B is given and of size [104 x 104], L is a lower tri-
 * angular matrix of size 104 stored in lower triangular 
 * storage format. Note the transpose of L AND B!
 *
 * Result: A in column major storage format.
 *
 */
void beliefPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_104_104(beliefPenaltyMPC_FLOAT *L, beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *A)
{
    int i,j,k,ii,di;
    beliefPenaltyMPC_FLOAT a;
    
    ii=0; di=0;
    for( j=0; j<104; j++ ){        
        for( i=0; i<104; i++ ){
            a = B[i*104+j];
            for( k=0; k<j; k++ ){
                a -= A[k*104+i]*L[ii+k];
            }    

			/* saturate for numerical stability */
			a = MIN(a, BIGM);
			a = MAX(a, -BIGM); 

			A[j*104+i] = a/L[ii+j];			
        }
        ii += ++di;
    }
}


/**
 * Compute L = L - A*A', where L is lower triangular of size 104
 * and A is a dense matrix of size [104 x 104] in column
 * storage format.
 * 
 * THIS ONE HAS THE WORST ACCES PATTERN POSSIBLE. 
 * POSSIBKE FIX: PUT A INTO ROW MAJOR FORMAT FIRST.
 * 
 */
void beliefPenaltyMPC_LA_DENSE_MMTSUB_104_104(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *L)
{
    int i, j, k, ii, di;
    beliefPenaltyMPC_FLOAT ltemp;
    
    ii = 0; di = 0;
    for( i=0; i<104; i++ ){        
        for( j=0; j<=i; j++ ){
            ltemp = 0; 
            for( k=0; k<104; k++ ){
                ltemp += A[k*104+i]*A[k*104+j];
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
void beliefPenaltyMPC_LA_DENSE_MVMSUB1_104_104(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *b, beliefPenaltyMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;

	for( i=0; i<104; i++ ){
		r[i] = b[i] - A[k++]*x[0];
	}	
	for( j=1; j<104; j++ ){		
		for( i=0; i<104; i++ ){
			r[i] -= A[k++]*x[j];
		}
	}
}


/**
 * Backward Substitution to solve L^T*x = y where L is a
 * lower triangular matrix in triangular storage format.
 * 
 * All involved dimensions are 104.
 */
void beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_104(beliefPenaltyMPC_FLOAT *L, beliefPenaltyMPC_FLOAT *y, beliefPenaltyMPC_FLOAT *x)
{
    int i, ii, di, j, jj, dj;
    beliefPenaltyMPC_FLOAT xel;    
	int start = 5356;
    
    /* now solve L^T*x = y by backward substitution */
    ii = start; di = 103;
    for( i=103; i>=0; i-- ){        
        xel = y[i];        
        jj = start; dj = 103;
        for( j=103; j>i; j-- ){
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
 * Matrix vector multiplication y = b - M'*x where M is of size [104 x 104]
 * and stored in column major format. Note the transpose of M!
 */
void beliefPenaltyMPC_LA_DENSE_MTVMSUB_104_104(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *b, beliefPenaltyMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0; 
	for( i=0; i<104; i++ ){
		r[i] = b[i];
		for( j=0; j<104; j++ ){
			r[i] -= A[k++]*x[j];
		}
	}
}


/*
 * Matrix vector multiplication y = b - M'*x where M is of size [104 x 208]
 * and stored in column major format. Note the transpose of M!
 */
void beliefPenaltyMPC_LA_DENSE_MTVMSUB_104_208(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *b, beliefPenaltyMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0; 
	for( i=0; i<208; i++ ){
		r[i] = b[i];
		for( j=0; j<104; j++ ){
			r[i] -= A[k++]*x[j];
		}
	}
}


/**
 * Backward Substitution to solve L^T*x = y where L is a
 * lower triangular matrix in triangular storage format.
 * 
 * All involved dimensions are 208.
 */
void beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_208(beliefPenaltyMPC_FLOAT *L, beliefPenaltyMPC_FLOAT *y, beliefPenaltyMPC_FLOAT *x)
{
    int i, ii, di, j, jj, dj;
    beliefPenaltyMPC_FLOAT xel;    
	int start = 21528;
    
    /* now solve L^T*x = y by backward substitution */
    ii = start; di = 207;
    for( i=207; i>=0; i-- ){        
        xel = y[i];        
        jj = start; dj = 207;
        for( j=207; j>i; j-- ){
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
 * Vector subtraction z = -x - y for vectors of length 3029.
 */
void beliefPenaltyMPC_LA_VSUB2_3029(beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *y, beliefPenaltyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<3029; i++){
		z[i] = -x[i] - y[i];
	}
}


/**
 * Forward-Backward-Substitution to solve L*L^T*x = b where L is a
 * diagonal matrix of size 325 in vector
 * storage format.
 */
void beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_325(beliefPenaltyMPC_FLOAT *L, beliefPenaltyMPC_FLOAT *b, beliefPenaltyMPC_FLOAT *x)
{
    int i;
            
    /* solve Ly = b by forward and backward substitution */
    for( i=0; i<325; i++ ){
		x[i] = b[i]/(L[i]*L[i]);
    }
    
}


/**
 * Forward-Backward-Substitution to solve L*L^T*x = b where L is a
 * diagonal matrix of size 104 in vector
 * storage format.
 */
void beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_104(beliefPenaltyMPC_FLOAT *L, beliefPenaltyMPC_FLOAT *b, beliefPenaltyMPC_FLOAT *x)
{
    int i;
            
    /* solve Ly = b by forward and backward substitution */
    for( i=0; i<104; i++ ){
		x[i] = b[i]/(L[i]*L[i]);
    }
    
}


/*
 * Vector subtraction z = x(xidx) - y where y, z and xidx are of length 325,
 * and x has length 325 and is indexed through yidx.
 */
void beliefPenaltyMPC_LA_VSUB_INDEXED_325(beliefPenaltyMPC_FLOAT *x, int* xidx, beliefPenaltyMPC_FLOAT *y, beliefPenaltyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<325; i++){
		z[i] = x[xidx[i]] - y[i];
	}
}


/*
 * Vector subtraction x = -u.*v - w for vectors of length 325.
 */
void beliefPenaltyMPC_LA_VSUB3_325(beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *v, beliefPenaltyMPC_FLOAT *w, beliefPenaltyMPC_FLOAT *x)
{
	int i;
	for( i=0; i<325; i++){
		x[i] = -u[i]*v[i] - w[i];
	}
}


/*
 * Vector subtraction z = -x - y(yidx) where y is of length 325
 * and z, x and yidx are of length 117.
 */
void beliefPenaltyMPC_LA_VSUB2_INDEXED_117(beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *y, int* yidx, beliefPenaltyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<117; i++){
		z[i] = -x[i] - y[yidx[i]];
	}
}


/*
 * Vector subtraction x = -u.*v - w for vectors of length 117.
 */
void beliefPenaltyMPC_LA_VSUB3_117(beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *v, beliefPenaltyMPC_FLOAT *w, beliefPenaltyMPC_FLOAT *x)
{
	int i;
	for( i=0; i<117; i++){
		x[i] = -u[i]*v[i] - w[i];
	}
}


/*
 * Vector subtraction z = x(xidx) - y where y, z and xidx are of length 104,
 * and x has length 104 and is indexed through yidx.
 */
void beliefPenaltyMPC_LA_VSUB_INDEXED_104(beliefPenaltyMPC_FLOAT *x, int* xidx, beliefPenaltyMPC_FLOAT *y, beliefPenaltyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<104; i++){
		z[i] = x[xidx[i]] - y[i];
	}
}


/*
 * Vector subtraction x = -u.*v - w for vectors of length 104.
 */
void beliefPenaltyMPC_LA_VSUB3_104(beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *v, beliefPenaltyMPC_FLOAT *w, beliefPenaltyMPC_FLOAT *x)
{
	int i;
	for( i=0; i<104; i++){
		x[i] = -u[i]*v[i] - w[i];
	}
}


/*
 * Vector subtraction z = -x - y(yidx) where y is of length 104
 * and z, x and yidx are of length 104.
 */
void beliefPenaltyMPC_LA_VSUB2_INDEXED_104(beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *y, int* yidx, beliefPenaltyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<104; i++){
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
        for( i=0; i<4186; i++ ){
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
        if( i == 4186 ){
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
    *mu_aff = mymu / (beliefPenaltyMPC_FLOAT)4186;
    return lsIt;
}


/*
 * Vector subtraction x = u.*v - a where a is a scalar
*  and x,u,v are vectors of length 4186.
 */
void beliefPenaltyMPC_LA_VSUB5_4186(beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *v, beliefPenaltyMPC_FLOAT a, beliefPenaltyMPC_FLOAT *x)
{
	int i;
	for( i=0; i<4186; i++){
		x[i] = u[i]*v[i] - a;
	}
}


/*
 * Computes x=0; x(uidx) += u/su; x(vidx) -= v/sv where x is of length 325,
 * u, su, uidx are of length 117 and v, sv, vidx are of length 325.
 */
void beliefPenaltyMPC_LA_VSUB6_INDEXED_325_117_325(beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *su, int* uidx, beliefPenaltyMPC_FLOAT *v, beliefPenaltyMPC_FLOAT *sv, int* vidx, beliefPenaltyMPC_FLOAT *x)
{
	int i;
	for( i=0; i<325; i++ ){
		x[i] = 0;
	}
	for( i=0; i<117; i++){
		x[uidx[i]] += u[i]/su[i];
	}
	for( i=0; i<325; i++){
		x[vidx[i]] -= v[i]/sv[i];
	}
}


/* 
 * Computes r = A*x + B*u
 * where A an B are stored in column major format
 */
void beliefPenaltyMPC_LA_DENSE_2MVMADD_208_325_325(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;
	int m = 0;
	int n;

	for( i=0; i<208; i++ ){
		r[i] = A[k++]*x[0] + B[m++]*u[0];
	}	

	for( j=1; j<325; j++ ){		
		for( i=0; i<208; i++ ){
			r[i] += A[k++]*x[j];
		}
	}
	
	for( n=1; n<325; n++ ){
		for( i=0; i<208; i++ ){
			r[i] += B[m++]*u[n];
		}		
	}
}


/* 
 * Computes r = A*x + B*u
 * where A an B are stored in column major format
 */
void beliefPenaltyMPC_LA_DENSE_2MVMADD_104_325_325(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;
	int m = 0;
	int n;

	for( i=0; i<104; i++ ){
		r[i] = A[k++]*x[0] + B[m++]*u[0];
	}	

	for( j=1; j<325; j++ ){		
		for( i=0; i<104; i++ ){
			r[i] += A[k++]*x[j];
		}
	}
	
	for( n=1; n<325; n++ ){
		for( i=0; i<104; i++ ){
			r[i] += B[m++]*u[n];
		}		
	}
}


/*
 * Computes x=0; x(uidx) += u/su; x(vidx) -= v/sv where x is of length 104,
 * u, su, uidx are of length 104 and v, sv, vidx are of length 104.
 */
void beliefPenaltyMPC_LA_VSUB6_INDEXED_104_104_104(beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *su, int* uidx, beliefPenaltyMPC_FLOAT *v, beliefPenaltyMPC_FLOAT *sv, int* vidx, beliefPenaltyMPC_FLOAT *x)
{
	int i;
	for( i=0; i<104; i++ ){
		x[i] = 0;
	}
	for( i=0; i<104; i++){
		x[uidx[i]] += u[i]/su[i];
	}
	for( i=0; i<104; i++){
		x[vidx[i]] -= v[i]/sv[i];
	}
}


/* 
 * Computes r = A*x + B*u
 * where A an B are stored in column major format
 */
void beliefPenaltyMPC_LA_DENSE_2MVMADD_104_325_104(beliefPenaltyMPC_FLOAT *A, beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *B, beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;
	int m = 0;
	int n;

	for( i=0; i<104; i++ ){
		r[i] = A[k++]*x[0] + B[m++]*u[0];
	}	

	for( j=1; j<325; j++ ){		
		for( i=0; i<104; i++ ){
			r[i] += A[k++]*x[j];
		}
	}
	
	for( n=1; n<104; n++ ){
		for( i=0; i<104; i++ ){
			r[i] += B[m++]*u[n];
		}		
	}
}


/*
 * Vector subtraction z = x - y for vectors of length 3029.
 */
void beliefPenaltyMPC_LA_VSUB_3029(beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *y, beliefPenaltyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<3029; i++){
		z[i] = x[i] - y[i];
	}
}


/** 
 * Computes z = -r./s - u.*y(y)
 * where all vectors except of y are of length 325 (length of y >= 325).
 */
void beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_325(beliefPenaltyMPC_FLOAT *r, beliefPenaltyMPC_FLOAT *s, beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *y, int* yidx, beliefPenaltyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<325; i++ ){
		z[i] = -r[i]/s[i] - u[i]*y[yidx[i]];
	}
}


/** 
 * Computes z = -r./s + u.*y(y)
 * where all vectors except of y are of length 117 (length of y >= 117).
 */
void beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_117(beliefPenaltyMPC_FLOAT *r, beliefPenaltyMPC_FLOAT *s, beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *y, int* yidx, beliefPenaltyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<117; i++ ){
		z[i] = -r[i]/s[i] + u[i]*y[yidx[i]];
	}
}


/** 
 * Computes z = -r./s - u.*y(y)
 * where all vectors except of y are of length 104 (length of y >= 104).
 */
void beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_104(beliefPenaltyMPC_FLOAT *r, beliefPenaltyMPC_FLOAT *s, beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *y, int* yidx, beliefPenaltyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<104; i++ ){
		z[i] = -r[i]/s[i] - u[i]*y[yidx[i]];
	}
}


/** 
 * Computes z = -r./s + u.*y(y)
 * where all vectors except of y are of length 104 (length of y >= 104).
 */
void beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_104(beliefPenaltyMPC_FLOAT *r, beliefPenaltyMPC_FLOAT *s, beliefPenaltyMPC_FLOAT *u, beliefPenaltyMPC_FLOAT *y, int* yidx, beliefPenaltyMPC_FLOAT *z)
{
	int i;
	for( i=0; i<104; i++ ){
		z[i] = -r[i]/s[i] + u[i]*y[yidx[i]];
	}
}


/*
 * Computes ds = -l.\(r + s.*dl) for vectors of length 4186.
 */
void beliefPenaltyMPC_LA_VSUB7_4186(beliefPenaltyMPC_FLOAT *l, beliefPenaltyMPC_FLOAT *r, beliefPenaltyMPC_FLOAT *s, beliefPenaltyMPC_FLOAT *dl, beliefPenaltyMPC_FLOAT *ds)
{
	int i;
	for( i=0; i<4186; i++){
		ds[i] = -(r[i] + s[i]*dl[i])/l[i];
	}
}


/*
 * Vector addition x = x + y for vectors of length 3029.
 */
void beliefPenaltyMPC_LA_VADD_3029(beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *y)
{
	int i;
	for( i=0; i<3029; i++){
		x[i] += y[i];
	}
}


/*
 * Vector addition x = x + y for vectors of length 1040.
 */
void beliefPenaltyMPC_LA_VADD_1040(beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *y)
{
	int i;
	for( i=0; i<1040; i++){
		x[i] += y[i];
	}
}


/*
 * Vector addition x = x + y for vectors of length 4186.
 */
void beliefPenaltyMPC_LA_VADD_4186(beliefPenaltyMPC_FLOAT *x, beliefPenaltyMPC_FLOAT *y)
{
	int i;
	for( i=0; i<4186; i++){
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
        for( i=0; i<4186; i++ ){
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
        if( i == 4186 ){
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
    for( i=0; i<3029; i++ ){
        z[i] += a_gamma*dz[i];
    }
    
    /* equality constraint multipliers */
    for( i=0; i<1040; i++ ){
        v[i] += a_gamma*dv[i];
    }
    
    /* inequality constraint multipliers & slacks, also update mu */
    *mu = 0;
    for( i=0; i<4186; i++ ){
        dltemp = l[i] + a_gamma*dl[i]; l[i] = dltemp;
        dstemp = s[i] + a_gamma*ds[i]; s[i] = dstemp;
        *mu += dltemp*dstemp;
    }
    
    *a = a_gamma;
    *mu /= (beliefPenaltyMPC_FLOAT)4186;
    return lsIt;
}




/* VARIABLE DEFINITIONS ------------------------------------------------ */
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_z[3029];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_v[1040];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_dz_aff[3029];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_dv_aff[1040];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_grad_cost[3029];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_grad_eq[3029];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_rd[3029];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_l[4186];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_s[4186];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_lbys[4186];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_dl_aff[4186];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ds_aff[4186];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_dz_cc[3029];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_dv_cc[1040];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_dl_cc[4186];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ds_cc[4186];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ccrhs[4186];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_grad_ineq[3029];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_z0 = beliefPenaltyMPC_z + 0;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzaff0 = beliefPenaltyMPC_dz_aff + 0;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzcc0 = beliefPenaltyMPC_dz_cc + 0;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_rd0 = beliefPenaltyMPC_rd + 0;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lbyrd0[325];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_cost0 = beliefPenaltyMPC_grad_cost + 0;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_eq0 = beliefPenaltyMPC_grad_eq + 0;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_ineq0 = beliefPenaltyMPC_grad_ineq + 0;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ctv0[325];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_v0 = beliefPenaltyMPC_v + 0;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_re0[208];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_beta0[208];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_betacc0[208];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvaff0 = beliefPenaltyMPC_dv_aff + 0;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvcc0 = beliefPenaltyMPC_dv_cc + 0;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_V0[67600];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Yd0[21736];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ld0[21736];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_yy0[208];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_bmy0[208];
int beliefPenaltyMPC_lbIdx0[325] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255, 256, 257, 258, 259, 260, 261, 262, 263, 264, 265, 266, 267, 268, 269, 270, 271, 272, 273, 274, 275, 276, 277, 278, 279, 280, 281, 282, 283, 284, 285, 286, 287, 288, 289, 290, 291, 292, 293, 294, 295, 296, 297, 298, 299, 300, 301, 302, 303, 304, 305, 306, 307, 308, 309, 310, 311, 312, 313, 314, 315, 316, 317, 318, 319, 320, 321, 322, 323, 324};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llb0 = beliefPenaltyMPC_l + 0;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_slb0 = beliefPenaltyMPC_s + 0;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llbbyslb0 = beliefPenaltyMPC_lbys + 0;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_rilb0[325];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbaff0 = beliefPenaltyMPC_dl_aff + 0;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbaff0 = beliefPenaltyMPC_ds_aff + 0;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbcc0 = beliefPenaltyMPC_dl_cc + 0;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbcc0 = beliefPenaltyMPC_ds_cc + 0;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsl0 = beliefPenaltyMPC_ccrhs + 0;
int beliefPenaltyMPC_ubIdx0[117] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lub0 = beliefPenaltyMPC_l + 325;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_sub0 = beliefPenaltyMPC_s + 325;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lubbysub0 = beliefPenaltyMPC_lbys + 325;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_riub0[117];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubaff0 = beliefPenaltyMPC_dl_aff + 325;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubaff0 = beliefPenaltyMPC_ds_aff + 325;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubcc0 = beliefPenaltyMPC_dl_cc + 325;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubcc0 = beliefPenaltyMPC_ds_cc + 325;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsub0 = beliefPenaltyMPC_ccrhs + 325;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Phi0[325];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_z1 = beliefPenaltyMPC_z + 325;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzaff1 = beliefPenaltyMPC_dz_aff + 325;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzcc1 = beliefPenaltyMPC_dz_cc + 325;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_rd1 = beliefPenaltyMPC_rd + 325;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lbyrd1[325];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_cost1 = beliefPenaltyMPC_grad_cost + 325;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_eq1 = beliefPenaltyMPC_grad_eq + 325;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_ineq1 = beliefPenaltyMPC_grad_ineq + 325;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ctv1[325];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_v1 = beliefPenaltyMPC_v + 208;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_re1[104];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_beta1[104];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_betacc1[104];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvaff1 = beliefPenaltyMPC_dv_aff + 208;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvcc1 = beliefPenaltyMPC_dv_cc + 208;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_V1[33800];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Yd1[5460];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ld1[5460];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_yy1[104];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_bmy1[104];
int beliefPenaltyMPC_lbIdx1[325] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255, 256, 257, 258, 259, 260, 261, 262, 263, 264, 265, 266, 267, 268, 269, 270, 271, 272, 273, 274, 275, 276, 277, 278, 279, 280, 281, 282, 283, 284, 285, 286, 287, 288, 289, 290, 291, 292, 293, 294, 295, 296, 297, 298, 299, 300, 301, 302, 303, 304, 305, 306, 307, 308, 309, 310, 311, 312, 313, 314, 315, 316, 317, 318, 319, 320, 321, 322, 323, 324};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llb1 = beliefPenaltyMPC_l + 442;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_slb1 = beliefPenaltyMPC_s + 442;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llbbyslb1 = beliefPenaltyMPC_lbys + 442;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_rilb1[325];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbaff1 = beliefPenaltyMPC_dl_aff + 442;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbaff1 = beliefPenaltyMPC_ds_aff + 442;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbcc1 = beliefPenaltyMPC_dl_cc + 442;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbcc1 = beliefPenaltyMPC_ds_cc + 442;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsl1 = beliefPenaltyMPC_ccrhs + 442;
int beliefPenaltyMPC_ubIdx1[117] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lub1 = beliefPenaltyMPC_l + 767;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_sub1 = beliefPenaltyMPC_s + 767;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lubbysub1 = beliefPenaltyMPC_lbys + 767;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_riub1[117];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubaff1 = beliefPenaltyMPC_dl_aff + 767;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubaff1 = beliefPenaltyMPC_ds_aff + 767;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubcc1 = beliefPenaltyMPC_dl_cc + 767;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubcc1 = beliefPenaltyMPC_ds_cc + 767;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsub1 = beliefPenaltyMPC_ccrhs + 767;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Phi1[325];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_W1[67600];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ysd1[21632];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lsd1[21632];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_z2 = beliefPenaltyMPC_z + 650;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzaff2 = beliefPenaltyMPC_dz_aff + 650;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzcc2 = beliefPenaltyMPC_dz_cc + 650;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_rd2 = beliefPenaltyMPC_rd + 650;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lbyrd2[325];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_cost2 = beliefPenaltyMPC_grad_cost + 650;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_eq2 = beliefPenaltyMPC_grad_eq + 650;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_ineq2 = beliefPenaltyMPC_grad_ineq + 650;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ctv2[325];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_v2 = beliefPenaltyMPC_v + 312;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_re2[104];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_beta2[104];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_betacc2[104];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvaff2 = beliefPenaltyMPC_dv_aff + 312;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvcc2 = beliefPenaltyMPC_dv_cc + 312;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_V2[33800];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Yd2[5460];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ld2[5460];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_yy2[104];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_bmy2[104];
int beliefPenaltyMPC_lbIdx2[325] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255, 256, 257, 258, 259, 260, 261, 262, 263, 264, 265, 266, 267, 268, 269, 270, 271, 272, 273, 274, 275, 276, 277, 278, 279, 280, 281, 282, 283, 284, 285, 286, 287, 288, 289, 290, 291, 292, 293, 294, 295, 296, 297, 298, 299, 300, 301, 302, 303, 304, 305, 306, 307, 308, 309, 310, 311, 312, 313, 314, 315, 316, 317, 318, 319, 320, 321, 322, 323, 324};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llb2 = beliefPenaltyMPC_l + 884;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_slb2 = beliefPenaltyMPC_s + 884;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llbbyslb2 = beliefPenaltyMPC_lbys + 884;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_rilb2[325];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbaff2 = beliefPenaltyMPC_dl_aff + 884;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbaff2 = beliefPenaltyMPC_ds_aff + 884;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbcc2 = beliefPenaltyMPC_dl_cc + 884;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbcc2 = beliefPenaltyMPC_ds_cc + 884;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsl2 = beliefPenaltyMPC_ccrhs + 884;
int beliefPenaltyMPC_ubIdx2[117] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lub2 = beliefPenaltyMPC_l + 1209;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_sub2 = beliefPenaltyMPC_s + 1209;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lubbysub2 = beliefPenaltyMPC_lbys + 1209;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_riub2[117];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubaff2 = beliefPenaltyMPC_dl_aff + 1209;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubaff2 = beliefPenaltyMPC_ds_aff + 1209;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubcc2 = beliefPenaltyMPC_dl_cc + 1209;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubcc2 = beliefPenaltyMPC_ds_cc + 1209;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsub2 = beliefPenaltyMPC_ccrhs + 1209;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Phi2[325];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_W2[33800];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ysd2[10816];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lsd2[10816];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_z3 = beliefPenaltyMPC_z + 975;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzaff3 = beliefPenaltyMPC_dz_aff + 975;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzcc3 = beliefPenaltyMPC_dz_cc + 975;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_rd3 = beliefPenaltyMPC_rd + 975;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lbyrd3[325];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_cost3 = beliefPenaltyMPC_grad_cost + 975;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_eq3 = beliefPenaltyMPC_grad_eq + 975;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_ineq3 = beliefPenaltyMPC_grad_ineq + 975;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ctv3[325];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_v3 = beliefPenaltyMPC_v + 416;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_re3[104];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_beta3[104];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_betacc3[104];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvaff3 = beliefPenaltyMPC_dv_aff + 416;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvcc3 = beliefPenaltyMPC_dv_cc + 416;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_V3[33800];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Yd3[5460];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ld3[5460];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_yy3[104];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_bmy3[104];
int beliefPenaltyMPC_lbIdx3[325] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255, 256, 257, 258, 259, 260, 261, 262, 263, 264, 265, 266, 267, 268, 269, 270, 271, 272, 273, 274, 275, 276, 277, 278, 279, 280, 281, 282, 283, 284, 285, 286, 287, 288, 289, 290, 291, 292, 293, 294, 295, 296, 297, 298, 299, 300, 301, 302, 303, 304, 305, 306, 307, 308, 309, 310, 311, 312, 313, 314, 315, 316, 317, 318, 319, 320, 321, 322, 323, 324};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llb3 = beliefPenaltyMPC_l + 1326;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_slb3 = beliefPenaltyMPC_s + 1326;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llbbyslb3 = beliefPenaltyMPC_lbys + 1326;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_rilb3[325];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbaff3 = beliefPenaltyMPC_dl_aff + 1326;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbaff3 = beliefPenaltyMPC_ds_aff + 1326;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbcc3 = beliefPenaltyMPC_dl_cc + 1326;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbcc3 = beliefPenaltyMPC_ds_cc + 1326;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsl3 = beliefPenaltyMPC_ccrhs + 1326;
int beliefPenaltyMPC_ubIdx3[117] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lub3 = beliefPenaltyMPC_l + 1651;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_sub3 = beliefPenaltyMPC_s + 1651;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lubbysub3 = beliefPenaltyMPC_lbys + 1651;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_riub3[117];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubaff3 = beliefPenaltyMPC_dl_aff + 1651;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubaff3 = beliefPenaltyMPC_ds_aff + 1651;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubcc3 = beliefPenaltyMPC_dl_cc + 1651;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubcc3 = beliefPenaltyMPC_ds_cc + 1651;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsub3 = beliefPenaltyMPC_ccrhs + 1651;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Phi3[325];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_W3[33800];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ysd3[10816];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lsd3[10816];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_z4 = beliefPenaltyMPC_z + 1300;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzaff4 = beliefPenaltyMPC_dz_aff + 1300;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzcc4 = beliefPenaltyMPC_dz_cc + 1300;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_rd4 = beliefPenaltyMPC_rd + 1300;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lbyrd4[325];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_cost4 = beliefPenaltyMPC_grad_cost + 1300;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_eq4 = beliefPenaltyMPC_grad_eq + 1300;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_ineq4 = beliefPenaltyMPC_grad_ineq + 1300;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ctv4[325];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_v4 = beliefPenaltyMPC_v + 520;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_re4[104];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_beta4[104];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_betacc4[104];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvaff4 = beliefPenaltyMPC_dv_aff + 520;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvcc4 = beliefPenaltyMPC_dv_cc + 520;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_V4[33800];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Yd4[5460];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ld4[5460];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_yy4[104];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_bmy4[104];
int beliefPenaltyMPC_lbIdx4[325] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255, 256, 257, 258, 259, 260, 261, 262, 263, 264, 265, 266, 267, 268, 269, 270, 271, 272, 273, 274, 275, 276, 277, 278, 279, 280, 281, 282, 283, 284, 285, 286, 287, 288, 289, 290, 291, 292, 293, 294, 295, 296, 297, 298, 299, 300, 301, 302, 303, 304, 305, 306, 307, 308, 309, 310, 311, 312, 313, 314, 315, 316, 317, 318, 319, 320, 321, 322, 323, 324};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llb4 = beliefPenaltyMPC_l + 1768;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_slb4 = beliefPenaltyMPC_s + 1768;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llbbyslb4 = beliefPenaltyMPC_lbys + 1768;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_rilb4[325];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbaff4 = beliefPenaltyMPC_dl_aff + 1768;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbaff4 = beliefPenaltyMPC_ds_aff + 1768;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbcc4 = beliefPenaltyMPC_dl_cc + 1768;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbcc4 = beliefPenaltyMPC_ds_cc + 1768;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsl4 = beliefPenaltyMPC_ccrhs + 1768;
int beliefPenaltyMPC_ubIdx4[117] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lub4 = beliefPenaltyMPC_l + 2093;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_sub4 = beliefPenaltyMPC_s + 2093;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lubbysub4 = beliefPenaltyMPC_lbys + 2093;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_riub4[117];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubaff4 = beliefPenaltyMPC_dl_aff + 2093;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubaff4 = beliefPenaltyMPC_ds_aff + 2093;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubcc4 = beliefPenaltyMPC_dl_cc + 2093;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubcc4 = beliefPenaltyMPC_ds_cc + 2093;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsub4 = beliefPenaltyMPC_ccrhs + 2093;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Phi4[325];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_W4[33800];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ysd4[10816];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lsd4[10816];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_z5 = beliefPenaltyMPC_z + 1625;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzaff5 = beliefPenaltyMPC_dz_aff + 1625;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzcc5 = beliefPenaltyMPC_dz_cc + 1625;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_rd5 = beliefPenaltyMPC_rd + 1625;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lbyrd5[325];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_cost5 = beliefPenaltyMPC_grad_cost + 1625;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_eq5 = beliefPenaltyMPC_grad_eq + 1625;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_ineq5 = beliefPenaltyMPC_grad_ineq + 1625;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ctv5[325];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_v5 = beliefPenaltyMPC_v + 624;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_re5[104];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_beta5[104];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_betacc5[104];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvaff5 = beliefPenaltyMPC_dv_aff + 624;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvcc5 = beliefPenaltyMPC_dv_cc + 624;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_V5[33800];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Yd5[5460];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ld5[5460];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_yy5[104];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_bmy5[104];
int beliefPenaltyMPC_lbIdx5[325] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255, 256, 257, 258, 259, 260, 261, 262, 263, 264, 265, 266, 267, 268, 269, 270, 271, 272, 273, 274, 275, 276, 277, 278, 279, 280, 281, 282, 283, 284, 285, 286, 287, 288, 289, 290, 291, 292, 293, 294, 295, 296, 297, 298, 299, 300, 301, 302, 303, 304, 305, 306, 307, 308, 309, 310, 311, 312, 313, 314, 315, 316, 317, 318, 319, 320, 321, 322, 323, 324};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llb5 = beliefPenaltyMPC_l + 2210;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_slb5 = beliefPenaltyMPC_s + 2210;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llbbyslb5 = beliefPenaltyMPC_lbys + 2210;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_rilb5[325];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbaff5 = beliefPenaltyMPC_dl_aff + 2210;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbaff5 = beliefPenaltyMPC_ds_aff + 2210;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbcc5 = beliefPenaltyMPC_dl_cc + 2210;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbcc5 = beliefPenaltyMPC_ds_cc + 2210;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsl5 = beliefPenaltyMPC_ccrhs + 2210;
int beliefPenaltyMPC_ubIdx5[117] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lub5 = beliefPenaltyMPC_l + 2535;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_sub5 = beliefPenaltyMPC_s + 2535;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lubbysub5 = beliefPenaltyMPC_lbys + 2535;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_riub5[117];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubaff5 = beliefPenaltyMPC_dl_aff + 2535;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubaff5 = beliefPenaltyMPC_ds_aff + 2535;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubcc5 = beliefPenaltyMPC_dl_cc + 2535;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubcc5 = beliefPenaltyMPC_ds_cc + 2535;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsub5 = beliefPenaltyMPC_ccrhs + 2535;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Phi5[325];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_W5[33800];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ysd5[10816];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lsd5[10816];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_z6 = beliefPenaltyMPC_z + 1950;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzaff6 = beliefPenaltyMPC_dz_aff + 1950;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzcc6 = beliefPenaltyMPC_dz_cc + 1950;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_rd6 = beliefPenaltyMPC_rd + 1950;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lbyrd6[325];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_cost6 = beliefPenaltyMPC_grad_cost + 1950;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_eq6 = beliefPenaltyMPC_grad_eq + 1950;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_ineq6 = beliefPenaltyMPC_grad_ineq + 1950;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ctv6[325];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_v6 = beliefPenaltyMPC_v + 728;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_re6[104];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_beta6[104];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_betacc6[104];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvaff6 = beliefPenaltyMPC_dv_aff + 728;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvcc6 = beliefPenaltyMPC_dv_cc + 728;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_V6[33800];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Yd6[5460];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ld6[5460];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_yy6[104];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_bmy6[104];
int beliefPenaltyMPC_lbIdx6[325] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255, 256, 257, 258, 259, 260, 261, 262, 263, 264, 265, 266, 267, 268, 269, 270, 271, 272, 273, 274, 275, 276, 277, 278, 279, 280, 281, 282, 283, 284, 285, 286, 287, 288, 289, 290, 291, 292, 293, 294, 295, 296, 297, 298, 299, 300, 301, 302, 303, 304, 305, 306, 307, 308, 309, 310, 311, 312, 313, 314, 315, 316, 317, 318, 319, 320, 321, 322, 323, 324};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llb6 = beliefPenaltyMPC_l + 2652;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_slb6 = beliefPenaltyMPC_s + 2652;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llbbyslb6 = beliefPenaltyMPC_lbys + 2652;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_rilb6[325];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbaff6 = beliefPenaltyMPC_dl_aff + 2652;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbaff6 = beliefPenaltyMPC_ds_aff + 2652;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbcc6 = beliefPenaltyMPC_dl_cc + 2652;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbcc6 = beliefPenaltyMPC_ds_cc + 2652;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsl6 = beliefPenaltyMPC_ccrhs + 2652;
int beliefPenaltyMPC_ubIdx6[117] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lub6 = beliefPenaltyMPC_l + 2977;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_sub6 = beliefPenaltyMPC_s + 2977;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lubbysub6 = beliefPenaltyMPC_lbys + 2977;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_riub6[117];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubaff6 = beliefPenaltyMPC_dl_aff + 2977;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubaff6 = beliefPenaltyMPC_ds_aff + 2977;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubcc6 = beliefPenaltyMPC_dl_cc + 2977;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubcc6 = beliefPenaltyMPC_ds_cc + 2977;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsub6 = beliefPenaltyMPC_ccrhs + 2977;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Phi6[325];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_W6[33800];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ysd6[10816];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lsd6[10816];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_z7 = beliefPenaltyMPC_z + 2275;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzaff7 = beliefPenaltyMPC_dz_aff + 2275;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzcc7 = beliefPenaltyMPC_dz_cc + 2275;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_rd7 = beliefPenaltyMPC_rd + 2275;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lbyrd7[325];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_cost7 = beliefPenaltyMPC_grad_cost + 2275;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_eq7 = beliefPenaltyMPC_grad_eq + 2275;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_ineq7 = beliefPenaltyMPC_grad_ineq + 2275;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ctv7[325];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_v7 = beliefPenaltyMPC_v + 832;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_re7[104];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_beta7[104];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_betacc7[104];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvaff7 = beliefPenaltyMPC_dv_aff + 832;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvcc7 = beliefPenaltyMPC_dv_cc + 832;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_V7[33800];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Yd7[5460];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ld7[5460];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_yy7[104];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_bmy7[104];
int beliefPenaltyMPC_lbIdx7[325] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255, 256, 257, 258, 259, 260, 261, 262, 263, 264, 265, 266, 267, 268, 269, 270, 271, 272, 273, 274, 275, 276, 277, 278, 279, 280, 281, 282, 283, 284, 285, 286, 287, 288, 289, 290, 291, 292, 293, 294, 295, 296, 297, 298, 299, 300, 301, 302, 303, 304, 305, 306, 307, 308, 309, 310, 311, 312, 313, 314, 315, 316, 317, 318, 319, 320, 321, 322, 323, 324};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llb7 = beliefPenaltyMPC_l + 3094;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_slb7 = beliefPenaltyMPC_s + 3094;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llbbyslb7 = beliefPenaltyMPC_lbys + 3094;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_rilb7[325];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbaff7 = beliefPenaltyMPC_dl_aff + 3094;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbaff7 = beliefPenaltyMPC_ds_aff + 3094;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbcc7 = beliefPenaltyMPC_dl_cc + 3094;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbcc7 = beliefPenaltyMPC_ds_cc + 3094;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsl7 = beliefPenaltyMPC_ccrhs + 3094;
int beliefPenaltyMPC_ubIdx7[117] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lub7 = beliefPenaltyMPC_l + 3419;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_sub7 = beliefPenaltyMPC_s + 3419;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lubbysub7 = beliefPenaltyMPC_lbys + 3419;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_riub7[117];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubaff7 = beliefPenaltyMPC_dl_aff + 3419;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubaff7 = beliefPenaltyMPC_ds_aff + 3419;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubcc7 = beliefPenaltyMPC_dl_cc + 3419;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubcc7 = beliefPenaltyMPC_ds_cc + 3419;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsub7 = beliefPenaltyMPC_ccrhs + 3419;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Phi7[325];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_W7[33800];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ysd7[10816];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lsd7[10816];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_z8 = beliefPenaltyMPC_z + 2600;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzaff8 = beliefPenaltyMPC_dz_aff + 2600;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzcc8 = beliefPenaltyMPC_dz_cc + 2600;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_rd8 = beliefPenaltyMPC_rd + 2600;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lbyrd8[325];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_cost8 = beliefPenaltyMPC_grad_cost + 2600;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_eq8 = beliefPenaltyMPC_grad_eq + 2600;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_ineq8 = beliefPenaltyMPC_grad_ineq + 2600;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ctv8[325];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_v8 = beliefPenaltyMPC_v + 936;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_re8[104];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_beta8[104];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_betacc8[104];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvaff8 = beliefPenaltyMPC_dv_aff + 936;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dvcc8 = beliefPenaltyMPC_dv_cc + 936;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_V8[33800];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Yd8[5460];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ld8[5460];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_yy8[104];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_bmy8[104];
int beliefPenaltyMPC_lbIdx8[325] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255, 256, 257, 258, 259, 260, 261, 262, 263, 264, 265, 266, 267, 268, 269, 270, 271, 272, 273, 274, 275, 276, 277, 278, 279, 280, 281, 282, 283, 284, 285, 286, 287, 288, 289, 290, 291, 292, 293, 294, 295, 296, 297, 298, 299, 300, 301, 302, 303, 304, 305, 306, 307, 308, 309, 310, 311, 312, 313, 314, 315, 316, 317, 318, 319, 320, 321, 322, 323, 324};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llb8 = beliefPenaltyMPC_l + 3536;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_slb8 = beliefPenaltyMPC_s + 3536;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llbbyslb8 = beliefPenaltyMPC_lbys + 3536;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_rilb8[325];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbaff8 = beliefPenaltyMPC_dl_aff + 3536;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbaff8 = beliefPenaltyMPC_ds_aff + 3536;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbcc8 = beliefPenaltyMPC_dl_cc + 3536;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbcc8 = beliefPenaltyMPC_ds_cc + 3536;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsl8 = beliefPenaltyMPC_ccrhs + 3536;
int beliefPenaltyMPC_ubIdx8[117] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lub8 = beliefPenaltyMPC_l + 3861;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_sub8 = beliefPenaltyMPC_s + 3861;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lubbysub8 = beliefPenaltyMPC_lbys + 3861;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_riub8[117];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubaff8 = beliefPenaltyMPC_dl_aff + 3861;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubaff8 = beliefPenaltyMPC_ds_aff + 3861;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubcc8 = beliefPenaltyMPC_dl_cc + 3861;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubcc8 = beliefPenaltyMPC_ds_cc + 3861;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsub8 = beliefPenaltyMPC_ccrhs + 3861;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Phi8[325];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_W8[33800];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Ysd8[10816];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lsd8[10816];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_f9[104] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_z9 = beliefPenaltyMPC_z + 2925;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzaff9 = beliefPenaltyMPC_dz_aff + 2925;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dzcc9 = beliefPenaltyMPC_dz_cc + 2925;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_rd9 = beliefPenaltyMPC_rd + 2925;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Lbyrd9[104];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_cost9 = beliefPenaltyMPC_grad_cost + 2925;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_eq9 = beliefPenaltyMPC_grad_eq + 2925;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_grad_ineq9 = beliefPenaltyMPC_grad_ineq + 2925;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_ctv9[104];
int beliefPenaltyMPC_lbIdx9[104] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llb9 = beliefPenaltyMPC_l + 3978;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_slb9 = beliefPenaltyMPC_s + 3978;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_llbbyslb9 = beliefPenaltyMPC_lbys + 3978;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_rilb9[104];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbaff9 = beliefPenaltyMPC_dl_aff + 3978;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbaff9 = beliefPenaltyMPC_ds_aff + 3978;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dllbcc9 = beliefPenaltyMPC_dl_cc + 3978;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dslbcc9 = beliefPenaltyMPC_ds_cc + 3978;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsl9 = beliefPenaltyMPC_ccrhs + 3978;
int beliefPenaltyMPC_ubIdx9[104] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103};
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lub9 = beliefPenaltyMPC_l + 4082;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_sub9 = beliefPenaltyMPC_s + 4082;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_lubbysub9 = beliefPenaltyMPC_lbys + 4082;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_riub9[104];
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubaff9 = beliefPenaltyMPC_dl_aff + 4082;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubaff9 = beliefPenaltyMPC_ds_aff + 4082;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dlubcc9 = beliefPenaltyMPC_dl_cc + 4082;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_dsubcc9 = beliefPenaltyMPC_ds_cc + 4082;
beliefPenaltyMPC_FLOAT* beliefPenaltyMPC_ccrhsub9 = beliefPenaltyMPC_ccrhs + 4082;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Phi9[104];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_W9[10816];
beliefPenaltyMPC_FLOAT musigma;
beliefPenaltyMPC_FLOAT sigma_3rdroot;
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Diag1_0[325];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_Diag2_0[325];
beliefPenaltyMPC_FLOAT beliefPenaltyMPC_L_0[52650];




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
beliefPenaltyMPC_LA_INITIALIZEVECTOR_3029(beliefPenaltyMPC_z, 0);
beliefPenaltyMPC_LA_INITIALIZEVECTOR_1040(beliefPenaltyMPC_v, 1);
beliefPenaltyMPC_LA_INITIALIZEVECTOR_4186(beliefPenaltyMPC_l, 1);
beliefPenaltyMPC_LA_INITIALIZEVECTOR_4186(beliefPenaltyMPC_s, 1);
info->mu = 0;
beliefPenaltyMPC_LA_DOTACC_4186(beliefPenaltyMPC_l, beliefPenaltyMPC_s, &info->mu);
info->mu /= 4186;
while( 1 ){
info->pobj = 0;
beliefPenaltyMPC_LA_DIAG_QUADFCN_325(params->H1, params->f1, beliefPenaltyMPC_z0, beliefPenaltyMPC_grad_cost0, &info->pobj);
beliefPenaltyMPC_LA_DIAG_QUADFCN_325(params->H2, params->f2, beliefPenaltyMPC_z1, beliefPenaltyMPC_grad_cost1, &info->pobj);
beliefPenaltyMPC_LA_DIAG_QUADFCN_325(params->H3, params->f3, beliefPenaltyMPC_z2, beliefPenaltyMPC_grad_cost2, &info->pobj);
beliefPenaltyMPC_LA_DIAG_QUADFCN_325(params->H4, params->f4, beliefPenaltyMPC_z3, beliefPenaltyMPC_grad_cost3, &info->pobj);
beliefPenaltyMPC_LA_DIAG_QUADFCN_325(params->H5, params->f5, beliefPenaltyMPC_z4, beliefPenaltyMPC_grad_cost4, &info->pobj);
beliefPenaltyMPC_LA_DIAG_QUADFCN_325(params->H6, params->f6, beliefPenaltyMPC_z5, beliefPenaltyMPC_grad_cost5, &info->pobj);
beliefPenaltyMPC_LA_DIAG_QUADFCN_325(params->H7, params->f7, beliefPenaltyMPC_z6, beliefPenaltyMPC_grad_cost6, &info->pobj);
beliefPenaltyMPC_LA_DIAG_QUADFCN_325(params->H8, params->f8, beliefPenaltyMPC_z7, beliefPenaltyMPC_grad_cost7, &info->pobj);
beliefPenaltyMPC_LA_DIAG_QUADFCN_325(params->H9, params->f9, beliefPenaltyMPC_z8, beliefPenaltyMPC_grad_cost8, &info->pobj);
beliefPenaltyMPC_LA_DIAG_QUADFCN_104(params->H10, beliefPenaltyMPC_f9, beliefPenaltyMPC_z9, beliefPenaltyMPC_grad_cost9, &info->pobj);
info->res_eq = 0;
info->dgap = 0;
beliefPenaltyMPC_LA_DENSE_MVMSUB3_208_325_325(params->C1, beliefPenaltyMPC_z0, params->D2, beliefPenaltyMPC_z1, params->e1, beliefPenaltyMPC_v0, beliefPenaltyMPC_re0, &info->dgap, &info->res_eq);
beliefPenaltyMPC_LA_DENSE_MVMSUB3_104_325_325(params->C2, beliefPenaltyMPC_z1, params->D3, beliefPenaltyMPC_z2, params->e2, beliefPenaltyMPC_v1, beliefPenaltyMPC_re1, &info->dgap, &info->res_eq);
beliefPenaltyMPC_LA_DENSE_MVMSUB3_104_325_325(params->C3, beliefPenaltyMPC_z2, params->D4, beliefPenaltyMPC_z3, params->e3, beliefPenaltyMPC_v2, beliefPenaltyMPC_re2, &info->dgap, &info->res_eq);
beliefPenaltyMPC_LA_DENSE_MVMSUB3_104_325_325(params->C4, beliefPenaltyMPC_z3, params->D5, beliefPenaltyMPC_z4, params->e4, beliefPenaltyMPC_v3, beliefPenaltyMPC_re3, &info->dgap, &info->res_eq);
beliefPenaltyMPC_LA_DENSE_MVMSUB3_104_325_325(params->C5, beliefPenaltyMPC_z4, params->D6, beliefPenaltyMPC_z5, params->e5, beliefPenaltyMPC_v4, beliefPenaltyMPC_re4, &info->dgap, &info->res_eq);
beliefPenaltyMPC_LA_DENSE_MVMSUB3_104_325_325(params->C6, beliefPenaltyMPC_z5, params->D7, beliefPenaltyMPC_z6, params->e6, beliefPenaltyMPC_v5, beliefPenaltyMPC_re5, &info->dgap, &info->res_eq);
beliefPenaltyMPC_LA_DENSE_MVMSUB3_104_325_325(params->C7, beliefPenaltyMPC_z6, params->D8, beliefPenaltyMPC_z7, params->e7, beliefPenaltyMPC_v6, beliefPenaltyMPC_re6, &info->dgap, &info->res_eq);
beliefPenaltyMPC_LA_DENSE_MVMSUB3_104_325_325(params->C8, beliefPenaltyMPC_z7, params->D9, beliefPenaltyMPC_z8, params->e8, beliefPenaltyMPC_v7, beliefPenaltyMPC_re7, &info->dgap, &info->res_eq);
beliefPenaltyMPC_LA_DENSE_MVMSUB3_104_325_104(params->C9, beliefPenaltyMPC_z8, params->D10, beliefPenaltyMPC_z9, params->e9, beliefPenaltyMPC_v8, beliefPenaltyMPC_re8, &info->dgap, &info->res_eq);
beliefPenaltyMPC_LA_DENSE_MTVM_208_325(params->C1, beliefPenaltyMPC_v0, beliefPenaltyMPC_grad_eq0);
beliefPenaltyMPC_LA_DENSE_MTVM2_104_325_208(params->C2, beliefPenaltyMPC_v1, params->D2, beliefPenaltyMPC_v0, beliefPenaltyMPC_grad_eq1);
beliefPenaltyMPC_LA_DENSE_MTVM2_104_325_104(params->C3, beliefPenaltyMPC_v2, params->D3, beliefPenaltyMPC_v1, beliefPenaltyMPC_grad_eq2);
beliefPenaltyMPC_LA_DENSE_MTVM2_104_325_104(params->C4, beliefPenaltyMPC_v3, params->D4, beliefPenaltyMPC_v2, beliefPenaltyMPC_grad_eq3);
beliefPenaltyMPC_LA_DENSE_MTVM2_104_325_104(params->C5, beliefPenaltyMPC_v4, params->D5, beliefPenaltyMPC_v3, beliefPenaltyMPC_grad_eq4);
beliefPenaltyMPC_LA_DENSE_MTVM2_104_325_104(params->C6, beliefPenaltyMPC_v5, params->D6, beliefPenaltyMPC_v4, beliefPenaltyMPC_grad_eq5);
beliefPenaltyMPC_LA_DENSE_MTVM2_104_325_104(params->C7, beliefPenaltyMPC_v6, params->D7, beliefPenaltyMPC_v5, beliefPenaltyMPC_grad_eq6);
beliefPenaltyMPC_LA_DENSE_MTVM2_104_325_104(params->C8, beliefPenaltyMPC_v7, params->D8, beliefPenaltyMPC_v6, beliefPenaltyMPC_grad_eq7);
beliefPenaltyMPC_LA_DENSE_MTVM2_104_325_104(params->C9, beliefPenaltyMPC_v8, params->D9, beliefPenaltyMPC_v7, beliefPenaltyMPC_grad_eq8);
beliefPenaltyMPC_LA_DENSE_MTVM_104_104(params->D10, beliefPenaltyMPC_v8, beliefPenaltyMPC_grad_eq9);
info->res_ineq = 0;
beliefPenaltyMPC_LA_VSUBADD3_325(params->lb1, beliefPenaltyMPC_z0, beliefPenaltyMPC_lbIdx0, beliefPenaltyMPC_llb0, beliefPenaltyMPC_slb0, beliefPenaltyMPC_rilb0, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD2_117(beliefPenaltyMPC_z0, beliefPenaltyMPC_ubIdx0, params->ub1, beliefPenaltyMPC_lub0, beliefPenaltyMPC_sub0, beliefPenaltyMPC_riub0, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD3_325(params->lb2, beliefPenaltyMPC_z1, beliefPenaltyMPC_lbIdx1, beliefPenaltyMPC_llb1, beliefPenaltyMPC_slb1, beliefPenaltyMPC_rilb1, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD2_117(beliefPenaltyMPC_z1, beliefPenaltyMPC_ubIdx1, params->ub2, beliefPenaltyMPC_lub1, beliefPenaltyMPC_sub1, beliefPenaltyMPC_riub1, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD3_325(params->lb3, beliefPenaltyMPC_z2, beliefPenaltyMPC_lbIdx2, beliefPenaltyMPC_llb2, beliefPenaltyMPC_slb2, beliefPenaltyMPC_rilb2, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD2_117(beliefPenaltyMPC_z2, beliefPenaltyMPC_ubIdx2, params->ub3, beliefPenaltyMPC_lub2, beliefPenaltyMPC_sub2, beliefPenaltyMPC_riub2, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD3_325(params->lb4, beliefPenaltyMPC_z3, beliefPenaltyMPC_lbIdx3, beliefPenaltyMPC_llb3, beliefPenaltyMPC_slb3, beliefPenaltyMPC_rilb3, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD2_117(beliefPenaltyMPC_z3, beliefPenaltyMPC_ubIdx3, params->ub4, beliefPenaltyMPC_lub3, beliefPenaltyMPC_sub3, beliefPenaltyMPC_riub3, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD3_325(params->lb5, beliefPenaltyMPC_z4, beliefPenaltyMPC_lbIdx4, beliefPenaltyMPC_llb4, beliefPenaltyMPC_slb4, beliefPenaltyMPC_rilb4, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD2_117(beliefPenaltyMPC_z4, beliefPenaltyMPC_ubIdx4, params->ub5, beliefPenaltyMPC_lub4, beliefPenaltyMPC_sub4, beliefPenaltyMPC_riub4, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD3_325(params->lb6, beliefPenaltyMPC_z5, beliefPenaltyMPC_lbIdx5, beliefPenaltyMPC_llb5, beliefPenaltyMPC_slb5, beliefPenaltyMPC_rilb5, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD2_117(beliefPenaltyMPC_z5, beliefPenaltyMPC_ubIdx5, params->ub6, beliefPenaltyMPC_lub5, beliefPenaltyMPC_sub5, beliefPenaltyMPC_riub5, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD3_325(params->lb7, beliefPenaltyMPC_z6, beliefPenaltyMPC_lbIdx6, beliefPenaltyMPC_llb6, beliefPenaltyMPC_slb6, beliefPenaltyMPC_rilb6, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD2_117(beliefPenaltyMPC_z6, beliefPenaltyMPC_ubIdx6, params->ub7, beliefPenaltyMPC_lub6, beliefPenaltyMPC_sub6, beliefPenaltyMPC_riub6, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD3_325(params->lb8, beliefPenaltyMPC_z7, beliefPenaltyMPC_lbIdx7, beliefPenaltyMPC_llb7, beliefPenaltyMPC_slb7, beliefPenaltyMPC_rilb7, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD2_117(beliefPenaltyMPC_z7, beliefPenaltyMPC_ubIdx7, params->ub8, beliefPenaltyMPC_lub7, beliefPenaltyMPC_sub7, beliefPenaltyMPC_riub7, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD3_325(params->lb9, beliefPenaltyMPC_z8, beliefPenaltyMPC_lbIdx8, beliefPenaltyMPC_llb8, beliefPenaltyMPC_slb8, beliefPenaltyMPC_rilb8, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD2_117(beliefPenaltyMPC_z8, beliefPenaltyMPC_ubIdx8, params->ub9, beliefPenaltyMPC_lub8, beliefPenaltyMPC_sub8, beliefPenaltyMPC_riub8, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD3_104(params->lb10, beliefPenaltyMPC_z9, beliefPenaltyMPC_lbIdx9, beliefPenaltyMPC_llb9, beliefPenaltyMPC_slb9, beliefPenaltyMPC_rilb9, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_VSUBADD2_104(beliefPenaltyMPC_z9, beliefPenaltyMPC_ubIdx9, params->ub10, beliefPenaltyMPC_lub9, beliefPenaltyMPC_sub9, beliefPenaltyMPC_riub9, &info->dgap, &info->res_ineq);
beliefPenaltyMPC_LA_INEQ_B_GRAD_325_325_117(beliefPenaltyMPC_lub0, beliefPenaltyMPC_sub0, beliefPenaltyMPC_riub0, beliefPenaltyMPC_llb0, beliefPenaltyMPC_slb0, beliefPenaltyMPC_rilb0, beliefPenaltyMPC_lbIdx0, beliefPenaltyMPC_ubIdx0, beliefPenaltyMPC_grad_ineq0, beliefPenaltyMPC_lubbysub0, beliefPenaltyMPC_llbbyslb0);
beliefPenaltyMPC_LA_INEQ_B_GRAD_325_325_117(beliefPenaltyMPC_lub1, beliefPenaltyMPC_sub1, beliefPenaltyMPC_riub1, beliefPenaltyMPC_llb1, beliefPenaltyMPC_slb1, beliefPenaltyMPC_rilb1, beliefPenaltyMPC_lbIdx1, beliefPenaltyMPC_ubIdx1, beliefPenaltyMPC_grad_ineq1, beliefPenaltyMPC_lubbysub1, beliefPenaltyMPC_llbbyslb1);
beliefPenaltyMPC_LA_INEQ_B_GRAD_325_325_117(beliefPenaltyMPC_lub2, beliefPenaltyMPC_sub2, beliefPenaltyMPC_riub2, beliefPenaltyMPC_llb2, beliefPenaltyMPC_slb2, beliefPenaltyMPC_rilb2, beliefPenaltyMPC_lbIdx2, beliefPenaltyMPC_ubIdx2, beliefPenaltyMPC_grad_ineq2, beliefPenaltyMPC_lubbysub2, beliefPenaltyMPC_llbbyslb2);
beliefPenaltyMPC_LA_INEQ_B_GRAD_325_325_117(beliefPenaltyMPC_lub3, beliefPenaltyMPC_sub3, beliefPenaltyMPC_riub3, beliefPenaltyMPC_llb3, beliefPenaltyMPC_slb3, beliefPenaltyMPC_rilb3, beliefPenaltyMPC_lbIdx3, beliefPenaltyMPC_ubIdx3, beliefPenaltyMPC_grad_ineq3, beliefPenaltyMPC_lubbysub3, beliefPenaltyMPC_llbbyslb3);
beliefPenaltyMPC_LA_INEQ_B_GRAD_325_325_117(beliefPenaltyMPC_lub4, beliefPenaltyMPC_sub4, beliefPenaltyMPC_riub4, beliefPenaltyMPC_llb4, beliefPenaltyMPC_slb4, beliefPenaltyMPC_rilb4, beliefPenaltyMPC_lbIdx4, beliefPenaltyMPC_ubIdx4, beliefPenaltyMPC_grad_ineq4, beliefPenaltyMPC_lubbysub4, beliefPenaltyMPC_llbbyslb4);
beliefPenaltyMPC_LA_INEQ_B_GRAD_325_325_117(beliefPenaltyMPC_lub5, beliefPenaltyMPC_sub5, beliefPenaltyMPC_riub5, beliefPenaltyMPC_llb5, beliefPenaltyMPC_slb5, beliefPenaltyMPC_rilb5, beliefPenaltyMPC_lbIdx5, beliefPenaltyMPC_ubIdx5, beliefPenaltyMPC_grad_ineq5, beliefPenaltyMPC_lubbysub5, beliefPenaltyMPC_llbbyslb5);
beliefPenaltyMPC_LA_INEQ_B_GRAD_325_325_117(beliefPenaltyMPC_lub6, beliefPenaltyMPC_sub6, beliefPenaltyMPC_riub6, beliefPenaltyMPC_llb6, beliefPenaltyMPC_slb6, beliefPenaltyMPC_rilb6, beliefPenaltyMPC_lbIdx6, beliefPenaltyMPC_ubIdx6, beliefPenaltyMPC_grad_ineq6, beliefPenaltyMPC_lubbysub6, beliefPenaltyMPC_llbbyslb6);
beliefPenaltyMPC_LA_INEQ_B_GRAD_325_325_117(beliefPenaltyMPC_lub7, beliefPenaltyMPC_sub7, beliefPenaltyMPC_riub7, beliefPenaltyMPC_llb7, beliefPenaltyMPC_slb7, beliefPenaltyMPC_rilb7, beliefPenaltyMPC_lbIdx7, beliefPenaltyMPC_ubIdx7, beliefPenaltyMPC_grad_ineq7, beliefPenaltyMPC_lubbysub7, beliefPenaltyMPC_llbbyslb7);
beliefPenaltyMPC_LA_INEQ_B_GRAD_325_325_117(beliefPenaltyMPC_lub8, beliefPenaltyMPC_sub8, beliefPenaltyMPC_riub8, beliefPenaltyMPC_llb8, beliefPenaltyMPC_slb8, beliefPenaltyMPC_rilb8, beliefPenaltyMPC_lbIdx8, beliefPenaltyMPC_ubIdx8, beliefPenaltyMPC_grad_ineq8, beliefPenaltyMPC_lubbysub8, beliefPenaltyMPC_llbbyslb8);
beliefPenaltyMPC_LA_INEQ_B_GRAD_104_104_104(beliefPenaltyMPC_lub9, beliefPenaltyMPC_sub9, beliefPenaltyMPC_riub9, beliefPenaltyMPC_llb9, beliefPenaltyMPC_slb9, beliefPenaltyMPC_rilb9, beliefPenaltyMPC_lbIdx9, beliefPenaltyMPC_ubIdx9, beliefPenaltyMPC_grad_ineq9, beliefPenaltyMPC_lubbysub9, beliefPenaltyMPC_llbbyslb9);
info->dobj = info->pobj - info->dgap;
info->rdgap = info->pobj ? info->dgap / info->pobj : 1e6;
if( info->rdgap < 0 ) info->rdgap = -info->rdgap;
if( info->mu < beliefPenaltyMPC_SET_ACC_KKTCOMPL
    && (info->rdgap < beliefPenaltyMPC_SET_ACC_RDGAP || info->dgap < beliefPenaltyMPC_SET_ACC_KKTCOMPL)
    && info->res_eq < beliefPenaltyMPC_SET_ACC_RESEQ
    && info->res_ineq < beliefPenaltyMPC_SET_ACC_RESINEQ ){
PRINTTEXT("OPTIMAL (within RESEQ=%2.1e, RESINEQ=%2.1e, (R)DGAP=(%2.1e)%2.1e, MU=%2.1e).\n",beliefPenaltyMPC_SET_ACC_RESEQ, beliefPenaltyMPC_SET_ACC_RESINEQ,beliefPenaltyMPC_SET_ACC_KKTCOMPL,beliefPenaltyMPC_SET_ACC_RDGAP,beliefPenaltyMPC_SET_ACC_KKTCOMPL);
exitcode = beliefPenaltyMPC_OPTIMAL; break; }
if( info->it == beliefPenaltyMPC_SET_MAXIT ){
PRINTTEXT("Maximum number of iterations reached, exiting.\n");
exitcode = beliefPenaltyMPC_MAXITREACHED; break; }
beliefPenaltyMPC_LA_VVADD3_3029(beliefPenaltyMPC_grad_cost, beliefPenaltyMPC_grad_eq, beliefPenaltyMPC_grad_ineq, beliefPenaltyMPC_rd);
beliefPenaltyMPC_LA_DIAG_CHOL_LBUB_325_325_117(params->H1, beliefPenaltyMPC_llbbyslb0, beliefPenaltyMPC_lbIdx0, beliefPenaltyMPC_lubbysub0, beliefPenaltyMPC_ubIdx0, beliefPenaltyMPC_Phi0);
beliefPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_208_325(beliefPenaltyMPC_Phi0, params->C1, beliefPenaltyMPC_V0);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_325(beliefPenaltyMPC_Phi0, beliefPenaltyMPC_rd0, beliefPenaltyMPC_Lbyrd0);
beliefPenaltyMPC_LA_DIAG_CHOL_LBUB_325_325_117(params->H2, beliefPenaltyMPC_llbbyslb1, beliefPenaltyMPC_lbIdx1, beliefPenaltyMPC_lubbysub1, beliefPenaltyMPC_ubIdx1, beliefPenaltyMPC_Phi1);
beliefPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_104_325(beliefPenaltyMPC_Phi1, params->C2, beliefPenaltyMPC_V1);
beliefPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_208_325(beliefPenaltyMPC_Phi1, params->D2, beliefPenaltyMPC_W1);
beliefPenaltyMPC_LA_DENSE_MMTM_208_325_104(beliefPenaltyMPC_W1, beliefPenaltyMPC_V1, beliefPenaltyMPC_Ysd1);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_325(beliefPenaltyMPC_Phi1, beliefPenaltyMPC_rd1, beliefPenaltyMPC_Lbyrd1);
beliefPenaltyMPC_LA_DIAG_CHOL_LBUB_325_325_117(params->H3, beliefPenaltyMPC_llbbyslb2, beliefPenaltyMPC_lbIdx2, beliefPenaltyMPC_lubbysub2, beliefPenaltyMPC_ubIdx2, beliefPenaltyMPC_Phi2);
beliefPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_104_325(beliefPenaltyMPC_Phi2, params->C3, beliefPenaltyMPC_V2);
beliefPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_104_325(beliefPenaltyMPC_Phi2, params->D3, beliefPenaltyMPC_W2);
beliefPenaltyMPC_LA_DENSE_MMTM_104_325_104(beliefPenaltyMPC_W2, beliefPenaltyMPC_V2, beliefPenaltyMPC_Ysd2);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_325(beliefPenaltyMPC_Phi2, beliefPenaltyMPC_rd2, beliefPenaltyMPC_Lbyrd2);
beliefPenaltyMPC_LA_DIAG_CHOL_LBUB_325_325_117(params->H4, beliefPenaltyMPC_llbbyslb3, beliefPenaltyMPC_lbIdx3, beliefPenaltyMPC_lubbysub3, beliefPenaltyMPC_ubIdx3, beliefPenaltyMPC_Phi3);
beliefPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_104_325(beliefPenaltyMPC_Phi3, params->C4, beliefPenaltyMPC_V3);
beliefPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_104_325(beliefPenaltyMPC_Phi3, params->D4, beliefPenaltyMPC_W3);
beliefPenaltyMPC_LA_DENSE_MMTM_104_325_104(beliefPenaltyMPC_W3, beliefPenaltyMPC_V3, beliefPenaltyMPC_Ysd3);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_325(beliefPenaltyMPC_Phi3, beliefPenaltyMPC_rd3, beliefPenaltyMPC_Lbyrd3);
beliefPenaltyMPC_LA_DIAG_CHOL_LBUB_325_325_117(params->H5, beliefPenaltyMPC_llbbyslb4, beliefPenaltyMPC_lbIdx4, beliefPenaltyMPC_lubbysub4, beliefPenaltyMPC_ubIdx4, beliefPenaltyMPC_Phi4);
beliefPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_104_325(beliefPenaltyMPC_Phi4, params->C5, beliefPenaltyMPC_V4);
beliefPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_104_325(beliefPenaltyMPC_Phi4, params->D5, beliefPenaltyMPC_W4);
beliefPenaltyMPC_LA_DENSE_MMTM_104_325_104(beliefPenaltyMPC_W4, beliefPenaltyMPC_V4, beliefPenaltyMPC_Ysd4);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_325(beliefPenaltyMPC_Phi4, beliefPenaltyMPC_rd4, beliefPenaltyMPC_Lbyrd4);
beliefPenaltyMPC_LA_DIAG_CHOL_LBUB_325_325_117(params->H6, beliefPenaltyMPC_llbbyslb5, beliefPenaltyMPC_lbIdx5, beliefPenaltyMPC_lubbysub5, beliefPenaltyMPC_ubIdx5, beliefPenaltyMPC_Phi5);
beliefPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_104_325(beliefPenaltyMPC_Phi5, params->C6, beliefPenaltyMPC_V5);
beliefPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_104_325(beliefPenaltyMPC_Phi5, params->D6, beliefPenaltyMPC_W5);
beliefPenaltyMPC_LA_DENSE_MMTM_104_325_104(beliefPenaltyMPC_W5, beliefPenaltyMPC_V5, beliefPenaltyMPC_Ysd5);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_325(beliefPenaltyMPC_Phi5, beliefPenaltyMPC_rd5, beliefPenaltyMPC_Lbyrd5);
beliefPenaltyMPC_LA_DIAG_CHOL_LBUB_325_325_117(params->H7, beliefPenaltyMPC_llbbyslb6, beliefPenaltyMPC_lbIdx6, beliefPenaltyMPC_lubbysub6, beliefPenaltyMPC_ubIdx6, beliefPenaltyMPC_Phi6);
beliefPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_104_325(beliefPenaltyMPC_Phi6, params->C7, beliefPenaltyMPC_V6);
beliefPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_104_325(beliefPenaltyMPC_Phi6, params->D7, beliefPenaltyMPC_W6);
beliefPenaltyMPC_LA_DENSE_MMTM_104_325_104(beliefPenaltyMPC_W6, beliefPenaltyMPC_V6, beliefPenaltyMPC_Ysd6);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_325(beliefPenaltyMPC_Phi6, beliefPenaltyMPC_rd6, beliefPenaltyMPC_Lbyrd6);
beliefPenaltyMPC_LA_DIAG_CHOL_LBUB_325_325_117(params->H8, beliefPenaltyMPC_llbbyslb7, beliefPenaltyMPC_lbIdx7, beliefPenaltyMPC_lubbysub7, beliefPenaltyMPC_ubIdx7, beliefPenaltyMPC_Phi7);
beliefPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_104_325(beliefPenaltyMPC_Phi7, params->C8, beliefPenaltyMPC_V7);
beliefPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_104_325(beliefPenaltyMPC_Phi7, params->D8, beliefPenaltyMPC_W7);
beliefPenaltyMPC_LA_DENSE_MMTM_104_325_104(beliefPenaltyMPC_W7, beliefPenaltyMPC_V7, beliefPenaltyMPC_Ysd7);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_325(beliefPenaltyMPC_Phi7, beliefPenaltyMPC_rd7, beliefPenaltyMPC_Lbyrd7);
beliefPenaltyMPC_LA_DIAG_CHOL_LBUB_325_325_117(params->H9, beliefPenaltyMPC_llbbyslb8, beliefPenaltyMPC_lbIdx8, beliefPenaltyMPC_lubbysub8, beliefPenaltyMPC_ubIdx8, beliefPenaltyMPC_Phi8);
beliefPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_104_325(beliefPenaltyMPC_Phi8, params->C9, beliefPenaltyMPC_V8);
beliefPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_104_325(beliefPenaltyMPC_Phi8, params->D9, beliefPenaltyMPC_W8);
beliefPenaltyMPC_LA_DENSE_MMTM_104_325_104(beliefPenaltyMPC_W8, beliefPenaltyMPC_V8, beliefPenaltyMPC_Ysd8);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_325(beliefPenaltyMPC_Phi8, beliefPenaltyMPC_rd8, beliefPenaltyMPC_Lbyrd8);
beliefPenaltyMPC_LA_DIAG_CHOL_ONELOOP_LBUB_104_104_104(params->H10, beliefPenaltyMPC_llbbyslb9, beliefPenaltyMPC_lbIdx9, beliefPenaltyMPC_lubbysub9, beliefPenaltyMPC_ubIdx9, beliefPenaltyMPC_Phi9);
beliefPenaltyMPC_LA_DIAG_MATRIXFORWARDSUB_104_104(beliefPenaltyMPC_Phi9, params->D10, beliefPenaltyMPC_W9);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_104(beliefPenaltyMPC_Phi9, beliefPenaltyMPC_rd9, beliefPenaltyMPC_Lbyrd9);
beliefPenaltyMPC_LA_DENSE_MMT2_208_325_325(beliefPenaltyMPC_V0, beliefPenaltyMPC_W1, beliefPenaltyMPC_Yd0);
beliefPenaltyMPC_LA_DENSE_MVMSUB2_208_325_325(beliefPenaltyMPC_V0, beliefPenaltyMPC_Lbyrd0, beliefPenaltyMPC_W1, beliefPenaltyMPC_Lbyrd1, beliefPenaltyMPC_re0, beliefPenaltyMPC_beta0);
beliefPenaltyMPC_LA_DENSE_MMT2_104_325_325(beliefPenaltyMPC_V1, beliefPenaltyMPC_W2, beliefPenaltyMPC_Yd1);
beliefPenaltyMPC_LA_DENSE_MVMSUB2_104_325_325(beliefPenaltyMPC_V1, beliefPenaltyMPC_Lbyrd1, beliefPenaltyMPC_W2, beliefPenaltyMPC_Lbyrd2, beliefPenaltyMPC_re1, beliefPenaltyMPC_beta1);
beliefPenaltyMPC_LA_DENSE_MMT2_104_325_325(beliefPenaltyMPC_V2, beliefPenaltyMPC_W3, beliefPenaltyMPC_Yd2);
beliefPenaltyMPC_LA_DENSE_MVMSUB2_104_325_325(beliefPenaltyMPC_V2, beliefPenaltyMPC_Lbyrd2, beliefPenaltyMPC_W3, beliefPenaltyMPC_Lbyrd3, beliefPenaltyMPC_re2, beliefPenaltyMPC_beta2);
beliefPenaltyMPC_LA_DENSE_MMT2_104_325_325(beliefPenaltyMPC_V3, beliefPenaltyMPC_W4, beliefPenaltyMPC_Yd3);
beliefPenaltyMPC_LA_DENSE_MVMSUB2_104_325_325(beliefPenaltyMPC_V3, beliefPenaltyMPC_Lbyrd3, beliefPenaltyMPC_W4, beliefPenaltyMPC_Lbyrd4, beliefPenaltyMPC_re3, beliefPenaltyMPC_beta3);
beliefPenaltyMPC_LA_DENSE_MMT2_104_325_325(beliefPenaltyMPC_V4, beliefPenaltyMPC_W5, beliefPenaltyMPC_Yd4);
beliefPenaltyMPC_LA_DENSE_MVMSUB2_104_325_325(beliefPenaltyMPC_V4, beliefPenaltyMPC_Lbyrd4, beliefPenaltyMPC_W5, beliefPenaltyMPC_Lbyrd5, beliefPenaltyMPC_re4, beliefPenaltyMPC_beta4);
beliefPenaltyMPC_LA_DENSE_MMT2_104_325_325(beliefPenaltyMPC_V5, beliefPenaltyMPC_W6, beliefPenaltyMPC_Yd5);
beliefPenaltyMPC_LA_DENSE_MVMSUB2_104_325_325(beliefPenaltyMPC_V5, beliefPenaltyMPC_Lbyrd5, beliefPenaltyMPC_W6, beliefPenaltyMPC_Lbyrd6, beliefPenaltyMPC_re5, beliefPenaltyMPC_beta5);
beliefPenaltyMPC_LA_DENSE_MMT2_104_325_325(beliefPenaltyMPC_V6, beliefPenaltyMPC_W7, beliefPenaltyMPC_Yd6);
beliefPenaltyMPC_LA_DENSE_MVMSUB2_104_325_325(beliefPenaltyMPC_V6, beliefPenaltyMPC_Lbyrd6, beliefPenaltyMPC_W7, beliefPenaltyMPC_Lbyrd7, beliefPenaltyMPC_re6, beliefPenaltyMPC_beta6);
beliefPenaltyMPC_LA_DENSE_MMT2_104_325_325(beliefPenaltyMPC_V7, beliefPenaltyMPC_W8, beliefPenaltyMPC_Yd7);
beliefPenaltyMPC_LA_DENSE_MVMSUB2_104_325_325(beliefPenaltyMPC_V7, beliefPenaltyMPC_Lbyrd7, beliefPenaltyMPC_W8, beliefPenaltyMPC_Lbyrd8, beliefPenaltyMPC_re7, beliefPenaltyMPC_beta7);
beliefPenaltyMPC_LA_DENSE_MMT2_104_325_104(beliefPenaltyMPC_V8, beliefPenaltyMPC_W9, beliefPenaltyMPC_Yd8);
beliefPenaltyMPC_LA_DENSE_MVMSUB2_104_325_104(beliefPenaltyMPC_V8, beliefPenaltyMPC_Lbyrd8, beliefPenaltyMPC_W9, beliefPenaltyMPC_Lbyrd9, beliefPenaltyMPC_re8, beliefPenaltyMPC_beta8);
beliefPenaltyMPC_LA_DENSE_CHOL_208(beliefPenaltyMPC_Yd0, beliefPenaltyMPC_Ld0);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_208(beliefPenaltyMPC_Ld0, beliefPenaltyMPC_beta0, beliefPenaltyMPC_yy0);
beliefPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_104_208(beliefPenaltyMPC_Ld0, beliefPenaltyMPC_Ysd1, beliefPenaltyMPC_Lsd1);
beliefPenaltyMPC_LA_DENSE_MMTSUB_104_208(beliefPenaltyMPC_Lsd1, beliefPenaltyMPC_Yd1);
beliefPenaltyMPC_LA_DENSE_CHOL_104(beliefPenaltyMPC_Yd1, beliefPenaltyMPC_Ld1);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_104_208(beliefPenaltyMPC_Lsd1, beliefPenaltyMPC_yy0, beliefPenaltyMPC_beta1, beliefPenaltyMPC_bmy1);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_104(beliefPenaltyMPC_Ld1, beliefPenaltyMPC_bmy1, beliefPenaltyMPC_yy1);
beliefPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_104_104(beliefPenaltyMPC_Ld1, beliefPenaltyMPC_Ysd2, beliefPenaltyMPC_Lsd2);
beliefPenaltyMPC_LA_DENSE_MMTSUB_104_104(beliefPenaltyMPC_Lsd2, beliefPenaltyMPC_Yd2);
beliefPenaltyMPC_LA_DENSE_CHOL_104(beliefPenaltyMPC_Yd2, beliefPenaltyMPC_Ld2);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_104_104(beliefPenaltyMPC_Lsd2, beliefPenaltyMPC_yy1, beliefPenaltyMPC_beta2, beliefPenaltyMPC_bmy2);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_104(beliefPenaltyMPC_Ld2, beliefPenaltyMPC_bmy2, beliefPenaltyMPC_yy2);
beliefPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_104_104(beliefPenaltyMPC_Ld2, beliefPenaltyMPC_Ysd3, beliefPenaltyMPC_Lsd3);
beliefPenaltyMPC_LA_DENSE_MMTSUB_104_104(beliefPenaltyMPC_Lsd3, beliefPenaltyMPC_Yd3);
beliefPenaltyMPC_LA_DENSE_CHOL_104(beliefPenaltyMPC_Yd3, beliefPenaltyMPC_Ld3);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_104_104(beliefPenaltyMPC_Lsd3, beliefPenaltyMPC_yy2, beliefPenaltyMPC_beta3, beliefPenaltyMPC_bmy3);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_104(beliefPenaltyMPC_Ld3, beliefPenaltyMPC_bmy3, beliefPenaltyMPC_yy3);
beliefPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_104_104(beliefPenaltyMPC_Ld3, beliefPenaltyMPC_Ysd4, beliefPenaltyMPC_Lsd4);
beliefPenaltyMPC_LA_DENSE_MMTSUB_104_104(beliefPenaltyMPC_Lsd4, beliefPenaltyMPC_Yd4);
beliefPenaltyMPC_LA_DENSE_CHOL_104(beliefPenaltyMPC_Yd4, beliefPenaltyMPC_Ld4);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_104_104(beliefPenaltyMPC_Lsd4, beliefPenaltyMPC_yy3, beliefPenaltyMPC_beta4, beliefPenaltyMPC_bmy4);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_104(beliefPenaltyMPC_Ld4, beliefPenaltyMPC_bmy4, beliefPenaltyMPC_yy4);
beliefPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_104_104(beliefPenaltyMPC_Ld4, beliefPenaltyMPC_Ysd5, beliefPenaltyMPC_Lsd5);
beliefPenaltyMPC_LA_DENSE_MMTSUB_104_104(beliefPenaltyMPC_Lsd5, beliefPenaltyMPC_Yd5);
beliefPenaltyMPC_LA_DENSE_CHOL_104(beliefPenaltyMPC_Yd5, beliefPenaltyMPC_Ld5);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_104_104(beliefPenaltyMPC_Lsd5, beliefPenaltyMPC_yy4, beliefPenaltyMPC_beta5, beliefPenaltyMPC_bmy5);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_104(beliefPenaltyMPC_Ld5, beliefPenaltyMPC_bmy5, beliefPenaltyMPC_yy5);
beliefPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_104_104(beliefPenaltyMPC_Ld5, beliefPenaltyMPC_Ysd6, beliefPenaltyMPC_Lsd6);
beliefPenaltyMPC_LA_DENSE_MMTSUB_104_104(beliefPenaltyMPC_Lsd6, beliefPenaltyMPC_Yd6);
beliefPenaltyMPC_LA_DENSE_CHOL_104(beliefPenaltyMPC_Yd6, beliefPenaltyMPC_Ld6);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_104_104(beliefPenaltyMPC_Lsd6, beliefPenaltyMPC_yy5, beliefPenaltyMPC_beta6, beliefPenaltyMPC_bmy6);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_104(beliefPenaltyMPC_Ld6, beliefPenaltyMPC_bmy6, beliefPenaltyMPC_yy6);
beliefPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_104_104(beliefPenaltyMPC_Ld6, beliefPenaltyMPC_Ysd7, beliefPenaltyMPC_Lsd7);
beliefPenaltyMPC_LA_DENSE_MMTSUB_104_104(beliefPenaltyMPC_Lsd7, beliefPenaltyMPC_Yd7);
beliefPenaltyMPC_LA_DENSE_CHOL_104(beliefPenaltyMPC_Yd7, beliefPenaltyMPC_Ld7);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_104_104(beliefPenaltyMPC_Lsd7, beliefPenaltyMPC_yy6, beliefPenaltyMPC_beta7, beliefPenaltyMPC_bmy7);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_104(beliefPenaltyMPC_Ld7, beliefPenaltyMPC_bmy7, beliefPenaltyMPC_yy7);
beliefPenaltyMPC_LA_DENSE_MATRIXTFORWARDSUB_104_104(beliefPenaltyMPC_Ld7, beliefPenaltyMPC_Ysd8, beliefPenaltyMPC_Lsd8);
beliefPenaltyMPC_LA_DENSE_MMTSUB_104_104(beliefPenaltyMPC_Lsd8, beliefPenaltyMPC_Yd8);
beliefPenaltyMPC_LA_DENSE_CHOL_104(beliefPenaltyMPC_Yd8, beliefPenaltyMPC_Ld8);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_104_104(beliefPenaltyMPC_Lsd8, beliefPenaltyMPC_yy7, beliefPenaltyMPC_beta8, beliefPenaltyMPC_bmy8);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_104(beliefPenaltyMPC_Ld8, beliefPenaltyMPC_bmy8, beliefPenaltyMPC_yy8);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_104(beliefPenaltyMPC_Ld8, beliefPenaltyMPC_yy8, beliefPenaltyMPC_dvaff8);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_104_104(beliefPenaltyMPC_Lsd8, beliefPenaltyMPC_dvaff8, beliefPenaltyMPC_yy7, beliefPenaltyMPC_bmy7);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_104(beliefPenaltyMPC_Ld7, beliefPenaltyMPC_bmy7, beliefPenaltyMPC_dvaff7);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_104_104(beliefPenaltyMPC_Lsd7, beliefPenaltyMPC_dvaff7, beliefPenaltyMPC_yy6, beliefPenaltyMPC_bmy6);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_104(beliefPenaltyMPC_Ld6, beliefPenaltyMPC_bmy6, beliefPenaltyMPC_dvaff6);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_104_104(beliefPenaltyMPC_Lsd6, beliefPenaltyMPC_dvaff6, beliefPenaltyMPC_yy5, beliefPenaltyMPC_bmy5);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_104(beliefPenaltyMPC_Ld5, beliefPenaltyMPC_bmy5, beliefPenaltyMPC_dvaff5);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_104_104(beliefPenaltyMPC_Lsd5, beliefPenaltyMPC_dvaff5, beliefPenaltyMPC_yy4, beliefPenaltyMPC_bmy4);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_104(beliefPenaltyMPC_Ld4, beliefPenaltyMPC_bmy4, beliefPenaltyMPC_dvaff4);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_104_104(beliefPenaltyMPC_Lsd4, beliefPenaltyMPC_dvaff4, beliefPenaltyMPC_yy3, beliefPenaltyMPC_bmy3);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_104(beliefPenaltyMPC_Ld3, beliefPenaltyMPC_bmy3, beliefPenaltyMPC_dvaff3);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_104_104(beliefPenaltyMPC_Lsd3, beliefPenaltyMPC_dvaff3, beliefPenaltyMPC_yy2, beliefPenaltyMPC_bmy2);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_104(beliefPenaltyMPC_Ld2, beliefPenaltyMPC_bmy2, beliefPenaltyMPC_dvaff2);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_104_104(beliefPenaltyMPC_Lsd2, beliefPenaltyMPC_dvaff2, beliefPenaltyMPC_yy1, beliefPenaltyMPC_bmy1);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_104(beliefPenaltyMPC_Ld1, beliefPenaltyMPC_bmy1, beliefPenaltyMPC_dvaff1);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_104_208(beliefPenaltyMPC_Lsd1, beliefPenaltyMPC_dvaff1, beliefPenaltyMPC_yy0, beliefPenaltyMPC_bmy0);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_208(beliefPenaltyMPC_Ld0, beliefPenaltyMPC_bmy0, beliefPenaltyMPC_dvaff0);
beliefPenaltyMPC_LA_DENSE_MTVM_208_325(params->C1, beliefPenaltyMPC_dvaff0, beliefPenaltyMPC_grad_eq0);
beliefPenaltyMPC_LA_DENSE_MTVM2_104_325_208(params->C2, beliefPenaltyMPC_dvaff1, params->D2, beliefPenaltyMPC_dvaff0, beliefPenaltyMPC_grad_eq1);
beliefPenaltyMPC_LA_DENSE_MTVM2_104_325_104(params->C3, beliefPenaltyMPC_dvaff2, params->D3, beliefPenaltyMPC_dvaff1, beliefPenaltyMPC_grad_eq2);
beliefPenaltyMPC_LA_DENSE_MTVM2_104_325_104(params->C4, beliefPenaltyMPC_dvaff3, params->D4, beliefPenaltyMPC_dvaff2, beliefPenaltyMPC_grad_eq3);
beliefPenaltyMPC_LA_DENSE_MTVM2_104_325_104(params->C5, beliefPenaltyMPC_dvaff4, params->D5, beliefPenaltyMPC_dvaff3, beliefPenaltyMPC_grad_eq4);
beliefPenaltyMPC_LA_DENSE_MTVM2_104_325_104(params->C6, beliefPenaltyMPC_dvaff5, params->D6, beliefPenaltyMPC_dvaff4, beliefPenaltyMPC_grad_eq5);
beliefPenaltyMPC_LA_DENSE_MTVM2_104_325_104(params->C7, beliefPenaltyMPC_dvaff6, params->D7, beliefPenaltyMPC_dvaff5, beliefPenaltyMPC_grad_eq6);
beliefPenaltyMPC_LA_DENSE_MTVM2_104_325_104(params->C8, beliefPenaltyMPC_dvaff7, params->D8, beliefPenaltyMPC_dvaff6, beliefPenaltyMPC_grad_eq7);
beliefPenaltyMPC_LA_DENSE_MTVM2_104_325_104(params->C9, beliefPenaltyMPC_dvaff8, params->D9, beliefPenaltyMPC_dvaff7, beliefPenaltyMPC_grad_eq8);
beliefPenaltyMPC_LA_DENSE_MTVM_104_104(params->D10, beliefPenaltyMPC_dvaff8, beliefPenaltyMPC_grad_eq9);
beliefPenaltyMPC_LA_VSUB2_3029(beliefPenaltyMPC_rd, beliefPenaltyMPC_grad_eq, beliefPenaltyMPC_rd);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_325(beliefPenaltyMPC_Phi0, beliefPenaltyMPC_rd0, beliefPenaltyMPC_dzaff0);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_325(beliefPenaltyMPC_Phi1, beliefPenaltyMPC_rd1, beliefPenaltyMPC_dzaff1);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_325(beliefPenaltyMPC_Phi2, beliefPenaltyMPC_rd2, beliefPenaltyMPC_dzaff2);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_325(beliefPenaltyMPC_Phi3, beliefPenaltyMPC_rd3, beliefPenaltyMPC_dzaff3);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_325(beliefPenaltyMPC_Phi4, beliefPenaltyMPC_rd4, beliefPenaltyMPC_dzaff4);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_325(beliefPenaltyMPC_Phi5, beliefPenaltyMPC_rd5, beliefPenaltyMPC_dzaff5);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_325(beliefPenaltyMPC_Phi6, beliefPenaltyMPC_rd6, beliefPenaltyMPC_dzaff6);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_325(beliefPenaltyMPC_Phi7, beliefPenaltyMPC_rd7, beliefPenaltyMPC_dzaff7);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_325(beliefPenaltyMPC_Phi8, beliefPenaltyMPC_rd8, beliefPenaltyMPC_dzaff8);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_104(beliefPenaltyMPC_Phi9, beliefPenaltyMPC_rd9, beliefPenaltyMPC_dzaff9);
beliefPenaltyMPC_LA_VSUB_INDEXED_325(beliefPenaltyMPC_dzaff0, beliefPenaltyMPC_lbIdx0, beliefPenaltyMPC_rilb0, beliefPenaltyMPC_dslbaff0);
beliefPenaltyMPC_LA_VSUB3_325(beliefPenaltyMPC_llbbyslb0, beliefPenaltyMPC_dslbaff0, beliefPenaltyMPC_llb0, beliefPenaltyMPC_dllbaff0);
beliefPenaltyMPC_LA_VSUB2_INDEXED_117(beliefPenaltyMPC_riub0, beliefPenaltyMPC_dzaff0, beliefPenaltyMPC_ubIdx0, beliefPenaltyMPC_dsubaff0);
beliefPenaltyMPC_LA_VSUB3_117(beliefPenaltyMPC_lubbysub0, beliefPenaltyMPC_dsubaff0, beliefPenaltyMPC_lub0, beliefPenaltyMPC_dlubaff0);
beliefPenaltyMPC_LA_VSUB_INDEXED_325(beliefPenaltyMPC_dzaff1, beliefPenaltyMPC_lbIdx1, beliefPenaltyMPC_rilb1, beliefPenaltyMPC_dslbaff1);
beliefPenaltyMPC_LA_VSUB3_325(beliefPenaltyMPC_llbbyslb1, beliefPenaltyMPC_dslbaff1, beliefPenaltyMPC_llb1, beliefPenaltyMPC_dllbaff1);
beliefPenaltyMPC_LA_VSUB2_INDEXED_117(beliefPenaltyMPC_riub1, beliefPenaltyMPC_dzaff1, beliefPenaltyMPC_ubIdx1, beliefPenaltyMPC_dsubaff1);
beliefPenaltyMPC_LA_VSUB3_117(beliefPenaltyMPC_lubbysub1, beliefPenaltyMPC_dsubaff1, beliefPenaltyMPC_lub1, beliefPenaltyMPC_dlubaff1);
beliefPenaltyMPC_LA_VSUB_INDEXED_325(beliefPenaltyMPC_dzaff2, beliefPenaltyMPC_lbIdx2, beliefPenaltyMPC_rilb2, beliefPenaltyMPC_dslbaff2);
beliefPenaltyMPC_LA_VSUB3_325(beliefPenaltyMPC_llbbyslb2, beliefPenaltyMPC_dslbaff2, beliefPenaltyMPC_llb2, beliefPenaltyMPC_dllbaff2);
beliefPenaltyMPC_LA_VSUB2_INDEXED_117(beliefPenaltyMPC_riub2, beliefPenaltyMPC_dzaff2, beliefPenaltyMPC_ubIdx2, beliefPenaltyMPC_dsubaff2);
beliefPenaltyMPC_LA_VSUB3_117(beliefPenaltyMPC_lubbysub2, beliefPenaltyMPC_dsubaff2, beliefPenaltyMPC_lub2, beliefPenaltyMPC_dlubaff2);
beliefPenaltyMPC_LA_VSUB_INDEXED_325(beliefPenaltyMPC_dzaff3, beliefPenaltyMPC_lbIdx3, beliefPenaltyMPC_rilb3, beliefPenaltyMPC_dslbaff3);
beliefPenaltyMPC_LA_VSUB3_325(beliefPenaltyMPC_llbbyslb3, beliefPenaltyMPC_dslbaff3, beliefPenaltyMPC_llb3, beliefPenaltyMPC_dllbaff3);
beliefPenaltyMPC_LA_VSUB2_INDEXED_117(beliefPenaltyMPC_riub3, beliefPenaltyMPC_dzaff3, beliefPenaltyMPC_ubIdx3, beliefPenaltyMPC_dsubaff3);
beliefPenaltyMPC_LA_VSUB3_117(beliefPenaltyMPC_lubbysub3, beliefPenaltyMPC_dsubaff3, beliefPenaltyMPC_lub3, beliefPenaltyMPC_dlubaff3);
beliefPenaltyMPC_LA_VSUB_INDEXED_325(beliefPenaltyMPC_dzaff4, beliefPenaltyMPC_lbIdx4, beliefPenaltyMPC_rilb4, beliefPenaltyMPC_dslbaff4);
beliefPenaltyMPC_LA_VSUB3_325(beliefPenaltyMPC_llbbyslb4, beliefPenaltyMPC_dslbaff4, beliefPenaltyMPC_llb4, beliefPenaltyMPC_dllbaff4);
beliefPenaltyMPC_LA_VSUB2_INDEXED_117(beliefPenaltyMPC_riub4, beliefPenaltyMPC_dzaff4, beliefPenaltyMPC_ubIdx4, beliefPenaltyMPC_dsubaff4);
beliefPenaltyMPC_LA_VSUB3_117(beliefPenaltyMPC_lubbysub4, beliefPenaltyMPC_dsubaff4, beliefPenaltyMPC_lub4, beliefPenaltyMPC_dlubaff4);
beliefPenaltyMPC_LA_VSUB_INDEXED_325(beliefPenaltyMPC_dzaff5, beliefPenaltyMPC_lbIdx5, beliefPenaltyMPC_rilb5, beliefPenaltyMPC_dslbaff5);
beliefPenaltyMPC_LA_VSUB3_325(beliefPenaltyMPC_llbbyslb5, beliefPenaltyMPC_dslbaff5, beliefPenaltyMPC_llb5, beliefPenaltyMPC_dllbaff5);
beliefPenaltyMPC_LA_VSUB2_INDEXED_117(beliefPenaltyMPC_riub5, beliefPenaltyMPC_dzaff5, beliefPenaltyMPC_ubIdx5, beliefPenaltyMPC_dsubaff5);
beliefPenaltyMPC_LA_VSUB3_117(beliefPenaltyMPC_lubbysub5, beliefPenaltyMPC_dsubaff5, beliefPenaltyMPC_lub5, beliefPenaltyMPC_dlubaff5);
beliefPenaltyMPC_LA_VSUB_INDEXED_325(beliefPenaltyMPC_dzaff6, beliefPenaltyMPC_lbIdx6, beliefPenaltyMPC_rilb6, beliefPenaltyMPC_dslbaff6);
beliefPenaltyMPC_LA_VSUB3_325(beliefPenaltyMPC_llbbyslb6, beliefPenaltyMPC_dslbaff6, beliefPenaltyMPC_llb6, beliefPenaltyMPC_dllbaff6);
beliefPenaltyMPC_LA_VSUB2_INDEXED_117(beliefPenaltyMPC_riub6, beliefPenaltyMPC_dzaff6, beliefPenaltyMPC_ubIdx6, beliefPenaltyMPC_dsubaff6);
beliefPenaltyMPC_LA_VSUB3_117(beliefPenaltyMPC_lubbysub6, beliefPenaltyMPC_dsubaff6, beliefPenaltyMPC_lub6, beliefPenaltyMPC_dlubaff6);
beliefPenaltyMPC_LA_VSUB_INDEXED_325(beliefPenaltyMPC_dzaff7, beliefPenaltyMPC_lbIdx7, beliefPenaltyMPC_rilb7, beliefPenaltyMPC_dslbaff7);
beliefPenaltyMPC_LA_VSUB3_325(beliefPenaltyMPC_llbbyslb7, beliefPenaltyMPC_dslbaff7, beliefPenaltyMPC_llb7, beliefPenaltyMPC_dllbaff7);
beliefPenaltyMPC_LA_VSUB2_INDEXED_117(beliefPenaltyMPC_riub7, beliefPenaltyMPC_dzaff7, beliefPenaltyMPC_ubIdx7, beliefPenaltyMPC_dsubaff7);
beliefPenaltyMPC_LA_VSUB3_117(beliefPenaltyMPC_lubbysub7, beliefPenaltyMPC_dsubaff7, beliefPenaltyMPC_lub7, beliefPenaltyMPC_dlubaff7);
beliefPenaltyMPC_LA_VSUB_INDEXED_325(beliefPenaltyMPC_dzaff8, beliefPenaltyMPC_lbIdx8, beliefPenaltyMPC_rilb8, beliefPenaltyMPC_dslbaff8);
beliefPenaltyMPC_LA_VSUB3_325(beliefPenaltyMPC_llbbyslb8, beliefPenaltyMPC_dslbaff8, beliefPenaltyMPC_llb8, beliefPenaltyMPC_dllbaff8);
beliefPenaltyMPC_LA_VSUB2_INDEXED_117(beliefPenaltyMPC_riub8, beliefPenaltyMPC_dzaff8, beliefPenaltyMPC_ubIdx8, beliefPenaltyMPC_dsubaff8);
beliefPenaltyMPC_LA_VSUB3_117(beliefPenaltyMPC_lubbysub8, beliefPenaltyMPC_dsubaff8, beliefPenaltyMPC_lub8, beliefPenaltyMPC_dlubaff8);
beliefPenaltyMPC_LA_VSUB_INDEXED_104(beliefPenaltyMPC_dzaff9, beliefPenaltyMPC_lbIdx9, beliefPenaltyMPC_rilb9, beliefPenaltyMPC_dslbaff9);
beliefPenaltyMPC_LA_VSUB3_104(beliefPenaltyMPC_llbbyslb9, beliefPenaltyMPC_dslbaff9, beliefPenaltyMPC_llb9, beliefPenaltyMPC_dllbaff9);
beliefPenaltyMPC_LA_VSUB2_INDEXED_104(beliefPenaltyMPC_riub9, beliefPenaltyMPC_dzaff9, beliefPenaltyMPC_ubIdx9, beliefPenaltyMPC_dsubaff9);
beliefPenaltyMPC_LA_VSUB3_104(beliefPenaltyMPC_lubbysub9, beliefPenaltyMPC_dsubaff9, beliefPenaltyMPC_lub9, beliefPenaltyMPC_dlubaff9);
info->lsit_aff = beliefPenaltyMPC_LINESEARCH_BACKTRACKING_AFFINE(beliefPenaltyMPC_l, beliefPenaltyMPC_s, beliefPenaltyMPC_dl_aff, beliefPenaltyMPC_ds_aff, &info->step_aff, &info->mu_aff);
if( info->lsit_aff == beliefPenaltyMPC_NOPROGRESS ){
PRINTTEXT("Affine line search could not proceed at iteration %d.\nThe problem might be infeasible -- exiting.\n",info->it+1);
exitcode = beliefPenaltyMPC_NOPROGRESS; break;
}
sigma_3rdroot = info->mu_aff / info->mu;
info->sigma = sigma_3rdroot*sigma_3rdroot*sigma_3rdroot;
musigma = info->mu * info->sigma;
beliefPenaltyMPC_LA_VSUB5_4186(beliefPenaltyMPC_ds_aff, beliefPenaltyMPC_dl_aff, musigma, beliefPenaltyMPC_ccrhs);
beliefPenaltyMPC_LA_VSUB6_INDEXED_325_117_325(beliefPenaltyMPC_ccrhsub0, beliefPenaltyMPC_sub0, beliefPenaltyMPC_ubIdx0, beliefPenaltyMPC_ccrhsl0, beliefPenaltyMPC_slb0, beliefPenaltyMPC_lbIdx0, beliefPenaltyMPC_rd0);
beliefPenaltyMPC_LA_VSUB6_INDEXED_325_117_325(beliefPenaltyMPC_ccrhsub1, beliefPenaltyMPC_sub1, beliefPenaltyMPC_ubIdx1, beliefPenaltyMPC_ccrhsl1, beliefPenaltyMPC_slb1, beliefPenaltyMPC_lbIdx1, beliefPenaltyMPC_rd1);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_325(beliefPenaltyMPC_Phi0, beliefPenaltyMPC_rd0, beliefPenaltyMPC_Lbyrd0);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_325(beliefPenaltyMPC_Phi1, beliefPenaltyMPC_rd1, beliefPenaltyMPC_Lbyrd1);
beliefPenaltyMPC_LA_DENSE_2MVMADD_208_325_325(beliefPenaltyMPC_V0, beliefPenaltyMPC_Lbyrd0, beliefPenaltyMPC_W1, beliefPenaltyMPC_Lbyrd1, beliefPenaltyMPC_beta0);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_208(beliefPenaltyMPC_Ld0, beliefPenaltyMPC_beta0, beliefPenaltyMPC_yy0);
beliefPenaltyMPC_LA_VSUB6_INDEXED_325_117_325(beliefPenaltyMPC_ccrhsub2, beliefPenaltyMPC_sub2, beliefPenaltyMPC_ubIdx2, beliefPenaltyMPC_ccrhsl2, beliefPenaltyMPC_slb2, beliefPenaltyMPC_lbIdx2, beliefPenaltyMPC_rd2);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_325(beliefPenaltyMPC_Phi2, beliefPenaltyMPC_rd2, beliefPenaltyMPC_Lbyrd2);
beliefPenaltyMPC_LA_DENSE_2MVMADD_104_325_325(beliefPenaltyMPC_V1, beliefPenaltyMPC_Lbyrd1, beliefPenaltyMPC_W2, beliefPenaltyMPC_Lbyrd2, beliefPenaltyMPC_beta1);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_104_208(beliefPenaltyMPC_Lsd1, beliefPenaltyMPC_yy0, beliefPenaltyMPC_beta1, beliefPenaltyMPC_bmy1);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_104(beliefPenaltyMPC_Ld1, beliefPenaltyMPC_bmy1, beliefPenaltyMPC_yy1);
beliefPenaltyMPC_LA_VSUB6_INDEXED_325_117_325(beliefPenaltyMPC_ccrhsub3, beliefPenaltyMPC_sub3, beliefPenaltyMPC_ubIdx3, beliefPenaltyMPC_ccrhsl3, beliefPenaltyMPC_slb3, beliefPenaltyMPC_lbIdx3, beliefPenaltyMPC_rd3);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_325(beliefPenaltyMPC_Phi3, beliefPenaltyMPC_rd3, beliefPenaltyMPC_Lbyrd3);
beliefPenaltyMPC_LA_DENSE_2MVMADD_104_325_325(beliefPenaltyMPC_V2, beliefPenaltyMPC_Lbyrd2, beliefPenaltyMPC_W3, beliefPenaltyMPC_Lbyrd3, beliefPenaltyMPC_beta2);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_104_104(beliefPenaltyMPC_Lsd2, beliefPenaltyMPC_yy1, beliefPenaltyMPC_beta2, beliefPenaltyMPC_bmy2);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_104(beliefPenaltyMPC_Ld2, beliefPenaltyMPC_bmy2, beliefPenaltyMPC_yy2);
beliefPenaltyMPC_LA_VSUB6_INDEXED_325_117_325(beliefPenaltyMPC_ccrhsub4, beliefPenaltyMPC_sub4, beliefPenaltyMPC_ubIdx4, beliefPenaltyMPC_ccrhsl4, beliefPenaltyMPC_slb4, beliefPenaltyMPC_lbIdx4, beliefPenaltyMPC_rd4);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_325(beliefPenaltyMPC_Phi4, beliefPenaltyMPC_rd4, beliefPenaltyMPC_Lbyrd4);
beliefPenaltyMPC_LA_DENSE_2MVMADD_104_325_325(beliefPenaltyMPC_V3, beliefPenaltyMPC_Lbyrd3, beliefPenaltyMPC_W4, beliefPenaltyMPC_Lbyrd4, beliefPenaltyMPC_beta3);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_104_104(beliefPenaltyMPC_Lsd3, beliefPenaltyMPC_yy2, beliefPenaltyMPC_beta3, beliefPenaltyMPC_bmy3);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_104(beliefPenaltyMPC_Ld3, beliefPenaltyMPC_bmy3, beliefPenaltyMPC_yy3);
beliefPenaltyMPC_LA_VSUB6_INDEXED_325_117_325(beliefPenaltyMPC_ccrhsub5, beliefPenaltyMPC_sub5, beliefPenaltyMPC_ubIdx5, beliefPenaltyMPC_ccrhsl5, beliefPenaltyMPC_slb5, beliefPenaltyMPC_lbIdx5, beliefPenaltyMPC_rd5);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_325(beliefPenaltyMPC_Phi5, beliefPenaltyMPC_rd5, beliefPenaltyMPC_Lbyrd5);
beliefPenaltyMPC_LA_DENSE_2MVMADD_104_325_325(beliefPenaltyMPC_V4, beliefPenaltyMPC_Lbyrd4, beliefPenaltyMPC_W5, beliefPenaltyMPC_Lbyrd5, beliefPenaltyMPC_beta4);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_104_104(beliefPenaltyMPC_Lsd4, beliefPenaltyMPC_yy3, beliefPenaltyMPC_beta4, beliefPenaltyMPC_bmy4);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_104(beliefPenaltyMPC_Ld4, beliefPenaltyMPC_bmy4, beliefPenaltyMPC_yy4);
beliefPenaltyMPC_LA_VSUB6_INDEXED_325_117_325(beliefPenaltyMPC_ccrhsub6, beliefPenaltyMPC_sub6, beliefPenaltyMPC_ubIdx6, beliefPenaltyMPC_ccrhsl6, beliefPenaltyMPC_slb6, beliefPenaltyMPC_lbIdx6, beliefPenaltyMPC_rd6);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_325(beliefPenaltyMPC_Phi6, beliefPenaltyMPC_rd6, beliefPenaltyMPC_Lbyrd6);
beliefPenaltyMPC_LA_DENSE_2MVMADD_104_325_325(beliefPenaltyMPC_V5, beliefPenaltyMPC_Lbyrd5, beliefPenaltyMPC_W6, beliefPenaltyMPC_Lbyrd6, beliefPenaltyMPC_beta5);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_104_104(beliefPenaltyMPC_Lsd5, beliefPenaltyMPC_yy4, beliefPenaltyMPC_beta5, beliefPenaltyMPC_bmy5);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_104(beliefPenaltyMPC_Ld5, beliefPenaltyMPC_bmy5, beliefPenaltyMPC_yy5);
beliefPenaltyMPC_LA_VSUB6_INDEXED_325_117_325(beliefPenaltyMPC_ccrhsub7, beliefPenaltyMPC_sub7, beliefPenaltyMPC_ubIdx7, beliefPenaltyMPC_ccrhsl7, beliefPenaltyMPC_slb7, beliefPenaltyMPC_lbIdx7, beliefPenaltyMPC_rd7);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_325(beliefPenaltyMPC_Phi7, beliefPenaltyMPC_rd7, beliefPenaltyMPC_Lbyrd7);
beliefPenaltyMPC_LA_DENSE_2MVMADD_104_325_325(beliefPenaltyMPC_V6, beliefPenaltyMPC_Lbyrd6, beliefPenaltyMPC_W7, beliefPenaltyMPC_Lbyrd7, beliefPenaltyMPC_beta6);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_104_104(beliefPenaltyMPC_Lsd6, beliefPenaltyMPC_yy5, beliefPenaltyMPC_beta6, beliefPenaltyMPC_bmy6);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_104(beliefPenaltyMPC_Ld6, beliefPenaltyMPC_bmy6, beliefPenaltyMPC_yy6);
beliefPenaltyMPC_LA_VSUB6_INDEXED_325_117_325(beliefPenaltyMPC_ccrhsub8, beliefPenaltyMPC_sub8, beliefPenaltyMPC_ubIdx8, beliefPenaltyMPC_ccrhsl8, beliefPenaltyMPC_slb8, beliefPenaltyMPC_lbIdx8, beliefPenaltyMPC_rd8);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_325(beliefPenaltyMPC_Phi8, beliefPenaltyMPC_rd8, beliefPenaltyMPC_Lbyrd8);
beliefPenaltyMPC_LA_DENSE_2MVMADD_104_325_325(beliefPenaltyMPC_V7, beliefPenaltyMPC_Lbyrd7, beliefPenaltyMPC_W8, beliefPenaltyMPC_Lbyrd8, beliefPenaltyMPC_beta7);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_104_104(beliefPenaltyMPC_Lsd7, beliefPenaltyMPC_yy6, beliefPenaltyMPC_beta7, beliefPenaltyMPC_bmy7);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_104(beliefPenaltyMPC_Ld7, beliefPenaltyMPC_bmy7, beliefPenaltyMPC_yy7);
beliefPenaltyMPC_LA_VSUB6_INDEXED_104_104_104(beliefPenaltyMPC_ccrhsub9, beliefPenaltyMPC_sub9, beliefPenaltyMPC_ubIdx9, beliefPenaltyMPC_ccrhsl9, beliefPenaltyMPC_slb9, beliefPenaltyMPC_lbIdx9, beliefPenaltyMPC_rd9);
beliefPenaltyMPC_LA_DIAG_FORWARDSUB_104(beliefPenaltyMPC_Phi9, beliefPenaltyMPC_rd9, beliefPenaltyMPC_Lbyrd9);
beliefPenaltyMPC_LA_DENSE_2MVMADD_104_325_104(beliefPenaltyMPC_V8, beliefPenaltyMPC_Lbyrd8, beliefPenaltyMPC_W9, beliefPenaltyMPC_Lbyrd9, beliefPenaltyMPC_beta8);
beliefPenaltyMPC_LA_DENSE_MVMSUB1_104_104(beliefPenaltyMPC_Lsd8, beliefPenaltyMPC_yy7, beliefPenaltyMPC_beta8, beliefPenaltyMPC_bmy8);
beliefPenaltyMPC_LA_DENSE_FORWARDSUB_104(beliefPenaltyMPC_Ld8, beliefPenaltyMPC_bmy8, beliefPenaltyMPC_yy8);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_104(beliefPenaltyMPC_Ld8, beliefPenaltyMPC_yy8, beliefPenaltyMPC_dvcc8);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_104_104(beliefPenaltyMPC_Lsd8, beliefPenaltyMPC_dvcc8, beliefPenaltyMPC_yy7, beliefPenaltyMPC_bmy7);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_104(beliefPenaltyMPC_Ld7, beliefPenaltyMPC_bmy7, beliefPenaltyMPC_dvcc7);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_104_104(beliefPenaltyMPC_Lsd7, beliefPenaltyMPC_dvcc7, beliefPenaltyMPC_yy6, beliefPenaltyMPC_bmy6);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_104(beliefPenaltyMPC_Ld6, beliefPenaltyMPC_bmy6, beliefPenaltyMPC_dvcc6);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_104_104(beliefPenaltyMPC_Lsd6, beliefPenaltyMPC_dvcc6, beliefPenaltyMPC_yy5, beliefPenaltyMPC_bmy5);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_104(beliefPenaltyMPC_Ld5, beliefPenaltyMPC_bmy5, beliefPenaltyMPC_dvcc5);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_104_104(beliefPenaltyMPC_Lsd5, beliefPenaltyMPC_dvcc5, beliefPenaltyMPC_yy4, beliefPenaltyMPC_bmy4);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_104(beliefPenaltyMPC_Ld4, beliefPenaltyMPC_bmy4, beliefPenaltyMPC_dvcc4);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_104_104(beliefPenaltyMPC_Lsd4, beliefPenaltyMPC_dvcc4, beliefPenaltyMPC_yy3, beliefPenaltyMPC_bmy3);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_104(beliefPenaltyMPC_Ld3, beliefPenaltyMPC_bmy3, beliefPenaltyMPC_dvcc3);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_104_104(beliefPenaltyMPC_Lsd3, beliefPenaltyMPC_dvcc3, beliefPenaltyMPC_yy2, beliefPenaltyMPC_bmy2);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_104(beliefPenaltyMPC_Ld2, beliefPenaltyMPC_bmy2, beliefPenaltyMPC_dvcc2);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_104_104(beliefPenaltyMPC_Lsd2, beliefPenaltyMPC_dvcc2, beliefPenaltyMPC_yy1, beliefPenaltyMPC_bmy1);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_104(beliefPenaltyMPC_Ld1, beliefPenaltyMPC_bmy1, beliefPenaltyMPC_dvcc1);
beliefPenaltyMPC_LA_DENSE_MTVMSUB_104_208(beliefPenaltyMPC_Lsd1, beliefPenaltyMPC_dvcc1, beliefPenaltyMPC_yy0, beliefPenaltyMPC_bmy0);
beliefPenaltyMPC_LA_DENSE_BACKWARDSUB_208(beliefPenaltyMPC_Ld0, beliefPenaltyMPC_bmy0, beliefPenaltyMPC_dvcc0);
beliefPenaltyMPC_LA_DENSE_MTVM_208_325(params->C1, beliefPenaltyMPC_dvcc0, beliefPenaltyMPC_grad_eq0);
beliefPenaltyMPC_LA_DENSE_MTVM2_104_325_208(params->C2, beliefPenaltyMPC_dvcc1, params->D2, beliefPenaltyMPC_dvcc0, beliefPenaltyMPC_grad_eq1);
beliefPenaltyMPC_LA_DENSE_MTVM2_104_325_104(params->C3, beliefPenaltyMPC_dvcc2, params->D3, beliefPenaltyMPC_dvcc1, beliefPenaltyMPC_grad_eq2);
beliefPenaltyMPC_LA_DENSE_MTVM2_104_325_104(params->C4, beliefPenaltyMPC_dvcc3, params->D4, beliefPenaltyMPC_dvcc2, beliefPenaltyMPC_grad_eq3);
beliefPenaltyMPC_LA_DENSE_MTVM2_104_325_104(params->C5, beliefPenaltyMPC_dvcc4, params->D5, beliefPenaltyMPC_dvcc3, beliefPenaltyMPC_grad_eq4);
beliefPenaltyMPC_LA_DENSE_MTVM2_104_325_104(params->C6, beliefPenaltyMPC_dvcc5, params->D6, beliefPenaltyMPC_dvcc4, beliefPenaltyMPC_grad_eq5);
beliefPenaltyMPC_LA_DENSE_MTVM2_104_325_104(params->C7, beliefPenaltyMPC_dvcc6, params->D7, beliefPenaltyMPC_dvcc5, beliefPenaltyMPC_grad_eq6);
beliefPenaltyMPC_LA_DENSE_MTVM2_104_325_104(params->C8, beliefPenaltyMPC_dvcc7, params->D8, beliefPenaltyMPC_dvcc6, beliefPenaltyMPC_grad_eq7);
beliefPenaltyMPC_LA_DENSE_MTVM2_104_325_104(params->C9, beliefPenaltyMPC_dvcc8, params->D9, beliefPenaltyMPC_dvcc7, beliefPenaltyMPC_grad_eq8);
beliefPenaltyMPC_LA_DENSE_MTVM_104_104(params->D10, beliefPenaltyMPC_dvcc8, beliefPenaltyMPC_grad_eq9);
beliefPenaltyMPC_LA_VSUB_3029(beliefPenaltyMPC_rd, beliefPenaltyMPC_grad_eq, beliefPenaltyMPC_rd);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_325(beliefPenaltyMPC_Phi0, beliefPenaltyMPC_rd0, beliefPenaltyMPC_dzcc0);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_325(beliefPenaltyMPC_Phi1, beliefPenaltyMPC_rd1, beliefPenaltyMPC_dzcc1);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_325(beliefPenaltyMPC_Phi2, beliefPenaltyMPC_rd2, beliefPenaltyMPC_dzcc2);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_325(beliefPenaltyMPC_Phi3, beliefPenaltyMPC_rd3, beliefPenaltyMPC_dzcc3);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_325(beliefPenaltyMPC_Phi4, beliefPenaltyMPC_rd4, beliefPenaltyMPC_dzcc4);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_325(beliefPenaltyMPC_Phi5, beliefPenaltyMPC_rd5, beliefPenaltyMPC_dzcc5);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_325(beliefPenaltyMPC_Phi6, beliefPenaltyMPC_rd6, beliefPenaltyMPC_dzcc6);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_325(beliefPenaltyMPC_Phi7, beliefPenaltyMPC_rd7, beliefPenaltyMPC_dzcc7);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_325(beliefPenaltyMPC_Phi8, beliefPenaltyMPC_rd8, beliefPenaltyMPC_dzcc8);
beliefPenaltyMPC_LA_DIAG_FORWARDBACKWARDSUB_104(beliefPenaltyMPC_Phi9, beliefPenaltyMPC_rd9, beliefPenaltyMPC_dzcc9);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_325(beliefPenaltyMPC_ccrhsl0, beliefPenaltyMPC_slb0, beliefPenaltyMPC_llbbyslb0, beliefPenaltyMPC_dzcc0, beliefPenaltyMPC_lbIdx0, beliefPenaltyMPC_dllbcc0);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_117(beliefPenaltyMPC_ccrhsub0, beliefPenaltyMPC_sub0, beliefPenaltyMPC_lubbysub0, beliefPenaltyMPC_dzcc0, beliefPenaltyMPC_ubIdx0, beliefPenaltyMPC_dlubcc0);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_325(beliefPenaltyMPC_ccrhsl1, beliefPenaltyMPC_slb1, beliefPenaltyMPC_llbbyslb1, beliefPenaltyMPC_dzcc1, beliefPenaltyMPC_lbIdx1, beliefPenaltyMPC_dllbcc1);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_117(beliefPenaltyMPC_ccrhsub1, beliefPenaltyMPC_sub1, beliefPenaltyMPC_lubbysub1, beliefPenaltyMPC_dzcc1, beliefPenaltyMPC_ubIdx1, beliefPenaltyMPC_dlubcc1);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_325(beliefPenaltyMPC_ccrhsl2, beliefPenaltyMPC_slb2, beliefPenaltyMPC_llbbyslb2, beliefPenaltyMPC_dzcc2, beliefPenaltyMPC_lbIdx2, beliefPenaltyMPC_dllbcc2);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_117(beliefPenaltyMPC_ccrhsub2, beliefPenaltyMPC_sub2, beliefPenaltyMPC_lubbysub2, beliefPenaltyMPC_dzcc2, beliefPenaltyMPC_ubIdx2, beliefPenaltyMPC_dlubcc2);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_325(beliefPenaltyMPC_ccrhsl3, beliefPenaltyMPC_slb3, beliefPenaltyMPC_llbbyslb3, beliefPenaltyMPC_dzcc3, beliefPenaltyMPC_lbIdx3, beliefPenaltyMPC_dllbcc3);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_117(beliefPenaltyMPC_ccrhsub3, beliefPenaltyMPC_sub3, beliefPenaltyMPC_lubbysub3, beliefPenaltyMPC_dzcc3, beliefPenaltyMPC_ubIdx3, beliefPenaltyMPC_dlubcc3);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_325(beliefPenaltyMPC_ccrhsl4, beliefPenaltyMPC_slb4, beliefPenaltyMPC_llbbyslb4, beliefPenaltyMPC_dzcc4, beliefPenaltyMPC_lbIdx4, beliefPenaltyMPC_dllbcc4);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_117(beliefPenaltyMPC_ccrhsub4, beliefPenaltyMPC_sub4, beliefPenaltyMPC_lubbysub4, beliefPenaltyMPC_dzcc4, beliefPenaltyMPC_ubIdx4, beliefPenaltyMPC_dlubcc4);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_325(beliefPenaltyMPC_ccrhsl5, beliefPenaltyMPC_slb5, beliefPenaltyMPC_llbbyslb5, beliefPenaltyMPC_dzcc5, beliefPenaltyMPC_lbIdx5, beliefPenaltyMPC_dllbcc5);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_117(beliefPenaltyMPC_ccrhsub5, beliefPenaltyMPC_sub5, beliefPenaltyMPC_lubbysub5, beliefPenaltyMPC_dzcc5, beliefPenaltyMPC_ubIdx5, beliefPenaltyMPC_dlubcc5);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_325(beliefPenaltyMPC_ccrhsl6, beliefPenaltyMPC_slb6, beliefPenaltyMPC_llbbyslb6, beliefPenaltyMPC_dzcc6, beliefPenaltyMPC_lbIdx6, beliefPenaltyMPC_dllbcc6);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_117(beliefPenaltyMPC_ccrhsub6, beliefPenaltyMPC_sub6, beliefPenaltyMPC_lubbysub6, beliefPenaltyMPC_dzcc6, beliefPenaltyMPC_ubIdx6, beliefPenaltyMPC_dlubcc6);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_325(beliefPenaltyMPC_ccrhsl7, beliefPenaltyMPC_slb7, beliefPenaltyMPC_llbbyslb7, beliefPenaltyMPC_dzcc7, beliefPenaltyMPC_lbIdx7, beliefPenaltyMPC_dllbcc7);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_117(beliefPenaltyMPC_ccrhsub7, beliefPenaltyMPC_sub7, beliefPenaltyMPC_lubbysub7, beliefPenaltyMPC_dzcc7, beliefPenaltyMPC_ubIdx7, beliefPenaltyMPC_dlubcc7);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_325(beliefPenaltyMPC_ccrhsl8, beliefPenaltyMPC_slb8, beliefPenaltyMPC_llbbyslb8, beliefPenaltyMPC_dzcc8, beliefPenaltyMPC_lbIdx8, beliefPenaltyMPC_dllbcc8);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_117(beliefPenaltyMPC_ccrhsub8, beliefPenaltyMPC_sub8, beliefPenaltyMPC_lubbysub8, beliefPenaltyMPC_dzcc8, beliefPenaltyMPC_ubIdx8, beliefPenaltyMPC_dlubcc8);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_104(beliefPenaltyMPC_ccrhsl9, beliefPenaltyMPC_slb9, beliefPenaltyMPC_llbbyslb9, beliefPenaltyMPC_dzcc9, beliefPenaltyMPC_lbIdx9, beliefPenaltyMPC_dllbcc9);
beliefPenaltyMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_104(beliefPenaltyMPC_ccrhsub9, beliefPenaltyMPC_sub9, beliefPenaltyMPC_lubbysub9, beliefPenaltyMPC_dzcc9, beliefPenaltyMPC_ubIdx9, beliefPenaltyMPC_dlubcc9);
beliefPenaltyMPC_LA_VSUB7_4186(beliefPenaltyMPC_l, beliefPenaltyMPC_ccrhs, beliefPenaltyMPC_s, beliefPenaltyMPC_dl_cc, beliefPenaltyMPC_ds_cc);
beliefPenaltyMPC_LA_VADD_3029(beliefPenaltyMPC_dz_cc, beliefPenaltyMPC_dz_aff);
beliefPenaltyMPC_LA_VADD_1040(beliefPenaltyMPC_dv_cc, beliefPenaltyMPC_dv_aff);
beliefPenaltyMPC_LA_VADD_4186(beliefPenaltyMPC_dl_cc, beliefPenaltyMPC_dl_aff);
beliefPenaltyMPC_LA_VADD_4186(beliefPenaltyMPC_ds_cc, beliefPenaltyMPC_ds_aff);
info->lsit_cc = beliefPenaltyMPC_LINESEARCH_BACKTRACKING_COMBINED(beliefPenaltyMPC_z, beliefPenaltyMPC_v, beliefPenaltyMPC_l, beliefPenaltyMPC_s, beliefPenaltyMPC_dz_cc, beliefPenaltyMPC_dv_cc, beliefPenaltyMPC_dl_cc, beliefPenaltyMPC_ds_cc, &info->step_cc, &info->mu);
if( info->lsit_cc == beliefPenaltyMPC_NOPROGRESS ){
PRINTTEXT("Line search could not proceed at iteration %d, exiting.\n",info->it+1);
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
output->z1[22] = beliefPenaltyMPC_z0[22];
output->z1[23] = beliefPenaltyMPC_z0[23];
output->z1[24] = beliefPenaltyMPC_z0[24];
output->z1[25] = beliefPenaltyMPC_z0[25];
output->z1[26] = beliefPenaltyMPC_z0[26];
output->z1[27] = beliefPenaltyMPC_z0[27];
output->z1[28] = beliefPenaltyMPC_z0[28];
output->z1[29] = beliefPenaltyMPC_z0[29];
output->z1[30] = beliefPenaltyMPC_z0[30];
output->z1[31] = beliefPenaltyMPC_z0[31];
output->z1[32] = beliefPenaltyMPC_z0[32];
output->z1[33] = beliefPenaltyMPC_z0[33];
output->z1[34] = beliefPenaltyMPC_z0[34];
output->z1[35] = beliefPenaltyMPC_z0[35];
output->z1[36] = beliefPenaltyMPC_z0[36];
output->z1[37] = beliefPenaltyMPC_z0[37];
output->z1[38] = beliefPenaltyMPC_z0[38];
output->z1[39] = beliefPenaltyMPC_z0[39];
output->z1[40] = beliefPenaltyMPC_z0[40];
output->z1[41] = beliefPenaltyMPC_z0[41];
output->z1[42] = beliefPenaltyMPC_z0[42];
output->z1[43] = beliefPenaltyMPC_z0[43];
output->z1[44] = beliefPenaltyMPC_z0[44];
output->z1[45] = beliefPenaltyMPC_z0[45];
output->z1[46] = beliefPenaltyMPC_z0[46];
output->z1[47] = beliefPenaltyMPC_z0[47];
output->z1[48] = beliefPenaltyMPC_z0[48];
output->z1[49] = beliefPenaltyMPC_z0[49];
output->z1[50] = beliefPenaltyMPC_z0[50];
output->z1[51] = beliefPenaltyMPC_z0[51];
output->z1[52] = beliefPenaltyMPC_z0[52];
output->z1[53] = beliefPenaltyMPC_z0[53];
output->z1[54] = beliefPenaltyMPC_z0[54];
output->z1[55] = beliefPenaltyMPC_z0[55];
output->z1[56] = beliefPenaltyMPC_z0[56];
output->z1[57] = beliefPenaltyMPC_z0[57];
output->z1[58] = beliefPenaltyMPC_z0[58];
output->z1[59] = beliefPenaltyMPC_z0[59];
output->z1[60] = beliefPenaltyMPC_z0[60];
output->z1[61] = beliefPenaltyMPC_z0[61];
output->z1[62] = beliefPenaltyMPC_z0[62];
output->z1[63] = beliefPenaltyMPC_z0[63];
output->z1[64] = beliefPenaltyMPC_z0[64];
output->z1[65] = beliefPenaltyMPC_z0[65];
output->z1[66] = beliefPenaltyMPC_z0[66];
output->z1[67] = beliefPenaltyMPC_z0[67];
output->z1[68] = beliefPenaltyMPC_z0[68];
output->z1[69] = beliefPenaltyMPC_z0[69];
output->z1[70] = beliefPenaltyMPC_z0[70];
output->z1[71] = beliefPenaltyMPC_z0[71];
output->z1[72] = beliefPenaltyMPC_z0[72];
output->z1[73] = beliefPenaltyMPC_z0[73];
output->z1[74] = beliefPenaltyMPC_z0[74];
output->z1[75] = beliefPenaltyMPC_z0[75];
output->z1[76] = beliefPenaltyMPC_z0[76];
output->z1[77] = beliefPenaltyMPC_z0[77];
output->z1[78] = beliefPenaltyMPC_z0[78];
output->z1[79] = beliefPenaltyMPC_z0[79];
output->z1[80] = beliefPenaltyMPC_z0[80];
output->z1[81] = beliefPenaltyMPC_z0[81];
output->z1[82] = beliefPenaltyMPC_z0[82];
output->z1[83] = beliefPenaltyMPC_z0[83];
output->z1[84] = beliefPenaltyMPC_z0[84];
output->z1[85] = beliefPenaltyMPC_z0[85];
output->z1[86] = beliefPenaltyMPC_z0[86];
output->z1[87] = beliefPenaltyMPC_z0[87];
output->z1[88] = beliefPenaltyMPC_z0[88];
output->z1[89] = beliefPenaltyMPC_z0[89];
output->z1[90] = beliefPenaltyMPC_z0[90];
output->z1[91] = beliefPenaltyMPC_z0[91];
output->z1[92] = beliefPenaltyMPC_z0[92];
output->z1[93] = beliefPenaltyMPC_z0[93];
output->z1[94] = beliefPenaltyMPC_z0[94];
output->z1[95] = beliefPenaltyMPC_z0[95];
output->z1[96] = beliefPenaltyMPC_z0[96];
output->z1[97] = beliefPenaltyMPC_z0[97];
output->z1[98] = beliefPenaltyMPC_z0[98];
output->z1[99] = beliefPenaltyMPC_z0[99];
output->z1[100] = beliefPenaltyMPC_z0[100];
output->z1[101] = beliefPenaltyMPC_z0[101];
output->z1[102] = beliefPenaltyMPC_z0[102];
output->z1[103] = beliefPenaltyMPC_z0[103];
output->z1[104] = beliefPenaltyMPC_z0[104];
output->z1[105] = beliefPenaltyMPC_z0[105];
output->z1[106] = beliefPenaltyMPC_z0[106];
output->z1[107] = beliefPenaltyMPC_z0[107];
output->z1[108] = beliefPenaltyMPC_z0[108];
output->z1[109] = beliefPenaltyMPC_z0[109];
output->z1[110] = beliefPenaltyMPC_z0[110];
output->z1[111] = beliefPenaltyMPC_z0[111];
output->z1[112] = beliefPenaltyMPC_z0[112];
output->z1[113] = beliefPenaltyMPC_z0[113];
output->z1[114] = beliefPenaltyMPC_z0[114];
output->z1[115] = beliefPenaltyMPC_z0[115];
output->z1[116] = beliefPenaltyMPC_z0[116];
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
output->z2[22] = beliefPenaltyMPC_z1[22];
output->z2[23] = beliefPenaltyMPC_z1[23];
output->z2[24] = beliefPenaltyMPC_z1[24];
output->z2[25] = beliefPenaltyMPC_z1[25];
output->z2[26] = beliefPenaltyMPC_z1[26];
output->z2[27] = beliefPenaltyMPC_z1[27];
output->z2[28] = beliefPenaltyMPC_z1[28];
output->z2[29] = beliefPenaltyMPC_z1[29];
output->z2[30] = beliefPenaltyMPC_z1[30];
output->z2[31] = beliefPenaltyMPC_z1[31];
output->z2[32] = beliefPenaltyMPC_z1[32];
output->z2[33] = beliefPenaltyMPC_z1[33];
output->z2[34] = beliefPenaltyMPC_z1[34];
output->z2[35] = beliefPenaltyMPC_z1[35];
output->z2[36] = beliefPenaltyMPC_z1[36];
output->z2[37] = beliefPenaltyMPC_z1[37];
output->z2[38] = beliefPenaltyMPC_z1[38];
output->z2[39] = beliefPenaltyMPC_z1[39];
output->z2[40] = beliefPenaltyMPC_z1[40];
output->z2[41] = beliefPenaltyMPC_z1[41];
output->z2[42] = beliefPenaltyMPC_z1[42];
output->z2[43] = beliefPenaltyMPC_z1[43];
output->z2[44] = beliefPenaltyMPC_z1[44];
output->z2[45] = beliefPenaltyMPC_z1[45];
output->z2[46] = beliefPenaltyMPC_z1[46];
output->z2[47] = beliefPenaltyMPC_z1[47];
output->z2[48] = beliefPenaltyMPC_z1[48];
output->z2[49] = beliefPenaltyMPC_z1[49];
output->z2[50] = beliefPenaltyMPC_z1[50];
output->z2[51] = beliefPenaltyMPC_z1[51];
output->z2[52] = beliefPenaltyMPC_z1[52];
output->z2[53] = beliefPenaltyMPC_z1[53];
output->z2[54] = beliefPenaltyMPC_z1[54];
output->z2[55] = beliefPenaltyMPC_z1[55];
output->z2[56] = beliefPenaltyMPC_z1[56];
output->z2[57] = beliefPenaltyMPC_z1[57];
output->z2[58] = beliefPenaltyMPC_z1[58];
output->z2[59] = beliefPenaltyMPC_z1[59];
output->z2[60] = beliefPenaltyMPC_z1[60];
output->z2[61] = beliefPenaltyMPC_z1[61];
output->z2[62] = beliefPenaltyMPC_z1[62];
output->z2[63] = beliefPenaltyMPC_z1[63];
output->z2[64] = beliefPenaltyMPC_z1[64];
output->z2[65] = beliefPenaltyMPC_z1[65];
output->z2[66] = beliefPenaltyMPC_z1[66];
output->z2[67] = beliefPenaltyMPC_z1[67];
output->z2[68] = beliefPenaltyMPC_z1[68];
output->z2[69] = beliefPenaltyMPC_z1[69];
output->z2[70] = beliefPenaltyMPC_z1[70];
output->z2[71] = beliefPenaltyMPC_z1[71];
output->z2[72] = beliefPenaltyMPC_z1[72];
output->z2[73] = beliefPenaltyMPC_z1[73];
output->z2[74] = beliefPenaltyMPC_z1[74];
output->z2[75] = beliefPenaltyMPC_z1[75];
output->z2[76] = beliefPenaltyMPC_z1[76];
output->z2[77] = beliefPenaltyMPC_z1[77];
output->z2[78] = beliefPenaltyMPC_z1[78];
output->z2[79] = beliefPenaltyMPC_z1[79];
output->z2[80] = beliefPenaltyMPC_z1[80];
output->z2[81] = beliefPenaltyMPC_z1[81];
output->z2[82] = beliefPenaltyMPC_z1[82];
output->z2[83] = beliefPenaltyMPC_z1[83];
output->z2[84] = beliefPenaltyMPC_z1[84];
output->z2[85] = beliefPenaltyMPC_z1[85];
output->z2[86] = beliefPenaltyMPC_z1[86];
output->z2[87] = beliefPenaltyMPC_z1[87];
output->z2[88] = beliefPenaltyMPC_z1[88];
output->z2[89] = beliefPenaltyMPC_z1[89];
output->z2[90] = beliefPenaltyMPC_z1[90];
output->z2[91] = beliefPenaltyMPC_z1[91];
output->z2[92] = beliefPenaltyMPC_z1[92];
output->z2[93] = beliefPenaltyMPC_z1[93];
output->z2[94] = beliefPenaltyMPC_z1[94];
output->z2[95] = beliefPenaltyMPC_z1[95];
output->z2[96] = beliefPenaltyMPC_z1[96];
output->z2[97] = beliefPenaltyMPC_z1[97];
output->z2[98] = beliefPenaltyMPC_z1[98];
output->z2[99] = beliefPenaltyMPC_z1[99];
output->z2[100] = beliefPenaltyMPC_z1[100];
output->z2[101] = beliefPenaltyMPC_z1[101];
output->z2[102] = beliefPenaltyMPC_z1[102];
output->z2[103] = beliefPenaltyMPC_z1[103];
output->z2[104] = beliefPenaltyMPC_z1[104];
output->z2[105] = beliefPenaltyMPC_z1[105];
output->z2[106] = beliefPenaltyMPC_z1[106];
output->z2[107] = beliefPenaltyMPC_z1[107];
output->z2[108] = beliefPenaltyMPC_z1[108];
output->z2[109] = beliefPenaltyMPC_z1[109];
output->z2[110] = beliefPenaltyMPC_z1[110];
output->z2[111] = beliefPenaltyMPC_z1[111];
output->z2[112] = beliefPenaltyMPC_z1[112];
output->z2[113] = beliefPenaltyMPC_z1[113];
output->z2[114] = beliefPenaltyMPC_z1[114];
output->z2[115] = beliefPenaltyMPC_z1[115];
output->z2[116] = beliefPenaltyMPC_z1[116];
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
output->z3[22] = beliefPenaltyMPC_z2[22];
output->z3[23] = beliefPenaltyMPC_z2[23];
output->z3[24] = beliefPenaltyMPC_z2[24];
output->z3[25] = beliefPenaltyMPC_z2[25];
output->z3[26] = beliefPenaltyMPC_z2[26];
output->z3[27] = beliefPenaltyMPC_z2[27];
output->z3[28] = beliefPenaltyMPC_z2[28];
output->z3[29] = beliefPenaltyMPC_z2[29];
output->z3[30] = beliefPenaltyMPC_z2[30];
output->z3[31] = beliefPenaltyMPC_z2[31];
output->z3[32] = beliefPenaltyMPC_z2[32];
output->z3[33] = beliefPenaltyMPC_z2[33];
output->z3[34] = beliefPenaltyMPC_z2[34];
output->z3[35] = beliefPenaltyMPC_z2[35];
output->z3[36] = beliefPenaltyMPC_z2[36];
output->z3[37] = beliefPenaltyMPC_z2[37];
output->z3[38] = beliefPenaltyMPC_z2[38];
output->z3[39] = beliefPenaltyMPC_z2[39];
output->z3[40] = beliefPenaltyMPC_z2[40];
output->z3[41] = beliefPenaltyMPC_z2[41];
output->z3[42] = beliefPenaltyMPC_z2[42];
output->z3[43] = beliefPenaltyMPC_z2[43];
output->z3[44] = beliefPenaltyMPC_z2[44];
output->z3[45] = beliefPenaltyMPC_z2[45];
output->z3[46] = beliefPenaltyMPC_z2[46];
output->z3[47] = beliefPenaltyMPC_z2[47];
output->z3[48] = beliefPenaltyMPC_z2[48];
output->z3[49] = beliefPenaltyMPC_z2[49];
output->z3[50] = beliefPenaltyMPC_z2[50];
output->z3[51] = beliefPenaltyMPC_z2[51];
output->z3[52] = beliefPenaltyMPC_z2[52];
output->z3[53] = beliefPenaltyMPC_z2[53];
output->z3[54] = beliefPenaltyMPC_z2[54];
output->z3[55] = beliefPenaltyMPC_z2[55];
output->z3[56] = beliefPenaltyMPC_z2[56];
output->z3[57] = beliefPenaltyMPC_z2[57];
output->z3[58] = beliefPenaltyMPC_z2[58];
output->z3[59] = beliefPenaltyMPC_z2[59];
output->z3[60] = beliefPenaltyMPC_z2[60];
output->z3[61] = beliefPenaltyMPC_z2[61];
output->z3[62] = beliefPenaltyMPC_z2[62];
output->z3[63] = beliefPenaltyMPC_z2[63];
output->z3[64] = beliefPenaltyMPC_z2[64];
output->z3[65] = beliefPenaltyMPC_z2[65];
output->z3[66] = beliefPenaltyMPC_z2[66];
output->z3[67] = beliefPenaltyMPC_z2[67];
output->z3[68] = beliefPenaltyMPC_z2[68];
output->z3[69] = beliefPenaltyMPC_z2[69];
output->z3[70] = beliefPenaltyMPC_z2[70];
output->z3[71] = beliefPenaltyMPC_z2[71];
output->z3[72] = beliefPenaltyMPC_z2[72];
output->z3[73] = beliefPenaltyMPC_z2[73];
output->z3[74] = beliefPenaltyMPC_z2[74];
output->z3[75] = beliefPenaltyMPC_z2[75];
output->z3[76] = beliefPenaltyMPC_z2[76];
output->z3[77] = beliefPenaltyMPC_z2[77];
output->z3[78] = beliefPenaltyMPC_z2[78];
output->z3[79] = beliefPenaltyMPC_z2[79];
output->z3[80] = beliefPenaltyMPC_z2[80];
output->z3[81] = beliefPenaltyMPC_z2[81];
output->z3[82] = beliefPenaltyMPC_z2[82];
output->z3[83] = beliefPenaltyMPC_z2[83];
output->z3[84] = beliefPenaltyMPC_z2[84];
output->z3[85] = beliefPenaltyMPC_z2[85];
output->z3[86] = beliefPenaltyMPC_z2[86];
output->z3[87] = beliefPenaltyMPC_z2[87];
output->z3[88] = beliefPenaltyMPC_z2[88];
output->z3[89] = beliefPenaltyMPC_z2[89];
output->z3[90] = beliefPenaltyMPC_z2[90];
output->z3[91] = beliefPenaltyMPC_z2[91];
output->z3[92] = beliefPenaltyMPC_z2[92];
output->z3[93] = beliefPenaltyMPC_z2[93];
output->z3[94] = beliefPenaltyMPC_z2[94];
output->z3[95] = beliefPenaltyMPC_z2[95];
output->z3[96] = beliefPenaltyMPC_z2[96];
output->z3[97] = beliefPenaltyMPC_z2[97];
output->z3[98] = beliefPenaltyMPC_z2[98];
output->z3[99] = beliefPenaltyMPC_z2[99];
output->z3[100] = beliefPenaltyMPC_z2[100];
output->z3[101] = beliefPenaltyMPC_z2[101];
output->z3[102] = beliefPenaltyMPC_z2[102];
output->z3[103] = beliefPenaltyMPC_z2[103];
output->z3[104] = beliefPenaltyMPC_z2[104];
output->z3[105] = beliefPenaltyMPC_z2[105];
output->z3[106] = beliefPenaltyMPC_z2[106];
output->z3[107] = beliefPenaltyMPC_z2[107];
output->z3[108] = beliefPenaltyMPC_z2[108];
output->z3[109] = beliefPenaltyMPC_z2[109];
output->z3[110] = beliefPenaltyMPC_z2[110];
output->z3[111] = beliefPenaltyMPC_z2[111];
output->z3[112] = beliefPenaltyMPC_z2[112];
output->z3[113] = beliefPenaltyMPC_z2[113];
output->z3[114] = beliefPenaltyMPC_z2[114];
output->z3[115] = beliefPenaltyMPC_z2[115];
output->z3[116] = beliefPenaltyMPC_z2[116];
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
output->z4[22] = beliefPenaltyMPC_z3[22];
output->z4[23] = beliefPenaltyMPC_z3[23];
output->z4[24] = beliefPenaltyMPC_z3[24];
output->z4[25] = beliefPenaltyMPC_z3[25];
output->z4[26] = beliefPenaltyMPC_z3[26];
output->z4[27] = beliefPenaltyMPC_z3[27];
output->z4[28] = beliefPenaltyMPC_z3[28];
output->z4[29] = beliefPenaltyMPC_z3[29];
output->z4[30] = beliefPenaltyMPC_z3[30];
output->z4[31] = beliefPenaltyMPC_z3[31];
output->z4[32] = beliefPenaltyMPC_z3[32];
output->z4[33] = beliefPenaltyMPC_z3[33];
output->z4[34] = beliefPenaltyMPC_z3[34];
output->z4[35] = beliefPenaltyMPC_z3[35];
output->z4[36] = beliefPenaltyMPC_z3[36];
output->z4[37] = beliefPenaltyMPC_z3[37];
output->z4[38] = beliefPenaltyMPC_z3[38];
output->z4[39] = beliefPenaltyMPC_z3[39];
output->z4[40] = beliefPenaltyMPC_z3[40];
output->z4[41] = beliefPenaltyMPC_z3[41];
output->z4[42] = beliefPenaltyMPC_z3[42];
output->z4[43] = beliefPenaltyMPC_z3[43];
output->z4[44] = beliefPenaltyMPC_z3[44];
output->z4[45] = beliefPenaltyMPC_z3[45];
output->z4[46] = beliefPenaltyMPC_z3[46];
output->z4[47] = beliefPenaltyMPC_z3[47];
output->z4[48] = beliefPenaltyMPC_z3[48];
output->z4[49] = beliefPenaltyMPC_z3[49];
output->z4[50] = beliefPenaltyMPC_z3[50];
output->z4[51] = beliefPenaltyMPC_z3[51];
output->z4[52] = beliefPenaltyMPC_z3[52];
output->z4[53] = beliefPenaltyMPC_z3[53];
output->z4[54] = beliefPenaltyMPC_z3[54];
output->z4[55] = beliefPenaltyMPC_z3[55];
output->z4[56] = beliefPenaltyMPC_z3[56];
output->z4[57] = beliefPenaltyMPC_z3[57];
output->z4[58] = beliefPenaltyMPC_z3[58];
output->z4[59] = beliefPenaltyMPC_z3[59];
output->z4[60] = beliefPenaltyMPC_z3[60];
output->z4[61] = beliefPenaltyMPC_z3[61];
output->z4[62] = beliefPenaltyMPC_z3[62];
output->z4[63] = beliefPenaltyMPC_z3[63];
output->z4[64] = beliefPenaltyMPC_z3[64];
output->z4[65] = beliefPenaltyMPC_z3[65];
output->z4[66] = beliefPenaltyMPC_z3[66];
output->z4[67] = beliefPenaltyMPC_z3[67];
output->z4[68] = beliefPenaltyMPC_z3[68];
output->z4[69] = beliefPenaltyMPC_z3[69];
output->z4[70] = beliefPenaltyMPC_z3[70];
output->z4[71] = beliefPenaltyMPC_z3[71];
output->z4[72] = beliefPenaltyMPC_z3[72];
output->z4[73] = beliefPenaltyMPC_z3[73];
output->z4[74] = beliefPenaltyMPC_z3[74];
output->z4[75] = beliefPenaltyMPC_z3[75];
output->z4[76] = beliefPenaltyMPC_z3[76];
output->z4[77] = beliefPenaltyMPC_z3[77];
output->z4[78] = beliefPenaltyMPC_z3[78];
output->z4[79] = beliefPenaltyMPC_z3[79];
output->z4[80] = beliefPenaltyMPC_z3[80];
output->z4[81] = beliefPenaltyMPC_z3[81];
output->z4[82] = beliefPenaltyMPC_z3[82];
output->z4[83] = beliefPenaltyMPC_z3[83];
output->z4[84] = beliefPenaltyMPC_z3[84];
output->z4[85] = beliefPenaltyMPC_z3[85];
output->z4[86] = beliefPenaltyMPC_z3[86];
output->z4[87] = beliefPenaltyMPC_z3[87];
output->z4[88] = beliefPenaltyMPC_z3[88];
output->z4[89] = beliefPenaltyMPC_z3[89];
output->z4[90] = beliefPenaltyMPC_z3[90];
output->z4[91] = beliefPenaltyMPC_z3[91];
output->z4[92] = beliefPenaltyMPC_z3[92];
output->z4[93] = beliefPenaltyMPC_z3[93];
output->z4[94] = beliefPenaltyMPC_z3[94];
output->z4[95] = beliefPenaltyMPC_z3[95];
output->z4[96] = beliefPenaltyMPC_z3[96];
output->z4[97] = beliefPenaltyMPC_z3[97];
output->z4[98] = beliefPenaltyMPC_z3[98];
output->z4[99] = beliefPenaltyMPC_z3[99];
output->z4[100] = beliefPenaltyMPC_z3[100];
output->z4[101] = beliefPenaltyMPC_z3[101];
output->z4[102] = beliefPenaltyMPC_z3[102];
output->z4[103] = beliefPenaltyMPC_z3[103];
output->z4[104] = beliefPenaltyMPC_z3[104];
output->z4[105] = beliefPenaltyMPC_z3[105];
output->z4[106] = beliefPenaltyMPC_z3[106];
output->z4[107] = beliefPenaltyMPC_z3[107];
output->z4[108] = beliefPenaltyMPC_z3[108];
output->z4[109] = beliefPenaltyMPC_z3[109];
output->z4[110] = beliefPenaltyMPC_z3[110];
output->z4[111] = beliefPenaltyMPC_z3[111];
output->z4[112] = beliefPenaltyMPC_z3[112];
output->z4[113] = beliefPenaltyMPC_z3[113];
output->z4[114] = beliefPenaltyMPC_z3[114];
output->z4[115] = beliefPenaltyMPC_z3[115];
output->z4[116] = beliefPenaltyMPC_z3[116];
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
output->z5[22] = beliefPenaltyMPC_z4[22];
output->z5[23] = beliefPenaltyMPC_z4[23];
output->z5[24] = beliefPenaltyMPC_z4[24];
output->z5[25] = beliefPenaltyMPC_z4[25];
output->z5[26] = beliefPenaltyMPC_z4[26];
output->z5[27] = beliefPenaltyMPC_z4[27];
output->z5[28] = beliefPenaltyMPC_z4[28];
output->z5[29] = beliefPenaltyMPC_z4[29];
output->z5[30] = beliefPenaltyMPC_z4[30];
output->z5[31] = beliefPenaltyMPC_z4[31];
output->z5[32] = beliefPenaltyMPC_z4[32];
output->z5[33] = beliefPenaltyMPC_z4[33];
output->z5[34] = beliefPenaltyMPC_z4[34];
output->z5[35] = beliefPenaltyMPC_z4[35];
output->z5[36] = beliefPenaltyMPC_z4[36];
output->z5[37] = beliefPenaltyMPC_z4[37];
output->z5[38] = beliefPenaltyMPC_z4[38];
output->z5[39] = beliefPenaltyMPC_z4[39];
output->z5[40] = beliefPenaltyMPC_z4[40];
output->z5[41] = beliefPenaltyMPC_z4[41];
output->z5[42] = beliefPenaltyMPC_z4[42];
output->z5[43] = beliefPenaltyMPC_z4[43];
output->z5[44] = beliefPenaltyMPC_z4[44];
output->z5[45] = beliefPenaltyMPC_z4[45];
output->z5[46] = beliefPenaltyMPC_z4[46];
output->z5[47] = beliefPenaltyMPC_z4[47];
output->z5[48] = beliefPenaltyMPC_z4[48];
output->z5[49] = beliefPenaltyMPC_z4[49];
output->z5[50] = beliefPenaltyMPC_z4[50];
output->z5[51] = beliefPenaltyMPC_z4[51];
output->z5[52] = beliefPenaltyMPC_z4[52];
output->z5[53] = beliefPenaltyMPC_z4[53];
output->z5[54] = beliefPenaltyMPC_z4[54];
output->z5[55] = beliefPenaltyMPC_z4[55];
output->z5[56] = beliefPenaltyMPC_z4[56];
output->z5[57] = beliefPenaltyMPC_z4[57];
output->z5[58] = beliefPenaltyMPC_z4[58];
output->z5[59] = beliefPenaltyMPC_z4[59];
output->z5[60] = beliefPenaltyMPC_z4[60];
output->z5[61] = beliefPenaltyMPC_z4[61];
output->z5[62] = beliefPenaltyMPC_z4[62];
output->z5[63] = beliefPenaltyMPC_z4[63];
output->z5[64] = beliefPenaltyMPC_z4[64];
output->z5[65] = beliefPenaltyMPC_z4[65];
output->z5[66] = beliefPenaltyMPC_z4[66];
output->z5[67] = beliefPenaltyMPC_z4[67];
output->z5[68] = beliefPenaltyMPC_z4[68];
output->z5[69] = beliefPenaltyMPC_z4[69];
output->z5[70] = beliefPenaltyMPC_z4[70];
output->z5[71] = beliefPenaltyMPC_z4[71];
output->z5[72] = beliefPenaltyMPC_z4[72];
output->z5[73] = beliefPenaltyMPC_z4[73];
output->z5[74] = beliefPenaltyMPC_z4[74];
output->z5[75] = beliefPenaltyMPC_z4[75];
output->z5[76] = beliefPenaltyMPC_z4[76];
output->z5[77] = beliefPenaltyMPC_z4[77];
output->z5[78] = beliefPenaltyMPC_z4[78];
output->z5[79] = beliefPenaltyMPC_z4[79];
output->z5[80] = beliefPenaltyMPC_z4[80];
output->z5[81] = beliefPenaltyMPC_z4[81];
output->z5[82] = beliefPenaltyMPC_z4[82];
output->z5[83] = beliefPenaltyMPC_z4[83];
output->z5[84] = beliefPenaltyMPC_z4[84];
output->z5[85] = beliefPenaltyMPC_z4[85];
output->z5[86] = beliefPenaltyMPC_z4[86];
output->z5[87] = beliefPenaltyMPC_z4[87];
output->z5[88] = beliefPenaltyMPC_z4[88];
output->z5[89] = beliefPenaltyMPC_z4[89];
output->z5[90] = beliefPenaltyMPC_z4[90];
output->z5[91] = beliefPenaltyMPC_z4[91];
output->z5[92] = beliefPenaltyMPC_z4[92];
output->z5[93] = beliefPenaltyMPC_z4[93];
output->z5[94] = beliefPenaltyMPC_z4[94];
output->z5[95] = beliefPenaltyMPC_z4[95];
output->z5[96] = beliefPenaltyMPC_z4[96];
output->z5[97] = beliefPenaltyMPC_z4[97];
output->z5[98] = beliefPenaltyMPC_z4[98];
output->z5[99] = beliefPenaltyMPC_z4[99];
output->z5[100] = beliefPenaltyMPC_z4[100];
output->z5[101] = beliefPenaltyMPC_z4[101];
output->z5[102] = beliefPenaltyMPC_z4[102];
output->z5[103] = beliefPenaltyMPC_z4[103];
output->z5[104] = beliefPenaltyMPC_z4[104];
output->z5[105] = beliefPenaltyMPC_z4[105];
output->z5[106] = beliefPenaltyMPC_z4[106];
output->z5[107] = beliefPenaltyMPC_z4[107];
output->z5[108] = beliefPenaltyMPC_z4[108];
output->z5[109] = beliefPenaltyMPC_z4[109];
output->z5[110] = beliefPenaltyMPC_z4[110];
output->z5[111] = beliefPenaltyMPC_z4[111];
output->z5[112] = beliefPenaltyMPC_z4[112];
output->z5[113] = beliefPenaltyMPC_z4[113];
output->z5[114] = beliefPenaltyMPC_z4[114];
output->z5[115] = beliefPenaltyMPC_z4[115];
output->z5[116] = beliefPenaltyMPC_z4[116];
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
output->z6[22] = beliefPenaltyMPC_z5[22];
output->z6[23] = beliefPenaltyMPC_z5[23];
output->z6[24] = beliefPenaltyMPC_z5[24];
output->z6[25] = beliefPenaltyMPC_z5[25];
output->z6[26] = beliefPenaltyMPC_z5[26];
output->z6[27] = beliefPenaltyMPC_z5[27];
output->z6[28] = beliefPenaltyMPC_z5[28];
output->z6[29] = beliefPenaltyMPC_z5[29];
output->z6[30] = beliefPenaltyMPC_z5[30];
output->z6[31] = beliefPenaltyMPC_z5[31];
output->z6[32] = beliefPenaltyMPC_z5[32];
output->z6[33] = beliefPenaltyMPC_z5[33];
output->z6[34] = beliefPenaltyMPC_z5[34];
output->z6[35] = beliefPenaltyMPC_z5[35];
output->z6[36] = beliefPenaltyMPC_z5[36];
output->z6[37] = beliefPenaltyMPC_z5[37];
output->z6[38] = beliefPenaltyMPC_z5[38];
output->z6[39] = beliefPenaltyMPC_z5[39];
output->z6[40] = beliefPenaltyMPC_z5[40];
output->z6[41] = beliefPenaltyMPC_z5[41];
output->z6[42] = beliefPenaltyMPC_z5[42];
output->z6[43] = beliefPenaltyMPC_z5[43];
output->z6[44] = beliefPenaltyMPC_z5[44];
output->z6[45] = beliefPenaltyMPC_z5[45];
output->z6[46] = beliefPenaltyMPC_z5[46];
output->z6[47] = beliefPenaltyMPC_z5[47];
output->z6[48] = beliefPenaltyMPC_z5[48];
output->z6[49] = beliefPenaltyMPC_z5[49];
output->z6[50] = beliefPenaltyMPC_z5[50];
output->z6[51] = beliefPenaltyMPC_z5[51];
output->z6[52] = beliefPenaltyMPC_z5[52];
output->z6[53] = beliefPenaltyMPC_z5[53];
output->z6[54] = beliefPenaltyMPC_z5[54];
output->z6[55] = beliefPenaltyMPC_z5[55];
output->z6[56] = beliefPenaltyMPC_z5[56];
output->z6[57] = beliefPenaltyMPC_z5[57];
output->z6[58] = beliefPenaltyMPC_z5[58];
output->z6[59] = beliefPenaltyMPC_z5[59];
output->z6[60] = beliefPenaltyMPC_z5[60];
output->z6[61] = beliefPenaltyMPC_z5[61];
output->z6[62] = beliefPenaltyMPC_z5[62];
output->z6[63] = beliefPenaltyMPC_z5[63];
output->z6[64] = beliefPenaltyMPC_z5[64];
output->z6[65] = beliefPenaltyMPC_z5[65];
output->z6[66] = beliefPenaltyMPC_z5[66];
output->z6[67] = beliefPenaltyMPC_z5[67];
output->z6[68] = beliefPenaltyMPC_z5[68];
output->z6[69] = beliefPenaltyMPC_z5[69];
output->z6[70] = beliefPenaltyMPC_z5[70];
output->z6[71] = beliefPenaltyMPC_z5[71];
output->z6[72] = beliefPenaltyMPC_z5[72];
output->z6[73] = beliefPenaltyMPC_z5[73];
output->z6[74] = beliefPenaltyMPC_z5[74];
output->z6[75] = beliefPenaltyMPC_z5[75];
output->z6[76] = beliefPenaltyMPC_z5[76];
output->z6[77] = beliefPenaltyMPC_z5[77];
output->z6[78] = beliefPenaltyMPC_z5[78];
output->z6[79] = beliefPenaltyMPC_z5[79];
output->z6[80] = beliefPenaltyMPC_z5[80];
output->z6[81] = beliefPenaltyMPC_z5[81];
output->z6[82] = beliefPenaltyMPC_z5[82];
output->z6[83] = beliefPenaltyMPC_z5[83];
output->z6[84] = beliefPenaltyMPC_z5[84];
output->z6[85] = beliefPenaltyMPC_z5[85];
output->z6[86] = beliefPenaltyMPC_z5[86];
output->z6[87] = beliefPenaltyMPC_z5[87];
output->z6[88] = beliefPenaltyMPC_z5[88];
output->z6[89] = beliefPenaltyMPC_z5[89];
output->z6[90] = beliefPenaltyMPC_z5[90];
output->z6[91] = beliefPenaltyMPC_z5[91];
output->z6[92] = beliefPenaltyMPC_z5[92];
output->z6[93] = beliefPenaltyMPC_z5[93];
output->z6[94] = beliefPenaltyMPC_z5[94];
output->z6[95] = beliefPenaltyMPC_z5[95];
output->z6[96] = beliefPenaltyMPC_z5[96];
output->z6[97] = beliefPenaltyMPC_z5[97];
output->z6[98] = beliefPenaltyMPC_z5[98];
output->z6[99] = beliefPenaltyMPC_z5[99];
output->z6[100] = beliefPenaltyMPC_z5[100];
output->z6[101] = beliefPenaltyMPC_z5[101];
output->z6[102] = beliefPenaltyMPC_z5[102];
output->z6[103] = beliefPenaltyMPC_z5[103];
output->z6[104] = beliefPenaltyMPC_z5[104];
output->z6[105] = beliefPenaltyMPC_z5[105];
output->z6[106] = beliefPenaltyMPC_z5[106];
output->z6[107] = beliefPenaltyMPC_z5[107];
output->z6[108] = beliefPenaltyMPC_z5[108];
output->z6[109] = beliefPenaltyMPC_z5[109];
output->z6[110] = beliefPenaltyMPC_z5[110];
output->z6[111] = beliefPenaltyMPC_z5[111];
output->z6[112] = beliefPenaltyMPC_z5[112];
output->z6[113] = beliefPenaltyMPC_z5[113];
output->z6[114] = beliefPenaltyMPC_z5[114];
output->z6[115] = beliefPenaltyMPC_z5[115];
output->z6[116] = beliefPenaltyMPC_z5[116];
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
output->z7[22] = beliefPenaltyMPC_z6[22];
output->z7[23] = beliefPenaltyMPC_z6[23];
output->z7[24] = beliefPenaltyMPC_z6[24];
output->z7[25] = beliefPenaltyMPC_z6[25];
output->z7[26] = beliefPenaltyMPC_z6[26];
output->z7[27] = beliefPenaltyMPC_z6[27];
output->z7[28] = beliefPenaltyMPC_z6[28];
output->z7[29] = beliefPenaltyMPC_z6[29];
output->z7[30] = beliefPenaltyMPC_z6[30];
output->z7[31] = beliefPenaltyMPC_z6[31];
output->z7[32] = beliefPenaltyMPC_z6[32];
output->z7[33] = beliefPenaltyMPC_z6[33];
output->z7[34] = beliefPenaltyMPC_z6[34];
output->z7[35] = beliefPenaltyMPC_z6[35];
output->z7[36] = beliefPenaltyMPC_z6[36];
output->z7[37] = beliefPenaltyMPC_z6[37];
output->z7[38] = beliefPenaltyMPC_z6[38];
output->z7[39] = beliefPenaltyMPC_z6[39];
output->z7[40] = beliefPenaltyMPC_z6[40];
output->z7[41] = beliefPenaltyMPC_z6[41];
output->z7[42] = beliefPenaltyMPC_z6[42];
output->z7[43] = beliefPenaltyMPC_z6[43];
output->z7[44] = beliefPenaltyMPC_z6[44];
output->z7[45] = beliefPenaltyMPC_z6[45];
output->z7[46] = beliefPenaltyMPC_z6[46];
output->z7[47] = beliefPenaltyMPC_z6[47];
output->z7[48] = beliefPenaltyMPC_z6[48];
output->z7[49] = beliefPenaltyMPC_z6[49];
output->z7[50] = beliefPenaltyMPC_z6[50];
output->z7[51] = beliefPenaltyMPC_z6[51];
output->z7[52] = beliefPenaltyMPC_z6[52];
output->z7[53] = beliefPenaltyMPC_z6[53];
output->z7[54] = beliefPenaltyMPC_z6[54];
output->z7[55] = beliefPenaltyMPC_z6[55];
output->z7[56] = beliefPenaltyMPC_z6[56];
output->z7[57] = beliefPenaltyMPC_z6[57];
output->z7[58] = beliefPenaltyMPC_z6[58];
output->z7[59] = beliefPenaltyMPC_z6[59];
output->z7[60] = beliefPenaltyMPC_z6[60];
output->z7[61] = beliefPenaltyMPC_z6[61];
output->z7[62] = beliefPenaltyMPC_z6[62];
output->z7[63] = beliefPenaltyMPC_z6[63];
output->z7[64] = beliefPenaltyMPC_z6[64];
output->z7[65] = beliefPenaltyMPC_z6[65];
output->z7[66] = beliefPenaltyMPC_z6[66];
output->z7[67] = beliefPenaltyMPC_z6[67];
output->z7[68] = beliefPenaltyMPC_z6[68];
output->z7[69] = beliefPenaltyMPC_z6[69];
output->z7[70] = beliefPenaltyMPC_z6[70];
output->z7[71] = beliefPenaltyMPC_z6[71];
output->z7[72] = beliefPenaltyMPC_z6[72];
output->z7[73] = beliefPenaltyMPC_z6[73];
output->z7[74] = beliefPenaltyMPC_z6[74];
output->z7[75] = beliefPenaltyMPC_z6[75];
output->z7[76] = beliefPenaltyMPC_z6[76];
output->z7[77] = beliefPenaltyMPC_z6[77];
output->z7[78] = beliefPenaltyMPC_z6[78];
output->z7[79] = beliefPenaltyMPC_z6[79];
output->z7[80] = beliefPenaltyMPC_z6[80];
output->z7[81] = beliefPenaltyMPC_z6[81];
output->z7[82] = beliefPenaltyMPC_z6[82];
output->z7[83] = beliefPenaltyMPC_z6[83];
output->z7[84] = beliefPenaltyMPC_z6[84];
output->z7[85] = beliefPenaltyMPC_z6[85];
output->z7[86] = beliefPenaltyMPC_z6[86];
output->z7[87] = beliefPenaltyMPC_z6[87];
output->z7[88] = beliefPenaltyMPC_z6[88];
output->z7[89] = beliefPenaltyMPC_z6[89];
output->z7[90] = beliefPenaltyMPC_z6[90];
output->z7[91] = beliefPenaltyMPC_z6[91];
output->z7[92] = beliefPenaltyMPC_z6[92];
output->z7[93] = beliefPenaltyMPC_z6[93];
output->z7[94] = beliefPenaltyMPC_z6[94];
output->z7[95] = beliefPenaltyMPC_z6[95];
output->z7[96] = beliefPenaltyMPC_z6[96];
output->z7[97] = beliefPenaltyMPC_z6[97];
output->z7[98] = beliefPenaltyMPC_z6[98];
output->z7[99] = beliefPenaltyMPC_z6[99];
output->z7[100] = beliefPenaltyMPC_z6[100];
output->z7[101] = beliefPenaltyMPC_z6[101];
output->z7[102] = beliefPenaltyMPC_z6[102];
output->z7[103] = beliefPenaltyMPC_z6[103];
output->z7[104] = beliefPenaltyMPC_z6[104];
output->z7[105] = beliefPenaltyMPC_z6[105];
output->z7[106] = beliefPenaltyMPC_z6[106];
output->z7[107] = beliefPenaltyMPC_z6[107];
output->z7[108] = beliefPenaltyMPC_z6[108];
output->z7[109] = beliefPenaltyMPC_z6[109];
output->z7[110] = beliefPenaltyMPC_z6[110];
output->z7[111] = beliefPenaltyMPC_z6[111];
output->z7[112] = beliefPenaltyMPC_z6[112];
output->z7[113] = beliefPenaltyMPC_z6[113];
output->z7[114] = beliefPenaltyMPC_z6[114];
output->z7[115] = beliefPenaltyMPC_z6[115];
output->z7[116] = beliefPenaltyMPC_z6[116];
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
output->z8[22] = beliefPenaltyMPC_z7[22];
output->z8[23] = beliefPenaltyMPC_z7[23];
output->z8[24] = beliefPenaltyMPC_z7[24];
output->z8[25] = beliefPenaltyMPC_z7[25];
output->z8[26] = beliefPenaltyMPC_z7[26];
output->z8[27] = beliefPenaltyMPC_z7[27];
output->z8[28] = beliefPenaltyMPC_z7[28];
output->z8[29] = beliefPenaltyMPC_z7[29];
output->z8[30] = beliefPenaltyMPC_z7[30];
output->z8[31] = beliefPenaltyMPC_z7[31];
output->z8[32] = beliefPenaltyMPC_z7[32];
output->z8[33] = beliefPenaltyMPC_z7[33];
output->z8[34] = beliefPenaltyMPC_z7[34];
output->z8[35] = beliefPenaltyMPC_z7[35];
output->z8[36] = beliefPenaltyMPC_z7[36];
output->z8[37] = beliefPenaltyMPC_z7[37];
output->z8[38] = beliefPenaltyMPC_z7[38];
output->z8[39] = beliefPenaltyMPC_z7[39];
output->z8[40] = beliefPenaltyMPC_z7[40];
output->z8[41] = beliefPenaltyMPC_z7[41];
output->z8[42] = beliefPenaltyMPC_z7[42];
output->z8[43] = beliefPenaltyMPC_z7[43];
output->z8[44] = beliefPenaltyMPC_z7[44];
output->z8[45] = beliefPenaltyMPC_z7[45];
output->z8[46] = beliefPenaltyMPC_z7[46];
output->z8[47] = beliefPenaltyMPC_z7[47];
output->z8[48] = beliefPenaltyMPC_z7[48];
output->z8[49] = beliefPenaltyMPC_z7[49];
output->z8[50] = beliefPenaltyMPC_z7[50];
output->z8[51] = beliefPenaltyMPC_z7[51];
output->z8[52] = beliefPenaltyMPC_z7[52];
output->z8[53] = beliefPenaltyMPC_z7[53];
output->z8[54] = beliefPenaltyMPC_z7[54];
output->z8[55] = beliefPenaltyMPC_z7[55];
output->z8[56] = beliefPenaltyMPC_z7[56];
output->z8[57] = beliefPenaltyMPC_z7[57];
output->z8[58] = beliefPenaltyMPC_z7[58];
output->z8[59] = beliefPenaltyMPC_z7[59];
output->z8[60] = beliefPenaltyMPC_z7[60];
output->z8[61] = beliefPenaltyMPC_z7[61];
output->z8[62] = beliefPenaltyMPC_z7[62];
output->z8[63] = beliefPenaltyMPC_z7[63];
output->z8[64] = beliefPenaltyMPC_z7[64];
output->z8[65] = beliefPenaltyMPC_z7[65];
output->z8[66] = beliefPenaltyMPC_z7[66];
output->z8[67] = beliefPenaltyMPC_z7[67];
output->z8[68] = beliefPenaltyMPC_z7[68];
output->z8[69] = beliefPenaltyMPC_z7[69];
output->z8[70] = beliefPenaltyMPC_z7[70];
output->z8[71] = beliefPenaltyMPC_z7[71];
output->z8[72] = beliefPenaltyMPC_z7[72];
output->z8[73] = beliefPenaltyMPC_z7[73];
output->z8[74] = beliefPenaltyMPC_z7[74];
output->z8[75] = beliefPenaltyMPC_z7[75];
output->z8[76] = beliefPenaltyMPC_z7[76];
output->z8[77] = beliefPenaltyMPC_z7[77];
output->z8[78] = beliefPenaltyMPC_z7[78];
output->z8[79] = beliefPenaltyMPC_z7[79];
output->z8[80] = beliefPenaltyMPC_z7[80];
output->z8[81] = beliefPenaltyMPC_z7[81];
output->z8[82] = beliefPenaltyMPC_z7[82];
output->z8[83] = beliefPenaltyMPC_z7[83];
output->z8[84] = beliefPenaltyMPC_z7[84];
output->z8[85] = beliefPenaltyMPC_z7[85];
output->z8[86] = beliefPenaltyMPC_z7[86];
output->z8[87] = beliefPenaltyMPC_z7[87];
output->z8[88] = beliefPenaltyMPC_z7[88];
output->z8[89] = beliefPenaltyMPC_z7[89];
output->z8[90] = beliefPenaltyMPC_z7[90];
output->z8[91] = beliefPenaltyMPC_z7[91];
output->z8[92] = beliefPenaltyMPC_z7[92];
output->z8[93] = beliefPenaltyMPC_z7[93];
output->z8[94] = beliefPenaltyMPC_z7[94];
output->z8[95] = beliefPenaltyMPC_z7[95];
output->z8[96] = beliefPenaltyMPC_z7[96];
output->z8[97] = beliefPenaltyMPC_z7[97];
output->z8[98] = beliefPenaltyMPC_z7[98];
output->z8[99] = beliefPenaltyMPC_z7[99];
output->z8[100] = beliefPenaltyMPC_z7[100];
output->z8[101] = beliefPenaltyMPC_z7[101];
output->z8[102] = beliefPenaltyMPC_z7[102];
output->z8[103] = beliefPenaltyMPC_z7[103];
output->z8[104] = beliefPenaltyMPC_z7[104];
output->z8[105] = beliefPenaltyMPC_z7[105];
output->z8[106] = beliefPenaltyMPC_z7[106];
output->z8[107] = beliefPenaltyMPC_z7[107];
output->z8[108] = beliefPenaltyMPC_z7[108];
output->z8[109] = beliefPenaltyMPC_z7[109];
output->z8[110] = beliefPenaltyMPC_z7[110];
output->z8[111] = beliefPenaltyMPC_z7[111];
output->z8[112] = beliefPenaltyMPC_z7[112];
output->z8[113] = beliefPenaltyMPC_z7[113];
output->z8[114] = beliefPenaltyMPC_z7[114];
output->z8[115] = beliefPenaltyMPC_z7[115];
output->z8[116] = beliefPenaltyMPC_z7[116];
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
output->z9[22] = beliefPenaltyMPC_z8[22];
output->z9[23] = beliefPenaltyMPC_z8[23];
output->z9[24] = beliefPenaltyMPC_z8[24];
output->z9[25] = beliefPenaltyMPC_z8[25];
output->z9[26] = beliefPenaltyMPC_z8[26];
output->z9[27] = beliefPenaltyMPC_z8[27];
output->z9[28] = beliefPenaltyMPC_z8[28];
output->z9[29] = beliefPenaltyMPC_z8[29];
output->z9[30] = beliefPenaltyMPC_z8[30];
output->z9[31] = beliefPenaltyMPC_z8[31];
output->z9[32] = beliefPenaltyMPC_z8[32];
output->z9[33] = beliefPenaltyMPC_z8[33];
output->z9[34] = beliefPenaltyMPC_z8[34];
output->z9[35] = beliefPenaltyMPC_z8[35];
output->z9[36] = beliefPenaltyMPC_z8[36];
output->z9[37] = beliefPenaltyMPC_z8[37];
output->z9[38] = beliefPenaltyMPC_z8[38];
output->z9[39] = beliefPenaltyMPC_z8[39];
output->z9[40] = beliefPenaltyMPC_z8[40];
output->z9[41] = beliefPenaltyMPC_z8[41];
output->z9[42] = beliefPenaltyMPC_z8[42];
output->z9[43] = beliefPenaltyMPC_z8[43];
output->z9[44] = beliefPenaltyMPC_z8[44];
output->z9[45] = beliefPenaltyMPC_z8[45];
output->z9[46] = beliefPenaltyMPC_z8[46];
output->z9[47] = beliefPenaltyMPC_z8[47];
output->z9[48] = beliefPenaltyMPC_z8[48];
output->z9[49] = beliefPenaltyMPC_z8[49];
output->z9[50] = beliefPenaltyMPC_z8[50];
output->z9[51] = beliefPenaltyMPC_z8[51];
output->z9[52] = beliefPenaltyMPC_z8[52];
output->z9[53] = beliefPenaltyMPC_z8[53];
output->z9[54] = beliefPenaltyMPC_z8[54];
output->z9[55] = beliefPenaltyMPC_z8[55];
output->z9[56] = beliefPenaltyMPC_z8[56];
output->z9[57] = beliefPenaltyMPC_z8[57];
output->z9[58] = beliefPenaltyMPC_z8[58];
output->z9[59] = beliefPenaltyMPC_z8[59];
output->z9[60] = beliefPenaltyMPC_z8[60];
output->z9[61] = beliefPenaltyMPC_z8[61];
output->z9[62] = beliefPenaltyMPC_z8[62];
output->z9[63] = beliefPenaltyMPC_z8[63];
output->z9[64] = beliefPenaltyMPC_z8[64];
output->z9[65] = beliefPenaltyMPC_z8[65];
output->z9[66] = beliefPenaltyMPC_z8[66];
output->z9[67] = beliefPenaltyMPC_z8[67];
output->z9[68] = beliefPenaltyMPC_z8[68];
output->z9[69] = beliefPenaltyMPC_z8[69];
output->z9[70] = beliefPenaltyMPC_z8[70];
output->z9[71] = beliefPenaltyMPC_z8[71];
output->z9[72] = beliefPenaltyMPC_z8[72];
output->z9[73] = beliefPenaltyMPC_z8[73];
output->z9[74] = beliefPenaltyMPC_z8[74];
output->z9[75] = beliefPenaltyMPC_z8[75];
output->z9[76] = beliefPenaltyMPC_z8[76];
output->z9[77] = beliefPenaltyMPC_z8[77];
output->z9[78] = beliefPenaltyMPC_z8[78];
output->z9[79] = beliefPenaltyMPC_z8[79];
output->z9[80] = beliefPenaltyMPC_z8[80];
output->z9[81] = beliefPenaltyMPC_z8[81];
output->z9[82] = beliefPenaltyMPC_z8[82];
output->z9[83] = beliefPenaltyMPC_z8[83];
output->z9[84] = beliefPenaltyMPC_z8[84];
output->z9[85] = beliefPenaltyMPC_z8[85];
output->z9[86] = beliefPenaltyMPC_z8[86];
output->z9[87] = beliefPenaltyMPC_z8[87];
output->z9[88] = beliefPenaltyMPC_z8[88];
output->z9[89] = beliefPenaltyMPC_z8[89];
output->z9[90] = beliefPenaltyMPC_z8[90];
output->z9[91] = beliefPenaltyMPC_z8[91];
output->z9[92] = beliefPenaltyMPC_z8[92];
output->z9[93] = beliefPenaltyMPC_z8[93];
output->z9[94] = beliefPenaltyMPC_z8[94];
output->z9[95] = beliefPenaltyMPC_z8[95];
output->z9[96] = beliefPenaltyMPC_z8[96];
output->z9[97] = beliefPenaltyMPC_z8[97];
output->z9[98] = beliefPenaltyMPC_z8[98];
output->z9[99] = beliefPenaltyMPC_z8[99];
output->z9[100] = beliefPenaltyMPC_z8[100];
output->z9[101] = beliefPenaltyMPC_z8[101];
output->z9[102] = beliefPenaltyMPC_z8[102];
output->z9[103] = beliefPenaltyMPC_z8[103];
output->z9[104] = beliefPenaltyMPC_z8[104];
output->z9[105] = beliefPenaltyMPC_z8[105];
output->z9[106] = beliefPenaltyMPC_z8[106];
output->z9[107] = beliefPenaltyMPC_z8[107];
output->z9[108] = beliefPenaltyMPC_z8[108];
output->z9[109] = beliefPenaltyMPC_z8[109];
output->z9[110] = beliefPenaltyMPC_z8[110];
output->z9[111] = beliefPenaltyMPC_z8[111];
output->z9[112] = beliefPenaltyMPC_z8[112];
output->z9[113] = beliefPenaltyMPC_z8[113];
output->z9[114] = beliefPenaltyMPC_z8[114];
output->z9[115] = beliefPenaltyMPC_z8[115];
output->z9[116] = beliefPenaltyMPC_z8[116];
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
output->z10[20] = beliefPenaltyMPC_z9[20];
output->z10[21] = beliefPenaltyMPC_z9[21];
output->z10[22] = beliefPenaltyMPC_z9[22];
output->z10[23] = beliefPenaltyMPC_z9[23];
output->z10[24] = beliefPenaltyMPC_z9[24];
output->z10[25] = beliefPenaltyMPC_z9[25];
output->z10[26] = beliefPenaltyMPC_z9[26];
output->z10[27] = beliefPenaltyMPC_z9[27];
output->z10[28] = beliefPenaltyMPC_z9[28];
output->z10[29] = beliefPenaltyMPC_z9[29];
output->z10[30] = beliefPenaltyMPC_z9[30];
output->z10[31] = beliefPenaltyMPC_z9[31];
output->z10[32] = beliefPenaltyMPC_z9[32];
output->z10[33] = beliefPenaltyMPC_z9[33];
output->z10[34] = beliefPenaltyMPC_z9[34];
output->z10[35] = beliefPenaltyMPC_z9[35];
output->z10[36] = beliefPenaltyMPC_z9[36];
output->z10[37] = beliefPenaltyMPC_z9[37];
output->z10[38] = beliefPenaltyMPC_z9[38];
output->z10[39] = beliefPenaltyMPC_z9[39];
output->z10[40] = beliefPenaltyMPC_z9[40];
output->z10[41] = beliefPenaltyMPC_z9[41];
output->z10[42] = beliefPenaltyMPC_z9[42];
output->z10[43] = beliefPenaltyMPC_z9[43];
output->z10[44] = beliefPenaltyMPC_z9[44];
output->z10[45] = beliefPenaltyMPC_z9[45];
output->z10[46] = beliefPenaltyMPC_z9[46];
output->z10[47] = beliefPenaltyMPC_z9[47];
output->z10[48] = beliefPenaltyMPC_z9[48];
output->z10[49] = beliefPenaltyMPC_z9[49];
output->z10[50] = beliefPenaltyMPC_z9[50];
output->z10[51] = beliefPenaltyMPC_z9[51];
output->z10[52] = beliefPenaltyMPC_z9[52];
output->z10[53] = beliefPenaltyMPC_z9[53];
output->z10[54] = beliefPenaltyMPC_z9[54];
output->z10[55] = beliefPenaltyMPC_z9[55];
output->z10[56] = beliefPenaltyMPC_z9[56];
output->z10[57] = beliefPenaltyMPC_z9[57];
output->z10[58] = beliefPenaltyMPC_z9[58];
output->z10[59] = beliefPenaltyMPC_z9[59];
output->z10[60] = beliefPenaltyMPC_z9[60];
output->z10[61] = beliefPenaltyMPC_z9[61];
output->z10[62] = beliefPenaltyMPC_z9[62];
output->z10[63] = beliefPenaltyMPC_z9[63];
output->z10[64] = beliefPenaltyMPC_z9[64];
output->z10[65] = beliefPenaltyMPC_z9[65];
output->z10[66] = beliefPenaltyMPC_z9[66];
output->z10[67] = beliefPenaltyMPC_z9[67];
output->z10[68] = beliefPenaltyMPC_z9[68];
output->z10[69] = beliefPenaltyMPC_z9[69];
output->z10[70] = beliefPenaltyMPC_z9[70];
output->z10[71] = beliefPenaltyMPC_z9[71];
output->z10[72] = beliefPenaltyMPC_z9[72];
output->z10[73] = beliefPenaltyMPC_z9[73];
output->z10[74] = beliefPenaltyMPC_z9[74];
output->z10[75] = beliefPenaltyMPC_z9[75];
output->z10[76] = beliefPenaltyMPC_z9[76];
output->z10[77] = beliefPenaltyMPC_z9[77];
output->z10[78] = beliefPenaltyMPC_z9[78];
output->z10[79] = beliefPenaltyMPC_z9[79];
output->z10[80] = beliefPenaltyMPC_z9[80];
output->z10[81] = beliefPenaltyMPC_z9[81];
output->z10[82] = beliefPenaltyMPC_z9[82];
output->z10[83] = beliefPenaltyMPC_z9[83];
output->z10[84] = beliefPenaltyMPC_z9[84];
output->z10[85] = beliefPenaltyMPC_z9[85];
output->z10[86] = beliefPenaltyMPC_z9[86];
output->z10[87] = beliefPenaltyMPC_z9[87];
output->z10[88] = beliefPenaltyMPC_z9[88];
output->z10[89] = beliefPenaltyMPC_z9[89];
output->z10[90] = beliefPenaltyMPC_z9[90];
output->z10[91] = beliefPenaltyMPC_z9[91];
output->z10[92] = beliefPenaltyMPC_z9[92];
output->z10[93] = beliefPenaltyMPC_z9[93];
output->z10[94] = beliefPenaltyMPC_z9[94];
output->z10[95] = beliefPenaltyMPC_z9[95];
output->z10[96] = beliefPenaltyMPC_z9[96];
output->z10[97] = beliefPenaltyMPC_z9[97];
output->z10[98] = beliefPenaltyMPC_z9[98];
output->z10[99] = beliefPenaltyMPC_z9[99];
output->z10[100] = beliefPenaltyMPC_z9[100];
output->z10[101] = beliefPenaltyMPC_z9[101];
output->z10[102] = beliefPenaltyMPC_z9[102];
output->z10[103] = beliefPenaltyMPC_z9[103];

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
