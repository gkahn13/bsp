#ifndef __UTILS_H__
#define __UTILS_H__

#define M_SQRT_PI_2 sqrt(M_PI/2.0)
#define M_1_SQRT2PI 1.0/sqrt(M_PI*2.0)
#define M_SQRT_2_PI sqrt(2.0/M_PI)
#define M_DEG2RAD M_PI/180.0

inline double sqr(double x) {
	return x*x;
}

inline double sigmoid(double x, double mean)
{
	double y = (x - mean);
	double s = (y/sqrt(1+y*y))+1.0;
	//double s = (y/(1.0 + abs(y)))+1.0;

	if (x < mean) 
		return s*0.1;
	else
		return s*10.0;
}

inline double random() {
	return ((double) rand()) / RAND_MAX;
}

/*************************************************************************
Cephes Math Library Release 2.8:  June, 2000
Copyright by Stephen L. Moshier

Contributors:
    * Sergey Bochkanov (ALGLIB project). Translation from C to pseudocode.

See subroutines comments for additional copyrights.

>>> SOURCE LICENSE >>>
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation (www.fsf.org); either version 2 of the 
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

A copy of the GNU General Public License is available at
http://www.fsf.org/licensing/licenses

>>> END OF LICENSE >>>
*************************************************************************/

/*************************************************************************
Natural logarithm of gamma function

Input parameters:
    X       -   argument

Result:
    logarithm of the absolute value of the Gamma(X).

Output parameters:
    SgnGam  -   sign(Gamma(X))

Domain:
    0 < X < 2.55e305
    -2.55e305 < X < 0, X is not an integer.

ACCURACY:
arithmetic      domain        # trials     peak         rms
   IEEE    0, 3                 28000     5.4e-16     1.1e-16
   IEEE    2.718, 2.556e305     40000     3.5e-16     8.3e-17
The error criterion was relative when the function magnitude
was greater than one but absolute when it was less than one.

The following test used the relative error criterion, though
at certain points the relative error could be much higher than
indicated.
   IEEE    -200, -4             10000     4.8e-16     1.3e-16
*************************************************************************/
double lngamma(double x, double& sgngam)
{
    double result;
    double a;
    double b;
    double c;
    double p;
    double q;
    double u;
    double w;
    double z;
    int i;
    double logpi;
    double ls2pi;
    double tmp;

    sgngam = 1;
    logpi = 1.14472988584940017414;
    ls2pi = 0.91893853320467274178;
    if( x < -34.0 )
    {
        q = -x;
        w = lngamma(q, tmp);
        p = floor(q);
        i = (int) floor(p+0.5);
        if( i%2==0 )
        {
            sgngam = -1;
        }
        else
        {
            sgngam = 1;
        }
        z = q-p;
        if( z > 0.5 )
        {
            p = p+1;
            z = p-q;
        }
        z = q*sin(M_PI*z);
        result = logpi-log(z)-w;
        return result;
    }
    if( x < 13 )
    {
        z = 1;
        p = 0;
        u = x;
        while( u >= 3 )
        {
            p = p-1;
            u = x+p;
            z = z*u;
        }
        while( u < 2 )
        {
            z = z/u;
            p = p+1;
            u = x+p;
        }
        if( z < 0 )
        {
            sgngam = -1;
            z = -z;
        }
        else
        {
            sgngam = 1;
        }
        if( u == 2 )
        {
            result = log(z);
            return result;
        }
        p = p-2;
        x = x+p;
        b = -1378.25152569120859100;
        b = -38801.6315134637840924+x*b;
        b = -331612.992738871184744+x*b;
        b = -1162370.97492762307383+x*b;
        b = -1721737.00820839662146+x*b;
        b = -853555.664245765465627+x*b;
        c = 1;
        c = -351.815701436523470549+x*c;
        c = -17064.2106651881159223+x*c;
        c = -220528.590553854454839+x*c;
        c = -1139334.44367982507207+x*c;
        c = -2532523.07177582951285+x*c;
        c = -2018891.41433532773231+x*c;
        p = x*b/c;
        result = log(z)+p;
        return result;
    }
    q = (x-0.5)*log(x)-x+ls2pi;
    if( x > 100000000 )
    {
        result = q;
        return result;
    }
    p = 1/(x*x);
    if( x >= 1000.0 )
    {
        q = q+((7.9365079365079365079365*0.0001*p-2.7777777777777777777778*0.001)*p+0.0833333333333333333333)/x;
    }
    else
    {
        a = 8.11614167470508450300*0.0001;
        a = -5.95061904284301438324*0.0001+p*a;
        a = 7.93650340457716943945*0.0001+p*a;
        a = -2.77777777730099687205*0.001+p*a;
        a = 8.33333333333331927722*0.01+p*a;
        q = q+a/x;
    }
    result = q;
    return result;
}

/*************************************************************************
Complemented incomplete gamma integral

The function is defined by


 igamc(a,x)   =   1 - igam(a,x)

                           inf.
                             -
                    1       | |  -t  a-1
              =   -----     |   e   t   dt.
                   -      | |
                  | (a)    -
                            x


In this implementation both arguments must be positive.
The integral is evaluated by either a power series or
continued fraction expansion, depending on the relative
values of a and x.

ACCURACY:

Tested at random a, x.
               a         x                      Relative error:
arithmetic   domain   domain     # trials      peak         rms
   IEEE     0.5,100   0,100      200000       1.9e-14     1.7e-15
   IEEE     0.01,0.5  0,100      200000       1.4e-13     1.6e-15
*************************************************************************/
// forward declarations
double incompletegamma(double a, double x);
double incompletegammac(double a, double x);

double incompletegammac(double a, double x)
{
    double result;
    double igammaepsilon;
    double igammabignumber;
    double igammabignumberinv;
    double ans;
    double ax;
    double c;
    double yc;
    double r;
    double t;
    double y;
    double z;
    double pk;
    double pkm1;
    double pkm2;
    double qk;
    double qkm1;
    double qkm2;
    double tmp;

    igammaepsilon = 0.000000000000001;
    igammabignumber = 4503599627370496.0;
    igammabignumberinv = 2.22044604925031308085*0.0000000000000001;
    if( x <= 0 || a <= 0 )
    {
        result = 1;
        return result;
    }
    if( x < 1 || x < a )
    {
        result = 1-incompletegamma(a, x);
        return result;
    }
    ax = a*log(x)-x-lngamma(a, tmp);
    if( ax < -709.78271289338399 )
    {
        result = 0;
        return result;
    }
    ax = exp(ax);
    y = 1-a;
    z = x+y+1;
    c = 0;
    pkm2 = 1;
    qkm2 = x;
    pkm1 = x+1;
    qkm1 = z*x;
    ans = pkm1/qkm1;
    do
    {
        c = c+1;
        y = y+1;
        z = z+2;
        yc = y*c;
        pk = pkm1*z-pkm2*yc;
        qk = qkm1*z-qkm2*yc;
        if( qk != 0 )
        {
            r = pk/qk;
            t = fabs((ans-r)/r);
            ans = r;
        }
        else
        {
            t = 1;
        }
        pkm2 = pkm1;
        pkm1 = pk;
        qkm2 = qkm1;
        qkm1 = qk;
        if( fabs(pk) > igammabignumber )
        {
            pkm2 = pkm2*igammabignumberinv;
            pkm1 = pkm1*igammabignumberinv;
            qkm2 = qkm2*igammabignumberinv;
            qkm1 = qkm1*igammabignumberinv;
        }
    }
    while( t > igammaepsilon );
    result = ans*ax;
    return result;
}

/*************************************************************************
Incomplete gamma integral

The function is defined by

                          x
                           -
                  1       | |  -t  a-1
 igam(a,x)  =   -----     |   e   t   dt.
                 -      | |
                | (a)    -
                          0


In this implementation both arguments must be positive.
The integral is evaluated by either a power series or
continued fraction expansion, depending on the relative
values of a and x.

ACCURACY:

                     Relative error:
arithmetic   domain     # trials      peak         rms
   IEEE      0,30       200000       3.6e-14     2.9e-15
   IEEE      0,100      300000       9.9e-14     1.5e-14
*************************************************************************/
double incompletegamma(double a, double x)
{
	//std::cout << "a: " << a << " x: " << x << std::endl;

    double result;
    double igammaepsilon;
    double ans;
    double ax;
    double c;
    double r;
    double tmp;

    igammaepsilon = 0.000000000000001;
    if( x <= 0 || a <= 0 )
    {
        result = 0;
        return result;
    }
    if( x > 1 && x > a )
    {
        result = 1-incompletegammac(a, x);
        return result;
    }
    ax = a*log(x)-x-lngamma(a, tmp);
    if( ax < -709.78271289338399 )
    {
        result = 0;
        return result;
    }
    ax = exp(ax);
    r = a;
    c = 1;
    ans = 1;
    do
    {
        r = r+1;
        c = c*x/r;
        ans = ans+c;
    }
    while( c/ans > igammaepsilon );
    result = ans*ax/a;
    return result;
}

inline double erf(double x)
{
	double t, z, retval;
	z = fabs( x );
	t = 1.0 / ( 1.0 + 0.5 * z );
	retval = t * exp( -z * z - 1.26551223 + t *
		( 1.00002368 + t *
		( 0.37409196 + t *
		( 0.09678418 + t *
		( -0.18628806 + t *
		( 0.27886807 + t *
		( -1.13520398 + t *
		( 1.48851587 + t *
		( -0.82215223 + t *
		0.1708727 ) ) ) ) ) ) ) ) );
	if( x < 0.0 ) return retval - 1.0;
	return 1.0 - retval;
}

inline double pdf(double x) {
	return exp(-0.5*x*x) * M_1_SQRT2PI;
}

inline double cdf(double x) {
	if (x < 0) {
		return 0.5 - 0.5*incompletegamma(0.5, 0.5*x*x);
	} else {
		return 0.5 + 0.5*incompletegamma(0.5, 0.5*x*x);
	}

	//return 0.5*(1.0 + erf(x*M_SQRT1_2));
}

#include "glut.h"
const double DEG2RAD = M_PI/180;

inline double mod2pi(double x) {
	// Returns a value 0 <= x < 2*PI
	double result = fmod(x, 2*M_PI);
	if (result < 0) {
		result += (2*M_PI);
	}
	return result;
}

inline void drawUnitCircle(void *userdef, float time) {
	glLineWidth(4.0f);
	glDisable(GL_LIGHTING);
	glBegin(GL_LINE_LOOP);
	glColor3f(0.0f, 0.0f, 0.0f);
	//glColor3f(0.5f, 0.0f, 1.0f);

	for (int i = 0; i < 36; i++)
	{
		float theta = (float)(i*10.0*DEG2RAD);
		glVertex3f(cosf(theta), sinf(theta), 0.0f);
	}
	glEnd();
	glEnable(GL_LIGHTING);
	glLineWidth(4.0f);
}

// Jacobian df/dx(x,u)
template <size_t _xDim, size_t _uDim>
inline Matrix<_xDim,_xDim> dfdx(double step, Matrix<_xDim> (*f)(double, const Matrix<_xDim>&, const Matrix<_uDim>&), const Matrix<_xDim>& x, const Matrix<_uDim>& u) 
{
	Matrix<_xDim,_xDim> A;
	Matrix<_xDim> xr(x), xl(x);
	for (size_t i = 0; i < _xDim; ++i) {
		xr[i] += step; xl[i] -= step;
		A.insert(0,i, (e.f(step, xr, u) - e.f(step, xl, u)) / (2.0*step));
		xr[i] = xl[i] = x[i];
	}
	return A;
}

// Jacobian df/du(x,u)
template <size_t _xDim, size_t _uDim>
inline Matrix<_xDim,_uDim> dfdu(double step, Matrix<_xDim> (*f)(double, const Matrix<_xDim>&, const Matrix<_uDim>&), const Matrix<_xDim>& x, const Matrix<_uDim>& u) 
{
	Matrix<_xDim,_uDim> B;
	Matrix<_uDim> ur(u), ul(u);
	for (size_t i = 0; i < _uDim; ++i) {
		ur[i] += step; ul[i] -= step;
		B.insert(0,i, (e.f(step, x, ur) - e.f(step, x, ul)) / (2.0*step));
		ur[i] = ul[i] = u[i];
	}
	return B;
}

// Jacobian dh/dx(x)
template <size_t _xDim, size_t _zDim>
inline Matrix<_zDim,_xDim> dhdx(double step, Matrix<_zDim> (*h)(const Matrix<_xDim>&), const Matrix<_xDim>& x) 
{
	Matrix<_zDim,_xDim> H;
	Matrix<_xDim> xr(x), xl(x);
	for (size_t i = 0; i < _xDim; ++i) {
		xr[i] += step; xl[i] -= step;
		H.insert(0,i, (e.h(xr) - e.h(xl)) / (2.0*step));
		xr[i] = xl[i] = x[i];
	}
	return H;
}

inline Matrix<4,1> quatFromRot(const Matrix<3,3>& R) 
{
	double x = R(2,1) - R(1,2);
	double y = R(0,2) - R(2,0);
	double z = R(1,0) - R(0,1);
	double r = sqrt(x*x+y*y+z*z);
	double t = R(0,0) + R(1,1) + R(2,2);
	double angle = atan2(r,t-1);
	if (angle != 0) {
		x /= r;
		y /= r;
		z /= r;
	} else {
		x = 0;
		y = 0;
		z = 0;
	}
	Matrix<4,1> q;
	q(0,0) = sin(angle/2)*x;
	q(1,0) = sin(angle/2)*y;
	q(2,0) = sin(angle/2)*z;
	q(3,0) = cos(angle/2);

	return q;
}

/*!
 *  @brief       Randomly samples from the uniform distribution over range [0,1].
 *  @returns     A uniform random number in [0,1].
 *  \ingroup globalfunc
 */
inline double random_highprecision() {
  return ((rand() << 15) + rand()) / 1073741823.0;
}

/*!
 *  @brief       Randomly samples from the univariate standard Gaussian distribution N(0,1).
 *  @returns     A random number from the univariate standard Gaussian distribution N(0,1).
 *  \ingroup globalfunc
 */
inline std::pair<double, double> normal() {
  double u, v, s;

  do {
    u = 2.0*random_highprecision()-1.0;
    v = 2.0*random_highprecision()-1.0;
    s = u*u + v*v;
  } while (s > 1.0 || s == 0.0);

  double r = sqrt(-2.0*log(s)/s);

  return std::make_pair(u*r, v*r);
}

/*!
 *  @brief       Randomly samples from the multivariate standard Gaussian distribution N(0,I).
 *  @tparam      dim     The dimension of the distribution.
 *  @returns     A random vector from the multivariate standard Gaussian distribution N(0,I).
 *  \ingroup globalfunc
 */
template <size_t dim>
inline Matrix<dim> sampleGaussian() {
  Matrix<dim> sample;
  for (int j = 0; j < dim / 2; ++j) {
	std::pair<double,double> n = normal();
    sample[j*2] = n.first;
	sample[j*2+1] = n.second;
  }
  if (dim % 2 == 1) {
    sample[dim - 1] = normal().first;
  }
  return sample;
}

/*!
 *  @brief       Randomly samples from the multivariate Gaussian distribution with specified 
 *               mean and variance.
 *  @tparam      dim     The dimension of the distribution.
 *  @param       mean    The mean of the distribution.
 *  @param       var     The variance (covariance matrix) of the distribution.
 *  @returns     A random vector from the specified Gaussian distribution.
 *  \ingroup globalfunc
 */
template <size_t dim>
inline Matrix<dim> sampleGaussian(const Matrix<dim>& mean, const SymmetricMatrix<dim>& var) {
  Matrix<dim> sample = sampleGaussian<dim>();
  Matrix<dim,dim> L;
  chol(var, L);
  return L * sample + mean;
}

template <typename T>
inline int sgn(T val) {
	return(T(0) < val) - (val < T(0));
}

#endif