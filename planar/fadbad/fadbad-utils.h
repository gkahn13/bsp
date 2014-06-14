#ifndef __FADBAD_UTILS_H__
#define __FADBAD_UTILS_H__

#include <Eigen/Eigen>
using namespace Eigen;

#include <assert.h>

#include "badiff.h"
using namespace fadbad;
typedef B<double> bdouble;

namespace Eigen {
template<> struct NumTraits<bdouble >
: NumTraits<double> // permits to get the epsilon, dummy_precision, lowest, highest functions
{
	typedef bdouble Real;
	typedef bdouble NonInteger;
	typedef bdouble Nested;
	enum {
		IsComplex = 0,
		IsInteger = 0,
		IsSigned = 1,
		RequireInitialization = 1,
		ReadCost = 1,
		AddCost = 3,
		MulCost = 3
	};
};
}

namespace fadbad {
inline const bdouble& conj(const bdouble& x) { return x; }
inline const bdouble& real(const bdouble& x) { return x; }
inline bdouble imag(const bdouble&) { return 0.; }
inline bdouble abs(const bdouble& x) { return (x < 0) ? -x : x; }
inline bdouble abs2(const bdouble& x) { return x*x; }
//inline bdouble sqrt(const bdouble& x)  { return sqrt(x); }
//inline bdouble exp(const bdouble&  x)  { return exp(x); }
//inline bdouble log(const bdouble&  x)  { return log(x); }
//inline bdouble sin(const bdouble&  x)  { return sin(x); }
//inline bdouble cos(const bdouble&  x)  { return cos(x); }
//inline bdouble pow(const bdouble& x, bdouble y)  { return pow(x, y); }
}


const bdouble bepsilon = 1e-5;

template <size_t _dim0, size_t _dim1>
using matb = Matrix<bdouble, _dim0, _dim1>;

template <size_t _dim>
using vecb = Matrix<bdouble, _dim, 1>;

#endif
