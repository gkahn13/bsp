#ifndef __FADBAD_UTILS_H__
#define __FADBAD_UTILS_H__

#include <Eigen/Eigen>
using namespace Eigen;

#include <assert.h>

#include "badiff.h"
using namespace fadbad;
typedef B<double> bdouble;
typedef B<int> bint;

typedef Matrix<bdouble,Eigen::Dynamic,1> VectorDynb;
typedef Matrix<bdouble,1,Eigen::Dynamic> RowVectorDynb;
typedef Matrix<bdouble,Eigen::Dynamic,Eigen::Dynamic> MatrixDynb;

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

inline bool fadbad_isnan(const bdouble& x) {
	return (abs(x) > 1e10);
}


#endif
