#ifndef ADOLCSUPPORT_H
#define ADOLCSUPPORT_H
#define ADOLC_TAPELESS

#include <adolc/adouble.h>
#include <Eigen/Core>

namespace Eigen {
template<> struct NumTraits<adtl::adouble>
: NumTraits<double> // permits to get the epsilon, dummy_precision, lowest, highest functions
{
	typedef adtl::adouble Real;
	typedef adtl::adouble NonInteger;
	typedef adtl::adouble Nested;

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
namespace adtl {

inline const adouble& conj(const adouble& x) { return x; }
inline const adouble& real(const adouble& x) { return x; }
inline adouble imag(const adouble&) { return 0.; }
inline adouble abs(const adouble& x) { return fabs(x); }
inline adouble abs2(const adouble& x) { return x*x; }

}

#endif // ADOLCSUPPORT_H
