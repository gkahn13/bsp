#include "../include/planar-system.h"

/**
 * Constructors and initializers
 */


/**
 * Public methods
 */

template< template <class> class VEC, class _xDim, class _uDim, class _qDim>
VEC<_xDim> PlanarSystem::dynfunc(const VEC<_xDim>& x, const VEC<_uDim>& u, const VEC<_qDim>& q, bool enforce_limits) {
	VEC<_xDim> x_new = x;

}

void test() {
	Matrix3d x;

}
