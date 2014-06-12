#ifndef __PLANAR_SYSTEM_H__
#define __PLANAR_SYSTEM_H__

#include <Python.h>

//#include "planar-utils.h"
//#include "geometry2d.h"

#include <boost/python.hpp>
#include <boost/python/numeric.hpp>
#include <boost/python/tuple.hpp>
#include <boost/numpy.hpp>
#include <boost/filesystem.hpp>

namespace py = boost::python;
namespace np = boost::numpy;

#include <Eigen/Eigen>
using namespace Eigen;

#include "badiff.h"
using namespace fadbad;
typedef B<double> bdouble;

#include "../../util/logging.h"

#define DT 1.0 // Note: if you change this, must change the FORCES matlab file
#define J_DIM 4 // joint dimension (three for robot, one for camera)
#define C_DIM 2 // object dimension
#define X_DIM 6
#define U_DIM 4
#define Q_DIM 4
#define Z_DIM 6
#define R_DIM 6

#define L_DIM 3 // number of links

template <size_t _dim0, size_t _dim1>
using mat = Matrix<double, _dim0, _dim1>;

template <size_t _dim>
using vec = Matrix<double, _dim, 1>;

template <size_t _dim0, size_t _dim1>
using matb = Matrix<bdouble, _dim0, _dim1>;

template <size_t _dim>
using vecb = Matrix<bdouble, _dim, 1>;

class PlanarSystem {
	const double step = 0.0078125*0.0078125;
	const double INFTY = 1e10;

	const double alpha_control = .01;
	const double alpha_belief = 1;
	const double alpha_final_belief = 1;
	const double alpha_goal = 10;

public:
	PlanarSystem(const vec<C_DIM>& camera_origin, const vec<C_DIM>& object, bool is_static);

//	template< template <class> class VEC, class _xDim, class _uDim, class _qDim>
	template<
	vec<_xDim> dynfunc(const VEC<_xDim>& x, const VEC<_uDim>& u, const VEC<_qDim>& q, bool enforce_limits=false);

private:
	bool is_static;
	vec<C_DIM> camera_origin;
	double camera_fov, camera_max_dist;
	vec<C_DIM> object;

	vec<C_DIM> robot_origin;
	vec<L_DIM> link_lengths;

	mat<Q_DIM, Q_DIM> Q;
	mat<R_DIM, R_DIM> R;
	vec<X_DIM> x_min, x_max;
	vec<U_DIM> u_min, u_max;
};

#endif
