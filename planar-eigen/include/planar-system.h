#ifndef __PLANAR_SYSTEM_H__
#define __PLANAR_SYSTEM_H__

#include <Python.h>

#include "planar-utils.h"
#include "geometry2d.h"

#include <boost/python.hpp>
#include <boost/python/numeric.hpp>
#include <boost/python/tuple.hpp>
#include <boost/numpy.hpp>
#include <boost/filesystem.hpp>

namespace py = boost::python;
namespace np = boost::numpy;

#include <Eigen/Eigen>
//#define EIGEN_USE_NEW_STDVECTOR
#include <Eigen/StdVector>
using namespace Eigen;

#include "badiff.h"
using namespace fadbad;
typedef B<double> bdouble;

#include "../../util/logging.h"

#define TIMESTEPS 10
#define DT 1.0 // Note: if you change this, must change the FORCES matlab file
#define J_DIM 4 // joint dimension (three for robot, one for camera)
#define C_DIM 2 // object dimension
#define X_DIM 6
#define U_DIM 4
#define Q_DIM 4
#define Z_DIM 6
#define R_DIM 6

#define TOTAL_VARS (TIMESTEPS*X_DIM + (TIMESTEPS-1)*U_DIM)

#define L_DIM 3 // number of links

template <size_t _dim0, size_t _dim1>
using mat = Matrix<double, _dim0, _dim1>;

template <size_t _dim>
using vec = Matrix<double, _dim, 1>;


//template <size_t _dim0, size_t _dim1>
//using matb = Matrix<bdouble, _dim0, _dim1>;
//
//template <size_t _dim>
//using vecb = Matrix<bdouble, _dim, 1>;

class PlanarSystem {
	const double step = 0.0078125*0.0078125;
	const double INFTY = 1e10;

	const double alpha_control = .01;
	const double alpha_belief = 1;
	const double alpha_final_belief = 1;
	const double alpha_goal = 10;

public:
	PlanarSystem(const vec<2>& camera_origin, const vec<2>& object, bool is_static);

//	template< template <class> class VEC, class _xDim, class _uDim, class _qDim>
	vec<X_DIM> dynfunc(const vec<X_DIM>& x, const vec<U_DIM>& u, const vec<Q_DIM>& q, bool enforce_limits=false);
	vec<Z_DIM> obsfunc(const vec<X_DIM>& x, const vec<2>& object, const vec<R_DIM>& r);

	mat<Z_DIM,Z_DIM> delta_matrix(const vec<X_DIM>& x, const vec<2>& object, const double alpha);

	void belief_dynamics(const vec<X_DIM>& x_t, const mat<X_DIM,X_DIM>& sigma_t, const vec<U_DIM>& u_t, const double alpha,
			vec<X_DIM>& x_tp1, mat<X_DIM,X_DIM>& sigma_tp1);
	void execute_control_step(const vec<X_DIM>& x_t_real, const vec<X_DIM>& x_t_t, const mat<X_DIM,X_DIM>& sigma_t_t, const vec<U_DIM>& u_t,
			vec<X_DIM>& x_tp1_real, vec<X_DIM>& x_tp1_tp1, mat<X_DIM,X_DIM>& sigma_tp1_tp1);

	std::vector<Beam> get_fov(const vec<X_DIM>& x);
	std::vector<Segment> get_link_segments(const vec<X_DIM>& x);

	void display(vec<X_DIM>& x, mat<X_DIM,X_DIM>& sigma, bool pause=true);
	void display(std::vector<vec<X_DIM>, aligned_allocator<vec<X_DIM>>>& X, mat<X_DIM,X_DIM>& sigma0, std::vector<vec<U_DIM>, aligned_allocator<vec<U_DIM>>>& U, const double alpha, bool pause=true);
	void display(std::vector<vec<X_DIM>, aligned_allocator<vec<X_DIM>>>& X, std::vector<mat<X_DIM,X_DIM>, aligned_allocator<mat<X_DIM,X_DIM>>>& S, bool pause=true);

	void get_limits(vec<X_DIM>& x_min, vec<X_DIM>& x_max, vec<U_DIM>& u_min, vec<U_DIM>& u_max);

	double cost(const std::vector<vec<X_DIM>, aligned_allocator<vec<X_DIM>>>& X, const mat<X_DIM,X_DIM>& sigma0, const std::vector<vec<U_DIM>, aligned_allocator<vec<U_DIM>>>& U, const double alpha);
	vec<TOTAL_VARS> cost_grad(std::vector<vec<X_DIM>, aligned_allocator<vec<X_DIM>>>& X, const mat<X_DIM,X_DIM>& sigma0, std::vector<vec<U_DIM>, aligned_allocator<vec<U_DIM>>>& U, const double alpha);

private:
	bool is_static;
	vec<2> camera_origin;
	double camera_fov, camera_max_dist;
	vec<2> object;

	vec<2> robot_origin;
	vec<L_DIM> link_lengths;

	mat<Q_DIM,Q_DIM> Q;
	mat<R_DIM,R_DIM> R;
	vec<X_DIM> x_min, x_max;
	vec<U_DIM> u_min, u_max;

	void init(const vec<2>& camera_origin, const vec<2>& object, bool is_static);

	void linearize_dynfunc(const vec<X_DIM>& x, const vec<U_DIM>& u, const vec<Q_DIM>& q, mat<X_DIM,X_DIM>& A, mat<X_DIM,Q_DIM>& M);
	void linearize_obsfunc(const vec<X_DIM>& x, const vec<R_DIM>& r, mat<Z_DIM,X_DIM>& H, mat<Z_DIM,R_DIM>& N);
};

#endif
