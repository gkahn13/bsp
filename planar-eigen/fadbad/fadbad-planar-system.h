#ifndef __FADBAD_PLANAR_SYSTEM_H__
#define __FADBAD_PLANAR_SYSTEM_H__

#include <Python.h>

#include "fadbad-utils.h"
#include "fadbad-geometry2d.h"

#include <boost/python.hpp>
#include <boost/python/numeric.hpp>
#include <boost/python/tuple.hpp>
#include <boost/numpy.hpp>
#include <boost/filesystem.hpp>

namespace py = boost::python;
namespace np = boost::numpy;

#include <Eigen/Eigen>
#include <Eigen/StdVector>
using namespace Eigen;

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


class FadbadPlanarSystem {
	const bdouble step = 0.0078125*0.0078125;
	const bdouble INFTY = 1e10;

	const bdouble alpha_control = .01;
	const bdouble alpha_belief = 1;
	const bdouble alpha_final_belief = 1;
	const bdouble alpha_goal = 10;

public:
	FadbadPlanarSystem() { };
	void init(const vecb<2>& camera_origin, const vecb<2>& object, bool is_static);
	FadbadPlanarSystem(const vecb<2>& camera_origin, const vecb<2>& object, bool is_static);

	vecb<X_DIM> dynfunc(const vecb<X_DIM>& x, const vecb<U_DIM>& u, const vecb<Q_DIM>& q, bool enforce_limits=false);
	vecb<Z_DIM> obsfunc(const vecb<X_DIM>& x, const vecb<2>& object, const vecb<R_DIM>& r);

	matb<Z_DIM,Z_DIM> delta_matrix(const vecb<X_DIM>& x, const vecb<2>& object, const bdouble alpha);

	void belief_dynamics(const vecb<X_DIM>& x_t, const matb<X_DIM,X_DIM>& sigma_t, const vecb<U_DIM>& u_t, const bdouble alpha,
			vecb<X_DIM>& x_tp1, matb<X_DIM,X_DIM>& sigma_tp1);

	std::vector<BeamB> get_fov(const vecb<X_DIM>& x);
	std::vector<SegmentB> get_link_segments(const vecb<X_DIM>& x);

	void get_limits(vecb<X_DIM>& x_min, vecb<X_DIM>& x_max, vecb<U_DIM>& u_min, vecb<U_DIM>& u_max);

	bdouble cost(const std::vector<vecb<X_DIM>, aligned_allocator<vecb<X_DIM>>>& X, const matb<X_DIM,X_DIM>& sigma0,
			const std::vector<vecb<U_DIM>, aligned_allocator<vecb<U_DIM>>>& U, const bdouble alpha);

private:
	bool is_static;
	vecb<2> camera_origin;
	bdouble camera_fov, camera_max_dist;
	vecb<2> object;

	vecb<2> robot_origin;
	vecb<L_DIM> link_lengths;

	matb<Q_DIM,Q_DIM> Q;
	matb<R_DIM,R_DIM> R;
	vecb<X_DIM> x_min, x_max;
	vecb<U_DIM> u_min, u_max;

	void linearize_dynfunc(const vecb<X_DIM>& x, const vecb<U_DIM>& u, const vecb<Q_DIM>& q, matb<X_DIM,X_DIM>& A, matb<X_DIM,Q_DIM>& M);
	void linearize_obsfunc(const vecb<X_DIM>& x, const vecb<R_DIM>& r, matb<Z_DIM,X_DIM>& H, matb<Z_DIM,R_DIM>& N);
};

#endif
