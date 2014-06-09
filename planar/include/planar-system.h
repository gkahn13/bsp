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

#include <armadillo>
using namespace arma;

#include "../../util/logging.h"


const double step = 0.0078125*0.0078125;
const double INFTY = 1e10;

class PlanarSystem {
public:
	PlanarSystem(const vec& camera_origin, const vec& object, bool is_static);

	vec dynfunc(const vec& x, const vec& u, const vec& q);
	vec obsfunc(const vec& x_robot, const vec& x_object, const vec& r);

	void belief_dynamics(const vec& x_t, const mat& sigma_t, const vec& u_t, vec& x_tp1, mat& sigma_tp1);
	void execute_control_step(const vec& x_t_real, const vec& x_t, const mat& sigma_t, const vec& u_t,
			vec& x_tp1_real, vec& x_tp1, mat& sigma_tp1);

	std::vector<Beam> get_fov(const vec& x);

	void display(std::vector<vec>& X, std::vector<mat>& S, bool pause=true);

	void get_limits(vec& x_min, vec& x_max, vec& u_min, vec& u_max);

private:
	bool is_static;
	vec camera_origin;
	double camera_fov, camera_max_dist;
	vec object;

	vec robot_origin;
	vec link_lengths;

	int X_DIM, U_DIM, Z_DIM, Q_DIM, R_DIM;
	int J_DIM, C_DIM; // joint dimension + camera dimension = X_DIM
	double DT;
	mat Q, R;
	vec x_min, x_max, u_min, u_max;

	void init(const vec& camera_origin, const vec& object, bool is_static);

	std::vector<Segment> get_link_segments(const vec& x);

	void linearize_dynfunc(const vec& x, const vec& u, const vec& q, mat& A, mat& M);
	void linearize_obsfunc(const vec& x, const vec& r, mat& H, mat& N);
	mat delta_matrix(const vec& x);

};

#endif
