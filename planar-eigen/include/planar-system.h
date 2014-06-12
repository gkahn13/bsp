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

//template <size_t _dim0, size_t _dim1>
//using mat = Matrix<double, _dim0, _dim1>;
//
//template <size_t _dim>
//using vec = Matrix<double, _dim, 1>;

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
	PlanarSystem(const VectorXd& camera_origin, const VectorXd& object, bool is_static);

//	template< template <class> class VEC, class _xDim, class _uDim, class _qDim>
	VectorXd dynfunc(const VectorXd& x, const VectorXd& u, const VectorXd& q, bool enforce_limits=false);
	VectorXd obsfunc(const VectorXd& x, const VectorXd& object, const VectorXd& r);

	MatrixXd delta_matrix(const VectorXd& x, const VectorXd& object, const double alpha);

	void belief_dynamics(const VectorXd& x_t, const MatrixXd& sigma_t, const VectorXd& u_t, const double alpha, VectorXd& x_tp1, MatrixXd& sigma_tp1);
	void execute_control_step(const VectorXd& x_t_real, const VectorXd& x_t_t, const MatrixXd& sigma_t_t, const VectorXd& u_t,
			VectorXd& x_tp1_real, VectorXd& x_tp1_tp1, MatrixXd& sigma_tp1_tp1);

	std::vector<Beam> get_fov(const VectorXd& x);
	std::vector<Segment> get_link_segments(const VectorXd& x);

	void display(VectorXd& x, MatrixXd& sigma, bool pause=true);
	void display(std::vector<VectorXd>& X, MatrixXd& sigma0, std::vector<VectorXd>& U, const double alpha, bool pause=true);
	void display(std::vector<VectorXd>& X, std::vector<MatrixXd>& S, bool pause=true);

	void get_limits(VectorXd& x_min, VectorXd& x_max, VectorXd& u_min, VectorXd& u_max);

	double cost(const std::vector<VectorXd>& X, const MatrixXd& sigma0, const std::vector<VectorXd>& U, const double alpha);
	VectorXd cost_grad(std::vector<VectorXd>& X, const MatrixXd& sigma0, std::vector<VectorXd>& U, const double alpha);

private:
	bool is_static;
	VectorXd camera_origin;
	double camera_fov, camera_max_dist;
	VectorXd object;

	VectorXd robot_origin;
	VectorXd link_lengths;

	MatrixXd Q;
	MatrixXd R;
	VectorXd x_min, x_max;
	VectorXd u_min, u_max;

	void init(const VectorXd& camera_origin, const VectorXd& object, bool is_static);

	void linearize_dynfunc(const VectorXd& x, const VectorXd& u, const VectorXd& q, MatrixXd& A, MatrixXd& M);
	void linearize_obsfunc(const VectorXd& x, const VectorXd& r, MatrixXd& H, MatrixXd& N);
};

#endif
