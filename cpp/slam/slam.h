#ifndef __SLAM_H__
#define __SLAM_H__

#include <fstream>
#include <math.h>

#define _USE_MATH_DEFINES

#include "../util/matrix.h"
#include "../util/utils.h"
//#include "util/Timer.h"
#include "../util/logging.h"

#include <Python.h>
//#include <pythonrun.h>
#include <boost/python.hpp>
#include <boost/filesystem.hpp>
#include "landmarks.h"

namespace py = boost::python;

#define TIMESTEPS 10
#define DT 1.0
#define NUM_LANDMARKS 3
#define NUM_WAYPOINTS 5

#define C_DIM 3 // car dimension [x, y, theta]
#define P_DIM 2 // Position dimension [x,y]
#define L_DIM 2*NUM_LANDMARKS

#define X_DIM C_DIM+L_DIM
#define U_DIM 2
#define Z_DIM L_DIM
#define Q_DIM 2
#define R_DIM L_DIM

#define S_DIM (((X_DIM+1)*(X_DIM))/2)
#define B_DIM (X_DIM+S_DIM)

const double length = 4;
const double camera_range = 3;
const double camera_view_angle = 3.1415926535/4;

const double step = 0.0078125*0.0078125;


Matrix<X_DIM> x0;
Matrix<X_DIM,X_DIM> SqrtSigma0;
Matrix<X_DIM> xGoal;
Matrix<X_DIM> xMin, xMax;
Matrix<U_DIM> uMin, uMax;
Matrix<Q_DIM, Q_DIM> Q;
Matrix<R_DIM, R_DIM> R;

const int T = TIMESTEPS;
const double INFTY = 1e10;

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

const double alpha_belief = 10, alpha_final_belief = 10, alpha_control = 1, alpha_goal_state = 1;

namespace config {
const double V = 3;
const double MAXG = 30*M_PI/180.;
const double RATEG = 20*M_PI/180.;
const double WHEELBASE = 4;
const double DT_CONTROLS = 0.025;

const double VELOCITY_NOISE = 0.3;
const double TURNING_NOISE = 3.0*M_PI/180.;

const double MAX_RANGE = 10.0;
const double DT_OBSERVE = 8*DT_CONTROLS;

const double OBS_DIST_NOISE = 0.1;
const double OBS_ANGLE_NOISE = 1.0*M_PI/180.;
}


double *inputVars, *vars;
std::vector<int> maskIndices;

Matrix<X_DIM> dynfunc(const Matrix<X_DIM>& x, const Matrix<U_DIM>& u, const Matrix<Q_DIM>& q)
{
	Matrix<X_DIM> xAdd = zeros<X_DIM,1>();
  
	xAdd[0] = (u[0]+q[0]) * DT * cos(x[2]+u[1]+q[0]);
	xAdd[1] = (u[0]+q[0]) * DT * sin(x[2]+u[1]+q[1]);
	xAdd[2] = (u[0]+q[0]) * DT * sin(u[1]+q[1])/config::WHEELBASE;

	Matrix<X_DIM> xNew = x + xAdd;
    return xNew;
}


Matrix<Z_DIM> obsfunc(const Matrix<X_DIM>& x, const Matrix<R_DIM>& r)
{
  //Matrix<L_DIM,1> landmarks = x.subMatrix<L_DIM,1>(P_DIM, 0);
	double xPos = x[0], yPos = x[1], angle = x[2];
	double l0, l1;
	double dx, dy;

	Matrix<Z_DIM> obs = zeros<Z_DIM,1>();

	for(int i = 0; i < L_DIM; i += 2) {
		dx = l0 - xPos;
		dy = l1 - yPos;

		if ((fabs(dx) < config::MAX_RANGE) &&
			(fabs(dy) < config::MAX_RANGE) &&
			(dx*cos(angle) + dy*sin(angle) > 0) &&
			(dx*dx + dy*dy < config::MAX_RANGE*config::MAX_RANGE))
		{
		  // TODO: is noise being added ocrrectly? (check add_observation_noise.m)
		  // This is correct
			obs[i] = sqrt(dx*dx + dy*dy) + r[i];
			obs[i+1] = atan2(dy, dx) - angle + r[i+1];
		}
	}

	return obs;
}


// Jacobians: df(x,u,q)/dx, df(x,u,q)/dq
void linearizeDynamics(const Matrix<X_DIM>& x, const Matrix<U_DIM>& u, const Matrix<Q_DIM>& q, Matrix<3,3>& A, Matrix<3,2>& M)
{
  //g is control input steer angle
  double s= sin(u[1]+x[2]); double c= cos(u[1]+x[2]);
  double vts= u[0]*DT*s; double vtc= u[0]*DT*c;

  M.reset();
  M(0, 0) = DT*c;
  M(0, 1) = -vts;
  M(1, 0) = DT*s;
  M(1, 1) = vtc;
  M(2, 0) = DT*sin(u[1])/length;
  M(2, 1) = u[0]*DT*cos(u[1])/length;

  A.reset();
  A(0,0) = 1;
  A(1,1) = 1;
  A(2,2) = 1;
  A(0,2) = -vts;
  A(1,2) = vtc;
  /*
	A.reset();
	Matrix<X_DIM> xr(x), xl(x);
	for (size_t i = 0; i < X_DIM; ++i) {
		xr[i] += step; xl[i] -= step;
		A.insert(0,i, (dynfunc(xr, u, q) - dynfunc(xl, u, q)) / (xr[i] - xl[i]));
		xr[i] = x[i]; xl[i] = x[i];
	}

	M.reset();
	Matrix<Q_DIM> qr(q), ql(q);
	for (size_t i = 0; i < Q_DIM; ++i) {
		qr[i] += step; ql[i] -= step;
		M.insert(0,i, (dynfunc(x, u, qr) - dynfunc(x, u, ql)) / (qr[i] - ql[i]));
		qr[i] = q[i]; ql[i] = q[i];
	}
  */
}

// Jacobians: df(x,u,q)/dx, df(x,u,q)/du
void linearizeDynamicsFunc(const Matrix<X_DIM>& x, const Matrix<U_DIM>& u, Matrix<X_DIM,X_DIM>& A, Matrix<X_DIM,Q_DIM>& M, const Matrix<X_DIM>& h)
{
	A.reset();
	Matrix<X_DIM> xr(x), xl(x);
	for (size_t i = 0; i < X_DIM; ++i) {
		xr[i] += step; xl[i] -= step;
		A.insert(0,i, (dynfunc(xr, u, zeros<Q_DIM,1>()) - dynfunc(xl, u,  zeros<Q_DIM,1>())) / (xr[i] - xl[i]));
		xr[i] = x[i]; xl[i] = x[i];
	}

	M.reset();
	Matrix<U_DIM> ur(u), ul(u);
	for (size_t i = 0; i < U_DIM; ++i) {
		ur[i] += step; ul[i] -= step;
		M.insert(0,i, (dynfunc(x, ur,  zeros<Q_DIM,1>()) - dynfunc(x, ul,  zeros<Q_DIM,1>())) / (ur[i] - ul[i]));
		ur[i] = u[i]; ul[i] = u[i];
	}

	h = dynfunc(x, u,  zeros<Q_DIM,1>());
}


// Jacobians: dh(x,r)/dx, dh(x,r)/dr
void linearizeObservation(const Matrix<X_DIM>& x, const Matrix<R_DIM>& r, Matrix<Z_DIM,X_DIM>& H, Matrix<Z_DIM,R_DIM>& N)
{
  //Approximate only setting H for 
  H.reset();
  for (int i=0; i<NUM_LANDMARKS; ++i) {
    double dx = x[3+2*i] - x[0];
    double dy = x[4+2*i] - x[1];
    double d2 = std::pow(dx,2) + std::pow(dy, 2);
    double d = std::sqrt(d2);
    double xd = dx/d;
    double yd = dy/d;
    double xd2 = dx/d2;
    double yd2 = dy/d2;
  
    //Approximate H being 0 for observations out of range
    //Currently full circle, should only be half circle
    double range_scale = (1+exp(10*(d-config::MAX_RANGE)));

    H(i, 0) = -xd / range_scale;
    H(i, 1) = -yd / range_scale;
    H(i, 2) = 0;
    H(i+1, 0) = yd2 / range_scale;
    H(i+1, 1) = -xd2 / range_scale;
    H(i+1, 2) = -1 / range_scale;
    H(i, 3+i) = xd / range_scale;
    H(i, 3+i+1) = yd / range_scale;
    H(i+1, 3+i) = -yd2 / range_scale;
    H(i+1, 3+i+1) = xd2 / range_scale;
  }
  
  /*
	H.reset();
	Matrix<X_DIM> xr(x), xl(x);
	for (size_t i = 0; i < X_DIM; ++i) {
	  xr[i] += step; xl[i] -= step;
	  H.insert(0,i, (obsfunc(xr, r) - obsfunc(xl, r)) / (xr[i] - xl[i]));
	  xr[i] = x[i]; xl[i] = x[i];
	}
  
	N.reset();
	Matrix<R_DIM> rr(r), rl(r);
	for (size_t i = 0; i < R_DIM; ++i) {
		rr[i] += step; rl[i] -= step;
		N.insert(0,i, (obsfunc(x, rr) - obsfunc(x, rl)) / (rr[i] - rl[i]));
		rr[i] = r[i]; rl[i] = r[i];

	}
  */
}

// Switch between belief vector and matrices
void unVec(const Matrix<B_DIM>& b, Matrix<X_DIM>& x, Matrix<X_DIM,X_DIM>& SqrtSigma) {
	x = b.subMatrix<X_DIM,1>(0,0);
	size_t idx = X_DIM;
	for (size_t j = 0; j < X_DIM; ++j) {
		for (size_t i = j; i < X_DIM; ++i) {
			SqrtSigma(i,j) = b[idx];
			SqrtSigma(j,i) = b[idx];
			++idx;
		}
	}
}

void vec(const Matrix<X_DIM>& x, const Matrix<X_DIM,X_DIM>& SqrtSigma, Matrix<B_DIM>& b) {
	b.insert(0,0,x);
	size_t idx = X_DIM;
	for (size_t j = 0; j < X_DIM; ++j) {
		for (size_t i = j; i < X_DIM; ++i) {
			b[idx] = 0.5 * (SqrtSigma(i,j) + SqrtSigma(j,i));
			++idx;
		}
	}
}

// Belief dynamics
Matrix<B_DIM> beliefDynamics(const Matrix<B_DIM>& b, const Matrix<U_DIM>& u) {
  //Need to approximate not recieving observations for landmarks out of range
  //Is it necessary to approximate augment somehow?
	Matrix<X_DIM> x;
	Matrix<X_DIM,X_DIM> SqrtSigma;
	unVec(b, x, SqrtSigma);

	Matrix<X_DIM,X_DIM> Sigma = SqrtSigma*SqrtSigma;

	Matrix<3,3> A;
	Matrix<3,2> M;
	linearizeDynamics(x, u, zeros<Q_DIM,1>(), A, M);
	Sigma.insert<3,3>(0,0, A*Sigma.subMatrix<3,3>(0,0)*~A + M*Q*~M);

	Sigma.insert<3,X_DIM-3>(0,3, A*(Sigma.subMatrix<3,X_DIM-3>(0,3)));
	x = dynfunc(x, u, zeros<Q_DIM,1>());

	//Should include a Q here
	//	Sigma = A*Sigma*~A + M*Q*~M;
	Matrix<Z_DIM,X_DIM> H;
	Matrix<Z_DIM,R_DIM> N;
	linearizeObservation(x, zeros<R_DIM,1>(), H, N);
	//Should include an R here
	Matrix<X_DIM,Z_DIM> K = Sigma*~H/(H*Sigma*~H + R);//N*R*~N);

	Sigma = (identity<X_DIM>() - K*H)*Sigma;

	Matrix<B_DIM> g;
	vec(x, sqrtm(Sigma), g);

	return g;
}


void setupDstarInterface(std::string mask) {
	std::stringstream ss(mask);
	int val, i=0;
	while (ss >> val) {
		if (val == 1) {
			maskIndices.push_back(i);
		}
		i++;
	}

	inputVars = new double[i];
}

void pythonDisplayTrajectory(std::vector< Matrix<B_DIM> >& B, std::vector< Matrix<P_DIM> >& waypoints, int time_steps)
{

	// B_vec is only for the robot, not the landmarks
	py::list B_vec;
	for(int j=0; j < P_DIM; j++) {
		for(int i=0; i < time_steps; i++) {
			B_vec.append(B[i][j]);
		}
	}

	Matrix<X_DIM> x;
	Matrix<X_DIM,X_DIM> SqrtSigma;
	for(int i=0; i < time_steps; i++) {
		unVec(B[i], x, SqrtSigma);
		B_vec.append(SqrtSigma[0,0]);
	}
	for(int i=0; i < time_steps; i++) {
		unVec(B[i], x, SqrtSigma);
		B_vec.append(SqrtSigma[1,1]);
	}
	for(int i=0; i < time_steps; i++) {
		unVec(B[i], x, SqrtSigma);
		B_vec.append(SqrtSigma[1,0]);
	}

	py::list waypoints_vec;
	for(int j=0; j < 2; j++) {
		for(int i=0; i < NUM_WAYPOINTS; i++) {
			waypoints_vec.append(waypoints[i][j]);
		}
	}

	std::string workingDir = boost::filesystem::current_path().normalize().string();

	try
	{
		Py_Initialize();
		py::object main_module = py::import("__main__");
		py::object main_namespace = main_module.attr("__dict__");
		py::exec("import sys, os", main_namespace);
		py::exec(py::str("sys.path.append('"+workingDir+"/point-slam')"), main_namespace);
		py::object plot_mod = py::import("plot_point_slam");
		py::object plot_traj = plot_mod.attr("plot_point_trajectory");

		plot_traj(B_vec, waypoints_vec, time_steps);
	}
	catch(py::error_already_set const &)
	{
		PyErr_Print();
	}

}



void pythonDisplayTrajectory2(std::vector< Matrix<B_DIM> >& B)
{
	py::list Bvec;
	for(int j=0; j < B_DIM; j++) {
		for(int i=0; i < T; i++) {
			Bvec.append(B[i][j]);
		}
	}

	/*
	py::list x0_list, xGoal_list;
	for(int i=0; i < X_DIM; i++) {
		x0_list.append(x0[i]);
		xGoal_list.append(xGoal[i]);
	}
	*/

	//std::string workingDir = boost::filesystem::current_path().normalize().string();
	//std::string bspDir = workingDir.substr(0,workingDir.find("bsp"));

	try
	{
		Py_Initialize();
		py::object main_module = py::import("__main__");
		py::object main_namespace = main_module.attr("__dict__");
		//py::exec("import sys, os", main_namespace);
		//py::exec(py::str("sys.path.append('"+bspDir+"bsp/python')"), main_namespace);
		//py::exec("from bsp_light_dark import LightDarkModel", main_namespace);
		//py::object model = py::eval("LightDarkModel()", main_namespace);
		py::object plot_mod = py::import("plot_point_slam");
		py::object plot_traj = plot_mod.attr("plot_belief_trajectory");

		plot_traj(Bvec, B_DIM, X_DIM, T, camera_range, camera_view_angle);
	}
	catch(py::error_already_set const &)
	{
		PyErr_Print();
	}

}

#endif
