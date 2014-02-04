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

#define TIMESTEPS 15
#define DT 1.0
#define NUM_LANDMARKS 1
#define NUM_WAYPOINTS 3

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


const double step = 0.0078125*0.0078125;


std::vector< Matrix<P_DIM> > waypoints(NUM_WAYPOINTS);
std::vector< Matrix<P_DIM> > landmarks(NUM_LANDMARKS);

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

const double alpha_belief = 10, alpha_final_belief = 100, alpha_control = .01, alpha_goal_state = 1;

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

const double OBS_DIST_NOISE = 1 * 0.1;
const double OBS_ANGLE_NOISE = 1 * 1.0*M_PI/180.;

const double ALPHA_OBS = .75;
}


void initProblemParams()
{
	Q = zeros<Q_DIM,Q_DIM>();
	Q(0,0) = 1*config::VELOCITY_NOISE*config::VELOCITY_NOISE; // 20
	Q(1,1) = 1*config::TURNING_NOISE*config::TURNING_NOISE; // 1e-2

	R = zeros<R_DIM, R_DIM>();
	for(int i = 0; i < R_DIM-1; i+=2) {
		R(i,i) = 1*config::OBS_DIST_NOISE*config::OBS_DIST_NOISE; // 1
		R(i+1,i+1) = 1*config::OBS_ANGLE_NOISE*config::OBS_ANGLE_NOISE; // 1e-5
	}

	waypoints[0][0] = 60; waypoints[0][1] = 0;
	waypoints[1][0] = 60; waypoints[1][1] = 20;
	waypoints[2][0] = 0; waypoints[2][1] = 20;

	landmarks[0][0] = 30; landmarks[0][1] = 12.5;
	//landmarks[0][0] = 0; landmarks[0][1] = 0;
	//landmarks[1][0] = 30; landmarks[1][1] = 0;
	//landmarks[2][0] = 60; landmarks[2][1] = 0;
	//landmarks[3][0] = 60; landmarks[3][1] = 20;
	//landmarks[4][0] = 30; landmarks[4][1] = 5;

	// start at (0, 0)
	// landmarks will be the same for all waypoint-waypoint optimizations
	x0.insert(0, 0, zeros<C_DIM,1>());
	for(int i = 0; i < NUM_LANDMARKS; ++i) {
		x0.insert(C_DIM+2*i, 0, landmarks[i]);
	}

	//This starts out at 0 for car, landmarks are set based on the car's sigma when first seen
	SqrtSigma0 = zeros<X_DIM, X_DIM>();
	for(int i = 0; i < C_DIM; ++i) { SqrtSigma0(i,i) = .1; }
	for(int i = 0; i < L_DIM; ++i) { SqrtSigma0(C_DIM+i,C_DIM+i) = 1; }

	// TODO: think of better values for these
	uMin[0] = -10;
	uMin[1] = -M_PI;
	uMax[0] = 10;
	uMax[1] = M_PI;

	xMin[0] = -20;
	xMin[1] = -20;
	xMin[2] = -M_PI;
	xMin[3] = landmarks[0][0] - 5;
	xMin[4] = landmarks[0][1] - 5;

	xMax[0] = 80;
	xMax[1] = 80;
	xMax[2] = M_PI;
	xMax[3] = landmarks[0][0] + 5;
	xMax[4] = landmarks[0][1] + 5;

}

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
	double dx, dy;

	Matrix<Z_DIM> obs = zeros<Z_DIM,1>();

	for(int i = 0; i < L_DIM; i += 2) {
		dx = x[C_DIM+i] - xPos;
		dy = x[C_DIM+i+1] - yPos;

		//if ((fabs(dx) < config::MAX_RANGE) &&
		//	(fabs(dy) < config::MAX_RANGE) &&
		//	(dx*cos(angle) + dy*sin(angle) > 0) &&
		//	(dx*dx + dy*dy < config::MAX_RANGE*config::MAX_RANGE))
		//{
			obs[i] = sqrt(dx*dx + dy*dy) + r[i];
			obs[i+1] = atan2(dy, dx) - angle + r[i+1];
		//}
	}

	return obs;
}

Matrix<Z_DIM,Z_DIM> deltaMatrix(const Matrix<X_DIM>& x) {
	Matrix<Z_DIM,Z_DIM> delta = zeros<Z_DIM,Z_DIM>();
	double l0, l1, dist;
	for(int i=C_DIM; i < X_DIM; i += 2) {
		l0 = x[i];
		l1 = x[i+1];

		dist = sqrt((x[0] - l0)*(x[0] - l0) + (x[1] - l1)*(x[1] - l1));

		double signed_dist = 1/(1+exp(-config::ALPHA_OBS*(config::MAX_RANGE-dist)));
		delta(i-C_DIM,i-C_DIM) = signed_dist;
		delta(i-C_DIM+1,i-C_DIM+1) = signed_dist;
	}

	return delta;
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
  M(2, 0) = DT*sin(u[1])/config::WHEELBASE;
  M(2, 1) = u[0]*DT*cos(u[1])/config::WHEELBASE;

  A.reset();
  A(0,0) = 1;
  A(1,1) = 1;
  A(2,2) = 1;
  A(0,2) = -vts;
  A(1,2) = vtc;

}

// Jacobians: df(x,u,q)/dx, df(x,u,q)/du
void linearizeDynamicsFiniteDiff(const Matrix<X_DIM>& x, const Matrix<U_DIM>& u, const Matrix<Q_DIM>& q, Matrix<X_DIM,X_DIM>& A, Matrix<X_DIM,Q_DIM>& M)
{
	A.reset();
	Matrix<X_DIM> xr(x), xl(x);
	for (size_t i = 0; i < X_DIM; ++i) {
		xr[i] += step; xl[i] -= step;
		A.insert(0,i, (dynfunc(xr, u, zeros<Q_DIM,1>()) - dynfunc(xl, u,  zeros<Q_DIM,1>())) / (xr[i] - xl[i]));
		xr[i] = x[i]; xl[i] = x[i];
	}

	M.reset();
	Matrix<Q_DIM> qr(q), ql(q);
	for (size_t i = 0; i < Q_DIM; ++i) {
		qr[i] += step; ql[i] -= step;
		M.insert(0,i, (dynfunc(x, u, qr) - dynfunc(x, u, ql)) / (qr[i] - ql[i]));
		qr[i] = q[i]; ql[i] = q[i];
	}
}


// Jacobians: dh(x,r)/dx, dh(x,r)/dr
void linearizeObservation(const Matrix<X_DIM>& x, const Matrix<R_DIM>& r, Matrix<Z_DIM,X_DIM>& H, Matrix<Z_DIM,R_DIM>& N)
{

  //Approximate only setting H for 
  H.reset();
  for (int i=0; i < L_DIM; i+=2) {
    double dx = x[C_DIM+i] - x[0];
    double dy = x[C_DIM+i+1] - x[1];
    double d2 = dx*dx + dy*dy + 1e-10;
    double d = sqrt(d2 + 1e-10);
    double xd = dx/d;
    double yd = dy/d;
    double xd2 = dx/d2;
    double yd2 = dy/d2;
  
    // Approximate H being 0 for observations out of range
    // Currently full circle, should only be half circle
    // double range_scale = (1+exp(10*(d-config::MAX_RANGE)));
    // GREG: trying delta matrix for now
    double range_scale = 1;

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
	*/
  
	N.reset();
	Matrix<R_DIM> rr(r), rl(r);
	for (size_t i = 0; i < R_DIM; ++i) {
		rr[i] += step; rl[i] -= step;
		N.insert(0,i, (obsfunc(x, rr) - obsfunc(x, rl)) / (rr[i] - rl[i]));
		rr[i] = r[i]; rl[i] = r[i];

	}

}

// Jacobians: dh(x,r)/dx, dh(x,r)/dr
void linearizeObservationFiniteDiff(const Matrix<X_DIM>& x, const Matrix<R_DIM>& r, Matrix<Z_DIM,X_DIM>& H, Matrix<Z_DIM,R_DIM>& N)
{

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
			//waypoints[4][0] = 0; waypoints[4][1] = 20;
				++idx;
		}
	}
}

// Belief dynamics
Matrix<B_DIM> beliefDynamics(const Matrix<B_DIM>& b, const Matrix<U_DIM>& u) {
  //Need to approximate not recieving observations for landmarks out of range
  //Is it necessary to approximate augment somehow?
	// we can try adding the binary matrix
	Matrix<X_DIM> x;
	Matrix<X_DIM,X_DIM> SqrtSigma;
	unVec(b, x, SqrtSigma);

	Matrix<X_DIM,X_DIM> Sigma = SqrtSigma*SqrtSigma;

	Matrix<C_DIM,C_DIM> Acar;
	Matrix<C_DIM,Q_DIM> Mcar;
	linearizeDynamics(x, u, zeros<Q_DIM,1>(), Acar, Mcar);

	Matrix<X_DIM,X_DIM> A = identity<X_DIM>();
	A.insert<C_DIM,C_DIM>(0, 0, Acar);
	Matrix<X_DIM,Q_DIM> M = zeros<X_DIM,Q_DIM>();
	M.insert<C_DIM, 2>(0, 0, Mcar);

	//Matrix<X_DIM,X_DIM> A = zeros<X_DIM,X_DIM>();
	//Matrix<X_DIM,Q_DIM> M = zeros<X_DIM,Q_DIM>();
	//linearizeDynamicsFiniteDiff(x, u, zeros<Q_DIM,1>(), A, M);

	Sigma = A*Sigma*~A + M*Q*~M;
	//Sigma.insert(0,C_DIM, Acar*(Sigma.subMatrix<C_DIM,X_DIM-C_DIM>(0,3)));

	x = dynfunc(x, u, zeros<Q_DIM,1>());

	Matrix<Z_DIM,X_DIM> H;
	Matrix<Z_DIM,R_DIM> N;
	linearizeObservation(x, zeros<R_DIM,1>(), H, N);
	//Should include an R here

	Matrix<Z_DIM,Z_DIM> delta = deltaMatrix(x);

	//std::cout << "A" << std::endl << A << std::endl;
	//std::cout << "M" << std::endl << M << std::endl;
	//std::cout << "Sigma" << std::endl << Sigma << std::endl;
	//std::cout << "x" << std::endl << x << std::endl;
	//std::cout << "H" << std::endl << H << std::endl;
	//std::cout << "N" << std::endl << N << std::endl;
	//std::cout << "deltaMatrix" << std::endl;
	//for(int i=0; i < Z_DIM; ++i) {
	//	std::cout << delta(i,i) << " ";
	//}
	//std::cout << std::endl;
	//std::cout << "Sigma*~H*delta" << std::endl << Sigma*~H*delta << std::endl;

	Matrix<X_DIM,Z_DIM> K = ((Sigma*~H*delta)/(delta*H*Sigma*~H*delta + R))*delta;//N*R*~N);

	Sigma = (identity<X_DIM>() - K*H)*Sigma;

	Matrix<B_DIM> g;
	vec(x, sqrtm(Sigma), g);

	//std::cout << "K" << std::endl << K << std::endl;
	//std::cout << "Sigma" << std::endl << Sigma << std::endl;

	return g;
}

// Jacobians: dg(b,u)/db, dg(b,u)/du
void linearizeBeliefDynamics(const Matrix<B_DIM>& b, const Matrix<U_DIM>& u, Matrix<B_DIM,B_DIM>& F, Matrix<B_DIM,U_DIM>& G, Matrix<B_DIM>& h)
{
	F.reset();
	Matrix<B_DIM> br(b), bl(b);
	for (size_t i = 0; i < B_DIM; ++i) {
		br[i] += step; bl[i] -= step;
		F.insert(0,i, (beliefDynamics(br, u) - beliefDynamics(bl, u)) / (br[i] - bl[i]));
		br[i] = b[i]; bl[i] = b[i];
	}

	G.reset();
	Matrix<U_DIM> ur(u), ul(u);
	for (size_t i = 0; i < U_DIM; ++i) {
		ur[i] += step; ul[i] -= step;
		G.insert(0,i, (beliefDynamics(b, ur) - beliefDynamics(b, ul)) / (ur[i] - ul[i]));
		ur[i] = u[i]; ul[i] = u[i];
	}

	h = beliefDynamics(b, u);
}


void pythonDisplayTrajectory(std::vector< Matrix<B_DIM> >& B, std::vector< Matrix<U_DIM> >& U, std::vector< Matrix<P_DIM> >& waypoints, std::vector< Matrix<P_DIM> >& landmarks, int time_steps)
{

	// B_vec is only for the robot, not the landmarks
	py::list B_vec;
	for(int j=0; j < P_DIM; j++) {
		for(int i=0; i < time_steps; i++) {
			B_vec.append(B[i][j]);
		}
	}

	py::list U_vec;
	for(int j=0; j < U_DIM; j++) {
		for(int i=0; i < time_steps-1; i++) {
			U_vec.append(U[i][j]);
		}
	}

	Matrix<X_DIM> x;
	Matrix<X_DIM,X_DIM> SqrtSigma;
	for(int i=0; i < time_steps; i++) {
		unVec(B[i], x, SqrtSigma);
		B_vec.append(SqrtSigma(0,0));
	}
	for(int i=0; i < time_steps; i++) {
		unVec(B[i], x, SqrtSigma);
		B_vec.append(SqrtSigma(1,0));
	}
	for(int i=0; i < time_steps; i++) {
		unVec(B[i], x, SqrtSigma);
		B_vec.append(SqrtSigma(1,1));
	}

	py::list waypoints_vec;
	for(int j=0; j < 2; j++) {
		for(int i=0; i < NUM_WAYPOINTS; i++) {
			waypoints_vec.append(waypoints[i][j]);
		}
	}

	py::list landmarks_vec;
	for(int j=0; j < 2; j++) {
		for(int i=0; i < NUM_LANDMARKS; i++) {
			landmarks_vec.append(landmarks[i][j]);
		}
	}

	std::string workingDir = boost::filesystem::current_path().normalize().string();

	try
	{
		Py_Initialize();
		py::object main_module = py::import("__main__");
		py::object main_namespace = main_module.attr("__dict__");
		py::exec("import sys, os", main_namespace);
		py::exec(py::str("sys.path.append('"+workingDir+"/slam')"), main_namespace);
		py::object plot_mod = py::import("plot_point_slam");
		py::object plot_traj = plot_mod.attr("plot_point_trajectory");

		plot_traj(B_vec, U_vec, waypoints_vec, landmarks_vec, config::MAX_RANGE, config::ALPHA_OBS, time_steps, DT);
	}
	catch(py::error_already_set const &)
	{
		PyErr_Print();
	}

}



#endif
