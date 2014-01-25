#ifndef __POINT_H__
#define __POINT_H__

#include <fstream>
#include <math.h>

#include "util/matrix.h"
//extern "C" {
#include "util/utils.h"
//}
//#include "util/Timer.h"
#include "util/logging.h"

#include <Python.h>
//#include <pythonrun.h>
#include <boost/python.hpp>
#include <boost/filesystem.hpp>
#include "landmarks.h"

namespace py = boost::python;

#define TIMESTEPS 15
#define DT 1.0
#define NUM_LANDMARKS 3
#define X_DIM 3+NUM_LANDMARKS*2
#define U_DIM 2
#define Z_DIM 2*NUM_LANDMARKS
#define Q_DIM 3
#define R_DIM 2*NUM_LANDMARKS

#define S_DIM (((X_DIM+1)*X_DIM)/2)
#define B_DIM (X_DIM+S_DIM)

const double length = 1;
const double camera_range = 500;
const double camera_view_angle = pi/4;

const double step = 0.0078125*0.0078125;


Matrix<X_DIM> x0;
Matrix<X_DIM,X_DIM> SqrtSigma0;
Matrix<X_DIM> xGoal;
Matrix<X_DIM> xMin, xMax;
Matrix<U_DIM> uMin, uMax;

//These are necessary for finding signed distance to approximate sensing regions
std::vector<Matrix<POS_DIM>> actualLandmarks;

Matrix<4,4> cameraParams;

const int T = TIMESTEPS;
const double INFTY = 1e10;

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

const double alpha_belief = 10, alpha_final_belief = 10, alpha_control = 1, alpha_goal_state = 1;

double *inputVars, *vars;
std::vector<int> maskIndices;


Matrix<X_DIM> dynfunc(const Matrix<X_DIM>& x, const Matrix<U_DIM>& u, const Matrix<Q_DIM>& q)
{
  Matrix<X_DIM> xAdd = zeros<X_DIM,1>();
    
  xAdd[0] = u[0] * DT * cos(x[2]) + q[0];//+ noise 
  xAdd[1] = u[0] * DT * sin(x[2]) + q[1];//+ noise 
  xAdd[2] = u[0] * DT * tan(u[1])/length + q[2];//+ noise 

  Matrix<X_DIM> xNew = x + xAdd;

  return xNew;
}

// Observation model
Matrix<Z_DIM> obsfunc(const Matrix<X_DIM>& x, const Matrix<R_DIM>& r)
{ 
  Matrix<Z_DIM> z;
  int idz=0;
  for (size_t i =0; i < NUM_LANDMARKS; ++i) {

    std::vector<Matrix<2>> camera_region(3);
    //Points form a triangle with height = sensing_range
    //First point is on the robots center position
    camera_region[0] = x<2, 1>(0,0);
    double side_length = camera_range / cos(camera_view_angle);

    camera_region[1][0] = x[0] + side_length*cos(x[2]+camera_view_angle);
    camera_region[1][1] = x[1] + side_length*sin(x[2]+camera_view_angle);

    camera_region[2][0] = x[0] + side_length*cos(x[2]-camera_view_angle);
    camera_region[2][1] = x[1] + side_length*sin(x[2]-camera_view_angle);

    double sd = signedDist(actualLandmarks[i], camera_region);
    z[2*i] = (actualLandmarks[i][0] + r[2*i])/ std::pow((1+exp(-sd)),2);
    z[2*i+1] = (actualLandmarks[i][1] + r[2*i+1])/ std::pow((1+exp(-sd)),2);
  }
  return z;
}

double signedDist(const Matrix<2> pt, const std::vector<Matrix<2>> camera_region)
{
  //camera_region = shrink_camera(camera_region);
  double d12 = line_point_signed_dist(camera_region[0],camera_region[1],pt);
  double d23 = line_point_signed_dist(camera_region[1],camera_region[2],pt);
  double d31 = line_point_signed_dist(camera_region[2],camera_region[0],pt);
  if (d12 <= 0.0 && d23 <= 0.0 && d31 <= 0.0) {
    return std::max({d12, d23, d31});
  } else {
    d12 = segment_point_dist(camera_region[0],camera_region[1], pt);
    d23 = segment_point_dist(camera_region[1],camera_region[2], pt);
    d31 = segment_point_dist(camera_region[2],camera_region[0], pt);
    return std::min({d12, d23, d31});
  }
}

double line_point_signed_dist(const Matrix<2> p1, const Matrix<2> p2, const Matrix<2> pt)
{
  double a = p2[1] - p1[1];
  double b = p1[0] - p2[0];
  double c = p2[0] * p1[1] - p1[0] * p2[1];
  return (a*pt[0] + b*pt[1] + c)/std::sqrt(a*a + b*b);
}

double segment_point_dist(const Matrix<2> p1, const Matrix<2> p2, const Matrix<2> pt)
{
  double len = std::pow((p2[0]-p1[0]), 2) + std::pow((p2[1] - p1[1]), 2);
  double t = ~(pt - p1) * (p2 - p1)/len;
  t = std::max(t,0.0);
  t = std::min(t,1.0);

  Matrix<2> pn = p1 + t*(p2-p1);
  double result = std::pow(pn[0] - pt[0], 2);
  result += std::pow(pn[1] - pt[1], 2);
  return std::sqrt(result);
}
/*
function T = shrink_camera(T)
n13 = [T(2,1) - T(2,3); T(1,3) - T(1,1)];
n32 = [T(2,3) - T(2,2); T(1,2) - T(1,3)];
n21 = [T(2,2) - T(2,1); T(1,1) - T(1,2)];
T(:,1) = T(:,1) + 0.75*n32/norm(n32);
T(:,2) = T(:,2) + 0.75*n13/norm(n13);
T(:,3) = T(:,3) + 0.75*n21/norm(n21);
end
*/

// Jacobians: df(x,u,q)/dx, df(x,u,q)/dq
void linearizeDynamics(const Matrix<X_DIM>& x, const Matrix<U_DIM>& u, const Matrix<Q_DIM>& q, Matrix<X_DIM,X_DIM>& A, Matrix<X_DIM,Q_DIM>& M)
{
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
}

// Jacobians: df(x,u,q)/dx, df(x,u,q)/du
void linearizeDynamicsFunc(const Matrix<X_DIM>& x, const Matrix<U_DIM>& u, Matrix<X_DIM,X_DIM>& A, Matrix<X_DIM,Q_DIM>& M, const Matrix<X_DIM>& h)
{
  
	A.reset();
	Matrix<X_DIM> xr(x), xl(x);
	for (size_t i = 0; i < X_DIM; ++i) {
		xr[i] += step; xl[i] -= step;
		A.insert(0,i, (dynfunc(xr, u, q) - dynfunc(xl, u, q)) / (xr[i] - xl[i]));
		xr[i] = x[i]; xl[i] = x[i];
	}

	M.reset();
	Matrix<U_DIM> ur(u), ul(u);
	for (size_t i = 0; i < U_DIM; ++i) {
		ur[i] += step; ul[i] -= step;
		M.insert(0,i, (dynfunc(x, ur, 0) - dynfunc(x, ul, 0)) / (ur[i] - ul[i]));
		ur[i] = u[i]; ul[i] = u[i];
	}

	h = dynfunc(x, u, 0);
}


// Jacobians: dh(x,r)/dx, dh(x,r)/dr
void linearizeObservation(const Matrix<X_DIM>& x, const Matrix<R_DIM>& r, Matrix<Z_DIM,X_DIM>& H, Matrix<Z_DIM,R_DIM>& N)
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
			++idx;
		}
	}
}

// Belief dynamics
Matrix<B_DIM> beliefDynamics(const Matrix<B_DIM>& b, const Matrix<U_DIM>& u) {
	Matrix<X_DIM> x;
	Matrix<X_DIM,X_DIM> SqrtSigma;
	unVec(b, x, SqrtSigma);

	Matrix<X_DIM,X_DIM> Sigma = SqrtSigma*SqrtSigma;

	Matrix<X_DIM,X_DIM> A;
	Matrix<X_DIM,Q_DIM> M;
	linearizeDynamics(x, u, zeros<Q_DIM,1>(), A, M);

	x = dynfunc(x, u, zeros<Q_DIM,1>());
	Sigma = A*Sigma*~A + M*~M;

	Matrix<Z_DIM,X_DIM> H;
	Matrix<Z_DIM,R_DIM> N;
	linearizeObservation(x, zeros<R_DIM,1>(), H, N);

	Matrix<X_DIM,Z_DIM> K = Sigma*~H/(H*Sigma*~H + N*~N);

	Sigma = (identity<X_DIM>() - K*H)*Sigma;

	Matrix<B_DIM> g;
	vec(x, sqrt(Sigma), g);

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

void pythonDisplayTrajectory(std::vector< Matrix<B_DIM> >& B, std::vector< Matrix<U_DIM> >& U)
{
	for (int t = 0; t < T-1; ++t) {
		B[t+1] = beliefDynamics(B[t], U[t]);
	}

	py::list Bvec;
	for(int j=0; j < B_DIM; j++) {
		for(int i=0; i < T; i++) {
			Bvec.append(B[i][j]);
		}
	}

	py::list Uvec;
		for(int j=0; j < U_DIM; j++) {
			for(int i=0; i < T-1; i++) {
			Uvec.append(U[i][j]);
		}
	}

	py::list x0_list, xGoal_list;
	for(int i=0; i < X_DIM; i++) {
		x0_list.append(x0[i]);
		xGoal_list.append(xGoal[i]);
	}

	std::string workingDir = boost::filesystem::current_path().normalize().string();
	std::string bspDir = workingDir.substr(0,workingDir.find("bsp"));

	try
	{
		Py_Initialize();
		py::object main_module = py::import("__main__");
		py::object main_namespace = main_module.attr("__dict__");
		py::exec("import sys, os", main_namespace);
		py::exec(py::str("sys.path.append('"+bspDir+"bsp/python')"), main_namespace);
		py::exec("from bsp_light_dark import LightDarkModel", main_namespace);
		py::object model = py::eval("LightDarkModel()", main_namespace);
		py::object plot_mod = py::import("plot");
		py::object plot_traj = plot_mod.attr("plot_belief_trajectory_cpp");

		plot_traj(Bvec, Uvec, model, x0_list, xGoal_list, T);
	}
	catch(py::error_already_set const &)
	{
		PyErr_Print();
	}

}

#endif
