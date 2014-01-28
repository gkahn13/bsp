#ifndef __POINT_H__
#define __POINT_H__

#include <fstream>
#include <math.h>

#include "../util/matrix.h"
//extern "C" {
#include "../util/utils.h"
//}
//#include "util/Timer.h"
#include "../util/logging.h"

#include <Python.h>
//#include <pythonrun.h>
#include <boost/python.hpp>
#include <boost/filesystem.hpp>
#include "landmarks.h"

namespace py = boost::python;

#define TIMESTEPS 30
#define DT 1.0
#define NUM_LANDMARKS 3
#define X_DIM 3+NUM_LANDMARKS*2
#define U_DIM 2
#define Z_DIM 2*NUM_LANDMARKS
#define Q_DIM 2
#define R_DIM 2*NUM_LANDMARKS

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

double *inputVars, *vars;
std::vector<int> maskIndices;

Matrix<X_DIM> dynfunc(const Matrix<X_DIM>& x, const Matrix<U_DIM>& u, const Matrix<Q_DIM>& q)
{
  Matrix<X_DIM> xAdd = zeros<X_DIM,1>();
  
  xAdd[0] = (u[0]+q[0]) * DT * cos(x[2]+u[1]+q[0]);//+ noise 
  xAdd[1] = (u[0]+q[0]) * DT * sin(x[2]+u[1]+q[1]);//+ noise 
  xAdd[2] = (u[0]+q[0]) * DT * sin(u[1]+q[1])/length;//+ noise 
  Matrix<X_DIM> xNew = x + xAdd;
  return xNew;
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
  double t = (~(pt - p1) * (p2 - p1)/len)[0];
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


// Observation model
Matrix<Z_DIM> obsfunc(const Matrix<X_DIM>& x, const Matrix<R_DIM>& r)
{ 
  Matrix<Z_DIM> z;
  int idz=0;
  for (size_t i =0; i < NUM_LANDMARKS; ++i) {

    std::vector<Matrix<2>> camera_region(3);
    //Points form a triangle with height = sensing_range
    //First point is on the robots center position
    camera_region[0][0] = x[0];
    camera_region[0][1] = x[1];

    double side_length = camera_range / cos(camera_view_angle);

    camera_region[1][0] = x[0] + side_length*cos(x[2]-camera_view_angle);
    camera_region[1][1] = x[1] + side_length*sin(x[2]-camera_view_angle);

    camera_region[2][0] = x[0] + side_length*cos(x[2]+camera_view_angle);
    camera_region[2][1] = x[1] + side_length*sin(x[2]+camera_view_angle);

    Matrix<2> landmarkPoint = x.subMatrix<2>(3+2*i,0);
    double sd = signedDist(landmarkPoint, camera_region);

    //double dist = std::sqrt(std::pow(landmarkPoint[0]-x[0], 2) + std::pow(landmarkPoint[1]-x[1],2));
    //double ang = atan2(landmarkPoint[1]-x[1], landmarkPoint[0]-x[0]);
    
    //z[2*i] = dist + r[2*i]/std::pow(1+(exp(-sd)),20);
    //z[2*i+1] = ang - x[2] + r[2*i+1]/ std::pow(1+(exp(-sd)),20);

    //jferguson -update this
    //New observation model -> observe a distance and an angle
    //Dynamics model stays more or less the same

    z[2*i] = (x[0] - x[3+2*i]) + r[2*i]/ std::pow(1+(exp(-sd)),20);
    z[2*i+1] = (x[1] - x[4+2*i]) + r[2*i+1]/ std::pow(1+(exp(-sd)),20);
  }
  return z;
}


// Jacobians: df(x,u,q)/dx, df(x,u,q)/dq
void linearizeDynamics(const Matrix<X_DIM>& x, const Matrix<U_DIM>& u, const Matrix<Q_DIM>& q, Matrix<3,3>& A, Matrix<3,2>& M)
{
  //jferguson -- hardcode this
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
  //jferguson -- hardcode this
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

  H.reset();
  for (int i=0; i<NUM_LANDMARKS; ++i) {
    double dx = x[3+2*i] - x[0];
    double dy = x[4+2*i] - x[1];
    double d2 = std::pow(dx,2) + std::pow(dy, 2);
    double d = std::sqrt(d2);
    double xd = dx/d;
    double yd = dy/d;
    double xd2 = dx/d2;
    double yd2 = dy/2;
  
    H(i, 0) = -xd;
    H(i, 1) = -yd;
    H(i, 2) = 0;
    H(i+1, 0) = yd2;
    H(i+1, 1) = -xd2;
    H(i+1, 2) = -1;
    H(i, 3+i) = xd;
    H(i, 3+i+1) = yd;
    H(i+1, 3+i) = -yd2;
    H(i+1, 3+i+1) = xd2;
  }
  
  //jferguson -- hardcode this
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

void pythonDisplayTrajectory(std::vector< Matrix<B_DIM> >& B)
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
		py::object plot_mod = py::import("plot_slam");
		py::object plot_traj = plot_mod.attr("plot_belief_trajectory");

		plot_traj(Bvec, B_DIM, X_DIM, T, camera_range, camera_view_angle);
	}
	catch(py::error_already_set const &)
	{
		PyErr_Print();
	}

}

#endif
