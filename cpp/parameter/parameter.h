#ifndef __POINT_H__
#define __POINT_H__

#include <fstream>

#include "util/matrix.h"


#include "util/utils.h"
#include "util/logging.h"

#include <Python.h>
//#include <pythonrun.h>
#include <boost/python.hpp>
#include <boost/filesystem.hpp>

namespace py = boost::python;

#define TIMESTEPS 15
#define DT 1.0

// J_DIM == SIM_X_DIM in original file
#define J_DIM 4 // number of joints in state (2 position and 2 velocity)
#define K_DIM 4 // number of parameters in state (2 masses and 2 lengths)

#define X_DIM 8
#define U_DIM 2
#define Z_DIM 4
#define Q_DIM 2
#define R_DIM 4

#define S_DIM (((X_DIM+1)*X_DIM)/2)
#define B_DIM (X_DIM+S_DIM)
#define XU_DIM (X_DIM*T+U_DIM*(T-1))

namespace dynamics {
const double gravity = 9.86;

const double mass1 = 0.092;
const double mass_center1 = 0.062;
const double length1 = 0.15;
const double damping1 = 0.015;
const double coulomb1 = 0.006;
const double rotational1_inertia = 0.00064;
const double motor1_inertia = 0.00000065;
const double motor1_torque_const = 0.0077;
const double gear1_ratio = 70;

const double mass2 = 0.077;
const double mass_center2 = 0.036;
const double length2 = 0.15;
const double damping2 = 0.015;
const double coulomb2 = 0.006;
const double rotational2_inertia = 0.0003;
const double motor2_inertia = 0.00000065;
const double motor2_torque_const = 0.0077;
const double gear2_ratio = 70;
}

const double step = 0.0078125 / 16;


const int T = TIMESTEPS;
const double INFTY = 1e10;

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

const double alpha_belief = 10, alpha_final_belief = 10, alpha_control = 1, alpha_goal_state = 1;

double *inputVars, *vars;
std::vector<int> maskIndices;


Matrix<J_DIM> jointdynfunc(const Matrix<J_DIM>& j, const Matrix<U_DIM>& u, double length1, double length2, double mass1, double mass2)
{
	Matrix<J_DIM> jNew;
	double j0 = j[0], j1 = j[1], j2 = j[2], j3 = j[3];

	//double mass2_sq = mass2*mass2;
	//double mass_center2_sq = dynamics::mass_center2*dynamics::mass_center2;
	double gear1_ratio_sq = dynamics::gear1_ratio*dynamics::gear1_ratio;
	double gear2_ratio_sq = dynamics::gear2_ratio*dynamics::gear2_ratio;
	double length1_sq = length1*length1;
	//double length2_sq = length2*length2;
	double cos2minus1 = cos(j1 - j0);
	//double cos2minus1_sq = cos2minus1*cos2minus1;

	double unknown_h = length1*dynamics::mass_center2*mass2*sin(j1 - j0);

	double H11 = dynamics::motor1_inertia*gear1_ratio_sq + dynamics::rotational1_inertia + mass2*length1_sq;
	double H12 = length1*dynamics::mass_center2*mass2*cos2minus1;
	double H22 = dynamics::motor2_inertia*gear2_ratio_sq + dynamics::rotational2_inertia;

	double D = H11*H22 - H12*H12;

	// sgn --> copysign since coulmb1 and coulomb 2 same
	double v1 = dynamics::gear1_ratio*dynamics::motor1_torque_const*u[0]
				+ unknown_h*j3*j3
				- (dynamics::mass_center1*mass1 + length1*mass2)*dynamics::gravity*cos(j0)
				- (dynamics::damping1*j2 + copysign(dynamics::coulomb1, j2));
	double v2 = dynamics::gear2_ratio*dynamics::motor2_torque_const*u[1]
				- unknown_h*j2*j2
				- dynamics::mass_center2*mass2*dynamics::gravity*cos(j1)
				- (dynamics::damping2*j3 + copysign(dynamics::coulomb2, j3));

	jNew[0] = j2;
	jNew[1] = j3;
	jNew[2] = (H22*v1 - H12*v2)/D;
	jNew[3] = (-H12*v1 + H11*v2)/D;

	return jNew;
}

// for both joints and params
Matrix<X_DIM> dynfunc(const Matrix<X_DIM>& x, const Matrix<U_DIM>& u, const Matrix<U_DIM>& q)
{
	// RK4 integration
	Matrix<J_DIM> k1, k2, k3, k4, jinit;
	Matrix<K_DIM> params;

	jinit = x.subMatrix<J_DIM>(0,0);
	params = x.subMatrix<K_DIM>(J_DIM,0);

	double length1 = (params[0] == 0 ? 0.0 : 1/params[0]);
	double length2 = (params[1] == 0 ? 0.0 : 1/params[1]);
	double mass1 = (params[2] == 0 ? 0.0 : 1/params[2]);
	double mass2 = (params[3] == 0 ? 0.0 : 1/params[3]);

	k1 = jointdynfunc(jinit, u, length1, length2, mass1, mass2);
	k2 = jointdynfunc(jinit + 0.5*step*k1, u, length1, length2, mass1, mass2);
	k3 = jointdynfunc(jinit + 0.5*step*k2, u, length1, length2, mass1, mass2);
	k4 = jointdynfunc(jinit + step*k3, u, length1, length2, mass1, mass2);

	Matrix<X_DIM> xNew = zeros<X_DIM,1>();
	xNew.insert(0, 0, jinit + step*(k1 + 2.0*(k2 + k3) + k4)/6.0);
	xNew.insert(J_DIM, 0, params);

	//std::cout << "x: " << ~x;
	//std::cout << "xNew: " << ~xNew << std::endl;

	return xNew;
}

// Observation model
Matrix<Z_DIM> obsfunc(const Matrix<X_DIM>& x, const Matrix<R_DIM>& r)
{
	Matrix<Z_DIM> z;
	double x0, x1, x2, x3, cosx0, sinx0, cosx1, sinx1;

	x0 = x[0]; x1 = x[1]; x2 = x[2]; x3 = x[3];

	cosx0 = cos(x0);
	sinx0 = sin(x0);
	cosx1 = cos(x1);
	sinx1 = sin(x1);

	z[0] = dynamics::length1*cosx0 + dynamics::length2*cosx1;
	z[1] = dynamics::length1*sinx0 + dynamics::length2*sinx1;
	z[2] = x2;
	z[3] = x3;

	return z;
}

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
void unVec(const Matrix<B_DIM>& b, Matrix<X_DIM>& x, Matrix<X_DIM,X_DIM>& S) {
	x = b.subMatrix<X_DIM,1>(0,0);
	size_t idx = X_DIM;
	for (size_t j = 0; j < X_DIM; ++j) {
		for (size_t i = j; i < X_DIM; ++i) {
			S(i,j) = b[idx];
			S(j,i) = b[idx];
			++idx;
		}
	}
}

void vec(const Matrix<X_DIM>& x, const Matrix<X_DIM,X_DIM>& S, Matrix<B_DIM>& b, bool isSqrtSigma = true) {
	b.insert(0,0,x);
	size_t idx = X_DIM;
	for (size_t j = 0; j < X_DIM; ++j) {
		for (size_t i = j; i < X_DIM; ++i) {
			if (isSqrtSigma) {
				b[idx] = 0.5 * (S(i,j) + S(j,i));
			} else {
				b[idx] = S(i,j);
			}
			++idx;
		}
	}
}


// Belief dynamics
Matrix<B_DIM> beliefDynamics(const Matrix<B_DIM>& b, const Matrix<U_DIM>& u, bool isSqrtSigma = true) {
	Matrix<X_DIM> x;
	Matrix<X_DIM,X_DIM> Sigma;
	unVec(b, x, Sigma);

	if (isSqrtSigma) {
		Sigma = Sigma*Sigma;
	}

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
	vec(x, Sigma, g, isSqrtSigma);

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

void pythonDisplayTrajectory(std::vector< Matrix<B_DIM> >& B, std::vector< Matrix<U_DIM> >& U, Matrix<X_DIM> x0, Matrix<X_DIM> xGoal)
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
