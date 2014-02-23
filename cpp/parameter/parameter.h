#ifndef __PARAMETER_H__
#define __PARAMETER_H__

#include <fstream>

#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include "util/matrix.h"


#include "util/utils.h"
#include "util/logging.h"

#include <Python.h>
//#include <pythonrun.h>
#include <boost/python.hpp>
#include <boost/filesystem.hpp>

namespace py = boost::python;

// horizon is total lifetime of planning
// timesteps is how far into future accounting for during MPC
#define HORIZON 500
#define TIMESTEPS 15
#define DT 1.0/5.0


// for ILQG
//#define TIMESTEPS 500
//#define DT 1.0/100.0

// J_DIM == SIM_X_DIM in original file
#define J_DIM 4 // number of joints in state (2 position and 2 velocity)
#define K_DIM 4 // number of parameters in state (2 masses and 2 lengths)

#define X_DIM 8
#define U_DIM 2
#define Z_DIM 4
#define Q_DIM 8
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


const double diffEps = 0.0078125 / 16;


const int T = TIMESTEPS;
const double INFTY = 1e10;

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

const double alpha_belief = 10, alpha_final_belief = 10, alpha_control = 0, alpha_goal_state = 1;
const double alpha_joint_belief = 0, alpha_param_belief = 10,
			  alpha_final_joint_belief = 0, alpha_final_param_belief = 10, alpha_goal_joint_state = 10, alpha_goal_param_state = 1;

double *inputVars, *vars;
std::vector<int> maskIndices;

boost::mt19937 rng; 
boost::normal_distribution<> nd(0.0, 1.0);
boost::variate_generator<boost::mt19937&, 
                           boost::normal_distribution<> > var_nor(rng, nd);

double sgn(double x) {
	double delta = 1e-10;
	if (x > delta) { return 1; }
	else if (x < delta) { return -1; }
	else { return 0; }
}


Matrix<J_DIM> jointdynfunc(const Matrix<J_DIM>& j, const Matrix<U_DIM>& u, double length1, double length2, double mass1, double mass2)
{

	Matrix<J_DIM> jNew;

	double j0 = j[0];
	double j1 = j[1];
	double j2 = j[2];
	double j3 = j[3];
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

	double v1 = dynamics::gear1_ratio*dynamics::motor1_torque_const*u[0]
				+ unknown_h*j3*j3
				- (dynamics::mass_center1*mass1 + length1*mass2)*dynamics::gravity*cos(j0)
				- (dynamics::damping1*j2 + dynamics::coulomb1*sgn(j2));
	double v2 = dynamics::gear2_ratio*dynamics::motor2_torque_const*u[1]
				- unknown_h*j2*j2
				- dynamics::mass_center2*mass2*dynamics::gravity*cos(j1)
				- (dynamics::damping2*j3 + dynamics::coulomb2*sgn(j3));

	jNew[0] = j2;
	jNew[1] = j3;
	jNew[2] = (H22*v1 - H12*v2)/D;
	jNew[3] = (-H12*v1 + H11*v2)/D;

	return jNew;

}

// for both joints and params
Matrix<X_DIM> dynfunc(const Matrix<X_DIM>& x, const Matrix<U_DIM>& u)
{
	// RK4 integration
	Matrix<J_DIM> k1, k2, k3, k4, jinit;

	jinit = x.subMatrix<J_DIM>(0,0);

	double x4 = x[4];
	double x5 = x[5];
	double x6 = x[6];
	double x7 = x[7];

	double length1, length2, mass1, mass2;
	length1 = (x4 == 0 ? 0.0 : 1/x4);
	length2 = (x5 == 0 ? 0.0 : 1/x5);
	mass1 = (x6 == 0 ? 0.0 : 1/x6);
	mass2 = (x7 == 0 ? 0.0 : 1/x7);

	k1 = jointdynfunc(jinit, u, length1, length2, mass1, mass2);
	k2 = jointdynfunc(jinit + 0.5*DT*k1, u, length1, length2, mass1, mass2);
	k3 = jointdynfunc(jinit + 0.5*DT*k2, u, length1, length2, mass1, mass2);
	k4 = jointdynfunc(jinit + DT*k3, u, length1, length2, mass1, mass2);

	Matrix<X_DIM> xNew = zeros<X_DIM,1>();
	xNew.insert(0, 0, jinit + DT*(k1 + 2.0*(k2 + k3) + k4)/6.0);
	xNew[4] = x4;
	xNew[5] = x5;
	xNew[6] = x6;
	xNew[7] = x7;

	return xNew;
}

// Observation model
Matrix<Z_DIM> obsfunc(const Matrix<X_DIM>& x)
{
	Matrix<Z_DIM> z;

	double x0 = x[0];
	double x1 = x[1];
	double x2 = x[2];
	double x3 = x[3];

	double length1, length2;
	length1 = (x[4] == 0 ? 0 : 1/x[4]);
	length2 = (x[5] == 0 ? 0 : 1/x[5]);

	double cosx0 = cos(x0);
	double sinx0 = sin(x0);
	double cosx1 = cos(x1);
	double sinx1 = sin(x1);

	z[0] = length1*cosx0 + length2*cosx1;
	z[1] = length1*sinx0 + length2*sinx1;
	z[2] = x2;
	z[3] = x3;

	return z;
}

// for M = df(x,u,q)/dq
// this returns M*~M
inline Matrix<X_DIM,X_DIM> getMMT(const Matrix<X_DIM>& x, const Matrix<U_DIM>& u)
{
	Matrix<X_DIM,X_DIM> S = identity<X_DIM>();
	S(0,0) = 0.25*0.001 + 0.0000000000001;
	S(1,1) = 0.25*0.001 + 0.0000000000001;
	S(2,2) = 1.0*0.001 + 0.0000000000001;
	S(3,3) = 1.0*0.001 + 0.0000000000001;
	S(4,4) = 0.0005;
	S(5,5) = 0.0005;
	S(6,6) = 0.0001;
	S(7,7) = 0.0001;
	return S;
}


SymmetricMatrix<Q_DIM> varQ() {

	SymmetricMatrix<Q_DIM> S = identity<Q_DIM>();
	S(0,0) = 0.25*0.001;
	S(1,1) = 0.25*0.001;
	S(2,2) = 1.0*0.001;
	S(3,3) = 1.0*0.001;
	S(4,4) = 0.0005;
	S(5,5) = 0.0005;
	S(6,6) = 0.0001;
	S(7,7) = 0.0001;
	return S;
}

SymmetricMatrix<R_DIM> varR(){ 
	SymmetricMatrix<R_DIM> S = identity<R_DIM>();	
	S(0,0) = 0.0001;	
	S(1,1) = 0.0001;	
	S(2,2) = 0.00001;	
	S(3,3) = 0.00001;	
	return S;	\
}


// Jacobians: df(x,u,q)/dx
void linearizeDynamics(const Matrix<X_DIM>& x, const Matrix<U_DIM>& u, Matrix<X_DIM,X_DIM>& A)
{
	A.reset();
	Matrix<X_DIM> xr(x), xl(x);
	for (size_t i = 0; i < X_DIM; ++i) {
		xr[i] += diffEps; xl[i] -= diffEps;
		A.insert(0,i, (dynfunc(xr, u) - dynfunc(xl, u)) / (xr[i] - xl[i]));
		xr[i] = x[i]; xl[i] = x[i];
	}

}

// for N = dh(x,r)/dr
// this returns N*~N
inline Matrix<Z_DIM,Z_DIM> getNNT(const Matrix<X_DIM>& x)
{
	Matrix<Z_DIM,Z_DIM> S = identity<Z_DIM>();
	S(0,0) = 0.0001;
	S(1,1) = 0.0001;
	S(2,2) = 0.00001;
	S(3,3) = 0.00001;
	return S;
}

// Jacobians: dh(x,r)/dx
void linearizeObservation(const Matrix<X_DIM>& x, Matrix<Z_DIM,X_DIM>& H)
{
	H.reset();
	Matrix<X_DIM> xr(x), xl(x);
	for (size_t i = 0; i < X_DIM; ++i) {
		xr[i] += diffEps; xl[i] -= diffEps;
		H.insert(0,i, (obsfunc(xr) - obsfunc(xl)) / (xr[i] - xl[i]));
		xr[i] = x[i]; xl[i] = x[i];
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

void vec(const Matrix<X_DIM>& x, const Matrix<X_DIM,X_DIM>& S, Matrix<B_DIM>& b) {
	b.insert(0,0,x);
	size_t idx = X_DIM;
	for (size_t j = 0; j < X_DIM; ++j) {
		for (size_t i = j; i < X_DIM; ++i) {
			b[idx] = 0.5 * (S(i,j) + S(j,i));
			++idx;
		}
	}
}

Matrix<Q_DIM> qNoise(){

	Matrix<Q_DIM> q; 

	for(int i = 0; i<Q_DIM; i++){
		q[i] =  var_nor();
	}

	return q;
}


Matrix<R_DIM> rNoise(){

	Matrix<R_DIM> r; 

	

	for(int i = 0; i<R_DIM; i++){
		r[i] =  var_nor();
	}

	return r;
}

// Belief dynamics
Matrix<B_DIM> beliefDynamics(const Matrix<B_DIM>& b, const Matrix<U_DIM>& u) {
	Matrix<X_DIM> x;
	Matrix<X_DIM,X_DIM> Sigma;
	unVec(b, x, Sigma);

	Sigma = Sigma*Sigma;

	Matrix<X_DIM,X_DIM> A;
	linearizeDynamics(x, u, A);
	Matrix<X_DIM,X_DIM> MMT = getMMT(x, u);

	x = dynfunc(x, u);
	Sigma = A*Sigma*~A + MMT;

	Matrix<Z_DIM,X_DIM> H = zeros<Z_DIM,X_DIM>();
	linearizeObservation(x, H);
	Matrix<Z_DIM,R_DIM> NNT = getNNT(x);


	Matrix<X_DIM,Z_DIM> K = Sigma*~H/(H*Sigma*~H + NNT);

	Sigma = (identity<X_DIM>() - K*H)*Sigma;

	Matrix<B_DIM> g;
	vec(x, sqrtm(Sigma), g);

	return g;
}

// returns updated belief based on real current state, estimated current belief, and input
Matrix<B_DIM> executeControlStep(Matrix<X_DIM>& x_t_real, const Matrix<B_DIM>& b_t, const Matrix<U_DIM>& u_t) {
	// useRealParams = true
	// update real state (maximum likelihood)
	Matrix<R_DIM,R_DIM> Rchol;
	Matrix<Q_DIM,Q_DIM> Qchol; 

	chol(varR(),Rchol);
	chol(varQ(),Qchol); 


	Matrix<X_DIM> x_tp1_real = dynfunc(x_t_real, u_t)+ ~Qchol*qNoise();
	x_t_real = x_tp1_real; 
	// sense real state (maximum likelihood)
	Matrix<Z_DIM> z_tp1_real = obsfunc(x_tp1_real) + ~Rchol*rNoise();

	// now do EKF on belief and incorporate discrepancy
	// using the kalman gain

	Matrix<X_DIM> x_t;
	Matrix<X_DIM,X_DIM> SqrtSigma_t;
	unVec(b_t, x_t, SqrtSigma_t);

	Matrix<X_DIM,X_DIM> Sigma_t = SqrtSigma_t * SqrtSigma_t;

	// calculate df/dx and df/dm
	Matrix<X_DIM,X_DIM> A;
	linearizeDynamics(x_t, u_t, A);
	Matrix<X_DIM,X_DIM> MMT = getMMT(x_t, u_t);

	// propagate estimated state
	Matrix<X_DIM> x_tp1 = dynfunc(x_t, u_t);
	Matrix<X_DIM,X_DIM> Sigma_tp1 = A*Sigma_t*~A + MMT;

	// calculate dh/dx and dh/dn
	Matrix<Z_DIM,X_DIM> H = zeros<Z_DIM,X_DIM>();
	linearizeObservation(x_tp1, H);
	Matrix<Z_DIM,R_DIM> NNT = getNNT(x_tp1);

	// calculate Kalman gain for estimated state
	Matrix<X_DIM,Z_DIM> K = Sigma_tp1*~H/(H*Sigma_tp1*~H + NNT);

	// correct the new state using Kalman gain and the observation
	Matrix<X_DIM> x_tp1_adj = x_tp1 + K*(z_tp1_real - obsfunc(x_tp1));

	// correct the new covariance
	Matrix<X_DIM,X_DIM> W = ~(H*Sigma_tp1)*(((H*Sigma_tp1*~H) + NNT) % (H*Sigma_tp1));
	Matrix<X_DIM,X_DIM> Sigma_tp1_adj = Sigma_tp1 - W;

	Matrix<B_DIM> b_tp1_adj;
	vec(x_tp1_adj, sqrtm(Sigma_tp1_adj), b_tp1_adj);
	return b_tp1_adj;
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

void pythonDisplayTrajectory(std::vector< Matrix<U_DIM> >& U, Matrix<X_DIM,X_DIM> SqrtSigma0, Matrix<X_DIM> x0, Matrix<X_DIM> xGoal)
{
	std::vector< Matrix<B_DIM> > B(T);
	vec(x0, SqrtSigma0, B[0]);
	for (int t = 0; t < T - 1; ++t)
	{
		B[t+1] = beliefDynamics(B[t], U[t]);
	}

	Py_Initialize();

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

	std::string workingDir = boost::filesystem::current_path().normalize().string();

	py::object main_module = py::import("__main__");
	py::object main_namespace = main_module.attr("__dict__");
	py::exec("import sys, os", main_namespace);
	py::exec(py::str("sys.path.append('"+workingDir+"/parameter')"), main_namespace);
	py::object plot_mod = py::import("plot_parameter");
	py::object plot_traj = plot_mod.attr("plot_parameter_trajectory");

	plot_traj(Bvec, Uvec, B_DIM, X_DIM, U_DIM, T);

}



void pythonDisplayHistory(std::vector< Matrix<U_DIM> >& U,std::vector< Matrix<B_DIM> >& B, Matrix<X_DIM,X_DIM> SqrtSigma0, Matrix<X_DIM> x0, int H)
{
	
	Py_Initialize();

	py::list Bvec;
	for(int j=0; j < B_DIM; j++) {
		for(int i=0; i < H; i++) {
			Bvec.append(B[i][j]);
		}
	}

	py::list Uvec;
	for(int j=0; j < U_DIM; j++) {
		for(int i=0; i < H-1; i++) {
			Uvec.append(U[i][j]);
		}
	}

	std::string workingDir = boost::filesystem::current_path().normalize().string();

	py::object main_module = py::import("__main__");
	py::object main_namespace = main_module.attr("__dict__");
	py::exec("import sys, os", main_namespace);
	py::exec(py::str("sys.path.append('"+workingDir+"/parameter')"), main_namespace);
	py::object plot_mod = py::import("plot_parameter");
	py::object plot_traj = plot_mod.attr("plot_parameter_trajectory");

	plot_traj(Bvec, Uvec, B_DIM, X_DIM, U_DIM, H);

}




void pythonPlotRobot(std::vector< Matrix<U_DIM> >& U, Matrix<X_DIM,X_DIM> SqrtSigma0, Matrix<X_DIM> x0, Matrix<X_DIM> xGoal)
{
	std::vector< Matrix<B_DIM> > B(T);
	vec(x0, SqrtSigma0, B[0]);
	for (int t = 0; t < T - 1; ++t)
	{
		B[t+1] = beliefDynamics(B[t], U[t]);
	}

	try {
		Py_Initialize();

		py::list Jvec;
		for(int j=0; j < 2; j++) {
			for(int i=0; i < T; i++) {
				Jvec.append(B[i][j]);
			}
		}

		std::string workingDir = boost::filesystem::current_path().normalize().string();

		py::object main_module = py::import("__main__");
		py::object main_namespace = main_module.attr("__dict__");
		py::exec("import sys, os", main_namespace);
		py::exec(py::str("sys.path.append('"+workingDir+"/parameter')"), main_namespace);
		py::object plot_mod = py::import("plot_parameter");
		py::object plot_robot = plot_mod.attr("plot_robot");

		plot_robot(Jvec, dynamics::length1, dynamics::length2, T);
	}
	catch(py::error_already_set const &)
	{
		PyErr_Print();
	}

}



#endif
