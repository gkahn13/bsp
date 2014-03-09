#ifndef __PARAMETER_H__
#define __PARAMETER_H__

#include <fstream>

#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include "util/matrix.h"
#include <cmath>

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
#define DT 0.1

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
#define UT_DIM U_DIM*(T-1)

#define S_DIM (((X_DIM+1)*X_DIM)/2)
#define B_DIM (X_DIM+S_DIM)
#define XU_DIM (X_DIM*T+U_DIM*(T-1))

namespace dynamics {

const double mass1 = 0.5;
const double mass2 = 0.5;


//coefficient of friction 
const double b1 = 0.0;
const double b2 = 0.0;

const double length1 = 0.5;
const double length2 = 0.5;

const double gravity = 9.82;

}


const double diffEps = 0.0078125 / 16;


const int T = TIMESTEPS;
const double INFTY = 1e10;

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

const double alpha_belief = 10, alpha_final_belief = 10, alpha_control = 0.01, alpha_goal_state = 1;
const double alpha_joint_belief = 0, alpha_param_belief = 10,
			  alpha_final_joint_belief = 0, alpha_final_param_belief = 10, alpha_goal_joint_state = 10, alpha_goal_param_state = 1;

double *inputVars, *vars;
std::vector<int> maskIndices;




Matrix<J_DIM> jointdynfunc(const Matrix<J_DIM>& x, const Matrix<U_DIM>& u, double l1, double l2, double m1, double m2){


	double l1_sqr = l1*l1; 
	double l2_sqr = l2*l2; 

	double I1 = m1*l1_sqr/12; 
	double I2 = m2*l2_sqr/12;



	Matrix<U_DIM,U_DIM> A; 

	A(0,0) = l1_sqr*(0.25*m2+m1)+I1; 

	A(0,1) = 0.5*m2*l1*l2*cos(x[0]-x[1]);

	A(1,0) = 0.5*m2*l1*l2*cos(x[0]-x[1]);

	A(1,1) = l2_sqr*0.25*m2+I2; 

	Matrix<U_DIM> B; 

	B[0] = dynamics::gravity*l1*sin(x[0])*(0.5*m1+m2) - 0.5*m2*l1*l2*(x[3]*x[3])*sin(x[0]-x[1])+u[0]-dynamics::b1*x[2];

	B[1] = 0.5*m2*l2*(l1*(x[2]*x[2])*sin(x[0]-x[1])+dynamics::gravity*sin(x[1])) + u[1]-dynamics::b2*x[3];
   

	Matrix<U_DIM> xd = A%B;


	Matrix<J_DIM> jNew = zeros<J_DIM,1>();

	jNew[0] = x[2];
	jNew[1] = x[3];
	jNew[2] = xd[0];
	jNew[3] = xd[1];

	return jNew;


}

// for both joints and params
Matrix<X_DIM> dynfunc(const Matrix<X_DIM>& x, const Matrix<U_DIM>& u, const Matrix<Q_DIM> q)
{

	Matrix<J_DIM> k1, k2, k3, k4, jinit;

	jinit = x.subMatrix<J_DIM>(0,0);

	double x4 = x[4];
	double x5 = x[5];
	double x6 = x[6];
	double x7 = x[7];
	
	double l1, l2, m1, m2;
	l1 = (x4 == 0 ? 0.0 : 1/x4);
	l2 = (x5 == 0 ? 0.0 : 1/x5);
	m1 = (x6 == 0 ? 0.0 : 1/x6);
	m2 = (x7 == 0 ? 0.0 : 1/x7);

	k1 = jointdynfunc(jinit, u, l1, l2, m1, m2);
	k2 = jointdynfunc(jinit + 0.5*DT*k1, u, l1, l2, m1, m2);
	k3 = jointdynfunc(jinit + 0.5*DT*k2, u, l1, l2, m1, m2);
	k4 = jointdynfunc(jinit + DT*k3, u, l1, l2, m1, m2);


	Matrix<X_DIM> xNew = zeros<X_DIM,1>();
	xNew.insert(0, 0, jinit + DT*(k1 + 2.0*(k2 + k3) + k4)/6.0);
	xNew[4] = x4;
	xNew[5] = x5;
	xNew[6] = x6;
	xNew[7] = x7;

	xNew += q*DT; 

	return xNew;
}

// Observation model
Matrix<Z_DIM> obsfunc(const Matrix<X_DIM>& x, const Matrix<R_DIM> & r)
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

	z += r; 

	return z;
}

SymmetricMatrix<Q_DIM> varQ() {

	SymmetricMatrix<Q_DIM> S = identity<Q_DIM>();
	S(0,0) = 0.01;
	S(1,1) = 0.01;
	S(2,2) = 0.01;
	S(3,3) = 0.01;
	S(4,4) = 1e-6;
	S(5,5) = 1e-6;
	S(6,6) = 1e-6;
	S(7,7) = 1e-6;
	return S;
}

SymmetricMatrix<R_DIM> varR(){
	SymmetricMatrix<R_DIM> S = identity<R_DIM>();
	S(0,0) = 1e-2;
	S(1,1) = 1e-2;
	S(2,2) = 1e-3;
	S(3,3) = 1e-3;
	return S;
}


// Jacobians: df(x,u,q)/dx, df(x,u,q)/dq
void linearizeDynamics(const Matrix<X_DIM>& x, const Matrix<U_DIM>& u, const Matrix<Q_DIM>& q, Matrix<X_DIM,X_DIM>& A, Matrix<X_DIM,Q_DIM>& M)
{
	A.reset();
	Matrix<X_DIM> xr(x), xl(x);
	for (size_t i = 0; i < X_DIM; ++i) {
		xr[i] += diffEps; xl[i] -= diffEps;
		A.insert(0,i, (dynfunc(xr, u, q) - dynfunc(xl, u, q)) / (2.0*diffEps));
		xr[i] = x[i]; xl[i] = x[i];
	}
	M.reset();
	Matrix<Q_DIM> qr(q), ql(q);
	for (size_t i = 0; i < Q_DIM; ++i) {
		qr[i] += diffEps; ql[i] -= diffEps;
		M.insert(0,i, (dynfunc(x, u, qr) - dynfunc(x, u, ql)) / (2.0*diffEps));
		qr[i] = q[i]; ql[i] = q[i];
	}
}

// for N = dh(x,r)/dr
// this returns N*~N

// Jacobians: dh(x,r)/dx, dh(x,r)/dr
void linearizeObservation(const Matrix<X_DIM>& x, const Matrix<R_DIM>& r, Matrix<Z_DIM,X_DIM>& H, Matrix<Z_DIM,R_DIM>& N)
{
	H.reset();
	Matrix<X_DIM> xr(x), xl(x);
	for (size_t i = 0; i < X_DIM; ++i) {

		xr[i] += diffEps; xl[i] -= diffEps;
		H.insert(0,i, (obsfunc(xr, r) - obsfunc(xl, r)) / (2.0*diffEps));
		
		xr[i] = x[i]; xl[i] = x[i];
	}

	N.reset();
	Matrix<R_DIM> rr(r), rl(r);
	for (size_t i = 0; i < R_DIM; ++i) {
		rr[i] += diffEps; rl[i] -= diffEps;
		N.insert(0,i, (obsfunc(x, rr) - obsfunc(x, rl)) / (2.0*diffEps));
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

// Belief dynamics
Matrix<B_DIM> beliefDynamics(const Matrix<B_DIM>& b, const Matrix<U_DIM>& u) {
	Matrix<X_DIM> x;
	Matrix<X_DIM,X_DIM> Sigma;
	unVec(b, x, Sigma);

	SymmetricMatrix<Q_DIM> Q = varQ();
	SymmetricMatrix<R_DIM> R = varR(); 
	Matrix<Q_DIM> q;
	Matrix<R_DIM> r; 

	Sigma = Sigma*Sigma;

	Matrix<X_DIM,X_DIM> A;
	Matrix<X_DIM,Q_DIM> M;
	 
	//LOG_INFO("Linearize dynamics");
	linearizeDynamics(x,u,q,A,M);

	//LOG_INFO("x_tp1 computation");
	//std::cout << ~q << std::endl;
	x = dynfunc(x, u,q);
	
	//std::cout << ~x << std::endl;

	Sigma = A*Sigma*~A + Q;

	Matrix<Z_DIM,X_DIM> H = zeros<Z_DIM,X_DIM>();
	Matrix<Z_DIM,R_DIM> N; 

	linearizeObservation(x,r, H,N);
	/*
	std::cout<<"X_BELIEF "<<x<<"\n";

	
	std::cout<<"X_BELIEF "<<(H*Sigma*~H + N*R*~N)<<"\n";

	*/

	//LOG_INFO("Compute K");
	Matrix<X_DIM,Z_DIM> K = Sigma*~H/(H*Sigma*~H + R);
	//LOG_INFO("Finished computing K");

	Sigma = (identity<X_DIM>() - K*H)*Sigma;
	
	Matrix<B_DIM> g;
	vec(x, sqrtm(Sigma), g);

	return g;

}





// returns updated belief based on real current state, estimated current belief, and input
Matrix<B_DIM> executeControlStep(Matrix<X_DIM>& x_t_real, const Matrix<B_DIM>& b_t, const Matrix<U_DIM>& u_t) {
	// useRealParams = true
	// update real state (maximum likelihood)
	SymmetricMatrix<Q_DIM> Q = varQ(); 
	SymmetricMatrix<R_DIM> R = varR(); 
	Matrix<Q_DIM> q =  sampleGaussian(zeros<Q_DIM,1>(),Q);
	Matrix<R_DIM> r = sampleGaussian(zeros<R_DIM,1>(),R);

	Matrix<X_DIM> x_tp1_real = dynfunc(x_t_real,u_t,q);
	x_t_real = x_tp1_real; 
	// sense real state (maximum likelihood)
	Matrix<Z_DIM> z_tp1_real = obsfunc(x_tp1_real,r);

	// now do EKF on belief and incorporate discrepancy
	// using the kalman gain

	q = zeros<Q_DIM,1>();
	r = zeros<R_DIM,1>();

	Matrix<X_DIM> x_t;
	Matrix<X_DIM,X_DIM> SqrtSigma_t;
	unVec(b_t, x_t, SqrtSigma_t);

	Matrix<X_DIM,X_DIM> Sigma_t = SqrtSigma_t * SqrtSigma_t;

	// calculate df/dx and df/dm
	Matrix<X_DIM,X_DIM> A;
	Matrix<X_DIM,Q_DIM> M;
	 
	linearizeDynamics(x_t,u_t,q,A,M);

	// propagate estimated state
	Matrix<X_DIM> x_tp1 = dynfunc(x_t, u_t,q);
	//Matrix<X_DIM,X_DIM> Sigma_tp1 = A*Sigma_t*~A + M*Q*~M;
	Matrix<X_DIM,X_DIM> Sigma_tp1 = A*Sigma_t*~A + Q;
	// calculate dh/dx and dh/dn
	Matrix<Z_DIM,X_DIM> H = zeros<Z_DIM,X_DIM>();
	Matrix<Z_DIM,R_DIM> N; 

	linearizeObservation(x_tp1,r, H,N);

	// calculate Kalman gain for estimated state
	//Matrix<X_DIM,Z_DIM> K = Sigma_tp1*~H/(H*Sigma_tp1*~H + N*R*~N);
	Matrix<X_DIM,Z_DIM> K = Sigma_tp1*~H/(H*Sigma_tp1*~H + R);
	//std::cout<<"x_tp1"<<x_tp1<<"\n";
	// correct the new state using Kalman gain and the observation

	std::cout << "Before update: " << std::endl;
	std::cout << ~x_tp1 << std::endl;

	Matrix<X_DIM> x_tp1_adj = x_tp1 + K*(z_tp1_real - obsfunc(x_tp1,r));

	std::cout << "After update: " << std::endl;
	std::cout<<"x_tp1: "<<~x_tp1_adj<<"\n";
	//Matrix<X_DIM,X_DIM> W = ~(H*Sigma_tp1)*(((H*Sigma_tp1*~H) + N*R*~N) % (H*Sigma_tp1));
	//Matrix<X_DIM,X_DIM> Sigma_tp1_adj = Sigma_tp1 - W;
	Matrix<X_DIM,X_DIM> Sigma_tp1_adj = Sigma_tp1 - K*H*Sigma_tp1;

	Matrix<B_DIM> b_tp1_adj;
	vec(x_tp1_adj, sqrtm(Sigma_tp1_adj), b_tp1_adj);
	return b_tp1_adj;
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
	std::cout<<"462\n";

	std::string workingDir = boost::filesystem::current_path().normalize().string();

	py::object main_module = py::import("__main__");
	py::object main_namespace = main_module.attr("__dict__");
	py::exec("import sys, os", main_namespace);
	py::exec(py::str("sys.path.append('"+workingDir+"/parameter')"), main_namespace);
	std::cout<<"472\n";
	py::object plot_mod = py::import("plot_parameter");
	std::cout<<"472\n";
	py::object plot_traj = plot_mod.attr("plot_parameter_trajectory");
	std::cout<<"472\n";
	plot_traj(Bvec, Uvec, B_DIM, X_DIM, U_DIM, H);

}


void pythonPaperPlot(std::vector< Matrix<U_DIM> >& U0,std::vector< Matrix<B_DIM> >& B0, std::vector< Matrix<U_DIM> >& U1,std::vector< Matrix<B_DIM> >& B1, Matrix<X_DIM,X_DIM> SqrtSigma0, Matrix<X_DIM> x0, int H)
{
	
	Py_Initialize();

	py::list Bvec0;
	for(int j=0; j < B_DIM; j++) {
		for(int i=0; i < H; i++) {
			Bvec0.append(B0[i][j]);
		}
	}

	py::list Uvec0;
	for(int j=0; j < U_DIM; j++) {
		for(int i=0; i < H-1; i++) {
			Uvec0.append(U0[i][j]);
		}
	}


	py::list Bvec1;
	for(int j=0; j < B_DIM; j++) {
		for(int i=0; i < H; i++) {
			Bvec1.append(B1[i][j]);
		}
	}

	py::list Uvec1;
	for(int j=0; j < U_DIM; j++) {
		for(int i=0; i < H-1; i++) {
			Uvec1.append(U1[i][j]);
		}
	}


	std::string workingDir = boost::filesystem::current_path().normalize().string();

	py::object main_module = py::import("__main__");
	py::object main_namespace = main_module.attr("__dict__");
	py::exec("import sys, os", main_namespace);
	py::exec(py::str("sys.path.append('"+workingDir+"/parameter')"), main_namespace);
	py::object plot_mod = py::import("plot_parameter");
	py::object plot_paper = plot_mod.attr("plot_for_paper");

	plot_paper(Bvec0, Uvec0,Bvec1, Uvec1, B_DIM, X_DIM, U_DIM, H);

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
