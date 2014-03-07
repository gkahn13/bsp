#include "slam-state.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <ctime>
#include <symbolic/casadi.hpp>
#include <symbolic/stl_vector_tools.hpp>
#include <cstdlib>

//#define TIMESTEPS 15
#define NUM_LANDMARKS 20
#define NUM_WAYPOINTS 4

#define C_DIM 3 // car dimension [x, y, theta]
#define P_DIM 2 // Position dimension [x,y]
#define L_DIM 2*NUM_LANDMARKS

#define X_DIM (C_DIM+L_DIM)
#define U_DIM 2
#define Z_DIM L_DIM
#define Q_DIM 2
#define R_DIM L_DIM

int T_CASADI;
const double DT = 1;

namespace config {
const double V = 3;
const double MAXG = 30*M_PI/180.;
const double RATEG = 20*M_PI/180.;
const double WHEELBASE = 4;
const double DT_CONTROLS = 0.025;

const double VELOCITY_NOISE = 0.3;
const double TURNING_NOISE = 3.0*M_PI/180.;

const double MAX_RANGE = 5.0;
const double DT_OBSERVE = 8*DT_CONTROLS;

const double OBS_DIST_NOISE = 1 * 0.1;
const double OBS_ANGLE_NOISE = 1 * 1.0*M_PI/180.;

const double ALPHA_OBS = .75;
}

using namespace CasADi;
using namespace std;


// dynfunc and obsfunc moved to arm-dynobsjac for analytical Jacobians
SXMatrix dynfunc(const SXMatrix& x_t, const SXMatrix& u_t)
{
	SXMatrix x_tp1 = x_t;

  	x_tp1(0) += u_t(0) * DT * cos(x_t(2) + u_t(1));
  	x_tp1(1) += u_t(0) * DT * sin(x_t(2) + u_t(1));
  	x_tp1(2) += u_t(0) * DT * sin(u_t(1))/config::WHEELBASE;

  	return x_tp1;
}

// Jacobians: df(x,u,q)/dx, df(x,u,q)/dq
void linearizeDynamics(const SXMatrix& x, const SXMatrix& u, SXMatrix& A, SXMatrix& M)
{
	//g is control input steer angle
	SXMatrix s = sin(u(1)+x(2));
	SXMatrix c= cos(u(1)+x(2));

	SXMatrix vts= u(0)*DT*s;
	SXMatrix vtc= u(0)*DT*c;

	M(0, 0) = DT*c;
	M(0, 1) = -vts;
	M(1, 0) = DT*s;
	M(1, 1) = vtc;
	M(2, 0) = DT*sin(u(1))/config::WHEELBASE;
	M(2, 1) = u(0)*DT*cos(u(1))/config::WHEELBASE;

	for(int i=0; i < X_DIM; ++i) {
		A(i,i) = 1;
	}

	A(0,0) = 1;
	A(1,1) = 1;
	A(2,2) = 1;
	A(0,2) = -vts;
	A(1,2) = vtc;
}

// Jacobians: dh(x,r)/dx, dh(x,r)/dr
void linearizeObservation(const SXMatrix& x, SXMatrix& H, SXMatrix& N)
{
	for (int i=0; i < L_DIM; i+=2) {
		SXMatrix dx = x(C_DIM+i) - x(0);
		SXMatrix dy = x(C_DIM+i+1) - x(1);
		SXMatrix d2 = dx*dx + dy*dy + 1e-10;
		SXMatrix d = sqrt(d2 + 1e-10);
		SXMatrix xd = dx/d;
		SXMatrix yd = dy/d;
		SXMatrix xd2 = dx/d2;
		SXMatrix yd2 = dy/d2;


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

	for(int i=0; i < Z_DIM; ++i) {
		N(i,i) = 1.0;
	}

}

SXMatrix deltaMatrix(const SXMatrix& x) {
	SXMatrix delta(Z_DIM, Z_DIM);
	SXMatrix l0, l1, dist;
	for(int i=C_DIM; i < X_DIM; i += 2) {
		l0 = x(i);
		l1 = x(i+1);

		dist = sqrt((x(0) - l0)*(x(0) - l0) + (x(1) - l1)*(x(1) - l1));

		SXMatrix signed_dist = 1/(1+exp(-config::ALPHA_OBS*(config::MAX_RANGE-dist)));
		delta(i-C_DIM,i-C_DIM) = signed_dist;
		delta(i-C_DIM+1,i-C_DIM+1) = signed_dist;
	}

	return delta;
}

void EKF(const SXMatrix& x_t, const SXMatrix& u_t, const SXMatrix& Sigma_t, SXMatrix& x_tp1, SXMatrix& Sigma_tp1)
{
	SXMatrix A(X_DIM,X_DIM), M(X_DIM,Q_DIM), QC(U_DIM,U_DIM);
	linearizeDynamics(x_t, u_t, A, M);


	QC(0,0) = config::VELOCITY_NOISE*config::VELOCITY_NOISE;
	QC(1,1) = config::TURNING_NOISE*config::TURNING_NOISE;

	Sigma_tp1 = mul(mul(A,Sigma_t),trans(A)) + mul(mul(M,QC),trans(M));

	x_tp1 = dynfunc(x_t, u_t);


	SXMatrix H(Z_DIM,X_DIM), N(Z_DIM,R_DIM), RC(Z_DIM,Z_DIM);

	linearizeObservation(x_tp1, H, N);

	SXMatrix delta = deltaMatrix(x_tp1);

	for(int i=0; i < R_DIM-1; i += 2) {
		RC(i,i) = config::OBS_DIST_NOISE*config::OBS_DIST_NOISE;
		RC(i+1,i+1) = config::OBS_ANGLE_NOISE*config::OBS_ANGLE_NOISE;
	}

	//K = ((Sigma_tp1*~H*delta)/(delta*H*Sigma_tp1*~H*delta + RC))*delta;
	SXMatrix K = mul(mul(mul(Sigma_tp1, mul(trans(H), delta)), solve(mul(delta, mul(H, mul(Sigma_tp1, mul(trans(H), delta)))) + RC, SXMatrix(DMatrix::eye(Z_DIM)))), delta);

	Sigma_tp1 = Sigma_tp1 - mul(K,mul(H,Sigma_tp1));

}

void generateCode(FX fcn, const std::string& name){
	cout << "Generating code for " << name << endl;

	fcn.generateCode(name + ".c");
}


// params[0] = alpha_belief
// params[1] = alpha_control
// params[2] = alpha_final_belief
SXMatrix costfunc(const SXMatrix& XU, const SXMatrix& Sigma_0, const SXMatrix& landmarks, const SXMatrix& params)
{
	SXMatrix cost = 0;

	SXMatrix x_tp1(X_DIM,1);
	SXMatrix Sigma_t = Sigma_0, Sigma_tp1(X_DIM,X_DIM);
	SXMatrix c_t(C_DIM,1), x_t(X_DIM,1), u_t(U_DIM,1);

	int offset = 0;

	for (int t = 0; t < (T_CASADI-1); ++t)
	{
		c_t = XU(Slice(offset,offset+C_DIM));
		offset += C_DIM;
		u_t = XU(Slice(offset,offset+U_DIM));
		offset += U_DIM;

		for(int i=0; i < X_DIM; ++i) {
			x_t(i) = (i < C_DIM) ? c_t(i) : landmarks(i-C_DIM);
		}

		cost += params[0]*trace(Sigma_t);
		cost += params[1]*inner_prod(u_t, u_t);

		EKF(x_t, u_t, Sigma_t, x_tp1, Sigma_tp1);
		Sigma_t = Sigma_tp1;
	}

	cost += params[2]*trace(Sigma_t);

	return cost;
}

//SXMatrix meritfunc(const SXMatrix& XU, const SXMatrix& Sigma_0, const SXMatrix& landmarks, const SXMatrix& params)
//{
//	SXMatrix cost = 0;
//
//	SXMatrix x_tp1(X_DIM,1);
//	SXMatrix Sigma_t = Sigma_0, Sigma_tp1(X_DIM,X_DIM);
//	SXMatrix c_t(C_DIM,1), x_t(X_DIM,1), u_t(U_DIM,1);
//
//	int offset = 0;
//
//	c_t = XU(Slice(offset,offset+C_DIM));
//	offset += C_DIM;
//
//	for(int i=0; i < X_DIM; ++i) {
//		x_t(i) = (i < C_DIM) ? c_t(i) : landmarks(i-C_DIM);
//	}
//
//	for (int t = 0; t < (T_CASADI-1); ++t)
//	{
//		u_t = XU(Slice(offset,offset+U_DIM));
//		offset += U_DIM;
//		c_t = XU(Slice(offset,offset+C_DIM));
//		offset += C_DIM;
//
//		for(int i=0; i < X_DIM; ++i) {
//			x_t(i) = (i < C_DIM) ? c_t(i) : landmarks(i-C_DIM);
//		}
//
//		cost += params[0]*trace(Sigma_t);
//		cost += params[1]*inner_prod(u_t, u_t);
//
//		EKF(x_t, u_t, Sigma_t, x_tp1, Sigma_tp1);
//		Sigma_t = Sigma_tp1;
//	}
//
//	cost += params[2]*trace(Sigma_t);
//
//	return cost;
//}

void test() {
	T_CASADI = 15;

	SXMatrix x_t(X_DIM,1), u_t(U_DIM,1), Sigma_t(X_DIM,X_DIM), x_tp1(X_DIM,1), Sigma_tp1(X_DIM,X_DIM);

	for(int i=0; i < X_DIM; ++i) { x_t(i) = 0; }
	x_t(3) = 30; x_t(4) = -10;
	x_t(5) = 70; x_t(6) = 12.5;
	x_t(7) = 20; x_t(8) = 10;

	u_t(0) = (60.0 - 0.0)/((double)(T_CASADI-1));
	u_t(1) = 0;

	cout << "x0: " << x_t << endl;
	cout << "u0: " << u_t << endl;

	for(int i = 0; i < C_DIM; ++i) { Sigma_t(i,i) = .1*.1; }
	for(int i = 0; i < L_DIM; ++i) { Sigma_t(C_DIM+i,C_DIM+i) = 1*1; }

	for(int t=0; t < T_CASADI-1; ++t) {
		cout << "\n\n\n\nt: " << t << endl;
		EKF(x_t, u_t, Sigma_t, x_tp1, Sigma_tp1);
		x_t = x_tp1;
		Sigma_t = Sigma_tp1;
	}

	//cout << "Sigma_t" << endl << Sigma_t << endl;
	//cout << "x_tp1" << endl << x_tp1 << endl;
	//cout << "Sigma_tp1" << endl << Sigma_tp1 << endl;
}

SXFunction casadiCostFunc(int timesteps) {
	T_CASADI = timesteps;

	int nXU = T_CASADI*C_DIM+(T_CASADI-1)*U_DIM;
	SXMatrix XU = ssym("XU",nXU,1);
	SXMatrix Sigma_0 = ssym("S0",X_DIM,X_DIM);
	SXMatrix landmarks = ssym("landmarks",L_DIM);
	SXMatrix params = ssym("params",3); // alpha_control, alpha_belief, alpha_final_belief

	// Objective
	SXMatrix f = costfunc(XU, Sigma_0, landmarks, params);

	// Create functions
	vector<SXMatrix> inp;
	inp.push_back(XU);
	inp.push_back(Sigma_0);
	inp.push_back(landmarks);
	inp.push_back(params);

	SXFunction f_fcn(inp,f);
	f_fcn.init();

	return f_fcn;
}

SXFunction casadiCostGradFunc(int timesteps) {
	T_CASADI = timesteps;

	int nXU = T_CASADI*C_DIM+(T_CASADI-1)*U_DIM;
	SXMatrix XU = ssym("XU",nXU,1);
	SXMatrix Sigma_0 = ssym("S0",X_DIM,X_DIM);
	SXMatrix landmarks = ssym("landmarks",L_DIM);
	SXMatrix params = ssym("params",3); // alpha_control, alpha_belief, alpha_final_belief

	// Objective
	SXMatrix f = costfunc(XU, Sigma_0, landmarks, params);

	SXMatrix grad_f = gradient(f,XU);

	// Create functions
	vector<SXMatrix> inp;
	inp.push_back(XU);
	inp.push_back(Sigma_0);
	inp.push_back(landmarks);
	inp.push_back(params);

	vector<SXMatrix> out;
	out.push_back(f);
	out.push_back(grad_f);
	SXFunction grad_f_fcn(inp,out);
	grad_f_fcn.init();

	return grad_f_fcn;
}

