#include <iostream>
#include <iomanip>
#include <fstream>
#include <ctime>
#include <symbolic/casadi.hpp>
#include <symbolic/stl_vector_tools.hpp>
#include <cstdlib>

#include "../slam.h"

//using namespace CasADi;
//using namespace std;


// dynfunc and obsfunc moved to arm-dynobsjac for analytical Jacobians
inline CasADi::SXMatrix dynfunc(const CasADi::SXMatrix& x_t, const CasADi::SXMatrix& u_t)
{
	CasADi::SXMatrix x_tp1 = x_t;

  	x_tp1(0) += u_t(0) * DT * cos(x_t(2) + u_t(1));
  	x_tp1(1) += u_t(0) * DT * sin(x_t(2) + u_t(1));
  	x_tp1(2) += u_t(0) * DT * sin(u_t(1))/config::WHEELBASE;

  	return x_tp1;
}

// Jacobians: df(x,u,q)/dx, df(x,u,q)/dq
inline void linearizeDynamics(const CasADi::SXMatrix& x, const CasADi::SXMatrix& u, CasADi::SXMatrix& A, CasADi::SXMatrix& M)
{
	//g is control input steer angle
	CasADi::SXMatrix s = sin(u(1)+x(2));
	CasADi::SXMatrix c= cos(u(1)+x(2));

	CasADi::SXMatrix vts= u(0)*DT*s;
	CasADi::SXMatrix vtc= u(0)*DT*c;

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
inline void linearizeObservation(const CasADi::SXMatrix& x, CasADi::SXMatrix& H, CasADi::SXMatrix& N)
{
	for (int i=0; i < L_DIM; i+=2) {
		CasADi::SXMatrix dx = x(C_DIM+i) - x(0);
		CasADi::SXMatrix dy = x(C_DIM+i+1) - x(1);
		CasADi::SXMatrix d2 = dx*dx + dy*dy + 1e-10;
		CasADi::SXMatrix d = sqrt(d2 + 1e-10);
		CasADi::SXMatrix xd = dx/d;
		CasADi::SXMatrix yd = dy/d;
		CasADi::SXMatrix xd2 = dx/d2;
		CasADi::SXMatrix yd2 = dy/d2;


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

inline CasADi::SXMatrix deltaMatrix(const CasADi::SXMatrix& x) {
	CasADi::SXMatrix delta(Z_DIM, Z_DIM);
	CasADi::SXMatrix l0, l1, dist;
	for(int i=C_DIM; i < X_DIM; i += 2) {
		l0 = x(i);
		l1 = x(i+1);

		dist = sqrt((x(0) - l0)*(x(0) - l0) + (x(1) - l1)*(x(1) - l1));

		CasADi::SXMatrix signed_dist = 1/(1+exp(-config::ALPHA_OBS*(config::MAX_RANGE-dist)));
		delta(i-C_DIM,i-C_DIM) = signed_dist;
		delta(i-C_DIM+1,i-C_DIM+1) = signed_dist;
	}

	return delta;
}

inline void EKF(const CasADi::SXMatrix& x_t, const CasADi::SXMatrix& u_t, const CasADi::SXMatrix& Sigma_t, CasADi::SXMatrix& x_tp1, CasADi::SXMatrix& Sigma_tp1)
{
	CasADi::SXMatrix A(X_DIM,X_DIM), M(X_DIM,Q_DIM), QC(U_DIM,U_DIM);
	linearizeDynamics(x_t, u_t, A, M);


	QC(0,0) = config::VELOCITY_NOISE*config::VELOCITY_NOISE;
	QC(1,1) = config::TURNING_NOISE*config::TURNING_NOISE;

	Sigma_tp1 = mul(mul(A,Sigma_t),trans(A)) + mul(mul(M,QC),trans(M));

	x_tp1 = dynfunc(x_t, u_t);


	CasADi::SXMatrix H(Z_DIM,X_DIM), N(Z_DIM,R_DIM), RC(Z_DIM,Z_DIM);

	linearizeObservation(x_tp1, H, N);

	CasADi::SXMatrix delta = deltaMatrix(x_tp1);

	for(int i=0; i < R_DIM-1; i += 2) {
		RC(i,i) = config::OBS_DIST_NOISE*config::OBS_DIST_NOISE;
		RC(i+1,i+1) = config::OBS_ANGLE_NOISE*config::OBS_ANGLE_NOISE;
	}

	//K = ((Sigma_tp1*~H*delta)/(delta*H*Sigma_tp1*~H*delta + RC))*delta;
	CasADi::SXMatrix K = mul(mul(mul(Sigma_tp1, mul(trans(H), delta)), solve(mul(delta, mul(H, mul(Sigma_tp1, mul(trans(H), delta)))) + RC, CasADi::SXMatrix(CasADi::DMatrix::eye(Z_DIM)))), delta);

	Sigma_tp1 = Sigma_tp1 - mul(K,mul(H,Sigma_tp1));

}

inline CasADi::SXMatrix beliefDynamics(const CasADi::SXMatrix& b_t, const CasADi::SXMatrix& u_t) {
	int index = 0;
	CasADi::SXMatrix x_t(X_DIM,1), SqrtSigma_t(X_DIM,X_DIM), Sigma_t(X_DIM,X_DIM);

	for(int i=0; i < X_DIM; ++i) { x_t(i,0) = b_t(i,0); }
	index = X_DIM;
	for (int j = 0; j < X_DIM; ++j) {
		for (int i = j; i < X_DIM; ++i) {
			SqrtSigma_t(i,j) = b_t(index,0);
			SqrtSigma_t(j,i) = b_t(index,0);
			++index;
		}
	}
	Sigma_t = mul(SqrtSigma_t,SqrtSigma_t);

	CasADi::SXMatrix x_tp1(X_DIM,1), Sigma_tp1(X_DIM,X_DIM), SqrtSigma_tp1(X_DIM,X_DIM);
	EKF(x_t, u_t, Sigma_t, x_tp1, Sigma_tp1);

	CasADi::SXMatrix b_fullsigma_tp1(X_DIM+(X_DIM*X_DIM),1);



	for(int i=0; i < X_DIM; ++i) { b_fullsigma_tp1(i,0) = x_tp1(i,0); }
	index = X_DIM;
	for (int j = 0; j < X_DIM; ++j) {
		for (int i = 0; i < X_DIM; ++i) {
			b_fullsigma_tp1(index++,0) = Sigma_tp1(i,j);
		}
	}
	return b_fullsigma_tp1;
}

inline CasADi::SXFunction casadiBeliefDynamicsFunc() {
	CasADi::SXMatrix b = CasADi::ssym("b",B_DIM,1);
	CasADi::SXMatrix u = CasADi::ssym("u",U_DIM,1);

	CasADi::SXMatrix bd = beliefDynamics(b, u);

	// Create functions
	std::vector<CasADi::SXMatrix> inp;
	inp.push_back(b);
	inp.push_back(u);

	CasADi::SXFunction bd_fcn(inp,bd);
	bd_fcn.init();

	return bd_fcn;
}


