#ifndef __POMDP_PENDULUM_H__
#define __POMDP_PENDULUM_H__

#define DIM 2
#define X_DIM 8
#define U_DIM 2
#define Z_DIM 4
#define K_DIM 4
#define S_DIM (((X_DIM+1)*X_DIM)/2)
#define B_DIM (X_DIM + S_DIM)
#define SIM_X_DIM (X_DIM-K_DIM)

#include "Dynamics-2dplanararm.hpp"

//Dynamics::Dynamics() : DynamicsBase("2dplanararm-lengths"), gravity(9.86), gravity_est(9.86),
//	mass1(0.092), mass1_est(0.105), mass_center1(0.062), mass_center1_est(0.062), length1(0.15), length1_est(0.05), damping1(0.015), damping1_est(0.015), coulomb1(0.006), coulomb1_est(0.006), rotational1_inertia(0.00064), rotational1_inertia_est(0.00064), motor1_inertia(0.00000065), motor1_inertia_est(0.00000065), motor1_torque_const(0.0077), motor1_torque_const_est(0.0077), gear1_ratio(70), gear1_ratio_est(70),
//	mass2(0.077), mass2_est(0.089), mass_center2(0.036), mass_center2_est(0.036), length2(0.15), length2_est(0.05), damping2(0.015), damping2_est(0.015), coulomb2(0.006), coulomb2_est(0.006), rotational2_inertia(0.0003), rotational2_inertia_est(0.0003), motor2_inertia(0.00000065), motor2_inertia_est(0.00000065), motor2_torque_const(0.0077), motor2_torque_const_est(0.0077), gear2_ratio(70), gear2_ratio_est(70)
Dynamics::Dynamics() : DynamicsBase("2dplanararm-lengths"), gravity(9.86), gravity_est(9.86),
	mass1(0.092), mass1_est(0.079), mass_center1(0.062), mass_center1_est(0.062), length1(0.15), length1_est(0.3), damping1(0.015), damping1_est(0.015), coulomb1(0.006), coulomb1_est(0.006), rotational1_inertia(0.00064), rotational1_inertia_est(0.00064), motor1_inertia(0.00000065), motor1_inertia_est(0.00000065), motor1_torque_const(0.0077), motor1_torque_const_est(0.0077), gear1_ratio(70), gear1_ratio_est(70),
	mass2(0.077), mass2_est(0.065), mass_center2(0.036), mass_center2_est(0.036), length2(0.15), length2_est(0.25), damping2(0.015), damping2_est(0.015), coulomb2(0.006), coulomb2_est(0.006), rotational2_inertia(0.0003), rotational2_inertia_est(0.0003), motor2_inertia(0.00000065), motor2_inertia_est(0.00000065), motor2_torque_const(0.0077), motor2_torque_const_est(0.0077), gear2_ratio(70), gear2_ratio_est(70)
	{
	horizon = 25;

	xParamReal(0,0) = 1/length1;
	xParamReal(1,0) = 1/length2	;
	xParamReal(2,0) = 1/mass1;
	xParamReal(3,0) = 1/mass2;

	Rint(0,0) = 1.0;
	Rint(1,1) = 1.0;

	Qint(0,0) = 1.0;
	Qint(1,1) = 1.0;
	Qint(2,2) = 1.0;
	Qint(3,3) = 1.0;
	Qint(4,4) = 0.0;
	Qint(5,5) = 0.0;
	Qint(6,6) = 0.0;
	Qint(7,7) = 0.0;

	QintVariance(0,0) = 0.0;
	QintVariance(1,1) = 0.0;
	QintVariance(2,2) = 0.0;
	QintVariance(3,3) = 0.0;
	QintVariance(4,4) = 100.0;
	QintVariance(5,5) = 100.0;
	QintVariance(6,6) = 100.0;
	QintVariance(7,7) = 100.0;

	QGoal(0,0) = 1.0;
	QGoal(1,1) = 1.0;
	QGoal(2,2) = 0.0;
	QGoal(3,3) = 0.0;
	QGoal(4,4) = 0.0;
	QGoal(5,5) = 0.0;
	QGoal(6,6) = 0.0;
	QGoal(7,7) = 0.0;

	QGoalVariance(0,0) = 0.0;
	QGoalVariance(1,1) = 0.0;
	QGoalVariance(2,2) = 0.0;
	QGoalVariance(3,3) = 0.0;
	QGoalVariance(4,4) = 100.0;
	QGoalVariance(5,5) = 100.0;
	QGoalVariance(6,6) = 100.0;
	QGoalVariance(7,7) = 100.0;

	experiment_t e = new double[2*K_DIM];
	e[0] = length1;
	e[1] = length1_est;
	e[2] = length2;
	e[3] = length2_est;
	e[4] = mass1;
	e[5] = mass1_est;
	e[6] = mass2;
	e[7] = mass2_est;
	initExperiment(e);

	Sigma0(4,4) = 0.5;
	Sigma0(5,5) = 0.5;
	Sigma0(6,6) = 1.0;
	Sigma0(7,7) = 1.0;

	control_bounds[0] = 1.0;
	control_bounds[1] = -1.0;
}

void Dynamics::setupExperiments(experiments_t& experiments) {
	double *e = new double[2*K_DIM];
	e[0] = length1;
	e[1] = length1_est;
	e[2] = length2;
	e[3] = length2_est;
	e[4] = mass1;
	e[5] = mass1_est;
	e[6] = mass2;
	e[7] = mass2_est;

	experiments.push_back(e);
}

inline SymmetricMatrix<X_DIM> Dynamics::varM(const Matrix<X_DIM>& x, const Matrix<U_DIM>& u) const {
	// TODO DWEBB -- Double check this and then ask Jur to verify it.
	SymmetricMatrix<X_DIM> S = identity<X_DIM>();
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

#define VARN(_dim) \
inline SymmetricMatrix<Z_DIM> Dynamics::varN(const Matrix<_dim>& x, const Matrix<U_DIM>& u) const { \
	SymmetricMatrix<Z_DIM> S = identity<Z_DIM>();	\
	S(0,0) = 0.0001;	\
	S(1,1) = 0.0001;	\
	S(2,2) = 0.00001;	\
	S(3,3) = 0.00001;	\
	return S;	\
}

/*
	S(4,4) = 0.00001;	\
	S(5,5) = 0.00001;	\
	S(6,6) = 0.00001;	\
	S(7,7) = 0.00001;	\
	S(8,8) = 0.00001;	\
	S(9,9) = 0.00001;	\
*/

VARN(X_DIM);
VARN(SIM_X_DIM);

inline Matrix<Z_DIM> h_custom(const Matrix<X_DIM>& x, double length1, double length2)
{
	Matrix<Z_DIM> z;

	float x0 = x[0];
	float x1 = x[1];
	float x2 = x[2];
	float x3 = x[3];
	float x4 = (x[4] == 0 ? 0 : 1/x[4]);
	float x5 = (x[5] == 0 ? 0 : 1/x[5]);
	float x0plus1 = x0 + x1;
	float cosx0 = cos(x0);
	float sinx0 = sin(x0);
	float cosx1 = cos(x1);
	float sinx1 = sin(x1);
	float cosx0plus1 = cos(x0plus1);
	float sinx0plus1 = sin(x0plus1);

	z[0] = x4*cosx0 + x5*cosx1;
	z[1] = x4*sinx0 + x5*sinx1;
	z[2] = x2;
	z[3] = x3;

	return z;
}

/*
	z[2] = x2;
	z[3] = x3;
	z[4] = x4*cosx0;
	z[5] = x4*sinx0;
	z[8] = x5*cosx0plus1;
	z[9] = x5*sinx0plus1;
*/

inline Matrix<Z_DIM> h_custom(const Matrix<SIM_X_DIM>& x, double length1, double length2)
{
	Matrix<Z_DIM> z;

	float x0 = x[0];
	float x1 = x[1];
	float x2 = x[2];
	float x3 = x[3];
	float x0plus1 = x0 + x1;
	float cosx0 = cos(x0);
	float sinx0 = sin(x0);
	float cosx1 = cos(x1);
	float sinx1 = sin(x1);
	float cosx0plus1 = cos(x0plus1);
	float sinx0plus1 = sin(x0plus1);

	z[0] = length1*cosx0 + length2*cosx1;
	z[1] = length1*sinx0 + length2*sinx1;
	z[2] = x2;
	z[3] = x3;

	return z;
}

/*
	z[2] = x2;
	z[3] = x3;
	z[4] = length1*cosx0;
	z[5] = length1*sinx0;
	z[8] = length2*cosx0plus1;
	z[9] = length2*sinx0plus1;
*/

#define HFN(_dim)	\
inline Matrix<Z_DIM> Dynamics::h(const Matrix<_dim>& x) const {	\
	/* This is a hack to ensure the correct observation function is called since I can't specialize the observation function at this point. */	\
	return h_custom(x, length1, length2);	\
}

HFN(X_DIM);
HFN(SIM_X_DIM);

std::string Dynamics::getDirectoryName() const {
	std::ostringstream oss;
	oss << name << "_length1_" << length1 << "_est_" << length1_est << "_length2_" << length2 << "_est_" << length2_est << "_mass1_" << mass1 << "_est_" << mass1_est << "_mass2_" << mass2 << "_est_" << mass2_est;
	return oss.str();
}

void Dynamics::initExperiment(const experiment_t& e) {
	length1 = e[0];
	length1_est = e[1];
	length2 = e[2];
	length2_est = e[3];
	mass1 = e[4];
	mass1_est = e[5];
	mass2 = e[6];
	mass2_est = e[7];

	x0[0] = -M_PI/2.0; x0[1] = -M_PI/2.0; x0[2] = 0; x0[3] = 0; x0[4] = 1/length1_est; x0[5] = 1/length2_est; x0[6] = 1/mass1_est; x0[7] = 1/mass2_est;
	xGoal[0] = -M_PI/2.0; xGoal[1] = -M_PI/2.0; xGoal[2] = 0.0; xGoal[3] = 0.0; xGoal[4] = 1/length1_est; xGoal[5] = 1/length2_est; xGoal[6] = 1/mass1_est; xGoal[7] = 1/mass2_est;
}

inline Matrix<X_DIM> Dynamics::f(double step, const Matrix<X_DIM>& x, const Matrix<U_DIM>& u) const {
	// RK4 integration
	Matrix<SIM_X_DIM> k1, k2, k3, k4, xinit;

	xinit = x.subMatrix<SIM_X_DIM>(0,0);

	double x4 = x[4];
	double x5 = x[5];
	double x6 = x[6];
	double x7 = x[7];

	double length1 = (x4 == 0 ? 0.0 : 1/x4);
	double length2 = (x5 == 0 ? 0.0 : 1/x5);
	double mass1 = (x6 == 0 ? 0.0 : 1/x6);
	double mass2 = (x7 == 0 ? 0.0 : 1/x7);

	k1 = g(xinit, u, mass1, mass2, mass_center1, mass_center2, length1, length2, damping1, damping2, coulomb1, coulomb2, rotational1_inertia, rotational2_inertia, motor1_inertia, motor2_inertia, motor1_torque_const, motor2_torque_const, gear1_ratio, gear2_ratio);
	k2 = g(xinit + 0.5*step*k1, u, mass1, mass2, mass_center1, mass_center2, length1, length2, damping1, damping2, coulomb1, coulomb2, rotational1_inertia, rotational2_inertia, motor1_inertia, motor2_inertia, motor1_torque_const, motor2_torque_const, gear1_ratio, gear2_ratio);
	k3 = g(xinit + 0.5*step*k2, u, mass1, mass2, mass_center1, mass_center2, length1, length2, damping1, damping2, coulomb1, coulomb2, rotational1_inertia, rotational2_inertia, motor1_inertia, motor2_inertia, motor1_torque_const, motor2_torque_const, gear1_ratio, gear2_ratio);
	k4 = g(xinit + step*k3, u, mass1, mass2, mass_center1, mass_center2, length1, length2, damping1, damping2, coulomb1, coulomb2, rotational1_inertia, rotational2_inertia, motor1_inertia, motor2_inertia, motor1_torque_const, motor2_torque_const, gear1_ratio, gear2_ratio);

	Matrix<X_DIM> xNew = zeros<X_DIM,1>();
	xNew.insert(0, 0, xinit + step*(k1 + 2.0*(k2 + k3) + k4)/6.0);
	xNew[4] = x4;
	xNew[5] = x5;
	xNew[6] = x6;
	xNew[7] = x7;

	return xNew;
}

inline void Dynamics::linearizeDynamics(double step, const Matrix<X_DIM>& xBar, const Matrix<U_DIM>& uBar, Matrix<X_DIM>& c, Matrix<X_DIM, X_DIM>& A, Matrix<X_DIM, U_DIM>& B, SymmetricMatrix<X_DIM>& M, unsigned int flag) const {
	if (flag & COMPUTE_c) c = f(step, xBar, uBar);
	if (flag & COMPUTE_A) {
		Matrix<X_DIM> xr(xBar), xl(xBar);
		for (size_t i = 0; i < X_DIM; ++i) {
			xr[i] += step; xl[i] -= step;
			A.insert(0,i, (f(step, xr, uBar) - f(step, xl, uBar)) / (2.0*step));
			xr[i] = xl[i] = xBar[i];
		}
		//A = dfdx(step, *this, xBar, uBar); //identity<X_DIM>();
	}
	if (flag & COMPUTE_B) {
		//B = dfdu(step, *this, xBar, uBar); //identity<U_DIM>();
		Matrix<U_DIM> ur(uBar), ul(uBar);
		for (size_t i = 0; i < U_DIM; ++i) {
			ur[i] += step; ul[i] -= step;
			B.insert(0,i, (f(step, xBar, ur) - f(step, xBar, ul)) / (2.0*step));
			ur[i] = ul[i] = uBar[i];
		}
	}
	if (flag & COMPUTE_M) M = varM(xBar, uBar);
}

inline void Dynamics::linearizeObservation(double step, const Matrix<X_DIM>& xBar, const Matrix<U_DIM>& uBar, Matrix<Z_DIM, X_DIM>& H, SymmetricMatrix<Z_DIM>& N) const {
	//H = dhdx<X_DIM>(step, h<X_DIM>, xBar); // identity<X_DIM>();
	//H = dhdx<X_DIM, Z_DIM>(*this, step, xBar); // identity<X_DIM>();
	Matrix<X_DIM> xr(xBar), xl(xBar);
	for (size_t i = 0; i < X_DIM; ++i) {
		xr[i] += step; xl[i] -= step;
		H.insert(0,i, (h(xr) - h(xl)) / (2.0*step));
		xr[i] = xl[i] = xBar[i];
	}
	N = varN(xBar, uBar);
}

// Compute closed form quadratic finalCost function around around b
inline void Dynamics::quadratizeFinalCost(const Matrix<X_DIM>& xBar, const SymmetricMatrix<X_DIM>& SigmaBar, double& s, SymmetricMatrix<X_DIM>& S, Matrix<1,X_DIM>& sT, Matrix<1,S_DIM>& tT, unsigned int flag) const {
	if (flag & COMPUTE_S) S = QGoal;
	if (flag & COMPUTE_sT) sT = ~(xBar - xGoal)*QGoal;
	if (flag & COMPUTE_s) s = 0.5*scalar(~(xBar - xGoal)*QGoal*(xBar - xGoal)) + scalar(vecTh(QGoalVariance)*vec(SigmaBar));
	if (flag & COMPUTE_tT) tT = vecTh(QGoalVariance);
}

// Compute closed form quadratic cost function around around b and u
inline bool Dynamics::quadratizeCost(double step, const Matrix<X_DIM>& xBar, const SymmetricMatrix<X_DIM>& SigmaBar, const Matrix<U_DIM>& uBar, double& q, SymmetricMatrix<X_DIM>& Q, 
						   SymmetricMatrix<U_DIM>& R, Matrix<U_DIM, X_DIM>& P, Matrix<1,X_DIM>& qT, Matrix<1,U_DIM>& rT, Matrix<1,S_DIM>& pT, unsigned int flag) const {
	if (flag & COMPUTE_Q) Q = Qint; // penalize uncertainty
	if (flag & COMPUTE_R) R = Rint;
	if (flag & COMPUTE_P) P = zeros<U_DIM,X_DIM>();
	if (flag & COMPUTE_qT) qT = ~(xBar - xGoal)*Qint;
	if (flag & COMPUTE_rT) rT = ~uBar*Rint;
	if (flag & COMPUTE_pT) pT = vecTh(QintVariance);
	if (flag & COMPUTE_q) q = 0.5*scalar(~(xBar - xGoal)*Qint*(xBar - xGoal)) + 0.5*scalar(~uBar*Rint*uBar) + scalar(vecTh(QintVariance)*vec(SigmaBar));
	return true; 
}

#endif //__POMDP_PENDULUM_H__