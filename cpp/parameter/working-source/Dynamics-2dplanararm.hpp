#ifndef __POMDP_2DPLANARARM_H__
#define __POMDP_2DPLANARARM_H__

#include <utility>
#include <vector>

#include "matrix.h"
#include "DynamicsBase.hpp"

// Forward declarations
inline std::pair<double, double> normal();
template <typename T>
inline int sgn(T val);

struct Dynamics : public DynamicsBase<DIM, X_DIM, SIM_X_DIM, U_DIM, Z_DIM, K_DIM, S_DIM, B_DIM> {
	Dynamics();

	virtual std::string getDirectoryName() const;

	virtual void initExperiment(const experiment_t& e);

	virtual inline Matrix<Z_DIM> h(const Matrix<X_DIM>& x) const;
	virtual inline Matrix<Z_DIM> h(const Matrix<SIM_X_DIM>& x) const;

	virtual inline SymmetricMatrix<Z_DIM> varN(const Matrix<X_DIM>& x, const Matrix<U_DIM>& u) const;
	virtual inline SymmetricMatrix<Z_DIM> varN(const Matrix<SIM_X_DIM>& x, const Matrix<U_DIM>& u) const;

	void saturateControl(Matrix<U_DIM>& u) const {
		if (u(0,0) > control_bounds(0,0)) u(0,0) = control_bounds(0,0);
		if (u(0,0) < control_bounds(0,1)) u(0,0) = control_bounds(0,1);
		if (u(1,0) > control_bounds(1,0)) u(1,0) = control_bounds(1,0);
		if (u(1,0) < control_bounds(1,1)) u(1,0) = control_bounds(1,1);
	}

	virtual void setupExperiments(experiments_t& experiments);
	void Dynamics::initControls(std::vector< Matrix<U_DIM> >& u) const {
		u.resize(horizon);
		for (int i = 0; i < horizon; ++i) {
			std::pair<double, double> samp;
			do {
				samp = normal();
			} while((abs(samp.first) > 1.0) || (abs(samp.second) > 1.0));
			u[i][0] = samp.first;
			u[i][1] = samp.second;
		}
	}

	virtual inline Matrix<X_DIM> f(double step, const Matrix<X_DIM>& x, const Matrix<U_DIM>& u) const;
	virtual inline SymmetricMatrix<X_DIM> varM(const Matrix<X_DIM>& x, const Matrix<U_DIM>& u) const;

	virtual inline void linearizeDynamics(double step, const Matrix<X_DIM>& xBar, const Matrix<U_DIM>& uBar, Matrix<X_DIM>& c, Matrix<X_DIM, X_DIM>& A, Matrix<X_DIM, U_DIM>& B, SymmetricMatrix<X_DIM>& M, unsigned int flag) const;

	virtual inline void linearizeObservation(double step, const Matrix<X_DIM>& xBar, const Matrix<U_DIM>& uBar, Matrix<Z_DIM, X_DIM>& H, SymmetricMatrix<Z_DIM>& N) const;

	virtual inline void quadratizeFinalCost(const Matrix<X_DIM>& xBar, const SymmetricMatrix<X_DIM>& SigmaBar, double& s, SymmetricMatrix<X_DIM>& S, Matrix<1,X_DIM>& sT, Matrix<1,S_DIM>& tT, unsigned int flag) const;
	virtual inline bool quadratizeCost(double step, const Matrix<X_DIM>& xBar, const SymmetricMatrix<X_DIM>& SigmaBar, const Matrix<U_DIM>& uBar, double& q, SymmetricMatrix<X_DIM>& Q, 
		SymmetricMatrix<U_DIM>& R, Matrix<U_DIM, X_DIM>& P, Matrix<1,X_DIM>& qT, Matrix<1,U_DIM>& rT, Matrix<1,S_DIM>& pT, unsigned int flag) const;

	inline Matrix<SIM_X_DIM> g(const Matrix<SIM_X_DIM>& x, const Matrix<U_DIM>& u, const double mass1, const double mass2, const double mass_center1, const double mass_center2, const double length1, const double length2, const double damping1, const double damping2, const double coulomb1, const double coulomb2, const double rotational1_inertia, const double rotational2_inertia, const double motor1_inertia, const double motor2_inertia, const double motor1_torque_const, const double motor2_torque_const, const double gear1_ratio, const double gear2_ratio) const {
		Matrix<SIM_X_DIM> xNew;

		double x0 = x[0];
		double x1 = x[1];
		double x2 = x[2];
		double x3 = x[3];
		double mass2_sq = mass2*mass2;
		double mass_center2_sq = mass_center2*mass_center2;
		double gear1_ratio_sq = gear1_ratio*gear1_ratio;
		double gear2_ratio_sq = gear2_ratio*gear2_ratio;
		double length1_sq = length1*length1;
		double length2_sq = length2*length2;
		double cos2minus1 = cos(x1 - x0);
		double cos2minus1_sq = cos2minus1*cos2minus1;

		double unknown_h = length1*mass_center2*mass2*sin(x1 - x0);

		double H11 = motor1_inertia*gear1_ratio_sq + rotational1_inertia + mass2*length1_sq;
		double H12 = length1*mass_center2*mass2*cos2minus1;
		double H22 = motor2_inertia*gear2_ratio_sq + rotational2_inertia;

		double D = H11*H22 - H12*H12;

		double v1 = gear1_ratio*motor1_torque_const*u[0]
		+ unknown_h*x3*x3
			- (mass_center1*mass1 + length1*mass2)*gravity*cos(x0)
			- (damping1*x2 + coulomb1*sgn(x2));
		double v2 = gear2_ratio*motor2_torque_const*u[1]
		- unknown_h*x2*x2
			- mass_center2*mass2*gravity*cos(x1)
			- (damping2*x3 + coulomb2*sgn(x3));

		xNew[0] = x2;
		xNew[1] = x3;
		xNew[2] = (H22*v1 - H12*v2)/D;
		xNew[3] = (-H12*v1 + H11*v2)/D;

		return xNew;
	}

	inline Matrix<SIM_X_DIM> freal(double step, const Matrix<SIM_X_DIM>& x, const Matrix<U_DIM>& u) const {
		// RK4 integration
		Matrix<SIM_X_DIM> k1, k2, k3, k4, xinit;

		xinit = x;

		k1 = g(xinit, u, mass1, mass2, mass_center1, mass_center2, length1, length2, damping1, damping2, coulomb1, coulomb2, rotational1_inertia, rotational2_inertia, motor1_inertia, motor2_inertia, motor1_torque_const, motor2_torque_const, gear1_ratio, gear2_ratio);
		k2 = g(xinit + 0.5*step*k1, u, mass1, mass2, mass_center1, mass_center2, length1, length2, damping1, damping2, coulomb1, coulomb2, rotational1_inertia, rotational2_inertia, motor1_inertia, motor2_inertia, motor1_torque_const, motor2_torque_const, gear1_ratio, gear2_ratio);
		k3 = g(xinit + 0.5*step*k2, u, mass1, mass2, mass_center1, mass_center2, length1, length2, damping1, damping2, coulomb1, coulomb2, rotational1_inertia, rotational2_inertia, motor1_inertia, motor2_inertia, motor1_torque_const, motor2_torque_const, gear1_ratio, gear2_ratio);
		k4 = g(xinit + step*k3, u, mass1, mass2, mass_center1, mass_center2, length1, length2, damping1, damping2, coulomb1, coulomb2, rotational1_inertia, rotational2_inertia, motor1_inertia, motor2_inertia, motor1_torque_const, motor2_torque_const, gear1_ratio, gear2_ratio);

		Matrix<SIM_X_DIM> xNew = xinit + step*(k1 + 2.0*(k2 + k3) + k4)/6.0;

		return xNew;
	}

	inline SymmetricMatrix<SIM_X_DIM> varMReal(const Matrix<SIM_X_DIM>& x, const Matrix<U_DIM>& u) const {
		Matrix<X_DIM> temp = zero<X_DIM>();
		temp.insert(0,0,x);
		SymmetricMatrix<SIM_X_DIM> temp2 = varM(temp, u).subSymmetricMatrix<SIM_X_DIM>(0);

		return temp2;
	}

	float gravity, gravity_est;

	double mass1, mass1_est;
	double mass_center1, mass_center1_est;
	double length1, length1_est;
	double damping1, damping1_est;
	double coulomb1, coulomb1_est;
	double rotational1_inertia, rotational1_inertia_est;
	double motor1_inertia, motor1_inertia_est;
	double motor1_torque_const, motor1_torque_const_est;
	double gear1_ratio, gear1_ratio_est;

	double mass2, mass2_est;
	double mass_center2, mass_center2_est;
	double length2, length2_est;
	double damping2, damping2_est;
	double coulomb2, coulomb2_est;
	double rotational2_inertia, rotational2_inertia_est;
	double motor2_inertia, motor2_inertia_est;
	double motor2_torque_const, motor2_torque_const_est;
	double gear2_ratio, gear2_ratio_est;
};

#include "common.h"

#endif // __POMDP_2DPLANARARM_H__
