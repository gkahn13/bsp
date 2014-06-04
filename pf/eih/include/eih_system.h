#pragma once

#include "../../system.h"
#include "pr2_sim.h"
#include "rave_utils.h"

#include <armadillo>
using namespace arma;

#include <openrave-core.h>
namespace rave = OpenRAVE;

const double UNKNOWN = -1e6;

template <class T>
inline T sigmoid(const T &x, const double alpha) {
	return 1/(1 + exp(-alpha*x));
}

class EihSystem : public virtual System {
public:
	EihSystem(rave::EnvironmentBasePtr e, Manipulator *m, KinectSensor *k);
	EihSystem(rave::EnvironmentBasePtr e, Manipulator *m, KinectSensor *k,
			int T, mat &R, mat &uMin, mat &uMax);

	mat dynfunc(const mat &x, const mat &u);
	mat obsfunc(const mat& x, const mat& t, const mat& r);

	void update_state_and_particles(const mat& x_t, const mat& P_t, const mat& u_t, mat& x_tp1, mat& P_tp1);
	double cost(const std::vector<mat>& X, const std::vector<mat>& U, const mat& P);

	void display_states_and_particles(const std::vector<mat>& X, const mat& P, bool pause=true);

	Manipulator* get_manip() { return manip; }
	KinectSensor* get_kinect() { return kinect; }

protected:
	void init(int Z_DIM, const mat &R, const mat &uMin, const mat &uMat, int T=10, double DT=1.0);

	double obsfunc_continuous_weight(const mat &particle, const cube &image, const mat &z_buffer);
	double obsfunc_discrete_weight(const mat &particle, const cube &image, const mat &z_buffer);

	rave::EnvironmentBasePtr env;
	Manipulator *manip;
	KinectSensor *kinect;

	std::vector<mat> desired_observations;

	// save handles. clear every time update_state_and_particles called
	std::vector<rave::GraphHandlePtr> handles;
};


//class System {
//public:
//	System();
//
//	virtual mat dynfunc(const mat& x, const mat& u) =0;
//	virtual mat obsfunc(const mat& x, const mat& t, const mat& r) =0;
//
//	virtual void update_state_and_particles(const mat& x_t, const mat& P_t, const mat& u_t, mat& x_tp1, mat& P_tp1) =0;
//	virtual double cost(const std::vector<mat>& X, const std::vector<mat>& U, const mat& P) =0;
//	mat cost_grad(std::vector<mat>& X, std::vector<mat>& U, const mat& P);
//
//	virtual void display_states_and_particles(const std::vector<mat>& X, const mat& P, bool pause=true) =0;
//
//	mat get_xMin() { return this->xMin; }
//	mat get_xMax() { return this->xMax; }
//	mat get_uMin() { return this->uMin; }
//	mat get_uMax() { return this->uMax; }
//
//protected:
//	int T, M, TOTAL_VARS, X_DIM, U_DIM, Z_DIM, Q_DIM, R_DIM;
//	double DT;
//
//	bool use_casadi;
//
//	mat R;
//	mat xMin, xMax, uMin, uMax;
//
//	CasadiSystem* casadi_sys;
//
//	double gauss_likelihood(const mat& v, const mat& S);
//	mat low_variance_sampler(const mat& P, const mat& W, double r);
//
//	double cost_entropy(const std::vector<mat>& X, const std::vector<mat>& U, const mat& P);
//	double cost_platt(const std::vector<mat>& X, const std::vector<mat>& U, const mat& P);
//};
