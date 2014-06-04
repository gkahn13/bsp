#ifndef _EIH_SYSTEM_H__
#define _EIH_SYSTEM_H__

#include "../../system.h"
#include "pr2_sim.h"
#include "rave_utils.h"
#include "utils.h"

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
	void update_state_and_particles(const mat& x_t, const mat& P_t, const mat& u_t, mat& x_tp1, mat& P_tp1, bool plot);

	double cost(const std::vector<mat>& X, const std::vector<mat>& U, const mat& P);
	double cost(const mat &x0, const std::vector<mat>& U, const mat& P);
	mat cost_grad(const mat &x0, std::vector<mat>& U, const mat& P);

	void display_states_and_particles(const std::vector<mat>& X, const mat& P, bool pause=true);

	Manipulator* get_manip() { return manip; }
	KinectSensor* get_kinect() { return kinect; }

protected:
	void init(int Z_DIM, const mat &R, const mat &uMin, const mat &uMat, int T=10, double DT=1.0);

	double obsfunc_continuous_weight(const mat &particle, const cube &image, const mat &z_buffer, bool plot=false);
	double obsfunc_discrete_weight(const mat &particle, const cube &image, const mat &z_buffer, bool plot=false);

	rave::EnvironmentBasePtr env;
	Manipulator *manip;
	KinectSensor *kinect;

	std::vector<mat> desired_observations;

	// save handles. clear every time update_state_and_particles called
	std::vector<rave::GraphHandlePtr> handles;
};

#endif
