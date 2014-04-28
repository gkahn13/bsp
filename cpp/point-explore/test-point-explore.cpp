//#include "point-explore.h"
#include "point-explore-system.h"

//#include <../util/matrix.h>
//#include <../util/utils.h>
#include <../util/logging.h>

#include <iostream>

#include <armadillo>

using namespace arma;

void test_update() {
//	srand(time(0));

	mat x0(N*X_DIM, 1, fill::zeros);
	x0 << 0 << endr << 0 << endr << .5 << endr << 0;

	mat target(X_DIM, 1);
	target << 3 << endr << 3;

	ObsType obs_type = ObsType::angle;
	CostType cost_type = CostType::platt;
	bool use_casadi = true;

	LOG_DEBUG("Initializing...");
	PointExploreSystem sys = PointExploreSystem(target, obs_type, cost_type, use_casadi);
	LOG_DEBUG("System initialized");

	mat xMin = sys.get_xMin();
	mat xMax = sys.get_xMax();

	const int M_FULL = 100;

	mat P0(X_DIM, M_FULL);
	for(int m=0; m < M_FULL; ++m) {
		for(int i=0; i < X_DIM; ++i) {
			P0(i,m) = uniform(0, 5);
		}
	}

	std::vector<mat> P(T, zeros<mat>(X_DIM, M_FULL));
	P[0] = P0;
	std::vector<mat> X(T, zeros<mat>(N*X_DIM,1));
	X[0] = x0;

//	LOG_DEBUG("Initial map");
//	sys.display_states_and_particles(std::vector<mat>(1, X[0]), P[0]);

	mat uinit(N*U_DIM, 1, fill::zeros);
	uinit << 0 << endr << .5 << endr << .5 << endr << 0;
//	uinit << .3 << endr << .3 << endr << .3 << endr << .3;
	std::vector<mat> U(T-1, uinit);

	double cost = sys.cost(X, U, P[0]);
//	LOG_DEBUG("Cost: %4.10f", cost);

	mat g = sys.cost_grad(X, U, P[0]);

	return;

	for(int t=0; t < T-1; ++t) {
		sys.update_state_and_particles(X[t], P[t], U[t], X[t+1], P[t+1]);

		sys.display_states_and_particles(std::vector<mat>(1,X[t+1]), P[t+1]);
	}

}

int main(int argc, char* argv[]) {
	test_update();
}
