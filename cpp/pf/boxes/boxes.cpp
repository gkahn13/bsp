#include "boxes-system.h"

#include <iostream>
#include <vector>

#include "../../util/Timer.h"
#include "../../util/logging.h"

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#define TIMESTEPS 10
#define BOXES 1
#define DT 1.0 // Note: if you change this, must change the FORCES matlab file
#define X_DIM 2
#define U_DIM 2
#define Z_DIM 1
#define Q_DIM 2
#define R_DIM 1

const int T = TIMESTEPS;
const int N = 1;
int M;
const int TOTAL_VARS = T*N*X_DIM + (T-1)*N*U_DIM;


mat::fixed<X_DIM, 1> x0;

namespace cfg {
const double improve_ratio_threshold = .1;
const double min_approx_improve = 1e-4;
const double min_trust_box_size = 1e-4;
const double trust_shrink_ratio = .5;
const double trust_expand_ratio = 1.5;
}

int main(int argc, char* argv[]) {
	srand(time(0));

	ObsType obs_type = ObsType::distance;
	CostType cost_type = CostType::entropy;
	bool use_casadi = true;
	M = 1000;

	mat box_centers(N*X_DIM, 1, fill::zeros);
	mat box_dims(N*X_DIM, 1, fill::zeros);
	box_centers << -0.5 << endr << -.25;
	box_dims << .5 << endr << 1;

	LOG_DEBUG("Initializing...");
	BoxesSystem sys = BoxesSystem(box_centers, box_dims, obs_type, cost_type, use_casadi,
			T, M, N, DT, X_DIM, U_DIM, Z_DIM, Q_DIM, R_DIM);
	LOG_DEBUG("System initialized");


	mat P0(X_DIM, M, fill::zeros);
	for(int m=0; m < M; ++m) {
		P0(0, m) = uniform(-2, 0.5);
		P0(1, m) = uniform(-2, 2);
	}

	x0 << 1.5 << endr << 1.5;
	std::vector<mat> X(T, zeros<mat>(X_DIM, 1));

	mat u(U_DIM, 1, fill::zeros);
	u(0) = 0;
	u(1) = -.2;

	std::vector<mat> P(T, zeros<mat>(N*X_DIM, M));

	X[0] = x0;
	P[0] = P0;
	for(int t=0; t < T-1; ++t) {
//		X[t+1] = sys.dynfunc(X[t], u);
		sys.update_state_and_particles(X[t], P[t], u, X[t+1], P[t+1]);

		sys.display_states_and_particles(std::vector<mat>(1, X[t]), P[t]);
	}

}
