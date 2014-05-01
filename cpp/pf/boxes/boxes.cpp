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
#define B_DIM 4
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
	M = 10;

	mat boxes(B_DIM, 1, fill::zeros);
	boxes << 0 << endr << 0 << endr << 1 << endr << 2;

	LOG_DEBUG("Initializing...");
	BoxesSystem sys = BoxesSystem(boxes, obs_type, cost_type, use_casadi,
			T, M, N, DT, B_DIM, X_DIM, U_DIM, Z_DIM, Q_DIM, R_DIM);
	LOG_DEBUG("System initialized");


	mat P(X_DIM, M, fill::zeros);
	for(int m=0; m < M; ++m) {
		for(int i=0; i < X_DIM; ++i) {
			P(i, m) = uniform(-1, 1);
		}
	}

	std::vector<mat> X(T, zeros<mat>(X_DIM, 1));

	sys.display_states_and_particles(X, P);
}
