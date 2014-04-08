#include "point-pf.h"
#include "casadi/casadi-point-pf.h"

#include <../util/matrix.h>
#include <../util/utils.h>
#include <../util/logging.h>

#include <iostream>

int main(int argc, char* argv[]) {
	srand(time(0));
	int M = 20; // number of particles
	point_pf::initialize(M);

	Matrix<X_DIM> x0, xGoal;
	x0[0] = -3.5; x0[1] = 2;
	xGoal[0] = -3.5; xGoal[1] = -2;

	SymmetricMatrix<X_DIM> Sigma0 = .01*identity<X_DIM>();

	Matrix<U_DIM> u = (xGoal - x0) / (DT*(T-1));

	std::vector<Matrix<X_DIM> > P0(M);
	for(int m=0; m < M; ++m) {
		P0[m] = sampleGaussian(x0, Sigma0);
	}

	std::vector<std::vector<Matrix<X_DIM>> > P_t(T);
	std::vector<Matrix<Q_DIM> > dyn_noise;
	std::vector<Matrix<R_DIM> > obs_noise;
	float sampling_noise;
	P_t[0] = P0;
	for(int t=0; t < T-1; ++t) {
		dyn_noise = sampleGaussianN(zeros<Q_DIM,1>(), Q, M);
		obs_noise = sampleGaussianN(zeros<R_DIM,1>(), R, M);
		sampling_noise = (1/float(M))*(rand() / float(RAND_MAX));
		//P_t[t+1] = point_pf::beliefDynamics(P_t[t], u, dyn_noise, obs_noise, sampling_noise);
		P_t[t+1] = point_pf::casadiBeliefDynamics(P_t[t], u, dyn_noise, obs_noise, sampling_noise);
	}

	for(int t=0; t < T; ++t) {
		std::cout << "\nt: " << t << "\n";
		for(int m=0; m < M; ++m) {
			std::cout << ~P_t[t][m];
		}
	}

	point_pf::pythonDisplayParticles(P_t);

	return 0;
}
