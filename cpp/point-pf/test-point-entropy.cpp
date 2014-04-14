#include "point-pf.h"
#include "point-entropy.h"

#include <../util/matrix.h>
#include <../util/utils.h>
#include <../util/logging.h>

#include <iostream>

int main(int argc, char* argv[]) {
	srand(time(0));
	point_pf::initialize();

	Matrix<X_DIM> x0, xGoal;
	x0[0] = -3.5; x0[1] = 2;
	xGoal[0] = -3.5; xGoal[1] = -2;

	SymmetricMatrix<X_DIM> Sigma0 = .01*identity<X_DIM>();

	Matrix<U_DIM> u = (xGoal - x0) / (DT*(T-1));
	std::vector<Matrix<U_DIM> > U(T-1, u);

	std::vector<Matrix<X_DIM> > P0(M);
	for(int m=0; m < M; ++m) {
		P0[m] = sampleGaussian(x0, Sigma0);
	}

	std::vector<std::vector<Matrix<X_DIM>> > P(T, std::vector<Matrix<X_DIM> >(M));
	std::vector<Matrix<Q_DIM> > dyn_noise;
	P[0] = P0;
	for(int t=0; t < T-1; ++t) {
		dyn_noise = sampleGaussianN(zeros<Q_DIM,1>(), 1e4*Q, M);
		for(int m=0; m < M; ++m) {
			P[t+1][m] = point_pf::dynfunc(P[t][m], u, dyn_noise[m]);
		}
	}

	float avg_entropy = 0;
	float max_iter = 1000;
	for(int iter=0; iter < max_iter; iter++) {
		float entropy = point_entropy::differential_entropy(P, U);
		std::cout << "entropy: " << entropy << "\n";
		avg_entropy += (1/float(max_iter))*entropy;
	}
	std::cout << "\nentropy: " << avg_entropy << "\n";

	point_pf::pythonDisplayParticles(P);

	return 0;
}
