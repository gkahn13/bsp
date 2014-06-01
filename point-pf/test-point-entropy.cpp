#include "point-pf.h"
#include "point-entropy.h"

#include <../util/matrix.h>
#include <../util/utils.h>
#include <../util/logging.h>

#include <iostream>

int main(int argc, char* argv[]) {
	//srand(time(0));
	point_pf::initialize();

	Matrix<X_DIM> x0, xGoal;
	x0[0] = -3.5; x0[1] = 2;
	xGoal[0] = -3.5; xGoal[1] = -2;

	SymmetricMatrix<X_DIM> Sigma0 = .01*identity<X_DIM>();

	Matrix<U_DIM> u = (xGoal - x0) / (DT*(T-1));
//	u[0] = .35;
	std::vector<Matrix<U_DIM> > U(T-1, u);

	std::vector<Matrix<X_DIM> > P0(M);
	for(int m=0; m < M; ++m) {
		P0[m] = sampleGaussian(x0, Sigma0);
	}

	std::vector<std::vector<Matrix<X_DIM>> > P(T, std::vector<Matrix<X_DIM> >(M));
	std::vector<Matrix<Q_DIM> > dyn_noise;
	P[0] = P0;
	for(int t=0; t < T-1; ++t) {
		dyn_noise = sampleGaussianN(zeros<Q_DIM,1>(), Q, M);
		for(int m=0; m < M; ++m) {
			P[t+1][m] = point_pf::dynfunc(P[t][m], u, dyn_noise[m]);
		}
	}

	float avg_entropy = point_entropy::differential_entropy(P, U);
	std::cout << "\nentropy: " << avg_entropy << "\n";

	Matrix<TOTAL_VARS> grad = point_entropy::grad_differential_entropy(P, U);
	int index = 0;
	for(int t=0; t < T; ++t) {
		std::cout << "\n\nt: " << t << "\n";
		std::cout << "x: ";
		for(int i=0; i < M*X_DIM; ++i) {
			std::cout << grad[index++] << " ";
		}
		if (t < T-1) {
			std::cout << "\nu: ";
			for(int i=0; i < U_DIM; ++i) {
				std::cout << grad[index++] << " ";
			}
		}
	}

	point_pf::pythonDisplayParticles(P);

	return 0;
}
