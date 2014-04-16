#include "point-explore.h"

#include <../util/matrix.h>
#include <../util/utils.h>
#include <../util/logging.h>

#include <iostream>

int main(int argc, char* argv[]) {
	srand(time(0));

	R = .01*identity<R_DIM>();

	xMin[0] = 0; xMin[1] = 0;
	xMax[0] = 5; xMax[1] = 5;

	x0[0] = 0; x0[1] = 0;
	target[0] = 3; target[1] = 3;

	std::vector<Matrix<X_DIM> > P0(M);
	for(int m=0; m < M; ++m) {
		for(int i=0; i < X_DIM; ++i) {
			P0[m][i] = (xMax[i] - xMin[i])*(rand() / float(RAND_MAX)) + xMin[i];
		}
	}

	std::cout << "Initial map\n";
	point_explore::pythonDisplayStateAndParticles(x0, P0, target);

	Matrix<U_DIM> u;
	u[0] = .25; u[1] = .25;

	std::vector<std::vector<Matrix<X_DIM>> > P(T);
	P[0] = P0;
	std::vector<Matrix<X_DIM> > X(T);
	X[0] = x0;
	for(int t=0; t < T-1; ++t) {
		point_explore::updateStateAndParticles(X[t], P[t], u, X[t+1], P[t+1]);
		point_explore::pythonDisplayStateAndParticles(X[t+1], P[t+1], target);
	}

	return 0;
}
