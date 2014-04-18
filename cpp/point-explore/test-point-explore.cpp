#include "point-explore.h"

#include <../util/matrix.h>
#include <../util/utils.h>
#include <../util/logging.h>

#include <iostream>

void test_update() {
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

	Matrix<U_DIM> uinit;
	uinit[0] = .25; uinit[1] = .25;
	std::vector<Matrix<U_DIM> > U(T-1, uinit);

	std::vector<std::vector<Matrix<X_DIM>> > P(T);
	P[0] = P0;
	std::vector<Matrix<X_DIM> > X(T);
	X[0] = x0;
	for(int t=0; t < T-1; ++t) {
		point_explore::updateStateAndParticles(X[t], P[t], U[t], X[t+1], P[t+1]);
		point_explore::pythonDisplayStateAndParticles(X[t+1], P[t+1], target);
	}

}

float uniform(float low, float high) {
	return (high - low)*(rand() / float(RAND_MAX)) + low;
}

void test_entropy() {
//	srand(time(0));

	R = .01*identity<R_DIM>();

	xMin[0] = .45; xMin[1] = .45;
	xMax[0] = .55; xMax[1] = .55;

	x0[0] = 0; x0[1] = 0;
	target[0] = 3; target[1] = 3;

	std::vector<Matrix<X_DIM> > P(M);
	for(int m=0; m < M; ++m) {
		if (m < M/4) {
			P[m][0] = uniform(.9, 1.1);
			P[m][1] = uniform(-.1, .1);
		} else {
			P[m][0] = uniform(-.1, .1);
			P[m][1] = uniform(.9, 1.1);
		}
	}

	Matrix<U_DIM> uinit;
	uinit[0] = 0; uinit[1] = 0.1;
	std::vector<Matrix<U_DIM> > U(T-1, uinit);

	std::vector<Matrix<X_DIM> > X(T);
	X[0] = x0;
	for(int t=0; t < T-1; ++t) {
		X[t+1] = point_explore::dynfunc(X[t], U[t]);
	}

	float entropy = point_explore::greg_entropy(X, U, P);

	std::cout << "entropy: " << entropy << "\n\n";
	point_explore::pythonDisplayStateAndParticles(x0, P, target);
}

int main(int argc, char* argv[]) {
//	test_update();
	test_entropy();
	return 0;
}
