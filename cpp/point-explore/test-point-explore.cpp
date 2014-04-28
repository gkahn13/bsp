//#include "point-explore.h"
#include "point-explore-system.h"

//#include <../util/matrix.h>
//#include <../util/utils.h>
#include <../util/logging.h>

#include <iostream>

#include <armadillo>

using namespace arma;

void test_update() {
	srand(time(0));

	mat x0(N*X_DIM, 1, fill::zeros);
	x0 << 0 << endr << 0 << endr << .5 << endr << 0;

	mat target(X_DIM, 1);
	target << 3 << endr << 3;

	ObsType obs_type = ObsType::distance;
	CostType cost_type = CostType::entropy;
	bool use_casadi = false;

	PointExploreSystem sys = PointExploreSystem(target, obs_type, cost_type, use_casadi);

	mat xMin = sys.get_xMin();
	mat xMax = sys.get_xMax();

	const int M_FULL = 1000;

	mat P0(X_DIM, M_FULL);
	for(int m=0; m < M_FULL; ++m) {
		for(int i=0; i < X_DIM; ++i) {
			P0(i,m) = uniform(xMin(i), xMax(i));
		}
	}

	std::vector<mat> P(T, zeros<mat>(X_DIM, M_FULL));
	P[0] = P0;
	std::vector<mat> X(T, zeros<mat>(N*X_DIM,1));
	X[0] = x0;

	std::cout << "Initial map\n";
	sys.display_states_and_particles(std::vector<mat>(1, X[0]), P[0]);

	mat uinit(N*U_DIM, 1, fill::zeros);
	//uinit << 0 << endr << .5 << endr << .5 << endr << 0;
	uinit << .3 << endr << .3 << endr << .3 << endr << .3;
	std::vector<mat> U(T-1, uinit);

	for(int t=0; t < T-1; ++t) {
		LOG_DEBUG("t: %d",t);
		sys.update_state_and_particles(X[t], P[t], U[t], X[t+1], P[t+1]);
		sys.display_states_and_particles(std::vector<mat>(1,X[t+1]), P[t+1]);
	}

}
//
//
//void test_entropy() {
////	srand(time(0));
//
//	R = .01*identity<N*R_DIM>();
//
//	xMin[0] = .45; xMin[1] = .45;
//	xMax[0] = .55; xMax[1] = .55;
//
//	x0[0] = 0; x0[1] = 0;
//	x0[2] = 4; x0[3] = 0;
//
//	target[0] = 3; target[1] = 3;
//
//	std::vector<Matrix<X_DIM> > P(M);
//	for(int m=0; m < M; ++m) {
//		if (m < 0) {
//			P[m][0] = uniform(.9, 1.1);
//			P[m][1] = uniform(-.1, .1);
//		} else {
//			P[m][0] = uniform(-.1, .1);
//			P[m][1] = uniform(.9, 1.1);
//		}
//	}
//
//	Matrix<N*U_DIM> uinit;
//	uinit[0] = 0.1; uinit[1] = 0;
//	uinit[2] = 0; uinit[3] = 0.1;
//	std::vector<Matrix<N*U_DIM> > U(T-1, uinit);
//
//	std::vector<Matrix<N*X_DIM> > X(T);
//	X[0] = x0;
//	for(int t=0; t < T-1; ++t) {
//		X[t+1] = point_explore::dynfunc(X[t], U[t]);
//	}
//
//	float entropy = point_explore::differential_entropy(X, U, P);
////	Matrix<TOTAL_VARS> grad_entropy = point_explore::grad_differential_entropy(X, U, P);
////
////	int index = 0;
////	for(int t=0; t < T; ++t) {
////		std::cout << "\nt: " << t << "\n";
////		std::cout << "x: ";
////		for(int i=0; i < X_DIM; ++i) {
////			std::cout << grad_entropy[index++] << " ";
////		}
////
////
////		if (t < T-1) {
////			std::cout << "\nu: ";
////			for(int i=0; i < U_DIM; ++i) {
////				std::cout << grad_entropy[index++] << " ";
////			}
////		}
////		std::cout << "\n";
////	}
//
//	std::cout << "entropy: " << entropy << "\n\n";
//	point_explore::pythonDisplayStatesAndParticles(std::vector<Matrix<N*X_DIM> >(1,x0), P, target);
//}

void test_system() {
	PointExploreSystem sys = PointExploreSystem();

	mat x(N*X_DIM,1,fill::ones), u(N*U_DIM,1,fill::zeros);
	mat x_tp1 = sys.dynfunc(x, u);

	x_tp1.rows(0,1) = 2*ones<mat>(2,1);
	x_tp1.print();
}

int main(int argc, char* argv[]) {
//	mat R = 1e-2*eye<mat>(N*R_DIM, N*R_DIM);
//	mat c = chol(R);
//	c.print();
	test_update();
//	test_entropy();
//	test_system();
}
