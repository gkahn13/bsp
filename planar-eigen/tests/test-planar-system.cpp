#include "../include/planar-system.h"


void test_belief_dynamics() {
//	vec camera = {0, 1};
//	vec object = {5, 8};
//	bool is_static = false;
//
//	PlanarSystem sys = PlanarSystem(camera, object, is_static);
//	int X_DIM = 6;
//	int U_DIM = 4;
//	int T = 10;
//
//	// initial uncertainty about object is large
//	mat sigma0 = .01*eye<mat>(X_DIM, X_DIM);
//	sigma0.submat(span(4,5), span(4,5)) = 20*eye<mat>(2, 2);
//
//	vec x0 = {M_PI/5, -M_PI/2, -M_PI/4, 0, 5, 5};
////	vec x0 = {M_PI/2, 0, 0, 0, 5, 5};
//
//	double alpha = 1;
//
//	// initialize state and controls
//	std::vector<vec> U(T-1, zeros<vec>(U_DIM));
//	std::vector<vec> X(T);
//	std::vector<mat> S(T);
//	X[0] = x0;
//	S[0] = sigma0;
//	for(int t=0; t < T-1; ++t) {
////		std::cout << "t: " << t << "\n";
////		X[t+1] = sys.dynfunc(X[t], U[t], zeros<vec>(U_DIM));
//		sys.belief_dynamics(X[t], S[t], U[t], .1, X[t+1], S[t+1]);
////		std::cout << "S:\n" << S[t] << "\n\n";
//	}
//
//	double cost = sys.cost(X, sigma0, U, alpha);
//	std::cout << "Cost: " << cost << "\n";
////
//	vec grad = sys.cost_grad(X, sigma0, U, alpha);
//	int index = 0;
//	for(int t=0; t < T; ++t) {
//		std::cout << "t: " << t << "\n";
//		std::cout << grad.subvec(index, index+X_DIM-1).t();
//		index += X_DIM;
//		if (t < T-1) {
//			std::cout << grad.subvec(index, index+U_DIM-1).t();
//			index += U_DIM;
//		}
//	}
//
//	sys.display(X, S);

}

int main(int argc, char* argv[]) {
	test_belief_dynamics();
}
