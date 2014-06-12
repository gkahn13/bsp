#include "../include/planar-system.h"


void test_belief_dynamics() {
	VectorXd camera(C_DIM), object(C_DIM);
	camera << 0, 1;
	object << 5, 8;
	bool is_static = false;

	PlanarSystem sys = PlanarSystem(camera, object, is_static);
	int T = 10;

	// initial uncertainty about object is large
	MatrixXd sigma0 = .01*MatrixXd::Identity(X_DIM, X_DIM);
	sigma0.block<C_DIM,C_DIM>(J_DIM,J_DIM) = 20*MatrixXd::Identity(C_DIM,C_DIM);

	VectorXd x0(X_DIM);
	x0 << M_PI/5, -M_PI/2, -M_PI/4, 0, 5, 5;
//	vec x0 = {M_PI/2, 0, 0, 0, 5, 5};

	double alpha = 1;

	// initialize state and controls
	std::vector<VectorXd> U(T-1, VectorXd::Zero(U_DIM));
	std::vector<VectorXd> X(T);
	std::vector<MatrixXd> S(T);
	X[0] = x0;
	S[0] = sigma0;
	for(int t=0; t < T-1; ++t) {
//		std::cout << "t: " << t << "\n";
//		X[t+1] = sys.dynfunc(X[t], U[t], zeros<vec>(U_DIM));
		sys.belief_dynamics(X[t], S[t], U[t], .1, X[t+1], S[t+1]);
//		std::cout << "S:\n" << S[t] << "\n\n";
	}

	double cost = sys.cost(X, sigma0, U, alpha);
	std::cout << "Cost: " << cost << "\n";
//
	VectorXd grad = sys.cost_grad(X, sigma0, U, alpha);
	int index = 0;
	for(int t=0; t < T; ++t) {
		std::cout << "t: " << t << "\n";
		std::cout << grad.segment<X_DIM>(index).transpose() << "\n";
		index += X_DIM;
		if (t < T-1) {
			std::cout << grad.segment<U_DIM>(index).transpose() << "\n";
			index += U_DIM;
		}
	}

	sys.display(X, S);
}

int main(int argc, char* argv[]) {
	test_belief_dynamics();
}
