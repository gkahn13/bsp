#include "../include/planar-system.h"


void test_belief_dynamics() {
	vec<2> camera, object;
	camera << 0, 1;
	object << 5, 8;
	bool is_static = false;

	PlanarSystem sys = PlanarSystem(camera, object, is_static);
	int T = 10;

	// initial uncertainty about object is large
	mat<X_DIM,X_DIM> sigma0 = .01*mat<X_DIM,X_DIM>::Identity();
	sigma0.block<C_DIM,C_DIM>(J_DIM,J_DIM) = 20*mat<2,2>::Identity();

	vec<X_DIM> x0;
	x0 << M_PI/5, -M_PI/2, -M_PI/4, 0, 5, 5;
//	vec x0 = {M_PI/2, 0, 0, 0, 5, 5};

	double alpha = 1;

	// initialize state and controls
	std::vector<vec<U_DIM>, aligned_allocator<vec<U_DIM>>> U(T-1, vec<U_DIM>::Zero());
	std::vector<vec<X_DIM>, aligned_allocator<vec<X_DIM>>> X(T);
	std::vector<mat<X_DIM,X_DIM>, aligned_allocator<mat<X_DIM,X_DIM>>> S(T);
	X[0] = x0;
	S[0] = sigma0;
	for(int t=0; t < T-1; ++t) {
//		std::cout << "t: " << t << "\n";
//		X[t+1] = sys.dynfunc(X[t], U[t], zeros<vec>(U_DIM));
		sys.belief_dynamics(X[t], S[t], U[t], .1, X[t+1], S[t+1]);
//		std::cout << "S:\n" << S[t] << "\n\n";
//		std::cout << "trace(S): " << S[t].trace() << "\n";
	}

	double cost = sys.cost(X, sigma0, U, alpha);
	std::cout << "Cost: " << cost << "\n";

//	return;

	vec<TOTAL_VARS> grad = sys.cost_grad(X, sigma0, U, alpha);
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

//	sys.display(X, S);
}


int main(int argc, char* argv[]) {
	test_belief_dynamics();
}
