#include "../include/planar-system.h"

int main(int argc, char* argv[]) {
	vec camera = {0, .05};
	vec object = {.5, .5};
	bool is_static = true;

	PlanarSystem sys = PlanarSystem(camera, object, is_static);

	int T = 3;
	int X_DIM = 6;
	std::vector<vec> X(T, zeros<vec>(X_DIM));

	X[0] = {M_PI/4, M_PI/4, M_PI/4, 0, .2, .2};
	X[1] = {M_PI/8, M_PI/8, M_PI/8, 0, .2, 0};
	X[2] = {M_PI/16, M_PI/16, M_PI/16, 0, .2, -.2};

	std::vector<mat> S(T, .05*eye<mat>(X_DIM, X_DIM));

	sys.display(X, S);
}
