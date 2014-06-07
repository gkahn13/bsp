#include "../include/planar-system.h"

int main(int argc, char* argv[]) {
	vec camera = {0, 0};
	vec object = {1, 1};
	bool is_static = true;

	PlanarSystem sys = PlanarSystem(camera, object, is_static);

	int T = 5;
	int X_DIM = 6;
	std::vector<vec> X(T, zeros<vec>(X_DIM));

	sys.display(X);
}
