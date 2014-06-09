#include "../include/planar-system.h"

void test_display() {
	vec camera = {0, .5};
	vec object = {5, 7};
	bool is_static = true;

	PlanarSystem sys = PlanarSystem(camera, object, is_static);

	int T = 3;
	int X_DIM = 6;
	std::vector<vec> X(T, zeros<vec>(X_DIM));

	X[0] = {M_PI/4, M_PI/4, M_PI/4, 0, 5, 5};
	X[1] = {M_PI/8, M_PI/8, M_PI/8, 0, 0, 5};
	X[2] = {M_PI/16, M_PI/16, M_PI/16, 0, -5, 5};

	std::vector<mat> S(T, 1*eye<mat>(X_DIM, X_DIM));

	sys.display(X, S);
}

void test_beams() {
	vec camera = {0, .5};
	vec object = {5, 7};
	bool is_static = false;

	PlanarSystem sys = PlanarSystem(camera, object, is_static);

	int T = 1;
	int X_DIM = 6;
	std::vector<vec> X(T);
	std::vector<mat> S(T, eye<mat>(X_DIM, X_DIM));

	X[0] = {M_PI/3, -M_PI/2, -M_PI/16, M_PI/8, 5, 5};

	sys.display(X, S);
}

void test_sd() {
	vec camera = {0, .5};
	vec object = {5, 7};
	bool is_static = false;

	PlanarSystem sys = PlanarSystem(camera, object, is_static);

	int T = 1;
	int X_DIM = 6;
	std::vector<vec> X(T);
	std::vector<mat> S(T, eye<mat>(X_DIM, X_DIM));

	X[0] = {M_PI/3, -M_PI/2, -M_PI/16, M_PI/8, 5, 5};

	std::vector<Beam> fov = sys.get_fov(X[0]);

	std::cout << "number of beams: " << fov.size() << "\n\n";

	std::vector<Segment> border = geometry2d::beams_border(fov);

//	std::cout << "border\n";
//	for(int i=0; i < border.size(); ++i) {
//		std::cout << border[i].p0.t() << border[i].p1.t();
//		std::cout << "Length: " << border[i].length() << "\n\n";
//	}

	double sd;
	vec p;

	std::cout << "Randomized tests\n";
	for(int i=0; i < 5; ++i) {
		p << planar_utils::uniform(-5, 5) << planar_utils::uniform(0, 5);
		sd = geometry2d::signed_distance(p, fov);
		std::cout << "sd to (" << p(0) << ", " << p(1) << ")\n";
		std::cout << "sd: " << sd << "\n\n";
	}

	p << 5 << 5;
	sd = geometry2d::signed_distance(p, fov);
	std::cout << "sd to (" << p(0) << ", " << p(1) << ") should be > 0\n";
	std::cout << "sd: " << sd << "\n\n";

	p << .05 << 4;
	sd = geometry2d::signed_distance(p, fov);
	std::cout << "sd to (" << p(0) << ", " << p(1) << ") should be < 0\n";
	std::cout << "sd: " << sd << "\n\n";

	sys.display(X, S);

}

int main(int argc, char* argv[]) {
//	test_display();
//	test_beams();
	test_sd();
}
