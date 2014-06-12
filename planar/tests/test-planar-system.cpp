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

//	X[0] = {M_PI/3, -M_PI/2, -M_PI/16, M_PI/8-M_PI/16, 5, 5};
	X[0] = {-M_PI/3, M_PI/2, M_PI/16, -M_PI/8+M_PI/16, 5, 5};

	std::vector<Beam> fov = sys.get_fov(X[0]);
	std::vector<Segment> link_segments = sys.get_link_segments(X[0]);

	std::cout << "Last link:\n";
	std::cout << link_segments.back().p0.t() << link_segments.back().p1.t() << "\n";

	for(int i=0; i < fov.size(); ++i) {
		std::cout << "i: " << i << "\n\n";

		std::cout << "original beam:\n";
		std::cout << fov[i].a.t();
		std::cout << fov[i].b.t() << "\n";

		std::vector<Beam> new_beams = fov[i].truncate(link_segments.back());
		std::cout << "new_beams:\n";
		for(int j=0; j < new_beams.size(); ++j) {
			std::cout << new_beams[i].a.t();
			std::cout << new_beams[i].b.t() << "\n";
		}
	}


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

void test_truncate_gaussian() {
//	vec camera = {0, .5};
	vec camera = {0, 0.01};
	vec object = {5, 7};
	bool is_static = false;

	PlanarSystem sys = PlanarSystem(camera, object, is_static);

	int X_DIM = 6;
	std::vector<vec> X_initial;
	std::vector<mat> S_initial;

	X_initial.push_back( vec{M_PI/3, -M_PI/2, -M_PI/16, M_PI/8, 2, 5} );
	S_initial.push_back(20*eye<mat>(X_DIM, X_DIM));

	X_initial.push_back( vec{M_PI/2.5, -M_PI/2, -M_PI/3, 0, 2, 5} );
	S_initial.push_back(20*eye<mat>(X_DIM, X_DIM));

	X_initial.push_back( vec{M_PI/3, -M_PI/2, -M_PI/16, -M_PI/8, 2, 5} );
	S_initial.push_back(20*eye<mat>(X_DIM, X_DIM));

	X_initial.push_back( vec{M_PI/2, 0, 0, 0, .05, 5} );
	S_initial.push_back(20*eye<mat>(X_DIM, X_DIM));

	X_initial.push_back( vec{M_PI/3, -M_PI/2, -M_PI/16, -M_PI/6, 2, 5} );
	S_initial.push_back(20*eye<mat>(X_DIM, X_DIM));


	for(int i=0; i < X_initial.size(); ++i) {
		std::cout << "Iteration: " << i << "\n";
		std::cout << "Before truncation\n";

		std::vector<vec> X(1, X_initial[i]);
		std::vector<mat> S(1, S_initial[i]);

		sys.display(X, S);

		std::vector<Beam> fov = sys.get_fov(X[0]);

		for(int t=0; t < 10; ++t) {
			vec cur_mean = X[0].subvec(4,5), out_mean;
			mat cur_cov = S[0].submat(span(4,5), span(4,5)), out_cov;

//			geometry2d::truncate_belief(fov, cur_mean, cur_cov, out_mean, out_cov);
			geometry2d::my_truncate_belief(fov, cur_mean, cur_cov, out_mean, out_cov);

			std::cout << "out_mean: " << out_mean.t();
			std::cout << "out_cov:\n" << out_cov;

			X[0].subvec(4,5) = out_mean;
			S[0].submat(span(4,5), span(4,5)) = out_cov;

			std::cout << "After truncation\n";
			sys.display(X, S);
		}
	}
}

void test_belief_dynamics() {


	vec camera = {0, 1};
	vec object = {5, 8};
	bool is_static = false;

	PlanarSystem sys = PlanarSystem(camera, object, is_static);
	int X_DIM = 6;
	int U_DIM = 4;
	int T = 10;

	// initial uncertainty about object is large
	mat sigma0 = .01*eye<mat>(X_DIM, X_DIM);
	sigma0.submat(span(4,5), span(4,5)) = 20*eye<mat>(2, 2);

	vec x0 = {M_PI/5, -M_PI/2, -M_PI/4, 0, 5, 5};
//	vec x0 = {M_PI/2, 0, 0, 0, 5, 5};

	double alpha = 1;

	// initialize state and controls
	std::vector<vec> U(T-1, zeros<vec>(U_DIM));
	std::vector<vec> X(T);
	std::vector<mat> S(T);
	X[0] = x0;
	S[0] = sigma0;
	for(int t=0; t < T-1; ++t) {
//		std::cout << "t: " << t << "\n";
//		X[t+1] = sys.dynfunc(X[t], U[t], zeros<vec>(U_DIM));
		sys.belief_dynamics(X[t], S[t], U[t], .1, X[t+1], S[t+1]);
//		std::cout << "S:\n" << S[t] << "\n\n";
//		std::cout << "trace(S): " << trace(S[t]) << "\n";
	}

	double cost = sys.cost(X, sigma0, U, alpha);
	std::cout << "Cost: " << cost << "\n";

//	return;

	vec grad = sys.cost_grad(X, sigma0, U, alpha);
	int index = 0;
	for(int t=0; t < T; ++t) {
		std::cout << "t: " << t << "\n";
		std::cout << grad.subvec(index, index+X_DIM-1).t();
		index += X_DIM;
		if (t < T-1) {
			std::cout << grad.subvec(index, index+U_DIM-1).t();
			index += U_DIM;
		}
	}

	sys.display(X, S);

}

int main(int argc, char* argv[]) {
//	test_display();
//	test_beams();
//	test_sd();
//	test_truncate_gaussian();
	test_belief_dynamics();
}
