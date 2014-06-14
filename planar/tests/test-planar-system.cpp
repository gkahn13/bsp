#include "../include/planar-system.h"

void test_display() {
	vec<2> camera, object;
	camera << 0, .5;
	object << 5, 7;
	bool is_static = true;

	PlanarSystem sys = PlanarSystem(camera, object, is_static);

	int T = 3;
	std::vector<vec<X_DIM>, aligned_allocator<vec<X_DIM>>> X(T, vec<X_DIM>::Zero());

	X[0] << M_PI/4, M_PI/4, M_PI/4, 0, 5, 5;
	X[1] << M_PI/8, M_PI/8, M_PI/8, 0, 0, 5;
	X[2] << M_PI/16, M_PI/16, M_PI/16, 0, -5, 5;

	std::vector<mat<X_DIM,X_DIM>, aligned_allocator<mat<X_DIM,X_DIM>>> S(T, 1*mat<X_DIM,X_DIM>::Identity());

	sys.display(X, S);
}

//void test_beams() {
//	vec<2> camera, object;
//	camera << 0, .5;
//	object << 5, 7;
//	bool is_static = false;
//
//	PlanarSystem sys = PlanarSystem(camera, object, is_static);
//
//	int T = 1;
//	std::vector<vec<X_DIM>, aligned_allocator<vec<X_DIM>>> X(T);
//	std::vector<mat<X_DIM,X_DIM>, aligned_allocator<mat<X_DIM,X_DIM>>> S(T, mat<X_DIM,X_DIM>::Identity());
//
//	X[0] << -M_PI/3, M_PI/2, M_PI/16, -M_PI/8+M_PI/16, 5, 5;
//
//	std::vector<Beam> fov = sys.get_fov(X[0]);
//	std::vector<Segment> link_segments = sys.get_link_segments(X[0]);
//
//	std::cout << "Last link:\n";
//	std::cout << link_segments.back().p0.transpose() << "\n" << link_segments.back().p1.transpose() << "\n\n";
//
//	for(int i=0; i < fov.size(); ++i) {
//		std::cout << "i: " << i << "\n\n";
//
//		std::cout << "original beam:\n";
//		std::cout << fov[i].a.transpose() << "\n";
//		std::cout << fov[i].b.transpose() << "\n\n";
//
//		std::vector<Beam> new_beams = fov[i].truncate(link_segments.back());
//		std::cout << "new_beams:\n";
//		for(int j=0; j < new_beams.size(); ++j) {
//			std::cout << new_beams[i].a.transpose() << "\n";
//			std::cout << new_beams[i].b.transpose() << "\n\n";
//		}
//	}
//
//	sys.display(X, S);
//}

void test_sd() {
	vec<2> camera, object;
	camera << 0, .5;
	object << 5, 7;
	bool is_static = false;

	PlanarSystem sys = PlanarSystem(camera, object, is_static);

	int T = 1;
	std::vector<vec<X_DIM>, aligned_allocator<vec<X_DIM>>> X(T);
	std::vector<mat<X_DIM,X_DIM>, aligned_allocator<mat<X_DIM,X_DIM>>> S(T, mat<X_DIM,X_DIM>::Identity());

	X[0] << M_PI/3, -M_PI/2, -M_PI/16, M_PI/8, 5, 5;

	std::vector<Beam> fov = sys.get_fov(X[0]);

	std::cout << "number of beams: " << fov.size() << "\n\n";

	std::vector<Segment> border = geometry2d::beams_border(fov);

//	std::cout << "border\n";
//	for(int i=0; i < border.size(); ++i) {
//		std::cout << border[i].p0.t() << border[i].p1.t();
//		std::cout << "Length: " << border[i].length() << "\n\n";
//	}

	double sd;
	vec<2> p;

	std::cout << "Randomized tests\n";
	for(int i=0; i < 5; ++i) {
		p << planar_utils::uniform(-5, 5), planar_utils::uniform(0, 5);
		sd = geometry2d::signed_distance(p, fov);
		std::cout << "sd to (" << p(0) << ", " << p(1) << ")\n";
		std::cout << "sd: " << sd << "\n\n";
	}

	p << 5, 5;
	sd = geometry2d::signed_distance(p, fov);
	std::cout << "sd to (" << p(0) << ", " << p(1) << ") should be > 0\n";
	std::cout << "sd: " << sd << "\n\n";

	p << .05, 4;
	sd = geometry2d::signed_distance(p, fov);
	std::cout << "sd to (" << p(0) << ", " << p(1) << ") should be < 0\n";
	std::cout << "sd: " << sd << "\n\n";

	sys.display(X, S);

}

void test_truncate_gaussian() {
	vec<2> camera, object;
	camera << 0, 0.01;
	object << 5, 7; // 5, 7
	bool is_static = false;

	PlanarSystem sys = PlanarSystem(camera, object, is_static);

	std::vector<vec<X_DIM>, aligned_allocator<vec<X_DIM>>> X_initial;
	std::vector<mat<X_DIM,X_DIM>, aligned_allocator<mat<X_DIM,X_DIM>>> S_initial;
	vec<X_DIM> tmp;

	tmp << M_PI/3, -M_PI/2, -M_PI/16, M_PI/8, 2, 5;
	X_initial.push_back(tmp);
	S_initial.push_back(20*mat<X_DIM,X_DIM>::Identity());

	tmp << M_PI/2.5, -M_PI/2, -M_PI/3, 0, 2, 5;
	X_initial.push_back(tmp);
	S_initial.push_back(20*mat<X_DIM,X_DIM>::Identity());

	tmp << M_PI/3, -M_PI/2, -M_PI/16, -M_PI/8, 2, 5;
	X_initial.push_back(tmp);
	S_initial.push_back(20*mat<X_DIM,X_DIM>::Identity());

	tmp << M_PI/2, 0, 0, M_PI/8, 4, 7;
	X_initial.push_back(tmp);
	S_initial.push_back(20*mat<X_DIM,X_DIM>::Identity());

	tmp << M_PI/3, -M_PI/2, -M_PI/16, -M_PI/6, 2, 5;
	X_initial.push_back(tmp);
	S_initial.push_back(20*mat<X_DIM,X_DIM>::Identity());


	for(int i=0; i < X_initial.size(); ++i) {
		std::cout << "Iteration: " << i << "\n";
		std::cout << "Before truncation\n";

		std::vector<vec<X_DIM>, aligned_allocator<vec<X_DIM>>> X(1, X_initial[i]);
		std::vector<mat<X_DIM,X_DIM>, aligned_allocator<mat<X_DIM,X_DIM>>> S(1, S_initial[i]);

		sys.display(X, S);

		std::vector<Beam> fov = sys.get_fov(X[0]);
		bool received_obs = geometry2d::is_inside(object, fov);
		std::cout << "Received observation: " << received_obs << "\n";

		for(int t=0; t < 10; ++t) {
			vec<2> cur_mean = X[0].segment<2>(4), out_mean;
			mat<2,2> cur_cov = S[0].block<2,2>(4,4), out_cov;

//			geometry2d::truncate_belief(fov, cur_mean, cur_cov, out_mean, out_cov);
			geometry2d::my_truncate_belief(fov, cur_mean, cur_cov, received_obs, out_mean, out_cov);

//			std::cout << "out_mean: " << out_mean.transpose() << "\n";
//			std::cout << "out_cov:\n" << out_cov << "\n";

			X[0].segment<2>(4) = out_mean;
			S[0].block<2,2>(4,4) = out_cov;

//			std::cout << "After truncation\n";
			sys.display(X, S);
		}
	}
}


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
	std::vector<vec<U_DIM>, aligned_allocator<vec<U_DIM>>> U(T-1, {.05, -.03, .01, -.02});
	U[T-2] = vec<U_DIM>::Zero();
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

void test_ik() {
	vec<C_DIM> camera, object;
	camera << 0, 1;
	object << 5, 8;
	bool is_static = false;

	PlanarSystem sys = PlanarSystem(camera, object, is_static);

	vec<C_DIM> goal_pos = {5, 5};
	vec<E_DIM> j;
	bool ik_success = sys.ik(goal_pos, j);

	std::cout << "ik_success: " << ik_success << "\n";
	std::cout << "joints: " << j.transpose() << "\n";
	std::cout << "ik pos soln: " << sys.get_ee_pos(j) << "\n";
}


int main(int argc, char* argv[]) {
//	test_display();
//	test_beams();
//	test_sd();
	test_truncate_gaussian();
//	test_belief_dynamics();
//	test_ik();
}
