#include "system/pr2_eih_system.h"
#include "../util/logging.h"
#include "../util/Timer.h"

#include <stdio.h>

#include <boost/preprocessor/iteration/local.hpp>

extern "C" {
#include "pr2eihMPC.h"
pr2eihMPC_FLOAT **H, **f, **lb, **ub, **z, *c;
}

const int T = TIMESTEPS;

namespace cfg {
const double alpha_init = .01; // .01
const double alpha_gain = 1.5; // 3
const double alpha_epsilon = .1; // .1
const double alpha_max_increases = 5; // 5

const double Xeps_initial = .1; // .1
const double Ueps_initial = .1; // .1
const double improve_ratio_threshold = 1; // .1
const double min_approx_improve = 1e-1; // .1
const double min_trust_box_size = .1; // .1
const double trust_shrink_ratio = .75; // .75
const double trust_expand_ratio = 1.1; // 1.25
}

void setup_mpc_vars(pr2eihMPC_params& problem, pr2eihMPC_output& output) {
	// inputs
	H = new pr2eihMPC_FLOAT*[T];
	f = new pr2eihMPC_FLOAT*[T];
	lb = new pr2eihMPC_FLOAT*[T];
	ub = new pr2eihMPC_FLOAT*[T];
	c = problem.c1;

	// output
	z = new pr2eihMPC_FLOAT*[T];

#define SET_VARS(n)    \
		H[ BOOST_PP_SUB(n,1) ] = problem.H##n ;  \
		f[ BOOST_PP_SUB(n,1) ] = problem.f##n ;  \
		lb[ BOOST_PP_SUB(n,1) ] = problem.lb##n ;	\
		ub[ BOOST_PP_SUB(n,1) ] = problem.ub##n ;	\
		z[ BOOST_PP_SUB(n,1) ] = output.z##n ;

#define BOOST_PP_LOCAL_MACRO(n) SET_VARS(n)
#define BOOST_PP_LOCAL_LIMITS (1, TIMESTEPS)
#include BOOST_PP_LOCAL_ITERATE()

	for(int i=0; i < J_DIM; ++i) { c[i] = INFINITY; }

	for(int t=0; t < T-1; ++t) {
		for(int i=0; i < (J_DIM+U_DIM); ++i) { H[t][i] = INFINITY; }
		for(int i=0; i < (J_DIM+U_DIM); ++i) { f[t][i] = INFINITY; }
		for(int i=0; i < (J_DIM+U_DIM); ++i) { lb[t][i] = INFINITY; }
		for(int i=0; i < (J_DIM+U_DIM); ++i) { ub[t][i] = INFINITY; }
		for(int i=0; i < (J_DIM+U_DIM); ++i) { z[t][i] = INFINITY; }
	}
	for(int i=0; i < J_DIM; ++i) { H[T-1][i] = INFINITY; }
	for(int i=0; i < J_DIM; ++i) { f[T-1][i] = INFINITY; }
	for(int i=0; i < J_DIM; ++i) { lb[T-1][i] = INFINITY; }
	for(int i=0; i < J_DIM; ++i) { ub[T-1][i] = INFINITY; }
	for(int i=0; i < J_DIM; ++i) { z[T-1][i] = INFINITY; }

}

void cleanup_mpc_vars() {
	delete[] H;
	delete[] f;
	delete[] lb;
	delete[] ub;
	delete[] z;
	delete c;
}

bool is_valid_inputs()
{
//	for(int t = 0; t < T-1; ++t) {
//		std::cout << "\n\nt: " << t << "\n";
//
//		if (t == 0) {
//			std::cout << "\nc[0]:\n";
//			for(int i=0; i < (J_DIM); ++i) {
//				std::cout << c[i] << " ";
//			}
//		}
//
//		std::cout << "\nH[" << t << "]: ";
//		for(int i=0; i < (J_DIM+U_DIM); ++i) {
//			std::cout << H[t][i] << " ";
//		}
//
//		std::cout << "\nf[" << t << "]: ";
//		for(int i=0; i < (J_DIM+U_DIM); ++i) {
//			std::cout << f[t][i] << " ";
//		}
//
//		std::cout << "\nlb[" << t << "]: ";
//		for(int i=0; i < (J_DIM+U_DIM); ++i) {
//			std::cout << lb[t][i] << " ";
//		}
//
//		std::cout << "\nub[" << t << "]: ";
//		for(int i=0; i < (J_DIM+U_DIM); ++i) {
//			std::cout << ub[t][i] << " ";
//		}
//	}
//	std::cout << "\n\nt: " << T-1 << "\n";
//
//	std::cout << "\nH[" << T-1 << "]: ";
//	for(int i=0; i < (J_DIM); ++i) {
//		std::cout << H[T-1][i] << " ";
//	}
//
//	std::cout << "\nf[" << T-1 << "]: ";
//	for(int i=0; i < (J_DIM); ++i) {
//		std::cout << f[T-1][i] << " ";
//	}
//
//	std::cout << "\nlb[" << T-1 << "]: ";
//	for(int i=0; i < (J_DIM); ++i) {
//		std::cout << lb[T-1][i] << " ";
//	}
//
//	std::cout << "\nub[" << T-1 << "]: ";
//	for(int i=0; i < (J_DIM); ++i) {
//		std::cout << ub[T-1][i] << " ";
//	}
//
//	std::cout << "\nA[" << T-1 << "]:\n";
//	for(int i=0; i < 2*G_DIM; ++i) {
//		for(int j=0; j < J_DIM; ++j) {
//			std::cout << A[i+j*(2*G_DIM)] << " ";
//		}
//		std::cout << "\n";
//	}
//	std::cout << "\n";
//
//	std::cout << "b[" << T-1 << "]: ";
//	for(int i=0; i < 2*G_DIM; ++i) {
//		std::cout << b[i] << " ";
//	}
//
//	std::cout << "\n\n";

	for(int i=0; i < (J_DIM); ++i) { if (c[i] > INFINITY/2) { LOG_ERROR("isValid0"); return false; } }
	for(int i=0; i < (J_DIM); ++i) { if (c[i] < lb[0][i]) { LOG_ERROR("isValid1"); return false; } }
	for(int i=0; i < (J_DIM); ++i) { if (c[i] > ub[0][i]) { LOG_ERROR("isValid2"); return false; } }

	for(int t = 0; t < T-1; ++t) {
		for(int i=0; i < (J_DIM+U_DIM); ++i) { if (H[t][i] > INFINITY/2) { LOG_ERROR("isValid5"); return false; } }
		for(int i=0; i < (J_DIM+U_DIM); ++i) { if (f[t][i] > INFINITY/2) { LOG_ERROR("isValid6"); return false; } }
		for(int i=0; i < (J_DIM+U_DIM); ++i) { if (lb[t][i] > INFINITY/2) { LOG_ERROR("isValid7"); return false; } }
		for(int i=0; i < (J_DIM+U_DIM); ++i) {if (ub[t][i] > INFINITY/2) { LOG_ERROR("isValid8"); return false; } }
		for(int i=0; i < (J_DIM+U_DIM); ++i) {if (ub[t][i] < lb[t][i]) { LOG_ERROR("isValid9"); return false; } }
	}
	for(int i=0; i < (J_DIM); ++i) { if (H[T-1][i] > INFINITY/2) { LOG_ERROR("isValid10"); return false; } }
	for(int i=0; i < (J_DIM); ++i) { if (f[T-1][i] > INFINITY/2) { LOG_ERROR("isValid11"); return false; } }
	for(int i=0; i < (J_DIM); ++i) { if (lb[T-1][i] > INFINITY/2) { LOG_ERROR("isValid12"); return false; } }
	for(int i=0; i < (J_DIM); ++i) { if (ub[T-1][i] > INFINITY/2) { LOG_ERROR("isValid13"); return false; } }
	for(int i=0; i < (J_DIM); ++i) {if (ub[T-1][i] < lb[T-1][i]) { LOG_ERROR("isValid14"); return false; } }

	return true;
}

void L_BFGS(const StdVectorJ& J, const StdVectorU& U, const VectorTOTAL &grad,
		const StdVectorJ &Jopt, const StdVectorU &Uopt, const VectorTOTAL &gradopt,
		MatrixTOTAL &hess) {
	VectorTOTAL s = VectorTOTAL::Zero();

	int index = 0;
	for(int t=0; t < T-1; ++t) {
		s.segment<J_DIM>(index) = Jopt[t] - J[t];
		index += J_DIM;
		s.segment<U_DIM>(index) = Uopt[t] - U[t];
		index += U_DIM;
	}
	s.segment<J_DIM>(index) = Jopt[T-1] - J[T-1];

	VectorTOTAL y = gradopt - grad;

	double theta;
	VectorTOTAL hess_s = hess*s;

	bool decision = s.dot(y) >= .2*s.dot(hess_s);
	if (decision) {
		theta = 1;
	} else {
		theta = (.8*s.dot(hess_s))/(s.dot(hess_s) - s.dot(y));
	}

	VectorTOTAL r = theta*y + (1-theta)*hess_s;

	hess = hess - (hess_s*hess_s.transpose())/(s.dot(hess_s)) + (r*r.transpose())/s.dot(r);
}

double pr2_eih_collocation(StdVectorJ& J, StdVectorU& U, const MatrixP& P,
		const std::vector<geometry3d::Triangle>& obstacles, PR2EihSystem& sys,
		pr2eihMPC_params &problem, pr2eihMPC_output &output, pr2eihMPC_info &info) {
	int max_iter = 100;
	double Xeps = cfg::Xeps_initial;
	double Ueps = cfg::Ueps_initial;

	double merit = 0, meritopt = 0;
	double constant_cost, hessian_constant, jac_constant;
	VectorTOTAL grad = VectorTOTAL::Zero();
	MatrixTOTAL hess = MatrixTOTAL::Identity();

	StdVectorJ Jopt(T, VectorJ::Zero());
	StdVectorU Uopt(T-1, VectorU::Zero());
	VectorTOTAL gradopt = VectorTOTAL::Zero();
	double optcost, model_merit, new_merit;
	double approx_merit_improve, exact_merit_improve, merit_improve_ratio;

	LOG_DEBUG("Initial trajectory entropy: %4.10f", sys.entropy(J, U, P, obstacles));

	int index = 0;
	bool solution_accepted = true;
	for(int it=0; it < max_iter; ++it) {

		LOG_DEBUG(" ");
		LOG_DEBUG(" ");
		LOG_DEBUG("Iter: %d", it);

		// only compute gradient/hessian if P/U has been changed
		if (solution_accepted) {

			if (it == 0) {
				merit = sys.entropy(J, U, P, obstacles);
				grad = sys.entropy_grad(J, U, P, obstacles);
			} else {
				merit = meritopt; // since L-BFGS calculation required it
				grad = gradopt;
			}

//			std::cout << "gradient:\n" << grad << "\n";

			VectorTOTAL diaghess = hess.diagonal();

			constant_cost = 0;
			hessian_constant = 0;
			jac_constant = 0;

			// fill in Hessian first so we can force it to be PSD
			index = 0;
			for(int t=0; t < T-1; ++t) {
				for(int i=0; i < (J_DIM+U_DIM); ++i) {
					double val = diaghess(index++);
					H[t][i] = (val < 0) ? 0 : val;
				}
			}
			for(int i=0; i < J_DIM; ++i) {
				double val = diaghess(index++);
				H[T-1][i] = (val < 0) ? 0 : val;
			}

			// fill in gradient
			index = 0;
			for(int t=0; t < T-1; ++t) {
				Matrix<double,(J_DIM+U_DIM),1> zbar;
				zbar.segment<J_DIM>(0) = J[t];
				zbar.segment<U_DIM>(J_DIM) = U[t];

				for(int i=0; i < (J_DIM+U_DIM); ++i) {
					hessian_constant += H[t][i]*zbar(i)*zbar(i);
					jac_constant -= grad(index)*zbar(i);
					f[t][i] = grad(index) - H[t][i]*zbar(i);
					index++;
				}
			}

			VectorJ zbar = J[T-1];

			for(int i=0; i < J_DIM; ++i) {
				hessian_constant += H[T-1][i]*zbar(i)*zbar(i);
				jac_constant -= grad(index)*zbar(i);
				f[T-1][i] = grad(index) - H[T-1][i]*zbar(i);
				index++;
			}

			for(int i=0; i < J_DIM; ++i) {
				c[i] = J[0](i);
			}

			constant_cost = 0.5*hessian_constant + jac_constant + merit;
		}


		VectorJ j_min, j_max;
		VectorU u_min, u_max;
		sys.get_limits(j_min, j_max, u_min, u_max);

		const double epsilon = 1e-5;
		// set trust region bounds based on current trust region size
		for(int t=0; t < T; ++t) {
			index = 0;
			for(int i=0; i < J_DIM; ++i) {
				lb[t][index] = std::max(j_min(i), J[t](i) - Xeps);
				ub[t][index] = std::min(j_max(i), J[t](i) + Xeps);
				if (ub[t][index] < lb[t][index]) {
					ub[t][index] = lb[t][index] + epsilon;
				}
				index++;
			}


			if (t < T-1) {
				// set each input lower/upper bound
				for(int i=0; i < U_DIM; ++i) {
					lb[t][index] = std::max(u_min(i), U[t](i) - Ueps);
					ub[t][index] = std::min(u_max(i), U[t](i) + Ueps);
					if (ub[t][index] < lb[t][index]) {
						ub[t][index] = lb[t][index] + epsilon;
					}
					index++;
				}
			}
		}

		// Verify problem inputs
		if (!is_valid_inputs()) {
			LOG_ERROR("Inputs are not valid!");
			exit(0);
		}

		// call FORCES
		int exitflag = pr2eihMPC_solve(&problem, &output, &info);
		if (exitflag == 1) {
			optcost = info.pobj;
			for(int t=0; t < T; ++t) {
				index = 0;
				for(int i=0; i < J_DIM; ++i) {
					Jopt[t](i) = z[t][index++];
				}

				if (t < T-1) {
					for(int i=0; i < U_DIM; ++i) {
						Uopt[t](i) = z[t][index++];
					}
				}
			}
		} else {
			LOG_FATAL("Some problem in solver");
			LOG_FATAL("Continuing");
			return INFINITY;
		}

//		LOG_INFO("Plotting Jopt");
//		sys.plot(Jopt, P, obstacles, false);

		model_merit = optcost + constant_cost; // need to add constant terms that were dropped

		new_merit = sys.entropy(Jopt, Uopt, P, obstacles);

		LOG_DEBUG("merit: %f", merit);
		LOG_DEBUG("model_merit: %f", model_merit);
		LOG_DEBUG("new_merit: %f", new_merit);
		LOG_DEBUG("constant cost term: %f", constant_cost);

		exact_merit_improve = merit - new_merit;
		approx_merit_improve = merit - model_merit;
		merit_improve_ratio = exact_merit_improve / approx_merit_improve;

		LOG_DEBUG("exact_merit_improve: %f", exact_merit_improve);
		LOG_DEBUG("approx_merit_improve: %f", approx_merit_improve);
		LOG_DEBUG("merit_improve_ratio: %f", merit_improve_ratio);

		if (approx_merit_improve < -1) {
			LOG_ERROR("Approximate merit function got worse: %f", approx_merit_improve);
			LOG_DEBUG("Shrinking trust region");

			Xeps *= cfg::trust_shrink_ratio;
			Ueps *= cfg::trust_shrink_ratio;
			solution_accepted = false;
//			return INFINITY;
		} else if (approx_merit_improve < cfg::min_approx_improve) {
			LOG_DEBUG("Converged: improvement small enough");
			J = Jopt; U = Uopt;
			solution_accepted = true;
			return sys.entropy(J, U, P, obstacles);
		} else if ((exact_merit_improve < 0) || (merit_improve_ratio < cfg::improve_ratio_threshold)) {
			Xeps *= cfg::trust_shrink_ratio;
			Ueps *= cfg::trust_shrink_ratio;
			LOG_DEBUG("Shrinking trust region size to: %2.6f %2.6f", Xeps, Ueps);
			solution_accepted = false;
		} else {
			// expand Xeps and Ueps
			Xeps *= cfg::trust_expand_ratio;
			Ueps *= cfg::trust_expand_ratio;
			LOG_DEBUG("Accepted, Increasing trust region size to:  %2.6f %2.6f", Xeps, Ueps);

			meritopt = sys.entropy(Jopt, Uopt, P, obstacles);
			gradopt = sys.entropy_grad(Jopt, Uopt, P, obstacles);

			L_BFGS(J, U, grad, Jopt, Uopt, gradopt, hess);

			J = Jopt; U = Uopt;
			solution_accepted = true;
		}

	}

	return sys.entropy(J, U, P, obstacles);
}


inline double uniform(double low, double high) {
	return (high - low)*(rand() / double(RAND_MAX)) + low;
}

void init_obstacles_and_particles(pr2_sim::Camera& cam,
		std::vector<geometry3d::Triangle>& obstacles, MatrixP& P) {
	Matrix4d cam_pose = cam.get_pose();
	Matrix3d cam_rot = cam_pose.block<3,3>(0,0);
	Vector3d cam_pos = cam_pose.block<3,1>(0,3);

//	std::vector<geometry3d::Triangle> obstacles_cam;
//
//	obstacles_cam.push_back(geometry3d::Triangle({0,0,.75}, {0,.1,.75}, {.05,.1,.75}));
//	obstacles_cam.push_back(geometry3d::Triangle({0,0,.75}, {.05,0,.75}, {.05,.1,.75}));
//
//	obstacles_cam.push_back(geometry3d::Triangle({-.2,0,.75}, {-.2,.1,.75}, {-.25,.1,.75}));
//	obstacles_cam.push_back(geometry3d::Triangle({-.2,0,.75}, {-.25,0,.75}, {-.25,.1,.75}));

	obstacles.clear();
//	for(const geometry3d::Triangle& obstacle_cam : obstacles_cam) {
//		obstacles.push_back(geometry3d::Triangle(cam_rot*obstacle_cam.a+cam_pos,
//				cam_rot*obstacle_cam.b+cam_pos,
//				cam_rot*obstacle_cam.c+cam_pos));
//	}

	obstacles.push_back(geometry3d::Triangle({.5,-.5,.55}, {.5,-.55,.55}, {.5,-.55,.75}));
	obstacles.push_back(geometry3d::Triangle({.5,-.5,.55}, {.5,-.5,.75}, {.5,-.55,.75}));

	// table
//	obstacles.push_back(geometry3d::Triangle({0,-1,.55}, {0,1,.55}, {1,1,.55}));
//	obstacles.push_back(geometry3d::Triangle({0,-1,.55}, {1,-1,.55}, {1,1,.55}));

	P.setZero();
	std::vector<geometry3d::TruncatedPyramid> truncated_frustum = cam.truncated_view_frustum(obstacles, true);
	int num_particles = 0;
	int particles_per_obstacle = double(M_DIM)/(obstacles.size());
	for(int i=0; i < obstacles.size(); i+=2) {
		geometry3d::Triangle& obstacle = obstacles[i];
		double x_min = std::min(obstacle.a(0), std::min(obstacle.b(0), obstacle.c(0)));
		double x_max = std::max(obstacle.a(0), std::max(obstacle.b(0), obstacle.c(0))) + .3;
		double y_min = std::min(obstacle.a(1), std::min(obstacle.b(1), obstacle.c(1))) - .1;
		double y_max = std::max(obstacle.a(1), std::max(obstacle.b(1), obstacle.c(1))) + .1;
		double z_min = std::min(obstacle.a(2), std::min(obstacle.b(2), obstacle.c(2)));
		double z_max = std::max(obstacle.a(2), std::max(obstacle.b(2), obstacle.c(2))) + .1;

		while(num_particles < (i+2)*particles_per_obstacle) {
			Vector3d p = Vector3d(uniform(x_min, x_max), uniform(y_min, y_max), uniform(z_min, z_max));
			if (!cam.is_in_fov(p, truncated_frustum)) {
				P.col(num_particles++) = p;
			}
		}
	}

}

void init_trajectory(StdVectorJ& J, StdVectorU& U, pr2_sim::Arm& arm, PR2EihSystem& sys) {
	Matrix4d start_pose = arm.fk(J[0]);
	Matrix4d next_pose = start_pose;
	VectorJ next_joints;
	for(int t=0; t < T-1; ++t) {
		U[t].setZero();
//		next_pose.block<3,1>(0,3) += Vector3d(0, 0, .005);
//		if (arm.ik(next_pose, next_joints)) {
//			U[t] = next_joints - J[t];
//		} else {
//			U[t] = VectorU::Zero();
//		}
	}

	for(int t=0; t < T-1; ++t) {
		J[t+1] = sys.dynfunc(J[t], U[t], VectorQ::Zero());
	}
}

int main(int argc, char* argv[]) {
	srand(time(0));

	// setup system
	pr2_sim::Simulator sim(true, false);
	pr2_sim::Arm arm(pr2_sim::Arm::right, &sim);
	pr2_sim::Camera cam(&arm, &sim);
	PR2EihSystem sys(&sim, &arm, &cam);

	pr2_sim::Arm other_arm(pr2_sim::Arm::left, &sim);
	other_arm.set_posture(pr2_sim::Arm::Posture::mantis);
	arm.set_posture(pr2_sim::Arm::Posture::mantis);

	// setup environment
	std::vector<geometry3d::Triangle> obstacles;
	MatrixP P_t, P_tp1;
	init_obstacles_and_particles(cam, obstacles, P_t);

	// setup initial state
	VectorJ j_t, j_tp1;
	j_t = arm.get_joints();

	// initialize state and controls
	StdVectorU U(T-1, VectorU::Zero());
	StdVectorJ J(T, j_t);

	// initialize FORCES
	pr2eihMPC_params problem;
	pr2eihMPC_output output;
	pr2eihMPC_info info;

	setup_mpc_vars(problem, output);
	util::Timer forces_timer, entropy_timer;

	while(true) {
		init_trajectory(J, U, arm, sys);

		LOG_INFO("Current state");
		sys.plot(J, P_t, obstacles);

		util::Timer_tic(&entropy_timer);
		double initial_entropy = sys.entropy(J, U, P_t, obstacles);
		double entropy_time = util::Timer_toc(&entropy_timer);
		LOG_INFO("Initial cost: %4.5f", initial_entropy);
		LOG_INFO("Entropy time: %5.3f ms", 1e3*entropy_time);

		// optimize
		util::Timer_tic(&forces_timer);
		double entropy = pr2_eih_collocation(J, U, P_t, obstacles, sys, problem, output, info);
		double forces_time = util::Timer_toc(&forces_timer);

		LOG_INFO("Optimized cost: %4.5f", entropy);
		LOG_INFO("Solve time: %5.3f ms", 1e3*forces_time);

		for(int t=0; t < T-1; ++t) {
			J[t+1] = sys.dynfunc(J[t], U[t], VectorQ::Zero());
		}

		LOG_INFO("Displaying optimized trajectory");
		sys.plot(J, P_t, obstacles);

		j_tp1 = sys.dynfunc(j_t, U[0], VectorQ::Zero());
		VectorP W_tp1 = sys.update_particle_weights(j_tp1, P_t, (1/double(M_DIM))*VectorP::Ones(), obstacles);
		P_tp1 = sys.low_variance_sampler(P_t, W_tp1);

		j_t = j_tp1;
		P_t = P_tp1;

		// TODO: replace with optimization vars from [1, ..., T-1]
		J = StdVectorJ(T, j_t);
		U = StdVectorU(T-1, VectorU::Zero());

	}

	return 0;
}
