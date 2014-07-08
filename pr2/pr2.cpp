#include "system/pr2-system.h"

#include <stdio.h>

#include <boost/preprocessor/iteration/local.hpp>

extern "C" {
#include "pr2MPC.h"
pr2MPC_FLOAT **H, **f, **lb, **ub, **z, *c;
}

const int T = TIMESTEPS;
const double INFTY = 1e10;

namespace cfg {
const double alpha_init = .01; // 1
const double alpha_gain = 3; // 3
const double alpha_epsilon = .1; // .001
const double alpha_max_increases = 5; // 10

const double improve_ratio_threshold = 1e-1; // .1
const double min_approx_improve = 1e-1; // 1
const double min_trust_box_size = .1; // .1
const double trust_shrink_ratio = .5; // .5
const double trust_expand_ratio = 1.5; // 1.5
}

void setup_mpc_vars(pr2MPC_params& problem, pr2MPC_output& output) {
	// inputs
	H = new pr2MPC_FLOAT*[T];
	f = new pr2MPC_FLOAT*[T];
	lb = new pr2MPC_FLOAT*[T];
	ub = new pr2MPC_FLOAT*[T];
	c = problem.c1;

	// output
	z = new pr2MPC_FLOAT*[T];

#define SET_VARS(n)    \
		H[ BOOST_PP_SUB(n,1) ] = problem.H##n ;  \
		f[ BOOST_PP_SUB(n,1) ] = problem.f##n ;  \
		lb[ BOOST_PP_SUB(n,1) ] = problem.lb##n ;	\
		ub[ BOOST_PP_SUB(n,1) ] = problem.ub##n ;	\
		z[ BOOST_PP_SUB(n,1) ] = output.z##n ;

#define BOOST_PP_LOCAL_MACRO(n) SET_VARS(n)
#define BOOST_PP_LOCAL_LIMITS (1, TIMESTEPS)
#include BOOST_PP_LOCAL_ITERATE()

	for(int i=0; i < J_DIM; ++i) { c[i] = INFTY; }

	for(int t=0; t < T-1; ++t) {
		for(int i=0; i < (J_DIM+U_DIM); ++i) { H[t][i] = INFTY; }
		for(int i=0; i < (J_DIM+U_DIM); ++i) { f[t][i] = INFTY; }
		for(int i=0; i < (J_DIM+U_DIM); ++i) { lb[t][i] = INFTY; }
		for(int i=0; i < (J_DIM+U_DIM); ++i) { ub[t][i] = INFTY; }
		for(int i=0; i < (J_DIM+U_DIM); ++i) { z[t][i] = INFTY; }
	}
	for(int i=0; i < J_DIM; ++i) { H[T-1][i] = INFTY; }
	for(int i=0; i < J_DIM; ++i) { f[T-1][i] = INFTY; }
	for(int i=0; i < J_DIM; ++i) { lb[T-1][i] = INFTY; }
	for(int i=0; i < J_DIM; ++i) { ub[T-1][i] = INFTY; }
	for(int i=0; i < J_DIM; ++i) { z[T-1][i] = INFTY; }

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
//
//	std::cout << "\n\n";

	for(int i=0; i < (J_DIM); ++i) { if (c[i] > INFTY/2) { LOG_ERROR("isValid0"); return false; } }
	for(int i=0; i < (J_DIM); ++i) { if (c[i] < lb[0][i]) { LOG_ERROR("isValid1"); return false; } }
	for(int i=0; i < (J_DIM); ++i) { if (c[i] > ub[0][i]) { LOG_ERROR("isValid2"); return false; } }

	for(int t = 0; t < T-1; ++t) {
		for(int i=0; i < (J_DIM+U_DIM); ++i) { if (H[t][i] > INFTY/2) { LOG_ERROR("isValid5"); return false; } }
		for(int i=0; i < (J_DIM+U_DIM); ++i) { if (f[t][i] > INFTY/2) { LOG_ERROR("isValid6"); return false; } }
		for(int i=0; i < (J_DIM+U_DIM); ++i) { if (lb[t][i] > INFTY/2) { LOG_ERROR("isValid7"); return false; } }
		for(int i=0; i < (J_DIM+U_DIM); ++i) {if (ub[t][i] > INFTY/2) { LOG_ERROR("isValid8"); return false; } }
		for(int i=0; i < (J_DIM+U_DIM); ++i) {if (ub[t][i] < lb[t][i]) { LOG_ERROR("isValid9"); return false; } }
	}
	for(int i=0; i < (J_DIM); ++i) { if (H[T-1][i] > INFTY/2) { LOG_ERROR("isValid10"); return false; } }
	for(int i=0; i < (J_DIM); ++i) { if (f[T-1][i] > INFTY/2) { LOG_ERROR("isValid11"); return false; } }
	for(int i=0; i < (J_DIM); ++i) { if (lb[T-1][i] > INFTY/2) { LOG_ERROR("isValid12"); return false; } }
	for(int i=0; i < (J_DIM); ++i) { if (ub[T-1][i] > INFTY/2) { LOG_ERROR("isValid13"); return false; } }
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

double pr2_approximate_collocation(StdVectorJ& J, StdVectorU& U, const MatrixJ& j_sigma0,
		const std::vector<ParticleGaussian>& particle_gmm, const double alpha,
		PR2System& sys, pr2MPC_params &problem, pr2MPC_output &output, pr2MPC_info &info) {
	int max_iter = 100;
	double Xeps = .1; // .5
	double Ueps = .1; // .5

	double merit = 0, meritopt = 0;
	double constant_cost, hessian_constant, jac_constant;
	VectorTOTAL grad = VectorTOTAL::Zero();
	MatrixTOTAL hess = MatrixTOTAL::Identity();

	StdVectorJ Jopt(T, VectorJ::Zero());
	StdVectorU Uopt(T-1, VectorU::Zero());
	VectorTOTAL gradopt = VectorTOTAL::Zero();
	double optcost, model_merit, new_merit;
	double approx_merit_improve, exact_merit_improve, merit_improve_ratio;

	LOG_DEBUG("Initial trajectory cost: %4.10f", sys.cost_gmm(J, j_sigma0, U, particle_gmm, alpha));

	int index = 0;
	bool solution_accepted = true;
	for(int it=0; it < max_iter; ++it) {

		LOG_DEBUG(" ");
		LOG_DEBUG(" ");
		LOG_DEBUG("Iter: %d", it);

		// only compute gradient/hessian if P/U has been changed
		if (solution_accepted) {

			if (it == 0) {
				merit = sys.cost_gmm(J, j_sigma0, U, particle_gmm, alpha);
				grad = sys.cost_gmm_grad_ripped(J, j_sigma0, U, particle_gmm, alpha);
			} else {
				merit = meritopt; // since L-BFGS calculation required it
				grad = gradopt;
			}

			std::cout << "gradient:\n" << grad << "\n";

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
		int exitflag = pr2MPC_solve(&problem, &output, &info);
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

//		LOG_DEBUG("Displaying Jopt");
//		sys.display(Jopt, particle_gmm);

		model_merit = optcost + constant_cost; // need to add constant terms that were dropped

		new_merit = sys.cost_gmm(Jopt, j_sigma0, Uopt, particle_gmm, alpha);

		LOG_DEBUG("merit: %f", merit);
		LOG_DEBUG("model_merit: %f", model_merit);
		LOG_DEBUG("new_merit: %f", new_merit);
		LOG_DEBUG("constant cost term: %f", constant_cost);

		approx_merit_improve = merit - model_merit;
		exact_merit_improve = merit - new_merit;
		merit_improve_ratio = exact_merit_improve / approx_merit_improve;

		LOG_DEBUG("approx_merit_improve: %f", approx_merit_improve);
		LOG_DEBUG("exact_merit_improve: %f", exact_merit_improve);
		LOG_DEBUG("merit_improve_ratio: %f", merit_improve_ratio);

		if (approx_merit_improve < -1) {
			LOG_ERROR("Approximate merit function got worse: %f", approx_merit_improve);
			LOG_DEBUG("Shrinking trust region");

			Xeps *= cfg::trust_shrink_ratio;
			Ueps *= cfg::trust_shrink_ratio;
			solution_accepted = false;
//			return INFTY;
		} else if (approx_merit_improve < cfg::min_approx_improve) {
			LOG_DEBUG("Converged: improvement small enough");
			J = Jopt; U = Uopt;
			solution_accepted = true;
			return sys.cost_gmm(J, j_sigma0, U, particle_gmm, alpha);
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

			meritopt = sys.cost_gmm(Jopt, j_sigma0, Uopt, particle_gmm, alpha);
			gradopt = sys.cost_gmm_grad_ripped(Jopt, j_sigma0, Uopt, particle_gmm, alpha);
			L_BFGS(J, U, grad, Jopt, Uopt, gradopt, hess);

			J = Jopt; U = Uopt;
			solution_accepted = true;
		}

	}

	return sys.cost_gmm(J, j_sigma0, U, particle_gmm, alpha);
}

double pr2_collocation(StdVectorJ& J, StdVectorU& U, const MatrixJ& j_sigma0,
		const std::vector<ParticleGaussian>& particle_gmm,
		PR2System& sys, pr2MPC_params &problem, pr2MPC_output &output, pr2MPC_info &info) {
	double alpha = cfg::alpha_init;
	double cost = INFINITY;

	for(int num_alpha_increases=0; num_alpha_increases < cfg::alpha_max_increases; ++num_alpha_increases) {
		LOG_DEBUG("Calling approximate collocation with alpha = %4.2f", alpha);
		cost = pr2_approximate_collocation(J, U, j_sigma0, particle_gmm, alpha, sys, problem, output, info);

		LOG_DEBUG("Reintegrating trajectory");
		for(int t=0; t < T-1; ++t) {
			J[t+1] = sys.dynfunc(J[t], U[t], VectorQ::Zero());
		}

//		double max_delta_diff = -INFINITY;
//		for(int t=0; t < T; ++t) {
//			double delta_alpha = sys.delta_matrix(J[t], particle_gmm[0].mean, alpha, particle_gmm[0].ODF)(J_DIM,J_DIM);
//			double delta_inf = sys.delta_matrix(J[t], particle_gmm[0].mean, INFTY, particle_gmm[0].ODF)(J_DIM,J_DIM);
//			max_delta_diff = std::max(max_delta_diff, fabs(delta_alpha - delta_inf));
//		}
//
//		LOG_DEBUG(" ");
//		LOG_DEBUG("Max delta difference: %4.2f", max_delta_diff);
//		if (max_delta_diff < cfg::alpha_epsilon) {
//			LOG_DEBUG("Max delta difference < %4.10f, exiting minimize merit", cfg::alpha_epsilon);
//			break;
//		}

		LOG_DEBUG("Increasing alpha by gain %4.5f", cfg::alpha_gain);
		alpha *= cfg::alpha_gain;
	}

	return cost;
}

void init_collocation(const VectorJ& j0, const MatrixP& P, PR2System& sys,
		StdVectorJ& J, StdVectorU& U, std::vector<ParticleGaussian>& particle_gmm) {
	// re-initialize GMM from PF
	sys.fit_gaussians_to_pf(P, particle_gmm);

	// only take max gaussian
	particle_gmm = std::vector<ParticleGaussian>(1, particle_gmm[0]);
//	particle_gmm[0].cov = .01*Matrix3d::Identity();
//	particle_gmm[0].cov *= 5000;

	for(int i=0; i < particle_gmm.size(); ++i) {
		particle_gmm[i].cov *= 1000; // 5000
		std::cout << particle_gmm[i].pct << "\n";
		std::cout << particle_gmm[i].cov << "\n\n";

		particle_gmm[i].ODF = sys.get_ODF(particle_gmm[i].mean);
	}

	// TODO: try a non-zero initialization
	VectorU uinit = VectorU::Zero();

	// integrate trajectory
	J[0] = j0;
	for(int t=0; t < T-1; ++t) {
		U[t] = uinit;
		J[t+1] = sys.dynfunc(J[t], U[t], VectorQ::Zero());
	}
}

MatrixP init_particles(rave::EnvironmentBasePtr env) {
	rave::KinBodyPtr table = env->GetKinBody("table");
	rave::KinBody::LinkPtr base = table->GetLink("base");
	rave::Vector extents = base->GetGeometry(0)->GetBoxExtents();

	rave::Vector table_pos = table->GetTransform().trans;
	double x_min, x_max, y_min, y_max, z_min, z_max;
//	x_min = table_pos.x - extents.x;
//	x_max = table_pos.x + extents.x;
//	y_min = table_pos.y - extents.y;
//	y_max = table_pos.y + extents.y;
//	z_min = table_pos.z + extents.z;
//	z_max = table_pos.z + extents.z + .2;

	x_min = table_pos.x - extents.x;
	x_max = table_pos.x + extents.x;
	y_min = table_pos.y - extents.y;
	y_max = table_pos.y;// + extents.y;
	z_min = table_pos.z + extents.z - .1;
	z_max = table_pos.z + extents.z - .2;// + .2;

	MatrixP P;

	// uniform
	for(int m=0; m < M_DIM; ++m) {
		P(0,m) = pr2_utils::uniform(x_min, x_max);
		P(1,m) = pr2_utils::uniform(y_min, y_max);
		P(2,m) = pr2_utils::uniform(z_min, z_max);
	}

	// two clumps
//	for(int m=0; m < M_DIM; ++m) {
//		if (m < M_DIM/2) {
//			P(0,m) = mm_utils::uniform(x_min, .7*x_min + .3*x_max);
//			P(1,m) = mm_utils::uniform(y_min, y_max);
//			P(2,m) = mm_utils::uniform(z_min, z_max);
//		} else {
//			P(0,m) = mm_utils::uniform(.3*x_min + .7*x_max, x_max);
//			P(1,m) = mm_utils::uniform(y_min, y_max);
//			P(2,m) = mm_utils::uniform(z_min, z_max);
//		}
//	}

	return P;
}

int main(int argc, char* argv[]) {
//	Vector3d object(3.35, -1.11, 0.8);
	Vector3d object = Vector3d(3.5+.1, -1.2-.1, .74-.1);
	Arm::ArmType arm_type = Arm::ArmType::right;
	bool view = true;
	PR2System sys(object, arm_type, view);

	PR2* brett = sys.get_brett();
	rave::EnvironmentBasePtr env = brett->get_env();
	Camera* cam = sys.get_camera();
	Arm* arm = sys.get_arm();
	arm->set_posture(Arm::Posture::mantis);

	// initialize starting state, belief, and pf
	VectorJ j_t, j_t_real, j_tp1, j_tp1_real;
	j_t = arm->get_joint_values();
	j_t_real = j_t; // TODO: have them be different

	MatrixJ j_sigma0 = .01*MatrixJ::Identity(); // TODO: never actually update it in MPC

	MatrixP P_t, P_tp1;
	P_t = init_particles(env);

	std::vector<ParticleGaussian> particle_gmm;

	// initialize state and controls
	StdVectorU U(T-1, VectorU::Zero());
	StdVectorJ J(T);

	// initialize FORCES
	pr2MPC_params problem;
	pr2MPC_output output;
	pr2MPC_info info;

	setup_mpc_vars(problem, output);
	util::Timer forces_timer;

	bool stop_condition = false;
	for(int iter=0; !stop_condition; iter++) {
		arm->set_joint_values(j_t);

		LOG_INFO("Updating internal TSDF and kinfu\n");
		StdVector3d pc = cam->get_pc(j_t_real);
		Matrix<double,HEIGHT_FULL,WIDTH_FULL> full_zbuffer = cam->get_full_zbuffer(j_t_real, pc);
		Matrix4d cam_pose = cam->get_pose(j_t_real);
		sys.update(pc, full_zbuffer, cam_pose);

		LOG_INFO("MPC iteration: %d",iter);
		init_collocation(j_t, P_t, sys, J, U, particle_gmm);

		LOG_INFO("Current state");
		sys.display(j_t, particle_gmm);

		LOG_INFO("Initialized trajectory");
		sys.display(J, particle_gmm);

		// optimize
		util::Timer_tic(&forces_timer);
		double cost = pr2_collocation(J, U, j_sigma0, particle_gmm, sys, problem, output, info);
		double forces_time = util::Timer_toc(&forces_timer);

		LOG_INFO("Optimized cost: %4.5f", cost);
		LOG_INFO("Solve time: %5.3f ms", forces_time*1000);

		for(int t=0; t < T-1; ++t) {
			J[t+1] = sys.dynfunc(J[t], U[t], VectorQ::Zero());
		}

		LOG_INFO("Post-optimization");
		sys.display(J, particle_gmm);

		sys.execute_control_step(j_t_real, j_t, U[0], P_t, particle_gmm[0].ODF, j_tp1_real, j_tp1, P_tp1);

		// TODO: stop condition

		j_t_real = j_tp1_real;
		j_t = j_tp1;
		P_t = P_tp1;

	}

	return 0;
}
