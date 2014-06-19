#include "include/planar-system.h"
#include "include/gmm.h"
#include "../util/Timer.h"

extern "C" {
#include "planarMPC.h"
planarMPC_FLOAT **H, **f, **lb, **ub, **z, *c, *A, *b;
}

const int T = TIMESTEPS;
const double INFTY = 1e10;

//const bool USE_FADBAD = true;

namespace cfg {
const double alpha_init = .01; // 1
const double alpha_gain = 3; // 3
const double alpha_epsilon = .1; // .001

const double improve_ratio_threshold = 1e-1; // .1
const double min_approx_improve = 1e-1; // 1
const double min_trust_box_size = .1; // .1
const double trust_shrink_ratio = .5; // .5
const double trust_expand_ratio = 1.5; // 1.5
}

void setup_mpc_vars(planarMPC_params& problem, planarMPC_output& output) {
	// inputs
	H = new planarMPC_FLOAT*[T];
	f = new planarMPC_FLOAT*[T];
	lb = new planarMPC_FLOAT*[T];
	ub = new planarMPC_FLOAT*[T];

	c = problem.c1;
	A = problem.A10;
	b = problem.b10;

	// output
	z = new planarMPC_FLOAT*[T];

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
	for(int i=0; i < (2*C_DIM)*J_DIM; ++i) { A[i] = INFTY; }
	for(int i=0; i < (2*C_DIM); ++i) { b[i] = INFTY; }

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

template<typename Derived>
inline void fill_col_major(double *X, const MatrixBase<Derived>& M) {
	int index = 0;
	for(int j=0; j < M.cols(); ++j) {
		for(int i=0; i < M.rows(); ++i) {
			X[index++] = M(i,j);
		}
	}
}

bool is_valid_inputs()
{
	for(int t = 0; t < T-1; ++t) {
		std::cout << "\n\nt: " << t << "\n";

		if (t == 0) {
			std::cout << "\nc[0]:\n";
			for(int i=0; i < (J_DIM); ++i) {
				std::cout << c[i] << " ";
			}
		}

		std::cout << "\nH[" << t << "]: ";
		for(int i=0; i < (J_DIM+U_DIM); ++i) {
			std::cout << H[t][i] << " ";
		}

		std::cout << "\nf[" << t << "]: ";
		for(int i=0; i < (J_DIM+U_DIM); ++i) {
			std::cout << f[t][i] << " ";
		}

		std::cout << "\nlb[" << t << "]: ";
		for(int i=0; i < (J_DIM+U_DIM); ++i) {
			std::cout << lb[t][i] << " ";
		}

		std::cout << "\nub[" << t << "]: ";
		for(int i=0; i < (J_DIM+U_DIM); ++i) {
			std::cout << ub[t][i] << " ";
		}
	}
	std::cout << "\n\nt: " << T-1 << "\n";

	std::cout << "\nH[" << T-1 << "]: ";
	for(int i=0; i < (J_DIM); ++i) {
		std::cout << H[T-1][i] << " ";
	}

	std::cout << "\nf[" << T-1 << "]: ";
	for(int i=0; i < (J_DIM); ++i) {
		std::cout << f[T-1][i] << " ";
	}

	std::cout << "\nlb[" << T-1 << "]: ";
	for(int i=0; i < (J_DIM); ++i) {
		std::cout << lb[T-1][i] << " ";
	}

	std::cout << "\nub[" << T-1 << "]: ";
	for(int i=0; i < (J_DIM); ++i) {
		std::cout << ub[T-1][i] << " ";
	}

	std::cout << "\nA[" << T-1 << "]:\n";
	for(int i=0; i < 2*C_DIM; ++i) {
		for(int j=0; j < J_DIM; ++j) {
			std::cout << A[i+j*(2*C_DIM)] << " ";
		}
		std::cout << "\n";
	}
	std::cout << "\n";

	std::cout << "b[" << T-1 << "]: ";
	for(int i=0; i < 2*C_DIM; ++i) {
		std::cout << b[i] << " ";
	}

	std::cout << "\n\n";

	for(int i=0; i < (J_DIM); ++i) { if (c[i] > INFTY/2) { return false; } }
	for(int i=0; i < (2*C_DIM)*J_DIM; ++i) { if (A[i] > INFTY/2) { return false; } }
	for(int i=0; i < (2*C_DIM); ++i) { if (b[i] > INFTY/2) { return false; } }

	for(int t = 0; t < T-1; ++t) {
		for(int i=0; i < (J_DIM+U_DIM); ++i) { if (H[t][i] > INFTY/2) { return false; } }
		for(int i=0; i < (J_DIM+U_DIM); ++i) { if (f[t][i] > INFTY/2) { return false; } }
		for(int i=0; i < (J_DIM+U_DIM); ++i) { if (lb[t][i] > INFTY/2) { return false; } }
		for(int i=0; i < (J_DIM+U_DIM); ++i) {if (ub[t][i] > INFTY/2) { return false; } }
		for(int i=0; i < (J_DIM+U_DIM); ++i) {if (ub[t][i] < lb[t][i]) { return false; } }
	}
	for(int i=0; i < (J_DIM); ++i) { if (H[T-1][i] > INFTY/2) { return false; } }
	for(int i=0; i < (J_DIM); ++i) { if (f[T-1][i] > INFTY/2) { return false; } }
	for(int i=0; i < (J_DIM); ++i) { if (lb[T-1][i] > INFTY/2) { return false; } }
	for(int i=0; i < (J_DIM); ++i) { if (ub[T-1][i] > INFTY/2) { return false; } }
	for(int i=0; i < (J_DIM); ++i) {if (ub[T-1][i] < lb[T-1][i]) { return false; } }

	return true;
}

void L_BFGS(const std::vector<vec<J_DIM>, aligned_allocator<vec<J_DIM>>>& J,
		const std::vector<vec<U_DIM>, aligned_allocator<vec<U_DIM>>>& U, const vec<TOTAL_VARS> &grad,
		const std::vector<vec<J_DIM>, aligned_allocator<vec<J_DIM>>> &Jopt,
		const std::vector<vec<U_DIM>, aligned_allocator<vec<U_DIM>>> &Uopt, const vec<TOTAL_VARS> &gradopt,
		mat<TOTAL_VARS, TOTAL_VARS> &hess) {
	vec<TOTAL_VARS> s = vec<TOTAL_VARS>::Zero();

	int index = 0;
	for(int t=0; t < T-1; ++t) {
		s.segment<J_DIM>(index) = Jopt[t] - J[t];
		index += J_DIM;
		s.segment<U_DIM>(index) = Uopt[t] - U[t];
		index += U_DIM;
	}
	s.segment<J_DIM>(index) = Jopt[T-1] - J[T-1];

	vec<TOTAL_VARS> y = gradopt - grad;

	double theta;
	vec<TOTAL_VARS> hess_s = hess*s;

	bool decision = s.dot(y) >= .2*s.dot(hess_s);
	if (decision) {
		theta = 1;
	} else {
		theta = (.8*s.dot(hess_s))/(s.dot(hess_s) - s.dot(y));
	}

	vec<TOTAL_VARS> r = theta*y + (1-theta)*hess_s;

	hess = hess - (hess_s*hess_s.transpose())/(s.dot(hess_s)) + (r*r.transpose())/s.dot(r);
}

double planar_collocation(std::vector<vec<J_DIM>, aligned_allocator<vec<J_DIM>>>& J,
		std::vector<vec<U_DIM>, aligned_allocator<vec<U_DIM>>>& U,
		const mat<J_DIM,J_DIM>& j_sigma0,
		const std::vector<PlanarGaussian>& planar_gmm,
		const double alpha,
		PlanarSystem& sys, planarMPC_params &problem, planarMPC_output &output, planarMPC_info &info) {
	int max_iter = 100;
	double Xeps = .5;
	double Ueps = .5;

	double merit = 0, meritopt = 0;
	double constant_cost, hessian_constant, jac_constant;
	vec<TOTAL_VARS> grad = vec<TOTAL_VARS>::Zero();
	mat<TOTAL_VARS,TOTAL_VARS> hess = mat<TOTAL_VARS,TOTAL_VARS>::Identity();

	std::vector<vec<J_DIM>, aligned_allocator<vec<J_DIM>>> Jopt(T, vec<J_DIM>::Zero());
	std::vector<vec<U_DIM>, aligned_allocator<vec<U_DIM>>> Uopt(T-1, vec<U_DIM>::Zero());
	vec<TOTAL_VARS> gradopt = vec<TOTAL_VARS>::Zero();
	double optcost, model_merit, new_merit;
	double approx_merit_improve, exact_merit_improve, merit_improve_ratio;

	LOG_DEBUG("Initial trajectory cost: %4.10f", sys.cost_gmm(J, j_sigma0, U, planar_gmm, alpha));

	int index = 0;
	bool solution_accepted = true;
	for(int it=0; it < max_iter; ++it) {

		LOG_DEBUG(" ");
		LOG_DEBUG(" ");
		LOG_DEBUG("Iter: %d", it);

		// only compute gradient/hessian if P/U has been changed
		if (solution_accepted) {

			if (it == 0) {
				merit = sys.cost_gmm(J, j_sigma0, U, planar_gmm, alpha);
				grad = sys.cost_gmm_grad(J, j_sigma0, U, planar_gmm, alpha);

//				sys.cost_and_cost_grad(X, sigma0, U, alpha, USE_FADBAD, merit, grad);
			} else {
				merit = meritopt; // since L-BFGS calculation required it
				grad = gradopt;
			}

			vec<TOTAL_VARS> diaghess = hess.diagonal();

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
				vec<(J_DIM+U_DIM)> zbar;
				zbar.segment<J_DIM>(0) = J[t];
				zbar.segment<U_DIM>(J_DIM) = U[t];

				for(int i=0; i < (J_DIM+U_DIM); ++i) {
					hessian_constant += H[t][i]*zbar(i)*zbar(i);
					jac_constant -= grad(index)*zbar(i);
					f[t][i] = grad(index) - H[t][i]*zbar(i);
					index++;
				}
			}

			vec<J_DIM> zbar = J[T-1];

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


		vec<X_DIM> x_min, x_max;
		vec<U_DIM> u_min, u_max;
		sys.get_limits(x_min, x_max, u_min, u_max);

		// set trust region bounds based on current trust region size
		for(int t=0; t < T; ++t) {
			index = 0;
			for(int i=0; i < J_DIM; ++i) {
				lb[t][index] = std::max(x_min(i), J[t](i) - Xeps);
				ub[t][index] = std::min(x_max(i), J[t](i) + Xeps);
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

		// end goal constraint on end effector
		vec<C_DIM> goal_pos = planar_gmm[0].obj_mean;
		double goal_delta = .01;

		vec<E_DIM> ee_jT = J[T-1].segment<E_DIM>(0);
		vec<C_DIM> ee_pos = sys.get_ee_pos(ee_jT);
		mat<C_DIM,E_DIM> ee_pos_jac;
		sys.get_ee_pos_jac(ee_jT, ee_pos_jac);
		mat<C_DIM,J_DIM> ee_pos_jac_full = mat<C_DIM,J_DIM>::Zero();
		ee_pos_jac_full.block<C_DIM,E_DIM>(0,0) = ee_pos_jac;

		mat<2*C_DIM,J_DIM> Amat;
		Amat.block<C_DIM,J_DIM>(0,0) = ee_pos_jac_full;
		Amat.block<C_DIM,J_DIM>(C_DIM,0) = -ee_pos_jac_full;
		fill_col_major(A, Amat);

		vec<2*C_DIM> bVec;
		bVec.segment<C_DIM>(0) = goal_pos - ee_pos + ee_pos_jac_full*J[T-1] + goal_delta*vec<C_DIM>::Ones();
		bVec.segment<C_DIM>(C_DIM) = -goal_pos + ee_pos - ee_pos_jac_full*J[T-1] + goal_delta*vec<C_DIM>::Ones();
		fill_col_major(b, bVec);

		// Verify problem inputs
//		if (!is_valid_inputs()) {
//			LOG_ERROR("Inputs are not valid!");
//			exit(0);
//		}

		// call FORCES
		int exitflag = planarMPC_solve(&problem, &output, &info);
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
			exit(-1);
		}

		model_merit = optcost + constant_cost; // need to add constant terms that were dropped

		new_merit = sys.cost_gmm(Jopt, j_sigma0, Uopt, planar_gmm, alpha);

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

		if (approx_merit_improve < -1e-5) {
			LOG_ERROR("Approximate merit function got worse: %f", approx_merit_improve);
			LOG_ERROR("Failure!");
			return INFTY;
		} else if (approx_merit_improve < cfg::min_approx_improve) {
			LOG_DEBUG("Converged: improvement small enough");
			J = Jopt; U = Uopt;
			solution_accepted = true;
			return sys.cost_gmm(J, j_sigma0, U, planar_gmm, alpha);
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

			meritopt = sys.cost_gmm(Jopt, j_sigma0, Uopt, planar_gmm, alpha);
			gradopt = sys.cost_gmm_grad(Jopt, j_sigma0, Uopt, planar_gmm, alpha);
//			sys.cost_and_cost_grad(Xopt, sigma0, Uopt, alpha, USE_FADBAD, meritopt, gradopt);
			L_BFGS(J, U, grad, Jopt, Uopt, gradopt, hess);

			J = Jopt; U = Uopt;
			solution_accepted = true;
		}

	}

	return sys.cost_gmm(J, j_sigma0, U, planar_gmm, alpha);
}

double planar_minimize_merit(std::vector<vec<J_DIM>, aligned_allocator<vec<J_DIM>>>& J,
		std::vector<vec<U_DIM>, aligned_allocator<vec<U_DIM>>>& U,
		const mat<J_DIM,J_DIM>& j_sigma0,
		const std::vector<PlanarGaussian>& planar_gmm,
		PlanarSystem& sys, planarMPC_params &problem, planarMPC_output &output, planarMPC_info &info) {
	double alpha = cfg::alpha_init;
	double cost = INFINITY;

	while(true) {
		LOG_DEBUG("Calling collocation with alpha = %4.2f", alpha);
		cost = planar_collocation(J, U, j_sigma0, planar_gmm, alpha, sys, problem, output, info);

		LOG_DEBUG("Reintegrating trajectory");
		for(int t=0; t < T-1; ++t) {
			J[t+1] = sys.dynfunc(J[t], U[t], vec<U_DIM>::Zero());
		}

		double max_delta_diff = -INFINITY;
		for(int t=0; t < T; ++t) {
			std::vector<Beam> fov = sys.get_fov(J[t]);
			double sd = geometry2d::signed_distance(planar_gmm[0].obj_mean, fov);
			double delta_alpha = 1.0 - 1.0/(1.0 + exp(-alpha*sd));
			double delta_inf = (sd > 0) ? 0 : 1;
			max_delta_diff = std::max(max_delta_diff, fabs(delta_alpha - delta_inf));
		}

		LOG_DEBUG(" ");
		LOG_DEBUG("Max delta difference: %4.2f", max_delta_diff);
		if (max_delta_diff < cfg::alpha_epsilon) {
			LOG_DEBUG("Max delta difference < %4.10f, exiting minimize merit", cfg::alpha_epsilon);
			break;
		}

		LOG_DEBUG("Increasing alpha by gain %4.5f", cfg::alpha_gain);
		alpha *= cfg::alpha_gain;
	}

	return cost;
}

void init_collocation(const vec<J_DIM>& j0, const mat<C_DIM,M_DIM>& P, PlanarSystem& sys,
		std::vector<vec<J_DIM>, aligned_allocator<vec<J_DIM>>>& J,
		std::vector<vec<U_DIM>, aligned_allocator<vec<U_DIM>>>& U,
		std::vector<PlanarGaussian>& planar_gmm) {
	// re-initialize GMM from PF
	sys.fit_gaussians_to_pf(P, planar_gmm);

	// find Gaussian with most particles
	// by construction of fit_gaussians_to_pf, is the first one
	vec<C_DIM> max_obj_mean = planar_gmm[0].obj_mean;

	// search for IK soln j_goal to get to max_obj_mean
	// starting with initial guess at j0
	vec<E_DIM> j_goal = j0.segment<E_DIM>(0);
	bool ik_success = sys.ik(max_obj_mean, j_goal);
	if (!ik_success) {
		LOG_ERROR("IK failed, exiting");
		exit(0);
	}

	vec<C_DIM> camera = sys.get_camera();
	// set straight-line trajectory (in joint space)
	vec<U_DIM> uinit = vec<U_DIM>::Zero();
	uinit.segment<E_DIM>(0) = (j_goal - j0.segment<E_DIM>(0)) / (double)((T-1)*DT);
	uinit(E_DIM) = atan((max_obj_mean(1) - camera(1))/(max_obj_mean(0) - camera(0))) / (double)((T-1)*DT);

	// integrate trajectory
	J[0] = j0;
	for(int t=0; t < T-1; ++t) {
		U[t] = uinit;
		J[t+1] = sys.dynfunc(J[t], U[t], vec<Q_DIM>::Zero());
	}
}

int main(int argc, char* argv[]) {
	vec<2> camera, object;
	camera << 0, 1;
	object << 5, 7; // 1, 12
	bool is_static = false;

	PlanarSystem sys = PlanarSystem(camera, object, is_static);

	// initialize starting state, belief, and pf
	mat<J_DIM,J_DIM> j_sigma0 = .01*mat<J_DIM,J_DIM>::Identity();

	vec<J_DIM> j0, j0_real;
	j0 << M_PI/5, -M_PI/2+M_PI/16, -M_PI/4, 0;
	j0_real = j0; // TODO: have them be different

	mat<C_DIM,M_DIM> P0;
	for(int m=0; m < M_DIM; ++m) {
		P0(0, m) = planar_utils::uniform(-10, 10);
		P0(1, m) = planar_utils::uniform(0, 10);
	}

	std::vector<PlanarGaussian> planar_gmm;

	// initialize state and controls
	std::vector<vec<U_DIM>, aligned_allocator<vec<U_DIM>>> U(T-1, vec<U_DIM>::Zero());
	std::vector<vec<J_DIM>, aligned_allocator<vec<J_DIM>>> J(T);

	// track real states
	std::vector<vec<J_DIM>, aligned_allocator<vec<J_DIM>>> J_real(1, j0_real);
	// particle filter
	std::vector<mat<C_DIM,M_DIM>, aligned_allocator<mat<C_DIM,M_DIM>>> pf_tracker(1, P0);

	// initialize FORCES variables
	planarMPC_params problem;
	planarMPC_output output;
	planarMPC_info info;

	setup_mpc_vars(problem, output);
	util::Timer forces_timer;

	bool stop_condition = false;
	while(!stop_condition) {
		init_collocation(j0, P0, sys,
				J, U, planar_gmm);

		LOG_INFO("Current state");
		sys.display(j0, planar_gmm);

		for(int i=0; i < planar_gmm.size(); ++i) {
			std::cout << "pct[" << i << "]: " << planar_gmm[i].pct << "\n";
		}

		LOG_INFO("Straight-line initial trajectory to Gaussian with most particles")
		sys.display(J, planar_gmm);

		// optimize
		util::Timer_tic(&forces_timer);
		double cost = planar_minimize_merit(J, U, j_sigma0, planar_gmm, sys, problem, output, info);
		double forces_time = util::Timer_toc(&forces_timer);

		LOG_INFO("Optimized cost: %4.5f", cost);
		LOG_INFO("Solve time: %5.3f ms", forces_time*1000);

		for(int t=0; t < T-1; ++t) {
			J[t+1] = sys.dynfunc(J[t], U[t], vec<Q_DIM>::Zero());
		}

		LOG_INFO("Post-optimization");
		sys.display(J, planar_gmm);

		vec<J_DIM> j_tp1_real, j_tp1;
		mat<C_DIM,M_DIM> P_tp1;
		sys.execute_control_step(J_real.back(), j0, U[0], P0,
				j_tp1_real, j_tp1, P_tp1);

		J_real.push_back(j_tp1_real);
		pf_tracker.push_back(P_tp1);

		// set start to the next time step
		j0 = j_tp1;
		P0 = P_tp1;
	}

}
