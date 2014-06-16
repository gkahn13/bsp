#include "include/planar-system.h"
#include "../util/Timer.h"

extern "C" {
#include "planarMPC.h"
planarMPC_FLOAT **H, **f, **lb, **ub, **z, *c, *A, *b;
}

const int T = TIMESTEPS;
const double INFTY = 1e10;

const bool USE_FADBAD = true;

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

	for(int i=0; i < X_DIM; ++i) { c[i] = INFTY; }
	for(int i=0; i < (2*C_DIM)*X_DIM; ++i) { A[i] = INFTY; }
	for(int i=0; i < (2*C_DIM); ++i) { b[i] = INFTY; }

	for(int t=0; t < T-1; ++t) {
		for(int i=0; i < (X_DIM+U_DIM); ++i) { H[t][i] = INFTY; }
		for(int i=0; i < (X_DIM+U_DIM); ++i) { f[t][i] = INFTY; }
		for(int i=0; i < (X_DIM+U_DIM); ++i) { lb[t][i] = INFTY; }
		for(int i=0; i < (X_DIM+U_DIM); ++i) { ub[t][i] = INFTY; }
		for(int i=0; i < (X_DIM+U_DIM); ++i) { z[t][i] = INFTY; }
	}
	for(int i=0; i < X_DIM; ++i) { H[T-1][i] = INFTY; }
	for(int i=0; i < X_DIM; ++i) { f[T-1][i] = INFTY; }
	for(int i=0; i < X_DIM; ++i) { lb[T-1][i] = INFTY; }
	for(int i=0; i < X_DIM; ++i) { ub[T-1][i] = INFTY; }
	for(int i=0; i < X_DIM; ++i) { z[T-1][i] = INFTY; }

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
			for(int i=0; i < (X_DIM); ++i) {
				std::cout << c[i] << " ";
			}
		}

		std::cout << "\nH[" << t << "]: ";
		for(int i=0; i < (X_DIM+U_DIM); ++i) {
			std::cout << H[t][i] << " ";
		}

		std::cout << "\nf[" << t << "]: ";
		for(int i=0; i < (X_DIM+U_DIM); ++i) {
			std::cout << f[t][i] << " ";
		}

		std::cout << "\nlb[" << t << "]: ";
		for(int i=0; i < (X_DIM+U_DIM); ++i) {
			std::cout << lb[t][i] << " ";
		}

		std::cout << "\nub[" << t << "]: ";
		for(int i=0; i < (X_DIM+U_DIM); ++i) {
			std::cout << ub[t][i] << " ";
		}
	}
	std::cout << "\n\nt: " << T-1 << "\n";

	std::cout << "\nH[" << T-1 << "]: ";
	for(int i=0; i < (X_DIM); ++i) {
		std::cout << H[T-1][i] << " ";
	}

	std::cout << "\nf[" << T-1 << "]: ";
	for(int i=0; i < (X_DIM); ++i) {
		std::cout << f[T-1][i] << " ";
	}

	std::cout << "\nlb[" << T-1 << "]: ";
	for(int i=0; i < (X_DIM); ++i) {
		std::cout << lb[T-1][i] << " ";
	}

	std::cout << "\nub[" << T-1 << "]: ";
	for(int i=0; i < (X_DIM); ++i) {
		std::cout << ub[T-1][i] << " ";
	}

	std::cout << "\nA[" << T-1 << "]:\n";
	for(int i=0; i < 2*C_DIM; ++i) {
		for(int j=0; j < X_DIM; ++j) {
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

	for(int i=0; i < (X_DIM); ++i) { if (c[i] > INFTY/2) { return false; } }
	for(int i=0; i < (2*C_DIM)*X_DIM; ++i) { if (A[i] > INFTY/2) { return false; } }
	for(int i=0; i < (2*C_DIM); ++i) { if (b[i] > INFTY/2) { return false; } }

	for(int t = 0; t < T-1; ++t) {
		for(int i=0; i < (X_DIM+U_DIM); ++i) { if (H[t][i] > INFTY/2) { return false; } }
		for(int i=0; i < (X_DIM+U_DIM); ++i) { if (f[t][i] > INFTY/2) { return false; } }
		for(int i=0; i < (X_DIM+U_DIM); ++i) { if (lb[t][i] > INFTY/2) { return false; } }
		for(int i=0; i < (X_DIM+U_DIM); ++i) {if (ub[t][i] > INFTY/2) { return false; } }
		for(int i=0; i < (X_DIM+U_DIM); ++i) {if (ub[t][i] < lb[t][i]) { return false; } }
	}
	for(int i=0; i < (X_DIM); ++i) { if (H[T-1][i] > INFTY/2) { return false; } }
	for(int i=0; i < (X_DIM); ++i) { if (f[T-1][i] > INFTY/2) { return false; } }
	for(int i=0; i < (X_DIM); ++i) { if (lb[T-1][i] > INFTY/2) { return false; } }
	for(int i=0; i < (X_DIM); ++i) { if (ub[T-1][i] > INFTY/2) { return false; } }
	for(int i=0; i < (X_DIM); ++i) {if (ub[T-1][i] < lb[T-1][i]) { return false; } }

	return true;
}

void L_BFGS(const std::vector<vec<X_DIM>, aligned_allocator<vec<X_DIM>>>& X,
		const std::vector<vec<U_DIM>, aligned_allocator<vec<U_DIM>>>& U, const vec<TOTAL_VARS> &grad,
		const std::vector<vec<X_DIM>, aligned_allocator<vec<X_DIM>>> &Xopt,
		const std::vector<vec<U_DIM>, aligned_allocator<vec<U_DIM>>> &Uopt, const vec<TOTAL_VARS> &gradopt,
		mat<TOTAL_VARS, TOTAL_VARS> &hess) {
	vec<TOTAL_VARS> s = vec<TOTAL_VARS>::Zero();

	int index = 0;
	for(int t=0; t < T-1; ++t) {
		s.segment<X_DIM>(index) = Xopt[t] - X[t];
		index += X_DIM;
		s.segment<U_DIM>(index) = Uopt[t] - U[t];
		index += U_DIM;
	}
	s.segment<X_DIM>(index) = Xopt[T-1] - X[T-1];

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

double planar_collocation(std::vector<vec<X_DIM>, aligned_allocator<vec<X_DIM>>>& X, mat<X_DIM,X_DIM>& sigma0,
		std::vector<vec<U_DIM>, aligned_allocator<vec<U_DIM>>>& U, const double alpha,
		PlanarSystem& sys, planarMPC_params &problem, planarMPC_output &output, planarMPC_info &info) {

	int max_iter = 100;
	double Xeps = .5;
	double Ueps = .5;

	double merit = 0, meritopt = 0;
	double constant_cost, hessian_constant, jac_constant;
	vec<TOTAL_VARS> grad = vec<TOTAL_VARS>::Zero();
	mat<TOTAL_VARS,TOTAL_VARS> hess = mat<TOTAL_VARS,TOTAL_VARS>::Identity();

	std::vector<vec<X_DIM>, aligned_allocator<vec<X_DIM>>> Xopt(T, vec<X_DIM>::Zero());
	std::vector<vec<U_DIM>, aligned_allocator<vec<U_DIM>>> Uopt(T-1, vec<U_DIM>::Zero());
	vec<TOTAL_VARS> gradopt = vec<TOTAL_VARS>::Zero();
	double optcost, model_merit, new_merit;
	double approx_merit_improve, exact_merit_improve, merit_improve_ratio;

	LOG_DEBUG("Initial trajectory cost: %4.10f", sys.cost(X, sigma0, U, alpha));

	int index = 0;
	bool solution_accepted = true;
	for(int it=0; it < max_iter; ++it) {

		LOG_DEBUG(" ");
		LOG_DEBUG(" ");
		LOG_DEBUG("Iter: %d", it);

		// only compute gradient/hessian if P/U has been changed
		if (solution_accepted) {

			if (it == 0) {
//				grad = sys.cost_grad(X, sigma0, U, alpha, USE_FADBAD);
				sys.cost_and_cost_grad(X, sigma0, U, alpha, USE_FADBAD, merit, grad);
			} else {
				merit = meritopt; // since L-BFGS calculation required it
				grad = gradopt;
			}

			vec<TOTAL_VARS> diaghess = hess.diagonal();
//			merit = sys.cost(X, sigma0, U, alpha);

			constant_cost = 0;
			hessian_constant = 0;
			jac_constant = 0;

			// fill in Hessian first so we can force it to be PSD
			index = 0;
			for(int t=0; t < T-1; ++t) {
				for(int i=0; i < (X_DIM+U_DIM); ++i) {
					double val = diaghess(index++);
					H[t][i] = (val < 0) ? 0 : val;
				}
			}
			for(int i=0; i < X_DIM; ++i) {
				double val = diaghess(index++);
				H[T-1][i] = (val < 0) ? 0 : val;
			}

			// fill in gradient
			index = 0;
			for(int t=0; t < T-1; ++t) {
				vec<(X_DIM+U_DIM)> zbar;
				zbar.segment<X_DIM>(0) = X[t];
				zbar.segment<U_DIM>(X_DIM) = U[t];

				for(int i=0; i < (X_DIM+U_DIM); ++i) {
					hessian_constant += H[t][i]*zbar(i)*zbar(i);
					jac_constant -= grad(index)*zbar(i);
					f[t][i] = grad(index) - H[t][i]*zbar(i);
					index++;
				}
			}

			vec<X_DIM> zbar = X[T-1];

			for(int i=0; i < X_DIM; ++i) {
				hessian_constant += H[T-1][i]*zbar(i)*zbar(i);
				jac_constant -= grad(index)*zbar(i);
				f[T-1][i] = grad(index) - H[T-1][i]*zbar(i);
				index++;
			}

			for(int i=0; i < X_DIM; ++i) {
				c[i] = X[0](i);
			}

			constant_cost = 0.5*hessian_constant + jac_constant + merit;
		}


		vec<X_DIM> x_min, x_max;
		vec<U_DIM> u_min, u_max;
		sys.get_limits(x_min, x_max, u_min, u_max);

		// set trust region bounds based on current trust region size
		for(int t=0; t < T; ++t) {
			index = 0;
			for(int i=0; i < X_DIM; ++i) {
				lb[t][index] = std::max(x_min(i), X[t](i) - Xeps);
				ub[t][index] = std::min(x_max(i), X[t](i) + Xeps);
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
		vec<C_DIM> goal_pos = X[T-1].segment<2>(J_DIM);
		double goal_delta = .01;

		vec<E_DIM> jT = X[T-1].segment<E_DIM>(0);
		vec<C_DIM> ee_pos = sys.get_ee_pos(jT);
		mat<C_DIM,E_DIM> ee_pos_jac;
		sys.get_ee_pos_jac(jT, ee_pos_jac);
		mat<C_DIM,X_DIM> ee_pos_jac_full = mat<C_DIM,X_DIM>::Zero();
		ee_pos_jac_full.block<C_DIM,E_DIM>(0,0) = ee_pos_jac;

		mat<2*C_DIM,X_DIM> Amat;
		Amat.block<C_DIM,X_DIM>(0,0) = ee_pos_jac_full;
		Amat.block<C_DIM,X_DIM>(C_DIM,0) = -ee_pos_jac_full;
		fill_col_major(A, Amat);

		vec<2*C_DIM> bVec;
		bVec.segment<C_DIM>(0) = goal_pos - ee_pos + ee_pos_jac_full*X[T-1] + goal_delta*vec<C_DIM>::Ones();
		bVec.segment<C_DIM>(C_DIM) = -goal_pos + ee_pos - ee_pos_jac_full*X[T-1] + goal_delta*vec<C_DIM>::Ones();
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
				for(int i=0; i < X_DIM; ++i) {
					Xopt[t](i) = z[t][index++];
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

		new_merit = sys.cost(Xopt, sigma0, Uopt, alpha);

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
			X = Xopt; U = Uopt;
			solution_accepted = true;
			return sys.cost(X, sigma0, U, alpha);
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

			sys.cost_and_cost_grad(Xopt, sigma0, Uopt, alpha, USE_FADBAD, meritopt, gradopt);
			L_BFGS(X, U, grad, Xopt, Uopt, gradopt, hess);

			X = Xopt; U = Uopt;
			solution_accepted = true;
		}

//		LOG_DEBUG("Displaying Xopt");
//		sys.display(Xopt, sigma0, Uopt, alpha);


	}

	return sys.cost(X, sigma0, U, alpha);
}

double planar_minimize_merit(std::vector<vec<X_DIM>, aligned_allocator<vec<X_DIM>>>& X, mat<X_DIM,X_DIM>& sigma0,
		std::vector<vec<U_DIM>, aligned_allocator<vec<U_DIM>>>& U, PlanarSystem& sys,
		planarMPC_params &problem, planarMPC_output &output, planarMPC_info &info) {
	double alpha = cfg::alpha_init;
	double cost = INFINITY;

	while(true) {
		LOG_DEBUG("Calling collocation with alpha = %4.2f", alpha);
		cost = planar_collocation(X, sigma0, U, alpha, sys, problem, output, info);

		LOG_DEBUG("Reintegrating trajectory");
		for(int t=0; t < T-1; ++t) {
			X[t+1] = sys.dynfunc(X[t], U[t], vec<U_DIM>::Zero());
		}

		double max_delta_diff = -INFINITY;
		for(int t=0; t < T; ++t) {
			mat<X_DIM,X_DIM> delta_alpha = sys.delta_matrix(X[t], X[t].segment<C_DIM>(J_DIM), alpha);
			mat<X_DIM,X_DIM> delta_inf = sys.delta_matrix(X[t], X[t].segment<C_DIM>(J_DIM), INFINITY);
			max_delta_diff = std::max(max_delta_diff, (delta_alpha - delta_inf).array().abs().maxCoeff());
		}

		LOG_DEBUG(" ");
		LOG_DEBUG("Max delta difference: %4.2f", max_delta_diff);
		if (max_delta_diff < cfg::alpha_epsilon) {
			LOG_DEBUG("Max delta difference < %4.10f, exiting minimize merit", cfg::alpha_epsilon);
			break;
		}

		LOG_DEBUG("Increasing alpha by gain %4.5f", cfg::alpha_gain);
		alpha *= cfg::alpha_gain;

//		LOG_DEBUG("Press enter to continue");
//		std::cin.ignore();
	}

	return cost;
}

int main(int argc, char* argv[]) {
	vec<2> camera, object;
	camera << 0, 1;
	object << 5, 7; // 1, 12
	bool is_static = false;

	PlanarSystem sys = PlanarSystem(camera, object, is_static);

	// initialize starting state, belief, and pf
	mat<X_DIM,X_DIM> sigma0 = .01*mat<X_DIM,X_DIM>::Identity();
	sigma0.block<C_DIM,C_DIM>(J_DIM, J_DIM) = 20*mat<C_DIM,C_DIM>::Identity(); // uncertainty about object large

	vec<X_DIM> x0, x0_real;
	x0 << M_PI/5, -M_PI/2+M_PI/16, -M_PI/4, 0, -3, 5;
	x0_real << x0.segment<J_DIM>(0) , object;

	mat<C_DIM,M_DIM> P0;
	for(int m=0; m < M_DIM; ++m) {
		P0(0, m) = planar_utils::uniform(-10, 10);
		P0(1, m) = planar_utils::uniform(0, 10);
	}

	vec<C_DIM> new_mean;
	mat<C_DIM,C_DIM> new_cov;
	sys.reinitialize(x0, P0, new_mean, new_cov);
	x0.segment<C_DIM>(J_DIM) = new_mean;
	sigma0.block<C_DIM,C_DIM>(J_DIM,J_DIM) = new_cov;

	LOG_INFO("Initial state");
	sys.display(x0, sigma0, P0);

	// initialize state and controls
	std::vector<vec<U_DIM>, aligned_allocator<vec<U_DIM>>> U(T-1, vec<U_DIM>::Zero());
	std::vector<vec<X_DIM>, aligned_allocator<vec<X_DIM>>> X(T);

	// track real states/beliefs
	std::vector<vec<X_DIM>, aligned_allocator<vec<X_DIM>>> X_real(1, x0_real);
	std::vector<mat<X_DIM,X_DIM>, aligned_allocator<mat<X_DIM,X_DIM>>> S_real(1, sigma0);
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
		// initialize straight-line trajectory
		// integrate dynamics
		vec<C_DIM> ee_goal = x0.segment<C_DIM>(J_DIM);
		vec<E_DIM> j_goal = x0.segment<E_DIM>(0);
		bool ik_success = sys.ik(ee_goal, j_goal);
		if (!ik_success) {
			LOG_ERROR("IK failed, exiting");
			exit(0);
		}

		vec<U_DIM> uinit = vec<U_DIM>::Zero();
		uinit.segment<E_DIM>(0) = (j_goal - x0.segment<E_DIM>(0)) / (double)((T-1)*DT);
		uinit(E_DIM) = atan((x0(X_DIM-1) - camera(1))/(x0(X_DIM-2) - camera(0))) / (double)((T-1)*DT);

		X[0] = x0;
		for(int t=0; t < T-1; ++t) {
			U[t] = uinit;
			X[t+1] = sys.dynfunc(X[t], U[t], vec<Q_DIM>::Zero());
		}

//		LOG_INFO("Initialized path");
//		sys.display(X, sigma0, U, pf_tracker.back(), 10);

		double init_cost = sys.cost(X, sigma0, U, INFINITY);
		LOG_INFO("Initial cost: %4.5f", init_cost);

		// optimize
		util::Timer_tic(&forces_timer);
		double cost = planar_minimize_merit(X, sigma0, U, sys, problem, output, info);
		double forces_time = util::Timer_toc(&forces_timer);

		LOG_INFO("Optimized cost: %4.5f", cost);
		LOG_INFO("Solve time: %5.3f ms", forces_time*1000);

		for(int t=0; t < T-1; ++t) {
			X[t+1] = sys.dynfunc(X[t], U[t], vec<Q_DIM>::Zero());
		}
//		LOG_INFO("Integrated trajectory");
//		sys.display(X, sigma0, U, pf_tracker.back(), INFINITY);

		// execute first control input
		vec<X_DIM> x_tp1_real, x_tp1_tp1;
		mat<X_DIM,X_DIM> sigma_tp1_tp1;
		mat<C_DIM,M_DIM> P_tp1;
		sys.execute_control_step(X_real.back(), x0, sigma0, U[0], pf_tracker.back(), x_tp1_real, x_tp1_tp1, sigma_tp1_tp1, P_tp1);

		X_real.push_back(x_tp1_real);
		S_real.push_back(sigma_tp1_tp1);
		pf_tracker.push_back(P_tp1);

//		LOG_DEBUG("Display after obtaining observation, but before truncating the belief");
//		sys.display(x_tp1_tp1, sigma_tp1_tp1, pf_tracker.back());

		// truncate belief
		std::vector<Beam> fov_tp1 = sys.get_fov(x_tp1_tp1);
		bool received_obs = geometry2d::is_inside(object, fov_tp1);
		vec<C_DIM> obj_pos_trunc;
		mat<C_DIM,C_DIM> obj_sigma_trunc;
		geometry2d::my_truncate_belief(fov_tp1, x_tp1_tp1.segment<C_DIM>(J_DIM),
				sigma_tp1_tp1.block<C_DIM,C_DIM>(J_DIM, J_DIM), received_obs,
				obj_pos_trunc, obj_sigma_trunc);

		vec<X_DIM> x_tp1_tp1_trunc;
		x_tp1_tp1_trunc << x_tp1_tp1.segment<J_DIM>(0), obj_pos_trunc;
		mat<X_DIM,X_DIM> sigma_tp1_tp1_trunc = sigma_tp1_tp1;
		sigma_tp1_tp1_trunc.block<C_DIM,C_DIM>(J_DIM, J_DIM) = obj_sigma_trunc;

		LOG_DEBUG("Display truncated belief");
		sys.display(x_tp1_tp1_trunc, sigma_tp1_tp1_trunc, pf_tracker.back());

//		stop_condition = sys.should_reinitialize(x_tp1_tp1_trunc, sigma_tp1_tp1_trunc, pf_tracker.back());
//		if (stop_condition) {
//			LOG_INFO("Stop condition met. Reinitializing");
//			vec<C_DIM> obj_mean_new;
//			mat<C_DIM,C_DIM> obj_cov_new;
//			sys.reinitialize(pf_tracker.back(), obj_mean_new, obj_cov_new);
//			x_tp1_tp1_trunc.segment<C_DIM>(J_DIM) = obj_mean_new;
//			sigma_tp1_tp1_trunc.block<C_DIM,C_DIM>(J_DIM,J_DIM) = obj_cov_new;
//			stop_condition = false;
//		}

		sys.reinitialize(x_tp1_tp1_trunc, pf_tracker.back(), new_mean, new_cov);
		x_tp1_tp1_trunc.segment<C_DIM>(J_DIM) = new_mean;
		sigma_tp1_tp1_trunc.block<C_DIM,C_DIM>(J_DIM,J_DIM) = new_cov;

		LOG_DEBUG("Display reinitialized belief");
		sys.display(x_tp1_tp1_trunc, sigma_tp1_tp1_trunc, pf_tracker.back());

		// set start to the next time step
		x0 = x_tp1_tp1_trunc;
		sigma0 = sigma_tp1_tp1_trunc;


		// no truncation
		//		x0 = x_tp1_tp1;
		//		sigma0 = sigma_tp1_tp1;

	}

	LOG_INFO("Stop condition met, exiting");

}
