#include "include/planar-system.h"
#include "../util/Timer.h"

extern "C" {
#include "planarMPC.h"
planarMPC_FLOAT **H, **f, **lb, **ub, **z, **c;
}

#define TIMESTEPS 10
#define DT 1.0 // Note: if you change this, must change the FORCES matlab file
#define X_DIM 6
#define U_DIM 4

const int T = TIMESTEPS;
const int TOTAL_VARS = T*X_DIM + (T-1)*U_DIM;

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))
const double INFTY = 1e10;

namespace cfg {
const double alpha_init = 1; // 1
const double alpha_gain = 3; // 3
const double alpha_epsilon = .001; // .001

const double improve_ratio_threshold = .1; // .1
const double min_approx_improve = 1; // 1
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
	c = new planarMPC_FLOAT*[1];

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

	c[0] = problem.c1;

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

	for(int i=0; i < X_DIM; ++i) { c[0][i] = INFTY; }
}

void cleanup_mpc_vars() {
	delete[] H;
	delete[] f;
	delete[] lb;
	delete[] ub;
	delete[] z;
	delete[] c;
}

bool is_valid_inputs()
{
	for(int t = 0; t < T-1; ++t) {
		std::cout << "\n\nt: " << t << "\n";

		if (t == 0) {
			std::cout << "\nc[0]:\n";
			for(int i=0; i < (X_DIM); ++i) {
				std::cout << c[0][i] << " ";
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

	std::cout << "\n";

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

	for(int i=0; i < (X_DIM); ++i) { if (c[0][i] > INFTY/2) { return false; } }

	return true;
}

void L_BFGS(const std::vector<vec>& X, const std::vector<vec>& U, const vec &grad,
		const std::vector<vec> &Xopt, const std::vector<vec> &Uopt, const vec &gradopt,
		mat &hess) {
	vec s(TOTAL_VARS, fill::zeros);

	int index = 0;
	for(int t=0; t < T-1; ++t) {
		s.rows(index, index+X_DIM-1) = Xopt[t] - X[t];
		index += X_DIM;
		s.rows(index, index+U_DIM-1) = Uopt[t] - U[t];
		index += U_DIM;
	}
	s.rows(index, index+X_DIM-1) = Xopt[T-1] - X[T-1];

	mat y = gradopt - grad;

	double theta;
	vec hess_s = hess*s;

	bool decision = dot(s, y) >= .2*dot(s, hess_s);
	if (decision) {
		theta = 1;
	} else {
		theta = (.8*dot(s, hess_s))/(dot(s, hess_s) - dot(s, y));
	}

	vec r = theta*y + (1-theta)*hess_s;

	hess = hess - (hess_s*hess_s.t())/(dot(s, hess_s)) + dot(r, r)/dot(s, r);
}

double eih_collocation(std::vector<vec>& X, mat& sigma0, std::vector<vec>& U, const double alpha,
		PlanarSystem& sys, planarMPC_params &problem, planarMPC_output &output, planarMPC_info &info) {

	int max_iter = 100;
	double Xeps = .5;
	double Ueps = .5;

	double merit = 0;
	double constant_cost, hessian_constant, jac_constant;
	vec grad;
	mat hess = eye<mat>(TOTAL_VARS, TOTAL_VARS);

	std::vector<vec> Xopt(T, zeros<vec>(X_DIM));
	std::vector<vec> Uopt(T-1, zeros<vec>(U_DIM));
	vec gradopt;
	double optcost, model_merit, new_merit;
	double approx_merit_improve, exact_merit_improve, merit_improve_ratio;

	LOG_DEBUG("Initial trajectory cost: %4.10f", sys.cost(X, sigma0, U, alpha));

	int index = 0;
	bool solution_accepted = true;
	for(int it=0; it < max_iter; ++it) {

		LOG_DEBUG("");
		LOG_DEBUG("");
		LOG_DEBUG("Iter: %d", it);

		// only compute gradient/hessian if P/U has been changed
		if (solution_accepted) {

			if (it == 0) {
				grad = sys.cost_grad(X, sigma0, U, alpha);
			} else {
				grad = gradopt; // since L-BFGS calculation required it
			}

			vec diaghess = diagvec(hess);

			merit = sys.cost(X, sigma0, U, alpha);

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
				vec::fixed<(X_DIM+U_DIM)> zbar;
				zbar.subvec(0, X_DIM-1) = X[t];
				zbar.subvec(X_DIM, (X_DIM+U_DIM)-1) = U[t];

				for(int i=0; i < (X_DIM+U_DIM); ++i) {
					hessian_constant += H[t][i]*zbar(i)*zbar(i);
					jac_constant -= grad(index)*zbar(i);
					f[t][i] = grad(index) - H[t][i]*zbar(i);
					index++;
				}
			}

			vec zbar = X[T-1];

			for(int i=0; i < X_DIM; ++i) {
				hessian_constant += H[T-1][i]*zbar(i)*zbar(i);
				jac_constant -= grad(index)*zbar(i);
				f[T-1][i] = grad(index) - H[T-1][i]*zbar(i);
				index++;
			}

			for(int i=0; i < X_DIM; ++i) {
				c[0][i] = X[0](i);
			}

			constant_cost = 0.5*hessian_constant + jac_constant + merit;
		}


		vec x_min, x_max, u_min, u_max;
		sys.get_limits(x_min, x_max, u_min, u_max);

		// set trust region bounds based on current trust region size
		for(int t=0; t < T; ++t) {
			index = 0;
			for(int i=0; i < X_DIM; ++i) {
				lb[t][index] = MAX(x_min(i), X[t](i) - Xeps);
				ub[t][index] = MIN(x_max(i), X[t](i) + Xeps);
				index++;
			}


			if (t < T-1) {
				// set each input lower/upper bound
				for(int i=0; i < U_DIM; ++i) {
					lb[t][index] = MAX(u_min(i), U[t](i) - Ueps);
					ub[t][index] = MIN(u_max(i), U[t](i) + Ueps);
					index++;
				}
			}

		}


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

			gradopt = sys.cost_grad(Xopt, sigma0, Uopt, alpha);
			L_BFGS(X, U, grad, Xopt, Uopt, gradopt, hess);

			X = Xopt; U = Uopt;
			solution_accepted = true;
		}

//		LOG_DEBUG("Displaying Xopt");
//		sys.display(Xopt, sigma0, Uopt, alpha);


	}

	return sys.cost(X, sigma0, U, alpha);
}

double eih_minimize_merit(std::vector<vec>& X, mat& sigma0, std::vector<vec>& U, PlanarSystem& sys,
		planarMPC_params &problem, planarMPC_output &output, planarMPC_info &info) {
	double alpha = cfg::alpha_init;
	double cost = INFINITY;

	while(true) {
		LOG_DEBUG("Calling collocation with alpha = %4.2f", alpha);
		cost = eih_collocation(X, sigma0, U, alpha, sys, problem, output, info);

		double max_delta_diff = -INFINITY;
		for(int t=0; t < T; ++t) {
			max_delta_diff = std::max(max_delta_diff,
					max(abs(diagvec(sys.delta_matrix(X[t], alpha)) - diagvec(sys.delta_matrix(X[t], INFINITY)))));
		}

		LOG_DEBUG("");
		LOG_DEBUG("Max delta difference: %4.2f", max_delta_diff);
		if (max_delta_diff < cfg::alpha_epsilon) {
			LOG_DEBUG("Max delta difference < %4.10f, exiting minimize merit", cfg::alpha_epsilon);
			break;
		}

		LOG_DEBUG("Increasing alpha by gain %4.10f and reintegrating the trajectory", cfg::alpha_gain);
		alpha *= cfg::alpha_gain;
		for(int t=0; t < T-1; ++t) {
			X[t+1] = sys.dynfunc(X[t], U[t], zeros<vec>(U_DIM));
		}

//		LOG_DEBUG("Press enter to continue");
//		std::cin.ignore();
	}

	return cost;
}

int main(int argc, char* argv[]) {
	vec camera = {0, 1};
	vec object = {5, 8};
	bool is_static = false;

	PlanarSystem sys = PlanarSystem(camera, object, is_static);

	// initial uncertainty about object is large
	mat sigma0 = .01*eye<mat>(X_DIM, X_DIM);
	sigma0.submat(span(4,5), span(4,5)) = 20*eye<mat>(2, 2);

	vec x0 = {M_PI/5, -M_PI/2+M_PI/16, -M_PI/4, 0, 5, 5};
//	vec x0 = {M_PI/2, 0, 0, M_PI/4, 5, 5};

	// initialize state and controls
	std::vector<vec> U(T-1, zeros<vec>(U_DIM));
	std::vector<vec> X(T);
	X[0] = x0;
	for(int t=0; t < T-1; ++t) {
		X[t+1] = sys.dynfunc(X[t], U[t], zeros<vec>(U_DIM));
	}

	// initialize FORCES variables
	planarMPC_params problem;
	planarMPC_output output;
	planarMPC_info info;

	setup_mpc_vars(problem, output);
	util::Timer forces_timer;

//	double init_cost = sys.cost(X, sigma0, U, INFTY);
//	LOG_DEBUG("Initial cost %4.10f", init_cost);

//	LOG_DEBUG("Initial setup");
//	sys.display(x0, sigma0);

	double cost = eih_minimize_merit(X, sigma0, U, sys, problem, output, info);

	LOG_DEBUG("Finished minimize merit");

	std::vector<mat> S(T);
	S[0] = sigma0;
	sys.display(X[0], S[0]);
	for(int t=0; t < T-1; ++t) {
		sys.belief_dynamics(X[t], S[t], U[t], INFINITY, X[t+1], S[t+1]);
		sys.display(X[t+1], S[t+1]);
	}

}
