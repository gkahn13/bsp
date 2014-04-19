#include "../point-explore.h"

#include <vector>

#include "../../util/matrix.h"
#include "../../util/Timer.h"
#include "../../util/logging.h"

extern "C" {
#include "exploreMPC.h"
exploreMPC_FLOAT **H, **f, **lb, **ub, **z, **c;
}

namespace cfg {
const double improve_ratio_threshold = .1;
const double min_approx_improve = 1e-4;
const double min_trust_box_size = 1e-4;
const double trust_shrink_ratio = .5;
const double trust_expand_ratio = 1.5;
}

void setupMPCVars(exploreMPC_params& problem, exploreMPC_output& output) {
	// inputs
	H = new exploreMPC_FLOAT*[T];
	f = new exploreMPC_FLOAT*[T];
	lb = new exploreMPC_FLOAT*[T];
	ub = new exploreMPC_FLOAT*[T];
	c = new exploreMPC_FLOAT*[1];

	// output
	z = new exploreMPC_FLOAT*[T];

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
		for(int i=0; i < X_DIM+U_DIM; ++i) { H[t][i] = INFTY; }
		for(int i=0; i < X_DIM+U_DIM; ++i) { f[t][i] = INFTY; }
		for(int i=0; i < X_DIM+U_DIM; ++i) { lb[t][i] = INFTY; }
		for(int i=0; i < X_DIM+U_DIM; ++i) { ub[t][i] = INFTY; }
		for(int i=0; i < X_DIM+U_DIM; ++i) { z[t][i] = INFTY; }
	}
	for(int i=0; i < X_DIM; ++i) { H[T-1][i] = INFTY; }
	for(int i=0; i < X_DIM; ++i) { f[T-1][i] = INFTY; }
	for(int i=0; i < X_DIM; ++i) { lb[T-1][i] = INFTY; }
	for(int i=0; i < X_DIM; ++i) { ub[T-1][i] = INFTY; }
	for(int i=0; i < X_DIM; ++i) { z[T-1][i] = INFTY; }

	for(int i=0; i < X_DIM; ++i) { c[0][i] = INFTY; }
}

void cleanupMPCVars() {
	delete[] H;
	delete[] f;
	delete[] lb;
	delete[] ub;
	delete[] z;
	delete[] c;
}

bool isValidInputs()
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
	}
	for(int i=0; i < (X_DIM); ++i) { if (H[T-1][i] > INFTY/2) { return false; } }
	for(int i=0; i < (X_DIM); ++i) { if (f[T-1][i] > INFTY/2) { return false; } }
	for(int i=0; i < (X_DIM); ++i) { if (lb[T-1][i] > INFTY/2) { return false; } }
	for(int i=0; i < (X_DIM); ++i) { if (ub[T-1][i] > INFTY/2) { return false; } }

	for(int i=0; i < (X_DIM); ++i) { if (c[0][i] > INFTY/2) { return false; } }

	return true;
}

double exploreCollocation(std::vector<Matrix<X_DIM> >& X, std::vector<Matrix<U_DIM> >& U, std::vector<Matrix<X_DIM> >& P,
		exploreMPC_params& problem, exploreMPC_output& output, exploreMPC_info& info) {

	int max_iter = 100;
	double Xeps = .5;
	double Ueps = .5;

	double merit = 0;
	double constant_cost, hessian_constant, jac_constant;
	// point_platt::costfunc derivative
	Matrix<TOTAL_VARS> d, diaghess;

	std::vector<Matrix<X_DIM> > Xopt(T);
	std::vector<Matrix<U_DIM> > Uopt(T-1);
	float optcost, model_merit, new_merit;
	float approx_merit_improve, exact_merit_improve, merit_improve_ratio;

	LOG_DEBUG("Initial trajectory cost: %4.10f", point_explore::differential_entropy(X, U, P));

	int index = 0;
	bool solution_accepted = true;
	for(int it=0; it < max_iter; ++it) {

		LOG_DEBUG("\nIter: %d", it);

		// only compute gradient/hessian if P/U has been changed
		if (solution_accepted) {
			d = point_explore::grad_differential_entropy(X, U, P);

//			diaghess = point_explore::diaghess_differential_entropy(X, U, P);
			diaghess.reset();
			merit = point_explore::differential_entropy(X, U, P);

			constant_cost = 0;
			hessian_constant = 0;
			jac_constant = 0;

			// compute Hessian first so we can force it to be PSD
			// TODO: use finite differences or BFGS to approximate. (else set to 0)
			index = 0;
			for(int t=0; t < T-1; ++t) {
				for(int i=0; i < X_DIM+U_DIM; ++i) {
					double val = diaghess[index++];
					H[t][i] = (val < 0) ? 0 : val;
				}
			}
			for(int i=0; i < X_DIM; ++i) {
				double val = diaghess[index++];
				H[T-1][i] = (val < 0) ? 0 : val;
			}

			// compute gradient
			index = 0;
			for(int t=0; t < T-1; ++t) {
				Matrix<X_DIM+U_DIM> zbar;
				zbar.insert(0, 0, X[t]);
				zbar.insert(X_DIM, 0, U[t]);

				for(int i=0; i < X_DIM+U_DIM; ++i) {
					hessian_constant += H[t][i]*zbar[i]*zbar[i];
					jac_constant -= d[index]*zbar[i];
					f[t][i] = d[index] - H[t][i]*zbar[i];
					//f[t][i] = (i < M*X_DIM) ? 0 : d[index] - H[t][i]*zbar[i];
					index++;
				}
			}

			Matrix<X_DIM> zbar;
			zbar.insert(0, 0, X[T]);

			for(int i=0; i < X_DIM; ++i) {
				hessian_constant += H[T-1][i]*zbar[i]*zbar[i];
				jac_constant -= d[index]*zbar[i];
				f[T-1][i] = d[index] - H[T-1][i]*zbar[i];
				index++;
			}

			for(int i=0; i < X_DIM; ++i) {
				c[0][i] = X[0][i];
			}

			constant_cost = 0.5*hessian_constant + jac_constant + merit;
		}


		// set trust region bounds based on current trust region size
		for(int t=0; t < T; ++t) {
			// set each particle lower/upper bound
			index = 0;
			for(int i=0; i < X_DIM; ++i) {
				lb[t][index] = MAX(xMin[i], X[t][i] - Xeps);
				ub[t][index] = MIN(xMax[i], X[t][i] + Xeps);
				index++;
			}

			if (t < T-1) {
				// set each input lower/upper bound
				for(int i=0; i < U_DIM; ++i) {
					lb[t][index] = MAX(uMin[i], U[t][i] - Ueps);
					ub[t][index] = MIN(uMax[i], U[t][i] + Ueps);
					index++;
				}
			}
		}


		// Verify problem inputs
//		if (!isValidInputs()) {
//			LOG_ERROR("Inputs are not valid!");
//			exit(0);
//		}

		// call FORCES
		int exitflag = exploreMPC_solve(&problem, &output, &info);
		if (exitflag == 1) {
			optcost = info.pobj;
			for(int t=0; t < T; ++t) {
				index = 0;
				for(int i=0; i < X_DIM; ++i) {
					Xopt[t][i] = z[t][index++];
				}

				if (t < T-1) {
					for(int i=0; i < U_DIM; ++i) {
						Uopt[t][i] = z[t][index++];
					}
				}
			}
		} else {
			LOG_FATAL("Some problem in solver");
			exit(-1);
		}

		model_merit = optcost + constant_cost; // need to add constant terms that were dropped

		new_merit = point_explore::differential_entropy(Xopt, Uopt, P);

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

//		std::cout << "Displaying FORCES\n";
//		point_explore::pythonDisplayStatesAndParticles(Xopt, P, target);

//		P = Popt; U = Uopt;
//		solution_accepted = true;

		if (approx_merit_improve < -1e-5) {
			LOG_ERROR("Approximate merit function got worse: %f", approx_merit_improve);
			LOG_ERROR("Failure!");
			return INFTY;
		} else if (approx_merit_improve < cfg::min_approx_improve) {
			LOG_DEBUG("Converged: improvement small enough");
			X = Xopt; U = Uopt;
			solution_accepted = true;
			break;
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
			X = Xopt; U = Uopt;
			solution_accepted = true;
		}

//		point_explore::pythonDisplayStatesAndParticles(X, P, target);

	}


	return point_explore::differential_entropy(X, U, P);
}

int main(int argc, char* argv[]) {
	//srand(time(0));
	point_explore::initialize();

	std::vector<Matrix<X_DIM> > P(M);
	for(int m=0; m < M; ++m) {
		if (m > M/2) {
			P[m][0] = uniform(1, 2);
			P[m][1] = uniform(3, 5);
		} else {
			P[m][0] = uniform(3, 5);
			P[m][1] = uniform(1, 2);
		}
	}
//	for(int m=0; m < M; ++m) {
//		P[m][0] = uniform(xMin[0], xMax[0]);
//		P[m][1] = uniform(xMin[1], xMax[1]);
//	}

	Matrix<U_DIM> uinit = (xMax - x0) / (DT*(T-1));
	std::vector<Matrix<U_DIM> > U(T-1, uinit);
//	for(int t=0; t < T-1; ++t) {
//		for(int i=0; i < U_DIM; ++i) {
//			U[t][i] = uniform(x0[i], xMax[i]) / (DT*(T-1));
//		}
//	}

	std::vector<Matrix<X_DIM> > X(T);

	for(int t=0; t < T-1; ++t) {
		X[t+1] = point_explore::dynfunc(X[t], U[t]);
	}

	LOG_DEBUG("Display initial trajectory");
	point_explore::pythonDisplayStatesAndParticles(X, P, target);

	double init_cost = point_explore::differential_entropy(X,U,P);
	LOG_DEBUG("Initial cost: %4.10f", init_cost);

	// initialize FORCES variables
	exploreMPC_params problem;
	exploreMPC_output output;
	exploreMPC_info info;

	setupMPCVars(problem, output);

	LOG_DEBUG("Calling exploreCollocation");

	double cost = exploreCollocation(X, U, P, problem, output, info);

	LOG_INFO("Initial cost: %4.10f", init_cost);
	LOG_INFO("Cost: %4.10f", cost);

	point_explore::pythonDisplayStatesAndParticles(X,P,target);
}

