#include "../point-platt.h"

#include <vector>

#include "../../util/matrix.h"
#include "../../util/Timer.h"

extern "C" {
#include "plattMPC.h"
plattMPC_FLOAT **H, **f, **lb, **ub, **z;
}

namespace cfg {
const double improve_ratio_threshold = .1;
const double min_approx_improve = 1e-2;
const double min_trust_box_size = 1e-2;
const double trust_shrink_ratio = .1;
const double trust_expand_ratio = 1.5;
}

void setupMPCVars(plattMPC_params& problem, plattMPC_output& output) {
	// inputs
	H = new plattMPC_FLOAT*[T];
	f = new plattMPC_FLOAT*[T];
	lb = new plattMPC_FLOAT*[T];
	ub = new plattMPC_FLOAT*[T];

	// output
	z = new plattMPC_FLOAT*[T];

#define SET_VARS(n)    \
		H[ BOOST_PP_SUB(n,1) ] = problem.H##n ;  \
		f[ BOOST_PP_SUB(n,1) ] = problem.f##n ;  \
		lb[ BOOST_PP_SUB(n,1) ] = problem.lb##n ;	\
		ub[ BOOST_PP_SUB(n,1) ] = problem.ub##n ;	\
		z[ BOOST_PP_SUB(n,1) ] = output.z##n ;

#define BOOST_PP_LOCAL_MACRO(n) SET_VARS(n)
#define BOOST_PP_LOCAL_LIMITS (1, TIMESTEPS)
#include BOOST_PP_LOCAL_ITERATE()

	for(int t=0; t < T-1; ++t) {
		for(int i=0; i < M*X_DIM+U_DIM; ++i) { H[t][i] = INFTY; }
		for(int i=0; i < M*X_DIM+U_DIM; ++i) { f[t][i] = INFTY; }
		for(int i=0; i < M*X_DIM+U_DIM; ++i) { lb[t][i] = INFTY; }
		for(int i=0; i < M*X_DIM+U_DIM; ++i) { ub[t][i] = INFTY; }
		for(int i=0; i < M*X_DIM+U_DIM; ++i) { z[t][i] = INFTY; }
	}
	for(int i=0; i < M*X_DIM; ++i) { H[T-1][i] = INFTY; }
	for(int i=0; i < M*X_DIM; ++i) { f[T-1][i] = INFTY; }
	for(int i=0; i < M*X_DIM; ++i) { lb[T-1][i] = INFTY; }
	for(int i=0; i < M*X_DIM; ++i) { ub[T-1][i] = INFTY; }
	for(int i=0; i < M*X_DIM; ++i) { z[T-1][i] = INFTY; }

}

void cleanupMPCVars() {
	delete[] H;
	delete[] f;
	delete[] lb;
	delete[] ub;
}


bool isValidInputs()
{
	for(int t = 0; t < T-1; ++t) {
		for(int i=0; i < (M*X_DIM+U_DIM); ++i) { if (H[t][i] > INFTY/2) { return false; } }
		for(int i=0; i < (M*X_DIM+U_DIM); ++i) { if (f[t][i] > INFTY/2) { return false; } }
		for(int i=0; i < (M*X_DIM+U_DIM); ++i) { if (lb[t][i] > INFTY/2) { return false; } }
		for(int i=0; i < (M*X_DIM+U_DIM); ++i) {if (ub[t][i] > INFTY/2) { return false; } }
	}
	for(int i=0; i < (M*X_DIM+U_DIM); ++i) { if (H[T-1][i] > INFTY/2) { return false; } }
	for(int i=0; i < (M*X_DIM+U_DIM); ++i) { if (f[T-1][i] > INFTY/2) { return false; } }
	for(int i=0; i < (M*X_DIM+U_DIM); ++i) { if (lb[T-1][i] > INFTY/2) { return false; } }
	for(int i=0; i < (M*X_DIM+U_DIM); ++i) { if (ub[T-1][i] > INFTY/2) { return false; } }

	return true;
}

double plattCollocation(std::vector<std::vector<Matrix<X_DIM> > >& P, std::vector<Matrix<U_DIM> >& U,
		plattMPC_params& problem, plattMPC_output& output, plattMPC_info& info) {

	int max_iter = 100;
	double Xeps = 1;
	double Ueps = 1;

	double delta = .01;

	double merit = 0;
	double constant_cost, hessian_constant, jac_constant;
	// costfunc derivative
	Matrix<TOTAL_VARS> d;

	std::vector<std::vector<Matrix<X_DIM> > > Popt(T);
	for(int t=0; t < T; ++t) {
		std::vector<Matrix<X_DIM> > P_t(M);
		Popt[t] = P_t;
	}
	std::vector<Matrix<U_DIM> > Uopt(T-1);
	float optcost, model_merit, new_merit;
	float approx_merit_improve, exact_merit_improve, merit_improve_ratio;

	LOG_DEBUG("Initial trajectory cost: %4.10f", costfunc(P,U));

	int index = 0;
	bool solution_accepted = true;
	for(int it=0; it < max_iter; ++it) {

		LOG_DEBUG("\nIter: %d", it);

		// only compute gradient/hessian if P/U has been changed
		if (solution_accepted) {
			d = deriv_costfunc(P, U);
			merit = costfunc(P, U);

			constant_cost = 0;
			hessian_constant = 0;
			jac_constant = 0;

			// compute Hessian first so we can force it to be PSD
			// TODO: use BFGS to approximate. setting to 0 for now
			for(int t=0; t < T-1; ++t) {
				for(int i=0; i < M*X_DIM+U_DIM; ++i) {
					H[t][i] = 0;
				}
			}
			for(int i=0; i < M*X_DIM; ++i) {
				H[T-1][i] = 0;
			}

			// compute gradient
			for(int t=0; t < T; ++t) {
				int stage_vars = (t < T-1) ? (M*X_DIM+U_DIM) : (M*X_DIM);

				Matrix<stage_vars> zbar;
				for(int m=0; m < M; ++m) { zbar.insert(m*X_DIM, 0, P[t][m]); }
				if (t < T-1) { zbar.insert(M*X_DIM, 0, U[t]); }

				index = 0;
				for(int i=0; i < stage_vars; ++i) {
					hessian_constant += H[t][i]*zbar[i]*zbar[i];
					jac_constant -= d[index]*zbar[i];
					f[t][i] = d[index] - H[t][i]*zbar[i];
					index++;
				}
			}

			constant_cost = 0.5*hessian_constant + jac_constant + merit;
		}


		// set trust region bounds based on current trust region size
		for(int t=0; t < T; ++t) {
			// set each particle lower/upper bound
			index = 0;
			for(int m=0; m < M; ++m) {
				for(int i=0; i < X_DIM; ++i) {
					lb[t][index] = MAX(xMin[i], P[t][m][i] - Xeps);
					ub[t][index] = MIN(xMax[i], P[t][m][i] + Xeps);
					index++;
				}
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

		// set goal constraint for MLE particle
		for(int i=0; i < X_DIM; ++i) {
			lb[T-1][i] = MAX(xGoal[i] - delta, P[T-1][0][i] - Xeps);
			ub[T-1][i] = MIN(xGoal[i] + delta, P[T-1][0][i] + Xeps);
		}

		// Verify problem inputs
		if (!isValidInputs()) {
			LOG_ERROR("Inputs are not valid!");
			exit(0);
		}

		// call FORCES
		int exitflag = plattMPC_solve(&problem, &output, &info);
		if (exitflag == 1) {
			optcost = info.pobj;
			for(int t=0; t < T; ++t) {
				index = 0;
				for(int m=0; m < M; ++m) {
					for(int i=0; i < X_DIM; ++i) {
						Popt[t][m][i] = z[t][index++];
					}
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

		new_merit = costfunc(Popt, Uopt);

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
			P = Popt; U = Uopt;
			solution_accepted = true;
			break;
		} else if ((exact_merit_improve < 0) || (merit_improve_ratio < cfg::improve_ratio_threshold)) {
			LOG_DEBUG("Shrinking trust region size to: %2.6f %2.6f", Xeps, Ueps);
			Xeps *= cfg::trust_shrink_ratio;
			Ueps *= cfg::trust_shrink_ratio;
			solution_accepted = false;
		} else {
			LOG_DEBUG("Accepted, Increasing trust region size to:  %2.6f %2.6f", Xeps, Ueps);
			// expand Xeps and Ueps
			Xeps *= cfg::trust_expand_ratio;
			Ueps *= cfg::trust_expand_ratio;
			P = Popt; U = Uopt;
			solution_accepted = true;
		}

	}


	return 0;
}


