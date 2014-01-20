#include <vector>
#include <iostream>
#include <cstdlib>
#include <iomanip>

#include "../point.h"

#include "util/matrix.h"
//#include "util/utils.h"
#include "util/Timer.h"
#include "util/logging.h"

#include <Python.h>
//#include <pythonrun.h>
#include <boost/python.hpp>
#include <numpy/ndarrayobject.h>
#include <boost/preprocessor/iteration/local.hpp>
#include <boost/preprocessor/arithmetic/sub.hpp>
#include <boost/filesystem.hpp>

namespace py = boost::python;

extern "C" {
#include "../sym/control-symeval.h"
#include "controlMPC.h"
}

namespace cfg {
const double improve_ratio_threshold = .1;
const double min_approx_improve = 1e-2;
const double min_trust_box_size = 1e-2;
const double trust_shrink_ratio = .1;
const double trust_expand_ratio = 1.5;
}

// controlMPC vars
controlMPC_FLOAT **H, **f, **lb, **ub, **z;


void setupcontrolMPCVars(controlMPC_params& problem, controlMPC_output& output)
{
	// inputs
	H = new controlMPC_FLOAT*[T-1];
	f = new controlMPC_FLOAT*[T-1];
	lb = new controlMPC_FLOAT*[T-1];
	ub = new controlMPC_FLOAT*[T-1];

	// output
	z = new controlMPC_FLOAT*[T-1];

#define SET_VARS(n)    \
		H[ BOOST_PP_SUB(n,1) ] = problem.H##n ;  \
		f[ BOOST_PP_SUB(n,1) ] = problem.f##n ;  \
		lb[ BOOST_PP_SUB(n,1) ] = problem.lb##n ;	\
		ub[ BOOST_PP_SUB(n,1) ] = problem.ub##n ;	\
		z[ BOOST_PP_SUB(n,1) ] = output.z##n ;

#define BOOST_PP_LOCAL_MACRO(n) SET_VARS(n)
#define BOOST_PP_LOCAL_LIMITS (1, TIMESTEPS-1)
#include BOOST_PP_LOCAL_ITERATE()
}


void cleanup()
{
	delete[] inputVars;
	delete[] vars;

	delete[] H;
	delete[] f;
	delete[] lb;
	delete[] ub;
	delete[] z;
}

void initVarVals(const std::vector< Matrix<U_DIM> >& U)
{
	int idx = 0;
	for (int t = 0; t < (T - 1); ++t) {
		for (int i = 0; i < U_DIM; ++i) {
			inputVars[idx++] = U[t][i];
		}
	}
	for (int i = 0; i < (Q_DIM+R_DIM); ++i) {
		inputVars[idx++] = 0;
	}
	for (int i = 0; i < X_DIM; ++i) {
		inputVars[idx++] = x0[i];
	}
	for (int i = 0; i < X_DIM; ++i) {
		inputVars[idx++] = xGoal[i];
	}
	for (int i = 0; i < (X_DIM+X_DIM); ++i) {
		inputVars[idx++] = SqrtSigma0[i];
	}
	inputVars[idx++] = alpha_belief; inputVars[idx++] = alpha_control; inputVars[idx++] = alpha_final_belief; inputVars[idx++] = alpha_goal_state;

	int nvars = (int)maskIndices.size();
	//vars = new double[nvars];

	// For evaluation
	idx = 0;
	for (int i = 0; i < nvars; ++i) {
		vars[idx++] = inputVars[maskIndices[i]];
	}
}

inline void forcePsdHessian() {
	// zero out negative diagonal entries
	for(int t = 0; t < T-1; ++t) {
		for(int i = 0; i < U_DIM; ++i) {
			H[t][i] = (H[t][i] < 0) ? 0 : H[t][i];
		}
	}
}

bool isValidInputs(double *result)
{
		for(int t = 0; t < T-1; ++t) {

		std::cout << "t: " << t << std::endl << std::endl;

		/*
		std::cout << "lb x: " << lb[t][0] << " " << lb[t][1] << std::endl;
		std::cout << "lb u: " << lb[t][2] << " " << lb[t][3] << std::endl;

		std::cout << "ub x: " << ub[t][0] << " " << ub[t][1] << std::endl;
		std::cout << "ub u: " << ub[t][2] << " " << ub[t][3] << std::endl;



		std::cout << "f: " << std::endl;
		for(int i = 0; i < 4; ++i) {
			std::cout << f[t][i] << std::endl;
		}
		 */

		std::cout << "H: " << std::endl;
		for(int i = 0; i < U_DIM; ++i) {
			std::cout << H[t][i] << std::endl;
		}

		std::cout << std::endl << std::endl;
	}
	return true;
}


double controlCollocation(std::vector< Matrix<U_DIM> >& U, controlMPC_params& problem, controlMPC_output& output, controlMPC_info& info)
{
	int maxIter = 100;
	double Ueps = 1;

	double prevcost = INFTY, optcost;
	double merit, model_merit, new_merit;
	double approx_merit_improve, exact_merit_improve, merit_improve_ratio;
	double constant_cost, hessian_constant, jac_constant;

	int dim = (T-1)*U_DIM;
	double* resultControlCost = new double;
	double* resultControlCostGradDiagHess = new double[2*dim + 1];

	double Hubar[4];

	int nvars = (int)maskIndices.size();
	vars = new double[nvars];

	initVarVals(U);

	evalCost(resultControlCost, vars);
	prevcost = resultControlCost[0];

	std::cout << "after evalCost" << std::endl;

	LOG_DEBUG("Initialization trajectory cost: %4.10f", prevcost);

	std::vector<Matrix<U_DIM> > Uopt(T-1);

	bool solution_accepted = true;

	for(int it = 0; it < maxIter; ++it)
	{
		LOG_DEBUG("\nIter: %d", it);

		if (solution_accepted) {
			initVarVals(U);
			evalCostGradDiagHess(resultControlCostGradDiagHess, vars);
			std::cout << "after evalCostGradDiagHess" << std::endl;
			merit = resultControlCostGradDiagHess[0];

			// evaluate constant cost term (omitted from optimization)
			// Need to compute:
			// ~z_bar * H * z_bar (easy to compute since H is diagonal)
			// -f * z_bar
			constant_cost = 0;
			hessian_constant = 0;
			jac_constant = 0;

			std::cout << "compute H" << std::endl;

			// compute Hessian first
			// so can force it to be PSD
			for (int t = 0; t < T-1; ++t) {
				std::cout << t << std::endl;
				H[t][0] = resultControlCostGradDiagHess[1+dim+t*U_DIM];
				H[t][1] = resultControlCostGradDiagHess[1+dim+t*U_DIM+1];
			}

			std::cout << "before forcePsdHessian" << std::endl;

			forcePsdHessian();

			std::cout << "after forcePsdHessian" << std::endl;

			for (int t = 0; t < T-1; ++t)
			{
				Matrix<U_DIM>& ut = U[t];

				for(int i = 0; i < U_DIM; ++i) {
					Hubar[i] = H[t][i]*ut[i];
				}

				f[t][0] = resultControlCostGradDiagHess[1+t*U_DIM];
				f[t][1] = resultControlCostGradDiagHess[1+t*U_DIM+1];

				for(int i = 0; i < U_DIM; ++i) {
					hessian_constant += H[t][i]*ut[i]*ut[i];
					jac_constant += -(f[t][i]*ut[i]);
				}

				for(int i = 0; i < U_DIM; ++i) {
					f[t][i] -= Hubar[i];
				}
			} //setting up problem

			constant_cost = 0.5*hessian_constant + jac_constant + merit;
		} // end solution_accepted

		// set trust region bounds based on current trust region size
		for (int t = 0; t < T-1; ++t)
		{
			Matrix<U_DIM>& ut = U[t];

			// Fill in lb, ub, C, e
			lb[t][0] = MAX(uMin[0], ut[0] - Ueps);
			lb[t][1] = MAX(uMin[1], ut[1] - Ueps);

			ub[t][0] = MIN(uMax[0], ut[0] + Ueps);
			ub[t][1] = MIN(uMax[1], ut[1] + Ueps);
		} //setting up problem

		// Verify problem inputs
		/*
		if (!isValidInputs(result)) {
			std::cout << "Inputs are not valid!" << std::endl;
			exit(0);
		}
		 */
		std::cout << "before solve" << std::endl;
		int exitflag = controlMPC_solve(&problem, &output, &info);
		std::cout << "after solve" << std::endl;
		if (exitflag == 1) {
			for(int t = 0; t < T-1; ++t) {
				Matrix<U_DIM>& ut = Uopt[t];

				for(int i = 0; i < U_DIM; ++i) {
					ut[i] = z[t][i];
				}
				optcost = info.pobj;
			}
		}
		else {
			LOG_FATAL("Some problem in solver");
			std::exit(-1);
		}

		model_merit = optcost + constant_cost; // need to add constant terms that were dropped

		initVarVals(Uopt);
		evalCost(resultControlCost, vars);
		new_merit = resultControlCost[0];

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
			delete resultControlCost;
			delete[] resultControlCostGradDiagHess;
			return INFTY;
		} else if (approx_merit_improve < cfg::min_approx_improve) {
			LOG_DEBUG("Converged: improvement small enough");
			U = Uopt;
			solution_accepted = true;
			break;
		} else if ((exact_merit_improve < 0) || (merit_improve_ratio < cfg::improve_ratio_threshold)) {
			LOG_DEBUG("Shrinking trust region size to: %2.6f", Ueps);
			Ueps *= cfg::trust_shrink_ratio;
			solution_accepted = false;
		} else {
			LOG_DEBUG("Accepted, Increasing trust region size to: %2.6f", Ueps);
			// expand Xeps and Ueps
			Ueps *= cfg::trust_expand_ratio;
			U = Uopt;
			prevcost = optcost;
			solution_accepted = true;
		}
	}

	initVarVals(U);
	evalCost(resultControlCost, vars);
	optcost = resultControlCost[0];

	delete resultControlCost;
	delete[] resultControlCostGradDiagHess;

	return optcost;
}


int main(int argc, char* argv[])
{
	x0[0] = -3.5; x0[1] = 2;
	SqrtSigma0 = identity<X_DIM>();
	xGoal[0] = -3.5; xGoal[1] = -2;

	xMin[0] = -5; xMin[1] = -3;
	xMax[0] = 5; xMax[1] = 3;
	uMin[0] = -1; uMin[1] = -1;
	uMax[0] = 1; uMax[1] = 1;

	Matrix<U_DIM> uinit;
	uinit[0] = (xGoal[0] - x0[0]) / (T-1);
	uinit[1] = (xGoal[1] - x0[1]) / (T-1);

	std::vector<Matrix<U_DIM> > U(T-1, uinit);
	std::vector<Matrix<X_DIM> > X(T);

	X[0] = x0;
	for (int t = 0; t < T-1; ++t) {
		X[t+1] = dynfunc(X[t], U[t], zeros<Q_DIM,1>());
	}

	int nparams = 4;
	int nvars = (T - 1) * U_DIM + Q_DIM + R_DIM + X_DIM + X_DIM + (X_DIM * X_DIM) + nparams;
	std::stringstream fileName;
	fileName << "control-masks-" << TIMESTEPS << ".txt";
	setupDstarInterface(fileName.str(), nparams, nvars);

	controlMPC_params problem;
	controlMPC_output output;
	controlMPC_info info;

	setupcontrolMPCVars(problem, output);

	util::Timer solveTimer;
	util::Timer_tic(&solveTimer);

	double cost = controlCollocation(U, problem, output, info);

	double solvetime = util::Timer_toc(&solveTimer);
	LOG_INFO("Cost: %4.10f", cost);
	LOG_INFO("Solve time: %5.3f ms", solvetime*1000);

	Matrix<B_DIM> binit = zeros<B_DIM>();
	std::vector<Matrix<B_DIM> > B(T, binit);

	vec(x0, SqrtSigma0, B[0]);
	for (size_t t = 0; t < T-1; ++t) {
		B[t+1] = beliefDynamics(B[t], U[t]);
	}

	pythonDisplayTrajectory(B, U);

	cleanup();

	int k;
	std::cin >> k;

	//CAL_End();
	return 0;
}
