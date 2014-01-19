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
#include <boost/python.hpp>
#include <boost/filesystem.hpp>

namespace py = boost::python;

extern "C" {
#include "../sym/state-symeval.h"
#include "lpMPC.h"
}

namespace cfg {
const double improve_ratio_threshold = .2;
const double min_approx_improve = 1e-2;
const double min_trust_box_size = 1e-2;
const double trust_shrink_ratio = .1;
const double trust_expand_ratio = 1.5;
}

// lpMPC vars
lpMPC_FLOAT **f, **lb, **ub, **C, **e, **z;



void setupLpMPCVars(lpMPC_params& problem, lpMPC_output& output)
{
	// inputs
	f = new lpMPC_FLOAT*[T];
	lb = new lpMPC_FLOAT*[T];
	ub = new lpMPC_FLOAT*[T];
	C = new lpMPC_FLOAT*[T-1];
	e = new lpMPC_FLOAT*[T];

	// output
	z = new lpMPC_FLOAT*[T];

#define SET_VARS(n)    \
		f[ BOOST_PP_SUB(n,1) ] = problem.f##n ;  \
		lb[ BOOST_PP_SUB(n,1) ] = problem.lb##n ;	\
		ub[ BOOST_PP_SUB(n,1) ] = problem.ub##n ;	\
		C[ BOOST_PP_SUB(n,1) ] = problem.C##n ; \
		e[ BOOST_PP_SUB(n,1) ] = problem.e##n ;  \
		z[ BOOST_PP_SUB(n,1) ] = output.z##n ;

#define BOOST_PP_LOCAL_MACRO(n) SET_VARS(n)
#define BOOST_PP_LOCAL_LIMITS (1, TIMESTEPS-1)
#include BOOST_PP_LOCAL_ITERATE()

#define SET_LAST_VARS(n)    \
		f[ BOOST_PP_SUB(n,1) ] = problem.f##n ;  \
		lb[ BOOST_PP_SUB(n,1) ] = problem.lb##n ;	\
		ub[ BOOST_PP_SUB(n,1) ] = problem.ub##n ;	\
		e[ BOOST_PP_SUB(n,1) ] = problem.e##n ;  \
		z[ BOOST_PP_SUB(n,1) ] = output.z##n ;

#define BOOST_PP_LOCAL_MACRO(n) SET_LAST_VARS(n)
#define BOOST_PP_LOCAL_LIMITS (TIMESTEPS, TIMESTEPS)
#include BOOST_PP_LOCAL_ITERATE()

}

void cleanup()
{
	delete[] inputVars;
	delete[] vars;

	delete[] f;
	delete[] lb; 
	delete[] ub; 
	delete[] C;
	delete[] e;
	delete[] z;
}

void initVarVals(const std::vector< Matrix<X_DIM> >& X, const std::vector< Matrix<U_DIM> >& U)
{
	int idx = 0;	
	for (int t = 0; t < T; ++t) {
		for (int i = 0; i < X_DIM; ++i) {
			inputVars[idx++] = X[t][i];
		}
	}
	for (int t = 0; t < (T - 1); ++t) {
		for (int i = 0; i < U_DIM; ++i) {
			inputVars[idx++] = U[t][i];
		}
	}
	for (int i = 0; i < (Q_DIM+R_DIM); ++i) {
		inputVars[idx++] = 0;
	}
	for (int i = 0; i < (X_DIM+X_DIM); ++i) {
		inputVars[idx++] = SqrtSigma0[i];
	}
	inputVars[idx++] = alpha_belief; inputVars[idx++] = alpha_control; inputVars[idx++] = alpha_final_belief;

	int nvars = (int)maskIndices.size();

	// For evaluation
	idx = 0;
	for (int i = 0; i < nvars; ++i) {
		vars[idx++] = inputVars[maskIndices[i]];
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

		std::cout << std::endl << std::endl;
	}
	return true;
}

double lpCollocation(std::vector< Matrix<X_DIM> >& X, std::vector< Matrix<U_DIM> >& U, lpMPC_params& problem, lpMPC_output& output, lpMPC_info& info)
{
	int maxIter = 100;
	double Xeps = 1;
	double Ueps = 1;

	// box constraint around goal
	double delta = 0.01;

	Matrix<X_DIM,1> x0 = X[0];

	double prevcost = INFTY, optcost;
	double merit, model_merit, new_merit;
	double approx_merit_improve, exact_merit_improve, merit_improve_ratio;
	double constant_cost, jac_constant;

	int dim = T*X_DIM + (T-1)*U_DIM;
	double* resultCost = new double;
	double* resultCostGrad = new double[dim + 1];

	int nvars = (int)maskIndices.size();
	vars = new double[nvars];

	initVarVals(X, U);
	evalCost(resultCost, vars);
	prevcost = resultCost[0];

	LOG_DEBUG("Initialization trajectory cost: %4.10f", prevcost);

	std::vector<Matrix<X_DIM> > Xopt(T);
	std::vector<Matrix<U_DIM> > Uopt(T-1);

	bool solution_accepted = true;

	for(int it = 0; it < maxIter; ++it)
	{
		LOG_DEBUG("\nIter: %d", it);

		if (solution_accepted)
		{
			initVarVals(X, U);
			evalCostGrad(resultCostGrad, vars);
			merit = resultCostGrad[0];

			// evaluate constant cost term (omitted from optimization)
			// Need to compute:
			// -f * z_bar
			constant_cost = 0;
			jac_constant = 0;

			for (int t = 0; t < T-1; ++t)
			{
				Matrix<X_DIM>& xt = X[t];
				Matrix<U_DIM>& ut = U[t];

				Matrix<X_DIM+U_DIM> zbar;
				zbar.insert(0,0,xt);
				zbar.insert(X_DIM,0,ut);

				f[t][0] = resultCostGrad[1+t*X_DIM];
				f[t][1] = resultCostGrad[1+t*X_DIM+1];
				f[t][2] = resultCostGrad[1+T*X_DIM+t*U_DIM];
				f[t][3] = resultCostGrad[1+T*X_DIM+t*U_DIM+1];

				for(int i = 0; i < (X_DIM+U_DIM); ++i) {
					jac_constant += -(f[t][i]*zbar[i]);
				}

				Matrix<X_DIM,X_DIM+U_DIM> CMat;
				Matrix<X_DIM> eVec;

				CMat.insert<X_DIM,X_DIM>(0,0,identity<X_DIM>());
				CMat.insert<X_DIM,U_DIM>(0,X_DIM,DT*identity<U_DIM>());
				int idx = 0;
				for(int c = 0; c < (X_DIM+U_DIM); ++c) {
					for(int r = 0; r < X_DIM; ++r) {
						C[t][idx++] = CMat[c + r*(X_DIM+U_DIM)];
					}
				}

				if (t == 0) {
					e[t][0] = x0[0]; e[t][1] = x0[1];
				} else {
					e[t][0] = 0; e[t][1] = 0;
				}
			} //setting up problem

			// Last stage
			Matrix<X_DIM>& xT = X[T-1];

			f[T-1][0] = resultCostGrad[1+(T-1)*X_DIM];
			f[T-1][1] = resultCostGrad[1+(T-1)*X_DIM+1];

			for(int i = 0; i < X_DIM; ++i) {
				jac_constant += -(f[T-1][i]*xT[i]);
			}

			e[T-1][0] = 0; e[T-1][1] = 0;

			constant_cost = jac_constant + resultCostGrad[0];
		} // end solution_accepted

		// set trust region bounds based on current trust region size
		for (int t = 0; t < T-1; ++t)
		{
			Matrix<X_DIM>& xt = X[t];
			Matrix<U_DIM>& ut = U[t];

			// Fill in lb, ub, C, e
			lb[t][0] = MAX(xMin[0], xt[0] - Xeps);
			lb[t][1] = MAX(xMin[1], xt[1] - Xeps);
			lb[t][2] = MAX(uMin[0], ut[0] - Ueps);
			lb[t][3] = MAX(uMin[1], ut[1] - Ueps);

			ub[t][0] = MIN(xMax[0], xt[0] + Xeps);
			ub[t][1] = MIN(xMax[1], xt[1] + Xeps);
			ub[t][2] = MIN(uMax[0], ut[0] + Ueps);
			ub[t][3] = MIN(uMax[1], ut[1] + Ueps);
		} //setting up problem

		// Fill in lb, ub, C, e
		Matrix<X_DIM>& xT = X[T-1];
		lb[T-1][0] = MAX(xGoal[0] - delta, xT[0] - Xeps);
		lb[T-1][1] = MAX(xGoal[1] - delta, xT[1] - Xeps);

		ub[T-1][0] = MIN(xGoal[0] + delta, xT[0] + Xeps);
		ub[T-1][1] = MIN(xGoal[1] + delta, xT[1] + Xeps);

		// Verify problem inputs
		/*
			if (!isValidInputs(result)) {
				std::cout << "Inputs are not valid!" << std::endl;
				exit(0);
			}
		 */

		int exitflag = lpMPC_solve(&problem, &output, &info);
		if (exitflag == 1) {
			for(int t = 0; t < T-1; ++t) {
				Matrix<X_DIM>& xt = Xopt[t];
				Matrix<U_DIM>& ut = Uopt[t];

				for(int i = 0; i < X_DIM; ++i) {
					xt[i] = z[t][i];
				}
				for(int i = 0; i < U_DIM; ++i) {
					ut[i] = z[t][X_DIM+i];
				}
				optcost = info.pobj;
			}
			Matrix<X_DIM>& xt = Xopt[T-1];
			xt[0] = z[T-1][0]; xt[1] = z[T-1][1];
		}
		else {
			LOG_FATAL("Some problem in solver");
			std::exit(-1);
		}

		model_merit = optcost + constant_cost; // need to add constant terms that were dropped

		initVarVals(Xopt, Uopt);
		evalCost(resultCost, vars);
		new_merit = resultCost[0];

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
			delete resultCost;
			delete[] resultCostGrad;
			return INFTY;
		} else if (approx_merit_improve < cfg::min_approx_improve) {
			LOG_DEBUG("Converged: improvement small enough");
			X = Xopt; U = Uopt;
			solution_accepted = true;
			break;
		} else if ((exact_merit_improve < 0) || (merit_improve_ratio < cfg::improve_ratio_threshold)) {
			LOG_DEBUG("Shrinking trust region size to: %2.6f %2.6f", Xeps, Ueps);
			Xeps *= cfg::trust_shrink_ratio;
			Ueps *= cfg::trust_shrink_ratio;
			solution_accepted = false;
		} else {
			LOG_DEBUG("Accepted, Increasing trust region size to:  %2.6f %2.6f", Xeps, Ueps);
			// expand Xeps and Ueps and break into outermost loop (which we don't have)
			Xeps *= cfg::trust_expand_ratio;
			Ueps *= cfg::trust_expand_ratio;
			X = Xopt; U = Uopt;
			prevcost = optcost;
			solution_accepted = true;
		}
	}

	initVarVals(X, U);
	evalCost(resultCost, vars);
	optcost = resultCost[0];

	delete resultCost;
	delete[] resultCostGrad;

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

	int nparams = 3;
	int nvars = T * X_DIM + (T - 1) * U_DIM + Q_DIM + R_DIM + (X_DIM * X_DIM) + nparams;
	std::stringstream fileName;
	fileName << "point/masks/state-masks-" << TIMESTEPS << ".txt";
	setupDstarInterface(fileName.str(), nparams, nvars);

	lpMPC_params problem;
	lpMPC_output output;
	lpMPC_info info;

	setupLpMPCVars(problem, output);

	util::Timer solveTimer;
	util::Timer_tic(&solveTimer);

	double cost = lpCollocation(X, U, problem, output, info);

	double solvetime = util::Timer_toc(&solveTimer);
	//LOG_INFO("Cost: %4.10f", cost);
	LOG_INFO("Solve time: %5.3f ms", solvetime*1000);

	Matrix<B_DIM> binit = zeros<B_DIM>();
	std::vector<Matrix<B_DIM> > B(T, binit);

	vec(X[0], SqrtSigma0, B[0]);
	for (size_t t = 0; t < T-1; ++t) {
		B[t+1] = beliefDynamics(B[t], U[t]);
	}

	pythonDisplayTrajectory(B, U);

	cleanup();

	/*
	std::cout << "Final X:" << std::endl;
	for (int t = 0; t < T; ++t) {
		std::cout << ~X[t] << std::endl;
	}
	 */

	int k;
	std::cin >> k;

	//CAL_End();
	return 0;
}
