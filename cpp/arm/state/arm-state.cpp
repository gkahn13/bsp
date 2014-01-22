#include <vector>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <iomanip>

#include "../arm.h"

#include "util/matrix.h"
#include "util/utils.h"
#include "util/Timer.h"
#include "util/logging.h"
//#include "util/adolcsupport.h"

#include <Python.h>
//#include <pythonrun.h>
#include <boost/python.hpp>
#include <boost/preprocessor/iteration/local.hpp>
#include <boost/preprocessor/arithmetic/sub.hpp>
#include <boost/filesystem.hpp>

namespace py = boost::python;

//#include <Eigen/Dense>
//#include <Eigen/Core>
//#include </usr/include/eigen3/unsupported/Eigen/AdolcForward>
//using namespace Eigen;

#include <boost/numeric/ublas/matrix.hpp>
using namespace boost::numeric::ublas;

#include <adolc/adouble.h>

extern "C" {
//#include "../sym/state-symeval.h"
#include "armStateMPC.h"
}

namespace cfg {
const double improve_ratio_threshold = .1;
const double min_approx_improve = 1e-2;
const double min_trust_box_size = 1e-2;
const double trust_shrink_ratio = .1;
const double trust_expand_ratio = 1.5;
}

// armStateMPC vars
armStateMPC_FLOAT **H, **f, **lb, **ub, **C, **e, **z;


void setupStateMPCVars(armStateMPC_params& problem, armStateMPC_output& output)
{
	// inputs
	H = new armStateMPC_FLOAT*[T];
	f = new armStateMPC_FLOAT*[T];
	lb = new armStateMPC_FLOAT*[T];
	ub = new armStateMPC_FLOAT*[T];
	C = new armStateMPC_FLOAT*[T-1];
	e = new armStateMPC_FLOAT*[T];

	// output
	z = new armStateMPC_FLOAT*[T];

#define SET_VARS(n)    \
		H[ BOOST_PP_SUB(n,1) ] = problem.H##n ;  \
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
		H[ BOOST_PP_SUB(n,1) ] = problem.H##n ;  \
		f[ BOOST_PP_SUB(n,1) ] = problem.f##n ;  \
		lb[ BOOST_PP_SUB(n,1) ] = problem.lb##n ;	\
		ub[ BOOST_PP_SUB(n,1) ] = problem.ub##n ;	\
		e[ BOOST_PP_SUB(n,1) ] = problem.e##n ;  \
		z[ BOOST_PP_SUB(n,1) ] = output.z##n ;

#define BOOST_PP_LOCAL_MACRO(n) SET_LAST_VARS(n)
#define BOOST_PP_LOCAL_LIMITS (TIMESTEPS, TIMESTEPS)
#include BOOST_PP_LOCAL_ITERATE()

}

/*
void cleanup()
{
	delete[] inputVars;
	delete[] vars;

	delete[] H;
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
	//vars = new double[nvars];

	// For evaluation
	idx = 0;
	for (int i = 0; i < nvars; ++i) {
		vars[idx++] = inputVars[maskIndices[i]];
	}
}

inline void forcePsdHessian() {
	// zero out negative diagonal entries
	if (force_type == 0) {
		for(int t = 0; t < T-1; ++t) {
			for(int i = 0; i < (X_DIM+U_DIM); ++i) {
				H[t][i] = (H[t][i] < 0) ? 0 : H[t][i];
			}
		}

		for(int i = 0; i < X_DIM; ++i) {
			H[T-1][i] = (H[T-1][i] < 0) ? 0 : H[T-1][i];
		}
	}
}

bool isValidInputs(double *result) {

	//armStateMPC_FLOAT **H, **f, **lb, **ub, **C, **e, **z;
	for(int t = 0; t < T-1; ++t) {
		std::cout << "t: " << t << std::endl << std::endl;



		std::cout << "H: " << std::endl;
		for(int i = 0; i < 4; ++i) {
			std::cout << H[t][i] << std::endl;
		}

		std::cout << std::endl << std::endl;
	}
	return true;
}


double stateCollocation(std::vector< Matrix<X_DIM> >& X, std::vector< Matrix<U_DIM> >& U, armStateMPC_params& problem, armStateMPC_output& output, armStateMPC_info& info)
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
	double constant_cost, hessian_constant, jac_constant;

	int dim = T*X_DIM + (T-1)*U_DIM;
	double* resultCost = new double;
	double* resultCostGradDiagHess = new double[2*dim + 1];

	double Hzbar[X_DIM+U_DIM];

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

		if (solution_accepted) {
			initVarVals(X, U);
			evalCostGradDiagHess(resultCostGradDiagHess, vars);
			merit = resultCostGradDiagHess[0];

			// evaluate constant cost term (omitted from optimization)
			// Need to compute:
			// ~z_bar * H * z_bar (easy to compute since H is diagonal)
			// -f * z_bar
			constant_cost = 0;
			hessian_constant = 0;
			jac_constant = 0;

			// compute Hessian first
			// so can force it to be PSD
			for (int t = 0; t < T-1; ++t) {
				H[t][0] = resultCostGradDiagHess[1+dim+t*X_DIM];
				H[t][1] = resultCostGradDiagHess[1+dim+t*X_DIM+1];
				H[t][2] = resultCostGradDiagHess[1+dim+T*X_DIM+t*U_DIM];
				H[t][3] = resultCostGradDiagHess[1+dim+T*X_DIM+t*U_DIM+1];
			}
			H[T-1][0] = resultCostGradDiagHess[1+dim+(T-1)*X_DIM];
			H[T-1][1] = resultCostGradDiagHess[1+dim+(T-1)*X_DIM+1];

			forcePsdHessian();

			for (int t = 0; t < T-1; ++t)
			{
				Matrix<X_DIM>& xt = X[t];
				Matrix<U_DIM>& ut = U[t];

				Matrix<X_DIM+U_DIM> zbar;
				zbar.insert(0,0,xt);
				zbar.insert(X_DIM,0,ut);

				for(int i = 0; i < (X_DIM+U_DIM); ++i) {
					Hzbar[i] = H[t][i]*zbar[i];
				}

				f[t][0] = resultCostGradDiagHess[1+t*X_DIM];
				f[t][1] = resultCostGradDiagHess[1+t*X_DIM+1];
				f[t][2] = resultCostGradDiagHess[1+T*X_DIM+t*U_DIM];
				f[t][3] = resultCostGradDiagHess[1+T*X_DIM+t*U_DIM+1];

				for(int i = 0; i < (X_DIM+U_DIM); ++i) {
					hessian_constant += H[t][i]*zbar[i]*zbar[i];
					jac_constant += -(f[t][i]*zbar[i]);
				}

				for(int i = 0; i < (X_DIM+U_DIM); ++i) {
					f[t][i] -= Hzbar[i];
				}

				// TODO: move following CMat code outside, is same every iteration
				Matrix<X_DIM,X_DIM+U_DIM> CMat;

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

			f[T-1][0] = resultCostGradDiagHess[1+(T-1)*X_DIM];
			f[T-1][1] = resultCostGradDiagHess[1+(T-1)*X_DIM+1];

			for(int i = 0; i < X_DIM; ++i) {
				hessian_constant += H[T-1][i]*xT[i]*xT[i];
				jac_constant += -(f[T-1][i]*xT[i]);
			}

			for(int i = 0; i < X_DIM; ++i) {
				Hzbar[i] = H[T-1][i]*xT[i];
				f[T-1][i] -= Hzbar[i];
			}

			e[T-1][0] = 0; e[T-1][1] = 0;

			constant_cost = 0.5*hessian_constant + jac_constant + resultCostGradDiagHess[0];
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

		//if (!isValidInputs(result)) {
		//	std::cout << "Inputs are not valid!" << std::endl;
		//	exit(0);
		//}


		int exitflag = armStateMPC_solve(&problem, &output, &info);
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
			delete[] resultCostGradDiagHess;
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
			// expand Xeps and Ueps
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
	delete[] resultCostGradDiagHess;

	return optcost;
}
*/

int main(int argc, char* argv[])
{
	//matrix<adouble> M (3, 3);
	//adouble adArr;

	//M(0,0) = adArr;

	x0[0] = .5*M_PI; x0[1] = -1.5431281995798991; x0[2] = -0.047595544887998331;
	x0[3] = 1.4423058659586809; x0[4] = 1.5334368368992011; x0[5] = -1.1431255223182604;

	SqrtSigma0 = identity<X_DIM>();

	xGoal[0] = -M_PI/2; xGoal[1] = .25*M_PI; xGoal[2] = -.5*M_PI;
	xGoal[3] = -M_PI/2; xGoal[4] = -.25*M_PI; xGoal[5] = .5*M_PI;

	xMin[0] = -M_PI; xMin[1] = -0.75*M_PI; xMin[2] = -0.75*M_PI; xMin[3] = -M_PI; xMin[4] = -0.75*M_PI; xMin[5] = -M_PI;
	xMax[0] = M_PI;  xMax[1] = 0.75*M_PI;  xMax[2] = 0.75*M_PI;  xMax[3] = M_PI;  xMax[4] =  0.75*M_PI; xMax[5] = M_PI;

	for(int i=0; i < U_DIM; ++i) { uMin[i] = -.5; uMax[i] = .5; }

	Matrix<U_DIM> uinit = (xGoal - x0) / (T-1);

	std::vector<Matrix<U_DIM> > U(T-1, uinit);
	std::vector<Matrix<X_DIM> > X(T);

	X[0] = x0;
	for (int t = 0; t < T-1; ++t) {
		X[t+1] = dynfunc(X[t], U[t], zeros<Q_DIM,1>());
	}

	//setupDstarInterface(std::string(getMask()));


	armStateMPC_params problem;
	armStateMPC_output output;
	armStateMPC_info info;

	setupStateMPCVars(problem, output);

	util::Timer solveTimer;
	util::Timer_tic(&solveTimer);

	// compute cost for the trajectory

	//for(int i = 0; i < 100000; ++i) {
	//	for (int t = 0; t < T-1; ++t) {
	//		X[t+1] = dynfunc(X[t], U[t], zeros<Q_DIM,1>());
	//	}
	//}


	std::cout << "before state collocation " << std::endl;

	double cost = stateCollocation(X, U, problem, output, info);

	std::cout << "after state collocation" << std::endl;

	double solvetime = util::Timer_toc(&solveTimer);
	LOG_INFO("Cost: %4.10f", cost);
	LOG_INFO("Solve time: %5.3f ms", solvetime*1000);

	Matrix<B_DIM> binit = zeros<B_DIM>();
	std::vector<Matrix<B_DIM> > B(T, binit);

	vec(X[0], SqrtSigma0, B[0]);
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
