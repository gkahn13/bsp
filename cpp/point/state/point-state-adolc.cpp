#include <vector>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <iomanip>

#include "../point-adolc.h"

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


#include <boost/numeric/ublas/matrix.hpp>
using namespace boost::numeric::ublas;

#include <adolc/adolc.h>

extern "C" {
#include "../sym/state-symeval.h"
#include "stateMPC.h"
}

namespace cfg {
const double improve_ratio_threshold = .1;
const double min_approx_improve = 1e-2;
const double min_trust_box_size = 1e-2;
const double trust_shrink_ratio = .1;
const double trust_expand_ratio = 1.5;
}

// stateMPC vars
stateMPC_FLOAT **H, **f, **lb, **ub, **C, **e, **z;


void setupStateMPCVars(stateMPC_params& problem, stateMPC_output& output)
{
	// inputs
	H = new stateMPC_FLOAT*[T];
	f = new stateMPC_FLOAT*[T];
	lb = new stateMPC_FLOAT*[T];
	ub = new stateMPC_FLOAT*[T];
	C = new stateMPC_FLOAT*[T-1];
	e = new stateMPC_FLOAT*[T];

	// output
	z = new stateMPC_FLOAT*[T];

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

void initVarVals(const matrix<double>& X, const matrix<double>& U)
{
	int idx = 0;
	for (int t = 0; t < T; ++t) {
		for (int i = 0; i < X_DIM; ++i) {
			inputVars[idx++] = X(i,t);
		}
	}
	for (int t = 0; t < (T - 1); ++t) {
		for (int i = 0; i < U_DIM; ++i) {
			inputVars[idx++] = U(i,t);
		}
	}
	for (int i = 0; i < (Q_DIM+R_DIM); ++i) {
		inputVars[idx++] = 0;
	}
	for (int i = 0; i < X_DIM; ++i) {
		for (int j = 0; j < X_DIM; ++j)
		inputVars[idx++] = SqrtSigma0(i,j);
	}
	inputVars[idx++] = alpha_belief; inputVars[idx++] = alpha_control; inputVars[idx++] = alpha_final_belief;

	int nvars = (int)maskIndices.size();
	//vars = new double[nvars];

	std::cout << "mask" << std::endl;
	// For evaluation
	idx = 0;
	for (int i = 0; i < nvars; ++i) {
		vars[idx++] = inputVars[maskIndices[i]];
	}

}

inline void forcePsdHessian() {
	// zero out negative diagonal entries
	for(int t = 0; t < T-1; ++t) {
		for(int i = 0; i < (X_DIM+U_DIM); ++i) {
			H[t][i] = (H[t][i] < 0) ? 0 : H[t][i];
		}
	}

	for(int i = 0; i < X_DIM; ++i) {
		H[T-1][i] = (H[T-1][i] < 0) ? 0 : H[T-1][i];
	}
}

bool isValidInputs(double *result) {

	//stateMPC_FLOAT **H, **f, **lb, **ub, **C, **e, **z;
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

void initXUAdolcMatrixes(matrix<double>& X, matrix<double>& U, matrix<adouble>& XAdolc, matrix<adouble>& UAdolc) {
	for(int t=0; t < T-1; ++t) {
		for(int i=0; i < X_DIM; ++i) {
			XAdolc(i,t) <<= X(i,t);
		}
		for(int i=0; i < U_DIM; ++i) {
			UAdolc(i,t) <<= U(i,t);
		}
	}
	for(int i=0; i < X_DIM; ++i) {
		XAdolc(i,T-1) <<= X(i,T-1);
	}
}

void initXUVectorArray(matrix<double>& X, matrix<double>& U, double* XU) {
	int idx = 0;
	for(int t=0; t < T-1; ++t) {
		for(int i=0; i < X_DIM; ++i) {
			XU[idx++] = X(i,t);
		}
		for(int i=0; i < U_DIM; ++i) {
			XU[idx++] = U(i,t);
		}
	}
	for(int i=0; i < X_DIM; ++i) {
		XU[idx++] = X(i,T-1);
	}
}

double stateCollocation(matrix<double> X, matrix<double> U, const stateMPC_params& problem, const stateMPC_output& output, const stateMPC_info& info)
{
	int maxIter = 100;
	double Xeps = 1;
	double Ueps = 1;

	// box constraint around goal
	double delta = 0.01;
	matrix<double> x0(X_DIM,1);
	column(x0,0) = column(X,0);

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


	//evalCost(resultCost, vars); // TODO: replace with cost function
	//prevcost = resultCost[0];
	matrix<double> q(Q_DIM,1), r(R_DIM,1);
	prevcost = costfunc<double>(X, U, q, r, SqrtSigma0);

	std::cout << "prevcost: " << prevcost << std::endl;


	LOG_DEBUG("Initialization trajectory cost: %4.10f", prevcost);

	matrix<adouble> Xopt(X_DIM, T);
	matrix<adouble> Uopt(U_DIM, T-1);

	/*
	double **aXU = new double*[T];
	for(int t=0; t < T-1; ++ t)
	{
		aXU[t] = new double[X_DIM+U_DIM];
		for(int x_idx=0; x_idx < X_DIM; ++x_idx) { aXU[t][x_idx] = X(t,x_idx); }
		for(int u_idx=0; u_idx < U_DIM; ++u_idx) { aXU[t][X_DIM+u_idx] = U(t,u_idx); }
	}
	*/


	bool solution_accepted = true;

	for(int it = 0; it < maxIter; ++it)
	{
		LOG_DEBUG("\nIter: %d", it);

		if (solution_accepted) {
			initVarVals(X, U);
			evalCostGradDiagHess(resultCostGradDiagHess, vars);
			merit = resultCostGradDiagHess[0];

			matrix<adouble> XAdolc(X_DIM,T);
			matrix<adouble> UAdolc(U_DIM,T-1);
			double merit;
			adouble meritAdolc;

			trace_on(GRAD_COST_TAG);
			initXUAdolcMatrixes(X,U,XAdolc,UAdolc);
			meritAdolc = costfunc<adouble>(XAdolc, UAdolc, q, r, SqrtSigma0);
			meritAdolc >>= merit;
			trace_off(GRAD_COST_TAG);

			double* grad_arr = new double[X_DIM*T + U_DIM*(T-1)];
			double* XU_arr = new double[X_DIM*T + U_DIM*(T-1)];

			initXUVectorArray(X,U,XU_arr);

			gradient(GRAD_COST_TAG, X_DIM*T + U_DIM*(T-1), XU_arr, grad_arr);

			std::cout << "adolc gradient" << std::endl;
			int idx = 0;
			for(int t=0; t < T-1; ++t) {
				for(int i=0; i < X_DIM; ++i) {
					std::cout << grad_arr[idx++] << " ";
				}
				for(int i=0; i < U_DIM; ++i) {
					std::cout << grad_arr[idx++] << " ";
				}
				std::cout << std::endl;
			}
			for(int i=0; i < X_DIM; ++i) {
				std::cout << grad_arr[idx++] << " ";
			}
			std::cout << std::endl;

			for(int t=0; t < T-1; ++t) {
				f[t][0] = resultCostGradDiagHess[1+t*X_DIM];
				f[t][1] = resultCostGradDiagHess[1+t*X_DIM+1];
				f[t][2] = resultCostGradDiagHess[1+T*X_DIM+t*U_DIM];
				f[t][3] = resultCostGradDiagHess[1+T*X_DIM+t*U_DIM+1];
			}
			f[T-1][0] = resultCostGradDiagHess[1+(T-1)*X_DIM];
			f[T-1][1] = resultCostGradDiagHess[1+(T-1)*X_DIM+1];

			std::cout << "symeval gradient" << std::endl;
			for(int t=0; t < T-1; ++t) {
				for(int i=0; i < X_DIM+U_DIM; ++i) {
					std::cout << f[t][i] << " ";
				}
				std::cout << std::endl;
			}
			for(int i=0; i < X_DIM; ++i) {
				std::cout << f[T-1][i] << " ";
			}
			std::cout << std::endl;

			exit(0);

			/*


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
				*/

			} //setting up problem

			/*
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
			*/
		} // end solution_accepted


		/*
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


		int exitflag = stateMPC_solve(&problem, &output, &info);
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
	*/
	return 0;
}

int main(int argc, char* argv[])
{
	x0(0,0) = -3.5; x0(1,0) = 2;
	SqrtSigma0 = identity_matrix<double>(X_DIM);
	xGoal(0,0) = -3.5; xGoal(1,0) = -2;

	xMin(0,0) = -5; xMin(1,0) = -3;
	xMax(0,0) = 5;  xMax(1,0) = 3;

	for(int i=0; i < U_DIM; ++i) { uMin(i,0) = -1; uMax(i,0) = 1; }

	matrix<double> uinit = (xGoal - x0) / (T-1);

	matrix<double> U(U_DIM,T-1);
	for(int t=0; t < T-1; ++t) { column(U,t) = column(uinit,0); }

	printMatrix("initial U", U);


	matrix<double> X = zeroMatrix<double>(X_DIM,T);

	printMatrix("pre initial X", X);

	column(X,0) = column(x0,0);
	matrix<double> Xcol(X_DIM,1), Ucol(U_DIM,1);
	for (int t = 0; t < T-1; ++t) {
		column(Xcol,0) = column(X,t);
		column(Ucol,0) = column(U,t);
		matrix<double>q = zeroMatrix<double>(Q_DIM,1);

		matrix<double> x_tp1 = dynfunc<double>(Xcol, Ucol, q);
		column(X,t+1) = column(x_tp1,0);
	}

	printMatrix("initial X", X);
	setupDstarInterface(std::string(getMask()));

	stateMPC_params problem;
	stateMPC_output output;
	stateMPC_info info;

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

	/*

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
	*/
	//CAL_End();
	return 0;


}
