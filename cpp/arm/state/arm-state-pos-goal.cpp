#include <vector>
#include <iomanip>

#include "util/matrix.h"
#include "util/Timer.h"

extern "C" {
#include "arm-state-pos-goal-MPC.h"
statePenaltyMPC_FLOAT **Q, **f, **lb, **ub, **z;
statePenaltyMPC_FLOAT *A, *b, *e;
#include "arm-state-casadi.h"
}

#include "boost/preprocessor.hpp"

#include "../arm.h"

namespace cfg {
const double improve_ratio_threshold = .1;
const double min_approx_improve = 1e-4;
const double min_trust_box_size = 1e-3;
const double trust_shrink_ratio = .5;
const double trust_expand_ratio = 1.5;
const double cnt_tolerance = 1e-4;
const double penalty_coeff_increase_ratio = 5;
const double initial_penalty_coeff = 5;
const double initial_trust_box_size = 1;
const int max_penalty_coeff_increases = 3;
const int max_sqp_iterations = 50;
}

double computeLQGMPcost(const std::vector<Matrix<X_DIM> >& X, const std::vector<Matrix<U_DIM> >& U)
{
	std::vector<Matrix<G_DIM, X_DIM> > J;
	std::vector<Matrix<X_DIM, X_DIM> > Sigma;
	preprocess(X, J, Sigma);

	double cost = 0;
	for(int t = 0; t < T-1; ++t) {
		cost += alpha_belief*tr(J[t]*Sigma[t]*~J[t]) + alpha_control*tr(~U[t]*U[t]);
	}
	cost += alpha_final_belief*tr(J[T-1]*Sigma[T-1]*~J[T-1]);
	return cost;
}

double casadiComputeCost(const std::vector< Matrix<X_DIM> >& X, const std::vector< Matrix<U_DIM> >& U)
{
	double cost = 0;

	int index = 0;
	double *XU_arr = new double[X_DIM*T + U_DIM*(T-1)];
	for(int t = 0; t < T-1; ++t) {
		for(int i=0; i < X_DIM; ++i) {
			XU_arr[index++] = X[t][i];
		}

		for(int i=0; i < U_DIM; ++i) {
			XU_arr[index++] = U[t][i];
		}
	}
	for(int i=0; i < X_DIM; ++i) {
		XU_arr[index++] = X[T-1][i];
	}

	Matrix<X_DIM,X_DIM> Sigma0 = SqrtSigma0*SqrtSigma0;
	double *Sigma0_arr = new double[X_DIM*X_DIM];
	index = 0;
	for(int i=0; i < X_DIM; ++i) {
		for(int j=0; j < X_DIM; ++j) {
			Sigma0_arr[index++] = Sigma0(i,j);
		}
	}

	double *params_arr = new double[3];
	params_arr[0] = alpha_belief;
	params_arr[1] = alpha_control;
	params_arr[2] = alpha_final_belief;

	double *cam0_arr = new double[3];
	cam0_arr[0] = cam0(0,0);
	cam0_arr[1] = cam0(1,0);
	cam0_arr[2] = cam0(2,0);

	double *cam1_arr = new double[3];
	cam1_arr[0] = cam1(0,0);
	cam1_arr[1] = cam1(1,0);
	cam1_arr[2] = cam1(2,0);

	const double **casadi_input = new double*[5];
	casadi_input[0] = XU_arr;
	casadi_input[1] = Sigma0_arr;
	casadi_input[2] = params_arr;
	casadi_input[3] = cam0_arr;
	casadi_input[4] = cam1_arr;

	double **cost_arr = new double*[1];
	cost_arr[0] = &cost;

	evaluateCostWrap(casadi_input, cost_arr);

	return cost;
}



double computeCost(const std::vector< Matrix<X_DIM> >& X, const std::vector< Matrix<U_DIM> >& U)
{
	double cost = 0;
	Matrix<B_DIM> b;
	Matrix<X_DIM> x;
	Matrix<X_DIM, X_DIM> SqrtSigma;
	vec(x0, SqrtSigma0, b);
	Matrix<G_DIM,X_DIM> J;

	for(int t = 0; t < T-1; ++t) {
		unVec(b, x, SqrtSigma);
		linearizeg(x, J);

		cost += alpha_belief*tr(J*SqrtSigma*SqrtSigma*~J) + alpha_control*tr(~U[t]*U[t]);
		b = beliefDynamics(b, U[t]);
	}
	unVec(b, x, SqrtSigma);
	linearizeg(x, J);
	cost += alpha_final_belief*tr(J*SqrtSigma*SqrtSigma*~J);
	return cost;
}

double computeMerit(const std::vector< Matrix<X_DIM> >& X, const std::vector< Matrix<U_DIM> >& U, double penalty_coeff)
{
	double merit = 0;
	Matrix<B_DIM> b;
	Matrix<X_DIM> x;
	Matrix<X_DIM,X_DIM> SqrtSigma;
	vec(x0, SqrtSigma0, b);
	Matrix<G_DIM,X_DIM> J;

	for(int t = 0; t < T-1; ++t) {
		unVec(b, x, SqrtSigma);
		linearizeg(x, J);

		merit += alpha_belief*tr(J*SqrtSigma*SqrtSigma*~J) + alpha_control*tr(~U[t]*U[t]);
		b = beliefDynamics(b, U[t]);
	}
	unVec(b, x, SqrtSigma);
	linearizeg(x, J);
	merit += alpha_final_belief*tr(J*SqrtSigma*SqrtSigma*~J);

	Matrix<G_DIM> delta;
	delta[0] = delta[1] = delta[2] = goaldelta;

	Matrix<2*G_DIM> goalposviol;
	goalposviol.insert<G_DIM,1>(0,0,g(x) - posGoal - delta);
	goalposviol.insert<G_DIM,1>(G_DIM,0, -g(x) + posGoal - delta);

	for(int i = 0; i < 2*G_DIM; ++i) {
		merit += penalty_coeff*MAX(goalposviol[i],0);
	}
	
	return merit;
}

void setupBeliefVars(statePenaltyMPC_params& problem, statePenaltyMPC_output& output)
{
	// problem inputs
	Q = new statePenaltyMPC_FLOAT*[T];
	f = new statePenaltyMPC_FLOAT*[T];
	lb = new statePenaltyMPC_FLOAT*[T];
	ub = new statePenaltyMPC_FLOAT*[T];

	// problem outputs
	z = new statePenaltyMPC_FLOAT*[T];

	// initial state
	e = problem.e1;

#define SET_VARS(n)    \
		Q[ BOOST_PP_SUB(n,1) ] = problem.Q##n ;  \
		f[ BOOST_PP_SUB(n,1) ] = problem.f##n ;  \
		lb[ BOOST_PP_SUB(n,1) ] = problem.lb##n ;	\
		ub[ BOOST_PP_SUB(n,1) ] = problem.ub##n ;	\
		z[ BOOST_PP_SUB(n,1) ] = output.z##n ;

#define BOOST_PP_LOCAL_MACRO(n) SET_VARS(n)
#define BOOST_PP_LOCAL_LIMITS (1, TIMESTEPS-1)
#include BOOST_PP_LOCAL_ITERATE()

#define SET_LAST_VARS(n)    \
		Q[ BOOST_PP_SUB(n,1) ] = problem.Q##n ;  \
		f[ BOOST_PP_SUB(n,1) ] = problem.f##n ;  \
		lb[ BOOST_PP_SUB(n,1) ] = problem.lb##n ;	\
		ub[ BOOST_PP_SUB(n,1) ] = problem.ub##n ;	\
		A = problem.A##n; \
		b = problem.b##n; \
		z[ BOOST_PP_SUB(n,1) ] = output.z##n ;

#define BOOST_PP_LOCAL_MACRO(n) SET_LAST_VARS(n)
#define BOOST_PP_LOCAL_LIMITS (TIMESTEPS, TIMESTEPS)
#include BOOST_PP_LOCAL_ITERATE()

}

void cleanupBeliefMPCVars()
{
	delete[] Q;
	delete[] f;
	delete[] lb;
	delete[] ub;
	delete[] z;
}

// TODO: Check if all inputs are valid, Q, f, lb, ub, A, b at last time step
bool isValidInputs()
{
	for(int t = 0; t < T-1; ++t) {

		// check if Q, f, lb, ub, C, e, b are valid!

		//std::cout << std::endl << std::endl;
	}

	for(int i = 0; i < 12; ++i) {
		std::cout << lb[T-1][i] << " ";
	}
	std::cout << "\n\n";

	for(int i = 0; i < 6; ++i) {
		std::cout << ub[T-1][i] << " ";
	}

	//for(int i = 0; i < 198; ++i) {
	//	std::cout << A[i] << " ";
	//}

	//for(int i = 0; i < 6; ++i) {
	//	std::cout << b[i] << "  ";
	//}
	std::cout << std::endl;

	return true;
}

/*
bool minimizeMeritFunction(std::vector< Matrix<X_DIM> >& X, std::vector< Matrix<U_DIM> >& U, statePenaltyMPC_params& problem, statePenaltyMPC_output& output, statePenaltyMPC_info& info, double penalty_coeff, double trust_box_size)
{
	//LOG_DEBUG("Solving sqp problem with penalty parameter: %2.4f", penalty_coeff);
	std::cout << "Solving sqp problem with penalty parameter: " << penalty_coeff << std::endl;

	Matrix<X_DIM,1> x0 = X[0];

	double Xeps = trust_box_size;
	double Ueps = trust_box_size;

	double prevcost, optcost;

	std::vector<Matrix<X_DIM> > Xopt(T);
	std::vector<Matrix<U_DIM> > Uopt(T-1);

	double merit, model_merit, new_merit;
	double approx_merit_improve, exact_merit_improve, merit_improve_ratio;

	int sqp_iter = 1, index = 0;
	bool success;

	Matrix<X_DIM+U_DIM, X_DIM+U_DIM> QMat;
	Matrix<X_DIM+2*G_DIM,X_DIM+2*G_DIM> QfMat;
	Matrix<X_DIM> eVec;
	Matrix<S_DIM,S_DIM> Hess;
	Matrix<2*G_DIM,X_DIM+2*G_DIM> AMat;
	Matrix<2*G_DIM,1> bVec;

	// sqp loop
	while(true)
	{
		// In this loop, we repeatedly construct a linear approximation to the nonlinear belief dynamics constraint
		//LOG_DEBUG("  sqp iter: %d", sqp_iter);
		std::cout << "  sqp iter: " << sqp_iter << std::endl;

		merit = computeMerit(X, U, penalty_coeff);
		
		//LOG_DEBUG("  merit: %4.10f", merit);
		std::cout << "  merit: " << merit << std::endl;

		// Problem linearization and definition
		// fill in Q, f, C, e
		
		for (int t = 0; t < T-1; ++t) 
		{
			Matrix<B_DIM>& xt = X[t];
			Matrix<U_DIM>& ut = U[t];

			constructHessian(bt, Hess);

			QMat.reset();
			QMat.insert<S_DIM,S_DIM>(X_DIM,X_DIM,2*alpha_belief*Hess);
			QMat.insert<U_DIM,U_DIM>(B_DIM,B_DIM,2*alpha_control*identity<U_DIM>());
			
			fillColMajor(Q[t], QMat);

			for(int i = 0; i < (B_DIM+U_DIM); ++i) {
				f[t][i] = 0;
			}
			for(int i = 0; i < 2*B_DIM; ++i) {
				f[t][B_DIM+U_DIM+i] = penalty_coeff;
			}
		}
		
		// For last stage, fill in Q, f, D, e
		constructHessian(B[T-1], Hess);

		QfMat.reset();
		QfMat.insert<S_DIM,S_DIM>(X_DIM,X_DIM,2*alpha_final_belief*Hess);
		
		fillColMajor(Q[T-1], QfMat);
				
		for(int i = 0; i < B_DIM; ++i) {
			f[T-1][i] = 0;
		}
		for(int i = 0; i < 2*G_DIM; ++i) {
			f[T-1][B_DIM+i] = penalty_coeff;
		}

		// fill in A and b
		Matrix<X_DIM> xT = X[T-1];
		
		Matrix<G_DIM,X_DIM> J;
		linearizeg(xT, J);

		AMat.reset();
		AMat.insert<G_DIM,X_DIM>(0,0,J);
		AMat.insert<G_DIM,X_DIM>(G_DIM,0,-J);

		fillColMajor(A, AMat);

		Matrix<G_DIM> delta;
		delta [0] = delta[1] = delta[2] = goaldelta;
		
		bVec.insert<G_DIM,1>(0,0,posGoal - g(xT) + J*xT + delta);
		bVec.insert<G_DIM,1>(G_DIM,0,-posGoal + g(xT) - J*xT + delta);

		fillColMajor(b, bVec);
		
		//std::cout << "PAUSED INSIDE MINIMIZEMERITFUNCTION" << std::endl;
		//int k;
		//std::cin >> k;


		// trust region size adjustment
		while(true)
		{
			//LOG_DEBUG("       trust region size: %2.6f %2.6f", Xeps, Ueps);
			std::cout << "       trust region size: " << Xeps << ", " << Ueps << std::endl;

			// solve the innermost QP here
			for(int t = 0; t < T-1; ++t)
			{
				Matrix<X_DIM>& xt = X[t];
				Matrix<U_DIM>& ut = U[t];

				// Fill in lb, ub

				index = 0;
				// x lower bound
				for(int i = 0; i < X_DIM; ++i) { lb[t][index++] = MAX(xMin[i], xt[i] - Xeps); }
				// u lower bound
				for(int i = 0; i < U_DIM; ++i) { lb[t][index++] = MAX(uMin[i], ut[i] - Ueps); }

				index = 0;
				// x upper bound
				for(int i = 0; i < X_DIM; ++i) { ub[t][index++] = MIN(xMax[i], xt[i] + Xeps); }
				// u upper bound
				for(int i = 0; i < U_DIM; ++i) { ub[t][index++] = MIN(uMax[i], ut[i] + Ueps); }
			}

			Matrix<X_DIM>& xT = X[T-1];

			// Fill in lb, ub, C, e
			index = 0;
			// xGoal lower bound
			for(int i = 0; i < X_DIM; ++i) { lb[T-1][index++] = MAX(xMin[i], xT[i] - Xeps); }
			
			// for lower bound on L1 slacks
			for(int i = 0; i < 2*G_DIM; ++i) { lb[T-1][index++] = 0; }

			index = 0;
			// xGoal upper bound
			for(int i = 0; i < X_DIM; ++i) { ub[T-1][index++] = MIN(xMax[i], xT[i] + Xeps); }
			
			// Verify problem inputs
			//if (!isValidInputs()) {
			//	std::cout << "Inputs are not valid!" << std::endl;
			//	exit(-1);
			//}

			//std::cerr << "PAUSING INSIDE MINIMIZE MERIT FUNCTION FOR INPUT VERIFICATION" << std::endl;
			//int num;
			//std::cin >> num;

			int exitflag = statePenaltyMPC_solve(&problem, &output, &info);
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
				for(int i = 0; i < X_DIM; ++i) {
					Xopt[T-1][i] = z[T-1][i];
				}
			}
			else {
				//LOG_ERROR("Some problem in solver");
				std::cerr << "Some problem in solver" << std::endl;
				std::exit(-1);
			}

			//LOG_DEBUG("Optimized cost: %4.10f", optcost);
			std::cout << "       Optimized cost: " << optcost << std::endl;

			model_merit = optcost;
			new_merit = computeMerit(Bopt, Uopt, penalty_coeff);

			//LOG_DEBUG("merit: %4.10f", merit);
			//LOG_DEBUG("model_merit: %4.10f", model_merit);
			//LOG_DEBUG("new_merit: %4.10f", new_merit);
			//std::cout << "       merit: " << merit << std::endl;
			//std::cout << "       model_merit: " << model_merit << std::endl;
			//std::cout << "       new_merit: " << new_merit << std::endl;
			
			approx_merit_improve = merit - model_merit;
			exact_merit_improve = merit - new_merit;
			merit_improve_ratio = exact_merit_improve / approx_merit_improve;

			//LOG_DEBUG("approx_merit_improve: %1.6f", approx_merit_improve);
			//LOG_DEBUG("exact_merit_improve: %1.6f", exact_merit_improve);
			//LOG_DEBUG("merit_improve_ratio: %1.6f", merit_improve_ratio);
			//std::cout << "       approx_merit_improve: " << approx_merit_improve << std::endl;
			//std::cout << "       exact_merit_improve: " << exact_merit_improve << std::endl;
			std::cout << "       merit_improve_ratio: " << merit_improve_ratio << std::endl;
			
			//std::cout << "PAUSED INSIDE minimizeMeritFunction" << std::endl;
			//int num;
			//std::cin >> num;

			if (approx_merit_improve < -1e-5) {
				//LOG_ERROR("Approximate merit function got worse: %1.6f", approx_merit_improve);
				//LOG_ERROR("Either convexification is wrong to zeroth order, or you are in numerical trouble");
				//LOG_ERROR("Failure!");
				std::cerr << "Approximate merit function got worse: " << approx_merit_improve << std::endl;
				std::cerr << "Either convexification is wrong to zeroth order, or you are in numerical trouble" << std::endl;
				std::cerr << "Failure!" << std::endl;

				success = false;
			} else if (approx_merit_improve < cfg::min_approx_improve) {
				//LOG_DEBUG("Converged: improvement small enough");
				std::cout << "  Converged: improvement small enough" << std::endl;
				B = Bopt; U = Uopt;
				return true;
			} else if ((exact_merit_improve < 0) || (merit_improve_ratio < cfg::improve_ratio_threshold)) {
				Xeps *= cfg::trust_shrink_ratio;
				Ueps *= cfg::trust_shrink_ratio;
				//LOG_DEBUG("Shrinking trust region size to: %2.6f %2.6f", Xeps, Ueps);
				std::cout << "  Shrinking trust region size to: " << Xeps << ", " << Ueps << std::endl;
			} else {
				Xeps *= cfg::trust_expand_ratio;
				Ueps *= cfg::trust_expand_ratio;
				B = Bopt; U = Uopt;
				prevcost = optcost;
				//LOG_DEBUG("Accepted, Increasing trust region size to:  %2.6f %2.6f", Xeps, Ueps);
				std::cout << "  Accepted, Increasing trust region size to:  " << Xeps << ", " << Ueps << std::endl;
				break;
			}

			if (Xeps < cfg::min_trust_box_size && Ueps < cfg::min_trust_box_size) {
			    //LOG_DEBUG("Converged: x tolerance");
				std::cout << "  Converged: x tolerance" << std::endl;
			    return true;
			}

		} // trust region loop
		sqp_iter++;
	} // sqp loop

	return success;
}


double statePenaltyCollocation(std::vector< Matrix<X_DIM> >& X, std::vector< Matrix<U_DIM> >& U, statePenaltyMPC_params& problem, statePenaltyMPC_output& output, statePenaltyMPC_info& info)
{
	double costTime = 0;

	double penalty_coeff = cfg::initial_penalty_coeff;
	double trust_box_size = cfg::initial_trust_box_size;

	int penalty_increases = 0;

	Matrix<2*G_DIM> goalposviol;

	// penalty loop
	while(penalty_increases < cfg::max_penalty_coeff_increases)
	{
		bool success = minimizeMeritFunction(B, U, problem, output, info, penalty_coeff, trust_box_size);

		double cntviol = 0;
	
		Matrix<G_DIM> delta;
		delta[0] = delta[1] = delta[2] = goaldelta;

		goalposviol.insert<G_DIM,1>(0,0,g(X[T-1]) - posGoal - delta);
		goalposviol.insert<G_DIM,1>(G_DIM,0, -g(X[T-1]) + posGoal - delta);

		for(int i = 0; i < 2*G_DIM; ++i) {
			cntviol += MAX(goalposviol[i],0);
		}
	
	    success = success && (cntviol < cfg::cnt_tolerance);
	    
		//LOG_DEBUG("Constraint violations: %2.10f",cntviol);
		std::cout << "Constraint violations: " << cntviol << std::endl;

	    if (!success) {
	        penalty_increases++;
	        penalty_coeff = penalty_coeff*cfg::penalty_coeff_increase_ratio;
	        trust_box_size = cfg::initial_trust_box_size;
	    }
	    else {
	    	return computeCost(X, U);
	    }
	}
	return computeCost(X, U);
}
*/

bool testInitializationFeasibility(const std::vector<Matrix<X_DIM> >& X, const std::vector<Matrix<U_DIM> >& U)
{
	LOG_DEBUG("X initial");
	for (int t = 0; t < T; ++t) { 
		const Matrix<X_DIM>& xt = X[t];
		for (int i = 0; i < X_DIM; ++i) {
			if (xt[i] > xMax[i] || xt[i] < xMin[i]) {
				LOG_ERROR("Joint angle limit violated at joint %d and time %d", i,t);
				return false;
			}
		}

		//std::cout << std::setprecision(8) << ~xt;
	}

	LOG_DEBUG("U initial");
	for (int t = 0; t < T-1; ++t) { 
		const Matrix<U_DIM>& ut = U[t];
		for (int i = 0; i < U_DIM; ++i) {
			if (ut[i] > uMax[i] || ut[i] < uMin[i]) {
				LOG_ERROR("Control limit violated at joint %d and time %d", i,t);
				return false;
			}
		}
		//std::cout << std::setprecision(8) << ~ut;
	}

	return true;
}

int main(int argc, char* argv[])
{
	
	initProblemParams();

	LOG_INFO("init problem params");

	Matrix<U_DIM> uinit = (xGoal - x0) / (double)((T-1)*DT);
	std::vector<Matrix<U_DIM> > U(T-1, uinit);

	std::vector<Matrix<X_DIM> > X(T);

	for (size_t t = 0; t < T-1; ++t) {
		X[t+1] = dynfunc(X[t], U[t], zeros<Q_DIM,1>());
		//std::cout << ~X[t] << std::endl;
	}

	bool feasible = testInitializationFeasibility(X, U);
	if (!feasible) {
		LOG_ERROR("Infeasible trajectory initialization detected");
		exit(-1);
	}



	double initTrajCost = computeCost(X, U);
	LOG_INFO("Initial trajectory cost: %4.10f", initTrajCost);

	double casadiInitCost = casadiComputeCost(X, U);
	LOG_INFO("Initial casadi cost: %4.10f", casadiInitCost);

	//readTrajectoryFromFile("data\\trajectory.txt", X);

	//double initLQGMPcost = computeLQGMPcost(X, U);
	//LOG_DEBUG("Initial trajectory LQG-MP cost: %4.10f", initLQGMPcost);

	//displayStateTrajectory(X, U, false);

	/*
	statePenaltyMPC_params problem;
	statePenaltyMPC_output output;
	statePenaltyMPC_info info;

	setupBeliefVars(problem, output);

	util::Timer solveTimer;
	Timer_tic(&solveTimer);

	double cost = statePenaltyCollocation(X, U, problem, output, info);

	double solvetime = util::Timer_toc(&solveTimer);

	LOG_INFO("Optimized cost: %4.10f", cost);
	LOG_INFO("Actual cost: %4.10f", computeCost(X, U));
	LOG_INFO("Solve time: %5.3f ms", solvetime*1000);
	
	std::cout << "Actual cost: " << computeCost(X, U) << std::endl;
	
	cleanupBeliefMPCVars();
	*/

	//double finalLQGMPcost = computeLQGMPcost(X, U);
	//LOG_DEBUG("Final trajectory LQG-MP cost: %4.10f",finalLQGMPcost);

	//saveOptimizedTrajectory(U);
	//readOptimizedTrajectory(U);
	
	/*
	vec(x0, SqrtSigma0, B[0]);
	for (size_t t = 0; t < T-1; ++t) {
		B[t+1] = beliefDynamics(B[t], U[t]);
	}
	*/
	
	LOG_INFO("Finished");
	int k;
	std::cin >> k;

	return 0;
}
