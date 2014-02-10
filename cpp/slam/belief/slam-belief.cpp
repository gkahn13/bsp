#include <vector>
#include <iomanip>
#include <sys/resource.h>

#include "util/matrix.h"
#include "util/Timer.h"
#include "util/logging.h"

extern "C" {
#include "beliefPenaltyMPC.h"
beliefPenaltyMPC_FLOAT **H, **f, **lb, **ub, **C, **D, **e, **z;
}

#include "boost/preprocessor.hpp"

#include "../slam.h"
#include "../traj/slam-traj.h"

const double alpha_belief = 10, alpha_final_belief = 50, alpha_control = .01;

namespace cfg {
const double improve_ratio_threshold = .1;
const double min_approx_improve = 1e-4;
const double min_trust_box_size = 1e-3;
const double trust_shrink_ratio = .5;
const double trust_expand_ratio = 1.5;
const double cnt_tolerance = 1e-4;
const double penalty_coeff_increase_ratio = 5;
const double initial_penalty_coeff = 5;
const double initial_trust_box_size = 10;
const int max_penalty_coeff_increases = 3;
const int max_sqp_iterations = 50;
}

// utility to fill Matrix in column major format in FORCES array
template <size_t _numRows>
inline void fillCol(double *X, const Matrix<_numRows>& XCol) {
	int idx = 0;
	for(size_t r = 0; r < _numRows; ++r) {
		X[idx++] = XCol[r];
	}
}

template <size_t _numRows, size_t _numColumns>
inline void fillColMajor(double *X, const Matrix<_numRows, _numColumns>& XMat) {
	int idx = 0;
	for(size_t c = 0; c < _numColumns; ++c) {
		for(size_t r = 0; r < _numRows; ++r) {
			X[idx++] = XMat[c + r*_numColumns];
		}
	}
}



double computeCost(const std::vector< Matrix<B_DIM> >& B, const std::vector< Matrix<U_DIM> >& U)
{
	double cost = 0;
	Matrix<X_DIM> x;
	Matrix<X_DIM, X_DIM> SqrtSigma;

	for(int t = 0; t < T-1; ++t) {
		unVec(B[t], x, SqrtSigma);
		cost += alpha_belief*tr(SqrtSigma*SqrtSigma) + alpha_control*tr(~U[t]*U[t]);
	}
	unVec(B[T-1], x, SqrtSigma);
	cost += alpha_final_belief*tr(SqrtSigma*SqrtSigma);
	return cost;
}


void setupBeliefVars(beliefPenaltyMPC_params &problem, beliefPenaltyMPC_output &output)
{
	// problem inputs
	H = new beliefPenaltyMPC_FLOAT*[T];
	f = new beliefPenaltyMPC_FLOAT*[T-1];
	lb = new beliefPenaltyMPC_FLOAT*[T];
	ub = new beliefPenaltyMPC_FLOAT*[T];
	C = new beliefPenaltyMPC_FLOAT*[T-1];
	//D = new beliefPenaltyMPC_FLOAT*[T];
	e = new beliefPenaltyMPC_FLOAT*[T];

	// problem outputs
	z = new beliefPenaltyMPC_FLOAT*[T];

#define SET_VARS(n)    \
		H[ BOOST_PP_SUB(n,1) ] = problem.H##n ; \
		f[ BOOST_PP_SUB(n,1) ] = problem.f##n ;  \
		lb[ BOOST_PP_SUB(n,1) ] = problem.lb##n ;	\
		ub[ BOOST_PP_SUB(n,1) ] = problem.ub##n ;	\
		C[ BOOST_PP_SUB(n,1) ] = problem.C##n ;  \
		e[ BOOST_PP_SUB(n,1) ] = problem.e##n ;  \
		z[ BOOST_PP_SUB(n,1) ] = output.z##n ;
#define BOOST_PP_LOCAL_MACRO(n) SET_VARS(n)
#define BOOST_PP_LOCAL_LIMITS (1, TIMESTEPS-1)
#include BOOST_PP_LOCAL_ITERATE()


#define SET_LAST_VARS(n)    \
		H[ BOOST_PP_SUB(n,1) ] = problem.H##n ; \
		lb[ BOOST_PP_SUB(n,1) ] = problem.lb##n ;	\
		ub[ BOOST_PP_SUB(n,1) ] = problem.ub##n ;	\
		e[ BOOST_PP_SUB(n,1) ] = problem.e##n ;  \
		z[ BOOST_PP_SUB(n,1) ] = output.z##n ;
#define BOOST_PP_LOCAL_MACRO(n) SET_LAST_VARS(n)
#define BOOST_PP_LOCAL_LIMITS (TIMESTEPS, TIMESTEPS)
#include BOOST_PP_LOCAL_ITERATE()


/*
#define SET_D(n)    \
		D[ BOOST_PP_SUB(n,1) ] = problem.D##n ;
#define BOOST_PP_LOCAL_MACRO(n) SET_D(n)
#define BOOST_PP_LOCAL_LIMITS (2, TIMESTEPS)
#include BOOST_PP_LOCAL_ITERATE()
*/


	// initialize H in x'*H*x to penalize covariance and controls
	// H is diagonal

	int index;
	for(int t=0; t < T-1; ++t) {
		index = 0;
		for(int i=0; i < X_DIM; ++i) { H[t][index++] = 0; }
		for(int i=0; i < S_DIM; ++i) { H[t][index++] = 2*alpha_belief; }
		for(int i=0; i < U_DIM; ++i) { H[t][index++] = 2*alpha_control; }
		for(int i=0; i < B_DIM; ++i) { H[t][index++] = 0; }
		for(int i=0; i < B_DIM; ++i) { H[t][index++] = 0; }
	}

	index = 0;
	for(int i=0; i < X_DIM; ++i) { H[T-1][index++] = 0; }
	for(int i=0; i < S_DIM; ++i) { H[T-1][index++] = 2*alpha_final_belief; }


	/*
	// set up D

	Matrix<2*B_DIM, 3*B_DIM+U_DIM> D1mat = zeros<2*B_DIM, 3*B_DIM+U_DIM>();
	D1mat.insert(B_DIM, 0, (Matrix<B_DIM,B_DIM>)identity<B_DIM>());
	fillColMajor(D[1], D1mat);

	Matrix<B_DIM, 3*B_DIM+U_DIM> Dmat = zeros<B_DIM, 3*B_DIM+U_DIM>();
	Dmat.insert(0, 0, (Matrix<B_DIM,B_DIM>)identity<B_DIM>());
	for(int t=2; t < T-1; ++t) {
		fillColMajor(D[t], Dmat);
	}

	fillColMajor(D[T-1], (Matrix<B_DIM,B_DIM>)identity<B_DIM>());
	*/

}

void cleanupBeliefMPCVars()
{
	//delete[] H;
	delete[] f;
	delete[] lb;
	delete[] ub;
	delete[] C;
	//delete[] D;
	delete[] e;
	delete[] z;
}

double computeMerit(const std::vector< Matrix<B_DIM> >& B, const std::vector< Matrix<U_DIM> >& U, double penalty_coeff)
{
	double merit = 0;
	Matrix<X_DIM> x;
	Matrix<X_DIM, X_DIM> SqrtSigma;
	Matrix<B_DIM> dynviol;
	for(int t = 0; t < T-1; ++t) {
		unVec(B[t], x, SqrtSigma);
		merit += alpha_belief*tr(SqrtSigma*SqrtSigma) + alpha_control*tr(~U[t]*U[t]);
		dynviol = (B[t+1] - beliefDynamics(B[t], U[t]) );
		for(int i = 0; i < B_DIM; ++i) {
			merit += penalty_coeff*fabs(dynviol[i]);
		}
	}
	unVec(B[T-1], x, SqrtSigma);
	merit += alpha_final_belief*tr(SqrtSigma*SqrtSigma);
	return merit;
}

// TODO: Check if all inputs are valid, H, f, lb, ub, C, e, D at last time step
bool isValidInputs()
{
	for(int t = 0; t < T-1; ++t) {

		std::cout << "t: " << t << std::endl << std::endl;
		/*
		std::cout << "f: ";
		for(int i = 0; i < (3*B_DIM+U_DIM); ++i) {
			std::cout << f[t][i] << " ";
		}
		std::cout << std::endl;
		*/

		std::cout << "lb c: ";
		for(int i = 0; i < C_DIM; ++i) {
			std::cout << lb[t][i] << " ";
		}
		std::cout << std::endl;

		std::cout << "ub c: ";
		for(int i = 0; i < C_DIM; ++i) {
			std::cout << ub[t][i] << " ";
		}
		std::cout << std::endl;

		std::cout << "lb l: ";
		for(int i = 0; i < L_DIM; ++i) {
			std::cout << lb[t][C_DIM+i] << " ";
		}
		std::cout << std::endl;

		std::cout << "ub l: ";
		for(int i = 0; i < L_DIM; ++i) {
			std::cout << ub[t][C_DIM+i] << " ";
		}
		std::cout << std::endl;


		std::cout << "lb u: ";
		for(int i = 0; i < U_DIM; ++i) {
			std::cout << lb[t][B_DIM+i] << " ";
		}
		std::cout << std::endl;

		std::cout << "ub u: ";
		for(int i = 0; i < U_DIM; ++i) {
			std::cout << ub[t][B_DIM+i] << " ";
		}
		std::cout << std::endl;

		/*
		//std::cout << "ub s, t: ";
		//for(int i = 0; i < 2*B_DIM; ++i) {
		//	std::cout << ub[t][B_DIM+U_DIM+i] << " ";
		//}
		//std::cout << std::endl;

		std::cout << "C:" << std::endl;
		if (t == 0) {
			for(int i = 0; i < 2*B_DIM*(3*B_DIM+U_DIM); ++i) {
				std::cout << C[t][i] << " ";
			}
		} else {
			for(int i = 0; i < B_DIM*(3*B_DIM+U_DIM); ++i) {
				std::cout << C[t][i] << " ";
			}
		}
		std::cout << std::endl;

		std::cout << "e:" << std::endl;
		if (t == 0) {
			for(int i = 0; i < 2*B_DIM; ++i) {
				std::cout << e[t][i] << " ";
			}
		} else {
			for(int i = 0; i < B_DIM; ++i) {
				std::cout << e[t][i] << " ";
			}
		}

		std::cout << "e:" << std::endl;
		for(int i = 0; i < B_DIM; ++i) {
			std::cout << e[t][i] << " ";
		}
		*/
		std::cout << std::endl << std::endl;
	}
	/*
	std::cout << "e:" << std::endl;
	for(int i = 0; i < B_DIM; ++i) {
		std::cout << e[T-1][i] << " ";
	}
	*/

	std::cout << "lb b: ";
	for(int i = 0; i < B_DIM; ++i) {
		std::cout << lb[T-1][i] << " ";
	}
	std::cout << std::endl;

	std::cout << "ub b: ";
	for(int i = 0; i < B_DIM; ++i) {
		std::cout << ub[T-1][i] << " ";
	}
	std::cout << std::endl;
	std::cout << std::endl;
	return true;
}


bool minimizeMeritFunction(std::vector< Matrix<B_DIM> >& B, std::vector< Matrix<U_DIM> >& U, beliefPenaltyMPC_params &problem, beliefPenaltyMPC_output &output, beliefPenaltyMPC_info &info, double penalty_coeff, double trust_box_size)
{
	LOG_DEBUG("Solving sqp problem with penalty parameter: %2.4f", penalty_coeff);

	//Matrix<B_DIM,1> b0 = B[0];

	std::vector< Matrix<B_DIM,B_DIM> > F(T-1);
	std::vector< Matrix<B_DIM,U_DIM> > G(T-1);
	std::vector< Matrix<B_DIM> > h(T-1);

	double Beps = trust_box_size;
	double Ueps = trust_box_size;

	double optcost;

	std::vector<Matrix<B_DIM> > Bopt(T);
	std::vector<Matrix<U_DIM> > Uopt(T-1);

	double merit, model_merit, new_merit;
	double approx_merit_improve, exact_merit_improve, merit_improve_ratio;

	int sqp_iter = 1, index = 0;
	bool success;

	Matrix<B_DIM,B_DIM> IB = identity<B_DIM>();
	Matrix<B_DIM,B_DIM> minusIB = IB;
	for(int i = 0; i < B_DIM; ++i) {
		minusIB(i,i) = -1;
	}
	
	//Matrix<3*B_DIM+U_DIM, 3*B_DIM+U_DIM> HMat;
	Matrix<B_DIM,3*B_DIM+U_DIM> CMat;
	Matrix<B_DIM> eVec;
	//Matrix<S_DIM,S_DIM> Hess;

	// sqp loop
	while(true)
	{
		// In this loop, we repeatedly construct a linear approximation to the nonlinear belief dynamics constraint
		LOG_DEBUG("  sqp iter: %d", sqp_iter);

		merit = computeMerit(B, U, penalty_coeff);
		
		LOG_DEBUG("  merit: %4.10f", merit);

		// Problem linearization and definition
		// fill in H, f, C, e
		
		util::Timer constructCETimer;
		util::Timer_tic(&constructCETimer);

		util::Timer linearizeTimer, fillCETimer;
		double linearizeTime = 0, fillCETime = 0;
		for (int t = 0; t < T-1; ++t) 
		{
			Matrix<B_DIM>& bt = B[t];
			Matrix<U_DIM>& ut = U[t];

			util::Timer_tic(&linearizeTimer);
			linearizeBeliefDynamics(bt, ut, F[t], G[t], h[t]);
			linearizeTime += util::Timer_toc(&linearizeTimer);

			//std::cout << "F[" << t << "]" << std::endl << F[t] << std::endl;
			//std::cout << "G[" << t << "]" << std::endl << G[t] << std::endl;
			//std::cout << "h[" << t << "]" << std::endl << h[t] << std::endl;
			//int k;
			//std::cin >> k;

			
			util::Timer_tic(&fillCETimer);
			// initialize f in cost function to penalize
			// belief dynamics slack variables
			index = 0;
			for(int i = 0; i < (B_DIM+U_DIM); ++i) { f[t][index++] = 0; }
			for(int i = 0; i < 2*B_DIM; ++i) { f[t][index++] = penalty_coeff; }

			CMat.reset();
			eVec.reset();

			CMat.insert<B_DIM,B_DIM>(0,0,F[t]);
			CMat.insert<B_DIM,U_DIM>(0,B_DIM,G[t]);
			CMat.insert<B_DIM,B_DIM>(0,B_DIM+U_DIM,IB);
			CMat.insert<B_DIM,B_DIM>(0,2*B_DIM+U_DIM,minusIB);

			fillColMajor(C[t], CMat);

			if (t == 0) {
				eVec.insert<B_DIM,1>(0,0,B[0]);
				fillCol(e[0], eVec);
			} 
			
			eVec = -h[t] + F[t]*bt + G[t]*ut;
			fillCol(e[t+1], eVec);

			fillCETime += util::Timer_toc(&fillCETimer);
		}

		double constructCETime = util::Timer_toc(&constructCETimer);
		//std::cout << "construct CE time = " << constructCETime*1000 << "ms" << std::endl;
		//std::cout << "linearize time = " << linearizeTime*1000 << "ms" << std::endl;
		//std::cout << "fill CE time = " << fillCETime*1000 << "ms" << std::endl;
		

		// trust region size adjustment
		while(true)
		{
			//std::cout << "PAUSED INSIDE MINIMIZEMERITFUNCTION" << std::endl;
			//std::cin.ignore();


			LOG_DEBUG("       trust region size: %2.6f %2.6f", Beps, Ueps);

			util::Timer fillVarsTimer;
			util::Timer_tic(&fillVarsTimer);
			// solve the innermost QP here
			for(int t = 0; t < T-1; ++t)
			{
				Matrix<B_DIM>& bt = B[t];
				Matrix<U_DIM>& ut = U[t];

				// Fill in lb, ub

				index = 0;
				// x lower bound
				for(int i = 0; i < X_DIM; ++i) { lb[t][index++] = MAX(xMin[i], bt[i] - Beps); }
				// sigma lower bound
				for(int i = 0; i < S_DIM; ++i) { lb[t][index] = bt[index] - Beps; index++; }
				// u lower bound
				for(int i = 0; i < U_DIM; ++i) { lb[t][index++] = MAX(uMin[i], ut[i] - Ueps); }

				// for lower bound on L1 slacks
				for(int i = 0; i < 2*B_DIM; ++i) { lb[t][index++] = 0; }

				index = 0;
				// x upper bound
				for(int i = 0; i < X_DIM; ++i) { ub[t][index++] = MIN(xMax[i], bt[i] + Beps); }
				// sigma upper bound
				for(int i = 0; i < S_DIM; ++i) { ub[t][index] = bt[index] + Beps; index++; }
				// u upper bound
				for(int i = 0; i < U_DIM; ++i) { ub[t][index++] = MIN(uMax[i], ut[i] + Ueps); }

				//for(int i = 0; i < 2*B_DIM; ++i) { ub[t][index++] = INFTY; }
			}

			Matrix<B_DIM>& bT = B[T-1];

			// Fill in lb, ub, C, e
			index = 0;
			double delta = .1;
			// xGoal lower bound
			for(int i = 0; i < P_DIM; ++i) { lb[T-1][index++] = xGoal[i] - delta; }
			// loose on landmarks
			for(int i = 0; i < X_DIM-P_DIM; ++i) { lb[T-1][index++] = MAX(xMin[i+P_DIM], bT[i+P_DIM] - Beps); }
			// sigma lower bound
			for(int i = 0; i < S_DIM; ++i) { lb[T-1][index] = bT[index] - Beps; index++;}

			index = 0;
			// xGoal upper bound
			for(int i = 0; i < P_DIM; ++i) { ub[T-1][index++] = xGoal[i] + delta; }
			// loose on landmarks
			for(int i = 0; i < X_DIM-P_DIM; ++i) { ub[T-1][index++] = MIN(xMax[i+P_DIM], bT[i+P_DIM] + Beps); }
			// sigma upper bound
			for(int i = 0; i < S_DIM; ++i) { ub[T-1][index] = bT[index] + Beps; index++;}

			// Verify problem inputs
			//if (!isValidInputs()) {
			//	std::cout << "Inputs are not valid!" << std::endl;
			//	exit(-1);
			//}

			//std::cerr << "PAUSING INSIDE MINIMIZE MERIT FUNCTION FOR INPUT VERIFICATION" << std::endl;
			//int num;
			//std::cin >> num;

			double fillVarsTime = util::Timer_toc(&fillVarsTimer);
			//std::cout << "fill variables time = " << fillVarsTime*1000 << "ms" << std::endl;

			int exitflag = beliefPenaltyMPC_solve(&problem, &output, &info);
			if (exitflag == 1) {
				for(int t = 0; t < T-1; ++t) {
					Matrix<B_DIM>& bt = Bopt[t];
					Matrix<U_DIM>& ut = Uopt[t];

					for(int i = 0; i < B_DIM; ++i) {
						bt[i] = z[t][i];
					}
					for(int i = 0; i < U_DIM; ++i) {
						ut[i] = z[t][B_DIM+i];
					}
					optcost = info.pobj;
				}
				for(int i = 0; i < B_DIM; ++i) {
					Bopt[T-1][i] = z[T-1][i];
				}
			}
			else {
				LOG_ERROR("Some problem in solver");
				exit(-1);
			}

			LOG_DEBUG("Optimized cost: %4.10f", optcost);

			model_merit = optcost;
			new_merit = computeMerit(Bopt, Uopt, penalty_coeff);

			LOG_DEBUG("merit: %4.10f", merit);
			LOG_DEBUG("model_merit: %4.10f", model_merit);
			LOG_DEBUG("new_merit: %4.10f", new_merit);
			
			approx_merit_improve = merit - model_merit;
			exact_merit_improve = merit - new_merit;
			merit_improve_ratio = exact_merit_improve / approx_merit_improve;

			LOG_DEBUG("approx_merit_improve: %1.6f", approx_merit_improve);
			LOG_DEBUG("exact_merit_improve: %1.6f", exact_merit_improve);
			LOG_DEBUG("merit_improve_ratio: %1.6f", merit_improve_ratio);
			
			//std::cout << "PAUSED INSIDE minimizeMeritFunction" << std::endl;
			//int num;
			//std::cin >> num;

			if (approx_merit_improve < -1e-5) {
				LOG_ERROR("Approximate merit function got worse: %1.6f", approx_merit_improve);
				LOG_ERROR("Either convexification is wrong to zeroth order, or you are in numerical trouble");
				LOG_ERROR("Failure!");

				success = false;
			} else if (approx_merit_improve < cfg::min_approx_improve) {
				LOG_DEBUG("Converged: improvement small enough");
				B = Bopt; U = Uopt;
				return true;
			} else if ((exact_merit_improve < 0) || (merit_improve_ratio < cfg::improve_ratio_threshold)) {
				Beps *= cfg::trust_shrink_ratio;
				Ueps *= cfg::trust_shrink_ratio;
				LOG_DEBUG("Shrinking trust region size to: %2.6f %2.6f", Beps, Ueps);
			} else {
				Beps *= cfg::trust_expand_ratio;
				Ueps *= cfg::trust_expand_ratio;
				B = Bopt; U = Uopt;
				LOG_DEBUG("Accepted, Increasing trust region size to:  %2.6f %2.6f", Beps, Ueps);
				break;
			}

			if (Beps < cfg::min_trust_box_size && Ueps < cfg::min_trust_box_size) {
			    LOG_DEBUG("Converged: x tolerance");
			    return true;
			}

		} // trust region loop
		sqp_iter++;
	} // sqp loop

	return success;
}

double beliefPenaltyCollocation(std::vector< Matrix<B_DIM> >& B, std::vector< Matrix<U_DIM> >& U, beliefPenaltyMPC_params &problem, beliefPenaltyMPC_output &output, beliefPenaltyMPC_info &info)
{
	double penalty_coeff = cfg::initial_penalty_coeff;
	double trust_box_size = cfg::initial_trust_box_size;

	int penalty_increases = 0;

	Matrix<B_DIM> dynviol;

	// penalty loop
	while(penalty_increases < cfg::max_penalty_coeff_increases)
	{
		bool success = minimizeMeritFunction(B, U, problem, output, info, penalty_coeff, trust_box_size);

		double cntviol = 0;
		for(int t = 0; t < T-1; ++t) {
			dynviol = (B[t+1] - beliefDynamics(B[t], U[t]) );
			for(int i = 0; i < B_DIM; ++i) {
				cntviol += fabs(dynviol[i]);
			}
		}
	    success = success && (cntviol < cfg::cnt_tolerance);
	    LOG_DEBUG("Constraint violations: %2.10f",cntviol);
	    if (!success) {
	        penalty_increases++;
	        penalty_coeff = penalty_coeff*cfg::penalty_coeff_increase_ratio;
	        trust_box_size = cfg::initial_trust_box_size;
	    }
	    else {
	    	return computeCost(B, U);
	    }
	}
	return computeCost(B, U);
}



int main(int argc, char* argv[])
{
	const rlim_t stackSize = 32 * 1024 * 1024;   // min stack size = 32 MB
	struct rlimit rl;
	int result;

	result = getrlimit(RLIMIT_STACK, &rl);
	if (result == 0)
	{
		if (rl.rlim_cur < stackSize)
		{
			rl.rlim_cur = stackSize;
			result = setrlimit(RLIMIT_STACK, &rl);
			if (result != 0)
			{
				LOG_ERROR("setrlimit returned result = %d\n", result);
				exit(-1);
			}
		}
	} else {
		LOG_ERROR("Could not get stack limit");
		exit(-1);
	}

	LOG_INFO("Initializing problem parameters");
	initProblemParams();

	LOG_INFO("Setting up belief variables");
	std::cout << "size of params: " << sizeof(struct beliefPenaltyMPC_params)*9.53674e-7 << "mb" << std::endl;
	std::cout << "size of output: " << sizeof(struct beliefPenaltyMPC_output)*9.53674e-7 << "mb" << std::endl;
	std::cout << "size of info: " << sizeof(struct beliefPenaltyMPC_info)*9.53674e-7 << "mb" << std::endl;
	//beliefPenaltyMPC_params &problem = malloc(sizeof(struct beliefPenaltyMPC_params));
	//beliefPenaltyMPC_output *output = malloc(sizeof(struct beliefPenaltyMPC_output));
	//beliefPenaltyMPC_info *info = malloc(sizeof(struct beliefPenaltyMPC_info));
	beliefPenaltyMPC_params problem;
	beliefPenaltyMPC_output output;
	beliefPenaltyMPC_info info;
	setupBeliefVars(problem, output);
	util::Timer solveTimer;

	std::vector<Matrix<B_DIM> > B_total(T*NUM_WAYPOINTS);
	std::vector<Matrix<B_DIM> > B(T);

	Matrix<U_DIM> uinit;

	Matrix<X_DIM,1> xtmp;
	Matrix<X_DIM,X_DIM> stmp;
	for(int i=0; i < NUM_WAYPOINTS; ++i) {
		LOG_INFO("Going to waypoint %d",i);
		// goal is waypoint position + direct angle + landmarks
		xGoal.insert(0, 0, waypoints[i]);

		//xGoal[2] = x0[2];
		//x0[2] = atan2(xGoal[1] - x0[1], xGoal[0] - x0[0]);
		xGoal[2] = atan2(xGoal[1] - x0[1], xGoal[0] - x0[0]);


		xGoal.insert(C_DIM, 0, x0.subMatrix<L_DIM,1>(C_DIM,0));


		/*
		// initialize velocity to dist / timesteps
		uinit[0] = sqrt((x0[0] - xGoal[0])*(x0[0] - xGoal[0]) + (x0[1] - xGoal[1])*(x0[1] - xGoal[1])) / (double)((T-1)*DT);
		// angle already pointed at goal, so is 0
		uinit[1] = 0;

		std::vector<Matrix<U_DIM> > U(T-1, uinit);
		 */

		std::vector<Matrix<U_DIM> > U(T-1);
		bool initTrajSuccess = initTraj(x0.subMatrix<C_DIM,1>(0,0), xGoal.subMatrix<C_DIM,1>(0,0), U);
		if (!initTrajSuccess) {
			LOG_ERROR("Failed to initialize trajectory, exiting slam-belief");
			exit(-1);
		}

		//std::cout << "X car initial" << std::endl;
		vec(x0, SqrtSigma0, B[0]);
		for(int t=0; t < T-1; ++t) {
			//std::cout << ~B[t].subMatrix<C_DIM,1>(0,0);
			B[t+1] = beliefDynamics(B[t], U[t]);
		}
		//std::cout << ~B[T-1].subMatrix<C_DIM,1>(0,0) << std::endl << std::endl;
		unVec(B[T-1], xtmp, stmp);
		std::cout << stmp.subMatrix<P_DIM,P_DIM>(0,0) << std::endl;


		/*
		std::cout << "U" << std::endl;
		for(int t=0; t < T-1; ++t) {
			std::cout << ~U[t];
		}
		std::cout << std::endl;
		*/


		pythonDisplayTrajectory(B, U, waypoints, landmarks, T);


		double initTrajCost = computeCost(B, U);
		LOG_INFO("Initial trajectory cost: %4.10f", initTrajCost);

		Timer_tic(&solveTimer);

		double cost = beliefPenaltyCollocation(B, U, problem, output, info);

		double solvetime = util::Timer_toc(&solveTimer);

		vec(x0, SqrtSigma0, B[0]);
		for (size_t t = 0; t < T-1; ++t) {
			B[t+1] = beliefDynamics(B[t], U[t]);
			unVec(B[t+1], xtmp, stmp);

			//std::cout << stmp.subMatrix<P_DIM,P_DIM>(0,0) << std::endl;
		}
		std::cout << stmp.subMatrix<P_DIM,P_DIM>(0,0) << std::endl;

		LOG_INFO("Initial cost: %4.10f", initTrajCost);
		LOG_INFO("Optimized cost: %4.10f", cost);
		LOG_INFO("Actual cost: %4.10f", computeCost(B,U));
		LOG_INFO("Solve time: %5.3f ms", solvetime*1000);

		pythonDisplayTrajectory(B, U, waypoints, landmarks, T);

		unVec(B[T-1], x0, SqrtSigma0);

	}
	
	cleanupBeliefMPCVars();
	
	//vec(x0, SqrtSigma0, B[0]);
	//for (size_t t = 0; t < T-1; ++t) {
	//	B[t+1] = beliefDynamics(B[t], U[t]);
	//}

	return 0;
}
