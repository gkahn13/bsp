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

const double alpha_belief = 10; // 10;
const double alpha_final_belief = 10; // 10;
const double alpha_control = 1; // .1


namespace cfg {
const double improve_ratio_threshold = .1; // .1
const double min_approx_improve = 1e-2; // 1e-3
const double min_trust_box_size = 1e-3; // 1e-2

const double trust_shrink_ratio = .5; // .75
const double trust_expand_ratio = 1.25; // 1.25

const double cnt_tolerance = .5;
const double penalty_coeff_increase_ratio = 5; // 2
const double initial_penalty_coeff = 10; // 15

const double initial_trust_box_size = 5; // 10
const double initial_Xpos_trust_box_size = 5; // 10;
const double initial_Xangle_trust_box_size = M_PI/6; // 10*M_PI/6;
const double initial_Uvel_trust_box_size = 5; // 10;
const double initial_Uangle_trust_box_size = M_PI/8; // 10*M_PI/8;

const int max_penalty_coeff_increases = 3; // 4
const int max_sqp_iterations = 50; // 50
}

struct forces_exception {
	forces_exception() { }
};

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

double nearestAngleFromTo(double from, double to) {

	while (to > from) {
		if (to - 2*M_PI < from) {
			if (fabs(to - from) < (fabs(to - 2*M_PI - from))) {
				return to;
			} else {
				return to - 2*M_PI;
			}
		}
		to -= 2*M_PI;
	}

	while (to < from) {
		if (to + 2*M_PI > from) {
			if (fabs(to - from) < (fabs(to + 2*M_PI - from))) {
				return to;
			} else {
				return to + 2*M_PI;
			}
		}
		to += 2*M_PI;
	}

	// should never reach this point
	return to;
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
	D = new beliefPenaltyMPC_FLOAT*[T];
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



#define SET_D(n)    \
		D[ BOOST_PP_SUB(n,1) ] = problem.D##n ;
#define BOOST_PP_LOCAL_MACRO(n) SET_D(n)
#define BOOST_PP_LOCAL_LIMITS (2, TIMESTEPS)
#include BOOST_PP_LOCAL_ITERATE()


	// initialize H in x'*H*x to penalize covariance and controls
	// H is diagonal

	int index;
	for(int t=0; t < T-1; ++t) {
		index = 0;
		for(int i=0; i < X_DIM; ++i) { H[t][index++] = 0; }
		for(int i=0; i < S_DIM; ++i) { H[t][index++] = 2*alpha_belief; }
		for(int i=0; i < U_DIM; ++i) { H[t][index++] = 2*alpha_control; }
		// TODO: why does this work?
		for(int i=0; i < 2*B_DIM; ++i) { H[t][index++] = 5e2; } // 1e4
	}

	index = 0;
	for(int i=0; i < X_DIM; ++i) { H[T-1][index++] = 0; }
	for(int i=0; i < S_DIM; ++i) { H[T-1][index++] = 2*alpha_final_belief; }


	// set up D

	Matrix<B_DIM, 3*B_DIM+U_DIM> Dmat = zeros<B_DIM, 3*B_DIM+U_DIM>();
	Dmat.insert(0, 0, -(Matrix<B_DIM,B_DIM>)identity<B_DIM>());
	for(int t=1; t < T-1; ++t) {
		fillColMajor(D[t], Dmat);
	}

	fillColMajor(D[T-1], -(Matrix<B_DIM,B_DIM>)identity<B_DIM>());


}

void cleanupBeliefMPCVars()
{
	delete[] H;
	delete[] f;
	delete[] lb;
	delete[] ub;
	delete[] C;
	delete[] D;
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


bool minimizeMeritFunction(std::vector< Matrix<B_DIM> >& B, std::vector< Matrix<U_DIM> >& U, beliefPenaltyMPC_params &problem, beliefPenaltyMPC_output &output, beliefPenaltyMPC_info &info, double penalty_coeff)
{
	LOG_DEBUG("Solving sqp problem with penalty parameter: %2.4f", penalty_coeff);

	//Matrix<B_DIM,1> b0 = B[0];

	std::vector< Matrix<B_DIM,B_DIM> > F(T-1);
	std::vector< Matrix<B_DIM,U_DIM> > G(T-1);
	std::vector< Matrix<B_DIM> > h(T-1);

	double Beps = cfg::initial_trust_box_size;
	double Ueps = cfg::initial_trust_box_size;

	double Xpos_eps = cfg::initial_Xpos_trust_box_size;
	double Xangle_eps = cfg::initial_Xangle_trust_box_size;
	double Uvel_eps = cfg::initial_Uvel_trust_box_size;
	double Uangle_eps = cfg::initial_Uangle_trust_box_size;

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
		
		for (int t = 0; t < T-1; ++t) 
		{
			Matrix<B_DIM>& bt = B[t];
			Matrix<U_DIM>& ut = U[t];

			linearizeBeliefDynamics(bt, ut, F[t], G[t], h[t]);

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

		}
		

		// trust region size adjustment
		while(true)
		{
			//std::cout << "PAUSED INSIDE MINIMIZEMERITFUNCTION" << std::endl;
			//std::cin.ignore();


			LOG_DEBUG("       trust region size: %2.6f %2.6f", Beps, Ueps);

			// solve the innermost QP here
			for(int t = 0; t < T-1; ++t)
			{
				Matrix<B_DIM>& bt = B[t];
				Matrix<U_DIM>& ut = U[t];

				// Fill in lb, ub

				index = 0;
				// x lower bound
				//for(int i = 0; i < X_DIM; ++i) { lb[t][index++] = MAX(xMin[i], bt[i] - Beps); }

				// car pos lower bound
				for(int i = 0; i < P_DIM; ++i) { lb[t][index++] = MAX(xMin[i], bt[i] - Xpos_eps); }
				// car angle lower bound
				lb[t][index++] = MAX(xMin[P_DIM], bt[P_DIM] - Xangle_eps);
				// landmark pos lower bound
				for(int i = C_DIM; i < X_DIM; ++i) { lb[t][index++] = MAX(xMin[i], bt[i] - Xpos_eps); }

				// sigma lower bound
				for(int i = 0; i < S_DIM; ++i) { lb[t][index] = bt[index] - Beps; index++; }

				// u lower bound
				//for(int i = 0; i < U_DIM; ++i) { lb[t][index++] = MAX(uMin[i], ut[i] - Ueps); }

				// u velocity lower bound
				lb[t][index++] = MAX(uMin[0], ut[0] - Uvel_eps);
				// u angle lower bound
				lb[t][index++] = MAX(uMin[1], ut[1] - Uangle_eps);

				// for lower bound on L1 slacks
				for(int i = 0; i < 2*B_DIM; ++i) { lb[t][index++] = 0; }

				index = 0;
				// x upper bound
				//for(int i = 0; i < X_DIM; ++i) { ub[t][index++] = MIN(xMax[i], bt[i] + Beps); }

				// car pos upper bound
				for(int i = 0; i < P_DIM; ++i) { ub[t][index++] = MIN(xMax[i], bt[i] + Xpos_eps); }
				// car angle upper bound
				ub[t][index++] = MIN(xMax[P_DIM], bt[P_DIM] + Xangle_eps);
				// landmark pos upper bound
				for(int i = C_DIM; i < X_DIM; ++i) { ub[t][index++] = MIN(xMax[i], bt[i] + Xpos_eps); }

				// sigma upper bound
				for(int i = 0; i < S_DIM; ++i) { ub[t][index] = bt[index] + Beps; index++; }

				// u upper bound
				//for(int i = 0; i < U_DIM; ++i) { ub[t][index++] = MIN(uMax[i], ut[i] + Ueps); }

				// u velocity upper bound
				ub[t][index++] = MIN(uMax[0], ut[0] + Uvel_eps);
				// u angle upper bound
				ub[t][index++] = MIN(uMax[1], ut[1] + Uangle_eps);
			}

			Matrix<B_DIM>& bT = B[T-1];

			// Fill in lb, ub, C, e
			index = 0;
			double finalPosDelta = .1;
			double finalAngleDelta = M_PI/4;

			// xGoal lower bound
			for(int i = 0; i < P_DIM; ++i) { lb[T-1][index++] = xGoal[i] - finalPosDelta; }
			// loose on car angles and landmarks
			lb[T-1][index++] = nearestAngleFromTo(bT[2], xGoal[2] - finalAngleDelta);
			//lb[T-1][index++] = xGoal[2] - finalAngleDelta;
			for(int i = C_DIM; i < X_DIM; ++i) { lb[T-1][index++] = MAX(xMin[i], bT[i] - Xpos_eps); }
			// sigma lower bound
			for(int i = 0; i < S_DIM; ++i) { lb[T-1][index] = bT[index] - Beps; index++;}

			index = 0;
			// xGoal upper bound
			for(int i = 0; i < P_DIM; ++i) { ub[T-1][index++] = xGoal[i] + finalPosDelta; }
			// loose on car angles and landmarks
			ub[T-1][index++] = nearestAngleFromTo(bT[2], xGoal[2] + finalAngleDelta);
			//ub[T-1][index++] = xGoal[2] + finalAngleDelta;
			for(int i = C_DIM; i < X_DIM; ++i) { ub[T-1][index++] = MIN(xMax[i], bT[i] + Xpos_eps); }
			// sigma upper bound
			for(int i = 0; i < S_DIM; ++i) { ub[T-1][index] = bT[index] + Beps; index++;}

			// Verify problem inputs
			//if (!isValidInputs()) {
			//	std::cout << "Inputs are not valid!" << std::endl;
			//	exit(-1);
			//}


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
				throw forces_exception();
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
				//LOG_ERROR("Either convexification is wrong to zeroth order, or you are in numerical trouble");
				//LOG_ERROR("Failure!");

				return false;
			} else if (approx_merit_improve < cfg::min_approx_improve) {
				LOG_DEBUG("Converged: improvement small enough");
				B = Bopt; U = Uopt;
				return true;
			} else if ((exact_merit_improve < 0) || (merit_improve_ratio < cfg::improve_ratio_threshold)) {
				Beps *= cfg::trust_shrink_ratio;
				Ueps *= cfg::trust_shrink_ratio;
				Xpos_eps *= cfg::trust_shrink_ratio;
				Xangle_eps *= cfg::trust_shrink_ratio;
				Uvel_eps *= cfg::trust_shrink_ratio;
				Uangle_eps *= cfg::trust_shrink_ratio;
				LOG_DEBUG("Shrinking trust region size to: %2.6f %2.6f %2.6f %2.6f", Xpos_eps, Xangle_eps, Uvel_eps, Uangle_eps);
			} else {
				Beps *= cfg::trust_expand_ratio;
				Ueps *= cfg::trust_expand_ratio;
				Xpos_eps *= cfg::trust_expand_ratio;
				Xangle_eps *= cfg::trust_expand_ratio;
				Uvel_eps *= cfg::trust_expand_ratio;
				Uangle_eps *= cfg::trust_expand_ratio;
				B = Bopt; U = Uopt;
				LOG_DEBUG("Accepted, Increasing trust region size to:  %2.6f %2.6f", Beps, Ueps);
				break;
			}

			if (Beps < cfg::min_trust_box_size && Ueps < cfg::min_trust_box_size &&
					Xpos_eps < cfg::min_trust_box_size && Xangle_eps < cfg::min_trust_box_size &&
					Uvel_eps < cfg::min_trust_box_size && Uangle_eps < cfg::min_trust_box_size) {
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
	//double trust_box_size = cfg::initial_trust_box_size;

	int penalty_increases = 0;

	Matrix<B_DIM> dynviol;

	// penalty loop
	while(penalty_increases < cfg::max_penalty_coeff_increases)
	{
		bool success = minimizeMeritFunction(B, U, problem, output, info, penalty_coeff);

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
	        //trust_box_size = cfg::initial_trust_box_size;
	    }
	    else {
	    	return computeCost(B, U);
	    }
	}
	return computeCost(B, U);
}

void planPath(std::vector<Matrix<P_DIM> > l, beliefPenaltyMPC_params& problem, beliefPenaltyMPC_output& output, beliefPenaltyMPC_info& info, std::ofstream& f) {
	initProblemParams(l);

	util::Timer solveTimer, trajTimer;
	double totalSolveTime = 0, trajTime = 0;

	double totalTrajCost = 0;

	std::vector<Matrix<B_DIM> > B_total(T*NUM_WAYPOINTS);
	std::vector<Matrix<U_DIM> > U_total((T-1)*NUM_WAYPOINTS);

	std::vector<Matrix<B_DIM> > B(T);

	Matrix<U_DIM> uinit;

	Matrix<X_DIM,1> xtmp;
	Matrix<X_DIM,X_DIM> stmp;
	for(int i=0; i < NUM_WAYPOINTS; ++i) {
		LOG_INFO("Going to waypoint %d",i);
		// goal is waypoint position + direct angle + landmarks
		xGoal.insert(0, 0, waypoints[i]);


		// want to be facing the next waypoint
		if (i < NUM_WAYPOINTS - 1) {
			xGoal[2] = atan2(waypoints[i+1][1] - waypoints[i][1], waypoints[i+1][0] - waypoints[i][0]);
		} else {
			xGoal[2] = atan2(xGoal[1] - x0[1], xGoal[0] - x0[0]);
		}

		xGoal.insert(C_DIM, 0, x0.subMatrix<L_DIM,1>(C_DIM,0));

		util::Timer_tic(&trajTimer);

		std::vector<Matrix<U_DIM> > U(T-1);
		bool initTrajSuccess = initTraj(x0.subMatrix<C_DIM,1>(0,0), xGoal.subMatrix<C_DIM,1>(0,0), U, T);
		if (!initTrajSuccess) {
			LOG_ERROR("Failed to initialize trajectory, continuing anyways");
			//exit(-1);
		}

		double initTrajTime = util::Timer_toc(&trajTimer);
		trajTime += initTrajTime;

		vec(x0, SqrtSigma0, B[0]);
		for(int t=0; t < T-1; ++t) {
			B[t+1] = beliefDynamics(B[t], U[t]);
		}

		//pythonDisplayTrajectory(B, U, waypoints, landmarks, T, false);


		//double initTrajCost = computeCost(B, U);
		//LOG_INFO("Initial trajectory cost: %4.10f", initTrajCost);

		Timer_tic(&solveTimer);

		double cost = 0;
		int iter = 0;
		while(true) {
			try {
				cost = beliefPenaltyCollocation(B, U, problem, output, info);
				break;
			}
			catch (forces_exception &e) {
				if (iter > 3) {
					LOG_ERROR("Tried too many times, giving up");
					pythonDisplayTrajectory(U, T, false);
					logDataToFile(f, B_total, INFTY, INFTY, 1);
					//exit(-1);
					return;
				}
				LOG_ERROR("Forces exception, trying again");
				iter++;
			}
		}

		double solvetime = util::Timer_toc(&solveTimer);
		totalSolveTime += solvetime;

		vec(x0, SqrtSigma0, B[0]);
		for (int t = 0; t < T-1; ++t) {
			B[t+1] = beliefDynamics(B[t], U[t]);
			unVec(B[t+1], xtmp, stmp);
		}

		for (int t = 0; t < T-1; ++t) {
			B_total[t+T*i] = B[t];
			U_total[t+(T-1)*i] = U[t];
		}
		B_total[T-1+T*i] = B[T-1];

		totalTrajCost += computeCost(B,U);

//		LOG_INFO("Initial cost: %4.10f", initTrajCost);
//		LOG_INFO("Optimized cost: %4.10f", cost);
//		LOG_INFO("Actual cost: %4.10f", computeCost(B,U));
//		LOG_INFO("Solve time: %5.3f ms", solvetime*1000);

		//pythonDisplayTrajectory(B, U, waypoints, landmarks, T, true);

		unVec(B[T-1], x0, SqrtSigma0);

	}

	LOG_INFO("Total trajectory cost: %4.10f", totalTrajCost);
	LOG_INFO("Total trajectory solve time: %5.3f ms", trajTime*1000);
	LOG_INFO("Total solve time: %5.3f ms", totalSolveTime*1000);

	logDataToFile(f, B_total, totalSolveTime*1000, trajTime*1000, 0);

	pythonDisplayTrajectory(B_total, U_total, waypoints, landmarks, T*NUM_WAYPOINTS, false);
}

int main(int argc, char* argv[])
{
	LOG_INFO("Setting up belief variables");
	std::cout << "size of params: " << sizeof(struct beliefPenaltyMPC_params)*9.53674e-7 << "mb" << std::endl;
	std::cout << "size of output: " << sizeof(struct beliefPenaltyMPC_output)*9.53674e-7 << "mb" << std::endl;
	std::cout << "size of info: " << sizeof(struct beliefPenaltyMPC_info)*9.53674e-7 << "mb" << std::endl;

	beliefPenaltyMPC_params problem;
	beliefPenaltyMPC_output output;
	beliefPenaltyMPC_info info;
	setupBeliefVars(problem, output);

	std::vector<std::vector<Matrix<P_DIM> > > l_list = landmarks_list();

	std::ofstream f;
	logDataHandle("slam/data/slam-belief", f);

	for(int i=0; i < l_list.size(); ++i) {
		planPath(l_list[i], problem, output, info, f);
	}

	cleanupBeliefMPCVars();

	return 0;
}
