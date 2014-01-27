#include <vector>
#include <iomanip>

#include "../parameter.h"

#include "util/matrix.h"
#include "util/Timer.h"
#include "util/logging.h"
//#include "util/utils.h"

#include <Python.h>
#include <boost/python.hpp>
#include <boost/filesystem.hpp>

namespace py = boost::python;

#define BELIEF_PENALTY_MPC
//#define BELIEF_MPC


extern "C" {
#include "beliefPenaltyMPC.h"
beliefPenaltyMPC_FLOAT **f, **lb, **ub, **C, **e, **z;
}

Matrix<X_DIM> x0;
Matrix<X_DIM,X_DIM> SqrtSigma0;
Matrix<X_DIM> xGoal;
Matrix<X_DIM> xMin, xMax;
Matrix<U_DIM> uMin, uMax;

namespace cfg {
const double improve_ratio_threshold = .1;
const double min_approx_improve = 1e-4;
const double min_trust_box_size = 1e-3;
const double trust_shrink_ratio = .1;
const double trust_expand_ratio = 1.5;
const double cnt_tolerance = 1e-4;
const double penalty_coeff_increase_ratio = 10;
const double initial_penalty_coeff = 10;
const double initial_trust_box_size = 1;
const int max_penalty_coeff_increases = 2;
const int max_sqp_iterations = 50;
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

// Jacobians: dg(b,u)/db, dg(b,u)/du
void linearizeBeliefDynamics(const Matrix<B_DIM>& b, const Matrix<U_DIM>& u, Matrix<B_DIM,B_DIM>& F, Matrix<B_DIM,U_DIM>& G, Matrix<B_DIM>& h)
{
	F.reset();
	Matrix<B_DIM> br(b), bl(b);
	for (size_t i = 0; i < B_DIM; ++i) {
		br[i] += step; bl[i] -= step;
		F.insert(0,i, (beliefDynamics(br, u) - beliefDynamics(bl, u)) / (br[i] - bl[i]));
		br[i] = b[i]; bl[i] = b[i];
	}

	G.reset();
	Matrix<U_DIM> ur(u), ul(u);
	for (size_t i = 0; i < U_DIM; ++i) {
		ur[i] += step; ul[i] -= step;
		G.insert(0,i, (beliefDynamics(b, ur) - beliefDynamics(b, ul)) / (ur[i] - ul[i]));
		ur[i] = u[i]; ul[i] = u[i];
	}

	h = beliefDynamics(b, u);
}


void setupBeliefVars(beliefPenaltyMPC_params& problem, beliefPenaltyMPC_output& output)
{
	f = new beliefPenaltyMPC_FLOAT*[T-1];
	lb = new beliefPenaltyMPC_FLOAT*[T];
	ub = new beliefPenaltyMPC_FLOAT*[T];
	C = new beliefPenaltyMPC_FLOAT*[T-1];
	e = new beliefPenaltyMPC_FLOAT*[T-1];
	z = new beliefPenaltyMPC_FLOAT*[T];

#define SET_VARS(n)    \
		f[ BOOST_PP_SUB(n,1) ] = problem.f##n ;  \
		C[ BOOST_PP_SUB(n,1) ] = problem.C##n ;  \
		e[ BOOST_PP_SUB(n,1) ] = problem.e##n ;  \
		lb[ BOOST_PP_SUB(n,1) ] = problem.lb##n ;	\
		ub[ BOOST_PP_SUB(n,1) ] = problem.ub##n ;	\
		z[ BOOST_PP_SUB(n,1) ] = output.z##n ;

#define BOOST_PP_LOCAL_MACRO(n) SET_VARS(n)
#define BOOST_PP_LOCAL_LIMITS (1, TIMESTEPS-1)
#include BOOST_PP_LOCAL_ITERATE()


#define SET_LAST_VARS(n)    \
		lb[ BOOST_PP_SUB(n,1) ] = problem.lb##n ;	\
		ub[ BOOST_PP_SUB(n,1) ] = problem.ub##n ;	\
		z[ BOOST_PP_SUB(n,1) ] = output.z##n ;

#define BOOST_PP_LOCAL_MACRO(n) SET_LAST_VARS(n)
#define BOOST_PP_LOCAL_LIMITS (TIMESTEPS, TIMESTEPS)
#include BOOST_PP_LOCAL_ITERATE()

}

void cleanupBeliefMPCVars()
{
	delete[] f;
	delete[] lb;
	delete[] ub;
	delete[] C;
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


		std::cout << "lb b: ";
		for(int i = 0; i < B_DIM; ++i) {
			std::cout << lb[t][i] << " ";
		}
		std::cout << std::endl;


		std::cout << "lb u: ";
		for(int i = 0; i < U_DIM; ++i) {
			std::cout << lb[t][B_DIM+i] << " ";
		}
		std::cout << std::endl;

		std::cout << "lb s, t: ";
		for(int i = 0; i < 2*B_DIM; ++i) {
			std::cout << lb[t][B_DIM+U_DIM+i] << " ";
		}
		std::cout << std::endl;

		std::cout << "ub b: ";
		for(int i = 0; i < B_DIM; ++i) {
			std::cout << ub[t][i] << " ";
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
			for(int i = 0; i < 170; ++i) {
				std::cout << C[t][i] << " ";
			}
		} else {
			for(int i = 0; i < 85; ++i) {
				std::cout << C[t][i] << " ";
			}
		}
		std::cout << std::endl;

		std::cout << "e:" << std::endl;
		if (t == 0) {
			for(int i = 0; i < 10; ++i) {
				std::cout << e[t][i] << " ";
			}
		} else {
			for(int i = 0; i < 5; ++i) {
				std::cout << e[t][i] << " ";
			}
		}
		*/
		std::cout << std::endl << std::endl;
	}
	return true;
}

// utility to fill Matrix in column major format in FORCES array

template <size_t _numRows, size_t _numColumns>
inline void fillColMajor(double *X, const Matrix<_numRows, _numColumns>& XMat) {
	int idx = 0;
	for(size_t c = 0; c < _numColumns; ++c) {
		for(size_t r = 0; r < _numRows; ++r) {
			X[idx++] = XMat[c + r*_numColumns];
		}
	}
}

bool minimizeMeritFunction(std::vector< Matrix<B_DIM> >& B, std::vector< Matrix<U_DIM> >& U, beliefPenaltyMPC_params& problem, beliefPenaltyMPC_output& output, beliefPenaltyMPC_info& info, double penalty_coeff, double trust_box_size)
{
	LOG_DEBUG("Solving sqp problem with penalty parameter: %2.4f", penalty_coeff);

	// box constraint around goal
	double delta = 0.01;

	Matrix<B_DIM,1> b0 = B[0];

	std::vector< Matrix<B_DIM,B_DIM> > F(T-1);
	std::vector< Matrix<B_DIM,U_DIM> > G(T-1);
	std::vector< Matrix<B_DIM> > h(T-1);

	double Beps = trust_box_size;
	double Ueps = trust_box_size;

	double prevcost, optcost;

	std::vector<Matrix<B_DIM> > Bopt(T);
	std::vector<Matrix<U_DIM> > Uopt(T-1);

	double merit, model_merit, new_merit;
	double approx_merit_improve, exact_merit_improve, merit_improve_ratio;

	int sqp_iter = 1, index = 0;
	bool success;

	Matrix<B_DIM,B_DIM> I = identity<B_DIM>();
	Matrix<B_DIM,B_DIM> minusI = I;
	for(int i = 0; i < B_DIM; ++i) {
		minusI(i,i) = -1;
	}

	// sqp loop
	while(true)
	{
		// In this loop, we repeatedly construct a linear approximation to the nonlinear belief dynamics constraint
		LOG_DEBUG("  sqp iter: %d", sqp_iter);

		merit = computeMerit(B, U, penalty_coeff);
		LOG_DEBUG("  merit: %4.10f", merit);

		for (int t = 0; t < T-1; ++t) {
			linearizeBeliefDynamics(B[t], U[t], F[t], G[t], h[t]);

			// Fill in f TODO: does this need H, like arm?
			for(int i = 0; i < (B_DIM+U_DIM); ++i) {
				f[t][i] = 0;
			}
			for(int i = 0; i < 2*B_DIM; ++i) {
				f[t][B_DIM+U_DIM+i] = penalty_coeff;
			}
		}

		// trust region size adjustment
		while(true)
		{
			LOG_DEBUG("       trust region size: %2.6f %2.6f", Beps, Ueps);

			// solve the innermost QP here
			for(int t = 0; t < T-1; ++t)
			{
				Matrix<B_DIM>& bt = B[t];
				Matrix<U_DIM>& ut = U[t];

				index = 0;
				// x lower bound
				for(int i = 0; i < X_DIM; ++i) { lb[t][index++] = MAX(xMin[i], bt[i] - Beps); }
				// sigma lower bound
				for(int i = 0; i < S_DIM; ++i) { lb[t][index] = bt[index] - Beps; index++; }
				// u lower bound
				for(int i = 0; i < U_DIM; ++i) { lb[t][index++] = MAX(uMin[i], ut[i] - Ueps); }

				// for lower bound on slacks
				for(int i = 0; i < 2*B_DIM; ++i) { lb[t][index++] = 0; }

				index = 0;
				// x upper bound
				for(int i = 0; i < X_DIM; ++i) { ub[t][index++] = MIN(xMax[i], bt[i] + Beps); }
				// sigma upper bound
				for(int i = 0; i < S_DIM; ++i) { ub[t][index] = bt[index] + Beps; index++; }
				// u upper bound
				for(int i = 0; i < U_DIM; ++i) { ub[t][index++] = MIN(uMax[i], ut[i] + Ueps); }

				//for(int i = 0; i < 2*B_DIM; ++i) { ub[t][index++] = INFTY; }

				if (t > 0) {
					Matrix<B_DIM,3*B_DIM+U_DIM> CMat;
					Matrix<B_DIM> eVec;

					CMat.insert<B_DIM,B_DIM>(0,0,F[t]);
					CMat.insert<B_DIM,U_DIM>(0,B_DIM,G[t]);
					CMat.insert<B_DIM,B_DIM>(0,B_DIM+U_DIM,I);
					CMat.insert<B_DIM,B_DIM>(0,2*B_DIM+U_DIM,minusI);

					fillColMajor(C[t], CMat);

					eVec = -h[t] + F[t]*bt + G[t]*ut;
					int nelems = eVec.numRows();
					for(int i = 0; i < nelems; ++i) {
						e[t][i] = eVec[i];
					}
				}
				else {
					Matrix<2*B_DIM,3*B_DIM+U_DIM> CMat;
					Matrix<2*B_DIM> eVec;

					CMat.insert<B_DIM,B_DIM>(0,0,I);
					CMat.insert<B_DIM,U_DIM>(0,B_DIM,zeros<B_DIM,U_DIM>());
					CMat.insert<B_DIM,2*B_DIM>(0,B_DIM+U_DIM,zeros<B_DIM,2*B_DIM>());

					CMat.insert<B_DIM,B_DIM>(B_DIM,0,F[t]);
					CMat.insert<B_DIM,U_DIM>(B_DIM,B_DIM,G[t]);
					CMat.insert<B_DIM,B_DIM>(B_DIM,B_DIM+U_DIM,I);
					CMat.insert<B_DIM,B_DIM>(B_DIM,2*B_DIM+U_DIM,minusI);

					fillColMajor(C[t], CMat);

					eVec.insert<B_DIM,1>(0,0,b0);
					eVec.insert<B_DIM,1>(B_DIM,0,zeros<B_DIM,1>());
					int nelems = eVec.numRows();
					for(int i = 0; i < nelems; ++i) {
						e[t][i] = eVec[i];
					}
				}
			}

			Matrix<B_DIM>& bT = B[T-1];

			// Fill in lb, ub, C, e
			index = 0;
			// xGoal lower bound
			for(int i = 0; i < X_DIM; ++i) { lb[T-1][index++] = MAX(xGoal[i] - delta, bT[i] - Beps); }
			// sigma lower bound
			for(int i = 0; i < S_DIM; ++i) { lb[T-1][index] = bT[index] - Beps; index++; }

			index = 0;
			// xGoal upper bound
			for(int i = 0; i < X_DIM; ++i) { ub[T-1][index++] = MIN(xGoal[i] + delta, bT[i] + Beps); }
			// sigma lower bound
			for(int i = 0; i < S_DIM; ++i) { ub[T-1][index] = bT[index] + Beps; index++; }

			// Verify problem inputs
			//if (!isValidInputs()) {
			//	std::cout << "Inputs are not valid!" << std::endl;
			//	exit(-1);
			//}


			//std::cerr << "PAUSING INSIDE MINIMIZE MERIT FUNCTION FOR INPUT VERIFICATION" << std::endl;
			//int num;
			//std::cin >> num;

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
				std::exit(-1);
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
				prevcost = optcost;
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

double beliefPenaltyCollocation(std::vector< Matrix<B_DIM> >& B, std::vector< Matrix<U_DIM> >& U, beliefPenaltyMPC_params& problem, beliefPenaltyMPC_output& output, beliefPenaltyMPC_info& info)
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
	double length1_est = .05, // inverse = 20
			length2_est = .05, // inverse = 20
			mass1_est = .105, // inverse = 9.52
			mass2_est = .089; // inverse = 11.24

	// position, then velocity
	x0[0] = -M_PI/2.0; x0[1] = -M_PI/2.0; x0[2] = 0; x0[3] = 0;
	// parameter start estimates (alphabetical, then numerical order)
	x0[4] = 1/length1_est; x0[5] = 1/length2_est; x0[6] = 1/mass1_est; x0[7] = 1/mass2_est;

	xGoal[0] = -M_PI/2.0; xGoal[1] = -M_PI/2.0; xGoal[2] = 0.0; xGoal[3] = 0.0;
	xGoal[4] = 1/length1_est; xGoal[5] = 1/length2_est; xGoal[6] = 1/mass1_est; xGoal[7] = 1/mass2_est;

	// from original file, possibly change
	SqrtSigma0(4,4) = sqrt(0.5);
	SqrtSigma0(5,5) = sqrt(0.5);
	SqrtSigma0(6,6) = 1.0;
	SqrtSigma0(7,7) = 1.0;

	xMin[0] = -M_PI; // joint pos 1
	xMin[1] = -M_PI; // joint pos 2
	xMin[2] = -M_PI/2; // joint vel 1
	xMin[3] = -M_PI/2; // joint vel 2
	xMin[4] = 0; // 1/length1
	xMin[5] = 0; // 1/length2
	xMin[6] = 0; // 1/mass1
	xMin[7] = 0; // 1/mass2

	xMax[0] = M_PI; // joint pos 1
	xMax[1] = M_PI; // joint pos 2
	xMax[2] = M_PI/2; // joint vel 1
	xMax[3] = M_PI/2; // joint vel 2
	xMax[4] = 1/.01; // 1/length1
	xMax[5] = 1/.01; // 1/length2
	xMax[6] = 1/.01; // 1/mass1
	xMax[7] = 1/.01; // 1/mass2

	for(int i = 0; i < U_DIM; ++i) {
		uMin[i] = -1;
		uMax[i] = 1;
	}

	Matrix<U_DIM> uinit = (xGoal.subMatrix<U_DIM,1>(0,0) - x0.subMatrix<U_DIM,1>(0,0))/(double)(T-1);
	
	std::vector<Matrix<U_DIM> > U(T-1, uinit); 

	std::vector<Matrix<B_DIM> > B(T);

	vec(x0, SqrtSigma0, B[0]);
	for (size_t t = 0; t < T-1; ++t) {
		B[t+1] = beliefDynamics(B[t], U[t]);
		//std::cout << ~B[t] << std::endl;
	}

	/*
	std::cout << "initial X:" << std::endl;
	for (int t = 0; t < T; ++t) {
		std::cout << ~B[t].subMatrix<X_DIM>(0,0);
	}

	std::cout << "initial U:" << std::endl;
	for (int t = 0; t < T; ++t) {
		std::cout << ~U[t];
	}
	*/


	beliefPenaltyMPC_params problem;
	beliefPenaltyMPC_output output;
	beliefPenaltyMPC_info info;

	setupBeliefVars(problem, output);

	util::Timer solveTimer;
	util::Timer_tic(&solveTimer);
	
	// B&U optimized in-place

	double cost = beliefPenaltyCollocation(B, U, problem, output, info);

	double solvetime = util::Timer_toc(&solveTimer);
	LOG_INFO("Optimized cost: %4.10f", cost);
	LOG_INFO("Solve time: %5.3f ms", solvetime*1000);
	
	cleanupBeliefMPCVars();

	for (size_t t = 0; t < T; ++t) {
		std::cout << ~B[t].subMatrix<X_DIM,1>(0,0) << std::endl;
	}


#define PLOT
#ifdef PLOT
	pythonDisplayTrajectory(U, SqrtSigma0, x0, xGoal);

	//int k;
	//std::cin >> k;
#endif

	return 0;
}
