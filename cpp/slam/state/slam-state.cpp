#include <vector>
#include <iomanip>

#include "../slam.h"

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
	Matrix<X_DIM> x_tp1;
	
	Matrix<B_DIM> b;
	Matrix<X_DIM, X_DIM> SqrtSigma = SqrtSigma0;
	vec(X[t], SqrtSigma, b);
	for(int t = 0; t < T-1; ++t) {
	  cost += alpha_belief*tr(SqrtSigma*SqrtSigma) + alpha_control*tr(~U[t]*U[t]);
	  b = beliefDynamics(b, U[t]);
	  unVec(b, x_tp1, SqrtSigma);
	}
	unVec(b, x, SqrtSigma);
	cost += alpha_final_belief*tr(SqrtSigma*SqrtSigma);
	return cost;
}

// Jacobians: dg(b,u)/db, dg(b,u)/du
void linearizeBeliefDynamics(const Matrix<X_DIM>& x, const Matrix<U_DIM>& u, Matrix<B_DIM,B_DIM>& F, Matrix<B_DIM,U_DIM>& G, Matrix<B_DIM>& h)
{
	F.reset();
	Matrix<B_DIM> b;
	Matrix<S_DIM> SqrtSigma = SqrtSigma0;
	Matrix<B_DIM> br(b), bl(b);
	for (size_t i = 0; i < B_DIM; ++i) {
		br[i] += step; bl[i] -= step;
		//std::cout << "bplus: " << ~(beliefDynamics(br, u)) << std::endl;
		//std::cout << "bminus: " << ~(beliefDynamics(bl, u)) << std::endl;
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

double computeMerit(const std::vector< Matrix<B_DIM> >& X, const std::vector< Matrix<U_DIM> >& U, double penalty_coeff)
{
	double merit = 0;
	Matrix<X_DIM> x;
	Matrix<X_DIM> x_tp1;
	
	Matrix<B_DIM> b;
	Matrix<X_DIM, X_DIM> SqrtSigma = SqrtSigma0;
	Matrix<X_DIM> dynviol;
	vec(X[t], SqrtSigma, b);
	for(int t = 0; t < T-1; ++t) {
	  merit += alpha_belief*tr(SqrtSigma*SqrtSigma) + alpha_control*tr(~U[t]*U[t]);
	  b = beliefDynamics(b, U[t]);
	  unVec(b, x_tp1, SqrtSigma);
	  dynviol = (X[t+1]-x_tp1);
	  for(int i = 0; i < X_DIM; ++i) {
	    merit += penalty_coeff*fabs(dynviol[i]);
	  }
	}
	unVec(b, x, SqrtSigma);
	merit += alpha_final_belief*tr(SqrtSigma*SqrtSigma);
	return merit;
}

bool isValidInputs()
{
	for(int t = 0; t < T-1; ++t) {

		std::cout << "t: " << t << std::endl << std::endl;

		std::cout << "f: ";
		for(int i = 0; i < (3*B_DIM+U_DIM); ++i) {
			std::cout << f[t][i] << " ";
		}
		std::cout << std::endl;

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

		std::cout << std::endl << std::endl;
	}
	return true;
}

bool minimizeMeritFunction(std::vector< Matrix<X_DIM> >& X, std::vector< Matrix<U_DIM> >& U, beliefPenaltyMPC_params& problem, beliefPenaltyMPC_output& output, beliefPenaltyMPC_info& info, double penalty_coeff, double trust_box_size)
{
	LOG_DEBUG("Solving sqp problem with penalty parameter: %2.4f", penalty_coeff);

	// box constraint around goal
	double delta = 0.01;

	Matrix<X_DIM,1> x0 = X[0];

	std::vector< Matrix<X_DIM,X_DIM> > F(T-1);
	std::vector< Matrix<X_DIM,U_DIM> > G(T-1);
	std::vector< Matrix<X_DIM> > h(T-1);

	double Beps = trust_box_size;
	double Ueps = trust_box_size;

	double prevcost, optcost;

	std::vector<Matrix<X_DIM> > Xopt(T);
	std::vector<Matrix<U_DIM> > Uopt(T-1);

	double merit, model_merit, new_merit;
	double approx_merit_improve, exact_merit_improve, merit_improve_ratio;

	int sqp_iter = 1;
	bool success;

	Matrix<X_DIM,X_DIM> I = identity<X_DIM>();
	Matrix<X_DIM,X_DIM> minusI = I;
	for(int i = 0; i < X_DIM; ++i) {
		minusI(i,i) = -1;
	}

	// sqp loop
	while(true)
	{
		// In this loop, we repeatedly construct a linear approximation to the nonlinear belief dynamics constraint
		LOG_DEBUG("  sqp iter: %d", sqp_iter);

		merit = computeMerit(X, U, penalty_coeff);
		LOG_DEBUG("  merit: %4.10f", merit);

		for (int t = 0; t < T-1; ++t) {
		  linearizeDynamicsFunc(X[t], U[t], F[t], G[t], h[t]);
		}

		// trust region size adjustment
		while(true)
		{
			LOG_DEBUG("       trust region size: %2.6f %2.6f", Beps, Ueps);

			// compute Hessian first
			// so can force it to be PSD
			//jferguson -- Have to add in H values for cost approximation
			for (int t = 0; t < T-1; ++t) {
			  for (int i = 0; i < X_DIM; ++i) {
			    H[t][i] = resultCostGradDiagHess[1+dim+t*X_DIM+i];
			  }

			  for (int i = 0; i < U_DIM; ++i) {
			    H[t][X_DIM+i] = resultCostGradDiagHess[1+dim+T*X_DIM+t*U_DIM+i];
			  }

			}
			for (int i = 0; i < X_DIM; ++i) {
			  H[T-1][i] = resultCostGradDiagHess[1+dim+(T-1)*X_DIM+i];
			}

			//jferguson -- Going to have to add in the hessian constraints that are 
			//currently hardcoded in the forces setup

			forcePsdHessian(0);

			// solve the innermost QP here
			for(int t = 0; t < T-1; ++t)
			{
				Matrix<B_DIM>& xt = X[t];
				Matrix<U_DIM>& ut = U[t];


				Matrix<X_DIM+U_DIM> zbar;
				zbar.insert(0,0,xt);
				zbar.insert(X_DIM,0,ut);

				for(int i = 0; i < (X_DIM+U_DIM); ++i) {
					Hzbar[i] = H[t][i]*zbar[i];
				}

				// Fill in f, lb, ub, C, e
				//jferguson  -- This changes for cost approximation
				for(int i = 0; i < X_DIM; ++i) {
					f[t][i] = resultCostGradDiagHess[1+t*X_DIM+i];
				}
				for(int i = 0; i < U_DIM; ++i) {
					f[t][X_DIM+i] = resultCostGradDiagHess[1+T*X_DIM+t*U_DIM+i];
				}				

				for(int i = 0; i < 2*X_DIM; ++i) {
					f[t][X_DIM+U_DIM+i] = penalty_coeff;
				}

				for(int i = 0; i < (X_DIM+U_DIM); ++i) {
					hessian_constant += H[t][i]*zbar[i]*zbar[i];
					jac_constant += -(f[t][i]*zbar[i]);
				}

				for(int i = 0; i < (X_DIM+U_DIM); ++i) {
					f[t][i] -= Hzbar[i];
				}



				for (int i = 0; i < X_DIM; ++i) {
				  lb[t][i] = MAX(xMin[i], xt[i]-Xeps);
				}

				for (int i = 0; i < U_DIM; ++i) {
				  lb[t][X_DIM+i] = MAX(uMin[i], ut[i]-Ueps);
				}

				for(int i = 0; i < 2*X_DIM; ++i) {
					lb[t][X_DIM+U_DIM+i] = 0;
				}

				for (int i = 0; i < X_DIM; ++i) {
				  ub[t][i] = MIN(xMax[i], xt[i]+Xeps);
				}

				for (int i = 0; i < U_DIM; ++i) {
				  ub[t][X_DIM+i] = MIN(uMax[i], ut[i]+Ueps);
				}
				

				if (t > 0) {
				  //jferguson -- Equality constraints remain the same
					Matrix<X_DIM,3*X_DIM+U_DIM> CMat;
					Matrix<X_DIM> eVec;

					CMat.insert<X_DIM,X_DIM>(0,0,F[t]);
					CMat.insert<X_DIM,U_DIM>(0,X_DIM,G[t]);
					CMat.insert<X_DIM,X_DIM>(0,X_DIM+U_DIM,I);
					CMat.insert<X_DIM,X_DIM>(0,2*X_DIM+U_DIM,minusI);

					//std::cout << CMat << std::endl;

					int idx = 0;
					int nrows = CMat.numRows(), ncols = CMat.numColumns();
					for(int c = 0; c < ncols; ++c) {
						for(int r = 0; r < nrows; ++r) {
							C[t][idx++] = CMat[c + r*ncols];
						}
					}
					eVec = -h[t] + F[t]*bt + G[t]*ut;
					int nelems = eVec.numRows();
					for(int i = 0; i < nelems; ++i) {
						e[t][i] = eVec[i];
					}
				}
				else {
					Matrix<2*X_DIM,3*X_DIM+U_DIM> CMat;
					Matrix<2*X_DIM> eVec;

					CMat.insert<X_DIM,X_DIM>(0,0,I);
					CMat.insert<X_DIM,U_DIM>(0,X_DIM,zeros<X_DIM,U_DIM>());
					CMat.insert<X_DIM,2*X_DIM>(0,X_DIM+U_DIM,zeros<X_DIM,2*X_DIM>());

					CMat.insert<X_DIM,X_DIM>(X_DIM,0,F[t]);
					CMat.insert<X_DIM,U_DIM>(X_DIM,X_DIM,G[t]);
					CMat.insert<X_DIM,X_DIM>(X_DIM,X_DIM+U_DIM,I);
					CMat.insert<X_DIM,X_DIM>(X_DIM,2*X_DIM+U_DIM,minusI);

					int idx = 0;
					int nrows = CMat.numRows(), ncols = CMat.numColumns();
					for(int c = 0; c < ncols; ++c) {
						for(int r = 0; r < nrows; ++r) {
							C[t][idx++] = CMat[c + r*ncols];
						}
					}
					eVec.insert<X_DIM,1>(0,0,x0);
					eVec.insert<X_DIM,1>(X_DIM,0,zeros<X_DIM,1>());
					int nelems = eVec.numRows();
					for(int i = 0; i < nelems; ++i) {
						e[t][i] = eVec[i];
					}
				}
			}

			} //setting up problem

			// Last stage
			Matrix<X_DIM>& xT = X[T-1];

			
			for (int i = 0; i < X_DIM; ++i) {
			  f[T-1][i] = resultCostGradDiagHess[1+(T-1)*X_DIM+i];
			}

			for(int i = 0; i < X_DIM; ++i) {
				hessian_constant += H[T-1][i]*xT[i]*xT[i];
				jac_constant += -(f[T-1][i]*xT[i]);
			}

			for(int i = 0; i < X_DIM; ++i) {
				Hzbar[i] = H[T-1][i]*xT[i];
				f[T-1][i] -= Hzbar[i];
			}

			//e[T-1][0] = 0; e[T-1][1] = 0;

			constant_cost = 0.5*hessian_constant + jac_constant + resultCostGradDiagHess[0];
			
			for (int i=0; i<X_DIM; ++i) {
			  lb[T-1][i] = MAX(xGoal[i] - delta, xT[i] - Beps);
			}

			for (int i=0; i<X_DIM; ++i) {
			  ub[T-1][i] = MIN(xGoal[i] + delta, xT[i] + Beps);
			}

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
					Matrix<B_DIM>& xt = Xopt[t];
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
				X = Xopt; U = Uopt;
				return true;
			} else if ((exact_merit_improve < 0) || (merit_improve_ratio < cfg::improve_ratio_threshold)) {
				Xeps *= cfg::trust_shrink_ratio;
				Ueps *= cfg::trust_shrink_ratio;
				LOG_DEBUG("Shrinking trust region size to: %2.6f %2.6f", Beps, Ueps);
			} else {
				Xeps *= cfg::trust_expand_ratio;
				Ueps *= cfg::trust_expand_ratio;
				X = Xopt; U = Uopt;
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

double statePenaltyCollocation(std::vector< Matrix<B_DIM> >& B, std::vector< Matrix<U_DIM> >& U, beliefPenaltyMPC_params& problem, beliefPenaltyMPC_output& output, beliefPenaltyMPC_info& info)
{
	double penalty_coeff = cfg::initial_penalty_coeff;
	double trust_box_size = cfg::initial_trust_box_size;

	int penalty_increases = 0;

	Matrix<X_DIM> x_tp1;
	
	Matrix<B_DIM> b;
	Matrix<B_DIM> b_tp1;
	Matrix<X_DIM, X_DIM> SqrtSigma = SqrtSigma0;	

	Matrix<X_DIM> dynviol;

	// penalty loop
	while(penalty_increases < cfg::max_penalty_coeff_increases)
	{
		bool success = minimizeMeritFunction(X, U, problem, output, info, penalty_coeff, trust_box_size);

		double cntviol = 0;
		vec(X[t], SqrtSigma, b);		
		for(int t = 0; t < T-1; ++t) {
		  b = beliefDynamics(b, U[t]);
		  unVec(b, x_tp1, SqrtSigma);
		  dynviol = (X[t+1]-x_tp1);
		  for(int i = 0; i < X_DIM; ++i) {
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
	    	return computeCost(X, U);
	    }
	}
	return computeCost(X, U);
}


int main(int argc, char* argv[])
{
  x0[0] = 0; x0[1] = 0; x0[2] = ;
	SqrtSigma0 = identity<X_DIM>();
	//What to initialize robot covariance with?
	//What to initialize landmarks with?
	//How to include multiple goals (separate trajectories for each one?)
	xGoal[0] = -3.5; xGoal[1] = -2;

	//Bounds to use? robot -> Bit beyond all landmarks+goals
	//               landmarks -> ??? (same?)
	//               Control bounds -> arbitrary
	xMin[0] = -5; xMin[1] = -3; 
	xMax[0] = 5; xMax[1] = 3;
	uMin[0] = -1; uMin[1] = -1;
	uMax[0] = 1; uMax[1] = 1;

	Matrix<U_DIM> uinit;
	uinit[0] = (xGoal[0] - x0[0]) / (T-1);
	uinit[1] = (xGoal[1] - x0[1]) / (T-1);
	
	std::vector<Matrix<U_DIM> > U(T-1, uinit); 

	std::vector<Matrix<X_DIM> > X(T);

	Matrix<B_DIM> b;
	Matrix<S_DIM> SqrtSigma;

	vec(x0, SqrtSigma0, b);
	for (size_t t = 0; t < T-1; ++t) {
		b = beliefDynamics(b, U[t]);
		unVec(b, X[t+1], SqrtSigma);
		//std::cout << ~B[t] << std::endl;
	}

	//for (size_t t = 0; t < T; ++t) {
	//	std::cout << ~B[t];
	//}

	beliefPenaltyMPC_params problem;
	beliefPenaltyMPC_output output;
	beliefPenaltyMPC_info info;

	setupBeliefVars(problem, output);

	util::Timer solveTimer;
	util::Timer_tic(&solveTimer);
	
	// B&U optimized in-place
	double cost = statePenaltyCollocation(X, U, problem, output, info);

	double solvetime = util::Timer_toc(&solveTimer);
	LOG_INFO("Optimized cost: %4.10f", cost);
	LOG_INFO("Solve time: %5.3f ms", solvetime*1000);
	
	cleanupBeliefMPCVars();

	pythonDisplayTrajectory(B, U);

	/*
	for (size_t t = 0; t < T; ++t) {
		std::cout << ~B[t] << std::endl;
	}
	*/

	int k;
	std::cin >> k;

	//CAL_End();
	return 0;
}
