#include <vector>
#include <iomanip>

#include "callisto.h"
#include "matrix.h"
#include "utils.h"

extern "C" {
#include "arm-belief-pos-goal-MPC.h"
beliefPenaltyMPC_FLOAT **Q, **f, **lb, **ub, **C, **e, **z;
beliefPenaltyMPC_FLOAT *A, *b;
}

#include "boost/preprocessor.hpp"

#include "arm.h"

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

void constructHessian(const Matrix<B_DIM>& b, Matrix<S_DIM,S_DIM>& Hess) {

	Matrix<X_DIM> x = b.subMatrix<X_DIM,1>(0,0);
	Matrix<S_DIM> SqrtSigma = b.subMatrix<S_DIM,1>(X_DIM,0);

	Matrix<G_DIM,X_DIM> J;
	linearizeg(x,J);
	Matrix<X_DIM,X_DIM> JTJ = ~J*J;

	// construct the full Hessian
	Matrix<X_DIM*X_DIM,X_DIM*X_DIM> Hfull;
	Hfull.reset();

	for(int i = 0; i < X_DIM; ++i) {
		Hfull.insert<X_DIM,X_DIM>(i*X_DIM,i*X_DIM,JTJ);
	}

	Matrix<X_DIM*X_DIM,S_DIM> A;
	A.reset();
	int idx = 0;
	for (int i = 0; i < X_DIM; ++i) {
	    for(int j = 0; j < X_DIM; ++j) {
	        if (i <= j) {
	            A(idx,(2*X_DIM-i+1)*i/2+j-i) = 1;
	        } else {
	            A(idx,(2*X_DIM-j+1)*j/2+i-j) = 1;
	        }
	        idx++;
	    }
	}

	Hess = ~A*Hfull*A;
	
	//std::cout << SqrtSigma << std::endl;
	//std::ofstream fptr("Hess.txt",std::ios::out);
	//fptr << Hess;
	//fptr.close();
	//std::cout << Hess << std::endl;
	//std::cout << "PAUSED INSIDE CONSTRUCT HESSIAN" << std::endl;
	//int k;
	//std::cin >> k;
	
	//tracecost = (~SqrtSigma*H*SqrtSigma)[0];
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

double computeCost(const std::vector< Matrix<B_DIM> >& B, const std::vector< Matrix<U_DIM> >& U)
{
	double cost = 0;
	Matrix<X_DIM> x;
	Matrix<X_DIM, X_DIM> SqrtSigma;
	Matrix<G_DIM,X_DIM> J;

	for(int t = 0; t < T-1; ++t) {
		unVec(B[t], x, SqrtSigma);
		linearizeg(x, J);
		cost += alpha_belief*tr(J*SqrtSigma*SqrtSigma*~J) + alpha_control*tr(~U[t]*U[t]);
	}
	unVec(B[T-1], x, SqrtSigma);
	linearizeg(x, J);
	cost += alpha_final_belief*tr(J*SqrtSigma*SqrtSigma*~J);
	return cost;
}

void setupBeliefVars(beliefPenaltyMPC_params& problem, beliefPenaltyMPC_output& output)
{
	// problem inputs
	Q = new beliefPenaltyMPC_FLOAT*[T];
	f = new beliefPenaltyMPC_FLOAT*[T];
	lb = new beliefPenaltyMPC_FLOAT*[T];
	ub = new beliefPenaltyMPC_FLOAT*[T];
	C = new beliefPenaltyMPC_FLOAT*[T-1];
	e = new beliefPenaltyMPC_FLOAT*[T];

	// problem outputs
	z = new beliefPenaltyMPC_FLOAT*[T];

#define SET_VARS(n)    \
		Q[ BOOST_PP_SUB(n,1) ] = problem.Q##n ;  \
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
		Q[ BOOST_PP_SUB(n,1) ] = problem.Q##n ;  \
		f[ BOOST_PP_SUB(n,1) ] = problem.f##n ;  \
		lb[ BOOST_PP_SUB(n,1) ] = problem.lb##n ;	\
		ub[ BOOST_PP_SUB(n,1) ] = problem.ub##n ;	\
		e[ BOOST_PP_SUB(n,1) ] = problem.e##n ;  \
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
	delete[] C;
	delete[] e;
	delete[] z;
}

double computeMerit(const std::vector< Matrix<B_DIM> >& B, const std::vector< Matrix<U_DIM> >& U, double penalty_coeff)
{
	double merit = 0;
	Matrix<X_DIM> x;
	Matrix<X_DIM, X_DIM> SqrtSigma;
	Matrix<G_DIM,X_DIM> J;

	Matrix<B_DIM> beliefdynviol;
	for(int t = 0; t < T-1; ++t) {
		unVec(B[t], x, SqrtSigma);
		linearizeg(x, J);

		merit += alpha_belief*tr(J*SqrtSigma*SqrtSigma*~J) + alpha_control*tr(~U[t]*U[t]);

		beliefdynviol = (B[t+1] - beliefDynamics(B[t], U[t]) );
		for(int i = 0; i < B_DIM; ++i) {
			merit += penalty_coeff*fabs(beliefdynviol[i]);
		}
	}
	unVec(B[T-1], x, SqrtSigma);
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

// TODO: Check if all inputs are valid, Q, f, lb, ub, C, e, D at last time step
bool isValidInputs()
{
	for(int t = 0; t < T-1; ++t) {

		// check if Q, f, lb, ub, C, e, b are valid!

		//std::cout << std::endl << std::endl;
	}

	for(int i = 0; i < 33; ++i) {
		std::cout << lb[T-1][i] << " ";
	}
	std::cout << "\n\n";

	for(int i = 0; i < 27; ++i) {
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

bool minimizeMeritFunction(std::vector< Matrix<B_DIM> >& B, std::vector< Matrix<U_DIM> >& U, beliefPenaltyMPC_params& problem, beliefPenaltyMPC_output& output, beliefPenaltyMPC_info& info, double penalty_coeff, double trust_box_size)
{
	//LOG_DEBUG("Solving sqp problem with penalty parameter: %2.4f", penalty_coeff);
	std::cout << "Solving sqp problem with penalty parameter: " << penalty_coeff << std::endl;

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

	Matrix<B_DIM,B_DIM> IB = identity<B_DIM>();
	Matrix<B_DIM,B_DIM> minusIB = IB;
	for(int i = 0; i < B_DIM; ++i) {
		minusIB(i,i) = -1;
	}
	
	Matrix<3*B_DIM+U_DIM, 3*B_DIM+U_DIM> QMat;
	Matrix<B_DIM+2*G_DIM,B_DIM+2*G_DIM> QfMat;
	Matrix<B_DIM,3*B_DIM+U_DIM> CMat;
	Matrix<B_DIM> eVec;
	Matrix<S_DIM,S_DIM> Hess;
	Matrix<2*G_DIM,B_DIM+2*G_DIM> AMat;
	Matrix<2*G_DIM,1> bVec;

	// sqp loop
	while(true)
	{
		// In this loop, we repeatedly construct a linear approximation to the nonlinear belief dynamics constraint
		//LOG_DEBUG("  sqp iter: %d", sqp_iter);
		std::cout << "  sqp iter: " << sqp_iter << std::endl;

		merit = computeMerit(B, U, penalty_coeff);
		
		//LOG_DEBUG("  merit: %4.10f", merit);
		std::cout << "  merit: " << merit << std::endl;

		// Problem linearization and definition
		// fill in Q, f, C, e
		
		for (int t = 0; t < T-1; ++t) 
		{
			Matrix<B_DIM>& bt = B[t];
			Matrix<U_DIM>& ut = U[t];

			linearizeBeliefDynamics(bt, ut, F[t], G[t], h[t]);

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

			CMat.reset();
			eVec.reset();

			CMat.insert<B_DIM,B_DIM>(0,0,F[t]);
			CMat.insert<B_DIM,U_DIM>(0,B_DIM,G[t]);
			CMat.insert<B_DIM,B_DIM>(0,B_DIM+U_DIM,IB);
			CMat.insert<B_DIM,B_DIM>(0,2*B_DIM+U_DIM,minusIB);

			fillColMajor(C[t], CMat);

			if (t == 0) {
				eVec.insert<B_DIM,1>(0,0,B[0]);
				//std::cout << "eVec: " << ~eVec << std::endl;
				fillCol(e[0], eVec);
			} 
			
			eVec = -h[t] + F[t]*bt + G[t]*ut;
			fillCol(e[t+1], eVec);
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
		Matrix<X_DIM> xT;
		Matrix<X_DIM,X_DIM> SqrtSigmaT;

		unVec(B[T-1], xT, SqrtSigmaT);

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
			//LOG_DEBUG("       trust region size: %2.6f %2.6f", Beps, Ueps);
			std::cout << "       trust region size: " << Beps << ", " << Ueps << std::endl;

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
			// xGoal lower bound
			for(int i = 0; i < X_DIM; ++i) { lb[T-1][index++] = MAX(xMin[i], bT[i] - Beps); }
			// sigma lower bound
			for(int i = 0; i < S_DIM; ++i) { lb[T-1][index] = bT[index] - Beps; index++;}
			
			// for lower bound on L1 slacks
			for(int i = 0; i < 2*G_DIM; ++i) { lb[T-1][index++] = 0; }

			index = 0;
			// xGoal upper bound
			for(int i = 0; i < X_DIM; ++i) { ub[T-1][index++] = MIN(xMax[i], bT[i] + Beps); }
			// sigma lower bound
			for(int i = 0; i < S_DIM; ++i) { ub[T-1][index] = bT[index] + Beps; index++;}

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
				Beps *= cfg::trust_shrink_ratio;
				Ueps *= cfg::trust_shrink_ratio;
				//LOG_DEBUG("Shrinking trust region size to: %2.6f %2.6f", Beps, Ueps);
				std::cout << "  Shrinking trust region size to: " << Beps << ", " << Ueps << std::endl;
			} else {
				Beps *= cfg::trust_expand_ratio;
				Ueps *= cfg::trust_expand_ratio;
				B = Bopt; U = Uopt;
				prevcost = optcost;
				//LOG_DEBUG("Accepted, Increasing trust region size to:  %2.6f %2.6f", Beps, Ueps);
				std::cout << "  Accepted, Increasing trust region size to:  " << Beps << ", " << Ueps << std::endl;
				break;
			}

			if (Beps < cfg::min_trust_box_size && Ueps < cfg::min_trust_box_size) {
			    //LOG_DEBUG("Converged: x tolerance");
				std::cout << "  Converged: x tolerance" << std::endl;
			    return true;
			}

		} // trust region loop
		sqp_iter++;
	} // sqp loop

	return success;
}

double beliefPenaltyCollocation(std::vector< Matrix<B_DIM> >& B, std::vector< Matrix<U_DIM> >& U, beliefPenaltyMPC_params& problem, beliefPenaltyMPC_output& output, beliefPenaltyMPC_info& info)
{
	double costTime = 0;

	double penalty_coeff = cfg::initial_penalty_coeff;
	double trust_box_size = cfg::initial_trust_box_size;

	int penalty_increases = 0;

	Matrix<B_DIM> beliefdynviol;
	Matrix<2*G_DIM> goalposviol;

	// penalty loop
	while(penalty_increases < cfg::max_penalty_coeff_increases)
	{
		bool success = minimizeMeritFunction(B, U, problem, output, info, penalty_coeff, trust_box_size);

		double cntviol = 0;
		for(int t = 0; t < T-1; ++t) {
			beliefdynviol = (B[t+1] - beliefDynamics(B[t], U[t]) );
			for(int i = 0; i < B_DIM; ++i) {
				cntviol += fabs(beliefdynviol[i]);
			}
		}
		Matrix<X_DIM> xT;
		Matrix<X_DIM, X_DIM> SqrtSigmaT;
		unVec(B[T-1], xT, SqrtSigmaT);

		Matrix<G_DIM> delta;
		delta[0] = delta[1] = delta[2] = goaldelta;

		goalposviol.insert<G_DIM,1>(0,0,g(xT) - posGoal - delta);
		goalposviol.insert<G_DIM,1>(G_DIM,0, -g(xT) + posGoal - delta);

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
	    	return computeCost(B, U);
	    }
	}
	return computeCost(B, U);
}

bool testInitializationFeasibility(std::vector<Matrix<B_DIM> >& B, std::vector<Matrix<U_DIM> >& U)
{
	std::cout << "X initial" << std::endl;
	for (int t = 0; t < T; ++t) { 
		Matrix<U_DIM>& xt = B[t].subMatrix<X_DIM,1>(0,0);
		for (int i = 0; i < X_DIM; ++i) {
			if (xt[i] > xMax[i] || xt[i] < xMin[i]) {
				std::cerr << "Joint angle limit violated at joint " << i << " and time " << t << std::endl;
				return false;
			}
		}

		std::cout << std::setprecision(8) << ~xt; 
	}

	std::cout << "U initial" << std::endl;
	for (int t = 0; t < T-1; ++t) { 
		Matrix<U_DIM>& ut = U[t];
		for (int i = 0; i < U_DIM; ++i) {
			if (ut[i] > uMax[i] || ut[i] < uMin[i]) {
				std::cerr << "Control limit violated at joint " << i << " and time " << t << std::endl;
				return false;
			}
		}
		std::cout << std::setprecision(8) << ~ut; 
	}

	return true;
}

int main(int argc, char* argv[])
{
	
	initProblemParams();

	std::cout << "init problem params" << std::endl;

	Matrix<U_DIM> uinit = (xGoal - x0) / (double)((T-1)*DT);
	std::vector<Matrix<U_DIM> > U(T-1, uinit);

	std::vector<Matrix<B_DIM> > B(T);

	vec(x0, SqrtSigma0, B[0]);
	for (size_t t = 0; t < T-1; ++t) {
		B[t+1] = beliefDynamics(B[t], U[t]);
		//std::cout << ~B[t] << std::endl;
	}

	bool feasible = testInitializationFeasibility(B, U);
	if (!feasible) {
		std::cerr << "Infeasible trajectory initialization detected" << std::endl;
		std::exit(-1);
	}

	double initTrajCost = computeCost(B, U);
	std::cout << "Initial trajectory cost: " << initTrajCost << std::endl;

	// Create Obstacles in Callisto
	CAL_Initialisation (true,true,true);
	initEnvironment();

	std::vector<Matrix<X_DIM> > X;
	
	//readTrajectoryFromFile("data\\trajectory.txt", X);

	for(int t = 0; t < T; ++t) {
		X.push_back(B[t].subMatrix<X_DIM,1>(0,0));
	}

	double initLQGMPcost = computeLQGMPcost(X, U);
	std::cout << "Initial trajectory LQG-MP cost: " << initLQGMPcost << std::endl;

	//displayTrajectory(X, U, B, false);

	beliefPenaltyMPC_params problem;
	beliefPenaltyMPC_output output;
	beliefPenaltyMPC_info info;

	setupBeliefVars(problem, output);

	Timer solveTimer;
	//Timer_tic(&solveTimer);
	solveTimer.start();

	double cost = beliefPenaltyCollocation(B, U, problem, output, info);

	solveTimer.stop();
	//double solvetime = util::Timer_toc(&solveTimer);

	//LOG_INFO("Optimized cost: %4.10f", cost);
	//LOG_INFO("Solve time: %5.3f ms", solvetime*1000);
	
	std::cout << "Optimized cost: " << cost << std::endl;
	std::cout << "Actual cost: " << computeCost(B, U) << std::endl;
	std::cout << "Optimization time: " << solveTimer.interval_mS() << " mS" << std::endl;
	
	cleanupBeliefMPCVars();
	
	double finalLQGMPcost = computeLQGMPcost(X, U);
	std::cout << "Final trajectory LQG-MP cost: " << finalLQGMPcost << std::endl;

	//saveOptimizedTrajectory(U);
	//readOptimizedTrajectory(U);
	
	//displayTrajectory(X, U, B);
	/*
	vec(x0, SqrtSigma0, B[0]);
	for (size_t t = 0; t < T-1; ++t) {
		B[t+1] = beliefDynamics(B[t], U[t]);
	}
	*/

	X.clear();
	for(int t = 0; t < T; ++t) {
		X.push_back(B[t].subMatrix<X_DIM,1>(0,0));
	}
	
	displayTrajectory(X, U, B);

	int k;
	std::cin >> k;

	CAL_End();

	return 0;
}
