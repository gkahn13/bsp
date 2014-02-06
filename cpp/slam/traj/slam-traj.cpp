#include <vector>
#include <iomanip>
#include <sys/resource.h>

#include "util/matrix.h"
#include "util/Timer.h"
#include "util/logging.h"

extern "C" {
#include "trajMPC.h"
trajMPC_FLOAT **f, **lb, **ub, **C, **e, **z;
}

#include "boost/preprocessor.hpp"

#include "../slam.h"

Matrix<C_DIM> c0;
Matrix<C_DIM> cGoal;

const double alpha_goal = 10;

namespace cfg {
const double improve_ratio_threshold = .1;
const double min_approx_improve = 1e-3;
const double min_trust_box_size = 1e-3;
const double trust_shrink_ratio = .1;
const double trust_expand_ratio = 1.5;
const double cnt_tolerance = 1e-4;
const double penalty_coeff_increase_ratio = 10;
const double initial_penalty_coeff = 25;
const double initial_trust_box_size = 1;
const int max_penalty_coeff_increases = 2;
const int max_sqp_iterations = 50;
}

Matrix<C_DIM> dynfunccar(const Matrix<C_DIM>& x, const Matrix<U_DIM>& u)
{
	Matrix<C_DIM> xAdd = zeros<C_DIM,1>();

	xAdd[0] = u[0] * DT * cos(x[2]+u[1]);
	xAdd[1] = u[0] * DT * sin(x[2]+u[1]);
	xAdd[2] = u[0] * DT * sin(u[1])/config::WHEELBASE;

	Matrix<C_DIM> xNew = x + xAdd;
    return xNew;
}

// Jacobians: df(x,u)/dx, df(x,u)/du
void linearizeCarDynamics(const Matrix<C_DIM>& x, const Matrix<U_DIM>& u, Matrix<C_DIM,C_DIM>& F, Matrix<C_DIM,U_DIM>& G, Matrix<C_DIM>& h)
{
	F.reset();
	Matrix<C_DIM> xr(x), xl(x);
	for (size_t i = 0; i < C_DIM; ++i) {
		xr[i] += step; xl[i] -= step;
		F.insert(0,i, (dynfunccar(xr, u) - dynfunccar(xl, u)) / (xr[i] - xl[i]));
		xr[i] = x[i]; xl[i] = x[i];
	}

	G.reset();
	Matrix<U_DIM> ur(u), ul(u);
	for (size_t i = 0; i < U_DIM; ++i) {
		ur[i] += step; ul[i] -= step;
		G.insert(0,i, (dynfunccar(x, ur) - dynfunccar(x, ul)) / (ur[i] - ul[i]));
		ur[i] = u[i]; ul[i] = u[i];
	}

	h = dynfunccar(x, u);
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



double computeCost(const std::vector< Matrix<C_DIM> >& X, const std::vector< Matrix<U_DIM> >& U)
{
	double cost = 0;
	for(int t = 0; t < T-1; ++t) {
		cost += alpha_control*tr(~U[t]*U[t]);
	}
	Matrix<P_DIM> pT = X[T-1].subMatrix<P_DIM,1>(0, 0);
	Matrix<P_DIM> pGoal = cGoal.subMatrix<P_DIM,1>(0, 0);
	cost += alpha_goal*tr(~(pT-pGoal)*(pT-pGoal));
	return cost;
}


void setupTrajVars(trajMPC_params &problem, trajMPC_output &output)
{
	// problem inputs
	f = new trajMPC_FLOAT*[T];
	lb = new trajMPC_FLOAT*[T];
	ub = new trajMPC_FLOAT*[T];
	C = new trajMPC_FLOAT*[T-1];
	e = new trajMPC_FLOAT*[T];

	// problem outputs
	z = new trajMPC_FLOAT*[T];

	// H[ BOOST_PP_SUB(n,1) ] = problem.H##n ;
#define SET_VARS(n)    \
		f[ BOOST_PP_SUB(n,1) ] = problem.f##n ;  \
		lb[ BOOST_PP_SUB(n,1) ] = problem.lb##n ;	\
		ub[ BOOST_PP_SUB(n,1) ] = problem.ub##n ;	\
		C[ BOOST_PP_SUB(n,1) ] = problem.C##n ;  \
		e[ BOOST_PP_SUB(n,1) ] = problem.e##n ;  \
		z[ BOOST_PP_SUB(n,1) ] = output.z##n ;
#define BOOST_PP_LOCAL_MACRO(n) SET_VARS(n)
#define BOOST_PP_LOCAL_LIMITS (1, TIMESTEPS-1)
#include BOOST_PP_LOCAL_ITERATE()

	// H[ BOOST_PP_SUB(n,1) ] = problem.H##n ;
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

void cleanupTrajMPCVars()
{
	delete[] f;
	delete[] lb;
	delete[] ub;
	delete[] C;
	delete[] e;
	delete[] z;
}

double computeMerit(const std::vector< Matrix<C_DIM> >& X, const std::vector< Matrix<U_DIM> >& U, double penalty_coeff)
{
	double merit = 0;
	Matrix<C_DIM> dynviol;
	for(int t = 0; t < T-1; ++t) {
		merit += alpha_control*tr(~U[t]*U[t]);
		dynviol = (X[t+1] - dynfunccar(X[t], U[t]) );
		for(int i = 0; i < C_DIM; ++i) {
			merit += penalty_coeff*fabs(dynviol[i]);
		}
	}
	Matrix<P_DIM> pT = X[T-1].subMatrix<P_DIM,1>(0, 0);
	Matrix<P_DIM> pGoal = cGoal.subMatrix<P_DIM,1>(0, 0);
	merit += alpha_goal*tr(~(pT-pGoal)*(pT-pGoal));
	return merit;
}

// TODO: Check if all inputs are valid, H, f, lb, ub, C, e, D at last time step
bool isValidInputs()
{
	for(int t = 0; t < T-1; ++t) {

		std::cout << "t: " << t << std::endl << std::endl;

		/*
		std::cout << "f: ";
		for(int i = 0; i < (2*C_DIM+U_DIM); ++i) {
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

		/*
		std::cout << "lb s, t: ";
		for(int i = 0; i < 2*C_DIM; ++i) {
			std::cout << lb[t][C_DIM+U_DIM+i] << " ";
		}
		std::cout << std::endl;
		*/

		std::cout << "lb u: ";
		for(int i = 0; i < U_DIM; ++i) {
			std::cout << lb[t][C_DIM+i] << " ";
		}
		std::cout << std::endl;

		std::cout << "ub u: ";
		for(int i = 0; i < U_DIM; ++i) {
			std::cout << ub[t][C_DIM+i] << " ";
		}
		std::cout << std::endl;

		/*

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
		for(int i = 0; i < C_DIM; ++i) {
			std::cout << e[t][i] << " ";
		}
		*/
		std::cout << std::endl << std::endl;

	}

	/*
	std::cout << "e:" << std::endl;
	for(int i = 0; i < C_DIM; ++i) {
		std::cout << e[T-1][i] << " ";
	}
	std::cout << std::endl << std::endl;
	*/

	std::cout << "lb c: ";
	for(int i = 0; i < C_DIM; ++i) {
		std::cout << lb[T-1][i] << " ";
	}
	std::cout << std::endl;

	std::cout << "ub c: ";
	for(int i = 0; i < C_DIM; ++i) {
		std::cout << ub[T-1][i] << " ";
	}
	std::cout << std::endl;
	std::cout << std::endl;
	return true;
}


bool minimizeMeritFunction(std::vector< Matrix<C_DIM> >& X, std::vector< Matrix<U_DIM> >& U, trajMPC_params &problem, trajMPC_output &output, trajMPC_info &info, double penalty_coeff, double trust_box_size)
{
	LOG_DEBUG("Solving sqp problem with penalty parameter: %2.4f", penalty_coeff);

	//Matrix<B_DIM,1> b0 = B[0];

	std::vector< Matrix<C_DIM,C_DIM> > F(T-1);
	std::vector< Matrix<C_DIM,U_DIM> > G(T-1);
	std::vector< Matrix<C_DIM> > h(T-1);

	double Xeps = trust_box_size;
	double Ueps = trust_box_size;

	double optcost;

	std::vector<Matrix<C_DIM> > Xopt(T);
	std::vector<Matrix<U_DIM> > Uopt(T-1);

	double merit, model_merit, new_merit;
	double approx_merit_improve, exact_merit_improve, merit_improve_ratio;

	int sqp_iter = 1, index = 0;
	bool success;

	Matrix<C_DIM,C_DIM> IB = identity<C_DIM>();
	Matrix<C_DIM,C_DIM> minusIB = -identity<C_DIM>();
	
	Matrix<C_DIM,3*C_DIM+U_DIM> CMat;
	Matrix<C_DIM> eVec;

	// sqp loop
	while(true)
	{
		// In this loop, we repeatedly construct a linear approximation to the nonlinear belief dynamics constraint
		LOG_DEBUG("  sqp iter: %d", sqp_iter);

		merit = computeMerit(X, U, penalty_coeff);
		
		LOG_DEBUG("  merit: %4.10f", merit);

		// Problem linearization and definition
		// fill in f, C, e
		
		for (int t = 0; t < T-1; ++t) 
		{
			Matrix<C_DIM>& xt = X[t];
			Matrix<U_DIM>& ut = U[t];

			linearizeCarDynamics(xt, ut, F[t], G[t], h[t]);

			//std::cout << "h[" << t << "]" << std::endl << h[t] << std::endl;
			//std::cout << "F[" << t << "]" << std::endl << F[t] << std::endl;
			//std::cout << "xt[" << t << "]" << ~xt;
			//std::cout << "G[" << t << "]" << std::endl << G[t] << std::endl;
			//std::cout << "ut[" << t << "]" << ~ut;

			// initialize f in cost function to penalize
			// belief dynamics slack variables
			index = 0;
			for(int i = 0; i < (C_DIM+U_DIM); ++i) { f[t][index++] = 0; }
			for(int i = 0; i < 2*C_DIM; ++i) { f[t][index++] = penalty_coeff; }

			CMat.reset();
			eVec.reset();

			CMat.insert<C_DIM,C_DIM>(0,0,F[t]);
			CMat.insert<C_DIM,U_DIM>(0,C_DIM,G[t]);
			CMat.insert<C_DIM,C_DIM>(0,C_DIM+U_DIM,IB);
			CMat.insert<C_DIM,C_DIM>(0,2*C_DIM+U_DIM,minusIB);

			fillColMajor(C[t], CMat);

			if (t == 0) {
				eVec.insert<C_DIM,1>(0,0,X[0]);
				fillCol(e[0], eVec);
			} 
			
			eVec = -h[t] + F[t]*xt + G[t]*ut;
			//std::cout << "eVec " << ~eVec << std::endl;
			fillCol(e[t+1], eVec);

		}
		
		for(int i=0; i < P_DIM; ++i) { f[T-1][i] = -2*alpha_goal*cGoal[i]; }
		f[T-1][2] = 0;

		//std::cout << "PAUSED INSIDE MINIMIZEMERITFUNCTION" << std::endl;
		//int k;
		//std::cin >> k;


		// trust region size adjustment
		while(true)
		{
			LOG_DEBUG("       trust region size: %2.6f %2.6f", Xeps, Ueps);

			util::Timer fillVarsTimer;
			util::Timer_tic(&fillVarsTimer);
			// solve the innermost QP here
			for(int t = 0; t < T-1; ++t)
			{
				Matrix<C_DIM>& xt = X[t];
				Matrix<U_DIM>& ut = U[t];

				// Fill in lb, ub

				index = 0;
				// x lower bound
				for(int i = 0; i < C_DIM; ++i) { lb[t][index++] = MAX(xMin[i], xt[i] - Xeps); }
				// u lower bound
				for(int i = 0; i < U_DIM; ++i) { lb[t][index++] = MAX(uMin[i], ut[i] - Ueps); }

				// for lower bound on L1 slacks
				for(int i = 0; i < 2*C_DIM; ++i) { lb[t][index++] = 0; }

				index = 0;
				// x upper bound
				for(int i = 0; i < C_DIM; ++i) { ub[t][index++] = MIN(xMax[i], xt[i] + Xeps); }
				// u upper bound
				for(int i = 0; i < U_DIM; ++i) { ub[t][index++] = MIN(uMax[i], ut[i] + Ueps); }

				//for(int i = 0; i < 2*B_DIM; ++i) { ub[t][index++] = INFTY; }
			}

			Matrix<C_DIM>& xT = X[T-1];

			// Fill in lb, ub, C, e
			index = 0;
			double delta = .5;
			// cGoal lower bound
			for(int i = 0; i < P_DIM; ++i) { lb[T-1][index++] = cGoal[i] - delta; }
			lb[T-1][index++] = MIN(xMin[2], xT[2] - Xeps); // none on angles

			index = 0;
			// cGoal upper bound
			for(int i = 0; i < P_DIM; ++i) { ub[T-1][index++] = cGoal[i] + delta; }
			ub[T-1][index++] = MAX(xMax[2], xT[2] + Xeps);

			double jac_constant = 0, hessian_constant = 0, constant_cost = 0;
			//for(int i = 0; i < C_DIM; ++i) { jac_constant += 0.5*(f[T-1][i]*cGoal[i]); }
			//jac_constant = alpha_goal*tr(~cGoal*cGoal);
			//hessian_constant = alpha_goal*tr(~xT*xT);

			Matrix<P_DIM> pGoal = cGoal.subMatrix<P_DIM,1>(0, 0);
			constant_cost = alpha_goal*tr(~pGoal*pGoal);

			//LOG_DEBUG("Hessian constant: %10.4f", hessian_constant);
			//LOG_DEBUG("Jacobian constant: %10.4f", jac_constant);
			LOG_DEBUG("Constant cost: %10.4f", constant_cost);

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

			int exitflag = trajMPC_solve(&problem, &output, &info);
			if (exitflag == 1) {
				for(int t = 0; t < T-1; ++t) {
					Matrix<C_DIM>& xt = Xopt[t];
					Matrix<U_DIM>& ut = Uopt[t];

					for(int i = 0; i < C_DIM; ++i) {
						xt[i] = z[t][i];
					}
					for(int i = 0; i < U_DIM; ++i) {
						ut[i] = z[t][C_DIM+i];
					}
					optcost = info.pobj;
				}
				for(int i = 0; i < C_DIM; ++i) {
					Xopt[T-1][i] = z[T-1][i];
				}
			}
			else {
				LOG_ERROR("Some problem in solver");
				exit(-1);
			}

			LOG_DEBUG("Optimized cost: %4.10f", optcost);

			model_merit = optcost + constant_cost;
			new_merit = computeMerit(Xopt, Uopt, penalty_coeff);

			//std::cout << "Xopt" << std::endl;
			//for(int t=0; t < T; ++t) { std::cout << ~Xopt[t]; }
			//std::cout << std::endl << std::endl;

			//std::cout << "Uopt" << std::endl;
			//for(int t=0; t < T-1; ++t) { std::cout << ~Uopt[t]; }
			//std::cout << std::endl << std::endl;

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
			//std::cin.ignore();

			if (approx_merit_improve < -1e-3) {
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
				LOG_DEBUG("Shrinking trust region size to: %2.6f %2.6f", Xeps, Ueps);
			} else {
				Xeps *= cfg::trust_expand_ratio;
				Ueps *= cfg::trust_expand_ratio;
				X = Xopt; U = Uopt;
				LOG_DEBUG("Accepted, Increasing trust region size to:  %2.6f %2.6f", Xeps, Ueps);
				break;
			}

			if (Xeps < cfg::min_trust_box_size && Ueps < cfg::min_trust_box_size) {
			    LOG_DEBUG("Converged: x tolerance");
			    return true;
			}

		} // trust region loop
		sqp_iter++;
	} // sqp loop

	return success;
}

double trajCollocation(std::vector< Matrix<C_DIM> >& X, std::vector< Matrix<U_DIM> >& U, trajMPC_params &problem, trajMPC_output &output, trajMPC_info &info)
{
	double penalty_coeff = cfg::initial_penalty_coeff;
	double trust_box_size = cfg::initial_trust_box_size;

	int penalty_increases = 0;

	Matrix<C_DIM> dynviol;

	// penalty loop
	while(penalty_increases < cfg::max_penalty_coeff_increases)
	{
		bool success = minimizeMeritFunction(X, U, problem, output, info, penalty_coeff, trust_box_size);

		double cntviol = 0;
		for(int t = 0; t < T-1; ++t) {
			dynviol = (X[t+1] - dynfunccar(X[t], U[t]) );
			for(int i = 0; i < C_DIM; ++i) {
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


void initTraj(const Matrix<C_DIM>& cStart, const Matrix<C_DIM>& cEnd, std::vector<Matrix<C_DIM> >& X, std::vector<Matrix<U_DIM> >& U) {
	trajMPC_params problem;
	trajMPC_output output;
	trajMPC_info info;
	setupTrajVars(problem, output);

	c0 = cStart;
	cGoal = cEnd;


	Matrix<U_DIM> uinit;
	uinit[0] = sqrt((c0[0] - cGoal[0])*(c0[0] - cGoal[0]) + (c0[1] - cGoal[1])*(c0[1] - cGoal[1])) / (double)((T-1)*DT);
	uinit[1] = 0;

	// TEMP
	uinit[0] *= .1;

	std::cout << "uinit: " << ~uinit << std::endl;

	for(int i=0; i < T-1; ++i) { U[i] = uinit; }

	std::cout << "X initial" << std::endl;
	X[0].insert(0,0,c0);
	for(int t=0; t < T-1; ++t) {
		std::cout << ~X[t];
		X[t+1] = dynfunccar(X[t],U[t]);
	}
	std::cout << ~X[T-1];

	double initTrajCost = computeCost(X, U);
	LOG_INFO("Initial trajectory cost: %4.10f", initTrajCost);

	util::Timer solveTimer;
	Timer_tic(&solveTimer);

	double cost = trajCollocation(X, U, problem, output, info);

	double solvetime = util::Timer_toc(&solveTimer);

	std::cout << "Final input" << std::endl;
	for(int t = 0; t < T-1; ++t) {
		std::cout << ~U[t];
	}
	std::cout << std::endl;

	std::cout << "Final trajectory" << std::endl;
	for (size_t t = 0; t < T-1; ++t) {
		std::cout << ~X[t];
		X[t+1] = dynfunccar(X[t], U[t]);
	}
	std::cout << ~X[T-1];

	LOG_INFO("Initial cost: %4.10f", initTrajCost);
	LOG_INFO("Optimized cost: %4.10f", cost);
	LOG_INFO("Actual cost: %4.10f", computeCost(X,U));
	LOG_INFO("Solve time: %5.3f ms", solvetime*1000);
}


int main(int argc, char* argv[])
{
	initProblemParams();

	Matrix<C_DIM> cStart, cEnd;
	std::vector<Matrix<C_DIM> > X(T, zeros<C_DIM,1>());
	std::vector<Matrix<U_DIM> > U(T-1, zeros<U_DIM,1>());

	cStart[0] = 0;
	cStart[1] = 0;

	cEnd[0] = 60;
	cEnd[1] = 0;

	//cStart[2] = atan2(cEnd[1] - cStart[1], cEnd[0] - cStart[0]);
	cStart[2] = .5*M_PI;
	cEnd[2] = atan2(cEnd[1] - cStart[1], cEnd[0] - cStart[0]);

	initTraj(cStart, cEnd, X, U);

	cleanupTrajMPCVars();
	
}
