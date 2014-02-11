#include "../slam.h"
#include "../traj/slam-traj.h"


#include <vector>
#include <iomanip>

#include "util/matrix.h"
#include "util/Timer.h"

extern "C" {
#include "stateMPC.h"
stateMPC_FLOAT **H, **f, **lb, **ub, **C, **e, **z;
#include "slam-state-casadi.h"
}

#include "boost/preprocessor.hpp"

const double alpha_belief = 10;//10;
const double alpha_final_belief = 50;//50;
const double alpha_control = .01;//.01

namespace cfg {
const double improve_ratio_threshold = .1; // .1
const double min_approx_improve = 1e-3; // 1e-4
const double min_trust_box_size = 1e-2; // 1e-3
const double trust_shrink_ratio = .5; // .1
const double trust_expand_ratio = 1.5; // 1.5
const double cnt_tolerance = 1e-4;
const double penalty_coeff_increase_ratio = 5; // 5
const double initial_penalty_coeff = 5; // 5
const double initial_trust_box_size = 1; // 5 // split up trust box size for X and U
const int max_penalty_coeff_increases = 3; // 3
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


void setupCasadiVars(const std::vector<Matrix<X_DIM> >& X, const std::vector<Matrix<U_DIM> >& U, double* XU_arr, double* Sigma0_arr, double* params_arr)
{
	int index = 0;
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
	index = 0;
	for(int i=0; i < X_DIM; ++i) {
		for(int j=0; j < X_DIM; ++j) {
			Sigma0_arr[index++] = Sigma0(i,j);
		}
	}

	params_arr[0] = alpha_belief;
	params_arr[1] = alpha_control;
	params_arr[2] = alpha_final_belief;

}

double casadiComputeCost(const std::vector< Matrix<X_DIM> >& X, const std::vector< Matrix<U_DIM> >& U)
{
	double XU_arr[XU_DIM];
	double Sigma0_arr[X_DIM*X_DIM];
	double params_arr[3];

	setupCasadiVars(X, U, XU_arr, Sigma0_arr, params_arr);

	const double **casadi_input = new const double*[3];
	casadi_input[0] = XU_arr;
	casadi_input[1] = Sigma0_arr;
	casadi_input[2] = params_arr;

	double cost = 0;
	double **cost_arr = new double*[1];
	cost_arr[0] = &cost;

	evaluateCostWrap(casadi_input, cost_arr);

	return cost;
}

double casadiComputeMerit(const std::vector< Matrix<X_DIM> >& X, const std::vector< Matrix<U_DIM> >& U, double penalty_coeff)
{
	double merit = 0;

	merit = casadiComputeCost(X, U);

	Matrix<X_DIM> x, dynviol;
	Matrix<X_DIM, X_DIM> SqrtSigma;
	Matrix<B_DIM> b, b_tp1;
	vec(x0, SqrtSigma0, b);

	for(int t = 0; t < T-1; ++t) {
		unVec(b, x, SqrtSigma);
		b_tp1 = beliefDynamics(b, U[t]);
		dynviol = (X[t+1] - b_tp1.subMatrix<X_DIM,1>(0,0) );
		for(int i = 0; i < X_DIM; ++i) {
			merit += penalty_coeff*fabs(dynviol[i]);
		}
		b = b_tp1;
	}
	return merit;
}

void casadiComputeCostGrad(const std::vector< Matrix<X_DIM> >& X, const std::vector< Matrix<U_DIM> >& U, double& cost, Matrix<XU_DIM>& Grad)
{
	double XU_arr[XU_DIM];
	double Sigma0_arr[X_DIM*X_DIM];
	double params_arr[3];

	setupCasadiVars(X, U, XU_arr, Sigma0_arr, params_arr);

	const double **casadi_input = new const double*[3];
	casadi_input[0] = XU_arr;
	casadi_input[1] = Sigma0_arr;
	casadi_input[2] = params_arr;

	double **costgrad_arr = new double*[2];
	costgrad_arr[0] = &cost;
	costgrad_arr[1] = Grad.getPtr();

	evaluateCostGradWrap(casadi_input, costgrad_arr);

	/*
	double jac_cost = 0;
	for(int i = 0; i < XU_DIM; ++i) {
		jac_cost += Grad[i]*XU_arr[i];
	}
	LOG_DEBUG("jac cost: %4.10f",jac_cost);
	*/
}

double computeCost(const std::vector< Matrix<X_DIM> >& X, const std::vector< Matrix<U_DIM> >& U)
{
	double cost = 0;
	Matrix<B_DIM> b;
	Matrix<X_DIM> x;
	Matrix<X_DIM, X_DIM> SqrtSigma;
	vec(x0, SqrtSigma0, b);

	for(int t = 0; t < T-1; ++t) {
		unVec(b, x, SqrtSigma);
		cost += alpha_belief*tr(SqrtSigma*SqrtSigma);
		cost += alpha_control*tr(~U[t]*U[t]);
		b = beliefDynamics(b, U[t]);
	}
	unVec(b, x, SqrtSigma);
	cost += alpha_final_belief*tr(SqrtSigma*SqrtSigma);
	return cost;
}

double computeMerit(const std::vector< Matrix<X_DIM> >& X, const std::vector< Matrix<U_DIM> >& U, double penalty_coeff)
{
	double merit = 0;
	Matrix<X_DIM> x, dynviol;
	Matrix<X_DIM, X_DIM> SqrtSigma;
	Matrix<B_DIM> b, b_tp1;
	vec(x0, SqrtSigma0, b);

	for(int t = 0; t < T-1; ++t) {
		unVec(b, x, SqrtSigma);
		merit += alpha_belief*tr(SqrtSigma*SqrtSigma) + alpha_control*tr(~U[t]*U[t]);
		b_tp1 = beliefDynamics(b, U[t]);
		dynviol = (X[t+1] - b_tp1.subMatrix<X_DIM,1>(0,0) );
		for(int i = 0; i < X_DIM; ++i) {
			merit += penalty_coeff*fabs(dynviol[i]);
		}
		b = b_tp1;
	}
	unVec(b, x, SqrtSigma);
	merit += alpha_final_belief*tr(SqrtSigma*SqrtSigma);
	return merit;
}

// Jacobians: dg(b,u)/db, dg(b,u)/du
void linearizeCarDynamics(const Matrix<X_DIM>& x, const Matrix<U_DIM>& u, Matrix<X_DIM,X_DIM>& F, Matrix<X_DIM,U_DIM>& G, Matrix<X_DIM>& h)
{
	F.reset();
	Matrix<X_DIM> xr(x), xl(x);
	for (size_t i = 0; i < X_DIM; ++i) {
		xr[i] += step; xl[i] -= step;
		F.insert(0,i, (dynfunc(xr, u, zeros<Q_DIM,1>()) - dynfunc(xl, u, zeros<Q_DIM,1>())) / (xr[i] - xl[i]));
		xr[i] = x[i]; xl[i] = x[i];
	}

	G.reset();
	Matrix<U_DIM> ur(u), ul(u);
	for (size_t i = 0; i < U_DIM; ++i) {
		ur[i] += step; ul[i] -= step;
		G.insert(0,i, (dynfunc(x, ur, zeros<Q_DIM,1>()) - dynfunc(x, ul, zeros<Q_DIM,1>())) / (ur[i] - ul[i]));
		ur[i] = u[i]; ul[i] = u[i];
	}

	h = dynfunc(x, u, zeros<Q_DIM,1>());
}


void setupStateVars(stateMPC_params& problem, stateMPC_output& output)
{
	// problem inputs
	H = new stateMPC_FLOAT*[T];
	f = new stateMPC_FLOAT*[T];
	lb = new stateMPC_FLOAT*[T];
	ub = new stateMPC_FLOAT*[T];
	C = new stateMPC_FLOAT*[T-1];
	e = new stateMPC_FLOAT*[T];

	// problem outputs
	z = new stateMPC_FLOAT*[T];

#define SET_VARS(n)    \
		H[ BOOST_PP_SUB(n,1) ] = problem.H##n ;  \
		f[ BOOST_PP_SUB(n,1) ] = problem.f##n ;  \
		lb[ BOOST_PP_SUB(n,1) ] = problem.lb##n ;	\
		ub[ BOOST_PP_SUB(n,1) ] = problem.ub##n ;	\
		C[ BOOST_PP_SUB(n,1) ] = problem.C##n ;	\
		e[ BOOST_PP_SUB(n,1) ] = problem.e##n ;	\
		z[ BOOST_PP_SUB(n,1) ] = output.z##n ;

#define BOOST_PP_LOCAL_MACRO(n) SET_VARS(n)
#define BOOST_PP_LOCAL_LIMITS (1, TIMESTEPS-1)
#include BOOST_PP_LOCAL_ITERATE()

#define SET_LAST_VARS(n)    \
		H[ BOOST_PP_SUB(n,1) ] = problem.H##n ;  \
		f[ BOOST_PP_SUB(n,1) ] = problem.f##n ;  \
		lb[ BOOST_PP_SUB(n,1) ] = problem.lb##n ;	\
		ub[ BOOST_PP_SUB(n,1) ] = problem.ub##n ;	\
		e[ BOOST_PP_SUB(n,1) ] = problem.e##n ;	\
		z[ BOOST_PP_SUB(n,1) ] = output.z##n ;

#define BOOST_PP_LOCAL_MACRO(n) SET_LAST_VARS(n)
#define BOOST_PP_LOCAL_LIMITS (TIMESTEPS, TIMESTEPS)
#include BOOST_PP_LOCAL_ITERATE()

}

void cleanupStateMPCVars()
{
	delete[] H;
	delete[] f;
	delete[] lb;
	delete[] ub;
	delete[] C;
	delete[] e;
	delete[] z;
}

// TODO: Check if all inputs are valid, Q, f, lb, ub, A, b at last time step
bool isValidInputs()
{
	// check if Q, f, lb, ub, e are valid!
	for(int t = 0; t < T-1; ++t)
	{
		//for(int i = 0; i < 144; ++i) {
		//	std::cout << Q[t][i] << " ";
		//}
		for(int i = 0; i < 12; ++i) {
			std::cout << lb[t][i] << " ";
		}
		std::cout << std::endl;
		for(int i = 0; i < 12; ++i) {
			std::cout << ub[t][i] << " ";
		}
		std::cout << "\n\n";
	}
	for(int i = 0; i < 12; ++i) {
		std::cout << lb[T-1][i] << " ";
	}
	std::cout << std::endl;
	for(int i = 0; i < 6; ++i) {
		std::cout << ub[T-1][i] << " ";
	}
	std::cout << "\n\n";


	int magic;
	std::cin >> magic;

	return true;
}

bool minimizeMeritFunction(std::vector< Matrix<X_DIM> >& X, std::vector< Matrix<U_DIM> >& U, stateMPC_params& problem, stateMPC_output& output, stateMPC_info& info, double penalty_coeff, double trust_box_size)
{
	LOG_DEBUG("Solving sqp problem with penalty parameter: %2.4f", penalty_coeff);

	//Matrix<X_DIM,1> x0 = X[0];

	// constrain initial state
	//for(int i = 0; i < X_DIM; ++i) {
	//	e[i] = x0[i];
	//}

	std::vector< Matrix<X_DIM,X_DIM> > F(T-1);
	std::vector< Matrix<X_DIM,U_DIM> > G(T-1);
	std::vector< Matrix<X_DIM> > h(T-1);

	double Xeps = trust_box_size;
	double Ueps = trust_box_size;

	double optcost;

	std::vector<Matrix<X_DIM> > Xopt(T);
	std::vector<Matrix<U_DIM> > Uopt(T-1);

	double merit, model_merit, new_merit;
	double approx_merit_improve, exact_merit_improve, merit_improve_ratio;
	double constant_cost, hessian_constant, jac_constant;

	int sqp_iter = 1, index = 0;
	bool success;

	Matrix<X_DIM+U_DIM, X_DIM+U_DIM> HMat;
	Matrix<X_DIM,X_DIM> HfMat;
	Matrix<X_DIM> eVec;
	Matrix<X_DIM,3*X_DIM+U_DIM> CMat;

	Matrix<X_DIM,X_DIM> IX = identity<X_DIM>();
	Matrix<X_DIM,X_DIM> minusIX = IX;
	for(int i = 0; i < X_DIM; ++i) {
		minusIX(i,i) = -1;
	}

	Matrix<X_DIM+U_DIM> zbar;

	// full Hessian from current timstep
	Matrix<XU_DIM,XU_DIM> B = identity<XU_DIM>();

	Matrix<XU_DIM> Grad, Gradopt;
	double cost;
	int idx = 0;

	// sqp loop
	while(true)
	{
		// In this loop, we repeatedly construct a linear approximation to the nonlinear belief dynamics constraint
		LOG_DEBUG("  sqp iter: %d", sqp_iter);

		merit = casadiComputeMerit(X, U, penalty_coeff);

		LOG_DEBUG("  merit: %4.10f", merit);

		// Compute gradients
		casadiComputeCostGrad(X, U, cost, Grad);

		// Problem linearization and definition
		// fill in Q, f

		hessian_constant = 0;
		jac_constant = 0;
		idx = 0;

		for (int t = 0; t < T-1; ++t)
		{
			Matrix<X_DIM>& xt = X[t];
			Matrix<U_DIM>& ut = U[t];

			idx = t*(X_DIM+U_DIM);
			//LOG_DEBUG("idx: %d",idx);

			// fill in gradients and Hessians

			HMat.reset();
			for(int i = 0; i < (X_DIM+U_DIM); ++i) {
				double val = B(idx+i,idx+i);
				HMat(i,i) = (val < 0) ? 0 : val;
			}

			// since diagonal, fill directly
			for(int i = 0; i < (X_DIM+U_DIM); ++i) { H[t][i] = HMat(i,i); }
			for(int i = 0; i < (2*X_DIM); ++i) { H[t][i + (X_DIM+U_DIM)] = 0; }

			zbar.insert(0,0,xt);
			zbar.insert(X_DIM,0,ut);

			for(int i = 0; i < (X_DIM+U_DIM); ++i) {
				hessian_constant += HMat(i,i)*zbar[i]*zbar[i];
				jac_constant -= Grad[idx+i]*zbar[i];
				f[t][i] = Grad[idx+i] - HMat(i,i)*zbar[i];
			}

			// penalize dynamics slack variables
			for(int i = X_DIM+U_DIM; i < 3*X_DIM+U_DIM; ++i) { f[t][i] = penalty_coeff; }

			// fill in linearizations
			linearizeCarDynamics(xt, ut, F[t], G[t], h[t]);


			CMat.reset();
			eVec.reset();

			CMat.insert<X_DIM,X_DIM>(0,0,F[t]);
			CMat.insert<X_DIM,U_DIM>(0,X_DIM,G[t]);
			CMat.insert<X_DIM,X_DIM>(0,X_DIM+U_DIM,IX);
			CMat.insert<X_DIM,X_DIM>(0,2*X_DIM+U_DIM,minusIX);

			fillColMajor(C[t], CMat);

			if (t == 0) {
				eVec.insert<X_DIM,1>(0,0,X[0]);
				fillCol(e[0], eVec);
			}

			eVec = -h[t] + F[t]*xt + G[t]*ut;
			fillCol(e[t+1], eVec);
		}

		// For last stage, fill in Q, f, A, b
		Matrix<X_DIM>& xT = X[T-1];

		idx = (T-1)*(X_DIM+U_DIM);
		//LOG_DEBUG("idx: %d",idx);

		HfMat.reset();
		for(int i = 0; i < X_DIM; ++i) {
			double val = B(idx+i,idx+i);
			HfMat(i,i) = (val < 0) ? 0 : val;
		}

		// since diagonal, fill directly
		for(int i = 0; i < X_DIM; ++i) { H[T-1][i] = HMat(i,i); }

		for(int i = 0; i < X_DIM; ++i) {
			hessian_constant += HfMat(i,i)*xT[i]*xT[i];
			jac_constant -= Grad[idx+i]*xT[i];
			f[T-1][i] = Grad[idx+i] - HfMat(i,i)*xT[i];
		}


		constant_cost = 0.5*hessian_constant + jac_constant + cost;
		LOG_DEBUG("  hessian cost: %4.10f", 0.5*hessian_constant);
		LOG_DEBUG("  jacobian cost: %4.10f", jac_constant);
		LOG_DEBUG("  constant cost: %4.10f", constant_cost);

		//std::cout << "PAUSED INSIDE MINIMIZEMERITFUNCTION" << std::endl;
		//int k;
		//std::cin >> k;


		// trust region size adjustment
		while(true)
		{
			LOG_DEBUG("       trust region size: %2.6f %2.6f", Xeps, Ueps);
			//std::cout << "       trust region size: " << Xeps << ", " << Ueps << std::endl;

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

				// for lower bound on L1 slacks
				for(int i = 0; i < 2*X_DIM; ++i) { lb[t][index++] = 0; }

				index = 0;
				// x upper bound
				for(int i = 0; i < X_DIM; ++i) { ub[t][index++] = MIN(xMax[i], xt[i] + Xeps); }
				// u upper bound
				for(int i = 0; i < U_DIM; ++i) { ub[t][index++] = MIN(uMax[i], ut[i] + Ueps); }
			}

			Matrix<X_DIM>& xT = X[T-1];

			// Fill in lb, ub, C, e
			index = 0;
			double delta = .1;
			// xGoal lower bound
			for(int i = 0; i < P_DIM; ++i) { lb[T-1][index++] = xGoal[i] - delta; }
			// loose on landmarks
			for(int i = 0; i < X_DIM-P_DIM; ++i) { lb[T-1][index++] = MAX(xMin[i+P_DIM], xT[i+P_DIM] - Xeps); }

			index = 0;
			// xGoal upper bound
			for(int i = 0; i < P_DIM; ++i) { ub[T-1][index++] = xGoal[i] + delta; }
			// loose on landmarks
			for(int i = 0; i < X_DIM-P_DIM; ++i) { ub[T-1][index++] = MIN(xMax[i+P_DIM], xT[i+P_DIM] + Xeps); }

			// Verify problem inputs
			//if (!isValidInputs()) {
			//	std::cout << "Inputs are not valid!" << std::endl;
			//	exit(-1);
			//}

			//std::cerr << "PAUSING INSIDE MINIMIZE MERIT FUNCTION FOR INPUT VERIFICATION" << std::endl;
			//int num;
			//std::cin >> num;

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
				for(int i = 0; i < X_DIM; ++i) {
					Xopt[T-1][i] = z[T-1][i];
				}
			}
			else {
				LOG_ERROR("Some problem in solver");
				pythonDisplayTrajectory(U, T, true);
				exit(-1);
			}

			LOG_DEBUG("       Optimized cost: %4.10f", optcost);

			model_merit = optcost + constant_cost;

			new_merit = casadiComputeMerit(Xopt, Uopt, penalty_coeff);

			LOG_DEBUG("       merit: %4.10f", merit);
			LOG_DEBUG("       model_merit: %4.10f", model_merit);
			LOG_DEBUG("       new_merit: %4.10f", new_merit);

			approx_merit_improve = merit - model_merit;
			exact_merit_improve = merit - new_merit;
			merit_improve_ratio = exact_merit_improve / approx_merit_improve;

			LOG_DEBUG("       approx_merit_improve: %1.6f", approx_merit_improve);
			LOG_DEBUG("       exact_merit_improve: %1.6f", exact_merit_improve);
			LOG_DEBUG("       merit_improve_ratio: %1.6f", merit_improve_ratio);

			//std::cout << "PAUSED INSIDE minimizeMeritFunction AFTER OPTIMIZATION" << std::endl;
			//int num;
			//std::cin >> num;

			if (approx_merit_improve < -1e-5) {
				//LOG_ERROR("Approximate merit function got worse: %1.6f", approx_merit_improve);
				//LOG_ERROR("Either convexification is wrong to zeroth order, or you are in numerical trouble");
				//LOG_ERROR("Failure!");

				return false;
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

				casadiComputeCostGrad(Xopt, Uopt, cost, Gradopt);

				Matrix<XU_DIM> s, y;

				idx = 0;
				for(int t = 0; t < T-1; ++t) {
					for(int i=0; i < X_DIM; ++i) {
						s[idx+i] = Xopt[t][i] - X[t][i];
						y[idx+i] = Gradopt[idx+i] - Grad[idx+i];
					}
					idx += X_DIM;

					for(int i=0; i < U_DIM; ++i) {
						s[idx+i] = Uopt[t][i] - U[t][i];
						y[idx+i] = Gradopt[idx+i] - Grad[idx+i];
					}
					idx += U_DIM;
				}
				for(int i=0; i < X_DIM; ++i) {
					s[idx+i] = Xopt[T-1][i] - X[T-1][i];
					y[idx+i] = Gradopt[idx+i] - Grad[idx+i];
				}

				double theta;
				Matrix<XU_DIM> Bs = B*s;

				bool decision = ((~s*y)[0] >= .2*(~s*Bs)[0]);
				if (decision) {
					theta = 1;
				} else {
					theta = (.8*(~s*Bs)[0])/((~s*Bs-~s*y)[0]);
				}

				//std::cout << "theta: " << theta << std::endl;

				Matrix<XU_DIM> r = theta*y + (1-theta)*Bs;
				//Matrix<XU_DIM> rBs = theta*(y -Bs);

				// SR1 update
				//B = B + (rBs*~rBs)/((~rBs*s)[0]);

				// L-BFGS update
				B = B - (Bs*~Bs)/((~s*Bs)[0]) + (r*~r)/((~s*r)[0]);

				// Do not update B
				//B = identity<XU_DIM>();

				X = Xopt; U = Uopt;

				LOG_DEBUG("Accepted, Increasing trust region size to:  %2.6f %2.6f", Xeps, Ueps);
				break;
			}

			if (Xeps < cfg::min_trust_box_size && Ueps < cfg::min_trust_box_size) {
			    LOG_DEBUG("Converged: x tolerance");
			    return true;
			}

			std::cout << "U" << std::endl;
			for(int t=0; t < T-1; ++t) {
				std::cout << ~U[t];
			}
			std::cout << std::endl << std::endl;
			//pythonDisplayTrajectory(U, T, false);
			//pythonDisplayTrajectory(X, T, true);

		} // trust region loop
		sqp_iter++;
	} // sqp loop

	return success;
}


double statePenaltyCollocation(std::vector< Matrix<X_DIM> >& X, std::vector< Matrix<U_DIM> >& U, stateMPC_params& problem, stateMPC_output& output, stateMPC_info& info)
{
	double penalty_coeff = cfg::initial_penalty_coeff;
	double trust_box_size = cfg::initial_trust_box_size;

	int penalty_increases = 0;

	Matrix<X_DIM> dynviol;

	// penalty loop
	while(penalty_increases < cfg::max_penalty_coeff_increases)
	{
		bool success = minimizeMeritFunction(X, U, problem, output, info, penalty_coeff, trust_box_size);

		double cntviol = 0;
		for(int t = 0; t < T-1; ++t) {
			dynviol = (X[t+1] - dynfunc(X[t], U[t], zeros<Q_DIM,1>()));
			for(int i = 0; i < X_DIM; ++i) {
				cntviol += fabs(dynviol[i]);
			}
		}

	    success = success && (cntviol < cfg::cnt_tolerance);

		LOG_DEBUG("Constraint violations: %2.10f",cntviol);
		//std::cout << "Constraint violations: " << cntviol << std::endl;

	    if (!success) {
	        penalty_increases++;
	        penalty_coeff = penalty_coeff*cfg::penalty_coeff_increase_ratio;
	        trust_box_size = cfg::initial_trust_box_size;
	    }
	    else {
	    	//return computeCost(X, U);
	    	return casadiComputeCost(X, U);
	    }
	}
	//return computeCost(X, U);
	return casadiComputeCost(X, U);
}



int main(int argc, char* argv[])
{

	LOG_INFO("Initializing problem parameters");
	initProblemParams();

	stateMPC_params problem;
	stateMPC_output output;
	stateMPC_info info;
	setupStateVars(problem, output);
	util::Timer solveTimer;

	std::vector<Matrix<B_DIM> > B_total(T*NUM_WAYPOINTS);
	std::vector<Matrix<B_DIM> > B(T);
	std::vector<Matrix<X_DIM> > X(T);

	Matrix<U_DIM> uinit;

	Matrix<X_DIM,1> x;
	Matrix<X_DIM,X_DIM> s;
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
			X[t] = B[t].subMatrix<X_DIM,1>(0,0);
			//std::cout << ~X[t];
			B[t+1] = beliefDynamics(B[t], U[t]);
		}
		X[T-1] = B[T-1].subMatrix<X_DIM,1>(0,0);
		//std::cout << ~X[T-1] << std::endl << std::endl;
		//unVec(B[T-1], x, s);
		//std::cout << s.subMatrix<P_DIM,P_DIM>(0,0) << std::endl;


		std::cout << "U" << std::endl;
		for(int t=0; t < T-1; ++t) {
			std::cout << ~U[t];
		}
		std::cout << std::endl;


		double initTrajCost = computeCost(X, U);
		LOG_INFO("Initial trajectory cost: %4.10f", initTrajCost);

		double initCasadiTrajCost = casadiComputeCost(X, U);
		LOG_INFO("Initial casadi trajectory cost: %4.10f", initCasadiTrajCost);

		pythonDisplayTrajectory(B, U, waypoints, landmarks, T, true);

		Timer_tic(&solveTimer);

		double cost = statePenaltyCollocation(X, U, problem, output, info);

		double solvetime = util::Timer_toc(&solveTimer);

		vec(x0, SqrtSigma0, B[0]);
		X[0] = x0;
		for (int t = 0; t < T-1; ++t) {
			B[t+1] = beliefDynamics(B[t], U[t]);
			unVec(B[t+1], x, s);
			X[t+1] = x;
			//std::cout << s.subMatrix<P_DIM,P_DIM>(0,0) << std::endl;
		}
		std::cout << s.subMatrix<P_DIM,P_DIM>(0,0) << std::endl;

		LOG_INFO("Initial cost: %4.10f", initTrajCost);
		LOG_INFO("Optimized cost: %4.10f", cost);
		LOG_INFO("Actual cost: %4.10f", computeCost(X,U));
		LOG_INFO("Solve time: %5.3f ms", solvetime*1000);

		//pythonDisplayTrajectory(B, U, waypoints, landmarks, T);

		unVec(B[T-1], x0, SqrtSigma0);

		pythonDisplayTrajectory(B, U, waypoints, landmarks, T, true);

	}

	cleanupStateMPCVars();

	//vec(x0, SqrtSigma0, B[0]);
	//for (size_t t = 0; t < T-1; ++t) {
	//	B[t+1] = beliefDynamics(B[t], U[t]);
	//}

	return 0;
}
