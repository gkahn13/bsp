#define USE_COST_HAM

#include "../slam.h"
#include "../traj/slam-traj.h"

#include "../casadi/casadi-slam-control.h"


#include <vector>
#include <iomanip>

#include "util/matrix.h"
#include "util/Timer.h"

extern "C" {
#include "controlMPC.h"
controlMPC_FLOAT **H, **f, **lb, **ub, **z;
}

#include "boost/preprocessor.hpp"

const double alpha_belief = 10; // 10;
const double alpha_final_belief = 10; // 10;
const double alpha_control = .1; // .1
const double alpha_goal_state = 10; // 10

CasADi::SXFunction casadi_cost_func, casadi_gradcost_func;
CasADi::SXFunction casadi_cost_func_ham, casadi_gradcost_func_ham;

namespace cfg {
const double improve_ratio_threshold = .1; // .1
const double min_approx_improve = 1e-2; // 1e-2
const double min_trust_box_size = 1e-3; // 1e-3

const double trust_shrink_ratio = .75; // .75
const double trust_expand_ratio = 1.25; // 1.25

const double cnt_tolerance = .5; // .5
const double penalty_coeff_increase_ratio = 5; // 5
const double initial_penalty_coeff = 10; // 10

const double initial_trust_box_size = .1; // .1 // split up trust box size for X and U
const double initial_Uvel_trust_box_size = 5; // 5;
const double initial_Uangle_trust_box_size = M_PI/8; // M_PI/8;

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

inline double wrapAngle(double angle) {
	return angle - 2*M_PI * floor(angle/(2*M_PI));
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


void setupCasadiVars(const std::vector<Matrix<U_DIM> >& U, double* U_arr, double* x0_arr, double* Sigma0_arr, double* xGoal_arr, double* params_arr)
{
	int index = 0;
	for(int t = 0; t < T-1; ++t) {
		for(int i=0; i < U_DIM; ++i) {
			U_arr[index++] = U[t][i];
		}
	}

	for(int i=0; i < X_DIM; ++i) {
		x0_arr[i] = x0[i];
	}

	Matrix<X_DIM,X_DIM> Sigma0 = SqrtSigma0*SqrtSigma0;
	index = 0;
	for(int i=0; i < X_DIM; ++i) {
		for(int j=0; j < X_DIM; ++j) {
			Sigma0_arr[index++] = Sigma0(i,j);
		}
	}

	for(int i=0; i < X_DIM; ++i) {
		xGoal_arr[i] = xGoal[i];
	}

	params_arr[0] = alpha_belief;
	params_arr[1] = alpha_control;
	params_arr[2] = alpha_final_belief;
	params_arr[3] = alpha_goal_state;
}

double casadiComputeCost(const std::vector< Matrix<U_DIM> >& U)
{
	double U_arr[TU_DIM];
	double x0_arr[X_DIM];
	double Sigma0_arr[X_DIM*X_DIM];
	double xGoal_arr[X_DIM];
	double params_arr[4];

	setupCasadiVars(U, U_arr, x0_arr, Sigma0_arr, xGoal_arr, params_arr);

	double cost = 0;

	casadi_cost_func.setInput(U_arr,0);
	casadi_cost_func.setInput(x0_arr,1);
	casadi_cost_func.setInput(Sigma0_arr,2);
	casadi_cost_func.setInput(xGoal_arr,3);
	casadi_cost_func.setInput(params_arr,4);

	casadi_cost_func.evaluate();

	casadi_cost_func.getOutput(&cost,0);

	return cost;
}

void casadiComputeCostGrad(const std::vector< Matrix<U_DIM> >& U, double& cost, Matrix<TU_DIM>& Grad)
{
	double U_arr[TU_DIM];
	double x0_arr[X_DIM];
	double Sigma0_arr[X_DIM*X_DIM];
	double xGoal_arr[X_DIM];
	double params_arr[4];

	setupCasadiVars(U, U_arr, x0_arr, Sigma0_arr, xGoal_arr, params_arr);

	casadi_gradcost_func.setInput(U_arr,0);
	casadi_gradcost_func.setInput(x0_arr,1);
	casadi_gradcost_func.setInput(Sigma0_arr,2);
	casadi_gradcost_func.setInput(xGoal_arr,3);
	casadi_gradcost_func.setInput(params_arr,4);

	casadi_gradcost_func.evaluate();

	casadi_gradcost_func.getOutput(&cost,0);
	casadi_gradcost_func.getOutput(Grad.getPtr(),1);

}

double computeCost(const std::vector< Matrix<U_DIM> >& U)
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
	cost += alpha_final_belief*tr(SqrtSigma*SqrtSigma) + alpha_goal_state*tr(~(x - xGoal)*(x - xGoal));
	return cost;
}

void computeCostGrad(std::vector< Matrix<U_DIM> >& U, double& cost, Matrix<TU_DIM>& Grad) {
	cost = computeCost(U);

	int index = 0;
	for(int t=0; t < T-1; ++t) {
		for(int i=0; i < U_DIM; ++i) {
			double u_orig = U[t][i];

			U[t][i] = u_orig + step;
			double cost_p = computeCost(U);
			U[t][i] = u_orig - step;
			double cost_m = computeCost(U);

			U[t][i] = u_orig;
			Grad[index++] = (cost_p - cost_m) / (2*step);
		}
	}
}

double computeCostHam(const std::vector<Matrix<U_DIM> >& U) {
	double cost = 0;
	// factored form, Sigma_t = B_ham_t * !C_ham_t
	std::vector<Matrix<X_DIM,X_DIM> > B_ham(T), C_ham(T);
	B_ham[0] = SqrtSigma0*SqrtSigma0;
	C_ham[0] = identity<X_DIM>();

	std::vector<Matrix<X_DIM> > X(T);
	X[0] = x0;
	for(int t=0; t < T-1; ++t) {
		Matrix<C_DIM,C_DIM> Acar;
		Matrix<C_DIM,Q_DIM> Mcar;
		linearizeDynamics(X[t], U[t], zeros<Q_DIM,1>(), Acar, Mcar);

		Matrix<X_DIM,X_DIM> A = identity<X_DIM>();
		A.insert<C_DIM,C_DIM>(0, 0, Acar);
		Matrix<X_DIM,Q_DIM> M = zeros<X_DIM,Q_DIM>();
		M.insert<C_DIM, 2>(0, 0, Mcar);

		X[t+1] = dynfunc(X[t], U[t], zeros<Q_DIM,1>());

		Matrix<Z_DIM,X_DIM> H;
		Matrix<Z_DIM,R_DIM> N;
		linearizeObservation(X[t+1], zeros<R_DIM,1>(), H, N);

		Matrix<Z_DIM,Z_DIM> delta = deltaMatrix(X[t+1]);

		// factored form update
		B_ham[t+1] = A*B_ham[t] + M*Q*~M*(!(~A)*C_ham[t]);
		C_ham[t+1] = !(~A)*C_ham[t] + (~H*delta)*!(Matrix<R_DIM,R_DIM>)R*delta*H*B_ham[t+1];

		cost += alpha_control*tr(~U[t]*U[t]);
	}

	Matrix<X_DIM,X_DIM> SigmaFinalInfo = B_ham[T-1]*!C_ham[T-1];

	cost += alpha_final_belief*tr(SigmaFinalInfo) + alpha_goal_state*tr(~(X[T-1] - xGoal)*(X[T-1] - xGoal));

	return cost;
}

void computeCostGradHam(std::vector< Matrix<U_DIM> >& U, double& cost, Matrix<TU_DIM>& Grad) {
	cost = computeCostHam(U);

	int index = 0;
	for(int t=0; t < T-1; ++t) {
		for(int i=0; i < U_DIM; ++i) {
			double u_orig = U[t][i];

			U[t][i] = u_orig + step;
			double cost_p = computeCostHam(U);
			U[t][i] = u_orig - step;
			double cost_m = computeCostHam(U);

			U[t][i] = u_orig;
			Grad[index++] = (cost_p - cost_m) / (2*step);
		}
	}
}


double casadiComputeCostHam(const std::vector< Matrix<U_DIM> >& U)
{
	double U_arr[TU_DIM];
	double x0_arr[X_DIM];
	double Sigma0_arr[X_DIM*X_DIM];
	double xGoal_arr[X_DIM];
	double params_arr[4];

	setupCasadiVars(U, U_arr, x0_arr, Sigma0_arr, xGoal_arr, params_arr);

	double cost = 0;

	casadi_cost_func_ham.setInput(U_arr,0);
	casadi_cost_func_ham.setInput(x0_arr,1);
	casadi_cost_func_ham.setInput(Sigma0_arr,2);
	casadi_cost_func_ham.setInput(xGoal_arr,3);
	casadi_cost_func_ham.setInput(params_arr,4);

	casadi_cost_func_ham.evaluate();

	casadi_cost_func_ham.getOutput(&cost,0);

	return cost;
}

void casadiComputeCostGradHam(const std::vector< Matrix<U_DIM> >& U, double& cost, Matrix<TU_DIM>& Grad)
{
	double U_arr[TU_DIM];
	double x0_arr[X_DIM];
	double Sigma0_arr[X_DIM*X_DIM];
	double xGoal_arr[X_DIM];
	double params_arr[4];

	setupCasadiVars(U, U_arr, x0_arr, Sigma0_arr, xGoal_arr, params_arr);

	casadi_gradcost_func_ham.setInput(U_arr,0);
	casadi_gradcost_func_ham.setInput(x0_arr,1);
	casadi_gradcost_func_ham.setInput(Sigma0_arr,2);
	casadi_gradcost_func_ham.setInput(xGoal_arr,3);
	casadi_gradcost_func_ham.setInput(params_arr,4);

	casadi_gradcost_func_ham.evaluate();

	casadi_gradcost_func_ham.getOutput(&cost,0);
	casadi_gradcost_func_ham.getOutput(Grad.getPtr(),1);

}


double computeCostInfo(const std::vector<Matrix<U_DIM> >& U) {
	double cost = 0;

	std::vector<Matrix<B_DIM> > B(T);
	vec(x0, SqrtSigma0, B[0]);
	for(int t=0; t < T-1; ++t) {
		B[t+1] = beliefDynamics(B[t], U[t]);
	}

	Matrix<X_DIM> xFinal;
	Matrix<X_DIM, X_DIM> SqrtSigmaFinal;
	unVec(B[T-1], xFinal, SqrtSigmaFinal);

	std::vector<Matrix<X_DIM> > X(T);
	std::vector<Matrix<X_DIM, X_DIM> > omega(T), omega_bar(T);
	X[0] = x0;
	omega[0] = !(SqrtSigma0*SqrtSigma0);
	for(int t=0; t < T-1; ++t) {
		Matrix<C_DIM,C_DIM> Acar;
		Matrix<C_DIM,Q_DIM> Mcar;
		linearizeDynamics(X[t], U[t], zeros<Q_DIM,1>(), Acar, Mcar);

		Matrix<X_DIM,X_DIM> A = identity<X_DIM>();
		A.insert<C_DIM,C_DIM>(0, 0, Acar);
		Matrix<X_DIM,Q_DIM> M = zeros<X_DIM,Q_DIM>();
		M.insert<C_DIM, 2>(0, 0, Mcar);

		X[t+1] = dynfunc(X[t], U[t], zeros<Q_DIM,1>());

		omega_bar[t+1] = !(A*!omega[t]*~A + M*Q*~M);
		std::cout << A*!omega[t]*~A << "\n";
//		Matrix<X_DIM,X_DIM> MQM = M*Q*~M;
//		std::cout << MQM << "\n";
//		Matrix<C_DIM,C_DIM> MQM_submat = MQM.subMatrix<C_DIM,C_DIM>(0,0);
//		Matrix<X_DIM,X_DIM> MQM_pinv = zeros<X_DIM,X_DIM>();
//		std::cout << !(MQM_submat + 1e-7*identity<C_DIM>()) << "\n";
//		MQM_pinv.insert(0,0, !(MQM_submat + 1e-7*identity<C_DIM>()));
//		omega_bar[t+1] = (!(~A))*omega[t]*!A + pseudoInverse(MQM);

		Matrix<Z_DIM,X_DIM> H;
		Matrix<Z_DIM,R_DIM> N;
		linearizeObservation(X[t+1], zeros<R_DIM,1>(), H, N);

		Matrix<Z_DIM,Z_DIM> delta = deltaMatrix(X[t+1]);
		omega[t+1] = omega_bar[t+1] + (~H*delta)*!(Matrix<R_DIM,R_DIM>)R*(delta*H);

		cost += alpha_belief*tr(!omega[t]);
		cost += alpha_control*tr(~U[t]*U[t]);
	}

	cost += alpha_final_belief*tr(!omega[T-1]) + alpha_goal_state*tr(~(X[T-1] - xGoal)*(X[T-1] - xGoal));

	return cost;
}


// Jacobians: dg(b,u)/db, dg(b,u)/du
void linearizeCarDynamics(const Matrix<C_DIM>& c, const Matrix<U_DIM>& u, Matrix<C_DIM,C_DIM>& F, Matrix<C_DIM,U_DIM>& G, Matrix<C_DIM>& h)
{
	Matrix<X_DIM,1> x;
	x.insert(0, 0, c);
	x.insert(C_DIM, 0, x0.subMatrix<L_DIM,1>(C_DIM,0));

	F.reset();
	Matrix<X_DIM> xr(x), xl(x), ddx;
	for (size_t i = 0; i < C_DIM; ++i) {
		xr[i] += step; xl[i] -= step;
		ddx = (dynfunc(xr, u, zeros<Q_DIM,1>()) - dynfunc(xl, u, zeros<Q_DIM,1>())) / (xr[i] - xl[i]);
		F.insert(0,i, ddx.subMatrix<C_DIM,1>(0, 0));
		xr[i] = x[i]; xl[i] = x[i];
	}

	G.reset();
	Matrix<U_DIM> ur(u), ul(u);
	Matrix<X_DIM> ddg;
	for (size_t i = 0; i < U_DIM; ++i) {
		ur[i] += step; ul[i] -= step;
		ddg = (dynfunc(x, ur, zeros<Q_DIM,1>()) - dynfunc(x, ul, zeros<Q_DIM,1>())) / (ur[i] - ul[i]);
		G.insert(0,i, ddg.subMatrix<C_DIM,1>(0, 0));
		ur[i] = u[i]; ul[i] = u[i];
	}

	h = dynfunc(x, u, zeros<Q_DIM,1>()).subMatrix<C_DIM,1>(0,0);
}


void setupControlVars(controlMPC_params& problem, controlMPC_output& output)
{
	// problem inputs
	H = new controlMPC_FLOAT*[T];
	f = new controlMPC_FLOAT*[T];
	lb = new controlMPC_FLOAT*[T];
	ub = new controlMPC_FLOAT*[T];

	// problem outputs
	z = new controlMPC_FLOAT*[T];

#define SET_VARS(n)    \
		H[ BOOST_PP_SUB(n,1) ] = problem.H##n ;  \
		f[ BOOST_PP_SUB(n,1) ] = problem.f##n ;  \
		lb[ BOOST_PP_SUB(n,1) ] = problem.lb##n ;	\
		ub[ BOOST_PP_SUB(n,1) ] = problem.ub##n ;	\
		z[ BOOST_PP_SUB(n,1) ] = output.z##n ;

#define BOOST_PP_LOCAL_MACRO(n) SET_VARS(n)
#define BOOST_PP_LOCAL_LIMITS (1, TIMESTEPS-1)
#include BOOST_PP_LOCAL_ITERATE()

	for(int t = 0; t < T-1; ++t) {
		for(int i=0; i < U_DIM; ++i) { H[t][i] = INFTY; }
		for(int i=0; i < U_DIM; ++i) { f[t][i] = INFTY; }
		for(int i=0; i < U_DIM; ++i) { lb[t][i] = INFTY; }
		for(int i=0; i < U_DIM; ++i) { ub[t][i] = INFTY; }
		for(int i=0; i < U_DIM; ++i) { z[t][i] = INFTY; }
	}
}

void cleanupControlMPCVars()
{
	delete[] H;
	delete[] f;
	delete[] lb;
	delete[] ub;
	delete[] z;
}

bool isValidInputs()
{

	for(int t = 0; t < T-1; ++t) {
		std::cout << "t: " << t << "\n";
		for(int i=0; i < U_DIM; ++i) { if (H[t][i] > INFTY/2) { std::cout << "H error: " << i << "\n"; } }
		for(int i=0; i < U_DIM; ++i) { if (f[t][i] > INFTY/2) { std::cout << "f error: " << i << "\n"; } }
		for(int i=0; i < U_DIM; ++i) { if (lb[t][i] > INFTY/2) { std::cout << "lb error: " << i << "\n"; } }
		for(int i=0; i < U_DIM; ++i) {if (lb[t][i] > INFTY/2) { std::cout << "ub error: " << i << "\n"; } }
	}

	for(int t = 0; t < T; ++t) {

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


		if (t < T-1) {
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
		}
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


	std::cout << std::endl;
	std::cout << std::endl;
	return true;
}

bool controlCollocation(std::vector< Matrix<U_DIM> >& U, controlMPC_params& problem, controlMPC_output& output, controlMPC_info& info)
{
	double Ueps = cfg::initial_trust_box_size;

	double Uvel_eps = cfg::initial_Uvel_trust_box_size;
	double Uangle_eps = cfg::initial_Uangle_trust_box_size;

	double optcost;

	std::vector<Matrix<U_DIM> > Uopt(T-1);

	double merit, model_merit, new_merit;
	double approx_merit_improve, exact_merit_improve, merit_improve_ratio;
	double constant_cost, hessian_constant, jac_constant;

	int index = 0;

	Matrix<C_DIM+U_DIM, C_DIM+U_DIM> HMat; // TODO: reduce
	Matrix<C_DIM,C_DIM> HfMat;



	// full Hessian from current timstep
	Matrix<TU_DIM,TU_DIM> B = identity<TU_DIM>();

	Matrix<TU_DIM> Grad, Gradopt;
	double cost;
	int idx = 0;

	// iter loop
	for(int it = 0; it < cfg::max_sqp_iterations; ++it)
	{
		LOG_DEBUG("  Iter: %d", it);

		// Compute gradients
#ifndef USE_COST_HAM
		casadiComputeCostGrad(U, cost, Grad);
//		computeCostGrad(U, cost, Grad);
#else
		casadiComputeCostGradHam(U, cost, Grad);
//		computeCostGradHam(U, cost, Grad);
#endif

		// Problem linearization and definition
		// fill in H, f

		hessian_constant = 0;
		jac_constant = 0;
		idx = 0;

		for (int t = 0; t < T-1; ++t)
		{
			Matrix<U_DIM>& ut = U[t];

			idx = t*(U_DIM);
			//LOG_DEBUG("idx: %d",idx);

			// fill in gradients and Hessians

			HMat.reset();
			for(int i = 0; i < U_DIM; ++i) {
				HMat(i,i) = B(idx+i,idx+i);
			}

			// since diagonal, fill directly
			for(int i = 0; i < U_DIM; ++i) { H[t][i] = HMat(i,i); }

			for(int i = 0; i < U_DIM; ++i) {
				hessian_constant += HMat(i,i)*ut[i]*ut[i];
				jac_constant -= Grad[idx+i]*ut[i];
				f[t][i] = Grad[idx+i] - HMat(i,i)*ut[i];
			}
		}

		//Matrix<C_DIM,1> cGoal = xGoal.subMatrix<C_DIM,1>(0,0);
		//double goal_cost = -alpha_goal_state*tr(~cGoal*cGoal);
		double goal_cost = 0;

		constant_cost = 0.5*hessian_constant + jac_constant + cost + goal_cost;
		LOG_DEBUG("  hessian cost: %4.10f", 0.5*hessian_constant);
		LOG_DEBUG("  jacobian cost: %4.10f", jac_constant);
		LOG_DEBUG("  goal cost: %4.10f", goal_cost);
		LOG_DEBUG("  constant cost: %4.10f", constant_cost);



		// trust region size adjustment
		while(true)
		{
			LOG_DEBUG("       trust region size: %2.6f", Ueps);

			// solve the innermost QP here
			for(int t = 0; t < T-1; ++t)
			{
				Matrix<U_DIM>& ut = U[t];

				index = 0;
				// u velocity lower bound
				lb[t][index++] = MAX(uMin[0], ut[0] - Uvel_eps);
				// u angle lower bound
				lb[t][index++] = MAX(uMin[1], ut[1] - Uangle_eps);

				index = 0;
				// u velocity upper bound
				ub[t][index++] = MIN(uMax[0], ut[0] + Uvel_eps);
				// u angle upper bound
				ub[t][index++] = MIN(uMax[1], ut[1] + Uangle_eps);
			}

			// Verify problem inputs
			//if (!isValidInputs()) {
			//	std::cout << "Inputs are not valid!" << std::endl;
			//	exit(-1);
			//}



			int exitflag = controlMPC_solve(&problem, &output, &info);
			if (exitflag == 1) {
				for(int t = 0; t < T-1; ++t) {
					Matrix<U_DIM>& ut = Uopt[t];

					for(int i = 0; i < U_DIM; ++i) {
						ut[i] = z[t][i];
					}
				}
				optcost = info.pobj;
			}
			else {
				LOG_ERROR("Some problem in solver");
				throw forces_exception();
			}

//			LOG_DEBUG("Optimized trajectory");
//			pythonDisplayTrajectory(Uopt, T, true);


			LOG_DEBUG("       Optimized cost: %4.10f", optcost);

			model_merit = optcost + constant_cost;

#ifndef USE_COST_HAM
			new_merit = casadiComputeCost(Uopt);
//			new_merit = computeCost(Uopt);
#else
			new_merit = casadiComputeCostHam(Uopt);
//			new_merit = computeCostHam(Uopt);
#endif

			merit = cost;
			LOG_DEBUG("       merit: %4.10f", merit);
			LOG_DEBUG("       model_merit: %4.10f", model_merit);
			LOG_DEBUG("       new_merit: %4.10f", new_merit);

			approx_merit_improve = merit - model_merit;
			exact_merit_improve = merit - new_merit;
			merit_improve_ratio = exact_merit_improve / approx_merit_improve;

			LOG_DEBUG("       approx_merit_improve: %1.6f", approx_merit_improve);
			LOG_DEBUG("       exact_merit_improve: %1.6f", exact_merit_improve);
			LOG_DEBUG("       merit_improve_ratio: %1.6f", merit_improve_ratio);


			if (approx_merit_improve < -1e-5) {
				LOG_ERROR("Approximate merit function got worse: %1.6f", approx_merit_improve);
				LOG_ERROR("Either convexification is wrong to zeroth order, or you are in numerical trouble");

				Uvel_eps *= cfg::trust_shrink_ratio;
				Uangle_eps *= cfg::trust_shrink_ratio;
				//break;
				//return false;
			} else if (approx_merit_improve < cfg::min_approx_improve) {
				LOG_DEBUG("Converged: improvement small enough");
				U = Uopt;
				return true;
			} else if ((exact_merit_improve < 0) || (merit_improve_ratio < cfg::improve_ratio_threshold)) {
				Uvel_eps *= cfg::trust_shrink_ratio;
				Uangle_eps *= cfg::trust_shrink_ratio;
				LOG_DEBUG("Shrinking trust region size to: %2.6f %2.6f", Uvel_eps, Uangle_eps);
			} else {
				Uvel_eps *= cfg::trust_expand_ratio;
				Uangle_eps *= cfg::trust_expand_ratio;

#ifndef USE_COST_HAM
				casadiComputeCostGrad(Uopt, cost, Gradopt);
//				computeCostGrad(Uopt, cost, Gradopt);
#else
				casadiComputeCostGradHam(Uopt, cost, Gradopt);
//				computeCostGradHam(Uopt, cost, Gradopt);
#endif

				Matrix<TU_DIM> s, y;

				idx = 0;
				for(int t = 0; t < T-1; ++t) {
					for(int i=0; i < U_DIM; ++i) {
						s[idx+i] = Uopt[t][i] - U[t][i];
						y[idx+i] = Gradopt[idx+i] - Grad[idx+i];
					}
					idx += U_DIM;
				}

				double theta;
				Matrix<TU_DIM> Bs = B*s;

				bool decision = ((~s*y)[0] >= .2*(~s*Bs)[0]);
				if (decision) {
					theta = 1;
				} else {
					theta = (.8*(~s*Bs)[0])/((~s*Bs-~s*y)[0]);
				}

				Matrix<TU_DIM> r = theta*y + (1-theta)*Bs;
				//Matrix<XU_DIM> rBs = theta*(y -Bs);

				// SR1 update
				//B = B + (rBs*~rBs)/((~rBs*s)[0]);

				// L-BFGS update
				B = B - (Bs*~Bs)/((~s*Bs)[0]) + (r*~r)/((~s*r)[0]);

				// Do not update B
				//B = identity<CU_DIM>();

				U = Uopt;

				LOG_DEBUG("Accepted, Increasing trust region size to:  %2.6f %2.6f", Uvel_eps, Uangle_eps);
				break;
			}

			//if (Xeps < cfg::min_trust_box_size && Ueps < cfg::min_trust_box_size) {
			if (Uvel_eps < cfg::min_trust_box_size && Uangle_eps < cfg::min_trust_box_size) {
			    LOG_DEBUG("Converged: x tolerance");
			    return true;
			}
		} // trust region loop
	} // iter loop

	return false;
}


void planPath(std::vector<Matrix<P_DIM> > l, controlMPC_params& problem, controlMPC_output& output, controlMPC_info& info, std::ofstream& f) {
	initProblemParams(l);

	util::Timer solveTimer, trajTimer;
	double totalSolveTime = 0, trajTime = 0;

	double totalTrajCost = 0;

	std::vector<Matrix<B_DIM> > B_total(T*NUM_WAYPOINTS);
	std::vector<Matrix<U_DIM> > U_total((T-1)*NUM_WAYPOINTS);
	int B_total_idx = 0, U_total_idx = 0;

	std::vector<Matrix<B_DIM> > B(T);

	Matrix<U_DIM> uinit;

	Matrix<X_DIM,1> x;
	Matrix<X_DIM,X_DIM> s;

	x0[2] = nearestAngleFromTo(0, x0[2]); // need to remod back to near 0
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
			LOG_ERROR("Failed to initialize trajectory for slam-control, continuing anyways");
			//exit(-1);
		}

		double initTrajTime = util::Timer_toc(&trajTimer);
		trajTime += initTrajTime;

		vec(x0, SqrtSigma0, B[0]);
		for(int t=0; t < T-1; ++t) {
			B[t+1] = beliefDynamics(B[t], U[t]);
		}

		//double initTrajCost = computeCost(U);
		//LOG_INFO("Initial trajectory cost: %4.10f", initTrajCost);

		//double initCasadiTrajCost = casadiComputeCost(U);
		//LOG_INFO("Initial casadi trajectory cost: %4.10f", initCasadiTrajCost);

//		pythonDisplayTrajectory(B, U, waypoints, landmarks, T, true);

		util::Timer_tic(&solveTimer);

		int iter = 0;
		while(true) {
			try {
				controlCollocation(U, problem, output, info);
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
		}

		for (int t = 0; t < T-1; ++t) {
			B_total[B_total_idx++] = B[t];
			U_total[U_total_idx++] = U[t];
		}
		B_total[B_total_idx++] = B[T-1];

		totalTrajCost += computeCost(U);


		unVec(B[T-1], x0, SqrtSigma0);

//		pythonDisplayTrajectory(B, U, waypoints, landmarks, T, true);
	}



	LOG_INFO("Total trajectory cost: %4.10f", totalTrajCost);
	LOG_INFO("Total trajectory solve time: %5.3f ms", trajTime*1000);
	LOG_INFO("Total solve time: %5.3f ms", totalSolveTime*1000);

	logDataToFile(f, B_total, totalSolveTime*1000, trajTime*1000, 0);

	pythonDisplayTrajectory(B_total, U_total, waypoints, landmarks, T*NUM_WAYPOINTS, true);
}

void test_hamiltonian(std::vector<Matrix<P_DIM> > l) {
	initProblemParams(l);

	xGoal.insert(0, 0, waypoints[0]);
	xGoal[2] = x0[2];

	Matrix<U_DIM> uinit;
	uinit[0] = (xGoal[0] - x0[0])/((float)(T-1));
	uinit[1] = 0;
	std::vector<Matrix<U_DIM> > U(T-1, uinit);
//	U[0][1] = M_PI/8;
//	U[7][1] = -M_PI/4;

	double cost_ham = computeCostHam(U);
	std::cout << "cost_ham: " << cost_ham << "\n";
	double casadi_cost_ham = casadiComputeCostHam(U);
	std::cout << "casadi_cost_ham: " << casadi_cost_ham << "\n";
//	pythonDisplayTrajectory(U, T, true);

	std::vector<Matrix<B_DIM> > B(T);
	vec(x0, SqrtSigma0, B[0]);
	for(int t=0; t < T-1; ++t) {
		B[t+1] = beliefDynamics(B[t], U[t]);
	}

	Matrix<X_DIM> xFinal;
	Matrix<X_DIM, X_DIM> SqrtSigmaFinal;
	unVec(B[T-1], xFinal, SqrtSigmaFinal);

//	pythonDisplayTrajectory(B, U, waypoints, l, T, true);

	// factored form, Sigma_t = B_ham_t * !C_ham_t
	std::vector<Matrix<X_DIM,X_DIM> > B_ham(T), C_ham(T);
	B_ham[0] = SqrtSigma0*SqrtSigma0;
	C_ham[0] = identity<X_DIM>();

	std::vector<Matrix<B_DIM> > B_info(T);
	std::vector<Matrix<X_DIM> > X(T);
	std::vector<Matrix<X_DIM, X_DIM> > omega(T), omega_bar(T);
	X[0] = x0;
	omega[0] = !(SqrtSigma0*SqrtSigma0);
	vec(X[0], sqrtm(!omega[0]), B_info[0]);
	for(int t=0; t < T-1; ++t) {
		Matrix<C_DIM,C_DIM> Acar;
		Matrix<C_DIM,Q_DIM> Mcar;
		linearizeDynamics(X[t], U[t], zeros<Q_DIM,1>(), Acar, Mcar);

		Matrix<X_DIM,X_DIM> A = identity<X_DIM>();
		A.insert<C_DIM,C_DIM>(0, 0, Acar);
		Matrix<X_DIM,Q_DIM> M = zeros<X_DIM,Q_DIM>();
		M.insert<C_DIM, 2>(0, 0, Mcar);

		X[t+1] = dynfunc(X[t], U[t], zeros<Q_DIM,1>());

		omega_bar[t+1] = !(A*!omega[t]*~A + M*Q*~M);

		Matrix<Z_DIM,X_DIM> H;
		Matrix<Z_DIM,R_DIM> N;
		linearizeObservation(X[t+1], zeros<R_DIM,1>(), H, N);

		Matrix<Z_DIM,Z_DIM> delta = deltaMatrix(X[t+1]);
		omega[t+1] = omega_bar[t+1] + (~H*delta)*!(Matrix<R_DIM,R_DIM>)R*(delta*H);

//		Matrix<X_DIM,Z_DIM> K = ((!omega_bar[t+1]*~H*delta)/(delta*H*!omega_bar[t+1]*~H*delta + R))*delta;
//		omega[t+1] = !((identity<X_DIM>() - K*H)*!omega_bar[t+1]);
//
//		vec(X[t+1], sqrtm(!omega[t+1]), B_info[t+1]);

		// factored form update
		B_ham[t+1] = A*B_ham[t] + M*Q*~M*(!(~A)*C_ham[t]);
		C_ham[t+1] = !(~A)*C_ham[t] + (~H*delta)*!(Matrix<R_DIM,R_DIM>)R*delta*H*B_ham[t+1];

		vec(X[t+1], sqrtm(B_ham[t+1]*!C_ham[t+1]), B_info[t+1]);

		Matrix<X_DIM,X_DIM> sigma_info = !omega[t+1];
		Matrix<X_DIM,X_DIM> sigma_factored = B_ham[t+1]*!C_ham[t+1];
		std::cout << "factored difference:\n" << sigma_info - sigma_factored << "\n";
	}

	// sum of square of singular values
	Matrix<X_DIM,X_DIM> SigmaFinal = SqrtSigmaFinal*SqrtSigmaFinal;
	double sigma_final_sum_sqr_singulars = tr(~SigmaFinal*SigmaFinal);
	Matrix<X_DIM,X_DIM> SigmaFinalInfo = !omega[T-1];
	double sigma_final_info_sum_sqr_singulars = tr(~SigmaFinalInfo*SigmaFinalInfo);

	std::cout << "SigmaFinal:\n" << SigmaFinal << "\n";
	std::cout << "SigmaFinal information:\n" << SigmaFinalInfo << "\n";

	pythonDisplayTrajectory(B_info, U, waypoints, l, T, true);

	std::cout << "SigmaFinal sum square singulars: " << sigma_final_sum_sqr_singulars << "\n";
	std::cout << "SigmaFinalInfo sum square singulars: " << sigma_final_info_sum_sqr_singulars << "\n";

	Matrix<X_DIM,X_DIM> diff = (SigmaFinal - SigmaFinalInfo);
	std::cout << "(SigmaFinal - SigmaFinalInfo):\n" << diff << "\n";

}

void test_info(std::vector<Matrix<P_DIM> > l) {
	initProblemParams(l);

	xGoal.insert(0, 0, waypoints[0]);
	xGoal[2] = x0[2];

	Matrix<U_DIM> uinit;
	uinit[0] = (xGoal[0] - x0[0])/((float)(T-1));
	uinit[1] = 0;
	std::vector<Matrix<U_DIM> > U(T-1, uinit);

	double cost_info = computeCostInfo(U);
	std::cout << "cost_info: " << cost_info << "\n";
	double cost = computeCost(U);
	std::cout << "cost: " << cost << "\n";
//	pythonDisplayTrajectory(U, T, true);
//
//	std::vector<Matrix<B_DIM> > B(T);
//	vec(x0, SqrtSigma0, B[0]);
//	for(int t=0; t < T-1; ++t) {
//		B[t+1] = beliefDynamics(B[t], U[t]);
//	}
}

int main(int argc, char* argv[])
{
	controlMPC_params problem;
	controlMPC_output output;
	controlMPC_info info;
	setupControlVars(problem, output);

	std::vector<std::vector<Matrix<P_DIM> > > l_list = landmarks_list();
	test_info(l_list[0]);
	return 0;

	LOG_INFO("initializing casadi functions...");

	std::ofstream f;
#ifndef USE_COST_HAM
	logDataHandle("slam/data/slam-control", f);
	casadi_cost_func = casadiCostFunc();
	casadi_gradcost_func = casadiCostGradFunc();
#else
	logDataHandle("slam/data/slam-control-ham", f);
	casadi_cost_func_ham = casadiCostFuncHam();
	casadi_gradcost_func_ham = casadiCostGradFuncHam();
#endif

	LOG_INFO("casadi functions initialized");

//	casadi_gradcost_func.generateCode("casadi_gradcost_func.c");
//	casadi_gradcost_func_ham.generateCode("casadi_gradcost_func_ham.c");

//	// TODO: temp
//	test_hamiltonian(l_list[0]);
//	return 0;


	for(size_t i=0; i < l_list.size(); ++i) {
		planPath(l_list[i], problem, output, info, f);
	}
	cleanupControlMPCVars();

	return 0;
}
