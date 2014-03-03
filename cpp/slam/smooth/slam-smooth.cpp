#include "slam-smooth.h"

#include <vector>
#include <iomanip>
#include <sys/resource.h>

#include "util/matrix.h"
#include "util/Timer.h"
#include "util/logging.h"

extern "C" {
#include "smoothMPC.h"
smoothMPC_FLOAT **H_smooth, **f_smooth, **lb_smooth, **ub_smooth, **C_smooth, **e_smooth, **z_smooth;
}

Matrix<U_DIM> uMinSmooth, uMaxSmooth;

const double alpha_control = .01; // .1
const double alpha_midway = 1; // 1
const double alpha_goal = 10; // 10

namespace cfg {
const double improve_ratio_threshold = .1; // .1
const double min_approx_improve = 1e-2; // 1e-2
const double min_trust_box_size = 1e-2; // 1e-2

const double trust_shrink_ratio = .5; // .5
const double trust_expand_ratio = 1.25; // 1.25

const double cnt_tolerance = 1e-1;
const double penalty_coeff_increase_ratio = 5; // 5
const double initial_penalty_coeff = 10; // 10

const double initial_trust_box_size = 1; // 1
const double initial_Xpos_trust_box_size = 1; // 1;
const double initial_Xangle_trust_box_size = M_PI/6; // M_PI/6;
const double initial_Uvel_trust_box_size = 1; // 1;
const double initial_Uangle_trust_box_size = M_PI/8; // M_PI/8;

const int max_penalty_coeff_increases = 4; // 4
const int max_sqp_iterations = 50; // 50
}

struct forces_exception {
	int c;
	forces_exception(int c):c(c) { }
};


// Jacobians: df(x,u)/dx, df(x,u)/du
void linearizeCarDynamicsSmooth(const Matrix<C_DIM>& x, const Matrix<U_DIM>& u, Matrix<C_DIM,C_DIM>& F, Matrix<C_DIM,U_DIM>& G, Matrix<C_DIM>& h)
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

inline double wrapAngle(double angle) {
	return angle - 2*M_PI * floor(angle/(2*M_PI));
}



double computeSmoothCost(const std::vector< Matrix<C_DIM> >& X, const std::vector< Matrix<U_DIM> >& U, const std::vector< Matrix<C_DIM> >& X_unsmooth)
{
	int T_unsmooth = X_unsmooth.size();
	int timestep_ratio = (int)T_SMOOTH / T_unsmooth;

	double cost = 0;
	for(int t = 0; t < T_SMOOTH-1; ++t) {
		// assuming controls don't wrap around
		cost += alpha_control*tr(~U[t]*U[t]);
	}
	for(int t=0; t < T_unsmooth-1; ++t) {
		Matrix<C_DIM> x = X[t*timestep_ratio];
		Matrix<C_DIM> x_unsmooth = X_unsmooth[t];
		Matrix<C_DIM> x_diff = x - x_unsmooth;
		x_diff[P_DIM] = 0; // don't care about angle
		cost += alpha_midway*tr(~x_diff*x_diff);
	}
	Matrix<C_DIM> x = X[T_SMOOTH-1];
	Matrix<C_DIM> x_unsmooth = X_unsmooth[T_unsmooth-1];
	cost += alpha_goal*tr(~(x-x_unsmooth)*(x-x_unsmooth));

	return cost;
}

double deviationCost(const std::vector< Matrix<C_DIM> >& X, const std::vector< Matrix<U_DIM> >& U, const std::vector< Matrix<C_DIM> >& X_unsmooth) {
	std::vector< Matrix<C_DIM> > X_integrated(X.size());
	X_integrated[0] = X[0];
	for(int t=0; t < T_SMOOTH-1; ++t) {
		X_integrated[t+1] = dynfunccar(X_integrated[t], U[t]);
	}

	double deviation = 0;

	int T_unsmooth = X_unsmooth.size();
	int timestep_ratio = (int)T_SMOOTH / T_unsmooth;

	for(int t=0; t < T_unsmooth-1; ++t) {
		Matrix<C_DIM> x = X_integrated[t*timestep_ratio];
		Matrix<C_DIM> x_unsmooth = X_unsmooth[t];
		Matrix<C_DIM> x_diff = x - x_unsmooth;
		x_diff[P_DIM] = 0; // don't care about angle
		deviation += sqrt(tr(~x_diff*x_diff));
	}
	Matrix<C_DIM> x = X_integrated[T_SMOOTH-1];
	Matrix<C_DIM> x_unsmooth = X_unsmooth[T_unsmooth-1];
	deviation += sqrt(tr(~(x-x_unsmooth)*(x-x_unsmooth)));

	return deviation;
}


void setupSmoothVars(smoothMPC_params &problem, smoothMPC_output &output)
{
	// problem inputs
	H_smooth = new smoothMPC_FLOAT*[T_SMOOTH];
	f_smooth = new smoothMPC_FLOAT*[T_SMOOTH];
	lb_smooth = new smoothMPC_FLOAT*[T_SMOOTH];
	ub_smooth = new smoothMPC_FLOAT*[T_SMOOTH];
	C_smooth = new smoothMPC_FLOAT*[T_SMOOTH-1];
	e_smooth = new smoothMPC_FLOAT*[T_SMOOTH];

	// problem outputs
	z_smooth = new smoothMPC_FLOAT*[T_SMOOTH];

#define SET_VARS(n)    \
		H_smooth[ BOOST_PP_SUB(n,1) ] = problem.H##n ;  \
		f_smooth[ BOOST_PP_SUB(n,1) ] = problem.f##n ;  \
		lb_smooth[ BOOST_PP_SUB(n,1) ] = problem.lb##n ;	\
		ub_smooth[ BOOST_PP_SUB(n,1) ] = problem.ub##n ;	\
		C_smooth[ BOOST_PP_SUB(n,1) ] = problem.C##n ;  \
		e_smooth[ BOOST_PP_SUB(n,1) ] = problem.e##n ;  \
		z_smooth[ BOOST_PP_SUB(n,1) ] = output.z##n ;
#define BOOST_PP_LOCAL_MACRO(n) SET_VARS(n)
#define BOOST_PP_LOCAL_LIMITS (1, TIMESTEPS_SMOOTH-1)
#include BOOST_PP_LOCAL_ITERATE()

#define SET_LAST_VARS(n)    \
		H_smooth[ BOOST_PP_SUB(n,1) ] = problem.H##n ;  \
		f_smooth[ BOOST_PP_SUB(n,1) ] = problem.f##n ; \
		lb_smooth[ BOOST_PP_SUB(n,1) ] = problem.lb##n ;	\
		ub_smooth[ BOOST_PP_SUB(n,1) ] = problem.ub##n ;	\
		e_smooth[ BOOST_PP_SUB(n,1) ] = problem.e##n ;  \
		z_smooth[ BOOST_PP_SUB(n,1) ] = output.z##n ;
#define BOOST_PP_LOCAL_MACRO(n) SET_LAST_VARS(n)
#define BOOST_PP_LOCAL_LIMITS (TIMESTEPS_SMOOTH, TIMESTEPS_SMOOTH)
#include BOOST_PP_LOCAL_ITERATE()

	for(int t = 0; t < T_SMOOTH-1; ++t) {
		for(int i=0; i < (3*C_DIM+U_DIM); ++i) { H_smooth[t][i] = INFTY; }
		for(int i=0; i < (3*C_DIM+U_DIM); ++i) { f_smooth[t][i] = INFTY; }
		for(int i=0; i < (3*C_DIM+U_DIM); ++i) { lb_smooth[t][i] = INFTY; }
		for(int i=0; i < (C_DIM+U_DIM); ++i) { ub_smooth[t][i] = INFTY; }
		for(int i=0; i < (C_DIM*(3*C_DIM+U_DIM)); ++i) { C_smooth[t][i] = INFTY; }
		for(int i=0; i < C_DIM; ++i) { e_smooth[t][i] = INFTY; }
		for(int i=0; i < (C_DIM+U_DIM); ++i) { z_smooth[t][i] = INFTY; }
	}
	for(int i=0; i < C_DIM; ++i) { H_smooth[T_SMOOTH-1][i] = INFTY; }
	for(int i=0; i < (C_DIM); ++i) { f_smooth[T_SMOOTH-1][i] = INFTY; }
	for(int i=0; i < (C_DIM); ++i) { lb_smooth[T_SMOOTH-1][i] = INFTY; }
	for(int i=0; i < (C_DIM); ++i) { ub_smooth[T_SMOOTH-1][i] = INFTY; }
	for(int i=0; i < C_DIM; ++i) { e_smooth[T_SMOOTH-1][i] = INFTY; }
	for(int i=0; i < (C_DIM); ++i) { z_smooth[T_SMOOTH-1][i] = INFTY; }

}

void cleanupSmoothMPCVars()
{
	delete[] H_smooth;
	delete[] f_smooth;
	delete[] lb_smooth;
	delete[] ub_smooth;
	delete[] C_smooth;
	delete[] e_smooth;
	delete[] z_smooth;
}

double computeSmoothMerit(const std::vector< Matrix<C_DIM> >& X, const std::vector< Matrix<U_DIM> >& U, const std::vector< Matrix<C_DIM> >& X_unsmooth, double penalty_coeff)
{
	double cost = computeSmoothCost(X, U, X_unsmooth);

	double merit = cost;

	Matrix<C_DIM> dynviol;
	for(int t = 0; t < T_SMOOTH-1; ++t) {
		dynviol = (X[t+1] - dynfunccar(X[t], U[t]) );
		merit += penalty_coeff*fabs(dynviol[0]);
		merit += penalty_coeff*fabs(dynviol[1]);
		merit += penalty_coeff*wrapAngle(fabs(dynviol[2])); // since angles wrap
	}
	return merit;
}

void pythonDisplaySmoothTrajectory(const std::vector< Matrix<C_DIM> >& X, std::vector< Matrix<U_DIM> >& U, bool integrate, bool pause, int timesteps = T_SMOOTH) {
	Matrix<X_DIM> x;
	x.insert(0, 0, X[0]);
	for(int i = 0; i < NUM_LANDMARKS; ++i) {
		x.insert(C_DIM+2*i, 0, landmarks[i]);
	}

	std::vector<Matrix<X_DIM> > X_full(timesteps, x);
	X_full[0] = x;
	for(int t=0; t < timesteps-1; ++t) {
		if (integrate) {
			x.insert(0, 0, dynfunccar(x.subMatrix<C_DIM,1>(0,0), U[t]));
			X_full[t+1] = x;
		} else{
			X_full[t+1].insert(0, 0, X[t]);
		}
	}

	pythonDisplayTrajectory(X_full, timesteps, pause);
}


bool isValidSmoothInputs()
{
	for(int t = 0; t < T_SMOOTH-1; ++t) {
		for(int i=0; i < (3*C_DIM+U_DIM); ++i) { if (H_smooth[t][i] > INFTY/2) { return false; } }
		for(int i=0; i < (3*C_DIM+U_DIM); ++i) { if (f_smooth[t][i] > INFTY/2) { return false; } }
		for(int i=0; i < (3*C_DIM+U_DIM); ++i) { if (lb_smooth[t][i] > INFTY/2) { return false; } }
		for(int i=0; i < (C_DIM+U_DIM); ++i) {if (lb_smooth[t][i] > INFTY/2) { return false; } }
		for(int i=0; i < (C_DIM*(3*C_DIM+U_DIM)); ++i) { if (C_smooth[t][i] > INFTY/2) { return false; } }
		for(int i=0; i < C_DIM; ++i) { if (e_smooth[t][i] > INFTY/2) { return false; } }
	}
	for(int i=0; i < (C_DIM); ++i) { if (H_smooth[T_SMOOTH-1][i] > INFTY/2) { return false; } }
	for(int i=0; i < (C_DIM); ++i) { if (f_smooth[T_SMOOTH-1][i] > INFTY/2) { return false; } }
	for(int i=0; i < (C_DIM); ++i) { if (lb_smooth[T_SMOOTH-1][i] > INFTY/2) { return false; } }
	for(int i=0; i < (C_DIM); ++i) { if (ub_smooth[T_SMOOTH-1][i] > INFTY/2) { return false; } }
	for(int i=0; i < C_DIM; ++i) { if (e_smooth[T_SMOOTH-1][i] > INFTY/2) { return false; } }

	return true;

	for(int t = 0; t < T_SMOOTH-1; ++t) {

		std::cout << "t: " << t << std::endl << std::endl;

		std::cout << "f_smooth: ";
		for(int i = 0; i < (3*C_DIM+U_DIM); ++i) {
			std::cout << f_smooth[t][i] << " ";
		}
		std::cout << std::endl;

		std::cout << "lb_smooth c: ";
		for(int i = 0; i < C_DIM; ++i) {
			std::cout << lb_smooth[t][i] << " ";
		}
		std::cout << std::endl;

		std::cout << "ub_smooth c: ";
		for(int i = 0; i < C_DIM; ++i) {
			std::cout << ub_smooth[t][i] << " ";
		}
		std::cout << std::endl;


		std::cout << "lb_smooth s, t: ";
		for(int i = 0; i < 2*C_DIM; ++i) {
			std::cout << lb_smooth[t][C_DIM+U_DIM+i] << " ";
		}
		std::cout << std::endl;


		std::cout << "lb_smooth u: ";
		for(int i = 0; i < U_DIM; ++i) {
			std::cout << lb_smooth[t][C_DIM+i] << " ";
		}
		std::cout << std::endl;

		std::cout << "ub_smooth u: ";
		for(int i = 0; i < U_DIM; ++i) {
			std::cout << ub_smooth[t][C_DIM+i] << " ";
		}
		std::cout << std::endl;


		/*
		std::cout << "C_smooth:" << std::endl;
		if (t == 0) {
			for(int i = 0; i < 2*B_DIM*(3*B_DIM+U_DIM); ++i) {
				std::cout << C_smooth[t][i] << " ";
			}
		} else {
			for(int i = 0; i < B_DIM*(3*B_DIM+U_DIM); ++i) {
				std::cout << C_smooth[t][i] << " ";
			}
		}
		std::cout << std::endl;

		std::cout << "e_smooth:" << std::endl;
		if (t == 0) {
			for(int i = 0; i < 2*B_DIM; ++i) {
				std::cout << e_smooth[t][i] << " ";
			}
		} else {
			for(int i = 0; i < B_DIM; ++i) {
				std::cout << e_smooth[t][i] << " ";
			}
		}


		std::cout << "e_smooth:" << std::endl;
		for(int i = 0; i < C_DIM; ++i) {
			std::cout << e_smooth[t][i] << " ";
		}

		std::cout << std::endl << std::endl;
		 */
	}

	/*
	std::cout << "e_smooth:" << std::endl;
	for(int i = 0; i < C_DIM; ++i) {
		std::cout << e_smooth[T-1][i] << " ";
	}
	std::cout << std::endl << std::endl;
	*/

	std::cout << "t: " << T_SMOOTH-1 << std::endl << std::endl;

	std::cout << "lb_smooth c: ";
	for(int i = 0; i < C_DIM; ++i) {
		std::cout << lb_smooth[T_SMOOTH-1][i] << " ";
	}
	std::cout << std::endl;

	std::cout << "ub_smooth c: ";
	for(int i = 0; i < C_DIM; ++i) {
		std::cout << ub_smooth[T_SMOOTH-1][i] << " ";
	}
	std::cout << std::endl;
	std::cout << std::endl;
	return true;
}



bool minimizeMeritFunction(std::vector< Matrix<C_DIM> >& X, std::vector< Matrix<U_DIM> >& U, const std::vector< Matrix<C_DIM> >& X_unsmooth, smoothMPC_params &problem, smoothMPC_output &output, smoothMPC_info &info, double penalty_coeff)
{
	LOG_DEBUG("Solving sqp problem with penalty parameter: %2.4f", penalty_coeff);

	int T_unsmooth = X_unsmooth.size();
	int timestep_ratio = (int)T_SMOOTH / T_unsmooth;

	std::vector< Matrix<C_DIM,C_DIM> > F(T_SMOOTH-1);
	std::vector< Matrix<C_DIM,U_DIM> > G(T_SMOOTH-1);
	std::vector< Matrix<C_DIM> > h(T_SMOOTH-1);

	double Xpos_eps = cfg::initial_Xpos_trust_box_size;
	double Xangle_eps = cfg::initial_Xangle_trust_box_size;
	double Uvel_eps = cfg::initial_Uvel_trust_box_size;
	double Uangle_eps = cfg::initial_Uangle_trust_box_size;

	double optcost = 0, cost = 0;

	std::vector<Matrix<C_DIM> > Xopt(T_SMOOTH);
	std::vector<Matrix<U_DIM> > Uopt(T_SMOOTH-1);

	double merit, model_merit, new_merit;
	double approx_merit_improve, exact_merit_improve, merit_improve_ratio;

	double hessian_constant;

	int sqp_iter = 1, index = 0;
	bool success = false;

	Matrix<C_DIM,C_DIM> IB = identity<C_DIM>();
	Matrix<C_DIM,C_DIM> minusIB = -identity<C_DIM>();
	
	Matrix<C_DIM,3*C_DIM+U_DIM> CMat;
	Matrix<C_DIM> eVec;

	// sqp loop
	while(true)
	{
		// In this loop, we repeatedly construct a linear approximation to the nonlinear belief dynamics constraint
		LOG_DEBUG("  sqp iter: %d", sqp_iter);

		merit = computeSmoothMerit(X, U, X_unsmooth, penalty_coeff);
		cost = computeSmoothCost(X, U, X_unsmooth);
		
		LOG_DEBUG("  merit: %4.10f", merit);

		// Problem linearization and definition
		// fill in f_smooth, C_smooth, e_smooth
		
		hessian_constant = 0;

		for (int t = 0; t < T_SMOOTH-1; ++t)
		{
			Matrix<C_DIM>& xt = X[t];
			Matrix<U_DIM>& ut = U[t];
			Matrix<(C_DIM+U_DIM)> zbar;

			linearizeCarDynamicsSmooth(xt, ut, F[t], G[t], h[t]);

			index = 0;
			for(int i=0; i < C_DIM; ++i) { H_smooth[t][index++] = 0; }
			for(int i=0; i < U_DIM; ++i) { H_smooth[t][index++] = 2*alpha_control; }
			for(int i=0; i < 2*C_DIM; ++i) { H_smooth[t][index++] = 0; } // TODO: add constant?

			// initialize f_smooth in cost function to penalize
			// belief dynamics slack variables
			index = 0;
			for(int i = 0; i < (C_DIM+U_DIM); ++i) { f_smooth[t][index++] = 0; }
			for(int i = 0; i < 2*C_DIM; ++i) { f_smooth[t][index++] = penalty_coeff; }

			if ((t % timestep_ratio) == 0) {
				for(int i=0; i < P_DIM; ++i) { H_smooth[t][i] = 2*alpha_midway; }
				for(int i=0; i < P_DIM; ++i) { f_smooth[t][i] = -2*alpha_midway*X_unsmooth[t / timestep_ratio][i]; }
			}

			zbar.insert(0, 0, xt);
			zbar.insert(C_DIM, 0, ut);

			for(int i = 0; i < (C_DIM+U_DIM); ++i) {
				hessian_constant += H_smooth[t][i]*zbar[i]*zbar[i];
			}

			CMat.reset();
			eVec.reset();

			CMat.insert<C_DIM,C_DIM>(0,0,F[t]);
			CMat.insert<C_DIM,U_DIM>(0,C_DIM,G[t]);
			CMat.insert<C_DIM,C_DIM>(0,C_DIM+U_DIM,IB);
			CMat.insert<C_DIM,C_DIM>(0,2*C_DIM+U_DIM,minusIB);

			fillColMajor(C_smooth[t], CMat);

			if (t == 0) {
				eVec.insert<C_DIM,1>(0,0,X[0]);
				fillCol(e_smooth[0], eVec);
			} 
			
			eVec = -h[t] + F[t]*xt + G[t]*ut;
			fillCol(e_smooth[t+1], eVec);

		}

		for(int i=0; i < C_DIM; ++i) { H_smooth[T_SMOOTH-1][i] = 2*alpha_goal; }
		for(int i=0; i < C_DIM; ++i) { f_smooth[T_SMOOTH-1][i] = -2*alpha_goal*X_unsmooth[T_unsmooth-1][i]; }

		for(int i = 0; i < (C_DIM); ++i) {
			hessian_constant += H_smooth[T_SMOOTH-1][i]*X[T_SMOOTH-1][i]*X[T_SMOOTH-1][i];
		}


		//std::cout << "PAUSED INSIDE MINIMIZEMERITFUNCTION" << std::endl;
		//int k;
		//std::cin >> k;


		// trust region size adjustment
		while(true)
		{
			LOG_DEBUG("       trust region size: %2.6f %2.6f %2.6f %2.6f", Xpos_eps, Xangle_eps, Uvel_eps, Uangle_eps);

			util::Timer fillVarsTimer;
			util::Timer_tic(&fillVarsTimer);
			// solve the innermost QP here
			for(int t = 0; t < T_SMOOTH-1; ++t)
			{
				Matrix<C_DIM>& xt = X[t];
				Matrix<U_DIM>& ut = U[t];

				// Fill in lb_smooth, ub_smooth

				index = 0;
				// x pos lower bound
				for(int i = 0; i < P_DIM; ++i) { lb_smooth[t][index++] = MAX(xMin[i], xt[i] - Xpos_eps); }
				// x angle lower bound
				lb_smooth[t][index++] = MAX(xMin[P_DIM], xt[P_DIM] - Xangle_eps);
				// u lower bound
				lb_smooth[t][index++] = MAX(uMinSmooth[0], ut[0] - Uvel_eps);
				lb_smooth[t][index++] = MAX(uMinSmooth[1], ut[1] - Uangle_eps);
				//for(int i = 0; i < U_DIM; ++i) { lb_smooth[t][index++] = MAX(uMinSmooth[i], ut[i] - Ueps); }

				// for lower bound on L1 slacks
				for(int i = 0; i < 2*C_DIM; ++i) { lb_smooth[t][index++] = 0; }

				index = 0;
				// x pos upper bound
				for(int i = 0; i < P_DIM; ++i) { ub_smooth[t][index++] = MIN(xMax[i], xt[i] + Xpos_eps); }
				// x angle upper bound
				ub_smooth[t][index++] = MIN(xMax[P_DIM], xt[P_DIM] + Xangle_eps);
				// u upper bound
				ub_smooth[t][index++] = MIN(uMaxSmooth[0], ut[0] + Uvel_eps);
				ub_smooth[t][index++] = MIN(uMaxSmooth[1], ut[1] + Uangle_eps);
				//for(int i = 0; i < U_DIM; ++i) { ub_smooth[t][index++] = MIN(uMaxSmooth[i], ut[i] + Ueps); }

				//for(int i = 0; i < 2*B_DIM; ++i) { ub_smooth[t][index++] = INFTY; }

			}

			double posDelta = .1;
			double angleDelta = M_PI/10;

			for(int t=0; t < T_unsmooth-1; ++t) {
				// x tight pos lower bound
				for(int i=0; i < P_DIM; ++i) { lb_smooth[t*timestep_ratio][i] = X_unsmooth[t][i] - posDelta; }

				// x tight pos upper bound
				for(int i=0; i < P_DIM; ++i) { ub_smooth[t*timestep_ratio][i] = X_unsmooth[t][i] + posDelta; }
			}


			index = 0;
			// lower bound last timestep. tight bounds on everything
			for(int i = 0; i < P_DIM; ++i) { lb_smooth[T_SMOOTH-1][index++] = X_unsmooth[T_unsmooth-1][i] - posDelta; }
			lb_smooth[T_SMOOTH-1][index++] = X_unsmooth[T_unsmooth-1][2] - angleDelta;

			index = 0;
			// upper bound last timestep. tight bounds on everything
			for(int i = 0; i < P_DIM; ++i) { ub_smooth[T_SMOOTH-1][index++] = X_unsmooth[T_unsmooth-1][i] + posDelta; }
			ub_smooth[T_SMOOTH-1][index++] = X_unsmooth[T_unsmooth-1][2] + angleDelta;

			double constant_cost = 0.5*hessian_constant + cost;

			LOG_DEBUG("  original cost: %4.10f", cost);
			LOG_DEBUG("  hessian cost: %4.10f", 0.5*hessian_constant);
			LOG_DEBUG("  constant cost: %4.10f", constant_cost);

			// Verify problem inputs
			//if (!isValidSmoothInputs()) {
			//	LOG_ERROR("Inputs are not valid!");
			//	exit(-1);
			//}

			int exitflag = smoothMPC_solve(&problem, &output, &info);
			if (exitflag == 1) {
				optcost = info.pobj;
				for(int t = 0; t < T_SMOOTH-1; ++t) {
					Matrix<C_DIM>& xt = Xopt[t];
					Matrix<U_DIM>& ut = Uopt[t];

					for(int i = 0; i < C_DIM; ++i) {
						xt[i] = z_smooth[t][i];
					}
					for(int i = 0; i < U_DIM; ++i) {
						ut[i] = z_smooth[t][C_DIM+i];
					}
				}
				for(int i = 0; i < C_DIM; ++i) {
					Xopt[T_SMOOTH-1][i] = z_smooth[T_SMOOTH-1][i];
				}
			}
			else {
				LOG_ERROR("Some problem in smooth solver");
				throw forces_exception(-1);
			}


			LOG_DEBUG("Optimized cost: %4.10f", optcost);

			model_merit = optcost + constant_cost;
			new_merit = computeSmoothMerit(Xopt, Uopt, X_unsmooth, penalty_coeff);

			LOG_DEBUG("merit: %4.10f", merit);
			LOG_DEBUG("model_merit: %4.10f", model_merit);
			LOG_DEBUG("new_merit: %4.10f", new_merit);
			
			double new_cost = computeSmoothCost(Xopt, Uopt, X_unsmooth);
			LOG_DEBUG("new_cost: %4.10f", new_cost);
			LOG_DEBUG("dynviol: %4.10f", new_merit-new_cost);
			double deviation = deviationCost(Xopt, Uopt, X_unsmooth);
			LOG_DEBUG("deviationCost: %4.10f", deviation);

			approx_merit_improve = merit - model_merit;
			exact_merit_improve = merit - new_merit;
			merit_improve_ratio = exact_merit_improve / approx_merit_improve;

			LOG_DEBUG("approx_merit_improve: %1.6f", approx_merit_improve);
			LOG_DEBUG("exact_merit_improve: %1.6f", exact_merit_improve);
			LOG_DEBUG("merit_improve_ratio: %1.6f", merit_improve_ratio);
			
			std::cout << "Paused in minimizeMeritFunction\n";
			std::cin.ignore();

			if (approx_merit_improve < -1) { // -1e-3
				LOG_ERROR("Approximate merit function got worse: %1.6f", approx_merit_improve);
				LOG_ERROR("Either convexification is wrong to zeroth order, or you are in numerical trouble");
				LOG_ERROR("Failure!");

				//for(int t=0; t < T_SMOOTH-1; ++t) {
				//	std::cout << ~U[t];
				//}
				//throw forces_exception(-1); // TODO: keep?
				pythonDisplaySmoothTrajectory(X, U, true, true);

				return false;
			} else if (approx_merit_improve < cfg::min_approx_improve) {
				LOG_DEBUG("Converged: improvement small enough");
				X = Xopt; U = Uopt;
				return true;
			} else if ((exact_merit_improve < 0) || (merit_improve_ratio < cfg::improve_ratio_threshold)) {
				Xpos_eps *= cfg::trust_shrink_ratio;
				Xangle_eps *= cfg::trust_shrink_ratio;
				Uvel_eps *= cfg::trust_shrink_ratio;
				Uangle_eps *= cfg::trust_shrink_ratio;
				LOG_DEBUG("Shrinking trust region size to: %2.6f %2.6f %2.6f %2.6f", Xpos_eps, Xangle_eps, Uvel_eps, Uangle_eps);
			} else {
				Xpos_eps *= cfg::trust_expand_ratio;
				Xangle_eps *= cfg::trust_expand_ratio;
				Uvel_eps *= cfg::trust_expand_ratio;
				Uangle_eps *= cfg::trust_expand_ratio;
				X = Xopt; U = Uopt;
				LOG_DEBUG("Accepted, Increasing trust region size to:  %2.6f %2.6f %2.6f %2.6f", Xpos_eps, Xangle_eps, Uvel_eps, Uangle_eps);
				break;
			}

			if (Xpos_eps < cfg::min_trust_box_size && Xangle_eps < cfg::min_trust_box_size &&
					Uvel_eps < cfg::min_trust_box_size && Uangle_eps < cfg::min_trust_box_size) {
				LOG_DEBUG("Converged: x tolerance");
				X = Xopt; U = Uopt; // TODO: added
				return true;
			}


		} // trust region loop
		sqp_iter++;
	} // sqp loop

	return success;
}

double smoothCollocation(std::vector< Matrix<C_DIM> >& X, std::vector< Matrix<U_DIM> >& U, const std::vector< Matrix<C_DIM> >& X_unsmooth, smoothMPC_params &problem, smoothMPC_output &output, smoothMPC_info &info)
{
	double penalty_coeff = cfg::initial_penalty_coeff;

	int penalty_increases = 0;

	Matrix<C_DIM> dynviol;

	// penalty loop
	while(penalty_increases < cfg::max_penalty_coeff_increases)
	{
		bool success = minimizeMeritFunction(X, U, X_unsmooth, problem, output, info, penalty_coeff);

		double cntviol = 0;
		for(int t = 0; t < T_SMOOTH-1; ++t) {
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
	    }
	    else {
	    	return computeSmoothCost(X, U, X_unsmooth);
	    }
	}
	return computeSmoothCost(X, U, X_unsmooth);
}

std::vector<Matrix<U_DIM> > smoothTraj(const std::vector<Matrix<C_DIM> >& X_unsmooth, const std::vector<Matrix<U_DIM> >& U_unsmooth) {
	int unsmooth_timesteps = X_unsmooth.size();

	uMinSmooth[0] = .01;
	uMinSmooth[1] = -M_PI/4;
	uMaxSmooth[0] = 10;
	uMaxSmooth[1] = M_PI/4;

	smoothMPC_params problem;
	smoothMPC_output output;
	smoothMPC_info info;
	setupSmoothVars(problem, output);

	std::vector<Matrix<C_DIM> > X(T_SMOOTH);

	int timestep_ratio = (int)T_SMOOTH / unsmooth_timesteps;

	int index = 0;
	for(int t=0; t < T-1; ++t) {
		for(int i=0; i < timestep_ratio; ++i) {
			//float next_weight = (i/((float)timestep_ratio));
			//X[index++] = (1-next_weight)*X_unsmooth[t] + next_weight*X_unsmooth[t+1];
			X[index++] = X_unsmooth[t]; // repeat X_unsmooth
		}
	}
	for(int i=index; i < T_SMOOTH; ++i) {
		X[i] = X_unsmooth[unsmooth_timesteps-1];
	}

	Matrix<U_DIM> uinit;
	uinit[0] = uMinSmooth[0];
	uinit[1] = 0;
	std::vector<Matrix<U_DIM> > U(T_SMOOTH-1, uinit);

	// U_unsmooth every timestep_ratio
	for(int t=0; t < T_SMOOTH-1; ++t) {
		if ((t % timestep_ratio) == 0) {
			U[t] = U_unsmooth[(int)t / timestep_ratio];
		}
	}

	std::cout << "X:\n";
	for(int t=0; t < T_SMOOTH; ++t) {
		std::cout << ~X[t];
	}

	std::cout << "\nU:\n";
	for(int t=0; t < T_SMOOTH-1; ++t) {
		std::cout << ~U[t];
	}


	util::Timer solveTimer;
	util::Timer_tic(&solveTimer);


	try {
		smoothCollocation(X, U, X_unsmooth, problem, output, info);
	} catch (forces_exception &e) {
		LOG_ERROR("Smooth solver failed!");
		pythonDisplaySmoothTrajectory(X, U, true, true);
		exit(-1);
	}


//	int iter = 0;
//	while(true) {
//		try {
//			smoothCollocation(X, U, X_unsmooth, problem, output, info);
//			break;
//		}
//		catch (forces_exception &e) {
//			if (iter > 3) {
//				LOG_ERROR("Tried too many times, giving up");
//				pythonDisplaySmoothTrajectory(X, U, true, true);
//				exit(-1);
//			}
//			LOG_ERROR("Forces exception, trying again");
//			X[0] = X_unsmooth[0];
//			for(int j=0; j < T_SMOOTH-1; ++j) {
//				std::cout << ~X[j];
//				X[j+1] = dynfunccar(X[j], U[j]);
//			}
//			std::cout << ~X[T_SMOOTH-1];
//			pythonDisplaySmoothTrajectory(X, U, true, true);
//			iter++;
//		}
//	}


	double solveTime = util::Timer_toc(&solveTimer);

	std::cout << "solveTime: " << solveTime*1000 << "ms\n";


	std::cout << "X[0]: " << ~X[0];
	std::cout << "X_unsmooth[0]: " << ~X_unsmooth[0];

	double deviationBefore = deviationCost(X, U, X_unsmooth);
	std::cout << "deviationCost before integration: " << deviationBefore << "\n";

	X[0] = X_unsmooth[0];
	for(int t=0; t < T_SMOOTH-1; ++t) {
		X[t+1] = dynfunccar(X[t], U[t]);
	}


	double cost = computeSmoothCost(X, U, X_unsmooth);
	std::cout << "Final cost: " << cost << "\n";

	double diffCost = cost;
	for(int t=0; t < T_SMOOTH-1; ++t) {
		diffCost -= alpha_control*tr(~U[t]*U[t]);
	}
	std::cout << "Difference between smoothed and unsmoothed: " << diffCost << "\n";

	double deviation = deviationCost(X, U, X_unsmooth);
	std::cout << "deviationCost after integration: " << deviation << "\n";

	pythonDisplaySmoothTrajectory(X, U, false, true);


	return U;
}


/*
int main(int argc, char* argv[])
{
	initProblemParams();

	std::vector<Matrix<C_DIM> > X(T, zeros<C_DIM,1>());
	std::vector<Matrix<U_DIM> > U(T-1, zeros<U_DIM,1>());


//	U[0][0] = 2.29551;      U[0][1] = 0.056562;
//	U[1][0] = 7.10411;      U[1][1] = 0.40421;
//	U[2][0] = 7.45238;      U[2][1] = -0.178322;
//	U[3][0] = 6.15177;      U[3][1] = -0.8003;
//	U[4][0] = 5.0965;       U[4][1] = -0.236125;
//	U[5][0] = 4.9362;       U[5][1] = 0.172297;
//	U[6][0] = 1.92251;      U[6][1] = -0.0218368;
//	U[7][0] = 2.3931;       U[7][1] = 0.526903;
//	U[8][0] = 1.86502;      U[8][1] = 0.556067;
//	U[9][0] = 1.75637;      U[9][1] = 0.322271;
//	U[10][0] = 2.71267;     U[10][1] = 0.254822;
//	U[11][0] = 3.30179;     U[11][1] = 0.122705;
//	U[12][0] = 9.99999;     U[12][1] = -0.0724144;
//	U[13][0] = 10;          U[13][1] = 0.329621;
//
//	X[0][0] = 0;
//	X[0][1] = 0;
//	X[0][2] = 0;

//	U[0][0] = 1.82213;       U[0][1] = 1.0472;
//	U[1][0] = 1.67913;       U[1][1] = 1.0472;
//	U[2][0] = 8.06062;       U[2][1] = -0.846693;
//	U[3][0] = 1.5;           U[3][1] = 0.0743642;
//	U[4][0] = 1.5;           U[4][1] = 0.52193;
//	U[5][0] = 1.5;           U[5][1] = 1.0472;
//	U[6][0] = 1.5;           U[6][1] = 1.0472;
//	U[7][0] = 1.5;           U[7][1] = 0.562566;
//	U[8][0] = 1.5;           U[8][1] = 0.358109;
//	U[9][0] = 1.5;           U[9][1] = 1.0472;
//	U[10][0] = 1.5;          U[10][1] = 1.0472;
//	U[11][0] = 1.5;          U[11][1] = 1.0472;
//	U[12][0] = 1.5;          U[12][1] = 0.967156;
//	U[13][0] = 3.61221;      U[13][1] = -0.172981;
//
//	X[0][0] = 59.7571;
//	X[0][1] = -0.0704439;
//	X[0][2] = 0.786252;

//	U[0][0] = 1.5;           U[0][1] = 0.646742;
//	U[1][0] = 1.5;           U[1][1] = 0.64126;
//	U[2][0] = 1.5;           U[2][1] = 0.579988;
//	U[3][0] = 2.18841;       U[3][1] = 0.75854;
//	U[4][0] = 1.5;           U[4][1] = 0.287933;
//	U[5][0] = 5.76078;       U[5][1] = 0.186944;
//	U[6][0] = 10;            U[6][1] = -0.46209;
//	U[7][0] = 10;            U[7][1] = 0.559336;
//	U[8][0] = 3.44668;       U[8][1] = -0.461086;
//	U[9][0] = 3.78677;       U[9][1] = -0.563789;
//	U[10][0] = 2.41875;      U[10][1] = 0.172149;
//	U[11][0] = 1.76995;      U[11][1] = -0.0564648;
//	U[12][0] = 8.97079;      U[12][1] = -0.676349;
//	U[13][0] = 10;           U[13][1] = 1.0472;
//
//	X[0][0] = 60.1005;
//	X[0][1] = 19.8891;
//	X[0][2] = 2.35812;


	U[0][0] = 4.34473;       U[0][1] = 1.0472;
	U[1][0] = 5.04883;       U[1][1] = 1.03323;
	U[2][0] = 4.16264;       U[2][1] = 0.225835;
	U[3][0] = 3.15039;       U[3][1] = -0.00538109;
	U[4][0] = 1.84172;       U[4][1] = 0.0459601;
	U[5][0] = 1.5;           U[5][1] = 0.0699177;
	U[6][0] = 1.5;           U[6][1] = -0.261927;
	U[7][0] = 1.5;           U[7][1] = -1.0472;
	U[8][0] = 1.5;           U[8][1] = -1.0472;
	U[9][0] = 1.5;           U[9][1] = -1.0472;
	U[10][0] = 1.5;          U[10][1] =  -1.0472;
	U[11][0] = 1.5;          U[11][1] =  -0.794009;
	U[12][0] = 7.27962;      U[12][1] =  -0.991759;
	U[13][0] = 10;           U[13][1] =  0.35943;

	X[0][0] = 0.1;
	X[0][1] = 19.9;
	X[0][2] = 3.927;

	for(int t=0; t < T-1; ++t) {
		std::cout << ~X[t];
		X[t+1] = dynfunccar(X[t],U[t]);
	}
	std::cout << ~X[T-1];

	std::cout << "initial unsmoothed trajectory\n";
	pythonDisplaySmoothTrajectory(X, U, false, true, T);

	smoothTraj(X, U);

	cleanupSmoothMPCVars();
	
}
*/
