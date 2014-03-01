#include "../slam.h"
#include "../traj/slam-traj.h"


#include <vector>
#include <iomanip>

#include "util/matrix.h"
#include "util/Timer.h"

extern "C" {
#include "controlMPC.h"
controlMPC_FLOAT **H, **f, **lb, **ub, **z;
#include "slam-control-casadi.h"
}

#include "boost/preprocessor.hpp"

int number = 0;

const double alpha_belief = 10; // 10;
const double alpha_final_belief = 10; // 10;
const double alpha_control = .1; // .1
const double alpha_goal_state = 10; // 10


namespace cfg {
const double improve_ratio_threshold = .1; // .1
const double min_approx_improve = 1e-2; // 1e-2
const double min_trust_box_size = 1e-3; // 1e-3

const double trust_shrink_ratio = .75; // .5
const double trust_expand_ratio = 1.25; // 1.25

const double cnt_tolerance = 1; // 1e-2
const double penalty_coeff_increase_ratio = 5; // 5
const double initial_penalty_coeff = 10; // 10

const double initial_trust_box_size = .1; // 5 // split up trust box size for X and U
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

	const double **casadi_input = new const double*[4];
	casadi_input[0] = U_arr;
	casadi_input[1] = x0_arr;
	casadi_input[2] = Sigma0_arr;
	casadi_input[3] = xGoal_arr;
	casadi_input[4] = params_arr;

	double cost = 0;
	double **cost_arr = new double*[1];
	cost_arr[0] = &cost;

	evaluateCostWrap(casadi_input, cost_arr);

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

	const double **casadi_input = new const double*[4];
	casadi_input[0] = U_arr;
	casadi_input[1] = x0_arr;
	casadi_input[2] = Sigma0_arr;
	casadi_input[3] = xGoal_arr;
	casadi_input[4] = params_arr;

	double **costgrad_arr = new double*[2];
	costgrad_arr[0] = &cost;
	costgrad_arr[1] = Grad.getPtr();

	evaluateCostGradWrap(casadi_input, costgrad_arr);

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
	bool success;

	Matrix<C_DIM+U_DIM, C_DIM+U_DIM> HMat;
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
		casadiComputeCostGrad(U, cost, Grad);

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

		Matrix<C_DIM,1> cGoal = xGoal.subMatrix<C_DIM,1>(0,0);
		//double goal_cost = -alpha_goal_state*tr(~cGoal*cGoal);
		double goal_cost = 0;

		constant_cost = 0.5*hessian_constant + jac_constant + cost + goal_cost;
		LOG_DEBUG("  hessian cost: %4.10f", 0.5*hessian_constant);
		LOG_DEBUG("  jacobian cost: %4.10f", jac_constant);
		LOG_DEBUG("  goal cost: %4.10f", goal_cost);
		LOG_DEBUG("  constant cost: %4.10f", constant_cost);

		//std::cout << "PAUSED INSIDE MINIMIZEMERITFUNCTION" << std::endl;
		//int k;
		//std::cin >> k;


		// trust region size adjustment
		while(true)
		{
			LOG_DEBUG("       trust region size: %2.6f", Ueps);

			// solve the innermost QP here
			for(int t = 0; t < T-1; ++t)
			{
				Matrix<U_DIM>& ut = U[t];

				// Fill in lb, ub

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


			LOG_DEBUG("       Optimized cost: %4.10f", optcost);

			model_merit = optcost + constant_cost;

			new_merit = casadiComputeCost(Uopt);

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

			//std::cout << "PAUSED INSIDE minimizeMeritFunction AFTER OPTIMIZATION" << std::endl;
			//int num;
			//std::cin >> num;

			if (number > 2) {
				//std::cout << "Displaying Uopt\n";
				//pythonDisplayTrajectory(Uopt, T, false);
			}

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


				casadiComputeCostGrad(Uopt, cost, Gradopt);

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

	return success;
}



void testCasadi() {
	initProblemParams();

	Matrix<U_DIM> uinit;

	xGoal = x0;
	xGoal.insert(0, 0, waypoints[0]);

	std::cout << "x0: " << ~x0;
	std::cout << "xGoal: " << ~xGoal << "\n";

	// initialize velocity to dist / timesteps
	uinit[0] = sqrt((x0[0] - xGoal[0])*(x0[0] - xGoal[0]) + (x0[1] - xGoal[1])*(x0[1] - xGoal[1])) / (double)((T-1)*DT);
	// angle already pointed at goal, so is 0
	uinit[1] = 0;

	std::vector<Matrix<U_DIM> > U(T-1, uinit);

	double initTrajCost = computeCost(U);
	LOG_INFO("Initial trajectory cost: %4.10f", initTrajCost);

	double initCasadiTrajCost = casadiComputeCost(U);
	LOG_INFO("Initial casadi trajectory cost: %4.10f", initCasadiTrajCost);


}


int main(int argc, char* argv[])
{

	LOG_INFO("Initializing problem parameters");
	initProblemParams();

	controlMPC_params problem;
	controlMPC_output output;
	controlMPC_info info;
	setupControlVars(problem, output);

	util::Timer solveTimer, trajTimer;
	double totalSolveTime = 0, trajTime = 0;

	double totalTrajCost = 0;

	int num_loops = 1;

	std::vector<Matrix<B_DIM> > B_total(T*NUM_WAYPOINTS*num_loops);
	std::vector<Matrix<U_DIM> > U_total((T-1)*NUM_WAYPOINTS*num_loops);
	int B_total_idx = 0, U_total_idx = 0;

	std::vector<Matrix<B_DIM> > B(T);

	Matrix<U_DIM> uinit;

	Matrix<X_DIM,1> x;
	Matrix<X_DIM,X_DIM> s;
	for(int loop=0; loop < num_loops; loop++) {
		LOG_INFO("Loop %d",loop);
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
				LOG_ERROR("Failed to initialize trajectory, exiting slam-state");
				exit(-1);
			}

			double initTrajTime = util::Timer_toc(&trajTimer);
			trajTime += initTrajTime;

			vec(x0, SqrtSigma0, B[0]);
			for(int t=0; t < T-1; ++t) {
				B[t+1] = beliefDynamics(B[t], U[t]);
			}

			//std::cout << ~X[0].subMatrix<C_DIM,1>(0,0);

			double initTrajCost = computeCost(U);
			//LOG_INFO("Initial trajectory cost: %4.10f", initTrajCost);

			double initCasadiTrajCost = casadiComputeCost(U);
			//LOG_INFO("Initial casadi trajectory cost: %4.10f", initCasadiTrajCost);

			//pythonDisplayTrajectory(B, U, waypoints, landmarks, T, true);

			util::Timer_tic(&solveTimer);

			double cost = 0;
			int iter = 0;
			while(true) {
				try {
					cost = controlCollocation(U, problem, output, info);
					break;
				}
				catch (forces_exception &e) {
					if (iter > 3) {
						LOG_ERROR("Tried too many times, giving up");
						pythonDisplayTrajectory(U, T, true);
						exit(-1);
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

			std::cout << SqrtSigma0 << "\n";

			pythonDisplayTrajectory(B, U, waypoints, landmarks, T, true);

			number++;

		}


	}



	LOG_INFO("Total trajectory cost: %4.10f", totalTrajCost);
	LOG_INFO("Total trajectory solve time: %5.3f ms", trajTime*1000);
	LOG_INFO("Total solve time: %5.3f ms", totalSolveTime*1000);


	pythonDisplayTrajectory(B_total, U_total, waypoints, landmarks, T*NUM_WAYPOINTS*num_loops, true);

	cleanupControlMPCVars();



	return 0;
}
