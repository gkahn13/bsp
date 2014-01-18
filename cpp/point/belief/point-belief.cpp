#include <vector>
#include <iomanip>

#include "util/matrix.h"
#include "util/utils.h"
#include "util/Timer.h"
#include "util/logging.h"

#include <Python.h>
#include <boost/python.hpp>
#include <boost/filesystem.hpp>

namespace py = boost::python;

#define BELIEF_PENALTY_MPC
//#define BELIEF_MPC

extern "C" {
#ifdef BELIEF_PENALTY_MPC
#include "beliefPenaltyMPC.h"
beliefPenaltyMPC_FLOAT **f, **lb, **ub, **C, **e, **z;
#endif

#ifdef BELIEF_MPC
#include "beliefMPC.h"
beliefMPC_FLOAT **lb, **ub, **C, **e, **z;
#endif
}

#define DT 1.0
#define TIMESTEPS 15
#define X_DIM 2
#define U_DIM 2
#define Z_DIM 2
#define Q_DIM 2
#define R_DIM 2

#define S_DIM (((X_DIM+1)*X_DIM)/2)
#define B_DIM (X_DIM+S_DIM)

const double step = 0.0078125*0.0078125;

Matrix<X_DIM> x0;
Matrix<X_DIM,X_DIM> SqrtSigma0;
Matrix<X_DIM> xGoal;
Matrix<X_DIM> xMin, xMax;
Matrix<U_DIM> uMin, uMax;

const int T = TIMESTEPS;
const double INFTY = 1e10;
const double alpha_belief = 10, alpha_final_belief = 10, alpha_control = 1;

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

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

inline Matrix<X_DIM> dynfunc(const Matrix<X_DIM>& x, const Matrix<U_DIM>& u, const Matrix<U_DIM>& q)
{  
	Matrix<X_DIM> xNew = x + u*DT + 0.01*q;
	return xNew;
}

// Observation model
inline Matrix<Z_DIM> obsfunc(const Matrix<X_DIM>& x, const Matrix<R_DIM>& r)
{
	double intensity = sqrt(sqr(0.5*x[0]) + 1e-6);
	Matrix<Z_DIM> z = x + intensity*r;
	return z;
}

// Jacobians: df(x,u,q)/dx, df(x,u,q)/dq
inline void linearizeDynamics(const Matrix<X_DIM>& x, const Matrix<U_DIM>& u, const Matrix<Q_DIM>& q, Matrix<X_DIM,X_DIM>& A, Matrix<X_DIM,Q_DIM>& M)
{
	A.reset();
	Matrix<X_DIM> xr(x), xl(x);
	for (size_t i = 0; i < X_DIM; ++i) {
		xr[i] += step; xl[i] -= step;
		A.insert(0,i, (dynfunc(xr, u, q) - dynfunc(xl, u, q)) / (xr[i] - xl[i]));
		xr[i] = x[i]; xl[i] = x[i];
	}

	M.reset();
	Matrix<Q_DIM> qr(q), ql(q);
	for (size_t i = 0; i < Q_DIM; ++i) {
		qr[i] += step; ql[i] -= step;
		M.insert(0,i, (dynfunc(x, u, qr) - dynfunc(x, u, ql)) / (qr[i] - ql[i]));
		qr[i] = q[i]; ql[i] = q[i];
	}
}

// Jacobians: dh(x,r)/dx, dh(x,r)/dr
inline void linearizeObservation(const Matrix<X_DIM>& x, const Matrix<R_DIM>& r, Matrix<Z_DIM,X_DIM>& H, Matrix<Z_DIM,R_DIM>& N)
{
	H.reset();
	Matrix<X_DIM> xr(x), xl(x);
	for (size_t i = 0; i < X_DIM; ++i) {
		xr[i] += step; xl[i] -= step;
		H.insert(0,i, (obsfunc(xr, r) - obsfunc(xl, r)) / (xr[i] - xl[i]));
		xr[i] = x[i]; xl[i] = x[i];
	}

	N.reset();
	Matrix<R_DIM> rr(r), rl(r);
	for (size_t i = 0; i < R_DIM; ++i) {
		rr[i] += step; rl[i] -= step;
		N.insert(0,i, (obsfunc(x, rr) - obsfunc(x, rl)) / (rr[i] - rl[i]));
		rr[i] = r[i]; rl[i] = r[i];
	}
}

// Switch between belief vector and matrices
inline void unVec(const Matrix<B_DIM>& b, Matrix<X_DIM>& x, Matrix<X_DIM,X_DIM>& SqrtSigma) {
	x = b.subMatrix<X_DIM,1>(0,0);
	size_t idx = X_DIM;
	for (size_t j = 0; j < X_DIM; ++j) {
		for (size_t i = j; i < X_DIM; ++i) {
			SqrtSigma(i,j) = b[idx];
			SqrtSigma(j,i) = b[idx];
			++idx;
		}
	}
}

inline void vec(const Matrix<X_DIM>& x, const Matrix<X_DIM,X_DIM>& SqrtSigma, Matrix<B_DIM>& b) {
	b.insert(0,0,x);
	size_t idx = X_DIM;
	for (size_t j = 0; j < X_DIM; ++j) {
		for (size_t i = j; i < X_DIM; ++i) {
			b[idx] = 0.5 * (SqrtSigma(i,j) + SqrtSigma(j,i));
			++idx;
		}
	}
}

// Belief dynamics
inline Matrix<B_DIM> beliefDynamics(const Matrix<B_DIM>& b, const Matrix<U_DIM>& u) {
	Matrix<X_DIM> x;
	Matrix<X_DIM,X_DIM> SqrtSigma;
	unVec(b, x, SqrtSigma);

	Matrix<X_DIM,X_DIM> Sigma = SqrtSigma*SqrtSigma;

	Matrix<X_DIM,X_DIM> A;
	Matrix<X_DIM,Q_DIM> M;
	linearizeDynamics(x, u, zeros<Q_DIM,1>(), A, M);

	x = dynfunc(x, u, zeros<Q_DIM,1>());
	Sigma = A*Sigma*~A + M*~M;

	Matrix<Z_DIM,X_DIM> H;
	Matrix<Z_DIM,R_DIM> N;
	linearizeObservation(x, zeros<R_DIM,1>(), H, N);

	Matrix<X_DIM,Z_DIM> K = Sigma*~H/(H*Sigma*~H + N*~N);

	Sigma = (identity<X_DIM>() - K*H)*Sigma;

	Matrix<B_DIM> g;
	vec(x, sqrt(Sigma), g);

	return g;
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

#ifdef BELIEF_MPC
void setupBeliefVars(beliefMPC_params& problem, beliefMPC_output& output)
{
	lb = new beliefMPC_FLOAT*[T];
	ub = new beliefMPC_FLOAT*[T];
	C = new beliefMPC_FLOAT*[T-1];
	e = new beliefMPC_FLOAT*[T-1];
	z = new beliefMPC_FLOAT*[T];

#define SET_VARS(n)    \
		C[ BOOST_PP_SUB(n,1) ] = problem.C##n ;  \
		e[ BOOST_PP_SUB(n,1) ] = problem.e##n ;  \
		lb[ BOOST_PP_SUB(n,1) ] = problem.lb##n ;	\
		ub[ BOOST_PP_SUB(n,1) ] = problem.ub##n ;	\
		z[ BOOST_PP_SUB(n,1) ] = output.z##n ;

#define BOOST_PP_LOCAL_MACRO(n) SET_VARS(n)
#define BOOST_PP_LOCAL_LIMITS (1, TIMESTEPS-1)
#include BOOST_PP_LOCAL_ITERATE()

#endif

#ifdef BELIEF_PENALTY_MPC
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

#endif

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
#ifdef BELIEF_PENALTY_MPC
	delete[] f;
#endif

	delete[] lb; 
	delete[] ub; 
	delete[] C;
	delete[] e;
	delete[] z;
}

#ifdef BELIEF_PENALTY_MPC

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

	int sqp_iter = 1;
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

				// Fill in f, lb, ub, C, e
				for(int i = 0; i < (B_DIM+U_DIM); ++i) {
					f[t][i] = 0;
				}
				for(int i = 0; i < 2*B_DIM; ++i) {
					f[t][B_DIM+U_DIM+i] = penalty_coeff;
				}
				lb[t][0] = MAX(xMin[0], bt[0] - Beps);
				lb[t][1] = MAX(xMin[1], bt[1] - Beps);
				lb[t][2] = bt[2] - Beps;
				lb[t][3] = bt[3] - Beps;
				lb[t][4] = bt[4] - Beps;
				lb[t][5] = MAX(uMin[0], ut[0] - Ueps);
				lb[t][6] = MAX(uMin[1], ut[1] - Ueps);
				for(int i = 0; i < 2*B_DIM; ++i) {
					lb[t][B_DIM+U_DIM+i] = 0;
				}

				ub[t][0] = MIN(xMax[0], bt[0] + Beps);
				ub[t][1] = MIN(xMax[1], bt[1] + Beps);
				ub[t][2] = bt[2] + Beps;
				ub[t][3] = bt[3] + Beps;
				ub[t][4] = bt[4] + Beps;
				ub[t][5] = MIN(uMax[0], ut[0] + Ueps);
				ub[t][6] = MIN(uMax[1], ut[1] + Ueps);
				//for(int i = 0; i < 2*B_DIM; ++i) {
				//	ub[t][B_DIM+U_DIM+i] = INFTY;
				//}

				if (t > 0) {
					Matrix<B_DIM,3*B_DIM+U_DIM> CMat;
					Matrix<B_DIM> eVec;

					CMat.insert<B_DIM,B_DIM>(0,0,F[t]);
					CMat.insert<B_DIM,U_DIM>(0,B_DIM,G[t]);
					CMat.insert<B_DIM,B_DIM>(0,B_DIM+U_DIM,I);
					CMat.insert<B_DIM,B_DIM>(0,2*B_DIM+U_DIM,minusI);

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
					Matrix<2*B_DIM,3*B_DIM+U_DIM> CMat;
					Matrix<2*B_DIM> eVec;

					CMat.insert<B_DIM,B_DIM>(0,0,I);
					CMat.insert<B_DIM,U_DIM>(0,B_DIM,zeros<B_DIM,U_DIM>());
					CMat.insert<B_DIM,2*B_DIM>(0,B_DIM+U_DIM,zeros<B_DIM,2*B_DIM>());

					CMat.insert<B_DIM,B_DIM>(B_DIM,0,F[t]);
					CMat.insert<B_DIM,U_DIM>(B_DIM,B_DIM,G[t]);
					CMat.insert<B_DIM,B_DIM>(B_DIM,B_DIM+U_DIM,I);
					CMat.insert<B_DIM,B_DIM>(B_DIM,2*B_DIM+U_DIM,minusI);

					int idx = 0;
					int nrows = CMat.numRows(), ncols = CMat.numColumns();
					for(int c = 0; c < ncols; ++c) {
						for(int r = 0; r < nrows; ++r) {
							C[t][idx++] = CMat[c + r*ncols];
						}
					}
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
			lb[T-1][0] = MAX(xGoal[0] - delta, bT[0] - Beps);
			lb[T-1][1] = MAX(xGoal[1] - delta, bT[1] - Beps);
			lb[T-1][2] = bT[2] - Beps;
			lb[T-1][3] = bT[3] - Beps;
			lb[T-1][4] = bT[4] - Beps;

			ub[T-1][0] = MIN(xGoal[0] + delta, bT[0] + Beps);
			ub[T-1][1] = MIN(xGoal[1] + delta, bT[1] + Beps);
			ub[T-1][2] = bT[2] + Beps;
			ub[T-1][3] = bT[3] + Beps;
			ub[T-1][4] = bT[4] + Beps;

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
#endif

#ifdef BELIEF_MPC
double beliefCollocation(std::vector< Matrix<B_DIM> >& B, std::vector< Matrix<U_DIM> >& U, beliefMPC_params& problem, beliefMPC_output& output, beliefMPC_info& info)
{
	int maxIter = 10;
	double Beps = 1;
	double Ueps = 1;

	// box constraint around goal
	double delta = 0.01;

	Matrix<B_DIM,1> b0 = B[0];

	std::vector< Matrix<B_DIM,B_DIM> > F(T-1);
	std::vector< Matrix<B_DIM,U_DIM> > G(T-1);
	std::vector< Matrix<B_DIM> > h(T-1);

	double prevcost = computeCost(B, U);
	double optcost;

	//std::cout << "Initialization trajectory cost: " << std::setprecision(10) << prevcost << std::endl;

	for(int it = 0; it < maxIter; ++it) 
	{
		//std::cout << "Iter: " << it << std::endl;

		// linearize belief dynamics constraint here
		for (int t = 0; t < T-1; ++t) 
		{
			Matrix<B_DIM>& bt = B[t];
			Matrix<U_DIM>& ut = U[t];

			linearizeBeliefDynamics(bt, ut, F[t], G[t], h[t]);

			// Fill in lb, ub, C, e
			lb[t][0] = MAX(xMin[0], bt[0] - Beps);
			lb[t][1] = MAX(xMin[1], bt[1] - Beps);
			lb[t][2] = bt[2] - Beps;
			lb[t][3] = bt[3] - Beps;
			lb[t][4] = bt[4] - Beps;
			lb[t][5] = MAX(uMin[0], ut[0] - Ueps);
			lb[t][6] = MAX(uMin[1], ut[1] - Ueps);

			ub[t][0] = MIN(xMax[0], bt[0] + Beps);
			ub[t][1] = MIN(xMax[1], bt[1] + Beps);
			ub[t][2] = bt[2] + Beps;
			ub[t][3] = bt[3] + Beps;
			ub[t][4] = bt[4] + Beps;
			ub[t][5] = MIN(uMax[0], ut[0] + Ueps);
			ub[t][6] = MIN(uMax[1], ut[1] + Ueps);

			if (t > 0) {
				Matrix<B_DIM,B_DIM+U_DIM> CMat;
				Matrix<B_DIM> eVec;

				CMat.insert<B_DIM,B_DIM>(0,0,F[t]);
				CMat.insert<B_DIM,U_DIM>(0,B_DIM,G[t]);
				Matrix<B_DIM+U_DIM, B_DIM> CMatT = ~CMat;
				int idx = 0;
				for(int c = 0; c < (B_DIM+U_DIM); ++c) {
					for(int r = 0; r < B_DIM; ++r) {
						C[t][idx++] = CMat[c + r*(B_DIM+U_DIM)];
					}
				}
				eVec = -h[t] + F[t]*bt + G[t]*ut;
				for(int i = 0; i < B_DIM; ++i) {
					e[t][i] = eVec[i];
				}
			}
			else {
				Matrix<2*B_DIM,B_DIM+U_DIM> CMat;
				Matrix<2*B_DIM> eVec;

				CMat.insert<B_DIM,B_DIM>(0,0,identity<B_DIM>());
				CMat.insert<B_DIM,U_DIM>(0,B_DIM,zeros<B_DIM,U_DIM>());
				CMat.insert<B_DIM,B_DIM>(B_DIM,0,F[t]);
				CMat.insert<B_DIM,U_DIM>(B_DIM,B_DIM,G[t]);
				int idx = 0;
				for(int c = 0; c < (B_DIM+U_DIM); ++c) {
					for(int r = 0; r < 2*B_DIM; ++r) {
						C[t][idx++] = CMat[c + r*(B_DIM+U_DIM)];
					}
				}
				eVec.insert<B_DIM,1>(0,0,bt);
				eVec.insert<B_DIM,1>(B_DIM,0,zeros<B_DIM,1>());
				for(int i = 0; i < 2*B_DIM; ++i) {
					e[t][i] = eVec[i];
				}
			}
		} //setting up problem

		Matrix<B_DIM>& bT = B[T-1];

		// Fill in lb, ub, C, e
		lb[T-1][0] = MAX(xGoal[0] - delta, bT[0] - Beps);
		lb[T-1][1] = MAX(xGoal[1] - delta, bT[1] - Beps);
		lb[T-1][2] = bT[2] - Beps;
		lb[T-1][3] = bT[3] - Beps;
		lb[T-1][4] = bT[4] - Beps;

		ub[T-1][0] = MIN(xGoal[0] + delta, bT[0] + Beps);
		ub[T-1][1] = MIN(xGoal[1] + delta, bT[1] + Beps);
		ub[T-1][2] = bT[2] + Beps;
		ub[T-1][3] = bT[3] + Beps;
		ub[T-1][4] = bT[4] + Beps;

		// Verify problem inputs

		//int num;
		//std::cin >> num;
		
		int exitflag = beliefMPC_solve(&problem, &output, &info);
		if (exitflag == 1) {
			for(int t = 0; t < T-1; ++t) {
				Matrix<B_DIM>& bt = B[t];
				Matrix<U_DIM>& ut = U[t];

				for(int i = 0; i < B_DIM; ++i) {
					bt[i] = z[t][i];
				}
				for(int i = 0; i < U_DIM; ++i) {
					ut[i] = z[t][B_DIM+i];
				}
				optcost = info.pobj;
			}
		}
		else {
			LOG_ERROR("Some problem in solver");
			std::exit(-1);
		}
		LOG_DEBUG("Optimized cost: %4.10f", optcost);

		if ((optcost > prevcost) || (fabs(optcost - prevcost)/prevcost < 0.01))
			break; 
		else {
			prevcost = optcost;
			// TODO: integrate trajectory?
			// TODO: plot trajectory
		}
		
		//int num;
		//std::cin >> num;
	}
	return computeCost(B, U);

}
#endif


void pythonDisplayTrajectory(std::vector< Matrix<B_DIM> >& B, std::vector< Matrix<U_DIM> >& U)
{
	for (int t = 0; t < T-1; ++t) {
		B[t+1] = beliefDynamics(B[t], U[t]);
	}

	py::list Bvec;
	for(int j=0; j < B_DIM; j++) {
		for(int i=0; i < T; i++) {
			Bvec.append(B[i][j]);
		}
	}

	py::list Uvec;
		for(int j=0; j < U_DIM; j++) {
			for(int i=0; i < T-1; i++) {
			Uvec.append(U[i][j]);
		}
	}

	py::list x0_list, xGoal_list;
	for(int i=0; i < X_DIM; i++) {
		x0_list.append(x0[0,i]);
		xGoal_list.append(xGoal[0,i]);
	}

	std::string workingDir = boost::filesystem::current_path().normalize().string();
	std::string bspDir = workingDir.substr(0,workingDir.find("bsp"));

	try
	{
		Py_Initialize();
		py::object main_module = py::import("__main__");
		py::object main_namespace = main_module.attr("__dict__");
		py::exec("import sys, os", main_namespace);
		py::exec(py::str("sys.path.append('"+bspDir+"bsp/python')"), main_namespace);
		py::exec("from bsp_light_dark import LightDarkModel", main_namespace);
		py::object model = py::eval("LightDarkModel()", main_namespace);
		py::object plot_mod = py::import("plot");
		py::object plot_traj = plot_mod.attr("plot_belief_trajectory_cpp");

		plot_traj(Bvec, Uvec, model, x0_list, xGoal_list, T);
	}
	catch(py::error_already_set const &)
	{
		PyErr_Print();
	}

}

int main(int argc, char* argv[])
{
	x0[0] = -3.5; x0[1] = 2;
	SqrtSigma0 = identity<X_DIM>();
	xGoal[0] = -3.5; xGoal[1] = -2;

	xMin[0] = -5; xMin[1] = -3; 
	xMax[0] = 5; xMax[1] = 3;
	uMin[0] = -1; uMin[1] = -1;
	uMax[0] = 1; uMax[1] = 1;

	Matrix<U_DIM> uinit;
	uinit[0] = (xGoal[0] - x0[0]) / (T-1);
	uinit[1] = (xGoal[1] - x0[1]) / (T-1);
	
	std::vector<Matrix<U_DIM> > U(T-1, uinit); 

	std::vector<Matrix<B_DIM> > B(T);

	vec(x0, SqrtSigma0, B[0]);
	for (size_t t = 0; t < T-1; ++t) {
		B[t+1] = beliefDynamics(B[t], U[t]);
		//std::cout << ~B[t] << std::endl;
	}

	//for (size_t t = 0; t < T; ++t) {
	//	std::cout << ~B[t];
	//}

#ifdef BELIEF_MPC
	beliefMPC_params problem;
	beliefMPC_output output;
	beliefMPC_info info;
#endif

#ifdef BELIEF_PENALTY_MPC
	beliefPenaltyMPC_params problem;
	beliefPenaltyMPC_output output;
	beliefPenaltyMPC_info info;
#endif

	setupBeliefVars(problem, output);

	util::Timer solveTimer;
	util::Timer_tic(&solveTimer);
	
	// B&U optimized in-place
#ifdef BELIEF_MPC
	double cost = beliefCollocation(B, U, problem, output, info);
#endif

#ifdef BELIEF_PENALTY_MPC
	double cost = beliefPenaltyCollocation(B, U, problem, output, info);
#endif

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
