#include <vector>
#include <iomanip>


#include "../parameter.h"

#include "util/matrix.h"
#include "util/Timer.h"
#include "util/logging.h"
#include "util/utils.h"

#include <Python.h>
#include <boost/python.hpp>

#include <boost/timer.hpp>


namespace py = boost::python;

#define STATE_PENALTY_MPC
//#define BELIEF_MPC


extern "C" {
#include "parameter-state-casadi.h"
#include "statePenaltyMPC.h"
statePenaltyMPC_FLOAT **f, **lb, **ub, **C, **e, **z, **H;
}

Matrix<X_DIM> x0;
Matrix<X_DIM,X_DIM> SqrtSigma0;
Matrix<X_DIM> xGoal;
Matrix<X_DIM> xMin, xMax;
Matrix<U_DIM> uMin, uMax;

/*
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
*/
namespace cfg {
const double improve_ratio_threshold = .25;
const double min_approx_improve = 1e-2;
const double min_trust_box_size = 1e-3;
const double trust_shrink_ratio = .9;
const double trust_expand_ratio = 1.1;
const double cnt_tolerance = 1e-4;
const double penalty_coeff_increase_ratio = 5;
const double initial_penalty_coeff = 5;
const double initial_trust_box_size = 1;
const int max_penalty_coeff_increases = 2;
const int max_sqp_iterations = 50;
}


double computeCost(const std::vector< Matrix<B_DIM> >& B, const std::vector< Matrix<U_DIM> >& U)
{
	double cost = 0;
	Matrix<X_DIM> x;
	Matrix<X_DIM, X_DIM> SqrtSigma, Sigma;

	for(int t = 0; t < T-1; ++t) {
		unVec(B[t], x, SqrtSigma);
		Sigma = SqrtSigma*SqrtSigma;
		cost += alpha_joint_belief*tr(Sigma.subMatrix<J_DIM,J_DIM>(0,0)) +
				 alpha_param_belief*tr(Sigma.subMatrix<K_DIM,K_DIM>(J_DIM,J_DIM)) +
				 alpha_control*tr(~U[t]*U[t]);
	}
	unVec(B[T-1], x, SqrtSigma);
	Sigma = SqrtSigma*SqrtSigma;
	cost += alpha_final_joint_belief*tr(Sigma.subMatrix<J_DIM,J_DIM>(0,0)) +
			 alpha_final_param_belief*tr(Sigma.subMatrix<K_DIM,K_DIM>(J_DIM,J_DIM));
	return cost;
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
	params_arr[2] = alpha_final_joint_belief;

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

double costfuncTest(const std::vector< Matrix<X_DIM> >& X, const std::vector< Matrix<U_DIM> >& U,  Matrix<X_DIM, X_DIM> SqrtSigm0)
{
	double merit = 0;
	Matrix<X_DIM> x_t;
	Matrix<X_DIM, X_DIM> Sigma;
	Matrix<B_DIM> B_t;

	for(int t = 0; t < T-1; ++t) {
		vec(X[t], SqrtSigm0, B_t);

		Sigma=SqrtSigm0*SqrtSigm0;

		merit += alpha_belief*tr(Sigma)+
				 alpha_control*tr(~U[t]*U[t]);

		B_t = beliefDynamics(B_t,U[t]);
		unVec(B_t, x_t, SqrtSigm0);
	}

	merit += alpha_final_joint_belief*tr(Sigma);
	return merit;
}

Matrix<XU_DIM> finiteGradTest(const std::vector< Matrix<X_DIM> >& X, const std::vector< Matrix<U_DIM> >& U, const Matrix<X_DIM, X_DIM>& Sigma0){

	Matrix<XU_DIM> grad; 
	Matrix<X_DIM> x_t,xr,xl; 
	Matrix<U_DIM> u_t,ur,ul; 

	std::vector< Matrix<X_DIM> > XR = X;
	std::vector< Matrix<X_DIM> > XL = X; 

	std::vector< Matrix<U_DIM> > UR = U; 
	std::vector< Matrix<U_DIM> > UL = U; 
	int index=0; 

	for(int t=0; t<T-1; t++){
		x_t = X[t]; u_t = U[t]; 
		for(int i=0; i<X_DIM; i++){
			xr = x_t; xl = x_t;
			XL = X; XR = X; 

			xr[i] += diffEps; xl[i] -= diffEps;
			XR[t]=xr; 
			XL[t]=xl; 

			grad[index++] = (costfuncTest(XR,U,Sigma0) - costfuncTest(XL,U,Sigma0)) / (2.0*diffEps);
			
		}
		for(int i=0; i<U_DIM; i++){
			ur = u_t; ul = u_t;
			UL = U; UR=U; 
			ur[i] += diffEps; ul[i] -= diffEps;
			UR[t]=ur; 
			UL[t]=ul; 
			
			grad[index++] = (costfuncTest(X,UR,Sigma0) - costfuncTest(X,UL,Sigma0)) / (2.0*diffEps);
			
		}
	}
	return grad; 

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

}





// Jacobians: dg(b,u)/db, dg(b,u)/du
void linearizeBeliefDynamics(const Matrix<B_DIM>& b, const Matrix<U_DIM>& u, Matrix<B_DIM,B_DIM>& F, Matrix<B_DIM,U_DIM>& G, Matrix<B_DIM>& h)
{
	F.reset();
	Matrix<B_DIM> br(b), bl(b);
	for (size_t i = 0; i < B_DIM; ++i) {
		br[i] += diffEps; bl[i] -= diffEps;
		F.insert(0,i, (beliefDynamics(br, u) - beliefDynamics(bl, u)) / (br[i] - bl[i]));
		br[i] = b[i]; bl[i] = b[i];
	}

	G.reset();
	Matrix<U_DIM> ur(u), ul(u);
	for (size_t i = 0; i < U_DIM; ++i) {
		ur[i] += diffEps; ul[i] -= diffEps;
		G.insert(0,i, (beliefDynamics(b, ur) - beliefDynamics(b, ul)) / (ur[i] - ul[i]));
		ur[i] = u[i]; ul[i] = u[i];
	}

	h = beliefDynamics(b, u);
}


// utility to fill Matrix in column major format in FORCES array
template <size_t _numRows>
inline void fillCol(double *X, const Matrix<_numRows>& XCol) {
	int idx = 0;
	for(size_t r = 0; r < _numRows; ++r) {
		X[idx++] = XCol[r];
	}
}


void setupBeliefVars(statePenaltyMPC_params& problem, statePenaltyMPC_output& output)
{
	H = new statePenaltyMPC_FLOAT*[T];
	f = new statePenaltyMPC_FLOAT*[T];
	lb = new statePenaltyMPC_FLOAT*[T];
	ub = new statePenaltyMPC_FLOAT*[T];
	C = new statePenaltyMPC_FLOAT*[T-1];
	e = new statePenaltyMPC_FLOAT*[T];
	z = new statePenaltyMPC_FLOAT*[T];

#define SET_VARS(n)    \
		f[ BOOST_PP_SUB(n,1) ] = problem.f##n ;  \
		C[ BOOST_PP_SUB(n,1) ] = problem.C##n ;  \
		e[ BOOST_PP_SUB(n,1) ] = problem.e##n ;  \
		lb[ BOOST_PP_SUB(n,1) ] = problem.lb##n ;	\
		ub[ BOOST_PP_SUB(n,1) ] = problem.ub##n ;	\
		z[ BOOST_PP_SUB(n,1) ] = output.z##n ; \
		H[ BOOST_PP_SUB(n,1) ] = problem.H##n ;

#define BOOST_PP_LOCAL_MACRO(n) SET_VARS(n)
#define BOOST_PP_LOCAL_LIMITS (1, TIMESTEPS-1)
#include BOOST_PP_LOCAL_ITERATE()


	//f[ BOOST_PP_SUB(n,1) ] = problem.f##n ;
#define SET_LAST_VARS(n)    \
		lb[ BOOST_PP_SUB(n,1) ] = problem.lb##n ;	\
		ub[ BOOST_PP_SUB(n,1) ] = problem.ub##n ;	\
		z[ BOOST_PP_SUB(n,1) ] = output.z##n ; \
		e[ BOOST_PP_SUB(n,1) ] = problem.e##n ; \
		f[ BOOST_PP_SUB(n,1) ] = problem.f##n ; \
		H[ BOOST_PP_SUB(n,1) ] = problem.H##n ;

#define BOOST_PP_LOCAL_MACRO(n) SET_LAST_VARS(n)
#define BOOST_PP_LOCAL_LIMITS (TIMESTEPS, TIMESTEPS)
#include BOOST_PP_LOCAL_ITERATE()

for(int t = 0; t < T-1; ++t) {
	for(int i=0; i < (3*X_DIM+U_DIM); ++i) { H[t][i] = INFTY; }
	for(int i=0; i < (3*X_DIM+U_DIM); ++i) { f[t][i] = INFTY; }
	for(int i=0; i < (3*X_DIM+U_DIM); ++i) { lb[t][i] = INFTY; }
 	for(int i=0; i < (X_DIM+U_DIM); ++i) { ub[t][i] = INFTY; }
 	for(int i=0; i < (X_DIM*(3*X_DIM+U_DIM)); ++i) { C[t][i] = INFTY; }
 	for(int i=0; i < X_DIM; ++i) { e[t][i] = INFTY; }
 	for(int i=0; i < (X_DIM+U_DIM); ++i) { z[t][i] = INFTY; }
 }
 	for(int i=0; i < (X_DIM); ++i) { H[T-1][i] = INFTY; }
 	for(int i=0; i < (X_DIM); ++i) { f[T-1][i] = INFTY; }
 	for(int i=0; i < (X_DIM); ++i) { lb[T-1][i] = INFTY; }
 	for(int i=0; i < (X_DIM); ++i) { ub[T-1][i] = INFTY; }
 	for(int i=0; i < X_DIM; ++i) { e[T-1][i] = INFTY; }
 	for(int i=0; i < (X_DIM); ++i) { z[T-1][i] = INFTY; }
		
}

void cleanupBeliefMPCVars()
{
	delete[] f;
	delete[] lb;
	delete[] ub;
	delete[] C;
	delete[] e;
	delete[] z;
	delete[] H;
}


double computeMerit(const std::vector< Matrix<B_DIM> >& B, const std::vector< Matrix<U_DIM> >& U, double penalty_coeff)
{
	double merit = 0;
	Matrix<X_DIM> x;
	Matrix<X_DIM, X_DIM> SqrtSigma, Sigma;
	Matrix<B_DIM> dynviol;
	for(int t = 0; t < T-1; ++t) {
		unVec(B[t], x, SqrtSigma);
		Sigma = SqrtSigma*SqrtSigma;
		merit += alpha_belief*tr(Sigma)+
				 alpha_control*tr(~U[t]*U[t]);
		dynviol = (B[t+1] - beliefDynamics(B[t], U[t]) );
		for(int i = 0; i < B_DIM; ++i) {
			merit += penalty_coeff*fabs(dynviol[i]);
		}
	}
	unVec(B[T-1], x, SqrtSigma);
	Sigma = SqrtSigma*SqrtSigma;
	merit += alpha_final_joint_belief*tr(Sigma);
	return merit;
}

bool isValidInputs()
{
	//TODO:FIX INDEXING

	bool boundsCorrect = true;

	for(int t = 0; t < T-1; ++t) {
		std::cout << "t: " << t << "\n";
		for(int i=0; i < (3*X_DIM+U_DIM); ++i) { if (H[t][i] > INFTY/2) { std::cout << "H error: " << i << "\n"; } }
		for(int i=0; i < (3*X_DIM+U_DIM); ++i) { if (f[t][i] > INFTY/2) { std::cout << "f error: " << i << "\n"; } }
		for(int i=0; i < (3*X_DIM+U_DIM); ++i) { if (lb[t][i] > INFTY/2) { std::cout << "lb error: " << i << "\n"; } }
		for(int i=0; i < (X_DIM+U_DIM); ++i) {if (lb[t][i] > INFTY/2) { std::cout << "ub error: " << i << "\n"; } }
		for(int i=0; i < (X_DIM*(3*X_DIM+U_DIM)); ++i) { if (C[t][i] > INFTY/2) { std::cout << "C error: " << i << "\n"; } }
		for(int i=0; i < X_DIM; ++i) { if (e[t][i] > INFTY/2) { std::cout << "e error: " << i; } }
	}
	std::cout << "t: " << T-1 << "\n";
	for(int i=0; i < (X_DIM); ++i) { if (H[T-1][i] > INFTY/2) { std::cout << "H error: " << i << "\n"; } }
	for(int i=0; i < (X_DIM); ++i) { if (lb[T-1][i] > INFTY/2) { std::cout << "lb error: " << i << "\n"; } }
	for(int i=0; i < (X_DIM); ++i) { if (ub[T-1][i] > INFTY/2) { std::cout << "ub error: " << i << "\n"; } }
	for(int i=0; i < X_DIM; ++i) { if (e[T-1][i] > INFTY/2) { std::cout << "e error: " << i << "\n"; } }
	
	

	for(int t = 0; t < T-1; ++t) {

		std::cout << "t: " << t << std::endl << std::endl;


		std::cout << "f: ";
		for(int i = 0; i < (3*X_DIM+U_DIM); ++i) {
			std::cout << f[t][i] << " ";
		}
		std::cout << std::endl;



		std::cout << "lb x: ";
		for(int i = 0; i < X_DIM; ++i) {
			std::cout << lb[t][i] << " ";
		}
		std::cout << std::endl;


		std::cout << "lb u: ";
		for(int i = 0; i < U_DIM; ++i) {
			std::cout << lb[t][X_DIM+i] << " ";
		}
		std::cout << std::endl;

		/*
		std::cout << "lb s, t: ";
		for(int i = 0; i < 2*X_DIM; ++i) {
			std::cout << lb[t][X_DIM+U_DIM+i] << " ";
		}
		std::cout << std::endl;
		*/

		std::cout << "ub x: ";
		for(int i = 0; i < X_DIM; ++i) {
			std::cout << ub[t][i] << " ";
		}
		std::cout << std::endl;

		std::cout << "ub u: ";
		for(int i = 0; i < U_DIM; ++i) {
			std::cout << ub[t][X_DIM+i] << " ";
		}
		std::cout << std::endl;
		

		/*
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
		*/
		std::cout << "e:" << std::endl;
		for(int i = 0; i < X_DIM; ++i) {
			std::cout << e[t][i] << " ";
		}
		std::cout << std::endl;
		
		std::cout << std::endl << std::endl;

		for(int i = 0; i < X_DIM; ++i) {
			boundsCorrect &= (lb[t][i] < ub[t][i]);
		}
	}

	for(int i = 0; i < X_DIM; ++i) {
		boundsCorrect &= (lb[T-1][i] < ub[T-1][i]);
	}

	std::cout << "boundsCorrect: " << boundsCorrect << std::endl;

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


// Jacobians: dg(b,u)/db, dg(b,u)/du
void linearizeArmDynamics(const Matrix<X_DIM>& x, const Matrix<U_DIM>& u, Matrix<X_DIM,X_DIM>& F, Matrix<X_DIM,U_DIM>& G, Matrix<X_DIM>& h)
{
	F.reset();
	Matrix<X_DIM> xr(x), xl(x);
	for (size_t i = 0; i < X_DIM; ++i) {
		xr[i] += diffEps; xl[i] -= diffEps;
		F.insert(0,i, (dynfunc(xr, u, zeros<Q_DIM,1>()) - dynfunc(xl, u, zeros<Q_DIM,1>())) / (xr[i] - xl[i]));
		xr[i] = x[i]; xl[i] = x[i];
	}

	G.reset();
	Matrix<U_DIM> ur(u), ul(u);
	for (size_t i = 0; i < U_DIM; ++i) {
		ur[i] += diffEps; ul[i] -= diffEps;
		G.insert(0,i, (dynfunc(x, ur, zeros<Q_DIM,1>()) - dynfunc(x, ul, zeros<Q_DIM,1>())) / (ur[i] - ul[i]));
		ur[i] = u[i]; ul[i] = u[i];

	}

	h = dynfunc(x, u, zeros<Q_DIM,1>());
}

bool minimizeMeritFunction(std::vector< Matrix<X_DIM> >& X, std::vector< Matrix<U_DIM> >& U, statePenaltyMPC_params& problem, statePenaltyMPC_output& output, statePenaltyMPC_info& info, double penalty_coeff)
{
	LOG_DEBUG("Solving sqp problem with penalty parameter: %2.4f", penalty_coeff);

	std::vector< Matrix<X_DIM,X_DIM> > F(T-1);
	std::vector< Matrix<X_DIM,U_DIM> > G(T-1);
	std::vector< Matrix<X_DIM> > h(T-1);

	double Xeps = cfg::initial_trust_box_size;
	double Ueps = cfg::initial_trust_box_size;



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

		//Grad = finiteGradTest(X, U,SqrtSigma0);

		/*for (int t=0; t<T-1; t++){
			std::cout<<"t "<<t<<"\n";
			for (int j=0; j<X_DIM; j++){
				std::cout<<Grad[index++]<<" ";
			}
			std::cout<<"\n";
			for (int j=0; j<U_DIM; j++){
				std::cout<<Grad[index++]<<" ";
			}
			std::cout<<"\n";
		}
		*/
		//std::cout<<"Gradcoo "<<Grad<<"\n";
		// Problem linearization and definition
		// fill in H, f

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
				HMat(i,i) = B(idx+i,idx+i);
			}

			// since diagonal, fill directly
			// TODO: check if H is 0 for the landmarks. if so, pull out of the state
			for(int i = 0; i < (X_DIM+U_DIM); ++i) { H[t][i] = HMat(i,i); }
			// TODO: why does this work???
			for(int i = 0; i < (2*X_DIM); ++i) { H[t][i + (X_DIM+U_DIM)] = 0; } //1e4

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
			linearizeArmDynamics(xt, ut, F[t], G[t], h[t]);


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

		// For last stage, fill in H, f
		Matrix<X_DIM>& xT = X[T-1];

		idx = (T-1)*(X_DIM+U_DIM);
		//LOG_DEBUG("idx: %d",idx);

		HfMat.reset();
		for(int i = 0; i < X_DIM; ++i) {
			double val = B(idx+i,idx+i);
			HfMat(i,i) = (val < 0) ? 0 : val;
		}

		// since diagonal, fill directly
		for(int i = 0; i < X_DIM; ++i) { H[T-1][i] = HfMat(i,i); }

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

			// Fill in lb, ub
			index = 0;
			for(int i = 0; i < X_DIM; ++i) { lb[T-1][index++] = MAX(xMin[i], xT[i] - Xeps); }
	
			// sigma lower bound

			index = 0;
			// xGoal upper bound
			for(int i = 0; i < X_DIM; ++i) { ub[T-1][index++] = MIN(xMax[i], xT[i] + Xeps); }

			// Verify problem inputs
			/*if (!isValidInputs()) {
				std::cout << "Inputs are not valid!" << std::endl;
				exit(-1);
			}
			*/



			int exitflag = statePenaltyMPC_solve(&problem, &output, &info);
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
				std::cout << "penalty coeff: " << penalty_coeff << "\n";
				std::exit(-1);
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
			//std::cin.ignore();
			//int num;
			//std::cin >> num;
			//exit(0);
			if (approx_merit_improve < -1e-5) {
				LOG_ERROR("Approximate merit function got worse: %1.6f", approx_merit_improve);
				LOG_ERROR("Either convexification is wrong to zeroth order, or you are in numerical trouble");
				LOG_ERROR("Failure!");
				

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


				// take in diagonal of B
				// find minimum value among diagonal elements
				// negate and add to other vals
				double minValue = INFTY;
				double maxValue = -INFTY;
				for(int i=0; i < XU_DIM; ++i) {
					minValue = MIN(minValue, B(i,i));
					maxValue = MAX(maxValue, B(i,i));
				}

				//std::cout << "minValue: " << minValue << "\n";
				//std::cout << "maxValue: " << maxValue << "\n\n";
				if (minValue < 0) {
					//std::cout << "negative minValue, press enter\n";
					//std::cin.ignore();
					B = B + fabs(minValue)*identity<XU_DIM>();
				}
				

				

				// Do not update B
				//B = identity<XU_DIM>();

				X = Xopt; U = Uopt;

				LOG_DEBUG("Accepted, Increasing trust region size to:  %2.6f %2.6f", Xeps,Ueps);
				break;
			}

			if (Xeps < cfg::min_trust_box_size && Ueps < cfg::min_trust_box_size) {
		 
			    LOG_DEBUG("Converged: x tolerance");
			    return true;
			}

			//std::cout << "U" << std::endl;
			//for(int t=0; t < T-1; ++t) {
			//	std::cout << ~U[t];
			//}
			//std::cout << std::endl << std::endl;
			//pythonDisplayTrajectory(U, T, false);
			//pythonDisplayTrajectory(X, T, true);

		} // trust region loop
		sqp_iter++;
	} // sqp loop

	return success;
}

double statePenaltyCollocation(std::vector< Matrix<X_DIM> >& X, std::vector< Matrix<U_DIM> >& U, statePenaltyMPC_params& problem, statePenaltyMPC_output& output, statePenaltyMPC_info& info)
{
	double penalty_coeff = cfg::initial_penalty_coeff;
	//double trust_box_size = cfg::initial_trust_box_size;

	int penalty_increases = 0;

	Matrix<X_DIM> dynviol;

	// penalty loop
	while(penalty_increases < cfg::max_penalty_coeff_increases)
	{
		bool success = minimizeMeritFunction(X, U, problem, output, info, penalty_coeff);

		double cntviol = 0;
		for(int t = 0; t < T-1; ++t) {
			dynviol = (X[t+1] - dynfunc(X[t], U[t], zeros<Q_DIM,1>()));
			for(int i = 0; i < X_DIM; ++i) {
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
	    	//return computeCost(X, U);
	    	return casadiComputeCost(X, U);
	    }
	}
	//return computeCost(X, U);
	return casadiComputeCost(X, U);
}



int main(int argc, char* argv[])
{

	// actual: 6.66, 6.66, 10.86, 13
	double length1_est = .3,
				length2_est = .7,
				mass1_est = .35,
				mass2_est = .35;

	// position, then velocity
	x0[0] = M_PI*0.5; x0[1] = M_PI*0.5; x0[2] = 0; x0[3] = 0;
	// parameter start estimates (alphabetical, then numerical order)
	x0[4] = 1/length1_est; x0[5] = 1/length2_est; x0[6] = 1/mass1_est; x0[7] = 1/mass2_est;


	Matrix<X_DIM> x_real;
	x_real[0] = M_PI*0.45; x_real[1] = M_PI*0.55; x_real[2] = -0.01; x_real[3] = 0.01;
	x_real[4] = 1/dynamics::length1; x_real[5] = 1/dynamics::length2; x_real[6] = 1/dynamics::mass1; x_real[7] = 1/dynamics::mass2;


	xGoal[0] = M_PI*0.5; xGoal[1] = M_PI*0.5; xGoal[2] = 0.0; xGoal[3] = 0.0;
	xGoal[4] = 1/length1_est; xGoal[5] = 1/length2_est; xGoal[6] = 1/mass1_est; xGoal[7] = 1/mass2_est;

	// from original file, possibly change
	SqrtSigma0(0,0) = 0.1;
	SqrtSigma0(1,1) = 0.1;
	SqrtSigma0(2,2) = 0.05;
	SqrtSigma0(3,3) = 0.05;
	SqrtSigma0(4,4) = 0.5;
	SqrtSigma0(5,5) = 0.5;
	SqrtSigma0(6,6) = 0.5;
	SqrtSigma0(7,7) = 0.5;


	xMin[0] = -1000; // joint pos 1
	xMin[1] = -1000; // joint pos 2
	xMin[2] = -1000; // joint vel 1
	xMin[3] = -1000; // joint vel 2
	xMin[4] = 0.01; // 1/length1
	xMin[5] = 0.01; // 1/length2
	xMin[6] = 0.01; // 1/mass1
	xMin[7] = 0.01; // 1/mass2

	xMax[0] = 1000; // joint pos 1
	xMax[1] = 1000; // joint pos 2
	xMax[2] = 1000; // joint vel 1
	xMax[3] = 1000; // joint vel 2
	xMax[4] = 100; // 1/length1
	xMax[5] = 100; // 1/length2
	xMax[6] = 100; // 1/mass1
	xMax[7] = 100; // 1/mass2

	for(int i = 0; i < U_DIM; ++i) {
		uMin[i] = -0.1;
		uMax[i] = 0.1;
	}

	//Matrix<U_DIM> uinit = (xGoal.subMatrix<U_DIM,1>(0,0) - x0.subMatrix<U_DIM,1>(0,0))/(double)(T-1);
	Matrix<U_DIM> uinit;
	uinit[0] = 0.0;
	uinit[1] = 0.0;
	
	std::vector<Matrix<U_DIM> > U(T-1, uinit); 
	std::vector<Matrix<X_DIM> > X(T);
	std::vector<Matrix<B_DIM> > B(T);

	std::vector<Matrix<U_DIM> > HistoryU(HORIZON);
	std::vector<Matrix<B_DIM> > HistoryB(HORIZON); 

	// pythonPlotRobot(U, SqrtSigma0, x0, xGoal);

	statePenaltyMPC_params problem;
	statePenaltyMPC_output output;
	statePenaltyMPC_info info;

	setupBeliefVars(problem, output);

	vec(x0, SqrtSigma0, B[0]);
	std::cout<<"HORIZON is "<<HORIZON<<'\n';

	for(int h = 0; h < HORIZON; ++h) {
		for (int t = 0; t < T-1; ++t) {
			X[t] = B[t].subMatrix<X_DIM,1>(0,0);
			B[t+1] = beliefDynamics(B[t], U[t]);
		}
		X[T-1] = B[T-1].subMatrix<X_DIM,1>(0,0);
		util::Timer solveTimer;
		util::Timer_tic(&solveTimer);
		
		
		double cost = statePenaltyCollocation(X, U, problem, output, info);
		
	

		//pythonDisplayTrajectory(U, SqrtSigma0, x0, xGoal);
		//pythonPlotRobot(U, SqrtSigma0, x0, xGoal);

		double solvetime = util::Timer_toc(&solveTimer);
		std::cout<<"Solve time: "<<solvetime*1000<<"\n";
		//LOG_INFO("Optimized cost: %4.10f", cost);
		//LOG_INFO("Solve time: %5.3f ms", solvetime*1000);
		
		unVec(B[0], x0, SqrtSigma0);

		//std::cout << "x0 before control step" << std::endl << ~x0;
		//std::cout << "control u: " << std::endl << ~U[0];

		
	
		HistoryU[h] = U[0];
		HistoryB[h] = B[0];


	
	
		B[0] = executeControlStep(x_real, B[0], U[0]);
		


		unVec(B[0], x0, SqrtSigma0);
		//std::cout << "x0 after control step" << std::endl << ~x0;
		#define SPEED_TEST
		#ifdef SPEED_TEST
		for(int t = 0; t < T-1; ++t) {
			for(int l=0; l<U_DIM; l++){
		
				U[t][l] = 0;
			}
		}
		#else
		for(int t = 0; t < T-2; ++t) {
		
			U[t] = U[t+1];
		}
		#endif


	}


	pythonDisplayHistory(HistoryU,HistoryB, SqrtSigma0, x0, HORIZON);
	cleanupBeliefMPCVars();




#define PLOT
#ifdef PLOT
	pythonDisplayHistory(HistoryU,HistoryB, SqrtSigma0, x0, HORIZON);

	//int k;
	//std::cin >> k;
#endif

	return 0;
}
