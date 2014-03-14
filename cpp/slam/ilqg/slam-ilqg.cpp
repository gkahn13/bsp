#include <vector>
#include <iomanip>

#include "util/Timer.h"

#include "../slam.h"
#include "../traj/slam-traj.h"
#include "ilqg.h"

const double alpha_belief = 10, alpha_final_belief = 10, alpha_control = 1, alpha_goal_state = 1;

// iLQG params
SymmetricMatrix<U_DIM> Rint;
SymmetricMatrix<X_DIM> Qint;
SymmetricMatrix<X_DIM> QGoal;
SymmetricMatrix<X_DIM> Sigma0;

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

Matrix<X_DIM> f(const Matrix<X_DIM>& x, const Matrix<U_DIM>& u)
{
	return dynfunc(x, u, zeros<Q_DIM,1>());
}

// Observation model
Matrix<Z_DIM> h(const Matrix<X_DIM>& x)
{
	return obsfunc(x, zeros<R_DIM,1>());
}

/*
double computeLQGMPcost(const std::vector<Matrix<X_DIM> >& X, const std::vector<Matrix<U_DIM> >& U)
{
	std::vector<Matrix<G_DIM, X_DIM> > J;
	std::vector<Matrix<X_DIM, X_DIM> > Sigma;
	preprocess(X, J, Sigma); // TODO: need to do same preprocess?

	double cost = 0;
	for(int t = 0; t < T-1; ++t) {
		cost += alpha_belief*tr(J[t]*Sigma[t]*~J[t]) + alpha_control*tr(~U[t]*U[t]);
	}
	cost += alpha_final_belief*tr(J[T-1]*Sigma[T-1]*~J[T-1]);
	return cost;
}
*/

//double computeCost(const std::vector< Matrix<B_DIM> >& B, const std::vector< Matrix<U_DIM> >& U)
//{
//	double cost = 0;
//	Matrix<X_DIM> x;
//	Matrix<X_DIM, X_DIM> SqrtSigma;
//
//	for(int t = 0; t < T-1; ++t) {
//		unVec(B[t], x, SqrtSigma);
//		cost += alpha_belief*tr(SqrtSigma*SqrtSigma) + alpha_control*tr(~U[t]*U[t]);
//	}
//	unVec(B[T-1], x, SqrtSigma);
//	cost += alpha_final_belief*tr(SqrtSigma*SqrtSigma);
//	return cost;
//}

// Jacobian df/dx(x,u)
template <size_t _xDim, size_t _uDim>
inline Matrix<_xDim,_xDim> dfdx(Matrix<_xDim> (*f)(const Matrix<_xDim>&, const Matrix<_uDim>&), const Matrix<_xDim>& x, const Matrix<_uDim>& u) 
{
	Matrix<_xDim,_xDim> A;
	Matrix<_xDim> xr(x), xl(x);
	for (size_t i = 0; i < _xDim; ++i) {
		xr[i] += step; xl[i] -= step;
		A.insert(0,i, (f(xr, u) - f(xl, u)) / (2.0*step));
		xr[i] = xl[i] = x[i];
	}
	return A;
}

// Jacobian df/du(x,u)
template <size_t _xDim, size_t _uDim>
inline Matrix<_xDim,_uDim> dfdu(Matrix<_xDim> (*f)(const Matrix<_xDim>&, const Matrix<_uDim>&),	const Matrix<_xDim>& x, const Matrix<_uDim>& u) 
{
	Matrix<_xDim,_uDim> B;
	Matrix<_uDim> ur(u), ul(u);
	for (size_t i = 0; i < _uDim; ++i) {
		ur[i] += step; ul[i] -= step;
		B.insert(0,i, (f(x, ur) - f(x, ul)) / (2.0*step));
		ur[i] = ul[i] = u[i];
	}
	return B;
}

// Jacobian dh/dx(x)
template <size_t _xDim, size_t _zDim>
inline Matrix<_zDim,_xDim> dhdx(Matrix<_zDim> (*h)(const Matrix<_xDim>&), const Matrix<_xDim>& x) 
{
	Matrix<_zDim,_xDim> H;
	Matrix<_xDim> xr(x), xl(x);
	for (size_t i = 0; i < _xDim; ++i) {
		xr[i] += step; xl[i] -= step;
		H.insert(0,i, (h(xr) - h(xl)) / (2.0*step));
		xr[i] = xl[i] = x[i];
	}
	return H;
}

inline SymmetricMatrix<X_DIM> varM(const Matrix<X_DIM>& x, const Matrix<U_DIM>& u)
{
	Matrix<C_DIM,C_DIM> Acar;
	Matrix<C_DIM,Q_DIM> Mcar;
	linearizeDynamics(x, u, zeros<Q_DIM,1>(), Acar, Mcar);

	Matrix<X_DIM,Q_DIM> M = zeros<X_DIM,Q_DIM>();
	M.insert<C_DIM, 2>(0, 0, Mcar);

	SymmetricMatrix<X_DIM> S = SymProd(M*Q,~M);

	return S;
}

inline SymmetricMatrix<Z_DIM> varN(const Matrix<X_DIM>& x)
{
	return R;
}

// TODO: needs to be changed?
// Compute closed form quadratic finalCost function around around b
inline void quadratizeFinalCost(const Matrix<X_DIM>& xBar, const SymmetricMatrix<X_DIM>& SigmaBar, double& s, SymmetricMatrix<X_DIM>& S, Matrix<1,X_DIM>& sT, Matrix<1,S_DIM>& tT, unsigned int flag) 
{
	if (flag & COMPUTE_S) S = QGoal;
	if (flag & COMPUTE_sT) sT = ~(xBar - xGoal)*QGoal;
	if (flag & COMPUTE_s) s = 0.5*scalar(~(xBar - xGoal)*QGoal*(xBar - xGoal)) + scalar(vecTh(QGoal)*vec(SigmaBar));
	//if (flag & COMPUTE_s) s = 2.0*(0.5*scalar(~(xBar - xGoal)*QGoal*(xBar - xGoal)) + scalar(vecTh(QGoal)*vec(SigmaBar)));
	if (flag & COMPUTE_tT) tT = vecTh(QGoal);
}

// TODO: needs to be changed?
// Compute closed form quadratic cost function around around b and u
inline bool quadratizeCost(const Matrix<X_DIM>& xBar, const SymmetricMatrix<X_DIM>& SigmaBar, const Matrix<U_DIM>& uBar, double& q, SymmetricMatrix<X_DIM>& Q, SymmetricMatrix<U_DIM>& R, 
						   Matrix<U_DIM, X_DIM>& P, Matrix<1,X_DIM>& qT, Matrix<1,U_DIM>& rT, Matrix<1,S_DIM>& pT, unsigned int flag) 
{
	if (flag & COMPUTE_Q) Q = zeros<X_DIM>(); // penalize uncertainty
	if (flag & COMPUTE_R) R = Rint;
	if (flag & COMPUTE_P) P = zeros<U_DIM,X_DIM>();
	if (flag & COMPUTE_qT) qT = zeros<1, X_DIM>();
	if (flag & COMPUTE_rT) rT = ~uBar*Rint;
	if (flag & COMPUTE_pT) pT = vecTh(Qint);
	if (flag & COMPUTE_q) q = 0.5*scalar(~uBar*Rint*uBar) + scalar(vecTh(Qint)*vec(SigmaBar));
	//if (flag & COMPUTE_q) q = 2*(0.5*scalar(~uBar*Rint*uBar) + scalar(vecTh(Qint)*vec(SigmaBar)));

	return true;
}

inline void linearizeDynamics(const Matrix<X_DIM>& xBar, const Matrix<U_DIM>& uBar, Matrix<X_DIM>& c, Matrix<X_DIM, X_DIM>& A, Matrix<X_DIM, U_DIM>& B, SymmetricMatrix<X_DIM>& M, unsigned int flag)
{
	if (flag & COMPUTE_c) c = f(xBar, uBar);
	if (flag & COMPUTE_A) A = dfdx(f, xBar, uBar);
	if (flag & COMPUTE_B) B = dfdu(f, xBar, uBar);
	if (flag & COMPUTE_M) M = varM(xBar, uBar);
}

inline void linearizeObservation(const Matrix<X_DIM>& xBar, Matrix<Z_DIM, X_DIM>& H, SymmetricMatrix<Z_DIM>& N)
{
	H = dhdx(h, xBar);
	N = varN(xBar);
}


void planPath(std::vector<Matrix<P_DIM> > l, std::ofstream& f) {
	LOG_INFO("Initializing problem parameters");
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
			LOG_ERROR("Failed to initialize trajectory, continuing anyways");
			//exit(-1);
		}

		double initTrajTime = util::Timer_toc(&trajTimer);
		trajTime += initTrajTime;

		vec(x0, SqrtSigma0, B[0]);
		for(int t=0; t < T-1; ++t) {
			B[t+1] = beliefDynamics(B[t], U[t]);
		}


		//pythonDisplayTrajectory(B, U, waypoints, landmarks, T, true);

		std::vector< Matrix<U_DIM, X_DIM> > L;
		std::vector< Matrix<X_DIM> > xBar;
		std::vector< SymmetricMatrix<X_DIM> > SigmaBar;

		xBar.push_back(x0);

		Sigma0 = SymProd(SqrtSigma0,SqrtSigma0);
		SigmaBar.push_back(Sigma0);

		util::Timer_tic(&solveTimer);

		solvePOMDP(linearizeDynamics, linearizeObservation, quadratizeFinalCost, quadratizeCost, xBar, SigmaBar, U, L);

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

		//totalTrajCost += computeCost(X,U);

		unVec(B[T-1], x0, SqrtSigma0);

		//pythonDisplayTrajectory(B, U, waypoints, landmarks, T, true);

	}

	//LOG_INFO("Total trajectory cost: %4.10f", totalTrajCost);
	LOG_INFO("Total trajectory solve time: %5.3f ms", trajTime*1000);
	LOG_INFO("Total solve time: %5.3f ms", totalSolveTime*1000);

	logDataToFile(f, B_total, totalSolveTime*1000, trajTime*1000, 0);


	pythonDisplayTrajectory(B_total, U_total, waypoints, landmarks, T*NUM_WAYPOINTS, false);
}

int main(int argc, char* argv[])
{
	Rint = alpha_control*identity<U_DIM>();
	Qint = alpha_belief*identity<X_DIM>();
	QGoal = alpha_final_belief*identity<X_DIM>();

	std::vector<std::vector<Matrix<P_DIM> > > l_list = landmarks_list();

	std::ofstream f;
	logDataHandle("slam/data/slam-ilqg", f);

	for(size_t i=0; i < l_list.size(); ++i) {
		planPath(l_list[i], f);
	}

	return 0;
}


