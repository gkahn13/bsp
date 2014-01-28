#include <vector>
#include <iomanip>

#include "callisto.h"
#include "matrix.h"
#include "utils.h"

#include "iLQG.h"

#include "arm.h"

Matrix<X_DIM> f(const Matrix<X_DIM>& x, const Matrix<U_DIM>& u)
{
	Matrix<X_DIM> xNew = x + u*DT;
	return xNew;
}

// Observation model
Matrix<Z_DIM> h(const Matrix<X_DIM>& x)
{
	Matrix<G_DIM> ee_pos = g(x);

	Matrix<Z_DIM> obs;
	obs[0] = (ee_pos[0] - cam0[0]) / (ee_pos[1] - cam0[1]);
	obs[1] = (ee_pos[2] - cam0[2]) / (ee_pos[1] - cam0[1]);
	obs[2] = (ee_pos[0] - cam1[0]) / (ee_pos[1] - cam1[1]);
    obs[3] = (ee_pos[2] - cam1[2]) / (ee_pos[1] - cam1[1]);

    return obs;
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
	SymmetricMatrix<U_DIM> S = identity<U_DIM>();
	for(int i = 0; i < U_DIM; ++i) {
		S(i,i) = 0.01;
	}
	return S;
}

inline SymmetricMatrix<Z_DIM> varN(const Matrix<X_DIM>& x)
{
	SymmetricMatrix<Z_DIM> S = identity<Z_DIM>();
	for(int i = 0; i < Z_DIM; ++i) {
		S(i,i) = 0.01;
	}
	return S;
}

// Compute closed form quadratic finalCost function around around b
inline void quadratizeFinalCost(const Matrix<X_DIM>& xBar, const SymmetricMatrix<X_DIM>& SigmaBar, double& s, SymmetricMatrix<X_DIM>& S, Matrix<1,X_DIM>& sT, Matrix<1,S_DIM>& tT, unsigned int flag) 
{
	if (flag & COMPUTE_S) S = QGoal;
	if (flag & COMPUTE_sT) sT = ~(xBar - xGoal)*QGoal;
	if (flag & COMPUTE_s) s = 0.5*scalar(~(xBar - xGoal)*QGoal*(xBar - xGoal)) + scalar(vecTh(QGoal)*vec(SigmaBar));
	//if (flag & COMPUTE_s) s = 2.0*(0.5*scalar(~(xBar - xGoal)*QGoal*(xBar - xGoal)) + scalar(vecTh(QGoal)*vec(SigmaBar)));
	if (flag & COMPUTE_tT) tT = vecTh(QGoal);
}

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
	if (flag & COMPUTE_A) A = dfdx(f, xBar, uBar); //identity<X_DIM>();
	if (flag & COMPUTE_B) B = dfdu(f, xBar, uBar); //identity<U_DIM>();
	if (flag & COMPUTE_M) M = varM(xBar, uBar);
}

inline void linearizeObservation(const Matrix<X_DIM>& xBar, Matrix<Z_DIM, X_DIM>& H, SymmetricMatrix<Z_DIM>& N)
{
	H = dhdx(h, xBar); // identity<X_DIM>();
	N = varN(xBar);
}

bool testInitializationFeasibility(std::vector<Matrix<B_DIM> >& B, std::vector<Matrix<U_DIM> >& U)
{
	//std::cout << "X initial" << std::endl;
	for (int t = 0; t < T; ++t) { 
		Matrix<U_DIM>& xt = B[t].subMatrix<X_DIM,1>(0,0);
		for (int i = 0; i < X_DIM; ++i) {
			if (xt[i] > xMax[i] || xt[i] < xMin[i]) {
				std::cerr << "Joint angle limit violated at joint " << i << " and time " << t << std::endl;
				return false;
			}
		}

		//std::cout << std::setprecision(8) << ~xt; 
	}

	//std::cout << "U initial" << std::endl;
	for (int t = 0; t < T-1; ++t) { 
		Matrix<U_DIM>& ut = U[t];
		for (int i = 0; i < U_DIM; ++i) {
			if (ut[i] > uMax[i] || ut[i] < uMin[i]) {
				std::cerr << "Control limit violated at joint " << i << " and time " << t << std::endl;
				return false;
			}
		}
		//std::cout << std::setprecision(8) << ~ut; 
	}

	return true;
}

void displayILQGTrajectory(std::vector<Matrix<U_DIM> >& U) 
{
	int nsteps = (int)U.size()+1;
	
	std::vector<Matrix<B_DIM> > B(T);
	vec(x0, SqrtSigma0, B[0]);
	for (size_t t = 0; t < T-1; ++t) {
		B[t+1] = beliefDynamics(B[t], U[t]);
		//std::cout << ~B[t] << std::endl;
	}

	Matrix<X_DIM, X_DIM> SqrtSigma;
	Matrix<X_DIM> x;
	
	int np[1] = {nsteps};
	float* p = new float[3*np[0]];
	int idx = 0;
	int npts = 0;
	Matrix<G_DIM> q = g(x0);
	p[idx] = (float) q[0]; p[idx+1] = (float)q[1]; p[idx+2] = (float)q[2];
	idx += 3;
	for (int t = 0; t < nsteps-1; ++t) {
		unVec(B[t], x, SqrtSigma);
		q = g(x);
		p[idx] = (float) q[0]; p[idx+1] = (float)q[1]; p[idx+2] = (float)q[2];
		idx += 3;
	}
	//CAL_CreatePolyline(cal_paths, 1, np, p, 0);
	drawPolyline3d(np[0], p, cal_paths);
	delete[] p;

	std::cout << "Last goal: " << ~g(x) << std::endl;

	//visualize last state
	//showState(x);
	
	Matrix<G_DIM, X_DIM> J;
	
	for (int t = 0; t < nsteps; ++t) 
	{
		Matrix<G_DIM,G_DIM> EVec, EVal;

		unVec(B[t], x, SqrtSigma);
		
		linearizeg(x, J);
		jacobi(J*SqrtSigma*SqrtSigma*~J, EVec, EVal);

		Matrix<4,1> q = quatFromRot(EVec);
		Matrix<G_DIM> p = g(B[t].subMatrix<X_DIM,1>(0,0));

		float pos[3] = {(float) p[0], (float) p[1], (float) p[2]};
		float quat[4] = {(float) q(0,0), (float) q(1,0), (float) q(2,0), (float) q(3,0)};
		float scale[3] = {(float) (2*sqrt(EVal(0,0))), (float) (2*sqrt(EVal(1,1))), (float) (2*sqrt(EVal(2,2)))};
		
		CAL_AddGroupKeyState(cal_ellipse, (float) (t * DT), pos, quat, scale, true);
		
		// Add joint angles for display along trajectory
		for (int j = 0; j < 6; ++j) {
			double angle = x[j] + null_joint[j];
			if (angle < 0) {
				angle += 2*M_PI;
			} else if (angle >= 2*M_PI) {
				angle -= 2*M_PI;
			}

			float rot[3] = {0.f, (float)angle, 0.f}; 
			CAL_AddGroupKeyState(joint_group[j], (float) (t * DT), 0, rot, 0, false);
		}
	}
}

int main(int argc, char* argv[])
{
	initProblemParams();

	Rint = alpha_control*identity<U_DIM>();
	Qint = alpha_belief*identity<X_DIM>();
	QGoal = alpha_final_belief*identity<X_DIM>();

	std::vector< Matrix<U_DIM, X_DIM> > L;
	std::vector< Matrix<X_DIM> > xBar;
	std::vector< SymmetricMatrix<X_DIM> > SigmaBar;
	
	Matrix<U_DIM> uinit = (xGoal - x0) / (double)((T-1)*DT);
	std::vector<Matrix<U_DIM> > uBar(T-1, uinit);

	Sigma0 = identity<X_DIM>();
	xBar.push_back(x0);
	SigmaBar.push_back(Sigma0);

	// Create Obstacles in Callisto
	CAL_Initialisation (true,true,true);
	initEnvironment();

	//displayILQGTrajectory(uBar);

	Timer solveTimer;
	//Timer_tic(&solveTimer);
	solveTimer.start();

	// solve ilqg here
	solvePOMDP(linearizeDynamics, linearizeObservation, quadratizeFinalCost, quadratizeCost, xBar, SigmaBar, uBar, L);

	solveTimer.stop();
	//double solvetime = util::Timer_toc(&solveTimer);

	//LOG_INFO("Optimized cost: %4.10f", cost);
	//LOG_INFO("Solve time: %5.3f ms", solvetime*1000);
	
	//std::cout << "Actual cost: " << computeCost(B, U) << std::endl;
	std::cout << "Optimization time: " << solveTimer.interval_mS() << " mS" << std::endl;
	
	//saveOptimizedTrajectory(U);
	//readOptimizedTrajectory(U);
	
	displayILQGTrajectory(uBar);
	/*
	vec(x0, SqrtSigma0, B[0]);
	for (size_t t = 0; t < T-1; ++t) {
		B[t+1] = beliefDynamics(B[t], U[t]);
	}
	*/

	//displayTrajectory(X, U, B);

	int k;
	std::cin >> k;

	CAL_End();

	return 0;
}
