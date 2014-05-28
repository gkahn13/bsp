#include "casadi-slam-control.h"
#include "casadi-slam.h"


// params[0] = alpha_belief
// params[1] = alpha_control
// params[2] = alpha_final_belief
// params[3] = alpha_goal_state
CasADi::SXMatrix costfunc(const CasADi::SXMatrix& U, const CasADi::SXMatrix& x0, const CasADi::SXMatrix& Sigma_0, const CasADi::SXMatrix& xGoal, const CasADi::SXMatrix& params)
{
	CasADi::SXMatrix cost = 0;

	CasADi::SXMatrix x_tp1(X_DIM,1);
	CasADi::SXMatrix Sigma_t = Sigma_0, Sigma_tp1(X_DIM,X_DIM);
	CasADi::SXMatrix x_t = x0, u_t(U_DIM,1);

	int offset = 0;

	for (int t = 0; t < (T-1); ++t)
	{
		u_t = U(CasADi::Slice(offset,offset+U_DIM));
		offset += U_DIM;

		cost += params[0]*trace(Sigma_t);
		cost += params[1]*inner_prod(u_t, u_t);

		EKF(x_t, u_t, Sigma_t, x_tp1, Sigma_tp1);
		Sigma_t = Sigma_tp1;
		x_t = x_tp1;
	}

	cost += params[2]*trace(Sigma_t) + params[3]*inner_prod(x_t - xGoal, x_t - xGoal);

	return cost;
}

// params[0] = alpha_belief NOT USED
// params[1] = alpha_control
// params[2] = alpha_final_belief
// params[3] = alpha_goal_state
CasADi::SXMatrix costfuncinfo(const CasADi::SXMatrix& U, const CasADi::SXMatrix& x0, const CasADi::SXMatrix& Sigma_0, const CasADi::SXMatrix& xGoal, const CasADi::SXMatrix& params)
{
	CasADi::SXMatrix cost = 0;

	std::vector<CasADi::SXMatrix> B_ham(T), C_ham(T);
	B_ham[0] = Sigma_0;
	C_ham[0] = CasADi::SXMatrix(CasADi::DMatrix::eye(X_DIM));

	std::vector<CasADi::SXMatrix> X(T);
	X[0] = x0;

	int offset = 0;
	for(int t=0; t < T-1; ++t) {
		const CasADi::SXMatrix u_t = U(CasADi::Slice(offset,offset+U_DIM));
		offset += U_DIM;

		CasADi::SXMatrix A(X_DIM,X_DIM), M(X_DIM,Q_DIM), QC(U_DIM,U_DIM), MC(X_DIM,Q_DIM), Mcar(C_DIM,Q_DIM);
		CasADi::SXMatrix Acar(C_DIM,C_DIM);

		QC(0,0) = config::VELOCITY_NOISE*config::VELOCITY_NOISE;
		QC(1,1) = config::TURNING_NOISE*config::TURNING_NOISE;

		CasADi::SXMatrix x_t = X[t];

		CasADi::SXMatrix s = sin(u_t(1)+x_t(2));
		CasADi::SXMatrix c = cos(u_t(1)+x_t(2));
		CasADi::SXMatrix vts = u_t(0)*DT*s;
		CasADi::SXMatrix vtc = u_t(0)*DT*c;

		Mcar(0, 0) = DT*c;
		Mcar(0, 1) = -vts;
		Mcar(1, 0) = DT*s;
		Mcar(1, 1) = vtc;
		Mcar(2, 0) = DT*sin(u_t(1))/config::WHEELBASE;
		Mcar(2, 1) = u_t(0)*DT*cos(u_t(1))/config::WHEELBASE;

		Acar(0,0) = 1;
		Acar(1,1) = 1;
		Acar(2,2) = 1;
		Acar(0,2) = -vts;
		Acar(1,2) = vtc;

		A = CasADi::SXMatrix(CasADi::DMatrix::eye(X_DIM));
		A(CasADi::Slice(0,C_DIM),CasADi::Slice(0,C_DIM)) = Acar;
		MC(CasADi::Slice(0,C_DIM), CasADi::Slice(0,2)) = Mcar;

//		linearizeCarDynamics(X[t], u_t, QC, Acar, MMTcar);
		CasADi::SXMatrix A_inv = solve(A, CasADi::SXMatrix(CasADi::DMatrix::eye(X_DIM)));

		X[t+1] = dynfunc(X[t], u_t);

		CasADi::SXMatrix H(Z_DIM,X_DIM), R_inv(Z_DIM,Z_DIM);

		linearizeObservation(X[t+1], H);

		CasADi::SXMatrix delta = deltaMatrix(X[t+1]);
		for(int i=0; i < R_DIM-1; i += 2) {
			R_inv(i,i) = 1/(config::OBS_DIST_NOISE*config::OBS_DIST_NOISE);
			R_inv(i+1,i+1) = 1/(config::OBS_ANGLE_NOISE*config::OBS_ANGLE_NOISE);
		}

		B_ham[t+1] = mul(A, B_ham[t]) + mul(mul(mul(mul(MC, QC), trans(MC)), trans(A_inv)), C_ham[t]);
		C_ham[t+1] = mul(trans(A_inv), C_ham[t]) + mul(mul(mul(mul(mul(trans(H), delta), R_inv), delta), H), B_ham[t+1]);

		cost += params[1]*inner_prod(u_t, u_t);
	}

	CasADi::SXMatrix SigmaFinalInfo = mul(B_ham[T-1], solve(C_ham[T-1], CasADi::SXMatrix(CasADi::DMatrix::eye(X_DIM))));

	cost += params[2]*trace(SigmaFinalInfo) + params[3]*inner_prod(X[T-1] - xGoal, X[T-1] - xGoal);

	return cost;
}


CasADi::SXFunction casadiCostFunc() {
	int nU = (T-1)*U_DIM;
	CasADi::SXMatrix U = CasADi::ssym("U",nU,1);
	CasADi::SXMatrix x0 = CasADi::ssym("x0",X_DIM,1);
	CasADi::SXMatrix Sigma_0 = CasADi::ssym("S0",X_DIM,X_DIM);
	CasADi::SXMatrix xGoal = CasADi::ssym("xGoal",X_DIM,1);
	CasADi::SXMatrix params = CasADi::ssym("params",4); // alpha_control, alpha_belief, alpha_final_belief, alpha_goal_state

	// Objective
	CasADi::SXMatrix f = costfunc(U, x0, Sigma_0, xGoal, params);

	// Create functions
	std::vector<CasADi::SXMatrix> inp;
	inp.push_back(U);
	inp.push_back(x0);
	inp.push_back(Sigma_0);
	inp.push_back(xGoal);
	inp.push_back(params);

	CasADi::SXFunction f_fcn(inp,f);
	f_fcn.init();

	return f_fcn;
}

CasADi::SXFunction casadiCostGradFunc() {
	int nU = (T-1)*U_DIM;
	CasADi::SXMatrix U = CasADi::ssym("U",nU,1);
	CasADi::SXMatrix x0 = CasADi::ssym("x0",X_DIM,1);
	CasADi::SXMatrix Sigma_0 = CasADi::ssym("S0",X_DIM,X_DIM);
	CasADi::SXMatrix xGoal = CasADi::ssym("xGoal",X_DIM,1);
	CasADi::SXMatrix params = CasADi::ssym("params",4); // alpha_control, alpha_belief, alpha_final_belief, alpha_goal_state

	// Objective
	CasADi::SXMatrix f = costfunc(U, x0, Sigma_0, xGoal, params);

	CasADi::SXMatrix grad_f = gradient(f,U);

	// Create functions
	std::vector<CasADi::SXMatrix> inp;
	inp.push_back(U);
	inp.push_back(x0);
	inp.push_back(Sigma_0);
	inp.push_back(xGoal);
	inp.push_back(params);

	std::vector<CasADi::SXMatrix> out;
	out.push_back(f);
	out.push_back(grad_f);
	CasADi::SXFunction grad_f_fcn(inp,out);
	grad_f_fcn.init();

	return grad_f_fcn;
}

CasADi::SXFunction casadiCostFuncInfo() {
	int nU = (T-1)*U_DIM;
	CasADi::SXMatrix U = CasADi::ssym("U",nU,1);
	CasADi::SXMatrix x0 = CasADi::ssym("x0",X_DIM,1);
	CasADi::SXMatrix Sigma_0 = CasADi::ssym("S0",X_DIM,X_DIM);
	CasADi::SXMatrix xGoal = CasADi::ssym("xGoal",X_DIM,1);
	CasADi::SXMatrix params = CasADi::ssym("params",4); // alpha_control, alpha_belief, alpha_final_belief, alpha_goal_state

	// Objective
	CasADi::SXMatrix f = costfuncinfo(U, x0, Sigma_0, xGoal, params);

	// Create functions
	std::vector<CasADi::SXMatrix> inp;
	inp.push_back(U);
	inp.push_back(x0);
	inp.push_back(Sigma_0);
	inp.push_back(xGoal);
	inp.push_back(params);

	CasADi::SXFunction f_fcn(inp,f);
	f_fcn.init();

	return f_fcn;
}

CasADi::SXFunction casadiCostGradFuncInfo() {
	int nU = (T-1)*U_DIM;
	CasADi::SXMatrix U = CasADi::ssym("U",nU,1);
	CasADi::SXMatrix x0 = CasADi::ssym("x0",X_DIM,1);
	CasADi::SXMatrix Sigma_0 = CasADi::ssym("S0",X_DIM,X_DIM);
	CasADi::SXMatrix xGoal = CasADi::ssym("xGoal",X_DIM,1);
	CasADi::SXMatrix params = CasADi::ssym("params",4); // alpha_control, alpha_belief, alpha_final_belief, alpha_goal_state

	// Objective
	CasADi::SXMatrix f = costfuncinfo(U, x0, Sigma_0, xGoal, params);

	CasADi::SXMatrix grad_f = gradient(f,U);

	// Create functions
	std::vector<CasADi::SXMatrix> inp;
	inp.push_back(U);
	inp.push_back(x0);
	inp.push_back(Sigma_0);
	inp.push_back(xGoal);
	inp.push_back(params);

	std::vector<CasADi::SXMatrix> out;
	out.push_back(f);
	out.push_back(grad_f);
	CasADi::SXFunction grad_f_fcn(inp,out);
	grad_f_fcn.init();

	return grad_f_fcn;
}


