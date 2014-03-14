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

