#include "casadi-slam-state.h"
#include "casadi-slam.h"


// params[0] = alpha_belief
// params[1] = alpha_control
// params[2] = alpha_final_belief
CasADi::SXMatrix costfunc(const CasADi::SXMatrix& XU, const CasADi::SXMatrix& Sigma_0, const CasADi::SXMatrix& landmarks, const CasADi::SXMatrix& params)
{
	CasADi::SXMatrix cost = 0;

	CasADi::SXMatrix x_tp1(X_DIM,1);
	CasADi::SXMatrix Sigma_t = Sigma_0, Sigma_tp1(X_DIM,X_DIM);
	CasADi::SXMatrix c_t(C_DIM,1), x_t(X_DIM,1), u_t(U_DIM,1);

	int offset = 0;

	for (int t = 0; t < (T-1); ++t)
	{
		c_t = XU(CasADi::Slice(offset,offset+C_DIM));
		offset += C_DIM;
		u_t = XU(CasADi::Slice(offset,offset+U_DIM));
		offset += U_DIM;

		for(int i=0; i < X_DIM; ++i) {
			x_t(i) = (i < C_DIM) ? c_t(i) : landmarks(i-C_DIM);
		}

		cost += params[0]*trace(Sigma_t);
		cost += params[1]*inner_prod(u_t, u_t);

		EKF(x_t, u_t, Sigma_t, x_tp1, Sigma_tp1);
		Sigma_t = Sigma_tp1;
	}

	cost += params[2]*trace(Sigma_t);

	return cost;
}


CasADi::SXFunction casadiCostFunc() {
	int nXU = T*C_DIM+(T-1)*U_DIM;
	CasADi::SXMatrix XU = CasADi::ssym("XU",nXU,1);
	CasADi::SXMatrix Sigma_0 = CasADi::ssym("S0",X_DIM,X_DIM);
	CasADi::SXMatrix landmarks = CasADi::ssym("landmarks",L_DIM);
	CasADi::SXMatrix params = CasADi::ssym("params",3); // alpha_control, alpha_belief, alpha_final_belief

	// Objective
	CasADi::SXMatrix f = costfunc(XU, Sigma_0, landmarks, params);

	// Create functions
	std::vector<CasADi::SXMatrix> inp;
	inp.push_back(XU);
	inp.push_back(Sigma_0);
	inp.push_back(landmarks);
	inp.push_back(params);

	CasADi::SXFunction f_fcn(inp,f);
	f_fcn.init();

	return f_fcn;
}

CasADi::SXFunction casadiCostGradFunc() {
	int nXU = T*C_DIM+(T-1)*U_DIM;
	CasADi::SXMatrix XU = CasADi::ssym("XU",nXU,1);
	CasADi::SXMatrix Sigma_0 = CasADi::ssym("S0",X_DIM,X_DIM);
	CasADi::SXMatrix landmarks = CasADi::ssym("landmarks",L_DIM);
	CasADi::SXMatrix params = CasADi::ssym("params",3); // alpha_control, alpha_belief, alpha_final_belief

	// Objective
	CasADi::SXMatrix f = costfunc(XU, Sigma_0, landmarks, params);

	CasADi::SXMatrix grad_f = gradient(f,XU);

	// Create functions
	std::vector<CasADi::SXMatrix> inp;
	inp.push_back(XU);
	inp.push_back(Sigma_0);
	inp.push_back(landmarks);
	inp.push_back(params);

	std::vector<CasADi::SXMatrix> out;
	out.push_back(f);
	out.push_back(grad_f);
	CasADi::SXFunction grad_f_fcn(inp,out);
	grad_f_fcn.init();

	return grad_f_fcn;
}

