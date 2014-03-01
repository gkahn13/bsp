#include <iostream>
#include <iomanip>
#include <fstream>
#include <ctime>
#include <symbolic/casadi.hpp>
#include <symbolic/stl_vector_tools.hpp>
#include <cstdlib>

#include "slam.h"

// params[0] = alpha_belief
// params[1] = alpha_control
// params[2] = alpha_final_belief
// params[3] = alpha_goal_state
SXMatrix costfunc(const SXMatrix& U, const SXMatrix& x0, const SXMatrix& Sigma_0, const SXMatrix& xGoal, const SXMatrix& params)
{
	SXMatrix cost = 0;

	SXMatrix x_tp1(X_DIM,1);
	SXMatrix Sigma_t = Sigma_0, Sigma_tp1(X_DIM,X_DIM);
	SXMatrix x_t = x0, u_t(U_DIM,1);

	int offset = 0;

	for (int t = 0; t < (T-1); ++t)
	{
		u_t = U(Slice(offset,offset+U_DIM));
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


int main(int argc, char* argv[])
{

	if (argc > 1) {
		T = atoi(argv[1]);
	} else {
		T = 15;
	}

	cout << "Creating casadi file for T = " << T << endl;

	int nU = (T-1)*U_DIM;
	SXMatrix U = ssym("U",nU,1);
	SXMatrix x0 = ssym("x0",X_DIM,1);
	SXMatrix Sigma_0 = ssym("S0",X_DIM,X_DIM);
	SXMatrix xGoal = ssym("xGoal",X_DIM,1);
	SXMatrix params = ssym("params",4); // alpha_control, alpha_belief, alpha_final_belief, alpha_goal_state

	// Objective
	SXMatrix f = costfunc(U, x0, Sigma_0, xGoal, params);

	SXMatrix grad_f = gradient(f,U);

	// Create functions
	vector<SXMatrix> inp;
	inp.push_back(U);
	inp.push_back(x0);
	inp.push_back(Sigma_0);
	inp.push_back(xGoal);
	inp.push_back(params);

	SXFunction f_fcn(inp,f);
	f_fcn.init();

	vector<SXMatrix> out;
	out.push_back(f);
	out.push_back(grad_f);
	SXFunction grad_f_fcn(inp,out);
	grad_f_fcn.init();


	// Generate code
	generateCode(f_fcn,"slam-control-cost");
	generateCode(grad_f_fcn,"slam-control-grad");


	return 0;
}
