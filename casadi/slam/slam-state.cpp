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
SXMatrix costfunc(const SXMatrix& XU, const SXMatrix& Sigma_0, const SXMatrix& landmarks, const SXMatrix& params)
{
	SXMatrix cost = 0;

	SXMatrix x_tp1(X_DIM,1);
	SXMatrix Sigma_t = Sigma_0, Sigma_tp1(X_DIM,X_DIM);
	SXMatrix c_t(C_DIM,1), x_t(X_DIM,1), u_t(U_DIM,1);

	int offset = 0;

	for (int t = 0; t < (T-1); ++t)
	{
		c_t = XU(Slice(offset,offset+C_DIM));
		offset += C_DIM;
		u_t = XU(Slice(offset,offset+U_DIM));
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

void test() {
	T = 15;

	SXMatrix x_t(X_DIM,1), u_t(U_DIM,1), Sigma_t(X_DIM,X_DIM), x_tp1(X_DIM,1), Sigma_tp1(X_DIM,X_DIM);

	for(int i=0; i < X_DIM; ++i) { x_t(i) = 0; }
	x_t(3) = 30; x_t(4) = -10;
	x_t(5) = 70; x_t(6) = 12.5;
	x_t(7) = 20; x_t(8) = 10;

	u_t(0) = (60.0 - 0.0)/((double)(T-1));
	u_t(1) = 0;

	cout << "x0: " << x_t << endl;
	cout << "u0: " << u_t << endl;

	for(int i = 0; i < C_DIM; ++i) { Sigma_t(i,i) = .1*.1; }
	for(int i = 0; i < L_DIM; ++i) { Sigma_t(C_DIM+i,C_DIM+i) = 1*1; }

	for(int t=0; t < T-1; ++t) {
		cout << "\n\n\n\nt: " << t << endl;
		EKF(x_t, u_t, Sigma_t, x_tp1, Sigma_tp1);
		x_t = x_tp1;
		Sigma_t = Sigma_tp1;
	}

	//cout << "Sigma_t" << endl << Sigma_t << endl;
	//cout << "x_tp1" << endl << x_tp1 << endl;
	//cout << "Sigma_tp1" << endl << Sigma_tp1 << endl;
}

int main(int argc, char* argv[])
{

	if (argc > 1) {
		T = atoi(argv[1]);
	} else {
		T = 15;
	}

	cout << "Creating casadi file for T = " << T << endl;

	vector<SXMatrix> X, U;
	int nXU = T*C_DIM+(T-1)*U_DIM;
	SXMatrix XU = ssym("XU",nXU,1);
	SXMatrix Sigma_0 = ssym("S0",X_DIM,X_DIM);
	SXMatrix landmarks = ssym("landmarks",L_DIM);
	SXMatrix params = ssym("params",3); // alpha_control, alpha_belief, alpha_final_belief

	// Objective
	SXMatrix f = costfunc(XU, Sigma_0, landmarks, params);

	SXMatrix grad_f = gradient(f,XU);

	// Create functions
	vector<SXMatrix> inp;
	inp.push_back(XU);
	inp.push_back(Sigma_0);
	inp.push_back(landmarks);
	inp.push_back(params);

	SXFunction f_fcn(inp,f);
	f_fcn.init();

	vector<SXMatrix> out;
	out.push_back(f);
	out.push_back(grad_f);
	SXFunction grad_f_fcn(inp,out);
	grad_f_fcn.init();


	// Generate code
	generateCode(f_fcn,"slam-state-cost");
	generateCode(grad_f_fcn,"slam-state-grad");

//#define TEST
#ifdef TEST
	// test evaluate function
	double x0[X_DIM], xGoal[X_DIM];
	double Sigma0[X_DIM][X_DIM];

	x0[0] = 0; x0[1] = 0; x0[2] = 0;
	x0[3] = 30; x0[4] = -10;
	x0[5] = 70; x0[6] = 12.5;
	x0[7] = 20; x0[8] = 10;

	xGoal[0] = 60; xGoal[1] = 0; xGoal[2] = 0;
	xGoal[3] = 30; xGoal[4] = -10;
	xGoal[5] = 70; xGoal[6] = 12.5;
	xGoal[7] = 20; xGoal[8] = 10;

	for(int i = 0; i < X_DIM; ++i) {
		for(int j = 0; j < X_DIM; ++j) {
			Sigma0[i][j] = 0;
		}
	}
	for(int i = 0; i < C_DIM; ++i) { Sigma0[i][i] = .1*.1; }
	for(int i = 0; i < L_DIM; ++i) { Sigma0[C_DIM+i][C_DIM+i] = 1*1; }


	double t_x0[nXU];
	double t_x1[X_DIM*X_DIM];
	double t_x2[3];

	//double t_r0[1];

	double u[X_DIM];
	for(int i = 0; i < X_DIM; ++i) {
		t_x0[i] = x0[i];
		if (i < U_DIM) {
			u[i] = (xGoal[i] - x0[i])/(double)(T-1);
		}
		//cout << u[i] << endl;
	}
	//cout << endl;

	int offset = X_DIM;
	for(int t = 0; t < (T-1); ++t) {
		for(int i = 0; i < U_DIM; ++i) {
			t_x0[offset++] = u[i];
		}
		offset += X_DIM;
	}

	offset = X_DIM+U_DIM;
	for(int i = 0; i < (T-1); ++i) {
		for(int i = 0; i < X_DIM; ++i) {
			// going in straight line horizontally
			if (i == 0) {
				t_x0[offset] = t_x0[offset-(X_DIM+U_DIM)] + u[i]*DT;
			}
			else {
				t_x0[offset] = t_x0[offset-(X_DIM+U_DIM)];
			}
			offset++;
		}
		offset += U_DIM;
	}

	for(int i = 0; i < X_DIM; ++i) {
		for(int j = 0; j < X_DIM; ++j) {
			t_x1[X_DIM*i+j] = Sigma0[i][j];
		}
	}

	t_x2[0] = 10; t_x2[1] = 1; t_x2[2] = 10;


	int index = 0;
	for(int t=0; t < T-1; ++t) {
		cout << "x[" << t << "]: ";
		for(int i=0; i < X_DIM; ++i) {
			cout << t_x0[index++] << " ";
		}
		cout << endl;

		cout << "u[" << t << "]: ";
		for(int i=0; i < U_DIM; ++i) {
			cout << t_x0[index++] << " ";
		}
		cout << endl;
	}
	cout << "x[" << T-1 << "]: ";
	for(int i=0; i < X_DIM; ++i) {
		cout << t_x0[index++] << " ";
	}
	cout << endl;

	/*
	index = 0;
	cout << "SqrtSigma0\n";
	for(int i=0; i < X_DIM; ++i) {
		for(int j=0; j < X_DIM; ++j) {
			cout << t_x1[index++] << " ";
		}
		cout << "\n";
	}
	*/

	f_fcn.setInput(t_x0,0);
	f_fcn.setInput(t_x1,1);
	f_fcn.setInput(t_x2,2);
	f_fcn.evaluate();

	double cost;
	f_fcn.getOutput(&cost,0);
	cout << "cost: " << setprecision(12) << cost << endl;
#endif


	return 0;
}
