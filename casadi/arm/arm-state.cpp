#include <iostream>
#include <iomanip>
#include <fstream>
#include <ctime>
#include <symbolic/casadi.hpp>
#include <symbolic/stl_vector_tools.hpp>

#define G_DIM 3
#define X_DIM 6
#define U_DIM 6
#define Z_DIM 4
#define Q_DIM 6
#define R_DIM 4
const int T = 15;
const double DT = 1;
const double l1 = 7.25;
const double l2 = 8;
const double l3 = 10.375;
const double l4 = 2.375;

using namespace CasADi;
using namespace std;

// dynfunc and obsfunc moved to arm-dynobsjac for analytical Jacobians
SXMatrix dynfunc(const SXMatrix& x_t, const SXMatrix& u_t)
{
  	SXMatrix x_tp1 = x_t + u_t*DT;
  	return x_tp1;
}

SXMatrix g(const SXMatrix& x)
{
	SXMatrix sx(X_DIM,1), cx(X_DIM,1);
	sx = sin(x);
	cx = cos(x);

	SXMatrix p(G_DIM,1);
	p(0) = sx(0)*(cx(1)*(sx(2)*(cx(4)*l4+l3)+cx(2)*cx(3)*sx(4)*l4)+sx(1)*(cx(2)*(cx(4)*l4+l3)-sx(2)*cx(3)*sx(4)*l4+l2))+cx(0)*sx(3)*sx(4)*l4;
	p(1) = -sx(1)*(sx(2)*(cx(4)*l4+l3)+cx(2)*cx(3)*sx(4)*l4)+cx(1)*(cx(2)*(cx(4)*l4+l3)-sx(2)*cx(3)*sx(4)*l4+l2)+l1;
	p(2) = cx(0)*(cx(1)*(sx(2)*(cx(4)*l4+l3)+cx(2)*cx(3)*sx(4)*l4)+sx(1)*(cx(2)*(cx(4)*l4+l3)-sx(2)*cx(3)*sx(4)*l4+l2))-sx(0)*sx(3)*sx(4)*l4;

	return p;
}

SXMatrix linearizeg(const SXMatrix& x)
{
	SXMatrix J(G_DIM,X_DIM);

	SXMatrix sx(X_DIM,1), cx(X_DIM,1);
	sx = sin(x);
	cx = cos(x);

    J(0,0) = cx(0)*(cx(1)*(sx(2)*(cx(4)*l4+l3)+cx(2)*cx(3)*sx(4)*l4)+sx(1)*(cx(2)*(cx(4)*l4+l3)-sx(2)*cx(3)*sx(4)*l4+l2))-sx(0)*sx(3)*sx(4)*l4;
    J(0,1) = sx(0)*(cx(1)*(cx(2)*(cx(4)*l4+l3)-sx(2)*cx(3)*sx(4)*l4+l2)-sx(1)*(sx(2)*(cx(4)*l4+l3)+cx(2)*cx(3)*sx(4)*l4));
    J(0,2) = sx(0)*(sx(1)*(-sx(2)*(cx(4)*l4+l3)-cx(2)*cx(3)*sx(4)*l4)+cx(1)*(cx(2)*(cx(4)*l4+l3)-sx(2)*cx(3)*sx(4)*l4));
    J(0,3) = sx(0)*(sx(1)*sx(2)*sx(3)*sx(4)*l4-cx(1)*cx(2)*sx(3)*sx(4)*l4)+cx(0)*cx(3)*sx(4)*l4;
    J(0,4) = sx(0)*(cx(1)*(cx(2)*cx(3)*cx(4)*l4-sx(2)*sx(4)*l4)+sx(1)*(-cx(2)*sx(4)*l4-sx(2)*cx(3)*cx(4)*l4))+cx(0)*sx(3)*cx(4)*l4;
    J(0,5) = 0;

    J(1,0) = 0;
    J(1,1) = -cx(1)*(sx(2)*(cx(4)*l4+l3)+cx(2)*cx(3)*sx(4)*l4)-sx(1)*(cx(2)*(cx(4)*l4+l3)-sx(2)*cx(3)*sx(4)*l4+l2);
    J(1,2) = cx(1)*(-sx(2)*(cx(4)*l4+l3)-cx(2)*cx(3)*sx(4)*l4)-sx(1)*(cx(2)*(cx(4)*l4+l3)-sx(2)*cx(3)*sx(4)*l4);
    J(1,3) = cx(1)*sx(2)*sx(3)*sx(4)*l4+sx(1)*cx(2)*sx(3)*sx(4)*l4;
    J(1,4) = cx(1)*(-cx(2)*sx(4)*l4-sx(2)*cx(3)*cx(4)*l4)-sx(1)*(cx(2)*cx(3)*cx(4)*l4-sx(2)*sx(4)*l4);
    J(1,5) = 0;

    J(2,0) = -sx(0)*(cx(1)*(sx(2)*(cx(4)*l4+l3)+cx(2)*cx(3)*sx(4)*l4)+sx(1)*(cx(2)*(cx(4)*l4+l3)-sx(2)*cx(3)*sx(4)*l4+l2))-cx(0)*sx(3)*sx(4)*l4;
    J(2,1) = cx(0)*(cx(1)*(cx(2)*(cx(4)*l4+l3)-sx(2)*cx(3)*sx(4)*l4+l2)-sx(1)*(sx(2)*(cx(4)*l4+l3)+cx(2)*cx(3)*sx(4)*l4));
    J(2,2) = cx(0)*(sx(1)*(-sx(2)*(cx(4)*l4+l3)-cx(2)*cx(3)*sx(4)*l4)+cx(1)*(cx(2)*(cx(4)*l4+l3)-sx(2)*cx(3)*sx(4)*l4));
    J(2,3) = cx(0)*(sx(1)*sx(2)*sx(3)*sx(4)*l4-cx(1)*cx(2)*sx(3)*sx(4)*l4)-sx(0)*cx(3)*sx(4)*l4;
    J(2,4) = cx(0)*(cx(1)*(cx(2)*cx(3)*cx(4)*l4-sx(2)*sx(4)*l4)+sx(1)*(-cx(2)*sx(4)*l4-sx(2)*cx(3)*cx(4)*l4))-sx(0)*sx(3)*cx(4)*l4;
    J(2,5) = 0;


	return J;
}

void EKF(const SXMatrix& x_t, const SXMatrix& u_t, const SXMatrix& Sigma_t, const SXMatrix& cam0, const SXMatrix& cam1, SXMatrix& x_tp1, SXMatrix& Sigma_tp1)
{
	// hand-coded derivatives (annoying!)
	// see point-dynobsjac.cpp

	//SXMatrix A = jacobian(dynfunc(X[t],U[t],q),X[t]);
	//SXMatrix M = jacobian(dynfunc(X[t],U[t],q),q);

	SXMatrix A(X_DIM,X_DIM), M(X_DIM,Q_DIM), QC(U_DIM,U_DIM);
	for(int i = 0; i < X_DIM; ++i) {
		A(i,i) = 1;
		M(i,i) = DT;
		QC(i,i) = 0.01;
	}

	//Sigma_tp1 = mul(mul(A,Sigma_t),trans(A)) + mul(M,trans(M));
	// because we know that A is identity
	Sigma_tp1 = A*Sigma_t*trans(A) + M*QC*trans(M);

	x_tp1 = dynfunc(x_t, u_t);

	//cout << "Observation Jacobians" << endl;
	//SXMatrix H = jacobian(obsfunc(x_tp1,r), x_tp1);
	//SXMatrix N = jacobian(obsfunc(x_tp1,r), x_tp1);

	SXMatrix H(Z_DIM,X_DIM), N(Z_DIM,R_DIM), RC(Z_DIM,Z_DIM);
	SXMatrix J(G_DIM,X_DIM);
	J = linearizeg(x_tp1);

    for (int i = 0; i < X_DIM; ++i) {
    	/*
        H(0,i) = (J(0,i) * (x_tp1(2) - cam0(2)) - (x_tp1(0) - cam0(0)) * J(2,i)) / ((x_tp1(2) - cam0(2)) * (x_tp1(2) - cam0(2)));
    	H(1,i) = (J(1,i) * (x_tp1(2) - cam0(2)) - (x_tp1(1) - cam0(1)) * J(2,i)) / ((x_tp1(2) - cam0(2)) * (x_tp1(2) - cam0(2)));
    	H(2,i) = (J(0,i) * (x_tp1(2) - cam1(2)) - (x_tp1(0) - cam1(0)) * J(2,i)) / ((x_tp1(2) - cam1(2)) * (x_tp1(2) - cam1(2)));
    	H(3,i) = (J(1,i) * (x_tp1(2) - cam1(2)) - (x_tp1(1) - cam1(1)) * J(2,i)) / ((x_tp1(2) - cam1(2)) * (x_tp1(2) - cam1(2)));
        */

    	H(0,i) = -(J(0,i) * (x_tp1(1) - cam0(1)) - (x_tp1(0) - cam0(0)) * J(1,i)) / ((x_tp1(1) - cam0(1)) * (x_tp1(1) - cam0(1)));
    	H(1,i) = -(J(2,i) * (x_tp1(1) - cam0(1)) - (x_tp1(2) - cam0(2)) * J(1,i)) / ((x_tp1(1) - cam0(1)) * (x_tp1(1) - cam0(1)));
    	H(2,i) = -(J(0,i) * (x_tp1(1) - cam1(1)) - (x_tp1(0) - cam1(0)) * J(1,i)) / ((x_tp1(1) - cam1(1)) * (x_tp1(1) - cam1(1)));
    	H(3,i) = -(J(2,i) * (x_tp1(1) - cam1(1)) - (x_tp1(2) - cam1(2)) * J(1,i)) / ((x_tp1(1) - cam1(1)) * (x_tp1(1) - cam1(1)));
    }

    for(int i = 0; i < Z_DIM; ++i) {
    	N(i,i) = 1;
    	RC(i,i) = 0.01;
    }

    //SXMatrix K = mul(mul(Sigma_tp1,trans(H)),solve(mul(mul(H,Sigma_tp1),trans(H)) + mul(N,trans(N)),SXMatrix(DMatrix::eye(Z_DIM))));
    SXMatrix K = mul(mul(Sigma_tp1,trans(H)),solve(mul(mul(H,Sigma_tp1),trans(H)) + N*RC*trans(N),SXMatrix(DMatrix::eye(Z_DIM))));
	Sigma_tp1 = Sigma_tp1 - mul(K,mul(H,Sigma_tp1));
}

// params[0] = alpha_belief
// params[1] = alpha_control
// params[2] = alpha_final_belief
SXMatrix costfunc(const SXMatrix& XU, const SXMatrix& Sigma_0, const SXMatrix& params, const SXMatrix& cam0, const SXMatrix& cam1)
{
	SXMatrix cost = 0;

	SXMatrix x_tp1(X_DIM,1);
	SXMatrix Sigma_t = Sigma_0, Sigma_tp1(X_DIM,X_DIM);
	SXMatrix x_t(X_DIM,1), u_t(U_DIM,1);
	SXMatrix J(G_DIM,X_DIM);

	int offset = 0;
	for (int t = 0; t < (T-1); ++t)
	{
		x_t = XU(Slice(offset,offset+X_DIM));
		offset += X_DIM;
		u_t = XU(Slice(offset,offset+U_DIM));
		offset += U_DIM;

		J = linearizeg(x_t);
		cost += params[0]*trace(mul(mul(J,Sigma_t),trans(J)));
		cost += params[1]*inner_prod(u_t, u_t);

		EKF(x_t, u_t, Sigma_t, cam0, cam1, x_tp1, Sigma_tp1);
		Sigma_t = Sigma_tp1;
	}

	x_t = XU(Slice(offset,offset+X_DIM));
	J = linearizeg(x_t);
	cost += params[2]*trace(mul(mul(J,Sigma_t),trans(J)));
	return cost;
}

void generateCode(FX fcn, const std::string& name){
	cout << "Generating code for " << name << endl;

	// Convert to an SXFunction (may or may not improve efficiency)
	//if(is_a<MXFunction>(fcn)){
	//	cout << "Casting as SXFunction" << endl;
	//	fcn = SXFunction(shared_cast<MXFunction>(fcn));
	//	fcn.init();
	//}

	// Generate C code
	fcn.generateCode(name + ".c");
}

int main()
{
	vector<SXMatrix> X, U;
	int nXU = T*X_DIM+(T-1)*U_DIM;
	SXMatrix XU = ssym("XU",nXU,1);
	SXMatrix Sigma_0 = ssym("S0",X_DIM,X_DIM);
	SXMatrix params = ssym("p",3);
	SXMatrix cam0 = ssym("cam0",G_DIM);
	SXMatrix cam1 = ssym("cam1",G_DIM);

	// Objective
	SXMatrix f = costfunc(XU, Sigma_0, params, cam0, cam1);

	SXMatrix grad_f = gradient(f,XU);

	//SXMatrix hess_f = hessian(f,XU);

	//SXMatrix diag_hess_f(nXU,1);
	//for(int i = 0; i < nXU; ++i) {
	//	diag_hess_f(i) = hess_f(i,i);
	//}

	// Create functions
	vector<SXMatrix> inp;
	inp.push_back(XU);
	inp.push_back(Sigma_0);
	inp.push_back(params);
	inp.push_back(cam0);
	inp.push_back(cam1);

	SXFunction f_fcn(inp,f);
	f_fcn.init();

	SXFunction grad_f_fcn(inp,grad_f);
	grad_f_fcn.init();

	//SXFunction hess_f_fcn(inp,diag_hess_f);
	//hess_f_fcn.init();

	// Generate code
	generateCode(f_fcn,"f");
	generateCode(grad_f_fcn,"grad_f");
	//generateCode(hess_f_fcn,"hess_f");

	return 0;
}
