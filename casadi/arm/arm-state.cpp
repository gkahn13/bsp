#include <iostream>
#include <iomanip>
#include <fstream>
#include <ctime>
#include <symbolic/casadi.hpp>
#include <symbolic/stl_vector_tools.hpp>

#define GDIM 3
#define XDIM 6
#define UDIM 6
#define ZDIM 4
#define QDIM 6
#define RDIM 4
const int T = 10;
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
	SXMatrix sx(XDIM,1), cx(XDIM,1);
	sx = sin(x);
	cx = cos(x);

	SXMatrix p(GDIM,1);
	p(0) = sx(0)*(cx(1)*(sx(2)*(cx(4)*l4+l3)+cx(2)*cx(3)*sx(4)*l4)+sx(1)*(cx(2)*(cx(4)*l4+l3)-sx(2)*cx(3)*sx(4)*l4+l2))+cx(0)*sx(3)*sx(4)*l4;
	p(1) = -sx(1)*(sx(2)*(cx(4)*l4+l3)+cx(2)*cx(3)*sx(4)*l4)+cx(1)*(cx(2)*(cx(4)*l4+l3)-sx(2)*cx(3)*sx(4)*l4+l2)+l1;
	p(2) = cx(0)*(cx(1)*(sx(2)*(cx(4)*l4+l3)+cx(2)*cx(3)*sx(4)*l4)+sx(1)*(cx(2)*(cx(4)*l4+l3)-sx(2)*cx(3)*sx(4)*l4+l2))-sx(0)*sx(3)*sx(4)*l4;

	return p;
}

void EKF(const SXMatrix& x_t, const SXMatrix& u_t, const SXMatrix& Sigma_t, const SXMatrix& cam0, const SXMatrix& cam1, SXMatrix& x_tp1, SXMatrix& Sigma_tp1)
{
	// hand-coded derivatives (annoying!)
	// see point-dynobsjac.cpp

	//SXMatrix A = jacobian(dynfunc(X[t],U[t],q),X[t]);
	//SXMatrix M = jacobian(dynfunc(X[t],U[t],q),q);

	SXMatrix A(XDIM,XDIM), M(XDIM,QDIM);
	for(int i = 0; i < XDIM; ++i) {
		A(i,i) = 1;
		M(i,i) = 0.5; // hard-coded noise scaling constant here for now
	}

	//Sigma_tp1 = mul(mul(A,Sigma_t),trans(A)) + mul(M,trans(M));
	// because we know that A is identity
	Sigma_tp1 = A*Sigma_t*trans(A) + M*trans(M);

	x_tp1 = dynfunc(x_t, u_t);

	//cout << "Observation Jacobians" << endl;
	//SXMatrix H = jacobian(obsfunc(x_tp1,r), x_tp1);
	//SXMatrix N = jacobian(obsfunc(x_tp1,r), x_tp1);

	SXMatrix H(ZDIM,XDIM), N(ZDIM,RDIM);

	SXMatrix J(GDIM,XDIM), sx(XDIM,1), cx(XDIM,1);
	sx = sin(x_tp1);
	cx = cos(x_tp1);

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

    for (int i = 0; i < XDIM; ++i) {
    	/*
        H(0,i) = (J(0,i) * (x_tp1(2) - cam0(2)) - (x_tp1(0) - cam0(0)) * J(2,i)) / ((x_tp1(2) - cam0(2)) * (x_tp1(2) - cam0(2)));
    	H(1,i) = (J(1,i) * (x_tp1(2) - cam0(2)) - (x_tp1(1) - cam0(1)) * J(2,i)) / ((x_tp1(2) - cam0(2)) * (x_tp1(2) - cam0(2)));
    	H(2,i) = (J(0,i) * (x_tp1(2) - cam1(2)) - (x_tp1(0) - cam1(0)) * J(2,i)) / ((x_tp1(2) - cam1(2)) * (x_tp1(2) - cam1(2)));
    	H(3,i) = (J(1,i) * (x_tp1(2) - cam1(2)) - (x_tp1(1) - cam1(1)) * J(2,i)) / ((x_tp1(2) - cam1(2)) * (x_tp1(2) - cam1(2)));
        */

    	H(0,i) = (J(0,i) * (x_tp1(1) - cam0(1)) - (x_tp1(0) - cam0(0)) * J(1,i)) / ((x_tp1(1) - cam0(1)) * (x_tp1(1) - cam0(1)));
    	H(1,i) = (J(2,i) * (x_tp1(1) - cam0(1)) - (x_tp1(2) - cam0(2)) * J(1,i)) / ((x_tp1(1) - cam0(1)) * (x_tp1(1) - cam0(1)));
    	H(2,i) = (J(0,i) * (x_tp1(1) - cam1(1)) - (x_tp1(0) - cam1(0)) * J(1,i)) / ((x_tp1(1) - cam1(1)) * (x_tp1(1) - cam1(1)));
    	H(3,i) = (J(2,i) * (x_tp1(1) - cam1(1)) - (x_tp1(2) - cam1(2)) * J(1,i)) / ((x_tp1(1) - cam1(1)) * (x_tp1(1) - cam1(1)));
    }

    for(int i = 0; i < ZDIM; ++i) {
    	N(i,i) = 1;
    }

    SXMatrix K = mul(mul(Sigma_tp1,trans(H)),solve(mul(mul(H,Sigma_tp1),trans(H)) + mul(N,trans(N)),SXMatrix(DMatrix::eye(ZDIM))));
	Sigma_tp1 = Sigma_tp1 - mul(K,mul(H,Sigma_tp1));
}

// params[0] = alpha_belief
// params[1] = alpha_control
// params[2] = alpha_final_belief
SXMatrix costfunc(const SXMatrix& XU, const SXMatrix& Sigma_0, const SXMatrix& params, const SXMatrix& cam0, const SXMatrix& cam1)
{
	SXMatrix cost = 0;

	SXMatrix x_tp1(XDIM,1);
	SXMatrix Sigma_t = Sigma_0, Sigma_tp1(XDIM,XDIM);
	SXMatrix x_t(XDIM,1), u_t(UDIM,1);
	int offset = 0;
	for (int t = 0; t < (T-1); ++t)
	{
		x_t = XU(Slice(offset,offset+XDIM));
		offset += XDIM;
		u_t = XU(Slice(offset,offset+UDIM));
		offset += UDIM;

		cost += params[0]*trace(Sigma_t);
		cost += params[1]*inner_prod(u_t, u_t);

		EKF(x_t, u_t, Sigma_t, cam0, cam1, x_tp1, Sigma_tp1);
		Sigma_t = Sigma_tp1;
	}

	cost += params[2]*trace(Sigma_t);
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
	int nXU = T*XDIM+(T-1)*UDIM;
	SXMatrix XU = ssym("XU",nXU,1);
	SXMatrix Sigma_0 = ssym("S0",XDIM,XDIM);
	SXMatrix params = ssym("p",3);
	SXMatrix cam0 = ssym("cam0",GDIM);
	SXMatrix cam1 = ssym("cam1",GDIM);

	// Objective
	SXMatrix f = costfunc(XU, Sigma_0, params, cam0, cam1);

	SXMatrix grad_f = gradient(f,XU);
	SXMatrix hess_f = hessian(f,XU);

	SXMatrix diag_hess_f(nXU,1);
	for(int i = 0; i < nXU; ++i) {
		diag_hess_f(i) = hess_f(i,i);
	}

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

	SXFunction hess_f_fcn(inp,diag_hess_f);
	hess_f_fcn.init();

	// Generate code
	generateCode(f_fcn,"f");
	generateCode(grad_f_fcn,"grad_f");
	generateCode(hess_f_fcn,"hess_f");

	return 0;
}
