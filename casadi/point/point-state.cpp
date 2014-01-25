#include <iostream>
#include <iomanip>
#include <fstream>
#include <ctime>
#include <symbolic/casadi.hpp>
#include <symbolic/stl_vector_tools.hpp>

#define XDIM 2
#define UDIM 2
#define ZDIM 2
#define QDIM 2
#define RDIM 2
const int T = 15;
const double DT = 1;

using namespace CasADi;
using namespace std;

// dynfunc and obsfunc moved to point-dynobsjac for analytical Jacobians
SXMatrix dynfunc(const SXMatrix& x_t, const SXMatrix& u_t)
{
  	SXMatrix x_tp1 = x_t + u_t*DT;
  	return x_tp1;
}

void EKF(const SXMatrix& x_t, const SXMatrix& u_t, const SXMatrix& Sigma_t, SXMatrix& x_tp1, SXMatrix& Sigma_tp1)
{
	// hand-coded derivatives (annoying!)
	// see point-dynobsjac.cpp

	//SXMatrix A = jacobian(dynfunc(X[t],U[t],q),X[t]);
	//SXMatrix M = jacobian(dynfunc(X[t],U[t],q),q);

	SXMatrix A(XDIM,XDIM), M(XDIM,QDIM);
	A(0,0) = 1; A(1,1) = 1;
	M(0,0) = 0.01; M(1,1) = 0.01;

	//Sigma_tp1 = mul(mul(A,Sigma_t),trans(A)) + mul(M,trans(M));
	// only because A and M are diagonal matrices
	Sigma_tp1 = A*Sigma_t*trans(A) + M*trans(M);

	x_tp1 = dynfunc(x_t, u_t);

	//cout << "Observation Jacobians" << endl;
	//SXMatrix H = jacobian(obsfunc(x_tp1,r), x_tp1);
	//SXMatrix N = jacobian(obsfunc(x_tp1,r), x_tp1);

	//SXMatrix H = jacobian(obsfunc(X[t+1],r), X[t+1]);
	//SXMatrix N = jacobian(obsfunc(X[t+1],r), r);

	SXMatrix H(ZDIM,XDIM), N(ZDIM,RDIM);
	//H(0,0) = (1+(r(0)*((0.5*(0.5*(x_tp1(0)+x_tp1(0))))/(sqrt(((0.5*(0.5*sq(x_tp1(0))))+1e-06))+sqrt(((0.5*(0.5*sq(x_tp1(0))))+1e-06))))));
	//H(1,0) = (r(1)*((0.5*(0.5*(x_tp1(0)+x_tp1(0))))/(sqrt(((0.5*(0.5*sq(x_tp1(0))))+1e-06))+sqrt(((0.5*(0.5*sq(x_tp1(0))))+1e-06)))));
	// Evaluated at r = 0, simplifies to:
	H(0,0) = 1;
	H(1,1) = 1;

	//N(0,0) = sqrt(((0.5*(0.5*sq(x_tp1(0))))+1e-06));
	//N(1,1) = sqrt(((0.5*(0.5*sq(x_tp1(0))))+1e-06));
	SXMatrix intensity = sqrt(x_tp1(0) * x_tp1(0) * 0.5 * 0.5 + 1e-6);
	N(0,0) = intensity;
	N(1,1) = intensity;

    SXMatrix K = mul(mul(Sigma_tp1,trans(H)),solve(mul(mul(H,Sigma_tp1),trans(H)) + mul(N,trans(N)),SXMatrix(DMatrix::eye(ZDIM))));
	Sigma_tp1 = Sigma_tp1 - mul(K,mul(H,Sigma_tp1));
}

// params[0] = alpha_belief
// params[1] = alpha_control
// params[2] = alpha_final_belief
SXMatrix costfunc(const std::vector<SXMatrix>& X, const std::vector<SXMatrix>& U, const SXMatrix& Sigma_0, const SXMatrix& params)
{
	SXMatrix cost = 0;

	SXMatrix x_tp1(XDIM,1);
	SXMatrix Sigma_t = Sigma_0, Sigma_tp1(XDIM,XDIM);

	for (int t = 0; t < (T-1); ++t)
	{
		cost += params[0]*trace(Sigma_t);
		cost += params[1]*inner_prod(U[t], U[t]);

		EKF(X[t], U[t], Sigma_t, x_tp1, Sigma_tp1);
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

int main(){

	// Optimization variables
	vector<SXMatrix> X, U;

	for(int t = 0; t < T-1; ++t) {
		stringstream xsstm;
		xsstm << "x" << t;
		X.push_back(ssym(xsstm.str(),XDIM,1));
		stringstream usstm;
		usstm << "u" << t;
		U.push_back(ssym(usstm.str(),UDIM,1));
		//cout << usstm.str() << endl;
	}
	stringstream xsstm;
	xsstm << "x" << (T-1);
	X.push_back(ssym(xsstm.str(),XDIM,1));

	SXMatrix Sigma_0 = ssym("S0",XDIM,XDIM);
	SXMatrix params = ssym("p",3);

	// Objective
	SXMatrix f = costfunc(X, U, Sigma_0, params);

	SXMatrix grad_f(T*XDIM+(T-1)*UDIM,1), hess_f(T*XDIM+(T-1)*UDIM,1);
	int offset = 0;
	for(int t = 0; t < T-1; ++t) {
		SXMatrix dfx = gradient(f,X[t]);
		SXMatrix d2fx2 = diag(jacobian(dfx,X[t]));
		if (dfx.size1() == XDIM) {
			for(int i = 0; i < XDIM; ++i) {
				grad_f(offset+i) = dfx(i);
			}
		}
		if (d2fx2.size1() == XDIM) {
			for(int i = 0; i < XDIM; ++i) {
				hess_f(offset+i) = d2fx2(i);
			}
		}
		offset += XDIM;

		SXMatrix dfu = gradient(f,U[t]);
		SXMatrix d2fu2 = diag(jacobian(dfu,U[t]));

		if (dfu.size1() == UDIM) {
			for(int i = 0; i < UDIM; ++i) {
				grad_f(offset+i) = dfu(i);
			}
		}
		if (d2fu2.size1() == UDIM) {
			for(int i = 0; i < UDIM; ++i) {
				hess_f(offset+i) = d2fu2(i);
			}
		}

		offset += UDIM;
	}
	SXMatrix dfx = gradient(f,X[T-1]);
	SXMatrix d2fx2 = diag(jacobian(dfx,X[T-1]));

	if (dfx.size1() == XDIM) {
		for(int i = 0; i < XDIM; ++i) {
			grad_f(offset+i) = dfx(i);
		}
	}
	if (d2fx2.size1() == XDIM) {
		for(int i = 0; i < XDIM; ++i) {
			hess_f(offset+i) = d2fx2(i);
		}
	}

	// Create functions
	vector<SXMatrix> inp;
	for(int t = 0; t < (T-1); ++t) {
		inp.push_back(X[t]);
		inp.push_back(U[t]);
	}
	inp.push_back(X[T-1]);
	inp.push_back(Sigma_0);
	inp.push_back(params);
	SXFunction f_fcn(inp,f);
	f_fcn.init();
	SXFunction grad_f_fcn(inp,grad_f);
	grad_f_fcn.init();
	SXFunction hess_f_fcn(inp,hess_f);
	hess_f_fcn.init();

	// Generate code
	generateCode(f_fcn,"f");
	generateCode(grad_f_fcn,"grad_f");
	generateCode(hess_f_fcn,"hess_f");

	return 0;
}
