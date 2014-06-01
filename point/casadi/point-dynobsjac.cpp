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
const double DT = 1;

using namespace CasADi;
using namespace std;

SXMatrix dynfunc(SXMatrix& x_t, SXMatrix& u_t, SXMatrix& q_t)
{
	SXMatrix x_tp1 = x_t + u_t*DT + q_t*0.01;
	return x_tp1;
}

SXMatrix obsfunc(SXMatrix& x_t, SXMatrix& r_t)
{
	SXMatrix intensity = sqrt(x_t(0) * x_t(0) * 0.5 * 0.5 + 1e-6);
	SXMatrix z_t = x_t + r_t*intensity;
	return z_t;
}

void generateCode(FX fcn, const std::string& name){
	cout << "Generating code for " << name << endl;

	// Convert to an SXFunction (may or may not improve efficiency)
	if(is_a<MXFunction>(fcn)){
		cout << "Casting as SXFunction" << endl;
		fcn = SXFunction(shared_cast<MXFunction>(fcn));
		fcn.init();
	}

	// Generate C code
	fcn.generateCode(name + ".c");
}

int main(){

	// Optimization variables
	SXMatrix x = ssym("x",XDIM,1);
	SXMatrix u = ssym("u",UDIM,1);
	SXMatrix q = ssym("q",QDIM,1);
	SXMatrix r = ssym("r",RDIM,1);

	// Objective
	SXMatrix A = jacobian(dynfunc(x,u,q),x);
	cout << "A:\n"; A.print(); cout << endl;
	SXMatrix M = jacobian(dynfunc(x,u,q),q);
	cout << "M:\n"; M.print(); cout << endl;

	SXMatrix H = jacobian(obsfunc(x,r),x);
	cout << "H:\n"; H.print(); cout << endl;
	SXMatrix N = jacobian(obsfunc(x,r),r);
	cout << "N:\n"; N.print(); cout << endl;

	vector<SXMatrix> inp;
	inp.push_back(x);
	inp.push_back(u);
	inp.push_back(q);
	inp.push_back(r);

	vector<SXMatrix> out;
	out.push_back(A);
	out.push_back(M);
	out.push_back(H);
	out.push_back(N);
	SXFunction f_fcn(inp,out);
	f_fcn.init();

	// Generate code
	generateCode(f_fcn,"f");

	return 0;
}
