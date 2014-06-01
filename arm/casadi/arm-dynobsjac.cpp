#include <iostream>
#include <iomanip>
#include <fstream>
#include <ctime>
#include <symbolic/casadi.hpp>
#include <symbolic/stl_vector_tools.hpp>

using namespace CasADi;
using namespace std;

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

SXMatrix dynfunc(SXMatrix& x_t, SXMatrix& u_t, SXMatrix& q_t)
{
	SXMatrix x_tp1 = x_t + (u_t + q_t*0.5)*DT;
	return x_tp1;
}

// Observation model
SXMatrix obsfunc(const SXMatrix& x, const SXMatrix& r, const SXMatrix& cam0, const SXMatrix& cam1)
{
	SXMatrix ee_pos = g(x);

	SXMatrix z_t(ZDIM,1);
	z_t(0) = (ee_pos(0) - cam0(0))/(ee_pos(1) - cam0(1)) + r(0);
	z_t(1) = (ee_pos(2) - cam0(2))/(ee_pos(1) - cam0(1)) + r(1);
	z_t(2) = (ee_pos(0) - cam1(0))/(ee_pos(1) - cam1(1)) + r(2);
	z_t(3) = (ee_pos(2) - cam1(2))/(ee_pos(1) - cam1(1)) + r(3);
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

int main()
{
	// Optimization variables
	SXMatrix x = ssym("x",XDIM,1);
	SXMatrix u = ssym("u",UDIM,1);
	SXMatrix q = ssym("q",QDIM,1);
	SXMatrix r = ssym("r",RDIM,1);
	SXMatrix cam0(GDIM,1), cam1(GDIM,1);

	// Objective
	SXMatrix A = jacobian(dynfunc(x,u,q),x);
	cout << "A:\n"; A.print(); cout << endl;
	SXMatrix M = jacobian(dynfunc(x,u,q),q);
	cout << "M:\n"; M.print(); cout << endl;

	SXMatrix H = jacobian(obsfunc(x,r,cam0,cam1),x);
	cout << "H:\n"; H.print(); cout << endl;
	SXMatrix N = jacobian(obsfunc(x,r,cam0,cam1),r);
	cout << "N:\n"; N.print(); cout << endl;

	vector<SXMatrix> inp;
	inp.push_back(x);
	inp.push_back(u);
	inp.push_back(q);
	inp.push_back(r);
	inp.push_back(cam0);
	inp.push_back(cam1);

	vector<SXMatrix> out;
	out.push_back(A);
	out.push_back(M);
	out.push_back(H);
	out.push_back(N);
	SXFunction f_fcn(inp,out);
	f_fcn.init();

	// Generate code
	//generateCode(f_fcn,"f");

	return 0;
}


