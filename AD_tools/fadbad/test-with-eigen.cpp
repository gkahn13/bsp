#include <iostream>

#include <Eigen/Eigen>
using namespace Eigen;

#include "badiff.h"
using namespace fadbad;

typedef B<double> bdouble;

typedef Eigen::Matrix<bdouble, 2, 1> VectorB;
typedef Eigen::Matrix<bdouble, 2, 2> MatrixB;
typedef Eigen::Matrix<bdouble, Dynamic, Dynamic > MatrixBd;
typedef Eigen::Matrix<bdouble, Dynamic, 1> VectorBd;

template <size_t _dim0, size_t _dim1>
using mat = Matrix<double, _dim0, _dim1>;

template <size_t _dim>
using vec = Matrix<double, _dim, 1>;

template <size_t _dim0, size_t _dim1>
using matb = Matrix<bdouble, _dim0, _dim1>;

template <size_t _dim>
using vecb = Matrix<bdouble, _dim, 1>;

namespace Eigen {
template<> struct NumTraits<bdouble >
: NumTraits<double> // permits to get the epsilon, dummy_precision, lowest, highest functions
{
	typedef bdouble Real;
	typedef bdouble NonInteger;
	typedef bdouble Nested;
	enum {
		IsComplex = 0,
		IsInteger = 0,
		IsSigned = 1,
		RequireInitialization = 1,
		ReadCost = 1,
		AddCost = 3,
		MulCost = 3
	};
};
}

namespace fadbad {
inline const bdouble& conj(const bdouble& x) { return x; }
inline const bdouble& real(const bdouble& x) { return x; }
inline bdouble imag(const bdouble&) { return 0.; }
inline bdouble abs(const bdouble& x) { return (x < 0) ? -x : x; }
inline bdouble abs2(const bdouble& x) { return x*x; }
inline bdouble sqrt(const bdouble& x)  { return sqrt(x); }
inline bdouble exp(const bdouble&  x)  { return exp(x); }
inline bdouble log(const bdouble&  x)  { return log(x); }
inline bdouble sin(const bdouble&  x)  { return sin(x); }
inline bdouble cos(const bdouble&  x)  { return cos(x); }
inline bdouble pow(const bdouble& x, bdouble y)  { return pow(x, y); }
}

bdouble func(const bdouble& x, const bdouble& y) {
	bdouble z=sqrt(x);
	return y*z+sin(z);
}

void test_backward() {
	bdouble x,y,f;    // Declare variables x,y,f
	x=1;                // Initialize variable x
	y=2;                // Initialize variable y
	f=func(x,y);        // Evaluate function and record DAG
	f.diff(0,1);        // Differentiate f (index 0 of 1)
	double fval=f.x();  // Value of function
	double dfdx=x.d(0); // Value of df/dx (index 0 of 1)
	double dfdy=y.d(0); // Value of df/dy (index 0 of 1)

	std::cout << "f(x,y)=" << fval << "\n";
	std::cout << "df/dx(x,y)=" << dfdx << "\n";
	std::cout << "df/dy(x,y)=" << dfdy << "\n";
}


bdouble func_eigen(const VectorB& m) {
	bdouble z=sqrt(m(0,0));
	return m(1,0)*z+sin(z);
}

void test_backward_eigen() {
	VectorB m;
	m(0,0) = 1;
	m(1,0) = 2;
	bdouble f = func_eigen(m);

	f.diff(0,1);        // Differentiate f (index 0 of 1)
	double fval=f.x();  // Value of function
	double dfdx=m(0,0).d(0); // Value of df/dx (index 0 of 1)
	double dfdy=m(1,0).d(0); // Value of df/dy (index 0 of 1)

	std::cout << "f(x,y)=" << fval << "\n";
	std::cout << "df/dx(x,y)=" << dfdx << "\n";
	std::cout << "df/dy(x,y)=" << dfdy << "\n";
}

bdouble func_matrix(const MatrixBd& m) {
//	MatrixBd m_inv = m.inverse();
	return (m*m).trace();
}

void test_matrix_func() {
	MatrixBd m(2,2);
//	m = MatrixBd::Ones(2, 2);
	m = 3*MatrixBd::Identity(2, 2);
	bdouble f = func_matrix(m);

	f.diff(0,1);        // Differentiate f (index 0 of 1)
	double fval=f.x();  // Value of function
	double dfd00 = m(0,0).d(0);
	double dfd10 = m(1,0).d(0);
	double dfd01 = m(0,1).d(0);
	double dfd11 = m(1,1).d(0);

	std::cout << "f(x,y)=" << fval << "\n";
	std::cout << dfd00 << "\t" << dfd10 << "\n";
	std::cout << dfd01 << "\t" << dfd11 << "\n";
}


// can't invert dynamic matrix?
void test_matrix_operations() {
	MatrixXd m = 3*MatrixXd::Identity(2,2);
	m.llt().solve(MatrixXd::Identity(2,2));
//	Eigen::LDLT<MatrixB> ldlt(m);
//	MatrixB m_inv = ldlt.solve(MatrixB::Identity());
//	std::cout << m_inv(0,0).x() << "\n";

//	MatrixB m = 3*MatrixB::Identity();
//	MatrixBd m_inv = m.inverse();
//	std::cout << m_inv(0,0).x() << "\n";

//	MatrixB m = 3*MatrixB::Identity();
//	MatrixB m_inv = m.inverse();
//	std::cout << m_inv(0,0).x() << "\n";

//	MatrixXd m = 2*MatrixXd::Identity(2, 2);
//	std::cout << "m:\n" << m << "\n";
//	MatrixXd m_inv = m.inverse();
//	std::cout << "m_inv:\n" << m_inv << "\n";
}

template <typename MAT, typename T>
void func_template(const MAT& m, T &f) {
	MAT m_inv = m.inverse();
	f = m_inv.sum();
}

void test_template() {
//	mat<2,2> m = 3*mat<2,2>::Identity();
//	double f;
//	func_template(m, f);
//	std::cout << f << "\n";

	matb<2,2> m = 3*matb<2,2>::Identity();
	bdouble f;
	func_template(m, f);

	f.diff(0,1);        // Differentiate f (index 0 of 1)
	double fval=f.x();  // Value of function
	double dfd00 = m(0,0).d(0);
	double dfd10 = m(1,0).d(0);
	double dfd01 = m(0,1).d(0);
	double dfd11 = m(1,1).d(0);

	std::cout << "f(x,y)=" << fval << "\n";
	std::cout << dfd00 << "\t" << dfd10 << "\n";
	std::cout << dfd01 << "\t" << dfd11 << "\n";
}

int main(int argc, char* argv[]) {
//	test_backward();
//	test_backward_eigen();
//	test_matrix_func();
//	test_matrix_operations();
	test_template();
}
