#include "util/ilqg-matrix.h"

#include "ilqg.h"

#include <time.h>


#define DT 1.0
#define X_DIM 2
#define U_DIM 2
#define Z_DIM 2
#define S_DIM (((X_DIM+1)*X_DIM)/2)
#define B_DIM (X_DIM + S_DIM)

SymmetricMatrix<U_DIM> Rint;
SymmetricMatrix<X_DIM> Qint;
SymmetricMatrix<X_DIM> QGoal;
SymmetricMatrix<X_DIM> Sigma0;

Matrix<X_DIM> x0, xGoal;

inline Matrix<X_DIM> f(const Matrix<X_DIM>& x, const Matrix<U_DIM>& u)
{  
	Matrix<X_DIM> xNew;

	xNew[0] = x[0] + u[0] * DT;
	xNew[1] = x[1] + u[1] * DT;

	return xNew;
}

// Observation model
inline Matrix<Z_DIM> h(const Matrix<X_DIM>& x)
{
	Matrix<Z_DIM> z;

	z[0] = x[0];
	z[1] = x[1];

	return z;
}

// Jacobian df/dx(x,u)
template <size_t _xDim, size_t _uDim>
inline Matrix<_xDim,_xDim> dfdx(Matrix<_xDim> (*f)(const Matrix<_xDim>&, const Matrix<_uDim>&), const Matrix<_xDim>& x, const Matrix<_uDim>& u) 
{
	Matrix<_xDim,_xDim> A;
	Matrix<_xDim> xr(x), xl(x);
	for (size_t i = 0; i < _xDim; ++i) {
		xr[i] += step; xl[i] -= step;
		A.insert(0,i, (f(xr, u) - f(xl, u)) / (2.0*step));
		xr[i] = xl[i] = x[i];
	}
	return A;
}

// Jacobian df/du(x,u)
template <size_t _xDim, size_t _uDim>
inline Matrix<_xDim,_uDim> dfdu(Matrix<_xDim> (*f)(const Matrix<_xDim>&, const Matrix<_uDim>&),	const Matrix<_xDim>& x, const Matrix<_uDim>& u) 
{
	Matrix<_xDim,_uDim> B;
	Matrix<_uDim> ur(u), ul(u);
	for (size_t i = 0; i < _uDim; ++i) {
		ur[i] += step; ul[i] -= step;
		B.insert(0,i, (f(x, ur) - f(x, ul)) / (2.0*step));
		ur[i] = ul[i] = u[i];
	}
	return B;
}

// Jacobian dh/dx(x)
template <size_t _xDim, size_t _zDim>
inline Matrix<_zDim,_xDim> dhdx(Matrix<_zDim> (*h)(const Matrix<_xDim>&), const Matrix<_xDim>& x) 
{
	Matrix<_zDim,_xDim> H;
	Matrix<_xDim> xr(x), xl(x);
	for (size_t i = 0; i < _xDim; ++i) {
		xr[i] += step; xl[i] -= step;
		H.insert(0,i, (h(xr) - h(xl)) / (2.0*step));
		xr[i] = xl[i] = x[i];
	}
	return H;
}

inline SymmetricMatrix<X_DIM> varM(const Matrix<X_DIM>& x, const Matrix<U_DIM>& u)
{
	SymmetricMatrix<X_DIM> S = identity<X_DIM>();
	S(0,0) = 0.01;
	S(1,1) = 0.01;
	return S;
}

inline SymmetricMatrix<Z_DIM> varN(const Matrix<X_DIM>& x)
{
	SymmetricMatrix<Z_DIM> S = identity<Z_DIM>();
	double intensity = sqrt(sqr(0.5*x[0]) + 1e-6);
	S(0,0) = intensity;
	S(1,1) = intensity;
	return S;
}

// Compute closed form quadratic finalCost function around around b
inline void quadratizeFinalCost(const Matrix<X_DIM>& xBar, const SymmetricMatrix<X_DIM>& SigmaBar, double& s, SymmetricMatrix<X_DIM>& S, Matrix<1,X_DIM>& sT, Matrix<1,S_DIM>& tT, unsigned int flag) 
{
	if (flag & COMPUTE_S) S = QGoal;
	if (flag & COMPUTE_sT) sT = ~(xBar - xGoal)*QGoal;
	//if (flag & COMPUTE_s) s = 0.5*scalar(~(xBar - xGoal)*QGoal*(xBar - xGoal)) + scalar(vecTh(QGoal)*vec(SigmaBar));
	if (flag & COMPUTE_s) s = 2*(0.5*scalar(~(xBar - xGoal)*QGoal*(xBar - xGoal)) + scalar(vecTh(QGoal)*vec(SigmaBar)));
	if (flag & COMPUTE_tT) tT = vecTh(QGoal);
}

// Compute closed form quadratic cost function around around b and u
inline bool quadratizeCost(const Matrix<X_DIM>& xBar, const SymmetricMatrix<X_DIM>& SigmaBar, const Matrix<U_DIM>& uBar, double& q, SymmetricMatrix<X_DIM>& Q, SymmetricMatrix<U_DIM>& R, 
						   Matrix<U_DIM, X_DIM>& P, Matrix<1,X_DIM>& qT, Matrix<1,U_DIM>& rT, Matrix<1,S_DIM>& pT, unsigned int flag) 
{
	if (flag & COMPUTE_Q) Q = zeros<X_DIM>(); // penalize uncertainty
	if (flag & COMPUTE_R) R = Rint;
	if (flag & COMPUTE_P) P = zeros<U_DIM,X_DIM>();
	if (flag & COMPUTE_qT) qT = zeros<1, X_DIM>();
	if (flag & COMPUTE_rT) rT = ~uBar*Rint;
	if (flag & COMPUTE_pT) pT = vecTh(Qint);
	//if (flag & COMPUTE_q) q = 0.5*scalar(~uBar*Rint*uBar) + scalar(vecTh(Qint)*vec(SigmaBar));
	if (flag & COMPUTE_q) q = 2*(0.5*scalar(~uBar*Rint*uBar) + scalar(vecTh(Qint)*vec(SigmaBar)));

	return true;
}

inline void linearizeDynamics(const Matrix<X_DIM>& xBar, const Matrix<U_DIM>& uBar, Matrix<X_DIM>& c, Matrix<X_DIM, X_DIM>& A, Matrix<X_DIM, U_DIM>& B, SymmetricMatrix<X_DIM>& M, unsigned int flag)
{
	if (flag & COMPUTE_c) c = f(xBar, uBar);
	if (flag & COMPUTE_A) A = dfdx(f, xBar, uBar); //identity<X_DIM>();
	if (flag & COMPUTE_B) B = dfdu(f, xBar, uBar); //identity<U_DIM>();
	if (flag & COMPUTE_M) M = varM(xBar, uBar);
}

inline void linearizeObservation(const Matrix<X_DIM>& xBar, Matrix<Z_DIM, X_DIM>& H, SymmetricMatrix<Z_DIM>& N)
{
	H = dhdx(h, xBar); // identity<X_DIM>();
	N = varN(xBar);
}

int main(int argc, char* argv[])
{
	size_t numSteps = 15;
	
	x0[0] = -3.5; x0[1] = 2;
	xGoal[0] = -3.5; xGoal[1] = -2;

	// init controls from start to goal -- straight line trajectory
	std::vector< Matrix<U_DIM> > uBar(numSteps, (xGoal - x0)/(double)numSteps); 

	Rint = identity<U_DIM>();
	Qint = 10.0*identity<X_DIM>();
	QGoal = 10.0*identity<X_DIM>();

	std::vector< Matrix<U_DIM, X_DIM> > L;
	std::vector< Matrix<X_DIM> > xBar;
	std::vector< SymmetricMatrix<X_DIM> > SigmaBar;
	
	Sigma0 = identity<X_DIM>();
	xBar.push_back(x0);
	SigmaBar.push_back(Sigma0);

	

	solvePOMDP(linearizeDynamics, linearizeObservation, quadratizeFinalCost, quadratizeCost, xBar, SigmaBar, uBar, L);

	//for (size_t i = 0; i < xBar.size(); ++i) {
	//	std::cout << ~(xBar[i]);
	//}

	int k;
	std::cin >> k;

	return 0;
}

