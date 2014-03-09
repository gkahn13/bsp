#include <iostream>
#include <iomanip>
#include <fstream>
#include <ctime>
#include <symbolic/casadi.hpp>
#include <symbolic/stl_vector_tools.hpp>
#include <cstdlib>

//#define TIMESTEPS 15
#define NUM_LANDMARKS 6
#define NUM_WAYPOINTS 4

#define C_DIM 3 // car dimension [x, y, theta]
#define P_DIM 2 // Position dimension [x,y]
#define L_DIM 2*NUM_LANDMARKS

#define X_DIM (C_DIM+L_DIM)
#define U_DIM 2
#define Z_DIM L_DIM
#define Q_DIM 2
#define R_DIM L_DIM

int T;
const double DT = 1;

namespace config {
const double V = 3;
const double MAXG = 30*M_PI/180.;
const double RATEG = 20*M_PI/180.;
const double WHEELBASE = 4;
const double DT_CONTROLS = 0.025;

const double VELOCITY_NOISE = 0.3;
const double TURNING_NOISE = 3.0*M_PI/180.;

const double MAX_RANGE = 5.0;
const double DT_OBSERVE = 8*DT_CONTROLS;

const double OBS_DIST_NOISE = 1 * 0.1;
const double OBS_ANGLE_NOISE = 1 * 1.0*M_PI/180.;

const double ALPHA_OBS = .75;
}

using namespace CasADi;
using namespace std;


// dynfunc and obsfunc moved to arm-dynobsjac for analytical Jacobians
SXMatrix dynfunc(const SXMatrix& x_t, const SXMatrix& u_t)
{
	SXMatrix x_tp1 = x_t;

  	x_tp1(0) += u_t(0) * DT * cos(x_t(2) + u_t(1));
  	x_tp1(1) += u_t(0) * DT * sin(x_t(2) + u_t(1));
  	x_tp1(2) += u_t(0) * DT * sin(u_t(1))/config::WHEELBASE;

  	return x_tp1;
}

// Jacobians: df(x,u,q)/dx, df(x,u,q)/dq
void linearizeDynamics(const SXMatrix& x, const SXMatrix& u, SXMatrix& A, SXMatrix& M)
{
	//g is control input steer angle
	SXMatrix s = sin(u(1)+x(2));
	SXMatrix c= cos(u(1)+x(2));

	SXMatrix vts= u(0)*DT*s;
	SXMatrix vtc= u(0)*DT*c;

	M(0, 0) = DT*c;
	M(0, 1) = -vts;
	M(1, 0) = DT*s;
	M(1, 1) = vtc;
	M(2, 0) = DT*sin(u(1))/config::WHEELBASE;
	M(2, 1) = u(0)*DT*cos(u(1))/config::WHEELBASE;

	for(int i=0; i < X_DIM; ++i) {
		A(i,i) = 1;
	}

	A(0,0) = 1;
	A(1,1) = 1;
	A(2,2) = 1;
	A(0,2) = -vts;
	A(1,2) = vtc;
}

inline void linearizeCarDynamics(const SXMatrix& x, const SXMatrix& u, const SXMatrix& QC, SXMatrix& Acar, SXMatrix& MMTcar) {
	//g is control input steer angle
	SXMatrix ssin = sin(u(1)+x(2));
	SXMatrix ccos = cos(u(1)+x(2));

	SXMatrix vts= u(0)*DT*ssin;
	SXMatrix vtc= u(0)*DT*ccos;

	SXMatrix a = DT*ccos;
	SXMatrix b = DT*ssin;
	SXMatrix c = DT*sin(u(1))/config::WHEELBASE;
	SXMatrix d = -vts;
	SXMatrix e = vtc;
	SXMatrix f = u(0)*DT*cos(u(1))/config::WHEELBASE;

	SXMatrix alpha = QC(0,0);
	SXMatrix beta = QC(1,1);

	MMTcar(0,0) = alpha*a*a + beta*d*d;
	MMTcar(0,1) = alpha*a*b + beta*d*e;
	MMTcar(0,2) = alpha*a*c + beta*d*f;
	MMTcar(1,1) = alpha*b*b + beta*e*e;
	MMTcar(1,2) = alpha*b*c + beta*e*f;
	MMTcar(2,2) = alpha*c*c + beta*f*f;

	MMTcar(1,0) = MMTcar(0,1);
	MMTcar(2,0) = MMTcar(0,2);
	MMTcar(2,1) = MMTcar(1,2);

	Acar(0,0) = 1;
	Acar(1,1) = 1;
	Acar(2,2) = 1;
	Acar(0,2) = -vts;
	Acar(1,2) = vtc;
}

// Jacobians: dh(x,r)/dx, dh(x,r)/dr
void linearizeObservation(const SXMatrix& x, SXMatrix& H, SXMatrix& N)
{
	for (int i=0; i < L_DIM; i+=2) {
		SXMatrix dx = x(C_DIM+i) - x(0);
		SXMatrix dy = x(C_DIM+i+1) - x(1);
		SXMatrix d2 = dx*dx + dy*dy + 1e-10;
		SXMatrix d = sqrt(d2 + 1e-10);
		SXMatrix xd = dx/d;
		SXMatrix yd = dy/d;
		SXMatrix xd2 = dx/d2;
		SXMatrix yd2 = dy/d2;


		H(i, 0) = -xd;
		H(i, 1) = -yd;
		H(i, 2) = 0;
		H(i+1, 0) = yd2;
		H(i+1, 1) = -xd2;
		H(i+1, 2) = -1;
		H(i, 3+i) = xd;
		H(i, 3+i+1) = yd;
		H(i+1, 3+i) = -yd2;
		H(i+1, 3+i+1) = xd2;
	}

	for(int i=0; i < Z_DIM; ++i) {
		N(i,i) = 1.0;
	}

}

SXMatrix deltaMatrix(const SXMatrix& x) {
	SXMatrix delta(Z_DIM, Z_DIM);
	SXMatrix l0, l1, dist;
	for(int i=C_DIM; i < X_DIM; i += 2) {
		l0 = x(i);
		l1 = x(i+1);

		dist = sqrt((x(0) - l0)*(x(0) - l0) + (x(1) - l1)*(x(1) - l1));

		SXMatrix signed_dist = 1/(1+exp(-config::ALPHA_OBS*(config::MAX_RANGE-dist)));
		delta(i-C_DIM,i-C_DIM) = signed_dist;
		delta(i-C_DIM+1,i-C_DIM+1) = signed_dist;
	}

	return delta;
}

template <size_t _size, size_t _numColumns>
inline CasADi::SXMatrix operator_percent(const CasADi::SXMatrix& p, const CasADi::SXMatrix& q) {
	// Cholesky factorization p = L*~L
	CasADi::SXMatrix L(_size,_size); // abuse SymmetricMatrix for triangular matrix
	for (size_t i = 0; i < _size; ++i) {
		for (size_t j = i; j < _size; ++j) {
			CasADi::SXMatrix sum = p(j,i);
			for (size_t k = 0; k < i; ++k) {
				sum -= L(j,k)*L(i,k);
			}
			if (i == j) {
				L(i,i) = sqrt(sum);
			} else {
				L(j,i) = sum / L(i,i);
			}
		}
	}

	// Backward and forward substitution
	CasADi::SXMatrix M(_size, _numColumns);
	for (size_t i = 0; i < _size; ++i) {
		for (size_t k = 0; k < _numColumns; ++k) {
			CasADi::SXMatrix sum = q(i,k);
			for (size_t j = 0; j < i; ++j) {
				sum -= L(i,j)*M(j,k);
			}
			M(i,k) = sum / L(i,i);
		}
	}
	for (size_t i = _size - 1; i != -1; --i) {
		for (size_t k = 0; k < _numColumns; ++k) {
			CasADi::SXMatrix sum = M(i,k);
			for (size_t j = i + 1; j < _size; ++j) {
				sum -= L(j,i)*M(j,k);
			}
			M(i,k) = sum / L(i,i);
		}
	}
	return M;
}


template <size_t _size, size_t _numRows>
inline CasADi::SXMatrix operator_divide(const CasADi::SXMatrix& p, const CasADi::SXMatrix& q) {
	return trans(operator_percent<_size,_numRows>(q, trans(p)));
}

// Compute product M*N of which one knows that the results is symmetric (and save half the computation)
template <size_t _size, size_t _numRows>
inline CasADi::SXMatrix SymProd(const CasADi::SXMatrix& M, const CasADi::SXMatrix& N) {
	CasADi::SXMatrix S(_size,_size);
	for (size_t j = 0; j < _size; ++j) {
		for (size_t i = j; i < _size; ++i) {
			//CasADi::SXMatrix temp(1,1);
			for (size_t k = 0; k < _numRows; ++k) {
				S(i,j) += M(i, k) * N(k, j);
			}
			S(j, i) = S(i,j);
		}
	}
	return S;
}

void EKF(const SXMatrix& x_t, const SXMatrix& u_t, const SXMatrix& Sigma_t, SXMatrix& x_tp1, SXMatrix& Sigma_tp1)
{
	SXMatrix A(X_DIM,X_DIM), M(X_DIM,Q_DIM), QC(U_DIM,U_DIM), MMTcar(C_DIM,C_DIM);
	SXMatrix Acar(C_DIM,C_DIM);

	QC(0,0) = config::VELOCITY_NOISE*config::VELOCITY_NOISE;
	QC(1,1) = config::TURNING_NOISE*config::TURNING_NOISE;

	//linearizeDynamics(x_t, u_t, A, M);
	linearizeCarDynamics(x_t, u_t, QC, Acar, MMTcar);

	SXMatrix SigmaCar(C_DIM,C_DIM), Sigma1(C_DIM,L_DIM), Sigma2(L_DIM,C_DIM), Sigma3(L_DIM,L_DIM);
	SigmaCar = Sigma_t(Slice(0,C_DIM),Slice(0,C_DIM));
	Sigma1 = Sigma_t(Slice(0,C_DIM),Slice(C_DIM,X_DIM));
	Sigma2 = Sigma_t(Slice(C_DIM,X_DIM),Slice(0,C_DIM));
	Sigma3 = Sigma_t(Slice(C_DIM,X_DIM),Slice(C_DIM,X_DIM));

	Sigma_tp1(CasADi::Slice(0,C_DIM),CasADi::Slice(0,C_DIM)) = SymProd<C_DIM,C_DIM>(Acar, mul(SigmaCar, trans(Acar)));//mul(Acar,mul(SigmaCar, trans(Acar)));
	Sigma_tp1(Slice(0,C_DIM),Slice(C_DIM,X_DIM)) = mul(Acar, Sigma1);
	Sigma_tp1(Slice(C_DIM,X_DIM),Slice(0,C_DIM)) = mul(Sigma2, trans(Acar));
	Sigma_tp1(Slice(C_DIM,X_DIM),Slice(C_DIM,X_DIM)) = Sigma3;

	//Sigma_tp1 = Sigma_tp1 + mul(mul(M,QC),trans(M));
	Sigma_tp1(Slice(0,C_DIM),Slice(0,C_DIM)) = Sigma_tp1(Slice(0,C_DIM),Slice(0,C_DIM)) + MMTcar;

	//Sigma_tp1 = mul(mul(A,Sigma_t),trans(A)) + mul(mul(M,QC),trans(M));

	x_tp1 = dynfunc(x_t, u_t);


	SXMatrix H(Z_DIM,X_DIM), N(Z_DIM,R_DIM), RC(Z_DIM,Z_DIM);

	linearizeObservation(x_tp1, H, N);

	SXMatrix delta = deltaMatrix(x_tp1);

	for(int i=0; i < R_DIM-1; i += 2) {
		RC(i,i) = config::OBS_DIST_NOISE*config::OBS_DIST_NOISE;
		RC(i+1,i+1) = config::OBS_ANGLE_NOISE*config::OBS_ANGLE_NOISE;
	}

	//K = ((Sigma_tp1*~H*delta)/(delta*H*Sigma_tp1*~H*delta + RC))*delta;
	//CasADi::SXMatrix K = mul(mul(mul(Sigma_tp1, mul(trans(H), delta)), solve(mul(delta, mul(H, mul(Sigma_tp1, mul(trans(H), delta)))) + RC, CasADi::SXMatrix(CasADi::DMatrix::eye(Z_DIM)))), delta);
	CasADi::SXMatrix HtransDelta = mul(trans(H), delta);
	CasADi::SXMatrix Sigma_tp1HtransDelta = mul(Sigma_tp1, HtransDelta);
	CasADi::SXMatrix K = mul(operator_divide<Z_DIM,X_DIM>(Sigma_tp1HtransDelta, SymProd<Z_DIM,X_DIM>(trans(HtransDelta), Sigma_tp1HtransDelta) + RC), delta);
	// SymProd<Z_DIM,X_DIM>(trans(HtransDelta), Sigma_tp1HtransDelta)
	// mul(trans(HtransDelta), Sigma_tp1HtransDelta)

	//Sigma_tp1 = Sigma_tp1 - mul(K,mul(H,Sigma_tp1));
	Sigma_tp1 = SymProd<X_DIM,X_DIM>(CasADi::SXMatrix(CasADi::DMatrix::eye(X_DIM)) - mul(K,H), Sigma_tp1);

}

void generateCode(FX fcn, const std::string& name){
	cout << "Generating code for " << name << endl;

	fcn.generateCode(name + ".c");
}

