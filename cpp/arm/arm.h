#ifndef __POINT_H__
#define __POINT_H__

#include <fstream>
#include <math.h>

#include "util/matrix.h"
#include "util/utils.h"
//#include "util/Timer.h"
#include "util/logging.h"
//#include "util/invertMatrix.hpp"

#include <Python.h>
//#include <pythonrun.h>
#include <boost/python.hpp>
#include <boost/filesystem.hpp>

namespace py = boost::python;

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
using namespace boost::numeric::ublas;

#include <adolc/adouble.h>

#define TIMESTEPS 5
#define DT 1
#define X_DIM 6
#define U_DIM 6
#define Z_DIM 4
#define Q_DIM 6
#define R_DIM 4
#define G_DIM 3 // dimension of g(x) function, which is 3d position

#define S_DIM (((X_DIM+1)*X_DIM)/2)
#define B_DIM (X_DIM+S_DIM)

const double l4 = 2.375;
const double l3 = 10.375;
const double l2 = 8;
const double l1 = 7.25;


const double step = .1;//0.0078125*0.0078125;

matrix<double> x0(X_DIM,1);
matrix<double> SqrtSigma0(X_DIM,X_DIM);
//matrix<adouble> posGoal;
matrix<double> xGoal(X_DIM,1); // TODO: temporary, since goal is a vector of joints
matrix<double> xMin(X_DIM,1), xMax(X_DIM,1);
matrix<double> uMin(U_DIM,1), uMax(U_DIM,1);

matrix<double> cam0(G_DIM,1), cam1(G_DIM,1);

matrix<double> Q = identity_matrix<double>(Q_DIM);
matrix<double> R = identity_matrix<double>(R_DIM);

const int T = TIMESTEPS;
const double INFTY = 1e10;

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

const double alpha_belief = 10, alpha_final_belief = 10, alpha_control = 1, alpha_goal_state = 1;

adouble *inputVars, *vars;
std::vector<int> maskIndices;

template<class TYPE>
matrix<TYPE> zeroMatrix(int rows, int cols)
{
	matrix<TYPE> M(rows,cols);
	for(int i=0; i < M.size1(); ++i)
	{
		for(int j=0; j < M.size2(); ++j)
		{
			M(i,j) = 0.0;
		}
	}
	return M;
}

template<class TYPE>
matrix<TYPE> identityMatrix(int size)
{
	matrix<TYPE> M = zeroMatrix<TYPE>(size, size);
	for(int i=0; i < size; ++i) { M(i,i) = 1.0; }
	return M;
}

template<class TYPE>
matrix<TYPE> copyMatrix(matrix<TYPE> M)
{
	matrix<TYPE> Mcopy(M.size1(),M.size2());
	for(int i=0; i < M.size1(); ++i)
	{
		for(int j=0; j < M.size2(); ++j)
		{
			Mcopy(i,j) = M(i,j);
		}
	}
	return M;
}

template<class TYPE>
void printMatrix(std::string name, matrix<TYPE> M)
{
	matrix<TYPE> transM = trans(M);
	std::cout << name << std::endl;
	for(int i=0; i < transM.size2(); ++i)
	{
		std::cout << column(transM,i) << std::endl;
	}
}

// converted
template<class TYPE>
matrix<TYPE> dynfunc(const matrix<TYPE>& x, const matrix<TYPE>& u, const matrix<TYPE>& q)
{
	matrix<TYPE> result = x + (u + q)*DT;
	return result;
}


// converted
// joint angles -> end effector position
template<class TYPE>
matrix<TYPE> g(const matrix<TYPE>& x)
{
    TYPE a0 = x(0,0), a1 = x(1,0), a2 = x(2,0), a3 = x(3,0), a4 = x(4,0), a5 = x(5,0);

    matrix<TYPE> p(G_DIM, 1);
    p(0,0) = sin(a0) * (cos(a1) * (sin(a2) * (cos(a4) * l4 + l3) + cos(a2) * cos(a3) * sin(a4) * l4) + sin(a1) * (cos(a2) * (cos(a4) * l4 + l3) - sin(a2) * cos(a3) * sin(a4) * l4 + l2)) + cos(a0) * sin(a3) * sin(a4) * l4;
    p(1,0) = -sin(a1) * (sin(a2) * (cos(a4) * l4 + l3) + cos(a2) * cos(a3) * sin(a4) * l4) + cos(a1) * (cos(a2) * (cos(a4) * l4 + l3) - sin(a2) * cos(a3) * sin(a4) * l4 + l2) + l1;
    p(2,0) = cos(a0) * (cos(a1) * (sin(a2) * (cos(a4) * l4 + l3) + cos(a2) * cos(a3) * sin(a4) * l4) + sin(a1) * (cos(a2) * (cos(a4) * l4 + l3) - sin(a2) * cos(a3) * sin(a4) * l4 + l2)) - sin(a0) * sin(a3) * sin(a4) * l4;

    return p;
}

// converted
// Observation model
template<class TYPE>
matrix<TYPE> obsfunc(const matrix<TYPE>& x, const matrix<TYPE>& r)
{
	matrix<TYPE> ee_pos = g<TYPE>(x);

	matrix<TYPE> obs(Z_DIM, 1);
	obs(0,0) = (ee_pos(0,0) - cam0(0,0)) / (ee_pos(1,0) - cam0(1,0)) + r(0,0);
	obs(1,0) = (ee_pos(2,0) - cam0(2,0)) / (ee_pos(1,0) - cam0(1,0)) + r(1,0);
	obs(2,0) = (ee_pos(0,0) - cam1(0,0)) / (ee_pos(1,0) - cam1(1,0)) + r(2,0);
    obs(3,0) = (ee_pos(2,0) - cam1(2,0)) / (ee_pos(1,0) - cam1(1,0)) + r(3,0);

    return obs;
}

// converted (assuming no work needs to be done)
// Jacobians: df(x,u,q)/dx, df(x,u,q)/dq
// TODO: remove!
void linearizeDynamics(const matrix<adouble>& x, const matrix<adouble>& u, const matrix<adouble>& q, matrix<adouble>& A, matrix<adouble>& M)
{
	//matrix<double> m(3,3), v(3,1);
	//column(m,0) = column(v,0);

	//A.reset();
	matrix<adouble> xr(x), xl(x);
	for (size_t i = 0; i < X_DIM; ++i) {
		xr(i,0) += step; xl(i,0) -= step;
		column(A, i) = column((dynfunc(xr, u, q) - dynfunc(xl, u, q)) / (xr(i,0) - xl(i,0)),0);
		xr(i,0) = x(i,0); xl(i,0) = x(i,0);
	}


	//M.reset();
	matrix<adouble> qr(q), ql(q);
	for (size_t i = 0; i < Q_DIM; ++i) {
		qr(i,0) += step; ql(i,0) -= step;
		column(M, i) = column((dynfunc(x, u, qr) - dynfunc(x, u, ql)) / (qr(i,0) - ql(i,0)),0);
		qr(i,0) = q(i,0); ql(i,0) = q(i,0);

	}

}

// converted (assuming no work needs to be done)
// Jacobians: dh(x,r)/dx, dh(x,r)/dr
// TODO: remove!
void linearizeObservation(const matrix<adouble>& x, const matrix<adouble>& r, matrix<adouble>& H, matrix<adouble>& N)
{
	LOG_DEBUG("linearize observation");
	/*
	//H.reset();
	printMatrix("x",x);
	matrix<adouble> xr(x), xl(x);
	for (size_t i = 0; i < X_DIM; ++i) {
		xr(i,0) += step; xl(i,0) -= step;
		column(H, i) = column((obsfunc(xr, r) - obsfunc(xl, r)) / (xr(i,0) - xl(i,0)),0);
		xr(i,0) = x(i,0); xl(i,0) = x(i,0);
	}
	*/
	//printMatrix("H",H);

	//for(int i=0; i < R_DIM; ++i) { r(i,0) = 0; } // TODO: remove
	printMatrix("r",r);
	printMatrix("obsfunc(x,r)",obsfunc(x,r));
	//N.reset();
	matrix<adouble> rr(r), rl(r);
	for (size_t i = 0; i < R_DIM; ++i) {
		rr(i,0) += step; rl(i,0) -= step;
		//printMatrix("rr",rr); printMatrix("rl",rl);
		LOG_DEBUG("%d",i);
		printMatrix("obsfunc(x,rr)",obsfunc(x,rr));
		printMatrix("obsfunc(x,rl)",obsfunc(x,rl));
		std::cout << "denom " << rr(i,0)-rl(i,0) << std::endl;
		//printMatrix("whole thing", (obsfunc(x, rr) - obsfunc(x, rl)) / (2*step));
		column(N, i) = column((obsfunc(x, rr) - obsfunc(x, rl)) / (rr(i,0) - rl(i,0)),0);
		rr(i,0) = r(i,0); rl(i,0) = r(i,0);
	}

	printMatrix("N",N);
	exit(0);

}


// converted (assuming no work needs to be done)
// Switch between belief vector and matrices
template<class TYPE>
void unVec(const matrix<TYPE>& b, matrix<TYPE>& x, matrix<TYPE>& SqrtSigma) {
	matrix<TYPE> x_tmp(X_DIM,1);
	for(int i=0; i < X_DIM; ++i) { x_tmp(i,0) = b(i,0); };
	x = x_tmp;

	size_t idx = X_DIM;
	for (size_t j = 0; j < X_DIM; ++j) {
		for (size_t i = j; i < X_DIM; ++i) {
			SqrtSigma(i,j) = b(idx,0);
			SqrtSigma(j,i) = b(idx,0);
			++idx;
		}
	}
}

// converted (assuming no work needs to be done)
template<class TYPE>
void vec(const matrix<TYPE>& x, const matrix<TYPE>& SqrtSigma, matrix<TYPE>& b) {
	for(int i=0; i < X_DIM; ++i) { b(i,0) = x(i,0); }
	size_t idx = X_DIM;
	for (size_t j = 0; j < X_DIM; ++j) {
		for (size_t i = j; i < X_DIM; ++i) {
			b(idx,0) = 0.5 * (SqrtSigma(i,j) + SqrtSigma(j,i));
			++idx;
		}
	}
}

template<class TYPE>
matrix<TYPE> mdivide(matrix<TYPE> A, matrix<TYPE> B)
{
	// Cholesky factorization A = L*~L
	// check if symmetric
	int size = A.size1();
	matrix<TYPE> L = zeroMatrix<TYPE>(size,size);
	for (int i = 0; i < size; ++i)
	{
		for (int j = i; j < size; ++j)
		{
			TYPE sum = A(j, i);
			for (int k = 0; k < i; ++k)
			{
				sum -= L(j, k) * L(i, k);
			}
			if (i == j)
			{
				L(i, i) = sqrt(sum);
			}
			else
			{
				L(j, i) = sum / L(i, i);
			}
		}
	}

	int ncols = B.size2();
	// Backward and forward substitution
	matrix<TYPE> M(size, ncols);
	for (int i = 0; i < size; ++i)
	{
		for (int k = 0; k < ncols; ++k)
		{
			TYPE sum = B(i, k);
			for (int j = 0; j < i; ++j)
			{
				sum -= L(i, j) * M(j, k);
			}
			M(i, k) = sum / L(i, i);
		}
	}
	for (int i = size - 1; i != -1; --i)
	{
		for (int k = 0; k < ncols; ++k)
		{
			TYPE sum = M(i, k);
			for (int j = i + 1; j < size; ++j)
			{
				sum -= L(j, i) * M(j, k);
			}
			M(i, k) = sum / L(i, i);
		}
	}
	return M;
}

/*
// TODO: can't do sqrt(Sigma), so can't comlete function
// converted (same except for A, M)
// Belief dynamics
matrix<adouble> beliefDynamics(const matrix<adouble>& b, const matrix<adouble>& u) {
	matrix<adouble> x(X_DIM,1);
	matrix<adouble> SqrtSigma(X_DIM,X_DIM);
	unVec(b, x, SqrtSigma);

	matrix<adouble> Sigma(X_DIM,X_DIM);
	Sigma = prod(SqrtSigma, SqrtSigma);


	//linearizeDynamics(x, u, zeros<Q_DIM,1>(), A, M);
	matrix<adouble> A = identityMatrix(X_DIM);
	matrix<adouble> M = DT*identityMatrix(U_DIM);

	x = dynfunc(x, u, zeroMatrix(Q_DIM,1));
	Sigma = prod(A, matrix<adouble>(prod(Sigma,trans(A)))) + prod(M,trans(M));

	matrix<adouble> H(Z_DIM,X_DIM);
	matrix<adouble> N(Z_DIM,R_DIM);
	matrix<adouble> r(R_DIM,1);
	linearizeObservation(x, r, H, N);

	matrix<adouble> prevInv = prod(H, matrix<adouble>(prod(Sigma,trans(H))));
	prevInv = prevInv + prod(N,trans(N));

	matrix<adouble> prevK = prod(Sigma, trans(H));
	matrix<adouble> K = mdivide(prevK, prevInv);
	//matrix<adouble> K = Sigma*~H/(H*Sigma*~H + N*~N);

	matrix<adouble> prevSigma = identityMatrix(X_DIM) - prod(K,H);
	Sigma = prod(prevSigma, Sigma);
	//Sigma = (identity<X_DIM>() - K*H)*Sigma;

	//matrix<adouble> b_g(B_DIM,1);
	//vec(x, sqrt(Sigma), b_g);

	return b;


	//return b_g;
}
*/

template<class TYPE>
matrix<TYPE> EKF(const matrix<TYPE>& x_t, const matrix<TYPE>& u_t, const matrix<TYPE>& q_t, const matrix<TYPE>& r_t, const matrix<TYPE>& Sigma_t)
{
	printMatrix("x_t",x_t);
	printMatrix("u_t",u_t);
	//printMatrix("q_t",q_t);
	//printMatrix("r_t",r_t);
	printMatrix("Sigma_t",Sigma_t);

	//computeDynJacobians(x_t, u_t, q_t, out A, out M);
	matrix<double> A = identityMatrix<double>(X_DIM); // d(dynfunc)/dx
	matrix<double> M = DT*identityMatrix<double>(U_DIM); // d(dynfunc)/dm (noise)


	printMatrix("A",A);
	printMatrix("M",M);

	matrix<TYPE> Sigma_tp1 = prod(((matrix<TYPE>)prod(A, Sigma_t)), trans(A)) + prod(((matrix<TYPE>)prod(M, Q)), trans(M));

	printMatrix("Sigma_tp1",Sigma_tp1);

	LOG_DEBUG("before dynfunc");
	matrix<TYPE> x_tp1 = dynfunc<TYPE>(x_t, u_t, zeroMatrix<double>(Q_DIM,1));

	printMatrix("x_tp1",x_tp1);

	matrix<TYPE> H(Z_DIM,X_DIM), N(Z_DIM,R_DIM);
	//matrix<adouble> r = zeroMatrix(R_DIM,1); // correct?
	LOG_DEBUG("before linearize observation");
	//linearizeObservation(x_tp1, zeroMatrix<double>(R_DIM,1), H, N); // TODO: use adol-c

	printMatrix("H",H);
	printMatrix("N",N);

	exit(0);

	LOG_DEBUG("before K1, K2");
	//printMatrix("Sigma_tp1",Sigma_tp1);
	//printMatrix("trans(H)",trans(H));
	matrix<TYPE> K1 = prod(Sigma_tp1, trans(H));
	matrix<TYPE> K2 = prod(((matrix<TYPE>)prod(H, Sigma_tp1)), trans(H)) + prod(((matrix<TYPE>)prod(N, R)), trans(N));

	LOG_DEBUG("before K");
	matrix<TYPE> K = trans(mdivide<TYPE>(trans(K2), trans(K1)));

	Sigma_tp1 = prod(((matrix<TYPE>)(identityMatrix<double>(X_DIM) - prod(K,H))), Sigma_tp1);

	LOG_DEBUG("before b_tp1");
	matrix<TYPE> b_tp1(B_DIM,1);
	vec<TYPE>(x_tp1, Sigma_tp1, b_tp1);
	return b_tp1;
}

template<class TYPE>
TYPE costfunc(const matrix<TYPE>& X, const matrix<TYPE>& U, const matrix<TYPE>& q, const matrix<TYPE>& r, const matrix<TYPE>& Sigma_0)
{
	LOG_DEBUG("costfunc");
	TYPE cost = 0;

	matrix<TYPE> x_tp1(X_DIM,1);
	column(x_tp1,0) = column(X,0);
	matrix<TYPE> Sigma_t(X_DIM,X_DIM), Sigma_tp1(X_DIM,X_DIM);
	Sigma_t = Sigma_0;

	matrix<TYPE> Ucol(U_DIM,1), Xcol(X_DIM,1);
	for (int t = 0; t < T - 1; ++t)
	{
		LOG_DEBUG("EKF iteration %d",t);
		for(int i=0; i < X_DIM; ++i) {
			cost = cost + alpha_belief * Sigma_t(i,i);
		}
		column(Ucol,0) = column(U,t);
		cost = cost + alpha_control * (prod(trans(Ucol), Ucol))(0,0);

		column(Xcol,0) = column(X,t);
		LOG_DEBUG("before EKF");
		matrix<TYPE> b_tp1 = EKF<TYPE>(Xcol, Ucol, q, r, Sigma_t);
		LOG_DEBUG("after EKF");
		unVec<TYPE>(b_tp1, x_tp1, Sigma_t);
	}

	for(int i=0; i < X_DIM; ++i) {
		cost = cost + alpha_final_belief * Sigma_t(i,i);
	}

	return cost;
}


/*
// converted (assuming no work needs to be done)
void setupDstarInterface(std::string mask) {
	std::stringstream ss(mask);
	int val, i=0;
	while (ss >> val) {
		if (val == 1) {
			maskIndices.push_back(i);
		}
		i++;
	}

	inputVars = new double[i];
}

// converted
void pythonDisplayTrajectory(std::vector< matrix<adouble> >& B, std::vector< matrix<adouble> >& U)
{
	// no python for this example
	// don't use
}
*/
#endif
