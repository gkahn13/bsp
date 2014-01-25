#ifndef __POINT_ADOLC_H__
#define __POINT_ADOLC_H__

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

#include <adolc/adolc.h>

#define TIMESTEPS 15
#define DT 1.0
#define X_DIM 2
#define U_DIM 2
#define Z_DIM 2
#define Q_DIM 2
#define R_DIM 2

#define S_DIM (((X_DIM+1)*X_DIM)/2)
#define B_DIM (X_DIM+S_DIM)


const double step = 0.0078125*0.0078125;

matrix<double> x0(X_DIM,1);
matrix<double> SqrtSigma0(X_DIM,X_DIM);
matrix<double> xGoal(X_DIM,1);
matrix<double> xMin(X_DIM,1), xMax(X_DIM,1);
matrix<double> uMin(U_DIM,1), uMax(U_DIM,1);


matrix<double> Q = identity_matrix<double>(Q_DIM);
matrix<double> R = identity_matrix<double>(R_DIM);

const int T = TIMESTEPS;
const double INFTY = 1e10;

const short int DOBS_DX_TAG = 1;
const short int DOBS_DR_TAG = 2;
const short int GRAD_COST_TAG = 3;
const short int HESS_COST_TAG = 4;

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

const double alpha_belief = 10, alpha_final_belief = 10, alpha_control = 1, alpha_goal_state = 1;

double *inputVars, *vars;
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

template<class TYPE>
void arrToMatrix(double** arr, matrix<TYPE>& M)
{
	for(int i=0; i < M.size1(); ++i)
	{
		for(int j=0; j < M.size2(); ++j)
		{
			M(i,j) = arr[i][j];
		}
	}
}

void matrixToArr(matrix<double>& M, double* arr)
{
	for(int i=0; i < M.size2(); ++i)
	{
		for(int j=0; j < M.size1(); ++j)
		{
			arr[i*M.size1()+j] = M(j,i);
		}
	}
}

// converted
template<class TYPE>
matrix<TYPE> dynfunc(const matrix<TYPE>& x, const matrix<TYPE>& u, const matrix<TYPE>& q)
{
	matrix<TYPE> xNew = x + u*DT + 0.01*q;
	return xNew;
}


// converted
// Observation model
template<class TYPE>
matrix<TYPE> obsfunc(const matrix<TYPE>& x, const matrix<TYPE>& r)
{
	TYPE intensity = sqrt(0.5*0.5*x(0,0)*x(0,0) + 1e-6);
	matrix<TYPE> z = x + intensity*r;
	return z;
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

void initAdolcMatrix(matrix<adouble>& mAdolc, matrix<double>& m)
{
	for(int i=0; i < mAdolc.size2(); ++i) {
		for(int j=0; j < mAdolc.size1(); ++j) {
			mAdolc(j,i) <<= m(j,i);
		}
	}
}

void retrieveAdolcMatrix(matrix<double>& m, matrix<adouble>& mAdolc)
{
	for(int i=0; i < m.size2(); ++i) {
		for(int j=0; j < m.size1(); ++j) {
			mAdolc(j,i) >>= m(j,i);
		}
	}
}

template<class TYPE>
void finiteDiffJac(matrix<TYPE>& x, matrix<TYPE>& r, matrix<TYPE>& H, matrix<TYPE>& N)
{
	matrix<TYPE> xr(x), xl(x);
	for (size_t i = 0; i < X_DIM; ++i) {
		xr(i,0) += step; xl(i,0) -= step;
		column(H, i) = column((obsfunc(xr, r) - obsfunc(xl, r)) / (xr(i,0) - xl(i,0)),0);
		xr(i,0) = x(i,0); xl(i,0) = x(i,0);
	}
	//printMatrix<double>("H finite diff", H);

	matrix<TYPE> rr(x), rl(x);
	for (size_t i = 0; i < R_DIM; ++i) {
		rr(i,0) += step; rl(i,0) -= step;
		column(N, i) = column((obsfunc(x, rr) - obsfunc(x, rl)) / (rr(i,0) - rl(i,0)),0);
		rr(i,0) = r(i,0); rl(i,0) = r(i,0);
	}
	//printMatrix<double>("N finite diff", N);

}

// converted (assuming no work needs to be done)
// Jacobians: dh(x,r)/dx, dh(x,r)/dr
template<class TYPE>
void linearizeObservation(matrix<TYPE>& x, matrix<TYPE>& r, matrix<TYPE>& H, matrix<TYPE>& N)
{
	LOG_DEBUG("inside linearizeObservation adolc");
	matrix<adouble> xAdolc(X_DIM,1);
	matrix<double> xObs(Z_DIM,1);

	trace_on(DOBS_DX_TAG);
	initAdolcMatrix(xAdolc, x);
	matrix<adouble> xObsAdolc = obsfunc<adouble>(xAdolc, r);
	retrieveAdolcMatrix(xObs, xObsAdolc);
	trace_off(DOBS_DX_TAG);

	double* x_arr = new double[X_DIM];
	std::copy(x.begin1(), x.end1(), x_arr);
	double** dobs_dx = new double*[Z_DIM];
	for(int i=0; i < Z_DIM; ++i) { dobs_dx[i] = new double[X_DIM]; }

	jacobian(DOBS_DX_TAG, Z_DIM, X_DIM, x_arr, dobs_dx);
	arrToMatrix(dobs_dx, H);

	matrix<adouble> rAdolc(R_DIM,1);
	matrix<double> rObs(Z_DIM,1);

	trace_on(DOBS_DR_TAG);
	initAdolcMatrix(rAdolc, r);
	matrix<adouble> rObsAdolc = obsfunc<adouble>(x, rAdolc);
	retrieveAdolcMatrix(rObs, rObsAdolc);
	trace_off(DOBS_DR_TAG);

	double* r_arr = new double[R_DIM];
	std::copy(r.begin1(), r.end1(), r_arr);
	double** dobs_dr = new double*[Z_DIM];
	for(int i=0; i < Z_DIM; ++i) { dobs_dr[i] = new double[R_DIM]; }

	jacobian(DOBS_DR_TAG, Z_DIM, R_DIM, r_arr, dobs_dr);
	arrToMatrix(dobs_dr, N);

	//finiteDiffJac(x,r);
}

/*
// adouble version b/c have to deal with embedded traces as a special case (I think...)
// Jacobians: dh(x,r)/dx, dh(x,r)/dr
template <>
void linearizeObservation<adouble>(matrix<adouble>& x, matrix<adouble>& r, matrix<adouble>& H, matrix<adouble>& N)
{
	//matrix<double> x_double(X_DIM,1), r_double(R_DIM,1);
	// extract value from x and r to compute jacs
	//for(int i=0; i < X_DIM; ++i) { x_double(i,0) = x(i,0).value(); }
	//for(int i=0; i < R_DIM; ++i) { r_double(i,0) = r(i,0).value(); }

	LOG_DEBUG("inside linearizeObservation adolc");
	matrix<adouble> xAdolc(X_DIM,1);
	matrix<double> xObs(Z_DIM,1);

	trace_on(DOBS_DX_TAG);
	//initAdolcMatrix(xAdolc, x_double);
	matrix<adouble> xObsAdolc = obsfunc<adouble>(x, r);
	//retrieveAdolcMatrix(xObs, xObsAdolc);
	trace_off(DOBS_DX_TAG);

	double* x_arr = new double[X_DIM];
	std::copy(x.begin1(), x.end1(), x_arr);
	double** dobs_dx = new double*[Z_DIM];
	for(int i=0; i < Z_DIM; ++i) { dobs_dx[i] = new double[X_DIM]; }

	jacobian(DOBS_DX_TAG, Z_DIM, X_DIM, x_arr, dobs_dx);
	arrToMatrix(dobs_dx, H);


	matrix<adouble> rAdolc(R_DIM,1);
	matrix<double> rObs(Z_DIM,1);

	trace_on(DOBS_DR_TAG);
	initAdolcMatrix(rAdolc, r_double);
	matrix<adouble> rObsAdolc = obsfunc<adouble>(x_double, rAdolc);
	retrieveAdolcMatrix(rObs, rObsAdolc);
	trace_off(DOBS_DR_TAG);

	double* r_arr = new double[R_DIM];
	std::copy(r_double.begin1(), r_double.end1(), r_arr);
	double** dobs_dr = new double*[Z_DIM];
	for(int i=0; i < Z_DIM; ++i) { dobs_dr[i] = new double[R_DIM]; }

	jacobian(DOBS_DR_TAG, Z_DIM, R_DIM, r_arr, dobs_dr);
	arrToMatrix(dobs_dr, N);


	//finiteDiffJac(x,r);
}
*/



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
	matrix<double> M = .01*identityMatrix<double>(U_DIM); // d(dynfunc)/dm (noise)


	printMatrix("A",A);
	printMatrix("M",M);

	matrix<TYPE> Sigma_tp1 = prod(((matrix<TYPE>)prod(A, Sigma_t)), trans(A)) + prod(((matrix<TYPE>)prod(M, Q)), trans(M));

	printMatrix("Sigma_tp1",Sigma_tp1);

	LOG_DEBUG("before dynfunc");
	matrix<TYPE> x_tp1 = dynfunc<TYPE>(x_t, u_t, zeroMatrix<double>(Q_DIM,1));

	printMatrix("x_tp1",x_tp1);

	//matrix<TYPE> H = zeroMatrix<TYPE>(Z_DIM,X_DIM);
	matrix<double> H = zeroMatrix<double>(Z_DIM,X_DIM);
	matrix<TYPE> N = zeroMatrix<TYPE>(Z_DIM,R_DIM);
	matrix<TYPE> r = zeroMatrix<TYPE>(R_DIM,1); // correct?
	LOG_DEBUG("before linearize observation");

	H(0,0) = 1; H(1,1) = 1;
	TYPE intensity = sqrt(x_tp1(0,0) * x_tp1(0,0) * 0.5 * 0.5 + 1e-6);
	N(0,0) = intensity;
	N(1,1) = intensity;
	//linearizeObservation<TYPE>(x_tp1, r, H, N);

	printMatrix("H",H);
	printMatrix("N",N);


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

/*
void pythonDisplayTrajectory(std::vector< Matrix<B_DIM> >& B, std::vector< Matrix<U_DIM> >& U)
{
	for (int t = 0; t < T-1; ++t) {
		B[t+1] = beliefDynamics(B[t], U[t]);
	}

	py::list Bvec;
	for(int j=0; j < B_DIM; j++) {
		for(int i=0; i < T; i++) {
			Bvec.append(B[i][j]);
		}
	}

	py::list Uvec;
		for(int j=0; j < U_DIM; j++) {
			for(int i=0; i < T-1; i++) {
			Uvec.append(U[i][j]);
		}
	}

	py::list x0_list, xGoal_list;
	for(int i=0; i < X_DIM; i++) {
		x0_list.append(x0[i]);
		xGoal_list.append(xGoal[i]);
	}

	std::string workingDir = boost::filesystem::current_path().normalize().string();
	std::string bspDir = workingDir.substr(0,workingDir.find("bsp"));

	try
	{
		Py_Initialize();
		py::object main_module = py::import("__main__");
		py::object main_namespace = main_module.attr("__dict__");
		py::exec("import sys, os", main_namespace);
		py::exec(py::str("sys.path.append('"+bspDir+"bsp/python')"), main_namespace);
		py::exec("from bsp_light_dark import LightDarkModel", main_namespace);
		py::object model = py::eval("LightDarkModel()", main_namespace);
		py::object plot_mod = py::import("plot");
		py::object plot_traj = plot_mod.attr("plot_belief_trajectory_cpp");

		plot_traj(Bvec, Uvec, model, x0_list, xGoal_list, T);
	}
	catch(py::error_already_set const &)
	{
		PyErr_Print();
	}

}
*/
#endif
