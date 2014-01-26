#pragma once
#ifndef __POMDP_On4_H__
#define __POMDP_On4_H__

#include "util/ilqg-matrix.h"
#include "util/utils.h"
#include <vector>

#define _sDim (((_xDim+1)*_xDim)/2)

const unsigned int COMPUTE_c = (1 << 0);
const unsigned int COMPUTE_A = (1 << 1);
const unsigned int COMPUTE_B = (1 << 2);
const unsigned int COMPUTE_M = (1 << 3);
const unsigned int COMPUTE_s = (1 << 4);
const unsigned int COMPUTE_S = (1 << 5); 
const unsigned int COMPUTE_sT = (1 << 6); 
const unsigned int COMPUTE_tT = (1 << 7);
const unsigned int COMPUTE_q = (1 << 8);
const unsigned int COMPUTE_Q = (1 << 9);
const unsigned int COMPUTE_R = (1 << 10);
const unsigned int COMPUTE_P = (1 << 11); 
const unsigned int COMPUTE_qT = (1 << 12);
const unsigned int COMPUTE_rT = (1 << 13);
const unsigned int COMPUTE_pT = (1 << 14);

//const double step = 0.0078125 / 16;
const double goldenRatio = 0.5*(3.0 - sqrt(5.0));
const double infCost = 1e10;

// Computes vector containing the elements (column by column) of the lower triangle of symmetric matrix S.
template <size_t _xDim>
inline Matrix<_sDim> vectorize(const SymmetricMatrix<_xDim>& S) {
	Matrix<_sDim> v;
	for (size_t i = 0; i < ((_xDim+1)*_xDim)/2; ++i) {
		v[i] = S[i];
	}
	return v;
}

// computes row vector such that 0.5*tr(S*Sigma) = vecTh(S)*vectorize(Sigma)
template <size_t _xDim>
inline Matrix<1,_sDim> vecTh(const SymmetricMatrix<_xDim>& S) {
	Matrix<1,_sDim> v;
	size_t idx = 0;
	for (size_t j = 0; j < _xDim; ++j) {
		for (size_t i = j; i < _xDim; ++i) {
			if (i == j) {
				v[idx] = 0.5 * S[idx];
			} else {
				v[idx] = S[idx];
			}			
			++idx;
		}
	}
	return v;
}

// Belief dynamics (extended Kalman filter)
template <size_t _xDim, size_t _zDim>
inline void beliefDynamics(const Matrix<_xDim, _xDim>& A, const SymmetricMatrix<_xDim>& M, const Matrix<_zDim, _xDim>& H, const SymmetricMatrix<_zDim>& N,
	const SymmetricMatrix<_xDim>& Sigma, SymmetricMatrix<_xDim>& SigmaNew, SymmetricMatrix<_xDim>& W) 
{
	SigmaNew = SymProd(A,Sigma*~A) + M;

	Matrix<_zDim,_xDim> HSigma = H*SigmaNew;
	W = SymProd(~HSigma,(SymProd(HSigma,~H) + N) % HSigma);
	SigmaNew -= W;
}

template <size_t _xDim, size_t _uDim, size_t _zDim>
inline void computeCEFJ(void (*linearizeDynamics)(const Matrix<_xDim>&, const Matrix<_uDim>&, Matrix<_xDim>&, Matrix<_xDim, _xDim>&, Matrix<_xDim,_uDim>&, SymmetricMatrix<_xDim>&, unsigned int), 
	void (*linearizeObservation)(const Matrix<_xDim>&, Matrix<_zDim, _xDim>&, SymmetricMatrix<_zDim>&),
	const Matrix<_xDim>& xBar, const SymmetricMatrix<_xDim>& SigmaBar, const Matrix<_uDim>& uBar,
	const Matrix<1,_sDim>& tT, const Matrix<1,_sDim>& hvecS,
	Matrix<1,_xDim>& tTC, Matrix<1,_uDim>& tTE, Matrix<1,_xDim>& vecSF, Matrix<1,_uDim>& vecSJ)
{
	Matrix<_xDim> xNext;
	Matrix<_xDim, _xDim> A;
	Matrix<_xDim, _uDim> B;
	SymmetricMatrix<_xDim> M;
	Matrix<_zDim, _xDim> H;
	SymmetricMatrix<_zDim> N;
	SymmetricMatrix<_xDim> GammaR, GammaL;
	SymmetricMatrix<_xDim> WR, WL;

	// Compute Jacobians C and F
	Matrix<_xDim> xR(xBar), xL(xBar);
	for (size_t i = 0; i < _xDim; ++i) {
		xR[i] += step; xL[i] -= step;

		linearizeDynamics(xR, uBar, xNext, A, B, M, COMPUTE_c|COMPUTE_A|COMPUTE_M);
		linearizeObservation(xNext, H, N);
		beliefDynamics(A, M, H, N, SigmaBar, GammaR, WR);

		linearizeDynamics(xL, uBar, xNext, A, B, M, COMPUTE_c|COMPUTE_A|COMPUTE_M);
		linearizeObservation(xNext, H, N);
		beliefDynamics(A, M, H, N, SigmaBar, GammaL, WL);

		tTC[i] = scalar(tT * vectorize((GammaR - GammaL) / (2.0*step)) );
		vecSF[i] = scalar(hvecS * vectorize((WR - WL) / (2.0*step)) );

		xR[i] = xL[i] = xBar[i];
	}

	// Compute E and J
	Matrix<_uDim> uR(uBar), uL(uBar);
	for (size_t i = 0; i < _uDim; ++i) {
		uR[i] += step; uL[i] -= step;

		linearizeDynamics(xBar, uR, xNext, A, B, M, COMPUTE_c|COMPUTE_A|COMPUTE_M);
		linearizeObservation(xNext, H, N);
		beliefDynamics(A, M, H, N, SigmaBar, GammaR, WR);

		linearizeDynamics(xBar, uL, xNext, A, B, M, COMPUTE_c|COMPUTE_A|COMPUTE_M);
		linearizeObservation(xNext, H, N);
		beliefDynamics(A, M, H, N, SigmaBar, GammaL, WL);

		tTE[i] = scalar(tT * vectorize((GammaR - GammaL) / (2.0*step)) );
		vecSJ[i] = scalar(hvecS * vectorize((WR - WL) / (2.0*step)) );

		uR[i] = uL[i] = uBar[i];
	}
}

// computes V dSigma_i V^T (is a symmetric matrix)
template <size_t _xDim>
inline SymmetricMatrix<_xDim> dSigma(const Matrix<_xDim, _xDim>& V, size_t row, size_t col)
{
	SymmetricMatrix<_xDim> S;
	if (row == col) {
		for (size_t j = 0; j < _xDim; ++j) {
			for (size_t i = j; i < _xDim; ++i) {
				S(i,j) = V(i,row)*V(j,row);
			}
		}
	} else {
		for (size_t j = 0; j < _xDim; ++j) {
			for (size_t i = j; i < _xDim; ++i) {
				S(i,j) = V(i,row)*V(j,col) + V(i,col)*V(j,row);
			}
		}
	}
	return S;
}

template <size_t _xDim, size_t _zDim>
inline void computeDG(const Matrix<_xDim, _xDim>& A, const SymmetricMatrix<_xDim>& M, const Matrix<_zDim, _xDim>& H, const SymmetricMatrix<_zDim>& N, 
	const SymmetricMatrix<_xDim>& SigmaBar, const Matrix<1,_sDim>& tT, const Matrix<1,_sDim>& hvecS,
	Matrix<1,_sDim>& tTD, Matrix<1,_sDim>& vecSG)
{
	Matrix<_xDim, _zDim> GammaH = (SymProd(A,SigmaBar*~A) + M)*~H;
	Matrix<_xDim, _xDim> V = A - GammaH*((SymProd(H,GammaH) + N)%(H*A));

	size_t idx = 0;
	for (size_t j = 0; j < _xDim; ++j) {
		for (size_t i = j; i < _xDim; ++i) {
			Matrix<_sDim> vecdSigmaVij = vectorize(dSigma(V, i, j));
			tTD[idx] = scalar( tT * vecdSigmaVij );
			vecSG[idx] = scalar( hvecS * (vectorize(dSigma(A, i, j)) - vecdSigmaVij) );
			++idx;
		}
	}
}


template <size_t _xDim, size_t _uDim, size_t _zDim>
inline void controlPolicy(void (*linearizeDynamics)(const Matrix<_xDim>&, const Matrix<_uDim>&, Matrix<_xDim>&, Matrix<_xDim, _xDim>&, Matrix<_xDim,_uDim>&, SymmetricMatrix<_xDim>&, unsigned int), 
	void (*linearizeObservation)(const Matrix<_xDim>&, Matrix<_zDim, _xDim>&, SymmetricMatrix<_zDim>&),
	bool (*quadratizeCost)(const Matrix<_xDim>&, const SymmetricMatrix<_xDim>&, const Matrix<_uDim>&, double&, SymmetricMatrix<_xDim>&, SymmetricMatrix<_uDim>&, Matrix<_uDim, _xDim>&, Matrix<1,_xDim>&, Matrix<1,_uDim>&, Matrix<1,_sDim>&, unsigned int),
	const Matrix<_xDim>& xBar, const SymmetricMatrix<_xDim>& SigmaBar, const Matrix<_uDim>& uBar,
	SymmetricMatrix<_xDim>& S, Matrix<1,_xDim>& sT, Matrix<1,_sDim>& tT, Matrix<_uDim,_xDim>& L, Matrix<_uDim>& l, double& gradient) 
{
	Matrix<_xDim, _xDim> A;
	Matrix<_xDim, _uDim> B;
	SymmetricMatrix<_xDim> M;
	Matrix<_xDim> xNext;

	linearizeDynamics(xBar, uBar, xNext, A, B, M, COMPUTE_c|COMPUTE_A|COMPUTE_B|COMPUTE_M); // O(n^3)

	Matrix<_zDim, _xDim> H;
	SymmetricMatrix<_zDim> N;

	linearizeObservation(xNext, H, N); // O(n^3);

	Matrix<1,_xDim> tTC;
	Matrix<1,_sDim> tTD;
	Matrix<1,_uDim> tTE;
	Matrix<1,_xDim> hvecSF;
	Matrix<1,_sDim> hvecSG;
	Matrix<1,_uDim> hvecSJ;

	Matrix<1,_sDim> hvecS = vecTh(S);

	computeCEFJ(linearizeDynamics, linearizeObservation, xBar, SigmaBar, uBar, tT, hvecS, tTC, tTE, hvecSF, hvecSJ); // O(n^4)
	computeDG(A, M, H, N, SigmaBar, tT, hvecS, tTD, hvecSG); // O(n^4)

	SymmetricMatrix<_xDim> Q;
	SymmetricMatrix<_uDim> R;
	Matrix<_uDim,_xDim> P;
	Matrix<1,_xDim> qT;
	Matrix<1,_uDim> rT;
	Matrix<1,_sDim> pT;
	double q;

	quadratizeCost(xBar, SigmaBar, uBar, q, Q, R, P, qT, rT, pT, COMPUTE_Q|COMPUTE_R|COMPUTE_P|COMPUTE_qT|COMPUTE_rT|COMPUTE_pT);

	Matrix<_xDim, _xDim> SA = S*A;

	Q += SymProd(~A,SA);
	R += SymProd(~B,S*B);
	P += ~B*SA;
	qT += sT*A + tTC + hvecSF;
	rT += sT*B + tTE + hvecSJ;
	pT += tTD + hvecSG;

	// control policy du = L dx + l
	L = -(R%P);
	l = -(R%~rT);

	// update value function
	S = Q + SymProd(~L,P);
	sT = qT + ~l*P;
	tT = pT;

	// update gradient
	gradient += scalar(rT*l);
}

template <size_t _xDim, size_t _uDim, size_t _zDim>
inline void backwardIteration(void (*linearizeDynamics)(const Matrix<_xDim>&, const Matrix<_uDim>&, Matrix<_xDim>&, Matrix<_xDim, _xDim>&, Matrix<_xDim,_uDim>&, SymmetricMatrix<_xDim>&, unsigned int), 
	void (*linearizeObservation)(const Matrix<_xDim>&, Matrix<_zDim, _xDim>&, SymmetricMatrix<_zDim>&),
	void (*quadratizeFinalCost)(const Matrix<_xDim>&, const SymmetricMatrix<_xDim>&, double&, SymmetricMatrix<_xDim>&, Matrix<1,_xDim>&, Matrix<1,_sDim>&, unsigned int),
	bool (*quadratizeCost)(const Matrix<_xDim>&, const SymmetricMatrix<_xDim>&, const Matrix<_uDim>&, double&, SymmetricMatrix<_xDim>&, SymmetricMatrix<_uDim>&, Matrix<_uDim, _xDim>&, Matrix<1,_xDim>&, Matrix<1,_uDim>&, Matrix<1,_sDim>&, unsigned int),
	const std::vector<Matrix<_xDim> >& xBar, const std::vector<SymmetricMatrix<_xDim> >& SigmaBar, const std::vector<Matrix<_uDim> >& uBar, 
	std::vector<Matrix<_uDim, _xDim> >& L, std::vector<Matrix<_uDim> >& l, double& gradient)
{
	// backward iteration: compute new control policy around nominal beliefs and control inputs
	SymmetricMatrix<_xDim> S;
	Matrix<1,_xDim> sT;
	Matrix<1,_sDim> tT;
	double s;

	quadratizeFinalCost(xBar.back(), SigmaBar.back(), s, S, sT, tT, COMPUTE_S|COMPUTE_sT|COMPUTE_tT);
	gradient = 0.0;

	for (int t = uBar.size() - 1; t != -1; --t) {
		controlPolicy(linearizeDynamics, linearizeObservation, quadratizeCost, xBar[t], SigmaBar[t], uBar[t], S, sT, tT, L[t], l[t], gradient);
	}
}

template <size_t _xDim, size_t _uDim>
inline void expectedCost(void (*linearizeDynamics)(const Matrix<_xDim>&, const Matrix<_uDim>&, Matrix<_xDim>&, Matrix<_xDim, _xDim>&, Matrix<_xDim,_uDim>&, SymmetricMatrix<_xDim>&, unsigned int), 
	bool (*quadratizeCost)(const Matrix<_xDim>&, const SymmetricMatrix<_xDim>&, const Matrix<_uDim>&, double&, SymmetricMatrix<_xDim>&, SymmetricMatrix<_uDim>&, Matrix<_uDim, _xDim>&, Matrix<1,_xDim>&, Matrix<1,_uDim>&, Matrix<1,_sDim>&, unsigned int),
	const Matrix<_xDim>& xBar, const SymmetricMatrix<_xDim>& SigmaBar, const Matrix<_uDim>& uBar, const SymmetricMatrix<_xDim>& WBar, const Matrix<_uDim,_xDim>& L,
	SymmetricMatrix<_xDim>& S, double& s) 
{
	Matrix<_xDim, _xDim> A;
	Matrix<_xDim, _uDim> B;
	SymmetricMatrix<_xDim> M;
	Matrix<_xDim> c;

	SymmetricMatrix<_xDim> Q;
	SymmetricMatrix<_uDim> R;
	Matrix<_uDim,_xDim> P;
	Matrix<1,_xDim> qT;
	Matrix<1,_uDim> rT;
	Matrix<1,_sDim> pT;
	double q;

	linearizeDynamics(xBar, uBar, c, A, B, M, COMPUTE_A|COMPUTE_B);
	bool retcode = quadratizeCost(xBar, SigmaBar, uBar, q, Q, R, P, qT, rT, pT, COMPUTE_q|COMPUTE_Q|COMPUTE_R|COMPUTE_P);

	if (retcode) {
		s = q + s; // + scalar(vecTh(S)*vectorize(WBar));
	} else {
		s = infCost;
	}

	Matrix<_xDim, _xDim> ApBL = A + B*L;
	S = Q + SymProd(~L,R*L) + SymSum(~L*P) + SymProd(~ApBL,S*ApBL);
}

template <size_t _xDim, size_t _uDim>
inline double computeExpectedCost(void (*linearizeDynamics)(const Matrix<_xDim>&, const Matrix<_uDim>&, Matrix<_xDim>&, Matrix<_xDim, _xDim>&, Matrix<_xDim,_uDim>&, SymmetricMatrix<_xDim>&, unsigned int), 
	void (*quadratizeFinalCost)(const Matrix<_xDim>&, const SymmetricMatrix<_xDim>&, double&, SymmetricMatrix<_xDim>&, Matrix<1,_xDim>&, Matrix<1,_sDim>&, unsigned int),
	bool (*quadratizeCost)(const Matrix<_xDim>&, const SymmetricMatrix<_xDim>&, const Matrix<_uDim>&, double&, SymmetricMatrix<_xDim>&, SymmetricMatrix<_uDim>&, Matrix<_uDim, _xDim>&, Matrix<1,_xDim>&, Matrix<1,_uDim>&, Matrix<1,_sDim>&, unsigned int),
	const std::vector<Matrix<_uDim, _xDim> >& L,
	const std::vector<Matrix<_xDim> >& xBar, const std::vector<SymmetricMatrix<_xDim> >& SigmaBar, const std::vector<Matrix<_uDim> >& uBar, const std::vector<SymmetricMatrix<_xDim> >& WBar) 
{
	SymmetricMatrix<_xDim> S;
	Matrix<1,_xDim> sT;
	Matrix<1,_sDim> tT;
	double s;

	size_t pathLen = uBar.size();

	// Compute expected cost
	quadratizeFinalCost(xBar[pathLen], SigmaBar[pathLen], s, S, sT, tT, COMPUTE_s|COMPUTE_S);
	for (int t = pathLen - 1; t != -1; --t) {
		expectedCost(linearizeDynamics, quadratizeCost, xBar[t], SigmaBar[t], uBar[t], WBar[t], L[t], S, s);
	}

	return s;
}

template <size_t _xDim, size_t _uDim, size_t _zDim>
inline void integrateControlPolicy(void (*linearizeDynamics)(const Matrix<_xDim>&, const Matrix<_uDim>&, Matrix<_xDim>&, Matrix<_xDim, _xDim>&, Matrix<_xDim,_uDim>&, SymmetricMatrix<_xDim>&, unsigned int), 
	void (*linearizeObservation)(const Matrix<_xDim>&, Matrix<_zDim, _xDim>&, SymmetricMatrix<_zDim>&),
	const std::vector<Matrix<_uDim, _xDim> >& L, const std::vector<Matrix<_uDim> >& l, double eps,
	const std::vector<Matrix<_xDim> >& xBar, const SymmetricMatrix<_xDim>& SigmaBar0, const std::vector<Matrix<_uDim> >& uBar, 
	std::vector<Matrix<_xDim> >& xNext, std::vector<SymmetricMatrix<_xDim> >& SigmaNext, std::vector<Matrix<_uDim> >& uNext, std::vector<SymmetricMatrix<_xDim> >& WNext)
{
	Matrix<_xDim, _xDim> A;
	Matrix<_xDim, _uDim> B;
	SymmetricMatrix<_xDim> M;
	Matrix<_zDim, _xDim> H;
	SymmetricMatrix<_zDim> N;

	size_t pathLen = uBar.size();

	int Umax = 1, Umin = -1;

	xNext[0] = xBar[0];
	SigmaNext[0] = SigmaBar0;
	for (size_t t = 0; t < pathLen; ++t) {
		uNext[t] = eps*l[t] + L[t]*(xNext[t] - xBar[t]) + uBar[t];

		for(size_t i = 0; i < _uDim; ++i) {
			uNext[t][i] = ((Umax-Umin)/2)*std::tanh(uNext[t][i]) + ((Umax+Umin)/2);
		}

		linearizeDynamics(xNext[t], uNext[t], xNext[t+1], A, B, M, COMPUTE_c|COMPUTE_A|COMPUTE_M);
		linearizeObservation(xNext[t+1], H, N);
		beliefDynamics(A, M, H, N, SigmaNext[t], SigmaNext[t+1], WNext[t]);
	}
}

template <size_t _xDim, size_t _uDim, size_t _zDim>
inline void forwardIteration(void (*linearizeDynamics)(const Matrix<_xDim>&, const Matrix<_uDim>&, Matrix<_xDim>&, Matrix<_xDim, _xDim>&, Matrix<_xDim,_uDim>&, SymmetricMatrix<_xDim>&, unsigned int), 
	void (*linearizeObservation)(const Matrix<_xDim>&, Matrix<_zDim, _xDim>&, SymmetricMatrix<_zDim>&),
	void (*quadratizeFinalCost)(const Matrix<_xDim>&, const SymmetricMatrix<_xDim>&, double&, SymmetricMatrix<_xDim>&, Matrix<1,_xDim>&, Matrix<1,_sDim>&, unsigned int),
	bool (*quadratizeCost)(const Matrix<_xDim>&, const SymmetricMatrix<_xDim>&, const Matrix<_uDim>&, double&, SymmetricMatrix<_xDim>&, SymmetricMatrix<_uDim>&, Matrix<_uDim, _xDim>&, Matrix<1,_xDim>&, Matrix<1,_uDim>&, Matrix<1,_sDim>&, unsigned int),
	const std::vector<Matrix<_uDim, _xDim> >& L, const std::vector<Matrix<_uDim> >& l, 
	std::vector<Matrix<_xDim> >& xBar, std::vector<SymmetricMatrix<_xDim> >& SigmaBar, std::vector<Matrix<_uDim> >& uBar, std::vector<SymmetricMatrix<_xDim> >& WBar, 
	double& bestCost, double& eps, double gradient)
{
	size_t pathLen = uBar.size();

	std::vector<Matrix<_xDim> > xNext(pathLen + 1);
	std::vector<SymmetricMatrix<_xDim> > SigmaNext(pathLen + 1);
	std::vector<Matrix<_uDim> > uNext(pathLen);
	std::vector<SymmetricMatrix<_xDim> > WNext(pathLen);

	double cost;

	/*for (eps *= 2.0; eps > 0.0; eps *= 0.5) {
	integrateControlPolicy(linearizeDynamics, linearizeObservation,
	L, l, eps, xBar, SigmaBar[0], uBar,
	xNext, SigmaNext, uNext, WNext);

	cost = computeExpectedCost(linearizeDynamics, quadratizeFinalCost, quadratizeCost,
	L, xNext, SigmaNext, uNext, WNext);

	if (cost < bestCost && abs(cost) < 1.0 / DBL_EPSILON) {
	bestCost = cost;
	xBar = xNext;
	SigmaBar = SigmaNext;
	uBar = uNext;
	WBar = WNext;
	return;
	}
	}
	if (eps == 0.0) {
	return;
	}*/


	// bracket minimum
	double leftEps, rightEps, bestEps, leftCost, rightCost;
	bool middleFound = false, rightFound = false;

	// Compute expected cost for eps = 0 (is always lower than current bestCost)
	bestCost = computeExpectedCost(linearizeDynamics, quadratizeFinalCost, quadratizeCost, L,
		xBar, SigmaBar, uBar, WBar);
	bestEps = 0.0;

	while (eps > 0.0 && (!middleFound || !rightFound)) {
		// Compute cost at current eps
		integrateControlPolicy(linearizeDynamics, linearizeObservation,
			L, l, eps, xBar, SigmaBar[0], uBar,
			xNext, SigmaNext, uNext, WNext);

		//std::cout << "#";

		//for (size_t i = 0; i < xBar.size(); ++i) {
		//	std::cout << "xBar: " << i << " " << xBar[i][0] << " " << xBar[i][1] << std::endl;
		//}

		cost = computeExpectedCost(linearizeDynamics, quadratizeFinalCost, quadratizeCost, L, xNext, SigmaNext, uNext, WNext);

		if (cost < bestCost && abs(cost) < 1.0 / DBL_EPSILON) {
			if (eps == 1.0) {
				bestCost = cost;
				xBar = xNext;
				SigmaBar = SigmaNext;
				uBar = uNext;
				WBar = WNext;
				return;
			}

			leftEps = bestEps;
			leftCost = bestCost;
			bestEps = eps;
			bestCost = cost;
			eps *= 2.0;
			if (eps > 1.0) {
				eps = 1.0;
			}
			middleFound = true;
		} else {
			rightEps = eps;
			rightCost = cost;
			eps *= 0.5;
			rightFound = true;
		}
	}
	if (eps == 0.0) {
		return;
	}
	//std::cout << std::endl;

	// Approximate minimum by taking parabola through three points.
	Matrix<3> costs; costs[0] = leftCost; costs[1] = bestCost; costs[2] = rightCost;
	Matrix<3> coeffs;
	Matrix<3,3> epss; 
	epss(0,0) = leftEps*leftEps; epss(0,1) = leftEps; epss(0,2) = 1.0;
	epss(1,0) = bestEps*bestEps; epss(1,1) = bestEps; epss(1,2) = 1.0;
	epss(2,0) = rightEps*rightEps; epss(2,1) = rightEps; epss(2,2) = 1.0;
	coeffs = epss % costs;
	eps = -0.5*coeffs[1]/coeffs[0];

	integrateControlPolicy(linearizeDynamics, linearizeObservation,
		L, l, eps, xBar, SigmaBar[0], uBar,
		xNext, SigmaNext, uNext, WNext);

	cost = computeExpectedCost(linearizeDynamics, quadratizeFinalCost, quadratizeCost, L, xNext, SigmaNext, uNext, WNext);

	if (cost < bestCost) {
		bestCost = cost;
	} else {
		eps = bestEps;
		integrateControlPolicy(linearizeDynamics, linearizeObservation,
			L, l, eps, xBar, SigmaBar[0], uBar,
			xNext, SigmaNext, uNext, WNext);
	}
	xBar = xNext;
	SigmaBar = SigmaNext;
	uBar = uNext;
	WBar = WNext;

	// golden section search
	/*eps = bestEps + goldenRatio*(rightEps - bestEps);

	for (size_t iter = 0; iter < 10; ++iter) {
	std::cout << "#";
	integrateControlPolicy(linearizeDynamics, linearizeObservation,
	L, l, eps, xBar, SigmaBar[0], uBar,
	xNext, SigmaNext, uNext, WNext);

	cost = computeExpectedCost(linearizeDynamics, quadratizeFinalCost, quadratizeCost, L, xNext, SigmaNext, uNext, WNext);

	if (cost < bestCost && abs(cost) < 1.0 / DBL_EPSILON) {
	if (eps > bestEps) {
	leftEps = bestEps;
	bestEps = eps;
	eps = bestEps + goldenRatio*(rightEps - bestEps);
	} else {
	rightEps = bestEps;
	bestEps = eps;
	eps = bestEps - goldenRatio*(bestEps - leftEps);
	}
	bestCost = cost;
	} else {
	if (eps > bestEps) {
	rightEps = eps;
	eps = bestEps - goldenRatio*(bestEps - leftEps);
	} else {
	leftEps = eps;
	eps = bestEps + goldenRatio*(rightEps - bestEps);
	} 
	}
	}

	eps = bestEps;
	integrateControlPolicy(linearizeDynamics, linearizeObservation,
	L, l, eps, xBar, SigmaBar[0], uBar,
	xNext, SigmaNext, uNext, WNext);

	xBar = xNext;
	SigmaBar = SigmaNext;
	uBar = uNext;
	WBar = WNext;*/
}

// nominal trajectory given as follows: xBar and SigmaBar must contain at least one item: the initial belief. uBar must contain the control inputs along the initial nominal trajectory.
template <size_t _xDim, size_t _uDim, size_t _zDim>
inline void solvePOMDP(void (*linearizeDynamics)(const Matrix<_xDim>&, const Matrix<_uDim>&, Matrix<_xDim>&, Matrix<_xDim, _xDim>&, Matrix<_xDim,_uDim>&, SymmetricMatrix<_xDim>&, unsigned int), 
	void (*linearizeObservation)(const Matrix<_xDim>&, Matrix<_zDim, _xDim>&, SymmetricMatrix<_zDim>&),
	void (*quadratizeFinalCost)(const Matrix<_xDim>&, const SymmetricMatrix<_xDim>&, double&, SymmetricMatrix<_xDim>&, Matrix<1,_xDim>&, Matrix<1,_sDim>&, unsigned int),
	bool (*quadratizeCost)(const Matrix<_xDim>&, const SymmetricMatrix<_xDim>&, const Matrix<_uDim>&, double&, SymmetricMatrix<_xDim>&, SymmetricMatrix<_uDim>&, 
						   Matrix<_uDim, _xDim>&, Matrix<1,_xDim>&, Matrix<1,_uDim>&, Matrix<1,_sDim>&, unsigned int),
	std::vector<Matrix<_xDim> >& xBar, std::vector<SymmetricMatrix<_xDim> >& SigmaBar, std::vector<Matrix<_uDim> >& uBar, std::vector<Matrix<_uDim, _xDim> >& L)
{
	double bestCost = DBL_MAX;
	bool terminate = false;

	double eps = 1.0;
	double gradient = 0.0;

	size_t pathLen = uBar.size();

	std::vector<Matrix<_uDim> > l(pathLen);
	std::vector<SymmetricMatrix<_xDim> > WBar(pathLen);
	L.resize(pathLen);
	xBar.resize(pathLen + 1);
	SigmaBar.resize(pathLen + 1);

	// integrate path
	Matrix<_xDim, _xDim> A;
	Matrix<_xDim, _uDim> B;
	SymmetricMatrix<_xDim> M;
	Matrix<_zDim, _xDim> H;
	SymmetricMatrix<_zDim> N;

	for (size_t t = 0; t < pathLen; ++t) {
		linearizeDynamics(xBar[t], uBar[t], xBar[t+1], A, B, M, COMPUTE_c|COMPUTE_A|COMPUTE_M);
		linearizeObservation(xBar[t+1], H, N);
		beliefDynamics(A, M, H, N, SigmaBar[t], SigmaBar[t+1], WBar[t]);
	}

	// compute expected cost 
	double initialCost = computeExpectedCost(linearizeDynamics, quadratizeFinalCost, quadratizeCost, L, xBar, SigmaBar, uBar, WBar);
	std::cout << "Cost to begin : " << initialCost << std::endl;
	


	size_t iter = 1;


	while(!terminate)
	{
		backwardIteration(linearizeDynamics, linearizeObservation, quadratizeFinalCost, quadratizeCost, xBar, SigmaBar, uBar, L, l, gradient);
		//double prevCost = bestCost;
		forwardIteration(linearizeDynamics, linearizeObservation, quadratizeFinalCost, quadratizeCost, L, l, xBar, SigmaBar, uBar, WBar, bestCost, eps, gradient);

		double absl = 0.0;
		double absu = 0.0;
		for (size_t i = 0; i < l.size(); ++i) {
			absl += scalar(~l[i]*l[i]);
			absu += scalar(~uBar[i]*uBar[i]);
		}
		terminate = (eps*eps*absl / absu < 1e-06);
		//terminate = (abs((prevCost - bestCost) / bestCost) < 1e-06);


		
		std::cout << "Iter: " << iter << " Grad: " << gradient << " RelErr: " << eps*eps*absl / absu << " Eps: " << eps << " Cost: " << bestCost << std::endl;

		++iter;
	}

}

#endif

