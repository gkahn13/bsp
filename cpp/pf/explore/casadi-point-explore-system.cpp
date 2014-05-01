#include "casadi-point-explore-system.h"

/**
 *
 * Constructors
 *
 */

CasadiPointExploreSystem::CasadiPointExploreSystem() {
	this->init_dims();
	mat R = 0*eye<mat>(N*R_DIM, N*R_DIM);
	this->init(ObsType::angle, CostType::entropy, R);
}

CasadiPointExploreSystem::CasadiPointExploreSystem(const ObsType obs_type, const CostType cost_type, mat& R) {
	this->init_dims();
	this->init(obs_type, cost_type, R);
}

CasadiPointExploreSystem::CasadiPointExploreSystem(const ObsType obs_type, const CostType cost_type, mat& R,
								int T, int M, int N, double DT, int X_DIM, int U_DIM, int Z_DIM, int Q_DIM, int R_DIM) {
	this->init_dims(T, M, N, DT, X_DIM, U_DIM, Z_DIM, Q_DIM, R_DIM);
	this->init(obs_type, cost_type, R);
}

void CasadiPointExploreSystem::init_dims(int T, int M, int N, double DT, int X_DIM, int U_DIM, int Z_DIM, int Q_DIM, int R_DIM) {
	this->T = T;
	this->M = M;
	this->N = N;
	this->DT = DT;
	this->X_DIM = X_DIM;
	this->U_DIM = U_DIM;
	this->Z_DIM = N*Z_DIM;
	this->Q_DIM = Q_DIM;
	this->R_DIM = R_DIM;

	this->TOTAL_VARS = T*N*X_DIM + (T-1)*N*U_DIM;
}

void CasadiPointExploreSystem::init(const ObsType obs_type, const CostType cost_type, mat& R) {
	this->obs_type = obs_type;
	this->cost_type = cost_type;
	this->R = R;

	this->cost_func = this->casadi_cost_func();
	this->cost_grad_func = this->casadi_cost_grad_func();
}


/**
 * Private methods
 */


AD::SXMatrix CasadiPointExploreSystem::dynfunc(const AD::SXMatrix& x_t, const AD::SXMatrix& u_t) {
	AD::SXMatrix x_tp1(N*X_DIM,1);

	x_tp1 = x_t + u_t*DT;

	return x_tp1;
}

AD::SXMatrix CasadiPointExploreSystem::obsfunc(const AD::SXMatrix& x, const AD::SXMatrix& t) {
	if (this->obs_type == ObsType::angle) {
		return this->obsfunc_angle(x, t);
	} else {
		return this->obsfunc_dist(x, t);
	}
}

AD::SXMatrix CasadiPointExploreSystem::obsfunc_dist(const AD::SXMatrix& x, const AD::SXMatrix& t) {
	AD::SXMatrix z(Z_DIM,1);

	for(int n=0; n < N; ++n) {
		AD::SXMatrix x_n = x(AD::Slice(n*X_DIM,(n+1)*X_DIM));
		AD::SXMatrix d = dist(x_n, t);
		z(n) = 1.0/(1.0+exp(-Constants::alpha*(Constants::max_range-d)));
	}

	return z;
}

AD::SXMatrix CasadiPointExploreSystem::obsfunc_angle(const AD::SXMatrix& x, const AD::SXMatrix& t) {
	AD::SXMatrix z(Z_DIM,1);

	for(int n=0; n < N; ++n) {
		AD::SXMatrix x_n = x(AD::Slice(n*X_DIM,(n+1)*X_DIM));
		z(n) = atan((x_n(1)-t(1))/(x_n(0)-t(0)));
	}

	return z;
}

AD::SXMatrix CasadiPointExploreSystem::cost_wrapper(const AD::SXMatrix& XU_vec, const AD::SXMatrix& P_vec) {
	std::vector<AD::SXMatrix> X(T), U(T-1), P(M);
	int index = 0;
	for(int t=0; t < T; ++t) {
		X[t] = XU_vec(AD::Slice(index,index+N*X_DIM));
		index += N*X_DIM;
		if (t < T-1) {
			U[t] = XU_vec(AD::Slice(index,index+N*U_DIM));
			index += N*U_DIM;
		}
	}

	index = 0;
	for(int m=0; m < M; ++m) {
		P[m] = P_vec(AD::Slice(index,index+X_DIM));
		index += X_DIM;
	}

	if (this->cost_type == CostType::entropy) {
		return this->cost_entropy(X, U, P);
	} else {
		return this->cost_platt(X, U, P);
	}
}
