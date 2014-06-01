#include "casadi-boxes-system.h"

/**
 *
 * Constructors/Initializers
 *
 */

CasadiBoxesSystem::CasadiBoxesSystem() {
	this->init_dims();
	mat R = 0*eye<mat>(N*R_DIM, N*R_DIM);
	mat box_dims(N*X_DIM, 1, fill::zeros);
	this->init(box_dims, ObsType::distance, CostType::entropy, R);
}

CasadiBoxesSystem::CasadiBoxesSystem(mat& box_dims, const ObsType obs_type, const CostType cost_type, mat& R) {
	this->init_dims();
	this->init(box_dims, obs_type, cost_type, R);
}

CasadiBoxesSystem::CasadiBoxesSystem(mat& box_dims, const ObsType obs_type, const CostType cost_type, mat& R,
		int T, int M, int N, double DT, int X_DIM, int U_DIM, int Z_DIM, int Q_DIM, int R_DIM) {
	this->init_dims(T, M, N, DT, X_DIM, U_DIM, Z_DIM, Q_DIM, R_DIM);
	this->init(box_dims, obs_type, cost_type, R);
}

void CasadiBoxesSystem::init_dims(int T, int M, int N, double DT, int X_DIM, int U_DIM, int Z_DIM, int Q_DIM, int R_DIM) {
	this->T = T;
	this->M = M;
	this->N = N;
	this->DT = DT;
	this->X_DIM = X_DIM;
	this->U_DIM = U_DIM;
	this->Z_DIM = N*Z_DIM;
	this->Q_DIM = Q_DIM;
	this->R_DIM = R_DIM;

	this->TOTAL_VARS = T*X_DIM + (T-1)*U_DIM;
}

void CasadiBoxesSystem::init(mat& box_dims, ObsType obs_type, CostType cost_type, mat& R) {
	this->box_dims = box_dims;
	this->obs_type = obs_type;
	this->cost_type = cost_type;
	this->R = R;

	this->cost_func = this->casadi_cost_func();
	this->cost_grad_func = this->casadi_cost_grad_func();
}


/**
 *
 * Private methods
 *
 */


AD::SXMatrix CasadiBoxesSystem::dynfunc(const AD::SXMatrix& x_t, const AD::SXMatrix& u_t) {
	AD::SXMatrix x_tp1(N*X_DIM,1);
	x_tp1 = x_t + u_t*DT;
	return x_tp1;
}


AD::SXMatrix CasadiBoxesSystem::obsfunc(const AD::SXMatrix& x, const AD::SXMatrix& b) {
	AD::SXMatrix z(Z_DIM, 1);
	AD::SXMatrix x_eye = x(0);
	AD::SXMatrix y_eye = x(1);
	for(int n=0; n < N; ++n) {
		AD::SXMatrix bx = b(n*X_DIM);
		AD::SXMatrix by = b(n*X_DIM+1);

		double width = this->box_dims(n*X_DIM);
		double height = this->box_dims(n*X_DIM+1);

		AD::SXMatrix vert_dist_to_center = fabs(y_eye-by);
		AD::SXMatrix horz_dist = fabs(x_eye - (bx+width/2.0));

		z(n) = (1.0/(1.0+exp(-Constants::alpha*(height/2.0 - vert_dist_to_center))));//*horz_dist;

//		z(n*X_DIM) = (1.0/(1.0+exp(-Constants::alpha*(height/2.0 - vert_dist_to_center))));
//		z(n*X_DIM+1) = horz_dist;
	}

	return z;
}

AD::SXMatrix CasadiBoxesSystem::cost_wrapper(const AD::SXMatrix& XU_vec, const AD::SXMatrix& P_vec) {
	std::vector<AD::SXMatrix> X(T), U(T-1), P(M);
	int index = 0;
	for(int t=0; t < T; ++t) {
		X[t] = XU_vec(AD::Slice(index,index+X_DIM));
		index += X_DIM;
		if (t < T-1) {
			U[t] = XU_vec(AD::Slice(index,index+U_DIM));
			index += U_DIM;
		}
	}

	index = 0;
	for(int m=0; m < M; ++m) {
		P[m] = P_vec(AD::Slice(index,index+N*X_DIM));
		index += N*X_DIM;
	}

	AD::SXMatrix cost;
	if (this->cost_type == CostType::entropy) {
		cost = this->cost_entropy(X, U, P);
	} else {
		cost = this->cost_platt(X, U, P);
	}

	for(int t=0; t < T-1; ++t) {
		cost += Constants::alpha_control_norm*mul(trans(U[t]), U[t]);
	}

	return cost;
}
