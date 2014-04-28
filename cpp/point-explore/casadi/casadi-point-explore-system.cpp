#include "casadi-point-explore-system.h"
#include "../point-explore-system.h"

/**
 *
 * Constructors
 *
 */

CasadiPointExploreSystem::CasadiPointExploreSystem() {
	R = 1e-2*eye<mat>(N*R_DIM, N*R_DIM);
	this->init(ObsType::angle, CostType::entropy, R);
}

CasadiPointExploreSystem::CasadiPointExploreSystem(const ObsType obs_type, const CostType cost_type, mat& R) {
	this->init(obs_type, cost_type, R);
}

void CasadiPointExploreSystem::init(const ObsType obs_type, const CostType cost_type, mat& R) {
	this->obs_type = obs_type;
	this->cost_type = cost_type;
	this->R = R;

	if(this->cost_type == CostType::entropy) {
		this->cost_func = this->casadi_cost_entropy_func();
		this->cost_grad_func = this->casadi_cost_entropy_grad_func();
	} else {
		// TODO: platt
	}
}


/**
 * Public methods
 */

double CasadiPointExploreSystem::casadi_cost(const std::vector<mat>& X, const std::vector<mat>& U,
														const mat& P) {
	double XU_arr[T*N*X_DIM+(T-1)*N*U_DIM];
	double P_arr[T*M*X_DIM];

	this->setup_casadi_vars(X, U, P, XU_arr, P_arr);

	this->cost_func.setInput(XU_arr, 0);
	this->cost_func.setInput(P_arr, 1);

	this->cost_func.evaluate();

	double cost = 0;
	this->cost_func.getOutput(&cost, 0);

	return cost;

}
mat CasadiPointExploreSystem::casadi_cost_grad(const std::vector<mat>& X, const std::vector<mat>& U,
															const mat& P) {
	double XU_arr[T*N*X_DIM+(T-1)*N*U_DIM];
	double P_arr[T*M*X_DIM];

	this->setup_casadi_vars(X, U, P, XU_arr, P_arr);

	this->cost_func.setInput(XU_arr, 0);
	this->cost_func.setInput(P_arr, 1);

	this->cost_func.evaluate();

	mat grad(TOTAL_VARS, 1, fill::zeros);
	this->cost_func.getOutput(grad.colptr(0), 0);

	return grad;
}

/**
 * Private methods
 */

AD::SXMatrix CasadiPointExploreSystem::dist(AD::SXMatrix a, AD::SXMatrix b) {
	return sqrt(trace(mul(trans(a-b),(a-b))));
}

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
	AD::SXMatrix z(N*Z_DIM,1);

	for(int n=0; n < N; ++n) {
		AD::SXMatrix x_n = x(AD::Slice(n*X_DIM,(n+1)*X_DIM));
		AD::SXMatrix d = dist(x_n, t);
		z(n) = 1.0/(1.0+exp(-Constants::alpha*(Constants::max_range-d)));
	}

	return z;
}

AD::SXMatrix CasadiPointExploreSystem::obsfunc_angle(const AD::SXMatrix& x, const AD::SXMatrix& t) {
	AD::SXMatrix z(N*Z_DIM,1);

	for(int n=0; n < N; ++n) {
		AD::SXMatrix x_n = x(AD::Slice(n*X_DIM,(n+1)*X_DIM));
		z(n) = atan((x_n(1)-t(1))/(x_n(0)-t(0)));
	}

	return z;
}

AD::SXMatrix CasadiPointExploreSystem::gauss_likelihood(const AD::SXMatrix& v) {
	mat Sf = chol(this->R);
	mat Sf_inv = inv(Sf);

//	float Sf_diag_prod = 1;
//	for(int i=0; i < N*R_DIM; ++i) { Sf_diag_prod *= Sf(i,i); }
//	float C = pow(2*M_PI, (N*Z_DIM)/2)*Sf_diag_prod;

	AD::SXMatrix Sf_inv_casadi(N*R_DIM,N*R_DIM);
	for(int i=0; i < N*R_DIM; ++i) {
		for(int j=0; j < N*R_DIM; ++j) {
			Sf_inv_casadi(i,j) = Sf_inv(i,j);
		}
	}

	AD::SXMatrix M = mul(Sf_inv_casadi, v);

	AD::SXMatrix E_exp_sum = exp(-0.5*mul(trans(M),M));
	AD::SXMatrix w = E_exp_sum; // / C; // no need to normalize here because normalized later on anyways
	return w;
}

AD::SXMatrix CasadiPointExploreSystem::cost_entropy(const std::vector<AD::SXMatrix>& X, const std::vector<AD::SXMatrix>& U,
										const std::vector<AD::SXMatrix>& P) {
	AD::SXMatrix entropy(1,1);

	std::vector<AD::SXMatrix> X_prop(T);
	std::vector<std::vector<AD::SXMatrix> > H(T, std::vector<AD::SXMatrix>(M));
	std::vector<std::vector<std::vector<AD::SXMatrix> > > GL_H(T, std::vector<std::vector<AD::SXMatrix> >(M, std::vector<AD::SXMatrix>(M)));
	for(int t=0; t < T-1; ++t) {
		X_prop[t+1] = this->dynfunc(X[t], U[t]);
		for(int m=0; m < M; ++m) {
			H[t+1][m] = this->obsfunc(X_prop[t+1], P[m]);
		}
	}

	for(int t=1; t < T; ++t) {
		for(int m=0; m < M; ++m) {
			for(int p=m; p < M; ++p) {
				GL_H[t][m][p] = this->gauss_likelihood(H[t][m] - H[t][p]);
			}
		}
	}

	std::vector<std::vector<AD::SXMatrix> > W(T, std::vector<AD::SXMatrix>(M, AD::SXMatrix(1,1)));
	for(int m=0; m < M; ++m) { W[0][m](0,0) = 1/float(M); }
	for(int t=1; t < T; ++t) {

		AD::SXMatrix W_sum(1,1);
		for(int m=0; m < M; ++m) {
			for(int p=0; p < M; ++p) {
				int a = MIN(m,p);
				int b = MAX(m,p);
				W[t][m] += GL_H[t][a][b];
			}
			W_sum += W[t][m];
		}
		for(int m=0; m < M; ++m) { W[t][m] = W[t][m] / W_sum; }

		// use skoglar version
		AD::SXMatrix entropy_t(1,1);
		for(int m=0; m < M; ++m) {
			entropy_t += -W[t][m]*log(W[t][m]);
		}

		// simplifies because zero particle dynamics
		for(int m=0; m < M; ++m) {
			entropy_t += -W[t][m]*log(W[t-1][m]);
		}

		AD::SXMatrix sum_cross_time_weights(1,1);
		for(int m=0; m < M; ++m) {
			sum_cross_time_weights += W[t-1][m]*W[t][m];
		}
		entropy_t += log(sum_cross_time_weights);

		entropy += entropy_t;
	}

	for(int t=0; t < T-2; ++t) {
		entropy += Constants::alpha_control_smooth*(mul(trans(U[t+1]-U[t]),(U[t+1]-U[t])));
	}

	for(int t=0; t < T-1; ++t) {
		entropy += Constants::alpha_control_norm*mul(trans(U[t]),U[t]);
	}

	for(int t=1; t < T; ++t) {
		for(int n=0; n < N; ++n) {
			AD::SXMatrix x_t_n = X_prop[t](AD::Slice(n*X_DIM,(n+1)*X_DIM));
			for(int n_other=0; n_other < N; n_other++) {
				AD::SXMatrix x_t_n_other = X_prop[t](AD::Slice(n_other*X_DIM,(n_other+1)*X_DIM));
				entropy += Constants::alpha_separation*mul(trans(x_t_n-x_t_n_other), x_t_n-x_t_n_other);
			}
		}
	}

	return entropy;
}

AD::SXMatrix CasadiPointExploreSystem::cost_entropy_wrapper(const AD::SXMatrix& XU_vec, const AD::SXMatrix& P_vec) {
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

	return this->cost_entropy(X, U, P);
}


AD::SXFunction CasadiPointExploreSystem::casadi_cost_entropy_func() {
	AD::SXMatrix XU_vec = AD::ssym("XU_vec", T*N*X_DIM + (T-1)*N*U_DIM);
	AD::SXMatrix P_vec = AD::ssym("P_vec", T*M*X_DIM);

	AD::SXMatrix entropy = this->cost_entropy_wrapper(XU_vec, P_vec);

	std::vector<AD::SXMatrix> inp;
	inp.push_back(XU_vec);
	inp.push_back(P_vec);

	AD::SXFunction entropy_fcn(inp, entropy);
	entropy_fcn.init();

	return entropy_fcn;
}

AD::SXFunction CasadiPointExploreSystem::casadi_cost_entropy_grad_func() {
	AD::SXMatrix XU_vec = AD::ssym("XU_vec", T*N*X_DIM + (T-1)*N*U_DIM);
	AD::SXMatrix P_vec = AD::ssym("P_vec", T*M*X_DIM);

	AD::SXMatrix entropy = this->cost_entropy_wrapper(XU_vec, P_vec);

	AD::SXMatrix grad_entropy = gradient(entropy,XU_vec);

	// Create functions
	std::vector<AD::SXMatrix> inp;
	inp.push_back(XU_vec);
	inp.push_back(P_vec);

	std::vector<AD::SXMatrix> out;
	out.push_back(grad_entropy);
	AD::SXFunction grad_entropy_fcn(inp,out);
	grad_entropy_fcn.init();

	return grad_entropy_fcn;
}

void CasadiPointExploreSystem::setup_casadi_vars(const std::vector<mat>& X, const std::vector<mat>& U,
			const mat& P, double* XU_arr, double* P_arr) {
	int index = 0;
	for(int t=0; t < T; ++t) {
		for(int i=0; i < N*X_DIM; ++i) {
			XU_arr[index++] = X[t](i);
		}

		if (t < T-1) {
			for(int i=0; i < N*U_DIM; ++i) {
				XU_arr[index++] = U[t](i);
			}
		}
	}

	index = 0;
	for(int m=0; m < M; ++m) {
		for(int i=0; i < X_DIM; ++i) {
			P_arr[index++] = P(i,m);
		}
	}
}
