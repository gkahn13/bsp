#include "casadi-system.h"

/**
 * Constructor
 */

CasadiSystem::CasadiSystem() {

}

/**
 * Public methods
 */

double CasadiSystem::casadi_cost(const std::vector<mat>& X, const std::vector<mat>& U,
														const mat& P) {
	double XU_arr[TOTAL_VARS];
	double P_arr[T*M*X_DIM];

	this->setup_casadi_vars(X, U, P, XU_arr, P_arr);

	this->cost_func.setInput(XU_arr, 0);
	this->cost_func.setInput(P_arr, 1);

	this->cost_func.evaluate();

	double cost = 0;
	this->cost_func.getOutput(&cost, 0);

	return cost;

}

mat CasadiSystem::casadi_cost_grad(const std::vector<mat>& X, const std::vector<mat>& U,
															const mat& P) {
	double XU_arr[TOTAL_VARS];
	double P_arr[T*M*X_DIM];

	this->setup_casadi_vars(X, U, P, XU_arr, P_arr);

	this->cost_grad_func.setInput(XU_arr, 0);
	this->cost_grad_func.setInput(P_arr, 1);

	this->cost_grad_func.evaluate();

	mat grad(TOTAL_VARS, 1, fill::zeros);
	this->cost_grad_func.getOutput(grad.colptr(0), 0);

	return grad;
}


/**
 * Private methods
 */

AD::SXMatrix CasadiSystem::dist(AD::SXMatrix a, AD::SXMatrix b) {
	return sqrt(trace(mul(trans(a-b),(a-b))));
}

AD::SXMatrix CasadiSystem::gauss_likelihood(const AD::SXMatrix& v) {
	mat Sf = chol(this->R);
	mat Sf_inv = inv(Sf);

	int r_size = this->R.n_cols;

	float C = pow(2*M_PI, r_size/2)*prod(diagvec(Sf));

	AD::SXMatrix Sf_inv_casadi(r_size,r_size);
	for(int i=0; i < r_size; ++i) {
		for(int j=0; j < r_size; ++j) {
			Sf_inv_casadi(i,j) = Sf_inv(i,j);
		}
	}

	AD::SXMatrix M = mul(Sf_inv_casadi, v);

	AD::SXMatrix E_exp_sum = exp(-0.5*mul(trans(M),M));
	AD::SXMatrix w = E_exp_sum / C; // no need to normalize here because normalized later on anyways

	return w;
}

AD::SXMatrix CasadiSystem::cost_entropy(const std::vector<AD::SXMatrix>& X, const std::vector<AD::SXMatrix>& U,
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

//	for(int t=0; t < T-2; ++t) {
//		entropy += Constants::alpha_control_smooth*(mul(trans(U[t+1]-U[t]),(U[t+1]-U[t])));
//	}
//
//	for(int t=0; t < T-1; ++t) {
//		entropy += Constants::alpha_control_norm*mul(trans(U[t]),U[t]);
//	}
//
//	for(int t=1; t < T; ++t) {
//		for(int n=0; n < N; ++n) {
//			AD::SXMatrix x_t_n = X_prop[t](AD::Slice(n*X_DIM,(n+1)*X_DIM));
//			for(int n_other=0; n_other < N; n_other++) {
//				AD::SXMatrix x_t_n_other = X_prop[t](AD::Slice(n_other*X_DIM,(n_other+1)*X_DIM));
//				entropy += Constants::alpha_separation*mul(trans(x_t_n-x_t_n_other), x_t_n-x_t_n_other);
//			}
//		}
//	}

	return entropy;
}


AD::SXMatrix CasadiSystem::cost_platt(const std::vector<AD::SXMatrix>& X, const std::vector<AD::SXMatrix>& U,
		const std::vector<AD::SXMatrix>& P) {
	std::vector<AD::SXMatrix> X_prop(T);
	X_prop[0] = X[0];
	for(int t=0; t < T-1; ++t) {
		X_prop[t+1] = this->dynfunc(X[t], U[t]);
	}

	std::vector<AD::SXMatrix> H(M, AD::SXMatrix(T*(Z_DIM),1));
	for(int m=0; m < M; ++m) {
		int index = 0;
		for(int t=0; t < T; ++t) {
			AD::SXMatrix Hm = this->obsfunc(X_prop[t], P[m]);
			for(int i=0; i < Z_DIM; ++i) {
				H[m](index++,0) = Hm(i,0);
			}
		}
	}

	AD::SXMatrix platt(1,1);
	for(int m=1; m < M; ++m) {
		AD::SXMatrix diff = H[m] - H[0];
		platt(0,0) += (1/float(M-1))*exp(-mul(trans(diff), diff));
	}

	return platt;
}

AD::SXFunction CasadiSystem::casadi_cost_func() {
	AD::SXMatrix XU_vec = AD::ssym("XU_vec", TOTAL_VARS);
	AD::SXMatrix P_vec = AD::ssym("P_vec", T*M*X_DIM);

	AD::SXMatrix cost = this->cost_wrapper(XU_vec, P_vec);

	std::vector<AD::SXMatrix> inp;
	inp.push_back(XU_vec);
	inp.push_back(P_vec);

	AD::SXFunction cost_fcn(inp, cost);
	cost_fcn.init();

	return cost_fcn;
}

AD::SXFunction CasadiSystem::casadi_cost_grad_func() {
	AD::SXMatrix XU_vec = AD::ssym("XU_vec", TOTAL_VARS);
	AD::SXMatrix P_vec = AD::ssym("P_vec", T*M*X_DIM);

	AD::SXMatrix cost = this->cost_wrapper(XU_vec, P_vec);

	AD::SXMatrix grad_cost = gradient(cost, XU_vec);

	// Create functions
	std::vector<AD::SXMatrix> inp;
	inp.push_back(XU_vec);
	inp.push_back(P_vec);

	std::vector<AD::SXMatrix> out;
	out.push_back(grad_cost);
	AD::SXFunction grad_cost_fcn(inp,out);
	grad_cost_fcn.init();

	return grad_cost_fcn;
}


void CasadiSystem::setup_casadi_vars(const std::vector<mat>& X, const std::vector<mat>& U,
			const mat& P, double* XU_arr, double* P_arr) {
	int index = 0;
	for(int t=0; t < T; ++t) {
		for(int i=0; i < X[t].n_rows; ++i) {
			XU_arr[index++] = X[t](i);
		}

		if (t < T-1) {
			for(int i=0; i < U[t].n_rows; ++i) {
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
