#include "system.h"

/**
 *
 * Constructors
 *
 */

System::System() { }

/**
 *
 * Public methods
 *
 */

mat System::cost_grad(std::vector<mat>& X, std::vector<mat>& U, const mat& P) {
	if (this->use_casadi) {
		return this->casadi_sys->casadi_cost_grad(X, U, P);
	}

	mat g(TOTAL_VARS, 1, fill::zeros);

	double orig, entropy_p, entropy_l;
	int index = 0;
	for(int t=0; t < T; ++t) {
		for(int i=0; i < X[t].n_rows; ++i) {
			orig = X[t](i);

			X[t](i) = orig + step;
			entropy_p = this->cost(X, U, P);

			X[t](i) = orig - step;
			entropy_l = this->cost(X, U, P);

			X[t][i] = orig;
			g(index++) = (entropy_p - entropy_l)/(2*step);
		}

		if (t < T-1) {
			for(int i=0; i < U[t].n_rows; ++i) {
				orig = U[t](i);

				U[t](i) = orig + step;
				entropy_p = this->cost(X, U, P);

				U[t](i) = orig - step;
				entropy_l = this->cost(X, U, P);

				U[t][i] = orig;
				g(index++) = (entropy_p - entropy_l)/(2*step);
			}
		}

	}

	return g;
}

/**
 *
 * Private methods
 *
 */

/**
 * S -- symmetric
 */
double System::gauss_likelihood(const mat& v, const mat& S) {
	mat Sf = chol(S);
	mat M = solve(Sf, v);

	double E = -0.5*accu(M % M);
	double C = pow(2*M_PI, S.n_cols/2) * prod(diagvec(Sf));
	double w = exp(E) / C;

	return w;
}

mat System::low_variance_sampler(const mat& P, const mat& W, double r) {
	int M = P.n_cols;
	mat P_sampled(P.n_rows, P.n_cols, fill::zeros);

	double c = W(0);
	int i = 0;
	for(int m=0; m < M; ++m) {
		double u = r + m * (1/double(M));
		while (u > c) {
			c += W(++i);
		}
		P_sampled.col(m) = P.col(i);
	}

	return P_sampled;
}


double System::cost_entropy(const std::vector<mat>& X, const std::vector<mat>& U, const mat& P) {
	int T = X.size();
	int M = P.n_cols;

	double entropy = 0;

	std::vector<mat> X_prop(T);
	std::vector<mat> H(T, zeros<mat>(Z_DIM, M));
	mat r(Z_DIM, 1, fill::zeros);
	for(int t=0; t < T-1; ++t) {
		X_prop[t+1] = this->dynfunc(X[t], U[t]);
		for(int m=0; m < M; ++m) {
			H[t+1].col(m) = this->obsfunc(X_prop[t+1], P.col(m), r);
		}
	}

	std::vector<mat> W(T, zeros<mat>(M, 1));
	W[0] = (1/double(M))*ones<mat>(M,1);
	for(int t=1; t < T; ++t) {
		for(int m=0; m < M; ++m) {
			for(int p=0; p < M; ++p) {
				mat diff = H[t].col(m) - H[t].col(p);
				W[t](m) += this->gauss_likelihood(diff, this->R);
			}
		}
		W[t] = W[t] / accu(W[t]);

		// skoglar version
		// second term simplifies because zero particle dynamics
		entropy += accu(-W[t] % log(W[t])) +
				   accu(-W[t] % log(W[t-1])) +
				   log(accu(W[t-1] % W[t]));
	}

	return entropy;
}

double System::cost_platt(const std::vector<mat>& X, const std::vector<mat>& U, const mat& P) {
	int T = X.size();
	int M = P.n_cols;

	std::vector<mat> X_prop(T, zeros<mat>(X[0].n_rows,1));
	X_prop[0] = X[0];
	for(int t=0; t < T-1; ++t) {
		X_prop[t+1] = this->dynfunc(X[t], U[t]);
	}

	std::vector<mat> H(M, zeros<mat>(T*(Z_DIM),1));
	mat r(Z_DIM, 1, fill::zeros);
	for(int m=0; m < M; ++m) {
		for(int t=0; t < T; ++t) {
			H[m].rows(t*(Z_DIM), (t+1)*(Z_DIM)-1) = this->obsfunc(X_prop[t], P.col(m), r);
		}
	}

	double platt = 0;
	for(int m=1; m < M; ++m) {
		mat diff = H[m] - H[0];
		platt += (1/float(M-1))*exp(-accu(diff % diff));
	}

	return platt;
}
