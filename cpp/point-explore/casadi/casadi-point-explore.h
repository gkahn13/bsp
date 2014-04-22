#ifndef __CASADI_POINT_EXPLORE_H__
#define __CASADI_POINT_EXPLORE_H__

#include "../point-explore.h"

#include <symbolic/casadi.hpp>
#include <symbolic/stl_vector_tools.hpp>
#include <cstdlib>

namespace casadi_point_explore {

namespace AD = CasADi;

AD::SXMatrix dist(AD::SXMatrix a, AD::SXMatrix b) {
	return sqrt(trace(mul(trans(a-b),(a-b))));
}

AD::SXMatrix dynfunc(const AD::SXMatrix& x_t, const AD::SXMatrix& u_t)
{
	AD::SXMatrix x_tp1(X_DIM,1);

	x_tp1(0) = x_t(0) + u_t(0)*DT;
	x_tp1(1) = x_t(1) + u_t(1)*DT;

  	return x_tp1;
}

AD::SXMatrix obsfunc(const AD::SXMatrix& x, const AD::SXMatrix& t)
{
	AD::SXMatrix z(Z_DIM,1);

	AD::SXMatrix d = dist(x, t);
	z(0) = 1.0/(1.0+exp(-alpha*(max_range-d)));

  	return z;
}


AD::SXMatrix gaussLikelihood(const AD::SXMatrix& v) {
	Matrix<R_DIM,R_DIM> Sf, Sf_inv;
	chol(R, Sf);
	Sf_inv = !Sf;

//	float Sf_diag_prod = 1;
//	for(int i=0; i < R_DIM; ++i) { Sf_diag_prod *= Sf(i,i); }
//	float C = pow(2*M_PI, Z_DIM/2)*Sf_diag_prod;

	AD::SXMatrix Sf_inv_casadi(R_DIM,R_DIM);
	for(int i=0; i < R_DIM; ++i) {
		for(int j=0; j < R_DIM; ++j) {
			Sf_inv_casadi(i,j) = Sf_inv(i,j);
		}
	}

//	int rows = Sf_inv.numRows();
//	int cols = Sf_inv.numColumns();
//	AD::SXMatrix M(rows,1);
//	for(int i=0; i < rows; ++i) {
//		for(int j=0; j < cols; ++j) {
//			M(i,0) += Sf_inv(i,j)*v(j,0);
//		}
//	}

	AD::SXMatrix M = mul(Sf_inv_casadi, v);

	AD::SXMatrix E_exp_sum = exp(-0.5*mul(trans(M),M));
	AD::SXMatrix w = E_exp_sum; // / C; // no need to normalize here because normalized later on anyways
	return w;
}

AD::SXMatrix differential_entropy(const std::vector<AD::SXMatrix>& X, const std::vector<AD::SXMatrix>& U,
									const std::vector<AD::SXMatrix>& P) {
	Matrix<R_DIM,R_DIM> Sf, Sf_inv;
	chol(R, Sf);
	Sf_inv = !Sf;

	float Sf_diag_prod = 1;
	for(int i=0; i < R_DIM; ++i) { Sf_diag_prod *= Sf(i,i); }

	AD::SXMatrix Sf_inv_casadi(R_DIM,R_DIM);
	for(int i=0; i < R_DIM; ++i) {
		for(int j=0; j < R_DIM; ++j) {
			Sf_inv_casadi(i,j) = Sf_inv(i,j);
		}
	}

	AD::SXMatrix entropy(1,1);

	std::vector<AD::SXMatrix> X_prop(T);
	std::vector<std::vector<AD::SXMatrix> > H(T, std::vector<AD::SXMatrix>(M));
	std::vector<std::vector<std::vector<AD::SXMatrix> > > GL_H(T, std::vector<std::vector<AD::SXMatrix> >(M, std::vector<AD::SXMatrix>(M)));
	for(int t=0; t < T-1; ++t) {
		X_prop[t+1] = dynfunc(X[t], U[t]);
		for(int m=0; m < M; ++m) {
			H[t+1][m] = obsfunc(X_prop[t+1], P[m]);
//			H[t+1][m] = mul(Sf_inv_casadi, obsfunc(X_prop[t+1], P[m]));
		}
	}

//	AD::SXMatrix v = AD::ssym("v", Z_DIM, 1);
//	AD::SXMatrix gl_output = gaussLikelihood(v);
//
//	std::vector<AD::SXMatrix> inp;
//	inp.push_back(v);
//
//	AD::SXFunction gl_fcn(inp, gl_output);
//	gl_fcn.init();
//
////	std::vector<AD::SXMatrix> c(1, H[0][0]);
//	const AD::MX val;
//	AD::DMatrix dmat(1,1);
//	dmat(0,0) = (AD::SX)(H[0][0](0,0));
//	gl_fcn.call(dmat);

	for(int t=1; t < T; ++t) {
		for(int m=0; m < M; ++m) {
			for(int n=m; n < M; ++n) {
				GL_H[t][m][n] = gaussLikelihood(H[t][m] - H[t][n]);
//				GL_H[t][m][n] = -0.5*mul(H[t][m]-H[t][n],H[t][m]-H[t][n]);
			}
		}
	}

	std::vector<std::vector<AD::SXMatrix> > W(T, std::vector<AD::SXMatrix>(M, AD::SXMatrix(1,1)));
	for(int m=0; m < M; ++m) { W[0][m](0,0) = 1/float(M); }
	for(int t=1; t < T; ++t) {

		AD::SXMatrix W_sum(1,1);
		for(int m=0; m < M; ++m) {
			for(int n=0; n < M; ++n) {
				int a = MIN(m,n);
				int b = MAX(m,n);
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
		entropy += alpha_control_smooth*(mul(trans(U[t+1]-U[t]),(U[t+1]-U[t])));
	}
//
//	for(int t=0; t < T-1; ++t) {
//		entropy += alpha_control_norm*tr(~U[t]*U[t]);
//	}

	return entropy;
}

AD::SXMatrix differential_entropy_wrapper(const AD::SXMatrix& XU_vec, const AD::SXMatrix& P_vec) {
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
		P[m] = P_vec(AD::Slice(index,index+X_DIM));
		index += X_DIM;
	}

	return differential_entropy(X, U, P);
}

void setup_casadi_vars(const std::vector<Matrix<X_DIM> >& X, const std::vector<Matrix<U_DIM> >& U,
					 const std::vector<Matrix<X_DIM> >& P, double* XU_arr, double* P_arr) {
	int index = 0;
	for(int t=0; t < T; ++t) {
		for(int i=0; i < X_DIM; ++i) {
			XU_arr[index++] = X[t][i];
		}

		if (t < T-1) {
			for(int i=0; i < U_DIM; ++i) {
				XU_arr[index++] = U[t][i];
			}
		}
	}

	index = 0;
	for(int m=0; m < M; ++m) {
		for(int i=0; i < X_DIM; ++i) {
			P_arr[index++] = P[m][i];
		}
	}
}

AD::SXFunction casadi_differential_entropy_func() {
	AD::SXMatrix XU_vec = AD::ssym("XU_vec", T*X_DIM + (T-1)*U_DIM);
	AD::SXMatrix P_vec = AD::ssym("P_vec", T*M*X_DIM);

	AD::SXMatrix entropy = differential_entropy_wrapper(XU_vec, P_vec);

	std::vector<AD::SXMatrix> inp;
	inp.push_back(XU_vec);
	inp.push_back(P_vec);

	AD::SXFunction entropy_fcn(inp, entropy);
	entropy_fcn.init();

	return entropy_fcn;
}

AD::SXFunction casadi_differential_entropy_gradfunc() {
	AD::SXMatrix XU_vec = AD::ssym("XU_vec", T*X_DIM + (T-1)*U_DIM);
	AD::SXMatrix P_vec = AD::ssym("P_vec", T*M*X_DIM);

	AD::SXMatrix entropy = differential_entropy_wrapper(XU_vec, P_vec);

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

AD::SXFunction casadi_differential_entropy_diaghessfunc() {
	AD::SXMatrix XU_vec = AD::ssym("XU_vec", T*X_DIM + (T-1)*U_DIM);
	AD::SXMatrix P_vec = AD::ssym("P_vec", T*M*X_DIM);

	AD::SXMatrix entropy = differential_entropy_wrapper(XU_vec, P_vec);

	AD::SXMatrix hess_entropy = hessian(entropy,XU_vec);

	AD::SXMatrix diaghess_entropy = diag(diag(hess_entropy));

	// Create functions
	std::vector<AD::SXMatrix> inp;
	inp.push_back(XU_vec);
	inp.push_back(P_vec);

	std::vector<AD::SXMatrix> out;
	out.push_back(diaghess_entropy);
	AD::SXFunction diaghess_entropy_fcn(inp,out);
	diaghess_entropy_fcn.init();

	return diaghess_entropy_fcn;
}


}

#endif
