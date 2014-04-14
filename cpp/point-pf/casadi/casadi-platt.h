#ifndef __CASADI_PLATT_H__
#define __CASADI_PLATT_H__

#include "casadi-point-pf.h"

namespace casadi_platt {

namespace AD = CasADi;

AD::SXMatrix costfunc_noise(const std::vector<std::vector<AD::SXMatrix> >& P, const std::vector<AD::SXMatrix >& U,
						const std::vector<AD::SXMatrix>& dyn_noise, const std::vector<AD::SXMatrix>& obs_noise) {

	AD::SXMatrix cost(1,1);
	cost(0,0) = 0;

	// need:
	// (T-1)*M dyn_noise
	// M*T obs_noise

	std::vector<std::vector<AD::SXMatrix> > P_prop(T, std::vector<AD::SXMatrix>(M));
	P_prop[0] = P[0];
	for(int t=0; t < T-1; ++t) {
		for(int m=0; m < M; ++m) {
			P_prop[t+1][m] = casadi_point_pf::dynfunc(P[t][m], U[t], dyn_noise[t*M+m]);
		}
	}

	AD::SXMatrix P_t_mle(X_DIM,1), P_t_other(X_DIM,1);
	AD::SXMatrix h_t_mle(Z_DIM,1), h_t_other(Z_DIM,1), diff(Z_DIM,1);
	AD::SXMatrix h_diff(T*Z_DIM,1);
	for(int m=1; m < M; ++m) {
		for(int t=0; t < T; ++t) {
			P_t_mle = P_prop[t][0];
			h_t_mle = casadi_point_pf::obsfunc(P_t_mle, obs_noise[t]);

			P_t_other = P_prop[t][m];
			h_t_other = casadi_point_pf::obsfunc(P_t_other, obs_noise[m*T+t]);
			diff = h_t_mle - h_t_other;
			for(int i=0; i < Z_DIM; ++i) {
				h_diff(t*X_DIM+i,0) = diff(i,0);
			}
		}
		AD::SXMatrix norm = mul(trans(h_diff), h_diff);
		cost(0,0) += (1/float(M))*exp(-norm(0,0));
	}

	for(int t=0; t < T-1; ++t) {
		cost += -alpha_control*mul(trans(U[t]),U[t]);
	}

	return -cost;
}

AD::SXMatrix costfunc(const AD::SXMatrix& PU_vec,
					  const AD::SXMatrix& dyn_noise_vec, const AD::SXMatrix& obs_noise_vec) {

	std::vector<std::vector<AD::SXMatrix> > P(T, std::vector<AD::SXMatrix>(M));
	std::vector<AD::SXMatrix> U(T-1);
	std::vector<AD::SXMatrix> dyn_noise((T-1)*M);
	std::vector<AD::SXMatrix> obs_noise(T*M);

	int index = 0;
	for(int t=0; t < T; ++t) {
		for(int m=0; m < M; ++m) {
			AD::SXMatrix p(X_DIM,1);
			for(int i=0; i < X_DIM; ++i) {
				p(i,0) = PU_vec(index++, 0);
			}
			P[t][m] = p;
		}
		if (t < T-1) {
			AD::SXMatrix u(U_DIM,1);
			for(int i=0; i < U_DIM; ++i) {
				u(i,0) = PU_vec(index++,0);
			}
			U[t] = u;
		}
	}
//
//	for(int t=0; t < T-1; ++t) {
//		AD::SXMatrix u(U_DIM,1);
//		for(int i=0; i < U_DIM; ++i) {
//			u(i,0) = PU_vec(index++,0);
//		}
//		U[t] = u;
//	}

	index = 0;
	for(int t=0; t < T-1; ++t) {
		for(int m=0; m < M; ++m) {
			AD::SXMatrix d(Q_DIM,1);
			for(int i=0; i < Q_DIM; ++i) {
				d(i,0) = dyn_noise_vec(index++, 0);
			}
			dyn_noise[t*M+m] = d;
		}
	}

	index = 0;
	for(int m=0; m < M; ++m) {
		for(int t=0; t < T; ++t) {
			AD::SXMatrix r(R_DIM,1);
			for(int i=0; i < R_DIM; ++i) {
				r(i,0) = obs_noise_vec(index++, 0);
			}
			obs_noise[m*T+t] = r;
		}
	}

	return costfunc_noise(P, U, dyn_noise, obs_noise);
}

AD::SXFunction casadi_costfunc() {
	AD::SXMatrix PU_vec = AD::ssym("PU_vec",T*M*X_DIM+(T-1)*U_DIM);
	AD::SXMatrix dyn_noise_vec = AD::ssym("dyn_noise_vec", (T-1)*M*Q_DIM);
	AD::SXMatrix obs_noise_vec = AD::ssym("obs_noise_vec", M*T*R_DIM);

	AD::SXMatrix cost = costfunc(PU_vec, dyn_noise_vec, obs_noise_vec);

	std::vector<AD::SXMatrix> inp;
	inp.push_back(PU_vec);
	inp.push_back(dyn_noise_vec);
	inp.push_back(obs_noise_vec);

	AD::SXFunction cost_fcn(inp, cost);
	cost_fcn.init();

	return cost_fcn;
}

AD::SXFunction casadi_costgradfunc() {
	AD::SXMatrix PU_vec = AD::ssym("PU_vec",T*M*X_DIM+(T-1)*U_DIM);
	AD::SXMatrix dyn_noise_vec = AD::ssym("dyn_noise_vec", (T-1)*M*Q_DIM);
	AD::SXMatrix obs_noise_vec = AD::ssym("obs_noise_vec", M*T*R_DIM);

	AD::SXMatrix cost = costfunc(PU_vec, dyn_noise_vec, obs_noise_vec);

	AD::SXMatrix grad_cost = gradient(cost,PU_vec);

	// Create functions
	std::vector<AD::SXMatrix> inp;
	inp.push_back(PU_vec);
	inp.push_back(dyn_noise_vec);
	inp.push_back(obs_noise_vec);

	std::vector<AD::SXMatrix> out;
	out.push_back(cost);
	out.push_back(grad_cost);
	AD::SXFunction grad_f_fcn(inp,out);
	grad_f_fcn.init();

	return grad_f_fcn;
}

AD::SXFunction casadi_costdiaghessfunc() {
	AD::SXMatrix PU_vec = AD::ssym("PU_vec",T*M*X_DIM+(T-1)*U_DIM);
	AD::SXMatrix dyn_noise_vec = AD::ssym("dyn_noise_vec", (T-1)*M*Q_DIM);
	AD::SXMatrix obs_noise_vec = AD::ssym("obs_noise_vec", M*T*R_DIM);

	AD::SXMatrix cost = costfunc(PU_vec, dyn_noise_vec, obs_noise_vec);

	AD::SXMatrix hess_cost = hessian(cost, PU_vec);

	AD::SXMatrix diaghess_cost = diag(diag(hess_cost));

	// Create functions
	std::vector<AD::SXMatrix> inp;
	inp.push_back(PU_vec);
	inp.push_back(dyn_noise_vec);
	inp.push_back(obs_noise_vec);

	std::vector<AD::SXMatrix> out;
	out.push_back(cost);
	out.push_back(diaghess_cost);
	AD::SXFunction diaghess_f_fcn(inp,out);
	diaghess_f_fcn.init();

	return diaghess_f_fcn;
}

void setup_casadi_vars(const std::vector<std::vector<Matrix<X_DIM> > >& P, const std::vector<Matrix<U_DIM> >& U,
						  const std::vector<Matrix<Q_DIM> >& dyn_noise, const std::vector<Matrix<R_DIM> >& obs_noise,
						  double* PU_arr, double* dyn_noise_arr, double* obs_noise_arr) {
	int index = 0;
	for(int t=0; t < T; ++t) {
		for(int m=0; m < M; ++m) {
			for(int i=0; i < X_DIM; ++i) {
				PU_arr[index++] = P[t][m][i];
			}
		}
		if (t < T-1) {
			for(int i=0; i < U_DIM; ++i) {
				PU_arr[index++] = U[t][i];
			}
		}
	}
//
//	for(int t=0; t < T-1; ++t) {
//		for(int i=0; i < U_DIM; ++i) {
//			PU_arr[index++] = U[t][i];
//		}
//	}

	index = 0;
	for(int t=0; t < T-1; ++t) {
		for(int m=0; m < M; ++m) {
			for(int i=0; i < Q_DIM; ++i) {
				dyn_noise_arr[index++] = dyn_noise[t*M+m][i];
			}
		}
	}

	index = 0;
	for(int m=0; m < M; ++m) {
		for(int t=0; t < T; ++t) {
			for(int i=0; i < R_DIM; ++i) {
				obs_noise_arr[index++] = obs_noise[m*T+t][i];
			}
		}
	}
}

}

#endif
