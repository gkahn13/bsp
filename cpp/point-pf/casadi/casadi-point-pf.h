#ifndef __CASADI_POINT_PF_H__
#define __CASADI_POINT_PF_H__

#include "../point-pf.h"
#include "../util/matrix.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <ctime>
#include <symbolic/casadi.hpp>
#include <symbolic/stl_vector_tools.hpp>
#include <cstdlib>

namespace casadi_point_pf {

//using namespace CasADi;
namespace AD = CasADi;

AD::SXMatrix dynfunc(const AD::SXMatrix& x_t, const AD::SXMatrix& u_t, const AD::SXMatrix& q_t)
{
	AD::SXMatrix x_tp1(X_DIM,1);

	x_tp1(0) = x_t(0) + u_t(0)*DT + Q(0,0)*q_t(0);
	x_tp1(1) = x_t(1) + u_t(1)*DT + Q(1,1)*q_t(1);

  	return x_tp1;
}

AD::SXMatrix obsfunc(const AD::SXMatrix& x_t, const AD::SXMatrix& r_t)
{
	AD::SXMatrix z(Z_DIM,1);

	z(0) = x_t(0) + sqrt(0.5*0.5*x_t(0)*x_t(0) + 1e-6)*r_t(0);
	z(1) = x_t(1) + sqrt(0.5*0.5*x_t(1)*x_t(1) + 1e-6)*r_t(1);

  	return z;
}

AD::SXMatrix gaussLikelihood(const AD::SXMatrix& v, const AD::SXMatrix& Sf_inv, const AD::SXMatrix& C) {
	AD::SXMatrix M = mul(Sf_inv, v);

	AD::SXMatrix E_exp_sum = exp(-0.5*trace(mul(trans(M),M)));
	AD::SXMatrix w = E_exp_sum / C;
	return w;
}

std::vector<AD::SXMatrix> lowVarianceSampler(const std::vector<AD::SXMatrix>& P, const std::vector<AD::SXMatrix>& W, AD::SXMatrix r) {
	int M = P.size();
	std::vector<AD::SXMatrix> P_sampled(M);

	AD::SXMatrix r_round = floor(r);
	isMinusOne(r);

	AD::SXMatrix c = W[0];
	int i = 0;
	for(int m=0; m < M; ++m) {
		AD::SXMatrix u = r + (m) * (1/float(M));
		AD::SXMatrix diff = u(0,0) - c(0,0);
		while (isOne(diff/fabs(diff))) {
		//while(u(0,0) > c(0,0)) {
			c += W[++i];
		}
		P_sampled[m] = P[i];
	}

	return P_sampled;
}

std::vector<AD::SXMatrix> beliefDynamics(const std::vector<AD::SXMatrix>& P_t, const AD::SXMatrix& u,
									   const std::vector<AD::SXMatrix>& dyn_noise, const std::vector<AD::SXMatrix>& obs_noise, const AD::SXMatrix& sampling_noise,
									   const AD::SXMatrix& Sf_inv, const AD::SXMatrix& C) {
	int M = P_t.size();
	std::vector<AD::SXMatrix> P_tp1_bar(M), P_tp1;
	std::vector<AD::SXMatrix> W(M);

	AD::SXMatrix W_sum = 0;
	for(int m=0; m < M; ++m) {
		P_tp1_bar[m] = dynfunc(P_t[m], u, dyn_noise[m]);
		AD::SXMatrix e = obsfunc(P_tp1_bar[m], obs_noise[m]) - obsfunc(P_tp1_bar[m], AD::SXMatrix::zeros(R_DIM,1));
		W[m] = gaussLikelihood(e, Sf_inv, C);
		W_sum += W[m];
	}
	for(int m=0; m < M; ++m) { W[m] = W[m] / W_sum; }

	P_tp1 = lowVarianceSampler(P_tp1_bar, W, sampling_noise);

	return P_tp1;
}

AD::SXMatrix beliefDynamicsWrapper(const AD::SXMatrix& P_t_u, const AD::SXMatrix& dyn_noise_mat, const AD::SXMatrix& obs_noise_mat,
									  const AD::SXMatrix& sampling_noise, const int M) {
	std::vector<AD::SXMatrix> P_t(M), dyn_noise(M), obs_noise(M);
	int index = 0;
	for(int m=0; m < M; ++m) {
		// below indexing only works since X_DIM == Q_DIM == R_DIM
		P_t[m] = P_t_u(AD::Slice(index,index+X_DIM));
		dyn_noise[m] = dyn_noise_mat(AD::Slice(index,index+Q_DIM));
		obs_noise[m] = obs_noise_mat(AD::Slice(index,index+R_DIM));
		index += X_DIM;
	}
	AD::SXMatrix u = P_t_u(AD::Slice(index,index+U_DIM));

	Matrix<R_DIM,R_DIM> Sf, Sf_inv;
	chol(R, Sf);
	Sf_inv = !Sf;

	float Sf_diag_prod = 1;
	for(int i=0; i < R_DIM; ++i) { Sf_diag_prod *= Sf(i,i); }
	float C = pow(2*M_PI, Z_DIM/2)*Sf_diag_prod;

	AD::SXMatrix Sf_inv_casadi(R_DIM,R_DIM), C_casadi(1,1);
	for(int i=0; i < R_DIM; ++i) {
		for(int j=0; j < R_DIM; ++j) {
			Sf_inv_casadi(i,j) = Sf_inv(i,j);
		}
	}
	C_casadi(0,0) = C;

	std::vector<AD::SXMatrix> P_tp1 = beliefDynamics(P_t, u, dyn_noise, obs_noise, sampling_noise, Sf_inv_casadi, C_casadi);
	AD::SXMatrix P_tp1_vec(T*M,1);
	index = 0;
	for(int m=0; m < M; ++m) {
		for(int i=0; i < X_DIM; ++i) {
			P_tp1_vec(index++) = P_tp1[m][i];
		}
	}

	return P_tp1_vec;
}

AD::SXFunction casadiBeliefDynamicsFunc(int M) {
	AD::SXMatrix P_t_u = AD::ssym("P_t_u",M*X_DIM+U_DIM);
	AD::SXMatrix dyn_noise_mat = AD::ssym("dyn_noise_mat", M*X_DIM);
	AD::SXMatrix obs_noise_mat = AD::ssym("obs_noise_mat", M*X_DIM);
	AD::SXMatrix sampling_noise = AD::ssym("sampling_noise", 1);

	AD::SXMatrix P_tp1 = beliefDynamicsWrapper(P_t_u, dyn_noise_mat, obs_noise_mat, sampling_noise, M);

	std::vector<AD::SXMatrix> inp;
	inp.push_back(P_t_u);
	inp.push_back(dyn_noise_mat);
	inp.push_back(obs_noise_mat);
	inp.push_back(sampling_noise);

	AD::SXFunction bd_fcn(inp, P_tp1);
	bd_fcn.init();

	return bd_fcn;
}


}

#endif
