#ifndef __POINT_ENTROPY_H__
#define __POINT_ENTROPY_H__

#include "point-pf.h"

namespace point_entropy {

float differential_entropy_noise(const std::vector<std::vector<Matrix<X_DIM> > >& P, const std::vector<Matrix<U_DIM> >& U,
							         const std::vector<Matrix<Q_DIM> >& dyn_noise, const std::vector<Matrix<R_DIM> >& obs_noise) {
	float entropy = 0;

	std::vector<std::vector<Matrix<X_DIM> > > P_prop(T, std::vector<Matrix<X_DIM> >(M));
	std::vector<std::vector<Matrix<X_DIM> > > P_prop_mle(T, std::vector<Matrix<X_DIM> >(M));
	std::vector<std::vector<Matrix<Z_DIM> > > H(T, std::vector<Matrix<Z_DIM> >(M));
	std::vector<std::vector<Matrix<Z_DIM> > > H_mle(T, std::vector<Matrix<Z_DIM> >(M));
	P_prop[0] = P[0];
	for(int t=0; t < T-1; ++t) {
		for(int m=0; m < M; ++m) {
			P_prop[t+1][m] = point_pf::dynfunc(P[t][m], U[t], dyn_noise[t*M+m]);
			P_prop_mle[t+1][m] = point_pf::dynfunc(P[t][m], U[t], zeros<Q_DIM,1>());
			H[t+1][m] = point_pf::obsfunc(P_prop[t+1][m], obs_noise[t*M+m]);
			H_mle[t+1][m] = point_pf::obsfunc(P_prop[t+1][m], zeros<R_DIM,1>());
		}
	}

	std::vector<std::vector<float> > W(T, std::vector<float>(M,1/float(M)));
	std::vector<std::vector<std::vector<float> > > dyn_prob(T, std::vector<std::vector<float> >(M, std::vector<float>(M)));
	for(int t=1; t < T; ++t) {
		// gauss likelihood of each particle at P[t]
		// set W[t]
		float W_sum = 0;
		for(int m=0; m < M; ++m) {
			W[t][m] = point_pf::gaussLikelihood<X_DIM, R_DIM>(H[t][m] - H_mle[t][m], R);
			W_sum += W[t][m];
		}
		for(int m=0; m < M; ++m) { W[t][m] = W[t][m] / W_sum; }

		// gauss likelihood of dynamics
		for(int m=0; m < M; ++m) {
			float dyn_prob_sum = 0;
			for(int n=0; n < M; ++n) {
				dyn_prob[t][m][n] = point_pf::gaussLikelihood<X_DIM, Q_DIM>(P_prop_mle[t][m] - P_prop_mle[t][n], Q);
				dyn_prob_sum += dyn_prob[t][m][n];
			}
			for(int n=0; n < M; ++n) { dyn_prob[t][m][n] = dyn_prob[t][m][n] / dyn_prob_sum; }
		}

		// use skoglar version
		float entropy_t = 0;
		for(int m=0; m < M; ++m) {
			entropy_t += -W[t][m]*log(W[t][m]);
		}
		for(int m=0; m < M; ++m) {
			float weighted_dyn_prob = 0;
			for(int n=0; n < M; ++n) {
				weighted_dyn_prob += W[t-1][n]*dyn_prob[t][m][n];
			}
			entropy_t += -W[t][m]*log(weighted_dyn_prob);
		}
		float sum_cross_time_weights = 0;
		for(int m=0; m < M; ++m) {
			sum_cross_time_weights += W[t-1][m]*W[t][m];
		}
		entropy_t += log(sum_cross_time_weights);

		entropy += entropy_t;
	}

	for(int t=0; t < T-1; ++t) {
		entropy -= alpha_control*tr(~U[t]*U[t]);
	}

	return -entropy;
}

float differential_entropy(std::vector<std::vector<Matrix<X_DIM> > >& P, std::vector<Matrix<U_DIM> >& U) {
	float avg_entropy = 0;
	float max_iter = 1000;
	for(int iter=0; iter < max_iter; iter++) {
		std::vector<Matrix<Q_DIM> > dyn_noise = sampleGaussianN(zeros<Q_DIM,1>(), Q, (T-1)*M);
		std::vector<Matrix<R_DIM> > obs_noise = sampleGaussianN(zeros<R_DIM,1>(), R, (T*M));
		float entropy = point_entropy::differential_entropy_noise(P, U, dyn_noise, obs_noise);
		avg_entropy += (1/float(max_iter))*entropy;
	}
	return avg_entropy;
}

Matrix<TOTAL_VARS> grad_differential_entropy(std::vector<std::vector<Matrix<X_DIM> > >& P, std::vector<Matrix<U_DIM> >& U) {
	Matrix<TOTAL_VARS> g_avg;

	int max_iter = 20;
	for(int iter=0; iter < max_iter; iter++) {
		Matrix<TOTAL_VARS> g;

		std::vector<Matrix<Q_DIM> > dyn_noise = sampleGaussianN(zeros<Q_DIM,1>(), Q, (T-1)*M);
		std::vector<Matrix<R_DIM> > obs_noise = sampleGaussianN(zeros<R_DIM,1>(), R, T*M);

		float orig, entropy_p, entropy_l;
		int index = 0;
		for(int t=0; t < T; ++t) {
			for(int m=0; m < M; ++m) {
				for(int i=0; i < X_DIM; ++i) {
					orig = P[t][m][i];

					P[t][m][i] = orig + step;
					entropy_p = differential_entropy_noise(P,U, dyn_noise, obs_noise);

					P[t][m][i] = orig - step;
					entropy_l = differential_entropy_noise(P,U, dyn_noise, obs_noise);

					P[t][m][i] = orig;
					g[index++] = (entropy_p - entropy_l)/(2*step);
				}
			}

			if (t < T-1) {
				for(int i=0; i < U_DIM; ++i) {
					orig = U[t][i];

					U[t][i] = orig + step;
					entropy_p = differential_entropy_noise(P,U, dyn_noise, obs_noise);

					U[t][i] = orig - step;
					entropy_l = differential_entropy_noise(P,U, dyn_noise, obs_noise);

					U[t][i] = orig;
					g[index++] = (entropy_p - entropy_l)/(2*step);
				}
			}
		}

		g_avg += (1/float(max_iter))*g;
	}

	return g_avg;
}

Matrix<TOTAL_VARS> diaghess_differential_entropy(std::vector<std::vector<Matrix<X_DIM> > >& P, std::vector<Matrix<U_DIM> >& U) {
	Matrix<TOTAL_VARS> diaghess_avg;

	int max_iter = 20;
	for(int iter=0; iter < max_iter; iter++) {
		Matrix<TOTAL_VARS> diaghess;

		std::vector<Matrix<Q_DIM> > dyn_noise = sampleGaussianN(zeros<Q_DIM,1>(), Q, (T-1)*M);
		std::vector<Matrix<R_DIM> > obs_noise = sampleGaussianN(zeros<R_DIM,1>(), R, T*M);
		std::vector<float> sampling_noise(T);
		for(int t=0; t < T; ++t) { sampling_noise[t] = (1/float(M))*(rand() / float(RAND_MAX)); }

		float orig, cost_p, cost, cost_l;
		int index = 0;
		for(int t=0; t < T; ++t) {
			for(int m=0; m < M; ++m) {
				for(int i=0; i < X_DIM; ++i) {
					cost = differential_entropy_noise(P, U, dyn_noise, obs_noise);
					orig = P[t][m][i];

					P[t][m][i] = orig + step;
					cost_p = differential_entropy_noise(P, U, dyn_noise, obs_noise);

					P[t][m][i] = orig - step;
					cost_l = differential_entropy_noise(P,U, dyn_noise, obs_noise);

					P[t][m][i] = orig;
					diaghess[index++] = (cost_p - 2*cost + cost_l)/(step*step);
				}
			}

			if (t < T-1) {
				for(int i=0; i < U_DIM; ++i) {
					cost = differential_entropy_noise(P, U, dyn_noise, obs_noise);
					orig = U[t][i];

					U[t][i] = orig + step;
					cost_p = differential_entropy_noise(P,U, dyn_noise, obs_noise);

					U[t][i] = orig - step;
					cost_l = differential_entropy_noise(P,U, dyn_noise, obs_noise);

					U[t][i] = orig;
					diaghess[index++] = (cost_p - 2*cost + cost_l)/(step*step);
				}
			}
		}

		diaghess_avg += (1/float(max_iter))*diaghess;
	}

	return diaghess_avg;
}


}

#endif
