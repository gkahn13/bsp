#ifndef __POINT_ENTROPY_H__
#define __POINT_ENTROPY_H__

#include "point-pf.h"

namespace point_entropy {

float differential_entropy(const std::vector<std::vector<Matrix<X_DIM> > >& P, const std::vector<Matrix<U_DIM> >& U) {
	float entropy = 0;

	std::vector<std::vector<Matrix<X_DIM> > > P_prop(T, std::vector<Matrix<X_DIM> >(M));
	std::vector<std::vector<Matrix<X_DIM> > > P_prop_mle(T, std::vector<Matrix<X_DIM> >(M));
	std::vector<std::vector<Matrix<Z_DIM> > > H(T, std::vector<Matrix<Z_DIM> >(M));
	std::vector<std::vector<Matrix<Z_DIM> > > H_mle(T, std::vector<Matrix<Z_DIM> >(M));
	std::vector<Matrix<Q_DIM> > dyn_noise = sampleGaussianN(zeros<Q_DIM,1>(), Q, (T-1)*M);
	std::vector<Matrix<R_DIM> > obs_noise = sampleGaussianN(zeros<R_DIM,1>(), R, (T*M));
	P_prop[0] = P[0];
	for(int t=0; t < T-1; ++t) {
		for(int m=0; m < M; ++m) {
			P_prop[t+1][m] = point_pf::dynfunc(P[t][m], U[t], dyn_noise[t*M+m]);
			P_prop_mle[t+1][m] = point_pf::dynfunc(P[t][m], U[t], zeros<Q_DIM,1>());
			H[t+1][m] = point_pf::obsfunc(P_prop[t+1][m], obs_noise[t*M+m]);
			H_mle[t+1][m] = point_pf::obsfunc(P_prop[t+1][m], zeros<R_DIM,1>());
		}
	}

//	std::cout << "Propagated everything\n";

	std::vector<std::vector<float> > W(T, std::vector<float>(M,1/float(M)));
	std::vector<std::vector<std::vector<float> > > dyn_prob(T, std::vector<std::vector<float> >(M, std::vector<float>(M)));
	for(int t=1; t < T; ++t) {
//		std::cout << "t: " << t << "\n";
		// gauss likelihood of each particle at P[t]
		// set W[t]
		float W_sum = 0;
		for(int m=0; m < M; ++m) {
			W[t][m] = point_pf::gaussLikelihood<X_DIM, R_DIM>(H[t][m] - H_mle[t][m], R);
			W_sum += W[t][m];
		}
		for(int m=0; m < M; ++m) { W[t][m] = W[t][m] / W_sum; }

//		std::cout << "Gauss likelihood dynamics\n";
		// gauss likelihood of dynamics
		for(int m=0; m < M; ++m) {
			float dyn_prob_sum = 0;
			for(int n=0; n < M; ++n) {
				dyn_prob[t][m][n] = point_pf::gaussLikelihood<X_DIM, Q_DIM>(P_prop_mle[t][m] - P_prop_mle[t][n], Q);
				dyn_prob_sum += dyn_prob[t][m][n];
			}
			for(int n=0; n < M; ++n) { dyn_prob[t][m][n] = dyn_prob[t][m][n] / dyn_prob_sum; }
		}

//		std::cout << "Calculate skoglar entropy\n";
		// use skoglar version
		float entropy_t = 0;
//		std::cout << "part 1\n";
		for(int m=0; m < M; ++m) {
			entropy_t += -W[t][m]*log(W[t][m]);
		}
//		std::cout << "part 2\n";
		for(int m=0; m < M; ++m) {
//			std::cout << "m: " << m << "\n";
			float weighted_dyn_prob = 0;
			for(int n=0; n < M; ++n) {
//				std::cout << "n: " << n << "\n";
				weighted_dyn_prob += W[t-1][n]*dyn_prob[t][m][n];
			}
			entropy_t += -W[t][m]*log(weighted_dyn_prob);
		}
//		std::cout << "part 3\n";
		float sum_cross_time_weights = 0;
		for(int m=0; m < M; ++m) {
			sum_cross_time_weights += W[t-1][m]*W[t][m];
		}
		entropy_t += log(sum_cross_time_weights);

//		std::cout << "entropy_" << t << ": " << entropy_t << "\n";
		entropy += entropy_t;
	}

	return entropy;
}

}

#endif
