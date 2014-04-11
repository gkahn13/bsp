#ifndef __POINT_PLATT_H__
#define __POINT_PLATT_H__

#include "point-pf.h"
#include "casadi/casadi-platt.h"

namespace AD = CasADi;

namespace point_platt {

// assume P[t][0] is most likely particle
float costfunc_noise(const std::vector<std::vector<Matrix<X_DIM> > >& P, const std::vector<Matrix<U_DIM> >& U,
						const std::vector<Matrix<Q_DIM> >& dyn_noise, const std::vector<Matrix<R_DIM> >& obs_noise,
						const std::vector<float> sampling_noise) {
	float cost = 0;

	// need:
	// (T-1)*M dyn_noise
	// M*T obs_noise
	// T sampling_noise

	// TODO: if no propagation, then all particle gradients = 0
	std::vector<std::vector<Matrix<X_DIM> > > P_prop(T, std::vector<Matrix<X_DIM> >(M));
	P_prop[0] = P[0];
	for(int t=0; t < T-1; ++t) {
//		std::vector<Matrix<Q_DIM> > dyn_noise_t(dyn_noise.begin()+t*M, dyn_noise.begin()+(t+1)*M);
//		std::vector<Matrix<R_DIM> > obs_noise_t(obs_noise.begin()+t*M, obs_noise.begin()+(t+1)*M);
//		float sampling_noise_t = sampling_noise[t];
//		P_prop[t+1] = point_pf::beliefDynamics(P[t], U[t], dyn_noise_t, obs_noise_t, sampling_noise_t);
		for(int m=0; m < M; ++m) {
			//P_prop[t+1][m] = point_pf::dynfunc(P[t][m], U[t], zeros<Q_DIM,1>());
			P_prop[t+1][m] = point_pf::dynfunc(P[t][m], U[t], dyn_noise[t*M+m]);
		}
	}

	Matrix<X_DIM> P_t_mle, P_t_other;
	Matrix<Z_DIM> h_t_mle, h_t_other;
	Matrix<T*Z_DIM> h_diff;
	for(int m=1; m < M; ++m) {
		for(int t=0; t < 1; ++t) { // T
			P_t_mle = P_prop[t][0];
			//h_t_mle = point_pf::obsfunc(P_t_mle, zeros<R_DIM,1>());
			h_t_mle = point_pf::obsfunc(P_t_mle, obs_noise[t]);

			P_t_other = P_prop[t][m];
			//h_t_other = point_pf::obsfunc(P_t_other, zeros<R_DIM,1>());
			h_t_other = point_pf::obsfunc(P_t_other, obs_noise[m*T+t]);
			h_diff.insert(t*X_DIM, 0, h_t_mle - h_t_other);
		}
		cost += (1/float(M))*exp(-tr(~h_diff*h_diff));
	}

	for(int t=0; t < T-1; ++t) {
		cost += alpha_control*tr(~U[t]*U[t]);
	}

	return cost;
}

float costfunc(const std::vector<std::vector<Matrix<X_DIM> > >& P, const std::vector<Matrix<U_DIM> >& U) {
	std::vector<Matrix<Q_DIM> > dyn_noise((T-1)*M);
	std::vector<Matrix<R_DIM> > obs_noise(T*M);
	std::vector<float> sampling_noise(T);

	float cost = 0;
	int max_iter = 100;
	for(int iter=0; iter < max_iter; ++iter) {
		dyn_noise = sampleGaussianN(zeros<Q_DIM,1>(), Q, (T-1)*M);
		obs_noise = sampleGaussianN(zeros<R_DIM,1>(), R, T*M);
		for(int t=0; t < T; ++t) { sampling_noise[t] = (1/float(M))*(rand() / float(RAND_MAX)); }
		cost += (1/float(max_iter))*costfunc_noise(P, U, dyn_noise, obs_noise, sampling_noise);
	}

	return cost;
	//return costfunc_noise(P, U, dyn_noise, obs_noise);
}

float casadi_costfunc(const std::vector<std::vector<Matrix<X_DIM> > >& P, const std::vector<Matrix<U_DIM> >& U) {

	double P_arr[T*M*X_DIM];
	double U_arr[(T-1)*U_DIM];
	double dyn_noise_arr[(T-1)*M*Q_DIM];
	double obs_noise_arr[M*T*R_DIM];

	int index = 0;
	for(int t=0; t < T; ++t) {
		for(int m=0; m < M; ++m) {
			for(int i=0; i < X_DIM; ++i) {
				P_arr[index++] = P[t][m][i];
			}
		}
	}

	index = 0;
	for(int t=0; t < T-1; ++t) {
		for(int i=0; i < U_DIM; ++i) {
			U_arr[index++] = U[t][i];
		}
	}

	std::vector<Matrix<Q_DIM> > dyn_noise = sampleGaussianN(zeros<Q_DIM,1>(), Q, (T-1)*M);
	index = 0;
	for(int t=0; t < T-1; ++t) {
		for(int m=0; m < M; ++m) {
			for(int i=0; i < Q_DIM; ++i) {
				dyn_noise_arr[index++] = dyn_noise[t*M+m][i];
			}
		}
	}

	std::vector<Matrix<Q_DIM> > obs_noise = sampleGaussianN(zeros<Q_DIM,1>(), Q, T*M);
	index = 0;
	for(int m=0; m < M; ++m) {
		for(int t=0; t < T; ++t) {
			for(int i=0; i < R_DIM; ++i) {
				obs_noise_arr[index++] = obs_noise[m*T+t][i];
			}
		}
	}


	double real_cost = costfunc_noise(P,U,dyn_noise,obs_noise,std::vector<float>(T));
	std::cout << "Actual cost: " << real_cost << "\n";

	AD::SXFunction func = casadi_platt::casadiCostfuncFunc();

	func.setInput(P_arr,0);
	func.setInput(U_arr,1);
	func.setInput(dyn_noise_arr,2);
	func.setInput(obs_noise_arr,3);

	func.evaluate();

	double cost = 0;
	func.getOutput(&cost,0);

	std::cout << "Casadi cost: " << cost << "\n";

	return cost;
}

Matrix<TOTAL_VARS> deriv_costfunc(std::vector<std::vector<Matrix<X_DIM> > >& P, std::vector<Matrix<U_DIM> >& U) {
	Matrix<TOTAL_VARS> d_avg;

	int max_iter = 100;
	for(int iter=0; iter < max_iter; iter++) {
		Matrix<TOTAL_VARS> d;

		std::vector<Matrix<Q_DIM> > dyn_noise = sampleGaussianN(zeros<Q_DIM,1>(), Q, (T-1)*M);
		std::vector<Matrix<R_DIM> > obs_noise = sampleGaussianN(zeros<R_DIM,1>(), R, T*M);
		std::vector<float> sampling_noise(T);
		for(int t=0; t < T; ++t) { sampling_noise[t] = (1/float(M))*(rand() / float(RAND_MAX)); }

		float orig, cost_p, cost_l;
		int index = 0;
		for(int t=0; t < T; ++t) {
			for(int m=0; m < M; ++m) {
				for(int i=0; i < X_DIM; ++i) {
					orig = P[t][m][i];

					P[t][m][i] = orig + step;
					cost_p = costfunc_noise(P,U, dyn_noise, obs_noise, sampling_noise);

					P[t][m][i] = orig - step;
					cost_l = costfunc_noise(P,U, dyn_noise, obs_noise, sampling_noise);

					P[t][m][i] = orig;
					d[index++] = (cost_p - cost_l)/(2*step);
				}
			}

			if (t < T-1) {
				for(int i=0; i < U_DIM; ++i) {
					orig = U[t][i];

					U[t][i] = orig + step;
					cost_p = costfunc_noise(P,U, dyn_noise, obs_noise, sampling_noise);

					U[t][i] = orig - step;
					cost_l = costfunc_noise(P,U, dyn_noise, obs_noise, sampling_noise);

					U[t][i] = orig;
					d[index++] = (cost_p - cost_l)/(2*step);
				}
			}
		}

		d_avg += (1/float(max_iter))*d;
	}

	return d_avg;
}

Matrix<TOTAL_VARS> diaghess_costfunc(std::vector<std::vector<Matrix<X_DIM> > >& P, std::vector<Matrix<U_DIM> >& U) {
	Matrix<TOTAL_VARS> diaghess_avg;

	int max_iter = 100;
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
					cost = costfunc_noise(P, U, dyn_noise, obs_noise, sampling_noise);
					orig = P[t][m][i];

					P[t][m][i] = orig + step;
					cost_p = costfunc_noise(P, U, dyn_noise, obs_noise, sampling_noise);

					P[t][m][i] = orig - step;
					cost_l = costfunc_noise(P,U, dyn_noise, obs_noise, sampling_noise);

					P[t][m][i] = orig;
					diaghess[index++] = (cost_p - 2*cost + cost_l)/(step*step);
				}
			}

			if (t < T-1) {
				for(int i=0; i < U_DIM; ++i) {
					cost = costfunc_noise(P, U, dyn_noise, obs_noise, sampling_noise);
					orig = U[t][i];

					U[t][i] = orig + step;
					cost_p = costfunc_noise(P,U, dyn_noise, obs_noise, sampling_noise);

					U[t][i] = orig - step;
					cost_l = costfunc_noise(P,U, dyn_noise, obs_noise, sampling_noise);

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
