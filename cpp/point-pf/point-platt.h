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

	// TODO: if no propagation, then all particle gradients = 0
	std::vector<std::vector<Matrix<X_DIM> > > P_prop(T, std::vector<Matrix<X_DIM> >(M));
	P_prop[0] = P[0];
	for(int t=0; t < T-1; ++t) {
		for(int m=0; m < M; ++m) {
			P_prop[t+1][m] = point_pf::dynfunc(P[t][m], U[t], dyn_noise[t*M+m]);
		}
	}

	Matrix<X_DIM> P_t_mle, P_t_other;
	Matrix<Z_DIM> h_t_mle, h_t_other;
	Matrix<T*Z_DIM> h_diff;
	for(int m=1; m < M; ++m) {
		for(int t=0; t < T; ++t) {
			P_t_mle = P_prop[t][0];
			h_t_mle = point_pf::obsfunc(P_t_mle, obs_noise[t]);

			P_t_other = P_prop[t][m];
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

	double PU_arr[T*M*X_DIM+(T-1)*U_DIM];
	double dyn_noise_arr[(T-1)*M*Q_DIM];
	double obs_noise_arr[M*T*R_DIM];

	std::vector<Matrix<Q_DIM> > dyn_noise = sampleGaussianN(zeros<Q_DIM,1>(), Q, (T-1)*M);
	std::vector<Matrix<Q_DIM> > obs_noise = sampleGaussianN(zeros<Q_DIM,1>(), Q, T*M);

	casadi_platt::setup_casadi_vars(P, U, dyn_noise, obs_noise, PU_arr, dyn_noise_arr, obs_noise_arr);


	double real_cost = costfunc_noise(P,U,dyn_noise,obs_noise,std::vector<float>(T));
	std::cout << "Actual cost: " << real_cost << "\n";

	AD::SXFunction func = casadi_platt::casadi_costfunc();

	func.setInput(PU_arr,0);
	func.setInput(dyn_noise_arr,1);
	func.setInput(obs_noise_arr,2);

	func.evaluate();

	double cost = 0;
	func.getOutput(&cost,0);

	std::cout << "Casadi cost: " << cost << "\n";

	return cost;
}

Matrix<TOTAL_VARS> grad_costfunc(std::vector<std::vector<Matrix<X_DIM> > >& P, std::vector<Matrix<U_DIM> >& U) {
	Matrix<TOTAL_VARS> g_avg;

	int max_iter = 100;
	for(int iter=0; iter < max_iter; iter++) {
		Matrix<TOTAL_VARS> g;

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
					g[index++] = (cost_p - cost_l)/(2*step);
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
					g[index++] = (cost_p - cost_l)/(2*step);
				}
			}
		}

		g_avg += (1/float(max_iter))*g;
	}

	return g_avg;
}

Matrix<TOTAL_VARS> casadi_grad_costfunc(std::vector<std::vector<Matrix<X_DIM> > >& P, std::vector<Matrix<U_DIM> >& U) {
	Matrix<TOTAL_VARS> g_avg, g;
	double cost;

	double PU_arr[T*M*X_DIM+(T-1)*U_DIM];
	double dyn_noise_arr[(T-1)*M*Q_DIM];
	double obs_noise_arr[M*T*R_DIM];

	AD::SXFunction costgrad_func = casadi_platt::casadi_costgradfunc();

	int max_iter = 1e3;
	for(int iter=0; iter < max_iter; iter++) {
		std::vector<Matrix<Q_DIM> > dyn_noise = sampleGaussianN(zeros<Q_DIM,1>(), Q, (T-1)*M);
		std::vector<Matrix<Q_DIM> > obs_noise = sampleGaussianN(zeros<Q_DIM,1>(), Q, T*M);

		casadi_platt::setup_casadi_vars(P, U, dyn_noise, obs_noise, PU_arr, dyn_noise_arr, obs_noise_arr);

		costgrad_func.setInput(PU_arr,0);
		costgrad_func.setInput(dyn_noise_arr,1);
		costgrad_func.setInput(obs_noise_arr,2);

		costgrad_func.evaluate();

		costgrad_func.getOutput(&cost,0);
		costgrad_func.getOutput(g.getPtr(),1);

		g_avg += (1/float(max_iter))*g;
	}

	return g_avg;
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

Matrix<TOTAL_VARS> casadi_diaghess_costfunc(std::vector<std::vector<Matrix<X_DIM> > >& P, std::vector<Matrix<U_DIM> >& U) {
	Matrix<TOTAL_VARS> diaghess_avg, diaghess;
	double cost;

	double PU_arr[T*M*X_DIM+(T-1)*U_DIM];
	double dyn_noise_arr[(T-1)*M*Q_DIM];
	double obs_noise_arr[M*T*R_DIM];

	AD::SXFunction costdiaghess_func = casadi_platt::casadi_costdiaghessfunc();

	int max_iter = 1e3;
	for(int iter=0; iter < max_iter; iter++) {
		std::vector<Matrix<Q_DIM> > dyn_noise = sampleGaussianN(zeros<Q_DIM,1>(), Q, (T-1)*M);
		std::vector<Matrix<Q_DIM> > obs_noise = sampleGaussianN(zeros<Q_DIM,1>(), Q, T*M);

		casadi_platt::setup_casadi_vars(P, U, dyn_noise, obs_noise, PU_arr, dyn_noise_arr, obs_noise_arr);

		costdiaghess_func.setInput(PU_arr,0);
		costdiaghess_func.setInput(dyn_noise_arr,1);
		costdiaghess_func.setInput(obs_noise_arr,2);

		costdiaghess_func.evaluate();

		costdiaghess_func.getOutput(&cost,0);
		costdiaghess_func.getOutput(diaghess.getPtr(),1);

		diaghess_avg += (1/float(max_iter))*diaghess;
	}

	std::cout << "diaghess_avg\n" << diaghess_avg;

	return diaghess_avg;
}


}

#endif
