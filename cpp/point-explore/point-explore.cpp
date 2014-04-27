#include "point-explore.h"
#include "casadi/casadi-point-explore.h"

SymmetricMatrix<N*R_DIM> R;

Matrix<N*X_DIM> x0;
Matrix<X_DIM> xMin, xMax;
Matrix<U_DIM> uMin, uMax;

Matrix<X_DIM> target;

AD::SXFunction casadi_differential_entropy_func;
AD::SXFunction casadi_grad_differential_entropy_func;
AD::SXFunction casadi_diaghess_differential_entropy_func;

namespace point_explore {

void initialize() {
	// 1e-2 for distance
	// 1e-3 for angle
	R = 1e-2*identity<N*R_DIM>();

	for(int i=0; i < N*X_DIM; ++i) {
		x0[i] = 0;
	}

	xMin[0] = -1; xMin[1] = -1;
	xMax[0] = 6; xMax[1] = 6;
	uMin[0] = -.25; uMin[1] = -.25;
	uMax[0] = .25; uMax[1] = .25;

	target[0] = 3; target[1] = 3;

	casadi_differential_entropy_func = casadi_point_explore::casadi_differential_entropy_func();
	casadi_grad_differential_entropy_func = casadi_point_explore::casadi_differential_entropy_gradfunc();
//	casadi_diaghess_differential_entropy_func = casadi_point_explore::casadi_differential_entropy_diaghessfunc();

}

// notice zero noise influence
Matrix<N*X_DIM> dynfunc(const Matrix<N*X_DIM>& x, const Matrix<N*U_DIM>& u) {
	Matrix<N*X_DIM> xNew = x + u*DT;
	return xNew;
}

Matrix<N*Z_DIM> obsfunc_dist(const Matrix<N*X_DIM>& x, const Matrix<X_DIM>& t, const Matrix<N*R_DIM>& r) {
	Matrix<N*Z_DIM> z;
	for(int n=0; n < N; ++n) {
		Matrix<X_DIM> x_n = x.subMatrix<X_DIM,1>(n*X_DIM,0);
		z[n] = 1.0/(1.0+exp(-alpha*(max_range-dist<X_DIM>(x_n, t)))) + r[n];
	}
	return z;
}

Matrix<N*Z_DIM> obsfunc_angle(const Matrix<N*X_DIM>& x, const Matrix<X_DIM>& t, const Matrix<N*R_DIM>& r) {
	Matrix<N*Z_DIM> z;
	for(int n=0; n < N; ++n) {
		Matrix<X_DIM> x_n = x.subMatrix<X_DIM,1>(n*X_DIM,0);
		z[n] = atan((x_n[1]-t[1])/(x_n[0]-t[0])) + r[n];
	}
	return z;
}

Matrix<N*Z_DIM> obsfunc(const Matrix<N*X_DIM>& x, const Matrix<X_DIM>& t, const Matrix<N*R_DIM>& r) {
//	return obsfunc_dist(x, t, r);
	return obsfunc_angle(x, t, r);
}


template <int _vDim, int _sDim>
double gaussLikelihood(const Matrix<_vDim>& v, const SymmetricMatrix<_sDim>& S) {
	Matrix<_sDim,_sDim> Sf;
	chol(S, Sf);
	Matrix<_sDim,1> M = (!Sf)*v;

	Matrix<_sDim> E_exp;
	double E_exp_sum = 0;
	E_exp_sum = exp(-0.5*tr(~M*M));

	double Sf_diag_prod = 1;
	for(int i=0; i < _sDim; ++i) { Sf_diag_prod *= Sf(i,i); }
	double C = pow(2*M_PI, _vDim/2)*Sf_diag_prod;

	double w = E_exp_sum / C;
	return w;
}

std::vector<Matrix<X_DIM> > lowVarianceSampler(const std::vector<Matrix<X_DIM> >& P, const std::vector<double>& W, double r) {
	int M = P.size(); // in case use high-res particle set
	std::vector<Matrix<X_DIM> > P_sampled(M);

	double c = W[0];
	int i = 0;
	for(int m=0; m < M; ++m) {
		double u = r + (m) * (1/double(M));
		while (u > c) {
			c += W[++i];
		}
		P_sampled[m] = P[i];
	}

	return P_sampled;
}


void updateStateAndParticles(const Matrix<N*X_DIM>& x_t, const std::vector<Matrix<X_DIM> >& P_t, const Matrix<N*U_DIM>& u_t,
								Matrix<N*X_DIM>& x_tp1, std::vector<Matrix<X_DIM> >& P_tp1) {
	int M = P_t.size(); // in case use high-res particle set
	x_tp1 = dynfunc(x_t, u_t);
	// receive noisy measurement
	Matrix<N*Z_DIM> z_tp1 = obsfunc(x_tp1, target, sampleGaussian(zeros<N*R_DIM,1>(), R));

	std::vector<double> W(M);
	double W_sum = 0;
	Matrix<N*Z_DIM> z_particle;
	// for each particle, weight by gaussLikelihood of that measurement given particle/robot distance
	for(int m=0; m < M; ++m) {
		z_particle = obsfunc(x_tp1, P_t[m], zeros<N*R_DIM,1>());
		Matrix<N*Z_DIM> e = z_particle - z_tp1;
		W[m] = gaussLikelihood<N*Z_DIM,N*Z_DIM>(e, R);
		W_sum += W[m];
	}
	for(int m=0; m < M; ++m) { W[m] = W[m] / W_sum; }

	double sampling_noise = (1/double(M))*(rand() / double(RAND_MAX));
	P_tp1 = lowVarianceSampler(P_t, W, sampling_noise);
}


double differential_entropy(const std::vector<Matrix<N*X_DIM> >& X, const std::vector<Matrix<N*U_DIM> >& U,
					 const std::vector<Matrix<X_DIM> >& P) {
	double entropy = 0;

	std::vector<Matrix<N*X_DIM>> X_prop(T);
	std::vector<std::vector<Matrix<N*Z_DIM> > > H(T, std::vector<Matrix<N*Z_DIM> >(M));
	for(int t=0; t < T-1; ++t) {
		X_prop[t+1] = dynfunc(X[t], U[t]);
		for(int m=0; m < M; ++m) {
			H[t+1][m] = obsfunc(X_prop[t+1], P[m], zeros<N*R_DIM,1>());
		}
	}

	std::vector<std::vector<double> > W(T, std::vector<double>(M,0));
	W[0] = std::vector<double>(M, 1/double(M));
	for(int t=1; t < T; ++t) {

		double W_sum = 0;
		for(int m=0; m < M; ++m) {
			for(int p=0; p < M; ++p) {
				W[t][m] += gaussLikelihood<N*Z_DIM, N*Z_DIM>(H[t][m] - H[t][p], R); // TODO: look into
			}
			W_sum += W[t][m];
		}
		for(int m=0; m < M; ++m) { W[t][m] = W[t][m] / W_sum; }

		// use skoglar version
		double entropy_t = 0;
		for(int m=0; m < M; ++m) {
			entropy_t += -W[t][m]*log(W[t][m]);
		}

		// simplifies because zero particle dynamics
		for(int m=0; m < M; ++m) {
			entropy_t += -W[t][m]*log(W[t-1][m]);
		}

		double sum_cross_time_weights = 0;
		for(int m=0; m < M; ++m) {
			sum_cross_time_weights += W[t-1][m]*W[t][m];
		}
		entropy_t += log(sum_cross_time_weights);

		entropy += entropy_t;

	}

	for(int t=0; t < T-2; ++t) {
		entropy += alpha_control_smooth*tr(~(U[t+1]-U[t])*(U[t+1]-U[t]));
	}

	for(int t=0; t < T-1; ++t) {
		entropy += alpha_control_norm*tr(~U[t]*U[t]);
	}

	return entropy;
}

double casadi_differential_entropy(const std::vector<Matrix<N*X_DIM> >& X, const std::vector<Matrix<N*U_DIM> >& U,
					 const std::vector<Matrix<X_DIM> >& P) {
	double XU_arr[T*N*X_DIM+(T-1)*N*U_DIM];
	double P_arr[T*M*X_DIM];

	casadi_point_explore::setup_casadi_vars(X, U, P, XU_arr, P_arr);

	casadi_differential_entropy_func.setInput(XU_arr,0);
	casadi_differential_entropy_func.setInput(P_arr,1);

	casadi_differential_entropy_func.evaluate();

	double cost = 0;
	casadi_differential_entropy_func.getOutput(&cost,0);

	return cost;
}

Matrix<TOTAL_VARS> grad_differential_entropy(std::vector<Matrix<N*X_DIM> >& X, std::vector<Matrix<N*U_DIM> >& U,
												const std::vector<Matrix<X_DIM> >& P) {
	Matrix<TOTAL_VARS> g;

	double orig, entropy_p, entropy_l;
	int index = 0;
	for(int t=0; t < T; ++t) {
		for(int i=0; i < N*X_DIM; ++i) {
			orig = X[t][i];

			X[t][i] = orig + step;
			entropy_p = differential_entropy(X, U, P);

			X[t][i] = orig - step;
			entropy_l = differential_entropy(X, U, P);

			X[t][i] = orig;
			g[index++] = (entropy_p - entropy_l)/(2*step);
		}

		if (t < T-1) {
			for(int i=0; i < N*U_DIM; ++i) {
				orig = U[t][i];

				U[t][i] = orig + step;
				entropy_p = differential_entropy(X, U, P);

				U[t][i] = orig - step;
				entropy_l = differential_entropy(X, U, P);

				U[t][i] = orig;
				g[index++] = (entropy_p - entropy_l)/(2*step);
			}
		}

	}

	return g;
}

Matrix<TOTAL_VARS> casadi_grad_differential_entropy(const std::vector<Matrix<N*X_DIM> >& X, const std::vector<Matrix<N*U_DIM> >& U,
					 const std::vector<Matrix<X_DIM> >& P) {
	double XU_arr[T*N*X_DIM+(T-1)*N*U_DIM];
	double P_arr[T*M*X_DIM];

	casadi_point_explore::setup_casadi_vars(X, U, P, XU_arr, P_arr);

	casadi_grad_differential_entropy_func.setInput(XU_arr,0);
	casadi_grad_differential_entropy_func.setInput(P_arr,1);

	casadi_grad_differential_entropy_func.evaluate();

	Matrix<TOTAL_VARS> grad;
	casadi_grad_differential_entropy_func.getOutput(grad.getPtr(),0);

	return grad;
}

Matrix<TOTAL_VARS> diaghess_differential_entropy(std::vector<Matrix<N*X_DIM> >& X, std::vector<Matrix<N*U_DIM> >& U,
												const std::vector<Matrix<X_DIM> >& P) {
	Matrix<TOTAL_VARS> g;

	double initial_entropy = differential_entropy(X, U, P);
	double orig, entropy_p, entropy_l;
	int index = 0;
	for(int t=0; t < T; ++t) {
		for(int i=0; i < N*X_DIM; ++i) {
			orig = X[t][i];

			X[t][i] = orig + step;
			entropy_p = differential_entropy(X, U, P);

			X[t][i] = orig - step;
			entropy_l = differential_entropy(X, U, P);

			X[t][i] = orig;
			g[index++] = (entropy_p - 2*initial_entropy + entropy_l)/(step*step);
		}

		if (t < T-1) {
			for(int i=0; i < N*U_DIM; ++i) {
				orig = U[t][i];

				U[t][i] = orig + step;
				entropy_p = differential_entropy(X, U, P);

				U[t][i] = orig - step;
				entropy_l = differential_entropy(X, U, P);

				U[t][i] = orig;
				g[index++] = (entropy_p - 2*initial_entropy + entropy_l)/(step*step);
			}
		}

	}

	return g;
}

Matrix<TOTAL_VARS> casadi_diaghess_differential_entropy(const std::vector<Matrix<N*X_DIM> >& X, const std::vector<Matrix<N*U_DIM> >& U,
					 const std::vector<Matrix<X_DIM> >& P) {
	double XU_arr[T*N*X_DIM+(T-1)*N*U_DIM];
	double P_arr[T*M*X_DIM];

	casadi_point_explore::setup_casadi_vars(X, U, P, XU_arr, P_arr);

	casadi_diaghess_differential_entropy_func.setInput(XU_arr,0);
	casadi_diaghess_differential_entropy_func.setInput(P_arr,1);

	casadi_diaghess_differential_entropy_func.evaluate();

	Matrix<TOTAL_VARS> diaghess;
	casadi_diaghess_differential_entropy_func.getOutput(diaghess.getPtr(),0);

	return diaghess;
}



void pythonDisplayStatesAndParticles(const std::vector<Matrix<N*X_DIM>>& X, const std::vector<Matrix<X_DIM> >& P, const Matrix<X_DIM>& targ) {
	int M = P.size(); // in case use high-res particle set

	py::list x_list;
	for(int i=0; i < N*X_DIM; ++i) {
		for(int t=0; t < X.size(); ++t) {
			x_list.append(X[t][i]);
		}
	}

	py::list targ_list;
	targ_list.append(targ[0]);
	targ_list.append(targ[1]);

	py::list particles_list;
	for(int i=0; i < X_DIM; ++i) {
		for(int m=0; m < M; ++m) {
			particles_list.append(P[m][i]);
		}
	}

	std::string workingDir = boost::filesystem::current_path().normalize().string();

	try
	{
		Py_Initialize();
		py::object main_module = py::import("__main__");
		py::object main_namespace = main_module.attr("__dict__");
		py::exec("import sys, os", main_namespace);
		py::exec(py::str("sys.path.append('"+workingDir+"/slam')"), main_namespace);
		py::object plot_module = py::import("plot_point_explore");
		py::object plot_state_and_particles = plot_module.attr("plot_state_and_particles");

		plot_state_and_particles(x_list, particles_list, targ_list, X_DIM, M, N);

		LOG_INFO("Press enter to continue");
		py::exec("raw_input()",main_namespace);
	}
	catch(py::error_already_set const &)
	{
		PyErr_Print();
	}


}

}
