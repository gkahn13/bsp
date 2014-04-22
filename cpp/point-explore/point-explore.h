#ifndef __POINT_EXPLORE_H__
#define __POINT_EXPLORE_H__

#include <vector>
#include <cmath>

#include "../util/matrix.h"
#include "../util/utils.h"
#include "../util/logging.h"

#include <Python.h>
#include <boost/python.hpp>
#include <boost/filesystem.hpp>

namespace py = boost::python;

#define TIMESTEPS 10
#define PARTICLES 100
#define DT 1.0
#define X_DIM 2
#define U_DIM 2
#define Z_DIM 1
#define Q_DIM 2
#define R_DIM 1

const int T = TIMESTEPS;
const int M = PARTICLES;
const int TOTAL_VARS = T*X_DIM + (T-1)*U_DIM;

SymmetricMatrix<R_DIM> R;

Matrix<X_DIM> x0;
Matrix<X_DIM> xMin, xMax;
Matrix<U_DIM> uMin, uMax;

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

const double step = 0.0078125*0.0078125;
const double INFTY = 1e10;

const double alpha = 2.5;
const double max_range = 0.25;

const double alpha_control_norm = 0; // 1e-5
const double alpha_control_smooth = 1e-2; // 1e-2

Matrix<X_DIM> target;

#include "casadi/casadi-point-explore.h"
namespace AD = CasADi;
AD::SXFunction casadi_differential_entropy_func;
AD::SXFunction casadi_grad_differential_entropy_func;
AD::SXFunction casadi_diaghess_differential_entropy_func;

inline float uniform(float low, float high) {
	return (high - low)*(rand() / float(RAND_MAX)) + low;
}

template<int _dim>
float dist(const Matrix<_dim>& a, const Matrix<_dim>& b) {
	return sqrt(tr(~(a-b)*(a-b)));
}

namespace point_explore {


void initialize() {
	R = 1e-2*identity<R_DIM>();

	x0[0] = 0; x0[1] = 0;

	xMin[0] = 0; xMin[1] = 0;
	xMax[0] = 5; xMax[1] = 5;
	uMin[0] = -.25; uMin[1] = -.25;
	uMax[0] = .25; uMax[1] = .25;

	target[0] = 3; target[1] = 3;

	casadi_differential_entropy_func = casadi_point_explore::casadi_differential_entropy_func();
	casadi_grad_differential_entropy_func = casadi_point_explore::casadi_differential_entropy_gradfunc();
//	casadi_diaghess_differential_entropy_func = casadi_point_explore::casadi_differential_entropy_diaghessfunc();

//	casadi_differential_entropy_func.generateCode("casadi_de_func.c");
//	casadi_grad_differential_entropy_func.generateCode("casadi_grad_de_func.c");
//	casadi_diaghess_differential_entropy_func.generateCode("casadi_diaghess_de_func.c");
}

// notice zero noise influence
Matrix<X_DIM> dynfunc(const Matrix<X_DIM>& x, const Matrix<U_DIM>& u) {
	Matrix<X_DIM> xNew = x + u*DT;
	return xNew;
}

Matrix<Z_DIM> obsfunc(const Matrix<X_DIM>& x, const Matrix<X_DIM>& t, const Matrix<R_DIM>& r) {
	Matrix<Z_DIM> z;
	z[0] = 1.0/(1.0+exp(-alpha*(max_range-dist<X_DIM>(x, t)))) + r[0];
	return z;
}

template <int _vDim, int _sDim>
float gaussLikelihood(const Matrix<_vDim>& v, const SymmetricMatrix<_sDim>& S) {
	Matrix<_sDim,_sDim> Sf;
	chol(S, Sf);
	Matrix<_sDim,1> M = (!Sf)*v;

	Matrix<_sDim> E_exp;
	float E_exp_sum = 0;
	E_exp_sum = exp(-0.5*tr(~M*M));

	float Sf_diag_prod = 1;
	for(int i=0; i < _sDim; ++i) { Sf_diag_prod *= Sf(i,i); }
	float C = pow(2*M_PI, _vDim/2)*Sf_diag_prod;

	float w = E_exp_sum / C;
	return w;
}

std::vector<Matrix<X_DIM> > lowVarianceSampler(const std::vector<Matrix<X_DIM> >& P, const std::vector<float>& W, float r) {
	std::vector<Matrix<X_DIM> > P_sampled(M);

	float c = W[0];
	int i = 0;
	for(int m=0; m < M; ++m) {
		float u = r + (m) * (1/float(M));
		while (u > c) {
			c += W[++i];
		}
		P_sampled[m] = P[i];
	}

	return P_sampled;
}


void updateStateAndParticles(const Matrix<X_DIM>& x_t, const std::vector<Matrix<X_DIM> >& P_t, const Matrix<U_DIM>& u_t,
								Matrix<X_DIM>& x_tp1, std::vector<Matrix<X_DIM> >& P_tp1) {
	x_tp1 = dynfunc(x_t, u_t);
	// receive noisy measurement
	Matrix<Z_DIM> z_tp1 = obsfunc(x_tp1, target, sampleGaussian(zeros<R_DIM,1>(), R));

	std::vector<float> W(M);
	float W_sum = 0;
	Matrix<Z_DIM> z_particle;
	// for each particle, weight by gaussLikelihood of that measurement given particle/robot distance
	for(int m=0; m < M; ++m) {
		z_particle = obsfunc(x_tp1, P_t[m], zeros<R_DIM,1>());
		Matrix<Z_DIM> e = z_particle - z_tp1;
		W[m] = gaussLikelihood<Z_DIM,Z_DIM>(e, R);
		W_sum += W[m];
	}
	for(int m=0; m < M; ++m) { W[m] = W[m] / W_sum; }

	float sampling_noise = (1/float(M))*(rand() / float(RAND_MAX));
	P_tp1 = lowVarianceSampler(P_t, W, sampling_noise);
}


float differential_entropy(const std::vector<Matrix<X_DIM> >& X, const std::vector<Matrix<U_DIM> >& U,
					 const std::vector<Matrix<X_DIM> >& P) {
	float entropy = 0;

	std::vector<Matrix<X_DIM>> X_prop(T);
	std::vector<std::vector<Matrix<Z_DIM> > > H(T, std::vector<Matrix<Z_DIM> >(M));
	for(int t=0; t < T-1; ++t) {
		X_prop[t+1] = dynfunc(X[t], U[t]);
		for(int m=0; m < M; ++m) {
			H[t+1][m] = obsfunc(X_prop[t+1], P[m], zeros<R_DIM,1>());
		}
	}

	std::vector<std::vector<float> > W(T, std::vector<float>(M,0));
	W[0] = std::vector<float>(M, 1/float(M));
	for(int t=1; t < T; ++t) {

		float W_sum = 0;
		for(int m=0; m < M; ++m) {
			for(int n=0; n < M; ++n) {
				W[t][m] += gaussLikelihood<Z_DIM, Z_DIM>(H[t][m] - H[t][n], R);
			}
			W_sum += W[t][m];
		}
		for(int m=0; m < M; ++m) { W[t][m] = W[t][m] / W_sum; }

		// use skoglar version
		float entropy_t = 0;
		for(int m=0; m < M; ++m) {
			entropy_t += -W[t][m]*log(W[t][m]);
		}

		// simplifies because zero particle dynamics
		for(int m=0; m < M; ++m) {
			entropy_t += -W[t][m]*log(W[t-1][m]);
		}

		float sum_cross_time_weights = 0;
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

double casadi_differential_entropy(const std::vector<Matrix<X_DIM> >& X, const std::vector<Matrix<U_DIM> >& U,
					 const std::vector<Matrix<X_DIM> >& P) {
	double XU_arr[T*X_DIM+(T-1)*U_DIM];
	double P_arr[T*M*X_DIM];

	casadi_point_explore::setup_casadi_vars(X, U, P, XU_arr, P_arr);

	casadi_differential_entropy_func.setInput(XU_arr,0);
	casadi_differential_entropy_func.setInput(P_arr,1);

	casadi_differential_entropy_func.evaluate();

	double cost = 0;
	casadi_differential_entropy_func.getOutput(&cost,0);

	return cost;
}

Matrix<TOTAL_VARS> grad_differential_entropy(std::vector<Matrix<X_DIM> >& X, std::vector<Matrix<U_DIM> >& U,
												const std::vector<Matrix<X_DIM> >& P) {
	Matrix<TOTAL_VARS> g;

	float orig, entropy_p, entropy_l;
	int index = 0;
	for(int t=0; t < T; ++t) {
		for(int i=0; i < X_DIM; ++i) {
			orig = X[t][i];

			X[t][i] = orig + step;
			entropy_p = differential_entropy(X, U, P);

			X[t][i] = orig - step;
			entropy_l = differential_entropy(X, U, P);

			X[t][i] = orig;
			g[index++] = (entropy_p - entropy_l)/(2*step);
		}

		if (t < T-1) {
			for(int i=0; i < U_DIM; ++i) {
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

Matrix<TOTAL_VARS> casadi_grad_differential_entropy(const std::vector<Matrix<X_DIM> >& X, const std::vector<Matrix<U_DIM> >& U,
					 const std::vector<Matrix<X_DIM> >& P) {
	double XU_arr[T*X_DIM+(T-1)*U_DIM];
	double P_arr[T*M*X_DIM];

	casadi_point_explore::setup_casadi_vars(X, U, P, XU_arr, P_arr);

	casadi_grad_differential_entropy_func.setInput(XU_arr,0);
	casadi_grad_differential_entropy_func.setInput(P_arr,1);

	casadi_grad_differential_entropy_func.evaluate();

	Matrix<TOTAL_VARS> grad;
	casadi_grad_differential_entropy_func.getOutput(grad.getPtr(),0);

	return grad;
}

Matrix<TOTAL_VARS> diaghess_differential_entropy(std::vector<Matrix<X_DIM> >& X, std::vector<Matrix<U_DIM> >& U,
												const std::vector<Matrix<X_DIM> >& P) {
	Matrix<TOTAL_VARS> g;

	float initial_entropy = differential_entropy(X, U, P);
	float orig, entropy_p, entropy_l;
	int index = 0;
	for(int t=0; t < T; ++t) {
		for(int i=0; i < X_DIM; ++i) {
			orig = X[t][i];

			X[t][i] = orig + step;
			entropy_p = differential_entropy(X, U, P);

			X[t][i] = orig - step;
			entropy_l = differential_entropy(X, U, P);

			X[t][i] = orig;
			g[index++] = (entropy_p - 2*initial_entropy + entropy_l)/(step*step);
		}

		if (t < T-1) {
			for(int i=0; i < U_DIM; ++i) {
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

Matrix<TOTAL_VARS> casadi_diaghess_differential_entropy(const std::vector<Matrix<X_DIM> >& X, const std::vector<Matrix<U_DIM> >& U,
					 const std::vector<Matrix<X_DIM> >& P) {
	double XU_arr[T*X_DIM+(T-1)*U_DIM];
	double P_arr[T*M*X_DIM];

	casadi_point_explore::setup_casadi_vars(X, U, P, XU_arr, P_arr);

	casadi_diaghess_differential_entropy_func.setInput(XU_arr,0);
	casadi_diaghess_differential_entropy_func.setInput(P_arr,1);

	casadi_diaghess_differential_entropy_func.evaluate();

	Matrix<TOTAL_VARS> diaghess;
	casadi_diaghess_differential_entropy_func.getOutput(diaghess.getPtr(),0);

	return diaghess;
}



void pythonDisplayStatesAndParticles(const std::vector<Matrix<X_DIM>>& X, const std::vector<Matrix<X_DIM> >& P, const Matrix<X_DIM>& targ) {
	py::list x_list, targ_list;
	for(int t=0; t < X.size(); ++t) {
		for(int i=0; i < X_DIM; ++i) {
			x_list.append(X[t][i]);
			targ_list.append(targ[i]);
		}
	}

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

		plot_state_and_particles(x_list, particles_list, targ_list, X_DIM, M);

		LOG_INFO("Press enter to continue");
		py::exec("raw_input()",main_namespace);
	}
	catch(py::error_already_set const &)
	{
		PyErr_Print();
	}


}

}


#endif
