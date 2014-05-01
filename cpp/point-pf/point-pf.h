#ifndef __POINT_PF_H__
#define __POINT_PF_H__

#include <vector>
#include <cmath>

#include "../util/matrix.h"
#include "../util/utils.h"
#include "../util/logging.h"

#include <Python.h>
#include <boost/python.hpp>
#include <boost/filesystem.hpp>

#include <symbolic/casadi.hpp>
#include <symbolic/stl_vector_tools.hpp>

namespace py = boost::python;
namespace AD = CasADi;

#define TIMESTEPS 15
#define PARTICLES 5
#define DT 1.0
#define X_DIM 2
#define U_DIM 2
#define Z_DIM 2
#define Q_DIM 2
#define R_DIM 2

const int T = TIMESTEPS;
const int M = PARTICLES;
const int TOTAL_VARS = T*M*X_DIM + (T-1)*U_DIM;

SymmetricMatrix<Q_DIM> Q;
SymmetricMatrix<R_DIM> R;

Matrix<X_DIM> x0, xGoal;
Matrix<X_DIM> xMin, xMax;
Matrix<U_DIM> uMin, uMax;

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

const double step = 0.0078125*0.0078125;
const double INFTY = 1e10;

// .001 for platt
// .1 for entropy
const double alpha_control = .1;

#include "casadi/casadi-point-pf.h"
AD::SXFunction casadi_belief_dynamics_func;

namespace point_pf {

void initialize() {
	Q = 1e-2*identity<Q_DIM>();
	R = 1e-2*identity<R_DIM>(); // not sure

	x0[0] = -3.5; x0[1] = 2;
	xGoal[0] = -3.5; xGoal[1] = -2;

	xMin[0] = -5; xMin[1] = -3;
	xMax[0] = 5; xMax[1] = 3;
	uMin[0] = -1; uMin[1] = -1;
	uMax[0] = 1; uMax[1] = 1;

	casadi_belief_dynamics_func = casadi_point_pf::casadiBeliefDynamicsFunc(M);
}


Matrix<X_DIM> dynfunc(const Matrix<X_DIM>& x, const Matrix<U_DIM>& u, const Matrix<Q_DIM>& q)
{
	Matrix<X_DIM> xNew = x + u*DT + .01*q;
	return xNew;
}

Matrix<Z_DIM> obsfunc(const Matrix<X_DIM>& x, const Matrix<R_DIM>& r)
{
	Matrix<Z_DIM> z = x + sqrt(0.5*0.5*x[0]*x[0] + 1e-6)*r;
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
//	for(int i=0; i < _sDim; ++i) {
//		Matrix<_sDim> M_col = M.subMatrix<_sDim,1>(0,0);
//		E_exp[i] = exp(-0.5 * tr(~M_col*M_col));
//		E_exp_sum += E_exp[i];
//	}

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

// other possible sampler
// http://www.stat.columbia.edu/~liam/teaching/neurostat-spr12/papers/EM/resampling.pdf



std::vector<Matrix<X_DIM> > beliefDynamics(const std::vector<Matrix<X_DIM> >& P_t, const Matrix<U_DIM>& u,
											 const std::vector<Matrix<Q_DIM> >& dyn_noise, const std::vector<Matrix<R_DIM> >& obs_noise, float sampling_noise) {
	std::vector<Matrix<X_DIM> > P_tp1_bar(M), P_tp1;
	std::vector<float> W(M);

	float W_sum = 0;
	for(int m=0; m < M; ++m) {
		P_tp1_bar[m] = dynfunc(P_t[m], u, dyn_noise[m]);
		Matrix<Z_DIM> e = obsfunc(P_tp1_bar[m], obs_noise[m]) - obsfunc(P_tp1_bar[m], zeros<R_DIM,1>());
		W[m] = gaussLikelihood<Z_DIM, R_DIM>(e, R);
		W_sum += W[m];
	}
	for(int m=0; m < M; ++m) { W[m] = W[m] / W_sum; }

	P_tp1 = lowVarianceSampler(P_tp1_bar, W, sampling_noise);

	return P_tp1;
}

std::vector<Matrix<X_DIM> > casadiBeliefDynamics(const std::vector<Matrix<X_DIM> >& P_t, const Matrix<U_DIM>& u,
													const std::vector<Matrix<Q_DIM> >& dyn_noise, const std::vector<Matrix<R_DIM> >& obs_noise, float sampling_noise) {
	double P_t_u_arr[M*X_DIM+U_DIM];
	double dyn_noise_arr[M*X_DIM];
	double obs_noise_arr[M*X_DIM];
	double sampling_noise_arr[1];

	int index = 0;
	for(int m=0; m < M; ++m) {
		for(int i=0; i < X_DIM; ++i) {
			P_t_u_arr[index] = P_t[m][i];
			dyn_noise_arr[index] = dyn_noise[m][i];
			obs_noise_arr[index] = obs_noise[m][i];
			index++;
		}
	}
	for(int i=0; i < U_DIM; ++i) { P_t_u_arr[index++] = u[i]; }

	casadi_belief_dynamics_func.setInput(P_t_u_arr,0);
	casadi_belief_dynamics_func.setInput(dyn_noise_arr,1);
	casadi_belief_dynamics_func.setInput(obs_noise_arr,2);
	casadi_belief_dynamics_func.setInput(sampling_noise_arr,3);

	casadi_belief_dynamics_func.evaluate();

	double P_tp1_arr[M*X_DIM];
	casadi_belief_dynamics_func.getOutput(P_tp1_arr,0);

	std::vector<Matrix<X_DIM> > P_tp1(M);
	index = 0;
	for(int m=0; m < M; ++m) {
		for(int i=0; i < X_DIM; ++i) {
			P_tp1[m][i] = P_tp1_arr[index++];
		}
	}

	return P_tp1;
}

void pythonDisplayParticles(std::vector<std::vector<Matrix<X_DIM> > >& particle_list) {
	py::list py_particle_list;
	for(int t=0; t < particle_list.size(); ++t) {
		for(int i=0; i < X_DIM; ++i) {
			for(int m=0; m < M; ++m) {
				py_particle_list.append(particle_list[t][m][i]);
			}
		}
	}

	std::string workingDir = boost::filesystem::current_path().normalize().string();

	try
	{
		Py_Initialize();
		py::object main_module = py::import("__main__");
		py::object main_namespace = main_module.attr("__dict__");
		py::exec("import sys, os", main_namespace);
		py::exec(py::str("sys.path.append('"+workingDir+"')"), main_namespace);
		py::object plot_module = py::import("plot_point_pf");
		py::object plot_particles = plot_module.attr("plot_particles");

		plot_particles(py_particle_list, T, X_DIM, M);

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
