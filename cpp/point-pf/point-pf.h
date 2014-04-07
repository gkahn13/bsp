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

namespace py = boost::python;

#define TIMESTEPS 3
#define DT 1.0
#define X_DIM 2
#define U_DIM 2
#define Z_DIM 2
#define Q_DIM 2
#define R_DIM 2

const int T = TIMESTEPS;

SymmetricMatrix<Q_DIM> Q;
SymmetricMatrix<R_DIM> R;

Matrix<X_DIM> dynfunc(const Matrix<X_DIM>& x, const Matrix<U_DIM>& u, const Matrix<Q_DIM>& q)
{
	Matrix<X_DIM> xNew = x + u*DT + 0.01*q;
	return xNew;
}

Matrix<Z_DIM> obsfunc(const Matrix<X_DIM>& x, const Matrix<R_DIM>& r)
{
	Matrix<Z_DIM> z = x + sqrt(0.5*0.5*x[0]*x[0] + 1e-6)*r;
	return z;
}

float gaussLikelihood(const Matrix<Z_DIM>& v, const SymmetricMatrix<R_DIM>& S) {
	Matrix<R_DIM,R_DIM> Sf;
	chol(S, Sf);
	Matrix<R_DIM,1> M = (!Sf)*v;

	Matrix<R_DIM> E_exp;
	float E_exp_sum = 0;
	E_exp_sum = exp(-0.5*tr(~M*M));
//	for(int i=0; i < R_DIM; ++i) {
//		Matrix<R_DIM> M_col = M.subMatrix<R_DIM,1>(0,0);
//		E_exp[i] = exp(-0.5 * tr(~M_col*M_col));
//		E_exp_sum += E_exp[i];
//	}

	float Sf_diag_prod = 1;
	for(int i=0; i < R_DIM; ++i) { Sf_diag_prod *= Sf(i,i); }
	float C = pow(2*M_PI, Z_DIM/2)*Sf_diag_prod;

	float w = E_exp_sum / C;
	return w;
}

std::vector<Matrix<X_DIM> > lowVarianceSampler(const std::vector<Matrix<X_DIM> >& P, const std::vector<float>& W, float r) {
	int M = P.size();
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

std::vector<Matrix<X_DIM> > beliefDynamics(const std::vector<Matrix<X_DIM> >& P_t, const Matrix<U_DIM>& u, const std::vector<Matrix<Q_DIM> > dyn_noise, const std::vector<Matrix<R_DIM> > obs_noise, float sampling_noise) {
	int M = P_t.size();
	std::vector<Matrix<X_DIM> > P_tp1_bar(M), P_tp1;
	std::vector<float> W(M);

	float W_sum = 0;
	for(int m=0; m < M; ++m) {
		P_tp1_bar[m] = dynfunc(P_t[m], u, dyn_noise[m]);
		Matrix<Z_DIM> e = obsfunc(P_tp1_bar[m], obs_noise[m]) - obsfunc(P_tp1_bar[m], zeros<R_DIM,1>());
		W[m] = gaussLikelihood(e, R);
		W_sum += W[m];
	}
	for(int m=0; m < M; ++m) { W[m] = W[m] / W_sum; }

	P_tp1 = lowVarianceSampler(P_tp1_bar, W, sampling_noise);

	return P_tp1;
}

void pythonDisplayParticles(std::vector<std::vector<Matrix<X_DIM> > >& particle_list) {
	int M = particle_list[0].size();

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
		py::exec(py::str("sys.path.append('"+workingDir+"/slam')"), main_namespace);
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


#endif
