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

#define TIMESTEPS 15
#define PARTICLES 100
#define DT 1.0
#define X_DIM 2
#define U_DIM 2
#define Z_DIM 1
#define Q_DIM 2
#define R_DIM 1

const int T = TIMESTEPS;
const int M = PARTICLES;

SymmetricMatrix<R_DIM> R;

Matrix<X_DIM> x0;
Matrix<X_DIM> xMin, xMax;

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

const double step = 0.0078125*0.0078125;
const double INFTY = 1e10;

const double alpha = 10;
const double max_range = 0.25;

Matrix<X_DIM> target;

namespace point_explore {

template<int _dim>
float dist(const Matrix<_dim>& a, const Matrix<_dim>& b) {
	return sqrt(tr(~(a-b)*(a-b)));
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

void pythonDisplayStateAndParticles(const Matrix<X_DIM>& x, const std::vector<Matrix<X_DIM> >& P, const Matrix<X_DIM>& targ) {
	py::list x_list, targ_list;
	for(int i=0; i < X_DIM; ++i) {
		x_list.append(x[i]);
		targ_list.append(targ[i]);
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
