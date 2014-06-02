#include <Python.h>

#include <boost/python.hpp>
#include <boost/filesystem.hpp>
namespace py = boost::python;

//#include "../kdpee/kdpee/kdpee.h"

#include <../util/logging.h>

#include <iostream>

#include <armadillo>

using namespace arma;

#define TIMESTEPS 50
#define PARTICLES 100

#define X_DIM 2
#define U_DIM 2
#define Z_DIM 2
#define Q_DIM 2
#define R_DIM 2

#define DT 1.0

const int T = TIMESTEPS;
const int M = PARTICLES;

mat::fixed<X_DIM, X_DIM> A;
mat::fixed<X_DIM, U_DIM> B;
mat::fixed<Z_DIM, X_DIM> C;
mat::fixed<Z_DIM, 1> D;

mat::fixed<Q_DIM, Q_DIM> Q;
mat::fixed<R_DIM, R_DIM> R;

inline mat sample_gaussian(mat mean, mat covariance) {
	mat sample(mean.n_rows, mean.n_cols, fill::randn);

	mat L = chol(covariance);
	sample = L*sample + mean;

	return sample;
}

mat dynfunc(const mat& x, const mat& u, const mat& q) {
	mat x_new = A*x + B*u + q;
	return x_new;
}

mat obsfunc(const mat& x, const mat& r) {
	mat z = C*x + D + r;
	return z;
}

void kalman_update(const mat& mu_tp1_p, const mat& sigma_t, const mat& z_tp1, mat& mu_tp1, mat& sigma_tp1) {
	mat sigma_tp1_p = A*sigma_t*A.t() + Q;

	mat K_tp1 = sigma_tp1_p*C.t()*inv(C*sigma_tp1_p*C.t() + R);
	mu_tp1 = mu_tp1_p + K_tp1*(z_tp1 - (C*mu_tp1_p + D));
	sigma_tp1 = (eye<mat>(X_DIM,X_DIM) - K_tp1*C)*sigma_tp1_p;
}

double gauss_likelihood(const mat& v, const mat& S) {
	mat Sf = chol(S);
	mat M = solve(Sf, v);

	double E = -0.5*accu(M % M);
	double C = pow(2*M_PI, S.n_cols/2) * prod(diagvec(Sf));
	double w = exp(E) / C;

	return w;
}

mat low_variance_sampler(const mat& P, const mat& W, double r) {
	int M = P.n_cols;
	mat P_sampled(P.n_rows, P.n_cols, fill::zeros);

	double c = W(0);
	int i = 0;
	for(int m=0; m < M; ++m) {
		double u = r + m * (1/double(M));
		while (u > c) {
			c += W(++i);
		}
		P_sampled.col(m) = P.col(i) + sample_gaussian(zeros<mat>(X_DIM,1), .1*eye<mat>(X_DIM,X_DIM));
	}

	return P_sampled;
}

void particle_update(const mat& P_tp1_p, const mat& z_tp1, mat& P_tp1) {
	int M = P_tp1_p.n_cols;

	mat W(M, 1, fill::zeros);
	mat r(R_DIM, 1, fill::zeros);
	for(int m=0; m < M; ++m) {
		mat z_particle = obsfunc(P_tp1_p.col(m), r);
		mat e = z_particle - z_tp1;
		W(m) = gauss_likelihood(e, R);
	}
	W = W / accu(W);

	double sampling_noise = (1/double(M))*(rand() / RAND_MAX);
	P_tp1 = low_variance_sampler(P_tp1_p, W, sampling_noise);
}

double particle_entropy(const mat& P) {
	int M = P.n_cols;

	double entropy = 0;

	mat H(Z_DIM, M);
	mat r(R_DIM, 1, fill::zeros);
	for(int m=0; m < M; ++m) {
		H.col(m) = obsfunc(P.col(m), r);
	}

	mat W_tm1 = (1/float(M))*ones<mat>(M, 1);
	mat W(M, 1, fill::zeros);
	mat U_inv = inv(chol(R));
	mat I = eye<mat>(R_DIM,R_DIM);
	for(int m=0; m < M; ++m) {
		for(int p=0; p < M; ++p) {
			mat diff = H.col(m) - H.col(p);
//			W(m) += gauss_likelihood(U_inv*diff,I);;
			W(m) += gauss_likelihood(diff, 2*R + I);
		}
	}
	W = W / accu(W);

	mat dyn_prob(M, M, fill::zeros);
	for(int m=0; m < M; ++m) {
		for(int p=0; p < M; ++p) {
			mat diff = P.col(m) - P.col(p);
			dyn_prob(m, p) = gauss_likelihood(diff, Q);
		}
	}

	entropy += accu(-W % log(W));

	for(int m=0; m < M; ++m) {
		double weighted_dyn_prob = 0;
		for(int p=0; p < M; ++p) {
			weighted_dyn_prob += W_tm1(m)*dyn_prob(m, p);
		}
		entropy += -W(m)*log(weighted_dyn_prob);
	}
//	entropy += accu(-W % log(W_tm1)); // zero dynamics approximation

	entropy += log(accu(W_tm1 % W));

	return entropy;
}

//floatval particle_entropy_kdpee(const mat& P) {
//	const int d = P.n_rows;
//	const int n = P.n_cols;
//
//	floatval **dimrefs = new floatval*[d];
//	for(int i=0; i < d; ++i) {
//		dimrefs[i] = new floatval[n];
//		for(int j=0; j < n; ++j) {
//			dimrefs[i][j] = P(i,j);
//		}
//	}
//
//	floatval *mins = new floatval[d];
//	floatval *maxs = new floatval[d];
//
//	for(int i=0; i < d; ++i) {
//		mat row = P.row(i);
//		mins[i] = row.min();
//		maxs[i] = row.max();
//	}
//
//	const floatval zcut = 1.96;
//
//	int *keys = new int[n];
//	for(int i=0; i < n; ++i) {
//		keys[i] = i+1;
//	}
//
//	const floatval **dimrefs_const = const_cast<const floatval**>(dimrefs);
//	floatval entropy_kdpee = kdpee(dimrefs_const, n, d, mins, maxs, zcut, keys);
//	return entropy_kdpee;
//}

void plot(const mat& mu, const mat& sigma, const mat& P) {
	int M = P.n_cols;

	py::list mu_list;
	for(int i=0; i < X_DIM; ++i) {
		mu_list.append(mu(i));
	}

	py::list sigma_list;
	for(int i=0; i < X_DIM; ++i) {
		for(int j=0; j < X_DIM; ++j) {
			sigma_list.append(sigma(i,j));
		}
	}

	py::list particles_list;
	for(int i=0; i < X_DIM; ++i) {
		for(int m=0; m < M; ++m) {
			particles_list.append(P(i,m));
		}
	}

	std::string working_dir = boost::filesystem::current_path().normalize().string();
	std::string bsp_dir = working_dir.substr(0,working_dir.find("bsp"));
	std::string explore_dir = bsp_dir + "bsp/pf/tests";

	try
	{
		Py_Initialize();
		py::object main_module = py::import("__main__");
		py::object main_namespace = main_module.attr("__dict__");
		py::exec("import sys, os", main_namespace);
		py::exec(py::str("sys.path.append('"+working_dir+"')"), main_namespace);
		py::object plot_module = py::import("plot_test_entropy");
		py::object plot = plot_module.attr("plot");

		plot(mu_list, sigma_list, particles_list, X_DIM, M);

		std::cout << "Press enter to continue\n";
		py::exec("raw_input()",main_namespace);
	}
	catch(py::error_already_set const &)
	{
		PyErr_Print();
	}
}

void test_entropy() {
	srand(time(0));

	A = eye<mat>(X_DIM, X_DIM);
	B = zeros<mat>(X_DIM, X_DIM);
	C = eye<mat>(Z_DIM, X_DIM);
	D = zeros<mat>(Z_DIM, 1);

	Q = 1*eye<mat>(Q_DIM, Q_DIM); // 1e-2
	R = 1*eye<mat>(R_DIM, R_DIM); // 1e-1

	int max_iter = 10;
	mat mean_err_entropy(1, 1, fill::zeros), std_err_entropy(1, 1, fill::zeros);
	mat mean_err_kdpee(1, 1, fill::zeros), std_err_kdpee(1, 1, fill::zeros);
	for(int iter=0; iter < max_iter; ++iter) {
		std::cout << "iter: " << iter << "\n";
		mat x0(X_DIM, 1, fill::zeros);
		mat sigma0 = 1*eye<mat>(X_DIM, X_DIM);

		mat P0(X_DIM, M, fill::zeros);
		for(int m=0; m < M; ++m) {
			P0.col(m) = sample_gaussian(x0, sigma0);
		}

		mat x_t = x0, sigma_t = sigma0;
		mat P_t = P0;

		mat x_tp1(X_DIM, 1, fill::zeros), sigma_tp1(X_DIM, X_DIM, fill::zeros);
		mat P_tp1(X_DIM, M, fill::zeros);

		mat k_entropys(T, 1, fill::zeros), p_entropys(T, 1, fill::zeros), kdpee_entropys(T, 1, fill::zeros);

		k_entropys(0) = 0.5*log(pow(2*M_PI*exp(1), X_DIM) * det(sigma0));
//		std::cout << "sigma0 entropy: " << k_entropys(0) << "\n";

		p_entropys(0) = particle_entropy(P0);
//		std::cout << "P0 entropy: " << p_entropys(0) << "\n";

//		kdpee_entropys(0) = particle_entropy_kdpee(P0);
//		std::cout << "kdpee0 entropy: " << kdpee_entropys(0) << "\n";

		for(int t=1; t < T; ++t) {
//			std::cout << "\nt: " << t << "\n";
			mat u = sample_gaussian(zeros<mat>(U_DIM,1), 1*eye<mat>(U_DIM,U_DIM));
			mat dyn_noise = sample_gaussian(zeros<mat>(Q_DIM,1), .1*Q);
			mat obs_noise = sample_gaussian(zeros<mat>(R_DIM,1), .1*R);

			// kalman update
			mat x_tp1_p = dynfunc(x_t, u, dyn_noise);
			mat z_tp1 = obsfunc(x_tp1_p, obs_noise);
			kalman_update(x_tp1_p, sigma_t, z_tp1, x_tp1, sigma_tp1);
			k_entropys(t) = 0.5*log(pow(2*M_PI*exp(1), X_DIM) * det(sigma_tp1));

//			std::cout << "kalman entropy: " << k_entropys(t) << "\n";
			x_t = x_tp1;
			sigma_t = sigma_tp1;

			// particle update
			mat P_tp1_p(X_DIM, M, fill::zeros);
			for(int m=0; m < M; ++m) {
				P_tp1_p.col(m) = dynfunc(P_t.col(m), u, dyn_noise);
			}
			particle_update(P_tp1_p, z_tp1, P_tp1);
			p_entropys(t) = particle_entropy(P_tp1);
//			std::cout << "particle entropy: " << p_entropys(t) << "\n";

//			kdpee_entropys(t) = particle_entropy_kdpee(P_tp1);
//			std::cout << "kdpee entropy: " << kdpee_entropys(t) << "\n";

			P_t = P_tp1;

			//		plot(x_t, sigma_t, P_t);
		}

		mat percent_diff_entropy = 100*abs((k_entropys - p_entropys)/p_entropys);
		mean_err_entropy += (1/double(max_iter))*mean(percent_diff_entropy);
		std_err_entropy += (1/double(max_iter))*stddev(percent_diff_entropy);

//		mat percent_diff_kdpee = 100*abs((k_entropys - kdpee_entropys)/kdpee_entropys);
//		mean_err_kdpee += (1/double(max_iter))*mean(percent_diff_kdpee);
//		std_err_kdpee+= (1/double(max_iter))*stddev(percent_diff_kdpee);
	}

	std::cout << "\nparticle total mean rel err: " << mean_err_entropy(0) << "\n";
	std::cout << "particle total std rel err: " << std_err_entropy(0) << "\n\n";

//	std::cout << "kdpee total mean rel err: " << mean_err_kdpee(0) << "\n";
//	std::cout << "kdpee total std rel err: " << std_err_kdpee(0) << "\n\n";
}

int main(int argc, char* argv[]) {
	test_entropy();
}
