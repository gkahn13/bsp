#ifndef __CASADI_SYSTEM_H__
#define __CASADI_SYSTEM_H__

#include <symbolic/casadi.hpp>
#include <symbolic/stl_vector_tools.hpp>
#include <cstdlib>

namespace AD = CasADi;

#include <armadillo>
using namespace arma;

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

class CasadiSystem {
public:
	CasadiSystem();

	double casadi_cost(const std::vector<mat>& X, const std::vector<mat>& U,
			const mat& P);
	mat casadi_cost_grad(const std::vector<mat>& X, const std::vector<mat>& U,
			const mat& P);

protected:
	int T, M, TOTAL_VARS, X_DIM, U_DIM, Z_DIM, Q_DIM, R_DIM;
	double DT;

	mat R;

	AD::SXFunction cost_func, cost_grad_func;

	AD::SXMatrix dist(AD::SXMatrix a, AD::SXMatrix b);
	virtual AD::SXMatrix dynfunc(const AD::SXMatrix& x_t, const AD::SXMatrix& u_t) =0;
	virtual AD::SXMatrix obsfunc(const AD::SXMatrix& x, const AD::SXMatrix& t) =0;

	AD::SXMatrix gauss_likelihood(const AD::SXMatrix& v);

	AD::SXMatrix cost_entropy(const std::vector<AD::SXMatrix>& X, const std::vector<AD::SXMatrix>& U,
										const std::vector<AD::SXMatrix>& P);
	AD::SXMatrix cost_platt(const std::vector<AD::SXMatrix>& X, const std::vector<AD::SXMatrix>& U,
									const std::vector<AD::SXMatrix>& P);

	virtual AD::SXMatrix cost_wrapper(const AD::SXMatrix& XU_vec, const AD::SXMatrix& P_vec) =0;
	AD::SXFunction casadi_cost_func();
	AD::SXFunction casadi_cost_grad_func();

	void setup_casadi_vars(const std::vector<mat>& X, const std::vector<mat>& U,
			const mat& P, double* XU_arr, double* P_arr);
};

#endif
