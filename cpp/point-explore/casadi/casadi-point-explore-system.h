#ifndef __CASADI_POINT_EXPLORE_SYSTEM_H__
#define __CASADI_POINT_EXPLORE_SYSTEM_H__

class CasadiPointExploreSystem;
#include "../point-explore-system.h"

#include <symbolic/casadi.hpp>
#include <symbolic/stl_vector_tools.hpp>
#include <cstdlib>

namespace AD = CasADi;

#include <armadillo>
using namespace arma;

class CasadiPointExploreSystem {
public:
	CasadiPointExploreSystem();
	CasadiPointExploreSystem(const ObsType obs_type, const CostType cost_type, mat& R);

	double casadi_cost(const std::vector<mat>& X, const std::vector<mat>& U,
			const mat& P);
	mat casadi_cost_grad(const std::vector<mat>& X, const std::vector<mat>& U,
			const mat& P);

private:
	ObsType obs_type;
	CostType cost_type;
	mat R;

	AD::SXFunction cost_func, cost_grad_func;

	void init(ObsType obs_type, CostType cost_type, mat& R);

	AD::SXMatrix dist(AD::SXMatrix a, AD::SXMatrix b);
	AD::SXMatrix dynfunc(const AD::SXMatrix& x_t, const AD::SXMatrix& u_t);
	AD::SXMatrix obsfunc(const AD::SXMatrix& x, const AD::SXMatrix& t);
	AD::SXMatrix obsfunc_dist(const AD::SXMatrix& x, const AD::SXMatrix& t);
	AD::SXMatrix obsfunc_angle(const AD::SXMatrix& x, const AD::SXMatrix& t);
	AD::SXMatrix gauss_likelihood(const AD::SXMatrix& v);

	AD::SXMatrix cost_entropy(const std::vector<AD::SXMatrix>& X, const std::vector<AD::SXMatrix>& U,
										const std::vector<AD::SXMatrix>& P);
	AD::SXMatrix cost_entropy_wrapper(const AD::SXMatrix& XU_vec, const AD::SXMatrix& P_vec);

	AD::SXFunction casadi_cost_entropy_func();
	AD::SXFunction casadi_cost_entropy_grad_func();

	void setup_casadi_vars(const std::vector<mat>& X, const std::vector<mat>& U,
			const mat& P, double* XU_arr, double* P_arr);
};

#endif
