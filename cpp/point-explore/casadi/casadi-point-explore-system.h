#ifndef __CASADI_POINT_EXPLORE_SYSTEM_H__
#define __CASADI_POINT_EXPLORE_SYSTEM_H__

#include "casadi-system.h"

class CasadiPointExploreSystem;
#include "../point-explore-system.h"

#include <symbolic/casadi.hpp>
#include <symbolic/stl_vector_tools.hpp>
#include <cstdlib>

namespace AD = CasADi;

#include <armadillo>
using namespace arma;

class CasadiPointExploreSystem : CasadiSystem {
public:
	CasadiPointExploreSystem();
	CasadiPointExploreSystem(const ObsType obs_type, const CostType cost_type, mat& R);
	CasadiPointExploreSystem(const ObsType obs_type, const CostType cost_type, mat& R,
								int T, int M, int N, double DT, int X_DIM, int U_DIM, int Z_DIM, int Q_DIM, int R_DIM);
	~CasadiPointExploreSystem() { };

	double casadi_cost(const std::vector<mat>& X, const std::vector<mat>& U,
			const mat& P);
	mat casadi_cost_grad(const std::vector<mat>& X, const std::vector<mat>& U,
			const mat& P);

protected:
	int T, M, N, TOTAL_VARS, X_DIM, U_DIM, Z_DIM, Q_DIM, R_DIM;
	double DT;

	ObsType obs_type;
	CostType cost_type;
	mat R;

	AD::SXFunction cost_func, cost_grad_func;

	void init_dims(int T=10, int M=100, int N=1, double DT=1.0, int X_DIM=2, int U_DIM=2, int Z_DIM=1, int Q_DIM=2, int R_DIM=1);
	void init(ObsType obs_type, CostType cost_type, mat& R);

	AD::SXMatrix dist(AD::SXMatrix a, AD::SXMatrix b);
	AD::SXMatrix dynfunc(const AD::SXMatrix& x_t, const AD::SXMatrix& u_t);
	AD::SXMatrix obsfunc(const AD::SXMatrix& x, const AD::SXMatrix& t);
	AD::SXMatrix obsfunc_dist(const AD::SXMatrix& x, const AD::SXMatrix& t);
	AD::SXMatrix obsfunc_angle(const AD::SXMatrix& x, const AD::SXMatrix& t);
	AD::SXMatrix gauss_likelihood(const AD::SXMatrix& v);

	AD::SXMatrix cost_entropy(const std::vector<AD::SXMatrix>& X, const std::vector<AD::SXMatrix>& U,
										const std::vector<AD::SXMatrix>& P);
	AD::SXMatrix cost_platt(const std::vector<AD::SXMatrix>& X, const std::vector<AD::SXMatrix>& U,
									const std::vector<AD::SXMatrix>& P);

	AD::SXMatrix cost_wrapper(const AD::SXMatrix& XU_vec, const AD::SXMatrix& P_vec);
	AD::SXFunction casadi_cost_func();
	AD::SXFunction casadi_cost_grad_func();

	void setup_casadi_vars(const std::vector<mat>& X, const std::vector<mat>& U,
			const mat& P, double* XU_arr, double* P_arr);
};

#endif
