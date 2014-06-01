#ifndef __CASADI_EXPLORE_SYSTEM_H__
#define __CASADI_EXPLORE_SYSTEM_H__

#include "../casadi-system.h"

class CasadiExploreSystem;
#include "explore-system.h"

#include <symbolic/casadi.hpp>
#include <symbolic/stl_vector_tools.hpp>
#include <cstdlib>

namespace AD = CasADi;

#include <armadillo>
using namespace arma;

class CasadiExploreSystem : public virtual CasadiSystem {
public:
	CasadiExploreSystem();
	CasadiExploreSystem(const ObsType obs_type, const CostType cost_type, mat& R);
	CasadiExploreSystem(const ObsType obs_type, const CostType cost_type, mat& R,
								int T, int M, int N, double DT, int X_DIM, int U_DIM, int Z_DIM, int Q_DIM, int R_DIM);
	~CasadiExploreSystem() { };

protected:
	int N;

	ObsType obs_type;
	CostType cost_type;

	void init_dims(int T=10, int M=100, int N=1, double DT=1.0, int X_DIM=2, int U_DIM=2, int Z_DIM=1, int Q_DIM=2, int R_DIM=1);
	void init(ObsType obs_type, CostType cost_type, mat& R);

	AD::SXMatrix dynfunc(const AD::SXMatrix& x_t, const AD::SXMatrix& u_t);
	AD::SXMatrix obsfunc(const AD::SXMatrix& x, const AD::SXMatrix& t);
	AD::SXMatrix obsfunc_dist(const AD::SXMatrix& x, const AD::SXMatrix& t);
	AD::SXMatrix obsfunc_angle(const AD::SXMatrix& x, const AD::SXMatrix& t);

	AD::SXMatrix cost_wrapper(const AD::SXMatrix& XU_vec, const AD::SXMatrix& P_vec);
};

#endif
