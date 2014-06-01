#ifndef __CASADI_BOXES_SYSTEM_H__
#define __CASADI_BOXES_SYSTEM_H__

#include "../casadi-system.h"

class CasadiBoxesSystem;
#include "boxes-system.h"

#include <symbolic/casadi.hpp>
#include <symbolic/stl_vector_tools.hpp>
#include <cstdlib>

namespace AD = CasADi;

#include <armadillo>
using namespace arma;

class CasadiBoxesSystem : public virtual CasadiSystem {
public:
	CasadiBoxesSystem();
	CasadiBoxesSystem(mat& box_dims, const ObsType obs_type, const CostType cost_type, mat& R);
	CasadiBoxesSystem(mat& box_dims, const ObsType obs_type, const CostType cost_type, mat& R,
								int T, int M, int N, double DT, int X_DIM, int U_DIM, int Z_DIM, int Q_DIM, int R_DIM);
	~CasadiBoxesSystem() { };

protected:
	int N;

	mat box_dims;
	ObsType obs_type;
	CostType cost_type;

	void init_dims(int T=10, int M=100, int N=1, double DT=1.0, int X_DIM=2, int U_DIM=2, int Z_DIM=1, int Q_DIM=2, int R_DIM=1);
	void init(mat& box_dims, ObsType obs_type, CostType cost_type, mat& R);

	AD::SXMatrix dynfunc(const AD::SXMatrix& x_t, const AD::SXMatrix& u_t);
	AD::SXMatrix obsfunc(const AD::SXMatrix& x, const AD::SXMatrix& b);

	AD::SXMatrix cost_wrapper(const AD::SXMatrix& XU_vec, const AD::SXMatrix& P_vec);
};

#endif
