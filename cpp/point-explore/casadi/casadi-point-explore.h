#ifndef __CASADI_POINT_EXPLORE_H__
#define __CASADI_POINT_EXPLORE_H__

#include "../point-explore-system.h"

#include <symbolic/casadi.hpp>
#include <symbolic/stl_vector_tools.hpp>
#include <cstdlib>

namespace casadi_point_explore {

namespace AD = CasADi;

AD::SXMatrix dist(AD::SXMatrix a, AD::SXMatrix b);

AD::SXMatrix dynfunc(const AD::SXMatrix& x_t, const AD::SXMatrix& u_t);

AD::SXMatrix obsfunc(const AD::SXMatrix& x, const AD::SXMatrix& t);

AD::SXMatrix gaussLikelihood(const AD::SXMatrix& v);

AD::SXMatrix differential_entropy(const std::vector<AD::SXMatrix>& X, const std::vector<AD::SXMatrix>& U,
									const std::vector<AD::SXMatrix>& P);

AD::SXMatrix differential_entropy_wrapper(const AD::SXMatrix& XU_vec, const AD::SXMatrix& P_vec);

//void setup_casadi_vars(const std::vector<Matrix<N*X_DIM> >& X, const std::vector<Matrix<N*U_DIM> >& U,
//					 const std::vector<Matrix<X_DIM> >& P, double* XU_arr, double* P_arr);

AD::SXFunction casadi_differential_entropy_func();

AD::SXFunction casadi_differential_entropy_gradfunc();

AD::SXFunction casadi_differential_entropy_diaghessfunc();

}

#endif
