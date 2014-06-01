#ifndef _SLAM_STATE_CASADI_H__
#define _SLAM_STATE_CASADI_H__

#include <symbolic/casadi.hpp>
#include <symbolic/stl_vector_tools.hpp>

CasADi::SXFunction casadiCostFunc();

CasADi::SXFunction casadiCostGradFunc();

#endif

