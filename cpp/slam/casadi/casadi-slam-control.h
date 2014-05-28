#ifndef _SLAM_CONTROL_CASADI_H__
#define _SLAM_CONTROL_CASADI_H__

#include <symbolic/casadi.hpp>
#include <symbolic/stl_vector_tools.hpp>

CasADi::SXFunction casadiCostFunc();

CasADi::SXFunction casadiCostGradFunc();

CasADi::SXFunction casadiCostFuncInfo();

CasADi::SXFunction casadiCostGradFuncInfo();

#endif
