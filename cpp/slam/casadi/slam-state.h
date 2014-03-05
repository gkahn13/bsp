#ifndef _SLAM_STATE_CASADI_H__
#define _SLAM_STATE_CASADI_H__

#include <iostream>
#include <iomanip>
#include <fstream>
#include <ctime>
#include <symbolic/casadi.hpp>
#include <symbolic/stl_vector_tools.hpp>
#include <cstdlib>

CasADi::SXFunction casadiCostFunc(int timesteps);

CasADi::SXFunction casadiCostGradFunc(int timesteps);

#endif

