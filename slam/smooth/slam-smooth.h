#ifndef _SLAM_SMOOTH_H__
#define _SLAM_SMOOTH_H__

#include "../slam.h"

#define TIMESTEPS_SMOOTH 60
const int T_SMOOTH = TIMESTEPS_SMOOTH;

std::vector<Matrix<U_DIM> > smoothTraj(const std::vector<Matrix<C_DIM> >& X_unsmooth);

#endif
