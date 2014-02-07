#ifndef _SLAM_TRAJ_H__
#define _SLAM_TRAJ_H__

#include "../slam.h"

bool initTraj(const Matrix<C_DIM>& cStart, const Matrix<C_DIM>& cEnd, std::vector<Matrix<U_DIM> >& U);

#endif
