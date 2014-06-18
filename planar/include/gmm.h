#ifndef __GMM_H__
#define __GMM_H__

#include "planar-utils.h"

#include <iostream>

#include <Eigen/Eigen>
#include <Eigen/StdVector>
using namespace Eigen;

template <size_t _dim0, size_t _dim1>
using mat = Matrix<double, _dim0, _dim1>;

template <size_t _dim>
using vec = Matrix<double, _dim, 1>;


namespace gmm {

void fit_gaussians_to_pf(const MatrixXd& P, std::vector<VectorXd>& obj_means, std::vector<MatrixXd>& obj_covs,
		std::vector<MatrixXd>& obj_particles);

void find_modes(const MatrixXd& P, std::vector<VectorXd>& modes, std::vector<std::vector<int>>& mode_particle_indices);

VectorXd find_nearest_mode(const VectorXd& p, const MatrixXd& P);

}

#endif
