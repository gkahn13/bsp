#include "../include/gmm.h"

namespace gmm {

void fit_gaussians_to_pf(const MatrixXd& P, std::vector<VectorXd>& obj_means, std::vector<MatrixXd>& obj_covs,
		std::vector<MatrixXd>& obj_particles) {
	// find modes and associated particles
	std::vector<VectorXd> modes;
	std::vector<std::vector<int>> mode_particle_indices;
	find_modes(P, modes, mode_particle_indices);

	std::cout << "number of modes found: " << obj_means.size() << "\n";
	std::cout << "mode_particle_indices.size(): " << mode_particle_indices.size() << "\n";

	// create matrices from associated mode_particle_indices
	// and then calculate covariance
	obj_means.clear();
	obj_covs.clear();
	obj_particles.clear();
	for(int i=0; i < mode_particle_indices.size(); ++i) {
		int num_mode_particles = mode_particle_indices[i].size();
		std::cout << "P_i size: " << num_mode_particles << "\n";
		MatrixXd P_i = MatrixXd::Zero(P.rows(), num_mode_particles);
		for(int m=0; m < num_mode_particles; ++m) {
			P_i.col(m) = P.col(mode_particle_indices[i][m]);
		}
		obj_means.push_back(P_i.rowwise().mean());
		obj_particles.push_back(P_i);
		MatrixXd P_i_centered = P_i.colwise() - P_i.rowwise().mean();
		MatrixXd obj_cov_i = (1/(double(num_mode_particles)-1))*(P_i_centered*P_i_centered.transpose());
		obj_covs.push_back(obj_cov_i);
	}
}

void find_modes(const MatrixXd& P,
		std::vector<VectorXd>& modes,
		std::vector<std::vector<int>>& mode_particle_indices) {
	for(int m=0; m < P.cols(); ++m) {
		VectorXd nearest_mode = find_nearest_mode(P.col(m), P);

		bool nearest_mode_exists = false;
		for(int i=0; i < modes.size(); ++i) {
			if ((modes[i] - nearest_mode).norm() < .05) {
				mode_particle_indices[i].push_back(m);
				nearest_mode_exists = true;
				break;
			}
		}

		if (!nearest_mode_exists) {
			modes.push_back(nearest_mode);
			mode_particle_indices.push_back(std::vector<int>());
		}
	}
}

VectorXd find_nearest_mode(const VectorXd& p, const MatrixXd& P) {
	VectorXd new_mean = p, mean = INFINITY*VectorXd::Ones(p.rows(), p.cols());

	while((mean - new_mean).norm() > 1e-3) {
		mean = new_mean;
		VectorXd kernel_weight = (-(1/1)*(P.colwise() - mean).colwise().norm()).array().exp(); // TODO: window size
		new_mean = kernel_weight.transpose().replicate(2,1).cwiseProduct(P).rowwise().sum() / kernel_weight.sum();
	}

	return mean;
}


}
