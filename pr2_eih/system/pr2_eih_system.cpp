#include "pr2_eih_system.h"

/**
 *
 * PR2EihSystem constructors
 *
 */

PR2EihSystem::PR2EihSystem(pr2_sim::Simulator *s, pr2_sim::Arm *a, pr2_sim::Camera *c) : sim(s), arm(a), cam(c) {
	if (sim == NULL) {
		sim = new pr2_sim::Simulator(true, false);
	}

	if (arm == NULL) {
		arm = new pr2_sim::Arm(pr2_sim::Arm::ArmType::right, sim);
	}

	if (cam == NULL) {
		cam = new pr2_sim::Camera(arm, sim);
	}

	arm->get_joint_limits(j_min, j_max);
	u_min = j_min;
	u_max = j_max;

	Q = (M_PI/4)*MatrixQ::Identity();
//	Q = 1e-2*MatrixQ::Identity();
	VectorR R_diag;
	R_diag << (M_PI/4)*VectorJ::Ones(), 5*Vector3d::Ones();
//	R_diag << 1e-2*VectorJ::Ones(), 5*Vector3d::Ones();
	R = R_diag.asDiagonal();

	arm->set_posture(pr2_sim::Arm::Posture::mantis);
}

/**
 *
 * PR2EihSystem public methods
 *
 */

VectorJ PR2EihSystem::dynfunc(const VectorJ& j, const VectorU& u, const VectorQ& q, bool enforce_limits) {
	if (!enforce_limits) {
		return (j + DT*(u + q));
	} else {
		VectorJ j_new = j + DT*(u+q);

		for(int i=0; i < J_DIM; ++i) {
			j_new(i) = std::max(j_new(i), j_min(i));
			j_new(i) = std::min(j_new(i), j_max(i));
		}

		return j_new;
	}
}

void PR2EihSystem::get_limits(VectorJ& j_min, VectorJ& j_max, VectorU& u_min, VectorU& u_max) {
	j_min = this->j_min;
	j_max = this->j_max;
	u_min = this->u_min;
	u_max = this->u_max;
}

/**
 * \brief Observation function is the joint positions and (object - camera)
 */
VectorZ PR2EihSystem::obsfunc(const VectorJ& j, const Vector3d& object, const VectorR& r) {
	VectorJ j_orig = arm->get_joints();
	arm->set_joints(j);

	VectorZ z;
	z.segment<J_DIM>(0) = j;
	z.segment<3>(J_DIM) = object - cam->get_pose().block<3,1>(0,3);

	arm->set_joints(j_orig);

	return z;
}

/**
 * \brief Delta for (object - camera) determined by signed distance of object from
 *        truncated view frustum of camera for arm joint values j
 */
MatrixZ PR2EihSystem::delta_matrix(const VectorJ& j, const Vector3d& object, const double alpha,
		const std::vector<geometry3d::Triangle>& obstacles) {
	VectorJ j_orig = arm->get_joints();
	arm->set_joints(j);

	MatrixZ delta = MatrixZ::Identity();

	for(int i=0; i < J_DIM; ++i) {
		delta(i,i) = 1; // always observe joints
	}

	std::vector<geometry3d::Pyramid> truncated_frustum = cam->truncated_view_frustum(obstacles, false);
	double sd = cam->signed_distance(object, truncated_frustum);

	double sd_sigmoid = 1.0 - 1.0/(1.0 + exp(-alpha*sd));
	for(int i=J_DIM; i < Z_DIM; ++i) {
		delta(i,i) = sd_sigmoid;
	}

	arm->set_joints(j_orig);

	return delta;
}

void PR2EihSystem::belief_dynamics(const VectorX& x_t, const MatrixX& sigma_t, const VectorU& u_t, const double alpha,
		const std::vector<geometry3d::Triangle>& obstacles, VectorX& x_tp1, MatrixX& sigma_tp1) {
	// propagate dynamics
	x_tp1 = x_t;
	x_tp1.segment<J_DIM>(0) = dynfunc(x_t.segment<J_DIM>(0), u_t, VectorQ::Zero(), true);

	// propagate belief through dynamics
	Matrix<double,X_DIM,X_DIM> A;
	Matrix<double,X_DIM,Q_DIM> M;
	linearize_dynfunc(x_t, u_t, VectorQ::Zero(), A, M);

	MatrixX sigma_tp1_bar = A*sigma_t*A.transpose() + M*Q*M.transpose();

	// propagate belief through observation
	Matrix<double,Z_DIM,X_DIM> H;
	linearize_obsfunc(x_tp1, VectorR::Zero(), H);

	MatrixZ delta = delta_matrix(x_tp1.segment<J_DIM>(0), x_tp1.segment<3>(J_DIM), alpha, obstacles);
	Matrix<double,X_DIM,Z_DIM> K = sigma_tp1_bar*H.transpose()*delta*(delta*H*sigma_tp1_bar*H.transpose()*delta + R).inverse()*delta;
	sigma_tp1 = (MatrixX::Identity() - K*H)*sigma_tp1_bar;
}

void PR2EihSystem::execute_control_step(const VectorJ& j_t_real, const VectorJ& j_t, const VectorU& u_t,
		const std::vector<Gaussian3d>& obj_gaussians_t,
		const std::vector<geometry3d::Triangle>& obstacles,
		VectorJ& j_tp1_real, VectorJ& j_tp1, std::vector<Gaussian3d>& obj_gaussians_tp1) {
	// find next real state from input + noise
	VectorQ control_noise = VectorQ::Zero();// + chol(.2*Q)*randn<vec>(Q_DIM);
	VectorR obs_noise = VectorR::Zero();// + chol(.2*R)*randn<vec>(R_DIM);
	j_tp1_real = dynfunc(j_t_real, u_t, control_noise, true);

	// max-likelihood dynfunc estimate
	j_tp1 = dynfunc(j_t, u_t, VectorQ::Zero(), true);

	arm->set_joints(j_tp1_real);
	std::vector<geometry3d::Pyramid> truncated_frustum = cam->truncated_view_frustum(obstacles, true);

	obj_gaussians_tp1.clear();
	for(const Gaussian3d& obj_gaussian_t : obj_gaussians_t) {
		const MatrixP& P_t = obj_gaussian_t.particles;
		VectorP W_tp1 = update_particle_weights(j_tp1, P_t, (1/double(M_DIM))*VectorP::Ones(), obstacles);
		MatrixP P_tp1 = low_variance_sampler(P_t, W_tp1);

		obj_gaussians_tp1.push_back(Gaussian3d(P_tp1));
	}
}

double PR2EihSystem::cost(const StdVectorJ& J, const MatrixJ& j_sigma0, const StdVectorU& U, const std::vector<Gaussian3d>& obj_gaussians,
		const double alpha, const std::vector<geometry3d::Triangle>& obstacles) {
	double cost = 0;

	for(const Gaussian3d& obj_gaussian : obj_gaussians) {
		VectorX x_t, x_tp1 = VectorX::Zero();
		MatrixX sigma_t = MatrixX::Zero(), sigma_tp1 = MatrixX::Zero();
		sigma_t.block<J_DIM,J_DIM>(0,0) = j_sigma0;
		sigma_t.block<3,3>(J_DIM,J_DIM) = obj_gaussian.cov;
		for(int t=0; t < TIMESTEPS-1; ++t) {
			x_t << J[t], obj_gaussian.mean;
			belief_dynamics(x_t, sigma_t, U[t], alpha, obstacles, x_tp1, sigma_tp1);

			cost += alpha_control*U[t].squaredNorm();
			// TODO: only penalize object?
			if (t < TIMESTEPS-2) {
				cost += alpha_belief*sigma_tp1.trace();
			} else {
				cost += alpha_final_belief*sigma_tp1.trace();
			}
			sigma_t = sigma_tp1;
		}
	}

	return cost;
}

VectorTOTAL PR2EihSystem::cost_grad(StdVectorJ& J, const MatrixJ& j_sigma0, StdVectorU& U, const std::vector<Gaussian3d>& obj_gaussians,
		const double alpha, const std::vector<geometry3d::Triangle>& obstacles) {
	VectorTOTAL grad;

	double orig, cost_p, cost_m;
	int index = 0;
	for(int t=0; t < TIMESTEPS; ++t) {
		for(int i=0; i < J_DIM; ++i) {
			orig = J[t][i];

			J[t](i) = orig + step;
			cost_p = cost(J, j_sigma0, U, obj_gaussians, alpha, obstacles);

			J[t](i) = orig - step;
			cost_m = cost(J, j_sigma0, U, obj_gaussians, alpha, obstacles);

			grad(index++) = (cost_p - cost_m) / (2*step);
		}

		if (t < TIMESTEPS-1) {
			for(int i=0; i < U_DIM; ++i) {
				orig = U[t][i];

				U[t](i) = orig + step;
				cost_p = cost(J, j_sigma0, U, obj_gaussians, alpha, obstacles);

				U[t](i) = orig - step;
				cost_m = cost(J, j_sigma0, U, obj_gaussians, alpha, obstacles);

				grad(index++) = (cost_p - cost_m) / (2*step);
			}
		}
	}
	return grad;
}

VectorP PR2EihSystem::update_particle_weights(const VectorJ& j_tp1, const MatrixP& P_t, const VectorP& W_t, const std::vector<geometry3d::Triangle>& obstacles) {
	arm->set_joints(j_tp1);
	std::vector<geometry3d::Pyramid> truncated_frustum = cam->truncated_view_frustum(obstacles, false);

	VectorP W_tp1;
	// for each particle, weight by sigmoid of signed distance
	for(int m=0; m < M_DIM; ++m) {
		double sd = cam->signed_distance(P_t.col(m), truncated_frustum);
		double sigmoid_sd = 1.0/(1.0 + exp(-alpha_particle_sd*sd));
		W_tp1(m) = W_t(m)*sigmoid_sd;
	}
	W_tp1 = W_tp1 / W_tp1.sum();

	return W_tp1;
}

MatrixP PR2EihSystem::low_variance_sampler(const MatrixP& P, const VectorP& W) {
	MatrixP P_sampled(3,M_DIM);

	double r = (1/double(M_DIM))*(rand() / double(RAND_MAX));
	double c = W(0);
	int i = 0;
	for(int m=0; m < M_DIM; ++m) {
		double u = r + m * (1/double(M_DIM));
		while (u > c) {
			c += W(++i);
		}
		P_sampled.col(m) = P.col(i) + .002*Vector3d::Random();
	}

	return P_sampled;
}

double PR2EihSystem::entropy(const StdVectorJ& J, const StdVectorU& U, const MatrixP& P, const std::vector<geometry3d::Triangle>& obstacles) {
	double entropy = 0;

	MatrixP P_t = P, P_tp1 = P;
	VectorP W_t = (1/double(M_DIM))*VectorP::Ones(), W_tp1;
	for(int t=0; t < TIMESTEPS-1; ++t) {
		VectorJ j_tp1 = dynfunc(J[t], U[t], VectorQ::Zero());

		W_tp1 = update_particle_weights(j_tp1, P_t, W_t, obstacles);

//		entropy += (-W_tp1.array() * W_tp1.array().log()).sum();
		for(int m=0; m < M_DIM; ++m) { // safe log
			if (W_tp1(m) > 1e-8) {
				entropy += -W_tp1(m)*log(W_tp1(m));
			}
		}

		W_t = W_tp1;

//		// TEMP: try resampling
//		P_tp1 = low_variance_sampler(P_t, W_tp1);
//
//		W_t = (1/double(M_DIM))*VectorP::Ones();
//		P_t = P_tp1;
	}
	arm->set_joints(J[0]);

	return entropy;
}

VectorTOTAL PR2EihSystem::entropy_grad(StdVectorJ& J, StdVectorU& U, const MatrixP& P, const std::vector<geometry3d::Triangle>& obstacles) {
	VectorTOTAL grad;

	double orig, cost_p, cost_m;
	int index = 0;
	for(int t=0; t < TIMESTEPS; ++t) {
		for(int i=0; i < J_DIM; ++i) {
			orig = J[t][i];

			J[t](i) = orig + step;
			cost_p = entropy(J, U, P, obstacles);

			J[t](i) = orig - step;
			cost_m = entropy(J, U, P, obstacles);

			grad(index++) = (cost_p - cost_m) / (2*step);
		}

		if (t < TIMESTEPS-1) {
			for(int i=0; i < U_DIM; ++i) {
				orig = U[t][i];

				U[t](i) = orig + step;
				cost_p = entropy(J, U, P, obstacles);

				U[t](i) = orig - step;
				cost_m = entropy(J, U, P, obstacles);

				grad(index++) = (cost_p - cost_m) / (2*step);
			}
		}
	}
	return grad;
}

void PR2EihSystem::plot(const StdVectorJ& J, const std::vector<Gaussian3d>& obj_gaussians,
		const std::vector<geometry3d::Triangle>& obstacles, bool pause) {
	sim->clear_plots();
	for(int t=0; t < J.size(); ++t) {
		arm->set_joints(J[t]);
		Matrix4d pose_world = sim->transform_from_to(cam->get_pose(), "base_link", "world");
		sim->plot_transform(pose_world);
	}

	for(const geometry3d::Triangle& obstacle : obstacles) {
		obstacle.plot(*sim, "base_link", {0,0,1}, true, 0.25);
	}

	for(const Gaussian3d& obj_gaussian : obj_gaussians) {
		Vector3d obj_world = sim->transform_from_to(obj_gaussian.mean, "base_link", "world");
		Matrix4d obj_sigma_T = Matrix4d::Zero();
		obj_sigma_T.block<3,3>(0,0) = obj_gaussian.cov;
		Matrix3d obj_sigma_world = sim->transform_from_to(obj_sigma_T, "base_link", "world").block<3,3>(0,0);
		sim->plot_gaussian(obj_world, obj_sigma_world, {0,1,0});

		for(int m=0; m < M_DIM; ++m) {
			Vector3d particle = obj_gaussian.particles.col(m);
			sim->plot_point(sim->transform_from_to(particle, "base_link", "world"), {0,1,0}, .001);
		}
	}

	arm->set_joints(J.back());
//	cam->plot({1,0,0});
	std::vector<geometry3d::Pyramid> truncated_frustum = cam->truncated_view_frustum(obstacles, false);
	for(const geometry3d::Pyramid& pyramid : truncated_frustum) {
		pyramid.plot(*sim, "base_link", {1,0,0}, true, true, 0.15);
	}

	arm->set_joints(J[0]);

	if (pause) {
		LOG_INFO("Plotted, press enter");
		std::cin.ignore();
	}
}

void PR2EihSystem::plot(const StdVectorJ& J, const MatrixP& P,
		const std::vector<geometry3d::Triangle>& obstacles, bool pause) {
	sim->clear_plots();
	for(int t=0; t < J.size(); ++t) {
		arm->set_joints(J[t]);
		Matrix4d pose_world = sim->transform_from_to(cam->get_pose(), "base_link", "world");
		sim->plot_transform(pose_world);
	}

	for(const geometry3d::Triangle& obstacle : obstacles) {
		obstacle.plot(*sim, "base_link", {0,0,1}, true, 0.25);
	}

	for(int m=0; m < M_DIM; ++m) {
		Vector3d particle = P.col(m);
		sim->plot_point(sim->transform_from_to(particle, "base_link", "world"), {0,1,0}, .005);
	}

	arm->set_joints(J.back());
//	cam->plot({1,0,0});
	std::vector<geometry3d::Pyramid> truncated_frustum = cam->truncated_view_frustum(obstacles, false);
	for(const geometry3d::Pyramid& pyramid : truncated_frustum) {
		pyramid.plot(*sim, "base_link", {1,0,0}, true, true, 0.15);
	}

	arm->set_joints(J[0]);

	if (pause) {
		LOG_INFO("Plotted, press enter");
		std::cin.ignore();
	}
}

/**
 *
 * PR2EihSystem private methods
 *
 */

void PR2EihSystem::linearize_dynfunc(const VectorX& x, const VectorU& u, const VectorQ& q,
		Matrix<double,X_DIM,X_DIM>& A, Matrix<double,X_DIM,Q_DIM>& M) {
	A = Matrix<double,X_DIM,X_DIM>::Identity();

	M.setZero();
	for(int j=0; j < Q_DIM; ++j) {
		for(int i=0; i < j+1; ++i) {
			M(i, j) = DT;
		}
	}

//	A.setZero();
//	VectorQ q_zero = VectorQ::Zero();
//	VectorX x_p = x, x_m = x;
//	for(int i=0; i < X_DIM; ++i) {
//		x_p(i) += step; x_m(i) -= step;
//		A.block<J_DIM,1>(0, i) = (dynfunc(x_p.segment<J_DIM>(0), u, q_zero) - dynfunc(x_m.segment<J_DIM>(0), u, q_zero)) / (2*step);
//		x_p(i) = x(i); x_m(i) = x(i);
//	}
//
//	M.setZero();
//	VectorQ q_p = q, q_m = q;
//	for(int i=0; i < Q_DIM; ++i) {
//		q_p(i) += step; q_m(i) -= step;
//		M.block<J_DIM,1>(0, i) = (dynfunc(x.segment<J_DIM>(0), u, q_p) - dynfunc(x.segment<J_DIM>(0), u, q_m)) / (2*step);
//	}
}

void PR2EihSystem::linearize_obsfunc(const VectorX& x, const VectorR& r,
		Matrix<double,Z_DIM,X_DIM>& H) {
	H = Matrix<double,Z_DIM,X_DIM>::Identity();

//	H.setZero();
//	VectorX x_p = x, x_m = x;
//	for(int i=0; i < X_DIM; ++i) {
//		x_p(i) += step; x_m(i) -= step;
//		H.block<Z_DIM,1>(0, i) = (obsfunc(x_p.segment<J_DIM>(0), x_p.segment<3>(J_DIM), r) -
//				obsfunc(x_m.segment<J_DIM>(0), x_m.segment<3>(J_DIM), r)) / (2*step);
//		x_p(i) = x(i); x_m(i) = x(i);
//	}
}

double PR2EihSystem::gauss_likelihood(const Vector3d& v, const Matrix3d& S) {
	Matrix3d Sf = S.llt().matrixL();
	Vector3d M = Sf.lu().solve(v);

	double E = -0.5*M.dot(M);
	double C = pow(2*M_PI, S.cols()/2) * Sf.diagonal().prod();
	double w = exp(E) / C;

	return w;
}

