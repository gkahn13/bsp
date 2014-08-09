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
	VectorR R_diag;
	R_diag << (M_PI/4)*VectorJ::Ones(), 5*Vector3d::Ones();
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

void PR2EihSystem::get_limits(VectorJ& j_min, VectorJ& j_max, VectorU& u_min, VectorU& u_max) {
	j_min = this->j_min;
	j_max = this->j_max;
	u_min = this->u_min;
	u_max = this->u_max;
}

double PR2EihSystem::cost(const StdVectorJ& J, const MatrixJ& j_sigma0, const StdVectorU& U, const Vector3d& obj, const Matrix3d& obj_sigma0,
		const double alpha, const std::vector<geometry3d::Triangle>& obstacles) {
	double cost = 0;

	VectorX x_t, x_tp1 = VectorX::Zero();
	MatrixX sigma_t = MatrixX::Zero(), sigma_tp1 = MatrixX::Zero();
	sigma_t.block<J_DIM,J_DIM>(0,0) = j_sigma0;
	sigma_t.block<3,3>(J_DIM,J_DIM) = obj_sigma0;
	for(int t=0; t < TIMESTEPS-1; ++t) {
		x_t << J[t], obj;
		belief_dynamics(x_t, sigma_t, U[t], alpha, obstacles, x_tp1, sigma_tp1);

		if (t < TIMESTEPS-2) {
			cost += alpha_belief*sigma_tp1.trace();
		} else {
			cost += alpha_final_belief*sigma_tp1.trace();
		}
		sigma_t = sigma_tp1;
	}

	return cost;
}

VectorTOTAL PR2EihSystem::cost_grad(StdVectorJ& J, const MatrixJ& j_sigma0, StdVectorU& U, const Vector3d& obj, const Matrix3d& obj_sigma0,
		const double alpha, const std::vector<geometry3d::Triangle>& obstacles) {
	VectorTOTAL grad;

	double orig, cost_p, cost_m;
	int index = 0;
	for(int t=0; t < TIMESTEPS; ++t) {
		for(int i=0; i < J_DIM; ++i) {
			orig = J[t][i];

			J[t](i) = orig + step;
			cost_p = cost(J, j_sigma0, U, obj, obj_sigma0, alpha, obstacles);

			J[t](i) = orig - step;
			cost_m = cost(J, j_sigma0, U, obj, obj_sigma0, alpha, obstacles);

			grad(index++) = (cost_p - cost_m) / (2*step);
		}

		if (t < TIMESTEPS-1) {
			for(int i=0; i < U_DIM; ++i) {
				orig = U[t][i];

				U[t](i) = orig + step;
				cost_p = cost(J, j_sigma0, U, obj, obj_sigma0, alpha, obstacles);

				U[t](i) = orig - step;
				cost_m = cost(J, j_sigma0, U, obj, obj_sigma0, alpha, obstacles);

				grad(index++) = (cost_p - cost_m) / (2*step);
			}
		}
	}
	return grad;
}

void PR2EihSystem::plot(const StdVectorJ& J, const Vector3d& obj, const Matrix3d& obj_sigma,
		const std::vector<geometry3d::Triangle>& obstacles, bool pause) {
	VectorJ curr_joints = arm->get_joints();

	sim->clear_plots();
	for(int t=0; t < J.size(); ++t) {
		arm->set_joints(J[t]);
		Matrix4d pose_world = sim->transform_from_to(cam->get_pose(), "base_link", "world");
		sim->plot_transform(pose_world);
	}

	Vector3d obj_world = sim->transform_from_to(obj, "base_link", "world");
	Matrix4d obj_sigma_T = Matrix4d::Zero();
	obj_sigma_T.block<3,3>(0,0) = obj_sigma;
	Matrix3d obj_sigma_world = sim->transform_from_to(obj_sigma_T, "base_link", "world").block<3,3>(0,0);
	sim->plot_gaussian(obj_world, obj_sigma_world, {0,1,0});

	for(const geometry3d::Triangle& obstacle : obstacles) {
		obstacle.plot(*sim, "base_link", {0,0,1}, true, 0.25);
	}

	cam->plot({1,0,0});

	arm->set_joints(curr_joints);

	if (pause) {
		std::cout << "Plotted, press enter\n";
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
