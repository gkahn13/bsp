#include "fadbad-pr2-system.h"

/**
 * FadbadPR2System Constructors
 */

FadbadPR2System::FadbadPR2System(rave::RobotBasePtr r, Arm::ArmType arm_type, PR2System& sys) {
	object = sys.object.cast<bdouble>();

	j_min = sys.j_min.cast<bdouble>();
	j_max = sys.j_max.cast<bdouble>();

	Q = sys.Q.cast<bdouble>();
	R = sys.R.cast<bdouble>();

	arm = new FadbadArm(r, arm_type);
	cam = new FadbadCamera(arm, sys.cam->get_gripper_tool_to_sensor());

	this->sys = &sys;
}

/**
 * FadbadPR2System public methods
 */

VectorJb FadbadPR2System::dynfunc(const VectorJb& j, const VectorUb& u, const VectorQb& q, bool enforce_limits) {
	VectorUb u_noisy = u + q;
	for(int i=0; i < U_DIM; ++i) { u_noisy(i) *= DT; }
	VectorJb j_new = j + u_noisy;

	if (!enforce_limits) {
		return j_new;
	}

	for(int i=0; i < J_DIM; ++i) {
		j_new(i) = (j_new(i) > j_min(i)) ? j_new(i) : j_min(i);
		j_new(i) = (j_new(i) < j_max(i)) ? j_new(i) : j_max(i);
	}

	return j_new;
}

VectorZb FadbadPR2System::obsfunc(const VectorJb& j, const Vector3b& object, const VectorRb& r) {
	VectorZb z;
	z.segment<J_DIM>(0) = j;
	z.segment<3>(J_DIM) = object - cam->get_position(j);
	return z;
}

MatrixZb FadbadPR2System::delta_matrix(const VectorJb& j, const Vector3b& object, const bdouble alpha) {
	MatrixZb delta = MatrixZb::Identity();

	for(int i=0; i < J_DIM; ++i) {
		delta(i,i) = 1;
	}

	std::vector<std::vector<FadbadBeam3d> > beams = cam->get_beams(j, pcl);
	std::vector<FadbadTriangle3d> border = cam->get_border(beams);
	bdouble sd = cam->signed_distance(object, beams, border);

	bdouble sd_sigmoid = 1.0 - 1.0/(1.0 + exp(-alpha*sd));
	for(int i=J_DIM; i < Z_DIM; ++i) {
		delta(i,i) = sd_sigmoid;
	}

	return delta;
}

void FadbadPR2System::belief_dynamics(const VectorXb& x_t, const MatrixXb& sigma_t, const VectorUb& u_t, const bdouble alpha,
		VectorXb& x_tp1, MatrixXb& sigma_tp1) {
	// propagate dynamics
	x_tp1 = x_t;
	x_tp1.segment<J_DIM>(0) = dynfunc(x_t.segment<J_DIM>(0), u_t, VectorQb::Zero(), true);

	// propagate belief through dynamics
	Matrix<bdouble,X_DIM,X_DIM> A;
	Matrix<bdouble,X_DIM,Q_DIM> M;
	linearize_dynfunc(x_t, u_t, VectorQb::Zero(), A, M);

	MatrixXb sigma_tp1_bar = A*sigma_t*A.transpose() + M*Q*M.transpose();

	// propagate belief through observation
	Matrix<bdouble,Z_DIM,X_DIM> H;
	linearize_obsfunc(x_tp1, VectorRb::Zero(), H);

	MatrixZb delta = delta_matrix(x_tp1.segment<J_DIM>(0), x_tp1.segment<3>(J_DIM), alpha);
	Matrix<bdouble,X_DIM,Z_DIM> K = sigma_tp1_bar*H.transpose()*delta*(delta*H*sigma_tp1_bar*H.transpose()*delta + R).fullPivLu().solve(MatrixZb::Identity())*delta;
	sigma_tp1 = (MatrixXb::Identity() - K*H)*sigma_tp1_bar;
}


bdouble FadbadPR2System::cost(const StdVectorJb& J, const Vector3b& obj, const MatrixXb& sigma0, const StdVectorUb& U, const bdouble alpha) {
	bdouble cost = 0;

	VectorXb x_t, x_tp1 = VectorXb::Zero();
	MatrixXb sigma_t = sigma0, sigma_tp1 = MatrixXb::Zero();
	for(int t=0; t < TIMESTEPS-1; ++t) {
		x_t << J[t], obj;
		belief_dynamics(x_t, sigma_t, U[t], alpha, x_tp1, sigma_tp1);

		if (t < TIMESTEPS-2) {
			cost += alpha_belief*sigma_tp1.trace();
		} else {
			cost += alpha_final_belief*sigma_tp1.trace();
		}
		sigma_t = sigma_tp1;
	}

	Vector3b final_pos = cam->get_position(J.back());
	Vector3b e = obj - final_pos;
	cost += alpha_goal*e.squaredNorm();

	return cost;
}

bdouble FadbadPR2System::cost_gmm(const StdVectorJb& J, const MatrixJb& j_sigma0, const StdVectorUb& U,
			const std::vector<FadbadParticleGaussian>& particle_gmm, const bdouble alpha) {
	bdouble cost_gmm = 0;

	MatrixXb sigma0 = MatrixXb::Zero();
	sigma0.block<J_DIM,J_DIM>(0,0) = j_sigma0;
	for(int i=0; i < particle_gmm.size(); ++i) {
		sigma0.block<3,3>(J_DIM,J_DIM) = particle_gmm[i].cov;
//		cost_gmm += particle_gmm[i].pct*cost(J, particle_gmm[i].mean, sigma0, U, alpha);
		cost_gmm += cost(J, particle_gmm[i].mean, sigma0, U, alpha);
	}

	return cost_gmm;
}

void FadbadPR2System::cost_gmm_and_grad(StdVectorJ& J, const MatrixJ& j_sigma0, StdVectorU& U,
			const std::vector<ParticleGaussian>& particle_gmm, const double alpha, const StdVector3d& pcl_d,
			double& cost, VectorTOTAL& grad) {
	// convert to bdouble

	StdVectorJb J_b;
	for(int i=0; i < J.size(); ++i) { J_b.push_back(J[i].cast<bdouble>()); }

	MatrixJb j_sigma0_b = j_sigma0.cast<bdouble>();

	StdVectorUb U_b;
	for(int i=0; i < U.size(); ++i) { U_b.push_back(U[i].cast<bdouble>()); }

	std::vector<FadbadParticleGaussian> particle_gmm_b;
	for(int i=0; i < particle_gmm.size(); ++i) {
		Vector3b mean_b = particle_gmm[i].mean.cast<bdouble>();
		Matrix3b cov_b = particle_gmm[i].cov.cast<bdouble>();
		bdouble pct_b = static_cast<bdouble>(particle_gmm[i].pct);
		particle_gmm_b.push_back(FadbadParticleGaussian(mean_b, cov_b, pct_b));
	}

	bdouble alpha_b = static_cast<bdouble>(alpha);

	pcl.clear();
	for(int i=0; i < pcl_d.size(); ++i) { pcl.push_back(pcl_d[i].cast<bdouble>()); }

	bdouble cost_b = cost_gmm(J_b, j_sigma0_b, U_b, particle_gmm_b, alpha_b);

	// get cost
	cost = cost_b.x();

	// get grad
	cost_b.diff(0,1);
	int index = 0;
	for(int t=0; t < TIMESTEPS; ++t) {
		for(int i=0; i < J_DIM; ++i) {
			grad(index++) = J_b[t](i).d(0);
		}

		if (t < TIMESTEPS-1) {
			for(int i=0; i < U_DIM; ++i) {
				grad(index++) = U_b[t](i).d(0);
			}
		}
	}
}

/**
 * FadbadPR2System private methods
 */

void FadbadPR2System::linearize_dynfunc(const VectorXb& x, const VectorUb& u, const VectorQb& q,
		Matrix<bdouble,X_DIM,X_DIM>& A, Matrix<bdouble,X_DIM,Q_DIM>& M) {
	A.setIdentity();

	M.setZero();
	for(int j=0; j < Q_DIM; ++j) {
		for(int i=0; i < j+1; ++i) {
			M(i, j) = DT;
		}
	}
}

void FadbadPR2System::linearize_obsfunc(const VectorXb& x, const VectorRb& r,
		Matrix<bdouble,Z_DIM,X_DIM>& H) {
	H.setIdentity();
}
