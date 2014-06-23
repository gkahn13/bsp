#include "pr2-system.h"

/**
 * PR2System Constructors and Initializers
 */

PR2System::PR2System(Vector3d& o) : object(o) {
	brett = new PR2();
	arm = brett->rarm;
	cam = brett->rcam;
	init();
}

PR2System::PR2System(Vector3d& o, Arm::ArmType arm_type, bool view) : object(o) {
	brett = new PR2(view);
	if (arm_type == Arm::ArmType::right) {
		arm = brett->rarm;
		cam = brett->rcam;
	} else {
		arm = brett->larm;
		cam = brett->lcam;
	}
	init();
}

PR2System::PR2System(Vector3d& o, Arm::ArmType arm_type, std::string env_file, std::string robot_name, bool view) : object(o) {
	brett = new PR2(env_file, robot_name, view);
	if (arm_type == Arm::ArmType::right) {
		arm = brett->rarm;
		cam = brett->rcam;
	} else {
		arm = brett->larm;
		cam = brett->lcam;
	}
	init();
}

void PR2System::init() {
	arm->get_limits(j_min, j_max);
	u_min = j_min;
	u_max = j_max;

	Q = (M_PI/4)*MatrixQ::Identity();
	VectorR R_diag;
	R_diag << (M_PI/4)*VectorJ::Ones(), 5*Vector3d::Ones();
	R = R_diag.asDiagonal();
}

/**
 * PR2System Public methods
 */

VectorJ PR2System::dynfunc(const VectorJ& j, const VectorU& u, const VectorQ& q, bool enforce_limits) {
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
VectorZ PR2System::obsfunc(const VectorJ& j, const Vector3d& object, const VectorR& r) {
	VectorZ z;
	z.segment<J_DIM>(0) = j;
	z.segment<3>(J_DIM) = object - cam->get_position();
	return z;
}

/**
 * \brief Delta for (object - camera) determined by FOV of camera for arm joint values j
 */
MatrixZ PR2System::delta_matrix(const VectorJ& j, const Vector3d& object, const double alpha) {
	MatrixZ delta = MatrixZ::Identity();

	for(int i=0; i < J_DIM; ++i) {
		delta(i,i) = 1; // TODO: should this depend on SD of joints?
	}

	VectorJ j_orig = arm->get_joint_values();
	arm->set_joint_values(j);

	std::vector<std::vector<Beam3d> > beams = cam->get_beams();
	std::vector<Triangle3d> border = cam->get_border(beams);
	double sd = cam->signed_distance(object, beams, border);

	arm->set_joint_values(j_orig);

	double sd_sigmoid = 1.0 - 1.0/(1.0 + exp(-alpha*sd));
	for(int i=J_DIM; i < Z_DIM; ++i) {
		delta(i,i) = sd_sigmoid;
	}

	return delta;
}

void PR2System::belief_dynamics(const VectorX& x_t, const MatrixX& sigma_t, const VectorU& u_t, const double alpha,
			VectorX& x_tp1, MatrixX& sigma_tp1) {
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
	Matrix<double,Z_DIM,R_DIM> N;
	linearize_obsfunc(x_tp1, VectorR::Zero(), H, N);

	MatrixZ delta = delta_matrix(x_tp1.segment<J_DIM>(0), x_tp1.segment<3>(J_DIM), alpha);
	Matrix<double,X_DIM,Z_DIM> K = sigma_tp1_bar*H.transpose()*delta*(delta*H*sigma_tp1_bar*H.transpose()*delta + R).inverse()*delta;
	sigma_tp1 = (MatrixX::Identity() - K*H)*sigma_tp1_bar;
}

void PR2System::execute_control_step(const VectorJ& j_t_real, const VectorJ& j_t, const VectorU& u_t, const MatrixP& P_t,
				VectorJ& j_tp1_real, VectorJ& j_tp1, MatrixP& P_tp1) {
	// find next real state from input + noise
	VectorQ control_noise = VectorQ::Zero();// + chol(.2*Q)*randn<vec>(Q_DIM);
	VectorR obs_noise = VectorR::Zero();// + chol(.2*R)*randn<vec>(R_DIM);
	j_tp1_real = dynfunc(j_t_real, u_t, control_noise, true);
	VectorZ z_tp1_real = obsfunc(j_tp1_real, this->object, obs_noise);

	// max-likelihood dynfunc estimate
	j_tp1 = dynfunc(j_t, u_t, VectorQ::Zero(), true);

	MatrixZ delta = delta_matrix(j_tp1_real, object, INFINITY);
	Vector3d delta_real = delta.diagonal().segment<3>(J_DIM);
	update_particles(j_tp1, delta_real(0), z_tp1_real, P_t, P_tp1);
}

void PR2System::get_limits(VectorJ& j_min, VectorJ& j_max, VectorU& u_min, VectorU& u_max) {
	j_min = this->j_min;
	j_max = this->j_max;
	u_min = this->u_min;
	u_max = this->u_max;
}

double PR2System::cost(const StdVectorJ& J, const Vector3d& obj, const MatrixX& sigma0, const StdVectorU& U, const double alpha) {
	double cost = 0;

	VectorX x_t, x_tp1 = VectorX::Zero();
	MatrixX sigma_t = sigma0, sigma_tp1 = MatrixX::Zero();
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

	VectorJ j_orig = arm->get_joint_values();
	arm->set_joint_values(J.back());

	Vector3d final_pos = cam->get_position();
	Vector3d e = obj - final_pos;
	cost += alpha_goal*e.squaredNorm();

	arm->set_joint_values(j_orig);


	return cost;

}

double PR2System::cost_gmm(const StdVectorJ& J, const MatrixJ& j_sigma0, const StdVectorU& U,
		const std::vector<ParticleGaussian>& particle_gmm, const double alpha) {
	double cost_gmm = 0;

	MatrixX sigma0 = MatrixX::Zero();
	sigma0.block<J_DIM,J_DIM>(0,0) = j_sigma0;
	for(int i=0; i < particle_gmm.size(); ++i) {
		sigma0.block<3,3>(J_DIM,J_DIM) = particle_gmm[i].cov;
//		cost_gmm += particle_gmm[i].pct*cost(J, particle_gmm[i].mean, sigma0, U, alpha);
		cost_gmm += cost(J, particle_gmm[i].mean, sigma0, U, alpha);
	}

	return cost_gmm;
}

VectorTOTAL PR2System::cost_gmm_grad(StdVectorJ& J, const MatrixJ& j_sigma0, StdVectorU& U,
			const std::vector<ParticleGaussian>& particle_gmm, const double alpha) {
	VectorTOTAL grad;

	double orig, cost_p, cost_m;
	int index = 0;
	for(int t=0; t < TIMESTEPS; ++t) {
		for(int i=0; i < J_DIM; ++i) {
			orig = J[t][i];

			J[t](i) = orig + step;
			cost_p = cost_gmm(J, j_sigma0, U, particle_gmm, alpha);

			J[t](i) = orig - step;
			cost_m = cost_gmm(J, j_sigma0, U, particle_gmm, alpha);

			grad(index++) = (cost_p - cost_m) / (2*step);
		}

		if (t < TIMESTEPS-1) {
			for(int i=0; i < U_DIM; ++i) {
				orig = U[t][i];

				U[t](i) = orig + step;
				cost_p = cost_gmm(J, j_sigma0, U, particle_gmm, alpha);

				U[t](i) = orig - step;
				cost_m = cost_gmm(J, j_sigma0, U, particle_gmm, alpha);

				grad(index++) = (cost_p - cost_m) / (2*step);
			}
		}
	}
	return grad;
}

void PR2System::fit_gaussians_to_pf(const MatrixP& P, std::vector<ParticleGaussian>& particle_gmm) {
	int d = 3;							 // dimension
	int M = M_DIM; 						 // number of targets
	int N = M_DIM; 						 // number of sources
	double h = sqrt(2)*.1; 				 // bandwith (h = sqrt(2)*sigma)
	double eps = .1; 				 // 1e-2
	double *x = new double[d*N]; 		 // source array
	double *y = new double[d*M]; 		 // target array
	int W = 3+1; 						 // 3 for weighted sum, +1 for the normalization
	double *q = new double[W*N]; 		 // weighting array
	double *output = new double[W*M]; 	 // output FGT

	// source and weights always constant

	for(int j=0; j < N; ++j) {
		for(int i=0; i < d; ++i) {
			x[j*d+i] = P(i,j);
		}
	}

	for(int i=0; i < d; ++i) {
		for(int j=0; j < N; ++j) {
			q[i*N+j] = P(i,j);
		}
	}
	for(int j=0; j < N; ++j) {
		q[d*N+j] = 1;
	}

	MatrixP means = P, new_means;
	double max_diff = INFINITY;
	while(max_diff > .001) { // TODO: determines convergence speed
		for(int j=0; j < M; ++j) {
			for(int i=0; i < d; ++i) {
				y[j*d+i] = means(i,j);
			}
		}

		memset(output, 0, sizeof(double)*W*M );

		figtree(d, N, M, W, x, h, q, y, eps, output);

		for(int j=0; j < M; ++j) {
			for(int i=0; i < d; ++i) {
				new_means(i,j) = output[i*M+j] / output[d*M+j];
			}
		}

		max_diff = (means - new_means).colwise().norm().maxCoeff();
		means = new_means;
	}

	StdVector3d modes;
	std::vector<std::vector<int>> mode_particle_indices;
	for(int m=0; m < M_DIM; ++m) {
		bool is_new_mode = true;
		for(int i=0; i < modes.size(); ++i) {
			if ((means.col(m) - modes[i]).norm() < .1) {
				is_new_mode = false;
				mode_particle_indices[i].push_back(m);
				break;
			}
		}

		if (is_new_mode) {
			modes.push_back(means.col(m));
			mode_particle_indices.push_back(std::vector<int>(1, m));
		}
	}

	// create ParticleGaussian vector
	particle_gmm.clear();

	Vector3d mean;
	Matrix3d cov;
	MatrixXd particles;
	for(int i=0; i < mode_particle_indices.size(); ++i) {
		int num_mode_particles = mode_particle_indices[i].size();
		particles = MatrixXd::Zero(3, num_mode_particles);
		for(int m=0; m < num_mode_particles; ++m) {
			particles.col(m) = P.col(mode_particle_indices[i][m]);
		}

		mean = particles.rowwise().mean();

		MatrixXd particles_centered = particles.colwise() - particles.rowwise().mean();
		cov = (1/(double(num_mode_particles)-1))*(particles_centered*particles_centered.transpose());

		particle_gmm.push_back(ParticleGaussian(mean, cov, particles, num_mode_particles/double(M_DIM)));
	}

	// sort with highest pct first
	std::sort(particle_gmm.begin(), particle_gmm.end(),
			[](const ParticleGaussian& a, const ParticleGaussian& b) {
		return a.pct > b.pct;
	});
}

/**
 * PR2System Display methods
 */

void PR2System::display(const VectorJ& j, bool pause) {
	display(StdVectorJ(1,j), pause);
}

void PR2System::display(const StdVectorJ& J, bool pause) {
	display(J, std::vector<ParticleGaussian>(), pause);
}

void PR2System::display(const VectorJ& j, const std::vector<ParticleGaussian>& particle_gmm, bool pause) {
	display(StdVectorJ(1,j), particle_gmm, pause);
}

void PR2System::display(const StdVectorJ& J, const std::vector<ParticleGaussian>& particle_gmm, bool pause) {
	if (pause) {
		rave_utils::clear_plots();
	}

	VectorJ j_orig = arm->get_joint_values();
	for(int t=0; t < J.size(); ++t) {
		arm->set_joint_values(J[t]);
		rave_utils::plot_transform(brett->get_env(), cam->get_pose());
	}
	std::vector<std::vector<Beam3d> > beams = cam->get_beams();
	cam->plot_fov(beams);
	arm->set_joint_values(j_orig);

	int num_gaussians = particle_gmm.size();
	for(int i=0; i < num_gaussians; ++i) {
		double hue = (2/3.0)*(i/double(num_gaussians));
		Vector3d color = utils::hsv_to_rgb({hue, 1, 1});

		rave_utils::plot_point(brett->get_env(), particle_gmm[i].mean, color, .02);
		for(int m=0; m < particle_gmm[i].particles.cols(); ++m) {
			rave_utils::plot_point(brett->get_env(), particle_gmm[i].particles.col(m), color, .005);
		}
	}

	rave_utils::plot_point(brett->get_env(), object, {0,0,1}, .03);

	if (pause) {
		LOG_INFO("Display: Press enter to continue");
		std::cin.ignore();
	}
}

/**
 * PR2System Private methods
 */

void PR2System::linearize_dynfunc(const VectorX& x, const VectorU& u, const VectorQ& q,
		Matrix<double,X_DIM,X_DIM>& A, Matrix<double,X_DIM,Q_DIM>& M) {
	A.setZero();
	M.setZero();

	VectorQ q_zero = VectorQ::Zero();
	VectorX x_p = x, x_m = x;
	for(int i=0; i < X_DIM; ++i) {
		x_p(i) += step; x_m(i) -= step;
		A.block<J_DIM,1>(0, i) = (dynfunc(x_p.segment<J_DIM>(0), u, q_zero) - dynfunc(x_m.segment<J_DIM>(0), u, q_zero)) / (2*step);
		x_p(i) = x(i); x_m(i) = x(i);
	}

	VectorQ q_p = q, q_m = q;
	for(int i=0; i < Q_DIM; ++i) {
		q_p(i) += step; q_m(i) -= step;
		M.block<J_DIM,1>(0, i) = (dynfunc(x.segment<J_DIM>(0), u, q_p) - dynfunc(x.segment<J_DIM>(0), u, q_m)) / (2*step);
	}
}

void PR2System::linearize_obsfunc(const VectorX& x, const VectorR& r,
		Matrix<double,Z_DIM,X_DIM>& H, Matrix<double,Z_DIM,R_DIM>& N) {
	H.setZero();
	N.setZero();

	VectorX x_p = x, x_m = x;
	for(int i=0; i < X_DIM; ++i) {
		x_p(i) += step; x_m(i) -= step;
		H.block<Z_DIM,1>(0, i) = (obsfunc(x_p.segment<J_DIM>(0), x_p.segment<3>(J_DIM), r) -
				obsfunc(x_m.segment<J_DIM>(0), x_m.segment<3>(J_DIM), r)) / (2*step);
		x_p(i) = x(i); x_m(i) = x(i);
	}

	VectorR r_p = r, r_m = r;
	for(int i=0; i < R_DIM; ++i) {
		r_p(i) += step; r_m(i) -= step;
		N.block<Z_DIM,1>(0, i) = (obsfunc(x.segment<J_DIM>(0), x.segment<3>(J_DIM), r_p) -
				obsfunc(x.segment<J_DIM>(0), x.segment<3>(J_DIM), r_m)) / (2*step);
		r_p(i) = x(i); r_m(i) = x(i);
	}
}


void PR2System::update_particles(const VectorJ& j_tp1_t, const double delta_fov_real, const VectorZ& z_tp1_real, const MatrixP& P_t,
		MatrixP& P_tp1) {
	Vector3d z_obj_real = z_tp1_real.segment<3>(J_DIM);

	VectorJ j_orig = arm->get_joint_values();
	arm->set_joint_values(j_tp1_t);
	std::vector<std::vector<Beam3d> > beams = cam->get_beams();

	VectorM W = VectorM::Zero();
	// for each particle, weight by gauss_likelihood of that measurement given particle/agent observation
	for(int m=0; m < M_DIM; ++m) {
		bool inside = cam->is_inside(P_t.col(m), beams);
		if (delta_fov_real < epsilon) {
			W(m) = (inside) ? 0 : 1;
		} else {
			if (inside) {
				Vector3d z_obj_m = obsfunc(j_tp1_t, P_t.col(m), VectorR::Zero()).segment<3>(J_DIM);
				Vector3d e = z_obj_real - z_obj_m;
				W(m) = gauss_likelihood(e, .1*R.block<3,3>(J_DIM,J_DIM));
			} else {
				W(m) = 0;
			}
		}
	}
	W = W / W.sum();

	low_variance_sampler(P_t, W, P_tp1);

	arm->set_joint_values(j_orig);
}

double PR2System::gauss_likelihood(const Vector3d& v, const Matrix3d& S) {
	Matrix3d Sf = S.llt().matrixL();
	Vector3d M = Sf.lu().solve(v);

	double E = -0.5*M.dot(M);
	double C = pow(2*M_PI, S.cols()/2) * Sf.diagonal().prod();
	double w = exp(E) / C;

	return w;
}

void PR2System::low_variance_sampler(const MatrixP& P, const VectorM& W, MatrixP& P_sampled) {
	double r = mm_utils::uniform(0, 1/double(M_DIM));
	double c = W(0);
	int i = 0;
	for(int m=0; m < M_DIM; ++m) {
		double u = r + m * (1/double(M_DIM));
		while (u > c) {
			c += W(++i);
		}
		P_sampled.col(m) = P.col(i);
	}
}
