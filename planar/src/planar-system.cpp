#include "../include/planar-system.h"

/**
 * Constructors and initializers
 */

PlanarSystem::PlanarSystem(const vec<C_DIM>& camera_origin, const vec<C_DIM>& object, bool is_static) {
	init(camera_origin, object, is_static);
	init_display();
}

void PlanarSystem::init(const vec<C_DIM>& camera_origin, const vec<C_DIM>& object, bool is_static) {
	this->camera_origin = camera_origin;
	camera_fov = M_PI/4;
	camera_max_dist = 15;
	this->object = object;
	this->is_static = is_static;

	robot_origin = Vector2d::Zero();
	link_lengths << 5, 5, 4;

	Q = (M_PI/4)*mat<U_DIM,U_DIM>::Identity();
//	R = 10*eye<mat>(Z_DIM, Z_DIM);
	vec<R_DIM> R_diag;
	R_diag << (M_PI/4)*vec<J_DIM>::Ones(),
			5*vec<C_DIM>::Ones();
	R = R_diag.asDiagonal();

	// x bound for object is sum of link_lengths
	double max_link_length = link_lengths.sum();
	x_min << -3*M_PI/2, -3*M_PI/2, -3*M_PI/2, -3*M_PI/2, -max_link_length, -max_link_length;
	x_max << 3*M_PI/2, 3*M_PI/2, 3*M_PI/2, 3*M_PI/2, max_link_length, max_link_length;

	double max_input = M_PI/12;
	u_min << -max_input, -max_input, -max_input, (is_static) ? 0 : -max_input;
	u_max << max_input, max_input, max_input, (is_static) ? 0 : max_input;

//	fps.init(camera_origin.cast<bdouble>(), object.cast<bdouble>(), is_static);
}

void PlanarSystem::init_display() {
	Py_Initialize();
	np::initialize();
	py::numeric::array::set_module_and_type("numpy", "ndarray");

	std::string working_dir = boost::filesystem::current_path().normalize().string();
	std::string bsp_dir = working_dir.substr(0,working_dir.find("bsp"));
	std::string planar_dir = bsp_dir + "bsp/planar";

	py::object main_module = py::import("__main__");
	py::object main_namespace = main_module.attr("__dict__");
	py::exec("import sys, os", main_namespace);
	py::exec(py::str("sys.path.append('"+planar_dir+"')"), main_namespace);
	py::object plot_mod = py::import("plot_planar");

	plot_planar = plot_mod.attr("plot_planar");
	plot_planar_gmm = plot_mod.attr("plot_planar_gmm");
}


/**
 * Public methods
 */

vec<J_DIM> PlanarSystem::dynfunc(const vec<J_DIM>& j, const vec<U_DIM>& u, const vec<Q_DIM>& q, bool enforce_limits) {
	return (j + DT*(u + q));
}

vec<Z_DIM> PlanarSystem::obsfunc(const vec<J_DIM>& j, const vec<C_DIM>& object, const vec<R_DIM>& r) {
	vec<Z_DIM> z;
	z(0) = j(0); // joint 0
	z(1) = j(1); // joint 1
	z(2) = j(2); // joint 2
	z(3) = j(3); // camera angle
	z(4) = object(0) - camera_origin(0); // delta x to object
	z(5) = object(1) - camera_origin(1); // delta y to object
	return z + r;
}


mat<Z_DIM,Z_DIM> PlanarSystem::delta_matrix(const vec<J_DIM>& j, const vec<C_DIM>& object, const double alpha) {
	mat<Z_DIM,Z_DIM> delta = mat<Z_DIM,Z_DIM>::Zero();

	for(int i=0; i < J_DIM; ++i) {
		delta(i, i) = 1; // TODO: should this depend on SD of link segments?
	}

	std::vector<Beam> fov = get_fov(j);
	double sd = geometry2d::signed_distance(object, fov);
	double sd_sigmoid = 1.0 - 1.0/(1.0 + exp(-alpha*sd));
	for(int i=J_DIM; i < X_DIM; ++i) {
		delta(i, i) = sd_sigmoid;
	}

	return delta;
}

/**
 * \brief Propagates belief through EKF with max-likelihood noise assumption
 */
void PlanarSystem::belief_dynamics(const vec<X_DIM>& x_t, const mat<X_DIM,X_DIM>& sigma_t, const vec<U_DIM>& u_t, const double alpha,
		vec<X_DIM>& x_tp1, mat<X_DIM,X_DIM>& sigma_tp1) {
	// propagate dynamics
	x_tp1 = x_t;
	x_tp1.segment<J_DIM>(0) = dynfunc(x_t.segment<J_DIM>(0), u_t, vec<Q_DIM>::Zero());

	// propagate belief through dynamics
	mat<X_DIM,X_DIM> A = mat<X_DIM,X_DIM>::Zero();
	mat<X_DIM,Q_DIM> M = mat<X_DIM,Q_DIM>::Zero();
	linearize_dynfunc(x_t, u_t, vec<Q_DIM>::Zero(), A, M);

	mat<X_DIM,X_DIM> sigma_tp1_bar = A*sigma_t*A.transpose() + M*Q*M.transpose();

	// propagate belief through observation
	mat<Z_DIM,X_DIM> H = mat<Z_DIM,X_DIM>::Zero();
	mat<Z_DIM,R_DIM> N = mat<Z_DIM,R_DIM>::Zero();
	linearize_obsfunc(x_tp1, vec<R_DIM>::Zero(), H, N);

	mat<Z_DIM,Z_DIM> delta = delta_matrix(x_tp1.segment<J_DIM>(0), x_tp1.segment<C_DIM>(J_DIM), alpha);
	mat<X_DIM,Z_DIM> K = sigma_tp1_bar*H.transpose()*delta*(delta*H*sigma_tp1_bar*H.transpose()*delta + R).inverse()*delta;
	sigma_tp1 = (mat<X_DIM,X_DIM>::Identity() - K*H)*sigma_tp1_bar;
}

void PlanarSystem::execute_control_step(const vec<X_DIM>& x_t_real, const vec<X_DIM>& x_t_t, const mat<X_DIM,X_DIM>& sigma_t_t, const vec<U_DIM>& u_t, const mat<C_DIM,M_DIM>& P_t,
			vec<X_DIM>& x_tp1_real, vec<X_DIM>& x_tp1_tp1, mat<X_DIM,X_DIM>& sigma_tp1_tp1, mat<C_DIM,M_DIM>& P_tp1) {
	// find next real state from input + noise
	vec<Q_DIM> control_noise = vec<Q_DIM>::Zero();// + chol(.2*Q)*randn<vec>(Q_DIM);
	vec<R_DIM> obs_noise = vec<R_DIM>::Zero();// + chol(.2*R)*randn<vec>(R_DIM);
	x_tp1_real = x_t_real;
	x_tp1_real.segment<J_DIM>(0) = dynfunc(x_t_real.segment<J_DIM>(0), u_t, control_noise, true);
	vec<Z_DIM> z_tp1_real = obsfunc(x_tp1_real.segment<J_DIM>(0), this->object, obs_noise);

	// now do update based on current belief (x_t, sigma_t)

	// propagate belief through dynamics
	mat<X_DIM,X_DIM> A = mat<X_DIM,X_DIM>::Zero();
	mat<X_DIM,Q_DIM> M = mat<X_DIM,Q_DIM>::Zero();
	linearize_dynfunc(x_t_t, u_t, vec<Q_DIM>::Zero(), A, M);

	mat<X_DIM,X_DIM> sigma_tp1_t = A*sigma_t_t*A.transpose() + M*Q*M.transpose();
	vec<X_DIM> x_tp1_t = x_t_t;
	x_tp1_t.segment<J_DIM>(0) = dynfunc(x_t_t.segment<J_DIM>(0), u_t, vec<Q_DIM>::Zero());

	// propagate belief through observation
	mat<Z_DIM,X_DIM> H = mat<Z_DIM,X_DIM>::Zero();
	mat<Z_DIM,R_DIM> N = mat<Z_DIM,R_DIM>::Zero();
	linearize_obsfunc(x_tp1_t, vec<R_DIM>::Zero(), H, N);

	mat<Z_DIM,Z_DIM> delta = delta_matrix(x_tp1_real.segment<J_DIM>(0), object, INFINITY);
	// calculate Kalman gain
	mat<X_DIM,Z_DIM> K = sigma_tp1_t*H.transpose()*delta*(delta*H*sigma_tp1_t*H.transpose()*delta + R).inverse()*delta;

	// update based on noisy measurement
	x_tp1_tp1 = x_tp1_t + K*(z_tp1_real - obsfunc(x_tp1_t.segment<J_DIM>(0), x_tp1_t.segment<C_DIM>(J_DIM), vec<R_DIM>::Zero()));
	sigma_tp1_tp1 = (mat<X_DIM,X_DIM>::Identity() - K*H)*sigma_tp1_t;

	// and update particle filter
	vec<C_DIM> delta_real = delta.diagonal().segment<C_DIM>(J_DIM);
	update_particles(x_tp1_t.segment<J_DIM>(0), delta_real(0), z_tp1_real, P_t, P_tp1);
}

void PlanarSystem::execute_control_step(const vec<J_DIM>& j_t_real, const vec<J_DIM>& j_t, const vec<U_DIM>& u_t, const mat<C_DIM,M_DIM>& P_t,
		vec<J_DIM>& j_tp1_real, vec<J_DIM>& j_tp1, mat<C_DIM,M_DIM>& P_tp1) {
	// find next real state from input + noise
	vec<Q_DIM> control_noise = vec<Q_DIM>::Zero();// + chol(.2*Q)*randn<vec>(Q_DIM);
	vec<R_DIM> obs_noise = vec<R_DIM>::Zero();// + chol(.2*R)*randn<vec>(R_DIM);
	j_tp1_real = dynfunc(j_t_real, u_t, control_noise, true);
	vec<Z_DIM> z_tp1_real = obsfunc(j_tp1_real, this->object, obs_noise);

	// max-likelihood dynfunc estimate
	j_tp1 = dynfunc(j_t, u_t, vec<Q_DIM>::Zero(), true);

	mat<Z_DIM,Z_DIM> delta = delta_matrix(j_tp1_real, object, INFINITY);
	vec<C_DIM> delta_real = delta.diagonal().segment<C_DIM>(J_DIM);
	update_particles(j_tp1, delta_real(0), z_tp1_real, P_t, P_tp1);
}

std::vector<Beam> PlanarSystem::get_fov(const vec<J_DIM>& j) {
	std::vector<Beam> beams, new_beams;

	// start out with full field of view
	double angle = j(J_DIM-1);
	vec<C_DIM> dir_right, dir_left;
	dir_right << sin(angle+camera_fov/2.0),
			cos(angle+camera_fov/2.0);
	dir_left << sin(angle-camera_fov/2.0),
			cos(angle-camera_fov/2.0);
	beams.push_back(Beam(camera_origin,
			camera_origin + camera_max_dist*dir_right,
			camera_origin + camera_max_dist*dir_left));

//	std::cout << "Untruncated FOV:\n" << beams[0].a.transpose() << "\n" << beams[0].b.transpose() << "\n\n";

	std::vector<Segment> link_segments = get_link_segments(j.segment<E_DIM>(0));

	for(int l=0; l < link_segments.size(); ++l) {
//		std::cout << "\ntruncating for link " << l << "\n";
		new_beams.clear();
		for(int i=0; i < beams.size(); ++i) {
//			std::cout << "beams[" << i << "]:\n" << beams[i].a.transpose() << "\n" << beams[i].b.transpose() << "\n";

			std::vector<Beam> new_beams_i = beams[i].truncate(link_segments[l]);
			for(int j=0; j < new_beams_i.size(); ++j) {
				// prevent small beams from being added
				if (new_beams_i[j].top_length() > epsilon) {
					new_beams.push_back(new_beams_i[j]);
//					std::cout << "new_beam:\n" << new_beams.back().a.transpose() << "\n" << new_beams.back().b.transpose() << "\n";
				}
			}
		}
		beams = new_beams;
//
//		std::cout << "\nbeams after truncation of link " << l << "\n";
//		for(int i=0; i < beams.size(); ++i) {
//			std::cout << "beams[" << i << "]:\n" << beams[i].a.transpose() << "\n" << beams[i].b.transpose() << "\n";
//		}
	}

	return beams;
}

std::vector<Segment> PlanarSystem::get_link_segments(const vec<E_DIM>& j) {
	std::vector<Segment> link_segments;

	vec<C_DIM> p0 = robot_origin, p1;
	double joint0 = 0, joint1;

	for(int i=0; i < link_lengths.rows(); ++i) {
		joint1 = joint0 + j[i];
		p1 << p0(0) + sin(joint1)*link_lengths(i),
				p0(1) + cos(joint1)*link_lengths(i);
		link_segments.push_back(Segment(p0, p1));

		joint0 = joint1;
		p0 = p1;
	}

	return link_segments;
}


vec<C_DIM> PlanarSystem::get_ee_pos(const vec<E_DIM>& j) {
	return get_link_segments(j).back().p1;
}

void PlanarSystem::get_ee_pos_jac(vec<E_DIM>& j, mat<C_DIM,E_DIM>& ee_jac) {
	vec<E_DIM> j_p = j, j_m = j;
	for(int i=0; i < E_DIM; ++i) {
		j_p(i) = j(i) + step;
		j_m(i) = j(i) - step;
		ee_jac.block<C_DIM,1>(0, i) = (get_ee_pos(j_p) - get_ee_pos(j_m)) / (2*step);
		j_p(i) = j(i);
		j_m(i) = j(i);
	}
}

bool PlanarSystem::ik(const vec<C_DIM>& ee_goal, vec<E_DIM>& j) {
	vec<C_DIM> ee_pos, error, jac_jacT_error;
	mat<C_DIM,E_DIM> ee_pos_jac;

	int iter = 0, max_iter = 1000;
	while((ee_goal - get_ee_pos(j)).norm() > epsilon) {
		ee_pos = get_ee_pos(j);
		get_ee_pos_jac(j, ee_pos_jac);

		error = ee_goal - ee_pos;

//		// x += alpha*jac.T*error
		jac_jacT_error = ee_pos_jac*ee_pos_jac.transpose()*error;
		double alpha = error.dot(jac_jacT_error) / jac_jacT_error.dot(jac_jacT_error);
		j += alpha*ee_pos_jac.transpose()*error;

		if (iter++ > max_iter) {
			LOG_ERROR("IK failed: reached max iterations");
			return false;
		}
	}

	// check limits
	for(int i=0; i < E_DIM; ++i) {
		if ((j(i) < x_min(i)) || (j(i) > x_max(i))) {
			LOG_ERROR("IK failed: solution outside of limits")
			return false;
		}
	}

	return true;
}

void PlanarSystem::get_limits(vec<X_DIM>& x_min, vec<X_DIM>& x_max, vec<U_DIM>& u_min, vec<U_DIM>& u_max) {
	x_min = this->x_min;
	x_max = this->x_max;
	u_min = this->u_min;
	u_max = this->u_max;
}

double PlanarSystem::cost(const std::vector<vec<J_DIM>, aligned_allocator<vec<J_DIM>>>& J, const vec<C_DIM>& obj,
		const mat<X_DIM,X_DIM>& sigma0, const std::vector<vec<U_DIM>, aligned_allocator<vec<U_DIM>>>& U, const double alpha) {
	double cost = 0;
	int T = J.size();

	vec<X_DIM> x_t, x_tp1 = vec<X_DIM>::Zero();
	mat<X_DIM,X_DIM> sigma_t = sigma0, sigma_tp1 = mat<X_DIM,X_DIM>::Zero();
	for(int t=0; t < T-1; ++t) {
		x_t << J[t], obj;
		belief_dynamics(x_t, sigma_t, U[t], alpha, x_tp1, sigma_tp1);
		cost += alpha_control*U[t].dot(U[t]);
		if (t < T-2) {
			cost += alpha_belief*sigma_tp1.trace();
		} else {
			cost += alpha_final_belief*sigma_tp1.trace();
		}
		sigma_t = sigma_tp1;
	}

	vec<C_DIM> final_pos = get_ee_pos(x_tp1.segment<E_DIM>(0));
	vec<C_DIM> e = obj - final_pos;
	cost += alpha_goal*e.squaredNorm();

//	cost += .1*(atan((obj(0) - camera_origin(0))/(obj(1) - camera_origin(1))) - x_tp1(E_DIM));

	return cost;
}

/**
 * \brief Propagate each Gaussian in planar_gmm through cost function
 * and add weighted sum
 */
double PlanarSystem::cost_gmm(const std::vector<vec<J_DIM>, aligned_allocator<vec<J_DIM>>>& J, const mat<J_DIM,J_DIM>& j_sigma0,
			const std::vector<vec<U_DIM>, aligned_allocator<vec<U_DIM>>>& U,
			const std::vector<PlanarGaussian>& planar_gmm, const double alpha) {
	double cost_gmm = 0;

	mat<X_DIM,X_DIM> sigma0 = mat<X_DIM,X_DIM>::Zero();
	sigma0.block<J_DIM,J_DIM>(0,0) = j_sigma0;
	for(int i=0; i < planar_gmm.size(); ++i) {
		sigma0.block<C_DIM,C_DIM>(J_DIM,J_DIM) = planar_gmm[i].obj_cov;
//		cost_gmm += planar_gmm[i].pct*cost(J, planar_gmm[i].obj_mean, sigma0, U, alpha);
		cost_gmm += cost(J, planar_gmm[i].obj_mean, sigma0, U, alpha);
	}

	return cost_gmm;
}

double PlanarSystem::cost_entropy(const std::vector<vec<J_DIM>, aligned_allocator<vec<J_DIM>>>& J,
		const std::vector<vec<U_DIM>, aligned_allocator<vec<U_DIM>>>& U,
		const mat<C_DIM,M_DIM>& P, const double alpha) {
	int T = J.size();
	double entropy = 0;
#define H_DIM 3

	vec<J_DIM> j_tp1;
	std::vector<mat<H_DIM,M_DIM>, aligned_allocator<mat<H_DIM,M_DIM>>> H(T);
	for(int t=0; t < T-1; ++t) {
		j_tp1 = dynfunc(J[t], U[t], vec<Q_DIM>::Zero());
		std::vector<Beam> fov = get_fov(j_tp1);
		for(int m=0; m < M_DIM; ++m) {
			double sd = geometry2d::signed_distance(P.col(m), fov);
			double delta = 1.0 - 1.0/(1.0 + exp(-alpha*sd));
			H[t+1].col(m) << delta, P.col(m) - camera_origin;
		}
	}

	vec<H_DIM> S_diag;
	S_diag << .01, 1, 1;
	mat<H_DIM,H_DIM> S = S_diag.asDiagonal();


	vec<M_DIM> W_t = (1/double(M_DIM))*vec<M_DIM>::Ones(), W_tp1;
	for(int t=1; t < T; ++t) {
		W_tp1.setZero();
		for(int m=0; m < M_DIM; ++m) {
			for(int p=0; p < M_DIM; ++p) {
				vec<H_DIM> v = H[t].col(m) - H[t].col(p);

				/* Gauss likelihood */
				mat<H_DIM,H_DIM> Sf = S.llt().matrixL();
				vec<H_DIM> M = Sf.lu().solve(v);

				double E = -0.5*M.dot(M);
				double C = pow(2*M_PI, S.cols()/2) * Sf.diagonal().prod();
				double w = exp(E) / C;
				/* end */

				W_tp1(m) += w;
			}
			W_tp1(m) *= W_t(m);
		}
		W_tp1 /= W_tp1.sum();

		entropy += (-W_tp1.array() * W_tp1.array().log()).sum();

		W_t = W_tp1;
	}

	return entropy;
}

vec<TOTAL_VARS> PlanarSystem::cost_grad(std::vector<vec<J_DIM>, aligned_allocator<vec<J_DIM>>>& J, const vec<C_DIM>& obj,
		const mat<X_DIM,X_DIM>& sigma0, std::vector<vec<U_DIM>, aligned_allocator<vec<U_DIM>>>& U, const double alpha) {
	int T = J.size();

	vec<TOTAL_VARS> grad;

	double orig, cost_p, cost_m;
	int index = 0;
	for(int t=0; t < T; ++t) {
		for(int i=0; i < J_DIM; ++i) {
			orig = J[t][i];

			J[t][i] = orig + step;
			cost_p = cost(J, obj, sigma0, U, alpha);

			J[t][i] = orig - step;
			cost_m = cost(J, obj, sigma0, U, alpha);

			grad(index++) = (cost_p - cost_m) / (2*step);
		}

		if (t < T-1) {
			for(int i=0; i < U_DIM; ++i) {
				orig = U[t][i];

				U[t][i] = orig + step;
				cost_p = cost(J, obj, sigma0, U, alpha);

				U[t][i] = orig - step;
				cost_m = cost(J, obj, sigma0, U, alpha);

				grad(index++) = (cost_p - cost_m) / (2*step);
			}
		}
	}
	return grad;
}

vec<TOTAL_VARS> PlanarSystem::cost_gmm_grad(std::vector<vec<J_DIM>, aligned_allocator<vec<J_DIM>>>& J, const mat<J_DIM,J_DIM>& j_sigma0,
		std::vector<vec<U_DIM>, aligned_allocator<vec<U_DIM>>>& U,
		const std::vector<PlanarGaussian>& planar_gmm, const double alpha) {
	int T = J.size();

	vec<TOTAL_VARS> grad;

	double orig, cost_p, cost_m;
	int index = 0;
	for(int t=0; t < T; ++t) {
		for(int i=0; i < J_DIM; ++i) {
			orig = J[t][i];

			J[t][i] = orig + step;
			cost_p = cost_gmm(J, j_sigma0, U, planar_gmm, alpha);

			J[t][i] = orig - step;
			cost_m = cost_gmm(J, j_sigma0, U, planar_gmm, alpha);

			grad(index++) = (cost_p - cost_m) / (2*step);
		}

		if (t < T-1) {
			for(int i=0; i < U_DIM; ++i) {
				orig = U[t][i];

				U[t][i] = orig + step;
				cost_p = cost_gmm(J, j_sigma0, U, planar_gmm, alpha);

				U[t][i] = orig - step;
				cost_m = cost_gmm(J, j_sigma0, U, planar_gmm, alpha);

				grad(index++) = (cost_p - cost_m) / (2*step);
			}
		}
	}
	return grad;
}

vec<TOTAL_VARS> PlanarSystem::cost_entropy_grad(std::vector<vec<J_DIM>, aligned_allocator<vec<J_DIM>>>& J,
		std::vector<vec<U_DIM>, aligned_allocator<vec<U_DIM>>>& U,
		const mat<C_DIM,M_DIM>& P, const double alpha) {
	int T = J.size();

	vec<TOTAL_VARS> grad;

	double orig, cost_p, cost_m;
	int index = 0;
	for(int t=0; t < T; ++t) {
		for(int i=0; i < J_DIM; ++i) {
			orig = J[t][i];

			J[t][i] = orig + step;
			cost_p = cost_entropy(J, U, P, alpha);

			J[t][i] = orig - step;
			cost_m = cost_entropy(J, U, P, alpha);

			grad(index++) = (cost_p - cost_m) / (2*step);
		}

		if (t < T-1) {
			for(int i=0; i < U_DIM; ++i) {
				orig = U[t][i];

				U[t][i] = orig + step;
				cost_p = cost_entropy(J, U, P, alpha);

				U[t][i] = orig - step;
				cost_m = cost_entropy(J, U, P, alpha);

				grad(index++) = (cost_p - cost_m) / (2*step);
			}
		}
	}
	return grad;

}

//void PlanarSystem::cost_and_cost_grad(std::vector<vec<X_DIM>, aligned_allocator<vec<X_DIM>>>& X, const mat<X_DIM,X_DIM>& sigma0,
//		std::vector<vec<U_DIM>, aligned_allocator<vec<U_DIM>>>& U, const double alpha, const bool use_fadbad,
//		double& cost, vec<TOTAL_VARS>& grad) {
//	if (!use_fadbad) {
//		cost = this->cost(X, sigma0, U, alpha);
//		grad = this->cost_grad(X, sigma0, U, alpha);
//	} else {
//		std::vector<vecb<X_DIM>, aligned_allocator<vecb<X_DIM>>> X_b(X.size());
//		matb<X_DIM,X_DIM> sigma0_b = sigma0.cast<bdouble>();
//		std::vector<vecb<U_DIM>, aligned_allocator<vecb<U_DIM>>> U_b(U.size());
//		bdouble alpha_b = (bdouble)alpha;
//
//		for(int t=0; t < X.size(); ++t) {
//			X_b[t] = X[t].cast<bdouble>();
//			if(t < X.size()-1) {
//				U_b[t] = U[t].cast<bdouble>();
//			}
//		}
//
//		bdouble cost_b = fps.cost(X_b, sigma0_b, U_b, alpha_b);
//		cost = cost_b.x();
//		cost_b.diff(0,1);
//		int index = 0;
//		for(int t=0; t < X.size(); ++t) {
//			for(int i=0; i < X_DIM; ++i) {
//				grad(index++) = X_b[t](i).d(0);
//			}
//
//			if (t < X.size()-1) {
//				for(int i=0; i < U_DIM; ++i) {
//					grad(index++) = U_b[t](i).d(0);
//				}
//			}
//		}
//	}
//}

/**
 * First index of planar_gmm contains the most particles
 */
void PlanarSystem::fit_gaussians_to_pf(const mat<C_DIM,M_DIM>& P,
		std::vector<PlanarGaussian>& planar_gmm) {
	planar_gmm.clear();

	// fit Gaussians
	std::vector<VectorXd> obj_means_dyn;
	std::vector<MatrixXd> obj_covs_dyn;
	std::vector<MatrixXd> obj_particles_tmp;
	gmm::fit_gaussians_to_pf(P, obj_means_dyn, obj_covs_dyn, obj_particles_tmp);

	// find Gaussian with most particles
	int max_obj_index = -1;
	int max_obj_num_particles = -1;
	for(int i=0; i < obj_particles_tmp.size(); ++i) {
		if (obj_particles_tmp[i].cols() > max_obj_num_particles) {
			max_obj_index = i;
			max_obj_num_particles = obj_particles_tmp[i].cols();
		}
	}

	// create structs and push
	vec<C_DIM> obj_mean;
	mat<C_DIM,C_DIM> obj_cov;

	obj_mean = obj_means_dyn[max_obj_index];
	obj_cov = obj_covs_dyn[max_obj_index];
	planar_gmm.push_back(PlanarGaussian(obj_mean, obj_cov,
			obj_particles_tmp[max_obj_index], obj_particles_tmp[max_obj_index].cols()/double(M_DIM)));

	for(int i=0; i < obj_means_dyn.size(); ++i) {
		int num_particles = obj_particles_tmp[i].cols();
		if ((i != max_obj_index) && (num_particles > 1)) {
			obj_mean = obj_means_dyn[i];
			obj_cov = obj_covs_dyn[i];
			planar_gmm.push_back(PlanarGaussian(obj_mean, obj_cov,
					obj_particles_tmp[i], num_particles/double(M_DIM)));
		}
	}
}

void PlanarSystem::fit_gaussians_to_pf_figtree(const mat<C_DIM,M_DIM>& P,
		std::vector<PlanarGaussian>& planar_gmm) {
	util::Timer figtree_timer;
	util::Timer_tic(&figtree_timer);

	int d = C_DIM;						 // dimension
	int M = M_DIM; 						 // number of targets
	int N = M_DIM; 						 // number of sources
	double h = sqrt(2)*3; 				 // bandwith (h = sqrt(2)*sigma)
	double epsilon = .1; 				 // 1e-2
	double *x = new double[d*N]; 		 // source array
	double *y = new double[d*M]; 		 // target array
	int W = C_DIM+1; 					 // C_DIM for weighted sum, +1 for the normalization
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

	mat<C_DIM,M_DIM> means = P, new_means;
	double max_diff = INFINITY;
	while(max_diff > 1e-3) {
		for(int j=0; j < M; ++j) {
			for(int i=0; i < d; ++i) {
				y[j*d+i] = means(i,j);
			}
		}

		memset(output, 0, sizeof(double)*W*M );

		figtree(d, N, M, W, x, h, q, y, epsilon, output);

		for(int j=0; j < M; ++j) {
			for(int i=0; i < d; ++i) {
				new_means(i,j) = output[i*M+j] / output[d*M+j];
			}
		}

		max_diff = (means - new_means).colwise().norm().maxCoeff();
		means = new_means;
	}

	std::vector<vec<C_DIM>, aligned_allocator<vec<C_DIM>>> modes;
	std::vector<std::vector<int>> mode_particle_indices;
	for(int m=0; m < M_DIM; ++m) {
		bool is_new_mode = true;
		for(int i=0; i < modes.size(); ++i) {
			if ((means.col(m) - modes[i]).norm() < .05) {
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

	// create PlanarGaussian vector
	planar_gmm.clear();

	vec<C_DIM> obj_mean;
	mat<C_DIM,C_DIM> obj_cov;
	MatrixXd obj_particles;
	for(int i=0; i < mode_particle_indices.size(); ++i) {
		int num_mode_particles = mode_particle_indices[i].size();
		obj_particles = MatrixXd::Zero(C_DIM, num_mode_particles);
		for(int m=0; m < num_mode_particles; ++m) {
			obj_particles.col(m) = P.col(mode_particle_indices[i][m]);
		}

		obj_mean = obj_particles.rowwise().mean();

		MatrixXd obj_particles_centered = obj_particles.colwise() - obj_particles.rowwise().mean();
		obj_cov = (1/(double(num_mode_particles)-1))*(obj_particles_centered*obj_particles_centered.transpose());

		planar_gmm.push_back(PlanarGaussian(obj_mean, obj_cov, obj_particles, num_mode_particles/double(M_DIM)));
	}

	// sort with highest pct first
	std::sort(planar_gmm.begin(), planar_gmm.end(),
	          [](const PlanarGaussian& a, const PlanarGaussian& b) {
	  return a.pct > b.pct;
	});


	double figtree_time = util::Timer_toc(&figtree_timer);
	LOG_INFO("Figtree time: %4.4f", figtree_time);
}

/**
 * Display methods
 */

void PlanarSystem::display(const vec<J_DIM>& j, bool pause) {
	std::vector<vec<J_DIM>, aligned_allocator<vec<J_DIM>>> J(1, j);
	display(J, pause);
}

void PlanarSystem::display(const std::vector<vec<J_DIM>, aligned_allocator<vec<J_DIM>>>& J, bool pause) {
	py::list J_pylist;
	for(int t=0; t < J.size(); ++t) {
		J_pylist.append(planar_utils::eigen_to_ndarray(J[t]));
	}

	np::ndarray robot_origin_ndarray = planar_utils::eigen_to_ndarray(robot_origin);
	np::ndarray link_lengths_ndarray = planar_utils::eigen_to_ndarray(link_lengths);
	np::ndarray camera_origin_ndarray = planar_utils::eigen_to_ndarray(camera_origin);
	np::ndarray object_ndarray = planar_utils::eigen_to_ndarray(object);

	py::list beams_pylist;
	std::vector<Beam> beams = get_fov(J.back());

	for(int i=0; i < beams.size(); ++i) {
		mat<2,3> m;
		m << beams[i].base, beams[i].a, beams[i].b;
		beams_pylist.append(planar_utils::eigen_to_ndarray(m));
	}

	try {
		plot_planar(J_pylist, robot_origin_ndarray, link_lengths_ndarray, camera_origin_ndarray, camera_fov, object_ndarray, beams_pylist, pause);
	}
	catch(py::error_already_set const &)
	{
		PyErr_Print();
	}
}

void PlanarSystem::display(const vec<J_DIM>& j,
		const std::vector<PlanarGaussian>& planar_gmm,
		bool pause) {
	std::vector<vec<J_DIM>, aligned_allocator<vec<J_DIM>>> J(1, j);
	display(J, planar_gmm, pause);
}

void PlanarSystem::display(const std::vector<vec<J_DIM>, aligned_allocator<vec<J_DIM>>>& J,
		const std::vector<PlanarGaussian>& planar_gmm,
		bool pause) {
	py::list J_pylist;
	for(int t=0; t < J.size(); ++t) {
		J_pylist.append(planar_utils::eigen_to_ndarray(J[t]));
	}

	np::ndarray robot_origin_ndarray = planar_utils::eigen_to_ndarray(robot_origin);
	np::ndarray link_lengths_ndarray = planar_utils::eigen_to_ndarray(link_lengths);
	np::ndarray camera_origin_ndarray = planar_utils::eigen_to_ndarray(camera_origin);
	np::ndarray object_ndarray = planar_utils::eigen_to_ndarray(object);

	py::list beams_pylist;
	std::vector<Beam> beams = get_fov(J.back());

	for(int i=0; i < beams.size(); ++i) {
		mat<2,3> m;
		m << beams[i].base, beams[i].a, beams[i].b;
		beams_pylist.append(planar_utils::eigen_to_ndarray(m));
	}

	py::list obj_means_list, obj_covs_list, obj_particles_list;
	for(int i=0; i < planar_gmm.size(); ++i) {
		obj_means_list.append(planar_utils::eigen_to_ndarray(planar_gmm[i].obj_mean));
		obj_covs_list.append(planar_utils::eigen_to_ndarray(planar_gmm[i].obj_cov));
		obj_particles_list.append(planar_utils::eigen_to_ndarray(planar_gmm[i].obj_particles));
	}

	try {
		plot_planar_gmm(J_pylist, obj_means_list, obj_covs_list, obj_particles_list,
				robot_origin_ndarray, link_lengths_ndarray, camera_origin_ndarray, camera_fov, object_ndarray, beams_pylist, pause);
	}
	catch(py::error_already_set const &)
	{
		PyErr_Print();
	}
}



//void PlanarSystem::display(const vec<X_DIM>& x, const mat<X_DIM,X_DIM>& sigma, bool pause) {
//	mat<C_DIM,M_DIM> P;
//	display(x, sigma, P, pause, false);
//}
//
//void PlanarSystem::display(const std::vector<vec<X_DIM>, aligned_allocator<vec<X_DIM>>>& X, const mat<X_DIM,X_DIM>& sigma0,
//		const std::vector<vec<U_DIM>, aligned_allocator<vec<U_DIM>>>& U, const double alpha, bool pause) {
//	mat<C_DIM,M_DIM> P;
//	display(X, sigma0, U, P, alpha, pause, false);
//}
//
//void PlanarSystem::display(const std::vector<vec<X_DIM>, aligned_allocator<vec<X_DIM>>>& X, const std::vector<mat<X_DIM,X_DIM>, aligned_allocator<mat<X_DIM,X_DIM>>>& S, bool pause) {
//	mat<C_DIM,M_DIM> P;
//	display(X, S, P, pause, false);
//}
//
//void PlanarSystem::display(const vec<X_DIM>& x, const mat<X_DIM,X_DIM>& sigma, const mat<C_DIM,M_DIM>& P, bool pause, bool plot_particles) {
//	std::vector<vec<X_DIM>, aligned_allocator<vec<X_DIM>>> X(1, x);
//	std::vector<mat<X_DIM,X_DIM>, aligned_allocator<mat<X_DIM,X_DIM>>> S(1, sigma);
//	display(X, S, P, pause, plot_particles);
//}
//
//void PlanarSystem::display(const std::vector<vec<X_DIM>, aligned_allocator<vec<X_DIM>>>& X,
//		const mat<X_DIM,X_DIM>& sigma0, const std::vector<vec<U_DIM>, aligned_allocator<vec<U_DIM>>>& U, const mat<C_DIM,M_DIM>& P, const double alpha, bool pause, bool plot_particles) {
//	int T = X.size();
//
//	std::vector<mat<X_DIM,X_DIM>, aligned_allocator<mat<X_DIM,X_DIM>>> S(T);
//	S[0] = sigma0;
//	vec<X_DIM> x_tp1;
//	for(int t=0; t < T-1; ++t) {
//		belief_dynamics(X[t], S[t], U[t], alpha, x_tp1, S[t+1]);
//	}
//	display(X, S, P, pause, plot_particles);
//}
//
//void PlanarSystem::display(const std::vector<vec<X_DIM>, aligned_allocator<vec<X_DIM>>>& X,
//		const std::vector<mat<X_DIM,X_DIM>, aligned_allocator<mat<X_DIM,X_DIM>>>& S, const mat<C_DIM,M_DIM>& P, bool pause, bool plot_particles) {
//	py::list X_pylist, S_pylist;
//	for(int t=0; t < X.size(); ++t) {
//		X_pylist.append(planar_utils::eigen_to_ndarray(X[t]));
//		S_pylist.append(planar_utils::eigen_to_ndarray(S[t]));
//	}
//
//	np::ndarray robot_origin_ndarray = planar_utils::eigen_to_ndarray(robot_origin);
//	np::ndarray link_lengths_ndarray = planar_utils::eigen_to_ndarray(link_lengths);
//	np::ndarray camera_origin_ndarray = planar_utils::eigen_to_ndarray(camera_origin);
//	np::ndarray object_ndarray = planar_utils::eigen_to_ndarray(object);
//
//	py::list beams_pylist;
//	std::vector<Beam> beams = get_fov(X.back().segment<J_DIM>(0));
//
//	for(int i=0; i < beams.size(); ++i) {
//		mat<2,3> m;
//		m << beams[i].base, beams[i].a, beams[i].b;
//		beams_pylist.append(planar_utils::eigen_to_ndarray(m));
//	}
//
//	try {
//		if (plot_particles) {
//			np::ndarray P_ndarray = planar_utils::eigen_to_ndarray(P);
//			vec<M_DIM> P_sd;
//			for(int m=0; m < M_DIM; ++m) {
//				P_sd(m) = geometry2d::signed_distance(P.col(m), beams);
//			}
//			np::ndarray P_sd_ndarray = planar_utils::eigen_to_ndarray(P_sd);
//
//			plot_planar(X_pylist, S_pylist, P_ndarray, P_sd_ndarray, robot_origin_ndarray, link_lengths_ndarray,
//							camera_origin_ndarray, camera_fov, object_ndarray, beams_pylist, pause);
//		} else {
//			plot_planar(X_pylist, S_pylist, NULL, NULL, robot_origin_ndarray, link_lengths_ndarray,
//					camera_origin_ndarray, camera_fov, object_ndarray, beams_pylist, pause);
//		}
//	}
//	catch(py::error_already_set const &)
//	{
//		PyErr_Print();
//	}
//}

/**
 * Private methods
 */

void PlanarSystem::linearize_dynfunc(const vec<X_DIM>& x, const vec<U_DIM>& u, const vec<Q_DIM>& q,
		mat<X_DIM,X_DIM>& A, mat<X_DIM,Q_DIM>& M) {
	A = mat<X_DIM,X_DIM>::Identity();

	for(int j=0; j < Q_DIM; ++j) {
		for(int i=0; i < j+1; ++i) {
			M(i, j) = 1;
		}
	}

//	VectorXd q_zero = VectorXd::Zero(Q_DIM);
//	VectorXd x_p = x, x_m = x;
//	for(int i=0; i < X_DIM; ++i) {
//		x_p(i) += step; x_m(i) -= step;
//		A.block<X_DIM,1>(0, i) = (dynfunc(x_p, u, q_zero) - dynfunc(x_m, u, q_zero)) / (2*step);
//		x_p(i) = x(i); x_m(i) = x(i);
//	}
//
//	VectorXd q_p = q, q_m = q;
//	for(int i=0; i < Q_DIM; ++i) {
//		q_p(i) += step; q_m(i) -= step;
//		M.block<X_DIM,1>(0, i) = (dynfunc(x, u, q_p) - dynfunc(x, u, q_m)) / (2*step);
//	}
}

void PlanarSystem::linearize_obsfunc(const vec<X_DIM>& x, const vec<R_DIM>& r, mat<Z_DIM,X_DIM>& H, mat<Z_DIM,R_DIM>& N) {
	H = mat<Z_DIM,X_DIM>::Identity();
	N = mat<Z_DIM,R_DIM>::Identity();

//	VectorXd x_p = x, x_m = x;
//	for(int i=0; i < X_DIM; ++i) {
//		x_p(i) += step; x_m(i) -= step;
//		H.block<Z_DIM,1>(0, i) = (obsfunc(x_p, x_p.segment<C_DIM>(J_DIM), r) -
//				obsfunc(x_m, x_m.segment<C_DIM>(J_DIM), r)) / (2*step);
//		x_p(i) = x(i); x_m(i) = x(i);
//	}
//
//	VectorXd r_p = r, r_m = r;
//	for(int i=0; i < R_DIM; ++i) {
//		r_p(i) += step; r_m(i) -= step;
//		N.block<Z_DIM,1>(0, i) = (obsfunc(x, x.segment<C_DIM>(J_DIM), r_p) -
//				obsfunc(x, x.segment<C_DIM>(J_DIM), r_m)) / (2*step);
//		r_p(i) = x(i); r_m(i) = x(i);
//	}
}

void PlanarSystem::update_particles(const vec<J_DIM>& j_tp1_t, const double delta_fov_real, const vec<Z_DIM>& z_tp1_real, const mat<C_DIM,M_DIM>& P_t,
		mat<C_DIM,M_DIM>& P_tp1) {

	vec<C_DIM> z_obj_real = z_tp1_real.segment<C_DIM>(J_DIM);
	std::vector<Beam> fov = get_fov(j_tp1_t);

	vec<M_DIM> W = vec<M_DIM>::Zero();
	// for each particle, weight by gauss_likelihood of that measurement given particle/agent observation
	for(int m=0; m < M_DIM; ++m) {
		bool inside = geometry2d::is_inside(P_t.col(m), fov);
		if (delta_fov_real < epsilon) {
			W(m) = (inside) ? 0 : 1;
		} else {
			if (inside) {
				vec<C_DIM> z_obj_m = obsfunc(j_tp1_t, P_t.col(m), vec<R_DIM>::Zero()).segment<C_DIM>(J_DIM);
				vec<C_DIM> e = z_obj_real - z_obj_m;
				W(m) = gauss_likelihood(e, .1*R.block<C_DIM,C_DIM>(J_DIM,J_DIM));
			} else {
				W(m) = 0;
			}
		}
	}
	W = W / W.sum();

	low_variance_sampler(P_t, W, P_tp1);
}


double PlanarSystem::gauss_likelihood(const vec<C_DIM>& v, const mat<C_DIM,C_DIM>& S) {
	mat<C_DIM,C_DIM> Sf = S.llt().matrixL();
	vec<C_DIM> M = Sf.lu().solve(v);

	double E = -0.5*M.dot(M);
	double C = pow(2*M_PI, S.cols()/2) * Sf.diagonal().prod();
	double w = exp(E) / C;

	return w;
}

void PlanarSystem::low_variance_sampler(const mat<C_DIM,M_DIM>& P, const vec<M_DIM>& W, mat<C_DIM,M_DIM>& P_sampled) {
	double r = planar_utils::uniform(0, 1/double(M_DIM));
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

