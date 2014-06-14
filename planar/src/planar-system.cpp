#include "../include/planar-system.h"

/**
 * Constructors and initializers
 */

PlanarSystem::PlanarSystem(const vec<C_DIM>& camera_origin, const vec<C_DIM>& object, bool is_static) {
	init(camera_origin, object, is_static);
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
	x_min << -3*M_PI/2, -3*M_PI/2, -3*M_PI/2, -M_PI/2, -max_link_length, -max_link_length;
	x_max << 3*M_PI/2, 3*M_PI/2, 3*M_PI/2, M_PI/2, max_link_length, max_link_length;

	double max_input = M_PI/12;
	u_min << -max_input, -max_input, -max_input, (is_static) ? 0 : -max_input;
	u_max << max_input, max_input, max_input, (is_static) ? 0 : max_input;

	fps.init(camera_origin.cast<bdouble>(), object.cast<bdouble>(), is_static);
}


/**
 * Public methods
 */

vec<X_DIM> PlanarSystem::dynfunc(const vec<X_DIM>& x, const vec<U_DIM>& u, const vec<Q_DIM>& q, bool enforce_limits) {
	vec<X_DIM> x_new = x;
	x_new.segment<J_DIM>(0) += DT*(u + q);
	return x_new;
}

vec<Z_DIM> PlanarSystem::obsfunc(const vec<X_DIM>& x, const vec<C_DIM>& object, const vec<R_DIM>& r) {
	vec<Z_DIM> z;
	z(0) = x(0); // joint 0
	z(1) = x(1); // joint 1
	z(2) = x(2); // joint 2
	z(3) = x(3); // camera angle
	z(4) = object(0) - camera_origin(0); // delta x to object
	z(5) = object(1) - camera_origin(1); // delta y to object
	return z + r;
}


mat<Z_DIM,Z_DIM> PlanarSystem::delta_matrix(const vec<X_DIM>& x, const vec<C_DIM>& object, const double alpha) {
	mat<Z_DIM,Z_DIM> delta = mat<Z_DIM,Z_DIM>::Zero();

	for(int i=0; i < J_DIM; ++i) {
		delta(i, i) = 1; // TODO: should this depend on SD of link segments?
	}

	std::vector<Beam> fov = get_fov(x);
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
	x_tp1 = dynfunc(x_t, u_t, vec<Q_DIM>::Zero());

	// propagate belief through dynamics
	mat<X_DIM,X_DIM> A = mat<X_DIM,X_DIM>::Zero();
	mat<X_DIM,Q_DIM> M = mat<X_DIM,Q_DIM>::Zero();
	linearize_dynfunc(x_t, u_t, vec<Q_DIM>::Zero(), A, M);

	mat<X_DIM,X_DIM> sigma_tp1_bar = A*sigma_t*A.transpose() + M*Q*M.transpose();

	// propagate belief through observation
	mat<Z_DIM,X_DIM> H = mat<Z_DIM,X_DIM>::Zero();
	mat<Z_DIM,R_DIM> N = mat<Z_DIM,R_DIM>::Zero();
	linearize_obsfunc(x_tp1, vec<R_DIM>::Zero(), H, N);

	mat<Z_DIM,Z_DIM> delta = delta_matrix(x_tp1, x_tp1.segment<C_DIM>(J_DIM), alpha);
	mat<X_DIM,Z_DIM> K = sigma_tp1_bar*H.transpose()*delta*(delta*H*sigma_tp1_bar*H.transpose()*delta + R).inverse()*delta;
	sigma_tp1 = (mat<X_DIM,X_DIM>::Identity() - K*H)*sigma_tp1_bar;
}

void PlanarSystem::execute_control_step(const vec<X_DIM>& x_t_real, const vec<X_DIM>& x_t_t, const mat<X_DIM,X_DIM>& sigma_t_t, const vec<U_DIM>& u_t,
		vec<X_DIM>& x_tp1_real, vec<X_DIM>& x_tp1_tp1, mat<X_DIM,X_DIM>& sigma_tp1_tp1) {
	// find next real state from input + noise
	vec<Q_DIM> control_noise = vec<Q_DIM>::Zero();// + chol(.2*Q)*randn<vec>(Q_DIM);
	vec<R_DIM> obs_noise = vec<R_DIM>::Zero();// + chol(.2*R)*randn<vec>(R_DIM);
	x_tp1_real = dynfunc(x_t_real, u_t, control_noise, true);
	vec<Z_DIM> z_tp1_real = obsfunc(x_tp1_real, this->object, obs_noise);

	// now do update based on current belief (x_t, sigma_t)

	// propagate belief through dynamics
	mat<X_DIM,X_DIM> A = mat<X_DIM,X_DIM>::Zero();
	mat<X_DIM,Q_DIM> M = mat<X_DIM,Q_DIM>::Zero();
	linearize_dynfunc(x_t_t, u_t, vec<Q_DIM>::Zero(), A, M);

	mat<X_DIM,X_DIM> sigma_tp1_t = A*sigma_t_t*A.transpose() + M*Q*M.transpose();
	vec<X_DIM> x_tp1_t = dynfunc(x_t_t, u_t, vec<Q_DIM>::Zero());

	// propagate belief through observation
	mat<Z_DIM,X_DIM> H = mat<Z_DIM,X_DIM>::Zero();
	mat<Z_DIM,R_DIM> N = mat<Z_DIM,R_DIM>::Zero();
	linearize_obsfunc(x_tp1_t, vec<R_DIM>::Zero(), H, N);

	mat<Z_DIM,Z_DIM> delta = delta_matrix(x_tp1_t, object, INFINITY);
	// calculate Kalman gain
	mat<X_DIM,Z_DIM> K = sigma_tp1_t*H.transpose()*delta*(delta*H*sigma_tp1_t*H.transpose()*delta + R).inverse()*delta;

	// update based on noisy measurement
	x_tp1_tp1 = x_tp1_t + K*(z_tp1_real - obsfunc(x_tp1_t, x_tp1_t.segment<C_DIM>(J_DIM), vec<R_DIM>::Zero()));
	sigma_tp1_tp1 = (mat<X_DIM,X_DIM>::Identity() - K*H)*sigma_tp1_t;
}

std::vector<Beam> PlanarSystem::get_fov(const vec<X_DIM>& x) {
	std::vector<Beam> beams, new_beams;

	// start out with full field of view
	double angle = x(J_DIM-1);
	vec<C_DIM> dir_right, dir_left;
	dir_right << sin(angle+camera_fov/2.0),
			cos(angle+camera_fov/2.0);
	dir_left << sin(angle-camera_fov/2.0),
			cos(angle-camera_fov/2.0);
	beams.push_back(Beam(camera_origin,
			camera_origin + camera_max_dist*dir_right,
			camera_origin + camera_max_dist*dir_left));

	std::vector<Segment> link_segments = get_link_segments(x);

	for(int l=0; l < link_segments.size(); ++l) {
		new_beams.clear();
		for(int i=0; i < beams.size(); ++i) {
			std::vector<Beam> new_beams_i = beams[i].truncate(link_segments[l]);
			for(int j=0; j < new_beams_i.size(); ++j) {
				// prevent small beams from being added
				if (new_beams_i[j].top_length() > epsilon) {
					new_beams.push_back(new_beams_i[j]);
				}
			}
		}
		beams = new_beams;
	}

	return beams;
}

vec<C_DIM> PlanarSystem::get_ee_pos(const vec<E_DIM>& j) {
	vec<X_DIM> x = vec<X_DIM>::Zero();
	x.segment<E_DIM>(0) = j;
	return get_link_segments(x).back().p1;
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

		if (iter > max_iter) {
			return false;
		}
		iter++;
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


void PlanarSystem::display(vec<X_DIM>& x, mat<X_DIM,X_DIM>& sigma, bool pause) {
	std::vector<vec<X_DIM>, aligned_allocator<vec<X_DIM>>> X(1, x);
	std::vector<mat<X_DIM,X_DIM>, aligned_allocator<mat<X_DIM,X_DIM>>> S(1, sigma);
	display(X, S);
}

void PlanarSystem::display(std::vector<vec<X_DIM>, aligned_allocator<vec<X_DIM>>>& X, mat<X_DIM,X_DIM>& sigma0, std::vector<vec<U_DIM>, aligned_allocator<vec<U_DIM>>>& U, const double alpha, bool pause) {
	int T = X.size();

	std::vector<mat<X_DIM,X_DIM>, aligned_allocator<mat<X_DIM,X_DIM>>> S(T);
	S[0] = sigma0;
	vec<X_DIM> x_tp1;
	for(int t=0; t < T-1; ++t) {
		belief_dynamics(X[t], S[t], U[t], alpha, x_tp1, S[t+1]);
	}
	display(X, S, pause);
}

void PlanarSystem::display(std::vector<vec<X_DIM>, aligned_allocator<vec<X_DIM>>>& X, std::vector<mat<X_DIM,X_DIM>, aligned_allocator<mat<X_DIM,X_DIM>>>& S, bool pause) {
	Py_Initialize();
	np::initialize();

	py::numeric::array::set_module_and_type("numpy", "ndarray");

	py::list X_pylist, S_pylist;
	for(int t=0; t < X.size(); ++t) {
		X_pylist.append(planar_utils::eigen_to_ndarray(X[t]));
		S_pylist.append(planar_utils::eigen_to_ndarray(S[t]));
	}

	np::ndarray robot_origin_ndarray = planar_utils::eigen_to_ndarray(robot_origin);
	np::ndarray link_lengths_ndarray = planar_utils::eigen_to_ndarray(link_lengths);
	np::ndarray camera_origin_ndarray = planar_utils::eigen_to_ndarray(camera_origin);
	np::ndarray object_ndarray = planar_utils::eigen_to_ndarray(object);

	py::list beams_pylist;
	std::vector<Beam> beams = get_fov(X.back());

	for(int i=0; i < beams.size(); ++i) {
		mat<2,3> m;
		m << beams[i].base, beams[i].a, beams[i].b;
		beams_pylist.append(planar_utils::eigen_to_ndarray(m));
	}

	std::string working_dir = boost::filesystem::current_path().normalize().string();
	std::string bsp_dir = working_dir.substr(0,working_dir.find("bsp"));
	std::string planar_dir = bsp_dir + "bsp/planar";

	try
	{
		py::object main_module = py::import("__main__");
		py::object main_namespace = main_module.attr("__dict__");
		py::exec("import sys, os", main_namespace);
		py::exec(py::str("sys.path.append('"+planar_dir+"')"), main_namespace);
		py::object plot_mod = py::import("plot_planar");
		py::object plot_planar = plot_mod.attr("plot_planar");

		plot_planar(X_pylist, S_pylist, robot_origin_ndarray, link_lengths_ndarray,
				camera_origin_ndarray, camera_fov, object_ndarray, beams_pylist);

		if (pause) {
			LOG_INFO("Press enter to continue");
			py::exec("raw_input()",main_namespace);
		}
	}
	catch(py::error_already_set const &)
	{
		PyErr_Print();
	}
}

void PlanarSystem::get_limits(vec<X_DIM>& x_min, vec<X_DIM>& x_max, vec<U_DIM>& u_min, vec<U_DIM>& u_max) {
	x_min = this->x_min;
	x_max = this->x_max;
	u_min = this->u_min;
	u_max = this->u_max;
}

double PlanarSystem::cost(const std::vector<vec<X_DIM>, aligned_allocator<vec<X_DIM>>>& X, const mat<X_DIM,X_DIM>& sigma0,
		const std::vector<vec<U_DIM>, aligned_allocator<vec<U_DIM>>>& U, const double alpha) {
	double cost = 0;
	int T = X.size();

	vec<X_DIM> x_tp1 = vec<X_DIM>::Zero();
	mat<X_DIM,X_DIM> sigma_t = sigma0, sigma_tp1 = mat<X_DIM,X_DIM>::Zero();
	for(int t=0; t < T-1; ++t) {
		belief_dynamics(X[t], sigma_t, U[t], alpha, x_tp1, sigma_tp1);
		cost += alpha_control*U[t].dot(U[t]);
		if (t < T-2) {
			cost += alpha_belief*sigma_tp1.trace();
		} else {
			cost += alpha_final_belief*sigma_tp1.trace();
		}
		sigma_t = sigma_tp1;
	}

	vec<C_DIM> final_ee_position = get_link_segments(x_tp1).back().p1;
	vec<C_DIM> goal_error = final_ee_position - x_tp1.segment<C_DIM>(J_DIM);
	cost += alpha_goal*goal_error.dot(goal_error);

	return cost;
}

vec<TOTAL_VARS> PlanarSystem::cost_grad(std::vector<vec<X_DIM>, aligned_allocator<vec<X_DIM>>>& X, const mat<X_DIM,X_DIM>& sigma0,
		std::vector<vec<U_DIM>, aligned_allocator<vec<U_DIM>>>& U, const double alpha) {
	int T = X.size();

	vec<TOTAL_VARS> grad;

	double orig, cost_p, cost_m;
	int index = 0;
	for(int t=0; t < T; ++t) {
		for(int i=0; i < X_DIM; ++i) {
			orig = X[t][i];

			X[t][i] = orig + step;
			cost_p = cost(X, sigma0, U, alpha);

			X[t][i] = orig - step;
			cost_m = cost(X, sigma0, U, alpha);

			grad(index++) = (cost_p - cost_m) / (2*step);
		}

		if (t < T-1) {
			for(int i=0; i < U_DIM; ++i) {
				orig = U[t][i];

				U[t][i] = orig + step;
				cost_p = cost(X, sigma0, U, alpha);

				U[t][i] = orig - step;
				cost_m = cost(X, sigma0, U, alpha);

				grad(index++) = (cost_p - cost_m) / (2*step);
			}
		}
	}
	return grad;
}

void PlanarSystem::cost_and_cost_grad(std::vector<vec<X_DIM>, aligned_allocator<vec<X_DIM>>>& X, const mat<X_DIM,X_DIM>& sigma0,
		std::vector<vec<U_DIM>, aligned_allocator<vec<U_DIM>>>& U, const double alpha, const bool use_fadbad,
		double& cost, vec<TOTAL_VARS>& grad) {
	if (!use_fadbad) {
		cost = this->cost(X, sigma0, U, alpha);
		grad = this->cost_grad(X, sigma0, U, alpha);
	} else {
		std::vector<vecb<X_DIM>, aligned_allocator<vecb<X_DIM>>> X_b(X.size());
		matb<X_DIM,X_DIM> sigma0_b = sigma0.cast<bdouble>();
		std::vector<vecb<U_DIM>, aligned_allocator<vecb<U_DIM>>> U_b(U.size());
		bdouble alpha_b = (bdouble)alpha;

		for(int t=0; t < X.size(); ++t) {
			X_b[t] = X[t].cast<bdouble>();
			if(t < X.size()-1) {
				U_b[t] = U[t].cast<bdouble>();
			}
		}

		bdouble cost_b = fps.cost(X_b, sigma0_b, U_b, alpha_b);
		cost = cost_b.x();
		cost_b.diff(0,1);
		int index = 0;
		for(int t=0; t < X.size(); ++t) {
			for(int i=0; i < X_DIM; ++i) {
				grad(index++) = X_b[t](i).d(0);
			}

			if (t < X.size()-1) {
				for(int i=0; i < U_DIM; ++i) {
					grad(index++) = U_b[t](i).d(0);
				}
			}
		}
	}
}

/**
 * Private methods
 */

std::vector<Segment> PlanarSystem::get_link_segments(const vec<X_DIM>& x) {
	std::vector<Segment> link_segments;

	vec<C_DIM> p0 = robot_origin, p1;
	double joint0 = 0, joint1;

	for(int i=0; i < link_lengths.rows(); ++i) {
		joint1 = joint0 + x[i];
		p1 << p0(0) + sin(joint1)*link_lengths(i),
				p0(1) + cos(joint1)*link_lengths(i);
		link_segments.push_back(Segment(p0, p1));

		joint0 = joint1;
		p0 = p1;
	}

	return link_segments;
}

void PlanarSystem::linearize_dynfunc(const vec<X_DIM>& x, const vec<U_DIM>& u, const vec<Q_DIM>& q, mat<X_DIM,X_DIM>& A, mat<X_DIM,Q_DIM>& M) {
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
