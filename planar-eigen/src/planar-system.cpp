#include "../include/planar-system.h"

/**
 * Constructors and initializers
 */

PlanarSystem::PlanarSystem(const VectorXd& camera_origin, const VectorXd& object, bool is_static) {
	init(camera_origin, object, is_static);
}

void PlanarSystem::init(const VectorXd& camera_origin, const VectorXd& object, bool is_static) {
	this->camera_origin = camera_origin;
	camera_fov = M_PI/4;
	camera_max_dist = 15;
	this->object = object;
	this->is_static = is_static;

	robot_origin = Vector2d::Zero();
	link_lengths = VectorXd(3);
	link_lengths << 5, 5, 4;

	Q = (M_PI/4)*MatrixXd::Identity(U_DIM, U_DIM);
//	R = 10*eye<mat>(Z_DIM, Z_DIM);
	VectorXd R_diag(R_DIM);
	R_diag << (M_PI/4)*VectorXd::Ones(J_DIM),
			5*VectorXd::Ones(C_DIM);
	R = R_diag.asDiagonal();

	// x bound for object is sum of link_lengths
	double max_link_length = link_lengths.sum();
	x_min = VectorXd::Zero(X_DIM);
	x_max = VectorXd::Zero(X_DIM);
	x_min << -M_PI/2, -M_PI/2, -M_PI/2, -M_PI/2, -max_link_length, -max_link_length;
	x_max << M_PI/2, M_PI/2, M_PI/2, M_PI/2, max_link_length, max_link_length;

	double max_input = M_PI/12;
	u_min = VectorXd::Zero(U_DIM);
	u_max = VectorXd::Zero(U_DIM);
	u_min << -max_input, -max_input, -max_input, (is_static) ? 0 : -max_input;
	u_max << max_input, max_input, max_input, (is_static) ? 0 : max_input;
}


/**
 * Public methods
 */

VectorXd PlanarSystem::dynfunc(const VectorXd& x, const VectorXd& u, const VectorXd& q, bool enforce_limits) {
	VectorXd x_new = x;
	x_new.segment<J_DIM>(0) += DT*(u + q);
	return x_new;
}

VectorXd PlanarSystem::obsfunc(const VectorXd& x, const VectorXd& object, const VectorXd& r) {
	VectorXd z(Z_DIM);
	z(0) = x(0); // joint 0
	z(1) = x(1); // joint 1
	z(2) = x(2); // joint 2
	z(3) = x(3); // camera angle
	z(4) = object(0) - camera_origin(0); // delta x to object
	z(5) = object(1) - camera_origin(1); // delta y to object
	return z + r;
}


MatrixXd PlanarSystem::delta_matrix(const VectorXd& x, const VectorXd& object, const double alpha) {
	MatrixXd delta = MatrixXd::Zero(Z_DIM, Z_DIM);

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
void PlanarSystem::belief_dynamics(const VectorXd& x_t, const MatrixXd& sigma_t, const VectorXd& u_t, const double alpha, VectorXd& x_tp1, MatrixXd& sigma_tp1) {
	// propagate dynamics
	x_tp1 = dynfunc(x_t, u_t, VectorXd::Zero(Q_DIM));

	// propagate belief through dynamics
	MatrixXd A = MatrixXd::Zero(X_DIM, X_DIM);
	MatrixXd M = MatrixXd::Zero(X_DIM, Q_DIM);
	linearize_dynfunc(x_t, u_t, VectorXd::Zero(Q_DIM), A, M);

	MatrixXd sigma_tp1_bar = A*sigma_t*A.transpose() + M*Q*M.transpose();

	// propagate belief through observation
	MatrixXd H = MatrixXd::Zero(Z_DIM, X_DIM);
	MatrixXd N = MatrixXd::Zero(Z_DIM, R_DIM);
	linearize_obsfunc(x_tp1, VectorXd::Zero(R_DIM), H, N);

	MatrixXd delta = delta_matrix(x_tp1, x_tp1.segment<C_DIM>(J_DIM), alpha);
	MatrixXd K = sigma_tp1_bar*H.transpose()*delta*(delta*H*sigma_tp1_bar*H.transpose()*delta + R).inverse()*delta;
	sigma_tp1 = (MatrixXd::Identity(X_DIM,X_DIM) - K*H)*sigma_tp1_bar;
}

void PlanarSystem::execute_control_step(const VectorXd& x_t_real, const VectorXd& x_t_t, const MatrixXd& sigma_t_t, const VectorXd& u_t,
		VectorXd& x_tp1_real, VectorXd& x_tp1_tp1, MatrixXd& sigma_tp1_tp1) {
	// find next real state from input + noise
	VectorXd control_noise = VectorXd::Zero(Q_DIM);// + chol(.2*Q)*randn<vec>(Q_DIM);
	VectorXd obs_noise = VectorXd::Zero(R_DIM);// + chol(.2*R)*randn<vec>(R_DIM);
	x_tp1_real = dynfunc(x_t_real, u_t, control_noise, true);
	VectorXd z_tp1_real = obsfunc(x_tp1_real, this->object, obs_noise);

	// now do update based on current belief (x_t, sigma_t)

	// propagate belief through dynamics
	MatrixXd A = MatrixXd::Zero(X_DIM, X_DIM);
	MatrixXd M = MatrixXd::Zero(X_DIM, Q_DIM);
	linearize_dynfunc(x_t_t, u_t, VectorXd::Zero(Q_DIM), A, M);

	MatrixXd sigma_tp1_t = A*sigma_t_t*A.transpose() + M*Q*M.transpose();
	VectorXd x_tp1_t = dynfunc(x_t_t, u_t, VectorXd::Zero(Q_DIM));

	// propagate belief through observation
	MatrixXd H = MatrixXd::Zero(Z_DIM, X_DIM);
	MatrixXd N = MatrixXd::Zero(Z_DIM, R_DIM);
	linearize_obsfunc(x_tp1_t, VectorXd::Zero(R_DIM), H, N);

	MatrixXd delta = delta_matrix(x_tp1_t, object, INFINITY);
	// calculate Kalman gain
	MatrixXd K = sigma_tp1_t*H.transpose()*delta*(delta*H*sigma_tp1_t*H.transpose()*delta + R).inverse()*delta;

	// update based on noisy measurement
	x_tp1_tp1 = x_tp1_t + K*(z_tp1_real - obsfunc(x_tp1_t, x_tp1_t.segment<C_DIM>(J_DIM), VectorXd::Zero(R_DIM)));
	sigma_tp1_tp1 = (MatrixXd::Identity(X_DIM,X_DIM) - K*H)*sigma_tp1_t;
}

std::vector<Beam> PlanarSystem::get_fov(const VectorXd& x) {
	std::vector<Beam> beams, new_beams;

	// start out with full field of view
	double angle = x(J_DIM-1);
	Vector2d dir_right, dir_left;
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

std::vector<Segment> PlanarSystem::get_link_segments(const VectorXd& x) {
	std::vector<Segment> link_segments;

	Vector2d p0 = robot_origin, p1;
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

void PlanarSystem::display(VectorXd& x, MatrixXd& sigma, bool pause) {
	std::vector<VectorXd> X(1, x);
	std::vector<MatrixXd> S(1, sigma);
	display(X, S);
}

void PlanarSystem::display(std::vector<VectorXd>& X, MatrixXd& sigma0, std::vector<VectorXd>& U, const double alpha, bool pause) {
	int T = X.size();

	std::vector<MatrixXd> S(T);
	S[0] = sigma0;
	VectorXd x_tp1(X_DIM);
	for(int t=0; t < T-1; ++t) {
		belief_dynamics(X[t], S[t], U[t], alpha, x_tp1, S[t+1]);
	}
	display(X, S, pause);
}

void PlanarSystem::display(std::vector<VectorXd>& X, std::vector<MatrixXd>& S, bool pause) {
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
		MatrixXd m(2, 3);
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

void PlanarSystem::get_limits(VectorXd& x_min, VectorXd& x_max, VectorXd& u_min, VectorXd& u_max) {
	x_min = this->x_min;
	x_max = this->x_max;
	u_min = this->u_min;
	u_max = this->u_max;
}

double PlanarSystem::cost(const std::vector<VectorXd>& X, const MatrixXd& sigma0, const std::vector<VectorXd>& U, const double alpha) {
	double cost = 0;
	int T = X.size();

	VectorXd x_tp1(X_DIM);
	MatrixXd sigma_t = sigma0, sigma_tp1(X_DIM, X_DIM);
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

	VectorXd final_ee_position = get_link_segments(x_tp1).back().p1;
	VectorXd goal_error = final_ee_position - x_tp1.segment<C_DIM>(J_DIM);
	cost += alpha_goal*goal_error.dot(goal_error);

	return cost;
}

VectorXd PlanarSystem::cost_grad(std::vector<VectorXd>& X, const MatrixXd& sigma0, std::vector<VectorXd>& U, const double alpha) {
	int T = X.size();

	VectorXd grad(T*X_DIM + (T-1)*U_DIM);
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

/**
 * Private methods
 */

/**
 * \param x is X_DIM by 1
 * \param u is U_DIM by 1
 * \param q is Q_DIM by 1
 * \param A is X_DIM by X_DIM
 * \param M is X_DIM by Q_DIM
 */
void PlanarSystem::linearize_dynfunc(const VectorXd& x, const VectorXd& u, const VectorXd& q, MatrixXd& A, MatrixXd& M) {
	A = MatrixXd::Identity(X_DIM, X_DIM);

	for(int j=0; j < M.cols(); ++j) {
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

/**
 * \param x is X_DIM by 1
 * \param r is R_DIM by 1
 * \param H is Z_DIM by X_DIM
 * \param N is Z_DIM by R_DIM
 */
void PlanarSystem::linearize_obsfunc(const VectorXd& x, const VectorXd& r, MatrixXd& H, MatrixXd& N) {
	H = MatrixXd::Identity(Z_DIM, X_DIM);
	N = MatrixXd::Identity(Z_DIM, R_DIM);

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
