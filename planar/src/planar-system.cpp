#include "../include/planar-system.h"

/**
 * Constructors and initializers
 */

PlanarSystem::PlanarSystem(const vec& camera_origin, const vec& object, bool is_static) {
	init(camera_origin, object, is_static);
}

void PlanarSystem::init(const vec& camera_origin, const vec& object, bool is_static) {
	this->camera_origin = camera_origin;
	camera_fov = M_PI/4;
	camera_max_dist = 10;
	this->object = object;
	this->is_static = is_static;

	robot_origin = zeros<vec>(2);
	link_lengths = {5, 5, 4};

	J_DIM = 4;
	C_DIM = 2;

	X_DIM = 6;
	U_DIM = 4;
	Z_DIM = 6;

	mat Q = eye<mat>(U_DIM, U_DIM);
	mat R = eye<mat>(Z_DIM, Z_DIM);

	x_min = {-M_PI/2, -M_PI/2, -M_PI/2, -M_PI/2, -INFTY, -INFTY};
	x_max = {M_PI/2, M_PI/2, M_PI/2, M_PI/2, INFTY, INFTY};

	u_min = {-M_PI/2, -M_PI/2, -M_PI/2, (is_static) ? 0 : -M_PI/2};
	u_max = {M_PI/2, M_PI/2, M_PI/2, (is_static) ? 0 : M_PI/2};

	Q_DIM = Q.n_rows;
	R_DIM = R.n_rows;

	DT = 1.0;
}

/**
 * Public methods
 */

vec PlanarSystem::dynfunc(const vec& x, const vec& u, const vec& q) {
	// TODO: should enforce joint limits?
	vec x_new = x;
	x_new.subvec(span(0, J_DIM-1)) += DT*(u + q);
	return x_new;
}

// For planning, object comes from state
// For execution update, object is this->object
vec PlanarSystem::obsfunc(const vec& x, const vec& object, const vec& r) {
	vec z(Z_DIM);
	z(0) = x(0); // joint 0
	z(1) = x(1); // joint 1
	z(2) = x(2); // joint 2
	z(3) = x(3); // camera angle
	z(4) = object(0) - camera_origin(0); // delta x to object
	z(5) = object(1) - camera_origin(1); // delta y to object
	return z;
}

/**
 * \brief Propagates belief through EKF with max-likelihood noise assumption
 */
void PlanarSystem::belief_dynamics(const vec& x_t, const mat& sigma_t, const vec& u_t, const double alpha, vec& x_tp1, mat& sigma_tp1) {
	// propagate dynamics
	x_tp1 = dynfunc(x_t, u_t, zeros<vec>(Q_DIM));

	// propagate belief through dynamics
	mat A(X_DIM, X_DIM), M(X_DIM, Q_DIM);
	linearize_dynfunc(x_t, u_t, zeros<vec>(Q_DIM), A, M);

	mat sigma_tp1_bar = A*sigma_t*A.t() + M*Q*M.t();

	// propagate belief through observation
	mat H(Z_DIM, X_DIM), N(Z_DIM, R_DIM);
	linearize_obsfunc(x_tp1, zeros<vec>(R_DIM), H, N);

	mat delta = delta_matrix(x_tp1, alpha);
	mat K = sigma_tp1_bar*H.t()*delta*inv(delta*H*sigma_tp1_bar*H.t()*delta + R)*delta;
	sigma_tp1 = (eye<mat>(X_DIM,X_DIM) - K*H)*sigma_tp1_bar;
}

void PlanarSystem::execute_control_step(const vec& x_t_real, const vec& x_t, const mat& sigma_t, const vec& u_t,
		vec& x_tp1_real, vec& x_tp1, mat& sigma_tp1) {

}

/**
 * \brief Creates full camera FOV, then adds occlusions
 * (in our case only the robot links) to truncate the FOV
 */
std::vector<Beam> PlanarSystem::get_fov(const vec& x) {
	std::vector<Beam> beams, new_beams;

	// start out with full field of view
	double angle = x(J_DIM-1);
	vec dir_right = {sin(angle+camera_fov/2.0), cos(angle+camera_fov/2.0)};
	vec dir_left = {sin(angle-camera_fov/2.0), cos(angle-camera_fov/2.0)};
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

void PlanarSystem::display(std::vector<vec>& X, std::vector<mat>& S, bool pause) {
	Py_Initialize();
	np::initialize();

	py::numeric::array::set_module_and_type("numpy", "ndarray");

	py::list X_pylist, S_pylist;
	for(int t=0; t < X.size(); ++t) {
		X_pylist.append(planar_utils::arma_to_ndarray(X[t]));
		S_pylist.append(planar_utils::arma_to_ndarray(S[t]));
	}

	np::ndarray robot_origin_ndarray = planar_utils::arma_to_ndarray(robot_origin);
	np::ndarray link_lengths_ndarray = planar_utils::arma_to_ndarray(link_lengths);
	np::ndarray camera_origin_ndarray = planar_utils::arma_to_ndarray(camera_origin);
	np::ndarray object_ndarray = planar_utils::arma_to_ndarray(object);

	py::list beams_pylist;
	std::vector<Beam> beams = get_fov(X.back());

	for(int i=0; i < beams.size(); ++i) {
		mat m = join_horiz(beams[i].base, join_horiz(beams[i].a, beams[i].b));
		beams_pylist.append(planar_utils::arma_to_ndarray(m));
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

void PlanarSystem::get_limits(vec& x_min, vec& x_max, vec& u_min, vec& u_max) {
	x_min = this->x_min;
	x_max = this->x_max;
	u_min = this->u_min;
	u_max = this->u_max;
}

double PlanarSystem::cost(const std::vector<vec>& X, const mat& sigma0, const std::vector<vec>& U, const double alpha) {
	double cost = 0;
	int T = X.size();
	int X_DIM = X[0].n_rows;

	vec x_tp1(X_DIM);
	mat sigma_t = sigma0, sigma_tp1(X_DIM, X_DIM);
	for(int t=0; t < T-1; ++t) {
		belief_dynamics(X[t], sigma_t, U[t], alpha, x_tp1, sigma_tp1);
		cost += alpha_control*dot(U[t], U[t]);
		if (t < T-2) {
			cost += alpha_belief*trace(sigma_tp1);
		} else {
			cost += alpha_final_belief*trace(sigma_tp1);
		}
		sigma_t = sigma_tp1;
	}

	vec final_ee_position = get_link_segments(x_tp1).back().p1;
	vec goal_error = final_ee_position - x_tp1.subvec(J_DIM, X_DIM-1);
	cost += alpha_goal*dot(goal_error, goal_error);

	return cost;
}

vec PlanarSystem::cost_grad(std::vector<vec>& X, const mat& sigma0, std::vector<vec>& U, const double alpha) {
	int T = X.size();
	int X_DIM = X[0].n_rows;
	int U_DIM = U[0].n_rows;

	vec grad(T*X_DIM + (T-1)*U_DIM);
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

std::vector<Segment> PlanarSystem::get_link_segments(const vec& x) {
	std::vector<Segment> link_segments;

	vec p0 = robot_origin, p1;
	double joint0 = 0, joint1;

	for(int i=0; i < link_lengths.n_rows; ++i) {
		joint1 = joint0 + x[i];
		p1 << p0(0) + sin(joint1)*link_lengths(i) <<
				p0(1) + cos(joint1)*link_lengths(i);
		link_segments.push_back(Segment(p0, p1));

		joint0 = joint1;
		p0 = p1;
	}

	return link_segments;
}

/**
 * \param x is X_DIM by 1
 * \param u is U_DIM by 1
 * \param q is Q_DIM by 1
 * \param A is X_DIM by X_DIM
 * \param M is X_DIM by Q_DIM
 */
void PlanarSystem::linearize_dynfunc(const vec& x, const vec& u, const vec& q, mat& A, mat& M) {
	vec q_zero = zeros<vec>(Q_DIM);
	vec x_p = x, x_m = x;
	for(int i=0; i < X_DIM; ++i) {
		x_p(i) += step; x_m(i) -= step;
		A.submat(span(0, X_DIM-1), span(i, i)) = (dynfunc(x_p, u, q_zero) - dynfunc(x_m, u, q_zero)) / (2*step);
		x_p(i) = x(i); x_m(i) = x(i);
	}

	vec q_p = q, q_m = q;
	for(int i=0; i < Q_DIM; ++i) {
		q_p(i) += step; q_m(i) -= step;
		M.submat(span(0, Q_DIM-1), span(i, i)) = (dynfunc(x, u, q_p) - dynfunc(x, u, q_m)) / (2*step);
	}
}

/**
 * \param x is X_DIM by 1
 * \param r is R_DIM by 1
 * \param H is Z_DIM by X_DIM
 * \param N is Z_DIM by R_DIM
 */
void PlanarSystem::linearize_obsfunc(const vec& x, const vec& r, mat& H, mat& N) {
	vec x_p = x, x_m = x;
	for(int i=0; i < X_DIM; ++i) {
		x_p(i) += step; x_m(i) -= step;
		H.submat(span(0, Z_DIM-1), span(i, i)) = (obsfunc(x_p.subvec(0, J_DIM-1), x_p.subvec(J_DIM, X_DIM-1), r) -
				obsfunc(x_p.subvec(0, J_DIM-1), x_p.subvec(J_DIM, X_DIM-1), r)) / (2*step);
		x_p(i) = x(i); x_m(i) = x(i);
	}

	vec r_p = r, r_m = r;
	for(int i=0; i < R_DIM; ++i) {
		r_p(i) += step; r_m(i) -= step;
		N.submat(span(0, Z_DIM-1), span(i, i)) = (obsfunc(x.subvec(0, J_DIM-1), x.subvec(J_DIM, X_DIM-1), r_p) -
				obsfunc(x.subvec(0, J_DIM-1), x.subvec(J_DIM, X_DIM-1), r_m)) / (2*step);
		r_p(i) = x(i); r_m(i) = x(i);
	}
}

mat PlanarSystem::delta_matrix(const vec& x, const double alpha) {
	mat delta(Z_DIM, Z_DIM, fill::zeros);

	for(int i=0; i < J_DIM; ++i) {
		delta(i, i) = 1; // TODO: should this depend on SD of link segments?
	}

	std::vector<Beam> fov = get_fov(x);
	vec object = x.subvec(J_DIM, X_DIM-1);
	double sd = geometry2d::signed_distance(object, fov);
	double sd_sigmoid = 1 - 1/(1 + exp(-alpha*sd));
	for(int i=J_DIM; i < X_DIM; ++i) {
		delta(i, i) = sd_sigmoid;
	}

	return delta;
}

