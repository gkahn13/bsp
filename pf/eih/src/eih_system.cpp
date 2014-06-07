#include "../include/eih_system.h"

/**
 * EihSystem constructors and initializers
 */

EihSystem::EihSystem(rave::EnvironmentBasePtr e, Manipulator *m, KinectSensor *k) : env(e), manip(m), kinect(k) {
	mat uMin, uMax;
	manip->get_limits(uMin, uMax);
	init(uMin, uMax);
}

EihSystem::EihSystem(rave::EnvironmentBasePtr e, Manipulator *m, KinectSensor *k,
		ObsType obs_type) : env(e), manip(m), kinect(k) {
	mat uMin, uMax;
	manip->get_limits(uMin, uMax);
	init(uMin, uMax, obs_type);
}

EihSystem::EihSystem(rave::EnvironmentBasePtr e, Manipulator *m, KinectSensor *k,
		ObsType obs_type, int T) : env(e), manip(m), kinect(k) {
	mat uMin, uMax;
	manip->get_limits(uMin, uMax);
	init(uMin, uMax, obs_type, T);
}

EihSystem::EihSystem(rave::EnvironmentBasePtr e, Manipulator *m, KinectSensor *k,
		ObsType obs_type, int T, const mat &uMin, const mat &uMax) : env(e), manip(m), kinect(k) {
	init(uMin, uMax, obs_type, T);
}

void EihSystem::init(const mat &uMin, const mat &uMax, ObsType obs_type, int T, double DT) {
	mat j = manip->get_joint_values();
	X_DIM = MAX(j.n_rows, j.n_cols);
	U_DIM = X_DIM;

	this->obs_type = obs_type;
	mat R_diag;
//	double height = kinect->get_height(), width = kinect->get_width();
//	double min_range = kinect->get_min_range(), max_range = kinect->get_max_range(),
//			optimal_range = kinect->get_optimal_range();
	switch(obs_type) {
	case ObsType::fov:
//		R_diag << 100*height << 100*width << 1*(max_range - min_range);
//		desired_observations = std::vector<mat>(1);
//		desired_observations[0] << height/2 << width/2 << optimal_range;
		break;
	case ObsType::fov_occluded:
//		R_diag << .5 << .01;
//		desired_observations = std::vector<mat>(2);
//		desired_observations[0] << 0 << UNKNOWN;
//		desired_observations[1] << 1 << 1;
		break;
	case ObsType::fov_occluded_color:
	default:
		R_diag << .5 << .01 << .05 << 1 << 1;
		desired_observations = std::vector<mat>(3);
		desired_observations[0] << 0 << UNKNOWN << UNKNOWN << UNKNOWN << UNKNOWN; // outside FOV
		desired_observations[1] << 1 << 1 << UNKNOWN << UNKNOWN << UNKNOWN; // inside FOV, occluded
		desired_observations[2] << 1 << .5 << 0 << 1 << 1; // inside FOV, not occluded, red
		break;
	}
	R = diagmat(R_diag);
	Z_DIM = R.n_rows;

	this->DT = DT;

	manip->get_limits(xMin, xMax);
	this->uMin = uMin;
	this->uMax = uMax;

	kinect->power_on();
	boost::this_thread::sleep(boost::posix_time::seconds(2));

	use_casadi = false;
}

/**
 * EihSystem public methods
 */

mat EihSystem::dynfunc(const mat &x, const mat &u) {
	int X_DIM = std::max(x.n_rows, x.n_cols);
	mat x_new = x + DT*u;

	for(int i=0; i < X_DIM; ++i) {
		x_new(i) = MAX(x_new(i), xMin(i));
		x_new(i) = MIN(x_new(i), xMax(i));
	}

	return x_new;
}

mat EihSystem::obsfunc(const mat& x, const mat& t, const mat& r) {

}

double EihSystem::obsfunc_continuous_weight(const mat &particle, const cube &image, const mat &z_buffer, bool plot) {
	switch(obs_type) {
	case ObsType::fov:
		return obsfunc_continuous_weight_fov(particle, image, z_buffer, plot);
	case ObsType::fov_occluded:
		return obsfunc_continuous_weight_fov_occluded(particle, image, z_buffer, plot);
	case ObsType::fov_occluded_color:
		return obsfunc_continuous_weight_fov_occluded_color(particle, image, z_buffer, plot);
	default:
		return obsfunc_continuous_weight_fov_occluded_color(particle, image, z_buffer, plot);
	}
}

double EihSystem::obsfunc_discrete_weight(const mat &particle, const cube &image, const mat &z_buffer, bool plot) {
	switch(obs_type) {
	case ObsType::fov:
		return obsfunc_discrete_weight_fov(particle, image, z_buffer, plot);
	case ObsType::fov_occluded:
		return obsfunc_discrete_weight_fov_occluded(particle, image, z_buffer, plot);
	case ObsType::fov_occluded_color:
		return obsfunc_discrete_weight_fov_occluded_color(particle, image, z_buffer, plot);
	default:
		return obsfunc_discrete_weight_fov_occluded_color(particle, image, z_buffer, plot);
	}
}

double EihSystem::obsfunc_continuous_weight_fov(const mat &particle, const cube &image, const mat &z_buffer, bool plot) {
	// compute comparison
	double height = kinect->get_height(), width = kinect->get_width();
	double min_range = kinect->get_min_range(), max_range = kinect->get_max_range(),
			optimal_range = kinect->get_optimal_range();

	mat R_diag;
	R_diag << 100*height << 100*width << 1*(max_range - min_range);
	mat R = diagmat(R_diag);

	mat z_low_weight;
	z_low_weight << height/2 << width/2 << optimal_range;

	// actual observation
	mat z;

	mat pixel = kinect->get_pixel_from_point(particle, false);
	double y = pixel(0), x = pixel(1);

	double d = norm(rave_utils::rave_vec_to_mat(kinect->get_pose().trans) - particle, 2);

	z << y << x << d;

	double max_likelihood = gauss_likelihood(zeros<mat>(Z_DIM,1), R);

	mat e = z - z_low_weight;
	double weight = (max_likelihood - gauss_likelihood(e.t(), R)) / max_likelihood;

	if (plot) {
		rave::Vector color = (weight > .5) ? rave::Vector(1, 1, 1) : rave::Vector(0, 1, 0);
		rave::Vector particle_vec = rave_utils::mat_to_rave_vec(particle);
		handles.push_back(rave_utils::plot_point(env, particle_vec, color));
	}

	return weight;
}

double EihSystem::obsfunc_discrete_weight_fov(const mat &particle, const cube &image, const mat &z_buffer, bool plot) {
	double weight;
	rave::Vector color;

	if (!kinect->is_visible(particle)) {
		weight = 1;
		color = rave::Vector(0, 1, 0);
	} else {
		weight = 0;
		color = rave::Vector(1, 1, 1);
	}

	if (plot) {
		rave::Vector particle_vec = rave_utils::mat_to_rave_vec(particle);
		handles.push_back(rave_utils::plot_point(env, particle_vec, color));
	}

	return weight;
}

double EihSystem::obsfunc_continuous_weight_fov_occluded(const mat &particle, const cube &image, const mat &z_buffer, bool plot) {
	double weight;
	rave::Vector color;

	double fov_weight = obsfunc_continuous_weight_fov(particle, image, z_buffer, false);

	if (kinect->is_visible(particle)) {
		mat kinect_pos = rave_utils::rave_vec_to_mat(kinect->get_pose().trans);

		double particle_dist = norm(particle - kinect_pos, 2);
		mat pixel = kinect->get_pixel_from_point(particle);
		int y = int(pixel(0)), x = int(pixel(1));
		double sd = particle_dist - z_buffer(y, x);

		double occluded_weight = sigmoid(sd, 5);

		weight = std::max(fov_weight, occluded_weight);

		color = (occluded_weight < .5) ? rave::Vector(1, 1, 1) : rave::Vector(0, 0, 0);

	} else {
		weight = fov_weight;
		color = rave::Vector(0, 1, 0);
	}

	if (plot) {
		rave::Vector particle_vec = rave_utils::mat_to_rave_vec(particle);
		handles.push_back(rave_utils::plot_point(env, particle_vec, color));
	}

	return weight;
}

double EihSystem::obsfunc_discrete_weight_fov_occluded(const mat &particle, const cube &image, const mat &z_buffer, bool plot) {
	double weight;
	rave::Vector color;

	if (!kinect->is_visible(particle)) {
		weight = 1;
		color = rave::Vector(0, 1, 0);
	} else {
		mat kinect_pos = rave_utils::rave_vec_to_mat(kinect->get_pose().trans);

		double particle_dist = norm(particle - kinect_pos, 2);
		mat pixel = kinect->get_pixel_from_point(particle);
		int y = int(pixel(0)), x = int(pixel(1));
		double sd = particle_dist - z_buffer(y, x);

		if (sd > -.01) {
			weight = 1;
			color = rave::Vector(0, 0, 0);
		} else {
			weight = 0;
			color = rave::Vector(1, 1, 1);
		}
	}

	if (plot) {
		rave::Vector particle_vec = rave_utils::mat_to_rave_vec(particle);
		handles.push_back(rave_utils::plot_point(env, particle_vec, color));
	}

	return weight;
}


double EihSystem::obsfunc_continuous_weight_fov_occluded_color(const mat &particle, const cube &image, const mat &z_buffer, bool plot) {
	double is_in_fov = UNKNOWN, sd_sigmoid = UNKNOWN;
	mat color = UNKNOWN*ones<mat>(3,1);

	mat pixel = kinect->get_pixel_from_point(particle);
	int y = int(pixel(0)), x = int(pixel(1));
	double alpha = 1e6;
	mat boundary;
	boundary << y << -y + kinect->get_height() << x << -x + kinect->get_width();
	is_in_fov = prod(prod(sigmoid(boundary, alpha)));

	is_in_fov = (1-.5)*is_in_fov + .5; // make being out of FOV not so advantageous for gauss likelihood

	if (is_in_fov > .75) {
		mat kinect_pos = rave_utils::rave_vec_to_mat(kinect->get_pose().trans);

		double particle_dist = norm(particle - kinect_pos, 2);
		mat pixel = kinect->get_pixel_from_point(particle);
		int y = int(pixel(0)), x = int(pixel(1));
		double sd = particle_dist - z_buffer(y, x);

		sd_sigmoid = sigmoid(sd, 10);

		if (fabs(sd) < .03) {
			mat rgb = image.subcube(y,x,0,y,x,2);
			color = utils::rgb_to_hsv(rgb);
		}

		sd_sigmoid = (1-.1)*sd_sigmoid + .1/2.0;
	}

	mat z(1, Z_DIM);
	z(0) = is_in_fov;
	z(1) = sd_sigmoid;
	z(2) = color(0);
	z(3) = color(1);
	z(4) = color(2);

	double weight = -INFINITY;
	for(int i=0; i < desired_observations.size(); ++i) {
		mat z_des = desired_observations[i];
		mat e = z - z_des;

		// deal with hue wrap-around
		double hue_small = MIN(z(2), z_des(2));
		double hue_large = MAX(z(2), z_des(2));
		e(2) = MIN(hue_large - hue_small, hue_small + (1-hue_large));

		weight = MAX(weight, gauss_likelihood(e.t(), R));
	}

	if (plot) {
		if (z(0) <= .5) { // outside FOV
			color << 0 << 1 << 0;
		} else if (z(1) < .25) { // free space
			color << 1 << 1 << 1;
		} else if (z(1) > .75) { // occluded
			color << 0 << 0 << 0;
		} else { // near surface and red
			color = image.subcube(y,x,0,y,x,2);
		}

		handles.push_back(rave_utils::plot_point(env, particle, color));
	}

	return weight;
}

double EihSystem::obsfunc_discrete_weight_fov_occluded_color(const mat &particle, const cube &image, const mat &z_buffer, bool plot) {
	double weight;
	rave::Vector color;

	if (!kinect->is_visible(particle)) {
		weight = 1;
		color = rave::Vector(0, 1, 0);
	} else {
		mat kinect_pos = rave_utils::rave_vec_to_mat(kinect->get_pose().trans);

		double particle_dist = norm(particle - kinect_pos, 2);
		mat pixel = kinect->get_pixel_from_point(particle);
		int y = int(pixel(0)), x = int(pixel(1));
		double sd = particle_dist - z_buffer(y, x);

		if (sd > .03) {
			weight = 1;
			color = rave::Vector(0, 0, 0);
		} else if (sd < -.03) {
			weight = 0;
			color = rave::Vector(1, 1, 1);
		} else {
			mat rgb = image.subcube(y,x,0,y,x,2);
			mat hsv = utils::rgb_to_hsv(rgb);

			double dist_to_red = MIN(hsv(0), 1-hsv(0));
			if (dist_to_red < .04) {
				weight = 1.1;
				color = rave_utils::mat_to_rave_vec(rgb);
			} else {
				weight = 0;
				color = rave::Vector(1, 1, 0);
			}
		}
	}

	if (plot) {
		rave::Vector particle_vec = rave_utils::mat_to_rave_vec(particle);
		handles.push_back(rave_utils::plot_point(env, particle_vec, color));
	}

	return weight;
}

void EihSystem::update_state_and_particles(const mat& x_t, const mat& P_t, const mat& u_t, mat& x_tp1, mat& P_tp1) {
	update_state_and_particles(x_t, P_t, u_t, x_tp1, P_tp1, false);
}

void EihSystem::update_state_and_particles(const mat& x_t, const mat& P_t, const mat& u_t, mat& x_tp1, mat& P_tp1, bool plot) {
	handles.clear();
	int M = P_t.n_cols;

	x_tp1 = dynfunc(x_t, u_t);
	manip->set_joint_values(x_tp1);

	cube image = kinect->get_image(true);
	mat z_buffer = kinect->get_z_buffer(false);

	mat W(M, 1, fill::zeros);
	for(int m=0; m < M; ++m) {
		W(m) = obsfunc_discrete_weight(P_t.col(m), image, z_buffer, plot);
//		W(m) = obsfunc_continuous_weight(P_t.col(m), image, z_buffer, plot);
	}

	W = W / accu(W);

	double sampling_noise = uniform(0, 1/double(M));
	P_tp1 = low_variance_sampler(P_t, W, sampling_noise);

	// add noise to each particle
	for(int m=0; m < M; ++m) {
		mat noise = .01*randu<mat>(3,1) - .005;
		P_tp1.col(m) += noise;
	}
}

double EihSystem::cost(const std::vector<mat>& X, const std::vector<mat>& U, const mat& P) {
	double entropy = 0;

	int M = P.n_cols;
	int T = U.size()+1;
	manip->set_joint_values(X[0]);

	mat W_t = (1/double(M))*ones<mat>(M, 1), W_tp1(M, 1, fill::zeros);
	for(int t=0; t < T-1; ++t) {
		mat x_tp1 = dynfunc(X[t], U[t]);
		manip->set_joint_values(x_tp1);

		cube image = kinect->get_image(true);
		mat z_buffer = kinect->get_z_buffer(false);

		for(int m=0; m < M; ++m) {
			W_tp1(m) = obsfunc_continuous_weight(P.col(m), image, z_buffer) * W_t(m);
		}
		W_tp1 = W_tp1 / accu(W_tp1);

		entropy += accu(-W_tp1 % log(W_tp1));

		W_t = W_tp1;
	}

	manip->set_joint_values(X[0]);

	return entropy;
}

double EihSystem::cost(const mat &x0, const std::vector<mat>& U, const mat& P) {
	int T = U.size()+1;
	std::vector<mat> X(T);
	X[0] = x0;
	for(int t=0; t < T-1; ++t) {
		X[t+1] = dynfunc(X[t], U[t]);
	}
	return cost(X, U, P);
}

mat EihSystem::cost_grad(std::vector<mat>& X, std::vector<mat>& U, const mat& P) {
	return System::cost_grad(X, U, P);
}

mat EihSystem::cost_grad(const mat &x0, std::vector<mat>& U, const mat& P) {
	int T = U.size()+1;
	mat g(U_DIM*(T-1), 1, fill::zeros);

	double orig, cost_p, cost_l;
	int index = 0;
	for(int t=0; t < T-1; ++t) {
		for(int i=0; i < U_DIM; ++i) {
			orig = U[t](i);

			U[t](i) = orig + step;
			cost_p = this->cost(x0, U, P);

			U[t](i) = orig - step;
			cost_l = this->cost(x0, U, P);

			U[t][i] = orig;
			g(index++) = (cost_p - cost_l)/(2*step);
		}
	}

	return g;
}

// plots kinect position
void EihSystem::display_states_and_particles(const std::vector<mat>& X, const mat& P, bool pause) {
	handles.clear();
	int M = P.n_cols;

	mat color;
	color << 0 << 0 << 1;
	for(int m=0; m < M; ++m) {
		handles.push_back(rave_utils::plot_point(env, P.col(m), color));
	}

	int T = X.size();
	double hue_start = 120.0/360, hue_end = 0.0/360;
	mat initial_joints = manip->get_joint_values();
	for(int t=0; t < T; ++t) {
		manip->set_joint_values(X[t]);
		mat hsv;
		hsv << hue_start + (hue_end - hue_start)*(t/double(T)) << 1 << 1;
		rave::Vector color = rave_utils::mat_to_rave_vec(utils::hsv_to_rgb(hsv));
		handles.push_back(rave_utils::plot_point(env, kinect->get_pose().trans, color));
		rave_utils::plot_transform(env, kinect->get_pose(), handles);
	}
	manip->set_joint_values(initial_joints);

	if (pause) {
		LOG_INFO("Press enter to continue");
		std::cin.ignore();
	}
}
