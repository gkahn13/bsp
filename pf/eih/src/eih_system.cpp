#include "../include/eih_system.h"

/**
 * EihSystem constructors and initializers
 */

EihSystem::EihSystem(rave::EnvironmentBasePtr e, Manipulator *m, KinectSensor *k) : env(e), manip(m), kinect(k) {
	mat R_diag;
	R_diag << .5 << .01 << .05 << 1 << 1;
	mat R = diagmat(R_diag);

	mat uMin, uMax;
	manip->get_limits(uMin, uMax);

	init(R.n_rows, R, uMin, uMax);
}

EihSystem::EihSystem(rave::EnvironmentBasePtr e, Manipulator *m, KinectSensor *k,
		int T, mat &R, mat &uMin, mat &uMax) : env(e), manip(m), kinect(k) {
	init(R.n_rows, R, uMin, uMax, T);
}

void EihSystem::init(int Z_DIM, const mat &R, const mat &uMin, const mat &uMax, int T, double DT) {
	mat j = manip->get_joint_values();
	X_DIM = MAX(j.n_rows, j.n_cols);
	U_DIM = X_DIM;
	this->Z_DIM = Z_DIM;
	this->Q_DIM = X_DIM;
	this->R_DIM = Z_DIM;
	this->TOTAL_VARS = T*X_DIM + (T-1)*U_DIM;
	this->DT = DT;
	this->R = R;

	manip->get_limits(xMin, xMax);
	this->uMin = uMin;
	this->uMax = uMax;

	desired_observations = std::vector<mat>(3);
	desired_observations[0] << 0 << UNKNOWN << UNKNOWN << UNKNOWN << UNKNOWN; // outside FOV
	desired_observations[1] << 1 << 1 << UNKNOWN << UNKNOWN << UNKNOWN; // inside FOV, occluded
	desired_observations[2] << 1 << .5 << 0 << 1 << 1; // inside FOV, not occluded, red

	kinect->power_on();
	boost::this_thread::sleep(boost::posix_time::seconds(2));
}

/**
 * EihSystem public methods
 */

mat EihSystem::dynfunc(const mat &x, const mat &u) {
	mat x_new = x + DT*u;

	for(int i=0; i < X_DIM; ++i) {
		x_new(i) = MAX(x_new(i), xMin(i));
		x_new(i) = MIN(x_new(i), xMax(i));
	}

	return x_new;
}

mat EihSystem::obsfunc(const mat& x, const mat& t, const mat& r) {

}

double EihSystem::obsfunc_continuous_weight(const mat &particle, const cube &image, const mat &z_buffer) {
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

	// plot
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

	return weight;
}

double EihSystem::obsfunc_discrete_weight(const mat &particle, const cube &image, const mat &z_buffer) {
	double weight;
	rave::Vector color;

	if (!kinect->is_in_fov(particle)) {
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

	rave::Vector particle_vec = rave_utils::mat_to_rave_vec(particle);
	handles.push_back(rave_utils::plot_point(env, particle_vec, color));

	return weight;
}

void EihSystem::update_state_and_particles(const mat& x_t, const mat& P_t, const mat& u_t, mat& x_tp1, mat& P_tp1) {
	handles.clear();
	int M = P_t.n_cols;

	x_tp1 = dynfunc(x_t, u_t);
	manip->set_joint_values(x_tp1);

	cube image = kinect->get_image();
	mat z_buffer = kinect->get_z_buffer();

	mat W(M, 1, fill::zeros);
	for(int m=0; m < M; ++m) {
//		W(m) = obsfunc_discrete_weight(P_t.col(m), image, z_buffer);
		W(m) = obsfunc_continuous_weight(P_t.col(m), image, z_buffer);
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

}

double EihSystem::cost(const mat &x0, const std::vector<mat>& U, const mat& P) {
	double entropy = 0;

	int M = P.n_cols;
	int T = U.size()+1;
	manip->set_joint_values(x0);

	mat x_t = x0, x_tp1;
	mat W_t = (1/double(M))*ones<mat>(M, 1), W_tp1(M, 1, fill::zeros);
	for(int t=0; t < T-1; ++t) {
		x_tp1 = dynfunc(x_t, U[t]);
		manip->set_joint_values(x_tp1);

		cube image = kinect->get_image();
		mat z_buffer = kinect->get_z_buffer();

		for(int m=0; m < M; ++m) {
			W_tp1(m) = obsfunc_continuous_weight(P.col(m), image, z_buffer) * W_t(m);
		}

		entropy += accu(-W_tp1 % log(W_tp1));

		W_t = W_tp1;
	}

	manip->set_joint_values(x0);

	return entropy;
}

void EihSystem::display_states_and_particles(const std::vector<mat>& X, const mat& P, bool pause) {
	handles.clear();
	int M = P.n_cols;

	mat color;
	color << 0 << 1 << 0;
	for(int m=0; m < M; ++m) {
		handles.push_back(rave_utils::plot_point(env, P.col(m), color));
	}

	int T = X.size();
	double hue_start = 180.0/360, hue_end = 240.0/360;
	mat initial_joints = manip->get_joint_values();
	for(int t=0; t < T; ++t) {
		manip->set_joint_values(X[t]);
		rave::Vector color = rave::Vector(hue_start + (hue_end - hue_start)*(t/T), 1, 1);
		handles.push_back(rave_utils::plot_point(env, manip->get_pose().trans, color));
	}
	manip->set_joint_values(initial_joints);
}
