#include "system/pr2_eih_system.h"

const double alpha = 100;

MatrixXd low_variance_sampler(const MatrixXd& P, const VectorXd& W) {
	int M = P.cols();
	MatrixXd P_sampled(3,M);

	double r = (1/double(M))*(rand() / double(RAND_MAX));
	double c = W(0);
	int i = 0;
	for(int m=0; m < M; ++m) {
		double u = r + m * (1/double(M));
		while (u > c) {
			c += W(++i);
		}
		P_sampled.col(m) = P.col(i);
	}

	return P_sampled;
}


VectorXd update_particle_weights(const VectorJ& j, const MatrixXd& P_t, const VectorXd& W_t, const std::vector<geometry3d::Triangle>& obstacles,
		pr2_sim::Arm& arm, pr2_sim::Camera& cam) {
	int M = P_t.cols();
	arm.set_joints(j);
	std::vector<geometry3d::Pyramid> truncated_frustum = cam.truncated_view_frustum(obstacles, false);

	VectorXd W_tp1(M);
	// for each particle, weight by sigmoid of signed distance
	for(int m=0; m < M; ++m) {
		double sd = cam.signed_distance(P_t.col(m), truncated_frustum);
		double sigmoid_sd = 1.0/(1.0 + exp(-alpha*sd));
		W_tp1(m) = W_t(m)*sigmoid_sd;
	}
	W_tp1 = W_tp1 / W_tp1.sum();

	return W_tp1;
}

double entropy(const VectorJ& j, const MatrixXd& P, const std::vector<geometry3d::Triangle>& obstacles,
		pr2_sim::Arm& arm, pr2_sim::Camera& cam) {
	int M = P.cols();

	VectorXd W_t = (1/double(M))*VectorXd::Ones(M);
	VectorXd W_tp1 = update_particle_weights(j, P, W_t, obstacles, arm, cam);

	double entropy = (-W_tp1.array() * W_tp1.array().log()).sum();

	return entropy;
}

VectorJ entropy_grad(const VectorJ& j, const MatrixXd& P, const std::vector<geometry3d::Triangle>& obstacles,
		pr2_sim::Arm& arm, pr2_sim::Camera& cam) {
	const double step = 1e-5;
	VectorJ grad;

	VectorU j_p = j, j_m = j;
	for(int i=0; i < J_DIM; ++i) {
		j_p(i) = j(i) + step;
		j_m(i) = j(i) - step;

		double entropy_p = entropy(j_p, P, obstacles, arm, cam);
		double entropy_m = entropy(j_m, P, obstacles, arm, cam);
		grad(i) = (entropy_p - entropy_m) / (2*step);

		j_p(i) = j(i);
		j_m(i) = j(i);
	}

	return grad;
}

inline double uniform(double low, double high) {
	return (high - low)*(rand() / double(RAND_MAX)) + low;
}

void init_obstacles_and_objects(pr2_sim::Camera& cam,
		std::vector<geometry3d::Triangle>& obstacles, std::vector<Gaussian3d>& obj_gaussians) {
	Matrix4d cam_pose = cam.get_pose();
	Matrix3d cam_rot = cam_pose.block<3,3>(0,0);
	Vector3d cam_pos = cam_pose.block<3,1>(0,3);

	std::vector<geometry3d::Triangle> obstacles_cam;

	obstacles_cam.push_back(geometry3d::Triangle({0,0,.75}, {0,.1,.75}, {.05,.1,.75}));
	obstacles_cam.push_back(geometry3d::Triangle({0,0,.75}, {.05,0,.75}, {.05,.1,.75}));

	obstacles_cam.push_back(geometry3d::Triangle({-.2,-.2,.75}, {-.2,-.1,.75}, {-.25,-.1,.75}));
	obstacles_cam.push_back(geometry3d::Triangle({-.2,-.2,.75}, {-.25,-.2,.75}, {-.25,-.1,.75}));

	obstacles.clear();
	for(const geometry3d::Triangle& obstacle_cam : obstacles_cam) {
		obstacles.push_back(geometry3d::Triangle(cam_rot*obstacle_cam.a+cam_pos,
				cam_rot*obstacle_cam.b+cam_pos,
				cam_rot*obstacle_cam.c+cam_pos));
	}

	std::vector<geometry3d::Pyramid> truncated_frustum = cam.truncated_view_frustum(obstacles, true);
	obj_gaussians.clear();
	for(int i=0; i < obstacles_cam.size(); i+=2) {
		geometry3d::Triangle& obstacle_cam = obstacles_cam[i];
		double x_min = std::min(obstacle_cam.a(0), std::min(obstacle_cam.b(0), obstacle_cam.c(0))) - .1;
		double x_max = std::max(obstacle_cam.a(0), std::max(obstacle_cam.b(0), obstacle_cam.c(0))) + .1;
		double y_min = std::min(obstacle_cam.a(1), std::min(obstacle_cam.b(1), obstacle_cam.c(1))) - .1;
		double y_max = std::max(obstacle_cam.a(1), std::max(obstacle_cam.b(1), obstacle_cam.c(1))) + .1;
		double z_min = std::max(obstacle_cam.a(2), std::max(obstacle_cam.b(2), obstacle_cam.c(2)));
		double z_max = z_min + 2*fabs(y_max - y_min);

		int num_particles = 0;
		MatrixP particles;
		while(num_particles < M_DIM) {
			Vector3d p_cam = Vector3d(uniform(x_min, x_max), uniform(y_min, y_max), uniform(z_min, z_max));
			Vector3d p_world = cam_rot*p_cam + cam_pos;
			if (!cam.is_in_fov(p_world, truncated_frustum)) {
				particles.col(num_particles++) = p_world;
			}
		}
		obj_gaussians.push_back(Gaussian3d(particles));
	}
}


int main() {
	// setup system
	pr2_sim::Simulator sim(true, false);
	pr2_sim::Arm arm(pr2_sim::Arm::right, &sim);
	pr2_sim::Camera cam(&arm, &sim);
	PR2EihSystem sys(&sim, &arm, &cam);

	pr2_sim::Arm other_arm(pr2_sim::Arm::left, &sim);
	other_arm.set_posture(pr2_sim::Arm::Posture::mantis);
	arm.set_posture(pr2_sim::Arm::Posture::mantis);

	// setup environment
	std::vector<geometry3d::Triangle> obstacles;
	std::vector<Gaussian3d> obj_gaussians;
	init_obstacles_and_objects(cam, obstacles, obj_gaussians);

	// get all particles into one batch
	int M_TOTAL = M_DIM*obj_gaussians.size();
	MatrixXd P_t(3,M_TOTAL), P_tp1(3,M_TOTAL);
	for(int i=0; i < obj_gaussians.size(); ++i) {
		P_t.block<3,M_DIM>(0, i*M_DIM) = obj_gaussians[i].particles;
	}

	VectorJ j_t, j_t_real, j_tp1, j_tp1_real;
	j_t = arm.get_joints();
	j_t_real = j_t;

	VectorU u_t;

	while(true) {
		LOG_INFO("Initial state");
		sys.plot(StdVectorJ(1, j_t), obj_gaussians, obstacles);

		VectorU grad = entropy_grad(j_t, P_t, obstacles, arm, cam);
		std::cout << "grad: " << grad.transpose() << "\n";

		j_tp1 = j_t - (M_PI/16.)*(grad/grad.norm());

		LOG_INFO("After grad step");
		sys.plot(StdVectorJ(1, j_tp1), obj_gaussians, obstacles);

		VectorXd W_tp1 = update_particle_weights(j_tp1, P_t, (1/double(M_TOTAL))*VectorXd::Ones(M_TOTAL), obstacles, arm, cam);
		P_tp1 = low_variance_sampler(P_t, W_tp1);

		for(int i=0; i < obj_gaussians.size(); ++i) {
			obj_gaussians[i] = Gaussian3d(P_tp1.block<3,M_DIM>(0,i*M_DIM));
		}

		j_t = j_tp1;
		P_t = P_tp1;
		arm.set_joints(j_t);
	}

}
