#include "system/pr2_eih_system.h"

double cost(const VectorJ& j, const VectorU& u, const std::vector<Gaussian3d>& obj_gaussians,
		const std::vector<geometry3d::Triangle>& obstacles, pr2_sim::Arm& arm, pr2_sim::Camera& cam, PR2EihSystem& sys) {
	double cost = 0;

	VectorJ j_tp1 = sys.dynfunc(j, u, VectorQ::Zero());
	arm.set_joints(j_tp1);

	std::vector<geometry3d::TruncatedPyramid> truncated_frustum = cam.truncated_view_frustum(obstacles, false);
	for(const Gaussian3d& obj_gaussian : obj_gaussians) {
		const MatrixP& particles = obj_gaussian.particles;
		for(int m=0; m < M_DIM; ++m) {
			cost += cam.signed_distance(particles.col(m), truncated_frustum);
		}
	}

	return cost;
}

VectorU cost_grad(const VectorJ& j, const VectorU& u, const std::vector<Gaussian3d>& obj_gaussians,
		const std::vector<geometry3d::Triangle>& obstacles, pr2_sim::Arm& arm, pr2_sim::Camera& cam, PR2EihSystem& sys) {
	const double step = 1e-5;
	VectorU grad;

	VectorU u_p = u, u_m = u;
	for(int i=0; i < j.rows(); ++i) {
		u_p(i) += step;
		u_m(i) -= step;

		double cost_p = cost(j, u_p, obj_gaussians, obstacles, arm, cam, sys);
		double cost_m = cost(j, u_m, obj_gaussians, obstacles, arm, cam, sys);
		grad(i) = (cost_p - cost_m) / (2*step);

		u_p(i) = u(i);
		u_m(i) = u(i);
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

	std::vector<geometry3d::TruncatedPyramid> truncated_frustum = cam.truncated_view_frustum(obstacles, true);
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
	std::vector<Gaussian3d> obj_gaussians_t, obj_gaussians_tp1;
	init_obstacles_and_objects(cam, obstacles, obj_gaussians_t);

	VectorJ j_t, j_t_real, j_tp1, j_tp1_real;
	j_t = arm.get_joints();
	j_t_real = j_t;

	VectorU u_t;

	while(true) {
		LOG_INFO("Initial state");
		sys.plot(StdVectorJ(1, j_t), obj_gaussians_t, obstacles);

		VectorU grad = cost_grad(j_t, VectorU::Zero(), obj_gaussians_t, obstacles, arm, cam, sys);
		std::cout << "grad: " << grad.transpose() << "\n";

		u_t = -(M_PI/16.)*(grad/grad.norm());

		LOG_INFO("After grad step");
		sys.plot(StdVectorJ(1, j_t+u_t), obj_gaussians_t, obstacles);

		sys.execute_control_step(j_t_real, j_t, u_t, obj_gaussians_t, obstacles, j_tp1_real, j_tp1, obj_gaussians_tp1);

		j_t = j_tp1;
		j_t_real = j_tp1_real;
		obj_gaussians_t = obj_gaussians_tp1;
	}
}
