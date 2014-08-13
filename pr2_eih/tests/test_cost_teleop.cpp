#include "system/pr2_eih_system.h"

const int T = TIMESTEPS;

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

	obstacles_cam.push_back(geometry3d::Triangle({-.2,0,.75}, {-.2,.1,.75}, {-.25,.1,.75}));
	obstacles_cam.push_back(geometry3d::Triangle({-.2,0,.75}, {-.25,0,.75}, {-.25,.1,.75}));

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


int main(int argc, char* argv[]) {
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

	// setup initial state
	VectorJ j_t, j_tp1;
	j_t = arm.get_joints();
	MatrixJ j_sigma0_t = (M_PI/4)*MatrixJ::Identity(); // TODO: never update it in MPC loop

	VectorU u_t;

	const double alpha = 100;
	while(true) {
		StdVectorJ J(2, j_t);
		StdVectorU U(1, VectorU::Zero());

		arm.set_joints(j_t);

		double cost_orig = sys.cost(J, j_sigma0_t, U, obj_gaussians_t, alpha, obstacles);
		LOG_INFO("cost orig: %10.10f",cost_orig);

		sys.plot(J, obj_gaussians_t, obstacles, false);
		arm.teleop();

		j_tp1 = arm.get_joints();
		u_t = (j_tp1 - j_t)/double(DT);

		std::cout << j_t.transpose() << "\n" << j_tp1.transpose() << "\n";

		J[1] = j_tp1;
		U[0] = u_t;
		double cost = sys.cost(J, j_sigma0_t, U, obj_gaussians_t, alpha, obstacles);
		LOG_INFO("cost orig: %10.10f",cost);
		std::cout << "\n";
	}


	return 0;
}
