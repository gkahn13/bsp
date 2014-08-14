//#include "pr2_eih_sqp.h"
#include "sqp/pr2_eih_sqp.h"
#include "../util/logging.h"

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

	obstacles_cam.push_back(geometry3d::Triangle({0,.1,.75}, {0,.2,.75}, {.05,.2,.75}));
	obstacles_cam.push_back(geometry3d::Triangle({0,.1,.75}, {.05,.1,.75}, {.05,.2,.75}));

	obstacles_cam.push_back(geometry3d::Triangle({-.2,.1,.75}, {-.2,.2,.75}, {-.25,.2,.75}));
	obstacles_cam.push_back(geometry3d::Triangle({-.2,.1,.75}, {-.25,.1,.75}, {-.25,.2,.75}));

	obstacles.clear();
	for(const geometry3d::Triangle& obstacle_cam : obstacles_cam) {
		obstacles.push_back(geometry3d::Triangle(cam_rot*obstacle_cam.a+cam_pos,
				cam_rot*obstacle_cam.b+cam_pos,
				cam_rot*obstacle_cam.c+cam_pos));
	}

//	obstacles.push_back(geometry3d::Triangle({0,2,.45}, {0,-2,.45}, {5,0,.45}));

	std::vector<geometry3d::TruncatedPyramid> truncated_frustum = cam.truncated_view_frustum(cam_pose, obstacles, true);
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

void init_obstacles_and_objects_from_data(pr2_sim::Camera& cam,
		std::vector<geometry3d::Triangle>& obstacles, std::vector<Gaussian3d>& obj_gaussians) {
	Matrix4d cam_pose = cam.get_pose();
	Matrix3d cam_rot = cam_pose.block<3,3>(0,0);
	Vector3d cam_pos = cam_pose.block<3,1>(0,3);
	std::vector<geometry3d::Triangle> obstacles_cam;

	Vector3d obj0_p0(-0.117518998682, -0.167983621359, 0.738430976868);
	Vector3d obj0_p1(0.0714884921908, -0.152321845293, 0.766414642334);
	Vector3d obj0_p2(0.0392220690846, 0.0534862950444, 0.869163274765);
	Vector3d obj0_p3(-0.149785429239, 0.0378245040774, 0.841179609299);
	obstacles_cam.push_back(geometry3d::Triangle(obj0_p0, obj0_p1, obj0_p2));
	obstacles_cam.push_back(geometry3d::Triangle(obj0_p0, obj0_p2, obj0_p3));

	Vector3d obj1_p0(-0.0622684955597, 0.0703296512365, 0.85747385025);
	Vector3d obj1_p1(-0.0780080407858, -0.0499060563743, 0.815235137939);
	Vector3d obj1_p2(-0.277430295944, 0.0255430340767, 0.674774885178);
	Vector3d obj1_p3(-0.261690735817, 0.145778745413, 0.717013597488);
	obstacles_cam.push_back(geometry3d::Triangle(obj1_p0, obj1_p1, obj1_p2));
	obstacles_cam.push_back(geometry3d::Triangle(obj1_p0, obj1_p2, obj1_p3));

	Vector3d obj2_p0(0.0894662365317, 0.0748614519835, 0.67870414257);
	Vector3d obj2_p1(0.0901472344995, 0.100558131933, 0.678196072578);
	Vector3d obj2_p2(-0.0529698356986, 0.104965731502, 0.709290742874);
	Vector3d obj2_p3(-0.0536508336663, 0.0792690515518, 0.709798812866);
	obstacles_cam.push_back(geometry3d::Triangle(obj2_p0, obj2_p1, obj2_p2));
	obstacles_cam.push_back(geometry3d::Triangle(obj2_p0, obj2_p2, obj2_p3));

	obstacles.clear();
	for(const geometry3d::Triangle& obstacle_cam : obstacles_cam) {
		obstacles.push_back(geometry3d::Triangle(cam_rot*obstacle_cam.a+cam_pos,
				cam_rot*obstacle_cam.b+cam_pos,
				cam_rot*obstacle_cam.c+cam_pos));
	}

	obj_gaussians.clear();
	Vector3d mean0(-0.0874505236686, -0.189941832512, 1.04621961403);
	Matrix3d cov0;
	cov0 << 0.0012533777221510694, 0.0011764044410806562, -0.00274903438138626,
			0.0011764044410806562, 0.008869640353479464, -0.012381557788010235,
			-0.00274903438138626, -0.012381557788010235, 0.02644106012457171;
	mean0 = cam_rot*mean0 + cam_pos;
	cov0 = cam_rot*cov0*cam_rot.transpose();
	obj_gaussians.push_back(Gaussian3d(mean0, cov0));

	Vector3d mean1(-0.293816802324, 0.0243655858351, 0.899660059232);
	Matrix3d cov1;
	cov1 << 0.004818769168691134, 0.0004966072557262238, -0.002970275982414134,
			0.0004966072557262238, 0.0010962437321780573, -0.0014799552310396417,
			-0.002970275982414134, -0.0014799552310396417, 0.005960974123895709;
	mean1 = cam_rot*mean1 + cam_pos;
	cov1 = cam_rot*cov1*cam_rot.transpose();
	obj_gaussians.push_back(Gaussian3d(mean1, cov1));


	Vector3d mean2(0.0228484183215, 0.0911101793249, 0.748269921541);
	Matrix3d cov2;
	cov2 << 0.0017403703895459083, -5.018270935997982e-05, -0.0002676494784628474,
			-5.018270935997982e-05, 4.81392402055504e-05, 1.7359733157569226e-05,
			-0.0002676494784628474, 1.7359733157569226e-05, 0.0005657157674245994;
	mean2 = cam_rot*mean2 + cam_pos;
	cov2 = cam_rot*cov2*cam_rot.transpose();
	obj_gaussians.push_back(Gaussian3d(mean2, cov2));
}

void init_trajectory(StdVectorJ& J, StdVectorU& U, const std::vector<Gaussian3d>& obj_gaussians,
		pr2_sim::Arm& arm, PR2EihSystem& sys) {
	Vector3d avg_obj_mean = Vector3d::Zero();
	double num_objs = obj_gaussians.size();
	for(const Gaussian3d& obj_gaussian : obj_gaussians) {
		avg_obj_mean += (1/num_objs)*obj_gaussian.mean;
	}

	Vector3d start_position = arm.fk(J[0]).block<3,1>(0,3);
	Vector3d next_position;
	VectorJ next_joints;
	for(int t=0; t < T-1; ++t) {
		next_position = start_position + (t+1)*Vector3d(.05,0,.05);
		if (arm.ik_lookat(next_position, avg_obj_mean, next_joints)) {
			U[t] = next_joints - J[t];
		} else {
			U[t] = VectorU::Zero();
		}
//		U[t].setZero(); // TODO :temp

		J[t+1] = sys.dynfunc(J[t], U[t], VectorQ::Zero(), true);
		U[t] = (J[t+1] - J[t])/double(DT);
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
//	init_obstacles_and_objects(cam, obstacles, obj_gaussians_t);
	init_obstacles_and_objects_from_data(cam, obstacles, obj_gaussians_t);

	// setup initial state
	VectorJ j_t, j_t_real, j_tp1, j_tp1_real;
	j_t = arm.get_joints();
	j_t_real = j_t; // TODO
	MatrixJ j_sigma0_t = (M_PI/4)*MatrixJ::Identity(); // TODO: never update it in MPC loop

	// initialize state and controls
	StdVectorU U(T-1, VectorU::Zero());
	StdVectorJ J(T, j_t);

	pr2_eih_sqp::PR2EihSqp pr2_eih_bsp;

	util::Timer forces_timer;

	while(true) {
		init_trajectory(J, U, obj_gaussians_t, arm, sys);

		double initial_cost = sys.cost(J, j_sigma0_t, U, obj_gaussians_t, INFINITY, obstacles);
		LOG_INFO("Initial cost: %4.5f", initial_cost);

		LOG_INFO("Current state");
		sys.plot(StdVectorJ(1, J[0]), obj_gaussians_t, obstacles);

		// optimize
		util::Timer_tic(&forces_timer);
		double cost = pr2_eih_bsp.collocation(J, U, j_sigma0_t, obj_gaussians_t, obstacles, sys);
		double forces_time = util::Timer_toc(&forces_timer);

		LOG_INFO("Optimized cost: %4.5f", cost);
		LOG_INFO("Solve time: %5.3f ms", 1e3*forces_time);

		for(int t=0; t < T-1; ++t) {
			J[t+1] = sys.dynfunc(J[t], U[t], VectorQ::Zero());
		}

		LOG_INFO("Displaying optimized trajectory");
		sys.plot(J, obj_gaussians_t, obstacles);

		sys.execute_control_step(j_t_real, j_t, U[0], obj_gaussians_t, obstacles, j_tp1_real, j_tp1, obj_gaussians_tp1);

		j_t_real = j_tp1_real;
		j_t = j_tp1;
		obj_gaussians_t = obj_gaussians_tp1;

		// TODO: replace with optimization vars from [1, ..., T-1]
		J = StdVectorJ(T, j_t);
		U = StdVectorU(T-1, VectorU::Zero());

	}

	return 0;
}
