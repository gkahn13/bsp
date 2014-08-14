#include "pr2_eih_sqp.h"
#include "../util/logging.h"

#include "pr2_utils/pr2/arm.h"
#include "pcl_utils/OccludedRegionArray.h"

class PR2EihBrett {
public:
	PR2EihBrett();

	void get_occluded_regions(std::vector<Gaussian3d>& obj_gaussians,
			std::vector<geometry3d::Triangle>& obstacles);
	void initialize_trajectory(StdVectorJ& J, StdVectorU& U, const std::vector<Gaussian3d>& obj_gaussians);
	void bsp(StdVectorJ& J, StdVectorU& U, const MatrixJ& j_sigma0,
			const std::vector<Gaussian3d>& obj_gaussians, const std::vector<geometry3d::Triangle>& obstacles);
	void execute_controls(const StdVectorU& U);

private:
	// ros
	ros::NodeHandle *nh_ptr;
	ros::Subscriber occluded_region_array_sub;

	bool received_occluded_region_array;
	pcl_utils::OccludedRegionArray current_occluded_region_array;

	// pr2_utils
	pr2_sim::Simulator *sim;
	pr2_sim::Arm *arm_sim;
	pr2_sim::Camera *cam_sim;
	PR2EihSystem *sys;
	pr2::Arm *arm;

	// FORCES
	pr2eihMPC_params problem;
	pr2eihMPC_output output;
	pr2eihMPC_info info;

	void _occluded_region_array_callback(const pcl_utils::OccludedRegionArrayConstPtr& msg);
};

PR2EihBrett::PR2EihBrett() {
	// setup system
	sim = new pr2_sim::Simulator(true, false);
	arm_sim = new pr2_sim::Arm(pr2_sim::Arm::right, sim);
	cam_sim = new pr2_sim::Camera(arm_sim, sim);
	sys = new PR2EihSystem(sim, arm_sim, cam_sim);

	// setup actual arm
	arm = new pr2::Arm(pr2::Arm::ArmType::right);

	// subscribe to OccludedRegionArray topic
	nh_ptr = new ros::NodeHandle();
	occluded_region_array_sub =
			nh_ptr->subscribe("/kinfu/occluded_region_array", 1, &PR2EihBrett::_occluded_region_array_callback, this);
	received_occluded_region_array = false;

	// setup FORCES
	setup_mpc_vars(problem, output);
}

void PR2EihBrett::get_occluded_regions(std::vector<Gaussian3d>& obj_gaussians,
		std::vector<geometry3d::Triangle>& obstacles) {
	obj_gaussians.clear();
	obstacles.clear();

	received_occluded_region_array = false;
	ROS_INFO("Waiting for new occluded region...");
	while (!received_occluded_region_array && !ros::isShuttingDown()) {
		ros::spinOnce();
	}

	sim->update();
	Matrix4d cam_pose = cam_sim->get_pose(arm->get_joints());
	Matrix3d cam_rot = cam_pose.block<3,3>(0,0);
	Vector3d cam_trans = cam_pose.block<3,1>(0,3);

	pcl_utils::OccludedRegionArray occluded_region_array = current_occluded_region_array;
	for(const pcl_utils::OccludedRegion& occluded_region : occluded_region_array.regions) {
		// add occlusions by forming triangles from front face polygon
		const std::vector<geometry_msgs::Point32>& points = occluded_region.front_face.points;
		Vector3d first_obstacle_point(points[0].x, points[0].y, points[0].z);
		for(int i=1; i < points.size()-1; ++i) {
			Vector3d second_obstacle_point(points[i].x, points[i].y, points[i].z);
			Vector3d third_obstacle_point(points[i+1].x, points[i+1].y, points[i+1].z);
			obstacles.push_back(geometry3d::Triangle(cam_rot*first_obstacle_point + cam_trans,
					cam_rot*second_obstacle_point + cam_trans, cam_rot*third_obstacle_point + cam_trans));
		}

		// add gaussian
		const pcl_utils::Gaussian& gaussian = occluded_region.gaussian;
		Vector3d mean(gaussian.mean.x, gaussian.mean.y, gaussian.mean.z);
		Matrix3d cov(gaussian.covariance.data());
		obj_gaussians.push_back(Gaussian3d(mean, cov));
	}
}

void PR2EihBrett::initialize_trajectory(StdVectorJ& J, StdVectorU& U, const std::vector<Gaussian3d>& obj_gaussians) {
	Vector3d avg_obj_mean = Vector3d::Zero();
	double num_objs = obj_gaussians.size();
	for(const Gaussian3d& obj_gaussian : obj_gaussians) {
		avg_obj_mean += (1/num_objs)*obj_gaussian.mean;
	}

	J[0] = arm->get_joints();
	Vector3d start_position = arm_sim->fk(J[0]).block<3,1>(0,3);
	Vector3d next_position;
	VectorJ next_joints;
	for(int t=0; t < T-1; ++t) {
		next_position = start_position + (t+1)*Vector3d(.05,0,.05);
		if (arm_sim->ik_lookat(next_position, avg_obj_mean, next_joints)) {
			U[t] = next_joints - J[t];
		} else {
			U[t] = VectorU::Zero();
		}
//		U[t].setZero();

		J[t+1] = sys->dynfunc(J[t], U[t], VectorQ::Zero(), true);
		U[t] = (J[t+1] - J[t])/double(DT);
	}
}

void PR2EihBrett::bsp(StdVectorJ& J, StdVectorU& U, const MatrixJ& j_sigma0,
		const std::vector<Gaussian3d>& obj_gaussians, const std::vector<geometry3d::Triangle>& obstacles) {
	// plan
	pr2_eih_collocation(J, U, j_sigma0, obj_gaussians, obstacles, *sys, problem, output, info);

	// reintegrate, just in case
	for(int t=0; t < T; ++t) {
		J[t+1] = sys->dynfunc(J[t], U[t], VectorQ::Zero(), true);
	}
}

void PR2EihBrett::execute_controls(const StdVectorU& U) {
	VectorJ current_joints = arm->get_joints();

	std::vector<VectorJ> joint_traj(U.size());
	for(int t=0; t < U.size(); ++t) {
		current_joints = sys->dynfunc(current_joints, U[t], VectorQ::Zero(), true);
		joint_traj[t] = current_joints;
	}

	arm->execute_joint_trajectory(joint_traj);
}

void PR2EihBrett::_occluded_region_array_callback(const pcl_utils::OccludedRegionArrayConstPtr& msg) {
	current_occluded_region_array = *msg;
	received_occluded_region_array = true;
}

int main(int argc, char* argv[]) {
	// initialize ros node
	ros::init(argc, argv, "robot_driver");
	log4cxx::LoggerPtr my_logger = log4cxx::Logger::getLogger(ROSCONSOLE_DEFAULT_NAME);
	my_logger->setLevel(ros::console::g_level_lookup[ros::console::levels::Info]);

	PR2EihBrett brett_bsp;

	std::vector<Gaussian3d> obj_gaussians;
	std::vector<geometry3d::Triangle> obstacles;

	MatrixJ j_sigma0 = (M_PI/4)*MatrixJ::Identity(); // TODO: never update in MPC
	StdVectorJ J(T);
	StdVectorU U(T-1);

	while(!ros::isShuttingDown()) {
		ROS_INFO("Getting occluded regions");
		brett_bsp.get_occluded_regions(obj_gaussians, obstacles);

		ROS_INFO("Initializing trajectory");
		brett_bsp.initialize_trajectory(J, U, obj_gaussians);

		ROS_INFO("Optimizing trajectory with bsp");
		brett_bsp.bsp(J, U, j_sigma0, obj_gaussians, obstacles);

		ROS_INFO("Executing first control");
		brett_bsp.execute_controls(StdVectorU(1, U[0]));
	}

}
