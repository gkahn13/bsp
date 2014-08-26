#include "sqp/pr2_eih_sqp.h"
#include "../util/logging.h"

#include "pr2_utils/pr2/arm.h"
#include "pcl_utils/OccludedRegionArray.h"

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <geometry_msgs/PoseArray.h>
#include <std_msgs/Empty.h>
#include <visualization_msgs/MarkerArray.h>

const int T = TIMESTEPS;

class OccludedRegionCompare {
public:
	OccludedRegionCompare() { }

	// assumes a and b are in the camera frame
	bool operator()(const pcl_utils::OccludedRegion& a, const pcl_utils::OccludedRegion& b) const {
		Vector3d a_gaussian_mean(a.gaussian.mean.x, a.gaussian.mean.y, a.gaussian.mean.z);
		Vector3d b_gaussian_mean(b.gaussian.mean.x, b.gaussian.mean.y, b.gaussian.mean.z);
		return (a_gaussian_mean.norm() < b_gaussian_mean.norm());
	}
};

class PR2EihBrett {
public:
	PR2EihBrett(int max_occluded_regions, const Vector3d& init_traj_control, double max_travel_distance);
	PR2EihBrett() : PR2EihBrett(1, Vector3d(0,0,0), 0.1) { }

	void get_occluded_regions(std::vector<Gaussian3d>& obj_gaussians,
			std::vector<geometry3d::Triangle>& obstacles);
	void display_gaussian_means(const std::vector<Gaussian3d>& obj_gaussians);
	void initialize_trajectory(StdVectorJ& J, StdVectorU& U, const std::vector<Gaussian3d>& obj_gaussians);
	void bsp(StdVectorJ& J, StdVectorU& U, const MatrixJ& j_sigma0,
			const std::vector<Gaussian3d>& obj_gaussians, const std::vector<geometry3d::Triangle>& obstacles);
	void display_trajectory(const StdVectorJ& J);
	void execute_controls(const StdVectorU& U);

private:
	// ros
	ros::NodeHandle *nh_ptr;
	ros::Subscriber occluded_region_array_sub;
	ros::Publisher display_trajectory_pub, get_occluded_pub, display_gaussian_pub;

	bool received_occluded_region_array;
	pcl_utils::OccludedRegionArray current_occluded_region_array;

	int max_occluded_regions;
	Vector3d init_traj_control;
	double max_travel_distance;

	// pr2_utils
	pr2_sim::Simulator *sim;
	pr2_sim::Arm *arm_sim;
	pr2_sim::Camera *cam_sim;
	PR2EihSystem *sys;
	pr2::Arm *arm;

	// sqp solver
	pr2_eih_sqp::PR2EihSqp pr2_eih_bsp;

	void _occluded_region_array_callback(const pcl_utils::OccludedRegionArrayConstPtr& msg);
};

PR2EihBrett::PR2EihBrett(int max_occluded_regions, const Vector3d& init_traj_control, double max_travel_distance) :
		max_occluded_regions(max_occluded_regions), init_traj_control(init_traj_control), max_travel_distance(max_travel_distance) {
	// setup system
	ROS_INFO("Setting up simulated system");
	sim = new pr2_sim::Simulator(true, true);
	arm_sim = new pr2_sim::Arm(pr2_sim::Arm::right, sim);
	cam_sim = new pr2_sim::Camera(arm_sim, sim);
	sys = new PR2EihSystem(sim, arm_sim, cam_sim);

	// setup actual arm
	ROS_INFO("Setting up actual arm");
	arm = new pr2::Arm(pr2::Arm::ArmType::right);

	// subscribe to OccludedRegionArray topic
	ROS_INFO("Creating publishers/subscribers");
	nh_ptr = new ros::NodeHandle();
	occluded_region_array_sub =
			nh_ptr->subscribe("/kinfu/occluded_regions", 1, &PR2EihBrett::_occluded_region_array_callback, this);
	received_occluded_region_array = false;
	display_gaussian_pub = nh_ptr->advertise<visualization_msgs::MarkerArray>("/bsp_object_means", 1);
	display_trajectory_pub = nh_ptr->advertise<geometry_msgs::PoseArray>("/bsp_trajectory", 1);
	get_occluded_pub = nh_ptr->advertise<std_msgs::Empty>("/get_occluded", 1);

	sim->update();
}

void PR2EihBrett::get_occluded_regions(std::vector<Gaussian3d>& obj_gaussians,
		std::vector<geometry3d::Triangle>& obstacles) {
	obj_gaussians.clear();
	obstacles.clear();

	received_occluded_region_array = false;
	get_occluded_pub.publish(std_msgs::Empty());
	ROS_INFO("Waiting for new occluded region...");
	while (!received_occluded_region_array && !ros::isShuttingDown()) {
		ros::spinOnce();
	}

	sim->update();
	Matrix4d cam_pose = cam_sim->get_pose(arm->get_joints());
	Matrix3d cam_rot = cam_pose.block<3,3>(0,0);
	Vector3d cam_trans = cam_pose.block<3,1>(0,3);

	std::vector<pcl_utils::OccludedRegion> occluded_regions = current_occluded_region_array.regions;
	std::sort(occluded_regions.begin(), occluded_regions.end(), OccludedRegionCompare());

	int num_occluded_regions = std::min(max_occluded_regions, int(occluded_regions.size()));
	ROS_INFO_STREAM("Getting the closest " << num_occluded_regions << " occluded regions");
	for(int i=0; i < num_occluded_regions; ++i) {
		const pcl_utils::OccludedRegion& occluded_region = occluded_regions[i];
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
		obj_gaussians.push_back(Gaussian3d(cam_rot*mean+cam_trans, cam_rot*cov*cam_rot.transpose()));

		ROS_INFO_STREAM("Gaussian " << i << "\nMean: " << obj_gaussians.back().mean.transpose()
				<< "\nCovariance:\n" << obj_gaussians.back().cov << "\n");
	}
}

void PR2EihBrett::display_gaussian_means(const std::vector<Gaussian3d>& obj_gaussians) {
	visualization_msgs::MarkerArray marker_array;

	marker_array.markers.resize(obj_gaussians.size());
	for(int i=0; i < obj_gaussians.size(); ++i) {
		marker_array.markers[i].header.frame_id = "/base_link";
		marker_array.markers[i].header.stamp = ros::Time::now();
		marker_array.markers[i].id = i + 1000;
		marker_array.markers[i].type = visualization_msgs::Marker::CUBE;
		marker_array.markers[i].action = visualization_msgs::Marker::ADD;
		marker_array.markers[i].pose.position.x = obj_gaussians[i].mean(0);
		marker_array.markers[i].pose.position.y = obj_gaussians[i].mean(1);
		marker_array.markers[i].pose.position.z = obj_gaussians[i].mean(2);
		marker_array.markers[i].pose.orientation.x = 0;
		marker_array.markers[i].pose.orientation.y = 0;
		marker_array.markers[i].pose.orientation.z = 0;
		marker_array.markers[i].pose.orientation.w = 1;
		marker_array.markers[i].scale.x = .02;
		marker_array.markers[i].scale.y = .02;
		marker_array.markers[i].scale.z = .02;
		marker_array.markers[i].color.a = 1;
		marker_array.markers[i].color.r = 1.0;
		marker_array.markers[i].color.g = 0.0;
		marker_array.markers[i].color.b = 1.0;
	}

	display_gaussian_pub.publish(marker_array);
}

void PR2EihBrett::initialize_trajectory(StdVectorJ& J, StdVectorU& U, const std::vector<Gaussian3d>& obj_gaussians) {
	ROS_INFO_STREAM("Delta position for each timestep is: " << init_traj_control.transpose());

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
		next_position = start_position + (t+1)*init_traj_control;
		if (arm_sim->ik_lookat(next_position, avg_obj_mean, next_joints)) {
			U[t] = next_joints - J[t];
		} else {
			U[t] = VectorU::Zero();
		}

		J[t+1] = sys->dynfunc(J[t], U[t], VectorQ::Zero(), true);
		U[t] = (J[t+1] - J[t])/double(DT);
	}
}

void PR2EihBrett::bsp(StdVectorJ& J, StdVectorU& U, const MatrixJ& j_sigma0,
		const std::vector<Gaussian3d>& obj_gaussians, const std::vector<geometry3d::Triangle>& obstacles) {
	ROS_INFO_STREAM("Number of gaussians: " << obj_gaussians.size());
	ROS_INFO_STREAM("Number of obstacles: " << obstacles.size());

	// plan
	pr2_eih_bsp.collocation(J, U, j_sigma0, obj_gaussians, obstacles, *sys, true);

	// reintegrate, just in case
	for(int t=0; t < T-1; ++t) {
		J[t+1] = sys->dynfunc(J[t], U[t], VectorQ::Zero(), true);
	}
}

void PR2EihBrett::display_trajectory(const StdVectorJ& J) {
	geometry_msgs::PoseArray pose_array;
	pose_array.poses.resize(T);
	for(int t=0; t < T; ++t) {
		Matrix4d pose = arm_sim->fk(J[t]);
		pose_array.poses[t].position.x = pose(0,3);
		pose_array.poses[t].position.y = pose(1,3);
		pose_array.poses[t].position.z = pose(2,3);

		Quaterniond quat(pose.block<3,3>(0,0));
		pose_array.poses[t].orientation.w = quat.w();
		pose_array.poses[t].orientation.x = quat.x();
		pose_array.poses[t].orientation.y = quat.y();
		pose_array.poses[t].orientation.z = quat.z();
	}

	pose_array.header.frame_id = "/base_link";
	pose_array.header.stamp = ros::Time::now();

	display_trajectory_pub.publish(pose_array);
}

void PR2EihBrett::execute_controls(const StdVectorU& U) {
	VectorJ current_joints = arm->get_joints();
	VectorJ next_joints;

	double travel_distance = 0;
	std::vector<VectorJ> joint_traj;
	for(int t=0; t < U.size(); ++t) {
		next_joints = sys->dynfunc(current_joints, U[t], VectorQ::Zero(), true);
		double dist = (arm_sim->fk(next_joints) - arm_sim->fk(current_joints)).block<3,1>(0,3).norm();
		ROS_INFO_STREAM("dist: " << dist);
		if ((travel_distance + dist < max_travel_distance) || (t == 0)) {
			travel_distance += dist;
			joint_traj.push_back(next_joints);
			current_joints = next_joints;
		} else {
			break;
		}
		ROS_INFO_STREAM("travel_distance: " << travel_distance);
	}

	arm->execute_joint_trajectory(joint_traj);

//	VectorJ current_joints = arm->get_joints();
//
//	std::vector<VectorJ> joint_traj(U.size());
//	for(int t=0; t < U.size(); ++t) {
//		current_joints = sys->dynfunc(current_joints, U[t], VectorQ::Zero(), true);
//		joint_traj[t] = current_joints;
//	}
//
//	arm->execute_joint_trajectory(joint_traj);
}


void PR2EihBrett::_occluded_region_array_callback(const pcl_utils::OccludedRegionArrayConstPtr& msg) {
	current_occluded_region_array = *msg;
	received_occluded_region_array = true;
}

void parse_args(int argc, char* argv[], int& max_occluded_regions, Vector3d& init_traj_control, double& max_travel_distance) {
	std::vector<double> init_traj_control_vec;
	// Declare the supported options.
	po::options_description desc("Allowed options");
	desc.add_options()
		("help", "produce help message")
		("m", po::value<int>(&max_occluded_regions)->default_value(1), "Maximum number of occluded regions during BSP")
		("init", po::value<std::vector<double> >(&init_traj_control_vec)->multitoken(), "Delta position (x,y,z) by which to initialize BSP trajectory")
		("d", po::value<double>(&max_travel_distance)->default_value(1.0), "Maximum distance traveled when executing BSP controls")
		;
	try {
		po::variables_map vm;
		po::store(po::parse_command_line(argc, argv, desc, po::command_line_style::unix_style ^ po::command_line_style::allow_short), vm);
		po::notify(vm);

		if (vm.count("init")) {
			init_traj_control = Vector3d(init_traj_control_vec.data());
		} else {
			init_traj_control << .05, 0, .05;
		}

		if (vm.count("help")) {
			std::cout << desc << "\n";
			exit(0);
		}
	} catch (std::exception &e) {
		std::cerr << "error: " << e.what() << "\n";
		exit(0);
	}

	ROS_INFO_STREAM("Max occluded regions: " << max_occluded_regions);
	ROS_INFO_STREAM("Init traj control: " << init_traj_control.transpose());
	ROS_INFO_STREAM("Max travel distance: " << max_travel_distance);
}

int main(int argc, char* argv[]) {
	// initialize ros node
	ros::init(argc, argv, "pr2_eih_brett");
	log4cxx::LoggerPtr my_logger = log4cxx::Logger::getLogger(ROSCONSOLE_DEFAULT_NAME);
	my_logger->setLevel(ros::console::g_level_lookup[ros::console::levels::Info]);
	ROS_INFO("Starting ROS node");

	int max_occluded_regions;
	Vector3d init_traj_control;
	double max_travel_distance;
	parse_args(argc, argv, max_occluded_regions, init_traj_control, max_travel_distance);

	PR2EihBrett brett_bsp(max_occluded_regions, init_traj_control, max_travel_distance);

	std::vector<Gaussian3d> obj_gaussians;
	std::vector<geometry3d::Triangle> obstacles;

	MatrixJ j_sigma0 = (M_PI/4)*MatrixJ::Identity(); // TODO: never update in MPC
	StdVectorJ J(T);
	StdVectorU U(T-1);

	while(!ros::isShuttingDown()) {
		ROS_INFO("Getting occluded regions, press enter");
		std::cin.ignore();
		brett_bsp.get_occluded_regions(obj_gaussians, obstacles);
		brett_bsp.display_gaussian_means(obj_gaussians);

		ROS_INFO("Initializing trajectory, press enter");
		std::cin.ignore();
		brett_bsp.initialize_trajectory(J, U, obj_gaussians);
		brett_bsp.display_trajectory(J);

		ROS_INFO("Optimizing trajectory with bsp, press enter");
		std::cin.ignore();
		brett_bsp.bsp(J, U, j_sigma0, obj_gaussians, obstacles);

		ROS_INFO("Displaying bsp trajectory, press enter");
		std::cin.ignore();
		brett_bsp.display_trajectory(J);

		ROS_INFO("Executing control, press enter");
		std::cin.ignore();
		brett_bsp.execute_controls(U);
	}

}
