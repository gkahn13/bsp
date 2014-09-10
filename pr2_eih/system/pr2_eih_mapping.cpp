#include "pr2_eih_mapping.h"

/**
 *
 * Helper classes
 *
 */

class OccludedRegionCompare {
public:
	OccludedRegionCompare() { }

	// assumes a and b are in the camera frame
	bool operator()(const pcl_utils::OccludedRegion& a, const pcl_utils::OccludedRegion& b) const {
//		Vector3d a_gaussian_mean(a.gaussian.mean.x, a.gaussian.mean.y, a.gaussian.mean.z);
//		Vector3d b_gaussian_mean(b.gaussian.mean.x, b.gaussian.mean.y, b.gaussian.mean.z);
//		return (a_gaussian_mean.norm() < b_gaussian_mean.norm());
		Matrix3d a_cov(a.gaussian.covariance.data());
		Matrix3d b_cov(b.gaussian.covariance.data());
		return (a_cov.trace() > b_cov.trace());
	}
};

/**
 *
 * Constructors
 *
 */

PR2EihMapping::PR2EihMapping(int max_occluded_regions, double max_travel_distance) :
		max_occluded_regions(max_occluded_regions), max_travel_distance(max_travel_distance) {
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
	ros::Duration(2.0).sleep();
	ros::spinOnce();

	occluded_region_array_sub =
			nh_ptr->subscribe("/kinfu/occluded_regions", 1, &PR2EihMapping::_occluded_region_array_callback, this);
	grasp_joint_traj_sub =
			nh_ptr->subscribe("/check_handle_grasps/grasp_joint_trajectory", 1,
					&PR2EihMapping::_grasp_joint_traj_callback, this);
	return_grasp_joint_traj_sub =
				nh_ptr->subscribe("/check_handle_grasps/return_grasp_trajectory", 1,
						&PR2EihMapping::_return_grasp_joint_traj_callback, this);
	received_occluded_region_array = false;
	display_gaussian_pub = nh_ptr->advertise<visualization_msgs::MarkerArray>("/bsp_object_means", 1);
	display_trajectory_pub = nh_ptr->advertise<geometry_msgs::PoseArray>("/bsp_trajectory", 1);
	get_occluded_pub = nh_ptr->advertise<std_msgs::Empty>("/get_occluded", 1);
	reset_kinfu_pub = nh_ptr->advertise<std_msgs::Empty>("/reset_kinfu", 1);
	logger_pub = nh_ptr->advertise<std_msgs::String>("/experiment_log", 1);
	home_pose_pub = nh_ptr->advertise<geometry_msgs::PoseStamped>("/bsp/home_pose", 1);

	sim->update();
	home_joints = arm->get_joints();

	reset_kinfu_pub.publish(std_msgs::Empty());
	ros::spinOnce();
	ros::Duration(1.0).sleep();
}

/**
 *
 * Public methods
 *
 */

void PR2EihMapping::get_occluded_regions(std::vector<Gaussian3d>& obj_gaussians,
		std::vector<geometry3d::Triangle>& obstacles) {
	publish_to_logger("start get_occluded_regions");
	publish_home_pose();

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
			obstacles.push_back(geometry3d::Triangle(cam_rot*first_obstacle_point+cam_trans,
					cam_rot*second_obstacle_point+cam_trans, cam_rot*third_obstacle_point+cam_trans));
		}

		// add gaussian
		const pcl_utils::Gaussian& gaussian = occluded_region.gaussian;
		Vector3d mean(gaussian.mean.x, gaussian.mean.y, gaussian.mean.z);
		Matrix3d cov(gaussian.covariance.data());
		obj_gaussians.push_back(Gaussian3d(cam_rot*mean+cam_trans, cam_rot*cov*cam_rot.transpose()));

		ROS_INFO_STREAM("Gaussian " << i << "\nMean: " << obj_gaussians.back().mean.transpose()
				<< "\nCovariance:\n" << obj_gaussians.back().cov << "\n");
	}

	publish_to_logger("end get_occluded_regions");
}

void PR2EihMapping::display_gaussian_means(const std::vector<Gaussian3d>& obj_gaussians) {
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

void PR2EihMapping::display_trajectory(const StdVectorJ& J) {
	geometry_msgs::PoseArray pose_array;
	pose_array.poses.resize(J.size());
	for(int t=0; t < J.size(); ++t) {
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

void PR2EihMapping::execute_controls(const StdVectorU& U) {
	publish_to_logger("start execute_bsp");

	VectorJ current_joints = arm->get_joints();
	VectorJ next_joints;

	double travel_distance = 0;
	StdVectorJ joint_traj;
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

	publish_to_logger("end execute_bsp");
}

bool PR2EihMapping::is_valid_grasp_trajectory(StdVectorJ& grasp_joint_traj, StdVectorJ& return_grasp_joint_traj) {
	bool received_grasp_traj = (grasp_joint_traj.size() > 0);
	grasp_joint_traj = this->grasp_joint_traj;
	return_grasp_joint_traj = this->return_grasp_joint_traj;
	return received_grasp_traj;
}

void PR2EihMapping::execute_grasp_trajectory(const StdVectorJ& grasp_joint_traj,
		const StdVectorJ& return_grasp_joint_traj) {
	publish_to_logger("start execute_grasp_trajectory");
	double speed = 0.06;

	ROS_INFO("Opening gripper and executing grasp trajectory");
	arm->open_gripper(0, false);
	arm->execute_joint_trajectory(grasp_joint_traj, speed);

	ROS_INFO("Closing the gripper");
	arm->close_gripper(0, true, 2.0);

	if (return_grasp_joint_traj.size() > 0) {
		arm->execute_joint_trajectory(return_grasp_joint_traj, speed);
	} else {
		ROS_INFO("Return grasp traj length is 0, returning via hard-coded behavior");

		ROS_INFO("Moving up");
		arm->go_to_pose(highest_pose_above(arm->get_pose()), speed);

//		ROS_INFO("Moving back along higher trajectory");
//		int traj_length = grasp_traj.size();
//		std::vector<Matrix4d> up_grasp_traj(traj_length);
//		for(int i=0; i < traj_length; ++i) {
//			up_grasp_traj[i] = highest_pose_above(grasp_traj[(traj_length-1)-i]);
//		}
//		arm->execute_pose_trajectory(up_grasp_traj);

		ROS_INFO("Moving to above home pose");
		arm->go_to_pose(highest_pose_above(arm->fk(home_joints)), speed);

		ROS_INFO("Moving to home pose and dropping off");
		arm->go_to_joints(home_joints, speed);
	}

	arm->open_gripper(0, false);
	arm->go_to_joints(home_joints);

	reset_kinfu_pub.publish(std_msgs::Empty());
	ros::spinOnce();
	ros::Duration(1.0).sleep(); // give kinfu some time
	this->grasp_joint_traj.clear();
	this->return_grasp_joint_traj.clear();


	publish_to_logger("end execute_grasp_trajectory");
}

/**
 *
 * Private methods
 *
 */

Matrix4d PR2EihMapping::highest_pose_above(const Matrix4d& pose) {
	Vector3d step(0,0,0.01);
	Matrix4d highest_above = pose;
	VectorJ joints;
	while (arm_sim->ik(highest_above, joints)) {
		highest_above.block<3,1>(0,3) += step;
	}
	highest_above.block<3,1>(0,3) -= step;
	return highest_above;
}

void PR2EihMapping::publish_to_logger(std::string str) {
	std_msgs::String msg;
	msg.data = str.c_str();
	logger_pub.publish(msg);
}

void PR2EihMapping::publish_home_pose() {
	// publish home pose for check_handle_grasps
	Matrix4d home_pose = arm->fk(home_joints);
	Affine3d home_pose_affine = Translation3d(home_pose.block<3,1>(0,3)) * AngleAxisd(home_pose.block<3,3>(0,0));
	geometry_msgs::PoseStamped home_pose_msg;
	tf::Pose tf_pose;
	tf::poseEigenToTF(home_pose_affine, tf_pose);
	tf::poseTFToMsg(tf_pose, home_pose_msg.pose);
	home_pose_msg.header.frame_id = "base_link";
	home_pose_msg.header.stamp = ros::Time::now();
	home_pose_pub.publish(home_pose_msg);
}

/**
 *
 * Callbacks
 *
 */

void PR2EihMapping::_occluded_region_array_callback(const pcl_utils::OccludedRegionArrayConstPtr& msg) {
	assert(msg->header.frame_id == "camera_rgb_optical_frame" || msg->header.frame_id == "/camera_rgb_optical_frame");
	current_occluded_region_array = *msg;
	received_occluded_region_array = true;
}

void PR2EihMapping::_grasp_joint_traj_callback(const trajectory_msgs::JointTrajectoryConstPtr& msg) {
	grasp_joint_traj.clear();
	for(const trajectory_msgs::JointTrajectoryPoint& jt_pt : msg->points) {
		VectorJ joints(jt_pt.positions.data());
		grasp_joint_traj.push_back(joints);
	}
}

void PR2EihMapping::_return_grasp_joint_traj_callback(const trajectory_msgs::JointTrajectoryConstPtr& msg) {
	return_grasp_joint_traj.clear();
	for(const trajectory_msgs::JointTrajectoryPoint& jt_pt : msg->points) {
		VectorJ joints(jt_pt.positions.data());
		return_grasp_joint_traj.push_back(joints);
	}
}
