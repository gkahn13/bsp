#ifndef __PR2_EIH_MAPPING_H__
#define __PR2_EIH_MAPPING_H__

#include "pr2_eih_system.h"
#include "../../util/logging.h"

#include "pr2_utils/pr2/arm.h"
#include "pcl_utils/OccludedRegionArray.h"
#include "pcl_utils/BoundingBox.h"

class PR2EihMapping {
public:
	PR2EihMapping(int max_occluded_regions, double max_travel_distance);
	PR2EihMapping() : PR2EihMapping(1, 0.1) { }

	void reset_kinfu();

	void get_occluded_regions(std::vector<Gaussian3d>& obj_gaussians,
			std::vector<geometry3d::Triangle>& obstacles);
	void display_gaussian_means(const std::vector<Gaussian3d>& obj_gaussians);

	// add planning in subclasses!

	void display_trajectory(const StdVectorJ& J);
	void execute_controls(const StdVectorU& U);

	bool is_valid_grasp_trajectory(StdVectorJ& grasp_joint_traj, StdVectorJ& return_grasp_joint_traj);
	void execute_grasp_trajectory(const StdVectorJ& grasp_joint_traj, const StdVectorJ& return_grasp_joint_traj);

protected:
	// ros
	ros::NodeHandle *nh_ptr;
	ros::Subscriber occluded_region_array_sub, grasp_joint_traj_sub, return_grasp_joint_traj_sub, table_pose_sub;
	ros::Publisher display_trajectory_pub, get_occluded_pub, display_gaussian_pub, reset_kinfu_pub, head_camera_time_pub,
		logger_pub, home_pose_pub;

	bool received_occluded_region_array;
	pcl_utils::OccludedRegionArray current_occluded_region_array;

	int max_occluded_regions;
	double max_travel_distance;

	// pr2_utils
	pr2_sim::Simulator *sim;
	pr2_sim::Arm *arm_sim;
	pr2_sim::Camera *cam_sim;
	PR2EihSystem *sys;
	pr2::Arm *arm;

	VectorJ home_joints;
	StdVectorJ grasp_joint_traj, return_grasp_joint_traj;
	geometry3d::Halfspace table_halfspace;

	Matrix4d highest_pose_above(const Matrix4d& pose);
	void publish_to_logger(std::string str);
	void publish_home_pose();
	void _occluded_region_array_callback(const pcl_utils::OccludedRegionArrayConstPtr& msg);
	void _grasp_joint_traj_callback(const trajectory_msgs::JointTrajectoryConstPtr& msg);
	void _return_grasp_joint_traj_callback(const trajectory_msgs::JointTrajectoryConstPtr& msg);
	void _table_pose_callback(const pcl_utils::BoundingBoxConstPtr& msg);
};

#endif
