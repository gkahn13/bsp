#pragma once

#include <armadillo>
using namespace arma;

#include <boost/thread/thread.hpp>
#include <boost/bind.hpp>
#include <boost/filesystem.hpp>

#include <openrave-core.h>
namespace rave = OpenRAVE;

#include "rave_utils.h"

// forward declarations
class Arm;

class PR2 {
public:
	Arm *larm, *rarm;

	PR2();
	PR2(std::string env_file, std::string robot_name, bool view=true);
	~PR2();

private:
	void init(std::string env_file, std::string robot_name, bool view);

	rave::EnvironmentBasePtr env;
	rave::ViewerBasePtr viewer;
	boost::shared_ptr<boost::thread> viewer_thread;

	rave::RobotBasePtr robot;
};

class Arm {
public:
	enum ArmType { left, right };
	enum Posture { untucked, tucked, up, side, mantis };

	Arm(rave::RobotBasePtr robot, ArmType arm_type);

	mat get_joint_values();
	mat get_pose();
	void get_lower_limits(mat &lower, mat &upper);

	void set_posture(Posture posture);
	void set_joint_values(mat &joint_values);
	void set_pose(mat &pose, std::string ref_frame="world");

private:
	void init(rave::RobotBasePtr robot, ArmType arm_type);

	ArmType arm_type;
	rave::RobotBasePtr robot;
	rave::RobotBase::ManipulatorPtr manip;
	std::vector<int> arm_indices;
	int num_joints;
};
