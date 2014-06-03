#pragma once

#include <map>
#include <ncurses.h>

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
class Head;
class Sensor;
class CameraSensor;

class PR2 {
public:
	Arm *larm, *rarm;
	Head *head;
	CameraSensor *camera;

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
	rave::Transform get_pose();
	void get_limits(mat &lower, mat &upper);

	void set_posture(Posture posture);
	void set_joint_values(const mat &joint_values);
	void set_pose(const rave::Transform &pose, std::string ref_frame="world");

	void teleop();

private:
	ArmType arm_type;
	std::string manip_name;
	rave::RobotBasePtr robot;
	rave::RobotBase::ManipulatorPtr manip;
	std::vector<int> arm_indices;
	int num_joints;
};

class Head {
public:
	Head(rave::RobotBasePtr robot);

	mat get_joint_values();
	void get_limits(mat &lower, mat &upper);

	void set_joint_values(const mat &joint_values);
	void look_at(const rave::Transform &pose,
			const std::string reference_frame="world", const std::string camera_frame="wide_stereo_link");

private:
	rave::RobotBasePtr robot;
	std::vector<int> head_indices;
	int num_joints;
};

class Sensor {
public:
	Sensor(rave::SensorBasePtr sensor);

	void power_on();
	void power_off();
	void render_on();
	void render_off();

	rave::SensorBase::SensorDataPtr get_data();

private:
	rave::SensorBasePtr sensor;
	bool is_powered, is_rendering;
	rave::SensorBase::SensorType type;
};

class CameraSensor : public Sensor {
public:
	CameraSensor(rave::SensorBasePtr sensor);

	cube get_image();
	std::vector<std::vector<mat> > get_pixels_and_colors(const std::vector<mat> &points);
	mat get_pixel_from_point(const mat &point);
	bool is_in_fov(const mat& point);

private:
	rave::SensorBasePtr sensor;
	int height, width;
	mat P;
};
