#ifndef _PR2_SIM_H__
#define _PR2_SIM_H__

#include <map>

#include <armadillo>
using namespace arma;

#include <boost/thread/thread.hpp>
#include <boost/bind.hpp>
#include <boost/filesystem.hpp>

#include <openrave-core.h>
namespace rave = OpenRAVE;

#include "rave_utils.h"
#include "utils.h"
#include "../../../util/Timer.h"

// forward declarations
class Arm;
class Head;
class Sensor;
class KinectSensor;

class PR2 {
public:
	Arm *larm, *rarm;
	Head *head;
	KinectSensor *h_kinect, *l_kinect, *r_kinect;

	PR2(bool view=true);
	PR2(std::string env_file, std::string robot_name, bool view=true);
	~PR2();

	rave::EnvironmentBasePtr get_env() { return env; }

private:
	void init(std::string env_file, std::string robot_name, bool view);

	rave::EnvironmentBasePtr env;
	rave::ViewerBasePtr viewer;
	boost::shared_ptr<boost::thread> viewer_thread;

	rave::RobotBasePtr robot;
};

class Manipulator {
public:
	mat get_joint_values();
	void get_limits(mat &lower, mat &upper);
	virtual rave::Transform get_pose() = 0;

	void set_joint_values(mat j);
	virtual void set_pose(const rave::Transform &pose, std::string ref_frame="world") =0;

	// DO NOT HOLD DOWN KEYS
	virtual void teleop() = 0;

protected:
	rave::RobotBasePtr robot;
	std::vector<int> joint_indices;
	int num_joints;
	mat lower, upper;

	void init();
};

class Arm : public virtual Manipulator {
public:
	enum ArmType { left, right };
	enum Posture { untucked, tucked, up, side, mantis };

	Arm(rave::RobotBasePtr robot, ArmType arm_type);

	rave::Transform get_pose();

	void set_posture(Posture posture);
	void set_pose(const rave::Transform &pose, std::string ref_frame="world");

	void teleop();

private:
	ArmType arm_type;
	std::string manip_name;
	rave::RobotBase::ManipulatorPtr manip;
};

class Head : public virtual Manipulator {
public:
	Head(rave::RobotBasePtr robot);

	rave::Transform get_pose();

	// look_at
	void set_pose(const rave::Transform &pose, const std::string ref_frame="world");

	void teleop();

private:
	rave::KinBody::LinkPtr pose_link;
};

class Sensor {
public:
	Sensor(rave::SensorBasePtr sensor);

	void power_on();
	void power_off();
	void render_on();
	void render_off();

	rave::SensorBase::SensorDataPtr get_data(bool wait_for_new=true);
	rave::Transform get_pose();

protected:
	rave::SensorBasePtr sensor;
	bool is_powered, is_rendering;
	rave::SensorBase::SensorType type;
	double data_timeout;
};

class DepthSensor : public Sensor {
public:
	DepthSensor(rave::SensorBasePtr sensor);

	std::vector<mat> get_points(bool wait_for_new=true);

	bool is_in_range(const mat& point);

	double get_min_range() { return min_range; }
	double get_max_range() { return max_range; }
	double get_optimal_range() { return optimal_range; }

private:
	double min_range, max_range, optimal_range;
};

class CameraSensor : public Sensor {
public:
	CameraSensor(rave::SensorBasePtr sensor);

	cube get_image(bool wait_for_new=true);
	std::vector<std::vector<mat> > get_pixels_and_colors(const std::vector<mat> &points, bool wait_for_new=true);
	mat get_pixel_from_point(const mat &point, bool is_round=true);
	bool is_in_fov(const mat& point);

	int get_height() { return height; }
	int get_width() { return width; }

private:
	int height, width;
	mat P;
};

class ColoredPoint {
public:
	mat point, color;
	ColoredPoint(mat &p, mat& c) : point(p), color(c) { };

	rave::GraphHandlePtr display(rave::EnvironmentBasePtr env) { return rave_utils::plot_point(env, point, color);	}
};

class KinectSensor {
public:
	KinectSensor(rave::RobotBasePtr robot, std::string depth_sensor_name, std::string camera_sensor_name);
	~KinectSensor();

	void power_on();
	void power_off();
	void render_on();
	void render_off();

	std::vector<ColoredPoint*> get_point_cloud(bool wait_for_new=true);
	mat get_z_buffer(bool wait_for_new=true);

	cube get_image(bool wait_for_new=true) { return camera_sensor->get_image(wait_for_new); }
	mat get_pixel_from_point(const mat &point, bool is_round=true) { return camera_sensor->get_pixel_from_point(point, is_round); }
	bool is_visible(const mat &point) { return (camera_sensor->is_in_fov(point) && depth_sensor->is_in_range(point)); }
	rave::Transform get_pose() { return camera_sensor->get_pose(); }

	void display_point_cloud(const std::vector<ColoredPoint*> &colored_points, std::vector<rave::GraphHandlePtr> &handles);

	int get_height() { return camera_sensor->get_height(); }
	int get_width() { return camera_sensor->get_width(); }
	double get_min_range() { return depth_sensor->get_min_range(); }
	double get_max_range() { return depth_sensor->get_max_range(); }
	double get_optimal_range() { return depth_sensor->get_optimal_range(); }

	// TODO: temp!
	DepthSensor* get_depth_sensor() { return depth_sensor; }
	CameraSensor* get_camera_sensor() { return camera_sensor; }

private:
	rave::RobotBasePtr robot;
	DepthSensor* depth_sensor;
	CameraSensor* camera_sensor;
};

#endif
