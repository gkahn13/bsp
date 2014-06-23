#ifndef _PR2_SIM_H__
#define _PR2_SIM_H__

#include <map>

#include <Eigen/Eigen>
using namespace Eigen;

#include <boost/thread/thread.hpp>
#include <boost/bind.hpp>
#include <boost/filesystem.hpp>

#include <openrave-core.h>
namespace rave = OpenRAVE;

#include "../geometry/geometry3d.h"
#include "../utils/rave-utils.h"
#include "../utils/utils.h"
#include "../../util/Timer.h"

#define ARM_DIM 7
#define HEAD_DIM 2

#define H_SUB 48
#define W_SUB 64
#define N_SUB (H_SUB*W_SUB)

// forward declarations
class Arm;
class Head;
class Camera;

class PR2 {
public:
	Arm *larm, *rarm;
	Head *head;
	Camera *hcam, *lcam, *rcam;

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

class Arm {
public:
	enum ArmType { left, right };
	enum Posture { untucked, tucked, up, side, mantis };

	Arm(rave::RobotBasePtr robot, ArmType arm_type);

	Matrix<double,ARM_DIM,1> get_joint_values();
	void get_limits(Matrix<double,ARM_DIM,1>& lower, Matrix<double,ARM_DIM,1>& upper);
	rave::Transform get_pose();

	void set_joint_values(const Matrix<double,ARM_DIM,1>& j);
	void set_pose(const rave::Transform &pose, std::string ref_frame="world");
	void set_posture(Posture posture);

	void teleop();

private:
	ArmType arm_type;
	std::string manip_name;
	rave::RobotBase::ManipulatorPtr manip;

	rave::RobotBasePtr robot;
	std::vector<int> joint_indices;
	int num_joints;
	Matrix<double,ARM_DIM,1> lower, upper;
};

class Head {
public:
	Head(rave::RobotBasePtr robot);

	Matrix<double,HEAD_DIM,1> get_joint_values();
	void get_limits(Matrix<double,HEAD_DIM,1>& lower, Matrix<double,HEAD_DIM,1>& upper);
	rave::Transform get_pose();

	void set_joint_values(Matrix<double,HEAD_DIM,1>& j);
	void look_at(const rave::Transform &pose, const std::string ref_frame="world");

	void teleop();

private:
	rave::RobotBasePtr robot;
	std::vector<int> joint_indices;
	int num_joints;
	Matrix<double,HEAD_DIM,1> lower, upper;

	rave::KinBody::LinkPtr pose_link;
};


class Camera {
public:
	Camera(rave::RobotBasePtr r, std::string camera_name, double mr);

	Matrix<double,N_SUB,3> get_directions();
	std::vector<std::vector<Beam3d> > get_beams();
	std::vector<Triangle3d> get_border(const std::vector<std::vector<Beam3d> >& beams);

	bool is_inside(const Vector3d& p, std::vector<std::vector<Beam3d> >& beams);
	double signed_distance(const Vector3d& p, std::vector<std::vector<Beam3d> >& beams, std::vector<Triangle3d>& border);

	void plot_fov(std::vector<std::vector<Beam3d> >& beams);

	inline Vector3d get_position() { return rave_utils::rave_to_eigen(sensor->GetTransform().trans); }
	inline rave::Transform get_pose() { return sensor->GetTransform(); }

private:
	rave::RobotBasePtr robot;
	rave::SensorBasePtr sensor;

	int height, width;
	double f, F, max_range;
	double H, W;
};

#endif
