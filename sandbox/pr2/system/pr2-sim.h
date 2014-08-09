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
#include "../utils/pr2-utils.h"
#include "../../util/Timer.h"

#define ARM_DIM 7
#define HEAD_DIM 2

// Camera constants for actual and subsampled

#define WIDTH_FULL	   640 // 256
#define HEIGHT_FULL    480 // 192

#define W_SUB 64 // 64
#define H_SUB 48 // 48
#define N_SUB (W_SUB*H_SUB)

namespace intrinsics {
const double FOCAL_LENGTH = .01;
const double MAX_RANGE = 0.75; // 0.75
const double MIN_RANGE = 0.2;

const double fx  = WIDTH_FULL*2.0;
const double fy  = HEIGHT_FULL*2.0;
const double cx  = double(WIDTH_FULL)/2.0 + 0.5;
const double cy  = double(HEIGHT_FULL)/2.0 + 0.5;
const double HEIGHT_M = FOCAL_LENGTH*(HEIGHT_FULL/fy);
const double WIDTH_M = FOCAL_LENGTH*(WIDTH_FULL/fx);

const double fx_sub  = W_SUB*2.0;
const double fy_sub  = H_SUB*2.0;
const double cx_sub  = double(W_SUB)/2.0 + 0.5;
const double cy_sub  = double(H_SUB)/2.0 + 0.5;
const double H_SUB_M = FOCAL_LENGTH*(H_SUB/fy_sub);
const double W_SUB_M = FOCAL_LENGTH*(W_SUB/fx_sub);
}

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
	rave::RobotBasePtr get_robot() { return robot; }

private:
	void init(std::string env_file, std::string robot_name, bool view);

	rave::EnvironmentBasePtr env;
	rave::ViewerBasePtr viewer;
	boost::shared_ptr<boost::thread> viewer_thread;

	rave::RobotBasePtr robot;
};

class Arm {
	const double step = 1e-5;

public:
	enum ArmType { left, right };
	enum Posture { untucked, tucked, up, side, mantis };

	Arm(rave::RobotBasePtr robot, ArmType arm_type);

	Matrix<double,ARM_DIM,1> get_joint_values();
	void get_limits(Matrix<double,ARM_DIM,1>& lower, Matrix<double,ARM_DIM,1>& upper);
	Matrix4d get_pose(const Matrix<double,ARM_DIM,1>& j);
	Vector3d get_position(const Matrix<double,ARM_DIM,1>& j) { return get_pose(j).block<3,1>(0,3); }

	Matrix<double,3,ARM_DIM> get_position_jacobian(const Matrix<double,ARM_DIM,1>& j);
	bool ik(const Vector3d& pos, Matrix<double,ARM_DIM,1>& j);

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

	rave::Transform origin;
	std::vector<rave::Vector> arm_joint_axes, arm_link_trans;
};

//class Head {
//public:
//	Head(rave::RobotBasePtr robot);
//
//	Matrix<double,HEAD_DIM,1> get_joint_values();
//	void get_limits(Matrix<double,HEAD_DIM,1>& lower, Matrix<double,HEAD_DIM,1>& upper);
//	rave::Transform get_pose();
//
//	void set_joint_values(Matrix<double,HEAD_DIM,1>& j);
//	void look_at(const rave::Transform &pose, const std::string ref_frame="world");
//
//	void teleop();
//
//private:
//	rave::RobotBasePtr robot;
//	std::vector<int> joint_indices;
//	int num_joints;
//	Matrix<double,HEAD_DIM,1> lower, upper;
//
//	rave::KinBody::LinkPtr pose_link;
//};


class Camera {
public:
	Camera(rave::RobotBasePtr r, std::string camera_name, Arm* a);

	// call once before collocation
	StdVector3d get_pc(const Matrix<double,ARM_DIM,1>& j);
	Matrix<double,HEIGHT_FULL,WIDTH_FULL> get_full_zbuffer(const Matrix<double,ARM_DIM,1>& j, const StdVector3d& obstacles);
	bool is_in_fov_full(const Vector3d& point, const Matrix<double,HEIGHT_FULL,WIDTH_FULL>& zbuffer, const Matrix4d& cam_pose);

	Matrix<double,H_SUB,W_SUB> get_zbuffer(const Matrix<double,ARM_DIM,1>& j, const StdVector3d& obstacles);

	Vector2i get_pixel_from_point(const Vector3d& point, const Matrix4d& cam_pose, const Matrix3d& P);
	Vector3d get_point_from_pixel_and_dist(const Vector2i& pixel, const double dist, const Matrix4d& cam_pose);
	bool is_in_fov(const Vector3d& point, const Matrix<double,H_SUB,W_SUB>& zbuffer, const Matrix4d& cam_pose);

	inline Matrix4d get_pose(const Matrix<double,ARM_DIM,1>& j) { return arm->get_pose(j)*gripper_tool_to_sensor; }
	inline Vector3d get_position(const Matrix<double,ARM_DIM,1>& j) { return get_pose(j).block<3,1>(0,3); }

	inline rave::SensorBasePtr get_sensor() { return sensor; }
	inline Matrix4d get_gripper_tool_to_sensor() { return gripper_tool_to_sensor; }

private:
	rave::RobotBasePtr robot;
	rave::SensorBasePtr sensor;
	Arm* arm;

	Matrix3d KK, KK_SUB;

	Matrix4d gripper_tool_to_sensor;

	Beam3d* fov;

	MatrixXd get_directions(const Matrix<double,ARM_DIM,1>& j, const int h, const int w, const double h_meters, const double w_meters);
};


#endif
