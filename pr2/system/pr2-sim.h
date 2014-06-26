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

#define FOCAL_LENGTH .01

#define WIDTH	   256 // 64
#define HEIGHT 	   192 // 48
const double fx  = WIDTH*2.0;
const double fy  = HEIGHT*2.0;
const double cx  = double(WIDTH)/2.0 + 0.5;
const double cy  = double(HEIGHT)/2.0 + 0.5;
const double HEIGHT_M = FOCAL_LENGTH*(HEIGHT/fy);
const double WIDTH_M = FOCAL_LENGTH*(WIDTH/fx);

#define W_SUB 64 // 64
#define H_SUB 48 // 48
const double fx_sub  = W_SUB*2.0;
const double fy_sub  = H_SUB*2.0;
const double cx_sub  = double(W_SUB)/2.0 + 0.5;
const double cy_sub  = double(H_SUB)/2.0 + 0.5;
const double H_SUB_M = FOCAL_LENGTH*(H_SUB/fy_sub);
const double W_SUB_M = FOCAL_LENGTH*(W_SUB/fx_sub);
#define N_SUB (W_SUB*H_SUB)

// forward declarations
class Arm;
class Head;
class Camera;
class PixelBucket;
class DepthMap;

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

	// call once before collocation
	void get_pcl();

	std::vector<std::vector<Beam3d> > get_beams();
	std::vector<Triangle3d> get_border(const std::vector<std::vector<Beam3d> >& beams, bool with_side_border=true);

	bool is_inside(const Vector3d& p, std::vector<std::vector<Beam3d> >& beams);
	double signed_distance(const Vector3d& p, std::vector<std::vector<Beam3d> >& beams, std::vector<Triangle3d>& border);

	void plot_fov(std::vector<std::vector<Beam3d> >& beams);
	void plot_pcl();

	inline Vector3d get_position() { return rave_utils::rave_to_eigen(sensor->GetTransform().trans); }
	inline rave::Transform get_pose() { return sensor->GetTransform(); }
	inline rave::SensorBasePtr get_sensor() { return sensor; }

private:
	rave::RobotBasePtr robot;
	rave::SensorBasePtr sensor;

//	int height, width;
//	double f, F, max_range;
//	double H, W;
	double max_range;

	Beam3d* fov = nullptr;

	StdVector3d env_points;
	DepthMap* depth_map;

	MatrixXd get_directions(const int h, const int w, const double h_meters, const double w_meters);
};

class PixelBucket {
public:
	PixelBucket(const Vector2d& pc) : pixel_center(pc) { };

	inline void add_point(const Vector2d& pixel, const Vector3d& point) { pixels.push_back(pixel); points.push_back(point); }
	inline bool is_empty() { return pixels.size() == 0; }
	inline Vector3d average_point() {
//		double dist_to_center = INFINITY;
//		Vector3d closest_point;
//		for(int i=0; i < pixels.size(); ++i) {
//			if ((pixels[i] - pixel_center).norm() < dist_to_center) {
//				dist_to_center = (pixels[i] - pixel_center).norm();
//				closest_point = points[i];
//			}
//		}
//		return closest_point;

		VectorXd weights(pixels.size());
		for(int i=0; i < pixels.size(); ++i) {
			weights(i) = (pixels[i] - pixel_center).norm();
		}
		weights /= weights.sum();

		Vector3d avg_point = Vector3d::Zero();
		for(int i=0; i < points.size(); ++i) {
			avg_point += weights(i)*points[i];
		}

		return avg_point;
	}
	inline void clear() { pixels.clear(); points.clear(); }

private:
	Vector2d pixel_center;
	StdVector2d pixels;
	StdVector3d points;
};

/**
 * DepthMap size is H_SUB x W_SUB
 */
class DepthMap {
public:
	DepthMap(rave::SensorBasePtr s, const Matrix3d& P_mat, double mr);

	void add_point(const Vector3d& point);
	Matrix<double,H_SUB,W_SUB> get_z_buffer();

	void clear();

private:
	rave::SensorBasePtr sensor;
	Matrix3d P;
	double max_range;

	std::vector<std::vector<PixelBucket*> > pixel_buckets;

	int num_neighbors_empty(int i, int j);
	Vector3d average_of_neighbors(int i, int j);


};

#endif
