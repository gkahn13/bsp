#ifndef _FADBAD_PR2_SIM_H__
#define _FADBAD_PR2_SIM_H__

#include <Eigen/Eigen>
using namespace Eigen;

#include <openrave-core.h>
namespace rave = OpenRAVE;

#include "../geometry/fadbad-geometry3d.h"
#include "../fadbad-utils.h"
//#include "../utils/pr2-utils.h"
//#include "../../util/Timer.h"

#include "../../system/pr2-sim.h"
#include "../../utils/rave-utils.h"

typedef Matrix<bdouble,4,1> Vector4b;
typedef Matrix<bdouble,4,4> Matrix4b;

// forward declarations
class FadbadCamera;
class FadbadPixelBucket;
class FadbadDepthMap;

class FadbadArm {
public:
	FadbadArm(rave::RobotBasePtr robot, Arm::ArmType arm_type);

	Matrix4b get_pose(const Matrix<bdouble,ARM_DIM,1>& j);

private:
	Matrix4b origin;
	std::vector<Vector3b> arm_joint_axes, arm_link_trans;
};

class FadbadCamera {
public:
	FadbadCamera(FadbadArm* a, const Matrix4d& g_t_to_s);

	std::vector<std::vector<FadbadBeam3d> > get_beams(const Matrix<bdouble,ARM_DIM,1>& j, const StdVector3b& pcl);
	std::vector<FadbadTriangle3d> get_border(const std::vector<std::vector<FadbadBeam3d> >& beams, bool with_side_border=true);

	bool is_inside(const Vector3b& p, std::vector<std::vector<FadbadBeam3d> >& beams);
	bdouble signed_distance(const Vector3b& p, std::vector<std::vector<FadbadBeam3d> >& beams, std::vector<FadbadTriangle3d>& border);

	inline Matrix4b get_pose(const Matrix<bdouble,ARM_DIM,1>& j) { return arm->get_pose(j)*gripper_tool_to_sensor; }
	inline Vector3b get_position(const Matrix<bdouble,ARM_DIM,1>& j) { return get_pose(j).block<3,1>(0,3); }

private:
	FadbadArm* arm;
	Matrix4b gripper_tool_to_sensor;

	FadbadDepthMap* depth_map;

	Matrix<bdouble,N_SUB,3> get_directions(const Matrix<bdouble,ARM_DIM,1>& j);
};

class FadbadPixelBucket {
public:
	FadbadPixelBucket(const Vector2b& pc) : pixel_center(pc) { };

	inline void add_point(const Vector2b& pixel, const Vector3b& point) { pixels.push_back(pixel); points.push_back(point); }
	inline bool is_empty() { return pixels.size() == 0; }
	inline Vector3b average_point() {
		std::cout << "pixels size: " << pixels.size() << "\n";
		assert(pixels.size() > 0);

		VectorDynb weights(pixels.size());
		for(int i=0; i < pixels.size(); ++i) {
			weights(i) = (pixels[i] - pixel_center).norm();
		}
		std::cout << "weights sum: " << weights.sum().x() << "\n";
		weights /= weights.sum();

		Vector3b avg_point = Vector3b::Zero();
		for(int i=0; i < points.size(); ++i) {
			for(int j=0; j < 3; ++j) { avg_point(j) += weights(i)*points[i](j); }
		}

		return avg_point;
	}
	inline void clear() { pixels.clear(); points.clear(); }

private:
	Vector2b pixel_center;
	StdVector2b pixels;
	StdVector3b points;
};

/**
 * FadbadDepthMap size is H_SUB x W_SUB
 */
class FadbadDepthMap {
public:
	FadbadDepthMap(const Matrix3b& P_mat);

	void add_point(const Vector3b& point, const Matrix4b& cam_pose);
	Matrix<bdouble,H_SUB,W_SUB> get_z_buffer(const Vector3b& cam_pos);

	void clear();

private:
	Matrix3b P;

	std::vector<std::vector<FadbadPixelBucket*> > pixel_buckets;

	int num_neighbors_empty(int i, int j);
	Vector3b average_of_neighbors(int i, int j);


};

#endif
