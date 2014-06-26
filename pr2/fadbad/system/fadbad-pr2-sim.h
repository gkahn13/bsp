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

class FadbadCamera {
public:
	FadbadCamera(rave::SensorBasePtr s, double mr);

	std::vector<std::vector<FadbadBeam3d> > get_beams(const StdVector3b& env_points);
	std::vector<FadbadTriangle3d> get_border(const std::vector<std::vector<FadbadBeam3d> >& beams, bool with_side_border=true);

	bool is_inside(const Vector3b& p, std::vector<std::vector<FadbadBeam3d> >& beams);
	bdouble signed_distance(const Vector3b& p, std::vector<std::vector<FadbadBeam3d> >& beams, std::vector<FadbadTriangle3d>& border);

	inline Vector3b get_position() { return rave_utils::rave_to_eigen(sensor->GetTransform().trans).cast<bdouble>(); }
	inline rave::Transform get_pose() { return sensor->GetTransform(); }

private:
	rave::SensorBasePtr sensor;

	bdouble max_range;

	FadbadDepthMap* depth_map;

	MatrixDynb get_directions(const int h, const int w, const bdouble h_meters, const bdouble w_meters);
};

class FadbadPixelBucket {
public:
	FadbadPixelBucket(const Vector2b& pc) : pixel_center(pc) { };

	inline void add_point(const Vector2b& pixel, const Vector3b& point) { pixels.push_back(pixel); points.push_back(point); }
	inline bool is_empty() { return pixels.size() == 0; }
	inline Vector3b average_point() {

		VectorDynb weights(pixels.size());
		for(int i=0; i < pixels.size(); ++i) {
			weights(i) = (pixels[i] - pixel_center).norm();
		}
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
	FadbadDepthMap(rave::SensorBasePtr s, const Matrix3b& P_mat, bdouble mr);

	void add_point(const Vector3b& point);
	Matrix<bdouble,H_SUB,W_SUB> get_z_buffer();

	void clear();

private:
	rave::SensorBasePtr sensor;
	Matrix3b P;
	bdouble max_range;

	std::vector<std::vector<FadbadPixelBucket*> > pixel_buckets;

	int num_neighbors_empty(int i, int j);
	Vector3b average_of_neighbors(int i, int j);


};

#endif
