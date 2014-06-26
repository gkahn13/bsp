#ifndef __FADBAD_GEOMETRY3D_H__
#define __FADBAD_GEOMETRY3D_H__

#include "../fadbad-utils.h"

#include "fadbad-geometry2d.h"

#include <math.h>
#include <queue>

#include <Eigen/Eigen>
using namespace Eigen;

#include <assert.h>

typedef Matrix<bdouble,3,1> Vector3b;
typedef Matrix<bdouble,3,3> Matrix3b;

typedef std::vector<Vector3b> StdVector3d;

class FadbadBeam3d {
public:
	/**
	 * Order goes
	 * b ---- a
	 * |      |
	 * |      |
	 * c ---- d
	 */
	Vector3b base, a, b, c, d;

	FadbadBeam3d() { };
	FadbadBeam3d(const Vector3b& base_tmp, const Vector3b& a_tmp, const Vector3b& b_tmp, const Vector3b& c_tmp, const Vector3b& d_tmp) :
		base(base_tmp), a(a_tmp), b(b_tmp), c(c_tmp), d(d_tmp) { };

	bool is_inside(const Vector3b& p);
	bool is_crossed_by(const Vector3b& p0, const Vector3b& p1);
};

class FadbadTriangle3d {
public:
	Vector3b a, b, c;

	FadbadTriangle3d() { };
	FadbadTriangle3d(const Vector3b& a_tmp, const Vector3b& b_tmp, const Vector3b& c_tmp) :
		a(a_tmp), b(b_tmp), c(c_tmp) { };

	Vector3b closest_point_to(const Vector3b& p);
	inline bdouble distance_to(const Vector3b& p) { return (closest_point_to(p) - p).norm(); }

	bool intersection_with_segment(const Vector3b& p0, const Vector3b& p1, Vector3b& intersection);

	bdouble area();

private:
	Matrix3b rotation_to_align_with(const Vector3b& target);
};


#endif
