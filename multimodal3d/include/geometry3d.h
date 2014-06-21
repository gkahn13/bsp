#ifndef __GEOMETRY3D_H__
#define __GEOMETRY3D_H__

#include "geometry2d.h"

#include <math.h>

#include <Eigen/Eigen>
using namespace Eigen;

#include <assert.h>

class Beam3d {
public:
	Vector3d base, a, b, c, d;

	Beam3d(const Vector3d& base_tmp, const Vector3d& a_tmp, const Vector3d& b_tmp, const Vector3d& c_tmp, const Vector3d& d_tmp) :
		base(base_tmp), a(a_tmp), b(b_tmp), c(c_tmp), d(d_tmp) { };

	bool is_inside(const Vector3d& p);
};

class Triangle3d {
	Vector3d a, b, c;

	Triangle3d(const Vector3d& a_tmp, const Vector3d& b_tmp, const Vector3d& c_tmp) :
		a(a_tmp), b(b_tmp), c(c_tmp) { };

	Vector3d closest_point_to(const Vector3d& p);
	inline double distance_to(const Vector3d& p) { return (closest_point_to(p) - p).norm(); }
	double area();

private:
	Matrix3d rotation_to_align_with(const Vector3d& target);
};

#endif
