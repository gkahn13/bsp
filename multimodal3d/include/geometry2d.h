#ifndef __GEOMETRY2D_H__
#define __GEOMETRY2D_H__

#include <Eigen/Eigen>
using namespace Eigen;

#include <assert.h>

#define epsilon (1e-5)

inline double triangle_area(const Vector2d& a, const Vector2d& b, const Vector2d& c) {
	return fabs((c(0)*(a(1) - b(1)) + a(0)*(b(1) - c(1)) + b(0)*(c(1) - a(1))) / 2.0);
}

class Triangle2d {
public:
	Vector2d a, b, c;

	Triangle2d(const Vector2d& a_pt, const Vector2d& b_pt, const Vector2d& c_pt) : a(a_pt), b(b_pt), c(c_pt) { };

	Vector2d closest_point_to(const Vector2d& x);
	bool is_inside(const Vector2d& x);
	double area();
};


#endif
