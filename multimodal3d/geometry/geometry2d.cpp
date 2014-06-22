#include "geometry2d.h"

/**
 * General functions
 */

inline Vector2d closest_point_to_segment(const Vector2d& p0, const Vector2d& p1, const Vector2d& x) {
	// min_{0<=t<=1} ||t*(p1-p0) + p0 - x||_{2}^{2}
	Vector2d v = p1 - p0;
	Vector2d b = p0 - x;

	double t = -v.dot(b) / v.squaredNorm();
	if ((0 <= t) && (t <= 1)) {
		return (t*(p1 - p0) + p0);
	} else {
		if ((x - p0).norm() < (x - p1).norm()) {
			return p0;
		} else {
			return p1;
		}
	}
}

/**
 * Triangle2d public methods
 */

/**
 * \brief Assumes x is not inside the triangle
 */
Vector2d Triangle2d::closest_point_to(const Vector2d& x) {
	assert(!is_inside(x));

	Vector2d p0, p1, p2;
	p0 = closest_point_to_segment(a, b, x);
	p1 = closest_point_to_segment(b, c, x);
	p2 = closest_point_to_segment(c, a, x);

	double dist0, dist1, dist2;
	dist0 = (p0 - x).norm();
	dist1 = (p1 - x).norm();
	dist2 = (p2 - x).norm();

	if ((dist0 <= dist1) && (dist0 <= dist2)) {
		return p0;
	} else if ((dist1 <= dist0) && (dist1 <= dist2)) {
		return p1;
	} else {
		return p2;
	}
}

bool Triangle2d::is_inside(const Vector2d& x) {
	double total_area, area0, area1, area2;
	total_area = triangle_area(a, b, c);
	area0 = triangle_area(a, b, x);
	area1 = triangle_area(b, c, x);
	area2 = triangle_area(c, a, x);

	return (fabs(total_area - (area0 + area1 + area2)) < epsilon);
}

double Triangle2d::area() {
	return triangle_area(a, b, c);
}
