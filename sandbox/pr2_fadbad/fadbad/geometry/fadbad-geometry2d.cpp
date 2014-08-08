#include "fadbad-geometry2d.h"

/**
 * General functions
 */

inline Vector2b fadbad_closest_point_to_segment(const Vector2b& p0, const Vector2b& p1, const Vector2b& x) {
	// min_{0<=t<=1} ||t*(p1-p0) + p0 - x||_{2}^{2}
	Vector2b v = p1 - p0;
	Vector2b b = p0 - x;

	bdouble t = -v.dot(b) / v.squaredNorm();
	if ((0 <= t) && (t <= 1)) {
//		return (t*(p1 - p0) + p0);
		Vector2b closest(t*(p1(0) - p0(0)) + p0(0), t*(p1(1) - p0(1)) + p0(1));
		return closest;
	} else {
		if ((x - p0).norm() < (x - p1).norm()) {
			return p0;
		} else {
			return p1;
		}
	}
}

/**
 * FadbadTriangle2d public methods
 */

/**
 * \brief Assumes x is not inside the triangle
 */
Vector2b FadbadTriangle2d::closest_point_to(const Vector2b& x) {
	assert(!is_inside(x));

	Vector2b p0, p1, p2;
	p0 = fadbad_closest_point_to_segment(a, b, x);
	p1 = fadbad_closest_point_to_segment(b, c, x);
	p2 = fadbad_closest_point_to_segment(c, a, x);

	bdouble dist0, dist1, dist2;
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

bool FadbadTriangle2d::is_inside(const Vector2b& x) {
	bdouble total_area, area0, area1, area2;
	total_area = fadbad_triangle_area(a, b, c);
	area0 = fadbad_triangle_area(a, b, x);
	area1 = fadbad_triangle_area(b, c, x);
	area2 = fadbad_triangle_area(c, a, x);

	return (abs(total_area - (area0 + area1 + area2)) < bepsilon);
}

bdouble FadbadTriangle2d::area() {
	return fadbad_triangle_area(a, b, c);
}
