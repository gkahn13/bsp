#include "fadbad-geometry3d.h"

/**
 * General functions
 */

inline bool halfplane_contains_point(const Vector3b& origin, const Vector3b& normal, const Vector3b& x) {
	return (normal.dot(x - origin) >= -bepsilon);
}

/**
 * FadbadBeam3d public methods
 */

bool FadbadBeam3d::is_inside(const Vector3b& p) {
	bool inside = true;

	inside &= halfplane_contains_point(base, (a-base).cross(d-base), p);
	inside &= halfplane_contains_point(base, (b-base).cross(a-base), p);
	inside &= halfplane_contains_point(base, (c-base).cross(b-base), p);
	inside &= halfplane_contains_point(base, (d-base).cross(c-base), p);
	inside &= halfplane_contains_point(a, (b-a).cross(d-a), p);

	return inside;
}

bool FadbadBeam3d::is_crossed_by(const Vector3b& p0, const Vector3b& p1) {
	bool doesnt_cross = false;

	doesnt_cross |= (halfplane_contains_point(base, -(a-base).cross(d-base), p0) && halfplane_contains_point(base, -(a-base).cross(d-base), p1));
	doesnt_cross |= (halfplane_contains_point(base, -(b-base).cross(a-base), p0) && halfplane_contains_point(base, -(b-base).cross(a-base), p1));
	doesnt_cross |= (halfplane_contains_point(base, -(c-base).cross(b-base), p0) && halfplane_contains_point(base, -(c-base).cross(b-base), p1));
	doesnt_cross |= (halfplane_contains_point(base, -(d-base).cross(c-base), p0) && halfplane_contains_point(base, -(d-base).cross(c-base), p1));

	return (!doesnt_cross);
}


/**
 * FadbadTriangle3d public methods
 */


/**
 * Rotates all points to be aligned with z-axis
 * Then projects to 2d, combining result back into 3d
 */
Vector3b FadbadTriangle3d::closest_point_to(const Vector3b& p) {
	Matrix3b rotation = rotation_to_align_with({0,0,1});

	Vector3b a_rot, b_rot, c_rot, p_rot;
	a_rot = rotation*a;
	b_rot = rotation*b;
	c_rot = rotation*c;
	p_rot = rotation*p;

	bdouble tri_rot_z = a_rot(2);

	Vector2b p_rot_2d(p_rot(0), p_rot(1));

	FadbadTriangle2d tri_2d = FadbadTriangle2d(a_rot.segment<2>(0), b_rot.segment<2>(0), c_rot.segment<2>(0));
	if (tri_2d.is_inside(p_rot_2d)) {
		// then distance is just the z of the rotated triangle, but with x,y from the point
		// then rotate back to original frame
		Vector3b p_rot_proj;
		p_rot_proj << p_rot_2d, tri_rot_z;
		return rotation.fullPivLu().solve(Matrix3b::Identity())*p_rot_proj;
	} else {
//		// find closest point on 2d triangle
//		// connect back in 3d, then rotate back to original frame
		Vector2b closest_pt_2d = tri_2d.closest_point_to(p_rot_2d);
		Vector3b closest_pt_3d;
		closest_pt_3d << closest_pt_2d, tri_rot_z;
		return rotation.fullPivLu().solve(Matrix3b::Identity())*closest_pt_3d;
	}
}

/**
 * Finds the intersection of line segment (p0, p1) with the triangle
 * Returns false if no intersection
 */
bool FadbadTriangle3d::intersection_with_segment(const Vector3b& p0, const Vector3b& p1, Vector3b& intersection) {
	// halfplane: <n, x - origin> = 0
	// segment: x = t*(p1 - p0) + p0 | 0 <= t <= 1
	// t = -<n, p0 - origin> / <n, p1 - p0>

	Vector3b n = (b-a).cross(c-a);
	Vector3b origin = a; // arbitrary
	bdouble t = -n.dot(p0 - origin) / n.dot(p1 - p0);

	if ((t < 0) || (t > 1)) {
		return false; // segment doesn't cross plane
	}

	for(int i=0; i < 3; ++i) { intersection(i) = t*(p1(i) - p0(i)) + p0(i); }

	if (distance_to(intersection) > 1e-3) {
		return false; // segment crosses plane, but not the triangle
	}

	return true;
}

bdouble FadbadTriangle3d::area() {
	Matrix3b rotation = rotation_to_align_with({0,0,1});
	return fadbad_triangle_area((rotation*a).segment<2>(0),
			(rotation*b).segment<2>(0),
			(rotation*c).segment<2>(0));
}

/**
 * FadbadTriangle3d private methods
 */


inline Matrix3b axis_angle_to_rotation_matrix(const bdouble angle, const Vector3b& axis) {
	Matrix3b res;
	bdouble s = sin(angle);
	Vector3b sin_axis(s*axis(0), s*axis(1), s*axis(2));
	bdouble c = cos(angle);
	Vector3b cos1_axis((1-c)*axis(0), (1-c)*axis(1), (1-c)*axis(2));

	bdouble tmp;
	tmp = cos1_axis(0) * axis(1);
	res(0,1) = tmp - sin_axis(2);
	res(1,0) = tmp + sin_axis(2);

	tmp = cos1_axis(0) * axis(2);
	res(0,2) = tmp + sin_axis(1);
	res(2,0) = tmp - sin_axis(1);

	tmp = cos1_axis(1) * axis(2);
	res(1,2) = tmp - sin_axis(0);
	res(2,1) = tmp + sin_axis(0);

	Vector3b diag;
	for(int i=0; i < 3; ++i) { diag(i) = cos1_axis(i)*axis(i) + c; }

	res.diagonal() = diag;

	return res;
}

Matrix3b FadbadTriangle3d::rotation_to_align_with(const Vector3b& target) {
	Vector3b source = (b-a).cross(c-a);
	source.normalize();

	Matrix3b rotation = Matrix3b::Identity();

	bdouble dot = source.dot(target);
	if (!fadbad_isnan(dot)) {
		bdouble angle = acos(dot);
		if (!fadbad_isnan(angle)) {
			Vector3b cross = source.cross(target);
			bdouble cross_norm = cross.norm();
			if ((!fadbad_isnan(cross_norm)) && (cross_norm > bepsilon)) {
				cross /= cross_norm;
				rotation = axis_angle_to_rotation_matrix(angle, cross);
			}
		}
	}

	return rotation;
}
