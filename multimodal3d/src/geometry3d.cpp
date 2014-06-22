#include "../include/geometry3d.h"

/**
 * General functions
 */

inline bool halfplane_contains_point(const Vector3d& origin, const Vector3d& normal, const Vector3d& x) {
	return (normal.dot(x - origin) >= epsilon);
}

/**
 * Beam3d public methods
 */

bool Beam3d::is_inside(const Vector3d& p) {
	bool inside = true;

	inside &= halfplane_contains_point(base, (a-base).cross(d-base), p);
	inside &= halfplane_contains_point(base, (b-base).cross(a-base), p);
	inside &= halfplane_contains_point(base, (c-base).cross(b-base), p);
	inside &= halfplane_contains_point(base, (d-base).cross(c-base), p);
	inside &= halfplane_contains_point(a, (b-a).cross(d-a), p);

	return inside;
}

void Beam3d::plot(rave::EnvironmentBasePtr env) {
	Vector3d color(1,0,0);
//	rave_utils::plot_segment(env, base, a, color);
//	rave_utils::plot_segment(env, base, b, color);
//	rave_utils::plot_segment(env, base, c, color);
//	rave_utils::plot_segment(env, base, d, color);

	rave_utils::plot_segment(env, a, b, color);
	rave_utils::plot_segment(env, b, c, color);
	rave_utils::plot_segment(env, c, d, color);
	rave_utils::plot_segment(env, d, a, color);
}

/**
 * Triangle3d public methods
 */


/**
 * Rotates all points to be aligned with z-axis
 * Then projects to 2d, combining result back into 3d
 */
Vector3d Triangle3d::closest_point_to(const Vector3d& p) {
	Matrix3d rotation = rotation_to_align_with({0,0,1});

	Vector3d a_rot, b_rot, c_rot, p_rot;
	a_rot = rotation*a;
	b_rot = rotation*b;
	c_rot = rotation*c;
	p_rot = rotation*p;

	double tri_rot_z = a_rot(2);

	Vector2d p_rot_2d;
	p_rot_2d << p_rot(0), p_rot(1);

	Triangle2d tri_2d = Triangle2d(a_rot.segment<2>(0), b_rot.segment<2>(0), c_rot.segment<2>(0));
	if (tri_2d.is_inside(p_rot_2d)) {
		// then distance is just the z of the rotated triangle, but with x,y from the point
		// then rotate back to original frame
		Vector3d p_rot_proj;
		p_rot_proj << p_rot_2d, tri_rot_z;
		return rotation.inverse()*p_rot_proj;
	} else {
		// find closest point on 2d triangle
		// connect back in 3d, then rotate back to original frame
		Vector2d closest_pt_2d = tri_2d.closest_point_to(p_rot_2d);
		Vector3d closest_pt_3d;
		closest_pt_3d << closest_pt_2d, tri_rot_z;
		return rotation.inverse()*closest_pt_3d;
	}
}

double Triangle3d::area() {
	Matrix3d rotation = rotation_to_align_with({0,0,1});
	return triangle_area((rotation*a).segment<2>(0),
			(rotation*b).segment<2>(0),
			(rotation*c).segment<2>(0));
}

void Triangle3d::plot(rave::EnvironmentBasePtr env) {
	Vector3d color(1,0,0);
	rave_utils::plot_segment(env, a, b, color);
	rave_utils::plot_segment(env, b, c, color);
	rave_utils::plot_segment(env, c, a, color);
}

/**
 * Triangle3d private methods
 */

Matrix3d Triangle3d::rotation_to_align_with(const Vector3d& target) {
	Vector3d source = (b-a).cross(c-a);
	source.normalize();

	Matrix3d rotation = Matrix3d::Identity();

	double dot = source.dot(target);
	if (!isnan(dot)) {
		double angle = acos(dot);
		if (!isnan(angle)) {
			Vector3d cross = source.cross(target);
			double cross_norm = cross.norm();
			if ((!isnan(cross_norm)) && (cross_norm > epsilon)) {
				cross /= cross_norm;
				rotation = Eigen::AngleAxis<double>(angle, cross).toRotationMatrix();
			}
		}
	}

	return rotation;
}
