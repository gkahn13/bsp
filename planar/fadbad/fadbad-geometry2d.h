#ifndef __FADBAD_GEOMETRY2D_H__
#define __FADBAD_GEOMETRY2D_H__

#include <Eigen/Eigen>
using namespace Eigen;

#include "fadbad-utils.h"

typedef Matrix<bdouble, 2, 1> Vector2b;
typedef Matrix<bdouble, 2, 2> Matrix2b;

class LineB;
class HalfspaceB;
class SegmentB;
class BeamB;

/**
 * Vector direction + origin
 */
class LineB {
public:
	LineB(const Vector2b& direction, const Vector2b& origin) : d(direction), o(origin) { };
	LineB(const SegmentB& seg);

	bool intersection(const LineB& other, Vector2b& intersection);
	bool intersection(const SegmentB& seg, Vector2b& intersection);
	bdouble distance_to(const Vector2b& x);
private:
	Vector2b d, o;
};

/**
 * Normal direction + origin
 */
class HalfspaceB {
public:
	HalfspaceB(const Vector2b& normal, const Vector2b& origin) : n(normal), o(origin) { };

	bool contains(const Vector2b& x);
	bool contains_part(const SegmentB& seg, SegmentB& seg_part);
private:
	Vector2b n, o;
};

class SegmentB {
public:
	Vector2b p0, p1;

	SegmentB(const Vector2b& point0, const Vector2b& point1) : p0(point0), p1(point1) { };

	bool intersection(const SegmentB& other, Vector2b& intersection);
	bool within_bounding_rect(const Vector2b& p) const;
	Vector2b closest_point_to(const Vector2b& x);
	bdouble distance_to(const Vector2b& x);

	bdouble length() { return (p1 - p0).norm(); }
};

/**
 * Triangle in which space inside the points
 * is in the FOV
 */
class BeamB {
public:
	Vector2b base, a, b;

	// base, a, b Counter-Clockwise
	BeamB(const Vector2b& base_pt, const Vector2b& a_pt, const Vector2b& b_pt);

	std::vector<BeamB> truncate(const SegmentB& s);

	bdouble signed_distance(const Vector2b& x);
	bool is_inside(const Vector2b& p);

	bdouble top_length() { return top_segment().length(); }

private:

	SegmentB right_segment() { return SegmentB(base, a); }
	SegmentB top_segment() { return SegmentB(a, b); }
	SegmentB left_segment() { return SegmentB(base, b); }

	bdouble area();
};

namespace fadbad_geometry2d {

bdouble signed_distance(const Vector2b& p, std::vector<BeamB>& beams);

std::vector<SegmentB> beams_border(const std::vector<BeamB>& beams);

}

#endif
