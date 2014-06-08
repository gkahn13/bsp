#ifndef __GEOMETRY2D_H__
#define __GEOMETRY2D_H__

#include <armadillo>
using namespace arma;

const double epsilon = 1e-4;

class Line;
class Segment;
class Beam;

/**
 * Vector direction + origin
 */
class Line {
public:
	Line(const vec& direction, const vec& origin) : d(direction), o(origin) { };

	bool intersection(const Segment& seg, vec& intersection);
private:
	vec d, o;
};

class Segment {
public:
	vec p0, p1;

	Segment(const vec& point0, const vec& point1) : p0(point0), p1(point1) { };

	bool intersection(const Segment& other, vec& intersection);
	bool within_bounding_rect(const vec& p) const;
};

/**
 * Triangle in which space inside the points
 * is in the FOV
 */
class Beam {
public:
	vec base, a, b;

	// base, a, b Counter-Clockwise
	Beam(const vec& base_pt, const vec& a_pt, const vec& b_pt);

	std::vector<Beam> truncate(const Segment& s);

private:

	Segment right_segment() { return Segment(base, a); }
	Segment top_segment() { return Segment(a, b); }
	Segment left_segment() { return Segment(base, b); }

	double area();
	bool is_inside(const vec& p);
};

#endif
