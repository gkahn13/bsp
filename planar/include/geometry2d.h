#ifndef __GEOMETRY2D_H__
#define __GEOMETRY2D_H__

#include <Python.h>

#include <boost/python.hpp>
#include <boost/python/numeric.hpp>
#include <boost/python/tuple.hpp>
#include <boost/numpy.hpp>
#include <boost/filesystem.hpp>

namespace py = boost::python;
namespace np = boost::numpy;


#include <armadillo>
using namespace arma;

#include <assert.h>

#include "planar-utils.h"

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
	double distance_from(const vec& x);

	double length() { return norm(p1 - p0, 2); }
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

	double signed_distance(const vec& x);
	bool is_inside(const vec& p);

	double top_length() { return top_segment().length(); }

private:

	Segment right_segment() { return Segment(base, a); }
	Segment top_segment() { return Segment(a, b); }
	Segment left_segment() { return Segment(base, b); }

	double area();
};

namespace geometry2d {

double signed_distance(const vec& p, std::vector<Beam>& beams);

std::vector<Segment> beams_border(const std::vector<Beam>& beams);

void plot_beams(std::vector<Beam>& beams);

}

#endif
