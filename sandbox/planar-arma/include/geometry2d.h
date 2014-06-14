#ifndef __GEOMETRY2D_H__
#define __GEOMETRY2D_H__

#include <Python.h>

#include <boost/python.hpp>
#include <boost/python/numeric.hpp>
#include <boost/python/tuple.hpp>
#include <boost/numpy.hpp>
#include <boost/filesystem.hpp>
#include <boost/math/distributions.hpp>
const boost::math::normal_distribution<> standard_normal;

namespace py = boost::python;
namespace np = boost::numpy;

#include <armadillo>
using namespace arma;

#include <assert.h>

#include "planar-utils.h"

const double epsilon = 1e-5;

class Line;
class Halfspace;
class Segment;
class Beam;

/**
 * Vector direction + origin
 */
class Line {
public:
	Line(const vec& direction, const vec& origin) : d(direction), o(origin) { };
	Line(const Segment& seg);

	bool intersection(const Line& other, vec& intersection);
	bool intersection(const Segment& seg, vec& intersection);
	double distance_to(const vec& x);
private:
	vec d, o;
};

/**
 * Normal direction + origin
 */
class Halfspace {
public:
	Halfspace(const vec& normal, const vec& origin) : n(normal), o(origin) { };

	bool contains(const vec& x);
	bool contains_part(const Segment& seg, Segment& seg_part);
private:
	vec n, o;
};

class Segment {
public:
	vec p0, p1;

	Segment(const vec& point0, const vec& point1) : p0(point0), p1(point1) { };

	bool intersection(const Segment& other, vec& intersection);
	bool within_bounding_rect(const vec& p) const;
	vec closest_point_to(const vec& x);
	double distance_to(const vec& x);

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

void truncate_belief(const std::vector<Beam>& beams, const vec& cur_mean, const mat& cur_cov,
		vec& out_mean, mat& out_cov);

void my_truncate_belief(const std::vector<Beam>& beams, const vec& cur_mean, const mat& cur_cov,
		vec& out_mean, mat& out_cov);

void truncate_gaussian(const vec& c, double d, const vec& cur_mean, const mat& cur_cov,
		vec& out_mean, mat& out_cov);

void my_truncate_gaussian(const vec& c, double d, const vec& cur_mean, const mat& cur_cov,
		vec& delta_mean_total, mat& delta_cov_total);

void truncate_univariate_gaussian(const double x, const double cur_mean, const double cur_var,
		double& out_mean, double& out_var);

void my_truncate_univariate_gaussian(const double x, const double cur_mean, const double cur_var,
		double& out_mean, double& out_var);

void plot_beams(std::vector<Beam>& beams);

}

#endif
