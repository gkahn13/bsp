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

#include <Eigen/Eigen>
using namespace Eigen;

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
	Line(const Vector2d& direction, const Vector2d& origin) : d(direction), o(origin) { };
	Line(const Segment& seg);

	bool intersection(const Line& other, Vector2d& intersection);
	bool intersection(const Segment& seg, Vector2d& intersection);
	double distance_to(const Vector2d& x);
private:
	Vector2d d, o;
};

/**
 * Normal direction + origin
 */
class Halfspace {
public:
	Halfspace(const Vector2d& normal, const Vector2d& origin) : n(normal), o(origin) { };

	bool contains(const Vector2d& x);
	bool contains_part(const Segment& seg, Segment& seg_part);
private:
	Vector2d n, o;
};

class Segment {
public:
	Vector2d p0, p1;

	Segment(const Vector2d& point0, const Vector2d& point1) : p0(point0), p1(point1) { };

	bool intersection(const Segment& other, Vector2d& intersection);
	bool within_bounding_rect(const Vector2d& p) const;
	Vector2d closest_point_to(const Vector2d& x);
	double distance_to(const Vector2d& x);

	double length() { return (p1 - p0).norm(); }
};

/**
 * Triangle in which space inside the points
 * is in the FOV
 */
class Beam {
public:
	Vector2d base, a, b;

	// base, a, b Counter-Clockwise
	Beam(const Vector2d& base_pt, const Vector2d& a_pt, const Vector2d& b_pt);

	std::vector<Beam> truncate(const Segment& s);

	double signed_distance(const Vector2d& x);
	bool is_inside(const Vector2d& p);

	double top_length() { return top_segment().length(); }

private:

	Segment right_segment() { return Segment(base, a); }
	Segment top_segment() { return Segment(a, b); }
	Segment left_segment() { return Segment(base, b); }

	double area();
};

namespace geometry2d {

double signed_distance(const Vector2d& p, std::vector<Beam>& beams);

std::vector<Segment> beams_border(const std::vector<Beam>& beams);

void truncate_belief(const std::vector<Beam>& beams, const Vector2d& cur_mean, const Matrix2d& cur_cov,
		Vector2d& out_mean, Matrix2d& out_cov);

void my_truncate_belief(const std::vector<Beam>& beams, const Vector2d& cur_mean, const Matrix2d& cur_cov,
		Vector2d& out_mean, Matrix2d& out_cov);

void truncate_gaussian(const Vector2d& c, double d, const Vector2d& cur_mean, const Matrix2d& cur_cov,
		Vector2d& out_mean, Matrix2d& out_cov);

void my_truncate_gaussian(const Vector2d& c, double d, const Vector2d& cur_mean, const Matrix2d& cur_cov,
		Vector2d& delta_mean_total, Matrix2d& delta_cov_total);

void truncate_univariate_gaussian(const double x, const double cur_mean, const double cur_var,
		double& out_mean, double& out_var);

void my_truncate_univariate_gaussian(const double x, const double cur_mean, const double cur_var,
		double& out_mean, double& out_var);

void plot_beams(std::vector<Beam>& beams);

}

#endif
