#include "../include/geometry2d.h"

/**
 * Line public methods
 */

bool Line::intersection(const Segment& seg, vec& intersection) {
	// v = seg.p1 - seg.p0
	// s*d + o = t*v + seg.p0

	vec v = seg.p1 - seg.p0;

	mat A = join_horiz(d, v);
	mat b = seg.p0 - o;

	if (fabs(det(A)) < epsilon) {
		return false;
	}

	vec soln = solve(A, b);
	double s = soln(0);

	intersection = s*d + o;

	return seg.within_bounding_rect(intersection);
}

/**
 * Segment public methods
 */

bool Segment::within_bounding_rect(const vec& p) const {
	double x_min, x_max, y_min, y_max;
	x_min = std::min(p0(0), p1(0));
	x_max = std::max(p0(0), p1(0));
	y_min = std::min(p0(1), p1(1));
	y_max = std::max(p0(1), p1(1));

	bool within = true;
	within &= p(0) > x_min - epsilon;
	within &= p(0) < x_max + epsilon;
	within &= p(1) > y_min - epsilon;
	within &= p(1) < y_max + epsilon;

	return within;
}

bool Segment::intersection(const Segment& other, vec& intersection) {
	vec p0_other = other.p0, p1_other = other.p1;

	// w = p1 - p0
	// v = p1_other - p0_other
	// s*w + p0 = t*v + p_other

	vec w = p1 - p0;
	vec v = p1_other - p0_other;

	mat A = join_horiz(w, v);
	mat b = p0_other - p0;

	if (fabs(det(A)) < epsilon) {
		return false;
	}

	vec soln = solve(A, b);
	double s = soln(0);

	intersection = s*w + p0;

	return (this->within_bounding_rect(intersection) && (other.within_bounding_rect(intersection)));
}

double Segment::distance_from(const vec& x) {
	// min_{0<=t<=1} ||t*(p1-p0) + p0 - x||_{2}^{2}
	vec v = p1 - p0;
	vec b = p0 - x;

	double t = -trace((v.t()*b)/(v.t()*v));
	if ((0 <= t) && (t <= 1)) {
		vec intersection = t*(p1 - p0) + p0;
		return norm(x - intersection, 2);
	} else {
		return std::min(norm(x - p0, 2), norm(x - p1, 2));
	}
}

/**
 * Beam constructor
 */

Beam::Beam(const vec& base_pt, const vec& a_pt, const vec& b_pt) : base(base_pt) {
	a = (a_pt(0) > b_pt(0)) ? a_pt : b_pt;
	b = (a_pt(0) <= b_pt(0)) ? a_pt : b_pt;
}

/**
 * Beam public methods
 */

std::vector<Beam> Beam::truncate(const Segment& s) {
	Segment left = left_segment();
	Segment top = top_segment();
	Segment right = right_segment();

	vec left_intersect, top_intersect, right_intersect;

	bool is_intersect_right = right.intersection(s, right_intersect);
	bool is_intersect_top = top.intersection(s, top_intersect);
	bool is_intersect_left = left.intersection(s, left_intersect);

	std::vector<Beam> new_beams;
	if (is_intersect_right && is_intersect_left) {
		new_beams.push_back(Beam(base, right_intersect, left_intersect));
	} else if (is_intersect_right && is_intersect_top) {
		new_beams.push_back(Beam(base, right_intersect, top_intersect));
		new_beams.push_back(Beam(base, top_intersect, b));
	} else if (is_intersect_top && is_intersect_left) {
		new_beams.push_back(Beam(base, a, top_intersect));
		new_beams.push_back(Beam(base, top_intersect, left_intersect));
	} else if (is_intersect_right) {
		vec p_inside = (is_inside(s.p0)) ? s.p0 : s.p1;
		vec top_projection_intersect;
		assert(Line(p_inside - base, base).intersection(top, top_projection_intersect));
		new_beams.push_back(Beam(base, right_intersect, p_inside));
		new_beams.push_back(Beam(base, top_projection_intersect, b));
	} else if (is_intersect_top) {
		vec p_inside = (is_inside(s.p0)) ? s.p0 : s.p1;
		vec top_projection_intersect;
		assert(Line(p_inside - base, base).intersection(top, top_projection_intersect));
		if (top_intersect(0) > p_inside(0)) {
			new_beams.push_back(Beam(base, a, top_intersect));
			new_beams.push_back(Beam(base, top_intersect, p_inside));
			new_beams.push_back(Beam(base, top_projection_intersect, b));
		} else {
			new_beams.push_back(Beam(base, a, top_projection_intersect));
			new_beams.push_back(Beam(base, p_inside, top_intersect));
			new_beams.push_back(Beam(base, top_intersect, b));
		}
	} else if (is_intersect_left) {
		vec p_inside = (is_inside(s.p0)) ? s.p0 : s.p1;
		vec top_projection_intersect;
		assert(Line(p_inside - base, base).intersection(top, top_projection_intersect));
		new_beams.push_back(Beam(base, a, top_projection_intersect));
		new_beams.push_back(Beam(base, p_inside, left_intersect));
	} else if (is_inside(s.p0) && is_inside(s.p1)) {
		vec right_pt = (s.p0(0) > s.p1(0)) ? s.p0 : s.p1;
		vec left_pt = (s.p0(0) <= s.p1(0)) ? s.p0 : s.p1;

		vec rtop_projection_intersect, ltop_projection_intersect;
		assert(Line(right_pt - base, base).intersection(top, rtop_projection_intersect));
		assert(Line(left_pt - base, base).intersection(top, ltop_projection_intersect));

		new_beams.push_back(Beam(base, a, rtop_projection_intersect));
		new_beams.push_back(Beam(base, right_pt, left_pt));
		new_beams.push_back(Beam(base, ltop_projection_intersect, b));
	} else {
		new_beams.push_back(*this);
	}

	return new_beams;
}

double Beam::signed_distance(const vec& x) {
	// sd positive if outside field of view
	bool inside = (is_inside(x)) ? -1 : 1;

	std::vector<Segment> segments = {right_segment(), top_segment(), left_segment()};
	double sd = (inside) ? -INFINITY : INFINITY;
	for(int i=0; i < segments.size(); ++i) {
		if (inside) {
			sd = std::max(sd, -segments[i].distance_from(x));
		} else {
			sd = std::min(sd, segments[i].distance_from(x));
		}
	}

	return sd;
}

bool Beam::is_inside(const vec& p) {
	double total_area = area();
	double area0 = Beam(base, a, p).area();
	double area1 = Beam(base, b, p).area();
	double area2 = Beam(a, b, p).area();

	return fabs(total_area - (area0 + area1 + area2)) < epsilon;
}

/**
 * Beam private methods
 */

double Beam::area() {
	return fabs((base(0)*(a(1) - b(1)) + a(0)*(b(1) - base(1)) + b(0)*(base(1) - a(1))) / 2.0);
}


/**
 * Functions
 */

namespace geometry2d {

double signed_distance(const vec& p, std::vector<Beam>& beams) {
	bool is_inside = false;
	for(int i=0; i < beams.size(); ++i) {
		is_inside |= beams[i].is_inside(p);
	}
	double sd_sign = (is_inside) ? -1 : 1;

	std::vector<Segment> border = beams_border(beams);
	double dist = INFINITY;
	for(int i=0; i < border.size(); ++i) {
		dist = std::min(dist, border[i].distance_from(p));
	}

	return sd_sign*dist;
}

// NOTE: assumes beams are sorted from right to left
std::vector<Segment> beams_border(const std::vector<Beam>& beams) {
	std::vector<Segment> segments;
	int num_beams = beams.size();

	segments.push_back(Segment(beams[0].base, beams[0].a));
	for(int i=0; i < num_beams; ++i) {
		segments.push_back(Segment(beams[i].a, beams[i].b));
		if (i < beams.size() - 1) {
			segments.push_back(Segment(beams[i].b, beams[i+1].a));
		}
	}
	segments.push_back(Segment(beams[num_beams-1].b, beams[num_beams-1].base));

	// filter out small segments
	std::vector<Segment> new_segments;
	for(int i=0; i < segments.size(); ++i) {
		if (segments[i].length() > epsilon) {
			new_segments.push_back(segments[i]);
		}
	}
	segments = new_segments;

	return segments;
}

void plot_beams(std::vector<Beam>& beams) {
	try {
		Py_Initialize();
		np::initialize();

		py::numeric::array::set_module_and_type("numpy", "ndarray");

		std::string working_dir = boost::filesystem::current_path().normalize().string();
		std::string bsp_dir = working_dir.substr(0,working_dir.find("bsp"));
		std::string planar_dir = bsp_dir + "bsp/planar";

		py::object main_module = py::import("__main__");
		py::object main_namespace = main_module.attr("__dict__");
		py::exec("import sys, os", main_namespace);
		py::exec(py::str("sys.path.append('"+planar_dir+"')"), main_namespace);
		py::object plot_mod = py::import("plot_planar");
		py::object plot_beams = plot_mod.attr("plot_beams");

		py::list beams_pylist;
		for(int i=0; i < beams.size(); ++i) {
			mat m = join_horiz(beams[i].base, beams[i].a);
			m = join_horiz(m, beams[i].b);
			beams_pylist.append(planar_utils::arma_to_ndarray(m));
		}

		plot_beams(beams_pylist, true);
	}
	catch(py::error_already_set const &)
	{
		PyErr_Print();
	}
}

}
