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

	bool within = true;
	within &= max((seg.p0 - epsilon) < intersection);
	within &= max(intersection < (seg.p1 + epsilon));

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

//	bool within = true;
//	within &= max((p0 - epsilon) < intersection);
//	within &= max(intersection < (p1 + epsilon));
//	within &= max((p0_other - epsilon) < intersection);
//	within &= max(intersection < (p1_other + epsilon));
//
//	return within;
	return (this->within_bounding_rect(intersection) && (other.within_bounding_rect(intersection)));
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

	std::cout << "is_intersect_right: " << is_intersect_right << "\n";
	std::cout << "is_intersect_top: " << is_intersect_top << "\n";
	std::cout << "is_intersect_left: " << is_intersect_left << "\n\n";

	std::cout << "right_intersect: " << right_intersect.t();
	std::cout << "top_intersect: " << top_intersect.t();
	std::cout << "left_intersect: " << left_intersect.t();


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
		vec p_inside = (s.p0(0) > s.p1(0)) ? s.p1 : s.p0;
		vec top_projection_intersect;
		if (!Line(p_inside - base, base).intersection(top, top_projection_intersect)) { std::cerr << "Beam truncate: should intersect\n"; exit(-1); }
		new_beams.push_back(Beam(base, right_intersect, p_inside));
		new_beams.push_back(Beam(base, top_projection_intersect, b));
	} else if (is_intersect_left) {
		vec p_inside = (s.p0(0) < s.p1(0)) ? s.p1 : s.p0;
		vec top_projection_intersect;
		if (!Line(p_inside - base, base).intersection(top, top_projection_intersect)) { std::cerr << "Beam truncate: should intersect\n"; exit(-1); }
		new_beams.push_back(Beam(base, a, top_projection_intersect));
		new_beams.push_back(Beam(base, p_inside, left_intersect));
	} else if (is_inside(s.p0) && is_inside(s.p1)) {
		std::cout << "both are inside\n";
		vec right_pt = (s.p0(0) > s.p1(0)) ? s.p0 : s.p1;
		vec left_pt = (s.p0(0) <= s.p1(0)) ? s.p0 : s.p1;

		vec rtop_projection_intersect, ltop_projection_intersect;
		if (!Line(right_pt - base, base).intersection(top, rtop_projection_intersect)) { std::cerr << "Beam truncate: should intersect\n"; exit(-1); }
		if (!Line(left_pt - base, base).intersection(top, ltop_projection_intersect)) { std::cerr << "Beam truncate: should intersect\n"; exit(-1); }

		new_beams.push_back(Beam(base, a, rtop_projection_intersect));
		new_beams.push_back(Beam(base, right_pt, left_pt));
		new_beams.push_back(Beam(base, ltop_projection_intersect, b));
	} else {
		new_beams.push_back(*this);
	}

	return new_beams;
}

/**
 * Beam private methods
 */

double Beam::area() {
	return fabs((base(0)*(a(1) - b(1)) + a(0)*(b(1) - base(1)) + b(0)*(base(1) - a(1))) / 2.0);
}

bool Beam::is_inside(const vec& p) {
	double total_area = area();
	double area0 = Beam(base, a, p).area();
	double area1 = Beam(base, b, p).area();
	double area2 = Beam(a, b, p).area();

	return fabs(total_area - (area0 + area1 + area2)) < epsilon;
}
