#include "fadbad-geometry2d.h"

/**
 * LineB constructors
 */

LineB::LineB(const SegmentB& seg) {
	d = seg.p1 - seg.p0;
	o = seg.p0;
}

/**
 * LineB public methods
 */

bool LineB::intersection(const LineB& other, Vector2b& intersection) {
	// s*d + o = t*d_other + o_other
	Matrix2b A;
	A << d, other.d; // [d d_other]
	Vector2b b = other.o - o;

	if (abs(A.determinant()) < bepsilon) {
		return false;
	}

	Vector2b soln = A.inverse()*b;

	bdouble s = soln(0);

	intersection << s*d(0) + o(0), s*d(1) + o(1);
	return true;
}

bool LineB::intersection(const SegmentB& seg, Vector2b& intersection) {
	// v = seg.p1 - seg.p0
	// s*d + o = t*v + seg.p0

	Vector2b v = seg.p1 - seg.p0;

	Matrix2b A;
	A << d, v;
	Vector2b b = seg.p0 - o;

	if (abs(A.determinant()) < bepsilon) {
		return false;
	}
	Vector2b soln = A.inverse()*b;

	bdouble s = soln(0);
	bdouble t = -soln(1);

	intersection << s*d(0) + o(0) , s*d(1) + o(1);

	return ((-bepsilon <= t) && (t <= 1+bepsilon));
//	return seg.within_bounding_rect(intersection);
}

bdouble LineB::distance_to(const Vector2b& x) {
	// y = o - x
	// min_{t} ||t*d + y||_{2}
	Vector2b y = o - x;
	bdouble t = -d.dot(y)/d.dot(d);

	Vector2b pt = {t*d(0) + y(0), t*d(1) + y(1)};

	return pt.norm();
}

/**
 * HalfspaceB public methods
 */

bool HalfspaceB::contains(const Vector2b& x) {
	return (n.dot(x - o) >= bepsilon);
}

bool HalfspaceB::contains_part(const SegmentB& seg, SegmentB& seg_part) {
	LineB split({-n(1), n(0)}, o);
	Vector2b intersection;

	if (split.intersection(seg, intersection)) {
		// segment crosses the half-space
		if (contains(seg.p0)) {// && norm(seg.p1 - intersection, 2) > bepsilon) {
			seg_part = SegmentB(seg.p0, intersection);
			return true;
		} else if (contains(seg.p1)) {// && norm(seg.p0 - intersection, 2) > bepsilon) {
			seg_part = SegmentB(seg.p1, intersection);
			return true;
		} else {
			return false;
		}

	} else {
		if (contains(seg.p0)) {
			seg_part = seg;
			return true;
		} else {
			return false;
		}
	}
}

/**
 * SegmentB public methods
 */

bool SegmentB::within_bounding_rect(const Vector2b& p) const {
	bdouble x_min, x_max, y_min, y_max;
	x_min = std::min(p0(0), p1(0));
	x_max = std::max(p0(0), p1(0));
	y_min = std::min(p0(1), p1(1));
	y_max = std::max(p0(1), p1(1));

	bool within = true;
	within &= p(0) > x_min - bepsilon;
	within &= p(0) < x_max + bepsilon;
	within &= p(1) > y_min - bepsilon;
	within &= p(1) < y_max + bepsilon;

	return within;
}

bool SegmentB::intersection(const SegmentB& other, Vector2b& intersection) {
	Vector2b p0_other = other.p0, p1_other = other.p1;

	// w = p1 - p0
	// v = p1_other - p0_other
	// s*w + p0 = t*v + p_other

	Vector2b w = p1 - p0;
	Vector2b v = p1_other - p0_other;

	Matrix2b A;
	A << w, v;
	Vector2b b = p0_other - p0;

	if (abs(A.determinant()) < bepsilon) {
		return false;
	}

	Vector2b soln = A.inverse()*b;
	bdouble s = soln(0);
	bdouble t = -soln(1);

	intersection << s*w(0) + p0(0), s*w(1) + p0(1);

	return ((-bepsilon <= s) && (s <= 1+bepsilon) && (-bepsilon <= t) && (t <= 1+bepsilon));
//	return (this->within_bounding_rect(intersection) && (other.within_bounding_rect(intersection)));
}

Vector2b SegmentB::closest_point_to(const Vector2b& x) {
	// min_{0<=t<=1} ||t*(p1-p0) + p0 - x||_{2}^{2}
	Vector2b v = p1 - p0;
	Vector2b b = p0 - x;

	bdouble t = -(v.dot(b)) / v.dot(v);
	if ((0 <= t) && (t <= 1)) {
		Vector2b intersection = {t*(p1(0) - p0(0)) + p0(0) , t*(p1(1) - p0(1)) + p0(1)};
		return intersection;
	} else {
		if ((x - p0).norm() < (x - p1).norm()) {
			return p0;
		} else {
			return p1;
		}
	}
}

bdouble SegmentB::distance_to(const Vector2b& x) {
	return (closest_point_to(x) - x).norm();
}


/**
 * BeamB constructor
 */

BeamB::BeamB(const Vector2b& base_pt, const Vector2b& a_pt, const Vector2b& b_pt) : base(base_pt) {
	a = (a_pt(0) > b_pt(0)) ? a_pt : b_pt;
	b = (a_pt(0) <= b_pt(0)) ? a_pt : b_pt;
}

/**
 * BeamB public methods
 */

std::vector<BeamB> BeamB::truncate(const SegmentB& s) {
	SegmentB left = left_segment();
	SegmentB top = top_segment();
	SegmentB right = right_segment();

	Vector2b left_intersect, top_intersect, right_intersect;

	bool is_intersect_right = right.intersection(s, right_intersect);
	bool is_intersect_top = top.intersection(s, top_intersect);
	bool is_intersect_left = left.intersection(s, left_intersect);

	std::vector<BeamB> new_beams;
	if (is_intersect_right && is_intersect_left) {
		new_beams.push_back(BeamB(base, right_intersect, left_intersect));
	} else if (is_intersect_right && is_intersect_top) {
		new_beams.push_back(BeamB(base, right_intersect, top_intersect));
		new_beams.push_back(BeamB(base, top_intersect, b));
	} else if (is_intersect_top && is_intersect_left) {
		new_beams.push_back(BeamB(base, a, top_intersect));
		new_beams.push_back(BeamB(base, top_intersect, left_intersect));
	} else if (is_intersect_right) {
		Vector2b p_inside = (is_inside(s.p0)) ? s.p0 : s.p1;
		Vector2b top_projection_intersect;
//		bool should_intersect = LineB(p_inside - base, base).intersection(top, top_projection_intersect);
		bool should_intersect = LineB(p_inside - base, base).intersection(LineB(top), top_projection_intersect);
		assert(should_intersect);
		new_beams.push_back(BeamB(base, right_intersect, p_inside));
		new_beams.push_back(BeamB(base, top_projection_intersect, b));
	} else if (is_intersect_top) {
		Vector2b p_inside = (is_inside(s.p0)) ? s.p0 : s.p1;
		Vector2b top_projection_intersect;
//		bool should_intersect = LineB(p_inside - base, base).intersection(top, top_projection_intersect);
		bool should_intersect = LineB(p_inside - base, base).intersection(LineB(top), top_projection_intersect);
		assert(should_intersect);
		if (top_intersect(0) > p_inside(0)) {
			new_beams.push_back(BeamB(base, a, top_intersect));
			new_beams.push_back(BeamB(base, top_intersect, p_inside));
			new_beams.push_back(BeamB(base, top_projection_intersect, b));
		} else {
			new_beams.push_back(BeamB(base, a, top_projection_intersect));
			new_beams.push_back(BeamB(base, p_inside, top_intersect));
			new_beams.push_back(BeamB(base, top_intersect, b));
		}
	} else if (is_intersect_left) {
		Vector2b p_inside = (is_inside(s.p0)) ? s.p0 : s.p1;
		Vector2b top_projection_intersect;
//		bool should_intersect = LineB(p_inside - base, base).intersection(top, top_projection_intersect);
		bool should_intersect = LineB(p_inside - base, base).intersection(LineB(top), top_projection_intersect);
		assert(should_intersect);
		new_beams.push_back(BeamB(base, a, top_projection_intersect));
		new_beams.push_back(BeamB(base, p_inside, left_intersect));
	} else if (is_inside(s.p0) && is_inside(s.p1)) {
		Vector2b right_pt = (s.p0(0) > s.p1(0)) ? s.p0 : s.p1;
		Vector2b left_pt = (s.p0(0) <= s.p1(0)) ? s.p0 : s.p1;

		Vector2b rtop_projection_intersect, ltop_projection_intersect;
//		bool should_intersect = LineB(right_pt - base, base).intersection(top, rtop_projection_intersect);
		bool should_intersect = LineB(right_pt - base, base).intersection(LineB(top), rtop_projection_intersect);
		assert(should_intersect);
//		should_intersect = LineB(left_pt - base, base).intersection(top, ltop_projection_intersect);
		should_intersect = LineB(left_pt - base, base).intersection(LineB(top), ltop_projection_intersect);
		assert(should_intersect);

		new_beams.push_back(BeamB(base, a, rtop_projection_intersect));
		new_beams.push_back(BeamB(base, right_pt, left_pt));
		new_beams.push_back(BeamB(base, ltop_projection_intersect, b));
	} else {
		new_beams.push_back(*this);
	}

	return new_beams;
}

bdouble BeamB::signed_distance(const Vector2b& x) {
	// sd positive if outside field of view
	bool inside = is_inside(x);

	std::vector<SegmentB> segments = {right_segment(), top_segment(), left_segment()};
	bdouble sd = (inside) ? -INFINITY : INFINITY;
	for(int i=0; i < segments.size(); ++i) {
		if (inside) {
			sd = std::max(sd, -segments[i].distance_to(x));
		} else {
			sd = std::min(sd, segments[i].distance_to(x));
		}
	}

	return sd;
}

bool BeamB::is_inside(const Vector2b& p) {
	bdouble total_area = area();
	bdouble area0 = BeamB(base, a, p).area();
	bdouble area1 = BeamB(base, b, p).area();
	bdouble area2 = BeamB(a, b, p).area();

	bool is_correct_area = abs(total_area - (area0 + area1 + area2)) < bepsilon;

	bdouble min_dist_to_side = std::min(right_segment().distance_to(p),
			std::min(top_segment().distance_to(p),
			left_segment().distance_to(p)));
	bool is_away_from_side = min_dist_to_side > bepsilon;

	return (is_correct_area && is_away_from_side);
}

/**
 * BeamB private methods
 */

bdouble BeamB::area() {
	return abs((base(0)*(a(1) - b(1)) + a(0)*(b(1) - base(1)) + b(0)*(base(1) - a(1))) / 2.0);
}

/**
 * Functions
 */

namespace fadbad_geometry2d {

bdouble signed_distance(const Vector2b& p, std::vector<BeamB>& beams) {
	bool is_inside = false;
	for(int i=0; i < beams.size(); ++i) {
		is_inside |= beams[i].is_inside(p);
	}
	bdouble sd_sign = (is_inside) ? -1 : 1;

	std::vector<SegmentB> border = beams_border(beams);
	bdouble dist = INFINITY;
	for(int i=0; i < border.size(); ++i) {
		dist = std::min(dist, border[i].distance_to(p));
	}

	return sd_sign*dist;
}

// NOTE: assumes beams are sorted from right to left
std::vector<SegmentB> beams_border(const std::vector<BeamB>& beams) {
	std::vector<SegmentB> segments;
	int num_beams = beams.size();

	segments.push_back(SegmentB(beams[0].base, beams[0].a));
	for(int i=0; i < num_beams; ++i) {
		segments.push_back(SegmentB(beams[i].a, beams[i].b));
		if (i < beams.size() - 1) {
			segments.push_back(SegmentB(beams[i].b, beams[i+1].a));
		}
	}
	segments.push_back(SegmentB(beams[num_beams-1].b, beams[num_beams-1].base));

	// filter out small segments
	std::vector<SegmentB> new_segments;
	for(int i=0; i < segments.size(); ++i) {
		if (segments[i].length() > bepsilon) {
			new_segments.push_back(segments[i]);
		}
	}
	segments = new_segments;

	return segments;
}


}
