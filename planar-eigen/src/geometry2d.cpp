#include "../include/geometry2d.h"

/**
 * Line constructors
 */

Line::Line(const Segment& seg) {
	d = seg.p1 - seg.p0;
	o = seg.p0;
}

/**
 * Line public methods
 */

bool Line::intersection(const Line& other, Vector2d& intersection) {
	// s*d + o = t*d_other + o_other
	Matrix2d A;
	A << d, other.d; // [d d_other]
	Vector2d b = other.o - o;

	if (fabs(A.determinant()) < epsilon) {
		return false;
	}

	Vector2d soln = A.lu().solve(b);
	double s = soln(0);

	intersection = s*d + o;
	return true;
}

bool Line::intersection(const Segment& seg, Vector2d& intersection) {
	// v = seg.p1 - seg.p0
	// s*d + o = t*v + seg.p0

	Vector2d v = seg.p1 - seg.p0;

	Matrix2d A;
	A << d, v;
	Vector2d b = seg.p0 - o;

	if (fabs(A.determinant()) < epsilon) {
		return false;
	}

	Vector2d soln = A.lu().solve(b);
	double s = soln(0);
	double t = -soln(1);

	intersection = s*d + o;

	return ((-epsilon <= t) && (t <= 1+epsilon));
//	return seg.within_bounding_rect(intersection);
}

double Line::distance_to(const Vector2d& x) {
	// y = o - x
	// min_{t} ||t*d + y||_{2}
	Vector2d y = o - x;
	double t = -d.dot(y)/d.dot(d);

	return (t*d + y).norm();
}

/**
 * Halfspace public methods
 */

bool Halfspace::contains(const Vector2d& x) {
	return (n.dot(x - o) >= epsilon);
}

bool Halfspace::contains_part(const Segment& seg, Segment& seg_part) {
	Line split({-n(1), n(0)}, o);
	Vector2d intersection;

	if (split.intersection(seg, intersection)) {
		// segment crosses the half-space
		if (contains(seg.p0)) {// && norm(seg.p1 - intersection, 2) > epsilon) {
			seg_part = Segment(seg.p0, intersection);
			return true;
		} else if (contains(seg.p1)) {// && norm(seg.p0 - intersection, 2) > epsilon) {
			seg_part = Segment(seg.p1, intersection);
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
 * Segment public methods
 */

bool Segment::within_bounding_rect(const Vector2d& p) const {
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

bool Segment::intersection(const Segment& other, Vector2d& intersection) {
	Vector2d p0_other = other.p0, p1_other = other.p1;

	// w = p1 - p0
	// v = p1_other - p0_other
	// s*w + p0 = t*v + p_other

	Vector2d w = p1 - p0;
	Vector2d v = p1_other - p0_other;

	Matrix2d A;
	A << w, v;
	Vector2d b = p0_other - p0;

	if (fabs(A.determinant()) < epsilon) {
		return false;
	}

	Vector2d soln = A.lu().solve(b);
	double s = soln(0);
	double t = -soln(1);

	intersection = s*w + p0;

	return ((-epsilon <= s) && (s <= 1+epsilon) && (-epsilon <= t) && (t <= 1+epsilon));
//	return (this->within_bounding_rect(intersection) && (other.within_bounding_rect(intersection)));
}

Vector2d Segment::closest_point_to(const Vector2d& x) {
	// min_{0<=t<=1} ||t*(p1-p0) + p0 - x||_{2}^{2}
	Vector2d v = p1 - p0;
	Vector2d b = p0 - x;

	double t = -(v.dot(b)) / v.dot(v);
	if ((0 <= t) && (t <= 1)) {
		Vector2d intersection = t*(p1 - p0) + p0;
		return intersection;
	} else {
		if ((x - p0).norm() < (x - p1).norm()) {
			return p0;
		} else {
			return p1;
		}
	}
}

double Segment::distance_to(const Vector2d& x) {
	return (closest_point_to(x) - x).norm();
}


/**
 * Beam constructor
 */

Beam::Beam(const Vector2d& base_pt, const Vector2d& a_pt, const Vector2d& b_pt) : base(base_pt) {
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

	Vector2d left_intersect, top_intersect, right_intersect;

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
		Vector2d p_inside = (is_inside(s.p0)) ? s.p0 : s.p1;
		Vector2d top_projection_intersect;
//		bool should_intersect = Line(p_inside - base, base).intersection(top, top_projection_intersect);
		bool should_intersect = Line(p_inside - base, base).intersection(Line(top), top_projection_intersect);
		assert(should_intersect);
		new_beams.push_back(Beam(base, right_intersect, p_inside));
		new_beams.push_back(Beam(base, top_projection_intersect, b));
	} else if (is_intersect_top) {
		Vector2d p_inside = (is_inside(s.p0)) ? s.p0 : s.p1;
		Vector2d top_projection_intersect;
//		bool should_intersect = Line(p_inside - base, base).intersection(top, top_projection_intersect);
		bool should_intersect = Line(p_inside - base, base).intersection(Line(top), top_projection_intersect);
		assert(should_intersect);
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
		Vector2d p_inside = (is_inside(s.p0)) ? s.p0 : s.p1;
		Vector2d top_projection_intersect;
//		bool should_intersect = Line(p_inside - base, base).intersection(top, top_projection_intersect);
		bool should_intersect = Line(p_inside - base, base).intersection(Line(top), top_projection_intersect);
		assert(should_intersect);
		new_beams.push_back(Beam(base, a, top_projection_intersect));
		new_beams.push_back(Beam(base, p_inside, left_intersect));
	} else if (is_inside(s.p0) && is_inside(s.p1)) {
		Vector2d right_pt = (s.p0(0) > s.p1(0)) ? s.p0 : s.p1;
		Vector2d left_pt = (s.p0(0) <= s.p1(0)) ? s.p0 : s.p1;

		Vector2d rtop_projection_intersect, ltop_projection_intersect;
//		bool should_intersect = Line(right_pt - base, base).intersection(top, rtop_projection_intersect);
		bool should_intersect = Line(right_pt - base, base).intersection(Line(top), rtop_projection_intersect);
		assert(should_intersect);
//		should_intersect = Line(left_pt - base, base).intersection(top, ltop_projection_intersect);
		should_intersect = Line(left_pt - base, base).intersection(Line(top), ltop_projection_intersect);
		assert(should_intersect);

		new_beams.push_back(Beam(base, a, rtop_projection_intersect));
		new_beams.push_back(Beam(base, right_pt, left_pt));
		new_beams.push_back(Beam(base, ltop_projection_intersect, b));
	} else {
		new_beams.push_back(*this);
	}

	return new_beams;
}

double Beam::signed_distance(const Vector2d& x) {
	// sd positive if outside field of view
	bool inside = is_inside(x);

	std::vector<Segment> segments = {right_segment(), top_segment(), left_segment()};
	double sd = (inside) ? -INFINITY : INFINITY;
	for(int i=0; i < segments.size(); ++i) {
		if (inside) {
			sd = std::max(sd, -segments[i].distance_to(x));
		} else {
			sd = std::min(sd, segments[i].distance_to(x));
		}
	}

	return sd;
}

bool Beam::is_inside(const Vector2d& p) {
	double total_area = area();
	double area0 = Beam(base, a, p).area();
	double area1 = Beam(base, b, p).area();
	double area2 = Beam(a, b, p).area();

	bool is_correct_area = fabs(total_area - (area0 + area1 + area2)) < epsilon;

	double min_dist_to_side = std::min(right_segment().distance_to(p),
			std::min(top_segment().distance_to(p),
			left_segment().distance_to(p)));
	bool is_away_from_side = min_dist_to_side > epsilon;

	return (is_correct_area && is_away_from_side);
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

double signed_distance(const Vector2d& p, std::vector<Beam>& beams) {
	bool is_inside = false;
	for(int i=0; i < beams.size(); ++i) {
		is_inside |= beams[i].is_inside(p);
	}
	double sd_sign = (is_inside) ? -1 : 1;

	std::vector<Segment> border = beams_border(beams);
	double dist = INFINITY;
	for(int i=0; i < border.size(); ++i) {
		dist = std::min(dist, border[i].distance_to(p));
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

/**
 * Find closest point to the top segment
 * that is on the bottom side of the top segment
 * (serves as a conservative approximation of convex region)
 */
void truncate_belief(const std::vector<Beam>& beams, const Vector2d& cur_mean, const Matrix2d& cur_cov,
		Vector2d& out_mean, Matrix2d& out_cov) {
	Vector2d top_right = beams[0].a;
	Vector2d top_left = beams.back().b;
	Segment top_seg(top_right, top_left);
	Line top_line(top_right, top_left);

	Vector2d max_point = top_right;
	double max_dist = 0.0;

	if (fabs(top_right(0) - top_left(0)) < epsilon) {
		out_mean = cur_mean;
		out_cov = cur_cov;
		return;
	}

	Vector2d intersection;
	double tmp_dist = -INFINITY;
	for(int i=0; i < beams.size(); ++i) {
		const Beam &beam = beams[i];

		if ((!top_seg.intersection(Segment(beam.base, beam.a), intersection)) &&
				((tmp_dist = top_line.distance_to(beam.a)) > max_dist)) {
			max_dist = tmp_dist;
			max_point = beam.a;
		}

		if ((!top_seg.intersection(Segment(beam.base, beam.b), intersection)) &&
				((tmp_dist = top_line.distance_to(beam.b)) > max_dist)) {
			max_dist = tmp_dist;
			max_point = beam.b;
		}
	}

	// normal vector to top
	Vector2d c = {top_right(1) - top_left(1), top_left(0) - top_right(0)};

	// TODO: MY ADDITION
//	std::vector<Beam> tmp_beams = beams;
//	double n_sign = (signed_distance(cur_mean, tmp_beams) > 0) ? 1 : -1;
//	c *= n_sign;

	double d = c.dot(max_point);
	truncate_gaussian(c, d, cur_mean, cur_cov, out_mean, out_cov);
}

void my_truncate_belief(const std::vector<Beam>& beams, const Vector2d& cur_mean, const Matrix2d& cur_cov,
		Vector2d& out_mean, Matrix2d& out_cov) {
	std::vector<Beam> shifted_beams(beams);
	// shift everything to be centered around cur_mean
	// transform all coordinates by inv(U), where cur_cov = U*U.t() (i.e. cholesky)
	Matrix2d U = cur_cov.llt().matrixL();
	Matrix2d U_inv = U.inverse();
	for(int i=0; i < shifted_beams.size(); ++i) {
		Beam& b = shifted_beams[i];
		b.base = U_inv*(b.base - cur_mean);
		b.a = U_inv*(b.a - cur_mean);
		b.b = U_inv*(b.b - cur_mean);
	}


	std::vector<Segment> border = beams_border(shifted_beams);

//	// remove side edges
//	std::vector<Segment> new_border;
//	for(int i=1; i < border.size()-1; ++i) {
//		new_border.push_back(border[i]);
//	}
//	border = new_border;

	Vector2d mean = Vector2d::Zero();
	Matrix2d cov = Matrix2d::Identity();

	Vector2d delta_mean_total = Vector2d::Zero();
	Matrix2d delta_cov_total = Matrix2d::Zero();

	while(border.size() > 0) {
		// find closest point p on geometry (i.e. border) to the origin
		// the valid space is defined by the hyperplane with normal n = -p
		Vector2d p;
		double min_dist = INFINITY;
		for(int i=0; i < border.size(); ++i) {
			if (border[i].distance_to(Vector2d::Zero()) < min_dist) {
				p = border[i].closest_point_to(Vector2d::Zero());
				min_dist = p.norm();
			}
		}
		Vector2d n = -p;

		// truncate gaussian w.r.t. to -n (the complement space that we want to truncate)
		truncate_gaussian(-n, n.dot(n), mean, cov, out_mean, out_cov);
		mean = out_mean;
		cov = out_cov;

//		my_truncate_gaussian(-n, dot(n, n), mean, cov, delta_mean_total, delta_cov_total);

		// prune all geometry in infeasible space
		Halfspace h(-p, p); // TODO: h(n, p) or h(-p, p)
		std::vector<Segment> new_border;
		for(int i=0; i < border.size(); ++i) {
			Segment seg_part(Vector2d::Zero(), Vector2d::Zero());
			if (h.contains_part(border[i], seg_part)) {
				new_border.push_back(seg_part);
			}
		}
		border = new_border;
		// repeat while still points left

//		if (border.size() > 0 ) {
//			std::cout << "multiple iterations!\n";
//			std::cin.ignore();
//		}
	}

//	out_mean = U*(delta_mean_total) + cur_mean;
//	out_cov = U*(eye<Matrix2d>(DIM,DIM) + delta_cov_total)*U.t();

	out_mean = U*mean + cur_mean;
	out_cov = U*cov*U.transpose();
}

void truncate_gaussian(const Vector2d& c, double d, const Vector2d& cur_mean, const Matrix2d& cur_cov,
		Vector2d& out_mean, Matrix2d& out_cov) {
	double y_mean = c.dot(cur_mean);
	double y_var = c.dot(cur_cov*c);

	double y_new_mean, y_new_var;
	truncate_univariate_gaussian(d, y_mean, y_var, y_new_mean, y_new_var);
//	my_truncate_univariate_gaussian(d, y_mean, y_var, y_new_mean, y_new_var);

	Vector2d xy_cov = cur_cov*c;
	Vector2d L = xy_cov / y_var;

	out_mean = cur_mean + L*(y_new_mean - y_mean);
	out_cov = cur_cov + (y_new_var/y_var - 1.0) * (L*xy_cov.transpose());
}

void my_truncate_gaussian(const Vector2d& c, double d, const Vector2d& cur_mean, const Matrix2d& cur_cov,
		Vector2d& delta_mean_total, Matrix2d& delta_cov_total) {
	double y_mean = c.dot(cur_mean);
	double y_var = c.dot(cur_cov*c);

	double y_new_mean, y_new_var;
//	truncate_univariate_gaussian(d, y_mean, y_var, y_new_mean, y_new_var);
	my_truncate_univariate_gaussian(d, y_mean, y_var, y_new_mean, y_new_var);

	Vector2d xy_cov = cur_cov*c;
	Vector2d L = xy_cov / y_var;

	delta_mean_total += L*(y_mean - y_new_mean);
//	delta_cov_total += (y_new_var/y_var - 1.0) * (L*xy_cov.t());
	delta_cov_total += ((cur_cov*c) / y_var)*(y_var - y_new_var)*((c.transpose()*cur_cov)/ y_var);
}


void truncate_univariate_gaussian(const double x, const double cur_mean, const double cur_var,
		double& out_mean, double& out_var) {
	double sd = sqrt(cur_var);
	double y = (x - cur_mean) / sd;
	double z = pdf(standard_normal, y) / cdf(standard_normal, y);

	out_mean = cur_mean - z*sd;
	out_var = cur_var*(1.0 - y*z - z*z);
}

void my_truncate_univariate_gaussian(const double x, const double cur_mean, const double cur_var,
		double& out_mean, double& out_var) {
	double sd = sqrt(cur_var);
	double alpha = (x - cur_mean) / sd;
	double lambda = pdf(standard_normal, alpha) / cdf(standard_normal, alpha);

	out_mean = cur_mean + lambda*sd;
	out_var = cur_var*(1 - lambda*lambda + alpha*lambda);
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
			Matrix<double,2,3> m;
			m << beams[i].base, beams[i].a, beams[i].b;
//			beams_pylist.append(planar_utils::eigen_to_ndarray(m));
		}

		plot_beams(beams_pylist, true);
	}
	catch(py::error_already_set const &)
	{
		PyErr_Print();
	}
}

}
