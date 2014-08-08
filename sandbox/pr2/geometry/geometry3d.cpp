#include "geometry3d.h"

/**
 * General functions
 */

inline bool halfplane_contains_point(const Vector3d& origin, const Vector3d& normal, const Vector3d& x) {
	return (normal.dot(x - origin) >= -epsilon);
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

bool Beam3d::is_crossed_by(const Vector3d& p0, const Vector3d& p1) {
	bool doesnt_cross = false;

	doesnt_cross |= (halfplane_contains_point(base, -(a-base).cross(d-base), p0) && halfplane_contains_point(base, -(a-base).cross(d-base), p1));
	doesnt_cross |= (halfplane_contains_point(base, -(b-base).cross(a-base), p0) && halfplane_contains_point(base, -(b-base).cross(a-base), p1));
	doesnt_cross |= (halfplane_contains_point(base, -(c-base).cross(b-base), p0) && halfplane_contains_point(base, -(c-base).cross(b-base), p1));
	doesnt_cross |= (halfplane_contains_point(base, -(d-base).cross(c-base), p0) && halfplane_contains_point(base, -(d-base).cross(c-base), p1));

	return (!doesnt_cross);
}

void Beam3d::plot(rave::EnvironmentBasePtr env, Vector3d color) {
//	Vector3d color(0,1,0);
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

/**
 * Finds the intersection of line segment (p0, p1) with the triangle
 * Returns false if no intersection
 */
bool Triangle3d::intersection_with_segment(const Vector3d& p0, const Vector3d& p1, Vector3d& intersection) {
	// halfplane: <n, x - origin> = 0
	// segment: x = t*(p1 - p0) + p0 | 0 <= t <= 1
	// t = -<n, p0 - origin> / <n, p1 - p0>

	Vector3d n = (b-a).cross(c-a);
	Vector3d origin = a; // arbitrary
	double t = -n.dot(p0 - origin) / n.dot(p1 - p0);

	if ((t < 0) || (t > 1)) {
		return false; // segment doesn't cross plane
	}

	intersection = t*(p1 - p0) + p0;

	if (distance_to(intersection) > 1e-3) {
		return false; // segment crosses plane, but not the triangle
	}

	return true;
}

double Triangle3d::area() {
	Matrix3d rotation = rotation_to_align_with({0,0,1});
	return triangle_area((rotation*a).segment<2>(0),
			(rotation*b).segment<2>(0),
			(rotation*c).segment<2>(0));
}

void Triangle3d::plot(rave::EnvironmentBasePtr env, Vector3d color) {
//	Vector3d color(0,0,1);
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

/**
 * MeshUnit Constructors
 */

MeshUnit::MeshUnit(const Vector3d& a, const Vector3d& b, const Vector3d& c, const Vector3d& d) : visited(false) {
	tri0 = new Triangle3d(a, b, c);
	tri1 = new Triangle3d(a, c, d);
}

/**
 * MeshUnit Public methods
 */

bool MeshUnit::intersection_with_segment(const Vector3d& p0, const Vector3d& p1, Vector3d& intersection) {
	if (tri0->intersection_with_segment(p0, p1, intersection)) {
		return true;
	}

	if (tri1->intersection_with_segment(p0, p1, intersection)) {
		return true;
	}

	return false;
}

/**
 * Mesh Constructors
 */

Mesh::Mesh(int r, int c) : rows(r), cols(c) {
	grid = std::vector<MeshUnit*>(r*c);
}

/**
 * Mesh Public methods
 */

void Mesh::connect() {
	// connect left/right column
	for(int i=0; i < rows; ++i) {
		// add above
		if (i > 0) {
			get(i,0)->add_neighbor(get(i-1,0));
			get(i,cols-1)->add_neighbor(get(i-1,cols-1));
		}
		// add to left/right
		get(i,0)->add_neighbor(get(i,1));
		get(i,cols-1)->add_neighbor(get(i,cols-2));
		// add bottom
		if (i < rows-1) {
			get(i,0)->add_neighbor(get(i+1,0));
			get(i,cols-1)->add_neighbor(get(i+1,cols-1));
		}
	}

	// connect top/bottom row
	for(int j=0; j < cols; ++j) {
		// add left
		if (j > 0) {
			get(0,j)->add_neighbor(get(0,j-1));
			get(rows-1,j)->add_neighbor(get(rows-1,j-1));
		}
		// add top/bottom
		get(0,j)->add_neighbor(get(1,j));
		get(rows-1,j)->add_neighbor(get(rows-2,j));
		// add right
		if (j < cols-1) {
			get(0,j)->add_neighbor(get(0,j+1));
			get(rows-1,j)->add_neighbor(get(rows-1,j+1));
		}
	}

	std::vector<std::vector<int> > offsets = { {0, 1}, {-1, 1}, {-1, 0}, {-1, -1}, {0, -1}, {1, -1}, {1, 0}, {1, 1} };
	// connect inside
	for(int i=1; i < rows-1; ++i) {
		for(int j=1; j < cols-1; ++j) {
			for(int k=0; k < offsets.size(); ++k) {
				get(i,j)->add_neighbor(get(i+offsets[k][0],j+offsets[k][1]));
			}
		}
	}
}

void Mesh::plot(rave::EnvironmentBasePtr env) {
	Vector3d color(0,0,1);
	for(int i=0; i < rows; ++i) {
		for(int j=0; j < cols; ++j) {
			get(i,j)->plot(env, color);
		}
	}
}

bool Mesh::find_intersection_starting_from(const Vector3d& p0, const Vector3d& p1, MeshUnit* start,
			Vector3d& intersection, MeshUnit** end) {
	for(int i=0; i < grid.size(); ++i) {
		grid[i]->set_visited(false);
	}

	std::queue<MeshUnit*> q;
	start->set_visited(true);
	q.push(start);
	MeshUnit* curr;
	int n = 0;
	while(q.size() > 0) {
		curr = q.front();
		q.pop();
		++n;

		if (curr->intersection_with_segment(p0, p1, intersection)) {
			*end = curr;
//			std::cin.ignore();
			if (n > 30) {
				std::cout << "num_visited: " << n << "\n";
			}
			return true;
		}

		std::vector<MeshUnit*> neighbors = curr->get_neighbors();
		for(int i=0; i < neighbors.size(); ++i) {
			if (!neighbors[i]->is_visited()) {
				neighbors[i]->set_visited(true);
				q.push(neighbors[i]);
			}
		}
	}

	*end = curr;
	return false;
}

int Mesh::num_not_visited() {
	int n = 0;
	for(int i=0; i < grid.size(); ++i) {
		if (!grid[i]->is_visited()) {
			n++;
		}
	}
	return n;
}
