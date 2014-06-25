#ifndef __GEOMETRY3D_H__
#define __GEOMETRY3D_H__

#include "../utils/rave-utils.h"
#include "geometry2d.h"

#include <math.h>
#include <queue>

#include <Eigen/Eigen>
using namespace Eigen;

#include <assert.h>

class Beam3d {
public:
	/**
	 * Order goes
	 * b ---- a
	 * |      |
	 * |      |
	 * c ---- d
	 */
	Vector3d base, a, b, c, d;

	Beam3d() { };
	Beam3d(const Vector3d& base_tmp, const Vector3d& a_tmp, const Vector3d& b_tmp, const Vector3d& c_tmp, const Vector3d& d_tmp) :
		base(base_tmp), a(a_tmp), b(b_tmp), c(c_tmp), d(d_tmp) { };

	bool is_inside(const Vector3d& p);
	bool is_crossed_by(const Vector3d& p0, const Vector3d& p1);

	void plot(rave::EnvironmentBasePtr env, Vector3d color);
};

class Triangle3d {
public:
	Vector3d a, b, c;

	Triangle3d() { };
	Triangle3d(const Vector3d& a_tmp, const Vector3d& b_tmp, const Vector3d& c_tmp) :
		a(a_tmp), b(b_tmp), c(c_tmp) { };

	Vector3d closest_point_to(const Vector3d& p);
	inline double distance_to(const Vector3d& p) { return (closest_point_to(p) - p).norm(); }

	bool intersection_with_segment(const Vector3d& p0, const Vector3d& p1, Vector3d& intersection);

	double area();

	void plot(rave::EnvironmentBasePtr env, Vector3d color);

private:
	Matrix3d rotation_to_align_with(const Vector3d& target);
};

class MeshUnit {
public:
	MeshUnit(const Vector3d& a, const Vector3d& b, const Vector3d& c, const Vector3d& d);

	bool intersection_with_segment(const Vector3d& p0, const Vector3d& p1, Vector3d& intersection);

	void add_neighbor(MeshUnit* neighbor) { neighbors.push_back(neighbor); }
	std::vector<MeshUnit*> get_neighbors() { return neighbors; }

	bool is_visited() { return visited; }
	void set_visited(bool v) { visited = v; }

	void plot(rave::EnvironmentBasePtr env, Vector3d color) { tri0->plot(env, color); tri1->plot(env, color); }

private:
	Triangle3d *tri0, *tri1;
	std::vector<MeshUnit*> neighbors;
	bool visited;
};

class Mesh {
public:
	Mesh(int r, int c);

	inline MeshUnit* get(int i, int j) { return grid[j*rows + i]; }
	inline void set(int i, int j, MeshUnit* m) { grid[j*rows + i] = m; }

	int num_rows() { return rows; }
	int num_cols() { return cols; }

	void connect();
	bool find_intersection_starting_from(const Vector3d& p0, const Vector3d& p1, MeshUnit* start,
			Vector3d& intersection, MeshUnit** end);
	void plot(rave::EnvironmentBasePtr env);

private:
	int rows, cols;
	std::vector<MeshUnit*> grid; // column major i.e. j*rows + i

	int num_visited() { return (rows*cols - num_not_visited()); }
	int num_not_visited();
};

#endif
