#include "voxel-grid.h"

/**
 * Local Helper classes
 */

struct VoxelDist {
	Vector3i voxel;
	double dist;

	VoxelDist(const Vector3i& v, const double d) : voxel(v), dist(d) { };

	bool operator<(VoxelDist const & rhs) const {
		return (dist > rhs.dist);
	}
};

/**
 * VoxelGrid Constructors
 */

VoxelGrid::VoxelGrid(const Vector3d& pos_center, double x, double y, double z, int res) : resolution(res) {
	bottom_corner = pos_center - Vector3d(x/2.0, y/2.0, z/2.0);
	top_corner = pos_center + Vector3d(x/2.0, y/2.0, z/2.0);
	dx = x / double(resolution);
	dy = y / double(resolution);
	dz = z / double(resolution);
	radius = std::min(dx, std::min(dy, dz))/10.0;

	TSDF = new Cube(resolution, resolution, resolution);
	TSDF->set_all(1);

	object = INFINITY*Vector3d::Ones();
	ODF = new Cube(resolution, resolution, resolution);
	ODF->set_all(-INFINITY);

	smallest_OD_neighbors = new CubeVector3i(resolution, resolution, resolution);

	offsets = {Vector3i(0, 0, -1),
			Vector3i(1, 0, -1),
			Vector3i(-1, 0, -1),
			Vector3i(0, 1, -1),
			Vector3i(0, -1, -1),
			Vector3i(1, 1, -1),
			Vector3i(1, -1, -1),
			Vector3i(-1, 1, -1),
			Vector3i(-1, -1, -1),

			Vector3i(1, 0, 0),
			Vector3i(-1, 0, 0),
			Vector3i(0, 1, 0),
			Vector3i(0, -1, 0),
			Vector3i(1, 1, 0),
			Vector3i(1, -1, 0),
			Vector3i(-1, 1, 0),
			Vector3i(-1, -1, 0),

			Vector3i(0, 0, 1),
			Vector3i(1, 0, 1),
			Vector3i(-1, 0, 1),
			Vector3i(0, 1, 1),
			Vector3i(0, -1, 1),
			Vector3i(1, 1, 1),
			Vector3i(1, -1, 1),
			Vector3i(-1, 1, 1),
			Vector3i(-1, -1, 1)};

	offset_dists = {dz,
                   Vector3d(dx,0,dz).norm(),
                   Vector3d(dx,0,dz).norm(),
                   Vector3d(0,dy,dz).norm(),
                   Vector3d(0,dy,dz).norm(),
                   Vector3d(dx,dy,dz).norm(),
                   Vector3d(dx,dy,dz).norm(),
                   Vector3d(dx,dy,dz).norm(),
                   Vector3d(dx,dy,dz).norm(),

                   dx,
                   dx,
                   dy,
                   dy,
                   Vector3d(dx,dy,0).norm(),
                   Vector3d(dx,dy,0).norm(),
                   Vector3d(dx,dy,0).norm(),
                   Vector3d(dx,dy,0).norm(),

                   dz,
                   Vector3d(dx,0,dz).norm(),
                   Vector3d(dx,0,dz).norm(),
                   Vector3d(0,dy,dz).norm(),
                   Vector3d(0,dy,dz).norm(),
                   Vector3d(dx,dy,dz).norm(),
                   Vector3d(dx,dy,dz).norm(),
                   Vector3d(dx,dy,dz).norm(),
                   Vector3d(dx,dy,dz).norm()};

}

/**
 * VoxelGrid Public methods
 */

void VoxelGrid::update_TSDF(const StdVector3d& pcl) {
	Vector3i voxel;
	StdVector3i neighbors;
	std::vector<double> dists;
	for(int i=0; i < pcl.size(); ++i) {
		if (is_valid_point(pcl[i])) {
			voxel = voxel_from_point(pcl[i]);
			TSDF->set(voxel, 0);

			// set neighbors also, fill in the gaps
			get_voxel_neighbors_and_dists(voxel, neighbors, dists);
			for(int i=0; i < neighbors.size(); ++i) {
				TSDF->set(neighbors[i], 0);
			}
		}
	}
}


void VoxelGrid::update_ODF(const Vector3d& obj, rave::EnvironmentBasePtr env) {
	typedef boost::heap::fibonacci_heap<VoxelDist>::handle_type handle_t;
	object = obj;
	ODF->set_all(INFINITY);

	Vector3i obj_voxel = voxel_from_point(obj);

	boost::heap::fibonacci_heap<VoxelDist> pq;
	std::map<std::vector<int>,handle_t> voxel_handle;

	for(int i=0; i < resolution; ++i) {
		for(int j=0; j < resolution; ++j) {
			for(int k=0; k < resolution; ++k) {
				if (TSDF->get(i,j,k) != 0) {
					const Vector3i voxel(i,j,k);
					double dist = (voxel == obj_voxel) ? 0 : INFINITY;
					handle_t handle = pq.push(VoxelDist(voxel, dist));
					voxel_handle[{i,j,k}] = handle;
				}
			}
		}
	}

//	VoxelDist* curr;
	StdVector3i neighbors;
	std::vector<double> neighbor_dists;
	while (pq.size() > 0) {
//		if (pq.size() % 1000 == 0) {
//			std::cout << "pq size: " << pq.size() << "\n";
//		}
		VoxelDist curr = pq.top();
		pq.pop();

//		rave_utils::plot_point(env, point_from_voxel(curr.voxel), Vector3d(1,0,0));
//		std::cin.ignore();

		ODF->set(curr.voxel, curr.dist);

		get_voxel_neighbors_and_dists(curr.voxel, neighbors, neighbor_dists);
		for(int i=0; i < neighbors.size(); ++i) {
			double dist = curr.dist + neighbor_dists[i];
			if (dist < ODF->get(neighbors[i])) {
				ODF->set(neighbors[i], dist);
				std::vector<int> n_vec = {neighbors[i](0), neighbors[i](1), neighbors[i](2)};
				handle_t handle = voxel_handle[n_vec];
				pq.update(handle, VoxelDist(neighbors[i], dist));
			}
		}
	}
}

double VoxelGrid::signed_distance_complete(Camera* cam, const Matrix<double,H_SUB,W_SUB>& zbuffer, const Matrix4d& cam_pose) {
	bool obj_in_fov = cam->is_in_fov(object, zbuffer, cam_pose);

	double min_dist = INFINITY;
	Vector3i min_voxel;
	for(int i=0; i < resolution; ++i) {
		for(int j=0; j < resolution; ++j) {
			for(int k=0; k < resolution; ++k) {
				if (TSDF->get(i,j,k) != 0) {
					Vector3d voxel_center = point_from_voxel(Vector3i(i,j,k));
					if (obj_in_fov) {
						if (!(cam->is_in_fov(voxel_center, zbuffer, cam_pose)) && (ODF->get(i,j,k) < min_dist)) {
							min_dist = ODF->get(i,j,k);
							min_voxel = {i,j,k};
						}
					} else {
						if ((cam->is_in_fov(voxel_center, zbuffer, cam_pose)) && (ODF->get(i,j,k) < min_dist)) {
							min_dist = ODF->get(i,j,k);
							min_voxel = {i,j,k};
						}
					}
				}
			}
		}
	}

	double sd = (obj_in_fov) ? min_dist : -min_dist;

	rave_utils::plot_point(cam->get_sensor()->GetEnv(), point_from_voxel(min_voxel), Vector3d(0,1,0), .03);

	return sd;
}

double VoxelGrid::signed_distance_greedy(Camera* cam, const Matrix<double,H_SUB,W_SUB>& zbuffer, const Matrix4d& cam_pose) {
	bool obj_in_fov = cam->is_in_fov(object, zbuffer, cam_pose);

//	double min_dist = INFINITY;
//	Vector3i min_voxel;

	StdVector3i neighbors;
	std::vector<double> neighbor_dists;
	if (obj_in_fov) {
		Vector3i curr_voxel(resolution-1, int(resolution/2), resolution-1), next_voxel = Vector3i::Zero();
		double curr_OD = ODF->get(curr_voxel), next_OD;

		while(true) {
			rave_utils::plot_point(cam->get_sensor()->GetEnv(), point_from_voxel(curr_voxel), Vector3d(0,0,1), .01);

			next_OD = INFINITY;
			get_voxel_neighbors_and_dists(curr_voxel, neighbors, neighbor_dists);
			for(int i=0; i < neighbors.size(); ++i) {
				if (!(cam->is_in_fov(point_from_voxel(neighbors[i]), zbuffer, cam_pose))) {
					double neighbor_OD = ODF->get(neighbors[i]);
					if (neighbor_OD < next_OD) {
						next_OD = neighbor_OD;
						next_voxel = neighbors[i];
					}
				}
			}

			if (next_OD >= curr_OD) {
				break;
			}

			curr_OD = next_OD;
			curr_voxel = next_voxel;
		}

		double sd = curr_OD;

		rave_utils::plot_point(cam->get_sensor()->GetEnv(), point_from_voxel(curr_voxel), Vector3d(0,0,1), .02);

		return sd;
	} else {
		Matrix4d start_cam_frame = Matrix4d::Identity();
		start_cam_frame(2,3) = MIN_RANGE; // zbuffer(H_SUB/2,W_SUB/2);
		Vector3d start_world = (cam_pose*start_cam_frame).block<3,1>(0,3);

//		int i, j;
//		for(i=0; i < H_SUB; ++i) {
//			for(j=0; j < W_SUB; ++j) {
////				double d = zbuffer(i,j);
////				if (cam->get_point_from_pixel_and_dist(Vector2i(i,j), d, cam_pose)) {
////
////				}
//				if ((zbuffer(i,j) > MIN_RANGE) && (zbuffer(i,j) < MAX_RANGE)) {
//					break;
//				}
//			}
//			if ((zbuffer(i,j) > MIN_RANGE) && (zbuffer(i,j) < MAX_RANGE)) {
//				break;
//			}
//		}
//		std::cout << "(i,j): (" << i << "," << j << ")\n";
//		std::cout << "zbuffer val: " << zbuffer(i,j) << "\n";
//
//		Vector2i pixel(i,j);
//		double dist = zbuffer(pixel(0), pixel(1));
//		Vector3d start_world = cam->get_point_from_pixel_and_dist(pixel, dist, cam_pose);

		rave_utils::plot_point(cam->get_sensor()->GetEnv(), start_world, Vector3d(1,1,0), .02);

		Vector3i curr_voxel = voxel_from_point(start_world), next_voxel = Vector3i::Zero();
		double curr_OD = ODF->get(curr_voxel), next_OD; // OD := object distance

		while(true) {
			rave_utils::plot_point(cam->get_sensor()->GetEnv(), point_from_voxel(curr_voxel), Vector3d(0,0,1), .01);

			next_OD = INFINITY;
			get_voxel_neighbors_and_dists(curr_voxel, neighbors, neighbor_dists);
			for(int i=0; i < neighbors.size(); ++i) {
				if (cam->is_in_fov(point_from_voxel(neighbors[i]), zbuffer, cam_pose)) {
					double neighbor_OD = ODF->get(neighbors[i]);
					if (neighbor_OD < next_OD) {
						next_OD = neighbor_OD;
						next_voxel = neighbors[i];
					}
				}
			}

			if (next_OD >= curr_OD) {
//				std::cout << "checking neighbors of neighbors\n";
//				std::cin.ignore();
//				// check neighbors of neighbors
//				StdVector3i n_ns;
//				std::vector<double> n_ns_dists;
//				for(int i=0; i < neighbors.size(); ++i) {
//					get_voxel_neighbors_and_dists(neighbors[i], n_ns, n_ns_dists);
//					for(int j=0; j < n_ns.size(); ++j) {
//						if (cam->is_in_fov(point_from_voxel(n_ns[i]), zbuffer, cam_pose)) {
//							double n_n_OD = ODF->get(n_ns[i]);
//							if (n_n_OD < next_OD) {
//								next_OD = n_n_OD;
//								next_voxel = n_ns[i];
//							}
//						}
//					}
//				}
//
//				if (next_OD >= curr_OD) {
//					break;
//				}

				break;
			}

			curr_OD = next_OD;
			curr_voxel = next_voxel;
		}

		double sd = -curr_OD;

		rave_utils::plot_point(cam->get_sensor()->GetEnv(), point_from_voxel(curr_voxel), Vector3d(0,0,1), .02);

		return sd;
	}

	return 0;
}

StdVector3d VoxelGrid::get_obstacles() {
	StdVector3d obstacles;
	for(int i=0; i < resolution; i++) {
		for(int j=0; j < resolution; j++) {
			for(int k=0; k < resolution; k++) {
				if (TSDF->get(i,j,k) == 0) {
					obstacles.push_back(point_from_voxel(Vector3i(i,j,k)));
				}
			}
		}
	}
	return obstacles;
}

/**
 * VoxelGrid Display methods
 */

void VoxelGrid::plot_TSDF(rave::EnvironmentBasePtr env) {
	Vector3d color(1,0,0);
	for(int i=0; i < resolution; ++i) {
		for(int j=0; j < resolution; ++j) {
			for(int k=0; k < resolution; ++k) {
				if (TSDF->get(i,j,k) == 0) {
					rave_utils::plot_point(env, point_from_voxel(Vector3i(i,j,k)), color, radius);
				}
			}
		}
	}
}

void VoxelGrid::plot_ODF(rave::EnvironmentBasePtr env) {
	double max_dist = -INFINITY;
	for(int i=0; i < resolution; ++i) {
		for(int j=0; j < resolution; ++j) {
			for(int k=0; k < resolution; ++k) {
				double dist = ODF->get(i,j,k);
				if ((dist > max_dist) && (dist < INFINITY)) {
					max_dist = dist;
				}
			}
		}
	}

	int step = 5;
	double size = .01;

	for(int i=0; i < resolution; i+=step) {
		for(int j=0; j < resolution; j+=step) {
			for(int k=0; k < resolution; k+=step) {
				double dist = ODF->get(i,j,k);
				if ((dist >= 0) && (dist < INFINITY)) {
					double dist_pct = dist / max_dist;
					Vector3d rgb = utils::hsv_to_rgb(Vector3d((2/3.)*dist_pct,1,1));
					rave_utils::plot_point(env, point_from_voxel(Vector3i(i,j,k)), rgb, size);
				} else {
					rave_utils::plot_point(env, point_from_voxel(Vector3i(i,j,k)), Vector3d(0,0,0), size);
				}
			}
		}
	}
}

void VoxelGrid::plot_FOV(rave::EnvironmentBasePtr env, Camera* cam, const Matrix<double,H_SUB,W_SUB>& zbuffer, const Matrix4d& cam_pose) {
	int step = 5;

	for(int i=0; i < resolution; i+=step) {
		for(int j=0; j < resolution; j+=step) {
			for(int k=0; k < resolution; k+=step) {
				if (TSDF->get(i,j,k) != 0) {
					Vector3d voxel_center = point_from_voxel(Vector3i(i,j,k));
					if (cam->is_in_fov(voxel_center, zbuffer, cam_pose)) {
						rave_utils::plot_point(env, voxel_center, Vector3d(0,1,0), .005);
					}
				}
			}
		}
	}
}

/**
 * VoxelGrid Private methods
 */

/**
 * \brief Called at the end of update_ODF.
 * 		  For each voxel, smallest_OD_neighbor to the nearest voxel
 * 		  that has the smallest object distance
 */
void VoxelGrid::update_smallest_OD_neighbors() {
	StdVector3i neighbors;
	std::vector<double> dists;
	for(int i=0; i < resolution; ++i) {
		for(int j=0; j < resolution; ++j) {
			for(int k=0; k < resolution; ++k) {
				Vector3i curr_voxel(i,j,k), min_neighbor(i,j,k);
				double min_neighbor_dist = INFINITY;

				get_voxel_neighbors_and_dists(curr_voxel, neighbors, dists);
				for(int l=0; l < neighbors.size(); ++l) {
					if (ODF->get(neighbors[i]) < ODF->get(min_neighbor)) {
						min_neighbor = neighbors[i];
					}
				}

				smallest_OD_neighbors->set(curr_voxel, min_neighbor);
			}
		}
	}
}

void VoxelGrid::get_voxel_neighbors_and_dists(const Vector3i& voxel, StdVector3i& neighbors, std::vector<double>& dists) {
	neighbors.clear();
	dists.clear();

	Vector3i neighbor;
	for(int i=0; i < offsets.size(); ++i) {
		neighbor = voxel + offsets[i];
		if ((neighbor.minCoeff() >= 0) && (neighbor.maxCoeff() < resolution)) {
			if (TSDF->get(neighbor) != 0) {
				neighbors.push_back(neighbor);
				dists.push_back(offset_dists[i]);
			}
		}
	}
}

Vector3i VoxelGrid::voxel_from_point(const Vector3d& point) {
	Vector3d relative = point - bottom_corner;
	Vector3i voxel(int(relative(0)/dx), int(relative(1)/dy), int(relative(2)/dz));

	if ((voxel.minCoeff() >= 0) && (voxel.maxCoeff() < resolution)) {
		return voxel;
	}

	LOG_WARN("Point is outside of VoxelGrid");
	return Vector3i(-1,-1,-1);
}

Vector3d VoxelGrid::point_from_voxel(const Vector3i& voxel) {
	Vector3d center(voxel(0)*dx, voxel(1)*dy, voxel(2)*dz);
	return (bottom_corner + center);
}
