#include "voxel-grid.h"

/**
 * Local Helper classes/functions
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

VoxelGrid::VoxelGrid(const Vector3d& pos_center, double x, double y, double z, int res, const Matrix4d& init_cam_pose) : resolution(res) {
	bottom_corner = pos_center - Vector3d(x/2.0, y/2.0, z/2.0);
	top_corner = pos_center + Vector3d(x/2.0, y/2.0, z/2.0);
	size = Vector3d(x,y,z);
	dx = x / double(resolution);
	dy = y / double(resolution);
	dz = z / double(resolution);
	radius = std::min(dx, std::min(dy, dz))/10.0;

	TSDF = new Cube(resolution, resolution, resolution);
	TSDF->set_all(1);

	offsets.push_back(Vector3i(0, 0, -1));
	offsets.push_back(Vector3i(1, 0, -1));
	offsets.push_back(Vector3i(-1, 0, -1));
	offsets.push_back(Vector3i(0, 1, -1));
	offsets.push_back(Vector3i(0, -1, -1));
	offsets.push_back(Vector3i(1, 1, -1));
	offsets.push_back(Vector3i(1, -1, -1));
	offsets.push_back(Vector3i(-1, 1, -1));
	offsets.push_back(Vector3i(-1, -1, -1));

	offsets.push_back(Vector3i(1, 0, 0));
	offsets.push_back(Vector3i(-1, 0, 0));
	offsets.push_back(Vector3i(0, 1, 0));
	offsets.push_back(Vector3i(0, -1, 0));
	offsets.push_back(Vector3i(1, 1, 0));
	offsets.push_back(Vector3i(1, -1, 0));
	offsets.push_back(Vector3i(-1, 1, 0));
	offsets.push_back(Vector3i(-1, -1, 0));

	offsets.push_back(Vector3i(0, 0, 1));
	offsets.push_back(Vector3i(1, 0, 1));
	offsets.push_back(Vector3i(-1, 0, 1));
	offsets.push_back(Vector3i(0, 1, 1));
	offsets.push_back(Vector3i(0, -1, 1));
	offsets.push_back(Vector3i(1, 1, 1));
	offsets.push_back(Vector3i(1, -1, 1));
	offsets.push_back(Vector3i(-1, 1, 1));
	offsets.push_back(Vector3i(-1, -1, 1));

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

	gpu_size = Vector3d(y,z,x);
	gpu_resolution = Vector3i(pcl::device::kinfuLS::VOLUME_X,
			pcl::device::kinfuLS::VOLUME_Y,
			pcl::device::kinfuLS::VOLUME_Z);

	pcl_tsdf = pcl::gpu::kinfuLS::TsdfVolume::Ptr(new pcl::gpu::kinfuLS::TsdfVolume(gpu_resolution));
	pcl_tsdf->setSize(gpu_size.cast<float>());
	pcl_tsdf->setTsdfTruncDist(.03);

	gpu_pcl_tsdf_origin = Matrix4d::Identity();
	gpu_pcl_tsdf_origin.block<3,3>(0,0) << -1, 0, 0,
			0, 0, -1,
			0, -1, 0;
	gpu_pcl_tsdf_origin.block<3,1>(0,3) = pos_center + Vector3d(x/2., y/2., z/2.);


	// setup kinfu tracker
	float shifting_distance = 5.0f;
	pcl_kinfu_tracker = new pcl::gpu::kinfuLS::KinfuTracker(gpu_size.cast<float>(), shifting_distance, HEIGHT_FULL, WIDTH_FULL);
	pcl_kinfu_tracker->setDepthIntrinsics(intrinsics::fx, intrinsics::fy, intrinsics::cx, intrinsics::cy);

	Matrix4d init_cam_pose_gpu = gpu_pcl_tsdf_origin.inverse()*init_cam_pose;
	Vector3f t = (init_cam_pose_gpu.block<3,1>(0,3)).cast<float>();
	Matrix3f R = (init_cam_pose_gpu.block<3,3>(0,0)).cast<float>();

	Affine3f affine_init_cam_pose = Translation3f(t) * AngleAxisf(R);

	pcl_kinfu_tracker->setInitialCameraPose(affine_init_cam_pose);
	pcl_kinfu_tracker->reset();
//	pcl_kinfu_tracker->setDisableICP();
}

/**
 * VoxelGrid Public methods
 */

void VoxelGrid::update_kinfu(const Matrix<double,HEIGHT_FULL,WIDTH_FULL>& zbuffer) {
	typedef pcl::gpu::DeviceArray2D<unsigned short> Depth;
	Depth depth(H_SUB,W_SUB);

	std::vector<unsigned short> data(HEIGHT_FULL*WIDTH_FULL);
	int cols = WIDTH_FULL;

	int index = 0;
	for(int i=0; i < HEIGHT_FULL; ++i) {
		for(int j=0; j < WIDTH_FULL; ++j) {
			data[index++] = static_cast<unsigned short>(1e3*zbuffer(i,j)); // in mm
		}
	}

	depth.upload(data, cols);

	std::cout << "prev pose:\n" << pcl_kinfu_tracker->getCameraPose().translation().transpose() << "\n";
	std::cout << pcl_kinfu_tracker->getCameraPose().rotation() << "\n";

	(*pcl_kinfu_tracker)(depth);

	std::cout << "pose pose: " << pcl_kinfu_tracker->getCameraPose().translation().transpose() << "\n";
	std::cout << pcl_kinfu_tracker->getCameraPose().rotation() << "\n\n";


}

void VoxelGrid::update_TSDF(const StdVector3d& pcl) {
	Vector3i voxel;
	StdVector3i neighbors;
	std::vector<double> dists;
	for(int i=0; i < pcl.size(); ++i) {
		if (is_valid_point(pcl[i])) {
			voxel = voxel_from_point(pcl[i]);
			TSDF->set(voxel, 0);

//			// set neighbors also, fill in the gaps
//			get_voxel_neighbors_and_dists(voxel, neighbors, dists);
//			for(int i=0; i < neighbors.size(); ++i) {
//				TSDF->set(neighbors[i], 0);
//			}
		}
	}

//	// TODO: temp
//	for(int y=0; y < resolution; ++y) {
//		for(int z=0; z < resolution; ++z) {
//			TSDF->set(Vector3i(resolution/2,y,z),0);
//		}
//	}

	upload_to_pcl_tsdf();
}

Cube VoxelGrid::get_ODF(const Vector3d& obj) {
	typedef boost::heap::fibonacci_heap<VoxelDist>::handle_type handle_t;

	Cube ODF(resolution, resolution, resolution);
	ODF.set_all(INFINITY);

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

	StdVector3i neighbors;
	std::vector<double> neighbor_dists;
	while (pq.size() > 0) {
		VoxelDist curr = pq.top();
		pq.pop();

//		rave_utils::plot_point(env, point_from_voxel(curr.voxel), Vector3d(1,0,0));
//		std::cin.ignore();

		ODF.set(curr.voxel, curr.dist);

		get_voxel_neighbors_and_dists(curr.voxel, neighbors, neighbor_dists);
		for(int i=0; i < neighbors.size(); ++i) {
			double dist = curr.dist + neighbor_dists[i];
			if (dist < ODF.get(neighbors[i])) {
				ODF.set(neighbors[i], dist);
				std::vector<int> n_vec = {neighbors[i](0), neighbors[i](1), neighbors[i](2)};
				handle_t handle = voxel_handle[n_vec];
				pq.update(handle, VoxelDist(neighbors[i], dist));
			}
		}
	}

	return ODF;
}

Vector3d VoxelGrid::signed_distance_complete_voxel_center(const Vector3d& object, const Cube& ODF,
		Camera* cam, const Matrix<double,H_SUB,W_SUB>& zbuffer, const Matrix4d& cam_pose) {
	bool obj_in_fov = cam->is_in_fov(object, zbuffer, cam_pose);

	double min_dist = INFINITY;
	Vector3i min_voxel;
	for(int i=0; i < resolution; ++i) {
		for(int j=0; j < resolution; ++j) {
			for(int k=0; k < resolution; ++k) {
				if (TSDF->get(i,j,k) != 0) {
					Vector3d voxel_center = point_from_voxel(Vector3i(i,j,k));
					if (obj_in_fov) {
						if (!(cam->is_in_fov(voxel_center, zbuffer, cam_pose)) && (ODF.get(i,j,k) < min_dist)) {
							min_dist = ODF.get(i,j,k);
							min_voxel = {i,j,k};
						}
					} else {
						if ((cam->is_in_fov(voxel_center, zbuffer, cam_pose)) && (ODF.get(i,j,k) < min_dist)) {
							min_dist = ODF.get(i,j,k);
							min_voxel = {i,j,k};
						}
					}
				}
			}
		}
	}

	rave_utils::plot_point(cam->get_sensor()->GetEnv(), point_from_voxel(min_voxel), Vector3d(0,1,0), .03);

	return point_from_voxel(min_voxel);
}

inline bool XOR(const bool& lhs, const bool& rhs) {
    return !( lhs && rhs ) && ( lhs || rhs );
}

Vector3d VoxelGrid::signed_distance_greedy_voxel_center(const Vector3d& object, const Cube& ODF,
			Camera* cam, const Matrix<double,H_SUB,W_SUB>& zbuffer, const Matrix4d& cam_pose) {
	bool obj_in_fov = cam->is_in_fov(object, zbuffer, cam_pose);

	StdVector3i neighbors;
	std::vector<double> neighbor_dists;

	std::vector<Vector3i> start_voxels;

	if (obj_in_fov) {
		start_voxels.push_back(Vector3i(resolution-1, int(resolution/2), resolution-1));
		start_voxels.push_back(Vector3i(0, resolution-1, resolution-1));
		start_voxels.push_back(Vector3i(0, 0, resolution-1));
	} else {
		Matrix4d start_cam_frame = Matrix4d::Identity();
		start_cam_frame(2,3) = .2; // intrinsics::MIN_RANGE;
		Vector3d start_world = (cam_pose*start_cam_frame).block<3,1>(0,3);

		start_voxels.push_back(voxel_from_point(start_world));
	}

	Vector3i min_voxel = start_voxels[0];
	for(int s=0; s < start_voxels.size(); ++s) {
		Vector3i curr_voxel = start_voxels[s], next_voxel = Vector3i::Zero();
		double curr_OD = ODF.get(curr_voxel), next_OD;

		while(true) {
			rave_utils::plot_point(cam->get_sensor()->GetEnv(), point_from_voxel(curr_voxel), Vector3d(0,0,1), .01);

			next_OD = INFINITY;
			get_voxel_neighbors_and_dists(curr_voxel, neighbors, neighbor_dists);
			for(int i=0; i < neighbors.size(); ++i) {
				if (XOR(obj_in_fov, cam->is_in_fov(point_from_voxel(neighbors[i]), zbuffer, cam_pose))) {
					double neighbor_OD = ODF.get(neighbors[i]);
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

		if (curr_OD < ODF.get(min_voxel)) {
			min_voxel = curr_voxel;
		}
	}

	rave_utils::plot_point(cam->get_sensor()->GetEnv(), point_from_voxel(min_voxel), Vector3d(0,0,1), .02);

	return point_from_voxel(min_voxel);
}

double VoxelGrid::signed_distance_complete(const Vector3d& object, Cube& ODF,
		Camera* cam, const Matrix<double,H_SUB,W_SUB>& zbuffer, const Matrix4d& cam_pose) {

	Vector3d voxel_center = signed_distance_complete_voxel_center(object, ODF, cam, zbuffer, cam_pose);
	double dist = (object - voxel_center).norm();

	bool obj_in_fov = cam->is_in_fov(object, zbuffer, cam_pose);
	double sd = (obj_in_fov) ? dist : -dist;

	return sd;
}

double VoxelGrid::signed_distance_greedy(const Vector3d& object, const Cube& ODF,
		Camera* cam, const Matrix<double,H_SUB,W_SUB>& zbuffer, const Matrix4d& cam_pose) {

	Vector3d voxel_center = signed_distance_greedy_voxel_center(object, ODF, cam, zbuffer, cam_pose);
	double dist = (object - voxel_center).norm();

	bool obj_in_fov = cam->is_in_fov(object, zbuffer, cam_pose);
	double sd = (obj_in_fov) ? dist : -dist;

	return sd;
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

Matrix<double,H_SUB,W_SUB> VoxelGrid::get_zbuffer(const Matrix4d& cam_pose) {
	pcl::gpu::kinfuLS::RayCaster ray_caster(H_SUB, W_SUB,
			intrinsics::fx_sub, intrinsics::fy_sub, intrinsics::cx_sub, intrinsics::cy_sub);

	Matrix4d cam_pose_gpu = gpu_pcl_tsdf_origin.inverse()*cam_pose;
	Vector3f t = (cam_pose_gpu.block<3,1>(0,3)).cast<float>();
	Matrix3f R = (cam_pose_gpu.block<3,3>(0,0)).cast<float>();

//	Matrix3f R = Matrix3f::Identity();  // 3X3 identity matrix
//	Vector3f t = gpu_size.cast<float>() * 0.5 - Vector3f (0,0, gpu_size(2) /2 * 1.2);


	Affine3f affine_cam_pose = Translation3f(t) * AngleAxisf(R);

	pcl::gpu::kinfuLS::tsdf_buffer *tsdf_buffer = new pcl::gpu::kinfuLS::tsdf_buffer();
	tsdf_buffer->tsdf_memory_start = 0;
	tsdf_buffer->tsdf_memory_end = 0;
	tsdf_buffer->tsdf_rolling_buff_origin = 0;
	tsdf_buffer->origin_GRID.x = 0;
	tsdf_buffer->origin_GRID.y = 0;
	tsdf_buffer->origin_GRID.z = 0;
	tsdf_buffer->origin_GRID_global.x = 0.f;
	tsdf_buffer->origin_GRID_global.y = 0.f;
	tsdf_buffer->origin_GRID_global.z = 0.f;
	tsdf_buffer->origin_metric.x = 0.f;
	tsdf_buffer->origin_metric.y = 0.f;
	tsdf_buffer->origin_metric.z = 0.f;
	tsdf_buffer->volume_size.x = gpu_size(0);
	tsdf_buffer->volume_size.y = gpu_size(1);
	tsdf_buffer->volume_size.z = gpu_size(2);
	tsdf_buffer->voxels_size.x = gpu_resolution(0);
	tsdf_buffer->voxels_size.y = gpu_resolution(1);
	tsdf_buffer->voxels_size.z = gpu_resolution(2);

	ray_caster.run(*pcl_tsdf, affine_cam_pose, tsdf_buffer);

	typedef pcl::gpu::DeviceArray2D<unsigned short> Depth;
	Depth depth(H_SUB,W_SUB);
	ray_caster.generateDepthImage(depth);

	int c;
	std::vector<unsigned short> data;
	depth.download(data, c);

	std::cout << "depth size: " << data.size() << "\n";
	std::cout << "should be: " << H_SUB*W_SUB << "\n";

	std::cout << "distance between rows: " << c << "\n";

	Matrix<double,H_SUB,W_SUB> zbuffer;
	int index = 0;
	for(int j=0; j < W_SUB; ++j) {
		for(int i=0; i < H_SUB; ++i) {
			std::cout << data[index] << " ";
			zbuffer(i,j) = data[index++];
		}
		std::cout << "\n";
	}
	std::cout << "\n\n";

	std::cout << "tsdf_memory_start: " << tsdf_buffer->tsdf_memory_start << "\n";
	std::cout << "tsdf_memory_end: " << tsdf_buffer->tsdf_memory_end << "\n";

	return zbuffer;
}

void VoxelGrid::test_gpu_conversions(rave::EnvironmentBasePtr env) {
	rave_utils::plot_transform(env, rave_utils::eigen_to_rave(gpu_pcl_tsdf_origin));

	Affine3f pose_affine = pcl_kinfu_tracker->getCameraPose();
	Matrix4d pose = Matrix4d::Identity();
	pose.block<3,1>(0,3) = pose_affine.translation().cast<double>();
	pose.block<3,3>(0,0) = pose_affine.rotation().cast<double>();
	Matrix4d pose_world = gpu_pcl_tsdf_origin*pose;
	rave_utils::plot_transform(env, rave_utils::eigen_to_rave(pose_world));

	for(int i=0; i < gpu_resolution(0); i+=5) {
		rave_utils::plot_point(env, point_from_gpu_voxel(Vector3i(i,0,0)), Vector3d(1,0,0));
	}

	for(int j=0; j < gpu_resolution(1); j+=5) {
		rave_utils::plot_point(env, point_from_gpu_voxel(Vector3i(0,j,0)), Vector3d(0,1,0));
	}

	for(int k=0; k < gpu_resolution(2); k+=5) {
		rave_utils::plot_point(env, point_from_gpu_voxel(Vector3i(0,0,k)), Vector3d(0,0,1));
	}

	std::vector<float> tsdf_vector;
//	pcl_tsdf->downloadTsdf(tsdf_vector);
	pcl_kinfu_tracker->volume().downloadTsdf(tsdf_vector);

	int index = 0;
	for(int i=0; i < gpu_resolution(0); i++) {
		for(int j=0; j < gpu_resolution(1); j++) {
			for(int k=0; k < gpu_resolution(2); k++) {
				float val = tsdf_vector[index++];
				if ((i % 15 == 0) && (j % 15 == 0) && (k % 15 == 0)) {
//					std::cout << val << "\n";
					if (val == 1) {
						rave_utils::plot_point(env, point_from_gpu_voxel(Vector3i(i,j,k)), Vector3d(1,0,0));
					}
				}
			}
		}
	}
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

void VoxelGrid::plot_ODF(Cube& ODF, rave::EnvironmentBasePtr env) {
	double max_dist = -INFINITY;
	for(int i=0; i < resolution; ++i) {
		for(int j=0; j < resolution; ++j) {
			for(int k=0; k < resolution; ++k) {
				double dist = ODF.get(i,j,k);
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
				double dist = ODF.get(i,j,k);
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
	int step = 2;

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

void VoxelGrid::plot_FOV_full(rave::EnvironmentBasePtr env, Camera* cam, const Matrix<double,HEIGHT_FULL,WIDTH_FULL>& zbuffer, const Matrix4d& cam_pose) {
	int step = 2;

	for(int i=0; i < resolution; i+=step) {
		for(int j=0; j < resolution; j+=step) {
			for(int k=0; k < resolution; k+=step) {
				if (TSDF->get(i,j,k) != 0) {
					Vector3d voxel_center = point_from_voxel(Vector3i(i,j,k));
					if (cam->is_in_fov_full(voxel_center, zbuffer, cam_pose)) {
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

void VoxelGrid::upload_to_pcl_tsdf() {
	pcl::device::DeviceArray2D<int> volume = pcl_tsdf->data();
	size_t step = volume.step();
	int rows = volume.rows();
	int cols = volume.cols();

	// Needed, otherwise vector 'tsdf' would be altered.
	std::vector<float> tsdf_tmp(gpu_resolution.prod());
	int x_scale = gpu_resolution(0) / resolution;
	int y_scale = gpu_resolution(1) / resolution;
	int z_scale = gpu_resolution(2) / resolution;

	std::cout << "x_scale: " << x_scale << "\n";
	std::cout << "y_scale: " << y_scale << "\n";
	std::cout << "z_scale: " << z_scale << "\n";

	int index = 0;
	for(int y=0; y < resolution; ++y) {
		for(int ys=0; ys < y_scale; ++ys) {
			for(int z=0; z < resolution; ++z) {
				for(int zs=0; zs < z_scale; ++zs) {
					for(int x=0; x < resolution; ++x) {
						for(int xs=0; xs < x_scale; ++xs) {
							short2 elem;

							elem.x = (short)(TSDF->get(x,y,z)*pcl::device::kinfuLS::DIVISOR);
							elem.y = (short)(1); // TODO: what should weight be?

							tsdf_tmp[index++] = *reinterpret_cast<float*>(&elem);
						}
					}
				}
			}
		}
	}

	pcl_tsdf->data().upload(&tsdf_tmp[0], step, rows, cols);
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

Vector3d VoxelGrid::point_from_gpu_voxel(const Vector3i& voxel) {
	double dx_gpu, dy_gpu, dz_gpu;
	dx_gpu = gpu_size(0) / gpu_resolution(0);
	dy_gpu = gpu_size(1) / gpu_resolution(1);
	dz_gpu = gpu_size(2) / gpu_resolution(2);

	Matrix4d gpu_pose = Matrix4d::Identity();
//	gpu_pose.block<3,1>(0,3) = Vector3d(voxel(0)*dx_gpu, voxel(1)*dy_gpu, voxel(2)*dz_gpu);
	gpu_pose.block<3,1>(0,3) = Vector3d(voxel(2)*dz_gpu, voxel(1)*dy_gpu, voxel(0)*dx_gpu); // swap x and z

	Vector3d world_center = (gpu_pcl_tsdf_origin*gpu_pose).block<3,1>(0,3);

	return world_center;
}
