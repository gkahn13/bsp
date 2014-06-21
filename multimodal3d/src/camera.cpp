#include "../include/camera.h"

/**
 * Camera constructors
 */

Camera::Camera(rave::RobotBasePtr r, rave::SensorBasePtr s, double mr) : robot(r), sensor(s), max_range(mr) {
	rave::SensorBase::SensorGeometryPtr geom = sensor->GetSensorGeometry(rave::SensorBase::ST_Camera);

	boost::shared_ptr<rave::SensorBase::CameraGeomData> cam_geom =
			boost::static_pointer_cast<rave::SensorBase::CameraGeomData>(geom);

	height = cam_geom->height;
	width = cam_geom->width;
	f = cam_geom->KK.fx;
	F = cam_geom->KK.focal_length;

	H = F*(HEIGHT/f);
	W = F*(WIDTH/f);
}

/**
 * Camera public methods
 */

Matrix<double,N,3> Camera::get_directions() {
	Matrix<double,H_SUB,W_SUB> height_grid = Matrix<double,H_SUB,1>::LinSpaced(H_SUB, -H/2.0, H/2.0).replicate(1,W_SUB);
	Matrix<double,H_SUB,W_SUB> width_grid = Matrix<double,1,W_SUB>::LinSpaced(W_SUB, -W/2.0, W/2.0).replicate(H_SUB,1);

	Matrix<double,N,1> height_grid_vec(height_grid.data());
	Matrix<double,N,1> width_grid_vec(width_grid.data());
	Matrix<double,N,1> z_grid = Matrix<double,N,1>::Zero();

	Matrix<double,N,3> offsets;
	offsets << width_grid_vec, height_grid_vec, z_grid;

	Matrix<double,N,3> points_cam = RowVector3d(0,0,max_range).replicate(N,1) + (max_range/F)*offsets;

	Matrix4d ref_from_world = rave_utils::rave_to_eigen(sensor->GetTransform());
	Vector3d origin_world_pos = ref_from_world.block<3,1>(0,3);

	Matrix<double,N,3> directions;

	Matrix4d point_cam = Matrix4d::Identity();
	Vector3d point_world;
	for(int i=0; i < N; ++i) {
		point_cam.block<3,1>(0,3) = points_cam.row(i);
		point_world = (ref_from_world*point_cam).block<3,1>(0,3);

		directions.row(i) = point_world - origin_world_pos;

//		Vector3d color(1,0,0);
//		rave_utils::plot_segment(robot->GetEnv(), origin_world_pos, origin_world_pos + directions.row(i).transpose(), color);
	}

	return directions;
}

std::vector<std::vector<Beam3d> > Camera::get_beams() {
	std::vector<std::vector<Beam3d> > beams(H_SUB-1, std::vector<Beam3d>(W_SUB-1));

	RowVector3d origin_pos = rave_utils::rave_to_eigen(sensor->GetTransform().trans);

	util::Timer dirs_timer;
	util::Timer_tic(&dirs_timer);
	Matrix<double,N,3> dirs = get_directions();
	double dirs_time = util::Timer_toc(&dirs_timer);
	std::cout << "dirs_time: " << dirs_time << "\n";

	Matrix<double,N,3> hits;

	rave::EnvironmentBasePtr env = robot->GetEnv();
	rave::RAY ray;
	ray.pos = sensor->GetTransform().trans;
	rave::CollisionReportPtr report(new rave::CollisionReport());
	for(int i=0; i < N; ++i) {
		ray.dir.x = dirs(i,0);
		ray.dir.y = dirs(i,1);
		ray.dir.z = dirs(i,2);
		if (env->CheckCollision(ray, report)) {
			hits.row(i) = rave_utils::rave_to_eigen(report->contacts[0].pos);
		} else {
			hits.row(i) = origin_pos + (max_range / sqrt(ray.dir.lengthsqr2()))*dirs.row(i);
		}
	}

//	for(int i=0; i < N; ++i) {
//		rave_utils::plot_point(env, hits.row(i), {0,1,0}, .01);
//		std::cout << hits.row(i) << "\n";
////		std::cin.ignore();
//	}

	for(int j=0; j < W_SUB-1; ++j) {
		for(int i=0; i < H_SUB-1; ++i) {
			beams[i][j].base = origin_pos;
			beams[i][j].a = hits.row((j+1)*H_SUB+i);
			beams[i][j].b = hits.row(j*H_SUB+i);
			beams[i][j].c = hits.row(j*H_SUB+i+1);
			beams[i][j].d = hits.row((j+1)*H_SUB+i+1);
		}
	}

	return beams;
}
