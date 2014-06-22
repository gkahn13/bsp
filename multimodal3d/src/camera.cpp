//#include "../include/camera.h"
//
///**
// * Camera constructors
// */
//
//Camera::Camera(rave::RobotBasePtr r, rave::SensorBasePtr s, double mr) : robot(r), sensor(s), max_range(mr) {
//	rave::SensorBase::SensorGeometryPtr geom = sensor->GetSensorGeometry(rave::SensorBase::ST_Camera);
//
//	boost::shared_ptr<rave::SensorBase::CameraGeomData> cam_geom =
//			boost::static_pointer_cast<rave::SensorBase::CameraGeomData>(geom);
//
//	height = cam_geom->height;
//	width = cam_geom->width;
//	f = cam_geom->KK.fx;
//	F = cam_geom->KK.focal_length;
//
//	H = F*(HEIGHT/f);
//	W = F*(WIDTH/f);
//}
//
///**
// * Camera public methods
// */
//
//Matrix<double,N,3> Camera::get_directions() {
//	Matrix<double,H_SUB,W_SUB> height_grid = Matrix<double,H_SUB,1>::LinSpaced(H_SUB, -H/2.0, H/2.0).replicate(1,W_SUB);
//	Matrix<double,H_SUB,W_SUB> width_grid = Matrix<double,1,W_SUB>::LinSpaced(W_SUB, -W/2.0, W/2.0).replicate(H_SUB,1);
//
//	Matrix<double,N,1> height_grid_vec(height_grid.data());
//	Matrix<double,N,1> width_grid_vec(width_grid.data());
//	Matrix<double,N,1> z_grid = Matrix<double,N,1>::Zero();
//
//	Matrix<double,N,3> offsets;
//	offsets << width_grid_vec, height_grid_vec, z_grid;
//
//	Matrix<double,N,3> points_cam = RowVector3d(0,0,max_range).replicate(N,1) + (max_range/F)*offsets;
//
//	Matrix4d ref_from_world = rave_utils::rave_to_eigen(sensor->GetTransform());
//	Vector3d origin_world_pos = ref_from_world.block<3,1>(0,3);
//
//	Matrix<double,N,3> directions;
//
//	Matrix4d point_cam = Matrix4d::Identity();
//	Vector3d point_world;
//	for(int i=0; i < N; ++i) {
//		point_cam.block<3,1>(0,3) = points_cam.row(i);
//		point_world = (ref_from_world*point_cam).block<3,1>(0,3);
//
//		directions.row(i) = point_world - origin_world_pos;
//	}
//
//	return directions;
//}
//
//std::vector<std::vector<Beam3d> > Camera::get_beams() {
//	std::vector<std::vector<Beam3d> > beams(H_SUB-1, std::vector<Beam3d>(W_SUB-1));
//
//	RowVector3d origin_pos = rave_utils::rave_to_eigen(sensor->GetTransform().trans);
//
//	Matrix<double,N,3> dirs = get_directions();
//
//	Matrix<double,N,3> hits;
//
//	rave::EnvironmentBasePtr env = robot->GetEnv();
//	rave::RAY ray;
//	ray.pos = sensor->GetTransform().trans;
//	rave::CollisionReportPtr report(new rave::CollisionReport());
//	for(int i=0; i < N; ++i) {
//		ray.dir.x = dirs(i,0);
//		ray.dir.y = dirs(i,1);
//		ray.dir.z = dirs(i,2);
//		if (env->CheckCollision(ray, report)) {
//			hits.row(i) = rave_utils::rave_to_eigen(report->contacts[0].pos);
//		} else {
//			hits.row(i) = origin_pos + (max_range / sqrt(ray.dir.lengthsqr2()))*dirs.row(i);
//		}
//	}
//
//	for(int j=0; j < W_SUB-1; ++j) {
//		for(int i=0; i < H_SUB-1; ++i) {
//			beams[i][j].base = origin_pos;
//			beams[i][j].a = hits.row((j+1)*H_SUB+i);
//			beams[i][j].b = hits.row(j*H_SUB+i);
//			beams[i][j].c = hits.row(j*H_SUB+i+1);
//			beams[i][j].d = hits.row((j+1)*H_SUB+i+1);
//		}
//	}
//
//	return beams;
//}
//
//std::vector<Triangle3d> Camera::get_border(const std::vector<std::vector<Beam3d> >& beams) {
//	std::vector<Triangle3d> border;
//
//	int rows = beams.size(), cols = beams[0].size();
//
//	// deal with left and right border columns
//	for(int i=0; i < rows; ++i) {
//		border.push_back(Triangle3d(beams[i][0].base, beams[i][0].b, beams[i][0].c));
//		border.push_back(Triangle3d(beams[i][cols-1].base, beams[i][cols-1].a, beams[i][cols-1].d));
//	}
//
//	// deal with top and bottom border rows
//	for(int j=0; j < cols; ++j) {
//		border.push_back(Triangle3d(beams[0][j].base, beams[0][j].a, beams[0][j].b));
//		border.push_back(Triangle3d(beams[rows-1][j].base, beams[rows-1][j].c, beams[rows-1][j].d));
//	}
//
//	// connect with left and top
//	for(int i=0; i < rows; ++i) {
//		for(int j=0; j < cols; ++j) {
//			// left
//			if (i > 0) {
//				border.push_back(Triangle3d(beams[i-1][j].a, beams[i-1][j].d, beams[i][j].b));
//				border.push_back(Triangle3d(beams[i-1][j].b, beams[i][j].b, beams[i][j].c));
//			}
//
//			if (j > 0) {
//				border.push_back(Triangle3d(beams[i][j-1].c, beams[i][j-1].d, beams[i][j].b));
//				border.push_back(Triangle3d(beams[i][j-1].b, beams[i][j].b, beams[i][j].a));
//			}
//		}
//	}
//
//	std::vector<Triangle3d> pruned_border;
//	for(int i=0; i < border.size(); ++i) {
//		if (border[i].area() > epsilon) {
//			pruned_border.push_back(border[i]);
//		}
//	}
//
//	return pruned_border;
//}
//
//double Camera::signed_distance(const Vector3d& p, std::vector<std::vector<Beam3d> >& beams, std::vector<Triangle3d>& border) {
//	bool is_inside = false;
//	for(int i=0; i < beams.size(); ++i) {
//		for(int j=0; j < beams[i].size(); ++j) {
//			if (beams[i][j].is_inside(p)) {
//				is_inside = true;
//				break;
//			}
//		}
//		if (is_inside) { break; }
//	}
//
//	double sd_sign = (is_inside) ? -1 : 1;
//
//	double sd = INFINITY;
//	for(int i=0; i < border.size(); ++i) {
//		sd = std::min(sd, border[i].distance_to(p));
//	}
//
//	return (sd_sign*sd);
//}
