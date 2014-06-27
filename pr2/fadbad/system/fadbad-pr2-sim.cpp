#include "fadbad-pr2-sim.h"

/**
 * FadbadArm constructors
 */

FadbadArm::FadbadArm(rave::RobotBasePtr robot, Arm::ArmType arm_type) {
	origin = rave_utils::rave_to_eigen(robot->GetLink("torso_lift_link")->GetTransform()).cast<bdouble>();

	arm_joint_axes = { Vector3b(0,0,1),
			Vector3b(0,1,0),
			Vector3b(1,0,0),
			Vector3b(0,1,0),
			Vector3b(1,0,0),
			Vector3b(0,1,0),
			Vector3b(1,0,0) };

	arm_link_trans = { Vector3b(0, (arm_type == Arm::ArmType::left) ? 0.188 : -0.188, 0),
			Vector3b(0.1, 0, 0),
			Vector3b(0, 0, 0),
			Vector3b(0.4, 0, 0),
			Vector3b(0, 0, 0),
			Vector3b(.321, 0, 0),
			Vector3b(.18, 0, 0) };
}

/**
 * FadbadArm public methods
 */

/**
 * \param angle: radians
 * \param axis: vector
 */
inline Matrix3b axis_angle_to_matrix(const bdouble& angle, const Vector3b& axis) {
	bdouble sina = sin(angle);
	bdouble cosa = cos(angle);

	Matrix3b R = Vector3b(cosa, cosa, cosa).asDiagonal();
	Matrix3b axis_outer = axis*axis.transpose();
	for(int i=0; i < 3; ++i) {
		for(int j=0; j < 3; ++j) {
			axis_outer(i,j) *= (1-cosa);
		}
	}
	R += axis_outer;

	Vector3b axis_sina = axis;
	for(int i=0; i < 3; ++i) { axis_sina(i) *= sina; }
	Matrix3b mat_sina;
	mat_sina << 0, -axis_sina(2), axis_sina(1),
			axis_sina(2), 0, -axis_sina(0),
			-axis_sina(1), axis_sina(0), 0;
	R += mat_sina;

	return R;
}

Matrix4b FadbadArm::get_pose(const Matrix<bdouble,ARM_DIM,1>& j) {
	Matrix4b pose_mat = origin;

	Matrix4b R;
	for(int i=0; i < ARM_DIM; ++i) {
		R.block<3,3>(0,0) = axis_angle_to_matrix(j(i), arm_joint_axes[i]);
		R.block<3,1>(0,3) = arm_link_trans[i];
		pose_mat = pose_mat*R;
	}

	return pose_mat;
}

/**
 * FadbadCamera constructors
 */

FadbadCamera::FadbadCamera(FadbadArm* a, const Matrix4d& g_t_to_s) : arm(a) {
	gripper_tool_to_sensor = g_t_to_s.cast<bdouble>();


	Matrix3b P = Matrix3b::Zero();
	P(0,0) = fx_sub;
	P(1,1) = fy_sub;
	P(2,2) = 1;
	P(0,2) = cx_sub;
	P(1,2) = cy_sub;

	depth_map = new FadbadDepthMap(P);
}


/**
 * FadbadCamera public methods
 */


std::vector<std::vector<FadbadBeam3d> > FadbadCamera::get_beams(const Matrix<bdouble,ARM_DIM,1>& j, const StdVector3b& pcl) {
	std::vector<std::vector<FadbadBeam3d> > beams(H_SUB-1, std::vector<FadbadBeam3d>(W_SUB-1));

	depth_map->clear();
	for(int i=0; i < pcl.size(); ++i) {
		depth_map->add_point(pcl[i], get_pose(j));
	}
	Matrix<bdouble,H_SUB,W_SUB> z_buffer = depth_map->get_z_buffer(get_position(j));

	RowVector3b origin_pos = get_position(j);

	Matrix<bdouble,N_SUB,3> dirs = get_directions(j, H_SUB, W_SUB, H_SUB_M, W_SUB_M);
	for(int i=0; i < N_SUB; ++i) {
		bdouble row_norm = dirs.row(i).norm();
		for(int j=0; j < 3; ++j) {
			dirs(i,j) /= row_norm;
		}
	}

	std::vector<std::vector<Vector3b> > hits(H_SUB, std::vector<Vector3b>(W_SUB));
	for(int i=0; i < H_SUB; ++i) {
		for(int j=0; j < W_SUB; ++j) {
//			hits[i][j] = origin_pos + z_buffer(i,j)*dirs.row(j*H_SUB+i);
			for(int k=0; k < 3; ++k) {
				hits[i][j](k) = origin_pos(k) + z_buffer(i,j)*dirs.row(j*H_SUB+i)(k);
			}
		}
	}

	for(int j=0; j < W_SUB-1; ++j) {
		for(int i=0; i < H_SUB-1; ++i) {
			beams[i][j].base = origin_pos;
			beams[i][j].a = hits[i][j+1];
			beams[i][j].b = hits[i][j];
			beams[i][j].c = hits[i+1][j];
			beams[i][j].d = hits[i+1][j+1];
		}
	}

	return beams;
}


std::vector<FadbadTriangle3d> FadbadCamera::get_border(const std::vector<std::vector<FadbadBeam3d> >& beams, bool with_side_border) {
	std::vector<FadbadTriangle3d> border;

	int rows = beams.size(), cols = beams[0].size();

	if (with_side_border) {
		// deal with left and right border columns
		for(int i=0; i < rows; ++i) {
			border.push_back(FadbadTriangle3d(beams[i][0].base, beams[i][0].b, beams[i][0].c));
			border.push_back(FadbadTriangle3d(beams[i][cols-1].base, beams[i][cols-1].a, beams[i][cols-1].d));
		}

		// deal with top and bottom border rows
		for(int j=0; j < cols; ++j) {
			border.push_back(FadbadTriangle3d(beams[0][j].base, beams[0][j].a, beams[0][j].b));
			border.push_back(FadbadTriangle3d(beams[rows-1][j].base, beams[rows-1][j].c, beams[rows-1][j].d));
		}
	}

	const bdouble min_sep = bepsilon;
	// connect with left and top and add center
	for(int i=0; i < rows; ++i) {
		for(int j=0; j < cols; ++j) {
			const FadbadBeam3d& curr = beams[i][j];
			// left
			if (i > 0) {
				const FadbadBeam3d& left = beams[i-1][j];
				if (((left.a - curr.b).norm() > min_sep) || ((left.d - curr.c).norm() > min_sep)) {
					// not touching, add linkage
					border.push_back(FadbadTriangle3d(left.a, left.d, curr.b));
					border.push_back(FadbadTriangle3d(left.b, curr.b, curr.c));
				}
			}

			// top
			if (j > 0) {
				const FadbadBeam3d& top = beams[i][j-1];
				if (((top.c - curr.b).norm() > min_sep) || ((top.d - curr.a).norm() > min_sep)) {
					// not touching, add linkage
					border.push_back(FadbadTriangle3d(top.c, top.d, curr.b));
					border.push_back(FadbadTriangle3d(top.b, curr.b, curr.a));
				}
			}

			border.push_back(FadbadTriangle3d(curr.a, curr.b, curr.c));
			border.push_back(FadbadTriangle3d(curr.a, curr.c, curr.d));
		}
	}

	std::vector<FadbadTriangle3d> pruned_border;
	for(int i=0; i < border.size(); ++i) {
		if (border[i].area() > bepsilon) {
			pruned_border.push_back(border[i]);
		}
	}

	return pruned_border;
}

bool FadbadCamera::is_inside(const Vector3b& p, std::vector<std::vector<FadbadBeam3d> >& beams) {
	bool inside = false;
	for(int i=0; i < beams.size(); ++i) {
		for(int j=0; j < beams[i].size(); ++j) {
			if (beams[i][j].is_inside(p)) {
				inside = true;
				break;
			}
		}
		if (inside) { break; }
	}

	return inside;
}

bdouble FadbadCamera::signed_distance(const Vector3b& p, std::vector<std::vector<FadbadBeam3d> >& beams, std::vector<FadbadTriangle3d>& border) {
	bdouble sd_sign = (is_inside(p, beams)) ? -1 : 1;

	bdouble sd = INFINITY;
	for(int i=0; i < border.size(); ++i) {
		bdouble dist = border[i].distance_to(p);
		sd = (dist < sd) ? dist : sd;
	}

	return (sd_sign*sd);
}


/**
 * FadbadCamera Private methods
 */


MatrixDynb FadbadCamera::get_directions(const Matrix<bdouble,ARM_DIM,1>& j, const int h, const int w, const bdouble h_meters, const bdouble w_meters) {
	const int n = h*w;
	MatrixDynb height_grid = VectorDynb::LinSpaced(h, -h_meters/2.0, h_meters/2.0).replicate(1,w);
	MatrixDynb width_grid = RowVectorDynb::LinSpaced(w, -w_meters/2.0, w_meters/2.0).replicate(h,1);

	MatrixDynb height_grid_vec(Map<VectorDynb>(height_grid.data(), n));
	MatrixDynb width_grid_vec(Map<VectorDynb>(width_grid.data(), n));
	VectorDynb z_grid = VectorDynb::Zero(n,1);

	MatrixDynb offsets(n,3);
	offsets << width_grid_vec, height_grid_vec, z_grid;

	for(int i=0; i < n; ++i) {
		for(int j=0; j < 3; ++j) {
			offsets *= (MAX_RANGE/FOCAL_LENGTH);
		}
	}

	MatrixDynb points_cam = RowVector3b(0,0,MAX_RANGE).replicate(n,1) + offsets;

	Matrix4b ref_from_world = get_pose(j);
	Vector3b origin_world_pos = ref_from_world.block<3,1>(0,3);

	MatrixDynb directions(n,3);

	Matrix4b point_cam = Matrix4b::Identity();
	Vector3b point_world;
	for(int i=0; i < n; ++i) {
		point_cam.block<3,1>(0,3) = points_cam.row(i);
		point_world = (ref_from_world*point_cam).block<3,1>(0,3);

		directions.row(i) = point_world - origin_world_pos;
	}

	return directions;
}


/**
 * FadbadDepthMap Constructors
 */

FadbadDepthMap::FadbadDepthMap(const Matrix3b& P_mat) : P(P_mat) {
	pixel_buckets = std::vector<std::vector<FadbadPixelBucket*> >(H_SUB, std::vector<FadbadPixelBucket*>(W_SUB));
	for(int i=0; i < H_SUB; ++i) {
		for(int j=0; j < W_SUB; ++j) {
			Vector2b center(i+0.5, j+0.5);
			pixel_buckets[i][j] = new FadbadPixelBucket(center);
		}
	}
};

/**
 * FadbadDepthMap Public methods
 */

// TODO: only add point to bucket if point is either (1) closer than points in bucket or (2) closet to points in bucket
void FadbadDepthMap::add_point(const Vector3b& point, const Matrix4b& cam_pose) {
	Matrix4b point_mat = Matrix4b::Identity();
	point_mat.block<3,1>(0,3) = point;

	Matrix4b point_mat_tilde = cam_pose.fullPivLu().solve(Matrix4b::Identity())*point_mat;

	Vector3b y = P*point_mat_tilde.block<3,1>(0,3);
	Vector2b pixel = {y(1)/y(2), y(0)/y(2)};

	if ((0 <= pixel(0)) && (pixel(0) < H_SUB) && (0 <= pixel(1)) && (pixel(1) < W_SUB) &&
			((cam_pose.block<3,1>(0,3) - point).norm() < MAX_RANGE)) { // TODO: should filter out points behind camera!
		for(int i=0; i < H_SUB; ++i) {
			for(int j=0; j < W_SUB; ++j) {
				bdouble i_b = i, j_b = j;
				if ((0 <= pixel(0) - i_b) && (pixel(0) - i_b < 1) &&
						(0 <= pixel(1) - j_b) && (pixel(1) - j_b < 1)) {
					pixel_buckets[i][j]->add_point(pixel, point);
				}
			}
		}
	}
}

Matrix<bdouble,H_SUB,W_SUB> FadbadDepthMap::get_z_buffer(const Vector3b& cam_pos) {
	Matrix<bdouble,H_SUB,W_SUB> z_buffer = Matrix<bdouble,H_SUB,W_SUB>::Ones();

	for(int i=0; i < H_SUB; ++i) {
		for(int j=0; j < W_SUB; ++j) {
			if (!(pixel_buckets[i][j]->is_empty())) {
				z_buffer(i,j) = (cam_pos - pixel_buckets[i][j]->average_point()).norm();
			} else if (num_neighbors_empty(i,j) >= 5 ) {
				z_buffer(i,j) = (cam_pos - average_of_neighbors(i, j)).norm();
			} else {
				z_buffer(i,j) = MAX_RANGE;
			}
		}
	}

	return z_buffer;
}

void FadbadDepthMap::clear() {
	for(int i=0; i < H_SUB; ++i) {
		for(int j=0; j < W_SUB; ++j) {
			pixel_buckets[i][j]->clear();
		}
	}
}

/**
 * FadbadDepthMap Private methods
 */

inline std::vector<std::vector<int> > offsets(int i, int j, int rows, int cols) {
	if ((i == 0) && (j == 0)) {
		return {{1,0}, {1,1}, {0,1}};
	} else if ((i == 0) && (j == cols-1)) {
		return {{0,-1}, {1,-1}, {1,0}};
	} else if ((i == rows-1) && (j == cols-1)) {
		return {{-1,0}, {-1,-1}, {0,-1}};
	} else if ((i == rows-1) && (j == 0)) {
		return {{-1,0}, {-1,1}, {0,1}};
	} else if (j == 0) {
		return {{-1,0}, {-1,1}, {0,1}, {1,1}, {1,0}};
	} else if (j == cols-1) {
		return {{-1,0}, {-1,-1}, {0,-1}, {1,-1}, {1,0}};
	} else if (i == 0) {
		return {{0,-1}, {1,-1}, {1,0}, {1,1}, {0,1}};
	} else if (i == rows-1) {
		return {{0,-1}, {-1,-1}, {-1,0}, {-1,1}, {0,1}};
	} else {
		return {{-1,-1}, {-1,0}, {-1,1}, {0,1}, {1,1}, {1,0}, {1,-1}, {0,-1}};
	}
}

int FadbadDepthMap::num_neighbors_empty(int i, int j) {
	int num_empty = 0;
	std::vector<std::vector<int> > o = offsets(i, j, H_SUB, W_SUB);
	for(int k=0; k < o.size(); ++k) {
		num_empty += (pixel_buckets[i+o[k][0]][j+o[k][1]]->is_empty()) ? 1 : 0;
	}

	return num_empty;
}

Vector3b FadbadDepthMap::average_of_neighbors(int i, int j) {
	std::vector<std::vector<int> > o = offsets(i, j, H_SUB, W_SUB);
	Vector3b avg_pt = Vector3b::Zero();
	bdouble num_neighbors = o.size();
	for(int k=0; k < o.size(); ++k) {
		Vector3b neighbor_avg = pixel_buckets[i+o[k][0]][j+o[k][1]]->average_point();
		for(int l=0; l < 3; ++l) { avg_pt(l) += (1/num_neighbors)*neighbor_avg(l); }
	}
	return avg_pt;
}
