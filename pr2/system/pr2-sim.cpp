#include "pr2-sim.h"

/**
 * PR2 Constructors/Destructors
 */

PR2::PR2(bool view) {
	std::string working_dir = boost::filesystem::current_path().normalize().string();
	std::string bsp_dir = working_dir.substr(0,working_dir.find("bsp"));
//	std::string env_file = bsp_dir + "bsp/pr2/envs/pr2-test.env.xml";
	std::string env_file = bsp_dir + "bsp/pr2/envs/pr2-empty.env.xml";

	std::string robot_name = "Brett";

	this->init(env_file, robot_name, view);
}

PR2::PR2(std::string env_file, std::string robot_name, bool view) {
	this->init(env_file, robot_name, view);
}

PR2::~PR2() {
	free(larm);
	free(rarm);
	free(head);
	free(hcam);
	free(lcam);
	free(rcam);
	env->Destroy();
}

/**
 * PR2 Initializers
 */

void SetViewer(rave::EnvironmentBasePtr penv, const std::string& viewername)
{
    rave::ViewerBasePtr viewer = RaveCreateViewer(penv,viewername);
    BOOST_ASSERT(!!viewer);

    // attach it to the environment:
    penv->Add(viewer);

    // finally call the viewer's infinite loop (this is why a separate thread is needed)
    bool showgui = true;
    viewer->main(showgui);
}

void PR2::init(std::string env_file, std::string robot_name, bool view) {
	RAVELOG_INFO("Initializing OpenRAVE\n");
	rave::RaveInitialize(true, rave::Level_Info);
	env = rave::RaveCreateEnvironment();
	RAVELOG_INFO("Loading environment: " + env_file + "\n");
	env->Load(env_file);

	robot = env->GetRobot(robot_name);

	if (view) {
		viewer_thread = boost::shared_ptr<boost::thread>(new boost::thread(boost::bind(SetViewer, env, "qtcoin")));
	}

	rarm = new Arm(robot, Arm::ArmType::right);
	larm = new Arm(robot, Arm::ArmType::left);
//	head = new Head(robot);

//	hcam = new Camera(robot, "head_cam", 5);
	rcam = new Camera(robot, "r_gripper_cam", rarm);
	lcam = new Camera(robot, "l_gripper_cam", larm);

	rave::CollisionCheckerBasePtr collision_checker = rave::RaveCreateCollisionChecker(env, "ode");
	env->SetCollisionChecker(collision_checker);
}


Arm::Arm(rave::RobotBasePtr robot, ArmType arm_type) {
	this->robot = robot;
	this->arm_type = arm_type;

	if (arm_type == ArmType::left) {
		manip_name = "leftarm";
	} else {
		manip_name = "rightarm";
	}

	manip = robot->GetManipulator(manip_name);
	joint_indices = manip->GetArmIndices();
	num_joints = joint_indices.size();

	std::vector<double> lower_vec(num_joints), upper_vec(num_joints);
	robot->GetDOFLimits(lower_vec, upper_vec, joint_indices);

	lower = Matrix<double,ARM_DIM,1>(lower_vec.data());
	upper = Matrix<double,ARM_DIM,1>(upper_vec.data());

	origin = robot->GetLink("torso_lift_link")->GetTransform();

	arm_joint_axes = { rave::Vector(0,0,1),
			rave::Vector(0,1,0),
			rave::Vector(1,0,0),
			rave::Vector(0,1,0),
			rave::Vector(1,0,0),
			rave::Vector(0,1,0),
			rave::Vector(1,0,0) };

	arm_link_trans = { rave::Vector(0, (arm_type == ArmType::left) ? 0.188 : -0.188, 0),
			rave::Vector(0.1, 0, 0),
			rave::Vector(0, 0, 0),
			rave::Vector(0.4, 0, 0),
			rave::Vector(0, 0, 0),
			rave::Vector(.321, 0, 0),
			rave::Vector(.18, 0, 0) };
}

/**
 * Arm public methods
 */

Matrix<double,ARM_DIM,1> Arm::get_joint_values() {
	std::vector<double> joint_values;
	robot->GetDOFValues(joint_values, joint_indices);
	return Matrix<double,ARM_DIM,1>(joint_values.data());
}

void Arm::get_limits(Matrix<double,ARM_DIM,1>& l, Matrix<double,ARM_DIM,1>& u) {
	l = lower;
	u = upper;
}

Matrix4d Arm::get_pose(const Matrix<double,ARM_DIM,1>& j) {
	rave::Transform pose_mat = origin;

	rave::Transform R;
	for(int i=0; i < ARM_DIM; ++i) {
		R.rot = rave::geometry::quatFromAxisAngle(arm_joint_axes[i], j[i]);
		R.trans = arm_link_trans[i];
		pose_mat = pose_mat*R;
	}

	return rave_utils::rave_to_eigen(pose_mat);
}

void Arm::set_joint_values(const Matrix<double,ARM_DIM,1>& j) {
	std::vector<double> joint_values_vec(ARM_DIM);
	for(int i=0; i < num_joints; ++i) {
		joint_values_vec[i] = std::min(j(i), upper(i));
		joint_values_vec[i] = std::max(j(i), lower(i));
	}

	robot->SetDOFValues(joint_values_vec, rave::KinBody::CheckLimitsAction::CLA_Nothing, joint_indices);
}

void Arm::set_pose(const rave::Transform& pose, std::string ref_frame) {
	std::vector<double> joint_values_vec(num_joints);
	rave_utils::cart_to_joint(manip, pose, ref_frame, "end_effector", joint_values_vec);

	Matrix<double,ARM_DIM,1> joint_values(joint_values_vec.data());
	set_joint_values(joint_values);
}

void Arm::set_posture(Posture posture) {
	std::vector<double> j;

	// joint values defined for left arm
	switch(posture) {
	case Posture::untucked:
		j = {0.4,  1.0,   0.0,  -2.05,  0.0,  -0.1,  0.0};
		break;
	case Posture::tucked:
		j = {0.06, 1.25, 1.79, -1.68, -1.73, -0.10, -0.09};
		break;
	case Posture::up:
		j = {0.33, -0.35,  2.59, -0.15,  0.59, -1.41, -0.27};
		break;
	case Posture::side:
		j = {1.832,  -0.332,   1.011,  -1.437,   1.1  ,  -2.106,  3.074};
		break;
	case Posture::mantis:
		j = {2.03018192, -0.05474993, 1.011, -1.47618716,  0.55995636, -1.42855926,  3.96467305};
		break;
	default:
		j = {0, 0, 0, 0, 0, 0, 0};
		break;
	}

	if (arm_type == ArmType::right) {
		j = {-j[0], j[1], -j[2], j[3], -j[4], j[5], -j[6]};
	}

	Matrix<double,ARM_DIM,1> j_vec(j.data());
	set_joint_values(j_vec);
}

void Arm::teleop() {
	double pos_step = .01;
	std::map<int,rave::Vector> delta_position =
	{
			{'a' , rave::Vector(0, pos_step, 0)},
			{'d' , rave::Vector(0, -pos_step, 0)},
			{'w' , rave::Vector(pos_step, 0, 0)},
			{'x' , rave::Vector(-pos_step, 0, 0)},
			{'+' , rave::Vector(0, 0, pos_step)},
			{'-' , rave::Vector(0, 0, -pos_step)},
	};

	double angle_step = 2.0*(M_PI/180);
	std::map<int,rave::Vector> delta_angle =
	{
			{'p', rave::Vector(angle_step, 0, 0)},
			{'o', rave::Vector(-angle_step, 0, 0)},
			{'k', rave::Vector(0, angle_step, 0)},
			{'l', rave::Vector(0, -angle_step, 0)},
			{'n', rave::Vector(0, 0, angle_step)},
			{'m', rave::Vector(0, 0, -angle_step)},
	};

	std::cout << manip_name << " teleop\n";

	char c;
	while ((c = utils::getch()) != 'q') {

		rave::Transform pose = manip->GetEndEffectorTransform();
		if (delta_position.count(c) > 0) {
			pose.trans += delta_position[c];
		} else if (delta_angle.count(c) > 0) {
			pose.rot = rave::geometry::quatFromAxisAngle(rave::geometry::axisAngleFromQuat(pose.rot) + delta_angle[c]);
		}

		set_pose(pose);
	}

	std::cout << manip_name << " end teleop\n";
}



///**
// * Head constructors
// */
//
//Head::Head(rave::RobotBasePtr robot) {
//	this->robot = robot;
//
//	std::vector<std::string> joint_names = {"head_pan_joint", "head_tilt_joint"};
//	for(int i=0; i < joint_names.size(); ++i) {
//		joint_indices.push_back(robot->GetJointIndex(joint_names[i]));
//	}
//	num_joints = joint_indices.size();
//
//	pose_link = robot->GetLink("wide_stereo_link");
//
//	std::vector<double> lower_vec(num_joints), upper_vec(num_joints);
//	robot->GetDOFLimits(lower_vec, upper_vec, joint_indices);
//
//	lower = Matrix<double,HEAD_DIM,1>(lower_vec.data());
//	upper = Matrix<double,HEAD_DIM,1>(upper_vec.data());
//}
//
///**
// * Head public methods
// */
//
//Matrix<double,HEAD_DIM,1> Head::get_joint_values() {
//	std::vector<double> joint_values;
//	robot->GetDOFValues(joint_values, joint_indices);
//	return Matrix<double,HEAD_DIM,1>(joint_values.data());
//}
//
//void Head::get_limits(Matrix<double,HEAD_DIM,1>& l, Matrix<double,HEAD_DIM,1>& u) {
//	l = lower;
//	u = upper;
//}
//
//rave::Transform Head::get_pose() {
//	return pose_link->GetTransform();
//}
//
//void Head::set_joint_values(Matrix<double,HEAD_DIM,1>& j) {
//	std::vector<double> joint_values_vec(HEAD_DIM);
//	for(int i=0; i < num_joints; ++i) {
//		joint_values_vec[i] = std::min(j(i), upper(i));
//		joint_values_vec[i] = std::max(j(i), lower(i));
//	}
//
//	robot->SetDOFValues(joint_values_vec, rave::KinBody::CheckLimitsAction::CLA_Nothing, joint_indices);
//}
//
//void Head::look_at(const rave::Transform &pose, const std::string ref_frame) {
//	rave::Transform world_from_ref, world_from_cam, ref_from_cam;
//
//	if (ref_frame == "world") {
//		world_from_ref.identity();
//	} else {
//		world_from_ref = robot->GetLink(ref_frame)->GetTransform();
//	}
//	world_from_cam = get_pose();
//	ref_from_cam = world_from_ref.inverse()*world_from_cam;
//
//	rave::Vector ax = pose.trans - ref_from_cam.trans;
//	double pan = atan(ax.y/ax.x);
//	double tilt = asin(-ax.z/sqrt(ax.lengthsqr3()));
//
//	Matrix<double,HEAD_DIM,1> joint_values(pan, tilt);
//	set_joint_values(joint_values);
//}
//
//void Head::teleop() {
//	double pos_step = .01;
//	std::map<int,std::vector<double> > delta_joints =
//	{
//			{'a' , {pos_step, 0}},
//			{'d' , {-pos_step, 0}},
//			{'w' , {0, -pos_step}},
//			{'x' , {0, pos_step}},
//	};
//
//	std::cout << "Head teleop\n";
//
//	char c;
//	while ((c = utils::getch()) != 'q') {
//
//		Matrix<double,HEAD_DIM,1> j = get_joint_values();
//		if (delta_joints.count(c) > 0) {
//			j += Matrix<double,HEAD_DIM,1>(delta_joints[c].data());
//		}
//
//		set_joint_values(j);
//
//	}
//
//	std::cout << "Head end teleop\n";
//}

/**
 * Camera constructors
 */

Camera::Camera(rave::RobotBasePtr r, std::string camera_name, Arm* a) : robot(r), arm(a), fov(nullptr) {
	bool found_sensor = false;
	std::vector<rave::RobotBase::AttachedSensorPtr> sensors = robot->GetAttachedSensors();
	for(int i=0; i < sensors.size(); ++i) {
		if (sensors[i]->GetName() == camera_name) {
			sensor = sensors[i]->GetSensor();
			found_sensor = true;
			break;
		}
	}

	assert(found_sensor);

	rave::SensorBase::SensorGeometryPtr geom = sensor->GetSensorGeometry(rave::SensorBase::ST_Camera);

	boost::shared_ptr<rave::SensorBase::CameraGeomData> cam_geom =
			boost::static_pointer_cast<rave::SensorBase::CameraGeomData>(geom);

	Matrix3d P = Matrix3d::Zero();
	P(0,0) = fx_sub;
	P(1,1) = fy_sub;
	P(2,2) = 1;
	P(0,2) = cx_sub;
	P(1,2) = cy_sub;

	depth_map = new DepthMap(sensor, P);

	Matrix4d sensor_pose = rave_utils::rave_to_eigen(sensor->GetTransform());
	Matrix4d ee_pose = a->get_pose(a->get_joint_values());

	gripper_tool_to_sensor = ee_pose.inverse()*sensor_pose;
}

/**
 * Camera public methods
 */

StdVector3d Camera::get_pcl(const Matrix<double,ARM_DIM,1>& j) {
	StdVector3d pcl;
	arm->set_joint_values(j);
	RowVector3d origin_pos = rave_utils::rave_to_eigen(sensor->GetTransform().trans);

	MatrixXd dirs = get_directions(j, HEIGHT, WIDTH, HEIGHT_M, WIDTH_M);

	if (fov != nullptr) {
		free(fov);
	}
	fov = new Beam3d(origin_pos, origin_pos + dirs.row(HEIGHT*(WIDTH-1)),
			origin_pos + dirs.row(0),
			origin_pos + dirs.row(HEIGHT-1),
			origin_pos + dirs.row(HEIGHT*WIDTH-1));

	rave::EnvironmentBasePtr env = robot->GetEnv();
	rave::RAY ray;
	ray.pos = sensor->GetTransform().trans;
	rave::CollisionReportPtr report(new rave::CollisionReport());
	for(int i=0; i < dirs.rows(); ++i) {
		ray.dir.x = dirs(i,0);
		ray.dir.y = dirs(i,1);
		ray.dir.z = dirs(i,2);
		if (env->CheckCollision(ray, report)) {
			pcl.push_back(rave_utils::rave_to_eigen(report->contacts[0].pos));

		}
	}

	return pcl;
}

std::vector<std::vector<Beam3d> > Camera::get_beams(const Matrix<double,ARM_DIM,1>& j, const StdVector3d& pcl) {
//	TimerCollection tc;
//	tc.start("get beams time");

	std::vector<std::vector<Beam3d> > beams(H_SUB-1, std::vector<Beam3d>(W_SUB-1));
	Matrix4d cam_pose = get_pose(j);

	depth_map->clear();
	for(int i=0; i < pcl.size(); ++i) {
		depth_map->add_point(pcl[i], cam_pose);
	}
	Matrix<double,H_SUB,W_SUB> z_buffer = depth_map->get_z_buffer(get_position(j));

	RowVector3d origin_pos = get_pose(j).block<3,1>(0,3);

	Matrix<double,N_SUB,3> dirs = get_directions(j, H_SUB, W_SUB, H_SUB_M, W_SUB_M);

	std::vector<std::vector<Vector3d> > hits(H_SUB, std::vector<Vector3d>(W_SUB));

	for(int i=0; i < H_SUB; ++i) {
		for(int j=0; j < W_SUB; ++j) {

			hits[i][j] = origin_pos + z_buffer(i,j)*(dirs.row(j*H_SUB+i) / dirs.row(j*H_SUB+i).norm());
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

//	tc.stop("get beams time");
//	tc.print_all_elapsed();

	return beams;
}


std::vector<Triangle3d> Camera::get_border(const std::vector<std::vector<Beam3d> >& beams, bool with_side_border) {
	std::vector<Triangle3d> border;

	int rows = beams.size(), cols = beams[0].size();

	if (with_side_border) {
		// deal with left and right border columns
		for(int i=0; i < rows; ++i) {
			border.push_back(Triangle3d(beams[i][0].base, beams[i][0].b, beams[i][0].c));
			border.push_back(Triangle3d(beams[i][cols-1].base, beams[i][cols-1].a, beams[i][cols-1].d));
		}

		// deal with top and bottom border rows
		for(int j=0; j < cols; ++j) {
			border.push_back(Triangle3d(beams[0][j].base, beams[0][j].a, beams[0][j].b));
			border.push_back(Triangle3d(beams[rows-1][j].base, beams[rows-1][j].c, beams[rows-1][j].d));
		}
	}

	const double min_sep = epsilon;
	// connect with left and top and add center
	for(int i=0; i < rows; ++i) {
		for(int j=0; j < cols; ++j) {
			const Beam3d& curr = beams[i][j];
			// left
			if (i > 0) {
				const Beam3d& left = beams[i-1][j];
				if (((left.a - curr.b).norm() > min_sep) || ((left.d - curr.c).norm() > min_sep)) {
					// not touching, add linkage
					border.push_back(Triangle3d(left.a, left.d, curr.b));
					border.push_back(Triangle3d(left.b, curr.b, curr.c));
				}
			}

			// top
			if (j > 0) {
				const Beam3d& top = beams[i][j-1];
				if (((top.c - curr.b).norm() > min_sep) || ((top.d - curr.a).norm() > min_sep)) {
					// not touching, add linkage
					border.push_back(Triangle3d(top.c, top.d, curr.b));
					border.push_back(Triangle3d(top.b, curr.b, curr.a));
				}
			}

			border.push_back(Triangle3d(curr.a, curr.b, curr.c));
			border.push_back(Triangle3d(curr.a, curr.c, curr.d));
		}
	}

//	std::cout << "before pruning size: " << border.size() << "\n";

	std::vector<Triangle3d> pruned_border;
	for(int i=0; i < border.size(); ++i) {
		if (border[i].area() > epsilon) {
			pruned_border.push_back(border[i]);
		}
	}

//	std::cout << "after pruning size: " << pruned_border.size() << "\n";

	return pruned_border;
}

bool Camera::is_inside(const Vector3d& p, std::vector<std::vector<Beam3d> >& beams) {
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

double Camera::signed_distance(const Vector3d& p, std::vector<std::vector<Beam3d> >& beams, std::vector<Triangle3d>& border) {
	double sd_sign = (is_inside(p, beams)) ? -1 : 1;

	double sd = INFINITY;
	for(int i=0; i < border.size(); ++i) {
		sd = std::min(sd, border[i].distance_to(p));
	}

	return (sd_sign*sd);
}

void Camera::plot_fov(std::vector<std::vector<Beam3d> >& beams) {
//	rave::EnvironmentBasePtr env = sensor->GetEnv();
//	rave::KinBodyPtr table = env->GetKinBody("table");
//	env->Remove(table);

	Vector3d color(0,1,0);
	// plot the ends
	for(int i=0; i < beams.size(); ++i) {
		for(int j=0; j < beams[i].size(); ++j) {
			beams[i][j].plot(sensor->GetEnv(), color);
		}
	}

//	int rows = beams.size(), cols = beams[0].size();
//	// plot the left and right
//	for(int i=0; i < rows; ++i) {
//		rave_utils::plot_segment(sensor->GetEnv(), beams[i][0].base, beams[i][0].b, color);
//		rave_utils::plot_segment(sensor->GetEnv(), beams[i][0].base, beams[i][0].c, color);
//
//		rave_utils::plot_segment(sensor->GetEnv(), beams[i][cols-1].base, beams[i][cols-1].a, color);
//		rave_utils::plot_segment(sensor->GetEnv(), beams[i][cols-1].base, beams[i][cols-1].d, color);
//	}
//	// plot the top and bottom
//	for(int j=0; j < cols; ++j) {
//		rave_utils::plot_segment(sensor->GetEnv(), beams[0][j].base, beams[0][j].a, color);
//		rave_utils::plot_segment(sensor->GetEnv(), beams[0][j].base, beams[0][j].b, color);
//
//		rave_utils::plot_segment(sensor->GetEnv(), beams[rows-1][j].base, beams[rows-1][j].c, color);
//		rave_utils::plot_segment(sensor->GetEnv(), beams[rows-1][j].base, beams[rows-1][j].d, color);
//	}
}

void Camera::plot_pcl(const StdVector3d& pcl) {
	rave::EnvironmentBasePtr env = sensor->GetEnv();

	Vector3d color(1,0,0);
	for(int i=0; i < pcl.size(); ++i) {
		rave_utils::plot_point(env, pcl[i], color, .005);
	}

	fov->plot(env, color);
	rave_utils::plot_segment(env, fov->base, fov->a, color);
	rave_utils::plot_segment(env, fov->base, fov->b, color);
	rave_utils::plot_segment(env, fov->base, fov->c, color);
	rave_utils::plot_segment(env, fov->base, fov->d, color);
}

/**
 * Camera Private methods
 */


MatrixXd Camera::get_directions(const Matrix<double,ARM_DIM,1>& j, const int h, const int w, const double h_meters, const double w_meters) {
	const int n = h*w;
	MatrixXd height_grid = VectorXd::LinSpaced(h, -h_meters/2.0, h_meters/2.0).replicate(1,w);
	MatrixXd width_grid = RowVectorXd::LinSpaced(w, -w_meters/2.0, w_meters/2.0).replicate(h,1);

	MatrixXd height_grid_vec(Map<VectorXd>(height_grid.data(), n));
	MatrixXd width_grid_vec(Map<VectorXd>(width_grid.data(), n));
	VectorXd z_grid = VectorXd::Zero(n,1);

	MatrixXd offsets(n,3);
	offsets << width_grid_vec, height_grid_vec, z_grid;

	MatrixXd points_cam = RowVector3d(0,0,MAX_RANGE).replicate(n,1) + (MAX_RANGE/FOCAL_LENGTH)*offsets;

	Matrix4d ref_from_world = get_pose(j);
	Vector3d origin_world_pos = ref_from_world.block<3,1>(0,3);

	MatrixXd directions(n,3);

	Matrix4d point_cam = Matrix4d::Identity();
	Vector3d point_world;
	for(int i=0; i < n; ++i) {
		point_cam.block<3,1>(0,3) = points_cam.row(i);
		point_world = (ref_from_world*point_cam).block<3,1>(0,3);

		directions.row(i) = point_world - origin_world_pos;
	}

	return directions;
}


/**
 * DepthMap Constructors
 */

DepthMap::DepthMap(rave::SensorBasePtr s, const Matrix3d& P_mat) : sensor(s), P(P_mat) {
	pixel_buckets = std::vector<std::vector<PixelBucket*> >(H_SUB, std::vector<PixelBucket*>(W_SUB));
	for(int i=0; i < H_SUB; ++i) {
		for(int j=0; j < W_SUB; ++j) {
			pixel_buckets[i][j] = new PixelBucket({i+0.5,j+0.5});
		}
	}
};

/**
 * DepthMap Public methods
 */

// TODO: only add point to bucket if point is either (1) closer than points in bucket or (2) closet to points in bucket
void DepthMap::add_point(const Vector3d& point, const Matrix4d& cam_pose) {
	Matrix4d point_mat = Matrix4d::Identity();
	point_mat.block<3,1>(0,3) = point;

	Matrix4d point_mat_tilde = cam_pose.inverse()*point_mat;

	Vector3d y = P*point_mat_tilde.block<3,1>(0,3);
	Vector2d pixel = {y(1)/y(2), y(0)/y(2)};
	int h_round = int(pixel(0));
	int w_round = int(pixel(1));

	if ((0 <= h_round) && (h_round < H_SUB) && (0 <= w_round) && (w_round < W_SUB) &&
			((cam_pose.block<3,1>(0,3) - point).norm() < MAX_RANGE)) { // TODO: should filter out points behind camera!
		pixel_buckets[h_round][w_round]->add_point(pixel, point);
	}
}

Matrix<double,H_SUB,W_SUB> DepthMap::get_z_buffer(const Vector3d& cam_pos) {
	Matrix<double,H_SUB,W_SUB> z_buffer = Matrix<double,H_SUB,W_SUB>::Ones();

	for(int i=0; i < H_SUB; ++i) {
		for(int j=0; j < W_SUB; ++j) {
			if (!(pixel_buckets[i][j]->is_empty())) {
				z_buffer(i,j) = (cam_pos - pixel_buckets[i][j]->average_point()).norm();
			} else if (num_neighbors_empty(i,j) >= 5 ) {
				z_buffer(i,j) = (cam_pos - average_of_neighbors(i, j)).norm();
			} else {
				z_buffer(i,j) = MAX_RANGE;
			}
//
//			if (pixel_buckets[i][j]->is_empty()) {
//				z_buffer(i,j) = MAX_RANGE;
//			} else {
//				z_buffer(i,j) = (cam_pos - pixel_buckets[i][j]->average_point()).norm();
//			}
		}
	}

	return z_buffer;
}

void DepthMap::clear() {
	for(int i=0; i < H_SUB; ++i) {
		for(int j=0; j < W_SUB; ++j) {
			pixel_buckets[i][j]->clear();
		}
	}
}

/**
 * DepthMap Private methods
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

int DepthMap::num_neighbors_empty(int i, int j) {
	int num_empty = 0;
	std::vector<std::vector<int> > o = offsets(i, j, H_SUB, W_SUB);
	for(int k=0; k < o.size(); ++k) {
		num_empty += (pixel_buckets[i+o[k][0]][j+o[k][1]]->is_empty()) ? 1 : 0;
	}

	return num_empty;
}

Vector3d DepthMap::average_of_neighbors(int i, int j) {
	std::vector<std::vector<int> > o = offsets(i, j, H_SUB, W_SUB);
	Vector3d avg_pt = Vector3d::Zero();
	double num_neighbors = o.size();
	for(int k=0; k < o.size(); ++k) {
		avg_pt += (1/num_neighbors)*(pixel_buckets[i+o[k][0]][j+o[k][1]]->average_point());
	}
	return avg_pt;
}

