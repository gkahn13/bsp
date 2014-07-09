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

//inline Matrix3d axis_angle_to_matrix(const double& angle, const Vector3d& axis) {
//	double sina = sin(angle);
//	double cosa = cos(angle);
//
//	Matrix3d R = Vector3d(cosa, cosa, cosa).asDiagonal();
//	Matrix3d axis_outer = axis*axis.transpose();
//	for(int i=0; i < 3; ++i) {
//		for(int j=0; j < 3; ++j) {
//			axis_outer(i,j) *= (1-cosa);
//		}
//	}
//	R += axis_outer;
//
//	Vector3d axis_sina = axis;
//	for(int i=0; i < 3; ++i) { axis_sina(i) *= sina; }
//	Matrix3d mat_sina;
//	mat_sina << 0, -axis_sina(2), axis_sina(1),
//			axis_sina(2), 0, -axis_sina(0),
//			-axis_sina(1), axis_sina(0), 0;
//	R += mat_sina;
//
//	return R;
//}

Matrix4d Arm::get_pose(const Matrix<double,ARM_DIM,1>& j) {
	rave::Transform pose_mat = origin;

	rave::Transform R;
	for(int i=0; i < ARM_DIM; ++i) {
//		Matrix3d rot = axis_angle_to_matrix(j(i), rave_utils::rave_to_eigen(arm_joint_axes[i]));
//		Vector3d trans = rave_utils::rave_to_eigen(arm_link_trans[i]);
//		Matrix4d R_eigen = Matrix4d::Identity();
//		R_eigen.block<3,3>(0,0) = rot;
//		R_eigen.block<3,1>(0,3) = trans;
//		std::cout << "R_eigen:\n" << R_eigen << "\n";
//		R = rave_utils::eigen_to_rave(R_eigen);

		R.rot = rave::geometry::quatFromAxisAngle(arm_joint_axes[i], j(i));
		R.trans = arm_link_trans[i];
		pose_mat = pose_mat*R;
	}

	return rave_utils::rave_to_eigen(pose_mat);
}

Matrix<double,3,ARM_DIM> Arm::get_position_jacobian(const Matrix<double,ARM_DIM,1>& j) {
	Matrix<double,3,ARM_DIM> jac;
	jac.setZero();

	Matrix<double,ARM_DIM,1> j_p = j, j_m = j;
	for(int i=0; i < ARM_DIM; ++i) {
		j_p(i) = j(i) + step;
		j_m(i) = j(i) - step;

		jac.block<3,1>(0,i) = (get_position(j_p) - get_position(j_m)) / (2*step);

		j_p(i) = j(i);
		j_m(i) = j(i);
	}

	return jac;
}

/**
 * \brief Find 3d translation IK
 * \param pos position to solve for
 * \param j   initial joints when searching for solution, and final solution
 * \return success/failure
 */
bool Arm::ik(const Vector3d& goal_pos, Matrix<double,ARM_DIM,1>& j) {
	Vector3d arm_pos, error, jac_jacT_error;
	Matrix<double,3,ARM_DIM> arm_pos_jac;

	int iter = 0, max_iter = 1000;
	while((goal_pos - get_position(j)).norm() > epsilon) {
		arm_pos = get_position(j);
		arm_pos_jac = get_position_jacobian(j);

		error = goal_pos - arm_pos;

		jac_jacT_error = arm_pos_jac*arm_pos_jac.transpose()*error;
		double alpha = error.dot(jac_jacT_error) / jac_jacT_error.dot(jac_jacT_error);
		j += alpha*arm_pos_jac.transpose()*error;

		if (iter++ > max_iter) {
			LOG_WARN("IK failed: reached max iterations");
			return false;
		}
	}

	if (((j - lower).minCoeff() < 0) || ((upper - j).minCoeff() < 0)) {
		return false;
	}

	return true;
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
	P(0,0) = intrinsics::fx_sub;
	P(1,1) = intrinsics::fy_sub;
	P(2,2) = 1;
	P(0,2) = intrinsics::cx_sub;
	P(1,2) = intrinsics::cy_sub;

	KK << intrinsics::fx, 0, intrinsics::cx,
			0, intrinsics::fy, intrinsics::cy,
			0, 0, 1;

	KK_SUB << intrinsics::fx_sub, 0, intrinsics::cx_sub,
				0, intrinsics::fy_sub, intrinsics::cy_sub,
				0, 0, 1;

//	depth_map = new DepthMap(sensor, P);

	Matrix4d sensor_pose = rave_utils::rave_to_eigen(sensor->GetTransform());
	Matrix4d ee_pose = a->get_pose(a->get_joint_values());

	gripper_tool_to_sensor = ee_pose.inverse()*sensor_pose;
}

/**
 * Camera public methods
 */

StdVector3d Camera::get_pc(const Matrix<double,ARM_DIM,1>& j) {
	StdVector3d pc;
	arm->set_joint_values(j);
	RowVector3d origin_pos = rave_utils::rave_to_eigen(sensor->GetTransform().trans);

	MatrixXd dirs = get_directions(j, HEIGHT_FULL, WIDTH_FULL, intrinsics::HEIGHT_M, intrinsics::WIDTH_M);

	if (fov != nullptr) {
		free(fov);
	}
	fov = new Beam3d(origin_pos, origin_pos + dirs.row(HEIGHT_FULL*(WIDTH_FULL-1)),
			origin_pos + dirs.row(0),
			origin_pos + dirs.row(HEIGHT_FULL-1),
			origin_pos + dirs.row(HEIGHT_FULL*WIDTH_FULL-1));

	rave::EnvironmentBasePtr env = robot->GetEnv();
	rave::RAY ray;
	ray.pos = sensor->GetTransform().trans;
	rave::CollisionReportPtr report(new rave::CollisionReport());
	for(int i=0; i < dirs.rows(); ++i) {
		ray.dir.x = dirs(i,0);
		ray.dir.y = dirs(i,1);
		ray.dir.z = dirs(i,2);
		if (env->CheckCollision(ray, report)) {
			pc.push_back(rave_utils::rave_to_eigen(report->contacts[0].pos));

		}
	}

	return pc;
}

Matrix<double,HEIGHT_FULL,WIDTH_FULL> Camera::get_full_zbuffer(const Matrix<double,ARM_DIM,1>& j, const StdVector3d& obstacles) {
	Matrix<double,HEIGHT_FULL,WIDTH_FULL> zbuffer = intrinsics::MAX_RANGE*Matrix<double,HEIGHT_FULL,WIDTH_FULL>::Ones();
	Matrix4d cam_pose = get_pose(j);

	for(int i=0; i < obstacles.size(); ++i) {
		Vector2i pixel = get_pixel_from_point(obstacles[i], cam_pose, KK);
		int h = pixel(0), w = pixel(1);
		if ((w >= 0) && (w < WIDTH_FULL) && (h >= 0) && (h < HEIGHT_FULL)) { // in frustrum
			Matrix4d obstacle_pose = Matrix4d::Identity();
			obstacle_pose.block<3,1>(0,3) = obstacles[i];
			if ((cam_pose.inverse()*obstacle_pose)(2,3) > 0) { // in front of camera pose
				double dist = (obstacles[i] - cam_pose.block<3,1>(0,3)).norm();
				zbuffer(h,w) = (zbuffer(h,w) < dist) ? zbuffer(h,w) : dist;
			}
		}
	}

	return zbuffer;
}

bool Camera::is_in_fov_full(const Vector3d& point, const Matrix<double,HEIGHT_FULL,WIDTH_FULL>& zbuffer, const Matrix4d& cam_pose) {
	Vector2i pixel = get_pixel_from_point(point, cam_pose, KK);
	int h = pixel(0), w = pixel(1);

	// out of frustrum
	if ((w < 0) || (w >= WIDTH_FULL) || (h < 0) || (h >= HEIGHT_FULL)) {
		return false;
	}

	// occluded
	if (zbuffer(h,w) < (cam_pose.block<3,1>(0,3) - point).norm()) {
		return false;
	}

	// behind camera pose
	Matrix4d point_pose = Matrix4d::Identity();
	point_pose.block<3,1>(0,3) = point;
	if ((cam_pose.inverse()*point_pose)(2,3) < 0) {
		return false;
	}

	return true;
}

inline std::vector<std::vector<int> > index_neighbors(int i, int j, int rows, int cols) {
	if ((i == 0) && (j == 0)) {
		return {{i+1,j+0}, {i+1,j+1}, {i+0,j+1}};
	} else if ((i == 0) && (j == cols-1)) {
		return {{i+0,j-1}, {i+1,j-1}, {i+1,j+0}};
	} else if ((i == rows-1) && (j == cols-1)) {
		return {{i-1,j+0}, {i-1,j-1}, {i+0,j-1}};
	} else if ((i == rows-1) && (j == 0)) {
		return {{i-1,j+0}, {i-1,j+1}, {i+0,j+1}};
	} else if (j == 0) {
		return {{i-1,j+0}, {i-1,j+1}, {i+0,j+1}, {i+1,j+1}, {i+1,j+0}};
	} else if (j == cols-1) {
		return {{i-1,j+0}, {i-1,j-1}, {i+0,j-1}, {i+1,j-1}, {i+1,j+0}};
	} else if (i == 0) {
		return {{i+0,j-1}, {i+1,j-1}, {i+1,j+0}, {i+1,j+1}, {i+0,j+1}};
	} else if (i == rows-1) {
		return {{i+0,j-1}, {i-1,j-1}, {i-1,j+0}, {i-1,j+1}, {i+0,j+1}};
	} else {
		return {{i-1,j-1}, {i-1,j+0}, {i-1,j+1}, {i+0,j+1}, {i+1,j+1}, {i+1,j+0}, {i+1,j-1}, {i+0,cols-1}};
	}
}

Matrix<double,H_SUB,W_SUB> Camera::get_zbuffer(const Matrix<double,ARM_DIM,1>& j, const StdVector3d& obstacles) {
	Matrix<double,H_SUB,W_SUB> zbuffer = intrinsics::MAX_RANGE*Matrix<double,H_SUB,W_SUB>::Ones();
	Matrix4d cam_pose = get_pose(j);

	for(int i=0; i < obstacles.size(); ++i) {
		Vector2i pixel = get_pixel_from_point(obstacles[i], cam_pose, KK_SUB);
		int h = pixel(0), w = pixel(1);
		if ((w >= 0) && (w < W_SUB) && (h >= 0) && (h < H_SUB)) { // in frustrum
			Matrix4d obstacle_pose = Matrix4d::Identity();
			obstacle_pose.block<3,1>(0,3) = obstacles[i];
			if ((cam_pose.inverse()*obstacle_pose)(2,3) > 0) { // in front of camera pose
				double dist = (obstacles[i] - cam_pose.block<3,1>(0,3)).norm();
				zbuffer(h,w) = (zbuffer(h,w) < dist) ? zbuffer(h,w) : dist;
			}
		}
	}

	Matrix<double,H_SUB,W_SUB> zbuffer_smoothed;
	zbuffer_smoothed.setZero();
	for(int i=0; i < H_SUB; ++i) {
		for(int j=0; j < W_SUB; ++j) {
			if (zbuffer(i,j) == intrinsics::MAX_RANGE) {
				// get rid of MAX_RANGEs surrounded by valids
				std::vector<std::vector<int> > n = index_neighbors(i, j, H_SUB, W_SUB);
				int num_not_max = 0;
				double depth_sum = 0;
				for(int k=0; k < n.size(); ++k) {
					if (zbuffer(n[k][0], n[k][1]) < intrinsics::MAX_RANGE) {
						num_not_max++;
						depth_sum += zbuffer(n[k][0], n[k][1]);
					}
				}

				if ((double(num_not_max)/double(n.size())) >= .5) {
					zbuffer_smoothed(i,j) = depth_sum / double(n.size());
				} else {
					zbuffer_smoothed(i,j) = zbuffer(i,j);
				}
			} else {
				// get rid of valids surrounded by MAX_RANGEs
				std::vector<std::vector<int> > n = index_neighbors(i, j, H_SUB, W_SUB);
				int num_max = 0;
				for(int k=0; k < n.size(); ++k) {
					if (zbuffer(n[k][0], n[k][1]) >= intrinsics::MAX_RANGE) {
						num_max++;
					}
				}

				if ((double(num_max)/double(n.size())) >=.75) {
					zbuffer_smoothed(i,j) = intrinsics::MAX_RANGE;
				} else {
					zbuffer_smoothed(i,j) = zbuffer(i,j);
				}
			}
		}
	}
	zbuffer = zbuffer_smoothed;

	return zbuffer;
}

Vector2i Camera::get_pixel_from_point(const Vector3d& point, const Matrix4d& cam_pose, const Matrix3d& P) {
	Matrix4d point_mat = Matrix4d::Identity();
	point_mat.block<3,1>(0,3) = point;

	Matrix4d point_mat_tilde = cam_pose.inverse()*point_mat;

	Vector3d y = P*point_mat_tilde.block<3,1>(0,3);

	Vector2i pixel = {int(y(1)/y(2)), int(y(0)/y(2))};
	return pixel;
}

Vector3d Camera::get_point_from_pixel_and_dist(const Vector2i& pixel, const double dist, const Matrix4d& cam_pose) {
	double x = intrinsics::FOCAL_LENGTH*(double(pixel(1)-W_SUB/2.0)/intrinsics::fx_sub);
	double y = intrinsics::FOCAL_LENGTH*(double(pixel(0)-H_SUB/2.0)/intrinsics::fy_sub);
	double z = dist;

	Matrix4d point_cam = Matrix4d::Identity();
	point_cam.block<3,1>(0,3) = Vector3d(x,y,z);
	Matrix4d point_world = cam_pose*point_cam;
	return point_world.block<3,1>(0,3);
}

bool Camera::is_in_fov(const Vector3d& point, const Matrix<double,H_SUB,W_SUB>& zbuffer, const Matrix4d& cam_pose) {
	Vector2i pixel = get_pixel_from_point(point, cam_pose, KK_SUB);
	int h = pixel(0), w = pixel(1);

	// out of frustrum
	if ((w < 0) || (w >= W_SUB) || (h < 0) || (h >= H_SUB)) {
		return false;
	}

	// occluded
	if (zbuffer(h,w) < (cam_pose.block<3,1>(0,3) - point).norm()) {
		return false;
	}

	// behind camera pose
	Matrix4d point_pose = Matrix4d::Identity();
	point_pose.block<3,1>(0,3) = point;
	if ((cam_pose.inverse()*point_pose)(2,3) < 0) {
		return false;
	}

	return true;
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

	MatrixXd points_cam = RowVector3d(0,0,intrinsics::MAX_RANGE).replicate(n,1) + (intrinsics::MAX_RANGE/intrinsics::FOCAL_LENGTH)*offsets;

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

