#include "../include/pr2_sim.h"

/**
 * PR2 Constructors/Destructors
 */

PR2::PR2(bool view) {
	std::string working_dir = boost::filesystem::current_path().normalize().string();
	std::string bsp_dir = working_dir.substr(0,working_dir.find("bsp"));
	std::string env_file = bsp_dir + "bsp/pf/eih/envs/pr2-test.env.xml";

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
	free(h_kinect);
	free(l_kinect);
	free(r_kinect);
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

	larm = new Arm(robot, Arm::ArmType::left);
	rarm = new Arm(robot, Arm::ArmType::right);
	head = new Head(robot);

	h_kinect = new KinectSensor(robot, "head_depth", "head_cam");
	l_kinect = new KinectSensor(robot, "l_gripper_depth", "l_gripper_cam");
	r_kinect = new KinectSensor(robot, "r_gripper_depth", "r_gripper_cam");
}

/**
 * Manipulator public methods
 */

mat Manipulator::get_joint_values() {
	std::vector<double> joint_values;
	robot->GetDOFValues(joint_values, joint_indices);
	return conv_to<mat>::from(joint_values);
}

void Manipulator::get_limits(mat &l, mat &u) {
	l = lower;
	u = upper;
}

void Manipulator::set_joint_values(mat j) {
	for(int i=0; i < num_joints; ++i) {
		j(i) = std::min(j(i), upper(i));
		j(i) = std::max(j(i), lower(i));
	}

	std::vector<double> joint_values_vec = conv_to<std::vector<double>>::from(j);

	robot->SetDOFValues(joint_values_vec, rave::KinBody::CheckLimitsAction::CLA_Nothing, joint_indices);
}

/**
 * Manipulator protected methods
 */

void Manipulator::init() {
	std::vector<double> lower_vec(num_joints), upper_vec(num_joints);
	robot->GetDOFLimits(lower_vec, upper_vec, joint_indices);

	lower = conv_to<mat>::from(lower_vec);
	upper = conv_to<mat>::from(upper_vec);
}

/**
 * Arm Constructors
 */

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

	init();
}

/**
 * Arm public methods
 */

rave::Transform Arm::get_pose() {
	return manip->GetEndEffectorTransform();
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

	mat j_mat = conv_to<mat>::from(j);
	set_joint_values(j_mat);
}

void Arm::set_pose(const rave::Transform &pose, std::string ref_frame) {
	std::vector<double> joint_values_vec(num_joints);
	rave_utils::cart_to_joint(manip, pose, ref_frame, "end_effector", joint_values_vec);

	mat joint_values = conv_to<mat>::from(joint_values_vec);
	set_joint_values(joint_values);
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

		rave::Transform pose = get_pose();
		if (delta_position.count(c) > 0) {
			pose.trans += delta_position[c];
//			std::cout << "\r " << c;
		} else if (delta_angle.count(c) > 0) {
			pose.rot = rave::geometry::quatFromAxisAngle(rave::geometry::axisAngleFromQuat(pose.rot) + delta_angle[c]);
//			std::cout << "\r " << c;
		}

		set_pose(pose);
	}

	std::cout << manip_name << " end teleop\n";
}



/**
 * Head constructors
 */

Head::Head(rave::RobotBasePtr robot) {
	this->robot = robot;

	std::vector<std::string> joint_names = {"head_pan_joint", "head_tilt_joint"};
	for(int i=0; i < joint_names.size(); ++i) {
		joint_indices.push_back(robot->GetJointIndex(joint_names[i]));
	}
	num_joints = joint_indices.size();

	pose_link = robot->GetLink("wide_stereo_link");

	init();
}

/**
 * Head public methods
 */

rave::Transform Head::get_pose() {
	return pose_link->GetTransform();
}

// should be called look_at (named set_pose b/c of Manipulator)
// points head to look at pose
void Head::set_pose(const rave::Transform &pose, const std::string ref_frame) {
	rave::Transform world_from_ref, world_from_cam, ref_from_cam;

	if (ref_frame == "world") {
		world_from_ref.identity();
	} else {
		world_from_ref = robot->GetLink(ref_frame)->GetTransform();
	}
	world_from_cam = get_pose();
	ref_from_cam = world_from_ref.inverse()*world_from_cam;

	rave::Vector ax = pose.trans - ref_from_cam.trans;
	double pan = atan(ax.y/ax.x);
	double tilt = asin(-ax.z/sqrt(ax.lengthsqr3()));

	mat joint_values;
	joint_values << pan << tilt;
	set_joint_values(joint_values);
}

void Head::teleop() {
	double pos_step = .01;
	std::map<int,std::vector<double> > delta_joints =
	{
			{'a' , {pos_step, 0}},
			{'d' , {-pos_step, 0}},
			{'w' , {0, -pos_step}},
			{'x' , {0, pos_step}},
	};

	std::cout << "Head teleop\n";

	char c;
	while ((c = utils::getch()) != 'q') {

		mat j = get_joint_values();
		if (delta_joints.count(c) > 0) {
			j += conv_to<mat>::from(delta_joints[c]);
		}

		set_joint_values(j);

	}

	std::cout << "Head end teleop\n";
}



/**
 * Sensor constructors
 */

Sensor::Sensor(rave::SensorBasePtr sensor) {
	this->sensor = sensor;
	is_powered = false;
	is_rendering = false;

	std::vector<rave::SensorBase::SensorType> types =
	{
			rave::SensorBase::SensorType::ST_Camera,
			rave::SensorBase::SensorType::ST_Laser,
	};

	int i = 0;
	for(; i < types.size(); ++i) {
		if (sensor->Supports(types[i])) {
			type = types[i];
			break;
		}
	}

	if (i == types.size()) {
		RAVELOG_ERROR("Invalid sensor type. Exiting.\n");
		exit(0);
	}


}

/**
 * Sensor public methods
 */

void Sensor::power_on() {
	if (!is_powered) {
		sensor->Configure(rave::SensorBase::ConfigureCommand::CC_PowerOn);
		is_powered = true;
	}
}

void Sensor::power_off() {
	if (is_powered) {
		sensor->Configure(rave::SensorBase::ConfigureCommand::CC_PowerOff);
		is_powered = false;
	}
}

void Sensor::render_on() {
	if (!is_rendering) {
		sensor->Configure(rave::SensorBase::ConfigureCommand::CC_RenderDataOn);
		is_rendering = true;
	}
}

void Sensor::render_off() {
	if (is_rendering) {
		sensor->Configure(rave::SensorBase::ConfigureCommand::CC_RenderDataOff);
		is_rendering = false;
	}
}

rave::SensorBase::SensorDataPtr Sensor::get_data() {
	rave::SensorBase::SensorDataPtr data = sensor->CreateSensorData(type);
	bool success = sensor->GetSensorData(data);

	if (!success) {
		RAVELOG_WARN("Couldn't retrieve data\n");
	}

	return data;
}

rave::Transform Sensor::get_pose() {
	return sensor->GetTransform();
}

/**
 * DepthSensor constructor
 */

DepthSensor::DepthSensor(rave::SensorBasePtr sensor) : Sensor(sensor) {

}

/**
 * DepthSensor public methods
 */

std::vector<mat> DepthSensor::get_points() {
	boost::shared_ptr<rave::SensorBase::LaserSensorData> l_data =
			boost::static_pointer_cast<rave::SensorBase::LaserSensorData>(get_data());
	std::vector<double> in_range = l_data->intensity;

	rave::Transform sensor_pose = sensor->GetTransform();
	std::vector<mat> points;
	for(int i=0; i < l_data->ranges.size(); ++i) {
		if (in_range[i] > .99) {
			rave::Vector p = sensor_pose.trans + l_data->ranges[i];
			mat point;
			point << p.x << p.y << p.z;
			points.push_back(point);
		}
	}

	return points;
}


/**
 * CameraSensor constructor
 */

CameraSensor::CameraSensor(rave::SensorBasePtr sensor) : Sensor(sensor) {
	rave::SensorBase::SensorGeometryPtr geom = sensor->GetSensorGeometry(rave::SensorBase::ST_Camera);

	boost::shared_ptr<rave::SensorBase::CameraGeomData> cam_geom =
			boost::static_pointer_cast<rave::SensorBase::CameraGeomData>(geom);

	P = zeros(3,3);
	P(0,0) = cam_geom->KK.fx;
	P(1,1) = cam_geom->KK.fy;
	P(2,2) = 1;
	P(0,2) = cam_geom->KK.cx;
	P(1,2) = cam_geom->KK.cy;

	height = cam_geom->height;
	width = cam_geom->width;
}

/**
 * Public methods
 */

// all values [0, 1]
cube CameraSensor::get_image() {
	boost::shared_ptr<rave::SensorBase::CameraSensorData> c_data =
				boost::static_pointer_cast<rave::SensorBase::CameraSensorData>(get_data());
	std::vector<uint8_t> raw_image = c_data->vimagedata;

	cube image(height, width, 3);
	int index = 0;
	for(int i=0; i < height; ++i) {
		for(int j=0; j < width; ++j) {
			for(int k=0; k < 3; ++k) {
				image(i,j,k) = double(raw_image[index++])/255;
			}
		}
	}

	return image;
}

std::vector<std::vector<mat> > CameraSensor::get_pixels_and_colors(const std::vector<mat> &points) {
	cube image = get_image();

	std::vector<std::vector<mat> > points_pixels_colors;
	for(int i=0; i < points.size(); ++i) {
		mat pixel = get_pixel_from_point(points[i]);
		if ((0 <= pixel(0)) && (pixel(0) < height) &&
				(0 <= pixel(1)) && (pixel(1) < width)) {
			mat color = image.subcube(int(pixel(0)), int(pixel(1)), 0, int(pixel(0)), int(pixel(1)), 2);
			std::vector<mat> point_pixel_color = {points[i], pixel, color};
			points_pixels_colors.push_back(point_pixel_color);
		}
	}
	return points_pixels_colors;
}

mat CameraSensor::get_pixel_from_point(const mat &point) {
	mat x_mat = zeros(4,4);
	if (point.is_colvec()) {
		x_mat.submat(0,3,2,3) = point;
	} else {
		x_mat.submat(0,3,2,3) = point.t();
	}
	x_mat(3,3) = 1;
	mat pose = rave_utils::rave_transform_to_mat(get_pose());
	mat xtilde = inv(rave_utils::rave_transform_to_mat(get_pose()))*x_mat;

//	double x1 = xtilde(0,3), x2 = xtilde(1,3), x3 = xtilde(2,3);
//	double f = P(0,0);
//	double y1 = (f/x3)*x1, y2 = (f/x3)*x2;

	mat y = P*xtilde.submat(0,3,2,3);
	mat pixel(2,1);
	pixel(0) = floor(y(1)/y(2));
	pixel(1) =  floor(y(0)/y(2));

	return pixel;
}

bool CameraSensor::is_in_fov(const mat& point) {
	mat pixel = get_pixel_from_point(point);
	return ((0 <= pixel(0)) && (pixel(0) < height) &&
			(0 <= pixel(1)) && (pixel(1) < width));
}


/**
 * KinectSensor constructors
 */

KinectSensor::KinectSensor(rave::RobotBasePtr robot, std::string depth_sensor_name, std::string camera_sensor_name) {
	this->robot = robot;
	camera_sensor = NULL;
	depth_sensor = NULL;

	std::vector<rave::RobotBase::AttachedSensorPtr> sensors = robot->GetAttachedSensors();
	for(int i=0; i < sensors.size(); ++i) {
		if (sensors[i]->GetName() == depth_sensor_name) {
			depth_sensor = new DepthSensor(sensors[i]->GetSensor());
		} else if (sensors[i]->GetName() == camera_sensor_name) {
			camera_sensor = new CameraSensor(sensors[i]->GetSensor());
		}
	}

	if ((camera_sensor == NULL) || (depth_sensor == NULL)) {
		RAVELOG_ERROR("Invalid sensor names. Exiting.\n");
		exit(0);
	}
}

KinectSensor::~KinectSensor() {
	free(depth_sensor);
	free(camera_sensor);
}

/**
 * KinectSensor public methods
 */

void KinectSensor::power_on() {
	depth_sensor->power_on();
	camera_sensor->power_on();
}

void KinectSensor::power_off() {
	depth_sensor->power_off();
	camera_sensor->power_off();
}

void KinectSensor::render_on() {
	depth_sensor->render_on();
	camera_sensor->render_on();
}

void KinectSensor::render_off() {
	depth_sensor->render_off();
	camera_sensor->render_off();
}


std::vector<ColoredPoint*> KinectSensor::get_point_cloud() {
	std::vector<mat> points = depth_sensor->get_points();
	std::vector<std::vector<mat> > points_pixels_colors = camera_sensor->get_pixels_and_colors(points);

	std::vector<ColoredPoint*> colored_points(points_pixels_colors.size());
	for(int i=0; i < points_pixels_colors.size(); ++i) {
		colored_points[i] = new ColoredPoint(points_pixels_colors[i][0], points_pixels_colors[i][2]);
	}

	return colored_points;
}

mat KinectSensor::get_z_buffer() {
	std::vector<mat> points = depth_sensor->get_points();
	std::vector<std::vector<mat> > points_pixels_colors = camera_sensor->get_pixels_and_colors(points);

	rave::Transform pose = camera_sensor->get_pose();
	mat pose_pos;
	pose_pos << pose.trans.x << pose.trans.y << pose.trans.z;

	int height = camera_sensor->get_height(), width = camera_sensor->get_width();
	mat z_buffer = INFINITY*ones(height,width);
	for(int i=0; i < points_pixels_colors.size(); ++i) {
		mat point = points_pixels_colors[i][0];
		mat pixel = points_pixels_colors[i][1];
		z_buffer(int(pixel(0)), int(pixel(1))) = norm(point - pose_pos, 2);
	}

	return z_buffer;
}

void KinectSensor::display_point_cloud(const std::vector<ColoredPoint*> &colored_points, std::vector<rave::GraphHandlePtr> &handles) {
	for(int i=0; i < colored_points.size(); ++i) {
		handles.push_back(colored_points[i]->display(robot->GetEnv()));
	}
}


