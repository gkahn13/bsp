#include "pr2_sim.h"

/**
 * PR2 Constructors/Destructors
 */

PR2::PR2() {
	std::string working_dir = boost::filesystem::current_path().normalize().string();
	std::string bsp_dir = working_dir.substr(0,working_dir.find("bsp"));
	std::string env_file = bsp_dir + "bsp/pf/eih/envs/pr2-test.env.xml";

	std::string robot_name = "Brett";

	this->init(env_file, robot_name, true);
}

PR2::PR2(std::string env_file, std::string robot_name, bool view) {
	this->init(env_file, robot_name, view);
}

PR2::~PR2() {
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
//		viewer_thread->join(); // TODO: temp
	}
//	env->Destroy(); // destroy

	larm = new Arm(robot, Arm::ArmType::left);
	rarm = new Arm(robot, Arm::ArmType::right);
	head = new Head(robot);

	std::vector<rave::RobotBase::AttachedSensorPtr> sensors = robot->GetAttachedSensors();
	for(int i=0; i < sensors.size(); ++i) {
		if (sensors[i]->GetName() == "r_gripper_cam") {
			std::cout << "Creating camera...\n";
			camera = new CameraSensor(sensors[i]->GetSensor());
		}
	}
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
	arm_indices = manip->GetArmIndices();
	num_joints = arm_indices.size();
}

/**
 * Arm public methods
 */

mat Arm::get_joint_values() {
	std::vector<double> joint_values;
	manip->GetArmDOFValues(joint_values);
	return conv_to<mat>::from(joint_values);
}

rave::Transform Arm::get_pose() {
	return manip->GetEndEffectorTransform();
}

void Arm::get_limits(mat& lower, mat& upper) {
	std::vector<double> lower_vec(num_joints), upper_vec(num_joints);
	robot->GetDOFLimits(lower_vec, upper_vec, arm_indices);

	lower = conv_to<mat>::from(lower_vec);
	upper = conv_to<mat>::from(upper_vec);
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

void Arm::set_joint_values(const mat &joint_values) {
	std::vector<double> joint_values_vec = conv_to<std::vector<double>>::from(joint_values);
	robot->SetDOFValues(joint_values_vec, rave::KinBody::CheckLimitsAction::CLA_CheckLimits, arm_indices);
}

void Arm::set_pose(const rave::Transform &pose, std::string ref_frame) {
	std::vector<double> joint_values_vec(num_joints);
	rave_utils::cart_to_joint(manip, pose, ref_frame, "end_effector", joint_values_vec);

	mat joint_values = conv_to<mat>::from(joint_values_vec);
	set_joint_values(joint_values);
}

void Arm::teleop() {
	initscr();
	raw();
	noecho();

	printw("%s teleop\n", manip_name.c_str());

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

	int c;
	do {
		boost::this_thread::sleep(boost::posix_time::seconds(.01));
		c = getch();

		if (c == ERR) {
			continue;
		}

		rave::Transform pose = get_pose();
		if (delta_position.count(c) > 0) {
			pose.trans += delta_position[c];
		} else if (delta_angle.count(c) > 0) {
			pose.rot = rave::geometry::quatFromAxisAngle(rave::geometry::axisAngleFromQuat(pose.rot) + delta_angle[c]);
		}

		set_pose(pose);

	} while (c != 'q');

	printw("%s end teleop\n", manip_name.c_str());

	refresh();
	endwin();
}



/**
 * Head constructors
 */

Head::Head(rave::RobotBasePtr robot) {
	this->robot = robot;

	std::vector<std::string> joint_names = {"head_pan_joint", "head_tilt_joint"};
	for(int i=0; i < joint_names.size(); ++i) {
		head_indices.push_back(robot->GetJointIndex(joint_names[i]));
	}
	num_joints = head_indices.size();
}

/**
 * Head public methods
 */

mat Head::get_joint_values() {
	std::vector<double> joint_values;
	robot->GetDOFValues(joint_values, head_indices);
	return conv_to<mat>::from(joint_values);
}

void Head::get_limits(mat &lower, mat &upper) {
	std::vector<double> lower_vec(num_joints), upper_vec(num_joints);
	robot->GetDOFLimits(lower_vec, upper_vec, head_indices);

	lower = conv_to<mat>::from(lower_vec);
	upper = conv_to<mat>::from(upper_vec);
}

void Head::set_joint_values(const mat &joint_values) {
	std::vector<double> joint_values_vec = conv_to<std::vector<double>>::from(joint_values);
	robot->SetDOFValues(joint_values_vec, rave::KinBody::CheckLimitsAction::CLA_CheckLimits, head_indices);
}

void Head::look_at(const rave::Transform &pose, const std::string reference_frame, const std::string camera_frame) {
	rave::Transform world_from_ref, world_from_cam, ref_from_cam;

	if (reference_frame == "world") {
		world_from_ref.identity();
	} else {
		world_from_ref = robot->GetLink(reference_frame)->GetTransform();
	}
	world_from_cam = robot->GetLink(camera_frame)->GetTransform();
	ref_from_cam = world_from_ref.inverse()*world_from_cam;

	rave::Vector ax = pose.trans - ref_from_cam.trans;
	double pan = atan(ax.y/ax.x);
	double tilt = asin(-ax.z/sqrt(ax.lengthsqr3()));

	mat joint_values;
	joint_values << pan << tilt;
	set_joint_values(joint_values);
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



/**
 * CameraSensor constructor
 */

CameraSensor::CameraSensor(rave::SensorBasePtr sensor) : Sensor(sensor) {
//	power_on();
//	boost::this_thread::sleep(boost::posix_time::seconds(1));
	std::cout << "CameraSensor constructor\n";

	boost::shared_ptr<rave::SensorBase::CameraSensorData> c_data =
			boost::static_pointer_cast<rave::SensorBase::CameraSensorData>(get_data());

	rave::SensorBase::SensorGeometryPtr geom = sensor->GetSensorGeometry(rave::SensorBase::ST_Camera);

	boost::shared_ptr<rave::SensorBase::CameraGeomData> cam_geom =
			boost::static_pointer_cast<rave::SensorBase::CameraGeomData>(geom);
	std::cout << "cam_geom\n";
	std::cout << "cx: " << cam_geom->KK.cx << "\n";
	std::cout << "cy: " << cam_geom->KK.cy << "\n";
	std::cout << "fx: " << cam_geom->KK.fx << "\n";
	std::cout << "fy: " << cam_geom->KK.fy << "\n";

	P = zeros(3,3);
	P(0,0) = cam_geom->KK.fx;
	P(1,1) = cam_geom->KK.fy;
	P(2,2) = 1;
	P(0,2) = cam_geom->KK.cx;
	P(1,2) = cam_geom->KK.cy;
	std::cout << P;

	height = cam_geom->height;
	width = cam_geom->width;
}

/**
 * Public methods
 */

cube CameraSensor::get_image() {
	boost::shared_ptr<rave::SensorBase::CameraSensorData> c_data =
				boost::static_pointer_cast<rave::SensorBase::CameraSensorData>(get_data());
	c_data->vimagedata;

	cube image(height, width, 3);
	int index = 0;
	for(int i=0; i < height; ++i) {
		for(int j=0; j < width; ++j) {
			for(int k=0; k < 3; ++k) {
				image(i,j,k) = c_data->vimagedata[index++];
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
		if ((0 <= pixel[0]) && (pixel[0] < height) &&
				(0 <= pixel[1]) && (pixel[1] < width)) {
			mat color = image.subcube(pixel[0], pixel[1], 0, pixel[0], pixel[1], 2);
			std::vector<mat> point_pixel_color = {points[i], pixel, color};
			points_pixels_colors.push_back(point_pixel_color);
		}
	}
	return points_pixels_colors;
}

mat CameraSensor::get_pixel_from_point(const mat &point) {
	mat x(point);
	x.submat(0,0,2,2) = zeros(3,3);
	mat xtilde = inv(rave_utils::rave_transform_to_mat(sensor->GetTransform()))*x;

//	double x1 = xtilde(0,3), x2 = xtilde(1,3), x3 = xtilde(2,3);
//	double f = P(0,0);
//	double y1 = (f/x3)*x1, y2 = (f/x3)*x2;

	mat y = P*xtilde.submat(0,2,2,2);
	mat pixel(2,1);
	pixel[0] = floor(y(1)/y(2));
	pixel[1] =  floor(y(0)/y(2));

	return pixel;
}

bool CameraSensor::is_in_fov(const mat& point) {
	mat pixel = get_pixel_from_point(point);
	return ((0 <= pixel[0]) && (pixel[0] < height) &&
			(0 <= pixel[1]) && (pixel[1] < width));
}




/**
 * TESTS
 */

void test_arm() {
	PR2 *brett = new PR2();
	boost::this_thread::sleep(boost::posix_time::seconds(1));

	brett->larm->set_posture(Arm::Posture::side);
	rave::Transform start_pose = brett->larm->get_pose();
	std::cout << "start joints:\n" << brett->larm->get_joint_values() << "\n";
	std::cout << "start pose:\n" << rave_utils::rave_transform_to_mat(start_pose) << "\n";

	brett->larm->set_posture(Arm::Posture::mantis);
	std::cout << "mantis joints:\n" << brett->larm->get_joint_values() << "\n";
	std::cout << "In mantis. Press enter to go back to start pose\n";
	std::cin.ignore();

	brett->larm->set_pose(start_pose);
	std::cout << "end joints:\n" << brett->larm->get_joint_values() << "\n";

	std::cout << "Press enter to quit\n";
	std::cin.ignore();
}

void test_teleop() {
	PR2 *brett = new PR2();
	boost::this_thread::sleep(boost::posix_time::seconds(1));

	brett->rarm->set_posture(Arm::Posture::mantis);
	brett->rarm->teleop();
}

void test_head() {
	PR2 *brett = new PR2();
	boost::this_thread::sleep(boost::posix_time::seconds(1));

	Arm *arm = brett->rarm;
	Head *head = brett->head;

	arm->set_posture(Arm::Posture::mantis);

	mat j = head->get_joint_values();
	std::cout << "head joints: " << j.t() << "\n";

	mat lower, upper;
	head->get_limits(lower, upper);
	std::cout << "lower: " << lower.t() << "\n";
	std::cout << "upper: " << upper.t() << "\n";

	arm->teleop();
	head->look_at(arm->get_pose());
}

void test_sensor() {
	PR2 *brett = new PR2();
	CameraSensor *camera = brett->camera;

	camera->power_on();
	boost::this_thread::sleep(boost::posix_time::seconds(1));
	std::cout << "Getting image...\n";
	cube image = camera->get_image();
	for(int i=0; i < image.n_rows; ++i) {
		for(int j=0; j < image.n_cols; ++j) {
			mat rgb = image.subcube(i,j,0,i,j,2);
//			std::cout << rgb;
		}
	}
	std::cout << "Image size: (" << image.n_rows << ", " << image.n_cols << ", " << image.n_slices << ")\n";

//	std::vector<std::vector<mat> > points_pixels_colors = camera->get_pixels_and_colors()

}

int main(int argc, char* argv[]) {
//	test_arm();
//	test_teleop();
//	test_head();
	test_sensor();
}
