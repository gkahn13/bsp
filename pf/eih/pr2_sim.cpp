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
	std::cout << "Initializing OpenRAVE\n";
	rave::RaveInitialize(true, rave::Level_Info);
	env = rave::RaveCreateEnvironment();
	std::cout << "Loading environment: " << env_file << "\n";
	env->Load(env_file);

	robot = env->GetRobot(robot_name);

//	robot->SetActiveManipulator("leftarm");
//	rave::RobotBase::ManipulatorPtr manip = robot->GetActiveManipulator();
//	if (!manip->GetIkSolver()) {
//		std::cout << "No IK solver!\n";
//		exit(0);
//	}
//
//	rave::ModuleBasePtr ikfast = rave::RaveCreateModule(env,"ikfast");
//	env->Add(ikfast,true,"");
//	std::stringstream ssin,ssout;
//	ssin << "LoadIKFastSolver " << robot->GetName() << " " << (int)rave::IKP_Transform6D;
//	// if necessary, add free inc for degrees of freedom
//	//ssin << " " << 0.04f;
//	// get the active manipulator
//	robot->SetActiveManipulator("leftarm");
//	rave::RobotBase::ManipulatorPtr pmanip = robot->GetActiveManipulator();
//	std::cout << "Active manipulator: " << pmanip->GetName() << "\n";
//	if( !ikfast->SendCommand(ssout,ssin) ) {
//		std::cout << ssin.str() << "\n";
//		RAVELOG_ERROR("failed to load iksolver\n");
//		exit(0);
//	}

	if (view) {
		viewer_thread = boost::shared_ptr<boost::thread>(new boost::thread(boost::bind(SetViewer, env, "qtcoin")));
//		viewer_thread->join(); // TODO: temp
	}
//	env->Destroy(); // destroy

	larm = new Arm(robot, Arm::ArmType::left);
}



/**
 * Arm Constructors
 */

Arm::Arm(rave::RobotBasePtr robot, ArmType arm_type) {
	init(robot, arm_type);
}

/**
 * Arm Initializers
 */

void Arm::init(rave::RobotBasePtr robot, ArmType arm_type) {
	this->robot = robot;
	this->arm_type = arm_type;

	std::string manip_name;
	if (arm_type == ArmType::left) {
		manip_name = "leftarm";
	} else {
		manip_name = "rightarm";
	}

	manip = robot->GetManipulator(manip_name);
	arm_indices = manip->GetArmIndices();
	num_joints = arm_indices.size();

//	// load inverse kinematics using ikfast
//	rave::EnvironmentBasePtr env = robot->GetEnv();
//	rave::ModuleBasePtr ikfast = rave::RaveCreateModule(env,"ikfast");
//	env->Add(ikfast,true,"");
//	std::stringstream ssin,ssout;
//	std::vector<double> solution;
//	ssin << "LoadIKFastSolver " << robot->GetName() << " " << (int)rave::IKP_Transform6D;
//	if( !ikfast->SendCommand(ssout,ssin) ) {
//		RAVELOG_ERROR("failed to load iksolver\n");
//	}
//	if( !manip->GetIkSolver()) {
//		env->Destroy();
//		exit(0);
//	}
}

/**
 * Arm public methods
 */

mat Arm::get_joint_values() {
	std::vector<double> joint_values;
	manip->GetArmDOFValues(joint_values);
	return conv_to<mat>::from(joint_values);
}

mat Arm::get_pose() {
	return rave_utils::rave_transform_to_mat(manip->GetEndEffectorTransform());
}

void Arm::get_lower_limits(mat& lower, mat& upper) {
	std::vector<double> lower_vec(num_joints), upper_vec(num_joints);
	robot->GetActiveDOFLimits(lower_vec, upper_vec);

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

void Arm::set_joint_values(mat &joint_values) {
	std::vector<double> joint_values_vec = conv_to<std::vector<double>>::from(joint_values);
	robot->SetDOFValues(joint_values_vec, rave::KinBody::CheckLimitsAction::CLA_CheckLimits, arm_indices);
//	robot->SetDOFValues(joint_values_vec, rave::KinBody::CheckLimitsAction::CLA_Nothing, arm_indices);
}

void Arm::set_pose(mat &pose, std::string ref_frame) {
	std::vector<double> joint_values_vec(num_joints);
	rave_utils::cart_to_joint(manip, pose, ref_frame, "end_effector", joint_values_vec);

	mat joint_values = conv_to<mat>::from(joint_values_vec);
	set_joint_values(joint_values);
}

/**
 * For testing only
 */

int main(int argc, char* argv[]) {
	PR2 *brett = new PR2();
	boost::this_thread::sleep(boost::posix_time::seconds(1));

	brett->larm->set_posture(Arm::Posture::side);
	mat start_pose = brett->larm->get_pose();
	std::cout << "start joints:\n" << brett->larm->get_joint_values() << "\n";
	std::cout << "start pose:\n" << start_pose << "\n";

	brett->larm->set_posture(Arm::Posture::mantis);
	std::cout << "mantis joints:\n" << brett->larm->get_joint_values() << "\n";
	std::cout << "In mantis. Press enter to go back to start pose\n";
	std::cin.ignore();

	brett->larm->set_pose(start_pose);
	std::cout << "end joints:\n" << brett->larm->get_joint_values() << "\n";

//	rave::Transform a, b, c;
//	a.identity();
//	b.identity();
//
//	c = a*b;
//	rave::RaveTransformMatrix<double> c_mat(c);
//	std::cout << "c:\n" << c_mat << "\n";

	std::cout << "Press enter to quit\n";
	std::cin.ignore();

	return 0;
}
