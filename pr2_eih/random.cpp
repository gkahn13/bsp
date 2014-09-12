#include "system/pr2_eih_mapping.h"
#include "sqp/pr2_eih_sqp.h"
#include "../util/logging.h"

#include <boost/program_options.hpp>
namespace po = boost::program_options;

const int T = TIMESTEPS;

class PR2EihMappingRandom : public PR2EihMapping {
public:
	PR2EihMappingRandom(double max_travel_distance) :
		PR2EihMapping(1, max_travel_distance) { };
	PR2EihMappingRandom() : PR2EihMappingRandom(0.1) { }

	void sample_joint_trajectory(StdVectorJ& J, StdVectorU& U);
};

void PR2EihMappingRandom::sample_joint_trajectory(StdVectorJ& J, StdVectorU& U) {
	VectorJ j0 = arm->get_joints();
	J = StdVectorJ(T);
	U = StdVectorU(T-1);

	J[0] = j0;
	for(int t=0; t < T-1; ) {
		U[t] = (1/float(T))*(M_PI/8.0)*VectorU::Random();
		J[t+1] = sys->dynfunc(J[t], U[t], VectorQ::Zero());

		Matrix4d pose = arm_sim->fk(J[t+1]);
		if (pose(2,3) >= arm_sim->fk(home_joints)(2,3)) {
			++t; // if z position above home pose, for safety
		}
	}
}

void parse_args(int argc, char* argv[],
		bool& pause, double& max_travel_distance) {
	std::vector<double> init_traj_control_vec;
	// Declare the supported options.
	po::options_description desc("Allowed options");
	desc.add_options()
		("help", "produce help message")
		("p", "Pause between stages")
		("d", po::value<double>(&max_travel_distance)->default_value(1.0), "Maximum distance traveled when executing BSP controls")
		;
	try {
		po::variables_map vm;
		po::store(po::parse_command_line(argc, argv, desc, po::command_line_style::unix_style ^ po::command_line_style::allow_short), vm);
		po::notify(vm);

		pause = vm.count("p");

		if (vm.count("help")) {
			std::cout << desc << "\n";
			exit(0);
		}
	} catch (std::exception &e) {
		std::cerr << "error: " << e.what() << "\n";
		exit(0);
	}

	ROS_INFO_STREAM("Max travel distance: " << max_travel_distance);
}

int main(int argc, char* argv[]) {
	// initialize ros node
	ros::init(argc, argv, "pr2_eih_random");
	log4cxx::LoggerPtr my_logger = log4cxx::Logger::getLogger(ROSCONSOLE_DEFAULT_NAME);
	my_logger->setLevel(ros::console::g_level_lookup[ros::console::levels::Info]);
	ROS_INFO("Starting ROS node");
	ros::spinOnce();

	bool pause;
	double max_travel_distance;
	parse_args(argc, argv, pause, max_travel_distance);

	PR2EihMappingRandom brett_random(max_travel_distance);

	StdVectorJ J;
	StdVectorU U;

	StdVectorJ grasp_joint_traj, return_grasp_joint_traj;

	ros::Duration(0.5).sleep();

	ROS_INFO("Resetting kinfu and turning on head");
	brett_nbv.reset_kinfu();

	while(!ros::isShuttingDown()) {
		ROS_INFO("Sampling random trajectory");
		if (pause) { ROS_INFO("Press enter"); std::cin.ignore(); }
		ros::spinOnce();
		brett_random.sample_joint_trajectory(J, U);

		ROS_INFO("Displaying random trajectory");
		if (pause) { ROS_INFO("Press enter"); std::cin.ignore(); }
		brett_random.display_trajectory(J);

		ROS_INFO("Executing controls");
		if (pause) { ROS_INFO("Press enter"); std::cin.ignore(); }
		ros::spinOnce();
		brett_random.execute_controls(U);

		ROS_INFO("Checking if there exists a valid grasp trajectory");
		if (brett_random.is_valid_grasp_trajectory(grasp_joint_traj, return_grasp_joint_traj)) {
			ROS_INFO("Valid grasp trajectory exists! Execute grasp");
			if (pause) { ROS_INFO("Press enter"); std::cin.ignore(); }

			brett_random.execute_grasp_trajectory(grasp_joint_traj, return_grasp_joint_traj);

			ROS_INFO("Resetting kinfu and turning on head");
			brett_nbv.reset_kinfu();
		}

		ros::spinOnce();
		ros::Duration(0.1).sleep();
	}
}
