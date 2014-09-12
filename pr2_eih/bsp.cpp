#include "system/pr2_eih_mapping.h"
#include "sqp/pr2_eih_sqp.h"
#include "../util/logging.h"

#include <boost/program_options.hpp>
namespace po = boost::program_options;

const int T = TIMESTEPS;

class PR2EihMappingBSP : public PR2EihMapping {
public:
	PR2EihMappingBSP(int max_occluded_regions, double max_travel_distance) :
		PR2EihMapping(max_occluded_regions, max_travel_distance) { };
	PR2EihMappingBSP() : PR2EihMappingBSP(1, 0.1) { }

	void initialize_trajectory(StdVectorJ& J, StdVectorU& U,
			const Vector3d& init_traj_control, const std::vector<Gaussian3d>& obj_gaussians);
	void bsp(StdVectorJ& J, StdVectorU& U, const MatrixJ& j_sigma0,
			const std::vector<Gaussian3d>& obj_gaussians, const std::vector<geometry3d::Triangle>& obstacles);

private:
	// sqp solver
	pr2_eih_sqp::PR2EihSqp pr2_eih_bsp;

};

void PR2EihMappingBSP::initialize_trajectory(StdVectorJ& J, StdVectorU& U,
		const Vector3d& init_traj_control, const std::vector<Gaussian3d>& obj_gaussians) {
	ROS_INFO_STREAM("Delta position for each timestep is: " << init_traj_control.transpose());

	Vector3d avg_obj_mean = Vector3d::Zero();
	double num_objs = obj_gaussians.size();
	for(const Gaussian3d& obj_gaussian : obj_gaussians) {
		avg_obj_mean += (1/num_objs)*obj_gaussian.mean;
	}

	J[0] = arm->get_joints();
	Vector3d start_position = arm_sim->fk(J[0]).block<3,1>(0,3);
	Vector3d next_position;
	VectorJ next_joints;
	for(int t=0; t < T-1; ++t) {
		next_position = start_position + (t+1)*init_traj_control;
		if (arm_sim->ik_lookat(next_position, avg_obj_mean, next_joints)) {
			U[t] = next_joints - J[t];
		} else {
			U[t] = VectorU::Zero();
		}

		J[t+1] = sys->dynfunc(J[t], U[t], VectorQ::Zero(), true);
		U[t] = (J[t+1] - J[t])/double(DT);
	}
}

void PR2EihMappingBSP::bsp(StdVectorJ& J, StdVectorU& U, const MatrixJ& j_sigma0,
		const std::vector<Gaussian3d>& obj_gaussians, const std::vector<geometry3d::Triangle>& obstacles) {
	publish_to_logger("start bsp");

	ROS_INFO_STREAM("Number of gaussians: " << obj_gaussians.size());
	ROS_INFO_STREAM("Number of obstacles: " << obstacles.size());

	// plan
	pr2_eih_bsp.collocation(J, U, j_sigma0, obj_gaussians, obstacles, *sys, true);

	// reintegrate, just in case
	for(int t=0; t < T-1; ++t) {
		J[t+1] = sys->dynfunc(J[t], U[t], VectorQ::Zero(), true);
	}

	publish_to_logger("end bsp");
}

void parse_args(int argc, char* argv[],
		bool& pause, int& max_occluded_regions, Vector3d& init_traj_control, double& max_travel_distance) {
	std::vector<double> init_traj_control_vec;
	// Declare the supported options.
	po::options_description desc("Allowed options");
	desc.add_options()
		("help", "produce help message")
		("p", "Pause between stages")
		("m", po::value<int>(&max_occluded_regions)->default_value(1), "Maximum number of occluded regions during BSP")
		("init", po::value<std::vector<double> >(&init_traj_control_vec)->multitoken(), "Delta position (x,y,z) by which to initialize BSP trajectory")
		("d", po::value<double>(&max_travel_distance)->default_value(1.0), "Maximum distance traveled when executing BSP controls")
		;
	try {
		po::variables_map vm;
		po::store(po::parse_command_line(argc, argv, desc, po::command_line_style::unix_style ^ po::command_line_style::allow_short), vm);
		po::notify(vm);

		pause = vm.count("p");

		if (vm.count("init")) {
			init_traj_control = Vector3d(init_traj_control_vec.data());
		} else {
			init_traj_control << .05, 0, .05;
		}

		if (vm.count("help")) {
			std::cout << desc << "\n";
			exit(0);
		}
	} catch (std::exception &e) {
		std::cerr << "error: " << e.what() << "\n";
		exit(0);
	}

	ROS_INFO_STREAM("Max occluded regions: " << max_occluded_regions);
	ROS_INFO_STREAM("Init traj control: " << init_traj_control.transpose());
	ROS_INFO_STREAM("Max travel distance: " << max_travel_distance);
}

int main(int argc, char* argv[]) {
	// initialize ros node
	ros::init(argc, argv, "pr2_eih_brett");
	log4cxx::LoggerPtr my_logger = log4cxx::Logger::getLogger(ROSCONSOLE_DEFAULT_NAME);
	my_logger->setLevel(ros::console::g_level_lookup[ros::console::levels::Info]);
	ROS_INFO("Starting ROS node");
	ros::spinOnce();

	bool pause;
	int max_occluded_regions;
	Vector3d init_traj_control;
	double max_travel_distance;
	parse_args(argc, argv, pause, max_occluded_regions, init_traj_control, max_travel_distance);

	PR2EihMappingBSP brett_bsp(max_occluded_regions, max_travel_distance);

	std::vector<Gaussian3d> obj_gaussians;
	std::vector<geometry3d::Triangle> obstacles;

	MatrixJ j_sigma0 = (M_PI/4)*MatrixJ::Identity(); // TODO: never update in MPC
	StdVectorJ J(T);
	StdVectorU U(T-1);

	StdVectorJ grasp_joint_traj, return_grasp_joint_traj;

	ros::Duration(0.5).sleep();

	ROS_INFO("Resetting kinfu and turning on head");
	brett_bsp.reset_kinfu();

	while(!ros::isShuttingDown()) {
		ROS_INFO("Getting occluded regions");
		if (pause) { ROS_INFO("Press enter"); std::cin.ignore(); }
		brett_bsp.get_occluded_regions(obj_gaussians, obstacles);
		brett_bsp.display_gaussian_means(obj_gaussians);
		ros::spinOnce();

		ROS_INFO("Initializing trajectory");
		if (pause) { ROS_INFO("Press enter"); std::cin.ignore(); }
		ros::spinOnce();
		brett_bsp.initialize_trajectory(J, U, init_traj_control, obj_gaussians);
		brett_bsp.display_trajectory(J);

		ROS_INFO("Optimizing trajectory with bsp");
		if (pause) { ROS_INFO("Press enter"); std::cin.ignore(); }
		ros::spinOnce();
		brett_bsp.bsp(J, U, j_sigma0, obj_gaussians, obstacles);

		ROS_INFO("Displaying bsp trajectory");
		if (pause) { ROS_INFO("Press enter"); std::cin.ignore(); }
		brett_bsp.display_trajectory(J);

		ROS_INFO("Executing control");
		if (pause) { ROS_INFO("Press enter"); std::cin.ignore(); }
		ros::spinOnce();
		brett_bsp.execute_controls(U);

		ROS_INFO("Checking if there exists a valid grasp trajectory");
		ros::spinOnce();
		ros::Duration(0.5).sleep();
		ros::spinOnce();
		if (brett_bsp.is_valid_grasp_trajectory(grasp_joint_traj, return_grasp_joint_traj)) {
			ROS_INFO("Valid grasp trajectory exists! Execute grasp");
			if (pause) { ROS_INFO("Press enter"); std::cin.ignore(); }

			brett_bsp.execute_grasp_trajectory(grasp_joint_traj, return_grasp_joint_traj);

			ROS_INFO("Resetting kinfu and turning on head");
			brett_bsp.reset_kinfu();
		}

		ros::spinOnce();
		ros::Duration(0.1).sleep();
	}

}
