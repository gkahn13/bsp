#include "system/pr2_eih_mapping.h"
#include "sqp/pr2_eih_sqp.h"
#include "../util/logging.h"

#include <boost/program_options.hpp>
namespace po = boost::program_options;

const int T = TIMESTEPS;

class PR2EihMappingNBV : public PR2EihMapping {
public:
	PR2EihMappingNBV(int max_occluded_regions, double max_travel_distance) :
		PR2EihMapping(max_occluded_regions, max_travel_distance) { };
	PR2EihMappingNBV() : PR2EihMappingNBV(1, 0.1) { }

	void next_best_view(int max_iters, const MatrixJ& j_sigma0,
			const std::vector<Gaussian3d>& obj_gaussians, const std::vector<geometry3d::Triangle>& obstacles,
			StdVectorJ& J, StdVectorU& U);

private:
	void sample_joint_trajectory(const VectorJ& j0, StdVectorJ& J, StdVectorU& U);
};

void PR2EihMappingNBV::next_best_view(int max_iters, const MatrixJ& j_sigma0,
			const std::vector<Gaussian3d>& obj_gaussians, const std::vector<geometry3d::Triangle>& obstacles,
			StdVectorJ& J, StdVectorU& U) {
	VectorJ j_t = arm->get_joints();
	const double alpha = 10;

	J = StdVectorJ(T);
	StdVectorJ J_i(T);
	U = StdVectorU(T-1);
	StdVectorU U_i(T-1);
	double min_cost = INFINITY;

	for(int i=0; i < max_iters; ++i) {
		sample_joint_trajectory(j_t, J_i, U_i);
		double cost = sys->cost(J_i, j_sigma0, U_i, obj_gaussians, alpha, obstacles);

		if (cost < min_cost) {
			min_cost = cost;
			J = J_i;
			U = U_i;
		}
	}

	for(const VectorJ& j : J) {
		sim->plot_transform(sim->transform_from_to(arm_sim->fk(j), "base_link", "world"));
	}
}

void PR2EihMappingNBV::sample_joint_trajectory(const VectorJ& j0, StdVectorJ& J, StdVectorU& U) {
	J = StdVectorJ(T);
	U = StdVectorU(T-1);

	J[0] = j0;
	for(int t=0; t < T-1; ) {
		U[t] = (1/float(T))*(M_PI/4.0)*VectorU::Random();
		J[t+1] = sys->dynfunc(J[t], U[t], VectorQ::Zero());

		Matrix4d pose = arm_sim->fk(J[t+1]);
		if (pose(2,3) >= arm_sim->fk(home_joints)(2,3)) {
			++t; // if z position above home pose, for safety
		}
	}
}

void parse_args(int argc, char* argv[],
		bool& pause, int& max_occluded_regions, double& max_travel_distance, int& max_iters) {
	std::vector<double> init_traj_control_vec;
	// Declare the supported options.
	po::options_description desc("Allowed options");
	desc.add_options()
		("help", "produce help message")
		("p", "Pause between stages")
		("m", po::value<int>(&max_occluded_regions)->default_value(1), "Maximum number of occluded regions during BSP")
		("d", po::value<double>(&max_travel_distance)->default_value(1.0), "Maximum distance traveled when executing BSP controls")
		("i", po::value<int>(&max_iters)->default_value(1000), "Maximum iterations for next-best-view")
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

	ROS_INFO_STREAM("Max occluded regions: " << max_occluded_regions);
	ROS_INFO_STREAM("Max travel distance: " << max_travel_distance);
}

int main(int argc, char* argv[]) {
	// initialize ros node
	ros::init(argc, argv, "pr2_eih_nbv");
	log4cxx::LoggerPtr my_logger = log4cxx::Logger::getLogger(ROSCONSOLE_DEFAULT_NAME);
	my_logger->setLevel(ros::console::g_level_lookup[ros::console::levels::Info]);
	ROS_INFO("Starting ROS node");
	ros::spinOnce();

	bool pause;
	int max_occluded_regions;
	double max_travel_distance;
	int max_iters;
	parse_args(argc, argv, pause, max_occluded_regions, max_travel_distance, max_iters);

	PR2EihMappingNBV brett_nbv(max_occluded_regions, max_travel_distance);

	std::vector<Gaussian3d> obj_gaussians;
	std::vector<geometry3d::Triangle> obstacles;

	MatrixJ j_sigma0 = (M_PI/4)*MatrixJ::Identity(); // TODO: never update in MPC
	StdVectorJ J;
	StdVectorU U;

	StdVectorJ grasp_joint_traj, return_grasp_joint_traj;

	ros::Duration(0.5).sleep();

	ROS_INFO("Resetting kinfu and turning on head");
	brett_nbv.reset_kinfu();

	while(!ros::isShuttingDown()) {
		ROS_INFO("Getting occluded regions");
		if (pause) { ROS_INFO("Press enter"); std::cin.ignore(); }
		brett_nbv.get_occluded_regions(obj_gaussians, obstacles);
		brett_nbv.display_gaussian_means(obj_gaussians);
		ros::spinOnce();

		ROS_INFO("Finding next-best-view");
		if (pause) { ROS_INFO("Press enter"); std::cin.ignore(); }
		ros::spinOnce();
		brett_nbv.next_best_view(max_iters, j_sigma0, obj_gaussians, obstacles, J, U);

		ROS_INFO("Displaying next-best-view");
		if (pause) { ROS_INFO("Press enter"); std::cin.ignore(); }
		brett_nbv.display_trajectory(J);

		ROS_INFO("Executing control");
		if (pause) { ROS_INFO("Press enter"); std::cin.ignore(); }
		ros::spinOnce();
		brett_nbv.execute_controls(U);

		ROS_INFO("Checking if there exists a valid grasp trajectory");
		if (brett_nbv.is_valid_grasp_trajectory(grasp_joint_traj, return_grasp_joint_traj)) {
			ROS_INFO("Valid grasp trajectory exists! Execute grasp");
			if (pause) { ROS_INFO("Press enter"); std::cin.ignore(); }

			brett_nbv.execute_grasp_trajectory(grasp_joint_traj, return_grasp_joint_traj);

			ROS_INFO("Resetting kinfu and turning on head");
			brett_nbv.reset_kinfu();
		}

		ros::spinOnce();
		ros::Duration(0.1).sleep();
	}
}
