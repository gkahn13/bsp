#include "../system/pr2-sim.h"
#include "../system/pr2-system.h"
#include "../system/voxel-grid.h"
#include "../geometry/geometry3d.h"
#include "../utils/rave-utils.h"
#include "../utils/pr2-utils.h"
#include "../../util/Timer.h"

#include <openrave-core.h>
namespace rave = OpenRAVE;

TimerCollection tc;

MatrixP setup_particles(rave::EnvironmentBasePtr env) {
	rave::KinBodyPtr table = env->GetKinBody("table");
	rave::KinBody::LinkPtr base = table->GetLink("base");
	rave::Vector extents = base->GetGeometry(0)->GetBoxExtents();

	rave::Vector table_pos = table->GetTransform().trans;
	double x_min, x_max, y_min, y_max, z_min, z_max;
	x_min = table_pos.x - extents.x;
	x_max = table_pos.x + extents.x;
	y_min = table_pos.y - extents.y;
	y_max = table_pos.y + extents.y;
	z_min = table_pos.z + extents.z;
	z_max = table_pos.z + extents.z + .2;

	MatrixP P;

	// uniform
	for(int m=0; m < M_DIM; ++m) {
		P(0,m) = mm_utils::uniform(x_min, x_max);
		P(1,m) = mm_utils::uniform(y_min, y_max);
		P(2,m) = mm_utils::uniform(z_min, z_max);
	}

	// two clumps
//	for(int m=0; m < M_DIM; ++m) {
//		if (m < M_DIM/2) {
//			P(0,m) = mm_utils::uniform(x_min, .7*x_min + .3*x_max);
//			P(1,m) = mm_utils::uniform(y_min, y_max);
//			P(2,m) = mm_utils::uniform(z_min, z_max);
//		} else {
//			P(0,m) = mm_utils::uniform(.3*x_min + .7*x_max, x_max);
//			P(1,m) = mm_utils::uniform(y_min, y_max);
//			P(2,m) = mm_utils::uniform(z_min, z_max);
//		}
//	}

	return P;
}
//
//void test_particle_update() {
//	Vector3d object(3.35, -1.11, 0.8);
//	Arm::ArmType arm_type = Arm::ArmType::right;
//	bool view = true;
//	PR2System sys(object, arm_type, view);
//
//	PR2* brett = sys.get_brett();
//	Arm* arm = sys.get_arm();
//	rave::EnvironmentBasePtr env = brett->get_env();
//	sleep(2);
//
//	arm->set_posture(Arm::Posture::mantis);
//
//	// setup scenario
//	VectorJ j_t_real = arm->get_joint_values(), j_tp1_real;
//	VectorJ j_t = j_t_real, j_tp1; // i.e. no belief == actual
//	VectorU u_t = VectorU::Zero();
//	MatrixP P_t = setup_particles(env), P_tp1;
//
//	LOG_INFO("Origin particles");
//	std::vector<ParticleGaussian> particle_gmm_t;
//	sys.fit_gaussians_to_pf(P_t, particle_gmm_t);
//	sys.display(j_t, particle_gmm_t);
//
//	sys.execute_control_step(j_t_real, j_t, u_t, P_t, j_tp1_real, j_tp1, P_tp1);
//
//	LOG_INFO("New particles")
//	std::vector<ParticleGaussian> particle_gmm_tp1;
//	sys.fit_gaussians_to_pf(P_tp1, particle_gmm_tp1);
//	sys.display(j_tp1, particle_gmm_tp1);
//}
//
//void test_figtree() {
//	Vector3d object(3.35, -1.11, 0.8);
//	Arm::ArmType arm_type = Arm::ArmType::right;
//	bool view = true;
//	PR2System sys(object, arm_type, view);
//
//	PR2* brett = sys.get_brett();
//	Arm* arm = sys.get_arm();
//	rave::EnvironmentBasePtr env = brett->get_env();
//	arm->set_posture(Arm::Posture::mantis);
//	sleep(1);
//
//	MatrixP P = setup_particles(env);
//	std::vector<ParticleGaussian> particle_gmm;
//
//	tc.start("figtree");
//	sys.fit_gaussians_to_pf(P, particle_gmm);
//	tc.stop("figtree");
//
//	std::cout << "particle_gmm size: " << particle_gmm.size() << "\n";
//	for(int i=0; i < particle_gmm.size(); ++i) {
//		std::cout << "mean: " << particle_gmm[i].mean.transpose() << "\n";
//		std::cout << "cov diag: " << particle_gmm[i].cov.diagonal().transpose() << "\n";
//		std::cout << "pct: " << particle_gmm[i].pct << "\n";
//	}
//
//	tc.print_all_elapsed();
//	sys.display(arm->get_joint_values(), particle_gmm);
//}
//
//void test_pr2_system() {
//	Vector3d object(3.35, -1.11, 0.8);
//	Arm::ArmType arm_type = Arm::ArmType::right;
//	bool view = true;
//	PR2System sys(object, arm_type, view);
//
//	PR2* brett = sys.get_brett();
//	Arm* arm = sys.get_arm();
//	rave::EnvironmentBasePtr env = brett->get_env();
//	sleep(2);
//
//	arm->set_posture(Arm::Posture::mantis);
//
//	// setup particles
//	MatrixP P = setup_particles(env);
//
//	// test plotting
//	VectorJ j = arm->get_joint_values();
//
//	std::vector<ParticleGaussian> particle_gmm;
//	particle_gmm.push_back(ParticleGaussian(Vector3d::Zero(), Matrix3d::Identity(), P.leftCols(M_DIM/2), 1));
//	particle_gmm.push_back(ParticleGaussian(Vector3d::Zero(), Matrix3d::Identity(), P.rightCols(M_DIM/2), 1));
//
//	sys.display(j, particle_gmm);
//}
//
//
//double cost(PR2System& sys, const VectorJ& j, const Vector3d& obj, const VectorU& u, const double alpha) {
//	VectorX x_t, x_tp1;
//	MatrixX sigma_t, sigma_tp1;
//
//	sigma_t = .01*MatrixX::Identity();
//	sigma_t.block<3,3>(J_DIM,J_DIM) = 20*Matrix3d::Identity();
//
//	x_t << j, obj;
//	sys.belief_dynamics(x_t, sigma_t, u, alpha, x_tp1, sigma_tp1);
//
//	return sigma_tp1.trace();
//}
//
//VectorU cost_grad(PR2System& sys, const VectorJ& j, const Vector3d& obj, const VectorU& u, const double alpha) {
//	VectorU grad;
//
//	MatrixX sigma = .01*MatrixX::Identity();
//	sigma.block<3,3>(J_DIM,J_DIM) = 20*Matrix3d::Identity();
//
//	double cost_p, cost_m;
//	VectorU u_plus = u, u_minus = u;
//	for(int i=0; i < U_DIM; ++i) {
//		u_plus = u; u_minus = u;
//		u_plus(i) = u(i) + epsilon;
//		u_minus(i) = u(i) - epsilon;
//
//		cost_p = cost(sys, j, obj, u_plus, alpha);
//
//		cost_m = cost(sys, j, obj, u_minus, alpha);
//
//		grad(i) = (cost_p - cost_m) / (2*epsilon);
//	}
//	return grad;
//}
//
//
//void test_camera() {
//	Vector3d object(3.35, -1.11, 0.8);
//	Arm::ArmType arm_type = Arm::ArmType::right;
//	bool view = true;
//	PR2System sys(object, arm_type, view);
//
//	PR2* brett = sys.get_brett();
//	Arm* arm = sys.get_arm();
//	Camera* cam = sys.get_camera();
//	rave::EnvironmentBasePtr env = brett->get_env();
//
//	MatrixP P = setup_particles(env);
//
//	arm->set_posture(Arm::Posture::mantis);
//	sleep(2);
//
//	VectorJ j = arm->get_joint_values();
//	StdVector3d pcl = cam->get_pcl(j);
//	cam->plot_pcl(pcl);
//
////	std::cout << "Displaying env mesh. Teleop and then get FOV\n";
//
//	arm->teleop();
//
//	std::vector<std::vector<Beam3d> > beams = cam->get_beams(j, pcl);
//	rave_utils::clear_plots();
//	cam->plot_fov(beams);
//
//	std::vector<Triangle3d> border = cam->get_border(beams);
//
////	tc.start("sd");
////	for(int m=0; m < M_DIM; ++m) {
////		double sd = cam->signed_distance(P.col(m), beams, border);
////	}
////	tc.stop("sd");
//
//	Vector3d p(3.35, -2.5, 0.8);
//	rave_utils::plot_point(env, p, Vector3d(1,0,0));
//
//	while(true) {
//		VectorU grad = cost_grad(sys, j, p, VectorU::Zero(), .01);
//		std::cout << "grad:\n" << grad << "\n";
//
//		if (grad.norm() < epsilon) {
//			std::cout << "crap\n";
//			exit(0);
//		}
//
//		VectorJ j_new = j - (M_PI/32)*grad/grad.norm();
////		arm->set_joint_values(j_new);
//
//		beams = cam->get_beams(j_new, pcl);
//		rave_utils::clear_plots();
//		rave_utils::plot_transform(env, rave_utils::eigen_to_rave(arm->get_pose(j_new)));
//		cam->plot_fov(beams);
//		rave_utils::plot_point(env, p, Vector3d(1,0,0));
//
//		j = j_new;
//
//		std::cin.ignore();
//	}
//
//	tc.print_all_elapsed();
//	std::cout << "Displaying current fov beams. Press enter to exit\n";
//	std::cin.ignore();
//}

void test_fk() {
	Vector3d object(3.35, -1.11, 0.8);
	Arm::ArmType arm_type = Arm::ArmType::right;
	bool view = true;
	PR2System sys(object, arm_type, view);

	PR2* brett = sys.get_brett();
	Arm* arm = sys.get_arm();
	Camera* cam = sys.get_camera();
	rave::EnvironmentBasePtr env = brett->get_env();

	arm->set_posture(Arm::Posture::mantis);
	sleep(1);

	VectorJ j = arm->get_joint_values();

	Matrix4d actual_arm_pose = rave_utils::transform_from_to(brett->get_robot(), Matrix4d::Identity(), "r_gripper_tool_frame", "world");
	Matrix4d fk_arm_pose = arm->get_pose(j);

	std::cout << "actual_arm_pose:\n" << actual_arm_pose << "\n\n";
	std::cout << "fk_arm_pose:\n" << fk_arm_pose << "\n\n";

	Matrix4d actual_cam_pose = rave_utils::rave_to_eigen(cam->get_sensor()->GetTransform());
	Matrix4d fk_cam_pose = cam->get_pose(j);

	std::cout << "actual_cam_pose:\n" << actual_cam_pose << "\n\n";
	std::cout << "fk_cam_pose:\n" << fk_cam_pose << "\n\n";
}

void test_voxel_grid() {
	Vector3d table(3.5, -1.2, 0.74);
	Vector3d object = table + Vector3d(0, 0, .4);
	Arm::ArmType arm_type = Arm::ArmType::right;
	bool view = true;
	PR2System sys(object, arm_type, view);

	PR2* brett = sys.get_brett();
	Arm* arm = sys.get_arm();
	Camera* cam = sys.get_camera();
	rave::EnvironmentBasePtr env = brett->get_env();

	arm->set_posture(Arm::Posture::mantis);
	sleep(1);

	rave_utils::plot_point(env, object, Vector3d(1,0,0), .05);

	Vector3d pos_center = table;
	double x_height = 1, y_height = 2, z_height = 1;
	int resolution = 100;
	VoxelGrid vgrid(pos_center, x_height, y_height, z_height, resolution);

	StdVector3d pcl = cam->get_pcl(arm->get_joint_values());

//	for(int i=0; i < pcl.size(); ++i) {
//		rave_utils::plot_point(env, pcl[i], Vector3d(1,0,0), .005);
//	}

	std::cout << "update_TSDF\n";
	vgrid.update_TSDF(pcl);

	std::cout << "update_ODF\n";
	tc.start("update_ODF");
	vgrid.update_ODF(object, env);
	tc.stop("update_ODF");

	StdVector3d obstacles = vgrid.get_obstacles();
	Matrix<double,H_SUB,W_SUB> zbuffer = cam->get_zbuffer(arm->get_joint_values(), pcl); // TODO: should be obstacles

	tc.start("sd");
	vgrid.signed_distance_complete(cam, zbuffer, cam->get_pose(arm->get_joint_values()));
	tc.stop("sd");

	vgrid.plot_TSDF(env);
//	vgrid.plot_ODF(env);
	vgrid.plot_FOV(env, cam, zbuffer, cam->get_pose(arm->get_joint_values()));

	tc.print_all_elapsed();

	std::cout << "Press enter to exit\n";
	std::cin.ignore();

}

int main() {
//	test_particle_update();
//	test_figtree();
//	test_pr2_system();
//	test_camera();
//	test_fk();
	test_voxel_grid();
	return 0;
}
