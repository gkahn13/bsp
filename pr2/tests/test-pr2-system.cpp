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
		P(0,m) = pr2_utils::uniform(x_min, x_max);
		P(1,m) = pr2_utils::uniform(y_min, y_max);
		P(2,m) = pr2_utils::uniform(z_min, z_max);
	}

	// two clumps
//	for(int m=0; m < M_DIM; ++m) {
//		if (m < M_DIM/2) {
//			P(0,m) = pr2_utils::uniform(x_min, .7*x_min + .3*x_max);
//			P(1,m) = pr2_utils::uniform(y_min, y_max);
//			P(2,m) = pr2_utils::uniform(z_min, z_max);
//		} else {
//			P(0,m) = pr2_utils::uniform(.3*x_min + .7*x_max, x_max);
//			P(1,m) = pr2_utils::uniform(y_min, y_max);
//			P(2,m) = pr2_utils::uniform(z_min, z_max);
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
//	StdVector3d pc = cam->get_pc(j);
//	cam->plot_pc(pc);
//
////	std::cout << "Displaying env mesh. Teleop and then get FOV\n";
//
//	arm->teleop();
//
//	std::vector<std::vector<Beam3d> > beams = cam->get_beams(j, pc);
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
//		beams = cam->get_beams(j_new, pc);
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

void test_ik() {
	Vector3d object(3.35, -1.11, 0.8);
	Arm::ArmType arm_type = Arm::ArmType::right;
	bool view = true;
	PR2System sys(object, arm_type, view);

	PR2* brett = sys.get_brett();
	Arm* arm = sys.get_arm();
	rave::EnvironmentBasePtr env = brett->get_env();

	arm->set_posture(Arm::Posture::mantis);
	sleep(1);

	for(int iter=0; iter < 10; ++iter) {
		std::cout << "iter: " << iter << "\n";
		Vector3d pos = object + .2*Vector3d::Random();
		rave_utils::plot_point(env, pos, Vector3d(1,0,0), .05);

		VectorJ j = arm->get_joint_values();
		bool success = arm->ik(pos, j);
		std::cout << "success: " << success << "\n";
		arm->set_joint_values(j);
		std::cin.ignore();
	}
}

void test_greedy() {
	srand(time(0));

	Vector3d table(3.5, -1.2, 0.74);
	Vector3d object = table + Vector3d(0, -.2, .05);
//	Vector3d object = table + Vector3d(0, .1, -.2);
	Arm::ArmType arm_type = Arm::ArmType::right;
	bool view = true;
	PR2System sys(object, arm_type, view);

	PR2* brett = sys.get_brett();
	Arm* arm = sys.get_arm();
	Camera* cam = sys.get_camera();
	rave::EnvironmentBasePtr env = brett->get_env();

	arm->set_posture(Arm::Posture::mantis);
	sleep(1);

	rave_utils::plot_transform(env, rave_utils::eigen_to_rave(cam->get_pose(arm->get_joint_values())));

	Vector3d pos_center = table;
	double x_height = 1.5, y_height = 2, z_height = 1;
	int resolution = 100;
	VoxelGrid vgrid(pos_center, x_height, y_height, z_height, resolution, cam->get_pose(arm->get_joint_values()));

	StdVector3d pc = cam->get_pc(arm->get_joint_values());
	Matrix<double,HEIGHT_FULL,WIDTH_FULL> full_zbuffer = cam->get_full_zbuffer(arm->get_joint_values(), pc);

	std::cout << "update_TSDF\n";
	vgrid.update(pc, full_zbuffer, cam->get_pose(arm->get_joint_values()));

	Matrix<double,H_SUB,W_SUB> zbuffer = cam->get_zbuffer(arm->get_joint_values(), pc); // TODO: should be obstacles

	Vector3d lower = pos_center - Vector3d(x_height/2., y_height/2., z_height/2.);
	Vector3d upper = pos_center + Vector3d(x_height/2., y_height/2., z_height/2.);
	for(int iter=0; iter < 10; ++iter) {
		rave_utils::clear_plots();
		vgrid.plot_TSDF(env);
		vgrid.plot_FOV(env, cam, zbuffer, cam->get_pose(arm->get_joint_values()));

		std::cout << "\niter: " << iter << "\n";
		for(int i=0; i < 3; ++i) {
			object(i) = pr2_utils::uniform(lower(i), upper(i));
		}

		std::cout << "update_ODF\n";
		tc.start("update_ODF");
		Cube ODF = vgrid.get_ODF(object);
		tc.stop("update_ODF");

		//	std::cout << "zbuffer\n" << zbuffer << "\n";

		tc.start("sd_complete");
		double sd_complete = vgrid.signed_distance_complete(object, ODF, cam, zbuffer, cam->get_pose(arm->get_joint_values()));
		tc.stop("sd_complete");

		tc.start("sd_greedy");
		double sd_greedy = vgrid.signed_distance_greedy(object, ODF, cam, zbuffer, cam->get_pose(arm->get_joint_values()));
		tc.stop("sd_greedy");

		std::cout << "sd_complete: " << sd_complete << "\n";
		std::cout << "sd_greedy: " << sd_greedy << "\n";

		tc.print_all_elapsed();
		tc.clear_all();

		rave_utils::plot_point(env, object, Vector3d(1,0,0), .05);
		std::cin.ignore();
	}


	std::cout << "Press enter to exit\n";
	std::cin.ignore();
}

VectorJ sd_gradient(VectorJ j, const Vector3d& object, const Cube& ODF, Camera* cam, VoxelGrid* vgrid) {
	VectorJ grad = VectorJ::Zero();
	const double step = 1e-5;
	const double alpha = 1;

	Matrix<double,H_SUB,W_SUB> zbuffer = vgrid->get_zbuffer(cam->get_pose(j));
	Matrix4d cam_pose = cam->get_pose(j), sd_pose_in_world = Matrix4d::Identity(), sd_pose_in_cam;
	sd_pose_in_world.block<3,1>(0,3) = vgrid->signed_distance_greedy_voxel_center(object, ODF, cam, zbuffer, cam_pose);
	sd_pose_in_cam = cam_pose.inverse()*sd_pose_in_world;
	bool obj_in_fov = cam->is_in_fov(object, zbuffer, cam_pose);
	double sd_sign = (obj_in_fov) ? -1 : 1;

	Vector3d v = object - sd_pose_in_world.block<3,1>(0,3);

	if (obj_in_fov) {
		std::cout << "obj is in fov\n";
	} else {
		std::cout << "obj is not in fov\n";
	}

	Matrix4d object_pose = Matrix4d::Identity();
	object_pose.block<3,1>(0,3) = object;
	Vector3d z(0,0,1);

	Vector3d object_in_cam_m, object_in_cam_p;
	Matrix4d cam_pose_m, cam_pose_p;
	Vector3d sd_in_world_m, sd_in_world_p, voxel_m, voxel_p;
	Vector3d u_m, u_p;
	double sd_m, sd_p, costheta_m, costheta_p;
	VectorJ j_m, j_p;
	double sd_sigmoid_m, sd_sigmoid_p;
	for(int i=0; i < J_DIM; ++i) {
		j_m = j;
		j_p = j;
		j_m(i) = j(i) - step;
		j_p(i) = j(i) + step;

		// find sd
		cam_pose_m = cam->get_pose(j_m);
		sd_in_world_m = (cam_pose_m*sd_pose_in_cam).block<3,1>(0,3);
		voxel_m = vgrid->exact_voxel_from_point(sd_in_world_m);
		sd_m = sd_sign*ODF.trilinear_interpolation(voxel_m);
//		// find cos(theta)
//		u_m = sd_in_world_m - cam_pose_m.block<3,1>(0,3);
//		costheta_m = (obj_in_fov) ? 0 : .1*(u_m.dot(v)) / (u_m.norm()*v.norm()) + .45;
		// find cos(theta) between camera center axis and object
		object_in_cam_m = (cam_pose_m.inverse()*object_pose).block<3,1>(0,3);
		costheta_m = (object_in_cam_m.dot(z)) / (object_in_cam_m.norm()*z.norm());
		// calculate delta
//		sd_sigmoid_m = 1.0 - 1.0/(1.0+exp(-alpha*sd_m));
		sd_sigmoid_m = costheta_m;

		// find sd
		cam_pose_p = cam->get_pose(j_p);
		sd_in_world_p = (cam_pose_p*sd_pose_in_cam).block<3,1>(0,3);
		voxel_p = vgrid->exact_voxel_from_point(sd_in_world_p);
		sd_p = sd_sign*ODF.trilinear_interpolation(voxel_p);
//		// find cos(theta)
//		u_p = sd_in_world_p - cam_pose_p.block<3,1>(0,3);
//		costheta_p = (obj_in_fov) ? 0 : .1*(u_p.dot(v)) / (u_p.norm()*v.norm()) + .45;
		// find cos(theta) between camera center axis and object
		object_in_cam_p = (cam_pose_p.inverse()*object_pose).block<3,1>(0,3);
		costheta_p = (object_in_cam_p.dot(z)) / (object_in_cam_p.norm()*z.norm());
		// calculate delta
//		sd_sigmoid_p = 1.0 - 1.0/(1.0+exp(-alpha*sd_p));
		sd_sigmoid_p = costheta_p;


		grad(i) = (sd_sigmoid_p - sd_sigmoid_m) / (2*step);
	}

	return grad;
}

void test_gradient() {
//	srand(time(0));

	Vector3d table(3.5, -1.2, 0.74);
	Vector3d object = table + Vector3d(0, .5, .05);
//	Vector3d object = table + Vector3d(.1, -.1, -.1);
	Arm::ArmType arm_type = Arm::ArmType::right;
	bool view = true;
	PR2System sys(object, arm_type, view);

	PR2* brett = sys.get_brett();
	Arm* arm = sys.get_arm();
	Camera* cam = sys.get_camera();
	VoxelGrid* vgrid = sys.get_voxel_grid();
	rave::EnvironmentBasePtr env = brett->get_env();

	arm->set_posture(Arm::Posture::mantis);
	sleep(1);

	VectorJ j = arm->get_joint_values();

	while(true) {
		std::cout << "Updating internal TSDF and kinfu\n";
		StdVector3d pc = cam->get_pc(j);
		Matrix<double,HEIGHT_FULL,WIDTH_FULL> full_zbuffer = cam->get_full_zbuffer(j, pc);
		Matrix4d cam_pose = cam->get_pose(j);
		sys.update(pc, full_zbuffer, cam_pose);

		std::cout << "Press enter:\n";
		std::cin.ignore();
		rave_utils::clear_plots();

		std::cout << "get_ODF\n";
		Cube ODF = sys.get_ODF(object);

		arm->set_joint_values(j);
		tc.start("sd_gradient");
		VectorJ grad = sd_gradient(j, object, ODF, cam, vgrid);
		tc.stop("sd_gradient");
		tc.print_all_elapsed();
		tc.clear_all();

		std::cout << "grad\n" << grad.transpose() << "\n";

		rave_utils::plot_point(env, object, Vector3d(1,0,0), .05);
		vgrid->plot_TSDF(env);
		vgrid->plot_FOV(env, cam, vgrid->get_zbuffer(cam->get_pose(j)), cam->get_pose(j));
		vgrid->plot_kinfu_tsdf(env);

		if (grad.norm() > epsilon) {
			grad /= grad.norm();
		}
		j = j + (M_PI/16)*(grad);
	}
}

void test_raycaster() {
	srand(time(0));

	Vector3d table(3.5, -1.2, 0.74);
	Vector3d object = table + Vector3d(0, .5, .05);
//	Vector3d object = table + Vector3d(.1, -.1, -.1);
	Arm::ArmType arm_type = Arm::ArmType::right;
	bool view = true;
	PR2System sys(object, arm_type, view);

	PR2* brett = sys.get_brett();
	Arm* arm = sys.get_arm();
	Camera* cam = sys.get_camera();
	VoxelGrid* vgrid = sys.get_voxel_grid();
	rave::EnvironmentBasePtr env = brett->get_env();

//	arm->set_posture(Arm::Posture::mantis);
	sleep(1);

	VectorJ j = arm->get_joint_values();

	std::cout << "Getting point cloud\n";
	StdVector3d pc = cam->get_pc(j);

	std::cout << "Getting full zbuffer\n";
	Matrix<double,HEIGHT_FULL,WIDTH_FULL> full_zbuffer = cam->get_full_zbuffer(j, pc);

	std::cout << "Updating\n";
	vgrid->update(pc, full_zbuffer, cam->get_pose(j));

	std::cout << "Visualize tsdf\n";
	vgrid->plot_kinfu_tsdf(env);

	arm->teleop();
	j = arm->get_joint_values();

	rave_utils::plot_transform(env, rave_utils::eigen_to_rave(cam->get_pose(j)));

	std::cout << "get_zbuffer using RayCaster\n";
	tc.start("get_zbuffer");
	Matrix<double,H_SUB,W_SUB> zbuffer = vgrid->get_zbuffer(cam->get_pose(j));
	tc.stop("get_zbuffer");

	std::cout << "plot_FOV_full\n";
	vgrid->plot_FOV(env, cam, zbuffer, cam->get_pose(j));

//	std::cout << "update_TSDF\n";
//	vgrid->update_TSDF(pc);
//
//	std::cout << "plot_TSDF\n";
//	vgrid->plot_TSDF(env);
//
//	std::cout << "plot_FOV_full\n";
//	vgrid->plot_FOV_full(env, cam, cam->get_full_zbuffer(j, pc), cam->get_pose(j));

//	std::cout << "get_zbuffer using projection\n";
//	tc.start("get_zbuffer projection");
//	cam->get_zbuffer(j, pc);
//	tc.stop("get_zbuffer projection");

	tc.print_all_elapsed();
	std::cout << "Press enter to exit\n";
	std::cin.ignore();
}

void test_tri_interp() {
	srand(time(0));

	Vector3d table(3.5, -1.2, 0.74);
//	Vector3d object = table + Vector3d(0, .5, .05);
	Vector3d object = table + Vector3d(.1, -.1, -.1);
	Arm::ArmType arm_type = Arm::ArmType::right;
	bool view = true;
	PR2System sys(object, arm_type, view);

	PR2* brett = sys.get_brett();
	Arm* arm = sys.get_arm();
	Camera* cam = sys.get_camera();
	VoxelGrid* vgrid = sys.get_voxel_grid();
	rave::EnvironmentBasePtr env = brett->get_env();

	arm->set_posture(Arm::Posture::mantis);
	sleep(1);

	VectorJ j = arm->get_joint_values();

	std::cout << "Updating internal TSDF and kinfu\n";
	StdVector3d pc = cam->get_pc(j);
	Matrix<double,HEIGHT_FULL,WIDTH_FULL> full_zbuffer = cam->get_full_zbuffer(j, pc);
	Matrix4d cam_pose = cam->get_pose(j);
	sys.update(pc, full_zbuffer, cam_pose);

	std::cout << "get_ODF\n";
	Cube ODF = sys.get_ODF(object);

	Vector3d voxel(1, 1, 1);
	double od;

	od = ODF.trilinear_interpolation(voxel);
	std::cout << od << "\n";

	od = ODF.trilinear_interpolation(voxel + Vector3d(.1,0,0));
	std::cout << od << "\n";

	od = ODF.trilinear_interpolation(voxel + Vector3d(0,.1,0));
	std::cout << od << "\n";

	od = ODF.trilinear_interpolation(voxel + Vector3d(0,0,.1));
	std::cout << od << "\n";

	od = ODF.trilinear_interpolation(Vector3d(20,20,20) + Vector3d(0,0,.1));
	std::cout << od << "\n";

	od = ODF.trilinear_interpolation(Vector3d(20,20,20) + Vector3d(0,0,-.1));
	std::cout << od << "\n";
}

void test_distance_to_TSDF() {
	Vector3d table(3.5, -1.2, 0.74);
//	Vector3d object = table + Vector3d(0, .5, .05);
	Vector3d object = table + Vector3d(.1, -.1, -.1);
	Arm::ArmType arm_type = Arm::ArmType::right;
	bool view = true;
	PR2System sys(object, arm_type, view);

	PR2* brett = sys.get_brett();
	Arm* arm = sys.get_arm();
	Camera* cam = sys.get_camera();
	VoxelGrid* vgrid = sys.get_voxel_grid();
	rave::EnvironmentBasePtr env = brett->get_env();

	arm->set_posture(Arm::Posture::mantis);
	sleep(1);

	VectorJ j = arm->get_joint_values();

	std::cout << "Updating internal TSDF and kinfu\n";
	StdVector3d pc = cam->get_pc(j);
	Matrix<double,HEIGHT_FULL,WIDTH_FULL> full_zbuffer = cam->get_full_zbuffer(j, pc);
	Matrix4d cam_pose = cam->get_pose(j);
	sys.update(pc, full_zbuffer, cam_pose);

	rave_utils::plot_point(env, object, Vector3d(0,0,1), .05);
	vgrid->plot_TSDF(env);

	while(true) {
		arm->teleop();

		double dist = vgrid->distance_to_TSDF(cam->get_position(arm->get_joint_values()));
		std::cout << "dist: " << dist << "\n";
	}
}

int main() {
//	test_particle_update();
//	test_figtree();
//	test_pr2_system();
//	test_camera();
//	test_fk();
//	test_ik();
//	test_voxel_grid();
//	test_greedy();
//	test_gradient();
//	test_raycaster();
//	test_tri_interp();
	test_distance_to_TSDF();
	return 0;
}
