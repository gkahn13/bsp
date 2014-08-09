#include "../system/pr2-sim.h"
#include "../system/pr2-system.h"
#include "../system/voxel-grid.h"
#include "../geometry/geometry3d.h"
#include "../utils/rave-utils.h"
#include "../utils/pr2-utils.h"
#include "../../util/Timer.h"

#include <openrave-core.h>
namespace rave = OpenRAVE;

StdVector3d get_voxel_centers_cam(const Cube& ODF, const VectorJ& j, VoxelGrid* vgrid, Camera* cam, rave::EnvironmentBasePtr env) {
	Matrix4d cam_pose = cam->get_pose(j);
	Matrix<double,H_SUB,W_SUB> zbuffer = vgrid->get_zbuffer(cam_pose);
//	std::cout << zbuffer << "\n";

//	std::cout << "get all voxel centers in the FOV (in the camera frame)\n";
	StdVector3d voxel_centers_cam;
	Vector3i size = ODF.size();
	for(int i=0; i < size(0); ++i) {
		for(int j=0; j < size(1); ++j) {
			for(int k=0; k < size(2); ++k) {
				Vector3d voxel_center_world = vgrid->exact_point_from_voxel(Vector3d(i,j,k));
				if (cam->is_in_fov(voxel_center_world, zbuffer, cam_pose)) {
					rave_utils::plot_point(env, voxel_center_world, Vector3d(0,1,0), .005);
					Matrix4d voxel_center_pose_world = Matrix4d::Identity();
					voxel_center_pose_world.block<3,1>(0,3) = voxel_center_world;
					Matrix4d voxel_center_pose_cam = cam_pose.inverse()*voxel_center_pose_world;
					voxel_centers_cam.push_back(voxel_center_pose_cam.block<3,1>(0,3));
				}
			}
		}
	}
//	std::cout << "total number of voxels: " << size.prod() << "\n";
//	std::cout << "num voxels in FOV: " << voxel_centers_cam.size() << "\n";

	return voxel_centers_cam;
}

double average_object_distance(const VectorJ& j, const StdVector3d& voxel_centers_cam, const Cube& ODF, VoxelGrid* vgrid, Camera* cam) {
	Matrix4d cam_pose = cam->get_pose(j);

//	std::cout << "convert all voxel centers in the FOV to world frame, then get voxels\n";
	StdVector3d voxels_world;
	for(int i=0; i < voxel_centers_cam.size(); ++i) {
		Matrix4d voxel_center_pose_cam = Matrix4d::Identity();
		voxel_center_pose_cam.block<3,1>(0,3) = voxel_centers_cam[i];
		Vector3d voxel_center_pos_world = (cam_pose*voxel_center_pose_cam).block<3,1>(0,3);
		if (vgrid->contains(voxel_center_pos_world)) {
			Vector3d voxel = vgrid->exact_voxel_from_point(voxel_center_pos_world);
			voxels_world.push_back(voxel);
		}
	}

//	std::cout << "calculate average object-distance for voxels in FOV\n";
	double total_od = 0;
	for(int i=0; i < voxels_world.size(); ++i) {
		total_od += ODF.trilinear_interpolation(voxels_world[i]);
	}
	double avg_od = total_od / double(voxels_world.size());

	return avg_od;
}

double average_object_distance_cost(const VectorJ& j, const StdVector3d& voxel_centers_cam, const Cube& ODF, VoxelGrid* vgrid, Camera* cam) {
	const double alpha_obj_dist = 1;
	const double alpha_collision = 1e3;

	double cost = 0;
	cost += alpha_obj_dist*average_object_distance(j, voxel_centers_cam, ODF, vgrid, cam);
	cost += alpha_collision*(1.0/(1.0+exp(1e3*(vgrid->distance_to_TSDF(cam->get_position(j)) - intrinsics::MIN_RANGE))));
	return cost;
}

VectorJ average_object_distance_grad(const VectorJ& j, const StdVector3d& voxel_centers_cam, const Cube& ODF, VoxelGrid* vgrid, Camera* cam) {
	const double step = 1e-5;
	VectorJ grad;

	VectorJ j_p, j_m;
	for(int i=0; i < J_DIM; ++i) {
		j_p = j;
		j_m = j;

		j_p(i) = j(i) + step;
		j_m(i) = j(i) - step;

		grad(i) = (average_object_distance_cost(j_p, voxel_centers_cam, ODF, vgrid, cam) - average_object_distance_cost(j_m, voxel_centers_cam, ODF, vgrid, cam)) / (2*step);
	}

	return grad;
}

void test_avg_dist() {
//	srand(time(0));

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

	while(true) {
//		arm->teleop();
		std::cout << "Press enter\n";
		std::cin.ignore();
		arm->set_joint_values(j);
		j = arm->get_joint_values();

//		rave_utils::clear_plots();
//		rave_utils::plot_point(env, object, Vector3d(0,0,1), .05);
//		vgrid->plot_TSDF(env);
//		vgrid->plot_kinfu_tsdf(env);
		StdVector3d voxel_centers_cam = get_voxel_centers_cam(ODF, j, vgrid, cam, env);

		VectorJ grad = average_object_distance_grad(j, voxel_centers_cam, ODF, vgrid, cam);
		std::cout << "avg_od grad: " << grad.transpose() << "\n";

		rave_utils::clear_plots();
		rave_utils::plot_point(env, object, Vector3d(0,0,1), .05);
		vgrid->plot_TSDF(env);
		vgrid->plot_FOV(env, cam, vgrid->get_zbuffer(cam->get_pose(j)), cam->get_pose(j));
		double avg_od = average_object_distance(j, voxel_centers_cam, ODF, vgrid, cam);
		std::cout << "avg_od: " << avg_od << "\n";

		j -= (M_PI/32)*(grad/grad.norm());
	}

	std::cout << "Press enter to exit\n";
	std::cin.ignore();
}

int main() {
	test_avg_dist();
	return 0;
}
