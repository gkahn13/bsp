#include "../system/pr2-sim.h"
#include "../system/pr2-system.h"
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

void test_particle_update() {
	Vector3d object(3.35, -1.11, 0.8);
	Arm::ArmType arm_type = Arm::ArmType::right;
	bool view = true;
	PR2System sys(object, arm_type, view);

	PR2* brett = sys.get_brett();
	Arm* arm = sys.get_arm();
	rave::EnvironmentBasePtr env = brett->get_env();
	sleep(2);

	arm->set_posture(Arm::Posture::mantis);

	// setup scenario
	VectorJ j_t_real = arm->get_joint_values(), j_tp1_real;
	VectorJ j_t = j_t_real, j_tp1; // i.e. no belief == actual
	VectorU u_t = VectorU::Zero();
	MatrixP P_t = setup_particles(env), P_tp1;

	LOG_INFO("Origin particles");
	std::vector<ParticleGaussian> particle_gmm_t;
	sys.fit_gaussians_to_pf(P_t, particle_gmm_t);
	sys.display(j_t, particle_gmm_t);

	sys.execute_control_step(j_t_real, j_t, u_t, P_t, j_tp1_real, j_tp1, P_tp1);

	LOG_INFO("New particles")
	std::vector<ParticleGaussian> particle_gmm_tp1;
	sys.fit_gaussians_to_pf(P_tp1, particle_gmm_tp1);
	sys.display(j_tp1, particle_gmm_tp1);
}

void test_figtree() {
	Vector3d object(3.35, -1.11, 0.8);
	Arm::ArmType arm_type = Arm::ArmType::right;
	bool view = true;
	PR2System sys(object, arm_type, view);

	PR2* brett = sys.get_brett();
	Arm* arm = sys.get_arm();
	rave::EnvironmentBasePtr env = brett->get_env();
	arm->set_posture(Arm::Posture::mantis);
	sleep(1);

	MatrixP P = setup_particles(env);
	std::vector<ParticleGaussian> particle_gmm;

	tc.start("figtree");
	sys.fit_gaussians_to_pf(P, particle_gmm);
	tc.stop("figtree");

	std::cout << "particle_gmm size: " << particle_gmm.size() << "\n";
	for(int i=0; i < particle_gmm.size(); ++i) {
		std::cout << "mean: " << particle_gmm[i].mean.transpose() << "\n";
		std::cout << "cov diag: " << particle_gmm[i].cov.diagonal().transpose() << "\n";
		std::cout << "pct: " << particle_gmm[i].pct << "\n";
	}

	tc.print_all_elapsed();
	sys.display(arm->get_joint_values(), particle_gmm);
}

void test_pr2_system() {
	Vector3d object(3.35, -1.11, 0.8);
	Arm::ArmType arm_type = Arm::ArmType::right;
	bool view = true;
	PR2System sys(object, arm_type, view);

	PR2* brett = sys.get_brett();
	Arm* arm = sys.get_arm();
	rave::EnvironmentBasePtr env = brett->get_env();
	sleep(2);

	arm->set_posture(Arm::Posture::mantis);

	// setup particles
	MatrixP P = setup_particles(env);

	// test plotting
	VectorJ j = arm->get_joint_values();

	std::vector<ParticleGaussian> particle_gmm;
	particle_gmm.push_back(ParticleGaussian(Vector3d::Zero(), Matrix3d::Identity(), P.leftCols(M_DIM/2), 1));
	particle_gmm.push_back(ParticleGaussian(Vector3d::Zero(), Matrix3d::Identity(), P.rightCols(M_DIM/2), 1));

	sys.display(j, particle_gmm);

}

int main() {
	test_particle_update();
//	test_figtree();
//	test_pr2_system();
	return 0;
}
