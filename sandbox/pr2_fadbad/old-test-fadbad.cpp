#include "../system/pr2-sim.h"
#include "../system/pr2-system.h"
#include "../geometry/geometry3d.h"
#include "../utils/rave-utils.h"
#include "../utils/pr2-utils.h"
#include "../../util/Timer.h"

#include "../fadbad/geometry/fadbad-geometry3d.h"
#include "../fadbad/fadbad-utils.h"
#include "../fadbad/system/fadbad-pr2-sim.h"
#include "../fadbad/system/fadbad-pr2-system.h"

#include <openrave-core.h>
namespace rave = OpenRAVE;

const int T = TIMESTEPS;

TimerCollection tc;

void init_cost(const VectorJ& j0, const MatrixP& P, PR2System& sys,
		StdVectorJ& J, StdVectorU& U, std::vector<ParticleGaussian>& particle_gmm) {
	// re-initialize GMM from PF
	sys.fit_gaussians_to_pf(P, particle_gmm);

	for(int i=0; i < particle_gmm.size(); ++i) {
		particle_gmm[i].cov *= 5000;
		std::cout << particle_gmm[i].pct << "\n";
		std::cout << particle_gmm[i].cov << "\n\n";
	}

	// only take max gaussian
	particle_gmm = std::vector<ParticleGaussian>(1, particle_gmm[0]);
//	particle_gmm[0].cov = .01*Matrix3d::Identity();
//	particle_gmm[0].cov *= 5000;

	// TODO: try a non-zero initialization
	VectorU uinit = VectorU::Zero();

	// integrate trajectory
	J[0] = j0;
	for(int t=0; t < T-1; ++t) {
		U[t] = uinit;
		J[t+1] = sys.dynfunc(J[t], U[t], VectorQ::Zero());
	}
}

MatrixP init_particles(rave::EnvironmentBasePtr env) {
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

	return P;
}

void test_cost() {
	Vector3d object(3.35, -1.11, 0.8);
	Arm::ArmType arm_type = Arm::ArmType::right;
	bool view = true;
	bool use_fadbad = true;
	PR2System sys(object, arm_type, view, use_fadbad);

	PR2* brett = sys.get_brett();
	Arm* arm = sys.get_arm();
	rave::EnvironmentBasePtr env = brett->get_env();
	arm->set_posture(Arm::Posture::mantis);
	sleep(2);

	MatrixP P = init_particles(env);

	VectorJ j0 = arm->get_joint_values();
	sys.get_pcl(j0);
	MatrixJ j_sigma0 = .01*MatrixJ::Identity(); // TODO: never actually update it in MPC

	StdVectorJ J(T);
	StdVectorU U(T-1);
	std::vector<ParticleGaussian> particle_gmm;

	init_cost(j0, P, sys, J, U, particle_gmm);

//	std::cout << "Displaying initial pose with particles\n";
//	sys.display(j0, particle_gmm);

	double alpha = 1;

	tc.start("cost_grad");
	double cost = sys.cost_gmm(J, j_sigma0, U, particle_gmm, alpha);
//	VectorTOTAL grad = sys.cost_gmm_grad(J, j_sigma0, U, particle_gmm, alpha);
	tc.stop("cost_grad");

	double cost_f;
	VectorTOTAL grad_f;
	tc.start("cost_grad_fadbad");
	sys.cost_gmm_and_grad(J, j_sigma0, U, particle_gmm, alpha, cost_f, grad_f);
	tc.stop("cost_grad_fadbad");

	std::cout << "cost: " << cost << "\n";
	std::cout << "cost_f: " << cost_f << "\n\n";

//	std::cout << "grad_f:\n" << grad_f << "\n\n";
//	for(int i=0; i < TOTAL_VARS; ++i) {
//		std::cout << grad(i) << "\t" << grad_f(i) << "\n";
//	}
//	std::cout << "\n";

	tc.print_all_elapsed();

}

int main() {
	test_cost();
	return 0;
}
