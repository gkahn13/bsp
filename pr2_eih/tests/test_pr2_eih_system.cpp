#include "system/pr2_eih_system.h"

double cost(const VectorU& u, const VectorJ& j, const MatrixJ& j_sigma0,
		const Vector3d& obj, const Matrix3d& obj_sigma0, const double alpha,
		const std::vector<geometry3d::Triangle>& obstacles, PR2EihSystem& sys) {
	double cost = 0;

	VectorX x_t, x_tp1 = VectorX::Zero();
	MatrixX sigma_t = MatrixX::Zero(), sigma_tp1 = MatrixX::Zero();
	sigma_t.block<J_DIM,J_DIM>(0,0) = j_sigma0;
	sigma_t.block<3,3>(J_DIM,J_DIM) = obj_sigma0;

	x_t << j, obj;
	sys.belief_dynamics(x_t, sigma_t, u, alpha, obstacles, x_tp1, sigma_tp1);

	cost += sigma_tp1.trace();

	return cost;
}

VectorU grad_cost(const VectorU& u, const VectorJ& j, const MatrixJ& j_sigma0,
		const Vector3d& obj, const Matrix3d& obj_sigma0, const double alpha,
		const std::vector<geometry3d::Triangle>& obstacles, PR2EihSystem& sys) {
	const double step = 1e-5;
	VectorU grad;

	VectorU u_p = u, u_m = u;
	double cost_p, cost_m;
	for(int i=0; i < U_DIM; ++i) {
		u_p(i) = u(i) + step;
		u_m(i) = u(i) - step;

		cost_p = cost(u_p, j, j_sigma0, obj, obj_sigma0, alpha, obstacles, sys);
		cost_m = cost(u_m, j, j_sigma0, obj, obj_sigma0, alpha, obstacles, sys);
		grad(i) = (cost_p - cost_m) / (2*step);

		u_p(i) = u(i);
		u_m(i) = u(i);
	}

	return grad;
}

void test_cost() {
	// setup system
	pr2_sim::Simulator sim(true, false);
	pr2_sim::Arm arm(pr2_sim::Arm::right, &sim);
	pr2_sim::Camera cam(&arm, &sim);
	PR2EihSystem sys(&sim, &arm, &cam);

	// setup environment
	Vector3d table_center(.2,.7,.5);
	std::vector<geometry3d::Triangle> obstacles = {
//			geometry3d::Triangle(table_center, table_center+Vector3d(.5,-1.4,0), table_center+Vector3d(.5,0,0)),
//			geometry3d::Triangle(table_center, table_center+Vector3d(0,-1.4,0), table_center+Vector3d(.5,-1.4,0)),
			geometry3d::Triangle(table_center+Vector3d(.25,-.7,0), table_center+Vector3d(.25,-.7,.2), table_center+Vector3d(.25,-1.2,0)),
			geometry3d::Triangle(table_center+Vector3d(.25,-1.2,0), table_center+Vector3d(.25,-.7,.2), table_center+Vector3d(.25,-1.2,.2))};

	// setup initial state
	arm.set_posture(pr2_sim::Arm::Posture::mantis);
	VectorJ j_t = arm.get_joints();
	MatrixJ j_sigma0 = .2*MatrixJ::Identity();

	Vector3d obj(0.640, -0.170, 0.560);
//	Vector3d obj(0.500, -0.170, 0.560);
	Matrix3d obj_sigma = Vector3d(.1,.1,.1).asDiagonal();

	VectorU u_t = VectorU::Zero();

	const double alpha = 1;
	double c = cost(u_t, j_t, j_sigma0, obj, obj_sigma, alpha, obstacles, sys);
	std::cout << "cost: " << c << "\n";

	VectorU grad = grad_cost(u_t, j_t, j_sigma0, obj, obj_sigma, alpha, obstacles, sys);
	std::cout << "grad: " << grad.transpose() << "\n\n";

	sys.plot(StdVectorJ(1, j_t), obj, obj_sigma, obstacles);

	VectorJ j_tp1 = j_t - (M_PI/16)*(grad / grad.norm());
	arm.set_joints(j_tp1);

	sys.plot(StdVectorJ(1, j_tp1), obj, obj_sigma, obstacles);

}

int main() {
	test_cost();
	return 0;
}
