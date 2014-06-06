#include "../include/eih_system.h"
#include "../include/pr2_sim.h"
#include "../include/rave_utils.h"
#include "../include/utils.h"
#include "../../../util/Timer.h"

/**
 * TESTS
 */

void test_arm() {
	PR2 *brett = new PR2();
	boost::this_thread::sleep(boost::posix_time::seconds(1));

	brett->larm->set_posture(Arm::Posture::side);
	rave::Transform start_pose = brett->larm->get_pose();
	std::cout << "start joints:\n" << brett->larm->get_joint_values() << "\n";
	std::cout << "start pose:\n" << rave_utils::rave_transform_to_mat(start_pose) << "\n";

	brett->larm->set_posture(Arm::Posture::mantis);
	std::cout << "mantis joints:\n" << brett->larm->get_joint_values() << "\n";
	std::cout << "In mantis. Press enter to go back to start pose\n";
	std::cin.ignore();

	brett->larm->set_pose(start_pose);
	std::cout << "end joints:\n" << brett->larm->get_joint_values() << "\n";

	std::cout << "Press enter to quit\n";
	std::cin.ignore();
}

void test_teleop() {
	PR2 *brett = new PR2();
	boost::this_thread::sleep(boost::posix_time::seconds(1));

	brett->larm->set_posture(Arm::Posture::mantis);
	brett->rarm->set_posture(Arm::Posture::mantis);

	brett->r_kinect->power_on();
	brett->r_kinect->render_on();
	while(true) {
		brett->rarm->teleop();

		std::cout << "\n" << brett->rarm->get_joint_values().t();
		std::cout << "Press enter to continue\n";
		std::cin.ignore();
	}
}

void test_head() {
	PR2 *brett = new PR2();
	boost::this_thread::sleep(boost::posix_time::seconds(1));

	Arm *arm = brett->rarm;
	Head *head = brett->head;

	arm->set_posture(Arm::Posture::mantis);

	mat j = head->get_joint_values();
	std::cout << "head joints: " << j.t() << "\n";

	mat lower, upper;
	head->get_limits(lower, upper);
	std::cout << "lower: " << lower.t() << "\n";
	std::cout << "upper: " << upper.t() << "\n";

	arm->teleop();
	head->set_pose(arm->get_pose());
}

void test_camera() {
	PR2 *brett = new PR2();
	KinectSensor *kinect = brett->h_kinect;

	kinect->power_on();
	boost::this_thread::sleep(boost::posix_time::seconds(1));
	std::cout << "Getting image...\n";
	cube image = kinect->get_image();
	for(int i=0; i < image.n_rows; ++i) {
		for(int j=0; j < image.n_cols; ++j) {
			mat rgb = image.subcube(i,j,0,i,j,2);
//			std::cout << rgb.n_rows << " " << rgb.n_cols << "\n";
//			std::cout << rgb;
		}
	}
	std::cout << "Image size: (" << image.n_rows << ", " << image.n_cols << ", " << image.n_slices << ")\n";

//	std::vector<std::vector<mat> > points_pixels_colors = camera->get_pixels_and_colors()

}

void test_plot() {
	PR2 *brett = new PR2();
	Arm *arm = brett->rarm;
	KinectSensor *kinect = brett->r_kinect;

	arm->set_posture(Arm::Posture::mantis);
	boost::this_thread::sleep(boost::posix_time::seconds(1));

	rave::Transform p = kinect->get_pose();

	std::vector<rave::GraphHandlePtr> h;
	rave_utils::plot_transform(brett->get_env(), p, h);

//	rave::Vector pos = p.trans;
//	rave::Vector color(0, 1, 0);
//	rave::GraphHandlePtr h = rave_utils::plot_point(brett->get_env(), pos, color);

	std::cout << "Press enter to exit\n";
	std::cin.ignore();
}

void test_kinect() {
	PR2 *brett = new PR2();
	Arm *arm = brett->rarm;
	KinectSensor *kinect = brett->r_kinect;

	kinect->power_on();
	kinect->render_on();
	arm->set_posture(Arm::Posture::mantis);
	boost::this_thread::sleep(boost::posix_time::seconds(4));

	std::vector<rave::GraphHandlePtr> handles;
	while(true) {
		arm->teleop();
		handles.clear();

		std::vector<ColoredPoint*> colored_points = kinect->get_point_cloud();
		kinect->display_point_cloud(colored_points, handles);
	}

	std::cout << "Press enter to exit\n";
	std::cin.ignore();
}

void setup_eih_environment(PR2 *brett, int M, Arm::ArmType arm_type, bool zero_seed, mat &P,
		EihSystem **sys, EihSystem::ObsType obs_type) {
	if (zero_seed) {
		srand(time(0));
	}

	rave::EnvironmentBasePtr env = brett->get_env();

	Arm *larm = brett->larm;
	Arm *rarm = brett->rarm;
	KinectSensor *l_kinect = brett->l_kinect;
	KinectSensor *r_kinect = brett->r_kinect;

	larm->set_posture(Arm::Posture::mantis);
	rarm->set_posture(Arm::Posture::mantis);

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

	for(int m=0; m < M; ++m) {
		P(0,m) = uniform(x_min, x_max);
		P(1,m) = uniform(y_min, y_max);
		P(2,m) = uniform(z_min, z_max);
	}

	Manipulator *manip;
	KinectSensor *kinect;
	if (arm_type == Arm::ArmType::left) {
		manip = larm;
		kinect = l_kinect;
	} else {
		manip = rarm;
		kinect = r_kinect;
	}
	*sys = new EihSystem(env, manip, kinect, obs_type);
	kinect->render_on();
	boost::this_thread::sleep(boost::posix_time::seconds(2));
}

void test_pf_update() {
	PR2 *brett = new PR2();
	int M = 1000;
	mat P(3,1000,fill::zeros);
	EihSystem *sys = NULL;
	Manipulator *manip = NULL;
	KinectSensor *kinect = NULL;
	setup_eih_environment(brett, 1000, Arm::ArmType::right, true, P, &sys, EihSystem::ObsType::fov);
	manip = sys->get_manip();
	kinect = sys->get_kinect();

	rave::EnvironmentBasePtr env = brett->get_env();

//	std::vector<rave::GraphHandlePtr> handles;
//	mat color;
//	color << 1 << 0 << 0;
//	for(int m=0; m < M; ++m) {
//		handles.push_back(rave_utils::plot_point(env, P.col(m), color));
//	}

	mat x_t = manip->get_joint_values();
	mat P_t = P;
	mat u_t(x_t.n_rows, x_t.n_cols, fill::zeros);
	mat x_tp1(x_t.n_rows, x_t.n_cols, fill::zeros), P_tp1(3, M, fill::zeros);

	std::cout << "x_t: " << x_t << "\n";
	std::cout << "u_t: " << u_t << "\n";

	while(true) {
		manip->teleop();
		x_t = manip->get_joint_values();
		sys->update_state_and_particles(x_t, P_t, u_t, x_tp1, P_tp1, true);

		P_t = P_tp1;
	}

	std::cout << "Press enter to exit\n";
	std::cin.ignore();
}

void test_cost() {
	PR2 *brett = new PR2();
	int M = 1000;
	mat P(3,1000,fill::zeros);
	EihSystem *sys = NULL;
	Manipulator *manip = NULL;
	KinectSensor *kinect = NULL;
	setup_eih_environment(brett, 1000, Arm::ArmType::right, false, P, &sys, EihSystem::ObsType::fov_occluded_color);
	manip = sys->get_manip();
	kinect = sys->get_kinect();

	rave::EnvironmentBasePtr env = brett->get_env();

	int T = 2;
	std::vector<mat> X_DES(T);

	for(int t=0; t < T; ++t) {
		manip->teleop();
		X_DES[t] = manip->get_joint_values().t();
	}

//	X_DES[0] << -2.0302  << -0.0547 << -1.0110 << -1.4762 << -0.5600 << -1.4286 << -3.9647;
//	X_DES[1] << -1.9594 << -0.0394 << -1.0110 << -1.3647 << -0.4439 << -1.6257 << -3.8557;
//	X_DES[2] << -1.9594 << -0.0394 << -1.0110 << -1.3647 << -0.4439 << -1.6257 << -3.8557;

	mat x0 = X_DES[0];

	std::vector<mat> U(T-1);
	for(int t=0; t < T-1; ++t) {
		U[t] = X_DES[t+1] - X_DES[t];
	}

	double cost = sys->cost(x0, U, P);
	std::cout << "cost: " << cost << "\n";

	mat cost_grad = sys->cost_grad(x0, U, P);
	std::cout << "cost grad:\n" << cost_grad;
}

void test_greedy() {
	PR2 *brett = new PR2();
	int M = 1000;
	mat P(3,1000,fill::zeros);
	EihSystem *sys = NULL;
	Manipulator *manip = NULL;
	KinectSensor *kinect = NULL;
	setup_eih_environment(brett, 1000, Arm::ArmType::right, false, P, &sys, EihSystem::ObsType::fov);
	manip = sys->get_manip();
	kinect = sys->get_kinect();

	rave::EnvironmentBasePtr env = brett->get_env();

	int T = 2, ch;
	mat x0 = manip->get_joint_values();
	std::vector<mat> U(1);
	while(true) {
//		std::cout << "\nPress 'q' to exit\n";
//		if (utils::getch() == 'q') { break; }
//		std::cout << "\n";

//		manip->teleop();

		U[0] = zeros<mat>(x0.n_rows, x0.n_cols);

		double cost = sys->cost(x0, U, P);
		std::cout << "cost: " << cost << "\n";

		mat cost_grad = sys->cost_grad(x0, U, P);
		std::cout << "cost grad: " << cost_grad.t();

		// step along negative gradient
//		x0 -= (2*(M_PI/180))*(cost_grad/norm(cost_grad,2));
//		manip->set_joint_values(x0);

		// step along negative gradient
		// (clip joint input)
//		mat u = -(2*(M_PI/180))*(cost_grad/norm(cost_grad,2));
//		mat x_tp1(x0.n_rows, x0.n_cols), P_tp1(3, M, fill::zeros);
//		sys->update_state_and_particles(x0, P, u, x_tp1, P_tp1, true);
//		x0 = x_tp1;
//		P = P_tp1;

		// step along negative gradient
		// (clip pos/rot change)
		mat u_grad = -(2*(M_PI/180))*(cost_grad/norm(cost_grad,2));
		mat x1 = x0 + u_grad;

		rave::Transform p0, p1;
		p0 = manip->get_pose();
		manip->set_joint_values(x1);
		p1 = manip->get_pose();

		rave::Vector trans_diff_dir = (p1.trans - p0.trans).normalize3();
		rave::Vector trans1_clipped = p0.trans + .03*(trans_diff_dir);
		rave::Vector rot1_clipped = rave::geometry::quatSlerp(p0.rot, p1.rot, 1.0);
		rave::Transform p1_clipped(rot1_clipped, trans1_clipped);

		manip->set_pose(p1_clipped);
		mat x1_clipped = manip->get_joint_values();
		mat u = x1_clipped - x0;

		manip->set_joint_values(x0);
		mat x_tp1(x0.n_rows, x0.n_cols), P_tp1(3, M, fill::zeros);
		sys->update_state_and_particles(x0, P, u, x_tp1, P_tp1, true);
		x0 = x_tp1;
		P = P_tp1;
	}
}

void test_stamps() {
	PR2 *brett = new PR2();

	Arm *arm = brett->rarm;
	KinectSensor *kinect = brett->r_kinect;
	DepthSensor *depth = kinect->get_depth_sensor();
	CameraSensor *cam = kinect->get_camera_sensor();

	arm->set_posture(Arm::Posture::mantis);
	kinect->power_on();
	kinect->render_on();
	boost::this_thread::sleep(boost::posix_time::seconds(2));


	double data_time_elapse;
	util::Timer data_timer;
	while(true) {
		util::Timer_tic(&data_timer);
		kinect->get_point_cloud(true);
		kinect->get_z_buffer(false);
		data_time_elapse = util::Timer_toc(&data_timer);
		std::cout << data_time_elapse << " (s)\n";
	}


}


int main(int argc, char* argv[]) {
//	test_arm();
//	test_teleop();
//	test_head();
//	test_camera();
//	test_plot();
//	test_kinect();
//	test_pf_update();
//	test_cost();
	test_greedy();
//	test_stamps();
}

