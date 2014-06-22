#include "../pr2/pr2-sim.h"
#include "../geometry/geometry3d.h"
#include "../utils/rave-utils.h"
#include "../utils/mm-utils.h"
#include "../../util/Timer.h"

#include <openrave-core.h>
namespace rave = OpenRAVE;

#define M 1000

TimerCollection tc;

void test_fov(bool view=true) {
	PR2 brett(view);
	rave::EnvironmentBasePtr env = brett.get_env();
	sleep(2);

	// choose any camera
	brett.rarm->set_posture(Arm::Posture::mantis);
	Camera* cam = brett.rcam;

	// setup particles
	rave::KinBodyPtr table = env->GetKinBody("table");
	rave::KinBody::LinkPtr base = table->GetLink("base");
	rave::Vector extents = base->GetGeometry(0)->GetBoxExtents();

	rave::Vector table_pos = table->GetTransform().trans;
	double x_min, x_max, y_min, y_max, z_min, z_max;
	x_min = table_pos.x - extents.x;
	x_max = table_pos.x + extents.x;
	y_min = table_pos.y - extents.y;
	y_max = table_pos.y + extents.y;
	z_min = table_pos.z + extents.z - .1;
	z_max = table_pos.z + extents.z + .1;

	Matrix<double, 3, M> P;
	for(int m=0; m < M; ++m) {
		P(0,m) = mm_utils::uniform(x_min, x_max);
		P(1,m) = mm_utils::uniform(y_min, y_max);
		P(2,m) = mm_utils::uniform(z_min, z_max);
//		rave_utils::plot_point(env, P.col(m), {0,1,0});
	}

	tc.start("beams");
	std::vector<std::vector<Beam3d> > beams = cam->get_beams();
	tc.stop("beams");

	tc.start("border");
	std::vector<Triangle3d> border = cam->get_border(beams);
	tc.stop("border");

//	for(int i=0; i < beams.size(); ++i) {
//		for(int j=0; j < beams[i].size(); ++j) {
//			beams[i][j].plot(env);
//		}
//	}

	for(int i=0; i < border.size(); ++i) {
		border[i].plot(env);
	}

	tc.start("sd");
	Matrix<double,M,1> sd;
	for(int m=0; m < M; ++m) {
		sd(m) = cam->signed_distance(P.col(m), beams, border);
//		std::cout << "sd: " << sd(m) << "\n";
//		rave_utils::plot_point(env, P.col(m), {0,1,0});
//		std::cin.ignore();
//		rave_utils::clear_plots(3);
	}
	tc.stop("sd");

	tc.print_all_elapsed();

	std::cin.ignore();
}

int main() {
	test_fov();
	return 0;
}
