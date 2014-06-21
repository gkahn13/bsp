#include "../include/camera.h"
#include "../include/geometry3d.h"
#include "../include/rave_utils.h"
#include "../../util/Timer.h"

#include <openrave-core.h>
namespace rave = OpenRAVE;

void SetViewer(rave::EnvironmentBasePtr penv, const std::string& viewername) {
    rave::ViewerBasePtr viewer = rave::RaveCreateViewer(penv,viewername);
    BOOST_ASSERT(!!viewer);

    // attach it to the environment:
    penv->Add(viewer);

    // finally call the viewer's infinite loop (this is why a separate thread is needed)
    bool showgui = true;
    viewer->main(showgui);
}

void test_fov() {
	std::cout << "Initializing OpenRAVE\n";
	rave::RaveInitialize(true, rave::Level_Info);
	rave::EnvironmentBasePtr env = rave::RaveCreateEnvironment();
	std::cout << "Loading environment";

	std::string env_file = "/home/gkahn/Research/bsp/multimodal3d/envs/pr2-test.env.xml";
	env->Load(env_file);

	rave::RobotBasePtr robot = env->GetRobot("Brett");

	boost::shared_ptr<boost::thread> viewer_thread(new boost::thread(boost::bind(SetViewer, env, "qtcoin")));
	sleep(2);

	rave::SensorBasePtr sensor;
	std::vector<rave::RobotBase::AttachedSensorPtr> sensors = robot->GetAttachedSensors();
	for(int i=0; i < sensors.size(); ++i) {
		if (sensors[i]->GetName() == "head_cam") {
			sensor = sensors[i]->GetSensor();
		}
	}

//	Vector3d color(1,0,0);
//	Vector3d start(0,0,0), end(1,0,0);
//	rave_utils::plot_segment(env, start, end, color);
//	rave::GraphHandlePtr h = rave_utils::plot_point(env, sensor->GetTransform().trans, rave::Vector(1,0,0), .1);

	double max_range = 5;
	Camera cam(robot, sensor, max_range);
//	cam.get_directions();

	util::Timer beams_timer;
	util::Timer_tic(&beams_timer);
	std::vector<std::vector<Beam3d> > beams = cam.get_beams();
	double beams_time = util::Timer_toc(&beams_timer);
	std::cout << "beams_time: " << beams_time << "\n";

	for(int i=0; i < beams.size(); ++i) {
		for(int j=0; j < beams[i].size(); ++j) {
			beams[i][j].plot(env);
		}
	}

	std::cin.ignore();
}

int main() {
	test_fov();
	return 0;
}
