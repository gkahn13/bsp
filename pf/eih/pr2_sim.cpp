#include "pr2_sim.h"

/**
 * Constructors
 */

PR2::PR2() {
	std::string working_dir = boost::filesystem::current_path().normalize().string();
	std::string bsp_dir = working_dir.substr(0,working_dir.find("bsp"));
	std::string env_file = bsp_dir + "bsp/pf/eih/envs/pr2-test.env.xml";
	this->init(env_file, true);
}

PR2::PR2(std::string env_file, bool view) {
	this->init(env_file, view);
}

/**
 * Initializer methods
 */

void SetViewer(rave::EnvironmentBasePtr penv, const std::string& viewername)
{
    rave::ViewerBasePtr viewer = RaveCreateViewer(penv,viewername);
    BOOST_ASSERT(!!viewer);

    // attach it to the environment:
    penv->Add(viewer);

    // finally call the viewer's infinite loop (this is why a separate thread is needed)
    bool showgui = true;
    viewer->main(showgui);
}

void PR2::init(std::string env_file, bool view) {
	std::cout << "Initializing OpenRAVE\n";
	rave::RaveInitialize(true, rave::Level_Info);
	env = rave::RaveCreateEnvironment();
	std::cout << "Loading environment: " << env_file << "\n";
	env->Load(env_file);

	if (view) {
		boost::thread viewer_thread(boost::bind(SetViewer, env, "qtcoin"));
		viewer_thread.join(); // TODO: temp
	}
	env->Destroy(); // destroy
}


/**
 * For testing only
 */

int main(int argc, char* argv[]) {
	PR2 *brett = new PR2();
	return 0;
}
