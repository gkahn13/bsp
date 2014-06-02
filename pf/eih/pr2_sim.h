#pragma once

#include <armadillo>
using namespace arma;

#include <boost/thread/thread.hpp>
#include <boost/bind.hpp>
#include <boost/filesystem.hpp>

#include <openrave-core.h>
namespace rave = OpenRAVE;

class PR2 {
public:
	PR2();
	PR2(std::string env_file, bool view=true);

private:
	void init(std::string env_file, bool view);

	rave::EnvironmentBasePtr env;
	rave::ViewerBasePtr viewer;
	boost::thread viewer_thread;
};
