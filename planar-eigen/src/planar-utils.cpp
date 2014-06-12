#include "../include/planar-utils.h"

/**
 * TimerCollection public methods
 */

void TimerCollection::start(const std::string& timer_name) {
	if (timer_dict.count(timer_name) == 0) {
		timer_dict[timer_name] = util::Timer();
		timer_elapsed[timer_name] = 0;
		timer_running[timer_name] = false;
	}

	if (!timer_running[timer_name]) {
		util::Timer_tic(&timer_dict[timer_name]);
		timer_running[timer_name] = true;
	}
}

void TimerCollection::stop(const std::string& timer_name) {
	if (timer_dict.count(timer_name) == 0) {
		LOG_WARN("Attempted to stop a timer that was never started. Ignoring");
		return;
	}

	timer_elapsed[timer_name] += util::Timer_toc(&timer_dict[timer_name]);
	timer_running[timer_name] = false;
}

void TimerCollection::print_elapsed(const std::string& timer_name) {
	if (timer_dict.count(timer_name) > 0) {
		std::cout << timer_name << " elapsed: " << timer_elapsed[timer_name] << " seconds\n";
	}
}

void TimerCollection::print_all_elapsed() {
	for(std::map<std::string,double>::iterator iter = timer_elapsed.begin(); iter != timer_elapsed.end(); ++iter) {
		std::cout << iter->first << " elapsed: " << iter->second << " seconds\n";
	}
}

void TimerCollection::clear_all() {
	timer_dict.clear();
	timer_elapsed.clear();
	timer_running.clear();
}

namespace planar_utils {

np::ndarray eigen_to_ndarray(VectorXd& m) {
	py::tuple shape = py::make_tuple(m.rows());
	np::dtype dtype = np::dtype::get_builtin<float>();
	np::ndarray n = np::zeros(shape, dtype);
	for(int i=0; i < m.rows(); ++i) {
		n[py::make_tuple(i)] = m(i);
	}

	return n;
}


np::ndarray eigen_to_ndarray(MatrixXd& m) {
	py::tuple shape = py::make_tuple(m.rows(), m.cols());
	np::dtype dtype = np::dtype::get_builtin<float>();
	np::ndarray n = np::zeros(shape, dtype);
	for(int i=0; i < m.rows(); ++i) {
		for(int j=0; j < m.cols(); ++j) {
			n[py::make_tuple(i,j)] = m(i,j);
		}
	}

	return n;
}

}
