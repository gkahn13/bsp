#ifndef __MM_UTILS_H__
#define __MM_UTILS_H__

#include <map>
#include <iostream>

#include <Eigen/Eigen>
using namespace Eigen;

#include "../../util/Timer.h"
#include "../../util/logging.h"

/**
 * Wrapper that tracks multiple timers
 */
class TimerCollection {
public:
	void start(const std::string& timer_name);
	void stop(const std::string& timer_name);

	void print_elapsed(const std::string& timer_name);
	void print_all_elapsed();

	void clear_all();
private:
	std::map<std::string,util::Timer> timer_dict;
	std::map<std::string,double> timer_elapsed;
	std::map<std::string,bool> timer_running;
};


namespace pr2_utils {

inline double uniform(double low, double high) {
	return (high - low)*(rand() / double(RAND_MAX)) + low;
}

}

#endif
