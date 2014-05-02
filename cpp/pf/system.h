#ifndef __SYSTEM_H__
#define __SYSTEM_H__

#include <Python.h>

#include <vector>
#include <cmath>

#include "../util/logging.h"

#include <boost/python.hpp>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
namespace py = boost::python;
namespace po = boost::program_options;


#include <symbolic/casadi.hpp>
#include <symbolic/stl_vector_tools.hpp>
#include <cstdlib>

#include <armadillo>
using namespace arma;

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

const double step = 0.0078125*0.0078125;
const double INFTY = 1e10;


inline double uniform(double low, double high) {
	return (high - low)*(rand() / double(RAND_MAX)) + low;
}

inline mat sample_gaussian(mat mean, mat covariance) {
	mat sample(mean.n_rows, mean.n_cols, fill::randn);

	mat L = chol(covariance);
	sample = L*sample + mean;

	return sample;
}

inline mat subsample(mat& P, int size) {
	mat P_shuffled = shuffle(P, 1);
	mat P_subsampled = P_shuffled.cols(0, size-1);
	return P_subsampled;
}

inline std::istream& operator>>(std::istream& in, mat& m)
{
    std::string token;
    in >> token;
    std::cout << token << "\n";
    try {
    	std::istringstream iss(token);
    	std::string word;
    	int index = 0;
    	while(iss >> word) {
    		m(index++) = std::atof(word.c_str());
    	}
    } catch (std::exception &e) {
    	throw po::validation_error(po::validation_error::invalid_option_value);
    }

    return in;
}

#include "casadi-system.h"
namespace AD = CasADi;

class System {
public:
	System();

	virtual mat dynfunc(const mat& x, const mat& u) =0;
	virtual mat obsfunc(const mat& x, const mat& t, const mat& r) =0;

	virtual void update_state_and_particles(const mat& x_t, const mat& P_t, const mat& u_t, mat& x_tp1, mat& P_tp1) =0;
	virtual double cost(const std::vector<mat>& X, const std::vector<mat>& U, const mat& P) =0;
	mat cost_grad(std::vector<mat>& X, std::vector<mat>& U, const mat& P);

	virtual void display_states_and_particles(const std::vector<mat>& X, const mat& P, bool pause=true) =0;

	mat get_xMin() { return this->xMin; }
	mat get_xMax() { return this->xMax; }
	mat get_uMin() { return this->uMin; }
	mat get_uMax() { return this->uMax; }

protected:
	int T, M, TOTAL_VARS, X_DIM, U_DIM, Z_DIM, Q_DIM, R_DIM;
	double DT;

	bool use_casadi;

	mat R;
	mat xMin, xMax, uMin, uMax;

	CasadiSystem* casadi_sys;

	double gauss_likelihood(const mat& v, const mat& S);
	mat low_variance_sampler(const mat& P, const mat& W, double r);

	double cost_entropy(const std::vector<mat>& X, const std::vector<mat>& U, const mat& P);
	double cost_platt(const std::vector<mat>& X, const std::vector<mat>& U, const mat& P);
};

#endif
