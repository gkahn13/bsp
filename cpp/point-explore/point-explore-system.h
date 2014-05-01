#ifndef __POINT_EXPLORE_SYSTEM_H__
#define __POINT_EXPLORE_SYSTEM_H__

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

namespace Constants {
	const double alpha = 2.5;
	const double max_range = 0.25;

	const double alpha_control_norm = 0; // 1e-2
	const double alpha_control_smooth = 0; // 1e-2
	const double alpha_separation = 0; // 1e-3
}


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

enum class ObsType { angle, distance};
enum class CostType { entropy, platt};

typedef std::vector<ObsType> ObsTypeList;
inline std::istream& operator>>(std::istream& in, ObsType& obs_type)
{
    std::string token;
    in >> token;
    if (token == "angle") {
        obs_type = ObsType::angle;
    } else if (token == "distance") {
        obs_type = ObsType::distance;
    } else {
    	throw po::validation_error(po::validation_error::invalid_option_value);
    }

    return in;
}

typedef std::vector<CostType> CostTypeList;
inline std::istream& operator>>(std::istream& in, CostType& cost_type)
{
    std::string token;
    in >> token;
    if (token == "entropy") {
        cost_type = CostType::entropy;
    } else if (token == "platt") {
        cost_type = CostType::platt;
    } else {
    	throw po::validation_error(po::validation_error::invalid_option_value);
    }

    return in;
}

#include "casadi/casadi-point-explore-system.h"
namespace AD = CasADi;

class PointExploreSystem {
public:
	PointExploreSystem();
	PointExploreSystem(mat& target, const ObsType obs_type, const CostType cost_type, bool use_casadi);
	PointExploreSystem(mat& target, const ObsType obs_type, const CostType cost_type, bool use_casadi,
			int T, int M, int N, double DT, int X_DIM, int U_DIM, int Z_DIM, int Q_DIM, int R_DIM);

	mat dynfunc(const mat& x, const mat& u);
	mat obsfunc(const mat& x, const mat& t, const mat& r);

	void update_state_and_particles(const mat& x_t, const mat& P_t, const mat& u_t, mat& x_tp1, mat& P_tp1);
	double cost(const std::vector<mat>& X, const std::vector<mat>& U, const mat& P);
	mat cost_grad(std::vector<mat>& X, std::vector<mat>& U, const mat& P);

	void display_states_and_particles(const std::vector<mat>& X, const mat& P);

	mat get_target() { return this->target; }
	mat get_xMin() { return this->xMin; }
	mat get_xMax() { return this->xMax; }
	mat get_uMin() { return this->uMin; }
	mat get_uMax() { return this->uMax; }

protected:
	int T, M, N, TOTAL_VARS, X_DIM, U_DIM, Z_DIM, Q_DIM, R_DIM;
	double DT;

	mat target;
	ObsType obs_type;
	CostType cost_type;
	bool use_casadi;

	mat R;
	mat xMin, xMax, uMin, uMax;

	CasadiPointExploreSystem* casadi_sys;

	void init_dims(int T=10, int M=100, int N=1, double DT=1.0, int X_DIM=2, int U_DIM=2, int Z_DIM=1, int Q_DIM=2, int R_DIM=1);
	void init(mat& target, const ObsType obs_type, const CostType cost_type, bool use_casadi,
			mat& xMin, mat& xMax, mat& uMin, mat& uMax, mat& R);

	mat obsfunc_dist(const mat& x, const mat& t, const mat& r);
	mat obsfunc_angle(const mat& x, const mat& t, const mat& r);

	double gauss_likelihood(const mat& v, const mat& S);
	mat low_variance_sampler(const mat& P, const mat& W, double r);

	double cost_entropy(const std::vector<mat>& X, const std::vector<mat>& U, const mat& P);
	double cost_platt(const std::vector<mat>& X, const std::vector<mat>& U, const mat& P);
};

#endif
