#ifndef __BOXES_SYSTEM_H__
#define __BOXES_SYSTEM_H__

#include "../system.h"

#include <boost/program_options.hpp>
namespace po = boost::program_options;

namespace Constants {

}

enum class ObsType { distance};
enum class CostType { entropy, platt};

typedef std::vector<ObsType> ObsTypeList;
inline std::istream& operator>>(std::istream& in, ObsType& obs_type)
{
    std::string token;
    in >> token;
    if (token == "distance") {
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

#include "casadi-boxes-system.h"

class BoxesSystem : public virtual System {
public:
	BoxesSystem();
	BoxesSystem(mat& boxes, const ObsType obs_type, const CostType cost_type, bool use_casadi);
	BoxesSystem(mat& boxes, const ObsType obs_type, const CostType cost_type, bool use_casadi,
			int T, int M, int N, double DT, int B_DIM, int X_DIM, int U_DIM, int Z_DIM, int Q_DIM, int R_DIM);

	mat dynfunc(const mat& x, const mat& u);
	mat obsfunc(const mat& x, const mat& b, const mat& r);

	void update_state_and_particles(const mat& x_t, const mat& P_t, const mat& u_t, mat& x_tp1, mat& P_tp1);
	double cost(const std::vector<mat>& X, const std::vector<mat>& U, const mat& P);

	void display_states_and_particles(const std::vector<mat>& X, const mat& P);

	mat get_boxes() { return this->boxes; }

protected:
	int B_DIM;
	int N;

	// boxes are stacked vertically, going <x, y, width, height>
	mat boxes;
	ObsType obs_type;
	CostType cost_type;

	void init_dims(int T=10, int M=100, int N=1, double DT=1.0, int B_DIM=4, int X_DIM=2, int U_DIM=2, int Z_DIM=1, int Q_DIM=2, int R_DIM=1);
	void init(mat& boxes, const ObsType obs_type, const CostType cost_type, bool use_casadi,
			mat& xMin, mat& xMax, mat& uMin, mat& uMax, mat& R);

};

#endif
