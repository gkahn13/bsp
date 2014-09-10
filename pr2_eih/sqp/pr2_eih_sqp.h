#ifndef _PR2_EIH_SQP_H__
#define _PR2_EIH_SQP_H__

#include "../system/pr2_eih_system.h"
#include "../../util/logging.h"
#include "../../util/Timer.h"

extern "C" {
#include "../mpc/pr2eihMPC.h"
}

#include <stdio.h>

#include <boost/preprocessor/iteration/local.hpp>

namespace pr2_eih_sqp {

const int T = TIMESTEPS;

namespace cfg {
const double alpha_init = 1; // .1
const double alpha_gain = 3; // 3
const double alpha_epsilon = .1; // .1
const double alpha_max_increases = 5; // 5

const double Xeps_initial = .1; // .1
const double Ueps_initial = .1; // .1
const double improve_ratio_threshold = .1; // .1
const double min_approx_improve = .05; // .05
const double min_trust_box_size = 1e-1; // 1e-2
const double trust_shrink_ratio = .75; // .75
const double trust_expand_ratio = 1.25; // 1.25

const int max_iters = 100;
}

class PR2EihSqp {
public:
	PR2EihSqp();
	~PR2EihSqp();

	double collocation(StdVectorJ& J, StdVectorU& U, const MatrixJ& j_sigma0,
				const std::vector<Gaussian3d>& obj_gaussians,
				const std::vector<geometry3d::Triangle>& obstacles, PR2EihSystem& sys, bool plot=false);

private:
	pr2eihMPC_params problem;
	pr2eihMPC_output output;
	pr2eihMPC_info info;
	pr2eihMPC_FLOAT **H, **f, **lb, **ub, **z, *c;

	void setup_mpc_vars();
	void cleanup_mpc_vars();

	void print_inputs() const;
	bool is_valid_inputs() const;
	void L_BFGS(const StdVectorJ& J, const StdVectorU& U, const VectorTOTAL &grad,
			const StdVectorJ &Jopt, const StdVectorU &Uopt, const VectorTOTAL &gradopt,
			MatrixTOTAL &hess) const;
	double approximate_collocation(StdVectorJ& J, StdVectorU& U, const MatrixJ& j_sigma0,
			const std::vector<Gaussian3d>& obj_gaussians, const double alpha,
			const std::vector<geometry3d::Triangle>& obstacles, PR2EihSystem& sys, bool plot=false);
};

}

#endif
