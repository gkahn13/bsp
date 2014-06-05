#include "./include/eih_system.h"
#include "./include/pr2_sim.h"
#include "./include/utils.h"

#include <armadillo>
using namespace arma;

#include "../../util/Timer.h"
#include "../../util/logging.h"

#include <boost/program_options.hpp>
namespace po = boost::program_options;

extern "C" {
#include "eihMPC.h"
eihMPC_FLOAT **H, **f, **lb, **ub, **z, **c;
}

#define TIMESTEPS 5
#define DT 1.0 // Note: if you change this, must change the FORCES matlab file
#define X_DIM 7
#define U_DIM 7

const int T = TIMESTEPS;
int M;
const int TOTAL_VARS = T*X_DIM + (T-1)*U_DIM;

mat::fixed<X_DIM, 1> x0;

namespace cfg {
const double improve_ratio_threshold = .1;
const double min_approx_improve = 1e-4;
const double min_trust_box_size = 1e-4;
const double trust_shrink_ratio = .5;
const double trust_expand_ratio = 1.5;
}

void setupMPCVars(eihMPC_params& problem, eihMPC_output& output) {
	// inputs
	H = new eihMPC_FLOAT*[T];
	f = new eihMPC_FLOAT*[T];
	lb = new eihMPC_FLOAT*[T];
	ub = new eihMPC_FLOAT*[T];
	c = new eihMPC_FLOAT*[1];

	// output
	z = new eihMPC_FLOAT*[T];

#define SET_VARS(n)    \
		H[ BOOST_PP_SUB(n,1) ] = problem.H##n ;  \
		f[ BOOST_PP_SUB(n,1) ] = problem.f##n ;  \
		lb[ BOOST_PP_SUB(n,1) ] = problem.lb##n ;	\
		ub[ BOOST_PP_SUB(n,1) ] = problem.ub##n ;	\
		z[ BOOST_PP_SUB(n,1) ] = output.z##n ;

#define BOOST_PP_LOCAL_MACRO(n) SET_VARS(n)
#define BOOST_PP_LOCAL_LIMITS (1, TIMESTEPS)
#include BOOST_PP_LOCAL_ITERATE()

	c[0] = problem.c1;

	for(int t=0; t < T-1; ++t) {
		for(int i=0; i < (X_DIM+U_DIM); ++i) { H[t][i] = INFTY; }
		for(int i=0; i < (X_DIM+U_DIM); ++i) { f[t][i] = INFTY; }
		for(int i=0; i < (X_DIM+U_DIM); ++i) { lb[t][i] = INFTY; }
		for(int i=0; i < (X_DIM+U_DIM); ++i) { ub[t][i] = INFTY; }
		for(int i=0; i < (X_DIM+U_DIM); ++i) { z[t][i] = INFTY; }
	}
	for(int i=0; i < X_DIM; ++i) { H[T-1][i] = INFTY; }
	for(int i=0; i < X_DIM; ++i) { f[T-1][i] = INFTY; }
	for(int i=0; i < X_DIM; ++i) { lb[T-1][i] = INFTY; }
	for(int i=0; i < X_DIM; ++i) { ub[T-1][i] = INFTY; }
	for(int i=0; i < X_DIM; ++i) { z[T-1][i] = INFTY; }

	for(int i=0; i < X_DIM; ++i) { c[0][i] = INFTY; }
}

void cleanupMPCVars() {
	delete[] H;
	delete[] f;
	delete[] lb;
	delete[] ub;
	delete[] z;
	delete[] c;
}

bool isValidInputs()
{
	for(int t = 0; t < T-1; ++t) {
		std::cout << "\n\nt: " << t << "\n";

		if (t == 0) {
			std::cout << "\nc[0]:\n";
			for(int i=0; i < (X_DIM); ++i) {
				std::cout << c[0][i] << " ";
			}
		}

		std::cout << "\nH[" << t << "]: ";
		for(int i=0; i < (X_DIM+U_DIM); ++i) {
			std::cout << H[t][i] << " ";
		}

		std::cout << "\nf[" << t << "]: ";
		for(int i=0; i < (X_DIM+U_DIM); ++i) {
			std::cout << f[t][i] << " ";
		}

		std::cout << "\nlb[" << t << "]: ";
		for(int i=0; i < (X_DIM+U_DIM); ++i) {
			std::cout << lb[t][i] << " ";
		}

		std::cout << "\nub[" << t << "]: ";
		for(int i=0; i < (X_DIM+U_DIM); ++i) {
			std::cout << ub[t][i] << " ";
		}
	}
	std::cout << "\n\nt: " << T-1 << "\n";

	std::cout << "\nH[" << T-1 << "]: ";
	for(int i=0; i < (X_DIM); ++i) {
		std::cout << H[T-1][i] << " ";
	}

	std::cout << "\nf[" << T-1 << "]: ";
	for(int i=0; i < (X_DIM); ++i) {
		std::cout << f[T-1][i] << " ";
	}

	std::cout << "\nlb[" << T-1 << "]: ";
	for(int i=0; i < (X_DIM); ++i) {
		std::cout << lb[T-1][i] << " ";
	}

	std::cout << "\nub[" << T-1 << "]: ";
	for(int i=0; i < (X_DIM); ++i) {
		std::cout << ub[T-1][i] << " ";
	}

	std::cout << "\n";

	for(int t = 0; t < T-1; ++t) {
		for(int i=0; i < (X_DIM+U_DIM); ++i) { if (H[t][i] > INFTY/2) { return false; } }
		for(int i=0; i < (X_DIM+U_DIM); ++i) { if (f[t][i] > INFTY/2) { return false; } }
		for(int i=0; i < (X_DIM+U_DIM); ++i) { if (lb[t][i] > INFTY/2) { return false; } }
		for(int i=0; i < (X_DIM+U_DIM); ++i) {if (ub[t][i] > INFTY/2) { return false; } }
	}
	for(int i=0; i < (X_DIM); ++i) { if (H[T-1][i] > INFTY/2) { return false; } }
	for(int i=0; i < (X_DIM); ++i) { if (f[T-1][i] > INFTY/2) { return false; } }
	for(int i=0; i < (X_DIM); ++i) { if (lb[T-1][i] > INFTY/2) { return false; } }
	for(int i=0; i < (X_DIM); ++i) { if (ub[T-1][i] > INFTY/2) { return false; } }

	for(int i=0; i < (X_DIM); ++i) { if (c[0][i] > INFTY/2) { return false; } }

	return true;
}

double eihCollocation(std::vector<mat> &X, std::vector<mat> &U, mat &P, EihSystem *sys,
		eihMPC_params &problem, eihMPC_output &output, eihMPC_info &info) {

	int max_iter = 100;
	double Xeps = .5;
	double Ueps = .5;

	double merit = 0;
	double constant_cost, hessian_constant, jac_constant;
	mat g(TOTAL_VARS, 1, fill::zeros), diaghess(TOTAL_VARS, 1, fill::zeros);

	std::vector<mat> Xopt(T, zeros<mat>(X_DIM, 1));
	std::vector<mat> Uopt(T-1, zeros<mat>(U_DIM, 1));
	double optcost, model_merit, new_merit;
	double approx_merit_improve, exact_merit_improve, merit_improve_ratio;

	LOG_DEBUG("Initial trajectory cost: %4.10f", sys->cost(X, U, P));

	int index = 0;
	bool solution_accepted = true;
	for(int it=0; it < max_iter; ++it) {

		LOG_DEBUG("\nIter: %d", it);

		// only compute gradient/hessian if P/U has been changed
		if (solution_accepted) {
			g = sys->cost_grad(X, U, P);

			diaghess.zeros();
			merit = sys->cost(X, U, P);

			constant_cost = 0;
			hessian_constant = 0;
			jac_constant = 0;

			// fill in Hessian first so we can force it to be PSD
			// TODO: use finite differences or BFGS to approximate. (else set to 0)
			index = 0;
			for(int t=0; t < T-1; ++t) {
				for(int i=0; i < (X_DIM+U_DIM); ++i) {
					double val = diaghess(index++);
					H[t][i] = (val < 0) ? 0 : val;
				}
			}
			for(int i=0; i < X_DIM; ++i) {
				double val = diaghess(index++);
				H[T-1][i] = (val < 0) ? 0 : val;
			}

			// fill in gradient
			index = 0;
			for(int t=0; t < T-1; ++t) {
				mat::fixed<(X_DIM+U_DIM), 1> zbar;
				zbar.rows(0, X_DIM-1) = X[t];
				zbar.rows(X_DIM, (X_DIM+U_DIM)-1) = U[t];

				for(int i=0; i < (X_DIM+U_DIM); ++i) {
					hessian_constant += H[t][i]*zbar(i)*zbar(i);
					jac_constant -= g[index]*zbar(i);
					f[t][i] = g[index] - H[t][i]*zbar(i);
					index++;
				}
			}

			mat zbar = X[T-1];

			for(int i=0; i < X_DIM; ++i) {
				hessian_constant += H[T-1][i]*zbar(i)*zbar(i);
				jac_constant -= g[index]*zbar(i);
				f[T-1][i] = g[index] - H[T-1][i]*zbar(i);
				index++;
			}

			for(int i=0; i < X_DIM; ++i) {
				c[0][i] = X[0](i);
			}

			constant_cost = 0.5*hessian_constant + jac_constant + merit;
		}


		mat xMin = sys->get_xMin();
		mat xMax = sys->get_xMax();
		mat uMin = sys->get_uMin();
		mat uMax = sys->get_uMax();

		// set trust region bounds based on current trust region size
		for(int t=0; t < T; ++t) {
			index = 0;
			for(int i=0; i < X_DIM; ++i) {
				lb[t][index] = MAX(xMin(i), X[t](i) - Xeps);
				ub[t][index] = MIN(xMax(i), X[t](i) + Xeps);
				index++;
			}


			if (t < T-1) {
				// set each input lower/upper bound
				for(int i=0; i < U_DIM; ++i) {
					lb[t][index] = MAX(uMin(i), U[t](i) - Ueps);
					ub[t][index] = MIN(uMax(i), U[t](i) + Ueps);
					index++;
				}
			}

		}


		// Verify problem inputs
//		if (!isValidInputs()) {
//			LOG_ERROR("Inputs are not valid!");
//			exit(0);
//		}


		// call FORCES
		int exitflag = eihMPC_solve(&problem, &output, &info);
		if (exitflag == 1) {
			optcost = info.pobj;
			for(int t=0; t < T; ++t) {
				index = 0;
				for(int i=0; i < X_DIM; ++i) {
					Xopt[t](i) = z[t][index++];
				}

				if (t < T-1) {
					for(int i=0; i < U_DIM; ++i) {
						Uopt[t](i) = z[t][index++];
					}
				}
			}
		} else {
			LOG_FATAL("Some problem in solver");
			exit(-1);
		}

		std::cout << "Displaying optimized path\n";
		sys->display_states_and_particles(Xopt, P, true);


		model_merit = optcost + constant_cost; // need to add constant terms that were dropped

		new_merit = sys->cost(Xopt, Uopt, P);

		LOG_DEBUG("merit: %f", merit);
		LOG_DEBUG("model_merit: %f", model_merit);
		LOG_DEBUG("new_merit: %f", new_merit);
		LOG_DEBUG("constant cost term: %f", constant_cost);

		approx_merit_improve = merit - model_merit;
		exact_merit_improve = merit - new_merit;
		merit_improve_ratio = exact_merit_improve / approx_merit_improve;

		LOG_DEBUG("approx_merit_improve: %f", approx_merit_improve);
		LOG_DEBUG("exact_merit_improve: %f", exact_merit_improve);
		LOG_DEBUG("merit_improve_ratio: %f", merit_improve_ratio);

		if (approx_merit_improve < -1e-5) {
			LOG_ERROR("Approximate merit function got worse: %f", approx_merit_improve);
			LOG_ERROR("Failure!");
			return INFTY;
		} else if (approx_merit_improve < cfg::min_approx_improve) {
			LOG_DEBUG("Converged: improvement small enough");
			X = Xopt; U = Uopt;
			solution_accepted = true;
			break;
		} else if ((exact_merit_improve < 0) || (merit_improve_ratio < cfg::improve_ratio_threshold)) {
			Xeps *= cfg::trust_shrink_ratio;
			Ueps *= cfg::trust_shrink_ratio;
			LOG_DEBUG("Shrinking trust region size to: %2.6f %2.6f", Xeps, Ueps);
			solution_accepted = false;
		} else {
			// expand Xeps and Ueps
			Xeps *= cfg::trust_expand_ratio;
			Ueps *= cfg::trust_expand_ratio;
			LOG_DEBUG("Accepted, Increasing trust region size to:  %2.6f %2.6f", Xeps, Ueps);
			X = Xopt; U = Uopt;
			solution_accepted = true;
		}


	}

	return sys->cost(X, U, P);
}

void setup_eih_environment(PR2 *brett, Arm::ArmType arm_type, bool zero_seed, mat &P, EihSystem **sys) {
	if (zero_seed) {
		srand(time(0));
	}

	rave::EnvironmentBasePtr env = brett->get_env();

	Arm *larm = brett->larm;
	Arm *rarm = brett->rarm;
	KinectSensor *l_kinect = brett->l_kinect;
	KinectSensor *r_kinect = brett->r_kinect;

	larm->set_posture(Arm::Posture::mantis);
	rarm->set_posture(Arm::Posture::mantis);

	rave::KinBodyPtr table = env->GetKinBody("table");
	rave::KinBody::LinkPtr base = table->GetLink("base");
	rave::Vector extents = base->GetGeometry(0)->GetBoxExtents();

	rave::Vector table_pos = table->GetTransform().trans;
	double x_min, x_max, y_min, y_max, z_min, z_max;
	x_min = table_pos.x - extents.x;
	x_max = table_pos.x + extents.x;
	y_min = table_pos.y - extents.y;
	y_max = table_pos.y + extents.y;
	z_min = table_pos.z + extents.z;
	z_max = table_pos.z + extents.z + .2;

	for(int m=0; m < M; ++m) {
		P(0,m) = uniform(x_min, x_max);
		P(1,m) = uniform(y_min, y_max);
		P(2,m) = uniform(z_min, z_max);
	}

	Manipulator *manip;
	KinectSensor *kinect;
	if (arm_type == Arm::ArmType::left) {
		manip = larm;
		kinect = l_kinect;
	} else {
		manip = rarm;
		kinect = r_kinect;
	}
	*sys = new EihSystem(env, manip, kinect);
	kinect->render_on();
	boost::this_thread::sleep(boost::posix_time::seconds(2));
}

int main(int argc, char* argv[]) {
	PR2* brett = new PR2();
	M = 1000;
	mat P(3, M);
	EihSystem *sys;
	setup_eih_environment(brett, Arm::ArmType::right, false, P, &sys);
	Arm *arm = brett->rarm;

	arm->set_posture(Arm::Posture::mantis);
	x0 = arm->get_joint_values();

	std::vector<mat> U(T-1, zeros<mat>(U_DIM, 1));
	std::vector<mat> X(T, x0);

	double init_cost = sys->cost(X, U, P);
	LOG_DEBUG("Initial cost: %4.10f", init_cost);

	// initialize FORCES variables
	eihMPC_params problem;
	eihMPC_output output;
	eihMPC_info info;

	setupMPCVars(problem, output);
	util::Timer forces_timer;
	std::vector<mat> X_actual(1, x0);

	do {
		init_cost = sys->cost(X, U, P);

		util::Timer_tic(&forces_timer);
		double cost = eihCollocation(X, U, P, sys, problem, output, info);
		double forces_time = util::Timer_toc(&forces_timer);

		LOG_DEBUG("Initial cost: %4.10f", init_cost);
		LOG_DEBUG("Cost: %4.10f", cost);
		LOG_DEBUG("Time: %4.10f ms", forces_time*1000);

		std::cout << "U\n";
		for(int t=0; t < T-1; ++t) {
			std::cout << U[t].t();
		}
		std::cout << "\n";

		X[0] = x0;
		std::cout << X[0].t();
		for(int t=0; t < T-1; ++t) {
			X[t+1] = sys->dynfunc(X[t], U[t]);
			std::cout << X[t+1].t();
		}
		std::cout << "\n";

		LOG_DEBUG("Optimized path");
		sys->display_states_and_particles(X, P, true);
		break;

		std::cout << "\nPress 'q' to exit\n";
	} while (utils::getch() != 'q');

}

