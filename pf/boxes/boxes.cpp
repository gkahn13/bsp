#include "boxes-system.h"

#include <iostream>
#include <vector>

#include "../../util/Timer.h"
#include "../../util/logging.h"

#include <boost/program_options.hpp>
namespace po = boost::program_options;

extern "C" {
#include "boxesMPC.h"
boxesMPC_FLOAT **H, **f, **lb, **ub, **z, **c;
}

#define TIMESTEPS 10
#define BOXES 1
#define DT 1.0 // Note: if you change this, must change the FORCES matlab file
#define X_DIM 2
#define U_DIM 2
#define Z_DIM 1
#define Q_DIM 2
#define R_DIM 1

const int T = TIMESTEPS;
const int N = 1;
int M;
const int TOTAL_VARS = T*X_DIM + (T-1)*U_DIM;


mat::fixed<X_DIM, 1> x0;

namespace cfg {
const double improve_ratio_threshold = .1;
const double min_approx_improve = 1e-4;
const double min_trust_box_size = 1e-4;
const double trust_shrink_ratio = .75;
const double trust_expand_ratio = 1.25;
}

void setupMPCVars(boxesMPC_params& problem, boxesMPC_output& output) {
	// inputs
	H = new boxesMPC_FLOAT*[T];
	f = new boxesMPC_FLOAT*[T];
	lb = new boxesMPC_FLOAT*[T];
	ub = new boxesMPC_FLOAT*[T];
	c = new boxesMPC_FLOAT*[1];

	// output
	z = new boxesMPC_FLOAT*[T];

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

double boxesCollocation(std::vector<mat>& X, std::vector<mat>& U, mat& P, BoxesSystem& sys,
		boxesMPC_params& problem, boxesMPC_output& output, boxesMPC_info& info) {

	int max_iter = 100;
	double Xeps = .5;
	double Ueps = .5;

	double merit = 0;
	double constant_cost, hessian_constant, jac_constant;
	mat d(TOTAL_VARS, 1, fill::zeros), diaghess(TOTAL_VARS, 1, fill::zeros);

	std::vector<mat> Xopt(T, zeros<mat>(X_DIM, 1));
	std::vector<mat> Uopt(T-1, zeros<mat>(U_DIM, 1));
	double optcost, model_merit, new_merit;
	double approx_merit_improve, exact_merit_improve, merit_improve_ratio;

	LOG_DEBUG("Initial trajectory cost: %4.10f", sys.cost(X, U, P));

	int index = 0;
	bool solution_accepted = true;
	for(int it=0; it < max_iter; ++it) {

		LOG_DEBUG("\nIter: %d", it);

		// only compute gradient/hessian if P/U has been changed
		if (solution_accepted) {
			d = sys.cost_grad(X, U, P);

//			diaghess = point_explore::casadi_diaghess_differential_entropy(X,U,P);
			diaghess.zeros();
			merit = sys.cost(X, U, P);

			constant_cost = 0;
			hessian_constant = 0;
			jac_constant = 0;

			// compute Hessian first so we can force it to be PSD
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

			// compute gradient
			index = 0;
			for(int t=0; t < T-1; ++t) {
				mat::fixed<(X_DIM+U_DIM), 1> zbar;
				zbar.rows(0, X_DIM-1) = X[t];
				zbar.rows(X_DIM, (X_DIM+U_DIM)-1) = U[t];

				for(int i=0; i < (X_DIM+U_DIM); ++i) {
					hessian_constant += H[t][i]*zbar(i)*zbar(i);
					jac_constant -= d[index]*zbar(i);
					f[t][i] = d[index] - H[t][i]*zbar(i);
					index++;
				}
			}

			mat zbar = X[T-1];

			for(int i=0; i < X_DIM; ++i) {
				hessian_constant += H[T-1][i]*zbar(i)*zbar(i);
				jac_constant -= d[index]*zbar(i);
				f[T-1][i] = d[index] - H[T-1][i]*zbar(i);
				index++;
			}

			for(int i=0; i < X_DIM; ++i) {
				c[0][i] = X[0](i);
			}

			constant_cost = 0.5*hessian_constant + jac_constant + merit;
		}


		mat xMin = sys.get_xMin();
		mat xMax = sys.get_xMax();
		mat uMin = sys.get_uMin();
		mat uMax = sys.get_uMax();

		// set trust region bounds based on current trust region size
		for(int t=0; t < T; ++t) {
			// set each particle lower/upper bound
			index = 0;
			for(int n=0; n < N; ++n) {
				for(int i=0; i < X_DIM; ++i) {
					lb[t][index] = MAX(xMin(i), X[t](n*X_DIM+i) - Xeps);
					ub[t][index] = MIN(xMax(i), X[t](n*X_DIM+i) + Xeps);
					index++;
				}
			}

			if (t < T-1) {
				// set each input lower/upper bound
				for(int n=0; n < N; ++n) {
					for(int i=0; i < U_DIM; ++i) {
						lb[t][index] = MAX(uMin(i), U[t](n*X_DIM+i) - Ueps);
						ub[t][index] = MIN(uMax(i), U[t](n*X_DIM+i) + Ueps);
						index++;
					}
				}
			}
		}


		// Verify problem inputs
//		if (!isValidInputs()) {
//			LOG_ERROR("Inputs are not valid!");
//			exit(0);
//		}

		// call FORCES
		int exitflag = boxesMPC_solve(&problem, &output, &info);
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

		model_merit = optcost + constant_cost; // need to add constant terms that were dropped

		new_merit = sys.cost(Xopt, Uopt, P);

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

	return sys.cost(X, U, P);
}

enum class InitType { average, zero };
// assume x0 is set before this
void initialize_trajectory(std::vector<mat>& X, std::vector<mat>& U, const mat& P, InitType init_type, BoxesSystem& sys) {
	int M = P.n_cols;

	mat uinit(U_DIM, 1, fill::zeros);

	if (init_type == InitType::average) {
		// go to average of particles
		mat avg_particle(X_DIM, 1, fill::zeros);
		for(int m=0; m < M; ++m) { avg_particle += (1/float(M))*P.col(m); }
		mat avg_particle_rep = repmat(avg_particle, N, 1);

		uinit = (avg_particle_rep - x0) / (DT*(T-1));
	}
//	else if (type == 1) {
//		// go to furthest heaviest particle
//		double eps = 1e-2;
//		std::vector<int> num_particles_nearby(M, 0);
//		int max_num_particles_nearby = -INFTY;
//		for(int m=0; m < M; ++m) {
//			for(int n=0; n < M; ++n) {
//				double d = dist<X_DIM>(P[m],P[n]);
//				if (d < eps) {
//					num_particles_nearby[m]++;
//				}
//			}
//			max_num_particles_nearby = MAX(max_num_particles_nearby, num_particles_nearby[m]);
//		}
//
//		Matrix<X_DIM> furthest_heaviest_particle;
//		double furthest = -INFTY;
//
//		for(int m=0; m < M; ++m) {
//			if (num_particles_nearby[m] == max_num_particles_nearby) {
//				double d = dist<X_DIM>(P[m], x0);
//				if (d > furthest) {
//					furthest = d;
//					furthest_heaviest_particle = P[m];
//				}
//			}
//		}
//
//		uinit = (furthest_heaviest_particle - x0) / (DT*(T-1));
	else if (init_type == InitType::zero) {
		// zero initialization
		uinit = zeros<mat>(U_DIM, 1);
	}

//	uinit << -.05 << endr << -.1;

	mat uMin = sys.get_uMin();
	mat uMax = sys.get_uMax();

	for(int i=0; i < U_DIM; ++i) {
		uinit(i) = (uinit(i) > uMax(i)) ? uMax(i) : uinit(i);
		uinit(i) = (uinit(i) < uMin(i)) ? uMin(i) : uinit(i);
	}


	X[0] = x0;
	for(int t=0; t < T-1; ++t) {
		U[t] = uinit;
		X[t+1] = sys.dynfunc(X[t], U[t]);
	}

}

// TODO: hardcode for one box
bool found(const mat& P, mat& dims) {
	int M = P.n_cols;

	double max_spread = -INFTY;
	for(int m=0; m < M; ++m) {
		for(int p=m; p < M; ++p) {
			max_spread = MAX(max_spread, norm(P.col(m) - P.col(p), 2));
		}
	}

	return (max_spread < dims(0) + dims(1));
}

void parse_boxes(int argc, char* argv[], ObsType& obs_type, CostType& cost_type, bool& use_casadi, mat& box_centers, mat& box_dims) {
	ObsTypeList obs_list;
	CostTypeList cost_list;
	std::vector<double> centers, dims;

	// Declare the supported options.
	po::options_description desc("Allowed options");
	desc.add_options()
		    				("help", "produce help message")
		    				("M", po::value<int>(&M), "Number of particles (default 100)")
		    				("obs", po::value<ObsTypeList>(&obs_list)->multitoken(), "Observation type <angle> or <distance> (default is <angle>)")
		    				("cost", po::value<CostTypeList>(&cost_list)->multitoken(), "Cost type <entropy> or <platt> (default is <entropy>)")
		    				("casadi", po::value<bool>(&use_casadi), "Use CasADi or not")
		    				("centers", po::value<std::vector<double> >(&centers)->multitoken(), "for N boxes, 2*N length of (x,y) coordinates")
		    				("dims", po::value<std::vector<double> >(&dims)->multitoken(), "for N boxes, 2*N length of (width, height)")
		    				;

	try {
		po::variables_map vm;
		po::store(po::parse_command_line(argc, argv, desc, po::command_line_style::unix_style ^ po::command_line_style::allow_short), vm);
		po::notify(vm);

		if (vm.count("help")) {
			cout << desc << "\n";
			exit(0);
		}

		if (vm.count("obs")) {
			obs_type = obs_list[0];
		}

		if (vm.count("cost")) {
			cost_type = cost_list[0];
		}

		if (vm.count("centers")) {
			if (centers.size() != N*X_DIM) {
				throw std::runtime_error(std::string("centers size incorrect"));
			}
			mat centers_mat(centers.data(), N*X_DIM, 1);
			box_centers = centers_mat;
		}

		if (vm.count("dims")) {
			if (dims.size() != N*X_DIM) {
				throw std::runtime_error(std::string("dims size incorrect"));
			}
			mat dims_mat(dims.data(), N*X_DIM, 1);
			box_dims = dims_mat;
		}


	} catch (std::exception &e) {
		std::cerr << "error: " << e.what() << "\n";
		exit(0);
	}

}

int main(int argc, char* argv[]) {
	srand(time(0));

	ObsType obs_type = ObsType::distance;
	CostType cost_type = CostType::entropy;
	bool use_casadi = true;
	M = 100;
	mat box_centers(N*X_DIM, 1, fill::zeros);
	mat box_dims(N*X_DIM, 1, fill::zeros);
	box_centers << -.25 << endr << 0;
	box_dims << .5 << endr << .25;

	parse_boxes(argc, argv, obs_type, cost_type, use_casadi, box_centers, box_dims);

	x0 << 1.5 << endr << 0;

	LOG_DEBUG("Initializing...");
	BoxesSystem sys = BoxesSystem(box_centers, box_dims, obs_type, cost_type, use_casadi,
			T, M, N, DT, X_DIM, U_DIM, Z_DIM, Q_DIM, R_DIM);
	LOG_DEBUG("System initialized");

	const int M_FULL = 1000;
	mat P_full(X_DIM, M_FULL, fill::zeros);
	for(int m=0; m < M_FULL; ++m) {
		P_full(0, m) = uniform(-2, 0.5);
		P_full(1, m) = uniform(-2, 2);
	}

	std::vector<mat> U(T-1);
	std::vector<mat> X(T);

	InitType init_type = InitType::zero;
	initialize_trajectory(X, U, P_full, init_type, sys);

	mat P = subsample(P_full, M);

	double init_cost = sys.cost(X,U,P);
	LOG_DEBUG("Initial cost: %4.10f", init_cost);

	LOG_DEBUG("Display initial trajectory");
	sys.display_states_and_particles(X, P_full);

	// initialize FORCES variables
	boxesMPC_params problem;
	boxesMPC_output output;
	boxesMPC_info info;

	setupMPCVars(problem, output);
	util::Timer forces_timer;
	std::vector<mat> X_actual(1, x0);
	while(true) {

		P = subsample(P_full, M);
		init_cost = sys.cost(X,U,P);

		util::Timer_tic(&forces_timer);
		double cost = boxesCollocation(X, U, P, sys, problem, output, info);
		double forces_time = util::Timer_toc(&forces_timer);

		LOG_DEBUG("Initial cost: %4.10f", init_cost);
		LOG_DEBUG("Cost: %4.10f", cost);
		LOG_DEBUG("Time: %4.10f ms", forces_time*1000);

		LOG_DEBUG("Optimized path");
		sys.display_states_and_particles(X, P, true);

		mat x = X[0], x_tp1(N*X_DIM, 1, fill::zeros);
		mat P_full_tp1(X_DIM, M_FULL, fill::zeros);
		int num_execute = 1;
		for(int t=0; t < num_execute; ++t) {
			sys.update_state_and_particles(x, P_full, U[t], x_tp1, P_full_tp1);
			P_full = P_full_tp1;
			x = x_tp1;

			X_actual.push_back(x);
		}

		if (found(P_full, box_dims)) {
			break;
		}

		x0 = x_tp1;
		initialize_trajectory(X, U, P_full, init_type, sys);

		LOG_DEBUG("Particle update step");
		sys.display_states_and_particles(X, P_full, true);
	}

	LOG_DEBUG("Found the landmark");

	double total_dist= 0;
	for(int t=0; t < X_actual.size()-1; ++t) {
		total_dist += norm(X_actual[t+1] - X_actual[t], 2);
	}
//	LOG_INFO("Total distance traveled: %4.10f", total_dist);
	std::cout << total_dist << ",\n";

	sys.display_states_and_particles(X_actual, P_full);


}
