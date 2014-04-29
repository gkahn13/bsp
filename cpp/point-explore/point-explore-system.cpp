#include "point-explore-system.h"

/**
 * Constructors and initializers
 */

PointExploreSystem::PointExploreSystem() {
	mat::fixed <X_DIM,1> target(fill::zeros);
	ObsType obs_type = ObsType::angle;
	CostType cost_type = CostType::entropy;
	bool use_casadi = false;

	mat::fixed <X_DIM,1> xMin(fill::ones), xMax(fill::ones);
	mat::fixed <U_DIM,1> uMin(fill::ones), uMax(fill::ones);

	xMin *= -1;
	xMax *= 6;
	uMin *= -.25;
	uMax *= .25;
	mat R = 1e-2*eye<mat>(N*R_DIM, N*R_DIM);

	this->init(target, obs_type, cost_type, use_casadi, xMin, xMax, uMin, uMax, R);
}

PointExploreSystem::PointExploreSystem(mat& target, const ObsType obs_type, const CostType cost_type, bool use_casadi) {
	mat::fixed <X_DIM,1> xMin(fill::ones), xMax(fill::ones);
	mat::fixed <U_DIM,1> uMin(fill::ones), uMax(fill::ones);

	xMin *= -1;
	xMax *= 6;
	uMin *= -.25;
	uMax *= .25;
	mat R = 1e-2*eye<mat>(N*R_DIM, N*R_DIM);

	this->init(target, obs_type, cost_type, use_casadi, xMin, xMax, uMin, uMax, R);
}

PointExploreSystem::PointExploreSystem(mat& target, const ObsType obs_type, const CostType cost_type, bool use_casadi,
											mat& xMin, mat& xMax, mat& uMin, mat& uMax, mat& R) {
	this->init(target, obs_type, cost_type, use_casadi, xMin, xMax, uMin, uMax, R);
}

void PointExploreSystem::init(mat& target, const ObsType obs_type, const CostType cost_type, bool use_casadi,
							  mat& xMin, mat& xMax, mat& uMin, mat& uMax, mat& R) {
	this->target = target;
	this->obs_type = obs_type;
	this->cost_type = cost_type;
	this->use_casadi = use_casadi;
	this->xMin = xMin;
	this->xMax = xMax;
	this->uMin = uMin;
	this->uMax = uMax;
	this->R = R;

	if (this->use_casadi) {
		this->casadi_sys = new CasadiPointExploreSystem(this->obs_type, this->cost_type, this->R);
	}
}

/**
 *
 * Public methods
 *
 */

/**
 * x -- N*X_DIM by 1
 * u -- N*U_DIM by 1
 */
mat PointExploreSystem::dynfunc(const mat& x, const mat& u) {
	mat x_new = x + DT*u;
	return x_new;
}

/**
 * x -- N*X_DIM by 1
 * t -- N*X_DIM by 1
 * r -- N*R_DIM by 1
 */
mat PointExploreSystem::obsfunc(const mat& x, const mat& t, const mat& r) {
	if (this->obs_type == ObsType::angle) {
		return this->obsfunc_angle(x, t, r);
	} else {
		return this->obsfunc_dist(x, t, r);
	}
}

/**
 * x_t -- N*X_DIM by 1
 * P_t -- N*X_DIM by m
 * u_t -- N*U_DIM by 1
 * x_tp1 -- N*X_DIM by 1
 * P_tp1 -- N*X_DIM by m
 */
void PointExploreSystem::update_state_and_particles(const mat& x_t, const mat& P_t, const mat& u_t, mat& x_tp1, mat& P_tp1) {
	int M = P_t.n_cols;
	x_tp1 = this->dynfunc(x_t, u_t);

	// receive noisy measurement
	mat z_tp1 = this->obsfunc(x_tp1, this->target, sample_gaussian(zeros<mat>(N*R_DIM,1), .01*this->R));

	mat W(M, 1, fill::zeros);
	mat r(N*R_DIM, 1, fill::zeros);
	// for each particle, weight by gauss_likelihood of that measurement given particle/agent observation
	for(int m=0; m < M; ++m) {
		mat z_particle = this->obsfunc(x_tp1, P_t.col(m), r);
		mat e = z_particle - z_tp1;
		W(m) = this->gauss_likelihood(e, this->R);
	}
	W = W / accu(W);

	double sampling_noise = uniform(0, 1/double(M));
	P_tp1 = this->low_variance_sampler(P_t, W, sampling_noise);
}

double PointExploreSystem::cost(const std::vector<mat>& X, const std::vector<mat>& U, const mat& P) {
	if (this->use_casadi) {
		return this->casadi_sys->casadi_cost(X, U, P);
	}

	double cost;
	if (this->cost_type == CostType::entropy) {
		cost = this->cost_entropy(X, U, P);
	} else {
		cost = this->cost_platt(X, U, P);
	}

	for(int t=0; t < T-2; ++t) {
		cost += Constants::alpha_control_smooth*norm(U[t+1] - U[t], 2);
		cost += Constants::alpha_control_norm*norm(U[t], 2);
	}
	cost += Constants::alpha_control_norm*norm(U[T-2], 2);

	return cost;
}

mat PointExploreSystem::cost_grad(std::vector<mat>& X, std::vector<mat>& U, const mat& P) {
	if (this->use_casadi) {
		return this->casadi_sys->casadi_cost_grad(X, U, P);
	}

	mat g(TOTAL_VARS, 1, fill::zeros);

	double orig, entropy_p, entropy_l;
	int index = 0;
	for(int t=0; t < T; ++t) {
		for(int i=0; i < N*X_DIM; ++i) {
			orig = X[t](i);

			X[t](i) = orig + step;
			entropy_p = this->cost(X, U, P);

			X[t](i) = orig - step;
			entropy_l = this->cost(X, U, P);

			X[t][i] = orig;
			g(index++) = (entropy_p - entropy_l)/(2*step);
		}

		if (t < T-1) {
			for(int i=0; i < N*U_DIM; ++i) {
				orig = U[t](i);

				U[t](i) = orig + step;
				entropy_p = this->cost(X, U, P);

				U[t](i) = orig - step;
				entropy_l = this->cost(X, U, P);

				U[t][i] = orig;
				g(index++) = (entropy_p - entropy_l)/(2*step);
			}
		}

	}

	return g;
}


/**
 *
 * Private methods
 *
 */

mat PointExploreSystem::obsfunc_angle(const mat& x, const mat& t, const mat& r) {
	mat z(N*Z_DIM, 1, fill::zeros);
	for(int n=0; n < N; ++n) {
		mat x_n = x.rows(n*X_DIM, (n+1)*X_DIM-1);
		z(n) = atan((x_n(1)-t(1))/(x_n(0)-t(0))) + r(n);
	}
	return z;
}

mat PointExploreSystem::obsfunc_dist(const mat& x, const mat& t, const mat& r) {
	mat z(N*Z_DIM, 1, fill::zeros);
	mat::fixed<X_DIM, 1> x_n;
	for(int n=0; n < N; ++n) {
		x_n = x.rows(n*X_DIM, (n+1)*X_DIM-1);
		z(n) = 1.0/(1.0 + exp(-Constants::alpha*(Constants::max_range - norm(x_n-t,2)))) + r(n);
	}
	return z;
}

/**
 * S -- symmetric
 */
double PointExploreSystem::gauss_likelihood(const mat& v, const mat& S) {
	mat Sf = chol(S);
	mat M = solve(Sf, v);

	double E = -0.5*accu(M % M);
	double C = pow(2*M_PI, S.n_cols/2) * prod(diagvec(Sf));
	double w = exp(E) / C;

	return w;
}

/*
 * P -- X_DIM by M
 * W -- M by 1
 */
mat PointExploreSystem::low_variance_sampler(const mat& P, const mat& W, double r) {
	int M = P.n_cols;
	mat P_sampled(P.n_rows, P.n_cols, fill::zeros);

	double c = W(0);
	int i = 0;
	for(int m=0; m < M; ++m) {
		double u = r + m * (1/double(M));
		while (u > c) {
			c += W(++i);
		}
		P_sampled.col(m) = P.col(i);
	}

	return P_sampled;
}

double PointExploreSystem::cost_entropy(const std::vector<mat>& X, const std::vector<mat>& U, const mat& P) {
	int T = X.size();
	int M = P.n_cols;

	double entropy = 0;

	std::vector<mat> X_prop(T);
	std::vector<mat> H(T, zeros<mat>(N*Z_DIM, M));
	mat r(N*R_DIM, 1, fill::zeros);
	for(int t=0; t < T-1; ++t) {
		X_prop[t+1] = this->dynfunc(X[t], U[t]);
		for(int m=0; m < M; ++m) {
			H[t+1].col(m) = this->obsfunc(X_prop[t+1], P.col(m), r);
		}
	}

	std::vector<mat> W(T, zeros<mat>(M, 1));
	W[0] = (1/double(M))*ones<mat>(M,1);
	for(int t=1; t < T; ++t) {
		for(int m=0; m < M; ++m) {
			for(int p=0; p < M; ++p) {
				mat diff = H[t].col(m) - H[t].col(p);
				W[t](m) += this->gauss_likelihood(diff, this->R);
			}
		}
		W[t] = W[t] / accu(W[t]);

		// skoglar version
		// second term simplifies because zero particle dynamics
		entropy += accu(-W[t] % log(W[t])) +
				   accu(-W[t] % log(W[t-1])) +
				   log(accu(W[t-1] % W[t]));
	}

	return entropy;
}

double PointExploreSystem::cost_platt(const std::vector<mat>& X, const std::vector<mat>& U, const mat& P) {
	int T = X.size();
	int M = P.n_cols;

	std::vector<mat> X_prop(T, zeros<mat>(N*X_DIM,1));
	X_prop[0] = X[0];
	for(int t=0; t < T-1; ++t) {
		X_prop[t+1] = this->dynfunc(X[t], U[t]);
	}

	std::vector<mat> H(M, zeros<mat>(T*(N*Z_DIM),1));
	mat r(N*Z_DIM, 1, fill::zeros);
	for(int m=0; m < M; ++m) {
		for(int t=0; t < T; ++t) {
			H[m].rows(t*(N*Z_DIM), (t+1)*(N*Z_DIM)-1) = this->obsfunc(X_prop[t], P.col(m), r);
		}
	}

	double platt = 0;
	for(int m=1; m < M; ++m) {
		mat diff = H[m] - H[0];
		platt += (1/float(M-1))*exp(-accu(diff % diff));
	}

	return platt;
}

void PointExploreSystem::display_states_and_particles(const std::vector<mat>& X, const mat& P) {
	int M = P.n_cols; // in case use high-res particle set

	py::list x_list;
	for(int i=0; i < N*X_DIM; ++i) {
		for(int t=0; t < X.size(); ++t) {
			x_list.append(X[t](i));
		}
	}

	py::list targ_list;
	for(int i=0; i < this->target.n_rows; ++i) {
		targ_list.append(this->target(i));
	}

	py::list particles_list;
	for(int i=0; i < X_DIM; ++i) {
		for(int m=0; m < M; ++m) {
			particles_list.append(P(i,m));
		}
	}

	std::string workingDir = boost::filesystem::current_path().normalize().string();

	try
	{
		Py_Initialize();
		py::object main_module = py::import("__main__");
		py::object main_namespace = main_module.attr("__dict__");
		py::exec("import sys, os", main_namespace);
		py::exec(py::str("sys.path.append('"+workingDir+"')"), main_namespace);
		py::object plot_module = py::import("plot_point_explore");
		py::object plot_state_and_particles = plot_module.attr("plot_state_and_particles");

		plot_state_and_particles(x_list, particles_list, targ_list, X_DIM, M, N);

		LOG_INFO("Press enter to continue");
		py::exec("raw_input()",main_namespace);
	}
	catch(py::error_already_set const &)
	{
		PyErr_Print();
	}


}

