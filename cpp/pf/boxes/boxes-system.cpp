#include "boxes-system.h"

/**
 *
 * Constructors
 *
 */


BoxesSystem::BoxesSystem() {
	this->init_dims();

	mat box_centers(N*X_DIM, 1, fill::zeros);
	mat box_dims(N*X_DIM, 1, fill::ones);
	ObsType obs_type = ObsType::distance;
	CostType cost_type = CostType::entropy;
	bool use_casadi = false;

	mat xMin(X_DIM, 1, fill::ones), xMax(X_DIM, 1, fill::ones);
	mat uMin(U_DIM, 1, fill::ones), uMax(U_DIM, 1, fill::ones);

	xMin *= -1;
	xMax *= 6;
	uMin *= -.25;
	uMax *= .25;
	mat R = 1e-2*eye<mat>(N*R_DIM, N*R_DIM);

	this->init(box_centers, box_dims, obs_type, cost_type, use_casadi, xMin, xMax, uMin, uMax, R);
}

BoxesSystem::BoxesSystem(mat& box_centers, mat& box_dims, const ObsType obs_type, const CostType cost_type, bool use_casadi) {
	this->init_dims();

	mat xMin(X_DIM, 1, fill::ones), xMax(X_DIM, 1, fill::ones);
	mat uMin(U_DIM, 1, fill::ones), uMax(U_DIM, 1, fill::ones);

	xMin *= -1;
	xMax *= 6;
	uMin *= -.25;
	uMax *= .25;
	mat R = 1e-2*eye<mat>(N*R_DIM, N*R_DIM);

	this->init(box_centers, box_dims, obs_type, cost_type, use_casadi, xMin, xMax, uMin, uMax, R);
}

BoxesSystem::BoxesSystem(mat& box_centers, mat& box_dims, const ObsType obs_type, const CostType cost_type, bool use_casadi,
								int T, int M, int N, double DT, int X_DIM, int U_DIM, int Z_DIM, int Q_DIM, int R_DIM) {
	this->init_dims(T, M, N, DT, X_DIM, U_DIM, Z_DIM, Q_DIM, R_DIM);

	mat xMin(X_DIM, 1, fill::ones), xMax(X_DIM, 1, fill::ones);
	mat uMin(U_DIM, 1, fill::ones), uMax(U_DIM, 1, fill::ones);

	xMin *= -1;
	xMax *= 6;
	uMin *= -.25;
	uMax *= .25;
	mat R = 1e-2*eye<mat>(N*R_DIM, N*R_DIM);

	this->init(box_centers, box_dims, obs_type, cost_type, use_casadi, xMin, xMax, uMin, uMax, R);
}


void BoxesSystem::init_dims(int T, int M, int N, double DT, int X_DIM, int U_DIM, int Z_DIM, int Q_DIM, int R_DIM) {
	this->T = T;
	this->M = M;
	this->N = N;
	this->DT = DT;
	this->X_DIM = X_DIM;
	this->U_DIM = U_DIM;
	this->Z_DIM = N*Z_DIM;
	this->Q_DIM = Q_DIM;
	this->R_DIM = R_DIM;

	this->TOTAL_VARS = T*N*X_DIM + (T-1)*N*U_DIM;
}

void BoxesSystem::init(mat& box_centers, mat& box_dims, const ObsType obs_type, const CostType cost_type, bool use_casadi,
							  mat& xMin, mat& xMax, mat& uMin, mat& uMax, mat& R) {
	this->box_centers = box_centers;
	this->box_dims = box_dims;
	this->obs_type = obs_type;
	this->cost_type = cost_type;
	this->use_casadi = use_casadi;
	this->xMin = xMin;
	this->xMax = xMax;
	this->uMin = uMin;
	this->uMax = uMax;
	this->R = R;

	if (this->use_casadi) {
//		this->casadi_sys = new CasadiExploreSystem(this->obs_type, this->cost_type, this->R,
//				T, M, N, DT, X_DIM, U_DIM, Z_DIM/N, Q_DIM, R_DIM);
	}
}

/**
 *
 * Public methods
 *
 */

mat BoxesSystem::dynfunc(const mat& x, const mat& u) {
	mat x_new = x + DT*u;
	return x_new;
}

mat BoxesSystem::obsfunc(const mat& x, const mat& b_centers, const mat& r) {
	mat z(Z_DIM, 1, fill::zeros);
	double x_eye = x(0);
	double y_eye = x(1);
	for(int n=0; n < N; ++n) {
		double bx = b_centers(n*X_DIM);
		double by = b_centers(n*X_DIM+1);

		double width = this->box_dims(n*X_DIM);
		double height = this->box_dims(n*X_DIM+1);

		if (fabs(y_eye-by) < height/2.0) {
			z(n) = sqrt((x_eye - (bx+width/2.0))*(x_eye - (bx+width/2.0))) + r(n);
		} else {
			z(n) = 0 + r(n);
		}
	}

	return z;
}

void BoxesSystem::update_state_and_particles(const mat& x_t, const mat& P_t, const mat& u_t, mat& x_tp1, mat& P_tp1) {
	int M = P_t.n_cols;
	x_tp1 = this->dynfunc(x_t, u_t);

	// receive noisy measurement
	mat z_tp1 = this->obsfunc(x_tp1, this->box_centers, sample_gaussian(zeros<mat>(N*R_DIM,1), this->R));

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

double BoxesSystem::cost(const std::vector<mat>& X, const std::vector<mat>& U, const mat& P) {
	if (this->use_casadi) {
		return this->casadi_sys->casadi_cost(X, U, P);
	}

	double cost;
	if (this->cost_type == CostType::entropy) {
		cost = this->cost_entropy(X, U, P);
	} else {
		cost = this->cost_entropy(X, U, P);
	}

	return cost;
}

void BoxesSystem::display_states_and_particles(const std::vector<mat>& X, const mat& P) {
	int M = P.n_cols; // in case use high-res particle set

	py::list x_list;
	for(int i=0; i < X_DIM; ++i) {
		for(int t=0; t < X.size(); ++t) {
			x_list.append(X[t](i));
		}
	}

	py::list box_centers_list, box_dims_list;
	for(int i=0; i < this->box_centers.n_rows; ++i) {
		box_centers_list.append(this->box_centers(i));
		box_dims_list.append(this->box_dims(i));
	}

	py::list particles_list;
	for(int i=0; i < P.n_rows; ++i) {
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
		py::exec(py::str("sys.path.append('"+workingDir+"/boxes')"), main_namespace);
		py::object plot_module = py::import("plot_boxes");
		py::object plot_state_and_particles = plot_module.attr("plot_state_and_particles");

		plot_state_and_particles(x_list, particles_list, box_centers_list, box_dims_list, X_DIM, M);

		LOG_INFO("Press enter to continue");
		py::exec("raw_input()",main_namespace);
	}
	catch(py::error_already_set const &)
	{
		PyErr_Print();
	}


}
