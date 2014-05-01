#include "boxes-system.h"

/**
 *
 * Constructors
 *
 */


BoxesSystem::BoxesSystem() {
	this->init_dims();

	mat boxes(B_DIM, 1, fill::zeros);
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

	this->init(boxes, obs_type, cost_type, use_casadi, xMin, xMax, uMin, uMax, R);
}

BoxesSystem::BoxesSystem(mat& boxes, const ObsType obs_type, const CostType cost_type, bool use_casadi) {
	this->init_dims();

	mat xMin(X_DIM, 1, fill::ones), xMax(X_DIM, 1, fill::ones);
	mat uMin(U_DIM, 1, fill::ones), uMax(U_DIM, 1, fill::ones);

	xMin *= -1;
	xMax *= 6;
	uMin *= -.25;
	uMax *= .25;
	mat R = 1e-2*eye<mat>(N*R_DIM, N*R_DIM);

	this->init(boxes, obs_type, cost_type, use_casadi, xMin, xMax, uMin, uMax, R);
}

BoxesSystem::BoxesSystem(mat& boxes, const ObsType obs_type, const CostType cost_type, bool use_casadi,
								int T, int M, int N, double DT, int B_DIM, int X_DIM, int U_DIM, int Z_DIM, int Q_DIM, int R_DIM) {
	this->init_dims(T, M, N, DT, B_DIM, X_DIM, U_DIM, Z_DIM, Q_DIM, R_DIM);

	mat xMin(X_DIM, 1, fill::ones), xMax(X_DIM, 1, fill::ones);
	mat uMin(U_DIM, 1, fill::ones), uMax(U_DIM, 1, fill::ones);

	xMin *= -1;
	xMax *= 6;
	uMin *= -.25;
	uMax *= .25;
	mat R = 1e-2*eye<mat>(N*R_DIM, N*R_DIM);

	this->init(boxes, obs_type, cost_type, use_casadi, xMin, xMax, uMin, uMax, R);
}


void BoxesSystem::init_dims(int T, int M, int N, double DT, int B_DIM, int X_DIM, int U_DIM, int Z_DIM, int Q_DIM, int R_DIM) {
	this->T = T;
	this->M = M;
	this->N = N;
	this->DT = DT;
	this->B_DIM = B_DIM;
	this->X_DIM = X_DIM;
	this->U_DIM = U_DIM;
	this->Z_DIM = N*Z_DIM;
	this->Q_DIM = Q_DIM;
	this->R_DIM = R_DIM;

	this->TOTAL_VARS = T*N*X_DIM + (T-1)*N*U_DIM;
}

void BoxesSystem::init(mat& boxes, const ObsType obs_type, const CostType cost_type, bool use_casadi,
							  mat& xMin, mat& xMax, mat& uMin, mat& uMax, mat& R) {
	this->boxes = boxes;
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
	return NULL;
}

mat BoxesSystem::obsfunc(const mat& x, const mat& t, const mat& r) {
	return NULL;
}

void BoxesSystem::update_state_and_particles(const mat& x_t, const mat& P_t, const mat& u_t, mat& x_tp1, mat& P_tp1) {

}

double BoxesSystem::cost(const std::vector<mat>& X, const std::vector<mat>& U, const mat& P) {
	return 0;
}

void BoxesSystem::display_states_and_particles(const std::vector<mat>& X, const mat& P) {
	int M = P.n_cols; // in case use high-res particle set

	py::list x_list;
	for(int i=0; i < X_DIM; ++i) {
		for(int t=0; t < X.size(); ++t) {
			x_list.append(X[t](i));
		}
	}

	py::list boxes_list;
	for(int i=0; i < this->boxes.n_rows; ++i) {
		boxes_list.append(this->boxes(i));
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

		plot_state_and_particles(x_list, particles_list, boxes_list, X_DIM, M);

		LOG_INFO("Press enter to continue");
		py::exec("raw_input()",main_namespace);
	}
	catch(py::error_already_set const &)
	{
		PyErr_Print();
	}


}
