#include "../point-explore.h"
#include "casadi-point-explore.h"

namespace AD = CasADi;

int main(int argc, char* argv[]) {
//	srand(time(0));
	point_explore::initialize();

	std::vector<Matrix<X_DIM> > P(M);
	for(int m=0; m < M; ++m) {
		P[m][0] = uniform(0, 5);
		P[m][1] = uniform(0, 5);
	}

	Matrix<X_DIM> avg_particle = zeros<X_DIM,1>();
	for(int m=0; m < M; ++m) { avg_particle += (1/float(M))*P[m]; }
	Matrix<N*X_DIM> avg_particle_rep;
	for(int n=0; n < N; ++n) {
		avg_particle_rep.insert(n*X_DIM, 0, avg_particle);
	}

	Matrix<N*U_DIM> uinit = (avg_particle_rep - x0) / (DT*(T-1));
	std::vector<Matrix<N*U_DIM> > U(T-1, uinit);

	std::vector<Matrix<N*X_DIM> > X(T);
	X[0] = x0;
	for(int t=0; t < T-1; ++t) {
		X[t+1] = point_explore::dynfunc(X[t], U[t]);
	}

	double entropy = point_explore::differential_entropy(X,U,P);
	std::cout << "entropy: " << entropy << "\n";

	double casadi_entropy = point_explore::casadi_differential_entropy(X,U,P);
	std::cout << "casadi entropy: " << casadi_entropy << "\n";

//	Matrix<TOTAL_VARS> grad = point_explore::grad_differential_entropy(X,U,P);
//	Matrix<TOTAL_VARS> casadi_grad = point_explore::casadi_grad_differential_entropy(X,U,P);
//
//	Matrix<TOTAL_VARS> diff = grad - casadi_grad;
//	std::cout << "diff norm: " << tr(~diff*diff) << "\n";

	casadi_differential_entropy_func.generateCode("casadi/casadi_de_func.c");
	casadi_grad_differential_entropy_func.generateCode("casadi/casadi_grad_de_func.c");

}
