#include "../point-explore.h"
#include "casadi-point-explore.h"

namespace AD = CasADi;

int main(int argc, char* argv[]) {
//	srand(time(0));
	point_explore::initialize();

	x0[0] = 0; x0[1] = 0;

	std::vector<Matrix<X_DIM> > P(M);
	for(int m=0; m < M; ++m) {
		P[m][0] = uniform(0, 5);
		P[m][1] = uniform(0, 5);
	}

	Matrix<X_DIM> avg_particle = zeros<X_DIM,1>();
	for(int m=0; m < M; ++m) { avg_particle += (1/float(M))*P[m]; }

	Matrix<U_DIM >uinit = (avg_particle - x0) / (DT*(T-1));
	std::vector<Matrix<U_DIM> > U(T-1, uinit);

	std::vector<Matrix<X_DIM> > X(T);
	X[0] = x0;
	for(int t=0; t < T-1; ++t) {
		X[t+1] = point_explore::dynfunc(X[t], U[t]);
	}

	double entropy = point_explore::casadi_differential_entropy(X,U,P);
	std::cout << "entropy: " << entropy << "\n";

	casadi_differential_entropy_func.generateCode("casadi/casadi_de_func.c");
	casadi_grad_differential_entropy_func.generateCode("casadi/casadi_grad_de_func.c");

}
