#include <vector>
#include <iomanip>

#include "../point-slam.h"

#include "util/matrix.h"
#include "util/Timer.h"
#include "util/logging.h"



int main(int argc, char* argv[])
{
	std::vector< Matrix<P_DIM> > waypoints(NUM_LANDMARKS);
	waypoints[0][0] = 2; waypoints[0][1] = 0;
	waypoints[1][2] = 4; waypoints[1][3] = 0;
	waypoints[2][4] = 4; waypoints[2][5] = 4;
	waypoints[3][6] = 0; waypoints[3][7] = 4;

	// for now, landmarks == waypoints
	Matrix<X_DIM> x0;
	for(int i = 0; i < L_DIM; i += 2) {
		x0.insert(P_DIM+i, 0, waypoints[i]);
	}

	Matrix<X_DIM,X_DIM> SqrtSigma0 = identity<X_DIM>();
	std::vector<Matrix<B_DIM> > B(T*NUM_LANDMARKS);

	Matrix<P_DIM> pGoal;
	pGoal[0] = 0; pGoal[1] = 0;
	x0.insert(0, 0, pGoal);

	vec(x0, SqrtSigma0, B[0]);

	int index = 1;
	for(int i = 0; i < NUM_LANDMARKS; ++i) {
		Matrix<P_DIM> p0 = pGoal;
		pGoal = waypoints[i];

		x0.insert(0, 0, p0);

		Matrix<U_DIM> uinit = (pGoal - p0[0]) / (T-1);
		std::vector<Matrix<U_DIM> > U(T-1, uinit);

		for(int t = 0; t < T - 1; ++t) {
			B[index+1] = beliefDynamics(B[index], U[t]);
			index++;
		}
	}

	Matrix<P_DIM> pInitial;
	pInitial[0] = 0; pInitial[1] = 0;
	pythonDisplayTrajectory(B, NULL, pInitial)

}
