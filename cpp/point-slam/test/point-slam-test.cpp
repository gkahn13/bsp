#include <vector>
#include <iomanip>

#include "../point-slam.h"

#include "util/matrix.h"
#include "util/Timer.h"
#include "util/logging.h"

std::vector<Matrix<P_DIM> > linearlyInterpolate(Matrix<P_DIM> start, Matrix<P_DIM> end, int steps) {
	std::vector<Matrix<P_DIM> > interp(steps);

	Matrix<P_DIM> u = (end-start)/((double)steps);
	for(int t = 0; t < steps; ++t) {
		interp[t] = start + t*u;
	}

	return interp;
}

void testPlotting(std::vector< Matrix<P_DIM> >& waypoints) {
	std::vector<Matrix<B_DIM> > B(T*NUM_WAYPOINTS, zeros<B_DIM,1>());

	Matrix<P_DIM> start = zeros<P_DIM,1>();
	int index = 0;
	for(int i = 0; i < NUM_WAYPOINTS; ++i) {
		std::vector<Matrix<P_DIM> > interp = linearlyInterpolate(start, waypoints[i], T);
		for(int t = 0; t < T; ++t) {
			B[index++].insert<P_DIM>(0, 0, interp[t]);
		}
		start = waypoints[i];
	}

	pythonDisplayTrajectory(B, waypoints);
}


int main(int argc, char* argv[])
{
	std::vector< Matrix<P_DIM> > waypoints(NUM_WAYPOINTS);
	waypoints[0][0] = 2; waypoints[0][1] = 0;
	waypoints[1][0] = 4; waypoints[1][1] = 0;
	waypoints[2][0] = 4; waypoints[2][1] = 4;
	waypoints[3][0] = 2; waypoints[3][1] = .2;
	waypoints[4][0] = 0; waypoints[4][1] = 4;

	//testPlotting(waypoints);


	std::vector< Matrix<P_DIM> > landmarks(NUM_LANDMARKS);
	landmarks[0][0] = 2; landmarks[0][1] = 0;
	landmarks[1][0] = 4; landmarks[1][1] = 0;
	landmarks[2][0] = 4; landmarks[2][1] = 4;
	landmarks[3][0] = 0; landmarks[3][1] = 4;


	// for now, landmarks == waypoints
	Matrix<X_DIM> x0;
	for(int i = 0; i < L_DIM; i += 2) {
		x0.insert(P_DIM+i, 0, landmarks[i]);
	}

	Matrix<X_DIM,X_DIM> SqrtSigma0 = identity<X_DIM>();
	std::vector<Matrix<B_DIM> > B(T*NUM_LANDMARKS);

	Matrix<P_DIM> pGoal;
	pGoal[0] = 0; pGoal[1] = 0;
	x0.insert(0, 0, pGoal);

	vec(x0, SqrtSigma0, B[0]);

	int index = 1;
	for(int i = 0; i < 1; ++i) {
		Matrix<P_DIM> p0 = pGoal;
		pGoal = waypoints[i];

		x0.insert(0, 0, p0);

		Matrix<U_DIM> uinit = (pGoal - p0) / (T-1);
		//std::vector<Matrix<U_DIM> > U(T-1, uinit);

		for(int t = 0; t < T - 1; ++t) {
			B[index+1] = beliefDynamics(B[index], uinit);
			std::cout << ~B[index+1].subMatrix<P_DIM>(0,0);
			index++;
		}
		//pythonDisplayTrajectory(B, waypoints);
	}

	pythonDisplayTrajectory(B, waypoints);



}

