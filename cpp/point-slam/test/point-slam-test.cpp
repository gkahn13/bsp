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

	pythonDisplayTrajectory(B, waypoints, T*NUM_WAYPOINTS);
}


int main(int argc, char* argv[])
{
	std::vector< Matrix<P_DIM> > waypoints(NUM_WAYPOINTS);
	waypoints[0][0] = 30; waypoints[0][1] = 0;
	waypoints[1][0] = 60; waypoints[1][1] = 0;
	waypoints[2][0] = 60; waypoints[2][1] = 20;
	//waypoints[3][0] = 30; waypoints[3][1] = 25;
	waypoints[3][0] = 30; waypoints[3][1] = 0;
	waypoints[4][0] = 0; waypoints[4][1] = 20;
	//waypoints[5][0] = 30; waypoints[5][1] = 0;

	//testPlotting(waypoints);


	std::vector< Matrix<P_DIM> > landmarks(NUM_LANDMARKS);
	landmarks[0][0] = 30; landmarks[0][1] = 0;
	landmarks[1][0] = 60; landmarks[1][1] = 0;
	landmarks[2][0] = 60; landmarks[2][1] = 20;
	//landmarks[3][0] = 30; landmarks[3][1] = 20;

	//landmarks[4][0] = 0; landmarks[4][1] = 40;


	// for now, landmarks == waypoints
	Matrix<X_DIM> x0;
	for(int i = 0; i < NUM_LANDMARKS; ++i) {
		x0.insert(P_DIM+2*i, 0, landmarks[i]);
	}

	Matrix<X_DIM,X_DIM> SqrtSigma0 = 10*identity<X_DIM>();
	for(int i = 0; i < L_DIM; ++i) { SqrtSigma0(P_DIM+i,P_DIM+i) = 100; }
	std::vector<Matrix<B_DIM> > B(T*NUM_WAYPOINTS);

	Matrix<P_DIM> pGoal;
	pGoal[0] = 0; pGoal[1] = 0;
	x0.insert(0, 0, pGoal);

	vec(x0, SqrtSigma0, B[0]);
	int index = 0;
	for(int i = 0; i < NUM_WAYPOINTS; ++i) {
		Matrix<P_DIM> p0 = pGoal;
		pGoal = waypoints[i];

		x0.insert(0, 0, p0);

		Matrix<U_DIM> uinit = (pGoal - p0) / (T-1);
		//std::vector<Matrix<U_DIM> > U(T-1, uinit);

		
		Matrix<X_DIM> xtmp;
		Matrix<X_DIM, X_DIM> stmp;
		std::cout << " " << i << std::endl;
		for(int t = 0; t < T - 1; ++t) {
			B[index+1] = beliefDynamics(B[index], uinit);
			unVec(B[index+1], xtmp, stmp);
			//std::cout << tr(stmp) << ": "<<xtmp[0] << ", "<<xtmp[1] << std::endl;
			index++;
		}
		//std::cout << "index " << index << std::endl;
		//std::cout << "----------" << std::endl;
		//pythonDisplayTrajectory(B, waypoints, (i+1)*T);

		// special last case
		if (i < NUM_WAYPOINTS - 1) {
			B[index+1] = B[index];
			index++;
		}
	}

	Matrix<X_DIM> xtmp;
	Matrix<X_DIM, X_DIM> stmp;

	for(int t = 0; t < T*NUM_WAYPOINTS; ++t) {
		unVec(B[t], xtmp, stmp);
		//std::cout << stmp.subMatrix<P_DIM,P_DIM>(0,0) << std::endl;
		//std::cout << stmp << std::endl;
		for(int i = 0; i < X_DIM; ++i) {
			std::cout << stmp(i,i) << " ";
		}
		std::cout << std::endl << std::endl;
	}


	//unVec(B[14], xtmp, stmp);
	//obsfunc(xtmp,zeros<R_DIM,1>());

	pythonDisplayTrajectory(B, waypoints, T*NUM_WAYPOINTS);



}

