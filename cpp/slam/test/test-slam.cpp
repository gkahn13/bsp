#include <vector>
#include <iomanip>

#include "../slam.h"

#include "util/matrix.h"
#include "util/Timer.h"
#include "util/logging.h"



int main(int argc, char* argv[])
{
	Q = zeros<Q_DIM,Q_DIM>();
	Q(0,0) = config::V*config::V;
	Q(1,1) = config::TURNING_NOISE*config::TURNING_NOISE;

	R = zeros<R_DIM, R_DIM>();
	for(int i = 0; i < R_DIM-1; i+=2) {
		R(i,i) = config::OBS_DIST_NOISE*config::OBS_DIST_NOISE;
		R(i+1,i+1) = config::OBS_ANGLE_NOISE*config::OBS_ANGLE_NOISE;
	}

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
		x0.insert(C_DIM+2*i, 0, landmarks[i]);
	}

	//This starts out at 0 for car, landmarks are set based on the car's sigma when first seen
	Matrix<X_DIM,X_DIM> SqrtSigma0 = zeros<X_DIM, X_DIM>();//10*identity<X_DIM>();
	for(int i = 0; i < L_DIM; ++i) { SqrtSigma0(C_DIM+i,C_DIM+i) = 100; }
	std::vector<Matrix<B_DIM> > B(T*NUM_WAYPOINTS);

	Matrix<P_DIM> pGoal;
	pGoal[0] = 0; pGoal[1] = 0;
	x0.insert(0, 0, pGoal);

	vec(x0, SqrtSigma0, B[0]);
	int index = 0;


	/*
Controls used when running openslam on our example set up and DT=1, V=3 (~10 timesteps per waypoint)

 Columns 1 through 12
         0         0         0         0         0         0         0         0         0         0         0         0
  Columns 13 through 24
         0         0         0         0         0         0         0         0    0.3491    0.5236    0.5236    0.5236
  Columns 25 through 36
    0.5236    0.3833    0.1028    0.0258    0.3749    0.5236    0.5236    0.5236    0.3964    0.1068    0.0269    0.0067
  Columns 37 through 48
    0.0017    0.0004    0.0001    0.0000    0.0000   -0.3491   -0.5236   -0.5236   -0.5114   -0.1623   -0.0209   -0.0052
  Columns 49 through 55
   -0.0013   -0.0003   -0.0001   -0.0000   -0.0000   -0.0000   -0.0000
	*/
	for(int i = 0; i < 2;++i){//NUM_WAYPOINTS; ++i) {
		Matrix<P_DIM> p0 = pGoal;
		pGoal = waypoints[i];

		x0.insert(0, 0, p0);

		


		//std::vector<Matrix<U_DIM> > U(T-1, uinit);


		Matrix<X_DIM> xtmp;
		Matrix<X_DIM, X_DIM> stmp;
		std::cout << " " << i << std::endl;
		for(int t = 0; t < T - 1; ++t) {
			B[index+1] = beliefDynamics(B[index], U[index]);
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

