#include <vector>
#include <iomanip>

#include "../slam.h"

#include "util/matrix.h"
#include "util/Timer.h"
#include "util/logging.h"


void test1()
{
	Q = zeros<Q_DIM,Q_DIM>();
	Q(0,0) = 20*config::VELOCITY_NOISE*config::VELOCITY_NOISE;
	Q(1,1) = 1e-2*config::TURNING_NOISE*config::TURNING_NOISE;

	R = zeros<R_DIM, R_DIM>();
	for(int i = 0; i < R_DIM-1; i+=2) {
		R(i,i) = 1*config::OBS_DIST_NOISE*config::OBS_DIST_NOISE;
		R(i+1,i+1) = 1e-5*config::OBS_ANGLE_NOISE*config::OBS_ANGLE_NOISE;
	}

	std::vector< Matrix<P_DIM> > waypoints(NUM_WAYPOINTS);
	waypoints[0][0] = 30; waypoints[0][1] = 0;
	waypoints[1][0] = 60; waypoints[1][1] = 0;
	waypoints[2][0] = 50; waypoints[2][1] = 20;
	//waypoints[3][0] = 30; waypoints[3][1] = 25;
	//waypoints[3][0] = 30; waypoints[3][1] = 0;
	//waypoints[4][0] = 0; waypoints[4][1] = 20;
	//waypoints[5][0] = 30; waypoints[5][1] = 0;

	//testPlotting(waypoints);

	std::vector< Matrix<P_DIM> > landmarks(NUM_LANDMARKS);
	landmarks[0][0] = 0; landmarks[0][1] = 0;
	landmarks[1][0] = 30; landmarks[1][1] = 0;
	landmarks[2][0] = 60; landmarks[2][1] = 0;
	landmarks[3][0] = 50; landmarks[3][1] = 20;
	//landmarks[4][0] = 0; landmarks[4][1] = 20;

	//landmarks[4][0] = 0; landmarks[4][1] = 40;

	Matrix<X_DIM> x0;
	// for now, landmarks == waypoints
	for(int i = 0; i < NUM_LANDMARKS; ++i) {
		x0.insert(C_DIM+2*i, 0, landmarks[i]);
	}

	//This starts out at 0 for car, landmarks are set based on the car's sigma when first seen
	Matrix<X_DIM,X_DIM> SqrtSigma0 = zeros<X_DIM, X_DIM>();//10*identity<X_DIM>();
	for(int i = 0; i < L_DIM; ++i) { SqrtSigma0(C_DIM+i,C_DIM+i) = 1; }
	//std::cout << "SqrtSigma0" << std::endl << SqrtSigma0 << std::endl;

	int num_controls = 55;
	/*
	double heading_control[55] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
					 0, 0, 0, 0, 0, 0, 0, 0, 0.3491, 0.5236, 0.5236, 0.5236,
				0.5236, 0.3833, 0.1028, 0.0258, 0.3749, 0.5236, 0.5236, 0, 0, 0,0,0,
				0, 0, 0, 0, 0.0000, 0, 0, 0, 0, 0, 0, 0,
							   0, 0, 0, 0, 0, 0, -0.0000 };
							   */
	double heading_control[55] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
							 0, 0, 0, 0, 0, 0, 0, 0, 0.3491, 0.5236, 0.5236, 0.5236,
							 0.5236, 0.3833, 0.1028, 0.0258, 0.3749, 0.5236, 0.5236, 0.5236, 0.3964, 0.1068, 0.0269, 0.0067,
							 0.0017, 0.0004, 0.0001, 0.0000, 0.0000, -0.3491, -0.5236, -0.5236, -0.5114, -0.1623, -0.0209, -0.0052,
							 -0.0013, -0.0003, -0.0001, -0.0000, -0.0000, -0.0000, -0.0000 };

	std::vector<Matrix<B_DIM> > B(num_controls+1);


	vec(x0, SqrtSigma0, B[0]);
	int index = 0;

	std::vector<Matrix<U_DIM> > U(num_controls);
	for(int t=0; t < num_controls; ++t) {
		Matrix<U_DIM> u = zeros<U_DIM,1>();
		u[0] = (t < num_controls - 24) ? config::V : config::V;///1.5;
		u[1] = heading_control[t];
		U[t] = u;
	}


	int k;
	Matrix<X_DIM> xtmp;
	Matrix<X_DIM, X_DIM> stmp;
	for(int t = 0; t <= num_controls-1; ++t) {
		B[index+1] = beliefDynamics(B[index], U[index]);
		unVec(B[index+1], xtmp, stmp);
		//std::cout << "x: " << ~xtmp << std::endl;
		std::cout << "car" << std::endl << tr(stmp.subMatrix<C_DIM,C_DIM>(0,0)) << ": " << stmp(0,0) << " " << stmp(1,1) << " " << stmp(2,2) << std::endl;
		std::cout << "landmarks" << std::endl;
		for(int i=C_DIM; i < X_DIM; ++i) { std::cout << stmp(i,i) << " "; }
		std::cout << std::endl << std::endl;
		index++;

		//std::cout << "press enter" << std::endl;
		//std::cin.ignore();
	}
	//std::cout << "index " << index << std::endl;
	//std::cout << "----------" << std::endl;
	pythonDisplayTrajectory(B, waypoints, num_controls);



	/*
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
	*/

}

void testRectangle()
{
	Q = zeros<Q_DIM,Q_DIM>();
	Q(0,0) = 20*config::VELOCITY_NOISE*config::VELOCITY_NOISE;
	Q(1,1) = 1e-2*config::TURNING_NOISE*config::TURNING_NOISE;

	R = zeros<R_DIM, R_DIM>();
	for(int i = 0; i < R_DIM-1; i+=2) {
		R(i,i) = 1*config::OBS_DIST_NOISE*config::OBS_DIST_NOISE;
		R(i+1,i+1) = 1e-5*config::OBS_ANGLE_NOISE*config::OBS_ANGLE_NOISE;
	}

	//std::vector< Matrix<P_DIM> > waypoints(NUM_WAYPOINTS);
	waypoints[0][0] = 60; waypoints[0][1] = 0;
	waypoints[1][0] = 60; waypoints[1][1] = 20;
	waypoints[2][0] = 0; waypoints[2][1] = 20;

	//std::vector< Matrix<P_DIM> > landmarks(NUM_LANDMARKS);
	landmarks[0][0] = 0; landmarks[0][1] = 0;
	landmarks[1][0] = 30; landmarks[1][1] = 0;
	landmarks[2][0] = 60; landmarks[2][1] = 0;
	landmarks[3][0] = 60; landmarks[3][1] = 20;
	landmarks[4][0] = 30; landmarks[4][1] = 5;

	Matrix<X_DIM> x0;
	for(int i = 0; i < NUM_LANDMARKS; ++i) {
		x0.insert(C_DIM+2*i, 0, landmarks[i]);
	}

	//This starts out at 0 for car, landmarks are set based on the car's sigma when first seen
	Matrix<X_DIM,X_DIM> SqrtSigma0 = zeros<X_DIM, X_DIM>();//10*identity<X_DIM>();
	for(int i = 0; i < L_DIM; ++i) { SqrtSigma0(C_DIM+i,C_DIM+i) = 1; }

	int num_controls = 47;

	/*
	double heading_control[47] = { 0,         0,         0,         0,         0,         0,         0,
	         0,         0,         0,         0,         0,         0,         0,
	         0,         0,         0,         0,         0,    0.3491,    0.5236,
	    0.5236,    0.5236,    0.4765,    0.1325,    0.0334,    0.3825,    0.5236,
	    0.5236,    0.2925,    0.0762,    0.0191,    0.0048,    0.0012,    0.0003,
	    0.0001,    0.0000,    0.0000,    0.0000,    0.0000,    0.0000,    0.0000,
	    0.0000,    0.0000,    0.0000,    0.0000,    0.0000 };
	 */

	double heading_control[47] = { 0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,
	         0,         0,         0,         0,         0,         0,         0,    0.3491,    0.5236,    0.5236,    0.5236,    0.4765,
	    0.1325,    0.0334,    0.3825,    0.5236,    0.5236,    0.5236,    0.3765,    0.1007,    0.0253,    0.0063,    0.0016,    0.0004,
	    0.0001,   -0.3490,   -0.5236,   -0.1771,   -0.0450,   -0.0113,   -0.0028,   -0.0007,   -0.0002,   -0.0000,   -0.0000 };


	std::vector<Matrix<B_DIM> > B(num_controls+1);


	vec(x0, SqrtSigma0, B[0]);
	int index = 0;

	std::vector<Matrix<U_DIM> > U(num_controls);
	for(int t=0; t < num_controls; ++t) {
		Matrix<U_DIM> u = zeros<U_DIM,1>();
		u[0] = config::V;
		u[1] = heading_control[t];
		U[t] = u;
	}


	int k;
	Matrix<X_DIM> xtmp;
	Matrix<X_DIM, X_DIM> stmp;
	for(int t = 0; t <= num_controls-1; ++t) {
		std::cout << "t: " << t << std::endl;
		B[index+1] = beliefDynamics(B[index], U[index]);
		unVec(B[index+1], xtmp, stmp);
		//std::cout << "x: " << ~xtmp << std::endl;
		std::cout << "car" << std::endl << tr(stmp.subMatrix<C_DIM,C_DIM>(0,0)) << ": " << stmp(0,0) << " " << stmp(1,1) << " " << stmp(2,2) << std::endl;
		std::cout << "landmarks" << std::endl;
		for(int i=C_DIM; i < X_DIM; ++i) { std::cout << stmp(i,i) << " "; }
		std::cout << std::endl << std::endl;
		index++;

		//std::cout << "press enter" << std::endl;
		//std::cin.ignore();
	}
	//std::cout << "index " << index << std::endl;
	//std::cout << "----------" << std::endl;
	pythonDisplayTrajectory(B, waypoints, num_controls);
}


int main(int argc, char* argv[])
{
	testRectangle();
}

