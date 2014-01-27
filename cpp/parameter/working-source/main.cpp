//#define MAX_LIKELIHOOD // Determines whether maximum likelihood measurements should be used. Good for debugging.
#define ITERATE_KALMAN 1 // Determines whether the Kalman filter in the simulation System should smooth and if so how many iterations it should perform.
#define INVERT_PARAMETER

//#include "Dynamics-2dplanararm-masses.hpp"
//#include "Dynamics-2dplanararm-damping.hpp"
//#include "Dynamics-2dplanararm-lengths.hpp"
#include "Dynamics-2dplanararm-masses-lengths.hpp"

int main(int argc, char* argv[])
{
	//unsigned int fp_control_word, new_fp_control_word;
	//_controlfp_s(&fp_control_word, 0, 0);
	////new_fp_control_word = fp_control_word & ~(_EM_INVALID | _EM_DENORMAL | _EM_ZERODIVIDE | _EM_OVERFLOW | _EM_UNDERFLOW /*| _EM_INEXACT*/);
	//new_fp_control_word = fp_control_word & ~(_EM_INVALID | _EM_ZERODIVIDE | _EM_OVERFLOW);
	//_controlfp_s(&fp_control_word, new_fp_control_word, _MCW_EM);

	/*
	std::vector< Matrix<U_DIM> > uBartemp;
	std::vector< Matrix<U_DIM, X_DIM> > Ltemp;
	readControlPolicy("control_policy.bin", uBartemp, Ltemp);
	*/
	Dynamics e;

	srand(0);

	Simulation sys(e);

	//test(e, sys);
	//gatherData(e, sys);
	//analyzePlanner(e, sys, 1);

	runMPC(e, sys);

	int k;
	std::cin >> k;

	return 0;
}