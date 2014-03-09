#ifndef __SLAM_H__
#define __SLAM_H__

#include <Python.h>

#include <iostream>
#include <fstream>
#include <string>
#include <math.h>

#define _USE_MATH_DEFINES

#include "../util/matrix.h"
#include "../util/logging.h"
#include "../util/utils.h"
#include "../util/Timer.h"

//#include <pythonrun.h>
#include <boost/python.hpp>
//#include <boost/filesystem.hpp>


#define TIMESTEPS 15
#define DT 1.0
#define NUM_LANDMARKS 50
#define NUM_WAYPOINTS 4

#define C_DIM 3 // car dimension [x, y, theta]
#define P_DIM 2 // Position dimension [x,y]
#define L_DIM 2*NUM_LANDMARKS

#define X_DIM (C_DIM+L_DIM)
#define U_DIM 2
#define Z_DIM L_DIM
#define Q_DIM 2
#define R_DIM L_DIM

#define S_DIM (((X_DIM+1)*(X_DIM))/2)
#define B_DIM (X_DIM+S_DIM)

#define XU_DIM (TIMESTEPS*X_DIM + (TIMESTEPS-1)*U_DIM)
#define CU_DIM (TIMESTEPS*C_DIM + (TIMESTEPS-1)*U_DIM)
#define TU_DIM ((TIMESTEPS-1)*U_DIM)

extern double sqrt_time;
extern double not_sqrt_time;


extern const double step;


extern std::vector< Matrix<P_DIM> > waypoints;
extern std::vector< Matrix<P_DIM> > landmarks;

extern Matrix<X_DIM> x0;
extern Matrix<X_DIM,X_DIM> SqrtSigma0;
extern Matrix<X_DIM> xGoal;
extern Matrix<X_DIM> xMin, xMax;
extern Matrix<U_DIM> uMin, uMax;
extern SymmetricMatrix<Q_DIM> Q;
extern SymmetricMatrix<R_DIM> R;

extern const int T;
extern const double INFTY;
extern const std::string landmarks_file;

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

//extern const double alpha_belief, alpha_final_belief, alpha_control, alpha_goal_state;

namespace config {
const double V = 3;
const double MAXG = 30*M_PI/180.;
const double RATEG = 20*M_PI/180.;
const double WHEELBASE = 4;
const double DT_CONTROLS = 0.025;

const double VELOCITY_NOISE = 0.3;
const double TURNING_NOISE = 3.0*M_PI/180.;

const double MAX_RANGE = 3.0; // 5.0
const double DT_OBSERVE = 8*DT_CONTROLS;

const double OBS_DIST_NOISE = 1 * 0.1;
const double OBS_ANGLE_NOISE = 1 * 1.0*M_PI/180.;

const double ALPHA_OBS = .75;
}


void initProblemParams(std::vector<Matrix<P_DIM> >& l);

std::vector<std::vector<Matrix<P_DIM>> > landmarks_list();

Matrix<X_DIM> dynfunc(const Matrix<X_DIM>& x, const Matrix<U_DIM>& u, const Matrix<Q_DIM>& q);

Matrix<C_DIM> dynfunccar(const Matrix<C_DIM>& x, const Matrix<U_DIM>& u);

Matrix<Z_DIM> obsfunc(const Matrix<X_DIM>& x, const Matrix<R_DIM>& r);

Matrix<Z_DIM,Z_DIM> deltaMatrix(const Matrix<X_DIM>& x);


// Jacobians: df(x,u,q)/dx, df(x,u,q)/dq
void linearizeDynamics(const Matrix<X_DIM>& x, const Matrix<U_DIM>& u, const Matrix<Q_DIM>& q, Matrix<3,3>& A, Matrix<3,2>& M);

// Jacobians: df(x,u,q)/dx, df(x,u,q)/du
void linearizeDynamicsFiniteDiff(const Matrix<X_DIM>& x, const Matrix<U_DIM>& u, const Matrix<Q_DIM>& q, Matrix<X_DIM,X_DIM>& A, Matrix<X_DIM,Q_DIM>& M);


// Jacobians: dh(x,r)/dx, dh(x,r)/dr
void linearizeObservation(const Matrix<X_DIM>& x, const Matrix<R_DIM>& r, Matrix<Z_DIM,X_DIM>& H, Matrix<Z_DIM,R_DIM>& N);

// Jacobians: dh(x,r)/dx, dh(x,r)/dr
void linearizeObservationFiniteDiff(const Matrix<X_DIM>& x, const Matrix<R_DIM>& r, Matrix<Z_DIM,X_DIM>& H, Matrix<Z_DIM,R_DIM>& N);

// Switch between belief vector and matrices
void unVec(const Matrix<B_DIM>& b, Matrix<X_DIM>& x, Matrix<X_DIM,X_DIM>& SqrtSigma);

void vec(const Matrix<X_DIM>& x, const Matrix<X_DIM,X_DIM>& SqrtSigma, Matrix<B_DIM>& b);

// Belief dynamics
Matrix<B_DIM> beliefDynamics(const Matrix<B_DIM>& b, const Matrix<U_DIM>& u);

Matrix<B_DIM> casadiBeliefDynamics(const Matrix<B_DIM>& b, const Matrix<U_DIM>& u);

Matrix<B_DIM> beliefDynamicsNoDelta(const Matrix<B_DIM>& b, const Matrix<U_DIM>& u);

void executeControlStep(const Matrix<X_DIM>& x_t_real, const Matrix<B_DIM>& b_t_t, const Matrix<U_DIM>& u_t, Matrix<X_DIM>& x_tp1_real, Matrix<B_DIM>& b_tp1_tp1);

// Jacobians: dg(b,u)/db, dg(b,u)/du
void linearizeBeliefDynamics(const Matrix<B_DIM>& b, const Matrix<U_DIM>& u, Matrix<B_DIM,B_DIM>& F, Matrix<B_DIM,U_DIM>& G, Matrix<B_DIM>& h);

void logDataHandle(std::string file_name, std::ofstream& f);

void logDataToFile(std::ofstream& f, const std::vector<Matrix<B_DIM> >& B, double solve_time, double initialization_time);

void pythonDisplayTrajectory(std::vector< Matrix<X_DIM> >& X, int time_steps, bool pause=false);

void pythonDisplayTrajectory(std::vector< Matrix<U_DIM> >& U, int time_steps, bool pause=false);

void pythonDisplayTrajectory(std::vector< Matrix<B_DIM> >& B, std::vector< Matrix<U_DIM> >& U, std::vector< Matrix<P_DIM> >& waypoints, std::vector< Matrix<P_DIM> >& landmarks, int time_steps, bool pause=false);



#endif
