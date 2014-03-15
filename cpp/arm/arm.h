#ifndef __ARM_H__
#define __ARM_H__

#include <fstream>

#include "util/matrix.h"

#include "util/logging.h"

#include <iostream>

#define TIMESTEPS 15
#define DT 1
#define X_DIM 6
#define U_DIM 6
#define Z_DIM 4
#define Q_DIM 6
#define R_DIM 4
#define G_DIM 3 // dimension of g(x) function, which is 3d position
#define XU_DIM (X_DIM*T+U_DIM*(T-1))

#define S_DIM (((X_DIM+1)*X_DIM)/2)
#define B_DIM (X_DIM+S_DIM)
#define XU_DIM (X_DIM*T+U_DIM*(T-1))
#define OPT_DIM (X_DIM*T+U_DIM*(T-1)+2*G_DIM)

const double l4 = 2.375;
const double l3 = 10.375;
const double l2 = 8;
const double l1 = 7.25;

const double step = 0.0078125*0.0078125;

Matrix<X_DIM> x0;
Matrix<X_DIM,X_DIM> SqrtSigma0;
Matrix<U_DIM, U_DIM> QC; // process noise covariance
Matrix<Z_DIM, Z_DIM> RC; // measurement noise covariance
	
Matrix<G_DIM> posGoal;
Matrix<X_DIM> xGoal; // TODO: temporary, since goal is a vector of joints
Matrix<X_DIM> xMin, xMax;
Matrix<U_DIM> uMin, uMax;

Matrix<G_DIM> cam0, cam1;

const double goaldelta = 0.1;

const int T = TIMESTEPS;
const double INFTY = 1e10;

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

// iLQG params
SymmetricMatrix<U_DIM> Rint;
SymmetricMatrix<X_DIM> Qint;
SymmetricMatrix<X_DIM> QGoal;
SymmetricMatrix<X_DIM> Sigma0;

const double alpha_belief = 10, alpha_final_belief = 10, alpha_control = 1, alpha_goal_state = 150;

// callisto visualization stuff
int cal_environment, cal_paths, cal_objects, cal_ellipse;
int joint_group[6];
float null_joint[6];
std::ifstream ifs; 

Matrix<X_DIM> dynfunc(const Matrix<X_DIM>& x, const Matrix<U_DIM>& u, const Matrix<U_DIM>& q)
{
	Matrix<X_DIM> xNew = x + (u + q)*DT;
	return xNew;
}

// joint angles -> end effector position
Matrix<G_DIM> g(const Matrix<X_DIM>& x)
{
	Matrix<X_DIM> sx, cx;
	for(int i = 0; i < X_DIM; ++i) {
		sx[i] = sin(x[i]);
		cx[i] = cos(x[i]);
	}

	Matrix<G_DIM> p;
	p[0] = sx[0]*(cx[1]*(sx[2]*(cx[4]*l4+l3)+cx[2]*cx[3]*sx[4]*l4)+sx[1]*(cx[2]*(cx[4]*l4+l3)-sx[2]*cx[3]*sx[4]*l4+l2))+cx[0]*sx[3]*sx[4]*l4;
	p[1] = -sx[1]*(sx[2]*(cx[4]*l4+l3)+cx[2]*cx[3]*sx[4]*l4)+cx[1]*(cx[2]*(cx[4]*l4+l3)-sx[2]*cx[3]*sx[4]*l4+l2)+l1;
	p[2] = cx[0]*(cx[1]*(sx[2]*(cx[4]*l4+l3)+cx[2]*cx[3]*sx[4]*l4)+sx[1]*(cx[2]*(cx[4]*l4+l3)-sx[2]*cx[3]*sx[4]*l4+l2))-sx[0]*sx[3]*sx[4]*l4;

	return p;
}

// Observation model
Matrix<Z_DIM> obsfunc(const Matrix<X_DIM>& x, const Matrix<R_DIM>& r)
{
	Matrix<G_DIM> ee_pos = g(x);

	Matrix<Z_DIM> obs;
	obs[0] = (ee_pos[0] - cam0[0]) / (ee_pos[1] - cam0[1]) + r[0];
	obs[1] = (ee_pos[2] - cam0[2]) / (ee_pos[1] - cam0[1]) + r[1];
	obs[2] = (ee_pos[0] - cam1[0]) / (ee_pos[1] - cam1[1]) + r[2];
    obs[3] = (ee_pos[2] - cam1[2]) / (ee_pos[1] - cam1[1]) + r[3];

    return obs;
}

// Jacobians: df(x,u,q)/dx, df(x,u,q)/dq
void linearizeDynamics(const Matrix<X_DIM>& x, const Matrix<U_DIM>& u, const Matrix<Q_DIM>& q, Matrix<X_DIM,X_DIM>& A, Matrix<X_DIM,Q_DIM>& M)
{
	A.reset();
	Matrix<X_DIM> xr(x), xl(x);
	for (size_t i = 0; i < X_DIM; ++i) {
		xr[i] += step; xl[i] -= step;
		A.insert(0,i, (dynfunc(xr, u, q) - dynfunc(xl, u, q)) / (xr[i] - xl[i]));
		xr[i] = x[i]; xl[i] = x[i];
	}

	M.reset();
	Matrix<Q_DIM> qr(q), ql(q);
	for (size_t i = 0; i < Q_DIM; ++i) {
		qr[i] += step; ql[i] -= step;
		M.insert(0,i, (dynfunc(x, u, qr) - dynfunc(x, u, ql)) / (qr[i] - ql[i]));
		qr[i] = q[i]; ql[i] = q[i];
	}
}

void linearizeg(const Matrix<X_DIM>& x, Matrix<G_DIM, X_DIM>& J)
{
	Matrix<X_DIM> sx, cx;
	for(int i = 0; i < X_DIM; ++i) {
		sx[i] = sin(x[i]);
		cx[i] = cos(x[i]);
	}

    J(0,0) = cx[0]*(cx[1]*(sx[2]*(cx[4]*l4+l3)+cx[2]*cx[3]*sx[4]*l4)+sx[1]*(cx[2]*(cx[4]*l4+l3)-sx[2]*cx[3]*sx[4]*l4+l2))-sx[0]*sx[3]*sx[4]*l4;
    J(0,1) = sx[0]*(cx[1]*(cx[2]*(cx[4]*l4+l3)-sx[2]*cx[3]*sx[4]*l4+l2)-sx[1]*(sx[2]*(cx[4]*l4+l3)+cx[2]*cx[3]*sx[4]*l4));
    J(0,2) = sx[0]*(sx[1]*(-sx[2]*(cx[4]*l4+l3)-cx[2]*cx[3]*sx[4]*l4)+cx[1]*(cx[2]*(cx[4]*l4+l3)-sx[2]*cx[3]*sx[4]*l4));
    J(0,3) = sx[0]*(sx[1]*sx[2]*sx[3]*sx[4]*l4-cx[1]*cx[2]*sx[3]*sx[4]*l4)+cx[0]*cx[3]*sx[4]*l4;
    J(0,4) = sx[0]*(cx[1]*(cx[2]*cx[3]*cx[4]*l4-sx[2]*sx[4]*l4)+sx[1]*(-cx[2]*sx[4]*l4-sx[2]*cx[3]*cx[4]*l4))+cx[0]*sx[3]*cx[4]*l4;
    J(0,5) = 0;

    J(1,0) = 0;
    J(1,1) = -cx[1]*(sx[2]*(cx[4]*l4+l3)+cx[2]*cx[3]*sx[4]*l4)-sx[1]*(cx[2]*(cx[4]*l4+l3)-sx[2]*cx[3]*sx[4]*l4+l2);
    J(1,2) = cx[1]*(-sx[2]*(cx[4]*l4+l3)-cx[2]*cx[3]*sx[4]*l4)-sx[1]*(cx[2]*(cx[4]*l4+l3)-sx[2]*cx[3]*sx[4]*l4);
    J(1,3) = cx[1]*sx[2]*sx[3]*sx[4]*l4+sx[1]*cx[2]*sx[3]*sx[4]*l4;
    J(1,4) = cx[1]*(-cx[2]*sx[4]*l4-sx[2]*cx[3]*cx[4]*l4)-sx[1]*(cx[2]*cx[3]*cx[4]*l4-sx[2]*sx[4]*l4);
    J(1,5) = 0;

    J(2,0) = -sx[0]*(cx[1]*(sx[2]*(cx[4]*l4+l3)+cx[2]*cx[3]*sx[4]*l4)+sx[1]*(cx[2]*(cx[4]*l4+l3)-sx[2]*cx[3]*sx[4]*l4+l2))-cx[0]*sx[3]*sx[4]*l4;
    J(2,1) = cx[0]*(cx[1]*(cx[2]*(cx[4]*l4+l3)-sx[2]*cx[3]*sx[4]*l4+l2)-sx[1]*(sx[2]*(cx[4]*l4+l3)+cx[2]*cx[3]*sx[4]*l4));
    J(2,2) = cx[0]*(sx[1]*(-sx[2]*(cx[4]*l4+l3)-cx[2]*cx[3]*sx[4]*l4)+cx[1]*(cx[2]*(cx[4]*l4+l3)-sx[2]*cx[3]*sx[4]*l4));
    J(2,3) = cx[0]*(sx[1]*sx[2]*sx[3]*sx[4]*l4-cx[1]*cx[2]*sx[3]*sx[4]*l4)-sx[0]*cx[3]*sx[4]*l4;
    J(2,4) = cx[0]*(cx[1]*(cx[2]*cx[3]*cx[4]*l4-sx[2]*sx[4]*l4)+sx[1]*(-cx[2]*sx[4]*l4-sx[2]*cx[3]*cx[4]*l4))-sx[0]*sx[3]*cx[4]*l4;
    J(2,5) = 0;
}

// Jacobians: dh(x,r)/dx, dh(x,r)/dr
void linearizeObservation(const Matrix<X_DIM>& x, const Matrix<R_DIM>& r, Matrix<Z_DIM,X_DIM>& H, Matrix<Z_DIM,R_DIM>& N)
{
	H.reset();
	Matrix<X_DIM> xr(x), xl(x);
	for (size_t i = 0; i < X_DIM; ++i) {
		xr[i] += step; xl[i] -= step;
		H.insert(0,i, (obsfunc(xr, r) - obsfunc(xl, r)) / (xr[i] - xl[i]));
		xr[i] = x[i]; xl[i] = x[i];
	}

	N.reset();
	Matrix<R_DIM> rr(r), rl(r);
	for (size_t i = 0; i < R_DIM; ++i) {
		rr[i] += step; rl[i] -= step;
		N.insert(0,i, (obsfunc(x, rr) - obsfunc(x, rl)) / (rr[i] - rl[i]));
		rr[i] = r[i]; rl[i] = r[i];
	}
}

// Switch between belief vector and matrices
void unVec(const Matrix<B_DIM>& b, Matrix<X_DIM>& x, Matrix<X_DIM,X_DIM>& SqrtSigma) {
	x = b.subMatrix<X_DIM,1>(0,0);
	size_t idx = X_DIM;
	for (size_t j = 0; j < X_DIM; ++j) {
		for (size_t i = j; i < X_DIM; ++i) {
			SqrtSigma(i,j) = b[idx];
			SqrtSigma(j,i) = b[idx];
			++idx;
		}
	}
}

void vec(const Matrix<X_DIM>& x, const Matrix<X_DIM,X_DIM>& SqrtSigma, Matrix<B_DIM>& b) {
	b.insert(0,0,x);
	size_t idx = X_DIM;
	for (size_t j = 0; j < X_DIM; ++j) {
		for (size_t i = j; i < X_DIM; ++i) {
			b[idx] = 0.5 * (SqrtSigma(i,j) + SqrtSigma(j,i));
			++idx;
		}
	}
}


// Belief dynamics
Matrix<B_DIM> beliefDynamics(const Matrix<B_DIM>& b, const Matrix<U_DIM>& u) {
	Matrix<X_DIM> x;
	Matrix<X_DIM,X_DIM> SqrtSigma;
	unVec(b, x, SqrtSigma);

	Matrix<X_DIM,X_DIM> Sigma = SqrtSigma*SqrtSigma;

	Matrix<X_DIM,X_DIM> A = identity<X_DIM>();
	Matrix<X_DIM,Q_DIM> M = DT*identity<U_DIM>();
	//linearizeDynamics(x, u, zeros<Q_DIM,1>(), A, M);

	x = dynfunc(x, u, zeros<Q_DIM,1>());
	Sigma = A*Sigma*~A + M*QC*~M;

	Matrix<Z_DIM,X_DIM> H = zeros<Z_DIM,X_DIM>();
	Matrix<Z_DIM,R_DIM> N = identity<Z_DIM>();
	//linearizeObservation(x, zeros<R_DIM>(), H, N);

	Matrix<G_DIM,X_DIM> J; 
	linearizeg(x, J);

	for (int i = 0; i < X_DIM; ++i) {
		/*
	        H(0,i) = (J(0,i) * (x[2] - cam0[2]) - (x[0] - cam0[0]) * J(2,i)) / ((x[2] - cam0[2]) * (x[2] - cam0[2]));
	    	H(1,i) = (J(1,i) * (x[2] - cam0[2]) - (x[1] - cam0[1]) * J(2,i)) / ((x[2] - cam0[2]) * (x[2] - cam0[2]));
	    	H(2,i) = (J(0,i) * (x[2] - cam1[2]) - (x[0] - cam1[0]) * J(2,i)) / ((x[2] - cam1[2]) * (x[2] - cam1[2]));
	    	H(3,i) = (J(1,i) * (x[2] - cam1[2]) - (x[1] - cam1[1]) * J(2,i)) / ((x[2] - cam1[2]) * (x[2] - cam1[2]));
		 */

		H(0,i) = -(J(0,i) * (x[1] - cam0[1]) - (x[0] - cam0[0]) * J(1,i)) / ((x[1] - cam0[1]) * (x[1] - cam0[1]));
		H(1,i) = -(J(2,i) * (x[1] - cam0[1]) - (x[2] - cam0[2]) * J(1,i)) / ((x[1] - cam0[1]) * (x[1] - cam0[1]));
		H(2,i) = -(J(0,i) * (x[1] - cam1[1]) - (x[0] - cam1[0]) * J(1,i)) / ((x[1] - cam1[1]) * (x[1] - cam1[1]));
		H(3,i) = -(J(2,i) * (x[1] - cam1[1]) - (x[2] - cam1[2]) * J(1,i)) / ((x[1] - cam1[1]) * (x[1] - cam1[1]));
	}


	Matrix<X_DIM,Z_DIM> K = Sigma*~H/(H*Sigma*~H + N*RC*~N);

	Sigma = (identity<X_DIM>() - K*H)*Sigma;

	Matrix<B_DIM> g;
	vec(x, sqrtm(Sigma), g);

	return g;
}

void initProblemParams(int iter)
{
	/*cam0[0] = -4;  cam0[1] = 11.5; cam0[2] = -23;
	cam1[0] = 4;  cam1[1] = 11.5; cam1[2] = -23;*/

	if(iter==0){
		int temp = 0; 
		
	}
	if(iter == 99){
		ifs.close();
	}

	cam0[0] = -4;  cam0[1] = 30; cam0[2] = 0;
	cam1[0] = 4;  cam1[1] = 30; cam1[2] = 0;
	for(int i = 0; i<6; i++){
		double temp =0; 
		ifs>>temp;
		
		x0[i] = temp; 
	}
	SqrtSigma0 = identity<X_DIM>();

	xGoal[0] = -1.4846950311433709; xGoal[1] = -2.2314918647565389; xGoal[2] = 1.4680882089972564;
	xGoal[3] = 0.37654505159140872; xGoal[4] = 1.660179950900027; xGoal[5] = -2.718448168983489;

	posGoal[0] = 11.5; posGoal[1] = 11.5; posGoal[2] = 0;
	//posGoal = g(xGoal);

	//xMin[0] = -M_PI; xMin[1] = -0.75*M_PI; xMin[2] = -0.75*M_PI; xMin[3] = -M_PI; xMin[4] = -0.75*M_PI; xMin[5] = -M_PI;
	//xMax[0] = M_PI;  xMax[1] = 0.75*M_PI;  xMax[2] = 0.75*M_PI;  xMax[3] = M_PI;  xMax[4] =  0.75*M_PI; xMax[5] = M_PI;

	xMin[0] = -M_PI; xMin[1] = -0.75*M_PI; xMin[2] = -0.75*M_PI; xMin[3] = -M_PI; xMin[4] = -0.75*M_PI; xMin[5] = -M_PI;
	xMax[0] = M_PI;  xMax[1] = 0.75*M_PI;  xMax[2] = 0.75*M_PI;  xMax[3] = M_PI;  xMax[4] =  0.75*M_PI; xMax[5] = M_PI;

	//uMin[0] = -0.5; uMin[1] = -0.5; uMin[2] = -0.5; uMin[3] = -0.5; uMin[4] = -0.5; uMin[5] = -0.5;
	//uMax[0] = 0.5; uMax[1] = 0.5; uMax[2] = 0.5; uMax[3] = 0.5; uMax[4] = 0.5; uMax[5] = 0.5;
	
	for(int i=0; i < U_DIM; ++i) { uMin[i] = -.5; uMax[i] = .5; }

	SqrtSigma0 = identity<X_DIM>() * 0.1;

	QC = identity<U_DIM>() * 0.01;
	RC = identity<Z_DIM>() * 0.01;
}

Matrix<4,1> quatFromAA(const Matrix<3>& axis, double angle);
void drawPolyline3d(int pathlen, float* pts, int groupid);
Matrix<4,1> quatFromRot(const Matrix<3,3>& R);
void preprocess(const std::vector<Matrix<X_DIM> >& path, std::vector<Matrix<G_DIM, X_DIM> >& J, std::vector<Matrix<X_DIM,X_DIM> >& Sigma);


// Jacobians: dg(b,u)/db, dg(b,u)/du
void linearizeBeliefDynamics(const Matrix<B_DIM>& b, const Matrix<U_DIM>& u, Matrix<B_DIM,B_DIM>& F, Matrix<B_DIM,U_DIM>& G, Matrix<B_DIM>& h)
{
	F.reset();
	Matrix<B_DIM> br(b), bl(b);
	for (size_t i = 0; i < B_DIM; ++i) {
		br[i] += step; bl[i] -= step;
		F.insert(0,i, (beliefDynamics(br, u) - beliefDynamics(bl, u)) / (br[i] - bl[i]));
		br[i] = b[i]; bl[i] = b[i];
	}

	G.reset();
	Matrix<U_DIM> ur(u), ul(u);
	for (size_t i = 0; i < U_DIM; ++i) {
		ur[i] += step; ul[i] -= step;
		G.insert(0,i, (beliefDynamics(b, ur) - beliefDynamics(b, ul)) / (ur[i] - ul[i]));
		ur[i] = u[i]; ul[i] = u[i];
	}

	h = beliefDynamics(b, u);
}





void preprocess(const std::vector<Matrix<X_DIM> >& path, std::vector<Matrix<G_DIM, X_DIM> >& J, std::vector<Matrix<X_DIM,X_DIM> >& Sigma) 
{
	std::vector<Matrix<X_DIM, X_DIM> > A; // process matrix
	std::vector<Matrix<X_DIM, U_DIM> > B; // input matrix
	Matrix<X_DIM, X_DIM> C; // LQR state deviation cost
	Matrix<U_DIM, U_DIM> D; // LQR input deviation cost
	std::vector<Matrix<X_DIM, U_DIM> > V; // process noise matrix
	Matrix<U_DIM, U_DIM> M; // process noise covariance
	std::vector<Matrix<Z_DIM, X_DIM> > H; // measurement matrix
	std::vector<Matrix<Z_DIM, Z_DIM> > W; // measurement noise matrix
	Matrix<Z_DIM, Z_DIM> N; // measurement noise covariance
	Matrix<X_DIM, X_DIM> P_0; // initial state covariance
	std::vector<Matrix<U_DIM, X_DIM> > L; // feedback matrix
	std::vector<Matrix<X_DIM, Z_DIM> > K; // Kalman-gain matrix

	std::vector<Matrix<2*X_DIM, 2*X_DIM> > R;

	int l = (int) path.size();

	A.resize(l); // process matrix
	B.resize(l); // input matrix
	V.resize(l); // process noise matrix
	H.resize(l); // measurement matrix
	W.resize(l); // measurement noise matrix
	L.resize(l); // feedback matrix
	K.resize(l); // Kalman-gain matrix
	R.resize(l);

	J.resize(l);
	Sigma.resize(l);

	// Initialization
	C = identity<X_DIM>();
	D = identity<U_DIM>();

	P_0 = identity<X_DIM>() * 0.01;

	M = identity<U_DIM>() * 0.01;
	N = identity<Z_DIM>() * 0.01;

	// Jacobians
	for (int k = 0; k < l; ++k) 
	{
		Matrix<G_DIM> g_path = g(path[k]);

		A[k] = identity<X_DIM>();
		B[k] = DT * identity<U_DIM>();
		V[k] = B[k];

		linearizeg(path[k], J[k]);

		for (int i = 0; i < X_DIM; ++i) {
			/*H[k](0,i) = (J[k](0,i)*(g_path[2] - cam0[2]) - (g_path[0] - cam0[0])*J[k](2,i)) / ((g_path[2] - cam0[2])*(g_path[2] - cam0[2]));
			H[k](1,i) = (J[k](1,i)*(g_path[2] - cam0[2]) - (g_path[1] - cam0[1])*J[k](2,i)) / ((g_path[2] - cam0[2])*(g_path[2] - cam0[2]));
			H[k](2,i) = (J[k](0,i)*(g_path[2] - cam1[2]) - (g_path[0] - cam1[0])*J[k](2,i)) / ((g_path[2] - cam1[2])*(g_path[2] - cam1[2]));
			H[k](3,i) = (J[k](1,i)*(g_path[2] - cam1[2]) - (g_path[1] - cam1[1])*J[k](2,i)) / ((g_path[2] - cam1[2])*(g_path[2] - cam1[2]));*/

			H[k](0,i) = -(J[k](0,i)*(g_path[1] - cam0[1]) - (g_path[0] - cam0[0])*J[k](1,i)) / ((g_path[1] - cam0[1])*(g_path[1] - cam0[1]));
			H[k](1,i) = -(J[k](2,i)*(g_path[1] - cam0[1]) - (g_path[2] - cam0[2])*J[k](1,i)) / ((g_path[1] - cam0[1])*(g_path[1] - cam0[1]));
			H[k](2,i) = -(J[k](0,i)*(g_path[1] - cam1[1]) - (g_path[0] - cam1[0])*J[k](1,i)) / ((g_path[1] - cam1[1])*(g_path[1] - cam1[1]));
			H[k](3,i) = -(J[k](2,i)*(g_path[1] - cam1[1]) - (g_path[2] - cam1[2])*J[k](1,i)) / ((g_path[1] - cam1[1])*(g_path[1] - cam1[1]));
		}

		W[k] = identity<Z_DIM>();
	}

	// LQR
	Matrix<X_DIM, X_DIM> S;
	S = C;
	L[l - 1] = zeros<U_DIM, X_DIM>();
	for (int k = l - 2; k >= 0; --k) {
		L[k] = -!(~B[k+1]*S*B[k+1] + D)*~B[k+1]*S*A[k+1];
		S = C + ~A[k+1]*S*A[k+1] + ~A[k+1]*S*B[k+1]*L[k];
	}

	// Kalman
	Matrix<X_DIM, X_DIM> P;
	P = P_0;
	K[0] = zeros<X_DIM, Z_DIM>();
	for (int k = 1; k < l; ++k) {
		P = A[k]*P*~A[k] + V[k]*M*~V[k];
		K[k] = P*~H[k]*!(H[k]*P*~H[k] + W[k]*N*~W[k]);
		P = (identity<X_DIM>() - K[k]*H[k])*P;
	}

	std::vector<Matrix<2*X_DIM, 2*X_DIM> > F(l);
	std::vector<Matrix<2*X_DIM, U_DIM+Z_DIM> > G(l);

	// Combination of LQR and Kalman
	Matrix<U_DIM+Z_DIM, U_DIM+Z_DIM> Q = zeros<U_DIM+Z_DIM, U_DIM+Z_DIM>();
	Q.insert(0,0, M); Q.insert(U_DIM, U_DIM, N);

	R[0].reset();
	R[0].insert(0,0,P_0);
	Sigma[0] = P_0; //R[0].subMatrix<X_DIM,X_DIM>(0,0);

	for (int k = 1; k < l; ++k) {
		F[k].insert(0,0,     A[k]);           F[k].insert(0,X_DIM,     B[k]*L[k-1]);
		F[k].insert(X_DIM,0, K[k]*H[k]*A[k]); F[k].insert(X_DIM,X_DIM, A[k] + B[k]*L[k-1] - K[k]*H[k]*A[k]);

		G[k].insert(0,0,     V[k]);           G[k].insert(0,U_DIM,     zeros<X_DIM,Z_DIM>());
		G[k].insert(X_DIM,0, K[k]*H[k]*V[k]); G[k].insert(X_DIM,U_DIM, K[k]*W[k]);

		R[k] = F[k]*R[k-1]*~F[k] + G[k]*Q*~G[k];
		Sigma[k] = R[k].subMatrix<X_DIM,X_DIM>(0,0);
	}
}


// helper functions

inline Matrix<4,1> quatFromAA(const Matrix<3>& axis, double angle)
{
	Matrix<4,1> q;
	double sa = sin(-angle*0.5);
	double ca = cos(angle*0.5);
	q(0,0) = sa*axis[0];
	q(1,0) = sa*axis[1];
	q(2,0) = sa*axis[2];
	q(3,0) = ca;

	return q;
}

Matrix<4,1> quatFromRot(const Matrix<3,3>& R) {
	double x = R(2,1) - R(1,2);
	double y = R(0,2) - R(2,0);
	double z = R(1,0) - R(0,1);
	double r = sqrt(x*x+y*y+z*z);
	double t = R(0,0) + R(1,1) + R(2,2);
	double angle = atan2(r,t-1);
	if (angle != 0) {
		x /= r;
		y /= r;
		z /= r;
	} else {
		x = 0;
		y = 0;
		z = 0;
	}
	Matrix<4,1> q;
	q(0,0) = sin(angle/2)*x;
	q(1,0) = sin(angle/2)*y;
	q(2,0) = sin(angle/2)*z;
	q(3,0) = cos(angle/2);

	return q;
}


inline void saveOptimizedTrajectory(std::vector<Matrix<U_DIM> >& Uopt) {
	std::ofstream fptr("optimized-controls.txt", std::ios::out);
	fptr << T << std::endl;
	for(int t = 0; t < T-1; ++t) {
		fptr << ~Uopt[t];
	}
	fptr.close();
}

inline void readOptimizedTrajectory(std::vector<Matrix<U_DIM> >& U) {
	std::ifstream fptr("optimized-controls.txt", std::ios::in);
	int nsteps;
	fptr >> nsteps;
	//std::cout << "nsteps: " << nsteps << std::endl;
	U.clear();
	U.resize(nsteps-1);
	for(int t = 0; t < nsteps-1; ++t) {
		fptr >> U[t];
		//std::cout << ~U[t];
	}
	fptr.close();
}

void readTrajectoryFromFile(const std::string& filename, std::vector<Matrix<X_DIM> >& path)
{
	std::ifstream fptr(filename.c_str(),std::ios::in);
	int nsteps;
	fptr >> nsteps;
	path.resize(nsteps);
	for(int i = 0; i < nsteps; ++i) {
		fptr >> path[i];
	}
	fptr.close();
}

// utility to fill Matrix in column major format in FORCES array
template <size_t _numRows>
inline void fillCol(double *X, const Matrix<_numRows>& XCol) {
	int idx = 0;
	for(int r = 0; r < _numRows; ++r) {
		X[idx++] = XCol[r];
	}
}

template <size_t _numRows, size_t _numColumns>
inline void fillColMajor(double *X, const Matrix<_numRows, _numColumns>& XMat) {
	int idx = 0;
	for(int c = 0; c < _numColumns; ++c) {
		for(int r = 0; r < _numRows; ++r) {
			X[idx++] = XMat[c + r*_numColumns];
		}
	}
}

#endif
