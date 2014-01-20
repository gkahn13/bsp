// PlanningControl.cpp : Defines the entry point for the console application.
//

#define _CRT_RAND_S

#define DIM 3
#define X_DIM 6
#define U_DIM 6
#define Z_DIM 4

#include "targetver.h"

#include <stdio.h>
#include <tchar.h>
#include <vector>
#include <map>
#include <queue>
#include <time.h>
#include <float.h>

#include "matrix.h"

struct RRTNode {
  Matrix<X_DIM> x;
  Matrix<U_DIM> u;
  int bp;
  Matrix<DIM> p;
};

struct PathNode {
  Matrix<X_DIM> x;
  Matrix<U_DIM> u;
};

int robot_nr;

double l4 = 2.375;
double l3 = 10.375;
double l2 = 8;
double l1 = 7.25;


int joint_group[6];
float null_joint[6];
//int cal_samples;
Matrix<U_DIM> u_min, u_max;
Matrix<X_DIM> x_min, x_max;
Matrix<X_DIM> start;
Matrix<DIM> goal;

double goalRadius;
double dt;

std::vector<Matrix<DIM> > cam(2);

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
std::vector<Matrix<DIM, X_DIM> > T;

std::vector<Matrix<2*X_DIM, 2*X_DIM> > R;

Matrix<X_DIM> f(const Matrix<X_DIM>& x, const Matrix<U_DIM>& u, const Matrix<U_DIM>& m) {
  Matrix<X_DIM> a = x + (u + m) * dt;
  return a;
}

Matrix<DIM> g(const Matrix<X_DIM>& x) {
  double a0 = x[0];
  double a1 = x[1];
  double a2 = x[2];
  double a3 = x[3];
  double a4 = x[4];
  double a5 = x[5];
    
  Matrix<DIM> p;
  p[0] = sin(a0)*(cos(a1)*(sin(a2)*(cos(a4)*l4+l3)+cos(a2)*cos(a3)*sin(a4)*l4)+sin(a1)*(cos(a2)*(cos(a4)*l4+l3)-sin(a2)*cos(a3)*sin(a4)*l4+l2))+cos(a0)*sin(a3)*sin(a4)*l4;
  p[1] = -sin(a1)*(sin(a2)*(cos(a4)*l4+l3)+cos(a2)*cos(a3)*sin(a4)*l4)+cos(a1)*(cos(a2)*(cos(a4)*l4+l3)-sin(a2)*cos(a3)*sin(a4)*l4+l2)+l1;
  p[2] = cos(a0)*(cos(a1)*(sin(a2)*(cos(a4)*l4+l3)+cos(a2)*cos(a3)*sin(a4)*l4)+sin(a1)*(cos(a2)*(cos(a4)*l4+l3)-sin(a2)*cos(a3)*sin(a4)*l4+l2))-sin(a0)*sin(a3)*sin(a4)*l4;

  return p;
}

void showState(const Matrix<X_DIM>& x) {
  
  //CAL_SetGroupPosition(cal_obstacles, 1, 0, 1);

  for (int k = 0; k < X_DIM; ++k) {
    double angle = x[k] + null_joint[k];
    if (angle >= 2*M_PI) {
      angle -= 2*M_PI;
    }
    if (angle < 0) {
      angle += 2*M_PI;
    }
    
    CAL_SetGroupOrientation(joint_group[k], 0, (float) angle, 0);
   
  }
  Matrix<DIM> p = g(x);
  //CAL_SetGroupPosition(cal_ellipse, p[0], p[1], p[2]); 

  
}

double random() {
  unsigned int i;
  rand_s(&i);
  return ((double) i) / UINT_MAX;
}

void initEnvironment() {
  x_min[0] = -M_PI; x_min[1] = -0.75*M_PI; x_min[2] = -0.75*M_PI; x_min[3] = -M_PI; x_min[4] = -0.75*M_PI; x_min[5] = -M_PI;
  x_max[0] = M_PI;  x_max[1] = 0.75*M_PI;  x_max[2] = 0.75*M_PI;  x_max[3] = M_PI;  x_max[4] =  0.75*M_PI; x_max[5] = M_PI;

  u_min[0] = -0.5; u_min[1] = -0.5; u_min[2] = -0.5; u_min[3] = -0.5; u_min[4] = -0.5; u_min[5] = -0.5;
  u_max[0] = 0.5; u_max[1] = 0.5; u_max[2] = 0.5; u_max[3] = 0.5; u_max[4] = 0.5; u_max[5] = 0.5;
  
  start[0] = M_PI/2; start[1] = -1.5431281995798991; start[2] = -0.047595544887998331;
  start[3] = 1.4423058659586809; start[4] = 1.5334368368992011; start[5] = -1.1431255223182604;
  
  goal[0] = 11.5; goal[1] = 11.5; goal[2] = 0;
  goalRadius = 1.5;

  dt = 0.5;

  /*cam[0][0] = -4;  cam[0][1] = 11.5; cam[0][2] = -23;
  cam[1][0] = 4;  cam[1][1] = 11.5; cam[1][2] = -23;*/

  cam[0][0] = -4;  cam[0][1] = 23; cam[0][2] = 0;
  cam[1][0] = 4;  cam[1][1] = 23; cam[1][2] = 0;

  //CAL_SetViewParams(0, 0, 0, 20, 0, 0, 0, 0, 1, 0); 
    
  CAL_CreateGroup(&cal_rrt, 0, false, "RRT");
  CAL_SetGroupColor(cal_rrt, 0, 1, 0);

  CAL_CreateGroup(&cal_paths, 0, false, "Paths");
  CAL_SetGroupColor(cal_paths, 0, 0, 1);

  CAL_CreateGroup(&cal_ellipse, 0, false, "Ellipse");
  CAL_SetGroupColor(cal_ellipse, 0, 1, 0, 0.7);
  CAL_CreateSphere(cal_ellipse, 1, 0, 0, 0);
  
  CAL_CreateGroup(&cal_environment, 0, false, "Environment");
  CAL_CreateGroup(&cal_obstacles, cal_environment, false, "Obstacles");

  CAL_SetGroupColor(cal_obstacles, 0, 0, 1);
  /*CAL_CreateBox(cal_obstacles, 0.4, 0.4, 1, cam[0][0], cam[0][1], cam[0][2]+0.5);
  CAL_CreateBox(cal_obstacles, 3, 3, 0.2, cam[0][0], cam[0][1], cam[0][2]+1);
  CAL_CreateBox(cal_obstacles, 0.4, 0.4, 1, cam[1][0], cam[1][1], cam[0][2]+0.5);
  CAL_CreateBox(cal_obstacles, 3, 3, 0.2, cam[1][0], cam[1][1], cam[0][2]+1);*/

  CAL_CreateBox(cal_obstacles, 0.4, 1, 0.4, cam[0][0], cam[0][1]-0.5, cam[0][2]);
  CAL_CreateBox(cal_obstacles, 3, 0.2, 3, cam[0][0], cam[0][1]-1, cam[0][2]);
  CAL_CreateBox(cal_obstacles, 0.4, 1, 0.4, cam[1][0], cam[1][1]-0.5, cam[0][2]);
  CAL_CreateBox(cal_obstacles, 3, 0.2, 3, cam[1][0], cam[1][1]-1, cam[0][2]);

    
  /*int np[1] = {5};
  float points[15] = {x_min[0][0], x_min[0][1], 0,
                      x_min[0][0], x_max[0][1], 0,
                      x_max[0][0], x_max[0][1], 0,
                      x_max[0][0], x_min[0][1], 0,
                      x_min[0][0], x_min[0][1], 0};
  CAL_CreatePolyline(cal_obstacles, 1, np, points);*/
  /*CAL_CreateBox(cal_obstacles, 21, 1, 0.5, 0, -10, 0);
  CAL_CreateBox(cal_obstacles, 1, 21, 0.5, -10, 0, 0);
  CAL_CreateBox(cal_obstacles, 1, 21, 0.5, 10, 0, 0);
  CAL_CreateBox(cal_obstacles, 21, 1, 0.5, 0, 10, 0);*/

  /*for (int i = 0; i < Z_DIM; ++i) {
    CAL_CreateSphere(cal_obstacles, 0.2, (float) beacon[i][0], (float) beacon[i][1], 1);
  }*/

  std::cerr << CAL_LoadScene("articulated.xml", cal_obstacles);
  
  CAL_GetID(&(joint_group[0]), "joint1");
  CAL_GetID(&(joint_group[1]), "joint2");
  CAL_GetID(&(joint_group[2]), "joint3");
  CAL_GetID(&(joint_group[3]), "joint4");
  CAL_GetID(&(joint_group[4]), "joint5");
  CAL_GetID(&(joint_group[5]), "joint6");

  null_joint[0] = 0;
  null_joint[1] = 5.78;
  null_joint[2] = 0;
  null_joint[3] = 0;
  null_joint[4] = 0;
  null_joint[5] = 0;


  showState(start);
  /*int k = 0;

  
  std::cin >> k;
  showState(goal);

  std::cin >> k;

  Matrix<X_DIM> x = start;

  while (k != -1) {
    std::cin >> k;
    double angle;
    std::cin >> angle;
    x[k] = angle;
    showState(x); 
  }*/

}

template <size_t size>
Matrix<size> sampleGaussian(const Matrix<size>& mean, const Matrix<size, size>& var) {
  Matrix<size> sample;
  for (int j = 0; j < size; ++j) {
    sample[j] = normal();
  }
  Matrix<size, size> SVec, SVal;
  jacobi(var, SVec, SVal);
  for (int i = 0; i < size; ++i) {
    SVal(i,i) = sqrt(SVal(i,i));
  }
  return SVec * SVal * sample + mean;
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



double normal() {
  double u_1 = 0;
  while (u_1 == 0) {
    u_1 = random();
  }
  double u_2 = 0;
  while (u_2 == 0) {
    u_2 = random();
  }
  return sqrt(-2*log(u_1)) * sin(2*M_PI*u_2);
}

unsigned int randomInt(unsigned int upper_bound) {
  unsigned int i;
  rand_s(&i);
  return i % upper_bound;
}

/*double computeProbability(const Matrix<DIM>& pos, int k) {
  double prob = 1;
  Matrix<DIM, DIM> P = R[robot_nr][k].subMatrix<DIM,DIM>(0,0);
  for (int r = 0; r < robot_nr; ++r) {
    int kr = k;
    if (kr >= (int) robot_paths[r].size()) {
      kr = (int) robot_paths[r].size() - 1;
    }
    Matrix<DIM,DIM> Pr = R[r][kr].subMatrix<DIM,DIM>(0,0);
    
    Matrix<DIM,DIM> invR = !(P + Pr);
    Matrix<DIM> p_ij = robot_paths[r][kr].x.subMatrix<DIM,1>(0,0) - pos;
    double colprob = 0;
    for (int s = 0; s < (int) disk_samples.size(); ++s) {
      colprob += exp(-0.5*tr( ~(disk_samples[s] - p_ij) * invR * (disk_samples[s] - p_ij) ));
    }
    colprob *= (car_l * car_l / ((double) disk_samples.size() * 2 * sqrt(det(P + Pr))));
    if (colprob > 1) {
      colprob = 1;
    }
    prob *= (1 - colprob);
  }

  return prob;
}*/

void preprocess(const std::vector<PathNode>& path) {
  int l = (int) path.size();

  A.resize(l); // process matrix
  B.resize(l); // input matrix
  V.resize(l); // process noise matrix
  H.resize(l); // measurement matrix
  W.resize(l); // measurement noise matrix
  L.resize(l); // feedback matrix
  K.resize(l); // Kalman-gain matrix
  T.resize(l);
  R.resize(l);

  // Initialization
  C = identity<X_DIM>();
  D = identity<U_DIM>();

  P_0 = identity<X_DIM>() * 0.01;

  M = identity<U_DIM>() * 0.01;
  N = identity<Z_DIM>() * 0.01;
  
  // Jacobians
  for (int k = 0; k < l; ++k) {
    double a0 = path[k].x[0]; double a1 = path[k].x[1]; double a2 = path[k].x[2]; double a3 = path[k].x[3]; double a4 = path[k].x[4]; double a5 = path[k].x[5];
    Matrix<DIM> g_path = g(path[k].x);

    A[k] = identity<X_DIM>(); // df/dx
    B[k] = dt * identity<U_DIM>(); // df/du
    V[k] = B[k]; // df/dm (noise)

    // T: Jacobian of end-effector position w.r.t joint angles (3x6 matrix)
    T[k](0,0) = cos(a0)*(cos(a1)*(sin(a2)*(cos(a4)*l4+l3)+cos(a2)*cos(a3)*sin(a4)*l4)+sin(a1)*(cos(a2)*(cos(a4)*l4+l3)-sin(a2)*cos(a3)*sin(a4)*l4+l2))-sin(a0)*sin(a3)*sin(a4)*l4;
    T[k](0,1) = sin(a0)*(cos(a1)*(cos(a2)*(cos(a4)*l4+l3)-sin(a2)*cos(a3)*sin(a4)*l4+l2)-sin(a1)*(sin(a2)*(cos(a4)*l4+l3)+cos(a2)*cos(a3)*sin(a4)*l4));
    T[k](0,2) = sin(a0)*(sin(a1)*(-sin(a2)*(cos(a4)*l4+l3)-cos(a2)*cos(a3)*sin(a4)*l4)+cos(a1)*(cos(a2)*(cos(a4)*l4+l3)-sin(a2)*cos(a3)*sin(a4)*l4));
    T[k](0,3) = sin(a0)*(sin(a1)*sin(a2)*sin(a3)*sin(a4)*l4-cos(a1)*cos(a2)*sin(a3)*sin(a4)*l4)+cos(a0)*cos(a3)*sin(a4)*l4;
    T[k](0,4) = sin(a0)*(cos(a1)*(cos(a2)*cos(a3)*cos(a4)*l4-sin(a2)*sin(a4)*l4)+sin(a1)*(-cos(a2)*sin(a4)*l4-sin(a2)*cos(a3)*cos(a4)*l4))+cos(a0)*sin(a3)*cos(a4)*l4;
    T[k](0,5) = 0;

    T[k](1,0) = 0;
    T[k](1,1) = -cos(a1)*(sin(a2)*(cos(a4)*l4+l3)+cos(a2)*cos(a3)*sin(a4)*l4)-sin(a1)*(cos(a2)*(cos(a4)*l4+l3)-sin(a2)*cos(a3)*sin(a4)*l4+l2);
    T[k](1,2) = cos(a1)*(-sin(a2)*(cos(a4)*l4+l3)-cos(a2)*cos(a3)*sin(a4)*l4)-sin(a1)*(cos(a2)*(cos(a4)*l4+l3)-sin(a2)*cos(a3)*sin(a4)*l4);
    T[k](1,3) = cos(a1)*sin(a2)*sin(a3)*sin(a4)*l4+sin(a1)*cos(a2)*sin(a3)*sin(a4)*l4;
    T[k](1,4) = cos(a1)*(-cos(a2)*sin(a4)*l4-sin(a2)*cos(a3)*cos(a4)*l4)-sin(a1)*(cos(a2)*cos(a3)*cos(a4)*l4-sin(a2)*sin(a4)*l4);
    T[k](1,5) = 0;
   
    T[k](2,0) = -sin(a0)*(cos(a1)*(sin(a2)*(cos(a4)*l4+l3)+cos(a2)*cos(a3)*sin(a4)*l4)+sin(a1)*(cos(a2)*(cos(a4)*l4+l3)-sin(a2)*cos(a3)*sin(a4)*l4+l2))-cos(a0)*sin(a3)*sin(a4)*l4;
    T[k](2,1) = cos(a0)*(cos(a1)*(cos(a2)*(cos(a4)*l4+l3)-sin(a2)*cos(a3)*sin(a4)*l4+l2)-sin(a1)*(sin(a2)*(cos(a4)*l4+l3)+cos(a2)*cos(a3)*sin(a4)*l4));
    T[k](2,2) = cos(a0)*(sin(a1)*(-sin(a2)*(cos(a4)*l4+l3)-cos(a2)*cos(a3)*sin(a4)*l4)+cos(a1)*(cos(a2)*(cos(a4)*l4+l3)-sin(a2)*cos(a3)*sin(a4)*l4));
    T[k](2,3) = cos(a0)*(sin(a1)*sin(a2)*sin(a3)*sin(a4)*l4-cos(a1)*cos(a2)*sin(a3)*sin(a4)*l4)-sin(a0)*cos(a3)*sin(a4)*l4;
    T[k](2,4) = cos(a0)*(cos(a1)*(cos(a2)*cos(a3)*cos(a4)*l4-sin(a2)*sin(a4)*l4)+sin(a1)*(-cos(a2)*sin(a4)*l4-sin(a2)*cos(a3)*cos(a4)*l4))-sin(a0)*sin(a3)*cos(a4)*l4;
    T[k](2,5) = 0;
    
    // H = dh/dx
    for (int i = 0; i < X_DIM; ++i) {
      /*H[k](0,i) = (T[k](0,i)*(g_path[2] - cam[0][2]) - (g_path[0] - cam[0][0])*T[k](2,i)) / ((g_path[2] - cam[0][2])*(g_path[2] - cam[0][2]));
      H[k](1,i) = (T[k](1,i)*(g_path[2] - cam[0][2]) - (g_path[1] - cam[0][1])*T[k](2,i)) / ((g_path[2] - cam[0][2])*(g_path[2] - cam[0][2]));
      H[k](2,i) = (T[k](0,i)*(g_path[2] - cam[1][2]) - (g_path[0] - cam[1][0])*T[k](2,i)) / ((g_path[2] - cam[1][2])*(g_path[2] - cam[1][2]));
      H[k](3,i) = (T[k](1,i)*(g_path[2] - cam[1][2]) - (g_path[1] - cam[1][1])*T[k](2,i)) / ((g_path[2] - cam[1][2])*(g_path[2] - cam[1][2]));*/
      
      H[k](0,i) = (T[k](0,i)*-(g_path[1] - cam[0][1]) - (g_path[0] - cam[0][0])*-T[k](1,i)) / ((g_path[1] - cam[0][1])*(g_path[1] - cam[0][1]));
      H[k](1,i) = (T[k](2,i)*-(g_path[1] - cam[0][1]) - (g_path[2] - cam[0][2])*-T[k](1,i)) / ((g_path[1] - cam[0][1])*(g_path[1] - cam[0][1]));
      H[k](2,i) = (T[k](0,i)*-(g_path[1] - cam[1][1]) - (g_path[0] - cam[1][0])*-T[k](1,i)) / ((g_path[1] - cam[1][1])*(g_path[1] - cam[1][1]));
      H[k](3,i) = (T[k](2,i)*-(g_path[1] - cam[1][1]) - (g_path[2] - cam[1][2])*-T[k](1,i)) / ((g_path[1] - cam[1][1])*(g_path[1] - cam[1][1]));
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

  for (int k = 1; k < l; ++k) {
    F[k].insert(0,0,     A[k]);           F[k].insert(0,X_DIM,     B[k]*L[k-1]);
    F[k].insert(X_DIM,0, K[k]*H[k]*A[k]); F[k].insert(X_DIM,X_DIM, A[k] + B[k]*L[k-1] - K[k]*H[k]*A[k]);
  
    G[k].insert(0,0,     V[k]);           G[k].insert(0,U_DIM,     zeros<X_DIM,Z_DIM>());
    G[k].insert(X_DIM,0, K[k]*H[k]*V[k]); G[k].insert(X_DIM,U_DIM, K[k]*W[k]);

    R[k] = F[k]*R[k-1]*~F[k] + G[k]*Q*~G[k];
  }


}

/*double simulate(const std::vector<PathNode>& path, int num_samples, bool vis = false) {
  int l = (int) path.size();

  preprocess(path);

  // Simulate
  Matrix<X_DIM> x_est;
  Matrix<X_DIM> x_true;
  Matrix<U_DIM> u;
  Matrix<U_DIM> m;
  Matrix<X_DIM> x_true_old;
  Matrix<Z_DIM> n;
  Matrix<Z_DIM> z;
  Matrix<X_DIM, X_DIM> P;
  
  int fail = 0;

  int* cal_samples;
  if (vis) {
    cal_samples = new int[num_samples];
  }

  for (int s = 0; s < num_samples; ++s) {
    if (vis) {
      CAL_CreateGroup(&(cal_samples[s]), 0, false);
      CAL_SetGroupColor(cal_samples[s], 0, 0, 1);
      CAL_SetGroupMotionOptions(cal_samples[s], CAL_NOINTERPOLATION);
      int obj;
      CAL_CreateCylinder(cal_samples[s], 0.075, 0.075, 0, 0, 0, &obj);
      CAL_SetObjectOrientation(obj, (float) (M_PI*0.5), 0, 0);
      //CAL_CreateSphere(cal_samples[s], 0.075, 0, 0, 0);
    }

    /*P = P_0;
    x_est = path[0].x;
    x_true = sampleGaussian(path[0].x, P_0);
    
    for (int k = 1; k < l; ++k) {
      u = path[k-1].u + L[k-1]*(x_est - path[k-1].x);
      m = sampleGaussian(zeros<U_DIM,1>(), M);

      x_true_old = x_true;

      x_true = f(x_true, u, m);

      if (vis) {
        float pos[3] = {(float)x_true_old[0], (float)x_true_old[1], 0};
        CAL_AddGroupKeyState(cal_samples[s], (k-1)*(float) dt, pos);
      } else {
        int col;
        CAL_CheckLineCollision(cal_environment, (float) x_true_old[0], (float) x_true_old[1], 0, (float) x_true[0], (float) x_true[1], 0, false, &col);
        if (col != 0) {
          ++fail;
          break;
        }
      }

      // Jacobian around current estimate
      A[k](0,0) = 1; A[k](0,1) = 0; A[k](0,2) = -dt*x_est[3]*sin(x_est[2]); A[k](0,3) = dt*cos(x_est[2]);       A[k](0,4) = 0;
      A[k](1,0) = 0; A[k](1,1) = 1; A[k](1,2) = dt*x_est[3]*cos(x_est[2]);  A[k](1,3) = dt*sin(x_est[2]);       A[k](1,4) = 0;
      A[k](2,0) = 0; A[k](2,1) = 0; A[k](2,2) = 1;                          A[k](2,3) = dt*tan(x_est[4])/car_l; A[k](2,4) = dt*x_est[3]*(1+tan(x_est[4])*tan(x_est[4]))/car_l;
      A[k](3,0) = 0; A[k](3,1) = 0; A[k](3,2) = 0;                          A[k](3,3) = 1;                      A[k](3,4) = 0;
      A[k](4,0) = 0; A[k](4,1) = 0; A[k](4,2) = 0;                          A[k](4,3) = 0;                      A[k](4,4) = 1;

      // Process update
      x_est = f(x_est, u, zeros<U_DIM,1>());
      P = A[k]*P*~A[k] + V[k]*M*~V[k];

      // Measurement update
      n = sampleGaussian(zeros<Z_DIM,1>(), N);
      z = H[k]*x_true + W[k]*n;

      K[k] = P*~H[k]*!(H[k]*P*~H[k] + W[k]*N*~W[k]);
      x_est = x_est + K[k] * (z - H[k]*x_est);
      P = (identity<X_DIM>() - K[k]*H[k])*P;
    }*/
  
/*
    P = P_0;
    x_est = zeros<X_DIM, 1>();
    x_true = sampleGaussian(zeros<X_DIM, 1>(), P_0);
    
    for (int k = 1; k < l; ++k) {
      u = L[k-1]*x_est;
      m = sampleGaussian(zeros<U_DIM,1>(), M);

      x_true_old = x_true;

      x_true = A[k] * x_true + B[k] * u + V[k] * m;

      if (vis) {
        float pos[3] = {x_true_old[0] + path[k-1].x[0], x_true_old[1] + path[k-1].x[1], 0};
        CAL_AddGroupKeyState(cal_samples[s], (k-1)*dt, pos);
      } else {
        int col;
        CAL_CheckLineCollision(cal_environment, (float) x_true_old[0] + path[k-1].x[0], (float) x_true_old[1] + path[k-1].x[1], 0, (float) x_true[0] + path[k].x[0], (float) x_true[1] + path[k].x[1], 0, false, &col);
        if (col != 0) {
          ++fail;
          break;
        }
      }

      // Process update
      x_est = A[k] * x_est + B[k] * u;
      
      // Measurement update
      n = sampleGaussian(zeros<Z_DIM,1>(), N);
      z = H[k]*x_true + W[k]*n;

      x_est = x_est + K[k] * (z - H[k]*x_est);
    }
  
  
  }

  return 1 - (double) fail / num_samples;
}*/

/*std::pair<double, double> computeQuality(const std::vector<PathNode>& path) {
  int l = (int) path.size();

  preprocess(path);

  std::vector<Matrix<2*X_DIM, 2*X_DIM> > F(l);
  std::vector<Matrix<2*X_DIM, U_DIM+Z_DIM> > G(l);
  R[robot_nr].resize(l);

  // Combination of LQR and Kalman
  Matrix<U_DIM+Z_DIM, U_DIM+Z_DIM> Q = zeros<U_DIM+Z_DIM, U_DIM+Z_DIM>();
  Q.insert(0,0, M); Q.insert(U_DIM, U_DIM, N);

  R[robot_nr][0].reset();
  R[robot_nr][0].insert(0,0,P_0);

  for (int k = 1; k < l; ++k) {
    F[k].insert(0,0,     A[k]);           F[k].insert(0,X_DIM,     B[k]*L[k-1]);
    F[k].insert(X_DIM,0, K[k]*H[k]*A[k]); F[k].insert(X_DIM,X_DIM, A[k] + B[k]*L[k-1] - K[k]*H[k]*A[k]);
  
    G[k].insert(0,0,     V[k]);           G[k].insert(0,U_DIM,     zeros<X_DIM,Z_DIM>());
    G[k].insert(X_DIM,0, K[k]*H[k]*V[k]); G[k].insert(X_DIM,U_DIM, K[k]*W[k]);

    R[robot_nr][k] = F[k]*R[robot_nr][k-1]*~F[k] + G[k]*Q*~G[k];
  }

  std::pair<double, double> quality;
  quality.first = 1;
  for (int k = 0; k < l; ++k) {
    double prob = computeProbability(path[k].x.subMatrix<DIM,1>(0,0), k);
    quality.first *= prob;
  }
  quality.second = 0;
  for (int k = 0; k < l; ++k) {
    quality.second -= tr(R[robot_nr][k].subMatrix<DIM,DIM>(0,0));
  }

  return quality;
    
}*/

std::vector<std::vector<PathNode> > smoothPaths(const std::vector<std::vector<PathNode> >& paths) {
  Matrix<X_DIM+U_DIM, X_DIM+U_DIM> smoothA, smoothM;
  Matrix<X_DIM, X_DIM> smoothN;
  Matrix<X_DIM, X_DIM+U_DIM> smoothH;

  smoothA.insert(0,0, identity<X_DIM>());
  smoothA.insert(X_DIM,0, zeros<X_DIM, X_DIM>());
  smoothA.insert(0, X_DIM, dt * identity<X_DIM>());
  smoothA.insert(X_DIM, X_DIM, identity<U_DIM>());

  smoothH.reset();
  smoothH.insert(0,0,identity<X_DIM>());
  
  smoothN = identity<X_DIM>();

  smoothM.reset();
  smoothM.insert(X_DIM, X_DIM, 1*identity<U_DIM>());

  std::vector<std::vector<PathNode> > smoothPaths(paths.size());

  for (int j = 0; j < paths.size(); ++j) {
    
    // N is number of examples, F is z_dim x z_dim, H is N*y_dim x z_dim, Q is z_dim x z_dim, R = N*y_dim x N*y_dim, y = [y_dim x 1], z = [z_dim x 1], P = [z_dim x z_dim]
    // Function returns most likely distribution over z's with mean z[i] and covariance P[i]
    // Also updates parameters Q and R

    int T = paths[j].size();

    std::vector<Matrix<X_DIM+U_DIM, 1> > apriori_z(T), aposteriori_z(T), final_z(T);
    std::vector<Matrix<X_DIM+U_DIM, X_DIM+U_DIM> > apriori_P(T), aposteriori_P(T), final_P(T);
    std::vector<PathNode> smoothPath(T);

    apriori_z[0].reset();
    apriori_z[0].insert(0,0, paths[j][0].x);

    apriori_P[0].reset();
    apriori_P[0].insert(X_DIM,X_DIM, 10000000 * identity<U_DIM>());

    // Kalman Smoother; Forward pass
    for (int t = 0; t < T; ++t) {
      // prediction step
      if (t != 0) {
        // Standard forward prediction
        apriori_z[t] = smoothA * aposteriori_z[t-1];
        apriori_P[t] = smoothA * aposteriori_P[t-1] * ~smoothA + smoothM;
      }

      

      // Update step
      Matrix<X_DIM+U_DIM, X_DIM> smoothK = apriori_P[t] * ~smoothH * !(smoothH * apriori_P[t] * ~smoothH + (t != T - 1 ? smoothN : 0.0*identity<X_DIM>()) );
      
      aposteriori_z[t] = apriori_z[t] + smoothK * (paths[j][t].x - smoothH * apriori_z[t]);
      aposteriori_P[t] = (identity<X_DIM+U_DIM>() - smoothK * smoothH) * apriori_P[t];   
    }

    
    // Backward pass and smoother in one
    final_z[T - 1] = aposteriori_z[T - 1];
    final_P[T - 1] = aposteriori_P[T - 1];

    smoothPath[T-1].x = final_z[T-1].subMatrix<X_DIM,1>(0,0);
    smoothPath[T-1].u = final_z[T-1].subMatrix<U_DIM,1>(X_DIM,0);

    for (int t = T - 2; t >= 0; --t) {
      Matrix<X_DIM+U_DIM,X_DIM+U_DIM> smoothL = aposteriori_P[t] * ~smoothA * !apriori_P[t+1];
      final_P[t] = aposteriori_P[t] - smoothL * (apriori_P[t+1] - final_P[t+1]) * ~smoothL;
      final_z[t] = aposteriori_z[t] + smoothL * (final_z[t+1] - apriori_z[t+1]);

      smoothPath[t].x = final_z[t].subMatrix<X_DIM,1>(0,0);
      smoothPath[t].u = final_z[t].subMatrix<U_DIM,1>(X_DIM,0);
    }

    smoothPaths[j] = smoothPath;
  }
  return smoothPaths;
}

void showPath(const std::vector<std::vector<PathNode> >& paths) {
  for (int j = 0; j < paths.size(); ++j) {
    int np[1] = {(int) paths[j].size()};
    float* p = new float[3*np[0]];
    for (int i = 0; i < (int) paths[j].size(); ++i) {
      Matrix<DIM> q = g(paths[j][i].x);
      p[3*i] = (float) q[0]; p[3*i+1] = q[1]; p[3*i+2] = q[2];
    }
    CAL_CreatePolyline(cal_paths, 1, np, p, 0);
    delete[] p;
  }

  //int k;
  //std::cin >> k;
  
  //
  /*CAL_SetGroupMotionOptions(joint_group[0], CAL_NOINTERPOLATION);
  CAL_SetGroupMotionOptions(joint_group[1], CAL_NOINTERPOLATION);
  CAL_SetGroupMotionOptions(joint_group[2], CAL_NOINTERPOLATION);
  CAL_SetGroupMotionOptions(joint_group[3], CAL_NOINTERPOLATION);
  CAL_SetGroupMotionOptions(joint_group[4], CAL_NOINTERPOLATION);
  CAL_SetGroupMotionOptions(joint_group[5], CAL_NOINTERPOLATION);*/

  preprocess(paths[0]);

  for (int i = 0; i < (int) paths[0].size(); ++i) {
    Matrix<DIM,DIM> EVec, EVal;

    jacobi(T[i]*R[i].subMatrix<X_DIM,X_DIM>(0,0)*~T[i], EVec, EVal);
    Matrix<4,1> q = quatFromRot(EVec);
    Matrix<DIM> p = g(paths[0][i].x);

    float pos[3] = {(float) p[0], (float) p[1], (float) p[2]};
    float quat[4] = {(float) q(0,0), (float) q(1,0), (float) q(2,0), (float) q(3,0)};
    float scale[3] = {(float) 2*sqrt(EVal(0,0)), (float) (2*sqrt(EVal(1,1))), (float) (2*sqrt(EVal(2,2)))};
    CAL_AddGroupKeyState(cal_ellipse, (float) (i * dt), pos, quat, scale, true);
    for (int j = 0; j < 6; ++j) {
      double angle = paths[0][i].x[j] + null_joint[j];
      if (angle < 0) {
        angle += 2*M_PI;
      } else if (angle >= 2*M_PI) {
        angle -= 2*M_PI;
      }

      float rot[3] = {0, angle, 0}; 
      CAL_AddGroupKeyState(joint_group[j], (float) (i * dt), 0, rot, 0, false);
    }
  }
  /*CAL_SetGroupMotionOptions(cal_ellipse, CAL_NOINTERPOLATION);
  CAL_SetGroupMotionOptions(joint_group[0], CAL_NOINTERPOLATION);
  CAL_SetGroupMotionOptions(joint_group[1], CAL_NOINTERPOLATION);
  CAL_SetGroupMotionOptions(joint_group[2], CAL_NOINTERPOLATION);
  CAL_SetGroupMotionOptions(joint_group[3], CAL_NOINTERPOLATION);
  CAL_SetGroupMotionOptions(joint_group[4], CAL_NOINTERPOLATION);
  CAL_SetGroupMotionOptions(joint_group[5], CAL_NOINTERPOLATION);*/
  
}

void readPaths(std::vector<std::vector<PathNode> >& paths) {
  int numPaths;
  std::cin >> numPaths;
  paths.resize(numPaths);
  for (int i = 0; i < numPaths; ++i) {
    int length;
    std::cin >> length;
    std::vector<PathNode> path(length);
    for (int j = 0; j < length; ++j) {
      PathNode node;
      std::cin >> node.x;
      std::cin >> node.u;
      path[j] = node;
    }
    paths[i] = path;
  }
}

void rrt() {
  std::vector<RRTNode> rrtTree;
  //CAL_EmptyGroup(cal_rrt);

  RRTNode startNode;
  startNode.x = start;
  startNode.p = g(start);
  
  rrtTree.push_back(startNode);

  while (tr(~(rrtTree.back().p - goal) * (rrtTree.back().p - goal)) > goalRadius*goalRadius) {
    Matrix<DIM> sample;
    sample[0] = random() * 44 - 22;
    sample[1] = random() * 22;
    sample[2] = random() * 44 - 22;
    if (tr(~sample * sample) > 22*22) {
      continue;
    }
                
    int node = -1;
    double mindist = 9e99;
    for (int j = 0; j < (int) rrtTree.size(); ++j) {
      double ddist = tr(~(rrtTree[j].p - sample) * (rrtTree[j].p - sample));
      if (ddist < mindist) {
        mindist = ddist;
        node = j;
      }
    }
    if (node == -1) {
      continue;
    }
        
    Matrix<U_DIM> input;
    for (int i = 0; i < U_DIM; ++i) {
     input[i] = random() * (u_max[i] - u_min[i]) + u_min[i];
    }
    Matrix<X_DIM> x_new = f(rrtTree[node].x, input, zeros<U_DIM,1>());
    
    bool valid = true;
    for (int i = 0; i < X_DIM; ++i) {
      if (x_new[i] < x_min[i] || x_new[i] > x_max[i]) {
        valid = false;
        break;
      }
    }
    if (!valid) {
      continue;
    }
    

    RRTNode newnode;
    newnode.x = x_new;
    newnode.p = g(x_new);
    newnode.u = input;
    newnode.bp = node;
    
    rrtTree.push_back(newnode);

    
    /*Matrix<DIM> last_p = g(rrtTree[node].x);
    int np[1] = {2};
    float p[6] = {(float) last_p[0], last_p[1], last_p[2], (float) q[0], (float) q[1], (float) q[2]};
    CAL_CreatePolyline(cal_rrt, 1, np, p);*/
  }

  //int k;
  //std::cin >> k;

  std::deque<PathNode> path;

  int i = (int) rrtTree.size() - 1;
  PathNode node;
  node.x = rrtTree[i].x;
  node.u = zeros<U_DIM, 1>();
  path.push_front(node);
  
  while (i != 0) {
    node.u = rrtTree[i].u;
    i = rrtTree[i].bp;
    node.x = rrtTree[i].x;
    
    path.push_front(node);
  }

  std::cout << path.size() << std::endl;
  for (int i = 0; i < (int) path.size(); ++i) {
    std::cout << ~path[i].x;
    std::cout << ~path[i].u;
  }
}

int _tmain(int argc, _TCHAR* argv[])
{
  clock_t startTime, endTime;

  // Create Obstacles in Callisto
  CAL_Initialisation (true,true,true);
  initEnvironment();

  // Construct Paths
  startTime = clock();
  int numPaths = 1000;
  std::cout << numPaths << std::endl;
  for (int i = 0; i < numPaths; ++i) {
    rrt();
    std::cerr << i << " ";
  }
  CAL_End();
  endTime = clock();
  std::cerr << (double) (endTime - startTime) / CLOCKS_PER_SEC << std::endl;
  /*return 0;*/


  // Read Paths
  //std::vector<std::vector<PathNode> > paths;
  //readPaths(paths);

  paths = smoothPaths(paths);

  // Show Path
  /*showPath(paths);
  int k;
  std::cin >> k;
  CAL_End();
	return 0;*/

  

  
  // Compute Optimal Path
  
  
  startTime = clock();

  int bestPath = 0;
  double bestQuality = DBL_MAX;

  for (int i = 0; i < (int) paths.size(); ++i) {
    preprocess(paths[i]);
    /*double quality = 0;
    for (int j = 0; j < paths[i].size(); ++j) {
      quality += tr(T[j]*R[j].subMatrix<X_DIM,X_DIM>(0,0)*~T[j]);
    }*/
    double quality = tr(T[paths[i].size()-1]*R[paths[i].size()-1].subMatrix<X_DIM,X_DIM>(0,0)*~T[paths[i].size()-1]);
    //double quality = simulate(paths[i], 10000);
    std::cerr << i << " " << quality << std::endl;
    
    if (quality < bestQuality) {
      bestPath = i;
      bestQuality = quality;
    }
  }

  endTime = clock();

  std::cerr << (double) (endTime - startTime) / CLOCKS_PER_SEC << std::endl;

  std::cout << "1" << std::endl << paths[bestPath].size() << std::endl;
  for (int i = 0; i < (int) paths[bestPath].size(); ++i) {
    std::cout << ~paths[bestPath][i].x;
    std::cout << ~paths[bestPath][i].u;
  }

  std::cerr << bestPath << " " << bestQuality << std::endl;
  
  CAL_End();
	return 0;
}

