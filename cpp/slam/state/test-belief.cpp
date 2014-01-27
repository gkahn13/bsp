#include <vector>
#include <iomanip>

#include "../slam.h"


int main() {
  //Q_DIM -> e
  //R_DIM -> 2*num_landmarks (2)
  //R_matrix -> block diagonal for each landmark
  Matrix<X_DIM, X_DIM> sqrtsigma0 = zeros<X_DIM, X_DIM>();
  R = zeros<R_DIM, R_DIM>();
  for (int i=0; i<NUM_LANDMARKS; ++i) {
    R(i*2, i*2) = 10;
    R(i*2+1, i*2+1) = 10;//std::pow(3.14159265358979/180, 2);
    sqrtsigma0(3+i*2, 3+i*2) = R(i*2, i*2);
    sqrtsigma0(3+i*2+1, 3+i*2+1) = R(i*2+1, i*2+1);
  }
  Q = zeros<Q_DIM, Q_DIM>();
  Q(0,0) = 3;
  Q(1,1) = 3;
  Q(2,2) = .3;//std::pow(3.141592653589793238462643383/60,2);

  init_landmarks();
  std::vector<Matrix<B_DIM> > B(T);
  Matrix<B_DIM> b;
  B[0] = b;
  std::vector< Matrix<2> > U(T-1);
  Matrix<2> u;
  u[0] = 1; u[1] = 3.14159265358979/4;
  Matrix<2> ustop;
  ustop[0] = 0; ustop[1] = 0;
  Matrix<2> uforward;
  uforward[0] = 1; uforward[1] = 0;
  Matrix<2> u2;
  u2[0] = 1; u2[1] = -3.14159265358979/4;


  for (int i=0;i<4;++i){
    U[i] = u;
  }
  U[4] = uforward;
  U[5] = uforward;
  for (int i=6;i<19;++i){
    U[i] = ustop;
  }
  
  Matrix<X_DIM> x0 = zeros<X_DIM,1>();
  x0[0]=0; x0[1]=0; x0[2]=0; 
  for (int i=0; i<NUM_LANDMARKS; ++i) {
    x0[3+2*i] = actualLandmarks[i][0];
    x0[4+2*i] = actualLandmarks[i][1];
  }
  vec(x0,sqrtsigma0, B[0]);

  for (int t=0; t<T-1; ++t) {
    B[t+1]=beliefDynamics(B[t],U[t]);
    unVec(B[t+1], x0, sqrtsigma0);
  }

  pythonDisplayTrajectory(B);

  return 0;
}
