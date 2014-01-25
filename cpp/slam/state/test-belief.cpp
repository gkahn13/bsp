#include <vector>
#include <iomanip>

#include "../slam.h"


int main() {
  Matrix<2> tmp;
  //tmp[0]=553.8; tmp[1]=634.3; 
  tmp[0] = 500; tmp[1] = 0;  

  actualLandmarks.push_back(tmp);
  //  init_landmarks();
  std::vector<Matrix<B_DIM> > B(10);
  Matrix<B_DIM> b;
  B[0] = b;
  std::vector< Matrix<2> > U(10);
  Matrix<2> u;
  u[0] = 1; u[1] = 0;

  for (int i=0;i<7;++i){
    U[i] = u;
  }
  u[0] = 3; u[1] = .1;
  U[7] = u;
  u[0] = 4; u[1] = .1;
  U[8] = u;
  u[0] = 5; u[1] = 0;
  U[9] = u;

  
  Matrix<X_DIM> x0 = zeros<X_DIM,1>();
  x0[0]=0; x0[1]=0; x0[2]=0; x0[3] = 10; x0[4] = 0;
  Matrix<X_DIM, X_DIM> sqrtsigma0 = 10*identity<X_DIM>();
  vec(x0,sqrtsigma0, B[0]);
  obsfunc(x0, zeros<R_DIM,1>());
  
  for (int t=0; t<9; ++t) {
    B[t+1]=beliefDynamics(B[t],U[t]);
    unVec(B[t+1], x0, sqrtsigma0);
  }
  
  //pythonDisplayTrajectory(B, U);

  return 0;
}
