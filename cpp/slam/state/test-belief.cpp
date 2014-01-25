#include <vector>
#include <iomanip>

#include "../slam.h"

#include "util/matrix.h"
#include "util/Timer.h"
#include "util/logging.h"
//#include "util/utils.h"

#include <Python.h>
#include <boost/python.hpp>
#include <boost/filesystem.hpp>

namespace py = boost::python;

int main() {
  std::vector< Matrix<2> >& U(10);
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
  u[0] = 6; u[1] = 0;
  U[10] = u;
  
  std::vector< Matrix<B_DIM> > B(T);

  Matrx<X_DIM> x0 = zeros<X_DIM,1>();
  x0[0]=0; x0[1]=0; x0[2]=0;
  Matrix<X_DIM, X_DIM> sqrtsigma0 = identity<X_DIM>();
  vec(x,sqrtsigma0, B[0]);
  for (int t=0; t<T; ++t) {
    B[t+1]=beliefDynamics(B[t],U[t])
  }
  pythonDisplayTrajectory(B, U);

  return 0;
}
