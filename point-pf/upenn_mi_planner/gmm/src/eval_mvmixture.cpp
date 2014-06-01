/** Copyright (C) 2012 Benjamin Charrow
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 */

#include <iostream>

#include <Eigen/Core>

#include <timeit/timeit.hpp>

#include "gmm/gmm.hpp"

using namespace Eigen;
using namespace std;
using namespace timeit;
using namespace rf;

VectorXd runTrial(double c_start, double c_stop, int nsteps, EntropyApproximator<2> *ea) {
  VectorXd results(nsteps);
  for (int step = 0; step < nsteps; ++step) {
    double c = c_start + static_cast<double>(step) / nsteps * (c_stop - c_start);
    
    GaussMM<2>::Ptr gmm = new GaussMM<2>();
    gmm->addComponent(Vector2d::Zero(), Vector2d(0.16, 1.0), 0.2);
    gmm->addComponent(Vector2d(3.0, 2.0), Vector2d(1.0, 0.16), 0.2);
    gmm->addComponent(Vector2d(1.0, -0.5), Vector2d(0.5, 0.5), 0.2);
    gmm->addComponent(Vector2d(2.5, 1.5), Vector2d(0.5, 0.5), 0.2);
    gmm->addComponent(Vector2d(c, c), Vector2d(0.5, 0.5), 0.2);

    results(step) = ea->evaluate(*gmm);

    delete gmm;
  }
  return results;
}

int main(int argc, char **argv) {
  double c_start = -3;
  double c_stop = 3;
  int nsteps = 100;

  // c in 0th column
  // Numerical in 1st column
  // taylor 2 in 2nd column
  // taylor 2 w/ splits in 3rd column
  MatrixXd results(nsteps, 4);

  for (int step = 0; step < nsteps; ++step) {
    double c = c_start + static_cast<double>(step) / nsteps * (c_stop - c_start);
    results(step, 0) = c;
  }

  {
    NumericalApproximator<2> na(3000);
    tic();
    results.col(1) = runTrial(c_start, c_stop, nsteps, &na);
    toc("Numerical:  ");
  }
  {
    TaylorApproximator<2> ta(2, 0);
    tic();
    results.col(2) = runTrial(c_start, c_stop, nsteps, &ta);
    toc("Taylor 2nd: ");
  }
  {
    TaylorApproximator<2> ta(2, 20);
    tic();
    results.col(3) = runTrial(c_start, c_stop, nsteps, &ta);
    toc("Taylor 2nd w/ 20 splits: ");
  }
  cout << results << endl;
}
