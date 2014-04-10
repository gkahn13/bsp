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

#include <ros_gsl/random.hpp>

#include "gmm/gmm.hpp"

using namespace Eigen;
using namespace timeit;
using namespace rf;

GaussMM<1>::Ptr buildParzen(VectorXd data) {
  GaussMM<1>::Ptr gmm = new GaussMM<1>();
  double mean = data.sum() / data.size();
  
  Vector1d cov = Vector1d::Zero();
  for (int i = 0; i < data.size(); ++i) {
    cov(0) += (data(i) - mean) * (data(i) - mean);
  }
  cov(0) = (cov(0) / (data.size() - 1)) * pow(4.0 / (3.0 * data.size()), 2.0 / 5.0);
  for (int i = 0; i < data.size(); ++i) {
    Vector1d mean = Vector1d::Constant(data(i));
    gmm->addComponent(mean, cov, 1.0 / data.size());
  }
  return gmm;
}

VectorXd runTrial(const VectorXd &inv_scales, const VectorXd &samples, const VectorXd &outputs,
                  EntropyApproximator<1> *ea) {
  VectorXd results;
  results.resize(inv_scales.size());
  for (int i = 0; i < inv_scales.size(); ++i) {
    GaussMM<1>::Ptr parzen = buildParzen(outputs * inv_scales(i) - samples);
    results(i) = ea->evaluate(*parzen);
    delete parzen;
  }
  return results;
}

int main(int argc, char **argv) {
  // GMM for data
  GaussMM<1> gmm;
  gmm.addComponent(Vector1d::Constant(-1), Matrix1d::Constant(0.025), 0.4);
  gmm.addComponent(Vector1d::Constant(1), Matrix1d::Constant(1), 0.6);
  
  // Param to estimate
  const double param = 2.0;
  // Samples to pass through model
  VectorXd samples = gmm.sample(100);

  double noise_sd = sqrt(0.04);
  VectorXd noises(samples.size());
  for (int i = 0; i < samples.size(); ++i) {
    noises(i) = gsl::normal(noise_sd);
  }

  VectorXd outputs = samples * param + noises;

  // inv_scale in 0th column
  // Numerical in 1st column
  // taylor 0 in 2nd column
  // taylor 2 in 3rd column
  int nevals = 100;
  const double inv_start = -2, inv_stop = 6;
  MatrixXd results(nevals, 4);
  
  for (int ind = 0; ind < nevals; ++ind) {
    results(ind, 0) = inv_start + (inv_stop - inv_start) * static_cast<double>(ind) / nevals;
  }

  {
    tic();
    NumericalApproximator<1> na(3000);
    results.col(1) = runTrial(results.col(0), samples, outputs, &na);
    toc("Numerical:   ");
  }
  {
    tic();
    TaylorApproximator<1> ta(0, 0);
    results.col(2) = runTrial(results.col(0), samples, outputs, &ta);
    toc("Taylor 0th:  ");
  }
  {
    tic();
    TaylorApproximator<1> ta(2, 0);
    results.col(3) = runTrial(results.col(0), samples, outputs, &ta);
    toc("Taylor 2nd:  ");
  }

  std::cout << results << std::endl;
}
