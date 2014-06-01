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

#include "gmm/gmm.hpp"
using namespace Eigen;
using namespace std;

namespace rf {
  double fast_exp(double val) {
    int tmp2 = (int)(1512775 * val + (1072693248 - 60801));
    double p = 0.0;
    *(1 + (int * )&p) = tmp2;
    return p;
  }
  
  template <>
  double GaussMM<1>::compLikelihood(const Vector1d &vec, int i) const {
    double diff = vec(0) - means_[i](0);
    double quad = diff * diff * invcovs_[i](0);
    return normalizers_[i] * exp(-0.5 * quad);
  }

  template <>
  double GaussMM<2>::compLikelihood(const Vector2d &vec, int i) const {
    double diff0 = vec(0) - means_[i](0);
    double diff1 = vec(1) - means_[i](1);
    double quad = diff0 * diff0 * invcovs_[i](0) + diff1 * diff1 * invcovs_[i](1);
    return normalizers_[i] * exp(-0.5 * quad);
  }
  
  static double eval_gmm_ent_1d(double x, void *g) {
    GaussMM<1>::ConstPtr gmm = static_cast<GaussMM<1>::ConstPtr>(g);
    GaussMM<1>::VectorD vec(1);
    vec(0) = x;
    double likelihood = gmm->likelihood(vec);
    if (likelihood < 1e-50) {
      return 0;
    } else {
      return -likelihood * log(likelihood);
    }
  }
  
  template <>
  Eigen::VectorXd GaussMM<1>::sample(int num) const {
    double *cdf = new double[numComponents()];
    for (size_t ind = 0; ind < numComponents(); ++ind) {
      cdf[ind] = weights_[ind] + (ind > 0 ? cdf[ind - 1] : 0.0);
    }

    std::vector<double> draws(num);
    for (int i = 0; i < num; ++i) {
      draws[i] = gsl::uniform(0, 1);
    }
    sort(draws.begin(), draws.end());

    Eigen::VectorXd samples(num);
    for (int draw = 0, cdf_ind = 0; draw < num; ++draw) {
      while (cdf[cdf_ind] < draws[draw]) {
        cdf_ind++;
      }
      
      samples(draw) = gsl::normal(sqrt(covs_[cdf_ind].value())) + means_[cdf_ind].value();
    }
  
    delete[] cdf;
    return samples;
  }


  static bool eigen_1dlt(const Vector1d &v1, const Vector1d &v2) {
    return v1(0) < v2(0);
  }
  
  template <>
  double NumericalApproximator<1>::evaluate(const GaussMM<1> &gmm) {
    gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(iters_);
    gsl_function F;
    F.function = &eval_gmm_ent_1d;
    F.params = static_cast<void*>(const_cast<GaussMM<1>::Ptr>(&gmm));
    
    double result, abserr;
    double mean_max = (*max_element(gmm.means().begin(), gmm.means().end(), eigen_1dlt))(0);
    double mean_min = (*min_element(gmm.means().begin(), gmm.means().end(), eigen_1dlt))(0);
    double cov_max = (*max_element(gmm.covs().begin(), gmm.covs().end(), eigen_1dlt))(0);

    gsl_integration_qag(&F, mean_min - 10 * sqrt(cov_max), mean_max + 10 * sqrt(cov_max),
                        1e-10, 1e-10, iters_, GSL_INTEG_GAUSS51,
                        workspace, &result, &abserr);
    gsl_integration_workspace_free(workspace);
    return result;
  }

  template<>
  int GaussMM<Eigen::Dynamic>::dim() const {
    return means_.size() != 0 ? means_[0].size() : -1;
  }
}

