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

#include "gtest/gtest.h"

#include "gmm/gmm.hpp"

#include <Eigen/Core>

#include <gsl/gsl_randist.h>

#include <ros_gsl/random.hpp>

#include <timeit/timeit.hpp>

using namespace Eigen;
using namespace std;
using namespace rf;

//=========================== Mixture model tests ===========================//

TEST(GMM, OneStandardGaussian) {
  double mu = 0, sd = 1.0;
  GaussMM<1> gmm;
  gmm.addComponent(VectorXd::Constant(1, mu), Vector1d::Constant(sd * sd), 1.0);

  for (double val = -3; val < 3; val += 0.2) {
    VectorXd vec = VectorXd::Constant(1, val);
    ASSERT_DOUBLE_EQ(gsl_ran_gaussian_pdf(val - mu, sd), gmm.likelihood(vec));
  }
}

TEST(GMM, OneNonstandardGaussian) {
  double mu = -1.0, sd = 4.0;
  
  GaussMM<1> gmm;
  gmm.addComponent(Vector1d::Constant(1, mu), Vector1d::Constant(sd * sd), 1.0);
  for (double val = -8; val < 8; val += 0.2) {
    VectorXd vec = Vector1d::Constant(1, val);
    ASSERT_DOUBLE_EQ(gsl_ran_gaussian_pdf(val - mu, sd), gmm.likelihood(vec));
  }
}

TEST(GMM, Sample) {
  double mu = 5.0, sd = 2.0;
  
  GaussMM<1> gmm;
  gmm.addComponent(Vector1d::Constant(mu), Vector1d::Constant(sd * sd), 1.0);
  int nsamp = 5000;
  VectorXd samples = gmm.sample(nsamp);
  double mean = samples.sum() / nsamp;
  ASSERT_NEAR(mu, mean, 1e-1);
  for (int i = 0; i < nsamp; ++i) {
    samples(i) -= mean;
  }
  ASSERT_NEAR(sd, sqrt(samples.squaredNorm() / (nsamp - 1)), 1e-1);
}

TEST(GMM, Two2DNonstandardGaussian) {
  Vector2d mean1(-0.5, 0.3), mean2(2.2, 0.1);
  Vector2d cov1(4.3, 0.2), cov2(10.0, 5.3);
  double w1 = 0.2, w2 = 0.8;
  
  GaussMM<2> gmm;
  gmm.addComponent(mean1, cov1, w1);
  gmm.addComponent(mean2, cov2, w2);
  for (double val1 = -10; val1 < 10; val1 += 0.5) {
    for (double val2 = -10; val2 < 10; val2 += 0.5) {
      double c1 = w1 * (gsl_ran_gaussian_pdf(mean1(0) - val1, sqrt(cov1(0))) *
                        gsl_ran_gaussian_pdf(mean1(1) - val2, sqrt(cov1(1))));
      double c2 = w2 * (gsl_ran_gaussian_pdf(mean2(0) - val1, sqrt(cov2(0))) *
                        gsl_ran_gaussian_pdf(mean2(1) - val2, sqrt(cov2(1))));
      double expec = c1 + c2;
      
      VectorXd vec(2);
      vec << val1, val2;
      ASSERT_NEAR(expec, gmm.likelihood(vec), 1e-15);
    }
  }
}

TEST(GMM, Split1DStandard) {
  GaussMM<1> gmm;
  gmm.addComponent(Vector1d::Zero(), Vector1d::Constant(1), 1.0);
  gmm.splitComponent(0);
  for (double val = -3.0; val < 3.0; val += 0.5) {
    ASSERT_NEAR(gsl_ran_gaussian_pdf(val, 1.0), gmm.likelihood(Vector1d::Constant(val)), 8e-3);
  }
}

TEST(GMM, Split2DNonStd) {
  Vector2d mu(3, -8);
  Vector2d cov(5, 2);
  double weight = 0.5;
  GaussMM<2> gmm;
  gmm.addComponent(mu, cov, weight);
  gmm.splitComponent(0);
  for (double x = -10; x < 10; x += 0.4) {
    for (double y = -20; y < 0; y += 0.4) {
      Vector2d data(x, y);
      double expect = (weight / (2.0 * M_PI * sqrt(cov.prod())) *
                       exp(-0.5 * (data.cwiseProduct(data).cwiseProduct(cov.cwiseInverse()).sum())));
      ASSERT_NEAR(expect, gmm.likelihood(data), 0.1);
    }
  }
}

TEST(GMM, CombineTwo1DGaussians) {
  GaussMM<1> gauss1, gauss2;
  gauss1.addComponent(Vector1d::Constant(-7.3), Vector1d::Constant(0.1), 1.0);
  gauss2.addComponent(Vector1d::Constant(2.23), Vector1d::Constant(0.4), 1.0);

  std::vector<GaussMM<1>::Ptr> mixtures;
  mixtures.push_back(&gauss1);
  mixtures.push_back(&gauss2);
  
  GaussMM<2>::Ptr combined = GaussMM<1>::Combine<2>(mixtures);

  ASSERT_TRUE(NULL != combined);
  ASSERT_TRUE(combined->numComponents() > 0);
  EXPECT_TRUE(combined->numComponents() == 1);
  EXPECT_EQ(Vector2d(-7.3, 2.23), combined->means()[0]);
  EXPECT_EQ(Vector2d(0.1, 0.4), combined->covs()[0]);
  EXPECT_EQ(1.0, combined->weights()[0]);
  delete combined;
}

TEST(GMM, Combine1DGaussianAndGMM) {
  GaussMM<1> gmm, gauss;

  gauss.addComponent(Vector1d::Constant(2.23), Vector1d::Constant(0.1), 1.0);
  gmm.addComponent(Vector1d::Constant(0.0), Vector1d::Constant(4.3), 0.2);
  gmm.addComponent(Vector1d::Constant(-4.5), Vector1d::Constant(1.7), 0.8);

  std::vector<GaussMM<1>::Ptr> mixtures;
  mixtures.push_back(&gauss);
  mixtures.push_back(&gmm);
  
  GaussMM2::Ptr combined = GaussMM1::Combine<2>(mixtures);

  ASSERT_TRUE(NULL != combined);
  ASSERT_TRUE(combined->numComponents() > 0);
  EXPECT_TRUE(combined->numComponents() == 2);  
  EXPECT_EQ(Vector2d(2.23, 0.0), combined->means()[0]);
  EXPECT_EQ(Vector2d(0.1, 4.3), combined->covs()[0]);
  EXPECT_EQ(0.2, combined->weights()[0]);

  EXPECT_EQ(Vector2d(2.23, -4.5), combined->means()[1]);
  EXPECT_EQ(Vector2d(0.1, 1.7), combined->covs()[1]);
  EXPECT_EQ(0.8, combined->weights()[1]);
  delete combined;
}

TEST(GMM, CombineDynamicGMM) {
  GaussMM<1> gauss1, gauss2;
  gauss1.addComponent(Vector1d::Constant(2.0), Vector1d::Constant(1), 1.0);
  gauss2.addComponent(Vector1d::Constant(-2.0), Vector1d::Constant(1), 1.0);
  
  std::vector<GaussMM<1>::Ptr> mixtures;
  mixtures.push_back(&gauss1);
  mixtures.push_back(&gauss2);

  GaussMM<Eigen::Dynamic>::Ptr result = GaussMM<1>::Combine<Eigen::Dynamic>(mixtures);
  
  ASSERT_TRUE(NULL != result);
  EXPECT_TRUE(result->numComponents() == 1);
  EXPECT_EQ(Vector2d(2.0, -2.0), result->means(0));
  EXPECT_EQ(Vector2d(1.0, 1.0), result->covs(0));  
}

TEST(GMM, SumTwo2DGaussians) {
  GaussMM<2> gmm1, gmm2;
  gmm1.addComponent(Vector2d(0.0, -8.3), Vector2d(1.2, 4.3), 1.0);
  gmm2.addComponent(Vector2d(-7.2, 1.0), Vector2d(4.2, 7.1), 1.0);

  std::vector<GaussMM<2>::Ptr> vec;
  vec.push_back(&gmm1);
  vec.push_back(&gmm2);
  std::vector<double> weights;
  weights.push_back(0.3);
  weights.push_back(0.7);
  GaussMM<2>::Ptr sum = GaussMM<2>::Sum(vec, weights);

  ASSERT_TRUE(sum != NULL);
  
  EXPECT_EQ(gmm1.means(0), sum->means(0));
  EXPECT_EQ(gmm1.covs(0), sum->covs(0));
  EXPECT_EQ(gmm1.weights(0) * weights[0], sum->weights(0));
  
  EXPECT_EQ(gmm2.means(0), sum->means(1));
  EXPECT_EQ(gmm2.covs(0), sum->covs(1));
  EXPECT_EQ(gmm2.weights(0) * weights[1], sum->weights(1));
  delete sum;
}

TEST(GMM, SumMixtureModels) {
  using namespace gsl;
  GaussMM<3> gmm1, gmm2, gmm3;
  gmm1.addComponent(Vector3d::Random(), Vector3d(uniform(0, 10), uniform(0, 10), uniform(0, 1)), 0.2);
  gmm1.addComponent(Vector3d::Random(), Vector3d(uniform(0, 10), uniform(0, 10), uniform(0, 1)), 0.5);
  gmm1.addComponent(Vector3d::Random(), Vector3d(uniform(0, 10), uniform(0, 10), uniform(0, 1)), 0.3);  
  
  gmm2.addComponent(Vector3d::Random(), Vector3d(uniform(0, 10), uniform(0, 10), uniform(0, 1)), 1.0);

  gmm3.addComponent(Vector3d::Random(), Vector3d(uniform(0, 10), uniform(0, 10), uniform(0, 1)), 0.8);
  gmm3.addComponent(Vector3d::Random(), Vector3d(uniform(0, 10), uniform(0, 10), uniform(0, 1)), 0.2);  
  
  std::vector<GaussMM<3>::Ptr> vec;
  vec.push_back(&gmm1);
  vec.push_back(&gmm2);
  vec.push_back(&gmm3);
  std::vector<double> weights;
  weights.push_back(0.3);
  weights.push_back(0.6);
  weights.push_back(0.1);
  GaussMM<3>::Ptr sum = GaussMM<3>::Sum(vec, weights);

  ASSERT_TRUE(sum != NULL);

  int index = 0;
  for (int i = 0; i < 3; ++i) {
    for (size_t j = 0; j < vec[i]->numComponents(); ++j) {
      ASSERT_EQ(vec[i]->means(j), sum->means(index));
      ASSERT_EQ(vec[i]->covs(j), sum->covs(index));
      ASSERT_EQ(vec[i]->weights(j) * weights[i], sum->weights(index));
      ++index;
    }
  }
  
  delete sum;
}

//============================== Entropy tests ==============================//
class EntropyTest : public testing::Test {
protected:
  EntropyTest() {

  }

  virtual void SetUp() {    
    mm_1d_std_ = new GaussMM<1>();
    mm_1d_std_->addComponent(Vector1d::Zero(1), Matrix1d::Identity(1, 1), 1.0);

    mm_2d_std_ = new GaussMM<2>();
    mm_2d_std_->addComponent(Vector2d::Zero(2), Vector2d(1, 1), 1.0);

    cov_2d_nonstd = Vector2d(1.2, 5);
    mm_2d_nonstd_ = new GaussMM<2>();
    mm_2d_nonstd_->addComponent((VectorXd(2) << 49.2, -50.3).finished(),
                                cov_2d_nonstd,
                                1.0);
  }

  virtual void TearDown() {
    delete mm_1d_std_;
    delete mm_2d_std_;
    delete mm_2d_nonstd_;
  }
  // A non-standard covariance matrix whose entropy is 3.733757
  Vector2d cov_2d_nonstd;
  // 1d & 2d mixture models with single gaussian component 
  GaussMM<1>::Ptr mm_1d_std_;
  GaussMM<2>::Ptr mm_2d_std_, mm_2d_nonstd_;
};

TEST_F(EntropyTest, gaussEntropy_exact) {
  ASSERT_NEAR(1.418939, EntropyApproximator<1>::gaussEntropy(Vector1d::Ones()), 1e-6);
  ASSERT_NEAR(2.837877, EntropyApproximator<2>::gaussEntropy(Vector2d::Ones()), 1e-6);
  ASSERT_NEAR(3.733757, EntropyApproximator<2>::gaussEntropy(cov_2d_nonstd), 1e-6);
}

TEST_F(EntropyTest, numerical_gaussian) {
  // Test numerical entropy for single component mixture model
  NumericalApproximator<1> numer1(2000);
  EXPECT_DOUBLE_EQ(EntropyApproximator<1>::gaussEntropy(Vector1d::Ones()),
                   numer1.evaluate(*mm_1d_std_));
  NumericalApproximator<2> numer2(2000);  
  EXPECT_NEAR(EntropyApproximator<2>::gaussEntropy(Vector2d::Ones()),
              numer2.evaluate(*mm_2d_std_), 1e-2);
  NumericalApproximator<2> numer3(2000);  
  EXPECT_NEAR(EntropyApproximator<2>::gaussEntropy(cov_2d_nonstd),
              numer2.evaluate(*mm_2d_nonstd_), 1e-2);
}

TEST_F(EntropyTest, taylor_0_gaussian) {
  TaylorApproximator<1> taylor1(0, 0);
  ASSERT_DOUBLE_EQ(EntropyApproximator<1>::gaussEntropy(Vector1d::Ones()),
                   taylor1.evaluate(*mm_1d_std_) + 0.5);
  TaylorApproximator<2> taylor2(0, 0);  
  ASSERT_DOUBLE_EQ(EntropyApproximator<2>::gaussEntropy(Vector2d::Ones()),
                   taylor2.evaluate(*mm_2d_std_) + 1.0);
  ASSERT_DOUBLE_EQ(EntropyApproximator<2>::gaussEntropy(cov_2d_nonstd),
                   taylor2.evaluate(*mm_2d_nonstd_) + 1.0);
}

TEST_F(EntropyTest, taylor_2_gaussian) {
  TaylorApproximator<1> taylor1(2, 0);
  ASSERT_DOUBLE_EQ(EntropyApproximator<1>::gaussEntropy(Vector1d::Ones()),
                   taylor1.evaluate(*mm_1d_std_));
  TaylorApproximator<2> taylor2(2, 0);  
  ASSERT_DOUBLE_EQ(EntropyApproximator<2>::gaussEntropy(Vector2d::Ones()),
                   taylor2.evaluate(*mm_2d_std_));
  ASSERT_DOUBLE_EQ(EntropyApproximator<2>::gaussEntropy(cov_2d_nonstd),
                   taylor2.evaluate(*mm_2d_nonstd_));
}

//============================ Performance Tests ============================//

class PerformanceTest : public testing::Test {
protected:
  void SetUp() {
    gmm1d = new GaussMM<1>();
    int num_comp = 500;
    for (int i = 0; i < num_comp; ++i) {
      gmm1d->addComponent(Vector1d::Random(1),
                          (Vector1d() << gsl::uniform(0, 10)).finished(),
                          1.0 / num_comp);
    }

    gmm2d = new GaussMM<2>();
    for (int i = 0; i < num_comp; ++i) {
      gmm2d->addComponent(Vector2d::Random(2),
                          Vector2d(gsl::uniform(0, 10), gsl::uniform(0, 1)),
                          1.0 / num_comp);
    }    
  }

  void TearDown() {
    delete gmm1d;
    delete gmm2d;
  }
  
  GaussMM<1>::Ptr gmm1d;
  GaussMM<2>::Ptr gmm2d;
};

TEST_F(PerformanceTest, GMM1DLikelihood) {
  for (int i = 0; i < 1000; ++i) {
    gmm1d->likelihood(Vector1d::Constant((i - 500.0) / 1000.0));
  }
}

TEST_F(PerformanceTest, GMM2DLikelihood) {
  for (int i = 0; i < 1000; ++i) {
    gmm2d->likelihood(Vector2d::Constant((i - 500.0) / 1000.0));
  }
}

TEST_F(PerformanceTest, Taylor0_1DComponents) {
  TaylorApproximator<1> ta(0, 0);
  ta.evaluate(*gmm1d);
}

TEST_F(PerformanceTest, Taylor0_2DComponents) {
  TaylorApproximator<2> ta(0, 0);
  ta.evaluate(*gmm2d);
}

TEST_F(PerformanceTest, Taylor2_1DComponents) {
  TaylorApproximator<1> ta(2, 0);
  ta.evaluate(*gmm1d);
}

TEST_F(PerformanceTest, Taylor2_2DComponents) {
  TaylorApproximator<2> ta(2, 0);
  ta.evaluate(*gmm2d);
}

TEST_F(PerformanceTest, Numerical1DComponents) {
  NumericalApproximator<1> na(500);
  na.evaluate(*gmm1d);
}

TEST_F(PerformanceTest, Numerical2DComponents) {
  NumericalApproximator<2> na(500);
  na.evaluate(*gmm2d);
}

// Run all the tests that were declared with TEST()
int main(int argc, char **argv){
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
