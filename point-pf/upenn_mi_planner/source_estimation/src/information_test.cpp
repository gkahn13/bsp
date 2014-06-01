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

#include <Eigen/Core>

#include "ros_gsl/random.hpp"

#include "information.hpp"
#include "planner.hpp"
#include "particles.hpp"

using namespace Eigen;
using namespace rf;

ParticleArray<>* sampleCircle(int nparticles,
                              double x, double y, double r, double sd) {
  ParticleArray<> *pa = new ParticleArray<>();
  pa->resize(nparticles);
  
  for (int i = 0; i < nparticles; ++i) {
    double theta = gsl::uniform(0, 2 * M_PI);
    double rad = gsl::normal(sd) + r;
    pa->ps[i].point.x = x + rad * cos(theta);
    pa->ps[i].point.y = y + rad * sin(theta);
    pa->ps[i].weight = 1.0 / nparticles;
  }
  return pa;
}


ParticleArray<>* sampleGaussians(double x1, double y1,
                                 double x2, double y2,
                                 int nparticles = 500, double sd = 0.2) {
  ParticleArray<> *pa = new ParticleArray<>();
  pa->resize(nparticles);
  
  for (int i = 0; i < nparticles / 2; ++i) {
    pa->ps[i].point.x = x1 + gsl::normal(sd);
    pa->ps[i].point.y = y1 + gsl::normal(sd);
    pa->ps[i].weight = 1.0 / nparticles;
  }
  for (int i = nparticles / 2; i < nparticles; ++i) {
    pa->ps[i].point.x = x2 + gsl::normal(sd);
    pa->ps[i].point.y = y2 + gsl::normal(sd);
    pa->ps[i].weight = 1.0 / nparticles;        
  }
  return pa;
}

ParticleArray<>* singleGaussian(double x, double y,
                                int nparticles = 500, double sd = 0.2) {
  ParticleArray<> *pa = new ParticleArray<>();
  pa->resize(nparticles);
  
  for (int i = 0; i < nparticles; ++i) {
    pa->ps[i].point.x = x + gsl::normal(sd);
    pa->ps[i].point.y = y + gsl::normal(sd);
    pa->ps[i].weight = 1.0 / nparticles;
  }
  return pa;
}


class InformationTest : public testing::Test {
protected:
    InformationTest() : up_(1), down_(1), left_(1), right_(1), origin_(1), checker_(NULL) {

    }
    virtual void SetUp() {
        up_.poses[0].setXYT(0, 1, 0);
        down_.poses[0].setXYT(0, -1, 0);
        right_.poses[0].setXYT(1, 0, 0);
        left_.poses[0].setXYT(-1, 0, 0);
        mi_params_.los.mu = 0;
        mi_params_.los.mu0 = 0;
        mi_params_.los.sigma2 = 1;
        mi_params_.nlos.mu = 0;
        mi_params_.nlos.mu0 = 0;
        mi_params_.nlos.sigma2 = 1;
    }

    RobotConfiguration up_, down_, left_, right_, origin_;
    LOSChecker checker_;
  
    MIEvaluator<1>::Params mi_params_;
    static int nparticles;
};

int InformationTest::nparticles = 500;

// Consider movement when one robot is midway between two distinct gaussians.
TEST_F(InformationTest, TwoSymmetricGaussiansMI) {
    double x1 = 0, y1 = 5;
    double x2 = 0, y2 = -5;

    ParticleArray<> *pa = sampleGaussians(x1, y1, x2, y2);
    
    // NumericalApproximator na(100);
    TaylorApproximator<1> ta(0, 0);
    MIEvaluator<1> eval(&ta, &checker_, mi_params_);

    double mi_up = eval.evaluate(up_, *pa);
    double mi_down = eval.evaluate(down_, *pa);
    double mi_right = eval.evaluate(right_, *pa);
    double mi_left = eval.evaluate(left_, *pa);

    // up/down & left/right pairs have similar rewards
    EXPECT_LT(abs(mi_up - mi_down), 1e-14);
    EXPECT_LT(abs(mi_left - mi_right), 1e-14);
    // up/down direction is better than left/right
    EXPECT_GT(mi_up - mi_left, 1e-4);
    
    delete pa;
}

TEST_F(InformationTest, TwoSymmetricGaussiansFindMax) {
    double x1 = 0, y1 = 4;
    double x2 = 0, y2 = -4;

    ParticleArray<> *pa = sampleGaussians(x1, y1, x2, y2);

    NumericalApproximator<1> na(10);
    MIEvaluator<1> eval(&na, &checker_, mi_params_);

    std::vector<RobotConfiguration *> configs;
    configs.push_back(&left_);
    
    RobotConfiguration *actual = eval.findMax(configs, *pa);
    EXPECT_TRUE(left_.samePosition(*actual));

    configs.push_back(&up_);
    actual = eval.findMax(configs, *pa);
    EXPECT_TRUE(up_.samePosition(*actual));

    configs.push_back(&left_);
    configs.push_back(&down_);
    actual = eval.findMax(configs, *pa);
    EXPECT_TRUE(up_.samePosition(*actual) || down_.samePosition(*actual));
    
    delete pa;
}

// Have two gaussians, and robot closer to one of them.
TEST_F(InformationTest, TwoAsymmetricGaussians) {
    double x1 = 2, y1 = 5;
    double x2 = 2, y2 = -2;

    ParticleArray<> *pa = sampleGaussians(x1, y1, x2, y2);

    NumericalApproximator<1> na(50);
    MIEvaluator<1> eval(&na, &checker_, mi_params_);    

    double mi_up = eval.evaluate(up_, *pa);
    double mi_down = eval.evaluate(down_, *pa);
    double mi_right = eval.evaluate(right_, *pa);
    double mi_left = eval.evaluate(left_, *pa);

    // Left/right are similar
    EXPECT_LT(abs(mi_left - mi_right), 1e-14);
    // Down is best, left/right is better than up
    EXPECT_GT(mi_down - mi_up, 1e-4);
    EXPECT_GT(mi_down - mi_right, 1e-4);
    EXPECT_GT(mi_left, mi_up);

    delete pa;
}

#define NBOTS 3
class PerformanceTest : public testing::Test {
protected:
  PerformanceTest() : checker_(NULL) {};
  
  virtual void SetUp() {
    double x1 = 2, y1 = 5;
    double x2 = 2, y2 = -2;
    pa_ = sampleGaussians(x1, y1, x2, y2);
    
    rcs_.resize(numConfigs);
    for (int i = 0; i < numConfigs; ++i) {
      rcs_[i] = new RobotConfiguration(NBOTS);
      for (int j = 0; j < NBOTS; ++j) {
        rcs_[i]->poses[j].setXYT(gsl::uniform(-20, 20), gsl::uniform(-20, 20), gsl::uniform(-20, 20));
      }
    }

    mi_params_.los.mu = 0;
    mi_params_.los.mu0 = 0;
    mi_params_.los.sigma2 = 5;
    mi_params_.nlos.mu = 0;
    mi_params_.nlos.mu0 = 0;
    mi_params_.nlos.sigma2 = 5;    
  }

  virtual void TearDown() {
    delete pa_;
    for (int i = 0; i < numConfigs; ++i) {
      delete rcs_[i];
    }
  }
  
  LOSChecker checker_;
  ParticleArray<> *pa_;
  std::vector<RobotConfiguration *> rcs_;
  static const int nparticles;
  static const int numConfigs;
  struct MIEvaluator<NBOTS>::Params mi_params_;  
};

const int PerformanceTest::nparticles = 1000;
const int PerformanceTest::numConfigs = 20;

TEST_F(PerformanceTest, TwoRobotsTwoAssymetricParticlesTaylor2) {
    TaylorApproximator<NBOTS> ta(2, 0);
    MIEvaluator<NBOTS> eval(&ta, &checker_, mi_params_);
    eval.findMax(rcs_, *pa_);
}

TEST_F(PerformanceTest, TwoRobotsTwoAssymetricParticlesTaylor0) {
    TaylorApproximator<NBOTS> ta(0, 0);
    MIEvaluator<NBOTS> eval(&ta, &checker_, mi_params_);
    eval.findMax(rcs_, *pa_);
}

TEST_F(PerformanceTest, TwoRobotsTwoAssymetricParticlesNumerical) {
    NumericalApproximator<NBOTS> na(1000);
    MIEvaluator<NBOTS> eval(&na, &checker_, mi_params_);
    eval.findMax(rcs_, *pa_);
}
