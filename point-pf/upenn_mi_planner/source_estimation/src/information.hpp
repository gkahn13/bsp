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

#ifndef INFORMATION_HPP
#define INFORMATION_HPP

#include <cstdio>

#include "ros/ros.h"
#include "ros/console.h"

#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_monte_miser.h>

#include <ros_gsl/random.hpp>

#include <player_map/rosmap.hpp>

#include <gmm/gmm.hpp>

#include <timeit/timeit.hpp>

#include "particles.hpp"
#include "planner.hpp"
#include "util.hpp"

// #define TIMING 

namespace rf {
  template <int D, typename T=rf::StationaryParticle>
  class MIEvaluator {
  public:
    struct Params {
      struct RangeGauss {
        double mu0, mu, sigma2;
      } los, nlos;
    };
    
    MIEvaluator(EntropyApproximator<D> *ea, LOSChecker *checker,
                const Params &p) :
      approx_(ea), checker_(checker), params_(p) {
      if (checker_ == NULL) {
        ROS_BREAK();
      }
    }
    ~MIEvaluator() {}

    // Evaluate MI for a single configuration.  Assumes that the
    // RobotConfiguration world coordinates are comparable to particle
    // coordinates.
    double evaluate(const RobotConfiguration &rc, const ParticleArray<T> &ps);

    // Return the configuration of robots that maximizes mutual information.
    // Assumes that the RobotConfiguration world coordinates are comparable to
    // particle coordinates.
    RobotConfiguration* findMax(const std::vector<RobotConfiguration *> &rcs,
                                const ParticleArray<T> &ps);


    static MIEvaluator<D, T>* FromROS(const ros::NodeHandle &nh,
                                      EntropyApproximator<D> *ea,
                                      LOSChecker *checker);
  private:
    EntropyApproximator<D> *approx_;
    LOSChecker *checker_;
    struct Params params_;
  };

  template <int D, typename T>
  RobotConfiguration* MIEvaluator<D, T>::findMax(const std::vector<RobotConfiguration *> &rcs,
                                                 const ParticleArray<T> &ps) {
    double max_mi = -1;
    RobotConfiguration *best_rc = NULL;
    for (size_t i = 0; i < rcs.size(); ++i) {
      RobotConfiguration *current = rcs[i];
      double mi = evaluate(*current, ps);
      ROS_DEBUG_STREAM("" << *current << mi);
      if (mi > max_mi) {
        best_rc = current;
        max_mi = mi;
      }
    }
    return best_rc;
  }

  template <int D, typename T>
  MIEvaluator<D, T>* MIEvaluator<D, T>::FromROS(const ros::NodeHandle &nh,
                                                EntropyApproximator<D> *ea,
                                                LOSChecker *checker) {
    Params p;
    nh.param("mieval_los_mu", p.los.mu, 0.0);
    nh.param("mieval_los_mu0", p.los.mu0, 0.0);
    nh.param("mieval_los_sigma2", p.los.sigma2, 1.0);
    nh.param("mieval_nlos_mu", p.nlos.mu, 0.0);
    nh.param("mieval_nlos_mu0", p.nlos.mu0, 0.0);
    nh.param("mieval_nlos_sigma2", p.nlos.sigma2, 1.0);
    ROS_INFO("Using MIEvaluator: ");
    ROS_INFO("  los:   mu=%6.2f mu0=%6.2f sigma2=%6.2f",
             p.los.mu, p.los.mu0, p.los.sigma2);
    ROS_INFO("  nlos:  mu=%6.2f mu0=%6.2f sigma2=%6.2f",
             p.nlos.mu, p.nlos.mu0, p.nlos.sigma2);
    return new MIEvaluator(ea, checker, p);
  }

  
  template <int D, typename T>
  double MIEvaluator<D, T>::evaluate(const RobotConfiguration &rc,
                                     const ParticleArray<T> &ps) {
    std::vector<typename GaussMM<D>::Ptr> mixtures(ps.size);
    std::vector<double> weights(ps.size);  
  
    // Build mixture for each particle
#ifdef TIMING
    timeit::tic();
#endif

    double cond_ent = 0.0;
    
    for(int pind = 0; pind < ps.size; ++pind) {
      const T *curr = ps.ps + pind;
      // Mixture of each 
      std::vector<GaussMM<1>::Ptr> bot_mixtures(rc.nbots);
    
      for (int bot = 0; bot < rc.nbots; ++bot) {

        Pose2D *bot_pose = &rc.poses[bot];
        bool los = checker_->LineOfSight(bot_pose->x(), bot_pose->y(),
                                         curr->point.x, curr->point.y);
        // Copy base GMM to new GMM and shift mean by current distance
        double distance = sqrt(l2_squared(curr->point, *bot_pose));
        bot_mixtures[bot] = new GaussMM<1>();

        struct Params::RangeGauss *p = NULL;
        
        if (distance < 1 || los) {
          p = &(params_.los);
        } else {
          p = &(params_.nlos);
        }
        bot_mixtures[bot]->addComponent(Vector1d::Constant(p->mu0 + p->mu * distance + distance),
                                        Vector1d::Constant(p->sigma2 * distance),
                                        1.0);
        
        // NOTE: This only works when the mixture model has one component
        // ROS_INFO_STREAM("cov: " << bot_mixtures[bot]->covs(0) << std::string(" ") << EntropyApproximator<1>::gaussEntropy(bot_mixtures[bot]->covs(0)));
          cond_ent += curr->weight * EntropyApproximator<1>::gaussEntropy(bot_mixtures[bot]->covs(0));
      }

      weights[pind] = curr->weight;
      mixtures[pind] = GaussMM<1>::Combine<D>(bot_mixtures);

      for (size_t i = 0; i < bot_mixtures.size(); ++i) {
        delete bot_mixtures[i];
      }
    }
#ifdef TIMING
    timeit::toc("Built mixture: ");
    timeit::tic();
#endif
    typename GaussMM<D>::Ptr combined_gmm = GaussMM<D>::Sum(mixtures, weights);

#ifdef TIMING
    std::stringstream ss;
    ss.str();
    ss << "Summed gmm (there are " << combined_gmm->numComponents() << " components): ";
    timeit::toc(ss.str().c_str());
    ROS_INFO("dim=%i, num components =%lu", combined_gmm->dim(), combined_gmm->numComponents());
    timeit::tic();
#endif
    double result = approx_->evaluate(*combined_gmm);
    result -= cond_ent;
    
#ifdef TIMING
    timeit::toc("Approximated: ");
#endif 
    delete combined_gmm;
    for (size_t i = 0; i < mixtures.size(); ++i) {
      delete mixtures[i];
    }

    return result;
  }
}
#endif
