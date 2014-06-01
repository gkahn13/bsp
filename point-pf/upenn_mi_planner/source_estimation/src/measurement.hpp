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

#ifndef MEASUREMENT_HPP
#define MEASUREMENT_HPP

#include <map>
#include <list>
#include <fstream>

#include <gmm/gmm.hpp>

#include <player_map/rosmap.hpp>

#include "types.hpp"

namespace rf {
  class Measurement {
  public:
    Measurement(double ts, int dest_id, int source_id, const Pose2D &source, double dist) :
      ts_(ts), dest_id_(dest_id), source_id_(source_id), dist_(dist), var_(0), dur_(0) {
      source_.setXYT(source);
    }
    
    Measurement(double ts, int dest_id, int source_id, const Pose2D &source, double dist,
                double variance, double duration) :
      ts_(ts), dest_id_(dest_id), source_id_(source_id), dist_(dist), var_(variance), dur_(duration) {
      source_.setXYT(source);
    }
  
    double timestamp() const { return ts_; }
    int destID() const { return dest_id_; }
    int sourceID() const { return source_id_; }
    double dist() const { return dist_; }
    double variance() const { return var_; }
    double duration() const { return dur_; }    
    const Pose2D & sourceLoc() const { return source_; }
  
    Measurement *Copy() const {
      Measurement *m = new Measurement(ts_, dest_id_, source_id_, source_,
                                       dist_, var_, dur_);
      return m;
    }
  
  private:
    double ts_;
    int dest_id_;
    int source_id_;
    double dist_, var_;
    double dur_;
    Pose2D source_;
    DISALLOW_COPY_AND_ASSIGN(Measurement);
  };

  std::ostream& operator<<(std::ostream &s, const Measurement &m);

  
  //======================== Measurement Aggregators ========================//
  class MeasurementAggregator {
  public:
    MeasurementAggregator() {}
    virtual ~MeasurementAggregator() {}
    virtual bool addMeasurement(const Measurement&) = 0;
    virtual const Measurement* getAggregate() const = 0;
  private:
    DISALLOW_COPY_AND_ASSIGN(MeasurementAggregator);
  };

  class NullAggregator : public MeasurementAggregator {
  public:
    NullAggregator();
    bool addMeasurement(const Measurement& m);
    const Measurement* getAggregate() const;
  private:
    const Measurement *msmt_;
    DISALLOW_COPY_AND_ASSIGN(NullAggregator);
  };

  // Aggregate distance measurements based on space and time.  
  //
  // When a new measurement is added, it is put into a queue.  When a new
  // measurement comes in that is 1) further than x-meters from the head and 2)
  // > t-seconds has elapsed since the measurement at thead head then an
  // aggregate measurement is returned.  This aggregate has the average of
  // distances in the queue, excluding the tail (which was the last
  // measurement).  
  //
  // This aggregator can be used to ensure that measurements are not added when
  // the robot is sitting still.
  // 
  // Measurements are grouped on a per sender & receiver pair.
  class STAggregator : public MeasurementAggregator {
  public:
    STAggregator(double secs, double dist);
    STAggregator(double secs, double dist, double max_age);

    ~STAggregator();

    bool addMeasurement(const Measurement& msmt);
    const Measurement* getAggregate() const;
    // Get original measurements that were used to make aggregate
    const std::vector<Measurement*>& getOrig() const;
  private:
    // Maps source ID to measurements that are pooled together
    std::map<std::pair<int, int>, std::list<Measurement*> *> source_msmts_;
    // The aggregated measurement to use for the measurement model
    Measurement *agg_msmt_;
    // The list of measurements used to create agg_msmt_
    std::vector<Measurement*> orig_msmts_;
    // Whether the aggregated measurement is valid
    bool valid_;
    // Cutoffs for determining how to group measurements
    double secs_, dist_;
    // Maximum age of measurement in the queue
    double max_age_;    
    DISALLOW_COPY_AND_ASSIGN(STAggregator);
  };  
  
  //============================== Base class ===============================//
  template<class P>
  class MeasurementModel {
  public:
    MeasurementModel() {};
    virtual ~MeasurementModel() {};
    virtual bool addMeasurement(const Measurement &m) = 0;
    virtual void applyMeasurement(P *p) = 0;
    virtual const Measurement* measurement() const = 0;
  private:
    DISALLOW_COPY_AND_ASSIGN(MeasurementModel<P>);
  };

  //============================ Gaussian Fading ============================//
  template<class P>
  class GaussFading : public MeasurementModel<P> {
  public:
    struct Params {
      double mu0, mu, sigma2;
    };
    
    GaussFading(const Params &pars, MeasurementAggregator *ma) :
      msmt_agg_(ma), params_(pars) {
    }
    
    ~GaussFading() { }
        
    bool addMeasurement(const Measurement &m) {
      return msmt_agg_->addMeasurement(m);
    }
        
    void applyMeasurement(P *particle) {
      const Measurement* msmt = measurement();
      double est_dist = sqrt(l2_squared(particle->point, msmt->sourceLoc()));

      GaussMM<1> gmm;
      gmm.addComponent(Vector1d::Constant(params_.mu0 +
                                          params_.mu * est_dist),
                       Vector1d::Constant(est_dist * params_.sigma2), 1.0);
      particle->weight *= gmm.likelihood(Vector1d::Constant(est_dist - msmt->dist()));
    }

    const Measurement* measurement() const {
      return msmt_agg_->getAggregate();
    }

    static GaussFading* FromROS(const ros::NodeHandle &nh, MeasurementAggregator *ma) {
      std::string type;
      nh.param<std::string>("measurement_model", type, "none");
      if (type != "gaussian") {
        ROS_ERROR("GaussFading::FromRos() Bad field 'measurement_model': %s",
                  type.c_str());
        ROS_BREAK();
      }
      Params p;
      nh.param("mu0", p.mu0, 2.1956);
      nh.param("mu", p.mu, -0.6286);
      nh.param("sigma2", p.sigma2, 1.9432);

      return new GaussFading(p, ma);
    }
    
  private:
    MeasurementAggregator *msmt_agg_;
    Params params_;
    DISALLOW_COPY_AND_ASSIGN(GaussFading<P>);
  };

  //====================== LOS / NLOS Gaussian Fading =======================//
  template<class P>
  class LosNlosFading : public MeasurementModel<P> {
  public:
    struct Params {
      struct RangeGauss {
        double mu0, mu, sigma2;
      } los, nlos;
      double close_mu, close_var;
      double range_cutoff;
    };
    
    LosNlosFading(const Params &ps, LOSChecker *checker,
                  MeasurementAggregator *ma) :
      checker_(checker), msmt_agg_(ma), params_(ps) {
    }

    virtual ~LosNlosFading() {};

    bool addMeasurement(const Measurement &m) {
      return msmt_agg_->getAggregate();
    }

    void applyMeasurement(P *particle) {
      // If weights are updated to frequently, the filter can converge incorrectly
      const Measurement *msmt = measurement();
      // if (msmt->dist() < 0) {
      //   ROS_WARN_THROTTLE(1, "LosNlosFading::setWeight() Measurement with negative distance");
      //   return;
      // }
      double est_dist = std::max(sqrt(l2_squared(particle->point, msmt->sourceLoc())),
                                 params_.range_cutoff);
      bool los = checker_->LineOfSight(msmt->sourceLoc().x(), msmt->sourceLoc().y(),
                                       particle->point.x, particle->point.y);
      particle->los = los;
      GaussMM<1> gmm;
      const typename LosNlosFading<P>::Params::RangeGauss *p;
      if (est_dist < 1) {
        gmm.addComponent(Vector1d::Constant(params_.close_mu),
                         Vector1d::Constant(params_.close_var),
                         1.0);
      } else if (los) {
        p = &(params_.los);
        gmm.addComponent(Vector1d::Constant(p->mu0 + p->mu * est_dist),
                         Vector1d::Constant(est_dist * p->sigma2),
                         1.0);
        // ROS_DEBUG_THROTTLE(10.0, "mu0: % 6.2f B: % 6.2f s2a: % 6.2f  est_dist: % 6.2f msmt_dist: % 6.2f new_weight = % 6.2f",
        //          p->mu0, p->mu, p->sigma2, est_dist, msmt->dist(), gmm.likelihood(Vector1d::Constant(est_dist - msmt->dist())));
      } else {
        p = &(params_.nlos);
        // ROS_DEBUG_THROTTLE(10.0, "mu0: % 6.2f B: % 6.2f s2a: % 6.2f  est_dist: % 6.2f msmt_dist: % 6.2f new_weight = % 6.2f",
        //          p->mu0, p->mu, p->sigma2, est_dist, msmt->dist(), gmm.likelihood(Vector1d::Constant(est_dist - msmt->dist())));        
        gmm.addComponent(Vector1d::Constant(p->mu0 + p->mu * est_dist),
                         Vector1d::Constant(est_dist * p->sigma2),
                         1.0);
        // ROS_INFO("mu0: % 6.2f B: % 6.2f s2a: % 6.2f  est_dist: % 6.2f msmt_dist: % 6.2f",
        //          p->mu0, p->B, p->sigma2alpha, est_dist, msmt->dist());
      }
      
      double new_weight = gmm.likelihood(Vector1d::Constant(est_dist - msmt->dist()));
      
      particle->weight *= new_weight;
    }
    
    const Measurement* measurement() const {
      return msmt_agg_->getAggregate();
    }

    static LosNlosFading<P>* FromROS(const ros::NodeHandle &nh,
                                     LOSChecker *checker,
                                     MeasurementAggregator *agg) {
      std::string type;
      nh.param<std::string>("measurement_model", type, "none");
      if (type != "fading") {
        ROS_ERROR("Wrong type");
        ROS_BREAK();
      }

      Params p;
      nh.param("los_mu0", p.los.mu0, 2.1956);
      nh.param("los_mu", p.los.mu, -0.6286);
      nh.param("los_sigma2", p.los.sigma2, 1.9432);
      nh.param("nlos_mu0", p.nlos.mu0, -1.3477);
      nh.param("nlos_mu", p.nlos.mu, -0.6445);
      nh.param("nlos_sigma2", p.nlos.sigma2, 3.0138);
      nh.param("close_mu", p.close_mu, -2.0);
      nh.param("close_var", p.close_var, 3.0);
      nh.param("range_cutoff", p.range_cutoff, 0.1);
      
      ROS_INFO("  range_cutoff: % 7.2f", p.range_cutoff);
      ROS_INFO("  close_mu: % 7.2f close_var: % 7.2f", p.close_mu, p.close_var);
      ROS_INFO("  los_mu0:  % 7.2f los_s2a:   % 7.2f los_mu:  % 7.2f",
               p.los.mu0, p.los.sigma2, p.los.mu);
      ROS_INFO("  nlos_mu0: % 7.2f nlos_s2a:  % 7.2f nlos_B: % 7.2f",
               p.nlos.mu0, p.nlos.sigma2, p.nlos.mu);
      return new LosNlosFading<P>(p, checker, agg);
    }

    const Params& getParams() const {
      return params_;
    }
	
  private:
    LOSChecker *checker_;
    MeasurementAggregator *msmt_agg_;
    Params params_;
  };
} 
#endif
