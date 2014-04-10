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

#ifndef PARTICLES_HPP
#define PARTICLES_HPP

#include "ros/ros.h"
#include "ros/console.h"

#include <vector>
#include <cmath>
#include <stdexcept>
#include <sstream>
#include <fstream>
#include <ios>
#include <cstring>
#include <iostream>
#include <map>
#include <list>
#include <cassert>

#include <Eigen/Dense>

#include <gsl/gsl_cdf.h>

#include <geometry_msgs/PoseArray.h>

#include <player_map/rosmap.hpp>

#include "types.hpp"

#include "measurement.hpp"

namespace rf {
  //============================ Construction =============================//
  template<class P>
  class ParticleBuilder {
  public:
    ParticleBuilder() : rsd_(10.0) {};

    virtual ~ParticleBuilder() {};

    void setRadSD(double rsd) {
      rsd_ = rsd;
    }
    
    ParticleArray<P>* NewParticles(int num, const Measurement &m);


  private:
    void repInv();
    double rsd_;
    DISALLOW_COPY_AND_ASSIGN(ParticleBuilder<P>);
  };

  //============================ Motion models ============================//
  template<class P>
  class MotionModel {
  public:
    MotionModel() {};
    virtual ~MotionModel() {} ;

    virtual void setMotion(const Pose2D &previous, const Pose2D &current, double duration) = 0;
    virtual void applyMotion(P *particle) = 0;
  private:
    DISALLOW_COPY_AND_ASSIGN(MotionModel<P>);
  };

  // Apply random noise to x, y component of particle
  class RandomXYMotion : public MotionModel<StationaryParticle> {
  public:
    RandomXYMotion(double sd_x, double sd_y, const gsl_rng *rng) :
      sd_x_(sd_x), sd_y_(sd_y), rng_(rng) {
    }

    void setMotion(const Pose2D&, const Pose2D &, double duration) {
    }

    void applyMotion(StationaryParticle *p) {
      p->point.x += gsl_ran_gaussian(rng_, sd_x_);
      p->point.y += gsl_ran_gaussian(rng_, sd_y_);
    }

  private:
    DISALLOW_COPY_AND_ASSIGN(RandomXYMotion);

    double sd_x_, sd_y_;
    const gsl_rng *rng_;
  };

  // Apply random noise to x, y component of particle
  class RandomRTMotion : public MotionModel<StationaryParticle> {
  public:
    RandomRTMotion(double r_sd_vel, double theta_sd_vel, double max_travel_dist, LOSChecker *checker) :
      r_sd_vel_(r_sd_vel), theta_sd_vel_(theta_sd_vel), checker_(checker),
      max_travel_dist_(max_travel_dist) {
    }

    void setMotion(const Pose2D& source_pose, const Pose2D &, double duration) {
      source_pose_.setXYT(source_pose);
      if (duration < 1e-6) {
        duration = 0.25;
      }
      dur_ = duration;
    }

    void applyMotion(StationaryParticle *p) {
      // Get polar representation of current particle using source as origin,
      // perturb state along theta and r
      double dx = p->point.x - source_pose_.x();
      double dy = p->point.y - source_pose_.y();
      double theta = atan2(dy, dx);
      double r = hypot(dx, dy);

      bool los;
      if (checker_ != NULL) {
        los = checker_->LineOfSight(p->point.x, p->point.y, source_pose_.x(), source_pose_.y());
      }

      double r_sd = r_sd_vel_ * sqrt(dur_), theta_sd = theta_sd_vel_ * sqrt(dur_);

      double new_r, new_theta, new_x, new_y, travel_dist;
      do {
        new_r = r + gsl::normal(r_sd);
        new_theta = theta + gsl::normal(theta_sd);
        new_x = new_r * cos(new_theta) + source_pose_.x();
        new_y = new_r * sin(new_theta) + source_pose_.y();
        r_sd /= 2.0;
        theta_sd /= 2.0;

        travel_dist = sqrt(l2_squared(Point2D(new_x, new_y),
                                      Point2D(p->point.x, p->point.y)));
      } while (travel_dist > max_travel_dist_ ||
               (checker_ != NULL &&
                los != checker_->LineOfSight(new_x, new_y,
                                             source_pose_.x(), source_pose_.y())));

      p->point.x = new_x;
      p->point.y = new_y;
    }

  private:
    DISALLOW_COPY_AND_ASSIGN(RandomRTMotion);
    Pose2D source_pose_;
    double r_sd_vel_, theta_sd_vel_;
    double dur_;
    LOSChecker *checker_;
    double max_travel_dist_;
  };

  //============================= Resampling ==============================//
  template<class P>
  class Resampler {
  public:
    virtual ~Resampler() {};

    virtual bool resample(ParticleArray<P>* ps, int updates) = 0;
  };

  template<class P>
  class PeriodicResampler : public Resampler<P> {
  public:
    PeriodicResampler(int period) {
      if (period < 1) {
        throw std::runtime_error("Period must be positive integer");
      }
      period_ = period;
    }

    bool resample(ParticleArray<P>* pa, int updates) {
      if (updates < period_) {
        return false;
      }

      double *cdf;
      double sum;
      pa->BuilldCDF(&cdf, &sum);

      // Copy particles to new location, so that current group of
      // particles can be mutated
      ParticleArray<P> *old_pa = pa->Copy();

      // Draw sample uniformly at random, find its value in CDF, and
      // copy over particle to ps from old_ps
      for (int i = 0; i < pa->size; ++i) {
        double sample = gsl::uniform(0, sum);
        int ind = findFirstGreater(cdf, pa->size, sample);
        // std::cout << "sample = " << sample << " ind = " << ind << std::endl;

        if (ind < 0 || ind >= pa->size) {
          std::cerr << "Bad index" << std::endl;
          delete old_pa;
          exit(-1);
        }

        pa->ps[i] = old_pa->ps[ind];
        pa->ps[i].weight = 1.0 / pa->size;
      }

      // std::cout << "Old: " << std::endl;
      // for (int i = 0; i < num_particles; ++i) {
      //     std::cout << old_ps[i] << std::endl;
      // }

      // std::cout << "New: " << std::endl;
      // for (int i = 0; i < num_particles; ++i) {
      //     std::cout << ps[i] << std::endl;
      // }
      delete old_pa;
      delete[] cdf;

      return true;
    }
  private:

    int period_;
  };


  // Implements Low_variance_sampler() from Probabilistic Robotics (PR), Section
  // 4.3, p. 110. Resamples only when neff is low
  template <typename P>
  class LowVarianceSampler : public Resampler<P> {
  public:
    LowVarianceSampler(double neff) {
      if (neff < 0 || neff > 1) {
        throw std::runtime_error("neff must be between 0 and 1");
      }
      neff_pct_ = neff;
    }

    bool resample(ParticleArray<P> *pa, int /* updates */) {
      // ROS_ERROR("sum: %f sq_sum: %f neff: %f", sum, sq_sum, neff);
      double neff = pa->neff();
      if (neff / pa->size >= neff_pct_) {
        return false;
      }

      // Copy original particle set to new storage and operate entirely
      // on that
      ParticleArray<P> *old_pa = pa->Copy();

      double normalizer = 0;
      for (int i = 0; i < pa->size; ++i) {
        normalizer += pa->ps[i].weight;
      }
      // DEBUG(normalizer);

      // "r" in PR
      double increment = gsl::uniform(0, 1.0 / pa->size);
      // "c" in PR
      double cdf = old_pa->ps[0].weight / normalizer;

      int i = 0;
      for (int particle_ind = 0; particle_ind < pa->size; ++particle_ind) {
        double U = increment + (particle_ind) * (1.0 / pa->size);
        while (U > cdf) {
          ++i;
          cdf += old_pa->ps[i].weight / normalizer;
        }
        assert(i < pa->size);
        pa->ps[particle_ind] = old_pa->ps[i];
        pa->ps[particle_ind].weight = 1.0 / pa->size;
      }

      // std::cout << "Old: " << std::endl;
      // for (int i = 0; i < num_particles; ++i) {
      //     std::cout << old_ps[i] << std::endl;
      // }

      // std::cout << "New: " << std::endl;
      // for (int i = 0; i < num_particles; ++i) {
      //     std::cout << ps[i] << std::endl;
      // }

      delete old_pa;
      return true;
    }
  private:
    double neff_pct_;
  };
  
  //=============================== Filter ================================//
  template<class P>
  class ParticleFilter {
  public:
    ParticleFilter(MotionModel<P> *mot, MeasurementModel<P> *msmt,
               Resampler<P> *r, ParticleBuilder<P> *pb) :
      particles_(NULL), weight_updates_(0), motion_(mot),
      measurement_(msmt), resampler_(r), builder_(pb) {
    }

    ~ParticleFilter() {
      delete particles_;
    }

    virtual bool addMeasurementAndMotion(const Measurement &m,
                                         const Pose2D &previous,
                                         const Pose2D &current) = 0;
    
    const Measurement* getMeasurement() const {
      return measurement_->measurement();
    }

    void initialize(int num_particles, const Measurement &measurement) {
      delete particles_;

      particles_ = builder_->NewParticles(num_particles, measurement);

      for (int i = 0; i < num_particles; ++i) {
        particles_->ps[i].weight = 1.0 / num_particles;
      }
    }

    const ParticleArray<P>* getParticles() {
      return particles_;
    }

  protected:
    DISALLOW_COPY_AND_ASSIGN(ParticleFilter<P>);
    ParticleArray<P> *particles_;

    int weight_updates_;

    MotionModel<P> *motion_;
    MeasurementModel<P> *measurement_;
    Resampler<P> *resampler_;
    ParticleBuilder<P> *builder_;
  };

  template<class P>
  class StandardPF : public ParticleFilter<P> {
  public:
    StandardPF(MotionModel<P> *mot, MeasurementModel<P> *msmt, Resampler<P> *r,
               ParticleBuilder<P> *pb) : ParticleFilter<P>(mot, msmt, r, pb) {
      ;
    }
    
    bool addMeasurementAndMotion(const Measurement &m,
                                 const Pose2D &previous, const Pose2D &current) {
      bool action = this->measurement_->addMeasurement(m);
      if (action) {
        ++this->weight_updates_;
        for (int i = 0; i < this->particles_->size; ++i) {
          this->measurement_->applyMeasurement(&(this->particles_->ps[i]));
        }

        bool resampled = this->resampler_->resample(this->particles_,
                                                    this->weight_updates_);
        if (resampled) {
          this->weight_updates_ = 0;
        }

        this->particles_->normalize();
        
        this->motion_->setMotion(previous, current,
                                 this->measurement_->measurement()->duration());
        for (int i = 0; i < this->particles_->size; ++i) {
          this->motion_->applyMotion(&(this->particles_->ps[i]));
        }
      }
      return action;
    }
  };

  //============================ Visualization ============================//
  template<typename P>
  class ExperimentSaver {
  public:
    ExperimentSaver(std::string fname) {
      out_ = new std::ofstream(fname.c_str(), std::ios::binary | std::ios::out);

    }

    ~ExperimentSaver() {
      out_->close();
      delete out_;
    }

    void add(const Measurement &m,
             const ParticleArray<P> *pa) {
      // Write measurement
      double ts = m.timestamp(), dist = m.dist();
      int32_t sourceId = m.sourceID(), destId = m.destID();
      double sourcex = m.sourceLoc().x(), sourcey = m.sourceLoc().y();
      double var = m.variance(), duration = m.duration();
      out_->write((char*) &ts, sizeof(ts));
      out_->write((char*) &dist, sizeof(dist));

      out_->write((char*) &var, sizeof(var));
      out_->write((char*) &duration, sizeof(duration));

      out_->write((char*) &destId, sizeof(destId));
      out_->write((char*) &sourceId, sizeof(sourceId));
      out_->write((char*) &sourcex, sizeof(sourcex));
      out_->write((char*) &sourcey, sizeof(sourcey));

      // Write particles
      int32_t diff = pa->size;
      out_->write((char*) &diff, sizeof(diff));
      for (int i = 0; i < diff; ++i) {
        double x = pa->ps[i].point.x, y = pa->ps[i].point.y, w = pa->ps[i].weight;
        out_->write((char*) &x, sizeof(x));
        out_->write((char*) &y, sizeof(y));
        out_->write((char*) &w, sizeof(w));
      }
    }

  private:
    std::ofstream *out_;
  };

  template<class P>
  class AdaptivePF : public ParticleFilter<P> {
  public:
    AdaptivePF(MotionModel<P> *mot, MeasurementModel<P> *msmt, Resampler<P> *r,
               ParticleBuilder<P> *pb) : ParticleFilter<P>(mot, msmt, r, pb) {
      ;
    }
    
    bool addMeasurementAndMotion(const Measurement &m,
                                 const Pose2D &previous, const Pose2D &current) {
      bool action = this->measurement_->addMeasurement(m);
      if (action) {
        this->motion_->setMotion(previous, current,
                                 this->measurement_->measurement()->duration());

        this->resampler_->resample(this->particles_, this->weight_updates_);

        this->particles_->normalize();
        
      }
      return action;
    }
  };
    
}
#endif
