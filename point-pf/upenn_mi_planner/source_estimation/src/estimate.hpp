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

#ifndef ESTIMATE_HPP
#define ESTIMATE_HPP

#include <iostream>
#include <fstream>
#include <sys/time.h>

#include <gsl/gsl_rng.h>

#include <Eigen/Core>

#include "ros/ros.h"
#include "ros/console.h"

#include "tf/tf.h"

#include "geometry_msgs/PoseArray.h"
#include "geometry_msgs/PoseWithCovarianceStamped.h"
#include "nav_msgs/GetMap.h"
#include "visualization_msgs/MarkerArray.h"

#include "nanotron/Range.h"
#include "nanotron/RangeWithPose.h"

#include <ros_gsl/random.hpp>

#include <player_map/rosmap.hpp>

#include <gmm/gmm.hpp>

#include "util.hpp"
#include "particles.hpp"

namespace rf { 
  // A node that estimates the locaiton of a static source
  class StaticSourceNode {
  public:
    StaticSourceNode();

    void measurement_callback(const nanotron::RangeWithPose::ConstPtr &msg);
    void updateMeasModel(const Measurement &msmt);
    // Publish new particles
    void pubParticles();

    geometry_msgs::PoseArray getParticles();
    bool updated();

    ~StaticSourceNode();
    void BuildFilter(map_t *map);
  
  private:
    ros::NodeHandle node_;
    ros::NodeHandle private_nh_;
    ros::Subscriber measurement_sub_;
    ros::Publisher particle_pub_, point_cloud_pub_;
    int target_id_; 
    ParticleFilter<StationaryParticle> *pf_;
    MeasurementModel<StationaryParticle> *meas_mod_;
    Resampler<StationaryParticle> *resampler_;
    MotionModel<StationaryParticle> *mot_mod_;
    map_t *map_;  
    ParticleBuilder<StationaryParticle> *builder_;
    MeasurementAggregator *agg_;
    LOSChecker *checker_;
    ExperimentSaver<StationaryParticle> *saver_;
    int save_every_, save_count_;
    // source id -> information about robots making measurements
    std::map<int, Point2D> bot_poses_;
    std::map<int, std::string> bot_names_;
    std::string logfile_;
    bool updated_;
    double nparticles_;
    bool pf_init_;
  };
}
#endif
