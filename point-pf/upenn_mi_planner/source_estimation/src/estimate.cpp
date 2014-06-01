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

#include "estimate.hpp"

namespace rf {
  // A node that estimates the locaiton of a static source
  StaticSourceNode::StaticSourceNode() :
    private_nh_("~"), pf_(NULL), meas_mod_(NULL), resampler_(NULL),
    mot_mod_(NULL), map_(NULL), builder_(NULL), checker_(NULL), saver_(NULL),
    save_every_(0), save_count_(0), updated_(false) {
    map_ = requestCSpaceMap("map_valid");
    // Initialize the particle filter, including measurement paramters
    BuildFilter(map_);
    particle_pub_ = node_.advertise<geometry_msgs::PoseArray>("source_estimate", 2, true);
    point_cloud_pub_ = node_.advertise<visualization_msgs::Marker>("/visualization_msgs", 2, true);
    measurement_sub_ = node_.subscribe("range_pose", 500,
                                       &StaticSourceNode::measurement_callback, this);

    private_nh_.getParam("target_id", target_id_);
    ROS_INFO("Tracking target %i", target_id_);
    
    if (private_nh_.hasParam("logfile")) {
      private_nh_.getParam("logfile", logfile_);
      saver_ = new ExperimentSaver<StationaryParticle>(logfile_);
      ROS_WARN("Saving log to %s", logfile_.c_str());
      private_nh_.param("save_every", save_every_, 50);
    } else {
      ROS_WARN("Not saving log");
    }

  }

  void StaticSourceNode::measurement_callback(const nanotron::RangeWithPose::ConstPtr &msg) {
    // Because this is the only callback, and the main loop doesn't do
    // anything, we don't need a mutex to guard anything.

    // Update distance
    std::map<int, std::string>::iterator f = bot_names_.find(msg->range.src_id);
    if (f == bot_names_.end()) {
      bot_names_[msg->range.src_id] = msg->id;
    }
    if (bot_names_[msg->range.src_id] != msg->id) {
      ROS_ERROR("Bot with ID=%i changed its name from %s to %s",
                msg->range.src_id, bot_names_[msg->range.src_id].c_str(), msg->id.c_str());
    }
        
    // Incorporate new measurement
    bot_poses_[msg->range.src_id] = Point2D(msg->pose.pose.pose.position.x,
                                            msg->pose.pose.pose.position.y);
    Pose2D source_pose;
    source_pose.setXYT(msg->pose.pose.pose.position.x, msg->pose.pose.pose.position.y, 0);    
    
    Measurement m(msg->range.header.stamp.toSec(), msg->range.dest_id, msg->range.src_id,
                  source_pose, msg->range.distance);

    if (msg->range.dest_id != target_id_) {
      return;
    }

    // Check if particle filter has been initialized yet
    if (!pf_init_) {
      pf_init_ = true;
      pf_->initialize(nparticles_, m);
      pubParticles();
      return;
    }

    if (!pf_->addMeasurementAndMotion(m, source_pose, source_pose)) {
      return;
    };
    
    if (saver_ != NULL) {
      if (save_count_ == 0) {
        saver_->add(*pf_->getMeasurement(), pf_->getParticles());
      }
      save_count_ = (save_count_ + 1) % save_every_;
    }

    pubParticles();
  }
  
  // Publish new particles
  void StaticSourceNode::pubParticles() {
    geometry_msgs::PoseArray parts = getParticles();
    particle_pub_.publish(parts);

    visualization_msgs::Marker m;
    m.header.stamp = ros::Time();
    m.header.frame_id = "/map";
    m.pose.orientation.w = 1.0;
    // m.action = visualization_msgs::Marker::ADD;
    // m.type = visualization_msgs::Marker::CUBE_LIST;
    // m.id = 0;
    // m.ns = "map";
    // m.pose.orientation.w = 1;
    // m.scale.x = 0.2;
    // m.scale.y = 0.2;
    // m.scale.z = 0.01;
    // m.color.a = 1.0;
    // m.color.r = 1.0;
    // m.points.resize(parts.poses.size());
    // m.colors.resize(parts.poses.size());

    const ParticleArray<StationaryParticle> *pa = pf_->getParticles();

    // double max_weight = 0.0;
    // for (int i =0; i < pa->size; ++i) {
    //   max_weight = std::max(pa->ps[i].weight, max_weight);
    // }

    double mean_x = 0.0, mean_y = 0.0;
    for (int i =0; i < pa->size; ++i) {
      // m.points[i].x = pa->ps[i].point.x;
      // m.points[i].y = pa->ps[i].point.y;

      mean_x += pa->ps[i].point.x * pa->ps[i].weight;
      mean_y += pa->ps[i].point.y * pa->ps[i].weight;
	}
    //   double scale = pa->ps[i].weight / max_weight;
    //   m.colors[i].a = 1.0;
    //   m.colors[i].r = 1.0 * scale;
    //   m.colors[i].g = 1.0 * scale;
    //   m.colors[i].b = 1.0 * scale;
    // }
    // point_cloud_pub_.publish(m);

    // Publish mean
    m.action = visualization_msgs::Marker::ADD;
    m.type = visualization_msgs::Marker::SPHERE;
    m.id = 1;
    m.ns = "map";
    m.pose.orientation.w = 1;
    m.scale.x = 0.4;
    m.scale.y = 0.4;
    m.scale.z = 0.4;
    m.color.a = 1.0;
    m.color.r = 1.0;
    m.pose.position.x = mean_x;
    m.pose.position.y = mean_y;
    m.pose.orientation.w = 1.0;
    point_cloud_pub_.publish(m);
  }

  geometry_msgs::PoseArray StaticSourceNode::getParticles() {
    const ParticleArray<StationaryParticle> *pa = pf_->getParticles();    
    geometry_msgs::PoseArray particle_msg;
    particle_msg.poses.resize(pa->size);
    particle_msg.header.stamp = ros::Time::now();
    particle_msg.header.frame_id = "/map"; // TODO: make this a ros param

    particle_msg.poses.resize(pa->size);

    for (int ind = 0; ind < pa->size; ++ind) {
      tf::poseTFToMsg(tf::Pose(tf::Quaternion::getIdentity(),
                               tf::Vector3(pa->ps[ind].point.x,
                                           pa->ps[ind].point.y,
                                           pa->ps[ind].weight)),
                      particle_msg.poses[ind]);

    }
    return particle_msg;
  }

  bool StaticSourceNode::updated() {
    return updated_;
  }

  StaticSourceNode::~StaticSourceNode() {
    std::cerr << "Saved log to " << logfile_ << std::endl;
    delete pf_;
    delete meas_mod_;
    delete resampler_;
    delete mot_mod_;
    delete builder_;
    delete checker_;
    delete saver_;
    map_free(map_);
  }

  void StaticSourceNode::BuildFilter(map_t *map) {
    //====================== Construct particle filter ======================//
    gsl_rng *rng = gsl::rng();
    double neff, init_range_sd;
    private_nh_.param("nparticles", nparticles_, 1000.0);
    private_nh_.param("neff", neff, 0.3);
    private_nh_.param("init_range_sd", init_range_sd, 15.0);
    
    double st_min_secs, st_min_dist, st_max_age;
    private_nh_.param("st_max_secs", st_min_secs, 0.0);
    private_nh_.param("st_min_dist", st_min_dist, 0.0);
    private_nh_.param("st_max_age", st_max_age, 30.0);
    agg_ = new STAggregator(st_min_secs, st_min_dist, st_max_age);

    checker_ = new LOSChecker(map);
    
    meas_mod_ = getMeasurementModel<StationaryParticle>(private_nh_, agg_, checker_);
    
    std::string mot_model_type;
    private_nh_.param("mot_model_type", mot_model_type, std::string("rt"));
    if (mot_model_type == "xy") {
      double x_sd, y_sd;
      ROS_INFO("Using random xy motion (cartesian)");
      private_nh_.getParam("x_sd", x_sd);
      private_nh_.getParam("y_sd", y_sd);
      mot_mod_ = new RandomXYMotion(x_sd, y_sd, rng);
    } else if (mot_model_type == "rt") {
      double r_sd, theta_sd, max_travel_dist;
      private_nh_.getParam("r_sd_vel", r_sd);
      private_nh_.getParam("theta_sd_vel", theta_sd);
      private_nh_.getParam("max_travel_dist", max_travel_dist);
      mot_mod_ = new RandomRTMotion(r_sd, theta_sd, max_travel_dist, checker_);
      ROS_INFO("Using random rt motion w/ map");
    } else if (mot_model_type == "rt_nomap") {
      double r_sd, theta_sd, max_travel_dist;
      private_nh_.getParam("r_sd_vel", r_sd);
      private_nh_.getParam("theta_sd_vel", theta_sd);
      private_nh_.getParam("max_travel_dist", max_travel_dist);
      mot_mod_ = new RandomRTMotion(r_sd, theta_sd, max_travel_dist, NULL);
      ROS_INFO("Using random rt motion w/out map");
    } else {
      ROS_ERROR("Unknown motion model type");
      ROS_BREAK();
    }

    builder_ = new ParticleBuilder<StationaryParticle>();
    builder_->setRadSD(init_range_sd);

    
    std::string resample_type;
    private_nh_.param("resample_type", resample_type, std::string("lowvar"));
    if (resample_type == "lowvar") {
      resampler_ = new LowVarianceSampler<StationaryParticle>(neff);
      pf_ = new StandardPF<StationaryParticle>(mot_mod_, meas_mod_,
                                               resampler_, builder_);
    } else {
      ROS_ERROR_STREAM("Unknown resample type " << resample_type);
      ROS_BREAK();
    }
    
    pf_init_ = false;
  }
}
