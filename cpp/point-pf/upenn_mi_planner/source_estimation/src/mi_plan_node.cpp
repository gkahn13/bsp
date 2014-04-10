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

#include <algorithm>
#include <string>

#include <boost/thread/mutex.hpp>

#include "ros/ros.h"
#include "ros/console.h"

#include "tf/tf.h"

#include "geometry_msgs/PoseArray.h"
#include "geometry_msgs/Pose.h"

#include "nav_msgs/Odometry.h"

#include <actionlib/client/simple_action_client.h>

#include <player_map/rosmap.hpp>

#include <rvo_move/bot_client.hpp>

#include <rvo_move/MoveAction.h>

#include "particles.hpp"
#include "information.hpp"
#include "planner.hpp"
#include "util.hpp"

using namespace rf;

typedef actionlib::SimpleActionClient<rvo_move::MoveAction> ActionClient;

template <int D>
class MIPlanNode {
public:
  MIPlanNode(EntropyApproximator<D> *ea,
             LatticePlanner<RobotConfiguration> *planner, LOSChecker *checker) :
    private_nh_("~"), grid_(NULL), have_particles_(false), approx_(ea), planner_(planner),
    all_arrived_(true) {

    if (planner_ == NULL) {
      ROS_ERROR("Planner is NULL");
      ros::shutdown();
      return;
    }
    evaluator_ = MIEvaluator<D>::FromROS(private_nh_, approx_, checker);

    private_nh_.param("filter_frac", filter_frac_, 0.05);
    private_nh_.param("grid_res", grid_res_, 0.1);
    ROS_INFO("Particle recduction:");
    ROS_INFO("  Grid resolution: %f  CDF: %0.3f%%",
             grid_res_, filter_frac_ * 100.0);
    
    // Subscribe to topics and setup action servers
    pf_sub_ = nh_.subscribe("source_estimate", 5, &MIPlanNode::particlesCallback, this);
    bots_ = BotClient::MakeBots(private_nh_);
    ROS_INFO("Connecting to action servers...");
    for (size_t i = 0; i < bots_.size(); ++i) {
      std::string name = ros::names::append(ros::names::resolve(bots_[i]->getName()),
                                           std::string("move_server"));
      act_clients_.push_back(new ActionClient(name, true));
      ROS_INFO_STREAM("  Waiting for action server " << name);
      act_clients_[i]->waitForServer(ros::Duration(5.0));
      if (!act_clients_[i]->isServerConnected()) {
        ROS_ERROR("  Couldn't connect to server");
        ros::shutdown();
      } else {
        ROS_INFO("  Action server started");
      }
      arrived_.push_back(true);
      std::string p = ros::names::append(ros::names::resolve(bots_[i]->getName()), "points");
      
      plan_pubs_.push_back(nh_.advertise<geometry_msgs::PoseArray>(p, 5, true));
    }
  }

  ~MIPlanNode() {
    for (size_t i = 0; i < bots_.size(); ++i) {
      delete bots_[i];
    }
    delete evaluator_;
  }
  
  void particlesCallback(const geometry_msgs::PoseArray::ConstPtr &msg) {
    boost::mutex::scoped_lock lock(pf_mutex_);
    have_particles_ = true;
    particles_.update(msg);
    delete grid_;
    grid_ = particles_.FilterCDF(filter_frac_);
    ROS_INFO("Reduced particles from %i to %i", particles_.size, grid_->size);
  }

  void plan() {
    // Check if bots are ready & pipes are connected
    {
      boost::mutex::scoped_lock lock(pf_mutex_);
      bool ready = have_particles_;
      if (!ready) {
        ROS_WARN("Haven't received particles yet");
        return;
      }
    }
    
    bool all_ready = true;
    for (size_t i = 0; i < bots_.size(); ++i) {
      if (!bots_[i]->havePose()) {
        ROS_WARN("MIPlanNode::plan(): %s does not have pose", bots_[i]->getName().c_str());
        all_ready = false;
      }
    }
    if (!all_ready) {
      ROS_WARN("MIPlanNode::plan(): Not all bots are ready");
      return;
    }
    
    // // Check if all bots have arrived at desired position, if not, return
    // {
    //   boost::mutex::scoped_lock lock(movement_mutex_);
    //   if (all_arrived_) {
    //     // ROS_INFO("MIPlanNode::plan(): All bots have arrived at destination");
    //   } else {
    //     ROS_WARN("MIPlanNode::plan(): Waiting for all bots to arrive at commanded positions");
    //     return;
    //   }
    // }
    // Everything is ready, let's find a new place to go!
    execute();
  }

  void execute() {
    // find next best position and send bots there
    RobotConfiguration rc(bots_.size());
    for (size_t i = 0; i < bots_.size(); ++i) {
      geometry_msgs::Pose bot_pose = bots_[i]->getPose();
      rc.poses[i].setXYT(bot_pose.position.x, bot_pose.position.y, 0);
    }
    ROS_INFO("Generating plan");
    std::vector<RobotConfiguration *> possible_plans = planner_->generate(rc);
    if (possible_plans.size() == 0) {
        ROS_ERROR("No valid moves exist!");
        return;
    }

    // Publish waypoints for each robot
    geometry_msgs::PoseArray poses;
    poses.header.frame_id = "/map";
    poses.header.stamp = ros::Time::now();
    poses.poses.resize(possible_plans.size());  
    for (size_t bot = 0; bot < bots_.size(); ++bot) {
      for (size_t j = 0; j < possible_plans.size(); ++j) {
        double x = possible_plans[j]->poses[bot].x();
        double y = possible_plans[j]->poses[bot].y();
        // Different orientations so they all show up in RViz
        double t = static_cast<double>(bot) / bots_.size() * 2.0 * M_PI;
        tf::poseTFToMsg(tf::Pose(tf::createQuaternionFromYaw(t),
                                 tf::Vector3(x, y, 0)),
                        poses.poses[j]);
      }
      plan_pubs_[bot].publish(poses);
    }
    
    RobotConfiguration* best_plan;
    {
      boost::mutex::scoped_lock lock(pf_mutex_);
      ROS_INFO("Getting next best move (%zu options)", possible_plans.size());
      best_plan = evaluator_->findMax(possible_plans, *grid_);
    }

    if (best_plan == NULL) {
      ROS_ERROR("Internal error in MIEvaluator::findMax()");
    } else {
      // send plan to each of the robots
      boost::mutex::scoped_lock lock(movement_mutex_);
      for (size_t i = 0; i < bots_.size(); ++i) {
        arrived_[i] = false;
        sendGoal(i, best_plan->poses[i].x(), best_plan->poses[i].y());
      }
    }
    planner_->freeConfigurations(&possible_plans);
  }

  void sendGoal(size_t id, double x, double y) {
    // Assumes the caller has acquired movement_mutex_
    rvo_move::MoveGoal goal;
    goal.target_pose.pose.position.x = x;
    goal.target_pose.pose.position.y = y;
    act_clients_[id]->sendGoal(goal, boost::bind(&MIPlanNode::arriveCb, this, _1, _2, id));
  }

  void arriveCb(const actionlib::SimpleClientGoalState& state,
                const rvo_move::MoveResultConstPtr& result, size_t id) {
    boost::mutex::scoped_lock lock(movement_mutex_);
    arrived_[id] = true;
  }
  
  
private:
  ros::NodeHandle nh_, private_nh_;
  ros::Subscriber pf_sub_;
  boost::mutex pf_mutex_;

  ParticleArray<StationaryParticle> particles_;
  ParticleArray<StationaryParticle> *grid_;
  bool have_particles_;
    
  std::vector<BotClient *> bots_;

  EntropyApproximator<D> *approx_;
  MIEvaluator<D> *evaluator_;
  LatticePlanner<RobotConfiguration> *planner_;

  // movement_mutex_ guards the arrived_ vector and the variable all_arrived_
  boost::mutex movement_mutex_;
  std::vector<ActionClient*> act_clients_;
  std::vector<bool> arrived_;
  std::vector<ros::Publisher> plan_pubs_;
  
  bool all_arrived_;
  
  double short_dist_min_, short_dist_max_;
  double filter_frac_;
  double grid_res_;
};

LatticePlanner<RobotConfiguration> *getPlanner(map_t *map) {
  ros::NodeHandle nh("~");

  std::string type;
  nh.param("planner_type", type, std::string("kinematic"));
  
  if (type == "kinematic") {
    KinematicPlanner *kin = new KinematicPlanner(map);
    double step_size, min_sep, min_wall_dist;
    nh.param("step_size", step_size, 0.3);
    nh.param("min_sep", min_sep, 0.25);
    nh.param("min_wall_dist", min_wall_dist, 0.2);
    kin->setStepSize(step_size);
    kin->setMinSep(min_sep);
    kin->setMinWallDist(min_wall_dist, map->max_occ_dist);
    return kin;
  } else if (type == "skeleton") {
    std::string fname;
    double min_dist;
    double short_dist_min, short_dist_max;
    double merge_dist;
    nh.param("filename", fname, std::string("graph.txt"));
    nh.param("min_travel_dist", min_dist, 2.0);
    SkeletonPlanner *skel = SkeletonPlanner::FromCSV(fname.c_str(), map);
    skel->setMinDist(min_dist);
    nh.param("short_dist_min", short_dist_min, 0.2);
    nh.param("short_dist_max", short_dist_max, 4.0);

    nh.param("merge_dist", merge_dist, 0.0);
    ROS_INFO("Short dist min: %0.2f Short dist max: %0.2f", short_dist_min, short_dist_max);
    skel->setShortDists(short_dist_min, short_dist_max);
    skel->setMergeDist(merge_dist);
    return skel;
  } else {
    ROS_ERROR("Unrecognized planner type: %s", type.c_str());
    return NULL;
  }
}

#define NUM_BOTS 2

int main(int argc, char **argv) {
  ros::init(argc, argv, "mi_plan");
  ros::NodeHandle n;
  ros::NodeHandle pnh("~");

  int nbots;
  pnh.param("nbots", nbots, 0);
  if (nbots != NUM_BOTS) {
    ROS_ERROR("mi_plan_node compiled for %i bots, but ~nbots=%i",
              NUM_BOTS, nbots);
    ROS_BREAK();
  }
  // Get the map for planning
  map_t *map = requestCSpaceMap("/static_map");  
  LatticePlanner<RobotConfiguration> *p = getPlanner(map);

  LOSChecker checker(map);

  TaylorApproximator<NUM_BOTS> ta(0, 0);
  MIPlanNode<NUM_BOTS> planner(&ta, p, &checker);


  double period;
  pnh.param("period", period, 4.0);
  for (ros::Duration d(period); ros::ok(); d.sleep()) {
    ros::spinOnce();
    planner.plan();
    ros::spinOnce();
  }
  map_free(map);
  delete p;
  
  return 0;
}
