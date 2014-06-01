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

#include <ros/ros.h>

#include <boost/thread/mutex.hpp>
#include <geometry_msgs/PoseWithCovarianceStamped.h>
#include <geometry_msgs/Pose.h>

#include <player_map/rosmap.hpp>

#include <rvo_move/bot_client.hpp>

#include "nanotron/Range.h"

#include "util.hpp"

using namespace geometry_msgs;
using namespace rf;

class BotSim : public BotClient {
public:
  BotSim(const ros::NodeHandle &parent, std::string prefix, std::string id) :
    BotClient(parent, prefix) {
    range_pub_ = nh_->advertise<nanotron::Range>("range", 2);
    if (!parent.hasParam(id)) {
      ROS_ERROR("No id=%s", id.c_str());
      ROS_BREAK();
    }
    parent.getParam(id, src_id_);
  }

  void sendRange(int dest_id, double d) {
    nanotron::Range range;
    range.src_id = src_id_;
    range.dest_id = dest_id;
    range.distance = d;
    range.header.stamp = ros::Time::now();
    
    range_pub_.publish(range);
  }
  
  int id() const { return src_id_; }

private:
  ros::Publisher range_pub_;

  int src_id_;
};

int main(int argc, char **argv) {
  ros::init(argc, argv, "range_sim");
  ros::NodeHandle nh("~");
  
  std::vector<BotSim *> bots;
  LOSChecker *los_checker;
  int nbots;
  map_t *map;  
  std::string map_name = "map_valid";
  double transmit_period;

  nh.param("period", transmit_period, 0.5);
  nh.getParam("nbots", nbots);

  map = requestCSpaceMap(map_name.c_str());
  los_checker = new LOSChecker(map);
  LosNlosFading<StationaryParticle> *model =
    LosNlosFading<StationaryParticle>::FromROS(nh, los_checker, NULL);
  const LosNlosFading<StationaryParticle>::Params &params = model->getParams();
  
  for (int i = 0; i < nbots; ++i) {
    std::stringstream ss_name, ss_id;
    ss_name.str();
    ss_name << "robot" << i;
    ss_id.str();
    ss_id << "id" << i;
    bots.push_back(new BotSim(nh, ss_name.str(), ss_id.str()));
  }
  
  ros::Rate r(1.0 / transmit_period);
  for (; ros::ok(); ros::spinOnce()) {
    // Send measurements
    Pose dest_pose, source_pose;    
    for (int dest_ind = 0; dest_ind < nbots; ++dest_ind) {
      if (!bots[dest_ind]->havePose()) {
        ROS_WARN("sim_range: Skipping measurements involving Bot %i (name=%s); no pose",
                 dest_ind, bots[dest_ind]->getName().c_str());
        continue;
      }
      int dest_id = bots[dest_ind]->id();
      
      dest_pose = bots[dest_ind]->getPose();      
      for (int source_ind = 0; source_ind < nbots; ++source_ind) {
        if (source_ind == dest_ind) {
          continue;
        }
        BotSim *bot = bots[source_ind];
        if (!bot->havePose()) {
          ROS_WARN("sim_range: Skipping measurements involving Bot %i (name=%s); no pose",
                   source_ind, bot->getName().c_str());
          continue;
        }
        source_pose = bot->getPose();
        double dist = hypot(dest_pose.position.x - source_pose.position.x,
                            dest_pose.position.y - source_pose.position.y);
        bool los = los_checker->LineOfSight(dest_pose.position.x, dest_pose.position.y,
                                            source_pose.position.x, source_pose.position.y);
        double samp;
        if (los) {
          samp = params.los.mu0 + params.los.mu * dist +
            gsl::normal(sqrt(dist * params.los.sigma2));
        } else {
          samp = params.nlos.mu0 + params.nlos.mu * dist +
            gsl::normal(sqrt(dist * params.nlos.sigma2));
        }
        double result = dist - samp;
        bot->sendRange(dest_id, result);
      }
      r.sleep();
    }
  }
}
