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

#include <rosbag/bag.h>
#include <rosbag/view.h>

#include <boost/foreach.hpp>

#include <std_msgs/String.h>
#include <nanotron/RangeWithPose.h>

#include <boost/foreach.hpp>

#include "estimate.hpp"

using namespace std;

int main(int argc, char **argv) {
  ros::init(argc, argv, "static_source");
  ros::NodeHandle n("~");
  if (!n.hasParam("target_id")) {
    ROS_ERROR("Must specify target_id");
    return -1;
  }
  if (argc != 2) {
    ROS_WARN("usage: static_source_bag bagfile");
    return -1;
  }
  int target_id;
  n.param("target_id", target_id, 0);
  std::stringstream ss;
  ss << target_id << "_";
  std::string target_str(ss.str());
  
  std::time_t t = std::time(NULL);
  char s[100];
  std::strftime(s, 80, "%Y-%m-%d--%H-%M-%S", std::localtime(&t));
  string bag_path = string(argv[1]);
  string folder = bag_path.substr(0, bag_path.rfind("/"));
  if (folder.size() != 0) {
    folder += "/";
  }
  string logfile = folder + target_str + string(s) + ".log";
  string bag_out_name = folder + target_str + string(s) + ".bag";
  
  rosbag::Bag bag_in, bag_out;
  bag_in.open(bag_path, rosbag::bagmode::Read);
  bag_out.open(bag_out_name.c_str(), rosbag::bagmode::Write);
  
  n.setParam("logfile", logfile);

  rf::StaticSourceNode node;

  vector<string> topics;
  string msmt_topic = "/measurements";
  topics.push_back(msmt_topic);
  for (int i = 20; i < 35; ++i) {
    stringstream ss;
    ss.str("");
    ss << "/scarab" << i << "/amcl_pose";
    string topic = ss.str();
    topics.push_back(topic);
  }
  rosbag::View view(bag_in, rosbag::TopicQuery(topics));

  bool added_metadata = false;
  BOOST_FOREACH(rosbag::MessageInstance const m, view) {
    if (!added_metadata) {
      std_msgs::String s;
      s.data = "Copied from " + string(argv[1]);
      bag_out.write("/metadata", m.getTime(), s);
      added_metadata = true;
    }
    
    if (m.getTopic() == msmt_topic) {
      // Apply message to filter and then write new estimate + measurement to bag
      nanotron::RangeWithPose::ConstPtr msg = m.instantiate<nanotron::RangeWithPose>();
      if (msg == NULL){
        ROS_ERROR("Found a non range with pose message in bag on /measurements");
        continue;
      }

      node.measurement_callback(msg);

      bag_out.write(m.getTopic(), msg->range.header.stamp, msg);
      // if (node.updated()) {
      //   bag_out.write("/laptop/source_estimate", msg->range.header.stamp, node.getParticles());
      // }
    } else {
      // Copy pose information
      geometry_msgs::PoseWithCovarianceStamped::ConstPtr pose =
        m.instantiate<geometry_msgs::PoseWithCovarianceStamped>();
      bag_out.write(m.getTopic(), pose->header.stamp, pose);
    }
  }
  return 0;
}
