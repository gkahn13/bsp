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

#include <player_map/rosmap.hpp>
#include <nanotron/RangeWithPose.h>
#include <boost/foreach.hpp>
#include <string>

#include "types.hpp"

using namespace std;
using namespace rf;

// list of (time, position)
typedef list<pair<double, Point2D> > path_t;

Point2D nearestPoint(std::map<int, path_t*> bot_paths, int id, double ts, int *status) {
  path_t* path = bot_paths[id];
  if (path == NULL) {
    *status = -1;
    return Point2D();
  }
  path_t::iterator prev = path->begin();
  for (path_t::iterator it = path->begin(); it != path->end(); ++it) {
    if (ts <= (*it).first) {
      if (fabs(((*it).first) - ts) > 5.0) {
        *status = -1;
      } else {
        *status = 0;
      }
      return Point2D(((*it).second.x), ((*it).second.y));      
    }
  }
  *status = -1;
  return Point2D();
}

Point2D targetLocation(int target_id, const rosbag::Bag &bag) {
  stringstream ss;
  ss << "/scarab" << target_id << "/amcl_pose";
  string target_topic(ss.str());
  vector<string> topics;
  topics.push_back(target_topic);
  rosbag::View view(bag, rosbag::TopicQuery(topics));

  vector<Point2D> points;
  BOOST_FOREACH(rosbag::MessageInstance const m, view) {
    geometry_msgs::PoseWithCovarianceStamped::ConstPtr s =
      m.instantiate<geometry_msgs::PoseWithCovarianceStamped>();
    if (s == NULL){
      ROS_ERROR_STREAM("Found a non pose message on " << target_topic);
      continue;
    }
    points.push_back(Point2D(s->pose.pose.position.x, s->pose.pose.position.y));    
  }

  if (points.size() == 0) {
    ROS_ERROR_STREAM("No pose information on " << target_topic);
	points.push_back(Point2D(53.5, 51.0));
    // ROS_BREAK();
  }

  Point2D base = points[0];
  for (size_t i = 0; i < points.size(); ++i) {
    double dist = sqrt(l2_squared(base, points[i]));
    if (dist > 0.1) {
      ROS_WARN("Target moved %0.2f from base location", dist);
    }
  }
  
  return points[0];
}


int main(int argc, char **argv) {
  ros::init(argc, argv, "bagread");
  ros::Time::init();
  rosbag::Bag bag;
  if (argc != 3) {
    ROS_WARN("usage: bagread bagfile target_id");
    return -1;
  }
  bag.open(argv[1], rosbag::bagmode::Read);
  double target_id = atoi(argv[2]);

  Point2D target_xy = targetLocation(target_id, bag);
  
  std::string fname = std::string(argv[1]) + ".txt";
  
  std::vector<std::string> topics;
  topics.push_back(std::string("/measurements"));

  rosbag::View view(bag, rosbag::TopicQuery(topics));

  map_t *map = requestCSpaceMap("/range_static_map");  
  LOSChecker los(map);
  // map from robot id -> list of (time, position)
  std::map<int, path_t*> bot_paths;

  // Gather all robot locations at each time
  BOOST_FOREACH(rosbag::MessageInstance const m, view) {
    nanotron::RangeWithPose::ConstPtr s = m.instantiate<nanotron::RangeWithPose>();
    if (s == NULL){
      ROS_ERROR("Found a non range with pose message in bag on /measurements");
      continue;
    }

    if (bot_paths.find(s->range.src_id) == bot_paths.end()) {
      bot_paths[s->range.src_id] = new path_t();
    }
    Point2D xy(s->pose.pose.pose.position.x,
               s->pose.pose.pose.position.y);
    bot_paths[s->range.src_id]->push_back(make_pair(s->pose.header.stamp.toSec(),
                                                    xy));
  }


  FILE *f = fopen(fname.c_str(), "w");
  
  BOOST_FOREACH(rosbag::MessageInstance const m, view) {
    nanotron::RangeWithPose::ConstPtr s = m.instantiate<nanotron::RangeWithPose>();
    if (s == NULL){
      ROS_ERROR("Found a non range with pose message in bag on /measurements");
      continue;
    }

    Point2D dest_xy;
    if (s->range.dest_id != target_id) {
      int status;
      Point2D temp = nearestPoint(bot_paths, s->range.dest_id, s->pose.header.stamp.toSec(),
                                     &status);
      dest_xy.x = temp.x;
      dest_xy.y = temp.y;
      if (status != 0) {
        std::cerr << "Nearest pose is too far away" << std::endl;
        continue;
      }
    } else {
      dest_xy.x = target_xy.x;
      dest_xy.y = target_xy.y;
    }
    Point2D source_xy(s->pose.pose.pose.position.x, s->pose.pose.pose.position.y);
    
    bool has_los = los.LineOfSight(dest_xy.x, dest_xy.y, source_xy.x, source_xy.y);

    fprintf(f, "%0.2f %i % 8.3f % 8.3f %i % 8.3f % 8.3f % 8.4f %i\n",
            s->pose.header.stamp.toSec(),
            s->range.src_id, source_xy.x, source_xy.y,
            s->range.dest_id, dest_xy.x, dest_xy.y,
            s->range.distance, has_los);
  }
  fclose(f);
  
  bag.close();
}
