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

#include <cstring>
#include <fstream>
#include <iostream>
#include <vector>
#include <cassert>
#include <cmath>
#include <cstdlib>

#include "planner.hpp"

using namespace std;
using namespace rf;

RobotConfiguration::RobotConfiguration(int nb) : nbots(nb) {
  poses = new Pose2D[nbots];    
}

RobotConfiguration::RobotConfiguration(const RobotConfiguration &rc) :
  nbots(rc.nbots) {
  poses = new Pose2D[nbots];
  for (int i = 0; i < nbots; ++i) {
    poses[i].setXYT(rc.poses[i]);
  }
}

bool RobotConfiguration::samePosition(const RobotConfiguration &rc, double dist /* = 0.1 */) {
  bool eq = rc.nbots == nbots;
  double sq_dist = dist * dist;
  for (int i = 0; i < nbots && eq; ++i) {
    eq &= rc.poses[i].distance_sq(poses[i]) < sq_dist;
  }
  return eq;
}

RobotConfiguration::~RobotConfiguration() {
  delete[] poses;
}
  
ostream &operator<<(ostream &os, const RobotConfiguration &r) {
  for (int i = 0; i < r.nbots; ++i) {
    os << "Bot " << i << ": " << r.poses[i] << endl;
  }
  return os;
}

DynamicConfiguration::DynamicConfiguration(int nb) : RobotConfiguration(nb) {
  velocities = new Velocity2D[nbots];
  for (int i = 0; i < nb; ++i) {
    velocities[i].v = 0;
    velocities[i].w = 0;
  }
}


DynamicConfiguration::DynamicConfiguration(const DynamicConfiguration &dc) :
  RobotConfiguration(dc) {
  velocities = new Velocity2D[nbots];
  memcpy(velocities, dc.velocities, dc.nbots * sizeof(velocities[0]));
}


DynamicConfiguration::~DynamicConfiguration() {
  delete[] velocities;
}

ostream& operator<<(ostream& os, const DynamicConfiguration &d) {
  for (int i = 0; i < d.nbots; ++i) {
    os << "Bot " << i;
    os << "  Pose:     " << d.poses[i];
    os << "  Velocity: " << "(" << d.velocities[i].v << ", " << d.velocities[i].w << ")";
    os << endl;
  }
  return os;
}

//============================ KinematicPlanner =============================//


KinematicPlanner::KinematicPlanner() :
  LatticePlanner<RobotConfiguration>(), step_(0.25)  {
    
}

KinematicPlanner::KinematicPlanner(map_t *c) :
  LatticePlanner<RobotConfiguration>(c), step_(0.25) {

}

vector<RobotConfiguration *>
KinematicPlanner::generate(const RobotConfiguration &rc) const {
  vector<RobotConfiguration *> configurations;
  RobotConfiguration start_config(rc);
  generate(rc, 0, &start_config, &configurations);
  return configurations;
}


void KinematicPlanner::generate(const RobotConfiguration &rc,
                                int bot_ind,
                                RobotConfiguration *current_config,
                                vector<RobotConfiguration *> *configs) const {
  if (bot_ind >= rc.nbots) {
    configs->push_back(new RobotConfiguration(*current_config));
    return;
  }
  for (int choice = 0; choice < 4; ++choice) {
    int sign;
    double *element;
    Pose2D *pose = &current_config->poses[bot_ind];
    double orig_x = pose->x(), orig_y = pose->y(), orig_t = pose->theta();
    double new_x = orig_x, new_y = orig_y, new_t = orig_t;
    switch(choice) {
    case UP:
      element = &new_y;
      sign = 1;
      break;
    case DOWN:
      element = &new_y;
      sign = -1;
      break;
    case LEFT:
      element = &new_x;
      sign = -1;
      break;
    case RIGHT:
      element = &new_x;
      sign = 1;
      break;
    default:
      abort();
    }
    (*element) += step_ * sign;
    pose->setXYT(new_x, new_y, new_t);
    bool map_free = true;
    if (map_ != NULL) {
      map_cell_t *cell = map_get_cell(map_, new_x, new_y, 0);
      assert(cell != NULL);
      // ROS_ERROR("occ_state = %i occ_dist: %f", cell->occ_state, cell->occ_dist);
      map_free = cell->occ_state == -1 && (cell->occ_dist > occ_dist_thresh_ || abs(cell->occ_dist - max_occ_dist_) < 1e-10);
    }
    if (map_free &&
        !collision(rc, *current_config, bot_ind) &&
        !overlap(*current_config, bot_ind)) {
      generate(rc, bot_ind + 1, current_config, configs);
    }
    pose->setXYT(orig_x, orig_y, orig_t);
  }
}

//============================= SkeletonPlanner =============================//

SkeletonPlanner::SkeletonPlanner() :
  LatticePlanner<RobotConfiguration>(), min_traversal_dist_(0.0),
  shortdist_min_(0.0), shortdist_max_(0.0),
  merge_dist_(0) {
  
}

SkeletonPlanner::SkeletonPlanner(map_t *map) :
  LatticePlanner<RobotConfiguration>(map) {
  SkeletonPlanner();
}

SkeletonPlanner *SkeletonPlanner::FromCSV(const char *fname, map_t *map) {
  SkeletonPlanner *sp = new SkeletonPlanner(map);
  std::string line;
  std::ifstream is(fname);
  if (!is) {
    ROS_ERROR("SkeletonPlanner::FromCSV()  Problem opening '%s'", fname);
    ROS_BREAK();
  }
  while (!is.eof() && getline(is, line)) {
    vector<string> tokens;
    stringstream ss;
    string t;

    ss.str(line);
    while (getline(ss, t, ' ')) {
      tokens.push_back(t);
    }
    double x, y;
    std::vector<int> inds;
    if (sscanf(tokens[0].c_str(), "%lf", &x) == 0) {
      ROS_ERROR("SkeletonPlanner::FromCSV() Problem parsing");
    }
    if (sscanf(tokens[1].c_str(), "%lf", &y) == 0) {
      ROS_ERROR("SkeletonPlanner::FromCSV() Problem parsing");
    }
    sp->vertices_.push_back(Point2D(x, y));
    if (map != NULL) {
      map_cell_t *cell = map_get_cell(map, x, y, 0);
      if ((cell == NULL) || cell->occ_state != -1) {
        ROS_ERROR("SkeletonPlanner::FromCSV() Bad coordinate: %f %f", x, y);
      }
    }
    stringstream debug;
    debug.str();
    for (size_t i = 2; i < tokens.size(); ++i) {
      int j;
      if (sscanf(tokens[i].c_str(), "%i ", &j) == 0) {
        ROS_ERROR("SkeletonPlanner::FromCSV() Problem parsing");
      }
      inds.push_back(j - 1);
      debug << j << " ";
    }

    ROS_DEBUG("Adding %zu with neighbors %s", 1 + sp->neighbors_.size(), debug.str().c_str());
    
    sp->neighbors_.push_back(inds);
  }

  // Consistency check
  bool inconsistent = false;
  for (size_t ver = 0; ver < sp->neighbors_.size(); ++ver) {
    vector<int> ver_ns = sp->neighbors_[ver];
    
    for (size_t j = 0; j < ver_ns.size(); ++j) {
      int neigh = ver_ns[j];
      if (neigh < 0 || neigh > (int)sp->neighbors_.size()) {
        ROS_ERROR("Invalid neighbor index %i for %zu", neigh + 1, ver + 1);
        continue;
      }
      vector<int> neigh_ns = sp->neighbors_[neigh];
      int num = count(neigh_ns.begin(), neigh_ns.end(), ver);
      if (num == 0) {
        ROS_ERROR("Edges are not symmetric %lu <-> %i", ver + 1, neigh + 1);
        inconsistent = true;
      } else if (num > 1) {
        ROS_ERROR("Same number appears multiple times in %i", neigh + 1);
        inconsistent = true;
      }
    }
  }
  if (inconsistent) {
    ROS_ERROR("Graph is inconsistent");
    ROS_BREAK();
  }
  return sp;
}

std::vector<RobotConfiguration *> SkeletonPlanner::generate(const RobotConfiguration &rc) const {
  vector<RobotConfiguration *> configs;

  vector<int> config(rc.nbots, 0), current(0, 0);
  // Get index for each position
  for (int bot = 0; bot < rc.nbots; ++bot) {
    int ind;
    ind = -1;
    double closest_dist = numeric_limits<double>::infinity();

    for (size_t vert = 0; vert < vertices_.size(); ++vert) {
      double dist = l2_squared(vertices_[vert], rc.poses[bot]);
      if (dist < closest_dist) {
        ind = vert;
        closest_dist = dist;
      }
    }
    if (ind == -1) {
      ROS_ERROR("SkeletonPlanner() Failed to find closest vertex");
    } else {
      config[bot] = ind;
    }
  }
  stringstream ss;
  for (size_t i = 0; i < config.size(); ++i) {
    ss << config[i] << " ";
  }

  vector<vector<int> > to_visit;
  for (size_t i = 0; i < config.size(); ++i) {
    // Get neighbors within distance threshold and merge them;
    vector<int> neighbors = getNeighbors(config[i]);
    vector<int> merged = mergeNeighbors(neighbors);
    
    to_visit.push_back(merged);

    std::stringstream ss;
    ss.str();
    for (size_t j = 0; j < merged.size(); ++j) {
      ss << merged[j] + 1 << " ";
    }
    ROS_DEBUG("%zu neighbors: %s", merged.size(), ss.str().c_str());
  }
  
  generate(config, to_visit, &current, &configs);

  return configs;
}

vector<int> SkeletonPlanner::mergeNeighbors(const vector<int> &ns) const {
  vector<bool> keep;
  keep.resize(ns.size());
  for (size_t i = 0; i < ns.size(); ++i) {
    keep[i] = true;
  }

  vector<int> compact;
  for (size_t i = 0; i < ns.size(); ++i) {
    if (keep[i]) {
      compact.push_back(ns[i]);
    }
    for (size_t j = i+1; j < ns.size(); ++j) {
      if (sqrt(l2_squared(vertices_[ns[i]], vertices_[ns[j]])) < merge_dist_) {
        keep[j] = false;
      }
    }
  }
  return compact;
}

void SkeletonPlanner::generate(const vector<int> &bots, const vector<vector<int> > &neighbors,
                               vector<int> *current,
                               vector<RobotConfiguration *> *configs) const {
  size_t curr_ind = current->size();
  // If at the end, build XY coordinates of configuration
  if (curr_ind == bots.size()) {
    RobotConfiguration *rc = new RobotConfiguration(bots.size());
    stringstream ss;
    for (size_t i = 0; i < curr_ind; ++i) {
      int vert_ind = (*current)[i];
      rc->poses[i].setXYT(vertices_[vert_ind].x, vertices_[vert_ind].y, 0);
      ss << vert_ind << " ";
    }
    configs->push_back(rc);
    return;
  }

  // Go over all edges for this bot, recursing to next bot for each one
  std::vector<int> edges = neighbors[curr_ind];
  
  // Go over all edges for this bot, recursing to next bot for each one
  for (size_t e_ind = 0; e_ind < edges.size(); e_ind++) {
    current->push_back(edges[e_ind]);

    generate(bots, neighbors, current, configs);
    assert(current->size() - 1 == curr_ind);
    
    current->pop_back();
  }
}


std::set<int> SkeletonPlanner::getNeighbors(int start_ind, double min_dist, double max_dist) const {
  // Generate neighbors.  Do BFS.  Take nodes that are leaves in the traversal,
  // or take a node if it is above a given distance.
  
  // Initialize queue
  list<int> queue;
  queue.push_back(start_ind);
  // Initialize visited; which nodes have been visited
  vector<bool> visited;
  visited.resize(neighbors_.size());  
  for (size_t i = 0; i < visited.size(); ++i) {
    visited[i] = false;
  }
  // in_queue tracks which elements of graph are in the queue
  vector<bool> in_queue;
  in_queue.resize(neighbors_.size());
  for (size_t i = 0; i < visited.size(); ++i) {
    in_queue[i] = false;
  }
  in_queue[start_ind] = true;

  Point2D base_xy = vertices_[start_ind];

  ROS_DEBUG("Starting at %i", start_ind + 1);
  set<int> edges;
  double node_dist;
  while (!queue.empty()) {
    // Get first element
    bool is_leaf;
    int node;
    vector<int> node_neighs;
    node = queue.front();
    queue.pop_front();
    ROS_DEBUG("  Visiting %2i", node + 1);
    if (visited[node] == true || in_queue[node] == false) {
      ROS_ERROR("Bad values for the queue!");
      ROS_BREAK();
    } else {
      visited[node] = true;
      in_queue[node] = false;
    }

    // Distance from node to bot
    node_dist = sqrt(l2_squared(base_xy, vertices_[node]));
    
    // Neighbors of node we're currently visitng
    node_neighs = neighbors_.at(node);
    std::stringstream ss;
    ss.str();
    for (size_t i = 0; i < node_neighs.size(); ++i) {
      ss << node_neighs[i] + 1 << " ";
    }
    ROS_DEBUG("  Neighbors of %i: %s", node + 1, ss.str().c_str());

    // Is node a leaf in the BFS traversal?
    is_leaf = true;
    for (size_t i = 0; i < node_neighs.size(); ++i) {
      is_leaf &= (visited[node_neighs[i]] || in_queue[node_neighs[i]]);
    }

    if (is_leaf || (min_dist < node_dist && node_dist < max_dist)) {
      edges.insert(node);
      ROS_DEBUG("    Adding %2i", node + 1);
      continue;
    } else if (max_dist < node_dist) {
      continue;
    } else {
      int neigh;
      for (size_t i = 0; i < node_neighs.size(); ++i) {
        neigh = node_neighs[i];
        if (!visited[neigh] && !in_queue[neigh]) {
          in_queue[neigh] = true;
          queue.push_back(neigh);
          ROS_DEBUG("    Inserting %2i into queue", neigh + 1);
        }
      }
    }
  }
  return edges;
}


std::vector<int> SkeletonPlanner::getNeighbors(int start_ind) const {
  // Generate neighbors.  Do BFS.  Take nodes that are leaves in the traversal,
  // or take a node if it is above a given distance.

  set<int> short_edges, far_edges;
  set<int> combined;
  vector<int> edges;
  far_edges = getNeighbors(start_ind, min_traversal_dist_, 50);
  short_edges = getNeighbors(start_ind, shortdist_min_, shortdist_max_);
  set<int>::iterator it;
  for (it = short_edges.begin(); it != short_edges.end(); ++it) {
    combined.insert(*it);
    // ROS_WARN("Short Edge %i", *it + 1);
  }
  for (it = far_edges.begin(); it != far_edges.end(); ++it) {
    combined.insert(*it);
    // ROS_WARN("Long Edge %i", *it + 1);
  }
  std::stringstream ss;
  for (it = combined.begin(); it != combined.end(); ++it) {
    ss << (*it) + 1 << " ";
    edges.push_back(*it);
  }

  return edges;
}
                                        
//============================= DynamicPlanner ==============================//

DynamicPlanner::DynamicPlanner() :
  LatticePlanner<DynamicConfiguration>(), dt_(0.5)  {
  setVel(-0.2, 0.2, 0.2, 0.4);
  setOmega(-0.1, 0.1, 0.1, 0.2);
}

DynamicPlanner::DynamicPlanner(map_t *c) :
  LatticePlanner<DynamicConfiguration>(c), dt_(0.5) {
  setVel(-0.2, 0.2, 0.2, 0.4);
  setOmega(-0.1, 0.1, 0.1, 0.2);
}

vector<DynamicConfiguration *> DynamicPlanner::generate(const DynamicConfiguration &rc) const {
  vector<DynamicConfiguration *> configurations;
  DynamicConfiguration start_config(rc);
  generate(rc, 0, &start_config, &configurations);
  return configurations;
}

void DynamicPlanner::generate(const DynamicConfiguration &rc,
                              int bot_ind,
                              DynamicConfiguration *current_config,
                              vector<DynamicConfiguration *> *configs) const {
  if (bot_ind >= rc.nbots) {
    configs->push_back(new DynamicConfiguration(*current_config));
    return;
  }
  for (double dv = vel_.step_min; dv <= vel_.step_max; dv += vel_.step_size) {
    for (double dw = omega_.step_min; dw <= omega_.step_max; dw += omega_.step_size) {
      // Calculate new position based on robot's previous velocities, time
      // step, and new command
      Pose2D *pose = &current_config->poses[bot_ind];
      Velocity2D *curr_vel = &current_config->velocities[bot_ind];
      double orig_v = curr_vel->v, orig_w = curr_vel->w;
      double orig_x = pose->x(), orig_y = pose->y(), orig_t = pose->theta();
      // Insure that robot does not drive backwards
      if (curr_vel->v + dv < 0) {
        continue;
      }
      
      double new_v = curr_vel->v + dv;
      double new_w = curr_vel->w + dw;
      curr_vel->v = new_v;
      curr_vel->w = new_w;
      ROS_ERROR("new_v: %6.2f new_w: %6.2f", curr_vel->v, curr_vel->w);
      // Calculate how state would change if we change velocity
      double dx, dy, dth;      
      dx = new_v * (dt_ - pow(new_w, 2)) * (pow(dt_, 3)) / 6.0;
      dy = new_v * ((new_w * pow(dt_, 2) / 2.0) -
                    (pow(new_w, 3) * pow(dt_, 4) / 24.0));
      dth = new_w * dt_;
      ROS_ERROR("dx: %6.2f dy: %6.2f dt: %6.2f", dx, dy, dth);
      // Calculate new state based on change
      double new_x, new_y, new_t;      
      new_x = orig_x + dx*cos(orig_t) - dy*sin(orig_t);
      new_y = orig_y + dx*sin(orig_t) - dy*cos(orig_t);
      new_t = orig_t + dth;
      ROS_ERROR("new_x: %6.2f new_y: %6.2f new_t: %6.2f", new_x, new_y, new_t);
      pose->setXYT(new_x, new_y, new_t);
      // Verify that new state wouldn't collide with other bots or walls.
      bool map_free = true;
      if (map_ != NULL) {
        map_cell_t *cell = map_get_cell(map_, new_x, new_y, 0);
        assert(cell != NULL);
        // ROS_ERROR("occ_state = %i occ_dist: %f", cell->occ_state, cell->occ_dist);
        map_free = (cell->occ_state == -1 &&
                    (cell->occ_dist > occ_dist_thresh_ || abs(cell->occ_dist - max_occ_dist_) < 1e-10));
      }
      if (map_free &&
          !collision(rc, *current_config, bot_ind) &&
          !overlap(*current_config, bot_ind)) {
        generate(rc, bot_ind + 1, current_config, configs);
      }
      // Reset state to original
      pose->setXYT(orig_x, orig_y, orig_t);
      curr_vel->v = orig_v;
      curr_vel->w = orig_w;
    }
  }
}
