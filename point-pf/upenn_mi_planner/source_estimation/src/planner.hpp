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

#ifndef PLANNER_HPP
#define PLANNER_HPP

#include <iostream>
#include <vector>
#include <map>

#include "ros/ros.h"
#include "ros/console.h"

#include <player_map/rosmap.hpp>

#include "types.hpp"

// Configuration of robots in world coordinates.

class RobotConfiguration {
public:
    RobotConfiguration(int nb);
    RobotConfiguration(const RobotConfiguration &rc);

    virtual ~RobotConfiguration();

    bool samePosition(const RobotConfiguration &rc, double dist = 0.1);
    const int nbots;
    rf::Pose2D *poses;
};

std::ostream &operator<<(std::ostream &os, const RobotConfiguration &r);

class DynamicConfiguration : public RobotConfiguration {
public:
    DynamicConfiguration(int nb);
    DynamicConfiguration(const DynamicConfiguration &dc);

    ~DynamicConfiguration();

    rf::Velocity2D *velocities;
};

std::ostream &operator<<(std::ostream &os, const DynamicConfiguration &dc);

//============================= LatticePlanner ==============================//

template <typename T>
class LatticePlanner {
public:
    LatticePlanner();
    LatticePlanner(map_t *map_);

    virtual ~LatticePlanner() {};

    // Get all valid configurations of robots that are one step away from the
    // current configuration.  If collision map is populated, all returned
    // configurations will not have any collisions with it.
    // 
    // Caller is responsible for de-allocating configurations, which you can do
    // via freeConfigurations().
    virtual std::vector<T *> generate(const T &rc) const = 0;
    // Free all configurations contained in the passed in vector.
    void freeConfigurations(std::vector<T *> *configs) const;

    void setMinSep(double sep) {
        assert(sep > 0);
        min_sep_seq_ = sep * sep;
    }

  // Keep at least this distance away from walls
  void setMinWallDist(double sep, double max_occ_dist) {
    assert(sep > 0);
    occ_dist_thresh_ = sep;
    max_occ_dist_ = max_occ_dist;
  }
    
protected:
    bool collision(const RobotConfiguration &prev,
                   const RobotConfiguration &curr,
                   int max_ind) const;
    bool overlap(const RobotConfiguration &current, int max_ind) const;
    map_t *map_;
    double min_sep_seq_;
    double occ_dist_thresh_;
    double max_occ_dist_;
};

template <typename T>
LatticePlanner<T>::LatticePlanner() :
  map_(NULL) , min_sep_seq_(0.5), occ_dist_thresh_(0.0), max_occ_dist_(0.0) {

}

template <typename T>
LatticePlanner<T>::LatticePlanner(map_t *c) :
  map_(c), min_sep_seq_(0.5), occ_dist_thresh_(0.1), max_occ_dist_(0.0) {

}


template <typename T>
void LatticePlanner<T>::freeConfigurations(std::vector<T *> *cs) const {
    for (size_t i = 0; i < cs->size(); ++i) {
        delete(cs->at(i));
        cs->at(i) = NULL;
    }
}

template <typename T>
bool LatticePlanner<T>::collision(const RobotConfiguration &prev,
                                  const RobotConfiguration &curr,
                                  int max_ind) const {
    for (int i = 0; i <= max_ind; ++i) {
        for (int j = i + 1; j <= max_ind; ++j) {
            bool i_to_j = prev.poses[i].distance_sq(curr.poses[j]) < min_sep_seq_;
            bool j_to_i = prev.poses[j].distance_sq(curr.poses[i]) < min_sep_seq_;
            if (i_to_j && j_to_i) {
                return true;
            }
        }
    }
    return false;
}

template <typename T>
bool LatticePlanner<T>::overlap(const RobotConfiguration &curr, int max_ind) const {
    for (int i = 0; i <= max_ind; ++i) {
        for (int j = i + 1; j <= max_ind; ++j) {
            if (curr.poses[i].distance_sq(curr.poses[j]) < min_sep_seq_) {
                return true;
            }
        }
    }
    return false;
}

//============================ KinematicPlanner =============================//

class KinematicPlanner : public LatticePlanner<RobotConfiguration> {
public:
    KinematicPlanner();
    KinematicPlanner(map_t *map_);

    std::vector<RobotConfiguration *> generate(const RobotConfiguration &rc) const;

    void setStepSize(double step) {
        assert(step > 0);
        step_ = step;
    }
    
private:
    enum NEIGHBOR_4 {UP, DOWN, LEFT, RIGHT};
    void generate(const RobotConfiguration &rc,
                  int bot_ind,
                  RobotConfiguration *current_config,
                  std::vector<RobotConfiguration *> *configs) const;
    double step_;
};

// Plan using a pre-defined set of points (i.e. a skeleton).  The points are in
// 2D, and have a user-specified connectivity structure.
class SkeletonPlanner : public LatticePlanner<RobotConfiguration> {
public:
  // Factory method
  static SkeletonPlanner* FromCSV(const char *s, map_t *map);

  void setMergeDist(double dist) { merge_dist_ = dist; }
  void setMinDist(double dist) { min_traversal_dist_ = dist; }
  void setShortDists(double min, double max) { shortdist_min_ = min;  shortdist_max_ = max; }
  std::vector<RobotConfiguration *> generate(const RobotConfiguration &rc) const;

private:
  SkeletonPlanner();
  SkeletonPlanner(map_t *map);

  void generate(const std::vector<int> &bots, const std::vector<std::vector<int> > &neighbors,
                std::vector<int> *current,
                std::vector<RobotConfiguration *> *configs) const;
  std::vector<int> getNeighbors(int bot_ind) const;  
  std::set<int> getNeighbors(int bot_ind, double min, double max) const;
  std::vector<int> mergeNeighbors(const std::vector<int> &ns) const;
  
  std::vector<rf::Point2D> vertices_;
  std::vector<std::vector<int > > neighbors_;
  double min_traversal_dist_;
  double shortdist_min_, shortdist_max_;
  double merge_dist_;
};

class DynamicPlanner : public LatticePlanner<DynamicConfiguration> {
public:
    DynamicPlanner();
    DynamicPlanner(map_t *map_);

    std::vector<DynamicConfiguration *> generate(const DynamicConfiguration &rc) const;

    void setOmega(double step_min, double step_max, double step_size, double max) {
        assert(step_size > 0);
        assert(step_max > step_min);
        assert(max > 0);        
        omega_.step_min = step_min;
        omega_.step_max = step_max;
        omega_.step_size = step_size;
        omega_.max = max;
    }

    void setVel(double step_min, double step_max, double step_size, double max) {
        assert(step_size > 0);
        assert(step_max > step_min);
        assert(max > 0);
        vel_.step_min = step_min;
        vel_.step_max = step_max;
        vel_.step_size = step_size;
        vel_.max = max;
    }

    void setTimeStep(double dt) {
        assert(dt > 0);
        dt_ = dt;
    }
    
private:
    void generate(const DynamicConfiguration &rc,
                  int bot_ind,
                  DynamicConfiguration *current_config,
                  std::vector<DynamicConfiguration *> *configs) const;
    struct {
      double step_min, step_max, step_size, max;
    } vel_, omega_;
    // Time step for forward calculations
    double dt_;
};

#endif
