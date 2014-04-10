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

#ifndef TYPES_HPP
#define TYPES_HPP
#include <iostream>
#include <cmath>

#include <geometry_msgs/PoseArray.h>

// A macro to disallow the copy constructor and operator= functions
// This should be used in the private: declarations for a class
#define DISALLOW_COPY_AND_ASSIGN(TypeName)      \
  TypeName(const TypeName&);                    \
  void operator=(const TypeName&)

#define DEBUG(x) std::cout << #x << "= " << x << std::endl;

namespace rf {
  class Point2D {
  public:
    Point2D() : x(0), y(0) {};
    Point2D(double xx, double yy) : x(xx), y(yy) { }
    Point2D(const Point2D &o) : x(o.x), y(o.y) { }
  
    double x, y;
  };

  struct Velocity2D {
    double v, w;
  };

  class Pose2D {
  public:
    Pose2D() {
      setXYT(0, 0, 0);
    };
        
    Pose2D(double x, double y, double theta) {
      setXYT(x, y, theta);
    }

    double distance_sq(const Pose2D &other) const {
      double x_diff = (x() - other.x());
      double y_diff = (y() - other.y());
      return x_diff * x_diff + y_diff * y_diff;
    }
    
    double distance(const Pose2D &other) const {
      return sqrt(distance_sq(other));
    }

    double angleDiff(const Pose2D &other) const {
      throw new std::exception();
    }

    void setXYT(double x, double y, double t) {
      m_[0] = x;
      m_[1] = y;
      m_[2] = t;
    }

    void setXYT(const Pose2D &p) {
      setXYT(p.x(), p.y(), p.theta());
    }
        
    double x() const {
      return m_[0];
    }

    double y() const {
      return m_[1];
    }

    double theta() const {
      return m_[2];
    }

    void x(double x) {
      m_[0] = x;
    }
    void y(double y) {
      m_[1] = y;
    }
    void theta(double t) {
      m_[2] = t;
    }
  private:
    DISALLOW_COPY_AND_ASSIGN(Pose2D);
        
    double m_[3];
  };

  std::ostream& operator<<(std::ostream &s, const Pose2D &p);

  double l2_squared(const Point2D &p1, const Pose2D &p2);
  double l2_squared(const Point2D &p1, const Point2D &p2);
  double l2_squared(const Point2D &p1, const Point2D &p2);

  double wrapToPi(double angle);  

  /** 
   * Return first index in cdf whose value is greather than val.
   * 
   * @param cdf The array.
   * @param sz Size of the array.
   * @param val Value to search for.
   * 
   * @return Index in array.
   */
  int findFirstGreater(double *cdf, int sz, double val);
  
  //=============================== Particles ===============================//
  struct StationaryParticle {
    Point2D point;
    bool los;
    double weight;
  };
  std::ostream& operator<<(std::ostream &s, const StationaryParticle &sp);

  template <typename T=rf::StationaryParticle>
  class ParticleArray {
  public:
    ParticleArray(int n) {
      size = 0;
      ps = NULL;
      resize(n);
    }

    ParticleArray* Copy() {
      ParticleArray<T> *pa = new ParticleArray();
      pa->resize(size);
      memcpy(pa->ps, ps, sizeof(T) * size);
      return pa;
    }

    ParticleArray() {
      size = 0;
      ps = NULL;
    }

    ~ParticleArray() {
      delete[] ps;
    }

    /** 
     * Normalize the weights to sum to 1
     * 
     */
    void normalize(void) {
      double sum = 0;
      for (int i = 0; i < size; ++i) {
        sum += ps[i].weight;
      }
      double new_sum = 0;
      // ROS_WARN("weight sum: %f", sum);
      for (int i = 0; i < size; ++i) {
        ps[i].weight /= sum;
        new_sum += ps[i].weight;
      }
      ROS_ASSERT(fabs(new_sum -1) < 1e-6);
    }
    
    void update(const geometry_msgs::PoseArray::ConstPtr &msg) {
      if (msg->poses.size() != static_cast<size_t>(size)) {
        delete[] ps;
        size = msg->poses.size();
        ps = new T[size];
      }
      double sum_weight = 0;
      for (size_t i = 0; i < msg->poses.size(); ++i) {
        ps[i].point.x = msg->poses.at(i).position.x;
        ps[i].point.y = msg->poses.at(i).position.y;
        ps[i].weight = msg->poses.at(i).position.z;
        sum_weight += ps[i].weight;
      }
      if (abs(sum_weight - 1.0) > 1e-10) {
        ROS_ERROR("Particle weights don't sum to one");
      }
    }

    void resize(int n) {
      size = n;
      delete[] ps;
      ps = new T[n];
    }

    ParticleArray *FilterCDF(double pct) {
      // Sort particles by weight
      typedef std::pair<double, std::pair<double, double> > w_pt;
      std::vector<w_pt> particles;
      for (int ind = 0; ind < size; ++ind) {
        particles.push_back(std::make_pair(ps[ind].weight,
                                           std::make_pair(ps[ind].point.x,
                                                          ps[ind].point.y)));
      }
      sort(particles.begin(), particles.end());

      double cdf = 0;
      // Find first point that is greater than pct
      int start = 0;
      for (int ind = 0; ind < size; ++ind) {
        if (cdf + particles[ind].first >= pct) {
          break;
        }
        cdf += particles[ind].first;
        ++start;
      }

      // Build new particle array
      ParticleArray *pa = new ParticleArray();
      pa->resize(size - start);
      for (int ind = start; ind < size; ++ind) {
        pa->ps[ind - start].weight = particles[ind].first / (1 - cdf);
        pa->ps[ind - start].point.x = particles[ind].second.first;
        pa->ps[ind - start].point.y = particles[ind].second.second;
      }

      return pa;
    }

    /** 
     * Calculate the number of effective particles.
     * 
     * @return 
     */
    double neff() {
      double sum = 0;
      double sq_sum = 0;
      for (int i = 0; i < size; ++i) {
        sum += ps[i].weight;
        sq_sum += ps[i].weight * ps[i].weight;
      }
      return sum * sum / sq_sum;
    }
    
    /** 
     * Build a CDF of the particles.
     *
     * Caller should free the array via delete[].
     *
     * @param cdf Where to store the CDF
     * @param sum Cumulative sum of weights (particles may not be normalized)
     */
    void BuildCDF(double **cdf, double *sum) {
      *cdf = new double[size];
      *sum = ps[0].weight;

      (*cdf)[0] = ps[0].weight;
      for (int i = 1; i < size; ++i) {
        (*cdf)[i] = (*cdf)[i - 1] + ps[i].weight;
        *sum += ps[i].weight;
      }
    }

    T *ps;
    int size;
  private:
    DISALLOW_COPY_AND_ASSIGN(ParticleArray<T>);
  };
}



#endif
