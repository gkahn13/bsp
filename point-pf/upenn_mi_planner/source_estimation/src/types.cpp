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

#include "types.hpp"

#include <cstdio>

namespace rf {
  double l2_squared(const Point2D &p1, const Pose2D &p2) {
    double dx = p1.x - p2.x();
    double dy = p1.y - p2.y();
    return dx * dx + dy * dy;
  }

  double l2_squared(const Point2D &p1, const Point2D &p2) {
    double dx = p1.x - p2.x;
    double dy = p1.y - p2.y;
    return dx * dx + dy * dy;
  }


  double wrapToPi(double angle) {
    while (angle >= M_PI) {
      angle -= 2 * M_PI;
    }
    while (angle <= -M_PI) {
      angle += 2 * M_PI;
    }
    return angle;
  }

  std::ostream& operator<<(std::ostream &s, const Pose2D &p) {
    int prev_prec = s.precision(2);
    s << "x: " << p.x() << " y: " << p.y() << " t: " << p.theta();
    s.precision(prev_prec);
    return s;
  }

  int findFirstGreater(double *cdf, int sz, double val) {
    int low = 1, high = sz;
    int curr;

    if (val <= cdf[0]) {
      return 0;
    }
    while (low <= high) {
      curr = (low + high) / 2;
      bool bigger_prev = cdf[curr - 1] < val;
      bool smaller_next = val <= cdf[curr];

      if (bigger_prev && smaller_next) {
        return curr;
      } else if (bigger_prev) { // Must be bigger than cdf[curr]
        low = curr + 1;
      } else { // Must be smaller than cdf[curr - 1]
        high = curr - 1;
      }
    }

    std::cerr << "PeriodicResampler::findBinFirstGreater() " <<
      "val = " << val << " cdf[end] = " << cdf[sz - 1] << std::endl;
    return sz;
  }
}
