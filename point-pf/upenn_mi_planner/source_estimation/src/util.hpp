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

#ifndef UTIL_HPP
#define UTIL_HPP

#include <ros/ros.h>

#include "gmm/gmm.hpp"

#include "particles.hpp"

namespace rf {
  template<typename P>
  MeasurementModel<P>* getMeasurementModel(const ros::NodeHandle &nh,
                                           MeasurementAggregator *agg,
                                           LOSChecker *checker) {
    std::string type;
    nh.param<std::string>("measurement_model", type, "none");

    ROS_INFO("Using '%s' measurement model", type.c_str());
    if (type == "gaussian") {
      return GaussFading<P>::FromROS(nh, agg);
    } else if (type == "fading") {
      return LosNlosFading<P>::FromROS(nh, checker, agg);
    }

    ROS_ERROR("getMeasurementModel(): Unrecognized model");
    ROS_BREAK();
    return NULL; // Prevent gcc from complaining about not returning anything
  }
}
#endif
