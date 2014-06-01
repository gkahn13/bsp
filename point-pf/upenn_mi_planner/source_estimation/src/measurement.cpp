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

#include "measurement.hpp"

namespace rf {
  std::ostream& operator<<(std::ostream &s, const Measurement &m) {
    s << std::setprecision(15);
    s << "Measurement: " << std::endl <<
      " timestamp: " << m.timestamp() << std::endl << 
      " source_id:  " << m.sourceID() << std::endl <<
      " dest_id: " << m.destID() << std::endl <<
      " recv_xyt: " << m.sourceLoc() << std::endl <<
      " dist: " << std::setprecision(2) << m.dist() <<
      " var: " << m.variance() << 
      " dur: " << m.duration();
    return s;
  }
  
  //=============================== Utilities ===============================//
  void meanAndVar(std::list<double> msmts, double *mean, double *var) {
    std::list<double>::iterator it;
    double sum = 0;
    for (it = msmts.begin(); it != msmts.end(); ++it) {
      sum += (*it);
    }
    *mean = sum / msmts.size();
    sum = 0;
    for (it = msmts.begin(); it != msmts.end(); ++it) {
      sum += pow((*it) - *mean, 2.0);
    }
    if (msmts.size() == 1) {
      *var = 0;
    } else {
      *var = sum / ( msmts.size() - 1);
    }
  }

  void meanAndVar(const std::list<Measurement*> &msmts, double *mean, double *var) {
    std::list<double> data;
    std::list<double>::iterator it;
    for (std::list<Measurement*>::const_iterator it = msmts.begin(); it != msmts.end(); ++it) {
      data.push_back((*it)->dist());
    }
    meanAndVar(data, mean, var);
  }  
  //============================ NullAggregator =============================//

  NullAggregator::NullAggregator() : msmt_(NULL) {}
        
  bool NullAggregator::addMeasurement(const Measurement& m) {
    msmt_ = &m;
    return true;
  }
        
  const Measurement* NullAggregator::getAggregate() const {
    return msmt_;
  }

  //============================= STAggregator ==============================//
  STAggregator::STAggregator(double secs, double dist) :
    agg_msmt_(NULL), valid_(false), secs_(secs), dist_(dist),
    max_age_(std::numeric_limits<double>::infinity())  {}
    
  STAggregator::STAggregator(double secs, double dist, double max_age) :
    agg_msmt_(NULL), valid_(false), secs_(secs), dist_(dist), max_age_(max_age)  {}

  STAggregator::~STAggregator() {
    std::map<std::pair<int, int>, std::list<Measurement*> *>::iterator it;
    for (it = source_msmts_.begin(); it != source_msmts_.end(); ++it) {
      std::list<Measurement *> *msmts = (*it).second;

      std::list<Measurement *>::iterator list_it;
      for (list_it = msmts->begin(); list_it != msmts->end(); ++list_it) {
        delete *list_it;
      }      
      delete (*it).second;
    }

    for (std::vector<Measurement *>::iterator list_it = orig_msmts_.begin();
         list_it != orig_msmts_.end(); ++list_it) {
      delete *list_it;
    }
    
           
    delete agg_msmt_;
  }

  bool STAggregator::addMeasurement(const Measurement& msmt) {
    std::pair<int, int> key = std::make_pair(msmt.sourceID(), msmt.destID());
    std::list<Measurement*> *msmts = source_msmts_[key];
    if (msmts == NULL) {
      source_msmts_[key] = new std::list<Measurement *>();
      msmts = source_msmts_[key];
    }

    if (!msmts->empty() && msmt.timestamp() < msmts->back()->timestamp()) {
      ROS_WARN("STAggregator(): Measurement arrived out of order, skipping");
      valid_ = false;
      return valid_;
    }
      
    // Clear old measurements
    while (!msmts->empty() &&
           (msmt.timestamp() - msmts->front()->timestamp() > max_age_)) {
      ROS_DEBUG_THROTTLE(1.0, "STAggregator(): Removing old entry");
      delete msmts->front();
      msmts->pop_front();
    }
    // ROS_WARN("size: %i", msmts->size());
      
    if (msmts->empty()) {
      msmts->push_front(msmt.Copy());
      valid_ = false;
    } else {
      const Measurement *front = msmts->front();
      bool old = msmt.timestamp() - front->timestamp() >= secs_;
      bool source_move = front->sourceLoc().distance(msmt.sourceLoc()) >= dist_;
      // ROS_WARN("ts_delta: %0.3f -> %i", msmt.timestamp() - front->timestamp());
      if (old && source_move) {
        valid_= true;
        // Calculate mean distance and pose to use for distance calculation
        double mean, var;
        meanAndVar(*msmts, &mean, &var);
        delete agg_msmt_;
        double elapsed = msmts->back()->timestamp() - front->timestamp();
        agg_msmt_ = new Measurement(front->timestamp(), front->destID(),
                                    front->sourceID(), front->sourceLoc(),
                                    mean, var, elapsed);
        // Release memory from old aggregated measurement and copy over new
        // measurements
        for (std::vector<Measurement*>::iterator it = orig_msmts_.begin();
             it != orig_msmts_.end(); ++it) {
          delete *it;
        }
        orig_msmts_.clear();
        for (std::list<Measurement*>::iterator it = msmts->begin();
             it != msmts->end(); ++it) {
          orig_msmts_.push_back(*it);
        }
        
        ROS_DEBUG_STREAM("Aggregated " << msmts->size() << " msmsts:" << *agg_msmt_);
        // Add new measurement which is not part of previous group
        msmts->clear();
        msmts->push_back(msmt.Copy());
      } else {
        valid_ = false;
        msmts->push_back(msmt.Copy());
      }
    }
    return valid_;
  }
  
  const Measurement* STAggregator::getAggregate() const {
    if (!valid_) {
      throw std::runtime_error("STAggregator::getAggregate() Invalid");
    }
    return agg_msmt_;
  }  

  const std::vector<Measurement*>& STAggregator::getOrig() const {
    if (!valid_) {
      throw std::runtime_error("STAggregator::getOrig() Invalid");
    }
    return orig_msmts_;
  }  

}
