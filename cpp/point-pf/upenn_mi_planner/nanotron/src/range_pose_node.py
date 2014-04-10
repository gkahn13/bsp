#!/usr/bin/python
# A node that synchronizes range and pose information from a robot.

# Copyright (C) 2012 Benjamin Charrow
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to
# deal in the Software without restriction, including without limitation the
# rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
# sell copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
# IN THE SOFTWARE.
#

import roslib; roslib.load_manifest('nanotron')
import rospy

import geometry_msgs.msg
import nanotron.msg

import threading

class Aggregator(object):
    def __init__(self):
        self._lock = threading.Lock()
        self._name = rospy.get_param("~name")
        self._pub = rospy.Publisher("range_pose", nanotron.msg.RangeWithPose)
        self._pose_sub = rospy.Subscriber("amcl_pose", geometry_msgs.msg.PoseWithCovarianceStamped,
                                          self._pose_callback)
        self._range_sub = rospy.Subscriber("range", nanotron.msg.Range,
                                           self._range_callback)
        self._pose = None
        self._prev_time = None
        
    def _pose_callback(self, pose_msg):
        with self._lock:
            curr_time = pose_msg.header.stamp
            if (self._prev_time is None) or (curr_time - self._prev_time >= rospy.Duration(0)):
                self._prev_time = curr_time
                self._pose = pose_msg
            else:
                rospy.logerr("Received pose from the past, ignoring")


    def _range_callback(self, range_msg):
        with self._lock:
            if self._pose is not None:
                msg = nanotron.msg.RangeWithPose(range = range_msg,
                                                 pose = self._pose,
                                                 id = self._name)
                self._pub.publish(msg)

def main():
    global pub
    rospy.init_node('range_pose_aggregator')

    agg = Aggregator()
    
    rospy.spin()

if __name__ == "__main__":
    main()
