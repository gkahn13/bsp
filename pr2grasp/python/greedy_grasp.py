#!/usr/bin/env python

import rospy
import roslib
roslib.load_manifest('tfx')
roslib.load_manifest('handle_detector') # make sure using hydro

import tfx
import handle_detector.msg as hd_msg
import geometry_msgs.msg as geom_msg



import IPython

class HandleDetector:
    topic_name = '/localization/handle_list'
    def __init__(self):
        self.handle_list_sub = rospy.Subscriber(HandleDetector.topic_name, hd_msg.HandleListMsg, self.callback)
        self.pose_pub = rospy.Publisher('/localization/handle_pose', geom_msg.PoseStamped)
        
        self.last_time_received = rospy.Time.now()
        self.last_msg = None
        self.last_pose = None
        
    def callback(self, msg):
        rospy.loginfo('Received handle list! Frame: {0}'.format(msg.header.frame_id))
        self.last_time_received = msg.header.stamp
        self.last_msg = msg
        self.last_pose = tfx.pose(msg.handles[0].cylinders[0].pose, frame=msg.header.frame_id)
        
        self.pose_pub.publish(self.last_pose.msg.PoseStamped())

def test_handle_detector():
    hd = HandleDetector()
    
    msg = hd.last_msg
    while not rospy.is_shutdown():
        if msg != hd.last_msg:
            msg = hd.last_msg
            #hd.plot_handle_list_msg(msg)
            
        rospy.sleep(.05)
        
    IPython.embed()

if __name__ == '__main__':
    rospy.init_node('greedy_grasp', anonymous=True)
    rospy.sleep(1)
    test_handle_detector()