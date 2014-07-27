#!/usr/bin/env python

import rospy
import roslib
import rosparam
roslib.load_manifest('sensor_msgs')
roslib.load_manifest('tfx')
import sensor_msgs.msg as sm
import tfx
import tf

import os
import subprocess
import signal
import argparse

import IPython

class SavePCs:
    def __init__(self, topic_name):
        rospy.loginfo('Listening for point cloud topic: {0}'.format(topic_name))
        self.topic_name = topic_name
        
        rospy.loginfo('Starting pcl_ros pointcloud_to_pcd')
        self.pub_name = 'save_pcs_points'
        self.pointcloud_to_pcd_pub = rospy.Publisher(self.pub_name, sm.PointCloud2)
        
        self.pc_sub = rospy.Subscriber(topic_name, sm.PointCloud2, self._pc_callback)
        self.last_pc = None
        
        self.num_saved = 0
        
        rospy.loginfo('Waiting for first point cloud...')
        while self.last_pc is None and not rospy.is_shutdown():
            rospy.sleep(.05)
        
        rospy.loginfo('Ready!')
        
        
    def run(self):
        rospy.loginfo('Press enter to save point cloud: {0}'.format(self.topic_name))
        
        while not rospy.is_shutdown():
            rospy.loginfo('Press enter')
            raw_input()
            
            pc = self.last_pc
            pc.header.stamp = rospy.Duration(self.num_saved)
            self.pointcloud_to_pcd_pub.publish(self.last_pc)
            
            rospy.loginfo('Publish point-cloud {0}'.format(self.num_saved))
            self.num_saved += 1            
        
        
    def _pc_callback(self, msg):
        self.last_pc = msg

if __name__ == '__main__':
    rospy.init_node('save_pcs', anonymous=True)
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--topic',type=str,default='/camera/depth_registered/points')
    args = parser.parse_args()
    
    save_pcs = SavePCs(args.topic)
    save_pcs.run()
    