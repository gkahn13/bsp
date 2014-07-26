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

#############
# IMPORTANT #
#############
""" Make sure to run 'rosrun pcl_ros pointcloud_to_pcd input:={file}' beforehand """

class SavePCs:
    def __init__(self, topic_name, file_name):
        rospy.loginfo('Listening for point cloud topic: {0}'.format(topic_name))
        self.topic_name = topic_name
        self.file_name = file_name
        
        taubin_dir = os.path.abspath('.')
        if not os.path.exists(self.file_name):
            os.mkdir(self.file_name)
            
        self.full_file_name = '{0}/{1}/{2}'.format(taubin_dir, self.file_name, self.file_name)
        
        self.pc_sub = rospy.Subscriber(topic_name, sm.PointCloud2, self._pc_callback)
        self.last_pc = None
        
        self.listener = tf.TransformListener()
        self.num_saved = 0
        
        rospy.loginfo('Waiting for first point cloud...')
        while self.last_pc is None and not rospy.is_shutdown():
            rospy.sleep(.05)
            
        rospy.loginfo('Starting pcl_ros pointcloud_to_pcd')
        self.pub_name = '/pointcloud_to_pcd_input'
        self.pointcloud_to_pcd_pub = rospy.Publisher(self.pub_name, sm.PointCloud2)
        #os.system('rosrun pcl_ros pointcloud_to_pcd input:=' + pub_name)
        #subprocess.Popen('rosrun pcl_ros pointcloud_to_pcd input:=' + self.pub_name, shell=True)
        
        rospy.loginfo('Ready!')
        
        
    def run(self):
        rospy.loginfo('Press enter to save point cloud: {0}'.format(self.topic_name))
        
        while not rospy.is_shutdown():
            rospy.loginfo('Press enter')
            raw_input()
            
            file_name_i = self.full_file_name + str(self.num_saved) + '_'
            self.num_saved += 1
            
            #rosparam.set_param('prefix',file_name_i)
            pro = subprocess.Popen('rosrun pcl_ros pointcloud_to_pcd input:={0} _prefix:={1}'.format(self.pub_name, file_name_i), stdout=subprocess.PIPE, 
                       shell=True, preexec_fn=os.setsid)
            rospy.sleep(1)
            self.pointcloud_to_pcd_pub.publish(self.last_pc)
            rospy.sleep(.5)
            
            os.killpg(pro.pid, signal.SIGTERM)
            
            rospy.loginfo('Saved to file {0}'.format(file_name_i))
            
        
        
    def _pc_callback(self, msg):
        self.last_pc = msg

if __name__ == '__main__':
    rospy.init_node('save_pcs', anonymous=True)
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--topic',type=str,default='/camera/depth_registered/points')
    parser.add_argument('--file',type=str,default='pc')
    args = parser.parse_args()
    
    save_pcs = SavePCs(args.topic, args.file)
    save_pcs.run()
    