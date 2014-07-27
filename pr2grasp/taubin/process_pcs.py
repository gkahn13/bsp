#!/usr/bin/env python

import rospy, roslib, rosparam, tf
roslib.load_manifest('sensor_msgs')
import sensor_msgs.msg as sm

import os
import sys
import subprocess
import signal
import argparse

import IPython

class HandleDetectorPCD:
    def __init__(self, file_path, **kwargs):
        self.params = dict()
        
        self.params['file'] = file_path
        self.params['target_radius'] = 0.025
        self.params['target_radius_error'] = 0.025
        self.params['affordance_gap'] = 0.008
        self.params['sample_size'] = 5000
        self.params['use_clearance_filter'] = 'true'
        self.params['use_occlusion_filter'] = 'true'
        self.params['curvature_estimator'] = 0
        self.params['point_cloud_source'] = 1
        self.params['update_interval'] = 1.0
        
        # RANSAC parameters
        self.params['ransac_runs'] = 5
        self.params['ransac_min_inliers'] = 8
        self.params['ransac_dist_radius'] = 0.02
        self.params['ransac_orient_radius'] = 0.4
        self.params['ransac_radius_radius'] = 0.003
        
        # workspace limits
        self.params['max_range'] = 0.9
        self.params['workspace_min_x'] = -0.25
        self.params['workspace_max_x'] = 0.45
        self.params['workspace_min_y'] = -0.2
        self.params['workspace_max_y'] = 0.4
        self.params['workspace_min_z'] = 0.3
        self.params['workspace_max_z'] = 1.0
        
        self.params['num_threads'] = 6
        
        for key, value in kwargs:
            if self.params.count(key) == 1:
                self.params[key] = value
        
        self.pro = None

        """        
        <param name="file" value="" />
        <param name="target_radius" value="0.025" /> <!--.03-->
        <param name="target_radius_error" value="0.025" /> <!--.025-->
        <param name="affordance_gap" value="0.008" />
        <param name="sample_size" value="5000" /> <!--20000-->
        <param name="use_clearance_filter" value="true" />
        <param name="use_occlusion_filter" value="true" />
        <param name="curvature_estimator" value="0" /> <!--0-->
        <param name="point_cloud_source" value="1" />
        <param name="update_interval" value="1.0" />
        
        <!-- RANSAC parameters -->
        <param name="ransac_runs" value="5" />
        <param name="ransac_min_inliers" value="8" />
        <param name="ransac_dist_radius" value="0.02" />
        <param name="ransac_orient_radius" value="0.4" />
        <param name="ransac_radius_radius" value="0.003" />
                
        <!-- workspace limits -->
        <param name="max_range" value="0.9" />
        <param name="workspace_min_x" value="-0.25" />
        <param name="workspace_max_x" value="0.45" />
        <param name="workspace_min_y" value="-0.2" />
        <param name="workspace_max_y" value="0.4" />
        <param name="workspace_min_z" value="0.3" />
        <param name="workspace_max_z" value="1.0" />
        
        <!-- number of threads to use -->
        <param name="num_threads" value="6" />
        """

    def run(self):
        rosrun_str = 'rosrun handle_detector handle_detector_localization'
        for param, value in self.params.items():
            rosrun_str += ' _{0}:={1}'.format(param, str(value))
            
        self.pro = subprocess.Popen(rosrun_str, stdout=subprocess.PIPE, shell=True, preexec_fn=os.setsid)
        
    def stop(self):
        os.killpg(self.pro.pid, signal.SIGTERM)

class ProcessPCs:
    def __init__(self, folder_name):
        taubin_dir = os.path.abspath('.')
        self.full_folder_name = '{0}/{1}'.format(taubin_dir, folder_name)
        
        if not os.path.exists(self.full_folder_name):
            rospy.logerr('Folder "{0}" does not exist, exiting'.format(self.full_folder_name))
            sys.exit()
        
        rospy.loginfo('Processing point-clouds in folder {0}'.format(self.full_folder_name))
        
        self.pcd_files = sorted([file for file in os.listdir(self.full_folder_name) if file[-4:]=='.pcd'])
        
        self.tf_broadcaster = tf.TransformBroadcaster()
        
        self.last_point_cloud = None
        self.cloud_sub = rospy.Subscriber('/localization/point_cloud', sm.PointCloud2, self._cloud_callback)
        self.cloud_pub = rospy.Publisher('/process_pcs/point_cloud', sm.PointCloud2)
        
        for i, pcd_file in enumerate(self.pcd_files):
            rospy.loginfo('Press enter to view pcd file {0}'.format(i))
            raw_input()
            
            hd_pcd = HandleDetectorPCD(pcd_file)
            hd_pcd.run()
            
            rospy.loginfo('Waiting for pcd file {0}'.format(i))
            while self.last_point_cloud is None and not rospy.is_shutdown():
                rospy.sleep(.01)
                
            #hd_pcd.stop()
            self.last_point_cloud = None
            break
            
        
        rospy.loginfo('Ready!')
        
    def _cloud_callback(self, msg):
        print('_cloud_callback')
        self.last_point_cloud = msg
        self.last_point_cloud.header.frame_id = '/base_link'
        self.cloud_pub.publish(self.last_point_cloud)
        
        self.tf_broadcaster.sendTransform((0, 0, 0),
                     [1, 0, 0, 0],
                     rospy.Time.now(),
                     '/camera_rgb_optical_frame',
                     '/base_link')


if __name__ == '__main__':
    rospy.init_node('process_pcs', anonymous=True)
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--folder',type=str,default='test')
    args = parser.parse_args()
    
    process_pcs = ProcessPCs(args.folder)
    rospy.spin()