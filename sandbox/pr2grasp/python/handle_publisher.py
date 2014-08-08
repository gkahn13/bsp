#!/usr/bin/env python

import rospy
import roslib

roslib.load_manifest('tfx')
roslib.load_manifest('handle_detector') # make sure using hydro

import tfx
import handle_detector.msg as hd_msg
import geometry_msgs.msg as geom_msg

import numpy as np

class HandlePublisher:
    topic_name = '/localization/handle_list'
    def __init__(self):
        self.handle_list_sub = rospy.Subscriber(HandlePublisher.topic_name, hd_msg.HandleListMsg, self.callback)
        self.pose_pub = rospy.Publisher('/localization/handle_pose', geom_msg.PoseStamped)
        rospy.loginfo('Publishing most recent handle to {0}'.format(HandlePublisher.topic_name))
        
    def callback(self, msg):
        handles = list()
        for handle in msg.handles:
            for cylinder in handle.cylinders:
                handle_pose = tfx.pose(cylinder.pose, stamp=msg.header.stamp, frame=msg.header.frame_id)
    
                if handle_pose.orientation.matrix[2,2] > 0:
                    continue
                
                
                g = np.array([cylinder.axis.x, cylinder.axis.y, cylinder.axis.z])
                r = np.array([cylinder.normal.x, cylinder.normal.y, cylinder.normal.z])
                b = np.cross(r, g)
                
                g = -g
                b = -b
                
                
                rot = np.zeros((3,3))
                rot[:,0] = r / np.linalg.norm(r)
                rot[:,1] = g / np.linalg.norm(g)
                rot[:,2] = b / np.linalg.norm(b)
                
                handle_pose.orientation = tfx.tb_angles(np.dot(tfx.tb_angles(-90,0,0).matrix,rot))
        
                """
                rot = np.array(handle_pose.orientation.matrix)
                if rot[2,2] > 0:
                    # facing wrong way
                    new_rot = np.zeros((3,3))
                    
                    #new_rot[:3,0] = rot[:3,1]
                    #new_rot[:3,1] = rot[:3,0]
                    #new_rot[:3,2] = -rot[:3,2]
                    
                    new_rot[:3,0] = rot[:3,0]
                    new_rot[:3,1] = -rot[:3,2]
                    new_rot[:3,2] = rot[:3,1]
                    
                    handle_pose.orientation = tfx.tb_angles(new_rot)
                """
                
                #if rot[2,2] > 0:
                #    self.pose_pub.publish(handle_pose.msg.PoseStamped())
                #    return
                
                handles.append(handle_pose)
                
                #self.pose_pub.publish(handle_pose.msg.PoseStamped())
                #return
                
        if len(handles) == 0:
            return
                
        #handles_base = [tfx.convertToFrame(handle, 'base_link') for handle in handles]
        handles_base = list()
        for handle in handles:
            cam_to_tool_frame = tfx.lookupTransform('r_gripper_tool_frame', handle.frame)
            tool_frame_to_base = tfx.lookupTransform('base_link', 'r_gripper_tool_frame')
            handles_base.append(tool_frame_to_base*(cam_to_tool_frame*handle))
                    
        sorted_handles_base = sorted(handles_base, key=lambda h: h.position.y)
        
        self.pose_pub.publish(sorted_handles_base[-1].msg.PoseStamped())
            
if __name__ == '__main__':
    rospy.init_node('handle_publisher', anonymous=True)
    handle_publisher = HandlePublisher()
    rospy.spin()