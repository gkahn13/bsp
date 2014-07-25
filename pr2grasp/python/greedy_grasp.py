#!/usr/bin/env python
import numpy as np

import rospy
import roslib
roslib.load_manifest('tfx')
roslib.load_manifest('handle_detector') # make sure using hydro

import tfx
import handle_detector.msg as hd_msg
import geometry_msgs.msg as geom_msg

from pr2 import arm, planner, simulator, utils

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
        #rospy.loginfo('Received handle list! Frame: {0}'.format(msg.header.frame_id))
        self.last_time_received = msg.header.stamp
        self.last_msg = msg
        self.last_pose = tfx.pose(msg.handles[0].cylinders[0].pose, stamp=msg.header.stamp, frame=msg.header.frame_id)
        
        rot = np.array(self.last_pose.orientation.matrix)
        if rot[2,2] > 0:
            # facing wrong way
            new_rot = np.zeros((3,3))
            new_rot[:3,0] = rot[:3,1]
            new_rot[:3,1] = rot[:3,0]
            new_rot[:3,2] = -rot[:3,2]
            self.last_pose.orientation = tfx.tb_angles(new_rot) 
        
        self.pose_pub.publish(self.last_pose.msg.PoseStamped())
        
    def get_new_handle(self):
        if self.last_pose is None:
            while self.last_pose is None:
                rospy.sleep(.01)
        else: 
            last_stamp = self.last_pose.stamp
            while last_stamp == self.last_pose.stamp:
                rospy.sleep(.01)
            
        return self.last_pose


class GreedyGrasp:
    """ Grasp pipeline which plans towards handles in state-space """
    def __init__(self):
        self.hd = HandleDetector()
        self.arm = arm.Arm('right')
        
        self.arm.go_to_posture('mantis')
        self.home_pose = self.arm.get_pose()
        
        self.handle_pose = None
        self.handle_pose_pub = rospy.Publisher('/debug_handle_pose', geom_msg.PoseStamped)
        
        rospy.loginfo('GreedyGrasp initialized')
        
    def start(self):
        utils.press_enter_to_continue('start')
        
        rospy.loginfo('Starting GreedyGrasp pipeline')
        self.go_to_home()
        
    def go_to_home(self):
        utils.press_enter_to_continue('go_to_home')
        
        rospy.loginfo('Moving to home pose')
        self.arm.go_to_pose(self.home_pose)
        
        rospy.loginfo('Opening gripper')
        self.arm.open_gripper()
        
        self.get_handle()
    
    def get_handle(self):
        utils.press_enter_to_continue('get_handle')
        
        rospy.loginfo('Waiting for handle pose...')
        handle_pose_cam = self.hd.get_new_handle()
        #handle_pose_cam = tfx.pose([-0.069, -0.008,  0.471],tfx.tb_angles(169.5, 14.5, 63.1), frame='/camera_rgb_optical_frame') # TEMP stub
        print('handle_pose_cam frame: {0}'.format(handle_pose_cam.frame))
        self.handle_pose = tfx.lookupTransform('base_link', handle_pose_cam.frame)*handle_pose_cam
        self.handle_pose_pub.publish(self.handle_pose.msg.PoseStamped())
        rospy.loginfo('Received new handle pose')
    
        self.grab_handle()
    
    def grab_handle(self):
        utils.press_enter_to_continue('grab_handle')
        
        rospy.loginfo('Going to handle')
        self.arm.go_to_pose(self.handle_pose)
        
        rospy.loginfo('Closing gripper')
        self.arm.close_gripper()
        
        self.drop_off_object()
        
    def drop_off_object(self):
        utils.press_enter_to_continue('drop_off_object')
        
        rospy.loginfo('Moving vertical')
        current_pose = self.arm.get_pose()
        self.arm.go_to_pose(current_pose + [0,0,.05])
        
        rospy.loginfo('Moving to home pose')
        self.arm.go_to_pose(self.home_pose)
        
        rospy.loginfo('Opening gripper')
        self.arm.open_gripper()


def test_greedy_grasp():
    gg = GreedyGrasp()
    gg.start()
    
    IPython.embed()
    
def test_handle_detector():
    hd = HandleDetector()
    
    while hd.last_msg is None:
        rospy.sleep(.05)
        
    msg = hd.last_msg
    c = msg.handles[0].cylinders[0]
    
    IPython.embed()

def test_planner():
    arm_name='right'
    sim = simulator.Simulator(view=False)
    a = arm.Arm(arm_name, sim=sim)
    p = planner.Planner(arm_name, sim=sim, interact=True)
    
    a.go_to_posture('untucked', speed=.2)
    
    current_joints = a.get_joints()
    current_pose = a.get_pose()
    target_pose = current_pose + [.2,0,0]
    
    print('Calling trajopt')
    joint_traj = p.get_joint_trajectory(current_joints, target_pose)
    
    print('\ndesired end pose:\n{0}'.format(target_pose.matrix))
    print('joint_traj end pose:\n{0}'.format(a.fk(joint_traj[-1]).matrix))
    
    IPython.embed()

if __name__ == '__main__':
    rospy.init_node('greedy_grasp', anonymous=True)
    #test_greedy_grasp()
    #test_handle_detector()
    test_planner()