#!/usr/bin/env python
import numpy as np

import rospy
import roslib
roslib.load_manifest('tfx')

import tfx
import geometry_msgs.msg as geom_msg

from pr2 import arm, planner, simulator, utils

import IPython

class GreedyGrasp:
    """ Grasp pipeline which plans towards handles in state-space """
    def __init__(self):
        self.arm_name = 'right'
        self.sim = simulator.Simulator(view=True)
        self.arm = arm.Arm(self.arm_name, sim=self.sim)
        self.planner = planner.Planner(self.arm_name, sim=self.sim, interact=False)
        
        self.sim.env.Load('../envs/SDHtable.xml')
        
        #self.arm.go_to_posture('mantis')
        self.home_pose = self.arm.get_pose()
        
        self.handle_pose_sub = rospy.Subscriber('/localization/handle_pose', geom_msg.PoseStamped, self._handle_pose_callback)
        self.handle_pose_callback, self.handle_pose = None, None
        self.handle_pose_pub = rospy.Publisher('/debug_handle_pose', geom_msg.PoseStamped)
        
        rospy.loginfo('GreedyGrasp initialized')
        
    def _handle_pose_callback(self, msg):
        self.handle_pose_callback = tfx.pose(msg)
        
    def start(self):
        utils.press_enter_to_continue('start')
        
        rospy.loginfo('Starting GreedyGrasp pipeline')
        self.go_to_home()
        
    def go_to_home(self):
        utils.press_enter_to_continue('go_to_home')
        
        rospy.loginfo('Moving to home pose')
        joint_traj = self.planner.get_joint_trajectory(self.arm.get_joints(), self.home_pose)
        self.arm.execute_joint_trajectory(joint_traj)
        #self.arm.go_to_pose(self.home_pose)
        
        rospy.loginfo('Opening gripper')
        self.arm.open_gripper()
        
        self.get_handle()
    
    def get_handle(self):
        utils.press_enter_to_continue('get_handle')
        
        rospy.loginfo('Waiting for handle pose...')
        
        last_handle_pose = self.handle_pose_callback
        while last_handle_pose.stamp == self.handle_pose_callback.stamp:
            rospy.sleep(.01)
        self.handle_pose = self.handle_pose_callback
        
        #cam_to_tool_frame = tfx.lookupTransform(self.arm.tool_frame, handle_pose_cam.frame)
        #tool_frame_to_base = tfx.lookupTransform('base_link', self.arm.tool_frame)
        #self.handle_pose = tool_frame_to_base*(cam_to_tool_frame*handle_pose_cam)
        self.handle_pose_pub.publish(self.handle_pose.msg.PoseStamped())
        rospy.loginfo('Received new handle pose')
        print(self.handle_pose)
    
        self.grab_handle()
    
    def grab_handle(self):
        utils.press_enter_to_continue('grab_handle')
        
        rospy.loginfo('Going to handle')
        #import trajoptpy
        #trajoptpy.SetInteractive(True)
        joint_traj = self.planner.get_joint_trajectory(self.arm.get_joints(), self.handle_pose)
        self.arm.execute_joint_trajectory(joint_traj)
        #self.arm.go_to_pose(self.handle_pose)
        
        rospy.loginfo('Closing gripper')
        self.arm.close_gripper()
        
        self.drop_off_object()
        
    def drop_off_object(self):
        utils.press_enter_to_continue('drop_off_object')
        
        rospy.loginfo('Moving vertical')
        current_pose = self.arm.get_pose()
        joint_traj = self.planner.get_joint_trajectory(self.arm.get_joints(), current_pose + [0,0,.05])
        self.arm.execute_joint_trajectory(joint_traj)
        #self.arm.go_to_pose(current_pose + [0,0,.05])
        
        rospy.loginfo('Moving to home pose')
        joint_traj = self.planner.get_joint_trajectory(self.arm.get_joints(), self.home_pose)
        self.arm.execute_joint_trajectory(joint_traj)
        #self.arm.go_to_pose(self.home_pose)
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
    p = planner.Planner(arm_name, sim=sim, interact=False)
    
    a.go_to_posture('untucked', speed=.2)
    
    current_joints = a.get_joints()
    current_pose = a.get_pose()
    target_pose = current_pose + [.2,0,0]
    
    print('Calling trajopt')
    joint_traj = p.get_joint_trajectory(current_joints, target_pose)
    
    print('\ndesired end pose:\n{0}'.format(target_pose.matrix))
    print('joint_traj end pose:\n{0}'.format(a.fk(joint_traj[-1]).matrix))
    
    IPython.embed()
    
def test_tf():
    tfx.lookupTransform('camera_rgb_optical_frame','r_gripper_tool_frame')
    IPython.embed()

if __name__ == '__main__':
    rospy.init_node('greedy_grasp', anonymous=True)
    test_greedy_grasp()
    #test_handle_detector()
    #test_planner()
    #test_tf()