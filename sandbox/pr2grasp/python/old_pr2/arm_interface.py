#!/usr/bin/env python

"""
Note: environment loads a xml file for the table. Must change depending on what table you are using
"""

import openravepy
import trajoptpy
import json
import numpy as np
import trajoptpy.kin_utils as ku
import os
import random

import rospy
import openravepy as rave
from pr2 import pr2

import tf
import tf.transformations as tft

from constants import Constants
#from Util import *
#from Pr2Debridement.srv import ReturnJointStates

from geometry_msgs.msg import PointStamped, PoseStamped
from trajectory_msgs.msg import JointTrajectory, JointTrajectoryPoint

import code

import IPython

class ArmInterface:
    """
    Interface class for controlling an arm of the PR2
    """

    joint_name_suffixes = ["_shoulder_pan_joint",
                           "_shoulder_lift_joint",
                           "_upper_arm_roll_joint",
                           "_elbow_flex_joint",
                           "_forearm_roll_joint",
                           "_wrist_flex_joint",
                           "_wrist_roll_joint"]

    tool_frame_suffix = '_gripper_tool_frame'

    

    def __init__ (self, pr2, armName):
        self.pr2 = pr2      
        if armName == Constants.ArmName.Left:
            self.arm = self.pr2.larm
            self.armName = 'leftarm'
        else:
            self.arm = self.pr2.rarm
            self.armName = 'rightarm'

        self.joint_names = [self.armName[0]+suffix for suffix in self.joint_name_suffixes]
        self.toolframe = self.armName[0] + self.tool_frame_suffix

        self.listener = tf.TransformListener()

        # number of iterations for trajopt
        self.n_steps = 60
        
        # used to slow down the pr2 motions
        slow_down_ratio = .5

        self.pr2.robot.SetDOFVelocityLimits(slow_down_ratio*self.pr2.robot.GetDOFVelocityLimits())
        self.pr2.robot.SetDOFAccelerationLimits(slow_down_ratio*self.pr2.robot.GetDOFAccelerationLimits())
 
        # slow down velocites
        self.arm.vel_limits = np.array([slow_down_ratio*limit for limit in self.arm.vel_limits])
        # slow down acceleration
        self.arm.acc_limits = np.array([slow_down_ratio*limit for limit in self.arm.acc_limits])

        rospy.sleep(1)

 
    def goToArmPose(self, pose, isPlanned, reqName=Constants.Request.noRequest):
        """
        Go to PoseStamped pose

        If isPlanned is True, then trajopt is used to plan the trajectory
        Otherwise, just IK is used to plan trajectory

        If IK fails, arm doesn't move
        
        The methods for executing the trajectories are in pr2/pr2.py
        """
        # must convert to BaseLink frame
        if self.listener != None:
            try:
                commonTime = self.listener.getLatestCommonTime(Constants.BaseLink,pose.header.frame_id)
                pose.header.stamp = commonTime
                pose = self.listener.transformPose(Constants.BaseLink,pose)
            except tf.Exception:
                return
            
        if pose.header.frame_id != Constants.BaseLink:
            return

        # get the current joint state for joint_names
        rospy.wait_for_service("return_joint_states")
        s = rospy.ServiceProxy("return_joint_states", ReturnJointStates)
        resp = s(self.joint_names)
        
        # set the start joint position
        joint_start = resp.position
        self.robot.SetDOFValues(joint_start, self.robot.GetManipulator(self.armName).GetArmIndices())

        # initialize trajopt inputs
        quat = pose.pose.orientation
        xyz = pose.pose.position
        quat_target = [quat.w, quat.x, quat.y, quat.z]
        xyz_target = [xyz.x, xyz.y, xyz.z]
        hmat_target = openravepy.matrixFromPose( np.r_[quat_target, xyz_target] )
        
        # if no planning
        if not isPlanned:
            self.arm.goto_pose_matrix(hmat_target, Constants.BaseLink, self.toolframe)
            return True

        # inverse kinematics
        manip = self.robot.GetManipulator(self.armName)
        init_joint_target = ku.ik_for_link(hmat_target, manip, self.toolframe, filter_options = openravepy.IkFilterOptions.CheckEnvCollisions)
        
        if init_joint_target == None:
            # inverse kinematics failed
            # will do nothing for now, may want to alter xyz_target a little
            rospy.loginfo('IK failed')
            return False




        request = self.getRequest(reqName, xyz_target, quat_target, init_joint_target)
        
        if request == None:
            return
               
        # convert dictionary into json-formatted string
        s = json.dumps(request) 
        # create object that stores optimization problem
        prob = trajoptpy.ConstructProblem(s, self.env)
        # do optimization
        result = trajoptpy.OptimizeProblem(prob)

        self.arm.follow_joint_trajectory(result.GetTraj())

        return True

    def getRequest(self, reqName, xyz_target, quat_target, init_joint_target):
        """
        Different request types for trajopt costs/constraints. See Constants.Request
        """
        if reqName == Constants.Request.goReceptacle:
            
            request = {
            "basic_info" : {
                "n_steps" : self.n_steps,
                "manip" : self.armName, 
                "start_fixed" : True 
                },
            "costs" : [
                {
                    "type" : "joint_vel",
                    "params": {"coeffs" : [1]} 
                    },
                {
                    "type" : "collision",
                    "params" : {
                        "coeffs" : [20],
                        "dist_pen" : [0.025] 
                        }
                    },
                {
                        "type" : "pose",
                        "name" : "target_pose",
                        "params" : {"xyz" : xyz_target, 
                                    "wxyz" : quat_target,
                                    "link": self.toolframe,
                                    "rot_coeffs" : [1,1,0],
                                    "pos_coeffs" : [0,0,0],
                                    "coeffs" : [20]
                                    }
                        },
                ],
            "constraints" : [
                # BEGIN pose_target
                {
                    "type" : "pose",
                    "name" : "target_pose",
                    "params" : {"xyz" : xyz_target, 
                                "wxyz" : quat_target,
                                "link": self.toolframe,
                                "rot_coeffs" : [0,0,0],
                                "pos_coeffs" : [1,1,1]
                                }
                    
                    },
                #END pose_target
                ],
            "init_info" : {
                "type" : "straight_line",
                "endpoint" : init_joint_target.tolist()
                }
            }

        else:
            request = {
                "basic_info" : {
                    "n_steps" : self.n_steps,
                    "manip" : self.armName, 
                    "start_fixed" : True 
                    },
                "costs" : [
                    {
                        "type" : "joint_vel",
                        "params": {"coeffs" : [1]} 
                        },
                    {
                        "type" : "collision",
                        "params" : {
                            "coeffs" : [20],
                            "dist_pen" : [0.025] 
                            }
                        },
                    ],
                "constraints" : [
                    {
                        "type" : "pose",
                        "name" : "target_pose",
                        "params" : {"xyz" : xyz_target, 
                                    "wxyz" : quat_target,
                                    "link": self.toolframe,
                                    "rot_coeffs" : [1,1,1],
                                    "pos_coeffs" : [1,1,1]
                                    }
                        
                        },
                    ],

                "init_info" : {
                    "type" : "straight_line",
                    "endpoint" : init_joint_target.tolist()
                    }
                }
            

        return request



    def goToSide(self):
        """
        Commands the arm to the side

        Useful for debugging
        """
        self.arm.goto_posture('side')

def test():
    rospy.init_node('arm_interface', anonymous=True)
    brett = pr2.PR2()
    arm = ArmInterface(brett, Constants.ArmName.Right)
    IPython.embed()

if __name__ == '__main__':    
    test()
