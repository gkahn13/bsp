#!/usr/bin/env python

"""
For controlling the gripper
"""

# Import required Python code.
import roslib
roslib.load_manifest('Pr2Debridement')
import rospy
import sys
import actionlib
from pr2_controllers_msgs.msg import Pr2GripperCommandAction, Pr2GripperCommandGoal, Pr2GripperCommand
from geometry_msgs.msg import PoseStamped

from Constants import ConstantsClass

class GripperControl():
    """
    Class for controlling the opening and closing of the gripper,
    as well as information related to the gripper
    """
    def __init__(self, gripperName):
        """
        gripperName must be from ConstantsClass.GripperName
        """

        if gripperName == ConstantsClass.GripperName.Left:
            self.toolframe = ConstantsClass.ToolFrame.Left
        else:
            self.toolframe = ConstantsClass.ToolFrame.Right

        # max gripper range spec is 90mm
        self.maxRange = .09
        # no effort limit (may change)
        self.effortLimit = -1
        # distance from center of gripper to tip
        self.gripperLength = .05

        self.client = actionlib.SimpleActionClient(gripperName + '_controller/gripper_action', Pr2GripperCommandAction)


    def openGripper(self):
        """
        Opens gripper all the way
        """
        return self.setGripper(1)

    def closeGripper(self):
        """
        Closes gripper all the way
        """
        return self.setGripper(0)

    def setGripper(self, openPercentage):
        """
        Opens the gripper openPercentage of the way
        """
        position = openPercentage * self.maxRange

        #self.client.wait_for_server()
        self.client.send_goal(Pr2GripperCommandGoal(Pr2GripperCommand(position, self.effortLimit)))
        self.client.wait_for_result()
        return True
        
        # below not guaranteed to work for grasping
        # result = self.client.get_result()        
        # return (result.reached_goal) or (result.stalled)

    def gripperLength(self):
        return self.gripperLength

    def gripperPose(self):
        """
        Returns a PoseStamped of current gripper pose
        according to tf
        """
        pose = PoseStamped()
        pose.header.stamp = rospy.Time.now()
        pose.header.frame_id = self.toolframe
        pose.pose.orientation.w = 1
        return pose





def test():
    """
    Opens and closes each gripper
    and checks for return status
    """

    success = True
    timeDelay = 3

    cgl = CommandGripperClass(ConstantsClass.GripperName.Left)
    success &= cgl.openGripper()
    rospy.sleep(timeDelay)
    success &= cgl.closeGripper()
    rospy.sleep(timeDelay)

    success &= cgl.setGripper(.5)
    rospy.sleep(timeDelay)
    success &= cgl.setGripper(.25)
    rospy.sleep(timeDelay)
    success &= cgl.setGripper(.75)
    rospy.sleep(timeDelay)

    success &= cgl.closeGripper()

    
    cgr = CommandGripperClass(ConstantsClass.GripperName.Right)    
    success &= cgr.openGripper()
    rospy.sleep(timeDelay)
    success &= cgr.closeGripper()
    rospy.sleep(timeDelay)

    success &= cgr.setGripper(.5)
    rospy.sleep(timeDelay)
    success &= cgr.setGripper(.25)
    rospy.sleep(timeDelay)
    success &= cgr.setGripper(.75)
    rospy.sleep(timeDelay)

    success &= cgr.closeGripper()


    if success:
        print('CommandGripperClass passed tests')
    else:
        print('CommandGripperClass failed!!!')
    

#for testing
if __name__ == '__main__':
    rospy.init_node('gripper_node')
    test()
