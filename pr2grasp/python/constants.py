#!/usr/bin/env python

"""
Contains general constants that
are frequently used
"""

class Constants:
    class ArmName:
        Left = 'l_arm'
        Right = 'r_arm'

    class GripperName:
        Left = 'l_gripper'
        Right = 'r_gripper'

    class ToolFrame:
        Left = 'l_gripper_tool_frame'
        Right = 'r_gripper_tool_frame'

    class Request:
        """
        Type of constraints/costs, for trajopt
        """
        noRequest = 0
        goNear = 1
        goReceptacle = 2

    BaseLink = 'base_link'
