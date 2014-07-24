import numpy as np

import rospy
import roslib
import actionlib
roslib.load_manifest("pr2_controllers_msgs")
roslib.load_manifest("move_base_msgs")
import trajectory_msgs.msg as tm
import geometry_msgs.msg as gm
import sensor_msgs.msg as sm
import pr2_controllers_msgs.msg as pcm


roslib.load_manifest('tfx')
import tfx
import tf.transformations as tft

import openravepy as rave # only for IK

import rave_utils
import utils

import IPython

class Arm:
    """
    Class for controlling the end effectors of the PR2
    """
    
    joint_name_suffixes = ["_shoulder_pan_joint",
                           "_shoulder_lift_joint",
                           "_upper_arm_roll_joint",
                           "_elbow_flex_joint",
                           "_forearm_roll_joint",
                           "_wrist_flex_joint",
                           "_wrist_roll_joint"]
    
    
    def __init__(self, arm_name, default_speed=.05):
        """
        :param arm_name: "left" or "right"
        :param default_speed: speed in meters/second
        """
        self.arm_name = arm_name
        self.default_speed = default_speed
        self.tool_frame = '{0}_gripper_tool_frame'.format(arm_name[0])
        self.joint_names = ['{0}{1}'.format(arm_name[0], joint_name) for joint_name in self.joint_name_suffixes]
        self.gripper_joint_name = '{0}_gripper_joint'.format(arm_name[0])
        self.min_grasp, self.max_grasp = -.01, .08
        self.default_max_effort = 40 # Newtons
        
        self.current_joints, self.current_grasp = None, None
        self.joint_state_sub = rospy.Subscriber('/joint_states', sm.JointState, self._joint_state_callback)
        
        rospy.loginfo('Waiting for /joint_states...')
        while self.current_joints is None or self.current_grasp is None:
            rospy.sleep(.01)
        
        self.joint_command_pub = rospy.Publisher('{0}_arm_controller/command'.format(arm_name[0]), tm.JointTrajectory)
        #self.gripper_command_pub = rospy.Publisher('{0}_gripper_controller/command'.format(arm_name[0]), pcm.Pr2GripperCommand)
        #self.gripper_command_pub = rospy.Publisher('{0}_gripper_traj_diagnostic'.format(arm_name[0]), tm.JointTrajectory)
        self.gripper_command_client = actionlib.SimpleActionClient('{0}_gripper_controller/gripper_action'.format(arm_name[0]), pcm.Pr2GripperCommandAction)
        
        # initialize openrave for IK
        self.env = rave.Environment()
        self.env.StopSimulation()
        self.env.Load("robots/pr2-beta-static.zae") # todo: use up-to-date urdf
        self.robot = self.env.GetRobots()[0]
        self.robot.SetActiveManipulator('{0}arm'.format(arm_name))
        self.manip = self.robot.GetActiveManipulator()
        ikmodel = rave.databases.inversekinematics.InverseKinematicsModel(self.robot,iktype=rave.IkParameterization.Type.Transform6D)
        if not ikmodel.load():
            ikmodel.autogenerate()
        
        rospy.sleep(1)
    
    #######################
    # trajectory commands #
    #######################
    
    def execute_pose_trajectory(self, pose_traj, block=True, speed=None):
        """
        :param pose_traj: list of tfx.pose
        :param block: if True, waits until trajectory is completed
        :param speed: execution speed in meters/sec
        """
        joint_traj = [self.ik(pose) for pose in pose_traj]
        print joint_traj
        self.execute_joint_trajectory(joint_traj, block=block, speed=speed)
    
    def execute_joint_trajectory(self, joint_traj, block=True, speed=None):
        """
        :param joint_traj: list of joint waypoints
        :param block: if True, waits until trajectory is completed
        :param speed: execution speed in meters/sec
        """
        for joints in joint_traj:
            assert joints is not None
            assert len(joints) == len(self.current_joints)
        speed = float(speed) if speed is not None else self.default_speed
        
        joint_trajectory = tm.JointTrajectory()
        joint_trajectory.joint_names = self.joint_names
        joint_trajectory.header.stamp = rospy.Time.now()
        
        curr_joints = self.current_joints
        time_from_start = 0.
        for next_joints in joint_traj:
            next_position = self.fk(next_joints).position
            curr_position = self.fk(curr_joints).position
            duration = ((next_position - curr_position).norm)/speed
            time_from_start += duration
            
            waypoint = tm.JointTrajectoryPoint()
            waypoint.positions = next_joints
            waypoint.velocities = [0 for _ in xrange(len(self.joint_names))]
            waypoint.time_from_start = rospy.Duration(time_from_start)
            
            joint_trajectory.points.append(waypoint)
            
            curr_joints = next_joints
            
        self.joint_command_pub.publish(joint_trajectory)
        
        if block:
            rospy.sleep(time_from_start)
        
            
    
    #####################
    # command methods   #
    #####################
    
    def go_to_pose(self, pose, block=True, speed=None):
        """
        :param pose: tfx.pose
        :param block: if True, waits until trajectory is completed
        :param speed: execution speed in meters/sec
        """
        self.execute_pose_trajectory([pose], block=block, speed=speed)
    
    def go_to_joints(self, joints, block=True, speed=None):
        """
        :param joints: list of joint values
        :param block: if True, waits until trajectory is completed
        :param speed: execution speed in meters/sec
        """
        self.execute_joint_trajectory([joints], block=block, speed=speed)
    
    def close_gripper(self, max_effort=None, block=True):
        """
        :param max_effort: max force in Newtons
        :param block: if True, waits until grasp completed
        """
        self.set_gripper(self.min_grasp, max_effort=max_effort, block=block)
    
    def open_gripper(self, max_effort=None, block=True):
        """
        :param max_effort: max force in Newtons
        :param block: if True, waits until grasp completed
        """
        self.set_gripper(self.max_grasp, max_effort=max_effort)
    
    def set_gripper(self, grasp, max_effort=None, block=True):
        """
        :param grasp: in meters
        :param max_effort: max force in Newtons
        :param block: if True, waits until grasp completed
        """
        assert self.min_grasp <= grasp <= self.max_grasp
        max_effort = max_effort if max_effort is not None else self.default_max_effort
        
        self.gripper_command_client.send_goal(pcm.Pr2GripperCommandGoal(pcm.Pr2GripperCommand(position=grasp,max_effort=max_effort)))
        
        if block:
            timeout = utils.Timeout(10)
            timeout.start()
            while(np.abs(self.current_grasp - max(grasp,0)) > .005 and not timeout.has_timed_out()):
                rospy.sleep(.005)
        
    
    #######################
    # state info methods  #
    #######################
    
    def get_pose(self):
        """
        :return current pose of tool_frame as tfx.pose
        """
        return self.fk(self.current_joints)
    
    def get_joints(self):
        """
        :return current joint values
        """
        return self.current_joints
    
    def fk(self, joints):
        """
        :return tfx.pose of tool_frame pose
        """
        self._update_rave()
        self.robot.SetDOFValues(joints, self.manip.GetArmIndices())
        pose_mat = rave_utils.transform_from_to(self.robot, np.eye(4), self.tool_frame, 'base_link')
        return tfx.pose(pose_mat, frame='base_link')
    
    def ik(self, pose):
        """
        :param pose: tfx.pose
        :return list of joints
        """
        assert pose.frame.count('base_link') == 1
        
        self._update_rave()
        goal_pose_mat = rave_utils.transform_relative_pose_for_ik(self.manip, pose.matrix, 'base_link', self.tool_frame)
        joints = self.manip.FindIKSolution(goal_pose_mat, 0) # rave.IkFilterOptions.CheckEnvCollisions
        
        return joints
        
    def _update_rave(self):
        msg = self.joint_state_msg
        
        indices, joint_values = list(), list()
        for name, joint_value in zip(msg.name, msg.position):
            for joint in self.robot.GetJoints():
                if joint.GetName() == name:
                    indices.append(joint.GetDOFIndex())
                    joint_values.append(joint_value)
                    break
                
        self.robot.SetDOFValues(joint_values, indices)
    
    #############
    # callbacks #
    #############
    
    def _joint_state_callback(self, msg):
        """
        Updates self.current_joints
        
        :param msg: sensor_msgs/JointState
        """
        self.joint_state_msg = msg
        
        j = list()
        msg_names_positions = zip(msg.name, msg.position)
        for joint_name in self.joint_names:
            for name, value in msg_names_positions:
                if name == joint_name:
                    j.append(value)
                    break
                
        self.current_joints = j
        
        for name, value in msg_names_positions:
            if name == self.gripper_joint_name:
                self.current_grasp = value
                break
        
        
def test_fk():
    arm = Arm('left')
    actual_pose = tfx.lookupTransform('base_link', 'l_gripper_tool_frame').as_pose()
    fk_pose = arm.fk(arm.current_joints)
    
    print('actual_pose: {0}'.format(actual_pose.matrix))
    print('fk_pose: {0}'.format(fk_pose.matrix))
    
    print('difference: {0}'.format(np.linalg.norm(actual_pose.matrix - fk_pose.matrix)))
    
    IPython.embed()
    
def test_ik():
    arm = Arm('left')
    curr_pose = arm.get_pose()
    
    curr_joints = arm.ik(curr_pose)
    fk_curr_joints = arm.fk(curr_joints)
    
    print('curr_pose: {0}'.format(curr_pose.matrix))
    print('fk(ik(curr_pose)): {0}'.format(fk_curr_joints.matrix))
    print('position difference: {0}'.format((curr_pose.position - fk_curr_joints.position).norm))
    
def test_commands():
    arm = Arm('left')
    curr_pose = arm.get_pose()
    
    #next_pose = curr_pose + [0,0,.2]
    #arm.go_to_pose(next_pose)
    
    pose_traj = [curr_pose + [0,0,-.2],
                 curr_pose + [-.1,.2,-.1],
                 curr_pose + [0,0,0]]
    
    print('pose_traj:')
    for pose in pose_traj:
        print(pose)
    
    arm.execute_pose_trajectory(pose_traj)
    
    IPython.embed()
    
def test_gripper():
    arm = Arm('left')
    
    #arm.open_gripper()
    arm.close_gripper()
    
    IPython.embed()
        
if __name__ == '__main__':
    rospy.init_node('test_arm', anonymous=True)
    #test_fk()
    #test_ik()
    #test_commands()
    test_gripper()