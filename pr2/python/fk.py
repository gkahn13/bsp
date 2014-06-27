import numpy as np
import time

import openravepy as rave
import roslib
roslib.load_manifest('tfx')
import tfx

import tf.transformations as tft

import utils

import IPython

arm_joint_names = ['r_shoulder_pan_joint',   # axis: [0,0,1]
                   'r_shoulder_lift_joint',  # axis: [0,1,0]
                   'r_upper_arm_roll_joint', # axis: [1,0,0]
                   'r_elbow_flex_joint',     # axis: [0,1,0]
                   'r_forearm_roll_joint',   # axis: [1,0,0]
                   'r_wrist_flex_joint',     # axis: [0,1,0]
                   'r_wrist_roll_joint']     # axis: [1,0,0]

arm_joint_axes = [[0,0,1],
                  [0,1,0],
                  [1,0,0],
                  [0,1,0],
                  [1,0,0],
                  [0,1,0],
                  [1,0,0]]

arm_link_trans = [[0., 0.188, 0.],
                 [ 0.1,  0.,  0.],
                 [ 0.,  0.,  0.],
                 [ 0.4,  0.,  0.],
                 [ 0.,  0.,  0.],
                 [ 0.321,  0.,  0.],
                 [0.18, 0, 0]]


arm_link_names = ['torso_lift_link',
                  'r_shoulder_pan_link',
                  'r_shoulder_lift_link',
                  'r_upper_arm_link',
                  'r_elbow_flex_link',
                  'r_forearm_link',
                  'r_wrist_flex_link',
                  'r_gripper_tool_frame']

class RobotFK:
    def __init__(self, robot, arm_name):
        """
        self.env = robot.GetEnv()
        self.robot = robot
        self.handles = list()
        
        self.arm = robot.GetManipulator(arm_name)
        
        self.arm_joints = robot.GetJoints(self.arm.GetArmIndices())
        self.num_joints = len(self.arm_joints)
        
        self.arm_links = [j.GetHierarchyParentLink() for j in self.arm_joints]
        arm_link_positions = [tfx.pose(l.GetTransform()).position for l in self.arm_links]
        self.link_trans = [(arm_link_positions[i+1] - arm_link_positions[i]).array for i in xrange(len(self.arm_links)-1)]
        self.link_trans.append([.18,0,0])
        #self.arm_link_lengths = [(arm_link_positions[i+1] - arm_link_positions[i]).norm for i in xrange(len(arm_links)-1)]
        
        """
        self.origin = robot.GetLink('torso_lift_link').GetTransform()
        
        self.arm_joint_axes =  [[0,0,1],
                                  [0,1,0],
                                  [1,0,0],
                                  [0,1,0],
                                  [1,0,0],
                                  [0,1,0],
                                  [1,0,0]]
        
        self.arm_link_trans = [[0., 0.188, 0.],
                             [ 0.1,  0.,  0.],
                             [ 0.,  0.,  0.],
                             [ 0.4,  0.,  0.],
                             [ 0.,  0.,  0.],
                             [ 0.321,  0.,  0.],
                             [0.18, 0, 0]]
        
        
        

    def fk(self, joint_values):
        pose_mat = self.origin
        
        R = np.eye(4)
        for i, j in enumerate(joint_values[:]):
            rot = tft.rotation_matrix(j, self.arm_joint_axes[i])
            trans = self.arm_link_trans[i]
            R[:3,:3] = rot[:3,:3]
            R[:3,3] = trans
            pose_mat = np.dot(pose_mat, R)
            
        return tfx.pose(pose_mat)

def test_fk():
    env = rave.Environment()
    env.Load('../robots/pr2-beta-sim.robot.xml')
    env.SetViewer('qtcoin')
    time.sleep(1)
    
    robot = env.GetRobots()[0]
    arm = robot.GetManipulator('leftarm')
    
    arm_joints = robot.GetJoints(arm.GetArmIndices())
    arm_links = [j.GetHierarchyParentLink() for j in arm_joints]
    arm_link_positions = [tfx.pose(l.GetTransform()).position for l in arm_links]
    arm_link_lengths = [(arm_link_positions[i+1] - arm_link_positions[i]).norm for i in xrange(len(arm_links)-1)]
    
    transforms = robot.GetLinkTransformations(0)
    indices = arm.GetArmIndices()
    arm_transforms = transforms[indices[0]:indices[-1]]
    
    robot_fk = RobotFK(robot, 'leftarm')
    
    joints = [np.random.random(7)-.5, np.random.random(7)-.5]
    
    for j in joints:
        robot.SetDOFValues(j, arm.GetArmIndices())
        
        first_link_pose = tfx.pose(arm_links[0].GetTransform())#[:3,3])
        last_link_fk = robot_fk.fk(arm.GetArmDOFValues())
        #last_link_actual = tfx.pose(arm_links[-1].GetTransform())
        last_link_actual = tfx.pose(utils.openraveTransformFromTo(robot,np.eye(4),'l_gripper_tool_frame','world'))
        
        print('first_link_pose: {0}'.format(first_link_pose.matrix))
        
        trans_err = (last_link_fk.position-last_link_actual.position).norm
        rot_err = np.linalg.norm(last_link_fk.matrix[:3,:3] - last_link_actual.matrix[:3,:3])
        print('trans_err: {0}'.format(trans_err))
        print('rot_err: {0}'.format(rot_err))
    
    IPython.embed()

if __name__ == '__main__':
    test_fk()