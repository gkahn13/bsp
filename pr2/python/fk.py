import numpy as np
import time

import openravepy as rave
import roslib
roslib.load_manifest('tfx')
import tfx

import tf.transformations as tft

import utils

import IPython

arm_joint_names = ['_shoulder_pan_joint',   # axis: [0,0,1]
                   '_shoulder_lift_joint',  # axis: [0,1,0]
                   '_upper_arm_roll_joint', # axis: [1,0,0]
                   '_elbow_flex_joint',     # axis: [0,1,0]
                   '_forearm_roll_joint',   # axis: [1,0,0]
                   '_wrist_flex_joint',     # axis: [0,1,0]
                   '_wrist_roll_joint']     # axis: [1,0,0]

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

start_link = 'torso_lift_link'

arm_link_names = ['_shoulder_pan_link',
                  '_shoulder_lift_link',
                  '_upper_arm_link',
                  '_elbow_flex_link',
                  '_forearm_link',
                  '_wrist_flex_link',
                  '_gripper_tool_frame']

class RobotFK:
    def __init__(self, robot, lr, gripper_tool_to_sensor):
        self.env = robot.GetEnv()
        self.robot = robot
        self.handles = list()
        """
        
        self.arm = robot.GetManipulator(arm_name)
        
        self.arm_joints = robot.GetJoints(self.arm.GetArmIndices())
        self.num_joints = len(self.arm_joints)
        
        self.arm_links = [j.GetHierarchyParentLink() for j in self.arm_joints]
        arm_link_positions = [tfx.pose(l.GetTransform()).position for l in self.arm_links]
        self.link_trans = [(arm_link_positions[i+1] - arm_link_positions[i]).array for i in xrange(len(self.arm_links)-1)]
        self.link_trans.append([.18,0,0])
        #self.arm_link_lengths = [(arm_link_positions[i+1] - arm_link_positions[i]).norm for i in xrange(len(arm_links)-1)]
        """
        
        self.joint_names = [lr+name for name in arm_joint_names]
        self.link_names = [start_link] + [lr+name for name in arm_link_names]
        
        links = [robot.GetLink(name) for name in self.link_names]
        link_positions = [tfx.pose(l.GetTransform()).position for l in links]
        self.link_trans = [(link_positions[i+1] - link_positions[i]).array for i in xrange(len(link_positions)-1)]
        
        self.gripper_tool_to_sensor = gripper_tool_to_sensor
        
        self.origin = robot.GetLink('torso_lift_link').GetTransform()
        
        self.arm_joint_axes =  [[0,0,1],
                                  [0,1,0],
                                  [1,0,0],
                                  [0,1,0],
                                  [1,0,0],
                                  [0,1,0],
                                  [1,0,0]]
        
        self.arm_link_trans = [[0., 0.188 if lr == 'l' else -0.188, 0.],
                             [ 0.1,  0.,  0.],
                             [ 0.,  0.,  0.],
                             [ 0.4,  0.,  0.],
                             [ 0.,  0.,  0.],
                             [ 0.321,  0.,  0.],
                             [0.18, 0, 0]]
        
        
        

    def fk(self, joint_values):
        self.handles = list()
        pose_mat = self.origin
        
        R = np.eye(4)
        for i, j in enumerate(joint_values[:]):
            rot = tft.rotation_matrix(j, self.arm_joint_axes[i])
            trans = self.arm_link_trans[i]
            R[:3,:3] = rot[:3,:3]
            R[:3,3] = trans
            print('{0}: {1}\n'.format(i,R))
            pose_mat = np.dot(pose_mat, R)
            
            self.handles += utils.plot_transform(self.env, pose_mat, .3)
          
        pose_mat = np.dot(pose_mat, self.gripper_tool_to_sensor)
        self.handles += utils.plot_transform(self.env, pose_mat, .3)
            
        return tfx.pose(pose_mat)

def test_fk():
    env = rave.Environment()
    env.Load('../robots/pr2-beta-sim.robot.xml')
    env.SetViewer('qtcoin')
    time.sleep(1)
    
    arm_name = 'left'
    lr = arm_name[0]
    
    robot = env.GetRobots()[0]
    arm = robot.GetManipulator(arm_name+'arm')
    sensor = robot.GetAttachedSensor(lr+'_gripper_cam')
    gripper_tool_to_sensor = utils.openraveTransformFromTo(robot, sensor.GetTransform(), 'world', lr+'_gripper_tool_frame')
    
    """
    arm_joints = robot.GetJoints(arm.GetArmIndices())
    arm_links = [j.GetHierarchyParentLink() for j in arm_joints]
    arm_link_positions = [tfx.pose(l.GetTransform()).position for l in arm_links]
    arm_link_lengths = [(arm_link_positions[i+1] - arm_link_positions[i]).norm for i in xrange(len(arm_links)-1)]
    
    transforms = robot.GetLinkTransformations(0)
    indices = arm.GetArmIndices()
    arm_transforms = transforms[indices[0]:indices[-1]]
    """
    
    robot_fk = RobotFK(robot, lr, gripper_tool_to_sensor)
    
    #joints = [np.random.random(7)-.5, np.random.random(7)-.5]
    #joints = [[-2.03018, -0.0547499   ,  -1.011 ,  -1.47619,  -0.559956 ,  -1.42856  , -3.96467]]
    joints = [[-0.0035595997252855227, -0.0001802885031949586, 0.009552644767137686, 0.0011127307030500688, -0.2960540594183074, -0.41989261289657076, 0.0005531658788795468]]

    for j in joints:
        robot.SetDOFValues(j, arm.GetArmIndices())
        
        last_link_fk = robot_fk.fk(arm.GetArmDOFValues())
        #last_link_actual = tfx.pose(utils.openraveTransformFromTo(robot,np.eye(4),arm_name[0]+'_gripper_tool_frame','world'))
        last_link_actual = tfx.pose(sensor.GetTransform())
        
        trans_err = (last_link_fk.position-last_link_actual.position).norm
        rot_err = np.linalg.norm(last_link_fk.matrix[:3,:3] - last_link_actual.matrix[:3,:3])
        print('trans_err: {0}'.format(trans_err))
        print('rot_err: {0}'.format(rot_err))
    
    IPython.embed()

if __name__ == '__main__':
    test_fk()