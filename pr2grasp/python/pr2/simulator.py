import rospy
import sensor_msgs.msg as sm

import openravepy as rave
import numpy as np

class Simulator:
    """ OpenRave simulator """
    def __init__(self, env_file='robots/pr2-beta-static.zae', view=False):
        self.joint_state_msg = None
        self.joint_state_sub = rospy.Subscriber('/joint_states', sm.JointState, self._joint_state_callback)
        
        self.env = rave.Environment()
        self.env.StopSimulation()
        self.env.Load(env_file)
        
        self.handles = list()
        self.view = view
        if view:
            self.env.SetViewer('qtcoin')
        
        self.robot = self.env.GetRobots()[0]
        self.larm = self.robot.GetManipulator('leftarm')
        self.rarm = self.robot.GetManipulator('rightarm')
        
        for arm in [self.larm, self.rarm]:
            self.robot.SetActiveManipulator(arm)
        
            ikmodel = rave.databases.inversekinematics.InverseKinematicsModel(self.robot,iktype=rave.IkParameterization.Type.Transform6D)
            if not ikmodel.load():
                ikmodel.autogenerate()
        
    def update(self):
        """
        Updates robot joints to match ROS /joint_states
        """
        msg = self.joint_state_msg
        if msg is None:
            return
        
        indices, joint_values = list(), list()
        for name, joint_value in zip(msg.name, msg.position):
            for joint in self.robot.GetJoints():
                if joint.GetName() == name:
                    indices.append(joint.GetDOFIndex())
                    joint_values.append(joint_value)
                    break
                
        self.robot.SetDOFValues(joint_values, indices)
        
    def _joint_state_callback(self, msg):
        """
        :type msg: sensor_msgs/JointState
        """
        self.joint_state_msg = msg
        
    #################
    # frame methods #
    #################
        
    def transform_from_to(self, pose_mat_in_ref, ref_link_name, targ_link_name):
        """
        :param pose_mat_in_ref: 4x4 np.array
        :param ref_link_name: string
        :param targ_link_name: string
        """
        pose_mat_in_ref = np.array(pose_mat_in_ref)
        
        # ref -> world
        if ref_link_name != 'world':
            ref_from_world = self.robot.GetLink(ref_link_name).GetTransform()
        else:
            ref_from_world = np.eye(4)
    
        # target -> world
        if targ_link_name != 'world':
            targ_from_world = self.robot.GetLink(targ_link_name).GetTransform()
        else:
            targ_from_world = np.eye(4)
    
        # target -> ref
        targ_from_ref = np.dot(np.linalg.inv(targ_from_world), ref_from_world)
    
        pose_mat_in_targ = np.dot(targ_from_ref, pose_mat_in_ref)
        return np.array(pose_mat_in_targ)
    
    def transform_relative_pose_for_ik(self, manip, matrix4, ref_frame, targ_frame):
        """
        Transforms the matrix to be used for IK
        (needed since last link in manipulator is not tool_frame
        
        :param manip: OpenRave manipulator
        :param matrix4: 4x4 np.array
        :param ref_frame: string reference frame
        :param targ_frame: string target frame
        """
        if ref_frame == 'world':
            world_from_ref = np.eye(4)
        else:
            ref = self.robot.GetLink(ref_frame)
            world_from_ref = ref.GetTransform()
    
        if targ_frame == 'end_effector':        
            targ_from_EE = np.eye(4)
        else:
            world_from_targ = self.robot.GetLink(targ_frame).GetTransform()
            world_from_EE = manip.GetEndEffectorTransform()    
            targ_from_EE = np.dot(np.linalg.inv(world_from_targ), world_from_EE)       
    
        ref_from_targ_new = matrix4
        world_from_EE_new = np.dot(np.dot(world_from_ref, ref_from_targ_new), targ_from_EE)    
    
        return np.array(world_from_EE_new)
    
    ############
    # plotting #
    ############
    
    def plot_point(self, pos_array, size=.01, color=(0,1,0)):
        """
        :param pos_array: 3d np.array
        :param size: radius in meters
        :param color: rgb [0,1], default green
        """
        self.handles += [self.env.plot3(points=pos_array,
                                        pointsize=size,
                                        colors=np.array(color),
                                        drawstyle=1)]
    
    def plot_segment(self, start, end, color=(1,0,0)):
        """
        :param start: 3d np.array
        :param end: 3d np.array
        :param color: rgb [0,1], default red
        """
        start = np.array(start)
        end = np.array(end)
        
        self.handles += [self.env.drawlinestrip(points=np.array([start, end]), linewidth=3.0, colors=np.array([color,color]))]
    
    def plot_transform(self, T, s=0.1):
        """
        :param T: 4x4 np.array
        :param s: length of axes in meters
        """
        T = np.array(T)
        h = []
        x = T[0:3,0]
        y = T[0:3,1]
        z = T[0:3,2]
        o = T[0:3,3]
        self.handles.append(self.env.drawlinestrip(points=np.array([o, o+s*x]), linewidth=3.0, colors=np.array([(1,0,0),(1,0,0)])))
        self.handles.append(self.env.drawlinestrip(points=np.array([o, o+s*y]), linewidth=3.0, colors=np.array(((0,1,0),(0,1,0)))))
        self.handles.append(self.env.drawlinestrip(points=np.array([o, o+s*z]), linewidth=3.0, colors=np.array(((0,0,1),(0,0,1)))))
    
    def save_view(self, file_name):
        """
        :param file_name: path string
        """
        self.env.GetViewer().SendCommand('SetFiguresInCamera 1') # also shows the figures in the image
        I = self.env.GetViewer().GetCameraImage(640,480,  self.env.GetViewer().GetCameraTransform(),[640,640,320,240])
        scipy.misc.imsave(file_name ,I)
        self.env.GetViewer().SendCommand('SetFiguresInCamera 0')
       
    