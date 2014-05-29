import IPython
import time

import numpy as np
from numpy import inf, zeros, dot, r_
from numpy.linalg import norm, inv
from matplotlib import pyplot as plt

import kinematics_utils as ku
import utils

import openravepy as rave

import roslib
roslib.load_manifest('tfx')
import tfx

def mirror_arm_joints(x):
    "mirror image of joints (r->l or l->r)"
    return r_[-x[0],x[1],-x[2],x[3],-x[4],x[5],-x[6]]

def transform_relative_pose_for_ik(manip, matrix4, ref_frame, targ_frame):
    robot = manip.GetRobot()


    if ref_frame == "world":
        worldFromRef = np.eye(4)
    else:
        ref = robot.GetLink(ref_frame)
        worldFromRef = ref.GetTransform()

    if targ_frame == "end_effector":        
        targFromEE = np.eye(4)
    else:
        targ = robot.GetLink(targ_frame)
        worldFromTarg = targ.GetTransform()
        worldFromEE = manip.GetEndEffectorTransform()    
        targFromEE = dot(inv(worldFromTarg), worldFromEE)       

    refFromTarg_new = matrix4
    worldFromEE_new = dot(dot(worldFromRef, refFromTarg_new), targFromEE)    

    return worldFromEE_new

def cart_to_joint(manip, matrix4, ref_frame, targ_frame, filter_options = 0):
    robot = manip.GetRobot()
    worldFromEE = transform_relative_pose_for_ik(manip, matrix4, ref_frame, targ_frame)        
    joint_positions = manip.FindIKSolution(worldFromEE, filter_options)
    if joint_positions is None: return joint_positions
    current_joints = robot.GetDOFValues(manip.GetArmIndices())
    joint_positions_unrolled = ku.closer_joint_angles(joint_positions, current_joints)
    return joint_positions_unrolled

class PR2:
    sensor_names = ['flashlidar3d','cam']
    def __init__(self, env_file='../robots/pr2-beta-static.zae'):
        self.env = rave.Environment()
        #self.env.StopSimulation()
        self.env.Load(env_file)
        self.robot = self.env.GetRobots()[0]
        
        self.larm = Arm(self, "l")
        self.rarm = Arm(self, "r")
        self.head = Head(self)
        
        self.depth_sensor = DepthSensor(self.robot.GetSensor('flashlidar3d').GetSensor())
        self.cam_sensor = CameraSensor(self.robot.GetSensor('cam').GetSensor())
        
class Arm:
    L_POSTURES = dict(        
        untucked = [0.4,  1.0,   0.0,  -2.05,  0.0,  -0.1,  0.0],
        tucked = [0.06, 1.25, 1.79, -1.68, -1.73, -0.10, -0.09],
        up = [ 0.33, -0.35,  2.59, -0.15,  0.59, -1.41, -0.27],
        side = [  1.832,  -0.332,   1.011,  -1.437,   1.1  ,  -2.106,  3.074]
    ) 
        
    def __init__(self, pr2, lr):
        self.lr = lr
        self.lrlong = {"r":"right", "l":"left"}[lr]
        self.tool_frame = "%s_gripper_tool_frame"%lr

        self.robot = pr2.robot
        self.manip = self.robot.GetManipulator("%sarm"%self.lrlong)
        self.arm_indices = self.manip.GetArmIndices()
        
    def get_joint_values(self):
        return self.manip.GetArmDOFValues()
    
    def get_pose(self):
        """ Returns pose of end effector """
        return tfx.pose(self.manip.GetEndEffectorTransform(), frame="world")
        
    def set_posture(self, name):
        l_joints = self.L_POSTURES[name]        
        joints = l_joints if self.lr == 'l' else mirror_arm_joints(l_joints)
        self.set_joint_values(joints)

    def set_joint_values(self, joint_values):
        self.robot.SetDOFValues(joint_values, self.arm_indices)
        
    def set_pose(self, pose):
        ref_frame = pose.frame if pose.frame is not None else "world"
        joint_values = cart_to_joint(self.manip, np.array(pose.matrix), ref_frame, "end_effector")
        if joint_values is not None:
            self.set_joint_values(joint_values)
        else:
            print('pose_matrix invalid, ignoring command')

class Head:
    joint_names = ['head_pan_joint', 'head_tilt_joint']
    def __init__(self, pr2):
        self.robot = pr2.robot
        self.joint_indices = [self.robot.GetJointIndex(joint_name) for joint_name in Head.joint_names]
        
    def set_pan_tilt(self, pan, tilt):
        self.robot.SetDOFValues([pan, tilt], self.joint_indices)
        
    def look_at(self, pose, camera_frame='wide_stereo_link'):
        reference_frame = pose.frame if pose.frame is not None else "world" 
        
        worldFromRef = self.robot.GetLink(reference_frame).GetTransform() if reference_frame != 'world' else np.eye(4)
        worldFromCam = self.robot.GetLink(camera_frame).GetTransform()
        refFromCam = dot(inv(worldFromRef), worldFromCam)        
    
        xyz_cam = refFromCam[:3,3]
        ax = pose.position - xyz_cam # pointing axis
        pan = np.arctan(ax[1]/ax[0])
        tilt = np.arcsin(-ax[2] / norm(ax))
        self.set_pan_tilt(pan,tilt)
        
class Sensor:
    def __init__(self, sensor):
        self.sensor = sensor
        self.is_powered = False
        self.is_rendering = False
        
        print(self.sensor.GetName())
        for type in rave.Sensor.Type.values.values():
            if self.sensor.Supports(type):
                self.type = type
                break
        else:
            # invalid sensor!
            self.type = -1
        
    def power_on(self):
        if not self.is_powered:
            self.sensor.Configure(rave.Sensor.ConfigureCommand.PowerOn)
            self.is_powered = True
        
    def power_off(self):
        if self.is_powered:
            self.sensor.Configure(rave.Sensor.ConfigureCommand.PowerOff)
            self.is_powered = False
        
    def get_data(self):
        return self.sensor.GetSensorData(self.type)
            
    def render_on(self):
        if not self.is_rendering:
            self.sensor.Configure(rave.Sensor.ConfigureCommand.RenderDataOn)
            self.is_rendering = True
        
    def render_off(self):
        if self.is_rendering:
            self.sensor.Configure(rave.Sensor.ConfigureCommand.RenderDataOff)
            self.is_rendering = False
            
class DepthSensor(Sensor):
    def __init__(self, sensor):
        Sensor.__init__(self, sensor)
        
    def get_depths(self):
        data = self.get_data()
        points = data.ranges
        depths = np.sqrt((points**2).sum(axis=1))
        return depths
        
    def get_points(self):
        data = self.get_data()
        point_mat = data.ranges
        in_range = np.array(data.intensity, dtype=int)
        
        sensor_pose = tfx.pose(self.sensor.GetTransform())
        points = list()
        for i in xrange(point_mat.shape[0]):
            if in_range[i] == 1:
                points.append(tfx.point(sensor_pose.position + point_mat[i,:]))
                
        return points
    
class CameraSensor(Sensor):
    def __init__(self, sensor):
        Sensor.__init__(self, sensor)
        
    def get_image(self):
        data = self.get_data()
        return data.imagedata
    
    def get_pixels_and_colors(self, points):
        """
        3d points --> 2d pixel coordinates
        Assumes all points within FOV
        
        Returns
        [(pixel, color), ...]
        """
        data = self.get_data()
        P = data.KK
        
        pixels_colors = list()
        for i in xrange(len(points)):
            x = points[i]
            
            x_mat = x.as_pose().matrix
            x_mat[0:3,0:3] = np.zeros((3,3))
            xtilde = np.dot(inv(self.sensor.GetTransform()), x_mat)
            
            x1, x2, x3 = xtilde[0,3], xtilde[1,3], xtilde[2,3]
            f = P[0,0]
            y1 = (f/x3)*x1
            y2 = (f/x3)*x2
            
            y = np.dot(P, xtilde[0:3,3])
            y_pixel = int(y[1]/y[2])
            x_pixel = int(y[0]/y[2])
            
            color = data.imagedata[y_pixel, x_pixel]
            pixels_colors.append((np.array([y_pixel, x_pixel]), np.array(color, dtype=float)))
            
        return pixels_colors
    
def test():
    brett = PR2('../envs/pr2-test.env.xml')
    env = brett.env
    env.SetViewer('qtcoin')
    
    larm = brett.larm
    rarm = brett.rarm
    head = brett.head
    depth_sensor = brett.depth_sensor
    cam_sensor = brett.cam_sensor
    
    depth_sensor.power_on()
    time.sleep(1)
    points = depth_sensor.get_points()
    depth_sensor.power_off()
    
    #handles = []
    #for p in points:
    #    handles.append(utils.plot_point(env, p.array))
        
    cam_sensor.power_on()
    time.sleep(1)
    cam_sensor.render_on()
    data = cam_sensor.get_data()
    
    
    pixels_colors = cam_sensor.get_pixels_and_colors(points)
    handles = []
    for point, pixel_color in zip(points, pixels_colors):
        pixel, color = pixel_color
        handles.append(utils.plot_point(env, point.array, color=color/255.))
    
    IPython.embed()
    
if __name__ == '__main__':
    test()