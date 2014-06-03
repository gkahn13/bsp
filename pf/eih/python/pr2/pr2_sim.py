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
    def __init__(self, env_file, view=True):
        self.env = rave.Environment()
        self.env.Load(env_file)
        if view:
            self.env.SetViewer('qtcoin')
        self.robot = self.env.GetRobots()[0]
        
        self.larm = Arm(self, "l")
        self.rarm = Arm(self, "r")
        self.head = Head(self)
        
        self.h_kinect = KinectSensor(self.robot, 'head_depth', 'head_cam')
        self.r_kinect = KinectSensor(self.robot, 'r_gripper_depth', 'r_gripper_cam')
        self.l_kinect = KinectSensor(self.robot, 'l_gripper_depth', 'l_gripper_cam')
        
class Arm:
    L_POSTURES = dict(        
        untucked = [0.4,  1.0,   0.0,  -2.05,  0.0,  -0.1,  0.0],
        tucked = [0.06, 1.25, 1.79, -1.68, -1.73, -0.10, -0.09],
        up = [ 0.33, -0.35,  2.59, -0.15,  0.59, -1.41, -0.27],
        side = [  1.832,  -0.332,   1.011,  -1.437,   1.1  ,  -2.106,  3.074],
        mantis = [ 2.03018192, -0.05474993, 1.011, -1.47618716,  0.55995636, -1.42855926,  3.96467305]
    ) 
        
    def __init__(self, pr2, lr):
        self.lr = lr
        self.lrlong = {"r":"right", "l":"left"}[lr]
        self.tool_frame = "%s_gripper_tool_frame"%lr

        self.pr2 = pr2
        self.robot = pr2.robot
        self.manip = self.robot.GetManipulator("%sarm"%self.lrlong)
        self.arm_indices = self.manip.GetArmIndices()
        
    def get_joint_values(self):
        return self.manip.GetArmDOFValues()
    
    def get_pose(self):
        """ Returns pose of end effector """
        return tfx.pose(self.manip.GetEndEffectorTransform(), frame='world')
    
    def get_limits(self):
        l_limits, u_limits = self.robot.GetDOFLimits()
        arm_l_limits = np.array([l_limits[i] for i in self.arm_indices])
        arm_u_limits = np.array([u_limits[i] for i in self.arm_indices])
        return arm_l_limits, arm_u_limits
        
    def set_posture(self, name):
        l_joints = self.L_POSTURES[name]        
        joints = l_joints if self.lr == 'l' else mirror_arm_joints(l_joints)
        self.set_joint_values(joints)

    def set_joint_values(self, joint_values):
        l_limits, u_limits = self.get_limits()
        valid_joint_values = np.array(joint_values)
        valid_joint_values = np.maximum(joint_values, l_limits)
        valid_joint_values = np.minimum(joint_values, u_limits)
        if not np.min(np.array(joint_values) == valid_joint_values):
            print('Invalid joints. Setting invalid values to limits')
        self.robot.SetDOFValues(joint_values, self.arm_indices)
        
    def set_pose(self, pose):
        ref_frame = pose.frame if pose.frame is not None else "world"
        joint_values = cart_to_joint(self.manip, np.array(pose.matrix), ref_frame, "end_effector")
        if joint_values is not None:
            self.set_joint_values(joint_values)
        else:
            print('pose_matrix invalid, ignoring command')
            
    def teleop(self):
        print('{0} arm teleop'.format(self.lrlong))
        
        pos_step = .01
        delta_position = {'a' : [0, pos_step, 0],
                          'd' : [0, -pos_step, 0],
                          'w' : [pos_step, 0, 0],
                          'x' : [-pos_step, 0, 0],
                          '+' : [0, 0, pos_step],
                          '-' : [0, 0, -pos_step]}
        
        angle_step = 2.0
        delta_angle = {'o' : [angle_step, 0, 0],
                       'p' : [-angle_step, 0, 0],
                       'k' : [0, angle_step, 0],
                       'l' : [0, -angle_step, 0],
                       'n' : [0, 0, angle_step],
                       'm' : [0, 0, -angle_step]}
        
        char = ''
        while char != 'q':
            char = utils.Getch.getch()
            pose = self.get_pose()
            new_pose = tfx.pose(pose)
            if delta_position.has_key(char):
                new_pose.position = pose.position.array + delta_position[char]
            ypr = np.array([pose.tb_angles.yaw_deg, pose.tb_angles.pitch_deg, pose.tb_angles.roll_deg])
            if delta_angle.has_key(char):
                ypr += np.array(delta_angle[char])
                new_pose = tfx.pose(pose.position, tfx.tb_angles(ypr[0], ypr[1], ypr[2]))
            self.set_pose(new_pose)    
            time.sleep(.01)
            
        print('{0} arm end teleop'.format(self.lrlong))
            
class Head:
    joint_names = ['head_pan_joint', 'head_tilt_joint']
    def __init__(self, pr2):
        self.robot = pr2.robot
        self.head_indices = [self.robot.GetJointIndex(joint_name) for joint_name in Head.joint_names]
    
    def get_joint_values(self):
        return self.robot.GetDOFValues(self.head_indices)
        
    def get_limits(self):
        l_limits, u_limits = self.robot.GetDOFLimits()
        head_l_limits = np.array([l_limits[i] for i in self.head_indices])
        head_u_limits = np.array([u_limits[i] for i in self.head_indices])
        return head_l_limits, head_u_limits
        
    def set_joint_values(self, joint_values):
        self.robot.SetDOFValues(joint_values, self.joint_indices)
        
    def look_at(self, pose, camera_frame='wide_stereo_link'):
        reference_frame = pose.frame if pose.frame is not None else "world" 
        
        worldFromRef = self.robot.GetLink(reference_frame).GetTransform() if reference_frame != 'world' else np.eye(4)
        worldFromCam = self.robot.GetLink(camera_frame).GetTransform()
        refFromCam = dot(inv(worldFromRef), worldFromCam)        
    
        xyz_cam = refFromCam[:3,3]
        ax = pose.position - xyz_cam # pointing axis
        pan = np.arctan(ax[1]/ax[0])
        tilt = np.arcsin(-ax[2] / norm(ax))
        self.set_joint_values([pan,tilt])
        
class Sensor:
    def __init__(self, sensor):
        self.sensor = sensor
        self.is_powered = False
        self.is_rendering = False
        
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
        
        
        self.power_on()
        time.sleep(1)
        data = self.get_data()
        
        self.P = data.KK
        self.height, self.width = data.imagedata.shape[0:2]
        
        self.power_off()
        
    def get_image(self):
        data = self.get_data()
        return data.imagedata
    
    def get_pixels_and_colors(self, points):
        """
        3d points --> 2d pixel coordinates
        Only for points in FOV
        
        Returns
        [(point, pixel, color), ...]
        """
        data = self.get_data()
        
        points_pixels_colors = list()
        for x in points:

            if self.is_in_fov(x):
                pixel = self.get_pixel_from_point(x)
                color = data.imagedata[pixel[0], pixel[1]]
                points_pixels_colors.append((x, np.array(pixel), np.array(color, dtype=float)))
            
        return points_pixels_colors
    
    def get_pixel_from_point(self, x):
        x_mat = x.as_pose().matrix
        x_mat[0:3,0:3] = np.zeros((3,3))
        xtilde = np.dot(inv(self.sensor.GetTransform()), x_mat)
        
        #x1, x2, x3 = xtilde[0,3], xtilde[1,3], xtilde[2,3]
        #f = self.P[0,0]
        #y1 = (f/x3)*x1
        #y2 = (f/x3)*x2
        
        y = np.dot(self.P, xtilde[0:3,3])
        y_pixel = int(y[1]/y[2])
        x_pixel = int(y[0]/y[2])
        
        return (y_pixel, x_pixel)
    
    def is_in_fov(self, x):
        """ True if x (tfx.point) is in the field-of-view """
        y_pixel, x_pixel = self.get_pixel_from_point(x)
        return ((0 <= y_pixel < self.height) and (0 <= x_pixel < self.width))
    
class KinectSensor:
    """
    Simulates a kinect by having a depth sensor + camera sensor
    """
    def __init__(self, robot, depth_sensor_name, camera_sensor_name):
        self.robot = robot
        self.depth_sensor = DepthSensor(robot.GetSensor(depth_sensor_name).GetSensor())
        self.camera_sensor = CameraSensor(robot.GetSensor(camera_sensor_name).GetSensor())
                
        geom = self.camera_sensor.sensor.GetSensorGeometry(rave.Sensor.Type.Camera)
        self.height, self.width = geom.imagedata.shape[0:2]
                
    def power_on(self):
        self.depth_sensor.power_on()
        self.camera_sensor.power_on()
        
    def power_off(self):
        self.depth_sensor.power_off()
        self.camera_sensor.power_off()
        
    def get_point_cloud(self):
        """ Returns a PointCloud, i.e. a list of ColoredPoint """
        points = self.depth_sensor.get_points()
        points_pixels_colors = self.camera_sensor.get_pixels_and_colors(points)
        
        colored_points = list()
        for point, pixel, color in points_pixels_colors:
            colored_points.append(ColoredPoint(point, color))
            
        return colored_points
    
    def get_z_buffer(self):
        points = self.depth_sensor.get_points()
        points_pixels_colors = self.camera_sensor.get_pixels_and_colors(points)
        
        origin_pos = self.get_pose().position
        
        from collections import defaultdict
        z_buffer = defaultdict(lambda: np.inf)
        for point, pixel, color in points_pixels_colors:
            z_buffer[tuple(pixel)] = (point - origin_pos).norm
            
        return z_buffer
    
    def get_image(self):
        return self.camera_sensor.get_image()
    
    def get_pixel_from_point(self, x):
        return self.camera_sensor.get_pixel_from_point(x)
    
    def is_in_fov(self, x):
        return self.camera_sensor.is_in_fov(x)
    
    def get_pose(self):
        return tfx.pose(self.camera_sensor.sensor.GetTransform())
    
    def distance_to(self, x):
        """ x -- tfx.point """
        origin_pos = tfx.pose(self.camera_sensor.sensor.GetTransform()).position
        return (x - origin_pos).norm
                
    def render_on(self):
        self.depth_sensor.render_on()
        self.camera_sensor.render_on()
        
    def render_off(self):
        self.depth_sensor.render_off()
        self.camera_sensor.render_off()
        
    def display_point_cloud(self, colored_points):
        env = self.robot.GetEnv()
        self.handles = list()
        for colored_point in colored_points:
            self.handles.append(colored_point.display(env))
        
        
class ColoredPoint:
    """
    Wrapper to simulate components of a point cloud
    Contains a point and a color
    """
    def __init__(self, point, color):
        """
        point -- tfx point
        color -- 3d array with values 0..255
        """
        self.point = point
        self.color = color
        
    def display(self, env):
        return utils.plot_point(env, self.point.array, color=self.color/255.)
    
def test():
    brett = PR2('../../envs/pr2-test.env.xml')
    env = brett.env
    
    larm = brett.larm
    rarm = brett.rarm
    head = brett.head
    h_kinect = brett.h_kinect
    l_kinect = brett.l_kinect
    r_kinect = brett.r_kinect
    
    kinect = r_kinect
    
    larm.set_posture('mantis')
    rarm.set_posture('mantis')
    
    kinect.power_on()
    kinect.render_on()
    time.sleep(1)
    colored_points = kinect.get_point_cloud()
    z_buffer = kinect.get_z_buffer()

    kinect.display_point_cloud(colored_points)
        
    image = kinect.camera_sensor.get_data().imagedata
    points_pixels_colors = kinect.camera_sensor.get_pixels_and_colors(kinect.depth_sensor.get_points())
    for point, pixel, color in points_pixels_colors:
        image[pixel[0],pixel[1]] = [0, 255, 0]
    
    plt.imshow(image)
    plt.show(block=False)
    
    
    
    
    IPython.embed()
    
def test_pose():
    brett = PR2('../../envs/pr2-test.env.xml')
    
    larm = brett.larm
    
    larm.set_posture('side')
    start_pose = larm.get_pose()
    print('start joints: {0}'.format(larm.get_joint_values()))
    print('start pose: {0}'.format(start_pose.matrix))
    
    larm.set_posture('mantis')
    print('mantis joints: {0}'.format(larm.get_joint_values()))
    print('Press enter to go back to start_pose')
    raw_input()
    
    larm.set_pose(start_pose)
    print('end joints: {0}'.format(larm.get_joint_values()))
    
    IPython.embed()
    
if __name__ == '__main__':
    test()
    #test_pose()