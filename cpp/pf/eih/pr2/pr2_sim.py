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
    def __init__(self, env_file):
        self.env = rave.Environment()
        self.env.Load(env_file)
        self.robot = self.env.GetRobots()[0]
        
        self.larm = Arm(self, "l")
        self.rarm = Arm(self, "r")
        self.head = Head(self)
        
        # 'head_depth', 'head_cam'
        # 'r_gripper_depth', 'r_gripper_cam'
        depth_sensor = DepthSensor(self.robot.GetSensor('r_gripper_depth').GetSensor())
        cam_sensor = CameraSensor(self.robot.GetSensor('r_gripper_cam').GetSensor())
        self.kinect_sensor = KinectSensor(depth_sensor, cam_sensor)
        
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

        self.pr2 = pr2
        self.robot = pr2.robot
        self.manip = self.robot.GetManipulator("%sarm"%self.lrlong)
        self.arm_indices = self.manip.GetArmIndices()
        
    def get_joint_values(self):
        return self.manip.GetArmDOFValues()
    
    def get_pose(self):
        """ Returns pose of tool_frame """
        #pose_mat = utils.openraveTransformFromTo(self.robot, np.eye(4), self.tool_frame, 'world')
        #return tfx.pose(pose_mat, frame='world')
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
        self.robot.SetDOFValues(joint_values, self.arm_indices)
        
    def set_pose(self, pose):
        ref_frame = pose.frame if pose.frame is not None else "world"
        joint_values = cart_to_joint(self.manip, np.array(pose.matrix), ref_frame, "end_effector")
        if joint_values is not None:
            self.set_joint_values(joint_values)
        else:
            print('pose_matrix invalid, ignoring command')
            
    def teleop(self):
        pos_step = .01
        delta_position = {'a' : [pos_step, 0, 0],
                          'd' : [-pos_step, 0, 0],
                          'w' : [0, pos_step, 0],
                          'x' : [0, -pos_step, 0],
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
                
                #p = tfx.pose(delta_position[char], frame=self.manip.GetEndEffector().GetName())
                #p = tfx.pose([0,0,0], frame=self.tool_frame)
                #new_mat = utils.openraveTransformFromTo(self.robot, p.matrix, self.tool_frame, 'world')
                #new_pose = tfx.pose(new_mat, frame='world')
                #print(pose)
                #print(new_pose)
                #IPython.embed()
                
                #p = tfx.pose(transform_relative_pose_for_ik(self.manip, np.array(new_pose.matrix), new_pose.frame, 'end_effector'))
                #h = utils.plot_transform(self.robot.GetEnv(), p.matrix)
                #IPython.embed()
                #pose += delta_position[char]
            ypr = np.array([pose.tb_angles.yaw_deg, pose.tb_angles.pitch_deg, pose.tb_angles.roll_deg])
            if delta_angle.has_key(char):
                ypr += np.array(delta_angle[char])
                new_pose = tfx.pose(pose.position, tfx.tb_angles(ypr[0], ypr[1], ypr[2]))
            #angle = tfx.tb_angles(ypr[0], ypr[1], ypr[2])
            #new_pose = tfx.pose(position, angle, frame=position.frame)
            self.set_pose(new_pose)    
            time.sleep(.01)
            
class Head:
    joint_names = ['head_pan_joint', 'head_tilt_joint']
    def __init__(self, pr2):
        self.robot = pr2.robot
        self.joint_indices = [self.robot.GetJointIndex(joint_name) for joint_name in Head.joint_names]
        
    def get_limits(self):
        l_limits, u_limits = self.robot.GetDOFLimits()
        head_l_limits = np.array([l_limits[i] for i in self.joint_indices])
        head_u_limits = np.array([u_limits[i] for i in self.joint_indices])
        return head_l_limits, head_u_limits
        
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
        Only for points in FOV
        
        Returns
        [(point, pixel, color), ...]
        """
        data = self.get_data()
        P = data.KK
        
        points_pixels_colors = list()
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

            y_max, x_max, _ = data.imagedata.shape
            if 0 <= y_pixel < y_max and 0 <= x_pixel < x_max:        
                color = data.imagedata[y_pixel, x_pixel]
                points_pixels_colors.append((x, np.array([y_pixel, x_pixel]), np.array(color, dtype=float)))
            
        return points_pixels_colors
    
class KinectSensor:
    """
    Simulates a kinect by having a depth sensor + camera sensor
    """
    def __init__(self, depth_sensor, camera_sensor):
        self.depth_sensor = depth_sensor
        self.camera_sensor = camera_sensor
        
    def power_on(self):
        self.depth_sensor.power_on()
        self.camera_sensor.power_on()
        
    def power_off(self):
        self.depth_sensor.power_off()
        self.camera_sensor.power_off()
        
    def get_data(self):
        """ Returns a PointCloud, i.e. a list of ColoredPoint """
        points = self.depth_sensor.get_points()
        points_pixels_colors = self.camera_sensor.get_pixels_and_colors(points)
        
        colored_points = list()
        for point, pixel, color in points_pixels_colors:
            colored_points.append(ColoredPoint(point, color))
            
        return colored_points
                
    def render_on(self):
        self.depth_sensor.render_on()
        self.camera_sensor.render_on()
        
    def render_off(self):
        self.depth_sensor.render_off()
        self.camera_sensor.render_off()
        
        
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
    brett = PR2('../envs/pr2-test.env.xml')
    env = brett.env
    env.SetViewer('qtcoin')
    
    larm = brett.larm
    rarm = brett.rarm
    head = brett.head
    kinect = brett.kinect_sensor
    
    larm.set_posture('side')
    rarm.set_posture('side')
    
    rarm.set_pose(tfx.pose((2.8011817158976595, -1.882240568190215, 0.9784708528906445),(-0.2167526477835155, 0.7996943738474431, 0.1764675951856764, 0.5313815822598347),frame='world'))
    p = rarm.get_pose()
    h = utils.plot_transform(env, p.matrix)
    time.sleep(1)
    
    ee_mat = rarm.manip.GetEndEffectorTransform()
    
    #name = 'r_gripper_tool_frame'
    #p_mat = utils.openraveTransformFromTo(brett.robot, np.eye(4), name, 'world')
    #h = utils.plot_transform(env, p_mat)
    IPython.embed()
    return
    
    
    kinect.power_on()
    kinect.render_on()
    time.sleep(1)
    colored_points = kinect.get_data()

    handles = list()    
    for colored_point in colored_points:
        handles.append(colored_point.display(env))
        
    image = kinect.camera_sensor.get_data().imagedata
    points_pixels_colors = kinect.camera_sensor.get_pixels_and_colors(kinect.depth_sensor.get_points())
    for point, pixel, color in points_pixels_colors:
        for i in xrange(-5,5):
            image[pixel[0],pixel[1]+i] = [0, 255, 0]
        for i in xrange(-5,5):
            image[pixel[0]+i,pixel[1]] = [0, 255, 0]
    
    plt.imshow(image)
    plt.show(block=False)
    IPython.embed()
    
if __name__ == '__main__':
    test()