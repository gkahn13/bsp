import IPython

import random

from pr2 import pr2_sim
from pr2 import utils

import roslib
roslib.load_manifest('tfx')
import tfx

UNKNOWN = -1

class EihSystem:
    def __init__(self, manip, kinect):
        """
        manip -- Head or Arm
           must have methods get_joint_values, set_joint_values, get_limits
        """
        self.manip = manip
        self.kinect = kinect
        
        self.kinect.power_on()
        
        self.DT = 1.0
        self.X_DIM = len(self.manip.get_joint_values())
        
        self.desired_observations = [(1, UNKNOWN, UNKNOWN), (1, 1, UNKNOWN), (1, 0, np.array([255.,0,0]))]
        
    def dynfunc(self, x, u):
        x_new = np.array(x + self.DT*u)
        
        l_limits, u_limits = self.manip.get_limits()
        x_new = np.maximum(x_new, l_limits)
        x_new = np.minimum(x_new, u_limits)
        
        return x_new
    
    def obsfunc(self, x, particle):
        """
        x -- current state
        particle -- tfx.point 
        """
        is_in_fov, is_occluded, color = UNKNOWN, UNKNOWN, UNKNOWN
        
        pixel = self.kinect.get_pixel_from_point(particle)
        is_in_fov = 1 if pixel is not None else 0
        
        if is_in_fov == 1:
            particle_dist = self.kinect.distance_to(particle)
            
            z_buffer = self.kinect.get_z_buffer()
            y, x = pixel
            
            is_occluded = 0 if particle_dist - .01 < z_buffer[y,x] else 1
            
            image = self.kinect.get_image()
            color = image[y, x]
        
        return (is_in_fov, is_occluded, color)

def random_within(lower, upper):
    return random.random()*(upper - lower) + lower

def test_eih_system():
    brett = pr2_sim.PR2('envs/pr2-test.env.xml')
    env = brett.env
    
    larm = brett.larm
    rarm = brett.rarm
    l_kinect = brett.l_kinect
    r_kinect = brett.r_kinect
    
    larm.set_posture('mantis')
    rarm.set_posture('mantis')
    
    table = env.GetKinBody('table')
    base = table.GetLink('base')
    extents = base.Geometry.GetBoxExtents(base.GetGeometries()[0])
    
    table_pose = tfx.pose(table.GetTransform())
    # assume table has orientation np.eye(3)
    x_min, x_max = table_pose.position.x - extents[0], table_pose.position.x + extents[0]
    y_min, y_max = table_pose.position.y - extents[1], table_pose.position.y + extents[1]
    z_min, z_max = table_pose.position.z + extents[2], table_pose.position.z + extents[2] + .2
    
    M = 1000
    particles = list()
    for i in xrange(M):
        x = random_within(x_min, x_max)
        y = random_within(y_min, y_max)
        z = random_within(z_min, z_max)
        particle = tfx.point([x,y,z])
        particles.append(particle)
        
    handles = list()
    for p in particles:
        handles.append(utils.plot_point(env, p.array, color=[1.,0,0]))
    
    IPython.embed()

if __name__ == '__main__':
    test_eih_system()