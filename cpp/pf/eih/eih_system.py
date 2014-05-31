import IPython

import random
import time
import colorsys

import numpy as np
from matplotlib import pyplot as plt

from pr2 import pr2_sim
from pr2 import utils

import roslib
roslib.load_manifest('tfx')
import tfx

UNKNOWN = -1

class EihSystem:
    def __init__(self, env, manip, kinect):
        """
        manip -- Head or Arm
           must have methods get_joint_values, set_joint_values, get_limits
        """
        self.env = env
        self.manip = manip
        self.kinect = kinect
        
        self.kinect.power_on()
        time.sleep(1)
        
        self.DT = 1.0
        self.X_DIM = len(self.manip.get_joint_values())
        self.U_DIM = self.X_DIM
        # is_in_fov
        # is_occluded
        # r, g, b
        self.Z_DIM = 5
        
        self.desired_observations = [np.array((0, UNKNOWN, UNKNOWN, UNKNOWN, UNKNOWN)),
                                     np.array((1, 1, UNKNOWN, UNKNOWN, UNKNOWN)),
                                     np.array((1, 0) + colorsys.rgb_to_hsv(1., 0, 0))]
        
        #self.desired_observations = [np.array((0.5, UNKNOWN, UNKNOWN, UNKNOWN, UNKNOWN)),
        #                             np.array((1, .5, UNKNOWN, UNKNOWN, UNKNOWN)),
        #                             np.array((1, 0) + colorsys.rgb_to_hsv(1., 0, 0))]
        
        self.R = np.diag([.5, .5, .2, 1, 1])
        
    def dynfunc(self, x, u):
        x_new = np.array(x + self.DT*u)
        
        l_limits, u_limits = self.manip.get_limits()
        x_new = np.maximum(x_new, l_limits)
        x_new = np.minimum(x_new, u_limits)
        
        return x_new
    
    def obsfunc(self, x, particle, z_buffer):
        """
        x -- current state
        particle -- tfx.point 
        """
        is_in_fov, is_occluded, color = UNKNOWN, UNKNOWN, (UNKNOWN, UNKNOWN, UNKNOWN)
        
        pixel = self.kinect.get_pixel_from_point(particle)
        is_in_fov = 1 if pixel is not None and z_buffer[pixel[0],pixel[1]] is not None else .5
        #if z_buffer[pixel[0],pixel[1]] is None:
        #    print 'z_buffer is None!'
        
        if is_in_fov == 1:
            particle_dist = self.kinect.distance_to(particle)
            
            #z_buffer = self.kinect.get_z_buffer()
            y, x = pixel
            
            is_occluded = 0 if particle_dist - .01 < z_buffer[y,x] else .5
            
            if is_occluded == 0:
                image = self.kinect.get_image()
                r, g, b = tuple(image[y, x]/255.)
                color = colorsys.rgb_to_hsv(r, g, b)
        
        return np.array((is_in_fov, is_occluded) + color)
    
    def update_state_and_particles(self, x_t, particles_t, u_t):
        M = len(particles_t)
        
        x_tp1 = self.dynfunc(x_t, u_t)
        self.manip.set_joint_values(x_tp1)
        
        z_buffer = self.kinect.get_z_buffer()
        
        handles = list()
        W = np.zeros(M)
        for m in xrange(M):
            z_m = self.obsfunc(x_tp1, particles_t[m], z_buffer)
            
            for z_d in self.desired_observations:
                e = z_m - z_d
                # special hue case, make it wrap around
                hue_small, hue_large = min(z_m[2], z_d[2]), max(z_m[2], z_d[2])
                e[2] = min(hue_large - hue_small, hue_small + (1-hue_large))
                
                W[m] = max(W[m], self.gauss_likelihood(e, self.R))
                #W[m] = max([self.gauss_likelihood(z_m - z_d, self.R) for z_d in self.desired_observations])
            
            h, s, v = z_m[2:]
            r, g, b = colorsys.hsv_to_rgb(h, s, v)
            if z_m[0] == 1:
                color = np.array((r,g,b))
            else:
                color = np.array((0,1,0))
            handles.append(utils.plot_point(self.env, particles_t[m].array, color=color))
            #fov_str = 'in fov' if z_m[0] == 1 else 'NOT in fov'
            #occ_str = 'occluded' if z_m[1] == 1 else 'NOT occluded'
            #color = (r,g,b)
            #print('m: {0}'.format(m))
            #print(fov_str)
            #print(occ_str)
            #print(str(color))
            #print('w: {0}\n'.format(W[m]))
            
            
            
            
        W = W / np.sum(W)
        
        sampling_noise = random_within(0, 1.0/float(M))
        particles_tp1 = self.low_variance_sampler(particles_t, W, sampling_noise)
        
        return x_tp1, particles_tp1, handles
            
    def gauss_likelihood(self, v, S):
        Sf = np.linalg.cholesky(S)
        M = np.linalg.solve(Sf, v)
        
        E = -0.5*np.sum(M*M)
        C = np.power(2*np.pi, S.shape[1]/2.) * np.prod(np.diag(Sf))
        w = np.exp(E) #since normalized anyways #np.exp(E) / C
        
        return w
    
    def low_variance_sampler(self, particles, W, r):
        M = len(particles)
        particles_sampled = list()
        
        c = W[0]
        i = 0
        for m in xrange(M):
            u = r + m* (1.0/float(M))
            while u > c:
                i += 1
                c += W[i]
            noise = [random_within(-.005, .005) for _ in xrange(3)]
            particles_sampled.append(particles[i].copy() + noise)
            
        return particles_sampled
            
        

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
    
    mug = env.GetKinBody('mug')
    mug_pose = tfx.pose(mug.GetTransform())
    #mug_pose.position.x -= .05
    #mug_pose.position.y -= .05
    h = utils.plot_transform(env, mug_pose.matrix)
    
    #x_min, x_max = mug_pose.position.x - .1, mug_pose.position.x + .1
    #y_min, y_max = mug_pose.position.y - .1, mug_pose.position.y + .1
    #z_min, z_max = mug_pose.position.z + .05, mug_pose.position.z + .2
    
    M = 1000
    particles = list()
    for i in xrange(M):
        x = random_within(x_min, x_max)
        y = random_within(y_min, y_max)
        z = random_within(z_min, z_max)
        particle = tfx.point([x,y,z])
        particles.append(particle)
        
    particles = sorted(particles, key=lambda x: (x-mug_pose.position).norm)
        
    #handles = list()
    #for p in particles:
    #    handles.append(utils.plot_point(env, p.array, color=[1.,0,0]))
        
    sys = EihSystem(env, rarm, r_kinect)
    r_kinect.render_on()
    rarm.set_pose(tfx.pose([2.901, -1.712,  0.978],tfx.tb_angles(-143.0, 67.9, 172.1)))
    time.sleep(1)
    
    """
    image = r_kinect.camera_sensor.get_data().imagedata
    points_pixels_colors = r_kinect.camera_sensor.get_pixels_and_colors(r_kinect.depth_sensor.get_points())
    for point, pixel, color in points_pixels_colors:
        image[pixel[0],pixel[1]] = [0, 255, 0]
    
    plt.imshow(image)
    plt.show(block=False)
    """
    
    x_t = rarm.get_joint_values()
    particles_t = particles
    u_t = np.zeros(x_t.shape[0])
    
    t = 0
    while True:
        x_tp1, particles_tp1, handles = sys.update_state_and_particles(x_t, particles_t, u_t)
    
        #handles = list()
        #for p in particles_tp1:
        #    handles.append(utils.plot_point(env, p.array, color=[1.,0,0]))
            
        particles_t = particles_tp1
        print('Iter: {0}'.format(t))
        rarm.teleop()
        handles = None
        x_t = rarm.get_joint_values()
        t += 1
    
    IPython.embed()

if __name__ == '__main__':
    test_eih_system()