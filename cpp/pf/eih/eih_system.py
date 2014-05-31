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

UNKNOWN = -1e6

"""
TODO LIST:
- occluded particles getting weighted too much
"""

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
        """
        is_in_fov
        sd_sigmoid
          sd = particle - z_buffer
          0 -- free space
          .5 -- near surface
          1 -- occluded
        h, s, v
        """
        self.Z_DIM = 5
        
        self.desired_observations = [np.array((0, UNKNOWN, UNKNOWN, UNKNOWN, UNKNOWN)),
                                     np.array((1, 1, UNKNOWN, UNKNOWN, UNKNOWN)),
                                     np.array((1, .5) + colorsys.rgb_to_hsv(1., 0, 0))]
        
        #self.desired_observations = [np.array((0.5, UNKNOWN, UNKNOWN, UNKNOWN, UNKNOWN)),
        #                             np.array((1, .5, UNKNOWN, UNKNOWN, UNKNOWN)),
        #                             np.array((1, 0) + colorsys.rgb_to_hsv(1., 0, 0))]
        
        self.R = np.diag([.5, .01, .1, 1, 1])
        
        # higher it is, more it weights particles inside of FOV
        # [.5, 1]
        self.exploitation = .75
        
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
        is_in_fov, sd_sigmoid, color = UNKNOWN, UNKNOWN, (UNKNOWN, UNKNOWN, UNKNOWN)
        
        exact_is_in_fov = 1 if self.kinect.is_in_fov(particle) else .5
        # sigmoid approximation to check if in FOV
        y, x = self.kinect.get_pixel_from_point(particle)
        alpha = 1e6
        h_l_fov = sigmoid(y, alpha)
        h_u_fov = sigmoid(-y + self.kinect.height, alpha)
        w_l_fov = sigmoid(x, alpha)
        w_u_fov = sigmoid(-x + self.kinect.width, alpha)
        is_in_fov = h_l_fov*h_u_fov*w_l_fov*w_u_fov
        
        is_in_fov = .5*is_in_fov + .5 # make being out of FOV not so advantageous for gauss likelihood
        
        # TEMP
        #is_in_fov = exact_is_in_fov
        
        if is_in_fov > self.exploitation:
            if exact_is_in_fov != 1:
                print('CRAPPPPPPPPPPPPPPP')
            particle_dist = self.kinect.distance_to(particle)
            
            y, x = self.kinect.get_pixel_from_point(particle)
            
            sd = particle_dist - z_buffer[y,x]
            sd_sigmoid = sigmoid(sd, 1e1)
            
            if abs(sd) < .03:
                image = self.kinect.get_image()
                r, g, b = tuple(image[y, x]/255.)
                color = colorsys.rgb_to_hsv(r, g, b)
                
            #sd = sigmoid(-sd, 1e4) * sd
        
        return np.array((is_in_fov, sd_sigmoid) + color)
    
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
            if z_m[0] <= .5: # out of fov
                color = np.array((0,1,0))
            elif z_m[1] < .25: # free space
                color = np.array((1,1,1))
            elif z_m[1] > .75: # occluded
                color = np.array((0,0,0))
            else: # near surface
                color = np.array((r,g,b))
            handles.append(utils.plot_point(self.env, particles_t[m].array, color=color))
            
            """
            fov_str = 'in fov' if z_m[0] == 1 else 'NOT in fov'
            sd_str = z_m[1]
            color = (r,g,b)
            print('m: {0}'.format(m))
            print(fov_str)
            print('signed_distance: {0}'.format(sd_str))
            print(str(color))
            print('w: {0}\n'.format(W[m]))
            raw_input()
            handles = list()
            """
            
            
            
            
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
            
def sigmoid(x, alpha):
    return 1.0/(1.0 + np.exp(-alpha*x))

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
    rarm.set_pose(tfx.pose([2.901, -1.712,  0.868],tfx.tb_angles(-143.0, 77.9, 172.1)))
    time.sleep(1)
    
    """
    image = r_kinect.camera_sensor.get_data().imagedata
    points_pixels_colors = r_kinect.camera_sensor.get_pixels_and_colors(r_kinect.depth_sensor.get_points())
    for point, pixel, color in points_pixels_colors:
        image[pixel[0],pixel[1]] = [0, 255, 0]
    
    plt.imshow(image)
    plt.show(block=False)
    """
    
    rarm.teleop()
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
        
        #rarm.teleop()
        x_t = rarm.get_joint_values()
        
        utils.save_view(env, 'figures/eih_pf_{0}.png'.format(t))
        
        t += 1
        handles = None
    
    IPython.embed()

if __name__ == '__main__':
    test_eih_system()