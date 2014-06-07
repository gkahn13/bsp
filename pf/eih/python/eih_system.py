import IPython

import sys
import random
import time
import colorsys

import numpy as np
from matplotlib import pyplot as plt

from pr2 import pr2_sim
from pr2 import utils

import openravepy as rave
import roslib
roslib.load_manifest('tfx')
import tfx
from tfx import transformations as tft

UNKNOWN = -1e6

class EihSystem:
    def __init__(self, env, manip, kinect, obs_type='fov_occluded_color'):
        """
        manip -- Head or Arm
           must have methods get_joint_values, set_joint_values, get_limits
        obs_type -- 'fov_occluded_color'
        """
        self.env = env
        self.manip = manip
        self.kinect = kinect
        
        height, width = self.kinect.height, self.kinect.width
        min_range, max_range, optimal_range = self.kinect.min_range, self.kinect.max_range, self.kinect.optimal_range
        self.obs_type = obs_type
        if obs_type == 'fov':
            self.obsfunc_discrete_weight = self.obsfunc_discrete_weight_fov
            self.obsfunc_continuous_weight = self.obsfunc_continuous_weight_fov
        elif obs_type == 'fov_occluded':
            self.obsfunc_discrete_weight = self.obsfunc_discrete_weight_fov_occluded
            self.obsfunc_continuous_weight = self.obsfunc_continuous_weight_fov_occluded
            
            self.fov_occluded_z_d = [np.array([height/2., width/2., optimal_range])]
            self.R = np.diag([0])
        else:
            self.obs_type = 'fov_occluded_color'
            self.obsfunc_discrete_weight = self.obsfunc_discrete_weight_fov_occluded_color
            self.obsfunc_continuous_weight = self.obsfunc_continuous_weight_fov_occluded_color
            
            self.desired_observations = [np.array((0, UNKNOWN, UNKNOWN, UNKNOWN, UNKNOWN)),
                                     np.array((1, 1, UNKNOWN, UNKNOWN, UNKNOWN)),
                                     np.array((1, .5) + colorsys.rgb_to_hsv(1., 0, 0))]
            
            self.R = np.diag([.5, .01, .05, 1, 1])
        
        self.kinect.power_on()
        time.sleep(1)
        
        self.DT = 1.0
        self.X_DIM = len(self.manip.get_joint_values())
        self.U_DIM = self.X_DIM
        
        self.handles = list()
        
    def dynfunc(self, x, u):
        x_new = np.array(x + self.DT*u)
        
        l_limits, u_limits = self.manip.get_limits()
        x_new = np.maximum(x_new, l_limits)
        x_new = np.minimum(x_new, u_limits)
        
        return x_new
    
    def obsfunc_discrete_weight_fov(self, particle, image, z_buffer):
        if not self.kinect.is_visible(particle):
            return 1, (0, 1, 0)
        else:
            return 0, (1, 1, 1)
    
    def obsfunc_continuous_weight_fov(self, particle, image, z_buffer):
        y, x = self.kinect.get_pixel_from_point(particle, False)
        d = self.kinect.distance_to(particle)
        
        z = np.array([y, x, d])
        
        height, width = self.kinect.height, self.kinect.width
        min_range, max_range, optimal_range = self.kinect.min_range, self.kinect.max_range, self.kinect.optimal_range
        
        R = np.diag([100*height, 100*width, 100*(max_range - min_range)])
        max_likelihood = self.gauss_likelihood(np.zeros((R.shape[0],1)), R)
        
        z_low_weight = np.array([height/2., width/2., optimal_range])
        
        weight = (max_likelihood - self.gauss_likelihood(z - z_low_weight, R)) / max_likelihood
        
        color = (1, 1, 1) if weight < .5 else (0, 1, 0)
        
        return weight, color
    
    def obsfunc_discrete_weight_fov_occluded(self, particle, image, z_buffer):
        if not self.kinect.is_visible(particle):
            return 1, (0, 1, 0)
        else:
            particle_dist = self.kinect.distance_to(particle)
            y, x = self.kinect.get_pixel_from_point(particle)
            sd = particle_dist - z_buffer[y,x]
            
            if sd > -.01:
                return 1, (0, 0, 0)
            else:
                return 0, (1, 1, 1)
    
    def obsfunc_continuous_weight_fov_occluded(self, particle, image, z_buffer):
        fov_weight, color = self.obsfunc_continuous_weight_fov(particle, image, z_buffer)
        
        if self.kinect.is_visible(particle):
            particle_dist = self.kinect.distance_to(particle)
            y, x = self.kinect.get_pixel_from_point(particle)
            sd = particle_dist - z_buffer[y,x]
        
            occluded_weight = sigmoid(sd, 5)
            
            weight = max(fov_weight, occluded_weight)
            
            if occluded_weight < .5:
                color = (1, 1, 1)
            else:
                color = (0, 0, 0)
        else:
            weight = fov_weight
            color = (0, 1, 0)
            
        return weight, color
        
    
    def obsfunc_discrete_weight_fov_occluded_color(self, particle, image, z_buffer):
        """ Returns weight, color to display """
        #if not self.kinect.is_visible(particle):
        if not self.kinect.camera_sensor.is_in_fov(particle):
            return 1.0, (0, 1, 0)
        
        particle_dist = self.kinect.distance_to(particle)
        y, x = self.kinect.get_pixel_from_point(particle)
        sd = particle_dist - z_buffer[y,x]
        
        if sd > .03:
            return 1.0, (0, 0, 0)
        elif sd < -.03:
            return 0.0, (1, 1, 1)
        
        r, g, b = tuple(image[y,x]/255.)
        h, s, v = colorsys.rgb_to_hsv(r, g, b)
        
        dist_to_red = min(h, 1 - h)
        if dist_to_red < .04:
            return 1.1, (r, g, b)
        else:
            return 0.0, (1, 1, 0)
        
        
    
    def obsfunc_continuous_weight_fov_occluded_color(self, particle, image, z_buffer):
        """
        x -- current state
        particle -- tfx.point 
        """
        is_in_fov, sd_sigmoid, color = UNKNOWN, UNKNOWN, (UNKNOWN, UNKNOWN, UNKNOWN)
        
        # exact_is_in_fov = 1 if self.kinect.is_in_fov(particle) else .5
        
        # sigmoid approximation to check if in FOV
        y, x = self.kinect.get_pixel_from_point(particle)
        boundary = np.array([y, -y + self.kinect.height, x, -x + self.kinect.width])
        is_in_fov = np.prod(sigmoid(boundary, 1e6))
        
        is_in_fov = (1-.5)*is_in_fov + .5 # make being out of FOV not so advantageous for gauss likelihood
        
        # higher it is, more it weights particles inside of FOV
        # [.5, 1]
        # .75 best so far
        if is_in_fov > .75:
            particle_dist = self.kinect.distance_to(particle)
            
            y, x = self.kinect.get_pixel_from_point(particle, False)
            
            sd = particle_dist - z_buffer[y,x]
            sd_sigmoid = sigmoid(sd, 10)
            
            if abs(sd) < .03:
                r, g, b = tuple(image[y, x]/255.)
                color = colorsys.rgb_to_hsv(r, g, b)
                
            sd_sigmoid = (1-.1)*sd_sigmoid + .1/2.0
            
        z = np.array((is_in_fov, sd_sigmoid) + color)
        
        weight = -np.inf
        for z_d in self.desired_observations:
            e = z - z_d
            # special hue case, make it wrap around
            hue_small, hue_large = min(z[2], z_d[2]), max(z[2], z_d[2])
            e[2] = min(hue_large - hue_small, hue_small + (1-hue_large))
            
            weight = max(weight, self.gauss_likelihood(e, self.R))
        
        h, s, v = z[2:]
        r, g, b = colorsys.hsv_to_rgb(h, s, v)
        if z[0] <= .5: # out of fov
            color = np.array((0,1,0))
        elif z[1] < .25: # free space
            color = np.array((1,1,1))
        elif z[1] > .75: # occluded
            color = np.array((0,0,0))
        else: # near surface
            color = np.array((r,g,b))
        
        return weight, color
    
    def update_state_and_particles(self, x_t, particles_t, u_t, plot=False, add_noise=True):
        self.handles = list()
        M = len(particles_t)
        
        x_tp1 = self.dynfunc(x_t, u_t)
        self.manip.set_joint_values(x_tp1)
        
        image = self.kinect.get_image()
        z_buffer = self.kinect.get_z_buffer()
        
        W = np.zeros(M)
        for m in xrange(M):    
            W[m], color = self.obsfunc_discrete_weight(particles_t[m], image, z_buffer)
            #W[m], color = self.obsfunc_continuous_weight(particles_t[m], image, z_buffer)
            if plot:
                self.handles.append(utils.plot_point(self.env, particles_t[m].array, color=color))
        
        W = W / np.sum(W)
        
        sampling_noise = random_within(0, 1.0/float(M))
        particles_tp1 = self.low_variance_sampler(particles_t, W, sampling_noise)
        
        if add_noise:
            for i in xrange(len(particles_tp1)):
                noise = [random_within(-.005, .005) for _ in xrange(3)]
                particles_tp1[i] = particles_tp1[i] + noise
        
        return x_tp1, particles_tp1
            
    def gauss_likelihood(self, v, S):
        Sf = np.linalg.cholesky(S)
        M = np.linalg.solve(Sf, v)
        
        E = -0.5*np.sum(M*M)
        C = np.power(2*np.pi, S.shape[1]/2.) * np.prod(np.diag(Sf))
        w = np.exp(E) / C
        
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
            particles_sampled.append(particles[i].copy())
            
        return particles_sampled
    
    def cost(self, x0, U, particles, use_discrete=False):
        """
        Starting at x0, propagates weights by each
        u in U, calculating entropy of the particle weights
        at each time step
        """
        M = len(particles)
        T = len(U)
        entropy = 0.
        
        x_t = x0
        W = (1./float(M))*np.ones(M)
        for t in xrange(T):
            u_t = U[t]
            x_tp1 = self.dynfunc(x_t, u_t)
            self.manip.set_joint_values(x_tp1)
            
            image = self.kinect.get_image()
            z_buffer = self.kinect.get_z_buffer()
            
            for m in xrange(M):
                if use_discrete:
                    W_m, color = self.obsfunc_discrete_weight(particles[m], image, z_buffer)
                else:
                    W_m, color = self.obsfunc_continuous_weight(particles[m], image, z_buffer)
                    
                W[m] *= W_m
            W = W / np.sum(W)
        
            for w in W.tolist():
                if w != 0:
                    entropy += -w*np.log(w)
                    
        self.manip.set_joint_values(x0)
                    
        return entropy
    
    def cost_grad(self, x0, U, particles, step=1e-5, use_discrete=False):
        T = len(U) + 1
        grad = [np.zeros(self.U_DIM) for _ in xrange(T-1)]
        
        for t in xrange(T-1):
            for i in xrange(self.U_DIM):
                u_orig = U[t][i]
                
                U[t][i] = u_orig + step
                cost_p = self.cost(x0, U, particles, use_discrete)
                
                U[t][i] = u_orig - step
                cost_m = self.cost(x0, U, particles, use_discrete)
                
                grad[t][i] = (cost_p - cost_m) / (2*step)
                
        return grad
                
    
    
    
            
def sigmoid(x, alpha):
    return 1.0/(1.0 + np.exp(-alpha*x))

def random_within(lower, upper):
    return random.random()*(upper - lower) + lower





"""  Test Functions """

def setup_environment(obs_type, M=1000, lr='r', zero_seed=True):
    if zero_seed:
        random.seed(0)
        
    brett = pr2_sim.PR2('../envs/pr2-test.env.xml')
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
    
    #x_min, x_max = mug_pose.position.x - .03, mug_pose.position.x + .03
    #y_min, y_max = mug_pose.position.y + .03, mug_pose.position.y + .03
    #z_min, z_max = mug_pose.position.z + extents[2], mug_pose.position.z + extents[2] + .2
    
    particles = list()
    for i in xrange(M):
        x = random_within(x_min, x_max)
        y = random_within(y_min, y_max)
        z = random_within(z_min, z_max)
        particle = tfx.point([x,y,z])
        particles.append(particle)
        
    particles = sorted(particles, key=lambda x: (x-mug_pose.position).norm)
        
    arm = larm if lr == 'l' else rarm
    kinect = l_kinect if lr == 'l' else r_kinect
        
    eih_sys = EihSystem(env, arm, kinect, obs_type)
    kinect.render_on()
    #arm.set_pose(tfx.pose([2.901, -1.712,  0.868],tfx.tb_angles(-143.0, 77.9, 172.1))) # FOR rarm ONLY
    time.sleep(1)
    
    return eih_sys, particles
    
def test_pf_update(M=1000):
    eih_sys, particles = setup_environment('fov_occluded_color', M=M, zero_seed=False)
    arm = eih_sys.manip
    kinect = eih_sys.kinect
    
    """
    arm.set_pose(tfx.pose(arm.get_pose().position, tfx.tb_angles(0, 90, 0)))
    p = arm.get_pose()
    center = [p.position.x + .75, p.position.y, p.position.z]
    
    particles = list()
    for i in xrange(M):
        pos = [0,0,0]
        for i in xrange(len(center)):
            pos[i] = random_within(center[i] - .025, center[i] + .025)
        particle = tfx.point(pos)
        particles.append(particle)
    """
    
    arm.set_pose(tfx.pose([2.901, -1.712,  0.868],tfx.tb_angles(-143.0, 77.9, 172.1))) # FOR rarm ONLY
    x_t = arm.get_joint_values()
    particles_t = particles
    u_t = np.zeros(x_t.shape[0])
    
    try:
        t = 0
        while True:
            x_tp1, particles_tp1 = eih_sys.update_state_and_particles(x_t, particles_t, u_t, True, False)
            
            particles_t = particles_tp1
            print('Iter: {0}'.format(t))
            t += 1
            
            arm.teleop()
            x_t = arm.get_joint_values()
            
            #utils.save_view(env, 'figures/eih_pf_{0}.png'.format(t))
            
    except KeyboardInterrupt:
        rave.RaveDestroy()
        print('Exiting')
        sys.exit()
    
    IPython.embed()
    
def test_cost(T = 2):
    eih_sys, particles = setup_environment('fov', zero_seed=False)
    arm = eih_sys.manip
    
    X_DIM, U_DIM = eih_sys.X_DIM, eih_sys.U_DIM

    """
    X_DES = [np.array([-1.8480293, 0.35748196, -1.011, -2.02331937, -0.84592055, -1.11009373, -3.77125935]),
             np.array([-1.77308591, 0.35872446, -1.011, -1.89298953, -0.63370048, -1.32043822, -3.97248228]),
             np.array([-1.9651356, 0.29791627, -1.111, -2.31242513, -1.59783918, -0.71440641, -3.02227944])]
    
    x0 = X_DES[0]
    U = [X_DES[t+1] - X_DES[t] for t in xrange(len(X_DES)-1)]
    """
    
    x0 = arm.get_joint_values()
    U = [np.zeros(U_DIM) for _ in xrange(T-1)]
    
    #cost = eih_sys.cost(x0, U, particles)
    #print('cost: {0}'.format(cost))
    
    print('computing grad...')
    grad = eih_sys.cost_grad(x0, U, particles, use_discrete=False)
    
    IPython.embed()
    
def test_greedy(M=1000):
    random.seed(0)
        
    brett = pr2_sim.PR2('../envs/pr2-test.env.xml')
    env = brett.env
    arm = brett.rarm
    kinect = brett.r_kinect
    
    arm.set_posture('mantis')
    p = arm.get_pose()
    arm.set_pose(tfx.pose(p.position + [.5, 0, 0], tfx.tb_angles(0, 90, 0)))
    
    # .23 just on edge of range
    # .25 out of range
    p = kinect.get_pose()
    center = [p.position.x + .3, p.position.y, p.position.z]
    
    P = list()
    for i in xrange(M):
        pos = [0,0,0]
        for i in xrange(len(center)):
            pos[i] = random_within(center[i] - .025, center[i] + .025)
        particle = tfx.point(pos)
        P.append(particle)
        
    handles = list()
    for p in P:
        handles.append(utils.plot_point(env, p.array,color=[1,0,0]))
    
                
    eih_sys = EihSystem(env, arm, kinect, 'fov')
    kinect.render_on()
    time.sleep(1)
        
    x0 = arm.get_joint_values()
    U = [np.zeros(x0.shape)]
    
    try:
        t = 0
        while True:
            print('Iter: {0}'.format(t))
            t += 1
            
            arm.set_joint_values(x0)
            
            grad = eih_sys.cost_grad(x0, U, P, step=1e-3, use_discrete=False)[0]
            print('grad: {0}'.format(list(grad)))
            
            u_grad = -(2*(np.pi/180.))*grad/np.linalg.norm(grad, 2)
            x1 = x0 + u_grad
            
            arm.set_joint_values(x0)
            p0 = arm.get_pose()
            arm.set_joint_values(x1)
            p1 = arm.get_pose()
            
            delta_pos = .05*(p1.position - p0.position)/(p1.position - p0.position).norm
            clipped_quat = tft.quaternion_slerp(p0.tb_angles.to_quaternion(), p1.tb_angles.to_quaternion(), .5)
            
            p1_clipped = tfx.pose(p0.position + delta_pos, clipped_quat)
            arm.set_pose(p1_clipped)
            x1_clipped = arm.get_joint_values()
            u_t = x1_clipped - x0
            arm.set_joint_values(x0)
            
            x_tp1, P_tp1 = eih_sys.update_state_and_particles(x0, P, u_t, plot=True, add_noise=False)
            
            P = P_tp1
            x0 = x_tp1
            
            
            #utils.save_view(env, 'figures/eih_pf_{0}.png'.format(t))
            
    except KeyboardInterrupt:
        rave.RaveDestroy()
        print('Exiting')
        sys.exit()
    

if __name__ == '__main__':
    test_pf_update()
    #test_cost()
    #test_greedy()


