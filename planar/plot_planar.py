import random
import os
import time

import numpy as np

import matplotlib.pyplot as plt

import IPython

def plot_planar_gmm(J, obj_means, obj_covs, obj_particles, robot_origin, link_lengths, camera_origin, camera_fov, object, beams, pause=True):
    """
    @param J                joints (robot + camera)
    @param obj_means        object means
    @param obj_covs         object covariances
    @param obj_particles    particle groups for each object
    @param robot_origin     (x,y) position of link0
    @param link_lengths
    @param camera_origin     (x,y) position of camera
    @param camera_fov angle  (in radians) of camera field-of-view
    @param object            actual position of the object
    @param beams             list of size (3,) arrays of triangles representing the FOV
    """
    plt.clf()
    plt.cla()
    
    fig = plt.figure(1)
    ax = fig.axes[0]
    ax.set_axis_bgcolor('black')
    ax.set_aspect('equal')
    
    T = len(J)
    
    max_length = np.sum(link_lengths)
    ax.axis([-max_length, max_length, -2, max_length])
    
    for t in xrange(T):
        """ plot robot arm """
        plot_arm(J[t][:3], robot_origin, link_lengths, alpha=(t+1)/float(T))
    
        """plot camera angle """
        plot_camera(J[t][3], camera_origin, camera_fov, alpha=(t+1)/float(T))
    
    """ plot object(s) mean, covariance, and particles """
    if obj_means is not None and obj_covs is not None and obj_particles is not None:
        num_objs = len(obj_means)
        colors = list(plt.cm.rainbow(np.linspace(0, 1, num_objs)))
        for obj_mean, obj_cov, obj_p, color in zip(obj_means, obj_covs, obj_particles, tuple(colors)):
            plt.plot(obj_mean[0], obj_mean[1], 's', markersize=10.0, color=color)
            plot_cov(obj_mean, obj_cov, color=color)
            
            plt.plot(obj_p[0,:], obj_p[1,:], 'x', color=color)
                
        
    """ plot object position """
    plt.plot(object[0], object[1], 's', markersize=10.0, color='green')
    
    """ plot beams for last timestep only """
    plot_beams(beams)
    
    plt.show(block=False)
    plt.pause(.1)
    
    save('../figures/full_runs/fig{0}'.format(time.time()), ext="png", close=False, verbose=False)
    
    if (pause):
        print('Press enter to continue')
        raw_input()

def plot_planar(J, robot_origin, link_lengths, camera_origin, camera_fov, object, beams, pause=True):
    """
    @param J                joints (robot + camera)
    @param robot_origin     (x,y) position of link0
    @param link_lengths
    @param camera_origin     (x,y) position of camera
    @param camera_fov angle  (in radians) of camera field-of-view
    @param object            actual position of the object
    @param beams             list of size (3,) arrays of triangles representing the FOV
    """
    plot_planar_gmm(J, None, None, None, robot_origin, link_lengths, camera_origin, camera_fov, object, beams, pause)


def plot_planar_old(X, S, P, P_sd, robot_origin, link_lengths, camera_origin, camera_fov, object, beams, pause=False):
    """
    @param X                 joints angles, camera angle, object position (for T timesteps)
    @param S                 covariance for X (for T timesteps)
    @param P                 particles (or NULL if no particles)
    @param P_sd              particles signed distance (or NULL if no particles)
    @param robot_origin      (x,y) position of link0
    @param link_lengths
    @param camera_origin     (x,y) position of camera
    @param camera_fov angle  (in radians) of camera field-of-view
    @param object            actual position of the object
    @param beams             list of size (3,) arrays of triangles representing the FOV
    """
    plt.clf()
    plt.cla()
    
    fig = plt.figure(1)
    ax = fig.axes[0]
    ax.set_axis_bgcolor('black')
    ax.set_aspect('equal')
    
    X_DIM, T = X[0].shape[0], len(X)
    
    max_length = np.sum(link_lengths)
    ax.axis([-max_length, max_length, -2, max_length])
    
    #colors = plt.cm.RdYlGn(np.linspace(0, 1, T))
    for t in xrange(T):
        """ plot robot arm """
        plot_arm(X[t][:3], robot_origin, link_lengths, alpha=(t+1)/float(T))
    
        """ plot object mean and covariance """
        mu = X[t][-2:]
        sigma = S[t][-2:, -2:]
        
        color = 'yellow'
        alpha = (t+1)/float(T)
        plt.plot(mu[0], mu[1], 's', markersize=10.0, color=color, alpha=alpha)
        plot_cov(mu, sigma, color=color, alpha=alpha)
        
        """plot camera angle """
        plot_camera(X[t][3], camera_origin, camera_fov, alpha=alpha)
        
    """ plot object position """
    plt.plot(object[0], object[1], 's', markersize=10.0, color='green')
    
    """ plot beams for last timestep only """
    plot_beams(beams)
    
    """ plot particles (if P is ndarray) """
    if type(P) is np.ndarray:
        M_DIM = P.shape[1]
        ee_pos = np.reshape(arm_joint_positions(X[-1][:3], robot_origin, link_lengths)[:,-1], (2,1))
        P_ee_dist = np.linalg.norm(P - np.repeat(ee_pos, M_DIM, axis=1), axis=0)
        #initialize_gaussian(P, P_sd, P_ee_dist)
        
        for m in xrange(P.shape[1]):
            plt.plot(P[0,m], P[1,m], 'x', color='yellow')
    
    
    plt.show(block=False)
    
    if (pause):
        print('Press enter to continue')
        raw_input()
    
def plot_arm(joints, robot_origin, link_lengths, color='red', alpha=1.0):
    joint_positions = arm_joint_positions(joints, robot_origin, link_lengths)
    for i in xrange(joint_positions.shape[1] - 1):
        plt.plot(joint_positions[0,i:i+2], joint_positions[1,i:i+2], '-o', color=color, alpha=alpha, markersize=10.0, linewidth=3.0)
        
def arm_joint_positions(joints, robot_origin, link_lengths):
    """ return 2d position of each link, starting from robot_origin """
    joint_positions = np.zeros((2, 1+link_lengths.shape[0]))
    joint_positions[:,0] = robot_origin
    
    x0, y0 = tuple(robot_origin)
    joint0 = 0
    
    index = 1
    for joint, link_length in zip(tuple(joints), tuple(link_lengths)):
        joint1 = joint0 + joint 
        x1 = x0 + np.sin(joint1)*link_length
        y1 = y0 + np.cos(joint1)*link_length
        
        joint_positions[:,index] = np.array([x1,y1])
        index += 1
        
        x0, y0, joint0 = x1, y1, joint1
        
    return joint_positions

    
    
def plot_cov(mu, sigma, color = 'yellow', alpha=1):
    t = np.linspace(-np.pi, np.pi, 2*np.pi/.1)
    x = np.sin(t)
    y = np.cos(t)

    D_diag, V = np.linalg.eigh(sigma)
    D = np.diag(D_diag)
    A = np.real((np.dot(V,np.sqrt(D))).T)

    z = np.dot(np.vstack((x.T,y.T)).T,A)

    plt.plot(z[:,0]+mu[0], z[:,1]+mu[1], color=color, alpha=alpha)
    
def plot_camera(angle, center, fov, color='blue', alpha=1.0):
    x_left = center[0] + 15*np.sin(angle - fov/2.0)
    y_left = center[1] + 15*np.cos(angle - fov/2.0)
    
    x_right = center[0] + 15*np.sin(angle + fov/2.0)
    y_right = center[1] + 15*np.cos(angle + fov/2.0)
    
    plt.plot([center[0], x_left], [center[1], y_left], linewidth=3.0, color=color, alpha=alpha)
    plt.plot([center[0], x_right], [center[1], y_right], linewidth=3.0, color=color, alpha=alpha)

def plot_beams(beams, debug=False):
    if debug:
        plt.clf()
        plt.cla()
        fig = plt.figure(1)
        ax = fig.add_subplot(111, axisbg='black')
    
    for beam in beams:
        plot_beam(beam)
        
    if debug:
        plt.show(block=False)
        print('Press enter to continue')
        raw_input()

def plot_beam(beam):
    """
    @param beam 2 by 3 ndarray of triangle points
    """
    #beam_closed = np.hstack((beam, beam[:,0].reshape(2,1)))
    #plt.plot(beam_closed[0,:], beam_closed[1,:])
    
    x = beam[0,2], beam[0,0], beam[0,1]
    y = beam[1,2], beam[1,0], beam[1,1]
    plt.fill(x, y, color='white', edgecolor='none', alpha=.5)
    
    
def save(path, ext='png', close=True, verbose=True):
    """Save a figure from pyplot.
     
    Parameters
    ----------
    path : string
    The path (and filename, without the extension) to save the
    figure to.
     
    ext : string (default='png')
    The file extension. This must be supported by the active
    matplotlib backend (see matplotlib.backends module). Most
    backends support 'png', 'pdf', 'ps', 'eps', and 'svg'.
     
    close : boolean (default=True)
    Whether to close the figure after saving. If you want to save
    the figure multiple times (e.g., to multiple formats), you
    should NOT close it in between saves or you will have to
    re-plot it.
     
    verbose : boolean (default=True)
    Whether to print information about when and where the image
    has been saved.
     
    """
    # Extract the directory and filename from the given path
    directory = os.path.split(path)[0]
    filename = "%s.%s" % (os.path.split(path)[1], ext)
    if directory == '':
        directory = '.'
     
    # If the directory does not exist, create it
    if not os.path.exists(directory):
        os.makedirs(directory)
     
    # The final path to save to
    savepath = os.path.join(directory, filename)
     
    if verbose:
        print("Saving figure to '%s'..." % savepath),
     
    # Actually save the figure
    plt.savefig(savepath)
    # Close it
    if close:
        plt.close()
     
    if verbose:
        print("Done")
    
