import numpy as np

import matplotlib.pyplot as plt

import IPython

def plot_planar(X, S, robot_origin, link_lengths, camera_origin, camera_fov, object):
    """
    @param X                 joints angles, camera angle, object position (for T timesteps)
    @param S                 covariance for X (for T timesteps)
    @param robot_origin      (x,y) position of link0
    @param link_lengths
    @param camera_origin     (x,y) position of camera
    @param camera_fov angle  (in radians) of camera field-of-view
    @param object            actual position of the object
    """
    plt.clf()
    plt.cla()
    
    fig = plt.figure(1)
    ax = fig.add_subplot(111, axisbg='black')
    
    X_DIM, T = X[0].shape[0], len(X)
    
    max_length = np.sum(link_lengths)
    ax.axis([-max_length, max_length, -max_length, max_length])
    
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
    
    plt.show(block=False)
    
    IPython.embed()
    
def plot_arm(joints, robot_origin, link_lengths, color='red', alpha=1.0):
    x0, y0 = tuple(robot_origin)
    joint0 = 0
    
    for joint, link_length in zip(tuple(joints), tuple(link_lengths)):
        joint1 = joint0 + joint 
        x1 = x0 + np.sin(joint1)*link_length
        y1 = y0 + np.cos(joint1)*link_length
        
        plt.plot([x0, x1], [y0, y1], '-o', color=color, alpha=alpha, markersize=10.0, linewidth=3.0)
        
        x0, y0, joint0 = x1, y1, joint1
    
    
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
    x_left = center[0] + 0.1*np.sin(angle - fov/2.0)
    y_left = center[1] + 0.1*np.cos(angle - fov/2.0)
    
    x_right = center[0] + 0.1*np.sin(angle + fov/2.0)
    y_right = center[1] + 0.1*np.cos(angle + fov/2.0)
    
    plt.plot([center[0], x_left], [center[1], y_left], linewidth=3.0, color=color, alpha=alpha)
    plt.plot([center[0], x_right], [center[1], y_right], linewidth=3.0, color=color, alpha=alpha)
    
def plot_beams(beams):
    # TODO: temp
    plt.clf()
    plt.cla()
    
    for beam in beams:
        print(beam)
        print('')
        plot_beam(beam)
        
    plt.show(block=False)
    print('Press enter to continue')
    raw_input()
    
def plot_beam(beam):
    """
    @param beam 2 by 3 ndarray of triangle points
    """
    beam_closed = np.hstack((beam, beam[:,0].reshape(2,1)))
    plt.plot(beam_closed[0,:], beam_closed[1,:])
    