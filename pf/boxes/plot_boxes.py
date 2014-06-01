import numpy as np
from numpy import matlib as ml
import math
import time

import matplotlib
import matplotlib.pyplot as plt
from scipy.interpolate import spline
import scipy.interpolate as interp


import IPython


def plot_state_and_particles(x, particles, box_centers, box_dims, xDim, M):
    plt.clf()
    plt.cla()
    
    #plt.axis([-2, 2, -2, 2])
    plt.axis([-5, 5, -5, 5])
    
    pDim = len(particles)/M
    P = np.matrix(particles)
    P = np.reshape(P, (pDim,M))
    
    X = np.matrix(x)
    X = np.reshape(X, (xDim, len(x)/(xDim)))
    
    N = int(len(box_centers)/xDim)
    colors = plt.cm.RdYlGn(np.linspace(0, 1, N))
    for n in xrange(N):
        plt.plot(P[n*xDim,:],P[n*xDim+1,:], color=colors[n], marker='x', markersize=5, mew=1)
        
        x, y = box_centers[n*xDim], box_centers[n*xDim+1]
        width = box_dims[n*xDim]
        height = box_dims[n*xDim+1]
        corners = np.zeros((4,2))
        
        corners = np.array([[x-width/2.0, y-height/2.0],
                            [x-width/2.0, y+height/2.0],
                            [x+width/2.0, y+height/2.0],
                            [x+width/2.0, y-height/2.0],
                            [x-width/2.0, y-height/2.0]])
        
        plt.plot(corners[:,0], corners[:,1], color=colors[n])
        plt.plot(x, y, color=colors[n], marker='^', markersize=10, mew=2.5)
        
    plt.plot(X[0,:], X[1,:], color='blue', marker='o', markersize=10, mew=2)
        
    plt.show(block=False)
    plt.pause(.05)
    