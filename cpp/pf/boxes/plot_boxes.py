import numpy as np
from numpy import matlib as ml
import math
import time

import matplotlib
import matplotlib.pyplot as plt
from scipy.interpolate import spline
import scipy.interpolate as interp


import IPython


def plot_state_and_particles(x, particles, boxes, xDim, M):
    plt.clf()
    plt.cla()
    
    plt.axis([-3, 3, -3, 3])
    
    pDim = len(particles)/M
    P = np.matrix(particles)
    P = np.reshape(P, (pDim,M))
    
    X = np.matrix(x)
    X = np.reshape(X, (xDim, len(x)/(xDim)))
    
    N = int(len(boxes)/4)
    colors = plt.cm.RdYlGn(np.linspace(0, 1, N))
    index = 0
    for n in xrange(N):
        plt.plot(P[index,:],P[index+1,:], color=colors[n], marker='o', markersize=10, mew=2)
        
        x, y = boxes[index], boxes[index+1]
        width = boxes[index+2]
        height = boxes[index+3]
        corners = np.zeros((4,2))
        
        corners = np.array([[x-width/2.0, y-height/2.0],
                            [x-width/2.0, y+height/2.0],
                            [x+width/2.0, y+height/2.0],
                            [x+width/2.0, y-height/2.0],
                            [x-width/2.0, y-height/2.0]])
        
        plt.plot(corners[:,0], corners[:,1], color=colors[n])
        plt.plot(x, y, color=colors[n], marker='^', markersize=10, mew=2.5)
        
        index += 4
        
    plt.plot(X[0,:], X[1,:], color='blue', marker='x', markersize=5, mew=1)
        
    plt.show(block=False)
    plt.pause(.05)
    