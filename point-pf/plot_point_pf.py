import numpy as np
from numpy import matlib as ml
import math
import time

import matplotlib
import matplotlib.pyplot as plt
from scipy.interpolate import spline
import scipy.interpolate as interp


import IPython


def plot_particles(particle_list, T, xDim, M):
    plt.clf()
    plt.cla()
    
    P_t = process_particle_list(particle_list, T, xDim, M)
    
    plot_domain()
    
    colors = plt.cm.RdYlGn(np.linspace(0, 1, len(P_t)))
    for P, c in zip(P_t, colors):
        plt.plot(P[0,:],P[1,:], ls='None', color=c, marker='x', markersize=10, mew=2)
    
    plt.show(block=False)
    plt.pause(.05)
    
    
    
def process_particle_list(particle_list, T, xDim, M):
    P_t = list()
    for t in xrange(T):
        P = np.matrix(particle_list[t*(xDim*M):(t+1)*(xDim*M)])
        P = P.reshape(xDim, M)
        P_t.append(P)
    return P_t

def plot_domain():
    xvec = np.linspace(-5,3,8/.025)
    yvec = np.linspace(-3,3,6/.025)
    imx, imy = np.meshgrid(xvec, yvec)

    sx = np.size(imx,0)
    sy = np.size(imy,1)
    imz = np.ones([sx,sy])

    for i in xrange(sx):
        for j in xrange(sy):
            imz[i,j] = (1.0/((imx[i,j]**2)+1))

    plt.imshow(imz,cmap=matplotlib.cm.Greys_r,extent=[-5,3,-3,3],aspect='equal')