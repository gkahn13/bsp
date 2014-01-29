import numpy as np
from numpy import matlib as ml
import math
import time

import matplotlib.pyplot as plt


import IPython



# belief is just belief of point robot, not waypoints
# T is time for total trajectory (i.e. to all the waypoints)
def plot_point_trajectory(B, waypoints, T):
    plt.clf()
    plt.cla()

    xDim = 2
    bDim = 5
    uDim = 2
    
    B = np.matrix(B)
    B = B.reshape(bDim, T)
    
    numWaypoints = len(waypoints)/2
    waypoints = np.matrix(waypoints)
    waypoints = waypoints.reshape(2, numWaypoints)
    
    plt.axis([-10,70,-10,50])
    #IPython.embed()
    #plt.axis('equal')
    
    # plot mean of trajectory
    plot_mean(B[0:2,:])
    
    
    for i in xrange(numWaypoints):
        plt.plot(waypoints[0,i],waypoints[1,i],color='purple',marker='s',markersize=8.0)

    Xt = ml.zeros([xDim,T])

    for t in xrange(0,T-1):
        Xt[:,t], SqrtSigma_t = decompose_belief(B[:,t], bDim, xDim)
        Sigma_t = SqrtSigma_t*SqrtSigma_t

        plot_cov(Xt[0:2,t], Sigma_t[0:2,0:2])

    Xt[:,T-1], SqrtSigma_T = decompose_belief(B[:,T-1], bDim, xDim)
    Sigma_T = SqrtSigma_T*SqrtSigma_T

    plot_cov(Xt[0:2,T-1], Sigma_T[0:2,0:2])
    
    plt.show(block=False)
    plt.pause(.05)
    
    raw_input()
    

def compose_belief(x, SqrtSigma, bDim, xDim):
    b = ml.zeros([bDim,1])
    b[0:xDim] = x
    idx = xDim
    for j in xrange(0,xDim):
        for i in xrange(j,xDim):
            b[idx] = 0.5 * (SqrtSigma.item(i,j) + SqrtSigma.item(j,i))
            idx = idx+1

    return b

def decompose_belief(b, bDim, xDim):
    x = b[0:xDim]
    idx = xDim
    
    SqrtSigma = ml.zeros([xDim, xDim])

    for j in xrange(0,xDim):
        for i in xrange(j,xDim):
            SqrtSigma[i,j] = b[idx]
            SqrtSigma[j,i] = SqrtSigma[i,j]
            idx = idx+1

    return x, SqrtSigma

    
def plot_mean(X):
    X = np.asarray(X)
    plt.plot(X[0,:],X[1,:],color='red',marker='s',markerfacecolor='yellow')
    
    plt.plot(X[0,0],X[1,0],'bs',markersize=10.0)
    plt.plot(X[0,-1],X[1,-1],'bs',markersize=10.0)
    

def plot_cov(mu, sigma):
    mu = np.asarray(mu)
    sigma = np.asarray(sigma)

    t = np.linspace(-math.pi, math.pi, 2*math.pi/.1)
    x = np.sin(t)
    y = np.cos(t)

    D_diag, V = ml.linalg.eigh(sigma)
    D = np.diag(D_diag)
    A = np.real((np.dot(V,ml.sqrt(D))).T)

    z = np.dot(np.vstack((x.T,y.T)).T,A)

    plt.plot(z[:,0]+mu[0], z[:,1]+mu[1], 'y-')
