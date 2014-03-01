import numpy as np
from numpy import matlib as ml
import math
import time

import matplotlib
import matplotlib.pyplot as plt
from scipy.interpolate import spline


import IPython



# belief is just belief of point robot, not waypoints
# T is time for total trajectory (i.e. to all the waypoints)
def plot_point_trajectory(B, U, waypoints, landmarks, max_range, alpha_obs, xDim, T, DT):
    plt.clf()
    plt.cla()

    sDim = (xDim*(xDim+1))/2.
    bDim = xDim + sDim
    uDim = 2
    cDim = 3
    
    B = np.matrix(B)
    B = B.reshape(bDim, T)
    
    U = np.matrix(U)
    U = U.reshape(uDim, T-1)
    
    numWaypoints = len(waypoints)/2
    waypoints = np.matrix(waypoints)
    waypoints = waypoints.reshape(2, numWaypoints)
    
    numLandmarks = len(landmarks)/2
    landmarks = np.array(landmarks)
    landmarks = landmarks.reshape(2, numLandmarks)
    
    extent = [-10,70,-25,40]
    plt.axis(extent)
    #plt.axis('equal')
    
    #smooth_trajectory(B, U, cDim, uDim, T, DT)
    
    #plot_domain(landmarks, max_range, alpha_obs, extent)
    plot_landmarks(B, landmarks, max_range, bDim, xDim, T)
    
    # plot mean of trajectory
    plot_mean(B[0:3,:], U, DT)
    
    plt.plot(landmarks[0,:], landmarks[1,:], ls='None', color='red', marker='x', markersize = 5.0)
    plt.plot(waypoints[0,:], waypoints[1,:], ls='None', color='purple', marker='s', markersize=8.0)

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
    
    #raw_input()
    
def plot_landmarks(B, landmarks, max_range, bDim, xDim, T):
    for i in xrange(landmarks.shape[1]):
        pos = landmarks[:,i]
        plot_cov(pos, max_range*max_range*ml.identity(2), 'b-')
        for t in xrange(T):
            x, SqrtSigma = decompose_belief(B[:,t], bDim, xDim)
            Sigma = SqrtSigma*SqrtSigma
            
            plot_cov(pos, Sigma[3+2*i:3+2*i+2,3+2*i:3+2*i+2], 'r-', alpha=(t+1)/float(T))
        
    
def plot_domain(landmarks, max_range, alpha_obs, extent):
    # note x and y are flipped for plotting!
    granularity = .25
    
    numLandmarks = landmarks.shape[1]
    
    xvec = np.linspace(extent[0],extent[1],(extent[1]-extent[0])/granularity)
    yvec = np.linspace(extent[2],extent[3],(extent[3]-extent[2])/granularity)
    imx, imy = np.meshgrid(xvec, yvec)

    sx = np.size(imx,0)
    sy = np.size(imy,1)

    imz = np.zeros([sx,sy])

    for i in xrange(sx):
        for j in xrange(sy):
            for k in xrange(numLandmarks):
                lx = landmarks[1,k]
                ly = landmarks[0,k]
                dist = math.sqrt((lx-imy[i,j])**2 + ((ly-imx[i,j])**2))
                try:
                    imz[i,j] += (1.0)/(1+math.exp(-alpha_obs*(max_range-dist)))
                except:
                    pass
                    
    imz = imz / imz.max()
    plt.imshow(imz,cmap=matplotlib.cm.Greys_r,extent=extent,aspect='equal',origin='lower')
    
def smooth_trajectory(B, U, cDim, uDim, T, DT):
    smooth_factor = 10
    dt = DT / float(smooth_factor)
    
    X = B[0:cDim,0]
    xidx = 0
    x = X[:,0]
    for t in xrange(T-1):
        u = U[:,t]
        x = B[0:cDim,t]
        X = np.hstack((X,x))
        for i in xrange(1,smooth_factor+1):
            x = dynfunc(x, u, dt)
            X = np.hstack((X,x))
            
    X = np.asarray(X)
    plt.plot(X[0,:],X[1,:])
    plt.show(block=False)
    plt.pause(.5)
    IPython.embed()
    
 
def dynfunc(x, u, dt, wheelbase=4.):
    xAdd = np.matrix(x)
    
    xAdd[0,0] = u[0,0] * dt * math.cos(x[2,0] + u[1,0]*dt)
    xAdd[1,0] = u[0,0] * dt * math.sin(x[2,0] + u[1,0]*dt)
    xAdd[2,0] = u[0,0] * dt * math.sin(u[1,0]*dt) / wheelbase
    
    return x + xAdd 
        
    
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

    
def plot_mean(X, U, DT, interp=False):
    
    X = np.asarray(X)
    
    if interp:
        xdata = X[0,:]
        ydata = X[1,:]
        
        xnew = np.linspace(xdata.min(), xdata.max(), 300)
        
        ysmooth = spline(xdata, ydata, xnew)
        
        plt.plot(xnew, ysmooth, color='red')
        plt.plot(X[0,:],X[1,:],ls='None',marker='s',markerfacecolor='yellow')
    else:
        plt.plot(X[0,:],X[1,:],color='red',marker='s',markerfacecolor='yellow')
    
    
    plt.plot(X[0,0],X[1,0],ls='None',marker='s',markersize=10.0)
    plt.plot(X[0,-1],X[1,-1],ls='None',marker='s',markersize=10.0)
    
    """
    X = np.asmatrix(X)
    
    T = X.shape[1]
    sampling = 2
    X_upsampled = ml.zeros([X.shape[0], T*sampling])
    X_upsampled[:,0] = X[:,0]
    index = 0
    for i in xrange(0, T-1):
        for j in xrange(0, sampling):
            X_upsampled[:,index+1] = dynfunc(X_upsampled[:,index], U[:,i], DT/float(sampling))
            index += 1
    
    X = np.asarray(X)
    X_upsampled = np.asarray(X_upsampled)
    plt.plot(X_upsampled[0,:],X_upsampled[1,:],color='red')
    plt.plot(X[0,:],X[1,:],ls='None',color='red',marker='s',markerfacecolor='yellow')
    plt.plot(X[0,0],X[1,0],ls='None',marker='s',markersize=10.0)
    plt.plot(X[0,-1],X[1,-1],ls='None',marker='s',markersize=10.0)
    """
    

def plot_cov(mu, sigma, plotType = 'y-', alpha=1):
    mu = np.asarray(mu)
    sigma = np.asarray(sigma)

    t = np.linspace(-math.pi, math.pi, 2*math.pi/.1)
    x = np.sin(t)
    y = np.cos(t)

    D_diag, V = ml.linalg.eigh(sigma)
    D = np.diag(D_diag)
    A = np.real((np.dot(V,ml.sqrt(D))).T)

    z = np.dot(np.vstack((x.T,y.T)).T,A)

    plt.plot(z[:,0]+mu[0], z[:,1]+mu[1], plotType, alpha=alpha)
