import numpy as np
from numpy import matlib as ml
import math
import time

import matplotlib.pyplot as plt

import belief

import IPython

# Compose belief vector as [state; square root of covariance matrix]
# Since the principal square root of the covariance matrix is symmetric, we only store the 
# lower (or upper) triangular portion of the square root to eliminate redundancy 
def compose_belief(x, SqrtSigma, bDim, xDim):
    b = ml.zeros([bDim,1])
    b[0:xDim] = x
    idx = xDim
    for j in xrange(0,xDim):
        for i in xrange(j,xDim):
            b[idx] = 0.5 * (SqrtSigma.item(i,j) + SqrtSigma.item(j,i))
            idx = idx+1

    return b

# Decompose belief vector into mean and square root of covariance
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

def plot_belief_trajectory(B, bDim, xDim, T):
    plt.clf()
    plt.cla()
    
    B = np.array(B)
    B.reshape(bDim, T)
    
    numLandmarks = float(xDim - 3) / 2.

    #plt.plot(start[0],start[1],'go',markersize=20.0)
    #plt.plot(goal[0],goal[1],'go',markersize=20.0)

    Xt = ml.zeros([xDim,T])

    for t in xrange(0,T):
        X[:,t], SqrtSigma_t = decompose_belief(B[1:,t], bDim, xDim)
        
        # plot car
        X_car = X[0:2, t]
        SqrtSigma_t_car = SqrtSigma_t[0:2, 0:2]
        
        plt.plot(X_car[0,0],X_car[1,0],color='red',marker='s',markerfacecolor='yellow')
        
        Sigma_t_car = SqrtSigma_t_car*SqrtSigma_t_car

        plot_cov(X_car, Sigma_t_car)
        
        # plot landmarks
        for i in xrange(0, numLandmarks):
            X_landmark = X[3+i*2:3+i*2+2, t]
            SqrtSigma_t_landmark = SqrtSigma_t[3+i*2:3+i*2+2, 3+i*2:3+i*2+2]
            
            plt.plot(X_landmark[0,0],X_landmark[1,0],color='blue',marker='o')
        
            Sigma_t_landmark = SqrtSigma_t_landmark*SqrtSigma_t_landmark
            
            plot_cov(X_landmark, Sigma_t_landmark)
            

        plt.show(block=False)
        plt.pause(.05)
    

    
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
