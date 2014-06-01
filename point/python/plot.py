import numpy as np
from numpy import matlib as ml
import math
import time

import matplotlib.pyplot as plt

import belief

import IPython

def plot_belief_trajectory(B, U, model, start=None, goal=None):
    plt.clf()
    plt.cla()

    xDim = model.xDim
    bDim = model.bDim
    uDim = model.uDim
    start = model.start if start is None else start
    goal = model.goal if goal is None else goal
    T = model.T

    plt.axis([-5,3,-3,3])

    model.plot_domain(B)

    if start is not None:
        plt.plot(start[0,0],start[1,0],'go',markersize=20.0)
    if goal is not None:
        plt.plot(goal[0,0],goal[1,0],'go',markersize=20.0)

    plot_mean(B[0:2,:])

    Xt = ml.zeros([xDim,T])

    for t in xrange(0,T-1):
        Xt[:,t], SqrtSigma_t = belief.decompose_belief(B[:,t], model)
        Sigma_t = SqrtSigma_t*SqrtSigma_t

        plot_cov(Xt[0:2,t], Sigma_t[0:2,0:2])

    Xt[:,T-1], SqrtSigma_T = belief.decompose_belief(B[:,T-1], model)
    Sigma_T = SqrtSigma_T*SqrtSigma_T

    plot_cov(Xt[0:2,T-1], Sigma_T[0:2,0:2])

    plt.show(block=False)
    plt.pause(.05)
    
def plot_belief_trajectory_cpp(B, U, model, start, goal, T):
    bDim = model.bDim
    uDim = model.uDim
    
    B = np.matrix(B)
    B = np.matrix(B.reshape(bDim, T))
    
    model.T = T
    model.setStartState(ml.matrix(start).T)
    model.setGoalState(ml.matrix(goal).T)

    plot_belief_trajectory(B, U, model)

    
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
