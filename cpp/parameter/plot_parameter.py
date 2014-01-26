import numpy as np
from numpy import matlib as ml
import math
import time

import matplotlib.pyplot as plt

import IPython

def compose_belief(x, Sigma, bDim, xDim):
    b = ml.zeros([bDim,1])
    b[0:xDim] = x
    idx = xDim
    for j in xrange(0,xDim):
        for i in xrange(j,xDim):
            b[idx] = Sigma.item(i,j)
            idx = idx+1

    return b

def decompose_belief(b, bDim, xDim):
    x = b[0:xDim]
    idx = xDim
    
    Sigma = ml.zeros([xDim, xDim])

    for j in xrange(0,xDim):
        for i in xrange(j,xDim):
            Sigma[i,j] = b[idx]
            Sigma[j,i] = Sigma[i,j]
            idx = idx+1

    return x, Sigma

def plot_parameter_trajectory(B, U, bDim, xDim, uDim, T):
    plt.clf()
    plt.cla()
    
    #plt.axes([0,T,-2,2])
    plt.autoscale(True,'both',True)
    
    B = np.matrix(B)
    U = np.matrix(U)
    
    B = B.reshape(bDim, T)
    U = U.reshape(uDim, T-1)
    
    X = ml.zeros([xDim, T])
    SigmaList = list()
    for i in xrange(T):
        X[:,i], Sigma = decompose_belief(B[:,i], bDim, xDim)
        SigmaList.append(Sigma)
        
    joint1pos = X[0,:].tolist()[0]
    joint2pos = X[1,:].tolist()[0]
    joint1vel = X[2,:].tolist()[0]
    joint2vel = X[3,:].tolist()[0]
    length1 = X[4,:].tolist()[0]
    length2 = X[5,:].tolist()[0]
    mass1 = X[6,:].tolist()[0]
    mass2 = X[7,:].tolist()[0]
        
    plt.figure(1)
    plt.subplot(8,1,1)
    plt.ylabel('joint1pos')
    plt.plot(joint1pos,'b-')
    
    plt.subplot(8,1,2)
    plt.ylabel('joint2pos')
    plt.plot(joint2pos,'b-')
    
    plt.subplot(8,1,3)
    plt.ylabel('joint1vel')
    plt.plot(joint1vel,'b-')
    
    plt.subplot(8,1,4)
    plt.ylabel('joint2vel')
    plt.plot(joint2vel,'b-')
    
    plt.subplot(8,1,5)
    plt.ylabel('length1')
    plt.plot(length1,'b-')
    
    plt.subplot(8,1,6)
    plt.ylabel('length2')
    plt.plot(length2,'b-')
    
    plt.subplot(8,1,7)
    plt.ylabel('mass1')
    plt.plot(mass1,'b-')
    
    plt.subplot(8,1,8)
    plt.ylabel('mass2')
    plt.plot(mass2,'b-')
    
    plt.show(block=False)
    plt.pause(.05)
    #print 'Press enter:'
    #raw_input()

