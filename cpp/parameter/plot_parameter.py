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

def plot_robot(jointAngles, length1, length2, T):
    plt.figure(1);
    plt.axis([-length1-length2, length1+length2,-length1-length2, length1+length2])
    
    jointAngles = np.array(jointAngles)
    jointAngles = jointAngles.reshape(2,T)
    
    for i in xrange(T):
        j0 = jointAngles[0,i]
        j1 = jointAngles[1,i]
        
        middleLinkPos = [length1*math.cos(j0), length1*math.sin(j0)]
        endEffectorPos = [length1*math.cos(j0) + length2*math.cos(j1), length1*math.sin(j0) + length2*math.sin(j1)]
        
        plt.plot([0, middleLinkPos[0], endEffectorPos[0]], [0, middleLinkPos[1], endEffectorPos[1]])
        
        plt.show(block=False)
        plt.pause(.05)
        raw_input()
        

def plot_parameter_trajectory(B, U, bDim, xDim, uDim, T):
    
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
    
    joint1poscov = [S[0,0] for S in SigmaList]
    joint2poscov = [S[1,1] for S in SigmaList]
    joint1velcov = [S[2,2] for S in SigmaList]
    joint2velcov = [S[3,3] for S in SigmaList]
    length1cov = [S[4,4] for S in SigmaList]
    length2cov = [S[5,5] for S in SigmaList]
    mass1cov = [S[6,6] for S in SigmaList]
    mass2cov = [S[7,7] for S in SigmaList]
        
    #IPython.embed()
        
    plt.figure(1)
    plt.clf()
    plt.cla()
    plt.autoscale(True,'both',True)
    
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
    
    
    plt.figure(2)
    plt.clf()
    plt.cla()
    plt.autoscale(True,'both',True)
    
    plt.subplot(8,1,1)
    plt.ylabel('joint1poscov')
    plt.plot(joint1poscov,'r-')
    
    plt.subplot(8,1,2)
    plt.ylabel('joint2poscov')
    plt.plot(joint2poscov,'r-')
    
    plt.subplot(8,1,3)
    plt.ylabel('joint1velcov')
    plt.plot(joint1velcov,'r-')
    
    plt.subplot(8,1,4)
    plt.ylabel('joint2velcov')
    plt.plot(joint2velcov,'r-')
    
    plt.subplot(8,1,5)
    plt.ylabel('length1cov')
    plt.plot(length1cov,'r-')
    
    plt.subplot(8,1,6)
    plt.ylabel('length2cov')
    plt.plot(length2cov,'r-')
    
    plt.subplot(8,1,7)
    plt.ylabel('mass1cov')
    plt.plot(mass1cov,'r-')
    
    plt.subplot(8,1,8)
    plt.ylabel('mass2cov')
    plt.plot(mass2cov,'r-')
    
    plt.show(block=False)
    plt.pause(.05)
    print 'Press enter:'
    raw_input()

