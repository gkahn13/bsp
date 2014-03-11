import numpy as np
from numpy import matlib as ml
import math
import time
import pickle 

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
    Base = X.copy(); 

    SigmaList = list()
    for i in xrange(T):
        Base[0,:] = 0.5; 
        X[:,i], Sigma = decompose_belief(B[:,i], bDim, xDim)
        SigmaList.append(Sigma)
   
    output =open('data.pkl','wb')

    pickle.dump([X, SigmaList] ,output)
    joint1pos = X[0,:].tolist()[0]
    joint2pos = X[1,:].tolist()[0]
    joint1vel = X[2,:].tolist()[0]
    joint2vel = X[3,:].tolist()[0]
    length1 = (1/X[4,:]).tolist()[0]
    length2 = (1/X[5,:]).tolist()[0]
    mass1 = (1/X[6,:]).tolist()[0]
    mass2 = (1/X[7,:]).tolist()[0]
    base = (Base[0,:]).tolist()[0]

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
   
    plt.autoscale(True,'x',True)
    
 
    
    plt.subplot(6,1,1)
    plt.ylim((0.2,0.7))
    plt.ylabel('length1')
    plt.plot(length1,'b-')
    plt.plot(base,'r-')
    
    plt.subplot(6,1,2)
    plt.ylabel('length2')
    plt.plot(length2,'b-')
    plt.plot(base,'r-')

    plt.subplot(6,1,3)
    plt.ylim((0.2,0.7))
    plt.ylabel('mass1')

    plt.plot(mass1,'b-')
    plt.plot(base,'r-')

    plt.subplot(6,1,4)
    plt.ylim((0.2,0.7))
    plt.ylabel('mass2')
    plt.plot(mass2,'b-')
    plt.plot(base,'r-')

    plt.subplot(6,1,5)
    plt.ylabel('j1')
    plt.plot(joint1pos,'b-')

    plt.subplot(6,1,6)
    plt.ylabel('j2')
    plt.plot(joint2pos,'b-')

    
    plt.figure(2)
    plt.clf()
    plt.cla()
    plt.autoscale(True,'x',True)
    
    
    plt.subplot(4,1,1)
    plt.ylabel('length1cov')
    plt.ylim((0.0,1.0))
    plt.plot(length1cov,'r-')
    
    plt.subplot(4,1,2)
    plt.ylabel('length2cov')
    plt.ylim((0.0,1.0))
    plt.plot(length2cov,'r-')
    
    plt.subplot(4,1,3)
    plt.ylabel('mass1cov')
    plt.ylim((0.0,1.0))
    plt.plot(mass1cov,'r-')
    
    plt.subplot(4,1,4)
    plt.ylabel('mass2cov')
    plt.ylim((0.0,1.0))
    plt.plot(mass2cov,'r-')
    
    plt.show(block=False)
    plt.pause(.05)
    print 'Press enter:'
    raw_input()

def plot_speed_bar(): 

    groups = 3; 
    index = np.arange(groups)

    Height_1l = (3.638762,151.98886,1594.754)
    std_1l = (2.20,109.3,67)

    Height_ml = (418.5,16901.0,206438.4)
    std_ml = (26.1,3897.4,36624.0)
    error_config = {'ecolor': '0.3'}
    bar_width = 0.4
    plt.figure(1)
    plt.clf()
    
    plt.cla()
    plt.bar(index,Height_1l,bar_width,yerr=std_1l,error_kw=error_config)
    plt.xticks(index+bar_width/2,('controls','state','belief'))
    plt.title('Average Speed Over One Step')
    plt.ylabel('Time (ms)')

    plt.figure(2)
    plt.clf()
    plt.cla()
    plt.bar(index,Height_ml,bar_width,yerr=std_ml,error_kw=error_config)
    plt.xticks(index+bar_width/2,('controls','state','belief'))
    plt.ylabel('Time (ms)')
    plt.title('Average Speed Over 300 Step')
    plt.show(block=False)
    plt.pause(.05)
    print 'Press enter:'



def plot_for_paper():
    
   
    bDim = 44
    uDim = 2
    T = 500
    xDim = 8

    pkl_file = open('data.pkl', 'rb')
    pkl_file2 = open('random.pkl', 'rb')

    [X, SigmaList] = pickle.load(pkl_file)
    [Xrand, SigmaListR] = pickle.load(pkl_file2)
    Base = X; 
   # IPython.embed()
    Base[0,T-1] = 0.5; 
    
         
        
  
    length1_rand = (1/Xrand[4,:]).tolist()[0]
    length1cov_rand = [S[4,4] for S in SigmaListR]

    length2_rand = (1/Xrand[5,:]).tolist()[0]
    length2cov_rand = [(S[5,5]) for S in SigmaListR]

    mass1_rand = (1/Xrand[6,:]).tolist()[0]
    mass1cov_rand = [(S[6,6]) for S in SigmaListR]

    mass2_rand = (1/Xrand[7,:]).tolist()[0]
    mass2cov_rand = [(S[7,7]) for S in SigmaListR]

    length1 = (1/X[4,:]).tolist()[0]
    length1cov = [S[4,4] for S in SigmaList]

    length2 = (1/X[5,:]).tolist()[0]
    length2cov = [S[5,5] for S in SigmaList]

    mass1 = (1/X[6,:]).tolist()[0]
    mass1cov = [S[6,6] for S in SigmaList]

    mass2 = (1/X[7,:]).tolist()[0]
    mass2cov = [S[7,7] for S in SigmaList]

    base = [0.5 for i in range(500)]

    Count = [i for i in range(500)]
        
    #IPython.embed()
        
    plt.figure(1)
    plt.clf()
    plt.cla()
   
    plt.autoscale(True,'x',True)
    lightgreen = '#00FF44'
    
    
    plt.subplot(4,1,1)
    plt.ylim((0.2,0.7))
    plt.ylabel('length1')
    plt.errorbar(Count,length1_rand,length1cov_rand,ecolor = 'c',errorevery=1)
    plt.errorbar(Count,length1,length1cov,ecolor=lightgreen,errorevery=1)
    plt.plot(length1_rand,'b',linewidth= 2)
    plt.plot(length1,'g-', linewidth = 2)
    plt.plot(base,'r-', linewidth = 2)
    
    plt.subplot(4,1,2)
    plt.ylim((0.2,0.7))
    plt.ylabel('length2')
    plt.errorbar(Count,length2_rand,length2cov_rand,ecolor ='c',errorevery=1)
    plt.errorbar(Count,length2,length2cov, ecolor=lightgreen,errorevery=1)
    plt.plot(length2_rand,'b',linewidth= 2)
    plt.plot(length2,'g-', linewidth = 2)
    plt.plot(base,'r-', linewidth = 2)


    plt.subplot(4,1,3)
    plt.ylim((0.2,0.7))
    plt.ylabel('mass1')
    plt.errorbar(Count,mass1_rand,mass1cov_rand,ecolor ='c',errorevery=1)
    plt.errorbar(Count,mass1,mass1cov,ecolor=lightgreen,errorevery=1)
    plt.plot(mass1_rand,'b',linewidth= 2)
    plt.plot(mass1,'g-', linewidth = 2)
    plt.plot(base,'r-', linewidth = 2)

    plt.subplot(4,1,4)
    plt.ylim((0.2,0.7))
    plt.ylabel('mass2')
    plt.xlabel('Timesteps')
    plt.errorbar(Count,mass2_rand,mass2cov_rand,ecolor ='c',errorevery=1)
    plt.errorbar(Count,mass2,mass2cov,ecolor=lightgreen,errorevery=1)
    plt.plot(mass2_rand,'b',linewidth= 2)
    plt.plot(mass2,'g-', linewidth = 2)
    plt.plot(base,'r-', linewidth = 2)
    
    
    
   
    plt.show(block=False)
    plt.pause(.05)
    print 'Press enter:'
    raw_input()


