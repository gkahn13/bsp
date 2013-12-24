import numpy as np
from numpy import matlib as ml
import scipy as sci
import cvxopt
import cvxpy

import belief
import math_util
import cvxpy_util
import plot
import util

import IPython


import matplotlib.pyplot as plt



def STOMP_BSP(B,model,plotter,profile):


    #Trajectory will be 3xT

    # set profiler to None if you want no profiling
    profiler = util.Profiler() if profile else None
    cost_obs,B,U = cost_func(B,model,profile,profiler)

    #Intialize variables
    K=9; 
    lmbda = 5e-2
    precost = -100 
    cost = 100
    iteration = 0 
    T = model.T; 

    traj = B[0:model.xDim,:]   
    dt = traj.copy();

    R = compute_R(T); 
    Rin = R.I; 
    M = compute_M(Rin.copy(),T); 

    best_costs = ml.zeros([5,1])+5e6;
    S = ml.zeros([T,5+K]) 
    eps = list()


    for i in range(model.xDim):
         eps.append([ml.zeros([5+K,model.T]),cost_obs]); 

    if profile:
        profiler.start('totalTime')

    #While it hasn't converged or has become invalid 
    while abs(cost - precost) > lmbda and cost < 1000:

       
        traj = B[0:model.xDim,:].copy();  
        precost = cost; 
        
        print 'NEW ITERATION'
        
       
        for k in range(K+5):
            ntraj = traj.copy(); 
            #print "Base Trajecroty",ntraj
            if k > 5:
                #Generate New Trajectories
                for i in range(model.xDim):

                   
                    e =  np.random.multivariate_normal(np.zeros(15),Rin)#*noise_decay**(iteration))
                    e[0] = 0
                    e[T-1] = 0
                    eps[i][0][k,:] = e.T
                    dt[i,:] = e.T    
                  
                
                B[0:model.xDim,:] = ntraj+dt; 

                cost_obs,B,U = cost_func(B,model,profile,profiler)
              
            
            else: 
                #Bring out the five previous best ones
                for i in range(model.xDim):
 
                    dt[i,:] = eps[i][0][k,:]
                cost_obs = eps[0][1]
            
            #Update five best ones
            eps,best_costs = top_trajs(sum(cost_obs),dt,best_costs,eps,model,cost_obs)

            S[:,k] = cost_obs

            #plot.plot_belief_trajectory(B, U, model)
     

        P = compute_probability(S,K+5,model);
         
        #Compute the expected gradient using probablities P
        for x in range(model.xDim):
            for i in range(T):
                dt[x,i] = P[i,:]*eps[x][0][:,i]

            dt[x,:] = (M*dt[x,:].T).T; 
            
      
        traj = traj+dt
       

        B[0:model.xDim,:] = traj


      
       
        iteration+=1
    #Compute Cost and see if convergence by a factor lambda
        cost_obs,B,U = cost_func(B,model,profile,profiler);
        cost = sum(cost_obs); 
     

        if plotter:
            plot.plot_belief_trajectory(B, U, model)
        #Check if valid   
        if(cost < 1000):
            Bopt = B.copy(); 
            Uopt = U.copy();

        #Update five best ones
        eps,best_costs = top_trajs(cost,dt,best_costs,eps,model,cost_obs);


    if profile:
        profiler.stop('totalTime')
        for name, time in profiler.allTimes():
            print('{0}: {1} seconds'.format(name, time))

    return Bopt,Uopt


#Compute Smoothing Matrix
def compute_M(Rinv,T):
    Tfl = float(T)
    M = Rinv

    for i in range(T):
        mx = np.max(Rinv[:,i]); 

        scale = (1/Tfl)/mx;
        M[:,i]= Rinv[:,i]*scale;

    
    M[0,:] = 0; 
    M[T-1,:] = 0; 
    return M; 


#Store five best trajectories
def top_trajs(cost,traj,best_costs,eps,model,cost_obs):
   

    if np.max(best_costs) > cost:

        index = np.argmax(best_costs); 

        for i in range(model.xDim):
            eps[i][0][index,:] = traj[i,:]; 
            eps[i][1] = cost_obs;
        best_costs[index] = cost;
        
    return eps,best_costs


#Compute cost by updating belief with ekf along the trajectory
def cost_func(B,model,profile,profiler):

    cost = ml.zeros([model.T,1])
    U = ml.zeros([model.uDim,model.T])
    T = model.T

    for t in range(T-1):
        U[:,t] = B[0:model.xDim,t+1]-B[0:model.xDim,t];

        B[:,t+1] = belief.belief_dynamics(B[:,t],U[:,t],None,model);
 
        if max(U[:,t])> 1:
            cost[t] = 1000

        elif abs(B[0,t]) > model.xMax[0]:
            cost[t] = 1000

        elif abs(B[1,t]) > model.xMax[1]:
            cost[t] = 1000
            
        else:  
            null, s = belief.decompose_belief(B[:,t], model)
            cost[t] = model.alpha_belief*ml.trace(s*s)+model.alpha_control*ml.sum(U[:,t].T*U[:,t])

    x, s = belief.decompose_belief(B[:,T-1], model)
    cost[T-1] = model.alpha_final_belief*ml.trace(s*s)



    return cost, B, U



#Compute probability of each random trajectory
def compute_probability(S,K,model):
    
    T = model.T; 
    h = 10; 

    Smin = np.zeros(T)
    Smax = np.zeros(T)
    P = ml.zeros([T,K])
    for i in range(T):
        Smin[i] = np.min(S[i,:]); 
        Smax[i] = np.max(S[i,:]); 
        #IPython.embed()
        if(Smin[i] == Smax[i]):
            Smax[i] = 1; 

    
    for i in range(T):
        P[i,:] = np.exp(-h*(S[i,:] - Smin[i])/(Smax[i]-Smin[i]));
      
        normalize = ml.sum(P[i,:]); 
        # if(normalize == K):
        #     P[i,:] = 0.005; 

        P[i,:] = P[i,:]/normalize;
    
     
    return P

#Compute inverse of covariance matrix for random trajectories 
def compute_R(T):

    A =  ml.identity(T)* -2; 
    A[0,0] = 1; 
    A[T-1,T-1] = 1; 

    for i in range(1,T-1):
        A[i,i+1] = 1; 
        A[i,i-1] = 1; 


    R = A.T*A; 
    R = R*200
    return R