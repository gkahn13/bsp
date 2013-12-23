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
    T = model.T;
    # set profiler to None if you want no profiling
    profiler = util.Profiler() if profile else None
    total_cost,B,U = cost_func(B,model)
    total_cost = sum(total_cost)
    #Intialize variables
    K=9; 
    lmbda = 5e-2
    precost = -100 
    cost = total_cost; 
    traj = B[0:model.xDim,:]
    R = compute_R(T); 
    Rin = R.I; 
    M = compute_M(Rin.copy(),T); 
    noise_decay = 7/8; 
    nd = 0; 
    dt = traj.copy();
    T = model.T; 
    best_costs = ml.zeros([5,1])+5e6;
    S = ml.zeros([T,5+K]) 
    eps = list()
    iteration = 0 
    for i in range(model.xDim):
         eps.append(ml.zeros([5+K,model.T])); 
    if profile:
        profiler.start('totalTime')
    while abs(cost - precost) > lmbda and cost != 0:

       
        traj = B[0:model.xDim,:].copy();  
        precost = cost; 
        
        print 'NEW ITERATION'
        
        for k in range(K+5):
            ntraj = traj.copy(); 
            #print "Base Trajecroty",ntraj
            if k > 5:
                for i in range(model.xDim):


                    e =  np.random.multivariate_normal(np.zeros(15),Rin)#*noise_decay**(iteration))
                    e[0] = 0
                    e[T-1] = 0
                    eps[i][k,:] = e.T
                    dt[i,:] = e.T    
                
                print dt
            
            else: 
                for i in range(model.xDim):
 
                    dt[i,:] = eps[i][k,:];  
               
            
            B[0:model.xDim,:] = ntraj+dt; 
            
            
            cost_obs,B,U = cost_func(B,model);
            

            eps,best_costs = top_trajs(sum(cost_obs),dt,best_costs,eps,model)

            S[:,k] = cost_obs;

           
            #plot.plot_belief_trajectory(B, U, model)
     

        P = compute_probability(S,K+5,model);
         
        for x in range(model.xDim):
            for i in range(T):
                dt[x,i] = P[i,:]*eps[x][:,i]

            dt[x,:] = (M*dt[x,:].T).T; 
            
      
        traj = traj+dt
       
   
        #Should smooth here 
        # gcost = cost
        # gprecost = cost+10
        # # Update Trajectory
        # while gprecost - gcost > 1e-2:
        #     traj = traj+3*dt
        #     gprecost = gcost 
        B[0:model.xDim,:] = traj
        #     cost_obs,B,U = cost_func(B,model)
    
        #     gcost = sum(cost_obs)
        #     print gprecost - gcost
            #IPython.embed()

        # if gcost > 1000:
        #     B[0:model.xDim,:] = B[0:model.xDim,:]-6*dt;

      
       
        iteration+=1
    #Compute Cost and see if convergence by a factor lambda
        cost_obs,B,U = cost_func(B,model);
        cost = sum(cost_obs); 
        print "COST LIST",cost_obs
        print "TOTAL COST",cost
 #       IPython.embed()
  #      plt.cla()
   #     plt.clf()
        if plotter:
            plot.plot_belief_trajectory(B, U, model)

        eps,best_costs = top_trajs(cost,dt,best_costs,eps,model);


    if profile:
        profiler.stop('totalTime')
        for name, time in profiler.allTimes():
            print('{0}: {1} seconds'.format(name, time))

    return B,U









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



def top_trajs(cost,traj,best_costs,eps,model):
   

    if np.max(best_costs) > cost:

        index = np.argmax(best_costs); 

        for i in range(model.xDim):
            eps[i][index,:] = traj[i,:]; 
        
        best_costs[index] = cost;
        
    return eps,best_costs


def cost_func(B,model):

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


def compute_R(T):

    A =  ml.identity(T)* -2; 
    A[0,0] = 1; 
    A[T-1,T-1] = 1; 

    for i in range(1,T-1):
        A[i,i+1] = 1; 
        A[i,i-1] = 1; 


    R = A.T*A; 
    R = R*100
    return R