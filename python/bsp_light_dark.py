import numpy as np
from numpy import matlib as ml
import matplotlib
import matplotlib.pyplot as plt

import model
import belief
import belief_opt
import belief_grad_free
import plot

import IPython

class LightDarkModel(model.Model):
    def __init__(self):
        # model dimensions
        self.xDim = 2 # state space dimension
        self.uDim = 2 # control input dimension
        self.qDim = 2 # dynamics noise dimension
        self.zDim = 2 # observation dimension
        self.rDim = 2 # observtion noise dimension

        # belief space dimension
        # note that we only store the lower (or upper) triangular portion
        # of the covariance matrix to eliminate redundancy
        self.bDim = int(self.xDim + self.xDim*((self.xDim+1)/2.))

        self.dT = 1. # time step for dynamics function
        self.T = 15 # number of time steps in trajectory

        self.alpha_belief = 10. # weighting factor for penalizing uncertainty at intermediate time steps
        self.alpha_final_belief = 10. # weighting factor for penalizing uncertainty at final time step
        self.alpha_control = 1. # weighting factor for penalizing control cost
        
        self.xMin = ml.vstack([-5,-3]) # minimum limits on state (xMin <= x)
        self.xMax = ml.vstack([5,3]) # maximum limits on state (x <= xMax)
        self.uMin = ml.vstack([-1,-1]) # minimum limits on control (uMin <= u)
        self.uMax = ml.vstack([1,1]) # maximum limits on control (u <= uMax)

        self.Q = ml.eye(self.qDim) # dynamics noise variance
        self.R = ml.eye(self.rDim) # observation noise variance

        self.start = ml.zeros([self.xDim,1]) # start state, OVERRIDE
        self.goal = ml.zeros([self.xDim,1]) # end state, OVERRIDE

        self.sqpParams = LightDarkSqpParams()

    def dynamics_func(self, x_t, u_t, q_t):
        return x_t + self.dT*u_t + .01*q_t

    def obs_func(self, x_t, r_t):
        intensity = 0.5*0.5*x_t.item(0,0)*x_t.item(0,0) + 1e-6
        return x_t + np.sqrt(intensity)*r_t

    def plot_domain(self, B):
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
        

class LightDarkSqpParams(model.SqpParams):
    def __init__(self):
        self.improve_ratio_threshold = .1
        self.min_trust_box_size = 1e-3
        self.min_approx_improve = 1e-4
        self.max_iter = 50.
        self.trust_shrink_ratio = .1
        self.trust_expand_ratio = 1.5
        self.cnt_tolerance = 1e-4
        self.max_penalty_coeff_increases = 3.
        self.penalty_coeff_increase_ratio = 10.
        self.initial_trust_box_size = 1.
        self.initial_penalty_coeff = 50.

def test_bsp_light_dark():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('--no-plotting',action='store_true',default=False)
    parser.add_argument('--profile',action='store_true',default=False)
    parser.add_argument('--gradient_free',action='store_true',default=False)
    args = parser.parse_args()

    plotting = not args.no_plotting
    profile = args.profile
    gradient_free = args.gradient_free
    model = LightDarkModel()

    X1 = ml.matrix([[-3.5,2],[-3.5,-2],[-4,0],[2,2],[-4,2]]).T
    G =  ml.matrix([[-3.5,-2],[-3.5,2],[-1,0],[2,-2],[-1,-2]]).T

    # [Reference] final belief trajectory costs for verification 
    # Allow for numerical inaccuracy
    verify_cost = [45.181701, 45.181643, 49.430339, 27.687003, 56.720314]

    for i_problem in xrange(0,5):
        # Setup initial conditions for problem
        x1 = X1[:,i_problem] # start (mean of initial belief)
        SqrtSigma1 = ml.eye(2) # initial covariance
        goal = G[:,i_problem] # goal
    
        model.setStartState(x1)
        model.setGoalState(goal)    
    
        # setup initial control vector -- straight line initialization from start to goal
        U = ml.tile(((model.goal - model.start)/(model.T-1)), (1,model.T-1))
    
        B = ml.zeros([model.bDim,model.T])
        B[:,0] = belief.compose_belief(x1, SqrtSigma1, model)
        for t in xrange(0,model.T-1):
            B[:,t+1] = belief.belief_dynamics(B[:,t], U[:,t], None, model,None,None)

        # display initialization
        if plotting:
            plot.plot_belief_trajectory(B, U, model)
    
        if gradient_free :
            [Bopt, Uopt] = belief_grad_free.STOMP_BSP(B,model,plotting,profile)
        else:
            [Bopt, Uopt] = belief_opt.belief_opt_penalty_sqp(B, U, model, plotting, profile)
        if plotting:
            plot.plot_belief_trajectory(Bopt, Uopt, model);
    
	# Forward simulated cost
        cost = belief.compute_forward_simulated_cost(B[:,0], Uopt, model)
        print('Total cost of optimized trajectory: %f' % cost)
    
        # save trajectory to png file 
        # saveas(gcf, sprintf('bsp-light-dark-plan-%i.png',i_problem));
    
        print('For verification (allow for numerical inaccuracy):')
        print('(Reference) Total cost of optimized trajectory: %2.6f' % verify_cost[i_problem])
    
        # Simulate execution of trajectory (test to see why the maximum likelihood observation assumption does not hold)
        #simulate_bsp_trajectory(B(:,1), Uopt, model);
    
        print('press enter to continue to the next problem')
        raw_input()
    

if __name__ == '__main__':
    test_bsp_light_dark()
