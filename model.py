import numpy as np

class Model:
    """
    Intended to be inherited for specific problems

    Below is an example of required variables
    and methods
    """
    def __init__(self):
        # model dimensions
        self.xDim = 0. # state space dimension
        self.uDim = 0. # control input dimension
        self.qDim = 0. # dynamics noise dimension
        self.zDim = 0. # observation dimension
        self.rDim = 0. # observtion noise dimension

        # belief space dimension
        # note that we only store the lower (or upper) triangular portion
        # of the covariance matrix to eliminate redundancy
        self.bDim = self.xDim + self.xDim*((self.xDim+1)/2.)

        self.dT = 0. # time step for dynamics function
        self.T = 0. # number of time steps in trajectory

        self.alpha_belief = 0. # weighting factor for penalizing uncertainty at intermediate time steps
        self.alpha_final_belief = 0. # weighting factor for penalizing uncertainty at final time step
        self.alpha_control = 0. # weighting factor for penalizing control cost
        
        self.xMin = np.zeros([self.xDim,1]) # minimum limits on state (xMin <= x)
        self.xMax = np.zeros([self.xDim,1]) # maximum limits on state (x <= xMax)
        self.uMin = np.zeros([self.uDim,1]) # minimum limits on control (uMin <= u)
        self.uMax = np.zeros([self.uDim,1]) # maximum limits on control (u <= uMax)

        self.Q = np.zeros([self.qDim,self.qDim]) # dynamics noise variance
        self.R = np.zeros([self.rDim,self.rDim]) # observation noise variance

        self.start = np.zeros([self.xDim,1]) # start state, OVERRIDE
        self.goal = np.zeros([self.xDim,1]) # end state, OVERRIDE

        self.sqpParams = SqpParams()

    def dynamics_func(self, x_t, u_t, q_t):
        return x_t + self.T*u_t + q_t

    def obs_func(self, x_t, r_t):
        return x_t + r_t

    def plot_domain(self, B):
        pass

    def setStartState(self, start):
        self.start = start

    def setGoalState(self, goal):
        self.goal = goal


class SqpParams:
    """
    Each Model has an instance of SqpParams
    Intended to be inherited for specific problems

    Below is an example of required variables
    """
    def __init__(self):
        self.improve_ratio_threshold = 0.
        self.min_trust_box_size = 0.
        self.min_approx_improve = 0.
        self.max_iter = 0.
        self.trust_shrink_ratio = 0.
        self.trust_expand_ratio = 0.
        self.cnt_tolerance = 0.
        self.max_penalty_coeff_increases = 0.
        self.penalty_coeff_increase_ratio = 0.
        self.initial_trust_box_size = 0.
        self.initial_penalty_coeff = 0.
