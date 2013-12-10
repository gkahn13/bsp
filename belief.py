import numpy as np
from numpy import matlib as ml
import cvxpy

from math_util import numerical_jac

import IPython

def ekf(x_t, Sigma_t, u_t, z_tp1, model):
# Extract function handles and useful definitions from model
    dynamics_func = model.dynamics_func
    obs_func = model.obs_func
    xDim = model.xDim
    qDim = model.qDim
    rDim = model.rDim
    Q = model.Q
    R = model.R


    dyn_varargin = [x_t, u_t, ml.zeros([qDim,1])]
    obs_varargin = [dynamics_func(x_t, u_t, ml.zeros([qDim,1])), ml.zeros([rDim,1])]

    # dynamics state jacobian
    A = numerical_jac(dynamics_func, 0, dyn_varargin)

    # dynamics noise jacobian
    M = numerical_jac(dynamics_func, 2, dyn_varargin)

    # observation state jacobian
    H = numerical_jac(obs_func, 0, obs_varargin)
    
    # observation noise jacobian
    N = numerical_jac(obs_func, 1, obs_varargin)

    Sigma_tp1_neg = A*Sigma_t*A.T + M*Q*M.T
    K = Sigma_tp1_neg*H.T*ml.linalg.inv(H*Sigma_tp1_neg*H.T + N*R*N.T)

    x = dynamics_func(x_t, u_t, ml.zeros([qDim,1]))
    x_tp1 = x + K*(z_tp1 - obs_func(x, ml.zeros([rDim,1])))
    Sigma_tp1 = (ml.eye(xDim) - K*H)*Sigma_tp1_neg

    return x_tp1, Sigma_tp1

# Belief dynamics: Given belief and control at time t, compute the belief
# at time (t+1) using EKF
def belief_dynamics(b_t, u_t, z_tp1, model):
    dynamics_func = model.dynamics_func
    obs_func = model.obs_func
    qDim = model.qDim
    rDim = model.rDim

    x_t, SqrtSigma_t = decompose_belief(b_t, model)
    Sigma_t = SqrtSigma_t*SqrtSigma_t

    if z_tp1 is None:
        # Maximum likelihood observation assumption
        z_tp1 = obs_func(dynamics_func(x_t, u_t, ml.zeros([qDim,1])), ml.zeros([rDim,1]))
        
    x_tp1, Sigma_tp1 = ekf(x_t, Sigma_t, u_t, z_tp1, model)
    
    # Compute square root for storage
    # Several different choices available -- we use the principal square root
    D_diag, V = ml.linalg.eig(Sigma_tp1)
    D = ml.diag(D_diag)
    SqrtSigma_tp1 = V*np.sqrt(D)*V.T

    b_tp1 = compose_belief(x_tp1, SqrtSigma_tp1, model)

    return b_tp1

# Compose belief vector as [state; square root of covariance matrix]
# Since the principal square root of the covariance matrix is symmetric, we only store the 
# lower (or upper) triangular portion of the square root to eliminate redundancy 
def compose_belief(x, SqrtSigma, model):
    xDim = model.xDim
    bDim = model.bDim

    b = ml.zeros([bDim,1])
    b[0:xDim] = x
    idx = xDim
    for j in xrange(0,xDim):
        for i in xrange(j,xDim):
            b[idx] = 0.5 * (SqrtSigma.item(i,j) + SqrtSigma.item(j,i))
            idx = idx+1

    return b

# Decompose belief vector into mean and square root of covariance
def decompose_belief(b, model):

    xDim = model.xDim

    x = b[0:xDim]
    idx = xDim
    
    # Matlab code
    # if (isa(b, 'double'))
    # SqrtSigma = zeros(xDim, xDim);
    # end
    # unsure what behavior is wanted
    SqrtSigma = ml.zeros([xDim, xDim])

    for j in xrange(0,xDim):
        for i in xrange(j,xDim):
            SqrtSigma[i,j] = b[idx]
            SqrtSigma[j,i] = SqrtSigma[i,j]
            idx = idx+1

    return x, SqrtSigma

# Decompose belief vector into mean and square root of covariance
def cvxpy_decompose_belief1(b, model):

    xDim = model.xDim

    x = b[0:xDim,0]
    idx = xDim
    
    # can create as an Expression??
    SqrtSigma = cvxpy.Variable(xDim, xDim)
    constraints = list()

    for j in xrange(0,xDim):
        for i in xrange(j,xDim):
            constraints.append(SqrtSigma[i,j] == b[idx,0])
            constraints.append(SqrtSigma[j,i] == SqrtSigma[i,j])
            idx = idx+1

    return x, SqrtSigma, constraints

# Decompose belief vector into mean and square root of covariance
def cvxpy_decompose_belief(b, model):

    xDim = model.xDim

    xVec, bVec = b[:xDim,0], b[xDim:,0]

    x = cvxpy.Variable(xDim, 1)
    SqrtSigma = cvxpy.Variable(xDim, xDim)
    constraints = list()

    constraints.append(x == xVec)

    bVec_offset = 0
    for i in xrange(xDim):
        constraints.append( SqrtSigma[i,i:xDim] == bVec[bVec_offset:bVec_offset+(xDim-i)].T )
        bVec_offset += xDim-i

    constraints.append(SqrtSigma == SqrtSigma.T)

    return x, SqrtSigma, constraints


def cvxpy_sigma_trace(b, model):
    xDim = model.xDim

    SqrtSigma = np.zeros([xDim,xDim]).tolist()

    idx = xDim
    for j in xrange(0,xDim):
        for i in xrange(j,xDim):
            SqrtSigma[i][j] = b[idx]
            SqrtSigma[j][i] = SqrtSigma[i][j]
            idx = idx+1

    trace = 0
    for i in xrange(0,xDim):
        trace += sum(cvxpy.square(SqrtSigma[i][j]) for j in xrange(0,xDim))

    return trace
    
            

# converts cvxpy matrix to vector
def cvxpy_vectorize1(A):
    rows, cols = A.size

    Avec = A[:,0]

    for col in xrange(1,cols):
        Avec = cvxpy.vstack(Avec,A[:,col])

    return Avec

# converts cvxpy matrix to vector
def cvxpy_vectorize(A):
    rows, cols = A.size
    constraints = list()

    Avec = cvxpy.Variable(rows*cols, 1)

    for col in xrange(0,cols):
         constraints.append( Avec[rows*col:rows*(col+1)] == A[:,col] )

    return Avec, constraints


if __name__ == '__main__':
    import model
    model = model.Model()
    model.xDim = 4
    x = np.matrix(np.random.rand(4,1))
    svec = np.matrix([0,1,2,3,4,5,6,7,8,9]).T
    b = np.vstack((x,svec))
    x, s, constraints = cvxpy_decompose_belief(b,model)

    x_actual, s_actual = decompose_belief(b,model)

    IPython.embed()
