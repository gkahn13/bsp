import numpy as np
from numpy import matlib as ml
import scipy as sci
import cvxopt
import cvxpy

import belief
import math_util
import cvxpy_util
import plot

import IPython

def belief_opt_penalty_sqp(B, U, model):
    cfg = model.sqpParams
      
    trust_box_size = cfg.initial_trust_box_size # The trust region will be a box around the current iterate.
    penalty_coeff = cfg.initial_penalty_coeff # Coefficient of l1 penalties 

    # TODO: The outer loop of the sqp algorithm, which repeatedly minimizes
    # the merit function --- Calls minimize_merit_function defined below

    # After this call, check to see if the
    # constraints are satisfied.
    # - If some constraint is violated, increase penalty_coeff by a factor of cfg.merit_coeff_increase_ratio
    # You should also reset the trust region size to be larger than cfg.min_trust_box_size,
    # which is used in the termination condition for the inner loop.
    # - If all constraints are satisfied (which in code means if they are satisfied up to tolerance cfg.cnt_tolerance), we're done.

    # cvx_quiet(true);
    # cvx_solver SDPT3;

    iterCount = 1
    while iterCount <= cfg.max_iter:    
    
        [B, U, success] = minimize_merit_function(B, U, model, cfg, penalty_coeff, trust_box_size)

        IPython.embed()
    
        if constraints_satisfied(B, U, model, cfg.cnt_tolerance):
            break
    
        if success is False:
            break
    
        penalty_coeff = cfg.initial_penalty_coeff
        trust_box_size = cfg.initial_trust_box_size

    return B, U

# Given a belief trajectory and a set of controls, compute merit function value -- Eq. 6 in the problem set
def compute_merit(B, U, model, penalty_coeff):
    xDim = model.xDim
    T = model.T
    merit = 0
    
    for t in xrange(0,T-1):
        x, s = belief.decompose_belief(B[:,t], model)
        merit += model.alpha_belief*ml.trace(s*s)
        merit += model.alpha_control*ml.sum(U[:,t].T*U[:,t])

        merit += penalty_coeff*ml.sum(np.abs(B[:,t+1]-g(B[:,t],U[:,t],model)))
    
    x, s = belief.decompose_belief(B[:,T-1], model)
    merit += model.alpha_final_belief*ml.trace(s*s)

    return merit

# Same as compute_merit, but linearized version for use with optimization
def linearized_compute_merit(B, Bcvx, U, Ucvx, gval, G, H, model, penalty_coeff):
    xDim = model.xDim
    uDim = model.uDim
    T = model.T

    trace_merits = list()
    control_merits = list()
    belief_penalty_merits = list()
    
    constraints = list()

    for t in xrange(0,T-1):
        x, s, decompose_constraint = belief.cvxpy_decompose_belief(Bcvx[:,t], model)
        constraints += decompose_constraint
        trace_merits.append(model.alpha_belief*cvxpy_util.sum_square(s))

        control_merits.append(model.alpha_control*cvxpy.quad_over_lin(Ucvx[:,t],1))
        
        belief_penalty_merits.append(penalty_coeff*cvxpy.sum(cvxpy.abs(B[:,t+1] - (gval[t] + G[t]*(Bcvx[:,t]-B[:,t]) + H[t]*(Ucvx[:,t]-U[:,t])))))
    
    x, s, decompose_constraint = belief.cvxpy_decompose_belief(Bcvx[:,T-1], model)
    constraints += decompose_constraint
    trace_merits.append(model.alpha_belief*cvxpy_util.sum_square(s))

    #merit = cvxpy.Variable()
    #constraints.append(merit == sum(trace_merits) + sum(control_merits) + sum(belief_penalty_merits))
    merit = sum(trace_merits) + sum(control_merits) + sum(belief_penalty_merits)

    return merit, constraints

def minimize_merit_function(B, U, model, cfg, penalty_coeff, trust_box_size):
    success = True
    sqp_iter = 1
    
    xDim = model.xDim
    bDim = model.bDim
    uDim = model.uDim
    T = model.T;

    Borig = B
    Uorig = U
    
    while True:
        # In this loop, we repeatedly construct a linear approximation to
        # the nonlinear belief dynamics constraint
        print('  sqp iter: %i' % sqp_iter)
        
        # Compute merit value using the current iterate (trajectory and set of controls)
        merit = compute_merit(B, U, model, penalty_coeff);
        print('merit: %f' % merit)
        
        print('     current merit: %g' % merit);
        
        # LINEARIZE THE BELIEF DYNAMICS CONSTRAINT (Eq. 5b)
	# Use the numerical_jac routine in q1_starter to compute the required Jacobians numerically
        gval, G, H = list(), list(), list()
        for t in xrange(0, T-1):
            gval.append(g(B[:,t], U[:,t], model))
            G.append(math_util.numerical_jac(g, 0, [B[:,t], U[:,t], model]))
            H.append(math_util.numerical_jac(g, 1, [B[:,t], U[:,t], model]))

        #IPython.embed()

        while True:
            # This is the trust region loop
            # Using the approximations computed above, this loop shrinks
            # the trust region until the progress on the approximate merit
            # function is a sufficiently large fraction of the progress on
            # the exact merit function.
            
            print('Press enter to continue with SQP')
            raw_input()

            print('    trust region size: %.3g' % trust_box_size)

            Bcvx = cvxpy.Variable(bDim, T)
            Ucvx = cvxpy.Variable(uDim, T-1)
            constraints = list()

            objective_merit, linearized_constraints = linearized_compute_merit(B,Bcvx,U,Ucvx,gval,G,H,model,penalty_coeff)
            objective = cvxpy.Minimize(objective_merit)
            

            constraints += linearized_constraints
            
            constraints += [Bcvx[:,0] == B[:,0], # Constraint to ensure that initial belief remains unchanged
                           Bcvx[0:xDim,T-1] == B[0:xDim,T-1], # reach goal at time T
                           Bcvx[0:xDim,:] <= np.tile(model.xMax, (1,T)), # upper x bound
                           Bcvx[0:xDim,:] >= np.tile(model.xMin, (1,T)), # lower x bound
                           Ucvx <= np.tile(model.uMax, (1,T-1)), # upper u bound
                           Ucvx >= np.tile(model.uMin, (1,T-1))] # lower u bound
            

            
            for t in xrange(0,T-1):
                constraints.append(cvxpy.norm(B[:,t]-Bcvx[:,t]) <= trust_box_size)
                constraints.append(cvxpy.norm(U[:,t]-Ucvx[:,t]) <= trust_box_size)
                
            constraints.append(cvxpy.norm(B[:,T-1]-Bcvx[:,T-1]) <= trust_box_size)
            

            problem = cvxpy.Problem(objective, constraints)

            try:
                cvx_optval = problem.solve()
            except Exception as e:
                print('Failed to solve QP subproblem.')
                print('Error %s' % e)
                IPython.embed()
                return B, U, False

            Bcvx = np.matrix(Bcvx.value)
            Ucvx = np.matrix(Ucvx.value)

            # TEMP!!!!
            B = Bcvx
            U = Ucvx
            plot.plot_belief_trajectory(B,U,model)

            model_merit = cvx_optval
			
            # Compute merit value using the optimized trajectory and set of controls
            new_merit = compute_merit(Bcvx, Ucvx, model, penalty_coeff)
			
            # line search
            approx_merit_improve = merit - model_merit
            exact_merit_improve = merit - new_merit
            merit_improve_ratio = exact_merit_improve / float(approx_merit_improve)
                        
            #IPython.embed()

            #info = struct('trust_box_size',trust_box_size);

            print('      approx improve: %.3g. exact improve: %.3g. ratio: %.3g' % (approx_merit_improve, exact_merit_improve, merit_improve_ratio))
            
            if approx_merit_improve < -1e-5:
                print('Approximate merit function got worse (%.3e).' % approx_merit_improve)
                print('Either convexification is wrong to zeroth order, or you are in numerical trouble')
                return B, U, False
            elif approx_merit_improve < cfg.min_approx_improve:
                print('Converged: y tolerance')
                B = Bcvx
                U = Ucvx
                plot.plot_belief_trajectory(B, U, model);
                #pause(0.01);
                return B, U, True
            elif (exact_merit_improve < 1e-2) or (merit_improve_ratio < cfg.improve_ratio_threshold):
                print('sqp 3rd clause')
                trust_box_size = trust_box_size * cfg.trust_shrink_ratio
                break # from trust region loop
            else:
                print('sqp else clause')
                trust_box_size = trust_box_size * cfg.trust_expand_ratio
                B = Bcvx
                U = Ucvx
                plot.plot_belief_trajectory(B, U, model);
                #pause(0.01);
                break # from trust region loop
            
            if trust_box_size < cfg.min_trust_box_size:
                print('Converged: x tolerance\n')
                return B, U, True
        
    sqp_iter += 1
    
# belief function
# b_tp1 = g(b_t, u_t)
def g(b, u, model):
    xDim = model.xDim
    qDim = model.qDim
    rDim = model.rDim
    
    Q = model.Q
    R = model.R
    
    dynamics_func = model.dynamics_func
    obs_func = model.obs_func
    
    x, s = belief.decompose_belief(b, model)
    
    x_tp1 = dynamics_func(x, u, ml.zeros([qDim,1]))

    dyn_varargin = [x, u, ml.zeros([qDim,1])]
    obs_varargin = [dynamics_func(x, u, ml.zeros([qDim,1])), ml.zeros([rDim,1])]

    # dynamics state jacobian
    A = math_util.numerical_jac(dynamics_func, 0, dyn_varargin)
    # dynamics noise jacobian
    M = math_util.numerical_jac(dynamics_func, 2, dyn_varargin)
    # observation state jacobian
    H = math_util.numerical_jac(obs_func, 0, obs_varargin)
    # observation noise jacobian
    N = math_util.numerical_jac(obs_func, 1, obs_varargin)
    
    s_tp1_neg = A*s*((A*s).T) + M*Q*M.T
    K = s_tp1_neg*H.T*ml.linalg.inv(H*s_tp1_neg*H.T + N*R*N.T)
    
    s_tp1 = np.asmatrix(sci.linalg.sqrtm((ml.eye(xDim) - K*H)*s_tp1_neg))
    s_tp1 = ml.real(s_tp1) # needed b/c s_tp1 will have unnecessary complex entries

    b_tp1 = belief.compose_belief(x_tp1, s_tp1, model)

    return b_tp1

# after each iteration in belief_penalty_opt_sqp
# check if done
def constraints_satisfied(B, U, model, tolerance):
    xDim = model.xDim
    xMax = model.xMax
    xMin = model.xMin
    uMax = model.uMax
    uMin = model.uMin
    T = model.T
    done = True

    for t in xrange(0, T-1):
        done &= np.min(B[:,t+1] - g(B[:,t],U[:,t],model) < tolerance)
        done &= np.min(B[1:xDim,t] <= xMax)
        done &= np.min(B[1:xDim,t] >= xMin)
        done &= np.min(U[:,t] <= uMax)
        done &= np.min(U[:,t] >= uMin)

        if not done:
            break
        
    return done
