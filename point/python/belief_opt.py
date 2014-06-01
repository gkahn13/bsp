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

def belief_opt_penalty_sqp(B, U, model, plotting=True, profile=False):
    cfg = model.sqpParams
      
    trust_box_size = cfg.initial_trust_box_size # The trust region will be a box around the current iterate.
    penalty_coeff = cfg.initial_penalty_coeff # Coefficient of l1 penalties 

    # set profiler to None if you want no profiling
    profiler = util.Profiler() if profile else None

    # The outer loop of the sqp algorithm, which repeatedly minimizes
    # the merit function --- Calls minimize_merit_function defined below

    # After this call, check to see if the
    # constraints are satisfied.
    # - If some constraint is violated, increase penalty_coeff by a factor of cfg.merit_coeff_increase_ratio
    # Reset the trust region size to be larger than cfg.min_trust_box_size, which is used in the termination condition for the inner loop.
    # - If all constraints are satisfied (which in code means if they are satisfied up to tolerance cfg.cnt_tolerance), we're done.

    if profiler:
        profiler.start('totalTime')

    iterCount = 1
    success = False
    while iterCount <= cfg.max_penalty_coeff_increases and not success:    
    
        [B, U, success] = minimize_merit_function(B, U, model, cfg, penalty_coeff, trust_box_size, plotting, profiler)
	
        if plotting:
            plot.plot_belief_trajectory(B, U, model);
    
        success = constraints_satisfied(B, U, model, cfg.cnt_tolerance)
	
        trust_box_size = cfg.initial_trust_box_size

	iterCount += 1

    if profiler:
        profiler.stop('totalTime')
        for name, time in profiler.allTimes():
            print('{0}: {1} seconds'.format(name, time))

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
        merit += penalty_coeff*ml.sum(np.abs(B[:,t+1]-belief.belief_dynamics(B[:,t],U[:,t],None,model)))
    
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
    
    for t in xrange(0,T-1):
        trace_merits.append(model.alpha_belief*belief.cvxpy_sigma_trace(Bcvx[:,t], model))
	#control_merits.append(model.alpha_control*cvxpy.quad_over_lin(Ucvx[:,t],1))
	control_merits.append(model.alpha_control*cvxpy.sum(cvxpy.square(Ucvx[:,t])))
        belief_penalty_merits.append(penalty_coeff*cvxpy.sum(cvxpy.abs(Bcvx[:,t+1] - (gval[t] + G[t]*(Bcvx[:,t]-B[:,t]) + H[t]*(Ucvx[:,t]-U[:,t])))))
    
    trace_merits.append(model.alpha_final_belief*belief.cvxpy_sigma_trace(Bcvx[:,T-1], model))

    merit = sum(trace_merits) + sum(control_merits) + sum(belief_penalty_merits)
    
    return merit

def minimize_merit_function(B, U, model, cfg, penalty_coeff, trust_box_size, plotting=True, profiler=None):
    success = True
    sqp_iter = 1
    
    xDim = model.xDim
    bDim = model.bDim
    uDim = model.uDim
    T = model.T;

    while True:

        # In this loop, we repeatedly construct a linear approximation to the nonlinear belief dynamics constraint
        print('  sqp iter: %i' % sqp_iter)
        
        # Compute merit value using the current iterate (trajectory and set of controls)
        merit = compute_merit(B, U, model, penalty_coeff);
        #print('merit: %f' % merit)
        #print('     current merit: %g' % merit);
        
        if profiler: profiler.start('jacobians')
        # Linearize the belief dynamics constraint 
	gval, G, H = list(), list(), list()
        for t in xrange(0, T-1):
            gval.append(belief.belief_dynamics(B[:,t], U[:,t], None, model))
	    G.append(math_util.numerical_jac(belief.belief_dynamics, 0, [B[:,t], U[:,t], None, model]))
	    H.append(math_util.numerical_jac(belief.belief_dynamics, 1, [B[:,t], U[:,t], None, model]))
        if profiler: profiler.stop('jacobians')

        while True:
            # This is the trust region loop
            # Using the approximations computed above, this loop shrinks
            # the trust region until the progress on the approximate merit
            # function is a sufficiently large fraction of the progress on
            # the exact merit function.
            
            #print('Press enter to continue with SQP')
            #raw_input()

            print('    trust region size: %.3g' % trust_box_size)

            if profiler: profiler.start('cvxpy')

            Bcvx = cvxpy.Variable(bDim, T)
            Ucvx = cvxpy.Variable(uDim, T-1)
            constraints = list()

            if profiler: profiler.start('linearized_compute_merit')
            objective_merit = linearized_compute_merit(B,Bcvx,U,Ucvx,gval,G,H,model,penalty_coeff)
            if profiler: profiler.stop('linearized_compute_merit')
            objective = cvxpy.Minimize(objective_merit)
            
            if profiler: profiler.start('cvxpy constraint creation')

            constraints += [Bcvx[:,0] == B[:,0], # Constraint to ensure that initial belief remains unchanged
                           Bcvx[0:xDim,T-1] == B[0:xDim,T-1], # reach goal at time T
                           Bcvx[0:xDim,:] <= np.tile(model.xMax, (1,T)), # upper x bound
                           Bcvx[0:xDim,:] >= np.tile(model.xMin, (1,T)), # lower x bound
                           Ucvx <= np.tile(model.uMax, (1,T-1)), # upper u bound
                           Ucvx >= np.tile(model.uMin, (1,T-1))] # lower u bound
            
            for t in xrange(0,T-1):
                constraints.append(cvxpy.abs(Bcvx[:,t]-B[:,t]) <= trust_box_size)
                constraints.append(cvxpy.abs(Ucvx[:,t]-U[:,t]) <= trust_box_size)
                
            constraints.append(cvxpy.abs(Bcvx[:,T-1]-B[:,T-1]) <= trust_box_size)

            if profiler: profiler.stop('cvxpy constraint creation')

            problem = cvxpy.Problem(objective, constraints)

            if profiler: profiler.start('cvxpy solve')
	    try:
                cvx_optval = problem.solve(verbose=False)
            except Exception as e:
                print('Failed to solve QP subproblem.')
                print('Error %s' % e)
                IPython.embed()
                return B, U, False
            if profiler: profiler.stop('cvxpy solve')

            if profiler: profiler.stop('cvxpy')
	    
	    Bcvx = np.matrix(Bcvx.value)
            Ucvx = np.matrix(Ucvx.value)

            model_merit = cvx_optval
	    	
            # Compute merit value using the optimized trajectory and set of controls
            new_merit = compute_merit(Bcvx, Ucvx, model, penalty_coeff)
	    #print('     new merit: %g' % new_merit);
			
            # line search
            approx_merit_improve = merit - model_merit
            exact_merit_improve = merit - new_merit
            merit_improve_ratio = exact_merit_improve / (approx_merit_improve)
                        
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
                if plotting:
                    plot.plot_belief_trajectory(B, U, model)
                return B, U, True
            elif (exact_merit_improve < 1e-2) or (merit_improve_ratio < cfg.improve_ratio_threshold):
                trust_box_size = trust_box_size * cfg.trust_shrink_ratio
            else:
                trust_box_size = trust_box_size * cfg.trust_expand_ratio
                B = Bcvx
                U = Ucvx
                if plotting:
                    plot.plot_belief_trajectory(B, U, model)
                break # from trust region loop
            
            if trust_box_size < cfg.min_trust_box_size:
                print('Converged: x tolerance\n')
                return B, U, True
        
    	sqp_iter += 1
    
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
    constraint_violations = 0
    for t in xrange(0, T-1):
	constraint_violations += ml.sum(np.abs(B[:,t+1] - belief.belief_dynamics(B[:,t],U[:,t],None,model)))
	done &= constraint_violations < tolerance
	done &= np.max(B[1:xDim,t] <= xMax)
	done &= np.min(B[1:xDim,t] >= xMin)
	done &= np.max(U[:,t] <= uMax)
	done &= np.min(U[:,t] >= uMin)

	if not done:
            break

    #print('Constraint violations: %g' % constraint_violations)
    return done
