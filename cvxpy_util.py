import numpy as np
from numpy import matlib as ml
import cvxpy

import IPython

# converts cvxpy matrix to vector
def vectorize1(A):
    rows, cols = A.size

    Avec = A[:,0]

    for col in xrange(1,cols):
        Avec = cvxpy.vstack(Avec,A[:,col])

    return Avec

# converts cvxpy matrix to vector
def vectorize(A):
    rows, cols = A.size
    constraints = list()

    Avec = cvxpy.Variable(rows*cols, 1)

    for col in xrange(0,cols):
         constraints.append( Avec[rows*col:rows*(col+1)] == A[:,col] )

    return Avec, constraints


def sum_square(A):
    rows, cols = A.size

    return sum([cvxpy.quad_over_lin(A[:,i],1) for i in xrange(cols)])
    
def vars_in_constraints_affine(constraints):
    variables = list()
    for c in constraints:
        rh_affine = c.rh_exp.curvature.AFFINE.is_affine()
        lh_affine = c.lh_exp.curvature.AFFINE.is_affine()

        if not rh_affine or not lh_affine:
            IPython.embed()
            return False
        
    return True
