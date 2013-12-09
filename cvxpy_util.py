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

    value = 0
    for i in xrange(cols):
        value += cvxpy.quad_over_lin(A[:,i],1)

    return value
