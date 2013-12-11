import numpy as np
from numpy import matlib as ml

import numdifftools as ndt
import theano
from theano import tensor

import IPython

# Numerical Jacobian of func 
# idx specifies the index of argument w.r.t which the Jacobian is computed
# varargin is of type list/tuple that encodes arguments passed to func

# For instance, for y = f(x1, x2, ..., xN)
# numerical_jac(f, 2, x1, x2, ..., xN) computes the Jacobian df/dx2

def numerical_jac(func, idx, varargin):

    step = 1e-6

    x = varargin[idx].copy() # make copy so not altered
    isTensor = type(x) == tensor.TensorVariable
    y = func(*varargin)
    lenx = x.shape[0]
    leny = y.shape[0]
    IPython.embed()
    J = ml.zeros([leny, lenx])
    if isTensor:
        J = tensor.shared(J)

    for i in xrange(0,lenx):
        xhi = x.item(i,0) + step
        xlo = x.item(i,0) - step

        varargin[idx][i,0] = xhi
        yhi = func(*varargin) 
        varargin[idx][i,0] = xlo
        ylo = func(*varargin)
        varargin[idx][i,0] = x[i,0]
    
        J[:,i] = (yhi - ylo)/(xhi - xlo)
        
    return J

# numerical jacobian using numdifftools
def numerical_jac1(func, idx, varargin):
    x_orig = varargin[idx]
    def input_func(x):
        varargin[idx] = x
        val = func(*varargin)
        varargin[idx] = x_orig
        return val
    Jfunc = ndt.Jacobian(input_func)
    IPython.embed()
    return np.asmatrix(Jfunc(np.array(x_orig)))

# numerical jacobian using Theano
def numerical_jac2(func, idx, varargin):
    x_orig = varargin[idx]
    def input_func(x):
        varargin[idx] = x
        val = func(*varargin)
        varargin[idx] = x_orig
        return val
    IPython.embed()
    return theano.gradient.jacobian(input_func, x_orig)


# for list of ndarrays [a1, ..., aN]
# returns a1*...*aN
def ndot(ndarray_list):
    prod = ndarray_list[0]
    
    for idx in xrange(1,len(ndarray_list)):
        prod = np.dot(prod, ndarray_list[idx])

    return prod
