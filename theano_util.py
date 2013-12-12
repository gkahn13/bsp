import numpy as np
from numpy import matlib as ml
import theano
import theano.tensor as T

import IPython

def gradient_flattened(func, x):
    """flat gradient of func wrt x"""
    return T.flatten(T.grad(func, x, disconnected_inputs='warn'))


def jacobian_flattened(func, x):
    """jacobian of func wrt x"""
    func = T.flatten(func)
    idx = T.arange(func.shape[0])

    def grad_i(i):
        return gradient_flattened(func[i], x)

    return theano.map(grad_i, idx)[0]


def jacobian(func, idx, varargin):
    x_orig = varargin[idx]
    def f(x):
        varargin[idx] = x
        val = func(*varargin)
        varargin[idx] = x_orig
        return val

    return jacobian_flattened(func(*varargin), varargin[idx])

# for list of ndarrays [a1, ..., aN]
# returns a1*...*aN
def ndot(ndarray_list):
    prod = ndarray_list[0]
    
    for idx in xrange(1,len(ndarray_list)):
        prod = T.dot(prod, ndarray_list[idx])

    return prod
