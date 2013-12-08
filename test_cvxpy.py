import numpy as np
import cvxpy

import IPython

"""
x = cvxpy.Variable(name='x')
y = cvxpy.Variable(name='y')

c0 = x+y <= 1
c1 = x >= 0
c2 = y >=0

c = [c0, c1, c2]

objective = cvxpy.Maximize(1.5*x+y)
p = cvxpy.Problem(objective, c)
"""

x = cvxpy.Variable(2,1,name='x')
y = cvxpy.Variable(2,1,name='y')

c0 = x[0,0] + x[1,0] <= 1
c1 = x >= np.zeros([2,1])
c2 = x[0,0] == .3

c = [c0, c1, c2]

objective = cvxpy.Maximize(cvxpy.sum(x))
p = cvxpy.Problem(objective, c)

IPython.embed()
