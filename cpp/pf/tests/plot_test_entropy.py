import numpy as np
import matplotlib.pyplot as plt
from numpy import matlib as ml

import IPython

def plot(mu, sigma, P, x_dim, M):
    mu = np.matrix(mu)
    mu = np.reshape(mu, (x_dim, 1))
    
    sigma = np.matrix(sigma)
    sigma = np.reshape(sigma, (x_dim, x_dim))
    
    P = np.matrix(P)
    P = np.reshape(P, (x_dim, M))
    
    plt.clf()
    plt.cla()
    
    plt.axis([-6, 6, -6, 6])
    
    for i in xrange(M):
        plt.plot(P[0,i], P[1,i], 'ro')
    
    
    plt.plot(mu[0,0], mu[1,0], 'bx', mew=2.0)
    plot_cov(mu, sigma)
    
    plt.show(block=False)
    plt.pause(.05)
    
    #IPython.embed()
    
def plot_cov(mu, sigma, plotType = 'b-', alpha=1):
    mu = np.asarray(mu)
    sigma = np.asarray(sigma)

    t = np.linspace(-np.pi, np.pi, 2*np.pi/.1)
    x = np.sin(t)
    y = np.cos(t)

    D_diag, V = ml.linalg.eigh(sigma)
    D = np.diag(D_diag)
    A = np.real((np.dot(V,ml.sqrt(D))).T)

    z = np.dot(np.vstack((x.T,y.T)).T,A)

    plt.plot(z[:,0]+mu[0], z[:,1]+mu[1], plotType, alpha=alpha, linewidth=4, mew=2.0)