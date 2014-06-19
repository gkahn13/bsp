import numpy as np

import matplotlib.pyplot as plt

import IPython

def plot_mean_shift(P, modes):
    plt.clf()
    plt.cla()
    
    fig = plt.figure(1)
    ax = fig.axes[0]
    ax.set_aspect('equal')
    
    ax.axis([-10, 10, -10, 10])
    
    plt.plot(P[0,:], P[1,:], 'x', color='red')
    
    for mode in modes:
        plt.plot(mode[0], mode[1], 'o', markersize=10.0, color='blue')
    
    plt.show(block=False)
    
    print('Press enter to exit')
    raw_input()