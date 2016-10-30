"""
This module contains a series of tests for functions in BSM.py.

The module was written for Python 3.* and has been tested on 
Python 3.5.2 and Python 2.7.12.

Authors:  Chris Kiernan, Eoin O'Driscoll, Sean Tully
Version:  1
Date:     30th October 2016
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
# The above line is required for compatibility with Python 2.7.

import sys
from BSM import *
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np

MAXITS = 1000
OMEGA = 1.3
MACHINE_EPSILON = sys.float_info.epsilon
X_TOLERANCE = 1e-06
R_TOLERANCE = 1e-06
  
T = 12.0
X = 100.0
Smax = 150.0
M = 12
N = 151
r = 0.04
sigma = 0.27

# Run model with above parameters.
A, b, fM = construct_BSM(T, X, Smax, M, N, r, sigma)
F = solve_BSM(A, b, fM, Smax, M, MAXITS, OMEGA, MACHINE_EPSILON, X_TOLERANCE,
              R_TOLERANCE)

# Write solutions to a file.
write_BSM_solution(F, filename='BSM_test.out')

# Create a surface plot showing option value as a function of the price
# of the underlying stock and the number of days until expiry.
Z = np.array(F)
x = np.arange(0,N-1,1)
y = np.arange(0,M+1,1)
X, Y = np.meshgrid(x, y)
fig = plt.figure()
ax = fig.gca(projection='3d')
surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)
xLabel = ax.set_xlabel('Stock price')
yLabel = ax.set_ylabel('Days to expiry')
zLabel = ax.set_zlabel('Strike price')
plt.savefig('BSM_test.pdf')
plt.close('all')

