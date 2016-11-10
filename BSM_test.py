"""
This module contains a series of tests for functions in BSM.py.

The module was written for Python 3.* and has been tested on 
Python 3.5.2 and Python 2.7.12.

Authors:  Chris Kiernan, Eoin O'Driscoll, Sean Tully
Version:  2
Date:     6th November 2016
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
from scipy.sparse.linalg import spsolve
from scipy.sparse import csr_matrix

MAXITS = 1000
OMEGA = 1.3
MACHINE_EPSILON = sys.float_info.epsilon
X_TOLERANCE = 1e-10
R_TOLERANCE = 1e-10
  
T = 90 / 365.0 # Option duration; T is in units of years.
X = 50.0       # Strike price.
Smax = 100.0   # Max stock price.
M = 50         # Number of time increments.
N = 50         # Number of stock price increments.
r = 0.02       # Annualised risk-free interest rate.
sigma = 0.3    # Volatility.


# Run model with above parameters.
A, fM, fmod = construct_BSM(T, X, Smax, M, N, r, sigma)
F = solve_BSM(A, fM, fmod, M, N, X, MAXITS, OMEGA, MACHINE_EPSILON,
              X_TOLERANCE, R_TOLERANCE)
              

# Solve using scipy.sparse.linalg.spsolve instead of sparse_sor
'''
def scipy_BSM(A, fM, fmod):
    Acsr = csr_matrix(A)
    F = [fM[:]]
    b = fM[1:-1]
    for m in range(1, M + 1):
        b[0] += fmod
        x = spsolve(Acsr, b).tolist()
        if x is None:
            return None
        F.append([X] + x[:] + [0.0])
        b = x[:]
    return F
    
F = scipy_solve(A, fM, fmod)
'''


# Write solutions to a file.
write_BSM_solution(F, filename='BSM_test.out')

# Create a surface plot showing option price as a function of the price
# of the underlying stock and the number of days until maturity.
h = Smax / float(N)
k = T / float(M)
Z = np.array(F)
x = np.linspace(0, Smax, N + 1)
y = np.linspace(0, 365 * T, M + 1)
X, Y = np.meshgrid(x, y)
fig = plt.figure()
ax = fig.gca(projection='3d', azim=-60., elev=20.)
#ax.plot_surface(X, Y, np.flipud(Z), rstride=1, cstride=1, cmap=cm.coolwarm,
#                linewidth=0, antialiased=False) # plot time to maturity
ax.plot_surface(X, Y, np.flipud(Z), rstride=1, cstride=1, cmap=cm.coolwarm,
                linewidth=0, antialiased=False) # plot time
ax.set_xlabel(r'stock price, $S$ (\$)')
#ax.set_ylabel(r'time to maturity, $T - t$ (days)')
ax.set_ylabel(r'time, $T$ (days)')
ax.set_zlabel(r'option price, $f$ (\$)')
ax.set_yticks(range(0, int(T * 365), 20))
params = r'$r$ = ' + str(r) + r'; $\sigma$ = ' + str(sigma)
ax.text(50, 50, 60, params)
plt.savefig('BSM_test.pdf')
plt.close('all')

