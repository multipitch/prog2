#!/usr/bin/env python
"""
This executable script solves the BSM problem, given some initial
parameters, by calling functions in BSM.py. It then writes the solution to
a file and also writes a 3d surface plot of the solution to file.

The script was written for Python 3.* and has been tested on 
Python 3.5.2 and Python 2.7.12.

Note that a seperate script, BSM.sh does all of the above except for
the creation of the plot.  It doesn't require additional libraries from
numpy, scipy and matplotlib and so exists seperately as a faster and
more portable alternative.

Authors:  Chris Kiernan, Eoin O'Driscoll, Sean Tully
Version:  1
Date:     12th November 2016

Args:
    T (float):  Total length of time.
    X (float):  Strike price of option at time T.
    Smax (float):   Maximum price of underlying stock.
    M (int):    Number of price increments.
    N (int):  Number of time increments.
    r (float):  Risk-free interest rate.
    sigma (float):  Volatility of stock price.
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


# Set some operating parameters for sparse_sor
MAXITS = 1000
OMEGA = 1.3
MACHINE_EPSILON = sys.float_info.epsilon
X_TOLERANCE = 1e-10
R_TOLERANCE = None

# Set default BSM parameters  
T = 90 / 365.0 # Option duration; T is in units of years.
X = 50.0       # Strike price.
Smax = 100.0   # Max stock price.
M = 50         # Number of time increments.
N = 50         # Number of stock price increments.
r = 0.02       # Annualised risk-free interest rate.
sigma = 0.3    # Volatility.

# Read in input parameters (use defaults otherwise).
argCount = len(sys.argv)
if argCount >= 2: 
    T = float(sys.argv[2])
    if argCount >= 3:
        X = float(sys.argv[3])
        if argCount >= 4:
            Smax = float(sys.argv[4])
            if argCount >= 5:
                M = int(sys.argv[5])
                if argCount >= 6:
                    N = int(sys.argv[6])
                    if argCount >= 7:
                        r = float(sys.argv[7])
                        if argCount >= 8:
                            sigma = float(sys.argv[8])  

# Construct BSM Matrix
A, fM, fmod = construct_BSM(T, X, Smax, M, N, r, sigma)

# Solve BSM Matrix
F = solve_BSM(A, fM, fmod, M, N, X, MAXITS, OMEGA, MACHINE_EPSILON,
              X_TOLERANCE, R_TOLERANCE)

# Write solutions to a file.
write_BSM_solution(F, filename='BSM.out')


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
ax.text(50, 60, 60, params)
plt.savefig('BSM.pdf')
plt.close('all')

