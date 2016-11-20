#!/usr/bin/env python
"""
This executable script solves the BSM problem, given some initial
parameters, by calling functions in BSM.py and writes the solution to
a file.

Note that a seperate script, BSM_plot.sh does all of the above and 
additionally creates a 3d surface plot.  It requires a number of
additional third-party libraries to create the plot.

The script was written for Python 3.* and has been tested on 
Python 3.5.2 and Python 2.7.12.

Authors:  Chris Kiernan, Eoin O'Driscoll, Sean Tully
Version:  2
Date:     12th November 2016

Args:
    T (float):  Total length of time (optional).
    X (float):  Strike price of option at time T (optional).
    Smax (float):   Maximum price of underlying stock (optional).
    M (int):    Number of price increments (optional).
    N (int):  Number of time increments (optional).
    r (float):  Risk-free interest rate (optional).
    sigma (float):  Volatility of stock price (optional).
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
# The above line is required for compatibility with Python 2.7.

import sys
from BSM import *

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

