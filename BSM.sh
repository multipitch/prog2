#!/usr/bin/env python
"""
This executable script solves the BSM problem, given some initial
parameters, by calling functions in BSM.py and writes the solution to
a file.

The script was written for Python 3.* and has been tested on 
Python 3.5.2 and Python 2.7.12.

Authors:  Chris Kiernan, Eoin O'Driscoll, Sean Tully
Version:  1
Date:     7th November 2016
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
R_TOLERANCE = 1e-10

# Set default BSM parameters  
r = 0.02       # Annualised risk-free interest rate.
sigma = 0.3    # Volatility.
T = 90 / 365.0 # Option duration; T is in units of years.
X = 50.0       # Strike price.
Smax = 100.0   # Max stock price.
M = 30         # Number of time increments.
N = 30         # Number of stock price increments.

# Read in input parameters (use defaults otherwise).
if len(sys.argv) >= 2: r = float(sys.argv[1])
if len(sys.argv) >= 3: sigma = float(sys.argv[2])
if len(sys.argv) >= 4: T = float(sys.argv[3])
if len(sys.argv) >= 5: X = float(sys.argv[4])
if len(sys.argv) >= 6: Smax = float(sys.argv[5])
if len(sys.argv) >= 7: M = float(sys.argv[6])
if len(sys.argv) >= 8: N = float(sys.argv[7])

# Construct BSM Matrix
A, fM, fmod = construct_BSM(T, X, Smax, M, N, r, sigma)

# Solve BSM Matrix
F = solve_BSM(A, fM, fmod, M, N, X, MAXITS, OMEGA, MACHINE_EPSILON,
              X_TOLERANCE, R_TOLERANCE)

# Write solutions to a file.
write_BSM_solution(F, filename='BSM.out')

