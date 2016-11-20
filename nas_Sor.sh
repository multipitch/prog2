#!/usr/bin/env python
"""
This executable script solves systems of linear equations of the form
Ax = b by calling functions in sparseSOR.py.

The script was written for Python 3.* and has been tested on 
Python 3.5.2 and Python 2.7.12.

Authors:  Chris Kiernan, Eoin O'Driscoll, Sean Tully
Version:  1
Date:     31st October 2016
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
# The above line is required for compatibility with Python 2.7.

import sys
from sparseSOR import *

# Set some fixed parameters.
MAXITS = 1000
OMEGA = 1.3
MACHINE_EPSILON = sys.float_info.epsilon
X_TOLERANCE = 1e-06
R_TOLERANCE = None

# Set some default parameters.
inputFilename = 'nas_Sor.in'
outputFilename = 'nas_Sor.out'
delimiter = ' '

# Read in input filename, or use default.
if len(sys.argv) >= 2:
    inputFilename = str(sys.argv[1])

# Read in output filename, or use default.
if len(sys.argv) >= 3:
    outputFilename = str(sys.argv[2])

# Read in delimiter or use default.
if len(sys.argv) >= 4:
    delimiter  = str(sys.argv[3])

# Read in problem data from file.
A, b, n = read_problem(inputFilename, delimiter)

# Solve the equations.
x, k, stopReason = sparse_sor(A, b, n, MAXITS, OMEGA, MACHINE_EPSILON,
                              X_TOLERANCE, R_TOLERANCE)

# Write result to file.
write_solution(x, stopReason, MAXITS, k, MACHINE_EPSILON, X_TOLERANCE,
                   R_TOLERANCE, outputFilename, delimiter, False)

