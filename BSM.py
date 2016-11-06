"""
This module solves the Black-Scholes-Merton option pricing problem for
European put options.

This module also includes functions for constructing the initial matrix
from some parameters, saving the initial matrix to a file and saving the
solution to a file.

The module was written for Python 3.* and has been tested on 
Python 3.5.2 and Python 2.7.12.

Authors:  Chris Kiernan, Eoin O'Driscoll, Sean Tully
Version:  2
Date:     6th November 2016
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
# The above line is required for compatibility with Python 2.7.

from sparseSOR import *


def solve_BSM(A, fM, fmod, M, N, X, maxits, omega, machine_epsilon, 
              x_tolerance, r_tolerance):
    """Solve an instance of the Black-Scholes-Merton problem.
    
    Args:
        A (list):  An (N-1) x (N-1) matrix represented as a list of
            N-1 rows, each containing a list of N-1 values as floats.
        fM (list):  A vector of N-1 floats, representing the solution
            at expiry.
        N (integer):  Number of time increments.
        M (int):    Number of price increments.
        maxits (int):  The maximum number of iterations to attempt.
        omega (float):  Relaxation parameter.
        machine_epsilon:  Machine error tolerance for float operations.
        x_tolerance (float):  User-specified tolerance in successive x
            approximations (optional).
        r_tolerance (float):  User-specified tolerance in successive
            residual approximations (optional).
        
    Returns:
        F (list):  An (M+1) x (N-1) matrix represented as a list of
            M+1 rows, each containing a list of N-1 values as floats.
            Each column contains a solution vector, f, with the first
            column containing the solutions at expiry and subsequent
            columns containing sulutions at preceding time increments.
    """
    F = [fM[:]]
    b = fM[1:-1]
    for m in range(1, M + 1):
        b[0] += fmod
        x, k, stopReason = sparse_sor(A, b, N - 1, maxits, omega,
                                       machine_epsilon, x_tolerance,
                                       r_tolerance)
        if x is None:
            print(stopReason + ' for iteration ' + str(k))
            return None
        F.append([X] + x[:] + [0.0])
        b = x[:]
    return F


def construct_BSM(T, X, Smax, M, N, r, sigma):
    """Construct a BMS matrix and vector from parameters.
    
    Args:
        T (float):  Total length of time.
        X (float):  Strike price of option at time T.
        Smax (float):   Maximum price of underlying stock.
        M (int):    Number of price increments.
        N (int):  Number of time increments.
        r (float):  Risk-free interest rate.
        sigma (float):  Volatility of stock price.
        
    Returns:
        A (list):  An (N-1) x (N-1) matrix represented as a list of
            N-1 rows, each containing a list of N-1 values as floats.
        fM (list):  A vector of N+1 floats, representing the solution
            at expiry.
        fmod (float):  Number to be added to the second element of the
            f(m+1) vector before using it to find fm vector.
    """
    # Ensure all float parameters are correctly cast
    T = float(T)
    X = float(X)
    Smax = float(Smax)
    r = float(r)
    sigma = float(sigma)
    
    h = Smax / float(N)
    k = T / float(M)

    # Initialise A and fM with zeros.
    fM = [0.0] * (N + 1)
    A = [[0.0] * (N - 1) for n in range(N - 1)]
    
    # Populate matrix A with values.
    A[0][0] = 1.0 + k * (r + sigma * sigma) # Set first diagonal
    for n in range(1, N - 1):
        # Set diagonal values.
        A[n][n] = 1.0 + k * (r + sigma * sigma * (n + 1) * (n + 1))
        # Set values above diagonal.
        A[n-1][n] = -0.5 * n * k * (n * sigma * sigma + r)
        # Set values below diagonal.
        A[n][n-1] = -0.5 * (n + 1) * k * ((n + 1) * sigma * sigma - r)
        
    # Populate vector fM with values.
    for n in range(N+1):
        fM[n] = max(X - n * h, 0.0) 
    
    # Calculate fmod
    fmod = 0.5 * k * (sigma * sigma - r) * X
    return A, fM, fmod

   
def write_BSM_problem(A, b, n, delimiter=' ', filename='nas_BSM.in'):
    """Writes initial BSM data to a file."""
    with open(filename, 'w') as f:
        data = str(n) + '\n'
        for i in range(n):
            data = data + delimiter.join(str(j) for j in A[i]) + '\n'
        data = data + delimiter.join(str(i) for i in b) + '\n'
        f.write(data)

def write_BSM_solution(F, delimiter=' ', filename='nas_BSM.out'):
    """Writes a matrix, F, to file."""
    with open(filename, 'w') as f:
        data = ''
        n = len(F)
        for i in range(n):
            data = data + delimiter.join(str(j) for j in F[i]) + '\n'
        f.write(data)


