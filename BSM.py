"""
This module solves the Black-Scholes-Merton option pricing problem for
European put options.

This module also includes functions for constructing the initial matrix
from some parameters, saving the initial matrix to a file and saving the
solution to a file.

The module was written for Python 3.* and has been tested on 
Python 3.5.2 and Python 2.7.12.

Authors:  Chris Kiernan, Eoin O'Driscoll, Sean Tully
Version:  1
Date:     30th October 2016
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
# The above line is required for compatibility with Python 2.7.

from sparseSOR import *

def solve_BSM(A, fM, N, Smax, M, maxits, omega, machine_epsilon, x_tolerance,
              r_tolerance):
    """Solve an instance of the Black-Scholes-Merton problem.
    
    Args:
        A (list):  An (N-1) x (N-1) matrix represented as a list of
            N-1 rows, each containing a list of N-1 values as floats.
        fM (list):  A vector of N-1 floats, representing the solution
            at expiry.
        N (integer):  Number of time increments.
        Smax (float):   Maximum price of underlying stock.
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
    F = [[None] for m in range(M + 1)]
    F[0] = fM[:]
    for m in range(M):
        fM = sparse_sor(A, fM, N - 1, maxits, omega, machine_epsilon,
                        x_tolerance, r_tolerance)[0]
        if fM is None:
            # TO DO: Raise error.
            print('Uh oh.')
            return None
        F[m+1] = fM[:]
    return F


def construct_BSM(T, X, Smax, M, N, r, sigma):
    """Construct a BMS matrix and vector from parameters.
    
    Args:
        T (float):  Total length of time.
        X (float):  Strike price of option at time T.
        Smax (float):   Maximum price of underlying stock.
        M (int):    Number of price increments.
        N (integer):  Number of time increments.
        r (float):  Risk-free interest rate.
        sigma (float):  Volatility of stock price.
        
    Returns:
        A (list):  An (N-1) x (N-1) matrix represented as a list of
            N-1 rows, each containing a list of N-1 values as floats.
        fM (list):  A vector of N-1 floats, representing the solution
            at expiry.
        N (integer):  Number of time increments.
    """
    h = float(Smax) / float(N)
    k = float(T) / float(M)

    # Initialise A and fM with zeros.
    fM = [0] * (N - 1)
    A = [fM[:] for n in range(N - 1)]
    
    # Populate matrix A with values.
    for n in range(N - 1):
        # Set diagonal values.
        A[n][n] = 1.0 + k * r + k * sigma * sigma * (n + 1) * (n + 1)
        if n > 0:
            # Set values above diagonal.
            A[n-1][n] = -0.5 * n * k * (n * sigma * sigma + r)
            # Set values below diagonal.
            A[n][n-1] = -0.5 * (n + 1) * k * ((n + 1) * sigma * sigma - r)
    
    # Populate matrix B with values.
    for n in range(N-1):
        fM[n] = max(X - n * h, 0.0) 
    return A, fM, N

   
def write_BSM_problem(A, b, n, delimiter=' ', filename='nas_Sor.in'):
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


