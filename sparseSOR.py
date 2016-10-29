"""
This module solves systems of linear equations of the form Ax = b.

It reads from an input file containing A and b and outputs x to an
output file, along with some other parameters

It was written for Python 3.* and has been tested on Python 3.5.2 and
Python 2.7.12.
"""
from __future__ import (absolute_import, division, print_function, 
                        unicode_literals)
# The above line is required for compatibility with Python 2.7.

import sys

MAXITS = 1000
MACHINE_EPSILON = sys.float_info.epsilon
USER_EPSILON = 0.0
OMEGA = 1.3

def sparse_sor(A, b, n, maxits, user_epsilon, machine_epsilon,
               omega, x=None):
    """Solves a system of linear equations of the form Ax = b.
        
    Args:
        A (list):  An n x n matrix represented as a list of n rows,
            each containing a list of n values as floats.
        b (list):  A vector of n floats.
        n (int):  Corresponds to the dimensions of A and b above.
        maxits (int):  The maximum number of iterations to attempt.
        user_epsilon:  A user-specified error tolerance.
        machine_epsilon:  Machine error tolerance for float operations.
        omega (float):  Relaxation parameter.
        x (list):  Initial guess for solution vector; a list of n floats
            (optional). If blank, x = [1, 1, ... , 1] is assumed.
    
    Returns:
        x (list):  Solution vector; a list of n floats.
        k (int):  The number of iterations taken.
    """
    val, col, rowStart = make_sparse(A)  # Convert A to CSR format.
    
    # TO DO: Implement some tests on the inputs to make sure they are
    #        suitable.
    
    if x is None:  # Construct initial x; elements must be non-zero.
        x = [1] * n
        
    k = 0
    converged = False  # TO DO: Implement convergence test.
    
    while not converged and k < maxits:
        for i in range(n):
            spam = 0.0
            for j in range(rowStart[i], rowStart[i + 1]):
                spam += val[j] * x[col[j]]
                if col[j] == i:
                    d = val[j]
            x[i] += omega * (b[i] - spam)/d
        k += 1
        
    # TO DO:  Return more information as per requirements of output
    #         file.
    
    return x, k


def make_sparse(A, m=None, n=None):
    """Converts a matrix to compressed sparse row (CSR) format.
    
    Args:
        A (list):  An n x n matrix represented as a list of n rows, each
            containing a list of n values as floats.
        m (int):  Number of rows in A (optional). If m is not
            specified, it is calculated.
        n (int):  Number of columns in A (optional).  If n is not
            specified and m is, n is set equal to m (the matrix is
            assumed to be square).  If neither m nor n are specified, n
            is calculated.
    
    Returns:
        val (list):  A list of the non-zero elements of matrix A in row
            order as floats.
        col (list):  A list of the column indices, corresponding to the
            elements in val, as integers.
        rowStart (list):  A list of integers, indicating where each row
            starts in val and col.
    """
    if n is None: 
        if m is None:  # Calculate dimensions if not specified.
            m = len(A)    
            n = len(A[0])            
        n = m  # Assume square if only one dimension is specified.
        
    val = []
    col = []
    rowStart = []
    rs = True  # Value of rs is True until we have a 'hit' in a row.
    colPos = -1  # The position of the last entry in col.

    # Convert to CSR format.
    for i in range(m):
        rs = True
        for j in range(n):
            if A[i][j] != 0:
                val.append(A[i][j])
                col.append(j)
                colPos += 1
                if rs:
                    rowStart.append(colPos)
                    rs = False
    
    # Note a 'hanging' column position is required:
    rowStart.append(colPos + 1)

    return val, col, rowStart


def read_problem(filename='nas_Sor.in', delimiter=' '):
    """Reads in a structured problem definition file.
    
    The file is assumed to be in the following format:
        - The first line contains an integer, n.
        - The next n lines each contain n numbers, seperated by a 
            delimiter; these lines represent the rows of a square n x n
            matrix, A.
        - The last line contains n numbers, seperated by a delimiter;
            these numbers represent the components of a vector, b, of 
            length n.
    Note that all elements of A and b are recast as floats.
    
    Args:
        filename (str):  The filename (optional)
        delimiter (str):  The delimiter between the numbers on each
            line (optional).
    
    Returns:
        A (list):  An n x n matrix represented as a list of n rows,
            each containing a list of n values as floats.
        b (list):  A vector of n floats.
        n (int):  Corresponds to the dimensions of A and b above.
    """
    A = []
    
    with open(filename, 'r') as f:
        n = int(f.readline())
        for line in range(n):
           A.append(list(float(i) for i in f.readline().split(delimiter)))
        b = list(float(i) for i in f.readline().split(delimiter))
        
    return A, b, n


def write_solution(x, k, filename='nas_Sor.out', delimiter=' ', 
                  printToConsole=None):
    """Write solution to a file.
    
    Args:
        x (list):  Solution vector; a list of n floats.
        k (int):  The number of iterations taken.
        filename (str):  The filename (optional).
        delimiter (str):  The delimiter to use between the values of x
            (optional).
        printToConsole (bool): If True, also print to the console 
            (optional).
    
    """
    includeRST = True # TO DO: Test for this.
    
    # TO DO: Pass in other parameters and replace all the 'spam' strings
    #        below with them.
    
    titles = ['Stopping reason','Max num of iterations','Number of iterations',
              'Machine epsilon','X seq tolerance']
    params = ['spam', str(MAXITS), str(k), str(MACHINE_EPSILON), 'spam']
    if includeRST:
        titles.append('Residual seq tolerance')
        params.append('spam')

    # Pad titles and parameters with spaces to justify text.
    offset = 2
    lines = ['', '']
    lt = len(titles)
    for i in range(lt - 1):
        lines[0] = lines[0] + titles[i].ljust(max(len(titles[i]), 
                                            len(params[i])) + offset)
        lines[1] = lines[1] + params[i].ljust(max(len(titles[i]), 
                                            len(params[i])) + offset)
    lines[0] = lines[0] + titles[lt - 1] + '\n'  # Don't pad last title.
    lines[1] = lines[1] + params[lt - 1] + '\n'  # Don't pad last param.
    
    # Construct results line using delimiter.
    lines.append(delimiter.join(str(i) for i in x) + '\n')
    
    # Write out the file and optionally print to screen if required.
    with open(filename, 'w') as f:
         for l in lines:
            f.write(l)
            if printToConsole: print(l, end="")
   

########################################################################  
#                                                                      #
#     TESTING                                                          #
#                                                                      #
########################################################################
def test():
    """Perform a range of general tests on module functions."""
    
    # Dummy data to validate make_sparse as per p.79 of Chapter 4 of the
    # lecture notes.
    A = [[21.0,  0.0,  0.0, 12.0,  0.0,  0.0],
         [ 0.0,  0.0, 49.0,  0.0,  0.0,  0.0],
         [31.0, 16.0,  0.0,  0.0,  0.0, 23.0],
         [ 0.0,  0.0,  0.0, 85.0,  0.0,  0.0],
         [55.0,  0.0,  0.0,  0.0, 91.0,  0.0],
         [ 0.0,  0.0,  0.0,  0.0,  0.0, 41.0]]

    # Test make_sparse.
    # Note that the integer values in col and rowStart are all one less
    # than those in the lecture notes as python enumerates list items
    # by starting from 0, whereas R starts from 1.
    print('Testing makeSparse()...')
    print('A\t\t:  ' + str(A))
    val, col, rowStart = make_sparse(A)
    print('val\t\t:  ' + str(val))
    print('col\t\t:  ' + str(col))
    print('rowStart\t:  ' + str(rowStart))
         
    # Dummy data to validate sparse_sor; solution should be [3.0, 4.0]
    # A = [[8.0, -1.0],[1.0, 2.0]]
    # b = [20.0, 11.0]
    # n = 2
    
    # Test read_problem.
    A, b, n = read_problem('test.in')

    # Test sparse_sor.
    print('\nTesting sparseSOR()...')
    print('A\t\t:  ' + str(A))
    print('b\t\t:  ' + str(b))
    print('n\t\t:  ' + str(n))
    x, k = sparse_sor(A, b, n, MAXITS, USER_EPSILON, MACHINE_EPSILON, OMEGA)
    print('x\t\t:  ' + str(x))
    print('k\t\t:  ' + str(k))
    
    # Test write_solution.
    write_solution(x, k, 'test.out', printToConsole=True)

# Run test.
test()

