"""
This module contains a function to solve systems of linear equations of
the form Ax = b as well as additional functions to I/O and for
performing some matrix operations.

The module was written for Python 3.* and has been tested on 
Python 3.5.2 and Python 2.7.12.

Authors:  Chris Kiernan, Eoin O'Driscoll, Sean Tully
Version:  4
Date:     31st October 2016
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
# The above line is required for compatibility with Python 2.7.

# TO DO:  Implement better error / exception handling


def sparse_sor(A, b, n, maxits, omega, machine_epsilon, x_tolerance=0.0,
               r_tolerance=None, x=None):
    """Solves a system of linear equations of the form Ax = b.
        
    Args:
        A (list):  An n x n matrix represented as a list of n rows,
            each containing a list of n values as floats.
        b (list):  A vector of n floats.
        n (int):  Corresponds to the dimensions of A and b above.
        maxits (int):  The maximum number of iterations to attempt.
        omega (float):  Relaxation parameter.
        machine_epsilon:  Machine error tolerance for float operations.
        x_tolerance (float):  User-specified tolerance in successive x
            approximations (optional).
        r_tolerance (float):  User-specified tolerance in successive
            residual approximations (optional).
        x (list):  Initial guess for solution vector; a list of n floats
            (optional). If blank, x = [1, 1, ... , 1] is assumed.
    
    Returns:
        x (list):  Solution vector; a list of n floats.
        k (int):  The number of iterations taken.
        stopReason (str):  Reason for halting the function
    """
    try:
        k = 0
        
        # Check matrix diagonal values for zeros.
        for i in range(n):
            if A[i][i] == 0:
                return None, k, 'Zero on Diagonal'
        
        # Check for strict row or column dominance.
        if not strict_dd_test(A):
                return None, k, 'Not strictly row or column dominant'        
                    
        val, col, rowStart = make_sparse(A)  # Convert A to CSR format.
                
        if x is None:  # Construct initial x; elements must be != 0.
            x = [1] * n

        px = 2  # Parameter for p-norm calculation on x.
        pr = 2  # Parameter for p-norm calculation on residual.
        
        # Main iteration.
        while k < maxits:
            xOld = x[:]
            for i in range(n):
                spam = 0.0
                for j in range(rowStart[i], rowStart[i + 1]):
                    spam += val[j] * x[col[j]]
                    if col[j] == i:
                        d = val[j]
                x[i] += omega * (b[i] - spam)/d
            k += 1
            
            # Calculate p-norms for x.
            xNorm = p_norm(x, px)
            deltax = []
            for i in range(n):
                deltax.append(x[i] - xOld[i])
            deltaxNorm = p_norm(deltax, px)
            
            # Calculate p-norm of residual.
            Ax = Ab(A, x)
            res = [bi - Axi for bi, Axi in zip(b, Ax)]
            resNorm = p_norm(res, pr)
            
            # Perform halting tests.
            if k > 1 and deltaxNorm > deltaxNormOld:
                return None, k, 'x Sequence divergence' 
            if deltaxNorm <= x_tolerance + 4.0 * machine_epsilon * xNorm:
                return x, k, 'x Sequence convergence'
            if r_tolerance is not None:
                if resNorm <= r_tolerance + 4.0 * machine_epsilon * xNorm:
                    return x, k, 'Residual convergence' 
                           
            deltaxNormOld = deltaxNorm            
        return x, k, 'Max Iterations reached'
        
    except:
        return None, k, 'Cannot proceed'


def make_sparse(A, m=None, n=None):
    """Converts a matrix to compressed sparse row (CSR) format.
    
    Args:
        A (list):  An m x n matrix represented as a list of m rows, each
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


def write_solution(x, stopReason, maxits, k, machine_epsilon, x_tolerance,
                   r_tolerance, filename='nas_Sor.out', delimiter=' ', 
                   printToConsole=None):
    """Write solution to a file.
    
    Args:
        x (list):  Solution vector; a list of n floats.
        stopReason (str):  Reason for halting the sparse_sor function.
        maxits (int):  The maximum number of iterations to attempt.
        k (int):  The number of iterations taken.
        machine_epsilon:  Machine error tolerance for float operations.
        x_tolerance (float):  User-specified tolerance in successive x
            approximations (optional).
        r_tolerance (float or None):  User-specified tolerance in 
            successive residual approximations (if None, this means the 
            residual tolerance test has not been applied.            
        filename (str):  The filename (optional).
        delimiter (str):  The delimiter to use between the values of x
            (optional).
        printToConsole (bool): If True, also print to the console 
            (optional).
    
    """
    titles = ['Stopping reason','Max num of iterations','Number of iterations',
              'Machine epsilon','X seq tolerance']
    params = [str(stopReason), str(maxits), str(k), str(machine_epsilon),
              str(x_tolerance)]
    if r_tolerance is not None:
        titles.append('Residual seq tolerance')
        params.append(str(r_tolerance))

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
    
    # Construct results line using delimiter (if solution obtained).
    if x is not None:
        lines.append(delimiter.join(str(i) for i in x) + '\n')
    
    # Write out the file and optionally print to screen if required.
    with open(filename, 'w') as f:
         for l in lines:
            f.write(l)
            if printToConsole is True: print(l, end="")


def p_norm(v, p=1):
    """Obtain the p-norm of a vector. Defaults to 1-norm.
    
    Args:
        v (list):  Vector (list of floats)
        p (number):  Power (float, or any number that can be converted
            to a float).
    
    Returns:
        norm (float):  The p-norm of the vector
    """
    if p == 1:
        return sum(v)
    p = float(p)
    pInv = 1/p
    spam = 0.0
    for i in range(len(v)):
        if p == 1:
            spam += abs(v[i]) ** p
        else:
            spam += abs(v[i])
    norm = spam ** pInv
    return norm


def strict_dd_test(A, n=None):
    """Checks if an n x n matrix, A is strictly diagonal dominant.
    
    Args:
        A (list):  An n x n matrix represented as a list of n rows,
            each containing a list of n values as floats.
        n (int):  Corresponds to the dimensions of A above (optional).
        
    Returns
        strictlyDiagonalDominant (bool):  True if matrix is strictly
            diagonal dominant (i.e. strictly row dominant or strictly
            column dominant.
    """
    if n is None:
        n = len(A)
    strictlyRowDominant = True
    strictlyDiagonalDominant = True
    for i in range(n):  # Check for strict row dominance.
        spam = 0
        for j in range(n):
            if j != i:
                spam += abs(A[i][j])
        if abs(A[i][i]) < spam:
            strictlyRowDominant = False
            break
    if not strictlyRowDominant: # Check for strict column dominance.
        for j in range(n):
            spam = 0
            for i in range(n):
                if i != j:
                    spam += abs(A[i][j])
            if abs(A[j][j]) < spam:
                strictlyDiagonalDominant = False
                break
    return strictlyDiagonalDominant
    
    
def Ab(A, b):
    """Multiplies an m x n matrix, A, by a vector, b, of length n.
    
    Args:
        A (list):  An l x m matrix represented as a list of l rows, each
            containing a list of m values as floats.
        b (list):  A vector of n floats.
            
    Returns:
        c (list):  A vector of n floats.
    """
    c = [sum(a*b for a, b in zip(aRow, b)) for aRow in A]
    return c

