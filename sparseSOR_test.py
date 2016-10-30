"""
This module contains a series of tests for functions in sparseSOR.py.

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
from sparseSOR import *

MAXITS = 1000
OMEGA = 1.3
MACHINE_EPSILON = sys.float_info.epsilon
X_TOLERANCE = 1e-06
R_TOLERANCE = 1e-06

def test_gen():
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
    print('Testing make_sparse()...')
    print('A\t\t:  ' + str(A))
    val, col, rowStart = make_sparse(A)
    print('val\t\t:  ' + str(val))
    print('col\t\t:  ' + str(col))
    print('rowStart\t:  ' + str(rowStart))
    
    # Test read_problem
    print('\nTesting read_problem...')
    A, b, n = read_problem('sparseSOR_test.in')
    print('A\t\t:  ' + str(A))
    print('b\t\t:  ' + str(b))
    print('n\t\t:  ' + str(n))
    
    # Test sparse_sor. Solution should be [3.0, 4.0]
    print('\nTesting sparse_sor...')
    x, k, stopReason = sparse_sor(A, b, n, MAXITS, OMEGA, MACHINE_EPSILON,
                                  X_TOLERANCE, R_TOLERANCE)
    print('x\t\t:  ' + str(x))
    print('k\t\t:  ' + str(k))
    
    # Test write_solution.
    print('\nTesting write_solution...')
    write_solution(x, stopReason, MAXITS, k, MACHINE_EPSILON, X_TOLERANCE,
                   R_TOLERANCE, 'sparseSOR_test.out', ' ', True)


def test_dd():
    """Test the diagonal dominance test function"""
    print('\nTesting dd...')
    A = [[8,2,1],[3,10,4],[1,1,6]]
    print('Both strictly row and strictly column dominant:')
    print('\tstrictlyDiagonalDomininant = ' + str(strict_dd_test(A)))

    A = [[8,2,1],[5,10,4],[4,1,6]]
    print('Strictly column dominant only:')
    print('\tstrictlyDiagonalDomininant = ' + str(strict_dd_test(A)))

    A = [[8,2,3],[3,10,4],[1,1,6]]
    print('Strictly row dominant only:')
    print('\tstrictlyDiagonalDomininant = ' + str(strict_dd_test(A)))

    A = [[1,3,4],[8,1,5],[6,7,2]]
    print('Neither strictly row nor strictly column dominant:')
    print('\tstrictlyDiagonalDomininant = ' + str(strict_dd_test(A)))

# Run tests
test_gen()
test_dd()

