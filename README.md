# prog2
MIS40530 Programming Assignment 2:  Linear Equations


###Main Objectives
* Read data from an input file containing A and b.
* Solve Ax = b using the Sparse-SOR algorithm and the sparse matrix data structure for SOR, as given in lectures.
* Write the computed solution vector x, together with the reason for stopping and other information, to an output file.
* Numerically solve the Black-Scholes-Merton (BSM) equation.


###Files
####BSM.py
Constructs a BSM problem from some parameters and solves the BSM equation by calling sparseSOR.py.
####BSM_test.py
Unit tests for BSM.py.
####sparseSOR.py
Solves Ax = b for sparse matrices using the Sparse-SOR algorithm. Also includes functions for I/O and matrix operations.
####sparseSOR_test.py
Unit tests for sparseSOR.py.
####nas_Sor.sh
Executable script which reads from an input file, solves Ax = b (using sparseSOR.py) and writes to an output file.
####makefile
Use 'make clean' to remove output files and temporary files.
