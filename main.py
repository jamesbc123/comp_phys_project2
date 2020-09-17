# -*- coding: utf-8 -*-
"""
Created on Thu Sep 17 13:06:56 2020

@author: jamesbc
"""

from solver import Solver
import numpy as np

tolerance = 1e-8
n = 100
a = 2*np.ones(n)
b = -1*np.ones(n-1)
A = np.zeros((n, n))

# Set the diagonal and off-diagonal elements. 
np.fill_diagonal(A, a)
np.fill_diagonal(A[:,1:], b)
np.fill_diagonal(A[1:,:], b)

solver = Solver(n, A, tolerance)
iterations, A, time_elap = solver.run()

print(iterations)
print(time_elap)
print(A.shape)