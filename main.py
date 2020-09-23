# -*- coding: utf-8 -*-
"""
Created on Thu Sep 17 13:06:56 2020

@author: jamesbc
"""

from solver import Solver
import numpy as np
import pandas as pd
from tqdm import tqdm

tolerance = 1e-8
pow = 3

toi = pd.DataFrame(columns=["n", "iter", "time"])

for j in tqdm(range(0,6)):
      for i in range(1, pow):
            n = 10**i
            a = 2*np.ones(n)
            b = -1*np.ones(n-1)
            A = np.zeros((n, n))
            
            # Set the diagonal and off-diagonal elements. 
            np.fill_diagonal(A, a)
            np.fill_diagonal(A[:,1:], b)
            np.fill_diagonal(A[1:,:], b)
            
            solver = Solver(n, A, tolerance)
            iterations, A, time_elap = solver.run()
            
            temp = pd.DataFrame({"n": n, "iter": iterations, "time": time_elap}, index=[0])
            toi = toi.append(temp)
      
toi.to_csv('./Results/toi.csv')
