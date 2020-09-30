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
n_schedule = [50, 100, 150, 200, 250]

toi = pd.DataFrame(columns=["n", "iter", "time"])

for j in range(0,1):
      for n in tqdm(n_schedule):
            h = 1/n
            hh = h*h
            a = 2/hh*np.ones(n)
            b = -1/hh*np.ones(n-1)
            A = np.zeros((n, n))
            
            # Set the diagonal and off-diagonal elements. 
            np.fill_diagonal(A, a)
            np.fill_diagonal(A[:,1:], b)
            np.fill_diagonal(A[1:,:], b)
            
            solver = Solver(n, A, tolerance)
            iterations, A, time_elap = solver.run()
            
            temp = pd.DataFrame({"n": n, "iter": iterations, "time": time_elap}, index=[0])
            toi = toi.append(temp)
      
toi.to_csv('../results/buckling_beam/toi.csv')
