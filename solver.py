# -*- coding: utf-8 -*-
"""
Created on Tue Sep 15 13:51:29 2020

@author: jamesbc
"""

import numpy as np
from tqdm import tqdm
import time

class Solver:
      def __init__(self, n, A, tolerance):
            self.tol = tolerance
            self.tol_reached = False
            self.n = n
            self.max_iter = n**3
            self.R = np.eye(n)
            self.A = A
            
      def max_off_diag(self):
            ''' Find the greatest off-diagonal element
            wiil need to add a np.where self.p and self.k
            are not equal. THen we will have it. 
            
            '''
            max = 0.0
            for i in range(0, self.n):
                  for j in range(i+1, self.n):
                        if np.abs(self.A[i,j]) > max:
                              max = np.abs(A[i,j])
                              self.p = i
                              self.k = j
           
            if max <= self.tol:
                  self.tol_reached = True
            return
      
      def calc_tau(self):
            p = self.p
            k = self.k
            
            return (self.A[p,p] - self.A[k,k])/(2*self.A[k,p])
      
      def calc_trig(self, tau):
            """return tan, cos and sin of theta.
            
            Not sure why in c++ code morten has an if tau 
            is less than 0. 
            """
            if self.A[self.k,self.p] != 0.0:
                  if tau > 0:
                        t = -tau + np.sqrt(tau**2 + 1.0)
                  else:
                        t = -tau - np.sqrt(tau**2 + 1.0)
                  
                  c = 1.0/(1.0 + t**2)
                  s = t*c
            else:
                  c = 1.0
                  s = 0.0
            return c, s
      
      def rotate(self, c, s):
            '''
            '''
            
            p = self.p
            k = self.k
            a_kk = self.A[k,k]
            a_pp = self.A[p,p]
            
            # Change the matrix elements with k and p indices
            self.A[k,k] = c*c*a_kk - 2.0*c*s*self.A[k,p] + s*s*a_pp
            self.A[p,p] = s*s*a_kk + 2.0*c*s*self.A[k,p] + c*c*a_pp
            self.A[k,p] = 0.0
            self.A[p,k] = 0.0
            
            # Change the rest of the elements
            for i in range(0, self.n):
                  if (i != k and i != p):
                        a_ik = self.A[i,k]
                        a_ip = self.A[i,p]
                        self.A[i,k] = c*a_ik - s*a_ip
                        self.A[k,i] = A[i,k]
                        self.A[i,p] = c*a_ip + s*a_ik
                        self.A[p,i] = A[i,p]
            
            r_ik = self.R[i,k]
            r_ip = self.R[i,p]
            self.R[i,k] = c*r_ik - s*r_ip
            self.R[i,p] = c*r_ip + s*r_ik
            
            return
            
      def run(self):
            start = time.time()
            for i in tqdm(range(0, self.max_iter)):
                  Solver.max_off_diag(self)
                  if self.tol_reached == True:
                        finish = time.time()
                        return i, self.A, (finish-start)
                  tau = Solver.calc_tau(self)
                  c, s = Solver.calc_trig(self, tau)
                  Solver.rotate(self, c, s)
            finish = time.time()
            return i, self.A, (finish-start)
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









