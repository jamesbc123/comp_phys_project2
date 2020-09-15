# -*- coding: utf-8 -*-
"""
Created on Tue Sep 15 13:51:29 2020

@author: jamesbc
"""

import numpy as np

class Solver:
      def __init__(self, n, tolerance):
            self.tol = tolerance
            self.tol_reached = False
            self.n = n
            self.max_iter = n**3
            self.A = np.eye((n, n))
            self.A[0,1]=1
            
      def norm(self):
            A_squ = self.A**2
            A_squ_max = A_squ.max()
            self.p, self.k = np.where(A_squ == A_squ_max)
            
            if A_squ_max <= self.tol:
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
            t1 = -tau + np.sqrt(tau**2 + 1)
            t2 = -tau - np.sqrt(tau**2 + 1)
            
            c = 1.0/(1.0 + t1**2)
            s = t1*c
            
            return t1, c, s
            
      def test(self):
            for i in range(0, self.max_iter):
                  self.norm(self)
                  if self.tol_reached == True:
                        break
                  tau = self.calc_tau(self)
                  t, c, s = self.cal_trig(self, tau)
                  
solver = Solver(10, 1e-8)
