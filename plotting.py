# -*- coding: utf-8 -*-
"""
Created on Wed Sep 23 13:20:40 2020

@author: jamesbc
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def load_timing_cpp():
      
      toi = np.loadtxt("./cpp_codes/timing.txt", skiprows=1, delimiter=",")
      toi = pd.DataFrame(toi, columns=["n","time"])
      toi = toi.groupby(["n"]).mean()
      print(toi.to_latex())
      
      return

def load_timing_py():
      toi = pd.read_csv("./results/toi.csv")
      
      toi = pd.DataFrame(toi, columns=["n", "time", "iter"])
      toi = toi.groupby(["n"]).mean()
      print(toi)
      print(toi.to_latex())
      
      return

def load_R_cpp():
      n = 10
      toi = np.loadtxt("./cpp_codes/eig_vec_10.txt")
      print(toi.shape)
      eig_vec = toi[0,:]
      h = 1/10
      x = np.arange(0,1,h)
      print(x)
      plt.plot(x, eig_vec)
      plt.xlabel("x")
      plt.ylabel("eigenvector value")
      plt.savefig("./results/eig_vec_10.png")
      return

load_R_cpp()