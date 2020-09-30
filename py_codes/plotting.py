# -*- coding: utf-8 -*-
"""
Created on Wed Sep 23 13:20:40 2020

@author: jamesbc
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re
import os

# Change the font size of all figures.
plt.rcParams.update({'font.size': 22, 'figure.figsize': (20,10)})


def plot_timing_table(directory):
      """ Prints the table of times to latex and plots arma vs jacobi
      """
      toi = np.loadtxt(directory + "timing.txt", skiprows=1, delimiter=",")
      toi = pd.DataFrame(toi, columns=["n","time used jacobi", "time used arma"])
      toi = toi.groupby(["n"], as_index=False).mean()
      
      
      plt.plot(toi["n"], toi["time used arma"], label = "armadillo solver")
      plt.plot(toi["n"], toi["time used jacobi"], label = "jacobi solver")
      plt.legend()
      plt.xlabel("n, dimensionality of matrix")
      plt.ylabel("time used (seconds)")
      plt.savefig(directory + "time_plot.png")
      plt.close()
      
      toi1 = plot_iter_table(directory, False)
      toi = pd.concat((toi, toi1['nbr of iterations']), axis=1)
      
      print(toi.to_latex(index=False))
      return

def plot_iter_table(directory, plot):
      """ Prints the table of iterations to latex and plots the number of 
      iterations and iterations divided by n.
      """
      toi = np.loadtxt(directory + "iterations.txt", skiprows=1, delimiter=",")
      
      toi = pd.DataFrame(toi, columns=["n", "nbr of iterations"])
      toi = toi.groupby(["n"], as_index=False).mean()
      
      if plot:
            plt.plot(toi["n"], toi["nbr of iterations"])
            plt.xlabel("n, dimensionality of matrix")
            plt.ylabel("number of iterations")
            plt.savefig(directory + "iterations.png")
            plt.close()
            
            plt.plot(toi["n"], np.log10(toi["nbr of iterations"])/np.log10(toi["n"]))
            plt.xlabel("n, dimensionality of matrix")
            plt.ylabel("log(iterations) base n")
            plt.savefig(directory + "log_base_n.png")
      plt.close()
      
      return toi
      
def plot_timing_table_py(directory):
      """ Prints to latex the results from python class, 
      i.e n, time used, and nbr of iterations.
      """
      toi = pd.read_csv(directory + "toi.csv")
      toi = pd.DataFrame(toi, columns=["n", "time", "iter"])
      toi = toi.groupby(["n"], as_index=False).mean()
      print(toi.to_latex(index=False))
      
      return

def plot_eig_vec(directory):
      ''' Plots the first eigenvector for each n.
      '''
      dict_nested = {'num':{}, 'ana':{}}
      
      dict = dict_nested['num']
      for filename in os.listdir(directory):
            if re.search('num_eigvec_[0-9]+', filename):
                  # Find n from the file name
                  n = int(re.findall('([0-9]+)', filename)[0])
                  toi = np.loadtxt(os.path.join(directory,filename))
                  
                  #extract the column vector and normalise them.
                  eig_vec = toi[:,0]
                  eig_vec = eig_vec / eig_vec.sum()
                  
                  rho_min = 0
                  rho_max = 1
                  h = (rho_max-rho_min)/n
                  x = np.linspace(h,rho_max,n-1)
                  
                  dict['{}'.format(n)] = (x, eig_vec)
    
                  plt.plot(x, eig_vec, label = "n={}".format(n))
                  plt.xlabel("\u03C1")
                  plt.ylabel("first eigenvector")
                  plt.legend()
      plt.savefig("../results/buckling_beam/numerical_eigenvectors.png")
      plt.close()
      
      dict = dict_nested['ana']
      for filename in os.listdir(directory):
            if re.search('analytic_eigvec_[0-9]+', filename):
                  # Find n from the file name
                  n = int(re.findall('([0-9]+)', filename)[0])
                  toi = np.loadtxt(os.path.join(directory,filename))
                  
                  eig_vec = toi[:,0]
                  eig_vec = eig_vec / eig_vec.sum()
                  
                  rho_min = 0
                  rho_max = 1
                  h = (rho_max-rho_min)/n
                  x = np.linspace(h,rho_max,n-1)
                  
                  dict['{}'.format(n)] = (x, eig_vec)
                  
                  plt.plot(x, eig_vec, label = "n={}".format(n))
                  plt.xlabel("\u03C1")
                  plt.ylabel("first eigenvector")
                  plt.legend()
      plt.savefig("../results/buckling_beam/analytic_eigvectors.png")
      plt.close()
      
      # Plots the first eigenvector for the analytic and numerical
      # just for one n, for comparison.
      
      x_ana, eigvec_ana = dict_nested['ana']['50']
      x_num, eigvec_num = dict_nested['num']['50']
      plt.plot(x_ana, eigvec_ana/eigvec_ana.sum(), label= "analytic eigenvector", marker = "v")
      plt.plot(x_num, eigvec_num/eigvec_num.sum(), label="numerical eigenvector", marker ="^") 
      plt.xlabel("\u03C1")
      plt.ylabel("eigenvector for n = 50")
      plt.legend(loc="upper right")
      plt.show()
      plt.savefig("../results/buckling_beam/numerical_vs_analytic.png")
     
      return

directory = "../results/buckling_beam/"

plot_eig_vec(directory)
plot_timing_table(directory)
toi = plot_iter_table(directory, plot=True)
plot_timing_table_py(directory)








