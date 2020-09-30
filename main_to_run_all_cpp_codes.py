# -*- coding: utf-8 -*-
"""
Created on Wed Sep 23 12:01:51 2020

@author: jamesbc
"""

''' Sor far I couldn't get this script to run as
expected. Instead please use the make file in the 
cpp_codes folder to execute the code and then 
the plotting script in py_codes.
'''
import os
import sys

# Make directory for plotting if it doesn't exist.
path = "./results/buckling_beam"
if not os.path.exists(path):
      os.makedirs(path)
      
# Compile all cpp files
all_cpp_codes = "./cpp_codes/*.cpp"
compiler_flags = "-larmadillo"

os.system("echo compiling...")
os.system("g++ -o main.out" + " " + all_cpp_codes + " " + compiler_flags)

os.system("echo executing")
os.system("./main.out")

os.system("echo creating plots...")
os.system("./py_codes/plotting.py")


