# -*- coding: utf-8 -*-
"""
Created on Wed Sep 23 12:01:51 2020

@author: jamesbc
"""

import os
import sys

# Compile all cpp files
all_cpp_codes = "./cpp_codes/*.cpp"
compiler_flags = "-larmadillo"

os.system("echo compiling...")
os.system("g++ -o main.out" + " " + all_cpp_codes + " " + compiler_flags)

os.system("echo executing")
os.system("./main.out")

# Make directory for plotting if it doesn't exist.
path = "./results"
if not os.path.exists(path):
      os.makedirs(path)

