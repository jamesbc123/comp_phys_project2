# -*- coding: utf-8 -*-
"""
Created on Wed Sep 23 13:20:40 2020

@author: jamesbc
"""

import pandas as pd
import numpy as np


toi = np.loadtxt("./cpp_codes/timing.txt", skiprows=1, delimiter=",")
toi = pd.DataFrame(toi, columns=["n","time"])
toi = toi.groupby(["n"]).mean()
print(toi.to_latex())

toi = pd.read_csv("./results/toi.csv")

toi = pd.DataFrame(toi, columns=["n", "time", "iter"])
toi = toi.groupby(["n"]).mean()
print(toi)
print(toi.to_latex())