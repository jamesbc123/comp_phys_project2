import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re
import os

omegaList = [0.01, 0.5, 1, 5]
# Read the comma-separated data files (two columns, x and y):
directory = "../results/problem_2e/"
data1 = np.loadtxt(directory + "2e_eigvecGS_omega0.01.txt", skiprows=1, delimiter=",")
data2 = np.loadtxt(directory + "2e_eigvecGS_omega0.50.txt", skiprows=1, delimiter=",")
data3 = np.loadtxt(directory + "2e_eigvecGS_omega1.00.txt", skiprows=1, delimiter=",")
data4 = np.loadtxt(directory + "2e_eigvecGS_omega5.00.txt", skiprows=1, delimiter=",")

# Get the columns as lists:
data1 = pd.DataFrame(data1, columns=["rhoList","eigenvectorGS"])
data2 = pd.DataFrame(data2, columns=["rhoList","eigenvectorGS"])
data3 = pd.DataFrame(data3, columns=["rhoList","eigenvectorGS"])
data4 = pd.DataFrame(data4, columns=["rhoList","eigenvectorGS"])

rhoList = data1["rhoList"] # Same for all values of omega_r
eigVecGS1 = data1["eigenvectorGS"]
eigVecGS2 = data2["eigenvectorGS"]
eigVecGS3 = data3["eigenvectorGS"]
eigVecGS4 = data4["eigenvectorGS"]


# Mal for selve plottet (2D plot):
#import numpy as np
#import matplotlib.pyplot as plt
plot1, = plt.plot(rhoList, eigVecGS1, label='omega_0.01')
plot2, = plt.plot(rhoList, eigVecGS2, label='omega_0.5')
plot3, = plt.plot(rhoList, eigVecGS3, label='omega_1')
plot4, = plt.plot(rhoList, eigVecGS4, label='omega_5')

#plt.legend([plot1, plot2, plot3, plot4], [r'$\rho_\text{max} = 0.01$', r'$\rho_\text{max} = 0.5$', r'$\rho_\text{max} = 1$', r'$\rho_\text{max} = 5$'])
# Hvis du får 'iterable'-feil (merk kommaene): plt.legend([(plot1,), (plot2,)], ['Plot 1', 'Plot 2'])
plt.grid()
#plt.xlabel(r'$\rho_\text{max}$')
#plt.ylabel(r'$\psi(\rho)$')
plt.suptitle('Eigenvector of the ground state')
plt.show()

# Eller bare klipp og lim:
"""
plot1, = plt.plot(#fyll inn lister her)
#plt.legend(plot1, )
plt.legend()
plt.grid()
plt.xlabel(r'')   # r means 'render'?
plt.ylabel(r'$$') # $$ for latex
#plt.xlim(0,max(xliste)*1.05)  # Sets the limits for the x axis.
#plt.ylim(0,max(yliste)*(1.05)) # Sets the limits for the y axis.
plt.suptitle('')
plt.show()
"""