import matplotlib.pyplot as plt
import numpy as np


# Eksempel, arange (kontrollerer størrelsen på dx)
x1 = np.arange(0, 5, 0.1)
y1 = x1**2

"""
# Eksempel, linspace (kontrollerer antall x-verdier i intervallet)
x2 = np.linspace(0,5,50)
y2 = np.log(x2)
"""

xListe, yListe, xListe2, yListe2 = [1],[1],[2],[2]  # Vilkårlige punkter

# Mal for selve plottet (2D plot):
#import numpy as np
#import matplotlib.pyplot as plt
plot1, = plt.plot(xListe, yListe, '-o', label='liste 1')
plot2, = plt.semilogy(xListe2, yListe2,'.', label='liste 2') # Hvis du vil ha logaritmisk plott på y-aksen. (eller evt. plt.semilogx eller plt.loglog)
plt.legend([plot1, plot2],['Kurve 1', 'Kurve 2'])
# Hvis du får 'iterable'-feil (merk kommaene): plt.legend([(plot1,), (plot2,)], ['Plot 1', 'Plot 2'])
plt.grid()
plt.xlabel(r'$x$ (enhet)')
plt.ylabel(r'$y$ (enhet)')
#plt.xlim(0,max(xliste)*1.05)  # Setter grensene for x-aksen.
#plt.ylim(0,max(yliste)*(1.05)) # Setter grensene for y-aksen.
plt.suptitle('Tittel')
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