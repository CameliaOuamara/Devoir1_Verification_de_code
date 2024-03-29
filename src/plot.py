# -*- coding: utf-8 -*-
"""
Classe qui trace le profil de concentration numerique dans le cylindre suivant le rayon pour une solution donnée
Trace le profil de concentration analytique du probleme
"""
import numpy as np
import matplotlib.pyplot as plt

class Plot_Concentration():
    def __init__(self, C, N):
        # Entrées :
        self.N = N
        self.C = C
        
        # Données :
        self.R = 0.5 # [m]
        self.Ce = 12 # [mol/m3]
        self.Deff = 10**-10 # [m2/s]
        self.S = 8*10**-9 # [mol/m3/s]
        
        self.r = np.linspace(0, self.R, self.N)
        
        self.C_exact = (0.25*self.S*self.R**2*(self.r**2/self.R**2-1)/self.Deff + self.Ce)
        
    def Plot_Numerique(self):
        
        # for i in range(len(self.C)):
        #     plt.scatter(self.r, self.C[i, :], label='i=%d'%i)
        plt.plot(self.r, self.C[-1, :], '.r', label = 'Solution numérique')        
        
    def Plot_Exact(self):
        plt.plot(self.r, self.C_exact, label='Solution exacte')
        
    def Save_plot(self,fileName,plotTitle):
        plt.grid()
        plt.xlabel("r")
        plt.ylabel("C")
        plt.legend()
        plt.title(plotTitle)
        plt.savefig(fileName+".png")
        plt.show()
            







































































