# -*- coding: utf-8 -*-
"""
Classe qui trace le profil de concentration numerique dans le cylindre suivant le rayon pour une solution donnée
Trace le profil de concentration analytique du probleme
"""
import numpy as np
import matplotlib.pyplot as plt

class Plot_Concentration():
    def __init__(self, C, N, sol_MNP, spline_MNP, t):
        # Entrées :
        self.N = N
        self.C = C
        
        # Données :
        self.R = 0.5 # [m]
        self.Ce = 12 # [mol/m3]
        self.Deff = 10**-10 # [m2/s]
        self.S = 8*10**-9 # [mol/m3/s]
        self.t = t
        
        self.r = np.linspace(0, self.R, self.N)
        
        self.C_exact = (0.25*self.S*self.R**2*(self.r**2/self.R**2-1)/self.Deff + self.Ce)
        
        self.sol_MNP = sol_MNP
        self.spline_MNP = spline_MNP
        
        
    def Plot_Numerique(self):
        
        # for i in range(len(self.C)):
        #     plt.scatter(self.r, self.C[i, :], label='i=%d'%i)
        plt.plot(self.r, self.C[-1, :], '.r', label = 'Solution numérique')        
        
    def Plot_MNP(self):
        
        # for i in range(len(self.C)):
        #     plt.scatter(self.r, self.C[i, :], label='i=%d'%i)
        plt.plot(np.linspace(0,max(self.r),len(self.sol_MNP)), self.sol_MNP, '.b', label = 'Solution MNP')    

    def Plot_spline_MNP(self):
        plot_r = np.linspace(0,max(self.r),15*len(self.sol_MNP))
        C_plot = []
        # for r_i in plot_r:
        #     C_plot_i = self.spline_MNP(self.t,r_i)
        #     C_plot.append(C_plot_i)
        C_plot = self.spline_MNP(self.t,plot_r)

        # for i in range(len(self.C)):
        #     plt.scatter(self.r, self.C[i, :], label='i=%d'%i)
        plt.plot(plot_r, C_plot[0], '-g', label = 'Solution spline MNP')    

        
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
            







































































