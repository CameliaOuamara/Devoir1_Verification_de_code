# -*- coding: utf-8 -*-
"""
Classe qui trace le profil de concentration numerique dans le cylindre suivant le rayon pour une solution donnée
Trace le profil de concentration analytique du probleme
"""
import numpy as np
import matplotlib.pyplot as plt

class Plot_Concentration():
    def __init__(self, C, N, spline_MNP, f_T_MMS, delta_t):
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
        
        self.spline_MNP = spline_MNP
        self.delta_t = delta_t
        self.t = np.linspace(0,(self.C.shape[0]-1)*delta_t, self.C.shape[0])
        self.r = np.linspace(0,self.R,C.shape[1])
        self.ri, self.ti = np.meshgrid(self.r, self.t)
        self.f_T_MMS = f_T_MMS
        
    def Plot_Numerique(self):
        
        # for i in range(len(self.C)):
        #     plt.scatter(self.r, self.C[i, :], label='i=%d'%i)
        for i in range(self.C.shape[0]):
            plt.plot(self.r, self.C[i, :], '1r', label = 'Solution numérique')        
        
    def Plot_MNP(self):
        
        # for i in range(len(self.C)):
        #     plt.scatter(self.r, self.C[i, :], label='i=%d'%i)
        for i in range(len(self.t)):
             plt.plot(self.r, self.spline_MNP(self.delta_t*i,self.r)[0,:], '2b', label = 'Solution MNP')
            
    def Plot_MMS(self):
        
        # for i in range(len(self.C)):
        #     plt.scatter(self.r, self.C[i, :], label='i=%d'%i)
        for i in range(len(self.t)):
            plt.plot(self.r, self.f_T_MMS(self.r,self.delta_t*i), '2b', label = 'Solution MMS')   
        
    def Plot_Exact(self):
        plt.plot(self.r, self.C_exact, label='Solution exacte')
        
    def Save_plot(self,fileName,plotTitle):
        plt.grid()
        plt.xlabel("r")
        plt.ylabel("C")
        # plt.legend()
        plt.title(plotTitle)
        plt.savefig(fileName+".png")
        plt.show()
            







































































