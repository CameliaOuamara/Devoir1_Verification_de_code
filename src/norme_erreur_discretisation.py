"""
Calcul les normes L1 L2 et Linf de l'erreur de discretisation de la solution numerique
"""

import numpy as np

class Norme_Erreur_Discretisation () :
    def __init__(self, solution_exacte, solution_numerique) :
        self.u_exact = solution_exacte
        self.u_numerique = solution_numerique

    def Calcul_Norme(self) :
        if len(self.u_numerique) != len(self.u_exact):
            # u_numerique = np.interp(np.arange(len(self.u_exact))/(len(self.u_exact)-1), np.arange(len(self.u_numerique))/(len(self.u_numerique)-1), self.u_numerique)
            u_exact = np.interp(np.arange(len(self.u_numerique))/(len(self.u_numerique)-1), np.arange(len(self.u_exact))/(len(self.u_exact)-1), self.u_exact)
            # self.Erreur_L1 = np.mean(abs(u_numerique - self.u_exact))
            # self.Erreur_L2 = np.sqrt(np.mean((u_numerique - self.u_exact)*(u_numerique - self.u_exact)))
            # self.Erreur_Linf = np.max(abs(u_numerique - self.u_exact))
            self.Erreur_L1 = np.mean(abs(self.u_numerique - u_exact))
            self.Erreur_L2 = np.sqrt(np.mean((self.u_numerique - u_exact)*(self.u_numerique - u_exact)))
            self.Erreur_Linf = np.max(abs(self.u_numerique - u_exact))            
        else:
            self.Erreur_L1 = np.mean(abs(self.u_numerique - self.u_exact))
            self.Erreur_L2 = np.sqrt(np.mean((self.u_numerique - self.u_exact)*(self.u_numerique - self.u_exact)))
            self.Erreur_Linf = np.max(abs(self.u_numerique - self.u_exact))
        return self.Erreur_L1, self.Erreur_L2, self.Erreur_Linf 
    
class Norme_Erreur_Discretisation_MMS () :
    def __init__(self, f_T_MMS, u_numerique, delta_t, R) :
        self.f_T_MMS = f_T_MMS
        self.u_numerique = u_numerique
        self.delta_t = delta_t
        self.R = R
        self.t = np.linspace(0,(u_numerique.shape[0]-1)*delta_t, u_numerique.shape[0])
        self.r = np.linspace(0,R,u_numerique.shape[1])
        self.ri, self.ti = np.meshgrid(self.r, self.t)

    def Calcul_Norme(self) :
        self.Erreur_matrice = abs(self.u_numerique - self.f_T_MMS(self.ri,self.ti))
        self.Erreur_L1 = np.mean(abs(self.u_numerique - self.f_T_MMS(self.ri,self.ti)))
        self.Erreur_L2 = np.sqrt(np.mean((self.u_numerique - self.f_T_MMS(self.ri,self.ti))*(self.u_numerique - self.f_T_MMS(self.ri,self.ti))))
        self.Erreur_Linf = np.max(abs(self.u_numerique - self.f_T_MMS(self.ri,self.ti)))
        
     
        return self.Erreur_L1, self.Erreur_L2, self.Erreur_Linf 