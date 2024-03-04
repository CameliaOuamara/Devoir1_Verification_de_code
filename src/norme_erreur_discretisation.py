"""
Calcul les normes L1 L2 et Linf de l'erreur de discretisation de la solution numerique
"""

import numpy as np

class Norme_Erreur_Discretisation () :
    def __init__(self, solution_exacte, solution_numerique) :
        self.u_exact = solution_exacte.flatten()
        self.u_numerique = solution_numerique.flatten()

    def Calcul_Norme(self) :
        if len(self.u_numerique) != len(self.u_exact):
            u_numerique = np.interp(np.arange(len(self.u_exact))/(len(self.u_exact)-1), np.arange(len(self.u_numerique))/(len(self.u_numerique)-1), self.u_numerique)
            self.Erreur_L1 = np.mean(abs(u_numerique - self.u_exact))
            self.Erreur_L2 = np.sqrt(np.mean((u_numerique - self.u_exact)*(u_numerique - self.u_exact)))
            self.Erreur_Linf = np.max(abs(u_numerique - self.u_exact))
        else:
            self.Erreur_L1 = np.mean(abs(self.u_numerique - self.u_exact))
            self.Erreur_L2 = np.sqrt(np.mean((self.u_numerique - self.u_exact)*(self.u_numerique - self.u_exact)))
            self.Erreur_Linf = np.max(abs(self.u_numerique - self.u_exact))
        return self.Erreur_L1, self.Erreur_L2, self.Erreur_Linf 
