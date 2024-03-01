"""
Calcul les normes L1 L2 et Linf de l'erreur de discretisation de la solution numerique
"""

import numpy as np

class Norme_Erreur_Discretisation () :
    def __init__(self, solution_exacte, solution_numerique) :
        self.u_exact = solution_exacte
        self.u_numerique = solution_numerique

    def Calcul_Norme(self) :
        self.Erreur_L1 = np.mean(abs(self.u_numerique - self.u_exact))
        self.Erreur_L2 = np.sqrt(np.mean((self.u_numerique - self.u_exact)*(self.u_numerique - self.u_exact)))
        self.Erreur_Linf = np.max(abs(self.u_numerique - self.u_exact))
        return self.Erreur_L1, self.Erreur_L2, self.Erreur_Linf 
