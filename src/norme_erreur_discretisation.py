"""
Calcul les normes L1 L2 et Linf de l'erreur de discretisation de la solution numerique
"""

import numpy as np

class Norme_Erreur_Discretisation () :
    def __init__(self, solution_exacte, solution_numerique) :
        self.u_exact = solution_exacte
        self.u_numerique = solution_numerique

    def Calcul_Norme(self) :
        print("Numerique")
        print(self.u_numerique.shape)
        print("Exact")
        print(self.u_exact.shape)
        matrice = self.u_numerique - self.u_exact
        # Erreur_L1 = np.mean(np.mean(abs(self.u_numerique - self.u_exact), axis=0),axis=1)
        # Erreur_L2 = np.sqrt(np.mean(np.mean((self.u_numerique - self.u_exact)*(self.u_numerique - self.u_exact)),axis=0),axis=1)
        # Erreur_Linf = np.mean(np.max(abs(self.u_numerique - self.u_exact)))
        Erreur_L1 = np.mean(np.mean(abs(matrice)))
        Erreur_L2 = np.sqrt(np.mean(np.multiply(np.mean(matrice,axis=0),np.mean(matrice,axis=0))))
        #Erreur_L2 = np.sqrt(np.mean(np.mean(np.multiply((matrice),(matrice)), axis=1)))
        #Erreur_Linf = np.mean(np.max(abs(matrice)))
        Erreur_Linf = np.max(np.mean(abs(matrice), axis=0))
        # print(Erreur_L1)
        # print(Erreur_L2)
        # print(Erreur_Linf)
        return Erreur_L1, Erreur_L2, Erreur_Linf 
