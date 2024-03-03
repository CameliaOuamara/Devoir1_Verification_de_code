
"""
Classe qui effectue une étude de convergence sur le probleme étudié
La classe calcule la solution numerique du probleme pour plusieurs discretisations geometriques
"""

# Libraries
import numpy as np
import matplotlib.pyplot as plt

# Functions
from profil_concentration import *
from plot import *
from norme_erreur_discretisation import *

class Etude_Convergence():
    def __init__(self, delta_r_vect, delta_t_vect, N_vect, R, critere_convergence, critere_max_iter, schema):
        self.delta_r_vect = delta_r_vect
        self.delta_t_vect = delta_t_vect
        self.N_vect = N_vect
        self.R = R
        self.critere_convergence = critere_convergence
        self.critere_max_iter = critere_max_iter
        self.schema = schema

        
    def Boucle_iterations(self,outputFolder):
        erreur_vect_L1 = np.zeros(len(self.N_vect))
        erreur_vect_L2 = np.zeros(len(self.N_vect))
        erreur_vect_L_inf = np.zeros(len(self.N_vect))
        
        
        for i in range(len(self.N_vect)):

            # Resolution
            if self.schema==1:
                Objet_Concentration = Profil_Concentration(self.delta_r_vect[i], self.delta_t_vect[i], self.N_vect[i], self.R, self.critere_convergence, self.critere_max_iter)
            elif self.schema == 2:
                Objet_Concentration = Profil_Concentration_Centree(self.delta_r_vect[i], self.delta_t_vect[i], self.N_vect[i], self.R, self.critere_convergence, self.critere_max_iter)
            Objet_Concentration.Algorithme_Resolution()

            # Plot
            Objet_Graphique = Plot_Concentration(Objet_Concentration.C, self.N_vect[i])
            Objet_Graphique.Plot_Numerique()
            #Objet_Graphique.Plot_Exact()
            #Objet_Graphique.Save_plot("schema1_"+str(N_vect[i]), "Comparaison de résultat premier schéma, "+str(N_vect[i])+" noeuds")
            #Objet_Graphique.Save_plot(outputFolder+"schema_%d_%d"%(self.schema,self.N_vect[i]), "Comparaison de résultat schéma %d ,%d noeuds"%(self.schema, self.N_vect[i]))
            
            # Erreur
            Objet_Norme_Erreur = Norme_Erreur_Discretisation(Objet_Graphique.C_exact, Objet_Concentration.C[-1,:])
            erreur_vect_L1[i], erreur_vect_L2[i], erreur_vect_L_inf[i] = Objet_Norme_Erreur.Calcul_Norme()

            del Objet_Concentration
            del Objet_Graphique

        return erreur_vect_L1, erreur_vect_L2, erreur_vect_L_inf
        

