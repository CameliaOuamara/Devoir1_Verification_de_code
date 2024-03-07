
"""
Classe qui effectue une étude de convergence sur le probleme étudié
La classe calcule la solution numerique du probleme pour plusieurs discretisations geometriques
"""

# Libraries
import numpy as np
import matplotlib.pyplot as plt
import time

# Functions
from profil_concentration import *
from plot import *
from norme_erreur_discretisation import *

class Etude_Convergence():
    def __init__(self, delta_r_vect, delta_t_vect, N_vect, R, critere_convergence, critere_max_iter, schema, sol_MNP, spline_bicubic, nombre_nodes_vect):
        self.delta_r_vect = delta_r_vect
        self.delta_t_vect = delta_t_vect
        self.N_vect = N_vect
        self.R = R
        self.critere_convergence = critere_convergence
        self.critere_max_iter = critere_max_iter
        self.schema = schema
        self.sol_MNP = sol_MNP
        self.spline_bicubic = spline_bicubic
        self.nombre_nodes_vect = nombre_nodes_vect

    def Boucle_iterations(self,outputFolder):
        erreur_vect_L1 = np.zeros(len(self.N_vect))
        erreur_vect_L2 = np.zeros(len(self.N_vect))
        erreur_vect_L_inf = np.zeros(len(self.N_vect))
        
        
        for i in range(len(self.N_vect)):

            # Resolution
            if self.schema==1:
                Objet_Concentration = Profil_Concentration_MNP(self.delta_r_vect[i], self.delta_t_vect[i], self.nombre_nodes_vect[i], self.R, self.critere_convergence, self.critere_max_iter[i], self.spline_bicubic)
            elif self.schema == 2:
                Objet_Concentration = Profil_Concentration_Centree_MNP(self.delta_r_vect[i], self.delta_t_vect[i], self.nombre_nodes_vect[i], self.R, self.critere_convergence, self.critere_max_iter[i], self.spline_bicubic)
            time_before_resolution = time.time()
            Objet_Concentration.Algorithme_Resolution()
            time_after_resolution = time.time()
            print("time for solving prob in convergence study:")
            print(time_after_resolution-time_before_resolution)
            '''
            # Plot
            Objet_Graphique = Plot_Concentration(Objet_Concentration.C, self.nombre_nodes_vect[i], self.sol_MNP[-1,:], self.spline_bicubic, self.delta_t_vect[i]*float(self.critere_max_iter[i]))
            Objet_Graphique.Plot_Numerique()
            # Objet_Graphique.Plot_Exact()
            #Objet_Graphique.Plot_MNP()
            Objet_Graphique.Plot_spline_MNP()
            Objet_Graphique.Plot_Numerique_decal(0)

            #Objet_Graphique.Save_plot("schema1_"+str(N_vect[i]), "Comparaison de résultat premier schéma, "+str(N_vect[i])+" noeuds")
            Objet_Graphique.Save_plot(outputFolder+"schema_%d_%d"%(self.schema,self.N_vect[i]), "Comparaison de résultat schéma %d ,%d noeuds"%(self.schema, self.N_vect[i]))
            '''
            # Erreur
            time_before_computing_error = time.time()
            Objet_Norme_Erreur = Norme_Erreur_Discretisation(Objet_Concentration.C_analytic[:,:], Objet_Concentration.C[:,:])
            erreur_vect_L1[i], erreur_vect_L2[i], erreur_vect_L_inf[i] = Objet_Norme_Erreur.Calcul_Norme()
            time_after_computing_error = time.time()
            print("time spent computing error of prob in convergence study:")
            print(time_after_computing_error-time_before_computing_error)
            del Objet_Concentration
            #del Objet_Graphique

        return erreur_vect_L1, erreur_vect_L2, erreur_vect_L_inf
        

