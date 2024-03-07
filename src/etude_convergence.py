
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
    def __init__(self, delta_r_vect, delta_t_vect, N_vect, R, critere_convergence, critere_max_iter, schema, sol_MNP, spline_bicubic):
        self.delta_r_vect = delta_r_vect
        self.delta_t_vect = delta_t_vect
        self.N_vect = N_vect
        self.R = R
        self.critere_convergence = critere_convergence
        self.critere_max_iter = critere_max_iter
        self.schema = schema
        self.sol_MNP = sol_MNP
        self.spline_bicubic = spline_bicubic

    def Boucle_iterations(self):
        erreur_vect_L1 = np.zeros(len(self.N_vect))
        erreur_vect_L2 = np.zeros(len(self.N_vect))
        erreur_vect_L_inf = np.zeros(len(self.N_vect))
        
        
        for i in range(len(self.N_vect)):

            # Resolution
            if self.schema==1:
                Objet_Concentration = Profil_Concentration_MNP(self.delta_r_vect[i], self.delta_t_vect[-1], self.N_vect[i], self.R, self.critere_convergence, self.critere_max_iter, self.spline_bicubic)
            elif self.schema == 2:
                Objet_Concentration = Profil_Concentration_Centree_MNP(self.delta_r_vect[i], self.delta_t_vect[-1], self.N_vect[i], self.R, self.critere_convergence, self.critere_max_iter, self.spline_bicubic)
            Objet_Concentration.Algorithme_Resolution()

            # Plot
            Objet_Graphique = Plot_Concentration(Objet_Concentration.C, self.N_vect[i], self.sol_MNP, 0.0, 0.0)
            Objet_Graphique.Plot_Numerique()
            # Objet_Graphique.Plot_Exact()
            Objet_Graphique.Plot_MNP()

            #Objet_Graphique.Save_plot("schema1_"+str(N_vect[i]), "Comparaison de résultat premier schéma, "+str(N_vect[i])+" noeuds")
            Objet_Graphique.Save_plot("schema %d_%d"%(self.schema,self.N_vect[i]), "Comparaison de résultat schéma %d ,%d noeuds"%(self.schema, self.N_vect[i]))
            
            # Erreur
            Objet_Norme_Erreur = Norme_Erreur_Discretisation(self.sol_MNP, Objet_Concentration.C[-1,:])
            erreur_vect_L1[i], erreur_vect_L2[i], erreur_vect_L_inf[i] = Objet_Norme_Erreur.Calcul_Norme()

            del Objet_Concentration
            del Objet_Graphique

        return erreur_vect_L1, erreur_vect_L2, erreur_vect_L_inf
    
    
    
class Etude_Convergence_MMS_spatial():
    def __init__(self, delta_r_vect, delta_t, N_vect, R, t_final, schema):
        self.delta_r_vect = delta_r_vect
        self.delta_t = delta_t
        self.N_vect = N_vect
        self.R = R
        self.t_final = t_final
        self.schema = schema


    def Boucle_iterations(self):
        erreur_vect_L1 = np.zeros(len(self.N_vect))
        erreur_vect_L2 = np.zeros(len(self.N_vect))
        erreur_vect_L_inf = np.zeros(len(self.N_vect))
        
        Deff = 1.0e3
        k = 4.0e3
        for i in range(len(self.N_vect)):

            # Resolution
            if self.schema==1:
                Objet_Concentration = Profil_Concentration_MMS(self.delta_r_vect[i], self.delta_t, self.N_vect[i], self.R, self.t_final, Deff, k)
            elif self.schema == 2:
                Objet_Concentration = Profil_Concentration_Centree_MMS(self.delta_r_vect[i], self.delta_t, self.N_vect[i], self.R, self.t_final, Deff, k)

            Objet_Concentration.Algorithme_Resolution()

            # # Plot
            Objet_Graphique = Plot_Concentration(Objet_Concentration.C, self.N_vect[i], 0.0, Objet_Concentration.f_T_MMS, Objet_Concentration.Delta_t)
            Objet_Graphique.Plot_MMS()
            Objet_Graphique.Plot_Numerique()

            Objet_Graphique.Save_plot("schema %d_%d"%(self.schema,self.N_vect[i]), "Comparaison de résultat schéma %d ,%d noeuds"%(self.schema, self.N_vect[i]))
            
            # Erreur
            # Objet_Norme_Erreur = Norme_Erreur_Discretisation_MMS(Objet_Concentration.f_T_MMS, Objet_Concentration.C[-1,:], Objet_Concentration.Delta_t, Objet_Concentration.R)
            Objet_Norme_Erreur = Norme_Erreur_Discretisation_MMS(Objet_Concentration.f_T_MMS, Objet_Concentration.C, self.delta_t, Objet_Concentration.R)
            erreur_vect_L1[i], erreur_vect_L2[i], erreur_vect_L_inf[i] = Objet_Norme_Erreur.Calcul_Norme()
            
            self.Erreur_matrice = Objet_Norme_Erreur.Erreur_matrice

            del Objet_Concentration
            del Objet_Graphique

        return erreur_vect_L1, erreur_vect_L2, erreur_vect_L_inf

        
class Etude_Convergence_MMS_temporel():
    def __init__(self, delta_r, delta_t_vect, N, R, t_final, schema):
        self.delta_r = delta_r
        self.delta_t_vect = delta_t_vect
        self.N_vect_t = t_final/delta_t_vect + 1
        self.N = N
        self.R = R
        self.t_final = t_final
        self.schema = schema


    def Boucle_iterations(self):
        erreur_vect_L1 = np.zeros(len(self.N_vect_t))
        erreur_vect_L2 = np.zeros(len(self.N_vect_t))
        erreur_vect_L_inf = np.zeros(len(self.N_vect_t))
        
        Deff = 1.0e-9
        k = 4.0e-9
        for i in range(len(self.N_vect_t)):

            # Resolution
            if self.schema==1:
                Objet_Concentration = Profil_Concentration_MMS(self.delta_r, self.delta_t_vect[i], self.N, self.R, self.t_final, Deff, k)
            elif self.schema == 2:
                Objet_Concentration = Profil_Concentration_Centree_MMS(self.delta_r, self.delta_t_vect[i], self.N, self.R, self.t_final, Deff, k)
            Objet_Concentration.Algorithme_Resolution()

    
            # # Plot
            Objet_Graphique = Plot_Concentration(Objet_Concentration.C, self.N, 0.0, Objet_Concentration.f_T_MMS, Objet_Concentration.Delta_t)
            Objet_Graphique.Plot_MMS()
            Objet_Graphique.Plot_Numerique()

            Objet_Graphique.Save_plot("etude temporel schema %d_%d"%(self.schema,self.delta_t_vect[i]), "Comparaison de résultat schéma %d ,dt = %d"%(self.schema, self.delta_t_vect[i]))
            
            # Erreur
            Objet_Norme_Erreur = Norme_Erreur_Discretisation_MMS(Objet_Concentration.f_T_MMS, Objet_Concentration.C, self.delta_t_vect[i], Objet_Concentration.R)
            erreur_vect_L1[i], erreur_vect_L2[i], erreur_vect_L_inf[i] = Objet_Norme_Erreur.Calcul_Norme()
            
            self.Erreur_matrice = Objet_Norme_Erreur.Erreur_matrice

            del Objet_Concentration
            del Objet_Graphique

        return erreur_vect_L1, erreur_vect_L2, erreur_vect_L_inf    
