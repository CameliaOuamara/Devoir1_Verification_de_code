
"""
Classe qui effectue une étude de convergence sur le probleme étudié
La classe calcule la solution numerique du probleme pour plusieurs discretisations geometriques
"""

# Libraries
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Functions
from profil_concentration import *
from plot import *
from norme_erreur_discretisation import *

class Etude_Convergence():
    def __init__(self, delta_r_vect, delta_t_vect, N_vect, R, critere_convergence, critere_max_iter, schema, T_final, Classe):
        self.delta_r_vect = delta_r_vect
        self.delta_t_vect = delta_t_vect
        self.N_vect = N_vect
        self.R = R
        self.critere_convergence = critere_convergence
        self.critere_max_iter = critere_max_iter
        self.schema = schema
        self.t_final = 0
        self.T_final = T_final
        self.Classe = Classe
        
        self.C_exact = 0


        
    def Boucle_iterations(self,outputFolder):
        erreur_vect_L1 = np.zeros(len(self.N_vect))
        erreur_vect_L2 = np.zeros(len(self.N_vect))
        erreur_vect_L_inf = np.zeros(len(self.N_vect))
        
        
        for i in range(len(self.N_vect)):

            # Solution exacte
            if self.Classe == 'MMS':
                # Solution mms
                Classe = MMS()
                r_vecteur = np.linspace(0,self.R, self.N_vect[i])
                self.C_exact = Classe.C(r_vecteur, self.T_final)
            else :
                data = pd.read_csv('Comsol_results_%d.txt' %self.N_vect[i],sep='\s+',header=None)
                data = pd.DataFrame(data)

                y_Comsol = data[1]
                
                self.C_exact = y_Comsol

            # Resolution
            if self.schema==1:
                Objet_Concentration = Profil_Concentration(self.delta_r_vect[i], self.delta_t_vect[i], self.N_vect[i], self.R, self.critere_convergence, self.critere_max_iter, self.T_final, self.Classe, self.schema)
            elif self.schema == 2:
                Objet_Concentration = Profil_Concentration_Centree(self.delta_r_vect[i], self.delta_t_vect[i], self.N_vect[i], self.R, self.critere_convergence, self.critere_max_iter, self.T_final, self.Classe, self.schema)
            Objet_Concentration.Algorithme_Resolution()
            self.t_final = Objet_Concentration.t

            # Plot
            Objet_Graphique = Plot_Concentration(Objet_Concentration.C, self.N_vect[i], self.T_final, self.C_exact)
            Objet_Graphique.Plot_Numerique()
            #Objet_Graphique.Plot_Exact()
            #Objet_Graphique.Save_plot("schema1_"+str(N_vect[i]), "Comparaison de résultat premier schéma, "+str(N_vect[i])+" noeuds")
            #Objet_Graphique.Save_plot(outputFolder+"schema_%d_%d"%(self.schema,self.N_vect[i]), "Comparaison de résultat schéma %d ,%d noeuds"%(self.schema, self.N_vect[i]))
            
            # Erreur
            # Objet_Norme_Erreur = Norme_Erreur_Discretisation(Objet_Graphique.C_exact, Objet_Concentration.C[-1,:])
            Objet_Norme_Erreur = Norme_Erreur_Discretisation(self.C_exact, Objet_Concentration.C[-1,:])
            erreur_vect_L1[i], erreur_vect_L2[i], erreur_vect_L_inf[i] = Objet_Norme_Erreur.Calcul_Norme()

            del Objet_Concentration
            del Objet_Graphique

        return erreur_vect_L1, erreur_vect_L2, erreur_vect_L_inf
        

class Etude_convergence_Comsol():
    def __init__(self, delta_r_vect, delta_t_vect, N_vect, R, critere_convergence, critere_max_iter, schema, t_final, Classe):
        self.delta_r_vect = delta_r_vect
        self.delta_t_vect = delta_t_vect
        self.N_vect = N_vect
        self.R = R
        self.critere_convergence = critere_convergence
        self.critere_max_iter = critere_max_iter
        self.schema = schema
        self.t_final = t_final
        self.Classe = Classe

        
        
    def Boucle_Iterations_Espace(self):
        
        # Etude de convergence du pas d'espace 
        
        # Initialisation vecteurs contenant la norme de l'erreur entre la solution
        # numerique et exact pour chaque pas d'espace
        erreur_vect_L1 = np.zeros(len(self.N_vect))
        erreur_vect_L2 = np.zeros(len(self.N_vect))
        erreur_vect_L_inf = np.zeros(len(self.N_vect))
        
        # Resolution numerique des equations
        for i in range(len(self.N_vect)):
            # La boucle for va passer a travers les differentes valeurs de delta_r
            # avec une valeur de delta_t fixée au pas de temps le plus petit
            
            Delta_t = 10**10
            
            # Solution numerique :
            # Initialisation des objets
            if self.schema==1:
                Objet_Concentration = Profil_Concentration(self.delta_r_vect[i], Delta_t, self.N_vect[i], self.R, self.critere_convergence, self.critere_max_iter, self.t_final, self.Classe, self.schema)
            elif self.schema == 2:
                Objet_Concentration = Profil_Concentration_Centree(self.delta_r_vect[i], Delta_t, self.N_vect[i], self.R, self.critere_convergence, self.critere_max_iter, self.t_final, self.Classe, self.schema)
            
            # Resolution
            Objet_Concentration.Algorithme_Resolution()
            
            # Solution exacte (Comsol):
            data = pd.read_csv('Comsol_results_%d_espace.txt'%(i+1),sep='\s+',header=None)
            data = pd.DataFrame(data)

            C_exact = data.iloc[1000:,:]
            # Calcul norme de l'erreur :
            Objet_Norme_Erreur = Norme_Erreur_Discretisation(C_exact, Objet_Concentration.C[1000:,:])
            erreur_vect_L1[i], erreur_vect_L2[i], erreur_vect_L_inf[i] = Objet_Norme_Erreur.Calcul_Norme()

            del Objet_Concentration

        return erreur_vect_L1, erreur_vect_L2, erreur_vect_L_inf
    
    def Boucle_Iterations_Temps(self):
        
        # Etude de convergence du pas de temps 
        
        # Initialisation vecteurs contenant la norme de l'erreur entre la solution
        # numerique et exact pour chaque pas de temps
        erreur_vect_L1 = np.zeros(len(self.delta_t_vect))
        erreur_vect_L2 = np.zeros(len(self.delta_t_vect))
        erreur_vect_L_inf = np.zeros(len(self.delta_t_vect))
        
        # Resolution numerique des equations
        for i in range(len(self.delta_t_vect)):
            # La boucle for va passer a travers les differentes valeurs de delta_r
            # avec une valeur de delta_r fixée au pas d'espace le plus petit
            print("--------------------------------------------")
            print(i)
            print("--------------------------------------------")

            N = 640
            Delta_r = self.R/(N-1)
            # Solution numerique :
            # Initialisation des objets
            if self.schema==1:
                Objet_Concentration = Profil_Concentration(Delta_r, self.delta_t_vect[i], N, self.R, self.critere_convergence, self.critere_max_iter, self.t_final, self.Classe, self.schema)
            elif self.schema == 2:
                Objet_Concentration = Profil_Concentration_Centree(Delta_r, self.delta_t_vect[i], N, self.R, self.critere_convergence, self.critere_max_iter, self.t_final, self.Classe, self.schema)
            
            # Resolution
            Objet_Concentration.Algorithme_Resolution()
            
            # Solution exacte (Comsol):
            data = pd.read_csv('Comsol_results_%d_temps.txt'%(i+1),sep='\s+',header=None)
            data = pd.DataFrame(data)

            C_exact = data.iloc[0:,:]
            # C_exact = data[:]
            
            # Calcul norme de l'erreur :
            Objet_Norme_Erreur = Norme_Erreur_Discretisation(C_exact, Objet_Concentration.C[0:,:])
            erreur_vect_L1[i], erreur_vect_L2[i], erreur_vect_L_inf[i] = Objet_Norme_Erreur.Calcul_Norme()

            del Objet_Concentration

        return erreur_vect_L1, erreur_vect_L2, erreur_vect_L_inf
        
        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    