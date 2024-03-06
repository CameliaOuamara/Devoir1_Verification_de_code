# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 12:33:38 2024

@author: ouama
"""

# Libraries
import numpy as np
#import sympy as sp
import matplotlib.pyplot as plt
import pandas as pd

# Functions
from profil_concentration import *
from plot import *
from norme_erreur_discretisation import *
from etude_convergence import *
from regression_loi_puissance import *

# -----------------------------------------------------------------------------------------------------------------
#                                               Definition parametres de l'étude
# -----------------------------------------------------------------------------------------------------------------
R = 0.5                                 # Coefficient de diffusion effectif [m2/s]
critere_convergence = 1.0e-14 
critere_max_iter = 100

# Pas d'espace
N_vect = np.arange(4,8,1, dtype=int)
N_vect = 5 * 2**N_vect
delta_r_vect = R/(N_vect-1)

# Pas de temps
delta_t_vect = np.linspace(11,9, 3)
delta_t_vect = 10**delta_t_vect

# Imposition du temps final de la simulation
t_final = 5*10**13

Classe = "Fick"

# -----------------------------------------------------------------------------------------------------------------
#                                         Schema 1 : Ordre 1 en temps et 1 en espace
# -----------------------------------------------------------------------------------------------------------------
schema = 1

# Initialisation de l'objet
Objet_Etude_Convergence = Etude_convergence_Comsol(delta_r_vect, delta_t_vect, N_vect, R, critere_convergence, critere_max_iter, schema, t_final, Classe)

# -----------------------------------------------------------------------------------------------------------------
#                                         Étude de convergence sur le pas d'espace
# -----------------------------------------------------------------------------------------------------------------
# # Calcul de l'erreur
# erreur_vect_L1, erreur_vect_L2, erreur_vect_L_inf = Objet_Etude_Convergence.Boucle_Iterations_Espace()

# # Regression en loi de puissance
# Objet_regression = Regression_Loi_Puissance(erreur_vect_L1, erreur_vect_L2, erreur_vect_L_inf, schema, delta_r_vect, delta_t_vect)
# Objet_regression.Plot_Regression_Espace()
# -----------------------------------------------------------------------------------------------------------------
#                                         Étude de convergence sur le pas de temps
# -----------------------------------------------------------------------------------------------------------------
# Calcul de l'erreur
erreur_vect_L1, erreur_vect_L2, erreur_vect_L_inf = Objet_Etude_Convergence.Boucle_Iterations_Temps()

#Regression en loi de puissance
Objet_regression = Regression_Loi_Puissance(erreur_vect_L1, erreur_vect_L2, erreur_vect_L_inf, schema, delta_r_vect, delta_t_vect)
Objet_regression.Plot_Regression_Temps()

#%% -----------------------------------------------------------------------------------------------------------------
#                                         Schema 2 : Ordre 1 en temps et 2 en espace
# -----------------------------------------------------------------------------------------------------------------
schema = 2

# Initialisation de l'objet
Objet_Etude_Convergence = Etude_convergence_Comsol(delta_r_vect, delta_t_vect, N_vect, R, critere_convergence, critere_max_iter, schema, t_final, Classe)

# -----------------------------------------------------------------------------------------------------------------
#                                         Étude de convergence sur le pas d'espace
# -----------------------------------------------------------------------------------------------------------------
# Calcul de l'erreur
# erreur_vect_L1, erreur_vect_L2, erreur_vect_L_inf = Objet_Etude_Convergence.Boucle_Iterations_Espace()

# # Regression en loi de puissance
# Objet_regression = Regression_Loi_Puissance(erreur_vect_L1, erreur_vect_L2, erreur_vect_L_inf, schema, delta_r_vect, delta_t_vect)
# Objet_regression.Plot_Regression_Espace()
# -----------------------------------------------------------------------------------------------------------------
#                                         Étude de convergence sur le pas de temps
# -----------------------------------------------------------------------------------------------------------------
# Calcul de l'erreur
erreur_vect_L1, erreur_vect_L2, erreur_vect_L_inf = Objet_Etude_Convergence.Boucle_Iterations_Temps()

# Regression en loi de puissance
Objet_regression = Regression_Loi_Puissance(erreur_vect_L1, erreur_vect_L2, erreur_vect_L_inf, schema, delta_r_vect, delta_t_vect)
Objet_regression.Plot_Regression_Temps()



































