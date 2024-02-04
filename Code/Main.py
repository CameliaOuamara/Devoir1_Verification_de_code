"""
Programme de difference finies pour le devoir 1 du cours MEC8811. A completer
"""

# Libraries
import numpy as np
import sympy as sp
import matplotlib.pyplot as plt

# Functions
from profil_concentration import *
from plot import *



# Debut du code
R = 0.5
N = 100
delta_r = R/(N-1)
# t_final = 1.0e10
# N_t = 16
# delta_t = t_final/(N_t-1)

#bonne règle de pouce est de prendre un dt à la limite de la stabilité pour le schéma explicite
D_eff = 1.0e-10
delta_t = 0.5 * delta_r*delta_r / D_eff

critere_conv = 1.0e-7

# Resolution
Objet_Concentration = Profil_Concentration(delta_r, delta_t, N, R, critere_conv)
Objet_Concentration.Algorithme_Resolution()

# Plot
Objet_Graphique = Plot_Concentration( Objet_Concentration.C, N)
Objet_Graphique.Plot_Numerique()
Objet_Graphique.Plot_Exact()

# Étude convergence
N_vect = np.arange(0.0,8.1,1.0)

delta_r_vect = 5.0 * 2**N_vect




    





















































