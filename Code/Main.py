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
N = 5
delta_r = R/(N-1)
t_final = 30 
N_t = 16
delta_t = t_final/(N_t-1)

# Resolution
Objet_Concentration = Profil_Concentration(delta_r, delta_t, t_final, N, R)
Objet_Concentration.Algorithme_Resolution()

# Plot
Objet_Graphique = Plot_Concentration( Objet_Concentration.C, N)
Objet_Graphique.Plot_Numerique()
Objet_Graphique.Plot_Exact()


    





















































