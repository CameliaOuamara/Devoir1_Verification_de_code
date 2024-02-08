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

critere_conv = 1.0e-30
critere_max_iter = 30

# # Resolution
# Objet_Concentration = Profil_Concentration(delta_r, delta_t, N, R, critere_conv)
# Objet_Concentration.Algorithme_Resolution()

# # Plot
# Objet_Graphique = Plot_Concentration( Objet_Concentration.C, N)
# Objet_Graphique.Plot_Numerique()
# Objet_Graphique.Plot_Exact()

# Étude convergence
N_vect = np.arange(0,5,1, dtype=int)

N_vect = 5 * 2**N_vect


delta_r_vect = R/(N_vect-1)

delta_t_vect = 1.0e10 * 0.5 * delta_r_vect*delta_r_vect / D_eff

erreur_vect_L1 = np.zeros(len(N_vect))
erreur_vect_L2 = np.zeros(len(N_vect))
erreur_vect_L_inf = np.zeros(len(N_vect))

plt.figure(0)
for i in range(len(N_vect)):
    # print("i: ", i)
    # Resolution
    Objet_Concentration = Profil_Concentration(delta_r_vect[i], delta_t_vect[i], N_vect[i], R, critere_conv, critere_max_iter)
    Objet_Concentration.Algorithme_Resolution()

    # Plot
    Objet_Graphique = Plot_Concentration(Objet_Concentration.C, N_vect[i])
    Objet_Graphique.Plot_Numerique()
    Objet_Graphique.Plot_Exact()
    
    #Erreur
    erreur_vect_L1[i] = np.mean(abs(Objet_Concentration.C[-1,:] - Objet_Graphique.C_exact))
    erreur_vect_L2[i] = np.sqrt(np.mean((Objet_Concentration.C[-1,:] - Objet_Graphique.C_exact)*(Objet_Concentration.C[-1,:] - Objet_Graphique.C_exact)))
    erreur_vect_L_inf[i] = np.max(abs(Objet_Concentration.C[-1,:] - Objet_Graphique.C_exact))

    del Objet_Concentration
    del Objet_Graphique
    


delta_r_plot = delta_r_vect[-1] - delta_r_vect[0]
ordre = 1.0
first_pt_plot_theo_L1 = -10.0**ordre * delta_r_plot + erreur_vect_L1[-1]
first_pt_plot_theo_L2 = -10.0**ordre * delta_r_plot + erreur_vect_L2[-1]
first_pt_plot_theo_L_inf = -10.0**ordre * delta_r_plot + erreur_vect_L_inf[-1]

plt.figure(1)
plt.loglog(delta_r_vect, erreur_vect_L1, '.r')
plt.loglog([delta_r_vect[0], delta_r_vect[-1]], [first_pt_plot_theo_L1, erreur_vect_L1[-1]])
plt.xlabel("delta_r")
plt.ylabel("erreur_L1")
plt.title("Premier schéma")

plt.figure(2)
plt.loglog(delta_r_vect, erreur_vect_L2, '.r')
plt.loglog([delta_r_vect[0], delta_r_vect[-1]], [first_pt_plot_theo_L2, erreur_vect_L2[-1]])
plt.xlabel("delta_r")
plt.ylabel("erreur_L2")
plt.title("Premier schéma")

plt.figure(3)
plt.loglog(delta_r_vect, erreur_vect_L_inf, '.r')
plt.loglog([delta_r_vect[0], delta_r_vect[-1]], [first_pt_plot_theo_L_inf, erreur_vect_L_inf[-1]])
plt.xlabel("delta_r")
plt.ylabel("erreur L_inf")
plt.title("Premier schéma")



erreur_vect_L1_Centree = np.zeros(len(N_vect))
erreur_vect_L2_Centree = np.zeros(len(N_vect))
erreur_vect_L_inf_Centree = np.zeros(len(N_vect))

plt.figure(4)
for i in range(len(N_vect)):
    # print("i: ", i)
    # Resolution
    Objet_Concentration = Profil_Concentration_Centree(delta_r_vect[i], delta_t_vect[i], N_vect[i], R, critere_conv, critere_max_iter)
    Objet_Concentration.Algorithme_Resolution()

    # Plot
    Objet_Graphique = Plot_Concentration(Objet_Concentration.C, N_vect[i])
    Objet_Graphique.Plot_Numerique()
    Objet_Graphique.Plot_Exact()
    
    #Erreur
    erreur_vect_L1_Centree[i] = np.mean(abs(Objet_Concentration.C[-1,:] - Objet_Graphique.C_exact))
    erreur_vect_L2_Centree[i] = np.sqrt(np.mean((Objet_Concentration.C[-1,:] - Objet_Graphique.C_exact)*(Objet_Concentration.C[-1,:] - Objet_Graphique.C_exact)))
    erreur_vect_L_inf_Centree[i] = np.max(abs(Objet_Concentration.C[-1,:] - Objet_Graphique.C_exact))

    del Objet_Concentration
    del Objet_Graphique
    

ordre_centree = 2.0
first_pt_plot_theo_L1_centree = -10.0**ordre_centree * delta_r_plot + erreur_vect_L1_Centree[-1]
first_pt_plot_theo_L2_centree = -10.0**ordre_centree * delta_r_plot + erreur_vect_L2_Centree[-1]
first_pt_plot_theo_L_inf_centree = -10.0**ordre_centree * delta_r_plot + erreur_vect_L_inf_Centree[-1]

plt.figure(5)
plt.loglog(delta_r_vect, erreur_vect_L1_Centree, '.r')
plt.loglog([delta_r_vect[0], delta_r_vect[-1]], [first_pt_plot_theo_L1_centree, erreur_vect_L1_Centree[-1]])
plt.xlabel("delta_r")
plt.ylabel("erreur_L1")
plt.title("Deuxième schéma")

plt.figure(6)
plt.loglog(delta_r_vect, erreur_vect_L2_Centree, '.r')
plt.loglog([delta_r_vect[0], delta_r_vect[-1]], [first_pt_plot_theo_L2_centree, erreur_vect_L2_Centree[-1]])
plt.xlabel("delta_r")
plt.ylabel("erreur_L2")
plt.title("Deuxième schéma")

plt.figure(7)
plt.loglog(delta_r_vect, erreur_vect_L_inf_Centree, '.r')
plt.loglog([delta_r_vect[0], delta_r_vect[-1]], [first_pt_plot_theo_L_inf_centree, erreur_vect_L_inf_Centree[-1]])
plt.xlabel("delta_r")
plt.ylabel("erreur L_inf")
plt.title("Deuxième schéma")












# # Debut du code
# R = 0.5
# N = 100
# delta_r = R/(N-1)
# # t_final = 1.0e10
# # N_t = 16
# # delta_t = t_final/(N_t-1)

# #bonne règle de pouce est de prendre un dt à la limite de la stabilité pour le schéma explicite
# D_eff = 1.0e-10
# delta_t = 0.5 * delta_r*delta_r / D_eff

# critere_conv = 1.0e-7

# # Resolution
# Objet_Concentration = Profil_Concentration_Centree(delta_r, delta_t, N, R, critere_conv)
# Objet_Concentration.Algorithme_Resolution()

# # Plot
# Objet_Graphique = Plot_Concentration( Objet_Concentration.C, N)
# Objet_Graphique.Plot_Numerique()
# Objet_Graphique.Plot_Exact()



