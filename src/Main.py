"""
------------------------------------------------------------------------------------------------------------------
MEC-8811 Vérification et Validation en Modélisation Numrérique
Devoir 1
------------------------------------------------------------------------------------------------------------------

Ce code utilise la methode des differences finies dans le but de resoudre un probleme de diffusion.
Le probleme étudié ici est le processus de diffusion du sel dans un pilier de béton poreux.
L'équation à résoudre est l'équation de Fick.
Le code resout le cas transitoire et permet d'atteindre le regime permanent.

"""

# Libraries
import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
import pandas as pd

# Functions
from profil_concentration import *
from plot import *
from norme_erreur_discretisation import *
from etude_convergence import *

# -----------------------------------------------------------------------------------------------------------------
#                                 Extraction solution COMSOL
# -----------------------------------------------------------------------------------------------------------------

# Solution Comsol

data = pd.read_csv('Comsol_results.txt',sep='\s+',header=None)
data = pd.DataFrame(data)

x_Comsol = data[0]
y_Comsol = data[1]

# -----------------------------------------------------------------------------------------------------------------
#                                               Debut du code
# -----------------------------------------------------------------------------------------------------------------

# Données 
R = 0.5                                              # Rayon du cylindre [m]
N = 100                                              # Nombre de points de discretisation [-]
Ce = 12
delta_r = R/(N-1)                                    # Taille de l'intervalle geometrique [m]
D_eff = 1.0e-10                                      # Coefficient de diffusion effectif [m2/s]
outputFolder = "y:\DOCuments\Verification et validation\Devoir1_Verification_de_code"

# Solution MMS
Terme_Source = 0
CL_R = Ce

# Règle de pouce : Prendre un dt à la limite de la stabilité pour le schéma explicite
delta_t = 0.5 * delta_r*delta_r / D_eff              # Pas de temps [s]

# Critere de convergence
critere_convergence = 1.0e-14                        # Critere sur la valeur de la concentration a l'iteration i vs i-1
critere_max_iter = 100                               # Nombre minimum d'iterations a realiser pour la convergence vers le regime permanent

# Étude convergence
#N_vect = np.arange(0,7,1, dtype=int)                 # Vecteur contenant les N utilisés dans l'étude de convergence
N_vect = np.arange(9,10,1, dtype=int)
N_vect = 5 * 2**N_vect
delta_r_vect = R/(N_vect-1)                          # Vecteur Delta r correspondant au vecteur N precedent
delta_t_vect = 1.0e5 * 0.5 * delta_r_vect*delta_r_vect / D_eff # Vecteur Delta t correspondant au vecteur Delta r precedent

# ----------------------------------------------------------------------------------------
#          Schéma 1 : Discretisation d'ordre 1 en temps et 1 en espace
# ----------------------------------------------------------------------------------------

# Solution numerique

plt.figure(0)

Objet_Etude_Convergence = Etude_Convergence(delta_r_vect, delta_t_vect, N_vect, R, critere_convergence, critere_max_iter, 1)
erreur_vect_L1, erreur_vect_L2, erreur_vect_L_inf = Objet_Etude_Convergence.Boucle_iterations(outputFolder)

# Plot solution Comsol

#plt.plot(x_Comsol,y_Comsol, label="Solution Comsol")
# Solution MMS (Temporaire)
Ce = 12
R = 0.5
D = 10**-10
k = 4*10**-9
r, t = sy.symbols('r t')

#%% Solution mms
C_MMS = sy.exp(sy.pi*r/R)*sy.sin(t)
C =  sy.lambdify([r, t], C_MMS, "numpy")
r = np.linspace(0, R, 5*2**5)
C_vec = np.zeros(len(r))
for i in range(len(C_vec)):
    C_vec[i] = C(r[i], 3*10**14)
plt.plot(r, C_vec, label="MMS")
plt.legend()
plt.grid()
plt.show()

# plt.figure(1)

# # Graphique log-log norme de l'erreur L1 vs delta_r
# plt.loglog(delta_r_vect, erreur_vect_L1, '.r', label = "Norme L1")

# # Ajuster une loi de puissance à toutes les valeurs (en utilisant np.polyfit avec logarithmes)
# coefficients = np.polyfit(np.log(delta_r_vect), np.log(erreur_vect_L1), 1)
# exponent_logreg = coefficients[0]
# constant_logreg = coefficients[1]

# # Fonction de régression en termes de logarithmes
# fit_function_log = lambda x: exponent_logreg * x + constant_logreg

# # Fonction de régression en termes originaux
# fit_function = lambda x: np.exp(fit_function_log(np.log(x)))

# # Extrapoler la valeur prédite pour la dernière valeur de h_values
# extrapolated_value = fit_function(delta_r_vect[-1])
# plt.loglog(delta_r_vect, fit_function(delta_r_vect), linestyle='--', color='r')

# # Afficher l'équation de la régression en loi de puissance pour la norme L1
# equation_text = f'$L_1 = {np.exp(constant_logreg):.4E} \\times Δr^{{{exponent_logreg:.4f}}}$'
# equation_text_obj = plt.text(0.5, 0.05, equation_text, fontsize=12, transform=plt.gca().transAxes, color='k')

# # Graphique log-log norme de l'erreur L2 vs delta_r
# plt.loglog(delta_r_vect, erreur_vect_L2, '.g', label = "Norme L2")

# # Ajuster une loi de puissance à toutes les valeurs (en utilisant np.polyfit avec logarithmes)
# coefficients = np.polyfit(np.log(delta_r_vect), np.log(erreur_vect_L2), 1)
# exponent_logreg = coefficients[0]
# constant_logreg = coefficients[1]

# # Fonction de régression en termes de logarithmes
# fit_function_log = lambda x: exponent_logreg * x + constant_logreg

# # Fonction de régression en termes originaux
# fit_function = lambda x: np.exp(fit_function_log(np.log(x)))

# # Extrapoler la valeur prédite pour la dernière valeur de h_values
# extrapolated_value = fit_function(delta_r_vect[-1])
# plt.loglog(delta_r_vect, fit_function(delta_r_vect), linestyle='--', color='g')

# # Afficher l'équation de la régression en loi de puissance pour la norme L2
# equation_text = f'$L_2 = {np.exp(constant_logreg):.4E} \\times Δr^{{{exponent_logreg:.4f}}}$'
# equation_text_obj = plt.text(0.5, 0.15, equation_text, fontsize=12, transform=plt.gca().transAxes, color='k')

# # Graphique log-log norme de l'erreur Linf vs delta_r
# plt.loglog(delta_r_vect, erreur_vect_L_inf, '.m', label='Norme $L_\infty$')

# # Ajuster une loi de puissance à toutes les valeurs (en utilisant np.polyfit avec logarithmes)
# coefficients = np.polyfit(np.log(delta_r_vect), np.log(erreur_vect_L_inf), 1)
# exponent_logreg = coefficients[0]
# constant_logreg = coefficients[1]

# # Fonction de régression en termes de logarithmes
# fit_function_log = lambda x: exponent_logreg * x + constant_logreg

# # Fonction de régression en termes originaux
# fit_function = lambda x: np.exp(fit_function_log(np.log(x)))

# # Extrapoler la valeur prédite pour la dernière valeur de h_values
# extrapolated_value = fit_function(delta_r_vect[-1])
# plt.loglog(delta_r_vect, fit_function(delta_r_vect), linestyle='--', color='m')

# # Afficher l'équation de la régression en loi de puissance pour la norme Linf
# equation_text = f'$L_\infty = {np.exp(constant_logreg):.4E} \\times Δr^{{{exponent_logreg:.4f}}}$'
# equation_text_obj = plt.text(0.5, 0.25, equation_text, fontsize=12, transform=plt.gca().transAxes, color='k')

# plt.xlabel("delta_r")
# plt.ylabel("Norme de l'erreur")
# plt.legend()
# plt.grid()
# plt.title("Normes des erreurs L1, L2 et $L_\infty$ schéma d'ordre 1 en fonction de N")
# plt.savefig(outputFolder+"Norme_des_erreurs_Schema_1.png")
# plt.show()

#%%
# ----------------------------------------------------------------------------------------
#          Schéma 2 : Discretisation d'ordre 1 en temps et 2 en espace
# ----------------------------------------------------------------------------------------
# erreur_vect_L1_Centree = np.zeros(len(N_vect))
# erreur_vect_L2_Centree = np.zeros(len(N_vect))
# erreur_vect_L_inf_Centree = np.zeros(len(N_vect))

# plt.figure(2)
# Objet_Etude_Convergence = Etude_Convergence(delta_r_vect, delta_t_vect, N_vect, R, critere_convergence, critere_max_iter, 2)
# erreur_vect_L1_Centree, erreur_vect_L2_Centree, erreur_vect_L_inf_Centree = Objet_Etude_Convergence.Boucle_iterations(outputFolder)

# # Plot solution Comsol

# plt.plot(x_Comsol,y_Comsol, label="Solution Comsol")
# plt.legend()
# plt.grid()
# plt.show()

# for i in range(len(N_vect)):
#     # print("i: ", i)
#     # Resolution
#     Objet_Concentration = Profil_Concentration_Centree(delta_r_vect[i], delta_t_vect[i], N_vect[i], R, critere_convergence, critere_max_iter)
#     Objet_Concentration.Algorithme_Resolution()

#     # Plot
#     Objet_Graphique = Plot_Concentration(Objet_Concentration.C, N_vect[i])
#     Objet_Graphique.Plot_Numerique()
#     Objet_Graphique.Plot_Exact()
#     Objet_Graphique.Save_plot("schema2_"+str(N_vect[i]), "Comparaison de résultat deuxième schéma, "+str(N_vect[i])+" noeuds")
    
#     # Erreur
#     Objet_Norme_Erreur = Norme_Erreur_Discretisation(Objet_Graphique.C_exact, Objet_Concentration.C[-1,:])
#     erreur_vect_L1_Centree[i], erreur_vect_L2_Centree[i], erreur_vect_L_inf_Centree[i] = Objet_Norme_Erreur.Calcul_Norme()

#     del Objet_Concentration
#     del Objet_Graphique
    

# plt.figure(3)

# # Graphique log-log norme de l'erreur L1 vs delta_r
# plt.loglog(delta_r_vect, erreur_vect_L1_Centree, '.r', label = "Norme L1")

# # Ajuster une loi de puissance à toutes les valeurs (en utilisant np.polyfit avec logarithmes)
# coefficients = np.polyfit(np.log(delta_r_vect), np.log(erreur_vect_L1_Centree), 1)
# exponent_logreg = coefficients[0]
# constant_logreg = coefficients[1]

# # Fonction de régression en termes de logarithmes
# fit_function_log = lambda x: exponent_logreg * x + constant_logreg

# # Fonction de régression en termes originaux
# fit_function = lambda x: np.exp(fit_function_log(np.log(x)))

# # Extrapoler la valeur prédite pour la dernière valeur de h_values
# extrapolated_value = fit_function(delta_r_vect[-1])
# plt.loglog(delta_r_vect, fit_function(delta_r_vect), linestyle='--', color='r')

# # Afficher l'équation de la régression en loi de puissance pour la norme L1
# equation_text = f'$L_1 = {np.exp(constant_logreg):.4E} \\times Δr^{{{exponent_logreg:.4f}}}$'
# equation_text_obj = plt.text(0.05, 0.05, equation_text, fontsize=12, transform=plt.gca().transAxes, color='k')

# # Graphique log-log norme de l'erreur L2 vs delta_r
# plt.loglog(delta_r_vect, erreur_vect_L2_Centree, '.g', label = "Nomre L2")

# # Ajuster une loi de puissance à toutes les valeurs (en utilisant np.polyfit avec logarithmes)
# coefficients = np.polyfit(np.log(delta_r_vect), np.log(erreur_vect_L2_Centree), 1)
# exponent_logreg = coefficients[0]
# constant_logreg = coefficients[1]

# # Fonction de régression en termes de logarithmes
# fit_function_log = lambda x: exponent_logreg * x + constant_logreg

# # Fonction de régression en termes originaux
# fit_function = lambda x: np.exp(fit_function_log(np.log(x)))

# # Extrapoler la valeur prédite pour la dernière valeur de h_values
# extrapolated_value = fit_function(delta_r_vect[-1])
# plt.loglog(delta_r_vect, fit_function(delta_r_vect), linestyle='--', color='g')

# # Afficher l'équation de la régression en loi de puissance pour la norme L2
# equation_text = f'$L_2 = {np.exp(constant_logreg):.4E} \\times Δr^{{{exponent_logreg:.4f}}}$'
# equation_text_obj = plt.text(0.05, 0.15, equation_text, fontsize=12, transform=plt.gca().transAxes, color='k')

# # Graphique log-log norme de l'erreur Linf vs delta_r
# plt.loglog(delta_r_vect, erreur_vect_L_inf_Centree, '.m', label = "Norme $L_\infty$")

# # Ajuster une loi de puissance à toutes les valeurs (en utilisant np.polyfit avec logarithmes)
# coefficients = np.polyfit(np.log(delta_r_vect), np.log(erreur_vect_L_inf_Centree), 1)
# exponent_logreg = coefficients[0]
# constant_logreg = coefficients[1]

# # Fonction de régression en termes de logarithmes
# fit_function_log = lambda x: exponent_logreg * x + constant_logreg

# # Fonction de régression en termes originaux
# fit_function = lambda x: np.exp(fit_function_log(np.log(x)))

# # Extrapoler la valeur prédite pour la dernière valeur de h_values
# extrapolated_value = fit_function(delta_r_vect[-1])
# plt.loglog(delta_r_vect, fit_function(delta_r_vect), linestyle='--', color='m')

# # Afficher l'équation de la régression en loi de puissance pour la norme Linf
# equation_text = f'$Linf = {np.exp(constant_logreg):.4E} \\times Δr^{{{exponent_logreg:.4f}}}$'
# equation_text_obj = plt.text(0.05, 0.25, equation_text, fontsize=12, transform=plt.gca().transAxes, color='k')

# plt.xlabel("delta_r")
# plt.ylabel("erreur L_inf")
# plt.legend()
# plt.grid()
# plt.title("Normes des erreurs L1, L2 et $L_\infty$ schéma d'ordre 2 en fonction de N")
# plt.savefig(outputFolder+"Norme_des_erreurs_Schema_2.png")
# plt.show()








