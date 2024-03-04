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
#import sympy as sp
import matplotlib.pyplot as plt
import pandas as pd

# Functions
from profil_concentration import *
from plot import *
from norme_erreur_discretisation import *
from etude_convergence import *



# -----------------------------------------------------------------------------------------------------------------
#                                               Debut du code
# -----------------------------------------------------------------------------------------------------------------
R = 0.5                                 # Coefficient de diffusion effectif [m2/s]
outputFolder = "y:\DOCuments\Verification et validation\Devoir1_Verification_de_code"

# Critere de convergence
critere_convergence = 1.0e-14                        # Critere sur la valeur de la concentration a l'iteration i vs i-1
critere_max_iter = 100                               # Nombre minimum d'iterations a realiser pour la convergence vers le regime permanent
N_vect = np.arange(5,8,1, dtype=int)
N_vect = 5 * 2**N_vect
delta_r_vect = R/(N_vect-1)                          # Vecteur Delta r correspondant au vecteur N precedent

# Le pas de temps va prendre les valeurs de : [2*10**12->5*10**12]
delta_t_values = np.linspace(0.005,0.1,5)
#delta_t_values = delta_t_values*10**5

# Variation du pas de temps et d'espace. Pour chaque pas de temps, faire varier le pas d'espace
# On impose le temps final de la simulation
# t_final = 5*10**13
t_final =10

Classe_name = 'MMS'

# Solution mms
Classe = MMS()
r_vecteur = np.linspace(0, R, 100)
C_exact = Classe.C(r_vecteur, t_final)

# ----------------------------------------------------------------------------------------
#          Schéma 1 : Discretisation d'ordre 1 en temps et 1 en espace
# ----------------------------------------------------------------------------------------

# Initialisation des matrices des normes des erreurs
# Chaque ligne correspond a un pas de temps
# Chaque colonne correspond a un pas d'espace
L1_matrix = np.zeros((len(delta_t_values),len(delta_r_vect)))
L2_matrix = np.zeros((len(delta_t_values),len(delta_r_vect)))
Linf_matrix = np.zeros((len(delta_t_values),len(delta_r_vect)))

for j in range(len(delta_t_values)):
    plt.figure()
    # Le vecteur delta_t_vec est juste utilisé pour utiliser le meme delta_t pendant qu'on fait varier delta_r
    delta_t_vect = np.ones(len(delta_r_vect))*delta_t_values[j]

    # Solution numerique    
    Objet_Etude_Convergence = Etude_Convergence(delta_r_vect, delta_t_vect, N_vect, R, critere_convergence, critere_max_iter, 1, t_final, Classe_name)
    
    # Calcul de l'erreur
    erreur_vect_L1, erreur_vect_L2, erreur_vect_L_inf = Objet_Etude_Convergence.Boucle_iterations(outputFolder)
    L1_matrix[j,:] = erreur_vect_L1
    L2_matrix[j,:] = erreur_vect_L2
    Linf_matrix[j,:] = erreur_vect_L_inf
    
    plt.plot(r_vecteur, C_exact, label="MMS")
    plt.legend()
    plt.grid()
    plt.show()
    
# L1_double_integrale = np.mean(L1_matrix, axis=0)
# L2_double_integrale = np.mean(L2_matrix, axis=0)
# Linf_double_integrale = np.mean(Linf_matrix, axis=0)
L1_double_integrale = L1_matrix[0,:]
L2_double_integrale = L2_matrix[0,:]
Linf_double_integrale = Linf_matrix[0,:]

# Plot solution Comsol
# plt.plot(x_Comsol,y_Comsol, label="Solution Comsol")

# plt.plot(r_vecteur, C_MMS, label="MMS")
# plt.legend()
# plt.grid()
#plt.show()

plt.figure(1)

# Graphique log-log norme de l'erreur L1 vs delta_r
plt.loglog(delta_r_vect, L1_double_integrale, '.r', label = "Norme L1")

# Ajuster une loi de puissance à toutes les valeurs (en utilisant np.polyfit avec logarithmes)
coefficients = np.polyfit(np.log(delta_r_vect), np.log(L1_double_integrale), 1)
exponent_logreg = coefficients[0]
constant_logreg = coefficients[1]

# Fonction de régression en termes de logarithmes
fit_function_log = lambda x: exponent_logreg * x + constant_logreg

# Fonction de régression en termes originaux
fit_function = lambda x: np.exp(fit_function_log(np.log(x)))

# Extrapoler la valeur prédite pour la dernière valeur de h_values
extrapolated_value = fit_function(delta_r_vect[-1])
plt.loglog(delta_r_vect, fit_function(delta_r_vect), linestyle='--', color='r')

# Afficher l'équation de la régression en loi de puissance pour la norme L1
equation_text = f'$L_1 = {np.exp(constant_logreg):.4E} \\times Δr^{{{exponent_logreg:.4f}}}$'
equation_text_obj = plt.text(0.5, 0.05, equation_text, fontsize=12, transform=plt.gca().transAxes, color='k')

# Graphique log-log norme de l'erreur L2 vs delta_r
plt.loglog(delta_r_vect, L2_double_integrale, '.g', label = "Norme L2")

# Ajuster une loi de puissance à toutes les valeurs (en utilisant np.polyfit avec logarithmes)
coefficients = np.polyfit(np.log(delta_r_vect), np.log(L2_double_integrale), 1)
exponent_logreg = coefficients[0]
constant_logreg = coefficients[1]

# Fonction de régression en termes de logarithmes
fit_function_log = lambda x: exponent_logreg * x + constant_logreg

# Fonction de régression en termes originaux
fit_function = lambda x: np.exp(fit_function_log(np.log(x)))

# Extrapoler la valeur prédite pour la dernière valeur de h_values
extrapolated_value = fit_function(delta_r_vect[-1])
plt.loglog(delta_r_vect, fit_function(delta_r_vect), linestyle='--', color='g')

# Afficher l'équation de la régression en loi de puissance pour la norme L2
equation_text = f'$L_2 = {np.exp(constant_logreg):.4E} \\times Δr^{{{exponent_logreg:.4f}}}$'
equation_text_obj = plt.text(0.5, 0.15, equation_text, fontsize=12, transform=plt.gca().transAxes, color='k')

# Graphique log-log norme de l'erreur Linf vs delta_r
plt.loglog(delta_r_vect, Linf_double_integrale, '.m', label='Norme $L_\infty$')

# Ajuster une loi de puissance à toutes les valeurs (en utilisant np.polyfit avec logarithmes)
coefficients = np.polyfit(np.log(delta_r_vect), np.log(Linf_double_integrale), 1)
exponent_logreg = coefficients[0]
constant_logreg = coefficients[1]

# Fonction de régression en termes de logarithmes
fit_function_log = lambda x: exponent_logreg * x + constant_logreg

# Fonction de régression en termes originaux
fit_function = lambda x: np.exp(fit_function_log(np.log(x)))

# Extrapoler la valeur prédite pour la dernière valeur de h_values
extrapolated_value = fit_function(delta_r_vect[-1])
plt.loglog(delta_r_vect, fit_function(delta_r_vect), linestyle='--', color='m')

# Afficher l'équation de la régression en loi de puissance pour la norme Linf
equation_text = f'$L_\infty = {np.exp(constant_logreg):.4E} \\times Δr^{{{exponent_logreg:.4f}}}$'
equation_text_obj = plt.text(0.5, 0.25, equation_text, fontsize=12, transform=plt.gca().transAxes, color='k')

plt.xlabel("$\Delta r$")
plt.ylabel("Norme de l'erreur")
plt.legend()
plt.grid()
plt.title("Normes des erreurs L1, L2 et $L_\infty$ schéma d'ordre 1 en fonction de $\Delta r$")
#plt.savefig(outputFolder+"Norme_des_erreurs_Schema_1.png")
plt.show()


#%% ----------------------------------------------------------------------------------------
#          Schéma 2 : Discretisation d'ordre 1 en temps et 2 en espace
# ----------------------------------------------------------------------------------------
erreur_vect_L1_Centree = np.zeros(len(N_vect))
erreur_vect_L2_Centree = np.zeros(len(N_vect))
erreur_vect_L_inf_Centree = np.zeros(len(N_vect))

L1_matrix = np.zeros((len(delta_t_values),len(delta_r_vect)))
L2_matrix = np.zeros((len(delta_t_values),len(delta_r_vect)))
Linf_matrix = np.zeros((len(delta_t_values),len(delta_r_vect)))

plt.figure(2)

for j in range(len(delta_t_values)):
    # Le vecteur delta_t_vec est juste utilisé pour utiliser le meme delta_t pendant qu'on fait varier delta_r
    delta_t_vect = np.ones(len(delta_r_vect))*delta_t_values[j]

    # Solution numerique    
    Objet_Etude_Convergence = Etude_Convergence(delta_r_vect, delta_t_vect, N_vect, R, critere_convergence, critere_max_iter, 2, t_final, Classe_name)
    
    # Calcul de l'erreur
    erreur_vect_L1_Centree, erreur_vect_L2_Centree, erreur_vect_L_inf_Centree = Objet_Etude_Convergence.Boucle_iterations(outputFolder)
    L1_matrix[j,:] = erreur_vect_L1_Centree
    L2_matrix[j,:] = erreur_vect_L2_Centree
    Linf_matrix[j,:] = erreur_vect_L_inf_Centree
plt.show()
 
L1_double_integrale = np.mean(L1_matrix, axis=0)
L2_double_integrale = np.mean(L2_matrix, axis=0)
Linf_double_integrale = np.mean(Linf_matrix, axis=0)

plt.figure(3)

# Graphique log-log norme de l'erreur L1 vs delta_r
plt.loglog(delta_r_vect, L1_double_integrale, '.r', label = "Norme L1")

# Ajuster une loi de puissance à toutes les valeurs (en utilisant np.polyfit avec logarithmes)
coefficients = np.polyfit(np.log(delta_r_vect), np.log(L1_double_integrale), 1)
exponent_logreg = coefficients[0]
constant_logreg = coefficients[1]

# Fonction de régression en termes de logarithmes
fit_function_log = lambda x: exponent_logreg * x + constant_logreg

# Fonction de régression en termes originaux
fit_function = lambda x: np.exp(fit_function_log(np.log(x)))

# Extrapoler la valeur prédite pour la dernière valeur de h_values
extrapolated_value = fit_function(delta_r_vect[-1])
plt.loglog(delta_r_vect, fit_function(delta_r_vect), linestyle='--', color='r')

# Afficher l'équation de la régression en loi de puissance pour la norme L1
equation_text = f'$L_1 = {np.exp(constant_logreg):.4E} \\times Δr^{{{exponent_logreg:.4f}}}$'
equation_text_obj = plt.text(0.05, 0.05, equation_text, fontsize=12, transform=plt.gca().transAxes, color='k')

# Graphique log-log norme de l'erreur L2 vs delta_r
plt.loglog(delta_r_vect, L2_double_integrale, '.g', label = "Nomre L2")

# Ajuster une loi de puissance à toutes les valeurs (en utilisant np.polyfit avec logarithmes)
coefficients = np.polyfit(np.log(delta_r_vect), np.log(L2_double_integrale), 1)
exponent_logreg = coefficients[0]
constant_logreg = coefficients[1]

# Fonction de régression en termes de logarithmes
fit_function_log = lambda x: exponent_logreg * x + constant_logreg

# Fonction de régression en termes originaux
fit_function = lambda x: np.exp(fit_function_log(np.log(x)))

# Extrapoler la valeur prédite pour la dernière valeur de h_values
extrapolated_value = fit_function(delta_r_vect[-1])
plt.loglog(delta_r_vect, fit_function(delta_r_vect), linestyle='--', color='g')

# Afficher l'équation de la régression en loi de puissance pour la norme L2
equation_text = f'$L_2 = {np.exp(constant_logreg):.4E} \\times Δr^{{{exponent_logreg:.4f}}}$'
equation_text_obj = plt.text(0.05, 0.15, equation_text, fontsize=12, transform=plt.gca().transAxes, color='k')

# Graphique log-log norme de l'erreur Linf vs delta_r
plt.loglog(delta_r_vect, Linf_double_integrale, '.m', label = "Norme $L_\infty$")

# Ajuster une loi de puissance à toutes les valeurs (en utilisant np.polyfit avec logarithmes)
coefficients = np.polyfit(np.log(delta_r_vect), np.log(Linf_double_integrale), 1)
exponent_logreg = coefficients[0]
constant_logreg = coefficients[1]

# Fonction de régression en termes de logarithmes
fit_function_log = lambda x: exponent_logreg * x + constant_logreg

# Fonction de régression en termes originaux
fit_function = lambda x: np.exp(fit_function_log(np.log(x)))

# Extrapoler la valeur prédite pour la dernière valeur de h_values
extrapolated_value = fit_function(delta_r_vect[-1])
plt.loglog(delta_r_vect, fit_function(delta_r_vect), linestyle='--', color='m')

# Afficher l'équation de la régression en loi de puissance pour la norme Linf
equation_text = f'$Linf = {np.exp(constant_logreg):.4E} \\times Δr^{{{exponent_logreg:.4f}}}$'
equation_text_obj = plt.text(0.05, 0.25, equation_text, fontsize=12, transform=plt.gca().transAxes, color='k')

plt.xlabel("$\Delta r$")
plt.ylabel("erreur L_inf")
plt.legend()
plt.grid()
plt.title("Normes des erreurs L1, L2 et $L_\infty$ schéma d'ordre 2 en fonction de $\Delta r$")
#plt.savefig(outputFolder+"Norme_des_erreurs_Schema_2.png")
plt.show()








