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

# Functions
from profil_concentration import *
from plot import *
from norme_erreur_discretisation import *
from etude_convergence import *

# -----------------------------------------------------------------------------------------------------------------
#                                               Debut du code
# -----------------------------------------------------------------------------------------------------------------

# Données 
R = 0.5                                              # Rayon du cylindre [m]
N = 100                                              # Nombre de points de discretisation [-]
delta_r = R/(N-1)                                    # Taille de l'intervalle geometrique [m]
D_eff = 1.0e-10                                      # Coefficient de diffusion effectif [m2/s]

# Règle de pouce : Prendre un dt à la limite de la stabilité pour le schéma explicite
delta_t = 0.5 * delta_r*delta_r / D_eff              # Pas de temps [s]

# Critere de convergence
critere_convergence = 1.0e-14                        # Critere sur la valeur de la concentration a l'iteration i vs i-1
critere_max_iter = 100000                               # Nombre minimum d'iterations a realiser pour la convergence vers le regime permanent

# Étude convergence
N_vect = np.arange(0,7,1, dtype=int)                 # Vecteur contenant les N utilisés dans l'étude de convergence
N_vect = 5 * 2**N_vect
delta_r_vect = R/(N_vect-1)                          # Vecteur Delta r correspondant au vecteur N precedent
delta_t_vect = 1.0e5 * 0.5 * delta_r_vect*delta_r_vect / D_eff # Vecteur Delta t correspondant au vecteur Delta r precedent
# delta_t_vect = 0.5 * delta_r_vect*delta_r_vect / D_eff # Vecteur Delta t correspondant au vecteur Delta r precedent


# -----------------------------------------------------------------------------------------------------------------
#                           Solution avec le maillage le plus fin (pour MNP)
# -----------------------------------------------------------------------------------------------------------------
sol_MNP = Profil_Concentration_Centree(delta_r_vect[-1], delta_t_vect[-1], N_vect[-1], R, critere_convergence, critere_max_iter)
sol_MNP.Algorithme_Resolution()
nb_time_step = sol_MNP.C.shape[0]
r_fine_mesh = np.linspace(0, R, N_vect[-1])
t_fine_mesh = np.linspace(0, delta_t_vect[-1]*(nb_time_step-1), nb_time_step)
spline_bicubic = sp.interpolate.RectBivariateSpline(t_fine_mesh, r_fine_mesh, sol_MNP.C[:,:])
# ----------------------------------------------------------------------------------------
#       Schéma 1 : Discretisation d'ordre 1 en temps et 1 en espace étude ordre spatial
# ----------------------------------------------------------------------------------------

plt.figure(0)

delta_t_MMS = 0.001
t_final = 1.0
N_vect = np.array([20, 40, 60, 80, 100])
delta_r_vect = R/(N_vect-1)
Objet_Etude_Convergence = Etude_Convergence_MMS_spatial(delta_r_vect, delta_t_MMS, N_vect, R, t_final, 1)
erreur_vect_L1, erreur_vect_L2, erreur_vect_L_inf = Objet_Etude_Convergence.Boucle_iterations()

plt.figure(1)

# Graphique log-log norme de l'erreur L1 vs delta_r
plt.loglog(delta_r_vect, erreur_vect_L1, '.r', label = "Norme L1")

# Ajuster une loi de puissance à toutes les valeurs (en utilisant np.polyfit avec logarithmes)
coefficients = np.polyfit(np.log(delta_r_vect), np.log(erreur_vect_L1), 1)
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
equation_text = f'$L_1 = {np.exp(constant_logreg):.4f} \\times Δr^{{{exponent_logreg:.4f}}}$'
equation_text_obj = plt.text(0.5, 0.05, equation_text, fontsize=12, transform=plt.gca().transAxes, color='k')

# Graphique log-log norme de l'erreur L2 vs delta_r
plt.loglog(delta_r_vect, erreur_vect_L2, '.g', label = "Norme L2")

# Ajuster une loi de puissance à toutes les valeurs (en utilisant np.polyfit avec logarithmes)
coefficients = np.polyfit(np.log(delta_r_vect), np.log(erreur_vect_L2), 1)
exponent_logreg = coefficients[0]
constant_logreg = coefficients[1]

# Fonction de régression en termes de logarithmes
fit_function_log = lambda x: exponent_logreg * x + constant_logreg

# Fonction de régression en termes originaux
fit_function = lambda x: np.exp(fit_function_log(np.log(x)))

# Extrapoler la valeur prédite pour la dernière valeur de h_values
extrapolated_value = fit_function(delta_r_vect[-1])
plt.loglog(delta_r_vect, fit_function(delta_r_vect), linestyle='--', color='g')

# Afficher l'équation de la régression en loi de puissance pour la norme L1
equation_text = f'$L_2 = {np.exp(constant_logreg):.4f} \\times Δr^{{{exponent_logreg:.4f}}}$'
equation_text_obj = plt.text(0.5, 0.15, equation_text, fontsize=12, transform=plt.gca().transAxes, color='k')

# Graphique log-log norme de l'erreur Linf vs delta_r
plt.loglog(delta_r_vect, erreur_vect_L_inf, '.m', label='Norme $L_\infty$')

# Ajuster une loi de puissance à toutes les valeurs (en utilisant np.polyfit avec logarithmes)
coefficients = np.polyfit(np.log(delta_r_vect), np.log(erreur_vect_L_inf), 1)
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
equation_text = f'$L_\infty = {np.exp(constant_logreg):.4f} \\times Δr^{{{exponent_logreg:.4f}}}$'
equation_text_obj = plt.text(0.5, 0.25, equation_text, fontsize=12, transform=plt.gca().transAxes, color='k')

plt.xlabel("$\Delta$r")
plt.ylabel("Norme de l'erreur")
plt.legend()
plt.grid()
plt.title("Normes des erreurs L1, L2 et $L_\infty$ avec schéma d'ordre 1 en \n fonction de $\Delta$r avec méthode MMS")
plt.savefig("Norme_des_erreurs_Schema_1_MMS_Spatial.png")
plt.show()


# -----------------------------------------------------------------------------------------------------------------
#                           Solution avec le maillage le plus fin (pour MNP)
# -----------------------------------------------------------------------------------------------------------------
sol_MNP_Centree = Profil_Concentration_Centree(delta_r_vect[-1], delta_t_vect[-1], N_vect[-1], R, critere_convergence, critere_max_iter)
sol_MNP_Centree.Algorithme_Resolution()
nb_time_step = sol_MNP_Centree.C.shape[0]
r_fine_mesh = np.linspace(0, R, N_vect[-1])
t_fine_mesh = np.linspace(0, delta_t_vect[-1]*(nb_time_step-1), nb_time_step)
spline_bicubic_Centree = sp.interpolate.RectBivariateSpline(t_fine_mesh, r_fine_mesh, sol_MNP_Centree.C[:,:])

#%%
# ----------------------------------------------------------------------------------------
#          Schéma 2 : Discretisation d'ordre 1 en temps et 2 en espace étude ordre spatial
# ----------------------------------------------------------------------------------------

plt.figure(2)

Objet_Etude_Convergence = Etude_Convergence_MMS_spatial(delta_r_vect, delta_t_MMS, N_vect, R, t_final, 2)
erreur_vect_L1_Centree, erreur_vect_L2_Centree, erreur_vect_L_inf_Centree = Objet_Etude_Convergence.Boucle_iterations()    

plt.figure(3)

# Graphique log-log norme de l'erreur L1 vs delta_r
plt.loglog(delta_r_vect, erreur_vect_L1_Centree, '.r', label = "Norme L1")

# Ajuster une loi de puissance à toutes les valeurs (en utilisant np.polyfit avec logarithmes)
coefficients = np.polyfit(np.log(delta_r_vect), np.log(erreur_vect_L1_Centree), 1)
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
equation_text = f'$L_1 = {np.exp(constant_logreg):.4f} \\times Δr^{{{exponent_logreg:.4f}}}$'
equation_text_obj = plt.text(0.5, 0.05, equation_text, fontsize=12, transform=plt.gca().transAxes, color='k')

# Graphique log-log norme de l'erreur L2 vs delta_r
plt.loglog(delta_r_vect, erreur_vect_L2_Centree, '.g', label = "Nomre L2")

# Ajuster une loi de puissance à toutes les valeurs (en utilisant np.polyfit avec logarithmes)
coefficients = np.polyfit(np.log(delta_r_vect), np.log(erreur_vect_L2_Centree), 1)
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
equation_text = f'$L_2 = {np.exp(constant_logreg):.4f} \\times Δr^{{{exponent_logreg:.4f}}}$'
equation_text_obj = plt.text(0.5, 0.15, equation_text, fontsize=12, transform=plt.gca().transAxes, color='k')

# Graphique log-log norme de l'erreur Linf vs delta_r
plt.loglog(delta_r_vect, erreur_vect_L_inf_Centree, '.m', label = "Norme $L_\infty$")

# Ajuster une loi de puissance à toutes les valeurs (en utilisant np.polyfit avec logarithmes)
coefficients = np.polyfit(np.log(delta_r_vect), np.log(erreur_vect_L_inf_Centree), 1)
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
equation_text = f'$L_\infty = {np.exp(constant_logreg):.4f} \\times Δr^{{{exponent_logreg:.4f}}}$'
equation_text_obj = plt.text(0.5, 0.25, equation_text, fontsize=12, transform=plt.gca().transAxes, color='k')


plt.xlabel("$\Delta$r")
plt.ylabel("Norme de l'erreur")
plt.legend()
plt.grid()
plt.title("Normes des erreurs L1, L2 et $L_\infty$ avec schéma d'ordre 2 en \n fonction de $\Delta$r avec méthode MMS")
plt.savefig("Norme_des_erreurs_Schema_2_MMS_Spatial.png")
plt.show()

# ----------------------------------------------------------------------------------------
#       Schéma 1 : Discretisation d'ordre 1 en temps et 1 en espace étude ordre temporel
# ----------------------------------------------------------------------------------------

plt.figure(4)

delta_t_vect = np.array([0.05, 0.025, 0.01, 0.005, 0.001])
t_final = 1.0
Objet_Etude_Convergence = Etude_Convergence_MMS_temporel(delta_r, delta_t_vect, N, R, t_final, 1)
erreur_vect_L1, erreur_vect_L2, erreur_vect_L_inf = Objet_Etude_Convergence.Boucle_iterations()

plt.figure(5)

# Graphique log-log norme de l'erreur L1 vs delta_t
plt.loglog(delta_t_vect, erreur_vect_L1, '.r', label = "Norme L1")

# Ajuster une loi de puissance à toutes les valeurs (en utilisant np.polyfit avec logarithmes)
coefficients = np.polyfit(np.log(delta_t_vect), np.log(erreur_vect_L1), 1)
exponent_logreg = coefficients[0]
constant_logreg = coefficients[1]

# Fonction de régression en termes de logarithmes
fit_function_log = lambda x: exponent_logreg * x + constant_logreg

# Fonction de régression en termes originaux
fit_function = lambda x: np.exp(fit_function_log(np.log(x)))

# Extrapoler la valeur prédite pour la dernière valeur de h_values
extrapolated_value = fit_function(delta_t_vect[-1])
plt.loglog(delta_t_vect, fit_function(delta_t_vect), linestyle='--', color='r')

# Afficher l'équation de la régression en loi de puissance pour la norme L1
equation_text = f'$L_1 = {np.exp(constant_logreg):.4f} \\times Δr^{{{exponent_logreg:.4f}}}$'
equation_text_obj = plt.text(0.5, 0.05, equation_text, fontsize=12, transform=plt.gca().transAxes, color='k')

# Graphique log-log norme de l'erreur L2 vs delta_t
plt.loglog(delta_t_vect, erreur_vect_L2, '.g', label = "Norme L2")

# Ajuster une loi de puissance à toutes les valeurs (en utilisant np.polyfit avec logarithmes)
coefficients = np.polyfit(np.log(delta_t_vect), np.log(erreur_vect_L2), 1)
exponent_logreg = coefficients[0]
constant_logreg = coefficients[1]

# Fonction de régression en termes de logarithmes
fit_function_log = lambda x: exponent_logreg * x + constant_logreg

# Fonction de régression en termes originaux
fit_function = lambda x: np.exp(fit_function_log(np.log(x)))

# Extrapoler la valeur prédite pour la dernière valeur de h_values
extrapolated_value = fit_function(delta_t_vect[-1])
plt.loglog(delta_t_vect, fit_function(delta_t_vect), linestyle='--', color='g')

# Afficher l'équation de la régression en loi de puissance pour la norme L1
equation_text = f'$L_2 = {np.exp(constant_logreg):.4f} \\times Δr^{{{exponent_logreg:.4f}}}$'
equation_text_obj = plt.text(0.5, 0.15, equation_text, fontsize=12, transform=plt.gca().transAxes, color='k')

# Graphique log-log norme de l'erreur Linf vs delta_t
plt.loglog(delta_t_vect, erreur_vect_L_inf, '.m', label='Norme $L_\infty$')

# Ajuster une loi de puissance à toutes les valeurs (en utilisant np.polyfit avec logarithmes)
coefficients = np.polyfit(np.log(delta_t_vect), np.log(erreur_vect_L_inf), 1)
exponent_logreg = coefficients[0]
constant_logreg = coefficients[1]

# Fonction de régression en termes de logarithmes
fit_function_log = lambda x: exponent_logreg * x + constant_logreg

# Fonction de régression en termes originaux
fit_function = lambda x: np.exp(fit_function_log(np.log(x)))

# Extrapoler la valeur prédite pour la dernière valeur de h_values
extrapolated_value = fit_function(delta_t_vect[-1])
plt.loglog(delta_t_vect, fit_function(delta_t_vect), linestyle='--', color='m')

# Afficher l'équation de la régression en loi de puissance pour la norme Linf
equation_text = f'$L_\infty = {np.exp(constant_logreg):.4f} \\times Δr^{{{exponent_logreg:.4f}}}$'
equation_text_obj = plt.text(0.5, 0.25, equation_text, fontsize=12, transform=plt.gca().transAxes, color='k')

plt.xlabel("$\Delta$t")
plt.ylabel("Norme de l'erreur")
plt.legend()
plt.grid()
plt.title("Normes des erreurs L1, L2 et $L_\infty$ avec schéma d'ordre 1 en \n fonction de $\Delta$t avec méthode MMS")
plt.savefig("Norme_des_erreurs_Schema_1_MMS_Temp.png")
plt.show()

#%%
# ----------------------------------------------------------------------------------------
#          Schéma 2 : Discretisation d'ordre 1 en temps et 2 en espace étude ordre temporel
# ----------------------------------------------------------------------------------------

plt.figure(6)

Objet_Etude_Convergence = Etude_Convergence_MMS_temporel(delta_r, delta_t_vect, N, R, t_final, 2)
erreur_vect_L1_Centree, erreur_vect_L2_Centree, erreur_vect_L_inf_Centree = Objet_Etude_Convergence.Boucle_iterations()    

plt.figure(7)

# Graphique log-log norme de l'erreur L1 vs delta_t
plt.loglog(delta_t_vect, erreur_vect_L1_Centree, '.r', label = "Norme L1")

# Ajuster une loi de puissance à toutes les valeurs (en utilisant np.polyfit avec logarithmes)
coefficients = np.polyfit(np.log(delta_t_vect), np.log(erreur_vect_L1_Centree), 1)
exponent_logreg = coefficients[0]
constant_logreg = coefficients[1]

# Fonction de régression en termes de logarithmes
fit_function_log = lambda x: exponent_logreg * x + constant_logreg

# Fonction de régression en termes originaux
fit_function = lambda x: np.exp(fit_function_log(np.log(x)))

# Extrapoler la valeur prédite pour la dernière valeur de h_values
extrapolated_value = fit_function(delta_t_vect[-1])
plt.loglog(delta_t_vect, fit_function(delta_t_vect), linestyle='--', color='r')

# Afficher l'équation de la régression en loi de puissance pour la norme L1
equation_text = f'$L_1 = {np.exp(constant_logreg):.4f} \\times Δr^{{{exponent_logreg:.4f}}}$'
equation_text_obj = plt.text(0.5, 0.05, equation_text, fontsize=12, transform=plt.gca().transAxes, color='k')

# Graphique log-log norme de l'erreur L2 vs delta_t
plt.loglog(delta_t_vect, erreur_vect_L2_Centree, '.g', label = "Nomre L2")

# Ajuster une loi de puissance à toutes les valeurs (en utilisant np.polyfit avec logarithmes)
coefficients = np.polyfit(np.log(delta_t_vect), np.log(erreur_vect_L2_Centree), 1)
exponent_logreg = coefficients[0]
constant_logreg = coefficients[1]

# Fonction de régression en termes de logarithmes
fit_function_log = lambda x: exponent_logreg * x + constant_logreg

# Fonction de régression en termes originaux
fit_function = lambda x: np.exp(fit_function_log(np.log(x)))

# Extrapoler la valeur prédite pour la dernière valeur de h_values
extrapolated_value = fit_function(delta_t_vect[-1])
plt.loglog(delta_t_vect, fit_function(delta_t_vect), linestyle='--', color='g')

# Afficher l'équation de la régression en loi de puissance pour la norme L2
equation_text = f'$L_2 = {np.exp(constant_logreg):.4f} \\times Δr^{{{exponent_logreg:.4f}}}$'
equation_text_obj = plt.text(0.5, 0.15, equation_text, fontsize=12, transform=plt.gca().transAxes, color='k')

# Graphique log-log norme de l'erreur Linf vs delta_t
plt.loglog(delta_t_vect, erreur_vect_L_inf_Centree, '.m', label = "Norme $L_\infty$")

# Ajuster une loi de puissance à toutes les valeurs (en utilisant np.polyfit avec logarithmes)
coefficients = np.polyfit(np.log(delta_t_vect), np.log(erreur_vect_L_inf_Centree), 1)
exponent_logreg = coefficients[0]
constant_logreg = coefficients[1]

# Fonction de régression en termes de logarithmes
fit_function_log = lambda x: exponent_logreg * x + constant_logreg

# Fonction de régression en termes originaux
fit_function = lambda x: np.exp(fit_function_log(np.log(x)))

# Extrapoler la valeur prédite pour la dernière valeur de h_values
extrapolated_value = fit_function(delta_t_vect[-1])
plt.loglog(delta_t_vect, fit_function(delta_t_vect), linestyle='--', color='m')

# Afficher l'équation de la régression en loi de puissance pour la norme Linf
equation_text = f'$L_\infty = {np.exp(constant_logreg):.4f} \\times Δr^{{{exponent_logreg:.4f}}}$'
equation_text_obj = plt.text(0.5, 0.25, equation_text, fontsize=12, transform=plt.gca().transAxes, color='k')

plt.xlabel("$\Delta$t")
plt.ylabel("Norme de l'erreur")
plt.legend()
plt.grid()
plt.title("Normes des erreurs L1, L2 et $L_\infty$ avec schéma d'ordre 2 en \n fonction de $\Delta$t avec méthode MMS")
plt.savefig("Norme_des_erreurs_Schema_2_MMS_Temp.png")
plt.show()







