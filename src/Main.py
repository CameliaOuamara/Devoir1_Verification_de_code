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
#import scipy as sp
import matplotlib.pyplot as plt
import scipy.interpolate

# Functions
from profil_concentration import *
from plot import *
from norme_erreur_discretisation import *
from etude_convergence import *
import time
def compute_dummy_MMS_sol(t_fine_mesh,r_fine_mesh):
    C_out = np.zeros((1, len(r_fine_mesh)))
    for ri in range(0,len(r_fine_mesh)):
        r = r_fine_mesh[ri]
        tempval = np.exp(-1.0*0.0*(1.0/1.0e9))*(0.5*np.cos(r*2*np.pi/r_fine_mesh[-1]))+(np.cos(r*2*np.pi/r_fine_mesh[-1]))
        C_out[0][ri]=tempval
    for ti in range(1,len(t_fine_mesh)):
        t = t_fine_mesh[ti]
        t_arr_temp = np.zeros((1, len(r_fine_mesh)))
        for ri in range(0,len(r_fine_mesh)):
            r = r_fine_mesh[ri]
            tempval = np.exp(-1.0*t*(1.0/1.0e9))*(0.5*np.cos(r*2*np.pi/r_fine_mesh[-1]))+(np.cos(r*2*np.pi/r_fine_mesh[-1]))
            t_arr_temp[0][ri]=tempval
        C_out=np.append(C_out,t_arr_temp,axis=0)
    return C_out
        

# -----------------------------------------------------------------------------------------------------------------
#                                               Debut du code
# -----------------------------------------------------------------------------------------------------------------
time_start = time.time()
# Données 
R = 0.5                                              # Rayon du cylindre [m]
N = 100                                              # Nombre de points de discretisation [-]
delta_r = R/(N-1)                                    # Taille de l'intervalle geometrique [m]
D_eff = 1.0e-10                                      # Coefficient de diffusion effectif [m2/s]
outputFolder = "/home/apollon/vilig/cours_VnV/devoir2/tracked_dir/Devoir1_Verification_de_code/results/"

# Règle de pouce : Prendre un dt à la limite de la stabilité pour le schéma explicite
delta_t = 0.5 * delta_r*delta_r / D_eff              # Pas de temps [s]

# Critere de convergence
critere_convergence = 1.0e-14                        # Critere sur la valeur de la concentration a l'iteration i vs i-1
critere_max_iter = 400                              # Nombre minimum d'iterations a realiser pour la convergence vers le regime permanent

##Input conditions for convergence study HW2:
convergence_dimension = "time" # switch to "time" for time, "space" for space
manufactured_sol_type = "MNP" # switch to "MNP" for near problem, "MMS" for manufactured solution
schema_check = "centered" # switch to "centered" for formulation of HW1QF, "notcentered" for HW1QA


# Étude convergence MNP
##### for space study
if convergence_dimension == "space" and manufactured_sol_type == "MNP":
    N_vect = np.arange(7,11,1, dtype=int)                 # Vecteur contenant les N utilisés dans l'étude de convergence
    N_vect = 5 * 2**N_vect
    delta_r_vect = R/(N_vect-1)                          # Vecteur Delta r correspondant au vecteur N precedent
    delta_t_vect = 1.0e1 * 2.0 +0.0*delta_r_vect*delta_r_vect / D_eff # Vecteur Delta t correspondant au vecteur Delta r precedent
    critere_max_iter = 1000+0.0*delta_r_vect
    nombre_nodes_vect = N_vect

##### for time study
if convergence_dimension == "time" and manufactured_sol_type == "MNP":
    N_vect = np.arange(0,5,1, dtype=int)                 # Vecteur contenant les N utilisés dans l'étude de convergence
    phys_time_simul_end = 2048002000
    N_vect = 160 * 2**N_vect
    critere_max_iter = N_vect
    delta_t_vect = phys_time_simul_end/(N_vect+1)                          # Vecteur Delta r correspondant au vecteur N precedent
    delta_r_vect = R/(5.0*(2.0**8.0)-1.0) +0.0*N_vect # Vecteur Delta t correspondant au vecteur Delta r precedent
    nombre_nodes_vect = (5 * 2**8)+0*N_vect


# Étude convergence MMS
##### for space study
if convergence_dimension == "space" and manufactured_sol_type == "MMS":
    N_vect = np.arange(1,7,1, dtype=int)                 # Vecteur contenant les N utilisés dans l'étude de convergence
    N_vect = 11 * 2**N_vect
    delta_r_vect = R/(N_vect-1)                          # Vecteur Delta r correspondant au vecteur N precedent
    delta_t_vect = 1.0e4 * 3.125 +0.0*delta_r_vect*delta_r_vect / D_eff # Vecteur Delta t correspondant au vecteur Delta r precedent
    critere_max_iter = 400+0.0*delta_r_vect
    nombre_nodes_vect = N_vect
##### for time study
if convergence_dimension == "time" and manufactured_sol_type == "MMS":
    N_vect = np.arange(0,5,1, dtype=int)                 # Vecteur contenant les N utilisés dans l'étude de convergence
    phys_time_simul_end = 1.0e7 * 3.125 * 401
    N_vect = 50 * 2**N_vect
    critere_max_iter = N_vect
    delta_t_vect = phys_time_simul_end/(N_vect+1)                          # Vecteur Delta r correspondant au vecteur N precedent
    delta_r_vect = R/(11.0*(2.0**8.0)-1.0) +0.0*N_vect # Vecteur Delta t correspondant au vecteur Delta r precedent
    nombre_nodes_vect = (11 * 2**8)+0*N_vect

print("convergence study values:")
print(delta_r_vect)
print(delta_t_vect)
print(critere_max_iter)
print(nombre_nodes_vect)
'''
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
#          Schéma 1 : Discretisation d'ordre 1 en temps et 1 en espace
# ----------------------------------------------------------------------------------------

plt.figure(0)

Objet_Etude_Convergence = Etude_Convergence(delta_r_vect, delta_t_vect, N_vect, R, critere_convergence, critere_max_iter, 1, sol_MNP.C[-1,:], spline_bicubic)
erreur_vect_L1, erreur_vect_L2, erreur_vect_L_inf = Objet_Etude_Convergence.Boucle_iterations(outputFolder)

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
equation_text = f'$L_1 = {np.exp(constant_logreg):.4E} \\times Δr^{{{exponent_logreg:.4f}}}$'
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

# Afficher l'équation de la régression en loi de puissance pour la norme L2
equation_text = f'$L_2 = {np.exp(constant_logreg):.4E} \\times Δr^{{{exponent_logreg:.4f}}}$'
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
equation_text = f'$L_\infty = {np.exp(constant_logreg):.4E} \\times Δr^{{{exponent_logreg:.4f}}}$'
equation_text_obj = plt.text(0.5, 0.25, equation_text, fontsize=12, transform=plt.gca().transAxes, color='k')

plt.xlabel("delta_r")
plt.ylabel("Norme de l'erreur")
plt.legend()
plt.grid()
plt.title("Normes des erreurs L1, L2 et $L_\infty$ schéma d'ordre 1 en fonction de N")
plt.savefig(outputFolder+"Norme_des_erreurs_Schema_1.png")
plt.show()
'''
# -----------------------------------------------------------------------------------------------------------------
#                           Solution avec le maillage le plus fin (pour MNP)
# -----------------------------------------------------------------------------------------------------------------
## input MNP delta_r_convergence
if convergence_dimension == "space" and manufactured_sol_type == "MNP":
    decalage_ignore = 20000
    decalage_fine_mesh = -4
    # decalage_ignore = 500
    # decalage_fine_mesh = -4
## input MNP delta_t_convergence
if convergence_dimension == "time" and manufactured_sol_type == "MNP":
    decalage_ignore = 0
    decalage_fine_mesh = -5
# input MMS convergence
if manufactured_sol_type == "MMS":
    decalage_ignore = 0
    decalage_fine_mesh = -1

sol_MNP_Centree = Profil_Concentration_Centree(delta_r_vect[decalage_fine_mesh], delta_t_vect[decalage_fine_mesh], nombre_nodes_vect[decalage_fine_mesh], R, critere_convergence, critere_max_iter[decalage_fine_mesh]+decalage_ignore)
sol_MNP_Centree.Algorithme_Resolution()
nb_time_step = sol_MNP_Centree.C.shape[0]-decalage_ignore
r_fine_mesh = np.linspace(0, R, nombre_nodes_vect[decalage_fine_mesh])
t_fine_mesh = np.linspace(0, delta_t_vect[decalage_fine_mesh]*(nb_time_step-1), nb_time_step)
time_before_compute_dummy = time.time()
C_dummy = compute_dummy_MMS_sol(t_fine_mesh,r_fine_mesh)
time_after_compute_dummy = time.time()
print("time to compute dummy:")
print(time_after_compute_dummy - time_before_compute_dummy)
if manufactured_sol_type == "MNP":
    spline_bicubic_Centree = sp.interpolate.RectBivariateSpline(t_fine_mesh, r_fine_mesh, sol_MNP_Centree.C[decalage_ignore:,:])
if manufactured_sol_type == "MMS":
    spline_bicubic_Centree = sp.interpolate.RectBivariateSpline(t_fine_mesh, r_fine_mesh, C_dummy[:,:])

print("Time after spline MNP:")
time_post_MNP_spline = time.time()
print(time_post_MNP_spline-time_start)
#%%
# ----------------------------------------------------------------------------------------
#          Schéma 2 : Discretisation d'ordre 1 en temps et 2 en espace
# ----------------------------------------------------------------------------------------
# erreur_vect_L1_Centree = np.zeros(len(N_vect))
# erreur_vect_L2_Centree = np.zeros(len(N_vect))
# erreur_vect_L_inf_Centree = np.zeros(len(N_vect))

plt.figure(2)
if schema_check == "centered":
    Objet_Etude_Convergence = Etude_Convergence(delta_r_vect, delta_t_vect, N_vect, R, critere_convergence, critere_max_iter, 2, sol_MNP_Centree.C, spline_bicubic_Centree, nombre_nodes_vect)
if schema_check == "notcentered":
    Objet_Etude_Convergence = Etude_Convergence(delta_r_vect, delta_t_vect, N_vect, R, critere_convergence, critere_max_iter, 1, sol_MNP_Centree.C, spline_bicubic_Centree, nombre_nodes_vect)
time_post_crea_etude_convergence = time.time()
print("Time after convergence etude cration:")
print(time_post_crea_etude_convergence-time_start)
erreur_vect_L1_Centree, erreur_vect_L2_Centree, erreur_vect_L_inf_Centree = Objet_Etude_Convergence.Boucle_iterations(outputFolder)
time_post_boucle_etude_convergence = time.time()
print("Time after convergence etude boucle:")
print(time_post_boucle_etude_convergence-time_start)
    

plt.figure(3)
if convergence_dimension == "space":
    array_x_conv = delta_r_vect
if convergence_dimension == "time":
    array_x_conv = delta_t_vect
# Graphique log-log norme de l'erreur L1 vs delta_r
plt.loglog(array_x_conv, erreur_vect_L1_Centree, '.r', label = "Norme L1")

# Ajuster une loi de puissance à toutes les valeurs (en utilisant np.polyfit avec logarithmes)
coefficients = np.polyfit(np.log(array_x_conv[0:-1]), np.log(erreur_vect_L1_Centree[0:-1]), 1)
exponent_logreg = coefficients[0]
constant_logreg = coefficients[1]

# Fonction de régression en termes de logarithmes
fit_function_log = lambda x: exponent_logreg * x + constant_logreg

# Fonction de régression en termes originaux
fit_function = lambda x: np.exp(fit_function_log(np.log(x)))

# Extrapoler la valeur prédite pour la dernière valeur de h_values
extrapolated_value = fit_function(array_x_conv[-1])
plt.loglog(array_x_conv, fit_function(array_x_conv), linestyle='--', color='r')

# Afficher l'équation de la régression en loi de puissance pour la norme L1
if convergence_dimension == "space":
    equation_text = f'$L_1 = {np.exp(constant_logreg):.4E} \\times Δr^{{{exponent_logreg:.4f}}}$'
    equation_text_obj = plt.text(0.05, 0.05, equation_text, fontsize=12, transform=plt.gca().transAxes, color='k')
if convergence_dimension == "time":
    equation_text = f'$L_1 = {np.exp(constant_logreg):.4E} \\times Δt^{{{exponent_logreg:.4f}}}$'
    equation_text_obj = plt.text(0.05, 0.05, equation_text, fontsize=12, transform=plt.gca().transAxes, color='k')

# Graphique log-log norme de l'erreur L2 vs delta_r
plt.loglog(array_x_conv, erreur_vect_L2_Centree, '.g', label = "Nomre L2")

# Ajuster une loi de puissance à toutes les valeurs (en utilisant np.polyfit avec logarithmes)
coefficients = np.polyfit(np.log(array_x_conv[0:-1]), np.log(erreur_vect_L2_Centree[0:-1]), 1)
exponent_logreg = coefficients[0]
constant_logreg = coefficients[1]

# Fonction de régression en termes de logarithmes
fit_function_log = lambda x: exponent_logreg * x + constant_logreg

# Fonction de régression en termes originaux
fit_function = lambda x: np.exp(fit_function_log(np.log(x)))

# Extrapoler la valeur prédite pour la dernière valeur de h_values
extrapolated_value = fit_function(array_x_conv[-1])
plt.loglog(array_x_conv, fit_function(array_x_conv), linestyle='--', color='g')

# Afficher l'équation de la régression en loi de puissance pour la norme L2
if convergence_dimension == "space":
    equation_text = f'$L_2 = {np.exp(constant_logreg):.4E} \\times Δr^{{{exponent_logreg:.4f}}}$'
    equation_text_obj = plt.text(0.05, 0.15, equation_text, fontsize=12, transform=plt.gca().transAxes, color='k')
if convergence_dimension == "time":
    equation_text = f'$L_2 = {np.exp(constant_logreg):.4E} \\times Δt^{{{exponent_logreg:.4f}}}$'
    equation_text_obj = plt.text(0.05, 0.15, equation_text, fontsize=12, transform=plt.gca().transAxes, color='k')


# Graphique log-log norme de l'erreur Linf vs delta_r
plt.loglog(array_x_conv, erreur_vect_L_inf_Centree, '.m', label = "Norme $L_\infty$")

# Ajuster une loi de puissance à toutes les valeurs (en utilisant np.polyfit avec logarithmes)
coefficients = np.polyfit(np.log(array_x_conv[0:-1]), np.log(erreur_vect_L_inf_Centree[0:-1]), 1)
exponent_logreg = coefficients[0]
constant_logreg = coefficients[1]

# Fonction de régression en termes de logarithmes
fit_function_log = lambda x: exponent_logreg * x + constant_logreg

# Fonction de régression en termes originaux
fit_function = lambda x: np.exp(fit_function_log(np.log(x)))

# Extrapoler la valeur prédite pour la dernière valeur de h_values
extrapolated_value = fit_function(array_x_conv[-1])
plt.loglog(array_x_conv, fit_function(array_x_conv), linestyle='--', color='m')

# Afficher l'équation de la régression en loi de puissance pour la norme Linf
if convergence_dimension == "space":
    equation_text = f'$Linf = {np.exp(constant_logreg):.4E} \\times Δr^{{{exponent_logreg:.4f}}}$'
    equation_text_obj = plt.text(0.05, 0.25, equation_text, fontsize=12, transform=plt.gca().transAxes, color='k')
if convergence_dimension == "time":
    equation_text = f'$Linf = {np.exp(constant_logreg):.4E} \\times Δt^{{{exponent_logreg:.4f}}}$'
    equation_text_obj = plt.text(0.05, 0.25, equation_text, fontsize=12, transform=plt.gca().transAxes, color='k')

if convergence_dimension == "time":
    plt.xlabel("delta_t")
if convergence_dimension == "space":
    plt.xlabel("delta_r")
plt.ylabel("erreur")
plt.legend()
plt.grid()
plt.title("Normes des erreurs L1, L2 et $L_\infty$ schéma d'ordre 2 en fonction de N")
if convergence_dimension == "time" and manufactured_sol_type == "MNP":
    plt.title("Normes des erreurs L1, L2 et $L_\infty$ schéma d'ordre 2 en fonction de Δt pour MNP")
    plt.savefig(outputFolder+"Norme_des_erreurs_Schema_2_MNP_time.png")
if convergence_dimension == "space" and manufactured_sol_type == "MNP":
    plt.title("Normes des erreurs L1, L2 et $L_\infty$ schéma d'ordre 2 en fonction de Δr pour MNP")
    plt.savefig(outputFolder+"Norme_des_erreurs_Schema_2_MNP_space.png")
if convergence_dimension == "time" and manufactured_sol_type == "MMS":
    plt.title("Normes des erreurs L1, L2 et $L_\infty$ schéma d'ordre 2 en fonction de Δt pour MMS")
    plt.savefig(outputFolder+"Norme_des_erreurs_Schema_2_MMS_time.png")
if convergence_dimension == "space" and manufactured_sol_type == "MMS":
    plt.title("Normes des erreurs L1, L2 et $L_\infty$ schéma d'ordre 2 en fonction de Δr pour MMS")
    plt.savefig(outputFolder+"Norme_des_erreurs_Schema_2_MMS_space.png")
plt.show()

time_end = time.time()
print("Time end:")
print(time_end-time_start)






