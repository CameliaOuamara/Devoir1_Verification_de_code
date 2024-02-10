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

critere_conv = 1.0e-14
critere_max_iter = 100

# # Resolution
# Objet_Concentration = Profil_Concentration(delta_r, delta_t, N, R, critere_conv)
# Objet_Concentration.Algorithme_Resolution()

# # Plot
# Objet_Graphique = Plot_Concentration( Objet_Concentration.C, N)
# Objet_Graphique.Plot_Numerique()
# Objet_Graphique.Plot_Exact()

# Étude convergence
N_vect = np.arange(0,7,1, dtype=int)

N_vect = 5 * 2**N_vect


delta_r_vect = R/(N_vect-1)

delta_t_vect = 1.0e5 * 0.5 * delta_r_vect*delta_r_vect / D_eff

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
    Objet_Graphique.Save_plot("schema1_"+str(N_vect[i]), "Comparaison de résultat premier schéma, "+str(N_vect[i])+" noeuds")
    
    #Erreur
    erreur_vect_L1[i] = np.mean(abs(Objet_Concentration.C[-1,:] - Objet_Graphique.C_exact))
    erreur_vect_L2[i] = np.sqrt(np.mean((Objet_Concentration.C[-1,:] - Objet_Graphique.C_exact)*(Objet_Concentration.C[-1,:] - Objet_Graphique.C_exact)))
    erreur_vect_L_inf[i] = np.max(abs(Objet_Concentration.C[-1,:] - Objet_Graphique.C_exact))

    del Objet_Concentration
    del Objet_Graphique
    


delta_r_plot = delta_r_vect[-1] - delta_r_vect[0]
ordre = 1.0

first_pt_plot_theo_L1 = erreur_vect_L1[-1] * (delta_r_vect[0]/delta_r_vect[-1])**ordre
first_pt_plot_theo_L2 = erreur_vect_L2[-1] * (delta_r_vect[0]/delta_r_vect[-1])**ordre
first_pt_plot_theo_L_inf = erreur_vect_L_inf[-1] * (delta_r_vect[0]/delta_r_vect[-1])**ordre

plt.figure(1)
plt.loglog(delta_r_vect, erreur_vect_L1, '.r', label = "résultat numérique")
#plt.loglog([delta_r_vect[0], delta_r_vect[-1]], [first_pt_plot_theo_L1, erreur_vect_L1[-1]], label = "résultat attendu selon ordre théorique")

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
plt.loglog(delta_r_vect, fit_function(delta_r_vect), linestyle='--', color='b', label='Régression en loi de puissance')
# Afficher l'équation de la régression en loi de puissance
equation_text = f'$L_1 = {np.exp(constant_logreg):.4f} \\times Δr^{{{exponent_logreg:.4f}}}$'
equation_text_obj = plt.text(0.05, 0.05, equation_text, fontsize=12, transform=plt.gca().transAxes, color='k')

plt.xlabel("delta_r")
plt.ylabel("erreur_L1")
plt.legend()
plt.title("Visualisation erreur L1 premier schéma contre solution analytique")
plt.savefig("erreur_L1_schema1.png")

plt.figure(2)
plt.loglog(delta_r_vect, erreur_vect_L2, '.r', label = "résultat numérique")
#plt.loglog([delta_r_vect[0], delta_r_vect[-1]], [first_pt_plot_theo_L2, erreur_vect_L2[-1]], label = "résultat attendu selon ordre théorique")

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
plt.loglog(delta_r_vect, fit_function(delta_r_vect), linestyle='--', color='b', label='Régression en loi de puissance')
# Afficher l'équation de la régression en loi de puissance
equation_text = f'$L_2 = {np.exp(constant_logreg):.4f} \\times Δr^{{{exponent_logreg:.4f}}}$'
equation_text_obj = plt.text(0.05, 0.05, equation_text, fontsize=12, transform=plt.gca().transAxes, color='k')

plt.xlabel("delta_r")
plt.ylabel("erreur_L2")
plt.legend()
plt.title("Visualisation erreur L2 premier schéma contre solution analytique")
plt.savefig("erreur_L2_schema1.png")

plt.figure(3)
plt.loglog(delta_r_vect, erreur_vect_L_inf, '.r', label = "résultat numérique")
#plt.loglog([delta_r_vect[0], delta_r_vect[-1]], [first_pt_plot_theo_L_inf, erreur_vect_L_inf[-1]], label = "résultat attendu selon ordre théorique")

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
plt.loglog(delta_r_vect, fit_function(delta_r_vect), linestyle='--', color='b', label='Régression en loi de puissance')
# Afficher l'équation de la régression en loi de puissance
equation_text = f'$Linf = {np.exp(constant_logreg):.4f} \\times Δr^{{{exponent_logreg:.4f}}}$'
equation_text_obj = plt.text(0.05, 0.05, equation_text, fontsize=12, transform=plt.gca().transAxes, color='k')

plt.xlabel("delta_r")
plt.ylabel("erreur L_inf")
plt.legend()
plt.title("Visualisation erreur Linf premier schéma contre solution analytique")
plt.savefig("erreur_Linf_schema1.png")



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
    Objet_Graphique.Save_plot("schema2_"+str(N_vect[i]), "Comparaison de résultat deuxième schéma, "+str(N_vect[i])+" noeuds")
    
    #Erreur
    erreur_vect_L1_Centree[i] = np.mean(abs(Objet_Concentration.C[-1,:] - Objet_Graphique.C_exact))
    erreur_vect_L2_Centree[i] = np.sqrt(np.mean((Objet_Concentration.C[-1,:] - Objet_Graphique.C_exact)*(Objet_Concentration.C[-1,:] - Objet_Graphique.C_exact)))
    erreur_vect_L_inf_Centree[i] = np.max(abs(Objet_Concentration.C[-1,:] - Objet_Graphique.C_exact))

    del Objet_Concentration
    del Objet_Graphique
    

ordre_centree = 2.0

first_pt_plot_theo_L1_centree = erreur_vect_L1_Centree[-1] * (delta_r_vect[0]/delta_r_vect[-1])**ordre_centree
first_pt_plot_theo_L2_centree = erreur_vect_L2_Centree[-1] * (delta_r_vect[0]/delta_r_vect[-1])**ordre_centree
first_pt_plot_theo_L_inf_centree = erreur_vect_L_inf_Centree[-1] * (delta_r_vect[0]/delta_r_vect[-1])**ordre_centree

plt.figure(5)
plt.loglog(delta_r_vect, erreur_vect_L1_Centree, '.r', label = "résultat numérique")
#plt.loglog([delta_r_vect[0], delta_r_vect[-1]], [first_pt_plot_theo_L1_centree, erreur_vect_L1_Centree[-1]], label = "résultat attendu selon ordre théorique")

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
plt.loglog(delta_r_vect, fit_function(delta_r_vect), linestyle='--', color='b', label='Régression en loi de puissance')
# Afficher l'équation de la régression en loi de puissance
equation_text = f'$L_1 = {np.exp(constant_logreg):.4f} \\times Δr^{{{exponent_logreg:.4f}}}$'
equation_text_obj = plt.text(0.05, 0.05, equation_text, fontsize=12, transform=plt.gca().transAxes, color='k')#

plt.xlabel("delta_r")
plt.ylabel("erreur_L1")
plt.legend()
plt.title("Visualisation erreur L1 deuxième schéma contre solution analytique")
plt.savefig("erreur_L1_schema2.png")

plt.figure(6)
plt.loglog(delta_r_vect, erreur_vect_L2_Centree, '.r', label = "résultat numérique")
#plt.loglog([delta_r_vect[0], delta_r_vect[-1]], [first_pt_plot_theo_L2_centree, erreur_vect_L2_Centree[-1]], label = "résultat attendu selon ordre théorique")

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
plt.loglog(delta_r_vect, fit_function(delta_r_vect), linestyle='--', color='b', label='Régression en loi de puissance')
# Afficher l'équation de la régression en loi de puissance
equation_text = f'$L_2 = {np.exp(constant_logreg):.4f} \\times Δr^{{{exponent_logreg:.4f}}}$'
equation_text_obj = plt.text(0.05, 0.05, equation_text, fontsize=12, transform=plt.gca().transAxes, color='k')

plt.xlabel("delta_r")
plt.ylabel("erreur_L2")
plt.legend()
plt.title("Visualisation erreur L2 deuxième schéma contre solution analytique")
plt.savefig("erreur_L2_schema2.png")

plt.figure(7)
plt.loglog(delta_r_vect, erreur_vect_L_inf_Centree, '.r', label = "résultat numérique")
#plt.loglog([delta_r_vect[0], delta_r_vect[-1]], [first_pt_plot_theo_L_inf_centree, erreur_vect_L_inf_Centree[-1]], label = "résultat attendu selon ordre théorique")

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
plt.loglog(delta_r_vect, fit_function(delta_r_vect), linestyle='--', color='b', label='Régression en loi de puissance')
# Afficher l'équation de la régression en loi de puissance
equation_text = f'$Linf = {np.exp(constant_logreg):.4f} \\times Δr^{{{exponent_logreg:.4f}}}$'
equation_text_obj = plt.text(0.05, 0.05, equation_text, fontsize=12, transform=plt.gca().transAxes, color='k')

plt.xlabel("delta_r")
plt.ylabel("erreur L_inf")
plt.legend()
plt.title("Visualisation erreur L_inf deuxième schéma contre solution analytique")
plt.savefig("erreur_Linf_schema2.png")
plt.show()












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



