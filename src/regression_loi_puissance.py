# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 14:45:50 2024

@author: ouama
"""

# Libraries
import numpy as np
#import sympy as sp
import matplotlib.pyplot as plt
import pandas as pd

class Regression_Loi_Puissance():
    def __init__(self, L1, L2, Linf, schema, delta_r_vect, delta_t_vect):
        self.L1 = L1
        self.L2 = L2
        self.Linf = Linf
        self.schema = schema
        self.delta_r_vect = delta_r_vect
        self.delta_t_vect = delta_t_vect
        
    def Plot_Regression_Espace(self):
        # Graphique log-log norme de l'erreur L1 vs delta_r
        plt.loglog(self.delta_r_vect, self.L1, '.r', label = "Norme L1")

        # Ajuster une loi de puissance à toutes les valeurs (en utilisant np.polyfit avec logarithmes)
        coefficients = np.polyfit(np.log(self.delta_r_vect), np.log(self.L1), 1)
        exponent_logreg = coefficients[0]
        constant_logreg = coefficients[1]

        # Fonction de régression en termes de logarithmes
        fit_function_log = lambda x: exponent_logreg * x + constant_logreg

        # Fonction de régression en termes originaux
        fit_function = lambda x: np.exp(fit_function_log(np.log(x)))

        # Extrapoler la valeur prédite pour la dernière valeur de h_values
        extrapolated_value = fit_function(self.delta_r_vect[-1])
        plt.loglog(self.delta_r_vect, fit_function(self.delta_r_vect), linestyle='--', color='r')

        # Afficher l'équation de la régression en loi de puissance pour la norme L1
        equation_text = f'$L_1 = {np.exp(constant_logreg):.4E} \\times Δr^{{{exponent_logreg:.4f}}}$'
        equation_text_obj = plt.text(0.5, 0.05, equation_text, fontsize=12, transform=plt.gca().transAxes, color='k')

        # Graphique log-log norme de l'erreur L2 vs delta_r
        plt.loglog(self.delta_r_vect, self.L2, '.g', label = "Norme L2")

        # Ajuster une loi de puissance à toutes les valeurs (en utilisant np.polyfit avec logarithmes)
        coefficients = np.polyfit(np.log(self.delta_r_vect), np.log(self.L2), 1)
        exponent_logreg = coefficients[0]
        constant_logreg = coefficients[1]

        # Fonction de régression en termes de logarithmes
        fit_function_log = lambda x: exponent_logreg * x + constant_logreg

        # Fonction de régression en termes originaux
        fit_function = lambda x: np.exp(fit_function_log(np.log(x)))

        # Extrapoler la valeur prédite pour la dernière valeur de h_values
        extrapolated_value = fit_function(self.delta_r_vect[-1])
        plt.loglog(self.delta_r_vect, fit_function(self.delta_r_vect), linestyle='--', color='g')

        # Afficher l'équation de la régression en loi de puissance pour la norme L2
        equation_text = f'$L_2 = {np.exp(constant_logreg):.4E} \\times Δr^{{{exponent_logreg:.4f}}}$'
        equation_text_obj = plt.text(0.5, 0.15, equation_text, fontsize=12, transform=plt.gca().transAxes, color='k')

        # Graphique log-log norme de l'erreur Linf vs delta_r
        plt.loglog(self.delta_r_vect, self.Linf, '.m', label='Norme $L_\infty$')

        # Ajuster une loi de puissance à toutes les valeurs (en utilisant np.polyfit avec logarithmes)
        coefficients = np.polyfit(np.log(self.delta_r_vect), np.log(self.Linf), 1)
        exponent_logreg = coefficients[0]
        constant_logreg = coefficients[1]

        # Fonction de régression en termes de logarithmes
        fit_function_log = lambda x: exponent_logreg * x + constant_logreg

        # Fonction de régression en termes originaux
        fit_function = lambda x: np.exp(fit_function_log(np.log(x)))

        # Extrapoler la valeur prédite pour la dernière valeur de h_values
        extrapolated_value = fit_function(self.delta_r_vect[-1])
        plt.loglog(self.delta_r_vect, fit_function(self.delta_r_vect), linestyle='--', color='m')

        # Afficher l'équation de la régression en loi de puissance pour la norme Linf
        equation_text = f'$L_\infty = {np.exp(constant_logreg):.4E} \\times Δr^{{{exponent_logreg:.4f}}}$'
        equation_text_obj = plt.text(0.5, 0.25, equation_text, fontsize=12, transform=plt.gca().transAxes, color='k')
        
        plt.xlabel("$\Delta r$")
        plt.ylabel("Norme de l'erreur")
        plt.legend(loc=2)
        plt.grid()
        plt.title("Normes des erreurs L1, L2 et $L_\infty$ en fonction de $\Delta r$ schéma %d " %self.schema)
        plt.show()

    def Plot_Regression_Temps(self):
        # Graphique log-log norme de l'erreur L1 vs delta_r
        plt.loglog(self.delta_t_vect, self.L1, '.r', label = "Norme L1")

        # Ajuster une loi de puissance à toutes les valeurs (en utilisant np.polyfit avec logarithmes)
        coefficients = np.polyfit(np.log(self.delta_t_vect), np.log(self.L1), 1)
        exponent_logreg = coefficients[0]
        constant_logreg = coefficients[1]

        # Fonction de régression en termes de logarithmes
        fit_function_log = lambda x: exponent_logreg * x + constant_logreg

        # Fonction de régression en termes originaux
        fit_function = lambda x: np.exp(fit_function_log(np.log(x)))

        # Extrapoler la valeur prédite pour la dernière valeur de h_values
        extrapolated_value = fit_function(self.delta_t_vect[-1])
        plt.loglog(self.delta_t_vect, fit_function(self.delta_t_vect), linestyle='--', color='r')

        # Afficher l'équation de la régression en loi de puissance pour la norme L1
        equation_text = f'$L_1 = {np.exp(constant_logreg):.4E} \\times Δt^{{{exponent_logreg:.4f}}}$'
        equation_text_obj = plt.text(0.5, 0.05, equation_text, fontsize=12, transform=plt.gca().transAxes, color='k')

        # Graphique log-log norme de l'erreur L2 vs delta_r
        plt.loglog(self.delta_t_vect, self.L2, '.g', label = "Norme L2")

        # Ajuster une loi de puissance à toutes les valeurs (en utilisant np.polyfit avec logarithmes)
        coefficients = np.polyfit(np.log(self.delta_t_vect), np.log(self.L2), 1)
        exponent_logreg = coefficients[0]
        constant_logreg = coefficients[1]

        # Fonction de régression en termes de logarithmes
        fit_function_log = lambda x: exponent_logreg * x + constant_logreg

        # Fonction de régression en termes originaux
        fit_function = lambda x: np.exp(fit_function_log(np.log(x)))

        # Extrapoler la valeur prédite pour la dernière valeur de h_values
        extrapolated_value = fit_function(self.delta_t_vect[-1])
        plt.loglog(self.delta_t_vect, fit_function(self.delta_t_vect), linestyle='--', color='g')

        # Afficher l'équation de la régression en loi de puissance pour la norme L2
        equation_text = f'$L_2 = {np.exp(constant_logreg):.4E} \\times Δt^{{{exponent_logreg:.4f}}}$'
        equation_text_obj = plt.text(0.5, 0.15, equation_text, fontsize=12, transform=plt.gca().transAxes, color='k')

        # Graphique log-log norme de l'erreur Linf vs delta_r
        plt.loglog(self.delta_t_vect, self.Linf, '.m', label='Norme $L_\infty$')

        # Ajuster une loi de puissance à toutes les valeurs (en utilisant np.polyfit avec logarithmes)
        coefficients = np.polyfit(np.log(self.delta_t_vect), np.log(self.Linf), 1)
        exponent_logreg = coefficients[0]
        constant_logreg = coefficients[1]

        # Fonction de régression en termes de logarithmes
        fit_function_log = lambda x: exponent_logreg * x + constant_logreg

        # Fonction de régression en termes originaux
        fit_function = lambda x: np.exp(fit_function_log(np.log(x)))

        # Extrapoler la valeur prédite pour la dernière valeur de h_values
        extrapolated_value = fit_function(self.delta_t_vect[-1])
        plt.loglog(self.delta_t_vect, fit_function(self.delta_t_vect), linestyle='--', color='m')

        # Afficher l'équation de la régression en loi de puissance pour la norme Linf
        equation_text = f'$L_\infty = {np.exp(constant_logreg):.4E} \\times Δt^{{{exponent_logreg:.4f}}}$'
        equation_text_obj = plt.text(0.5, 0.25, equation_text, fontsize=12, transform=plt.gca().transAxes, color='k')
        
        plt.xlabel("$\Delta t$")
        plt.ylabel("Norme de l'erreur")
        plt.legend(loc=2)
        plt.grid()
        plt.title("Normes des erreurs L1, L2 et $L_\infty$ en fonction de $\Delta t$ schéma %d " %self.schema)
        plt.show()