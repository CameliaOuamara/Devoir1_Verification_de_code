"""
------------------------------------------------------------------------------------------------------------------
MEC-8211 Vérification et Validation en Modélisation Numrérique
Devoir 3, QA
------------------------------------------------------------------------------------------------------------------

Ce code utilise un fichier de résultat précédemment obtenu pour réaliser une étude de convergence
sur la taille de maille.

"""

# Libraries
import time
import matplotlib.pyplot as plt
import numpy as np
        

# -----------------------------------------------------------------------------------------------------------------
#                                               Debut du code
# -----------------------------------------------------------------------------------------------------------------
time_start = time.time()
# Functions
def read_data_from_text_file(filename):
    data = {}
    with open(filename, 'r') as file:
        # Read the header line
        header = file.readline().strip().split('\t')
        
        # Initialize dictionary with empty lists for each variable
        for var in header:
            data[var] = []
        
        # Read data line by line
        for line in file:
            values = line.strip().split('\t')
            for i, var in enumerate(header):
                data[var].append(float(values[i]))
    
    return data

# Example usage:
filename = 'src/tempoutputfile.dat'
data = read_data_from_text_file(filename)
print(data)


# -----------------------------------------------------------------------------------------------------------------
#                           Solution avec le maillage le plus fin (pour MNP)
# -----------------------------------------------------------------------------------------------------------------


# Graphique log-log norme de l'erreur L1 vs delta_r
outputFolder="results"
array_x_conv = data['delta_x']
print(array_x_conv)
erreur_vect_L1 = data['k']
print(erreur_vect_L1)
print(erreur_vect_L1[0])
for i in range(0,len(erreur_vect_L1)):
    erreur_vect_L1[i] = np.abs(erreur_vect_L1[i]-erreur_vect_L1[-1])
print(erreur_vect_L1)

plt.figure(1)
plt.loglog(array_x_conv[:-2], erreur_vect_L1[:-2], '.r', label = "Norme L1")

# Ajuster une loi de puissance à toutes les valeurs (en utilisant np.polyfit avec logarithmes)
coefficients = np.polyfit(np.log(array_x_conv[0:-2]), np.log(erreur_vect_L1[0:-2]), 1)
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
equation_text = f'$L_1 = {np.exp(constant_logreg):.4E} \\times Δr^{{{exponent_logreg:.4f}}}$'
equation_text_obj = plt.text(0.05, 0.05, equation_text, fontsize=12, transform=plt.gca().transAxes, color='k')

plt.xlabel("delta_x")
plt.ylabel("erreur")
plt.legend()
plt.grid()
plt.title("Normes de l'écart avec maillage le plus fin en fonction de Δx")
plt.savefig(outputFolder+"etude_convergence_L1.png")

plt.show()

time_end = time.time()
print("Time end:")
print(time_end-time_start)
