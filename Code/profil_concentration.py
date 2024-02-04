"""
Fonction de resolution d'un probleme stationnaire
"""

# Libraries
import numpy as np
import sympy as sp

class Profil_Concentration:

    def __init__(self, delta_r, delta_t, t_final, N, R):
        """
        Parameters
        ----------
        delta_x : float
            Pas geometrique en [m].
        delta_t : float
            Pas de temps en [s].
        t_final : float
            Temps final de la simulation en [s].
        N : int
            Nombre de noeuds du maillage.

        Returns
        -------
        C : numpy array
            Profil de concentration en chaque noeud du maillage.
        """
        # Entrees :
        self.Delta_r = delta_r
        self.Delta_t = delta_t
        self.t_final = t_final
        self.N = N
        self.R = R
        
        # Données :
        self.Ce = 12 # [mol/m3]
        self.Deff = 10**-10 # [m2/s]
        self.S = 8*10**-9 # [mol/m3/s]
        
        self.a = self.Deff/(self.Delta_r**2)
        self.b = self.Deff/self.Delta_r
        self.e = 1/self.Delta_t
        self.r = np.linspace(0, self.R, self.N)
        
    def Matrice_A(self):
        """

        Returns
        -------
        Membre de gauche. Matrice A des coefficients

        """
        self.A = np.zeros((self.N,self.N))
        
        self.A[0,0]   = -1
        self.A[0,1]   =  1
        self.A[-1,-1] =  1
        
        for i in range(1,self.N-1):
            self.A[i, i-1]  = self.a
            
            self.A[i, i]   = -2*self.a - self.b/self.r[i] - self.e

            self.A[i, i+1] = self.a + self.b/self.r[i]
            
        
    def Matrice_B(self, C_t):
        """
        Parameters
        ----------
        C_t : numpy array. [1,N]
            Vecteur contenant la valeur de concentration au temps t.

        Returns
        -------
        Membre de droite. Matrice des resultats B.
        """
        
        self.B = np.zeros(self.N)
        
        self.B[0]  = 0
        self.B[-1] = self.Ce
        
        for i in range(1,self.N-1):
            self.B[i] = self.S - self.e * C_t[i]
            
    def Algorithme_Resolution(self):
        
        # Calcul matrice A
        self.Matrice_A()
        
        # Initialisation du temps
        t = 0 + self.Delta_t
        i = 0
        
        # Nombre de points dans le temps
        N_t = self.t_final / self.Delta_t + 1
    
        # Concentration au temps t0
        C_t = np.zeros(self.N)
        
        # Initailisation matrice concentration a chaque temps
        self.C = np.zeros((int(N_t), self.N))
        self.C[0,:] = C_t

        while t <= self.t_final:
            # Construction matrice B
            self.Matrice_B(self.C[i])
            
            # Resolution du systeme matriciel
            C_t_plus_1 = np.linalg.solve(self.A, self.B)
            
            # Remplissage de la matrice de concentration
            self.C[i+1, :] = C_t_plus_1
            
            # Avancer au temps suivant
            t += self.Delta_t
            i+= 1
            
            
            
        
            
                
        
            
            
                
                
        
        

        
        















































