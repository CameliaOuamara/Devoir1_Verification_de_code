"""
Classe qui calcule la solution numerique du probleme étudié pour une discretisation geometrique donnée
"""

# Libraries
import numpy as np

# import sympy as sp
import scipy.sparse as sp
import scipy.sparse.linalg
import sympy as sy
from Classes_Termes_Sources import *

class Profil_Concentration:

    def __init__(self, delta_r, delta_t, N, R, critere_conv,critere_max_iter, t_final, Classe, schema):
        """
        Parameters
        ----------
        delta_r : float
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
        self.critere_conv = critere_conv
        self.critere_max_iter = critere_max_iter
        self.N = N
        self.R = R
        self.t_final = t_final
        self.Classe = Classe
        self.schema = schema
        
        # Données :
        self.Ce = 12 # [mol/m3]
        # self.Deff = 10**-10 # [m2/s]
        self.Deff = 1*10**-10        # [m2/s]
        self.S = 8*10**-9 # [mol/m3/s]
        self.k = 4*10**-9 # []
        
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
        self.A = sp.lil_matrix((self.N,self.N))
        self.A_inverse = np.zeros((self.N,self.N))

        self.A[0,0]   = -1
        self.A[0,1]   =  1
        self.A[-1,-1] =  1
        
        for i in range(1,self.N-1):
            self.A[i, i-1]  = self.a
            
            self.A[i, i]   = -2*self.a - self.b/self.r[i] - self.e - self.k

            self.A[i, i+1] = self.a + self.b/self.r[i]
            
        #storer l'inverse puisqu'elle ne change pas
        self.A = self.A.tocsr()
        self.A_inverse = sp.linalg.inv(self.A)
        
            
        
    def Matrice_B(self, C_t, t_simulation):
        """
        Parameters
        ----------
        C_t : numpy array. [1,N]
            Vecteur contenant la valeur de concentration au temps t.

        Returns
        -------
        Membre de droite. Matrice des resultats B.
        """
        
        self.B = sp.lil_matrix((self.N, 1))
        #t_simulation = self.t_final
        # Recuperation du terme source et conditions limites
        if self.Classe == 'MMS':
            Classe = MMS()
            S = Classe.Terme_source
            dCdr_r0 = Classe.dCdr_r0(0, t_simulation)
            C_rR = Classe.C(self.R, t_simulation)
            
            for i in range(1,self.N-1):
                self.B[i,0] = - self.e * C_t[i] + S(self.r[i], t_simulation)
        else :
            Classe = Fick(self.N)
            S = Classe.Terme_source
            dCdr_r0 = Classe.dCdr_r0
            C_rR = Classe.C_rR 
            
            for i in range(1,self.N-1):
                self.B[i,0] = - self.e * C_t[i] + S
        
        if self.schema == 1:
            self.B[0,0]  = dCdr_r0 *self.Delta_r # 0.0 
        elif self.schema == 2:
            self.B[0,0]  = dCdr_r0 *2*self.Delta_r # 0.0
        self.B[-1,0] = C_rR     # self.Ce
        
        self.B = self.B.tocsc()
            
    def Algorithme_Resolution(self):
        
        # Calcul matrice A
        self.Matrice_A()
        # Initialisation du temps
        self.t = 0 + self.Delta_t
        i = 0
    
        # Initailisation matrice concentration a chaque temps
        if self.Classe=='MMS':
            Classe = Fick(self.N)
            C_0 = Classe.C0
        else:
            Classe = MMS()
            C_0 = Classe.C(self.r, 0)
            
        # Concentration au temps t0
        C_t = C_0
        C_t_plus_1 = C_t
        
        self.C = np.zeros((1, self.N))
        self.C[0,:] = C_0

        #initialisation de diff_temporelle pour s'assurer qu'on soit steady à la dernière itération
        diff_temporelle = 1.0e10
        
        i=0
        if self.t_final ==0:
            
            while diff_temporelle>=self.critere_conv and i < self.critere_max_iter:
                # Construction matrice B
                #self.Matrice_B(self.C[i], self.t)
                self.Matrice_B(C_t_plus_1, self.t)
                
                # Resolution du systeme matriciel
                # C_t_plus_1[0,:] = sp.linalg.spsolve(self.A, self.B)
                C_t_plus_1 = self.A_inverse.dot(self.B).toarray()
    
                # Remplissage de la matrice de concentration
                self.C = np.append(self.C, C_t_plus_1.T, axis=0)
                
                # Avancer au temps suivant
                self.t += self.Delta_t
                i+= 1
            
                # calcul de la différence temporelle
                diff_temporelle = np.linalg.norm(self.C[i, :] - self.C[i-1, :])/np.linalg.norm(self.C[i, :])
        else:
            while self.t<=self.t_final:
                # Construction matrice B
                #self.Matrice_B(self.C[i], self.t)
                self.Matrice_B(C_t_plus_1, self.t)
                
                # Resolution du systeme matriciel
                # C_t_plus_1[0,:] = sp.linalg.spsolve(self.A, self.B)
                C_t_plus_1 = self.A_inverse.dot(self.B).toarray()
    
                # Remplissage de la matrice de concentration
                self.C = np.append(self.C, C_t_plus_1.T, axis=0)
                
                # Avancer au temps suivant
                self.t += self.Delta_t
                i+= 1
            
                # calcul de la différence temporelle
                diff_temporelle = np.linalg.norm(self.C[i, :] - self.C[i-1, :])/np.linalg.norm(self.C[i, :])
                
            
            

# Nouvelle classe qui est pareille que Profil_Concentration sauf le schéma utilisé pour la dérivé dC/dr qui est centrée plutôt qu'avant          
class Profil_Concentration_Centree(Profil_Concentration):
    def Matrice_A(self):
        """

        Returns
        -------
        Membre de gauche. Matrice A des coefficients

        """
        self.A = sp.lil_matrix((self.N,self.N))
        self.A_inverse = np.zeros((self.N,self.N))

        # self.A[0,0]   = -1
        # self.A[0,1]   =  1
        
        self.A[0,0]   = -3.0
        self.A[0,1]   =  4.0
        self.A[0,2]   =  -1.0
        
        self.A[-1,-1] =  1
        
        for i in range(1,self.N-1):
            self.A[i, i-1]  = self.a - 0.5 * self.b/self.r[i]
            
            self.A[i, i]   = -2*self.a - self.e - self.k

            self.A[i, i+1] = self.a + 0.5 * self.b/self.r[i]
            
        #storer l'inverse puisqu'elle ne change pas
        self.A = self.A.tocsr()
        self.A_inverse = sp.linalg.inv(self.A)
                        
            
        
            
                
        
            
            
                
                
        
        

        
        















































