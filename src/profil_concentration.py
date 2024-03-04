"""
Classe qui calcule la solution numerique du probleme étudié pour une discretisation geometrique donnée
"""

# Libraries
import numpy as np

# import sympy as sp
# import scipy.sparse as sp
# import scipy.sparse.linalg
import scipy as sp
import sympy
from scipy.linalg import svd

class Profil_Concentration:

    def __init__(self, delta_r, delta_t, N, R, critere_conv,critere_max_iter):
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
        
        # Données :
        self.Ce = 12 # [mol/m3]
        self.Deff = 10**-10 # [m2/s]
        self.S = 8*10**-9 # [mol/m3/s]
        self.k = 4.0e-9
        
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
        self.A = sp.sparse.lil_matrix((self.N,self.N))
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
        self.A_inverse = sp.sparse.linalg.inv(self.A)
            
        
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
        
        self.B = sp.sparse.lil_matrix((self.N, 1))
        
        self.B[0,0]  = 0.0
        self.B[-1,0] = self.Ce
        
        for i in range(1,self.N-1):
            self.B[i,0] = - self.e * C_t[i]
            
        self.B = self.B.tocsc()
            
    def Algorithme_Resolution(self):
        
        # Calcul matrice A
        self.Matrice_A()
        
        # Initialisation du temps
        self.t = 0 + self.Delta_t
        i = 0
    
        # Concentration au temps t0
        C_t = np.zeros((1, self.N))
        C_t_plus_1 = np.zeros((1, self.N))
        
        # Initailisation matrice concentration a chaque temps
        self.C = np.zeros((1, self.N))
        # self.C[0,:] = C_t
            
        #initialisation de diff_temporelle pour s'assurer qu'on soit steady à la dernière itération
        diff_temporelle = 1.0e10
        
        while diff_temporelle>=self.critere_conv and i < self.critere_max_iter:
            # Construction matrice B
            self.Matrice_B(self.C[i])
            
            # Resolution du systeme matriciel
            # C_t_plus_1[0,:] = sp.sparse.linalg.spsolve(self.A, self.B)
            C_t_plus_1 = self.A_inverse.dot(self.B).toarray()

            # Remplissage de la matrice de concentration
            self.C = np.append(self.C, C_t_plus_1.T, axis=0)
            
            # Avancer au temps suivant
            self.t += self.Delta_t
            i+= 1
            
            
            # calcul de la différence temporelle
            diff_temporelle = np.linalg.norm(self.C[i, :] - self.C[i-1, :])/np.linalg.norm(self.C[i, :])
        
        # for i in range(self.critere_max_iter):
            # # Construction matrice B
            # self.Matrice_B(self.C[i])
            
            # # Resolution du systeme matriciel
            # # C_t_plus_1[0,:] = sp.sparse.linalg.spsolve(self.A, self.B)
            # C_t_plus_1 = self.A_inverse.dot(self.B).toarray()

            # # Remplissage de la matrice de concentration
            # self.C = np.append(self.C, C_t_plus_1.T, axis=0)
            
            # # Avancer au temps suivant
            # self.t += self.Delta_t            

# Nouvelle classe qui est pareille que Profil_Concentration sauf le schéma utilisé pour la dérivé dC/dr qui est centrée plutôt qu'avant          
class Profil_Concentration_Centree(Profil_Concentration):
    def Matrice_A(self):
        """

        Returns
        -------
        Membre de gauche. Matrice A des coefficients

        """
        self.A = sp.sparse.lil_matrix((self.N,self.N))
        self.A_inverse = np.zeros((self.N,self.N))
        
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
        self.A_inverse = sp.sparse.linalg.inv(self.A)
        
        
class Profil_Concentration_Centree_MNP(Profil_Concentration_Centree):
    def __init__(self, delta_r, delta_t, N, R, critere_conv,critere_max_iter, spline_bicubic):
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
        
        # Données :
        self.Ce = 12 # [mol/m3]
        self.Deff = 10**-10 # [m2/s]
        self.S = 8*10**-9 # [mol/m3/s]
        self.k = 4.0e-9
        
        self.a = self.Deff/(self.Delta_r**2)
        self.b = self.Deff/self.Delta_r
        self.e = 1/self.Delta_t
        self.r = np.linspace(0, self.R, self.N)   
        self.spline_bicubic = spline_bicubic
        
            
                
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
        dc_dt = self.spline_bicubic.partial_derivative(1,0)
        dc_dr = self.spline_bicubic.partial_derivative(0,1)
        d2c_dr2 = self.spline_bicubic.partial_derivative(0,2)
        c = self.spline_bicubic

        # S_MNP = dc_dt(self.t,self.r[0]) - self.Deff*(1/self.r[0] * dc_dr(self.t,self.r[0]) + d2c_dr2(self.t,self.r[0])) + self.k * c(self.t,self.r[0])
        
        self.B = sp.sparse.lil_matrix((self.N, 1))
        
        self.B[0,0]  = 0.0
        self.B[-1,0] = self.Ce
        
        for i in range(1,self.N-1):
            S_MNP = dc_dt(self.t,self.r[i]) - self.Deff*(1/self.r[i] * dc_dr(self.t,self.r[i]) + d2c_dr2(self.t,self.r[i])) + self.k * c(self.t,self.r[i])
            # S_MNP = -self.Deff*(1/self.r[i] * dc_dr(self.t,self.r[i]) + d2c_dr2(self.t,self.r[i])) + self.k * c(self.t,self.r[i])
            self.B[i,0] = - self.e * C_t[i] + S_MNP           
            # self.B[i,0] = - self.e * C_t[i]         
            
                
class Profil_Concentration_MNP(Profil_Concentration):
    def __init__(self, delta_r, delta_t, N, R, critere_conv,critere_max_iter, spline_bicubic):
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
        
        # Données :
        self.Ce = 12 # [mol/m3]
        self.Deff = 10**-10 # [m2/s]
        self.S = 8*10**-9 # [mol/m3/s]
        self.k = 4.0e-9
        
        self.a = self.Deff/(self.Delta_r**2)
        self.b = self.Deff/self.Delta_r
        self.e = 1/self.Delta_t
        self.r = np.linspace(0, self.R, self.N)   
        self.spline_bicubic = spline_bicubic
        
            
                
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
        dc_dt = self.spline_bicubic.partial_derivative(1,0)
        dc_dr = self.spline_bicubic.partial_derivative(0,1)
        d2c_dr2 = self.spline_bicubic.partial_derivative(0,2)
        c = self.spline_bicubic

        # S_MNP = dc_dt(self.t,self.r[0]) - self.Deff*(1/self.r[0] * dc_dr(self.t,self.r[0]) + d2c_dr2(self.t,self.r[0])) + self.k * c(self.t,self.r[0])
        
        self.B = sp.sparse.lil_matrix((self.N, 1))
        
        self.B[0,0]  = 0.0
        self.B[-1,0] = self.Ce
        
        for i in range(1,self.N-1):
            S_MNP = dc_dt(self.t,self.r[i]) - self.Deff*(1/self.r[i] * dc_dr(self.t,self.r[i]) + d2c_dr2(self.t,self.r[i])) + self.k * c(self.t,self.r[i])
            # S_MNP = -self.Deff*(1/self.r[i] * dc_dr(self.t,self.r[i]) + d2c_dr2(self.t,self.r[i])) + self.k * c(self.t,self.r[i])
            self.B[i,0] = - self.e * C_t[i] + S_MNP           
            # self.B[i,0] = - self.e * C_t[i]  
            
            
class Profil_Concentration_MMS(Profil_Concentration):
    def __init__(self, delta_r, delta_t, N, R, t_final):
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
        self.t_final = t_final
        self.N = N
        self.R = R
        
        # Données :
        self.Ce = 12 # [mol/m3]
        self.Deff = 10**-10 # [m2/s]
        self.S = 8*10**-9 # [mol/m3/s]
        self.k = 4.0e-9
        
        self.a = self.Deff/(self.Delta_r**2)
        self.b = self.Deff/self.Delta_r
        self.e = 1/self.Delta_t
        self.r = np.linspace(0, self.R, self.N)
        
        r_MMS,t_MMS = sympy.symbols('r_MMS t_MMS')

        #solution MMS
        self.C_MMS = sympy.cos(sympy.pi*r_MMS/(2.*self.R))*sympy.exp(-t_MMS) + sympy.cos(sympy.pi*r_MMS/(2.*self.R))
        # self.C_MMS = sympy.cos(sympy.pi*r_MMS/(2.*self.R))/(t_MMS+1.0) + sympy.cos(sympy.pi*r_MMS/(2.*self.R))
        # create callable function for symbolic expression
        self.f_T_MMS = sympy.lambdify([r_MMS,t_MMS], self.C_MMS, "numpy")

        # Appliquer l'opérateur sur la solution MMS
        self.source = sympy.diff(self.C_MMS,t_MMS) - self.Deff * (sympy.diff(self.C_MMS,r_MMS)/r_MMS + sympy.diff(sympy.diff(self.C_MMS,r_MMS),r_MMS)) + self.k*self.C_MMS
        
        # print("self.source: ", self.source)
        # create callable function for symbolic expression
        self.f_source = sympy.lambdify([r_MMS,t_MMS], self.source, "numpy")
        
        print("Profil_Concentration_MMS")
        

    def Matrice_A(self):
        """

        Returns
        -------
        Membre de gauche. Matrice A des coefficients

        """
        self.A = sp.sparse.lil_matrix((self.N,self.N))
        self.A_inverse = np.zeros((self.N,self.N))

        self.A[0,0]   = -1
        self.A[0,1]   =  1
        # self.A[0,0] =  1

        
        self.A[-1,-1] =  1
        
        for i in range(1,self.N-1):
            self.A[i, i-1]  = self.a
            
            self.A[i, i]   = -2*self.a - self.b/self.r[i] - self.e - self.k

            self.A[i, i+1] = self.a + self.b/self.r[i]
            
        #storer l'inverse puisqu'elle ne change pas
        self.A = self.A.tocsr()
        self.A_inverse = sp.sparse.linalg.inv(self.A)
        # # Perform singular value decomposition
        # u, s, vh = svd(self.A.toarray())
        
        # # Calculate the condition number
        # condition_number = np.max(s) / np.min(s)
        
        # print("condition_number: ", condition_number)

                
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

        
        self.B = sp.sparse.lil_matrix((self.N, 1))
        
        # self.B[0,0]  = self.f_T_MMS(0.0,self.t)
        self.B[0,0]  = 0.0
        self.B[-1,0] = self.f_T_MMS(self.R,self.t)
        
        for i in range(1,self.N-1):
            # S_MNP = -self.Deff*(1/self.r[i] * dc_dr(self.t,self.r[i]) + d2c_dr2(self.t,self.r[i])) + self.k * c(self.t,self.r[i])
            self.B[i,0] = - self.e * C_t[i] - self.f_source(self.r[i],self.t)           
            # self.B[i,0] = - self.e * C_t[i]              


    def Algorithme_Resolution(self):
        
        # Calcul matrice A
        self.Matrice_A()
        
        # Initialisation du temps
        self.t = 0 + self.Delta_t
        i = 0
    
        # Concentration au temps t0
        C_t = self.f_T_MMS(self.r,0.0)
        C_t_plus_1 = np.zeros((1, self.N))
        
        # Initailisation matrice concentration a chaque temps
        self.C = np.zeros((1, self.N))
        self.C[0,:] = C_t
        
        while self.t < self.t_final:

            # Construction matrice B
            self.Matrice_B(self.C[i])
            
            # Resolution du systeme matriciel
            # C_t_plus_1[0,:] = sp.sparse.linalg.spsolve(self.A, self.B)
            C_t_plus_1 = self.A_inverse.dot(self.B).toarray()

            # Remplissage de la matrice de concentration
            self.C = np.append(self.C, C_t_plus_1.T, axis=0)
            
            # Avancer au temps suivant
            self.t += self.Delta_t
            i+= 1
            
        # print("self.t: ", self.t)
        # print("self.Delta_r: ", self.Delta_r)

class Profil_Concentration_Centree_MMS(Profil_Concentration_Centree):
    def __init__(self, delta_r, delta_t, N, R, t_final):
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
        self.t_final = t_final
        self.N = N
        self.R = R
        
        # Données :
        self.Ce = 12 # [mol/m3]
        self.Deff = 10**-10 # [m2/s]
        self.S = 8*10**-9 # [mol/m3/s]
        self.k = 4.0e-9
        
        self.a = self.Deff/(self.Delta_r**2)
        self.b = self.Deff/self.Delta_r
        self.e = 1/self.Delta_t
        self.r = np.linspace(0, self.R, self.N)
        
        r_MMS,t_MMS = sympy.symbols('r_MMS t_MMS')

        #solution MMS
        self.C_MMS = sympy.cos(sympy.pi*r_MMS/(2.*self.R))*sympy.exp(-t_MMS) + sympy.cos(sympy.pi*r_MMS/(2.*self.R))
        # self.C_MMS = sympy.cos(sympy.pi*r_MMS/(2.*self.R))/(t_MMS+1.0) + sympy.cos(sympy.pi*r_MMS/(2.*self.R))
        # create callable function for symbolic expression
        self.f_T_MMS = sympy.lambdify([r_MMS,t_MMS], self.C_MMS, "numpy")

        # Appliquer l'opérateur sur la solution MMS
        self.source = sympy.diff(self.C_MMS,t_MMS) - self.Deff * (sympy.diff(self.C_MMS,r_MMS)/r_MMS + sympy.diff(sympy.diff(self.C_MMS,r_MMS),r_MMS)) + self.k*self.C_MMS
        # create callable function for symbolic expression
        self.f_source = sympy.lambdify([r_MMS,t_MMS], self.source, "numpy")
        
        print("Profil_Concentration_Centree_MMS")
        
        
    def Matrice_A(self):
        """

        Returns
        -------
        Membre de gauche. Matrice A des coefficients

        """
        self.A = sp.sparse.lil_matrix((self.N,self.N))
        self.A_inverse = np.zeros((self.N,self.N))
        
        self.A[0,0]   = -3.0
        self.A[0,1]   =  4.0
        self.A[0,2]   =  -1.0
        
        # self.A[0,0] =  1
        self.A[-1,-1] =  1
        
        for i in range(1,self.N-1):
            self.A[i, i-1]  = self.a - 0.5 * self.b/self.r[i]
            
            self.A[i, i]   = -2*self.a - self.e - self.k

            self.A[i, i+1] = self.a + 0.5 * self.b/self.r[i]
            
        #storer l'inverse puisqu'elle ne change pas
        self.A = self.A.tocsr()            
        self.A_inverse = sp.sparse.linalg.inv(self.A)
        
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

        # S_MNP = dc_dt(self.t,self.r[0]) - self.Deff*(1/self.r[0] * dc_dr(self.t,self.r[0]) + d2c_dr2(self.t,self.r[0])) + self.k * c(self.t,self.r[0])
        
        self.B = sp.sparse.lil_matrix((self.N, 1))
        
        self.B[0,0]  = 0.0
        self.B[-1,0] = 0.0
        
        for i in range(1,self.N-1):
            # S_MNP = -self.Deff*(1/self.r[i] * dc_dr(self.t,self.r[i]) + d2c_dr2(self.t,self.r[i])) + self.k * c(self.t,self.r[i])
            self.B[i,0] = - self.e * C_t[i] - self.f_source(self.r[i],self.t)           
            # self.B[i,0] = - self.e * C_t[i]              


    # def Algorithme_Resolution(self):
        
    #     # Calcul matrice A
    #     self.Matrice_A()
        
    #     # Initialisation du temps
    #     self.t = 0
    #     i = 0
    
    #     # Concentration au temps t0
    #     C_t = self.f_T_MMS(self.r,0.0)
    #     C_t_plus_1 = np.zeros((1, self.N))
        
    #     # Initailisation matrice concentration a chaque temps
    #     self.C = np.zeros((1, self.N))
    #     self.C[0,:] = C_t
        
    #     while self.t < self.t_final:
    #         # Avancer au temps suivant
    #         self.t += self.Delta_t
    #         # Construction matrice B
    #         self.Matrice_B(self.C[i])
            
    #         # Resolution du systeme matriciel
    #         # C_t_plus_1[0,:] = sp.sparse.linalg.spsolve(self.A, self.B)
    #         C_t_plus_1 = self.A_inverse.dot(self.B).toarray()

    #         # Remplissage de la matrice de concentration
    #         self.C = np.append(self.C, C_t_plus_1.T, axis=0)
        
    #         i+= 1
            
    #     # print("self.t: ", self.t)
    #     print("self.Delta_r: ", self.Delta_r)

    def Algorithme_Resolution(self):
        
        # Calcul matrice A
        self.Matrice_A()
        
        # Initialisation du temps
        self.t = 0 + self.Delta_t
        i = 0
    
        # Concentration au temps t0
        C_t = self.f_T_MMS(self.r,0.0)
        C_t_plus_1 = np.zeros((1, self.N))
        
        # Initailisation matrice concentration a chaque temps
        self.C = np.zeros((1, self.N))
        self.C[0,:] = C_t
        
        while self.t < self.t_final:

            # Construction matrice B
            self.Matrice_B(self.C[i])
            
            # Resolution du systeme matriciel
            # C_t_plus_1[0,:] = sp.sparse.linalg.spsolve(self.A, self.B)
            C_t_plus_1 = self.A_inverse.dot(self.B).toarray()

            # Remplissage de la matrice de concentration
            self.C = np.append(self.C, C_t_plus_1.T, axis=0)
            
            # Avancer au temps suivant
            self.t += self.Delta_t
            i+= 1









































