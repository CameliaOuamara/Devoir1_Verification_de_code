"""
Classe qui calcule la solution numerique du probleme étudié pour une discretisation geometrique donnée
"""

# Libraries
import numpy as np

# import sympy as sp
import scipy.sparse as sp
import scipy.sparse.linalg
import sympy as sy

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
        
        # Solution MMS (Temporaire)
        Ce = 12
        R = 0.5
        D = 10**-10
        k = 4*10**-9

        r, t = sy.symbols('r t')

        # Solution mms
        #C_MMS = Ce * sy.sin(t) * sy.sin(sy.pi*r/R)
        C_MMS = sy.exp(sy.pi*r/R)*sy.sin(t)

        # Derivées
        dCdt   = sy.diff(C_MMS, t)
        dCdr   = sy.diff(C_MMS, r)
        d2Cdr2 = sy.diff(sy.diff(C_MMS, r), r)

        # Terme source
        S_MMS = D*d2Cdr2 + D*dCdr/r - k*C_MMS - dCdt

        # Callable functions
        C =  sy.lambdify([r, t], C_MMS, "numpy")
        S =  sy.lambdify([r, t], S_MMS, "numpy")
        dCdt = sy.lambdify([r,t], dCdt, "numpy")
        
        self.B = sp.lil_matrix((self.N, 1))
        
        self.B[0,0]  = dCdt(0, t_simulation) #0.0
        self.B[-1,0] = C(self.R, t_simulation) #self.Ce
        
        for i in range(1,self.N-1):
            self.B[i,0] = - self.e * C_t[i] + S(self.r[i], t_simulation) #- self.e * C_t[i] 
            # self.S
            
        self.B = self.B.tocsc()
            
    def Algorithme_Resolution(self):
        
        # Calcul matrice A
        self.Matrice_A()
        
        # Initialisation du temps
        self.t = 0 + self.Delta_t
        i = 0
    
        # # Concentration au temps t0
        # C_t = np.zeros((1, self.N))
        # C_t_plus_1 = np.zeros((1, self.N))
        
        # Solution MMS (Temporaire)
        Ce = 12
        R = 0.5
        D = 10**-10
        k = 4*10**-9
        r, t = sy.symbols('r t')
        
        # Solution mms
        C_MMS = sy.exp(sy.pi*r/R)*sy.sin(t)
        C =  sy.lambdify([r, t], C_MMS, "numpy")
        C_vec = np.zeros(len(self.r))
        for i in range(len(C_vec)):
            C_vec[i] = C(self.r[i], 0)
        
        C_t = C(self.r, 0)
        C_t_plus_1 = C(self.r, 0)
        
        # Initailisation matrice concentration a chaque temps
        #self.C = np.zeros((1, self.N))
        # self.C[0,:] = C_t
        
        self.C = np.zeros((1, self.N))
        self.C[0,:] = C_vec
        

        # while t <= self.t_final:
        #     # Construction matrice B
        #     self.Matrice_B(self.C[i])
            
        #     # Resolution du systeme matriciel
        #     C_t_plus_1 = np.linalg.solve(self.A, self.B)
            
        #     # Remplissage de la matrice de concentration
        #     self.C[i+1, :] = C_t_plus_1
            
        #     # Avancer au temps suivant
        #     t += self.Delta_t
        #     i+= 1
            
        #initialisation de diff_temporelle pour s'assurer qu'on soit steady à la dernière itération
        diff_temporelle = 1.0e10
        i=0
        while diff_temporelle>=self.critere_conv and i < self.critere_max_iter:
            # Construction matrice B
            self.Matrice_B(self.C[i], self.t)
            
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
                        
            
        
            
                
        
            
            
                
                
        
        

        
        















































