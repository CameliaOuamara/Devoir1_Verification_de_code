# -*- coding: utf-8 -*-
"""
Created on Sun Mar  3 13:50:51 2024

@author: ouama
"""
import numpy as np
import sympy as sy

"""
Fichier contenant les classes utilisées dans le code
Les classes definissent les termes sources, conditions frontieres et condition initiale 
qui seront utilisées pour definir le probleme
"""

class Fick():
    def __init__(self, N):
        self.Terme_source = 0
        self.dCdr_r0 = 0
        self.C_rR = 12 # Ce
        self.C0 = np.zeros(N)
        
class MMS():
    def __init__(self):
        Ce = 12
        R = 0.5
        D = 1
        k = 4*10**-9

        r, t = sy.symbols('r t')

        # Solution mms
        C_MMS = Ce * sy.sin(t) * sy.exp(sy.pi*r/R)

        # Derivées
        dCdt   = sy.diff(C_MMS, t)
        dCdr   = sy.diff(C_MMS, r)
        d2Cdr2 = sy.diff(sy.diff(C_MMS, r), r)

        # Terme source
        S_MMS = D*d2Cdr2 + D*dCdr/r - k*C_MMS - dCdt

        # Callable functions
        C =  sy.lambdify([r, t], C_MMS, "numpy")
        S =  sy.lambdify([r, t], S_MMS, "numpy")
        
        self.Terme_source = S
        self.dCdr_r0 = sy.lambdify([r,t], dCdr, "numpy")
        self.C = C
    
    