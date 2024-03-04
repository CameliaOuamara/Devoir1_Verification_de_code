# -*- coding: utf-8 -*-
"""
Created on Sat Mar  2 21:04:08 2024

@author: ouama
"""

# Libraries
import numpy as np
import sympy as sy

# Functions
# from profil_concentration import *
# from plot import *
# from norme_erreur_discretisation import *
# from etude_convergence import *

Ce = 12
R = 0.5
D = 10**-10
k = 4*10**-9

r, t = sy.symbols('r t')

# Solution mms
C_MMS = Ce * sy.sin(t) * sy.sin(sy.pi*r/R)

# Derivées
dCdt   = sy.diff(C_MMS, t)
dCdr   = sy.diff(C_MMS, r)
d2Cdr2 = sy.diff(sy.diff(C_MMS, r), r)

# Terme source
S_MMS = D*d2Cdr2 + D*dCdr/r - k*C_MMS - dCdt

# Callable functions
C =  sy.lambdify([r, t], C_MMS, "numpy")
S =  sy.lambdify([r, t], S_MMS, "numpy")

























































