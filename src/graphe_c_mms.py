import numpy as np
import matplotlib.pyplot as plt
import sympy

class YourClass:
    def __init__(self, R):
        self.R = R
        r_MMS, t_MMS, Deff, k = sympy.symbols('r_MMS t_MMS Deff k')
        self.C_MMS = sympy.cos(sympy.pi*r_MMS/(2.*self.R))*sympy.exp(-t_MMS) + sympy.cos(sympy.pi*r_MMS/(2.*self.R))
        self.f_T_MMS = sympy.lambdify([r_MMS, t_MMS], self.C_MMS, "numpy")
        # Appliquer l'opérateur sur la solution MMS
        self.source = sympy.diff(self.C_MMS,t_MMS) - Deff * (sympy.diff(self.C_MMS,r_MMS)/r_MMS + sympy.diff(sympy.diff(self.C_MMS,r_MMS),r_MMS)) + k*self.C_MMS
        print("self.source: ", self.source)
        # print("self.source: ", self.source)
        # create callable function for symbolic expression
        self.f_source = sympy.lambdify([r_MMS,t_MMS], self.source, "numpy")
        
    def plot_function(self):
        t_values = np.arange(0, 2.1, 0.1)
        r_values = np.linspace(0, 0.5, 200)

        fig, ax = plt.subplots()
        
        for i, t_MMS in enumerate(t_values):
            color = plt.cm.viridis(i / len(t_values))  # Gradual color change
            ax.plot(r_values, self.f_T_MMS(r_values, t_MMS), color=color, label=f't_MMS={t_MMS:.1f}')

        ax.set_title('$C_{MMS}(r,t)$ en fonction de $r$ pour différent $t$')
        ax.set_xlabel('$r$')
        ax.set_ylabel('$C_{MMS}(r,t)$')
        # ax.legend()
        ax.grid(True)  # Add grid

        plt.savefig('c_mms.png',dpi=700)  # Save plot as a .png file
        plt.show()

# Example usage:
your_instance = YourClass(R=0.5)
your_instance.plot_function()
