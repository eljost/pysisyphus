import numpy as np
from sympy import symbols, exp

from pysisyphus.calculators.AnaPotBase import AnaPotBase

# [1] G. Henkelman, G. Jóhannesson, and H. Jónsson
# Theoretical Methods in Condensed Phase Chemistry
# (Springer, Netherlands, 2002), pp. 287

class LEPSBase(AnaPotBase):

    def __init__(self): 
        self.abc = (0.05, 0.80, 0.05)
        self.ds = (4.746, 4.746, 3.445)
        self.rac = 3.473
        self.r0 = 0.742
        self.alpha = 1.942
        # Used for LEPS + harmonic oscillator
        self.kc = 0.2025
        self.c_ = 1.154

        x, y = symbols("x y")

        V_str = str(self.V_LEPS(self.abc, self.ds, self.alpha, self.r0,
                                (x, y, self.rac))) 
        # xlim = (0, 6)
        # ylim = (0, 4)
        # levels = np.linspace(-5, 5, 250)

        # V_str = str(self.V_harmonic())
        # levels = np.linspace(-20, 20, 250)

        V_str = str(self.V_tot())
        levels = np.linspace(-5, 2, 75)

        xlim = (0.4, 3.2)
        ylim = (-2, 2)

        super().__init__(V_str=V_str, xlim=xlim, ylim=ylim, levels=levels)

    def Q(self, d, alpha, r0, r):
        """Coulomb interactions."""
        return d/2* (3/2*exp(-2*alpha*(r-r0)) - exp(-alpha*(r-r0)))

    def J(self, d, alpha, r0, r):
        """Quantum mechanical exchange interaction."""
        return d/4 * (exp(-2*alpha*(r-r0))-6*exp(-alpha*(r-r0)))

    def G(self, a, b):
        "Gaussian function."""
        return exp(-0.5*((a/0.1)**2 +(b/0.35)**2))

    def V_LEPS(self, abc, ds, alpha, r0, rs):
        """Equation (A.1) in [1].
        Mimics reaction involving three atoms confined to motion along
        a line."""
        a, b, c = abc
        dab, dbc, dac = ds
        rab, rbc, rac = rs
        rs = (rab, rbc, rac)
        
        Qab, Qbc, Qac = [self.Q(d, alpha, r0, r) for d, r in zip(ds, rs)]
        Jab, Jbc, Jac = [self.J(d, alpha, r0, r) for d, r in zip(ds, rs)]
        
        first_term = Qab/(1+a) + Qbc/(1+b) + Qac/(1+c)
        second_term = Jab**2/(1+a)**2 + Jbc**2/(1+b)**2 + Jac**2/(1+c)**2
        third_term = (-Jab*Jbc / ((1+a)*(1+b))
                      -Jbc*Jac / ((1+b)*(1+c))
                      -Jab*Jac / ((1+a)*(1+c))
        )
        return first_term - (second_term + third_term)**0.5

    def V_harmonic(self):
        """Equation (A.2) in [1].
        A and C are fixed, only B can move. A condensed phase environment
        is represented by adding a harmonic oscillator degree of freedom."""
        x, y = symbols("x y")
        return (self.V_LEPS(self.abc, self.ds, self.alpha, self.r0,
                            (x, self.rac-x, self.rac))
               + 2*self.kc*(x-(self.rac/2 - y/self.c_))**2)

    def V_tot(self):
        """Equation (A.3) in [1].
        Additional saddle point."""
        x, y = symbols("x y")
        return self.V_harmonic() + 1.5 * self.G(x-2.02083, y+0.272881)


    def __str__(self):
        return "LEPSBase calculator"


if __name__ == "__main__":
    lp = LEPSBase()
    lp.plot()
    import matplotlib.pyplot as plt
    plt.show()
