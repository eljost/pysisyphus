import numpy as np
from sympy import symbols, exp

# [1]   10.1142/9789812839664_0016
#       G. Henkelman, G. Jóhannesson, and H. Jónsson
#       Theoretical Methods in Condensed Phase Chemistry
#       (Springer, Netherlands, 2002), pp. 287
# [2]   https://aip.scitation.org/doi/pdf/10.1063/1.480097?class=pdf
#       I'm very sure that Fig. 6 and Fig. 5 are not the potentials that are
#       described in the paper. There are two errors in the paper. A_1 is
#       actually -1.5 and not 1.5 and s_x_2 is 0.5 and not 5.0. You also
#       have to divide by (2*s_x_i**2) and not only by (2*s_x_i) as given
#       in Eq. (22).

class LEPSExpr:


    def __init__(self):
        """Generates sympy expression for various LEPS potentials."""
        # V_harmonic uses b = 0.80
        self.abc = (0.05, 0.30, 0.05)
        self.ds = (4.746, 4.746, 3.445)
        self.r0 = 0.742
        self.alpha = 1.942

        # Used for LEPS + harmonic oscillator
        self.rac = 3.742
        self.kc = 0.2025
        self.c_ = 1.154

        self.choices = {
            "leps": self.get_leps,
            "harmonic": self.get_harmonic,
            "tot": self.get_tot,
            "dimer": self.get_dimer,
        }

    def get_leps(self):
        V_expr = self.V_LEPS()
        xlim = (0, 6)
        ylim = (0, 12)
        levels = np.linspace(-5, 5, 250)
        return V_expr, xlim, ylim, levels

    def get_harmonic(self):
        V_expr = self.V_harmonic()
        levels = np.linspace(-20, 20, 250)
        xlim = (0.4, 3.2)
        ylim = (-2, 2)
        return V_expr, xlim, ylim, levels

    def get_tot(self):
        V_expr = self.V_tot()
        levels = np.linspace(-5, 2, 75)
        xlim = (0.4, 3.2)
        ylim = (-2, 2)
        return V_expr, xlim, ylim, levels

    def get_dimer(self):
        V_expr = self.V_dimer()
        levels = None
        xlim = (0.25, 3.25)
        ylim = (-3.5, 3.5)
        return V_expr, xlim, ylim, levels

    def get_expr(self, pot_type="leps"):
        assert pot_type in self.choices
        return self.choices[pot_type]()


    def Q(self, d, alpha, r0, r):
        """Coulomb interactions."""
        return d/2* (3/2*exp(-2*alpha*(r-r0)) - exp(-alpha*(r-r0)))

    def J(self, d, alpha, r0, r):
        """Quantum mechanical exchange interaction."""
        return d/4 * (exp(-2*alpha*(r-r0))-6*exp(-alpha*(r-r0)))

    def G(self, a, b):
        "Gaussian function."""
        return exp(-0.5*((a/0.1)**2 +(b/0.35)**2))

    def Gdimer(self, x, y, A, x0, y0, sx, sy):
        return A * exp(-(x-x0)**2/(2*sx**2)) * exp(-(y-y0)**2/(2*sy**2))

    def V_LEPS(self, x=None, y=None, abc=None):
        """Equation (A.1) in [1].
        Mimics reaction involving three atoms confined to motion along
        a line."""
        if abc is None:
            abc = self.abc
        a, b, c = abc
        if x is None:
            x = symbols("x")
        if y is None:
            y = symbols("y")
        rs = (x, y, self.rac)
        
        Qab, Qbc, Qac = [self.Q(d, self.alpha, self.r0, r)
                         for d, r in zip(self.ds, rs)]
        Jab, Jbc, Jac = [self.J(d, self.alpha, self.r0, r)
                         for d, r in zip(self.ds, rs)]
        
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
        abc = (0.05, 0.80, 0.05)
        return (self.V_LEPS(x, self.rac-x, abc)
                + 2*self.kc*(x-(self.rac/2 - y/self.c_))**2
        )

    def V_tot(self):
        """Equation (A.3) in [1].
        Additional saddle point."""
        x, y = symbols("x y")
        return self.V_harmonic() + 1.5 * self.G(x-2.02083, y+0.272881)

    def V_dimer(self):
        """III. Results Section A in [3]. Two additional saddle points
        from two added gaussians."""
        x, y = symbols("x y")
        #           Ai , x0i    ,  y0i     , sxi, syi
        params = ((-1.5, 2.02083, -0.172881, 0.1, 0.35),
                  ( 6.0, 0.8    ,  2.0     , 0.5, 0.7 ))
        G0, G1 = [self.Gdimer(x, y, *parms) for parms in params]
        return self.V_harmonic() + G0 + G1

    def __str__(self):
        return "LEPSExpr"
