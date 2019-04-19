import numpy as np

from pysisyphus.calculators.AnaPotBase import AnaPotBase

# Taken from [1] 10.1021/ct9005147
#  J. Chem. Theory Comput., 2010, 6 (4), pp 1136–1144

class AnaPotCBM(AnaPotBase):

    def __init__(self): 
        V_str = "x**4 + 4*x**2*y**2 - 2*x**2 + 2*y**2"
        xlim = (-1.25, 1.25)
        ylim = (-0.75, 0.75)
        levels = np.linspace(-1, 1.5, 80)
        super().__init__(V_str=V_str, xlim=xlim, ylim=ylim, levels=levels)

    def __str__(self):
        return "AnaPotCBM calculator"


if __name__ == "__main__":
    ap = AnaPotCBM()
    ap.plot()
    import matplotlib.pyplot as plt
    plt.show()
