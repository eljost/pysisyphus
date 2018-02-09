from pysisyphus.calculators.AnaPotBase import AnaPotBase


class AnaPot4(AnaPotBase):

    def __init__(self): 
        V_str = "-10*x**2 + 10*y**2 + 4*sin(x*y)-2*x+x**4"
        xlim = (-3.5, 3.5)
        ylim = (-3.5, 3.5)
        super().__init__(V_str=V_str, xlim=xlim, ylim=ylim)

    def __str__(self):
        return "AnaPot4 calculator"


if __name__ == "__main__":
    ap4 = AnaPot4()
    ap4.plot()
    import matplotlib.pyplot as plt
    plt.show()
