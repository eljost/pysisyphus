from pysisyphus.calculators.AnaPotBase import AnaPotBase

# [1] https://aip.scitation.org/doi/abs/10.1063/1.442352
#     Cerjan, Miller


class CerjanMiller(AnaPotBase):

    def __init__(self, a=1, b=1, c=1): 
        # Eq. (3.1) in [1]
        V_str = f"({a}-{b}*y**2)*x**2*exp(-x**2)+{c}/2*y**2"
        xlim = (-1.3, 1.3)
        ylim = (-0.7, 1.9)
        super().__init__(V_str=V_str, xlim=xlim, ylim=ylim)

    def __str__(self):
        return "CerjanMiller calculator"


if __name__ == "__main__":
    cj = CerjanMiller()
    cj.plot()
    import matplotlib.pyplot as plt
    plt.show()
