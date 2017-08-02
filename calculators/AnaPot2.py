import numpy as np

from pysisyphus.calculators.Calculator import Calculator

# https://www.wolframalpha.com/input/?i=plot+arccot(-exp(y)*cot(x/2-pi/4))+-+2*exp(-(y-sin(x))^2/2)
# [1] http://aip.scitation.org/doi/abs/10.1063/1.461606
# https://www.wolframalpha.com/input/?i=derivative+of+(arccot(-exp(y)*cot(x/2-pi/4))+-+2*exp(-(y-sin(x))^2/2))

class AnaPot2(Calculator):

    def __init__(self): 
        super(AnaPot2, self).__init__()

    def get_energy(self, atoms, coords):
        x, y, z = coords
        cot = 1 / np.tan(x/2 - np.pi/4)
        arccot = np.arctan(
                1 / (-np.exp(y) * cot)
        )
        energy = (
            arccot - 2 * np.exp(-(y-np.sin(x))**2 / 2)
        )
        return {"energy": energy}

    def __str__(self):
        return "AnaPot2 calculator"


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    xlim = (-np.pi/2, np.pi)
    ylim = (-2, 2)
    x = np.linspace(*xlim, 100)
    y = np.linspace(*ylim, 100)
    X, Y = np.meshgrid(x, y)
    Z = np.full_like(X, 0)
    fake_atoms = ("H", )
    pot_coords = np.stack((X, Y, Z))
    pot = AnaPot2().get_energy(fake_atoms, pot_coords)["energy"]

    levels=(-2, 1, 40)
    levels = np.linspace(*levels)
    fig, ax = plt.subplots()
    contours = ax.contour(X, Y, pot, levels)
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    fig.colorbar(contours, cax=cbar_ax)
    plt.show()
