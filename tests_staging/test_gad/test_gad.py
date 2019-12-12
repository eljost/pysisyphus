import matplotlib.pyplot as plt
import numpy as np

from pysisyphus.calculators.WolfeQuapp import WolfeQuapp
from pysisyphus.tsoptimizers.GAD import GAD


def test():
    x0 = (1.15, -1.4, 0.)

    geom = WolfeQuapp.get_geom(x0)

    gad_kwargs = {
        "hessian_init": "calc",
        "dt": 0.15,
    }
    gad = GAD(geom, **gad_kwargs)
    gad.run()

    coords = np.array(gad.coords)

    calc = geom.calculator
    calc.plot()
    ax = calc.ax
    ax.plot(*coords.T[:2], "o-")
    plt.show()


if __name__ == "__main__":
    test()
