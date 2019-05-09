#!/usr/bin/env python3

from pysisyphus.helpers import geom_from_library
from pysisyphus.optimizers.guess_hessians import (lindh_guess, fischer_guess,
                                                  simple_guess, swart_guess)

import matplotlib.pyplot as plt
import numpy as np


def plot_hessian(H, title=""):
    fig, ax = plt.subplots()
    im = ax.imshow(H, vmin=0, vmax=1)
    fig.colorbar(im)
    fig.suptitle(title)


def test_guess_hessians():
    geom = geom_from_library("birkholz/vitamin_c.xyz", coord_type="redund")

    H_lindh = lindh_guess(geom)
    print("Lindh")
    print(np.diag(H_lindh))
    plot_hessian(H_lindh, "Lindh")

    H_fischer = fischer_guess(geom)
    print("Fischer")
    print(np.diag(H_fischer))
    plot_hessian(H_fischer, "Fischer")

    H_simple = simple_guess(geom)
    print("Simple")
    print(np.diag(H_simple))
    plot_hessian(H_simple, "Simple 0.5/0.2/0.1")

    H_swart = swart_guess(geom)
    print("Swart")
    print(np.diag(H_swart))
    plot_hessian(H_swart, "Swart")

    plt.show()


if __name__ == "__main__":
    test_guess_hessians()
