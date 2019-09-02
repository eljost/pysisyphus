#!/usr/bin/env python3

import numpy as np

from pysisyphus.calculators import MOPAC
from pysisyphus.calculators import Gaussian16
from pysisyphus.helpers import geom_from_library


def test_mopac():
    geom = geom_from_library("h2o.xyz")
    geom2 = geom.copy()
    calc = MOPAC()
    # calc = MOPAC(mult=3)
    geom.set_calculator(calc)
    grad = geom.gradient
    energy = geom.energy
    print("mop, energy", energy)
    print("mop, gradient", grad)
    hessian = geom.hessian

    g16calc = Gaussian16("PM7")
    geom2.set_calculator(g16calc)
    ref_gradient = geom2.gradient
    ref_energy = geom2.energy
    ref_hessian = geom2.hessian

    print("g16, energy", ref_energy)
    print("g16, gradient", ref_gradient)
    np.testing.assert_allclose(energy, ref_energy, atol=2e-3)
    np.testing.assert_allclose(grad, ref_gradient, atol=1e-4)
    np.testing.assert_allclose(hessian, ref_hessian, atol=1e-3)


if __name__ == "__main__":
    test_mopac()
