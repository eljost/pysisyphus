#!/usr/bin/env python3


from pysisyphus.helpers import geom_from_library
from pysisyphus.calculators.PySCF import PySCF

import numpy as np
from pytest import approx


def get_calc():
    return PySCF(basis="321g", pal=2)


def test_inthess():
    geom = geom_from_library("h2o2_rot2.xyz", coord_type="redund")
    geom.set_calculator(get_calc())
    ref_hess = geom.hessian
    # np.savetxt("ref_hess", ref_hess)
    print("ref_hess")
    print(ref_hess)
    print("forces")
    print(geom.forces)

    int_ = geom.internal
    coords0 = geom.coords.copy()
    step_size = 0.001
    num_hess = list()
    num_hess2 = list()
    for i in range(geom.coords.size):
        step = np.zeros_like(coords0)
        step[i] = step_size
        geom_ = geom.copy()
        geom_.set_calculator(get_calc())
        geom_.coords = coords0 + step
        plus_grad = geom_.gradient
        geom_.coords = coords0 - step
        minus_grad = geom_.gradient
        # Using one step size
        fd = (plus_grad - minus_grad) / (2*step_size)
        print(f"{i:02d}:")
        print("\tone step")
        print("\t", fd)
        num_hess.append(fd)

        # # Using two step sizes
        # step2 = step*2
        # geom_.coords = coords0 + step2
        # plus2_grad = geom_.gradient
        # geom_.coords = coords0 - step2
        # minus2_grad = geom_.gradient
        # fd2 = (2/3 * (plus_grad - minus_grad) - 1/12*(plus2_grad - minus2_grad))/step_size
        # print("\ttwo steps")
        # print("\t", fd2)
        # num_hess2.append(fd2)

    num_hess = np.array(num_hess)
    # num_hess2 = np.array(num_hess2)
    # Symmetrize
    num_hess = (num_hess + num_hess.T) / 2
    # num_hess2 = (num_hess2 + num_hess2.T) / 2
    print("one step")
    print(num_hess)
    # np.savetxt("num_hess", num_hess)
    # print("two steps")
    # print(num_hess2)
    # np.savetxt("num_hess2", num_hess2)

    diff = ref_hess - num_hess
    print("diff")
    print(diff)
    ds = np.abs(diff).sum()
    print(f"sum(abs(diff)) {ds:.8f}")
    rms = np.sqrt(np.mean(diff**2))
    print(f"rms {rms:.8f}")
    assert rms == approx(0.00077671)

    # diff2 = ref_hess - num_hess2
    # print("diff2")
    # print(diff2)
    # ds2 = np.abs(diff2).sum()
    # print(f"sum(abs(diff2)) {ds2:.8f}")
    # rms2 = np.sqrt(np.mean(diff2**2))
    # print(f"rms2 {rms2:.8f}")


if __name__ == "__main__":
    test_inthess()
