#!/usr/bin/env python3

import numpy as np
import pytest
from pytest import approx


from pysisyphus.helpers import geom_from_library, eigval_to_wavenumber
from pysisyphus.calculators.AnaPot import AnaPot
from pysisyphus.calculators.PySCF import PySCF
from pysisyphus.calculators.ORCA import ORCA
from pysisyphus.calculators.XTB import XTB


def numhess(geom, step_size=0.0001):
    coords = geom.coords
    cnum = len(coords)
    H = list()
    for i in range(cnum):
        print(f"Step {i+1}/{cnum}")
        # Plus
        step = np.zeros_like(coords)
        step[i] += step_size
        pl_res = geom.get_energy_and_forces_at(coords+step)
        # print(pl_res)
        min_res = geom.get_energy_and_forces_at(coords-step)
        # print(min_res)

        gp = -pl_res["forces"]
        gm = -min_res["forces"]
        fd = (gp - gm) / (2*step_size)
        # print(i, fd)
        H.append(fd)
    H = np.array(H)
    # Symmetrize
    H = (H+H.T)/2
    return H


# geom = geom_from_library("h2o2_hf_321g_opt.xyz")
# geom = geom_from_library("h2o.xyz")
geom = geom_from_library("hcn_bent.xyz")
# calc = XTB(pal=2)
# calc = ORCA(keywords="BP86 def2-SVP", pal=4)
calc = ORCA(keywords="hf sto-3g tightscf", pal=4)
geom.set_calculator(calc)
# geom = AnaPot.get_geom((-0.428, 0.981, 0.))
H = geom.hessian
print("ref hessian")
print(H)
nH = numhess(geom)
print("num hessian")
print(nH)
diff = np.abs(nH-H).sum()
print(diff)
rms = np.sqrt(np.mean((nH-H)**2))
print("rms", rms)


mw_nH = geom.mass_weigh_hessian(nH)
mwp_nH = geom.eckart_projection(mw_nH)
w, v = np.linalg.eigh(mwp_nH)
nus = eigval_to_wavenumber(w)
print(nus)


assert np.allclose(nH, H) 
# geom.calculator.plot(show=True)
# import matplotlib.pyplot as plt
# plt.show()
# numhes
# np.testing.assert_allclose(nH, H) 
# import pdb; pdb.set_trace()

# for i in range(1, 6):
    # step_size = 10**(-i)
    # nH = numhess(geom, step_size)
    # diff = np.abs(nH-H).sum()
    # try:
        # assert np.allclose(nH, H) 
        # print(f"{i:02d} ss={step_size:.4e} diff={diff:.4e}")
    # except AssertionError:
        # print(f"step_size={step_size:.4e} is not enough!")


def get_calc():
    # Interestingly enough the test will fail with keep_chk == True ...
    # as the RMS value will be much higher
    return PySCF(basis="321g", pal=2, keep_chk=False)


@pytest.mark.pyscf
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
