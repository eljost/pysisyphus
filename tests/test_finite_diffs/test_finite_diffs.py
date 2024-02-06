import numpy as np
import pytest

from pysisyphus.calculators.AnaPot import AnaPot
from pysisyphus.calculators import XTB
from pysisyphus.calculators.PySCF import PySCF
from pysisyphus.constants import ANG2BOHR
from pysisyphus.finite_diffs import finite_difference_hessian
from pysisyphus.helpers import geom_loader
from pysisyphus.testing import using


def assert_hessians(geom):
    hessian = geom.hessian

    def grad_func(coords):
        return -geom.get_energy_and_forces_at(coords)["forces"]

    fd_hessian = finite_difference_hessian(geom.coords, grad_func, acc=4)
    np.testing.assert_allclose(hessian, fd_hessian, atol=2.5e-4)


def test_fd_hessian():
    geom = AnaPot().get_saddles(i=0)
    assert_hessians(geom)


@using("pyscf")
def test_water_fd_hessian():
    geom = geom_loader("lib:h2o.xyz")
    calc = PySCF(basis="321g", pal=2)
    geom.set_calculator(calc)
    assert_hessians(geom)


@using("xtb")
@pytest.mark.parametrize(
    "pal",
    (1, 2),
)
def test_num_hessian(pal):
    geom = geom_loader("lib:benzene_xtb_opt.xyz")
    calc = XTB(acc=0.3, pal=pal)
    geom.set_calculator(calc)

    num_res = geom.calculator.get_num_hessian(geom.atoms, geom.cart_coords)
    num_hess = num_res["hessian"]
    hess = geom.hessian
    np.testing.assert_allclose(num_hess, hess, atol=5e-5)


from pysisyphus.calculators import ORCA


@using("orca")
@pytest.mark.skip_ci
@pytest.mark.parametrize(
    "serial",
    (False, True),
)
def test_orca_num_hessian(serial):
    geom = geom_loader("lib:benzene_xtb_opt.xyz")
    num_hess_kwargs = {
        "serial": serial,
        "step_size": 0.005 * ANG2BOHR,
    }
    calc = ORCA(
        keywords="tpss def2-svp tightscf",
        pal=8,
        # Compare against finite diff. Hessian calculated by ORCA itself
        numfreq=True,
        num_hess_kwargs=num_hess_kwargs,
    )
    geom.set_calculator(calc)

    num_res = geom.calculator.get_num_hessian(geom.atoms, geom.cart_coords)
    num_hess = num_res["hessian"]
    hess = geom.hessian
    np.testing.assert_allclose(num_hess, hess, atol=5e-4)
