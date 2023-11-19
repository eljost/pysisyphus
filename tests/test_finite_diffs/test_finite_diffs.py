import numpy as np

from pysisyphus.calculators.AnaPot import AnaPot
from pysisyphus.calculators import XTB
from pysisyphus.calculators.PySCF import PySCF
from pysisyphus.finite_diffs import (
    finite_difference_hessian,
    finite_difference_hessian_pal,
)
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
def test_fd_hessian_pal():
    geom = geom_loader("lib:benzene_xtb_opt.xyz")

    def calc_getter(**calc_kwargs):
        # By default, XTB uses acc=0.03 for numerical Hessians
        return XTB(acc=0.3, **calc_kwargs)

    geom.set_calculator(calc_getter())
    org_name = geom.calculator.base_name
    hessian = geom.hessian

    def grad_func(calc_number, coords):
        displ_geom = geom.copy()
        displ_geom.coords = coords
        base_name = f"{org_name}_num_hess"
        calc = calc_getter(base_name=base_name, calc_number=calc_number, pal=1)
        displ_geom.set_calculator(calc)
        return -displ_geom.get_energy_and_forces_at(coords)["forces"]

    fd_hessian = finite_difference_hessian_pal(geom.coords, grad_func, acc=4)
    np.testing.assert_allclose(hessian, fd_hessian, atol=6e-5)
