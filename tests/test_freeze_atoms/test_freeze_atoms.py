import numpy as np
import pytest

from pysisyphus.calculators.PySCF import PySCF
from pysisyphus.helpers import geom_loader
from pysisyphus.run import run_from_dict
from pysisyphus.testing import using
from pysisyphus.optimizers.RFOptimizer import RFOptimizer


def run_opt(geom):
    geom.set_calculator(PySCF(pal=1, basis="sto3g", verbose=0))
    opt = RFOptimizer(geom)
    opt.run()

    return geom


@using("pyscf")
@pytest.mark.parametrize("coord_type", ("cart", "redund"))
@pytest.mark.parametrize("ind", range(3))
def test_freeze_one_atom(coord_type, ind):
    """Freeze every atom once."""
    geom = geom_loader("lib:h2o.xyz", freeze_atoms=[ind], coord_type=coord_type)
    org_coords = geom.coords3d[ind].copy()
    geom = run_opt(geom)

    np.testing.assert_allclose(geom.coords3d[ind], org_coords)
    assert geom.energy == pytest.approx(-74.96589666)


@using("pyscf")
def test_freeze_two_atom():
    inds = [0, 1]
    geom = geom_loader("lib:h2o.xyz", freeze_atoms=inds)
    org_coords = geom.coords3d[inds].copy()
    geom = run_opt(geom)

    np.testing.assert_allclose(geom.coords3d[inds], org_coords)
    # Energy is higher, compared to the test with only 1 frozen atom.
    assert geom.energy == pytest.approx(-74.96484452)


@using("pyscf")
def test_run_dict_freeze():
    run_dict = {
        "geom": {
            "type": "cart",
            "fn": "lib:h2o.xyz",
            "freeze_atoms": [0, 1],
        },
        "calc": {
            "type": "pyscf",
            "basis": "sto3g",
            "pal": 1,
            "verbose": 0,
        },
        "opt": {},
    }
    result = run_from_dict(run_dict)

    geom = result.opt_geom
    assert geom.energy == pytest.approx(-74.96484452)


@using("pyscf")
def test_frozen_atoms_get_forces_at():
    ind = 0
    geom = geom_loader("lib:h2o.xyz", freeze_atoms=[ind])
    geom.set_calculator(PySCF(pal=1, basis="sto3g", verbose=0))

    funcs = (geom.get_energy_and_cart_forces_at, geom.get_energy_and_forces_at)
    for func in funcs:
        results = func(geom.coords)
        forces = results["forces"]
        np.testing.assert_allclose(forces.reshape(-1, 3)[ind], np.zeros(3))
