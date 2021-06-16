import pytest

from pysisyphus.helpers import geom_loader
from pysisyphus.optimizers import MicroOptimizer
from pysisyphus.testing import using
from pysisyphus.calculators.PySCF import PySCF


def get_coords(atoms, cart_coords):
    mod_coords = cart_coords.reshape(-1, 3).copy()
    mod_coords[:, 0] += 1
    results = {"coords": mod_coords}
    return results


@using("pyscf")
@pytest.mark.parametrize(
    "coord_type",
    [
        "cart",
        "redund",
    ],
)
@pytest.mark.parametrize(
    "step",
    [
        "sd",
        "cg",
        "lbfgs",
    ],
)
def test_micro_optimizer(coord_type, step):
    geom = geom_loader("lib:methane.xyz", coord_type=coord_type)
    calc = PySCF(basis="sto3g", verbose=0)
    # Patch calculator so it can reparametrize
    calc.get_coords = get_coords
    geom.set_calculator(calc)
    opt = MicroOptimizer(geom, step=step, rms_force=1e-5)
    opt.run()

    assert opt.is_converged
    assert opt.cur_cycle in (3, 4)  # Depends on step and coord_type...
