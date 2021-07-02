import numpy as np
import pytest

from pysisyphus.calculators.PySCF import PySCF
from pysisyphus.helpers import geom_loader
from pysisyphus.intcoords.PrimTypes import (
    PrimTypes,
    normalize_prim_inputs,
    normalize_prim_input,
)
from pysisyphus.optimizers.guess_hessians import get_guess_hessian
from pysisyphus.optimizers.RFOptimizer import RFOptimizer
from pysisyphus.testing import using


def test_normalize_prim_inputs():
    # ref_prim_inp = [PrimTypes(0), 1, 2]
    ref_prim_inp = (PrimTypes.BOND, 1, 2)
    prim_inps = [["BOND", 1, 2], [0, 1, 2], [0.0, 1, 2]]
    normalize_prim_inputs(prim_inps)
    for pi in prim_inps:
        assert normalize_prim_input(pi)[0] == ref_prim_inp


def test_normalize_prim_inputs_shortcuts():
    prim_inps = [["ATOM", 0], ["XYZ", 0], ["atom", 0], ["xyz", 0]]
    ref_prim_inp = [
        (PrimTypes.CARTESIAN_X, 0),
        (PrimTypes.CARTESIAN_Y, 0),
        (PrimTypes.CARTESIAN_Z, 0),
    ]
    for pi in prim_inps:
        assert normalize_prim_input(pi) == ref_prim_inp


@pytest.mark.xfail
def test_normalize_fail():
    prim_inp = ["NOT_DEFINED!", 1]
    normalize_prim_input(prim_inp)


def h2o_fix(fixed_atoms):
    constrain_prims = [["atom", ind] for ind in fixed_atoms]
    kwargs = {
        "constrain_prims": constrain_prims,
    }
    geom = geom_loader("lib:h2o.xyz", coord_type="redund", coord_kwargs=kwargs)
    return geom


@pytest.fixture
def h2o_full():
    return h2o_fix((0, 1, 2))


@pytest.fixture
def h2o_ofix():
    return h2o_fix((0,))


@pytest.mark.parametrize("hessian_init", ["fischer", "lindh", "simple", "swart"])
def test_constrain_guess_hessians(hessian_init, h2o_ofix):
    """Checks that the different Hessian don't fail"""

    get_guess_hessian(h2o_ofix, hessian_init)


@using("pyscf")
def test_constrain_h2o_full(h2o_full):
    """Immediate convergence, as all 3 atoms are constrained."""
    start = h2o_full.coords3d.copy()
    calc = PySCF(basis="sto3g")
    h2o_full.set_calculator(calc)
    opt = RFOptimizer(h2o_full)
    opt.run()

    assert opt.is_converged
    assert opt.cur_cycle == 0
    end = h2o_full.coords3d
    np.testing.assert_allclose(end, start)


@using("pyscf")
def test_constrain_h2o_ofix(h2o_ofix):
    """Immediate convergence, as all 3 atoms are constrained."""
    start = h2o_ofix.coords3d.copy()
    calc = PySCF(basis="sto3g")
    h2o_ofix.set_calculator(calc)
    opt = RFOptimizer(h2o_ofix, thresh="gau", dump=True)
    opt.run()

    assert opt.is_converged
    assert h2o_ofix.energy == pytest.approx(-74.96590119)
    end = h2o_ofix.coords3d
    np.testing.assert_allclose(end[0], start[0], atol=1e-18)
