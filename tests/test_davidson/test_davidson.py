import numpy as np
import pytest

from pysisyphus.calculators import XTB
from pysisyphus.helpers import geom_loader
from pysisyphus.helpers_pure import eigval_to_wavenumber
from pysisyphus.modefollow import geom_davidson
from pysisyphus.modefollow.NormalMode import NormalMode
from pysisyphus.tsoptimizers.RSPRFOptimizer import RSPRFOptimizer
from pysisyphus.testing import using


@using("xtb")
@pytest.mark.parametrize(
    "precon, ref_cyc, ref_nu",
    [
        (True, 1, 1690.976973),
        (False, 3, 1690.970311),
    ],
)
def test_block_davidson_acet(precon, ref_cyc, ref_nu, this_dir):
    geom = geom_loader("lib:acet_tm.xyz")
    calc = XTB(pal=2)
    geom.set_calculator(calc)

    *_, cart_displs = geom.get_normal_modes()
    inds = [10, 2]

    rg = np.random.default_rng(20180325)

    def get_guess(vec, masses_rep, scale=5e-3):
        return NormalMode(vec + scale * rg.random(*vec.shape), masses_rep)

    guess_modes = [get_guess(cart_displs[:, ind], geom.masses_rep) for ind in inds]

    hessian_precon = None
    if precon:
        hessian_precon = np.loadtxt(this_dir / "hessian_precon")

    result = geom_davidson(
        geom,
        guess_modes,
        hessian_precon=hessian_precon,
        start_precon=5,
        print_level=1,
    )

    nu = result.nus[result.mode_inds[0]]

    assert result.converged
    assert result.cur_cycle == ref_cyc
    assert nu == pytest.approx(ref_nu)


@using("pyscf")
def test_block_davidson_hcn():
    geom = geom_loader("lib:hcn_iso_pm6_near_ts.xyz", coord_type="redund")
    calc = XTB(pal=2)
    geom.set_calculator(calc)

    opt = RSPRFOptimizer(geom, thresh="gau")
    opt.run()

    *_, cart_displs = geom.get_normal_modes()
    rg = np.random.default_rng(20180325)

    def get_guess(vec, masses_rep, scale=1e-1):
        return NormalMode(vec + scale * rg.random(*vec.shape), masses_rep)

    guess_modes = [get_guess(cart_displs[:, 0], geom.masses_rep), ]

    hessian_precon = None
    result = geom_davidson(
        geom,
        guess_modes,
        hessian_precon=None,
        start_precon=5,
        print_level=1,
    )

    nu = result.nus[result.mode_inds[0]]

    assert result.converged
    assert result.cur_cycle == 2
