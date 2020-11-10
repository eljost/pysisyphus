import numpy as np
import pytest

from pysisyphus.calculators import XTB
from pysisyphus.helpers import geom_loader
from pysisyphus.helpers_pure import eigval_to_wavenumber
from pysisyphus.modefollow import geom_davidson
from pysisyphus.modefollow.NormalMode import NormalMode
from pysisyphus.testing import using


@using("xtb")
@pytest.mark.parametrize(
    "precon, ref_cyc, ref_nu",
    [
        (True, 3, 1691.0526969),
        (False, 5, 1691.0311213),
    ],
)
def test_block_davidson_acet(precon, ref_cyc, ref_nu, this_dir):
    geom = geom_loader("lib:acet_tm.xyz")
    calc = XTB(pal=2)
    geom.set_calculator(calc)

    mw_H = geom.mw_hessian
    H = geom.eckart_projection(mw_H)
    w, v = np.linalg.eigh(H)
    inds = [16, 8]

    rg = np.random.default_rng(20180325)

    def get_guess(vec, masses_rep, scale=5e-3):
        return NormalMode(vec + scale * rg.random(*vec.shape), masses_rep)

    guess_modes = [get_guess(v[:, ind], geom.masses_rep) for ind in inds]

    hessian_precon = None
    if precon:
        hessian_precon = np.loadtxt(this_dir / "hessian_precon")

    result = geom_davidson(
        geom, guess_modes, hessian_precon=hessian_precon, start_precon=5, print_level=1,
    )

    nu = result.nus[result.mode_inds[0]]

    assert result.converged
    assert result.cur_cycle == ref_cyc
    assert nu == pytest.approx(ref_nu)
