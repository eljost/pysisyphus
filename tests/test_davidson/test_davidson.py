import os
from pathlib import Path

import numpy as np
import pytest

from pysisyphus.calculators import ORCA, XTB
from pysisyphus.calculators.PySCF import PySCF
from pysisyphus.helpers import geom_loader
from pysisyphus.helpers_pure import eigval_to_wavenumber
from pysisyphus.init_logging import init_logging
from pysisyphus.modefollow import davidson
from pysisyphus.modefollow.davidson import block_davidson
from pysisyphus.modefollow.NormalMode import NormalMode
from pysisyphus.testing import using


THIS_DIR = Path(os.path.abspath(os.path.dirname(__file__)))


@pytest.mark.parametrize(
    "precon", [
        True,
        False,
    ]
)
@pytest.mark.parametrize(
    "calc_cls, calc_kwargs, ref_nu_pre, ref_cyc_pre, ref_nu, ref_cyc", [
        pytest.param(
            # 1791.60 cm⁻¹ from ORCA itself
            ORCA, {"keywords": "BP86 def2-SVP"},
            1793.274093, 3, 1793.293512, 7, marks=using("orca")),
        pytest.param(
            # 1691.888 cm⁻¹ from XTB itself
            XTB, {},
            1691.052689, 3, 1691.005160, 7, marks=using("xtb")),
        # pytest.param(
            # PySCF, {"basis": "321g",},
            # 1889.2027727, 16, 1888.9849030442317, 16, marks=using("pyscf")),
        pytest.param(
            PySCF, {"basis": "def2svp", "xc": "bp86"},
            1790.744, 3, 1790.342, 7, marks=using("pyscf")),
    ],
)
def test_davidson_acet(precon, calc_cls, calc_kwargs,
                       # ref_nu, ref_cyc, ref_nu_pre, ref_cyc_pre):
                       ref_nu_pre, ref_cyc_pre, ref_nu, ref_cyc):
    geom = geom_loader("lib:acet_tm.xyz")

    calc_kwargs.update({
        "pal": 2,
    })
    calc = calc_cls(**calc_kwargs)
    geom.set_calculator(calc)

    # Initial guess
    l = np.zeros_like(geom.coords).reshape(-1, 3)
    l[0][2] = 0.8
    l[1][2] = -0.6
    q = NormalMode(l.flatten(), geom.masses_rep)
    
    if precon:
        # Hessian for preconditioner
        # precon_calc = XTB(pal=2)
        # precon_geom = geom.copy()
        # precon_geom.set_calculator(precon_calc)
        # hessian_precon = precon_geom.eckart_projection(precon_geom.mw_hessian)
        # np.savetxt("hessian_precon", hessian_precon)
        hessian_precon = np.loadtxt(THIS_DIR / "hessian_precon")
    else:
        hessian_precon = None

    davidson_kwargs = {
        "hessian_precon": hessian_precon,
    }
    result = davidson(geom, q, **davidson_kwargs)

    nu = result.nus[result.mode_ind]
    print("nu", nu)
    print("cur_cycle", result.cur_cycle)

    ref_cycles = ref_cyc_pre if precon else ref_cyc
    nu_ref = ref_nu_pre if precon else ref_nu

    assert nu == pytest.approx(nu_ref, abs=0.5)
    assert result.cur_cycle+1 == ref_cycles



# def test_block_davidson():
    # geom = geom_loader("block_davidson/xtbopt.xyz")
    # # geom.jmol()

    # H = np.loadtxt("block_davidson/xtbhess")
    # mw_H = geom.eckart_projection(geom.mass_weigh_hessian(H))
    # w, v = np.linalg.eigh(mw_H)
    # # nus = eigval_to_wavenumber(w)
    # inds = [20, 28]
    # b = get_guess(v, inds)
    # import pdb; pdb.set_trace()
    # # w, v = np.linalg.eigh(H)
    # # ind = [6, 7]
    # # init_guess = v[:,ind]
    # # init_guess += np.random.rand(*init_guess.shape)
    # # init_guess /= np.linalg.norm(init_guess, axis=0)


def get_guess(vec, masses_rep, scale=0.0):
    return NormalMode(-vec + scale * np.random.rand(*vec.shape), masses_rep)


@pytest.mark.parametrize(
    "precon", [True, False]
)
def test_block_davidson_acet(precon):
    geom = geom_loader("lib:acet_tm.xyz")
    calc = XTB(pal=2)
    geom.set_calculator(calc)

    mw_H = geom.mw_hessian
    H = geom.eckart_projection(mw_H)
    w, v = np.linalg.eigh(H)
    nus = eigval_to_wavenumber(w)
    # inds = [16, 8]
    inds = [16, ]
    guess_modes = [get_guess(v[:,ind], geom.masses_rep) for ind in inds]

    # l = np.zeros_like(geom.coords).reshape(-1, 3)
    # l[0][2] = 0.8
    # l[1][2] = -0.6
    # q = NormalMode(l.flatten(), geom.masses_rep)
    # guess_modes = (q, )

    if precon:
        hessian_precon = np.loadtxt(THIS_DIR / "hessian_precon")
        nu_ref = 1691.052689
    else:
        hessian_precon = None
        nu_ref = 1691.005160

    result = block_davidson(geom, guess_modes, hessian_precon=hessian_precon, max_cycles=5)

    nu = result.nus[result.mode_inds[0]]
    print("nu", nu)
    print("cur_cycle", result.cur_cycle)

    assert nu == pytest.approx(nu_ref, abs=0.5)
