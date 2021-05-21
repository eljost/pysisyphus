from pathlib import Path

import numpy as np
import pytest

from pysisyphus.helpers import geom_loader
from pysisyphus.init_logging import init_logging
from pysisyphus.calculators import ORCA, Gaussian16, Turbomole
from pysisyphus.calculators.PySCF import PySCF
from pysisyphus.optimizers.ConjugateGradient import ConjugateGradient
from pysisyphus.optimizers.FIRE import FIRE
from pysisyphus.optimizers.LBFGS import LBFGS
from pysisyphus.optimizers.QuickMin import QuickMin
from pysisyphus.optimizers.RFOptimizer import RFOptimizer
from pysisyphus.optimizers.SteepestDescent import SteepestDescent
from pysisyphus.optimizers.PreconLBFGS import PreconLBFGS
from pysisyphus.testing import using


init_logging()


@pytest.fixture
def this_dir(request):
    return Path(request.module.__file__).parents[0]


@pytest.mark.parametrize(
    "calc_cls, calc_kwargs, chk_exts",
    [
        pytest.param(
            ORCA, {"keywords": "HF def2-SVP", }, ("gbw", ),
            marks=using("orca")),
        pytest.param(
            Gaussian16, {"route": "HF/def2SVP", }, ("fchk", ),
            marks=using("gaussian16")),
        pytest.param(
            PySCF, {"method": "scf", "basis": "def2svp"}, ("chkfile", ),
            marks=using("pyscf")),
        pytest.param(
            Turbomole, {"control_path": "benzene_control_path"}, ("mos", ),
            marks=using("turbomole")),
        pytest.param(
            Turbomole, {"control_path": "benzene_control_path_uhf"}, ("alpha", "beta"),
            marks=using("turbomole")),
    ]
)
def test_restart(calc_cls, calc_kwargs, chk_exts, this_dir):
    geom = geom_loader("lib:benzene.xyz")

    if calc_cls == Turbomole:
        calc_kwargs["control_path"] = this_dir / calc_kwargs["control_path"]
    calc = calc_cls(**calc_kwargs)
    geom.set_calculator(calc)

    forces = geom.forces

    ref_energy = -230.534572588
    assert geom.energy == pytest.approx(ref_energy)

    restart_info = calc.get_restart_info()
    chkfiles = restart_info["chkfiles"]
    assert len(chkfiles) >= 1
    assert all([ext in chkfiles for ext in chk_exts])

    calc2 = calc_cls(**calc_kwargs)
    calc2.set_restart_info(restart_info)
    geom2 = geom.copy()
    geom2.set_calculator(calc2)
    forces2 = geom2.forces

    np.testing.assert_allclose(forces2, forces, atol=1e-5)


@using("pyscf")
def test_geometry_get_restart_info():
    geom = geom_loader("lib:benzene.xyz")
    calc = PySCF(method="scf", basis="def2svp")

    geom.set_calculator(calc)
    restart = geom.get_restart_info()

    atoms = restart["atoms"]
    coords = restart["cart_coords"]

    assert atoms == geom.atoms
    assert len(coords) == len(geom.atoms * 3)
    assert "calc_info" in restart


@using("pyscf")
@pytest.mark.parametrize(
    "opt_cls, opt_kwargs_, ref_norm",
    [
        pytest.param(ConjugateGradient, {}, 0.01693523, marks=using("pyscf")),
        pytest.param(FIRE, {}, 0.50285483, marks=using("pyscf")),
        pytest.param(LBFGS, {"gamma_mult": True, }, 2.2002337e-6, marks=using("pyscf")),
        pytest.param(LBFGS, {"gamma_mult": False, }, 1.36271012e-5, marks=using("pyscf")),
        pytest.param(PreconLBFGS, {}, 9.11439241e-6, marks=using("pyscf")),
        pytest.param(QuickMin, {}, 0.02305389, marks=using("pyscf")),
        pytest.param(RFOptimizer, {}, 0.00189796616, marks=using("pyscf")),
        pytest.param(SteepestDescent, {}, 0.05535400, marks=using("pyscf")),
    ]
)
def test_opt_restart(opt_cls, opt_kwargs_, ref_norm):
    def get_calc():
        return PySCF(method="scf", basis="def2svp")

    def get_geom():
        geom = geom_loader("lib:h2o_shaken.xyz")
        geom.set_calculator(get_calc())
        return geom

    def get_opt(geom, max_cycles, restart_info=None):
        opt_kwargs = {
            "max_cycles": max_cycles,
            "restart_info": restart_info,
            "dump": True,
            "thresh": "gau_tight",
        }
        opt_kwargs.update(opt_kwargs_)
        opt = opt_cls(geom, **opt_kwargs)
        return opt

    max_cycles = 4

    # Reference run
    ref_geom = get_geom()
    ref_opt = get_opt(ref_geom, max_cycles=2*max_cycles)
    ref_opt.run()

    ref_energy = ref_geom.energy
    ref_forces = ref_geom.forces

    assert ref_opt.cur_cycle == 2*max_cycles - 1
    assert np.linalg.norm(ref_forces) == pytest.approx(ref_norm)

    # Try two runs, each with 1*max_cycles and restart inbetween
    # First run
    first_geom = get_geom()
    first_opt = get_opt(first_geom, max_cycles)
    first_opt.run()
    # Get restart info
    restart_info = first_opt.get_restart_info()

    # Second run, restarted from the first one
    re_geom = get_geom()
    re_opt = get_opt(re_geom, max_cycles, restart_info)
    re_opt.run()

    assert re_opt.cur_cycle == 2*max_cycles - 1
    assert re_geom.energy == pytest.approx(ref_energy)
    re_forces = re_geom.forces
    np.testing.assert_allclose(re_forces, ref_forces, atol=1e-6)
