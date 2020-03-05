from pathlib import Path

import numpy as np
import pytest

from pysisyphus.helpers import geom_from_library
from pysisyphus.init_logging import init_logging
from pysisyphus.calculators import ORCA, Gaussian16, Turbomole
from pysisyphus.calculators.PySCF import PySCF
from pysisyphus.optimizers.QuickMin import QuickMin
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
    geom = geom_from_library("benzene.xyz")

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


def test_geometry_get_restart_info():
    geom = geom_from_library("benzene.xyz")
    calc = PySCF(method="scf", basis="def2svp")

    geom.set_calculator(calc)
    restart = geom.get_restart_info()

    atoms = restart["atoms"]
    coords = restart["cart_coords"]

    assert atoms == geom.atoms
    assert len(coords) == len(geom.atoms * 3)
    assert "calc_info" in restart


def test_opt_restart():
    def get_calc():
        return ORCA("HF sto-3g")

    def get_geom():
        geom = geom_from_library("h2o_shaken.xyz")
        geom.set_calculator(get_calc())
        return geom

    def get_opt(geom, restart_info=None):
        opt_kwargs = {
            "max_cycles": 2,
            "restart_info": restart_info,
            "dump": True,
        }
        opt = QuickMin(geom, **opt_kwargs)
        return opt

    # Initial run
    geom = get_geom()
    opt = get_opt(geom)
    opt.run()
    restart_info = opt.get_restart_info()

    # Restarted run
    re_geom = get_geom()
    re_opt = get_opt(re_geom, restart_info)
    re_opt.run()

    # import pdb; pdb.set_trace()
    pass
