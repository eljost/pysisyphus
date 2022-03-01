from pathlib import Path

import pytest

from pysisyphus.calculators.PySCF import PySCF
from pysisyphus.calculators import Turbomole, ORCA, ORCA5, Gaussian16, Psi4
from pysisyphus.helpers import geom_loader
from pysisyphus.init_logging import init_logging
from pysisyphus.optimizers.RFOptimizer import RFOptimizer
from pysisyphus.testing import using


init_logging()


@pytest.mark.parametrize(
    "calc_cls, calc_kwargs_",
    [
        pytest.param(Gaussian16,
            {"route": "HF/STO-3G"},
            marks=using("gaussian16")),
        pytest.param(ORCA,
            {"keywords": "HF STO-3G tightscf"},
            marks=using("orca")),
        # pytest.param(Psi4,
            # {"method": "scf", "basis": "sto-3g", "to_set": {"scf_type": "pk"}},
            # marks=using("psi4"),
        # ),
        pytest.param(PySCF,
            {"basis": "sto-3g"},
            marks=using("pyscf")),
        pytest.param(Turbomole,
            {"control_path": "./control_path_hf_sto3g_gs"},
            marks=using("turbomole")),
])
def test_cytosin_gs_opt(calc_cls, calc_kwargs_, this_dir):
    geom = geom_loader("lib:cytosin.xyz", coord_type="redund")

    print("@Using", calc_cls)
    calc_kwargs = {
        "pal": 2,
        "mem": 1000,
    }
    calc_kwargs.update(calc_kwargs_)

    if "control_path" in calc_kwargs:
        calc_kwargs["control_path"] = this_dir / calc_kwargs["control_path"]

    calc = calc_cls(**calc_kwargs)
    geom.set_calculator(calc)

    opt_kwargs = {
        "thresh": "gau_tight",
        "overachieve_factor": 2.,
        "line_search": True,
        "gdiis": True,
    }
    opt = RFOptimizer(geom, **opt_kwargs)
    opt.run()

    assert opt.is_converged
    assert geom.energy == pytest.approx(-387.54925356)


@pytest.mark.parametrize(
    "calc_cls, calc_kwargs",
    [
        pytest.param(Gaussian16,
            {"route": "PBE1PBE/def2SVP TD=(nstates=2,root=1)"},
            marks=using("gaussian16")
        ),
        pytest.param(Turbomole,
            {"control_path": "./control_path_pbe0_def2svp_s1"},
            marks=using("turbomole")
        ),
        pytest.param(ORCA5,
            {"keywords": "PBE0 def2-SVP tightscf",
             "blocks": "%tddft nroots 2 iroot 1 tda false end"},
            marks=using("orca5")
        ),
        pytest.param(PySCF,
            {"xc": "pbe0", "method": "tddft", "basis": "def2SVP",
             "nstates": 2, "root": 1},
            # Skip for now, as this takes 30 min in the CI
            marks=(using("pyscf"), pytest.mark.skip),
        )
])
def test_cytosin_s1_opt(calc_cls, calc_kwargs, this_dir):
    geom = geom_loader("lib:cytosin.xyz", coord_type="redund")

    if "control_path" in calc_kwargs:
        calc_kwargs["control_path"] = this_dir / calc_kwargs["control_path"]

    calc_kwargs.update({
        "mem": 2000,
        "pal": 4,
        "ovlp_type": "tden",
        "mos_renorm": True,
        "track": True,
    })
    calc = calc_cls(**calc_kwargs)
    geom.set_calculator(calc)

    opt_kwargs = {
        "thresh": "gau",
        "overachieve_factor": 2.,
        # "trust_radius": 0.3,
        # "trust_max": 0.3,
        "line_search": True,
        "gdiis": True,
    }
    opt = RFOptimizer(geom, **opt_kwargs)
    opt.run()

    assert opt.is_converged
    # assert geom.energy == pytest.approx(-394.06081796)
    assert calc.root_flips[2]
    assert calc.root == 2
