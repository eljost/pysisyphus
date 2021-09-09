import numpy as np
import pytest

from pysisyphus.calculators import Gaussian16, ORCA5, Turbomole
from pysisyphus.calculators.PySCF import PySCF
from pysisyphus.helpers import geom_loader
from pysisyphus.testing import using


np.set_printoptions(suppress=True, precision=8)


@pytest.mark.parametrize(
    "XY", ["X", "X+Y", "X-Y"]
)
@pytest.mark.parametrize(
    "calc_cls, calc_kwargs",
    [
            # ORCA TDA, Singlet -> Singlet
            pytest.param(
                ORCA5,
                {
                    "keywords": "hf def2-svp nocosx tightscf",
                    "blocks": "%tddft nroots 3 iroot 1 tda true end",
                    "retry_calc": 0,
                },
                marks=using("orca5"),
            ),
            # ORCA TDA, Singlet -> Triplet
            pytest.param(
                ORCA5,
                {
                    "keywords": "hf def2-svp nocosx tightscf",
                    "blocks": "%tddft nroots 3 iroot 1 tda true triplets true end",
                    "retry_calc": 0,
                },
                marks=using("orca5"),
            ),
            pytest.param(
                Turbomole,
                {
                    "control_path": "control_path",
                },
                marks=using("turbomole"),
            ),
            pytest.param(
                PySCF,
                {
                    "basis": "def2svp",
                    "method": "tda",
                    "nstates": 3,
                    "root": 1,
                    "verbose": 4,
                },
                marks=using("pyscf"),
            ),
            pytest.param(
                Gaussian16,
                {
                    "route": "hf def2svp cis=(nstates=3,root=1)"
                },
                marks=using("gaussian16")
            ),
            # Singlet->Triplet
            pytest.param(
                Gaussian16,
                {
                    "route": "hf def2svp cis=(nstates=3,root=1,triplets)"
                },
                marks=using("gaussian16")
            )
    ],
)
def test_h2o2_overlaps(XY, calc_cls, calc_kwargs, this_dir):
    geom = geom_loader("lib:h2o2_def2svp_hf_opt.xyz")

    try:
        calc_kwargs["control_path"] = this_dir / calc_kwargs["control_path"]
    except KeyError:
        pass

    calc_kwargs.update(
        {
            "pal": 4,
            "mem": 1000,
            "track": True,
            "ovlp_type": "tden",
            "XY": XY,
        }
    )
    calc = calc_cls(**calc_kwargs)
    geom.set_calculator(calc)

    energy = geom.energy
    # Self overlap ... should yield identity matrix as states are orthogonal
    ovlps = calc.get_tden_overlaps(indices=(0, 0))
    ref_ovlps = np.eye(ovlps.shape[0])
    np.testing.assert_allclose(ovlps, ref_ovlps, atol=1e-8)
    print("\n", ovlps)
