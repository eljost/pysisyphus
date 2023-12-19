import numpy as np
import pytest

from pysisyphus.calculators import ORCA, Gaussian16, Turbomole, DFTBp
from pysisyphus.calculators.PySCF import PySCF
from pysisyphus.helpers import geom_loader
from pysisyphus.testing import using


@pytest.mark.parametrize(
    "mult, ref_energies",
    (
        (
            1,
            # Reference energies from Turbomole
            (
                -75.96023978184000,
                -75.62222276222431,
                -75.55696018049534,
                -75.53779576649994,
            ),
        ),
        (
            3,
            # Reference energies from Turbomole
            (
                -75.71293821803000,
                -75.64358728543596,
                -75.62823894103637,
                -75.46436197245151,
            ),
        ),
    ),
)
@pytest.mark.parametrize(
    "calc_cls, calc_kwargs",
    (
        pytest.param(
            Turbomole,
            {},
            marks=using("turbomole"),
        ),
        pytest.param(
            ORCA,
            {
                "keywords": "hf def2-svp tightscf norijcosx",
                "blocks": "%tddft tda false nroots 3 end",
            },
            marks=using("orca"),
        ),
        pytest.param(
            PySCF,
            {
                "method": "tdhf",
                "basis": "def2svp",
                "nroots": 3,
            },
            marks=using("PySCF"),
        ),
        pytest.param(
            Gaussian16,
            {
                "route": "hf def2svp td(nstates=3)",
            },
            marks=using("gaussian16"),
        ),
    ),
)
def test_h2o_all_energies(mult, ref_energies, calc_cls, calc_kwargs, this_dir):
    geom = geom_loader("lib:h2o.xyz")
    # Fix control path for Turbomole
    if calc_cls == Turbomole:
        calc_kwargs["control_path"] = this_dir / f"control_path_mult_{mult}"
    calc_kwargs["mult"] = mult
    calc_kwargs["base_name"] = f"calc_mult_{mult}"
    calc = calc_cls(**calc_kwargs)
    geom.set_calculator(calc)
    all_energies = geom.all_energies
    # PySCF and Turbomole agree extermely well, at least in the restricted calculation.
    # ORCA deviates up to 5e-5 Eh.
    np.testing.assert_allclose(all_energies, ref_energies, atol=5e-5)

    # As we did not set any root the geometries energy should correspond to the GS energy
    energy = geom.energy
    assert energy == pytest.approx(ref_energies[0])

    for root in range(4):
        root_en = geom.get_root_energy(root)
        assert root_en == pytest.approx(ref_energies[root])
