import numpy as np
import pytest

from pysisyphus.calculators import ORCA, Gaussian16, Turbomole, DFTBp
from pysisyphus.calculators.PySCF import PySCF
from pysisyphus.wavefunction.pop_analysis import mulliken_charges
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


@using("dftbp")
def test_dftbp_h2o_all_energies():
    geom = geom_loader("lib:h2o.xyz")
    nroots = 4
    calc_kwargs = {
        "mult": 1,
        "parameter": "mio-ext",
        "nroots": nroots,
    }
    calc = DFTBp(**calc_kwargs)
    geom.set_calculator(calc)
    all_energies = geom.all_energies
    ref_energies = (
        -4.07775075,
        -3.39575683,
        -3.33313599,
        -3.24453337,
        -3.21035650,
    )
    np.testing.assert_allclose(all_energies, ref_energies)

    # As we did not set any root the geometries energy should correspond to the GS energy
    energy = geom.energy
    assert energy == pytest.approx(ref_energies[0])

    for root in range(nroots + 1):
        root_en = geom.get_root_energy(root)
        assert root_en == pytest.approx(ref_energies[root])


@pytest.mark.parametrize(
    "mult, ref_energy, ref_charges",
    (
        # Ref charges from ORCA calculation
        (1, -74.960702, (-0.372544, 0.186272, 0.186272)),
        (3, -74.590125, (0.043753, -0.021876, -0.021876)),
    ),
)
@pytest.mark.parametrize(
    "calc_cls, calc_kwargs",
    (
        pytest.param(
            ORCA,
            {
                "keywords": "hf sto-3g",
            },
            marks=using("orca"),
        ),
        pytest.param(
            Turbomole,
            {},
            marks=using("turbomole"),
        ),
        pytest.param(
            PySCF,
            {
                "method": "scf",
                "basis": "sto3g",
            },
            marks=using("PySCF"),
        ),
        #
        # DFTB+ seems to use Slater-type-orbitals and does not
        # seeom to support a wavefunction export in a known/supported
        # format. So there is no DFTB+ test here.
        #
        pytest.param(
            Gaussian16,
            {
                "route": "hf sto-3g",
            },
            marks=using("gaussian16"),
        ),
    ),
)
def test_get_wavefunction(
    mult, ref_energy, ref_charges, calc_cls, calc_kwargs, this_dir
):
    geom = geom_loader("lib:h2o.xyz")
    if calc_cls == Turbomole:
        calc_kwargs["control_path"] = this_dir / f"control_path_h2o_mult_{mult}"
    calc_kwargs["mult"] = mult
    calc_kwargs["base_name"] = f"calc_mult_{mult}"
    calc = calc_cls(**calc_kwargs)
    geom.set_calculator(calc)
    wf = geom.get_wavefunction()["wavefunction"]
    pop_ana = mulliken_charges(wf)
    np.testing.assert_allclose(pop_ana.charges, ref_charges, atol=4e-6)
    assert geom.energy == pytest.approx(ref_energy)
