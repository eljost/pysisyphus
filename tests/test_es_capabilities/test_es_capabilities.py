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
    wf = geom.calc_wavefunction()["wavefunction"]
    pop_ana = mulliken_charges(wf)
    np.testing.assert_allclose(pop_ana.charges, ref_charges, atol=4e-6)
    # Access geom._energy to ensure that energy is actually set
    assert geom._energy == pytest.approx(ref_energy)


@pytest.mark.parametrize(
    "root, ref_energy, ref_dpm",
    (
        # Ref dipole moment from ORCA calculation
        (0, -98.572847, (0.0, 0.0, -0.49257)),
        (2, -98.107181, (0.0, 0.0, 0.44944)),
    ),
)
@pytest.mark.parametrize(
    "calc_cls, calc_kwargs",
    (
        pytest.param(
            ORCA,
            {
                "keywords": "hf sto-3g tightscf",
                "blocks": "%tddft tda true nroots 2 end",
                "root": 1,
            },
            marks=using("orca"),
        ),
        pytest.param(
            Turbomole,
            {},
            marks=using("turbomole"),
        ),
        pytest.param(
            Gaussian16,
            {
                "route": "hf sto-3g cis(nstates=2)",
                "root": 1,
            },
            marks=using("gaussian16"),
        ),
    ),
)
def test_relaxed_density(root, ref_energy, ref_dpm, calc_cls, calc_kwargs, this_dir):
    """
    PySCF buries the calculation of relaxed densities in the
    respective gred_elec() functions in the pyscf.grad.{tdrhf,tduhf} modules.
    I'm not in the mood to extract the respective parts into separate functions.

    Even though Turbomole keeps relaxed densities to itself, it still dumps
    the OV-correction to the density matrix to a file, so we can reconstruct
    the relaxed density matrix.

    DFTB+ does not seem to allow wavefunction export, so there is no point in
    trying to get the relaxed density.
    """
    geom = geom_loader("lib:hf_hf_sto3g_opt.xyz")
    calc_kwargs.update(
        {
            "base_name": f"root_{root}",
        }
    )
    if calc_cls == Turbomole:
        calc_kwargs["control_path"] = this_dir / f"control_path_hf"
    calc = calc_cls(**calc_kwargs)
    geom.set_calculator(calc)

    result = geom.calc_relaxed_density(root)
    dens = result["density"]
    wf = geom.wavefunction
    dens_key = wf.set_relaxed_density(root, dens)
    P_tot = dens.sum(axis=0)
    dpm = wf.get_dipole_moment(P_tot)

    # Also test context manager and shortcut
    with wf.current_density(dens_key):
        # wf.get_dipole_moment() and wf.dipole_moment will use the currently selected density.
        context_dpm = wf.get_dipole_moment()
        context_dpm2 = wf.dipole_moment

    energy = result["energy"]
    assert energy == pytest.approx(ref_energy)

    atol = 1e-5
    np.testing.assert_allclose(dpm, ref_dpm, atol=atol)
    np.testing.assert_allclose(context_dpm, ref_dpm, atol=atol)
    np.testing.assert_allclose(context_dpm2, ref_dpm, atol=atol)


@pytest.mark.parametrize(
    "calc_cls, calc_kwargs",
    (
        pytest.param(Turbomole, {}, marks=using("turbomole")),
        pytest.param(
            ORCA,
            {
                "keywords": "hf def2-svp",
                "blocks": "%tddft tda false nroots 3 end",
            },
            marks=using("orca"),
        ),
        pytest.param(
            Gaussian16,
            {
                "route": "hf def2svp td(nstates=3)",
                "iop9_40": 4,
            },
            marks=using("gaussian16"),
        ),
    ),
)
def test_td_1tdms(calc_cls, calc_kwargs, this_dir):
    mult = 1
    geom = geom_loader("lib:h2o.xyz")
    # Fix control path for Turbomole
    if calc_cls == Turbomole:
        calc_kwargs["control_path"] = this_dir / f"control_path_mult_{mult}"
    calc_kwargs["mult"] = mult
    calc_kwargs["base_name"] = f"calc_mult_{mult}"
    calc = calc_cls(**calc_kwargs)
    geom.set_calculator(calc)

    Xa, Ya, Xb, Yb = geom.td_1tdms
    XpYa = Xa + Ya
    XpYb = Xb + Yb

    wf = calc.get_stored_wavefunction()

    tdms = wf.get_transition_dipole_moment(XpYa, XpYb)
    # From Turbomole
    ref_tdms = np.array(
        (
            (0.0, 0.0, 0.319042),
            (0.0, 0.0, 0.0),
            (0.0, -0.567422, 0.0),
        )
    )
    # Compare squares of TDMs, as the sign of each state's TDM is undetermined
    np.testing.assert_allclose(tdms**2, ref_tdms**2, atol=6e-5)
