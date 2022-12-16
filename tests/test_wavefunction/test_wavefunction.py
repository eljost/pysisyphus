import json

import h5py
import numpy as np
import pytest

from pysisyphus.calculators.ORCA import parse_orca_cis
from pysisyphus.config import WF_LIB_DIR
from pysisyphus.wavefunction import Wavefunction
from pysisyphus.wavefunction.pop_analysis import (
    mulliken_charges,
    mulliken_charges_from_wf,
    iao_charges_from_wf,
)
from pysisyphus.wavefunction.excited_states import norm_ci_coeffs


@pytest.mark.parametrize(
    "fn, unrestricted",
    (
        ("orca_ch4_sto3g.json", False),
        ("orca_ch4_sto3g_uhf.json", True),
    ),
)
def test_orca_json(fn, unrestricted):
    wf = Wavefunction.from_orca_json(WF_LIB_DIR / fn)
    assert wf.unrestricted == unrestricted
    assert wf.occ == (5, 5)


@pytest.mark.parametrize(
    "fn",
    ("orca_ch4_sto3g.json", "orca_ch4_sto3g_uhf.json"),
)
def test_orca_mulliken(fn):
    wf = Wavefunction.from_orca_json(WF_LIB_DIR / fn)
    charges = mulliken_charges_from_wf(wf)
    np.testing.assert_allclose(
        charges,
        (-0.27748112, 0.06937027, 0.06937027, 0.06937027, 0.06937027),
        atol=2e-5,
    )


@pytest.mark.parametrize(
    "fn",
    (
        "iao_ref.json",
        "iao_ref_uhf.json",
    ),
)
def test_orca_iao(fn):
    wf = Wavefunction.from_orca_json(WF_LIB_DIR / fn)
    charges = iao_charges_from_wf(wf)
    np.testing.assert_allclose(
        charges,
        (-0.523749, 0.130937, 0.130937, 0.130937, 0.130937),
        atol=2e-5,
    )


@pytest.mark.parametrize(
    "fn, ref_dip_mom",
    (
        ("orca_dipmom.json", (0.65100, -0.15774, -0.67608)),
        ("orca_dipmom_uhf.json", (0.00000, 0.00000, -0.34612)),
    ),
)
def test_orca_dipole_moments(fn, ref_dip_mom):
    wf = Wavefunction.from_orca_json(WF_LIB_DIR / fn)
    dip_mom = wf.dipole_moment
    np.testing.assert_allclose(dip_mom, ref_dip_mom, atol=1e-5)


def test_mulliken_exc(this_dir):
    jfn = this_dir / "h2o2" / "rhf_000.000.orca.json"
    ovfn = this_dir / "h2o2" / "overlap_data.h5"

    wf = Wavefunction.from_orca_json(jfn)
    with h5py.File(ovfn) as handle:
        CI = handle["ci_coeffs"][0]
    print()
    ref_charges = (0.399031, 0.025996, 0.474676, 0.369144, 0.009265)
    for i, (tden, ref_charge) in enumerate(zip(CI, ref_charges), 1):
        P_exc = wf.P_exc(tden)
        _ = wf.get_dipole_moment(P=P_exc + P_exc)
        charges = mulliken_charges(
            P=(P_exc / 2, P_exc / 2),
            S=wf.S,
            nuc_charges=wf.nuc_charges,
            ao_centers=wf.ao_centers,
        )
        print(f"\tq: {charges}")
        ref_charges = np.full(4, ref_charge)
        ref_charges[[1, 2]] *= -1  # Oxygens have negative charge
        np.testing.assert_allclose(charges, ref_charges, atol=1e-5)


def test_orca_unrestricted(this_dir):
    """H2O2, HF/3-21G, charge=-1, mult=2"""

    base_fn = this_dir / "h2o2_anion" / "00_h2o2_anion"
    json_fn = base_fn.with_suffix(".json")
    cis_fn = base_fn.with_suffix(".cis")

    wf = Wavefunction.from_orca_json(json_fn)

    P_a, P_b = wf.P
    assert P_a.shape == P_b.shape

    occ_a, occ_b = wf.occ
    assert occ_a == 10
    assert occ_b == 9
    # ORCA uses center of mass by default, which differs a bit from the one
    # that pysisyphus calculates.
    origin = (-0.864214, -1.994637,  1.223355)
    dpm = wf.get_dipole_moment(origin=origin)
    ref_dpm = (0.00047,       0.00016,      -0.00029)
    np.testing.assert_allclose(dpm, ref_dpm, atol=1e-5)

    ref_mulliken = (0.077267, 0.077267, -0.577267, -0.577267)
    mulliken = mulliken_charges_from_wf(wf)
    np.testing.assert_allclose(mulliken, ref_mulliken, atol=1e-5)

    Xa, Ya, Xb, Yb = parse_orca_cis(cis_fn)
    # Xa, Ya = norm_ci_coeffs(Xa, Ya)
    # Xb, Yb = norm_ci_coeffs(Xb, Yb)
    Xa, Ya, Xb, Yb = norm_ci_coeffs(Xa, Ya, Xb, Yb)
    assert Xa.shape == Ya.shape == (3, 10, 12)
    assert Xb.shape == Yb.shape == (3, 9, 13)
    norms_a = np.linalg.norm(Xa, axis=(1, 2))
    norms_b = np.linalg.norm(Xb, axis=(1, 2))
    norms = norms_a**2 + norms_b**2
    ones = np.ones_like(norms)
    np.testing.assert_allclose(norms, ones)

    tdms = wf.get_transition_dipole_moment(Xa, Xb, origin=origin)
    exc_ens = (0.101429, 0.137824, 0.186742)
    fosc = wf.oscillator_strength(exc_ens, tdms)
    ref_fosc = (0.061637792, 0.000000021, 0.000435488)
    np.testing.assert_allclose(fosc, ref_fosc, atol=1e-6)


@pytest.mark.parametrize(
    "charge, ecp_electrons",
    (
        (0, (28,)),
        (
            0,
            {
                0: 28,
            },
        ),
        pytest.param(
            10,
            {
                0: 28,
            },
            marks=pytest.mark.xfail,
        ),
        pytest.param(0, (0,), marks=pytest.mark.xfail),
        pytest.param(0, None, marks=pytest.mark.xfail),
    ),
)
def test_orca_ecp(charge, ecp_electrons):
    fn = WF_LIB_DIR / "orca_yttrium_ecp.molden.input"
    Wavefunction.from_orca_molden(fn, charge=charge, ecp_electrons=ecp_electrons)


def test_quadrupole_moment():
    wf = Wavefunction.from_file(WF_LIB_DIR / "orca_benzene_quad.json")
    quad_moms = wf.quadrupole_moment
    diag = np.diag(quad_moms)
    np.testing.assert_allclose(diag, (-23.80829, -23.80856, -30.13957), atol=1e-5)


def test_one_el_terms():
    fn = WF_LIB_DIR / "orca_ch4_sto3g.json"
    wf = Wavefunction.from_orca_json(fn)
    shells = wf.shells
    T = shells.T_sph
    V = shells.V_sph
    H = T + V

    with open(fn) as handle:
        ref_data = json.load(handle)
    H_ref = ref_data["Molecule"]["H-Matrix"]
    np.testing.assert_allclose(H, H_ref, atol=1e-12)


def test_eval_esp():
    wf = Wavefunction.from_file(WF_LIB_DIR / "orca_lih_ccpvdz.molden.input")
    coords3d = np.array((
        (0.0, 0.0, 1.0),
        (0.0, 0.0, 2.0),
        (1.0, 0.0, 3.0),
    ))
    print()
    esp = wf.eval_esp(coords3d)
    # From orca_vpot
    esp_ref = (0.7086898598813272, 0.0984935998538505, -0.0830229403844891)
    np.testing.assert_allclose(esp, esp_ref, atol=1e-10)
