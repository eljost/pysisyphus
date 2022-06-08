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
    dip_mom = wf.dipole_moment()
    np.testing.assert_allclose(dip_mom, ref_dip_mom, atol=1e-5)


def test_transition_dipole_moments(this_dir):
    def get_fn(fn):
        return this_dir / "ch4_3states" / fn

    wf = Wavefunction.from_orca_json(get_fn("tda.json"))

    step = 0
    fn = get_fn("tda_overlap_data.h5")
    with h5py.File(fn) as handle:
        all_ens = handle["all_energies"][step]
        CI = handle["ci_coeffs"][step]
    exc_ens = all_ens[1:] - all_ens[0]

    tdms = wf.transition_dipole_moment(CI)
    fosc = wf.oscillator_strength(exc_ens, tdms)
    fmt = " .6f"
    for i, (f, tdm) in enumerate(zip(fosc, tdms)):
        x, y, z = tdm
        tot = (tdm ** 2).sum()
        print(f"{i}: {f=:{fmt}} {x:{fmt}} {y:{fmt}} {z:{fmt}}, tot={tot:{fmt}}")

    np.testing.assert_allclose(fosc, (0.506729325, 0.506497773, 0.506627503), atol=1e-6)


def test_P_exc(this_dir):
    def get_fn(fn):
        return this_dir / "ch4_3states" / fn

    wf = Wavefunction.from_orca_json(get_fn("tda.json"))
    fn = get_fn("tda_overlap_data.h5")
    with h5py.File(fn) as handle:
        CI = handle["ci_coeffs"][0]
    S1 = CI[0]
    P_exc = wf.P_exc(S1)
    exc_dpm = wf.dipole_moment(P_exc[None, :, :])
    print(exc_dpm)


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
        _ = wf.dipole_moment((P_exc, P_exc))
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

    fn = this_dir / "h2o2_anion" / "05_h2o2_anion.json"
    wf = Wavefunction.from_orca_json(fn)

    P_a, P_b = wf.P
    assert P_a.shape == P_b.shape

    occ_a, occ_b = wf.occ
    assert occ_a == 10
    assert occ_b == 9
    # ORCA uses center of mass by default, which differs a bit from the one
    # that pysisyphus calculates.
    dpm = wf.dipole_moment(origin=(-1.487808, 3.034774, 0.292895))
    ref_dpm = (-0.00017, -0.00085, 0.00235)
    np.testing.assert_allclose(dpm, ref_dpm, atol=1e-5)

    ref_mulliken = (0.165653, -0.665653, -0.665653, 0.165653)
    mulliken = mulliken_charges_from_wf(wf)
    np.testing.assert_allclose(mulliken, ref_mulliken, atol=1e-5)

    # cis = "/home/johannes/tmp/284_marcus/05_h2o2_anion/05_h2o2_anion.cis"
    cis = this_dir / "h2o2_anion" / "05_h2o2_anion.cis"
    Xa, Ya, Xb, Yb = parse_orca_cis(cis)
    assert Xa.shape == Ya.shape == (5, 10, 12)
    assert Xb.shape == Yb.shape == (5, 9, 13)
    norms_a = np.linalg.norm(Xa, axis=(1, 2))
    norms_b = np.linalg.norm(Xb, axis=(1, 2))
    norms = norms_a ** 2 + norms_b ** 2
    ones = np.ones_like(norms)
    np.testing.assert_allclose(norms, ones)

    tdms = wf.transition_dipole_moment(Xa, Xb)
    exc_ens = (0.080127, 0.128713, 0.214026, 0.230042, 0.246588)
    fosc = wf.oscillator_strength(exc_ens, tdms)
    ref_fosc = (0.030505412, 0.000000181, 0.000000011, 0.000072701, 0.000000226)
    np.testing.assert_allclose(fosc, ref_fosc, atol=1e-6)
