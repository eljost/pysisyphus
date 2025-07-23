import numpy as np
import pytest

from pysisyphus.calculators.ORCA import parse_orca_cis
from pysisyphus.wavefunction import Wavefunction
from pysisyphus.wavefunction.excited_states import (
    norm_ci_coeffs,
    make_density_matrices_for_root,
    get_state_to_state_transition_density,
)

np.set_printoptions(suppress=True, precision=6, linewidth=180)


@pytest.mark.parametrize(
    "base, dpm_ref",
    (
        ("00_h2o_restricted", (0.0, -0.49594, 0.0)),
        ("01_h2o_unrestricted", (0.0, -0.43172, 0.0)),
    ),
)
def test_exc_dens(base, dpm_ref, this_dir):
    base = this_dir / "data" / base
    wf_fn = base.with_suffix(".bson")
    wf = Wavefunction.from_file(wf_fn)
    print(wf)
    cis_fn = base.with_suffix(".cis")
    restricted = wf.restricted
    Xa, Ya, Xb, Yb = parse_orca_cis(cis_fn, restricted_same_ab=True)
    Xa, Ya, Xb, Yb = norm_ci_coeffs(Xa, Ya, Xb, Yb)

    root = 1
    rootm1 = root - 1
    Ca, Cb = wf.C
    Pa_exc, Pb_exc = make_density_matrices_for_root(
        rootm1, restricted, Xa, Ya, Xb, Yb, Ca=Ca, Cb=Cb
    )
    P_exc = Pa_exc + Pb_exc
    dpm = wf.get_dipole_moment(P=P_exc)
    np.testing.assert_allclose(dpm, dpm_ref, atol=1e-5)


@pytest.mark.parametrize(
    "base, ref_fn",
    (
        ("02_h2o_dotrans_restricted", "02_trans_moms_ref"),
        ("03_h2o_dotrans_unrestricted", "03_trans_moms_ref"),
    ),
)
def test_trans_dens(base, ref_fn, this_dir):
    tdpm_ref = np.loadtxt(this_dir / "data" / ref_fn)[:, 6:]

    base = this_dir / "data" / base
    wf_fn = base.with_suffix(".bson")
    wf = Wavefunction.from_file(wf_fn)
    print(wf)
    cis_fn = base.with_suffix(".cis")
    Xa, Ya, Xb, Yb = parse_orca_cis(cis_fn, restricted_same_ab=True)
    Xa, Ya, Xb, Yb = norm_ci_coeffs(Xa, Ya, Xb, Yb)

    XpYa = Xa + Ya
    XpYb = Xb + Yb

    tdpm = np.zeros_like(tdpm_ref)
    nroots = Xa.shape[0]
    i = 0
    for init in range(nroots - 1):
        for final in range(init + 1, nroots):
            Pa = get_state_to_state_transition_density(XpYa[init], XpYa[final])
            Pb = get_state_to_state_transition_density(XpYb[init], XpYb[final])
            tdpm[i] = wf.get_transition_dipole_moment(Pa, Pb, full=True)
            i += 1
    np.testing.assert_allclose(tdpm, tdpm_ref, atol=1e-5)
