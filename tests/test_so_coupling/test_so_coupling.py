import numpy as np
import pytest

from pysisyphus.calculators.Gaussian16 import parse_ci_coeffs
from pysisyphus.constants import NU2AU
from pysisyphus.wavefunction import Wavefunction
from pysisyphus.wavefunction.excited_states import norm_ci_coeffs
from pysisyphus.wavefunction.so_coupling import singlet_triplet_so_couplings


# Multiplicative conversion factor from atomic units in Hartree
# to wavenumbers in cm⁻¹.
_AU2NU = 1 / NU2AU


@pytest.mark.parametrize(
    "base_name",
    ("00_h2co_cart", "01_h2co_sph"),
)
def test_socs(base_name, this_dir):
    data_dir = this_dir / "data"
    wf_fn = (data_dir / base_name).with_suffix(".fchk")
    log_fn = wf_fn.with_suffix(".log")
    ref_fn = wf_fn.with_suffix(".npy")

    wf = Wavefunction.from_file(wf_fn)
    Xa, Ya, Xb, Yb = parse_ci_coeffs(log_fn)

    # Singlets and triplets are mixed in one calculation
    singlets = [1, 5, 6, 8, 9]
    triplets = [0, 2, 3, 4, 7]
    Xas = Xa[singlets]
    Yas = Ya[singlets]
    Xat = Xa[triplets]
    Yat = Ya[triplets]

    # Singlet-singlet excitations
    nsings, *_ = Xas.shape
    Xas, Yas = norm_ci_coeffs(Xas, Yas)
    XpYs = Xas + Yas

    # Singlet-triplet excitations
    ntrips, *_ = Xat.shape
    Xat, Yat = norm_ci_coeffs(Xat, Yat)
    XpYt = Xat + Yat

    # The way CI coefficients are normalized in pysisyphus mandates multiplication by
    # sqrt(2). This is also discussed in the PySOC paper.
    # TODO: move this somehwere else?!
    XpYs = np.sqrt(2) * XpYs
    XpYt = np.sqrt(2) * XpYt

    socs = singlet_triplet_so_couplings(wf, XpYs, XpYt)
    socs2 = np.abs(socs) ** 2
    tot_socs = np.sqrt(socs2.sum(axis=1))

    socs_nu = socs * _AU2NU
    tot_socs_nu = tot_socs * _AU2NU

    cur_st = 0
    fmt = "> 10.5f"
    print()
    for sstate in range(nsings + 1):
        for tstate in range(ntrips):
            ts = tot_socs_nu[cur_st]
            sm1, s0, sp1 = np.abs(socs_nu[cur_st])
            print(
                f"<S{sstate}|Hso|T{tstate+1}>: {ts:{fmt}}  {sm1:{fmt}} {s0:{fmt}} {sp1:{fmt}}"
            )
            cur_st += 1
        print()

    if not ref_fn.exists():
        np.save(ref_fn, socs)
    ref_socs = np.load(ref_fn)
    np.testing.assert_allclose(socs, ref_socs)
