import numpy as np
import pytest

from pysisyphus.calculators.Gaussian16 import parse_ci_coeffs
from pysisyphus.constants import AU2EV
from pysisyphus.wavefunction import Wavefunction
from pysisyphus.wavefunction.so_coupling import run


@pytest.mark.parametrize(
    "base_name",
    (
        "00_h2co_cart",
        "01_h2co_sph",
    ),
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

    sing_ens = np.array((0.0, 4.0272, 8.4726, 9.2789, 9.5754, 9.7038)) / AU2EV
    trip_ens = np.array((3.3510, 5.6721, 8.0749, 8.1309, 9.3176)) / AU2EV

    ################
    # SO-couplings #
    ################

    socs = run(wf, Xas, Yas, Xat, Yat, sing_ens, trip_ens)

    if not ref_fn.exists():
        np.save(ref_fn, socs)
    ref_socs = np.load(ref_fn)
    np.testing.assert_allclose(socs, ref_socs.reshape(*socs.shape))
