import numpy as np
import pytest

from pysisyphus.wavefunction import Wavefunction

from pysisyphus.calculators.Turbomole import get_density_matrices_for_root


@pytest.mark.parametrize(
    "path, vec_fn, root, ref_dpm",
    (
        ("rtda", "ciss_a", 2, (0.0, 0.0, 0.449441)),
        ("utda", "ucis_a", 2, (0.0, 0.0, 0.271004)),
        ("rtd", "sing_a", 2, (0.0, 0.0, 0.45343)),
        ("utd", "unrs_a", 2, (-0.101564, 0.911588, -0.000159)),
    ),
)
def test(path, vec_fn, root, ref_dpm, this_dir):
    path = this_dir / "make_dens" / path

    wf_fn = path / "wavefunction.molden"
    wf = Wavefunction.from_file(wf_fn)
    log_fn = path / "turbomole.out"
    vec_fn = path / vec_fn
    rlx_vec_fn = path / "dipl_a"
    Ca, Cb = wf.C
    Pa, Pb = get_density_matrices_for_root(log_fn, vec_fn, root, rlx_vec_fn, Ca, Cb)
    P_tot = Pa + Pb
    dpm = wf.get_dipole_moment(P_tot)
    np.testing.assert_allclose(dpm, ref_dpm, atol=3e-5)
