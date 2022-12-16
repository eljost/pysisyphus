import numpy as np

from pysisyphus.calculators.ORCA import parse_orca_cis
from pysisyphus.wavefunction import Wavefunction
from pysisyphus.wavefunction.excited_states import top_differences, norm_ci_coeffs


def test_top_differences(this_dir):
    base = this_dir / "benzene"

    # Transition densities
    cis_ref = base / "calculator_000.000.orca.cis"
    Xra, Yra, *_ = parse_orca_cis(cis_ref)
    Xra, Yra = norm_ci_coeffs(Xra, Yra)

    cis_cur = base / "calculator_000.001.orca.cis"
    Xca, Yca, *_ = parse_orca_cis(cis_cur)
    Xca, Yca = norm_ci_coeffs(Xca, Yca)

    # MO-Overlaps
    jsr = base / "calculator_000.000.orca.json"
    wfr = Wavefunction.from_orca_json(jsr)
    jsc = base / "calculator_000.001.orca.json"
    wfc = Wavefunction.from_orca_json(jsc)
    S_MO = wfr.S_MO_with(wfc)

    # Closed shell case
    diffs = 2 * top_differences(Xra, Yra, Xca, Yca, S_MO)

    # Open shell case
    # Renormalize to get the appropriate norms for alpha and beta.
    Xra, Yra, Xrb, Yrb = norm_ci_coeffs(Xra, Yra, Xra, Yra)
    Xca, Yca, Xcb, Ycb = norm_ci_coeffs(Xca, Yca, Xca, Yca)
    alpha_diffs = top_differences(Xra, Yra, Xca, Yca, S_MO)
    beta_diffs = top_differences(Xrb, Yrb, Xcb, Ycb, S_MO)
    np.testing.assert_allclose(alpha_diffs, beta_diffs)

    open_shell_diffs = alpha_diffs + beta_diffs
    np.testing.assert_allclose(diffs, open_shell_diffs)

    ref_diffs = np.array(
        (
            (0.00425723, 0.80192988, 0.91658875, 0.90855886),
            (0.79511846, 0.00439037, 0.88472107, 0.89259826),
            (0.80270685, 0.78982435, 0.64191080, 0.03065335),
            (0.81173349, 0.78518363, 0.03065289, 0.64191028),
        )
    )
    np.testing.assert_allclose(diffs, ref_diffs, atol=1e-7)
