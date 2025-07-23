# [1] https://doi.org/10.1063/1.3463717
#     Comparisons of classical and Wigner sampling of transition state energy levels
#     for quasiclassical trajectory chemical dynamics simulations
#     Sun, Hase, 2010

import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

from pysisyphus.dynamics import wigner


def test_h2_wigner_sampling(this_dir):
    h5_fn = this_dir / "h2_hfdef2svp_final_hessian.h5"
    seed = 20182020
    sampler, _ = wigner.get_wigner_sampler(str(h5_fn), temperature=300, seed=seed)
    n = 500
    displ_coords3d = np.zeros((n, 2, 3))
    for i in range(n):
        wsample = sampler()
        displ_coords3d[i] = wsample.coords3d
    dists = np.linalg.norm(displ_coords3d[:, 0] - displ_coords3d[:, 1], axis=1)
    descr = stats.describe(dists)
    min_, max_ = descr.minmax
    mean = descr.mean
    var = descr.variance
    np.testing.assert_allclose(
        (min_, max_, mean, var), (1.01207, 1.88892, 1.41359, 0.02420), atol=1e-5
    )


def test_c2h5f_ts_wigner_sampling(this_dir):
    """Recreation of Fig. 2 from [1].

    The transition state was optimized using ORCA 5.0.3 at the RI-MP2/6-31(d)
    level of theory using the (tightscf, def2-tzvp/c, nofrozencore) keywords.
    """

    h5_fn = this_dir / "c2h5f_631gdrimp2_ts_final_hessian.h5"
    seed = 20182020
    T = 300
    sampler, _ = wigner.get_wigner_sampler(str(h5_fn), temperature=T, seed=seed)
    n = 400
    displ_coords3d = np.zeros((n, 8, 3))
    qs_nodim = np.zeros((n, 18))
    ps_nodim = np.zeros((n, 18))
    qs = np.zeros((n, 18))
    for i in range(n):
        wsample = sampler()
        displ_coords3d[i] = wsample.coords3d
        qs[i] = wsample.qs
        qs_nodim[i] = wsample.qs_nodim
        ps_nodim[i] = wsample.ps_nodim

    nu470 = 1

    """
    fig, ax = plt.subplots()
    ax.violinplot(qs, quantiles=[[0.1, 0.9]] * 18)
    plt.show()

    nu3340 = 17
    fig, ((ax0, ax1), (ax2, ax3)) = plt.subplots(nrows=2, ncols=2, figsize=(12, 12))
    lim = 4
    bins = np.linspace(-lim, lim, num=25)
    ax0.hist(qs_nodim[:, nu470], bins=bins)
    ax0.set_xlabel("Q'")
    ax0.set_title("$C_2H_4$ rotation (470 cm⁻¹)")
    ax1.hist(ps_nodim[:, nu470], bins=bins)
    ax1.set_xlabel("P'")
    ax1.set_title("C$_2$H$_4$ rotation (470 cm⁻¹)")

    ax2.hist(qs_nodim[:, nu3340], bins=bins)
    ax2.set_xlabel("Q'")
    ax2.set_title("CH stretch (3340 cm⁻¹)")
    ax3.hist(ps_nodim[:, nu3340], bins=bins)
    ax3.set_title("CH stretch (3340 cm⁻¹)")
    ax3.set_xlabel("P'")

    for ax in (ax0, ax1, ax2, ax3):
        ax.set_xlim(-lim, lim)
        ax.set_ylabel("Count")
        ax.axvline(0.0, c="k", ls=":")
    fig.tight_layout()
    plt.show()
    """

    descr = stats.describe(qs_nodim[:, nu470])
    min_, max_ = descr.minmax
    mean = descr.mean
    var = descr.variance
    np.testing.assert_allclose(
        (min_, max_, mean, var), (-2.14725, 2.12394, -0.06041, 0.56758), atol=1e-5
    )
