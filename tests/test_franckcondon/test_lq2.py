import numpy as np
import pytest

from pysisyphus.constants import NU2AU
from pysisyphus.franckcondon.helpers import nu2angfreq_au
from pysisyphus.franckcondon import lq2_abs_cross_sec, unitless_displs_from_eigensystem
from pysisyphus.io import geom_from_hessian
from pysisyphus.helpers import eigval_to_wavenumber


def test_lq2(this_dir):
    fosc = 0.075871684
    dE_vert = 38864.1  # S1
    geom = geom_from_hessian(this_dir / "naphthalene_hessian.h5")
    grad = np.loadtxt(this_dir / "naphthalene_gradient")

    Hmw = geom.mw_hessian
    Hproj = geom.eckart_projection(Hmw, full=True)
    eigvals, eigvecs = np.linalg.eigh(Hproj)
    mask = np.abs(eigvals) > 1e-5

    significant_eigvals = eigvals[mask]
    significant_eigvecs = eigvecs[:, mask]
    # Wavenumbers in cm⁻¹
    nus = eigval_to_wavenumber(significant_eigvals)
    # Angular frequencies in au
    angfreqs = nu2angfreq_au(nus)
    mw_grad = grad / np.sqrt(geom.masses_rep)
    deltas = unitless_displs_from_eigensystem(
        mw_grad, significant_eigvals, significant_eigvecs
    )
    # HWHM in au
    hwhm = np.sqrt(np.log(2) * np.sum(deltas**2 * angfreqs**2))
    # HWHM in wavenumbers
    hwhm = hwhm / NU2AU

    Es = np.linspace(34_000, 50_000, 200)
    Es_au = Es * NU2AU
    dE_vert_au = dE_vert * NU2AU
    sigmas = lq2_abs_cross_sec(Es_au, fosc, dE_vert_au, angfreqs, deltas)
    sigmas /= sigmas.max()

    # xref, yref = np.loadtxt("/home/johannes/tmp/431_fc/naph_asa.abs.dat").T
    # import matplotlib.pyplot as plt
    # fig, ax = plt.subplots()
    # # ax.plot(Es, sigmas)
    # ax.plot(Es + hwhm, sigmas)
    # ax.axvline(dE_vert, c="k")
    # ax.plot(xref, yref)
    # plt.show()
    assert sigmas[10] == pytest.approx(0.127434001)
    assert sigmas[100] == pytest.approx(0.28347098)
