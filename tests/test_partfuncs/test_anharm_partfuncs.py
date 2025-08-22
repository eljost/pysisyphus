import matplotlib.pyplot as plt
import numpy as np
import pytest

from pysisyphus.calculators import XTB
from pysisyphus.constants import AU2KJPERMOL
from pysisyphus.io import geom_from_hessian
import pysisyphus.partfuncs.driver as pf_driver
from pysisyphus.TablePrinter import TablePrinter
from pysisyphus.testing import using


def plot_fit(dists, energies):
    poly = pf_driver.fit_quartic_poly(dists, energies)

    dists_fit = np.linspace(dists[0], dists[-1], num=100)
    energies_fit = poly(dists_fit)

    fig, ax = plt.subplots(figsize=(10, 5))
    energy0 = energies.min()

    energies -= energy0
    energies_fit -= energy0
    energies *= AU2KJPERMOL
    energies_fit *= AU2KJPERMOL

    ax.plot(dists, energies, "o-", label="calc")
    ax.plot(dists_fit, energies_fit, ":", label="Fit")
    with np.printoptions(precision=3):
        ax.set_title(poly)
    ax.set_xlabel("Displacement / Bohr")
    ax.set_ylabel("ΔE / kJ mol⁻¹")
    ax.legend()
    ax.set_ylim(0, 60)
    fig.tight_layout()
    # plt.show()
    return fig


@using("xtb")
def test_anharm_partfuncs(this_dir):
    fn = this_dir / "ts_final_hessian_xtb.h5"
    calc = XTB(acc=0.01, mult=2)

    geom = geom_from_hessian(fn, calculator=calc)

    npoints = 10
    step_size = 0.4
    temperature = 648
    # nu_thresh = 100
    nu_thresh = 25
    all_part_funcs = pf_driver.run(
        geom,
        temperature,
        npoints=npoints,
        step_size=step_size,
        nu_thresh=nu_thresh,
        out_dir=this_dir,
    )

    print()
    header = (
        "#",
        "nu / cm⁻¹",
        "nu_aq / cm⁻¹",
        "q_aq",
        "q_hq",
        "U_aq / E_h",
        "S_sq / E_h·K⁻¹",
    )
    fmts = (
        "{:2d}",
        "{:6.2f}",
        "{:6.2f}",
        "{:10.6f}",
        "{:10.6f}",
        "{:12.4e}",
        "{:12.4e}",
    )
    table = TablePrinter(header, fmts, width=14, sub_underline=False)
    table.print_header()

    for pf in all_part_funcs:
        table.print_row(
            (
                pf.nu_ind,
                pf.nu,
                pf.nu_aq,
                pf.q_aq,
                pf.q_hq,
                pf.U_aq,
                pf.S_aq,
            )
        )
    ref1 = (15.99, 26.42, 10.144578, 28.158026)
    ref2 = (21.18, 35.95, 8.085407, 21.259538)
    apf1, apf2, *_ = all_part_funcs
    for apf, ref in ((apf1, ref1), (apf2, ref2)):
        nu_ref, nu_aq_ref, q_aq_ref, q_hq_ref = ref
        assert apf.nu == pytest.approx(nu_ref, abs=1e-2)
        assert apf.nu_aq == pytest.approx(nu_aq_ref, abs=1e-2)
        assert apf.q_aq == pytest.approx(q_aq_ref, abs=1e-6)
        assert apf.q_hq == pytest.approx(q_hq_ref, abs=1e-6)


def test_from_poly():
    coeffs = (-60.543, 3.651e-6, 4.953e-4, -7.419e-5, 4.203e-3)
    poly = np.polynomial.Polynomial(coeffs, domain=[-4, 4], window=[-1, 1])
    # with np.printoptions(precision=3):
    # print(poly)
    pf = pf_driver.calculate_partfuncs(
        nu=15.99,
        nu_ind=1,
        red_mass_au=9590.436648051831,
        poly=poly,
        temperature=648,
        natoms=41,
    )
    assert pf.q_wk == pytest.approx(10.1445, abs=1e-4)
    assert pf.q_aq == pytest.approx(10.1445, abs=1e-4)
    assert pf.q_ac == pytest.approx(10.1525, abs=1e-4)
    assert pf.q_hc == pytest.approx(28.1665, abs=1e-4)
    assert pf.q_hq == pytest.approx(28.1650, abs=1e-4)
    assert pf.dlnq_aq_dT == pytest.approx(0.0012003, abs=1e-7)
    assert pf.U_aq == pytest.approx(0.0015962, abs=1e-7)
    assert pf.S_aq == pytest.approx(9.800537e-6)
