import numpy as np
import pytest

from pysisyphus.calculators import XTB
from pysisyphus.constants import AU2NU, AU2SEC
from pysisyphus.finite_diffs import periodic_fd_2_8
from pysisyphus.helpers import geom_loader
from pysisyphus.hindered_rotor import torsion_driver
from pysisyphus.hindered_rotor.opt import spline_closure
from pysisyphus.testing import using


@using("xtb")
def test_torsion_gpr_ethane(this_dir):
    geom = geom_loader(this_dir / "ethane_xtbopt.xyz")
    torsion = (4, 0, 1, 6)

    def calc_getter(**kwargs):
        calc_kwargs = {
            "keep_kind": "latest",
            "charge": 0,
            "mult": 1,
        }
        calc_kwargs.update(**kwargs)

        calc = XTB(acc=0.001, **calc_kwargs)
        return calc

    result = torsion_driver.run(
        geom,
        torsion,
        calc_getter=calc_getter,
        plot=True,
        out_dir=this_dir / "ethane",
        partfunc_thresh=1e-3,
    )
    np.testing.assert_allclose(
        (result.hr_partfunc, result.cancel_partfunc),
        # For partfunc_thresh=1e-3; values for xtb 6.7.1
        (4.283, 0.670),
        atol=1e-2,
    )


def test_h2o2_g16_analytical(this_dir):
    geom = geom_loader(this_dir / "h2o2_b3lyp_tzvp_opt.xyz")
    fn = this_dir / "h2o2_g16_b3_tz_relaxed_scan.dat"
    data = np.loadtxt(fn)
    energy_getter = spline_closure(
        4,
        *data.T,
    )
    result = torsion_driver.run(
        geom,
        [1, 3, 2, 0],
        energy_getter=energy_getter,
        out_dir=this_dir / "h2o2_spline",
    )
    w = result.eigvals
    tunnel_splitting = (w[1] - w[0]) * AU2NU
    assert tunnel_splitting == pytest.approx(12.76, abs=5e-2)


def test_cancel_partfunc(this_dir):
    fn = this_dir / "h2o2_g16_b3_tz_relaxed_scan.dat"
    rads, energies = np.loadtxt(fn).T
    energies -= energies.min()
    dx = rads[1] - rads[0]
    force_constant = periodic_fd_2_8(0, energies, dx)
    imom = 2730.6925715477555
    ho_freq = np.sqrt(force_constant / imom) / (2 * np.pi)
    ho_freq_si = ho_freq / AU2SEC
    temperature = 298.15
    from pysisyphus.partfuncs import partfuncs as pf

    cancel_partfunc = pf.harmonic_quantum_partfunc(ho_freq_si, temperature)
    print(f"{cancel_partfunc=: >10.5f}")
    assert cancel_partfunc == pytest.approx(0.05015, abs=1e-5)
