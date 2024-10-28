import numpy as np
import pytest

from pysisyphus.calculators import XTB
from pysisyphus.constants import AU2NU
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
        en_thresh=1e-6,
        plot=True,
        out_dir=this_dir / "ethane",
    )
    np.testing.assert_allclose(
        (result.hr_partfunc, result.cancel_partfunc),
        (4.285, 1.051),
        atol=2e-2,
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
    print(f"{tunnel_splitting: >6.2f} cm⁻¹")
    assert tunnel_splitting == pytest.approx(12.71, abs=1e-1)
