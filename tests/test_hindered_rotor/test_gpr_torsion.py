import numpy as np

from pysisyphus.calculators import XTB
from pysisyphus.helpers import geom_loader
from pysisyphus.hindered_rotor import torsion_driver
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
        geom, calc_getter, torsion, plot=True, out_dir=this_dir / "ethane"
    )
    # print(f"{result.hr_partfunc}")
    # print(f"{result.cancel_partfunc}")
    np.testing.assert_allclose(
        (result.hr_partfunc, result.cancel_partfunc),
        (4.284727, 1.0424495),
        atol=1e-4,
    )
