import pytest

from pysisyphus.calculators import XTB
from pysisyphus.dynamics.helpers import get_mb_velocities_for_geom
from pysisyphus.helpers import geom_loader
from pysisyphus.testing import using


@using("xtb")
def test_run_md():
    geom = geom_loader("lib:so3hcl_diss_ts_opt.xyz")
    calc = XTB(pal=2)

    T = 298.15
    seed = 20182503
    energy_ref = -20.46920283843
    velocities = get_mb_velocities_for_geom(geom, T=T, seed=seed)
    geoms = calc.run_md(geom.atoms, geom.cart_coords, t=200, dt=0.1,
                        velocities=velocities)

    assert len(geoms) == 200
    last_geom = geoms[-1]
    energy = calc.get_energy(last_geom.atoms, last_geom.cart_coords)["energy"]
    assert energy == pytest.approx(energy_ref)
