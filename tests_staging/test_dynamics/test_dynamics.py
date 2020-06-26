import numpy as np

from pysisyphus.constants import BOHR2ANG
from pysisyphus.dynamics.velocity_verlet import md
from pysisyphus.helpers import geom_loader


def test_oniom_md():
    calc_dict = {
        "high": {
            "type": "pypsi4",
            "method": "scf",
            "basis": "sto-3g",
        },
        "low": {
            "type": "pyxtb",
        },
    }
    high_inds = (4,5,6)
    from pysisyphus.calculators.ONIOM import ONIOM
    oniom = ONIOM(calc_dict, high_inds)

    geom = geom_loader("lib:acetaldehyd_oniom.xyz")
    geom.set_calculator(oniom)

    v0 = .005 * np.random.rand(*geom.coords.shape)
    md_kwargs = {
        "v0": v0,
        "t": 40,
        "dt": 0.5,
    }
    md_result = md(geom, **md_kwargs)
    from pysisyphus.xyzloader import make_trj_str

    coords = md_result.coords.reshape(-1, len(geom.atoms), 3) * BOHR2ANG
    trj_str = make_trj_str(geom.atoms, coords)
    with open("md.trj", "w") as handle:
        handle.write(trj_str)
