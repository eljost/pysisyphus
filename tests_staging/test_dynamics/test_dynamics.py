import numpy as np

from pysisyphus.constants import BOHR2ANG
from pysisyphus.calculators.AnaPot import AnaPot
from pysisyphus.calculators.MullerBrownSympyPot import MullerBrownPot
from pysisyphus.calculators.XTB import XTB
from pysisyphus.dynamics import mdp
from pysisyphus.dynamics.velocity_verlet import md
from pysisyphus.helpers import geom_loader


def test_so3hcl_diss():
    """See [1]"""
    def get_geom():
        geom = geom_loader("lib:so3hcl_diss_ts_opt.xyz")
        geom.set_calculator(XTB(pal=4))
        return geom

    geom = get_geom()
    mdp_kwargs = {
        # About 5 kcal/mol
        "E_excess": 0.0079,
        "term_funcs": None,
        "epsilon": 5e-4,
        "ascent_alpha": 0.05,
        "t_init": 20,
        # Paper uses 200
        "t": 100,
        "dt": .5,
        "seed": 25032018,
        # "external_md": True,
        "max_init_trajs": 1,
    }
    res = mdp(geom, **mdp_kwargs)

    import pdb; pdb.set_trace()
    # geom = get_geom()
    # mdp_kwargs["E_excess"] = 0
    # res_ee = mdp(geom, **mdp_kwargs)


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
