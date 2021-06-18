import itertools as it

from pysisyphus.calculators import TIP3P, ExternalPotential
from pysisyphus.dynamics.helpers import get_mb_velocities_for_geom
from pysisyphus.dynamics import md
from pysisyphus.helpers import geom_loader


def get_water_constraints(index):
    i = index*3
    return [[i, i+1], [i, i+2], [i+1, i+2]]


def test_rattle_tip3p(this_dir):
    geom = geom_loader(this_dir / "output_10.xyz")
    # geom = geom_loader(this_dir / "output_2.xyz")
    # geom.jmol()

    T = 298.15
    calc = TIP3P()
    potentials = [
        {
            "type": "logfermi",
            "beta": 6,
            "T": T,
            "radius": 10,
        },
    ]
    ext_calc = ExternalPotential(calc, potentials=potentials)
    geom.set_calculator(ext_calc)

    constraints = list(it.chain(*[get_water_constraints(i)
                                  for i in range(len(geom.atoms) // 3)]))

    seed = 20200626
    v0 = get_mb_velocities_for_geom(geom, T, seed=seed).flatten()

    md_kwargs = {
        "v0": v0,
        # "steps": 100,
        # "dt": 1,
        "steps": 250,
        "dt": 1.5,
        # "steps": 500,
        # "dt": 0.25,
        "constraints": constraints,
        "constraint_kwargs": {
            # "remove_com_v": False,
        },
        # "thermostat": "csvr",
    }

    # import pdb; pdb.set_trace()
    md_result = md(geom, **md_kwargs)

    from pysisyphus.xyzloader import coords_to_trj
    coords = md_result.coords
    trj_fn = "md.trj"
    atoms = geom.atoms
    coords_to_trj(trj_fn, atoms, coords)
