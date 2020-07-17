import itertools as it

import pytest

from pysisyphus.calculators import TIP3P
from pysisyphus.dynamics.helpers import get_mb_velocities_for_geom
from pysisyphus.dynamics.velocity_verlet import md
from pysisyphus.dynamics import rattle_closure
from pysisyphus.helpers import geom_loader

from pysisyphus.dynamics.lincs import lincs_closure


def get_water_constraints(index):
    i = index*3
    return [[i, i+1], [i, i+2], [i+1, i+2]]


def test_rattle_tip3p(this_dir):
    # geom = geom_loader(this_dir / "output_10.xyz")
    geom = geom_loader(this_dir / "output_2.xyz")
    geom.set_calculator(TIP3P())
    # geom.jmol()

    constraints = list(it.chain(*[get_water_constraints(i) for i in range(len(geom.atoms) // 3)]))
    # constraints.append([2, 3])

    forces = geom.forces

    T = 298.15
    seed = 20200626
    velocities = get_mb_velocities_for_geom(geom, T, seed=seed).flatten()

    rattle = rattle_closure(geom, constraints, dt=0.5, tol=1e-2,
                            max_cycles=10)
    c, v, f = rattle(geom.coords, velocities, forces)
    print("coords")
    print(c)
    print("velocities")
    print(v)
    print("forces")
    print(f)
