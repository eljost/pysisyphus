import numpy as np
import pytest

from pysisyphus.constants import KB
from pysisyphus.calculators.ExternalPotential import ExternalPotential
from pysisyphus.calculators.LennardJones import LennardJones
from pysisyphus.Geometry import Geometry
from pysisyphus.helpers import geom_from_library


def test_external_potential():
    atoms = ("X", )
    coords = (0., 0., 0.)
    geom = Geometry(atoms, coords)

    calc = LennardJones()
    geom.set_calculator(calc)

    ref_energy = geom.energy
    potentials = [
        {
            "type": "logfermi",
            "beta": 6,
            "T": 1/KB,  # will result in kT = 1
            "radius": 10,
        },
    ]
    ext_calc = ExternalPotential(calc, potentials=potentials)
    geom.set_calculator(ext_calc)
    ext_energy = geom.energy
    assert ext_energy == pytest.approx(ref_energy)

    displ_coords = (9., 0., 0)
    geom.coords = displ_coords
    ext_energy = geom.energy
    # from math import exp, log; log(1 + exp(6*(9-10)))
    assert ext_energy == pytest.approx(0.002475685)

    forces = geom.forces
    ref_forces = np.array((-0.01483574, 0., 0.))
    np.testing.assert_allclose(forces, ref_forces)


def test_geometry_sphere_radius():
    geom = geom_from_library("benzene.xyz")
    radius = geom.get_sphere_radius(offset=4)

    assert radius == pytest.approx(8.7142660)
