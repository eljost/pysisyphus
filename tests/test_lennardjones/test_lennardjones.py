#!/usr/bin/env python3

from ase.cluster.icosahedron import Icosahedron
from ase.calculators.lj import LennardJones as ase_LJ
import numpy as np
import pytest

from pysisyphus.calculators.LennardJones import LennardJones
from pysisyphus.Geometry import Geometry


def test_lennard_jones():
    atoms = Icosahedron("Ar", noshells=2, latticeconstant=3)
    atoms.set_calculator(ase_LJ())
    ase_forces = atoms.get_forces()

    coords = atoms.positions.flatten()
    geom = Geometry(atoms.get_chemical_symbols(), coords)
    geom.set_calculator(LennardJones())
    np.testing.assert_allclose(geom.forces, ase_forces.flatten(), atol=1e-15)
