#!/usr/bin/env python3

from ase.cluster.icosahedron import Icosahedron
from ase.calculators.lj import LennardJones as ase_LJ
import numpy as np
import pytest

from pysisyphus.calculators.LennardJones import LennardJones
from pysisyphus.constants import BOHR2ANG
from pysisyphus.Geometry import Geometry
from pysisyphus.helpers import geom_loader
from pysisyphus.optimizers.RFOptimizer import RFOptimizer


def test_lennard_jones():
    atoms = Icosahedron("Ar", noshells=2, latticeconstant=3)
    atoms.calc = ase_LJ()
    ase_forces = atoms.get_forces()
    ase_energy = atoms.get_potential_energy()

    coords = atoms.positions.flatten()
    geom = Geometry(atoms.get_chemical_symbols(), coords / BOHR2ANG)
    geom.set_calculator(LennardJones())

    pysis_energy = geom.energy
    assert pysis_energy == pytest.approx(ase_energy)

    pysis_forces = geom.forces / BOHR2ANG
    np.testing.assert_allclose(pysis_forces, ase_forces.flatten(), atol=1e-15)


def test_ar_cluster():
    geom = geom_loader("lib:ar14cluster.xyz")
    geom.set_calculator(LennardJones())

    opt_kwargs = {
        "max_cycles": 150,
        "gediis": True,
        "thresh": "gau_vtight",
    }
    opt = RFOptimizer(geom, **opt_kwargs)
    opt.run()

    assert opt.is_converged
    assert opt.cur_cycle == 110
