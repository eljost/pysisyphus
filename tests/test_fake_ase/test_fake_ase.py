from ase.calculators.lj import LennardJones as ase_LJ
from ase.cluster.icosahedron import Icosahedron
from ase.optimize import BFGS
import numpy as np
import pytest

from pysisyphus.calculators import FakeASE, LennardJones


def test_fake_ase():
    # pysisyphus
    pysis_calc = LennardJones()
    ase_calc = FakeASE(pysis_calc)

    # ASE atoms, pysisyphus calculator
    atoms_pysis = Icosahedron("Ar", noshells=2, latticeconstant=3)
    # ASE atoms, ASE calculator
    atoms_ase = atoms_pysis.copy()

    atoms_pysis.calc = ase_calc
    atoms_ase.calc = ase_LJ()

    pysis_forces = atoms_pysis.get_forces()
    ase_forces = atoms_ase.get_forces()
    np.testing.assert_allclose(pysis_forces, ase_forces, atol=1e-15)


def test_fake_ase_opt():
    atoms = Icosahedron("Ar", noshells=2, latticeconstant=3)
    atoms.calc = FakeASE(LennardJones())

    dyn = BFGS(atoms)
    dyn.run(fmax=0.0005)

    assert dyn.converged()
    assert dyn.get_number_of_steps() == 14
    assert np.linalg.norm(dyn.f0) == pytest.approx(0.0041871980)
