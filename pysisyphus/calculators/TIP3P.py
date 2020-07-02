import numpy as np

from pysisyphus.constants import ANG2BOHR, AU2KJPERMOL
from pysisyphus.calculators.Calculator import Calculator
from pysisyphus.calculators.LennardJones import LennardJones


class TIP3P(Calculator):
    rOH = 0.9572 * ANG2BOHR
    aHOH = 104.52  # deg
    # Charges
    qO = -0.8340
    qH = +0.4170
    sigma = 3.15061 * ANG2BOHR
    epsilon = 0.6364 / AU2KJPERMOL

    # rc = 5 Å in Bohr
    def __init__(self, rc=9.44863062728914):
        super().__init__()

        # Cutoff distance
        self.rc = rc
        self.lj = LennardJones(sigma=self.sigma, epsilon=self.epsilon, rc=self.rc)

    def calculate(self, coords3d):
        # Lennard-Jones interaction between Oxygens only
        o_coords3d = coords3d[::3]
        lj_energy, lj_forces = self.lj.calculate(o_coords3d)

        forces = np.zeros_like(coords3d)
        forces[::3] = lj_forces
        energy = lj_energy

        return energy, forces

    def get_energy(self, atoms, coords):
        energy, _ = self.calculate(coords.reshape(-1, 3))
        return {"energy": energy}

    def get_forces(self, atoms, coords):
        energy, forces = self.calculate(coords.reshape(-1, 3))
        return {"energy": energy,
                "forces": forces.flatten(),
        }
