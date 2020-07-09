import numpy as np

from pysisyphus.constants import ANG2BOHR, AU2KJPERMOL
from pysisyphus.calculators.Calculator import Calculator
from pysisyphus.calculators.LennardJones import LennardJones

# [1] https://aip.scitation.org/doi/abs/10.1063/1.445869
#     Jorgensen, 1983
# [2] https://pubs.acs.org/doi/pdf/10.1021/ja00344a001
#     Jorgensen, 1983
# Not yet implemented
#   OPC3, https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4991989/


class TIP3P(Calculator):
    """Transferable Intermolecular Potential 3 Point"""
    rOH = 0.9572 * ANG2BOHR
    aHOH = 104.52  # deg
    # Charges
    qH = +0.4170
    qO = -2*qH
    sigma = 3.15061 * ANG2BOHR
    epsilon = 0.6364 / AU2KJPERMOL

    # rc = 5 Å in Bohr
    def __init__(self, rc=9.44863062728914):
        super().__init__()

        # Cutoff distance
        self.rc = rc
        self.lj = LennardJones(sigma=self.sigma, epsilon=self.epsilon, rc=self.rc)

        self.charges = np.array((self.qO, self.qH, self.qH))
        """
        coulomb_energy = (multiple of elem. charge * multiple of elem. charge)
                         / (distance in bohr) * 1 / (4 * pi * vacuum permittivity)

        coulomb_prefactor converts everything to atomic units and it is ... drum
        roll ... 1.
        from scipy.constants import value as pcval
        self.coulomb_prefactor = (1 / (4 * np.pi) * pcval("elementary charge")**2
                                  / pcval("Hartree energy") / pcval("Bohr radius")
                                  / pcval("vacuum electric permittivity")
        )
        """
        self.coulomb_prefactor = 1

    def coulomb(self, coords3d):
        waters = coords3d.shape[0] // 3
        if waters == 1:
            return 0., np.zeros((3, 3))
        stencil = np.array((0, 1, 2))
        # Pair indices of the interacting water molecules
        a, b = np.triu_indices(waters, 1)

        # Pair indices of the respective interacting atoms
        # a: 0, 1, 2, 0, 1, 2, 0, 1, 2, ...
        a = (3*np.repeat(a, 3)[:, None] + stencil).flatten()
        # b: 3, 3, 3, 4, 4, 4, 5, 5, 5, ... 
        b = np.repeat((3 * b[:, None] + stencil).flatten(), 3)

        # Pair coordinate differences
        diffs = coords3d[a] - coords3d[b]  # Shape: (Interacting pairs, 3)
        # Distances
        rs = np.linalg.norm(diffs, axis=1)

        # Assemble water charges for all atoms/coordinates
        charges = np.tile(self.charges, coords3d.shape[0] // 3)
        # Keep the uncontracted pair-energies, as we still need them for the
        # gradient.
        pair_energies = self.coulomb_prefactor * charges[a] * charges[b] / rs
        energy = pair_energies.sum()

        # See the LennardJones calculator for the explicit derivation of the
        # derivative of 1/r**n.
        products = (pair_energies / rs**2)[:, None] * diffs

        gradient = np.zeros_like(coords3d)
        # Every pair (a, b) contributes to the total gradient of atoms a and b.
        for i, prod in enumerate(products):
            gradient[a[i]] -= prod
            gradient[b[i]] += prod

        return energy, -gradient

    def calculate(self, coords3d):
        coulomb_energy, coulomb_forces = self.coulomb(coords3d)

        # Lennard-Jones interaction between Oxygens only
        o_coords3d = coords3d[::3]
        lj_energy, lj_forces = self.lj.calculate(o_coords3d)

        forces = np.zeros_like(coords3d)
        forces[::3] += lj_forces

        forces = coulomb_forces.copy()
        # Add LennardJones forces to oxygens
        forces[::3] += lj_forces
        energy = coulomb_energy + lj_energy

        return energy, forces

    def get_energy(self, atoms, coords):
        energy, _ = self.calculate(coords.reshape(-1, 3))
        return {"energy": energy}

    def get_forces(self, atoms, coords):
        energy, forces = self.calculate(coords.reshape(-1, 3))
        return {"energy": energy,
                "forces": forces.flatten(),
        }
