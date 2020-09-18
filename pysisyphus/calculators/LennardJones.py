import numpy as np

from pysisyphus.calculators.Calculator import Calculator


class LennardJones(Calculator):

    # Corresponds to σ = 1 Å, as the default value in ASE, but
    # pysisyphus uses au/Bohr.
    def __init__(self, sigma=1.8897261251, epsilon=1, rc=None):
        super().__init__()

        self.sigma = sigma
        self.epsilon = epsilon
        # Cutoff distance
        if rc is None:
            rc = 3 * self.sigma
        self.rc = rc
        # Shift energy
        self.e0 = (4 * self.epsilon *
                   ((self.sigma/self.rc)**12 - (self.sigma/self.rc)**6)
        )

    def calculate(self, coords3d):
        # Index pairs
        a, b = np.triu_indices(len(coords3d), 1)

        # Distances
        diffs = coords3d[a] - coords3d[b]  # Shape: (N_atoms, 3)
        rs = np.linalg.norm(diffs, axis=1)
        c6 = np.where(rs <= self.rc, (self.sigma/rs)**6, np.zeros_like(rs))
        energy = -self.e0 * (c6 != 0.0).sum()
        c12 = c6**2
        energy += np.sum(4*self.epsilon * (c12 - c6))

        """
        Gradient related remarks:

        Lennard-Jones-potential:
            LJ(r) =
                4*ε*[(σ/r)**12 - (σ/r)**6]

        Derivative of quotient appearing in LJ potential w.r.t first cartesian
        coordinate:
            d(σ/r)**n/dx_1
                = σ**n * d/dx_1 r**(-n)
                = σ**n * (-n/2) * (2x1-2x2) * r**(-n) * r**(-2)
                = σ**n * r**(-n) * (-n) * (x1-x2) * r**(-2)
                = (σ/r)**n * (-n) * (x1-x2) / r**2

        Derivate w.r.t to cartesian x coordinate of atom A (x_1):
            dLJ(r)/dx_1 = 
                24*ε*[-2*(σ/r)**12 + (σ/r)**6]*(x1-x2)/r**2

        Derivate w.r.t to cartesian x coordinate of atom B (x_2):
            dLJ(r)/dx_2 = 
                -24*ε*[-2*(σ/r)**12 + (σ/r)**6]*(x1-x2)/r**2

        The derivate w.r.t to x_2 differs only by a factor of -1!
        """
        prefactors = 24*self.epsilon * (c6 - 2*c12) / rs**2
        products = prefactors[:,None] * diffs

        gradient = np.zeros_like(coords3d)
        # Every pair (a, b) contributes to the total gradient of atoms a and b.
        for i, prod in enumerate(products):
            gradient[a[i]] += prod
            gradient[b[i]] -= prod

        return energy, -gradient

    def get_energy(self, atoms, coords):
        energy, _ = self.calculate(coords.reshape(-1, 3))
        return {"energy": energy}

    def get_forces(self, atoms, coords):
        energy, forces = self.calculate(coords.reshape(-1, 3))
        return {"energy": energy,
                "forces": forces.flatten(),
        }
