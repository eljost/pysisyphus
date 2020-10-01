import numpy as np
from scipy.spatial.distance import pdist, squareform

from pysisyphus.calculators.Calculator import Calculator


class IDPPCalculator(Calculator):

    def __init__(self, target): 
        self.target = squareform(target)

        super().__init__(base_name="idpp")

    def get_forces(self, atoms, coords):
        coords_reshaped = coords.reshape((-1, 3))

        D = []
        for c in coords_reshaped:
            Di = coords_reshaped - c
            D.append(Di)
        D = np.array(D)

        curr_pdist = pdist(coords_reshaped)
        curr_square = squareform(curr_pdist)
        curr_diff = curr_square - self.target

        curr_square = curr_square + np.eye(curr_square.shape[0])

        # The bigger the differences 'curr_diff', the bigger the energy.
        # The smaller the current distances 'current_pdist', the bigger
        # the energy.
        energy = 0.5 * (curr_diff**2 / curr_square**4).sum()

        # Adapted from ASE IDPP calculator
        # https://gitlab.com/ase/ase/blob/master/ase/neb.py, GPL2
        forces = -2 * ((curr_diff *
                       (1 - 2 * curr_diff / curr_square) /
                        curr_square**5)[...,np.newaxis] * D).sum(0)

        results = {
            "energy" : energy,
            "forces": forces.flatten()
        }
        return results

    def __str__(self):
        return "IDPP calculator"


