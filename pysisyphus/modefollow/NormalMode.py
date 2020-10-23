from math import sqrt

import numpy as np


class NormalMode:
    """See http://gaussian.com/vib/"""

    def __init__(self, l, masses):
        """NormalMode class.

        Cartesian displacements are normalized to 1.

        Parameters
        ----------
        l : np.array
            Cartesian, non-mass-weighted displacements.
        masses : np.array
            Atomic masses.
        """

        self.l = np.array(l.flatten())
        self.l /= np.linalg.norm(l)
        self.masses = masses
        assert self.l.shape == self.masses.shape

    def __len__(self):
        return self.l.size

    @property
    def red_mass(self):
        return 1 / np.sum(np.square(self.l_mw) / self.masses)

    def mw_norm_for_norm(self, norm=0.01):
        return norm * sqrt(self.red_mass)

    @property
    def l_mw(self):
        l_mw = self.l * np.sqrt(self.masses)
        return l_mw / np.linalg.norm(l_mw)
