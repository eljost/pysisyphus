import dataclasses
from typing import Optional

import numpy as np


def occ_from_occs(occs):
    occ = occs.sum()
    iocc = int(occ)
    # Don't allow fractional total occupations
    assert float(iocc) == occ
    return iocc


@dataclasses.dataclass
class MOCoeffs:
    Ca: np.ndarray
    ensa: np.ndarray
    occsa: int
    Cb: Optional[np.ndarray] = None
    ensb: Optional[np.ndarray] = None
    occsb: Optional[int] = None

    def __post_init__(self):
        self._restricted = self.Cb is None

        # All beta quantities must be given together
        if self.unrestricted:
            assert self.ensb is not None and self.occsb is not None
        else:
            self.Cb = self.Ca.copy()
            self.ensb = self.ensa.copy()
            self.occsa = self.occsa / 2.0
            self.occsb = self.occsa.copy()

        self._occa = occ_from_occs(self.occsa)
        self._occb = occ_from_occs(self.occsb)

    @property
    def restricted(self) -> bool:
        return self._restricted

    @property
    def unrestricted(self) -> bool:
        return not self.restricted

    @property
    def occa(self) -> int:
        return self._occa

    @property
    def occb(self) -> int:
        return self._occb

    def swap_inplace(self, C: np.ndarray, ens: np.ndarray, ind1: int, ind2: int):
        """Swap a pair of MO coeffs and energies inplace."""
        tmp = C[:, ind1].copy()
        C[:, ind1] = C[:, ind2]
        C[:, ind2] = tmp
        ens[ind1], ens[ind2] = ens[ind2], ens[ind1]

    def swap_mos(self, alpha_pairs, beta_pairs=None):
        if beta_pairs is None:
            beta_pairs = list()

        new_kwargs = dataclasses.asdict(self)

        Ca = new_kwargs["Ca"]
        ensa = new_kwargs["ensa"]
        for ind1, ind2 in alpha_pairs:
            self.swap_inplace(Ca, ensa, ind1, ind2)

        Cb = new_kwargs["Cb"]
        ensb = new_kwargs["ensb"]
        for ind1, ind2 in beta_pairs:
            self.swap_inplace(Cb, ensb, ind1, ind2)
        return MOCoeffs(**new_kwargs)
