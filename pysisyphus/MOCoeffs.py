import dataclasses
from typing import Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np

from pysisyphus.constants import AU2EV


def occ_from_occs(occs):
    occ = occs.sum()
    iocc = int(occ)
    # Don't allow fractional total occupations
    assert float(iocc) == occ
    return iocc


def filter_ens(
    ens: np.ndarray, en_homo: float, below_homo: float, above_homo: float
) -> Tuple[np.ndarray, np.ndarray]:
    below_thresh = en_homo - below_homo
    above_thresh = en_homo + above_homo
    ens_between = list()
    indices_between = list()
    for i, en in enumerate(ens):
        if below_thresh <= en <= above_thresh:
            ens_between.append(en)
            indices_between.append(i)
    ens_between = np.array(ens_between)
    indices_between = np.array(indices_between, dtype=int)
    return ens_between, indices_between


@dataclasses.dataclass
class MOCoeffs:
    # Ca/Cb: 2d array of floats w/ shape (naos, nmos) containing MO coefficients
    # ensa/ensb: 1d array of floats w/ shape (nmos, ) containing MO energies
    # occsa/occs: 1d array of floats w/ shape (nmos, ) containing MO occupation numbers
    Ca: np.ndarray
    ensa: np.ndarray
    occsa: np.ndarray
    # Beta electron/MO quantities are optional, as they are not present in restricted
    # calculations. When only alpha quantities are given all beta quantities will be
    # derived/copied from them.
    Cb: Optional[np.ndarray] = None
    ensb: Optional[np.ndarray] = None
    occsb: Optional[np.ndarray] = None

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

    def _homo(self, occ):
        return occ - 1 if occ else None

    @property
    def homoa(self) -> Optional[int]:
        return self._homo(self.occa)

    @property
    def homob(self) -> Optional[int]:
        return self._homo(self.occb)

    @property
    def lumoa(self):
        return self.occa

    @property
    def lumob(self):
        return self.occb

    def _virt_inds(self, occs):
        return np.arange(occs.size)[occs == 0.0]

    @property
    def virt_indsa(self):
        return self._virt_inds(self.occsa)

    @property
    def virt_indsb(self):
        return self._virt_inds(self.occsb)

    def swap_inplace(
        self,
        C: np.ndarray,
        ens: np.ndarray,
        occs: np.ndarray,
        ind1: int,
        ind2: int,
        swap_energies: bool = True,
        swap_occs: bool = True,
    ):
        """Swap a pair of MO coeffs and energies inplace."""
        tmp = C[:, ind1].copy()
        C[:, ind1] = C[:, ind2]
        C[:, ind2] = tmp
        if swap_energies:
            ens[ind1], ens[ind2] = ens[ind2], ens[ind1]
        if swap_occs:
            occs[ind1], occs[ind2] = occs[ind2], occs[ind1]

    def swap_mos(self, alpha_pairs, beta_pairs=None, **kwargs):
        if beta_pairs is None:
            beta_pairs = list()

        new_kwargs = dataclasses.asdict(self)

        Ca = new_kwargs["Ca"]
        ensa = new_kwargs["ensa"]
        occsa = new_kwargs["occsa"]
        for ind1, ind2 in alpha_pairs:
            self.swap_inplace(Ca, ensa, occsa, ind1, ind2, **kwargs)

        Cb = new_kwargs["Cb"]
        ensb = new_kwargs["ensb"]
        occsb = new_kwargs["occsa"]
        for ind1, ind2 in beta_pairs:
            self.swap_inplace(Cb, ensb, occsb, ind1, ind2, **kwargs)
        return MOCoeffs(**new_kwargs)

    def plot_mo_energies(self, below_homo=0.5, above_homo=0.5, show=False):
        ensa = self.ensa
        homoa = self.homoa
        en_homoa = ensa[homoa]
        ensb = self.ensb
        homob = self.homob
        en_homob = ensa[homob]

        ensa, indsa = filter_ens(ensa, en_homoa, below_homo, above_homo)
        ensb, indsb = filter_ens(ensb, en_homob, below_homo, above_homo)
        occsa = self.occsa[indsa]
        occsb = self.occsb[indsb]

        def colors(occs):
            return ["red" if occ else "blue" for occ in occs]

        colorsa = colors(occsa)
        colorsb = colors(occsb)

        en_homoa = en_homoa * AU2EV
        en_homob = en_homob * AU2EV
        ensa = ensa * AU2EV
        ensb = ensb * AU2EV
        xa = np.ones_like(ensa)
        xb = np.ones_like(ensb) + 1

        fig, ax = plt.subplots()

        def annot(xs, ens, inds):
            for x, en, ind in zip(xs, ens, inds):
                text = str(ind)
                xy = (x + 0.125, en)
                ax.annotate(text, xy)

        kwargs = {
            "s": 200,
            "marker": "_",
            "zorder": 3,
        }
        ax.scatter(xa, ensa, c=colorsa, label="α", **kwargs)
        ax.scatter(xb, ensb, c=colorsb, label="β", **kwargs)
        annot(xa, ensa, indsa)
        annot(xb, ensb, indsb)
        ax.axhline(en_homoa, c="k", ls="--", label="HOMO α")
        ax.axhline(en_homob, c="k", ls="--", label="HOMO β")
        ax.legend()
        ax.set_ylabel("E / eV")
        ax.set_xlim(0, 3)
        fig.tight_layout()
        if show:
            plt.show()
        return fig, ax
