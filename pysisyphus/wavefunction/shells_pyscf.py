import numpy as np

from pysisyphus.helpers_pure import argsort
from pysisyphus.wavefunction.helpers import permut_for_order
from pysisyphus.wavefunction.shells import Shells


def pyscf_cart_order(l):
    order = list()
    for lx in reversed(range(l + 1)):
        for ly in reversed(range(l + 1 - lx)):
            lz = l - lx - ly
            order.append("x" * lx + "y" * ly + "z" * lz)
    return tuple(order)


_CART_ORDER = [pyscf_cart_order(l) for l in range(5)]
_CART_PS = permut_for_order(_CART_ORDER)
_SPH_PS = {
    0: [[1]],  # s
    1: [[1, 0, 0], [0, 0, 1], [0, 1, 0]],  # px py pz
    2: [
        [0, 0, 0, 0, 1],  # dxy
        [0, 0, 0, 1, 0],  # dyz
        [0, 0, 1, 0, 0],  # dz²
        [0, 1, 0, 0, 0],  # dxz
        [1, 0, 0, 0, 0],  # dx² - y²
    ],
    3: [
        [0, 0, 0, 0, 0, 0, 1],
        [0, 0, 0, 0, 0, 1, 0],
        [0, 0, 0, 0, 1, 0, 0],
        [0, 0, 0, 1, 0, 0, 0],
        [0, 0, 1, 0, 0, 0, 0],
        [0, 1, 0, 0, 0, 0, 0],
        [1, 0, 0, 0, 0, 0, 0],
    ],
    4: [
        [0, 0, 0, 0, 0, 0, 0, 0, 1],
        [0, 0, 0, 0, 0, 0, 0, 1, 0],
        [0, 0, 0, 0, 0, 0, 1, 0, 0],
        [0, 0, 0, 0, 0, 1, 0, 0, 0],
        [0, 0, 0, 0, 1, 0, 0, 0, 0],
        [0, 0, 0, 1, 0, 0, 0, 0, 0],
        [0, 0, 1, 0, 0, 0, 0, 0, 0],
        [0, 1, 0, 0, 0, 0, 0, 0, 0],
        [1, 0, 0, 0, 0, 0, 0, 0, 0],
    ],
}


def get_pyscf_P(shells: Shells, cartesian: bool = True) -> np.ndarray:
    """Get Cartesian/spherical Pyscf-permutation matrix P.

      --- PySCF --->
      |
      |
    pysis    P
      |
      v

    Parameters
    ----------
    shells
        Pysisyphus Shells object.

    Returns
    -------
    P.T
        Permutation matrix. PySCF-order will be along columns, pysisyphus-order
        along the rows.
    """
    if cartesian:
        nbfs = shells.cart_size
        _PS = _CART_PS

        def size_func(L):
            return (L + 2) * (L + 1) // 2

    # Spherical basis functions
    else:
        nbfs = shells.sph_size
        _PS = _SPH_PS

        def size_func(L):
            return 2 * L + 1

    # Empty permutation matrix
    P = np.zeros((nbfs, nbfs), dtype=int)
    pysis_offset = 0
    pyscf_offset = 0
    for _, center_shells in shells.center_shell_iter():
        # Angular momenta of shells on a center
        center_Ls = [shell.L for shell in center_shells]
        # Their shell sizes
        sizes = list(map(size_func, center_Ls))
        # Total shell size on a center
        center_tot_size = sum(sizes)
        # Determine starting indices on pysisyphus' side as cumulative sum
        # of the respective shell sizes.
        pysis_starts = pysis_offset + np.cumsum(np.array([0] + sizes[:-1], dtype=int))
        # Determine the shell ordering in PySCF; there, shells on a center
        # are ordered by increasing angular momentum.
        pyscf_inds = argsort(center_Ls)

        # Set blocks in the permutatin matrix P
        #
        # Loop over shells in PySCF-order
        for i in pyscf_inds:
            # Determine, where shell i starts in pysisyphus and its angular momentum L
            pysis_start = pysis_starts[i]
            L = center_Ls[i]
            size = size_func(L)
            # Pick permutation-subblock w/ PySCF's ordering ...
            p = _PS[L]
            # ... and put it in the correct place in the permutation matrix.
            P[
                pysis_start : pysis_start + size,
                pyscf_offset : pyscf_offset + size,
            ] = p
            # Increase offsets and advance to the next shell
            pyscf_offset += size
        pysis_offset += center_tot_size
    # Permutation matrix
    return P.T


class PySCFShells(Shells):
    """
    Cartesian bfs >= d angular momentum are not normalized!
    S_ref = mol.intor("int1e_ovlp_cart")
    N = 1 / np.diag(S_ref)**0.5
    ao *= N
    """

    cart_Ps = _CART_PS
    sph_Ps = _SPH_PS

    @property
    def P_cart(self):
        """Permutation matrix for Cartesian basis functions."""
        if self.ordering == "pysis":
            return np.eye(self.cart_size)
        elif self._P_cart is None:
            self._P_cart = get_pyscf_P(self, cartesian=True)
        return self._P_cart

    @property
    def P_sph(self):
        """Permutation matrix for Spherical basis functions."""
        if self.ordering == "pysis":
            return np.eye(self.sph_size)
        elif self._P_sph is None:
            self._P_sph = get_pyscf_P(self, cartesian=False)
        return self._P_sph
