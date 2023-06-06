from typing import Dict, Optional

import numpy as np


ArrDict = Dict[str, np.ndarray]


class DIIS:
    def __init__(self, arr_specs: ArrDict):
        """DIIS container.

        Parameters
        ----------
        arr_specs
            Dict with string-keys and array-values. The key 'err_vecs' must always be
            present. Additional keys can be provided, depending on what must be stored
            next to the error vectors.
            The arrays should be at least 2d, with the first dimension having size 'nvecs'.
        """
        self.nvecs = arr_specs["err_vecs"].shape[0]
        for k, v in arr_specs.items():
            assert v.shape[0] == self.nvecs
            setattr(self, k, v)

        self.ref_keys = set(arr_specs.keys())
        self.cur_nvecs = 0
        self.cur_index = 0

    def store(self, arrs: ArrDict) -> int:
        """Store arrays.

        Returns the index, where the arrays were stored in their respective arrays.
        """
        assert set(arrs.keys()) == self.ref_keys
        index = self.cur_index % self.nvecs

        for key in self.ref_keys:
            getattr(self, key)[index] = arrs[key].copy()
        self.cur_nvecs = min(self.cur_nvecs + 1, self.nvecs)
        self.cur_index += 1
        # print(f"\tstored at {index=}, {self.cur_nvecs=}, {self.cur_index=}")
        return index

    def reset(self):
        """Reset DIIS container."""
        self.cur_index = 0
        self.cur_nvecs = 0
        # Also zero the arrays?!

    def get_coeffs(self) -> Optional[np.ndarray]:
        """Try to calculate DIIS coefficients from the error matrix."""
        if not self.can_get_coeffs():
            return None

        err_vecs = self.err_vecs[: self.cur_nvecs]
        err_mat = np.einsum("ij,kj->ik", err_vecs, err_vecs)
        try:
            coeffs = np.linalg.solve(err_mat, np.ones(self.cur_nvecs))
            coeffs /= coeffs.sum()
        except np.linalg.LinAlgError:
            coeffs = None
            self.reset()
        return coeffs

    def get(self, key: str) -> np.ndarray:
        return getattr(self, key)[: self.cur_nvecs]

    def can_get_coeffs(self) -> bool:
        return self.cur_nvecs >= 2

    def __str__(self):
        return f"DIIS(nvecs={self.nvecs}, {self.cur_nvecs} stored)"

    def __repr__(self):
        return self.__str__()
