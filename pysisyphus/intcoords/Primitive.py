import abc
import logging
from math import exp

import numpy as np

from pysisyphus.elem_data import COVALENT_RADII as CR
from pysisyphus.helpers_pure import hash_arr


class Primitive(metaclass=abc.ABCMeta):
    def __init__(self, indices, periodic=False, calc_kwargs=None, cache=False):
        self.indices = list(indices)
        self.periodic = periodic
        if calc_kwargs is None:
            calc_kwargs = ()
        self.calc_kwargs = calc_kwargs
        self.cache = cache

        self.logger = logging.getLogger("internal_coords")
        self.val_cache = {}
        self.grad_cache = {}

    def log(self, msg, lvl=logging.DEBUG):
        self.logger.log(lvl, msg)

    def log_dbg(self, msg):
        self.log(msg, lvl=logging.DEBUG)

    @staticmethod
    def parallel(u, v, thresh=1e-6):
        dot = u.dot(v) / (np.linalg.norm(u) * np.linalg.norm(v))
        return (1 - abs(dot)) < thresh

    @staticmethod
    def _get_cross_vec(coords3d, indices):
        m, o, n = indices
        # Select initial vector for cross product, similar to
        # geomeTRIC. It must NOT be parallel to u and/or v.
        x_dash = coords3d[n] - coords3d[m]
        x = x_dash / np.linalg.norm(x_dash)
        cross_vecs = np.eye(3)
        min_ind = np.argmin([np.dot(cv, x) ** 2 for cv in cross_vecs])
        return cross_vecs[min_ind]

    def set_cross_vec(self, coords3d, indices):
        self.cross_vec = self._get_cross_vec(coords3d, self.indices)
        self.log(f"Cross vector for {self} set to {self.cross_vec}")

    @abc.abstractmethod
    def _calculate(*, coords3d, indices, gradient, **kwargs):
        pass

    @abc.abstractmethod
    def _weight(self, atoms, coords3d, indices, f_damping):
        pass

    def weight(self, atoms, coords3d, f_damping=0.12):
        return self._weight(atoms, coords3d, self.indices, f_damping)

    @staticmethod
    def rho(atoms, coords3d, indices):
        i, j = indices
        distance = np.linalg.norm(coords3d[i] - coords3d[j])
        cov_rad_sum = CR[atoms[i].lower()] + CR[atoms[j].lower()]
        return exp(-(distance / cov_rad_sum - 1))

    # def calculate(self, coords3d, indices=None, gradient=False):
        # if indices is None:
            # indices = self.indices

        # # Gather calc_kwargs
        # calc_kwargs = {key: getattr(self, key) for key in self.calc_kwargs}

        # return self._calculate(
            # coords3d=coords3d,
            # indices=indices,
            # gradient=gradient,
            # **calc_kwargs,
        # )

    def calculate(self, coords3d, indices=None, gradient=False):
        if indices is None:
            indices = self.indices

        if self.cache:
            cur_hash = hash_arr(coords3d[indices], precision=8)
            try:
                val = self.val_cache[cur_hash]
                if gradient:
                    grad = self.grad_cache[cur_hash]
                    self.log_dbg(f"Returning cached value & gradient for hash '{cur_hash}'.")
                    return val, grad
                else:
                    self.log_dbg(f"Returning cached value for hash '{cur_hash}'.")
                    return val
            except KeyError:
                self.log_dbg(f"Hash '{cur_hash}' is not yet cached.")

        # Gather calc_kwargs
        calc_kwargs = {key: getattr(self, key) for key in self.calc_kwargs}

        results =  self._calculate(
            coords3d=coords3d,
            indices=indices,
            gradient=gradient,
            **calc_kwargs,
        )

        if self.cache:
            if gradient:
                val, grad = results
                self.val_cache[cur_hash] = val
                self.grad_cache[cur_hash] = grad
            else:
                self.val_cache[cur_hash] = results
        return results

    def jacobian(self, coords3d, indices=None):
        if indices is None:
            indices = self.indices

        # Gather calc_kwargs
        calc_kwargs = {key: getattr(self, key) for key in self.calc_kwargs}

        return self._jacobian(
            coords3d=coords3d,
            indices=indices,
            **calc_kwargs,
        )

    def __str__(self):
        return f"{self.__class__.__name__}({self.indices})"

    def __repr__(self):
        return self.__str__()
