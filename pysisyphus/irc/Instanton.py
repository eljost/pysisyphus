# [1] https://doi.org/10.1021/acs.jctc.8b00068
#     Dual-Level Approach to Instanton Theory
#     Meisner, Kästner 2018
# [2] https://doi.org/10.1021/ct100658y
#     Locating Instantons in Many Degrees of Freedom
#     Rommel, Goumans, Kästner, 2011
import numpy as np

from pysisyphus.config import T_DEFAULT
from pysisyphus.constants import AU2SEC, AU2J, HBAR, KBAU
from pysisyphus.Geometry import Geometry


class Instanton:
    def __init__(
        self, calc_getter, P, T=T_DEFAULT, ts_geom=None, other_instanton=None, dr=0.4
    ):
        self.calc_getter = calc_getter
        self.P = P
        self.T = T
        self.ts_geom = ts_geom
        self.other_instanton = other_instanton
        assert self.ts_geom or self.other_instanton
        self.dr = 0.4

        if self.ts_geom:
            assert self.ts_geom.coord_type == "mwcartesian"

        hbar_au = HBAR / AU2SEC / AU2J
        self.beta_hbar = 1 / (self.T * KBAU) * hbar_au
        # Instanton is periodic; so the first image is connected to the
        # last. We start at k=0, so k-1 will be -1, which points to the
        # last image.
        self.ks = np.arange(self.P)
        self.ksm1 = self.ks - 1
        self.ksp1 = self.ks + 1
        self.ksp1[-1] = 0
        self.P_over_beta_hbar = self.P / self.beta_hbar
        self.P_bh = self.P_over_beta_hbar

    def prepare_from_ts_geom(self):
        self.ts_atoms = self.ts_geom.atoms
        self.ts_coords = self.ts_geom.coords
        mw_hessian = self.ts_geom.mw_hessian
        if self.ts_geom.is_analytical_2d:
            proj_hessian = mw_hessian
        else:
            proj_hessian = self.ts_geom.eckart_projection(mw_hessian, full=True)
        eigvals, eigvecs = np.linalg.eigh(proj_hessian)
        imag_mode = eigvecs[:, 0]
        # TODO: calculate crossover Temperature here, if no T is given
        is_ = np.arange(self.P) + 1
        cos_ = np.cos((is_ - 0.5) / self.P * np.pi)
        all_mw_coords = self.ts_coords + self.dr * cos_[:, None] * imag_mode
        all_coords = all_mw_coords / np.sqrt(self.ts_geom.masses_rep)
        geoms = list()
        # self.ts_geom.jmol()
        for coords in all_coords:
            geom = Geometry(self.ts_atoms, coords, coord_type="mwcartesian")
            calc = self.calc_getter()
            geom.set_calculator(calc)
            # geom.jmol()
            geoms.append(geom)
        # TODO: align geometries?!
        return geoms

    def prepare(self):
        if self.ts_geom:
            geoms = self.prepare_from_ts_geom()
        self.images = geoms

    def as_xyz(self):
        return "\n".join([geom.as_xyz() for geom in self.images])

    def kin_grad(self, image_coords, inds1, inds2):
        """d ||image_coords[inds1]-image_coords[inds2]||_2 / d inds1"""
        diffs = image_coords[inds1] - image_coords[inds2]
        dists = np.linalg.norm(diffs, axis=1)
        grad = (diffs / dists[:, None]).flatten()
        return grad

    def action(self, energies=None):
        if energies is None:
            energies = [image.energy for image in self.images]

        all_coords = np.array([image.coords for image in self.images])
        diffs = all_coords[self.ks] - all_coords[self.ksm1]
        dists = np.linalg.norm(diffs, axis=1)
        S_0 = self.P_bh * dists.sum()

        S_pot = 1 / self.P_bh * sum(energies)
        S_E = S_0 / 2 + S_pot
        results = {
            "action": S_E,
        }
        return results

    def action_gradient(self):
        """
        kin_grad corresponds to the gradient of S_0 in (Eq. 1 in [1], or
        first term in Eq. 6 in [2].) It boils down to the derivative of a sum
        of vector norms

            d     sum_k||y_k - y_(k-1)||_2
            ---
            d y_k

        The derivative of a norm of a vector difference is quite simple, but
        care has to be taken to recognize, that y_k appears two times in the sum.
        It appears in the first summand for k and in the second summand for k+1.

                                           1. term                 2. term
            sum_k||y_k - y_(k-1)||_2 = ||y_k - y_(k-1)||_2 + ||y_(k+1) - y_k||_2
                                       + ... and so on

        The derivative of the first term is

            (y_k - y_(k-1)) / ||y_k - y_(k-1)||_2

        and the derivative of the second term is

            -(y_(k+1) - y_k)) / ||y_(k+1) - y_k)||_2

        which is equal to

            (y_k - y_(k+1)) / ||y_(k+1) - y_k)||_2 .

        By using the fact that ||a-b||_2 == ||b-a||_2 we can use self.kin_grad
        for both terms by calling it with ks and ks+1 for the 2nd term, as this
        switches the signs.
        """
        all_coords = np.array([image.coords for image in self.images])
        # Grad of 1. term
        kin_grad = self.kin_grad(all_coords, self.ks, self.ksm1)
        # Grad of 2. term
        kin_grad += self.kin_grad(all_coords, self.ks, self.ksp1)
        kin_grad *= self.P_bh
        pot_grad = np.array([image.gradient for image in self.images]).flatten()
        pot_grad /= self.P_bh
        gradient = kin_grad / 2 + pot_grad
        results = {
            "gradient": gradient,
        }
        energies = [image.energy for image in self.images]
        results.update(self.action(energies=energies))
        return results

    @property
    def coords(self):
        return np.array(([image.coords for image in self.images]))

    @coords.setter
    def coords(self, coords):
        coords = coords.reshape(len(self.images), -1)
        for img_coords, image in zip(coords, self.images):
            image.coords = img_coords
