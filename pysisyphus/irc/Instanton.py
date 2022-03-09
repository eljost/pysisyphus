# [1] https://doi.org/10.1021/acs.jctc.8b00068
#     Dual-Level Approach to Instanton Theory
#     Meisner, Kästner 2018
# [2] https://doi.org/10.1021/ct100658y
#     Locating Instantons in Many Degrees of Freedom
#     Rommel, Goumans, Kästner, 2011
# [3] https://aip.scitation.org/doi/full/10.1063/1.4932362
#     Ring-polymer instanton theory of electron transofer in the
#     nonadiabatic limit
#     Richardson, 2015
import logging

import numpy as np
import scipy as sp

from pysisyphus import logger as pysis_logger
from pysisyphus.Geometry import Geometry
from pysisyphus.helpers import align_coords, get_coords_diffs
from pysisyphus.helpers_pure import eigval_to_wavenumber


def T_crossover_from_eigval(eigval):
    nu = eigval_to_wavenumber(eigval)  # in cm⁻¹
    nu_m = 100 * abs(nu)  # in m⁻¹
    freq = nu_m * sp.constants.speed_of_light
    # In papers the crossover temperature is usually defined with the angular
    # frequency nu*2*pi. When using the "normal" frequency we don't have to
    # divide by 2*pi. Both approaches obviously yield the same T_c.
    # ang_freq = freq * 2 * np.pi
    # T_c_ = sp.constants.hbar * ang_freq / (sp.constants.Boltzmann * 2 * np.pi)
    T_c = sp.constants.hbar * freq / (sp.constants.Boltzmann)
    return T_c


def T_crossover_from_ts(ts_geom):
    mw_hessian = ts_geom.mw_hessian
    proj_hessian = ts_geom.eckart_projection(mw_hessian, full=True)
    eigvals, eigvecs = np.linalg.eigh(proj_hessian)
    T_c = T_crossover_from_eigval(eigvals[0])
    return T_c


logger = pysis_logger.getChild("instanton")
logger.setLevel(logging.DEBUG)
file_handler = logging.FileHandler("instanton.log", mode="w", delay=True)
logger.addHandler(file_handler)


def log_progress(val, key, i):
    logger.debug(f"Calculated {key} at image {i}")
    return val


class Instanton:
    def __init__(self, images, calc_getter, T):
        self.images = images
        self.calc_getter = calc_getter
        for image in self.images:
            image.set_calculator(calc_getter())
        self.T = T

        # Pre-calculate action prefactors for the given temperature
        beta = 1 / (sp.constants.Boltzmann * self.T)  # Joule
        beta_hbar = beta * sp.constants.hbar  # seconds
        beta_hbar_fs = beta_hbar * 1e15  # fs
        self.beta_hbar = beta_hbar_fs
        self.P_over_beta_hbar = self.P / self.beta_hbar
        self.P_bh = self.P_over_beta_hbar
        """The Instanton is periodic: the first image is connected to the
        last image.
        At k=0 the index k-1 will be -1, which points to the last image.

        Below we pre-calculate some indices (assuming N images).
            unshifted indices ks: k = {0, 1, .. , N-1}
            shifted indices ksm1: k-1 = {-1, 0, 1, .. , N-2}
            shifted indices ksp1: k+1 = {1, 2, .. , N-1, 0}
        """
        self.ks = np.arange(self.P)
        self.ksm1 = self.ks - 1
        self.ksp1 = self.ks + 1
        self.ksp1[-1] = 0  # k+1 for the last image points to the first image

        self.coord_type = "mwcartesian"
        self.internal = None

    @property
    def P(self):
        return len(self.images)

    @classmethod
    def from_ts(cls, ts_geom, P, dr=0.4, delta_T=25, cart_hessian=None, **kwargs):
        assert ts_geom.coord_type == "mwcartesian"
        atoms = ts_geom.atoms
        ts_coords = ts_geom.coords

        if cart_hessian is None:
            mw_hessian = ts_geom.mw_hessian
        else:
            mw_hessian = ts_geom.mass_weigh_hessian(cart_hessian)

        proj_hessian = ts_geom.eckart_projection(mw_hessian, full=True)
        eigvals, eigvecs = np.linalg.eigh(proj_hessian)
        # Use crossover temperature with a little offset (delta_T) if no T is given.
        try:
            kwargs["T"]
        except KeyError:
            T_c = T_crossover_from_eigval(eigvals[0])
            kwargs["T"] = T_c - delta_T
        imag_mode = eigvecs[:, 0]
        cosines = np.cos((np.arange(P) + 1 - 0.5) / P * np.pi)
        image_mw_coords = ts_coords + dr * cosines[:, None] * imag_mode
        image_coords = image_mw_coords / np.sqrt(ts_geom.masses_rep)
        if not ts_geom.is_analytical_2d:
            align_coords(image_coords)
        images = [
            Geometry(atoms, coords, coord_type="mwcartesian") for coords in image_coords
        ]
        instanton = Instanton(images, **kwargs)
        return instanton

    @classmethod
    def from_instanton(cls, other, **kwargs):
        images = other.images
        instanton = Instanton(images, **kwargs)
        return instanton

    def as_xyz(self):
        return "\n".join([geom.as_xyz() for geom in self.images])

    @property
    def coords(self):
        return np.ravel([image.coords for image in self.images])

    @coords.setter
    def coords(self, coords):
        coords = coords.reshape(len(self.images), -1)
        for img_coords, image in zip(coords, self.images):
            image.coords = img_coords

    def action(self):
        """Action in au / fs, Hartree per femtosecond."""
        all_coords = np.array([image.coords for image in self.images])
        diffs = all_coords[self.ks] - all_coords[self.ksm1]
        dists = np.linalg.norm(diffs, axis=1)
        S_0 = self.P_bh * (dists**2).sum()

        energies = [image.energy for image in self.images]
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

            d     sum_k (||y_k - y_(k-1)||_2)²
            ---
            d y_k

        The derivative of a norm of a vector difference is quite simple, but
        care has to be taken to recognize, that y_k appears two times in the sum.
        It appears in the first summand for k and in the second summand for k+1.

            sum_k (||y_k - y_(k-1)||_2)²
                        1. term                 2. term
                = (||y_k - y_(k-1)||_2)² + (||y_(k+1) - y_k||_2)² + ... and so on

        The derivative of the first term is

            2 * (y_k - y_(k-1))

        and the derivative of the second term is

            -2 * (y_(k+1) - y_k))

        which is equal to

            2 * (y_k - y_(k+1)) .

        To summarize:

            d     sum_k(||y_k - y_(k-1)||_2)²
            ---
            d y_k

            =   2 * (2 * y_k - y_(k-1) - y_(k+1)) .
        """
        image_coords = np.array([image.coords for image in self.images])
        kin_grad = (
            2
            * (
                2 * image_coords  # y_k
                - image_coords[self.ksm1]  # y_(k-1)
                - image_coords[self.ksp1]  # y_(k+1)
            ).flatten()
        )
        kin_grad *= self.P_bh
        pot_grad = np.array(
            [
                log_progress(image.gradient, "gradient", i)
                for i, image in enumerate(self.images)
            ]
        ).flatten()
        pot_grad /= self.P_bh
        gradient = kin_grad / 2 + pot_grad
        # gradient = pot_grad
        # gradient = kin_grad
        # gradient = kin_grad / 2
        results = {
            "gradient": gradient,
        }
        results.update(self.action())
        return results

    def action_hessian(self):
        image_hessians = [
            log_progress(image.hessian, "hessian", i)
            for i, image in enumerate(self.images)
        ]
        pot_hess = sp.linalg.block_diag(*image_hessians)
        pot_hess /= self.P_bh
        coord_num = pot_hess.shape[0]
        zeroes = np.zeros((coord_num, coord_num))
        image_coord_num = self.images[0].coords.size
        inds = np.arange(coord_num).reshape(-1, image_coord_num)
        ks = inds[self.ks].flatten()
        ksm1 = inds[self.ksm1].flatten()
        ksp1 = inds[self.ksp1].flatten()
        km = zeroes.copy()
        km[ks, ksm1] = 1.0
        kp = zeroes.copy()
        kp[ks, ksp1] = 1.0
        kin_hess = 4 * np.eye(coord_num) - 2 * km - 2 * kp  # y_k  # y_k-1  # y_k+1
        kin_hess *= self.P_bh
        hessian = kin_hess / 2 + pot_hess
        # hessian = pot_hess
        # hessian = kin_hess / 2
        results = {
            "hessian": hessian,
        }
        results.update(self.action())
        return results

    @property
    def energy(self):
        return self.action()["action"]

    @property
    def gradient(self):
        return self.action_gradient()["gradient"]

    @property
    def forces(self):
        return -self.gradient

    @property
    def hessian(self):
        return self.action_hessian()["hessian"]

    @property
    def cart_hessian(self):
        return sp.linalg.block_diag(*[image.cart_hessian for image in self.images])

    @property
    def cart_coords(self):
        return np.ravel([image.cart_coords for image in self.images])

    @property
    def cart_forces(self):
        return np.ravel([image.cart_forces for image in self.images])

    @property
    def cart_hessian(self):
        return sp.linalg.block_diag(*[image.cart_forces for image in self.images])

    def is_analytical_2d(self):
        return self.images[0].is_analytical_2d

    @property
    def path_length(self):
        image_coords3d = [image.coords3d for image in self.images]
        coord_diffs = get_coords_diffs(image_coords3d, align=False, normalize=False)
        length = coord_diffs.sum()
        return length

    def get_additional_print(self):
        length = self.path_length
        return f"\t\tInstanton length={length:.2f} √a̅m̅u̅·au"
