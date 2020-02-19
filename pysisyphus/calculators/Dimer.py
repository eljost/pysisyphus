import logging

import numpy as np

from pysisyphus.calculators.Calculator import Calculator
from pysisyphus.linalg import perp_comp, make_unit_vec
from pysisyphus.optimizers.closures import small_lbfgs_closure
from pysisyphus.optimizers.restrict_step import get_scale_max
from pysisyphus.helpers import rms


class RotationConverged(Exception):
    pass


class Dimer(Calculator):

    def __init__(self, calculator, *args, N_init=None, length=0.0189, rotation_max_cycles=15,
                 rotation_method="fourier", rotation_thresh=1e-4, rotation_tol=1,
                 rotation_max_element=0.001, rotation_interpolate=True, **kwargs):
        super().__init__(*args, **kwargs)

        self.logger = logging.getLogger("dimer")

        self.calculator = calculator
        self.length = float(length)

        # Rotation parameters
        self.rotation_max_cycles = int(rotation_max_cycles)

        rotation_methods = {
            "direct": self.direct_rotation,
            "fourier": self.fourier_rotation,
        }
        try:
            self.rotation_method = rotation_methods[rotation_method]
        except KeyError as err:
            print(f"Invalid rotation_method={rotation_method}! Valid types are: "
                  f"{tuple(self.rotation_methods.keys())}"
            )
            raise err
        self.rotation_thresh = float(rotation_thresh)
        self.rotation_tol = np.deg2rad(rotation_tol)
        self.rotation_max_element = float(rotation_max_element)
        self.rotation_interpolate = bool(rotation_interpolate)

        restrict_steps = {
            "direct": get_scale_max(self.rotation_max_element),
            "fourier": None,
        }
        self.restrict_step = restrict_steps[rotation_method]

        self._N = None
        self._coords0 = None
        self._energy0 = None
        self._f0 = None
        self._f1 = None

        self.force_evals = 0

        # Set dimer direction if given
        if N_init is not None:
            self.log("Setting initial orientation from given 'N_init'.")
            self.N = N_init

    def log(self, message=""):
        self.logger.debug(message)

    @property
    def N(self):
        return self._N

    @N.setter
    def N(self, N_new):
        N_new = np.array(N_new, dtype=float) / np.linalg.norm(N_new)
        self._N = N_new

    @property
    def coords0(self):
        return self._coords0

    @coords0.setter
    def coords0(self, coords0_new):
        self._coords0 = coords0_new
        self._energy0 = None
        self._f0 = None
        self._f1 = None

    @property
    def coords1(self):
        return self.coords0 + self.length * self.N

    @coords1.setter
    def coords1(self, coords1_new):
        N_new = coords1_new - self.coords0
        self.N = N_new
        self._f1 = None

    @property
    def energy0(self):
        if self._energy0 is None:
            results = self.calculator.get_energy(self.atoms, self.coords0)["energy"]
            self._energy0 = results["energy"]
        return self._energy0

    @property
    def f0(self):
        if self._f0 is None:
            results = self.calculator.get_forces(self.atoms, self.coords0)
            self.force_evals += 1
            self._f0 = results["forces"]
            self._energy0 = results["energy"]
        return self._f0

    @property
    def f1(self):
        if self._f1 is None:
            results = self.calculator.get_forces(self.atoms, self.coords1)
            self.force_evals += 1
            self._f1 = results["forces"]
        return self._f1

    @f1.setter
    def f1(self, f1_new):
        self._f1 = f1_new

    @property
    def f2(self):
        """Never calculated explicitly, but estimated from f0 and f1."""
        return 2 * self.f0 - self.f1

    @property
    def rot_force(self):
        f1_perp = perp_comp(self.f1, self.N)
        f2_perp = perp_comp(self.f2, self.N)

        f_perp = f1_perp - f2_perp
        return f_perp

    def curvature(self, f1, f2, N):
        """Curvature of the mode represented by the dimer."""
        return (f2 - f1).dot(N) / (2 * self.length)

    @property
    def C(self):
        """Shortcut for the curvature."""
        return self.curvature(self.f1, self.f2, self.N)

    def rotate(self, rad, theta):
        """Rotate dimer and produce new coords1."""
        return self.coords0 + (self.N*np.cos(rad) + theta*np.sin(rad)) * self.length

    def direct_rotation(self, optimizer, prev_step):
        rot_step = optimizer(self.rot_force, prev_step)
        rot_step = self.restrict_step(rot_step)
        # Strictly speaking rot_step should be constrained to conserve the desired
        # dimer length (coords1 - coords0)*2. This step is unconstrained.
        # Later on we calculate the actual step between the old coords1 and the new
        # coords1 that have been reconstrained.
        self.coords1 = self.coords1 + rot_step

    def fourier_rotation(self, optimizer, prev_step):
        theta_dir = optimizer(self.rot_force, prev_step)
        # Remove component that is parallel to N
        theta_dir = theta_dir - theta_dir.dot(self.N)*self.N
        theta = theta_dir / np.linalg.norm(theta_dir)

        # Get rotated endpoint geometries. The rotation takes place in a plane
        # spanned by N and theta. Theta is a unit vector perpendicular to N that
        # can be formed from the perpendicular components of the forces at the
        # endpoints.

        C = self.C
        # Derivative of the curvature, Eq. (29) in [2]
        # (f2 - f1) or -(f1 - f2)
        dC = 2*(self.f0 - self.f1).dot(theta) / self.length
        rad_trial = -0.5 * np.arctan2(dC, 2*abs(C))
        # logger.debug(f"rad_trial={rad_trial:.2f}")
        if np.abs(rad_trial) < self.rotation_tol:
            # logger.debug(f"rad_trial={rad_trial:.2f} below threshold. Breaking.")
            raise RotationConverged

        # Trial rotation for finite difference calculation of rotational force
        # and rotational curvature.
        coords1_trial = self.rotate(rad_trial, theta)
        f1_trial = self.calculator.get_forces(self.atoms, coords1_trial)["forces"]
        self.force_evals += 1
        f2_trial = 2*self.f0 - f1_trial
        N_trial = make_unit_vec(coords1_trial, self.coords0)
        C_trial = self.curvature(f1_trial, f2_trial, N_trial)

        b1 = 0.5 * dC
        a1 = (C - C_trial + b1 * np.sin(2 * rad_trial)) / (1 - np.cos(2 *  rad_trial))
        a0 = 2 * (C - a1)

        rad_min = 0.5 * np.arctan(b1/a1)
        # logger.debug(f"rad_min={rad_min:.2f}")
        def get_C(theta_rad):
            return a0/2 + a1*np.cos(2*theta_rad) + b1*np.sin(2*theta_rad)
        C_min = get_C(rad_min)  # lgtm [py/multiple-definition]
        if C_min > C:
            rad_min += np.deg2rad(90)
            # C_min_new = get_C(rad_min)
            # logger.debug( "Predicted theta_min lead us to a curvature maximum "
                         # f"(C(theta)={C_min:.6f}). Adding pi/2 to theta_min. "
                         # f"(C(theta+pi/2)={C_min_new:.6f})"
            # )

        # TODO: handle cases where the curvature is still positive, but
        # the angle is small, so the rotation is skipped.
        # Don't do rotation for small angles
        if np.abs(rad_min) < self.rotation_tol:
            # logger.debug(f"rad_min={rad_min:.2f} below threshold. Breaking.")
            raise RotationConverged

        f1 = None
        # Interpolate force at coords1_rot; see Eq. (12) in [4]
        if self.rotation_interpolate:
            f1 = (np.sin(rad_trial - rad_min) / np.sin(rad_trial) * self.f1
                   + np.sin(rad_min) / np.sin(rad_trial) * f1_trial
                   + (1 - np.cos(rad_min) - np.sin(rad_min)
                      * np.tan(rad_trial / 2)) * self.f0
            )

        self.coords1 = self.rotate(rad_min, theta)
        self.f1 = f1

    def get_forces(self, atoms, coords):
        # Generate random guess for the dimer orientation if not yet set
        if self.N is None:
            self.log("No initial orientation given. Using random guess.")
            self.N = np.random.rand(coords.size)
        self.atoms = atoms
        self.coords0 = coords

        lbfgs = small_lbfgs_closure()
        try:
            self.log("Starting dimer rotation")
            prev_step = None
            for i in range(self.rotation_max_cycles):
                rot_force = self.rot_force
                rms_rot_force = rms(rot_force)
                self.log(f"\t{i:02d}: rms(rot_force)={rms_rot_force:.6f}")
                if rms_rot_force <= self.rotation_thresh:
                    self.log("\trms(rot_force) is below threshold!")
                    raise RotationConverged
                coords1_old = self.coords1
                self.rotation_method(lbfgs, prev_step)
                actual_step = self.coords1 - coords1_old
                prev_step = actual_step
            else:
                msg =  "\tDimer rotation dit not converge in " \
                      f"{self.rotation_max_cycles}"
        except RotationConverged:
            msg = f"\tDimer rotation converged in {i+1} cycle(s)."
        self.log(msg)
        self.log("N after rotation:\n\t" + str(self.N))
        self.log()

        energy = self.energy0

        f0 = self.f0
        N = self.N

        f_parallel = f0.dot(N)*N
        f_perp = f0 - f_parallel
        f_tran = f_perp - f_parallel
        results = {
            "energy": energy,
            "forces": f_tran
        }
        return results
