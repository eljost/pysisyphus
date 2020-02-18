import numpy as np

from pysisyphus.calculators.Calculator import Calculator
from pysisyphus.linalg import perp_comp, make_unit_vec
from pysisyphus.optimizers.closures import small_lbfgs_closure
from pysisyphus.helpers import rms


class RotationConverged(Exception):
    pass


class DimerMethod(Calculator):

    rotation_types = ("fourier", "direct")

    def __init__(self, calculator, geometry, N_init, length=0.01, max_rotations=10,
                 rotation_type="direct", rotation_thresh=1e-3, rotation_tol=1.0,
                 interpolate=True):
        super().__init__(self)

        self.calculator = calculator
        self.geometry = geometry
        N_init = np.array(N_init, dtype=float)
        self.length = float(length)

        # Rotation parameters
        self.max_rotations = int(max_rotations)
        assert rotation_type in self.rotation_types, \
            f"Invalid rotation_type={rotation_type}! Valid types are: " \
            f"{self.rotation_types}"
        self.rotation_type = rotation_type
        self.rotation_thresh = float(rotation_thresh)
        self.rotation_tol = np.deg2rad(rotation_tol)
        self.interpolate = bool(interpolate)

        # Set normalized dimer direction
        if N_init is None:
            N_init = np.random.rand(self.geometry.coords.size)
        N_init /= np.linalg.norm(N_init)

        self.N = N_init.copy()
        self.atoms = self.geometry.atoms

        self._f0 = None
        self._energy0 = None
        self._f1 = None

    @property
    def coords0(self):
        return self.geometry.coords

    @coords0.setter
    def coords0(self, coords0_new):
        self.geometry.coords = coords0_new
        self._energy0 = None
        self._f0 = None
        self._f1 = None

    @property
    def coords1(self):
        return self.geometry.coords + self.length * self.N

    @coords1.setter
    def coords1(self, coords1_new):
        N_new = coords1_new - self.coords0
        N_new /= np.linalg.norm(N_new)
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
            self._f0 = results["forces"]
            self._energy0 = results["energy"]
        return self._f0

    @property
    def f1(self):
        if self._f1 is None:
            results = self.calculator.get_forces(self.atoms, self.coords1)
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
        return (f2 - f1).dot(N) / (2 * self.length)

    @property
    def C(self):
        """Curvature"""
        return self.curvature(self.f1, self.f2, self.N)

    def rotate(self, rad, theta):
        """Return new coords1"""
        step = (self.N*np.cos(rad) + theta*np.sin(rad)) * self.length
        return self.coords0 + step

    def direct_rotation_step(self):
        step = self.rot_force
        norm = np.linalg.norm(step)
        print(f"norm(rot_force)={norm:.6f}")
        trust = 0.0005
        if norm > trust:
            step = trust * step / norm

        new_coords1 = self.coords1 + step
        self.coords1  = new_coords1

    def fourier_rotation_step(self, optimizer):
        theta_dir = optimizer(self.rot_force)
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
            C_min_new = get_C(rad_min)
            # logger.debug( "Predicted theta_min lead us to a curvature maximum "
                         # f"(C(theta)={C_min:.6f}). Adding pi/2 to theta_min. "
                         # f"(C(theta+pi/2)={C_min_new:.6f})"
            # )
            C_min = C_min_new

        # TODO: handle cases where the curvature is still positive, but
        # the angle is small, so the rotation is skipped.
        # Don't do rotation for small angles
        if np.abs(rad_min) < self.rotation_tol:
            # logger.debug(f"rad_min={rad_min:.2f} below threshold. Breaking.")
            raise RotationConverged

        f1 = None
        # Interpolate force at coords1_rot; see Eq. (12) in [4]
        if self.interpolate:
            f1 = (np.sin(rad_trial - rad_min) / np.sin(rad_trial) * self.f1
                   + np.sin(rad_min) / np.sin(rad_trial) * f1_trial
                   + (1 - np.cos(rad_min) - np.sin(rad_min)
                      * np.tan(rad_trial / 2)) * self.f0
            )
        self.coords1 = self.rotate(rad_min, theta)
        self.f1 = f1

    # def get_energy(self, atoms, coords):
        # return self.calculator.get_energy(atoms, coords)

    def get_forces(self, atoms, coords):
        self.coords0 = coords
        lbfgs = small_lbfgs_closure()

        try:
            for i in range(self.max_rotations):
                self.fourier_rotation_step(lbfgs)
        except RotationConverged:
            print(f"Rotation converged in {i+1} cycles.")

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
        print("N", self.N)
        return results
