import logging
from pathlib import Path

import h5py
import numpy as np

from pysisyphus.calculators.Calculator import Calculator
from pysisyphus.helpers import rms, get_tangent_trj_str
from pysisyphus.intcoords.helpers import get_weighted_bond_mode
from pysisyphus.linalg import perp_comp, make_unit_vec
from pysisyphus.optimizers.closures import small_lbfgs_closure
from pysisyphus.optimizers.restrict_step import get_scale_max


class RotationConverged(Exception):
    pass


class Gaussian:
    def __init__(self, height, center, std, N):
        self.center = np.array(center)
        self.height = float(height)
        self.std = float(std)
        self.N = np.array(N)

        self.s0 = self.center.dot(self.N)

    def energy(self, R, height=None):
        if height is None:
            height = self.height

        return height * np.exp(-((R.dot(self.N) - self.s0) ** 2) / (2 * self.std ** 2))

    def forces(self, R, height=None):
        if height is None:
            height = self.height
        s_diff = R.dot(self.N) - self.s0

        return (
            height
            * np.exp(-(s_diff ** 2) / (2 * self.std ** 2))
            * s_diff
            / self.std ** 2
            * self.N
        )

    def __str__(self):
        return (
            f"Gaussian(height={self.height:.4f}, center={self.center}, "
            f"s0={self.s0:.4f}, std={self.std:.4f}"
        )


# [1] https://aip.scitation.org/doi/abs/10.1063/1.480097
#     Original Dimer paper
#     Henkelmann, 1999
# [2] https://doi.org/10.1063/1.1809574
#     Comparison of TS search methods
#     Olsen, 2004
# [3] https://doi.org/10.1063/1.2104507
#     Comparison of Dimer and P-RFO
#     Heyden, 2005
# [4] https://aip.scitation.org/doi/abs/10.1063/1.2815812
#     Superlinear Dimer method
#     Kästner, 2008


class Dimer(Calculator):
    def __init__(
        self,
        calculator,
        *args,
        N_raw=None,
        length=0.0189,
        rotation_max_cycles=15,
        rotation_method="fourier",
        rotation_thresh=1e-4,
        rotation_tol=1,
        rotation_max_element=0.001,
        rotation_interpolate=True,
        rotation_disable=False,
        rotation_disable_pos_curv=True,
        rotation_remove_trans=True,
        trans_force_f_perp=True,
        bonds=None,
        N_hessian=None,
        bias_rotation=False,
        bias_translation=False,
        bias_gaussian_dot=0.1,
        seed=None,
        write_orientations=True,
        forward_hessian=True,
        **kwargs,
    ):
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
            print(
                f"Invalid rotation_method={rotation_method}! Valid types are: "
                f"{tuple(self.rotation_methods.keys())}"
            )
            raise err
        self.rotation_thresh = float(rotation_thresh)
        self.rotation_tol = np.deg2rad(rotation_tol)
        self.rotation_max_element = float(rotation_max_element)
        self.rotation_interpolate = bool(rotation_interpolate)
        self.rotation_disable = bool(rotation_disable)
        self.rotation_disable_pos_curv = bool(rotation_disable_pos_curv)
        self.rotation_remove_trans = bool(rotation_remove_trans)
        self.trans_force_f_perp = trans_force_f_perp
        self.forward_hessian = forward_hessian

        # Regarding generation of initial orientation
        self.bonds = bonds
        self.N_hessian = N_hessian
        # Bias
        self.bias_rotation = bool(bias_rotation)
        self.bias_rotation_a = 0.0
        self.bias_translation = bool(bias_translation)
        self.bias_gaussian_dot = float(bias_gaussian_dot)

        self.write_orientations = write_orientations

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
        self.gaussians = list()

        # Set dimer direction if given
        self.N_raw = N_raw
        if self.N_raw is not None:
            msg = "Setting initial orientation from given 'N_raw'"
            if isinstance(self.N_raw, str) and Path(self.N_raw).exists():
                fn = self.N_raw
                self.N_raw = np.loadtxt(fn)
                msg = f"Read initial orientation from file '{fn}'"
                N_raw = self.N_raw.copy()
            self.N = N_raw
            self.log(msg)
        self.N_init = None

        if seed is None:
            # 2**32 - 1
            seed = np.random.randint(4294967295)
        np.random.seed(seed)
        msg = f"Using seed {seed} to initialize the random number generator.\n"
        print(msg)
        self.log(msg)

    @property
    def N(self):
        return self._N

    @N.setter
    def N(self, N_new):
        N_new = np.array(N_new, dtype=float).flatten()
        if self.rotation_remove_trans:
            N_new = self.remove_translation(N_new)
        N_new /= np.linalg.norm(N_new)
        self._N = N_new
        self._f1 = None

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
    def can_bias_f1(self):
        return (
            self.bias_rotation
            and (self.N_init is not None)
            and (self.bias_rotation_a is not None)
            and (self.bias_rotation_a > 0.0)
        )

    @property
    def should_bias_f1(self):
        """May lead to calculation of f0 and/or f1 if present!"""
        return self.can_bias_f1 and self.C > 0.0

    @property
    def can_bias_f0(self):
        return self.bias_translation and (self.N is not None)

    @property
    def should_bias_f0(self):
        """May lead to calculation of f0 and/or f1 if present!"""
        return self.can_bias_f0 and self.C > 0.0

    @property
    def f1_bias(self):
        # Apply bias force to f1 if desired. Dont apply bias if N_init is not (yet)
        # set. When N_raw was converged to a reasonable N_init we can add
        # the bias.
        assert self.bias_rotation_a >= 0.0, "This should not be negative!"

        fN = self.bias_rotation_a * self.length * self.N.dot(self.N_init) * self.N_init

        return fN

    @property
    def rot_force(self):
        f1_perp = perp_comp(self.f1, self.N)
        f2_perp = perp_comp(self.f2, self.N)
        f_perp = f1_perp - f2_perp

        # Don't bias rotational force if curvature is already negative
        if self.should_bias_f1:
            f1_bias_perp = perp_comp(self.f1_bias, self.N)
            f_perp += f1_bias_perp

        return f_perp

    def curvature(self, f1, f2, N):
        """Curvature of the mode represented by the dimer."""
        return (f2 - f1).dot(N) / (2 * self.length)

    @property
    def C(self):
        """Shortcut for the curvature."""
        return self.curvature(self.f1, self.f2, self.N)

    def get_gaussian_energies(self, coords, sum_=True):
        energies = [gauss.energy(coords) for gauss in self.gaussians]
        if sum_:
            energies = sum(energies)
        return energies

    def get_gaussian_forces(self, coords, sum_=True):
        forces = [gauss.forces(coords) for gauss in self.gaussians]
        if sum_:
            forces = np.sum(forces, axis=0)
        return forces

    def add_gaussian(
        self, atoms, center, N, height=0.1, std=0.0529, max_cycles=50, dot_ref=None
    ):
        # Create new gaussian object with default height that will be
        # refined later.
        gaussian = Gaussian(height=height, center=center, std=std, N=N)

        if dot_ref is None:
            dot_ref = self.bias_gaussian_dot

        # Calculate real forces at inflection point of new gaussian
        infl_coords = center + gaussian.std * N
        infl_results = self.calculator.get_forces(atoms, infl_coords)
        self.force_evals += 1
        infl_forces = infl_results["forces"]

        assert (
            infl_forces.dot(N) < 0
        ), "We probably overstepped the TS. See Section 2.3 in the paper."

        forces = self.get_gaussian_forces(infl_coords) + infl_forces

        def get_dot(height):
            """Dot product of forces for a given height and orientation N."""
            tmp_forces = forces.copy()
            tmp_forces += gaussian.forces(infl_coords, height=height)
            dot = tmp_forces.dot(N)
            return dot

        def can_break(dot):
            """Convergence indicator."""
            return abs(dot - dot_ref) <= 1e-3

        def bisect(
            min_,
            max_,
        ):
            for i in range(max_cycles):
                if abs(min_ - max_) <= 1e-10:
                    raise Exception("min_ and max_ became too similar!")

                # Determie value at half of the internval
                height = min_ + (max_ - min_) / 2
                dot = get_dot(height)

                if can_break(dot):
                    break

                if dot > dot_ref:
                    max_ = height
                elif dot < dot_ref:
                    min_ = height
            return height

        # Determine appropriate height
        grow = 2
        min_height = 0
        assert get_dot(0) < dot_ref
        for i in range(max_cycles):
            dot = get_dot(height)

            if can_break(dot):
                break

            if 0 < dot < dot_ref:
                min_height = height
                height *= grow
            elif dot < dot_ref:
                height *= grow
            else:
                height = bisect(min_height, height)
                break
        dot = get_dot(height)
        gaussian.height = height
        self.gaussians.append(gaussian)

        self.log(f"Added gaussian with height={height:.6f}")
        self.log(f"There are now {len(self.gaussians)} gaussians.")

        return gaussian

    def get_N_raw_from_hessian(self, h5_fn, root=0):
        with h5py.File(h5_fn, "r") as handle:
            hessian = handle["hessian"][:]
        w, v = np.linalg.eigh(hessian)
        assert w[root] < -1e-3
        N_raw = v[:, root]
        return N_raw

    def set_N_raw(self, coords):
        self.log("No initial orientation given. Generating one.")
        if self.bonds is not None:
            N_raw = get_weighted_bond_mode(self.bonds, coords.reshape(-1, 3))
            msg = "weighted bond mode"
        elif self.N_hessian is not None:
            N_raw = self.get_N_raw_from_hessian(self.N_hessian)
            msg = "first imaginary mode of HDF5 Hessian"
        else:
            msg = "random guess"
            N_raw = np.random.rand(coords.size)
        self.log(f"Obtained initial orientation from {msg}.")
        # Make N_raw translationally invariant and normalize
        self.N = N_raw
        # Now we keep the normalized dimer orientation
        self.N_raw = self.N

    def remove_translation(self, displacement):
        # Average vector components over cartesian directions (x,y,z)
        # Sum over each direction should equal zero if translationally invariant
        average = displacement.reshape(-1, 3).mean(axis=0)

        if max(abs(average)) > 1e-8:
            self.log(
                f"N-vector not translationally invariant. Removing average before normalization."
            )
        else:
            return displacement
        # Subtract the average component along each direction to make sum zero
        invariant_displacement = (
            displacement.reshape(-1, 3) - average[None, :]
        ).flatten()
        return invariant_displacement

    def rotate_coords1(self, rad, theta):
        """Rotate dimer and produce new coords1."""
        return self.coords0 + (self.N * np.cos(rad) + theta * np.sin(rad)) * self.length

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
        theta_dir = theta_dir - theta_dir.dot(self.N) * self.N
        theta = theta_dir / np.linalg.norm(theta_dir)

        # Get rotated endpoint geometries. The rotation takes place in a plane
        # spanned by N and theta. Theta is a unit vector perpendicular to N that
        # can be formed from the perpendicular components of the forces at the
        # endpoints.

        C = self.C
        # Derivative of the curvature, Eq. (29) in [2]
        # (f2 - f1) or -(f1 - f2)
        dC = 2 * (self.f0 - self.f1).dot(theta) / self.length
        rad_trial = -0.5 * np.arctan2(dC, 2 * abs(C))
        # self.log(f"rad_trial={rad_trial:.2f}")
        if np.abs(rad_trial) < self.rotation_tol:
            self.log(f"rad_trial={rad_trial:.2f} below threshold. Breaking.")
            raise RotationConverged

        # Trial rotation for finite difference calculation of rotational force
        # and rotational curvature.
        coords1_trial = self.rotate_coords1(rad_trial, theta)
        f1_trial = self.calculator.get_forces(self.atoms, coords1_trial)["forces"]
        self.force_evals += 1
        f2_trial = 2 * self.f0 - f1_trial
        N_trial = make_unit_vec(coords1_trial, self.coords0)
        C_trial = self.curvature(f1_trial, f2_trial, N_trial)

        b1 = 0.5 * dC
        a1 = (C - C_trial + b1 * np.sin(2 * rad_trial)) / (1 - np.cos(2 * rad_trial))
        a0 = 2 * (C - a1)

        rad_min = 0.5 * np.arctan(b1 / a1)
        # self.log(f"rad_min={rad_min:.2f}")
        def get_C(theta_rad):
            return a0 / 2 + a1 * np.cos(2 * theta_rad) + b1 * np.sin(2 * theta_rad)

        C_min = get_C(rad_min)  # lgtm [py/multiple-definition]
        if C_min > C:
            rad_min += np.deg2rad(90)
            # C_min_new = get_C(rad_min)
            # self.log("Predicted theta_min lead us to a curvature maximum "
            # f"(C(theta)={C_min:.6f}). Adding pi/2 to theta_min. "
            # f"(C(theta+pi/2)={C_min_new:.6f})"
            # )

        # TODO: handle cases where the curvature is still positive, but
        # the angle is small, so the rotation is skipped.
        # Don't do rotation for small angles
        if np.abs(rad_min) < self.rotation_tol:
            # self.log(f"rad_min={rad_min:.2f} below threshold. Breaking.")
            raise RotationConverged

        f1 = None
        # Interpolate force at coords1_rot; see Eq. (12) in [4]
        if self.rotation_interpolate:
            f1 = (
                np.sin(rad_trial - rad_min) / np.sin(rad_trial) * self.f1
                + np.sin(rad_min) / np.sin(rad_trial) * f1_trial
                + (1 - np.cos(rad_min) - np.sin(rad_min) * np.tan(rad_trial / 2))
                * self.f0
            )

        self.coords1 = self.rotate_coords1(rad_min, theta)
        self.f1 = f1

    def do_dimer_rotations(self, rotation_thresh=None):
        self.log("Doing dimer rotations")
        if rotation_thresh is None:
            rotation_thresh = self.rotation_thresh
            self.log(f"\tThreshold norm(rot_force)={rotation_thresh:.6f}")

        lbfgs = small_lbfgs_closure(gamma_mult=True)
        try:
            N_first = self.N
            prev_step = None
            for i in range(self.rotation_max_cycles):  # lgtm [py/redundant-else]
                N_cur = self.N
                rot_force = self.rot_force
                rms_rot_force = rms(rot_force)
                if self.should_bias_f1:
                    C_real = self.C
                    C_bias = -self.bias_rotation_a * (self.N.dot(self.N_init)) ** 2
                    C = C_real + C_bias
                    C_str = f"C={C: .6f}, C_real={C_real: .6f}, C_bias={C_bias: .6f}"
                else:
                    C_str = f"C={self.C: .6f}"
                self.log(f"\t{i:02d}: rms(rot_force)={rms_rot_force:.6f} {C_str}")
                if rms_rot_force <= rotation_thresh:
                    self.log("\trms(rot_force) is below threshold!")
                    raise RotationConverged
                coords1_old = self.coords1
                self.rotation_method(lbfgs, prev_step)
                actual_step = self.coords1 - coords1_old
                prev_step = actual_step
                rot_deg = np.rad2deg(np.arccos(N_cur.dot(self.N)))
                self.log(f"\t\tRotated by {rot_deg:.1f}°")
            else:
                msg = (
                    "\tDimer rotation did not converge in "
                    f"{self.rotation_max_cycles}"
                )
        except RotationConverged:
            msg = f"\tDimer rotation converged in {i+1} cycle(s)."
        self.log(msg)
        # self.log("\tN after rotation:\n\t" + str(self.N))
        self.log()
        # Restrict to interval [-1,1] where arccos is defined
        rot_deg = np.rad2deg(np.arccos(max(min(N_first.dot(self.N), 1.0), -1.0)))
        self.log(
            f"\tRotated by {rot_deg:.1f}° w.r.t. the orientation " "before rotation(s)."
        )

    def update_orientation(self, coords):
        # Generate random guess for the dimer orientation if not yet set
        if self.N is None:
            self.set_N_raw(coords)

        # Refine N_raw to N_init if not yet done
        if self.bias_rotation and self.N_init is None and self.C > 0.0:
            # Run initial sweep with a much softer convergence threshold
            self.log("Initial sweep to refine N_raw to N_init.")
            self.do_dimer_rotations(10 * self.rotation_thresh)
            if self.C > 0:
                self.N_init = self.N
                rot_rad = np.arccos(self.N_raw.dot(self.N_init))
                rot_deg = np.rad2deg(rot_rad)
                self.log(f"N_raw:\n\t{self.N_raw}")
                self.log(f"Rotated N_raw by {rot_deg:.1f}° to N_init")
                self.log(f"N_init:\n\t{self.N_init}")
                C = self.C
                self.log(f"Curvature after intial sweep is C={C:.6f}")
                self.log("Determining proper scaling factor for bias potential.")
                # Determine proper scaling factor for the quadratic bias potential
                # from the current curvature.
                scale_fact = 1.5
                self.bias_rotation_a = scale_fact * self.C
                self.log(f"Using a={scale_fact}*C={self.bias_rotation_a:.6f}")
            else:
                self.log(
                    f"Curvature after initial sweep C={self.C:.6f} is "
                    "already negative. Not setting N_init and bias_rotation_a!"
                )

        self.do_dimer_rotations()

    def get_forces(self, atoms, coords):
        self.atoms = atoms
        self.coords0 = coords

        try:
            N_backup = self.N.copy()
        except AttributeError:
            N_backup = None
        if not self.rotation_disable:
            self.update_orientation(coords)
        if (N_backup is not None) and self.rotation_disable_pos_curv and self.C > 0:
            self.log("Rotation did not yield a negative curvature. "
                     "Restoring previous unrotated N."
            )
            self.N = N_backup
        # Now we (have an updated self.N and) can do the force projections
        N = self.N
        # self.log(f"Orientation N:\n\t{N}")
        # Save orientation N in human-readable format, aka .trj
        # file in Angstrom.
        if self.write_orientations:
            trj_str = get_tangent_trj_str(atoms, coords, N)
            trj_fn = self.make_fn("N.trj")
            with open(trj_fn, "w") as handle:
                handle.write(trj_str)
            self.log(f"Wrote current orientation animation to '{trj_fn}'")
        # Always save orientation in Bohr
        N_fn = self.make_fn("N")
        np.savetxt(N_fn, N)

        energy = self.energy0
        self.log(f"\tenergy={self.energy0:.8f} au")

        f0 = self.f0

        if self.should_bias_f0:
            self.log("Biasing translation forces")
            self.log(f"There are currently {len(self.gaussians)} gaussians present.")
            bias_energy = self.get_gaussian_energies(coords)
            energy += bias_energy
            bias_forces = self.get_gaussian_forces(coords, sum_=False)
            try:
                bias_norms = np.linalg.norm(bias_forces, axis=1)
                bias_norms_str = np.array2string(bias_norms, precision=4)
                self.log(f"\tnorm(bias_forces)={bias_norms_str}")
            except np.AxisError:
                self.log("Skipping calculation of norm(bias_forces)")
            f0 += np.sum(bias_forces, axis=0)

        norm_f0 = np.linalg.norm(f0)
        self.log(f"\tnorm(forces)={norm_f0:.6f}")

        f_parallel = f0.dot(N) * N
        norm_parallel = np.linalg.norm(f_parallel)
        self.log(f"\tnorm(forces_parallel)={norm_parallel:.6f}")
        # self.log(f"\tforce_parallel:\n\t{f_parallel}")

        f_perp = f0 - f_parallel
        norm_perp = np.linalg.norm(f_perp)
        self.log(f"\tnorm(forces_perp)={norm_perp:.6f}")
        # self.log(f"\tforce_perp:\n\t{f_perp}")

        # Only return perpendicular component when curvature is negative
        f_tran = -f_parallel

        curv_str = "positive" if self.C > 0 else "negative"

        force_str = "reversed parallel component of"
        if (self.C < 0) or self.trans_force_f_perp:
            f_tran += f_perp
            force_str = "full"
        self.log(f"Curvature is {curv_str}. Returning {force_str} f_tran.")
        # self.log(f"\tf_tran:\n\t{f_tran}")
        self.log()

        self.calc_counter += 1

        results = {"energy": energy, "forces": f_tran}

        return results

    def get_hessian(self, atoms, coords):
        if not self.forward_hessian:
            raise Exception("Actual Hessian method not forwarded by Dimer calculator!")
        results = self.calculator.get_hessian(atoms, coords)
        return results
