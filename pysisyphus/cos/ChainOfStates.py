# [1]   http://dx.doi.org/10.1063/1.1323224
#       Improved tangent estimate in the nudged elastic band method
#       for finding minimum energy paths and saddle points
#       Henkelman, Jonsson, 2000

from copy import copy
import logging
import sys
from typing import Literal, Optional

from distributed import Client, LocalCluster
import numpy as np
import psutil
from scipy.interpolate import interp1d, splprep, splev

from pysisyphus.helpers import align_coords, get_coords_diffs
from pysisyphus.helpers_pure import hash_arr, rms
from pysisyphus.modefollow import geom_lanczos

from pysisyphus.cos.distributed import distributed_calculations


class ClusterDummy:
    def close(self):
        pass


class ChainOfStates:
    logger = logging.getLogger("cos")
    valid_coord_types = ("cart", "cartesian", "dlc")

    def __init__(
        self,
        images,
        fix_first=True,
        fix_last=True,
        align_fixed=True,
        climb=False,
        climb_rms=5e-3,
        climb_lanczos=False,
        climb_lanczos_rms=5e-3,
        climb_fixed=False,
        ts_opt=False,
        ts_opt_rms=2.5e-3,
        energy_min_mix=False,
        scheduler=None,
        cluster=False,
        cluster_kwargs=None,
        progress=False,
    ):
        assert len(images) >= 2, "Need at least 2 images!"
        self.images = list(images)
        self.fix_first = fix_first
        self.fix_last = fix_last
        self.align_fixed = align_fixed
        self.climb = climb
        self.climb_rms = climb_rms
        self.climb_lanczos = climb_lanczos
        # Must not be bigger than climb_rms, so Lanczos tangent is NOT activated
        # before a climbing image is assigned.
        self.climb_lanczos_rms = min(self.climb_rms, climb_lanczos_rms)
        self.ts_opt = ts_opt
        self.ts_opt_rms = ts_opt_rms
        # I really wonder what made me pick 'climb_fixed = True' in e9f039a0 ...
        # ... now i know! Sometimes it hampers the convergence, when the CI index
        # flips between multiple images.
        self.climb_fixed = climb_fixed
        self.energy_min_mix = energy_min_mix
        self.scheduler = scheduler
        # Providing only cluster_kwargs also enables dask
        self._cluster = bool(cluster) or (cluster_kwargs is not None)
        if cluster_kwargs is None:
            cluster_kwargs = dict()
        # _cluster_kwargs = {
        # # Only one calculation / worker
        # "threads_per_worker": 1,
        # "n_workers": psutil.cpu_count(logical=False),
        # }
        ncores = psutil.cpu_count(logical=False)
        _cluster_kwargs = {
            # One worker per cluster with multiple threads
            "threads_per_worker": ncores,
            "n_workers": 1,
            # Register available cores as resource
            "resources": {
                "CPU": ncores,
            },
        }
        _cluster_kwargs.update(cluster_kwargs)
        self.cluster_kwargs = _cluster_kwargs
        self.progress = progress

        self._external_scheduler = scheduler is not None
        # Can't call self.log before eval counter are initialized ...
        self.image_energy_evals = 0
        self.image_force_evals = 0

        self._coords = None
        self._image_energies = None
        self._energy = None
        self._image_forces = None
        self._perpendicular_forces = None
        self._forces = None
        self._tangents = None

        self.coords_length = self.images[0].coords.size
        self.cart_coords_length = self.images[0].cart_coords.size
        self.zero_vec = np.zeros(self.coords_length)

        """
        It was a rather unfortunate choice to fill/grow self.coords_list and
        self.all_true_forces in different methods. self.prepare_opt_cycle()
        appends to self.coords_list, while self.calculate_forces() appends to
        self.all_true_forces().

        After a succsessful COS optimization both lists differ in length;
        self.all_true_forces has 1 additional item, compared to self.coords_list,
        as self.calculate_forces() is called in Optimizer.prepare_opt() once.
        Afterwards, self.coords_list and self.all_true_forces grow in a consistent
        manner.

        Two choices can be made: keep this discrepancy in mind and omit/neglect
        the first item in self.coords_list, or grow another list in
        self.calculate_forces(). For now, we will go with the latter option.
        """
        # coords_list & force_list are updated in prepare_opt_cycle
        self.coords_list = list()
        self.forces_list = list()
        self.perp_forces_list = list()
        # all_energies & all_true_forces are updated in calculate_forces()
        self.all_energies = list()
        self.all_true_forces = list()
        # See the multiline comment above
        self.all_cart_coords = list()
        self.lanczos_tangents = dict()
        self.prev_lanczos_hash = None

        # Start climbing immediateley with climb_rms == -1
        self.started_climbing = self.climb_rms == -1
        if self.started_climbing:
            self.log("Will start climbing immediately.")
        self.started_climbing_lanczos = self.climb_lanczos_rms == -1
        if self.started_climbing_lanczos:
            self.log("Will use Lanczos algorithm immediately.")
        self.started_ts_opt = self.ts_opt_rms == -1
        if self.started_ts_opt:
            self.log("Will use TS image(s) immediately.")
        self.fixed_climb_indices = None
        self.started_energy_min_mixing = False
        self.fixed_ts_indices = None
        # Use original forces for these images
        self.org_forces_indices = list()

        img0 = self.images[0]
        self.image_atoms = copy(img0.atoms)
        self.coord_type = img0.coord_type
        assert (
            self.coord_type in self.valid_coord_types
        ), f"Invalid coord_type! Supported types are: {self.valid_coord_types}"
        assert all(
            [img.coord_type == self.coord_type for img in self.images]
        ), "coord_type of images differ!"
        try:
            self.typed_prims = img0.internal.typed_prims
        except AttributeError:
            self.typed_prims = None

    def has_tangents(self):
        return self._tangents != None

    def has_image_forces(self):
        return self._image_forces != None

    @property
    def nimages(self) -> int:
        return len(self.images)

    @property
    def calculator(self):
        try:
            calc = self.images[0].calculator
        except IndexError:
            calc = None
        return calc

    def propagate(self, force_unique_overlap_data_fns=True):
        """Propagate chkfiles and root information along COS.

        Does an energy calculation at every image and tries to propagate
        the converged wavefunction to the next image.

        When excited states should be tracked it is assumed that the correct root is
        set at the first image. In most cases, e.g. in COS calculations started from
        YAML input the initial root will be the same for all images. When the correct
        root switches between the first image (with the correct root) and the last image,
        then using the same root for all images will be wrong. This is also corrected
        here. From the second image on the excited state (ES) overlaps are calculated
        between the current image and the image before it and the root with the highest
        overlap to the root on the previous image is set on the current image. As we
        start from the correct root at the first image this will ensure that the correct
        root is selected at all images.

        If requested, unique names for the dumped overlap data HDF5 will be picked.
        """
        try:
            track = self.images[0].calculator.track
        except AttributeError:
            track = False

        # To avoid cyclic import
        from pysisyphus.calculators.OverlapCalculator import (
            track_root_between_ovlp_cals,
        )

        # If requested, check if all calculators dump their overlap data into unique
        # HDF5 files. If they don't, update the filenames.
        # Won't really work for Growing-String calculations ...
        if track and force_unique_overlap_data_fns:
            nimages = len(self.images)
            cur_dump_fns = set([image.calculator.dump_fn for image in self.images])
            if len(cur_dump_fns) != nimages:
                for i, image in enumerate(self.images):
                    calc = image.calculator
                    new_dump_fn = calc.dump_fn.with_name(
                        f"overlap_data_{calc.calc_number:03d}.h5"
                    )
                    calc.dump_fn = new_dump_fn
                    self.log(f"Updated dump_fn of image {i} to {new_dump_fn.name}")
            assert (
                len(set([image.calculator.dump_fn for image in self.images])) == nimages
            )
        if track:
            for image in self.images:
                image.calculator.ovlp_with = "previous"

        # Run an energy calculation on all calculators sequentially,
        # then propagate chkfiles and root-information along them.
        for i, image in enumerate(self.images):
            # Do an energy calculation. Don't worry when the calculated/selected root
            # is actually the wrong one. This will be correct below and from the second
            # iteration on, the roots will be correct.
            image.energy
            self.log(f"Calculated energy for image {i}.")

            calc = image.calculator
            # Try to set chkfiles on the next calculator
            try:
                next_calc = self.images[i + 1].calculator
                next_calc.set_chkfiles(calc.get_chkfiles())
                self.log(f"Set chkfiles of image {i} at image {i+1}.")
            except IndexError:
                self.log("Last image. No further image to set chkfiles on!")
            except AttributeError:
                self.log("Calculator does not implement set_chkfiles/get_chkfiles!")

            # Update root information when it is a tracking calculation
            if track and i > 0:
                prev_calc = self.images[i - 1].calculator
                # prev_root = prev_calc.root
                # root = calc.root
                new_root = track_root_between_ovlp_cals(prev_calc, calc)
                calc.root = new_root
                self.log(
                    f"Calculated ES overlaps between image {i-1} and image {i}. "
                    f"New root at image {i} is {new_root}."
                )
                # self.log(
                # f"Root at image {i-1}: {prev_root}, previous root at image {i}: {root}, "
                # f"new root at image {i}: {new_root}."
                # )

        # Clear the lists storing infomration about the calculated roots.
        # If we don't reset them then the first reference cycle will be a cycle with a
        # potentially erroneous root, leading to undesired root flips.
        if track:
            for i, image in enumerate(self.images):
                image.calculator.clear_stored_calculations()
                image.clear()
                print(f"Using root {image.calculator.root} for image {i}.")
        self.log("Finished propagating chkfiles/roots along COS.")
        # TODO: return roots?!

    @property
    def is_analytical_2d(self):
        try:
            ia2d = self.images[0].calculator.analytical_2d
        except AttributeError:
            ia2d = False
        return ia2d

    def log(self, message):
        self.logger.debug(message)

    def get_fixed_indices(self):
        fixed = list()
        if self.fix_first:
            fixed.append(0)
        if self.fix_last:
            fixed.append(len(self.images) - 1)
        return fixed

    def get_dask_local_cluster(self):
        cluster = LocalCluster(**self.cluster_kwargs)
        link = cluster.dashboard_link
        self.log(f"Created {cluster}. Dashboard is available at {link}")
        return cluster

    @property
    def use_dask(self):
        return self.scheduler is not None

    def init_dask(self):
        # Return dummy when scheduler address is already present or dask is disabled.
        if self._external_scheduler or not self._cluster:
            cluster = ClusterDummy()
        # Create LocalCluster and return it if dask is enabled but no scheduler address
        # is set.
        # It seems like we can't save the LocalCluster object in the ChainOfStates object,
        # because this leads to strange errors while object pickling.
        elif self.scheduler is None:
            cluster = self.get_dask_local_cluster()
            self.scheduler = cluster.scheduler_address
        else:
            raise Exception("How did I end up here?")
        return cluster

    def exit_dask(self, cluster):
        # Don't do anything is cluster/scheduler is externally managed
        if not self._external_scheduler:
            cluster.close()
            self.scheduler = None
            self._external_scheduler = False

    @property
    def moving_indices(self):
        """Returns the indices of the images that aren't fixed and can be
        optimized."""
        fixed = self.get_fixed_indices()
        return [i for i in range(len(self.images)) if i not in fixed]

    @property
    def last_index(self):
        return len(self.images) - 1

    @property
    def moving_images(self):
        return [self.images[i] for i in self.moving_indices]

    @property
    def max_image_num(self):
        return len(self.images)

    @property
    def image_inds(self):
        return list(range(self.max_image_num))

    def zero_fixed_vector(self, vector):
        fixed = self.get_fixed_indices()
        for i in fixed:
            vector[i] = self.zero_vec
        return vector

    def clear(self):
        self._image_energies = None
        self._energy = None
        self._image_forces = None
        self._perpendicular_forces = None
        self._forces = None
        self._hessian = None
        self._tangents = None

    @property
    def atoms(self):
        atoms_ = self.images[0].atoms
        return len(self.images) * atoms_

    def set_vector(self, name, vector, clear=False):
        vec_per_image = vector.reshape(-1, self.coords_length)
        assert len(self.images) == len(vec_per_image)
        for i in self.moving_indices:
            setattr(self.images[i], name, vec_per_image[i])
        if clear:
            self.clear()

    @property
    def coords(self):
        """Return a flat 1d array containing the coordinates of all images."""
        all_coords = [image.coords for image in self.images]
        # Note: why does this getter set self._coords? ... I wrote this line 6 years ago.
        self._coords = np.concatenate(all_coords)
        return self._coords

    @coords.setter
    def coords(self, coords):
        """Distribute the flat 1d coords array over all images."""
        self.set_vector("coords", coords, clear=True)

    @property
    def cart_coords(self):
        """Return a flat 1d array containing the cartesian coordinates of all
        images."""
        return np.concatenate([image.cart_coords for image in self.images])

    @property
    def coords3d(self):
        assert self.images[0].coord_type == "cart"
        return self.coords.reshape(-1, 3)

    @property
    def image_coords(self):
        return np.array([image.coords for image in self.images])

    def set_coords_at(self, i, coords):
        """Called from helpers.procrustes with cartesian coordinates.
        Then tries to set cartesian coordinate as self.images[i].coords
        which will raise an error when coord_type != "cart".
        """
        assert self.images[i].coord_type in ("cart", "cartesian"), (
            "ChainOfStates.set_coords_at() has to be reworked to support "
            "internal coordiantes. Try to set 'align: False' in the 'opt' "
            "section of the .yaml input file."
        )
        if i in self.moving_indices:
            self.images[i].coords = coords
        # When dealing with a fixed image don't set coords through the
        # property, which would result in resetting the image's calculated
        # data. Instead, assign coords directly. This only occurs when
        # aligning the fixed images.
        elif self.align_fixed:
            self.images[i]._coords = coords

    @property
    def energy(self):
        """Currently, ChainOfStates.energy and ChainOfStates.image_energies are the same."""
        if self._energy is None:
            self._energy = self.image_energies
        return self._energy

    @energy.setter
    def energy(self, energies):
        """This is needed for some optimizers like CG and BFGS, when they skip a step
        and reset the Geometry."""
        assert len(self.images) == len(energies)
        # We do not allow setting the energy of a fixed image
        for i in self.moving_indices:
            self.images[i].energy = energies[i]

        self._energy = energies

    def calc_image_energy(self, image):
        image.calc_energy()
        return image

    def calc_image_energy_and_forces(self, image):
        image.calc_energy_and_forces()
        return image

    def concurrent_image_calcs(
        self, method, images_to_calculate, image_indices, calc_kind: str
    ):
        client = self.get_dask_client()
        self.log(f"Doing concurrent {calc_kind} calculations using {client}")
        calculated_images = distributed_calculations(
            client,
            images_to_calculate,
            method,
            logger=self.logger,
        )
        # Set calculated images
        for ind, image in zip(image_indices, calculated_images):
            self.images[ind] = image
        client.close()

    def concurrent_energy_calcs(self, images_to_calculate, image_indices):
        return self.concurrent_image_calcs(
            self.calc_image_energy,
            images_to_calculate,
            image_indices,
            "energy",
        )

    def concurrent_force_calcs(self, images_to_calculate, image_indices):
        return self.concurrent_image_calcs(
            self.calc_image_energy_and_forces,
            images_to_calculate,
            image_indices,
            "forces",
        )

    def get_images_to_calculate(self, image_indices=None):
        # Determine the number of images for which we have to do calculations.
        # There may also be calculations for fixed images, as they need an
        # energy value for an upwinding tangent calculation.
        # But in principle, every fixed image needs only an energy calculation once.
        if image_indices is None:
            images_to_calculate = self.moving_images
            image_indices = self.moving_indices
            if self.fix_first and (self.images[0]._energy is None):
                images_to_calculate = [self.images[0]] + images_to_calculate
                image_indices = [0] + list(image_indices)
            if self.fix_last and (self.images[-1]._energy is None):
                images_to_calculate = images_to_calculate + [self.images[-1]]
                image_indices = list(image_indices) + [-1]
            assert len(images_to_calculate) <= len(self.images)
        else:
            images_to_calculate = [self.images[i] for i in image_indices]
        return images_to_calculate, image_indices

    def calculate_image_energies(self, image_indices=None):
        images_to_calculate, image_indices = self.get_images_to_calculate(image_indices)

        # Parallel calculation with dask
        if self.use_dask:
            self.concurrent_energy_calcs(images_to_calculate, image_indices)
        # Serial calculation
        else:
            for image in images_to_calculate:
                image.calc_energy()
                # Poor mans progress bar ;)
                if self.progress:
                    print(".", end="")
                    sys.stdout.flush()
            if self.progress:
                print("\r", end="")
        self.image_energy_evals += 1

        energies = np.array([image.energy for image in self.images])
        return {
            "energies": energies,
        }

    @property
    def image_energies(self):
        if self._image_energies is None:
            image_results = self.calculate_image_energies()
            self._image_energies = image_results["energies"]
        return self._image_energies

    # We probably should not allow setting the energies of the underlying images.
    # @image_energies.setter
    # def image_energies(self, image_energies):
    # assert image_energies.size == self.nimages
    # # TODO: handle setting of fixed images
    # # At the least the calculation of fixed image energies seems to be already
    # # handled in calculate_image_energies/calculate_image_forces.
    # self._image_energies = image_energies

    def calculate_image_forces(self, image_indices=None):
        images_to_calculate, image_indices = self.get_images_to_calculate(image_indices)

        if self.use_dask:
            self.concurrent_force_calcs(images_to_calculate, image_indices)
        # Serial calculation
        else:
            for image in images_to_calculate:
                image.calc_energy_and_forces()
                # Poor mans progress bar ;)
                if self.progress:
                    print(".", end="")
                    sys.stdout.flush()
            if self.progress:
                print("\r", end="")
        self.set_zero_forces_for_fixed_images()
        self.image_force_evals += 1

        energies = np.array([image.energy for image in self.images])
        forces = np.array([image.forces for image in self.images])
        self.all_energies.append(energies)
        self.all_true_forces.append(forces.copy())
        self.all_cart_coords.append(self.cart_coords.copy())

        return {
            "energies": energies,
            "forces": forces,
        }

    @property
    def image_forces(self):
        if self._image_forces is None:
            image_results = self.calculate_image_forces()
            self._image_energies = image_results["energies"]
            self._image_forces = image_results["forces"]
        return self._image_forces

    # We probably should not allow setting the forces of the underlying images.
    # @image_forces.setter
    # def image_forces(self, image_forces):
    # # TODO: handle setting of fixed images
    # # At the least the calculation of fixed image energies seems to be already
    # # handled in calculate_image_energies/calculate_image_forces.
    # self._image_forces = image_forces

    def prepare_forces(self):
        image_forces = self.image_forces

        # Tangents are required for the projection
        tangents = self.tangents
        self.perpendicular_forces = self.calculate_perpendicular_forces(
            image_forces, tangents
        )
        self.perp_forces_list.append(self.perpendicular_forces.copy())
        forces = self.perpendicular_forces.copy()
        return forces

    def finalize_forces(self, forces):
        image_forces = self.image_forces
        image_energies = self.image_energies
        tangents = self.tangents
        self.update_with_climbing_forces(forces, tangents, image_energies, image_forces)
        self.update_with_org_forces(forces, self.org_forces_indices, image_forces)
        self.forces = forces.flatten()

    @property
    def forces(self):
        if self._forces is None:
            forces = self.prepare_forces()
            self.finalize_forces(forces)
        return self._forces

    @forces.setter
    def forces(self, forces):
        # TODO: zero forces for fixed images
        self._forces = forces

    def calculate_perpendicular_forces(self, image_forces, tangents):
        """[1] Eq. 12"""
        perp_forces = np.zeros((self.nimages, self.coords_length))
        for j, _ in enumerate(self.image_inds):
            if j not in self.moving_indices:
                continue
            fi = image_forces[j]
            ti = tangents[j]
            perp_forces[j] = fi - fi.dot(ti) * ti
        return perp_forces

    @property
    def perpendicular_forces(self):
        return self._perpendicular_forces

    @perpendicular_forces.setter
    def perpendicular_forces(self, perpendicular_forces):
        self._perpendicular_forces = perpendicular_forces

    @property
    def pure_perpendicular_forces(self):
        assert self.has_image_forces and self.has_tangents
        return self.perpendicular_forces

    @property
    def gradient(self):
        return -self.forces

    @gradient.setter
    def gradient(self, gradient):
        self.forces = -gradient

    @property
    def masses_rep(self):
        return np.array([image.masses_rep for image in self.images]).flatten()

    @property
    def results(self):
        tmp_results = list()
        for image in self.images:
            res = image.results
            res["coords"] = image.coords
            res["cart_coords"] = image.cart_coords
            tmp_results.append(res)
        return tmp_results

    def set_zero_forces_for_fixed_images(self):
        """This is always done in cartesian coordinates, independent
        of the actual coord_type of the images as setting forces only
        work with cartesian forces."""
        zero_forces = np.zeros_like(self.images[0].cart_coords)
        if self.fix_first:
            self.images[0].cart_forces = zero_forces
            self.log("Zeroed forces on fixed first image.")
        if self.fix_last:
            self.images[-1].cart_forces = zero_forces
            self.log("Zeroed forces on fixed last image.")

    def get_tangent(
        self,
        i: int,
        kind: Literal["upwinding", "bisect", "simple", "lanczos"] = "upwinding",
        lanczos_guess: Optional[np.ndarray] = None,
        disable_lanczos: bool = False,
        energies: Optional[np.ndarray] = None,
    ):
        """[1] Equations (8) - (11).

        To avoid any calculations when using upwinding tangents energies must be passed
        explicitly for this kind. The Lanczos tangent will lead to additional calculations
        nonetheless.
        """
        if kind == "upwinding":
            assert energies is not None and (
                len(energies) == self.nimages
            ), "'upwinding'-tangents require energies!"

        # Converge to lowest curvature mode at the climbing image.
        # In the current implementation the given kind may be overwritten when
        # Lanczos iterations are enabled and there are climbing images. By
        # setting 'disable_lanczos=True' the provided kind is never overwritten.
        if (
            not disable_lanczos
            and self.started_climbing_lanczos
            # and (i in self.get_climbing_indices())
            and (i == self.get_hei_index())
        ):
            kind = "lanczos"

        tangent_kinds = ("upwinding", "simple", "bisect", "lanczos")
        assert kind in tangent_kinds, "Invalid kind! Valid kinds are: {tangent_kinds}"
        prev_index = max(i - 1, 0)
        next_index = min(i + 1, len(self.images) - 1)

        prev_image = self.images[prev_index]
        ith_image = self.images[i]
        next_image = self.images[next_index]

        # If (i == 0) or (i == len(self.images)-1) then one
        # of this tangents is zero.
        tangent_plus = next_image - ith_image
        tangent_minus = ith_image - prev_image

        # Handle first and last image
        if i == 0:
            return tangent_plus / np.linalg.norm(tangent_plus)
        elif i == (len(self.images) - 1):
            return tangent_minus / np.linalg.norm(tangent_minus)

        # [1], Eq. (1)
        if kind == "simple":
            tangent = next_image - prev_image
        # [1], Eq. (2)
        elif kind == "bisect":
            first_term = tangent_minus / np.linalg.norm(tangent_minus)
            sec_term = tangent_plus / np.linalg.norm(tangent_plus)
            tangent = first_term + sec_term
        # Upwinding tangent from [1] Eq. (8) and so on
        elif kind == "upwinding":
            # If no energies are present at the energies accessing the energy
            # attribute will initiate an energy calculation. As
            prev_energy = energies[prev_index]
            ith_energy = energies[i]
            next_energy = energies[next_index]

            next_energy_diff = abs(next_energy - ith_energy)
            prev_energy_diff = abs(prev_energy - ith_energy)
            delta_energy_max = max(next_energy_diff, prev_energy_diff)
            delta_energy_min = min(next_energy_diff, prev_energy_diff)

            # Uphill
            if next_energy > ith_energy > prev_energy:
                tangent = tangent_plus
            # Downhill
            elif next_energy < ith_energy < prev_energy:
                tangent = tangent_minus
            # Minimum or Maximum
            else:
                if next_energy >= prev_energy:
                    tangent = (
                        tangent_plus * delta_energy_max
                        + tangent_minus * delta_energy_min
                    )
                # next_energy < prev_energy
                else:
                    tangent = (
                        tangent_plus * delta_energy_min
                        + tangent_minus * delta_energy_max
                    )
        elif kind == "lanczos":
            # Calculating a lanczos tangent is costly, so we store the
            # tangent in a dictionary. The current coordinates are
            # stringified with precision=4 and then hashed. The tangent
            # is stored/looked up with this hash.
            cur_hash = hash_arr(ith_image.coords, precision=4)
            try:
                tangent = self.lanczos_tangents[cur_hash]
                self.log(
                    "Returning previously calculated Lanczos tangent with "
                    f"hash={cur_hash}"
                )
            except KeyError:
                # Try to use previous Lanczos tangent
                guess = lanczos_guess
                if (guess is None) and (self.prev_lanczos_hash is not None):
                    guess = self.lanczos_tangents[self.prev_lanczos_hash]
                    self.log(
                        f"Using tangent with hash={self.prev_lanczos_hash} "
                        "as initial guess for Lanczos algorithm."
                    )
                # Drop the eigenvalue
                _, tangent, _ = geom_lanczos(ith_image, guess=guess, logger=self.logger)
                self.lanczos_tangents[cur_hash] = tangent
                # Update hash
                self.prev_lanczos_hash = cur_hash
        else:
            raise Exception(f"Unknown tangent kind '{kind}' requested!")

        tangent /= np.linalg.norm(tangent)
        return tangent

    @property
    def tangents(self):
        if self._tangents is None:
            image_energies = self.image_energies
            # Besides kind="lanczos" tangent caluculation is basically free
            self.tangents = np.array(
                [
                    self.get_tangent(i, energies=image_energies)
                    for i in range(len(self.images))
                ]
            )
        return self._tangents

    @tangents.setter
    def tangents(self, tangents):
        assert tangents.shape == (self.nimages, self.coords_length)
        self._tangents = tangents

    def as_xyz(self, comments=None):
        if comments is None:
            comments = [""] * self.nimages
        assert len(comments) == self.nimages
        return "\n".join(
            [
                image.as_xyz(comment=comment)
                for image, comment in zip(self.images, comments)
            ]
        )

    def get_dask_client(self):
        return Client(self.scheduler)

    def get_hei_index(self, energies=None):
        """Return index of highest energy image."""
        if energies is None:
            energies = [image.energy for image in self.images]
        return np.argmax(energies)

    def get_full_cycles(self) -> np.ndarray:
        """Return array of integers that indexes self.all_true_forces/self.all_cart_coords.

        When the ChainOfStates is not yet fully grown this list will be empty.
        The items of this list can be used to index self.all_true_forces and
        related lists, to extract image coordinate & forces data for all
        cycles when the COS was already fully grown.

        This data can then be used to, e.g., construct a (TS)-Hessian for an
        selected image."""

        try:
            fully_grown = self.fully_grown
        except AttributeError:
            fully_grown = True

        if not fully_grown:
            return np.empty(0)

        full_size = self.cart_coords_length * self.nimages
        all_sizes = np.array([true_forces.size for true_forces in self.all_true_forces])
        mask = all_sizes == full_size
        indices = np.arange(len(all_sizes))
        full_cycles = indices[mask]
        return full_cycles

    def prepare_opt_cycle(self, last_coords, last_energies, last_forces):
        """Implements additional logic in preparation of the next
        optimization cycle.

        Should be called by the optimizer at the beginning of a new
        optimization cycle. Can be used to implement additional logic
        as needed for AdaptiveNEB etc.
        """
        self.coords_list.append(last_coords)
        self.forces_list.append(last_forces)

        messages = list()
        # Check if we can start climbing.
        already_climbing = self.started_climbing
        # TODO: use provided forces here to decided whether we want to climb
        # Currently, new forces are calculated here when deciding if we want to climb ...
        if self.climb and not already_climbing:
            # last_image_forces = last_forces.reshape(self.nimages, -1)
            self.started_climbing = self.compare_image_rms_forces(
                self.climb_rms,  # last_image_forces
            )
            if self.started_climbing:
                msg = "Will use climbing image(s) in next cycle."
                self.log(msg)
                messages.append(msg)

        # Fix climbing index/indices if not already set, but requested.
        if already_climbing and self.climb_fixed and (self.fixed_climb_indices is None):
            self.fixed_climb_indices = self.get_climbing_indices()

        # Check if we can start to converge the HEI tangent using the Lanczos algorithm.
        already_climbing_lanczos = self.started_climbing_lanczos
        if (
            self.climb_lanczos
            and self.started_climbing
            and not already_climbing_lanczos
        ):
            self.started_climbing_lanczos = self.compare_image_rms_forces(
                self.climb_lanczos_rms
            )
            if self.started_climbing_lanczos:
                msg = "Will use Lanczos algorithm for HEI tangent in next cycle."
                self.log(msg)
                messages.append(msg)

        # Check for mixing images
        if self.energy_min_mix and not self.started_energy_min_mixing:
            assert all([image.has_all_energies for image in self.images])
            all_energies = np.array([image.all_energies for image in self.images])
            assert all_energies.shape[1] == 2  # Force two states/energies
            energy_diffs = np.diff(all_energies, axis=1).flatten()
            calc_inds = all_energies.argmin(axis=1)

            mix_at = []
            for i, calc_ind in enumerate(calc_inds[:-1]):
                next_ind = calc_inds[i + 1]
                if (
                    (calc_ind != next_ind)
                    and (i not in self.org_forces_indices)
                    and (i + 1 not in self.org_forces_indices)
                ):
                    min_diff_offset = energy_diffs[[i, i + 1]].argmin()
                    mix_at.append(i + min_diff_offset)

            for ind in mix_at:
                self.images[ind].calculator.mix = True
                self.org_forces_indices.append(ind)
                self.log(f"Enable energy mixing for  calculator of image {ind}.")
            self.started_energy_min_mixing = True

        # Starting TS optimization does not influence the return value. Expand on this?!
        if self.ts_opt and not self.started_ts_opt:
            self.started_ts_opt = self.compare_image_rms_forces(self.ts_opt_rms)
            if self.started_ts_opt:
                msg = "Will use TS image(s) in next cycle. Disabled any climbing/Lanczos image."
                self.climb_lanczos = False
                self.climb = False
                self.log(msg)
                messages.append(msg)

        reset_flag = not already_climbing and self.started_climbing
        return reset_flag, messages

    def compare_image_rms_forces(self, ref_rms):
        """Compare rms(forces) value of an image against a reference value.

        Used to decide if we designate an image as climbing image or a
        TS node.

        Only initiate climbing on a sufficiently converged MEP.
        This can be determined from a supplied threshold for the
        RMS force (rms_force) or from a multiple of the
        RMS force convergence threshold (rms_multiple, default).
        """
        # TODO:: del me; this starts force calculations if forces are not set
        # TODO: self.perpendicular_forces is not a property anymore ...
        rms_forces = rms(self.perpendicular_forces)
        # Only start climbing when the COS is fully grown. This
        # attribute may not be defined in all subclasses, so it
        # defaults to True here.
        try:
            fully_grown = self.fully_grown
        except AttributeError:
            fully_grown = True
        start_climbing = (rms_forces <= ref_rms) and fully_grown
        return start_climbing

    def get_climbing_indices(self, image_energies):
        # Index of the highest energy image (HEI)
        hei_index = self.get_hei_index(image_energies)

        move_inds = self.moving_indices
        # Don't climb if not yet enabled or requested.
        if not (self.climb and self.started_climbing):
            climb_indices = tuple()
        elif self.fixed_climb_indices is not None:
            climb_indices = self.fixed_climb_indices
            _ = "index" if len(climb_indices) == 1 else "indices"
            self.log(f"Returning fixed climbing {_}.")
        # Do one image climbing (C1) neb if explicitly requested or
        # the HEI is the first or last item in moving_indices.
        elif self.climb == "one" or ((hei_index == 1) or (hei_index == move_inds[-1])):
            climb_indices = (hei_index,)
        # We can do two climbing (C2) neb if the highest energy image (HEI)
        # is in moving_indices but not the first or last item in this list.
        # elif self.climb != "one" and hei_index in move_inds[1:-1]:
        elif hei_index in move_inds[1:-1]:
            climb_indices = (hei_index - 1, hei_index + 1)
        # Don't climb when the HEI is the first or last image of the whole NEB.
        else:
            climb_indices = tuple()
            self.log("Want to climb but can't. HEI is first or last image!")
        # self.log(f"Climbing indices: {climb_indices}")
        return climb_indices

    def update_with_climbing_forces(
        self, forces, tangents, image_energies, image_forces
    ):
        if not (self.climb and self.started_climbing):
            return forces

        for i in self.get_climbing_indices(image_energies):
            # Energy and forces vector of climbing image
            ienergy = image_energies[i]
            iforces = image_forces[i]
            tangent = tangents[i]
            climbing_forces = iforces - 2 * iforces.dot(tangent) * tangent
            # Updated original forces w/ forces reversed along tangent
            forces[i] = climbing_forces
            norm = np.linalg.norm(climbing_forces)
            self.log(
                f"Climbing with image {i}, E = {ienergy:.6f} au, norm(forces)={norm:.6f}"
            )
        return forces

    def update_with_org_forces(self, forces, indices, image_forces):
        for ind in indices:
            forces[ind] = image_forces[ind]

    def get_ts_image_indices(self):
        if not self.started_ts_opt:
            ts_images = tuple()
        else:
            if self.fixed_ts_indices is None:
                hei_index = self.get_hei_index()
                self.fixed_ts_indices = (hei_index,)
                self.log(f"Fixed TS image to index {hei_index}.")
            ts_images = self.fixed_ts_indices
        return ts_images

    def get_splined_hei(self):
        self.log("Splining HEI")
        # Interpolate energies
        cart_coords = np.array([image.cart_coords for image in self.images])
        if not self.is_analytical_2d:
            cart_coords = align_coords(cart_coords)
        coord_diffs = get_coords_diffs(cart_coords)
        self.log(f"\tCoordinate differences: {coord_diffs}")
        energies = np.array(self.energy)
        energies_spline = interp1d(coord_diffs, energies, kind="cubic")
        x_fine = np.linspace(0, 1, 500)
        energies_fine = energies_spline(x_fine)
        # Determine index that yields the highest energy
        hei_ind = energies_fine.argmax()
        hei_x = x_fine[hei_ind]
        self.log(f"Found splined HEI at x={hei_x:.4f}")
        hei_frac_index = hei_x * (len(self.images) - 1)
        hei_energy = energies_fine[hei_ind]

        reshaped = cart_coords.reshape(-1, self.cart_coords_length)
        # To use splprep we have to transpose the coords.
        transp_coords = reshaped.transpose()
        tcks, us = zip(
            *[
                splprep(transp_coords[i : i + 9], s=0, k=3, u=coord_diffs)
                for i in range(0, len(transp_coords), 9)
            ]
        )

        # Reparametrize mesh
        hei_coords = np.vstack(
            [
                # WTF, Black? This looks horrible.
                splev(
                    [
                        hei_x,
                    ],
                    tck,
                )
                for tck in tcks
            ]
        )
        hei_coords = hei_coords.flatten()

        # Actually it looks like that splined tangents are really bad approximations
        # to the actual imaginary mode. The Cartesian upwinding tangent is usually
        # much much better. In 'run_tsopt_from_cos' we actually mix two "normal" tangents
        # to obtain the HEI tangent.
        hei_tangent = np.vstack(
            [
                # WTF, Black? This looks horrible.
                splev(
                    [
                        hei_x,
                    ],
                    tck,
                    der=1,
                )
                for tck in tcks
            ]
        ).T
        hei_tangent = hei_tangent.flatten()
        hei_tangent /= np.linalg.norm(hei_tangent)
        return hei_coords, hei_energy, hei_tangent, hei_frac_index

    def get_image_calc_counter_sum(self):
        return sum([image.calculator.calc_counter for image in self.images])

    def describe(self):
        imgs = self.images
        img = imgs[0]
        return f"ChainOfStates, {len(imgs)} images, ({img.sum_formula}, {len(img.atoms)} atoms) per image"

    def __str__(self):
        name = self.__class__.__name__
        nimages = len(self.images)
        return f"{name}({nimages} images)"
