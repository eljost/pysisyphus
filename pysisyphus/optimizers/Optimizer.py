import abc
from dataclasses import dataclass
import logging
import os
from pathlib import Path
import sys
import textwrap
import time
from typing import Literal, Optional, Tuple

import numpy as np
import yaml

from pysisyphus.cos.ChainOfStates import ChainOfStates
from pysisyphus.Geometry import Geometry
from pysisyphus.helpers import (
    check_for_end_sign,
    fit_rigid,
    get_coords_diffs,
    procrustes,
)
from pysisyphus.helpers_pure import highlight_text
from pysisyphus.intcoords.exceptions import RebuiltInternalsException
from pysisyphus.intcoords.helpers import interfragment_distance
from pysisyphus.io.hdf5 import get_h5_group, resize_h5_group
from pysisyphus.optimizers.exceptions import ZeroStepLength
from pysisyphus.TablePrinter import TablePrinter


def get_optimizer(geom, index, opt_kwargs, prefix):
    # Import here, to avoid cyclic import error
    from pysisyphus.optimizers.cls_map import get_opt_cls

    opt_kwargs = opt_kwargs.copy()
    opt_type = opt_kwargs.pop("type")

    prefix = f"{prefix}_img_{index:02d}"
    _opt_kwargs = {
        "prefix": prefix,
        "h5_group_name": f"{prefix}_opt",
        # Child optimizers don't print to STDOUT/pysisyphus.log by default,
        # only to their respective log files.
        "logging_level": logging.DEBUG,
    }
    opt_kwargs.update(_opt_kwargs)

    opt = get_opt_cls(opt_type)(geom, **opt_kwargs)
    return opt


def configure_opt_logger(logger, prefix):
    logger.setLevel(logging.DEBUG)
    handler = logging.FileHandler(f"{prefix}optimizer.log", mode="w", delay=True)
    fmt_str = "%(message)s"
    formatter = logging.Formatter(fmt_str)
    handler.setFormatter(formatter)
    handler.setLevel(logging.DEBUG)
    logger.addHandler(handler)


def get_data_model(geometry, is_cos, max_cycles):
    try:
        # Attribute is only present in COS classes
        image_num = geometry.max_image_num
        dummy_geom = geometry.images[0]
    except AttributeError:
        image_num = 1
        dummy_geom = geometry

    # Define dataset shapes. As pysisyphus offers growing COS methods where
    # the number of images changes along the optimization we have to define
    # the shapes accordingly by considering the maximum number of images.
    _1d = (max_cycles,)
    _2d = (max_cycles, image_num * dummy_geom.coords.size)
    _image_inds = (max_cycles, image_num)
    # Number of cartesian coordinates is probably different from the number
    # of internal coordinates.
    _2d_cart = (max_cycles, image_num * dummy_geom.cart_coords.size)
    # The dimensionality of energies depends on whether a COS is optimized or
    # not. I know this is probably not the best idea...
    _energy = _1d if (not is_cos) else (max_cycles, geometry.max_image_num)

    data_model = {
        "image_nums": _1d,
        "image_inds": _image_inds,
        "cart_coords": _2d_cart,
        "coords": _2d,
        "energies": _energy,
        "forces": _2d,
        # AFIR related
        "true_energies": _energy,
        "true_forces": _2d_cart,
        "steps": _2d,
        # Convergence related
        "max_forces": _1d,
        "rms_forces": _1d,
        "max_steps": _1d,
        "rms_steps": _1d,
        # Misc
        "cycle_times": _1d,
        "modified_forces": _2d,
        # COS specific
        "tangents": _2d,
    }

    return data_model


CONV_THRESHS = {
    #                max_force, rms_force, max_step, rms_step
    "nwchem_loose": (4.5e-3, 3.0e-3, 5.4e-3, 3.6e-3),
    "gau_loose": (2.5e-3, 1.7e-3, 1.0e-2, 6.7e-3),
    "gau": (4.5e-4, 3.0e-4, 1.8e-3, 1.2e-3),
    "gau_tight": (1.5e-5, 1.0e-5, 6.0e-5, 4.0e-5),
    "gau_vtight": (2.0e-6, 1.0e-6, 6.0e-6, 4.0e-6),
    "baker": (3.0e-4, 2.0e-4, 3.0e-4, 2.0e-4),
    # Dummy thresholds
    "never": (2.0e-6, 1.0e-6, 6.0e-6, 4.0e-6),
}
Thresh = Literal["gau_loose", "gau", "gau_tight", "gau_vtight", "baker", "never"]


@dataclass(frozen=True)
class ConvInfo:
    cur_cycle: int
    energy_converged: bool
    max_force_converged: bool
    rms_force_converged: bool
    max_step_converged: bool
    rms_step_converged: bool
    desired_eigval_structure: bool

    def get_convergence(self):
        return (
            self.energy_converged,
            self.max_force_converged,
            self.rms_force_converged,
            self.max_step_converged,
            self.rms_step_converged,
            self.desired_eigval_structure,
        )

    def is_converged(self):
        return all(self.get_convergence())


class Optimizer(metaclass=abc.ABCMeta):
    def __init__(
        self,
        geometry: Geometry,
        thresh: Thresh = "gau_loose",
        max_step: float = 0.04,
        max_cycles: int = 150,
        min_step_norm: float = 1e-8,
        assert_min_step: bool = True,
        rms_force: Optional[float] = None,
        rms_force_only: bool = False,
        max_force_only: bool = False,
        force_only: bool = False,
        converge_to_geom_rms_thresh: float = 0.05,
        align: bool = False,
        align_factor: float = 1.0,
        dump: bool = False,
        dump_restart: bool = False,
        print_every: int = 1,
        dump_ascii: bool = False,
        prefix: str = "",
        reparam_thresh: float = 1e-3,
        reparam_check_rms: bool = True,
        # reparam_when: Optional[Literal["before", "after"]] = "after",
        reparam_when: Optional[Literal["after", None]] = "after",
        overachieve_factor: float = 0.0,
        check_eigval_structure: bool = False,
        restart_info=None,
        check_coord_diffs: bool = True,
        coord_diff_thresh: float = 0.01,
        fragments: Optional[Tuple] = None,
        monitor_frag_dists: int = 0,
        out_dir: str = ".",
        h5_fn: str = "optimization.h5",
        h5_group_name: str = "opt",
        image_opt_kwargs: Optional[dict] = None,
        logging_level: int = logging.INFO,
    ) -> None:
        """Optimizer baseclass. Meant to be subclassed.

        Parameters
        ----------
        geometry
            Geometry to be optimized.
        thresh
            Convergence threshold.
        max_step
            Maximum absolute component of the allowed step vector. Utilized in
            optimizers that don't support a trust region or line search.
        max_cycles
            Maximum number of allowed optimization cycles.
        min_step_norm
            Minimum norm of an allowed step. If the step norm drops below
            this value a ZeroStepLength-exception is raised. The unit depends
            on the coordinate system of the supplied geometry.
        assert_min_step
            Flag that controls whether the norm of the proposed step is check
            for being too small.
        rms_force
            Root-mean-square of the force from which user-defined thresholds
            are derived. When 'rms_force' is given 'thresh' is ignored.
        rms_force_only
            When set, convergence is signalled only based on rms(forces).
        max_force_only
            When set, convergence is signalled only based on max(|forces|).
        force_only
            When set, convergence is signalled only based on max(|forces|) and rms(forces).
        converge_to_geom_rms_thresh
            Threshold for the RMSD with another geometry. When the RMSD drops
            below this threshold convergence is signalled. Only used with
            Growing Newton trajectories.
        align
            Flag that controls whether the geometry is aligned in every step
            onto the coordinates of the previous step. Must not be used with
            internal coordinates.
        align_factor
            Factor that controls the strength of the alignment. 1.0 means
            full alignment, 0.0 means no alignment. The factor mixes the
            rotation matrix of the alignment with the identity matrix.
        dump
            Flag to control dumping/writing of optimization progress to the
            filesystem
        dump_restart
            Flag to control whether restart information is dumped to the
            filesystem.
        print_every
            Report optimization progress every nth cycle.
        dump_ascii
            Dump an ASCII representation of the optimized geometry to the
            filesystem. Has no effect for COS optimizations.
        prefix
            Short string that is prepended to several files created by
            the optimizer. Allows distinguishing several optimizations carried
            out in the same directory.
        reparam_thresh
            Controls the minimal allowed similarity between coordinates
            after two successive reparametrizations. Convergence is signalled
            if the coordinates did not change significantly.
        reparam_check_rms
            Whether to check for (too) similar coordinates after reparametrization.
        reparam_when
            Reparametrize before or after calculating the step. Can also be turned
            off by setting it to None.
        overachieve_factor
            Signal convergence when max(forces) and rms(forces) fall below the
            chosen threshold, divided by this factor. Convergence of max(step) and
            rms(step) is ignored.
        check_eigval_structure
            Check the eigenvalues of the modes we maximize along. Convergence requires
            them to be negative. Useful if TS searches are started from geometries close
            to a minimum.
        restart_info
            Restart information. Undocumented.
        check_coord_diffs
            Whether coordinates of chain-of-sates images are checked for being
            too similar.
        coord_diff_thresh
            Unitless threshold for similary checking of COS image coordinates.
            The first image is assigned 0, the last image is assigned to 1.
        fragments
            Tuple of lists containing atom indices, defining two fragments.
        monitor_frag_dists
            Monitor fragment distances for N cycles. The optimization is terminated
            when the interfragment distances falls below the initial value after N
            cycles.
        out_dir
            String poiting to a directory where optimization progress is
            dumped.
        h5_fn
            Basename of the HDF5 file used for dumping.
        h5_group_name
            Groupname used for dumping of this optimization.
        image_opt_kwargs
            Optimizer kwargs that are used to construct child-optimizers. Currently
            only used to obtain TS-optimizers in COS-optimizations to converge an
            image to an actual TS.
        logging_level
            Controls logging level of TablePrinter. Setting it to logging.DEBUG will
            suppress printing to STDOUT and writing to pysisyphus.log. Messages will
            be recorded into a respective f'{prefix}optimizers.log' file, regardless
            of this argument.
        """
        assert thresh in CONV_THRESHS.keys()

        self.geometry = geometry
        self.thresh = thresh
        self.max_step = max_step
        self.min_step_norm = min_step_norm
        self.assert_min_step = assert_min_step
        self.rms_force_only = rms_force_only
        self.max_force_only = max_force_only
        self.force_only = force_only
        self.converge_to_geom_rms_thresh = converge_to_geom_rms_thresh
        self.align = align
        self.align_factor = align_factor
        self.dump = dump
        self.dump_restart = dump_restart
        print_every = int(print_every)
        assert print_every >= 1
        self.print_every = print_every
        self.dump_ascii = dump_ascii
        self.prefix = f"{prefix}_" if prefix else prefix
        self.reparam_thresh = reparam_thresh
        self.reparam_check_rms = reparam_check_rms
        self.reparam_when = reparam_when
        assert self.reparam_when in ("after", None)
        self.overachieve_factor = float(overachieve_factor)
        self.check_eigval_structure = check_eigval_structure
        self.check_coord_diffs = check_coord_diffs
        self.coord_diff_thresh = float(coord_diff_thresh)
        if image_opt_kwargs is None:
            image_opt_kwargs = dict()
        self.image_opt_kwargs = image_opt_kwargs
        self.logging_level = logging_level

        self.logger = logging.getLogger(f"pysis.{self.prefix}opt")
        configure_opt_logger(self.logger, self.prefix)
        self.is_cos = issubclass(type(self.geometry), ChainOfStates)

        # Set up convergence thresholds
        self.convergence = self.make_conv_dict(
            thresh, rms_force, rms_force_only, max_force_only, force_only
        )
        for key, value in self.convergence.items():
            setattr(self, key, value)

        if self.thresh == "never":
            max_cycles = 1_000_000_000
            self.dump = False
            self.log(
                f"Got threshold {self.thresh}, set 'max_cycles' to {max_cycles} "
                "and disabled dumping!"
            )
            self.conv_dict = {}
        try:
            self.converge_to_geom = self.geometry.converge_to_geom
        except AttributeError:
            # TODO: log that attribute is not present at debug level
            self.converge_to_geom = None
        self.max_cycles = max_cycles

        self.fragments = fragments
        self.monitor_frag_dists = monitor_frag_dists
        if self.monitor_frag_dists:
            assert (
                len(self.fragments) == 2
            ), "Interfragment monitoring requires two fragments!"
            assert all(
                [len(frag) > 0 for frag in self.fragments]
            ), "Fragments must not be empty!"
        # Setting some default values
        self.monitor_frag_dists_counter = self.monitor_frag_dists
        self.interfrag_dists = list()
        self.resetted = False
        try:
            out_dir = Path(out_dir)
        except TypeError:
            out_dir = Path(".")
        self.out_dir = out_dir.resolve()
        self.out_dir.mkdir(parents=True, exist_ok=True)

        if self.is_cos:
            moving_image_num = len(self.geometry.moving_indices)
            self.print(
                f"Path with {moving_image_num} moving images "
                f"(first fixed: {self.geometry.fix_first}, last fixed: {self.geometry.fix_last})."
            )

        # Don't use prefix for this fn, as different optimizations
        # can be distinguished according to their group in the HDF5 file.
        self.h5_fn = self.get_path_for_fn(h5_fn, with_prefix=False)
        # TODO: make this also depend on 'prefix'?!
        self.h5_group_name = h5_group_name

        current_fn = "current_geometries.trj" if self.is_cos else "current_geometry.xyz"
        self.current_fn = self.get_path_for_fn(current_fn)
        final_fn = "final_geometries.trj" if self.is_cos else "final_geometry.xyz"
        self.final_fn = self.get_path_for_fn(final_fn)
        self.hei_trj_fn = self.get_path_for_fn("cos_hei.trj")
        try:
            os.remove(self.hei_trj_fn)
        except FileNotFoundError:
            pass

        # Setting some empty lists as default. The actual shape of the respective
        # entries is not considered, which gives us some flexibility.
        self.data_model = get_data_model(self.geometry, self.is_cos, self.max_cycles)
        for la in self.data_model.keys():
            setattr(self, la, list())

        if self.dump:
            out_trj_fn = self.get_path_for_fn("optimization.trj")
            self.out_trj_handle = open(out_trj_fn, "w")
            # Call with reset=True to delete remnants of previous calculations, unless
            # the optimizer was restarted. Given a previous optimization with, e.g. 30
            # cycles and a second restarted optimization with 20 cycles the last 10 cycles
            # of the previous optimization would still be present.
            reset = restart_info is None
            h5_group = get_h5_group(
                self.h5_fn, self.h5_group_name, self.data_model, reset=reset
            )
            h5_group.file.close()
        if self.prefix:
            self.log(f"Created optimizer with prefix {self.prefix}")

        self.restarted = False
        self.last_cycle = 0
        self.cur_cycle = 0
        if restart_info is not None:
            if isinstance(restart_info, str):
                restart_info = yaml.load(restart_info, Loader=yaml.SafeLoader)
            self.set_restart_info(restart_info)
            self.restarted = True
        self.image_opts = dict()

        header = "cycle Δ(energy) max(|force|) rms(force) max(|step|) rms(step) s/cycle".split()
        col_fmts = "int float float float float float float_short".split()
        self.table = TablePrinter(
            header, col_fmts, width=12, logger=self.logger, level=self.logging_level
        )
        self.is_converged = False

    def get_path_for_fn(self, fn, with_prefix=True):
        prefix = self.prefix if with_prefix else ""
        return self.out_dir / (prefix + fn)

    def make_conv_dict(
        self, key, rms_force=None, rms_force_only=False, max_force_only=False, force_only=False
    ):
        if not rms_force:
            threshs = CONV_THRESHS[key]
        else:
            self.print(
                "Deriving convergence threshold from supplied "
                f"rms_force={rms_force}."
            )
            threshs = (
                1.5 * rms_force,
                rms_force,
                6 * rms_force,
                4 * rms_force,
            )
        keys = keep_keys = [
            "max_force_thresh",
            "rms_force_thresh",
            "max_step_thresh",
            "rms_step_thresh",
        ]
        conv_dict = {k: v for k, v in zip(keys, threshs)}

        # Only used gradient information for COS optimizations
        if self.is_cos:
            keep_keys = ["max_force_thresh", "rms_force_thresh"]

        if rms_force_only:
            self.log("Checking convergence with rms(forces) only!")
            keep_keys = ["rms_force_thresh"]
        elif max_force_only:
            self.log("Checking convergence with max(forces) only!")
            keep_keys = ["max_force_thresh"]
        elif force_only:
            self.log("Checking convergence with max(forces) and rms(forces) only!")
            keep_keys = ["max_force_thresh", "rms_force_thresh"]

        # The dictionary should only contain pairs that are needed
        conv_dict = {key: value for key, value in conv_dict.items() if key in keep_keys}
        return conv_dict

    def report_conv_thresholds(self):
        oaf = self.overachieve_factor

        # Overachieved
        def oa(val):
            return f", ({val/oaf:.6f})" if oaf > 0.0 else ""

        internal_coords = self.geometry.coord_type not in (
            "cart",
            "cartesian",
            "mwcartesian",
        )
        fu = "E_h a_0⁻¹" + (" (rad⁻¹)" if internal_coords else "")  # forces unit
        su = "a_0" + (" (rad)" if internal_coords else "")  # step unit

        try:
            rms_thresh = f"\tmax(|force|) <= {self.max_force_thresh:.6f}{oa(self.max_force_thresh)} {fu}"
        except AttributeError:
            rms_thresh = None
        try:
            max_thresh = f"\t  rms(force) <= {self.rms_force_thresh:.6f}{oa(self.rms_force_thresh)} {fu}"
        except AttributeError:
            max_thresh = None
        threshs = (rms_thresh, max_thresh)

        if self.rms_force_only:
            use_threshs = (threshs[1],)
        elif self.max_force_only:
            use_threshs = (threshs[0],)
        elif self.force_only:
            use_threshs = threshs
        elif self.is_cos:
            use_threshs = threshs
        else:
            use_threshs = threshs + (
                f"\t max(|step|) <= {self.max_step_thresh:.6f} {su}",
                f"\t   rms(step) <= {self.rms_step_thresh:.6f} {su}",
            )
        self.print(
            "Convergence thresholds"
            + (", (overachieved when)" if oaf > 0.0 else "")
            + ":\n"
            + "\n".join(use_threshs)
            + f"\n\t'Superscript {self.table.mark}' indicates convergence"
            + "\n"
        )

    def log(self, message, level=logging.DEBUG):
        """Log to optimizer-specific logfile.

        By default, messages passed to Optimizer.log will only appear in the
        optimizer-specific logfile 'f{self.prefix}optimizer.log'."""
        self.logger.log(level, message)

    def print(self, message):
        """Print/write to STDOUT, pysisyphus.log and {prefix}optimizer.log.

        Whether this method actually prints/writes depends on self.logging_level.
        If the level is below logging.INFO (20) 'msg' will not appear in STDOUT and
        pysisyphus.log. Until logging.DEBUG (10) it will appear in the optimizer-
        specific log file.

        This method should only be used in the Optimizer baseclass (this module).
        Optimizers inheriting from this class should use Optimizer.log() instead.
        """
        self.logger.log(self.logging_level, message)

    def check_convergence(self, step=None, multiple=1.0, overachieve_factor=None):
        """Check if the current convergence of the optimization
        is equal to or below the required thresholds, or a multiple
        thereof. The latter may be used in initiating the climbing image.
        """

        if step is None:
            step = self.steps[-1]
        if overachieve_factor is None:
            overachieve_factor = self.overachieve_factor

        # When using a ChainOfStates method convergence depends on the
        # perpendicular component of the force along the MEP.
        if self.is_cos:
            forces = self.geometry.perpendicular_forces
        elif len(self.modified_forces) == len(self.forces):
            self.log("Using modified forces to determine convergence!")
            forces = self.modified_forces[-1]
        else:
            forces = self.forces[-1]

        # The forces of fixed images may be zero and this may distort the RMS
        # values. So we take into account the number of moving images with
        # non-zero forces vectors.
        if self.is_cos:
            non_zero_elements = (
                len(self.geometry.moving_indices) * self.geometry.coords_length
            )
            rms_force = np.sqrt(np.sum(np.square(forces)) / non_zero_elements)
            rms_step = np.sqrt(np.sum(np.square(step)) / non_zero_elements)
        else:
            rms_force = np.sqrt(np.mean(np.square(forces)))
            rms_step = np.sqrt(np.mean(np.square(step)))

        max_force = np.abs(forces).max()
        max_step = np.abs(step).max()

        self.max_forces.append(max_force)
        self.rms_forces.append(rms_force)
        self.max_steps.append(max_step)
        self.rms_steps.append(rms_step)

        # Give geometry a chance to signal convergence, e.g. GrowingNT that
        # is supposed to stop when a TS was passed.
        try:
            geom_converged = self.geometry.check_convergence()
        except AttributeError:
            geom_converged = False

        converged_to_geom = False
        if self.converge_to_geom is not None:
            rmsd = np.sqrt(
                np.mean((self.converge_to_geom.coords - self.geometry.coords) ** 2)
            )
            converged_to_geom = rmsd < self.converge_to_geom_rms_thresh

        this_cycle = {
            "max_force_thresh": max_force,
            "rms_force_thresh": rms_force,
            "max_step_thresh": max_step,
            "rms_step_thresh": rms_step,
        }

        def check(key):
            # Always return True if given key is not checked
            key_is_checked = key in self.convergence
            if key_is_checked:
                result = this_cycle[key] <= getattr(self, key) * multiple
            else:
                result = True
            return result

        convergence = {
            "energy_converged": True,
            "max_force_converged": check("max_force_thresh"),
            "rms_force_converged": check("rms_force_thresh"),
            "max_step_converged": check("max_step_thresh"),
            "rms_step_converged": check("rms_step_thresh"),
        }
        # For TS optimizations we also try to check the eigenvalue structure of the
        # Hessian. A saddle point of order n must have exatcly only n significant negative
        # eigenvalues. We try to check this below.
        #
        # Currently, this is not totally strict,
        # as only the values in self.ts_mode_eigvals are checked but actually all eigenvalues
        # would have to be checked.
        desired_eigval_structure = True
        if self.check_eigval_structure:
            try:
                desired_eigval_structure = (
                    # Acutally all eigenvalues would have to be checked, but
                    # currently they are not stored anywhere.
                    self.ts_mode_eigvals < self.small_eigval_thresh
                ).sum() == len(self.roots)
            except AttributeError:
                self.log(
                    "Skipping check of eigenvalue structure, as information is unavailable."
                )
        convergence["desired_eigval_structure"] = desired_eigval_structure
        conv_info = ConvInfo(self.cur_cycle, **convergence)

        # Check if force convergence is overachieved. If the eigenvalue-structure is
        # checked, a wrong structure will prevent overachieved convergence.
        overachieved = False
        if overachieve_factor > 0 and desired_eigval_structure:
            max_thresh = self.convergence.get("max_force_thresh", 0) / overachieve_factor
            rms_thresh = self.convergence.get("rms_force_thresh", 0) / overachieve_factor
            max_ = max_force < max_thresh
            rms_ = rms_force < rms_thresh
            overachieved = max_ and rms_
            if max_:
                self.log("max(force) is overachieved")
            if rms_:
                self.log("rms(force) is overachieved")
            if max_ and rms_:
                self.log("Force convergence overachieved!")

        converged = all([val for val in convergence.values()])
        not_never = self.thresh != "never"

        if self.thresh == "baker":
            energy_converged = False
            if self.cur_cycle > 0:
                cur_energy = self.energies[-1]
                prev_energy = self.energies[-2]
                energy_converged = abs(cur_energy - prev_energy) < 1e-6
            converged = (max_force < 3e-4) and (energy_converged or (max_step < 3e-4))
        return (
            any((converged_to_geom, converged, overachieved, geom_converged))
            and not_never,
            conv_info,
        )

    def print_opt_progress(self, conv_info):
        try:
            energy_diff = self.energies[-1] - self.energies[-2]
        # ValueError: maybe raised when the number of images differ in cycles
        # IndexError: raised in first cycle when only one energy is present
        except (ValueError, IndexError):
            energy_diff = float("nan")

        # Try to sum COS energies
        try:
            energy_diff = sum(energy_diff)
        except TypeError:
            pass

        if (self.cur_cycle > 1) and (self.cur_cycle % 10 == 0):
            self.table.print_sep()

        # desired_eigval_structure is also returned, but currently not reported.
        marks = [False, *conv_info.get_convergence()[:-1], False]
        try:
            cycle_time = self.cycle_times[-1]
        except IndexError:
            cycle_time = 0.0
        self.table.print_row(
            (
                self.cur_cycle,
                energy_diff,
                self.max_forces[-1],
                self.rms_forces[-1],
                self.max_steps[-1],
                self.rms_steps[-1],
                cycle_time,
            ),
            marks=marks,
        )
        try:
            # Geometries/ChainOfStates objects can also do some printing.
            add_info = self.geometry.get_additional_print()
            if add_info:
                self.table.print(add_info)
        except AttributeError:
            pass

    @property
    def is_cart_opt(self):
        return self.geometry.coord_type in ("cart", "cartesian", "mw_cartesian")

    def _fit_rigid(self, *, vectors=None, vector_lists=None, hessian=None):
        if vector_lists is None:
            vector_lists = list()

        # Always rotate cartesian coordinates and forces
        base_vector_lists = [self.cart_coords, self.forces, self.steps]
        vector_lists = base_vector_lists + vector_lists
        (
            rot_vecs,
            (rot_cart_coords, rot_forces, rot_steps, *rot_vec_lists),
            rot_hessian,
        ) = fit_rigid(
            self.geometry,
            vectors=vectors,
            vector_lists=vector_lists,
            hessian=hessian,
            align_factor=self.align_factor,
        )
        self.cart_coords = rot_cart_coords
        self.forces = rot_forces
        self.steps = rot_steps
        return rot_vecs, rot_vec_lists, rot_hessian

    def _procrustes(self):
        """Wrapper for procrustes that passes additional arguments along."""
        procrustes(self.geometry, self.align_factor)

    def scale_by_max_step(self, steps):
        steps_max = np.abs(steps).max()
        if steps_max > self.max_step:
            steps *= self.max_step / steps_max
        return steps

    def prepare_opt(self):
        pass

    def postprocess_opt(self):
        pass

    @abc.abstractmethod
    def get_step(
        self,
        energy,
        forces: np.ndarray,
        hessian: Optional[np.ndarray] = None,
        eigvals: Optional[np.ndarray] = None,
        eigvecs: Optional[np.ndarray] = None,
        resetted: Optional[bool] = None,
    ):
        pass

    def housekeeping(self):
        # Calculate energy and forces
        forces = self.geometry.forces
        energy = self.geometry.energy
        self.forces.append(forces)
        self.energies.append(energy)

        if not self.is_cos:
            self.log(f"      Energy: {energy: >12.6f} au")
        self.log(f"norm(forces): {np.linalg.norm(forces): >12.6f} au / bohr (rad)")
        self.log(f" rms(forces): {np.sqrt(np.mean(forces**2)): >12.6f} au / bohr (rad)")

        return {
            "energy": energy,
            "forces": forces,
        }

    def update_image_step(self, image_index, prefix=""):
        """Calculate a new step for an image in a COS using another Optimizer.

        This method is useful to update an image step, e.g, by a TS optimizer.
        """
        if prefix != "":
            prefix = f"{prefix}-"
        geom = self.geometry.images[image_index]
        # Try to look up the Optimizer ...
        try:
            opt = self.image_opts[image_index]
        # ... if it is not already present we create it and store it.
        except KeyError:
            opt = get_optimizer(
                geom, image_index, opt_kwargs=self.image_opt_kwargs, prefix="ts"
            )
            self.image_opts[image_index] = opt
            self.print(f"Created new {prefix}optimizer for image {image_index}")

        try:
            # Advance the generator one cycle and fetch the new step
            next(opt.run_generator)
            updated_step = opt.steps[-1]
        except StopIteration:
            # When the optimzation already converged the generator is exhausted
            # and we return a zero-step.
            updated_step = np.zeros_like(self.geometry.images[0].coords)

        return updated_step

    def write_to_out_dir(self, out_fn, content, mode="w"):
        out_path = self.out_dir / out_fn
        with open(out_path, mode) as handle:
            handle.write(content)

    def write_image_trjs(self):
        base_name = "image_{:03d}.trj"
        for i, image in enumerate(self.geometry.images):
            image_fn = base_name.format(i)
            comment = f"cycle {self.cur_cycle}"
            as_xyz = image.as_xyz(comment)
            self.write_to_out_dir(image_fn, as_xyz + "\n", "a")

    def write_results(self):
        # Save results from the Optimizer to HDF5 file if requested
        h5_group = get_h5_group(self.h5_fn, self.h5_group_name)

        # Some attributes never change and are only set in the first cycle
        if self.cur_cycle == 0:
            h5_group.attrs["is_cos"] = self.is_cos
            try:
                atoms = self.geometry.images[0].atoms
                coord_size = self.geometry.images[0].coords.size
            except AttributeError:
                atoms = self.geometry.atoms
                coord_size = self.geometry.coords.size
            h5_group.attrs["atoms"] = np.bytes_(atoms)
            h5_group.attrs["coord_type"] = self.geometry.coord_type
            h5_group.attrs["coord_size"] = coord_size
            h5_group.attrs["overachieve_factor"] = self.overachieve_factor
            for key in (
                "max_force_thresh",
                "rms_force_thresh",
                "max_step_thresh",
                "rms_step_thresh",
            ):
                try:
                    h5_group.attrs[key] = self.convergence[key]
                # Step threshold may not be present
                except KeyError:
                    pass

        # Update changing attributes
        h5_group.attrs["cur_cycle"] = self.cur_cycle
        h5_group.attrs["is_converged"] = self.is_converged

        for key, shape in self.data_model.items():
            value = getattr(self, key)
            # Don't try to set empty values, e.g. 'tangents' are only present
            # for COS methods. 'modified_forces' are only present for NCOptimizer.
            if not value:
                continue
            if len(shape) > 1:
                h5_group[key][self.cur_cycle, : len(value[-1])] = value[-1]
            else:
                h5_group[key][self.cur_cycle] = value[-1]

        h5_group.file.close()

    def write_cycle_to_file(self):
        as_xyz_str = self.geometry.as_xyz()

        if self.is_cos:
            out_fn = "cycle_{:03d}.trj".format(self.cur_cycle)
            self.write_to_out_dir(out_fn, as_xyz_str)
            # Also write separate .trj files for every image in the cos
            self.write_image_trjs()

            # Dump current HEI
            max_ind = np.argmax(self.energies[-1])
            with open(self.hei_trj_fn, "a") as handle:
                handle.write(self.geometry.images[max_ind].as_xyz() + "\n")

        else:
            # Append to .trj file
            self.out_trj_handle.write(as_xyz_str + "\n")
            self.out_trj_handle.flush()
        # Dump to HDF5
        self.write_results()

    def final_summary(self):
        # If the optimization was stopped _forces may not be set, so
        # then we force a calculation if it was not already set.
        _ = self.geometry.forces
        cart_forces = self.geometry.cart_forces
        max_cart_forces = np.abs(cart_forces).max()
        rms_cart_forces = np.sqrt(np.mean(cart_forces**2))
        int_str = ""
        if self.geometry.coord_type not in ("cart", "cartesian", "mwcartesian"):
            int_forces = self.geometry.forces
            max_int_forces = np.abs(int_forces).max()
            rms_int_forces = np.sqrt(np.mean(int_forces**2))
            int_str = f"""
            \tmax(forces, internal): {max_int_forces:.6f} hartree/(bohr,rad)
            \trms(forces, internal): {rms_int_forces:.6f} hartree/(bohr,rad)"""
        energy = self.geometry.energy
        final_summary = f"""
        Final summary:{int_str}
        \tmax(forces,cartesian): {max_cart_forces:.6f} hartree/bohr
        \trms(forces,cartesian): {rms_cart_forces:.6f} hartree/bohr
        \tenergy: {energy:.8f} hartree
        """
        return textwrap.dedent(final_summary.strip())

    def get_run_generator(self):
        self.print("If not specified otherwise, all quantities are given in au.\n")

        if not self.restarted:
            prep_start_time = time.time()
            # Set up the optimizer in prepare_opt, e.g., prepare initial Hessiane etc.
            print("@ Calling prepare_opt()")
            self.prepare_opt()
            print("@ Done w/ prepare_opt()")
            self.log(f"{self.geometry.coords.size} degrees of freedom.")
            prep_end_time = time.time()
            prep_time = prep_end_time - prep_start_time
            self.report_conv_thresholds()
            self.print(f"Spent {prep_time:.1f} s preparing the first cycle.")

        self.table.print_header()
        self.stopped = False
        # Actual optimization loop
        for self.cur_cycle in range(self.last_cycle, self.max_cycles):
            start_time = time.time()
            self.log("\n" + highlight_text(f"Cycle {self.cur_cycle:03d}"))

            if self.is_cos and self.check_coord_diffs:
                image_coords = [image.cart_coords for image in self.geometry.images]
                align = len(image_coords[0]) > 3
                cds = get_coords_diffs(image_coords, align=align)
                # Differences of coordinate differences ;)
                cds_diffs = np.diff(cds)
                min_ind = cds_diffs.argmin()
                if cds_diffs[min_ind] < self.coord_diff_thresh:
                    similar_inds = min_ind, min_ind + 1
                    msg = (
                        f"Cartesian coordinates of images {similar_inds} are "
                        "too similar. Stopping optimization!"
                    )
                    # I should improve my logging :)
                    self.print(msg)
                    self.log(msg)
                    sim_fn = "too_similar.trj"
                    with open(sim_fn, "w") as handle:
                        handle.write(self.geometry.as_xyz())
                    self.print(f"Dumped latest coordinates to '{sim_fn}'.")
                    break

            # Check if something considerably changed in the optimization, e.g.,
            # new images were added/interpolated. Then the optimizer should be reset.
            reset_flag = False
            if self.cur_cycle > 0 and self.is_cos:
                reset_flag, messages = self.geometry.prepare_opt_cycle(
                    self.coords[-1], self.energies[-1], self.forces[-1]
                )
                if messages:
                    self.print("\n".join(messages))
            # Reset when number of coordinates changed, compared to the previous cycle.
            elif self.cur_cycle > 0:
                reset_flag = reset_flag or (
                    self.geometry.coords.size != self.coords[-1].size
                )

            if reset_flag:
                self.reset()

            """
            # Coordinates may be updated here.
            if self.reparam_when == "before" and hasattr(
                self.geometry, "reparametrize"
            ):
                # This call actually returns a bool, but right now we just drop it.
                self.geometry.reparametrize()
            """

            if self.is_cos and self.align and self.is_cart_opt:
                # Try to use fit_rigid() of the optimizer class, if implemented;
                # fall back to the generic _fit_rigid() method else.
                try:
                    self.fit_rigid()
                except AttributeError:
                    self._fit_rigid()

            self.coords.append(self.geometry.coords.copy())
            self.cart_coords.append(self.geometry.cart_coords.copy())

            # Determine and store number of currenctly actively optimized images.
            # In non-COS optimizations we store some dummy values [0] and 1, resembling
            # a COS with one image.
            try:
                image_inds = self.geometry.image_inds
                image_num = len(image_inds)
            except AttributeError:
                image_inds = [0]
                image_num = 1
            self.image_inds.append(image_inds)
            self.image_nums.append(image_num)

            #####################################################################
            #  Actual step calculation in the Optimizer subclass is done here!  #
            #  The step is not yet taken in the underlying Geometry/COS object. #
            #####################################################################

            # Calculate step
            get_step_kwargs = self.housekeeping()
            step = self.get_step(**get_step_kwargs)

            try:
                ts_image_indices = self.geometry.get_ts_image_indices()
                # Create view on step with per image steps
                image_steps = step.reshape(self.geometry.nimages, -1)
            except AttributeError:
                ts_image_indices = list()
            for ts_image_index in ts_image_indices:
                # Update view and underlying step array
                image_steps[ts_image_index] = self.update_image_step(
                    ts_image_index, prefix="TS"
                )

            step_norm = np.linalg.norm(step)
            self.log(f"norm(step)={step_norm:.6f} au (rad)")
            for source, target in (
                ("true_energy", "true_energies"),
                ("true_forces", "true_forces"),
            ):
                try:
                    if (value := getattr(self.geometry, source)) is not None:
                        getattr(self, target).append(value)
                except AttributeError:
                    pass

            if step is None:
                # Remove the previously added coords
                self.coords.pop(-1)
                self.cart_coords.pop(-1)
                continue

            if self.is_cos:
                self.tangents.append(self.geometry.get_tangents().flatten())

            self.steps.append(step)

            # Convergence check
            self.is_converged, conv_info = self.check_convergence()

            end_time = time.time()
            elapsed_seconds = end_time - start_time
            self.cycle_times.append(elapsed_seconds)

            if self.dump:
                self.write_cycle_to_file()
                with open(self.current_fn, "w") as handle:
                    handle.write(self.geometry.as_xyz())

            if (
                self.dump
                and self.dump_restart
                and (self.cur_cycle % self.dump_restart) == 0
            ):
                self.dump_restart_info()

            if (self.cur_cycle % self.print_every) == 0 or self.is_converged:
                self.print_opt_progress(conv_info)

            if self.dump_ascii and not self.is_cos:
                asciiart = self.geometry.as_ascii_art()
                if asciiart != "":
                    ascii_fn = self.get_path_for_fn("geom_ascii")
                    with open(ascii_fn, "w") as handle:
                        handle.write(asciiart)

            if self.is_converged:
                self.table.print("Converged!")
                break
            # Allow convergence, before checking for too small steps
            elif self.assert_min_step and (step_norm <= self.min_step_norm):
                raise ZeroStepLength

            # Now actually apply the step to the stored Geometry/ChainOfStates object.
            new_coords = self.geometry.coords.copy() + step
            try:
                self.geometry.coords = new_coords
                # Use the actual step. It may differ from the proposed step
                # when internal coordinates are used, as the internal-Cartesian
                # transformation is done iteratively.
                self.steps[-1] = self.geometry.coords - self.coords[-1]
            except RebuiltInternalsException as exception:
                self.print("Rebuilt internal coordinates!")
                rebuilt_fn = self.get_path_for_fn("rebuilt_primitives.xyz")
                with open(rebuilt_fn, "w") as handle:
                    handle.write(self.geometry.as_xyz())
                if self.is_cos:
                    for image in self.geometry.images:
                        image.reset_coords(exception.typed_prims)
                self.reset()

                if self.dump:
                    self.data_model = get_data_model(
                        self.geometry, self.is_cos, self.max_cycles
                    )
                    self.h5_group_name += "_rebuilt"
                    h5_group = get_h5_group(
                        self.h5_fn,
                        self.h5_group_name,
                        self.data_model,
                        reset=True,
                    )
                    h5_group.file.close()

            # Coordinates may be updated here.
            if (self.reparam_when == "after") and hasattr(
                self.geometry, "reparametrize"
            ):
                reparametrized = self.geometry.reparametrize()
            else:
                reparametrized = False

            cur_coords = self.geometry.coords
            prev_coords = self.coords[-1]
            if (
                self.is_cos
                and self.reparam_check_rms
                and reparametrized
                and (cur_coords.size == prev_coords.size)
            ):
                self.log("Did reparametrization")

                rms = np.sqrt(np.mean((prev_coords - cur_coords) ** 2))
                self.log(f"rms of coordinates after reparametrization={rms:.6f}")
                self.is_converged = rms < self.reparam_thresh
                if self.is_converged:
                    self.table.print(
                        "Insignificant coordinate change after "
                        "reparametrization. Signalling convergence!"
                    )
                    break

            # Alternative: calcualte overlap of AFIR force and step. If this
            # overlap is negative the step is taken against the AFIR force.
            if self.monitor_frag_dists_counter > 0:
                interfrag_dist = interfragment_distance(
                    *self.fragments, self.geometry.coords3d
                )
                try:
                    prev_interfrag_dist = self.interfrag_dists[-1]
                    if interfrag_dist > prev_interfrag_dist:
                        self.print("Interfragment distances increased!")
                        self.stopped = True  # Can't use := in if clause
                        break
                except IndexError:
                    pass
                self.interfrag_dists.append(interfrag_dist)
                self.monitor_frag_dists_counter -= 1

            sys.stdout.flush()
            sign = check_for_end_sign()
            if sign == "stop":
                self.stopped = True
                break
            elif sign == "converged":
                self.converged = True
                self.table.print("Operator indicated convergence!")
                break

            self.log("")
            yield self.cur_cycle
        else:
            self.table.print("Number of cycles exceeded!")

        #############################################
        # Outside big loop driving the optimization #
        #############################################

        self.print("")
        print("@ Finished optimization; doing final summary")
        if self.dump:
            self.out_trj_handle.close()

        if (not self.is_cos) and (not self.stopped):
            self.print(self.final_summary())
            # Remove 'current_geometry.xyz' file
            try:
                os.remove(self.current_fn)
            except FileNotFoundError:
                self.log(f"Tried to delete '{self.current_fn}'. Couldn't find it.")
        with open(self.final_fn, "w") as handle:
            handle.write(self.geometry.as_xyz())
        self.table.print(
            f"Wrote final, hopefully optimized, geometry to '{self.final_fn.name}'"
        )
        self.postprocess_opt()
        sys.stdout.flush()

    @property
    def run_generator(self):
        try:
            self._run_generator
        except AttributeError:
            self._run_generator = self.get_run_generator()
        return self._run_generator

    def run(self):
        while True:
            try:
                next(self.run_generator)
            except StopIteration:
                return

    def _get_opt_restart_info(self):
        """To be re-implemented in the derived classes."""
        return dict()

    def _set_opt_restart_info(self, opt_restart_info):
        """To be re-implemented in the derived classes."""
        return

    def get_restart_info(self):
        restart_info = {
            "geom_info": self.geometry.get_restart_info(),
            "last_cycle": self.cur_cycle,
            "max_cycles": self.max_cycles,
            "energies": self.energies,
            "coords": self.coords,
            "forces": [forces.tolist() for forces in self.forces],
            "steps": [step.tolist() for step in self.steps],
        }
        restart_info.update(self._get_opt_restart_info())
        return restart_info

    def set_restart_info(self, restart_info):
        # Set restart information general to all optimizers
        self.last_cycle = restart_info["last_cycle"] + 1

        must_resize = self.last_cycle >= self.max_cycles
        if must_resize:
            self.max_cycles += restart_info["max_cycles"]
            # Resize HDF5
            if self.dump:
                h5_group = get_h5_group(self.h5_fn, self.h5_group_name)
                resize_h5_group(h5_group, self.max_cycles)
                h5_group.file.close()

        self.coords = [np.array(coords) for coords in restart_info["coords"]]
        self.energies = restart_info["energies"]
        self.forces = [np.array(forces) for forces in restart_info["forces"]]
        self.steps = [np.array(step) for step in restart_info["steps"]]

        # Set subclass specific information
        self._set_opt_restart_info(restart_info)

        # Propagate restart information downwards to the geometry
        self.geometry.set_restart_info(restart_info["geom_info"])

    def dump_restart_info(self):
        restart_info = self.get_restart_info()

        restart_fn = f"restart_{self.cur_cycle:03d}.yaml"
        restart_yaml = yaml.dump(restart_info)
        self.write_to_out_dir(restart_fn, restart_yaml)
