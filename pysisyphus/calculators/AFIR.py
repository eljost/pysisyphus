# [1] https://pubs.acs.org/doi/pdf/10.1021/ct200290m?rand=dcfwsf09
# [2] https://onlinelibrary.wiley.com/doi/epdf/10.1002/jcc.23481
# [3] https://onlinelibrary.wiley.com/doi/epdf/10.1002/tcr.201600043

from dataclasses import dataclass
import itertools as it
from typing import List

import autograd
import autograd.numpy as anp
import numpy as np
import numpy.typing as npt

from pysisyphus.calculators.Calculator import Calculator
from pysisyphus.elem_data import COVALENT_RADII
from pysisyphus.Geometry import Geometry
from pysisyphus.helpers import complete_fragments
from pysisyphus.helpers_pure import log
from pysisyphus.io.hdf5 import get_h5_group, resize_h5_group
from pysisyphus.linalg import finite_difference_hessian


@dataclass
class AFIRPath:
    atoms: tuple
    cart_coords: np.ndarray
    energies: np.ndarray
    forces: np.ndarray
    opt_is_converged: bool


def get_data_model(atoms, max_cycles):
    coord_size = 3 * len(atoms)
    _1d = (max_cycles,)
    _2d = (max_cycles, coord_size)
    _3d = (max_cycles, coord_size, coord_size)

    data_model = {
        "cart_coords": _2d,
        "energy": _1d,
        "forces": _2d,
        "hessian": _3d,
        "true_energy": _1d,
        "true_forces": _2d,
        "true_hessian": _3d,
    }

    return data_model


def afir_closure(
    fragment_indices, cov_radii, gamma, rho=1, p=6, prefactor=1.0, logger=None
):
    """rho=1 pushes fragments together, rho=-1 pulls fragments apart."""

    # See https://onlinelibrary.wiley.com/doi/full/10.1002/qua.24757
    # Eq. (9) for extension to 3 fragments.
    assert len(fragment_indices) == 2

    inds = np.array(list(it.product(*fragment_indices)))
    cov_rad_sums = cov_radii[inds].sum(axis=1)

    # 3.8164 Angstrom in Bohr
    R0 = 7.21195
    # 1.0061 kJ/mol to Hartree
    epsilon = 0.000383203368

    # Avoid division by zero for gamma = 0.
    if gamma == 0.0:
        alpha = 0.0
    else:
        alpha = gamma / (
            (2 ** (-1 / 6) - (1 + (1 + gamma / epsilon) ** 0.5) ** (-1 / 6)) * R0
        )

    rho_verbose = {1: ("pushing", "together"), -1: ("pulling", "apart")}
    w1, w2 = rho_verbose[rho]
    log(
        logger,
        f"Creating AFIR closure with α={alpha:.6f}, prefactor {prefactor:.6f}, "
        f"rho={rho}, {w1} framgents {w2}",
    )

    def afir_func(coords3d):
        diffs = anp.diff(coords3d[inds], axis=1).reshape(-1, 3)
        rs = anp.linalg.norm(diffs, axis=1)

        omegas = (cov_rad_sums / rs) ** p

        f = prefactor * alpha * rho * (omegas * rs).sum() / omegas.sum()
        return f

    return afir_func


class AFIR(Calculator):
    def __init__(
        self,
        calculator: Calculator,
        fragment_indices: List[List[int]],
        gamma: npt.ArrayLike,
        rho: npt.ArrayLike = 1,
        p: int = 6,
        ignore_hydrogen: bool = True,
        complete_fragments: bool = True,
        dump: bool = True,
        h5_fn: str = "afir.h5",
        h5_group_name: str = "afir",
        **kwargs,
    ):
        """
        Artifical Force Induced Reaction calculator.

        Currently, there are no automated drivers to run large-scale AFIR calculations
        with many different initial orientations and/or increasing collision energy
        parameter γ. Nontheless, selected AFIR calculations can be carried out by hand.
        After convergence, artificial potential & forces, as well as real energies
        and forces can be plotted with 'pysisplot --afir'. The highest energy point
        along the AFIR path can then be selected for a subsequent TS-optimization,
        e.g. via 'pysistrj --get [index] optimzation.trj'.

        Future versions of pysisyphus may provide drivers for more automatted
        AFIR calculations.

        Parameters
        ----------
        calculator
            Actual QC calculator that provides energies and its derivatives,
            that are modified by the AFIR calculator, e.g., ORCA or Psi4.
        fragment_indices
            List of lists of integers, specifying the separate fragments. If
            the indices in theses lists don't comprise all atoms in the molecule,
            the reamining indices will be added as a separate fragment. If a AFIR
            calculation is carried out with 2 fragments and 'complete_fragments' is
            True (see below) it is enough to specify only
            the indices of one fragment, e.g., for a system of 10 atoms
            'fragment_indices=[[0,1,2,3]]' is enough. The second system will be set
            up automatically with indices [4,5,6,7,8,9].
        gamma
            Collision energy parameter γ in au. For 2 fragments it can be a single
            integer, while for > 2 fragments a list of gammas must be given, specifying
            the pair-wise collision energy parameters. For 3 fragments 3 gammas must be
            given [γ_01, γ_02, γ_12], for 4 fragments 6 gammas would be required
            [γ_01, γ_02, γ_03, γ_12, γ_13, γ_23] and so on.
        rho
            Direction of the artificial force, either 1 or -1. The same comments
            as for gamma apply. For 2 fragments a single integer is enough, for
            > 2 fragments a list of rhos must be given (see above).  For rho=1
            fragments are pushed together, for rho=-1 fragments are pulled apart.
        p
            Exponent p used in the calculation of the weight function ω. Defaults
            to 6 and probably does not have to be changed.
        ignore_hydrogen
            Whether hydrogens are ignored in the calculation of the artificial force.
        complete_fragments
            Whether an incomplete specification in 'fragment_indices' is automatically
            completed.
        dump
            Whether an HDF5 file is created.
        h5_fn
            Filename of the HDF5 file used for dumping.
        h5_group_name
            HDF5 group name used for dumping.

        Other Parameters
        ----------------
        **kwargs
            Keyword arguments passed to the Calculator baseclass.
        """
        super().__init__(**kwargs)

        self.calculator = calculator
        # Update own charge and multiplicity with values from the wrapped
        # calculator.
        self.charge = calculator.charge
        self.mult = calculator.mult
        self.fragment_indices = fragment_indices
        assert len(self.fragment_indices) > 0
        self.gamma = np.atleast_1d(gamma).astype(float)
        self.rho = np.atleast_1d(rho).astype(int)
        np.testing.assert_allclose(np.abs(self.rho), np.ones_like(self.rho))
        self.p = p
        self.ignore_hydrogen = ignore_hydrogen
        self.complete_fragments = complete_fragments
        self.dump = dump
        self.h5_fn = self.out_dir / h5_fn
        self.h5_group_name = h5_group_name
        # We can't initialize the HDF5 group as we don't know the shape of
        # atoms/coords yet. So we wait until after the first calculation.
        self.h5_group = None
        self.h5_cycles = 50

        if self.ignore_hydrogen:
            self.log("No artificial force contribution from hydrogens!")
        self.atoms = None
        self.calc_counter = 0

    """We try to keep charge and multiplicity consistent between the AFIR
    calculator and the actual wrapped calculator. But we will always return
    the charge and multiplicity of the underlying calculator."""

    @property
    def charge(self):
        return self.calculator.charge

    @charge.setter
    def charge(self, charge):
        try:
            self.calculator.charge = charge
        except AttributeError:
            pass
        self._charge = charge

    @property
    def mult(self):
        return self.calculator.mult

    @mult.setter
    def mult(self, mult):
        try:
            self.calculator.mult = mult
        except AttributeError:
            pass
        self._mult = mult

    def init_h5_group(self, atoms, max_cycles=None):
        if max_cycles is None:
            max_cycles = self.h5_cycles
        self.data_model = get_data_model(atoms, max_cycles)
        self.h5_group = get_h5_group(
            self.h5_fn, self.h5_group_name, self.data_model, reset=True
        )

    def dump_h5(self, atoms, coords, results):
        # Initialize if not done yet
        if self.h5_group is None:
            self.init_h5_group(atoms)
            # Write atoms once
            self.h5_group.attrs["atoms"] = atoms

        # Check if HDF5 datasets have to be resized
        cur_max_cycles = self.h5_group["cart_coords"].shape[0]
        need_resize = self.calc_counter > cur_max_cycles - 5
        if need_resize:
            new_max_cycles = cur_max_cycles + self.h5_cycles
            resize_h5_group(self.h5_group, max_cycles=new_max_cycles)

        for k, v in results.items():
            self.h5_group[k][self.calc_counter] = v
        self.h5_group["cart_coords"][self.calc_counter] = coords
        self.h5_group.attrs["cur_cycle"] = self.calc_counter

    def log_fragments(self):
        self.log(f"Using {len(self.fragment_indices)} fragments")
        for i, frag in enumerate(self.fragment_indices):
            self.log(f"Fragment {i:02d}, {len(frag)} atoms:")
            self.log(f"\t{frag}")

    def write_fragment_geoms(self, atoms, coords):
        geom = Geometry(atoms, coords)
        for i, frag in enumerate(self.fragment_indices):
            frag_geom = geom.get_subgeom(frag)
            fn = self.make_fn(f"frag_geom_{i:02d}.xyz")
            with open(fn, "w") as handle:
                handle.write(frag_geom.as_xyz())
            self.log(f"Wrote geometry of fragment {i:02d} to {fn}.")

    def set_atoms_and_funcs(self, atoms, coords):
        """Initially atoms was also an argument to the constructor of AFIR.
        I removed it so creation becomes easier.
        The first time a calculation is requested with a proper atom set
        everything is set up (cov. radii, afir function and corresponding
        gradient). Afterwards there is only a check if atoms != None and it
        is expected that all functions are properly set.

        fragment_indices can also be incomplete w.r.t. to the number of
        atoms. If the sum of the specified fragment atoms is less than the
        number of atoms present then all remaining unspecified atoms will
        be gathered in one fragment.
        """

        if self.atoms is not None:
            assert self.atoms == atoms
            return

        self.log("Setting atoms on AFIR calculator")
        self.atoms = atoms
        if self.complete_fragments:
            self.fragment_indices = complete_fragments(
                self.atoms, self.fragment_indices
            )
        self.log_fragments()
        self.write_fragment_geoms(atoms, coords)

        if self.ignore_hydrogen:
            hydrogen_inds = [i for i, atom in enumerate(atoms) if atom.lower() == "h"]
            fragment_indices = list()
            for i, frag_inds in enumerate(self.fragment_indices):
                frag_inds_no_h = [j for j in frag_inds if j not in hydrogen_inds]
                fragment_indices.append(frag_inds_no_h)
                dropped_hydrogens = len(frag_inds) - len(frag_inds_no_h)
                self.log(f"Ignoring {dropped_hydrogens} hydrogen(s) from fragment {i}.")
            self.fragment_indices = fragment_indices

        self.cov_radii = np.array([COVALENT_RADII[atom.lower()] for atom in atoms])
        self.log("Set covalent radii")
        self.afir_funcs = list()
        self.afir_grad_funcs = list()
        pairs = list(it.combinations(self.fragment_indices, 2))
        prefactor = 1 / len(pairs)
        self.log(
            f"Doing AFIR with {len(pairs)} fragment pairs and prefactor={prefactor:.4f}"
        )
        try:
            self.gamma = [float(self.gamma)] * len(pairs)
        except TypeError:
            assert len(self.gamma) == len(pairs)
        try:
            self.rho = [float(self.rho)] * len(pairs)
        except TypeError:
            assert len(self.rho) == len(pairs)
        self.log(f"Using gamma(s): {self.gamma}")
        self.log(f"  Using rho(s): {self.rho}")

        for (frag1, frag2), gamma, rho in zip(pairs, self.gamma, self.rho):
            afir_func = afir_closure(
                (frag1, frag2),
                self.cov_radii,
                gamma=gamma,
                rho=rho,
                p=self.p,
                prefactor=prefactor,
                logger=self.logger,
            )
            afir_grad_func = autograd.grad(afir_func)
            self.afir_funcs.append(afir_func)
            self.afir_grad_funcs.append(afir_grad_func)
        self.log("Created and set AFIR function & gradient function.")

    def get_energy(self, atoms, coords, **prepare_kwargs):
        self.set_atoms_and_funcs(atoms, coords)

        true_results = self.calculator.get_energy(atoms, coords, **prepare_kwargs)
        true_energy = true_results["energy"]
        coords3d = coords.reshape(-1, 3)
        # Iterate over all fragment pairs
        afir_energy = sum([afir_func(coords3d) for afir_func in self.afir_funcs])
        self.log()

        results = {
            "energy": true_energy + afir_energy,
            "true_energy": true_energy,
        }
        if self.dump:
            self.dump_h5(atoms, coords, results)
        self.calc_counter += 1
        return results

    def get_forces(self, atoms, coords, **prepare_kwargs):
        self.set_atoms_and_funcs(atoms, coords)

        true_results = self.calculator.get_forces(atoms, coords, **prepare_kwargs)
        true_energy = true_results["energy"]
        true_forces = true_results["forces"]

        coords3d = coords.reshape(-1, 3)
        afir_energy = 0.0
        afir_forces = np.zeros_like(coords)
        # Iterate over all fragment pairs
        for afir_func, afir_grad_func in zip(self.afir_funcs, self.afir_grad_funcs):
            afir_energy += afir_func(coords3d)
            # Add negative of the gradient (the force)
            afir_forces += -afir_grad_func(coords3d).flatten()

        true_norm = np.linalg.norm(true_forces)
        afir_norm = np.linalg.norm(afir_forces)
        self.log(
            f"\ntrue_energy={true_energy:.6f} au\n"
            f"afir_energy={afir_energy:.6f} au\n"
            f" sum_energy={true_energy+afir_energy:.6f} au\n"
            f"norm(true_forces)={true_norm:.6f} au/bohr\n"
            f"norm(afir_forces)={afir_norm:.6f} au/bohr\n"
        )

        results = {
            "energy": true_energy + afir_energy,
            "forces": true_forces + afir_forces,
            "true_forces": true_forces,
            "true_energy": true_energy,
        }
        if self.dump:
            self.dump_h5(atoms, coords, results)
        self.calc_counter += 1
        return results

    def afir_fd_hessian_wrapper(self, coords3d, afir_grad_func):
        coords = coords3d.flatten()

        def grad_func(coords):
            afir_grad = afir_grad_func(coords.reshape(-1, 3))
            return afir_grad.flatten()

        fd_hessian = finite_difference_hessian(coords, grad_func, acc=4)
        return fd_hessian

    def get_hessian(self, atoms, coords, **prepare_kwargs):
        self.set_atoms_and_funcs(atoms, coords)

        true_results = self.calculator.get_hessian(atoms, coords, **prepare_kwargs)
        true_energy = true_results["energy"]
        true_hessian = true_results["hessian"]

        coords3d = coords.reshape(-1, 3)
        afir_energy = 0.0
        afir_hessian = np.zeros_like(true_hessian)
        # Iterate over all fragment pairs
        for afir_func, afir_grad_func in zip(self.afir_funcs, self.afir_grad_funcs):
            afir_energy += afir_func(coords3d)
            # AFIR Hessian from finite differences
            afir_hessian += self.afir_fd_hessian_wrapper(coords3d, afir_grad_func)

        results = {
            "energy": true_energy + afir_energy,
            "hessian": true_hessian + afir_hessian,
            "true_hessian": true_hessian,
            "true_energy": true_energy,
        }
        if self.dump:
            self.dump_h5(atoms, coords, results)
        self.calc_counter += 1
        return results
