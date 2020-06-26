#!/usr/bin/env python3

# [1] https://pubs.acs.org/doi/pdf/10.1021/ct200290m?rand=dcfwsf09
# [2] https://onlinelibrary.wiley.com/doi/epdf/10.1002/jcc.23481
# [3] https://onlinelibrary.wiley.com/doi/epdf/10.1002/tcr.201600043

import itertools as it

import autograd
import autograd.numpy as anp
import numpy as np

from pysisyphus.calculators.Calculator import Calculator
from pysisyphus.constants import AU2KJPERMOL
from pysisyphus.elem_data import COVALENT_RADII
from pysisyphus.Geometry import Geometry
from pysisyphus.helpers import complete_fragments
from pysisyphus.io.hdf5 import get_h5_group, resize_h5_group


def get_data_model(atoms, max_cycles):
    coord_size = 3 * len(atoms)
    _1d = (max_cycles, )
    _2d = (max_cycles, coord_size)

    data_model = {
        "cart_coords": _2d,
        "energy": _1d,
        "forces": _2d,
        "true_energy": _1d,
        "true_forces": _2d,
    }

    return data_model


def afir_closure(fragment_indices, cov_radii, gamma, rho=1, p=6):
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
    if gamma == 0.:
        alpha = 0.
    else:
        alpha = gamma / ((2**(-1/6) - (1 + (1 + gamma/epsilon)**0.5)**(-1/6)) * R0)
    
    def afir_func(coords3d):
        diffs = anp.diff(coords3d[inds], axis=1).reshape(-1, 3)
        rs = anp.linalg.norm(diffs, axis=1)

        omegas = (cov_rad_sums / rs)**p

        f = alpha * rho * (omegas*rs).sum() / omegas.sum()
        return f
    return afir_func


class AFIR(Calculator):

    def __init__(self, calculator, fragment_indices, gamma, rho=1, p=6,
                 dump=True, h5_fn="afir.h5", h5_group_name="afir", **kwargs):
        super().__init__(**kwargs)

        self.calculator = calculator
        self.fragment_indices = fragment_indices
        assert len(self.fragment_indices) in (1, 2)
        # gamma is expected to be given in kJ/mol. convert it to au.
        self.gamma = gamma / AU2KJPERMOL
        assert self.gamma > 0
        self.rho = int(rho)
        assert self.rho in (-1, 1)
        self.p = p
        self.dump = dump
        self.h5_fn = h5_fn
        self.h5_group_name = h5_group_name
        # We can initialize the HDF5 group as we don't know the shape of
        # atoms/coords yet. So we wait until after the first calculation.
        self.h5_group = None
        self.h5_cycles = 50

        rho_verbose = { 1: ("pushing", "together"),
                       -1: ("pulling", "apart")
        }
        w1, w2 = rho_verbose[self.rho]
        self.log(f"rho={self.rho}, {w1} framgents {w2}")

        self.atoms = None
        self.calc_counter = 0

    def init_h5_group(self, atoms, max_cycles=None):
        if max_cycles is None:
            max_cycles = self.h5_cycles
        self.data_model = get_data_model(atoms, max_cycles)
        self.h5_group = get_h5_group(self.h5_fn, self.h5_group_name, self.data_model,
                                     reset=True)

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
            fn = f"frag_geom_{i:02d}.xyz"
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
        self.fragment_indices = complete_fragments(self.atoms, self.fragment_indices)
        self.log_fragments()
        self.write_fragment_geoms(atoms, coords)
        self.cov_radii = np.array([COVALENT_RADII[atom.lower()] for atom in atoms]) 
        self.log("Set covalent radii")
        self.afir_func = afir_closure(self.fragment_indices,
                                      self.cov_radii,
                                      self.gamma,
                                      rho=self.rho,
                                      p=self.p)
        self.log("Created and set AFIR function.")
        self.afir_grad_func = autograd.grad(self.afir_func)
        self.log("Created and set AFIR gradient function.")

    def get_energy(self, atoms, coords):
        self.set_atoms_and_funcs(atoms, coords)

        true_energy = self.calculator.get_energy(atoms, coords)["energy"]
        afir_energy = self.afir_func(coords.reshape(-1, 3))
        self.log()

        results = {
            "energy": true_energy+afir_energy,
            "true_energy": true_energy,
        }
        if self.dump:
            self.dump_h5(atoms, coords, results)
        self.calc_counter += 1
        return results

    def get_forces(self, atoms, coords):
        self.set_atoms_and_funcs(atoms, coords)

        coords3d = coords.reshape(-1, 3)
        results = self.calculator.get_forces(atoms, coords)
        true_energy = results["energy"]
        true_forces = results["forces"]

        afir_energy = self.afir_func(coords3d)
        afir_forces = -self.afir_grad_func(coords3d).flatten()

        true_norm = np.linalg.norm(true_forces)
        afir_norm = np.linalg.norm(afir_forces)
        self.log(f"norm(true_forces)={true_norm:.6f} au/bohr")
        self.log(f"norm(afir_forces)={afir_norm:.6f} au/bohr")
        self.log()

        results = {
            "energy": true_energy+afir_energy,
            "forces": true_forces+afir_forces,
            "true_forces": true_forces,
            "true_energy": true_energy,
        }
        if self.dump:
            self.dump_h5(atoms, coords, results)
        self.calc_counter += 1
        return results
