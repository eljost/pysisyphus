from collections import Counter, namedtuple
import copy
import itertools as it
from pathlib import Path
import re
import subprocess
import tempfile
import sys
import warnings

import h5py
import numpy as np
from scipy.spatial.distance import pdist
from scipy.spatial.transform import Rotation
import rmsd

try:
    from thermoanalysis.QCData import QCData
    from thermoanalysis.thermo import thermochemistry
except ModuleNotFoundError:
    pass

from pysisyphus.config import p_DEFAULT, T_DEFAULT
from pysisyphus.constants import BOHR2ANG
from pysisyphus.exceptions import DifferentAtomOrdering

from pysisyphus.wavefunction.excited_states import norm_ci_coeffs
from pysisyphus.hessian_proj import get_hessian_projector, inertia_tensor
from pysisyphus.linalg import are_collinear
from pysisyphus.elem_data import (
    ATOMIC_NUMBERS,
    COVALENT_RADII as CR,
    INV_ATOMIC_NUMBERS,
    ISOTOPE_DICT,
    MASS_DICT,
    VDW_RADII as VDWR,
)
from pysisyphus.helpers_pure import (
    eigval_to_wavenumber,
    full_expand,
    molecular_volume,
    to_subscript_num,
)
from pysisyphus.intcoords import (
    DLC,
    HDLC,
    RedundantCoords,
    TRIC,
    TMTRIC,
    HybridRedundantCoords,
    CartesianCoords,
    MWCartesianCoords,
)
from pysisyphus.intcoords.exceptions import (
    NeedNewInternalsException,
    RebuiltInternalsException,
    DifferentCoordLengthsException,
)
from pysisyphus.intcoords.helpers import get_tangent
from pysisyphus.intcoords.setup import BOND_FACTOR
from pysisyphus.intcoords.setup_fast import find_bonds
from pysisyphus.plot_ascii import plot_wrapper
from pysisyphus.xyzloader import make_xyz_str


def normalize_atoms(atoms) -> tuple[str]:
    atomic_numbers = set(INV_ATOMIC_NUMBERS.keys())
    _atoms = list()
    for atom in atoms:
        try:
            atom_int = int(atom)
            if atom_int in atomic_numbers:
                atom = INV_ATOMIC_NUMBERS[atom_int]
        except ValueError:
            pass
        # Was atom.capitalize() before ...
        atom = atom.lower()
        _atoms.append(atom)
    return tuple(_atoms)


class Geometry:
    coord_types = {
        "cart": None,
        "redund": RedundantCoords,
        "hredund": HybridRedundantCoords,
        "dlc": DLC,
        "hdlc": HDLC,
        "tric": TRIC,
        "tmtric": TMTRIC,
        "cartesian": CartesianCoords,
        "mwcartesian": MWCartesianCoords,
    }

    def __init__(
        self,
        atoms,
        coords,
        fragments=None,
        coord_type="cart",
        coord_kwargs=None,
        isotopes=None,
        freeze_atoms=None,
        remove_com=False,
        remove_centroid=False,
        comment="",
        name="",
    ):
        """Object representing atoms in a coordinate system.

        The Geometry represents atoms and their positions in coordinate
        system. By default cartesian coordinates are used, but internal
        coordinates are also possible.

        Parameters
        ----------
        atoms : iterable
            Iterable of length N, containing  element symbols.
        coords : 1d iterable
            1d iterable of length 3N, containing the cartesian coordinates
            of N atoms.
        fragments : dict, optional
            Dict with different keys denoting different fragments. The values
            contain lists of atom indices.
        coord_type : {"cart", "redund"}, optional
            Type of coordinate system to use. Right now cartesian (cart)
            and redundand (redund) are supported.
        coord_kwargs : dict, optional
            Dictionary containing additional arguments that get passed
            to the constructor of the internal coordinate class.
        isotopes : iterable of pairs, optional
            Iterable of pairs consisting of 0-based atom index and either an integer
            or a float. If an integer is given the closest isotope mass will be selected.
            Given a float, this float will be directly used as mass.
        freeze_atoms : iterable of integers
            Specifies which atoms should remain fixed at their initial positions.
        remove_com : bool, optional
            Move center of mass to the origin.
        remove_centroid : bool, optional
            Move centroid to the origin.
        comment : str, optional
            Comment string.
        name : str, optional
            Verbose name of the geometry, e.g. methanal or water. Used for printing
        """
        self.atoms = normalize_atoms(atoms)
        # self._coords always holds cartesian coordinates.
        self._coords = np.array(coords, dtype=float).flatten()
        assert self._coords.size == (3 * len(self.atoms)), (
            f"Expected 3N={3*len(self.atoms)} cartesian coordinates but got "
            f"{self._coords.size}. Did you accidentally supply internal "
            "coordinates?"
        )
        if fragments is None:
            fragments = dict()
        self.fragments = fragments
        self.coord_type = coord_type
        if coord_kwargs is None:
            coord_kwargs = dict()
        self.coord_kwargs = coord_kwargs
        if isotopes is None:
            isotopes = list()
        self.isotopes = isotopes
        if freeze_atoms is None:
            freeze_atoms = list()
        elif type(freeze_atoms) is str:
            freeze_atoms = full_expand(freeze_atoms)
        self.freeze_atoms = np.array(freeze_atoms, dtype=int)
        self.remove_com = bool(remove_com)
        self.remove_centroid = bool(remove_centroid)
        self.comment = comment
        self.name = name

        self._masses = None
        if self.remove_com:
            self.coords3d = self.coords3d - self.center_of_mass[None, :]
        if self.remove_centroid:
            self.coords3d = self.coords3d - self.centroid[None, :]

        self._energy = None
        self._forces = None
        self._hessian = None
        self._all_energies = None
        self.calculator = None

        assert (
            # Negative atom indices are not allowed.
            all(self.freeze_atoms >= 0)
            and (
                # Allow an empty array, no frozen atoms.
                (self.freeze_atoms.size == 0)
                # Or check that the biggest index is still in the valid range
                or (self.freeze_atoms.max() < len(self.atoms))
            )
        ), f"'freeze_atoms' must all be >= 0 and < {len(self.atoms)}!"

        # Disallow any coord_kwargs with coord_type == 'cart'
        if (coord_type == "cart") and not (coord_kwargs is None or coord_kwargs == {}):
            print(
                "coord_type is set to 'cart' but coord_kwargs were given. "
                "This is probably not intended. Exiting!"
            )
            sys.exit()

        # Coordinate systems are handled below
        coord_class = self.coord_types[self.coord_type]
        if coord_class:
            if (len(self.freeze_atoms) > 0) and ("freeze_atoms" not in coord_kwargs):
                coord_kwargs["freeze_atoms"] = freeze_atoms
            self.internal = coord_class(
                atoms,
                self.coords3d.copy(),
                masses=self.masses,
                **coord_kwargs,
            )
        else:
            self.internal = None

    @property
    def moving_atoms(self):
        return [atom for i, atom in enumerate(self.atoms) if i not in self.freeze_atoms]

    def moving_atoms_jmol(self):
        atoms = list()
        freeze_atoms = self.freeze_atoms
        for i, atom in enumerate(self.atoms):
            atom = atom if i not in freeze_atoms else "X"
            atoms.append(atom)
        self.jmol(atoms=atoms)

    @property
    def sum_formula(self):
        atoms = self.atoms
        atoms = [atom.capitalize() for atom in atoms]
        unique_atoms = sorted(set(atoms))
        counter = Counter(atoms)
        atoms = list()
        num_strs = list()

        def set_atom(atom):
            atoms.append(atom)
            num = counter[atom]
            if num == 1:
                num_str = ""
            else:
                num_str = to_subscript_num(num)
            num_strs.append(num_str)

        # Hill-System
        for atom in ("C", "H"):
            try:
                unique_atoms.remove(atom)
                set_atom(atom)
            except ValueError:
                pass
        for atom in unique_atoms:
            set_atom(atom)

        return "".join([f"{atom}{num_str}" for atom, num_str in zip(atoms, num_strs)])

    def assert_compatibility(self, other):
        """Assert that two Geometries can be substracted from each other.

        Parameters
        ----------
        other : Geometry
            Geometry for comparison.
        """
        same_atoms = self.atoms == other.atoms
        same_coord_type = self.coord_type == other.coord_type
        same_coord_length = len(self.coords) == len(other.coords)
        assert same_atoms, "Atom number/ordering is incompatible!"
        assert same_coord_type, "coord_types are incompatible!"
        try:
            assert same_coord_length, "Different length of coordinate vectors!"
        except AssertionError:
            raise DifferentCoordLengthsException

    def __eq__(self, other):
        return (self.atoms == other.atoms) and np.allclose(
            self.coords, other.coords, atol=1e-8
        )

    def __sub__(self, other):
        self.assert_compatibility(other)
        if self.coord_type in ("cart", "cartesian"):
            diff = self.coords - other.coords
        elif self.coord_type in ("redund", "dlc"):
            # Take periodicity of dihedrals into account by calling
            # get_tangent(). Care has to be taken regarding the orientation
            # of the returned tangent vector. It points from self to other.
            #
            # As we want to return the difference between two vectors we
            # have to reverse the direction of the tangent by multiplying it
            # with -1 to be consistent with basic subtraction laws ...
            # A - B = C, where C is a vector pointing from B to A (B + C = A)
            # In our case get_tangent returns B - A, that is a vector pointing
            # from A to B.
            diff = -get_tangent(
                self.internal.prim_coords,
                other.internal.prim_coords,
                self.internal.dihedral_indices,
            )
        else:
            raise Exception("Invalid coord_type!")

        # Convert to DLC
        if self.coord_type == "dlc":
            diff = self.internal.U.T.dot(diff)
        return diff

    def __add__(self, other):
        atoms = tuple(self.atoms) + tuple(other.atoms)
        coords = np.concatenate((self.cart_coords, other.cart_coords))
        return Geometry(atoms, coords)

    @property
    def is_linear(self):
        return are_collinear(self.coords3d)

    def atom_xyz_iter(self):
        return iter(zip(self.atoms, self.coords3d))

    def copy(self, coord_type=None, coord_kwargs=None):
        """Returns a new Geometry object with same atoms and coordinates.

        Parameters
        ----------
        coord_type : str
            Desired coord_type, defaults to current coord_type.

        coord_kwargs : dict, optional
            Any desired coord_kwargs that will be passed to the RedundantCoords
            object.
        Returns
        -------
        geom : Geometry
            New Geometry object with the same atoms and coordinates.
        """
        if coord_type is None:
            coord_type = self.coord_type

        if coord_kwargs is None:
            coord_kwargs = dict()

        # Geometry constructor will exit when coord_kwargs are given
        # with coord_type == 'cart'. So we only supply it when we are
        # NOT using cartesian coordinates.
        _coord_kwargs = None
        if coord_type != "cart":
            try:
                typed_prims = self.internal.typed_prims
            # Will be raised if the current coord_type is 'cart'
            except AttributeError:
                typed_prims = None
            _coord_kwargs = {
                "typed_prims": typed_prims,
                "check_bends": True,
            }
            _coord_kwargs.update(coord_kwargs)
        return Geometry(
            self.atoms,
            self._coords.copy(),
            coord_type=coord_type,
            coord_kwargs=_coord_kwargs,
            isotopes=copy.deepcopy(self.isotopes),
            freeze_atoms=self.freeze_atoms.copy(),
        )

    def copy_all(self, coord_type=None, coord_kwargs=None):
        new_geom = self.copy(coord_type, coord_kwargs)
        new_geom.set_calculator(self.calculator)
        new_geom.energy = self._energy
        if self._forces is not None:
            new_geom.cart_forces = self._forces
        if self._hessian is not None:
            new_geom.cart_hessian = self._hessian
        return new_geom

    def atom_indices(self):
        """Dict with atom types as key and corresponding indices as values.

        Returns
        -------
        inds_dict : dict
            Unique atom types as keys, corresponding indices as values.
        """
        inds_dict = {}
        for atom_type in set(self.atoms):
            inds_dict[atom_type] = [
                i for i, atom in enumerate(self.atoms) if atom == atom_type
            ]
        return inds_dict

    @property
    def atom_types(self):
        return set(self.atoms)

    @property
    def atomic_numbers(self):
        return [ATOMIC_NUMBERS[a.lower()] for a in self.atoms]

    def get_fragments(self, regex):
        regex = re.compile(regex)
        frags = [frag for frag in self.fragments.keys() if regex.search(frag)]
        org_indices = list(it.chain(*[self.fragments[frag] for frag in frags]))

        new_atoms = [self.atoms[ind] for ind in org_indices]
        new_coords = self.coords3d[org_indices].copy()
        new_fragments = dict()
        i = 0
        for frag in frags:
            frag_atoms = len(self.fragments[frag])
            new_fragments[frag] = list(range(i, i + frag_atoms))
            i += frag_atoms
        return Geometry(new_atoms, new_coords, fragments=new_fragments)

    @property
    def layers(self):
        try:
            layers = self.calculator.layers
        except AttributeError:
            layers = (None,)
        return layers

    def del_atoms(self, inds, **kwargs):
        atoms = [atom for i, atom in enumerate(self.atoms) if i not in inds]
        c3d = self.coords3d
        coords3d = np.array([c3d[i] for i, _ in enumerate(self.atoms) if i not in inds])
        return Geometry(atoms, coords3d.flatten(), **kwargs)

    def set_calculator(self, calculator, clear=True):
        """Reset the object and set a calculator."""
        if clear:
            self.clear()
        self.calculator = calculator

    @property
    def is_analytical_2d(self):
        try:
            return self.calculator.analytical_2d
        except AttributeError:
            return False

    @property
    def mm_inv(self):
        """Inverted mass matrix.

        Returns a diagonal matrix containing the inverted atomic
        masses.
        """
        return np.diag(1 / self.masses_rep)

    @property
    def mm_sqrt_inv(self):
        """Inverted square root of the mass matrix."""
        return np.diag(1 / (self.masses_rep**0.5))

    @property
    def coords(self):
        """1d vector of atomic coordinates.

        Returns
        -------
        coords : np.array
            1d array holding the current coordinates.
        """
        if self.internal:
            coords = self.internal.coords
        else:
            # self._coords will always hold Cartesian coordinates.
            coords = self._coords
        return coords

    def set_coord(self, ind, coord):
        """Set a coordinate by index.

        Parameters
        ----------
        ind : int
            Index in of the coordinate to set in the self.coords array.
        coord : float
            Coordinate value.
        """
        assert (
            self.coord_type == "cart" and len(self.freeze_atoms) == 0
        ), "set_coord was not yet tested with coord_type != 'cart' and frozen atoms!"
        self.coords[ind] = coord
        self.clear()

    def set_coords(self, coords, cartesian=False, update_constraints=False):
        coords = np.array(coords).flatten()

        # Do Internal->Cartesian backtransformation if internal coordinates are used.
        if hasattr(self, "internal") and self.internal:
            # When internal coordinates are employed it may happen, that the underlying
            # Cartesian coordinates are updated, e.g. from the IPIServer calculator, which
            # may yield different internal coordinates.
            #
            # Here we update the Cartesians of the internal coordinate object to the new
            # values and calculate new internal coordinates, from which we can derive a step
            # in internals.
            if cartesian:
                self.assert_cart_coords(coords)
                cart_coords = coords.copy()
                # Update Cartesians of internal coordinate object and calculate
                # new internals.
                self.internal.coords3d = coords
                # Determine new internal coordinates, so we can later calculate a
                # step in internal coordinates.
                coords = self.internal.coords
                # Finally we also update the Cartesian coordinates of the Geometry object,
                # so the subsequent sanity check does not fail. This also allows updating
                # the coordiantes of atoms that are frozen. We set Geometry._coords directly,
                # instead of Geometry.cart_coords or Geometry.coords3d, to avoid an infinite
                # recursion.
                self._coords = cart_coords

            # Sanity check, asserting that the cartesian coordinates of the
            # Geometry object and the internal coordinate object are the same.
            np.testing.assert_allclose(self.coords3d, self.internal.coords3d)

            try:
                int_step = coords - self.internal.coords
                cart_step = self.internal.transform_int_step(
                    int_step, update_constraints=update_constraints
                )
                # From now on coords will always hold Cartesian coordinates!
                coords = self._coords + cart_step
            except NeedNewInternalsException as exception:
                invalid_inds = exception.invalid_inds
                # Check if the remaining internal coordinates are valid
                valid_typed_prims = [
                    typed_prim
                    for i, typed_prim in enumerate(self.internal.typed_prims)
                    if i not in invalid_inds
                ]
                coords3d = exception.coords3d.copy()
                coord_class = self.coord_types[self.coord_type]
                coord_kwargs = self.coord_kwargs.copy()
                """Instead of using only the remaining, valid typed_prims
                we could look for an entirely new set of typed_prims.
                
                But when we do this and we end up with more coordinates
                than before, this will lead to problems with the HDF5 dump.
                
                No problems arise when fewer coordinates are used
                (valid_typed_prims <= self.internal.typed_prims).
                With typed prims, only the remaining, valid typed_prims
                will be defined for the new geometry.
                
                coord_kwargs["typed_prims"] = valid_typed_prims # Currently disabled
                
                With 'define_prims' the remaining, valid typed_prims
                will be used, together with newly determined internal
                coordinates. This supports, e.g., the switch from a simple
                bend to a linear bend and its complement.
                
                Currently the default."""
                coord_kwargs.update(
                    {
                        "define_prims": valid_typed_prims,
                        "constrain_prims": self.internal.constrain_prims,
                    }
                )

                self.internal = coord_class(self.atoms, coords3d, **coord_kwargs)
                self._coords = coords3d.flatten()
                raise RebuiltInternalsException(
                    typed_prims=self.internal.typed_prims.copy()
                )

        # Restore original coordinates of frozen atoms. Right now this should
        # be redundant, as the Cartesian step is also constrainted in the
        # Internal->Cartesian backtransformation. But we keep it for now.
        coords.reshape(-1, 3)[self.freeze_atoms] = self.coords3d[self.freeze_atoms]
        # Set new Cartesian coordinates
        self._coords = coords
        # Reset all values because no calculations with the new coords
        # have been performed yet.
        self.clear()

    def reset_coords(self, new_typed_prims=None):
        if self.coord_type == "cart":
            return

        coord_class = self.coord_types[self.coord_type]
        self.internal = coord_class(
            self.atoms, self.coords3d, typed_prims=new_typed_prims
        )

    @coords.setter
    def coords(self, coords):
        """Wrapper for saving coordinates internally.

        Parameters
        ----------
        coords : np.array
            1d array containing atomic coordiantes. It's length
            depends on the coordinate system.
        """
        self.set_coords(coords)

    @property
    def coords3d(self):
        """Coordinates in 3d.

        Returns
        -------
        coords3d : np.array
            Coordinates of the Geometry as 2D array.
        """
        return self._coords.reshape(-1, 3)

    @coords3d.setter
    def coords3d(self, coords3d):
        self.set_coords(coords3d, cartesian=True)

    @property
    def cart_coords(self):
        return self._coords

    @cart_coords.setter
    def cart_coords(self, coords):
        self.set_coords(coords, cartesian=True)

    @property
    def coords_by_type(self):
        """Coordinates in 3d by atom type and their corresponding indices.

        Returns
        -------
        cbt : dict
            Dictionary with the unique atom types of the Geometry as keys.
            It's values are the 3d coordinates of the corresponding atom type.
        inds : dict
            Dictionary with the unique atom types of the Geometry as keys.
            It's values are the original indices of the 3d coordinates in the
            whole coords3d array.
        """
        cbt = dict()
        inds = dict()
        # for i, (atom, c3d) in enumerate(zip(self.atoms, self.coords3d)):
        # cbt.setdefault(atom, list()).append((i, c3d.tolist()))
        for i, (atom, c3d) in enumerate(zip(self.atoms, self.coords3d)):
            cbt.setdefault(atom, list()).append((c3d))
            inds.setdefault(atom, list()).append(i)
        for atom, c3d in cbt.items():
            cbt[atom] = np.array(c3d)
            inds[atom] = np.array(inds[atom])
        return cbt, inds

    @property
    def comment(self):
        en_width = 20
        # Check if we have to drop an (old) energy entry
        try:
            _ = float(self._comment[:en_width])
            # Drop old energy entry
            self._comment = self._comment[en_width + 2 :]
        except (ValueError, IndexError):
            pass

        # Prepend (new) energy, if present
        if self._energy:
            en_str = f"{self._energy: >{en_width}.8f} , "
        else:
            en_str = ""
        return f"{en_str}{self._comment}"

    @comment.setter
    def comment(self, new_comment):
        self._comment = new_comment

    @property
    def masses(self) -> np.ndarray:
        if self._masses is None:
            # Lookup tabuled masses in internal database
            masses = np.array([MASS_DICT[atom.lower()] for atom in self.atoms])
            # Use (different) isotope masses if requested
            for atom_index, iso_mass in self.isotopes:
                if "." not in str(iso_mass):
                    atom = self.atoms[atom_index].lower()
                    key = (atom, iso_mass)
                    try:
                        iso_mass = ISOTOPE_DICT[key]
                    except KeyError as err:
                        print(
                            f"Found no suitable mass for '{atom.capitalize()}' with approx. "
                            f"mass of ~{iso_mass} au!"
                        )
                        raise err
                masses[atom_index] = float(iso_mass)
            self.masses = masses
        return self._masses

    @masses.setter
    def masses(self, masses):
        assert len(masses) == len(self.atoms)
        masses = np.array(masses, dtype=float)
        self._masses = masses
        # Also try to propagate updated masses to the internal coordiante object
        try:
            self.internal.masses = masses
        except AttributeError:
            pass

    @property
    def masses_rep(self):
        # Some of the analytical potentials are only 2D
        repeat_masses = 2 if (self._coords.size == 2) else 3
        return np.repeat(self.masses, repeat_masses)

    @property
    def total_mass(self):
        return sum(self.masses)

    def center_of_mass_at(self, coords3d):
        """Returns the center of mass at given coords3d.

        Parameters
        ----------
        coords3d : np.array, shape(N, 3)
            Cartesian coordiantes.

        Returns
        -------
        R : np.array, shape(3, )
            Center of mass.
        """
        return 1 / self.total_mass * np.sum(coords3d * self.masses[:, None], axis=0)

    @property
    def center_of_mass(self):
        """Returns the center of mass.

        Returns
        -------
        R : np.array, shape(3, )
            Center of mass.
        """
        return self.center_of_mass_at(self.coords3d)

    @property
    def centroid(self):
        """Geometric center of the Geometry.

        Returns
        -------
        R : np.array, shape(3, )
            Geometric center of the Geometry.
        """
        return self.coords3d.mean(axis=0)

    def center(self):
        self.coords3d -= self.centroid[None, :]

    @property
    def mw_coords(self):
        """Mass-weighted coordinates.

        Returns
        -------
        mw_coords : np.array
            1d array containing the mass-weighted cartesian coordiantes.
        """
        return np.sqrt(self.masses_rep) * self._coords

    @mw_coords.setter
    def mw_coords(self, mw_coords):
        """Set mass-weighted coordinates."""
        self.coords = mw_coords / np.sqrt(self.masses_rep)

    def fd_coords3d_gen(self, step_size=1e-3):
        """Iterator returning 3d Cartesians for finite-differences."""
        coords3d = self.coords3d.copy()
        zeros = np.zeros_like(coords3d)
        for i, _ in enumerate(self.coords3d):
            for j in (0, 1, 2):
                step = zeros.copy()
                step[i, j] = step_size
                yield i, j, coords3d + step, coords3d - step

    @property
    def covalent_radii(self):
        return np.array([CR[a.lower()] for a in self.atoms])

    @property
    def vdw_radii(self):
        return np.array([VDWR[a.lower()] for a in self.atoms])

    def vdw_volume(self, **kwargs):
        V_au, *_ = molecular_volume(self.coords3d, self.vdw_radii, **kwargs)
        return V_au

    @property
    def inertia_tensor(self):
        return inertia_tensor(self.coords3d, self.masses)

    def principal_axes_are_aligned(self):
        """Check if the principal axes are aligned with the cartesian axes.

        Returns
        -------
        aligned : bool
            Wether the principal axes are aligned or not.
        """
        w, v = np.linalg.eigh(self.inertia_tensor)
        return np.allclose(v, np.eye(3)), v

    def align_principal_axes(self):
        """Align the principal axes to the cartesian axes.

        https://math.stackexchange.com/questions/145023
        """
        I = self.inertia_tensor
        w, v = np.linalg.eigh(I)
        # rot = np.linalg.solve(v, np.eye(3))
        # self.coords3d = rot.dot(self.coords3d.T).T
        self.coords3d = v.T.dot(self.coords3d.T).T

    def standard_orientation(self):
        # Translate center of mass to cartesian origin
        self.coords3d -= self.center_of_mass
        # Try to rotate the principal axes onto the cartesian axes
        for i in range(5):
            self.align_principal_axes()
            aligned, vecs = self.principal_axes_are_aligned()
            if aligned:
                break

    def reparametrize(self, energy, forces):
        # Currently, self.calculator.get_coords is only implemented by the
        # IPServer, but it is deactivated there.
        try:
            # TODO: allow skipping the update
            results = self.calculator.get_coords(self.atoms, self.cart_coords)
            self.set_coords(results["coords"], cartesian=True)
            reparametrized = True
        except AttributeError:
            reparametrized = False
        return reparametrized

    @property
    def energy(self):
        """Energy of the current atomic configuration.

        Returns
        -------
        energy : float
            Energy of the current atomic configuration.
        """
        if self._energy is None:
            results = self.calculator.get_energy(self.atoms, self._coords)
            self.set_results(results)
        return self._energy

    @energy.setter
    def energy(self, energy):
        """Internal wrapper for setting the energy.

        Parameters
        ----------
        energy : float
        """
        self._energy = energy

    @property
    def has_energy(self):
        return self._energy is not None

    @property
    def all_energies(self):
        """Return energies of all states that were calculated.

        This will also set self.energy, which may NOT be the ground state,
        but the state correspondig to the 'root' attribute of the calculator."""
        if self._all_energies is None:
            results = self.calculator.get_all_energies(self.atoms, self._coords)
            self.set_results(results)
        return self._all_energies

    def get_root_energy(self, root):
        return self.all_energies[root]

    def has_all_energies(self):
        return self._all_energies is not None

    @all_energies.setter
    def all_energies(self, all_energies):
        """Internal wrapper for setting all energies.

        Parameters
        ----------
        all_energies : np.array
        """
        self._all_energies = all_energies

    @property
    def cart_forces(self):
        if self._forces is None:
            results = self.calculator.get_forces(self.atoms, self._coords)
            self.set_results(results)
        return self._forces

    @cart_forces.setter
    def cart_forces(self, cart_forces):
        cart_forces = np.array(cart_forces)
        assert cart_forces.shape == self.cart_coords.shape
        self._forces = cart_forces

    @property
    def forces(self):
        """Energy of the current atomic configuration.

        Returns
        -------
        force : np.array
            1d array containing the forces acting on the atoms. Negative
            of the gradient.
        """
        forces = self.cart_forces
        if self.internal:
            forces = self.internal.transform_forces(forces)
        return forces

    @forces.setter
    def forces(self, forces):
        """Internal wrapper for setting the forces.

        Parameters
        ----------
        forces : np.array
        """
        forces = np.array(forces)
        assert forces.shape == self.cart_coords.shape
        self._forces = forces

    @property
    def has_forces(self):
        return self._forces is not None

    @property
    def cart_gradient(self):
        return -self.cart_forces

    @cart_gradient.setter
    def cart_gradient(self, cart_gradient):
        self.cart_forces = -cart_gradient

    @property
    def gradient(self):
        """Negative of the force.

        Returns
        -------
        gradient : np.array
            1d array containing the negative of the current forces.
        """
        return -self.forces

    # @gradient.setter
    # def gradient(self, gradient):
    # """Internal wrapper for setting the gradient."""
    # # No check here as this is handled by in the forces.setter.
    # self.forces = -gradient

    @property
    def mw_gradient(self):
        """Mass-weighted gradient.

        Returns
        -------
        mw_gradient : np.array
            Returns the mass-weighted gradient.
        """
        return -self.forces / np.sqrt(self.masses_rep)

    @property
    def cart_hessian(self):
        if self._hessian is None:
            results = self.calculator.get_hessian(self.atoms, self._coords)
            self.set_results(results)
        return self._hessian

    @cart_hessian.setter
    def cart_hessian(self, cart_hessian):
        if cart_hessian is not None:
            cart_hessian = np.array(cart_hessian)
            assert cart_hessian.shape == (self.cart_coords.size, self.cart_coords.size)
        self._hessian = cart_hessian

    @property
    def hessian(self):
        """Matrix of second derivatives of the energy in respect to atomic
        displacements.

        Returns
        -------
        hessian : np.array
            2d array containing the second derivatives of the energy with respect
            to atomic/coordinate displacements depending on the type of
            coordiante system.
        """
        hessian = self.cart_hessian
        if self.internal:
            int_gradient = self.gradient
            return self.internal.transform_hessian(hessian, int_gradient)
        return hessian

    def has_hessian(self):
        return self._hessian is not None

    # @hessian.setter
    # def hessian(self, hessian):
    # """Internal wrapper for setting the hessian."""
    # assert hessian.shape == (self.coords.size, self.coords.size)
    # self._hessian = hessian

    def mass_weigh_hessian(self, hessian):
        return self.mm_sqrt_inv.dot(hessian).dot(self.mm_sqrt_inv)

    @property
    def mw_hessian(self):
        """Mass-weighted hessian.

        Returns
        -------
        mw_hessian : np.array
            2d array containing the mass-weighted hessian M^(-1/2) H M^(-1/2).
        """
        # M^(-1/2) H M^(-1/2)
        # TODO: Do the right thing here when the hessian is not yet calculated.
        #       this would probably involve figuring out how to mass-weigh and
        #       internal coordinat hessian... I think this is described in one
        #       of the Gonzalez-Schlegel-papers about the GS2 algorithm.
        return self.mass_weigh_hessian(self.cart_hessian)

    def unweight_mw_hessian(self, mw_hessian):
        """Unweight a mass-weighted hessian.

        Parameters
        ----------
        mw_hessian : np.array
            Mass-weighted hessian to be unweighted.

        Returns
        -------
        hessian : np.array
            2d array containing the hessian.
        """
        mm_sqrt = np.diag(self.masses_rep**0.5)
        return mm_sqrt.dot(mw_hessian).dot(mm_sqrt)

    def set_h5_hessian(self, fn):
        with h5py.File(fn, "r") as handle:
            atoms = handle.attrs["atoms"]
            hessian = handle["hessian"][:]

        # Also check lengths, as zip would lead to trunction for
        # different lenghts of self.atoms and atoms.
        valid = (len(atoms) == len(self.atoms)) and all(
            [ga.lower() == a.lower() for ga, a in zip(self.atoms, atoms)]
        )
        if valid:
            self.cart_hessian = hessian

    def get_normal_modes(
        self, cart_hessian=None, cart_gradient=None, proj_gradient=False, full=False
    ):
        """Normal mode wavenumbers, eigenvalues and Cartesian displacements Hessian."""
        if cart_hessian is None:
            cart_hessian = self.cart_hessian
        if not proj_gradient:
            mw_gradient = None
        elif cart_gradient is None:
            mw_gradient = self.mw_gradient
        else:
            mw_gradient = cart_gradient / np.sqrt(self.masses_rep)

        mw_hessian = self.mass_weigh_hessian(cart_hessian)
        proj_hessian, P = self.eckart_projection(
            mw_hessian, return_P=True, mw_gradient=mw_gradient, full=full
        )
        eigvals, eigvecs = np.linalg.eigh(proj_hessian)
        mw_cart_displs = P.T.dot(eigvecs)
        cart_displs = self.mm_sqrt_inv.dot(mw_cart_displs)
        cart_displs /= np.linalg.norm(cart_displs, axis=0)
        nus = eigval_to_wavenumber(eigvals)
        return nus, eigvals, mw_cart_displs, cart_displs

    def get_imag_frequencies(self, hessian=None, thresh=1e-6):
        vibfreqs, eigvals, *_ = self.get_normal_modes(hessian)
        return vibfreqs[eigvals < thresh]

    def get_thermoanalysis(
        self, energy=None, cart_hessian=None, T=T_DEFAULT, p=p_DEFAULT, point_group="c1"
    ):
        if cart_hessian is None:
            cart_hessian = self.cart_hessian
            # Delte any supplied energy value when a Hessian calculation is carried out
            energy = None

        if energy is None:
            energy = self.energy

        vibfreqs, *_ = self.get_normal_modes(cart_hessian)
        try:
            mult = self.calculator.mult
        except AttributeError:
            mult = 1
            warnings.warn(
                "Multiplicity for electronic entropy could not be determined! "
                f"Falling back to Using 2S+1 = {mult}."
            )

        thermo_dict = {
            "masses": self.masses,
            "wavenumbers": vibfreqs,
            "coords3d": self.coords3d,
            "scf_energy": energy,
            "mult": mult,
        }

        qcd = QCData(thermo_dict, point_group=point_group)
        thermo = thermochemistry(
            qcd, temperature=T, pressure=p, invert_imags=-15.0, cutoff=25.0
        )

        return thermo

    def get_trans_rot_projector(self, full=False):
        warnings.warn(
            "'Geometry.get_trans_rot_projector()' is deprecated. Please use "
            "'Geometry.get_hessian_projector() instead.",
            DeprecationWarning,
        )
        return get_hessian_projector(self.cart_coords, masses=self.masses, full=full)

    def get_hessian_projector(self, cart_gradient=None, full=False):
        return get_hessian_projector(
            self.cart_coords, masses=self.masses, cart_gradient=cart_gradient, full=full
        )

    def eckart_projection(
        self, mw_hessian, return_P=False, mw_gradient=None, full=False
    ):
        # Must not project analytical 2d potentials.
        if self.is_analytical_2d:
            return mw_hessian

        if mw_gradient is not None:
            cart_gradient = mw_gradient * np.sqrt(self.masses_rep)
        else:
            cart_gradient = None
        P = self.get_hessian_projector(cart_gradient=cart_gradient, full=full)
        proj_hessian = P.dot(mw_hessian).dot(P.T)
        # Projection seems to slightly break symmetry (sometimes?). Resymmetrize.
        proj_hessian = (proj_hessian + proj_hessian.T) / 2
        if return_P:
            return proj_hessian, P
        else:
            return proj_hessian

    def calc_energy(self):
        """Force energy calculation at the current coordinates."""
        results = self.calculator.get_energy(self.atoms, self.cart_coords)
        self.set_results(results)

    def calc_energy_and_forces(self):
        """Force energy and forces calculation at the current coordinates."""
        results = self.calculator.get_forces(self.atoms, self.cart_coords)
        self.set_results(results)

    def assert_cart_coords(self, coords):
        assert coords.size == self.cart_coords.size, (
            "This method only works with cartesian coordinate input. "
            "Did you accidentally provide internal coordinates?"
        )

    def get_temporary_coords(self, coords):
        if self.coord_type != "cart":
            int_step = coords - self.internal.coords
            cart_step = self.internal.transform_int_step(int_step, pure=True)
            coords = self.cart_coords + cart_step
        self.assert_cart_coords(coords)
        return coords

    def get_energy_at(self, coords):
        coords = self.get_temporary_coords(coords)
        return self.calculator.get_energy(self.atoms, coords)["energy"]

    def get_energy_at_cart_coords(self, cart_coords):
        self.assert_cart_coords(cart_coords)
        return self.calculator.get_energy(self.atoms, cart_coords)["energy"]

    def get_energy_and_forces_at(self, coords):
        """Calculate forces and energies at the given coordinates.

        The results are not saved in the Geometry object."""
        coords = self.get_temporary_coords(coords)
        results = self.calculator.get_forces(self.atoms, coords)
        self.zero_frozen_forces(results["forces"])

        if self.coord_type != "cart":
            results["forces"] = self.internal.transform_forces(results["forces"])

        return results

    def get_energy_and_cart_forces_at(self, cart_coords):
        self.assert_cart_coords(cart_coords)
        results = self.calculator.get_forces(self.atoms, cart_coords)
        self.zero_frozen_forces(results["forces"])
        return results

    def get_energy_and_cart_hessian_at(self, cart_coords):
        self.assert_cart_coords(cart_coords)
        results = self.calculator.get_hessian(self.atoms, cart_coords)
        return results

    def calc_double_ao_overlap(self, geom2):
        return self.calculator.run_double_mol_calculation(
            self.atoms, self.coords, geom2.coords
        )

    def zero_frozen_forces(self, cart_forces):
        cart_forces.reshape(-1, 3)[self.freeze_atoms] = 0.0

    def calc_wavefunction(self, **prepare_kwargs):
        # TODO: support wf (kw)-args?
        results = self.calculator.get_wavefunction(
            self.atoms, self.cart_coords, **prepare_kwargs
        )
        self.set_results(results)
        return results

    @property
    def wavefunction(self):
        if self._wavefunction is None:
            self.calc_wavefunction()
        return self._wavefunction

    @wavefunction.setter
    def wavefunction(self, wavefunction):
        self._wavefunction = wavefunction

    @property
    def td_1tdms(self):
        """1-particle transition density matrices from TD-DFT/TDA.

        Returns list of Xa, Ya, Xb and Yb in MO basis."""
        if self._td_1tdms is None:
            self.all_energies
        return self._td_1tdms

    @td_1tdms.setter
    def td_1tdms(self, td_1tdms):
        self._td_1tdms = td_1tdms

    def calc_relaxed_density(self, root, **prepare_kwargs):
        """Calculate a relaxed excited state density via an ES gradient calculation.

        The question is, if this method should set the wavefunction property
        at the current Geometry. On one hand, staying in pure python w/o numba
        the wavefunction sanity-check can become costly, even though it shouldn't
        be.
        On the other hand, setting the wavefunction would ensure consistency
        between the levels of theory used for density and wavefunction.

        For now, calculating an ES density does not set a wavefunction on the
        Geometry, whereas requesting the relaxed density for the GS does.

        TODO: add flag that allows setting the wavefunction (WF). Then,
        calculators should also include the WF in their results."""
        if root == 0:
            results = self.calc_wavefunction(**prepare_kwargs)
            wf = results["wavefunction"]
            density = wf.get_relaxed_density(0)
        else:
            results = self.calculator.get_relaxed_density(
                self.atoms, self.cart_coords, root, **prepare_kwargs
            )
            # Don't set density on Geometry
            density = results.pop("density")
            self.set_results(results)
        results["density"] = density
        if self.internal:
            results["forces"] = self.internal.transform_forces(results["forces"])
        return results

    def clear(self):
        """Reset the object state."""

        self._energy = None
        self._forces = None
        self._hessian = None
        self.true_energy = None
        self.true_forces = None
        self.true_hessian = None
        self._all_energies = None
        self._wavefunction = None
        self._td_1tdms = None

    def set_results(self, results):
        """Save the results from a dictionary.

        Parameters
        ----------
        results : dict
            The keys in this dict will be set as attributes in the current
            object, with the corresponding item as value.
        """

        trans = {
            "energy": "energy",
            "forces": "cart_forces",
            "hessian": "cart_hessian",
            # True properties in AFIR calculations
            "true_energy": "true_energy",
            "true_forces": "true_forces",
            "true_hessian": "true_hessian",
            # Overlap calculator; includes excited states
            "all_energies": "all_energies",
            "td_1tdms": "td_1tdms",
            # Wavefunction related
            "wavefunction": "wavefunction",
        }

        for key, value in results.items():
            # Zero forces of frozen atoms
            if key == "forces":
                self.zero_frozen_forces(value)
            elif key == "td_1tdms":
                value = norm_ci_coeffs(*value)
                results[key] = value

            setattr(self, trans[key], value)
        self.results = results

    def as_xyz(self, comment="", atoms=None, cart_coords=None):
        """Current geometry as a string in XYZ-format.

        Parameters
        ----------
        comment : str, optional
            Will be written in the second line (comment line) of the
            XYZ-string.
        cart_coords : np.array, 1d, shape (3 * atoms.size, )
            Cartesians for dumping instead of self._coords.

        Returns
        -------
        xyz_str : str
            Current geometry as string in XYZ-format.
        """
        if atoms is None:
            atoms = self.atoms
        if cart_coords is None:
            cart_coords = self._coords
        cart_coords = cart_coords.copy()
        cart_coords *= BOHR2ANG
        if comment is None:
            comment = ""
        elif comment == "":
            comment = self.comment
        return make_xyz_str(atoms, cart_coords.reshape((-1, 3)), comment)

    def dump_xyz(self, fn, cart_coords=None, **kwargs):
        fn = str(fn)
        if not fn.lower().endswith(".xyz"):
            fn = fn + ".xyz"
        with open(fn, "w") as handle:
            handle.write(self.as_xyz(cart_coords=cart_coords, **kwargs))

    def dump_trj(self, fn, trj_cart_coords, **kwargs):
        fn = Path(fn).with_suffix(".trj")
        xyzs = list()
        for cart_coords in trj_cart_coords:
            xyz = self.as_xyz(cart_coords=cart_coords, **kwargs)
            xyzs.append(xyz)
        trj = "\n".join(xyzs)
        with open(fn, "w") as handle:
            handle.write(trj)
        return fn

    def get_subgeom(self, indices, coord_type="cart", sort=False, cart_coords=None):
        """Return a Geometry containing a subset of the current Geometry.

        Parameters
        ----------
        indices : iterable of ints
            Atomic indices that the define the subset of the current Geometry.
        coord_type : str, ("cart", "redund"), optional
            Coordinate system of the new Geometry.
        cart_coords
            Optional 1d array of Cartesian coordinates of shape (3*natoms, ).

        Returns
        -------
        sub_geom : Geometry
            Subset of the current Geometry.
        """
        if cart_coords is not None:
            coords3d = cart_coords.reshape(-1, 3)
        else:
            coords3d = self.coords3d
        if sort:
            indices = sorted(indices)
        ind_list = list(indices)
        sub_atoms = [self.atoms[i] for i in ind_list]
        sub_coords = coords3d[ind_list]
        sub_geom = Geometry(sub_atoms, sub_coords.flatten(), coord_type=coord_type)
        return sub_geom

    def get_subgeom_without(self, indices, **kwargs):
        with_indices = [ind for ind, _ in enumerate(self.atoms) if ind not in indices]
        return self.get_subgeom(with_indices, **kwargs)

    def rmsd(self, geom, align=True):
        if not self.atoms == geom.atoms:
            raise DifferentAtomOrdering
        if align:
            return rmsd.kabsch_rmsd(
                self.coords3d - self.centroid, geom.coords3d - geom.centroid
            )
        else:
            return np.sqrt(np.mean((self.cart_coords - geom.cart_coords) ** 2))

    def as_g98_list(self):
        """Returns data for fake Gaussian98 standard orientation output.

        Returns
        -------
        g98_list : list
            List with one row per atom. Every row contains [center number,
            atomic number, atomic type (always 0 for now), X Y Z coordinates
            in Angstrom.
        """
        Atom = namedtuple("Atom", "center_num atom_num atom_type x y z")
        atoms = list()
        for i, (a, c) in enumerate(zip(self.atoms, self.coords3d), 1):
            x, y, z = c * BOHR2ANG
            atom = Atom(i, ATOMIC_NUMBERS[a.lower()], 0, x, y, z)
            atoms.append(atom)
        return atoms

    def tmp_xyz_handle(self, atoms=None, cart_coords=None):
        tmp_xyz = tempfile.NamedTemporaryFile(suffix=".xyz")
        tmp_xyz.write(self.as_xyz(atoms=atoms, cart_coords=cart_coords).encode("utf-8"))
        tmp_xyz.flush()
        return tmp_xyz

    def jmol(self, atoms=None, cart_coords=None, stdin=None, jmol_cmd="jmol"):
        """Show geometry in jmol.

        TODO: read jmol command from .pysisyphusrc ?!"""
        tmp_xyz = self.tmp_xyz_handle(atoms, cart_coords)
        cmd = [jmol_cmd, tmp_xyz.name]
        if stdin is not None:
            cmd = cmd + ["-s", "-"]
        try:
            subprocess.run(
                cmd,
                input=stdin,
                text=True,
            )
        except FileNotFoundError:
            print(f"'{jmol_cmd}' does not seem to be on your $PATH!")
        tmp_xyz.close()

    def modes3d(self):
        try:
            bonds = self.internal.bond_atom_indices
            bonds_str = " --bonds " + " ".join(map(str, it.chain(*bonds)))
        except AttributeError:
            bonds_str = ""

        tmp_xyz = self.tmp_xyz_handle()
        subprocess.run(f"modes3d.py {tmp_xyz.name}{bonds_str}", shell=True)
        tmp_xyz.close()

    def as_ase_atoms(self, vacuum=None):
        try:
            import ase
        except ImportError:
            print("Please install the 'ase' package!")
            return None

        # ASE coordinates are in Angstrom
        atoms = ase.Atoms(symbols=self.atoms, positions=self.coords3d * BOHR2ANG)
        if vacuum is not None:
            atoms.center(vacuum=vacuum)

        if self.calculator is not None:
            from pysisyphus.calculators import FakeASE

            ase_calc = FakeASE(self.calculator)
            atoms.set_calculator(ase_calc)
        return atoms

    def as_ascii_art(self) -> str:
        """ASCII-art representation of the Geometry.

        Using code from gpaw. Requires an ase installation."""
        return plot_wrapper(self)

    def get_restart_info(self):
        # Geometry restart information
        restart_info = {
            "atoms": self.atoms,
            "cart_coords": self.cart_coords.tolist(),
            "coord_type": self.coord_type,
            "comment": self.comment,
        }
        try:
            typed_prims = self.internal.typed_prims
        except AttributeError:
            typed_prims = None
        restart_info["typed_prims"] = typed_prims

        # Calculator restart information
        try:
            calc_restart_info = self.calculator.get_restart_info()
        except AttributeError:
            calc_restart_info = dict()
        restart_info["calc_info"] = calc_restart_info

        return restart_info

    def set_restart_info(self, restart_info):
        assert self.atoms == restart_info["atoms"]
        self.cart_coords = np.array(restart_info["cart_coords"], dtype=float)

        try:
            self.calculator.set_restart_info(restart_info["calc_info"])
        except KeyError:
            print("No calculator restart information found!")
        except AttributeError:
            print("Could not restart calculator, as no calculator is set!")

    def get_sphere_radius(self, offset=4):
        distances = pdist(self.coords3d)

        radius = (distances.max() / 2) + offset
        return radius

    def without_hydrogens(self):
        atoms_no_h, coords3d_no_h = zip(
            *[
                (atom, coords)
                for atom, coords in zip(self.atoms, self.coords3d)
                if atom.lower() != "h"
            ]
        )
        return Geometry(atoms_no_h, np.array(coords3d_no_h).flatten())

    def describe(self):
        return f"Geometry({self.sum_formula}, {len(self.atoms)} atoms)"

    def approximate_radius(self):
        """Approximate molecule radius from the biggest atom distance along an axis."""
        coords3d = self.coords3d - self.centroid[None, :]
        mins = coords3d.min(axis=0)
        maxs = coords3d.max(axis=0)
        dists = maxs - mins
        max_dist = dists.max()
        return max_dist

    def rotate(self, copy=False, rng=None):
        if copy:
            geom = self.copy()
        else:
            geom = self

        rot = Rotation.random(random_state=rng)
        geom.coords3d = rot.apply(geom.coords3d)
        return geom

    @property
    def bond_sets(self, bond_factor=BOND_FACTOR):
        bonds = find_bonds(
            self.atoms, self.coords3d, self.covalent_radii, bond_factor=bond_factor
        )
        bond_sets = set([frozenset(b) for b in bonds])
        return bond_sets

    def __str__(self):
        name = ""
        if self.name:
            name = f"{self.name}, "
        return f"Geometry({name}{self.sum_formula})"

    def __repr__(self):
        return self.__str__()
