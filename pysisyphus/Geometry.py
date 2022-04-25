from collections import Counter, namedtuple
import copy
import itertools as it
import re
import subprocess
import tempfile
import sys

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

from pysisyphus import logger
from pysisyphus.config import p_DEFAULT, T_DEFAULT
from pysisyphus.constants import BOHR2ANG
from pysisyphus.elem_data import (
    MASS_DICT,
    ISOTOPE_DICT,
    ATOMIC_NUMBERS,
    COVALENT_RADII as CR,
)
from pysisyphus.helpers_pure import eigval_to_wavenumber, full_expand
from pysisyphus.intcoords import (
    DLC,
    HDLC,
    RedundantCoords,
    TRIC,
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
from pysisyphus.xyzloader import make_xyz_str


def inertia_tensor(coords3d, masses):
    """Inertita tensor.

                          | x² xy xz |
    (x y z)^T . (x y z) = | xy y² yz |
                          | xz yz z² |
    """
    x, y, z = coords3d.T
    squares = np.sum(coords3d ** 2 * masses[:, None], axis=0)
    I_xx = squares[1] + squares[2]
    I_yy = squares[0] + squares[2]
    I_zz = squares[0] + squares[1]
    I_xy = -np.sum(masses * x * y)
    I_xz = -np.sum(masses * x * z)
    I_yz = -np.sum(masses * y * z)
    I = np.array(((I_xx, I_xy, I_xz), (I_xy, I_yy, I_yz), (I_xz, I_yz, I_zz)))
    return I


def get_trans_rot_vectors(cart_coords, masses, rot_thresh=1e-6):
    """Vectors describing translation and rotation.

    These vectors are used for the Eckart projection by constructing
    a projector from them.

    See Martin J. Field - A Pratcial Introduction to the simulation
    of Molecular Systems, 2007, Cambridge University Press, Eq. (8.23),
    (8.24) and (8.26) for the actual projection.

    See also https://chemistry.stackexchange.com/a/74923.

    Parameters
    ----------
    cart_coords : np.array, 1d, shape (3 * atoms.size, )
        Atomic masses in amu.
    masses : iterable, 1d, shape (atoms.size, )
        Atomic masses in amu.

    Returns
    -------
    ortho_vecs : np.array(6, 3*atoms.size)
        2d array containing row vectors describing translations
        and rotations.
    """

    coords3d = np.reshape(cart_coords, (-1, 3))
    total_mass = masses.sum()
    com = 1 / total_mass * np.sum(coords3d * masses[:, None], axis=0)
    coords3d_centered = coords3d - com[None, :]

    I = inertia_tensor(coords3d, masses)
    _, Iv = np.linalg.eigh(I)
    Iv = Iv.T

    masses_rep = np.repeat(masses, 3)
    sqrt_masses = np.sqrt(masses_rep)
    num = len(masses)

    def get_trans_vecs():
        """Mass-weighted unit vectors of the three cartesian axes."""

        for vec in ((1, 0, 0), (0, 1, 0), (0, 0, 1)):
            _ = sqrt_masses * np.tile(vec, num)
            yield _ / np.linalg.norm(_)

    def get_rot_vecs():
        """As done in geomeTRIC."""

        rot_vecs = np.zeros((3, cart_coords.size))
        # p_vecs = Iv.dot(coords3d_centered.T).T
        for i in range(masses.size):
            p_vec = Iv.dot(coords3d_centered[i])
            for ix in range(3):
                rot_vecs[0, 3 * i + ix] = Iv[2, ix] * p_vec[1] - Iv[1, ix] * p_vec[2]
                rot_vecs[1, 3 * i + ix] = Iv[2, ix] * p_vec[0] - Iv[0, ix] * p_vec[2]
                rot_vecs[2, 3 * i + ix] = Iv[0, ix] * p_vec[1] - Iv[1, ix] * p_vec[0]
        rot_vecs *= sqrt_masses[None, :]
        return rot_vecs

    trans_vecs = list(get_trans_vecs())
    rot_vecs = np.array(get_rot_vecs())
    # Drop vectors with vanishing norms
    rot_vecs = rot_vecs[np.linalg.norm(rot_vecs, axis=1) > rot_thresh]
    tr_vecs = np.concatenate((trans_vecs, rot_vecs), axis=0)
    tr_vecs = np.linalg.qr(tr_vecs.T)[0].T
    return tr_vecs


def get_trans_rot_projector(cart_coords, masses, full=False):
    tr_vecs = get_trans_rot_vectors(cart_coords, masses=masses)
    U, s, _ = np.linalg.svd(tr_vecs.T)
    if full:
        P = np.eye(cart_coords.size)
        for tr_vec in tr_vecs:
            P -= np.outer(tr_vec, tr_vec)
    else:
        P = U[:, s.size :].T
    return P


class Geometry:

    coord_types = {
        "cart": None,
        "redund": RedundantCoords,
        "hredund": HybridRedundantCoords,
        "dlc": DLC,
        "hdlc": HDLC,
        "tric": TRIC,
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
        comment : str, optional
            Comment string.
        name : str, optional
            Verbose name of the geometry, e.g. methanal or water. Used for printing
        """
        self.atoms = atoms
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
        self.comment = comment
        self.name = name

        self._masses = None
        self._energy = None
        self._forces = None
        self._hessian = None
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
        return "_".join(
            [f"{atom.title()}{num}" for atom, num in Counter(self.atoms).items()]
        )

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
        return (self.atoms == other.atoms) and all(self.coords == other.coords)

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
        atoms = [atom for i, atom in enumerate(self.atoms) if not (i in inds)]
        c3d = self.coords3d
        coords3d = np.array(
            [c3d[i] for i, _ in enumerate(self.atoms) if not (i in inds)]
        )
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
        return np.diag(1 / (self.masses_rep ** 0.5))

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
        if self.internal:
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
                self.internal = coord_class(
                    # Instead of using only the remaining, valid typed_prims
                    # we could look for an entirely new set of typed_prims.
                    #
                    # But when we do this and we end up with more coordinates
                    # than before, this will lead to problems with the HDF5 dump.
                    # No problems arise when fewer coordinates are used
                    # (valid_typed_prims <= self.internal.typed_prims).
                    self.atoms,
                    coords3d,
                    # With typed prims, only the remaining, valid typed_prims
                    # will be defined for the new geometry.
                    #
                    # typed_prims=valid_typed_prims,
                    #
                    # With 'define_prims' the remaining, valid typed_prims
                    # will be used, together with newly determined internal
                    # coordinates. This supports, e.g., the switch from a simple
                    # bend to a linear bend and its complement.
                    define_prims=valid_typed_prims,
                )
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

    def reset_coords(self, new_typed_prims):
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
    def masses(self):
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

    def reparametrize(self):
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
        mm_sqrt = np.diag(self.masses_rep ** 0.5)
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

    def get_normal_modes(self, cart_hessian=None, full=False):
        """Normal mode wavenumbers, eigenvalues and Cartesian displacements Hessian."""
        if cart_hessian is None:
            cart_hessian = self.cart_hessian

        mw_hessian = self.mass_weigh_hessian(cart_hessian)
        proj_hessian, P = self.eckart_projection(mw_hessian, return_P=True, full=full)
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
            logger.debug(
                "Multiplicity for electronic entropy could not be determined! "
                f"Using 2S+1 = {mult}."
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
        return get_trans_rot_projector(self.cart_coords, masses=self.masses, full=full)

    def eckart_projection(self, mw_hessian, return_P=False, full=False):
        # Must not project analytical 2d potentials.
        if self.is_analytical_2d:
            return mw_hessian

        P = self.get_trans_rot_projector(full=full)
        proj_hessian = P.dot(mw_hessian).dot(P.T)
        if return_P:
            return proj_hessian, P
        else:
            return proj_hessian

    def calc_energy_and_forces(self):
        """Force a calculation of the current energy and forces."""
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

    def clear(self):
        """Reset the object state."""

        self._energy = None
        self._forces = None
        self._hessian = None
        self.true_energy = None
        self.true_forces = None
        self.true_hessian = None
        self.all_energies = None

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
        }

        for key in results:
            # Zero forces of frozen atoms
            if key == "forces":
                self.zero_frozen_forces(results[key])

            setattr(self, trans[key], results[key])
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
        if comment == "":
            comment = self.comment
        return make_xyz_str(atoms, cart_coords.reshape((-1, 3)), comment)

    def dump_xyz(self, fn):
        if not fn.lower().endswith(".xyz"):
            fn = fn + ".xyz"
        with open(fn, "w") as handle:
            handle.write(self.as_xyz())

    def get_subgeom(self, indices, coord_type="cart", sort=False):
        """Return a Geometry containing a subset of the current Geometry.

        Parameters
        ----------
        indices : iterable of ints
            Atomic indices that the define the subset of the current Geometry.
        coord_type : str, ("cart", "redund"), optional
            Coordinate system of the new Geometry.

        Returns
        -------
        sub_geom : Geometry
            Subset of the current Geometry.
        """
        if sort:
            indices = sorted(indices)
        ind_list = list(indices)
        sub_atoms = [self.atoms[i] for i in ind_list]
        sub_coords = self.coords3d[ind_list]
        sub_geom = Geometry(sub_atoms, sub_coords.flatten(), coord_type=coord_type)
        return sub_geom

    def get_subgeom_without(self, indices, **kwargs):
        with_indices = [ind for ind, _ in enumerate(self.atoms) if ind not in indices]
        return self.get_subgeom(with_indices, **kwargs)

    def rmsd(self, geom):
        return rmsd.kabsch_rmsd(
            self.coords3d - self.centroid, geom.coords3d - geom.centroid
        )

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

    def jmol(self, atoms=None, cart_coords=None):
        """Show geometry in jmol."""
        tmp_xyz = tempfile.NamedTemporaryFile(suffix=".xyz")
        tmp_xyz.write(self.as_xyz(atoms=atoms, cart_coords=cart_coords).encode("utf-8"))
        tmp_xyz.flush()
        jmol_cmd = "jmol"
        try:
            subprocess.run([jmol_cmd, tmp_xyz.name])
        except FileNotFoundError:
            print(f"'{jmol_cmd}' seems not to be on your path!")
        tmp_xyz.close()

    def as_ase_atoms(self):
        try:
            import ase
        except ImportError:
            print("Please install the 'ase' package!")
            return None

        # ASE coordinates are in Angstrom
        atoms = ase.Atoms(symbols=self.atoms, positions=self.coords3d * BOHR2ANG)

        if self.calculator is not None:
            from pysisyphus.calculators import FakeASE

            ase_calc = FakeASE(self.calculator)
            atoms.set_calculator(ase_calc)
        return atoms

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

    def rotate(self, copy=False):
        if copy:
            geom = self.copy()
        else:
            geom = self

        rot = Rotation.random()
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
