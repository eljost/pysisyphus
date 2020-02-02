from collections import Counter, namedtuple
import subprocess
import tempfile
import warnings

import numpy as np
import rmsd

from pysisyphus.constants import BOHR2ANG
from pysisyphus.elem_data import MASS_DICT, ATOMIC_NUMBERS
from pysisyphus.InternalCoordinates import RedundantCoords
from pysisyphus.intcoords.DLC import DLC
from pysisyphus.intcoords.helpers import get_tangent
from pysisyphus.linalg import gram_schmidt
from pysisyphus.xyzloader import make_xyz_str


class Geometry:

    coord_types = {
        "cart": None,
        "redund": RedundantCoords,
        "dlc": DLC,
    }

    def __init__(self, atoms, coords, coord_type="cart", comment="",
                 prim_indices=None, define_prims=None, bond_factor=1.3):
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
        coord_type : {"cart", "redund"}, optional
            Type of coordinate system to use. Right now cartesian (cart)
            and redundand (redund) are supported.
        comment : str, optional
            Comment string.
        prim_indices : iterable of shape (3, (bonds, bends, dihedrals))
            Iterable containing definitions for primitive internal coordinates.
            Three items are expected in the iterable: the first one should
            contain integer pairs, defining bonds between atoms, the next one
            should contain integer triples containing bends, and finally
            integer quadrupel for dihedrals.
        """
        self.atoms = atoms
        # self._coords always holds cartesian coordinates.
        self._coords = np.array(coords, dtype=np.float).flatten()
        assert self._coords.size == (3*len(self.atoms)), \
            f"Expected 3N={3*len(self.atoms)} cartesian coordinates but got " \
            f"{self._coords.size}. Did you accidentally supply internal " \
             "coordinates?"

        if (prim_indices is not None) and coord_type == "cart":
            coord_type = "redund"
            print("coord_type is set to 'cart' but primitive indices were "
                  "provided. Using 'redund' coord_type.")
        self.coord_type = coord_type
        coord_class = self.coord_types[self.coord_type]
        if coord_class:
            self.internal = coord_class(atoms, self._coords,
                                        prim_indices=prim_indices,
                                        define_prims=define_prims,
                                        bond_factor=bond_factor,
            )
        else:
            self.internal = None
        self.comment = comment

        self._energy = None
        self._forces = None
        self._hessian = None
        self.calculator = None

        self.masses = np.array([MASS_DICT[atom.lower()] for atom in self.atoms])
        self.total_mass = sum(self.masses)
        # Some of the analytical potentials are only 2D
        repeat_masses = 2 if (self._coords.size == 2) else 3
        self.masses_rep = np.repeat(self.masses, repeat_masses)

        atom_counter = Counter(self.atoms)
        self.sum_formula = "_".join(
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
        assert same_coord_length, "Different length of coordinate vectors!"

    def __eq__(self, other):
        return (self.atoms == other.atoms) and all(self.coords == other.coords)

    def __sub__(self, other):
        self.assert_compatibility(other)
        if self.coord_type == "cart":
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
            diff = -get_tangent(self.internal.prim_coords, other.internal.prim_coords,
                                self.internal.dihed_start)
        else:
            raise Exception("Invalid coord_type!")

        # Convert to DLC
        if self.coord_type == "dlc":
            diff = self.internal.U.T.dot(diff)
        return diff

    def copy(self, coord_type=None):
        """Returns a new Geometry object with same atoms and coordinates.

        Parameters
        ----------
        coord_type : str
            Desired coord_type, defaults to current coord_type.

        Returns
        -------
        geom : Geometry
            New Geometry object with the same atoms and coordinates.
        """
        if coord_type is None:
            coord_type = self.coord_type
        prim_indices = None
        if coord_type != "cart":
            prim_indices = self.internal.prim_indices
        return Geometry(self.atoms, self._coords, coord_type=coord_type,
                        prim_indices=prim_indices)

    def atom_indices(self):
        """Dict with atom types as key and corresponding indices as values.

        Returns
        -------
        inds_dict : dict
            Unique atom types as keys, corresponding indices as values.
        """
        unique_atoms = set(self.atoms)
        inds_dict = {}
        for atom_type in set(self.atoms):
            inds_dict[atom_type] = [i for i, atom in enumerate(self.atoms)
                                    if atom == atom_type]
        return inds_dict

    @property
    def atom_types(self):
        return set(self.atoms)

    def clear(self):
        """Reset the object state."""

        self._energy = None
        self._forces = None
        self._hessian = None

    def set_calculator(self, calculator):
        """Reset the object and set a calculator."""
        self.clear()
        self.calculator = calculator

    @property
    def mm_inv(self):
        """Inverted mass matrix.

        Returns a diagonal matrix containing the inverted atomic
        masses.
        """
        return np.diag(1/self.masses_rep)

    @property
    def mm_sqrt_inv(self):
        """Inverted square root of the mass matrix."""
        return np.diag(1/(self.masses_rep**0.5))

    @property
    def coords(self):
        """1d vector of atomic coordinates.

        Returns
        -------
        coords : np.array
            1d array holding the current coordinates.
        """
        if self.internal:
            return self.internal.coords
        return self._coords

    @coords.setter
    def coords(self, coords):
        """Wrapper for saving coordinates internally.

        Parameters
        ----------
        coords : np.array
            1d array containing atomic coordiantes. It's length
            may vary depending on the chosen coordinate system.
        """
        # Do the backtransformation from internal to cartesian.
        coords = np.array(coords).flatten()
        if self.internal:
            int_step = coords - self.internal.coords
            cart_diff = self.internal.transform_int_step(int_step)
            coords = self._coords + cart_diff
            self.internal.cart_coords = coords
        self._coords = coords
        # Reset all values because no calculations with the new coords
        # have been performed yet.
        self.clear()

    def set_coord(self, ind, coord):
        """Set a coordinate by index.

        Parameters
        ----------
        ind : int
            Index in of the coordinate to set in the self.coords array.
        coord : float
            Coordinate value.
        """
        self.coords[ind] = coord
        self.clear()

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
        self._coords = coords3d.flatten()

    @property
    def cart_coords(self):
        return self._coords

    @cart_coords.setter
    def cart_coords(self, coords):
        self._coords = coords

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
        energy_str = f"{self._energy} , " if self._energy else ""
        return f"{energy_str}{self._comment}"

    @comment.setter
    def comment(self, new_comment):
        self._comment = new_comment

    @property
    def center_of_mass(self):
        """Returns the center of mass.

        Returns
        -------
        R : np.array, shape(3, )
            Center of mass.
        """
        return 1/self.total_mass * np.sum(self.coords3d*self.masses[:,None],
                                          axis=0)

    @property
    def centroid(self):
        """Geometric center of the Geometry.

        Returns
        -------
        R : np.array, shape(3, )
            Geometric center of the Geometry.
        """
        return self.coords3d.mean(axis=0)

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

    @property
    def inertia_tensor(self):
        """Inertita tensor.

                              | x² xy xz |
        (x y z)^T . (x y z) = | xy y² yz |
                              | xz yz z² |
        """
        x, y, z = self.coords3d.T
        squares = np.sum(self.coords3d**2 * self.masses[:,None], axis=0)
        I_xx = squares[1] + squares[2]
        I_yy = squares[0] + squares[2]
        I_zz = squares[0] + squares[1]
        I_xy = -np.sum(self.masses*x*y)
        I_xz = -np.sum(self.masses*x*z)
        I_yz = -np.sum(self.masses*y*z)
        I = np.array((
                (I_xx, I_xy, I_xz),
                (I_xy, I_yy, I_yz),
                (I_xz, I_yz, I_zz)
        ))
        return I

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
        """
        # else:
            # print(vecs)
        # return aligned

        I = self.inertia_tensor
        w, v = np.linalg.eigh(I)
        rot = np.linalg.solve(v, np.eye(3))
        rot_c3d = rot.dot(c3d.T).T

        # print(c3d)
        # print()
        # print(rot_c3d)
        # print()
        # print()
        # import pdb; pdb.set_trace()
        self.coords3d = rot_c3d
        assert self.principal_axes_are_aligned()
        return self.coords3d
        """

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

    # @forces.setter
    # def forces(self, forces):
        # """Internal wrapper for setting the forces.

        # Parameters
        # ----------
        # forces : np.array
        # """
        # forces = np.array(forces)
        # assert forces.shape == self.coords.shape
        # self._forces = forces

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
        # TODO: Do the right thin here when the hessian is not yet calculated.
        #       this would probably involve figuring out how to mass-weigh and
        #       internal coordinat hessian... I think this is described in one
        #       of the Gonzales-Schlegel-papers about the GS2 algorithm.
        return self.mass_weigh_hessian(self.cart_hessian)

    def get_initial_hessian(self):
        """Return and initial guess for the hessian."""
        warnings.warn(
                "This method will be removed in the future. Get hessians from "
                "'pysisyphus.optimizers.guess_hessians' instead.",
                DeprecationWarning
        )
        if self.internal:
            H = self.internal.get_initial_hessian()
        if self.coord_type == "dlc":
            U = self.internal.U
            H = U.T.dot(H).dot(U)
        else:
            H = np.eye(self.coords.size)
        return H

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

    def get_trans_rot_vectors(self):
        """Orthonormal vectors describing translation and rotation.

        These vectors are used in the Eckart projection.

        See Martin J. Field - A Pratcial Introduction to the simulation
        of Molecular Systems, 2007, Cambridge University Press, Eq. (8.23),
        (8.24) and (8.26) for the actual projection.

        Returns
        -------
        ortho_vecs : np.array(6, atoms*3)
            2d array containing row vectors describing translations
            and rotations.
        """
        M_sqrt = np.sqrt(self.masses_rep)
        num = len(self.atoms)
        def get_trans_vecs():
            """Mass-weighted unit vectors of the three cartesian axes."""
            for vec in ((1, 0, 0), (0, 1, 0), (0, 0, 1)):
                _ = M_sqrt * np.tile(vec, num)
                yield _ / np.linalg.norm(_)
        trans_vecs = list(get_trans_vecs())

        x, y, z = self.coords3d.T
        zeros = np.zeros(x.size)

        def get_rot_vecs():
            """Mass-weighted unit vectors Rx ~ (0, -z, y, 0, -z, y, ...)
            Ry ~ (z, 0, -x, z, 0, -x, ...) and Rz ~ (-y, x, 0, -y, x, 0)."""
            for c3d in ((zeros, -z, y), (z, zeros, -x), (-y, x, zeros)):
                _ = np.array(c3d).T.flatten()
                _ *= M_sqrt
                yield _ / np.linalg.norm(_)
        rot_vecs = list(get_rot_vecs())
        ortho_vecs = np.array(gram_schmidt(trans_vecs + rot_vecs))

        return ortho_vecs

    def eckart_projection(self, mw_hessian):
        P = np.eye(self.cart_coords.size)
        for vec in self.get_trans_rot_vectors():
            P -= np.outer(vec, vec)
        return P.T.dot(mw_hessian).dot(P)

    def calc_energy_and_forces(self):
        """Force a calculation of the current energy and forces."""
        results = self.calculator.get_forces(self.atoms, self.cart_coords)
        self.set_results(results)

    def assert_cart_coords(self, coords):
        assert coords.size == self.cart_coords.size, \
            "This method only works with cartesian coordinate input. " \
            "Did you accidentally provide internal coordinates?"

    def get_energy_at(self, coords):
        self.assert_cart_coords(coords)
        return self.calculator.get_energy(self.atoms, coords)

    def get_energy_and_forces_at(self, coords):
        """Calculate forces and energies at the given coordinates.
        
        The results are not saved in the Geometry object."""
        self.assert_cart_coords(coords)
        return self.calculator.get_forces(self.atoms, coords)

    def calc_double_ao_overlap(self, geom2):
        return self.calculator.run_double_mol_calculation(self.atoms,
                                                          self.coords,
                                                          geom2.coords
        )

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
            "true_forces": "true_forces",
            "true_energy": "true_energy",
            "all_energies": "all_energies",
        }

        for key in results:
            setattr(self, trans[key], results[key])
        self.results = results

    def as_xyz(self, comment=""):
        """Current geometry as a string in XYZ-format.

        Parameters
        ----------
        comment : str, optional
            Will be written in the second line (comment line) of the
            XYZ-string.

        Returns
        -------
        xyz_str : str
            Current geometry as string in XYZ-format.
        """
        coords = self._coords * BOHR2ANG
        return make_xyz_str(self.atoms, coords.reshape((-1,3)), self.comment)

    def get_subgeom(self, indices, coord_type="cart"):
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
        ind_list = list(indices)
        sub_atoms = [self.atoms[i] for i in ind_list]
        sub_coords = self.coords3d[ind_list]
        sub_geom = Geometry(sub_atoms, sub_coords.flatten(), coord_type=coord_type)
        return sub_geom

    def rmsd(self, geom):
        return rmsd.kabsch_rmsd(self.coords3d-self.centroid,
                                geom.coords3d-geom.centroid)

    def as_g98_list(self):
        """Returns data for fake Gaussian98 standard orientation output.

        Returns
        -------
        g98_list : list
            List with one row per atom. Every row contains [center number,
            atomic number, atomic type (always 0 for now), X Y Z coordinates
            in Angstrom.
        """
        Atom = namedtuple("Atom",
                          "center_num atom_num atom_type x y z")
        atoms = list()
        for i, (a, c) in enumerate(zip(self.atoms, self.coords3d), 1):
            x, y, z = c*BOHR2ANG
            atom = Atom(i, ATOMIC_NUMBERS[a.lower()], 0, x, y, z)
            atoms.append(atom)
        return atoms

    def jmol(self):
        """Show geometry in jmol."""
        tmp_xyz = tempfile.NamedTemporaryFile(suffix=".xyz")
        tmp_xyz.write(self.as_xyz().encode("utf-8"))
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
            return None

        # ASE coordinates are in Angstrom
        atoms = ase.Atoms(symbols=self.atoms, positions=self.coords3d*BOHR2ANG)

        if self.calculator is not None:
            from pysisyphus.calculators import FakeASE

            ase_calc = FakeASE(self.calculator)
            atoms.set_calculator(ase_calc)
        return atoms

    def __str__(self):
        return f"Geometry({self.sum_formula})"

    def __repr__(self):
        return self.__str__()
