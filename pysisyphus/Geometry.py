from collections import Counter
import logging

import numpy as np
import rmsd

from pysisyphus.constants import BOHR2ANG
from pysisyphus.xyzloader import make_xyz_str
from pysisyphus.elem_data import MASS_DICT
from pysisyphus.InternalCoordinates import RedundantCoords

class Geometry:

    coord_types = {
        "cart": None,
        "redund": RedundantCoords,
    }

    def __init__(self, atoms, coords, coord_type="cart", comment=""):
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
        """
        self.atoms = atoms
        # self._coords always holds cartesian coordinates.
        self._coords = np.array(coords)

        self.coord_type = coord_type
        coord_class = self.coord_types[self.coord_type]
        if coord_class:
            self.internal = coord_class(atoms, coords)
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

    def __eq__(self, other):
        return (self.atoms == other.atoms) and all(self.coords == other.coords)

    def copy(self):
        """Returns a new Geometry object with same atoms and coordinates.

        Returns
        -------
        geom : Geometry
            New Geometry object with the same atoms and coordinates.
        """
        return Geometry(self.atoms, self._coords, self.coord_type)

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
        return np.linalg.inv(np.diag(self.masses_rep))

    @property
    def mm_sqrt_inv(self):
        """Inverted square root of the mass matrix."""
        return np.linalg.inv(np.sqrt(np.diag(self.masses_rep)))

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
        if self.internal:
            int_step = coords - self.internal.coords
            cart_diff = self.internal.transform_int_step(int_step)
            coords = self._coords + cart_diff
            self.internal.cart_coords = coords
        self._coords = coords
        # Reset all values because no calculations with the new coords
        # have been performed yet.
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
        self.coords = coords3d.flatten()

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
    def forces(self):
        """Energy of the current atomic configuration.

        Returns
        -------
        force : np.array
            1d array containing the forces acting on the atoms. Negative
            of the gradient.
        """
        if self._forces is None:
            results = self.calculator.get_forces(self.atoms, self._coords)
            self.set_results(results)
        if self.internal:
            return self.internal.transform_forces(self._forces)
        return self._forces

    @forces.setter
    def forces(self, forces):
        """Internal wrapper for setting the forces.

        Parameters
        ----------
        forces : np.array
        """
        #if self.internal:
        #    raise Exception("Setting forces in internal coordinates not "
        #                    "yet implemented!")
        self._forces = forces

    @property
    def gradient(self):
        """Negative of the force.

        Returns
        -------
        gradient : np.array
            1d array containing the negative of the current forces.
        """
        return -self.forces

    @gradient.setter
    def gradient(self, gradient):
        """Internal wrapper for setting the gradient."""
        self._forces = -gradient

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
    def hessian(self):
        """Matrix of second derivatives of the energy in respect to atomic
        displacements.


        Returns
        -------
        hessian : np.array
            2d array containing the second derivatives of the energy in respect
            to atomic displacements.
        """
        if self._hessian is None:
            results = self.calculator.get_hessian(self.atoms, self._coords)
            self.set_results(results)
        if self.internal:
            raise Exception("Hessian in internal coordinates not implemented!")
        return self._hessian

    @property
    def mw_hessian(self):
        """Mass-weighted hessian.

        Returns
        -------
        mw_hessian : np.array
            2d array containing the mass-weighted hessian M^(-1/2) H M^(-1/2).
        """
        # M^(-1/2) H M^(-1/2)
        return self.mm_sqrt_inv.dot(self.hessian).dot(self.mm_sqrt_inv)

    @hessian.setter
    def hessian(self, hessian):
        """Internal wrapper for setting the hessian."""
        self._hessian = hessian

    def get_initial_hessian(self):
        """Return and initial guess for the hessian."""
        if self.internal:
            return self.internal.get_initial_hessian()
        return np.eye(self.coords.size)

    def calc_energy_and_forces(self):
        """Force a calculation of the current energy and forces."""
        results = self.calculator.get_forces(self.atoms, self.coords)
        self.set_results(results)

    def get_energy_and_forces_at(self, coords):
        """Calculate forces and energies at the given coordinates.
        
        The results are not saved in the Geometry object."""
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

        for key in results:
            setattr(self, key, results[key])
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

    def __str__(self):
        return f"Geometry({self.sum_formula})"

    def __repr__(self):
        return self.__str__()
