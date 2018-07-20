import logging

import numpy as np

from pysisyphus.constants import BOHR2ANG
from pysisyphus.xyzloader import make_xyz_str
from pysisyphus.elem_data import MASS_DICT
from pysisyphus.InternalCoordinates import RedundantCoords

class Geometry:

    coord_types = {
        "cart": None,
        "redund": RedundantCoords,
    }

    def __init__(self, atoms, coords, coord_type="cart"):
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
        """
        self.atoms = atoms
        # self._coords always holds cartesian coordinates.
        self._coords = np.array(coords)

        coord_class = self.coord_types[coord_type]
        if coord_class:
            self.internal = coord_class(atoms, coords)
        else:
            self.internal = None

        self._energy = None
        self._forces = None
        self._hessian = None
        self.calculator = None

        self.masses = [MASS_DICT[atom.lower()] for atom in self.atoms]
        # Some of the analytical potentials are only 2D
        repeat_masses = 2 if (self._coords.size == 2) else 3
        self.masses_rep = np.repeat(self.masses, repeat_masses)

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
        return self._coords.reshape(-1, 3)

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
        if self._energy and (comment == ""):
            comment = f"{comment} {self._energy}"
        return make_xyz_str(self.atoms, coords.reshape((-1,3)), comment)

    def __str__(self):
        return "Geometry, {} atoms".format(len(self.atoms))
