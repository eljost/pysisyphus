#!/usr/bin/env python3

import numpy as np

from pysisyphus.xyzloader import make_trj_str
from pysisyphus.constants import AMU2KG, BOHR2ANG, KB, MPERSEC2AU, AU2J


def dump_coords(atoms, coords, trj_fn):
    coords = np.array(coords)
    coords = coords.reshape(-1, len(atoms), 3) * BOHR2ANG
    trj_str = make_trj_str(atoms, coords)

    with open(trj_fn, "w") as handle:
        handle.write(trj_str)


def kinetic_energy_from_velocities(masses, velocities):
    """Kinetic energy for given velocities and masses.

    Parameters
    ----------
    masses : 2d array, (number of atoms, 1)
        Masses are expected in amu.
    velocities : 2d array, (number of atoms, 3)
        Velocities in m/s.

    Returns
    -------
    E_kin : float
        Kinetic energy in Joule.
    """
    return np.sum(masses*AMU2KG / 2 * velocities**2)


def kinetic_energy_for_temperature(atom_number, T):
    return 3/2 * atom_number * T * KB


def temperature_for_kinetic_energy(atom_number, E_kin):
    return 2/3 * E_kin / (atom_number * KB)


def temperature_for_kinetic_energy_au(atom_number, E_kin):
    return temperature_for_kinetic_energy(atom_number, E_kin*AU2J)


def get_velocities(geom, T=298.15):
    """Velocities in atomic units for the given temperature."""
    # Calculate E_kin in joule for a given temperature
    E_kin_ref = kinetic_energy_for_temperature(len(geom.atoms), T)

    velocities = np.random.rand(*geom.coords3d.shape) - 0.5
    # Calculate E_kin for the given velocities
    masses = geom.masses[:,None]
    E_kin_cur = kinetic_energy_from_velocities(masses, velocities)

    # Scale velocities to achieve the desired kinetic energy.
    # We derive the factor from E_kin_ref / E_kin_cur, but we want to
    # scale the velocities that enter quadratically into to the energy
    # expression. So we have to take the square root of the quotient.
    factor = (E_kin_ref / E_kin_cur)**0.5
    velocities_scaled = factor * velocities
    E_kin_scaled = kinetic_energy_from_velocities(masses, velocities_scaled)
    np.testing.assert_allclose(E_kin_scaled, E_kin_ref)
    return velocities_scaled * MPERSEC2AU
