import numpy as np

from pysisyphus.xyzloader import make_trj_str
from pysisyphus.Geometry import get_trans_rot_projector
from pysisyphus.constants import BOHR2ANG, KBAU, VELO2E


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
    masses : 1d array, shape (number of atoms, )
        Atomic masses in amu.
    velocities : 2d array, (number of atoms, 3)
        Atomic velocities in Bohr/fs.

    Returns
    -------
    E_kin : float
        Kinetic energy in Hartree.
    """
    return np.sum(masses[:,None] / 2 * velocities**2) * VELO2E


def kinetic_energy_for_temperature(atom_number, T, fixed_dof=0):
    """Kinetic energy for given temperature and number of atoms.

    Each atom has three degrees of freedom (1/2 * 3 == 3/2).

    Parameters
    ----------
    atom_number : int
        Number of atoms. Each atom has three degrees of freedom.
    T : float
        Temperature in Kelvin.
    fixed_dof : int, optional, default=0
        Number of fixed degrees of freedom, e.g. 3 when the center-of-mass
        velocity is removed.

    Returns
    -------
    E_kin : float
        Kinetic energy in Hartree.
    """
    return (3*atom_number - fixed_dof) / 2 * T * KBAU


def temperature_for_kinetic_energy(atom_number, E_kin, fixed_dof=0):
    """Temperature for given kinetic energy and atom number.

    Each atom has three degrees of freedom (1/2 * 3 == 3/2).

    Parameters
    ----------
    atom_number : int
        Number of atoms. Each atom has three degrees of freedom.
    E_kin : float
        Kinetic energy in Hartree.
    fixed_dof : int, optional, default=0
        Number of fixed degrees of freedom, e.g. 3 when the center-of-mass
        velocity is removed.

    Returns
    -------
    Temperature : float
        Temperature in Kelvin.
    """
    return 2 * E_kin / ((3*atom_number - fixed_dof) * KBAU)


def remove_com_velocity(v, masses, keep_norm=True):
    """Remove center-of-mass velocity.

    Returned units vary with the given input units.

    Parameters
    ----------
    v : np.array, 2d, shape (number of atoms, 3)
        Velocities.
    masses : np.array, 1d, shape (number of atoms, )
        Atomic masses.
    keep_norm : bool, default=True
        Whether to rescale v to its original norm, after removal of v_com.

    Returns
    -------
    v : np.array
        Velocities without center-of-mass velocity.
    """
    org_norm = np.linalg.norm(v)
    v_com = (v * masses[:,None] / np.sum(masses)).sum(axis=0)
    v -= v_com
    if keep_norm:
        v /= np.linalg.norm(v) / org_norm
    return v


def scale_velocities_to_temperatue(masses, v, T_desired, fixed_dof=0):
    """Scale velocities to a given temperature.

    Parameters
    ----------
    masses : np.array, 1d, shape (number of atoms, )
        Atomic masses in amu.
    v : np.array, 2d, shape (number of atoms, 3)
        (Unscaled) velocities in Bohr/fs.
    T_desired : float
        Desired temperature in Kelvin.
    fixed_dof : int, optional, default=0
        Number of fixed degrees of freedom, e.g. 3 when the center-of-mass
        velocity is removed.

    Returns
    -------
    v : np.array, 2d, shape (number of atoms, 3)
        Scaled velocities in Bohr/fs.
    """

    E_ref = kinetic_energy_for_temperature(len(masses), T_desired,
                                           fixed_dof=fixed_dof)  # in Hartree
    E_cur = kinetic_energy_from_velocities(masses, v)  # in Hartree
    scale = (E_ref / E_cur)**0.5
    v *= scale
    # E_now = kinetic_energy_from_velocities(masses, v)
    # np.testing.assert_allclose(E_now, E_ref)
    return v


def scale_velocities_to_energy(masses, v, E_desired):
    """Scale velocities to a given temperature.

    Parameters
    ----------
    masses : np.array, 1d, shape (number of atoms, )
        Atomic masses in amu.
    v : np.array, 2d, shape (number of atoms, 3)
        (Unscaled) velocities in Bohr/fs.
    E_desired : float
        Desired kinetic energy in Hartree.

    Returns
    -------
    v : np.array, 2d, shape (number of atoms, 3)
        Scaled velocities in Bohr/fs.
    """

    E_cur = kinetic_energy_from_velocities(masses, v)  # in Hartree
    scale = (E_desired / E_cur)**0.5
    v *= scale
    return v


def unscaled_velocity_distribution(masses, T, seed=None):
    """
    ρ ∝ exp(- v² * m / (2 * kT))
    ln(ρ) ∝ - 1/2 * v² * m / kT
    v² ∝ -2 * ln(ρ) * kT / m
    v ∝ ±(-2 * ln(ρ) * kT / m)**0.5
    v ∝ ±(-2 * k)**0.5 * (ln(ρ) * T / m)**0.5

    The first term of the RHS is constant. As we later scale the
    velocities we neglect it. Don't use these velocities unscaled!

    v ∝ ±(ln(ρ) * T / m)**0.5
    """

    if seed is not None:
        np.random.seed(seed)

    ps = np.random.standard_normal((len(masses), 3))
    v = ps * np.sqrt(T / masses[:,None])
    return v


def get_mb_velocities(masses, cart_coords, T, remove_com_v=True, remove_rot_v=True,
                      seed=None):
    """Initial velocities from Maxwell-Boltzmann distribution.

    Parameters
    ----------
    masses : np.array, 1d, shape (number of atoms, )
        Atomic masses in amu.
    cart_coords : iterable, 1d, shape (3 * number of atoms, )
        Atomic cartesian coordinates. Needed for removal of rotation.
    T : float
        Temperature in Kelvin.
    remove_com_v : bool, default=True, optional
        Whether to remove center-of-mass velocity.
    remove_rot_v : bool, default=True, optional
        Whether to remove rotational velocity.
    seed : int, default=None, optional
        Seed for the random-number-generator.

    Returns
    -------
    v : np.array, 2d, shape (number of atoms, 3)
        Initial velocities in Bohr/fs.
    """

    masses = np.array(masses)

    # Initial velocities
    v = unscaled_velocity_distribution(masses, T, seed=seed)

    if (len(masses) == 1) and remove_com_v:
        raise Exception("Removing COM velocity with only 1 atom is a bad idea!")

    fixed_dof = 0
    # Remove center-of-mass velocity
    if remove_com_v:
        v = remove_com_velocity(v, masses)
        fixed_dof += 3

    if remove_rot_v:
        P = get_trans_rot_projector(cart_coords, masses)
        v = P.dot(v.flatten()).reshape(-1, 3)
        # Right now this also removes the translational components
        fixed_dof = 6

    # In Bohr/fs
    v = scale_velocities_to_temperatue(masses, v, T, fixed_dof=fixed_dof)

    return v


def get_mb_velocities_for_geom(geom, T, remove_com_v=True, remove_rot_v=True,
                               seed=None):
    """Initial velocities from Maxwell-Boltzmann distribution.

    See 'get_mb_velocities' for explanation.
    """

    return get_mb_velocities(geom.masses, geom.cart_coords, T,
                             remove_com_v=remove_com_v,
                             remove_rot_v=remove_rot_v,
                             seed=seed,
    )


def energy_forces_getter_closure(geom):
    def energy_forces_getter(coords):
        results = geom.get_energy_and_forces_at(coords)
        energy = results["energy"]
        forces = results["forces"]
        return energy, forces
    return energy_forces_getter
