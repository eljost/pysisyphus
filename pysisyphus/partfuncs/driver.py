# [1] https://doi.org/10.1021/jp075861+
#     Kinetic Analysis of the Pyrolysis of Phenethyl Phenyl Ether:
#     Computational Prediction of α/β-Selectivities
#     Beste, Buchanan, Britt, Hathorn, Harrison, 2007

import dataclasses
import logging
from pathlib import Path

import numpy as np

from pysisyphus.constants import AMU2AU, AU2J, C, KBAU, PLANCK
from pysisyphus.Geometry import Geometry
import pysisyphus.numerov as numerov
from pysisyphus.partfuncs import partfuncs


pf_logger = logging.getLogger("partfunc_logger")


@dataclasses.dataclass
class Displacement:
    # Index of displacement in a mode
    ind: int
    # Index of the parent mode in the set of all normal modes
    nu_ind: int
    # Distance from starting geometry in Bohr
    dist: float
    # Displaced Cartesian coordinates in Bohr; 1d array with shape (3*natoms, )
    coords: np.ndarray


@dataclasses.dataclass
class Mode:
    ind: int
    nu: float
    cart_displs: np.ndarray
    coords: np.ndarray
    masses_rep: np.ndarray

    def __post_init__(self):
        l_mw = self.cart_displs * np.sqrt(self.masses_rep)
        self.l_mw = l_mw / np.linalg.norm(l_mw)

    def displace(self, npoints, step_size, forward=True) -> list[Displacement]:
        """Cartesian displacements along normal mode.

        Parameters
        ----------
        npoints
            Positive integer, giving the number of displacements in the selected
            direction along the normal mode.
        step_size
            Positive floating point number giving a step size in Bohr.
        forward
            Boolean flag that controls the direction of the displacement.

        Returns
        -------
        displs
            List of Displacement objects.
        """
        step = step_size * self.cart_displs
        factor = 1 if forward else -1
        displs: list[Displacement] = list()
        for i in range(1, npoints + 1):
            i *= factor
            displ_coords = i * step + self.coords
            dist = i * step_size
            displ = Displacement(ind=i, nu_ind=self.ind, dist=dist, coords=displ_coords)
            displs.append(displ)
        return displs

    def displace_full(self, npoints, step_size) -> list[Displacement]:
        """Full displacement in both directions along normal mode.

        See Mode.displace() for a full docstring.
        """
        backward = self.displace(npoints, step_size, forward=False)
        forward = self.displace(npoints, step_size, forward=True)
        displ0 = Displacement(
            ind=0, nu_ind=self.ind, dist=0.0, coords=self.coords.copy()
        )
        return backward[::-1] + [displ0] + forward

    @property
    def freq_si(self) -> float:
        """Harmonic frequency of oscillator in s⁻¹."""
        return C * self.nu * 100

    @property
    def red_mass(self) -> float:
        """Reduced mass of oscillator in amu."""
        return 1 / np.sum(np.square(self.l_mw) / self.masses_rep)

    @property
    def red_mass_au(self) -> float:
        """Reduced mass of oscillator in atomic units, not amu!."""
        return self.red_mass * AMU2AU


def get_modes_for_displacements(
    coords: np.ndarray,
    nus: np.ndarray,
    cart_displs: np.ndarray,
    masses_rep: np.ndarray,
    nu_thresh=100.0,
    logger=pf_logger,
) -> list[Mode]:
    """Determine normal modes for calculation of corrected partition functions.

    Parameters
    ----------
    coords
        1d array of shape (3*natoms, ) containing Cartesian coordinates in Bohr.
    nus
        1d array of shape (nmodes, ) containing normal mode wavenumbers in cm⁻¹.
    cart_displs
        2d array of shape (3*natoms, nmodes) containing Cartesian displacement
        vectors for all normal modes.
    masses_rep
        1d array of shape (3*natmos, ) containing triplets of atomic masses.
    nu_thresh
        Wavenumber threshold in cm⁻¹ for filtering the normal modes.
    logger
        logging.Logger object. Defaults to 'pf_logger' (see top of the module).

    Returns
    -------
    modes
        List of Mode objects.
    """
    coords = coords.flatten()
    modes: list[Mode] = list()
    for i, nu in enumerate(nus):
        if nu < 0.0 or nu > nu_thresh:
            continue
        logger.info(f"Create Mode {i:02d} at {nu: >8.2f} cm⁻¹")
        mode = Mode(
            ind=i,
            nu=nu,
            cart_displs=cart_displs[:, i],
            coords=coords,
            masses_rep=masses_rep,
        )
        modes.append(mode)
    return modes


def fit_quartic_poly(
    dists: np.ndarray, energies: np.ndarray
) -> np.polynomial.Polynomial:
    """Determine quartic potential V(x) = a0 + a1*x + a2*x**2 + a3*x**3 + a4*x**4.

    Parameters
    ----------
    dists
        1d arary of shape (ndists, ) containing distances at which the energies
        were evaluated.
    energies
        1d array of shape (ndists, ) containing energies.

    Return
    ------
    poly
        Quartic np.poly.Polynomial object containing the fitted polynomial. The
        unit(s) depend on the units of dists and energies.
    """
    poly = np.polynomial.Polynomial.fit(dists, energies, deg=4)
    return poly


def determine_xlims(
    poly: np.polynomial.Polynomial,
    dq_thresh: float,
    temperature: float,
    natoms: int,
    xgrow: float = 1.1,
    displ_per_atom: float = 0.15,
):
    """Determine suitable limits for the potential for a converged partition function.

    Parameters
    ----------
    poly
        1d polynomial of arbitrary order describing the potential.
    dq_thresh
        Threshold for the sum-over-states partition function.
    temperature
        Temperature in Kelvin.
    natoms
        Number of atoms.
    xgrow
        Floating point number > 1.0 used for upscaling the limits.
    displ_per_atom
        Mean displacment per atom; used to determine an initial value of xlim.

    Returns
    -------
    xlim
        Coordinate limit in Bohr to yield converged sos partition function for a potential
        evalulated in the range np.linspace(-xlim, xlim, num=enough_points).
    """
    kbT = KBAU * temperature
    en_max = -kbT * np.log(dq_thresh)
    en_eq = poly(0.0)
    xlim = displ_per_atom * natoms
    j = 0
    while True:
        xscale = xgrow**j
        en_left = poly(-xlim * xscale) - en_eq
        en_right = poly(xlim * xscale) - en_eq
        if (en_left > en_max) or (en_right > en_max):
            # Grow an additional time to have some margin
            xlim = xgrow * xscale * xlim
            break
        j += 1
    return xlim


def run_displs(
    geom: Geometry, displs: list[Displacement], chkfiles: str | Path
) -> tuple[np.ndarray, np.ndarray]:
    """Calculation driver for displacements along a normal mode.

    Parameters
    ----------
    geom
        Geometry object.
    displs
        List with ndists Displacement objects.
    chkfiles
        Optional chkfiles as calculated at the coordinates of Geometry. Will
        be set on the calculator if provided and supported by it.

    Returns
    -------
    dists
        1d array of shape (ndists, ) containing the distances from the initial
        geometry in Bohr.
    energies
        1d array of shape (ndists, ) containing the absolute electronic energies
        in Hartree.
    """

    # Try to set chkfiles to accelerate convergence
    if chkfiles is not None:
        try:
            geom.calculator.set_chkfile(chkfiles)
        except AttributeError:
            pass
    # Create arrays to keep the results
    dists = np.zeros(len(displs))
    energies = np.zeros(len(displs))
    # Loop over all displacements and calculate the energy at each
    for i, displ in enumerate(displs):
        energy = geom.get_energy_at(displ.coords)
        dists[i] = displ.dist
        energies[i] = energy
    return dists, energies


def run_full_displs(
    mode: Mode,
    geom: Geometry,
    energy0: float,
    step_size: float,
    chkfiles: str | Path,
    npoints: int,
) -> np.ndarray:
    """Driver to run all calculations along a normal mode.

    Parameters
    ----------
    mode
        Mode object to scan along.
    geom
        Geometry object w/ a calculator.
    energy0
        Energy at the initial/undisplaced geometry in Hartree.
    step_size
        Positive floating point number denoting the distance between two points along
        the scan.
    chkfiles
        Chkfiles that will be set before the calculations to hopefully
        increase convergence.
    npoints
        Positive integer; number of displacments to calculate in each direction.

    Returns
    -------
    potential
        2d array of shape (2 * npoints + 1, 2) containing distances in Bohr in the first
        column and energies in the second column.
    """
    # Create arrays to keep distances and energies along the Normal mode.
    dists = np.zeros(2 * npoints + 1)
    energies = np.zeros(2 * npoints + 1)
    # Set the energy of the equilibrium structure. We don't have to set the
    # displacement for the equilibrium geometry, as the arrays was already
    # initialized with zero.
    energies[npoints] = energy0

    # Backward (negative) direction
    back_displs = mode.displace(npoints, step_size, forward=False)
    back_dists, back_ens = run_displs(geom, back_displs, chkfiles)
    # Reverse energies and distances
    dists[:npoints] = back_dists[::-1]
    energies[:npoints] = back_ens[::-1]

    # Foward (positive) direction
    forward_displs = mode.displace(npoints, step_size, forward=True)
    forward_dists, forward_ens = run_displs(geom, forward_displs, chkfiles)
    dists[npoints + 1 :] = forward_dists
    energies[npoints + 1 :] = forward_ens

    potential = np.stack((dists, energies), axis=1)
    return potential


def load_or_calc(
    fn: str | Path, calc_func, io_funcs=(np.savetxt, np.loadtxt)
) -> tuple[np.ndarray, str]:
    """If 'fn' exists, load it, otherwise calculate it using 'calc_func'."""
    save_func, load_func = io_funcs
    fn = Path(fn)
    if fn.exists():
        data = load_func(fn)
        msg = f"Loaded data from '{fn}'"
    else:
        data = calc_func()
        save_func(fn, data)
        msg = f"Calculated data and saved it to '{fn}'"
    return data, msg


@dataclasses.dataclass
class PartitionFunctions:
    nu: float
    nu_ind: int
    # Semi-classical Wigner-Kirkwood partition function
    q_wk: float
    # Sum-over-states partition function
    q_aq: float
    # Anharmonic classical partition function
    q_ac: float
    # Harmonic classical partition function
    q_hc: float
    # Harmonic quantum partition function
    q_hq: float
    # Derivative of ln(sum-over-states partition function) w.r.t. temperature
    dlnq_aq_dT: float
    # Temperature in Kelvin
    temperature: float
    # Reduced mass in atomic units NOT in amu!
    red_mass_au: float
    # Eigenvalues
    eigvals: np.ndarray

    def __post_init__(self):
        if len(self.eigvals) > 1:
            # Excitation energy of fundamental in Joule
            fundamental_si = (self.eigvals[1] - self.eigvals[0]) * AU2J
            freq_fundamental_si = fundamental_si / PLANCK
            # Corrected, potentially anharmonic wavenumber of the fundamental in cm⁻¹
            self.nu_aq = freq_fundamental_si / (C * 1e2)

    """
    def render_report(self):
        report = (
            f"nu_aq={self.nu_aq: >8.2f} cm⁻¹, "
            f"q_aq={self.q_aq: >12.6f}, q_ac={self.q_ac: >12.6f}, "
            f"q_wk={self.q_wk: >12.6f}, "
            f"q_hc={self.q_hc: >12.6f}, q_hq={self.q_hq: >12.6f}"
        )
        return report
    """

    def _internal_energy(self, dlnqdT: float) -> float:
        U = KBAU * self.temperature**2 * dlnqdT
        return U

    @property
    def U_aq(self) -> float:
        """Internal energy from sum-over-states partition function in Hartree."""
        return self._internal_energy(self.dlnq_aq_dT)

    def _entropy(self, q: float, dlnqdT: float) -> float:
        S = KBAU * np.log(q) + KBAU * self.temperature * dlnqdT
        return S

    @property
    def S_aq(self) -> float:
        """Entropy from sum-over-states partition function in Hartree/Kelvin."""
        return self._entropy(self.q_aq, self.dlnq_aq_dT)


def calculate_partfuncs(
    nu: float,
    nu_ind: float,
    red_mass_au: float,
    poly: np.polynomial.Polynomial,
    temperature: float,
    natoms: int,
    # As outlined by Piccini, one could also check for the convergence of the fundamental
    dq_thresh: float = 1e-4,
    npoints_grid: int = 1000,
    nstates: int = 200,
) -> PartitionFunctions:
    """Driver that calculates partition functions for a given Mode and its potential.

    Parameters
    ----------
    nu
        Harmonic wavenumber in cm⁻¹.
    nu_ind
        Integer index of the normal modes in the list of all normal modes.
    red_mass_au
        Reduced mass of the oscillator in atomic units, not amu!
    poly
        1d np.polynomial.Polynomial representing the potential along the normal mode.
    temperature
        Temperature in Kelvin.
    dq_thresh
        Convergence threshold for the sum-over-states partition function.
    npoints_grid
        Positive integer denoting the number of grid points for the Numerov method.
    nstates
        Number of states to determine w/ the Numerov method. Must be less or equal
        than npoints_grid.

    Returns
    -------
    part_funcs
        PartitionFunctions object.
    """
    assert nstates > 0
    assert nstates <= npoints_grid

    # red_mass_au = mode.red_mass_au
    # Frequency in 1/s from speed of light C and wavenumber in m⁻¹
    freq_si = C * nu * 100

    # Numverov method start
    #
    # Energy at initial geometry
    en_eq = poly(0.0)

    def energy_getter(i, x):
        return poly(x) - en_eq

    # Pick grid wide enough
    xlim = determine_xlims(poly, dq_thresh, temperature, natoms)
    grid = np.linspace(-xlim, xlim, num=npoints_grid)
    eigvals, _ = numerov.run(
        grid=grid,
        energy_getter=energy_getter,
        mass=red_mass_au,
        nstates=200,
    )
    # Numverov method end

    # Semi-classical Wigner-Kirkwood partition function
    q_wk, _ = partfuncs.wigner_kirkwood_partfunc(poly, temperature, red_mass_au)
    # Sum-over-states partitin function
    q_aq = partfuncs.sos_partfunc(eigvals, temperature)
    # Anharmonic classical partition function
    q_ac = partfuncs.anharmonic_classic_partfunc(poly, temperature, red_mass_au)
    # Harmonic classical partition functions
    q_hc = partfuncs.harmonic_classic_partfunc(freq_si, temperature)
    # Harmonic quantum partition function
    q_hq = partfuncs.harmonic_quantum_partfunc(freq_si, temperature)

    # Derivative (d ln(q) / dT) for calculation of internal energy and entropy
    dlnq_aq_dT = partfuncs.dln_sos_partfunc_dT(eigvals, temperature)

    part_funcs = PartitionFunctions(
        nu=nu,
        nu_ind=nu_ind,
        q_wk=q_wk,
        q_aq=q_aq,
        q_ac=q_ac,
        q_hc=q_hc,
        q_hq=q_hq,
        dlnq_aq_dT=dlnq_aq_dT,
        temperature=temperature,
        red_mass_au=red_mass_au,
        eigvals=eigvals,
    )
    return part_funcs


def run(
    geom: Geometry,
    temperature: float,
    nu_thresh: float = 100,
    npoints: int = 4,
    step_size: float = 0.4,
    logger=pf_logger,
    out_dir: str | Path = ".",
) -> list[PartitionFunctions]:
    """Driver for partition function calculation.

    Parameters
    ----------
    geom
        Geometry object.
    temperature
        Temperature in Kelvin.
    nu_thresh
        Wavenumber threshold in cm⁻¹ that determines for which normal modes
        corrected partition functions will be calculated.
    npoints
        Number of points per side for the displacements. Overall (2 * npoints * nmodes + 1)
        calculations will be done. The + 1 comes from one calculation at the initial geomtry.
    step_size
        Step size in Bohr between two points along the 1d displacement.
    logger
        Optional logger.

    Returns
    -------
    all_part_funcs
        List of PartitionFunctions objects, containing corrected partition functions
        and additional information.
    """
    out_dir = Path(out_dir)
    # Calculation at equilibrium geometry to determine energy and chkfile
    energy0 = geom.energy
    try:
        chkfiles = geom.calculator.get_chkfiles()
    except AttributeError:
        chkfiles = None

    nus, *_, cart_displs = geom.get_normal_modes()
    modes = get_modes_for_displacements(
        geom.coords3d, nus, cart_displs, geom.masses_rep, nu_thresh=nu_thresh
    )

    polys = list()
    # Loop over all normal modes, calculate energies along them and fit the polynomials
    for mode in modes:
        out_fn = out_dir / f"mode_energies_{mode.ind:03d}.dat"

        potential, msg = load_or_calc(
            out_fn,
            lambda: run_full_displs(mode, geom, energy0, step_size, chkfiles, npoints),
        )
        shape_matches = potential.shape == (2 * npoints + 1)
        step_size_matches = abs(potential[0, 0] - potential[1, 0] - step_size) <= 1e-10
        # Recalcuate when arguments are different
        # TODO: make this more robust
        if not (shape_matches and step_size_matches):
            out_fn.unlink()
            potential, msg = load_or_calc(
                out_fn,
                lambda: run_full_displs(
                    mode, geom, energy0, step_size, chkfiles, npoints
                ),
            )
        dists, energies = potential.T
        logger.info(msg)
        poly = fit_quartic_poly(dists, energies)
        polys.append(poly)
    # End of loop over all normal modes

    natoms = len(geom.atoms)

    all_part_funcs = list()
    # Calculate partition functions
    for mode, poly in zip(modes, polys, strict=True):
        part_funcs = calculate_partfuncs(
            mode.nu, mode.ind, mode.red_mass_au, poly, temperature, natoms
        )
        all_part_funcs.append(part_funcs)
    return all_part_funcs
