# [1] https://doi.org/10.1063/1.3463717
#     Comparisons of classical and Wigner sampling of transition state
#     energy levels for quasiclassical trajectory chemical dynamics simulations
#     Sun, Hase, 2010
# [2] https://doi.org/10.1002/qua.25049
#     Effects of different initial condition samplings on photodynamics and
#     spectrum of pyrrole
#     Barbatti, Sen, 2015
# [3] https://doi.org/10.1063/1.453761
#     The Morse oscillator in position space, momentum space, and phase space
#     Dahl, Springborg, 1988
#
# An interesting paper also seems to be (but i never read it ...)
#
# [4] https://doi.org/10.1063/5.0039592
#     Sampling initial positions and momenta for nuclear trajectories
#     from quantum mechanical distributions
#     Yao, Hase, Granucci, Persico

import functools
from math import exp, pi
from typing import Callable, Optional

import numpy as np
from numpy.typing import NDArray
from numpy.polynomial.laguerre import Laguerre

from pysisyphus.constants import AMU2AU, AMU2KG, BOHR2M, C, HBAR, KB, P_AU, PLANCK
from pysisyphus.helpers_pure import eigval_to_wavenumber
from pysisyphus.Geometry import Geometry, get_trans_rot_projector


def get_vib_state(wavenumber: float, temperature: Optional[float] = None) -> int:
    """Return random vibrational state n for given wavenumber and temperature."""
    if temperature is None:
        return 0  # Ground state

    freq = wavenumber * 100 * C  # from cm⁻¹ to s⁻¹
    vib_en_J = PLANCK * freq  # Energy in J
    quot = vib_en_J / (KB * temperature)
    Z = np.exp(-quot / 2) / (1 - np.exp(-quot))  # Partition function
    _1overZ = 1 / Z
    # TODO: check for sensible values?

    def get_p(n):
        """Probability of vibrational state with quantum number n.

        Given by a Boltzmann distribution:

                   exp(-e_n / kT)
        p_n = ----------------------------
              sum_{n = 0}^N exp(-e_n / kT)

        for a harmonic oscillator e_n = h * freq * (n + 1/2)

        See also https://chemistry.stackexchange.com/a/61120.

        """
        return np.exp(-quot * (n + 0.5)) * _1overZ

    # Determine sensible maximum vibrational state n.
    n = 0
    probabilities = list()
    probability_sum = 0.0
    # 1.0 may never be reached, so we stop earlier.
    thresh = 0.999999
    while True:
        p = get_p(n)
        probabilities.append(p)
        probability_sum += p
        if probability_sum >= thresh:
            break
        n += 1

    # Generate random number that is smaller than the current sum.
    while True:
        # Sample from the possible interval
        random_state = np.random.random() * thresh
        if random_state < probability_sum:
            break

    cumsum = np.cumsum(probabilities)
    for n, cs in enumerate(cumsum):
        if cs >= random_state:
            break
    return n


def normal_mode_reduced_masses(masses_rep, normal_modes):
    return 1 / (np.square(normal_modes) / masses_rep[:, None]).sum(axis=0)


@functools.singledispatch
def get_wigner_sampler(
    coords3d: NDArray,
    masses: NDArray,
    hessian: NDArray,
    temperature: Optional[float] = None,
    nu_thresh: float = 5.0,
    stddevs: float = 6.0,
) -> Callable:
    assert coords3d.shape == (len(masses), 3)
    assert hessian.shape == (coords3d.size, coords3d.size)

    # Projector to remove translation & rotation
    Proj = get_trans_rot_projector(coords3d, masses, full=True)
    masses_rep = np.repeat(masses, 3)
    mm_sqrt = np.sqrt(masses_rep)
    M = np.diag(mm_sqrt)
    M_inv = np.diag(1 / mm_sqrt)
    PM = Proj @ M_inv

    # Diagonalize projected, mass-weighted Hessian.
    w, v = np.linalg.eigh(PM @ hessian @ PM.T)
    nus = eigval_to_wavenumber(w)
    small_nu_mask = np.abs(nus) < nu_thresh
    w = w[~small_nu_mask]
    v = v[:, ~small_nu_mask]
    nus = nus[~small_nu_mask]
    assert (nus >= nu_thresh).all(), "Imaginary wavenumbers are not yet handled!"
    nnus = len(nus)  # Number of non-zero wavenumbers

    # Convert wavenumbers (cm⁻¹) to angular frequency in 1/s
    ang_freqs = 2 * pi * nus * 100 * C

    # Reduced masses in amu and then in kg
    mus_amu = normal_mode_reduced_masses(masses_rep, v)
    mus = mus_amu * AMU2KG

    # Standard deviations. Eq. (3) in [2].
    sigma2_q = HBAR / (2 * mus * ang_freqs)
    sigma2_p = (HBAR * mus * ang_freqs) / 2
    # Sigmas are now in atomic units!
    sigma2_q = sigma2_q / (BOHR2M**2)
    sigma2_p = sigma2_p / (P_AU**2)
    q_quot = 1 / (2 * sigma2_q)
    p_quot = 1 / (2 * sigma2_p)

    span = 2 * stddevs

    def sampler():
        Qs = np.zeros(nnus)
        Ps = np.zeros(nnus)
        for i in range(nnus):
            qi_quot = q_quot[i]
            pi_quot = p_quot[i]
            n = get_vib_state(nus[i], temperature)
            lag_coeffs = np.zeros(n + 1)
            lag_coeffs[n] = 1.0
            lag = Laguerre(lag_coeffs)
            prefact = (-1)**n / pi
            while True:
                q, p, ref = np.random.random_sample(3)
                # Map q and p from [0, 1) onto chosen interval [-stddevs, stddevs)
                q = q * span - stddevs
                p = p * span - stddevs
                r2 = q**2 * qi_quot + p**2 * pi_quot
                # Wigner function. Eq. (2) in [2] or (28) in [3]. Both equations
                # differ in the presence or absence of a n! term. Here and in [3]
                # the n! is absorbed into the Laguerre polynomial.
                p_wig = prefact * exp(-r2) * lag(2 * r2)
                if p_wig >= ref:
                    break
            Qs[i] = q
            Ps[i] = p

        # Convert to Cartesian coordinates
        displ = M_inv @ v @ Qs
        velocities = M @ v @ Ps

        # The COM remains unaffected, as we displace along vibrations,
        # not translations.
        displ_coords3d = coords3d.copy() + displ.reshape(-1, 3)

        # dp is still in atomic units with electron_mass as mass unit
        velocities /= masses_rep * AMU2AU
        P = get_trans_rot_projector(displ_coords3d, masses)
        velocities = P.dot(velocities)
        return displ_coords3d, velocities

    return sampler


@get_wigner_sampler.register
def _(geom: Geometry, **kwargs):
    return get_wigner_sampler(geom.coords3d, geom.masses, geom.hessian, **kwargs)
