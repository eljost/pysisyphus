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
# [5] https://doi.org/10.1039/C8CP03273D
#     Finite-temperature Wigner phase-space sampling and temperature effects on the
#     excited-state dynamics of 2-nitronaphthalene.
#     Zobel, Nogueira, Gonzalez

import argparse
import functools
from math import exp, pi
from typing import Callable, Optional
import secrets
import sys

import matplotlib.pyplot as plt
import numpy as np
from numpy.polynomial.laguerre import Laguerre

from pysisyphus.constants import AMU2AU, AU2EV, AU2SEC, C, KB, PLANCK
from pysisyphus.io import geom_from_hessian
from pysisyphus.Geometry import Geometry


# From cm⁻¹ to angular frequency in atomic units
# cm⁻¹ * 100 -> m⁻¹
# m⁻¹ * C -> s⁻¹
NU2ANGFREQAU = 2 * pi * 100 * C * AU2SEC


def get_vib_state(
    wavenumber: float,
    rng: Optional[np.random.Generator] = None,
    temperature: Optional[float] = None,
) -> int:
    """Return random vibrational state n for given wavenumber and temperature."""
    if rng is None:
        rng = np.random.default_rng()
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
        random_state = rng.random() * thresh
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
    atoms,
    coords3d: np.ndarray,
    masses: np.ndarray,
    hessian: np.ndarray,
    temperature: Optional[float] = None,
    nu_thresh: float = 20.0,
    stddevs: float = 6.0,
    seed: Optional[int] = None,
) -> tuple[Callable, int]:
    assert coords3d.shape == (len(masses), 3)
    assert hessian.shape == (coords3d.size, coords3d.size)

    if seed is None:
        seed = secrets.randbits(128)
    rng = np.random.default_rng(seed)

    tmp_geom = Geometry(atoms, coords3d)
    tmp_geom.masses = masses
    tmp_geom.cart_hessian = hessian
    # nus, eigvals, mw_cart_displs (v), cart_displs
    nus, _, v, _ = tmp_geom.get_normal_modes()

    # Square root of angular frequencies in atomic units. Required to convert the
    # dimensionless Q and P values into atomic units.
    ang_freqs_au_sqrt = np.sqrt(nus * NU2ANGFREQAU)

    assert (nus >= nu_thresh).all(), "Imaginary wavenumbers are not yet handled!"
    nnus = len(nus)  # Number of non-zero wavenumbers

    span = 2 * stddevs
    # Pre-calculate some of the Laguerre polynomials
    # We use some shortcuts for n == 0 and n == 1
    laguerres = {
        0: lambda _: 1.0,
        1: lambda x: 1.0 - x,
    }

    def get_laguerre(n):
        lag_coeffs = np.zeros(n + 1)
        lag_coeffs[n] = 1.0
        return Laguerre(lag_coeffs)

    # Precompute some often accessed Laguere polynomials
    for i in range(2, 11):
        laguerres[i] = get_laguerre(i)

    mm_sqrt_au = np.sqrt(tmp_geom.masses_rep * AMU2AU)
    M_inv_au = np.diag(1 / mm_sqrt_au)

    def sampler():
        Qs = np.zeros(nnus)
        Ps = np.zeros(nnus)
        for i in range(nnus):
            n = get_vib_state(nus[i], rng, temperature=temperature)
            try:
                lag = laguerres[n]
            except KeyError:
                lag = get_laguerre(n)
                laguerres[n] = lag
            # According to eq. (31) in [3], the absolute value of the Wigner function is
            # bound between 0 and 1/π (see also Figure 3 in [3]). By omitting 1/π in the
            # Wigner function (eq. (28) in [3]) its absolute value will be between 0 and 1.
            # This saves us some multiplications/division as with 1/π it would also be sensible
            # to map 'ref' from [0, 1) to [0, 1/π)].
            prefact = (-1) ** n
            while True:
                q, p, ref = rng.random(3)
                # Map q and p from [0, 1) onto chosen interval [-stddevs, stddevs)
                q = q * span - stddevs
                p = p * span - stddevs
                r2 = q**2 + p**2
                # Wigner function. Eq. (2) in [2] or (28) in [3]. Both equations
                # differ in the presence or absence of a n! term. Here and in [3]
                # the n! is absorbed into the Laguerre polynomial.
                p_wig = prefact * exp(-r2) * lag(2 * r2)
                # We don't use p and q that correspond to a negative value of the
                # Wigner function. See the SI of [5] for a discussion.
                if 0.0 < p_wig >= ref:
                    break
            Qs[i] = q
            Ps[i] = p

        # The actual displacements/momenta depend on the vibrational frequencies.
        # Now we the dimensionless units to atomic units. See eq. (5) in [1].
        Qs /= ang_freqs_au_sqrt
        Ps *= ang_freqs_au_sqrt

        # Convert to Cartesian coordinates
        displ = M_inv_au @ v @ Qs
        velocities = M_inv_au @ v @ Ps

        # The COM remains unaffected, as we displace along vibrations,
        # not translations.
        displ_coords3d = coords3d.copy() + displ.reshape(-1, 3)

        # Remove rotation & translation from velocities at new coordinates.
        displ_geom = Geometry(atoms, displ_coords3d)
        displ_geom.masses = masses
        P = displ_geom.get_hessian_projector(full=True)
        velocities = P.dot(velocities)
        velocities = velocities.reshape(-1, 3)
        return displ_coords3d, velocities

    return sampler, seed


@get_wigner_sampler.register
def _(geom: Geometry, **kwargs):
    return get_wigner_sampler(
        geom.atoms, geom.coords3d, geom.masses, geom.hessian, **kwargs
    )


@get_wigner_sampler.register
def _(h5_fn: str, **kwargs):
    geom = geom_from_hessian(h5_fn)
    return get_wigner_sampler(geom, **kwargs)


def plot_normal_coords(normal_coords):
    ncoords, nnormal_coords = normal_coords.shape
    fig, ax = plt.subplots()
    ax.set_title(
        f"Normal coordinate distribution from Wigner sampling ({ncoords} samples)"
    )
    ax.axhline(0.0, c="k", ls="--", zorder=0)
    ax.violinplot(normal_coords)
    ax.set_xlabel("Normal mode")
    ax.set_ylabel("Normal coordinate")
    ax.set_xlim(0, nnormal_coords + 1)
    fig.tight_layout()
    return fig, ax


def to_normal_coords_getter(geom):
    coords_eq = geom.coords
    sqrt_masses = np.repeat(np.sqrt(geom.masses), 3)
    M_sqrt = np.diag(sqrt_masses)
    _, eigvecs = np.linalg.eigh(geom.eckart_projection(geom.mw_hessian, full=True))
    drop_first = 5 if geom.is_linear else 6
    eigvecs = eigvecs[:, drop_first:]

    def to_normal_coords(coords):
        displ = coords - coords_eq
        displ_mw = M_sqrt @ displ
        # Transform mass-weighted displacements to normal coordinates
        norm_coords = eigvecs.T @ displ_mw
        return norm_coords

    return to_normal_coords


def parse_args(args):
    parser = argparse.ArgumentParser()
    parser.add_argument("h5_fn", type=str, help="Filename of pysisyphus HDF5 Hessian.")
    parser.add_argument("-n", type=int, default=100, help="Number of samples.")
    parser.add_argument("-T", type=float, default=None, help="Temperature in K.")
    parser.add_argument("--seed", type=int)
    parser.add_argument("--plotekin", action="store_true", help="Plot kinetic energy.")
    parser.add_argument(
        "--plotnormalcoords",
        action="store_true",
        help="Plot distribution of normal coordinates",
    )
    return parser.parse_args(args)


def run():
    args = parse_args(sys.argv[1:])

    h5_fn = args.h5_fn
    n = args.n
    temperature = args.T
    seed = args.seed
    plotekin = args.plotekin
    plotnormalcoords = args.plotnormalcoords

    geom = geom_from_hessian(h5_fn)

    sampler, seed = get_wigner_sampler(geom, temperature=temperature, seed=seed)
    print(f"Seed: {seed}")
    xyzs = list()
    coords3d = np.empty((n, len(geom.atoms), 3))
    velocities = np.empty_like(coords3d)
    for i in range(n):
        coords3d[i], velocities[i] = sampler()
        xyzs.append(geom.as_xyz(cart_coords=coords3d[i]))

    trj_fn = f"samples_{n}.trj"
    with open(trj_fn, "w") as handle:
        handle.write("\n".join(xyzs))

    if plotekin:
        half_masses_au = geom.masses * AMU2AU / 2

        def E_kin(v_au):
            return (half_masses_au[:, None] * (v_au**2)).sum()

        E_kins = np.array([E_kin(v) for v in velocities])
        E_kins_eV = E_kins * AU2EV
        _, ax = plt.subplots()
        ax.hist(E_kins_eV, bins=100)
        plt.show()
    if plotnormalcoords:
        to_normal_coords = to_normal_coords_getter(geom)
        normal_coords = np.array([to_normal_coords(c3d.flatten()) for c3d in coords3d])
        fig, ax = plot_normal_coords(normal_coords)
        fig.tight_layout()
        fig.savefig("normal_coords.pdf")
