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
import dataclasses
import functools
from math import ceil, exp, pi
import secrets
import sys
from typing import Callable, Optional
import warnings

import matplotlib.pyplot as plt
import numpy as np
from numpy.polynomial.laguerre import Laguerre
import rmsd
import scipy as sp

from pysisyphus.constants import AMU2AU, AU2EV, BOHR2ANG, AU2SEC, C, KB, PLANCK
from pysisyphus.Geometry import Geometry
from pysisyphus.io import geom_from_hessian
from pysisyphus.helpers_pure import render_sp_stats


# From cm⁻¹ to angular frequency in atomic units
# cm⁻¹ * 100 -> m⁻¹
# m⁻¹ * C -> s⁻¹
NU2ANGFREQAU = 2 * pi * 100 * C * AU2SEC


@dataclasses.dataclass
class WignerSample:
    # Normal coordinates
    qs: np.ndarray
    # Momenta
    ps: np.ndarray
    # Dimensionless displacement
    qs_nodim: np.ndarray
    # Dimensionless momentum
    ps_nodim: np.ndarray
    # Displaced coordinates
    coords3d: np.ndarray
    # Equilibrium coordinates
    coords3d_eq: np.ndarray
    # Velocities
    velocities: np.ndarray


def get_vib_state(
    wavenumber: float,
    rng: Optional[np.random.Generator] = None,
    temperature: float = 0.0,
) -> int:
    """Return random vibrational state n for given wavenumber and temperature."""
    if rng is None:
        rng = np.random.default_rng()
    temperature_thresh = 1e-8
    if abs(temperature) <= temperature_thresh:
        return 0  # Ground state
    assert temperature > temperature_thresh, f"Got negative {temperature=:8.4f}!"

    # The code below is only executed when temperature is > 0 K.

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
    temperature: float = 0.0,
    nu_thresh: float = 0.0,
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
    nnus = len(nus)  # Number of non-zero wavenumbers

    imag_nus = nus <= nu_thresh
    nimag_nus = imag_nus.sum()
    use_nu_inds = np.arange(nnus)[~imag_nus]
    if nimag_nus:
        warnings.warn(
            f"Detected {nimag_nus} imaginary wavenumber(s)! They will be ignored in the sampling."
        )

    # Square root of angular frequencies in atomic units. Required to convert the
    # dimensionless Q and P values into atomic units.
    ang_freqs_au_sqrt = np.sqrt(nus * NU2ANGFREQAU)
    ang_freqs_au_sqrt[imag_nus] = 1.0

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

    def sampler() -> WignerSample:
        qs_nodim = np.zeros(nnus)
        ps_nodim = np.zeros(nnus)
        # Modes that are ignored will not be part of 'use_nu_inds' and their associated
        # sampled normal coordinates and momenta will stay at 0.0 throughout.
        for i in use_nu_inds:
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
            qs_nodim[i] = q
            ps_nodim[i] = p

        # The actual displacements/momenta depend on the vibrational frequencies.
        # Now we convert the dimensionless units to atomic units. See eq. (5) in [1].
        qs = qs_nodim / ang_freqs_au_sqrt
        ps = ps_nodim * ang_freqs_au_sqrt

        # Convert to Cartesian coordinates
        displ = M_inv_au @ v @ qs
        velocities = M_inv_au @ v @ ps

        # The COM remains unaffected, as we displace along vibrations,
        # not translations.
        displ_coords3d = coords3d.copy() + displ.reshape(-1, 3)

        # Remove rotation & translation from velocities at new coordinates.
        displ_geom = Geometry(atoms, displ_coords3d)
        displ_geom.masses = masses
        P = displ_geom.get_hessian_projector(full=True)
        velocities = P.dot(velocities)
        velocities = velocities.reshape(-1, 3)
        return WignerSample(
            qs=qs,
            ps=ps,
            qs_nodim=qs_nodim,
            ps_nodim=ps_nodim,
            coords3d=displ_coords3d,
            coords3d_eq=coords3d.copy(),
            velocities=velocities,
        )

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
    ax.axhline(0.0, c="k", ls="--", zorder=0)
    quants = [0.1, 0.9]
    quantiles = [quants] * nnormal_coords
    ax.violinplot(normal_coords, quantiles=quantiles)
    ax.set_xlabel("Normal mode")
    ax.set_ylabel("Normal coordinate")
    ax.set_title(
        f"Normal coordinates of {ncoords} geomtries.\n"
        f"Markers at quantiles [0, {', '.join(map(str, quants))}, 1.0]."
    )
    ax.set_xlim(0, nnormal_coords + 1)
    fig.tight_layout()
    return fig, ax


def plot_distances(all_coords3d):
    ncoords, natoms, _ = all_coords3d.shape
    ndists = sum(range(natoms))
    dists = list()
    for c3d in all_coords3d:
        for i, ci in enumerate(c3d):
            for cj in c3d[i + 1 :]:
                dists.append(np.linalg.norm(ci - cj))
    dists = np.array(dists) * BOHR2ANG

    fig, ax = plt.subplots()
    max_ = dists.max()
    bins = np.linspace(0.0, max_, num=50)
    ax.hist(dists, bins=bins)

    ax.set_xlim(0.0, ceil(max_))
    ax.set_xlabel("Distances / au")
    ax.set_ylabel("Count")
    ax.set_title(f"{natoms} atoms, {ncoords} geometries, {ndists} distances per geom")
    return fig, ax


def parse_args(args):
    parser = argparse.ArgumentParser()
    parser.add_argument("h5_fn", type=str, help="Filename of pysisyphus HDF5 Hessian.")
    parser.add_argument("-n", type=int, default=100, help="Number of samples.")
    parser.add_argument(
        "-T",
        type=float,
        default=0.0,
        help="Temperature in K. For T = 0 K only the vibrational ground-state will be sampled.",
    )
    parser.add_argument("--seed", type=int)
    parser.add_argument("--plotekin", action="store_true", help="Plot kinetic energy.")
    parser.add_argument(
        "--plotnormalcoords",
        action="store_true",
        help="Plot distribution of normal coordinates",
    )
    parser.add_argument(
        "--plotdistances",
        action="store_true",
        help="Plot histogram of pairwise atomic distances.",
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
    plotdistances = args.plotdistances

    geom = geom_from_hessian(h5_fn)

    sampler, seed = get_wigner_sampler(geom, temperature=temperature, seed=seed)
    print(f"Seed: {seed}, Temperature = {temperature:8.4f} K")
    xyzs = list()
    coords3d = np.empty((n, len(geom.atoms), 3))
    velocities = np.empty_like(coords3d)
    samples = list()
    for i in range(n):
        sample = sampler()
        samples.append(sample)
        coords3d[i] = sample.coords3d
        velocities[i] = sample.velocities
        xyzs.append(geom.as_xyz(cart_coords=coords3d[i]))

    c3d_ref = geom.coords3d
    rmsds = np.array([rmsd.kabsch_rmsd(c3d_ref, c3d) for c3d in coords3d])
    rmsd_stats = sp.stats.describe(rmsds)
    print("RMSD statistics:")
    print(render_sp_stats(rmsd_stats, unit="au", unit2="Å", conv2=BOHR2ANG))

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
        normal_coords = np.array([sample.qs for sample in samples])
        fig, ax = plot_normal_coords(normal_coords)
        fig.tight_layout()
        fig.savefig("normal_coords.pdf")
    if plotdistances:
        fig, ax = plot_distances(coords3d)
        fig.tight_layout()
        fig.savefig("distances.pdf")
