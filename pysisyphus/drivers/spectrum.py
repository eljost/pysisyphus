# [1] https://gaussian.com/uvvisplot/

import dataclasses
from pathlib import Path
from typing import Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
from numpy.typing import NDArray

from pysisyphus.constants import _1OVER_AU2NM, C, M_E, NA, AU2EV
from pysisyphus.config import T_DEFAULT
from pysisyphus.drivers import boltzmann


# Computation of prefactor from Gaussian whitepaper
Q_E_ESU = 4.803204e-10  # Crazy unit, cm**1.5 g**0.5 s⁻¹
C_CM = C * 100  # Speed of light in cm/s
M_E_G = M_E * 1000  # Electron mass in gram
NM2CM = 1e7  # Nanometer to centimeter
PREFACTOR = (  # Prefactor in eq. (5) of [1]
    np.sqrt(np.pi) * NA * Q_E_ESU**2 / (1e3 * np.log(10) * C_CM**2 * M_E_G)
) / NM2CM
# in Hartree
_04EV = 0.4 / AU2EV
_1EV = 1.0 / AU2EV
_2EV = 2.0 / AU2EV


@dataclasses.dataclass
class Spectrum:
    all_energies: np.ndarray
    fosc: np.ndarray
    nm: np.ndarray
    epsilon: np.ndarray

    def __post_init__(self):
        self.gs_en, *exc_ens = self.all_energies

        self.exc_ens = np.array(exc_ens) - self.gs_en
        self.exc_ens_nm = _1OVER_AU2NM / self.exc_ens

    def plot(self, **kwargs):
        return plot_spectrum(
            self.nm, self.epsilon, self.exc_ens_nm, self.fosc, **kwargs
        )

    def save(self, fn):
        act_fn = Path(fn).with_suffix(".npz")
        data = dataclasses.asdict(self)
        np.savez(act_fn, **data)
        return act_fn

    @staticmethod
    def load(fn):
        handle = np.load(fn)
        kwargs = {key: handle[key] for key in handle.files}
        return Spectrum(**kwargs)

    @staticmethod
    def from_orca(wf_fn, cis_fn, log_fn, **kwargs):
        from pysisyphus.calculators.ORCA import get_exc_ens_fosc

        exc_ens, fosc = get_exc_ens_fosc(wf_fn, cis_fn, log_fn)
        return spectrum_from_ens_fosc(exc_ens, fosc, **kwargs)

    def report(self):
        for i, (f, en) in enumerate(zip(self.fosc, self.exc_ens_nm)):
            print(f"{i: >4d}: {f=: >12.6f} ΔE={en: >8.2f} nm")


def get_grid(resolution, exc_ens, padding, min_, max_):
    from_ = max(min_, int(exc_ens.min() - padding))
    to_ = min(max_, int(exc_ens.max() + padding))
    return np.arange(from_, to_ + resolution, step=resolution)


def homogeneous_broadening(
    exc_ens: np.ndarray,
    osc_strengths: np.ndarray,
    resolution: float = 0.5,
    stddev: float = _04EV,
    from_to=None,
    grid_kwargs=None,
) -> tuple[np.ndarray, np.ndarray]:
    """Homogeneous broadening of stick spectra as outlined in Gaussian
    whitepaper [1].

    σ = 0.0147 au corresponds to about 0.4 eV. The function yields molar
    extinction coefficients in l mol cm⁻¹."""
    assert len(exc_ens) == len(osc_strengths)
    exc_ens_nm = _1OVER_AU2NM / exc_ens
    if from_to is None:
        exc_en_lowest = exc_ens[0]
        lower_energy_bound = max(exc_en_lowest - _1EV, 0.03)
        exc_en_highest = exc_ens[-1]
        upper_energy_bound = exc_en_highest + _2EV
        from_to = _1OVER_AU2NM / np.array((lower_energy_bound, upper_energy_bound))
    if grid_kwargs is None:
        grid_kwargs = {}
    _grid_kwargs = {
        "padding": 100,
    }
    _grid_kwargs.update(grid_kwargs)
    from_to = np.array(from_to)
    nm = get_grid(resolution, from_to, min_=25, max_=1500, **_grid_kwargs)
    stddev_nm = _1OVER_AU2NM / stddev

    quot = stddev_nm * (1 / nm[None, :] - (1 / exc_ens_nm[:, None]))
    exp_ = np.exp(-(quot**2))
    epsilon = PREFACTOR * osc_strengths[:, None] * stddev_nm * exp_
    epsilon = epsilon.sum(axis=0)
    return nm, epsilon


def homogeneous_broadening_on_grid(
    exc_ens: np.ndarray,
    osc_strengths: np.ndarray,
    grid: np.ndarray,
    stddev: float = _04EV,
) -> tuple[np.ndarray, np.ndarray]:
    """Homogeneous broadening of stick spectra as outlined in Gaussian
    whitepaper [1].

    σ = 0.0147 au corresponds to about 0.4 eV. The function yields molar
    extinction coefficients in l mol cm⁻¹."""
    exc_ens_nm = _1OVER_AU2NM / exc_ens
    stddev_nm = _1OVER_AU2NM / stddev

    quot = stddev_nm * (1 / grid[None, :] - (1 / exc_ens_nm[:, None]))
    exp_ = np.exp(-(quot**2))
    epsilon = PREFACTOR * osc_strengths[:, None] * stddev_nm * exp_
    epsilon = epsilon.sum(axis=0)
    return epsilon


def homogeneous_broadening_eV(
    exc_ens: np.ndarray,
    osc_strengths: np.ndarray,
    resolution: float = 0.05,
    stddev: float = 0.4,
    from_to=None,
) -> np.ndarray:
    """
    Homogeneous broadening of stick spectra as outlined in eq. (5) and (6)
    of https://doi.org/10.1063/1.4948471.

    Parameters
    ----------
    exc_ens
        Excitation energies in Hartree.
    osc_strengths
        Unitless oscillator strengths.
    resolution
        Resolution of broadened spectra in eV.
    stddev
        Standard deviation, sigma in eq. (6).
    from_to
        Limits for abscissa of spectrum in eV.

    Returns
    -------
    eV
        Spectral grid in eV.
    epsilon
        Molar extinction coefficient.
    """
    exc_ens_eV = exc_ens * AU2EV
    if from_to is None:
        from_to = np.array((exc_ens_eV[0], exc_ens_eV[-1]))
    eV = get_grid(resolution, from_to, padding=2, min_=0, max_=5)
    epsilon = (osc_strengths / (3.7922e33 * 4 * 2.926e-39 * np.sqrt(np.pi) * stddev))[
        None, :
    ] * np.exp(-(((eV[None, :] - exc_ens_eV) / stddev) ** 2))
    epsilon = epsilon.sum(axis=0)
    return eV, epsilon


def spectrum_from_ens_fosc(all_energies, fosc, **kwargs) -> Spectrum:
    gs_en, *exc_ens = all_energies
    exc_ens = np.array(exc_ens) - gs_en
    nm, epsilon = homogeneous_broadening(exc_ens, fosc, **kwargs)
    spectrum = Spectrum(
        all_energies=all_energies,
        fosc=fosc,
        nm=nm,
        epsilon=epsilon,
    )
    return spectrum


def spectrum_from_wf_tdms_ens(wf, Xa, Ya, Xb, Yb, all_energies, **kwargs) -> Spectrum:
    gs_en, *exc_ens = all_energies
    exc_ens = exc_ens - gs_en
    trans_dip_moms = wf.get_transition_dipole_moment(Xa + Ya, Xb + Yb)
    fosc = wf.oscillator_strength(exc_ens, trans_dip_moms)
    return spectrum_from_ens_fosc(all_energies, fosc, **kwargs)


def plot_spectrum(nm, epsilon, exc_ens_nm=None, fosc=None, show=False):
    fig, ax = plt.subplots()
    ax.plot(nm, epsilon, label="broadened")
    ax.set_ylabel(r"$\epsilon$")
    ax.set_xlabel(r"$\lambda$ / nm")
    # ax.legend(loc="upper left")
    axs = [
        ax,
    ]

    if (exc_ens_nm is not None) and (fosc is not None):
        ax2 = ax.twinx()
        ax2.stem(exc_ens_nm, fosc, "r", markerfmt=" ", basefmt=" ", label="fosc")
        ax2.set_ylabel("fosc")
        axs.append(ax2)
        # ax2.legend(loc="upper right")

    fig.tight_layout()
    if show:
        plt.show()

    return fig, axs


def plot_ensemble_spectrum(
    all_energies,
    fosc,
    labels,
    temperature=T_DEFAULT,
    annotate_states=True,
    annotate_thresh=0.75,
    show=False,
):
    nspectra, nstates = all_energies.shape
    assert fosc.shape == (nspectra, nstates - 1)
    # all_energies0 = spectra[0].all_energies
    gs_energies = all_energies[:, 0]
    weights = boltzmann.boltzmann_weights(gs_energies, T=temperature)
    max_weight = weights.max()
    if max_weight < annotate_thresh:
        annotate_thresh = 0.95 * max_weight
    fosc_weighted = fosc * weights[:, None]

    # exc_ens = np.array([spectrum.exc_ens for spectrum in spectra])
    exc_ens = all_energies[:, 1:] - gs_energies[:, None]
    exc_ens_nm = _1OVER_AU2NM / exc_ens

    exc_en_lowest = exc_ens.min()
    exc_en_highest = exc_ens.max()

    lower_energy_bound = max(exc_en_lowest - _1EV, 0.03)
    upper_energy_bound = exc_en_highest + _2EV
    from_to = _1OVER_AU2NM / np.array((upper_energy_bound, lower_energy_bound))
    grid = np.arange(*from_to.astype(int), step=0.05)

    broadened = np.zeros((nspectra, grid.size))
    for i in range(nspectra):
        broadened[i] = homogeneous_broadening_on_grid(exc_ens[i], fosc[i], grid=grid)
    broadened_weighted = broadened * weights[:, None]
    fin = broadened_weighted.sum(axis=0)

    fig, ax = plt.subplots()
    ax.plot(grid, fin)
    ax2 = ax.twinx()

    inds_sorted = weights.argsort()[::-1]
    # Go from conformers with high weights to conformers with low weights
    for i in inds_sorted:
        bw = broadened_weighted[i]
        kws = {}
        weight = weights[i]
        foscw = fosc_weighted[i]
        if (label := labels[i]) is not None and weight > 0.03:
            kws["label"] = f"{label} ({weight: >5.1%})"
        alpha = weight
        (lines,) = ax.plot(grid, bw, ":", alpha=alpha, **kws)
        color = lines.get_color()
        _, stem_lines, _ = ax2.stem(
            exc_ens_nm[i],
            foscw,
            linefmt=color,
            markerfmt=" ",
            basefmt=" ",
        )
        if annotate_states and weight > annotate_thresh:
            # TODO: support threshold for fosc?
            for j, foscw_j in enumerate(foscw):
                shift = min(0.1 * foscw_j, 0.025)
                xy = (exc_ens_nm[i, j], foscw_j + shift)
                text = str(j + 1)
                ax2.annotate(text, xy, ha="center", fontsize=6, color=color)
        plt.setp(stem_lines, "alpha", alpha)

    ax.set_xlim(grid[0], grid[-1])
    for ax_ in (ax, ax2):
        ax_.set_ylim(0)

    ax.set_ylabel(r"$\epsilon$")
    ax.set_xlabel(r"$\lambda$ / nm")
    ax2.set_ylabel("weighted oscillator strength")
    ax.legend()

    if show:
        plt.show()
    return fig, ax
