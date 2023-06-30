# [1] https://gaussian.com/uvvisplot/

import dataclasses
from pathlib import Path
from typing import Tuple

import matplotlib.pyplot as plt
import numpy as np
from numpy.typing import NDArray

from pysisyphus.constants import AU2J, C, M_E, NA, PLANCK, AU2EV


# Computation of prefactor from Gaussian whitepaper
Q_E_ESU = 4.803204e-10  # Crazy unit, cm**1.5 g**0.5 s⁻¹
C_CM = C * 100  # Speed of light in cm/s
M_E_G = M_E * 1000  # Electron mass in gram
NM2CM = 1e7  # Nanometer to centimeter
PREFACTOR = (  # Prefactor in eq. (5) of [1]
    np.sqrt(np.pi) * NA * Q_E_ESU**2 / (1e3 * np.log(10) * C_CM**2 * M_E_G)
) / NM2CM
_04EV = 0.4 / AU2EV  # in Hartree
# Factor used in converting energy in Hartree to wavelength
_AU2NM = PLANCK * C * 1e9 / AU2J


@dataclasses.dataclass
class Spectrum:
    exc_ens: NDArray[float]
    exc_ens_nm: NDArray[float]
    fosc: NDArray[float]
    nm: NDArray[float]
    epsilon: NDArray[float]

    def plot(self, **kwargs):
        plot_spectrum(self.nm, self.epsilon, self.exc_ens_nm, self.fosc, **kwargs)

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


def au2nm(au):
    return _AU2NM / au


def get_grid(resolution, exc_ens, padding, min_, max_):
    from_ = max(min_, int(exc_ens.min() - padding))
    to_ = min(max_, int(exc_ens.max() + padding))
    return np.arange(from_, to_ + resolution, step=resolution)


def homogeneous_broadening(
    exc_ens: NDArray[float],
    osc_strengths: NDArray[float],
    resolution: float = 0.5,
    stddev: float = _04EV,
    from_to=None,
) -> Tuple[NDArray[float], NDArray[float]]:
    """Homogeneous broadening of stick spectra as outlined in Gaussian
    whitepaper [1].

    σ = 0.0147 au corresponds to about 0.4 eV. The function yields molar
    extinction coefficients in l mol cm⁻¹."""
    exc_ens_nm = au2nm(exc_ens)
    if from_to is None:
        from_to = np.array((exc_ens_nm[0], exc_ens_nm[-1]))
    nm = get_grid(resolution, from_to, padding=100, min_=100, max_=900)
    stddev_nm = au2nm(stddev)

    quot = stddev_nm * (1 / nm[None, :] - (1 / exc_ens_nm[:, None]))
    exp_ = np.exp(-(quot**2))
    epsilon = PREFACTOR * osc_strengths[:, None] * stddev_nm * exp_
    epsilon = epsilon.sum(axis=0)
    return nm, epsilon


def spectrum_from_ens_fosc(exc_ens, fosc, **kwargs):
    exc_ens_nm = au2nm(exc_ens)
    nm, epsilon = homogeneous_broadening(exc_ens, fosc, **kwargs)
    spectrum = Spectrum(
        exc_ens=exc_ens,
        exc_ens_nm=exc_ens_nm,
        fosc=fosc,
        nm=nm,
        epsilon=epsilon,
    )
    return spectrum


def homogeneous_broadening_eV(
    exc_ens: NDArray[float],
    osc_strengths: NDArray[float],
    resolution: float = 0.05,
    stddev: float = 0.4,
    from_to=None,
) -> NDArray[float]:
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


def plot_spectrum(nm, epsilon, exc_ens_nm=None, fosc=None, show=False):
    fig, ax = plt.subplots()
    ax.plot(nm, epsilon)
    ax.set_ylabel(r"$\epsilon$")
    ax.set_xlabel(r"$\lambda$ / nm")
    axs = [
        ax,
    ]

    if (exc_ens_nm is not None) and (fosc is not None):
        ax2 = ax.twinx()
        ax2.stem(exc_ens_nm, fosc, "r", markerfmt=" ", basefmt=" ")
        ax2.set_ylabel("fosc")
        axs.append(ax2)
    fig.tight_layout()
    if show:
        plt.show()

    return fig, axs
