from typing import Tuple
import warnings

import matplotlib.pyplot as plt
import numpy as np
from numpy.typing import NDArray

from pysisyphus.constants import AU2J, C, M_E, NA, PLANCK, AU2EV
from pysisyphus.wavefunction import Wavefunction
from pysisyphus.calculators.ORCA import parse_orca_cis, parse_orca_all_energies


# Computation of prefactor from Gaussian whitepaper
Q_E_ESU = 4.803204e-10  # Crazy unit, cm**1.5 g**0.5 s⁻¹
C_CM = C * 100  # Speed of light in cm/s
M_E_G = M_E * 1000  # Electron mass in gram
NM2CM = 1e7  # Nanometer to centimeter
PREFACTOR = (  # Prefactor in eq. (5) of [1]
    np.sqrt(np.pi) * NA * Q_E_ESU ** 2 / (1e3 * np.log(10) * C_CM ** 2 * M_E_G)
) / NM2CM
_04EV = 0.4 / AU2EV
# Factor used in converting energy in Hartree to wavelength
_AU2NM = PLANCK * C * 1e9 / AU2J


def au2nm(au):
    return _AU2NM / au


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
        from_ = min(300, int(exc_ens_nm.min() - 100))
        to_ = max(900, int(exc_ens_nm.max() + 100))
    else:
        from_, to_ = from_to

    nm = np.arange(from_, to_ + resolution, step=resolution)
    stddev_nm = au2nm(stddev)

    quot = stddev_nm * (1 / nm[None, :] - (1 / exc_ens_nm[:, None]))
    exp_ = np.exp(-(quot ** 2))
    epsilon = PREFACTOR * osc_strengths[:, None] * stddev_nm * exp_
    epsilon = epsilon.sum(axis=0)
    return nm, epsilon


def from_orca(json_fn, cis_fn, log_fn):
    wf = Wavefunction.from_orca_json(json_fn)
    Xa, Ya, Xb, Yb = parse_orca_cis(cis_fn)
    all_energies = parse_orca_all_energies(log_fn, do_tddft=True)

    exc_ens = all_energies[1:] - all_energies[0]
    tdens = wf.transition_dipole_moment(Xa + Ya)
    warnings.warn("Only alpha TDM is currently taken into account!")
    fosc = wf.oscillator_strength(exc_ens, tdens)
    return exc_ens, fosc


def plot_spectrum(nm, epsilon, exc_ens_nm=None, fosc=None, show=False):
    fig, ax = plt.subplots()
    ax.plot(nm, epsilon)
    ax.set_ylabel("$\epsilon$")
    ax.set_xlabel("$\lambda$ / nm")
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
