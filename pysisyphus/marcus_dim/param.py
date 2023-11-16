# [1] https://doi.org/10.1039/D3SC01402A
#     Identifying the Marcus dimension of electron transfer from
#     ab initio calculations
#     Šrut, Lear, Krewald, 2023


from dataclasses import dataclass
import functools
import math
import os
import warnings

import matplotlib.pyplot as plt
import numpy as np
import sympy as sym

from pysisyphus.constants import BOHR2ANG, NU2AU

# import pysisyphus.drivers.marcusdim_types as mdtypes


@dataclass
class MarcusModel:
    reorg_en: float  # Lambda, Hartree
    dG: float  # ΔG, barrier height, in Hartree
    coupling: float  # Vab, in Hartree
    R: float  # Distance of adiabatic minimum to top of barrier in Bohr
    f: float  # Force constant in Hartree / Bohr**2
    d: float  # Separation of diabatic states in Bohr

    def as_wavenums_and_ang_tuple(self):
        return (
            self.reorg_en / NU2AU,
            self.dG / NU2AU,
            self.coupling / NU2AU,
            self.R * BOHR2ANG,
            self.f / NU2AU / BOHR2ANG**2,
            self.d * BOHR2ANG,
        )

    def G_diabatic(self, x):
        Ga = self.f * x**2
        Gb = self.f * (x - self.d) ** 2
        return Ga, Gb

    def G_adiabatic(self, x):
        Ga, Gb = self.G_diabatic(x)
        plus = Ga + Gb
        minus2 = (Ga - Gb) ** 2
        sqrt = np.sqrt(minus2 + 4 * self.coupling)
        G1 = 0.5 * (plus - sqrt)
        G2 = 0.5 * (plus + sqrt)
        return G1, G2

    def plot(self, x, show=False):
        Ga, Gb = self.G_diabatic(x)
        G1, G2 = self.G_adiabatic(x)
        fig, ax = plt.subplots()
        for label, state in (("$G_a$", Ga), ("$G_b$", Gb)):
            ax.plot(x, state, label=label)
        for label, state in (("$G_1$", G1), ("$G_2$", G2)):
            ax.plot(x, state, label=label)
        ax.legend()
        if show:
            plt.show()
        return fig, ax

    def pretty(self):
        reorg_en, dG, coupling, *_, d = self.as_wavenums_and_ang_tuple()
        reorg_en = f"{reorg_en: >5.0f} cm⁻¹"
        dG = f"{dG: >5.0f} cm⁻¹"
        _2coupling = f"{2*coupling: >5.0f} cm⁻¹"
        d = f"{d: >6.4f} Å"
        return f"MarcusModel(λ={reorg_en}, ΔG={dG}, 2Vab={_2coupling}, d={d})"


def solve_marcus(R, Vab, dG=None, en_exc_min=None):
    """Fit Marcus-Hush model to given parameters R, Vab, dG and en_exc_min.
    Either dG or en_exc_min is mandatory.

    All argument should be given in atomic units (Hartree and Bohr).

    Quartic extension is described here:
        https://www.science.org/doi/epdf/10.1126/science.278.5339.846
    This is not yet implemented!"""

    assert (dG is not None) or (
        en_exc_min is not None
    ), "Either 'dG' or 'en_exc_min' must be given!"
    f, d = sym.symbols("f d", real=True)
    # Equation 1 is the same for both parametrizations; either eq. (5b) or (6b)
    eq1 = (
        d / 2
        - 1 / 2 * ((d**2 * f - sym.sqrt(d**4 * f**2 - 4 * Vab**2)) / (d * f))
        - R
    )

    nan = math.nan

    def inner(eq2):
        # Solve equation system using sympy
        results = sym.solve([eq1, eq2], d, f, dict=True)
        # Class III systems with only one minimum can't be solved, but we
        # return a model nonetheless.
        if len(results) == 0:
            fval = nan
            dval = nan
            reorg_en = nan
            dG = nan
        elif len(results) == 1:
            res = results[0]
            fval = res[f]
            dval = res[d]
            reorg_en = fval * dval**2
            dG = (dval**2 * fval - 2 * Vab) ** 2 / (4 * dval**2 * fval)
        else:
            raise Exception("Solving equations yielded more than one result!")
        return MarcusModel(
            reorg_en=float(reorg_en),
            dG=float(dG),
            coupling=float(Vab),
            R=float(R),
            f=float(fval),
            d=float(dval),
        )

    results = dict()
    # Depending on whether dG and/or en_exc_min we can try to run both parametrizations.
    # Use parametrization A
    if dG is not None:
        # Eq. (5c) in SI of [1]
        eq2_a = (d**2 * f - 2 * Vab) ** 2 / (4 * d**2 * f) - dG
        results["a"] = inner(eq2_a)
    # Use parametrization B
    if en_exc_min is not None:
        # Eq. (6c) in SI of [1]
        eq2_b = d**2 * f - en_exc_min
        results["b"] = inner(eq2_b)
    return results


def solve_marcus_wavenums_and_ang(R, Vab, dG=None, en_exc_min=None):
    """Wrapper for solve_marcus that accepts wavenumbers and Ångstrom."""
    kwargs = {
        "R": R / BOHR2ANG,
        "Vab": Vab * NU2AU,
        "dG": dG * NU2AU if dG is not None else None,
        "en_exc_min": en_exc_min * NU2AU if en_exc_min is not None else None,
    }
    return solve_marcus(**kwargs)


def find_minima(arr):
    assert (arr.ndim == 1) and (len(arr) >= 3)
    first_min = [0] if arr[0] < arr[1] else []
    last_min = [len(arr) - 1] if arr[-1] < arr[-2] else []
    inner_mins = [i for i in range(1, arr.size) if arr[i - 1] > arr[i] < arr[i + 1]]
    return first_min + inner_mins + last_min


@functools.singledispatch
def param_marcus(coordinate: np.ndarray, energies: np.ndarray):
    """Parametrize Marcus model with results from scan along Marcus dimension."""
    assert coordinate.ndim == 1, (
        "Parametrization requires a 1d coordinate, e.g. an "
        "array collecting displacement along the Marcus dimension!"
    )
    assert (
        energies.ndim == 2
    ), "Parametrization requires potential energy curves of 2 states!"

    energies = energies - energies.min()
    # Excitation energy at adiabatic minimum
    min_inds = find_minima(energies[:, 0])  # Search minima in lower state
    if len(min_inds) == 2:
        ind0, ind1 = min_inds
        adia_min_ind = ind0 if energies[ind0, 0] < energies[ind1, 0] else ind1
        energies_between_mins = energies[ind0 : ind1 + 1]
        # Determine index of barrier in lower state
        barr_ind = energies_between_mins[:, 0].argmax()
        barr_ind += ind0
    elif len(min_inds) == 1:
        warnings.warn(
            "Found class III system. Parametrizing Marcus model is not possible!"
        )
        adia_min_ind = barr_ind = min_inds[0]
    else:
        raise Exception("How did I get here?!")

    adia_min_ens = energies[adia_min_ind]
    en_exc_min = adia_min_ens[1] - adia_min_ens[0]
    barr_ens = energies[barr_ind]
    en_exc_barr = barr_ens[1] - barr_ens[0]

    # Electronic coupling
    Vab = en_exc_barr / 2
    # Distance R between adiabatic minimum and top of barrier
    R = abs(coordinate[barr_ind] - coordinate[adia_min_ind])
    # Barrier height ΔG
    dG = energies[barr_ind, 0]

    return solve_marcus(R, Vab, dG=dG, en_exc_min=en_exc_min)


@param_marcus.register
def _(path: os.PathLike):
    data = np.load(path)
    factors = data["factors"]
    energies = data["energies"]
    return param_marcus(factors, energies)
