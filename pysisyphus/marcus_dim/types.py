from dataclasses import dataclass
from enum import Enum
import functools

import matplotlib.pyplot as plt
import numpy as np

from pysisyphus.constants import BOHR2ANG, NU2AU

Property = Enum("Property", ["EPOS_IAO", "EPOS_MULLIKEN", "EEXC"])
RobinDay = Enum("RobinDay", ["CLASS1", "CLASS2", "CLASS3"])
PopKind = Enum("PopulationKind", ["IAO", "MULLIKEN"])


@functools.singledispatch
def get_rd_class(rd_class: RobinDay) -> RobinDay:
    return rd_class


@get_rd_class.register
def _(rd_class: str) -> RobinDay:
    return RobinDay[rd_class.upper()]


@get_rd_class.register
def _(rd_class: int) -> RobinDay:
    return {
        1: RobinDay.CLASS1,
        2: RobinDay.CLASS2,
        3: RobinDay.CLASS3,
    }[rd_class]


@dataclass
class MarcusDimScanResult:
    factors: np.ndarray
    coords: np.ndarray
    energies: np.ndarray
    properties: np.ndarray


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

    def G_adiabatic(self, x, shift=True):
        Ga, Gb = self.G_diabatic(x)
        plus = Ga + Gb
        minus2 = (Ga - Gb) ** 2
        sqrt = np.sqrt(minus2 + 4 * self.coupling**2)
        G1 = 0.5 * (plus - sqrt)
        G2 = 0.5 * (plus + sqrt)
        if shift:
            shift = G1.min()
            G1 -= shift
            G2 -= shift
        return G1, G2

    def plot(self, x, adiabatic=True, diabatic=False, show=False):
        fig, ax = plt.subplots()
        if adiabatic:
            G1, G2 = self.G_adiabatic(x)
            for label, state in (
                ("$G_1$ adia., model", G1),
                ("$G_2$ adia., model", G2),
            ):
                ax.plot(x, state, label=label)
        if diabatic:
            Ga, Gb = self.G_diabatic(x)
            for label, state in (("$G_a$", Ga), ("$G_b$", Gb)):
                ax.plot(x, state, label=label)
        ax.set_xlabel("Marcus dimension / $a_0$")
        ax.set_ylabel(r"Energy / $E_\mathrm{H}$")
        ax.legend()
        ax.set_title("Marcus Model")
        if show:
            plt.show()
        return fig, ax

    def pretty(self):
        reorg_en, dG, coupling, _, f, d = self.as_wavenums_and_ang_tuple()
        reorg_en = f"{reorg_en: >5.0f} cm⁻¹"
        dG = f"{dG: >5.0f} cm⁻¹"
        _2coupling = f"{2*coupling: >5.0f} cm⁻¹"
        d = f"{d: >6.4f} Å"
        f = f"{f: >10.1f} cm⁻¹ Å⁻²"
        return f"MarcusModel(λ={reorg_en}, ΔG={dG}, 2Vab={_2coupling}, d={d}, f={f})"
