import dataclasses
from typing import Optional

import numpy as np

from pysisyphus.constants import AU2KJPERMOL, AU2SEC, AU2NU, AU2EV
from pysisyphus.finite_diffs import periodic_fd_2_8
from pysisyphus.helpers_pure import highlight_text
from pysisyphus.partfuncs import partfuncs as pf
from pysisyphus.TablePrinter import TablePrinter
from pysisyphus.hindered_rotor import RotorInfo


@dataclasses.dataclass
class TorsionGPRResult:
    rotor_info: RotorInfo
    grid: np.ndarray
    energies: np.ndarray
    std: np.ndarray
    grid_train: np.ndarray
    energies_train: np.ndarray
    eigvals: np.ndarray
    eigvecs: np.ndarray
    temperature: float
    # Boltzmann weights
    weights: np.ndarray
    # TODO: Store optimized geometries?!

    @property
    def dx(self):
        return abs(self.grid[1] - self.grid[0])

    def calc_hr_partfunc(self, temperature: Optional[float] = None):
        if temperature is None:
            temperature = self.temperature
        eigvals = self.eigvals - self.eigvals.min()
        return pf.sos_partfunc(eigvals, temperature)

    def calc_cancel_partfunc(self, temperature: Optional[float] = None):
        if temperature is None:
            temperature = self.temperature
        # TODO: make calculation of force constant more robust; currently
        # index 0 is hardcoded, but this is not necessarily the minimum.
        # If an index != 0 is supposed to be used we would also have to update
        # the stored moments of inertia, as they are currently calculated for
        # the initial geometry (with index 0) only.
        #
        # Idea: add index as argument; if not provided use actual minimum and
        # recalculate moments of inertia?
        energies_cut = self.energies[:-1]
        force_constant = periodic_fd_2_8(0, energies_cut, self.dx)
        ho_freq = np.sqrt(force_constant / self.rotor_info.imom_left) / (2 * np.pi)
        ho_freq_si = ho_freq / AU2SEC
        return pf.harmonic_quantum_partfunc(ho_freq_si, temperature)

    def calc_corr_factor(self, temperature: Optional[float] = None):
        if temperature is None:
            temperature = self.temperature
        hr_partfunc = self.calc_hr_partfunc(temperature)
        cancel_partfunc = self.calc_cancel_partfunc(temperature)
        return hr_partfunc / cancel_partfunc

    @property
    def B(self):
        return 1 / (2 * self.rotor_info.imom_left)

    @property
    def BeV(self):
        return self.B * AU2EV

    def __post_init__(self):
        # Calculation of partition function correction factor (HR / HO)
        #
        # Currently, we evaluate the partition function at the initial geometry.
        # TODO: make evaluation more flexible; Maybe pick the global minimum of the scan
        # and don't always stay at the initial geometry.
        # TODO: move partfunc stuff into separate function?!
        self.hr_partfunc = self.calc_hr_partfunc()
        self.cancel_partfunc = self.calc_cancel_partfunc()
        self.corr_factor = self.calc_corr_factor()

    def report(self, boltzmann_thresh: float = 0.9999):
        print(highlight_text("Hindered Rotor Scan Summary"))
        print(
            f"             Rotor moment of inertia: {self.rotor_info.imom_left: >14.6f} au"
        )
        print(f"                        B = ħ²/(2·I): {self.BeV*1e6:>14.2f} ଼μeV")
        print(f"1d hindered rotor partition function: {self.hr_partfunc: >14.6f}")
        print(f"        HO cancel partition function: {self.cancel_partfunc: >14.6f}")
        print(f"           Correction factor (HR/HO): {self.corr_factor: >14.6f}")
        print(f"Reporting states until Σp_Boltz. >= {boltzmann_thresh}")
        print(f"                         Temperature: {self.temperature: >14.4f} K")
        header = ("# State", "E/Eh", "E/kJ mol⁻¹", "ΔE/cm⁻¹", "p_Boltz.", "Σp_Boltz.")
        col_fmts = ("{:3d}", "float", "{: 12.6f}", "{: >8.2f}", "float", "float")
        table = TablePrinter(header, col_fmts, width=12, sub_underline=False)

        weights_cum = np.cumsum(self.weights)
        nstates = np.argmax(weights_cum > boltzmann_thresh) + 1

        eigvals = self.eigvals - self.eigvals.min()
        eigvals_kJ = eigvals * AU2KJPERMOL
        eigvals_nu = eigvals * AU2NU
        nrows = len(eigvals)

        table.print_header()
        table.print_rows(
            (range(nrows), eigvals, eigvals_kJ, eigvals_nu, self.weights, weights_cum),
            first_n=nstates,
        )
