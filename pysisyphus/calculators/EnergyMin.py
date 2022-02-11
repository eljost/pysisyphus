# [1] https://doi.org/10.1063/5.0021923
#     Extending NEB method to reaction pathways involving multiple spin states
#     Zhao et al, 2020
import numpy as np

from pysisyphus.calculators.Calculator import Calculator
from pysisyphus.constants import AU2KCALPERMOL, AU2KJPERMOL


class EnergyMin(Calculator):
    def __init__(self, calculator1, calculator2, **kwargs):
        super().__init__(**kwargs)
        self.calc1 = calculator1
        self.calc2 = calculator2

        self.mix = False

    def do_calculations(self, name, atoms, coords, **prepare_kwargs):
        results1 = getattr(self.calc1, name)(atoms, coords, **prepare_kwargs)
        results2 = getattr(self.calc2, name)(atoms, coords, **prepare_kwargs)
        energy1 = results1["energy"]
        energy2 = results2["energy"]
        all_energies = np.array((energy1, energy2))

        if self.mix:
            alpha = 0.02 / AU2KCALPERMOL
            sigma = 3.5 / AU2KCALPERMOL
            # Must be positive, so substract lower energy from higher energy.
            # I is the higher state, j the lower one.
            en_i, en_j = (
                (energy2, energy1) if (energy1 < energy2) else (energy1, energy2)
            )
            energy_diff = en_i - en_j
            assert energy_diff > 0.0
            self.log(
                f"Mix mode, ΔE={energy_diff:.6f} au ({energy_diff*AU2KJPERMOL:.2f} kJ mol⁻¹)"
            )
            energy_diff_sq = energy_diff ** 2
            denom = energy_diff + alpha
            energy = (en_i + en_j) / 2 + sigma * energy_diff_sq / denom
            results = {
                "energy": energy,
                "all_energies": all_energies,
            }
            self.log(f"Mixed energy: {energy:.6f} au")
            if name == "get_forces":
                forces1 = results1["forces"]
                forces2 = results2["forces"]
                forces_i, forces_j = (
                    (forces2, forces1) if (energy1 < energy2) else (forces1, forces2)
                )
                forces = (forces_i + forces_j) / 2 + sigma * (
                    energy_diff_sq + 2 * alpha * energy_diff
                ) / denom ** 2 * (forces_j - forces_i)
                results["forces"] = forces
                self.log(f"norm(mixed forces)={np.linalg.norm(forces):.6f} au a0⁻¹")
            return results

        min_ind = [1, 0][int(energy1 < energy2)]
        en1_or_en2 = ("calc1", "calc2")[min_ind]
        results = (results1, results2)[min_ind]
        results["all_energies"] = all_energies
        en_diff_kJ = abs(energy1 - energy2) * AU2KJPERMOL

        self.log(
            f"energy_calc1={energy1:.6f} au, energy_calc2={energy2:.6f} au, returning "
            f"results for {en1_or_en2}, {en_diff_kJ:.2f} kJ mol⁻¹ lower."
        )
        return results

    def get_energy(self, atoms, coords, **prepare_kwargs):
        return self.do_calculations("get_energy", atoms, coords, **prepare_kwargs)

    def get_forces(self, atoms, coords, **prepare_kwargs):
        return self.do_calculations("get_forces", atoms, coords, **prepare_kwargs)

    def get_hessian(self, atoms, coords, **prepare_kwargs):
        return self.do_calculations("get_hessian", atoms, coords, **prepare_kwargs)
