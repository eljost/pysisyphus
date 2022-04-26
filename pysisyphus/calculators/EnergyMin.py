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
        def run_calculation(calc, name):
            results = getattr(calc, name)(atoms, coords, **prepare_kwargs)
            return results

        # Avoid unnecessary costly Hessian calculations; decide which Hessian
        # to calculate from previous energy calculation.
        if name == "get_hessian":
            tmp_energy1 = run_calculation(self.calc1, "get_energy")["energy"]
            tmp_energy2 = run_calculation(self.calc2, "get_energy")["energy"]
            if tmp_energy1 <= tmp_energy2:
                results1 = run_calculation(self.calc1, name)
                results2 = {"energy": tmp_energy2}
                self.log("Skipped Hessian calculation for calc2.")
            else:
                results1 = {"energy": tmp_energy1}
                results2 = run_calculation(self.calc2, name)
                self.log("Skipped Hessian calculation for calc1.")
        # Do both full calculation otherwise
        else:
            results1 = run_calculation(self.calc1, name)
            results2 = run_calculation(self.calc2, name)
        energy1 = results1["energy"]
        energy2 = results2["energy"]
        all_energies = np.array((energy1, energy2))

        # Mixed forces to optimize crossing points
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
            energy_diff_sq = energy_diff**2
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
                ) / denom**2 * (forces_j - forces_i)
                results["forces"] = forces
                self.log(f"norm(mixed forces)={np.linalg.norm(forces):.6f} au a0⁻¹")
            return results
        # Mixed forces end

        min_ind = [1, 0][int(energy1 < energy2)]
        en1_or_en2 = ("calc1", "calc2")[min_ind]
        results = (results1, results2)[min_ind]
        results["all_energies"] = all_energies
        en_diff_kJ = abs(energy1 - energy2) * AU2KJPERMOL

        self.log(
            f"energy_calc1={energy1:.6f} au, energy_calc2={energy2:.6f} au, returning "
            f"results for {en1_or_en2}, {en_diff_kJ: >12.2f} kJ mol⁻¹ lower."
        )
        self.calc_counter += 1
        return results

    def get_energy(self, atoms, coords, **prepare_kwargs):
        return self.do_calculations("get_energy", atoms, coords, **prepare_kwargs)

    def get_forces(self, atoms, coords, **prepare_kwargs):
        return self.do_calculations("get_forces", atoms, coords, **prepare_kwargs)

    def get_hessian(self, atoms, coords, **prepare_kwargs):
        return self.do_calculations("get_hessian", atoms, coords, **prepare_kwargs)

    def get_chkfiles(self):
        chkfiles = {}
        for key in ("calc1", "calc2"):
            try:
                chkfiles[key] = getattr(self, key).get_chkfiles()
            except AttributeError:
                pass
        return chkfiles

    def set_chkfiles(self, chkfiles):
        for key in ("calc1", "calc2"):
            try:
                getattr(self, key).set_chkfiles(chkfiles[key])
                msg = f"Set chkfile on '{key}'"
            except KeyError:
                msg = f"Found no information for '{key}' in chkfiles!"
            except AttributeError:
                msg = f"Setting chkfiles is not supported by '{key}'"
            self.log(msg)
