# [1] https://doi.org/10.1063/5.0021923
#     Extending NEB method to reaction pathways involving multiple spin states
#     Zhao et al, 2020
# [2] https://doi.org/10.1021/jp0761618
#     Optimizing Conical Intersections without Derivative Coupling Vectors:
#     Application to Multistate Multireference Second-Order Perturbation Theory
#     Levine, Coe, Martinez 2008

import numpy as np

from pysisyphus.calculators.Calculator import Calculator
from pysisyphus.constants import AU2KJPERMOL
from pysisyphus.helpers import norm_max_rms


class EnergyMin(Calculator):
    def __init__(
        self,
        calculator1: Calculator,
        calculator2: Calculator,
        mix: bool = False,
        alpha: float = 0.02,  # Hartree
        sigma: float = 3.5,  # Unitless; default for ethene case in [1]
        min_energy_diff: float = 0.0,
        check_after: int = 0,
        **kwargs,
    ):
        """
        Use energy and derivatives of the calculator with lower energy.

        This calculators carries out two calculations with different settings
        and returns the results of the lower energy one. This can be used
        to consider flips between a singlet and a triplet PES etc.

        Parameters
        ----------
        calculator1
            Wrapped QC calculator that provides energies and its derivatives.
        calculator2
            Wrapped QC calculator that provides energies and its derivatives.
        mix
            Enable mixing of both forces, according to the approach outlined
            in [2]. Can be used to optimize guesses for MECPs.
            Pass
        alpha
            Smoothing parameter in Hartree. See [2] for a discussion.
        sigma
            Unitless gap size parameter. The final gap becomes
            smaller for bigga sigmas. Has to be adapted for each case. See
            [2] for a discussion (p. 407 right column and p. 408 left column.)
        min_energy_diff
            Energy difference in Hartree. When set to a value != 0 and the
            energy difference between both
            calculators drops below this value, execution of both calculations
            is diabled for 'check_after' cycles. In these cycles the calculator choice
            remains fixed. After 'check_after' cycles, both energies
            will be calculated and it is checked, if the previous calculator
            choice remains valid.
            In conjunction with 'check_after' both arguments can be used to
            save computational ressources.
        check_after
            Amount of cycles in which the calculator choice remains fixed.

        Other Parameters
        ----------------
        **kwargs
            Keyword arguments passed to the Calculator baseclass.
        """
        super().__init__(**kwargs)
        self.calc1 = calculator1
        self.calc2 = calculator2
        self.alpha = alpha
        self.sigma = sigma
        self.min_energy_diff = float(min_energy_diff)
        self.check_after = int(check_after)
        assert self.check_after >= 0, "'check_after' must not be negative!"

        self.mix = mix
        self.recalc_in = self.check_after
        self.fixed_calc = None

    def do_calculations(self, name, atoms, coords, **prepare_kwargs):
        def run_calculation(calc, name=name):
            results = getattr(calc, name)(atoms, coords, **prepare_kwargs)
            return results

        if (self.fixed_calc is not None) and self.recalc_in > 0:
            self.log(
                f"Used fixed calculator '{self.fixed_calc}'. Re-checking both "
                f"calculators in {self.recalc_in} cycles."
            )
            results = run_calculation(self.fixed_calc)
            self.recalc_in -= 1
            self.calc_counter += 1
            return results
        elif (self.fixed_calc is not None) and self.recalc_in == 0:
            self.log(f"Unset fixed calculator {self.calc1}.")
            self.fixed_calc = None

        # Avoid unnecessary costly Hessian calculations; decide which Hessian
        # to calculate from previous energy calculation.
        if name == "get_hessian":
            tmp_energy1 = run_calculation(self.calc1, "get_energy")["energy"]
            tmp_energy2 = run_calculation(self.calc2, "get_energy")["energy"]
            if tmp_energy1 <= tmp_energy2:
                results1 = run_calculation(self.calc1)
                results2 = {"energy": tmp_energy2}
                self.log("Calculated Hessian for calc1, skipped it for calc2.")
            else:
                results1 = {"energy": tmp_energy1}
                results2 = run_calculation(self.calc2)
                self.log("Calculated Hessian for calc2, skipped it for calc1.")
        # Do both full calculation otherwise
        else:
            results1 = run_calculation(self.calc1, name)
            results2 = run_calculation(self.calc2, name)
        energy1 = results1["energy"]
        energy2 = results2["energy"]
        all_energies = np.array((energy1, energy2))

        # Mixed forces to optimize crossing points
        if self.mix:
            # Must be positive, so substract lower energy from higher energy.
            # I is the higher state, j the lower one.
            en_i, en_j = (
                (energy2, energy1) if (energy1 < energy2) else (energy1, energy2)
            )
            energy_diff = en_i - en_j
            energy_diff_kJ = energy_diff * AU2KJPERMOL
            assert energy_diff > 0.0
            self.log(
                f"Mix mode, ΔE={energy_diff:.6f} au ({energy_diff_kJ:.2f} kJ mol⁻¹)"
            )
            energy_diff_sq = energy_diff**2
            denom = energy_diff + self.alpha
            energy = (en_i + en_j) / 2 + self.sigma * energy_diff_sq / denom
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
                # There seems to be a typo in Eq. (8) in [1]; the final term
                # should be (dV_J - dV_I) instead of the (dV_I - dV_J).
                # The formula is correct though, in the original publication
                # (Eq. (7) in [2]).
                forces = (forces_i + forces_j) / 2 + self.sigma * (
                    energy_diff_sq + 2 * self.alpha * energy_diff
                ) / denom**2 * (forces_i - forces_j)
                results["forces"] = forces
                for key, forces in (
                    ("mixed forces", forces),
                    ("forces1", forces1),
                    ("forces2", forces2),
                ):
                    norm, max_, rms_ = norm_max_rms(forces)
                    self.log(
                        f"{key: >14s}: (norm={norm:.4f}, max(|{key: >14s}|)={max_:.4f}, "
                        f"rms={rms_:.4f}) au a0⁻¹"
                    )
            self.calc_counter += 1
            return results
        # Mixed forces end

        min_ind = [1, 0][int(energy1 < energy2)]
        en1_or_en2 = ("calc1", "calc2")[min_ind]
        energy_diff = energy1 - energy2
        # Try to fix calculator, if requested
        if (self.min_energy_diff and self.check_after) and (
            # When the actual difference is above to minimum differences
            # or
            # no calculator is fixed yet
            (energy_diff > self.min_energy_diff)
            or (self.fixed_calc is None)
        ):
            self.fixed_calc = (self.calc1, self.calc2)[min_ind]
            self.recalc_in = self.check_after
            self.log(
                f"Fixed calculator choice ({en1_or_en2}) for {self.recalc_in} cycles."
            )
        results = (results1, results2)[min_ind]
        results["all_energies"] = all_energies
        energy_diff_kJ = abs(energy_diff) * AU2KJPERMOL

        self.log(
            f"energy_calc1={energy1:.6f} au, energy_calc2={energy2:.6f} au, returning "
            f"results for {en1_or_en2}, {energy_diff_kJ: >10.2f} kJ mol⁻¹ lower."
        )
        self.calc_counter += 1
        return results

    def get_energy(self, atoms, coords, **prepare_kwargs):
        return self.do_calculations("get_energy", atoms, coords, **prepare_kwargs)

    def get_forces(self, atoms, coords, **prepare_kwargs):
        return self.do_calculations("get_forces", atoms, coords, **prepare_kwargs)

    def get_hessian(self, atoms, coords, **prepare_kwargs):
        return self.do_calculations("get_hessian", atoms, coords, **prepare_kwargs)

    def get_chkfiles(self) -> dict:
        chkfiles = {}
        for key in ("calc1", "calc2"):
            try:
                chkfiles[key] = getattr(self, key).get_chkfiles()
            except AttributeError:
                pass
        return chkfiles

    def set_chkfiles(self, chkfiles: dict):
        for key in ("calc1", "calc2"):
            try:
                getattr(self, key).set_chkfiles(chkfiles[key])
                msg = f"Set chkfile on '{key}'"
            except KeyError:
                msg = f"Found no information for '{key}' in chkfiles!"
            except AttributeError:
                msg = f"Setting chkfiles is not supported by '{key}'"
            self.log(msg)
