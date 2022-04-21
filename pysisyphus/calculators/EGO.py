# [1] http://dx.doi.org/10.1021/acs.jctc.8b00885
#     Exploring Potential Energy Surface with External Forces

import numpy as np

from pysisyphus.calculators.Calculator import Calculator


class EGO(Calculator):
    def __init__(
        self,
        calculator,
        ref_geom,
        max_force=0.175,
        **kwargs,
    ):
        super().__init__(**kwargs)

        self.calculator = calculator
        self.ref_geom = ref_geom
        self.max_force = max_force
        assert (
            self.max_force > 0.0
        ), f"Maximum force must be positve bug to max_force={self.max_force:.4f} au"

        self._ref_hessian = None
        self._s = None

    @property
    def ref_hessian(self):
        if self._ref_hessian is None:
            geom = self.ref_geom
            results = self.calculator.get_hessian(geom.atoms, geom.cart_coords)
            self._ref_hessian = results["hessian"]
            Hr0 = self._ref_hessian @ self.ref_geom.cart_coords[:, None]
            self._s = self.max_force / np.abs(Hr0).max()
            self.log(f"Set EGO s={self._s:.6f}")
        return self._ref_hessian

    @property
    def s(self):
        return self._s

    def get_mods(self, atoms, coords):
        assert atoms == self.ref_geom.atoms
        Hr = self.ref_hessian @ coords[:, None]
        s = self.s
        energy_mod = float(coords[None, :] @ Hr)
        forces_mod = -s * Hr.flatten()
        return energy_mod, forces_mod

    def get_energy(self, atoms, coords, **prepare_kwargs):
        true_energy = self.calculator.get_energy(atoms, coords)["energy"]
        energy_mod, _ = self.get_mods(atoms, coords)

        results = {
            "energy": true_energy + energy_mod,
            "true_energy": true_energy,
        }
        self.calc_counter += 1
        return results

    def get_forces(self, atoms, coords, **prepare_kwargs):
        true_results = self.calculator.get_forces(atoms, coords, **prepare_kwargs)
        true_energy = true_results["energy"]
        true_forces = true_results["forces"]
        energy_mod, forces_mod = self.get_mods(atoms, coords)

        results = {
            "energy": true_energy + energy_mod,
            "forces": true_forces + forces_mod,
            "true_forces": true_forces,
            "true_energy": true_energy,
        }
        self.calc_counter += 1
        return results
