# [1] http://dx.doi.org/10.1021/acs.jctc.8b00885
#     Exploring Potential Energy Surface with External Forces
#     Wolinski, 2018

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

    def set_ref_hessian_for_geom(self, ref_geom):
        # Calculate actual Hessian with internal calculator
        results = self.calculator.get_hessian(ref_geom.atoms, ref_geom.cart_coords)
        self._ref_hessian = results["hessian"]
        Hr0 = self._ref_hessian @ self.ref_geom.cart_coords[:, None]
        # As shown in left column on p. 6309 of [1]
        self.s = self.max_force / np.abs(Hr0).max()
        self.log(f"Set EGO s={self._s:.6f}")

    @property
    def ref_hessian(self):
        if self._ref_hessian is None:
            self.set_ref_hessian_for_geom(self.ref_geom)
        return self._ref_hessian

    @property
    def s(self):
        return self._s

    @s.setter
    def s(self, s):
        self._s = s

    def get_external_forces(self, atoms, coords):
        assert atoms == self.ref_geom.atoms
        Hr = self.ref_hessian @ coords[:, None]
        energy_mod = float(coords[None, :] @ Hr)
        forces_mod = -self.s * Hr.flatten()
        return energy_mod, forces_mod

    def get_energy(self, atoms, coords, **prepare_kwargs):
        true_energy = self.calculator.get_energy(atoms, coords)["energy"]
        energy_mod, _ = self.get_external_forces(atoms, coords)

        results = {
            "energy": true_energy + energy_mod,
            "true_energy": true_energy,
        }
        # Manually increase calculation counter of EGO calculator
        self.calc_counter += 1
        return results

    def get_forces(self, atoms, coords, **prepare_kwargs):
        true_results = self.calculator.get_forces(atoms, coords, **prepare_kwargs)
        true_energy = true_results["energy"]
        true_forces = true_results["forces"]
        energy_mod, forces_mod = self.get_external_forces(atoms, coords)

        results = {
            "energy": true_energy + energy_mod,
            "forces": true_forces + forces_mod,
            "true_forces": true_forces,
            "true_energy": true_energy,
        }
        # Manually increase calculation counter of EGO calculator
        self.calc_counter += 1
        return results

    def get_hessian(self, atoms, coords, **prepare_kwargs):
        raise NotImplementedError("EGO-Hessian is not implemented!")
