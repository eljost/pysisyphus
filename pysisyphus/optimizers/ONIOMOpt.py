import logging

import numpy as np

from pysisyphus.constants import AU2KJPERMOL
from pysisyphus.optimizers.Optimizer import Optimizer
from pysisyphus.optimizers import poly_fit


logger = logging.getLogger("optimizer")


class ONIOMOpt(Optimizer):

    def __init__(self, geometry, *args, micro_cycles=None, **kwargs):
        print("The ONIOMOpt optimizer is not really ready yet!")
        super().__init__(geometry, *args, **kwargs)

        layers = self.geometry.layers
        print(f"found {len(layers)} layers: {layers}")

        if micro_cycles is None:
            micro_cycles = np.ones(len(layers), dtype=int)
        try:
            micro_cycles[0] = 5
        except IndexError:
            micro_cycles = None
        self.micro_cycles = micro_cycles
        self.log(f"Micro cycles: {self.micro_cycles}")

        self.calc = self.geometry.calculator
        if len(layers) > 1:
            self.layer_indices = [self.calc.atom_inds_in_layer(i, exclude_inner=True)
                                  for i, _ in enumerate(layers)
            ]
        else:
            self.layer_indices = [[i for i, atom in enumerate(self.geometry.atoms)] ]

        # Conjugate gradient, previous search direction
        self.prev_directions = [None for layer in layers]
        # Initial step length for line search
        self.trial_lengths = [0.1 for layer in layers]

    def cg_step(self, atoms, coords, index, beta_formula="HS", full=False):
        if full:
            res = self.geometry.get_energy_and_forces_at(coords)
            forces = res["forces"]
            energy = res["energy"]
        else:
            energy, forces = self.calc.calc_layer(atoms, coords, index,
                                                  parent_correction=False)
        def stat(forces):
            f3d = forces.reshape(-1, 3)
            if not full:
                f3d = f3d[self.layer_indices[index]]
            # max_ = np.abs(forces).max()
            # rms_ = np.sqrt(np.mean(forces**2))
            max_ = np.abs(f3d).max()
            rms_ = np.sqrt(np.mean(f3d**2))
            self.log(f"\tStart: max={max_:.6f}, rms={rms_:.6f}")
        stat(forces)

        prev_grad = -forces
        prev_energy = energy

        # Direction of steepst descent in the first cycle
        prev_direction = self.prev_directions[index]
        if prev_direction is None:
            prev_direction = forces
        atom_indices = self.layer_indices[index]
        if not full:
            if atom_indices == [10, 11, 12, 13, 14]:
                atom_indices = [7] + atom_indices
            _ = np.zeros_like(prev_direction).reshape(-1, 3)
            _[atom_indices] = prev_direction.reshape(-1, 3)[atom_indices]
            prev_direction = _.flatten()

        trial_length = self.trial_lengths[index]

        norm = np.linalg.norm(prev_direction)
        for i in range(3):
            self.log(f"Linesearch with trial step length {trial_length:.6f}")
            trial_step = trial_length * prev_direction / norm
            trial_coords = coords + trial_step
            if full:
                res = self.geometry.get_energy_and_forces_at(trial_coords)
                trial_forces = res["forces"]
                trial_energy = res["energy"]
            else:
                trial_energy, trial_forces = self.calc.calc_layer(atoms, trial_coords, index,
                                                                  parent_correction=False)
            ls_kwargs = {
                "cur_energy": trial_energy,
                "cur_grad": -trial_forces,
                "prev_energy": prev_energy,
                "prev_grad": prev_grad,
                "prev_step": trial_step,
                "prev_coords": coords,
                "allow_cubic": True,
                "allow_none": False,
                "cubic_max": 2.,
                "quartic_max": 4.,
            }
            ls_result = poly_fit.poly_line_search(**ls_kwargs)
            if ls_result[0] is not None:
                energy, grad, new_coords, step = ls_result
                trial_length = np.linalg.norm(step)
                break
            else:
                trial_length *= 2
                self.log("Linesearch did not produced a result. Trying longer "
                         "trial step length.")
        else:
            self.trial_lengths[index] = 0.1
            self.prev_directions[index] = forces
            step = forces
            norm = np.linalg.norm(step)
            if norm > 0.5:
                step = step/norm * 0.5
            self.log("Steepest descent FALLBACK")
            return step
            # ls_result = poly_fit.poly_line_search(**ls_kwargs)
            # raise Exception("Linesearchfailed")

        # Hestensen-Stiefel
        if beta_formula == "HS":
            beta = grad.dot(grad - prev_grad) / (grad - prev_grad).dot(prev_direction)
        # Polak-Ribiere
        elif beta_formula == "PR":
            beta = grad.dot(grad - prev_grad) / prev_grad.dot(prev_grad)
        else:
            raise Exception("Invalid 'beta_formula'. Use 'PR' or 'HS'!")
        beta = max(0, beta)
        self.log(f"beta = {beta:.4f}")
        self.prev_directions[index] = -grad + beta*prev_direction
        self.trial_lengths[index] = trial_length

        return step

    def sd_step(self, atoms, coords, index, full=True):
        if full:
            res = self.geometry.get_energy_and_forces_at(coords)
            forces = res["forces"]
            energy = res["energy"]
        else:
            energy, forces = self.calc.calc_layer(atoms, coords, index,
                                                  parent_correction=False)
        def stat(forces):
            f3d = forces.reshape(-1, 3)
            f3d = f3d[self.layer_indices[index]]
            # max_ = np.abs(forces).max()
            # rms_ = np.sqrt(np.mean(forces**2))
            max_ = np.abs(f3d).max()
            rms_ = np.sqrt(np.mean(f3d**2))
            self.log(f"\tStart: max={max_:.6f}, rms={rms_:.6f}")
        stat(forces)

        prev_grad = -forces
        prev_energy = energy

        # Steepest descent direction
        prev_direction = forces
        # atom_indices = self.layer_indices[index]
        # if atom_indices == [10, 11, 12, 13, 14]:
            # atom_indices = [7] + atom_indices
        # self.log(f"\tatom_indices={atom_indices}")
        # _ = np.zeros_like(prev_direction).reshape(-1, 3)
        # _[atom_indices] = prev_direction.reshape(-1, 3)[atom_indices]
        # prev_direction = _.flatten()

        trial_length = self.trial_lengths[index]

        norm = np.linalg.norm(prev_direction)
        for i in range(3):
            self.log(f"\tLinesearch with trial step length {trial_length:.6f}")
            trial_step = trial_length * prev_direction / norm
            trial_coords = coords + trial_step
            res = self.geometry.get_energy_and_forces_at(trial_coords)
            if full:
                res = self.geometry.get_energy_and_forces_at(trial_coords)
                trial_forces = res["forces"]
                trial_energy = res["energy"]
            else:
                trial_forces = res["forces"]
                trial_energy = res["energy"]
            ls_kwargs = {
                "cur_energy": trial_energy,
                "cur_grad": -trial_forces,
                "prev_energy": prev_energy,
                "prev_grad": prev_grad,
                "prev_step": trial_step,
                "prev_coords": coords,
                "allow_cubic": True,
                "allow_none": False,
                "cubic_max": 2.,
                "quartic_max": 4.,
            }
            ls_result = poly_fit.poly_line_search(**ls_kwargs)
            if ls_result[0] is not None:
                energy, grad, new_coords, step = ls_result
                trial_length = np.linalg.norm(step)
                break
            else:
                trial_length *= 2
                self.log("Linesearch did not produced a result. Trying longer "
                         "trial step length.")
        else:
            # Reset everything
            self.trial_lengths[index] = 0.1
            self.prev_directions[index] = forces
            step = forces
            self.log("Steepest descent FALLBACK")

        self.trial_lengths[index] = trial_length

        return step

    # def optimize(self):
        # atoms = self.geometry.atoms

        # forces = self.geometry.forces
        # energy = self.geometry.energy
        # self.forces.append(forces)
        # self.energies.append(energy)

        # if self.cur_cycle > 0:
            # self.log(f"Current energy={energy:.6f}")
            # dE = energy - self.energies[-2]
            # dE_str = "raised" if dE >= 0 else "lowered"
            # dEkj = dE*AU2KJPERMOL
            # self.log(f"Current energy: {energy:.6f}, energy {dE_str} "
                     # f"by {dEkj:.2f} kJ mol⁻¹")
            # if dE_str == "raised":
                # print("Raised!")

        # # Calling copy is important, otherwise we would modify the geomtries
        # # coordinates.
        # coords3d = self.geometry.coords3d.copy()
        # li0 = self.layer_indices[0]
        # li1 = self.layer_indices[1]
        # coords3d_1 = coords3d[li1].copy()

        # # step_func = self.cg_step
        # step_func = self.sd_step

        # layer = 0
        # l0_step = np.zeros_like(coords3d)
        # # Microcycles
        # self.log(f"Starting microcycles for layer {layer}")
        # for i in range(self.micro_cycles[0]):
            # self.log(f"Microcycle {i}")
            # step = step_func(atoms, coords3d.flatten(), layer)
            # # Only set coordinates of atoms in layer 0
            # _ = l0_step.copy()
            # _[li0] = step.reshape(-1, 3)[li0]
            # coords3d += _
        # np.testing.assert_allclose(coords3d[li1], coords3d_1)

        # # # Step for inner layer
        # layer = 1
        # self.log(f"\n\nStarting cycle for inner layer {layer}")
        # # step = self.cg_step(atoms, coords3d.flatten(), layer, full=True)
        # step = step_func(atoms, coords3d.flatten(), layer)
        # # print(step.reshape(-1,3))
        # coords3d += step.reshape(-1,3)

        # step = (coords3d - self.geometry.coords3d).flatten()
        # return step

    def optimize(self):
        atoms = self.geometry.atoms
        coords = self.geometry.coords
        forces = self.geometry.forces
        energy = self.geometry.energy
        self.forces.append(forces)
        self.energies.append(energy)
        if self.cur_cycle > 0:
            self.log(f"Current energy={energy:.6f}")
            dE = energy - self.energies[-2]
            dE_str = "raised" if dE >= 0 else "lowered"
            dEkj = dE*AU2KJPERMOL
            self.log(f"Current energy: {energy:.6f}, energy {dE_str} "
                     f"by {dEkj:.2f} kJ mol⁻¹")
            if dE_str == "raised":
                print("Raised!")
        return self.cg_step(atoms, coords, 0, full=True, beta_formula="PR")
        # return self.sd_step(atoms, coords, 0, full=True)
