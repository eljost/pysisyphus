import logging

import numpy as np

from pysisyphus.optimizers.Optimizer import Optimizer
from pysisyphus.optimizers import poly_fit


logger = logging.getLogger("optimizer")


class ONIOMOpt(Optimizer):

    def __init__(self, geometry, *args, micro_cycles=None, **kwargs):
        super().__init__(geometry, *args, **kwargs)

        layers = self.geometry.layers
        print(f"found {len(layers)} layers: {layers}")

        if micro_cycles is None:
            micro_cycles = np.ones(len(layers), dtype=int)
            micro_cycles[0] = 5
        self.micro_cycles = micro_cycles
        self.log("Micro cycles: {self.micro_cycles}")

        self.calc = self.geometry.calculator
        self.layer_indices = [self.calc.atom_inds_in_layer(i, exclude_inner=True)
                              for i, _ in enumerate(layers)
        ]

        # Conjugate gradient, previous search direction
        self.prev_direction = None
        # Initial step length for line search
        self.trial_length = 0.1

    def cg_step(self, atoms, coords, index, beta_formula="HS"):
        energy, forces = self.calc.calc_layer(atoms, coords, index,
                                              parent_correction=False)
        self.energies.append(energy)
        self.forces.append(forces)

        prev_grad = -forces
        prev_energy = energy

        # Direction of steepst descent in the first cycle
        if self.prev_direction is None:
            self.prev_direction = forces

        norm = np.linalg.norm(self.prev_direction)
        for i in range(3):
            self.log(f"Linesearch with trial step length {self.trial_length:.4f}")
            trial_step = self.trial_length * self.prev_direction / norm
            trial_coords = coords + trial_step
            trial_result = self.calc.calc_layer(atoms, trial_coords, index,
                                                parent_correction=False)
            trial_energy, trial_forces = trial_result
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
                self.trial_length = np.linalg.norm(step)
                break
            else:
                self.trial_length *= 2
                self.log("Linesearch did not produced a result. Trying longer "
                         "trial step length.")
        else:
            raise Exception("Linesearchfailed")

        # Hestensen-Stiefel
        if beta_formula == "HS":
            beta = grad.dot(grad - prev_grad) / (grad - prev_grad).dot(self.prev_direction)
        # Polak-Ribiere
        elif beta_formula == "PR":
            beta = grad.dot(grad - prev_grad) / prev_grad.dot(prev_grad)
        else:
            raise Exception("Invalid 'beta_formula'. Use 'PR' or 'HS'!")
        beta = max(0, beta)
        self.log(f"beta = {beta:.4f}")
        self.prev_direction = -grad + beta*self.prev_direction
        return step

    def optimize(self):
        atoms = self.geometry.atoms
        coords = self.geometry.cart_coords

        # Microcycles:
        for i in range(self.micro_cycles[0]):
            print(i)

        return self.cg_step(atoms, coords, 0, beta_formula="HS")
