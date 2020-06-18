import numpy as np

from pysisyphus.optimizers.Optimizer import Optimizer
from pysisyphus.optimizers import poly_fit


def poly_line_search(cur_energy, prev_energy, cur_grad, prev_grad, prev_step, prev_coords):

    # Generate directional gradients by projecting them on the previous step.
    prev_grad_proj = prev_step @ prev_grad
    cur_grad_proj =  prev_step @ cur_grad
    cubic_result = poly_fit.cubic_fit(prev_energy, cur_energy,
                                          prev_grad_proj, cur_grad_proj)
    quartic_result = poly_fit.quartic_fit(prev_energy, cur_energy,
                                          prev_grad_proj, cur_grad_proj)
    accept = {
        "cubic": lambda x: (x > 0.) and (x <= 2),
        # "cubic": lambda x: False,
        "quartic": lambda x: (x > 0.) and (x <= 4),
    }

    fit_result = None
    if quartic_result and accept["quartic"](quartic_result.x):
        fit_result = quartic_result
        deg = "quartic"
    elif cubic_result and accept["cubic"](cubic_result.x):
        fit_result = cubic_result
        deg = "cubic"
    else:
        # Midpoint fallback as described by gaussian?
        if cur_energy < prev_energy:
            x = 1
            y = cur_energy
            deg = "current point"
        else:
            x = 0.5
            y = (cur_energy + prev_energy) / 2
            deg = "midpoint"
        fit_result = poly_fit.FitResult(x, y, polys=None)

    fit_energy = None
    fit_grad = None
    fit_coords = None
    fit_step = None
    if fit_result and fit_result.y < prev_energy:
        x = fit_result.x
        fit_energy = fit_result.y
        print(f"Did {deg} interpolation with x={x:.6f}.")
        # Interpolate coordinates and gradient
        fit_step = x * prev_step
        fit_coords = prev_coords + fit_step
        fit_grad = (1-x)*prev_grad + x*cur_grad
    return fit_energy, fit_grad, fit_coords, fit_step


class ONIOMOpt(Optimizer):

    def __init__(self, geometry, *args, max_micro_cycles=None, **kwargs):
        super().__init__(geometry, *args, **kwargs)

        layers = self.geometry.layers
        print(f"found {len(layers)} layers: {layers}")

        if max_micro_cycles is None:
            max_micro_cycles = np.ones(len(layers), dtype=int)
            max_micro_cycles[0] = 5
        self.max_micro_cycles = max_micro_cycles
        print("max_micro_cycles", self.max_micro_cycles)

        self.calc = self.geometry.calculator
        # l0 = calc.atom_inds_in_layer(0)
        # l0e = calc.atom_inds_in_layer(0, exclude_inner=True)
        # l1 = calc.atom_inds_in_layer(1)
        # l1e = calc.atom_inds_in_layer(1, exclude_inner=True)
        self.layer_indices = [self.calc.atom_inds_in_layer(i, exclude_inner=True)
                              for i, _ in enumerate(layers)
        ]
        # atoms = self.geometry.atoms
        # coords = self.geometry.cart_coords
        # calc.calc_layer(atoms, coords, 0)
        # calc.calc_layer(atoms, coords, 1)

        self.prev_direction = None
        self.trial_length = 0.1

    def cg_step(self, index):
        atoms, coords = self.geometry.atoms, self.geometry.cart_coords
        # res = self.calc.calc_layer(atoms, coords, index, parent_correction=False)
        energy, forces = self.calc.calc_layer(atoms, coords, index, parent_correction=False)

        # First cycle, steepest descent
        if self.prev_step is None:
            self.prev_direction = forces
            self.prev_energy = energy
            self.prev_forces = forces
            self.prev_coords = coords

            step = forces
            norm = np.linalg.norm(step)
            if norm > 0.5:
                direction = forces / np.linalg.norm(forces)
                step = 0.5 * direction
            self.prev_step = step
            self.geometry.coords = coords + step
            energy, forces = self.calc.calc_layer(atoms, coords, index, parent_correction=False)

        ls_kwargs = {
            "cur_energy": energy,
            "prev_energy": self.prev_energy,
            "cur_grad": -forces,
            "prev_grad": -self.prev_forces,
            "prev_step": self.prev_step,
            "prev_coords": self.prev_coords,
        }
        res = poly_line_search(**ls_kwargs)
        if res[0] is None:
            res = (energy, -forces, coords, self.prev_step)
        ienergy, igrad, icoords, istep = res
        beta = igrad.dot(igrad + self.prev_forces) / self.prev_forces.dot(self.prev_forces)
        print(f"beta = {beta:.4f}")
        self.prev_direction = -igrad + beta*self.prev_direction#self.prev_direc
        step = istep

        self.prev_energy = ienergy
        self.prev_forces = -igrad
        self.prev_step = step
        self.prev_coords = icoords.copy()

        return istep
        """
        Calc energy and forces at new geometry
        Linesearch
            Get new interpolated geometry, energy & gradient
        Calculate beta
        Calculate new step direction

        Store:
            prev_energy
            prev_forces
            prev_step
        """

    def cg_step(self, index):
        atoms, coords = self.geometry.atoms, self.geometry.cart_coords

        energy, forces = self.calc.calc_layer(atoms, coords, index,
                                              parent_correction=False)
        self.energies.append(energy)
        self.forces.append(forces)


        prev_grad = -forces
        prev_energy = energy

        if self.prev_direction is None:
            self.prev_direction = forces

        norm = np.linalg.norm(self.prev_direction)
        for i in range(3):
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
            }
            ls_result = poly_line_search(**ls_kwargs)
            if ls_result[0] is not None:
                energy, grad, new_coords, step = ls_result
                # self.trial_length = max(0.1, np.linalg.norm(step))
                self.trial_length = np.linalg.norm(step)
                print(f"\t trial_length={self.trial_length:.6f}")
                self.log("Haha")
                break
            else:
                # step = forces
                # snorm = np.linalg.norm(forces)
                # if snorm > 0.5:
                    # step = step/snorm*0.5
                # return step
                # import pdb; pdb.set_trace()
                # ls_result = poly_line_search(**ls_kwargs)
                self.trial_length *= 2
                print(f"\tHARHAR, new trial_length {self.trial_length:.4f}")
                continue
        else:
            raise Exception("Linesearchfailed")

        # beta = grad.dot(grad - prev_grad) / prev_grad.dot(prev_grad)  # PR
        beta = grad.dot(grad - prev_grad) / (grad - prev_grad).dot(self.prev_direction)  # HS
        print(f"beta {beta:.4f}")
        # self.log("beta", beta)
        # beta = max(0, beta)
        # # print(f"\talpha={alpha:.4f}, beta={beta:.4f}")
        self.prev_direction = -grad + beta*self.prev_direction
        # self.prev_direction = -grad
        return step

    def optimize(self):
        return self.cg_step(0)
        # forces = self.geometry.forces
        # energy = self.geometry.energy
        atoms, coords = self.geometry.atoms, self.geometry.cart_coords
        index = 0
        energy, forces = self.calc.calc_layer(atoms, coords, index, parent_correction=False)
        self.forces.append(forces)
        self.energies.append(energy)


        self.forces.append(forces)
        self.energies.append(energy)

        step = forces
        norm = np.linalg.norm(step)
        if norm > 0.5:
            direction = forces / np.linalg.norm(forces)
            step = 0.5 * direction
        return step
