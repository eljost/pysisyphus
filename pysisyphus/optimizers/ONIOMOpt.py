import numpy as np

from pysisyphus.optimizers.Optimizer import Optimizer
from pysisyphus.optimizers.LBFGS import LBFGS
from pysisyphus.optimizers.closures import small_lbfgs_closure
from pysisyphus.optimizers.restrict_step import scale_by_max_step
from pysisyphus.helpers_pure import highlight_text


class ONIOMOpt(Optimizer):
    def __init__(
        self,
        geometry,
        *args,
        micro_cycles=None,
        max_micro_cycles=50,
        control_step=False,
        max_step=0.2,
        step="full",
        **kwargs,
    ):
        super().__init__(geometry, max_step=max_step, **kwargs)

        self.max_micro_cycles = max_micro_cycles
        self.control_step = control_step
        self.step = step
        assert self.step in ("full", "high")

        layers = self.geometry.layers
        assert len(layers) == 2, "Only ONIOM2 supported yet!"

        # Set up micro cycles for every layer
        if micro_cycles is None:
            micro_cycles = np.ones(len(layers), dtype=int)
            micro_cycles[0] = 0
        self.micro_cycles = micro_cycles
        self.log(f"Micro cycles: {self.micro_cycles}")

        self.lbfgs_closure = small_lbfgs_closure(history=10, gamma_mult=False)
        self.high_steps = list()

        (self.real_model,), (self.high_model,) = self.geometry.layers
        # Freeze high layer, but also freeze atoms in real layer that become
        # link atoms in the high layer.
        link_inds = [link.parent_ind for link in self.high_model.links]
        self.freeze_in_real = self.high_model.atom_inds + link_inds
        self.freeze_in_high = [2, 1, 3]
        print("!!! HARDCODED self.freeze_in_high !!!")

        self.micros_converged = 0
        self.micro_opt_cycles = list()

    def optimize(self):

        #######################
        # Relax real geometry #
        #######################

        coords3d_org = self.geometry.coords3d.copy()
        real_geom = self.real_model.as_geom(self.geometry.atoms, coords3d_org.copy())
        real_geom.freeze_atoms = self.freeze_in_real

        key = "real"
        micro_cycles = self.micro_cycles[0]
        if micro_cycles == 0:
            micro_cycles = self.max_micro_cycles
        real_opt_kwargs = {
            "prefix": key,
            "h5_group_name": f"{key}_opt",
            "max_cycles": micro_cycles,
            "thresh": self.thresh,  # Respect parents convergence threshold
            "line_search": True,
            "align": False,
        }
        real_opt = LBFGS(real_geom, **real_opt_kwargs)
        print("\n" + highlight_text(f"Opt Cycle {self.cur_cycle}, Micro cycles") + "\n")
        real_opt.run()
        real_step = real_geom.coords3d - coords3d_org
        self.micros_converged += real_opt.is_converged
        self.micro_opt_cycles.append(real_opt.cur_cycle + 1)
        print("\n" + highlight_text(f"Micro cycles finished") + "\n")

        #######################
        # Relax full geometry #
        #######################

        # Calculate full ONIOM forces with previously releaxed, real coordinates
        results = self.geometry.get_energy_and_forces_at(real_geom.coords3d.flatten())
        energy = results["energy"]
        forces = results["forces"]
        self.energies.append(energy)
        self.forces.append(forces.copy())

        try:
            prev_step = self.steps[-1] + real_step.flatten()
        except IndexError:
            prev_step = None
        if self.step == "high":
            forces.reshape(-1, 3)[self.freeze_in_high] = 0.0
        step = self.lbfgs_closure(forces, prev_step=prev_step)

        if self.control_step:
            step = scale_by_max_step(step, self.max_step)
        step += real_step.flatten()
        return step

    def postprocess_opt(self):
        tot_micro_cycs = sum(self.micro_opt_cycles)
        msg = (
            f"\nMicro-cycle optimizations:\n"
            f"\t        Attempted: {self.cur_cycle+1}\n"
            f"\t        Converged: {self.micros_converged}\n"
            f"\t     Total cycles: {tot_micro_cycs}\n"
        )
        self.log(msg)
