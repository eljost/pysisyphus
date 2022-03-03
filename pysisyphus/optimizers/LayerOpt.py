import numpy as np

from pysisyphus.calculators import IPIServer
from pysisyphus.Geometry import Geometry
from pysisyphus.helpers_pure import full_expand
from pysisyphus.optimizers.Optimizer import Optimizer
from pysisyphus.optimizers.LBFGS import LBFGS
from pysisyphus.optimizers.PreconLBFGS import PreconLBFGS
from pysisyphus.optimizers.RFOptimizer import RFOptimizer

"""
Require all layers to be defined?
TODO: allow setting user-chosen optimizer.
"""


class LayerOpt(Optimizer):
    def __init__(
        self,
        geometry,
        layers,
        **kwargs,
    ):
        super().__init__(geometry, **kwargs)
        assert geometry.coord_type == "cart"

        self.layers = layers

        atoms = geometry.atoms
        all_indices = np.arange(len(geometry.atoms))
        freeze_mask = np.full_like(all_indices, True, dtype=bool)

        self.layer_indices = list()
        self.layer_get_geoms = list()
        self.layer_get_opts = list()

        # Iterate in reverse order from smallest (lowest) layer to biggest (highest) layer.
        indices_below = set()
        self.handle = open("model.trj", "w")
        for i, layer in enumerate(layers[::-1]):
            try:
                indices = full_expand(layer["indices"])
            # Allow missing indices. Then the indices for the whole system are assumed.
            except KeyError:
                indices = all_indices
            # Drop indices from layer layers
            indices = sorted(set(indices) - indices_below)
            indices_below.update(set(indices))
            self.layer_indices.append(indices)
            layer_mask = freeze_mask.copy()
            # Don't freeze current layer
            layer_mask[indices] = False

            try:
                address = layer["address"]
                calc = IPIServer(address=address)
                if i == 0:
                    self.geometry.set_calculator(calc)
            except KeyError as err:
                print("Currently, a socket address is mandatory!")
                raise err

            # Geometry
            def get_geom_getter():
                coord_type = "redund" if (i == 0) else "cart"
                freeze_atoms = all_indices[layer_mask] if i != 0 else None
                layer_calc = calc
                coord_kwargs = (
                    {
                        "define_for": indices,
                    }
                    if i == 0
                    else {}
                )
                # freeze_atoms = None
                # coord_type = "redund"
                # coord_kwargs = {"define_for": indices}

                def get_geom(coords3d):
                    geom = Geometry(
                        atoms,
                        coords3d.copy(),
                        freeze_atoms=freeze_atoms,
                        coord_type=coord_type,
                        coord_kwargs=coord_kwargs,
                    )
                    geom.set_calculator(layer_calc)
                    return geom

                return get_geom

            get_geom = get_geom_getter()
            if i == 0:
                geom = get_geom(self.geometry.coords3d)

                def get_geom(coords3d):
                    self.handle.write(geom.as_xyz() + "\n")
                    self.handle.flush()
                    geom.coords3d = coords3d
                    return geom

            self.layer_get_geoms.append(get_geom)

            # Optimizer
            def get_opt_getter():
                layer_indices = indices
                rms_force_thresh = self.convergence["rms_force_thresh"]
                key = f"layer{i}"
                _opt_kwargs = {
                    "prefix": key,
                    "max_cycles": 150,
                    # "thresh": "gau",
                    "thresh": self.thresh,
                    "line_search": True,
                    "align": False,
                    "overachieve_factor": 2,
                    "keep_last": 15,
                    # "max_cycles": 50,
                    # "mu_reg": 0.1,
                }
                # opt_kwargs = {
                # "max_cycles": 150,
                # "thresh": self.thresh,
                # "overachieve_factor": 2,
                # "hessian_init": "unit",
                # }
                def get_opt(geom, forces=None):
                    opt_kwargs = _opt_kwargs.copy()
                    # if forces is not None:
                        # forces_thresh = forces.reshape(-1, 3)[layer_mask].flatten()
                        # rms = np.sqrt(np.mean(forces_thresh ** 2))
                        # fact = 10
                        # if (rms / fact) > rms_force_thresh:
                            # opt_kwargs.update(
                                # {
                                    # "rms_force": 1 * rms_force_thresh,
                                    # "rms_force_only": True,
                                # }
                            # )
                    # else:
                        # opt_kwargs.update(
                            # {
                                # "rms_force": rms_force_thresh,
                                # "rms_force_only": True,
                            # }
                        # )
                    opt = LBFGS(geom, **opt_kwargs)
                    # opt = RFOptimizer(geom, **opt_kwargs)
                    return opt

                return get_opt

            get_opt = get_opt_getter()
            if i == 0:
                model_opt = RFOptimizer(geom, thresh="never")
                model_opt.prepare_opt()  # TODO: do this outside of constructor

                def get_opt(geom, forces):
                    return model_opt

            self.layer_get_opts.append(get_opt)

        self.layer_indices = self.layer_indices[::-1]
        self.layer_get_geoms = self.layer_get_geoms[::-1]
        self.layer_get_opts = self.layer_get_opts[::-1]

    @property
    def layer_num(self):
        return len(self.layers)

    def optimize(self):
        coords3d_org = self.geometry.coords3d.copy()
        coords3d_cur = coords3d_org.copy()
        for i, (indices, get_geom, get_opt) in enumerate(
            zip(self.layer_indices, self.layer_get_geoms, self.layer_get_opts)
        ):
            geom = get_geom(coords3d_cur)
            get_opt_kwargs = {
                "forces": self.forces[-1] if self.cur_cycle > 0 else None,
            }
            opt = get_opt(geom, **get_opt_kwargs)
            if i == self.layer_num - 1:
                break
            opt.run()
            coords3d_cur[indices] = geom.coords3d[indices]
            # geom.jmol()

        ####################
        # Relax last layer #
        ####################

        # 'geom' and 'indices' for the last layer were defined in the for-loop
        # above, before breaking out.
        cart_forces = geom.cart_forces
        # print(cart_forces.reshape(-1, 3))
        energy = geom.energy
        self.energies.append(energy)
        self.forces.append(cart_forces.copy())

        opt.coords.append(geom.coords.copy())
        opt.cart_coords.append(geom.cart_coords.copy())

        # Calculate one step
        int_step = opt.optimize()
        opt.steps.append(int_step)
        geom.coords = geom.coords + int_step
        coords3d_cur[indices] = geom.coords3d[indices]

        full_step = coords3d_cur - coords3d_org
        return full_step.flatten()
