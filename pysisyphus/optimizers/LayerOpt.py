# [1] https://doi.org/10.1002/jcc.10156
#     Geometry optimization with QM/MM, ONIOM, and other combined methods.
#     I. Microiterations and constraints
#     Vreven, Morokuma, Farkas, Schlegel, Frisch, 2003
# [2] https://doi.org/10.1039/A909486E
#     Linear scaling geometry optimisation and transition state search
#     in hybrid delocalised internal coordinates
#     Billeter, Turner, Thiel, 2000
import numpy as np

from pysisyphus.calculators import IPIServer
from pysisyphus.Geometry import Geometry
from pysisyphus.helpers_pure import full_expand
from pysisyphus.optimizers.Optimizer import Optimizer
from pysisyphus.optimizers.LBFGS import LBFGS
from pysisyphus.optimizers.RFOptimizer import RFOptimizer

"""
Require all layers to be defined?
TODO: allow setting user-chosen optimizer.
"""


class Layers:
    def __init__(self, geometry, opt_thresh, layers=None):
        self.geometry = geometry
        self.layers = layers
        self.opt_thresh = opt_thresh

        atoms = geometry.atoms
        all_indices = np.arange(len(geometry.atoms))
        # Boolean array; every item corresponds to an atom in the total system.
        freeze_mask = np.full_like(all_indices, True, dtype=bool)

        self.indices = list()
        self.geom_getters = list()
        self.opt_getters = list()

        indices_below = set()
        # Iterate in reverse order from smallest (lowest) layer to biggest (highest) layer.
        # to setup geom_getter and opt_getter.
        for i, layer in enumerate(layers[::-1]):
            try:
                indices = full_expand(layer["indices"])
            # Allow missing 'indices' key. Then it is assumed that this layer contains the
            # whole system.
            except KeyError:
                assert i != 0, "Found whole system in layer 0. I don't like that!"
                indices = all_indices
            # Drop indices from layers below ...
            indices = sorted(set(indices) - indices_below)
            # ... and update 'indices_below' for the next layer
            indices_below.update(set(indices))
            self.indices.append(indices)
            # Set up mask that indicates which atoms to freeze in the current layer (
            # all atoms that are not in the current layer.)
            layer_mask = freeze_mask.copy()
            layer_mask[indices] = False

            try:
                address = layer["address"]
                calc = IPIServer(address=address)
                # When calculating layer 0 we have access to the true energy of the system.
                # So we assign the calculator of layer0 to the actual geometry containing
                # the whole system.
                if i == 0:
                    geometry.set_calculator(calc)
            # If no address is given, we assume that pysisyphus' ONIOM calculator
            # is used.
            except KeyError as err:
                try:
                    calc = layer["layer_calc"]
                except KeyError:
                    print(
                        "Currently, a socket address for an IPI-protol client is mandatory!"
                    )
                    raise err

            # Geometry
            def get_geom_getter():
                # coord_type = "redund" if (i == 0) else "cart"
                coord_type = "redund" if (i == 0) else "cartesian"
                # Don't freeze anything in layer 0; but use only the
                # mobile atoms to define internal coordinates.
                freeze_atoms = all_indices[layer_mask] if i != 0 else None
                coord_kwargs = (
                    {
                        "define_for": indices,
                    }
                    if i == 0
                    else {}  # no coord_kwargs for higher layers
                )
                layer_calc = calc

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
            # Use a persistent Geometry for the layer 0. Overwrite the function
            # above with a definition that always returns the same geometry.
            if i == 0:
                geom = get_geom(geometry.coords3d)

                def get_geom(coords3d):
                    geom.coords3d = coords3d
                    return geom

            self.geom_getters.append(get_geom)

            # Optimizer
            def get_opt_getter():
                _opt_kwargs = {
                    "max_cycles": 250,
                    "thresh": self.opt_thresh,
                    "overachieve_factor": 2,
                    # LBFGS
                    "line_search": True,
                    "keep_last": 15,
                    "align": False,
                    # RFO
                    # "hessian_init": "unit",
                }

                def get_opt(geom, forces=None):
                    opt_kwargs = _opt_kwargs.copy()
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

            self.opt_getters.append(get_opt)

        self.indices = self.indices[::-1]
        self.geom_getters = self.geom_getters[::-1]
        self.opt_getters = self.opt_getters[::-1]

    @classmethod
    def from_oniom_calculator(cls, geometry, oniom_calc=None, **kwargs):
        calc = geometry.calculator
        if calc is None:
            calc = oniom_calc
        layers = list()
        for i, layer in enumerate(calc.layers):
            assert len(layer) == 1, "Multicenter-ONIOM is not yet supported!"
            model = layer[0]
            link_hosts = [link.parent_ind for link in model.links]
            indices = model.atom_inds + link_hosts
            layer_calc = calc.get_layer_calc(i)
            layer = {
                    "layer_calc": layer_calc,
                    "indices": indices,
            }
            layers.append(layer)
        return Layers(geometry, layers=layers, **kwargs)

    def __len__(self):
        return len(self.layers)


class LayerOpt(Optimizer):
    def __init__(
        self,
        geometry,
        layers=None,
        **kwargs,
    ):
        super().__init__(geometry, **kwargs)
        assert geometry.coord_type == "cart"

        layers_kwargs = {
            "geometry": self.geometry,
            "opt_thresh": self.thresh,
        }
        if layers is not None:
            layers_kwargs["layers"] = layers
            layers = Layers(**layers_kwargs)
        # Construct layers from ONIOM calculator
        else:
            layers = Layers.from_oniom_calculator(**layers_kwargs)
        self.layers = layers

    @property
    def layer_num(self):
        return len(self.layers)

    def optimize(self):
        coords3d_org = self.geometry.coords3d.copy()
        coords3d_cur = coords3d_org.copy()
        for i, (indices, get_geom, get_opt) in enumerate(
            zip(self.layers.indices, self.layers.geom_getters, self.layers.opt_getters)
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

        ####################
        # Relax last layer #
        ####################

        # 'geom' and 'indices' for the last layer were defined in the for-loop
        # above, before breaking from the loop.
        #
        # Calculate forces and energy of the last layer. This have to be the "true"
        # ONIOM forces of the system, containing all contributions. That's why we
        # also save them as the true forces in the optimizer.
        cart_forces = geom.cart_forces
        energy = geom.energy
        self.energies.append(energy)
        self.forces.append(cart_forces.copy())

        # Also store relevant quantities in the optimizer of the last layer, so
        # stuff like Hessian updates are possible. These quantities are usually
        # stored in the big optimization-loop in Optimizer.run(). As run() is
        # never called for the last optimizer we have to store them manually.
        opt.coords.append(geom.coords.copy())
        opt.cart_coords.append(geom.cart_coords.copy())

        # Calculate one step
        int_step = opt.optimize()
        opt.steps.append(int_step)
        geom.coords = geom.coords + int_step
        coords3d_cur[indices] = geom.coords3d[indices]

        full_step = coords3d_cur - coords3d_org
        return full_step.flatten()
