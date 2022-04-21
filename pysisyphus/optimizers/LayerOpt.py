# [1] https://doi.org/10.1002/jcc.10156
#     Geometry optimization with QM/MM, ONIOM, and other combined methods.
#     I. Microiterations and constraints
#     Vreven, Morokuma, Farkas, Schlegel, Frisch, 2003
# [2] https://doi.org/10.1039/A909486E
#     Linear scaling geometry optimisation and transition state search
#     in hybrid delocalised internal coordinates
#     Billeter, Turner, Thiel, 2000
import logging
from pprint import pprint

import numpy as np

from pysisyphus.calculators import IPIServer, ONIOM
from pysisyphus.Geometry import Geometry
from pysisyphus.helpers_pure import full_expand, highlight_text, recursive_update
from pysisyphus.optimizers.Optimizer import Optimizer


def get_geom_kwargs(layer_ind, layer_mask):
    coord_type = "tric" if layer_ind == 0 else "cartesian"
    all_indices = np.arange(layer_mask.size)
    freeze_atoms = all_indices[layer_mask]
    geom_kwargs = {
        "type": coord_type,
        "freeze_atoms": freeze_atoms,
        "coord_kwargs": {
            "freeze_atoms_exclude": layer_ind == 0,
        },
    }
    return geom_kwargs


def get_opt_kwargs(opt_key, layer_ind, thresh):
    # Some defaults tailored to LayerOpt
    opt_defaults = {
        "lbfgs": {
            "mu_reg": 0.1,
        },
        "plbfgs": {
            "precon_kind": "full_fast",
            "precon_update": 50,
        },
    }
    if layer_ind == 0:
        opt_kwargs = {
            "type": "rfo",
            "thresh": "never",
        }
    else:
        if opt_key is None:
            opt_key = "lbfgs"

        opt_kwargs = {
            "type": opt_key,
            "max_cycles": 250,
            "thresh": thresh,
            "overachieve_factor": 5,
            "prefix": f"layer_{layer_ind:02d}",
        }
        try:
            opt_kwargs.update(opt_defaults[opt_key])
        except KeyError:
            pass
    return opt_kwargs


logger = logging.getLogger("optimizer")


class Layers:
    def __init__(self, geometry, opt_thresh, layers=None):
        self.geometry = geometry
        self.layers = layers
        self.opt_thresh = opt_thresh

        print(highlight_text("Layers", level=1))
        pprint(self.layers, compact=True)
        print("")

        atoms = geometry.atoms
        all_indices = np.arange(len(geometry.atoms))
        # Boolean array; every item corresponds to an atom in the total system.
        freeze_mask = np.full_like(all_indices, True, dtype=bool)

        self.indices = list()
        self.geom_getters = list()
        self.opt_getters = list()

        # We import 'get_opt_cls' in the constructor, as it relies on dicts, that
        # contain references to this class (LayerOpt). Trying to import 'get_opt_cls'
        # at the top of the module will result in an circular import.
        from pysisyphus.optimizers.cls_map import get_opt_cls

        indices_below = set()
        # Iterate in reverse order from smallest (lowest) layer to biggest (highest) layer.
        # to setup geom_getter and opt_getter.
        for i, layer in enumerate(layers[::-1]):
            try:
                indices = full_expand(layer["indices"])
            # Allow missing 'indices' key. Then it is assumed that this layer contains the
            # whole system.
            except KeyError:
                assert (
                    i != 0
                ), "Found whole system in highest level layer. I don't like that!"
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
                try:
                    calc_kwargs = {
                        "address": layer["address"],
                    }
                except KeyError:
                    _calc_kwargs = layer["calc"]
                    calc_kwargs = recursive_update({}, _calc_kwargs)
                    # The popped calc key is currently unused as using an IPIServer is
                    # mandatory.
                    _ = calc_kwargs.pop("type", None)
                calc = IPIServer(**calc_kwargs)
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

            ####################
            #  Geometry setup  #
            ####################

            _geom_kwargs = layer.get("geom", dict())
            _geom_kwargs = recursive_update({}, _geom_kwargs)
            geom_kwargs = get_geom_kwargs(i, layer_mask=layer_mask)
            try:
                geom_kwargs.update(_geom_kwargs)
            # Allow empty "geom:" block
            except TypeError:
                pass
            coord_type = geom_kwargs.pop("type")

            # Geometry
            def get_geom_getter():
                # Rebind the variables here, otherwise the wrong geom_kwargs
                # and calc will be used, as they are redefined in the next loop cycle.
                layer_geom_kwargs = geom_kwargs.copy()
                layer_calc = calc

                def get_geom(coords3d):
                    geom = Geometry(
                        atoms,
                        coords3d.copy(),
                        coord_type=coord_type,
                        **layer_geom_kwargs,
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

            #####################
            #  Optimizer setup  #
            #####################

            _opt_kwargs = layer.get("opt", dict())
            _opt_kwargs = recursive_update({}, _opt_kwargs)
            opt_key = _opt_kwargs.get("type", None)
            opt_kwargs = get_opt_kwargs(opt_key, i, thresh=self.opt_thresh)
            try:
                opt_kwargs.update(_opt_kwargs)
            # Allow empty "opt:" block
            except TypeError:
                pass
            opt_key = opt_kwargs.pop("type")
            opt_cls = get_opt_cls(opt_key)

            def get_opt_getter():
                layer_opt_kwargs = opt_kwargs
                layer_opt_cls = opt_cls

                def get_opt(geom):
                    opt = layer_opt_cls(geom, **layer_opt_kwargs)
                    return opt

                return get_opt

            get_opt = get_opt_getter()
            if i == 0:
                model_opt = get_opt(geom)

                def get_opt(geom):
                    return model_opt

            self.opt_getters.append(get_opt)

        self.indices = self.indices[::-1]
        self.geom_getters = self.geom_getters[::-1]
        self.opt_getters = self.opt_getters[::-1]

    @classmethod
    def from_oniom_calculator(cls, geometry, oniom_calc=None, layers=None, **kwargs):
        calc = geometry.calculator
        if calc is None:
            calc = oniom_calc

        if layers is not None:
            assert len(layers) == len(calc.layers), (
                f"ONIOM calculator has {len(calc.layers)} layers, but only "
                f"{len(layers)} layer were defined in 'layers:'!"
            )

            layers = list(layers)
            for i, layer in enumerate(layers):
                if layer is None:
                    layers[i] = dict()
                elif isinstance(layer, dict):
                    pass
                else:
                    raise Exception("Layer definition must be empty or dict-like!")
        else:
            layers = [dict() for _ in calc.layers]

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
            layers[i].update(layer)
        return Layers(geometry, layers=layers, **kwargs)

    def __len__(self):
        return len(self.layers)


class LayerOpt(Optimizer):
    def __init__(
        self,
        geometry: Geometry,
        layers: dict = None,
        **kwargs,
    ) -> None:
        super().__init__(geometry, **kwargs)
        assert geometry.coord_type in ("cart", "cartesian")

        layers_kwargs = {
            "geometry": self.geometry,
            "opt_thresh": self.thresh,
            "layers": layers,
        }
        # Construct layers from ONIOM calculator
        if isinstance(self.geometry.calculator, ONIOM):
            layers = Layers.from_oniom_calculator(**layers_kwargs)
        else:
            layers = Layers(**layers_kwargs)
        self.layers = layers

        self.micro_cycles = list()
        self.micro_cycles_converged = list()

    @property
    def layer_num(self) -> int:
        return len(self.layers)

    def optimize(self) -> None:
        coords3d_org = self.geometry.coords3d.copy()
        coords3d_cur = coords3d_org.copy()
        cur_micro_cycles = list()
        cur_micro_cycles_converged = list()
        for i, (indices, get_geom, get_opt) in enumerate(
            zip(self.layers.indices, self.layers.geom_getters, self.layers.opt_getters)
        ):
            print(highlight_text(f"Layer {i}", level=1))
            is_last_layer = i == self.layer_num - 1
            geom = get_geom(coords3d_cur)
            opt = get_opt(geom)
            if is_last_layer:
                if self.cur_cycle == 0:
                    opt.prepare_opt()
                break
            opt.run()
            cur_micro_cycles.append(opt.cur_cycle + 1)
            cur_micro_cycles_converged.append(opt.is_converged)
            coords3d_cur[indices] = geom.coords3d[indices]
        self.micro_cycles.append(cur_micro_cycles)
        self.micro_cycles_converged.append(cur_micro_cycles_converged)

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

    def postprocess_opt(self) -> None:
        coord_types = list()
        for layer in self.layers.layers:
            try:
                coord_type = layer["geom"]["type"]
                coord_types.append(coord_type)
            except KeyError:
                pass
        micro_sum = np.array(self.micro_cycles).sum()
        print("\nMicrocycles:")
        print("\t", end="")
        pprint(self.micro_cycles)
        print(f"\t@@@ Σ {micro_sum},", ", ".join(coord_types))
        print(f"\t@@@ Macrocycles: {self.cur_cycle+1}, converged? {self.is_converged}")
        print("\t@@@")
