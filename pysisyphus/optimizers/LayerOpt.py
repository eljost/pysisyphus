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


class LayerOpt(Optimizer):
    def __init__(
        self,
        geometry,
        layers,
        **kwargs,
    ):
        super().__init__(geometry, **kwargs)

        self.layers = layers

        atoms = geometry.atoms
        self.layer_indices = list()
        self.layer_calcs = list()
        self.layer_get_geoms = list()
        self.layer_get_opts = list()
        indices_above = set(range(len(self.geometry.atoms)))
        for i, layer in enumerate(layers):
            try:
                indices = full_expand(layer["indices"])
            except KeyError:
                indices = sorted(indices_above)
            try:
                next_indices = full_expand(layers[i + 1]["indices"])
            except IndexError:
                next_indices = None

            # Geometry
            def get_geom_getter():
                layer_indices = indices
                layer_atoms = [atoms[i] for i in layer_indices]
                coord_type = "cart"  # allow user choice
                freeze_atoms = next_indices

                # Calculator
                try:
                    address = layer["address"]
                    calc = IPIServer(address=address)
                except KeyError:
                    pass

                def get_geom(coords3d):
                    layer_coords = coords3d[layer_indices].flatten()
                    geom = Geometry(
                        layer_atoms,
                        layer_coords,
                        freeze_atoms=freeze_atoms,
                        coord_type=coord_type,
                    )
                    geom.set_calculator(calc)
                    return geom

                return get_geom

            geom_getter = get_geom_getter()
            if i == (len(layers) - 1):
                model_geom = geom_getter(self.geometry.coords3d)

                def geom_getter(coords3d):
                    return model_geom

            self.layer_get_geoms.append(geom_getter)

            # Optimizer
            def get_opt_getter():
                key = f"layer{i}"
                opt_kwargs = {
                    "prefix": key,
                    "h5_group_name": f"{key}_opt",
                    "max_cycles": 150,
                    "thresh": "gau",
                    "line_search": True,
                    "align": False,
                    "overachieve_factor": 2,
                }

                def get_opt(geom):
                    opt = LBFGS(geom, **opt_kwargs)
                    return opt

                return get_opt

            get_opt = get_opt_getter()
            if i == (len(layers) - 1):
                model_opt = RFOptimizer(
                    geom_getter(self.geometry.coords3d), thresh="never"
                )
                model_opt.prepare_opt()

                def get_opt(geom):
                    return model_opt

            self.layer_get_opts.append(get_opt)

        # We may not need this check because it could also be possible to only relax
        # an active center and not the whole molecule.
        # all_indices = set(it.chain(*self.layer_indices))
        # assert all_indices == set(range(len(self.geometry.atoms)))

    def optimize(self):
        # np.testing.assert_allclose(
        # self.layer_geoms[0].coords3d, self.geometry.coords3d[self.freeze_in_real]
        # )

        #####################################
        # Microiterations for real geometry #
        #####################################

        coords3d_org = self.geometry.coords3d.copy()
        for get_geom, get_opt in zip(self.layer_get_geoms, self.layer_get_opts):
            geom = get_geom(coords3d_org)
            opt = get_opt(geom)
            opt.run()
            break
        real_step = real_geom.coords3d - coords3d_org
        import pdb; pdb.set_trace()  # fmt: skip

        #######################
        # Relax full geometry #
        #######################

        # Calculate full ONIOM forces with previously releaxed, real coordinates
        results = self.geometry.get_energy_and_forces_at(real_geom.coords3d.flatten())
        energy = results["energy"]
        forces = results["forces"]
        self.energies.append(energy)
        self.forces.append(forces.copy())

        high_forces = forces.reshape(-1, 3)[self.freeze_in_real].flatten()
        self.model_geom._forces = high_forces
        self.model_geom._energy = energy
        model_coords3d_org = self.model_geom.coords3d.copy()
        self.model_opt.coords.append(self.model_geom.coords.copy())
        self.model_opt.cart_coords.append(self.model_geom.cart_coords.copy())
        # Calculate one step
        int_step = self.model_opt.optimize()
        self.model_opt.steps.append(int_step)
        self.model_geom.coords = self.model_geom.coords + int_step
        model_step = self.model_geom.coords3d - model_coords3d_org

        cart_step = real_step
        cart_step.reshape(-1, 3)[self.freeze_in_real] = model_step

        return cart_step.flatten()
