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
        assert geometry.coord_type == "cart"

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
            self.layer_indices.append(indices)
            try:
                next_indices = full_expand(layers[i + 1]["indices"])
            except IndexError:
                next_indices = None

            # Calculator
            try:
                address = layer["address"]
                calc = IPIServer(address=address)
            except KeyError:
                pass

            # Geometry
            def get_geom_getter():
                layer_indices = indices
                layer_atoms = [atoms[i] for i in layer_indices]
                freeze_atoms = next_indices
                layer_calc = calc

                def get_geom(coords3d, coord_type="cart"):
                    layer_coords = coords3d[layer_indices].flatten()
                    geom = Geometry(
                        layer_atoms,
                        layer_coords,
                        freeze_atoms=freeze_atoms,
                        coord_type=coord_type,
                    )
                    geom.set_calculator(layer_calc)
                    return geom

                return get_geom

            geom_getter = get_geom_getter()
            if i == (len(layers) - 1):
                model_geom = geom_getter(self.geometry.coords3d, "redund")

                def geom_getter(coords3d):
                    return model_geom

            self.layer_get_geoms.append(geom_getter)

            # Optimizer
            def get_opt_getter():
                key = f"layer{i}"
                opt_kwargs = {
                    "prefix": key,
                    "max_cycles": 150,
                    "thresh": self.thresh,
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
                model_opt.prepare_opt()  # TODO: do this outside of constructor

                def get_opt(geom):
                    return model_opt

            self.layer_get_opts.append(get_opt)
        self.geometry.set_calculator(calc)

        # We may not need this check because it could also be possible to only relax
        # an active center and not the whole molecule.
        # all_indices = set(it.chain(*self.layer_indices))
        # assert all_indices == set(range(len(self.geometry.atoms)))

    @property
    def layer_num(self):
        return len(self.layers)

    def optimize(self):
        # np.testing.assert_allclose(
        # self.layer_geoms[0].coords3d, self.geometry.coords3d[self.freeze_in_real]
        # )

        #####################################
        # Microiterations for real geometry #
        #####################################

        coords3d_org = self.geometry.coords3d.copy()
        coords3d_cur = coords3d_org.copy()
        # for indices, get_geom, get_opt in zip(
        # self.layer_indices, self.layer_get_geoms, self.layer_get_opts
        # ):
        for i, (indices, get_geom, get_opt) in enumerate(
            zip(self.layer_indices, self.layer_get_geoms, self.layer_get_opts)
        ):
            geom = get_geom(coords3d_cur)
            opt = get_opt(geom)
            if i == self.layer_num - 1:
                break
            opt.run()
            coords3d_cur[indices] = geom.coords3d
            # geom.jmol()

        ####################
        # Relax last layer #
        ####################

        # Calculate full ONIOM forces.
        results = self.geometry.get_energy_and_forces_at(coords3d_cur.flatten())
        forces = results["forces"]
        energy = results["energy"]
        # forces = self.geometry.forces
        # energy = self.geometry.energy
        self.energies.append(energy)
        self.forces.append(forces.copy())

        # 'geom' and 'indices' for the last layer were defined in the for-loop
        # above, before breaking out.
        last_layer_forces = forces.reshape(-1, 3)[indices].flatten()
        # geom.jmol()
        geom._forces = last_layer_forces
        geom._energy = energy
        geom_coords3d_org = geom.coords3d.copy()

        opt.coords.append(geom.coords.copy())
        opt.cart_coords.append(geom.cart_coords.copy())

        # Calculate one step
        int_step = opt.optimize()
        opt.steps.append(int_step)
        geom.coords = geom.coords + int_step
        coords3d_cur[indices] = geom.coords3d

        full_step = coords3d_cur - coords3d_org
        return full_step.flatten()


class LayerOptEven(Optimizer):
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
            self.layer_indices.append(indices)
            try:
                next_indices = full_expand(layers[i + 1]["indices"])
            except IndexError:
                next_indices = None

            # Calculator
            try:
                address = layer["address"]
                calc = IPIServer(address=address)
            except KeyError:
                pass

            # Geometry
            def get_geom_getter():
                layer_indices = indices
                layer_atoms = [atoms[i] for i in layer_indices]
                freeze_atoms = next_indices
                layer_calc = calc

                def get_geom(coords3d, coord_type="cart"):
                    layer_coords = coords3d[layer_indices].flatten()
                    geom = Geometry(
                        layer_atoms,
                        layer_coords,
                        freeze_atoms=freeze_atoms,
                        coord_type=coord_type,
                    )
                    geom.set_calculator(layer_calc)
                    return geom

                return get_geom

            geom_getter = get_geom_getter()
            if i == (len(layers) - 1):
                model_geom = geom_getter(self.geometry.coords3d, "redund")

                def geom_getter(coords3d):
                    return model_geom

            self.layer_get_geoms.append(geom_getter)

            # Optimizer
            def get_opt_getter():
                key = f"layer{i}"
                opt_kwargs = {
                    "prefix": key,
                    "max_cycles": 150,
                    "thresh": self.thresh,
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
                model_opt.prepare_opt()  # TODO: do this outside of constructor

                def get_opt(geom):
                    return model_opt

            self.layer_get_opts.append(get_opt)
        self.geometry.set_calculator(calc)

        # We may not need this check because it could also be possible to only relax
        # an active center and not the whole molecule.
        # all_indices = set(it.chain(*self.layer_indices))
        # assert all_indices == set(range(len(self.geometry.atoms)))

    @property
    def layer_num(self):
        return len(self.layers)

    def optimize(self):
        # np.testing.assert_allclose(
        # self.layer_geoms[0].coords3d, self.geometry.coords3d[self.freeze_in_real]
        # )

        #####################################
        # Microiterations for real geometry #
        #####################################

        coords3d_org = self.geometry.coords3d.copy()
        coords3d_cur = coords3d_org.copy()
        # for indices, get_geom, get_opt in zip(
        # self.layer_indices, self.layer_get_geoms, self.layer_get_opts
        # ):
        for i, (indices, get_geom, get_opt) in enumerate(
            zip(self.layer_indices, self.layer_get_geoms, self.layer_get_opts)
        ):
            geom = get_geom(coords3d_cur)
            opt = get_opt(geom)
            if i == self.layer_num - 1:
                break
            opt.run()
            coords3d_cur[indices] = geom.coords3d
            # geom.jmol()

        ####################
        # Relax last layer #
        ####################

        # Calculate full ONIOM forces.
        results = self.geometry.get_energy_and_forces_at(coords3d_cur.flatten())
        forces = results["forces"]
        energy = results["energy"]
        # forces = self.geometry.forces
        # energy = self.geometry.energy
        self.energies.append(energy)
        self.forces.append(forces.copy())

        # 'geom' and 'indices' for the last layer were defined in the for-loop
        # above, before breaking out.
        last_layer_forces = forces.reshape(-1, 3)[indices].flatten()
        # geom.jmol()
        geom._forces = last_layer_forces
        geom._energy = energy
        geom_coords3d_org = geom.coords3d.copy()

        opt.coords.append(geom.coords.copy())
        opt.cart_coords.append(geom.cart_coords.copy())

        # Calculate one step
        int_step = opt.optimize()
        opt.steps.append(int_step)
        geom.coords = geom.coords + int_step
        coords3d_cur[indices] = geom.coords3d

        full_step = coords3d_cur - coords3d_org
        return full_step.flatten()
