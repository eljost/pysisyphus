# [1] https://www.sciencedirect.com/science/article/pii/S0166128098004758
#     https://doi.org/10.1016/S0166-1280(98)00475-8
#     Dapprich, Frisch, 1998
# [2] https://onlinelibrary.wiley.com/doi/abs/10.1002/9783527629213.ch2
#     Clemente, Frisch, 2010
#
# Not implemented in pysisyphus
#
# [2] https://aip.scitation.org/doi/pdf/10.1063/1.2814164?class=pdf
#     QM/QM ONIOM EE based on Mulliken charges
#     Hratchian, Raghavachari, 2008
# [3] https://aip.scitation.org/doi/full/10.1063/1.3315417<Paste>
#     QM/QM ONIOM EE based on Löwdin charges
#     Mayhall, Hratchian, 2010
# [4] https://www.frontiersin.org/articles/10.3389/fchem.2018.00089/full
#     Overview on hybrid methods
#
# Excited state ONIOM
# [5] https://aip.scitation.org/doi/pdf/10.1063/1.4972000?class=pdf
# [6] https://pubs.rsc.org/en/content/articlehtml/2012/pc/c2pc90007f


import itertools as it
import logging
from collections import namedtuple

import numpy as np

from pysisyphus.calculators import Gaussian16, OpenMolcas, ORCA, Psi4, Turbomole, XTB
from pysisyphus.calculators.Calculator import Calculator
from pysisyphus.Geometry import Geometry
from pysisyphus.elem_data import COVALENT_RADII as CR
from pysisyphus.intcoords.setup import get_bond_sets


CALC_DICT = {
    "g16": Gaussian16,
    "openmolcas": OpenMolcas.OpenMolcas,
    "orca": ORCA,
    "psi4": Psi4,
    "turbomole": Turbomole,
    "xtb": XTB.XTB,
    # "pypsi4": PyPsi4,
    # "pyxtb": PyXTB,
}
try:
    from pysisyphus.calculators.PySCF import PySCF
    CALC_DICT["pyscf"] = PySCF
except ImportError:
    print("Error importing PySCF in ONIOMv2")


Link = namedtuple("Link", "ind parent_ind atom g")


def get_g_value(atom, parent_atom, link_atom):
    cr, pcr, lcr = [CR[a.lower()] for a in (atom, parent_atom, link_atom)]

    # Ratio between sum of CR_atom and CR_link with sum of CR_atom CR_parent_atom.
    # See [1] Sect. 2.2 page 5.
    g = (cr + lcr) / (cr + pcr)
    return g


def cap_fragment(atoms, coords, fragment, link_atom="H", g=0.709):
    coords3d = coords.reshape(-1, 3)

    fragment_set = set(fragment)
    
    # Determine bond(s) that connect fragment with the rest
    bonds = get_bond_sets(atoms, coords3d)
    bond_sets = [set(b) for b in bonds]
    
    # Find all bonds that involve one atom of model. These bonds
    # connect the model to the real geometry. We want to cap these
    # bonds.
    break_bonds = [b for b in bond_sets if len(b & fragment_set) == 1]
    
    # Put capping atoms at every bond to break.
    # The model fragment size corresponds to the length of the union of
    # the model set and the atoms in break_bonds.
    capped_frag = fragment_set.union(*break_bonds)
    capped_inds = list(sorted(capped_frag))

    # Index map between the new model geometry and the original indices
    # in the real geometry.
    atom_map = {model_ind: real_ind for real_ind, model_ind
                in zip(capped_inds, range(len(capped_inds)))}
    
    links = list()
    for bb in break_bonds:
        to_cap = bb - fragment_set
        assert len(to_cap) == 1
        ind = list(bb - to_cap)[0]
        parent_ind = tuple(to_cap)[0]
        if g is None:
            g = get_g_value(atoms[ind], atoms[parent_ind], link_atom)
        link = Link(ind=ind, parent_ind=parent_ind, atom=link_atom, g=g)
        links.append(link)
    
    return atom_map, links


def atom_inds_to_cart_inds(atom_inds):
    stencil = np.array((0, 1, 2), dtype=int)
    size_ = len(atom_inds)
    cart_inds = np.tile(stencil, size_) + np.repeat(atom_inds, 3)*3
    return cart_inds


class Model():

    def __init__(self, name, calc_level, calc,
                 parent_name, parent_calc_level, parent_calc,
                 atom_inds, parent_atom_inds, use_link_atoms=True):

        self.name = name
        self.calc_level = calc_level
        self.calc = calc

        self.parent_name = parent_name
        self.parent_calc_level = parent_calc_level
        self.parent_calc = parent_calc

        self.atom_inds = atom_inds
        self.parent_atom_inds = parent_atom_inds

        self.use_link_atoms = use_link_atoms

        self.links = list()
        self.capped = False
        self.J = None

    def log(self, message=""):
        logger = logging.getLogger("calculator")
        logger.debug(self.__str__() + " " + message)

    def create_links(self, atoms, coords, debug=False):
        self.capped = True

        if self.use_link_atoms and self.parent_name is not None:
            _, self.links = cap_fragment(atoms, coords, self.atom_inds)
        self.capped_atom_num = len(self.atom_inds) + len(self.links)
        for i, link in enumerate(self.links):
            ind, parent_ind = link.ind, link.parent_ind
            self.log(f"\tCreated Link atom ({link.atom}) between {atoms[ind]}{ind} "
                     f"and {atoms[parent_ind]}{parent_ind} (g={link.g:.6f})")

        if len(self.links) == 0:
            self.log("Didn't create any link atoms!\n")

        try:
            self.J = self.get_jacobian()
        except TypeError:
            self.log("Skipping definition of jacobian shape")

        if debug:
            catoms, ccoords = self.capped_atoms_coords(atoms, coords)
            geom = Geometry(catoms, ccoords)
            geom.jmol()

    def capped_atoms_coords(self, all_atoms, all_coords):
        assert self.capped, "Did you forget to call create_links()?"

        org_atom_num = len(self.atom_inds)
        c3d = all_coords.reshape(-1, 3)

        capped_atoms = [all_atoms[i] for i in self.atom_inds]
        # Initialize empty coordinate array
        capped_coords = np.zeros((self.capped_atom_num, 3))
        # Copy non-capped coordinates
        capped_coords[:org_atom_num] = c3d[self.atom_inds]

        for i, link in enumerate(self.links):
            capped_atoms.append(link.atom)

            r1 = c3d[link.ind]
            r3 = c3d[link.parent_ind]
            r2 = r1 + link.g*(r3-r1)
            capped_coords[org_atom_num + i] = r2
        return capped_atoms, capped_coords

    # def coords(self, atoms, coords):
        # """Wrapper for self.capped_atoms_coords."""
        # return self.capped_atoms_coords(atoms, coords)

    def get_jacobian(self):
        try:
            # Shape of Jacobian is (model + link, real)
            jac_shape = (len(self.atom_inds)*3 + len(self.links)*3,
                         len(self.parent_atom_inds)*3)
        except TypeError:
            # Return array instead of 1, because the integer 1 can't be transposed
            # as we have to do for the hessian calculation.
            return np.array(1)

        # assert self.parent_atom_inds == list(range(self.parent_atom_inds[0],
                                                   # self.parent_atom_inds[-1]+1))

        J = np.zeros(jac_shape)
        # Stencil for diagonal elements of 3x3 submatrix
        stencil = np.array((0, 1, 2), dtype=int)
        size_ = len(self.atom_inds)

        model_rows = np.arange(size_ * 3)
        # When more than two layers are present the inner layers aren't directly
        # embedded in the outermost layer. This means parent_inds does not begin
        # at 0, but with a higher index. So we need a map of the actual indices
        # (not starting at 0) to the indices in the Jacobian which start at 0.
        atom_inds = [self.parent_atom_inds.index(ind) for ind in self.atom_inds]
        ind_map = {k: v for k, v in zip(self.atom_inds, atom_inds)}
        model_cols = atom_inds_to_cart_inds(atom_inds)
        J[model_rows, model_cols] = 1

        # Link atoms
        link_start = model_rows.max() + 1
        for i, (ind, parent_ind, atom, g) in enumerate(self.links):
            rows = link_start + i*3 + stencil
            cols = ind_map[ind]*3 + stencil

            J[rows, cols] = 1 - g
            try:
                parent_cols = self.parent_atom_inds.index(parent_ind)*3 + stencil
                J[rows, parent_cols] = g
            # Raised when link atom is not coupled to layer above, but
            # to a layer higher above.
            except ValueError:
                pass

        return J

    def get_energy(self, atoms, coords, point_charges=None,
                   parent_correction=True):
        self.log("Energy calculation")
        catoms, ccoords = self.capped_atoms_coords(atoms, coords)

        prepare_kwargs = {
            "point_charges": point_charges,
        }

        self.log("Calculation at layer level")
        results = self.calc.get_energy(catoms, ccoords, prepare_kwargs)
        energy = results["energy"]

        # Calculate correction if parent layer is present and it is requested
        if (self.parent_calc is not None) and parent_correction:
            self.log("Calculation at parent layer level")
            parent_results= self.parent_calc.get_energy(catoms, ccoords, prepare_kwargs)
            parent_energy = parent_results["energy"]
            energy -= parent_energy
        elif not parent_correction:
            self.log("No parent correction!")

        return energy

    def get_forces(self, atoms, coords, point_charges=None,
                   parent_correction=True):
        self.log("Force calculation")
        catoms, ccoords = self.capped_atoms_coords(atoms, coords)

        prepare_kwargs = {
            "point_charges": point_charges,
        }

        self.log("Calculation at layer level")
        results = self.calc.get_forces(catoms, ccoords, prepare_kwargs)
        forces = results["forces"]
        energy = results["energy"]
        forces = forces.dot(self.J)

        # Calculate correction if parent layer is present and it is requested
        if (self.parent_calc is not None) and parent_correction:
            self.log("Calculation at parent layer level")
            parent_results = self.parent_calc.get_forces(catoms, ccoords, prepare_kwargs)
            parent_forces = parent_results["forces"]
            parent_energy = parent_results["energy"]

            # Correct energy and forces
            energy -= parent_energy
            forces -= parent_forces.dot(self.J)
        elif not parent_correction:
            self.log("No parent correction!")

        return energy, forces

    def get_hessian(self, atoms, coords, point_charges=None,
                    parent_correction=True):
        self.log("Hessian calculation")
        catoms, ccoords = self.capped_atoms_coords(atoms, coords)

        # prepare_kwargs = {
            # "point_charges": point_charges,
        # }
        if point_charges is not None:
            raise Exception("point_charges & hessian is not yet implemented")

        self.log("Calculation at layer level")
        # results = self.calc.get_hessian(catoms, ccoords, prepare_kwargs)
        results = self.calc.get_hessian(catoms, ccoords)
        hessian = results["hessian"]
        energy = results["energy"]
        hessian = self.J.T.dot(hessian.dot(self.J))

        # Calculate correction if parent layer is present and it is requested
        if (self.parent_calc is not None) and parent_correction:
            self.log("Calculation at parent layer level")
            parent_results = self.parent_calc.get_hessian(catoms, ccoords)
            parent_hessian = parent_results["hessian"]
            parent_energy = parent_results["energy"]

            # Correct energy and hessian
            energy -= parent_energy
            hessian -= self.J.T.dot(parent_hessian.dot(self.J))
        elif not parent_correction:
            self.log("No parent correction!")

        return energy, hessian

    def parse_charges(self):
        charges = self.calc.parse_charges()
        try:
            parent_charges = self.parent_calc.parse_charges()
        except AttributeError:
            parent_charges = None

        return charges, parent_charges
    
    def __str__(self):
        return f"Model({self.name}, {len(self.atom_inds)} atoms, " \
               f"level={self.calc_level}, parent_level={self.parent_calc_level})"

    def __repr__(self):
        return self.__str__()


class ONIOM(Calculator):

    def __init__(self, calcs, models, geom, layers=None, embedding=None,
                 real_key="real", use_link_atoms=True, *args, **kwargs):
        """
        layer: list of models
            len(layer) == 1: normal ONIOM, len(layer) >= 1: multicenter ONIOM.
        model:
            (sub)set of all atoms that resides in a certain layer and has
            a certain calculator.
        """

        super().__init__(*args, **kwargs)

        valid_embeddings = (None, "electronic")
        assert embedding in valid_embeddings, f"Valid embeddings are: {valid_embeddings}"
        self.embedding = embedding

        assert real_key not in models, \
            f'"{real_key}" must not be defined in "models"!'
        assert real_key in calcs, \
            f'"{real_key}" must be defined in "calcs"!'

        self.use_link_atoms = use_link_atoms

        # When no ordering of layers is given we try to guess it from
        # the size of the respective models. It's probably a better idea
        # to always specify the layer ordering though ;)
        if layers is None:
            self.log("No explicit layer ordering specified! Determining layer "
                     "hierarchy from model sizes. This does not support multi-"
                     "center ONIOM!")
            as_list = [(key, val) for key, val in models.items()]
            # Determine hierarchy of models, from biggest to smallest model
            layers = [
                key for key, val
                in sorted(as_list, key=lambda model: -len(model[1]["inds"]))
            ]

        assert real_key not in layers, \
            f'"{real_key}" must not be defined in "layers"!'

        ############
        #          #
        #  LAYERS  #
        #          #
        ############

        # Add real model and layer as they are missing right now. The real
        # layer is always the last layer. The real layer is always calculated
        # by the 'realkey'-calculator.
        layers = [real_key] + layers
        models[real_key] = {
            "calc": real_key,
            "inds": list(range(len(geom.atoms))),
        }
        self.log(f"Layer-ordering from big to small: {layers}")

        # Single-model layers will be given as strings. As we also support
        # multicenter-ONIOM there may also be layers that are given as lists
        # that contain multiple models per layer.
        # Now we convert the single-model layers to lists of length 1, so
        # every layer is a list.
        layers = [
            [layer, ] if isinstance(layer, str) else layer for layer in layers
        ]
        self.layer_num = len(layers)
        assert self.layer_num > 1, "ONIOM with only 1 layer requested. Aborting!"

        ############
        #          #
        #  MODELS  #
        #          #
        ############

        # Create mapping between model and its parent layer. Actually
        # this is a bit hacky right now, as the mapping should not be between
        # model and parent layer, but between model and parent model.
        # This way we expect the parent layer to have the same calculator
        # throughout, so multicenter ONIOM with different calculators
        # in all but the smallest layer (highest level) is not well defined.
        #
        # If multicenter ONIOM in an intermediate layer is useful may
        # be another question to be answered ;).
        self.model_parent_layers = dict()
        for i, layer in enumerate(layers[1:]):
            self.model_parent_layers.update(
                {model: i for model in layer}
            )
        model_keys = list(it.chain(*layers))

        cur_calc_num = 0
        def get_calc(calc_key, base_name=None):
            """Helper function for easier generation of calculators
            with incrementing calc_number."""
            nonlocal cur_calc_num

            kwargs = calcs[calc_key].copy()
            type_ = kwargs.pop("type")

            kwargs["calc_number"] = cur_calc_num
            if base_name is not None:
                kwargs["base_name"] = base_name

            calc = CALC_DICT[type_](**kwargs)
            cur_calc_num += 1
            return calc

        # Create models and required calculators
        self.models = list()
        self.layers = [list() for _ in layers]
        for model in model_keys[1:]:
            parent_layer_ind = self.model_parent_layers[model]
            parent_layer = layers[parent_layer_ind]
            parent_calc_keys = set([models[model]["calc"] for model in parent_layer])
            assert len(parent_calc_keys) == 1, \
                "It seems you are trying to run a multicenter ONIOM setup in " \
                "an intermediate layer with different calculators. This is " \
                "not supported right now."

            parent = parent_layer[0]
            model_calc_key = models[model]["calc"]
            parent_calc_key = models[parent]["calc"]

            model_base_name = f"{model}_{model_calc_key}"
            model_calc = get_calc(model_calc_key, base_name=model_base_name)
            parent_base_name = f"{model}_parent"
            parent_calc = get_calc(parent_calc_key, base_name=parent_base_name)

            model = Model(
                name=model,
                calc_level=model_calc_key,
                calc=model_calc,
                parent_name=parent,
                parent_calc_level=parent_calc_key,
                parent_calc=parent_calc,
                atom_inds=models[model]["inds"],
                parent_atom_inds=models[parent]["inds"],
                use_link_atoms=self.use_link_atoms,
            )
            self.models.append(model)
            self.layers[parent_layer_ind+1].append(model)

        # All real model
        real_calc = get_calc(real_key)
        real_model = Model(
                name=real_key,
                calc_level=real_key,
                calc=real_calc,
                parent_name=None, parent_calc_level=None, parent_calc=None,
                atom_inds=list(range(len(geom.atoms))),
                parent_atom_inds=None,
            )
        self.models.insert(0, real_model)
        self.layers[0].append(real_model)

        # Reverse order of models so the first model is the real system
        # self.models = self.models[::-1]

        self.log("Created all ONIOM layers:")
        for model in self.models:
            self.log("\t"+str(model))

        # Create link atoms
        [model.create_links(geom.atoms, geom.cart_coords) for model in self.models]

        # And do a quick sanity check
        assert len(self.models[0].links) == 0, \
            "There must not be any links in the 'real' layer!"
        # Look for link atoms that appear in two adjacent layers. In such situations
        # the higher layer is coupled to a layer two levels below. This may be a bad
        # idea.
        for i, (lower_model, model) in enumerate(zip(self.models[:-1], self.models[1:])):
            lower_links = lower_model.links
            links = model.links
            same_links = [link for link in links if link in lower_links]
            if same_links:
                print(f"Found {len(same_links)} link(s) that appear(s) in two layers!")
                for j, link in enumerate(same_links):
                    print(f"\t{j:02d}: {link}")
                print(f"Your current setup couples layer '{model.name}' to "
                      f"to layer '{self.models[i-1].name}' two levels below! "
                       "This is probably a bad idea!")

        self.log(f"Created ONIOM calculator with {self.layer_num} layers and "
                 f"{len(self.models)} models.")

    def run_calculations(self, atoms, coords, method):
        titles = {
            None: "",
            "electronic": "Electronic embedding",
        }
        self.log(f"{titles[self.embedding]} ONIOM calculation")

        all_results = list()
        point_charges = None
        for i, layer in enumerate(self.layers):
            # Only consider charges that belong to atoms in the parent
            # layer. Otherwise this would result in additonal charges at
            # the same positions as the atoms we would like to calculate.
            if self.embedding == "electronic" and (i > 0):
                parent_layer = self.layers[i-1]
                assert len(parent_layer) == 1, \
                    "Multicenter ONIOM in intermediate layer is not supported!"
                parent_model = parent_layer[0]
                parent_inds = set(parent_model.atom_inds)
                charges, _ = parent_model.parse_charges()

                layer_inds = set(*it.chain([model.atom_inds for model in layer]))
                # Determine indices of atoms that are in the parent layer, but
                # not in the current layer
                only_parent_inds = list(parent_inds - layer_inds)
                ee_charges = charges[only_parent_inds]
                ee_charge_sum = sum(ee_charges)

                point_charges = np.zeros((ee_charges.size, 4))
                point_charges[:,:3] = coords.reshape(-1, 3)[only_parent_inds]
                point_charges[:,3] = ee_charges

                self.log(f"Polarizing calculation in layer {i} ({layer}) by "
                         f"charges from layer {i-1} ({self.layers[i-1]})."
                )
                self.log(f"sum(charges)={ee_charge_sum:.4f}")

            results = [getattr(model, method)(atoms, coords, point_charges=point_charges)
                       for model in layer]
            all_results.extend(results)
        self.calc_counter += 1

        return all_results

    def run_calculation(self, atoms, coords):
        self.log("run_calculation() called. Doing simple energy calculation!")
        return self.get_energy(atoms, coords)

    def get_energy(self, atoms, coords):
        all_energies = self.run_calculations(atoms, coords, "get_energy")

        energy = sum(all_energies)

        return {
                "energy": energy,
        }

    def get_forces(self, atoms, coords):
        all_results = self.run_calculations(atoms, coords, "get_forces")

        energies, forces_ = zip(*all_results)
        forces_ = [np.array(f).reshape(-1, 3) for f in forces_]
        energy = sum(energies)

        forces = forces_[0]
        for mdl, f in zip(self.models[1:], forces_[1:]):
            forces[mdl.parent_atom_inds] += f

        return {
                "energy": energy,
                "forces": forces.flatten(),
        }

    def get_hessian(self, atoms, coords):
        all_results = self.run_calculations(atoms, coords, "get_hessian")

        energies, hessians = zip(*all_results)
        energy = sum(energies)

        hessian = hessians[0]
        for mdl, h in zip(self.models[1:], hessians[1:]):
            inds = atom_inds_to_cart_inds(mdl.parent_atom_inds)
            # Keep in mind that we modify hessians[0] in place
            hessian[inds[:,None], inds[None,:]] += h

        return {
                "energy": energy,
                "hessian": hessian,
        }

    def atom_inds_in_layer(self, index, exclude_inner=False):
        """Returns list of atom indices in layer at index.

        Atoms that also appear in inner layer can be excluded on request.

        Parameters
        ----------
        index : int
            pasd
        exclude_inner : bool, default=False, optional
            Whether to exclude atom indices that also appear in inner layers.

        Returns
        -------
        atom_indices : list
            List containing the atom indices in the selected layer.
        """

        layer = self.layers[index]
        atom_inds = list(it.chain(*[model.atom_inds for model in layer]))
        if exclude_inner and (index < len(self.layers)-1):
            lower_inds = self.atom_inds_in_layer(index+1)
            # Drop indices that appear in inner layers
            atom_inds = [i for i in atom_inds if i not in lower_inds]
        return atom_inds

    def calc_layer(self, atoms, coords, index, parent_correction=True):
        layer = self.layers[index]
        assert len(layer) == 1, \
            "Multicenter not yet supported!"
        model, = layer
        result = model.get_forces(atoms, coords, parent_correction=parent_correction)
        return result
