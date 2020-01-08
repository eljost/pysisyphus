#!/usr/bin/env python3

from collections import namedtuple

import numpy as np

from pysisyphus.calculators import Gaussian16, OpenMolcas, ORCA, Psi4, Turbomole, XTB
from pysisyphus.calculators.Calculator import Calculator
from pysisyphus.Geometry import Geometry
from pysisyphus.elem_data import COVALENT_RADII as CR
from pysisyphus.intcoords.findbonds import get_bond_sets


CALC_DICT = {
    "g16": Gaussian16,
    "openmolcas": OpenMolcas.OpenMolcas,
    "orca": ORCA.ORCA,
    "psi4": Psi4,
    "turbomole": Turbomole.Turbomole,
    "xtb": XTB.XTB,
    # "pyscf": PySCF,
    # "pypsi4": PyPsi4,
    # "pyxtb": PyXTB,
}


Link = namedtuple("Link", "ind parent_ind atom g")


def cap_fragment(atoms, coords, fragment, link_atom="H"):
    coords3d = coords.reshape(-1, 3)

    high_set = set(fragment)
    ind_set = set(range(len(atoms)))
    rest = ind_set - high_set
    
    # Determine bond(s) that connect fragment with the rest
    bonds = get_bond_sets(atoms, coords3d)
    bond_sets = [set(b) for b in bonds]
    
    # Find all bonds that involve one atom of model. These bonds
    # connect the model to the real geometry. We want to cap these
    # bonds.
    break_bonds = [b for b in bond_sets if len(b & high_set) == 1]
    
    # Put capping atoms at every bond to break.
    # The model fragment size corresponds to the length of the union of
    # the model set and the atoms in break_bonds.
    capped_frag = high_set.union(*break_bonds)
    capped_inds = list(sorted(capped_frag))

    # Index map between the new model geometry and the original indices
    # in the real geometry.
    atom_map = {model_ind: real_ind for real_ind, model_ind
                in zip(capped_inds, range(len(capped_inds)))}
    
    # g = 0.723886 # Gaussian16
    g = 0.709 # Paper g98-ONIOM-implementation
    c3d = coords3d.copy()
    new_atoms = list(atoms)
    # links = dict()
    links = list()
    for bb in break_bonds:
        to_cap = bb - high_set
        assert len(to_cap) == 1
        ind = list(bb - to_cap)[0]
        parent_ind = tuple(to_cap)[0]
        new_ind = np.sum(np.array(fragment) < parent_ind)
        link = Link(ind=ind, parent_ind=parent_ind, atom=link_atom, g=g)
        links.append(link)
    
    return atom_map, links


class Model():

    def __init__(self, name, calc_level, calc,
                 parent_name, parent_calc_level, parent_calc,
                 atom_inds, parent_atom_inds, all_atom_inds):

        self.name = name
        self.calc_level = calc_level
        self.calc = calc

        self.parent_name = parent_name
        self.parent_calc_level = parent_calc_level
        self.parent_calc = parent_calc

        self.atom_inds = atom_inds
        self.parent_atom_inds = parent_atom_inds
        self.all_atom_inds = atom_inds

        self.links = list()
        self.capped = False

    def create_links(self, atoms, coords, debug=False):
        self.capped = True

        if self.parent_name is not None:
            _, self.links = cap_fragment(atoms, coords, self.atom_inds)
        self.capped_atom_num = len(self.atom_inds) + len(self.links)

        if debug:
            catoms, ccords = self.capped_atoms_coords(atoms, coords)
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

    def get_energy(self, atoms, coords):
        print("calc energy", self.__str__())
        catoms, ccoords = self.capped_atoms_coords(atoms, coords)
        energy = self.calc.get_energy(catoms, ccoords)["energy"]
        try:
            parent_energy = self.parent_calc.get_energy(catoms, ccoords)["energy"]
        except AttributeError:
            parent_energy = 0.
        return energy - parent_energy

    def get_forces(self, atoms, coords):
        print("calc energy", self.__str__())
        catoms, ccoords = self.capped_atoms_coords(atoms, coords)
        energy = self.calc.get_energy(catoms, ccoords)["energy"]
        try:
            parent_energy = self.parent_calc.get_energy(catoms, ccoords)["energy"]
        except AttributeError:
            parent_energy = 0.
        return {"energy": energy - parent_energy}
    
    def __str__(self):
        return f"Model({self.name}, {len(self.atom_inds)} atoms, " \
               f"level={self.calc_level}, parent_level={self.parent_calc_level})"


class ONIOMext(Calculator):

    def __init__(self, calcs, models, layers, geom):
        super().__init__()

        real_key = "real"
        assert real_key not in models, \
            f'"{real_key}" must not be defined in "models"!'
        assert real_key in calcs, \
            f'"{real_key}" must be defined in "calcs"!'

        # Handle layers
        self.layer_num = len(layers)
        layers = [
            [layer, ] if isinstance(layer, str) else layer for layer in layers
        ]
        # Convert single level layers into lists of length 1.
        self.layers = {
            i: layer for i, layer in enumerate(layers)
        }

        models[real_key] = {
            "calc": real_key,
            "inds": list(range(len(geom.atoms))),
        }

        # Handle models
        model_keys, model_sizes = zip(
                *{key: len(model["inds"]) for key, model in models.items()}.items()
        )
        model_keys = list(models.keys())
        model_sizes = list()
        model_atom_sets = list()
        for key in model_keys:
            model = models[key]
            inds = model["inds"]
            model_atom_sets.append(set(inds))
            model_sizes.append(len(inds))
        # Determine hierarchy of models, from smallest to biggest model
        sort_args = np.argsort(model_sizes)
        all_atom_inds = range(len(geom.atoms))

        print(model_keys)
        print(model_sizes)
        print(model_atom_sets)
        print(sort_args)

        # The models are sorted from smaller to bigger. Now we check every
        # model but the last (biggest) to which bigger model it belongs.
        # The last (biggest) model is assumed to be embedded in the real system.
        model_parents = dict()
        for i, model_ind in enumerate(sort_args[:-1]):
            # Exclude the current model and the real model that contains all atoms
            rest_inds = sort_args[i+1:]
            atom_set = model_atom_sets[model_ind]
            is_subset = [atom_set.issubset(model_atom_sets[ind]) for ind in rest_inds]

            # Each model should belong to all models of higher layers
            parent_ind = rest_inds[is_subset.index(True)]
            
            model_key = model_keys[model_ind]
            parent_key = model_keys[parent_ind]
            model_parents[model_key] = parent_key
            print(i, model_key, parent_key)

        print(model_parents)
        cur_calc_num = 0
        def get_calc(calc_key):
            nonlocal cur_calc_num

            kwargs = calcs[calc_key].copy()
            type_ = kwargs.pop("type")
            calc = CALC_DICT[type_](**kwargs, calc_number=cur_calc_num)
            cur_calc_num += 1
            return calc

        # Create models and required calculators
        self.models = list()
        for model, parent in model_parents.items():
            print("\t", model, parent)
            model_calc_key = models[model]["calc"]
            try:
                parent_calc_key = models[parent]["calc"]
            except KeyError:
                parent_calc_key = "real"

            print(model_calc_key, parent_calc_key)
            model_calc = get_calc(model_calc_key)
            parent_calc = get_calc(parent_calc_key)

            model = Model(
                name=model,
                calc_level=model_calc_key,
                calc=model_calc,
                parent_name=parent,
                parent_calc_level=parent_calc_key,
                parent_calc=parent_calc,
                atom_inds=models[model]["inds"],
                parent_atom_inds=models[parent]["inds"],
                all_atom_inds=all_atom_inds,
            )
            self.models.append(model)

        # All real model
        real_calc = get_calc(real_key)
        self.models.append(
            Model(
                name=real_key,
                calc_level=real_key,
                calc=real_calc,
                parent_name=None, parent_calc_level=None, parent_calc=None,
                atom_inds=all_atom_inds,
                parent_atom_inds=None,
                all_atom_inds=all_atom_inds,
            )
        )

        for m in self.models:
            print(m)

        # Create link atoms
        [model.create_links(geom.atoms, geom.coords) for model in self.models]

    def get_energy(self, atoms, coords):
        energy = sum(
            [model.get_energy(atoms, coords) for model in self.models]
        )
        return {"energy": energy}
