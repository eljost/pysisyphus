#!/usr/bin/env python3

from collections import namedtuple
import itertools as it
from pprint import pprint

import numpy as np
from scipy.spatial.distance import pdist

# from pysisyphus.calculators import *
from pysisyphus.calculators import Gaussian16, OpenMolcas, ORCA, Psi4, Turbomole, XTB
from pysisyphus.calculators.Calculator import Calculator
from pysisyphus.Geometry import Geometry
from pysisyphus.elem_data import COVALENT_RADII as CR


CALC_DICT = {
    # "g09": Gaussian09.Gaussian09,
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


LinkMap = namedtuple("LinkMap", "r1_ind r3_ind atom g")


def get_cov_radii_sum_array(atoms, coords):
    coords3d = coords.reshape(-1, 3)
    atom_indices = list(it.combinations(range(len(coords3d)),2))
    atom_indices = np.array(atom_indices, dtype=int)
    cov_rad_sums = list()
    for i, j in atom_indices:
        atom1 = atoms[i].lower()
        atom2 = atoms[j].lower()
        cov_rad_sums.append(CR[atom1] + CR[atom2])
    cov_rad_sums = np.array(cov_rad_sums)
    return cov_rad_sums


def get_bond_sets(atoms, coords3d, bond_factor=1.3):
    cdm = pdist(coords3d)
    # Generate indices corresponding to the atom pairs in the
    # condensed distance matrix cdm.
    atom_indices = list(it.combinations(range(len(coords3d)),2))
    atom_indices = np.array(atom_indices, dtype=int)
    cov_rad_sums = get_cov_radii_sum_array(atoms, coords3d.flatten())
    cov_rad_sums *= bond_factor
    bond_flags = cdm <= cov_rad_sums
    bond_indices = atom_indices[bond_flags]
    return bond_indices


def cap(atoms, coords, high_frag, link_atom="H"):
    coords3d = coords.reshape(-1, 3)

    high_set = set(high_frag)
    ind_set = set(range(len(atoms)))
    rest = ind_set - high_set
    
    # Determine bond(s) that connect high_frag with the rest
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
    link_maps = dict()
    for bb in break_bonds:
        to_cap = bb - high_set
        assert len(to_cap) == 1
        r1_ind = list(bb - to_cap)[0]
        r3_ind = tuple(to_cap)[0]
        r1 = c3d[r1_ind]
        r3 = c3d[r3_ind]
        r2 = r1 + g*(r3-r1)
        c3d[r3_ind] = r2
        new_atoms[r3_ind] = link_atom
        new_ind = np.sum(np.array(high_frag) < r3_ind)
        link_map = LinkMap(r1_ind=r1_ind, r3_ind=r3_ind, atom=link_atom, g=g)
        link_maps[new_ind] = link_map
    
    capped_atoms = [new_atoms[i] for i in capped_inds]
    capped_coords = c3d[capped_inds].flatten()
    capped_geom = Geometry(capped_atoms, capped_coords)
    
    return capped_geom, atom_map, link_maps


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

        self.links = None

    def create_links(self, atoms, coords):
        if self.parent_name is None:
            self.links = list()
            return

        pinds = self.parent_atom_inds
        parent_atoms = [atoms[i] for i in pinds]
        parent_coords = coords.reshape(-1, 3)[pinds]
        # capped_geom, atom_map, link_maps
        _, __, links = cap(parent_atoms, parent_coords, self.atom_inds)
        import pdb; pdb.set_trace()
        pass

    def capped_atoms_coords(self, atoms, coords):
        assert self.links is not None, "Did you forget to call create_links()?"

    def energy(atoms, coords):
        catoms, ccoords = self.capped_atoms_coords(atoms, coords)
        energy = self.calc.get_energy(catoms, ccoords)
        try:
            parent_energy = self.parent_calc.get_energy(catoms, ccoords)
        except AttributeError:
            parent_energy = 0.
        return energy - parent_energy
    
    def __str__(self):
        return f"Model({self.name}, {len(self.atom_inds)} atoms, " \
               f"level={self.calc_level}, parent_level={self.parent_calc_level}"


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
        def get_calc(calc_key):
            kwargs = calcs[calc_key].copy()
            type_ = kwargs.pop("type")
            return CALC_DICT[type_](**kwargs)

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

        pprint([str(m) for m in self.models])

        # Create link atoms
        [model.create_links(geom.atoms, geom.coords) for model in self.models]

    def get_energy(self, atoms, coords):
        return sum([model.get_energy(atoms, coords) for model in self.models])
