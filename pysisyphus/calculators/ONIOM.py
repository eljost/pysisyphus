#!/usr/bin/env python3

import itertools as it
from collections import namedtuple

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
    "orca": ORCA,
    "psi4": Psi4,
    "turbomole": Turbomole,
    "xtb": XTB,
    # "pyscf": PySCF,
    # "pypsi4": PyPsi4,
    # "pyxtb": PyXTB,
}


LinkMap = namedtuple("LinkMap", "r1_ind r3_ind g")


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


def cap(geom, high_frag):
    high_set = set(high_frag)
    ind_set = set(range(len(geom.atoms)))
    rest = ind_set - high_set
    
    # Determine bond(s) that connect high_frag with the rest
    # bonds, _, _ = geom.internal.prim_indices
    bonds = get_bond_sets(geom.atoms, geom.coords3d)
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
    c3d = geom.coords3d.copy()
    new_atoms = list(geom.atoms)
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
        new_atoms[r3_ind] = "H"
        new_ind = np.sum(np.array(high_frag) < r3_ind)
        link_map = LinkMap(r1_ind=r1_ind, r3_ind=r3_ind, g=g)
        link_maps[new_ind] = link_map
    
    capped_atoms = [new_atoms[i] for i in capped_inds]
    capped_coords = c3d[capped_inds].flatten()
    capped_geom = Geometry(capped_atoms, capped_coords)
    
    return capped_geom, atom_map, link_maps


class ONIOM(Calculator):

    def __init__(self, calc_dict, model_inds):
        super().__init__()

        self.calc_dict = calc_dict
        self.model_inds = model_inds

        high_level = self.calc_dict["high"]
        high_type = high_level.pop("type")
        high_cls = CALC_DICT[high_type]

        low_level = self.calc_dict["low"]
        low_type = low_level.pop("type")
        low_cls = CALC_DICT[low_type]

        self.model_high_calc = high_cls(calc_number=0, **high_level)
        self.real_low_calc = low_cls(calc_number=1, **low_level)
        self.model_low_calc = low_cls(calc_number=2, **low_level)

    def get_forces(self, atoms, coords):
        tmp_geom = Geometry(atoms, coords)
        results = self.oniom2(tmp_geom)
        return results
    
    def oniom2(self, real_geom):
        model_geom, atom_map, links = cap(real_geom, self.model_inds)
            
        results_3 = self.real_low_calc.get_forces(real_geom.atoms, real_geom.cart_coords)
        results_1 = self.model_low_calc.get_forces(model_geom.atoms, model_geom.cart_coords)
        results_2 = self.model_high_calc.get_forces(model_geom.atoms, model_geom.cart_coords)

        E3 = results_3["energy"]
        E1 = results_1["energy"]
        E2 = results_2["energy"]
        g3 = -results_3["forces"].reshape(-1, 3)
        g1 = -results_1["forces"].reshape(-1, 3)
        g2 = -results_2["forces"].reshape(-1, 3)
        
        # ONIOM energy
        energy_oniom = E3 - E1 + E2
        
        # ONIOM gradient
        gradient_oniom = g3 
        for model_ind, real_ind in atom_map.items():
            qm_correction = -g1[model_ind] + g2[model_ind]
            if model_ind in links:
                model_sum = -g1[model_ind] + g3[model_ind]
                r1_ind, r3_ind, g = links[model_ind]
                gradient_oniom[r1_ind] += (1-g) * qm_correction
                gradient_oniom[r3_ind] += g * qm_correction
            else:
                gradient_oniom[real_ind] += qm_correction

        results = {
            "energy": energy_oniom,
            "forces": -gradient_oniom.flatten(),
        }
        return results
