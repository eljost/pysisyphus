from dataclasses import dataclass
import itertools as it
from typing import List

import numpy as np
from numpy.typing import NDArray

from pysisyphus.Geometry import Geometry
from pysisyphus.io.pdb import parse_pdb
from pysisyphus.io.psf import parse_psf


@dataclass
class Atom:
    id: int
    segment: str
    resid: int
    resname: str
    name: str
    type: str
    charge: float
    mass: float
    coords: NDArray[float]
    element: str

    @staticmethod
    def from_psf_line(line, coords, element):
        return Atom(
            id=line["id"],
            segment=line["segment"],
            resid=line["resid"],
            resname=line["resname"],
            name=line["atom_name"],
            type=line["atom_type"],
            charge=line["charge"],
            mass=line["mass"],
            coords=coords,
            element=element,
        )


@dataclass
class Residue:
    id: int
    name: str
    segment: str
    atoms: List[Atom]

    @staticmethod
    def from_psf_lines(lines, coords3d, elements):
        assert len(lines) == len(coords3d) == len(elements)
        atoms = [
            Atom.from_psf_line(line, coords, element)
            for line, coords, element in zip(lines, coords3d, elements)
        ]
        atom0 = atoms[0]
        resid0 = atom0.resid
        resname0 = atom0.resname
        segment0 = atom0.segment
        if len(lines) > 1:
            resids, resnames = zip(*[(atom.resid, atom.resname) for atom in atoms[1:]])
            assert all([resname == resname0 for resname in resnames])
            assert all([resid == resid0 for resid in resids])
        return Residue(resid0, resname0, segment0, atoms)

    @property
    def key(self):
        return self.segment, self.id

    @property
    def charge(self):
        charge = sum([atom.charge for atom in self.atoms])
        # assert abs(charge % 1) <= 1e-10  # I guess charges must not be integer ...
        return charge

    @property
    def masses(self):
        return np.array([atom.mass for atom in self.atoms])

    @property
    def total_mass(self):
        return sum(self.masses)

    @property
    def atom_indices(self):
        return [atom.id for atom in self.atoms]

    @property
    def elements_coords3d(self):
        elements = [atom.element for atom in self.atoms]
        return elements, self.coords3d

    @property
    def coords3d(self):
        return np.array([atom.coords for atom in self.atoms])

    @property
    def center_of_mass(self):
        return (
            1 / self.total_mass * np.sum(self.coords3d * self.masses[:, None], axis=0)
        )

    @property
    def com(self):
        return self.center_of_mass

    def __len__(self):
        return len(self.atoms)

    def __str__(self):
        return f"{self.resname}{self.resid}"


def residues_from_psf(psf_data, atoms, coords3d, atom_map):
    psf_atoms = psf_data["atoms"]
    res_key = lambda atom: (atom["segment"], atom["resid"])
    psf_atoms = sorted(psf_atoms, key=res_key)
    psf_atoms_by_res = it.groupby(psf_atoms, key=res_key)
    residues = dict()
    for _, g in psf_atoms_by_res:
        res_atoms = sorted(g, key=lambda atom: atom["id"])
        atom_inds = [atom_map[atom["id"]] for atom in res_atoms]
        res_coords3d = coords3d[atom_inds]
        res_elements = [atoms[i] for i in atom_inds]
        res = Residue.from_psf_lines(res_atoms, res_coords3d, res_elements)
        residues[res.key] = res
    return residues


def residues_within_com_dist(
    atoms, coords3d, atom_map, psf_data, within_resid, within_dist
):

    residuesd = residues_from_psf(psf_data, atoms, coords3d, atom_map)

    coms = {key: res.com for key, res in residuesd.items()}

    def within(resid, com_dist):
        ref_com = coms[resid]
        res_ids_below_dist = [
            key
            for key, res in residuesd.items()
            if np.linalg.norm(coms[res.key] - ref_com) <= com_dist
        ]
        return res_ids_below_dist

    within_cys = within(within_resid, within_dist)
    wresidues = [residuesd[resid] for resid in within_cys]
    return wresidues


def geom_from_residues(residues):
    atoms = list()
    coords3d = np.zeros((sum([len(res) for res in residues]), 3))
    i = 0
    for res in residues:
        res_elements, res_coords3d = res.elements_coords3d
        len_res = len(res)
        coords3d[i : i + len_res] = res_coords3d
        atoms.extend(res_elements)
        i += len_res
    return Geometry(atoms, coords3d)


def link_atoms_for_residues(
    residues, bonds, coords3d, atom_map, link_element="H", g=0.709
):
    atom_inds = list(it.chain(*[res.atom_indices for res in residues]))

    bond_dict = dict()
    for bond in bonds:
        from_, to_ = bond
        bond_dict.setdefault(from_, set()).add(to_)
        bond_dict.setdefault(to_, set()).add(from_)
    atom_set = set(atom_inds)

    cut_bonds = list()
    # Check all bonds from residue-atoms. Determine which bonds are cut.
    for atom in atom_inds:
        try:
            cut_bonds_with = bond_dict[atom] - atom_set
        # Single ions/atoms may not have any bonds.
        except KeyError:
            cut_bonds_with = list()
        for cbw in cut_bonds_with:
            cut_bonds.append((atom, cbw))

    link_atoms = list()
    link_hosts = np.zeros(len(cut_bonds))
    link_coords3d = np.zeros((len(cut_bonds), 3))
    for i, cut_bond in enumerate(cut_bonds):
        from_, to_ = cut_bond
        link_atoms.append(link_element)
        link_hosts[i] = to_
        from_coords = coords3d[atom_map[from_]]
        to_coords = coords3d[atom_map[to_]]
        bond_vec = to_coords - from_coords
        link_coords3d[i] = from_coords + g * bond_vec
    link_atoms = tuple(link_atoms)
    return link_hosts, link_atoms, link_coords3d


def cluster_from_psf_pdb(psf_fn, pdb_fn, within_resid, within_dist):
    atoms, coords, _, atom_map = parse_pdb(pdb_fn)
    coords3d = coords.reshape(-1, 3)
    print("Loaded PDB data.")

    # Topology etc. from PSF
    psf_data = parse_psf(psf_fn)
    # import pickle
    # with open("d.pickle", "wb") as handle:
    # pickle.dump(psf_data, handle)
    # with open("d.pickle", "rb") as handle:
    # psf_data = pickle.load(handle)
    print("Loaded PSF data.")

    residues = residues_within_com_dist(
        atoms, coords3d, atom_map, psf_data, within_resid, within_dist
    )
    geom = geom_from_residues(residues)
    print("Create cluster geometry.")
    bonds = np.array(psf_data["nbond"]["inds"], dtype=int).reshape(-1, 2)
    link_hosts, link_atoms, link_coords3d = link_atoms_for_residues(
        residues, bonds, coords3d, atom_map
    )
    sat_geom = Geometry(
        geom.atoms + link_atoms, np.concatenate((geom.coords3d, link_coords3d), axis=0)
    )
    print("Created satured geometry with link atoms.")
    return sat_geom, sum([res.charge for res in residues])
