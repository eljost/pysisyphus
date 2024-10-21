from dataclasses import dataclass
import itertools as it
import pickle
from typing import List

import numpy as np

from pysisyphus.Geometry import Geometry
from pysisyphus.io.pdb import parse_pdb
from pysisyphus.io.psf import parse_psf
from pysisyphus.TablePrinter import TablePrinter


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
    coords: np.ndarray
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


"""
@dataclass
class Residues:
    residues: Dict
    psf_data: Dict

    def as_geom(self, with_link_atoms=True):
        geom = geom_from_residues(self.residues)
        atoms = geom.atoms
        coords3d = geom.coords3d
        bonds = np.array(self.psf_data["nbond"]["inds"], dtype=int).reshape(-1, 2)
        if with_link_atoms:
            link_hosts, link_atoms, link_coords3d = link_atoms_for_residues(
                self.residues, bonds, coords3d, atom_map
            )
        atoms += link_atoms
        coords3d = np.concatenate((coords3d, link_coords3d), axis=0)
        geom = Geometry(atoms, coords3d)
        return geom
"""


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


def residues_within_dist(
    residues,
    within_resid,
    within_dist,
    kind="com",
):
    def com_within():
        coms = {key: res.com for key, res in residues.items()}
        ref_com = coms[within_resid]
        res_ids_below_dist = [
            key
            for key, res in residues.items()
            if np.linalg.norm(coms[res.key] - ref_com) <= within_dist
        ]
        return res_ids_below_dist

    def atom_within():
        ref_res = residues[within_resid]
        ref_coords3d = ref_res.coords3d
        res_ids_below_dist = list()
        for key, res in residues.items():
            dist_vecs = ref_coords3d[:, None, :] - res.coords3d
            dists = np.linalg.norm(dist_vecs, axis=2)
            if (dists <= within_dist).any():
                res_ids_below_dist.append(key)
        return res_ids_below_dist

    if kind == "com":
        within_resid = com_within()
    elif kind == "atom":
        within_resid = atom_within()
        pass
    else:
        raise Exception(f"{kind=} is not supported!")

    residues_within = [residues[resid] for resid in within_resid]
    return residues_within


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
    residues,
    bonds,
    coords3d,
    atom_map,
    link_element="H",
    g=0.709,
    ignore_bonds=None,
):
    if ignore_bonds is None:
        ignore_bonds = list()
    atom_inds = list(it.chain(*[res.atom_indices for res in residues]))

    bond_dict = dict()
    for bond in bonds:
        from_, to_ = bond
        bond_dict.setdefault(from_, set()).add(to_)
        bond_dict.setdefault(to_, set()).add(from_)
    atom_set = set(atom_inds)

    for ibond in ignore_bonds:
        from_, to_ = ibond
        try:
            bond_dict[from_].remove(to_)
            bond_dict[to_].remove(from_)
        except KeyError:
            pass

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


def load_psf(psf_fn):
    if str(psf_fn).lower().endswith(".psf"):
        psf_data = parse_psf(psf_fn)
    else:
        with open(psf_fn, "rb") as handle:
            psf_data = pickle.load(handle)
    return psf_data


def cluster_from_psf_pdb(
    # psf_fn, pdb_fn, within_resid=None, within_dist=0.0, ref_residues=None, kind="atom,"
    psf_data,
    pdb_fn,
    within_resid=None,
    within_dist=0.0,
    ref_residues=None,
    kind="atom",
    ignore_bonds=None,
):
    atoms, coords, _, atom_map = parse_pdb(pdb_fn)
    coords3d = coords.reshape(-1, 3)
    print("Loaded PDB data.")

    residues = residues_from_psf(psf_data, atoms, coords3d, atom_map)

    # Select according to COM distance
    if within_resid and within_dist:
        # sel_residues = residues_within_com_dist(
        sel_residues = residues_within_dist(
            residues,
            within_resid,
            within_dist,
            kind=kind,
        )
    # Select residues according to provided ref_residues.
    elif ref_residues:
        sel_residues = [residues[ref_res] for ref_res in ref_residues]
    # Or complain!
    else:
        raise Exception(
            "Either 'within_resid' and 'within_dist' or 'residues' must be given!"
        )

    residue_charges = [res.charge for res in sel_residues]
    tot_charge = sum(residue_charges)
    charge_table = TablePrinter("# Residue Charge".split(), ("int", "str", "float"))
    charge_table.print_header()
    for i, res in enumerate(sel_residues):
        charge_table.print_row((i, f"{res.id}{res.name}", residue_charges[i]))
    print(f"Total charge of selected residues: {tot_charge: >8.6f} e")
    print()

    geom = geom_from_residues(sel_residues)
    print("Created cluster geometry.")
    bonds = np.array(psf_data["nbond"]["inds"], dtype=int).reshape(-1, 2)
    link_hosts, link_atoms, link_coords3d = link_atoms_for_residues(
        sel_residues,
        bonds,
        coords3d,
        atom_map,
        ignore_bonds=ignore_bonds,
    )
    sat_geom = Geometry(
        geom.atoms + link_atoms, np.concatenate((geom.coords3d, link_coords3d), axis=0)
    )
    nlink_atoms = len(link_atoms)
    header = ("#", "Org Host (1b)", "Link atom (0b)", "x", "y", "z")
    col_fmts = "int3 str str float float float".split()
    link_table = TablePrinter(header, col_fmts, width=14)
    print(f"Created satured geometry with {nlink_atoms} link atoms.")
    link_table.print_header()
    loffset = len(geom.atoms)
    for i in range(nlink_atoms):
        lx, ly, lz = link_coords3d[i]
        la = f"{i+loffset}{link_atoms[i]}"
        lh = int(link_hosts[i])
        ha = atoms[lh - 1]
        lha = f"{lh}{ha}"
        link_table.print_row((i, lha, la, lx, ly, lz))

    # Determine backbone atoms and/or link atoms and report their indices, so they can
    # be restrained in subsequent RMSD-biased optimizations. Also report the distance
    # between their positions and the COM of the reference residue.
    backbone_names = {"CA", "C", "O", "N"}  # HN, HA not included
    i = 0
    backbone_inds = list()
    backbone_com_dists = list()
    ref_com = None
    if within_resid is not None:
        ref_com = residues[within_resid].com
    for res in sel_residues:
        for atom in res.atoms:
            if atom.name in backbone_names:
                backbone_inds.append(i)
                if ref_com is not None:
                    bb_dist = np.linalg.norm(ref_com - atom.coords)
                    backbone_com_dists.append(bb_dist)
                    print(
                        f"{res.name: >4s}{res.id: <5d}, {atom.element: >2s}{atom.id: <5d}, "
                        f"type={atom.name: >4s}, id={i: >5d}, {bb_dist: >8.4f} au"
                    )
            i += 1
    backbone_inds = np.array(backbone_inds, dtype=int)
    backbone_com_dists = np.array(backbone_com_dists, dtype=float)

    return sat_geom, sel_residues, backbone_inds, backbone_com_dists
