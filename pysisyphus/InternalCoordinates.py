#!/usr/bin/env python3

# [1] https://doi.org/10.1063/1.1515483
# [2] https://doi.org/10.1063/1.471864 delocalized internal coordinates

import itertools
import logging

import numpy as np
from scipy.spatial.distance import pdist, squareform

from pysisyphus.helpers import geom_from_library
from pysisyphus.elem_data import COVALENT_RADII as CR


class InternalCoordinates:

    def __init__(self, geom):
        self.geom = geom

        self.bond_indices = list()
        self.bending_indices = list()
        self.dihedral_indices = list()

        self.set_bond_indices()
        self.set_bending_indices()
        self.set_dihedral_indices()
        B_prim = self.make_B_prim()
        self.set_B_inv(B_prim)

    def merge_fragments(self, fragments):
        """Merge a list of sets recursively. Pop the first element
        of the list and check if it intersects with one of the remaining
        elements. If yes, delete the intersecting set from the list, form
        the union of both sets and append it at the end of the list. If
        it doesn't intersect with any of the remaining sets append the
        popped set at the end of the list."""
        if len(fragments) == 1:
            return fragments
        popped = fragments.pop(0)
        for i, frag in enumerate(fragments):
            merged = popped & frag
            if merged:
                fragments.remove(frag)
                fragments.append(popped | frag)
                break
        if merged:
            return self.merge_fragments(fragments)
        fragments.append(popped)
        return fragments

    def connect_fragments(self, cdm, fragments):
        """Determine the smallest interfragment bond for a list
        of fragments and a condensed distance matrix."""
        dist_mat = squareform(cdm)
        interfragment_indices = list()
        for frag1, frag2 in itertools.combinations(fragments, 2):
            arr1 = np.array(list(frag1))[None,:]
            arr2 = np.array(list(frag2))[:,None]
            indices = [(i1, i2) for i1, i2 in itertools.product(frag1, frag2)]
            distances = np.array([dist_mat[ind] for ind in indices])
            min_index = indices[distances.argmin()]
            interfragment_indices.append(min_index)
        # Or as Philipp proposed: two loops over the fragments and only
        # generate interfragment distances. So we get a full matrix with
        # the original indices but only the required distances.
        return interfragment_indices

    def set_bond_indices(self, factor=1.3):
        """
        Default factor taken from [1] A.1.
        """
        coords = self.geom.coords.reshape(-1, 3)
        # Condensed distance matrix
        cdm = pdist(coords)
        # Generate indices corresponding to the atom pairs in
        # condensed distance matrix cdm.
        atom_indices = list(itertools.combinations(range(len(coords)),2))
        atom_indices = np.array(atom_indices, dtype=int)
        cov_rad_sums = list()
        for i, j in atom_indices:
            atom1 = self.geom.atoms[i].lower()
            atom2 = self.geom.atoms[j].lower()
            cov_rad1 = CR[atom1]
            cov_rad2 = CR[atom2]
            cov_rad_sum = factor * (cov_rad1 + cov_rad2)
            cov_rad_sums.append(cov_rad_sum)
        cov_rad_sums = np.array(cov_rad_sums)
        bond_flags = cdm <= cov_rad_sums
        bond_indices = atom_indices[bond_flags]

        # Check if there are any disconnected fragments
        bond_ind_sets = [frozenset(bi) for bi in bond_indices]
        fragments = self.merge_fragments(bond_ind_sets)
        if len(fragments) != 1:
            interfragment_inds = connect_fragments(cdm, fragments)
            bond_indices = np.concatenate((bond_indices, interfragment_inds))

        logging.warning("No check for single atoms!")
        logging.warning("No check for hydrogen bonds!")
        self.bond_indices = bond_indices

    def sort_by_central(self, set1, set2):
        """Determines a common index in two sets and returns a length 3
        tuple with the central index at the middle position and the two
        terminal indices as first and last indices."""
        central_set = set1 & set2
        union = set1 | set2
        assert len(central_set) == 1
        terminal1, terminal2 = union - central_set
        (central, ) = central_set
        return (terminal1, central, terminal2), central

    def set_bending_indices(self):
        bond_sets = {frozenset(bi) for bi in self.bond_indices}
        for bond_set1, bond_set2 in itertools.combinations(bond_sets, 2):
            union = bond_set1 | bond_set2
            if len(union) == 3:
                as_tpl, _ = self.sort_by_central(bond_set1, bond_set2)
                self.bending_indices.append(as_tpl)
        logging.warning("No check for (nearly) linear angles!"
                        "No additional orthogonal bending coordinates.")
        self.bending_indices = np.array(self.bending_indices, dtype=int)

    def set_dihedral_indices(self):
        dihedral_sets = list()
        for bond, bend in itertools.product(self.bond_indices,
                                            self.bending_indices):
            central = bend[1]
            bes = set((bend[0], bend[2]))
            bois = set(bond)
            if (len(bes & bois) == 1) and (central not in bois):
                (intersect,)  = set(bond) & set(bend)
                intersect_ind = list(bond).index(intersect)
                term_ind = 1 - intersect_ind
                terminal = bond[term_ind]
                if intersect == bend[0]:
                    dihedral_ind = [terminal] + list(bend)
                else:
                    dihedral_ind = list(bend) + [terminal]
                dihedral_set = set(dihedral_ind)
                if dihedral_set not in dihedral_sets:
                    self.dihedral_indices.append(dihedral_ind)
                    dihedral_sets.append(dihedral_set)
        logging.warning("No brute force method for molecules that should have a "
                        "dihedral but where no one can be found.")
        self.dihedral_indices = np.array(self.dihedral_indices)

    def get_bond_B(self, bond_ind):
        coords = self.geom.coords.reshape(-1, 3)
        n, m = bond_ind
        bond = coords[m] - coords[n]
        bond_length = np.linalg.norm(bond)
        bond_normed = bond / bond_length
        row = np.zeros_like(coords)
        # 1 / -1 correspond to the sign factor [1] Eq. 18
        row[m,:] = 1 * bond_normed
        row[n,:] = -1 * bond_normed
        row = row.flatten()
        return row

    def get_angle_B(self, angle_ind):
        def are_parallel(vec1, vec2, thresh=1e-6):
            rad = np.arccos(vec1.dot(vec2))
            return abs(rad) > (np.pi - thresh)
        coords = self.geom.coords.reshape(-1, 3)
        m, o, n = angle_ind
        u_dash = coords[m] - coords[o]
        v_dash = coords[n] - coords[o]
        u_norm = np.linalg.norm(u_dash)
        v_norm = np.linalg.norm(v_dash)
        u = u_dash / u_norm
        v = v_dash / v_norm
        rad = np.arccos(u.dot(v))
        # Eq. (24) in [1]
        if are_parallel(u, v):
            tmp_vec = np.array((1, -1, 1))
            par = are_parallel(u, tmp_vec) and are_parallel(v, tmp_vec)
            tmp_vec = np.array((-1, 1, 1)) if par else tmp_vec
            w_dash = np.cross(u, tmp_vec)
        else:
            w_dash = np.cross(u, v)
        w_norm = np.linalg.norm(w_dash)
        w = w_dash / w_norm
        uxw = np.cross(u, w)
        wxv = np.cross(w, v)

        row = np.zeros_like(coords)
        #                  |  m  |  n  |  o  |
        # -----------------------------------
        # sign_factor(amo) |  1  |  0  | -1  | first_term
        # sign_factor(ano) |  0  |  1  | -1  | second_term
        first_term = uxw / u_norm
        second_term = wxv / v_norm
        row[m,:] = first_term
        row[o,:] = -first_term - second_term
        row[n,:] = second_term
        row = row.flatten()
        return row

    def get_dihedral_B(self, dihedral_ind):
        coords = self.geom.coords.reshape(-1, 3)
        m, o, p, n = dihedral_ind
        u_dash = coords[m] - coords[o]
        v_dash = coords[n] - coords[p]
        w_dash = coords[p] - coords[o]
        u_norm = np.linalg.norm(u_dash)
        v_norm = np.linalg.norm(v_dash)
        w_norm = np.linalg.norm(w_dash)
        u = u_dash / u_norm
        v = v_dash / v_norm
        w = w_dash / w_norm
        phi_u = np.arccos(u.dot(w))
        phi_v = np.arccos(w.dot(v))
        uxw = np.cross(u, w)
        vxw = np.cross(v, w)
        cos_dihed = uxw.dot(vxw)/(np.sin(phi_u)*np.sin(phi_v))
        # Restrict cos_dihed to [-1, 1]
        cos_dihed = min(cos_dihed, 1)
        cos_dihed = max(cos_dihed, -1)
        dihedral = np.arccos(cos_dihed)

        row = np.zeros_like(coords)
        #                  |  m  |  n  |  o  |  p  |
        # ------------------------------------------
        # sign_factor(amo) |  1  |  0  | -1  |  0  | 1st term
        # sign_factor(apn) |  0  | -1  |  0  |  1  | 2nd term
        # sign_factor(aop) |  0  |  0  |  1  | -1  | 3rd term
        sin2_u = np.sin(phi_u)**2
        sin2_v = np.sin(phi_v)**2
        first_term  = uxw/(u_norm*sin2_u)
        second_term = vxw/(v_norm*sin2_v)
        third_term  = (uxw*np.cos(phi_u)/(w_norm*sin2_u)
                      -vxw*np.cos(phi_v)/(w_norm*sin2_v)
        )
        row[m,:] = first_term
        row[n,:] = -second_term
        row[o,:] = -first_term + third_term
        row[p,:] = second_term - third_term
        row = row.flatten()
        return row

    def make_B_prim(self, save=None, tm_format=False):
        bond_rows = [self.get_bond_B(ind) for ind in self.bond_indices]
        bend_rows = [self.get_angle_B(ind) for ind in self.bending_indices]
        dihed_rows = [self.get_dihedral_B(ind) for ind in self.dihedral_indices]

        B = np.array((*bond_rows, *bend_rows, *dihed_rows))
        if save:
            fmt = "% 1.4f"
            np.savetxt(save, B_mat, fmt=fmt)
            if tm_format:
                atm_fmt = " {: .04f}"
                tm_str = ""
                col = 1
                for row in B_mat:
                    per_atom = row.reshape(-1, 3)
                    tm_str += f"\ncolumn{col:3d}\n"
                    col += 1
                    atom = 1
                    for per_atom in row.reshape(-1, 3):
                        tm_str += f"atom {atom:3d}"

                        tm_str += (atm_fmt*3+"\n").format(*per_atom)
                        atom += 1
                with open(save + ".tm", "w") as handle:
                    handle.write(tm_str)
        return B

    def set_B_inv(self, B_prim):
        G = B_prim.dot(B_prim.T)
        w, v = np.linalg.eigh(G)
        #print(w)
        #print(w.shape)
        #print(v.T)
        import pdb; pdb.set_trace()
        non_zero_inds = np.where(abs(w) > 1e-6)
        degrees_of_freedom = 3*len(self.geom.atoms)-6
        assert(len(non_zero_inds[0]) == degrees_of_freedom)
        eigenvecs = v.T[non_zero_inds]
        # Eq. 3 in [2], transformation of B to the active coordinate set
        B = eigenvecs.dot(B_prim)
        self.B_inv = np.linalg.pinv(B.dot(B.T)).dot(B)

    # internal indices entkoppeln und in objekt speichern
    # bei bedarf b-matrix regenerieren

    def backtransform(self, int_step):
        def rms(coords1, coords2):
            return np.sqrt(np.mean((coords1-coords2)**2))
        #x_new = coords + B_mat_inv.dot(int_step)
        tmp_coords = self.geom.coords.copy()
        last_int_step = int_step
        while True:
            break
            cart_step =  self.B_inv.T.dot(last_int_step)
            new_coords = tmp_coords + cart_step
            #import pdb; pdb.set_trace()
        pass
