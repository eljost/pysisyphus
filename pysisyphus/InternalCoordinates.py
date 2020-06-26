#!/usr/bin/env python3

# [1] https://doi.org/10.1063/1.1515483 optimization review
# [2] https://doi.org/10.1063/1.471864 delocalized internal coordinates
# [3] https://doi.org/10.1016/0009-2614(95)00646-L lindh model hessian
# [4] 10.1002/(SICI)1096-987X(19990730)20:10<1067::AID-JCC9>3.0.CO;2-V
#     Handling of corner cases
# [5] https://doi.org/10.1063/1.462844

from collections import namedtuple
from functools import reduce
import itertools as it
import logging

import numpy as np
from scipy.spatial.distance import pdist, squareform

from pysisyphus.constants import BOHR2ANG
from pysisyphus.elem_data import VDW_RADII, COVALENT_RADII as CR
from pysisyphus.intcoords.derivatives import d2q_b, d2q_a, d2q_d
from pysisyphus.intcoords.exceptions import NeedNewInternalsException
from pysisyphus.intcoords.findbonds import get_pair_covalent_radii
from pysisyphus.intcoords.fragments import merge_fragments


PrimitiveCoord = namedtuple(
                    "PrimitiveCoord",
                    "inds val grad",
)


class RedundantCoords:

    RAD_175 = 3.05432619
    BEND_MIN_DEG = 15
    BEND_MAX_DEG = 180

    def __init__(self, atoms, cart_coords, bond_factor=1.3,
                 prim_indices=None, define_prims=None, bonds_only=False,
                 check_bends=True, check_dihedrals=False):
        self.atoms = atoms
        self._cart_coords = cart_coords
        self.bond_factor = bond_factor
        self.define_prims = define_prims
        self.bonds_only = bonds_only
        self.check_bends = check_bends
        self.check_dihedrals = check_dihedrals

        self._B_prim = None
        self.bond_indices = list()
        self.bending_indices = list()
        self.dihedral_indices = list()
        self.hydrogen_bond_indices = list()

        if prim_indices is None:
            self.set_primitive_indices(self.define_prims)
        else:
            to_arr = lambda _: np.array(list(_), dtype=int)
            bonds, bends, dihedrals = prim_indices
            # We accept all bond indices. What could possibly go wrong?! :)
            self.bond_indices = to_arr(bonds)
            valid_bends = [inds for inds in bends
                           if self.is_valid_bend(inds)]
            self.bending_indices = to_arr(valid_bends)
            valid_dihedrals = [inds for inds in dihedrals if
                               self.is_valid_dihedral(inds)]
            self.dihedral_indices = to_arr(valid_dihedrals)

        if self.bonds_only:
            self.bending_indices = list()
            self.dihedral_indices = list()
        self._prim_internals = self.calculate(self.cart_coords)
        self._prim_coords = np.array([pc.val for pc in self._prim_internals])

    def log(self, message):
        logger = logging.getLogger("internal_coords")
        logger.debug(message)

    @property
    def prim_indices(self):
        return [self.bond_indices, self.bending_indices, self.dihedral_indices]

    @property
    def prim_indices_set(self):
        return set([tuple(prim_ind) for prim_ind in it.chain(*self.prim_indices)])

    @property
    def prim_coords(self):
        if self._prim_coords is None:
           self._prim_coords = np.array(
                            [pc.val for pc in self.calculate(self.cart_coords)]
            )
        return self._prim_coords

    @property
    def cart_coords(self):
        return self._cart_coords

    @cart_coords.setter
    def cart_coords(self, cart_coords):
        self._cart_coords = cart_coords
        self._B_prim = None

    @property
    def coords(self):
        return self.prim_coords

    @property
    def coord_indices(self):
        ic_ind_tuples = [tuple(ic.inds) for ic in self._prim_internals]
        return {ic_inds: i for i, ic_inds in enumerate(ic_ind_tuples)}

    @property
    def dihed_start(self):
        return len(self.bond_indices) + len(self.bending_indices)

    def get_index_of_prim_coord(self, prim_ind):
        """Index of primitive internal for the given atom indices.

        TODO: simplify this so when we get a prim_ind of len 2
        (bond) we don't have to check the bending and dihedral indices."""
        prim_ind_set = set(prim_ind)
        indices = [i for i, pi in enumerate(it.chain(*self.prim_indices))
                 if set(pi) == prim_ind_set]
        index = None
        try:
            index = indices[0]
        except IndexError:
            self.log(f"Primitive internal with indices {prim_ind} "
                      "is not defined!")
        return index

    @property
    def c3d(self):
        return self.cart_coords.reshape(-1, 3)

    @property
    def B_prim(self):
        """Wilson B-Matrix"""
        if self._B_prim is None:
            self._B_prim = np.array([c.grad for c in self.calculate(self.cart_coords)])

        return self._B_prim

    @property
    def B(self):
        """Wilson B-Matrix"""
        return self.B_prim

    @property
    def Bt_inv(self):
        """Transposed generalized inverse of the Wilson B-Matrix."""
        B = self.B
        return np.linalg.pinv(B.dot(B.T)).dot(B)

    @property
    def B_inv(self):
        """Generalized inverse of the Wilson B-Matrix."""
        B = self.B
        return B.T.dot(np.linalg.pinv(B.dot(B.T)))

    @property
    def P(self):
        """Projection matrix onto B. See [1] Eq. (4)."""
        return self.B.dot(self.B_inv)

    def transform_forces(self, cart_forces):
        """Combination of Eq. (9) and (11) in [1]."""
        return self.Bt_inv.dot(cart_forces)

    def get_K_matrix(self, int_gradient=None):
        if int_gradient is not None:
            assert len(int_gradient) == len(self._prim_internals)
        size_ = self.cart_coords.size
        if int_gradient is None:
            return np.zeros((size_, size_))

        dg_funcs = {
            2: d2q_b,
            3: d2q_a,
            4: d2q_d,
        }
        def grad_deriv_wrapper(inds):
            coords_flat = self.c3d[inds].flatten()
            dgrad = dg_funcs[len(inds)](*coords_flat)
            return dgrad

        K_flat = np.zeros(size_ * size_)
        for pc, int_grad_item in zip(self._prim_internals, int_gradient):
            # Contract with gradient
            try:
                dg = int_grad_item * grad_deriv_wrapper(pc.inds)
            except (ValueError, ZeroDivisionError) as err:
                self.log( "Error in calculation of 2nd derivative of primitive "
                         f"internal {pc.inds}."
                )
                continue
            # Depending on the type of internal coordinate dg is a flat array
            # of size 36 (stretch), 81 (bend) or 144 (torsion).
            #
            # An internal coordinate contributes to an element K[j, k] of the
            # K matrix if the cartesian coordinate indices j and k belong to an
            # atom that contributes to the respective internal coordinate.
            #
            # As for now we build up the K matrix as flat array. To add the dg
            # entries at the appropriate places in K_flat we have to calculate
            # the corresponding flat indices of dg in K_flat.
            cart_inds = list(it.chain(*[range(3*i,3*i+3) for i in pc.inds]))
            flat_inds = [row*size_ + col for row, col in it.product(cart_inds, cart_inds)]
            K_flat[flat_inds] += dg
        K = K_flat.reshape(size_, size_)
        return K

    def transform_hessian(self, cart_hessian, int_gradient=None):
        """Transform Cartesian Hessian to internal coordinates."""
        if int_gradient is None:
            self.log("Supplied 'int_gradient' is None. K matrix will be zero, "
                     "so derivatives of the Wilson-B-matrix are neglected in "
                     "the hessian transformation."
            )
        K = self.get_K_matrix(int_gradient)
        return self.Bt_inv.dot(cart_hessian-K).dot(self.B_inv)

    def backtransform_hessian(self, redund_hessian, int_gradient=None):
        """Transform Hessian in internal coordinates to Cartesians."""
        if int_gradient is None:
            self.log("Supplied 'int_gradient' is None. K matrix will be zero, "
                     "so derivatives of the Wilson-B-matrix are neglected in "
                     "the hessian transformation."
            )
        K = self.get_K_matrix(int_gradient)
        return self.B.T.dot(redund_hessian).dot(self.B) + K

    def project_hessian(self, H, shift=1000):
        """Expects a hessian in internal coordinates. See Eq. (11) in [1]."""
        P = self.P
        return P.dot(H).dot(P) + shift*(np.eye(P.shape[0]) - P)

    def project_vector(self, vector):
        """Project supplied vector onto range of B."""
        return self.P.dot(vector)

    def connect_fragments(self, cdm, fragments):
        """Determine the smallest interfragment bond for a list
        of fragments and a condensed distance matrix."""
        dist_mat = squareform(cdm)
        interfragment_indices = list()
        for frag1, frag2 in it.combinations(fragments, 2):
            arr1 = np.array(list(frag1))[None,:]
            arr2 = np.array(list(frag2))[:,None]
            indices = [(i1, i2) for i1, i2 in it.product(frag1, frag2)]
            distances = np.array([dist_mat[ind] for ind in indices])
            min_index = indices[distances.argmin()]
            interfragment_indices.append(min_index)
        # Or as Philipp proposed: two loops over the fragments and only
        # generate interfragment distances. So we get a full matrix with
        # the original indices but only the required distances.
        return interfragment_indices

    def set_hydrogen_bond_indices(self, bond_indices):
        coords3d = self.cart_coords.reshape(-1, 3)
        tmp_sets = [frozenset(bi) for bi in bond_indices]
        # Check for hydrogen bonds as described in [1] A.1 .
        # Find hydrogens bonded to small electronegative atoms X = (N, O
        # F, P, S, Cl).
        hydrogen_inds = [i for i, a in enumerate(self.atoms)
                         if a.lower() == "h"]
        x_inds = [i for i, a in enumerate(self.atoms)
                  if a.lower() in "n o f p s cl".split()]
        hydrogen_bond_inds = list()
        for h_ind, x_ind in it.product(hydrogen_inds, x_inds):
            as_set = set((h_ind, x_ind))
            if not as_set in tmp_sets:
                continue
            # Check if distance of H to another electronegative atom Y is
            # greater than the sum of their covalent radii but smaller than
            # the 0.9 times the sum of their van der Waals radii. If the
            # angle X-H-Y is greater than 90° a hydrogen bond is asigned.
            y_inds = set(x_inds) - set((x_ind, ))
            for y_ind in y_inds:
                y_atom = self.atoms[y_ind].lower()
                cov_rad_sum = CR["h"] + CR[y_atom]
                distance = self.calc_stretch(coords3d, (h_ind, y_ind))
                vdw = 0.9 * (VDW_RADII["h"] + VDW_RADII[y_atom])
                angle = self.calc_bend(coords3d, (x_ind, h_ind, y_ind))
                if (cov_rad_sum < distance < vdw) and (angle > np.pi/2):
                    self.hydrogen_bond_indices.append((h_ind, y_ind))
                    self.log(f"Added hydrogen bond between atoms {h_ind} "
                             f"({self.atoms[h_ind]}) and {y_ind} ({self.atoms[y_ind]})")
        self.hydrogen_bond_indices = np.array(self.hydrogen_bond_indices)

    def set_bond_indices(self, define_bonds=None, factor=None):
        """
        Default factor of 1.3 taken from [1] A.1.
        Gaussian uses somewhat less, like 1.2, or different radii than we do.
        """
        bond_factor = factor if factor else self.bond_factor
        coords3d = self.cart_coords.reshape(-1, 3)
        # Condensed distance matrix
        cdm = pdist(coords3d)
        # Generate indices corresponding to the atom pairs in the
        # condensed distance matrix cdm.
        atom_indices = list(it.combinations(range(len(coords3d)),2))
        atom_indices = np.array(atom_indices, dtype=int)
        cov_rad_sums = get_pair_covalent_radii(self.atoms)
        cov_rad_sums *= bond_factor
        bond_flags = cdm <= cov_rad_sums
        bond_indices = atom_indices[bond_flags]

        if define_bonds:
            bond_indices = np.concatenate(((bond_indices, define_bonds)), axis=0)

        self.bare_bond_indices = bond_indices

        # Look for hydrogen bonds
        self.set_hydrogen_bond_indices(bond_indices)
        if self.hydrogen_bond_indices.size > 0:
            bond_indices = np.concatenate((bond_indices,
                                           self.hydrogen_bond_indices))

        # Merge bond index sets into fragments
        bond_ind_sets = [frozenset(bi) for bi in bond_indices]
        fragments = merge_fragments(bond_ind_sets)

        # Look for unbonded single atoms and create fragments for them.
        bonded_set = set(tuple(bond_indices.flatten()))
        unbonded_set = set(range(len(self.atoms))) - bonded_set
        fragments.extend(
            [frozenset((atom, )) for atom in unbonded_set]
        )
        self.fragments = fragments

        # Check if there are any disconnected fragments. If there are some
        # create interfragment bonds between all of them.
        if len(fragments) != 1:
            interfragment_inds = self.connect_fragments(cdm, fragments)
            bond_indices = np.concatenate((bond_indices, interfragment_inds))

        self.bond_indices = bond_indices

    def are_parallel(self, vec1, vec2, angle_ind=None, thresh=1e-6):
        dot = max(min(vec1.dot(vec2), 1), -1)
        rad = np.arccos(dot)#vec1.dot(vec2))
        # angle > 175°
        if abs(rad) > self.RAD_175:
            # self.log(f"Nearly linear angle {angle_ind}: {np.rad2deg(rad)}")
            ind_str = f" ({angle_ind})" if (angle_ind is not None) else ""
            self.log(f"Nearly linear angle{ind_str}: {np.rad2deg(rad)}")
        return abs(rad) > (np.pi - thresh)

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

    def is_valid_bend(self, bend_ind):
        val = self.calc_bend(self.c3d, bend_ind)
        deg = np.rad2deg(val)
        # Always return true if bends should not be checked
        return (
            not self.check_bends) or (self.BEND_MIN_DEG <= deg <= self.BEND_MAX_DEG
        )

    def set_bending_indices(self, define_bends=None):
        bond_sets = {frozenset(bi) for bi in self.bond_indices}
        for bond_set1, bond_set2 in it.combinations(bond_sets, 2):
            union = bond_set1 | bond_set2
            if len(union) == 3:
                as_tpl, _ = self.sort_by_central(bond_set1, bond_set2)
                if not self.is_valid_bend(as_tpl):
                    self.log(f"Didn't create bend {list(as_tpl)}")
                             # f" with value of {deg:.3f}°")
                    continue
                self.bending_indices.append(as_tpl)
        self.bending_indices = np.array(self.bending_indices, dtype=int)

        if define_bends:
            bis = np.concatenate(( (self.bending_indices, define_bends)), axis=0)
            self.bending_indices = bis

    def is_valid_dihedral(self, dihedral_ind, thresh=1e-6):
        # Check for linear atoms
        first_angle = self.calc_bend(self.c3d, dihedral_ind[:3])
        second_angle = self.calc_bend(self.c3d, dihedral_ind[1:])
        pi_thresh = np.pi - thresh
        return ((abs(first_angle) < pi_thresh)
                and (abs(second_angle) < pi_thresh)
        )

    def set_dihedral_indices(self, define_dihedrals=None):
        dihedrals = list()
        def set_dihedral_index(dihedral_ind):
            dihed = tuple(dihedral_ind)
            # Check if this dihedral is already present
            if (dihed in dihedrals) or (dihed[::-1] in dihedrals):
                return
            # Assure that the angles are below 175° (3.054326 rad)
            if not self.is_valid_dihedral(dihedral_ind, thresh=0.0873):
                self.log(f"Skipping generation of dihedral {dihedral_ind} "
                          "as some of the the atoms are (nearly) linear."
                )
                return
            self.dihedral_indices.append(dihedral_ind)
            dihedrals.append(dihed)

        improper_dihedrals = list()
        coords3d = self.cart_coords.reshape(-1, 3)
        for bond, bend in it.product(self.bond_indices, self.bending_indices):
            central = bend[1]
            bend_set = set(bend)
            bond_set = set(bond)
            # Check if the two sets share one common atom. If not continue.
            intersect = bend_set & bond_set
            if len(intersect) != 1:
                continue
            # When the common atom is a terminal atom of the bend, that is
            # it's not the central atom of the bend, we create a
            # proper dihedral. Before we create any improper dihedrals we
            # create these proper dihedrals.
            if central not in bond_set:
                # The new terminal atom in the dihedral is the one that
                # doesn' intersect.
                terminal = tuple(bond_set - intersect)[0]
                intersecting_atom = tuple(intersect)[0]
                if intersecting_atom == bend[0]:
                    dihedral_ind = [terminal] + bend.tolist()
                else:
                    dihedral_ind = bend.tolist() + [terminal]
                set_dihedral_index(dihedral_ind)
            # If the common atom is the central atom we try to form an out
            # of plane bend / improper torsion. They may be created later on.
            else:
                fourth_atom = list(bond_set - intersect)
                dihedral_ind = bend.tolist() + fourth_atom
                # This way dihedrals may be generated that contain linear
                # atoms and these would be undefinied. So we check for this.
                dihed = self.calc_dihedral(coords3d, dihedral_ind)
                if not np.isnan(dihed):
                    improper_dihedrals.append(dihedral_ind)
                else:
                    self.log(f"Dihedral {dihedral_ind} is undefinied. Skipping it!")

        # Now try to create the remaining improper dihedrals.
        if (len(self.atoms) >= 4) and (len(self.dihedral_indices) == 0):
            for improp in improper_dihedrals:
                set_dihedral_index(improp)
            self.log("Permutational symmetry not considerd in "
                            "generation of improper dihedrals.")

        self.dihedral_indices = np.array(self.dihedral_indices)

        if define_dihedrals:
            dis = np.concatenate(((self.dihedral_indices, define_dihedrals)), axis=0)
            self.dihedral_indices = dis

    def sort_by_prim_type(self, to_sort):
        by_prim_type = [[], [], []]
        if to_sort is None:
            to_sort = list()
        for item in to_sort:
            len_ = len(item)
            by_prim_type[len_-2].append(item)
        return by_prim_type

    def set_primitive_indices(self, define_prims=None):
        stretches, bends, dihedrals = self.sort_by_prim_type(define_prims)
        self.set_bond_indices(stretches)
        self.set_bending_indices(bends)
        self.set_dihedral_indices(dihedrals)

    def calculate(self, coords, attr=None):
        coords3d = coords.reshape(-1, 3)
        def per_type(func, ind):
            val, grad = func(coords3d, ind, True)
            return PrimitiveCoord(ind, val, grad)
        self.bonds = list()
        self.bends = list()
        self.dihedrals = list()
        for ind in self.bond_indices:
            bonds = per_type(self.calc_stretch, ind)
            self.bonds.append(bonds)
        for ind in self.bending_indices:
            bend = per_type(self.calc_bend, ind)
            self.bends.append(bend)
        for ind in self.dihedral_indices:
            dihedral = per_type(self.calc_dihedral, ind)
            self.dihedrals.append(dihedral)
        int_coords = self.bonds + self.bends + self.dihedrals
        if attr:
            return np.array([getattr(ic,attr) for ic in int_coords])
        return int_coords

    def calc_stretch(self, coords3d, bond_ind, grad=False):
        n, m = bond_ind
        bond = coords3d[m] - coords3d[n]
        bond_length = np.linalg.norm(bond)
        if grad:
            bond_normed = bond / bond_length
            row = np.zeros_like(coords3d)
            # 1 / -1 correspond to the sign factor [1] Eq. 18
            row[m,:] =  bond_normed
            row[n,:] = -bond_normed
            row = row.flatten()
            return bond_length, row
        return bond_length

    def calc_bend(self, coords3d, angle_ind, grad=False):
        m, o, n = angle_ind
        u_dash = coords3d[m] - coords3d[o]
        v_dash = coords3d[n] - coords3d[o]
        u_norm = np.linalg.norm(u_dash)
        v_norm = np.linalg.norm(v_dash)
        u = u_dash / u_norm
        v = v_dash / v_norm
        angle_rad = np.arccos(u.dot(v))
        if grad:
            # Eq. (24) in [1]
            if self.are_parallel(u, v, angle_ind):
                tmp_vec = np.array((1, -1, 1))
                par = self.are_parallel(u, tmp_vec) and self.are_parallel(v, tmp_vec)
                tmp_vec = np.array((-1, 1, 1)) if par else tmp_vec
                w_dash = np.cross(u, tmp_vec)
            else:
                w_dash = np.cross(u, v)
            w_norm = np.linalg.norm(w_dash)
            w = w_dash / w_norm
            uxw = np.cross(u, w)
            wxv = np.cross(w, v)

            row = np.zeros_like(coords3d)
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
            return angle_rad, row
        return angle_rad

    def calc_dihedral(self, coords3d, dihedral_ind, grad=False, cos_tol=1e-9):
        m, o, p, n = dihedral_ind
        u_dash = coords3d[m] - coords3d[o]
        v_dash = coords3d[n] - coords3d[p]
        w_dash = coords3d[p] - coords3d[o]
        u_norm = np.linalg.norm(u_dash)
        v_norm = np.linalg.norm(v_dash)
        w_norm = np.linalg.norm(w_dash)
        u = u_dash / u_norm
        v = v_dash / v_norm
        w = w_dash / w_norm
        phi_u = np.arccos(u.dot(w))
        phi_v = np.arccos(-w.dot(v))
        uxw = np.cross(u, w)
        vxw = np.cross(v, w)
        cos_dihed = uxw.dot(vxw)/(np.sin(phi_u)*np.sin(phi_v))

        # Restrict cos_dihed to [-1, 1]
        if cos_dihed >= 1 - cos_tol:
            dihedral_rad = 0
        elif cos_dihed <= -1 + cos_tol:
            dihedral_rad = np.arccos(-1)
        else:
            dihedral_rad = np.arccos(cos_dihed)

        if dihedral_rad != np.pi:
            # wxv = np.cross(w, v)
            # if wxv.dot(u) < 0:
            if vxw.dot(u) < 0:
                dihedral_rad *= -1
        if grad:
            row = np.zeros_like(coords3d)
            #                  |  m  |  n  |  o  |  p  |
            # ------------------------------------------
            # sign_factor(amo) |  1  |  0  | -1  |  0  | 1st term
            # sign_factor(apn) |  0  | -1  |  0  |  1  | 2nd term
            # sign_factor(aop) |  0  |  0  |  1  | -1  | 3rd term
            # sign_factor(apo) |  0  |  0  | -1  |  1  | 4th term
            sin2_u = np.sin(phi_u)**2
            sin2_v = np.sin(phi_v)**2
            first_term  = uxw/(u_norm*sin2_u)
            second_term = vxw/(v_norm*sin2_v)
            third_term  = uxw*np.cos(phi_u)/(w_norm*sin2_u)
            fourth_term = -vxw*np.cos(phi_v)/(w_norm*sin2_v)
            row[m,:] = first_term
            row[n,:] = -second_term
            row[o,:] = -first_term + third_term - fourth_term
            row[p,:] = second_term - third_term + fourth_term
            row = row.flatten()
            return dihedral_rad, row
        return dihedral_rad

    def update_internals(self, new_cartesians, prev_internals):
        new_internals = self.calculate(new_cartesians, attr="val")
        internal_diffs = np.array(new_internals - prev_internals)
        _, _, dihedrals = self.prim_indices
        dihedral_diffs = internal_diffs[-len(dihedrals):]
        # Find differences that are shifted by 2*pi
        shifted_by_2pi = np.abs(np.abs(dihedral_diffs) - 2*np.pi) < np.pi/2
        org = dihedral_diffs.copy()
        new_dihedrals = new_internals[-len(dihedrals):]
        new_dihedrals[shifted_by_2pi] -= 2*np.pi * np.sign(dihedral_diffs[shifted_by_2pi])
        new_internals[-len(dihedrals):] = new_dihedrals
        return new_internals

    def dihedrals_are_valid(self, cart_coords):
        _, _, dihedrals = self.prim_indices

        def collinear(v1, v2, thresh=1e-4):
            # ~4e-5 corresponds to 179.5°
            return 1 - abs(v1.dot(v2)) <= thresh

        coords3d = cart_coords.reshape(-1, 3)
        def check(indices):
            m, o, p, n = indices
            u_dash = coords3d[m] - coords3d[o]
            v_dash = coords3d[n] - coords3d[p]
            w_dash = coords3d[p] - coords3d[o]
            u_norm = np.linalg.norm(u_dash)
            v_norm = np.linalg.norm(v_dash)
            w_norm = np.linalg.norm(w_dash)
            u = u_dash / u_norm
            v = v_dash / v_norm
            w = w_dash / w_norm

            valid = not (collinear(u, w) or collinear(v, w))
            return valid

        all_valid = all([check(indices) for indices in dihedrals])
        return all_valid

    def transform_int_step(self, step, cart_rms_thresh=1e-6):
        """This is always done in primitive internal coordinates so care
        has to be taken that the supplied step is given in primitive internal
        coordinates."""

        remaining_int_step = step
        cur_cart_coords = self.cart_coords.copy()
        cur_internals = self.prim_coords
        target_internals = cur_internals + step
        B_prim = self.B_prim

        # Bt_inv may be overriden in other coordiante systems so we
        # calculate it 'manually' here.
        Bt_inv_prim = np.linalg.pinv(B_prim.dot(B_prim.T)).dot(B_prim)

        last_rms = 9999
        prev_internals = cur_internals
        self.backtransform_failed = True
        for i in range(25):
            cart_step = Bt_inv_prim.T.dot(remaining_int_step)
            # Recalculate exact Bt_inv every cycle. Costly.
            # cart_step = self.Bt_inv.T.dot(remaining_int_step)
            cart_rms = np.sqrt(np.mean(cart_step**2))
            # Update cartesian coordinates
            cur_cart_coords += cart_step
            # Determine new internal coordinates
            new_internals = self.update_internals(cur_cart_coords, prev_internals)
            remaining_int_step = target_internals - new_internals
            internal_rms = np.sqrt(np.mean(remaining_int_step**2))
            self.log(f"Cycle {i}: rms(Δcart)={cart_rms:1.4e}, "
                     f"rms(Δinternal) = {internal_rms:1.5e}"
            )

            # This assumes the first cart_rms won't be > 9999 ;)
            if (cart_rms < last_rms):
                # Store results of the conversion cycle for laster use, if
                # the internal-cartesian-transformation goes bad.
                best_cycle = (cur_cart_coords.copy(), new_internals.copy())
                best_cycle_ind = i
            elif i != 0:
                # If the conversion somehow fails we return the step
                # saved above.
                self.log( "Internal to cartesian failed! Using from step "
                         f"from cycle {best_cycle_ind}."
                )
                cur_cart_coords, new_internals = best_cycle
                break
            else:
                raise Exception("Internal-cartesian back-transformation already "
                                "failed in the first step. Aborting!"
                )
            prev_internals = new_internals

            last_rms = cart_rms
            if cart_rms < cart_rms_thresh:
                self.log("Internal to cartesian transformation converged!")
                self.backtransform_failed = False
                break
            self._prim_coords = np.array(new_internals)

        if self.check_dihedrals and (not self.dihedrals_are_valid(cur_cart_coords)):
            raise NeedNewInternalsException(cur_cart_coords)

        self.log("")
        # Return the difference between the new cartesian coordinates that yield
        # the desired internal coordinates and the old cartesian coordinates.
        return cur_cart_coords - self.cart_coords

    def __str__(self):
        bonds = len(self.bond_indices)
        bends = len(self.bending_indices)
        dihedrals = len(self.dihedral_indices)
        name = self.__class__.__name__
        return f"{name}({bonds} bonds, {bends} bends, {dihedrals} dihedrals)"
