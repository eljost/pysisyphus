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

from pysisyphus.elem_data import VDW_RADII, COVALENT_RADII as CR

#logging.basicConfig(level=logging.INFO)

PrimitiveCoord = namedtuple("PrimitiveCoord", "inds val grad")


class RedundantCoords:

    def __init__(self, atoms, cart_coords):
        self.atoms = atoms
        self.cart_coords = cart_coords

        self.bond_indices = list()
        self.bending_indices = list()
        self.dihedral_indices = list()
        self.hydrogen_bond_indices = list()

        self.set_primitive_indices()
        self._prim_coords = self.calculate(self.cart_coords)
        self.set_rho()

    @property
    def prim_indices(self):
        return [self.bond_indices, self.bending_indices, self.dihedral_indices]

    @property
    def B(self):
        """Wilson B-Matrix"""
        return np.array([c.grad for c in self.calculate(self.cart_coords)])

    @property
    def Bt_inv(self):
        """Generalized inverse of the transposed Wilson B-Matrix."""
        B = self.B
        return np.linalg.pinv(B.dot(B.T)).dot(B)

    @property
    def P(self):
        """Projection matrix onto B. See [1] Eq. (4)."""
        return self.B.dot(self.Bt_inv.T)

    def transform_forces(self, cart_forces):
        """Combination of Eq. (9) and (11) in [1]."""
        return self.P.dot(self.Bt_inv.dot(cart_forces))

    def project_hessian(self, H):
        """Expects a hessian in internal coordinates. See Eq. (11) in [1]."""
        P = self.P
        return self.P.dot(H).dot(P)

    @property
    def coords(self):
        return np.array([pc.val for pc in self.calculate(self.cart_coords)])

    def set_rho(self):
        """Calculated rho values as required for the Lindh model hessian
        as described in [3], similar to pyberny.
        Instead of using the tabulated r_ref,ij values we will use the covalent
        radii as in pyberny. The tabulated r_ref,ij value for two carbons
        (2nd period) is 2.87 Bohr. Carbons covalent radius is ~ 1.44 Bohr,
        so two times it is 2.88 Bohr which fits nicely with the tabulate value.
        Hydrogens covalent radius is 0.59 bohr, so C-H gives 2.03 Bohr
        (tabulated 2.10). If values for elements > 3rd are requested the alpha
        values for the 3rd period will be (re)used.
        """
        first_period = "h he".split()
        def get_alpha(atom1, atom2):
            if (atom1 in first_period) and (atom2 in first_period):
                return 1.
            elif (atom1 in first_period) or (atom2 in first_period):
                return 0.3949
            else:
                return 0.28
        atoms = [a.lower() for a in self.atoms]
        alphas = [get_alpha(a1, a2)
                  for a1, a2 in it.combinations(atoms, 2)]
        cov_radii = np.array([CR[a.lower()] for a in atoms])
        rref = np.array([r1+r2
                         for r1, r2 in it.combinations(cov_radii, 2)])
        coords3d = self.cart_coords.reshape(-1, 3)
        cdm = pdist(coords3d)
        # It shouldn't be a problem that the diagonal is 0 because
        # no primitive internal coordinates will ever access a diagonal
        # element.
        self.rho = squareform(np.exp(alphas*(rref**2-cdm**2)))

    def get_initial_hessian(self):
        """
        k_dict = {
            2: 0.45,
            3: 0.15,
            4: 0.005,
        }
        k_diag = list()
        for primitive in self._prim_coords:
            rho_product = 1
            for i1, i2 in it.combinations(primitive.inds, 2):
                rho_product *= self.rho[i1, i2]
            k_diag.append(k_dict[len(primitive.inds)] * rho_product)
        return np.diagflat(k_diag)
        """
        k_dict = {
            2: 0.5,
            3: 0.2,
            4: 0.1,
        }
        k_diag = [k_dict[len(prim.inds)] for prim in self._prim_coords]
        logging.warning("Using simple 0.5/0.2/0.1 model hessian!")
        return np.diagflat(k_diag)

    def merge_fragments(self, fragments):
        """Merge a list of sets recursively. Pop the first element
        of the list and check if it intersects with one of the remaining
        elements. If yes, delete the intersecting set from the list, form
        the union of both sets and append it at the end of the list.
        If the popped set doesn't intersect with any of the remaining sets
        append the it at the end of the list."""
        if len(fragments) == 1:
            return fragments
        popped = fragments.pop(0)
        for frag in fragments:
            if popped & frag:
                fragments.remove(frag)
                fragments.append(popped | frag)
                return self.merge_fragments(fragments)
        fragments.append(popped)
        return fragments

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
                    logging.debug("Added hydrogen bond between {h_ind} and {y_ind}")
        self.hydrogen_bond_indices = np.array(self.hydrogen_bond_indices)

    def set_bond_indices(self, factor=1.3):
        """
        Default factor of 1.3 taken from [1] A.1.
        Gaussian uses somewhat less, like 1.2, or different radii than we do.
        """
        coords3d = self.cart_coords.reshape(-1, 3)
        # Condensed distance matrix
        cdm = pdist(coords3d)
        # Generate indices corresponding to the atom pairs in the
        # condensed distance matrix cdm.
        atom_indices = list(it.combinations(range(len(coords3d)),2))
        atom_indices = np.array(atom_indices, dtype=int)
        cov_rad_sums = list()
        for i, j in atom_indices:
            atom1 = self.atoms[i].lower()
            atom2 = self.atoms[j].lower()
            cov_rad1 = CR[atom1]
            cov_rad2 = CR[atom2]
            cov_rad_sum = factor * (cov_rad1 + cov_rad2)
            cov_rad_sums.append(cov_rad_sum)
        cov_rad_sums = np.array(cov_rad_sums)
        bond_flags = cdm <= cov_rad_sums
        bond_indices = atom_indices[bond_flags]

        # Look for hydrogen bonds
        self.set_hydrogen_bond_indices(bond_indices)
        if self.hydrogen_bond_indices.size > 0:
            bond_indices = np.concatenate((bond_indices,
                                           self.hydrogen_bond_indices))

        # Merge bond index sets into fragments
        bond_ind_sets = [frozenset(bi) for bi in bond_indices]
        fragments = self.merge_fragments(bond_ind_sets)

        # Look for unbonded single atoms and create fragments for them.
        bonded_set = reduce(lambda x, y: set(x) | set(y), bond_indices)
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
        for bond_set1, bond_set2 in it.combinations(bond_sets, 2):
            union = bond_set1 | bond_set2
            if len(union) == 3:
                as_tpl, _ = self.sort_by_central(bond_set1, bond_set2)
                self.bending_indices.append(as_tpl)
        self.bending_indices = np.array(self.bending_indices, dtype=int)

    def is_valid_dihedral(self, dihedral_ind, thresh=1e-6):
        coords3d = self.cart_coords.reshape(-1, 3)
        # Check for linear atoms
        first_angle = self.calc_bend(coords3d, dihedral_ind[:3])
        second_angle = self.calc_bend(coords3d, dihedral_ind[1:])
        pi_thresh = np.pi - thresh
        return ((abs(first_angle) < pi_thresh)
                and (abs(second_angle) < pi_thresh)
        )

    def set_dihedral_indices(self):
        dihedral_sets = list()
        def set_dihedral_index(dihedral_ind):
            dihedral_set = set(dihedral_ind)
            # Check if this dihedral is already present
            if dihedral_set in dihedral_sets:
                return
            # Assure that the angles are below 175° (3.054326 rad)
            if not self.is_valid_dihedral(dihedral_ind, thresh=0.0873):
                logging.warning("Skipping generation of dihedral "
                               f"{dihedral_ind} as some of the the atoms "
                                "are linear."
                )
                return
            self.dihedral_indices.append(dihedral_ind)
            dihedral_sets.append(dihedral_set)

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
                improper_dihedrals.append(dihedral_ind)

        # Now try to create the remaining improper dihedrals.
        if (len(self.atoms) >= 4) and (len(self.dihedral_indices) == 0):
            for improp in improper_dihedrals:
                set_dihedral_index(improp)
            logging.warning("Permutational symmetry not considerd in "
                            "generation of improper dihedrals.")

        self.dihedral_indices = np.array(self.dihedral_indices)

    def set_primitive_indices(self):
        self.set_bond_indices()
        self.set_bending_indices()
        self.set_dihedral_indices()

    def calculate(self, coords, attr=None):
        coords3d = coords.reshape(-1, 3)
        def per_type(func, ind):
            val, grad = func(coords3d, ind, True)
            return PrimitiveCoord(ind, val, grad)
        int_coords = list()
        for ind in self.bond_indices:
            int_coords.append(per_type(self.calc_stretch, ind))
        for ind in self.bending_indices:
            int_coords.append(per_type(self.calc_bend, ind))
        for ind in self.dihedral_indices:
            int_coords.append(per_type(self.calc_dihedral, ind))
        if attr:
            return np.array([getattr(ic,attr) for ic in int_coords])
        return int_coords

    def calculate_val_diffs(self, coords1, coords2):
        vals1 = np.array(self.calculate(coords1, attr="val"))
        vals2 = np.array(self.calculate(coords2, attr="val"))
        return vals1-vals2

    def calc_stretch(self, coords, bond_ind, grad=False):
        n, m = bond_ind
        bond = coords[m] - coords[n]
        bond_length = np.linalg.norm(bond)
        if grad:
            bond_normed = bond / bond_length
            row = np.zeros_like(coords)
            # 1 / -1 correspond to the sign factor [1] Eq. 18
            row[m,:] =  bond_normed
            row[n,:] = -bond_normed
            row = row.flatten()
            return bond_length, row
        return bond_length

    def calc_bend(self, coords, angle_ind, grad=False):
        def are_parallel(vec1, vec2, thresh=1e-6):
            dot = max(min(vec1.dot(vec2), 1), -1)
            rad = np.arccos(dot)#vec1.dot(vec2))
            # angle > 175°
            if abs(rad) > (np.pi - 0.088):
                logging.warning(f"Nearly linear angle {angle_ind}: {np.rad2deg(rad)}")
            return abs(rad) > (np.pi - thresh)
        m, o, n = angle_ind
        u_dash = coords[m] - coords[o]
        v_dash = coords[n] - coords[o]
        u_norm = np.linalg.norm(u_dash)
        v_norm = np.linalg.norm(v_dash)
        u = u_dash / u_norm
        v = v_dash / v_norm
        angle_rad = np.arccos(u.dot(v))
        if grad:
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
            return angle_rad, row
        return angle_rad

    def calc_dihedral(self, coords, dihedral_ind, grad=False):
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
        dihedral_rad = np.arccos(cos_dihed)
        if grad:
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
            return dihedral_rad, row
        return dihedral_rad

    def get_internal_diffs(self, new_coords, last_vals):
        new_internals = self.calculate(new_coords)
        diffs = list()
        for new_i, last_val in zip(new_internals, last_vals):
            diff = new_i.val - last_val
            # Remove multiples of 2pi from angles and dihedrals
            if new_i.inds.size > 2:
                pre_diff = diff
                diff = np.sign(diff) * (abs(diff) % (2*np.pi))
            diffs.append(diff)
        return np.array(diffs)

    def transform_int_step(self, step, cart_rms_thresh=1e-6):
        def rms(coords1, coords2):
            return np.sqrt(np.mean((coords1-coords2)**2))
        last_step = step
        last_coords = self.cart_coords.copy()
        last_vals = self.calculate(last_coords, attr="val")
        last_rms = None
        full_cart_step = np.zeros_like(self.cart_coords)
        Bt_inv = self.Bt_inv

        first_cart_step = None
        for i in range(25):
            cart_step = Bt_inv.T.dot(last_step)
            new_coords = last_coords + cart_step
            cart_rms = rms(last_coords, new_coords)
            if i == 0:
                # Store the first converted cartesian step if the
                # transformation goes wrong.
                first_cart_step = cart_step
            elif cart_rms > last_rms:
                # If the conversion somehow fails we return the step
                # saved above.
                full_cart_step = first_cart_step
                break
            full_cart_step += cart_step
            new_vals = self.calculate(new_coords)
            int_diffs = self.get_internal_diffs(new_coords, last_vals)
            tmp_vals = self.calculate(last_coords)

            last_step -= int_diffs
            last_coords = new_coords
            last_vals += int_diffs
            last_rms = cart_rms
            logging.info(f"Cycle {i}: rms(ΔCart) = {cart_rms:1.4e}")
            if cart_rms < cart_rms_thresh:
                logging.info("Internal to cartesian transformation converged!")
                break
        return full_cart_step


class DelocalizedCoords(RedundantCoords):
    def __init__(self, geom):
        super().__init__(geom)

    def set_delocalized_vectors(self, thresh=1e-6):
        """See [5] between Eq. (7) and Eq. (8) for advice regarding
        the threshold."""
        G = self.B_prim.dot(self.B_prim.T)
        w, v = np.linalg.eigh(G)
        #print(w)
        #print(w.shape)
        #print(v.T)
        #import pdb; pdb.set_trace()
        non_zero_inds = np.where(abs(w) > thresh)
        degrees_of_freedom = 3*len(self.atoms)-6
        assert(len(non_zero_inds[0]) == degrees_of_freedom)
        self.delocalized_vectors = v[:,non_zero_inds[0]]
        # Eq. 3 in [2], transformation of B to the active coordinate set
        self.B = self.delocalized_vectors.T.dot(self.B_prim)
        self.Bt_inv = np.linalg.pinv(self.B.dot(self.B.T)).dot(self.B)

    def get_delocalized(self):
        primitives = self.get_primitives()
        return primitives.dot(self.delocalized_vectors)

    """
    @property
    def G(self):
        B = self.B
        return B.dot(B.T)
    """

