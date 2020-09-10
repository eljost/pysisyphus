# [1] https://doi.org/10.1063/1.1515483 optimization review
# [2] https://doi.org/10.1063/1.471864 delocalized internal coordinates
# [3] https://doi.org/10.1016/0009-2614(95)00646-L lindh model hessian
# [4] 10.1002/(SICI)1096-987X(19990730)20:10<1067::AID-JCC9>3.0.CO;2-V
#     Handling of corner cases
# [5] https://doi.org/10.1063/1.462844

import itertools as it
import logging

import numpy as np
from scipy.spatial.distance import pdist, squareform

from pysisyphus.elem_data import VDW_RADII, COVALENT_RADII as CR
from pysisyphus.intcoords.derivatives import d2q_b, d2q_a, d2q_d
from pysisyphus.intcoords.findbonds import get_pair_covalent_radii
from pysisyphus.intcoords.fragments import merge_fragments
from pysisyphus.intcoords.backconversion import transform_int_step
from pysisyphus.intcoords.eval import eval_prim_internals, \
                                      calc_stretch, calc_bend, calc_dihedral
        

class RedundantCoords:

    RAD_175 = 3.05432619
    BEND_MIN_DEG = 15
    BEND_MAX_DEG = 180

    def __init__(
        self,
        atoms,
        cart_coords,
        bond_factor=1.3,
        prim_indices=None,
        define_prims=None,
        bonds_only=False,
        check_bends=True,
        check_dihedrals=False,
    ):
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

        self.logger = logging.getLogger("internal_coords")

        if prim_indices is None:
            self.set_primitive_indices(self.define_prims)
        else:
            to_arr = lambda _: np.array(list(_), dtype=int)
            bonds, bends, dihedrals = prim_indices
            # We accept all bond indices. What could possibly go wrong?! :)
            self.bond_indices = to_arr(bonds)
            valid_bends = [inds for inds in bends if self.is_valid_bend(inds)]
            self.bending_indices = to_arr(valid_bends)
            valid_dihedrals = [
                inds for inds in dihedrals if self.is_valid_dihedral(inds)
            ]
            self.dihedral_indices = to_arr(valid_dihedrals)

        if self.bonds_only:
            self.bending_indices = list()
            self.dihedral_indices = list()

        prim_ints = self.calculate(self.cart_coords)
        self._prim_internals = prim_ints
        self._prim_coords = np.array([prim.val for prim in self._prim_internals])

    def log(self, message):
        self.logger.debug(message)

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
        indices = [
            i
            for i, pi in enumerate(it.chain(*self.prim_indices))
            if set(pi) == prim_ind_set
        ]
        index = None
        try:
            index = indices[0]
        except IndexError:
            self.log(f"Primitive internal with indices {prim_ind} " "is not defined!")
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
                self.log(
                    "Error in calculation of 2nd derivative of primitive "
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
            cart_inds = list(it.chain(*[range(3 * i, 3 * i + 3) for i in pc.inds]))
            flat_inds = [
                row * size_ + col for row, col in it.product(cart_inds, cart_inds)
            ]
            K_flat[flat_inds] += dg
        K = K_flat.reshape(size_, size_)
        return K

    def transform_hessian(self, cart_hessian, int_gradient=None):
        """Transform Cartesian Hessian to internal coordinates."""
        if int_gradient is None:
            self.log(
                "Supplied 'int_gradient' is None. K matrix will be zero, "
                "so derivatives of the Wilson-B-matrix are neglected in "
                "the hessian transformation."
            )
        K = self.get_K_matrix(int_gradient)
        return self.Bt_inv.dot(cart_hessian - K).dot(self.B_inv)

    def backtransform_hessian(self, redund_hessian, int_gradient=None):
        """Transform Hessian in internal coordinates to Cartesians."""
        if int_gradient is None:
            self.log(
                "Supplied 'int_gradient' is None. K matrix will be zero, "
                "so derivatives of the Wilson-B-matrix are neglected in "
                "the hessian transformation."
            )
        K = self.get_K_matrix(int_gradient)
        return self.B.T.dot(redund_hessian).dot(self.B) + K

    def project_hessian(self, H, shift=1000):
        """Expects a hessian in internal coordinates. See Eq. (11) in [1]."""
        P = self.P
        return P.dot(H).dot(P) + shift * (np.eye(P.shape[0]) - P)

    def project_vector(self, vector):
        """Project supplied vector onto range of B."""
        return self.P.dot(vector)

    def connect_fragments(self, cdm, fragments):
        """Determine the smallest interfragment bond for a list
        of fragments and a condensed distance matrix."""
        dist_mat = squareform(cdm)
        interfragment_indices = list()
        for frag1, frag2 in it.combinations(fragments, 2):
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
        hydrogen_inds = [i for i, a in enumerate(self.atoms) if a.lower() == "h"]
        x_inds = [
            i for i, a in enumerate(self.atoms) if a.lower() in "n o f p s cl".split()
        ]
        for h_ind, x_ind in it.product(hydrogen_inds, x_inds):
            as_set = set((h_ind, x_ind))
            if as_set not in tmp_sets:
                continue
            # Check if distance of H to another electronegative atom Y is
            # greater than the sum of their covalent radii but smaller than
            # the 0.9 times the sum of their van der Waals radii. If the
            # angle X-H-Y is greater than 90° a hydrogen bond is asigned.
            y_inds = set(x_inds) - set((x_ind,))
            for y_ind in y_inds:
                y_atom = self.atoms[y_ind].lower()
                cov_rad_sum = CR["h"] + CR[y_atom]
                distance = calc_stretch(coords3d, (h_ind, y_ind))
                vdw = 0.9 * (VDW_RADII["h"] + VDW_RADII[y_atom])
                angle = calc_bend(coords3d, (x_ind, h_ind, y_ind))
                if (cov_rad_sum < distance < vdw) and (angle > np.pi / 2):
                    self.hydrogen_bond_indices.append((h_ind, y_ind))
                    self.log(
                        f"Added hydrogen bond between atoms {h_ind} "
                        f"({self.atoms[h_ind]}) and {y_ind} ({self.atoms[y_ind]})"
                    )
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
        atom_indices = list(it.combinations(range(len(coords3d)), 2))
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
            bond_indices = np.concatenate((bond_indices, self.hydrogen_bond_indices))

        # Merge bond index sets into fragments
        bond_ind_sets = [frozenset(bi) for bi in bond_indices]
        fragments = merge_fragments(bond_ind_sets)

        # Look for unbonded single atoms and create fragments for them.
        bonded_set = set(tuple(bond_indices.flatten()))
        unbonded_set = set(range(len(self.atoms))) - bonded_set
        fragments.extend([frozenset((atom,)) for atom in unbonded_set])
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
        (central,) = central_set
        return (terminal1, central, terminal2), central

    def is_valid_bend(self, bend_ind):
        val = calc_bend(self.c3d, bend_ind)
        deg = np.rad2deg(val)
        # Always return true if bends should not be checked
        return (not self.check_bends) or (self.BEND_MIN_DEG <= deg <= self.BEND_MAX_DEG)

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
            bis = np.concatenate(((self.bending_indices, define_bends)), axis=0)
            self.bending_indices = bis

    def is_valid_dihedral(self, dihedral_ind, thresh=1e-6):
        # Check for linear atoms
        first_angle = calc_bend(self.c3d, dihedral_ind[:3])
        second_angle = calc_bend(self.c3d, dihedral_ind[1:])
        pi_thresh = np.pi - thresh
        return (abs(first_angle) < pi_thresh) and (abs(second_angle) < pi_thresh)

    def set_dihedral_indices(self, define_dihedrals=None):
        dihedrals = list()

        def set_dihedral_index(dihedral_ind):
            dihed = tuple(dihedral_ind)
            # Check if this dihedral is already present
            if (dihed in dihedrals) or (dihed[::-1] in dihedrals):
                return
            # Assure that the angles are below 175° (3.054326 rad)
            if not self.is_valid_dihedral(dihedral_ind, thresh=0.0873):
                self.log(
                    f"Skipping generation of dihedral {dihedral_ind} "
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
                dihed = calc_dihedral(coords3d, dihedral_ind)
                if not np.isnan(dihed):
                    improper_dihedrals.append(dihedral_ind)
                else:
                    self.log(f"Dihedral {dihedral_ind} is undefinied. Skipping it!")

        # Now try to create the remaining improper dihedrals.
        if (len(self.atoms) >= 4) and (len(self.dihedral_indices) == 0):
            for improp in improper_dihedrals:
                set_dihedral_index(improp)
            self.log(
                "Permutational symmetry not considerd in "
                "generation of improper dihedrals."
            )

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
            by_prim_type[len_ - 2].append(item)
        return by_prim_type

    def set_primitive_indices(self, define_prims=None):
        stretches, bends, dihedrals = self.sort_by_prim_type(define_prims)
        self.set_bond_indices(stretches)
        self.set_bending_indices(bends)
        self.set_dihedral_indices(dihedrals)

    def set_prim_internals(self, prim_internals):
        self._prim_internals = list(it.chain(*prim_internals))
        self._prim_coords = np.array([prim.val for prim in self._prim_internals])

    def calculate(self, coords, attr=None):
        prim_internals = eval_prim_internals(coords, self.prim_indices)
        self.bonds, self.bends, self.dihedrals = prim_internals

        if attr:
            return np.array([getattr(ic, attr) for ic in prim_internals])

        return list(it.chain(*prim_internals))

    def transform_int_step(self, int_step, pure=False):
        new_prim_ints, cart_step, failed = transform_int_step(
            int_step,
            self.cart_coords,
            self._prim_coords,
            self.B_prim,
            self.prim_indices,
            logger=self.logger,
        )
        if not pure:
            self.set_prim_internals(new_prim_ints)
        return cart_step

    def __str__(self):
        bonds = len(self.bond_indices)
        bends = len(self.bending_indices)
        dihedrals = len(self.dihedral_indices)
        name = self.__class__.__name__
        return f"{name}({bonds} bonds, {bends} bends, {dihedrals} dihedrals)"
