#!/usr/bin/env python3

# [1] https://doi.org/10.1063/1.1515483 optimization review
# [2] https://doi.org/10.1063/1.471864 delocalized internal coordinates
# [3] https://doi.org/10.1016/0009-2614(95)00646-L lindh model hessian
# [4] 10.1002/(SICI)1096-987X(19990730)20:10<1067::AID-JCC9>3.0.CO;2-V
#     Handling of corner cases
# [5] https://doi.org/10.1063/1.462844

import itertools as it
import logging
import typing

import attr
import numpy as np
from scipy.spatial.distance import squareform

from pysisyphus.constants import BOHR2ANG
from pysisyphus.elem_data import VDW_RADII, COVALENT_RADII as CR
from pysisyphus.intcoords import Bend, LinearBend, Stretch, Torsion
from pysisyphus.intcoords.derivatives import d2q_b, d2q_a, d2q_d
from pysisyphus.intcoords.findbonds import get_bond_sets
from pysisyphus.intcoords.fragments import merge_fragments


@attr.s(auto_attribs=True)
class PrimitiveCoord:
    inds : typing.List[int]
    val : float
    grad : np.ndarray


class RedundantCoords:

    RAD_175 = 3.05432619
    BEND_MIN_DEG = 15
    LIN_BEND_DEG = 170

    def __init__(self, atoms, cart_coords, bond_factor=1.3,
                 prim_indices=None, define_prims=None, bonds_only=False,
                 check_bends=True, linear_bend_deg=170, make_complement=True):
        """Redundant internal coordinate handling.

        Parameters
        ----------
        atoms : iterable
            Iterable of length N, containing  element symbols.
        coords : 1d iterable
            1d iterable of length 3N, containing the cartesian coordinates
            of N atoms.
        bond_factor : float, optional
            Scaling factor for the sum of the covalent radii. Used for
            bond/stretch setup. Bigger values will lead to increased number
            of bond definitions.
        prim_indices : iterable of shape (3, (bonds, bends, dihedrals))
            Iterable containing definitions for primitive internal coordinates.
            Three items are expected in the iterable: the first one should
            contain integer pairs, defining bonds between atoms, the next one
            should contain integer triples containing bends, and finally
            integer quadrupel for dihedrals.
        define_prims : list of lists, optional
            Define additional primitives.
        bonds_only : bool, optional
            Whether only bonds/stretches shall be defined. Default is False.
        bond_factor : float, default=1.3
            Scaling factor of the the pair covalent radii used for bond detection.
        check_bends : bool, default=True
            Whether to check for invalid bends. If disabled these bends will
            be defined nonetheless.
        """
        self.atoms = atoms
        self.cart_coords = cart_coords
        self.bond_factor = bond_factor
        self.define_prims = define_prims
        self.bonds_only = bonds_only
        self.check_bends = check_bends

        self.linear_bend_deg = float(linear_bend_deg)
        self.make_complement = bool(make_complement)

        self.stretch_indices = list()
        self.bend_indices = list()
        self.torsion_indices = list()
        self.hydrogen_bond_indices = list()

        if prim_indices is None:
            self.set_primitive_indices(self.define_prims)
        else:
            to_arr = lambda _: np.array(list(_), dtype=int)
            bonds, bends, dihedrals = prim_indices
            # We accept all bond indices. What could possibly go wrong?! :)
            self.stretch_indices = to_arr(bonds)
            valid_bends = [inds for inds in bends
                           if self.is_valid_bend(inds)]
            self.bend_indices = to_arr(valid_bends)
            valid_dihedrals = [inds for inds in dihedrals if
                               self.is_valid_dihedral(inds)]
            self.torsion_indices = to_arr(valid_dihedrals)

        if self.bonds_only:
            self.bend_indices = list()
            self.torsion_indices = list()

        # Create the actual primitive internal coordinate (PIC) objects
        self._primitives = self.get_primitives(
                                self.stretch_indices,
                                self.bend_indices,
                                self.torsion_indices,
                                self.cart_coords,
        )
        self._prim_internals = self.calculate(self.cart_coords)
        self._prim_coords = np.array([pc.val for pc in self._prim_internals])

    def log(self, message):
        logger = logging.getLogger("internal_coords")
        logger.debug(message)

    @property
    def prim_indices(self):
        return [self.stretch_indices, self.bend_indices, self.torsion_indices]

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
    def coords(self):
        return self.prim_coords

    @property
    def coord_indices(self):
        ic_ind_tuples = [tuple(ic.inds) for ic in self._prim_internals]
        return {ic_inds: i for i, ic_inds in enumerate(ic_ind_tuples)}

    @property
    def dihed_start(self):
        return len(self.stretch_indices) + len(self.bend_indices)

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
        return np.array([c.grad for c in self.calculate(self.cart_coords)])

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
        if int_gradient is None:
            self.log("Supplied 'int_gradient' is None. K matrix will be zero, "
                     "so derivatives of the Wilson-B-matrix are neglected in "
                     "the hessian transformation."
            )
        K = self.get_K_matrix(int_gradient)
        return self.Bt_inv.dot(cart_hessian-K).dot(self.B_inv)

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
        self.log(f"Connecting fragments")
        for i, frag in enumerate(fragments):
            self.log(f"\t{i: >5d}: {tuple(frag)}")
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

    def set_hydrogen_bond_indices(self, stretch_indices):
        coords3d = self.cart_coords.reshape(-1, 3)
        tmp_sets = [frozenset(bi) for bi in stretch_indices]
        # Check for hydrogen bonds as described in [1] A.1 .
        # Find hydrogens bonded to small electronegative atoms X = (N, O
        # F, P, S, Cl).
        hydrogen_inds = [i for i, a in enumerate(self.atoms)
                         if a.lower() == "h"]
        x_inds = [i for i, a in enumerate(self.atoms)
                  if a.lower() in "n o f p s cl".split()]
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
                distance = np.linalg.norm(coords3d[h_ind] - coords3d[y_ind])
                vdw = 0.9 * (VDW_RADII["h"] + VDW_RADII[y_atom])
                angle = Bend._calculate(coords3d, (x_ind,  h_ind, y_ind))
                if (cov_rad_sum < distance < vdw) and (angle > np.pi/2):
                    self.hydrogen_bond_indices.append((h_ind, y_ind))
                    self.log(f"Added hydrogen bond between atoms {h_ind} "
                             f"({self.atoms[h_ind]}) and {y_ind} ({self.atoms[y_ind]})")
        self.hydrogen_bond_indices = np.array(self.hydrogen_bond_indices)

    def get_stretch_indices(self, define_bonds=None, bond_factor=None):
        """ Default factor of 1.3 taken from [1] A.1."""

        bond_factor = float(bond_factor) if bond_factor else self.bond_factor
        coords3d = self.cart_coords.reshape(-1, 3)

        # Set up bond indices
        stretch_indices, cdm, cbm = get_bond_sets(
                                        self.atoms,
                                        coords3d,
                                        bond_factor=bond_factor,
                                        return_cdm=True,
                                        return_cbm=True,
        )
        if define_bonds:
            stretch_indices = np.concatenate(((stretch_indices, define_bonds)), axis=0)
            for from_, to_ in define_bonds:
                cdm[from_, to_] = 1
        self.bond_matrix = squareform(cbm)

        # Bond indices without interfragment bonds and/or hydrogen bonds
        self.bare_stretch_indices = stretch_indices

        # Look for hydrogen bonds
        self.set_hydrogen_bond_indices(stretch_indices)
        if self.hydrogen_bond_indices.size > 0:
            stretch_indices = np.concatenate((stretch_indices,
                                              self.hydrogen_bond_indices))

        # Merge bond index sets into fragments
        bond_ind_sets = [frozenset(bi) for bi in stretch_indices]
        fragments = merge_fragments(bond_ind_sets)

        # Look for unbonded single atoms and create fragments for them.
        bonded_set = set(tuple(stretch_indices.flatten()))
        unbonded_set = set(range(len(self.atoms))) - bonded_set
        fragments.extend(
            [frozenset((atom, )) for atom in unbonded_set]
        )
        self.fragments = fragments

        # Check if there are any disconnected fragments. If there are some
        # create interfragment bonds between all of them.
        if len(fragments) != 1:
            interfragment_inds = self.connect_fragments(cdm, fragments)
            stretch_indices = np.concatenate((stretch_indices, interfragment_inds))

        self.stretch_indices = stretch_indices

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
        val = Bend._calculate(self.c3d, bend_ind)
        deg = np.rad2deg(val)
        # Always return true if bends should not be checked
        return (not self.check_bends) or (self.BEND_MIN_DEG <= deg)

    def get_bend_indices(self, define_bends=None):
        bond_sets = {frozenset(bi) for bi in self.stretch_indices}
        for bond_set1, bond_set2 in it.combinations(bond_sets, 2):
            union = bond_set1 | bond_set2
            if len(union) == 3:
                as_tpl, _ = self.sort_by_central(bond_set1, bond_set2)
                if not self.is_valid_bend(as_tpl):
                    self.log(f"Didn't create bend {list(as_tpl)}")
                    continue
                self.bend_indices.append(as_tpl)
        self.bend_indices = np.array(self.bend_indices, dtype=int)

        if define_bends:
            self.bend_indices = np.concatenate(
                                    ((self.bend_indices, define_bends)),
                                    axis=0
            )

    def is_valid_dihedral(self, dihedral_ind, thresh=1e-6):
        # Check for linear atoms
        first_angle = Bend._calculate(self.c3d, dihedral_ind[:3])
        second_angle = Bend._calculate(self.c3d, dihedral_ind[1:])
        pi_thresh = np.pi - thresh
        return ((abs(first_angle) < pi_thresh)
                and (abs(second_angle) < pi_thresh)
        )

    def get_torsion_indices(self, define_dihedrals=None):
        dihedral_sets = list()
        def set_dihedral_index(dihedral_ind):
            dihedral_set = set(dihedral_ind)
            # Check if this dihedral is already present
            if dihedral_set in dihedral_sets:
                return
            # Assure that the angles are below 175° (3.054326 rad)
            if not self.is_valid_dihedral(dihedral_ind, thresh=0.0873):
                self.log(f"Did not create dihedral {dihedral_ind} as some "
                          "vectors are (nearly) colinear."
                )
                return
            self.torsion_indices.append(dihedral_ind)
            dihedral_sets.append(dihedral_set)

        improper_dihedrals = list()
        coords3d = self.cart_coords.reshape(-1, 3)
        for bond, bend in it.product(self.stretch_indices, self.bend_indices):
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
                dihed = Torsion._calculate(coords3d=coords3d, indices=dihedral_ind)
                if not np.isnan(dihed):
                    improper_dihedrals.append(dihedral_ind)
                else:
                    self.log(f"Dihedral {dihedral_ind} is undefinied. Skipping it!")

        # Now try to create the remaining improper dihedrals.
        if (len(self.atoms) >= 4) and (len(self.torsion_indices) == 0):
            for improp in improper_dihedrals:
                set_dihedral_index(improp)
            self.log("Permutational symmetry not considerd in "
                            "generation of improper dihedrals.")

        self.torsion_indices = np.array(self.torsion_indices)

        if define_dihedrals:
            dis = np.concatenate(((self.torsion_indices, define_dihedrals)), axis=0)
            self.torsion_indices = dis

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
        self.get_stretch_indices(stretches)
        self.get_bend_indices(bends)
        self.get_torsion_indices(dihedrals)

    def get_primitives(self, stretch_indices, bend_indices, torsion_indices,
                       coords):
        coords3d = coords.reshape(-1, 3)
        classes = {
            2: Stretch,
            3: Bend,
            # 3: LinearBend,  # Will be handled explicitly
            4: Torsion,
        }

        primitives = list()
        for inds in it.chain(stretch_indices, bend_indices, torsion_indices):
            prim_kwargs = {
                "indices": inds,
                "periodic": len(inds) == 4,
            }
            cls = classes[len(inds)]
            val = cls._calculate(coords3d, inds)

            # Check for linear angles
            linear = (
                len(inds) == 3
                and self.linear_bend_deg > 0
                and np.rad2deg(val) >= self.linear_bend_deg
                # No need for linear bend if already enough bonds present
                and sum(self.bond_matrix[inds[1]]) < 5
            )
            if linear:
                self.log(f"Bend {inds}={np.rad2deg(val):.1f}° is (close to) linear. "
                          "Creating linear bend & complement.")
                # Create LinearBend instead of regular Bend
                cls = LinearBend

            # Create primitive coordinate and append
            prim = cls(**prim_kwargs)
            primitives.append(prim)

            if linear and self.make_complement:
                self.log(f"Created complement for Bend {inds}")
                prim_kwargs["complement"] = True
                prim = cls(**prim_kwargs)
                primitives.append(prim)

        self.log("Defined primitives")
        for i, p in enumerate(primitives):
            self.log(f"\t{i}: {p.indices}")
        return primitives

    def calculate(self, coords, attr=None):
        coords3d = coords.reshape(-1, 3)
        self.bonds = list()
        self.bends = list()
        self.dihedrals = list()

        lists = {
            2: self.bonds,
            3: self.bends,
            4: self.dihedrals,
        }

        pcs = list()
        for prim in self._primitives:
            val, grad = prim.calculate(coords3d, gradient=True)
            pc = PrimitiveCoord(prim.indices, val, grad)
            pcs.append(pc)
            lists[len(prim.indices)].append(pc)

        if attr:
            return np.array([getattr(pc, attr) for pc in pcs])
        return pcs

    def update_internals(self, new_cartesians, prev_internals):
        new_internals = self.calculate(new_cartesians, attr="val")
        internal_diffs = np.array(new_internals - prev_internals)
        _, _, dihedrals = self.prim_indices
        dihedral_diffs = internal_diffs[-len(dihedrals):]
        # Find differences that are shifted by 2*pi
        shifted_by_2pi = np.abs(np.abs(dihedral_diffs) - 2*np.pi) < np.pi/2
        new_dihedrals = new_internals[-len(dihedrals):]
        new_dihedrals[shifted_by_2pi] -= 2*np.pi * np.sign(dihedral_diffs[shifted_by_2pi])
        new_internals[-len(dihedrals):] = new_dihedrals
        return new_internals

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
        self.log("")
        return cur_cart_coords - self.cart_coords

    def __str__(self):
        bonds = len(self.stretch_indices)
        bends = len(self.bend_indices)
        dihedrals = len(self.torsion_indices)
        name = self.__class__.__name__
        return f"{name}({bonds} bonds, {bends} bends, {dihedrals} dihedrals)"
