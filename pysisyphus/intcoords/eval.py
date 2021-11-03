import logging
import random

import numpy as np


class PrimInternal:
    def __init__(self, inds, val, grad=None):
        self.inds = inds
        self.val = val
        self.grad = grad

    def __str__(self):
        return f"PrimInternal({self.inds})"

    def __repr__(self):
        return f"PrimInternal({self.inds}, {self.val:.4f})"


def eval_primitives(coords3d, primitives):
    prim_internals = list()
    for primitive in primitives:
        value, gradient = primitive.calculate(coords3d, gradient=True)
        prim_internal = PrimInternal(primitive.indices, value, gradient)
        prim_internals.append(prim_internal)
    return prim_internals


def eval_B(coords3d, primitives):
    prim_internals = eval_primitives(coords3d, primitives)
    return np.array([prim_int.grad for prim_int in prim_internals])


def check_primitives(coords3d, primitives, B=None, thresh=1e-6, logger=None):
    assert len(primitives) > 0

    def log(msg, level=logging.DEBUG):
        if logger is not None:
            logger.log(level, msg)

    if B is None:
        B = eval_B(coords3d, primitives)
    G = B.T.dot(B)
    w, v = np.linalg.eigh(G)
    nonzero_inds = np.abs(w) > thresh
    # More coordinates may be expected when collinear atoms are present.
    expected = coords3d.size - 6
    nonzero_num = sum(nonzero_inds)
    missing = expected - nonzero_num
    if missing > 0:
        log(
            "Not enough internal coordinates defined! Expected at least "
            f"{expected} nonzero eigenvalues. Found only {nonzero_num}!"
        )
    nonzero_w = w[nonzero_inds]
    # Condition number
    kappa = abs(nonzero_w.max() / nonzero_w.min())
    log(f"Condition number of B^T.B=G: {kappa:.4e}")
    return missing + 1, kappa


def augment_primitives(missing_prims, coords3d, prim_indices, fragments):
    add_bonds = list()
    add_bends = list()
    add_dihedrals = list()

    fragment_tpls = [tuple(fragment) for fragment in fragments]
    if len(fragments) > 1:
        bond_inds = prim_indices[0]
        bond_sets = [frozenset(bond) for bond in bond_inds]
        while missing_prims > 0:
            random.shuffle(fragment_tpls)
            frag1, frag2 = fragment_tpls[:2]
            atom1 = random.choice(frag1)
            atom2 = random.choice(frag2)
            bond_set = frozenset((atom1, atom2))
            if (bond_set not in bond_sets) and (bond_set not in add_bonds):
                add_bonds.append(list(bond_set))
                missing_prims -= 1
    return add_bonds, add_bends, add_dihedrals
